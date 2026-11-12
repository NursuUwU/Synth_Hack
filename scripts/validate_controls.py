"""Offline stress controls, deliberately separate from external clinical validation."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
from collections import Counter
from contextlib import nullcontext
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from genomic_detective.annotations import AnnotationCache, VepProvider
from genomic_detective.dataset import Dataset
from genomic_detective.demo_data import ensure_demo_data
from genomic_detective.evaluation import read_development_labels
from genomic_detective.knowledge import Knowledge
from genomic_detective.models import utc_now
from genomic_detective.molecular import MolecularKnowledge, default_knowledge
from genomic_detective.pipeline import Pipeline
from genomic_detective.ranking import RankingConfig

SEED = 'controls-v1-2026-11-12'
SCENARIOS = ('baseline', 'remove_known_causal', 'empty_positive_hpo',
             'half_positive_hpo', 'add_five_hpo', 'frequency_unavailable',
             'clinical_unavailable', 'causal_annotation_unavailable')
NOISE = ('HP:0001250', 'HP:0001513', 'HP:0001631', 'HP:0000098',
         'HP:0001252', 'HP:0000505', 'HP:0000077', 'HP:0000365',
         'HP:0001249', 'HP:0000407')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def transform(case, annotations, expected, scenario, knowledge):
    """Only the validation harness knows the label used for targeted deletion/masking."""
    if scenario not in SCENARIOS:
        raise ValueError('Unknown control scenario')
    if expected not in {v.id for v in case.variants}:
        raise ValueError('Known causal variant absent from original input')
    modified, changed = copy.deepcopy(case), copy.deepcopy(annotations)
    edits = {}
    if scenario == 'remove_known_causal':
        modified.variants = [v for v in modified.variants if v.id != expected]
        if not modified.variants:
            raise ValueError('Control requires at least one background variant')
        edits['removed_variant'] = expected
    elif scenario == 'empty_positive_hpo':
        edits['removed_hpo'] = sorted(modified.positive_hpo)
        modified.positive_hpo = {}
    elif scenario == 'half_positive_hpo':
        ordered = sorted(modified.positive_hpo,
                         key=lambda term: hashlib.sha256((SEED + case.id + term).encode()).hexdigest())
        keep = set(ordered[:math.ceil(len(ordered) / 2)])
        edits['removed_hpo'] = sorted(set(ordered) - keep)
        modified.positive_hpo = {k: v for k, v in modified.positive_hpo.items() if k in keep}
    elif scenario == 'add_five_hpo':
        extra = [t for t in NOISE if t in knowledge.parents and knowledge.ic(t) > 0
                 and t not in modified.positive_hpo and t not in modified.excluded_hpo][:5]
        if len(extra) != 5:
            raise ValueError('Cannot construct five non-overlapping HPO additions')
        edits['added_hpo'] = extra
        modified.positive_hpo.update({t: knowledge.labels[t] for t in extra})
    elif scenario == 'frequency_unavailable':
        for item in changed.values():
            item['gnomad_max_af'], item['frequencies'] = None, []
        edits['masked'] = 'All VEP gnomAD frequencies; no change to variants or consequences'
    elif scenario == 'clinical_unavailable':
        for item in changed.values():
            item['clinical'] = []
        edits['masked'] = 'All VEP clinical labels and bundled ClinVar records; literature retained'
    elif scenario == 'causal_annotation_unavailable':
        changed[expected] = {'variant_id': expected, 'assembly': 'GRCh38', 'status': 'unavailable',
                             'transcripts': [], 'gnomad_max_af': None, 'clinical': [], 'frequencies': []}
        edits['unavailable_variant'] = expected
    if scenario != 'baseline':
        modified.metadata['validation_control'] = scenario
        if scenario in {'empty_positive_hpo', 'half_positive_hpo', 'add_five_hpo'}:
            modified.metadata['phenotype_override'] = 'validation_control'
    modified.__post_init__()
    return modified, {v.id: changed[v.id] for v in modified.variants}, edits


class FrozenProvider:
    def __init__(self, annotations, source):
        self.annotations, self.source = annotations, source

    def annotate(self, variants, **kwargs):
        if not kwargs.get('offline', True):
            raise ValueError('Control runner never permits online annotation')
        return {v.id: copy.deepcopy(self.annotations[v.id]) for v in variants}, self.source


def summarize(rows):
    groups = []
    for model in sorted({r['model'] for r in rows}):
        for scenario in SCENARIOS:
            selected = [r for r in rows if r['model'] == model and r['scenario'] == scenario]
            if not selected:
                continue
            n = len(selected)
            positive = scenario != 'remove_known_causal'
            usable = [r for r in selected if r.get('status') not in {'error', 'incomplete'}]
            groups.append({'model': model, 'scenario': scenario, 'attempted': n,
                'completed': len(usable), 'statuses': dict(Counter(r['status'] for r in selected)),
                'candidate_count': sum(r['status'] == 'candidate' for r in selected),
                'candidate_wrong_first': sum(r['status'] == 'candidate' and r['expected_rank'] != 1 for r in selected),
                'top_1_all_attempted': sum(r.get('expected_rank') == 1 for r in usable) / n if positive else None,
                'top_3_all_attempted': sum(r.get('expected_rank') is not None and r['expected_rank'] <= 3 for r in usable) / n if positive else None,
                'top_5_all_attempted': sum(r.get('expected_rank') is not None and r['expected_rank'] <= 5 for r in usable) / n if positive else None,
                'mrr_all_attempted': sum(1 / r['expected_rank'] for r in usable if r.get('expected_rank')) / n if positive else None})
    return groups


def execute(cases, annotations, source, knowledge, configs, labels, output, progress=print):
    output.mkdir(parents=True, exist_ok=False)
    store = default_knowledge()
    masked_snapshot = copy.deepcopy(store.snapshot)
    masked_snapshot.update(records={}, scope='Artificial clinical masking control; not a new database release')
    masked = MolecularKnowledge(copy.deepcopy(store.rules), masked_snapshot)
    manifest = {'schema_version': 1, 'created_at': utc_now(), 'seed': SEED,
        'scope': 'Stress controls on known development cases; NOT independent clinical validation',
        'frozen_model_commit': 'd1df5c05c1796482c65306d6ede7745b9a877fb5',
        'scenarios': list(SCENARIOS), 'noise_hpo_candidates': list(NOISE),
        'configs': {name: config.dictionary() for name, config in configs.items()},
        'original_inputs_sha256': digest([asdict(case) for case in cases]),
        'original_annotations_sha256': digest(annotations),
        'original_sources': {'vep': source, 'hpo': knowledge.manifest, 'molecular': store.metadata},
        'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitations': ['Known labels are used only to construct targeted controls and evaluate results.',
                       'Added HPO terms are artificial perturbations, not verified patient findings.',
                       'Deleted causal variant does not establish a clinically negative case.',
                       'No threshold or evidence curation is changed in response to this run.']}
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    rows = []
    for model, config in configs.items():
        for scenario in SCENARIOS:
            for case in cases:
                expected = labels[case.id]['expected_answer']
                row = {'model': model, 'scenario': scenario, 'case_id': case.id, 'expected': expected}
                try:
                    changed_case, changed_annotations, edits = transform(case, annotations, expected, scenario, knowledge)
                    context = patch('genomic_detective.molecular.default_knowledge', return_value=masked) if scenario == 'clinical_unavailable' else nullcontext()
                    with context:
                        report = Pipeline(FrozenProvider(changed_annotations, source), knowledge, config).run([changed_case])[0]
                    target = next((c for c in report['candidates'] if c['variant_id'] == expected), None)
                    best = report['candidates'][0]
                    row.update(status=report['status'], prediction=report['prediction'],
                        expected_rank=target['rank'] if target else None,
                        expected_score=target['score'] if target else None,
                        top_variant=best['variant_id'], top_gene=best['gene'], top_score=best['score'],
                        annotation_coverage=report['annotation_coverage'], variant_count=report['variant_count'],
                        positive_hpo_count=len(changed_case.positive_hpo), edits=edits,
                        uncertainty_reasons=report['uncertainty_reasons'],
                        ranking_guard=report['ranking_guard']['status'],
                        presentation_guard=report['presentation']['guard']['status'],
                        input_sha256=report['input_sha256'], annotation_sha256=report['annotation_sha256'],
                        implementation_sha256=report['implementation_sha256'], config_sha256=report['config_sha256'],
                        molecular_sources=report['sources'].get('molecular'))
                except Exception as error:
                    row.update(status='error', error=f'{type(error).__name__}: {error}')
                rows.append(row)
                with (output / 'results.jsonl').open('a', encoding='utf-8') as stream:
                    stream.write(json.dumps(row, ensure_ascii=False) + '\n')
            if progress:
                counts = Counter(r['status'] for r in rows if r['model'] == model and r['scenario'] == scenario)
                progress(f'{model} / {scenario}: {dict(counts)}')
    result = {'manifest': manifest, 'summary': summarize(rows), 'cases': rows}
    (output / 'controls.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    lines = ['# Контрольные сценарии', '', manifest['scope'], '',
        'Места считаются среди всех попыток: неполные/ошибочные расчёты не засчитываются как успех.', '',
        '| Модель | Сценарий | Попытки | Статусы | candidate с неверным первым местом | Top-1 |',
        '|---|---|---:|---|---:|---:|']
    for group in result['summary']:
        metric = '—' if group['top_1_all_attempted'] is None else f"{group['top_1_all_attempted']:.1%}"
        lines.append(f"| {group['model']} | {group['scenario']} | {group['attempted']} | {group['statuses']} | {group['candidate_wrong_first']} | {metric} |")
    lines += ['', 'Удаление известного варианта — искусственный контроль, не клиническая специфичность.',
              'Статус candidate — прохождение исследовательских порогов, не подтверждённый диагноз.', '']
    (output / 'controls.md').write_text('\n'.join(lines), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/validation_controls_v1')
    args = parser.parse_args()
    ensure_demo_data(ROOT, progress=print)
    dataset = Dataset(ROOT / 'Project_3_genomic_detecive_for_students.zip')
    cache = AnnotationCache(ROOT / 'data/cache/vep.sqlite3')
    try:
        cases = dataset.load_cases('development')
        annotations, source = VepProvider(cache).annotate([v for c in cases for v in c.variants], offline=True)
    finally:
        dataset.close()
        cache.close()
    labels = {row['case_id']: row for row in read_development_labels(ROOT / 'Project_3_genomic_detecive_for_students.zip')}
    configs = {'v1.1': RankingConfig.load(ROOT / 'config/ranking-v1.1.json'),
               'v1.2': RankingConfig.load(ROOT / 'config/ranking.json')}
    result = execute(cases, annotations, source, Knowledge.load(ROOT / 'data/knowledge'), configs, labels, args.output)
    return int(any(r['status'] == 'error' for r in result['cases']))


if __name__ == '__main__':
    sys.exit(main())
