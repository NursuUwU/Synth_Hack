"""The only production module allowed to read development labels. Never imported by Pipeline."""
from __future__ import annotations

import csv
import io
import json
import hashlib
from dataclasses import asdict, replace
from pathlib import Path

from .dataset import Dataset


def read_development_labels(dataset_path: Path) -> list[dict]:
    dataset = Dataset(dataset_path)
    try:
        return list(csv.DictReader(io.StringIO(dataset.read('development_labels.tsv')), delimiter='\t'))
    finally:
        dataset.close()


def evaluate(dataset_path: Path, reports_directory: Path) -> dict:
    labels = read_development_labels(dataset_path)
    run = json.loads((reports_directory / 'run.json').read_text(encoding='utf-8'))
    active_cases = set(run['case_ids'])
    rows = []
    for label in labels:
        path = reports_directory / f'{label["case_id"]}.json'
        if label['case_id'] not in active_cases or not path.exists():
            rows.append({'case_id': label['case_id'], 'status': 'not_run', 'rank': None})
            continue
        report = json.loads(path.read_text(encoding='utf-8'))
        if report['input']['metadata'].get('phenotype_override'):
            raise ValueError('Evaluate original cases separately from phenotype overrides')
        expected = label['expected_answer']
        position = next((c['rank'] for c in report['candidates'] if c['variant_id'] == expected), None)
        rows.append({'case_id': label['case_id'], 'status': report['status'],
                     'expected': expected, 'rank': position, 'prediction': report['prediction'],
                     'annotation_coverage': report['annotation_coverage'],
                     'evidence_status': label.get('evidence'), 'source_url': label.get('source_url')})
    evaluable = [row for row in rows if row['status'] not in {'not_run', 'incomplete'}]
    count = len(evaluable)
    return {'scope': 'Open development set; not held-out clinical accuracy',
        'labeled_cases': len(labels), 'evaluated_complete_cases': count,
        'top_1': sum(row['rank'] == 1 for row in evaluable) / count if count else None,
        'top_3': sum(row['rank'] is not None and row['rank'] <= 3 for row in evaluable) / count if count else None,
        'top_5': sum(row['rank'] is not None and row['rank'] <= 5 for row in evaluable) / count if count else None,
        'mean_reciprocal_rank': sum(1 / row['rank'] for row in evaluable if row['rank']) / count if count else None,
        'prediction_exact_match': sum(row['prediction'] == row['expected'] for row in evaluable) / count if count else None,
        'negative_controls': 0,
        'limitations': ['No supplied negative controls; abstention is not validated by this metric.',
                       'Published labels were not independently clinically curated.',
                       'Incomplete runs are counted separately, never as successful evaluations.'],
        'cases': rows}


METHODS = {
    'combined': None,
    'molecular_only': {'phenotype': 0.0, 'molecular': 1.0, 'inheritance': 0.0},
    'phenotype_only': {'phenotype': 1.0, 'molecular': 0.0, 'inheritance': 0.0},
    'without_phenotype': {'phenotype': 0.0, 'molecular': 0.35 / 0.45, 'inheritance': 0.10 / 0.45},
}


def compare_experts(cases, annotations: dict, knowledge, config, dataset_path: Path, progress=None):
    """Recompute every gene/disease hypothesis before reading evaluation labels.

    This compares rankings, not calibrated decision thresholds. Zeroing a column in
    saved winning hypotheses would miss alternative genes/diseases and is avoided.
    """
    from .guards import check_ranking
    from .literature import default_registry
    from .models import CONTEXT_MODES, utc_now
    from .ranking import rank

    literature = default_registry()

    rankings, configs, included, excluded = {}, {}, [], []
    for case in cases:
        if case.metadata.get('split') != 'development' or case.metadata.get('phenotype_override'):
            raise ValueError('Benchmark requires original development cases only')
        missing = [v.id for v in case.variants if annotations.get(v.id, {}).get('status') != 'annotated']
        informative = any(knowledge.resolve(term) in knowledge.parents and knowledge.ic(term) > 0
                          for term in case.positive_hpo)
        if missing or not knowledge.diseases or not informative:
            excluded.append({'case_id': case.id, 'missing_annotations': len(missing),
                             'informative_phenotype': informative})
        else:
            included.append(case)
    for method, weights in METHODS.items():
        # Inheritance remaining in without_phenotype can still use disease assignments
        # from HPO. Molecular-only is the strict phenotype-independent comparator.
        method_config = replace(config, weights=dict(config.weights if weights is None else weights))
        method_config.validate()
        configs[method] = method_config.dictionary() | {'weights': dict(method_config.weights)}
        rankings[method] = {}
        for case in included:
            applied = literature if method_config.molecular_mode in CONTEXT_MODES else None
            candidates = rank(case, annotations, knowledge, method_config, applied)
            check_ranking(case, candidates, annotations, method_config, applied)
            rankings[method][case.id] = candidates
        if progress:
            progress(f'Benchmark: {method}, {len(included)} cases')

    # This is the first access to labels; none were available to the ranker.
    labels = {row['case_id']: row for row in read_development_labels(dataset_path)}
    rows = []
    for case in included:
        if case.id not in labels:
            raise ValueError(f'Missing development label for {case.id}')
        expected = labels[case.id]['expected_answer']
        if expected == 'no_causal_variant_in_vcf':
            raise ValueError('Rank-only benchmark does not evaluate negative labels')
        positions, ties, target_genes = {}, {}, {}
        for method in METHODS:
            candidates = rankings[method][case.id]
            target = next((candidate for candidate in candidates if candidate['variant_id'] == expected), None)
            positions[method] = target['rank'] if target else None
            target_genes[method] = target['gene'] if target else None
            ties[method] = sum(abs(candidate['score'] - target['score']) <= 1e-12
                               for candidate in candidates) if target else 0
        rows.append({'case_id': case.id, 'expected_variant': expected,
                     'expected_gene_from_label': labels[case.id].get('causal_gene'),
                     'ranks': positions, 'equal_score_group_size': ties,
                     'target_selected_gene': target_genes})
    metrics = {}
    for method in METHODS:
        positions = [row['ranks'][method] for row in rows]
        n = len(positions)
        metrics[method] = {
            'cases': n,
            **{f'top_{k}': sum(p is not None and p <= k for p in positions) / n if n else None for k in (1, 3, 5)},
            'mrr': sum(1 / p for p in positions if p) / n if n else None,
        }
    inputs = {case.id: {'variants': [asdict(v) for v in case.variants],
                       'positive_hpo': case.positive_hpo, 'excluded_hpo': case.excluded_hpo}
              for case in included}
    output = {'created_at': utc_now(), 'scope': 'Open development ranking comparison; not independent clinical accuracy',
        'configs': configs, 'metrics': metrics, 'cases': rows, 'excluded_cases': excluded,
        'input_sha256': hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest(),
        'config_sha256': config.digest,
        'implementation_sha256': hashlib.sha256(b''.join(p.name.encode() + b'\0' + p.read_bytes()
            for p in sorted(Path(__file__).parent.glob('*.py')))).hexdigest(),
        'limitations': [
            'Only rankings are compared; score thresholds and abstention are not calibrated here.',
            'The same cases, annotation snapshot and ontology are used for every method.',
            'All gene/disease hypotheses are recomputed under each set of weights.',
            'without_phenotype still uses disease-level inheritance; molecular_only is the strict comparator.',
            'Ties use deterministic variant-ID ordering; equal-score group sizes are reported.',
            'No open negative controls; no claims about clinical accuracy or causal ablation effects.',
        ]}
    return output


def write_comparison(result: dict, directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'benchmark.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    methods = list(METHODS)
    with (directory / 'ranks.tsv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['case_id', 'expected_gene', *methods])
        for row in result['cases']:
            writer.writerow([row['case_id'], row['expected_gene_from_label'], *[row['ranks'][m] for m in methods]])
    def percentage(value):
        return f'{100 * value:.1f}%' if value is not None else 'нет данных'
    lines = ['# Сравнение экспертных модулей', '',
        'Открытый учебный development-набор. Это сравнение вариантов ранжирования, а не независимая клиническая точность.', '',
        '| Метод | Случаи | Top-1 | Top-3 | Top-5 | MRR |', '|---|---:|---:|---:|---:|---:|']
    for method, metric in result['metrics'].items():
        mrr = f'{metric["mrr"]:.4f}' if metric['mrr'] is not None else 'нет данных'
        lines.append(f'| {method} | {metric["cases"]} | {percentage(metric["top_1"])} | '
                     f'{percentage(metric["top_3"])} | {percentage(metric["top_5"])} | {mrr} |')
    lines.extend(['', '| Случай | Ожидаемый ген | ' + ' | '.join(methods) + ' |',
                  '|---|---|' + '---:|' * len(methods)])
    for row in result['cases']:
        lines.append('| ' + row['case_id'] + ' | ' + str(row['expected_gene_from_label']) + ' | '
                     + ' | '.join(str(row['ranks'][m]) for m in methods) + ' |')
    lines.extend(['', 'Веса, размеры групп равных баллов, источники и контрольные суммы находятся в JSON.',
        'molecular_only использует только молекулярные сведения; without_phenotype сохраняет канал наследования.',
        'Сортировка при равных баллах зависит от ID варианта; сравнение на шести случаях имеет ограниченную силу.',
        'Точки принятия решения и качество отказа этой командой не проверяются.', ''])
    (directory / 'benchmark.md').write_text('\n'.join(lines), encoding='utf-8')
