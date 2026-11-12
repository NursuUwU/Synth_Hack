"""Offline phenotype sensitivity sweep for case_01 (not a production module).

Each positive HPO term is moved to excluded one at a time. The script writes a
TSV with the resulting top hypothesis and FBN1 rank. It demonstrates ranking
sensitivity to phenotype; it does not create a new clinical case.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from genomic_detective.annotations import AnnotationCache, VepProvider
from genomic_detective.dataset import Dataset
from genomic_detective.knowledge import Knowledge
from genomic_detective.pipeline import Pipeline
from genomic_detective.ranking import RankingConfig


def parse_args():
    parser = argparse.ArgumentParser(description='case_01 phenotype sensitivity sweep')
    parser.add_argument('--dataset', type=Path,
                        default=Path('Project_3_genomic_detecive_for_students.zip'))
    parser.add_argument('--cache', type=Path, default=Path('data/cache/vep.sqlite3'))
    parser.add_argument('--knowledge', type=Path, default=Path('data/knowledge'))
    parser.add_argument('--config', type=Path, default=Path('config/ranking.json'))
    parser.add_argument('--case', default='case_01')
    parser.add_argument('--gene', default='FBN1')
    parser.add_argument('--output', type=Path, default=Path('outputs/sensitivity_case01.tsv'))
    return parser.parse_args()


def rank_map(report):
    return {candidate['variant_id']: candidate['rank'] for candidate in report['candidates']}


def changed_positions(baseline, override):
    keys = set(baseline) | set(override)
    return sum(1 for key in keys if baseline.get(key) != override.get(key))


def main():
    args = parse_args()
    dataset = Dataset(args.dataset)
    cache = AnnotationCache(args.cache)
    try:
        case = next(item for item in dataset.load_cases('development', [args.case]))
        knowledge = Knowledge.load(args.knowledge)
        config = RankingConfig.load(args.config)
        provider = VepProvider(cache)
        pipeline = Pipeline(provider, knowledge, config)
        baseline = pipeline.run([copy.deepcopy(case)], offline=True)[0]
        baseline_ranks = rank_map(baseline)
        baseline_gene = baseline['candidates'][0]['gene']
        rows = []
        for term, label in sorted(case.positive_hpo.items()):
            modified = copy.deepcopy(case)
            modified.positive_hpo.pop(term)
            modified.excluded_hpo[term] = label
            modified.metadata['phenotype_override'] = f'sensitivity_case01:{term}'
            modified.__post_init__()
            report = pipeline.run([modified], offline=True)[0]
            top = report['candidates'][0]
            fbn1 = next((c for c in report['candidates'] if c['gene'] == args.gene), None)
            rows.append({
                'term_id': term,
                'label': label,
                'top1_gene': top['gene'],
                'top1_score': round(top['score'], 6),
                'FBN1_rank': fbn1['rank'] if fbn1 else '',
                'FBN1_score': round(fbn1['score'], 6) if fbn1 else '',
                'n_changed_positions': changed_positions(baseline_ranks, rank_map(report)),
                'leader_changed': top['gene'] != baseline_gene,
            })
            print(f"{term} {label}: top={top['gene']} FBN1_rank={fbn1['rank'] if fbn1 else '-'} "
                  f"changed={rows[-1]['n_changed_positions']}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        print(json.dumps({'baseline_top1_gene': baseline_gene,
                          'baseline_top1_score': round(baseline['candidates'][0]['score'], 6),
                          'terms': len(rows), 'output': str(args.output)}, ensure_ascii=False))
        return 0
    finally:
        cache.close()
        dataset.close()


if __name__ == '__main__':
    raise SystemExit(main())
