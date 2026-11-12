"""Read-only structural audit of the supplied teaching dataset.

Does not annotate variants, infer diagnoses, or use labels to rank candidates.
Run with the bundled Python runtime from the workspace root.
"""

import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from itertools import combinations
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parent.parent
ARCHIVE = WORKSPACE / 'Project_3_genomic_detecive_for_students.zip'
ROOT = 'Project_3_genomic_detecive_for_students/'


def key(row):
    return f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"


def main():
    problems = []
    cases = []
    variant_sets = {}
    annotation_missing = Counter()
    annotation_fields = ['gene', 'consequence', 'population_AF', 'clinical_significance']
    with zipfile.ZipFile(ARCHIVE) as archive:
        def read(name):
            return archive.read(ROOT + name).decode('utf-8-sig')

        def table(name):
            return list(csv.DictReader(io.StringIO(read(name)), delimiter='\t'))

        objects = table('objects.tsv')
        labels = {row['case_id']: row for row in table('development_labels.tsv')}
        for obj in objects:
            cid, split = obj['case_id'], obj['split']
            prefix = f'{split}/{cid}'
            rows = table(prefix + '_variants.tsv')
            meta = json.loads(read(prefix + '_metadata.json'))
            pheno = json.loads(read(prefix + '_phenotype.json'))
            features = pheno['phenotypicFeatures']
            vcf = [line.split('\t') for line in read(prefix + '.vcf').splitlines()
                   if line and not line.startswith('#')]
            ids = [key(row) for row in rows]
            variant_sets[cid] = set(ids)
            tsv_records = Counter((key(row), row['genotype']) for row in rows)
            vcf_records = Counter((f'{row[0]}:{row[1]}:{row[3]}:{row[4]}', row[9]) for row in vcf)
            positive = [f['type']['id'] for f in features if not f.get('excluded', False)]
            negative = [f['type']['id'] for f in features if f.get('excluded', False)]
            invalid_rows = [i + 2 for i, row in enumerate(rows)
                            if None in row or any(value is None for value in row.values())
                            or not row['pos'].isdigit() or int(row['pos']) < 1
                            or not re.fullmatch('[ACGT]', row['ref'])
                            or not re.fullmatch('[ACGT]', row['alt'])
                            or row['ref'] == row['alt']]
            missing = {field: sum(row[field] in ('NA', '', '.') for row in rows)
                       for field in annotation_fields}
            annotation_missing.update(missing)
            label = labels.get(cid)
            result = {
                'case_id': cid,
                'split': split,
                'chromosomes': sorted({row['chrom'] for row in rows}),
                'variants': len(rows),
                'unique_variants': len(set(ids)),
                'duplicate_variants': len(ids) - len(set(ids)),
                'variant_types': dict(Counter(row['variant_type'] for row in rows)),
                'genotypes': dict(Counter(row['genotype'] for row in rows)),
                'vcf_records': len(vcf),
                'vcf_tsv_match_including_genotype': tsv_records == vcf_records,
                'vcf_unavailable_fields': {
                    name: sum(row[index] == '.' for row in vcf)
                    for name, index in [('QUAL', 5), ('FILTER', 6), ('INFO', 7)]
                },
                'annotation_missing': missing,
                'positive_hpo': len(positive),
                'excluded_hpo': len(negative),
                'features_with_onset': sum('onset' in f for f in features),
                'duplicate_hpo_ids': len(positive + negative) - len(set(positive + negative)),
                'contradictory_hpo_ids': sorted(set(positive) & set(negative)),
                'invalid_hpo_ids': [term for term in positive + negative
                                    if not re.fullmatch(r'HP:\d{7}', term)],
                'invalid_variant_rows': invalid_rows,
                'sex': pheno.get('sex'),
                'phenotype_truncated': meta['phenotype_truncated'],
                'metadata': meta,
                'phenotype_top_level_fields': sorted(pheno),
                'variant_columns': list(rows[0]),
                'object_counts_match': len(rows) == int(obj['variant_count'])
                                       and len(features) == int(obj['phenotype_count']),
                'case_ids_match': meta.get('case_id') == cid and pheno.get('id') == cid,
                'object_chromosome_matches': {row['chrom'] for row in rows} == {obj['chromosome']},
                'label': label,
                'label_variant_present': label['causal_variant'] in set(ids) if label else None,
                'label_genotypes': [row['genotype'] for row in rows
                                    if label and key(row) == label['causal_variant']],
            }
            for field in ['vcf_tsv_match_including_genotype', 'object_counts_match',
                          'case_ids_match', 'object_chromosome_matches']:
                if not result[field]:
                    problems.append(f'{cid}: failed {field}')
            for field in ['duplicate_variants', 'duplicate_hpo_ids', 'contradictory_hpo_ids',
                          'invalid_hpo_ids', 'invalid_variant_rows']:
                if result[field]:
                    problems.append(f'{cid}: {field} = {result[field]}')
            if label and not result['label_variant_present']:
                problems.append(f'{cid}: supplied label variant absent')
            cases.append(result)

    overlaps = [{'a': a, 'b': b, 'shared_variants': len(variant_sets[a] & variant_sets[b]),
                 'identical_variant_sets': variant_sets[a] == variant_sets[b]}
                for a, b in combinations(variant_sets, 2) if variant_sets[a] & variant_sets[b]]
    report = {
        'archive': ARCHIVE.name,
        'archive_sha256': hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
        'method': 'Structural audit only; no external annotation or biological validation.',
        'case_count': len(cases),
        'total_variant_rows': sum(case['variants'] for case in cases),
        'unique_variant_keys_across_cases': len(set().union(*variant_sets.values())),
        'total_positive_hpo': sum(case['positive_hpo'] for case in cases),
        'total_excluded_hpo': sum(case['excluded_hpo'] for case in cases),
        'vcf_unavailable_fields': {
            field: sum(case['vcf_unavailable_fields'][field] for case in cases)
            for field in ['QUAL', 'FILTER', 'INFO']
        },
        'development_label_genotypes': {
            case['case_id']: case['label_genotypes'] for case in cases if case['label']
        },
        'development_label_presence_in_same_chromosome_cases': {
            case['case_id']: {
                other['case_id']: case['label']['causal_variant'] in variant_sets[other['case_id']]
                for other in cases if other['chromosomes'] == case['chromosomes']
            }
            for case in cases if case['label']
        },
        'annotation_missing': dict(annotation_missing),
        'structural_problems': problems,
        'cases': cases,
        'overlaps': overlaps,
    }
    (Path(__file__).parent / 'dataset_audit.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'cases'}, ensure_ascii=False, indent=2))
    print('CASE\tSPLIT\tN\tHPO+\tHPO-\tTRUNCATED\tGT\tLABEL_PRESENT')
    for case in cases:
        print('\t'.join(str(case[field]) for field in
                        ['case_id', 'split', 'variants', 'positive_hpo', 'excluded_hpo',
                         'phenotype_truncated', 'genotypes', 'label_variant_present']))


if __name__ == '__main__':
    main()
