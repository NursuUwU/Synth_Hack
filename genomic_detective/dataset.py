"""Read cases only. Deliberately never reads development_labels.tsv or worked_example.md."""
from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path

from .models import Case, Variant


class Dataset:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.archive = zipfile.ZipFile(self.path) if self.path.is_file() else None
        if self.archive:
            matches = [name for name in self.archive.namelist()
                       if name.endswith('/objects.tsv') and not name.startswith('__MACOSX/')]
            if len(matches) != 1:
                raise ValueError('Archive must contain exactly one objects.tsv')
            self.prefix = matches[0][:-len('objects.tsv')]
        else:
            self.prefix = ''

    def close(self):
        if self.archive:
            self.archive.close()

    def read(self, relative: str) -> str:
        if '..' in Path(relative).parts or Path(relative).is_absolute():
            raise ValueError('Unsafe dataset path')
        if self.archive:
            return self.archive.read(self.prefix + relative).decode('utf-8-sig')
        return (self.path / relative).read_text(encoding='utf-8-sig')

    def objects(self) -> list[dict]:
        return list(csv.DictReader(io.StringIO(self.read('objects.tsv')), delimiter='\t'))

    def load_cases(self, split: str = 'all', case_ids: list[str] | None = None) -> list[Case]:
        objects = self.objects()
        selected = [obj for obj in objects if (split == 'all' or obj['split'] == split)
                    and (not case_ids or obj['case_id'] in case_ids)]
        if case_ids and set(case_ids) - {obj['case_id'] for obj in selected}:
            raise ValueError('Requested cases not found in selected split')
        return [self.load_case(obj) for obj in selected]

    def load_case(self, obj: dict) -> Case:
        cid = obj['case_id']
        prefix = f"{obj['split']}/{cid}"
        meta = json.loads(self.read(prefix + '_metadata.json'))
        phenotype = json.loads(self.read(prefix + '_phenotype.json'))
        if meta.get('case_id') != cid or phenotype.get('id') != cid:
            raise ValueError(f'{cid}: inconsistent case identifiers')
        rows = list(csv.DictReader(io.StringIO(self.read(prefix + '_variants.tsv')), delimiter='\t'))
        required = {'chrom', 'pos', 'ref', 'alt', 'genotype'}
        if not rows or not required.issubset(rows[0]):
            raise ValueError(f'{cid}: missing variant columns')
        variants = [Variant(row['chrom'], int(row['pos']), row['ref'], row['alt'],
                            row['genotype'], meta['genome_assembly']) for row in rows]
        vcf = [line.split('\t') for line in self.read(prefix + '.vcf').splitlines()
               if line and not line.startswith('#')]
        vcf_records = [(Variant(row[0], int(row[1]), row[3], row[4], row[9],
                                meta['genome_assembly']).id, row[9]) for row in vcf]
        if sorted(vcf_records) != sorted((v.id, v.genotype) for v in variants):
            raise ValueError(f'{cid}: VCF/TSV mismatch')
        features = phenotype['phenotypicFeatures']
        if len(variants) != int(obj['variant_count']) or len(features) != int(obj['phenotype_count']):
            raise ValueError(f'{cid}: counts do not match objects.tsv')
        if any(v.chrom != obj['chromosome'] for v in variants):
            raise ValueError(f'{cid}: chromosome mismatch')
        positive, excluded = {}, {}
        seen = set()
        for feature in features:
            term = feature['type']['id']
            if term in seen:
                raise ValueError(f'{cid}: repeated HPO term {term}')
            seen.add(term)
            target = excluded if feature.get('excluded', False) else positive
            target[term] = feature['type'].get('label', term)
        meta['split'] = obj['split']
        return Case(cid, variants, positive, excluded, meta, phenotype.get('sex', 'UNKNOWN_SEX'))
