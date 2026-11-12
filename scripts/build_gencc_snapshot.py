"""Build a compact GenCC gene-disease validity snapshot for genes in the local VEP cache.

Input is the GenCC submissions TSV export (CC0, https://thegencc.org/), downloaded separately:
https://search.thegencc.org/download/action/submissions-export-tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://search.thegencc.org/download/action/submissions-export-tsv'


def cached_genes(cache: Path) -> set[str]:
    genes = set()
    with sqlite3.connect(cache) as connection:
        for (payload,) in connection.execute('SELECT payload FROM annotations'):
            raw = json.loads(payload).get('raw')
            for entry in raw if isinstance(raw, list) else [raw]:
                for transcript in (entry or {}).get('transcript_consequences', []):
                    if transcript.get('gene_symbol'):
                        genes.add(transcript['gene_symbol'])
    return genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'resources/gencc-snapshot.json')
    args = parser.parse_args()
    genes = cached_genes(ROOT / 'data/cache/vep.sqlite3')
    records = {}
    with args.tsv.open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['gene_symbol'] not in genes:
                continue
            ids = sorted({row['disease_curie'], row['disease_original_curie'], row['submitted_as_disease_id']} - {''})
            records.setdefault(row['gene_symbol'], []).append({
                'disease_ids': ids, 'disease': row['disease_title'],
                'classification': row['classification_title'], 'moi': row['moi_title'],
                'submitter': row['submitter_title'], 'date': row['submitted_as_date'][:10]})
    snapshot = {'schema_version': 1, 'source': 'GenCC submissions export', 'url': URL, 'license': 'CC0 1.0',
                'fetched_at': datetime.fromtimestamp(args.tsv.stat().st_mtime, timezone.utc).isoformat(),
                'input_sha256': hashlib.sha256(args.tsv.read_bytes()).hexdigest(),
                'scope': 'Genes present in the local VEP annotation cache',
                'records': dict(sorted(records.items()))}
    args.output.write_text(json.dumps(snapshot, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
    print(f'{len(records)} genes, {sum(map(len, records.values()))} submissions -> {args.output}')


if __name__ == '__main__':
    main()
