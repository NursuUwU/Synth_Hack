"""Maintainer-only public ClinVar refresh; uses all VEP clinical hits, never labels."""
import hashlib
import json
import sys
import time
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from genomic_detective.annotations import AnnotationCache, VepProvider
from genomic_detective.clinvar import parse_vcv
from genomic_detective.dataset import Dataset
from genomic_detective.models import utc_now
from genomic_detective.network import HttpClient


def main():
    dataset = Dataset(ROOT / 'Project_3_genomic_detecive_for_students.zip')
    cache = AnnotationCache(ROOT / 'data/cache/vep.sqlite3')
    try:
        cases = dataset.load_cases()
        annotations, _ = VepProvider(cache).annotate([v for c in cases for v in c.variants], offline=True)
    finally:
        dataset.close()
        cache.close()
    requested = {key: value for key, value in annotations.items() if value.get('clinical')}
    rsids = sorted({c['record_id'] for a in requested.values() for c in a['clinical']
                    if c['record_id'] and c['record_id'].startswith('rs') and c['record_id'][2:].isdigit()})
    client = HttpClient(timeout=60)
    raw_dir = ROOT / 'data/cache/clinvar/requests'
    raw_dir.mkdir(parents=True, exist_ok=True)
    found, raw_sources = {}, []
    for start in range(0, len(rsids), 20):
        group = rsids[start:start + 20]
        query = urllib.parse.urlencode({'db': 'clinvar', 'term': ' OR '.join(group), 'retmode': 'json', 'retmax': 500})
        base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        search = client.json(base + 'esearch.fcgi?' + query)['esearchresult']
        if int(search['count']) > len(search['idlist']):
            raise ValueError('ClinVar search truncated; use smaller batches')
        time.sleep(0.4)
        if not search['idlist']:
            continue
        url = base + 'efetch.fcgi?' + urllib.parse.urlencode({
            'db': 'clinvar', 'rettype': 'vcv', 'is_variationid': 'true', 'id': ','.join(search['idlist'])})
        payload = client.request(url, headers={'Accept': 'application/xml'})
        records = parse_vcv(payload)
        if not records:
            raise ValueError('ClinVar returned no classified simple alleles')
        sha = hashlib.sha256(payload).hexdigest()
        (raw_dir / (sha + '.xml')).write_bytes(payload)
        raw_sources.append({'url': url, 'sha256': sha, 'bytes': len(payload)})
        for record in records:
            if record['variant_id'] in requested:
                found.setdefault(record['variant_id'], {})[record['accession']] = record
        print(f'ClinVar: {min(start + 20, len(rsids))}/{len(rsids)} rsIDs, {len(found)} exact allele matches', flush=True)
        time.sleep(0.4)
    result = {'schema_version': 1, 'fetched_at': utc_now(), 'assembly': 'GRCh38',
        'scope': 'All variants with an allele-specific VEP clinical annotation; no phenotype or label-based selection',
        'requested_variants': sorted(requested), 'requested_rsids': rsids,
        'unmatched_variants': sorted(set(requested) - set(found)),
        'sources': raw_sources,
        'records': {key: list(value.values()) for key, value in sorted(found.items())}}
    target = ROOT / 'resources/clinvar-snapshot.json'
    temporary = target.with_suffix('.json.part')
    temporary.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    temporary.replace(target)
    print(f'Saved {len(found)}/{len(requested)} variants to {target.name}', flush=True)


if __name__ == '__main__':
    main()
