"""Versioned VEP raw cache and conservative, allele-specific normalization."""
from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path

from .models import Variant, chromosome, utc_now
from .network import HttpClient, NetworkError


SERVER = 'https://rest.ensembl.org'
OPTIONS = 'canonical=1&mane=1&hgvs=1&numbers=1&transcript_version=1'
NORMALIZER_VERSION = 1


class AnnotationCache:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path)
        self.db.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)')
        self.db.execute('CREATE TABLE IF NOT EXISTS annotations '
                        '(namespace TEXT, variant_key TEXT, payload TEXT NOT NULL, '
                        'PRIMARY KEY(namespace, variant_key))')

    def metadata(self) -> dict | None:
        row = self.db.execute("SELECT value FROM metadata WHERE key='active'").fetchone()
        return json.loads(row[0]) if row else None

    def set_metadata(self, data: dict):
        with self.db:
            self.db.execute("INSERT OR REPLACE INTO metadata VALUES ('active', ?)",
                            (json.dumps(data),))

    def get(self, variant: Variant, namespace: str) -> dict | None:
        row = self.db.execute('SELECT payload FROM annotations WHERE namespace=? AND variant_key=?',
                              (namespace, variant.cache_key)).fetchone()
        return json.loads(row[0]) if row else None

    def put(self, variant: Variant, namespace: str, record: dict):
        with self.db:
            self.db.execute('INSERT OR REPLACE INTO annotations VALUES (?, ?, ?)',
                            (namespace, variant.cache_key, json.dumps(record)))

    def close(self):
        self.db.close()


def normalize(raw: dict, variant: Variant, source: dict, fetched_at: str) -> dict:
    if raw.get('assembly_name') != variant.assembly:
        raise ValueError('VEP assembly mismatch')
    if chromosome(raw.get('seq_region_name', variant.chrom)) != variant.chrom:
        raise ValueError('VEP chromosome mismatch')
    if int(raw.get('start', 0)) != variant.pos or int(raw.get('end', 0)) != variant.pos:
        raise ValueError('VEP coordinate mismatch')
    if raw.get('allele_string') != f'{variant.ref}/{variant.alt}':
        raise ValueError('VEP reference/alternate allele mismatch')
    transcripts = []
    for entry in raw.get('transcript_consequences', []):
        if entry.get('variant_allele') != variant.alt:
            continue
        transcripts.append({key: entry.get(key) for key in (
            'gene_id', 'gene_symbol', 'transcript_id', 'biotype', 'impact', 'consequence_terms',
            'canonical', 'mane_select', 'mane_plus_clinical', 'hgvsc', 'hgvsp')})
    clinical = []
    frequencies = []
    ambiguous_clinical = False
    for entry in raw.get('colocated_variants', []):
        strand = entry.get('strand', 1)
        alt = variant.alt if strand == 1 else variant.alt.translate(str.maketrans('ACGT', 'TGCA'))
        alleles = entry.get('allele_string', '').split('/')
        if strand == -1:
            alleles = [a.translate(str.maketrans('ACGT', 'TGCA')) for a in alleles]
        coordinate_match = (int(entry.get('start', -1)) == variant.pos
                            and int(entry.get('end', -1)) == variant.pos)
        allele_match = variant.ref in alleles and variant.alt in alleles
        if not coordinate_match or not allele_match:
            continue
        for name, value in entry.get('frequencies', {}).get(alt, {}).items():
            if isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value) and 0 <= value <= 1:
                frequencies.append({'dataset': name, 'af': float(value), 'record_id': entry.get('id')})
        expected_alleles = {variant.ref, variant.alt}
        exact = len(alleles) == 2 and set(alleles) == expected_alleles
        significance = entry.get('clin_sig', [])
        if isinstance(significance, str):
            significance = [significance]
        allele_specific = entry.get('clin_sig_allele')
        selected_significance = []
        if isinstance(allele_specific, str):
            for assertion in allele_specific.split(';'):
                allele, separator, classification = assertion.partition(':')
                if separator and allele == alt and classification:
                    selected_significance.extend(classification.split(','))
        elif not allele_specific and exact:
            selected_significance = significance
        if selected_significance:
            clinical.append({'record_id': entry.get('id'), 'significance': sorted(set(selected_significance)),
                             'allele_specific_field': bool(allele_specific),
                             'review_status': None, 'condition': None,
                             'source': 'Ensembl co-located clinical annotation',
                             'caveat': 'Condition and ClinVar review status not independently retrieved'})
        elif significance:
            ambiguous_clinical = True
    gnomad = [f['af'] for f in frequencies if f['dataset'].lower().startswith('gnomad')]
    return {
        'schema_version': NORMALIZER_VERSION,
        'variant_id': variant.id, 'assembly': variant.assembly, 'status': 'annotated',
        'fetched_at': fetched_at, 'source': source,
        'transcripts': transcripts,
        'intergenic': not transcripts,
        'most_severe_consequence': raw.get('most_severe_consequence'),
        'frequencies': frequencies,
        'gnomad_max_af': max(gnomad) if gnomad else None,
        'clinical': clinical,
        'ambiguous_clinical_skipped': ambiguous_clinical,
    }


class VepProvider:
    def __init__(self, cache: AnnotationCache, client: HttpClient | None = None):
        self.cache = cache
        self.client = client or HttpClient()

    def metadata(self, offline: bool) -> dict:
        if offline:
            return self.cache.metadata() or {'namespace': 'empty', 'assembly': 'GRCh38', 'offline': True}
        assembly = self.client.json(SERVER + '/info/assembly/homo_sapiens?content-type=application/json')
        if assembly.get('default_coord_system_version') != 'GRCh38':
            raise ValueError('Remote default assembly is not GRCh38; refusing to annotate')
        releases = self.client.json(SERVER + '/info/data?content-type=application/json').get('releases', [])
        if not releases:
            raise ValueError('Ensembl release could not be determined')
        data = {'namespace': f'vep-{max(releases)}-GRCh38-n{NORMALIZER_VERSION}-{OPTIONS}',
                'assembly': 'GRCh38', 'assembly_patch': assembly.get('assembly_name'),
                'ensembl_releases': releases,
                'endpoint': SERVER + '/vep/homo_sapiens/region',
                'options': OPTIONS, 'verified_at': utc_now()}
        self.cache.set_metadata(data)
        return data

    def annotate(self, variants: list[Variant], offline: bool = False,
                 batch_size: int = 100, progress=None) -> tuple[dict, dict]:
        if not 1 <= batch_size <= 200:
            raise ValueError('VEP batch size must be 1..200')
        unique = {variant.cache_key: variant for variant in variants}
        source = self.metadata(offline)
        annotations, pending = {}, []
        for variant in unique.values():
            cached = self.cache.get(variant, source['namespace'])
            if cached:
                annotations[variant.id] = normalize(cached['raw'], variant, source, cached['fetched_at'])
            else:
                pending.append(variant)
        if not offline:
            for start in range(0, len(pending), batch_size):
                batch = pending[start:start + batch_size]
                records = self._batch(batch, source)
                annotations.update(records)
                if progress:
                    progress(f'Annotation: {min(start + batch_size, len(pending))}/{len(pending)} uncached variants')
        for variant in pending:
            if variant.id not in annotations:
                annotations[variant.id] = {
                    'variant_id': variant.id, 'assembly': variant.assembly, 'status': 'unavailable',
                    'reason': 'not_in_cache' if offline else 'annotation_request_failed',
                    'source': source, 'transcripts': [], 'clinical': [], 'frequencies': [],
                    'gnomad_max_af': None,
                }
        return annotations, source

    def _batch(self, batch: list[Variant], source: dict) -> dict:
        try:
            response = self.client.json(source['endpoint'] + '?' + OPTIONS,
                                        {'variants': [v.vcf_line() for v in batch]})
        except NetworkError as error:
            # A malformed record must not poison the whole batch. Do not fan out on outages.
            if 'HTTP 400' in str(error) and len(batch) > 1:
                middle = len(batch) // 2
                return self._batch(batch[:middle], source) | self._batch(batch[middle:], source)
            return {}
        if not isinstance(response, list):
            return {}
        expected = {v.vcf_line(): v for v in batch}
        records = {}
        for raw in response:
            if not isinstance(raw, dict):
                continue
            variant = expected.get(' '.join(raw.get('input', '').split()))
            if not variant or 'error' in raw:
                continue
            fetched_at = utc_now()
            try:
                record = normalize(raw, variant, source, fetched_at)
            except (ValueError, TypeError):
                continue
            self.cache.put(variant, source['namespace'], {'raw': raw, 'fetched_at': fetched_at})
            records[variant.id] = record
        return records
