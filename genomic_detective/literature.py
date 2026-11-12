"""Curated literature evidence, separate from VEP consequence prediction.

The registry is a small, reviewed subset. An exact-allele functional record may
raise the consequence component, but the VEP consequence itself is never
rewritten. Non-exact scopes are displayed as context only and never change a
score. Support values are heuristic, not ACMG codes, likelihood ratios or
diagnosis probabilities.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATH = ROOT / 'config/molecular_evidence.json'
ALLOWED_STATUS = 'approved_for_research_ranking'
EVIDENCE_SCOPES = frozenset({'exact_variant', 'same_site_other_allele',
                             'codon_cluster', 'gene_mechanism'})
SUPPORTING_SCOPES = frozenset({'exact_variant'})
VARIANT_PATTERN = re.compile(r'^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)):(\d+):([ACGT]):([ACGT])$')


class LiteratureError(ValueError):
    pass


def _variant_fields(value: str) -> tuple[str, int, str, str]:
    match = VARIANT_PATTERN.fullmatch(str(value))
    if not match:
        raise LiteratureError(f'Invalid curated variant key: {value!r}')
    chrom, position, ref, alt = match.groups()
    if int(position) < 1 or ref == alt:
        raise LiteratureError(f'Invalid curated variant key: {value!r}')
    return chrom, int(position), ref, alt


def _support(value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise LiteratureError('Invalid literature support')
    if not 0 <= value <= 1:
        raise LiteratureError('Literature support must be inside [0,1]')
    return float(value)


class LiteratureRegistry:
    """Validated, immutable-in-practice view of the curated literature file."""

    def __init__(self, payload: dict):
        if not isinstance(payload, dict) or payload.get('schema_version') != 1:
            raise LiteratureError('Unsupported literature evidence schema')
        self.version = str(payload.get('version') or '')
        self.status = str(payload.get('status') or '')
        self.reviewed_on = str(payload.get('reviewed_on') or '')
        self.warnings: list[str] = []
        self.entries: list[dict] = []
        if self.status != ALLOWED_STATUS:
            self.warnings.append('Literature registry status is not approved_for_research_ranking; no entries loaded')
        else:
            seen = set()
            for raw in payload.get('entries', []):
                entry = self._validate(raw)
                if entry is None:
                    continue
                if entry['id'] in seen:
                    self.warnings.append(f"Duplicate literature entry ignored: {entry['id']}")
                    continue
                seen.add(entry['id'])
                self.entries.append(entry)
        self.digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def _validate(self, raw) -> dict | None:
        if not isinstance(raw, dict):
            self.warnings.append('Non-object literature entry ignored')
            return None
        if raw.get('status', self.status) != ALLOWED_STATUS:
            self.warnings.append('Literature entry without research approval ignored')
            return None
        required = ('id', 'assembly', 'variant', 'gene', 'disease_ids', 'evidence_scope',
                    'evidence_type', 'fact', 'source', 'approved_by')
        missing = [name for name in required if not raw.get(name)]
        if missing:
            self.warnings.append('Literature entry without ' + ', '.join(missing) + ' ignored')
            return None
        if not re.fullmatch(r'[a-z0-9_]+', str(raw['id'])):
            self.warnings.append('Literature entry with unsafe id ignored')
            return None
        if raw['assembly'] != 'GRCh38':
            self.warnings.append(f"Literature entry {raw['id']} is not GRCh38; ignored")
            return None
        try:
            _variant_fields(raw['variant'])
            support = _support(raw.get('support', 0))
        except LiteratureError as error:
            self.warnings.append(f"Literature entry {raw['id']} ignored: {error}")
            return None
        if raw['evidence_scope'] not in EVIDENCE_SCOPES:
            self.warnings.append(f"Literature entry {raw['id']} has unknown evidence_scope; ignored")
            return None
        if not isinstance(raw['disease_ids'], list) or not all(
                isinstance(item, str) and item for item in raw['disease_ids']):
            self.warnings.append(f"Literature entry {raw['id']} has invalid disease_ids; ignored")
            return None
        source = raw['source']
        if not isinstance(source, dict) or not (source.get('pmid') or source.get('url')):
            self.warnings.append(f"Literature entry {raw['id']} has no source; ignored")
            return None
        if not isinstance(raw['approved_by'], str):
            self.warnings.append(f"Literature entry {raw['id']} has no approval; ignored")
            return None
        entry = dict(raw)
        entry['support'] = support if raw['evidence_scope'] in SUPPORTING_SCOPES else 0.0
        entry['disease_ids'] = list(dict.fromkeys(raw['disease_ids']))
        entry['not_established'] = list(raw.get('not_established') or [])
        return entry

    @classmethod
    def load(cls, path=DEFAULT_PATH):
        return cls(json.loads(Path(path).read_text(encoding='utf-8')))

    def match(self, variant, gene, disease_id) -> list[dict]:
        """Return registry entries applicable to this allele, gene and hypothesis."""
        if variant is None or not gene or not disease_id:
            return []
        if isinstance(variant, str):
            chrom, position, ref, alt = _variant_fields(variant)
            assembly = 'GRCh38'
            variant_key = variant
        else:
            chrom, position, ref, alt, assembly = variant.chrom, variant.pos, variant.ref, variant.alt, variant.assembly
            variant_key = variant.id
        matched = []
        for entry in self.entries:
            if entry['gene'] != gene or disease_id not in entry['disease_ids']:
                continue
            if entry['assembly'] != assembly:
                continue
            scope = entry['evidence_scope']
            if scope == 'exact_variant':
                if entry['variant'] == variant_key:
                    matched.append(entry)
            elif scope == 'same_site_other_allele':
                e_chrom, e_pos, e_ref, e_alt = _variant_fields(entry['variant'])
                if (e_chrom, e_pos) == (chrom, position) and e_alt != alt:
                    matched.append(entry)
            else:
                # codon_cluster and gene_mechanism are context only and cannot
                # be tied to this allele without additional coordinates.
                matched.append(entry)
        return matched

    def exact(self, variant, gene, disease_id) -> list[dict]:
        return [entry for entry in self.match(variant, gene, disease_id)
                if entry['evidence_scope'] in SUPPORTING_SCOPES]

    def sources(self, entries) -> set[str]:
        return {'LIT:' + entry['id'] for entry in entries}

    def metadata(self) -> dict:
        return {'version': self.version, 'status': self.status, 'reviewed_on': self.reviewed_on,
                'entries': len(self.entries), 'digest': self.digest,
                'scope': 'Curated research evidence; absence of a record is not benign evidence'}


@lru_cache(maxsize=1)
def default_registry() -> LiteratureRegistry:
    return LiteratureRegistry.load(DEFAULT_PATH)
