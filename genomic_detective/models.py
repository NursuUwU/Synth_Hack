from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def chromosome(value: str) -> str:
    value = str(value).removeprefix('chr')
    if value == 'M':
        value = 'MT'
    if not re.fullmatch(r'(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)', value):
        raise ValueError(f'Unsupported chromosome: {value!r}')
    return 'chr' + value


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str
    assembly: str = 'GRCh38'

    def __post_init__(self):
        object.__setattr__(self, 'chrom', chromosome(self.chrom))
        if self.assembly != 'GRCh38':
            raise ValueError('The first implementation supports GRCh38 only')
        if self.pos < 1 or self.ref not in 'ACGT' or self.alt not in 'ACGT':
            raise ValueError(f'Invalid SNV: {self}')
        if len(self.ref) != 1 or len(self.alt) != 1 or self.ref == self.alt:
            raise ValueError('Only single alternate SNVs are supported')
        if self.genotype not in {'0/1', '1/0', '1/1', '0|1', '1|0', '1|1', './.', '.|.'}:
            raise ValueError(f'Unsupported genotype: {self.genotype}')

    @property
    def id(self) -> str:
        return f'{self.chrom}:{self.pos}:{self.ref}:{self.alt}'

    @property
    def cache_key(self) -> str:
        return f'{self.assembly}:{self.id}'

    @property
    def heterozygous(self) -> bool:
        return self.genotype in {'0/1', '1/0', '0|1', '1|0'}

    @property
    def homozygous_alt(self) -> bool:
        return self.genotype in {'1/1', '1|1'}

    def vcf_line(self) -> str:
        return f'{self.chrom[3:]} {self.pos} . {self.ref} {self.alt} . . .'


@dataclass
class Case:
    id: str
    variants: list[Variant]
    positive_hpo: dict[str, str]
    excluded_hpo: dict[str, str]
    metadata: dict = field(default_factory=dict)
    sex: str = 'UNKNOWN_SEX'

    def __post_init__(self):
        if not re.fullmatch(r'[A-Za-z0-9_-]+', self.id):
            raise ValueError('Case id must be a safe file identifier')
        if not self.variants:
            raise ValueError(f'{self.id}: no variants')
        if len({v.id for v in self.variants}) != len(self.variants):
            raise ValueError(f'{self.id}: duplicate variants')
        if set(self.positive_hpo) & set(self.excluded_hpo):
            raise ValueError(f'{self.id}: contradictory HPO features')
        for term in self.positive_hpo | self.excluded_hpo:
            if not re.fullmatch(r'HP:\d{7}', term):
                raise ValueError(f'Invalid HPO id: {term}')


# Disease-context molecular modes; v3 adds the v1.3 rarity, ClinVar and MOI rules.
CONTEXT_MODES = frozenset({'context_v2', 'context_v3'})


def finite_score(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
