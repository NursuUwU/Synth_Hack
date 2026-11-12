"""Local HPO ontology and disease-specific gene/phenotype associations."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from .models import utc_now
from .network import HttpClient


BASE = 'https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/'
FILES = ('hp.obo', 'genes_to_phenotype.txt', 'phenotype.hpoa')
INHERITANCE = {
    'HP:0000006': 'AD', 'HP:0000007': 'AR', 'HP:0001417': 'XLD',
    'HP:0001419': 'XLR', 'HP:0001423': 'XL', 'HP:0001427': 'MITO',
}


def fetch_knowledge(directory: Path, client: HttpClient | None = None, progress=None):
    directory.mkdir(parents=True, exist_ok=True)
    client = client or HttpClient(timeout=60)
    manifest = {'downloaded_at': utc_now(), 'files': {}}
    for name in FILES:
        payload = client.request(BASE + name, headers={'Accept': '*/*'})
        if len(payload) < 1000 or payload[:200].lstrip().lower().startswith(b'<!doctype html'):
            raise ValueError(f'Invalid knowledge file: {name}')
        temporary = directory / (name + '.part')
        temporary.write_bytes(payload)
        temporary.replace(directory / name)
        manifest['files'][name] = {'url': BASE + name, 'bytes': len(payload),
                                    'sha256': hashlib.sha256(payload).hexdigest()}
        if progress:
            progress(f'HPO: {name} ({len(payload):,} bytes)')
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return manifest


def read_tsv(path: Path) -> list[dict]:
    lines = []
    with path.open(encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('#database_id\t'):
                lines.append(line[1:])
            elif line.strip() and not line.startswith('#'):
                lines.append(line)
    return list(csv.DictReader(io.StringIO(''.join(lines)), delimiter='\t'))


class Knowledge:
    def __init__(self):
        self.parents: dict[str, set[str]] = {}
        self.labels: dict[str, str] = {}
        self.aliases: dict[str, str] = {}
        self.diseases: dict[str, dict] = {}
        self.gene_diseases: dict[str, set[str]] = defaultdict(set)
        self.term_disease_count = defaultdict(int)
        self.manifest = {}
        self.warnings = []

    @classmethod
    def load(cls, directory: Path):
        kb = cls()
        if not all((directory / name).exists() for name in FILES):
            kb.warnings.append('HPO knowledge files are missing; phenotype/inheritance evidence unavailable')
            return kb
        manifest_path = directory / 'manifest.json'
        if manifest_path.exists():
            kb.manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
            for name, entry in kb.manifest.get('files', {}).items():
                if Path(name).name != name or name not in FILES:
                    raise ValueError('Unsafe knowledge manifest entry')
                if hashlib.sha256((directory / name).read_bytes()).hexdigest() != entry['sha256']:
                    raise ValueError(f'Knowledge checksum mismatch: {name}')
        else:
            kb.warnings.append('HPO knowledge manifest missing; source version is not recorded')
        kb._read_ontology(directory / 'hp.obo')
        rows = read_tsv(directory / 'genes_to_phenotype.txt')
        if not rows or not {'gene_symbol', 'hpo_id', 'disease_id'}.issubset(rows[0]):
            raise ValueError('Unsupported genes_to_phenotype schema')
        for row in rows:
            gene, disease_id, term = row['gene_symbol'], row['disease_id'], row['hpo_id']
            disease = kb.diseases.setdefault(disease_id, {
                'id': disease_id, 'name': disease_id, 'positive': set(), 'negative': set(),
                'inheritance': set(), 'genes': set(), 'references': set(),
            })
            disease['genes'].add(gene)
            kb.gene_diseases[gene].add(disease_id)
            if term in INHERITANCE:
                disease['inheritance'].add(INHERITANCE[term])
            else:
                disease['positive'].add(kb.resolve(term))
        for row in read_tsv(directory / 'phenotype.hpoa'):
            disease = kb.diseases.get(row.get('database_id'))
            if not disease:
                continue
            disease['name'] = row.get('disease_name') or disease['id']
            term = row.get('hpo_id')
            if not term:
                continue
            negated = row.get('qualifier') == 'NOT'
            if term in INHERITANCE:
                if not negated:
                    disease['inheritance'].add(INHERITANCE[term])
            elif negated:
                disease['negative'].add(kb.resolve(term))
                disease['positive'].discard(kb.resolve(term))
            if row.get('reference'):
                disease['references'].add(row['reference'])
        kb.build_index()
        return kb

    def _read_ontology(self, path: Path):
        current = None
        version = None
        with path.open(encoding='utf-8') as stream:
            for line in stream:
                line = line.strip()
                if line.startswith('data-version: '):
                    version = line[14:]
                if line.startswith('['):
                    current = None
                elif line.startswith('id: HP:'):
                    current = line[4:]
                    self.parents[current] = set()
                elif current and line.startswith('name: '):
                    self.labels[current] = line[6:]
                elif current and line.startswith('is_a: '):
                    self.parents[current].add(line[6:].split()[0])
                elif current and line.startswith('alt_id: '):
                    self.aliases[line[8:]] = current
        self.manifest['ontology_version'] = version

    def resolve(self, term: str) -> str:
        return self.aliases.get(term, term)

    @lru_cache(maxsize=40000)
    def ancestors(self, term: str) -> frozenset[str]:
        term = self.resolve(term)
        # Iterative traversal also tolerates an accidental ontology cycle.
        seen, pending = set(), [term]
        while pending:
            node = pending.pop()
            if node not in seen:
                seen.add(node)
                pending.extend(self.parents.get(node, set()))
        return frozenset(seen)

    def build_index(self):
        self.term_disease_count.clear()
        for disease in self.diseases.values():
            expanded = set()
            for term in disease['positive']:
                expanded.update(self.ancestors(term))
            disease['expanded'] = expanded
            for term in expanded:
                self.term_disease_count[term] += 1

    def ic(self, term: str) -> float:
        return -math.log((self.term_disease_count.get(self.resolve(term), 0) + 1)
                         / (len(self.diseases) + 1)) if self.diseases else 0.0
