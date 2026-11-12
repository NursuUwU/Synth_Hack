"""Install the versioned offline demonstration snapshot shipped with the repository."""
from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import zipfile
from pathlib import Path


MEMBERS = {'cache/vep.sqlite3', 'knowledge/hp.obo', 'knowledge/genes_to_phenotype.txt',
           'knowledge/phenotype.hpoa', 'knowledge/manifest.json'}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def ensure_demo_data(root: Path, progress=None):
    """Restore missing files only; never mix a partial custom HPO with our snapshot."""
    root = root.resolve()
    destinations = {name: (root / 'data' / name).resolve() for name in MEMBERS}
    if any(not path.is_relative_to(root / 'data') for path in destinations.values()):
        raise ValueError('Demo data destination escapes the project data directory')
    missing = {name for name, path in destinations.items() if not path.is_file()}
    if not missing:
        return []
    resources = root / 'resources'
    manifest = json.loads((resources / 'demo-data.manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema_version') != 1 or set(manifest['files']) != MEMBERS:
        raise ValueError('Unsupported demo snapshot manifest')
    archive = resources / 'demo-data.zip'
    if digest(archive) != manifest['bundle_sha256']:
        raise ValueError('Demo snapshot checksum mismatch. Download the repository again.')
    if any(name.startswith('knowledge/') for name in missing):
        for name in MEMBERS - missing:
            if name.startswith('knowledge/') and digest(destinations[name]) != manifest['files'][name]['sha256']:
                raise ValueError('Partial custom HPO differs from the bundled snapshot. '
                                 'Complete it with prepare or select another --knowledge directory.')
    if progress:
        progress('Preparing bundled offline data (HPO + Ensembl VEP); no downloads...')
    data_dir = root / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.demo-staging-', dir=data_dir) as temporary:
        staged = {}
        with zipfile.ZipFile(archive) as source:
            if set(source.namelist()) != MEMBERS or len(source.infolist()) != len(MEMBERS):
                raise ValueError('Unexpected files in demo snapshot')
            for number, name in enumerate(sorted(missing)):
                entry = source.getinfo(name)
                expected = manifest['files'][name]
                if entry.file_size != expected['bytes'] or entry.file_size > 150_000_000:
                    raise ValueError('Invalid demo snapshot member size')
                target = Path(temporary) / str(number)
                with source.open(name) as incoming, target.open('wb') as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                if digest(target) != expected['sha256']:
                    raise ValueError('Demo snapshot member checksum mismatch: ' + name)
                staged[name] = target
        # Verify every staged file before installing any of them.
        restored = []
        for name, target in staged.items():
            destination = destinations[name]
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                target.replace(destination)
                restored.append(name)
    if progress:
        progress(f'Offline data ready ({len(restored)} files restored).')
    return restored
