"""Package the frozen external pilot for offline reproduction; maintainer only."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = ('manifest.json', 'inputs.json', 'labels.json', 'selection-audit.json',
         'annotations.json', 'clinical-snapshot.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'outputs/validation_external_v1')
    args = parser.parse_args()
    archive = ROOT / 'resources/validation-pilot-v1.zip'
    sources = [args.source / name for name in FILES]
    sources.extend(sorted(args.source.glob('clinvar-*.xml')))
    files = {}
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as target:
        for source in sources:
            data = source.read_bytes()
            target.writestr(source.name, data)
            files[source.name] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    manifest = {'schema_version': 1, 'archive': archive.name, 'bytes': archive.stat().st_size,
        'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'files': files,
        'scope': 'Frozen public external semi-synthetic pilot; inputs and annotations, not precomputed rankings',
        'license': 'phenopacket-store-LICENSE', 'phenopacket_release': '0.1.27'}
    (ROOT / 'resources/validation-pilot-v1.manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(f'Bundled {len(files)} files, {manifest["bytes"]} bytes')


if __name__ == '__main__':
    main()
