"""Maintainer tool: package public reference data and the synthetic variant cache."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from genomic_detective.demo_data import MEMBERS, digest


def main():
    resources = ROOT / 'resources'
    resources.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as temporary:
        snapshot = Path(temporary) / 'vep.sqlite3'
        original = sqlite3.connect((ROOT / 'data/cache/vep.sqlite3').as_uri() + '?mode=ro', uri=True)
        backup = sqlite3.connect(snapshot)
        try:
            original.backup(backup)
            count = backup.execute('SELECT COUNT(*) FROM annotations').fetchone()[0]
            metadata = json.loads(backup.execute("SELECT value FROM metadata WHERE key='active'").fetchone()[0])
        finally:
            backup.close()
            original.close()
        files = {}
        archive = resources / 'demo-data.zip'
        with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
            for name in sorted(MEMBERS):
                source = snapshot if name == 'cache/vep.sqlite3' else ROOT / 'data' / name
                payload = source.read_bytes()
                entry = zipfile.ZipInfo(name, date_time=(2026, 11, 12, 0, 0, 0))
                entry.compress_type = zipfile.ZIP_DEFLATED
                entry.external_attr = 0o100644 << 16
                output.writestr(entry, payload)
                files[name] = {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest()}
        manifest = {'schema_version': 1, 'snapshot_date': '2026-11-12',
                    'bundle_sha256': digest(archive), 'bundle_bytes': archive.stat().st_size,
                    'hpo_version': '2026-09-01', 'vep': metadata, 'cached_variants': count,
                    'dataset_sha256': digest(ROOT / 'Project_3_genomic_detecive_for_students.zip'),
                    'files': files, 'attribution': 'See resources/README.md and docs/data_sources.md'}
        (resources / 'demo-data.manifest.json').write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(f'Packaged {count} cached variants; {archive.stat().st_size:,} bytes.')


if __name__ == '__main__':
    main()
