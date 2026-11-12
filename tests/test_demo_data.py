import hashlib
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

from genomic_detective.demo_data import MEMBERS, digest, ensure_demo_data


class DemoDataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        resources = self.root / 'resources'
        resources.mkdir()
        self.payloads = {name: ('fixture ' + name).encode() for name in MEMBERS}
        archive = resources / 'demo-data.zip'
        with zipfile.ZipFile(archive, 'w') as output:
            for name, payload in self.payloads.items():
                output.writestr(name, payload)
        self.manifest = {'schema_version': 1, 'bundle_sha256': digest(archive),
                         'files': {name: {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                                   for name, data in self.payloads.items()}}
        self.write_manifest()

    def write_manifest(self):
        (self.root / 'resources/demo-data.manifest.json').write_text(json.dumps(self.manifest), encoding='utf-8')

    def test_clean_repository_restores_all_data_offline_and_is_idempotent(self):
        self.assertEqual(set(ensure_demo_data(self.root)), MEMBERS)
        for name, payload in self.payloads.items():
            self.assertEqual((self.root / 'data' / name).read_bytes(), payload)
        self.assertEqual(ensure_demo_data(self.root), [])

    def test_complete_custom_data_is_preserved(self):
        ensure_demo_data(self.root)
        target = self.root / 'data/knowledge/hp.obo'
        target.write_bytes(b'custom version')
        self.assertEqual(ensure_demo_data(self.root), [])
        self.assertEqual(target.read_bytes(), b'custom version')

    def test_partial_custom_knowledge_is_not_silently_mixed(self):
        target = self.root / 'data/knowledge/hp.obo'
        target.parent.mkdir(parents=True)
        target.write_bytes(b'custom version')
        with self.assertRaisesRegex(ValueError, 'Partial custom HPO'):
            ensure_demo_data(self.root)
        self.assertFalse((self.root / 'data/cache/vep.sqlite3').exists())

    def test_partial_matching_snapshot_is_completed(self):
        target = self.root / 'data/knowledge/hp.obo'
        target.parent.mkdir(parents=True)
        target.write_bytes(self.payloads['knowledge/hp.obo'])
        self.assertEqual(len(ensure_demo_data(self.root)), 4)

    def test_corrupt_archive_is_rejected_before_installation(self):
        (self.root / 'resources/demo-data.zip').write_bytes(b'bad download')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            ensure_demo_data(self.root)
        self.assertFalse((self.root / 'data').exists())

    def test_invalid_member_digest_installs_nothing(self):
        self.manifest['files']['knowledge/hp.obo']['sha256'] = '0' * 64
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, 'member checksum'):
            ensure_demo_data(self.root)
        self.assertFalse((self.root / 'data/cache/vep.sqlite3').exists())

    def test_unlisted_archive_members_are_rejected(self):
        archive = self.root / 'resources/demo-data.zip'
        with zipfile.ZipFile(archive, 'a') as output:
            output.writestr('../outside.txt', b'no')
        self.manifest['bundle_sha256'] = digest(archive)
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, 'Unexpected files'):
            ensure_demo_data(self.root)
        self.assertFalse((self.root / 'outside.txt').exists())


class DistributionIntegrityTests(unittest.TestCase):
    def test_real_demo_archive_matches_manifest(self):
        root = Path(__file__).resolve().parents[1]
        manifest = json.loads((root / 'resources/demo-data.manifest.json').read_text(encoding='utf-8'))
        archive = root / 'resources/demo-data.zip'
        self.assertEqual(digest(archive), manifest['bundle_sha256'])
        self.assertEqual(manifest['cached_variants'], 3690)
        self.assertEqual(digest(root / 'Project_3_genomic_detecive_for_students.zip'), manifest['dataset_sha256'])
        with zipfile.ZipFile(archive) as source:
            self.assertEqual(set(source.namelist()), MEMBERS)

    def test_official_windows_runtime_is_intact_and_includes_license(self):
        root = Path(__file__).resolve().parents[1]
        settings = json.loads((root / 'config/python-runtime.json').read_text(encoding='utf-8'))
        name = settings['windows']['amd64']['url'].rsplit('/', 1)[-1]
        archive = root / 'resources' / name
        self.assertEqual(digest(archive), settings['windows']['amd64']['sha256'])
        with zipfile.ZipFile(archive) as source:
            self.assertIn('python.exe', source.namelist())
            self.assertTrue(any('license' in name.casefold() for name in source.namelist()))


if __name__ == '__main__':
    unittest.main()
