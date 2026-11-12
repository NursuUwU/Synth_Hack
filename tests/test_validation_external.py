import copy
import hashlib
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from validate_external import digest, metrics, restore_bundle, unpack_case, verify_reference


def packet():
    return {'phenotypicFeatures': [{'type': {'id': 'HP:0001250', 'label': 'Seizure'}}],
        'metaData': {'externalReferences': [{'id': 'PMID:123456'}]},
        'interpretations': [{'progressStatus': 'SOLVED', 'diagnosis': {'genomicInterpretations': [
            {'interpretationStatus': 'CAUSATIVE', 'variantInterpretation': {'variationDescriptor': {
                'vcfRecord': {'genomeAssembly': 'hg38', 'chrom': 'chr1', 'pos': 100, 'ref': 'A', 'alt': 'G'},
                'geneContext': {'symbol': 'GENEA'}, 'allelicState': {'id': 'GENO:0000135'}}}}]}}]}


class ExternalTests(unittest.TestCase):
    def test_replay_detects_rank_and_config_regression(self):
        result = {'cases': [{'model': 'v1.2', 'case_id': 'fixture', 'expected_rank': 1,
                            'status': 'candidate', 'prediction': 'chr1:100:A:G',
                            'annotation_coverage': 1, 'config_sha256': 'original'}]}
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / 'reference.json'
            reference.write_text(json.dumps(result), encoding='utf-8')
            verify_reference(result, reference)
            for field, replacement in [('expected_rank', 2), ('config_sha256', 'changed')]:
                altered = copy.deepcopy(result)
                altered['cases'][0][field] = replacement
                with self.assertRaises(ValueError):
                    verify_reference(altered, reference)

    def test_frozen_bundle_checksums_and_cohort_independence(self):
        root = Path(__file__).resolve().parents[1]
        info = json.loads((root / 'resources/validation-pilot-v1.manifest.json').read_text())
        archive = root / 'resources' / info['archive']
        self.assertEqual(hashlib.sha256(archive.read_bytes()).hexdigest(), info['sha256'])
        with zipfile.ZipFile(archive) as source:
            for name, entry in info['files'].items():
                self.assertEqual(hashlib.sha256(source.read(name)).hexdigest(), entry['sha256'])
            manifest = json.loads(source.read('manifest.json'))
            inputs = json.loads(source.read('inputs.json'))
            labels = json.loads(source.read('labels.json'))
        self.assertEqual(digest(inputs), manifest['inputs_sha256'])
        self.assertEqual(digest(labels), manifest['labels_sha256'])
        self.assertEqual(len(inputs), 20)
        self.assertEqual(len({d['expected_gene'] for d in labels}), 20)
        pmids = [p for d in labels for p in d['pmids']]
        self.assertEqual(len(set(pmids)), len(pmids))
        self.assertFalse(set(manifest['excluded_pmids']) & set(pmids))
        self.assertFalse(set(manifest['excluded_development_genes']) & {d['expected_gene'] for d in labels})
        for case, label in zip(inputs, labels):
            self.assertEqual(len(case['variants']), 400)
            self.assertNotIn('expected_variant', case['metadata'])
            self.assertEqual(sum('chr' + v['chrom'].removeprefix('chr') + ':' + str(v['pos']) + ':' + v['ref'] + ':' + v['alt'] == label['expected_variant'] for v in case['variants']), 1)

    def test_bundle_restoration_is_verified_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'frozen'
            restore_bundle(output)
            self.assertTrue((output / 'clinical-snapshot.json').is_file())
            with self.assertRaises(FileExistsError):
                restore_bundle(output)

    def test_explicit_grch38_single_snv_and_zygosity(self):
        variant, gene, pmids, positive, excluded = unpack_case(packet())
        self.assertEqual(variant.id, 'chr1:100:A:G')
        self.assertEqual(variant.genotype, '0/1')
        self.assertEqual(gene, 'GENEA')
        self.assertEqual(pmids, ['PMID:123456'])
        self.assertEqual(set(positive), {'HP:0001250'})
        self.assertEqual(excluded, {})

    def test_no_implicit_liftover_or_zygosity_guess(self):
        for key, value in [('genomeAssembly', 'hg19'), ('alt', 'GG'), ('chrom', 'chrX')]:
            data = packet()
            data['interpretations'][0]['diagnosis']['genomicInterpretations'][0]['variantInterpretation']['variationDescriptor']['vcfRecord'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                unpack_case(data)
        data = packet()
        data['interpretations'][0]['diagnosis']['genomicInterpretations'][0]['variantInterpretation']['variationDescriptor']['allelicState'] = {}
        with self.assertRaises(ValueError):
            unpack_case(data)

    def test_multiallele_and_contradictory_phenotype_excluded(self):
        data = packet()
        genomic = data['interpretations'][0]['diagnosis']['genomicInterpretations']
        genomic.append(copy.deepcopy(genomic[0]))
        with self.assertRaises(ValueError):
            unpack_case(data)
        data = packet()
        data['phenotypicFeatures'].append({'type': {'id': 'HP:0001250'}, 'excluded': True})
        with self.assertRaises(ValueError):
            unpack_case(data)

    def test_unfinished_cases_remain_in_denominator(self):
        rows = [{'status': status, 'expected_rank': rank, 'case_id': str(i)}
                for i, (status, rank) in enumerate([('candidate', 1), ('candidate', 10), ('incomplete', 1), ('error', None)])]
        result = metrics(rows)
        self.assertEqual(result['attempted'], 4)
        self.assertEqual(result['completed'], 2)
        self.assertEqual(result['top_1_all_attempted'], .25)
        self.assertEqual(result['candidate_wrong_first'], 1)
        self.assertEqual(result['errors'], ['3'])


if __name__ == '__main__':
    unittest.main()
