from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from genomic_detective.annotations import AnnotationCache, VepProvider, normalize
from genomic_detective.dataset import Dataset
from genomic_detective.experts import InheritanceExpert, MolecularExpert, clinical_support
from genomic_detective.guards import GuardError, check_presentation, check_ranking
from genomic_detective.knowledge import Knowledge
from genomic_detective.models import Case, Variant
from genomic_detective.pipeline import Pipeline, write_run
from genomic_detective.ranking import RankingConfig, rank
from genomic_detective.reporting import OpenAIReporter, html_report, make_presentation, template_selection


ROOT = Path(__file__).resolve().parent.parent


def knowledge():
    kb = Knowledge()
    kb.parents = {'HP:0000001': set(), 'HP:0000002': {'HP:0000001'},
                  'HP:0000003': {'HP:0000001'}}
    kb.labels = {'HP:0000001': 'root', 'HP:0000002': 'feature A', 'HP:0000003': 'feature B'}
    for gene, disease, term in [('GENEA', 'TEST:A', 'HP:0000002'), ('GENEB', 'TEST:B', 'HP:0000003')]:
        kb.diseases[disease] = {'id': disease, 'name': disease, 'positive': {term},
            'negative': set(), 'inheritance': {'AD'}, 'genes': {gene}, 'references': set()}
        kb.gene_diseases[gene].add(disease)
    kb.build_index()
    return kb


def annotation(variant, gene, consequence='missense_variant', af=0.00001):
    return {'variant_id': variant.id, 'assembly': 'GRCh38', 'status': 'annotated',
        'transcripts': [{'gene_symbol': gene, 'transcript_id': 'TEST:transcript:' + gene,
                         'biotype': 'protein_coding', 'canonical': 1,
                         'consequence_terms': [consequence]}],
        'gnomad_max_af': af, 'clinical': [], 'frequencies': [],
        'source': {'namespace': 'test-fixture-only'}}


class FakeProvider:
    def __init__(self, annotations):
        self.annotations = annotations

    def annotate(self, variants, **kwargs):
        return self.annotations, {'namespace': 'test-fixture-only'}


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.a = Variant('chr1', 100, 'A', 'T', '0/1')
        self.b = Variant('chr1', 200, 'G', 'C', '0/1')
        self.case = Case('fixture', [self.a, self.b], {'HP:0000002': 'feature A'}, {})
        self.kb = knowledge()
        self.config = RankingConfig.load(ROOT / 'config/ranking.json')
        self.annotations = {self.a.id: annotation(self.a, 'GENEA'), self.b.id: annotation(self.b, 'GENEB')}

    def pipeline(self, annotations=None, kb=None):
        return Pipeline(FakeProvider(annotations or self.annotations), kb or self.kb, self.config)

    def test_phenotype_change_changes_ranking(self):
        first = self.pipeline().run([self.case])[0]
        self.case.positive_hpo = {'HP:0000003': 'feature B'}
        second = self.pipeline().run([self.case])[0]
        self.assertEqual(first['candidates'][0]['gene'], 'GENEA')
        self.assertEqual(second['candidates'][0]['gene'], 'GENEB')
        self.assertNotEqual(first['input_sha256'], second['input_sha256'])

    def test_exclusion_reduces_support_but_does_not_veto(self):
        self.case.positive_hpo = {'HP:0000001': 'root', 'HP:0000003': 'feature B'}
        self.case.excluded_hpo = {'HP:0000002': 'feature A'}
        report = self.pipeline().run([self.case])[0]
        a = next(c for c in report['candidates'] if c['gene'] == 'GENEA')
        self.assertEqual(len(a['experts']['phenotype']['conflicts']), 1)
        self.assertEqual(len(report['candidates']), 2)

    def test_missing_frequency_is_not_zero_frequency(self):
        expert = MolecularExpert()
        missing = expert.assess(self.a, annotation(self.a, 'GENEA', af=None), 'GENEA')
        rare = expert.assess(self.a, annotation(self.a, 'GENEA', af=0), 'GENEA')
        self.assertEqual(missing['components']['frequency'], 0)
        self.assertGreater(rare['score'], missing['score'])
        self.assertTrue(any('unavailable' in entry for entry in missing['missing']))

    def test_vus_and_conflicts_add_no_pathogenic_support(self):
        cases = [(['uncertain_significance'], 'uncertain_significance'),
                 (['pathogenic', 'benign'], 'conflicting'),
                 (['conflicting interpretations of pathogenicity'], 'conflicting'),
                 (['risk_factor'], 'not_actionable'), (['not pathogenic'], 'not_actionable')]
        for labels, expected_status in cases:
            with self.subTest(labels=labels):
                self.assertEqual(clinical_support(labels), (expected_status, 0.0))
                data = annotation(self.a, 'GENEA')
                data['clinical'] = [{'significance': labels}]
                result = MolecularExpert().assess(self.a, data, 'GENEA')
                self.assertEqual(result['components']['clinical'], 0.0)
                self.assertEqual(result['clinical_status'], expected_status)
                self.assertTrue(any(ev['criterion'] == 'clinical_annotation' for ev in result['evidence']))

    def test_recognized_pathogenic_labels_are_separate_from_vus(self):
        self.assertEqual(clinical_support(['pathogenic/likely pathogenic']), ('pathogenic', 0.8))
        self.assertEqual(clinical_support(['likely_pathogenic']), ('likely_pathogenic', 0.8))
        self.assertEqual(clinical_support(['likely_benign']), ('likely_benign', 0.0))

    def test_unphased_allele_order_is_equivalent(self):
        expert = InheritanceExpert(self.kb)
        reverse = Variant('chr1', 100, 'A', 'T', '1/0')
        self.assertEqual(expert.assess(self.a, 'TEST:A')['score'], expert.assess(reverse, 'TEST:A')['score'])

    def test_single_heterozygote_does_not_confirm_recessive_model(self):
        self.kb.diseases['TEST:A']['inheritance'] = {'AR'}
        expert = InheritanceExpert(self.kb)
        hetero = expert.assess(self.a, 'TEST:A')
        homo = expert.assess(Variant('chr1', 100, 'A', 'T', '1/1'), 'TEST:A')
        self.assertLess(hetero['score'], homo['score'])
        self.assertIn('не подтверждает', hetero['evidence'][0]['text'])

    def test_disease_hypotheses_combined_before_selection(self):
        self.kb.diseases['TEST:Z'] = {'id': 'TEST:Z', 'name': 'TEST:Z',
            'positive': {'HP:0000002'}, 'negative': set(), 'inheritance': set(),
            'genes': {'GENEA'}, 'references': set()}
        self.kb.gene_diseases['GENEA'].add('TEST:Z')
        self.kb.build_index()
        report = self.pipeline().run([self.case])[0]
        best = report['candidates'][0]
        self.assertEqual(best['disease_id'], 'TEST:A')
        self.assertIn('TEST:Z', [item['disease_id'] for item in best['other_disease_hypotheses']])

    def test_missing_annotations_withhold_prediction(self):
        self.annotations[self.b.id] = {'variant_id': self.b.id, 'status': 'unavailable',
            'assembly': 'GRCh38', 'transcripts': [], 'gnomad_max_af': None}
        report = self.pipeline().run([self.case])[0]
        self.assertEqual(report['status'], 'incomplete')
        self.assertIsNone(report['prediction'])
        self.assertNotEqual(report['prediction'], 'no_causal_variant_in_vcf')

    def test_missing_knowledge_does_not_generate_negative_prediction(self):
        report = self.pipeline(kb=Knowledge()).run([self.case])[0]
        self.assertEqual(report['status'], 'incomplete')
        self.assertIsNone(report['prediction'])

    def test_uninformative_root_phenotype_withholds_prediction(self):
        self.case.positive_hpo = {'HP:0000001': 'root'}
        report = self.pipeline().run([self.case])[0]
        self.assertEqual(report['status'], 'incomplete')
        self.assertIsNone(report['prediction'])

    def test_full_data_weak_evidence_can_abstain(self):
        for variant, gene in [(self.a, 'GENEA'), (self.b, 'GENEB')]:
            self.annotations[variant.id] = annotation(variant, gene, 'intron_variant', 0.1)
        report = self.pipeline().run([self.case])[0]
        self.assertEqual(report['status'], 'no_convincing_candidate')
        self.assertEqual(report['prediction'], 'no_causal_variant_in_vcf')

    def test_guard_rejects_altered_score_and_non_input_variant(self):
        candidates = rank(self.case, self.annotations, self.kb, self.config)
        bad = copy.deepcopy(candidates)
        bad[0]['score'] += 0.01
        with self.assertRaises(GuardError):
            check_ranking(self.case, bad, self.annotations, self.config)
        bad = copy.deepcopy(candidates)
        bad[0]['variant_id'] = 'chr1:999:A:T'
        with self.assertRaises(GuardError):
            check_ranking(self.case, bad, self.annotations, self.config)

    def test_guard_rejects_llm_rank_change_and_unknown_evidence(self):
        report = self.pipeline().run([self.case])[0]
        selection = template_selection(report)
        selection['candidates'].reverse()
        with self.assertRaises(GuardError):
            check_presentation(selection, report)

    def test_guard_rejects_malformed_presentation(self):
        report = self.pipeline().run([self.case])[0]
        for bad in [None, [], {'case_id': report['case_id'], 'status': report['status'], 'candidates': [None]}]:
            with self.assertRaises(GuardError):
                check_presentation(bad, report)
        selection = template_selection(report)
        selection['candidates'][0]['evidence_ids'].append('fabricated_source')
        with self.assertRaises(GuardError):
            check_presentation(selection, report)

    def test_llm_invalid_output_falls_back(self):
        report = self.pipeline().run([self.case])[0]

        class BadReporter:
            def select(self, data):
                result = template_selection(data)
                result['candidates'][0]['variant_id'] = 'invented'
                return result

        presentation = make_presentation(report, True, BadReporter())
        self.assertEqual(presentation['provider'], 'template')
        self.assertIn('GuardError', presentation['warning'])

    def test_missing_llm_credentials_falls_back_without_network(self):
        report = self.pipeline().run([self.case])[0]
        with patch.dict('os.environ', {}, clear=True):
            presentation = make_presentation(report, True)
        self.assertEqual(presentation['provider'], 'template')
        self.assertEqual(presentation['guard']['status'], 'passed')

    def test_llm_adapter_requests_only_constrained_fact_selection(self):
        report = self.pipeline().run([self.case])[0]
        expected = template_selection(report)

        class FakeHttp:
            def json(inner, url, payload, headers):
                self.assertEqual(url, 'https://api.openai.com/v1/responses')
                self.assertFalse(payload['store'])
                self.assertTrue(payload['text']['format']['strict'])
                self.assertEqual(payload['model'], 'fixture-model')
                self.assertNotIn('expected_answer', payload['input'])
                self.assertNotIn('fixture-key', payload['input'])
                return {'status': 'completed', 'output': [{'type': 'message', 'content': [
                    {'type': 'output_text', 'text': json.dumps(expected)}]}]}

        with patch.dict('os.environ', {'OPENAI_API_KEY': 'fixture-key', 'GENOMIC_LLM_MODEL': 'fixture-model'}):
            result = make_presentation(report, True, OpenAIReporter(FakeHttp()))
        self.assertEqual(result['provider'], 'openai_evidence_selection')
        self.assertEqual(result['selection'], expected)

    def test_external_html_is_escaped(self):
        report = self.pipeline().run([self.case])[0]
        report['candidates'][0]['disease_name'] = '<script>alert(1)</script>'
        result = html_report(report)
        self.assertNotIn('<script>', result)
        self.assertIn('&lt;script&gt;', result)

    def test_partial_export_never_leaves_stale_complete_prediction(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            report = self.pipeline().run([self.case])[0]
            write_run([report], directory)
            self.assertTrue((directory / 'predictions.tsv').exists())
            report['prediction'], report['status'] = None, 'incomplete'
            write_run([report], directory)
            self.assertFalse((directory / 'predictions.tsv').exists())
            self.assertTrue((directory / 'predictions.partial.tsv').exists())

    def test_open_labels_are_not_read_during_loading(self):
        dataset = Dataset(ROOT / 'Project_3_genomic_detecive_for_students.zip')
        original = dataset.read

        def guarded_read(name):
            if 'labels' in name or 'worked_example' in name:
                raise AssertionError('Label leakage')
            return original(name)

        try:
            with patch.object(dataset, 'read', guarded_read):
                cases = dataset.load_cases()
            self.assertEqual(len(cases), 10)
            self.assertEqual(sum(len(case.variants) for case in cases), 4096)
            self.assertEqual(sum(len(case.positive_hpo) for case in cases), 179)
        finally:
            dataset.close()


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant('chr1', 100, 'A', 'T', '0/1')
        self.raw = {'assembly_name': 'GRCh38', 'seq_region_name': '1', 'start': 100, 'end': 100,
            'allele_string': 'A/T', 'input': self.variant.vcf_line(),
            'most_severe_consequence': 'missense_variant',
            'transcript_consequences': [{'gene_symbol': 'GENEA', 'variant_allele': 'T',
                                        'consequence_terms': ['missense_variant']}],
            'colocated_variants': [{'id': 'test', 'start': 100, 'end': 100, 'strand': 1,
                'allele_string': 'A/T', 'clin_sig': ['pathogenic'],
                'frequencies': {'A': {'gnomadg': 0.8}, 'T': {'gnomadg': 0.002}}}]}
        self.source = {'namespace': 'fixture'}

    def test_frequency_uses_alternate_allele(self):
        result = normalize(self.raw, self.variant, self.source, 'fixture')
        self.assertEqual(result['gnomad_max_af'], 0.002)
        self.assertEqual(len(result['clinical']), 1)

    def test_multiallelic_clinical_assertion_not_attributed_blindly(self):
        self.raw['colocated_variants'][0]['allele_string'] = 'A/T/C'
        result = normalize(self.raw, self.variant, self.source, 'fixture')
        self.assertEqual(result['clinical'], [])
        self.assertTrue(result['ambiguous_clinical_skipped'])

    def test_multiallelic_allele_specific_field_separates_classifications(self):
        entry = self.raw['colocated_variants'][0]
        entry['allele_string'] = 'A/T/C'
        entry['clin_sig'] = ['pathogenic', 'benign']
        entry['clin_sig_allele'] = 'T:pathogenic;C:benign;T:risk_factor'
        result = normalize(self.raw, self.variant, self.source, 'fixture')
        self.assertEqual(result['clinical'][0]['significance'], ['pathogenic', 'risk_factor'])
        self.assertNotIn('benign', result['clinical'][0]['significance'])

    def test_overlapping_mane_select_precedes_downstream_mane_plus(self):
        data = annotation(self.variant, 'GENEA', 'synonymous_variant')
        data['transcripts'][0]['mane_select'] = 'TEST:select'
        data['transcripts'].append({'gene_symbol': 'GENEA', 'transcript_id': 'TEST:downstream',
            'biotype': 'protein_coding', 'mane_plus_clinical': 'TEST:plus',
            'consequence_terms': ['downstream_gene_variant']})
        result = MolecularExpert().assess(self.variant, data, 'GENEA')
        self.assertNotEqual(result['selected_transcript'], 'TEST:downstream')
        self.assertEqual(result['consequences'], ['synonymous_variant'])

    def test_wrong_assembly_or_alleles_rejected(self):
        self.raw['assembly_name'] = 'GRCh37'
        with self.assertRaises(ValueError):
            normalize(self.raw, self.variant, self.source, 'fixture')
        self.raw['assembly_name'] = 'GRCh38'
        self.raw['allele_string'] = 'G/T'
        with self.assertRaises(ValueError):
            normalize(self.raw, self.variant, self.source, 'fixture')

    def test_offline_cache_is_release_specific_and_does_not_call_network(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = AnnotationCache(Path(temp) / 'cache.sqlite3')
            try:
                cache.set_metadata(self.source)
                cache.put(self.variant, 'old-release', {'raw': self.raw, 'fetched_at': 'fixture'})

                class NoNetwork:
                    def json(self, *args, **kwargs):
                        raise AssertionError('Offline made network request')

                provider = VepProvider(cache, NoNetwork())
                result, _ = provider.annotate([self.variant], offline=True)
                self.assertEqual(result[self.variant.id]['status'], 'unavailable')
                cache.put(self.variant, 'fixture', {'raw': self.raw, 'fetched_at': 'fixture'})
                result, _ = provider.annotate([self.variant], offline=True)
                self.assertEqual(result[self.variant.id]['status'], 'annotated')
            finally:
                cache.close()


if __name__ == '__main__':
    unittest.main()
