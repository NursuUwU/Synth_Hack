import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from genomic_detective.cli import parser
from genomic_detective.evaluation import METHODS, compare_experts, write_comparison
from genomic_detective.knowledge import Knowledge
from genomic_detective.models import Case, Variant
from genomic_detective.ranking import RankingConfig, rank


ROOT = Path(__file__).resolve().parent.parent


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.config = RankingConfig.load(ROOT / 'config/ranking.json')
        self.kb = Knowledge()
        self.kb.parents = {'HP:0000001': set(), 'HP:0000002': {'HP:0000001'}, 'HP:0000003': {'HP:0000001'}}
        self.kb.labels = {term: term for term in self.kb.parents}
        self.a = Variant('chr1', 100, 'A', 'T', '0/1')
        self.b = Variant('chr1', 200, 'G', 'C', '0/1')
        self.annotations = {}
        for variant, gene, disease, term in [(self.a, 'TESTA', 'TEST:A', 'HP:0000002'),
                                             (self.b, 'TESTB', 'TEST:B', 'HP:0000003')]:
            self.kb.diseases[disease] = {'id': disease, 'name': disease, 'positive': {term},
                'negative': set(), 'inheritance': {'AD'}, 'genes': {gene}, 'references': set()}
            self.kb.gene_diseases[gene].add(disease)
            self.annotations[variant.id] = {'status': 'annotated', 'gnomad_max_af': 0.00001, 'clinical': [],
                'transcripts': [{'gene_symbol': gene, 'transcript_id': 'TEST:' + gene,
                    'biotype': 'protein_coding', 'canonical': 1, 'consequence_terms': ['missense_variant']}]}
        self.kb.build_index()
        self.cases = [Case('fixture_A', [self.a, self.b], {'HP:0000002': 'A'}, {}, {'split': 'development'}),
                      Case('fixture_B', [self.a, self.b], {'HP:0000003': 'B'}, {}, {'split': 'development'})]
        self.labels = [{'case_id': 'fixture_A', 'expected_answer': self.a.id, 'causal_gene': 'TESTA'},
                       {'case_id': 'fixture_B', 'expected_answer': self.b.id, 'causal_gene': 'TESTB'}]

    def compare(self):
        with patch('genomic_detective.evaluation.read_development_labels', return_value=self.labels):
            return compare_experts(self.cases, self.annotations, self.kb, self.config, Path('fixture.zip'))

    def test_comparison_recomputes_rankings_before_reading_labels(self):
        with patch('genomic_detective.ranking.rank', wraps=rank) as ranker:
            def read_labels(_):
                self.assertEqual(ranker.call_count, len(METHODS) * len(self.cases))
                return self.labels
            with patch('genomic_detective.evaluation.read_development_labels', side_effect=read_labels):
                result = compare_experts(self.cases, self.annotations, self.kb, self.config, Path('fixture.zip'))
        self.assertEqual(result['metrics']['combined']['top_1'], 1.0)
        self.assertEqual(result['metrics']['molecular_only']['top_1'], 0.5)
        self.assertEqual(result['cases'][0]['equal_score_group_size']['molecular_only'], 2)

    def test_config_and_cases_are_not_mutated(self):
        weights = dict(self.config.weights)
        cases = copy.deepcopy(self.cases)
        self.compare()
        self.assertEqual(self.config.weights, weights)
        self.assertEqual(self.cases, cases)

    def test_incomplete_cases_are_not_counted_as_successes(self):
        self.annotations[self.b.id]['status'] = 'unavailable'
        result = self.compare()
        self.assertEqual(len(result['excluded_cases']), 2)
        self.assertIsNone(result['metrics']['combined']['top_1'])
        self.assertEqual(result['metrics']['combined']['cases'], 0)

    def test_original_development_only(self):
        self.cases[0].metadata['phenotype_override'] = 'fixture.json'
        with self.assertRaises(ValueError):
            self.compare()
        self.cases[0].metadata = {'split': 'hidden'}
        with self.assertRaises(ValueError):
            self.compare()

    def test_exports_and_default_split(self):
        result = self.compare()
        with tempfile.TemporaryDirectory() as temp:
            write_comparison(result, Path(temp))
            self.assertTrue((Path(temp) / 'benchmark.json').exists())
            self.assertIn('molecular_only', (Path(temp) / 'benchmark.md').read_text(encoding='utf-8'))
            self.assertIn('fixture_B', (Path(temp) / 'ranks.tsv').read_text(encoding='utf-8'))
        self.assertEqual(parser().parse_args(['benchmark']).split, 'development')


if __name__ == '__main__':
    unittest.main()
