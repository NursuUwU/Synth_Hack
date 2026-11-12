import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from validate_controls import FrozenProvider, execute, summarize, transform
from genomic_detective.models import Case, Variant
from genomic_detective.ranking import RankingConfig
from test_workflow import annotation, knowledge


class ControlTests(unittest.TestCase):
    def setUp(self):
        self.kb = knowledge()
        self.a = Variant('1', 100, 'A', 'T', '0/1')
        self.b = Variant('1', 200, 'G', 'C', '0/1')
        self.case = Case('fixture', [self.a, self.b], {'HP:0000002': 'feature A'}, {})
        self.annotations = {self.a.id: annotation(self.a, 'GENEA'), self.b.id: annotation(self.b, 'GENEB')}

    def test_targeted_deletion_and_masks_never_modify_original(self):
        original = copy.deepcopy((self.case, self.annotations))
        removed, anns, _ = transform(self.case, self.annotations, self.a.id, 'remove_known_causal', self.kb)
        self.assertEqual([v.id for v in removed.variants], [self.b.id])
        self.assertEqual(set(anns), {self.b.id})
        _, anns, _ = transform(self.case, self.annotations, self.a.id, 'frequency_unavailable', self.kb)
        self.assertIsNone(anns[self.a.id]['gnomad_max_af'])
        self.assertEqual((self.case, self.annotations), original)

    def test_masked_annotation_is_not_a_benign_annotation(self):
        case, anns, _ = transform(self.case, self.annotations, self.a.id, 'causal_annotation_unavailable', self.kb)
        self.assertEqual(anns[self.a.id]['status'], 'unavailable')
        self.assertIn(self.a, case.variants)
        self.assertEqual(anns[self.b.id], self.annotations[self.b.id])

    def test_half_profile_is_deterministic_and_preserves_exclusions(self):
        self.case.positive_hpo.update({'HP:0000003': 'feature B', 'HP:0000001': 'root'})
        first, _, _ = transform(self.case, self.annotations, self.a.id, 'half_positive_hpo', self.kb)
        second, _, _ = transform(self.case, self.annotations, self.a.id, 'half_positive_hpo', self.kb)
        self.assertEqual(first, second)
        self.assertEqual(len(first.positive_hpo), 2)
        self.assertEqual(first.excluded_hpo, self.case.excluded_hpo)

    def test_all_attempted_denominator_counts_failures_and_no_rank_for_deletion(self):
        rows = [{'model': 'v1.2', 'scenario': 'baseline', 'status': status, 'expected_rank': 1}
                for status in ['candidate', 'incomplete', 'error']]
        rows.append({'model': 'v1.2', 'scenario': 'remove_known_causal', 'status': 'candidate', 'expected_rank': None})
        groups = summarize(rows)
        self.assertEqual(groups[0]['attempted'], 3)
        self.assertEqual(groups[0]['top_1_all_attempted'], 1 / 3)
        self.assertIsNone(groups[1]['top_1_all_attempted'])
        self.assertEqual(groups[1]['candidate_wrong_first'], 1)

    def test_frozen_provider_copies_and_refuses_network(self):
        provider = FrozenProvider(self.annotations, {})
        anns, _ = provider.annotate([self.a], offline=True)
        anns[self.a.id]['transcripts'].clear()
        self.assertTrue(self.annotations[self.a.id]['transcripts'])
        with self.assertRaises(ValueError):
            provider.annotate([self.a], offline=False)

    def test_real_pipeline_controls_withhold_on_missing_annotation_or_hpo(self):
        config = RankingConfig.load(Path(__file__).resolve().parents[1] / 'config/ranking.json')
        from genomic_detective.pipeline import Pipeline
        for scenario in ['causal_annotation_unavailable', 'empty_positive_hpo']:
            case, anns, _ = transform(self.case, self.annotations, self.a.id, scenario, self.kb)
            report = Pipeline(FrozenProvider(anns, {}), self.kb, config).run([case])[0]
            self.assertEqual(report['status'], 'incomplete')
            self.assertIsNone(report['prediction'])


if __name__ == '__main__':
    unittest.main()
