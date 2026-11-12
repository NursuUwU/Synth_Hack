"""v1.3 rules: MOI compatibility, gnomAD absence and gene-level ClinVar support."""
import json
import unittest
from dataclasses import replace
from pathlib import Path

from genomic_detective.experts import InheritanceExpert, MolecularExpert
from genomic_detective.guards import check_ranking
from genomic_detective.models import Case, Variant
from genomic_detective.molecular import MolecularKnowledge
from genomic_detective.ranking import RankingConfig, decide, rank
from test_molecular_context import assertion
from test_workflow import annotation, knowledge

ROOT = Path(__file__).resolve().parents[1]


def recessive_knowledge():
    kb = knowledge()
    for disease in kb.diseases.values():
        disease['inheritance'] = {'AR'}
    return kb


class InheritanceV13Tests(unittest.TestCase):
    def setUp(self):
        self.dominant = InheritanceExpert(knowledge())
        self.recessive = InheritanceExpert(recessive_knowledge())
        self.het = Variant('chr1', 100, 'A', 'T', '0/1')
        self.hom = Variant('chr1', 100, 'A', 'T', '1/1')

    def test_single_recessive_heterozygote_is_filtered(self):
        result = self.recessive.assess(self.het, 'TEST:A', True, None, False)
        self.assertFalse(result['moi_compatible'])
        self.assertEqual(result['score'], 0.0)

    def test_second_rare_allele_allows_possible_compound_heterozygosity(self):
        result = self.recessive.assess(self.het, 'TEST:A', True, 0.005, True)
        self.assertTrue(result['moi_compatible'])
        self.assertEqual(result['score'], 0.6)
        self.assertIn('фаза', result['evidence'][0]['text'])

    def test_mode_specific_frequency_limits(self):
        self.assertFalse(self.dominant.assess(self.het, 'TEST:A', True, 0.002)['moi_compatible'])
        self.assertTrue(self.dominant.assess(self.het, 'TEST:A', True, 0.0005)['moi_compatible'])
        self.assertFalse(self.recessive.assess(self.hom, 'TEST:A', True, 0.005)['moi_compatible'])
        self.assertFalse(self.recessive.assess(self.het, 'TEST:A', True, 0.03, True)['moi_compatible'])

    def test_clinvar_whitelist_bypasses_frequency_but_not_genotype(self):
        dominant = self.dominant.assess(self.het, 'TEST:A', True, 0.005, False, True)
        self.assertTrue(dominant['moi_compatible'])
        self.assertIn('whitelist', dominant['evidence'][0]['text'])
        lone = self.recessive.assess(self.het, 'TEST:A', True, 0.005, False, True)
        self.assertFalse(lone['moi_compatible'])

    def test_unknown_frequency_is_not_treated_as_common(self):
        result = self.dominant.assess(self.het, 'TEST:A', True, None)
        self.assertTrue(result['moi_compatible'])
        self.assertEqual(result['score'], 1.0)

    def test_legacy_assessment_never_filters(self):
        result = self.recessive.assess(self.het, 'TEST:A')
        self.assertNotIn('moi_compatible', result)
        self.assertEqual(result['score'], 0.2)


class RankingV13Tests(unittest.TestCase):
    def setUp(self):
        self.kb = recessive_knowledge()
        self.config = RankingConfig.load(ROOT / 'config/ranking.json')
        self.assertEqual(self.config.molecular_mode, 'context_v3')
        # The fixture genes have no curated molecular rules; legacy scoring keeps it small.
        self.config = replace(self.config, molecular_mode='legacy_v1')

    def run_case(self, variants, annotations, mode='context_v3'):
        config = replace(self.config, molecular_mode=mode)
        case = Case('fixture', variants, {'HP:0000002': 'feature A'}, {})
        candidates = rank(case, annotations, self.kb, config)
        return candidates, config, case

    def test_v13_mode_uses_filter_and_abstains_on_lone_recessive_heterozygote(self):
        lone = Variant('chr1', 100, 'A', 'T', '0/1')
        annotations = {lone.id: annotation(lone, 'GENEA', 'stop_gained')}
        from unittest.mock import patch
        with patch('genomic_detective.experts.MolecularExpert.assess_facts',
                   lambda self, facts, disease_id=None: facts['baseline']):
            candidates, config, case = self.run_case([lone], annotations)
            check_ranking(case, candidates, annotations, config)
        self.assertFalse(candidates[0]['experts']['inheritance']['moi_compatible'])
        self.assertEqual(decide(candidates, 1.0, True, config),
                         ('no_convincing_candidate', 'no_causal_variant_in_vcf'))

    def test_compatible_hypothesis_is_ranked_before_filtered_one(self):
        lone = Variant('chr1', 100, 'A', 'T', '0/1')
        hom = Variant('chr1', 300, 'C', 'G', '1/1')
        annotations = {lone.id: annotation(lone, 'GENEA', 'stop_gained'),
                       hom.id: annotation(hom, 'GENEB', 'missense_variant')}
        from unittest.mock import patch
        with patch('genomic_detective.experts.MolecularExpert.assess_facts',
                   lambda self, facts, disease_id=None: facts['baseline']):
            candidates, config, case = self.run_case([lone, hom], annotations)
            check_ranking(case, candidates, annotations, config)
        self.assertEqual(candidates[0]['variant_id'], hom.id)
        self.assertEqual(candidates[1]['variant_id'], lone.id)


class MolecularV13Tests(unittest.TestCase):
    def setUp(self):
        self.rules = json.loads((ROOT / 'config/molecular_knowledge.json').read_text(encoding='utf-8'))
        self.variant = Variant('1', 156138613, 'C', 'T', '0/1')
        self.annotation = annotation(self.variant, 'LMNA', 'missense_variant', af=None)

    def store(self, assertions=()):
        snapshot = {'schema_version': 1, 'assembly': 'GRCh38', 'fetched_at': 'fixture', 'scope': 'test', 'records': {}}
        if assertions:
            snapshot['records'][self.variant.id] = [{'variant_id': self.variant.id, 'assembly': 'GRCh38',
                'genes': ['LMNA'], 'accession': 'VCV000000001.1', 'assertions': list(assertions)}]
        return MolecularKnowledge(self.rules, snapshot)

    def assess(self, mode, store=None, disease='OMIM:176670'):
        return MolecularExpert(mode, store or self.store()).assess(
            self.variant, self.annotation, 'LMNA', disease, literature=None)

    def test_absence_from_gnomad_is_weak_rarity_only_in_v13(self):
        self.assertEqual(self.assess('context_v2')['components']['frequency'], 0)
        result = self.assess('context_v3')
        self.assertEqual(result['components']['frequency'], 0.5)
        self.assertTrue(any(e['criterion'] == 'population_frequency' for e in result['evidence']))

    def test_failed_annotation_never_becomes_rarity(self):
        self.annotation['status'] = 'unavailable'
        self.assertIsNone(self.assess('context_v3')['score'])

    def test_reviewed_unspecified_condition_gives_gene_level_support(self):
        generic = assertion(['MEDGEN:C3661900'], review='criteria provided, single submitter')
        generic['conditions'] = ['not provided']
        store = self.store([generic])
        self.assertEqual(self.assess('context_v2', store)['components']['clinical'], 0)
        result = self.assess('context_v3', store)
        self.assertEqual(result['components']['clinical'], 0.5)
        self.assertEqual(result['clinical_context_status'], 'gene_level')

    def test_unreviewed_or_benign_unspecified_records_give_no_support(self):
        unreviewed = assertion(['MEDGEN:C1'], review='no assertion criteria provided')
        unreviewed['conditions'] = ['not provided']
        self.assertEqual(self.assess('context_v3', self.store([unreviewed]))['components']['clinical'], 0)
        benign = assertion(['MEDGEN:C1'], classification='Benign', review='criteria provided, single submitter')
        benign['conditions'] = ['not specified']
        self.assertEqual(self.assess('context_v3', self.store([benign]))['components']['clinical'], 0)

    def test_other_named_disease_does_not_count_as_unspecified(self):
        other = assertion(['OMIM:999999'], review='criteria provided, single submitter')
        other['conditions'] = ['Some other syndrome']
        self.assertEqual(self.assess('context_v3', self.store([other]))['components']['clinical'], 0)


if __name__ == '__main__':
    unittest.main()
