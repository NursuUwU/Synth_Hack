from __future__ import annotations

import copy
import unittest
from dataclasses import replace
from pathlib import Path

from genomic_detective.experts import MolecularExpert
from genomic_detective.guards import GuardError, check_presentation, check_ranking
from genomic_detective.knowledge import Knowledge
from genomic_detective.literature import LiteratureError, LiteratureRegistry
from genomic_detective.models import Case, Variant
from genomic_detective.ranking import RankingConfig, rank
from genomic_detective.reporting import template_selection
from test_workflow import annotation


ROOT = Path(__file__).resolve().parents[1]
LMNA = 'chr1:156138613:C:T'


def entry(**changes):
    record = {
        'id': 'lmna_g608g_splice', 'assembly': 'GRCh38', 'variant': LMNA, 'gene': 'LMNA',
        'transcript': 'NM_170707.4', 'hgvsc': 'c.1824C>T', 'hgvsp': 'p.Gly608=',
        'disease_ids': ['OMIM:176670', 'ORPHA:740'],
        'crosswalk_source': {'url': 'https://api.orphadata.com/rd-cross-referencing/orphacodes/740',
                             'checked_on': '2026-11-12'},
        'evidence_scope': 'exact_variant', 'evidence_type': 'functional_experiment',
        'support': 0.85,
        'fact': 'Для G608G описаны аномальная обработка РНК и укороченный белковый продукт.',
        'source': {'pmid': '12714972', 'material': 'Figure 2, Table 1'},
        'not_established': ['RNA-эксперимент для входного образца не проводился'],
        'approved_by': 'participant_4',
    }
    record.update(changes)
    return record


def registry(*entries):
    return LiteratureRegistry({'schema_version': 1, 'version': '1.2',
                               'status': 'approved_for_research_ranking',
                               'reviewed_on': '2026-11-12',
                               'entries': list(entries or (entry(),))})


def lmna_knowledge():
    kb = Knowledge()
    kb.parents = {'HP:0000001': set(), 'HP:0000002': {'HP:0000001'}}
    kb.labels = {'HP:0000001': 'root', 'HP:0000002': 'Progeroid feature'}
    disease = 'ORPHA:740'
    kb.diseases[disease] = {'id': disease, 'name': 'Hutchinson-Gilford progeria',
        'positive': {'HP:0000002'}, 'negative': set(), 'inheritance': {'AD'},
        'genes': {'LMNA'}, 'references': set()}
    kb.gene_diseases['LMNA'].add(disease)
    kb.build_index()
    return kb


def lmna_annotation(variant, consequence='synonymous_variant', af=0.00001):
    result = annotation(variant, 'LMNA', consequence, af)
    result['transcripts'][0].update(transcript_id='NM_170707.4', mane_select='NM_170707.4',
                                    hgvsc='NM_170707.4:c.1824C>T')
    return result


class LiteratureTests(unittest.TestCase):
    def setUp(self):
        self.variant = Variant('1', 156138613, 'C', 'T', '0/1')
        self.annotation = lmna_annotation(self.variant)
        self.kb = lmna_knowledge()
        self.case = Case('literature_fixture', [self.variant],
                         {'HP:0000002': 'Progeroid feature'}, {})
        self.config = replace(RankingConfig.load(ROOT / 'config/ranking.json'),
                              molecular_mode='legacy_v1')
        self.registry = registry()

    def assess(self, variant=None, ann=None, disease='ORPHA:740', literature=None):
        return MolecularExpert('legacy_v1').assess(
            variant or self.variant, ann or self.annotation, 'LMNA', disease,
            literature or self.registry)

    def baseline_score(self):
        return 0.55 * 0.1 + 0.15 * 1.0

    # ТЗ-1.13
    def test_other_synonymous_position_gets_no_support(self):
        variant = Variant('1', 156138614, 'C', 'T', '0/1')
        ann = lmna_annotation(variant)
        result = self.assess(variant, ann)
        self.assertEqual(result['components']['consequence'], 0.1)
        self.assertNotIn('literature_applied', result)
        self.assertFalse(any(ev['criterion'] == 'literature_functional_effect' for ev in result['evidence']))

    # ТЗ-1.14
    def test_same_site_other_alt_gets_no_support(self):
        variant = Variant('1', 156138613, 'C', 'A', '0/1')
        ann = lmna_annotation(variant)
        result = self.assess(variant, ann)
        self.assertEqual(result['components']['consequence'], 0.1)
        self.assertNotIn('literature_applied', result)

    # ТЗ-1.15
    def test_same_variant_other_disease_gets_no_support(self):
        result = self.assess(disease='OMIM:999999')
        self.assertEqual(result['components']['consequence'], 0.1)
        self.assertNotIn('literature_applied', result)

    # ТЗ-1.16
    def test_entry_without_source_or_approval_is_not_loaded(self):
        payload = registry(entry(id='no_source', source={}),
                           entry(id='no_approval', approved_by='')).entries
        self.assertEqual(payload, [])

    # ТЗ-1.17
    def test_same_site_other_allele_is_context_only(self):
        context = registry(entry(id='lmna_other_alt', variant='chr1:156138613:C:A',
                                 evidence_scope='same_site_other_allele', support=0.85))
        result = self.assess(literature=context)
        self.assertEqual(result['components']['consequence'], 0.1)
        self.assertNotIn('literature_applied', result)
        self.assertTrue(any(ev['criterion'] == 'literature_context' for ev in result['evidence']))

    # ТЗ-1.18
    def test_exact_support_uses_max_not_sum(self):
        result = self.assess()
        self.assertEqual(result['components']['consequence'], 0.85)
        self.assertAlmostEqual(result['score'], 0.55 * 0.85 + 0.15 * 1.0)
        self.assertNotAlmostEqual(result['score'], 0.55 * 0.85 + 0.85 + 0.15 * 1.0)
        self.assertNotAlmostEqual(result['score'], self.baseline_score() + 0.85)

    def ranked_candidate(self):
        annotations = {self.variant.id: self.annotation}
        candidates = rank(self.case, annotations, self.kb, self.config, self.registry)
        candidate = next(item for item in candidates if item['gene'] == 'LMNA')
        return annotations, candidates, candidate

    # ТЗ-1.19
    def test_evidence_from_outside_registry_is_rejected(self):
        annotations, candidates, candidate = self.ranked_candidate()
        foreign = {'id': 'ev_foreign', 'criterion': 'literature_context',
                   'text': 'Invented foreign source', 'source': 'LIT:not_in_registry',
                   'basis': 'literature'}
        candidate['evidence'].append(copy.deepcopy(foreign))
        candidate['experts']['phenotype']['evidence'].append(copy.deepcopy(foreign))
        with self.assertRaises(GuardError):
            check_ranking(self.case, candidates, annotations, self.config, self.registry)

    # ТЗ-1.20
    def test_presentation_must_keep_literature_evidence(self):
        _, _, candidate = self.ranked_candidate()
        report = {'case_id': self.case.id, 'status': 'candidate', 'candidates': [candidate]}
        selection = template_selection(report)
        literature_ids = [ev['id'] for ev in candidate['evidence']
                          if ev['criterion'] == 'literature_functional_effect']
        self.assertTrue(literature_ids)
        selection['candidates'][0]['evidence_ids'] = [
            eid for eid in selection['candidates'][0]['evidence_ids'] if eid not in literature_ids]
        with self.assertRaises(GuardError):
            check_presentation(selection, report)

    def test_crosswalk_allows_orphanet_hypothesis(self):
        result = self.assess(disease='ORPHA:740')
        self.assertEqual(result['literature_applied'], ['lmna_g608g_splice'])

    def test_invalid_schema_rejected(self):
        with self.assertRaises(LiteratureError):
            LiteratureRegistry({'schema_version': 99, 'entries': []})


if __name__ == '__main__':
    unittest.main()
