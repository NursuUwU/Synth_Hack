import copy
import json
import unittest
from pathlib import Path
from unittest.mock import patch

from genomic_detective.clinvar import parse_vcv, review_stars
from genomic_detective.experts import MolecularExpert
from genomic_detective.guards import GuardError, check_presentation, check_ranking
from genomic_detective.models import Case, Variant
from genomic_detective.molecular import MolecularKnowledge
from genomic_detective.pipeline import Pipeline
from genomic_detective.ranking import RankingConfig, aggregate, decide, rank, uncertainty_reasons
from genomic_detective.reporting import template_selection
from test_workflow import FakeProvider, annotation, knowledge

ROOT = Path(__file__).resolve().parents[1]


def assertion(ids, classification='Pathogenic', review='criteria provided, multiple submitters, no conflicts', compound=False, number=1):
    accession = f'RCV{number:09}.1'
    return {'accession': accession, 'classification': classification, 'review_status': review,
            'review_stars': review_stars(review), 'last_evaluated': '2026-01-01',
            'conditions': ['fixture condition'], 'condition_ids': ids,
            'compound_condition': compound, 'source': 'ClinVar:' + accession,
            'url': 'https://www.ncbi.nlm.nih.gov/clinvar/' + accession + '/'}


class ContextTests(unittest.TestCase):
    def setUp(self):
        self.rules = json.loads((ROOT / 'config/molecular_knowledge.json').read_text(encoding='utf-8'))
        self.variant = Variant('1', 156138613, 'C', 'T', '0/1')
        self.annotation = annotation(self.variant, 'LMNA', 'synonymous_variant', af=None)
        self.annotation['transcripts'][0].update(mane_select='NM_170707.4', hgvsc='TEST:c.1824C>T')
        self.config = RankingConfig.load(ROOT / 'config/ranking.json')

    def store(self, assertions=()):
        snapshot = {'schema_version': 1, 'assembly': 'GRCh38', 'fetched_at': 'fixture', 'scope': 'test', 'records': {}}
        if assertions:
            snapshot['records'][self.variant.id] = [{'variant_id': self.variant.id, 'assembly': 'GRCh38',
                'genes': ['LMNA'], 'accession': 'VCV000000001.1', 'assertions': list(assertions)}]
        return MolecularKnowledge(self.rules, snapshot)

    def assess(self, store=None, disease='OMIM:176670', variant=None, ann=None):
        return MolecularExpert('context_v2', store or self.store()).assess(
            variant or self.variant, ann or self.annotation, 'LMNA', disease)

    def test_synonymous_preserved_and_exact_functional_effect_added(self):
        result = self.assess()
        self.assertEqual(result['consequences'], ['synonymous_variant'])
        self.assertEqual(result['components']['functional'], 0.85)
        self.assertTrue(any(e['criterion'] == 'functional_effect' for e in result['evidence']))
        self.assertEqual(self.annotation['transcripts'][0]['consequence_terms'], ['synonymous_variant'])

    def test_other_alt_or_position_receives_no_known_function(self):
        for variant in [Variant('1', self.variant.pos, 'C', 'A', '0/1'),
                        Variant('1', self.variant.pos + 1, 'C', 'T', '0/1')]:
            ann = copy.deepcopy(self.annotation)
            ann['variant_id'] = variant.id
            self.assertEqual(self.assess(variant=variant, ann=ann)['components']['functional'], 0)

    def test_wrong_transcript_or_hgvs_receives_no_known_function(self):
        for changes in [{'mane_select': 'NM_OTHER.1'}, {'hgvsc': 'TEST:c.1825C>T'}]:
            ann = copy.deepcopy(self.annotation)
            ann['transcripts'][0].update(changes)
            self.assertEqual(self.assess(ann=ann)['components']['functional'], 0)

    def test_other_disease_does_not_inherit_progeria_evidence(self):
        first = self.assess()
        other = self.assess(disease='OMIM:115200')
        self.assertGreater(first['score'], other['score'])
        self.assertEqual(other['components']['functional'], 0)

    def test_single_condition_xrefs_match_but_compound_does_not(self):
        ids = ['OMIM:176670', 'ORPHA:740']
        matched = self.assess(store=self.store([assertion(ids)]), disease='ORPHA:740')
        self.assertEqual(matched['components']['functional'], 0.85)
        compound = self.assess(store=self.store([assertion(ids, compound=True)]), disease='ORPHA:740')
        self.assertEqual(compound['components']['functional'], 0)
        self.assertEqual(compound['components']['clinical'], 0)

    def test_dependent_evidence_uses_max_not_sum(self):
        result = self.assess(store=self.store([assertion(['OMIM:176670'])]))
        self.assertEqual(result['components']['clinical'], 0.8)
        self.assertEqual(result['score'], 0.85)
        self.assertEqual(result['components']['frequency'], 0)

    def test_unscoped_vep_pathogenic_does_not_supply_condition_points(self):
        ann = copy.deepcopy(self.annotation)
        ann['clinical'] = [{'significance': ['pathogenic']}]
        result = self.assess(ann=ann, disease='TEST:OTHER')
        self.assertEqual(result['components']['clinical'], 0)
        self.assertEqual(result['clinical_context_status'], 'unresolved')

    def test_other_condition_and_no_criteria_classification_are_not_positive(self):
        for record in [assertion(['TEST:OTHER']),
                       assertion(['OMIM:176670'], review='no assertion criteria provided')]:
            result = self.assess(store=self.store([record]))
            self.assertEqual(result['components']['clinical'], 0)
        result = self.assess(store=self.store([assertion(['OMIM:176670'])]))
        self.assertEqual(result['components']['clinical'], 0.8)

    def test_conflict_keeps_both_facts_and_forces_uncertainty(self):
        store = self.store([assertion(['OMIM:176670'], 'Benign')])
        result = self.assess(store=store)
        candidate = {'variant_id': self.variant.id, 'score': 0.8, 'experts': {
            'molecular': result, 'phenotype': {'score': 0.9}, 'inheritance': {'score': 1}}}
        self.assertTrue(result['functional_clinical_conflict'])
        self.assertTrue(any(e['criterion'] == 'evidence_conflict' for e in result['evidence']))
        self.assertEqual(decide([candidate], 1, True, self.config)[0], 'uncertain')

    def test_pathogenic_and_vus_across_matching_records_are_not_silently_merged(self):
        store = self.store([assertion(['OMIM:176670']),
                            assertion(['OMIM:176670'], 'Uncertain significance', number=2)])
        result = self.assess(store=store)
        self.assertEqual(result['clinical_status'], 'conflicting')
        self.assertEqual(result['components']['clinical'], 0)
        self.assertTrue(result['clinical_conflict'])

    def test_limited_evidence_is_not_confused_with_small_score_gap(self):
        result = self.assess(disease='TEST:OTHER')
        # Keep the example above the molecular decision threshold without creating
        # a reviewed classification or functional observation.
        result['score'] = 0.4
        candidate = {'variant_id': 'fixture', 'score': 0.8, 'experts': {
            'molecular': result, 'phenotype': {'score': 0.9}, 'inheritance': {'score': 1}}}
        self.assertEqual(uncertainty_reasons([candidate], self.config), ['limited_condition_evidence'])
        self.assertEqual(decide([candidate], 1, True, self.config)[0], 'uncertain')

    def test_null_prediction_requires_nmd_and_disease_context(self):
        variant = Variant('1', 100, 'A', 'T', '0/1')
        expert = MolecularExpert('context_v2', self.store())
        ann = annotation(variant, 'PTPN11', 'stop_gained', af=None)
        result = expert.assess(variant, ann, 'PTPN11', 'OMIM:163950')
        self.assertEqual(result['components']['consequence'], 0)
        ann = annotation(variant, 'KMT2D', 'stop_gained', af=None)
        result = expert.assess(variant, ann, 'KMT2D', 'OMIM:147920')
        self.assertLessEqual(result['components']['consequence'], 0.55)
        self.assertTrue(any('NMD' in s for s in result['missing']))

    def test_snapshot_rejects_wrong_allele_and_forged_review_level(self):
        store = self.store([assertion(['OMIM:176670'])])
        for field, value in [('variant_id', 'chr1:1:A:T'), ('assembly', 'GRCh37')]:
            snapshot = copy.deepcopy(store.snapshot)
            snapshot['records'][self.variant.id][0][field] = value
            with self.assertRaises(ValueError):
                MolecularKnowledge(self.rules, snapshot)
        snapshot = copy.deepcopy(store.snapshot)
        snapshot['records'][self.variant.id][0]['assertions'][0]['review_stars'] = 4
        with self.assertRaises(ValueError):
            MolecularKnowledge(self.rules, snapshot)

    def test_wrong_annotation_assembly_is_rejected(self):
        ann = copy.deepcopy(self.annotation)
        ann['assembly'] = 'GRCh37'
        with self.assertRaises(ValueError):
            self.assess(ann=ann)

    def test_facts_extracted_once_and_all_disease_hypotheses_evaluated(self):
        case = Case('fixture', [self.variant], {'HP:0000002': 'feature A'}, {})
        kb = knowledge()
        kb.gene_diseases['LMNA'] = set(kb.diseases)
        original = MolecularExpert.extract_facts
        with patch.object(MolecularExpert, 'extract_facts', autospec=True, side_effect=original) as extract:
            candidates = rank(case, {self.variant.id: self.annotation}, kb, self.config)
        self.assertEqual(extract.call_count, 1)
        self.assertEqual(len(candidates[0]['other_disease_hypotheses']), 1)

    def test_case_identifier_cannot_change_ranking(self):
        kb = knowledge()
        kb.gene_diseases['LMNA'] = set(kb.diseases)
        first = Case('case_05', [self.variant], {'HP:0000002': 'feature A'}, {})
        renamed = copy.deepcopy(first)
        renamed.id = 'unseen_identifier'
        anns = {self.variant.id: self.annotation}
        self.assertEqual(rank(first, anns, kb, self.config), rank(renamed, anns, kb, self.config))

    def test_guard_replays_molecular_context_not_only_aggregate(self):
        variant = Variant('1', 100, 'A', 'T', '0/1')
        case = Case('fixture', [variant], {'HP:0000002': 'feature A'}, {})
        anns = {variant.id: annotation(variant, 'GENEA')}
        candidates = rank(case, anns, knowledge(), self.config)
        candidates[0]['experts']['molecular']['components']['functional'] = 0.99
        candidates[0]['experts']['molecular']['score'] = 0.99
        candidates[0]['score'] = aggregate(candidates[0]['experts'], self.config)
        with self.assertRaises(GuardError):
            check_ranking(case, candidates, anns, self.config)

    def test_presentation_cannot_drop_functional_caveats(self):
        store = self.store()
        case = Case('fixture', [self.variant], {'HP:0000002': 'feature A'}, {})
        kb = knowledge()
        kb.diseases['OMIM:176670'] = copy.deepcopy(kb.diseases['TEST:A'])
        kb.diseases['OMIM:176670']['id'] = 'OMIM:176670'
        kb.gene_diseases['LMNA'] = {'OMIM:176670'}
        with patch('genomic_detective.molecular.default_knowledge', return_value=store):
            report = Pipeline(FakeProvider({self.variant.id: self.annotation}), kb, self.config).run([case])[0]
        selection = template_selection(report)
        functional = {e['id'] for e in report['candidates'][0]['evidence'] if e['criterion'] == 'functional_effect'}
        selection['candidates'][0]['evidence_ids'] = [e for e in selection['candidates'][0]['evidence_ids'] if e not in functional]
        with self.assertRaises(GuardError):
            check_presentation(selection, report)


class ClinVarParserTests(unittest.TestCase):
    def test_only_variant_location_and_germline_single_condition_xrefs(self):
        payload = b'''<ClinVarResult-Set><VariationArchive Accession="VCV000000001" Version="2" RecordType="classified">
          <ClassifiedRecord><SimpleAllele><GeneList><Gene Symbol="TEST"><Location><SequenceLocation Assembly="GRCh38" Chr="1" positionVCF="999" referenceAlleleVCF="A" alternateAlleleVCF="T"/></Location></Gene></GeneList>
          <Location><SequenceLocation Assembly="GRCh38" Accession="NC_000001.11" Chr="1" positionVCF="100" referenceAlleleVCF="A" alternateAlleleVCF="T"/>
          <SequenceLocation Assembly="GRCh38" Accession="NW_12345.1" Chr="1" positionVCF="200" referenceAlleleVCF="A" alternateAlleleVCF="T"/></Location></SimpleAllele>
          <RCVList><RCVAccession Accession="RCV000000001" Version="3"><ClassifiedConditionList><ClassifiedCondition DB="MedGen" ID="C1">condition A</ClassifiedCondition></ClassifiedConditionList>
          <RCVClassifications><GermlineClassification><ReviewStatus>reviewed by expert panel</ReviewStatus><Description>Pathogenic</Description></GermlineClassification><OncogenicityClassification><Description>Oncogenic</Description></OncogenicityClassification></RCVClassifications></RCVAccession></RCVList>
          <ClinicalAssertionList><TraitSet><Trait><XRef DB="MedGen" ID="C1"/><XRef DB="OMIM" ID="123456"/></Trait><Trait><XRef DB="MedGen" ID="C2"/><XRef DB="OMIM" ID="654321"/></Trait></TraitSet></ClinicalAssertionList>
          </ClassifiedRecord></VariationArchive></ClinVarResult-Set>'''
        records = parse_vcv(payload)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record['variant_id'], 'chr1:100:A:T')
        self.assertEqual(record['assertions'][0]['condition_ids'], ['MEDGEN:C1', 'OMIM:123456'])
        self.assertEqual(record['assertions'][0]['classification'], 'Pathogenic')
        self.assertEqual(record['assertions'][0]['review_stars'], 3)

    def test_haplotype_classification_is_not_inherited_by_component(self):
        payload = b'<ClinVarResult-Set><VariationArchive RecordType="classified"><ClassifiedRecord><Haplotype><SimpleAllele/></Haplotype></ClassifiedRecord></VariationArchive></ClinVarResult-Set>'
        self.assertEqual(parse_vcv(payload), [])


if __name__ == '__main__':
    unittest.main()
