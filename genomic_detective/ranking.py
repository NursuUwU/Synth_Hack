from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from .experts import InheritanceExpert, MolecularExpert, PhenotypeExpert
from .knowledge import Knowledge
from .models import CONTEXT_MODES, Case, finite_score


# Research heuristic shared with the UI; not a calibrated confidence interval.
CANDIDATE_GAP = 0.025


@dataclass
class RankingConfig:
    version: str
    weights: dict
    candidate_min_score: float
    candidate_min_phenotype: float
    candidate_min_molecular: float
    min_annotation_coverage: float
    note: str = ''
    molecular_mode: str = 'legacy_v1'

    @classmethod
    def load(cls, path: Path):
        result = cls(**json.loads(path.read_text(encoding='utf-8')))
        result.validate()
        return result

    def validate(self):
        if self.molecular_mode not in {'legacy_v1'} | CONTEXT_MODES:
            raise ValueError('Unknown molecular mode')
        if set(self.weights) != {'phenotype', 'molecular', 'inheritance'}:
            raise ValueError('Exactly three expert weights are required')
        numbers = list(self.weights.values()) + [self.candidate_min_score, self.candidate_min_phenotype,
                  self.candidate_min_molecular, self.min_annotation_coverage]
        if not all(finite_score(x) and 0 <= x <= 1 for x in numbers):
            raise ValueError('Weights and thresholds must be finite values in [0,1]')
        if abs(sum(self.weights.values()) - 1.0) > 1e-9:
            raise ValueError('Expert weights must sum to 1')

    def dictionary(self):
        return vars(self)

    @property
    def digest(self):
        return hashlib.sha256(json.dumps(self.dictionary(), sort_keys=True).encode()).hexdigest()


def moi_compatible(candidate: dict) -> bool:
    # Only v1.3 inheritance assessments carry the flag; earlier versions never filter.
    return candidate['experts']['inheritance'].get('moi_compatible', True)


def order_key(candidate: dict):
    # v1.3 places MOI-incompatible hypotheses after compatible ones (Exomiser-style filter).
    return (not moi_compatible(candidate), -candidate['score'])


def rare_alleles_by_gene(case: Case, annotations: dict) -> dict:
    from .experts import MOI_MAX_AF
    genes = {}
    for variant in case.variants:
        annotation = annotations[variant.id]
        af = annotation.get('gnomad_max_af')
        if annotation.get('status') != 'annotated' or af is not None and af > MOI_MAX_AF['AR_comphet']:
            continue
        for entry in annotation.get('transcripts', []):
            if entry.get('gene_symbol'):
                genes.setdefault(entry['gene_symbol'], set()).add(variant.id)
    return genes


def clinvar_whitelisted(molecular: dict) -> bool:
    # Disease-matched, reviewed (>=2 stars) P/LP ClinVar assertion; never gene-level records.
    return molecular.get('clinical_context_status') == 'matched' and any(
        a['review_stars'] >= 2 and 'pathogenic' in a['classification'].lower()
        and 'conflict' not in a['classification'].lower() and 'benign' not in a['classification'].lower()
        for a in molecular.get('clinical_assertions', []))


def aggregate(experts: dict, config: RankingConfig) -> float:
    return sum(config.weights[name] * (assessment['score'] if assessment['score'] is not None else 0.0)
               for name, assessment in experts.items())


def rank(case: Case, annotations: dict, knowledge: Knowledge, config: RankingConfig,
         literature=None) -> list[dict]:
    if literature is None and config.molecular_mode in CONTEXT_MODES:
        from .literature import default_registry
        literature = default_registry()
    phenotype = PhenotypeExpert(knowledge)
    molecular = MolecularExpert(config.molecular_mode)
    inheritance = InheritanceExpert(knowledge)
    moi_filter = config.molecular_mode == 'context_v3'
    rare_alleles = rare_alleles_by_gene(case, annotations) if moi_filter else {}
    candidates = []
    for variant in case.variants:
        annotation = annotations[variant.id]
        genes = sorted({entry['gene_symbol'] for entry in annotation.get('transcripts', [])
                        if entry.get('gene_symbol')})
        alternatives = []
        for gene in genes or [None]:
            hypotheses = phenotype.hypotheses(case, gene) if gene else [{
                'score': None, 'disease_id': None, 'disease_name': None,
                'evidence': [], 'missing': ['No annotated gene'], 'matches': [], 'conflicts': []}]
            for pheno in hypotheses:
                mol = molecular.assess(variant, annotation, gene, pheno['disease_id'], literature)
                experts = {'phenotype': pheno, 'molecular': mol,
                           'inheritance': inheritance.assess(
                               variant, pheno['disease_id'], moi_filter, annotation.get('gnomad_max_af'),
                               bool(rare_alleles.get(gene, set()) - {variant.id}), clinvar_whitelisted(mol))}
                alternatives.append({'variant_id': variant.id, 'assembly': variant.assembly,
                    'genotype': variant.genotype, 'gene': gene, 'disease_id': pheno['disease_id'],
                    'disease_name': pheno['disease_name'], 'score': aggregate(experts, config),
                    'experts': experts,
                    'evidence': [item for assessment in experts.values() for item in assessment['evidence']],
                    'missing': sorted({item for assessment in experts.values() for item in assessment['missing']}),
                    'annotation_status': annotation['status']})
        ordered = sorted(alternatives, key=lambda item: (*order_key(item), item['gene'] or '', item['disease_id'] or ''))
        best = ordered[0]
        gene_alternatives = {}
        for item in ordered:
            if item['gene'] != best['gene']:
                gene_alternatives.setdefault(item['gene'], item)
        best['other_gene_hypotheses'] = [
            {'gene': item['gene'], 'disease_id': item['disease_id'], 'score': item['score']}
            for item in gene_alternatives.values()]
        best['other_disease_hypotheses'] = [
            {'disease_id': item['disease_id'], 'score': item['score']}
            for item in ordered if item is not best and item['gene'] == best['gene']][:5]
        candidates.append(best)
    candidates.sort(key=lambda item: (*order_key(item), item['variant_id']))
    for position, candidate in enumerate(candidates, 1):
        candidate['rank'] = position
    return candidates


def uncertainty_reasons(candidates, config):
    if not candidates:
        return []
    reasons = []
    rivals = [c for c in candidates[1:] if moi_compatible(c)]
    if rivals and candidates[0]['score'] - rivals[0]['score'] < CANDIDATE_GAP:
        reasons.append('close_scores')
    if config.molecular_mode in CONTEXT_MODES:
        assessment = candidates[0]['experts']['molecular']
        if assessment.get('clinical_conflict'):
            reasons.append('clinical_conflict')
        if assessment.get('functional_clinical_conflict'):
            reasons.append('functional_clinical_conflict')
        components = assessment.get('components', {})
        if (not components.get('clinical', 0) and not components.get('functional', 0)
                and not assessment.get('literature_applied')):
            reasons.append('limited_condition_evidence')
        elif assessment.get('clinical_context_status') == 'gene_level' and not components.get('functional', 0):
            reasons.append('gene_level_clinical_only')
    return reasons


def decide(candidates: list[dict], coverage: float, knowledge_ready: bool, config: RankingConfig):
    if coverage < config.min_annotation_coverage or not knowledge_ready:
        return 'incomplete', None
    if candidates:
        best = candidates[0]
        if not moi_compatible(best):
            return 'no_convincing_candidate', 'no_causal_variant_in_vcf'
        phenotype = best['experts']['phenotype']['score']
        molecular = best['experts']['molecular']['score']
        if (best['score'] >= config.candidate_min_score
                and phenotype is not None and phenotype >= config.candidate_min_phenotype
                and molecular is not None and molecular >= config.candidate_min_molecular):
            return ('uncertain' if uncertainty_reasons(candidates, config) else 'candidate'), best['variant_id']
    # This is the benchmark's abstention token, not proof that no biological cause exists.
    return 'no_convincing_candidate', 'no_causal_variant_in_vcf'
