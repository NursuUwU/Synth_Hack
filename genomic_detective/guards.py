"""Mechanical integrity guards; not an independent validation of medical truth."""
from __future__ import annotations

from .models import Case, finite_score
from .ranking import RankingConfig, aggregate


class GuardError(ValueError):
    pass


def check_ranking(case: Case, candidates: list[dict], annotations: dict, config: RankingConfig):
    expected = {variant.id: variant for variant in case.variants}
    if len(candidates) != len(expected) or {c['variant_id'] for c in candidates} != set(expected):
        raise GuardError('Candidate set does not match input variants')
    previous = float('inf')
    seen_evidence = {}
    for position, candidate in enumerate(candidates, 1):
        variant = expected[candidate['variant_id']]
        if candidate['assembly'] != variant.assembly or candidate['genotype'] != variant.genotype:
            raise GuardError('Candidate assembly/genotype changed')
        score = candidate['score']
        if not finite_score(score) or not 0 <= score <= 1 or score > previous + 1e-12:
            raise GuardError('Invalid or unsorted ranking score')
        if candidate['rank'] != position:
            raise GuardError('Incorrect rank number')
        previous = score
        if set(candidate['experts']) != set(config.weights):
            raise GuardError('Missing expert result')
        for assessment in candidate['experts'].values():
            value = assessment['score']
            if value is not None and (not finite_score(value) or not 0 <= value <= 1):
                raise GuardError('Invalid expert score')
        if abs(aggregate(candidate['experts'], config) - score) > 1e-12:
            raise GuardError('Aggregate does not match expert scores')
        annotation = annotations[candidate['variant_id']]
        if annotation.get('gnomad_max_af') is None:
            if candidate['experts']['molecular'].get('components', {}).get('frequency', 0) != 0:
                raise GuardError('Missing AF was converted into rarity evidence')
        gene = candidate['gene']
        if gene and gene not in {t.get('gene_symbol') for t in annotation.get('transcripts', [])}:
            raise GuardError('Gene does not come from variant annotation')
        allowed_sources = {f'VEP:{variant.cache_key}'}
        if candidate['disease_id']:
            allowed_sources.add(f'HPO:{candidate["disease_id"]}')
        expected_evidence = [ev for expert in candidate['experts'].values() for ev in expert['evidence']]
        if candidate['evidence'] != expected_evidence:
            raise GuardError('Evidence differs from expert output')
        for item in candidate['evidence']:
            if item['source'] not in allowed_sources:
                raise GuardError('Evidence refers to an unknown source')
            if item['id'] in seen_evidence and seen_evidence[item['id']] != item:
                raise GuardError('Evidence identifier collision')
            seen_evidence[item['id']] = item
    return {'status': 'passed', 'checks': [
        'input_variant_identity', 'assembly_and_genotype', 'finite_scores', 'score_recalculation',
        'rank_order', 'missing_frequency', 'gene_provenance', 'evidence_references'],
        'scope': 'Mechanical integrity; source correctness and clinical validity require expert review'}


def check_presentation(selection: dict, report: dict, count: int = 3):
    expected = report['candidates'][:count]
    if not isinstance(selection, dict):
        raise GuardError('Presentation must be a JSON object')
    if selection.get('case_id') != report['case_id'] or selection.get('status') != report['status']:
        raise GuardError('Presentation changed case or status')
    items = selection.get('candidates', [])
    if (not isinstance(items, list) or not all(isinstance(item, dict) for item in items)
            or [item.get('variant_id') for item in items] != [c['variant_id'] for c in expected]):
        raise GuardError('Presentation changed candidate set or rank order')
    for item, candidate in zip(items, expected):
        selected = item.get('evidence_ids')
        allowed = {ev['id'] for ev in candidate['evidence']}
        if (not isinstance(selected, list) or not all(isinstance(eid, str) for eid in selected)
                or len(selected) != len(set(selected)) or not set(selected).issubset(allowed)):
            raise GuardError('Presentation contains unknown or repeated evidence')
        mandatory = {ev['id'] for ev in candidate['evidence']
                     if ev['criterion'] in {'phenotype_conflict', 'clinical_annotation', 'inheritance',
                                             'consequence', 'population_frequency'}}
        if not mandatory.issubset(selected):
            raise GuardError('Presentation omitted conflicting or qualified evidence')
    return {'status': 'passed', 'checks': ['case_status', 'candidate_order', 'evidence_ids', 'required_caveats'],
            'scope': 'Facts rendered verbatim from evidence; LLM cannot add free-form claims'}
