"""Mechanical integrity guards; not an independent validation of medical truth."""
from __future__ import annotations

from .models import CONTEXT_MODES, Case, finite_score
from .ranking import RankingConfig, aggregate, order_key


REQUIRED_CRITERIA = frozenset({'phenotype_conflict', 'clinical_annotation', 'inheritance',
    'consequence', 'population_frequency', 'clinical_assertion', 'clinical_context',
    'mechanism_context', 'functional_effect', 'literature_context', 'evidence_conflict',
    'literature_functional_effect'})


class GuardError(ValueError):
    pass


def check_ranking(case: Case, candidates: list[dict], annotations: dict, config: RankingConfig,
                  literature=None):
    from .experts import MolecularExpert
    if literature is None:
        from .literature import default_registry
        literature = default_registry()
    molecular = MolecularExpert(config.molecular_mode)
    expected = {variant.id: variant for variant in case.variants}
    if len(candidates) != len(expected) or {c['variant_id'] for c in candidates} != set(expected):
        raise GuardError('Candidate set does not match input variants')
    previous = None
    seen_evidence = {}
    for position, candidate in enumerate(candidates, 1):
        variant = expected[candidate['variant_id']]
        if candidate['assembly'] != variant.assembly or candidate['genotype'] != variant.genotype:
            raise GuardError('Candidate assembly/genotype changed')
        score = candidate['score']
        if not finite_score(score) or not 0 <= score <= 1:
            raise GuardError('Invalid ranking score')
        key = order_key(candidate)
        if previous is not None and (key[0] < previous[0] or key[0] == previous[0] and score > -previous[1] + 1e-12):
            raise GuardError('Unsorted ranking')
        if candidate['rank'] != position:
            raise GuardError('Incorrect rank number')
        previous = key
        if set(candidate['experts']) != set(config.weights):
            raise GuardError('Missing expert result')
        for assessment in candidate['experts'].values():
            value = assessment['score']
            if value is not None and (not finite_score(value) or not 0 <= value <= 1):
                raise GuardError('Invalid expert score')
        if abs(aggregate(candidate['experts'], config) - score) > 1e-12:
            raise GuardError('Aggregate does not match expert scores')
        annotation = annotations[candidate['variant_id']]
        absent_rule = (config.molecular_mode == 'context_v3' and annotation.get('status') == 'annotated'
                       and candidate['experts']['molecular'].get('policy') == 'v1.3')
        if annotation.get('gnomad_max_af') is None and not absent_rule:
            if candidate['experts']['molecular'].get('components', {}).get('frequency', 0) != 0:
                raise GuardError('Missing AF was converted into rarity evidence')
        gene = candidate['gene']
        if gene and gene not in {t.get('gene_symbol') for t in annotation.get('transcripts', [])}:
            raise GuardError('Gene does not come from variant annotation')
        allowed_sources = {f'VEP:{variant.cache_key}'}
        if config.molecular_mode in CONTEXT_MODES:
            try:
                recalculated = molecular.assess(variant, annotation, gene, candidate['disease_id'], literature)
            except (ValueError, KeyError, TypeError) as error:
                raise GuardError('Invalid molecular evidence input') from error
            if recalculated != candidate['experts']['molecular']:
                raise GuardError('Molecular evidence applicability or score changed')
            allowed_sources.update(item['source'] for item in recalculated['evidence'])
        matched = literature.match(variant, gene, candidate['disease_id'])
        allowed_sources.update(literature.sources(matched))
        exact = [entry for entry in matched if entry['evidence_scope'] == 'exact_variant']
        if exact and candidate['experts']['molecular'].get('literature_applied'):
            result = candidate['experts']['molecular']
            base = result.get('literature_base_consequence')
            support = max(entry['support'] for entry in exact)
            consequence = result.get('components', {}).get('consequence')
            if base is None or consequence is None or abs(consequence - max(base, support)) > 1e-12:
                raise GuardError('Literature support was summed or did not use max()')
            if not set('LIT:' + entry['id'] for entry in exact).issubset(
                    {item['source'] for item in result['evidence']}):
                raise GuardError('Literature evidence is not attributable to the registry')
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
    checks = [
        'input_variant_identity', 'assembly_and_genotype', 'finite_scores', 'score_recalculation',
        'rank_order', 'missing_frequency', 'gene_provenance', 'evidence_references']
    if config.molecular_mode in CONTEXT_MODES:
        checks.extend(['contextual_molecular_recalculation', 'evidence_applicability'])
    checks.append('literature_max_not_sum')
    return {'status': 'passed', 'checks': checks,
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
        mandatory = {ev['id'] for ev in candidate['evidence'] if ev['criterion'] in REQUIRED_CRITERIA}
        if not mandatory.issubset(selected):
            raise GuardError('Presentation omitted conflicting or qualified evidence')
    return {'status': 'passed', 'checks': ['case_status', 'candidate_order', 'evidence_ids', 'required_caveats'],
            'scope': 'Facts rendered verbatim from evidence; LLM cannot add free-form claims'}
