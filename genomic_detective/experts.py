"""Deterministic experts. Scores are heuristic support, never disease probabilities."""
from __future__ import annotations

import copy
import hashlib
import json

from .knowledge import Knowledge
from .models import CONTEXT_MODES, Case, Variant


CONSEQUENCE_SUPPORT = {
    'stop_gained': 0.9, 'frameshift_variant': 0.9,
    'splice_acceptor_variant': 0.85, 'splice_donor_variant': 0.85,
    'start_lost': 0.75, 'stop_lost': 0.7,
    'missense_variant': 0.55, 'protein_altering_variant': 0.5,
    'splice_region_variant': 0.35, 'synonymous_variant': 0.1,
    'intron_variant': 0.08, 'intergenic_variant': 0.02,
    '3_prime_UTR_variant': 0.1, '5_prime_UTR_variant': 0.1,
}


def evidence(subject: str, criterion: str, text: str, source: str, basis: str) -> dict:
    identity = f'{subject}|{criterion}|{text}|{source}'
    return {'id': 'ev_' + hashlib.sha256(identity.encode()).hexdigest()[:16],
            'criterion': criterion, 'text': text, 'source': source, 'basis': basis}


def _literature_url(entry: dict) -> str:
    source = entry.get('source') or {}
    if source.get('pmid'):
        return 'https://pubmed.ncbi.nlm.nih.gov/' + str(source['pmid']) + '/'
    return str(source.get('url') or '')


def clinical_support(significance: list[str]) -> tuple[str, float]:
    """Interpret supplied labels; never derive ACMG criteria from an aggregate label.

    VUS, conflicts and unrecognized classifications supply no positive pathogenicity
    support. Neutrality does not imply a benign classification.
    """
    terms = {part.strip().lower().replace('_', ' ')
             for value in significance for part in value.split('/')}
    pathogenic = bool(terms & {'pathogenic', 'likely pathogenic'})
    benign = bool(terms & {'benign', 'likely benign'})
    if pathogenic and benign or any('conflict' in term for term in terms):
        return 'conflicting', 0.0
    if pathogenic:
        return ('pathogenic' if 'pathogenic' in terms else 'likely_pathogenic'), 0.8
    if benign:
        return ('benign' if 'benign' in terms else 'likely_benign'), 0.0
    if 'uncertain significance' in terms:
        return 'uncertain_significance', 0.0
    return ('not_actionable' if terms else 'unavailable'), 0.0


class PhenotypeExpert:
    def __init__(self, knowledge: Knowledge):
        self.kb = knowledge
        self.cache = {}

    def assess(self, case: Case, gene: str) -> dict:
        return self.hypotheses(case, gene)[0]

    def hypotheses(self, case: Case, gene: str) -> list[dict]:
        # Cache includes profile content: changing HPO terms cannot reuse a stale score.
        key = (gene, tuple(sorted(case.positive_hpo)), tuple(sorted(case.excluded_hpo)))
        if key in self.cache:
            return self.cache[key]
        diseases = self.kb.gene_diseases.get(gene, set())
        if not diseases or not case.positive_hpo:
            return [{'score': None, 'disease_id': None, 'disease_name': None,
                    'evidence': [], 'missing': ['No disease associations or positive phenotype'],
                    'matches': [], 'conflicts': []}]
        possibilities = [self._disease(case, gene, self.kb.diseases[d]) for d in sorted(diseases)]
        possibilities.sort(key=lambda result: (-(result['score'] if result['score'] is not None else -1),
                                                result['disease_id']))
        self.cache[key] = possibilities
        return possibilities

    def _disease(self, case: Case, gene: str, disease: dict) -> dict:
        source = f"HPO:{disease['id']}"
        subject = f"{gene}:{disease['id']}"
        weights, total, matches, conflicts, items = 0.0, 0.0, [], [], []
        unknown = []
        for term, label in case.positive_hpo.items():
            resolved = self.kb.resolve(term)
            if resolved not in self.kb.parents or self.kb.ic(resolved) <= 0:
                unknown.append(term)
                continue
            weight = self.kb.ic(resolved)
            ancestors = self.kb.ancestors(resolved)
            common = ancestors & disease['expanded']
            shared = max(common, key=lambda t: (self.kb.ic(t), t)) if common else None
            similarity = min(1.0, self.kb.ic(shared) / weight) if shared else 0.0
            total += weight * similarity
            weights += weight
            if shared and similarity > 0:
                matches.append({'query': term, 'shared_ancestor': shared, 'similarity': round(similarity, 6)})
                items.append(evidence(subject, 'phenotype_match',
                    f'{label} ({term}): соответствие профилю {disease["id"]} '
                    f'через {self.kb.labels.get(shared, shared)} ({shared}).', source, 'database'))
            if any(negative in ancestors for negative in disease.get('negative', set())):
                conflicts.append({'term': term, 'kind': 'disease_excludes_observed_feature'})
        for term, label in case.excluded_hpo.items():
            resolved = self.kb.resolve(term)
            # A patient-level exclusion is a warning: variable penetrance prevents hard vetoes.
            if resolved in disease['expanded']:
                conflicts.append({'term': term, 'kind': 'patient_excludes_disease_feature'})
                items.append(evidence(subject, 'phenotype_conflict',
                    f'{label} ({term}) явно отсутствует в случае, но указан в профиле '
                    f'{disease["id"]}; возможна вариабельность проявлений.', source, 'database'))
        score = max(0.0, total / weights - min(0.2, 0.04 * len(conflicts))) if weights else None
        return {'score': score, 'disease_id': disease['id'], 'disease_name': disease['name'],
                'matches': matches, 'conflicts': conflicts, 'evidence': items,
                'missing': [f'HPO term unavailable or uninformative: {term}' for term in unknown]}


class MolecularExpert:
    def __init__(self, mode='legacy_v1', knowledge=None):
        if mode not in {'legacy_v1'} | CONTEXT_MODES:
            raise ValueError('Unknown molecular mode')
        self.mode = mode
        if mode in CONTEXT_MODES and knowledge is None:
            from .molecular import default_knowledge
            knowledge = default_knowledge()
        self.knowledge = knowledge
        self.cache = {}
        self.facts_cache = {}

    def extract_facts(self, variant, annotation, gene):
        baseline = self._baseline(variant, annotation, gene)
        if self.mode in CONTEXT_MODES:
            return self.knowledge.extract(variant, annotation, gene, baseline)
        return {'baseline': baseline}

    def assess_facts(self, facts, disease_id=None):
        if self.mode not in CONTEXT_MODES:
            return facts['baseline']
        return self.knowledge.assess(facts, disease_id, v3=self.mode == 'context_v3')

    def extract_facts_cached(self, variant, annotation, gene):
        identity = hashlib.sha256(json.dumps(annotation, sort_keys=True, default=str).encode()).hexdigest()
        key = (variant.id, gene, identity)
        if key not in self.facts_cache:
            self.facts_cache[key] = self.extract_facts(variant, annotation, gene)
        return self.facts_cache[key]

    def assess(self, variant: Variant, annotation: dict, gene: str | None, disease_id=None,
               literature=None) -> dict:
        # The curated layer belongs to v1.2. It is never auto-loaded for the
        # frozen legacy mode; explicit registries remain supported for tests.
        if literature is None and self.mode in CONTEXT_MODES:
            from .literature import default_registry
            literature = default_registry()
        identity = hashlib.sha256(json.dumps(annotation, sort_keys=True, default=str).encode()).hexdigest()
        key = (variant.id, gene, disease_id, identity,
               literature.digest if literature is not None else None)
        if key in self.cache:
            return copy.deepcopy(self.cache[key])
        result = self.assess_facts(self.extract_facts_cached(variant, annotation, gene), disease_id)
        if literature is not None:
            result = self.apply_literature(result, variant, gene, disease_id, literature)
        self.cache[key] = result
        return copy.deepcopy(result)

    def apply_literature(self, result: dict, variant: Variant, gene: str | None, disease_id,
                         literature) -> dict:
        """Exact-allele curation raises consequence support with max(), never by summing."""
        if result.get('score') is None or not gene or not disease_id:
            return result
        matched = literature.match(variant, gene, disease_id)
        if not matched:
            return result
        subject = f'{variant.id}:{gene}:{disease_id}'
        exact = [entry for entry in matched if entry['evidence_scope'] == 'exact_variant']
        context = [entry for entry in matched if entry['evidence_scope'] != 'exact_variant']
        items = result['evidence']
        for entry in context:
            text = (f"{entry['fact']} Область применимости: {entry['evidence_scope']}; "
                    'свидетельство показано как контекст и не добавляет поддержку точному аллелю.')
            item = evidence(subject, 'literature_context', text,
                                 'LIT:' + entry['id'], 'literature')
            item['url'] = _literature_url(entry)
            items.append(item)
        if not exact:
            result['literature_context'] = [entry['id'] for entry in context]
            return result
        support = max(entry['support'] for entry in exact)
        base = result['components']['consequence']
        result['literature_base_consequence'] = base
        result['literature_support'] = support
        result['components']['consequence'] = max(base, support)
        for entry in exact:
            text = (f"{entry['fact']} Источник: PMID {entry['source'].get('pmid', 'не указан')}, "
                    f"{entry['source'].get('material', 'материал не указан')}. "
                    'Это установленный функциональный эффект для точного аллеля, '
                    'не ACMG-класс и не вероятность диагноза.')
            boundaries = entry.get('not_established') or []
            if boundaries:
                text += ' Границы: ' + '; '.join(boundaries) + '.'
            item = evidence(subject, 'literature_functional_effect', text,
                                 'LIT:' + entry['id'], 'literature')
            item['url'] = _literature_url(entry)
            items.append(item)
        result['literature_applied'] = [entry['id'] for entry in exact]
        if result.get('scoring_policy') == 'context_v2_max_not_sum':
            components = result['components']
            support_block = max(0.55 * components['consequence'], components.get('functional', 0.0),
                                components.get('clinical', 0.0))
            result['support_block'] = support_block
            result['score'] = support_block + 0.15 * components['frequency']
        else:
            components = result['components']
            result['score'] = (0.55 * components['consequence'] + 0.15 * components['frequency']
                               + 0.30 * components.get('clinical', 0.0))
        return result

    def _baseline(self, variant: Variant, annotation: dict, gene: str | None) -> dict:
        if annotation.get('status') != 'annotated':
            return {'score': None, 'evidence': [], 'missing': ['Variant annotation unavailable'],
                    'selected_transcript': None, 'consequences': [], 'components': {}}
        subject = f'{variant.id}:{gene or "intergenic"}'
        source = f'VEP:{variant.cache_key}'
        items, missing = [], []
        transcripts = [t for t in annotation['transcripts'] if t.get('gene_symbol') == gene] if gene else []
        overlapping = [t for t in transcripts if set(t.get('consequence_terms') or [])
                       - {'upstream_gene_variant', 'downstream_gene_variant'}]
        relevant = overlapping or transcripts
        preferred = [t for t in relevant if t.get('mane_select')]
        selection = 'overlapping MANE Select'
        if not preferred:
            preferred = [t for t in relevant if t.get('mane_plus_clinical')]
            selection = 'overlapping MANE Plus Clinical'
        if not preferred:
            preferred = [t for t in relevant if t.get('canonical') and t.get('biotype') == 'protein_coding']
            selection = 'canonical protein-coding'
        if not preferred:
            preferred = [t for t in relevant if t.get('biotype') == 'protein_coding'] or relevant
            selection = 'available transcript'
        selected = max(preferred, key=lambda t: max(
            (CONSEQUENCE_SUPPORT.get(c, 0.15) for c in t.get('consequence_terms') or []), default=0.0)) if preferred else None
        consequences = selected.get('consequence_terms') or [] if selected else [annotation.get('most_severe_consequence') or 'unknown']
        consequence_score = max((CONSEQUENCE_SUPPORT.get(c, 0.15) for c in consequences), default=0.0)
        items.append(evidence(subject, 'consequence',
            f'Предполагаемое последствие: {", ".join(consequences)}; '
            f'транскрипт: {selected.get("transcript_id") if selected else "не определён"} '
            f'({selection}). Само по себе не устанавливает патогенность.', source, 'prediction'))
        frequency = annotation.get('gnomad_max_af')
        frequency_score = 0.0
        if frequency is None:
            missing.append('gnomAD frequency unavailable; not assumed zero')
        else:
            frequency_score = 1.0 if frequency < 0.0001 else 0.5 if frequency < 0.001 else 0.0
            items.append(evidence(subject, 'population_frequency',
                f'Максимальная доступная частота gnomAD: {frequency:.6g}. '
                'Редкость не является самостоятельным доказательством патогенности.', source, 'database'))
        significance = sorted({value.lower().replace('_', ' ')
                               for record in annotation.get('clinical', []) for value in record['significance']})
        clinical_status, clinical_score = clinical_support(significance)
        conflict = clinical_status == 'conflicting'
        if significance:
            neutrality = (' VUS, конфликт или непатогенная/неподдерживаемая категория '
                          'не добавляют положительной поддержки патогенности.'
                          if clinical_score == 0 else '')
            items.append(evidence(subject, 'clinical_annotation',
                f'Клиническая аннотация Ensembl: {", ".join(significance)}. '
                'Заболевание и статус экспертной проверки ClinVar отдельно не получены.'
                + neutrality, source, 'database'))
            missing.append('Condition-specific ClinVar assertion and review status not retrieved')
        else:
            missing.append('Allele-specific clinical interpretation unavailable')
        components = {'consequence': consequence_score, 'frequency': frequency_score,
                      'clinical': clinical_score}
        # Missing evidence contributes no support. This is not a benign classification.
        score = 0.55 * consequence_score + 0.15 * frequency_score + 0.30 * clinical_score
        return {'score': score, 'evidence': items, 'missing': missing,
                'selected_transcript': selected.get('transcript_id') if selected else None,
                'transcript_selection': selection, 'consequences': consequences,
                'components': components, 'clinical_conflict': bool(conflict),
                'clinical_status': clinical_status}


# Exomiser default maximum allele frequencies per mode of inheritance
# (Jacobsen et al. 2022, doi:10.1093/bib/bbac188): 0.1% dominant or
# homozygous recessive, 2% compound-heterozygous recessive.
MOI_MAX_AF = {'AD': 0.001, 'AR_hom': 0.001, 'AR_comphet': 0.02}


class InheritanceExpert:
    def __init__(self, knowledge: Knowledge):
        self.kb = knowledge

    def assess(self, variant: Variant, disease_id: str | None, moi_filter: bool = False,
               allele_frequency: float | None = None, second_allele: bool = False,
               clinvar_whitelisted: bool = False, stated_mode: str | None = None) -> dict:
        disease = self.kb.diseases.get(disease_id)
        modes = sorted(disease['inheritance']) if disease else []
        if moi_filter and stated_mode and not modes and disease:
            # A disease record without modes must not bypass the stated mode: use the modes
            # described for the same gene's other diseases (e.g. Orphanet vs OMIM records).
            gene_modes = sorted({mode for gene in disease.get('genes', ())
                                 for other in self.kb.gene_diseases.get(gene, ())
                                 for mode in self.kb.diseases[other]['inheritance']} & {'AD', 'AR'})
            if gene_modes and stated_mode not in gene_modes:
                text = (f'Для этой болезни тип наследования не указан; для других болезней гена описан только '
                        f'{", ".join(gene_modes)}. Заданный тип {stated_mode} с ним не совпадает.')
                return {'score': 0.0, 'modes': [], 'gene_modes': gene_modes, 'moi_compatible': False,
                        'stated_mode': stated_mode, 'evidence': [evidence(
                            variant.id + ':' + str(disease_id), 'inheritance', text,
                            f'HPO:{disease_id}', 'compatibility')],
                        'missing': ['Inheritance mode of this disease record unavailable']}
        if not modes or variant.genotype in {'./.', '.|.'}:
            return {'score': None, 'modes': modes, 'evidence': [],
                    'missing': ['Inheritance mode or genotype unavailable']}
        if moi_filter:
            return self._filtered(variant, disease_id, modes, allele_frequency, second_allele,
                                  clinvar_whitelisted, stated_mode)
        scores, explanations = [], []
        if 'AD' in modes:
            scores.append(1.0 if variant.heterozygous else 0.5)
            explanations.append('Генотип совместим с гетерозиготной моделью AD.' if variant.heterozygous
                                else 'Гомозиготный генотип требует отдельной проверки для AD-модели.')
        if 'AR' in modes:
            scores.append(1.0 if variant.homozygous_alt else 0.2)
            explanations.append('Гомозиготный альтернативный генотип совместим с AR-моделью.'
                                if variant.homozygous_alt else
                                'Одна гетерозиготная запись не подтверждает AR-модель; второй аллель и фаза не установлены.')
        if not scores:
            return {'score': None, 'modes': modes, 'evidence': [],
                    'missing': ['Inheritance model is outside the autosomal baseline']}
        text = ' '.join(explanations) + ' Семейная сегрегация и происхождение de novo не подтверждены.'
        return {'score': max(scores), 'modes': modes, 'evidence': [evidence(
            variant.id + ':' + str(disease_id), 'inheritance', text,
            f'HPO:{disease_id}', 'compatibility')],
            'missing': ['Family segregation and phase unavailable']}

    def _filtered(self, variant, disease_id, modes, af, second_allele, whitelisted=False, stated=None):
        """v1.3: Exomiser-style MOI compatibility, including allele frequency per mode."""
        def rare(limit):
            # Unknown or absent frequency is not treated as common. Well-supported ClinVar
            # P/LP alleles bypass frequency limits, as in Exomiser (Jacobsen et al. 2022),
            # but still need a genotype compatible with the mode.
            return af is None or af <= limit or whitelisted
        options, explanations = [], []
        if stated and {'AD', 'AR'} & set(modes) and stated not in modes:
            # The clinician's stated mode excludes hypotheses described only with other modes.
            text = (f'Заданный тип наследования {stated} не описан для этой болезни '
                    f'({", ".join(modes)}). Гипотеза не может быть выдана как прогноз.')
            return {'score': 0.0, 'modes': modes, 'moi_compatible': False, 'stated_mode': stated,
                    'evidence': [evidence(variant.id + ':' + str(disease_id), 'inheritance', text,
                                          f'HPO:{disease_id}', 'compatibility')],
                    'missing': ['Family segregation and phase unavailable']}
        if stated:
            modes = [stated]
        if 'AD' in modes:
            ok = rare(MOI_MAX_AF['AD'])
            options.append(((1.0 if variant.heterozygous else 0.5) if ok else 0.0, ok))
            explanations.append(
                ('Генотип совместим с AD.' if variant.heterozygous else 'Гомозигота для AD требует проверки.')
                if ok else f'Частота {af:.4g} выше порога {MOI_MAX_AF["AD"]} для AD.')
        if 'AR' in modes:
            if variant.homozygous_alt:
                ok = rare(MOI_MAX_AF['AR_hom'])
                options.append((1.0 if ok else 0.0, ok))
                explanations.append('Гомозигота совместима с AR.' if ok else
                                    f'Частота {af:.4g} выше порога {MOI_MAX_AF["AR_hom"]} для гомозиготы AR.')
            elif second_allele and rare(MOI_MAX_AF['AR_comphet']):
                options.append((0.6, True))
                explanations.append('В гене есть второй редкий вариант: компаунд-гетерозиготность AR возможна, '
                                    'но фаза (trans) не установлена.')
            else:
                options.append((0.0, False))
                explanations.append('Одиночная гетерозигота без второго редкого аллеля в гене не объясняет AR-болезнь.')
        if not options:
            return {'score': None, 'modes': modes, 'evidence': [], 'moi_compatible': True,
                    'missing': ['Inheritance model is outside the autosomal baseline']}
        compatible = any(ok for _, ok in options)
        score = max(value for value, ok in options if ok) if compatible else 0.0
        text = ' '.join(explanations) + ' Семейная сегрегация и происхождение de novo не подтверждены.'
        if whitelisted and af is not None and af > MOI_MAX_AF['AD']:
            text += (f' Частота {af:.4g} не исключает гипотезу: аллель имеет сопоставленную запись ClinVar '
                     'P/LP с ≥2 звёздами (исключение ClinVar whitelist, как в Exomiser).')
        if not compatible:
            text += ' Гипотеза не проходит фильтр совместимости с типом наследования и не может быть выдана как прогноз.'
        return {'score': score, 'modes': modes, 'moi_compatible': compatible, 'evidence': [evidence(
            variant.id + ':' + str(disease_id), 'inheritance', text, f'HPO:{disease_id}', 'compatibility')],
            'missing': ['Family segregation and phase unavailable']}
