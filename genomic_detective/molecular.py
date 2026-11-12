"""Versioned disease-context evidence. Internal support scores are not ACMG/LR."""
from __future__ import annotations

import copy
import hashlib
import json
from functools import lru_cache
from pathlib import Path

from .models import Variant, finite_score


DISRUPTIVE = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
              'splice_donor_variant', 'start_lost', 'stop_lost'}
ROOT = Path(__file__).resolve().parents[1]


class MolecularKnowledge:
    def __init__(self, rules, snapshot):
        if rules.get('schema_version') != 1 or snapshot.get('schema_version') != 1:
            raise ValueError('Unsupported molecular evidence schema')
        if snapshot.get('assembly') != 'GRCh38':
            raise ValueError('Molecular evidence requires GRCh38')
        self.rules, self.snapshot = rules, snapshot
        self.sources = {s['id']: s for s in rules['sources']}
        self.records = snapshot['records']
        from .clinvar import review_stars
        for key, records in self.records.items():
            for record in records:
                if record['variant_id'] != key or record['assembly'] != 'GRCh38':
                    raise ValueError('ClinVar snapshot allele/assembly mismatch')
                for assertion in record['assertions']:
                    if (assertion['source'] != 'ClinVar:' + assertion['accession']
                            or assertion['review_stars'] != review_stars(assertion['review_status'])
                            or not isinstance(assertion['compound_condition'], bool)):
                        raise ValueError('Invalid ClinVar assertion provenance/review status')
        self.by_variant = {}
        for record in rules['variant_evidence']:
            fields = record['variant_id'].split(':')
            if len(fields) != 4:
                raise ValueError('Invalid curated allele')
            variant = Variant(fields[0], int(fields[1]), fields[2], fields[3], '0/1', record['assembly'])
            if variant.id != record['variant_id'] or not record['disease_ids']:
                raise ValueError('Invalid curated applicability')
            if (not finite_score(record['support']) or not 0 <= record['support'] <= 1
                    or record['kind'] not in {'established_functional_effect', 'literature_context'}
                    or record['kind'] == 'literature_context' and record['support'] != 0
                    or not set(record['source_ids']) <= set(self.sources)):
                raise ValueError('Invalid curated support/source')
            self.by_variant.setdefault(variant.id, []).append(record)
        for mechanism in rules['mechanisms']:
            if mechanism['source_id'] not in self.sources:
                raise ValueError('Unknown mechanism source')
        self.metadata = {'version': rules['version'], 'clinical_snapshot_date': snapshot['fetched_at'],
            'clinical_lookup_scope': snapshot['scope'], 'clinical_variants': len(self.records),
            'curated_alleles': len(self.by_variant),
            'rules_sha256': hashlib.sha256(json.dumps(rules, sort_keys=True).encode()).hexdigest(),
            'clinical_snapshot_sha256': hashlib.sha256(json.dumps(snapshot, sort_keys=True).encode()).hexdigest(),
            'scoring': 'max(0.55 * prediction, established_function, condition_classification) + 0.15 * rarity',
            'scope': 'Heuristic prioritization; no calibrated probability, LR or automatic ACMG codes'}

    @classmethod
    def load(cls, rules_path, snapshot_path):
        return cls(json.loads(Path(rules_path).read_text(encoding='utf-8')),
                   json.loads(Path(snapshot_path).read_text(encoding='utf-8')))

    def extract(self, variant, annotation, gene, baseline):
        if (annotation.get('variant_id', variant.id) != variant.id
                or annotation.get('assembly', variant.assembly) != variant.assembly):
            raise ValueError('Molecular annotation allele/assembly mismatch')
        transcript = next((t for t in annotation.get('transcripts', [])
                           if t.get('transcript_id') == baseline['selected_transcript']), {})
        hgvs = (transcript.get('hgvsc') or '').partition(':')[2]
        curated = [r for r in self.by_variant.get(variant.id, [])
                   if r['assembly'] == variant.assembly and r['gene'] == gene
                   and r['refseq'] == transcript.get('mane_select') and r['hgvsc'] == hgvs]
        clinical = [r for r in self.records.get(variant.id, [])
                    if r['assembly'] == variant.assembly and gene in r['genes']]
        return {'baseline': baseline, 'curated': curated, 'clinical': clinical,
                'gene': gene, 'variant': variant, 'transcript': transcript}

    @staticmethod
    def context_match(disease_id, ids, assertions):
        if disease_id in ids:
            return True
        # Only one-condition records provide cross-references; compound lists cannot
        # make two different diseases aliases of one another.
        return any(not a['compound_condition'] and disease_id in a['condition_ids']
                   and set(ids).intersection(a['condition_ids']) for a in assertions)

    def assess(self, facts, disease_id):
        from .experts import clinical_support, evidence
        base = copy.deepcopy(facts['baseline'])
        if base['score'] is None:
            return base
        variant, gene = facts['variant'], facts['gene']
        subject = f'{variant.id}:{gene}:{disease_id}'
        # Keep original VEP consequence and frequency facts; aggregate clinical labels
        # remain visible, but cannot supply condition-specific positive support.
        items = base['evidence']
        for item in items:
            if item['criterion'] == 'clinical_annotation':
                item['text'] = item['text'].replace(
                    'Заболевание и статус экспертной проверки ClinVar отдельно не получены.',
                    'Агрегированная метка VEP не добавляет поддержку конкретной болезни без отдельного сопоставления ClinVar.')
                # Text changes require a new deterministic identity.
                replacement = evidence(subject, item['criterion'], item['text'], item['source'], item['basis'])
                item.update(replacement)
        missing = [text for text in base['missing'] if text not in {
            'Condition-specific ClinVar assertion and review status not retrieved',
            'Allele-specific clinical interpretation unavailable'}]
        assertions = [a for record in facts['clinical'] for a in record['assertions']]
        matched = [a for a in assertions if not a['compound_condition'] and disease_id in a['condition_ids']]
        reviewed = [a for a in matched if a['review_stars'] >= 1]
        status, clinical = clinical_support([a['classification'] for a in reviewed])
        categories = set()
        for assertion in reviewed:
            category, _ = clinical_support([assertion['classification']])
            category = {'likely_pathogenic': 'pathogenic', 'likely_benign': 'benign'}.get(category, category)
            if category in {'pathogenic', 'benign', 'uncertain_significance', 'conflicting'}:
                categories.add(category)
        if len(categories) > 1:
            status, clinical = 'conflicting', 0.0
        conflict = status == 'conflicting'
        for assertion in matched:
            text = (f"ClinVar {assertion['accession']}: {assertion['classification']}; "
                    f"заболевание: {'; '.join(assertion['conditions'])}; "
                    f"review status: {assertion['review_status']} ({assertion['review_stars']} зв.); "
                    f"оценка: {assertion['last_evaluated'] or 'дата неизвестна'}. "
                    'Число звёзд не является вероятностью патогенности.')
            if assertion['review_stars'] == 0:
                text += ' Без предоставленных критериев классификация не добавляет баллы в этой версии.'
            item = evidence(subject, 'clinical_assertion', text, assertion['source'], 'database')
            item['url'] = assertion['url']
            items.append(item)
        for record in facts['clinical']:
            skipped = [a for a in record['assertions'] if a not in matched]
            if skipped:
                source = 'ClinVar:' + record['accession']
                items.append(evidence(subject, 'clinical_context',
                    f'{len(skipped)} записей ClinVar для другого, составного или несопоставленного '
                    'заболевания не добавляют поддержку выбранной гипотезе.', source, 'database'))
        if not reviewed:
            missing.append('Reviewed condition-specific ClinVar classification unavailable; no positive support assigned')
        if not facts['clinical']:
            missing.append('Exact-allele ClinVar record unavailable in the bundled snapshot; not assumed benign')
        consequence = base['components']['consequence']
        mechanism_name = None
        if set(base['consequences']) & DISRUPTIVE:
            consequence = min(consequence, 0.55)
            missing.append('NMD, exon relevance and disease-specific PVS1 requirements not established')
        for mechanism in self.rules['mechanisms']:
            if mechanism['gene'] != gene or not self.context_match(disease_id, mechanism['disease_ids'], assertions):
                continue
            mechanism_name = mechanism['mechanism']
            item = evidence(subject, 'mechanism_context', mechanism['text'],
                            'Literature:' + self.source_key(mechanism['source_id']), 'literature_context')
            item['url'] = self.sources[mechanism['source_id']]['url']
            items.append(item)
            # A generic null consequence does not demonstrate an altered-product mechanism.
            if not mechanism['generic_null_support'] and set(base['consequences']) & {'stop_gained', 'frameshift_variant'}:
                consequence = 0.0
                missing.append('Generic null consequence does not establish the selected disease mechanism')
        functional = 0.0
        applied = []
        for record in facts['curated']:
            if not self.context_match(disease_id, record['disease_ids'], assertions):
                continue
            criterion = 'functional_effect' if record['kind'] == 'established_functional_effect' else 'literature_context'
            for source_id in record['source_ids']:
                item = evidence(subject, criterion, record['text'],
                    'Literature:' + self.source_key(source_id),
                    'experimental_literature' if criterion == 'functional_effect' else 'literature_context')
                item['url'] = self.sources[source_id]['url']
                items.append(item)
            functional = max(functional, record['support'])
            applied.append({'variant_id': record['variant_id'], 'refseq': record['refseq'],
                            'hgvsc': record['hgvsc'], 'kind': record['kind'], 'source_ids': record['source_ids']})
        contradiction = functional > 0 and status in {'benign', 'likely_benign', 'conflicting'}
        if conflict:
            items.append(evidence(subject, 'evidence_conflict',
                'Для выбранной болезни клинические классификации конфликтуют или расходятся '
                'между записями; положительная клиническая поддержка не назначена.',
                'Policy:molecular-context-v1.2', 'compatibility'))
        if contradiction:
            items.append(evidence(subject, 'evidence_conflict',
                'Функциональное литературное свидетельство и клиническая интерпретация требуют '
                'совместной проверки; функциональное влияние не равно патогенности.',
                'Policy:molecular-context-v1.2', 'compatibility'))
        frequency = base['components']['frequency']
        support_block = max(0.55 * consequence, functional, clinical)
        base.update(score=support_block + 0.15 * frequency,
            components={'consequence': consequence, 'functional': functional,
                        'clinical': clinical, 'frequency': frequency},
            support_block=support_block, scoring_policy='context_v2_max_not_sum',
            context_disease_id=disease_id, mechanism=mechanism_name,
            clinical_status=status, clinical_conflict=conflict,
            clinical_context_status='matched' if reviewed else 'unreviewed' if matched else 'unresolved',
            functional_clinical_conflict=contradiction,
            clinical_assertions=matched, curated_evidence=applied,
            evidence=items, missing=sorted(set(missing)))
        return base

    def source_key(self, source_id):
        source = self.sources[source_id]
        return 'PMID:' + source['pmid'] if source.get('pmid') else source_id


@lru_cache(maxsize=1)
def default_knowledge():
    return MolecularKnowledge.load(ROOT / 'config/molecular_knowledge.json',
                                   ROOT / 'resources/clinvar-snapshot.json')
