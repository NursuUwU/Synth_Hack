"""Fixed orchestrator. No access to evaluation labels."""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from .guards import check_ranking
from .models import Case, utc_now
from .ranking import decide, rank, uncertainty_reasons
from .reporting import get_reporter, html_report, make_presentation, markdown


class Pipeline:
    def __init__(self, provider, knowledge, config, llm=False, literature=None):
        self.provider, self.knowledge, self.config, self.llm = provider, knowledge, config, llm
        if literature is None:
            from .literature import default_registry
            literature = default_registry()
        self.literature = literature
        self.reporter = get_reporter() if llm else None

    def run(self, cases: list[Case], offline=True, progress=None) -> list[dict]:
        annotations, source = self.provider.annotate([v for case in cases for v in case.variants],
                                                      offline=offline, progress=progress)
        reports = []
        for case in cases:
            literature = self.literature if self.config.molecular_mode == 'context_v2' else None
            candidates = rank(case, annotations, self.knowledge, self.config, literature)
            guard = check_ranking(case, candidates, annotations, self.config, literature)
            count = sum(annotations[v.id]['status'] == 'annotated' for v in case.variants)
            coverage = count / len(case.variants)
            knowledge_ready = bool(self.knowledge.diseases) and any(
                self.knowledge.resolve(term) in self.knowledge.parents and self.knowledge.ic(term) > 0
                for term in case.positive_hpo)
            status, prediction = decide(candidates, coverage, knowledge_ready, self.config)
            warnings = list(self.knowledge.warnings) + [
                'Synthetic composite genotype; not a full patient genome or clinical benchmark.',
                'Heuristic weights/thresholds are uncalibrated; ranking score is not a probability.',
                'Only GRCh38 SNVs are supported in this baseline.',
                'No family segregation, phase confirmation or independent reference FASTA validation.',
                'VEP clinical annotations do not include independently fetched condition-specific ClinVar review status.',
                'Absence of a convincing candidate does not prove absence of a genetic cause.',
                'Curated literature evidence is a small research sample; absence of a record is not evidence of benignity.',
            ]
            if case.metadata.get('dataset_modification'):
                warnings.append('Dataset modification: ' + str(case.metadata['dataset_modification'])
                                + '. This is not an organizer case and must not be used for benchmark scoring.')
            if coverage < 1:
                warnings.append(f'{len(case.variants) - count} variants lack annotations; missingness is not benign evidence.')
            if case.metadata.get('phenotype_truncated'):
                warnings.append('The organizer marked the phenotype profile as truncated.')
            if case.metadata.get('phenotype_override'):
                warnings.append('Phenotype was explicitly overridden for this run; not the original case profile.')
            unknown = sorted(term for term in case.positive_hpo | case.excluded_hpo
                             if self.knowledge.resolve(term) not in self.knowledge.parents)
            if unknown:
                warnings.append('HPO terms absent from the loaded ontology: ' + ', '.join(unknown))
            input_data = {'case_id': case.id, 'variants': [asdict(v) for v in case.variants],
                          'positive_hpo': case.positive_hpo, 'excluded_hpo': case.excluded_hpo,
                          'metadata': case.metadata, 'sex': case.sex}
            annotation_data = {v.id: annotations[v.id] for v in case.variants}
            genes = {}
            for candidate in candidates:
                if candidate['gene'] and candidate['gene'] not in genes:
                    genes[candidate['gene']] = {'gene': candidate['gene'], 'score': candidate['score'],
                        'best_variant': candidate['variant_id'], 'disease_id': candidate['disease_id']}
            report = {'schema_version': 1, 'created_at': utc_now(), 'case_id': case.id,
                'status': status, 'prediction': prediction,
                'uncertainty_reasons': uncertainty_reasons(candidates, self.config) if status == 'uncertain' else [],
                'variant_count': len(case.variants), 'annotation_count': count, 'annotation_coverage': coverage,
                'input': input_data,
                'input_sha256': hashlib.sha256(json.dumps(input_data, sort_keys=True).encode()).hexdigest(),
                'annotations': annotation_data,
                'annotation_sha256': hashlib.sha256(json.dumps(annotation_data, sort_keys=True).encode()).hexdigest(),
                'config': self.config.dictionary(), 'config_sha256': self.config.digest,
                'implementation_sha256': hashlib.sha256(b''.join(
                    path.name.encode() + b'\0' + path.read_bytes()
                    for path in sorted(Path(__file__).parent.glob('*.py')))).hexdigest(),
                'sources': {'vep': source, 'hpo': self.knowledge.manifest,
                            'literature': self.literature.metadata()},
                'candidates': candidates, 'genes': list(genes.values()), 'warnings': warnings,
                'ranking_guard': guard}
            report['presentation'] = make_presentation(report, use_llm=self.llm, reporter=self.reporter)
            if self.reporter and report['presentation']['provider'] != 'template':
                report['presentation']['model_metadata'] = dict(getattr(self.reporter, 'last_metadata', {}))
            if report['presentation']['warning']:
                report['warnings'].append(report['presentation']['warning'])
            if self.config.molecular_mode == 'context_v2':
                from .molecular import default_knowledge
                report['sources']['molecular'] = default_knowledge().metadata
                # Keep this limitation for VEP; independent ClinVar records are now separate.
                report['warnings'].append('ClinVar snapshot covers VEP clinical hits only; missing or '
                    'unmatched condition evidence supplies no positive support, not a benign classification.')
                report['warnings'].append('Literature curation was developed after inspection of open cases; '
                    'development improvement is not independent validation.')
            reports.append(report)
            if progress:
                progress(f'{case.id}: {status}; ranking and presentation guards passed')
        return reports


def write_run(reports: list[dict], directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    for report in reports:
        stem = directory / report['case_id']
        stem.with_suffix('.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        stem.with_suffix('.md').write_text(markdown(report), encoding='utf-8')
        stem.with_suffix('.html').write_text(html_report(report), encoding='utf-8')
        with (directory / f'{report["case_id"]}_ranking.tsv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['rank', 'variant_id', 'gene', 'disease_id', 'score', 'phenotype', 'molecular', 'inheritance'])
            for candidate in report['candidates']:
                writer.writerow([candidate['rank'], candidate['variant_id'], candidate['gene'], candidate['disease_id'],
                    candidate['score'], *[candidate['experts'][name]['score'] for name in ['phenotype', 'molecular', 'inheritance']]])
    complete = all(report['prediction'] is not None for report in reports)
    target = directory / ('predictions.tsv' if complete else 'predictions.partial.tsv')
    with target.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['case_id', 'prediction'])
        writer.writerows((r['case_id'], r['prediction']) for r in reports if r['prediction'] is not None)
    # Remove only our own stale prediction counterpart in the explicitly selected output directory.
    stale = directory / ('predictions.partial.tsv' if complete else 'predictions.tsv')
    if stale.is_file():
        stale.unlink()
    summary = {'created_at': utc_now(), 'case_ids': [r['case_id'] for r in reports],
        'predictions_complete': complete, 'prediction_file': target.name,
        'withheld_cases': [r['case_id'] for r in reports if r['prediction'] is None],
        'cases': [{'case_id': r['case_id'], 'status': r['status'], 'prediction': r['prediction'],
                   'annotation_coverage': r['annotation_coverage'], 'input_sha256': r['input_sha256']}
                  for r in reports]}
    (directory / 'run.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
