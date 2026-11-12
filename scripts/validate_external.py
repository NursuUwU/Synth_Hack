"""Pinned external, semi-synthetic pilot. Model rules remain unchanged."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shutil
import sys
import time
import urllib.parse
import zipfile
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from validate_controls import FrozenProvider, digest
from genomic_detective.annotations import AnnotationCache, VepProvider
from genomic_detective.clinvar import parse_vcv
from genomic_detective.dataset import Dataset
from genomic_detective.evaluation import read_development_labels
from genomic_detective.knowledge import Knowledge
from genomic_detective.models import Case, Variant, utc_now
from genomic_detective.molecular import MolecularKnowledge, default_knowledge
from genomic_detective.network import HttpClient
from genomic_detective.pipeline import Pipeline
from genomic_detective.ranking import RankingConfig

RELEASE = '0.1.27'
ARCHIVE_SHA = 'd0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a'
SEED = 'external-pilot-v1-2026-11-12'
SIZE = 20
BACKGROUND_SIZE = 399


def unpack_case(data):
    solved = [i for i in data.get('interpretations', []) if i.get('progressStatus') == 'SOLVED']
    causal = [g for i in solved for g in i.get('diagnosis', {}).get('genomicInterpretations', [])
              if g.get('interpretationStatus') == 'CAUSATIVE']
    if len(causal) != 1:
        raise ValueError('not_exactly_one_causative_interpretation')
    descriptor = causal[0].get('variantInterpretation', {}).get('variationDescriptor', {})
    vcf = descriptor.get('vcfRecord', {})
    if vcf.get('genomeAssembly') not in {'GRCh38', 'hg38'}:
        raise ValueError('no_explicit_GRCh38_vcf_record')
    state = descriptor.get('allelicState', {}).get('id')
    genotype = {'GENO:0000135': '0/1', 'GENO:0000136': '1/1'}.get(state)
    if not genotype:
        raise ValueError('unsupported_or_unknown_allelic_state')
    try:
        variant = Variant(vcf['chrom'], int(vcf['pos']), vcf['ref'], vcf['alt'], genotype)
    except (ValueError, KeyError, TypeError):
        raise ValueError('unsupported_or_invalid_SNV') from None
    if variant.chrom in {'chrX', 'chrY', 'chrMT'}:
        raise ValueError('non_autosomal_pilot_scope')
    gene = descriptor.get('geneContext', {}).get('symbol')
    if not gene:
        raise ValueError('no_causal_gene')
    pmids = sorted({r['id'] for r in data.get('metaData', {}).get('externalReferences', [])
                    if r.get('id', '').startswith('PMID:')})
    if not pmids:
        raise ValueError('no_publication_identifier')
    positive, excluded = {}, {}
    for feature in data.get('phenotypicFeatures', []):
        term = feature.get('type', {}).get('id', '')
        if not re.fullmatch(r'HP:\d{7}', term):
            raise ValueError('invalid_HPO')
        target = excluded if feature.get('excluded', False) else positive
        if term in positive or term in excluded:
            raise ValueError('duplicate_or_contradictory_HPO')
        target[term] = feature['type'].get('label', term)
    if not positive:
        raise ValueError('no_positive_HPO')
    return variant, gene, pmids, positive, excluded


def prepare(archive, output):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA:
        raise ValueError('Phenopacket archive checksum mismatch')
    original = Dataset(ROOT / 'Project_3_genomic_detecive_for_students.zip')
    try:
        cases = original.load_cases()
    finally:
        original.close()
    labels = read_development_labels(ROOT / 'Project_3_genomic_detecive_for_students.zip')
    known_variants = {r['expected_answer'] for r in labels}
    store = default_knowledge()
    known_genes = {m['gene'] for m in store.rules['mechanisms']}
    text = json.dumps([c.metadata for c in cases]) + json.dumps(labels) + json.dumps(store.rules)
    known_pmids = {'PMID:' + p for p in re.findall(r'(?:PMID:|pubmed.ncbi.nlm.nih.gov/|"pmid":\s*")(\d+)', text)}
    output.mkdir(parents=True, exist_ok=False)
    selected, audit, genes, publications, alleles = [], [], set(), set(), set()
    with zipfile.ZipFile(archive) as source:
        paths = sorted([p for p in source.namelist() if p.endswith('.json')],
                       key=lambda p: hashlib.sha256((SEED + p).encode()).hexdigest())
        for path in paths:
            row = {'source_path': path}
            try:
                payload = source.read(path)
                data = json.loads(payload)
                variant, gene, pmids, positive, excluded = unpack_case(data)
                if gene in known_genes or variant.id in known_variants or set(pmids) & known_pmids:
                    raise ValueError('overlap_with_development_or_curated_registry')
                if gene in genes or variant.id in alleles or set(pmids) & publications:
                    raise ValueError('repeat_gene_allele_or_publication_in_pilot')
                if len(selected) >= SIZE:
                    raise ValueError('eligible_after_prespecified_sample_limit')
                record = {'case_id': f'external_{len(selected) + 1:02}', 'source_path': path,
                    'source_sha256': hashlib.sha256(payload).hexdigest(), 'source_id': data.get('id'),
                    'expected_variant': variant.id, 'expected_gene': gene, 'pmids': pmids,
                    'disease': [i.get('diagnosis', {}).get('disease') for i in data.get('interpretations', [])],
                    'variant': asdict(variant), 'positive_hpo': positive, 'excluded_hpo': excluded,
                    'sex': data.get('subject', {}).get('sex', 'UNKNOWN_SEX')}
                selected.append(record)
                genes.add(gene); publications.update(pmids); alleles.add(variant.id)
                row['status'] = 'selected'
            except (ValueError, KeyError, TypeError) as error:
                row.update(status='excluded', reason=str(error))
            audit.append(row)
    if len(selected) != SIZE:
        raise ValueError(f'Only {len(selected)} eligible cases; no adaptive relaxation permitted')
    unique = {v.id: v for c in reversed(cases) for v in c.variants}
    pool = [v for key, v in unique.items() if key not in known_variants | alleles]
    background = sorted(pool, key=lambda v: hashlib.sha256((SEED + v.id).encode()).hexdigest())[:BACKGROUND_SIZE]
    inputs = []
    for record in selected:
        inputs.append(asdict(Case(record['case_id'], background + [Variant(**record['variant'])],
            record['positive_hpo'], record['excluded_hpo'],
            {'validation_scope': 'external_semi_synthetic', 'phenopacket_release': RELEASE}, record['sex'])))
    manifest = {'schema_version': 1, 'prepared_at': utc_now(), 'seed': SEED,
        'frozen_model_commit': 'd1df5c05c1796482c65306d6ede7745b9a877fb5',
        'scope': 'External published-case pilot with shared synthetic background; not clinical validation',
        'release': RELEASE, 'archive_sha256': ARCHIVE_SHA,
        'archive_url': f'https://github.com/monarch-initiative/phenopacket-store/releases/download/{RELEASE}/all_phenopackets.zip',
        'sample_size': SIZE, 'background_variants': BACKGROUND_SIZE, 'total_variants_per_case': BACKGROUND_SIZE + 1,
        'excluded_development_genes': sorted(known_genes), 'excluded_pmids': sorted(known_pmids),
        'inputs_sha256': digest(inputs), 'labels_sha256': digest(selected),
        'selection_rule': 'SHA256(seed+archive path); first 20 eligible distinct genes/publications/alleles; no ranking-based selection',
        'background_rule': '399 SHA256(seed+variant ID) selected original background SNVs; six open causal alleles and pilot targets removed',
        'audit_counts': dict(Counter(r.get('reason', r['status']) for r in audit)),
        'limitations': ['Public interpretations are accepted from the source; not independently clinically re-curated.',
                       'HPO/ClinVar can contain evidence from these publications; no temporal or knowledge-source holdout.',
                       'Shared synthetic background is not a full patient exome and limits generalization.',
                       'Autosomal single-SNV cases only; unsupported cases are explicitly excluded.',
                       'No score thresholds, weights or curated rules are tuned on pilot outcomes.']}
    for name, obj in [('manifest.json', manifest), ('inputs.json', inputs), ('labels.json', selected), ('selection-audit.json', audit)]:
        (output / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'Frozen {len(inputs)} cases with {BACKGROUND_SIZE} background SNVs each.', flush=True)


def fetch_clinical(annotations, output, offline):
    target = output / 'clinical-snapshot.json'
    ann_sha = digest(annotations)
    if target.exists():
        snapshot = json.loads(target.read_text(encoding='utf-8'))
        if snapshot['annotation_sha256'] != ann_sha:
            raise ValueError('External clinical snapshot annotation mismatch')
        return snapshot
    if offline:
        raise ValueError('No frozen external clinical snapshot for offline replay')
    requested = {k: a for k, a in annotations.items() if a.get('clinical')}
    rsids = sorted({c['record_id'] for a in requested.values() for c in a['clinical']
                    if c.get('record_id', '').startswith('rs') and c['record_id'][2:].isdigit()})
    client = HttpClient(timeout=60)
    found, sources = {}, []
    for start in range(0, len(rsids), 20):
        base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        query = urllib.parse.urlencode({'db': 'clinvar', 'term': ' OR '.join(rsids[start:start + 20]), 'retmode': 'json', 'retmax': 500})
        search = client.json(base + 'esearch.fcgi?' + query)['esearchresult']
        if int(search['count']) > len(search['idlist']):
            raise ValueError('Truncated clinical query')
        time.sleep(0.4)
        if not search['idlist']:
            continue
        url = base + 'efetch.fcgi?' + urllib.parse.urlencode({'db': 'clinvar', 'rettype': 'vcv', 'is_variationid': 'true', 'id': ','.join(search['idlist'])})
        payload = client.request(url, headers={'Accept': 'application/xml'})
        sha = hashlib.sha256(payload).hexdigest()
        (output / f'clinvar-{sha}.xml').write_bytes(payload)
        sources.append({'url': url, 'sha256': sha, 'bytes': len(payload)})
        for record in parse_vcv(payload):
            if record['variant_id'] in requested:
                found.setdefault(record['variant_id'], {})[record['accession']] = record
        time.sleep(0.4)
    snapshot = {'schema_version': 1, 'assembly': 'GRCh38', 'fetched_at': utc_now(),
        'annotation_sha256': ann_sha, 'scope': 'Uniform lookup of ALL VEP clinical hits across background and pilot targets; no label filtering',
        'requested_variants': sorted(requested), 'requested_rsids': rsids,
        'unmatched_variants': sorted(set(requested) - set(found)), 'sources': sources,
        'records': {k: list(v.values()) for k, v in sorted(found.items())}}
    target.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'ClinVar exact allele matches: {len(found)}/{len(requested)} VEP clinical hits.', flush=True)
    return snapshot


def metrics(rows):
    n = len(rows)
    completed = [r for r in rows if r['status'] not in {'incomplete', 'error'}]
    return {'attempted': n, 'completed': len(completed),
        'statuses': dict(Counter(r['status'] for r in rows)),
        **{f'top_{k}_all_attempted': sum(r['expected_rank'] is not None and r['expected_rank'] <= k for r in completed) / n if n else None for k in [1, 3, 5]},
        'mrr_all_attempted': sum(1 / r['expected_rank'] for r in completed if r['expected_rank']) / n if n else None,
        'candidate_wrong_first': sum(r['status'] == 'candidate' and r['expected_rank'] != 1 for r in rows),
        'errors': [r['case_id'] for r in rows if r['status'] == 'error']}


def restore_bundle(output):
    info = json.loads((ROOT / 'resources/validation-pilot-v1.manifest.json').read_text(encoding='utf-8'))
    archive = ROOT / 'resources' / info['archive']
    if hashlib.sha256(archive.read_bytes()).hexdigest() != info['sha256']:
        raise ValueError('External validation bundle checksum mismatch')
    output.mkdir(parents=True, exist_ok=False)
    with zipfile.ZipFile(archive) as source:
        if set(source.namelist()) != set(info['files']):
            raise ValueError('Unexpected files in external validation bundle')
        for name, entry in info['files'].items():
            if Path(name).name != name:
                raise ValueError('Unsafe bundle member')
            data = source.read(name)
            if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise ValueError('External validation member checksum mismatch')
            (output / name).write_bytes(data)


def verify_reference(result, reference):
    """Regression replay of a published pilot, not another independent evaluation."""
    expected = json.loads(reference.read_text(encoding='utf-8'))
    key = lambda row: (row['model'], row['case_id'])
    actual = {key(row): row for row in result['cases']}
    saved = {key(row): row for row in expected['cases']}
    if set(actual) != set(saved):
        raise ValueError('Replay case/model set differs from frozen reference')
    for identifier, row in actual.items():
        for field in ['expected_rank', 'status', 'prediction', 'annotation_coverage', 'config_sha256']:
            if row.get(field) != saved[identifier].get(field):
                raise ValueError(f'Frozen replay changed: {identifier} / {field}')
    print('PASS: all frozen ranks/statuses/predictions/coverage reproduced.', flush=True)


def run(output, offline, reference=None):
    manifest = json.loads((output / 'manifest.json').read_text(encoding='utf-8'))
    inputs = json.loads((output / 'inputs.json').read_text(encoding='utf-8'))
    labels = json.loads((output / 'labels.json').read_text(encoding='utf-8'))
    if digest(inputs) != manifest['inputs_sha256'] or digest(labels) != manifest['labels_sha256']:
        raise ValueError('Frozen external input/label checksum mismatch')
    cases = [Case(d['id'], [Variant(**v) for v in d['variants']], d['positive_hpo'], d['excluded_hpo'], d['metadata'], d['sex']) for d in inputs]
    ann_file = output / 'annotations.json'
    if ann_file.exists():
        saved = json.loads(ann_file.read_text(encoding='utf-8'))
        annotations, source = saved['annotations'], saved['source']
    else:
        cache = AnnotationCache(ROOT / 'data/cache/vep.sqlite3')
        try:
            annotations, source = VepProvider(cache).annotate([v for c in cases for v in c.variants], offline=offline, progress=print)
        finally:
            cache.close()
        ann_file.write_text(json.dumps({'annotations': annotations, 'source': source}, ensure_ascii=False) + '\n', encoding='utf-8')
    clinical = fetch_clinical(annotations, output, offline)
    molecular = MolecularKnowledge(copy.deepcopy(default_knowledge().rules), clinical)
    knowledge = Knowledge.load(ROOT / 'data/knowledge')
    base = RankingConfig.load(ROOT / 'config/ranking.json')
    configs = {'v1.1': RankingConfig.load(ROOT / 'config/ranking-v1.1.json'), 'v1.2': base,
               'phenotype_only': replace(base, weights={'phenotype': 1., 'molecular': 0., 'inheritance': 0.})}
    # Labels are consulted below only AFTER all model rankings have been computed.
    reports = {}
    with patch('genomic_detective.molecular.default_knowledge', return_value=molecular):
        for name, config in configs.items():
            reports[name] = {}
            for case in cases:
                try:
                    report = Pipeline(FrozenProvider(annotations, source), knowledge, config).run([case])[0]
                    reports[name][case.id] = report
                    print(f'{name}: {case.id} {report["status"]}', flush=True)
                except Exception as error:
                    reports[name][case.id] = {'status': 'error', 'error': f'{type(error).__name__}: {error}'}
    answers = {d['case_id']: d for d in labels}
    publication_overlap = {}
    for cid, answer in answers.items():
        references = {ref for disease in knowledge.gene_diseases.get(answer['expected_gene'], set())
                      for ref in knowledge.diseases[disease]['references']}
        cited_pmids = {p for ref in references for p in re.findall(r'PMID:\d+', ref)}
        publication_overlap[cid] = sorted(set(answer['pmids']) & cited_pmids)
    rows = []
    for name, results in reports.items():
        for cid, report in results.items():
            answer = answers[cid]
            row = {'model': name, 'case_id': cid, 'expected': answer['expected_variant'], 'gene': answer['expected_gene'],
                   'pmids': answer['pmids'], 'source_path': answer['source_path'], 'status': report['status'],
                   'hpo_gene_profile_publication_overlap': publication_overlap[cid]}
            if report['status'] == 'error':
                row.update(error=report['error'], expected_rank=None)
            else:
                target = next((c for c in report['candidates'] if c['variant_id'] == answer['expected_variant']), None)
                row.update(expected_rank=target['rank'] if target else None, prediction=report['prediction'],
                    target_gene=target['gene'] if target else None, target_disease=target['disease_id'] if target else None,
                    target_score=target['score'] if target else None,
                    target_tie_count=sum(abs(c['score'] - target['score']) <= 1e-12 for c in report['candidates']) if target else 0,
                    target_molecular=target['experts']['molecular'] if target else None,
                    top_three=[{k: c[k] for k in ['rank', 'variant_id', 'gene', 'disease_id', 'score']} for c in report['candidates'][:3]],
                    annotation_coverage=report['annotation_coverage'], uncertainty_reasons=report['uncertainty_reasons'],
                    ranking_guard=report['ranking_guard']['status'], presentation_guard=report['presentation']['guard']['status'],
                    input_sha256=report['input_sha256'], annotation_sha256=report['annotation_sha256'],
                    implementation_sha256=report['implementation_sha256'], config_sha256=report['config_sha256'])
            rows.append(row)
    result = {'created_at': utc_now(), 'manifest': manifest,
        'sources': {'vep': source, 'hpo': knowledge.manifest, 'molecular': molecular.metadata},
        'annotations_sha256': digest(annotations),
        'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'hpo_publication_overlap_cases': sum(bool(v) for v in publication_overlap.values()),
        'metrics': {name: metrics([r for r in rows if r['model'] == name]) for name in configs},
        'descriptive_no_direct_hpo_publication_overlap': {
            name: metrics([r for r in rows if r['model'] == name and not r['hpo_gene_profile_publication_overlap']]) for name in configs},
        'cases': rows, 'notes': ['phenotype_only compares RANKS; its decision statuses use combined-model gates and are not calibrated for this mode.',
                               'All attempted cases remain in main metric denominators; incomplete/error cases count as unsuccessful.']}
    (output / 'external-results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    if reference:
        verify_reference(result, reference)
    print(json.dumps(result['metrics'], indent=2), flush=True)
    return int(any(r['status'] == 'error' for r in rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/validation_external_v1')
    parser.add_argument('--archive', type=Path, default=ROOT / f'data/cache/validation/phenopacket-store-{RELEASE}.zip')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--bundled', action='store_true', help='Restore the included frozen pilot and recalculate offline')
    parser.add_argument('--expect', type=Path, help='Compare recalculated ranks/statuses with a saved pilot regression reference')
    parser.add_argument('--replay-from', type=Path, help='Copy frozen artifacts to a new output directory and replay offline')
    args = parser.parse_args()
    if args.bundled and args.replay_from:
        parser.error('Choose either --bundled or --replay-from')
    if args.bundled:
        restore_bundle(args.output)
        args.offline = True
    if args.replay_from:
        args.output.mkdir(parents=True, exist_ok=False)
        for name in ['manifest.json', 'inputs.json', 'labels.json', 'selection-audit.json', 'annotations.json', 'clinical-snapshot.json']:
            shutil.copy2(args.replay_from / name, args.output / name)
        args.offline = True
    if not args.output.exists():
        prepare(args.archive, args.output)
    if args.prepare_only:
        return 0
    if (args.output / 'external-results.json').exists():
        raise ValueError('Completed pilot is immutable; use a new output directory for replay')
    from genomic_detective.demo_data import ensure_demo_data
    ensure_demo_data(ROOT, progress=print)
    return run(args.output, args.offline, args.expect)


if __name__ == '__main__':
    sys.exit(main())
