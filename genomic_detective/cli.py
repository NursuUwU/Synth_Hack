from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .annotations import AnnotationCache, VepProvider
from .dataset import Dataset
from .knowledge import Knowledge, fetch_knowledge


DEFAULT_ARCHIVE = 'Project_3_genomic_detecive_for_students.zip'


def parser():
    root = argparse.ArgumentParser(description='Genomic Detective: educational variant prioritization')
    commands = root.add_subparsers(dest='command', required=True)
    for name in ['prepare', 'run', 'benchmark', 'serve']:
        sub = commands.add_parser(name)
        sub.add_argument('--dataset', type=Path, default=Path(DEFAULT_ARCHIVE))
        sub.add_argument('--cache', type=Path, default=Path('data/cache/vep.sqlite3'))
        sub.add_argument('--knowledge', type=Path, default=Path('data/knowledge'))
        if name != 'serve':
            sub.add_argument('--split', choices=['all', 'development', 'hidden'], default='all')
            sub.add_argument('--case', action='append', dest='cases')
    prep = commands.choices['prepare']
    prep.add_argument('--batch-size', type=int, default=100)
    prep.add_argument('--skip-knowledge', action='store_true')
    run = commands.choices['run']
    run.add_argument('--output', type=Path, default=Path('outputs/latest'))
    run.add_argument('--config', type=Path, default=Path('config/ranking.json'))
    run.add_argument('--online', action='store_true', help='Allow requests for uncached VEP annotations')
    run.add_argument('--llm', action='store_true', help='Use optional LLM fact selection (GigaChat-2-Pro by default)')
    run.add_argument('--phenotype', type=Path, help='Override HPO profile for exactly one selected case')
    ui = commands.choices['serve']
    ui.add_argument('--port', type=int, default=8765)
    ui.add_argument('--output', type=Path, default=Path('outputs/ui'))
    ui.add_argument('--config', type=Path, default=Path('config/ranking.json'))
    ui.add_argument('--open', action='store_true', help='Open the local interface in the default browser')
    benchmark = commands.choices['benchmark']
    benchmark.set_defaults(split='development')
    benchmark.add_argument('--output', type=Path, default=Path('outputs/benchmark'))
    benchmark.add_argument('--config', type=Path, default=Path('config/ranking.json'))
    evaluate = commands.add_parser('evaluate')
    evaluate.add_argument('--dataset', type=Path, default=Path(DEFAULT_ARCHIVE))
    evaluate.add_argument('--reports', type=Path, default=Path('outputs/latest'))
    evaluate.add_argument('--output', type=Path, help='Defaults to <reports>/evaluation.json')
    check = commands.add_parser('llm-check', help='Check GigaChat access and strict JSON on an existing report')
    check_input = check.add_mutually_exclusive_group(required=True)
    check_input.add_argument('--report', type=Path)
    check_input.add_argument('--synthetic', action='store_true', help='Use an invented fixture without genomic/phenotype data')
    check.add_argument('--output', type=Path, default=Path('outputs/llm_check.json'))
    return root


def progress(message):
    print(message, file=sys.stderr, flush=True)


def main(argv=None):
    args = parser().parse_args(argv)
    dataset = cache = None
    try:
        if args.command == 'serve':
            from .web import serve
            return serve(args)
        if args.command == 'llm-check':
            from .gigachat import DEFAULT_MODEL, GigaChatReporter
            from .guards import check_presentation
            import os
            reporter = GigaChatReporter()
            model = os.environ.get('GENOMIC_LLM_MODEL', DEFAULT_MODEL)
            models = reporter.models()
            if model not in models:
                raise ValueError(f'The selected model {model} is not available to this API project')
            if args.synthetic:
                report = {'case_id': 'connection_check', 'status': 'uncertain', 'candidates': [
                    {'variant_id': 'demo_item_1', 'gene': 'DEMO', 'missing': [], 'evidence': [
                        {'id': 'ev_demo_1', 'criterion': 'consequence',
                         'text': 'Invented JSON connectivity fixture. No genomic or patient data.',
                         'source': 'synthetic_fixture', 'basis': 'test'}]}]}
            else:
                report = json.loads(args.report.read_text(encoding='utf-8'))
            selection = reporter.select(report)
            guard = check_presentation(selection, report)
            result = {'connection': 'passed', 'selected_model': model,
                      'input_mode': 'synthetic_fixture' if args.synthetic else 'report_top_candidates',
                      'guard': guard, 'model_metadata': reporter.last_metadata, 'selection': selection}
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
            print(json.dumps({'connection': 'passed', 'model': model,
                              'guard': guard['status'], 'output': str(args.output)}, ensure_ascii=False))
            return 0
        if args.command == 'evaluate':
            from .evaluation import evaluate
            result = evaluate(args.dataset, args.reports)
            output = args.output or args.reports / 'evaluation.json'
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0
        dataset = Dataset(args.dataset)
        cases = dataset.load_cases(args.split, args.cases)
        if not cases:
            raise ValueError('No cases selected')
        cache = AnnotationCache(args.cache)
        provider = VepProvider(cache)
        if args.command == 'prepare':
            if not args.skip_knowledge:
                fetch_knowledge(args.knowledge, progress=progress)
            annotations, source = provider.annotate(
                [v for case in cases for v in case.variants], batch_size=args.batch_size, progress=progress)
            complete = sum(a['status'] == 'annotated' for a in annotations.values())
            print(json.dumps({'annotated': complete, 'requested': len(annotations),
                              'source': source}, ensure_ascii=False, indent=2))
            return 0 if complete == len(annotations) else 2
        if args.command == 'benchmark':
            from .evaluation import compare_experts, write_comparison
            from .ranking import RankingConfig
            knowledge = Knowledge.load(args.knowledge)
            config = RankingConfig.load(args.config)
            annotations, source = provider.annotate([v for case in cases for v in case.variants], offline=True)
            result = compare_experts(cases, annotations, knowledge, config, args.dataset, progress)
            result['sources'] = {'vep': source, 'hpo': knowledge.manifest}
            write_comparison(result, args.output)
            print(json.dumps(result['metrics'], ensure_ascii=False, indent=2))
            return 0 if result['cases'] and not result['excluded_cases'] else 2
        from .pipeline import Pipeline, write_run
        from .ranking import RankingConfig
        if args.phenotype:
            if len(cases) != 1:
                raise ValueError('--phenotype requires exactly one --case')
            phenotype = json.loads(args.phenotype.read_text(encoding='utf-8-sig'))
            positive, excluded = {}, {}
            for feature in phenotype['phenotypicFeatures']:
                target = excluded if feature.get('excluded', False) else positive
                term = feature['type']['id']
                if term in positive or term in excluded:
                    raise ValueError('Repeated HPO term in override')
                target[term] = feature['type'].get('label', term)
            cases[0].positive_hpo, cases[0].excluded_hpo = positive, excluded
            cases[0].__post_init__()
            cases[0].metadata['phenotype_override'] = str(args.phenotype.name)
        kb = Knowledge.load(args.knowledge)
        config = RankingConfig.load(args.config)
        pipeline = Pipeline(provider, kb, config, llm=args.llm)
        reports = pipeline.run(cases, offline=not args.online, progress=progress)
        write_run(reports, args.output)
        for report in reports:
            print(f"{report['case_id']}: {report['status']}; prediction={report['prediction']}; "
                  f"coverage={report['annotation_coverage']:.1%}")
        return 0 if all(r['status'] != 'incomplete' for r in reports) else 2
    except (ValueError, OSError, RuntimeError, KeyError, json.JSONDecodeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    finally:
        if cache:
            cache.close()
        if dataset:
            dataset.close()
