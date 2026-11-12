"""Start a repository-only copy and exercise the real offline HTTP pipeline."""
from __future__ import annotations

import argparse
import json
import os
import queue
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--case', default='all', help='Case ID or all')
    args = parser.parse_args()
    root = args.root.resolve()
    # A clean smoke check must not inherit a developer's API credentials.
    environment = {key: value for key, value in os.environ.items()
                   if not any(word in key.upper() for word in ('GIGACHAT', 'OPENAI', 'GENOMIC_LLM'))}
    environment['PYTHONUTF8'] = '1'
    environment['PYTHONUNBUFFERED'] = '1'
    bootstrap = "import sys,runpy;sys.path.insert(0,sys.argv.pop(1));runpy.run_module('genomic_detective',run_name='__main__')"
    process = subprocess.Popen([sys.executable, '-I', '-c', bootstrap, str(root),
                                'serve', '--port', '0'], cwd=root, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding='utf-8')
    lines = queue.Queue()
    def consume():
        for line in process.stdout:
            lines.put(line.strip())
    reader = threading.Thread(target=consume, daemon=True)
    reader.start()
    try:
        deadline, origin = time.monotonic() + 180, None
        while time.monotonic() < deadline:
            try:
                line = lines.get(timeout=1)
            except queue.Empty:
                if process.poll() is not None:
                    raise RuntimeError('Server stopped before becoming ready')
                continue
            print(line, flush=True)
            if line.startswith('Local UI: '):
                origin = line.removeprefix('Local UI: ')
                break
        if not origin:
            raise RuntimeError('Server startup timed out')
        def request(path, payload=None, token=None):
            headers = {'Content-Type': 'application/json', 'Origin': origin}
            if token:
                headers['X-Workspace-Token'] = token
            data = json.dumps(payload).encode() if payload is not None else None
            with urllib.request.urlopen(urllib.request.Request(origin + path, data=data, headers=headers), timeout=60) as response:
                return response.read()
        initial = json.loads(request('/api/init'))
        assert len(initial['cases']) == 10
        assert initial['hpo_terms'] > 20_000
        selected = [c for c in initial['cases'] if args.case == 'all' or c['id'] == args.case]
        assert selected, 'No cases selected'
        for case in selected:
            payload = {'case_id': case['id'], 'positive': list(case['positive']), 'excluded': list(case['excluded'])}
            report = json.loads(request('/api/analyze', payload, initial['token']))
            assert report['annotation_coverage'] == 1
            assert report['status'] != 'incomplete'
            assert report['ranking_guard']['status'] == report['presentation_guard']['status'] == 'passed'
            assert len(report['candidates']) == 50
            if 'close_scores' in report.get('uncertainty_reasons', []):
                assert len(report['ambiguity']['variant_ids']) >= 2
            export = json.loads(request(report['downloads']['json']))
            assert export['case_id'] == case['id']
            assert len(export['candidates']) == case['variant_count']
            print(f"PASS {case['id']}: {case['variant_count']} variants, {report['status']}", flush=True)
        print(f'PASS: repository-only startup and {len(selected)} case(s), without API credentials.')
    finally:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=10)
        process.stdout.close()


if __name__ == '__main__':
    main()
