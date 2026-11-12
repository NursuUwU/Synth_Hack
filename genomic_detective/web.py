"""Loopback-only, offline UI over the same ranking pipeline as the CLI."""
from __future__ import annotations

import copy
import json
import secrets
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

from .annotations import AnnotationCache, VepProvider
from .dataset import Dataset
from .knowledge import Knowledge
from .models import CONTEXT_MODES
from .pipeline import Pipeline, write_run
from .ranking import CANDIDATE_GAP, RankingConfig


ASSETS = Path(__file__).parent / 'static'


class BusyError(ValueError):
    pass


class Workspace:
    def __init__(self, dataset, cache, knowledge, config, output):
        source = Dataset(dataset)
        try:
            self.cases = {case.id: case for case in source.load_cases()}
        finally:
            source.close()
        self.cache_path = Path(cache)
        self.knowledge = Knowledge.load(Path(knowledge))
        self.config = RankingConfig.load(Path(config))
        from .literature import default_registry
        self.literature = default_registry()
        self.output = Path(output)
        self.token = secrets.token_urlsafe(32)
        self.lock = threading.Lock()
        self.exports = {}

    def case_info(self, case):
        return {'id': case.id, 'variant_count': len(case.variants), 'sex': case.sex,
                'split': case.metadata.get('split'),
                'truncated': bool(case.metadata.get('phenotype_truncated')),
                'positive': case.positive_hpo, 'excluded': case.excluded_hpo}

    def bootstrap(self):
        clinical_date = None
        if self.config.molecular_mode in CONTEXT_MODES:
            from .molecular import default_knowledge
            clinical_date = default_knowledge().metadata['clinical_snapshot_date'][:10]
        return {'token': self.token, 'cases': [self.case_info(c) for c in self.cases.values()],
                'clinical_snapshot_date': clinical_date,
                'hpo_version': self.knowledge.manifest.get('ontology_version', 'version unavailable'),
                'hpo_terms': len(self.knowledge.labels), 'diseases': len(self.knowledge.diseases),
                'weights': self.config.dictionary()['weights']}

    def search(self, query):
        query = query.strip().casefold()
        if len(query) < 2:
            return []
        terms = [(term, label) for term, label in self.knowledge.labels.items()
                 if query in term.casefold() or query in label.casefold()]
        terms.sort(key=lambda item: (item[0].casefold() != query,
                                   not item[1].casefold().startswith(query), item[1]))
        return [{'id': term, 'label': label} for term, label in terms[:20]]

    def profile(self, payload):
        if not isinstance(payload, dict) or payload.get('case_id') not in self.cases:
            raise ValueError('Выберите существующий случай.')
        case = copy.deepcopy(self.cases[payload['case_id']])
        for key, attribute in [('positive', 'positive_hpo'), ('excluded', 'excluded_hpo')]:
            terms = payload.get(key)
            if not isinstance(terms, list) or len(terms) > 250 or any(not isinstance(t, str) for t in terms):
                raise ValueError('HPO-профиль должен содержать список терминов (до 250).')
            if len(set(terms)) != len(terms):
                raise ValueError('В HPO-профиле есть повторяющиеся термины.')
            original = self.cases[case.id]
            labels = original.positive_hpo | original.excluded_hpo
            if any(t not in self.knowledge.labels and t not in labels for t in terms):
                raise ValueError('Неизвестный HPO-термин. Выберите его через поиск.')
            setattr(case, attribute, {t: self.knowledge.labels.get(t, labels.get(t, t)) for t in terms})
        case.__post_init__()
        original = self.cases[case.id]
        modified = (set(case.positive_hpo) != set(original.positive_hpo)
                    or set(case.excluded_hpo) != set(original.excluded_hpo))
        if modified:
            case.metadata['phenotype_override'] = 'interactive_ui'
        return case, modified

    def analyze(self, payload):
        case, modified = self.profile(payload)
        if not self.lock.acquire(blocking=False):
            raise BusyError('Другой анализ ещё выполняется. Повторите запуск через несколько секунд.')
        cache = None
        try:
            # SQLite connection belongs to the request thread; no shared connection.
            cache = AnnotationCache(self.cache_path)
            report = Pipeline(VepProvider(cache), self.knowledge, self.config,
                              literature=getattr(self, 'literature', None)).run([case], offline=True)[0]
            run_id = secrets.token_hex(12)
            directory = self.output / run_id
            write_run([report], directory)
            self.exports[run_id] = (directory, case.id)
            near = [c['variant_id'] for c in report['candidates']
                    if report['candidates'][0]['score'] - c['score'] < CANDIDATE_GAP
                    ] if report['status'] == 'uncertain' else []
            if len(near) < 2:
                near = []
            return {'case_id': case.id, 'status': report['status'], 'prediction': report['prediction'],
                    'uncertainty_reasons': report['uncertainty_reasons'],
                    'ambiguity': {'gap_threshold': CANDIDATE_GAP, 'variant_ids': near},
                    'created_at': report['created_at'], 'modified': modified,
                    'variant_count': report['variant_count'], 'annotation_count': report['annotation_count'],
                    'annotation_coverage': report['annotation_coverage'],
                    'candidates': report['candidates'][:50], 'warnings': report['warnings'],
                    'ranking_guard': report['ranking_guard'],
                    'presentation_guard': report['presentation']['guard'],
                    'config': report['config'], 'input_sha256': report['input_sha256'],
                    'downloads': {kind: f'/download/{run_id}/{kind}' for kind in ['html', 'json', 'tsv']}}
        finally:
            if cache:
                cache.close()
            self.lock.release()

    def download(self, run_id, kind):
        if run_id not in self.exports or kind not in {'html', 'json', 'tsv'}:
            raise KeyError('Unknown export')
        directory, case_id = self.exports[run_id]
        name = f'{case_id}_ranking.tsv' if kind == 'tsv' else f'{case_id}.{kind}'
        return directory / name


DEFAULT_PORT = 8765


class LocalServer(ThreadingHTTPServer):
    daemon_threads = True
    # SO_REUSEADDR permits multiple simultaneous listeners on Windows.
    allow_reuse_address = sys.platform != 'win32'
    allow_reuse_port = False

    def __init__(self, port, workspace):
        super().__init__(('127.0.0.1', port), Handler)
        self.workspace = workspace
        self.authority = f'127.0.0.1:{self.server_port}'
        self.origin = f'http://{self.authority}'


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        pass  # No query/profile logging.

    def send(self, status, body, content_type='application/json; charset=utf-8', filename=None):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; form-action 'self'")
        if filename:
            self.send_header('Content-Disposition', f'attachment; filename="{filename}"')
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def trusted(self, mutation=False):
        # Reject DNS rebinding and cross-origin writes; no CORS access.
        if self.headers.get('Host') != self.server.authority:
            return False
        origin = self.headers.get('Origin')
        if origin and origin != self.server.origin:
            return False
        return not mutation or secrets.compare_digest(
            self.headers.get('X-Workspace-Token', ''), self.server.workspace.token)

    def do_GET(self):
        if not self.trusted():
            return self.send(403, {'error': 'Недопустимый источник запроса.'})
        url = urlsplit(self.path)
        assets = {'/': ('index.html', 'text/html; charset=utf-8'),
                  '/app.js': ('app.js', 'text/javascript; charset=utf-8'),
                  '/style.css': ('style.css', 'text/css; charset=utf-8')}
        if url.path in assets:
            name, kind = assets[url.path]
            return self.send(200, (ASSETS / name).read_bytes(), kind)
        if url.path == '/api/init':
            return self.send(200, self.server.workspace.bootstrap())
        if url.path == '/api/hpo':
            query = parse_qs(url.query).get('q', [''])[0][:200]
            return self.send(200, self.server.workspace.search(query))
        if url.path.startswith('/download/'):
            parts = url.path.split('/')
            try:
                if len(parts) != 4:
                    raise KeyError('Unknown export')
                path = self.server.workspace.download(parts[2], parts[3])
                kinds = {'html': 'text/html; charset=utf-8', 'json': 'application/json; charset=utf-8',
                         'tsv': 'text/tab-separated-values; charset=utf-8'}
                return self.send(200, path.read_bytes(), kinds[parts[3]], path.name)
            except (KeyError, OSError):
                pass
        self.send(404, {'error': 'Страница не найдена.'})

    def do_POST(self):
        if not self.trusted(mutation=True):
            return self.send(403, {'error': 'Недопустимый источник запроса. Обновите страницу.'})
        if urlsplit(self.path).path != '/api/analyze':
            return self.send(404, {'error': 'Страница не найдена.'})
        try:
            size = int(self.headers.get('Content-Length', '0'))
            if size <= 0 or size > 65536 or self.headers.get('Content-Type') != 'application/json':
                return self.send(400, {'error': 'Ожидается JSON размером до 64 КБ.'})
            payload = json.loads(self.rfile.read(size).decode('utf-8'))
            return self.send(200, self.server.workspace.analyze(payload))
        except BusyError as error:
            self.send(409, {'error': str(error)})
        except (ValueError, TypeError, UnicodeError) as error:
            message = str(error) if isinstance(error, ValueError) and not isinstance(error, json.JSONDecodeError) else 'Некорректный запрос.'
            self.send(400, {'error': message})
        except Exception:
            self.send(500, {'error': 'Анализ не завершён. Проверьте локальный кэш и справочники; повторите запуск через CLI для диагностики.'})


def port_in_use(error):
    return getattr(error, 'winerror', None) == 10048 or error.errno in {48, 98, 10048}


def bind(port):
    # Bind first, before expensive HPO loading. An explicit busy port fails immediately;
    # the default moves on, since a stale earlier launch often still holds 8765.
    if port is not None:
        if not 0 <= port <= 65535:
            raise ValueError('Port must be between 0 and 65535')
        try:
            return LocalServer(port, None)
        except OSError as error:
            if port_in_use(error):
                raise ValueError(f'Port {port} is already in use. Open http://127.0.0.1:{port}/ '
                                 'or start with --port 8766.') from error
            raise
    for candidate in range(DEFAULT_PORT, DEFAULT_PORT + 11):
        try:
            server = LocalServer(candidate, None)
        except OSError as error:
            if port_in_use(error):
                continue
            raise
        if candidate != DEFAULT_PORT:
            print(f'Port {DEFAULT_PORT} is busy (probably an earlier launch is still running); '
                  f'using {candidate}.', flush=True)
        return server
    raise ValueError(f'Ports {DEFAULT_PORT}-{DEFAULT_PORT + 10} are all in use. Close earlier launches '
                     '(python.exe in Task Manager) or start with --port <number>.')


def serve(args):
    server = bind(args.port)
    with server:
        if (args.cache.resolve() == Path('data/cache/vep.sqlite3').resolve()
                and args.knowledge.resolve() == Path('data/knowledge').resolve()):
            from .demo_data import ensure_demo_data
            ensure_demo_data(Path.cwd(), progress=lambda text: print(text, flush=True))
        print('Loading local HPO reference data. Please wait...', flush=True)
        server.workspace = Workspace(args.dataset, args.cache, args.knowledge, args.config, args.output)
        print(f'Local UI: {server.origin}', flush=True)
        print('Offline annotations; no external API calls. Stop: Ctrl+C.', flush=True)
        if args.open:
            webbrowser.open(server.origin)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
    return 0
