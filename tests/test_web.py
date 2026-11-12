"""HTTP contract and phenotype recalculation through the real orchestrator."""
import copy
import http.client
import json
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

from genomic_detective.models import Case, Variant
from genomic_detective.ranking import RankingConfig
from genomic_detective.web import LocalServer, Workspace
from test_workflow import FakeProvider, annotation, knowledge


class WebTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        a, b = Variant('1', 100, 'A', 'T', '0/1'), Variant('1', 200, 'G', 'C', '0/1')
        workspace = Workspace.__new__(Workspace)
        workspace.cases = {'fixture': Case('fixture', [a, b], {'HP:0000002': 'feature A'}, {})}
        workspace.knowledge = knowledge()
        workspace.config = RankingConfig.load(Path(__file__).resolve().parents[1] / 'config/ranking.json')
        workspace.cache_path = Path(self.temp.name) / 'cache.sqlite'
        workspace.output = Path(self.temp.name) / 'reports'
        workspace.token = 'test-token'
        workspace.lock = threading.Lock()
        workspace.exports = {}
        self.workspace = workspace
        self.provider = patch('genomic_detective.web.VepProvider', return_value=FakeProvider(
            {a.id: annotation(a, 'GENEA'), b.id: annotation(b, 'GENEB')}))
        self.provider.start()
        self.addCleanup(self.provider.stop)
        self.server = LocalServer(0, workspace)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.addCleanup(self.stop_server)

    def stop_server(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=3)

    def request(self, path, payload=None, headers=None):
        connection = http.client.HTTPConnection('127.0.0.1', self.server.server_port, timeout=10)
        request_headers = {'Host': self.server.authority}
        if payload is not None:
            request_headers.update({'Content-Type': 'application/json', 'X-Workspace-Token': 'test-token'})
        request_headers.update(headers or {})
        body = None if payload is None else json.dumps(payload)
        connection.request('GET' if payload is None else 'POST', path, body, request_headers)
        response = connection.getresponse()
        status, response_headers, data = response.status, dict(response.getheaders()), response.read()
        connection.close()
        return status, response_headers, data

    def analyze(self, positive, excluded=None):
        status, _, body = self.request('/api/analyze', {'case_id': 'fixture', 'positive': positive, 'excluded': excluded or []})
        self.assertEqual(status, 200, body)
        return json.loads(body)

    def test_recalculate_with_edited_profile_and_restore(self):
        original = copy.deepcopy(self.workspace.cases['fixture'])
        first = self.analyze(['HP:0000002'])
        second = self.analyze(['HP:0000003'])
        restored = self.analyze(['HP:0000002'])
        self.assertEqual(first['candidates'][0]['gene'], 'GENEA')
        self.assertEqual(second['candidates'][0]['gene'], 'GENEB')
        self.assertEqual(restored['candidates'][0]['gene'], 'GENEA')
        self.assertFalse(first['modified'])
        self.assertTrue(second['modified'])
        self.assertNotEqual(first['input_sha256'], second['input_sha256'])
        self.assertEqual(original, self.workspace.cases['fixture'])
        self.assertEqual(second['ranking_guard']['status'], 'passed')

    def test_without_informative_phenotype_withholds_prediction(self):
        report = self.analyze([])
        self.assertEqual(report['status'], 'incomplete')
        self.assertIsNone(report['prediction'])

    def test_score_gap_group_is_separate_from_evidence_uncertainty(self):
        a, b = self.workspace.cases['fixture'].variants
        with patch('genomic_detective.web.VepProvider', return_value=FakeProvider(
                {a.id: annotation(a, 'GENEA'), b.id: annotation(b, 'GENEA')})):
            report = self.analyze(['HP:0000002'])
        self.assertEqual(report['status'], 'uncertain')
        self.assertEqual(report['ambiguity']['variant_ids'], [a.id, b.id])
        self.assertEqual(report['ambiguity']['gap_threshold'], 0.025)
        clear = self.analyze(['HP:0000002'])
        self.assertEqual(clear['status'], 'uncertain')
        self.assertEqual(clear['uncertainty_reasons'], ['limited_condition_evidence'])
        self.assertEqual(clear['ambiguity']['variant_ids'], [])

    def test_exports_preserve_separate_runs(self):
        first = self.analyze(['HP:0000002'])
        second = self.analyze(['HP:0000003'])
        self.assertNotEqual(first['downloads']['json'], second['downloads']['json'])
        status, headers, data = self.request(first['downloads']['json'])
        self.assertEqual(status, 200)
        self.assertIn('attachment', headers['Content-Disposition'])
        self.assertEqual(json.loads(data)['candidates'][0]['gene'], 'GENEA')
        self.assertNotIn('annotations', first)  # Full raw report available only in explicit export.
        self.assertEqual(self.request('/download/../../.secrets/gigachat.dpapi')[0], 404)

    def test_cross_origin_and_wrong_token_rejected(self):
        payload = {'case_id': 'fixture', 'positive': [], 'excluded': []}
        for headers in [{'Origin': 'https://example.org'}, {'X-Workspace-Token': ''}, {'Host': 'attacker.example'}]:
            self.assertEqual(self.request('/api/analyze', payload, headers)[0], 403)
        self.assertEqual(self.request('/api/init', headers={'Host': 'attacker.example'})[0], 403)

    def test_profile_validation(self):
        for positive, excluded in [(['HP:9999999'], []), (['HP:0000002'] * 2, []), (['HP:0000002'], ['HP:0000002'])]:
            self.assertEqual(self.request('/api/analyze', {'case_id': 'fixture', 'positive': positive, 'excluded': excluded})[0], 400)
        self.assertEqual(self.request('/api/analyze', [1, 2])[0], 400)

    def test_busy_request_does_not_overwrite(self):
        self.workspace.lock.acquire()
        try:
            self.assertEqual(self.request('/api/analyze', {'case_id': 'fixture', 'positive': [], 'excluded': []})[0], 409)
            self.assertFalse(self.workspace.exports)
        finally:
            self.workspace.lock.release()

    def test_second_server_cannot_share_listening_port(self):
        with self.assertRaises(OSError):
            LocalServer(self.server.server_port, self.workspace)

    def test_bootstrap_search_and_assets(self):
        status, _, body = self.request('/api/init')
        self.assertEqual(status, 200)
        self.assertEqual(json.loads(body)['cases'][0]['id'], 'fixture')
        status, _, body = self.request('/api/hpo?q=feature%20B')
        self.assertEqual(json.loads(body), [{'id': 'HP:0000003', 'label': 'feature B'}])
        for asset in ['/', '/app.js', '/style.css']:
            status, headers, body = self.request(asset)
            self.assertEqual(status, 200)
            self.assertIn("frame-ancestors 'none'", headers['Content-Security-Policy'])
            self.assertGreater(len(body), 100)


if __name__ == '__main__':
    unittest.main()


class PortTests(unittest.TestCase):
    def test_default_port_moves_on_when_busy_but_explicit_port_fails(self):
        from genomic_detective import web
        try:
            busy = LocalServer(web.DEFAULT_PORT, None)
        except OSError:
            self.skipTest('Default port is held by another process')
        try:
            moved = web.bind(None)
            self.assertNotEqual(moved.server_port, web.DEFAULT_PORT)
            moved.server_close()
            with self.assertRaises(ValueError):
                web.bind(web.DEFAULT_PORT)
        finally:
            busy.server_close()

