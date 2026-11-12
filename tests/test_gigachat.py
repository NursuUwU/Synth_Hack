import json
import os
import unittest
import uuid
from unittest.mock import patch

from genomic_detective.gigachat import API, OAUTH, DEFAULT_MODEL, GigaChatReporter
from genomic_detective.reporting import get_reporter, make_presentation, template_selection


class GigaChatTests(unittest.TestCase):
    def setUp(self):
        self.env = patch.dict(os.environ, {'GIGACHAT_CREDENTIALS': 'fixture-key'}, clear=True)
        self.env.start()
        self.addCleanup(self.env.stop)

    def report(self):
        return {'case_id': 'fixture', 'status': 'uncertain', 'candidates': [
            {'variant_id': 'chr1:100:A:T', 'gene': 'FIXTURE', 'missing': [], 'evidence': [
                {'id': 'ev_fixture', 'criterion': 'consequence', 'text': 'Fixture only',
                 'source': 'fixture', 'basis': 'prediction'}]}]}

    def test_default_provider_is_gigachat_pro(self):
        self.assertIsInstance(get_reporter(), GigaChatReporter)
        self.assertEqual(DEFAULT_MODEL, 'GigaChat-2-Pro')

    def test_oauth_token_is_reused_and_refreshed(self):
        requests = []
        clock = [1_700_000_000.0]

        class FakeHttp:
            def form_json(inner, url, fields, headers):
                self.assertEqual(url, OAUTH)
                self.assertEqual(fields, {'scope': 'GIGACHAT_API_PERS'})
                self.assertEqual(headers['Authorization'], 'Basic fixture-key')
                self.assertEqual(uuid.UUID(headers['RqUID']).version, 4)
                requests.append(url)
                return {'access_token': 'fixture-token', 'expires_at': (clock[0] + 1800) * 1000}

        reporter = GigaChatReporter(FakeHttp(), clock=lambda: clock[0])
        self.assertEqual(reporter.token(), 'fixture-token')
        self.assertEqual(reporter.token(), 'fixture-token')
        self.assertEqual(len(requests), 1)
        clock[0] += 1750
        reporter.token()
        self.assertEqual(len(requests), 2)

    def test_manual_access_token_does_not_request_oauth(self):
        with patch.dict(os.environ, {'GIGACHAT_ACCESS_TOKEN': 'fixture-token'}):
            self.assertEqual(GigaChatReporter(object()).token(), 'fixture-token')

    def test_gigachat_schema_request_and_guard(self):
        report = self.report()
        expected = template_selection(report)

        class FakeHttp:
            def form_json(inner, *args):
                return {'access_token': 'fixture-token', 'expires_at': 2800}

            def json(inner, url, payload=None, headers=None):
                self.assertEqual(headers['Authorization'], 'Bearer fixture-token')
                if url == API + '/models':
                    return {'data': [{'id': DEFAULT_MODEL}]}
                self.assertEqual(url, API + '/chat/completions')
                self.assertEqual(payload['model'], DEFAULT_MODEL)
                self.assertTrue(payload['response_format']['strict'])
                self.assertEqual(payload['response_format']['type'], 'json_schema')
                self.assertIn('properties', payload['response_format']['schema'])
                self.assertNotIn('fixture-key', json.dumps(payload))
                return {'model': DEFAULT_MODEL + ':fixture',
                    'choices': [{'finish_reason': 'stop', 'message': {'content': json.dumps(expected)}}],
                    'usage': {'prompt_tokens': 10, 'completion_tokens': 5, 'access_token': 'not-retained'}}

        reporter = GigaChatReporter(FakeHttp(), clock=lambda: 1000)
        self.assertIn(DEFAULT_MODEL, reporter.models())
        result = make_presentation(report, True, reporter)
        self.assertEqual(result['provider'], 'gigachat_evidence_selection')
        self.assertEqual(result['guard']['status'], 'passed')
        self.assertNotIn('access_token', reporter.last_metadata['usage'])

    def test_truncated_generation_falls_back(self):
        class FakeHttp:
            def json(inner, *args):
                return {'choices': [{'finish_reason': 'length', 'message': {'content': '{}'}}]}

        with patch.dict(os.environ, {'GIGACHAT_ACCESS_TOKEN': 'fixture-token'}):
            result = make_presentation(self.report(), True, GigaChatReporter(FakeHttp()))
        self.assertEqual(result['provider'], 'template')
        self.assertIn('ValueError', result['warning'])

    def test_missing_credentials_are_not_sent(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(ValueError):
                GigaChatReporter(object()).token()


if __name__ == '__main__':
    unittest.main()
