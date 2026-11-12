import socket
import ssl
import unittest
import urllib.error
from unittest.mock import patch

from genomic_detective.network import HttpClient, NetworkError, connection_error_message


class NetworkTests(unittest.TestCase):
    def test_tls_error_explains_ca_without_leaking_query(self):
        reason = ssl.SSLCertVerificationError(1, 'certificate verification failed')
        reason.verify_code = 19
        error = urllib.error.URLError(reason)
        message = connection_error_message('https://api.giga.chat/v1?secret=do-not-log', error)
        self.assertIn('TLS certificate verification failed', message)
        self.assertIn('api.giga.chat', message)
        self.assertIn('19', message)
        self.assertNotIn('do-not-log', message)

    def test_dns_error_is_separate_from_tls(self):
        message = connection_error_message('https://api.giga.chat',
                                            urllib.error.URLError(socket.gaierror(11001, 'fixture')))
        self.assertIn('DNS lookup failed', message)

    def test_form_encoding_and_tls_context_are_passed(self):
        class Response:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

            def read(self):
                return b'{"ok": true}'

        def fake_open(request, **kwargs):
            self.assertEqual(request.data, b'scope=GIGACHAT_API_PERS')
            self.assertEqual(request.get_header('Content-type'), 'application/x-www-form-urlencoded')
            self.assertIn('context', kwargs)
            return Response()

        with patch('urllib.request.urlopen', side_effect=fake_open):
            self.assertEqual(HttpClient(retries=0).form_json('https://example.org',
                {'scope': 'GIGACHAT_API_PERS'}), {'ok': True})

    def test_failed_transport_raises_sanitized_error(self):
        error = urllib.error.URLError(socket.gaierror(11001, 'fixture-secret'))
        with patch('urllib.request.urlopen', side_effect=error):
            with self.assertRaises(NetworkError) as raised:
                HttpClient(retries=0).json('https://api.giga.chat/?secret=fixture-secret')
        self.assertIn('DNS lookup failed', str(raised.exception))
        self.assertNotIn('fixture-secret', str(raised.exception))


if __name__ == '__main__':
    unittest.main()
