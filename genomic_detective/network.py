from __future__ import annotations

import json
import ssl
import socket
import time
import urllib.error
import urllib.parse
import urllib.request


class NetworkError(RuntimeError):
    pass


def connection_error_message(url: str, error: Exception) -> str:
    # Only the hostname, error class/code and fixed advice are exposed. No headers or bodies.
    host = urllib.parse.urlsplit(url).hostname or 'remote service'
    reason = error.reason if isinstance(error, urllib.error.URLError) else error
    if isinstance(reason, ssl.SSLCertVerificationError):
        return (f'TLS certificate verification failed for {host} '
                f'(verify code {reason.verify_code}). Configure the trusted application CA bundle; '
                'see docs/gigachat_setup.md.')
    if isinstance(reason, socket.gaierror):
        return f'DNS lookup failed for {host} (code {reason.errno})'
    if isinstance(reason, (TimeoutError, socket.timeout)):
        return f'Connection timed out for {host}'
    code = getattr(reason, 'winerror', None) or getattr(reason, 'errno', None)
    return f'Connection failed for {host} ({type(reason).__name__}, code {code})'


class HttpClient:
    def __init__(self, timeout: float = 45, retries: int = 2, ca_file: str | None = None):
        self.timeout = timeout
        self.retries = retries
        self.ca_file = ca_file

    def request(self, url: str, payload: dict | None = None, headers: dict | None = None,
                raw_data: bytes | None = None) -> bytes:
        if payload is not None and raw_data is not None:
            raise ValueError('JSON and form bodies cannot be combined')
        data = json.dumps(payload).encode() if payload is not None else raw_data
        context = None
        if self.ca_file:
            context = ssl.create_default_context()
            context.load_verify_locations(cafile=self.ca_file)
        request_headers = {'Accept': 'application/json', 'User-Agent': 'GenomicDetective/0.1'}
        if data is not None:
            request_headers['Content-Type'] = 'application/json'
        request_headers.update(headers or {})
        for attempt in range(self.retries + 1):
            request = urllib.request.Request(url, data=data, headers=request_headers)
            try:
                with urllib.request.urlopen(request, timeout=self.timeout, context=context) as response:
                    return response.read()
            except urllib.error.HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504} or attempt == self.retries:
                    # Do not log request headers, credentials or potentially sensitive bodies.
                    raise NetworkError(f'HTTP {error.code} from {urllib.parse.urlsplit(url).hostname}') from error
                try:
                    delay = min(30.0, float(error.headers.get('Retry-After', 2 ** attempt)))
                except ValueError:
                    delay = 2 ** attempt
            except (urllib.error.URLError, TimeoutError, OSError) as error:
                if attempt == self.retries:
                    raise NetworkError(connection_error_message(url, error)) from error
                delay = 2 ** attempt
            time.sleep(delay)
        raise NetworkError('Request failed')

    def json(self, url: str, payload: dict | None = None, headers: dict | None = None):
        try:
            return json.loads(self.request(url, payload, headers))
        except (json.JSONDecodeError, UnicodeError) as error:
            raise NetworkError('Remote response is not valid JSON') from error

    def form_json(self, url: str, fields: dict, headers: dict | None = None):
        request_headers = dict(headers or {})
        request_headers['Content-Type'] = 'application/x-www-form-urlencoded'
        body = urllib.parse.urlencode(fields).encode('utf-8')
        try:
            return json.loads(self.request(url, headers=request_headers, raw_data=body))
        except (json.JSONDecodeError, UnicodeError) as error:
            raise NetworkError('OAuth response is not valid JSON') from error
