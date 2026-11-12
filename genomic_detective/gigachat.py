"""GigaChat-2-Pro with PERS OAuth, strict JSON and an in-memory access token."""
from __future__ import annotations

import json
import math
import os
import time
import uuid
from pathlib import Path

from .network import HttpClient
from .reporting import SELECTION_INSTRUCTIONS, build_selection_request


API = 'https://api.giga.chat/v1'
OAUTH = 'https://ngw.devices.sberbank.ru:9443/api/v2/oauth'
DEFAULT_MODEL = 'GigaChat-2-Pro'


class GigaChatReporter:
    provider_name = 'gigachat_evidence_selection'

    def __init__(self, client=None, clock=None):
        bundle = Path(__file__).resolve().parent.parent / 'config' / 'certs' / 'gigachat-ca.pem'
        ca_file = os.environ.get('GIGACHAT_CA_FILE') or (str(bundle) if bundle.is_file() else None)
        self.client = client or HttpClient(timeout=60, retries=0, ca_file=ca_file)
        self.clock = clock or time.time
        self._access_token = None
        self._expires_at = 0.0
        self.last_metadata = {}

    def token(self):
        supplied_token = os.environ.get('GIGACHAT_ACCESS_TOKEN', '').strip()
        if supplied_token:
            return supplied_token
        if self._access_token and self.clock() < self._expires_at - 60:
            return self._access_token
        credentials = os.environ.get('GIGACHAT_CREDENTIALS', '').strip()
        if credentials.startswith('Basic '):
            credentials = credentials[6:].strip()
        if not credentials or any(char.isspace() for char in credentials):
            raise ValueError('Set the local GIGACHAT_CREDENTIALS authorization key')
        scope = os.environ.get('GIGACHAT_SCOPE', 'GIGACHAT_API_PERS')
        if scope != 'GIGACHAT_API_PERS':
            raise ValueError('This free personal-project adapter requires GIGACHAT_API_PERS')
        result = self.client.form_json(OAUTH, {'scope': scope}, {
            'Authorization': f'Basic {credentials}', 'RqUID': str(uuid.uuid4())})
        if not isinstance(result, dict) or not isinstance(result.get('access_token'), str) or not result['access_token']:
            raise ValueError('OAuth did not return an access token')
        expires = float(result.get('expires_at', 0))
        if expires > 100_000_000_000:
            expires /= 1000  # API responses can use Unix milliseconds.
        if not math.isfinite(expires) or expires <= self.clock() + 60:
            raise ValueError('OAuth returned an expired or unusable token')
        self._access_token, self._expires_at = result['access_token'], expires
        return self._access_token

    def models(self):
        result = self.client.json(API + '/models', headers={'Authorization': f'Bearer {self.token()}'})
        if not isinstance(result, dict) or not isinstance(result.get('data'), list):
            raise ValueError('Unexpected models response')
        return [row['id'] for row in result['data'] if isinstance(row, dict) and isinstance(row.get('id'), str)]

    def select(self, report: dict):
        model = os.environ.get('GENOMIC_LLM_MODEL', DEFAULT_MODEL)
        schema, facts = build_selection_request(report)
        payload = {'model': model, 'stream': False, 'temperature': 0, 'max_tokens': 4096,
            'messages': [{'role': 'system', 'content': SELECTION_INSTRUCTIONS},
                         {'role': 'user', 'content': json.dumps(facts, ensure_ascii=False)}],
            'response_format': {'type': 'json_schema', 'schema': schema, 'strict': True}}
        result = self.client.json(API + '/chat/completions', payload,
                                   {'Authorization': f'Bearer {self.token()}'})
        if not isinstance(result, dict):
            raise ValueError('GigaChat response must be a JSON object')
        choices = result.get('choices', [])
        if len(choices) != 1 or choices[0].get('finish_reason') != 'stop':
            raise ValueError('GigaChat response was truncated, refused or incomplete')
        content = choices[0].get('message', {}).get('content')
        if not isinstance(content, str):
            raise ValueError('GigaChat returned no JSON content')
        selection = json.loads(content)
        # Only operational metadata is retained; credentials and access tokens stay in memory.
        usage = result.get('usage', {})
        safe_usage = {key: usage[key] for key in ('prompt_tokens', 'completion_tokens', 'total_tokens', 'precached_prompt_tokens')
                      if isinstance(usage, dict) and isinstance(usage.get(key), int)
                      and not isinstance(usage[key], bool) and usage[key] >= 0}
        self.last_metadata = {'provider': 'gigachat', 'requested_model': model,
                              'returned_model': result.get('model'), 'usage': safe_usage}
        return selection
