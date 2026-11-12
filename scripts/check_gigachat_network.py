"""Check GigaChat transport without reading or transmitting any credentials."""
import json
import os
import ssl
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def main():
    bundle = Path(__file__).resolve().parent.parent / 'config' / 'certs' / 'gigachat-ca.pem'
    ca_file = os.environ.get('GIGACHAT_CA_FILE') or (str(bundle) if bundle.is_file() else None)
    context = ssl.create_default_context()
    if ca_file:
        context.load_verify_locations(cafile=ca_file)
    results = []
    for url in ['https://ngw.devices.sberbank.ru:9443/api/v2/oauth',
                'https://api.giga.chat/v1/models']:
        row = {'host': urllib.parse.urlsplit(url).netloc, 'extra_ca_configured': bool(ca_file)}
        try:
            with urllib.request.urlopen(url, context=context, timeout=20) as response:
                row.update(tls='ok', http=response.status)
        except urllib.error.HTTPError as error:
            row.update(tls='ok', http=error.code)
        except urllib.error.URLError as error:
            reason = error.reason
            row.update(tls='failed', reason_type=type(reason).__name__,
                       verify_code=getattr(reason, 'verify_code', None),
                       verify_message=getattr(reason, 'verify_message', None),
                       errno=getattr(reason, 'errno', None))
        except (OSError, TimeoutError) as error:
            row.update(tls='failed', reason_type=type(error).__name__, errno=getattr(error, 'errno', None))
        results.append(row)
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
