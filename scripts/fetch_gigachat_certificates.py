"""Download official public CA files over verified HTTPS for app-scoped GigaChat TLS."""
import hashlib
import json
import ssl
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def main():
    directory = Path(__file__).resolve().parent.parent / 'config' / 'certs'
    directory.mkdir(parents=True, exist_ok=True)
    names = ['russian_trusted_root_ca_pem.crt', 'russian_trusted_sub_ca_pem.crt']
    records, certificates = [], []
    for name in names:
        url = 'https://gu-st.ru/content/lending/' + name
        with urllib.request.urlopen(url, timeout=30) as response:
            data = response.read()
        text = data.decode('ascii').strip() + '\n'
        if '-----BEGIN CERTIFICATE-----' not in text or '-----END CERTIFICATE-----' not in text:
            raise ValueError('Downloaded content is not a PEM certificate')
        # Validate PEM syntax before using anything as a CA bundle.
        context = ssl.create_default_context()
        context.load_verify_locations(cadata=text)
        certificates.append(text)
        records.append({'url': url, 'sha256': hashlib.sha256(data).hexdigest()})
        print('Downloaded and parsed:', name)
    bundle = directory / 'gigachat-ca.pem'
    bundle.write_text('\n'.join(certificates), encoding='ascii')
    (directory / 'sources.json').write_text(json.dumps({
        'downloaded_at': datetime.now(timezone.utc).isoformat(),
        'documentation': 'https://developers.sber.ru/docs/ru/gigachat/certificates',
        'sources': records, 'bundle_sha256': hashlib.sha256(bundle.read_bytes()).hexdigest(),
    }, indent=2) + '\n', encoding='utf-8')
    print('Application CA bundle created:', bundle.name)


if __name__ == '__main__':
    main()
