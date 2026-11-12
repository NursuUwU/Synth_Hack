"""User-supplied cases: VCF + optional Phenopacket, online annotation of uncached variants."""
from __future__ import annotations

import gzip
import json
import re
import time
import urllib.parse
from pathlib import Path

from .models import Case, Variant, utc_now
from .network import HttpClient, NetworkError

MAX_VARIANTS = 5000
GENOTYPES = {'0/1', '1/0', '1/1', '0|1', '1|0', '1|1'}
OTHER_BUILDS = re.compile(r'(GRCh37|hg19|b37|hs37)', re.IGNORECASE)
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def decode(data: bytes) -> str:
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    return data.decode('utf-8-sig')


def parse_vcf(text: str) -> tuple[list[Variant], dict]:
    """Single-sample GRCh38 SNVs. Everything skipped is counted and reported, never silent."""
    variants, seen = [], set()
    skipped = {'not_snv': 0, 'multiallelic': 0, 'reference_or_missing_genotype': 0,
               'unsupported_chromosome': 0, 'duplicate': 0, 'filtered': 0}
    header_seen = False
    for number, line in enumerate(text.splitlines(), 1):
        if line.startswith('##'):
            if line.lower().startswith(('##reference', '##assembly', '##contig')) and OTHER_BUILDS.search(line):
                raise ValueError('VCF указывает на сборку GRCh37/hg19. Поддерживается только GRCh38.')
            continue
        if line.startswith('#'):
            columns = line.lstrip('#').split('\t')
            if len(columns) < 10:
                raise ValueError('В VCF нет колонки FORMAT и образца с генотипом.')
            if len(columns) > 10:
                raise ValueError('Поддерживается VCF с одним образцом (пробанд).')
            header_seen = True
            continue
        if not line.strip():
            continue
        if not header_seen:
            raise ValueError('Не найдена строка заголовка #CHROM.')
        fields = line.split('\t')
        if len(fields) < 10:
            raise ValueError(f'Строка {number}: ожидается 10 колонок VCF.')
        chrom, pos, _, ref, alt, _, status, _, fmt, sample = fields[:10]
        if status not in {'.', 'PASS', ''}:
            skipped['filtered'] += 1
            continue
        if ',' in alt:
            skipped['multiallelic'] += 1
            continue
        if len(ref) != 1 or len(alt) != 1 or ref.upper() not in 'ACGT' or alt.upper() not in 'ACGT':
            skipped['not_snv'] += 1
            continue
        keys = fmt.split(':')
        if 'GT' not in keys:
            raise ValueError(f'Строка {number}: нет поля GT.')
        values = sample.split(':')
        genotype = values[keys.index('GT')] if keys.index('GT') < len(values) else '.'
        if genotype not in GENOTYPES:
            skipped['reference_or_missing_genotype'] += 1
            continue
        try:
            variant = Variant(chrom, int(pos), ref.upper(), alt.upper(), genotype)
        except ValueError:
            skipped['unsupported_chromosome'] += 1
            continue
        if variant.id in seen:
            skipped['duplicate'] += 1
            continue
        seen.add(variant.id)
        variants.append(variant)
        if len(variants) > MAX_VARIANTS:
            raise ValueError(f'Больше {MAX_VARIANTS} вариантов. Отфильтруйте VCF до панели или региона '
                             '(онлайн-аннотация всего экзома в этом прототипе не поддерживается).')
    if not variants:
        raise ValueError('В VCF нет подходящих вариантов: нужны однонуклеотидные замены с генотипом 0/1 или 1/1.')
    return variants, {key: value for key, value in skipped.items() if value}


def parse_table(text: str) -> tuple[list[Variant], dict]:
    """Prepared variant table (TSV/CSV) with chrom, pos, ref, alt, genotype — the organizers' format."""
    import csv
    import io
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith('##')]
    if not lines:
        raise ValueError('Таблица вариантов пуста.')
    delimiter = '\t' if '\t' in lines[0] else ','
    rows = list(csv.DictReader(io.StringIO('\n'.join(lines)), delimiter=delimiter))
    header = {name.strip().lstrip('#').lower(): name for name in rows[0]} if rows else {}
    required = ['chrom', 'pos', 'ref', 'alt', 'genotype']
    if any(name not in header for name in required):
        raise ValueError('В таблице нужны колонки chrom, pos, ref, alt, genotype (как в *_variants.tsv организаторов).')
    vcf = ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND']
    for row in rows:
        values = [str(row[header[name]] or '').strip() for name in required]
        vcf.append(f'{values[0]}\t{values[1]}\t.\t{values[2]}\t{values[3]}\t.\t.\t.\tGT\t{values[4]}')
    # One validation path for both formats: the table is checked exactly like a VCF.
    return parse_vcf('\n'.join(vcf))


def parse_variants(text: str) -> tuple[list[Variant], dict]:
    first = next((line for line in text.splitlines() if line.strip()), '')
    return parse_vcf(text) if first.startswith('#') else parse_table(text)


def parse_phenopacket(text: str, labels: dict) -> tuple[dict, dict, str, str | None]:
    """Phenopacket v2 phenotypicFeatures; unknown HPO IDs are rejected with a list."""
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError('Файл фенотипа не является корректным JSON.') from error
    features = data.get('phenotypicFeatures') if isinstance(data, dict) else None
    if not isinstance(features, list):
        raise ValueError('В JSON нет списка phenotypicFeatures (формат Phenopacket).')
    positive, excluded, unknown = {}, {}, []
    for feature in features:
        term = (feature.get('type') or {}).get('id') if isinstance(feature, dict) else None
        if not isinstance(term, str) or not re.fullmatch(r'HP:\d{7}', term):
            raise ValueError('Каждый признак должен иметь type.id в формате HP:NNNNNNN.')
        if term not in labels:
            unknown.append(term)
            continue
        label = feature['type'].get('label') or labels[term]
        (excluded if feature.get('excluded') else positive)[term] = label
    if unknown:
        raise ValueError('Неизвестные HPO-термины: ' + ', '.join(sorted(unknown)[:10]))
    subject = data.get('subject') or {}
    return positive, excluded, subject.get('sex') or data.get('sex') or 'UNKNOWN_SEX', data.get('id')


def make_case(identifier: str, name: str, variants, positive, excluded, sex, skipped) -> Case:
    metadata = {'case_id': identifier, 'split': 'uploaded', 'uploaded': True, 'source_name': name,
                'genome_assembly': 'GRCh38', 'uploaded_at': utc_now(), 'skipped_records': skipped}
    return Case(identifier, variants, positive, excluded, metadata, sex)


def annotate_missing(provider, variants, batch_size=100) -> dict:
    """Fetch uncached variants into the existing cache namespace so offline runs reuse them."""
    annotations, source = provider.annotate(variants, offline=True)
    missing = [v for v in variants if annotations[v.id]['status'] != 'annotated']
    summary = {'cached': len(variants) - len(missing), 'fetched': 0, 'failed': 0}
    if not missing:
        return summary
    if source.get('namespace') == 'empty':
        raise ValueError('Локальный кэш аннотаций не инициализирован. Запустите prepare или Start.cmd ещё раз.')
    for start in range(0, len(missing), batch_size):
        records = provider._batch(missing[start:start + batch_size], source)
        summary['fetched'] += len(records)
    summary['failed'] = len(missing) - summary['fetched']
    if summary['fetched'] == 0:
        raise ValueError('Не удалось получить аннотации Ensembl VEP для новых вариантов. Проверьте интернет; '
                         'встроенные случаи работают офлайн.')
    return summary


def fetch_clinvar(annotations: dict, known: set, client: HttpClient | None = None, pause=0.4) -> list[dict]:
    """Same exact-allele ClinVar lookup as the bundled snapshot, for new VEP clinical hits only."""
    from .clinvar import parse_vcv
    requested = {key: value for key, value in annotations.items()
                 if value.get('clinical') and key not in known}
    rsids = sorted({c['record_id'] for a in requested.values() for c in a['clinical']
                    if (c.get('record_id') or '').startswith('rs') and c['record_id'][2:].isdigit()})
    client = client or HttpClient(timeout=60)
    found = {}
    for start in range(0, len(rsids), 20):
        query = urllib.parse.urlencode({'db': 'clinvar', 'term': ' OR '.join(rsids[start:start + 20]),
                                        'retmode': 'json', 'retmax': 500})
        search = client.json(EUTILS + 'esearch.fcgi?' + query)['esearchresult']
        time.sleep(pause)
        if not search['idlist']:
            continue
        url = EUTILS + 'efetch.fcgi?' + urllib.parse.urlencode(
            {'db': 'clinvar', 'rettype': 'vcv', 'is_variationid': 'true', 'id': ','.join(search['idlist'])})
        for record in parse_vcv(client.request(url, headers={'Accept': 'application/xml'})):
            if record['variant_id'] in requested:
                found.setdefault(record['variant_id'], {})[record['accession']] = record
        time.sleep(pause)
    return [record for records in found.values() for record in records.values()]


def save_supplement(path: Path, records: list[dict]):
    existing = json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {'records': {}}
    for record in records:
        bucket = existing['records'].setdefault(record['variant_id'], [])
        if record['accession'] not in {r['accession'] for r in bucket}:
            bucket.append(record)
    existing.update(schema_version=1, assembly='GRCh38', updated_at=utc_now(),
                    scope='Exact-allele ClinVar records for variants of user-uploaded cases')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(existing, ensure_ascii=False, indent=1) + '\n', encoding='utf-8')
