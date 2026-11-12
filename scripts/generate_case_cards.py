"""Generate case cards for the defense from the final v1_3 reports.

Numbers are read from outputs/v1_3 so cards cannot silently drift from the run.
Expected answers come from the organizer label file and are marked as such.
"""
from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CASES = {
    'case_01': {
        'title': 'FBN1 — marfanoid-progeroid-lipodystrophy spectrum',
        'positive_hpo': ['HP:0009064', 'HP:0005328', 'HP:0003758', 'HP:0001166', 'HP:0000426'],
        'negative_hpo': ['HP:0001083', 'HP:0004942', 'HP:0000545'],
        'negative_role': 'исключённая эктопия хрусталика и аневризма аорты спорят с классическим Марфаном; миопия неспецифична.',
        'competitor': 'SMAD3 / ORPHA:91387 (ближайший другой ген); FBN1 / OMIM:154700 классический Марфан (та же болезнь-спектр, ниже по баллу)',
        'competitor_text': '«за» редкий липодистрофический/прогероидный спектр: генерализованная липодистрофия, прогероидное лицо, снижение подкожного жира; «против» классического Марфана: эктопия хрусталика и аневризма аорты явно исключены.',
        'limitation': 'Синтетический составной генотип, нет семьи/фазы; condition-specific ClinVar с review status и функциональный эффект точного аллеля не подтверждены, поэтому статус `uncertain`, а не `candidate`.',
        'research': ['fbn1_2011', 'fbn1_2014'],
    },
    'case_03': {
        'title': 'KMT2D — Kabuki syndrome 1 / расширенный фенотип',
        'positive_hpo': ['HP:0000453', 'HP:0011381', 'HP:0008551', 'HP:0000430', 'HP:0000407'],
        'negative_hpo': ['HP:0001212', 'HP:0004209', 'HP:0000687'],
        'negative_role': 'исключённые prominent fingertip pads и клинодактилия V пальца нетипичны для классического Kabuki; это согласуется с расширенным фенотипом (PMID 32083401).',
        'competitor': 'PTPN11 / OMIM:163950 (ближайший другой ген); KMT2D / OMIM:620186 (вторая форма Kabuki, ниже по баллу)',
        'competitor_text': '«за» KMT2D: хоанальная атрезия, аплазия полукружных каналов, аномалии ушей, гипотиреоз; CHARGE-подобный спектр обсуждается клинически, но ближайший балл в этом прогоне — PTPN11 (RASopathy).',
        'limitation': 'Синтетический составной генотип; миссенс-эффект KMT2D не установлен, condition-specific ClinVar с review status отсутствует, поэтому `uncertain`.',
        'research': ['kmt2d_2010', 'kmt2d_2020'],
    },
    'case_05': {
        'title': 'LMNA — Hutchinson–Gilford progeria / прогероидный спектр',
        'positive_hpo': ['HP:0009125', 'HP:0003758', 'HP:0000905', 'HP:0001596', 'HP:0004322'],
        'negative_hpo': ['HP:0012632'],
        'negative_role': 'исключённое повышение внутриглазного давления не поддерживает альтернативный глазной спектр и не отменяет системные признаки.',
        'competitor': 'NRAS / OMIM:613224 (ближайший другой ген); TGFB2 / OMIM:614816 — третий по баллу',
        'competitor_text': 'LMNA поднимается точным литературным свидетельством для G608G: аномальный сплайсинг и укороченный продукт (PMID 12714972). VEP-класс остаётся `synonymous_variant`, поддержка последствия берётся как `max(0.1, 0.85)`, а не сумма.',
        'limitation': 'Нет RNA-эксперимента на входном образце, семьи и фазы; литературное de novo не является de novo входного случая. Балл — эвристическая поддержка, не вероятность диагноза и не ACMG-класс.',
        'research': ['lmna_2003'],
    },
    'case_06': {
        'title': 'TGFBR2 — Loeys–Dietz syndrome 2',
        'positive_hpo': ['HP:0005116', 'HP:0000193', 'HP:0001363', 'HP:0002616', 'HP:0100775'],
        'negative_hpo': ['HP:0001083', 'HP:0010648', 'HP:0001519'],
        'negative_role': 'исключённые эктопия хрусталика и прозрачность кожи спорят с классическим Марфаном и сосудистым EDS; пропорционально высокий рост также исключён.',
        'competitor': 'MRAS / OMIM:618499 (ближайший другой ген); клинически обсуждаются Марфан и сосудистый EDS',
        'competitor_text': '«за» LDS: извитость артерий, расщеплённый язычок, краниосиностоз, аневризма корня аорты, дуральная эктазия; «против» Марфана: эктопия хрусталика исключена; «против» vEDS: прозрачность кожи исключена.',
        'limitation': 'ClinVar-запись для G357W имеет 0 звёзд и не даёт положительной condition-specific поддержки; функциональные assay G357W извлечены не полностью. Поэтому `uncertain`.',
        'research': ['tgfbr2_2005'],
    },
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reports', type=Path, default=Path('outputs/v1_3'))
    parser.add_argument('--dataset', type=Path,
                        default=Path('Project_3_genomic_detecive_for_students.zip'))
    parser.add_argument('--output', type=Path, default=Path('demo/cards'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(args.dataset) as archive:
        prefix = next(name for name in archive.namelist() if name.endswith('development_labels.tsv'))
        text = archive.read(prefix).decode('utf-8-sig').splitlines()
    header = text[0].split('\t')
    labels = {row['case_id']: row for row in
              (dict(zip(header, line.split('\t'))) for line in text[1:] if line.strip())}

    cards = []
    for case_id, spec in CASES.items():
        report = json.loads((args.reports / f'{case_id}.json').read_text(encoding='utf-8'))
        top = report['candidates'][0]
        label = labels[case_id]
        positives = report['input']['positive_hpo']
        excluded = report['input']['excluded_hpo']
        key_positive = [f'`{term}` {positives.get(term, "")}' for term in spec['positive_hpo']]
        key_negative = [f'`{term}` {excluded.get(term, "")}' for term in spec['negative_hpo']]
        same_gene = top.get('other_disease_hypotheses') or []
        competitor = report['candidates'][1]
        phenotype = top['experts']['phenotype']
        molecular = top['experts']['molecular']
        inheritance = top['experts']['inheritance']
        lit = [ev for ev in top['evidence'] if ev['source'].startswith('LIT:')]
        func = [ev for ev in top['evidence'] if ev['criterion'] in {'functional_effect', 'literature_functional_effect'}]
        matches = '; '.join(f"{m['query']} → {m['shared_ancestor']}"
                            for m in phenotype['matches'] if m['query'] in spec['positive_hpo'][:8])
        if not matches:
            matches = '; '.join(f"{m['query']} → {m['shared_ancestor']}" for m in phenotype['matches'][:5])
        reasons = ', '.join(report.get('uncertainty_reasons') or []) or 'нет'
        card = [
            f'# Карточка {case_id}: {spec["title"]}',
            '',
            '## 1. Вход',
            f'- {report["variant_count"]} вариантов, хромосома {report["input"]["variants"][0]["chrom"]}.',
            '- Ключевые «+» HPO: ' + ', '.join(key_positive) + '.',
            '- Ключевые «−» HPO: ' + ', '.join(key_negative) + '.',
            f'- Дискриминирующая роль исключений: {spec["negative_role"]}',
            '',
            '## 2. Ожидаемый ответ',
            f'- {label["causal_gene"]} `{label["causal_variant"]}`, {label["disease_id"]}, PMID {label["source_url"].rstrip("/").split("/")[-1]}.',
            '- **Из открытой разметки организаторов, не наш вывод.** Полная независимая клиническая курация не проводилась.',
            '',
            '## 3. Результат системы',
            f'- Ранг {top["rank"]}, балл {top["score"]:.3f}, статус `{report["status"]}`.',
            f'- Ближайший конкурент: {competitor["gene"]} / {competitor["disease_id"]}, балл {competitor["score"]:.3f} (ранг {competitor["rank"]}).',
            f'- Причины неопределённости: {reasons}.',
            f'- Конкурирующая гипотеза: {spec["competitor"]}.',
            '',
            '## 4. Объяснение',
            f'- Фенотип: {phenotype["score"]:.3f}; совпадения: {matches}.',
            f'- Молекулярный канал: consequence {molecular["components"].get("consequence")}, '
            f'functional {molecular["components"].get("functional")}, clinical {molecular["components"].get("clinical")}, '
            f'frequency {molecular["components"].get("frequency")}.',
        ]
        if func:
            card.append('- Функциональное/литературное свидетельство: ' + ' | '.join(
                f'{ev["criterion"]} ({ev["source"]}): {ev["text"]}' for ev in func))
        if lit:
            card.append('- Литературный канал: ' + ', '.join(ev['source'] for ev in lit) + '.')
        card += [
            f'- Наследование: score {inheritance["score"]}, модели {inheritance.get("modes")}; '
            + ('модель наследования для выбранного disease_id не загружена — это ограничение, не отрицание AD. '
               if inheritance["score"] is None else '')
            + 'сегрегация и фаза недоступны.',
            f'- Конкуренция: {spec["competitor_text"]}',
            '',
            '## 5. Ограничение',
            spec['limitation'],
            '',
            '## 6. Источник',
            f'- Разметка: {label["source_url"]}.',
            '- Исследовательские карточки: `config/molecular_evidence.research.json` ('
            + ', '.join(spec['research']) + ').',
            '- Финальный прогон: `outputs/v1_3/' + case_id + '.json` (единая версия чисел).',
            '',
        ]
        cards.append((case_id, '\n'.join(card)))

    for case_id, card in cards:
        (args.output / f'{case_id}.md').write_text(card, encoding='utf-8')
    (args.output / 'all_cards.md').write_text(
        '# Карточки случаев (финальный прогон v1.2)\n\n' +
        '\n\n---\n\n'.join(card for _, card in cards), encoding='utf-8')
    print(json.dumps({'output': str(args.output), 'cards': [c for c, _ in cards]}, ensure_ascii=False))


if __name__ == '__main__':
    main()
