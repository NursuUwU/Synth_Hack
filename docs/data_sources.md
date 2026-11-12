# Источники и благодарности

## Human Phenotype Ontology

Приложение использует **Human Phenotype Ontology, HPO Consortium и Monarch Initiative**, версия онтологии **2026-09-01**. Версия отображается в интерфейсе после загрузки справочника. Исходные файлы в локальном снимке не изменены; ссылки и хеши первоначального получения находятся внутри `knowledge/manifest.json`, версия — в заголовках самих файлов.

- [Проект HPO](https://www.human-phenotype-ontology.org/).
- [Релиз 2026-09-01](https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01).
- [LICENSE.md этого релиза](https://github.com/obophenotype/human-phenotype-ontology/blob/v2026-09-01/LICENSE.md), [опубликованные условия HPO](https://human-phenotype-ontology.github.io/license.html).
- Gargano et al., *The Human Phenotype Ontology in 2024: phenotypes around the world*: [публикация](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767975/).

В снимок входят hp.obo, genes_to_phenotype.txt и phenotype.hpoa целиком. Для phenotype.hpoa дополнительно сохранён исходный заголовок с датой 2026-09-02. Идентификаторы OMIM, Orphanet и ссылки на публикации остаются атрибутами исходных аннотаций; отдельные полные базы OMIM/Orphanet не распространяются.

## Ensembl Variant Effect Predictor

В снимке сохранены исходные ответы **Ensembl VEP release 116**, координатная система GRCh38, сборка источника GRCh38.p14, только для вариантов учебного набора. Пространство имён, параметры запросов и время проверки источника записаны в `demo-data.manifest.json` и SQLite metadata.

- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).
- [Условия предоставления данных Ensembl](https://www.ensembl.org/info/about/legal/disclaimer.html).

Аннотации клинической значимости из VEP не представлены как новая самостоятельная экспертиза ClinVar. Их ограничения остаются в отчётах.

## ClinVar

В v1.2 добавлен отдельный нормализованный снимок **NCBI ClinVar от 12 ноября 2026 года** — `resources/clinvar-snapshot.json`. Запрошены все 114 уникальных вариантов с клиническими отметками в VEP-кэше учебного набора; фенотип и открытые ответы не использовались для отбора запросов. Это ограниченный снимок, не полный импорт базы.

Получение: официальные Entrez ESearch/EFetch, формат VCV XML. В снимке сохранены URL запросов, размеры и SHA-256 исходных XML. Берутся только простые точные SNV на основной хромосоме GRCh38 и germline RCV; записи гаплотипов, соматические классификации и альтернативные последовательности исключаются. Идентификаторы разных заболеваний в составном списке не объединяются в синонимы. Классификации без предоставленных критериев остаются видимыми, но не добавляют клиническую поддержку по правилам прототипа.

- [ClinVar: доступ к данным](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).
- [ClinVar: review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
- [ClinVar: RCV и VCV](https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/).

Снимок базы и литературный реестр имеют отдельные контрольные суммы в отчёте. Публикации перечислены в `config/molecular_knowledge.json` и [предметном разборе](molecular_evidence_review.md). Литературные сведения подготовлены после знакомства с открытыми случаями; это ограничивает оценку обобщающей способности модели.

## Происхождение учебного набора

`Project_3_genomic_detecive_for_students.zip` — предоставленный для кейса учебный синтетический составной набор. Исходный архив не изменяется; его SHA-256 записан в манифесте снимка. Исходные VCF/TSV и открытые development-ответы остаются внутри архива организаторов. Модуль ранжирования не читает ответы; отдельная команда evaluate использует их для проверки.

## Python для Windows

При отсутствии установленного Python Windows-лаунчер распаковывает официальный embeddable ZIP **Python 3.14.10** из resources/ (Windows x64); для ARM64/32-bit получает подходящий ZIP с python.org. URL, архитектуры и опубликованные SHA-256 зафиксированы в `config/python-runtime.json`; архив проверяется перед запуском. Официальный Windows x64 ZIP распространяется без изменений вместе с оригинальной лицензией внутри. Рабочие файлы runtime сохраняются в `.runtime/`, за пределами Git. Источник: [страница выпуска Python](https://www.python.org/downloads/release/python-31410/). Установка в системные каталоги и изменение PATH не выполняются.
