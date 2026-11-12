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

## Учебный набор

`Project_3_genomic_detecive_for_students.zip` — предоставленный для кейса учебный синтетический составной набор. Исходный архив не изменяется; его SHA-256 записан в манифесте снимка. Исходные VCF/TSV и открытые development-ответы остаются внутри архива организаторов. Модуль ранжирования не читает ответы; отдельная команда evaluate использует их для проверки.

## Python для Windows

При отсутствии установленного Python Windows-лаунчер распаковывает официальный embeddable ZIP **Python 3.14.10** из resources/ (Windows x64); для ARM64/32-bit получает подходящий ZIP с python.org. URL, архитектуры и опубликованные SHA-256 зафиксированы в `config/python-runtime.json`; архив проверяется перед запуском. Официальный Windows x64 ZIP распространяется без изменений вместе с оригинальной лицензией внутри. Рабочие файлы runtime сохраняются в `.runtime/`, за пределами Git. Источник: [страница выпуска Python](https://www.python.org/downloads/release/python-31410/). Установка в системные каталоги и изменение PATH не выполняются.
