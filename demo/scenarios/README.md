# Четыре статуса: фактические артефакты

Единая версия чисел — финальный прогон `outputs/v1_2` (правила `research-baseline-v1.2-literature-channel`).

| Статус | Условие | Артефакт | Факт |
|---|---|---|---|
| `candidate` | пороги пройдены, разрыв ≥ 0.025, нет ограничивающих причин | `outputs/v1_2/case_05.*`, `outputs/v1_2/case_08.*` | LMNA 0.703; COQ4 0.683 |
| `uncertain` | пороги пройдены, но разрыв < 0.025 **или** ограниченные condition-specific свидетельства | `outputs/v1_2/case_09.*` | GLI3 0.756, разрыв 0.136, причина `limited_condition_evidence` |
| `incomplete` | покрытие < 100 % либо нет информативного фенотипа/справочников | `outputs/incomplete_demo/*` | coverage 99.75 %, prediction None, exit 2, `predictions.partial.tsv` |
| `no_convincing_candidate` | обработка завершена, пороги не пройдены | в текущем наборе нет чистого примера | удаление причинного варианта дало фоновый `candidate` SMAD3; это ограничение, не доказательство отсутствия причины |

## Важное расхождение с ожиданиями ТЗ-2.2

- ТЗ ожидало `case_08 = uncertain`, `case_09 = candidate`; фактически `case_08 = candidate`, `case_09 = uncertain`.
- Причина: у case_09 сработала `limited_condition_evidence`, у case_08 разрыв 0.086 превысил порог. Сценарии исправлены по факту, не подогнаны.
- Оба случая имеют по 5 HPO и `phenotype_truncated = true`; разные статусы показывают, что сокращённый профиль сам по себе не предопределяет отказ.

## Формулировки

- `incomplete`: «данных для завершённого вывода недостаточно», **не** «причины нет».
- `no_convincing_candidate`: «в пределах обработанного VCF убедительного кандидата нет»; отсутствие кандидата не равно отсутствию болезни.
- `candidate`/`uncertain`: «внутренний исследовательский балл», не вероятность диагноза.

## Сценарные листы

- `demo/scenarios/case_08.md`
- `demo/scenarios/case_09.md`
- `demo/scenarios/case_01_no_causal.md`
- `demo/scenarios/incomplete_case_01.md`
- What-if: `outputs/sensitivity_case01.tsv`, `outputs/sensitivity_case01.md`, `outputs/whatif_case01_lipodystrophy/`.
