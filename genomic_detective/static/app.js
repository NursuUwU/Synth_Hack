'use strict';
const $ = id => document.getElementById(id);
const API_VERSION = 2;  // must match web.API_VERSION
const state = {cases: [], token: '', current: null, features: new Map(), report: null,
  previous: new Map(), selected: null, visible: 10, busy: false, searchVersion: 0};
const statusTitles = {candidate: 'Кандидат найден', uncertain: 'Ведущая гипотеза',
  incomplete: 'Недостаточно данных', no_convincing_candidate: 'Убедительного кандидата нет'};
const statusNotes = {
  candidate: 'Ведущий вариант прошёл все пороги: фенотип, молекулярные данные и совместимость с типом наследования.',
  uncertain: 'Ведущий вариант выбран; ниже указано, какого подтверждения ему не хватает.',
  incomplete: 'Для вывода не хватает данных — причина в предупреждениях ниже. Рейтинг показан для просмотра.',
  no_convincing_candidate: 'Ни один вариант не прошёл пороги. В таблице — варианты, ближе всего подошедшие к ним.'};
const missingTranslations = {
  'Reviewed condition-specific ClinVar classification unavailable; no positive support assigned': 'Нет классификации ClinVar с предоставленными критериями для выбранного заболевания; положительная поддержка не назначена.',
  'Exact-allele ClinVar record unavailable in the bundled snapshot; not assumed benign': 'Точный аллель не найден в доступном снимке ClinVar; это не означает доброкачественность.',
  'NMD, exon relevance and disease-specific PVS1 requirements not established': 'Не установлены NMD, значимость экзона и условия применения PVS1 для этой болезни.',
  'Generic null consequence does not establish the selected disease mechanism': 'Общая метка потери функции не подтверждает механизм выбранного заболевания.',
  'gnomAD frequency unavailable; not assumed zero': 'Частота gnomAD неизвестна; она не приравнивается к нулю.',
  'Condition-specific ClinVar assertion and review status not retrieved': 'Не получены утверждение ClinVar для конкретного заболевания и статус экспертной проверки.',
  'Family segregation and phase unavailable': 'Нет данных о семейной сегрегации и фазе вариантов.',
  'Allele absent from VEP gnomAD data; absence is weak rarity evidence only': 'Аллель не найден в gnomAD; это лишь слабый довод редкости.',
  'ClinVar condition not specified; gene-level support only': 'В ClinVar заболевание не указано; поддержка только на уровне гена.',
  'GenCC gene-disease validity for the selected disease ID not found': 'В GenCC нет оценки для этого идентификатора болезни.',
  'Inheritance mode of this disease record unavailable': 'Для этой записи о болезни тип наследования не указан.',
  'Inheritance mode or genotype unavailable': 'Тип наследования болезни или генотип неизвестен.'};
const warningTranslations = {
  'ClinVar snapshot covers VEP clinical hits only; missing or unmatched condition evidence supplies no positive support, not a benign classification.': 'Снимок ClinVar охватывает варианты с клиническими аннотациями VEP. Отсутствие записи или совпадения по заболеванию не означает доброкачественность.',
  'Literature curation was developed after inspection of open cases; development improvement is not independent validation.': 'Литературные сведения собраны после знакомства с открытыми случаями. Улучшение на них не является независимой валидацией.',
  'Synthetic composite genotype; not a full patient genome or clinical benchmark.': 'Синтетический составной генотип: набор не представляет полный геном пациента или клиническую выборку.',
  'Heuristic weights/thresholds are uncalibrated; ranking score is not a probability.': 'Веса и пороги эвристические, без калибровки. Итоговый балл не является вероятностью.',
  'Only GRCh38 SNVs are supported in this baseline.': 'В этой версии поддерживаются только однонуклеотидные варианты GRCh38.',
  'No family segregation, phase confirmation or independent reference FASTA validation.': 'Не выполнены семейная сегрегация, подтверждение фазы и независимая проверка по референсной FASTA.',
  'VEP clinical annotations do not include independently fetched condition-specific ClinVar review status.': 'Аннотации VEP не содержат отдельно проверенного статуса ClinVar для конкретного заболевания.',
  'Absence of a convincing candidate does not prove absence of a genetic cause.': 'Отсутствие убедительного кандидата не доказывает отсутствие генетической причины.',
  'Phenotype was explicitly overridden for this run; not the original case profile.': 'Для этого расчёта использован изменённый HPO-профиль.',
  'The organizer marked the phenotype profile as truncated.': 'Организаторы пометили фенотипический профиль как сокращённый.',
  'User-uploaded case: variants absent from the local cache were annotated online (Ensembl VEP, ClinVar); input is not organizer data.': 'Загруженный случай: новые варианты аннотированы онлайн (Ensembl VEP, ClinVar).'};
function translateWarning(text) {
  if (warningTranslations[text]) return warningTranslations[text];
  let match = text.match(/^Inheritance mode set by the user: (\w+);/);
  if (match) return `Задан тип наследования ${match[1]}: гипотезы с другим типом не выдаются как ответ.`;
  match = text.match(/^Uploaded VCF records skipped: (.*)$/);
  if (match) return `Пропущены строки VCF: ${match[1].replace('not_snv', 'не SNV').replace('multiallelic', 'мультиаллельные')
    .replace('reference_or_missing_genotype', 'без альтернативного генотипа').replace('filtered', 'не PASS')
    .replace('unsupported_chromosome', 'неподдерживаемая хромосома').replace('duplicate', 'дубликаты')}.`;
  return text;
}

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined && text !== null) node.textContent = text;
  return node;
}
async function api(url, payload) {
  const options = payload === undefined ? {} : {method: 'POST',
    headers: {'Content-Type': 'application/json', 'X-Workspace-Token': state.token}, body: JSON.stringify(payload)};
  let response, data;
  try {
    response = await fetch(url, options);
    data = await response.json();
  } catch {
    throw new Error('Не удалось получить ответ локального сервера. Запустите .\\run.ps1 serve и дождитесь строки Local UI. Оставьте терминал открытым, затем обновите эту страницу.');
  }
  if (!response.ok) throw new Error(data.error || 'Не удалось выполнить запрос.');
  return data;
}
function showError(message) { $('error').textContent = message; $('error').hidden = !message; }
function changed() {
  const original = state.current;
  const equal = (mode, source) => {
    const selected = [...state.features].filter(([, f]) => f.mode === mode).map(([id]) => id);
    return selected.length === Object.keys(source).length && selected.every(id => id in source);
  };
  return !equal('positive', original.positive) || !equal('excluded', original.excluded);
}
function markEdited() {
  $('profile-state').textContent = changed() ? 'Профиль изменён · требуется пересчёт' : 'Исходный профиль';
  $('result').hidden = true;
  $('empty').hidden = false;
  const active = [...state.features.values()].filter(f => f.mode !== 'ignored').length;
  $('feature-count').textContent = active;
  $('run-label').textContent = state.report ? 'Пересчитать рейтинг' : 'Рассчитать рейтинг';
}
function renderFeatures() {
  $('features').replaceChildren();
  for (const [id, feature] of state.features) {
    const item = el('div', `feature is-${feature.mode}`);
    item.append(el('span', 'feature-name', feature.label));
    const bottom = el('div', 'feature-bottom');
    bottom.append(el('code', '', id));
    const group = el('div', 'segmented');
    group.setAttribute('role', 'group'); group.setAttribute('aria-label', feature.label);
    for (const [mode, label] of [['positive', 'Есть'], ['excluded', 'Нет'], ['ignored', '−']]) {
      const button = el('button', '', label);
      button.type = 'button'; button.disabled = state.busy;
      button.setAttribute('aria-pressed', String(feature.mode === mode));
      button.setAttribute('aria-label', `${feature.label}: ${mode === 'ignored' ? 'не учитывать' : mode === 'excluded' ? 'явно отсутствует' : 'есть'}`);
      button.addEventListener('click', () => { feature.mode = mode; renderFeatures(); markEdited(); });
      group.append(button);
    }
    bottom.append(group); item.append(bottom); $('features').append(item);
  }
  $('feature-count').textContent = [...state.features.values()].filter(f => f.mode !== 'ignored').length;
}
function resetProfile() {
  state.features.clear();
  for (const mode of ['positive', 'excluded']) {
    for (const [id, label] of Object.entries(state.current[mode])) state.features.set(id, {label, mode});
  }
  state.searchVersion++; $('hpo-search').value = ''; $('search-results').hidden = true;
  renderFeatures(); markEdited();
}
function selectCase(id) {
  state.current = state.cases.find(c => c.id === id);
  state.report = null; state.previous.clear(); state.selected = null;
  const kind = {development: 'Открытый случай', hidden: 'Скрытый случай', uploaded: 'Загруженный VCF'}[state.current.split] || '';
  $('case-meta').replaceChildren(el('span', '', `${state.current.variant_count} вариантов`), el('span', '', kind), el('span', '', 'GRCh38'));
  const skipped = Object.values(state.current.skipped || {}).reduce((a, b) => a + b, 0);
  if (skipped) $('case-meta').append(el('span', '', `пропущено строк: ${skipped}`));
  $('inheritance').value = '';
  $('truncated').hidden = !state.current.truncated;
  resetProfile(); showError('');
}
function setBusy(busy) {
  state.busy = busy;
  for (const id of ['case-select', 'reset', 'analyze', 'hpo-search', 'inheritance', 'upload-submit']) $(id).disabled = busy;
  renderFeatures();
  $('loading').hidden = !busy;
  $('results-panel')?.setAttribute('aria-busy', String(busy));
  if (busy) { $('result').hidden = true; $('empty').hidden = true; $('run-label').textContent = 'Идёт расчёт…'; }
}
function score(value) { return value === null || value === undefined ? '—' : value.toFixed(3); }
function scoreCell(value, total = false) {
  const cell = el('td', total ? 'total-score' : '', score(value));
  if (!total && value !== null && value !== undefined) {
    const track = el('span', 'score-track'); const bar = el('i');
    bar.style.width = `${Math.max(0, Math.min(1, value)) * 100}%`; track.append(bar); cell.append(track);
  }
  return cell;
}
function rowTone(candidate) {
  // Green: the answer passed every threshold. Yellow: contested leaders (uncertain status
  // or near-equal scores). Everything else stays neutral.
  // Near-equal groups are capped: with a sparse phenotype dozens of variants can tie.
  const report = state.report; const near = new Set(report.ambiguity.variant_ids);
  if (report.status === 'candidate' && candidate.rank === 1) return 'tone-best';
  if (report.status === 'uncertain' && (candidate.rank === 1 || (near.has(candidate.variant_id) && candidate.rank <= 5)))
    return 'tone-contested';
  return '';
}
function inheritanceCell(expert) {
  // Compatibility check, not a graded score: show its meaning instead of 1 / — .
  const value = expert.score;
  const [mark, title, cls] = value === null || value === undefined
    ? ['?', 'Тип наследования болезни не указан в HPO — не проверялось', 'inh-unknown']
    : expert.moi_compatible === false || value === 0 ? ['✗', 'Не совместимо с типом наследования болезни', 'inh-no']
    : value >= 1 ? ['✓', 'Генотип совместим с типом наследования болезни', 'inh-yes']
    : ['≈', `Совместимо частично (${score(value)}): например, вторая копия аллеля или фаза не установлены`, 'inh-partial'];
  const cell = el('td', `inh-cell ${cls}`, mark); cell.title = title; return cell;
}
function renderRanking() {
  $('ranking').replaceChildren();
  const showDelta = state.previous.size > 0;
  $('delta-head').hidden = !showDelta;
  for (const candidate of state.report.candidates.slice(0, state.visible)) {
    const selected = candidate.variant_id === state.selected;
    const row = el('tr', `${rowTone(candidate)} ${selected ? 'selected' : ''}`.trim());
    row.append(el('td', 'rank-number', String(candidate.rank).padStart(2, '0')));
    const name = el('td'); const button = el('button', 'gene-button', candidate.gene || 'Ген не определён');
    button.setAttribute('aria-label', `Доказательства: ${candidate.gene || candidate.variant_id}, место ${candidate.rank}`);
    button.setAttribute('aria-pressed', String(selected));
    name.append(button, el('span', 'variant-id', candidate.variant_id)); row.append(name, scoreCell(candidate.score, true));
    for (const expert of ['phenotype', 'molecular']) row.append(scoreCell(candidate.experts[expert].score));
    row.append(inheritanceCell(candidate.experts.inheritance));
    if (showDelta) {
      const old = state.previous.get(candidate.variant_id); const delta = old === undefined ? null : old - candidate.rank;
      const cell = el('td', delta > 0 ? 'delta-up' : delta < 0 ? 'delta-down' : '',
        delta === null ? 'новый' : delta === 0 ? '·' : `${delta > 0 ? '↑' : '↓'}${Math.abs(delta)}`);
      cell.title = 'Изменение места по сравнению с предыдущим расчётом этого случая'; row.append(cell);
    }
    row.addEventListener('click', () => { state.selected = candidate.variant_id; renderRanking(); renderDetails(candidate); });
    $('ranking').append(row);
  }
  $('shown-count').textContent = `${Math.min(state.visible, state.report.candidates.length)} из ${state.report.variant_count}`;
  $('more').hidden = state.visible >= state.report.candidates.length;
}
function evidenceSection(title, expert, weight) {
  const section = el('section', 'detail-section'); const heading = el('h4', '', title);
  heading.append(el('small', '', `${score(expert.score)} · вес ${weight}`)); section.append(heading);
  const list = el('ul', 'evidence-list');
  for (const evidence of expert.evidence || []) {
    const item = el('li', evidence.criterion === 'phenotype_conflict' ? 'conflict' : '');
    item.append(el('span', '', evidence.text));
    const source = el('span', 'evidence-source');
    const basis = {database: 'Справочник', prediction: 'Предсказание', compatibility: 'Совместимость',
      experimental_literature: 'Опубликованный эксперимент', literature_context: 'Литература', literature: 'Литература'}[evidence.basis] || evidence.basis;
    source.append(el('span', 'basis', basis));
    if (evidence.url?.startsWith('https://')) {
      const link = el('a', '', evidence.source); link.href = evidence.url;
      link.target = '_blank'; link.rel = 'noopener noreferrer'; source.append(link);
    } else source.append(document.createTextNode(evidence.source));
    item.append(source); list.append(item);
  }
  if (!list.children.length) list.append(el('li', '', 'Подтверждающие данные отсутствуют.'));
  section.append(list); return section;
}
function renderDetails(candidate) {
  const container = $('details'); container.replaceChildren();
  const heading = el('div', 'candidate-heading'); const name = el('div');
  name.append(el('h3', '', candidate.gene || 'Ген не определён'), el('p', 'helper', candidate.variant_id));
  heading.append(name, el('span', 'rank-badge', `${candidate.rank}-е место`));
  const disease = el('p', 'disease-name', candidate.disease_name || 'Заболевание не определено');
  disease.append(el('small', '', candidate.disease_id || ''));
  const facts = el('div', 'detail-facts');
  for (const value of [`Генотип ${candidate.genotype}`, candidate.assembly,
    candidate.experts.molecular.selected_transcript,
    ...(candidate.experts.molecular.consequences || [])]) if (value) facts.append(el('span', '', value));
  container.append(heading, disease, facts);
  for (const [key, label] of [['phenotype', 'Фенотип и противоречия'], ['molecular', 'Молекулярные данные'], ['inheritance', 'Наследование']])
    container.append(evidenceSection(label, candidate.experts[key], state.report.config.weights[key]));
  if (candidate.experts.molecular.scoring_policy === 'context_v2_max_not_sum') {
    container.append(el('p', 'helper', 'Молекулярная оценка использует максимум поддержки от прогноза, опубликованного функционального эффекта и классификации для конкретной болезни. Связанные свидетельства не суммируются; редкость учитывается отдельно. Балл эвристический.'));
  }
  if (candidate.missing.length) {
    const section = el('section', 'detail-section'); section.append(el('h4', '', 'Чего не хватает для проверки'));
    const list = el('ul', 'missing-list');
    candidate.missing.forEach(text => list.append(el('li', '', missingTranslations[text] || text)));
    section.append(list); container.append(section);
  }
}
function renderReport() {
  const report = state.report;
  $('empty').hidden = true; $('result').hidden = false;
  $('result-case').textContent = `${report.case_id.toUpperCase()} / РЕЗУЛЬТАТ АНАЛИЗА`;
  $('status-title').textContent = statusTitles[report.status] || report.status;
  $('result-profile').textContent = `${report.modified ? 'Изменённый' : 'Исходный'} HPO-профиль · ${new Date(report.created_at).toLocaleTimeString('ru-RU', {hour: '2-digit', minute: '2-digit'})}`;
  if (state.inheritanceUsed) $('result-profile').textContent += ` · наследование ${state.inheritanceUsed}`;
  const reasons = report.uncertainty_reasons || [];
  if (report.status === 'uncertain' && reasons.length === 1 && reasons[0] === 'close_scores')
    $('status-title').textContent = 'Конкурирующие гипотезы';
  const descriptions = {
    close_scores: 'У нескольких вариантов почти равные баллы — это равноправные гипотезы.',
    clinical_conflict: 'Клинические интерпретации в ClinVar расходятся.',
    functional_clinical_conflict: 'Функциональные данные и клиническая классификация расходятся.',
    limited_condition_evidence: 'Фенотип указывает на этот вариант, но записи ClinVar или функциональных данных для этой болезни нет.',
    gene_level_clinical_only: 'Вариант патогенный по ClinVar, но запись не привязана к конкретной болезни.'
  };
  $('status-note').textContent = report.status === 'uncertain' && reasons.length
    ? reasons.map(reason => descriptions[reason] || reason).join(' ') : statusNotes[report.status] || '';
  $('status-note').className = `status-note${report.status === 'candidate' ? '' : ' amber'}`;
  $('variant-count').textContent = report.variant_count;
  $('coverage').textContent = `${(report.annotation_coverage * 100).toFixed(0)}%`;
  const uncertain = report.status === 'uncertain' && report.ambiguity.variant_ids.length >= 2;
  const gap = report.ambiguity.gap_threshold.toLocaleString('ru-RU');
  const hasHypothesis = ['candidate', 'uncertain'].includes(report.status);
  $('priority-label').textContent = uncertain ? 'Близкие по баллу варианты' : hasHypothesis ? 'Первый в рейтинге' : 'Итоговый кандидат';
  $('leading-gene').textContent = uncertain ? report.ambiguity.variant_ids.length : hasHypothesis ? report.candidates[0]?.gene || '—' : '—';
  $('priority-note').textContent = uncertain ? `Разрыв с первым < ${gap}` : '';
  const legend = [];
  if (report.status === 'candidate') legend.push(['tone-best', 'ответ: прошёл все пороги']);
  if (report.status === 'uncertain') legend.push(['tone-contested', uncertain ? `спорные: разрыв с первым меньше ${gap} (показаны до 5)` : 'ведущая гипотеза с неполным подтверждением']);
  $('ranking-legend').replaceChildren(...legend.map(([cls, text]) => { const item = el('span', 'legend-item'); item.append(el('i', cls), text); return item; }));
  $('ranking-legend').append(el('span', 'legend-item', 'Наследование: ✓ совместимо · ≈ частично · ✗ нет · ? не указано'));
  for (const kind of ['html', 'tsv', 'json']) { $(`download-${kind}`).href = report.downloads[kind]; $(`download-${kind}`).setAttribute('download', ''); }
  $('guard-status').textContent = report.ranking_guard.status === 'passed' && report.presentation_guard.status === 'passed'
    ? 'Проверка целостности расчёта и ссылок на доказательства пройдена. Клиническую достоверность эта проверка не устанавливает.'
    : 'Проверка целостности не пройдена.';
  $('warnings').replaceChildren(...report.warnings.map(w => el('li', '', translateWarning(w))));
  state.visible = 10; state.selected = report.candidates[0]?.variant_id;
  renderRanking(); if (report.candidates[0]) renderDetails(report.candidates[0]);
  $('details').hidden = !report.candidates.length;
  $('profile-state').textContent = report.modified ? 'Изменённый профиль · рассчитано' : 'Исходный профиль · рассчитано';
  $('run-label').textContent = 'Пересчитать рейтинг';
}
async function analyze() {
  if (state.busy) return;
  showError(''); setBusy(true);
  const payload = {case_id: state.current.id, positive: [], excluded: [], inheritance: $('inheritance').value || null};
  state.inheritanceUsed = $('inheritance').value;
  for (const [id, f] of state.features) if (f.mode !== 'ignored') payload[f.mode].push(id);
  try {
    const report = await api('/api/analyze', payload);
    state.previous = new Map((state.report?.candidates || []).map(c => [c.variant_id, c.rank]));
    state.report = report; setBusy(false); renderReport();
  } catch (error) { setBusy(false); $('empty').hidden = false; showError(error.message); markEdited(); }
}
let searchTimer;
$('hpo-search').addEventListener('input', () => {
  clearTimeout(searchTimer); const version = ++state.searchVersion; const query = $('hpo-search').value.trim();
  $('search-results').hidden = true;
  if (query.length < 2) return;
  searchTimer = setTimeout(async () => {
    try {
      const terms = await api(`/api/hpo?q=${encodeURIComponent(query)}`);
      if (version !== state.searchVersion) return;
      const results = $('search-results'); results.replaceChildren(); results.hidden = false;
      for (const term of terms) {
        const button = el('button', '', term.label); button.append(el('small', '', term.id));
        button.addEventListener('click', () => {
          if (state.busy) return;
          state.features.set(term.id, {label: term.label, mode: 'positive'});
          $('hpo-search').value = ''; results.hidden = true; state.searchVersion++; renderFeatures(); markEdited();
        }); results.append(button);
      }
      if (!terms.length) results.append(el('p', '', 'Ничего не найдено. Попробуйте HPO ID или английское название.'));
    } catch (error) { if (version === state.searchVersion) showError(error.message); }
  }, 220);
});
$('case-select').addEventListener('change', event => selectCase(event.target.value));
$('inheritance').addEventListener('change', markEdited);
$('upload-toggle').addEventListener('click', () => { $('upload-form').hidden = !$('upload-form').hidden; });
function caseOption(c) {
  const option = el('option', '', c.split === 'uploaded' ? `${c.name} (загружен)` : c.id.replace('case_', 'Случай '));
  option.value = c.id; return option;
}
async function fileBase64(file) {
  const bytes = new Uint8Array(await file.arrayBuffer()); let binary = '';
  for (let i = 0; i < bytes.length; i += 0x8000) binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(binary);
}
$('upload-form').addEventListener('submit', async event => {
  event.preventDefault();
  if (state.busy) return;
  const vcf = $('upload-vcf').files[0]; const phenotype = $('upload-phenotype').files[0];
  if (!vcf) { showError('Выберите VCF-файл.'); return; }
  if (vcf.size > 12 * 1024 * 1024) { showError('VCF больше 12 МБ: отфильтруйте его до панели генов или региона.'); return; }
  showError(''); setBusy(true); $('upload-label').textContent = 'Аннотируем варианты…';
  $('loading').querySelector('h2').textContent = 'Загружаем и аннотируем VCF';
  try {
    const result = await api('/api/upload', {vcf: await fileBase64(vcf), name: $('upload-name').value,
      phenotype: phenotype ? await phenotype.text() : null});
    state.cases.push(result.case); $('case-select').append(caseOption(result.case));
    $('case-select').value = result.case.id; selectCase(result.case.id);
    $('upload-form').reset(); $('upload-form').hidden = true;
    const a = result.annotation;
    $('profile-state').textContent = `Загружено: ${result.case.variant_count} вариантов (из кэша ${a.cached}, онлайн ${a.fetched}`
      + `${a.failed ? `, не аннотировано ${a.failed}` : ''}${a.clinvar_lookup === 'failed' ? ', ClinVar недоступен' : ''}).`;
  } catch (error) { showError(error.message); }
  finally {
    setBusy(false); $('upload-label').textContent = 'Загрузить и аннотировать';
    $('loading').querySelector('h2').textContent = 'Сопоставляем доказательства';
    if (state.report) $('result').hidden = false; else $('empty').hidden = false;
  }
});
$('reset').addEventListener('click', resetProfile);
$('analyze').addEventListener('click', analyze);
$('more').addEventListener('click', () => { state.visible += 10; renderRanking(); });
(async () => {
  try {
    const initial = await api('/api/init'); state.cases = initial.cases; state.token = initial.token;
    if (initial.api_version !== API_VERSION)
      showError('Запущен сервер старой версии: страница обновилась, а программа — нет. Закройте окно Start.cmd '
        + '(и лишние python.exe в Диспетчере задач) и запустите Start.cmd снова.');
    $('hpo-credit').textContent = `Human Phenotype Ontology · ${initial.hpo_version || 'версия не передана'}`;
    $('clinvar-credit').hidden = !initial.clinical_snapshot_date;
    if (initial.clinical_snapshot_date) $('clinvar-credit').textContent = `ClinVar · снимок ${initial.clinical_snapshot_date}`;
    if (!state.cases.length) throw new Error('В наборе нет случаев.');
    $('case-select').replaceChildren(...state.cases.map(caseOption));
    selectCase(state.cases[0].id); setBusy(false);
  } catch (error) { showError(`Не удалось открыть рабочую область: ${error.message}`); }
})();
