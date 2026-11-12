'use strict';
const $ = id => document.getElementById(id);
const API_VERSION = 3;  // must match web.API_VERSION
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
const CONSEQUENCES = {
  missense_variant: 'Замена аминокислоты', synonymous_variant: 'Синонимичная замена', stop_gained: 'Преждевременный стоп',
  stop_lost: 'Потеря стоп-кодона', start_lost: 'Потеря старт-кодона', frameshift_variant: 'Сдвиг рамки',
  splice_donor_variant: 'Нарушение сплайсинга', splice_acceptor_variant: 'Нарушение сплайсинга',
  splice_region_variant: 'Область сплайсинга', intron_variant: 'Интрон', '5_prime_UTR_variant': "5'-UTR",
  '3_prime_UTR_variant': "3'-UTR", upstream_gene_variant: 'Рядом с геном', downstream_gene_variant: 'Рядом с геном',
  intergenic_variant: 'Межгенный', non_coding_transcript_exon_variant: 'Некодирующий транскрипт'};
const SEVERITY = ['stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant', 'start_lost',
  'stop_lost', 'missense_variant', 'splice_region_variant', 'synonymous_variant', '5_prime_UTR_variant',
  '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'intron_variant', 'upstream_gene_variant',
  'downstream_gene_variant', 'intergenic_variant'];
const CLINVAR_CLASS = {pathogenic: 'Патогенный', likely_pathogenic: 'Вероятно патогенный', benign: 'Доброкачественный',
  likely_benign: 'Вероятно доброкачественный', uncertain_significance: 'VUS — значение неизвестно', conflicting: 'Противоречивые оценки'};
const MODES = {AD: 'доминантное', AR: 'рецессивное', XLD: 'X-сцепл. доминантное', XLR: 'X-сцепл. рецессивное'};
const stars = n => n > 0 ? '★'.repeat(n) : '';
function consequence(candidate) {
  const terms = candidate.experts.molecular.consequences || [];
  const main = [...terms].sort((a, b) => (SEVERITY.indexOf(a) + 99) % 99 - (SEVERITY.indexOf(b) + 99) % 99)[0];
  return {text: main ? CONSEQUENCES[main] || main.replaceAll('_', ' ') : 'Не определено', term: main || ''};
}
function frequency(candidate) {
  const facts = candidate.facts || {};
  if (facts.annotation_status !== 'annotated') return {text: 'нет данных', note: 'Аннотация не получена'};
  const af = facts.gnomad_af;
  if (af === null || af === undefined) return {text: 'нет в gnomAD', note: 'Аллель не встречается в gnomAD'};
  const percent = af * 100;
  const text = percent >= 0.1 ? `${percent.toFixed(2)} %` : percent >= 0.001 ? `${percent.toFixed(3)} %` : `${percent.toExponential(1)} %`;
  return {text: text.replace('.', ','), note: af > 0.001 ? 'Выше 0,1 % — слишком часто для доминантной болезни' : 'Редкий вариант'};
}
function clinvar(candidate) {
  const m = candidate.experts.molecular;
  const best = list => list.reduce((a, b) => (b.review_stars > (a?.review_stars ?? -1) ? b : a), null);
  const matched = best(m.clinical_assertions || []);
  if (m.clinical_context_status === 'matched' && matched)
    return {text: `${CLINVAR_CLASS[m.clinical_status] || matched.classification} ${stars(matched.review_stars)}`, note: 'для этой болезни', tone: m.clinical_status};
  if (m.clinical_context_status === 'gene_level')
    return {text: `${CLINVAR_CLASS[m.clinical_status] || ''} ${stars(best(m.gene_level_assertions || [])?.review_stars || 0)}`, note: 'болезнь не указана', tone: m.clinical_status};
  if (matched) return {text: matched.classification, note: 'без критериев (0 ★)', tone: ''};
  if ((m.evidence || []).some(e => e.criterion === 'clinical_context')) return {text: 'Есть запись', note: 'для другой болезни', tone: ''};
  return {text: 'Нет записи', note: '', tone: ''};
}
function symptoms(candidate) {
  const matches = candidate.experts.phenotype.matches || [];
  // Exact = same or more specific term in the disease profile; near = a close parent term (similarity >= 0.5).
  const exact = matches.filter(m => m.similarity >= 0.999);
  const near = matches.filter(m => m.similarity >= 0.5 && m.similarity < 0.999);
  const label = id => state.features.get(id)?.label || id;
  const names = [...exact.map(m => label(m.query)), ...near.map(m => `~${label(m.query)}`)];
  return {count: exact.length, near: near.length, total: state.positiveCount || 0, names};
}
function inheritance(candidate) {
  const expert = candidate.experts.inheritance; const value = expert.score;
  const modes = (expert.modes || []).map(m => MODES[m] || m).join(' / ');
  if (value === null || value === undefined) return {mark: '?', cls: 'inh-unknown', text: 'не указано', title: 'Тип наследования болезни не указан в HPO'};
  if (expert.moi_compatible === false || value === 0) return {mark: '✗', cls: 'inh-no', text: modes || 'не совместимо', title: expert.evidence?.[0]?.text || 'Не совместимо'};
  if (value >= 1) return {mark: '✓', cls: 'inh-yes', text: modes, title: 'Генотип совместим с типом наследования'};
  return {mark: '≈', cls: 'inh-partial', text: modes, title: expert.evidence?.[0]?.text || 'Совместимо частично'};
}
function clipCell(text, target, candidate, title) {
  // Truncated text: clicking it opens this candidate and scrolls smoothly to the matching section.
  const cell = el('td', 'clip-cell'); const span = el('button', 'clip', text);
  span.type = 'button'; span.title = title || text;
  span.addEventListener('click', event => { event.stopPropagation(); openCandidate(candidate, target); });
  cell.append(span); return cell;
}
function openCandidate(candidate, target) {
  state.selected = candidate.variant_id; renderRanking(); renderDetails(candidate);
  if (!target) return;
  const section = document.getElementById(target) || $('details');
  section.scrollIntoView({behavior: 'smooth', block: 'start'});
  section.classList.remove('flash'); void section.offsetWidth; section.classList.add('flash');
}
function renderRanking() {
  $('ranking').replaceChildren();
  for (const candidate of state.report.candidates.slice(0, state.visible)) {
    const selected = candidate.variant_id === state.selected;
    const row = el('tr', `${rowTone(candidate)} ${selected ? 'selected' : ''}`.trim());
    row.append(el('td', 'rank-number', String(candidate.rank).padStart(2, '0')));
    const name = el('td'); const button = el('button', 'gene-button', candidate.gene || 'Ген не определён');
    button.setAttribute('aria-label', `Доказательства: ${candidate.gene || candidate.variant_id}, место ${candidate.rank}`);
    button.setAttribute('aria-pressed', String(selected));
    name.append(button, el('span', 'variant-id', candidate.variant_id)); row.append(name);
    const hom = ['1/1', '1|1'].includes(candidate.genotype);
    const genotype = el('td', 'small-cell', hom ? 'гомо' : 'гетеро'); genotype.title = `${candidate.genotype}: ${hom ? 'вариант в обеих копиях' : 'вариант в одной копии'}`;
    row.append(genotype);
    const effect = consequence(candidate); const effectCell = el('td', 'small-cell', effect.text); effectCell.title = effect.term; row.append(effectCell);
    const freq = frequency(candidate); const freqCell = el('td', 'small-cell', freq.text); freqCell.title = freq.note; row.append(freqCell);
    const cv = clinvar(candidate); const cvCell = el('td', `small-cell cv-${cv.tone || 'none'}`);
    cvCell.append(el('span', 'cv-class', cv.text.trim())); if (cv.note) cvCell.append(el('span', 'cell-note', cv.note));
    row.append(cvCell);
    row.append(clipCell(candidate.disease_name || '—', 'sec-why', candidate,
      `${candidate.disease_name || ''} ${candidate.disease_id || ''}`.trim()));
    const sym = symptoms(candidate);
    row.append(clipCell(`${sym.count} из ${sym.total}${sym.names.length ? ': ' + sym.names.join(', ') : ''}`, 'sec-phenotype', candidate,
      sym.names.join('\n') || 'Совпадений нет'));
    const inh = inheritance(candidate); const inhCell = el('td', `small-cell ${inh.cls}`);
    inhCell.append(el('span', 'inh-mark', inh.mark), document.createTextNode(` ${inh.text}`)); inhCell.title = inh.title;
    row.append(inhCell, scoreCell(candidate.score, true));
    row.addEventListener('click', () => openCandidate(candidate));
    $('ranking').append(row);
  }
  $('shown-count').textContent = `${Math.min(state.visible, state.report.candidates.length)} из ${state.report.variant_count}`;
  $('more').hidden = state.visible >= state.report.candidates.length;
}
function whyHere(candidate) {
  // One or two arguments for the place and the main limitation, from the same evidence.
  const pros = [], cons = []; const e = candidate.experts; const sym = symptoms(candidate);
  const cv = clinvar(candidate); const freq = frequency(candidate);
  if ((e.phenotype.score ?? 0) >= 0.5) pros.push(`совпало ${sym.count} из ${sym.total} признаков${sym.near ? ` (+${sym.near} близких)` : ''}`);
  else cons.push(`слабое совпадение симптомов (${sym.count} из ${sym.total})`);
  if (['pathogenic', 'likely_pathogenic'].includes(cv.tone)) pros.push(`ClinVar: ${cv.text.trim().toLowerCase()} — ${cv.note}`);
  if ((e.molecular.curated_evidence || []).length || (e.molecular.literature_applied || []).length) pros.push('опубликован функциональный эффект этого аллеля');
  if (freq.text === 'нет в gnomAD' || (candidate.facts?.gnomad_af ?? 1) < 0.0001) pros.push(`вариант редкий (${freq.text})`);
  if (e.inheritance.moi_compatible === false) cons.push('не совместим с типом наследования');
  if (cv.text === 'Нет записи' || cv.note === 'для другой болезни') cons.push('нет записи ClinVar для этой болезни');
  else if (cv.tone === 'uncertain_significance') cons.push('значение варианта по ClinVar неизвестно');
  if (!cons.length) cons.push('семья и фаза вариантов не проверены');
  return {pros: pros.slice(0, 2), cons: cons.slice(0, 1)};
}
function sourceLinks(candidate) {
  const links = new Map(); const [chrom, pos, ref, alt] = candidate.variant_id.split(':');
  const add = (label, url) => { if (!links.has(label)) links.set(label, url); };
  for (const evidence of candidate.evidence || []) {
    if (evidence.source.startsWith('ClinVar:') && evidence.url) add(evidence.source.replace('ClinVar:', 'ClinVar '), evidence.url);
    if (evidence.source.startsWith('LIT:') || evidence.source.startsWith('Literature:')) add(evidence.source.replace(/^(LIT|Literature):/, ''), evidence.url || '');
  }
  add('gnomAD', `https://gnomad.broadinstitute.org/variant/${chrom.replace('chr', '')}-${pos}-${ref}-${alt}?dataset=gnomad_r4`);
  const id = candidate.disease_id || '';
  if (id.startsWith('OMIM:')) add(id, `https://omim.org/entry/${id.slice(5)}`);
  if (id.startsWith('ORPHA:')) add(id, `https://www.orpha.net/en/disease/detail/${id.slice(6)}`);
  if (id) add('HPO', `https://hpo.jax.org/browse/disease/${id}`);
  add('Ensembl VEP', `https://rest.ensembl.org/vep/human/region/${chrom.replace('chr', '')}:${pos}-${pos}/${alt}?content-type=application/json`);
  return [...links];
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
  const why = whyHere(candidate); const whySection = el('section', 'detail-section why'); whySection.id = 'sec-why';
  whySection.append(el('h4', '', 'Почему на этом месте'));
  const list = el('ul', 'why-list');
  why.pros.forEach(text => list.append(el('li', 'pro', text))); why.cons.forEach(text => list.append(el('li', 'con', text)));
  const e = candidate.experts; const inh = inheritance(candidate);
  whySection.append(list, el('p', 'helper sub-scores',
    `Итог ${score(candidate.score)} = фенотип ${score(e.phenotype.score)} × ${state.report.config.weights.phenotype}`
    + ` + молекулярные ${score(e.molecular.score)} × ${state.report.config.weights.molecular}`
    + ` + наследование ${score(e.inheritance.score)} × ${state.report.config.weights.inheritance} (${inh.mark})`));
  const sources = el('div', 'source-links'); sources.append(el('span', 'source-title', 'Источники:'));
  for (const [label, url] of sourceLinks(candidate)) {
    if (url) { const link = el('a', '', label); link.href = url; link.target = '_blank'; link.rel = 'noopener noreferrer'; sources.append(link); }
    else sources.append(el('span', '', label));
  }
  whySection.append(sources); container.append(whySection);
  for (const [key, label, id] of [['phenotype', 'Фенотип и противоречия', 'sec-phenotype'], ['molecular', 'Молекулярные данные', 'sec-molecular'], ['inheritance', 'Наследование', 'sec-inheritance']]) {
    const section = evidenceSection(label, candidate.experts[key], state.report.config.weights[key]); section.id = id; container.append(section);
  }
  if (candidate.missing.length) {
    const section = el('section', 'detail-section'); section.append(el('h4', '', 'Чего не хватает для проверки'));
    const missing = el('ul', 'missing-list');
    candidate.missing.forEach(text => missing.append(el('li', '', missingTranslations[text] || text)));
    section.append(missing); container.append(section);
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
  $('status-note').className = `status-note${{candidate: '', uncertain: ' soft'}[report.status] ?? ' neutral'}`;
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
  if (report.status === 'uncertain') legend.push(['tone-contested', uncertain ? `близкие гипотезы: разрыв с первым меньше ${gap}` : 'ведущая гипотеза, подтверждение неполное']);
  $('ranking-legend').replaceChildren(...legend.map(([cls, text]) => { const item = el('span', 'legend-item'); item.append(el('i', cls), text); return item; }));
  $('ranking-legend').append(el('span', 'legend-item', 'Наследование: ✓ совместимо · ≈ частично · ✗ нет · ? не указано'),
    el('span', 'legend-item', 'Обрезанный текст кликабелен — откроет доказательства'));
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
  state.positiveCount = payload.positive.length;
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
