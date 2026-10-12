// Suite workspace: inputs → simulation → explained outputs, equations, checks, recommendations,
// sensitivity / uncertainty, specification knowledge and live reference data.
import { h, icon, toast, download, tabs, table, toCSV, parseCSV, fileBytes, debounce } from '../core/ui.js';
import { fmt } from '../core/num.js';
import { lineChart, heatmap, barChart, histChart } from '../core/plot.js';
import { viewer } from '../core/view3d.js';
import { loadSuite, loadSpec, resolveInputs, runSuite, defaultInputs, finishRun, validateInputs, benchStatus } from '../core/engine.js';
import { runSuiteBenchmarks, loadBenchmarks, summarise, KIND as BENCH_KIND } from '../core/benchmarks.js';
import { kv } from '../core/store.js';
import { project, setInputs, saveProject, onProject } from '../core/store.js';
import { visible, oat, morris, sobol, monteCarlo, numericInputs } from '../core/sensitivity.js';
import * as live from '../core/live.js';
import { buildReport } from '../core/report.js';

const EQ_KINDS = { chemical: 'Chemical equations', ionic: 'Ionic equations', 'half-cell': 'Half-cell reactions', governing: 'Governing equations', constitutive: 'Constitutive relations', formula: 'Formulae' };
const PRI = { high: 'Act first', medium: 'Recommended', low: 'Consider' };

export async function suitePage(root, id) {
  const spec = await loadSpec(), s = spec.suites[id - 1];
  if (!s) { root.append(h('p', 'Unknown suite.')); return; }
  document.title = `Suite ${id} · ${s.name} — Designer Materials`;
  let suite; try { suite = await loadSuite(id); } catch (e) { root.append(h('div.card.error', h('h2', 'This suite could not be loaded'), h('p', e.message))); return; }
  let { values, origin } = resolveInputs(suite), last = null;
  const state = { running: false };

  root.append(h('header.pagehead',
    h('div.crumbs', h('a', { href: '#/suites' }, 'Suites'), ' / ', `Suite ${id}`),
    h('div.headrow', h('span.badge.num', String(id).padStart(2, '0')), h('div', h('h1', s.name), h('p.lead', suite.tagline), h('p.muted', s.resp + '.')),
      h('div.pager', id > 1 ? h('a.btn.ghost.sm', { href: `#/suite/${id - 1}`, title: spec.suites[id - 2].name }, icon('back', 16), ` ${id - 1}`) : null, id < 29 ? h('a.btn.ghost.sm', { href: `#/suite/${id + 1}`, title: spec.suites[id].name }, `${id + 1} `, icon('fwd', 16)) : null))));

  // ---------- Simulate tab ----------
  const formBox = h('form.inputs', { onsubmit: e => { e.preventDefault(); run(); } }), resBox = h('div.results');
  function field(i) {
    const idn = `in-${id}-${i.key}`, type = i.type || 'number'; let ctl;
    const set = v => { values[i.key] = v; origin[i.key] = 'user'; if (type === 'select' || type === 'bool') renderForm(); else chip.replaceChildren(); };
    const chip = h('span.chips');
    if (origin[i.key] === 'bus' && i.bus && project.bus[i.bus]) chip.append(h('a.chip.link', { href: `#/suite/${project.bus[i.bus].from}`, title: `Value published by Suite ${project.bus[i.bus].from} (${i.bus}). Edit to override.` }, icon('link', 12), ` from Suite ${project.bus[i.bus].from}`));
    if (type === 'select') ctl = h('select', { id: idn, onchange: e => set(e.target.value) }, i.options.map(o => { const v = typeof o === 'object' ? o.value : o; return h('option', { value: v, selected: v === values[i.key] }, typeof o === 'object' ? o.label : o); }));
    else if (type === 'bool') ctl = h('label.switch', h('input', { type: 'checkbox', id: idn, checked: !!values[i.key], onchange: e => set(e.target.checked) }), h('span'));
    else if (type === 'text') ctl = h('input', { id: idn, type: 'text', value: values[i.key] ?? '', oninput: e => set(e.target.value), spellcheck: false });
    else if (type === 'table') {
      const ta = h('textarea', { id: idn, rows: 6, spellcheck: false, value: values[i.key] ?? '', oninput: e => set(e.target.value) });
      ctl = h('div.tablein', ta, h('div.row', h('label.btn.sm.ghost', icon('upload', 14), ' Upload data file', h('input', { type: 'file', hidden: true, accept: '.csv,.tsv,.txt,.dat,.xy,.json,.xlsx', onchange: async e => { const f = e.target.files[0]; if (!f) return; try { const { readTabular } = await import('../core/dataio.js'); const t = await readTabular(f.name, await fileBytes(f)); ta.value = toCSV(t.columns, t.rows); set(ta.value); toast(`Loaded ${t.rows.length} rows from ${f.name}`, 'ok'); } catch (err) { toast(err.message, 'bad', 6000); } } })), h('span.muted', i.columns ? `Columns: ${i.columns.join(', ')}` : '')));
    } else {
      const inp = h('input', { id: idn, type: 'number', value: values[i.key], min: i.min, max: i.max, step: i.step || 'any', inputmode: 'decimal', oninput: e => { const v = e.target.value === '' ? NaN : Number(e.target.value); set(v); inp.classList.toggle('bad', !(v >= i.min && v <= i.max)); } });
      ctl = h('div.numrow', inp, i.unit ? h('span.unit', i.unit) : null);
      if (i.live && live.FEEDS[i.live]) chip.append(h('button.chip.live', { type: 'button', title: `Fill from ${live.FEEDS[i.live].source}`, onclick: async () => { try { const r = await live.feed(i.live, { force: live.online() }); const v = Math.min(i.max, Math.max(i.min, r.data.value)); inp.value = v; set(v); toast(`${live.FEEDS[i.live].label}: ${fmt(r.data.value)} ${live.FEEDS[i.live].unit} (${r.data.period}${r.stale ? ', cached' : ''})`, 'ok'); } catch (e) { toast('Live value unavailable: ' + e.message, 'bad'); } } }, icon('wifi', 12), ' live'));
    }
    return h('div.field', { class: type === 'table' ? 'wide' : '' }, h('label', { for: idn }, i.label, chip), ctl, h('p.help', i.help, type === 'number' ? h('span.range', ` Range ${fmt(i.min)} – ${fmt(i.max)}${i.unit ? ' ' + i.unit : ''}.`) : null));
  }
  function renderForm() {
    const groups = new Map(); for (const i of suite.inputs) { if (!visible(i, values)) continue; const g = i.group || 'Inputs'; if (!groups.has(g)) groups.set(g, []); groups.get(g).push(i); }
    const preset = h('select', { 'aria-label': 'Load a worked case', onchange: e => { const p = suite.presets[+e.target.value]; if (!p) return; values = { ...defaultInputs(suite), ...p.values }; for (const k in values) origin[k] = k in p.values ? 'user' : 'default'; renderForm(); run(); toast(`Loaded case: ${p.name}`, 'ok'); } }, h('option', { value: '' }, 'Load a worked case…'), suite.presets.map((p, k) => h('option', { value: k, title: p.note || '' }, p.name)));
    formBox.replaceChildren(
      h('div.formtop', preset, h('button.btn.ghost.sm', { type: 'button', title: 'Reset to defaults and re-link values published by other suites', onclick: () => { delete project.inputs[id]; saveProject(); ({ values, origin } = resolveInputs(suite)); renderForm(); } }, icon('refresh', 14), ' Reset')),
      ...[...groups].map(([g, list], k) => h('details.group', { open: k < 3 }, h('summary', g, h('span.count', list.length)), h('div.fields', list.map(field)))),
      h('div.runbar', h('button.btn.primary', { type: 'submit' }, icon('play', 16), ' Run simulation'), h('span.muted', 'Runs on this device — works offline.')),
      hifiPanel());
  }
  // Full-resolution background run (Web Worker): lifts the interactive grid / time-step caps.
  const hv = {}; let worker = null;
  function hifiPanel() {
    const H = suite.hifi; if (!H || (H.modes && !H.modes.includes(values.mode))) return null;
    for (const i of H.inputs) if (!(i.key in hv)) hv[i.key] = i.value;
    const plan = h('p.hint'), prog = h('progress', { max: 1, value: 0, hidden: true }), stat = h('p.hint'), start = h('button.btn', { type: 'button' }, icon('play', 15), ' Start full-resolution run'), stop = h('button.btn.ghost.sm', { type: 'button', hidden: true }, 'Stop and keep the result so far');
    const fmtT = s => s < 90 ? `${Math.round(s)} s` : s < 5400 ? `${Math.round(s / 60)} min` : `${(s / 3600).toFixed(1)} h`;
    const showPlan = async () => { try { const m = await import(`../suites/s${String(id).padStart(2, '0')}.js`); if (m.hifiPlan) { const p = m.hifiPlan({ ...values, ...hv }); plan.textContent = `${p.cells.toLocaleString()} cells × ${p.steps.toLocaleString()} steps. Estimated ${fmtT(p.estimatedSeconds)} on a typical laptop (longer on a phone), about ${Math.round(p.memoryMB)} MB of memory. The estimate is refined once the run starts.`; } } catch { plan.textContent = ''; } };
    const busy = b => { start.hidden = b; stop.hidden = !b; prog.hidden = !b; };
    start.onclick = () => {
      const errs = validateInputs(suite, values); if (errs.length) { toast(errs[0], 'bad', 6000); return; }
      if (typeof Worker === 'undefined') { toast('Background computation is not available in this browser.', 'bad'); return; }
      const wurl = new URL('../core/worker.js', import.meta.url); try { worker = wurl.origin === location.origin ? new Worker(wurl, { type: 'module' }) : new Worker(URL.createObjectURL(new Blob([`import ${JSON.stringify(wurl.href)};`], { type: 'text/javascript' })), { type: 'module' }); } catch (e) { toast('Background computation could not start: ' + e.message, 'bad', 7000); return; } busy(true); prog.value = 0; stat.textContent = 'Starting…'; const t0 = performance.now();
      worker.onmessage = ev => { const m = ev.data;
        if (m.type === 'progress') { prog.value = m.fraction || 0; stat.textContent = `Step ${m.step?.toLocaleString() ?? '—'} of ${m.steps?.toLocaleString() ?? '—'} · simulated ${fmt(m.simTime, 4)} s · about ${fmtT(m.etaSeconds ?? 0)} left` + (m.maxT ? ` · T max ${fmt(m.maxT, 4)} K` : '') + (m.meanLength ? ` · mean length ${fmt(m.meanLength * 1e9, 4)} nm` : ''); }
        else if (m.type === 'done') { busy(false); worker.terminate(); worker = null; last = finishRun(suite, id, values, m.result, m.runtime, { origin: { ...origin } }); renderResults(last); stat.textContent = `Full-resolution run finished in ${fmtT((performance.now() - t0) / 1000)}. Results are shown on the right.`; toast('Full-resolution run finished', 'ok', 6000); }
        else if (m.type === 'error') { busy(false); worker.terminate(); worker = null; stat.textContent = 'The run stopped: ' + m.message; toast(m.message, 'bad', 8000); } };
      worker.onerror = ev => { busy(false); stat.textContent = 'The background solver could not start: ' + (ev.message || 'unknown error'); };
      worker.postMessage({ type: 'run', id, inputs: { ...values }, hifi: { ...hv }, bus: Object.fromEntries(Object.entries(project.bus).map(([k, v]) => [k, v.value])) });
    };
    stop.onclick = () => { if (worker) { worker.postMessage({ type: 'cancel' }); stat.textContent = 'Stopping at the next time step…'; } };
    showPlan();
    return h('details.group.hifi', h('summary', icon('layers', 15), ' ', H.label), h('p.help', 'Runs the simulation at the full specified resolution in the background on this device. You can keep using the app while it computes; keep this page open until it finishes.'),
      h('div.fields', H.inputs.map(i => h('div.field', h('label', i.label), h('div.numrow', h('input', { type: 'number', value: hv[i.key], min: i.min, max: i.max, step: 'any', oninput: e => { hv[i.key] = Math.min(i.max, Math.max(i.min, +e.target.value || i.value)); showPlan(); } }), i.unit ? h('span.unit', i.unit) : null)))), plan, h('div.runbar', start, stop, prog), stat);
  }
  async function run() {
    if (state.running) return; state.running = true; resBox.classList.add('busy'); if (!resBox.children.length) resBox.replaceChildren(h('div.card', h('p.muted', 'Solving…')));
    await new Promise(r => setTimeout(r, 20));
    try { last = await runSuite(id, values, { origin: { ...origin } }); renderResults(last); }
    catch (e) { resBox.replaceChildren(h('div.card.error', h('h3', icon('warn', 18), ' The simulation could not run'), h('p', 'Check the highlighted inputs and try again.'), h('pre', e.message))); }
    state.running = false; resBox.classList.remove('busy');
  }
  const fig = (title, body, explain) => h('figure.card.fig', h('figcaption', h('h4', title)), body, explain ? h('p.explain', icon('info', 14), ' ', explain) : null);
  function geometryCard(g, r) {
    const sel = h('select', { 'aria-label': 'Export format' }), box = h('div.row');
    import('../geometry/importer.js').then(m => { const ex = Object.entries(m.EXPORTERS).filter(([, e]) => e.accepts.includes(g.kind) || m.convert); ex.sort((a, b) => (b[1].accepts.includes(g.kind) ? 1 : 0) - (a[1].accepts.includes(g.kind) ? 1 : 0)); sel.append(...ex.map(([k, e]) => h('option', { value: k }, `${e.name} (${e.ext})${e.accepts.includes(g.kind) ? '' : ' — converted'}`))); box.append(sel, h('button.btn.sm', { onclick: () => { try { const o = m.exportGeometry(g, sel.value, {}); download(o.name || `suite${id}_geometry${m.EXPORTERS[sel.value].ext}`, o.bytes || o.text, o.mime); if (o.lost && o.lost.length) toast('Exported. Not preserved: ' + o.lost.join('; '), 'info', 7000); else toast('Geometry exported', 'ok'); } catch (e) { toast('Export failed: ' + e.message, 'bad', 6000); } } }, icon('download', 14), ' Export geometry'), h('button.btn.sm.ghost', { onclick: () => { project.geometry.unshift({ name: `Suite ${id}: ${g.title || g.kind}`, at: r.at, geometry: JSON.parse(JSON.stringify(g, (k, v) => ArrayBuffer.isView(v) ? Array.from(v) : v)) }); project.geometry = project.geometry.slice(0, 12); saveProject(); toast('Sent to Geometry Studio', 'ok'); } }, icon('cube', 14), ' Send to Geometry Studio')); }).catch(() => box.append(h('button.btn.sm', { onclick: () => download(`suite${id}_geometry.json`, JSON.stringify(g, (k, v) => ArrayBuffer.isView(v) ? Array.from(v) : v), 'application/json') }, 'Export geometry (JSON)')));
    return fig(g.title || 'Generated geometry', h('div', viewer(g), box), g.frames ? 'Geometry computed by the model. Use “Play growth” to watch it evolve over simulation time — it is not a fixed drawing.' : 'Geometry generated from the model parameters. It can be exported for meshing, fabrication or further simulation.');
  }
  function renderResults(r) {
    const lvl = r.credibility.level, bad = /Review needed/.test(lvl), eqBy = {}; for (const e of r.equations || []) (eqBy[e.kind] ||= []).push(e);
    const recs = [...(r.recommendations || [])].sort((a, b) => 'hml'.indexOf(a.priority[0]) - 'hml'.indexOf(b.priority[0]));
    const practice = h('div.practice', h('p.muted', 'Looking up current published practice…'));
    live.literature(live.SUITE_TOPICS[id] + ' ' + (values.material || project.bus['material.formula']?.value || '')).then(res => practice.replaceChildren(h('h4', 'Current published practice ', h('span.muted', `(${res.cached ? 'cached ' + live.age(res.at) : 'fetched now'} · OpenAlex)`)), h('ul.refs', res.data.slice(0, 5).map(w => h('li', h('a', { href: w.url, target: '_blank', rel: 'noopener' }, w.title), h('span.muted', ` — ${w.venue || 'source n/a'}, ${w.year || ''}`)))), h('p.hint', 'Recent literature for this topic, retrieved live so you can check each recommendation against the newest evidence.'))).catch(() => practice.replaceChildren(h('p.hint', 'Offline: live literature for this topic will load when a connection is available. Recommendations above are computed locally from this run.')));
    resBox.replaceChildren(
      h('div.card.summary', h('div.row.between', h('h3', 'What this run shows'), h('span.muted', `${fmt(r.runtime_ms, 3)} ms · ${new Date(r.at).toLocaleString()}`)), h('p', r.summary),
        h('div.cred', { class: bad ? 'bad' : '' }, icon(bad ? 'warn' : 'info', 16), h('div', h('b', 'Evidence status: ', lvl), h('span.muted', ' — a converged calculation is not proof of physical accuracy. See the Credibility tab for what is still missing.')))),
      (r.warnings || []).length ? h('div.card.warnbox', h('h4', icon('warn', 16), ' Validity limits for this run'), h('ul', r.warnings.map(w => h('li', w)))) : null,
      h('section', h('h3.sec', 'Key results'), h('div.kpis', (r.kpis || []).map(k => h('div.kpi', { class: k.good === true ? 'good' : k.good === false ? 'poor' : '' }, h('div.kl', k.label, k.good === true ? h('span.st.good', icon('check', 12), ' within target') : k.good === false ? h('span.st.poor', icon('warn', 12), ' needs attention') : null), h('div.kv', typeof k.value === 'number' ? fmt(k.value) : k.value, k.unit ? h('span.ku', ' ' + k.unit) : null), h('p.explain', k.explain))))),
      r.geometry ? h('section', h('h3.sec', 'Geometry'), geometryCard(r.geometry, r)) : null,
      (r.series || []).length ? h('section', h('h3.sec', 'Charts'), h('div.figs', r.series.map(sr => fig(sr.title, lineChart(sr), sr.explain)))) : null,
      (r.maps || []).length ? h('section', h('h3.sec', 'Field maps'), h('div.figs', r.maps.map(m => fig(m.title, heatmap(m), m.explain)))) : null,
      (r.tables || []).length ? h('section', h('h3.sec', 'Tables'), r.tables.map(t => fig(t.title, h('div', table(t.columns, t.rows), h('button.btn.sm.ghost', { onclick: () => download(`suite${id}_${t.key || 'table'}.csv`, toCSV(t.columns, t.rows), 'text/csv') }, icon('download', 14), ' CSV')), t.explain))) : null,
      h('section', h('h3.sec', 'Equations for this case'), h('div.card.eqs', Object.entries(EQ_KINDS).filter(([k]) => eqBy[k]).map(([k, title]) => h('div.eqgroup', h('h4', title), eqBy[k].map(e => h('div.eq', h('span.eql', e.label), h('code', e.text), e.note ? h('span.eqn', e.note) : null)))))),
      h('section', h('h3.sec', 'Recommendations'), h('div.recs', recs.map(c => h('div.rec.' + c.priority, h('div.rh', h('span.pri', PRI[c.priority]), c.sustainability ? h('span.sus', icon('leaf', 13), ' Sustainability') : null), h('p', c.text), h('p.basis', h('b', 'Why: '), c.basis)))), h('div.card', practice)),
      h('section', h('h3.sec', 'Physical and numerical checks'), h('div.card', h('ul.checks', (r.checks || []).map(c => h('li', { class: c.pass ? 'ok' : 'no' }, icon(c.pass ? 'check' : 'x', 14), h('b', c.name), ' — ', c.detail))), r.numerics ? h('p.numerics', h('b', 'Numerical method: '), r.numerics.method, ...['dof', 'steps', 'iterations', 'residual'].filter(k => r.numerics[k] !== undefined).map(k => ` · ${k} ${fmt(r.numerics[k], 3)}`), r.numerics.notes ? ` · ${r.numerics.notes}` : '') : null)),
      h('section', h('h3.sec', 'Shared with other suites'), h('div.card', h('p.muted', 'These values were published to the project data bus; suites that list them as inputs pick them up automatically.'), h('div.pubs', Object.entries(r.publish || {}).map(([k, v]) => h('span.pub', h('code', k), ' = ', typeof v === 'number' ? fmt(v) : v))))),
      h('section', h('h3.sec', 'Export'), h('div.card.row.wrap',
        h('button.btn.sm', { onclick: () => download(`suite${id}_report.html`, buildReport(s, suite, r, resBox), 'text/html') }, icon('doc', 14), ' Report (HTML)'),
        h('button.btn.sm', { onclick: () => download(`suite${id}_report.doc`, buildReport(s, suite, r, resBox, { word: true }), 'application/msword') }, icon('doc', 14), ' Report (Word)'),
        h('button.btn.sm', { onclick: () => { document.body.classList.add('printing'); window.print(); document.body.classList.remove('printing'); } }, icon('doc', 14), ' Print / PDF'),
        h('button.btn.sm.ghost', { onclick: () => download(`suite${id}_results.json`, JSON.stringify({ suite: id, name: s.name, ...r }, (k, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 1), 'application/json') }, 'Results (JSON)'),
        h('button.btn.sm.ghost', { onclick: () => download(`suite${id}_kpis.csv`, toCSV(['key', 'label', 'value', 'unit', 'explanation'], r.kpis.map(k => [k.key, k.label, k.value, k.unit || '', k.explain])), 'text/csv') }, 'Key results (CSV)'),
        (r.series || []).length ? h('button.btn.sm.ghost', { onclick: () => { const rows = []; for (const sr of r.series) for (const l of sr.lines) for (let k = 0; k < l.x.length; k++) rows.push([sr.title, l.name, l.x[k], l.y[k]]); download(`suite${id}_series.csv`, toCSV(['chart', 'series', 'x', 'y'], rows), 'text/csv'); } }, 'Chart data (CSV)') : null,
        (r.exports || []).map(x => h('button.btn.sm.ghost', { title: x.description || '', onclick: () => download(x.name, x.bytes || x.content, x.mime || 'text/plain') }, icon('download', 14), ' ' + x.name)))));
  }

  // ---------- Sensitivity tab ----------
  function sensTab() {
    const out = h('div'), kpiSel = h('select'), meth = h('select', h('option', { value: 'oat' }, 'Tornado — one input at a time (fast, local)'), h('option', { value: 'morris' }, 'Morris screening (global ranking, interactions flagged)'), h('option', { value: 'sobol' }, 'Sobol indices (variance shares, most thorough)'), h('option', { value: 'mc' }, 'Uncertainty propagation (Monte Carlo)')), frac = h('input', { type: 'number', value: 10, min: 1, max: 50, step: 1 }), prog = h('progress', { max: 1, value: 0, hidden: true });
    const fill = () => { const ks = (last?.kpis || []).filter(k => typeof k.value === 'number'); kpiSel.replaceChildren(...ks.map(k => h('option', { value: k.key }, k.label))); };
    const go = async () => {
      if (!last) await run(); fill(); const key = kpiSel.value, label = kpiSel.selectedOptions[0]?.textContent, f = Math.max(0.01, +frac.value / 100), ctx = { bus: Object.fromEntries(Object.entries(project.bus).map(([k, v]) => [k, v.value])) }, n = numericInputs(suite, values).length;
      if (!key) return; prog.hidden = false; prog.value = 0; out.replaceChildren(h('p.muted', `Evaluating the model around the current case (${n} numeric inputs)…`)); const onProgress = p => prog.value = p;
      try {
        if (meth.value === 'oat') { const r = await oat(suite, values, key, ctx, { frac: f, onProgress }), top = r.items.filter(i => Number.isFinite(i.swing)); out.replaceChildren(h('h4', `Effect of ±${frac.value} % in each input on “${label}” (base ${fmt(r.base)})`), barChart(top.map(i => ({ label: i.label, lo: i.dLo, hi: i.dHi }))), h('div.legend', h('span.lg', h('i', { style: { background: 'var(--s1)' } }), 'input decreased'), h('span.lg', h('i', { style: { background: 'var(--s2)' } }), 'input increased')), h('p.explain', icon('info', 14), ` Bars show how far the result moves when one input changes by ±${frac.value} % and the rest stay fixed. `, top[0] ? h('b', `“${top[0].label}” dominates` + (top[1] ? `, followed by “${top[1].label}”.` : '.')) : '', ' An elasticity of 1 means a 1 % input change gives a 1 % output change; inputs near zero can be fixed at nominal values, the top ones deserve tight control and measurement.'), table(['Input', 'Low', 'High', 'Result at low', 'Result at high', 'Elasticity'], top.map(i => [i.label + (i.unit ? ` (${i.unit})` : ''), i.lo, i.hi, i.yLo, i.yHi, i.elasticity])), h('p.muted', `${r.evaluations} model evaluations.`)); record('oat', key, top.slice(0, 5).map(i => i.label)); }
        else if (meth.value === 'morris') { const r = await morris(suite, values, key, ctx, { frac: f * 2, onProgress }); out.replaceChildren(h('h4', `Morris screening for “${label}” over ±${(+frac.value * 2)} % windows`), barChart(r.items.map(i => ({ label: i.label, value: i.muStar })), { signed: false }), h('p.explain', icon('info', 14), ' μ* is the mean absolute elementary effect: the larger it is, the more the input matters overall. A large σ relative to μ* means the effect depends on other inputs (interaction) or is non-linear — do not tune that input in isolation.'), table(['Input', 'μ* (importance)', 'μ (signed)', 'σ (interaction / non-linearity)'], r.items.map(i => [i.label, i.muStar, i.mu, i.sigma])), h('p.muted', `${r.evaluations} model evaluations.`)); record('morris', key, r.items.slice(0, 5).map(i => i.label)); }
        else if (meth.value === 'sobol') { const r = await sobol(suite, values, key, ctx, { frac: f * 2, onProgress }); out.replaceChildren(h('h4', `Sobol variance decomposition for “${label}”`), barChart(r.items.map(i => ({ label: i.label, value: i.ST })), { signed: false }), h('p.explain', icon('info', 14), ` Total index Sₜ is the share of output variance removed if the input were known exactly; first-order S₁ is its share acting alone. Roughly ${(r.interaction * 100).toFixed(0)} % of the variance comes from interactions between inputs. Estimates from ${r.evaluations} evaluations carry sampling error of a few hundredths.`), table(['Input', 'S₁ (alone)', 'Sₜ (total)'], r.items.map(i => [i.label, i.S1, i.ST]))); record('sobol', key, r.items.slice(0, 5).map(i => i.label)); }
        else { const keys = (last.kpis || []).filter(k => typeof k.value === 'number').map(k => k.key), r = await monteCarlo(suite, values, keys, ctx, { cv: f, onProgress }), st = r.stats[key]; out.replaceChildren(h('h4', `Spread of “${label}” when every numeric input scatters by ±${frac.value} % (1σ)`), st ? histChart(st.values, { xlabel: label }) : h('p', 'No finite samples.'), st ? h('p.explain', icon('info', 14), ` Across ${r.n} sampled cases the result has mean ${fmt(st.mean)}, standard deviation ${fmt(st.std)} (${(st.cov * 100).toFixed(1)} % of the mean), and a 90 % interval of ${fmt(st.p5)} to ${fmt(st.p95)}. If that interval crosses a specification limit, tighten the dominant inputs identified by the tornado or Sobol analysis.`) : null, table(['Result', 'Mean', 'Std dev', '5th pct', 'Median', '95th pct', 'CoV %'], (last.kpis || []).filter(k => r.stats[k.key]).map(k => { const q = r.stats[k.key]; return [k.label + (k.unit ? ` (${k.unit})` : ''), q.mean, q.std, q.p5, q.p50, q.p95, q.cov * 100]; })), r.failed ? h('p.muted', `${r.failed} samples failed to solve and were excluded.`) : null); record('monte-carlo', key, []); }
      } catch (e) { out.replaceChildren(h('div.card.error', e.message)); }
      prog.hidden = true;
    };
    const record = (method, key, top) => { if (project.runs[id]) { project.runs[id].sensitivity = { method, key, top, at: new Date().toISOString() }; saveProject(); } };
    fill();
    return h('div', h('div.card', h('p', 'Find out which inputs control each result, how they interact, and how measurement or process scatter propagates into the prediction. The analysis re-runs this suite’s model around the current case.'), h('div.fields', h('div.field', h('label', 'Result to analyse'), kpiSel), h('div.field', h('label', 'Method'), meth), h('div.field', h('label', 'Variation of each input'), h('div.numrow', frac, h('span.unit', '%')))), h('div.runbar', h('button.btn.primary', { onclick: go }, icon('chart', 16), ' Run analysis'), prog)), h('div.card', out));
  }

  // ---------- knowledge tabs ----------
  const chips = list => h('div.taglist', (list || []).map(x => h('span.tag', x)));
  const modelsTab = () => h('div', h('div.card', h('h3', 'Governing equations and mathematical formulations'), chips(s.eq.formulations)), h('div.card', h('h3', 'Classical and contemporary models'), chips(s.eq.models)), h('div.card', h('h3', 'Hybrid and computational methods'), chips(s.eq.hybrid)), h('div.card', h('h3', 'Engineering outputs'), chips(s.eq.outputs)), s.capabilities ? h('div.card', h('h3', 'Core capabilities'), h('ol.caps', s.capabilities.map(c => h('li', h('b', c.name + ': '), c.text)))) : null, h('div.card.note', icon('info', 16), h('p', 'Only the models physically relevant to the selected analysis are activated in a run; the equations actually solved are listed with every result. Reduced-order models implemented in-browser are named explicitly, and input decks can be exported for high-fidelity external solvers where offered.')));
  const condTab = () => h('div', h('div.card', h('h3', 'Initial conditions'), h('p', s.ibc.ic)), h('div.card', h('h3', s.ibc.bcLabel || 'Boundary conditions'), h('p', s.ibc.bc)), s.ibc.note ? h('div.card.note', icon('info', 16), h('p', s.ibc.note)) : null, h('div.card', h('h3', 'Boundary-condition classes available platform-wide'), table(['Class', 'Typical application'], spec.cross.bcClasses.map(b => [b.cls, b.use]))));
  const numTab = () => h('div', h('div.card', h('h3', 'Numerical settings that matter most for this suite'), h('p', s.num + '.')), h('div.card', h('h3', 'Configuration levels'), table(['Mode', 'How it operates'], spec.cross.numModes.map(m => [m.mode, m.text]))), h('div.card', h('h3', 'Numerical parameter groups'), spec.cross.numGroups.map(g => h('details', h('summary', `${g.id}. ${g.title}`, h('span.muted', '  ' + g.part)), g.text.map(t => h('p', t))))));
  function benchCard() {
    const box = h('div', h('p.muted', 'Loading reference benchmarks…')), prog = h('progress', { max: 1, value: 0, hidden: true });
    const show = res => { const sm = summarise(res); benchStatus[id] = sm; kv.set('bench:' + id, { ...sm, at: Date.now() }); box.replaceChildren(h('p', h('b', `${sm.passed} of ${sm.total} benchmarks reproduced`), ` on this device; worst relative deviation ${(sm.worst * 100).toFixed(2)} %.`), table(['Benchmark', 'Type', 'Reference value', 'This solver', 'Deviation', 'Allowed', 'Result', 'Source'], res.map(r => [h('span', { title: r.note || '' }, r.title), BENCH_KIND[r.kind] || r.kind, r.expected, r.error ? 'error' : r.got, r.relative ? `${(r.dev * 100).toFixed(2)} %` : fmt(r.dev, 3), r.abs !== undefined ? `± ${r.abs}` : `${(r.tol * 100).toFixed(1)} %`, h('span.tag', { class: r.ok ? 'sup-full' : 'sup-inventory' }, r.ok ? 'pass' : 'differs'), r.reference])), h('p.hint', 'Each reference value comes from an exact solution, handbook data or a published benchmark and is independent of this software. Agreement shows the solver reproduces established results for those cases; your own case still needs comparison with your own measurements (Suite 24).')); };
    const go = async () => { prog.hidden = false; const res = await runSuiteBenchmarks(id, p => prog.value = p); prog.hidden = true; show(res); };
    loadBenchmarks(id).then(list => { if (!list.length) { box.replaceChildren(h('p.muted', 'No reference benchmarks are packaged for this suite yet.')); return; } box.replaceChildren(h('p', `${list.length} independent reference results are packaged for this suite.`)); go(); });
    return h('div.card', h('div.row.between.wrap', h('h3', 'Reference benchmarks'), h('button.btn.sm', { onclick: go }, icon('refresh', 14), ' Run again')), prog, box);
  }
  const cvvTab = () => h('div', benchCard(), last ? h('div.card', h('h3', 'Qualification gates for the current run'), h('ul.gates', last.credibility.gates.map(g => h('li.' + g.status, h('span.gs', { pass: 'Passed', partial: 'Partly evidenced', open: 'Open', fail: 'Failed' }[g.status]), h('b', `Gate ${g.id} — ${g.name}`), h('p', g.detail)))), h('p', h('b', 'Evidence status: '), last.credibility.level), h('a.btn.sm', { href: '#/suite/24' }, 'Open Suite 24 to add verification, calibration and validation evidence')) : h('div.card', h('p.muted', 'Run the simulation to see its qualification gates.')),
    h('div.card', h('h3', 'Calibration'), h('p', s.cvv.cal)), h('div.card', h('h3', 'Verification'), h('p', s.cvv.ver)), h('div.card', h('h3', 'Validation'), h('p', s.cvv.val)), h('div.card', h('h3', 'Expected evidence'), h('p', s.cvv.exp)), s.cvv.note ? h('div.card.note', icon('info', 16), h('p', s.cvv.note)) : null);
  const ioTab = () => h('div.two', h('div', h('h3', 'Accepted input data'), s.io.inputs.map(g => h('div.card', h('h4', g.group), h('ul.cols', g.items.map(x => h('li', x))))), h('a.btn.sm', { href: '#/data' }, icon('upload', 14), ' Open the data portal to import case data')), h('div', h('h3', 'Generated output data'), s.io.outputs.map(g => h('div.card', h('h4', g.group), h('ul.cols', g.items.map(x => h('li', x)))))));
  function liveTab() {
    const box = h('div'), q = h('input', { type: 'search', value: live.SUITE_TOPICS[id], 'aria-label': 'Topic' }), src = h('select', ['openalex', 'crossref', 'europepmc', 'zenodo', 'datacite'].map(k => h('option', { value: k }, `${live.CONNECTORS[k].name} — ${live.CONNECTORS[k].kind}`)));
    const go = async (force) => { box.replaceChildren(h('p.muted', 'Querying…')); try { const r = await live.literature(q.value, src.value, { force }); box.replaceChildren(h('p.muted', `${r.data.length} records · ${r.cached ? 'cached ' + live.age(r.at) + (r.stale ? ' (offline copy)' : '') : 'fetched just now'}`), h('ul.refs', r.data.map(w => h('li', h('a', { href: w.url, target: '_blank', rel: 'noopener' }, w.title || '(untitled)'), h('span.muted', ` — ${w.venue || ''} ${w.year || ''}${w.cited != null ? ` · cited ${w.cited}` : ''}`))))); } catch (e) { box.replaceChildren(h('p.muted', live.online() ? `The source did not respond (${e.message}).` : 'You are offline and nothing is cached for this query yet.')); } };
    setTimeout(() => go(false), 0);
    return h('div', h('div.card', h('h3', 'Live reference data for this suite'), h('p', 'Queries go straight from your browser to the source, so results are current whenever you are online and the last answers stay available offline.'), h('div.row.wrap', q, src, h('button.btn.sm', { onclick: () => go(true) }, icon('refresh', 14), ' Refresh')), box), h('div.card', h('h3', 'More trusted sources'), h('div.taglist', live.PORTALS.map(p => h('a.tag', { href: p.url(project.bus['material.formula']?.value || 'ZnO'), target: '_blank', rel: 'noopener', title: p.kind }, p.name))), h('a.btn.sm.ghost', { href: '#/live' }, 'All live connections')));
  }

  renderForm();
  root.append(tabs([
    { label: 'Simulate', render: () => h('div.workspace', h('aside.card.formcard', formBox), resBox) },
    { label: 'Sensitivity & uncertainty', render: sensTab },
    { label: 'Models & equations', render: modelsTab },
    { label: 'Initial & boundary conditions', render: condTab },
    { label: 'Numerics', render: numTab },
    { label: 'Credibility (CVV)', render: cvvTab, fresh: true },
    { label: 'Inputs & outputs', render: ioTab },
    { label: 'Live data', render: liveTab },
  ]));
  run();
}
