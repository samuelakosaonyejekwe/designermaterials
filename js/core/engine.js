// Suite orchestration: lazy module loading, input resolution against the shared bus, execution,
// credibility gating, cross-suite publication and the full-case pipeline.
import { project, busValues, publish, recordRun, setInputs } from './store.js';
import { visible } from './sensitivity.js';

let specP = null; export let SPEC = null;
export const loadSpec = () => specP || (specP = fetch('data/spec.json').then(r => r.json()).then(s => (SPEC = s)));
const mods = new Map();
export async function loadSuite(id) { if (!mods.has(id)) mods.set(id, import(`../suites/s${String(id).padStart(2, '0')}.js`).then(m => m.default)); return mods.get(id); }

// Order in which a complete case flows through the platform (forward engineering), then reverse engineering.
export const PIPELINE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29];
export const STAGES = [
  { name: 'Define material and architecture', suites: [1, 2, 3, 4, 5], text: 'Composition, crystal structure, nano-geometry and unit-cell architecture.' },
  { name: 'Fundamental physics and behaviour', suites: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], text: 'Chemistry, transport, growth, thermal, electrical, optical, magnetic, mechanical, wave and atomistic physics.' },
  { name: 'Multiscale prediction and control', suites: [18, 19, 20, 22], text: 'Effective properties, component-scale performance, tuning, degradation and lifetime.' },
  { name: 'Fabrication and characterisation', suites: [21, 23], text: 'Manufacturing conditions, predicted production outcomes and measured properties.' },
  { name: 'Validation, optimisation and deployment', suites: [24, 25, 26, 27, 28], text: 'Credibility, inverse design, industrialisation, safety, economics and lifecycle.' },
  { name: 'Reverse engineering', suites: [29], text: 'Start from a measured material; reconstruct geometry, composition, properties and plausible fabrication history.' },
];

export function defaultInputs(suite) { const o = {}; for (const i of suite.inputs) o[i.key] = i.value; return o; }
// Saved values > bus-linked values > defaults. Returns values plus where each came from.
export function resolveInputs(suite, { useBus = true } = {}) {
  const saved = project.inputs[suite.id] || {}, bus = project.bus, values = {}, origin = {};
  for (const i of suite.inputs) {
    if (i.key in saved) { values[i.key] = saved[i.key]; origin[i.key] = 'user'; }
    else if (useBus && i.bus && bus[i.bus] && bus[i.bus].from !== suite.id && fits(i, bus[i.bus].value)) { values[i.key] = coerce(i, bus[i.bus].value); origin[i.key] = 'bus'; }
    else { values[i.key] = i.value; origin[i.key] = 'default'; }
  }
  return { values, origin };
}
const fits = (i, v) => (i.type || 'number') === 'number' ? typeof v === 'number' && Number.isFinite(v) : typeof v === 'string';
const coerce = (i, v) => (i.type || 'number') === 'number' ? Math.min(i.max, Math.max(i.min, v)) : v;
export function validateInputs(suite, values) {
  const errs = [];
  for (const i of suite.inputs) { if (!visible(i, values)) continue; const v = values[i.key]; if ((i.type || 'number') === 'number') { if (typeof v !== 'number' || !Number.isFinite(v)) errs.push(`${i.label}: enter a number`); else if (v < i.min || v > i.max) errs.push(`${i.label}: ${v} is outside the supported range ${i.min} – ${i.max} ${i.unit || ''}`); } }
  return errs;
}

// Eight qualification gates. Only evidence the platform can actually observe is scored automatically; gates that need
// experimental evidence stay open until the user records it in Suite 24.
export const benchStatus = {};
// Nearest published experiment for the current inputs of each suite: { inputs: JSON, title, source, distance, ok, dev }
export const expStatus = {};
export function assessGates(suite, values, result, origin) {
  const checks = result.checks || [], failed = checks.filter(c => !c.pass), assumed = Object.values(origin || {}).filter(o => o === 'default').length, total = Object.keys(origin || {}).length || 1;
  const cred = project.bus, g = [];
  g.push({ id: 1, name: 'Input data quality', status: assumed === 0 ? 'pass' : assumed / total > 0.6 ? 'open' : 'partial', detail: `${total - assumed} of ${total} inputs supplied by the user or linked from other suites; ${assumed} still at library defaults. Units are fixed per field.` });
  g.push({ id: 2, name: 'Mathematical and model consistency', status: failed.length ? 'fail' : 'pass', detail: failed.length ? `${failed.length} physical/admissibility check(s) failed: ${failed.map(c => c.name).join('; ')}` : `${checks.length} conservation / admissibility checks passed for this run.` });
  const bm = benchStatus[suite.id];
  g.push({ id: 3, name: 'Code verification', status: bm ? (bm.passed === bm.total ? 'pass' : 'partial') : 'partial', detail: bm ? `${bm.passed} of ${bm.total} independent reference benchmarks reproduced on this device (worst deviation ${(bm.worst * 100).toFixed(2)} %). They cover the suite's models, not your specific case.` : 'This suite ships independent reference benchmarks (analytical solutions and handbook data). Run them from the Credibility tab to verify the solvers on this device.' });
  const gci = cred['cred.gci'];
  g.push({ id: 4, name: 'Solution verification', status: gci ? (gci.value < 0.05 ? 'pass' : 'partial') : (result.numerics && result.numerics.residual !== undefined ? 'partial' : 'open'), detail: gci ? `Grid-convergence index ${(gci.value * 100).toFixed(2)} % recorded in Suite 24.` : 'Solver residuals are reported, but no mesh / time-step convergence study has been recorded for this case.' });
  const rec = project.runs[suite.id] || {}, cal = rec.calibration, val = rec.validation && rec.validation.inputs === JSON.stringify(values) ? rec.validation : null;
  g.push({ id: 5, name: 'Parameter calibration', status: cal || cred['cred.calibrated'] || cred['cred.cal_rmse'] ? 'pass' : 'open', detail: cal ? `“${cal.inputLabel}” was fitted to ${cal.value} so that “${cal.kpiLabel}” reproduces your measured ${cal.measured} (${new Date(cal.at).toLocaleDateString()}).` : (cred['cred.calibrated'] || cred['cred.cal_rmse']) ? 'Calibrated parameters recorded in Suite 24.' : 'Model parameters are handbook / literature-typical values, not calibrated to your measurements. Use “Calibrate to a measurement” below.' });
  const lit = expStatus[suite.id] && expStatus[suite.id].inputs === JSON.stringify(values) && expStatus[suite.id].ok && expStatus[suite.id].distance <= 0.15 ? expStatus[suite.id] : null;
  g.push({ id: 6, name: 'Independent experimental validation', status: val ? (val.passed === val.n ? 'pass' : 'fail') : lit ? 'pass' : cred['cred.rmse'] ? 'partial' : 'open', detail: val ? `${val.passed} of ${val.n} predicted quantities agree with your measurements within the stated acceptance (${val.acceptance} % or 2 × measurement uncertainty), recorded ${new Date(val.at).toLocaleDateString()} for exactly these inputs.` : lit ? `This case lies within ${(lit.distance * 100).toFixed(0)} % of a published experiment that the solver reproduces to ${(lit.dev * 100).toFixed(1)} %: ${lit.title} (${lit.source}). Your own measurements would strengthen this further.` : cred['cred.rmse'] ? `Validation error recorded in Suite 24 (RMSE ${Number(cred['cred.rmse'].value).toPrecision(3)}). Check it covers this operating range.` : 'No measurements have been compared with this prediction yet. Use “Validate with my measurements” below.' });
  g.push({ id: 7, name: 'Predictive robustness', status: project.runs[suite.id]?.sensitivity ? 'partial' : 'open', detail: project.runs[suite.id]?.sensitivity ? 'A sensitivity / uncertainty study has been run for this suite.' : 'Run the sensitivity and uncertainty tabs to see how input scatter moves the result.' });
  g.push({ id: 8, name: 'Application-specific qualification', status: 'open', detail: 'Requires documented review for the intended engineering use; never granted automatically.' });
  let level = 'Theoretical prediction — not verified';
  if (bm && bm.passed === bm.total && g[1].status === 'pass') level = `Solver benchmarked against ${bm.total} reference results — case not yet validated`;
  if (g[1].status === 'pass' && g[3].status !== 'open') level = 'Numerically checked for this case';
  if (g[1].status === 'pass' && g[4].status === 'pass') level = 'Calibrated against identified datasets';
  if (g[1].status === 'pass' && g[4].status === 'pass' && g[5].status !== 'open') level = 'Compared with independent measurements';
  if (lit && g[1].status === 'pass') level = `Validated against a published experiment near this operating point (${lit.source})`;
  if (val && val.passed === val.n && g[1].status === 'pass') level = `Experimentally validated at this operating point (${val.n} measured ${val.n === 1 ? 'quantity' : 'quantities'})`;
  if (val && val.passed < val.n) level = 'Disagrees with your measurements — calibrate or revise the model';
  if (g[1].status === 'fail') level = 'Review needed — a physical or validity check failed';
  return { gates: g, level };
}

export async function runSuite(id, values, { origin = {}, record = true } = {}) {
  const suite = await loadSuite(id); await loadSpec();
  const errs = validateInputs(suite, values); if (errs.length) { const e = new Error(errs.join('\n')); e.validation = errs; throw e; }
  const ctx = { bus: busValues(), spec: SPEC.suites[id - 1], project: { name: project.name } }, t0 = performance.now();
  const result = await suite.run({ ...values }, ctx);
  return finishRun(suite, id, values, result, performance.now() - t0, { origin, record });
}
// Shared post-processing for interactive and background (worker) runs.
export function finishRun(suite, id, values, result, runtime, { origin = {}, record = true } = {}) {
  const at = new Date().toISOString();
  result.runtime_ms = runtime; result.at = at; result.inputs = { ...values }; result.origin = origin;
  result.credibility = assessGates(suite, values, result, origin);
  if (record) {
    setInputs(id, values); publish(id, result.publish, at);
    recordRun(id, { at, runtime, inputs: values, summary: result.summary, kpis: (result.kpis || []).map(k => ({ key: k.key, label: k.label, value: k.value, unit: k.unit, explain: k.explain, good: k.good })), equations: result.equations, recommendations: result.recommendations, warnings: result.warnings, checks: result.checks, level: result.credibility.level, numerics: result.numerics, sensitivity: project.runs[id]?.sensitivity || null, calibration: project.runs[id]?.calibration || null, validation: project.runs[id]?.validation && project.runs[id].validation.inputs === JSON.stringify(values) ? project.runs[id].validation : null });
  }
  return result;
}
// Run every suite in pipeline order so each consumes what earlier suites published.
export async function runCase({ suites = PIPELINE, onStep, passes = 1 } = {}) {
  const log = [];
  for (let p = 0; p < passes; p++) for (const id of suites) {
    const suite = await loadSuite(id), { values, origin } = resolveInputs(suite); let entry;
    try { const r = await runSuite(id, values, { origin }); entry = { id, ok: true, runtime: r.runtime_ms, summary: r.summary, linked: Object.values(origin).filter(o => o === 'bus').length, published: Object.keys(r.publish || {}).length, warnings: (r.warnings || []).length, failedChecks: (r.checks || []).filter(c => !c.pass).length, level: r.credibility.level }; }
    catch (e) { entry = { id, ok: false, error: e.message }; }
    log.push(entry); if (onStep) onStep(entry, log.length, suites.length * passes); await new Promise(r => setTimeout(r, 0));
  }
  return log;
}
// Which suites feed / consume which bus keys (for the connectivity map)
export async function wiring() {
  const edges = [], prod = {}, cons = {};
  for (const id of PIPELINE) { let s; try { s = await loadSuite(id); } catch { continue; } for (const i of s.inputs) if (i.bus) (cons[i.bus] ||= new Set()).add(id); }
  for (const k in project.bus) (prod[k] ||= new Set()).add(project.bus[k].from);
  for (const k in cons) for (const to of cons[k]) for (const from of prod[k] || []) if (from && from !== to) edges.push({ from, to, key: k });
  return { edges, consumes: cons, produced: prod };
}
