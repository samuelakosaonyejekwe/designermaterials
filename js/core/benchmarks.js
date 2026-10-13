// Reference benchmarks: independent analytical results, handbook data and published numerical references that each
// suite is compared with on the user's device.
import { loadSuite, defaultInputs } from './engine.js';

const cache = new Map();
export async function loadBenchmarks(id) {
  if (!cache.has(id)) cache.set(id, import(`../../data/benchmarks/s${String(id).padStart(2, '0')}.js`).then(m => m.default).catch(() => []));
  return cache.get(id);
}
export const KIND = { experiment: 'Published experiment', analytical: 'Exact / analytical solution', 'reference-data': 'Handbook or measured reference data', 'numerical-reference': 'Published numerical benchmark', conservation: 'Exact invariant', calibration: 'Calibration case (plausibility range; tighten with your own data)' };
export function compare(b, got) {
  if (typeof got !== 'number' || !Number.isFinite(got)) return { ok: false, dev: NaN, got };
  const abs = Math.abs(got - b.expected), rel = b.expected !== 0 ? abs / Math.abs(b.expected) : abs;
  return { ok: b.abs !== undefined ? abs <= b.abs : rel <= b.tol, dev: b.abs !== undefined ? abs : rel, relative: b.abs === undefined, got };
}
export async function runBenchmark(b, suite) {
  suite = suite || await loadSuite(b.suite);
  try {
    const r = await suite.run({ ...defaultInputs(suite), ...b.inputs }, { bus: {}, benchmark: true });
    const got = b.pub ? r.publish?.[b.pub] : r.kpis.find(k => k.key === b.kpi)?.value;
    return { ...b, ...compare(b, got), unit: b.pub ? '' : r.kpis.find(k => k.key === b.kpi)?.unit || '' };
  } catch (e) { return { ...b, ok: false, dev: NaN, got: NaN, error: e.message }; }
}
export async function runSuiteBenchmarks(id, onProgress) {
  const list = await loadBenchmarks(id), suite = list.length ? await loadSuite(id) : null, out = [];
  for (const b of list) { out.push(await runBenchmark(b, suite)); if (onProgress) onProgress(out.length / list.length); await new Promise(r => setTimeout(r, 0)); }
  return out;
}
export const summarise = res => ({ total: res.length, passed: res.filter(r => r.ok).length, worst: res.reduce((m, r) => r.relative && Number.isFinite(r.dev) ? Math.max(m, r.dev) : m, 0) });

// How close a user's case is to the published experiments packaged for a suite: for each experiment of the same analysis
// mode, the largest relative difference over the inputs that define its operating point.
export async function nearestExperiments(id, values, suite) {
  const list = (await loadBenchmarks(id)).filter(b => b.kind === 'experiment' && (!b.inputs?.mode || b.inputs.mode === values.mode)), def = suite ? Object.fromEntries(suite.inputs.map(i => [i.key, i])) : {};
  return list.map(b => { const keys = (b.match && b.match.length ? b.match : Object.keys(b.inputs || {})).filter(k => k !== 'mode'); let worst = 0, diff = [];
    for (const k of keys) { const e = b.inputs?.[k] ?? def[k]?.value, v = values[k]; if (typeof e === 'number' && typeof v === 'number') { const span = def[k] ? Math.abs(def[k].max - def[k].min) : 0, d = e !== 0 ? Math.abs(v - e) / Math.abs(e) : span ? Math.abs(v - e) / span : Math.abs(v - e); worst = Math.max(worst, d); if (d > 0.02) diff.push({ key: k, label: def[k]?.label || k, yours: v, theirs: e, unit: def[k]?.unit || '' }); } else if (e !== undefined && v !== e) { worst = Infinity; diff.push({ key: k, label: def[k]?.label || k, yours: v, theirs: e, unit: '' }); } }
    return { ...b, distance: worst, diff }; }).sort((a, b) => a.distance - b.distance);
}
