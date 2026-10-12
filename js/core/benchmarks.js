// Reference benchmarks: independent analytical results, handbook data and published numerical references that each
// suite is compared with on the user's device.
import { loadSuite, defaultInputs } from './engine.js';

const cache = new Map();
export async function loadBenchmarks(id) {
  if (!cache.has(id)) cache.set(id, import(`../../data/benchmarks/s${String(id).padStart(2, '0')}.js`).then(m => m.default).catch(() => []));
  return cache.get(id);
}
export const KIND = { analytical: 'Exact / analytical solution', 'reference-data': 'Handbook or measured reference data', 'numerical-reference': 'Published numerical benchmark', conservation: 'Exact invariant', calibration: 'Calibration case (plausibility range; tighten with your own data)' };
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
