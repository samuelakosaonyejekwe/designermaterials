// Sensitivity and uncertainty analysis that works on any suite: one-at-a-time (tornado), Morris screening,
// Sobol variance decomposition (Saltelli estimator) and Monte Carlo / Latin-hypercube uncertainty propagation.
import { rng, lhs, mean, std, quantile, normInv, clamp } from './num.js';

export const numericInputs = (suite, values) => suite.inputs.filter(i => (i.type || 'number') === 'number' && i.sens !== false && i.max > i.min && visible(i, values));
export const visible = (i, values) => !i.showIf || Object.entries(i.showIf).every(([k, v]) => (Array.isArray(v) ? v : [v]).includes(values[k]));
const kpiOf = (r, key) => { const k = r && r.kpis && r.kpis.find(x => x.key === key); return k && typeof k.value === 'number' && Number.isFinite(k.value) ? k.value : NaN; };
async function evalAt(suite, base, over, key, ctx) { try { return kpiOf(await suite.run({ ...base, ...over }, ctx), key); } catch { return NaN; } }
const span = (i, base, frac) => { const v = base[i.key]; if (i.log && v > 0) return [clamp(v / (1 + frac), i.min, i.max), clamp(v * (1 + frac), i.min, i.max)]; const d = Math.abs(v) > 1e-300 ? Math.abs(v) * frac : (i.max - i.min) * frac * 0.5; return [clamp(v - d, i.min, i.max), clamp(v + d, i.min, i.max)]; };
const fromUnit = (i, u, base, frac) => { const [lo, hi] = frac ? span(i, base, frac) : [i.min, i.max]; return i.log && lo > 0 ? lo * (hi / lo) ** u : lo + (hi - lo) * u; };
const tick = () => new Promise(r => setTimeout(r, 0));

// Local one-at-a-time: ±frac around the base point. Returns normalised sensitivity (elasticity) and swing.
export async function oat(suite, base, key, ctx, { frac = 0.1, onProgress } = {}) {
  const ins = numericInputs(suite, base), y0 = await evalAt(suite, base, {}, key, ctx), out = [];
  for (let n = 0; n < ins.length; n++) {
    const i = ins[n], [lo, hi] = span(i, base, frac), yl = await evalAt(suite, base, { [i.key]: lo }, key, ctx), yh = await evalAt(suite, base, { [i.key]: hi }, key, ctx);
    const dx = (hi - lo) / (Math.abs(base[i.key]) || 1), el = Number.isFinite(yl) && Number.isFinite(yh) && dx > 0 && y0 !== 0 ? ((yh - yl) / y0) / dx : NaN;
    out.push({ key: i.key, label: i.label, unit: i.unit, lo, hi, yLo: yl, yHi: yh, dLo: yl - y0, dHi: yh - y0, swing: Math.abs(yh - yl), elasticity: el });
    if (onProgress) onProgress((n + 1) / ins.length); if (n % 3 === 2) await tick();
  }
  out.sort((a, b) => (b.swing || 0) - (a.swing || 0)); return { base: y0, items: out, evaluations: 1 + 2 * ins.length, frac };
}
// Morris elementary effects over the full declared ranges (or ±frac window)
export async function morris(suite, base, key, ctx, { r = 8, levels = 6, frac = 0, seed = 5, onProgress } = {}) {
  const ins = numericInputs(suite, base), k = ins.length, rand = rng(seed), delta = levels / (2 * (levels - 1)), ee = ins.map(() => []);
  for (let t = 0; t < r; t++) {
    let u = ins.map(() => rand.int(levels / 2) / (levels - 1)); const order = ins.map((_, i) => i).sort(() => rand() - 0.5);
    const point = uu => Object.fromEntries(ins.map((i, j) => [i.key, fromUnit(i, uu[j], base, frac)])); let y = await evalAt(suite, base, point(u), key, ctx);
    for (const j of order) { const u2 = [...u]; u2[j] = u[j] + delta > 1 ? u[j] - delta : u[j] + delta; const y2 = await evalAt(suite, base, point(u2), key, ctx); if (Number.isFinite(y) && Number.isFinite(y2)) ee[j].push((y2 - y) / (u2[j] - u[j])); u = u2; y = y2; }
    if (onProgress) onProgress((t + 1) / r); await tick();
  }
  const items = ins.map((i, j) => ({ key: i.key, label: i.label, muStar: mean(ee[j].map(Math.abs)) || 0, mu: mean(ee[j]) || 0, sigma: ee[j].length > 1 ? std(ee[j]) : 0, n: ee[j].length })).sort((a, b) => b.muStar - a.muStar);
  return { items, evaluations: r * (k + 1) };
}
// Sobol first-order and total indices (Saltelli 2010 / Jansen estimators)
export async function sobol(suite, base, key, ctx, { n = 96, frac = 0.2, seed = 9, onProgress } = {}) {
  const ins = numericInputs(suite, base), k = ins.length, A = lhs(n, k, rng(seed)), B = lhs(n, k, rng(seed + 101));
  const run = async row => evalAt(suite, base, Object.fromEntries(ins.map((i, j) => [i.key, fromUnit(i, row[j], base, frac)])), key, ctx);
  const yA = [], yB = []; for (let s = 0; s < n; s++) { yA.push(await run(A[s])); yB.push(await run(B[s])); if (s % 8 === 7) await tick(); }
  const okRows = yA.map((v, s) => Number.isFinite(v) && Number.isFinite(yB[s])), all = yA.concat(yB).filter(Number.isFinite), V = std(all) ** 2 || 1e-300, items = [];
  for (let j = 0; j < k; j++) {
    let s1 = 0, st = 0, cnt = 0;
    for (let s = 0; s < n; s++) { if (!okRows[s]) continue; const row = [...A[s]]; row[j] = B[s][j]; const y = await run(row); if (!Number.isFinite(y)) continue; s1 += yB[s] * (y - yA[s]); st += (yA[s] - y) ** 2; cnt++; }
    items.push({ key: ins[j].key, label: ins[j].label, S1: cnt ? clamp(s1 / cnt / V, -0.2, 1.2) : NaN, ST: cnt ? clamp(st / (2 * cnt) / V, 0, 1.5) : NaN });
    if (onProgress) onProgress((j + 1) / k); await tick();
  }
  items.sort((a, b) => (b.ST || 0) - (a.ST || 0)); return { items, variance: V, mean: mean(all), evaluations: n * (k + 2), frac, interaction: Math.max(0, 1 - items.reduce((a, i) => a + Math.max(0, i.S1 || 0), 0)) };
}
// Uncertainty propagation: each numeric input ~ Normal(base, cv*|base|) truncated to its bounds
export async function monteCarlo(suite, base, keys, ctx, { n = 200, cv = 0.05, seed = 21, onProgress } = {}) {
  const ins = numericInputs(suite, base), U = lhs(n, ins.length, rng(seed)), out = Object.fromEntries(keys.map(k => [k, []])); let failed = 0;
  for (let s = 0; s < n; s++) {
    const over = Object.fromEntries(ins.map((i, j) => [i.key, clamp(base[i.key] + normInv(clamp(U[s][j], 1e-6, 1 - 1e-6)) * cv * (Math.abs(base[i.key]) || (i.max - i.min) * 0.1), i.min, i.max)]));
    let r = null; try { r = await suite.run({ ...base, ...over }, ctx); } catch { failed++; }
    for (const k of keys) { const v = kpiOf(r, k); if (Number.isFinite(v)) out[k].push(v); }
    if (s % 10 === 9) { if (onProgress) onProgress((s + 1) / n); await tick(); }
  }
  return { n, cv, failed, stats: Object.fromEntries(keys.map(k => { const a = out[k]; return [k, a.length ? { values: a, mean: mean(a), std: std(a), p5: quantile(a, 0.05), p50: quantile(a, 0.5), p95: quantile(a, 0.95), cov: Math.abs(std(a) / (mean(a) || 1)) } : null]; })) };
}
