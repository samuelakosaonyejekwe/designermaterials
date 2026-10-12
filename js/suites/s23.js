// Suite 23 — Experimental Characterisation and Metrology.
// Reduction of measured data to material properties with uncertainty: XRD (peak fitting, phase identification,
// lattice refinement, Scherrer / Williamson–Hall), particle sizing (DLS cumulants, counting statistics), gas
// adsorption (BET, t-plot, Kelvin pore size), optical and vibrational spectroscopy (Tauc, Beer–Lambert, LSPR,
// Raman mode matching), indentation (Oliver–Pharr, Hertz), electrical (four-point probe, van der Pauw, Hall,
// impedance), thermal / mechanical / magnetic reductions and GUM uncertainty budgets.
// Every default data table is seeded synthetic data from the matching forward model plus noise.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import { CRYSTALS, powderXRD } from '../geometry/atomic.js';

const F = N.fmt, PI = Math.PI, D2R = PI / 180, LN2 = Math.LN2;
const fin = (v, d = 0) => (Number.isFinite(v) ? v : d);
const SYN = new Set();
export const SYN_MSG = 'synthetic demonstration data — replace with measurements';
export const registerSynthetic = t => { SYN.add(String(t).trim()); return t; };
export const isSynthetic = t => SYN.has(String(t ?? '').trim());
const syn = registerSynthetic;
export const csv = (head, rows) => head.join(',') + '\n' + rows.map(r => r.map(v => (typeof v === 'number' ? +v.toPrecision(6) : v)).join(',')).join('\n');

// ---------- defensive table parsing ----------
export function parseTable(text, minCols, minRows, label) {
  const lines = String(text ?? '').split(/\r?\n/).map(l => l.trim()).filter(l => l && !/^(#|\/\/|%)/.test(l)).slice(0, 20001);
  const rows = []; let head = null;
  for (const l of lines) {
    let c = l.split(/[,;\t]/).map(s => s.trim()); if (c.length === 1 && /\s/.test(l)) c = l.split(/\s+/);
    const v = c.map(s => (s === '' ? NaN : Number(s)));
    if (c.length >= minCols && v.slice(0, minCols).every(Number.isFinite)) rows.push(v); else if (!head && !rows.length) head = c;
  }
  if (rows.length < minRows) throw new Error(`${label}: at least ${minRows} numeric rows with ${minCols} column(s) are needed (found ${rows.length}). Paste CSV/TSV text with one header line.`);
  return { head: head || [], rows };
}
function parseNamed(text) { // name,value,u rows
  const out = {};
  for (const l of String(text ?? '').split(/\r?\n/)) { const c = l.trim().split(/[,;\t]/).map(s => s.trim()); if (c.length >= 2 && c[0] && Number.isFinite(Number(c[1])) && c[1] !== '') out[c[0].toLowerCase()] = [Number(c[1]), Math.abs(fin(Number(c[2]), 0))]; }
  return out;
}

// ---------- small statistics ----------
export function linfit(x, y) {
  const n = x.length, mx = N.mean(x), my = N.mean(y); let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxx += (x[i] - mx) ** 2; sxy += (x[i] - mx) * (y[i] - my); syy += (y[i] - my) ** 2; }
  const slope = sxx > 0 ? sxy / sxx : 0, intercept = my - slope * mx; let ss = 0; for (let i = 0; i < n; i++) ss += (y[i] - intercept - slope * x[i]) ** 2;
  const s2 = n > 2 ? ss / (n - 2) : 0, seSlope = sxx > 0 ? Math.sqrt(s2 / sxx) : 0, seInt = Math.sqrt(s2 * (1 / Math.max(1, n) + (sxx > 0 ? mx * mx / sxx : 0)));
  return { slope, intercept, seSlope, seInt, r2: syy > 0 ? Math.max(0, 1 - ss / syy) : 1, s: Math.sqrt(s2), n, mx, sxx };
}
// Student-t quantile (closed form for ν = 1, 2; Cornish–Fisher expansion in 1/ν otherwise)
export function tq(p, nu) { p = N.clamp(p, 1e-9, 1 - 1e-9); if (nu <= 1) return Math.tan(PI * (p - 0.5)); if (nu === 2) return (2 * p - 1) / Math.sqrt(2 * p * (1 - p));   // exact for 1 and 2 degrees of freedom, where the expansion is poor
  const z = N.normInv(p), v = Math.max(1, nu), z3 = z ** 3, z5 = z ** 5, z7 = z ** 7, z9 = z ** 9; return z + (z3 + z) / (4 * v) + (5 * z5 + 16 * z3 + 3 * z) / (96 * v * v) + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * v ** 3) + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * v ** 4); }
export function smooth(y, h = 2) { const n = y.length, o = new Array(n); for (let i = 0; i < n; i++) { let s = 0, c = 0; for (let j = Math.max(0, i - h); j <= Math.min(n - 1, i + h); j++) { s += y[j]; c++; } o[i] = s / c; } return o; }
// SNIP background: iterative peak clipping with a shrinking half-window of m channels
export function snip(y, m) { const n = y.length, b = [...y]; for (let p = Math.min(m, (n - 1) >> 1); p >= 1; p--) { const t = [...b]; for (let i = p; i < n - p; i++) { const a = 0.5 * (t[i - p] + t[i + p]); if (a < b[i]) b[i] = a; } } return b; }
function noiseOf(y) { const d = []; for (let i = 1; i < y.length - 1; i++) d.push(Math.abs(y[i] - 0.5 * (y[i - 1] + y[i + 1]))); return d.length ? 1.4826 * N.quantile(d, 0.5) / Math.sqrt(1.5) : 0; }
// min/max-preserving decimation to ≤ 400 points
export function thin(x, y, cap = 400) {
  const X = [], Y = []; for (let i = 0; i < x.length; i++) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) { X.push(x[i]); Y.push(y[i]); }
  if (X.length <= cap) return { x: X, y: Y };
  const w = Math.ceil(X.length / (cap / 2)), ox = [], oy = [];
  for (let i = 0; i < X.length; i += w) { let a = i, b = i; for (let j = i; j < Math.min(X.length, i + w); j++) { if (Y[j] < Y[a]) a = j; if (Y[j] > Y[b]) b = j; } for (const k of a === b ? [a] : [Math.min(a, b), Math.max(a, b)]) { ox.push(X[k]); oy.push(Y[k]); } }
  return { x: ox, y: oy };
}
export function makeOut() {
  const o = { summary: '', kpis: [], series: [], maps: [], tables: [], equations: [], checks: [], warnings: [], recommendations: [], publish: {}, exports: [], numerics: { method: '' }, synthetic: false };
  o.K = (key, label, value, unit, explain, good = null) => o.kpis.push({ key, label, value, unit, explain, good });
  o.S = (key, title, xlabel, ylabel, lines, explain, opt = {}) => o.series.push({ key, title, xlabel, ylabel, logx: !!opt.logx, logy: !!opt.logy, lines: lines.map(l => ({ name: l[0], ...thin(l[1], l[2]) })), explain });
  o.T = (key, title, columns, rows, explain) => o.tables.push({ key, title, columns, rows: rows.map(r => r.map(v => (typeof v === 'number' ? (Number.isFinite(v) ? +v.toPrecision(5) : '—') : String(v ?? '—')))), explain });
  o.E = (kind, label, text, note) => o.equations.push(note ? { kind, label, text, note } : { kind, label, text });
  o.C = (name, pass, detail) => o.checks.push({ name, pass: !!pass, detail: String(detail) });
  o.W = t => { if (!o.warnings.includes(t)) o.warnings.push(t); };
  o.R = (priority, text, basis, sustainability) => o.recommendations.push(sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
  o.tbl = (text, minCols, minRows, label) => { if (isSynthetic(text)) o.synthetic = true; return parseTable(text, minCols, minRows, label); };
  return o;
}
export function sealOut(o) {
  for (const k of o.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; k.good = null; k.explain += ' Not computable for this dataset (reported as 0) — see warnings.'; o.W(`"${k.label}" could not be computed from the supplied data and is reported as 0.`); }
  for (const k of Object.keys(o.publish)) if (typeof o.publish[k] === 'number' && !Number.isFinite(o.publish[k])) delete o.publish[k];
  if (o.synthetic) o.warnings.unshift(`The input table is ${SYN_MSG}.`);
  const { K, S, T, E, C, W, R, tbl, synthetic, ...r } = o; if (!r.maps.length) delete r.maps; return r;
}

// ============================================================================================
// X-ray diffraction kernel (exported: reused by the reverse-engineering suite)
// ============================================================================================
const pv = (x, x0, h, w, eta) => { const u = (x - x0) / w, u2 = 4 * u * u; return h * (eta / (1 + u2) + (1 - eta) * Math.exp(-LN2 * u2)); };
const pvArea = (h, w, eta) => h * w * (eta * PI / 2 + (1 - eta) * Math.sqrt(PI / LN2) / 2);
/** Find and fit peaks of a background-free profile with pseudo-Voigt functions (clusters fitted jointly, local linear baseline). */
export function fitPeaks(x, y, noise, { maxPeaks = 24, minRel = 0.012, minSep = 0.15 } = {}) {
  const n = x.length, step = (x[n - 1] - x[0]) / Math.max(1, n - 1) || 1, ys = smooth(y, 2), ymax = N.maxOf(ys);
  if (!(ymax > 0) || n < 8) return { peaks: [], fit: new Array(n).fill(0), iterations: 0 };
  let pk = N.findPeaks(ys, { minProminence: Math.max(4 * noise, minRel * ymax), minDistance: Math.max(2, Math.round(minSep / step)) });
  pk = pk.sort((a, b) => b.y - a.y).slice(0, maxPeaks).sort((a, b) => a.i - b.i).map(p => { const h = ys[p.i]; let l = p.i, r = p.i; while (l > 0 && ys[l] > h / 2) l--; while (r < n - 1 && ys[r] > h / 2) r++; return { pos: x[p.i], h, w: Math.max(2 * step, (r - l - 1) * step) }; });
  pk.forEach((p, i) => { const d = Math.min(i ? p.pos - pk[i - 1].pos : Infinity, i < pk.length - 1 ? pk[i + 1].pos - p.pos : Infinity); if (p.w > d) p.w = Math.max(2 * step, d * 0.6); });
  const out = [], fit = new Array(n).fill(0); let cl = [], iters = 0;
  const flush = () => {
    if (!cl.length) return; const k = cl.length, lo = cl[0].pos - 5 * cl[0].w, hi = cl[k - 1].pos + 5 * cl[k - 1].w;
    const i0 = Math.max(0, Math.floor((lo - x[0]) / step)), i1 = Math.min(n - 1, Math.ceil((hi - x[0]) / step)), xs = x.slice(i0, i1 + 1), yv = y.slice(i0, i1 + 1), xm = xs[0];
    let p = [], b = []; for (const q of cl) { p.push(q.pos, q.h, q.w, 0.5); b.push([q.pos - q.w, q.pos + q.w], [0, 4 * q.h + 1e-12], [Math.max(step * 0.5, q.w / 5), q.w * 5], [0, 1]); }
    p.push(0, 0); b.push([-Infinity, Infinity], [-Infinity, Infinity]); let se = null, ok = false;
    if (xs.length >= p.length + 2) {
      const model = (q, xx) => { let s = q[4 * k] + q[4 * k + 1] * (xx - xm); for (let j = 0; j < k; j++) s += pv(xx, q[4 * j], q[4 * j + 1], q[4 * j + 2], q[4 * j + 3]); return s; };
      try { const lm = N.levenbergMarquardt(q => xs.map((xx, i) => { const v = model(q, xx) - yv[i]; return Number.isFinite(v) ? v : 1e30; }), p, { maxIter: 25, tol: 1e-7, bounds: b }); if (lm.p.every(Number.isFinite)) { p = lm.p; se = lm.stderr; ok = true; iters += lm.iterations; } } catch { /* keep raw estimates */ }
    }
    for (let j = 0; j < k; j++) {
      const [pos, h, w, eta] = p.slice(4 * j, 4 * j + 4), s = (i, d) => (se && Number.isFinite(se[4 * j + i]) && se[4 * j + i] > 0 ? se[4 * j + i] : d);
      if (!(h > 2 * noise) || !(w > 0)) continue;
      out.push({ pos, sePos: s(0, step / 2), h, fwhm: w, seFwhm: s(2, step), eta, area: pvArea(h, w, eta), fitted: ok });
      for (let i = 0; i < n; i++) { const d = Math.abs(x[i] - pos); if (d < 30 * w) fit[i] += pv(x[i], pos, h, w, eta); }
    }
    cl = [];
  };
  for (const q of pk) { if (cl.length && (q.pos - cl[cl.length - 1].pos > 2.5 * (q.w + cl[cl.length - 1].w) || cl.length >= 4)) flush(); cl.push(q); }
  flush();
  return { peaks: out, fit, iterations: iters };
}
const REF = new Map();
/** Reference reflections (kinematic powder pattern) of a CRYSTALS entry for the given wavelength and angular range. */
export const refLimit = (lam, tmax) => Math.min(tmax, Math.floor(2 * Math.asin(Math.min(1, lam / 1.8)) / D2R * 10) / 10);   // reference reflections are generated for d ≥ 0.9 Å
export function refPeaks(key, lam, tmin, tmax) {
  tmax = refLimit(lam, tmax); if (!(tmax > tmin + 1)) return [];
  const k = `${key}|${lam.toFixed(4)}|${tmin.toFixed(1)}|${tmax.toFixed(1)}`; if (REF.has(k)) return REF.get(k);
  let pk = []; try { pk = powderXRD(key, { wavelength_A: lam, twoThetaMin: tmin, twoThetaMax: tmax, step: 0.1, crystallite_nm: 40 }).peaks.filter(p => p.I_rel >= 2).map(p => ({ tt: p.twoTheta, I: p.I_rel, hkl: p.hkl, d: p.d_A })); } catch { pk = []; }
  if (REF.size > 500) REF.clear(); REF.set(k, pk); return pk;
}
/** Figure-of-merit ranking of all library phases against observed peaks [{pos, area}]. */
export function rankPhases(obs, lam, tmin, tmax, sigPos = 0.08) {
  const tol = 3 * sigPos + 0.12, out = []; obs = obs.filter(p => p.pos <= refLimit(lam, tmax) + tol); const Io = N.sum(obs.map(p => p.area)) || 1;
  for (const key of Object.keys(CRYSTALS)) {
    if (key === 'graphene') continue;                                  // vacuum-padded cell: c is not a lattice spacing
    const ref = refPeaks(key, lam, tmin, tmax); if (!ref.length) continue;
    const Ir = N.sum(ref.map(p => p.I)), used = new Set(), pairs = [];
    for (const r of [...ref].sort((a, b) => b.I - a.I)) { let best = -1, bd = tol; obs.forEach((p, j) => { if (used.has(j)) return; const d = Math.abs(p.pos - r.tt); if (d < bd) { bd = d; best = j; } }); if (best >= 0) { used.add(best); pairs.push({ ref: r, obs: obs[best] }); } }
    const def = CRYSTALS[key], row = { key, name: def.name, formula: def.formula, system: def.system, nRef: ref.length, nMatch: pairs.length, pairs, fom: 0, R: 0, O: 0, scale: 0, meanAbs: 0 };
    if (pairs.length) {
      const wI = N.sum(pairs.map(q => q.ref.I)), s = N.sum(pairs.map(q => q.ref.I * (Math.sin(q.ref.tt / 2 * D2R) / Math.sin(q.obs.pos / 2 * D2R) - 1))) / wI;
      let z2 = 0, ma = 0; for (const q of pairs) { const tts = 2 * Math.asin(N.clamp(Math.sin(q.ref.tt / 2 * D2R) / (1 + s), 0, 1)) / D2R; z2 += ((q.obs.pos - tts) / sigPos) ** 2; ma += Math.abs(q.obs.pos - q.ref.tt); }
      z2 /= pairs.length; row.meanAbs = ma / pairs.length; row.scale = s;
      const so = N.sum(pairs.map(q => q.obs.area)); row.R = wI / Ir; row.O = so / Io;
      const rint = 1 - 0.5 * N.sum(pairs.map(q => Math.abs(q.obs.area / so - q.ref.I / wI)));
      row.rint = rint; row.fom = row.R * row.O * Math.exp(-0.5 * Math.min(z2, 60)) * Math.exp(-0.5 * (s / 0.003) ** 2) * (0.6 + 0.4 * rint);
    }
    out.push(row);
  }
  const tot = N.sum(out.map(r => r.fom)); for (const r of out) r.score = tot > 0 ? r.fom / tot : 0;
  return out.sort((a, b) => b.fom - a.fom);
}
const cellKind = sys => (/cubic/i.test(sys) ? 'cubic' : /hex|trig/i.test(sys) ? 'hex' : /tetra/i.test(sys) ? 'tet' : null);
/** Weighted linear least squares on Q = 1/d² for cubic / hexagonal / tetragonal cells. */
export function refineLattice(key, pairs, lam, sigPos = 0.08) {
  const def = CRYSTALS[key], kind = cellKind(def.system), a0 = def.a, c0 = def.c ?? def.a;
  const res = { kind, a: a0, c: c0, seA: 0, seC: 0, n: pairs.length, refined: false, chi2: 0, a0, c0 };
  if (!kind || !pairs.length) return res;
  const rows = pairs.map(q => { const [h, k, l] = q.ref.hkl, th = q.obs.pos / 2 * D2R, Q = 4 * Math.sin(th) ** 2 / (lam * lam), sQ = Math.max(1e-12, Q * 2 / Math.tan(th) * (Math.hypot(sigPos, q.obs.sePos || 0) / 2) * D2R); return { A: kind === 'cubic' ? [h * h + k * k + l * l] : kind === 'hex' ? [4 / 3 * (h * h + h * k + k * k), l * l] : [h * h + k * k, l * l], Q, sQ }; });
  const np = kind === 'cubic' ? 1 : 2, M = N.zeros(np, np), g = N.zeros(np);
  for (const r of rows) for (let i = 0; i < np; i++) { g[i] += r.A[i] * r.Q / r.sQ ** 2; for (let j = 0; j < np; j++) M[i][j] += r.A[i] * r.A[j] / r.sQ ** 2; }
  let p, cov;
  const ok = np === 1 ? M[0][0] > 0 : M[0][0] * M[1][1] - M[0][1] * M[1][0] > 1e-9 * M[0][0] * M[1][1];
  if (!ok) { // only one family of reflections: keep the library c/a and refine the scale
    if (!(M[0][0] + M[1][1] > 0)) return res;
    const ra = (c0 / a0) ** 2; let sw = 0, sq = 0; for (const r of rows) { const u = r.A[0] + r.A[1] / ra; sw += u * u / r.sQ ** 2; sq += u * r.Q / r.sQ ** 2; }
    if (!(sw > 0 && sq > 0)) return res; const A = sq / sw; res.a = 1 / Math.sqrt(A); res.c = res.a * c0 / a0; res.seA = res.a / 2 / Math.sqrt(sw) / A; res.seC = res.seA * c0 / a0; res.refined = true; res.constrained = true; return res;
  }
  try { p = N.solveLinear(M, g); cov = N.inverse(M); } catch { return res; }
  if (!p.every(v => v > 0)) return res;
  let chi2 = 0; for (const r of rows) chi2 += ((r.Q - N.dot(r.A, p)) / r.sQ) ** 2; const dof = Math.max(1, rows.length - np), sc = Math.max(1, chi2 / dof);
  res.a = 1 / Math.sqrt(p[0]); res.seA = res.a / 2 * Math.sqrt(Math.abs(cov[0][0]) * sc) / p[0];
  if (np === 2) { res.c = 1 / Math.sqrt(p[1]); res.seC = res.c / 2 * Math.sqrt(Math.abs(cov[1][1]) * sc) / p[1]; } else { res.c = res.a; res.seC = res.seA; }
  res.chi2 = chi2 / dof; res.refined = true; return res;
}
export function cellVolume(kind, a, c) { return kind === 'hex' ? Math.sqrt(3) / 2 * a * a * c : kind === 'tet' ? a * a * c : a ** 3; }
export function xrayDensity(key, a, c) { const def = CRYSTALS[key], kind = cellKind(def.system); if (!kind) return NaN; let M = 0; for (const s of def.basis) M += (chem.ELEMENTS[s.el]?.mass || 0) * (s.occ ?? 1); return M * 1e-3 / (N.N_A * cellVolume(kind, a, c) * 1e-30); }
/** Complete reduction of a powder pattern: background, peaks, Scherrer, Williamson–Hall, phase ranking, lattice refinement, texture. */
export function xrdAnalyse(X, Y, { lam = 1.5406, inst = 0.08, K = 0.9, sigPos = 0.08, corr = 'quadrature' } = {}) {
  let pts = X.map((v, i) => [v, Y[i]]).filter(p => p[0] > 0.5 && p[0] < 178 && Number.isFinite(p[1])).sort((a, b) => a[0] - b[0]);
  if (pts.length > 6000) { const g = Math.ceil(pts.length / 6000), q = []; for (let i = 0; i + g <= pts.length; i += g) { let sx = 0, sy = 0; for (let j = i; j < i + g; j++) { sx += pts[j][0]; sy += pts[j][1]; } q.push([sx / g, sy / g]); } pts = q; }
  const x = pts.map(p => p[0]), y = pts.map(p => p[1]), n = x.length, step = (x[n - 1] - x[0]) / Math.max(1, n - 1) || 0.02;
  const bg = snip(smooth(y, 2), N.clamp(Math.round(2.5 / step), 6, 300)), net = y.map((v, i) => v - bg[i]), noise = noiseOf(y);
  const fp = fitPeaks(x, net, noise, { minSep: 0.15 }), peaks = fp.peaks;
  for (const p of peaks) {
    const th = p.pos / 2 * D2R; p.d = lam / (2 * Math.sin(th));
    const b = corr === 'linear' ? p.fwhm - inst : Math.sqrt(Math.max(0, p.fwhm ** 2 - inst ** 2));
    p.beta = b * D2R; p.resolved = b > 0.15 * Math.max(inst, step); p.size = p.resolved ? K * lam / (p.beta * Math.cos(th)) / 10 : NaN;      // nm
    p.seSize = p.resolved ? p.size * p.seFwhm * p.fwhm / Math.max(1e-12, b * (corr === 'linear' ? p.fwhm : b)) : NaN;
  }
  const res = peaks.filter(p => p.resolved), top = [...res].sort((a, b) => b.area - a.area).slice(0, 5);
  const size = top.length ? N.mean(top.map(p => p.size)) : NaN, sizeSd = top.length > 1 ? N.std(top.map(p => p.size)) : (top[0]?.seSize ?? NaN);
  let wh = null;
  if (res.length >= 3) { const xs = res.map(p => 4 * Math.sin(p.pos / 2 * D2R)), yv = res.map(p => p.beta * Math.cos(p.pos / 2 * D2R)), f = linfit(xs, yv); wh = { xs, yv, ...f, strain: f.slope, seStrain: f.seSlope, size: f.intercept > 0 ? K * lam / f.intercept / 10 : NaN, seSize: f.intercept > 0 ? K * lam / f.intercept / 10 * f.seInt / f.intercept : NaN }; }
  const tmin = Math.floor(x[0] * 10) / 10, tmax = Math.ceil(x[n - 1] * 10) / 10, phases = peaks.length ? rankPhases(peaks, lam, tmin, tmax, sigPos) : [], best = phases[0] && phases[0].fom > 0 ? phases[0] : null;
  let lattice = null, texture = null, density = NaN;
  if (best) {
    lattice = refineLattice(best.key, best.pairs, lam, sigPos); density = xrayDensity(best.key, lattice.a, lattice.c);
    for (const q of best.pairs) q.obs.hkl = q.ref.hkl;
    if (best.pairs.length >= 3) { const so = N.sum(best.pairs.map(q => q.obs.area)), sr = N.sum(best.pairs.map(q => q.ref.I)), r = best.pairs.map(q => (q.obs.area / so) / (q.ref.I / sr)), m = N.mean(r), tc = r.map(v => v / m), k = N.argmax(tc); texture = { tcMax: tc[k], hkl: best.pairs[k].ref.hkl, tc }; }
  }
  let ss = 0, st = 0; for (let i = 0; i < n; i++) { ss += (net[i] - fp.fit[i]) ** 2; st += net[i] ** 2; }
  return { x, y, bg, net, fit: fp.fit, noise, step, peaks, size, sizeSd, wh, phases, best, lattice, texture, density, rwp: st > 0 ? Math.sqrt(ss / st) : 0, iterations: fp.iterations, lam, tmin, tmax };
}
export function synthXRD(key, { size = 25, strain = 5e-4, seed = 11, peak = 6000, tmin = 20, tmax = 80, step = 0.05 } = {}) {
  const p = powderXRD(key, { crystallite_nm: size, strain, twoThetaMin: tmin, twoThetaMax: tmax, step }), r = N.rng(seed);
  return syn(csv(['two_theta_deg', 'intensity'], p.twoTheta.map((t, i) => { const v = peak * p.intensity[i] / 100 + 180 + 120 * Math.exp(-(t - tmin) / 25); return [t, Math.max(0, Math.round(v + r.normal(0, Math.sqrt(v))))]; })));
}
const hklS = h => (h ? `(${h.join(h.some(v => Math.abs(v) > 9) ? ' ' : '')})` : '—');

// ============================================================================================
// Seeded synthetic demonstration datasets (forward model + noise)
// ============================================================================================
const XRD_ZNO = synthXRD('ZnO_wurtzite', { size: 25, seed: 11 }), XRD_AU = synthXRD('Au', { size: 12, seed: 12, strain: 1e-3 });
const DLS_Q = (n, lamNm, angDeg) => 4 * PI * n / (lamNm * 1e-9) * Math.sin(angDeg / 2 * D2R);
function synthDLS(d_nm, pdi, seed) { const r = N.rng(seed), D = N.K_B * 298.15 / (3 * PI * 0.89e-3 * d_nm * 1e-9), G = D * DLS_Q(1.33, 633, 173) ** 2, mu2 = pdi * G * G; return syn(csv(['tau_s', 'g2'], N.logspace(1e-6, 1e-1, 64).map(t => { const g1 = Math.exp(-G * t) * (1 + mu2 * t * t / 2); return [t, 1 + 0.85 * g1 * g1 + r.normal(0, 0.003)]; }))); }
const DLS_DEF = synthDLS(60, 0.08, 3), DLS_AU = synthDLS(24, 0.12, 4);
const synthTEM = (med, gsd, n, seed) => { const r = N.rng(seed); return syn(csv(['diameter_nm'], Array.from({ length: n }, () => [r.lognormal(Math.log(med), Math.log(gsd))]))); };
const TEM_DEF = synthTEM(42, 1.22, 200, 5), TEM_AU = synthTEM(13, 1.12, 150, 6);
const N2 = { rK: 2 * 8.85e-3 * 34.7e-6 / (N.R_GAS * 77.35) * 1e9, liq: 1.547e-3, t: x => 0.354 * Math.cbrt(5 / -Math.log(x)), tHJ: x => 0.1 * Math.sqrt(13.99 / (0.034 - Math.log10(x))) };   // Kelvin constant (nm), cm³ liquid per cm³ STP, Halsey film thickness (nm)
const BET_DEF = (() => { const r = N.rng(7), vm = 23, Cc = 110, nL = 5, xs = [0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99];
  return syn(csv(['p_rel', 'v_ads_cm3_g'], xs.map(x => { const v = vm * Cc * x / (1 - x) * (1 - (nL + 1) * x ** nL + nL * x ** (nL + 1)) / (1 + (Cc - 1) * x - Cc * x ** (nL + 1)), rp = N2.rK / -Math.log(x) + N2.t(x); return [x, (v + 120 / (1 + Math.exp(-(Math.log(rp) - Math.log(4.5)) / 0.12))) * (1 + r.normal(0, 0.004))]; }))); })();
const UV_DEF = (() => { const r = N.rng(8), Eg = 3.30; return syn(csv(['wavelength_nm', 'absorbance'], N.linspace(300, 520, 111).map(l => { const E = 1239.842 / l; return [l, 0.9 * Math.sqrt(Math.max(0, E - Eg)) / E * 3.3 + 0.05 * Math.exp(Math.min(0, E - Eg) / 0.06) + 0.03 + r.normal(0, 0.003)]; }))); })();
const synthLSPR = (peakNm, ratio, seed) => { const r = N.rng(seed), base = l => 0.5 * (450 / l) ** 3.2, w = 42; return syn(csv(['wavelength_nm', 'absorbance'], N.linspace(380, 800, 141).map(l => [l, base(l) + (ratio * 0.5 - base(peakNm)) / (1 - ratio / (1 + ((450 - peakNm) / w) ** 2)) / (1 + ((l - peakNm) / w) ** 2) + r.normal(0, 0.002)]))); };
const LSPR_DEF = synthLSPR(524, (Math.log(20) + 2.2) / 3, 9);
const synthRaman = (modes, seed, lo, hi) => { const r = N.rng(seed); return syn(csv(['shift_cm1', 'intensity'], N.linspace(lo, hi, Math.round((hi - lo) / 2) + 1).map(s => { let v = 60 + 0.01 * s; for (const [p, h, w] of modes) v += h / (1 + ((s - p) / (w / 2)) ** 2); return [s, v + r.normal(0, 4)]; }))); };
const RAMAN_DEF = synthRaman([[520.2, 900, 5], [1350, 180, 40], [1582, 600, 22], [2690, 950, 34]], 10, 100, 3000);
const RAMAN_ZNO = synthRaman([[437.6, 800, 9], [574, 120, 24], [332, 90, 20]], 13, 200, 1200);
const DIAMOND_TIP = { E: 1141e9, nu: 0.07 };
function synthOP(E, nu, H, Pmax_mN, seed) {
  const r = N.rng(seed), Er = 1 / ((1 - nu * nu) / E + (1 - DIAMOND_TIP.nu ** 2) / DIAMOND_TIP.E), P = Pmax_mN * 1e-3, A = P / H, hc = Math.sqrt(A / 24.5), S = 2 * 1.034 / Math.sqrt(PI) * Er * Math.sqrt(A), hmax = hc + 0.75 * P / S, m = 1.3, hf = hmax - m * P / S, rows = [];
  for (let i = 1; i <= 40; i++) { const h = hmax * i / 40; rows.push([h * 1e9 + r.normal(0, 0.15), Pmax_mN * (h / hmax) ** 2 * (1 + r.normal(0, 0.002))]); }
  for (let i = 1; i <= 45; i++) { const h = hmax - (hmax - hf) * i / 46; rows.push([h * 1e9 + r.normal(0, 0.15), Pmax_mN * ((h - hf) / (hmax - hf)) ** m * (1 + r.normal(0, 0.002))]); }
  return syn(csv(['depth_nm', 'load_mN'], rows));
}
const OP_DEF = synthOP(72e9, 0.17, 9.3e9, 10, 14);
const HZ_DEF = (() => { const r = N.rng(15), Er = 2e6 / (1 - 0.25), Rt = 2.5e-6; return syn(csv(['indentation_nm', 'force_nN'], N.linspace(-40, 240, 57).map(d => [d, 4 / 3 * Er * Math.sqrt(Rt) * Math.max(0, d * 1e-9) ** 1.5 * 1e9 + r.normal(0, 3)]))); })();
const FPP_DEF = (() => { const r = N.rng(16), Rm = 2e-6 / 150e-9 * LN2 / PI; return syn(csv(['current_mA', 'voltage_mV'], N.linspace(-1, 1, 21).map(i => [i, i * Rm + r.normal(0, 0.01)]))); })();
const HALL_DEF = (() => { const r = N.rng(17); return syn(csv(['B_T', 'V_hall_mV'], N.linspace(-1, 1, 11).map(b => [b, -1e-3 * b / (N.Q_E * 3e26 * 150e-9) * 1e3 + 0.02 + r.normal(0, 0.002)]))); })();
function randles(f, Rs, Rct, Q, al, sw) { const w = 2 * PI * f, yr = Q * w ** al * Math.cos(al * PI / 2), yi = Q * w ** al * Math.sin(al * PI / 2), zr = Rct + sw / Math.sqrt(w), zi = -sw / Math.sqrt(w), d = zr * zr + zi * zi, Yr = yr + zr / d, Yi = yi - zi / d, D = Yr * Yr + Yi * Yi; return [Rs + Yr / D, -Yi / D]; }
const EIS_DEF = (() => { const r = N.rng(18); return syn(csv(['f_Hz', 'Z_re_ohm', 'Z_im_ohm'], N.logspace(1e5, 1e-2, 50).map(f => { const z = randles(f, 20, 150, 3e-5, 0.9, 30); return [f, z[0] * (1 + r.normal(0, 0.006)), z[1] * (1 + r.normal(0, 0.006))]; }))); })();
const parker = w => { let s = 1; for (let k = 1; k < 60; k++) s += 2 * (k % 2 ? -1 : 1) * Math.exp(-k * k * w); return w < 0.05 ? 0 : s; };   // rear-face rise, w = π²αt/L²
const parkerK = fr => N.bisect(w => parker(w) - fr, 0.05, 20, 1e-10) / (PI * PI);
const FLASH_DEF = (() => { const r = N.rng(19), al = 1.8e-5, L = 1e-3; return syn(csv(['time_ms', 'signal'], N.linspace(0, 60, 121).map(t => [t, parker(PI * PI * al * t * 1e-3 / (L * L)) + r.normal(0, 0.008)]))); })();
function synthTens(E_MPa, sy, n, smax, nu, seed, noise) { const r = N.rng(seed); return syn(csv(['strain', 'stress_MPa', 'transverse_strain'], N.linspace(0, smax, 70).map(s => { const e = s / E_MPa + 0.002 * (s / sy) ** n; return [e, s + r.normal(0, noise), -nu * e + r.normal(0, 2e-6 + 2e-3 * Math.abs(nu * e))]; }))); }
const TENS_DEF = synthTens(69000, 276, 12, 318, 0.33, 20, 0.4), TENS_AUX = synthTens(45, 1.2, 8, 1.45, -0.55, 21, 0.003);
const VSM_DEF = (() => { const r = N.rng(22), Ms = 350, Hc = 8, k = Math.tan(PI * 0.3 / 2) / Hc, Hs = [...N.linspace(800, -800, 81), ...N.linspace(-800, 800, 81)]; return syn(csv(['H_kA_m', 'M_kA_m'], Hs.map((H, i) => [H, 2 * Ms / PI * Math.atan((i < 81 ? H + Hc : H - Hc) * k) + r.normal(0, 1.2)]))); })();
const REP_DEF = (() => { const r = N.rng(23); return syn(csv(['value'], Array.from({ length: 10 }, () => [r.normal(100.9, 1.1)]))); })();

// GUM reductions: value = f(named inputs). Defaults are typical magnitudes; the budget table overrides them by name.
const GUM = {
  scherrer: { label: 'Scherrer crystallite size D = Kλ/(β cosθ)', unit: 'nm', eq: 'D = K·λ / (√(B² − b²)·cos θ)', vars: [['K', 0.9, 0.05], ['lambda_A', 1.5406, 0.0001], ['fwhm_deg', 0.345, 0.008], ['inst_deg', 0.08, 0.01], ['two_theta_deg', 36.25, 0.02]], f: v => v[0] * v[1] / (Math.sqrt(Math.max(1e-12, v[2] ** 2 - v[3] ** 2)) * D2R * Math.cos(v[4] / 2 * D2R)) / 10 },
  bragg: { label: 'Bragg d-spacing d = λ/(2 sinθ)', unit: 'Å', eq: 'd = λ / (2·sin θ)', vars: [['lambda_A', 1.5406, 0.0001], ['two_theta_deg', 36.25, 0.02]], f: v => v[0] / (2 * Math.sin(Math.max(1e-6, v[1] / 2 * D2R))) },
  stokes: { label: 'Stokes–Einstein diameter d = k_B·T/(3πηD)', unit: 'nm', eq: 'd_h = k_B·T / (3π·η·D)', vars: [['T_K', 298.15, 0.2], ['eta_mPas', 0.89, 0.02], ['D_m2_s', 8.2e-12, 2e-13]], f: v => N.K_B * v[0] / (3 * PI * v[1] * 1e-3 * Math.max(1e-30, v[2])) * 1e9 },
  bet_size: { label: 'BET equivalent diameter d = 6/(ρ·SSA)', unit: 'nm', eq: 'd = 6 / (ρ·S_BET)', vars: [['rho_kg_m3', 5606, 30], ['ssa_m2_g', 100, 3]], f: v => 6 / Math.max(1e-30, v[0] * v[1] * 1e3) * 1e9 },
  fpp: { label: 'Four-point-probe resistivity ρ = (π/ln 2)·(V/I)·t', unit: 'Ω·m', eq: 'ρ = (π/ln 2)·(V/I)·t', vars: [['V_mV', 2.94, 0.01], ['I_mA', 1, 0.002], ['t_nm', 150, 5]], f: v => PI / LN2 * v[0] / Math.max(1e-30, v[1]) * v[2] * 1e-9 },
  flash: { label: 'Laser-flash conductivity k = 0.1388·L²/t½·ρ·c_p', unit: 'W/m/K', eq: 'k = 0.1388·L²/t½ · ρ·c_p', vars: [['L_mm', 1, 0.005], ['t_half_ms', 7.7, 0.15], ['rho_kg_m3', 5606, 30], ['cp_J_kgK', 494, 15]], f: v => 0.1388 * (v[0] * 1e-3) ** 2 / Math.max(1e-30, v[1] * 1e-3) * v[2] * v[3] },
  hardness: { label: 'Indentation hardness H = P/(24.5·h_c²)', unit: 'GPa', eq: 'H = P_max / (24.5·h_c²)', vars: [['P_mN', 10, 0.02], ['hc_nm', 209.5, 1.5]], f: v => v[0] * 1e-3 / (24.5 * Math.max(1e-30, (v[1] * 1e-9) ** 2)) / 1e9 },
};
const gumTable = k => csv(['quantity', 'value', 'u_standard'], GUM[k].vars);
const GUM_DEF = syn(gumTable('scherrer'));

// Well-known vibrational mode positions (cm⁻¹) used for peak assignment
const RAMAN_MODES = [['ZnO (wurtzite)', 'E2(high)', 437], ['Si (crystalline)', 'F2g optical phonon', 520.7], ['Graphene / graphitic carbon', 'D (defect-activated)', 1350], ['Graphene / graphitic carbon', 'G (E2g)', 1580], ['Graphene / graphitic carbon', '2D (G′)', 2700],
  ['TiO2 anatase', 'Eg', 144], ['TiO2 anatase', 'B1g', 399], ['TiO2 anatase', 'A1g + B1g', 516], ['TiO2 anatase', 'Eg', 639], ['TiO2 rutile', 'Eg', 447], ['TiO2 rutile', 'A1g', 612], ['Diamond', 'F2g', 1332]];
// Indicative LSPR peak positions of citrate-type silver spheres in water (nm diameter → nm peak); order-of-magnitude guide only
const AG_LSPR = [[10, 395], [20, 400], [40, 412], [60, 430], [80, 455], [100, 490]];

// ============================================================================================
// Analyses
// ============================================================================================
function mXRD(inp, o) {
  const t = o.tbl(inp.xrd_table, 2, 30, 'XRD pattern'), lam = inp.xrd_lambda, A = xrdAnalyse(t.rows.map(r => r[0]), t.rows.map(r => r[1]), { lam, inst: inp.xrd_inst, K: inp.xrd_K, sigPos: inp.xrd_pos_tol, corr: inp.xrd_corr });
  const b = A.best, L = A.lattice, second = A.phases[1], strongest = A.peaks.length ? A.peaks[N.argmax(A.peaks.map(p => p.area))] : null;
  o.K('n_peaks', 'Fitted diffraction peaks', A.peaks.length, '', `Number of peaks above 4σ noise (σ = ${F(A.noise, 3)} counts) after SNIP background removal, each fitted with a pseudo-Voigt profile by Levenberg–Marquardt.`, A.peaks.length >= 3);
  o.K('phase', 'Best-matching phase', b ? `${b.name} [${b.formula}]` : 'no match', '', b ? `Highest figure of merit among ${A.phases.length} library structures: ${b.nMatch} of ${b.nRef} reference reflections matched, mean |Δ2θ| = ${F(b.meanAbs, 3)}°. It is a ranked candidate, not a proof of identity.` : 'No library structure reproduces the observed peak positions at this wavelength; check the wavelength and the 2θ calibration.', b ? b.score > 0.5 : false);
  o.K('score', 'Evidence score of best phase', b ? b.score : 0, '', `Figure of merit normalised over all library candidates (sums to 1). Runner-up: ${second && second.fom > 0 ? `${second.name} with ${F(second.score, 3)}` : 'none'}. Values below about 0.6 mean the pattern alone does not discriminate.`, b ? b.score > 0.6 : false);
  o.K('a', 'Refined lattice parameter a', L && L.refined ? L.a : 0, 'Å', L && L.refined ? `Weighted least squares on 1/d² of ${L.n} indexed reflections (${L.kind} cell${L.constrained ? ', c/a held at the library value' : ''}); standard error ${F(L.seA, 2)} Å, library value ${F(L.a0, 5)} Å.` : 'Not refined: the best phase is not cubic, hexagonal or tetragonal, or no reflection could be indexed.', L && L.refined ? Math.abs(L.a / L.a0 - 1) < 0.005 : null);
  o.K('c', 'Refined lattice parameter c', L && L.refined ? L.c : 0, 'Å', L && L.refined ? `Second cell edge from the same refinement (equal to a for cubic cells); standard error ${F(L.seC, 2)} Å, library value ${F(L.c0, 5)} Å. No zero-shift or sample-displacement term is refined.` : 'Not refined for this dataset — see the lattice parameter a entry for the reason.');
  o.K('size', 'Scherrer crystallite size', fin(A.size), 'nm', `Mean of the ${Math.min(5, A.peaks.filter(p => p.resolved).length)} strongest resolved peaks, D = Kλ/(β cos θ) with K = ${inp.xrd_K} and instrument breadth ${inp.xrd_inst}° removed (${inp.xrd_corr}); scatter between peaks ±${F(fin(A.sizeSd), 2)} nm. It is a volume-weighted column length, not a particle diameter.`, Number.isFinite(A.size) ? A.size < 120 : false);
  o.K('wh_size', 'Williamson–Hall size', A.wh ? fin(A.wh.size) : 0, 'nm', A.wh ? `From the intercept of β cos θ against 4 sin θ over ${A.wh.n} peaks (R² = ${F(A.wh.r2, 3)}); ± ${F(fin(A.wh.seSize), 2)} nm. A non-positive intercept is reported as 0 and means size broadening is not separable.` : 'At least three resolved peaks are needed for a Williamson–Hall line.');
  o.K('strain', 'Microstrain (Williamson–Hall slope)', A.wh ? A.wh.strain : 0, '', A.wh ? `Slope of the Williamson–Hall line, ε ± ${F(A.wh.seStrain, 2)}; a slope within two standard errors of zero means strain broadening is not resolved from size broadening.` : 'Not available with fewer than three resolved peaks; strain and size cannot be separated from one peak.', A.wh ? Math.abs(A.wh.strain) < 5e-3 : null);
  o.K('d_main', 'd-spacing of strongest peak', strongest ? strongest.d : 0, 'Å', strongest ? `Bragg's law for the peak at 2θ = ${F(strongest.pos, 5)}° (± ${F(strongest.sePos, 2)}° fit error) with λ = ${lam} Å${strongest.hkl ? `; indexed as ${hklS(strongest.hkl)}` : ''}.` : 'No peak was found above the noise threshold in this pattern.');
  o.K('texture', 'Largest texture coefficient', A.texture ? A.texture.tcMax : 0, '', A.texture ? `Harris texture coefficient TC = (I/I₀)/mean(I/I₀) against the simulated random-powder intensities; maximum for ${hklS(A.texture.hkl)}. TC ≫ 1 hints at preferred orientation (e.g. c-axis aligned rods) but also absorbs errors of the kinematic reference.` : 'Fewer than three matched reflections — texture cannot be assessed.', A.texture ? A.texture.tcMax < 1.6 : null);
  o.K('rho_x', 'X-ray density', fin(A.density), 'kg/m³', 'Mass of the unit-cell contents divided by the refined cell volume; an upper bound of the bulk density (no pores, vacancies or amorphous fraction).');
  o.K('rwp', 'Profile residual', A.rwp, '', 'Root of the summed squared difference between the background-free pattern and the sum of fitted peaks, relative to the pattern norm; values above about 0.15 mean missed or badly shaped peaks.', A.rwp < 0.15);
  o.S('pattern', 'Diffraction pattern, background and fit', '2θ (°)', 'intensity (counts)', [['measured', A.x, A.y], ['SNIP background', A.x, A.bg], ['background + fitted peaks', A.x, A.fit.map((v, i) => v + A.bg[i])]], 'Measured pattern with the clipped background and the sum of fitted pseudo-Voigt peaks; unexplained intensity indicates a second phase or a poor profile shape.');
  if (A.wh) o.S('wh', 'Williamson–Hall plot', '4 sin θ', 'β cos θ (rad)', [['peaks', A.wh.xs, A.wh.yv], ['fit', [0, N.maxOf(A.wh.xs)], [A.wh.intercept, A.wh.intercept + A.wh.slope * N.maxOf(A.wh.xs)]]], 'Intercept gives Kλ/D (size), slope gives microstrain. Scatter about the line reflects anisotropic broadening and fit errors.');
  if (b) o.S('sticks', 'Observed peak areas against the reference pattern', '2θ (°)', 'relative intensity (%)', [['observed (normalised)', A.peaks.map(p => p.pos), A.peaks.map(p => 100 * p.area / N.maxOf(A.peaks.map(q => q.area)))], [`${b.formula} reference`, refPeaks(b.key, lam, A.tmin, A.tmax).map(p => p.tt), refPeaks(b.key, lam, A.tmin, A.tmax).map(p => p.I)]], 'Integrated peak areas compared with the simulated random-powder intensities of the top candidate; systematic ratios point to texture.');
  o.T('peaks', 'Fitted peaks', ['2θ (°)', '± (°)', 'd (Å)', 'FWHM (°)', 'η (Lorentz fraction)', 'area', 'hkl', 'size (nm)'], A.peaks.map(p => [p.pos, p.sePos, p.d, p.fwhm, p.eta, p.area, hklS(p.hkl), p.resolved ? p.size : 'beyond resolution']), 'Pseudo-Voigt parameters of each peak with the Bragg d-spacing, the index assigned from the best phase and the single-peak Scherrer size.');
  o.T('cands', 'Ranked candidate phases', ['rank', 'phase', 'formula', 'system', 'evidence score', 'FOM', 'matched/expected', 'mean |Δ2θ| (°)', 'cell scale (%)'], A.phases.slice(0, 8).map((p, i) => [i + 1, p.name, p.formula, p.system, p.score, p.fom, `${p.nMatch}/${p.nRef}`, p.meanAbs, 100 * p.scale]), 'FOM = (matched reference intensity fraction) × (explained observed intensity fraction) × position likelihood × cell-scale prior × intensity agreement. Close scores mean the phases cannot be told apart by this pattern.');
  o.E('governing', "Bragg's law", 'n·λ = 2·d·sin θ'); o.E('formula', 'Pseudo-Voigt peak profile', 'I(2θ) = h·[η/(1 + 4u²) + (1 − η)·exp(−4 ln2·u²)],  u = (2θ − 2θ₀)/FWHM');
  o.E('formula', 'Scherrer equation with instrument correction', inp.xrd_corr === 'linear' ? 'D = K·λ/(β·cos θ),  β = B_obs − b_inst' : 'D = K·λ/(β·cos θ),  β = √(B_obs² − b_inst²)'); o.E('formula', 'Williamson–Hall relation', 'β·cos θ = K·λ/D + 4·ε·sin θ');
  o.E('governing', 'Lattice refinement (weighted least squares)', 'min Σ wᵢ·(1/dᵢ² − Qᵢ(a, c))²;  cubic Q = (h²+k²+l²)/a²;  hexagonal Q = 4(h²+hk+k²)/(3a²) + l²/c²;  tetragonal Q = (h²+k²)/a² + l²/c²');
  o.E('formula', 'Phase-identification figure of merit', 'FOM = R_ref·R_obs·exp(−½⟨(Δ2θ/σ)²⟩)·exp(−½(s/0.3 %)²)·(0.6 + 0.4·R_int);  score = FOM/ΣFOM');
  if (b) { let M = NaN; try { M = chem.molarMass(b.formula); } catch { /* formula not parsable */ } o.E('formula', 'Identified compound', `${chem.prettyFormula(b.formula)} — ${b.name}, ${b.system}${Number.isFinite(M) ? `, M = ${F(M, 5)} g/mol` : ''}`); }
  o.C('Peaks found above noise', A.peaks.length >= 3, `${A.peaks.length} peaks; at least three are needed for phase identification and Williamson–Hall analysis`);
  o.C('Phase identification is discriminating', !!b && b.score > 0.6, b ? `top score ${F(b.score, 3)}${second && second.fom > 0 ? `, runner-up ${second.name} ${F(second.score, 3)}` : ''}` : 'no candidate matched');
  o.C('Peaks broader than the instrument', A.peaks.some(p => p.resolved), `instrument breadth ${inp.xrd_inst}°, narrowest fitted FWHM ${A.peaks.length ? F(N.minOf(A.peaks.map(p => p.fwhm)), 3) : '—'}°`);
  if (L && L.refined) o.C('Refined cell close to the library cell', Math.abs(L.a / L.a0 - 1) < 0.01, `Δa/a = ${F(100 * (L.a / L.a0 - 1), 3)} %, Δc/c = ${F(100 * (L.c / L.c0 - 1), 3)} %`);
  if (!b) o.W('No library phase matches: the wavelength may be wrong for this pattern, or the phase is not in the built-in structure library.');
  if (b && second && second.fom > 0.5 * b.fom) o.W(`${b.name} and ${second.name} fit almost equally well (isostructural cells within the position tolerance); composition analysis (EDS/XPS) is needed to decide.`);
  if (Number.isFinite(A.size) && A.size > 100) o.W('Crystallite size above about 100 nm: broadening is close to the instrument breadth and the Scherrer value is only a lower bound.');
  if (A.wh && A.wh.intercept <= 0) o.W('Williamson–Hall intercept is not positive: size broadening cannot be separated from strain with these peaks.');
  o.W('Kinematic single-wavelength reference patterns (no Kα₂, absorption or preferred orientation); no zero-shift refinement. Use a full Rietveld refinement for quantitative phase fractions.');
  o.R('high', b && b.score > 0.6 ? `Confirm ${b.formula} with a second method (Raman modes or EDS composition) before using a = ${F(L?.a ?? 0, 5)} Å downstream.` : 'Collect composition data (EDS/XPS) and a wider 2θ range: the diffraction pattern alone does not single out one phase.', 'Phase ranking score and runner-up margin.');
  o.R('medium', `Measure a line-profile standard (e.g. LaB₆ or Si powder) to replace the assumed instrument breadth of ${inp.xrd_inst}°; at ${F(fin(A.size), 3)} nm the size error scales with the breadth correction.`, 'Scherrer correction sensitivity.');
  if (A.texture && A.texture.tcMax > 1.5) o.R('medium', `Texture coefficient ${F(A.texture.tcMax, 3)} for ${hklS(A.texture.hkl)}: verify orientation with a rocking curve or pole figure, or grind/spin the sample if a random powder is intended.`, 'Harris texture coefficient above 1.5.');
  Object.assign(o.publish, b ? { 'material.formula': b.formula, 'material.name': b.name } : {}, L && L.refined ? { 'prop.lattice_a': L.a * 1e-10, 'prop.lattice_c': L.c * 1e-10, 'cred.u_lattice_a': L.seA * 1e-10 } : {}, Number.isFinite(A.size) ? { 'geom.crystallite': A.size * 1e-9, 'cred.u_rel': fin(A.sizeSd / A.size, 0) } : {}, { 'cred.rmse': A.rwp }, Number.isFinite(A.density) ? { 'prop.density_xray': A.density } : {});
  o.numerics = { method: 'SNIP background + pseudo-Voigt cluster fits (Levenberg–Marquardt) + weighted linear least-squares cell refinement', dof: A.x.length, iterations: A.iterations, residual: A.rwp, notes: `${A.phases.length} library phases screened; pattern step ${F(A.step, 3)}°` };
  o.summary = b ? `The pattern is best explained by ${b.name} (evidence score ${F(b.score, 3)}, ${b.nMatch}/${b.nRef} reflections) with a = ${F(L.a, 5)} Å${L.kind !== 'cubic' ? ` and c = ${F(L.c, 5)} Å` : ''}. Scherrer crystallite size is ${F(fin(A.size), 3)} nm${A.wh ? ` and the Williamson–Hall microstrain ${F(A.wh.strain, 2)}` : ''}.` : `${A.peaks.length} peaks were fitted but no library phase matches them at λ = ${lam} Å; no lattice parameters were refined.`;
}

function mDLS(inp, o) {
  const t = o.tbl(inp.dls_table, 2, 8, 'DLS autocorrelation'), rows = t.rows.filter(r => r[0] > 0).sort((a, b) => a[0] - b[0]), tau = rows.map(r => r[0]), g2 = rows.map(r => r[1]), n = tau.length;
  if (n < 8) throw new Error('DLS autocorrelation: at least 8 rows with positive delay time are needed.');
  const T = inp.dls_T + 273.15, eta = inp.dls_visc * 1e-3, q = DLS_Q(inp.dls_n, inp.dls_lambda, inp.dls_angle), nb = Math.max(3, Math.round(n * 0.1)), B = N.mean(g2.slice(n - nb)), beta0 = N.mean(g2.slice(0, 3)) - B;
  let m = 0; while (m < n && g2[m] - B > 0.1 * beta0) m++;
  const valid = beta0 > 0 && m >= 5; let G = NaN, mu2 = 0, beta = beta0, seG = NaN, its = 0;
  if (valid) {
    const ts = tau[m - 1], xs = tau.slice(0, m).map(v => v / ts), c = N.polyfit(xs, g2.slice(0, m).map(v => Math.log(v - B)), 2); G = -c[1] / 2 / ts; mu2 = c[2] / ts / ts; beta = Math.exp(c[0]);
    if (G > 0) { try { const lm = N.levenbergMarquardt(p => xs.map((x, i) => B + p[0] * Math.exp(-2 * p[1] * x) * (1 + p[2] * x * x / 2) ** 2 - g2[i]), [beta, G * ts, Math.max(0, mu2) * ts * ts], { maxIter: 60, bounds: [[0, 2], [1e-9, 1e9], [0, 1e9]] }); if (lm.p[1] > 0) { beta = lm.p[0]; G = lm.p[1] / ts; mu2 = lm.p[2] / ts / ts; seG = fin(lm.stderr?.[1], NaN) / ts; its = lm.iterations; } } catch { /* keep the linearised estimate */ } }
  }
  const ok = valid && G > 0, D = ok ? G / (q * q) : NaN, d = ok ? N.K_B * T / (3 * PI * eta * D) : NaN, pdi = ok ? Math.max(0, mu2 / (G * G)) : NaN, urel = ok && seG > 0 ? seG / G : NaN;
  o.K('d_z', 'z-average hydrodynamic diameter', fin(d * 1e9), 'nm', `Stokes–Einstein diameter from the first cumulant Γ = ${F(fin(G), 4)} s⁻¹ at q = ${F(q, 4)} m⁻¹; intensity-weighted (∝ d⁶), so a few large particles or aggregates dominate. Fit uncertainty ± ${F(fin(urel * d * 1e9), 2)} nm.`, ok);
  o.K('pdi', 'Polydispersity index', fin(pdi), '', 'Second cumulant normalised by Γ²; below 0.1 is near-monodisperse, above about 0.3 the cumulant method is no longer meaningful and a distribution inversion is needed.', ok ? pdi < 0.3 : false);
  o.K('D', 'Translational diffusion coefficient', fin(D), 'm²/s', 'D = Γ/q² from the decay rate of the field correlation function; valid for dilute, non-interacting spheres in single-scattering conditions.');
  o.K('width', 'Distribution width (σ ≈ d·√PDI)', fin(d * Math.sqrt(pdi) * 1e9), 'nm', 'Standard deviation of the equivalent Gaussian intensity distribution implied by the polydispersity index; an estimate, not a measured distribution.');
  o.K('beta', 'Intercept (coherence factor β)', fin(beta), '', 'Amplitude of g₂ − B at zero delay; values well below about 0.5 indicate poor coherence, multiple scattering or a wrong baseline.', beta > 0.3 && beta <= 1.05);
  o.K('baseline', 'Baseline B', B, '', `Mean of the last ${nb} delay channels; should equal 1 for a normalised correlation function. A drifting baseline signals dust or sedimentation.`, Math.abs(B - 1) < 0.02);
  o.K('gamma', 'Decay rate Γ', fin(G), '1/s', 'First cumulant of ln g₁(τ): the mean decay rate of the scattered-field correlation, proportional to the diffusion coefficient.');
  const gf = tau.map(x => (ok ? B + beta * Math.exp(-2 * G * x) * (1 + mu2 * x * x / 2) ** 2 : B));
  o.S('g2', 'Intensity autocorrelation and cumulant fit', 'delay τ (s)', 'g₂(τ)', [['measured', tau, g2], ['cumulant fit', tau.slice(0, Math.max(m, 2)), gf.slice(0, Math.max(m, 2))]], 'The fit covers the decay down to 10 % of the intercept; deviations at long delays indicate a second, slower population.', { logx: true });
  o.S('lng', 'Linearised decay ln(g₂ − B)', 'delay τ (s)', 'ln(g₂ − B)', [['measured', tau.slice(0, m), g2.slice(0, m).map(v => Math.log(Math.max(1e-12, v - B)))], ['quadratic cumulant', tau.slice(0, m), tau.slice(0, m).map(x => Math.log(Math.max(1e-12, beta)) - 2 * fin(G) * x + mu2 * x * x)]], 'A straight line means a single decay rate; upward curvature is the polydispersity term μ₂τ².');
  o.E('governing', 'Siegert relation and cumulant expansion', 'g₂(τ) = B + β·exp(−2Γτ)·(1 + μ₂τ²/2)²'); o.E('formula', 'Scattering vector', 'q = (4π·n/λ)·sin(θ/2)'); o.E('formula', 'Stokes–Einstein', 'd_h = k_B·T/(3π·η·D),  D = Γ/q²,  PDI = μ₂/Γ²');
  o.C('Decay resolved', ok, valid ? `${m} channels above 10 % of the intercept` : 'fewer than 5 channels before the correlation decays — the delay range does not cover the decay');
  o.C('Baseline normalised', Math.abs(B - 1) < 0.02, `B = ${F(B, 5)}`);
  if (!ok) o.W('The correlation function does not show a usable decay; diameter and PDI are reported as 0.');
  if (ok && pdi > 0.3) o.W('PDI above 0.3: the cumulant z-average is not a reliable size; use a regularised distribution inversion (CONTIN/NNLS) and microscopy.');
  o.W('Assumes dilute non-interacting spheres and a known dispersant viscosity at the measurement temperature; viscosity error maps one-to-one into size error.');
  o.R('high', `Verify the dispersant viscosity (${inp.dls_visc} mPa·s at ${inp.dls_T} °C): a 2 % viscosity error shifts the ${F(fin(d * 1e9), 3)} nm z-average by 2 %.`, 'Stokes–Einstein sensitivity ∂d/∂η = −d/η.');
  o.R('medium', pdi > 0.15 ? 'Compare with TEM/SEM number statistics: at this polydispersity the intensity-weighted z-average overstates the number-mean size.' : 'Run at least three repeat measurements and report the repeatability standard deviation of the z-average alongside the fit error.', 'Cumulant PDI and intensity weighting.');
  Object.assign(o.publish, ok ? { 'geom.diameter': d, 'prop.D_diff': D, 'geom.pdi': pdi } : {}, Number.isFinite(urel) ? { 'cred.u_rel': urel } : {}, { 'cred.rmse': N.rmse(gf.slice(0, Math.max(m, 1)), g2.slice(0, Math.max(m, 1))) });
  o.numerics = { method: 'Cumulant analysis: quadratic fit of ln(g₂ − B) refined by Levenberg–Marquardt on g₂', dof: m, iterations: its, notes: 'baseline from the last 10 % of channels' };
  o.summary = ok ? `Cumulant analysis gives a z-average hydrodynamic diameter of ${F(d * 1e9, 4)} nm with PDI ${F(pdi, 3)} (D = ${F(D, 3)} m²/s at ${inp.dls_T} °C). The value is intensity-weighted and assumes spherical, non-interacting particles.` : 'The autocorrelation data do not contain a resolvable decay, so no hydrodynamic size could be derived; extend the delay range or check the baseline.';
}

function mCounts(inp, o) {
  const t = o.tbl(inp.tem_table, 1, 5, 'Particle diameters'), d = t.rows.map(r => r[0]).filter(v => v > 0), n = d.length;
  if (n < 5) throw new Error('Particle diameters: at least 5 positive values are needed.');
  const m = N.mean(d), s = N.std(d), cv = s / m, tt = tq(0.975, n - 1), ci = tt * s / Math.sqrt(n), ln = d.map(Math.log), mu = N.mean(ln), sg = N.std(ln), s3 = N.sum(d.map(v => v ** 3)), d43 = N.sum(d.map(v => v ** 4)) / s3, d32 = s3 / N.sum(d.map(v => v * v)), dv = Math.cbrt(s3 / n);
  const need = Math.ceil((1.96 * cv / (inp.tem_target / 100)) ** 2), sorted = [...d].sort((a, b) => a - b), ks = f => { let k = 0; sorted.forEach((v, i) => { const c = f(v); k = Math.max(k, Math.abs(c - i / n), Math.abs(c - (i + 1) / n)); }); return k; };
  const ksL = ks(v => N.normCdf(Math.log(v), mu, sg || 1e-12)), ksN = ks(v => N.normCdf(v, m, s || 1e-12)), h = N.histogram(d, Math.max(6, Math.min(30, Math.round(Math.sqrt(n)))));
  o.K('mean', 'Number-mean diameter', m, 'nm', `Arithmetic mean of ${n} counted particles; 95 % confidence interval of the mean ± ${F(ci, 3)} nm (Student-t, ${n - 1} degrees of freedom).`, ci / m < 0.05);
  o.K('sd', 'Standard deviation', s, 'nm', `Sample standard deviation of the diameters; coefficient of variation ${F(100 * cv, 3)} % describes the width of the population, not the measurement error.`);
  o.K('d43', 'Volume-weighted mean D[4,3]', d43, 'nm', 'Σd⁴/Σd³ — the mean a volume- or mass-sensitive method (laser diffraction, XRD volume weighting) would report; always at or above the number mean.');
  o.K('d32', 'Sauter mean D[3,2]', d32, 'nm', 'Σd³/Σd² — the diameter with the same surface-to-volume ratio as the population; compare with the BET equivalent diameter.');
  o.K('gmean', 'Geometric mean (log-normal median)', Math.exp(mu), 'nm', 'exp(mean ln d): the median of the fitted log-normal distribution obtained by the method of moments in log space.');
  o.K('gsd', 'Geometric standard deviation', Math.exp(sg), '', 'exp(std ln d); 1.0 is monodisperse, 1.1–1.3 is typical of good colloidal syntheses.', Math.exp(sg) < 1.3);
  o.K('need', 'Count needed for target precision', need, '', `Number of particles for a 95 % confidence half-width of ${inp.tem_target} % of the mean, n = (1.96·CV/ε)². ${n >= need ? 'The present count is sufficient.' : `The present count of ${n} is not sufficient.`}`, n >= need);
  o.K('ks', 'Log-normal goodness of fit (KS distance)', ksL, '', `Kolmogorov–Smirnov distance to the fitted log-normal; normal distribution gives ${F(ksN, 3)}. Approximate 5 % critical value 0.886/√n = ${F(0.886 / Math.sqrt(n), 3)} (parameters estimated from the data).`, ksL < 0.886 / Math.sqrt(n));
  const xs = N.linspace(sorted[0] * 0.8, sorted[n - 1] * 1.1, 120);
  o.S('hist', 'Size histogram and log-normal fit', 'diameter (nm)', 'count', [['counts', h.x, h.y], ['log-normal fit', xs, xs.map(v => n * h.width * Math.exp(-((Math.log(v) - mu) ** 2) / (2 * sg * sg || 1e-12)) / (v * (sg || 1e-6) * Math.sqrt(2 * PI)))]], 'Counted diameters binned by √n classes with the log-normal density scaled to the same area.');
  o.S('cdf', 'Cumulative number and volume distributions', 'diameter (nm)', 'cumulative fraction', [['number', sorted, sorted.map((_, i) => (i + 1) / n)], ['volume', sorted, (() => { let c = 0; return sorted.map(v => (c += v ** 3) / s3); })()]], 'The volume curve lies to the right of the number curve; the gap grows with polydispersity and explains why DLS and XRD report larger sizes than counting.');
  o.E('formula', 'Moment means', 'D[1,0] = Σd/n;  D[3,2] = Σd³/Σd²;  D[4,3] = Σd⁴/Σd³'); o.E('formula', 'Log-normal density', 'f(d) = exp(−(ln d − μ)²/(2σ²)) / (d·σ·√(2π))'); o.E('formula', 'Confidence interval and required count', 'CI₉₅ = t₀.₉₇₅,ₙ₋₁·s/√n;  n_req = (1.96·CV/ε)²');
  o.C('Enough particles counted', n >= need, `${n} counted, ${need} needed for ±${inp.tem_target} %`); o.C('Log-normal model acceptable', ksL < 0.886 / Math.sqrt(n), `KS distance ${F(ksL, 3)} against critical ${F(0.886 / Math.sqrt(n), 3)}`);
  if (n < 100) o.W('Fewer than 100 particles: the tails of the distribution (and D[4,3]) are poorly determined.');
  o.W('Projected-area diameters from 2-D images; overlapping or truncated particles and the image-scale calibration are not included in the confidence interval.');
  o.R(n >= need ? 'low' : 'high', n >= need ? `The ${n} counted particles already give ±${F(100 * ci / m, 2)} % on the mean; further counting is not needed for the ${inp.tem_target} % target.` : `Count at least ${need} particles (currently ${n}) from several grid areas to reach ±${inp.tem_target} % on the mean.`, 'Required-count formula with the observed CV.');
  o.R('medium', `Calibrate the image scale with a traceable grating: a 1 % pixel-size error shifts the ${F(m, 4)} nm mean by the same 1 %, larger than the statistical ± ${F(100 * ci / m, 2)} % here.`, 'Type-B uncertainty of the magnification.');
  Object.assign(o.publish, { 'geom.diameter': m * 1e-9, 'geom.diameter_sd': s * 1e-9, 'cred.u_rel': ci / 1.96 / m, 'geom.ssa_number': 6 / (d32 * 1e-9) });
  o.numerics = { method: 'Moment statistics, log-space method of moments, Kolmogorov–Smirnov distance', dof: n, notes: 'Student-t interval of the mean' };
  o.summary = `${n} particles give a number-mean diameter of ${F(m, 4)} ± ${F(ci, 2)} nm (95 % CI) with a standard deviation of ${F(s, 3)} nm; the log-normal fit has a median of ${F(Math.exp(mu), 4)} nm and geometric standard deviation ${F(Math.exp(sg), 3)}. The volume-weighted mean is ${F(d43, 4)} nm.`;
}

function mBET(inp, o) {
  const t = o.tbl(inp.bet_table, 2, 6, 'Adsorption isotherm'), rows = t.rows.filter(r => r[0] > 0 && r[0] < 1 && r[1] > 0).sort((a, b) => a[0] - b[0]);
  if (rows.length < 6) throw new Error('Adsorption isotherm: at least 6 rows with 0 < p/p0 < 1 and positive uptake are needed.');
  const x = rows.map(r => r[0]), v = rows.map(r => r[1]), n = x.length, rq = x.map((p, i) => v[i] * (1 - p));
  let iMax = 0; while (iMax < n - 1 && rq[iMax + 1] > rq[iMax] && x[iMax + 1] <= 0.5) iMax++;
  let idx = x.map((_, i) => i).filter(i => x[i] >= 0.05 && x[i] <= Math.min(0.35, x[iMax])); let widened = false;
  if (idx.length < 3) { idx = x.map((_, i) => i).filter(i => x[i] <= 0.35); widened = true; } if (idx.length < 3) idx = [0, 1, 2];
  const bx = idx.map(i => x[i]), by = idx.map(i => x[i] / (v[i] * (1 - x[i]))), f = linfit(bx, by), sum = f.slope + f.intercept, vm = sum > 0 ? 1 / sum : NaN, Cc = f.intercept !== 0 ? f.slope / f.intercept + 1 : NaN;
  const ssa = vm / 22414 * N.N_A * inp.bet_sigma * 1e-18, seVm = vm * vm * Math.sqrt(Math.max(0, f.seSlope ** 2 + f.seInt ** 2 - 2 * f.mx * f.seSlope ** 2)), urel = fin(seVm / vm, 0), xm = Cc > 0 ? 1 / (Math.sqrt(Cc) + 1) : NaN, xAt = N.interp1(v, x, vm);
  const dEq = 6 / (inp.rho * ssa * 1e3), vTot = v[n - 1] * N2.liq, dPore = 4 * vTot * 1e-6 / ssa, por = vTot * 1e-6 * inp.rho * 1e3 / (1 + vTot * 1e-6 * inp.rho * 1e3);
  // t-plot (Halsey thickness) and Kelvin/BJH-style distribution on the adsorption branch
  const ti = x.map((_, i) => i).filter(i => N2.tHJ(x[i]) >= 0.35 && N2.tHJ(x[i]) <= 0.5), tf = ti.length >= 3 ? linfit(ti.map(i => N2.tHJ(x[i])), ti.map(i => v[i])) : null, sExt = tf ? tf.slope * 1.547 : NaN, vMicro = tf ? Math.max(0, tf.intercept * N2.liq) : NaN;
  const pr = [], pd = []; let vw = 0, rw = 0;
  for (let i = 1; i < n; i++) if (x[i] > 0.35 && v[i] > v[i - 1]) { const xm2 = Math.sqrt(x[i] * x[i - 1]), r = N2.rK / -Math.log(xm2) + N2.t(xm2), r0 = N2.rK / -Math.log(x[i - 1]) + N2.t(x[i - 1]), r1 = N2.rK / -Math.log(x[i]) + N2.t(x[i]), dV = (v[i] - v[i - 1]) * N2.liq; if (r < 60) { pr.push(2 * r); pd.push(dV / Math.max(1e-9, 2 * (r1 - r0))); vw += dV; rw += dV * 2 * r; } }
  const dBJH = vw > 0 ? rw / vw : NaN, okC = Cc > 0, okX = Number.isFinite(xm) && xm >= bx[0] - 1e-9 && xm <= bx[bx.length - 1] + 1e-9, okM = Number.isFinite(xAt) && Math.abs(xAt - xm) / xm < 0.2;
  o.K('ssa', 'BET specific surface area', fin(ssa), 'm²/g', `Monolayer capacity v_m = ${F(fin(vm), 4)} cm³(STP)/g × molecular cross-section ${inp.bet_sigma} nm²; linear BET fit over ${idx.length} points in p/p₀ = ${F(bx[0], 2)}–${F(bx[bx.length - 1], 2)} (R² = ${F(f.r2, 5)}), ± ${F(100 * urel, 2)} % from the regression.`, okC && f.r2 > 0.999);
  o.K('C', 'BET constant C', fin(Cc), '', 'C = slope/intercept + 1 reflects the adsorbate–surface interaction energy; 50–200 is typical for N₂ on oxides. Negative C means the pressure window is invalid; C above ~500 suggests micropore filling.', okC && Cc < 500);
  o.K('vm', 'Monolayer capacity', fin(vm), 'cm³(STP)/g', 'Gas volume completing one statistical monolayer, 1/(slope + intercept) of the BET line.');
  o.K('d_eq', 'Equivalent particle diameter', fin(dEq * 1e9), 'nm', `d = 6/(ρ·S) for dense, smooth, monodisperse spheres with ρ = ${inp.rho} kg/m³; internal porosity or roughness makes it smaller than the real particle size, agglomerate necks make it larger.`);
  o.K('v_pore', 'Total pore volume', vTot, 'cm³/g', `Liquid-equivalent uptake at the highest relative pressure (p/p₀ = ${F(x[n - 1], 3)}), Gurvich rule; includes interparticle condensation.`);
  o.K('d_pore', 'Mean pore diameter 4V/S', fin(dPore * 1e9), 'nm', 'Hydraulic diameter of cylindrical pores having the measured total pore volume and BET area; shape-dependent and only an average.');
  o.K('d_bjh', 'Kelvin (BJH-style) mean pore diameter', fin(dBJH), 'nm', 'Volume-weighted mean of 2(r_K + t) from increments of the adsorption branch above p/p₀ = 0.35, without the wall-film correction of the full BJH recursion; meaningful only for 2–50 nm mesopores.');
  o.K('s_ext', 't-plot external area', fin(sExt), 'm²/g', tf ? `Slope of uptake against the Harkins–Jura film thickness between 0.35 and 0.50 nm (${ti.length} points); micropore volume from the intercept ${F(fin(vMicro), 3)} cm³/g. Agreement with the BET area indicates the absence of micropores.` : 'Fewer than three points fall in the 0.35–0.50 nm thickness window, so no t-plot was fitted.');
  o.K('porosity', 'Porosity implied by pore volume', por, '', 'ε = V_p·ρ/(1 + V_p·ρ) using the skeletal density; a bed or particle porosity only if all pores were filled at the highest pressure.');
  o.S('iso', 'Adsorption isotherm', 'p/p₀', 'uptake (cm³ STP/g)', [['measured', x, v], ['BET model', x.filter(p => p < 0.5), x.filter(p => p < 0.5).map(p => (okC ? vm * Cc * p / ((1 - p) * (1 - p + Cc * p)) : 0))]], 'Measured uptake with the BET isotherm from the fitted v_m and C; the model is only expected to hold in the fitted window.');
  o.S('bet', 'BET transform', 'p/p₀', '(p/p₀)/[v(1 − p/p₀)]', [['data in window', bx, by], ['linear fit', [bx[0], bx[bx.length - 1]], [f.intercept + f.slope * bx[0], f.intercept + f.slope * bx[bx.length - 1]]]], 'Linearised BET plot over the Rouquerol-consistent window; curvature means the window is too wide.');
  if (pr.length >= 2) o.S('psd', 'Pore-size distribution (Kelvin, adsorption branch)', 'pore diameter (nm)', 'dV/dD (cm³/g/nm)', [['dV/dD', pr, pd]], 'Incremental liquid volume against Kelvin diameter plus film thickness; a simplified BJH estimate.');
  o.E('governing', 'BET equation (linear form)', '(p/p₀)/[v(1 − p/p₀)] = 1/(v_m·C) + (C − 1)/(v_m·C)·(p/p₀)'); o.E('formula', 'Specific surface area', 'S_BET = v_m·N_A·σ/V_molar,  V_molar = 22 414 cm³/mol'); o.E('formula', 'Kelvin equation and film thickness (N₂, 77 K)', `r_K = −2γV_m/(R·T·ln(p/p₀)) = ${F(N2.rK, 4)} nm/(−ln(p/p₀));  t_Halsey = 0.354·[5/(−ln(p/p₀))]^(1/3) nm (Kelvin correction);  t_HJ = 0.1·[13.99/(0.034 − log₁₀(p/p₀))]^½ nm (t-plot)`); o.E('formula', 'Equivalent diameter', 'd = 6/(ρ·S_BET)');
  o.C('Rouquerol 1: v(1 − p/p₀) increases over the window', x[iMax] >= bx[bx.length - 1] - 1e-12 && !widened, `v(1 − p/p₀) rises up to p/p₀ = ${F(x[iMax], 3)}`); o.C('Rouquerol 2: C positive', okC, `C = ${F(fin(Cc), 4)}`);
  o.C('Rouquerol 3: monolayer pressure inside the window', okX, `1/(√C + 1) = ${F(fin(xm), 3)}`); o.C('Rouquerol 4: monolayer pressure consistent', okM, `p/p₀ at v_m read from the isotherm ${F(fin(xAt), 3)} against ${F(fin(xm), 3)}`);
  if (widened) o.W('Fewer than three points lie in the consistent 0.05–0.35 window; the fit was widened and the surface area is unreliable.');
  if (Cc > 500) o.W('Very high C constant: microporosity is likely, the BET area is then an apparent value.');
  o.W('Assumes N₂ at 77 K (Kelvin constant, film thickness and liquid density are for nitrogen). The pore-size estimate neglects the BJH wall-film recursion and uses the adsorption branch only.');
  o.R('high', okC && okX && okM ? `Report S_BET = ${F(ssa, 4)} m²/g with the fitted window (${F(bx[0], 2)}–${F(bx[bx.length - 1], 2)}) and C = ${F(Cc, 3)}; the consistency criteria are met.` : 'Re-select the pressure range: at least one Rouquerol consistency criterion fails, so the reported area depends on the window.', 'Rouquerol consistency checks.');
  o.R('medium', `Compare the equivalent diameter ${F(fin(dEq * 1e9), 3)} nm with TEM and XRD sizes: a much smaller BET diameter indicates internal porosity or roughness, a larger one hard agglomeration.`, 'd = 6/(ρ·S) cross-check.');
  o.R('low', 'Recover or minimise liquid-nitrogen use by batching degassed samples in one dewar fill and re-using the free-space calibration where the sample tube is unchanged.', 'Cryogen consumption per isotherm.', true);
  Object.assign(o.publish, { 'geom.ssa': fin(ssa) * 1e3, 'geom.diameter': fin(dEq), 'geom.porosity': por, 'geom.pore_diameter': fin(dPore), 'cred.u_rel': urel, 'cred.rmse': f.s });
  o.numerics = { method: 'Linear BET regression over a Rouquerol-consistent window; Harkins–Jura t-plot; Kelvin increments with Halsey film thickness', dof: idx.length, residual: f.s, notes: `${n} isotherm points` };
  o.summary = `The BET fit over p/p₀ = ${F(bx[0], 2)}–${F(bx[bx.length - 1], 2)} gives a specific surface area of ${F(fin(ssa), 4)} m²/g (C = ${F(fin(Cc), 3)}), equivalent to dense spheres of ${F(fin(dEq * 1e9), 3)} nm. Total pore volume is ${F(vTot, 3)} cm³/g with a mean pore diameter of ${F(fin(dPore * 1e9), 3)} nm (4V/S).`;
}

function mTauc(inp, o) {
  const t = o.tbl(inp.uv_table, 2, 12, 'Absorbance spectrum'), rows = t.rows.filter(r => r[0] > 0).map(r => [1239.842 / r[0], r[1]]).sort((a, b) => a[0] - b[0]), E = rows.map(r => r[0]), A = rows.map(r => r[1]), n = E.length, ex = inp.tauc_type === 'indirect' ? 0.5 : 2;
  const base = N.quantile(A, 0.03), raw = A.map((a, i) => Math.pow(Math.max(0, (a - (inp.tauc_base === 'subtract' ? base : 0)) * E[i]), ex)), ymax = N.maxOf(raw) || 1, y = raw.map(v => v / ymax), w = Math.max(6, Math.round(n * 0.08));
  const wins = []; for (let i = 0; i + w <= n; i++) { const f = linfit(E.slice(i, i + w), y.slice(i, i + w)), mid = N.mean(y.slice(i, i + w)); if (f.slope > 0 && mid > 0.1 && mid < 0.9) wins.push({ i, f, Eg: -f.intercept / f.slope }); }
  const good = wins.filter(q => q.f.r2 > 0.985), pool = good.length ? good : wins, bestW = pool.length ? pool.reduce((a, b) => (b.f.slope > a.f.slope ? b : a)) : null;
  let Eg = NaN, uFit = NaN, uSel = NaN, EgB = NaN, Eu = NaN;
  if (bestW) {
    const f = bestW.f; Eg = bestW.Eg; uFit = f.s / f.slope * Math.sqrt(1 / f.n + (f.sxx > 0 ? (Eg - f.mx) ** 2 / f.sxx : 0));
    const near = pool.filter(q => q.f.slope > 0.8 * f.slope).map(q => q.Eg); uSel = near.length > 1 ? N.std(near) : 0;
    const lo = E.map((_, i) => i).filter(i => E[i] < Eg - 0.25); if (lo.length >= 4) { const fb = linfit(lo.map(i => E[i]), lo.map(i => y[i])); if (f.slope - fb.slope > 0) EgB = (fb.intercept - f.intercept) / (f.slope - fb.slope); }
    const nz = noiseOf(A), ui = E.map((_, i) => i).filter(i => E[i] > Eg - 0.3 && E[i] < Eg - 0.02 && A[i] - base > 4 * nz); if (ui.length >= 4) { const fu = linfit(ui.map(i => E[i]), ui.map(i => Math.log(A[i] - base))); if (fu.slope > 0) Eu = 1 / fu.slope; }
  }
  const u = Math.hypot(fin(uFit), fin(uSel), Number.isFinite(EgB) ? Math.abs(EgB - Eg) / Math.sqrt(3) : 0), ok = Number.isFinite(Eg) && Eg > E[0] - 1 && Eg < E[n - 1];
  o.K('Eg', 'Optical band gap (Tauc)', fin(Eg), 'eV', `x-axis intercept of the steepest linear segment of (A·hν)^${ex} (${inp.tauc_type} allowed transition), chosen automatically among ${wins.length} sliding windows of ${w} points with R² > 0.985. Combined standard uncertainty ± ${F(u, 2)} eV.`, ok);
  o.K('u_Eg', 'Band-gap standard uncertainty', u, 'eV', `Root-sum-square of the regression error (${F(fin(uFit), 2)} eV), the spread between equally acceptable windows (${F(fin(uSel), 2)} eV) and the baseline-method difference; it excludes the choice of transition type, which usually dominates.`, u < 0.05);
  o.K('Eg_base', 'Band gap with baseline intersection', fin(EgB, fin(Eg)), 'eV', 'Intersection of the Tauc line with a linear fit to the sub-gap region instead of the zero line; the preferred estimate when sub-gap absorption or scattering lifts the baseline.');
  o.K('lambda_g', 'Absorption edge wavelength', ok ? 1239.842 / Eg : 0, 'nm', 'λ = hc/E_g = 1239.84 eV·nm / E_g: the wavelength below which the material absorbs strongly.');
  o.K('r2', 'Linearity of the selected segment', bestW ? bestW.f.r2 : 0, '', `Coefficient of determination of the chosen Tauc segment (${bestW ? `${F(E[bestW.i], 4)}–${F(E[bestW.i + w - 1], 4)} eV` : 'none found'}); below about 0.99 the transition-type assumption or the baseline is doubtful.`, bestW ? bestW.f.r2 > 0.99 : false);
  o.K('Eu', 'Urbach energy', fin(Eu * 1000), 'meV', 'Inverse slope of ln(A) just below the gap; a measure of band-tail disorder (defects, dopants, thermal broadening). Zero when the tail is not resolved.');
  o.S('tauc', 'Tauc plot', 'photon energy (eV)', `(A·hν)^${ex} (normalised)`, [['data', E, y], ...(bestW ? [['linear segment', [Eg, E[bestW.i + w - 1]], [0, bestW.f.intercept + bestW.f.slope * E[bestW.i + w - 1]]]] : [])], 'The straight segment extrapolated to zero gives the band gap; the fit window is selected for maximum slope with high linearity.');
  o.S('abs', 'Absorbance spectrum', 'wavelength (nm)', 'absorbance', [['measured', rows.map(r => 1239.842 / r[0]).reverse(), [...A].reverse()]], 'Raw absorbance against wavelength; a sloping long-wavelength background indicates scattering that should be removed before the Tauc analysis.');
  o.E('governing', 'Tauc relation', '(α·hν)^(1/r) = B·(hν − E_g);  r = ½ direct allowed, r = 2 indirect allowed'); o.E('formula', 'Photon energy', 'hν (eV) = 1239.842/λ (nm)'); o.E('formula', 'Urbach tail', 'α = α₀·exp((hν − E₀)/E_U)');
  o.C('Linear Tauc segment found', !!bestW && bestW.f.r2 > 0.985, bestW ? `R² = ${F(bestW.f.r2, 5)} over ${w} points` : 'no rising segment'); o.C('Band gap inside the measured range', ok, `E_g = ${F(fin(Eg), 4)} eV; data span ${F(E[0], 3)}–${F(E[n - 1], 3)} eV`);
  if (!ok) o.W('No band edge was found inside the measured range; the band gap is reported as 0.');
  o.W('The Tauc method assumes parabolic bands and a known transition type; absorbance is used in place of the absorption coefficient, which is valid for a uniform non-scattering film or dilute dispersion. For powders use Kubelka–Munk transformed diffuse reflectance.');
  o.R('high', `State the transition type with the result: E_g = ${F(fin(Eg), 4)} ± ${F(u, 1)} eV holds for a ${inp.tauc_type} allowed transition only; re-run with the other exponent to see the model dependence.`, 'Tauc exponent is an assumption, not an output.');
  o.R('medium', Number.isFinite(EgB) && Math.abs(EgB - Eg) > 0.03 ? `Baseline and zero-line intercepts differ by ${F(Math.abs(EgB - Eg), 2)} eV: remove the scattering background (integrating sphere) before quoting the gap.` : 'Measure with an integrating sphere if the sample scatters; the present baseline influence is small.', 'Difference between the two intercept constructions.');
  Object.assign(o.publish, ok ? { 'prop.bandgap_eV': Eg, 'cred.u_bandgap_eV': u, 'cred.u_rel': u / Eg } : {}, { 'cred.rmse': bestW ? bestW.f.s : 0 });
  o.numerics = { method: 'Sliding-window linear regression on the Tauc transform with automatic segment selection', dof: n, iterations: wins.length, residual: bestW ? bestW.f.s : 0 };
  o.summary = ok ? `The Tauc construction for a ${inp.tauc_type} allowed transition gives an optical band gap of ${F(Eg, 4)} ± ${F(u, 1)} eV (absorption edge ${F(1239.842 / Eg, 4)} nm). The Urbach tail energy is ${F(fin(Eu * 1000), 3)} meV.` : 'No linear rising Tauc segment with an intercept inside the measured range was found, so no band gap is reported; extend the spectrum across the absorption edge.';
}

function mLSPR(inp, o) {
  const t = o.tbl(inp.lspr_table, 2, 10, 'Extinction spectrum'), rows = t.rows.sort((a, b) => a[0] - b[0]), l = rows.map(r => r[0]), A = rows.map(r => r[1]), n = l.length, au = inp.lspr_metal !== 'Ag', As = smooth(A, 2);
  const lo = au ? 480 : 370, hi = au ? 650 : 520, cand = l.map((_, i) => i).filter(i => l[i] >= lo && l[i] <= hi), ip = cand.length ? cand.reduce((a, b) => (As[b] > As[a] ? b : a)) : N.argmax(As);
  let lp = l[ip], Ap = As[ip]; const i0 = Math.max(0, ip - 4), i1 = Math.min(n - 1, ip + 4);
  if (i1 - i0 >= 4) { const c = N.polyfit(l.slice(i0, i1 + 1).map(v => v - l[ip]), A.slice(i0, i1 + 1), 2); if (c[2] < 0) { const dx = N.clamp(-c[1] / (2 * c[2]), l[i0] - l[ip], l[i1] - l[ip]); lp = l[ip] + dx; Ap = c[0] + c[1] * dx + c[2] * dx * dx; } }
  const isPeak = ip > 0 && ip < n - 1 && As[ip] > As[Math.max(0, ip - 5)] && As[ip] > As[Math.min(n - 1, ip + 5)], A450 = N.interp1(l, As, 450), A650 = N.interp1(l, As, au ? 650 : 550), ratio = A450 > 0 ? Ap / A450 : NaN;
  let r = ip; while (r < n - 1 && As[r] > Ap / 2) r++; const hwhm = l[r] - lp;
  const dAu = Math.exp(3.0 * ratio - 2.2), dAg = (() => { const p = AG_LSPR; if (lp <= p[0][1]) return p[0][0]; if (lp >= p[p.length - 1][1]) return p[p.length - 1][0]; return N.interp1(p.map(q => q[1]), p.map(q => q[0]), lp); })(), d = au ? dAu : dAg, dOK = isPeak && Number.isFinite(d) && d > 3 && d < 120;
  const c = Ap / (inp.lspr_eps * inp.lspr_path), rho = au ? 19300 : 10490, mP = rho * PI * (d * 1e-9) ** 3 / 6, cMass = c * 1e3 * N.N_A * mP;
  o.K('lambda', 'LSPR peak wavelength', lp, 'nm', `Maximum of the extinction band located by a parabola through nine points around the highest value in ${lo}–${hi} nm. ${isPeak ? '' : 'No distinct maximum was found in this window — the value is the window maximum only.'}`, isPeak);
  o.K('A_peak', 'Peak absorbance', Ap, '', 'Extinction at the plasmon maximum; keep it between about 0.1 and 1.5 for a linear detector response.', Ap > 0.05 && Ap < 1.5);
  o.K('d_est', `${au ? 'Au' : 'Ag'} particle diameter (empirical estimate)`, dOK ? d : 0, 'nm', au ? `Empirical correlation for citrate-type gold spheres in water, d ≈ exp(3.00·A_spr/A_450 − 2.20) with A_spr/A_450 = ${F(fin(ratio), 4)}; typical scatter 10–20 %, intended for about 5–80 nm. Treat as a screening value, not a measurement.` : 'Interpolated from indicative peak positions of silver spheres in water (≈ 395 nm at 10 nm to ≈ 490 nm at 100 nm); uncertainty of several tens of percent because shape, coating and medium shift the peak as much as size does.', dOK ? null : false);
  o.K('conc', 'Particle molar concentration', c, 'mol/L', `Beer–Lambert law c = A/(ε·l) with ε = ${F(inp.lspr_eps, 3)} L/mol/cm and l = ${inp.lspr_path} cm; ε of metal nanoparticles scales roughly with d³, so the value is only as good as the assumed ε.`);
  o.K('n_conc', 'Number concentration', c * N.N_A * 1e3, '1/m³', 'Molar particle concentration times Avogadro constant, per cubic metre of dispersion.');
  o.K('mass_conc', 'Metal mass concentration', dOK ? cMass : 0, 'kg/m³', `Number concentration × mass of one sphere of the estimated diameter with bulk density ${rho} kg/m³; inherits the cube of the size error.`);
  o.K('hwhm', 'Band half-width (red side)', fin(hwhm), 'nm', 'Half width at half maximum on the long-wavelength side; broadening beyond about 45 nm for small gold spheres indicates polydispersity or aggregation.');
  o.K('agg', 'Aggregation indicator', Ap > 0 ? A650 / Ap : 0, '', `Ratio of extinction at ${au ? 650 : 550} nm to the peak; rising values during storage or after salt addition signal flocculation (coupled plasmons absorb at longer wavelength).`, Ap > 0 ? A650 / Ap < 0.35 : null);
  o.S('spec', 'Extinction spectrum', 'wavelength (nm)', 'absorbance', [['measured', l, A], ['smoothed', l, As]], 'Extinction of the colloid with the plasmon band; a shoulder or tail to the red indicates anisotropic particles or aggregates.');
  o.S('norm', 'Spectrum normalised to the peak', 'wavelength (nm)', 'A/A_peak', [['normalised', l, As.map(v => v / (Ap || 1))]], 'Normalisation removes concentration so that band position and width can be compared between batches.');
  o.E('governing', 'Beer–Lambert law', 'A = ε·c·l'); o.E('formula', au ? 'Empirical size correlation (gold spheres in water)' : 'Indicative size–peak relation (silver spheres in water)', au ? 'd (nm) ≈ exp(3.00·A_spr/A_450 − 2.20)' : 'λ_spr ≈ 395 nm (10 nm) … 490 nm (100 nm), interpolated', 'Empirical; no traceable calibration is implied.'); o.E('formula', 'Quasi-static dipole resonance condition', 'Re ε_metal(λ_spr) = −2·ε_medium');
  try { const eq = au ? chem.ionicEquations([{ f: 'AuCl4-', phase: 'aq' }, { f: 'e-', phase: null }], [{ f: 'Au', phase: 's' }, { f: 'Cl-', phase: 'aq' }]) : chem.ionicEquations([{ f: 'Ag+', phase: 'aq' }, { f: 'e-', phase: null }], [{ f: 'Ag', phase: 's' }]); o.E('half-cell', 'Reduction half-reaction forming the colloid', eq.netIonic); } catch { /* optional */ }
  o.C('Distinct plasmon maximum', isPeak, `peak at ${F(lp, 4)} nm`); o.C('Absorbance in the linear range', Ap > 0.05 && Ap < 1.5, `A_peak = ${F(Ap, 3)}`); o.C('Size inside the correlation range', dOK, `d = ${F(fin(d), 3)} nm`);
  o.W('The size estimate is an empirical correlation for near-spherical, non-aggregated particles in water; ligand shell, shape and refractive index of the medium shift the band by amounts comparable to the size effect. Confirm by TEM or DLS.');
  o.R('high', `Use the optical size (${F(dOK ? d : 0, 3)} nm) only for batch-to-batch screening and calibrate it once against TEM counts of the same synthesis.`, 'Empirical nature of the LSPR size relation.');
  o.R('medium', Ap > 1.2 ? 'Dilute the sample: absorbance above about 1.2 leaves the linear range of most spectrometers and biases the A_spr/A_450 ratio.' : 'Record the spectrum of the dispersant in the same cuvette and subtract it; the ratio method is sensitive to baseline offsets at 450 nm.', 'Photometric linearity and baseline.');
  o.R('low', 'Use micro-volume cuvettes so that only a few tens of microlitres of precious-metal colloid are consumed per measurement, and return the sample to the batch.', 'Gold and silver are critical, high-footprint raw materials.', true);
  Object.assign(o.publish, dOK ? { 'geom.diameter': d * 1e-9, 'cred.u_rel': au ? 0.15 : 0.4 } : {}, { 'res.lspr_nm': lp, 'proc.c_particles': c * 1e3, 'material.formula': au ? 'Au' : 'Ag' });
  o.numerics = { method: 'Parabolic peak location, Beer–Lambert reduction and empirical size correlation', dof: n };
  o.summary = `The plasmon band peaks at ${F(lp, 4)} nm with absorbance ${F(Ap, 3)}; the empirical ${au ? 'A_spr/A_450 ratio' : 'peak-position'} correlation suggests ${au ? 'gold' : 'silver'} spheres of about ${F(dOK ? d : 0, 3)} nm. Beer–Lambert gives ${F(c, 3)} mol/L of particles for the stated extinction coefficient; both numbers are estimates to be confirmed by microscopy.`;
}

function mRaman(inp, o) {
  const t = o.tbl(inp.raman_table, 2, 20, 'Raman/FTIR spectrum'), rows = t.rows.sort((a, b) => a[0] - b[0]), x = rows.map(r => r[0]), y = rows.map(r => r[1]), n = x.length, step = (x[n - 1] - x[0]) / Math.max(1, n - 1) || 1;
  const bg = snip(smooth(y, 2), N.clamp(Math.round(80 / step), 6, 300)), net = y.map((v, i) => v - bg[i]), noise = noiseOf(y), fp = fitPeaks(x, net, noise, { minSep: 6, minRel: 0.02, maxPeaks: 20 }), pk = fp.peaks, tol = inp.raman_tol;
  const asg = pk.map(p => RAMAN_MODES.filter(m => Math.abs(m[2] - p.pos) <= tol * (m[2] > 2000 ? 3 : 1)).sort((a, b) => Math.abs(a[2] - p.pos) - Math.abs(b[2] - p.pos)));
  const find = name => { let b = null; pk.forEach((p, i) => { const m = asg[i].find(q => q[1].startsWith(name)); if (m && (!b || p.h > b.h)) b = p; }); return b; };
  const D = find('D'), G = find('G'), D2 = find('2D'), Si = pk.find((p, i) => asg[i].some(m => m[0].startsWith('Si'))), idig = D && G ? D.h / G.h : G ? 0 : NaN, i2dg = D2 && G ? D2.h / G.h : NaN, La = idig > 0 ? 2.4e-10 * inp.raman_laser ** 4 / idig : NaN;
  const dSi = Si ? Si.pos - 520.7 : NaN, stress = Si ? -inp.raman_stress * dSi : NaN, mats = {}; pk.forEach((p, i) => { if (asg[i][0]) mats[asg[i][0][0]] = (mats[asg[i][0][0]] || 0) + 1; });
  const ranked = Object.entries(mats).map(([m, c]) => [m, c, RAMAN_MODES.filter(q => q[0] === m && q[2] >= x[0] && q[2] <= x[n - 1]).length]).sort((a, b) => b[1] / b[2] - a[1] / a[2] || b[1] - a[1]), amb = pk.filter((p, i) => new Set(asg[i].map(m => m[0])).size > 1).length;
  o.K('n_peaks', 'Fitted bands', pk.length, '', `Bands above 4σ noise after SNIP baseline removal, fitted with pseudo-Voigt profiles; ${pk.filter((p, i) => asg[i].length).length} of them lie within ${tol} cm⁻¹ of a library mode.`, pk.length > 0);
  o.K('top', 'Most consistent material', ranked[0] ? ranked[0][0] : 'none', '', ranked[0] ? `${ranked[0][1]} of ${ranked[0][2]} library modes of this material inside the measured range were found. The built-in list holds only a dozen well-known modes, so this is a hint and not an identification.` : 'No fitted band coincides with a built-in mode position within the tolerance.', ranked[0] ? ranked[0][1] / ranked[0][2] > 0.6 : false);
  o.K('idig', 'Graphitic defect ratio I_D/I_G', fin(idig), '', D && G ? `Height ratio of the D band (${F(D.pos, 5)} cm⁻¹) to the G band (${F(G.pos, 5)} cm⁻¹); grows with defect density in the low-defect regime and falls again in highly disordered carbon, so it is not monotonic.` : 'D and/or G band not present — reported as 0; not a carbon spectrum or no defect band above noise.', D && G ? idig < 0.5 : null);
  o.K('La', 'In-plane crystallite size L_a', fin(La), 'nm', `L_a = 2.4×10⁻¹⁰·λ⁴/(I_D/I_G) with λ = ${inp.raman_laser} nm (empirical relation for nanographite, valid for L_a above a few nm); 0 when no D band is resolved.`);
  o.K('i2dg', 'I_2D/I_G ratio', fin(i2dg), '', 'Height ratio of the 2D to the G band; above about 2 with a narrow single-Lorentzian 2D band is characteristic of monolayer graphene, below 1 of multilayers — substrate and doping also change it.');
  o.K('si_shift', 'Si peak shift from 520.7 cm⁻¹', fin(dSi), 'cm⁻¹', Si ? `Fitted position ${F(Si.pos, 6)} ± ${F(Si.sePos, 1)} cm⁻¹ relative to stress-free crystalline silicon; requires a wavenumber calibration better than the shift itself.` : 'No silicon band found within the tolerance — reported as 0.');
  o.K('si_stress', 'Si stress from peak shift', fin(stress), 'MPa', `σ = −${inp.raman_stress} MPa per cm⁻¹ × Δω (user coefficient: about 500 for uniaxial and 250 for equi-biaxial stress in (100) silicon). Positive = tensile. Laser heating also shifts the peak down and mimics tension.`);
  o.K('noise', 'Noise level (1σ)', noise, 'counts', 'Robust estimate from second differences of the raw spectrum; sets the 4σ detection threshold for bands.');
  o.S('spec', 'Spectrum, baseline and fitted bands', 'Raman shift (cm⁻¹)', 'intensity', [['measured', x, y], ['baseline', x, bg], ['baseline + fit', x, fp.fit.map((v, i) => v + bg[i])]], 'Raw spectrum with the clipped baseline and the sum of fitted bands.');
  o.S('net', 'Baseline-corrected spectrum', 'Raman shift (cm⁻¹)', 'net intensity', [['net', x, net], ['fit', x, fp.fit]], 'Net band intensities used for the ratios; residual structure indicates overlapping or asymmetric bands.');
  o.T('bands', 'Band table and assignments', ['position (cm⁻¹)', '± (cm⁻¹)', 'FWHM (cm⁻¹)', 'height', 'area', 'assignment', 'Δ to library (cm⁻¹)', 'alternatives'], pk.map((p, i) => [p.pos, p.sePos, p.fwhm, p.h, p.area, asg[i][0] ? `${asg[i][0][0]} ${asg[i][0][1]}` : 'unassigned', asg[i][0] ? p.pos - asg[i][0][2] : '—', asg[i].slice(1).map(m => `${m[0]} ${m[1]}`).join('; ') || '—']), 'Each fitted band with the nearest built-in mode; alternatives are listed where two materials have modes within the tolerance.');
  o.E('formula', 'Raman shift', 'Δν̃ = 1/λ_laser − 1/λ_scattered'); o.E('formula', 'Graphitic crystallite size (empirical)', 'L_a (nm) = 2.4×10⁻¹⁰·λ_laser⁴ (nm⁴) / (I_D/I_G)'); o.E('constitutive', 'Stress from the silicon phonon shift', `σ (MPa) = −${inp.raman_stress}·Δω (cm⁻¹)`, 'Coefficient depends on stress state and orientation.'); o.E('formula', 'Lorentzian band shape', 'I(ν̃) = h/(1 + 4(ν̃ − ν̃₀)²/Γ²)');
  o.C('Bands detected', pk.length > 0, `${pk.length} bands, noise ${F(noise, 3)}`); o.C('Assignments unambiguous', amb === 0, amb ? `${amb} band(s) match modes of more than one material (e.g. ZnO 437 and rutile 447 cm⁻¹)` : 'every assigned band matches a single library material');
  if (!pk.length) o.W('No band above the noise threshold was found.');
  o.W('Mode matching uses a short built-in list of well-known positions; real identification needs a full reference spectrum, polarisation and the complete set of modes.');
  o.R('high', Si ? `Calibrate the wavenumber axis on a stress-free Si reference before interpreting the ${F(dSi, 2)} cm⁻¹ shift: a 0.5 cm⁻¹ calibration error equals ${F(0.5 * inp.raman_stress, 3)} MPa.` : 'Calibrate the wavenumber axis with a silicon wafer (520.7 cm⁻¹) before comparing band positions with reference values.', 'Stress coefficient × calibration uncertainty.');
  o.R('medium', D && G ? `I_D/I_G = ${F(idig, 3)}: map at least 10 points to separate local defects from the batch average, and keep the laser power low enough that the ratio does not change with exposure.` : 'Reduce laser power or defocus if band positions drift with time — heating shifts and broadens phonon bands.', 'Defect ratio variability and laser heating.');
  Object.assign(o.publish, { 'cred.rmse': noise }, D && G ? { 'defect.id_ig': idig } : {}, Si ? { 'mech.residual_stress': stress * 1e6 } : {}, ranked[0] ? { 'material.name': ranked[0][0] } : {});
  o.numerics = { method: 'SNIP baseline + pseudo-Voigt band fits (Levenberg–Marquardt) + nearest-mode assignment', dof: n, iterations: fp.iterations };
  o.summary = `${pk.length} bands were fitted; the built-in mode list is most consistent with ${ranked[0] ? ranked[0][0] : 'no library material'}.${D && G ? ` The graphitic defect ratio I_D/I_G is ${F(idig, 3)}${Number.isFinite(i2dg) ? ` and I_2D/I_G is ${F(i2dg, 3)}` : ''}.` : ''}${Si ? ` The silicon band sits ${F(dSi, 2)} cm⁻¹ from 520.7 cm⁻¹, equivalent to ${F(stress, 3)} MPa with the stated coefficient.` : ''} Assignments are hints that need a full reference spectrum.`;
}

function mOP(inp, o) {
  const t = o.tbl(inp.op_table, 2, 12, 'Load–displacement curve'), h = t.rows.map(r => r[0]), P = t.rows.map(r => r[1]), n = h.length, iP = N.argmax(P), Pmax = P[iP];
  let iU = iP; for (let i = iP; i < n; i++) if (P[i] >= 0.995 * Pmax) iU = i;                       // end of any hold segment
  const hmax = N.maxOf(h.slice(0, iU + 1)), un = []; for (let i = iU; i < n; i++) if (P[i] > 0.2 * Pmax && P[i] <= Pmax && h[i] > 0) un.push([h[i], P[i]]);
  if (un.length < 5 || !(Pmax > 0)) throw new Error('Load–displacement curve: an unloading segment with at least 5 points between 20 % and 100 % of the peak load is needed after the maximum.');
  const hminU = N.minOf(un.map(q => q[0])), hLast = Math.min(h[n - 1], hminU), hf0 = N.clamp(hLast - 0.3 * (hmax - hLast), 0, hminU * 0.98);
  let al, hf = hf0, m = 1.5, se = null, its = 0;
  try { const lm = N.levenbergMarquardt(p => un.map(q => Pmax * Math.pow(Math.max(1e-9, (q[0] - p[0]) / Math.max(1e-9, hmax - p[0])), p[1]) * p[2] - q[1]), [hf0, 1.5, 1], { maxIter: 100, bounds: [[0, hminU * 0.999], [1, 3], [0.5, 1.5]] }); hf = lm.p[0]; m = lm.p[1]; al = lm.p[2]; se = lm.stderr; its = lm.iterations; } catch { al = 1; }
  const S = al * Pmax * m / Math.max(1e-9, hmax - hf), hc = hmax - inp.op_eps * Pmax / S, hcOK = hc > 0, hcU = Math.max(hc, 0.05 * hmax), A = 24.5 * hcU * hcU + inp.op_C1 * hcU;   // mN/nm, nm, nm²
  const H = Pmax * 1e-3 / (A * 1e-18), Er = Math.sqrt(PI) / (2 * inp.op_beta) * (S * 1e6) / Math.sqrt(A * 1e-18), inv = 1 / Er - (1 - DIAMOND_TIP.nu ** 2) / DIAMOND_TIP.E, E = inv > 0 ? (1 - inp.op_nu ** 2) / inv : NaN;
  const fitP = q => al * Pmax * Math.pow(Math.max(0, (q - hf) / Math.max(1e-9, hmax - hf)), m), rm = N.rmse(un.map(q => fitP(q[0])), un.map(q => q[1])), seS = se && se.every(Number.isFinite) ? S * Math.hypot(se[1] / m, se[0] / Math.max(1e-9, hmax - hf), se[2] / al) : NaN, urel = fin(seS / S, rm / Pmax);
  let Wt = 0; for (let i = 1; i <= iU; i++) Wt += 0.5 * (P[i] + P[i - 1]) * (h[i] - h[i - 1]); let We = 0; for (let i = iU + 1; i < n; i++) We += 0.5 * (P[i] + P[i - 1]) * (h[i - 1] - h[i]);
  o.K('E', 'Sample elastic modulus', fin(E / 1e9), 'GPa', `From the reduced modulus with ν = ${inp.op_nu} and a diamond indenter (E = 1141 GPa, ν = 0.07): 1/E_r = (1 − ν²)/E + (1 − ν_i²)/E_i. Relative standard uncertainty from the unloading fit about ${F(100 * urel, 2)} %.`, Number.isFinite(E));
  o.K('H', 'Indentation hardness', H / 1e9, 'GPa', `Peak load ${F(Pmax, 4)} mN divided by the projected contact area A(h_c) = 24.5·h_c²${inp.op_C1 > 0 ? ` + ${inp.op_C1}·h_c` : ''} = ${F(A, 4)} nm² (ideal Berkovich geometry unless a tip-rounding term is given).`, hcOK);
  o.K('Er', 'Reduced modulus', Er / 1e9, 'GPa', `E_r = √π/(2β)·S/√A with β = ${inp.op_beta}; combines the elastic compliance of sample and indenter.`);
  o.K('S', 'Unloading contact stiffness', S * 1e6, 'N/m', `Slope dP/dh at peak load of the power-law fit P = α(h − h_f)^m to the upper 80 % of the unloading curve (${un.length} points, m = ${F(m, 4)}, h_f = ${F(hf, 4)} nm).`);
  o.K('hc', 'Contact depth', hc, 'nm', `h_c = h_max − ε·P_max/S with ε = ${inp.op_eps} and h_max = ${F(hmax, 4)} nm; sink-in is assumed, so pile-up makes the true contact area larger.`, hcOK);
  o.K('m', 'Unloading exponent m', m, '', 'Power-law exponent of the unloading curve; 1.2–1.6 is typical for a Berkovich tip. Values at the bounds (1 or 3) indicate a poor fit or creep during unloading.', m > 1.1 && m < 1.8);
  o.K('hf_ratio', 'Residual depth ratio h_f/h_max', hf / hmax, '', 'Below about 0.7 the Oliver–Pharr sink-in assumption holds; above 0.7 in a non-work-hardening material pile-up is likely and hardness and modulus are overestimated.', hf / hmax < 0.7);
  o.K('we', 'Elastic work fraction', Wt > 0 ? N.clamp(We / Wt, 0, 1) : 0, '', 'Recovered (unloading) work divided by total indentation work from trapezoidal integration; complements h_f/h_max as a plasticity indicator.');
  const hs = N.linspace(Math.max(hf, hminU * 0.9), hmax, 60);
  o.S('ph', 'Load–displacement curve with unloading fit', 'displacement (nm)', 'load (mN)', [['measured', h, P], ['power-law unloading fit', hs, hs.map(fitP)], ['stiffness tangent', [hmax - Pmax / S, hmax], [0, Pmax]]], 'The tangent at peak load is the contact stiffness; its intercept with the depth axis lies between h_f and h_c.');
  o.S('resid', 'Unloading fit residuals', 'displacement (nm)', 'residual (mN)', [['measured − fit', un.map(q => q[0]), un.map(q => q[1] - fitP(q[0]))]], 'Structure in the residuals (rather than scatter) reveals thermal drift, creep or an unsuitable fit range.');
  o.E('governing', 'Unloading power law', 'P = α·(h − h_f)^m,  S = dP/dh at h_max'); o.E('formula', 'Contact depth and area (Berkovich)', 'h_c = h_max − ε·P_max/S;  A = 24.5·h_c² + C₁·h_c'); o.E('formula', 'Hardness and reduced modulus', 'H = P_max/A;  E_r = (√π/(2β))·S/√A'); o.E('constitutive', 'Contact compliance', '1/E_r = (1 − ν²)/E + (1 − ν_i²)/E_i');
  o.C('Contact depth positive', hcOK, `h_c = ${F(hc, 4)} nm`); o.C('Sink-in assumption (h_f/h_max < 0.7)', hf / hmax < 0.7, `h_f/h_max = ${F(hf / hmax, 3)}`);
  if (inp.op_film > 0) o.C('10 % rule for coatings', hmax < 0.1 * inp.op_film, `h_max = ${F(hmax, 4)} nm against film thickness ${inp.op_film} nm (limit ${F(0.1 * inp.op_film, 4)} nm)`);
  o.C('Sample compliance physical', inv > 0, inv > 0 ? `E = ${F(E / 1e9, 4)} GPa` : 'reduced modulus exceeds the indenter modulus — check the area function or frame compliance');
  if (!hcOK) o.W('Computed contact depth is not positive; the unloading stiffness is too low for this ε — hardness and modulus are not meaningful.');
  if (inp.op_film > 0 && hmax >= 0.1 * inp.op_film) o.W(`Indentation depth exceeds 10 % of the ${inp.op_film} nm film: the substrate contributes to the measured hardness, and for modulus the influence starts even earlier.`);
  if (hf / hmax >= 0.7) o.W('h_f/h_max ≥ 0.7: pile-up is possible; verify the contact area by AFM/SEM imaging of the residual impression.');
  o.W('Frame compliance and thermal drift are assumed to be already corrected in the supplied curve; the ideal Berkovich area function underestimates area below roughly 100–200 nm contact depth unless C₁ is calibrated.');
  o.R('high', `Calibrate the tip area function on fused silica (E ≈ 72 GPa, H ≈ 9 GPa) at contact depths around ${F(hc, 3)} nm before quoting H = ${F(H / 1e9, 3)} GPa.`, 'Area-function error enters H linearly and E as its square root.');
  o.R('medium', inp.op_film > 0 ? `Keep h_max below ${F(0.1 * inp.op_film, 3)} nm (10 % rule) or use a film/substrate model; the present depth is ${F(hmax, 4)} nm.` : 'Use at least 10 indents and report mean and standard deviation; single-curve fit errors understate site-to-site variability.', 'Substrate influence and sampling statistics.');
  Object.assign(o.publish, Number.isFinite(E) ? { 'prop.E': E } : {}, { 'prop.hardness': H, 'cred.u_rel': urel, 'cred.rmse': rm / Pmax });
  o.numerics = { method: 'Oliver–Pharr: power-law unloading fit by Levenberg–Marquardt', dof: un.length, iterations: its, residual: rm, notes: 'upper 80 % of unloading used' };
  o.summary = `Oliver–Pharr analysis gives a contact stiffness of ${F(S * 1e6, 4)} N/m at h_max = ${F(hmax, 4)} nm, a contact depth of ${F(hc, 4)} nm, hardness ${F(H / 1e9, 3)} GPa and reduced modulus ${F(Er / 1e9, 4)} GPa, i.e. a sample modulus of ${F(fin(E / 1e9), 4)} GPa for ν = ${inp.op_nu}.`;
}

function mHertz(inp, o) {
  const t = o.tbl(inp.hz_table, 2, 8, 'Force–indentation curve'), rows = t.rows.sort((a, b) => a[0] - b[0]), d = rows.map(r => r[0]), Fn = rows.map(r => r[1]), n = d.length, Rt = inp.hz_R * 1e-9, dmax = d[n - 1], Fmax = N.maxOf(Fn);
  if (!(Fmax > 0) || !(dmax > d[0])) throw new Error('Force–indentation curve: increasing indentation with positive force is needed.');
  const span = dmax - d[0]; let k = Fmax / Math.pow(Math.max(1e-9, span), 1.5), d0 = d[0] + 0.1 * span, se = null, its = 0;
  try { const lm = N.levenbergMarquardt(p => d.map((x, i) => p[0] * k * Math.pow(Math.max(0, x - p[1]), 1.5) - Fn[i]), [1, d0], { maxIter: 100, bounds: [[1e-6, 1e6], [d[0] - span, dmax - 0.05 * span]] }); k *= lm.p[0]; d0 = lm.p[1]; se = lm.stderr ? [lm.stderr[0] / lm.p[0], lm.stderr[1]] : null; its = lm.iterations; } catch { /* initial estimate kept */ }
  const kSI = k * 1e-9 / Math.pow(1e-9, 1.5), Er = 3 * kSI / (4 * Math.sqrt(Rt)), E = Er * (1 - inp.hz_nu ** 2), dc = Math.max(0, dmax - d0) * 1e-9, a = Math.sqrt(Rt * dc), p0 = a > 0 ? 1.5 * Fmax * 1e-9 / (PI * a * a) : 0, fit = d.map(x => k * Math.pow(Math.max(0, x - d0), 1.5));
  const cz = d.map((x, i) => i).filter(i => d[i] - d0 > 0.1 * (dmax - d0) && Fn[i] > 0), fl = cz.length >= 4 ? linfit(cz.map(i => Math.log(d[i] - d0)), cz.map(i => Math.log(Fn[i]))) : null, urel = fin(se?.[0], N.rmse(fit, Fn) / Fmax);
  o.K('E', 'Sample elastic modulus', E, 'Pa', `E = E_r·(1 − ν²) with ν = ${inp.hz_nu} for a rigid spherical tip of radius ${inp.hz_R} nm; reduced modulus E_r = 3k/(4√R) = ${F(Er, 4)} Pa from the fitted Hertz prefactor (± ${F(100 * urel, 2)} %).`, true);
  o.K('d0', 'Contact point', d0, 'nm', `Fitted offset of the indentation axis (± ${F(fin(se?.[1]), 2)} nm); an error of the contact point is the main source of modulus error in soft-sample AFM indentation.`);
  o.K('a', 'Contact radius at maximum load', a * 1e9, 'nm', 'a = √(R·δ) from Hertz theory at the deepest point; sets the lateral resolution of the measurement.');
  o.K('p0', 'Peak contact pressure', p0, 'Pa', 'p₀ = 3F/(2πa²), the maximum Hertz pressure under the centre of the contact; compare with the yield strength to judge whether the contact stayed elastic.');
  o.K('ratio', 'Indentation-to-radius ratio δ/R', dc / Rt, '', 'Hertz theory is a small-strain model: keep δ/R below about 0.1 (a/R below 0.3); beyond that use a Sneddon or finite-strain model.', dc / Rt < 0.1);
  o.K('n_exp', 'Free power-law exponent', fl ? fl.slope : 1.5, '', 'Slope of log F against log(δ − δ₀); 1.5 is expected for a sphere, 2 for a cone/pyramid, 1 for a flat punch. A departure indicates the wrong tip model, adhesion or substrate stiffening.', fl ? Math.abs(fl.slope - 1.5) < 0.15 : null);
  o.K('rmse', 'Fit residual (RMSE)', N.rmse(fit, Fn), 'nN', 'Root-mean-square deviation between the measured force and the Hertz fit over the whole curve.');
  o.S('fd', 'Force–indentation curve with Hertz fit', 'indentation (nm)', 'force (nN)', [['measured', d, Fn], ['Hertz fit', d, fit]], 'The fit includes the non-contact baseline before the contact point; systematic deviation at large depth indicates substrate effect or non-linear elasticity.');
  o.S('lin', 'Linearised Hertz plot', 'indentation − δ₀ (nm)', 'F^(2/3) (nN^(2/3))', [['measured', d.map(x => x - d0), Fn.map(v => Math.pow(Math.max(0, v), 2 / 3))], ['Hertz', d.map(x => x - d0), fit.map(v => Math.pow(v, 2 / 3))]], 'F^(2/3) is linear in indentation for a Hertzian contact; curvature exposes adhesion or a wrong contact point.');
  o.E('governing', 'Hertz contact of a sphere on a half-space', 'F = (4/3)·E_r·√R·δ^(3/2);  a = √(R·δ);  p₀ = 3F/(2πa²)'); o.E('constitutive', 'Reduced modulus (rigid tip)', 'E_r = E/(1 − ν²)');
  o.C('Small-strain validity δ/R < 0.1', dc / Rt < 0.1, `δ/R = ${F(dc / Rt, 3)}`); o.C('Exponent consistent with a sphere', fl ? Math.abs(fl.slope - 1.5) < 0.15 : false, fl ? `n = ${F(fl.slope, 3)}` : 'too few contact points');
  if (dc / Rt >= 0.1) o.W('Indentation exceeds 10 % of the tip radius: the Hertz small-strain assumption is violated and the modulus is biased.');
  o.W('Adhesion (JKR/DMT), viscoelasticity, finite sample thickness and tip-shape error are neglected; the tip is treated as rigid.');
  o.R('high', `Verify the tip radius (${inp.hz_R} nm) by imaging a tip characteriser: E scales with R^(−1/2), so a 20 % radius error gives a 10 % modulus error on ${F(E, 3)} Pa.`, 'Hertz prefactor sensitivity.');
  o.R('medium', 'Acquire a force map (≥ 64 curves) and report the modulus distribution; a single curve does not represent a heterogeneous surface.', 'Site-to-site variability.');
  Object.assign(o.publish, { 'prop.E': E, 'cred.u_rel': urel, 'cred.rmse': N.rmse(fit, Fn) / Fmax });
  o.numerics = { method: 'Hertz fit with free contact point (Levenberg–Marquardt)', dof: n, iterations: its, residual: N.rmse(fit, Fn) };
  o.summary = `The Hertz fit gives a reduced modulus of ${F(Er, 4)} Pa and a sample modulus of ${F(E, 4)} Pa for ν = ${inp.hz_nu}, with the contact point at ${F(d0, 3)} nm. At maximum load the contact radius is ${F(a * 1e9, 3)} nm and δ/R = ${F(dc / Rt, 3)}.`;
}

function mFPP(inp, o) {
  const t = o.tbl(inp.fpp_table, 2, 3, 'Four-point-probe I–V data'), I = t.rows.map(r => r[0] * 1e-3), V = t.rows.map(r => r[1] * 1e-3), f = linfit(I, V), Rm = f.slope, th = inp.el_thick * 1e-9, s = inp.fpp_s * 1e-3, ds = Math.max(2, inp.fpp_d / inp.fpp_s), ts = th / s;
  const F2 = LN2 / (LN2 + Math.log((ds * ds + 3) / (ds * ds - 3))), F1 = ts < 1e-3 ? 1 : LN2 / Math.log(Math.sinh(Math.min(ts, 300)) / Math.sinh(Math.min(ts, 300) / 2)), Rs = PI / LN2 * Math.abs(Rm) * F2, rho = Rs * th * F1, sig = 1 / rho, urel = Rm !== 0 ? Math.abs(f.seSlope / Rm) : 0;
  o.K('Rs', 'Sheet resistance', Rs, 'Ω/sq', `R_s = (π/ln 2)·(V/I)·F₂ with V/I = ${F(Rm, 5)} Ω from the slope of ${I.length} I–V points and the lateral correction F₂ = ${F(F2, 4)} for a circular sample ${F(ds, 3)} probe spacings across. Valid for a film much thinner than the probe spacing.`, true);
  o.K('rho', 'Resistivity', rho, 'Ω·m', `ρ = R_s·t·F₁ with t = ${inp.el_thick} nm and thickness correction F₁ = ${F(F1, 4)} (t/s = ${F(ts, 2)}); thickness uncertainty enters linearly and is not included in the fit error.`);
  o.K('sigma', 'Electrical conductivity', sig, 'S/m', 'Reciprocal of the resistivity; compare with the bulk value of the material to judge grain-boundary and surface scattering in thin films.');
  o.K('u', 'Relative standard uncertainty (fit)', urel, '', 'Standard error of the I–V slope divided by the slope; covers electrical noise only. Add probe-spacing (≈ 1 %), thickness and temperature contributions in the uncertainty-budget mode.', urel < 0.01);
  o.K('offset', 'Voltage offset at zero current', f.intercept * 1e3, 'mV', 'Intercept of the I–V line; thermoelectric voltages and amplifier offsets appear here and are removed by using the slope rather than a single-point ratio.');
  o.K('r2', 'Ohmic linearity R²', f.r2, '', 'Coefficient of determination of the I–V line; values below 0.999 point to non-ohmic contacts, self-heating or injection effects.', f.r2 > 0.999);
  o.K('F', 'Geometric correction F₁·F₂', F1 * F2, '', 'Product of the finite-thickness and finite-diameter correction factors applied to the infinite-sheet formula (collinear equally spaced probes at the centre of a circular sample).');
  o.S('iv', 'Current–voltage line', 'current (mA)', 'voltage (mV)', [['measured', I.map(v => v * 1e3), V.map(v => v * 1e3)], ['linear fit', [N.minOf(I) * 1e3, N.maxOf(I) * 1e3], [(f.intercept + f.slope * N.minOf(I)) * 1e3, (f.intercept + f.slope * N.maxOf(I)) * 1e3]]], 'A straight line through both polarities confirms ohmic behaviour; its slope is the measured resistance V/I.');
  o.S('res', 'Residuals of the linear fit', 'current (mA)', 'residual (µV)', [['residual', I.map(v => v * 1e3), V.map((v, i) => (v - f.intercept - f.slope * I[i]) * 1e6)]], 'A systematic S-shape reveals heating or non-ohmic contacts; random scatter is electrical noise.');
  o.E('governing', 'Collinear four-point probe on a thin sheet', 'R_s = (π/ln 2)·(V/I)·F₂,  ρ = R_s·t·F₁'); o.E('formula', 'Finite-size corrections', 'F₁ = ln 2 / ln[sinh(t/s)/sinh(t/2s)];  F₂ = ln 2 / {ln 2 + ln[((d/s)² + 3)/((d/s)² − 3)]}'); o.E('constitutive', "Ohm's law", 'J = σ·E,  σ = 1/ρ');
  o.C('Ohmic response', f.r2 > 0.999, `R² = ${F(f.r2, 6)}`); o.C('Sample much larger than probe spacing', ds > 10, `d/s = ${F(ds, 3)} (correction ${F(100 * (1 - F2), 2)} %)`); o.C('Positive resistance', Rm > 0, `V/I = ${F(Rm, 5)} Ω`);
  if (inp.fpp_d / inp.fpp_s < 2) o.W('Sample diameter below two probe spacings: the correction formula is not valid and d/s was limited to 2.');
  if (ts > 0.4) o.W('Thickness exceeds 40 % of the probe spacing: the thin-sheet assumption is weak and the thickness correction F₁ is significant.');
  o.R('high', `Measure the film thickness independently (profilometry or ellipsometry): ρ = ${F(rho, 3)} Ω·m inherits its relative error one-to-one.`, 'ρ = R_s·t.');
  o.R('medium', f.r2 > 0.999 ? 'Reverse current polarity and average, as done here with a bipolar sweep, in every routine measurement to cancel thermoelectric offsets.' : 'Lower the current or improve the contacts: the I–V line is not linear enough for a resistance to be defined.', 'Linearity of the I–V sweep.');
  Object.assign(o.publish, { 'prop.sigma_elec': sig, 'elec.sheet_resistance': Rs, 'cred.u_rel': urel, 'cred.rmse': f.s });
  o.numerics = { method: 'Linear regression of V on I with Smits-type geometric corrections', dof: I.length, residual: f.s };
  o.summary = `The I–V slope of ${F(Rm, 5)} Ω gives a sheet resistance of ${F(Rs, 4)} Ω/sq and, with a thickness of ${inp.el_thick} nm, a resistivity of ${F(rho, 3)} Ω·m (conductivity ${F(sig, 3)} S/m). The fit contributes ${F(100 * urel, 2)} % relative standard uncertainty.`;
}

function mHall(inp, o) {
  const t = o.tbl(inp.hall_table, 2, 3, 'Hall data'), B = t.rows.map(r => r[0]), V = t.rows.map(r => r[1] * 1e-3), f = linfit(B, V), th = inp.el_thick * 1e-9, I = inp.hall_I * 1e-3, RH = f.slope * th / I, RA = inp.hall_RA, RB = inp.hall_RB;
  const Rs0 = PI * (RA + RB) / (2 * LN2); let Rs = Rs0; try { Rs = N.brent(r => Math.exp(-PI * RA / r) + Math.exp(-PI * RB / r) - 1, Rs0 * 1e-3, Rs0 * 1.0001, 1e-14 * Rs0); } catch { /* symmetric estimate kept */ }
  const rho = Rs * th, nC = RH !== 0 ? 1 / (N.Q_E * Math.abs(RH)) : NaN, mu = Math.abs(RH) / rho, type = RH < 0 ? 'n-type (electrons)' : RH > 0 ? 'p-type (holes)' : 'undetermined', urel = f.slope !== 0 ? Math.abs(f.seSlope / f.slope) : 1, sig = Math.abs(f.slope) > 2 * f.seSlope;
  o.K('type', 'Majority carrier type', sig ? type : 'undetermined', '', `Sign of the Hall coefficient R_H = ${F(RH, 3)} m³/C from the slope of the Hall voltage against field (${B.length} points). ${sig ? '' : 'The slope is within two standard errors of zero, so the sign is not significant.'}`, sig);
  o.K('n', 'Carrier density', fin(nC), '1/m³', `n = 1/(q·|R_H|) assuming a single carrier type and Hall scattering factor r_H = 1 (= ${F(fin(nC) * 1e-6, 3)} cm⁻³); relative standard uncertainty from the slope ${F(100 * urel, 2)} %.`, sig);
  o.K('mu', 'Hall mobility', mu, 'm²/V/s', `µ_H = |R_H|/ρ = ${F(mu * 1e4, 4)} cm²/V/s using the van der Pauw resistivity; it differs from the drift mobility by the scattering factor (1–2).`);
  o.K('Rs', 'Sheet resistance (van der Pauw)', Rs, 'Ω/sq', `Numerical solution of exp(−πR_A/R_s) + exp(−πR_B/R_s) = 1 for R_A = ${RA} Ω and R_B = ${RB} Ω (Brent's method); correction factor f = ${F(Rs / Rs0, 4)} relative to the symmetric formula.`);
  o.K('rho', 'Resistivity', rho, 'Ω·m', `ρ = R_s·t with t = ${inp.el_thick} nm; requires a uniform, hole-free film with small contacts on the perimeter.`);
  o.K('RH', 'Hall coefficient', RH, 'm³/C', 'R_H = (dV_H/dB)·t/I; negative for electrons, positive for holes. Mixed conduction makes it field-dependent and smaller in magnitude.');
  o.K('offset', 'Misalignment offset voltage', f.intercept * 1e3, 'mV', 'Zero-field intercept of the Hall voltage caused by contact misalignment; removed by fitting the slope over both field polarities.');
  o.K('sym', 'van der Pauw asymmetry R_A/R_B', RB > 0 ? RA / RB : 0, '', 'Ratio of the two four-terminal resistances; values far from 1 indicate an anisotropic or non-uniform sample or unequal contact placement.', RB > 0 ? RA / RB > 0.5 && RA / RB < 2 : false);
  o.S('vh', 'Hall voltage against magnetic field', 'B (T)', 'V_H (mV)', [['measured', B, V.map(v => v * 1e3)], ['linear fit', [N.minOf(B), N.maxOf(B)], [(f.intercept + f.slope * N.minOf(B)) * 1e3, (f.intercept + f.slope * N.maxOf(B)) * 1e3]]], 'The slope gives the Hall coefficient; curvature would indicate two-carrier conduction or magnetoresistive pickup.');
  const rr = N.logspace(Rs0 * 0.2, Rs0 * 2, 60);
  o.S('vdp', 'van der Pauw equation residual', 'trial R_s (Ω/sq)', 'exp(−πR_A/R_s) + exp(−πR_B/R_s) − 1', [['residual', rr, rr.map(r => Math.exp(-PI * RA / r) + Math.exp(-PI * RB / r) - 1)]], 'The sheet resistance is the zero crossing of this monotonic function.', { logx: true });
  o.E('governing', 'Hall effect (single carrier)', 'V_H = R_H·I·B/t,  R_H = −1/(q·n) (electrons) or +1/(q·p) (holes)'); o.E('governing', 'van der Pauw equation', 'exp(−π·R_A/R_s) + exp(−π·R_B/R_s) = 1'); o.E('constitutive', 'Drude conductivity', 'σ = q·n·µ,  µ_H = |R_H|/ρ');
  o.C('Hall slope significant', sig, `slope ${F(f.slope, 3)} ± ${F(f.seSlope, 2)} V/T`); o.C('Hall voltage linear in field', f.r2 > 0.99, `R² = ${F(f.r2, 5)}`); o.C('van der Pauw solution consistent', Math.abs(Math.exp(-PI * RA / Rs) + Math.exp(-PI * RB / Rs) - 1) < 1e-6, `residual ${F(Math.exp(-PI * RA / Rs) + Math.exp(-PI * RB / Rs) - 1, 2)}`);
  if (!sig) o.W('The Hall slope is not statistically different from zero; carrier type and density are not determined (compensated material, too low field or too thick a sample).');
  o.W('Single-carrier, isotropic model with Hall factor 1; degenerate or mixed-conduction films, grain-boundary barriers and finite contact size bias density and mobility.');
  o.R('high', `Repeat with both current polarities and all contact permutations; with ${F(100 * urel, 2)} % slope uncertainty the carrier density ${F(fin(nC), 3)} m⁻³ is limited by ${urel > 0.02 ? 'electrical noise' : 'thickness and geometry rather than noise'}.`, 'Slope standard error against type-B contributions.');
  o.R('medium', 'Keep contacts below about 5 % of the sample edge length and on the perimeter; finite contacts inside the sample area cause errors of several percent in van der Pauw resistivity.', 'van der Pauw contact-size condition.');
  Object.assign(o.publish, { 'prop.sigma_elec': 1 / rho, 'elec.sheet_resistance': Rs, 'elec.mobility': mu, 'cred.u_rel': urel, 'cred.rmse': f.s }, Number.isFinite(nC) ? { 'elec.carrier_density': nC } : {});
  o.numerics = { method: 'Linear regression of Hall voltage; Brent root of the van der Pauw equation', dof: B.length, residual: f.s };
  o.summary = `The Hall coefficient of ${F(RH, 3)} m³/C indicates ${sig ? type : 'no significant carrier sign'} conduction with a carrier density of ${F(fin(nC), 3)} m⁻³. With the van der Pauw sheet resistance of ${F(Rs, 4)} Ω/sq the resistivity is ${F(rho, 3)} Ω·m and the Hall mobility ${F(mu * 1e4, 4)} cm²/V/s.`;
}

function mEIS(inp, o) {
  const t = o.tbl(inp.eis_table, 3, 8, 'Impedance spectrum'), rows = t.rows.filter(r => r[0] > 0).sort((a, b) => b[0] - a[0]), f = rows.map(r => r[0]), zr = rows.map(r => r[1]), zi = rows.map(r => r[2]), n = f.length;
  if (n < 8) throw new Error('Impedance spectrum: at least 8 rows with positive frequency are needed.');
  const mod = zr.map((v, i) => Math.hypot(v, zi[i]) || 1), Rs0 = Math.max(1e-6, N.minOf(zr)), span = Math.max(1e-6, N.maxOf(zr) - Rs0), wLow = 2 * PI * f[n - 1];
  const resid = (model, p) => { const r = new Array(2 * n); for (let i = 0; i < n; i++) { const z = model(f[i], p); r[2 * i] = fin((z[0] - zr[i]) / mod[i], 1e6); r[2 * i + 1] = fin((z[1] - zi[i]) / mod[i], 1e6); } return r; };
  const full = (ff, p) => randles(ff, Math.exp(p[0]), Math.exp(p[1]), Math.exp(p[2]), p[3], Math.exp(p[4])), simple = (ff, p) => randles(ff, Math.exp(p[0]), Math.exp(p[1]), Math.exp(p[2]), 1, 0);
  let best = null, bestS = null, its = 0;
  for (const fr of [0.3, 0.6, 0.9]) {
    const Rct = fr * span, sub = f.map((_, i) => i).filter(i => zr[i] < Rs0 + Rct), ipk = sub.length ? sub.reduce((a, b) => (-zi[b] > -zi[a] ? b : a)) : 0, C0 = 1 / (2 * PI * f[ipk] * Rct), sw0 = Math.max(1e-6, Math.abs(zi[n - 1]) * Math.sqrt(wLow) * 0.5);
    try { const lm = N.levenbergMarquardt(p => resid(full, p), [Math.log(Rs0), Math.log(Rct), Math.log(C0), 0.9, Math.log(sw0)], { maxIter: 80, bounds: [[-14, 20], [-14, 20], [-35, 5], [0.5, 1], [-25, 15]] }); its += lm.iterations; if (!best || lm.cost < best.cost) best = lm; } catch { /* next start */ }
    try { const lm = N.levenbergMarquardt(p => resid(simple, p), [Math.log(Rs0), Math.log(Rct), Math.log(C0)], { maxIter: 60, bounds: [[-14, 20], [-14, 20], [-35, 5]] }); if (!bestS || lm.cost < bestS.cost) bestS = lm; } catch { /* next start */ }
  }
  if (!best) throw new Error('Impedance spectrum: the equivalent-circuit fit did not converge; check the sign convention (capacitive Z_im negative).');
  const p = best.p, Rs = Math.exp(p[0]), Rct = Math.exp(p[1]), Q = Math.exp(p[2]), al = p[3], sw = Math.exp(p[4]), se = best.stderr && best.stderr.every(Number.isFinite) ? best.stderr : null, rel = i => (se ? Math.min(se[i], 10) : NaN);
  const Ceff = Math.pow(Q, 1 / al) * Math.pow(1 / Rs + 1 / Rct, (al - 1) / al), aic = (c, k) => 2 * n * Math.log(Math.max(1e-300, c) / (2 * n)) + 2 * k, aF = aic(best.cost, 5), aS = bestS ? aic(bestS.cost, 3) : Infinity, chi = Math.sqrt(best.cost / (2 * n)), wSig = sw / Math.sqrt(wLow) > 0.03 * Rct;
  const zf = f.map(ff => full(ff, p)), zs = bestS ? f.map(ff => simple(ff, bestS.p)) : null;
  o.K('Rs', 'Series (solution) resistance R_s', Rs, 'Ω', `High-frequency intercept of the Nyquist arc from the complex non-linear least-squares fit; relative standard error ${F(100 * fin(rel(0)), 2)} %.`, true);
  o.K('Rct', 'Charge-transfer resistance R_ct', Rct, 'Ω', `Diameter of the semicircle; inversely proportional to the exchange current (R_ct = RT/(nF·i₀) at small overpotential). Relative standard error ${F(100 * fin(rel(1)), 2)} %.`, fin(rel(1), 1) < 0.1);
  o.K('Cdl', 'Effective double-layer capacitance', Ceff, 'F', `Brug conversion of the constant-phase element, C = Q^(1/α)·(1/R_s + 1/R_ct)^((α−1)/α), with Q = ${F(Q, 3)} S·s^α; divide by the electrode area for a specific capacitance.`);
  o.K('alpha', 'CPE exponent α', al, '', `1 is an ideal capacitor; 0.8–0.95 is usual for rough or porous electrodes; 0.5 would be diffusion-like. Standard error ${F(fin(se?.[3]), 2)}.`, al > 0.75);
  o.K('sigma_w', 'Warburg coefficient σ', sw, 'Ω·s^(-1/2)', `Semi-infinite diffusion element Z_W = σ(1 − j)/√ω; relative standard error ${F(100 * fin(rel(4)), 2)} %. ${wSig ? '' : 'Its contribution at the lowest frequency is below 3 % of R_ct, so it is not identifiable from this spectrum.'}`, wSig && fin(rel(4), 1) < 0.3);
  o.K('tau', 'Interfacial time constant R_ct·C', Rct * Ceff, 's', `Characteristic relaxation time of the charge-transfer arc; apex frequency ${F(1 / (2 * PI * Rct * Ceff), 3)} Hz.`);
  o.K('chi', 'Fit residual (modulus-weighted RMS)', chi, '', 'Root-mean-square relative deviation of real and imaginary parts; about the noise level (≤ 1–2 %) for an adequate circuit.', chi < 0.03);
  o.K('daic', 'ΔAIC (ideal R–C circuit − full Randles)', Number.isFinite(aS) ? aS - aF : 0, '', 'Akaike information criterion difference between the 3-parameter R_s + R_ct‖C circuit and the 5-parameter CPE + Warburg circuit; above about 10 the extra elements are justified by the data.', aS - aF > 10);
  o.S('nyq', 'Nyquist plot', 'Z′ (Ω)', '−Z″ (Ω)', [['measured', zr, zi.map(v => -v)], ['Randles (CPE + Warburg) fit', zf.map(z => z[0]), zf.map(z => -z[1])], ...(zs ? [['ideal R–C circuit', zs.map(z => z[0]), zs.map(z => -z[1])]] : [])], 'Semicircle from charge transfer in parallel with the double layer, followed by the 45° Warburg line at low frequency.');
  o.S('bode', 'Bode magnitude', 'frequency (Hz)', '|Z| (Ω)', [['measured', f, mod], ['fit', f, zf.map(z => Math.hypot(z[0], z[1]))]], 'Plateaus at R_s (high frequency) and R_s + R_ct (intermediate) with the diffusion rise at low frequency.', { logx: true, logy: true });
  o.T('par', 'Equivalent-circuit parameters', ['parameter', 'value', 'unit', 'relative standard error (%)'], [['R_s', Rs, 'Ω', 100 * rel(0)], ['R_ct', Rct, 'Ω', 100 * rel(1)], ['Q (CPE)', Q, 'S·s^α', 100 * rel(2)], ['α', al, '—', se ? 100 * se[3] / al : NaN], ['σ (Warburg)', sw, 'Ω·s^-½', 100 * rel(4)]], 'Standard errors come from the Jacobian of the modulus-weighted fit at the optimum (linearised covariance); they assume the circuit is the right model.');
  o.E('governing', 'Randles circuit impedance', 'Z(ω) = R_s + 1/[Q·(jω)^α + 1/(R_ct + σ(1 − j)/√ω)]'); o.E('governing', 'Complex non-linear least squares', 'min Σᵢ [(Z′ᵢ − Z′(ωᵢ))² + (Z″ᵢ − Z″(ωᵢ))²]/|Zᵢ|²'); o.E('formula', 'Effective capacitance of a CPE (Brug)', 'C = Q^(1/α)·(1/R_s + 1/R_ct)^((α − 1)/α)'); o.E('formula', 'Model selection', 'AIC = N·ln(RSS/N) + 2k');
  o.C('Fit residual at noise level', chi < 0.03, `χ = ${F(chi, 3)}`); o.C('Passive circuit (all elements positive)', Rs > 0 && Rct > 0 && Q > 0, 'R_s, R_ct, Q, σ > 0 enforced by log-parameterisation'); o.C('Warburg element identifiable', wSig, `σ/√ω_min = ${F(sw / Math.sqrt(wLow), 3)} Ω against R_ct = ${F(Rct, 3)} Ω`);
  if (!se) o.W('Parameter covariance could not be computed (singular Jacobian): at least one circuit element is not identifiable from this spectrum.');
  o.W('Equivalent circuits are not unique: different circuits can fit the same spectrum. Linearity, stability and causality (Kramers–Kronig) of the measurement are assumed, not tested.');
  o.R('high', wSig ? `Extend the spectrum below ${F(f[n - 1], 2)} Hz only if the diffusion coefficient is needed; R_ct = ${F(Rct, 4)} Ω is already determined to ${F(100 * fin(rel(1)), 2)} %.` : `Measure to lower frequency (below ${F(f[n - 1], 2)} Hz): the Warburg branch is not resolved, so diffusion parameters cannot be reported.`, 'Identifiability of the Warburg coefficient.');
  o.R('medium', `Normalise C_dl = ${F(Ceff, 3)} F by the geometric electrode area and compare with 10–40 µF/cm² for a smooth surface to estimate the electrochemically active area.`, 'Double-layer capacitance per unit area.');
  Object.assign(o.publish, { 'elec.R_ct': Rct, 'elec.R_s': Rs, 'elec.C_dl': Ceff, 'cred.rmse': chi, 'cred.u_rel': fin(rel(1), chi) });
  o.numerics = { method: 'Complex non-linear least squares (Levenberg–Marquardt, modulus weighting, log-parameters, three starts)', dof: 2 * n - 5, iterations: its, residual: chi };
  o.summary = `The Randles fit gives R_s = ${F(Rs, 4)} Ω, R_ct = ${F(Rct, 4)} Ω, an effective double-layer capacitance of ${F(Ceff, 3)} F (CPE exponent ${F(al, 3)}) and a Warburg coefficient of ${F(sw, 3)} Ω·s^-½, with a modulus-weighted residual of ${F(100 * chi, 2)} %. ΔAIC against an ideal R–C circuit is ${F(Number.isFinite(aS) ? aS - aF : 0, 3)}.`;
}

function mFlash(inp, o) {
  const t = o.tbl(inp.flash_table, 2, 10, 'Laser-flash signal'), rows = t.rows.sort((a, b) => a[0] - b[0]), tm = rows.map(r => r[0] * 1e-3), s = rows.map(r => r[1]), n = tm.length, ss = smooth(s, 2), base = N.mean(s.slice(0, Math.max(2, Math.round(n * 0.03)))), im = N.argmax(ss), top = ss[im], amp = top - base, L = inp.fl_L * 1e-3;
  if (!(amp > 0)) throw new Error('Laser-flash signal: the detector signal must rise above its initial baseline.');
  const cross = fr => { for (let i = 1; i <= im; i++) if ((ss[i] - base) / amp >= fr && (ss[i - 1] - base) / amp < fr) return tm[i - 1] + (fr - (ss[i - 1] - base) / amp) / ((ss[i] - ss[i - 1]) / amp) * (tm[i] - tm[i - 1]); return NaN; };
  const th = cross(0.5), fr = [0.25, 0.5, 0.75], al = fr.map(q => parkerK(q) * L * L / cross(q)), a = al[1], spread = N.std(al.filter(Number.isFinite)) / a, decay = (top - ss[n - 1]) / amp, k = a * inp.rho * inp.fl_cp, ok = Number.isFinite(a) && a > 0;
  const model = tm.map(x => base + amp * parker(PI * PI * fin(a) * x / (L * L)));
  o.K('alpha', 'Thermal diffusivity', fin(a), 'm²/s', `Parker half-time method α = 0.1388·L²/t½ with L = ${inp.fl_L} mm and t½ = ${F(fin(th) * 1e3, 4)} ms read from the baseline-corrected, smoothed rear-face signal. Adiabatic, instantaneous-pulse model.`, ok);
  o.K('k', 'Thermal conductivity', fin(k), 'W/m/K', `k = α·ρ·c_p with ρ = ${inp.rho} kg/m³ and c_p = ${inp.fl_cp} J/kg/K; density and heat capacity errors add in quadrature to the diffusivity error.`);
  o.K('t_half', 'Half-rise time', fin(th) * 1e3, 'ms', 'Time for the rear face to reach half of its maximum temperature rise; should be at least ~10 times the laser pulse width for the instantaneous-pulse assumption.');
  o.K('spread', 'Consistency of 25/50/75 % rise times', fin(spread), '', `Relative standard deviation of the diffusivity evaluated from the 25 %, 50 % and 75 % rise times with their own Parker constants (${fr.map(q => F(parkerK(q), 4)).join(', ')}); above about 5 % indicates heat loss, finite pulse or non-uniform heating.`, spread < 0.05);
  o.K('decay', 'Signal decay after the maximum', decay, '', 'Fraction of the rise lost by the end of the record; a clear decay means radiative/convective heat loss, for which a Cowan or Cape–Lehman correction should replace the adiabatic Parker formula.', decay < 0.05);
  o.K('effus', 'Thermal effusivity', ok ? Math.sqrt(k * inp.rho * inp.fl_cp) : 0, 'W·s^½/m²/K', 'e = √(k·ρ·c_p): governs contact temperature and transient heat exchange with adjoining materials.');
  o.K('rmse', 'Model residual', N.rmse(model, s) / amp, '', 'RMS difference between the measured curve and the adiabatic Parker curve for the fitted diffusivity, relative to the signal amplitude.', N.rmse(model, s) / amp < 0.05);
  o.S('rise', 'Rear-face temperature rise', 'time (ms)', 'detector signal', [['measured', tm.map(v => v * 1e3), s], ['Parker model', tm.map(v => v * 1e3), model]], 'Measured transient against the adiabatic one-dimensional solution; a measured curve falling below the model after the maximum shows heat loss.');
  o.S('norm', 'Normalised rise against dimensionless time', 'π²·α·t/L²', '(V − V₀)/(V_max − V₀)', [['measured', tm.map(x => PI * PI * fin(a) * x / (L * L)), s.map(v => (v - base) / amp)], ['Parker', N.linspace(0, 10, 100), N.linspace(0, 10, 100).map(parker)]], 'In dimensionless form every adiabatic measurement collapses on one master curve; the half-rise occurs at 1.37.');
  o.E('governing', 'One-dimensional heat conduction', '∂T/∂t = α·∂²T/∂x²,  adiabatic faces, instantaneous surface pulse'); o.E('formula', 'Parker rear-face solution', 'V(t) = 1 + 2·Σₙ (−1)ⁿ·exp(−n²π²αt/L²);  α = 0.1388·L²/t½'); o.E('constitutive', 'Conductivity', 'k = α·ρ·c_p');
  o.C('Half-rise time found', ok, `t½ = ${F(fin(th) * 1e3, 4)} ms`); o.C('Adiabatic (no heat-loss) assumption', decay < 0.05, `signal decay ${F(100 * decay, 2)} %`); o.C('Rise-time consistency', spread < 0.05, `spread ${F(100 * fin(spread), 2)} %`);
  if (decay >= 0.05) o.W('The signal decays after its maximum: apply a heat-loss correction (Cowan) — the Parker value overestimates diffusivity.');
  o.W('Assumes a homogeneous, opaque, isotropic disc; porous, translucent or layered samples need a coating and a multilayer or radiation model.');
  o.R('high', `Measure the thickness to ±1 µm: α ∝ L², so a 1 % thickness error gives 2 % on α = ${F(fin(a), 3)} m²/s.`, 'Quadratic thickness sensitivity.');
  o.R('medium', 'Determine c_p by DSC on the same batch and density by Archimedes; literature values of dense material overestimate k of porous nanostructured samples.', 'k = α·ρ·c_p error propagation.');
  Object.assign(o.publish, ok ? { 'prop.k_thermal': k, 'prop.alpha_thermal': a, 'cred.u_rel': Math.max(fin(spread), 0.02) } : {}, { 'cred.rmse': N.rmse(model, s) / amp });
  o.numerics = { method: 'Parker half-time reduction with constants from the series solution (bisection)', dof: n, residual: N.rmse(model, s) / amp };
  o.summary = `The half-rise time of ${F(fin(th) * 1e3, 4)} ms through ${inp.fl_L} mm gives a thermal diffusivity of ${F(fin(a), 3)} m²/s and a conductivity of ${F(fin(k), 4)} W/m/K with the stated density and heat capacity. Rise-time consistency is ${F(100 * fin(spread), 2)} % and the post-maximum decay ${F(100 * decay, 2)} %.`;
}

function mTensile(inp, o) {
  const t = o.tbl(inp.tens_table, 2, 10, 'Stress–strain curve'), e = t.rows.map(r => r[0]), s = t.rows.map(r => r[1]), n = e.length, hasT = t.rows.every(r => Number.isFinite(r[2])), et = hasT ? t.rows.map(r => r[2]) : null, iu = N.argmax(s), uts = s[iu], off = inp.tens_offset / 100;
  if (!(uts > 0)) throw new Error('Stress–strain curve: positive stresses are needed.');
  let best = null;                                                            // elastic segment: 5 % to an upper stress fraction of the maximum, shrunk until linear
  for (const up of [0.4, 0.33, 0.27, 0.22, 0.18, 0.15]) { const id = []; for (let i = 0; i <= iu; i++) if (s[i] >= 0.05 * uts && s[i] <= up * uts) id.push(i); if (id.length < 4) continue; const f = linfit(id.map(i => e[i]), id.map(i => s[i])); if (!best || f.r2 > best.f.r2) best = { i: id[0], j: id[id.length - 1], f, up, m: id.length }; if (f.r2 > 0.9995) break; }
  if (!best || !(best.f.slope > 0)) { const w0 = Math.max(4, Math.round(n * 0.15)); best = { i: 0, j: w0 - 1, f: linfit(e.slice(0, w0), s.slice(0, w0)), up: 1, m: w0 }; }
  const w = best.j - best.i + 1;
  const E = best.f.slope, e0 = -best.f.intercept / (E || 1); let sy = NaN, ey = NaN;
  for (let i = best.j + 1; i < n; i++) { const g0 = s[i - 1] - E * (e[i - 1] - e0 - off), g1 = s[i] - E * (e[i] - e0 - off); if (g0 > 0 && g1 <= 0) { const q = g0 / (g0 - g1); sy = s[i - 1] + q * (s[i] - s[i - 1]); ey = e[i - 1] + q * (e[i] - e[i - 1]); break; } }
  const yielded = Number.isFinite(sy), tough = N.trapz(e, s), resil = yielded ? sy * sy / (2 * E) : uts * uts / (2 * E), fT = hasT ? linfit(e.slice(best.i, best.i + w), et.slice(best.i, best.i + w)) : null, nu = fT ? -fT.slope : NaN, urel = E !== 0 ? Math.abs(best.f.seSlope / E) : 0;
  o.K('E', "Young's modulus", E, 'MPa', `Least-squares slope of the segment between 5 % and ${F(100 * best.up, 2)} % of the maximum stress (${best.m} points, R² = ${F(best.f.r2, 5)}), the upper limit being lowered automatically until the segment is linear; standard error ${F(100 * urel, 2)} %. Machine compliance must already be removed (extensometer strain).`, best.f.r2 > 0.99);
  o.K('sy', `${inp.tens_offset} % offset yield strength`, yielded ? sy : 0, 'MPa', yielded ? `Intersection of the curve with a line of slope E offset by ${inp.tens_offset} % strain, at total strain ${F(ey, 3)}.` : 'The curve never falls below the offset line — no yield point before the end of the record (brittle failure or test stopped early); reported as 0.', yielded);
  o.K('uts', 'Ultimate tensile strength', uts, 'MPa', `Maximum engineering stress, reached at strain ${F(e[iu], 3)}; beyond it deformation localises (necking) in ductile materials.`);
  o.K('nu', "Poisson's ratio", fin(nu), '', hasT ? `Negative slope of transverse against axial strain over the elastic segment (± ${F(fT.seSlope, 2)}). Negative values identify auxetic behaviour; isotropic continua are bounded by −1 < ν < 0.5, architected lattices need not be isotropic.` : 'No transverse-strain column supplied — reported as 0.', hasT ? null : false);
  o.K('ef', 'Strain at end of record', e[n - 1], '', 'Total engineering strain at the last data point; equals the elongation at fracture only if the test ran to failure.');
  o.K('tough', 'Strain energy density to end of record', tough, 'MJ/m³', 'Area under the engineering stress–strain curve by trapezoidal integration (MPa × strain = MJ/m³): the energy absorbed per unit volume.');
  o.K('resil', 'Modulus of resilience', resil, 'MJ/m³', 'σ_y²/(2E): elastic energy density stored at yield (at the maximum stress when no yield point was found).');
  o.K('ratio', 'Yield-to-tensile ratio', yielded ? sy / uts : 1, '', 'σ_y/UTS: low values indicate large work-hardening capacity; values near 1 little reserve between yielding and failure.');
  o.S('ss', 'Engineering stress–strain curve', 'strain (–)', 'stress (MPa)', [['measured', e, s], ['elastic line', [e0, e0 + uts / E], [0, uts]], ['offset line', [e0 + off, e0 + off + uts / E], [0, uts]]], 'The elastic line is the fitted modulus; where the offset line crosses the curve is the proof stress.');
  o.S(hasT ? 'trans' : 'tangent', hasT ? 'Transverse against axial strain' : 'Tangent modulus', 'axial strain (–)', hasT ? 'transverse strain (–)' : 'dσ/dε (MPa)', hasT ? [['measured', e, et], ['elastic fit', [e[best.i], e[best.i + w - 1]], [fT.intercept + fT.slope * e[best.i], fT.intercept + fT.slope * e[best.i + w - 1]]]] : [['tangent modulus', e.slice(1), e.slice(1).map((v, i) => (s[i + 1] - s[i]) / Math.max(1e-12, v - e[i]))]], hasT ? "The slope in the elastic range is −ν; a positive slope means the specimen widens when stretched (negative Poisson's ratio)." : 'Numerical derivative of the curve; the plateau at small strain is the elastic modulus and its drop marks yielding.');
  o.E('constitutive', "Hooke's law (uniaxial)", 'σ = E·ε,  ε_transverse = −ν·ε_axial'); o.E('formula', 'Offset yield construction', `σ_y: σ(ε) = E·(ε − ${F(off, 2)})`); o.E('formula', 'Energy densities', 'U = ∫σ dε;  U_r = σ_y²/(2E)');
  o.C('Linear elastic segment found', best.f.r2 > 0.99, `R² = ${F(best.f.r2, 5)}`); o.C('Yield point reached', yielded, yielded ? `σ_y = ${F(sy, 4)} MPa` : 'no intersection with the offset line');
  if (hasT) o.C("Poisson's ratio within isotropic stability bounds", nu > -1 && nu < 0.5, `ν = ${F(nu, 3)} (−1 < ν < 0.5 for an isotropic solid; anisotropic lattices may exceed these)`);
  if (hasT && nu < 0) o.W(`Negative Poisson's ratio (${F(nu, 3)}): auxetic response. Check that transverse strain is measured over several unit cells and away from the grips, where the lattice is constrained.`);
  o.W('Engineering (not true) stress and strain; strain is taken as given — crosshead-derived strain includes machine compliance and underestimates the modulus.');
  o.R('high', `Use an extensometer or DIC for strain: the modulus (${F(E, 4)} MPa) from crosshead displacement is typically low by tens of percent.`, 'Machine compliance in series with the specimen.');
  o.R('medium', hasT && nu < 0 ? 'For the auxetic lattice, test at least 5 × 5 unit cells and report ν over the central cells; edge effects make small specimens appear less auxetic.' : 'Test at least five specimens and report mean and standard deviation of σ_y and UTS; single curves do not capture batch scatter.', 'Size effects and specimen-to-specimen variability.');
  Object.assign(o.publish, { 'prop.E': E * 1e6, 'eff.E': E * 1e6, 'cred.u_rel': urel, 'cred.rmse': best.f.s }, yielded ? { 'prop.yield': sy * 1e6 } : {}, hasT ? { 'prop.nu': nu, 'eff.nu': nu } : {}, { 'mech.uts': uts * 1e6 });
  o.numerics = { method: 'Linear regression over an automatically shrunk elastic stress window; offset-line intersection; trapezoidal integration', dof: n, residual: best.f.s };
  o.summary = `The elastic segment gives E = ${F(E, 4)} MPa; the ${inp.tens_offset} % offset yield strength is ${yielded ? F(sy, 4) + ' MPa' : 'not reached'} and the maximum stress ${F(uts, 4)} MPa.${hasT ? ` Poisson's ratio from the transverse strain is ${F(nu, 3)}${nu < 0 ? ' (auxetic)' : ''}.` : ''} Energy absorbed to the end of the record is ${F(tough, 3)} MJ/m³.`;
}

function mVSM(inp, o) {
  const t = o.tbl(inp.vsm_table, 2, 10, 'Hysteresis loop'), H = t.rows.map(r => r[0]), M = t.rows.map(r => r[1]), n = H.length, Hmax = N.maxOf(H.map(Math.abs)), hc = [], mr = [];
  for (let i = 1; i < n; i++) { if (M[i - 1] * M[i] < 0 || (M[i] === 0 && M[i - 1] !== 0)) hc.push(H[i - 1] - M[i - 1] * (H[i] - H[i - 1]) / (M[i] - M[i - 1])); if (H[i - 1] * H[i] < 0 || (H[i] === 0 && H[i - 1] !== 0)) mr.push(H[i] === 0 ? M[i] : M[i - 1] - H[i - 1] * (M[i] - M[i - 1]) / (H[i] - H[i - 1])); }
  const hi = H.map((_, i) => i).filter(i => Math.abs(H[i]) > 0.7 * Hmax), fa = hi.length >= 4 ? linfit(hi.map(i => 1 / Math.abs(H[i])), hi.map(i => Math.abs(M[i]))) : null, Mmax = N.maxOf(M.map(Math.abs)), Ms = fa && fa.intercept > 0 && fa.intercept < 2 * Mmax ? fa.intercept : Mmax;
  const Hc = hc.length ? N.mean(hc.map(Math.abs)) : 0, Mr = mr.length ? N.mean(mr.map(Math.abs)) : 0, bias = hc.length >= 2 ? N.mean(hc) : 0, sq = Ms > 0 ? Mr / Ms : 0; let area = 0; for (let i = 0; i < n; i++) { const j = (i + 1) % n; area += H[i] * M[j] - H[j] * M[i]; }
  const loss = N.MU0 * Math.abs(area) / 2 * 1e6, sat = Ms > 0 ? Mmax / Ms : 0, closed = Math.hypot(H[0] - H[n - 1], 0) < 0.1 * Hmax;
  o.K('Ms', 'Saturation magnetisation', Ms, 'kA/m', fa && Ms === fa.intercept ? `Extrapolation of |M| against 1/|H| to infinite field over the ${hi.length} points above 70 % of the maximum field (law of approach to saturation, first term); the largest measured value is ${F(Mmax, 4)} kA/m.` : `Largest measured magnetisation (no reliable high-field extrapolation was possible); a lower bound if the loop is not saturated.`, sat > 0.95);
  o.K('Mr', 'Remanent magnetisation', Mr, 'kA/m', `Mean |M| at the ${mr.length} zero-field crossings, by linear interpolation between neighbouring points.`, mr.length >= 2);
  o.K('Hc', 'Coercive field', Hc, 'kA/m', `Mean |H| at the ${hc.length} zero-magnetisation crossings (= ${F(N.MU0 * Hc * 1e6, 4)} mT). Near-zero coercivity with zero remanence is the signature of superparamagnetic nanoparticles at the measurement temperature.`, hc.length >= 2);
  o.K('sq', 'Squareness M_r/M_s', sq, '', '0.5 is expected for randomly oriented uniaxial single-domain particles (Stoner–Wohlfarth), about 0.83 for cubic anisotropy; lower values indicate multi-domain or superparamagnetic behaviour, higher values texture.');
  o.K('sigma_s', 'Specific saturation magnetisation', Ms * 1e3 / inp.rho, 'A·m²/kg', `M_s divided by the density ${inp.rho} kg/m³ (numerically equal to emu/g); the usual way to compare nanoparticle powders with bulk.`);
  o.K('loss', 'Hysteresis loss per cycle', loss, 'J/m³', 'µ₀·∮H dM evaluated with the shoelace formula over the supplied loop; multiplied by frequency it gives the hysteresis heating power density.');
  o.K('bias', 'Loop shift (exchange-bias field)', bias, 'kA/m', 'Mean of the signed coercive fields; a non-zero value indicates exchange bias (e.g. an oxide shell on a metallic core) or an uncorrected field offset.');
  o.K('sat', 'Degree of saturation reached', sat, '', 'Largest measured magnetisation divided by the extrapolated M_s; below about 0.95 a minor loop may have been measured and H_c, M_r are then underestimated.', sat > 0.95);
  o.S('loop', 'Magnetisation loop', 'H (kA/m)', 'M (kA/m)', [['measured', H, M]], 'Full hysteresis loop; the intercepts with the axes are the coercive field and the remanence.');
  o.S('approach', 'Approach to saturation', '1/|H| (m/kA)', '|M| (kA/m)', [['high-field data', hi.map(i => 1 / Math.abs(H[i])), hi.map(i => Math.abs(M[i]))], ...(fa ? [['linear extrapolation', [0, 1 / (0.7 * Hmax)], [fa.intercept, fa.intercept + fa.slope / (0.7 * Hmax)]]] : [])], 'The intercept at 1/H → 0 is the saturation magnetisation; a slope that does not flatten suggests a paramagnetic or disordered-surface contribution.');
  o.E('constitutive', 'Magnetic constitutive relation', 'B = µ₀·(H + M)'); o.E('formula', 'Law of approach to saturation (first term)', 'M(H) = M_s·(1 − a/H)'); o.E('formula', 'Hysteresis loss', 'W = µ₀·∮ H dM'); o.E('formula', 'Squareness', 'S = M_r/M_s');
  o.C('Both branches present', hc.length >= 2 && mr.length >= 2, `${hc.length} coercive and ${mr.length} remanence crossings found`); o.C('Loop saturated', sat > 0.95, `M_max/M_s = ${F(sat, 3)}`); o.C('Loop closed', closed, `start and end fields differ by ${F(Math.abs(H[0] - H[n - 1]), 3)} kA/m`);
  o.C('Squareness physical', sq >= 0 && sq <= 1, `M_r/M_s = ${F(sq, 3)}`);
  if (hc.length < 2) o.W('Fewer than two zero crossings of M: the data are not a complete loop, so coercivity and remanence are incomplete.');
  o.W('No demagnetising-field correction is applied (H is the applied field); for dense or plate-like samples the internal field is lower and the loop is sheared. Volume magnetisation requires an accurate sample volume.');
  o.R('high', sat > 0.95 ? `M_s = ${F(Ms, 4)} kA/m is well defined; compare σ_s = ${F(Ms * 1e3 / inp.rho, 3)} A·m²/kg with the bulk value to estimate the magnetically dead surface layer of the nanoparticles.` : `Increase the maximum field beyond ${F(Hmax, 3)} kA/m: the loop reaches only ${F(100 * sat, 3)} % of the extrapolated saturation.`, 'Degree of saturation.');
  o.R('medium', 'Correct for the demagnetising factor of the sample shape and for the diamagnetic holder signal before comparing coercivity and susceptibility between samples.', 'Internal field H_i = H − N_d·M.');
  Object.assign(o.publish, { 'prop.Ms': Ms * 1e3, 'mag.Hc': Hc * 1e3, 'mag.Mr': Mr * 1e3, 'cred.u_rel': fa ? Math.abs(fa.seInt / Ms) : 0.05, 'cred.rmse': fa ? fa.s / Ms : 0 });
  o.numerics = { method: 'Linear interpolation of axis crossings; high-field linear extrapolation in 1/H; shoelace loop area', dof: n };
  o.summary = `The loop gives M_s = ${F(Ms, 4)} kA/m, M_r = ${F(Mr, 4)} kA/m and H_c = ${F(Hc, 4)} kA/m (squareness ${F(sq, 3)}). The hysteresis loss is ${F(loss, 3)} J/m³ per cycle and the loop reaches ${F(100 * sat, 3)} % of the extrapolated saturation.`;
}

function mGUM(inp, o) {
  if (isSynthetic(inp.gum_table)) o.synthetic = true;
  const red = GUM[inp.gum_red] || GUM.scherrer, tb = parseNamed(inp.gum_table), miss = [], vars = red.vars.map(v => { const q = tb[v[0].toLowerCase()]; if (!q) miss.push(v[0]); return q ? [v[0], q[0], q[1]] : [...v]; });
  const x = vars.map(v => v[1]), u = vars.map(v => v[2]), f = z => fin(red.f(z), NaN), y = f(x), kf = inp.gum_k;
  const c = x.map((xi, i) => { const h = Math.max(Math.abs(xi) * 1e-5, u[i] * 1e-2, 1e-12), a = [...x], b = [...x]; a[i] += h; b[i] -= h; return fin((f(a) - f(b)) / (2 * h), 0); });
  const contrib = c.map((ci, i) => (ci * u[i]) ** 2), uc = Math.sqrt(N.sum(contrib)), dom = N.argmax(contrib), r = N.rng(inp.seed), mc = [];
  for (let k = 0; k < 4000; k++) { const v = f(x.map((xi, i) => xi + r.normal(0, u[i]))); if (Number.isFinite(v)) mc.push(v); }
  const mcM = mc.length ? N.mean(mc) : y, mcS = mc.length > 1 ? N.std(mc) : 0, lo = mc.length ? N.quantile(mc, 0.025) : y, hi = mc.length ? N.quantile(mc, 0.975) : y, lin = uc > 0 ? Math.abs(mcS - uc) / uc : 0;
  o.K('y', `Result: ${red.label.split(' ')[0]} reduction`, fin(y), red.unit, `Value of ${red.label} at the tabulated input estimates.`, Number.isFinite(y));
  o.K('uc', 'Combined standard uncertainty', uc, red.unit, `Law of propagation of uncertainty for uncorrelated inputs, u_c² = Σ(cᵢ·uᵢ)², with sensitivity coefficients from central finite differences; ${F(100 * uc / Math.abs(y || 1), 3)} % of the result.`, uc / Math.abs(y || 1) < 0.1);
  o.K('U', `Expanded uncertainty (k = ${kf})`, kf * uc, red.unit, `U = k·u_c; k = 2 corresponds to about 95 % coverage for a normal output distribution. Report as ${F(fin(y), 4)} ± ${F(kf * uc, 2)} ${red.unit}.`);
  o.K('dom', 'Dominant contributor', vars[dom][0], '', `${vars[dom][0]} accounts for ${F(100 * contrib[dom] / (uc * uc || 1), 3)} % of the combined variance; reducing any other input uncertainty has little effect until this one is improved.`);
  o.K('dom_share', 'Variance share of dominant input', uc > 0 ? contrib[dom] / (uc * uc) : 0, '', 'Fraction of u_c² carried by the largest contribution (cᵢuᵢ)²; a value near 1 means a single quantity controls the measurement quality.');
  o.K('mc_sd', 'Monte Carlo standard deviation', mcS, red.unit, `Standard deviation of ${mc.length} seeded samples propagated through the full non-linear reduction (normal inputs); differs from the linear u_c by ${F(100 * lin, 3)} %.`, lin < 0.1);
  o.K('mc_lo', 'Monte Carlo 2.5 % quantile', lo, red.unit, 'Lower end of the 95 % coverage interval from the Monte Carlo output distribution; asymmetric intervals reveal non-linearity that k·u_c hides.');
  o.K('mc_hi', 'Monte Carlo 97.5 % quantile', hi, red.unit, 'Upper end of the 95 % coverage interval from the Monte Carlo output distribution.');
  o.T('budget', 'Uncertainty budget', ['quantity', 'estimate', 'standard uncertainty', 'sensitivity cᵢ', 'contribution |cᵢ|·uᵢ', 'share of variance (%)'], vars.map((v, i) => [v[0], v[1], v[2], c[i], Math.abs(c[i]) * u[i], 100 * contrib[i] / (uc * uc || 1)]), 'One row per input quantity; contributions add in quadrature to the combined standard uncertainty (inputs assumed uncorrelated).');
  o.S('share', 'Variance contributions', 'input index', 'share of variance (%)', [['share', vars.map((_, i) => i + 1), contrib.map(v => 100 * v / (uc * uc || 1))]], `Share of the combined variance per input in table order (${vars.map((v, i) => `${i + 1} = ${v[0]}`).join(', ')}).`);
  if (mc.length > 10) { const h = N.histogram(mc, 30); o.S('mc', 'Monte Carlo output distribution', `result (${red.unit})`, 'count', [['samples', h.x, h.y], ['normal with u_c', h.x, h.x.map(v => mc.length * h.width * Math.exp(-((v - y) ** 2) / (2 * uc * uc || 1e-300)) / (Math.sqrt(2 * PI) * (uc || 1e-300)))]], 'Histogram of the propagated samples against the normal distribution implied by the linear budget; skew signals non-linearity.'); }
  o.E('governing', 'Law of propagation of uncertainty (GUM)', 'u_c²(y) = Σᵢ (∂f/∂xᵢ)²·u²(xᵢ)  (uncorrelated inputs);  U = k·u_c'); o.E('formula', 'Measurement model', red.eq); o.E('formula', 'Finite-difference sensitivity', 'cᵢ = [f(xᵢ + h) − f(xᵢ − h)]/(2h)');
  o.C('Linear budget agrees with Monte Carlo', lin < 0.1, `u_c = ${F(uc, 3)}, Monte Carlo ${F(mcS, 3)} ${red.unit}`); o.C('Result finite at the estimates', Number.isFinite(y), `y = ${F(fin(y), 5)} ${red.unit}`);
  if (miss.length) o.W(`Quantities not found in the budget table (${miss.join(', ')}): typical default estimates and uncertainties were used for them. Table rows must be named ${red.vars.map(v => v[0]).join(', ')}.`);
  if (lin >= 0.1) o.W('The model is noticeably non-linear over the input uncertainties: use the Monte Carlo interval instead of k·u_c.');
  o.W('Inputs are treated as uncorrelated and normally distributed; correlated inputs (e.g. two readings from the same instrument) need covariance terms.');
  o.R('high', `Reduce the uncertainty of ${vars[dom][0]} first: it carries ${F(100 * contrib[dom] / (uc * uc || 1), 3)} % of the variance, so halving it lowers u_c from ${F(uc, 3)} to about ${F(Math.sqrt(Math.max(0, uc * uc - 0.75 * contrib[dom])), 3)} ${red.unit}.`, 'Variance decomposition of the budget.');
  o.R('medium', 'Document the origin of each standard uncertainty (type A from repeats, type B from certificates or resolution with the assumed distribution) alongside the budget.', 'Traceability of the uncertainty statement.');
  Object.assign(o.publish, { 'cred.u_rel': uc / Math.abs(y || 1), 'cred.u_combined': uc, 'cred.U_expanded': kf * uc });
  o.numerics = { method: 'GUM first-order propagation with central finite-difference sensitivities; Monte Carlo cross-check', dof: vars.length, steps: mc.length };
  o.summary = `${red.label}: y = ${F(fin(y), 5)} ${red.unit} with combined standard uncertainty ${F(uc, 3)} ${red.unit} and expanded uncertainty ${F(kf * uc, 3)} ${red.unit} (k = ${kf}). The dominant contributor is ${vars[dom][0]} with ${F(100 * contrib[dom] / (uc * uc || 1), 3)} % of the variance; Monte Carlo gives a 95 % interval of ${F(lo, 4)}–${F(hi, 4)} ${red.unit}.`;
}

function mRepeat(inp, o) {
  const t = o.tbl(inp.rep_table, 1, 3, 'Replicate measurements'), v = t.rows.map(r => r[0]), n = v.length, m = N.mean(v), s = N.std(v), se = s / Math.sqrt(n), tt = tq(0.975, n - 1), uB = inp.rep_uB, uc = Math.hypot(se, uB), U = 2 * uc, bias = m - inp.rep_ref, En = bias / Math.max(1e-300, Math.hypot(U, inp.rep_Uref));
  const dev = v.map(x => Math.abs(x - m)), ig = N.argmax(dev), G = s > 0 ? dev[ig] / s : 0, tc = tq(1 - 0.05 / (2 * n), Math.max(1, n - 2)), Gc = (n - 1) / Math.sqrt(n) * Math.sqrt(tc * tc / (n - 2 + tc * tc)), need = Math.ceil((1.96 * s / Math.max(1e-300, inp.rep_target / 100 * Math.abs(m))) ** 2);
  o.K('mean', 'Mean of replicates', m, '', `Arithmetic mean of ${n} replicate measurements; 95 % confidence interval ± ${F(tt * se, 3)} (Student-t).`);
  o.K('s', 'Repeatability standard deviation', s, '', `Sample standard deviation under repeatability conditions; relative ${F(100 * s / Math.abs(m || 1), 3)} %. The repeatability limit r = 2.8·s = ${F(2.8 * s, 3)} is the largest difference expected between two single results (95 %).`);
  o.K('se', 'Standard uncertainty of the mean (type A)', se, '', 's/√n: the part of the uncertainty that further repetition can reduce.');
  o.K('U', 'Expanded uncertainty of the laboratory value (k = 2)', U, '', `Combination of the type-A term with the instrument type-B standard uncertainty ${uB}: U = 2·√(s²/n + u_B²).`);
  o.K('bias', 'Bias against the reference value', bias, '', `Mean minus the reference value ${inp.rep_ref} (${F(100 * bias / Math.abs(inp.rep_ref || 1), 3)} %); significant only if it exceeds the combined expanded uncertainty.`, Math.abs(bias) <= Math.hypot(U, inp.rep_Uref));
  o.K('En', 'E_n number', En, '', 'E_n = (x_lab − X_ref)/√(U_lab² + U_ref²) with expanded (k = 2) uncertainties; |E_n| ≤ 1 means the result agrees with the reference within the stated uncertainties, |E_n| > 1 calls for investigation.', Math.abs(En) <= 1);
  o.K('G', 'Grubbs outlier statistic', G, '', `Largest absolute deviation divided by s (value ${F(v[ig], 5)}); critical value at 5 % significance ${F(Gc, 4)} for n = ${n}. Exceeding it flags a possible outlier, to be investigated, not silently removed.`, G < Gc);
  o.K('need', 'Replicates needed for target precision', need, '', `Number of repeats for a 95 % half-width of ${inp.rep_target} % of the mean from type-A scatter alone, n = (1.96·s/(ε·x̄))².`, n >= need);
  o.S('run', 'Run chart of replicates', 'measurement number', 'value', [['replicates', v.map((_, i) => i + 1), v], ['mean', [1, n], [m, m]], ['mean + 2s', [1, n], [m + 2 * s, m + 2 * s]], ['mean − 2s', [1, n], [m - 2 * s, m - 2 * s]], ['reference', [1, n], [inp.rep_ref, inp.rep_ref]]], 'Sequence of results with warning limits; a trend indicates drift, which violates repeatability conditions.');
  const h = N.histogram(v, Math.max(4, Math.round(Math.sqrt(n)) + 2));
  o.S('hist', 'Distribution of replicates', 'value', 'count', [['counts', h.x, h.y]], 'Histogram of the replicate values; with few points it only reveals gross asymmetry or isolated outliers.');
  o.E('formula', 'Type-A evaluation', 'x̄ = Σxᵢ/n;  s² = Σ(xᵢ − x̄)²/(n − 1);  u_A = s/√n'); o.E('formula', 'E_n number', 'E_n = (x_lab − X_ref)/√(U_lab² + U_ref²)'); o.E('formula', 'Grubbs test', 'G = max|xᵢ − x̄|/s;  G_crit = ((n − 1)/√n)·√(t²/(n − 2 + t²)),  t = t(α/(2n), n − 2)');
  o.C('Calibration check |E_n| ≤ 1', Math.abs(En) <= 1, `E_n = ${F(En, 3)}`); o.C('No outlier by Grubbs test', G < Gc, `G = ${F(G, 3)} against ${F(Gc, 3)}`); o.C('At least 5 replicates', n >= 5, `${n} values`);
  if (n < 5) o.W('Fewer than five replicates: the standard deviation itself is uncertain by more than 35 %.');
  if (Math.abs(En) > 1) o.W('|E_n| > 1: the instrument reading disagrees with the reference beyond the stated uncertainties — recalibrate or revise the uncertainty budget before using the results.');
  o.R('high', Math.abs(En) <= 1 ? `The calibration check passes (E_n = ${F(En, 2)}); keep the reference measurement in the routine at a fixed interval and plot it on a control chart.` : `Investigate the bias of ${F(bias, 3)}: check calibration date, reference handling and environmental conditions before any further measurement.`, 'E_n criterion of proficiency testing.');
  o.R('medium', se > uB ? `Type-A scatter dominates (u_A = ${F(se, 2)} against u_B = ${uB}): more replicates will reduce the uncertainty.` : `The instrument type-B term dominates (u_B = ${uB} against u_A = ${F(se, 2)}): more replicates will not help — improve the calibration instead.`, 'Comparison of type-A and type-B contributions.');
  Object.assign(o.publish, { 'cred.u_rel': uc / Math.abs(m || 1), 'cred.En': En, 'cred.repeatability': s, 'cred.rmse': s });
  o.numerics = { method: 'Type-A statistics, Grubbs test and E_n number', dof: n - 1 };
  o.summary = `${n} replicates give ${F(m, 5)} ± ${F(U, 2)} (k = 2) with a repeatability standard deviation of ${F(s, 3)}. Against the reference ${inp.rep_ref} ± ${inp.rep_Uref} the bias is ${F(bias, 3)} and E_n = ${F(En, 3)}, so the calibration check ${Math.abs(En) <= 1 ? 'passes' : 'fails'}.`;
}

const MODES = {
  xrd: ['XRD — peak fitting, phase identification, lattice refinement', mXRD], dls: ['Particle sizing — DLS cumulant analysis', mDLS], counts: ['Particle sizing — counted diameters (TEM/SEM)', mCounts], bet: ['Gas adsorption — BET area and pore size', mBET],
  tauc: ['Optical — Tauc band gap', mTauc], lspr: ['Optical — Beer–Lambert and LSPR size estimate', mLSPR], raman: ['Raman / FTIR — band table and mode matching', mRaman], op: ['Nanoindentation — Oliver–Pharr', mOP], hertz: ['Indentation — Hertz (spherical / AFM)', mHertz],
  fpp: ['Electrical — four-point probe', mFPP], hall: ['Electrical — Hall effect and van der Pauw', mHall], eis: ['Electrical — impedance spectrum (Randles fit)', mEIS], flash: ['Thermal — laser-flash diffusivity', mFlash], tensile: ['Mechanical — tensile curve', mTensile], vsm: ['Magnetic — VSM hysteresis loop', mVSM],
  gum: ['Metrology — GUM uncertainty budget', mGUM], repeat: ['Metrology — repeatability and calibration check (E_n)', mRepeat],
};
const num = (key, label, unit, value, min, max, modes, help, extra = {}) => ({ key, label, unit, value, min, max, group: extra.group || 'Analysis settings', showIf: { mode: modes }, help, ...extra });
const tab = (key, label, columns, value, mode, help) => ({ key, label, type: 'table', columns, value, group: 'Measured data', showIf: { mode: [mode] }, help });

export default {
  id: 23,
  tagline: 'Reduces measured diffraction, sizing, adsorption, spectroscopy, indentation, electrical, thermal, mechanical and magnetic data to material properties with uncertainty.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'xrd', options: Object.entries(MODES).map(([value, m]) => ({ value, label: m[0] })), group: 'Analysis', help: 'Which measurement to reduce. Each analysis has its own data table (CSV text) and settings; the default tables are seeded synthetic demonstration data.' },
    tab('xrd_table', 'XRD pattern (2θ, intensity)', ['two_theta_deg', 'intensity'], XRD_ZNO, 'xrd', 'Powder diffraction pattern on a uniform 2θ grid. The default is a simulated ZnO pattern (25 nm crystallites) with background and counting noise.'),
    num('xrd_lambda', 'X-ray wavelength', 'Å', 1.5406, 0.5, 2.3, ['xrd'], 'Wavelength of the radiation (Cu Kα₁ 1.5406, Co Kα 1.7890, Mo Kα 0.7107 Å). All d-spacings and the phase match depend on it.'),
    num('xrd_inst', 'Instrumental breadth (FWHM)', '°', 0.08, 0, 0.5, ['xrd'], 'Peak width of a strain-free, coarse-grained standard on the same instrument; removed from the observed width before the Scherrer and Williamson–Hall analyses.'),
    num('xrd_K', 'Scherrer constant K', '', 0.9, 0.6, 1.3, ['xrd'], 'Shape factor: about 0.9 for FWHM of roughly equiaxed crystallites. The size scales linearly with it.'),
    num('xrd_pos_tol', '2θ position uncertainty', '°', 0.08, 0.01, 0.5, ['xrd'], 'Standard uncertainty of peak positions including zero shift and sample displacement; sets the matching tolerance for phase identification and the weights of the lattice refinement.'),
    { key: 'xrd_corr', label: 'Instrument-broadening correction', type: 'select', value: 'quadrature', options: [{ value: 'quadrature', label: 'Quadrature (Gaussian profiles): β² = B² − b²' }, { value: 'linear', label: 'Linear (Lorentzian profiles): β = B − b' }], group: 'Analysis settings', showIf: { mode: ['xrd'] }, help: 'How the instrumental breadth is removed. Real profiles lie between the two limits; the difference is a measure of the model uncertainty of the size.' },
    tab('dls_table', 'DLS autocorrelation (τ, g₂)', ['tau_s', 'g2'], DLS_DEF, 'dls', 'Normalised intensity autocorrelation function against delay time. Default: simulated 60 nm particles with PDI 0.08 in water at 25 °C, 633 nm, 173° detection.'),
    num('dls_T', 'Measurement temperature', '°C', 25, 5, 90, ['dls'], 'Sample temperature; enters Stokes–Einstein directly and through the viscosity, which must be given for this temperature.'),
    num('dls_visc', 'Dispersant viscosity', 'mPa·s', 0.89, 0.2, 100, ['dls'], 'Dynamic viscosity of the liquid at the measurement temperature (water: 0.89 at 25 °C). The size is proportional to 1/η.', { log: true }),
    num('dls_n', 'Dispersant refractive index', '', 1.33, 1, 1.8, ['dls'], 'Refractive index of the liquid at the laser wavelength; enters the scattering vector.'),
    num('dls_lambda', 'Laser wavelength', 'nm', 633, 400, 1100, ['dls'], 'Vacuum wavelength of the DLS laser.'),
    num('dls_angle', 'Scattering angle', '°', 173, 10, 175, ['dls'], 'Detection angle (173° backscatter or 90° in most instruments).'),
    tab('tem_table', 'Counted particle diameters', ['diameter_nm'], TEM_DEF, 'counts', 'One measured diameter per row from TEM/SEM image analysis. Default: 200 values drawn from a log-normal population (median 42 nm, geometric standard deviation 1.22).'),
    num('tem_target', 'Target precision of the mean', '%', 2, 0.1, 20, ['counts'], 'Desired 95 % confidence half-width of the mean diameter; used to compute how many particles must be counted.'),
    tab('bet_table', 'Adsorption isotherm (p/p₀, uptake)', ['p_rel', 'v_ads_cm3_g'], BET_DEF, 'bet', 'Nitrogen adsorption branch at 77 K: relative pressure and adsorbed volume in cm³(STP)/g. Default: simulated mesoporous oxide (≈ 100 m²/g).'),
    num('bet_sigma', 'Adsorbate cross-section', 'nm²', 0.162, 0.1, 0.3, ['bet'], 'Area occupied by one adsorbed molecule (N₂ at 77 K: 0.162 nm²).'),
    num('rho', 'Sample (skeletal) density', 'kg/m³', 5606, 500, 23000, ['bet', 'flash', 'vsm'], 'Solid density: converts surface area to an equivalent particle size, diffusivity to conductivity and volume to mass magnetisation.', { bus: 'prop.density', group: 'Sample' }),
    tab('uv_table', 'Absorbance spectrum (λ, A)', ['wavelength_nm', 'absorbance'], UV_DEF, 'tauc', 'UV–vis absorbance against wavelength across the absorption edge. Default: simulated direct-gap semiconductor (3.30 eV) with an Urbach tail and offset.'),
    { key: 'tauc_type', label: 'Transition type', type: 'select', value: 'direct', options: [{ value: 'direct', label: 'Direct allowed: (αhν)²' }, { value: 'indirect', label: 'Indirect allowed: (αhν)^½' }], group: 'Analysis settings', showIf: { mode: ['tauc'] }, help: 'Exponent of the Tauc transform. It is an assumption about the band structure (ZnO, CdSe: direct; anatase TiO₂, Si: indirect) and changes the result.' },
    { key: 'tauc_base', label: 'Absorbance baseline', type: 'select', value: 'subtract', options: [{ value: 'subtract', label: 'Subtract the sub-gap baseline level' }, { value: 'none', label: 'Use absorbance as supplied' }], group: 'Analysis settings', showIf: { mode: ['tauc'] }, help: 'Removes a constant offset (reflection, scattering) estimated from the lowest absorbance values before the Tauc transform.' },
    tab('lspr_table', 'Extinction spectrum (λ, A)', ['wavelength_nm', 'absorbance'], LSPR_DEF, 'lspr', 'Extinction of a metal colloid. Default: simulated spectrum of about 20 nm gold spheres in water.'),
    { key: 'lspr_metal', label: 'Metal', type: 'select', value: 'Au', options: [{ value: 'Au', label: 'Gold spheres (A_spr/A_450 correlation)' }, { value: 'Ag', label: 'Silver spheres (indicative peak position)' }], group: 'Analysis settings', showIf: { mode: ['lspr'] }, help: 'Selects the empirical size relation and the wavelength window in which the plasmon band is searched.' },
    num('lspr_eps', 'Molar extinction coefficient', 'L/mol/cm', 8.8e8, 1e3, 1e12, ['lspr'], 'Extinction coefficient at the band maximum per mole of particles (about 10⁹ L/mol/cm for 20 nm gold; it scales roughly with d³).', { log: true }),
    num('lspr_path', 'Optical path length', 'cm', 1, 0.01, 10, ['lspr'], 'Cuvette path length for the Beer–Lambert concentration.'),
    tab('raman_table', 'Raman/FTIR spectrum (shift, intensity)', ['shift_cm1', 'intensity'], RAMAN_DEF, 'raman', 'Spectrum against wavenumber. Default: simulated graphene on silicon (D, G, 2D bands and the Si phonon).'),
    num('raman_tol', 'Mode-matching tolerance', 'cm⁻¹', 8, 1, 30, ['raman'], 'Maximum distance between a fitted band and a library mode for an assignment (three times wider for second-order bands above 2000 cm⁻¹).'),
    num('raman_laser', 'Excitation wavelength', 'nm', 532, 325, 1064, ['raman'], 'Laser wavelength; enters the empirical graphitic crystallite-size relation as λ⁴.'),
    num('raman_stress', 'Si stress coefficient', 'MPa per cm⁻¹', 500, 100, 1000, ['raman'], 'Stress per unit shift of the 520.7 cm⁻¹ silicon phonon: about 500 for uniaxial and 250 for equi-biaxial stress in (100) Si.'),
    tab('op_table', 'Load–displacement curve (h, P)', ['depth_nm', 'load_mN'], OP_DEF, 'op', 'Loading and unloading segments of an instrumented indentation, in acquisition order. Default: simulated Berkovich indent in a fused-silica-like solid (E = 72 GPa, H = 9.3 GPa).'),
    num('op_nu', "Sample Poisson's ratio", '', 0.17, 0, 0.5, ['op'], "Needed to convert the reduced modulus into Young's modulus; an error of 0.1 changes E by only a few percent.", { bus: 'prop.nu' }),
    num('op_beta', 'Indenter geometry factor β', '', 1.034, 1, 1.1, ['op'], 'Correction for the non-axisymmetric Berkovich shape (1.034); 1.0 for an axisymmetric indenter.'),
    num('op_eps', 'Contact-depth factor ε', '', 0.75, 0.7, 1, ['op'], 'Geometry constant of the Oliver–Pharr sink-in model: 0.75 for a paraboloid/Berkovich tip, 0.72 for a cone, 1 for a flat punch.'),
    num('op_C1', 'Tip-rounding area term C₁', 'nm', 0, 0, 5000, ['op'], 'First correction term of the calibrated area function A = 24.5·h_c² + C₁·h_c; zero for an ideally sharp Berkovich tip.'),
    num('op_film', 'Film thickness (0 = bulk)', 'nm', 0, 0, 1e5, ['op'], 'If the sample is a coating, its thickness is used for the 10 % depth rule that limits substrate influence.'),
    tab('hz_table', 'Force–indentation curve (δ, F)', ['indentation_nm', 'force_nN'], HZ_DEF, 'hertz', 'Approach segment of a spherical-tip or colloidal-probe AFM indentation including some non-contact baseline. Default: simulated 2 MPa elastomer with a 2.5 µm probe.'),
    num('hz_R', 'Tip radius', 'nm', 2500, 2, 1e5, ['hertz'], 'Radius of the spherical indenter or colloidal probe; the modulus scales with R^(−1/2).', { log: true }),
    num('hz_nu', "Sample Poisson's ratio", '', 0.5, 0, 0.5, ['hertz'], 'About 0.5 for elastomers and hydrated soft matter, 0.2–0.35 for stiff solids.'),
    tab('fpp_table', 'Four-point-probe I–V data', ['current_mA', 'voltage_mV'], FPP_DEF, 'fpp', 'Current through the outer probes and voltage across the inner probes, both polarities. Default: simulated 150 nm transparent-conductor film (ρ = 2×10⁻⁶ Ω·m).'),
    num('el_thick', 'Film thickness', 'nm', 150, 1, 1e6, ['fpp', 'hall'], 'Thickness of the conducting layer; converts sheet resistance into resistivity and Hall voltage into carrier density.', { log: true }),
    num('fpp_s', 'Probe spacing', 'mm', 1, 0.1, 5, ['fpp'], 'Distance between adjacent collinear probes.'),
    num('fpp_d', 'Sample diameter', 'mm', 50, 5, 300, ['fpp'], 'Lateral size of the (circular) sample for the finite-size correction; the correction is negligible above about 40 probe spacings.'),
    tab('hall_table', 'Hall data (B, V_H)', ['B_T', 'V_hall_mV'], HALL_DEF, 'hall', 'Hall voltage against magnetic field at constant current, both field polarities. Default: simulated n-type doped ZnO film (3×10²⁰ cm⁻³) with a small misalignment offset.'),
    num('hall_I', 'Hall measurement current', 'mA', 1, 0.001, 100, ['hall'], 'Constant current during the Hall sweep.', { log: true }),
    num('hall_RA', 'van der Pauw resistance R_A', 'Ω', 12.0, 1e-4, 1e6, ['hall'], 'Four-terminal resistance V₄₃/I₁₂ along one edge of the sample.', { log: true }),
    num('hall_RB', 'van der Pauw resistance R_B', 'Ω', 12.5, 1e-4, 1e6, ['hall'], 'Four-terminal resistance V₁₄/I₂₃ along the perpendicular edge.', { log: true }),
    tab('eis_table', 'Impedance spectrum (f, Z′, Z″)', ['f_Hz', 'Z_re_ohm', 'Z_im_ohm'], EIS_DEF, 'eis', 'Frequency, real and imaginary impedance (capacitive Z″ negative). Default: simulated Randles cell with CPE and Warburg diffusion plus 0.6 % noise.'),
    tab('flash_table', 'Laser-flash rear-face signal (t, V)', ['time_ms', 'signal'], FLASH_DEF, 'flash', 'Detector signal of the rear face after the flash at t = 0. Default: adiabatic Parker curve for α = 1.8×10⁻⁵ m²/s through 1 mm with noise.'),
    num('fl_L', 'Sample thickness', 'mm', 1, 0.05, 10, ['flash'], 'Thickness of the disc along the heat-flow direction; enters the diffusivity squared.'),
    num('fl_cp', 'Specific heat capacity', 'J/kg/K', 494, 50, 5000, ['flash'], 'Heat capacity at the measurement temperature for the conversion to thermal conductivity.', { bus: 'prop.cp', group: 'Sample' }),
    tab('tens_table', 'Stress–strain curve (ε, σ, optional ε_transverse)', ['strain', 'stress_MPa', 'transverse_strain'], TENS_DEF, 'tensile', 'Engineering strain, stress and optionally transverse strain. Default: simulated aluminium-alloy-like curve (E = 69 GPa, proof stress 276 MPa, ν = 0.33).'),
    num('tens_offset', 'Offset strain for proof stress', '%', 0.2, 0.01, 1, ['tensile'], 'Plastic strain offset of the proof-stress construction (0.2 % by convention for metals).'),
    tab('vsm_table', 'Hysteresis loop (H, M)', ['H_kA_m', 'M_kA_m'], VSM_DEF, 'vsm', 'Applied field and volume magnetisation around a full loop, in acquisition order. Default: simulated ferrimagnetic nanoparticle powder (M_s = 350 kA/m, H_c = 8 kA/m).'),
    { key: 'gum_red', label: 'Measurement model', type: 'select', value: 'scherrer', options: Object.entries(GUM).map(([value, g]) => ({ value, label: g.label })), group: 'Analysis settings', showIf: { mode: ['gum'] }, help: 'Reduction equation whose uncertainty budget is evaluated. The budget table must name its input quantities; missing ones fall back to typical values.' },
    { key: 'gum_table', label: 'Budget inputs (quantity, estimate, standard uncertainty)', type: 'table', columns: ['quantity', 'value', 'u_standard'], value: GUM_DEF, group: 'Measured data', showIf: { mode: ['gum'] }, help: `One row per input quantity with its estimate and standard uncertainty. Names per model — ${Object.entries(GUM).map(([k, g]) => `${k}: ${g.vars.map(v => v[0]).join(', ')}`).join('; ')}.` },
    num('gum_k', 'Coverage factor k', '', 2, 1, 3, ['gum'], 'Multiplier of the combined standard uncertainty for the expanded uncertainty (2 for about 95 % coverage).'),
    tab('rep_table', 'Replicate measurements', ['value'], REP_DEF, 'repeat', 'Repeated results of the same measurand under repeatability conditions, one per row. Default: ten simulated readings of a 100-unit reference artefact.'),
    num('rep_ref', 'Reference value', '', 100, 0, 1e6, ['repeat'], 'Certified or accepted value of the reference artefact in the same unit as the replicates.'),
    num('rep_Uref', 'Expanded uncertainty of the reference (k = 2)', '', 1.2, 0, 1e5, ['repeat'], 'Expanded uncertainty stated on the reference certificate.'),
    num('rep_uB', 'Instrument type-B standard uncertainty', '', 0.4, 0, 1e5, ['repeat'], 'Standard uncertainty from calibration, resolution and environment that repetition does not reduce.'),
    num('rep_target', 'Target precision of the mean', '%', 1, 0.01, 50, ['repeat'], 'Desired 95 % half-width from type-A scatter, used to compute the number of replicates needed.'),
    { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Analysis', help: 'Seed of the Monte Carlo cross-check in the uncertainty-budget mode; identical inputs and seed give identical results.' },
  ],
  presets: [
    { name: 'ZnO nanorods — XRD phase, lattice and crystallite size', note: 'Wurtzite ZnO pattern with 25 nm crystallites: phase ranking, a and c refinement, Scherrer and Williamson–Hall analysis.', values: { mode: 'xrd', xrd_table: XRD_ZNO } },
    { name: 'Au nanoparticles — XRD (Au/Ag ambiguity)', note: 'A 12 nm gold pattern: gold and silver are isostructural with nearly equal cells, so the ranking shows a genuinely ambiguous identification.', values: { mode: 'xrd', xrd_table: XRD_AU } },
    { name: 'Au nanoparticle colloid — LSPR and Beer–Lambert', note: 'Extinction spectrum of a citrate gold sol: peak position, empirical size estimate and particle concentration.', values: { mode: 'lspr', lspr_metal: 'Au' } },
    { name: 'Auxetic re-entrant lattice coupon — tensile test', note: 'Polymer re-entrant honeycomb with transverse strain: effective modulus, proof stress and a negative Poisson’s ratio.', values: { mode: 'tensile', tens_table: TENS_AUX } },
    { name: 'Fused-silica reference — Oliver–Pharr nanoindentation', note: 'Standard calibration material: modulus near 72 GPa and hardness near 9 GPa from the unloading curve.', values: { mode: 'op', op_nu: 0.17 } },
    { name: 'ZnO nanorod film — Raman modes', note: 'E2(high) band of wurtzite ZnO near 437 cm⁻¹ with weaker defect-related bands; shows the ZnO/rutile assignment ambiguity.', values: { mode: 'raman', raman_table: RAMAN_ZNO } },
    { name: 'Au nanoparticles — DLS z-average', note: 'Cumulant analysis of a 24 nm gold sol (hydrodynamic, intensity-weighted size).', values: { mode: 'dls', dls_table: DLS_AU } },
    { name: 'Au nanoparticles — counted TEM diameters', note: 'Number statistics of 150 counted gold particles: means, log-normal fit and the count needed for a target precision.', values: { mode: 'counts', tem_table: TEM_AU } },
  ],
  run(inp) {
    const o = makeOut(), mode = MODES[inp.mode] ? inp.mode : 'xrd';
    MODES[mode][1](inp, o);
    if (!o.checks.length) o.C('Data table parsed', true, 'numeric rows were read');
    o.R('medium', 'Record sample identification, preparation history, instrument, calibration date and environmental conditions with this reduction so that the result can be traced and repeated.', 'Metrological traceability of characterisation records.');
    if (!o.recommendations.some(r => r.sustainability)) o.R('low', 'Archive the raw data with this analysis and re-reduce existing files when settings change instead of re-measuring; batch samples per instrument session to cut standby energy, cryogen and consumable use.', 'Resource efficiency of laboratory characterisation.', true);
    const rep = [`Suite 23 characterisation report — ${MODES[mode][0]}`, '', o.summary, '', 'Results:', ...o.kpis.map(k => `  ${k.label}: ${typeof k.value === 'number' ? F(k.value, 5) : k.value} ${k.unit || ''}`), '', 'Checks:', ...o.checks.map(c => `  [${c.pass ? 'pass' : 'FAIL'}] ${c.name}: ${c.detail}`), '', 'Limits of validity:', ...o.warnings.map(w => '  - ' + w), ...(o.synthetic ? [`  - Input table is ${SYN_MSG}.`] : [])].join('\n');
    o.exports.push({ name: `characterisation_${mode}.csv`, mime: 'text/csv', content: 'key,label,value,unit\n' + o.kpis.map(k => `${k.key},"${k.label.replace(/"/g, "'")}",${typeof k.value === 'number' ? k.value : `"${k.value}"`},${k.unit || ''}`).join('\n'), description: 'Structured property record (key results) for simulation calibration.' }, { name: `characterisation_${mode}.txt`, mime: 'text/plain', content: rep, description: 'Plain-text characterisation report with results, checks and validity limits.' });
    return sealOut(o);
  },
};
