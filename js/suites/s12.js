// Suite 12 — Electromagnetic, Photonic and Optical Engineering.
// Complex transfer-matrix multilayers and 1-D photonic crystals, Drude/wire negative permittivity, split-ring negative
// permeability, negative-index composites with S-parameter retrieval, Mie/Gans nanoparticle optics, effective-medium and
// hyperbolic multilayers, and phase-gradient metasurfaces / metalenses.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { PI, C0, EPS0, MU0, C } = N;
const f3 = v => N.fmt(v, 3), f4 = v => N.fmt(v, 4);
const cx = z => `${f3(z.re)} ${z.im < 0 ? '−' : '+'} ${f3(Math.abs(z.im))}i`;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ser = (key, title, xlabel, ylabel, lines, explain, o = {}) => ({ key, title, xlabel, ylabel, logx: !!o.logx, logy: !!o.logy, lines, explain });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const pos = v => (v > 1e-300 ? v : 1e-300);
function finish(o) {
  o.warnings = o.warnings || []; o.tables = o.tables || []; let bad = 0;
  for (const k of o.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; bad++; }
  for (const s of o.series) for (const l of s.lines) {
    if (l.x.length > 400) { const st = (l.x.length - 1) / 399, ix = Array.from({ length: 400 }, (_, i) => Math.round(i * st)); l.x = ix.map(i => l.x[i]); l.y = ix.map(i => l.y[i]); }
    l.x = Array.from(l.x); l.y = Array.from(l.y, v => { if (!Number.isFinite(v)) { bad++; return 0; } return v; });
  }
  for (const k of Object.keys(o.publish)) if (typeof o.publish[k] === 'number' && !Number.isFinite(o.publish[k])) delete o.publish[k];
  if (bad) o.warnings.push('Some outputs were not computable for this input combination and are reported as 0 (non-finite values suppressed).');
  return o;
}
const up = z => (z.im < 0 || (z.im === 0 && z.re < 0) ? { re: -z.re, im: -z.im } : z); // passive branch, Im ≥ 0
const csqrt = z => up(C.sqrt(z)), ONE = { re: 1, im: 0 }, J = { re: 0, im: 1 };
const band = (xs, flag) => { // widest contiguous run where flag is true → [lo, hi] or null
  let best = null, s = -1; for (let i = 0; i <= xs.length; i++) { const on = i < xs.length && flag[i]; if (on && s < 0) s = i; if (!on && s >= 0) { if (!best || xs[i - 1] - xs[s] > best[1] - best[0]) best = [xs[s], xs[i - 1]]; s = -1; } } return best;
};

// ---------------------------------------------------------------- dispersive materials (λ in nm → ε = (n + ik)²)
const LD = { // Lorentz–Drude fits (Rakić et al. 1998 parameter set), energies in eV
  Au: { wp: 9.03, f: [0.760, 0.024, 0.010, 0.071, 0.601, 4.384], G: [0.053, 0.241, 0.345, 0.870, 2.494, 2.214], w: [0, 0.415, 0.830, 2.969, 4.304, 13.32] },
  Ag: { wp: 9.01, f: [0.845, 0.065, 0.124, 0.011, 0.840, 5.646], G: [0.048, 3.886, 0.452, 0.065, 0.916, 2.419], w: [0, 0.816, 4.481, 8.185, 9.083, 20.29] },
  Al: { wp: 14.98, f: [0.523, 0.227, 0.050, 0.166, 0.030], G: [0.047, 0.333, 0.312, 1.351, 3.382], w: [0, 0.162, 1.544, 1.808, 3.473] },
};
// Measured optical constants of silver, [photon energy (eV), n, k]: Johnson & Christy, "Optical constants of the noble metals",
// Phys. Rev. B 6, 4370 (1972), room-temperature evaporated films, 0.64–6.60 eV. n and k are interpolated linearly in energy.
const JC_AG = [[0.64,0.24,14.08],[0.77,0.15,11.85],[0.89,0.13,10.1],[1.02,0.09,8.828],[1.14,0.04,7.795],[1.26,0.04,6.992],[1.39,0.04,6.312],[1.51,0.04,5.727],[1.64,0.03,5.242],[1.76,0.04,4.838],[1.88,0.05,4.483],[2.01,0.06,4.152],[2.13,0.05,3.858],[2.26,0.06,3.586],[2.38,0.05,3.324],[2.5,0.05,3.093],[2.63,0.05,2.869],[2.75,0.04,2.657],[2.88,0.04,2.462],[3,0.05,2.275],[3.12,0.05,2.07],[3.25,0.05,1.864],[3.37,0.07,1.657],[3.5,0.1,1.419],[3.62,0.14,1.142],[3.74,0.17,0.829],[3.87,0.81,0.392],[3.99,1.13,0.616],[4.12,1.34,0.964],[4.24,1.39,1.161],[4.36,1.41,1.264],[4.49,1.41,1.331],[4.61,1.38,1.372],[4.74,1.35,1.387],[4.86,1.33,1.393],[4.98,1.31,1.389],[5.11,1.3,1.378],[5.23,1.28,1.367],[5.36,1.28,1.357],[5.48,1.26,1.344],[5.6,1.25,1.342],[5.73,1.22,1.336],[5.85,1.2,1.325],[5.98,1.18,1.312],[6.1,1.15,1.296],[6.22,1.14,1.277],[6.35,1.12,1.255],[6.47,1.1,1.232],[6.6,1.07,1.212]];
const jcAg = E => { let i = 1; while (i < JC_AG.length - 1 && JC_AG[i][0] < E) i++; const a = JC_AG[i - 1], b = JC_AG[i], t = (E - a[0]) / (b[0] - a[0]), n = a[1] + t * (b[1] - a[1]), k = a[2] + t * (b[2] - a[2]); return { re: n * n - k * k, im: 2 * n * k }; };
const sell = (l2, terms, a = 1) => terms.reduce((s, [B, Cc]) => s + B * l2 / (l2 - Cc), a);
const MATS = {
  Au: { label: 'Gold (Lorentz–Drude)', metal: true }, Ag: { label: 'Silver (Johnson & Christy table)', metal: true }, Al: { label: 'Aluminium (Lorentz–Drude)', metal: true },
  ITO: { label: 'Indium tin oxide (Drude)', metal: true, formula: 'In2O3' }, Si: { label: 'Silicon' }, SiO2: { label: 'Fused silica (Sellmeier)' }, TiO2: { label: 'Titanium dioxide, rutile (ordinary)' },
  ZnO: { label: 'Zinc oxide (ordinary)' }, glass: { label: 'Borosilicate crown glass (BK7 Sellmeier)', formula: 'SiO2' }, air: { label: 'Air / vacuum', formula: 'N2' }, custom: { label: 'Custom constant n + ik', formula: null },
};
function epsOf(name, lam, inp) {
  const E = 1239.84193 / lam, um = lam / 1000;
  if (LD[name]) { const m = LD[name], ld = (() => { let e = C.sub(ONE, C.div(C.of(m.f[0] * m.wp * m.wp), C.mul(C.of(E), C.of(E, m.G[0])))); for (let j = 1; j < m.f.length; j++) e = C.add(e, C.div(C.of(m.f[j] * m.wp * m.wp), C.of(m.w[j] * m.w[j] - E * E, -E * m.G[j]))); return e; })();
    // silver: tabulated Johnson & Christy constants inside their range (the Lorentz–Drude fit is too smooth at the 3.8 eV interband edge and
    // places Re ε = −2 at 368 nm instead of 354 nm); blended into the Lorentz–Drude fit over the outermost table interval at each end
    if (name === 'Ag' && E > 0.64 && E < 6.6) { const jc = jcAg(E), w = Math.min(1, (E - 0.64) / 0.13, (6.6 - E) / 0.13); return { re: w * jc.re + (1 - w) * ld.re, im: w * jc.im + (1 - w) * ld.im }; }
    return ld; }
  if (name === 'ITO') return C.sub(C.of(3.9), C.div(C.of(1.9 * 1.9), C.mul(C.of(E), C.of(E, 0.12))));
  let n = 1, k = 0; const u = Math.max(um, 0.4), l2 = u * u;
  if (name === 'Si') { n = Math.sqrt(sell(l2, [[10.6684293, 0.301516485 ** 2], [1.54133408, 1104 ** 2]])); k = E > 1.12 ? 4.6e5 * (E - 1.12) ** 2 * lam * 1e-9 / (4 * PI) : 0; }
  else if (name === 'SiO2') { const s2 = Math.max(um, 0.21) ** 2; n = Math.sqrt(sell(s2, [[0.6961663, 0.0684043 ** 2], [0.4079426, 0.1162414 ** 2], [0.8974794, 9.896161 ** 2]])); }
  else if (name === 'TiO2') n = Math.sqrt(5.913 + 0.2441 / (l2 - 0.0803));
  else if (name === 'ZnO') n = Math.sqrt(2.81418 + 0.87968 * l2 / (l2 - 0.3042 ** 2) - 0.00711 * l2);
  else if (name === 'glass') { const s2 = Math.max(um, 0.3) ** 2; n = Math.sqrt(sell(s2, [[1.03961212, 0.00600069867], [0.231792344, 0.0200179144], [1.01046945, 103.560653]])); }
  else if (name === 'custom') { n = inp.n_user; k = inp.k_user; }
  return { re: n * n - k * k, im: 2 * n * k };
}
const nOf = (name, lam, inp) => csqrt(epsOf(name, lam, inp));
const matFormula = k => (MATS[k].formula === undefined ? k : MATS[k].formula);
function matEquations(keys) { const out = []; for (const k of [...new Set(keys)]) { const f = matFormula(k); if (f) { try { out.push(eqn('formula', `Material: ${MATS[k].label}`, `${chem.prettyFormula(f)}, M = ${f4(chem.molarMass(f))} g/mol`)); } catch { /* formula not parsed */ } } } return out; }
const dispWarn = (keys, lmin, W) => { if (lmin < 400 && keys.some(k => ['Si', 'TiO2', 'ZnO'].includes(k))) W.push('Below 400 nm the Si, TiO₂ and ZnO dispersion formulae are held at their 400 nm values and interband absorption is not included; results in the ultraviolet are indicative only.'); };

// ---------------------------------------------------------------- transfer-matrix kernel
function layerMat(n, d, lam, kx, pol) {
  const q = csqrt(C.sub(C.mul(n, n), C.of(kx * kx))), eta = pol === 'p' ? C.div(C.mul(n, n), q) : q, del = C.scale(q, 2 * PI * d / lam); del.im = Math.min(del.im, 300);
  const c = C.cos(del), s = C.sin(del), is = C.mul({ re: 0, im: -1 }, s);
  return [c, C.div(is, eta), C.mul(is, eta), c];
}
const mul2 = (a, b) => [C.add(C.mul(a[0], b[0]), C.mul(a[1], b[2])), C.add(C.mul(a[0], b[1]), C.mul(a[1], b[3])), C.add(C.mul(a[2], b[0]), C.mul(a[3], b[2])), C.add(C.mul(a[2], b[1]), C.mul(a[3], b[3]))];
function tmm(ns, ds, n0, nsub, lam, th, pol, field) {
  const kx = n0 * Math.sin(th), qs = csqrt(C.sub(C.mul(nsub, nsub), C.of(kx * kx))), es = pol === 'p' ? C.div(C.mul(nsub, nsub), qs) : qs, e0 = pol === 'p' ? n0 / Math.cos(th) : n0 * Math.cos(th);
  let B = ONE, Cc = es, L = 0; const fl = [];
  for (let j = ns.length - 1; j >= 0; j--) {
    const m = layerMat(ns[j], ds[j], lam, kx, pol), B2 = C.add(C.mul(m[0], B), C.mul(m[1], Cc)), C2 = C.add(C.mul(m[2], B), C.mul(m[3], Cc)); B = B2; Cc = C2;
    const mg = Math.max(C.abs(B), C.abs(Cc)); if (mg > 1e50) { B = C.scale(B, 1 / mg); Cc = C.scale(Cc, 1 / mg); L += Math.log(mg); }
    if (field) fl.push([C.abs(B), L]);
  }
  const den = C.add(C.scale(B, e0), Cc), r = C.div(C.sub(C.scale(B, e0), Cc), den), R = Math.min(C.abs2(r), 1), T = Math.min(4 * e0 * Math.max(es.re, 0) / pos(C.abs2(den)) * Math.exp(-2 * L), 1 - R);
  const out = { R, T, A: Math.max(1 - R - T, 0), r };
  if (field) out.E2 = fl.reverse().map(([b, l]) => (b * 2 * e0 / pos(C.abs(den)) * Math.exp(l - L)) ** 2);
  return out;
}
const tmmPol = (ns, ds, n0, nsub, lam, th, pol) => { if (pol !== 'u') return tmm(ns, ds, n0, nsub, lam, th, pol); const a = tmm(ns, ds, n0, nsub, lam, th, 's'), b = tmm(ns, ds, n0, nsub, lam, th, 'p'); return { R: (a.R + b.R) / 2, T: (a.T + b.T) / 2, A: (a.A + b.A) / 2 }; };
const cacos = z => { if (!(C.abs(z) < 1e7)) return { re: Math.abs(Math.atan2(z.im, z.re)), im: Math.log(Math.min(2 * C.abs(z), 1e300)) }; const u = C.add(z, C.mul(J, C.sqrt(C.sub(ONE, C.mul(z, z))))); return { re: Math.abs(Math.atan2(u.im, u.re)), im: Math.abs(Math.log(pos(C.abs(u)))) }; };

function multilayer(inp) {
  const W = [], l0 = inp.lam0, th = Math.min(inp.theta, 89) * PI / 180, n0 = inp.n_inc, pol = inp.pol, H = inp.matH, Lm = inp.matL, Np = Math.round(inp.Nper);
  const nH0 = nOf(H, l0, inp), nL0 = nOf(Lm, l0, inp), dH = inp.dH > 0 ? inp.dH : l0 / (4 * Math.max(nH0.re, 0.05)), dL = inp.dL > 0 ? inp.dL : l0 / (4 * Math.max(nL0.re, 0.05));
  const seq = []; // [material, thickness]
  if (inp.stack === 'film') seq.push([H, dH]);
  else if (inp.stack === 'cavity') { for (let i = 0; i < Np; i++) seq.push([H, dH], [Lm, dL]); seq.push([Lm, dL]); for (let i = 0; i < Np; i++) seq.push([Lm, dL], [H, dH]); }
  else for (let i = 0; i < Np; i++) seq.push([H, dH], [Lm, dL]);
  const build = lam => ({ ns: seq.map(s => nOf(s[0], lam, inp)), ds: seq.map(s => s[1]) }), lmin = Math.min(inp.lam_min, inp.lam_max - 10), lams = N.linspace(lmin, inp.lam_max, 321), pols = pol === 'u' ? 's' : pol;
  dispWarn([H, Lm, inp.substrate], lmin, W);
  const spec = lams.map(l => { const b = build(l); return tmmPol(b.ns, b.ds, n0, nOf(inp.substrate, l, inp), l, th, pol); });
  const at0 = (() => { const b = build(l0); return tmmPol(b.ns, b.ds, n0, nOf(inp.substrate, l0, inp), l0, th, pol); })();
  // Bloch wavenumber of the infinite periodic stack
  const Lam = dH + dL, kx = n0 * Math.sin(th), bloch = l => { const m = mul2(layerMat(nOf(H, l, inp), dH, l, kx, pols), layerMat(nOf(Lm, l, inp), dL, l, kx, pols)), hf = C.scale(C.add(m[0], m[3]), 0.5); return { K: cacos(hf), tr: hf.re }; };
  const bl = lams.map(bloch), periodic = inp.stack !== 'film';
  let gap = null;
  if (periodic) { const flag = bl.map(b => Math.abs(b.tr) > 1); let i0 = N.argmin(lams.map(l => Math.abs(l - l0))); if (!flag[i0]) { const g = band(lams, flag); if (g) i0 = lams.indexOf(g[0]); } if (flag[i0]) { let a = i0, b = i0; while (a > 0 && flag[a - 1]) a--; while (b < lams.length - 1 && flag[b + 1]) b++; const edge = (i, j) => (j < 0 || j >= lams.length ? lams[i] : N.brent(l => Math.abs(bloch(l).tr) - 1, lams[i], lams[j], 1e-6)); try { gap = [edge(a, a - 1), edge(b, b + 1)]; } catch { gap = [lams[a], lams[b]]; } } }
  const K0 = bloch(l0).K, pen = periodic && K0.im > 1e-6 ? Lam / (2 * K0.im) : 0, fc = gap ? 0.5 * (C0 / (gap[0] * 1e-9) + C0 / (gap[1] * 1e-9)) : C0 / (l0 * 1e-9);
  const nHr = nH0.re, nLr = nL0.re, bwAn = (4 / PI) * Math.asin(Math.abs(nHr - nLr) / (nHr + nLr)), nsub0 = nOf(inp.substrate, l0, inp);
  const qw = inp.stack === 'bragg' && inp.dH === 0 && inp.dL === 0 && th === 0 && nH0.im < 1e-6 && nL0.im < 1e-6 && nsub0.im < 1e-6;
  // analytic quarter-wave reflectance: Y = (n_H/n_L)^(2N) n_s
  const Yq = (nHr / nLr) ** (2 * Np) * nsub0.re, Ran = ((n0 - Yq) / (n0 + Yq)) ** 2;
  // field profile at λ0 (sub-sliced)
  const sub = 6, fs = build(l0), nsF = [], dsF = [], zs = [0]; fs.ns.forEach((n, i) => { for (let k = 0; k < sub; k++) { nsF.push(n); dsF.push(fs.ds[i] / sub); zs.push(zs[zs.length - 1] + fs.ds[i] / sub); } });
  const fld = tmm(nsF, dsF, n0, nsub0, l0, th, pols, true), E2 = [...fld.E2, fld.E2.length ? fld.E2[fld.E2.length - 1] : 1].slice(0, zs.length);
  // cavity mode
  let Qc = 0, lc = l0; if (inp.stack === 'cavity' && gap) { const i1 = lams.findIndex(l => l >= gap[0]), i2 = lams.findIndex(l => l > gap[1]), seg = spec.slice(i1, i2 < 0 ? undefined : i2).map(s => s.T); if (seg.length > 2) { const ip = i1 + N.argmax(seg), Tf = l => { const b = build(l); return tmmPol(b.ns, b.ds, n0, nOf(inp.substrate, l, inp), l, th, pol).T; }; lc = N.goldenSection(l => -Tf(l), lams[Math.max(ip - 1, 0)], lams[Math.min(ip + 1, lams.length - 1)], 1e-7); const Tp = Tf(lc); let hw = 1e-6; for (let k = 0; k < 60 && Tf(lc + hw) > Tp / 2; k++) hw *= 1.5; Qc = lc / (2 * hw); } }
  const Rmax = N.maxOf(spec.map(s => s.R)), maxA = N.maxOf(spec.map(s => s.A)), cons = N.maxOf(spec.map(s => Math.abs(s.R + s.T + s.A - 1)));
  if (inp.theta > 89) W.push('Incidence angle limited to 89° (grazing incidence is singular for the plane-wave formulation).');
  if (periodic && !gap) W.push('No photonic stop band was found in the wavelength window: the index contrast is too low or the window does not contain the Bragg wavelength.');
  if (seq.some(s => MATS[s[0]].metal) && inp.stack !== 'film') W.push('Periodic metal–dielectric stacks are strongly absorbing; the Bloch stop band of a lossy stack is not sharply defined.');
  const fTHz = lams.map(l => C0 / (l * 1e-9) / 1e12);
  return {
    summary: `${inp.stack === 'film' ? `A ${f3(dH)} nm ${H} film` : inp.stack === 'cavity' ? `A ${H}/${Lm} microcavity with ${Np} mirror pairs per side` : `A ${Np}-period ${H}/${Lm} stack (${f3(dH)}/${f3(dL)} nm)`} on ${inp.substrate} reflects ${f3(100 * at0.R)} %, transmits ${f3(100 * at0.T)} % and absorbs ${f3(100 * at0.A)} % at ${l0} nm (${inp.theta}°, ${pol === 'u' ? 'unpolarised' : pol + '-polarised'}). ${gap ? `The stop band of the periodic stack spans ${f3(gap[0])}–${f3(gap[1])} nm (mid-gap frequency at ${f3(2 * gap[0] * gap[1] / (gap[0] + gap[1]))} nm).` : 'No stop band lies in the plotted window.'}${Qc ? ` The cavity mode at ${f4(lc)} nm has Q ≈ ${f3(Qc)}.` : ''}`,
    kpis: [
      kpi('R0', 'Reflectance at design wavelength', at0.R, '–', `Fraction of incident power reflected at ${l0} nm from the complex transfer-matrix solution for the finite stack on its substrate.`),
      kpi('T0', 'Transmittance at design wavelength', at0.T, '–', 'Fraction of incident power delivered into the substrate (back-surface reflection of the substrate is not included).'),
      kpi('A0', 'Absorptance at design wavelength', at0.A, '–', 'Power absorbed in the layers, 1 − R − T; zero for lossless dielectrics and large for thin metal films.'),
      kpi('Rmax', 'Peak reflectance in window', Rmax, '–', 'Highest reflectance over the plotted wavelength range; for a Bragg mirror it approaches 1 as the number of periods grows.'),
      kpi('gap_lo', 'Stop-band lower edge', gap ? gap[0] : 0, 'nm', 'Short-wavelength edge of the photonic band gap of the infinite periodic stack, where |½ tr M| crosses 1 (zero if no gap).'),
      kpi('gap_hi', 'Stop-band upper edge', gap ? gap[1] : 0, 'nm', 'Long-wavelength edge of the photonic band gap from the Bloch condition cos KΛ = ½ tr M.'),
      kpi('gap_c', 'Stop-band centre', gap ? 2 * gap[0] * gap[1] / (gap[0] + gap[1]) : 0, 'nm', `Wavelength of the mid-gap frequency; for a quarter-wave design at normal incidence it sits close to the design wavelength ${l0} nm and blue-shifts with angle.`),
      kpi('bw', 'Fractional bandwidth', gap ? 2 * (gap[1] - gap[0]) / (gap[1] + gap[0]) : 0, '–', `Gap width relative to centre frequency, Δω/ω₀; the quarter-wave formula (4/π)·asin((n_H − n_L)/(n_H + n_L)) gives ${f3(bwAn)}.`),
      kpi('pen', 'Field penetration depth', pen, 'nm', 'Intensity 1/e decay length of the Bloch wave inside the stop band, Λ/(2 Im KΛ); it sets how many periods are needed and the effective cavity length.'),
      kpi('dH', 'High-index layer thickness', dH, 'nm', `Physical thickness of the ${H} layer (quarter-wave at ${l0} nm when the input is zero).`),
      kpi('dL', 'Low-index layer thickness', dL, 'nm', `Physical thickness of the ${Lm} layer (quarter-wave at ${l0} nm when the input is zero).`),
      kpi('Q', 'Cavity quality factor', Qc, '–', 'λ/Δλ of the defect-mode transmission peak (zero unless the microcavity design is selected).'),
    ],
    series: [
      ser('rta', 'Reflectance, transmittance and absorptance spectra', 'wavelength (nm)', 'fraction (–)', [{ name: 'R', x: lams, y: spec.map(s => s.R) }, { name: 'T', x: lams, y: spec.map(s => s.T) }, { name: 'A', x: lams, y: spec.map(s => s.A) }], 'Power coefficients of the finite stack. A flat-topped reflectance plateau is the photonic stop band; the ripples outside it are Fabry–Perot fringes of the whole stack.'),
      ser('band', 'Bloch band diagram of the periodic stack', 'K·Λ/π (–)', 'frequency (THz)', [{ name: 'Re K (propagating)', x: bl.map(b => b.K.re / PI), y: fTHz }, { name: 'Im K (evanescent)', x: bl.map(b => b.K.im / PI), y: fTHz }], 'Dispersion of the infinite periodic stack from cos KΛ = ½ tr M. Where Re K is pinned at the zone edge and Im K is non-zero the wave is evanescent: this is the band gap.'),
      ser('field', 'Field intensity profile at the design wavelength', 'depth from the top surface (nm)', '|E_t|²/|E_inc|² (–)', [{ name: '|E|²', x: zs, y: E2 }], 'Tangential electric-field intensity through the stack at the design wavelength. In a stop band the envelope decays exponentially with the penetration depth; in a cavity it peaks in the spacer.'),
    ],
    equations: [...matEquations([H, Lm, inp.substrate]),
      eqn('governing', 'Helmholtz equation in each layer', 'd²E/dz² + (k₀² n_j² − k_x²) E = 0'), eqn('governing', 'Characteristic matrix', 'M_j = [[cos δ_j, −i sin δ_j/η_j], [−i η_j sin δ_j, cos δ_j]],  δ_j = 2π n_j d_j cos θ_j/λ'),
      eqn('formula', 'Fresnel coefficients of the stack', 'r = (η₀B − C)/(η₀B + C),  T = 4η₀ Re(η_s)/|η₀B + C|²,  [B; C] = ΠM_j [1; η_s]'), eqn('governing', 'Bloch condition', 'cos(KΛ) = ½ (M₁₁ + M₂₂)'),
      eqn('formula', 'Quarter-wave mirror', 'R = [(n₀ − Y)/(n₀ + Y)]²,  Y = (n_H/n_L)^(2N) n_s;  Δω/ω₀ = (4/π) asin[(n_H − n_L)/(n_H + n_L)]')],
    checks: [
      chk('Energy conservation R + T + A = 1', cons < 1e-9, `max residual ${f3(cons)}`), chk('Passivity: absorptance non-negative', spec.every(s => s.A >= 0), `max A = ${f3(maxA)}`),
      chk('Analytical quarter-wave reflectance reproduced', !qw || Math.abs(at0.R - Ran) < 2e-3, qw ? `TMM ${f4(at0.R)} against closed form ${f4(Ran)}` : 'not a lossless quarter-wave stack at normal incidence; check skipped'),
      chk('Stop-band width consistent with quarter-wave formula', !qw || !gap || Math.abs(2 * (1 / gap[0] - 1 / gap[1]) / (1 / gap[0] + 1 / gap[1]) - bwAn) < 0.15 * bwAn + 0.01, gap ? `Bloch Δω/ω₀ = ${f3(2 * (1 / gap[0] - 1 / gap[1]) / (1 / gap[0] + 1 / gap[1]))}, formula ${f3(bwAn)} (index dispersion causes the difference)` : 'no gap'),
    ],
    warnings: W,
    recommendations: [
      rec('high', periodic ? `With ${Np} periods the mirror reaches R = ${f4(at0.R)}; each added pair multiplies the transmitted leakage by about (n_L/n_H)² = ${f3((nLr / nHr) ** 2)}. ${at0.R > 0.99 ? 'Further pairs add stress and deposition time for little gain.' : 'Add pairs or raise the index contrast if higher reflectance is needed.'}` : `For an antireflection film choose n ≈ √(n₀·n_s) = ${f3(Math.sqrt(n0 * nsub0.re))} at quarter-wave thickness; the present film gives R = ${f4(at0.R)}.`, 'Quarter-wave stack scaling'),
      rec('medium', `Hold layer thickness to about ±1 % (±${f3(0.01 * dH)} nm for ${H}): the stop-band centre shifts in proportion to optical thickness, i.e. ±${f3(0.01 * l0)} nm.`, 'Linear dependence of Bragg wavelength on optical thickness'),
      rec('medium', `Design for the working angle: at ${inp.theta}° the band centre is at ${gap ? f3(2 * gap[0] * gap[1] / (gap[0] + gap[1])) : '—'} nm; s- and p-polarised gaps separate with angle, so specify polarisation.`, 'Angular blue shift of Bragg reflection'),
      rec('low', 'Measure the real n and k of deposited films by ellipsometry; sputtered and sol–gel oxides are typically 2–10 % less dense than the handbook dispersion used here.', 'Built-in dispersion models are bulk reference data'),
      rec('medium', `Dielectric mirrors (A = ${f3(100 * at0.A)} %) avoid the ohmic loss of metal mirrors and need no precious metals; prefer abundant oxides (SiO₂, TiO₂, ZnO) and minimise the pair count to cut deposition energy.`, 'Material and energy efficiency of coating design', true),
    ],
    publish: { 'res.transmission': at0.T, 'res.absorption': at0.A, 'res.f0': fc, ...(gap ? { 'res.bandgap_lo': C0 / (gap[1] * 1e-9), 'res.bandgap_hi': C0 / (gap[0] * 1e-9) } : {}), 'geom.thickness': zs[zs.length - 1] * 1e-9, 'prop.n_refr': nHr },
    numerics: { method: 'Complex characteristic (transfer) matrix method with rescaling; Bloch eigenvalue of the period matrix', dof: seq.length, steps: lams.length, notes: 'Plane waves, planar homogeneous isotropic layers, coherent interference; semi-infinite substrate.' },
  };
}

// ---------------------------------------------------------------- Kramers–Kronig numerical check
function kkMismatch(chi, w, ws) { // Re χ(ω) = (2/π) P∫ ω′ Im χ/(ω′² − ω²) dω′, singularity-subtracted, tan-mapped grid
  const n = 20000, a = w * chi(w).im; let s = 0;
  for (let i = 0; i < n; i++) { const t = (i + 0.5) * PI / (2 * n), wp = ws * Math.tan(t), d = wp * wp - w * w; if (Math.abs(d) < 1e-12 * w * w) continue; s += (wp * chi(wp).im - a) / d * ws / (Math.cos(t) ** 2); }
  const re = (2 / PI) * s * PI / (2 * n), ex = chi(w).re; return Math.abs(re - ex) / Math.max(Math.abs(ex), 1e-12);
}
const kkCheck = (name, chi, w, ws, gRel) => (gRel < 3e-4 ? chk(name, true, 'resonance too narrow for the numerical quadrature; causality follows analytically because all poles lie in the lower half-plane (damping > 0)') : (() => { const m = kkMismatch(chi, w, ws); return chk(name, m < 0.03, `numerical Hilbert transform of Im χ reproduces Re χ within ${f3(100 * m)} % at the probe frequency`); })());
const fUnit = f => (f >= 3e11 ? [1e12, 'THz'] : f >= 3e8 ? [1e9, 'GHz'] : [1e6, 'MHz']);

// ---------------------------------------------------------------- negative permittivity
const DRUDE = { Au: { einf: 9.84, wp: 9.01, g: 0.072 }, Ag: { einf: 3.7, wp: 9.2, g: 0.021 }, Al: { einf: 1, wp: 14.75, g: 0.0818 }, ITO: { einf: 3.9, wp: 1.9, g: 0.12 } };
const EV = 1.602176634e-19 / N.HBAR; // rad/s per eV
// Plasma (cut-off) frequency of a square lattice of parallel conducting wires, E along the wires.
//  'pendry'  : ω_p² = 2πc²/(a² ln(a/r))                      — original logarithmic thin-wire estimate
//  'belov'   : ω_p² = 2πc²/(a² [ln(a/(2πr)) + 0.5275])        — quasi-static lattice sum (Belov, Tretyakov & Maslovski), r ≪ a and ω_p a/c ≪ 1
//  'lattice' : dynamic lattice sum, 1/κ² = Σ_{G≠0} J₀²(Gr)/(G² − κ²), κ = ω_p/c — zero of the cell-periodic field of a uniform surface current
//              on each wire, averaged over the wire surface (monopole order); valid for thick wires (r/a up to about 0.2) and it reduces to
//              the 'belov' expression when κa ≪ 1. No fitted constant.
const WIRE_SUM = {};
function wirePlasma(a, r, model) {
  const ra = r / a, pend = Math.sqrt(2 * PI / Math.log(1 / ra)) / a, belArg = Math.log(1 / (2 * PI * ra)) + 0.5275;
  if (model === 'pendry') return { wp: C0 * pend, model };
  if (model === 'belov') return belArg > 0.05 ? { wp: C0 * Math.sqrt(2 * PI / belArg) / a, model } : { wp: C0 * pend, model: 'pendry', fallback: true };
  // dimensionless: lengths in units of a, G = 2π(m, n)
  const key = ra.toPrecision(6); if (WIRE_SUM[key]) return { wp: C0 * WIRE_SUM[key] / a, model };
  const Mc = Math.min(60, Math.max(16, Math.ceil(6 / (2 * PI * ra)))), Gc = 2 * PI * Mc, G2 = [], J2 = [], mult = [];
  for (let m = 1; m <= Mc; m++) for (let n = 0; n <= m; n++) { const g = 2 * PI * Math.hypot(m, n); if (g > Gc) continue; G2.push(g * g); const j = N.besselJ(0, g * ra); J2.push(j * j); mult.push(n === 0 || n === m ? 4 : 8); }
  // static part Σ J₀²/G²: numerical with an integral tail when the cut-off resolves the wire (G_c r ≥ 6), otherwise the thin-wire closed form
  let stat; if (Gc * ra >= 5.9) { stat = 1 / (2 * PI * PI * ra * Gc); for (let i = 0; i < G2.length; i++) stat += mult[i] * J2[i] / G2[i]; } else stat = belArg / (2 * PI);
  const f = k2 => { let d = 0; for (let i = 0; i < G2.length; i++) d += mult[i] * J2[i] * k2 / (G2[i] * (G2[i] - k2)); return 1 / k2 - stat - d; };
  let lo = 1e-8, hi = 4 * PI * PI * (1 - 1e-9); for (let i = 0; i < 70; i++) { const m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; }
  const ka = Math.sqrt(0.5 * (lo + hi)); WIRE_SUM[key] = ka; return { wp: C0 * ka / a, model };
}
function wireParams(inp) { const a = inp.wire_a * 1e-3, r = Math.min(inp.wire_r * 1e-6, a / 4), q = wirePlasma(a, r, ['pendry', 'belov'].includes(inp.wire_model) ? inp.wire_model : 'lattice'), wp = q.wp; return { a, r, wp, model: q.model, fallback: !!q.fallback, wpPendry: wirePlasma(a, r, 'pendry').wp, g: EPS0 * a * a * wp * wp / (PI * r * r * inp.metal_sigma) }; }
const WIRE_MODEL = { pendry: ['Pendry\'s thin-wire formula, ω_p² = 2πc²/(a² ln(a/r))', 'ω_p² = 2πc₀²/(a² ln(a/r))'], belov: ['the quasi-static lattice-sum formula, ω_p² = 2πc²/(a² [ln(a/(2πr)) + 0.5275])', 'ω_p² = 2πc₀²/(a² [ln(a/(2πr)) + 0.5275])'], lattice: ['the dynamic lattice sum 1/κ² = Σ J₀²(Gr)/(G² − κ²), κ = ω_p/c, valid also for thick wires', '(c₀/ω_p)² = Σ_{G≠0} J₀²(Gr)/(G² − ω_p²/c₀²),  G = 2π(m, n)/a'] };
function negEps(inp) {
  const W = [], wire = inp.medium === 'wire'; let einf, wp, g, wp_ = null;
  if (wire) { wp_ = wireParams(inp); einf = 1; wp = wp_.wp; g = wp_.g; if (inp.wire_r * 1e-6 > wp_.a / 4) W.push('Wire radius limited to a quarter of the period.');
    const ra = wp_.r / wp_.a;
    if (wp_.fallback) W.push('The quasi-static formula has no solution for wires this thick (ln(a/2πr) + 0.5275 ≤ 0); Pendry\'s expression was used. Select the dynamic lattice sum.');
    if (wp_.model === 'pendry' && ra > 0.01) W.push(`Pendry's logarithmic formula is a thin-wire estimate: at r/a = ${f3(ra)} it lies ${f3(100 * (1 - wp_.wpPendry / wirePlasma(wp_.a, wp_.r, 'lattice').wp))} % below the lattice-sum plasma frequency.`);
    if (wp_.model === 'belov' && ra > 0.05) W.push(`The quasi-static formula assumes r ≪ a and ω_p a/c ≪ 1; at r/a = ${f3(ra)} it lies ${f3(100 * (wp / wirePlasma(wp_.a, wp_.r, 'lattice').wp - 1))} % above the dynamic lattice sum.`);
    if (ra > 0.2) W.push('For r/a above about 0.2 the wires are no longer thin on the scale of the cell: higher multipoles of the wire current matter and the monopole lattice sum is approximate (a full-wave band calculation is needed).');
    const dsk = Math.sqrt(2 / (wp * N.MU0 * inp.metal_sigma)); if (dsk < 0.5 * wp_.r) W.push(`The skin depth at the plasma frequency (${f3(dsk * 1e6)} µm) is well below the wire radius: the current flows in a surface layer, so the d.c.-resistance damping rate underestimates the loss by a factor of about ${f3(wp_.r / (2 * dsk))}.`); }
  else { const d = DRUDE[inp.medium] || DRUDE.Ag; einf = d.einf; wp = d.wp * EV; g = d.g * EV; }
  const eps = w => C.sub(C.of(einf), C.div(C.of(wp * wp), C.mul(C.of(w), C.of(w, g)))), wenz2 = wp * wp / einf - g * g, wenz = wenz2 > 0 ? Math.sqrt(wenz2) : 0;
  const fp = wp / (2 * PI), [sc, un] = fUnit(fp), ws = N.linspace(0.03 * wp, Math.max(1.5 * wp / Math.sqrt(einf), 0.2 * wp), 300), E = ws.map(eps), wpr = inp.f_rel * (wenz || wp / Math.sqrt(einf)), ep = eps(wpr), np = csqrt(ep);
  const skin = w => { const n = csqrt(eps(w)); return C0 / (w * Math.max(n.im, 1e-12)); }, Rf = w => { const n = csqrt(eps(w)); return C.abs2(C.div(C.sub(ONE, n), C.add(ONE, n))); };
  const lamP = 2 * PI * C0 / (wenz || wp), homog = wire ? lamP / wp_.a : Infinity;
  if (wenz === 0) W.push('Damping exceeds the screened plasma frequency: Re ε never becomes negative (over-damped plasma).');
  if (wire && homog < 4) W.push(`The plasma wavelength is only ${f3(homog)} lattice periods: the wire array diffracts rather than behaving as an effective medium near its plasma frequency.`);
  if (!wire) W.push('A single Drude term is used: interband transitions (below about 520 nm for Au, 320 nm for Ag) are represented only by the constant ε∞.');
  const fl = ws.map(w => w / (2 * PI) / sc);
  return {
    summary: `${wire ? `A wire lattice (radius ${inp.wire_r} µm, period ${inp.wire_a} mm)` : `${inp.medium} (Drude model)`} has a plasma frequency of ${f3(fp / sc)} ${un} and Re ε < 0 from DC up to the epsilon-near-zero frequency ${f3(wenz / (2 * PI) / sc)} ${un}. At ${f3(wpr / (2 * PI) / sc)} ${un} the permittivity is ${cx(ep)}, the skin depth ${f3(skin(wpr) * (wire ? 1e3 : 1e9))} ${wire ? 'mm' : 'nm'} and a thick slab reflects ${f3(100 * Rf(wpr))} %.`,
    kpis: [
      kpi('fp', 'Plasma frequency', fp / sc, un, wire ? `Effective plasma frequency of the wire lattice from ${WIRE_MODEL[wp_.model][0]}; it depends only on geometry. Pendry's thin-wire estimate gives ${f3(wp_.wpPendry / (2 * PI) / sc)} ${un}.` : 'Unscreened Drude plasma frequency of the free-electron gas.'),
      kpi('f_enz', 'Epsilon-near-zero frequency', wenz / (2 * PI) / sc, un, 'Frequency where Re ε crosses zero, √(ω_p²/ε∞ − γ²); below it the medium is opaque and reflective, above it transparent.', wenz > 0),
      kpi('lam_enz', 'ENZ free-space wavelength', wenz ? 2 * PI * C0 / wenz * (wire ? 1e3 : 1e9) : 0, wire ? 'mm' : 'nm', 'Free-space wavelength at the epsilon-near-zero point; field enhancement and wavelength stretching occur around it.'),
      kpi('gamma', 'Damping rate γ/2π', g / (2 * PI) / sc, un, wire ? 'Collision rate from the finite wire conductivity, γ = ε₀a²ω_p²/(πr²σ); thinner or more resistive wires are lossier.' : 'Electron collision rate of the Drude fit; it sets ohmic loss and plasmon linewidth.'),
      kpi('eps_re', 'Re ε at probe frequency', ep.re, '–', `Real permittivity at ${f3(inp.f_rel)} × the ENZ frequency; negative values mean the field decays evanescently inside the medium.`, ep.re < 0),
      kpi('eps_im', 'Im ε at probe frequency', ep.im, '–', 'Imaginary permittivity (loss) at the probe frequency; it must be non-negative for a passive medium.', ep.im >= 0),
      kpi('skin', 'Skin depth at probe', skin(wpr), 'm', 'Field 1/e penetration depth c/(ω·Im n); a slab several skin depths thick is opaque.'),
      kpi('R', 'Half-space reflectance at probe', Rf(wpr), '–', 'Normal-incidence power reflectance |(1 − n)/(1 + n)|² of a thick slab at the probe frequency.'),
      kpi('Qp', 'Plasma quality factor ω_p/γ', wp / g, '–', 'Ratio of plasma frequency to damping; high values give a sharp ENZ transition and low loss in the negative-ε band.', wp / g > 10),
      kpi('homog', 'Plasma wavelength per lattice period', Number.isFinite(homog) ? homog : 0, '–', 'λ_p/a for the wire medium (zero for bulk metals); effective-medium behaviour needs this ratio well above about 4.', !wire || homog > 4),
    ],
    series: [
      ser('eps', 'Permittivity dispersion', `frequency (${un})`, 'ε (–)', [{ name: 'Re ε', x: fl, y: E.map(e => Math.max(e.re, -50)) }, { name: 'Im ε', x: fl, y: E.map(e => Math.min(e.im, 50)) }, { name: 'zero', x: [fl[0], fl[299]], y: [0, 0] }], 'Drude permittivity against frequency (clipped at ±50). The zero crossing of Re ε is the ENZ point; the negative-ε band lies to its left.'),
      ser('skin', 'Skin depth and reflectance', `frequency (${un})`, 'value', [{ name: wire ? 'skin depth (mm)' : 'skin depth (nm)', x: fl, y: ws.map(w => pos(Math.min(skin(w), 1e3) * (wire ? 1e3 : 1e9))) }, { name: 'reflectance (–)', x: fl, y: ws.map(w => pos(Rf(w))) }], 'Penetration depth and half-space reflectance. Below the ENZ frequency reflectance is high and the skin depth small; above it the medium becomes transparent and the skin depth diverges.', { logy: true }),
    ],
    equations: [
      eqn('constitutive', 'Drude permittivity', 'ε(ω) = ε∞ − ω_p²/(ω(ω + iγ))', 'time convention e^(−iωt), so Im ε ≥ 0 for a passive medium'),
      ...(wire ? [eqn('formula', 'Wire-lattice plasma frequency', WIRE_MODEL[wp_.model][1], wp_.model === 'lattice' ? 'cell-periodic field of a uniform surface current on each wire vanishing on the wire surface; reduces to ω_p² = 2πc₀²/(a² [ln(a/(2πr)) + 0.5275]) for thin wires' : undefined), eqn('formula', 'Damping of the wire medium', 'γ = ε₀ a² ω_p²/(π r² σ)', 'd.c. resistance of the wire over its inductance per unit length')] : matEquations([inp.medium])),
      eqn('formula', 'ENZ condition', 'Re ε(ω) = 0  ⇒  ω_ENZ = √(ω_p²/ε∞ − γ²)'), eqn('formula', 'Skin depth', 'δ = c₀/(ω Im √ε)'),
    ],
    checks: [
      chk('Passivity: Im ε ≥ 0 at all frequencies', E.every(e => e.im >= 0), `min Im ε = ${f3(N.minOf(E.map(e => e.im)))}`),
      chk('Negative-permittivity band exists', wenz > 0, wenz > 0 ? `Re ε < 0 for 0 < f < ${f3(wenz / (2 * PI) / sc)} ${un}` : 'over-damped: no band'),
      kkCheck('Causality (Kramers–Kronig) of the Drude form', w => C.sub(eps(w), C.of(einf)), wpr, wp, g / wp),
      chk('Effective-medium validity (λ_p ≫ a)', !wire || homog > 4, wire ? `λ_p/a = ${f3(homog)}` : 'bulk material'),
    ],
    warnings: W,
    recommendations: [
      rec('high', wire ? `To move the plasma frequency, change the period: f_p scales almost as 1/a (now ${f3(fp / 1e9)} GHz at a = ${inp.wire_a} mm); the logarithmic dependence on radius is weak.` : `Use ${inp.medium} between about ${f3(0.3 * wenz / (2 * PI) / sc)} and ${f3(0.9 * wenz / (2 * PI) / sc)} ${un} for plasmonic or shielding functions; near the ENZ point loss (Im ε = ${f3(eps(wenz || wp).im)}) limits field enhancement.`, 'Dispersion of the Drude response'),
      rec('medium', `Loss figure ω_p/γ = ${f3(wp / g)}: ${wire ? 'use thicker or higher-conductivity wires to raise it' : 'silver and aluminium outperform gold and ITO at their respective bands'}.`, 'Damping relative to plasma frequency'),
      rec('medium', 'Verify the retrieved ε from measured S-parameters of a slab at least three periods thick; a single layer of wires is not a medium.', 'Homogenisation limit of periodic structures'),
      rec('low', wire ? 'Printed copper or aluminium wire grids replace bulk metal with a metal fill fraction below 10⁻⁴, saving material and mass in shields and radomes.' : 'Where performance allows, substitute aluminium or transparent conducting oxides for gold and silver to reduce reliance on precious metals.', 'Material efficiency and critical raw materials', true),
    ],
    publish: { 'eff.eps': ep.re, 'res.f0': wenz / (2 * PI), 'res.bandgap_lo': 0, 'res.bandgap_hi': wenz / (2 * PI), 'eff.n': np.re },
    numerics: { method: 'Closed-form Drude dispersion; wire-lattice plasma frequency by bisection on the reciprocal-lattice sum (or the selected closed form); numerical Kramers–Kronig quadrature', dof: ws.length, notes: 'Local, isotropic (for field along the wires) effective medium; spatial dispersion of wire media not included.' },
  };
}

// ---------------------------------------------------------------- negative permeability (split-ring resonators)
// Two ring models:
//  'cylinder' — Pendry's formula for rings stacked densely along their axis (or split cylinders), ω₀² = 3 l c₀²/(π ln(2c/d) r³), F = πr²/a².
//               It assumes l ≪ r; applied to a single printed ring pair it overestimates the resonance several-fold.
//  'planar'   — printed (edge-coupled) double split ring as a quasi-static LC circuit (Marqués, Mesa, Martel & Medina): L is the inductance of one
//               flat ring of mean radius r₀ = r + c + d/2 and width c carrying a uniform sheet current; the two ring halves on either side of the
//               splits form capacitors π r₀ C′ in series, C = π r₀ C′/2, with C′ the per-unit-length capacitance of coplanar strips of width c and
//               gap d on a substrate of thickness h (conformal mapping). ω₀ = 1/√(LC). No fitted constant.
//               Smith et al. (Phys. Rev. Lett. 84, 4184, 2000): r = 1.5, c = 0.8, d = 0.2 mm, measured 4.845 GHz; this model gives 4.66 GHz on a
//               0.25 mm FR-4-type board (ε_r = 4.4, assumed — the source does not specify the laminate), the cylinder formula 49.8 GHz.
const ellipK = k => { let a = 1, b = Math.sqrt(Math.max(1 - k * k, 1e-300)); for (let i = 0; i < 40 && Math.abs(a - b) > 1e-15 * a; i++) { const t = 0.5 * (a + b); b = Math.sqrt(a * b); a = t; } return PI / (2 * a); };
const kRatio = k => (k < 1e-8 ? PI / (2 * Math.log(4 / Math.max(k, 1e-300))) : ellipK(k) / ellipK(Math.sqrt(1 - k * k)));   // K(k)/K(k′)
function cpsCapacitance(c, d, h, er) {   // F/m, coplanar strips of width c separated by d on a substrate (thickness h, permittivity er) in air
  const k = d / (d + 2 * c), x1 = PI * d / (4 * h), x2 = PI * (d + 2 * c) / (4 * h), k1 = x2 > 30 ? Math.exp(x1 - x2) * (x1 > 30 ? 1 : 1 - Math.exp(-2 * x1)) : Math.sinh(x1) / Math.sinh(x2);
  const q = 0.5 * kRatio(k1) / kRatio(k), ee = 1 + (er - 1) * N.clamp(q, 0, 0.5); return { C: EPS0 * ee / kRatio(k), ee, Cair: EPS0 / kRatio(k) };
}
// Self-inductance of a flat annulus r0 ± c/2 carrying a uniform azimuthal sheet current: L = (1/c²) ∫∫ M(ρ₁, ρ₂) dρ₁ dρ₂ with Maxwell's mutual
// inductance of coaxial coplanar filaments, M = µ₀√(ρ₁ρ₂)[(2/κ − κ)K(κ) − (2/κ)E(κ)], κ² = 4ρ₁ρ₂/(ρ₁ + ρ₂)². Midpoint rule on n × n cells; the
// logarithmic singularity of the diagonal cells is integrated analytically (cell average of ln(1/|x − y|) = 3/2 − ln h).
function ellipKE(k) { let a = 1, b = Math.sqrt(Math.max(1 - k * k, 1e-300)), cc = k, sum = 0.5 * cc * cc, pw = 0.5; for (let i = 0; i < 40 && Math.abs(cc) > 1e-16; i++) { const t = 0.5 * (a + b); cc = 0.5 * (a - b); b = Math.sqrt(a * b); a = t; pw *= 2; sum += pw * cc * cc; } const K = PI / (2 * a); return [K, K * (1 - sum)]; }
function ringInductance(r0, c) {   // H
  const n = 24, h = c / n; let L = 0;
  for (let i = 0; i < n; i++) { const p1 = r0 - c / 2 + (i + 0.5) * h; L += MU0 * p1 * (Math.log(8 * p1 / h) - 0.5);
    for (let j = 0; j < i; j++) { const p2 = r0 - c / 2 + (j + 0.5) * h, kap = 2 * Math.sqrt(p1 * p2) / (p1 + p2), [K, E] = ellipKE(kap); L += 2 * MU0 * Math.sqrt(p1 * p2) * ((2 / kap - kap) * K - (2 / kap) * E); } }
  return L / (n * n);
}
function srrParams(inp) {
  const planar = inp.srr_type !== 'cylinder', r = inp.srr_r * 1e-3, c = inp.srr_c * 1e-3, d = Math.min(inp.srr_d * 1e-3, 1.9 * c), l = inp.srr_l * 1e-3, rOut = r + 2 * c + d, r0 = r + c + d / 2;
  const aMin = planar ? 2.02 * rOut : 2.05 * r, a = Math.max(inp.srr_a * 1e-3, aMin), clampA = inp.srr_a * 1e-3 < aMin, clampD = inp.srr_d * 1e-3 > 1.9 * c;
  if (!planar) {
    const w0 = Math.sqrt(3 * l * C0 * C0 / (PI * Math.log(2 * c / d) * r ** 3)), Fm = PI * r * r / (a * a), dskin = Math.sqrt(2 / (w0 * MU0 * inp.metal_sigma)), teff = Math.min(inp.srr_t * 1e-6, dskin);
    const R1 = 1 / (inp.metal_sigma * c * teff), G = 2 * l * R1 / (r * MU0), Lr = MU0 * PI * r * r / l, Cr = 1 / (w0 * w0 * Lr);
    return { planar, r, a, c, d, l, w0, Fm, G, Lr, Cr, dskin, clampA, clampD, r0, rOut };
  }
  const cap = cpsCapacitance(c, d, Math.max(inp.srr_h * 1e-3, 1e-9), inp.srr_epsr), Lr = ringInductance(r0, c), Cr = PI * r0 * cap.C / 2, w0 = 1 / Math.sqrt(Lr * Cr);
  const dskin = Math.sqrt(2 / (w0 * MU0 * inp.metal_sigma)), teff = Math.min(inp.srr_t * 1e-6, dskin), R = (4 * PI * r0 / 3) / (inp.metal_sigma * c * teff), G = R / Lr;
  const Fraw = MU0 * PI * PI * r0 ** 4 / (Lr * a * a * l), Fm = Math.min(Fraw, 0.9);
  return { planar, r, a, c, d, l, w0, Fm, Fraw, G, Lr, Cr, dskin, clampA, clampD, r0, rOut, R, ee: cap.ee, Cpul: cap.C, f0air: 1 / Math.sqrt(Lr * PI * r0 * cap.Cair / 2) / (2 * PI) };
}
const muPendry = (w, w0, F, G) => C.sub(ONE, C.div(C.of(F * w * w), C.of(w * w - w0 * w0, G * w)));
function negMu(inp) {
  const W = [], p = srrParams(inp), { w0, Fm, G } = p, wmp = w0 / Math.sqrt(1 - Fm), f0 = w0 / (2 * PI), [sc, un] = fUnit(f0);
  const mu = w => muPendry(w, w0, Fm, G), ws = N.linspace(0.6 * w0, Math.max(1.6 * w0, 1.25 * wmp), 400), Mu = ws.map(mu), neg = band(ws, Mu.map(m => m.re < 0)), wpr = inp.f_rel * w0, mp = mu(wpr), imin = N.argmin(Mu.map(m => m.re));
  const lamA = 2 * PI * C0 / w0 / p.a;
  if (p.clampA) W.push(p.planar ? `Lattice period raised to ${f3(p.a * 1e3)} mm: the outer ring (radius r + 2c + d = ${f3(p.rOut * 1e3)} mm) must fit inside the cell.` : 'Lattice period raised to 2.05 ring radii so that neighbouring rings do not overlap.');
  if (p.planar && p.Fraw > 0.9) W.push('The magnetic oscillator strength of this dense stack exceeds 0.9 and was limited: rings this closely stacked couple strongly and the dilute mixing rule no longer holds (use the cylinder formula for l ≪ r).');
  if (p.planar && p.c > 0.6 * p.r0) W.push('The ring width exceeds 60 % of the mean radius: the uniform-current ring inductance and the coplanar-strip capacitance (which ignores curvature) become rough estimates.');
  if (p.planar && p.d < inp.srr_t * 1e-6) W.push('The gap between the rings is smaller than the metal thickness: the parallel-plate capacitance of the facing ring edges, not included in the zero-thickness coplanar-strip formula, lowers the real resonance.');
  if (p.planar && p.l < 0.5 * p.r0) W.push('Rings stacked closer than about half their radius couple inductively: the single-ring inductance used here is too large and the resonance of the stack lies higher (towards the cylinder formula).');
  if (!p.planar && p.l > p.r) W.push(`The cylinder formula assumes rings stacked much closer than their radius (l ≪ r); here l/r = ${f3(p.l / p.r)}. For printed ring pairs that are not densely stacked it overestimates the resonance several-fold (49.8 GHz instead of the measured 4.845 GHz for the rings of Smith et al., 2000): select the planar model.`);
  if (p.clampD) W.push('Ring gap limited to 1.9 ring widths: the capacitance formula needs ln(2c/d) > 0.');
  if (!neg) W.push('Loss is too high for Re µ to become negative: use thicker or more conductive rings, or a larger filling factor.');
  if (lamA < 5) W.push(`The free-space wavelength at resonance is only ${f3(lamA)} lattice periods; the effective-medium description is marginal.`);
  const fl = ws.map(w => w / (2 * PI) / sc);
  return {
    summary: `The split-ring medium resonates at ${f3(f0 / sc)} ${un} (${p.planar ? `planar edge-coupled rings, quasi-static LC model: L = ${f3(p.Lr * 1e9)} nH, C = ${f3(p.Cr * 1e12)} pF, substrate ε_eff = ${f3(p.ee)}; ` : 'Pendry cylinder/stack formula; '}inner ring radius ${inp.srr_r} mm, oscillator strength ${f3(Fm)}). ${neg ? `Re µ is negative from ${f3(neg[0] / (2 * PI) / sc)} to ${f3(neg[1] / (2 * PI) / sc)} ${un}, reaching ${f3(Mu[imin].re)}` : 'Re µ stays positive'}; the lossless magnetic plasma frequency is ${f3(wmp / (2 * PI) / sc)} ${un}. At the probe frequency µ = ${cx(mp)}.`,
    kpis: [
      kpi('f0', 'Magnetic resonance frequency', f0 / sc, un, p.planar ? `Quasi-static LC resonance of the printed ring pair, ω₀ = 1/√(LC) with the flat-ring inductance and the series capacitance of the two ring halves across the inter-ring gap; without the substrate it would be ${f3(p.f0air / sc)} ${un}.` : 'LC resonance of the split rings from Pendry\'s model for densely stacked rings, ω₀² = 3l·c₀²/(π ln(2c/d) r³); it scales as r^(−3/2).'),
      kpi('fmp', 'Magnetic plasma frequency', wmp / (2 * PI) / sc, un, 'Upper edge of the lossless negative-µ band, ω₀/√(1 − F); the band widens with filling factor.'),
      kpi('F', 'Filling factor', Fm, '–', p.planar ? 'Oscillator strength of the magnetic resonance, µ₀π²r₀⁴/(L a² l): the magnetic polarisability of one ring pair times the ring density; it plays the role of the filling factor πr²/a² of the cylinder model.' : 'Fraction of the unit-cell area enclosed by the ring, πr²/a²; it is the strength of the magnetic resonance.'),
      kpi('neg_lo', 'Negative-µ band lower edge', neg ? neg[0] / (2 * PI) / sc : 0, un, 'Lowest frequency in the sweep with Re µ < 0 including loss (zero if there is no negative band).', !!neg),
      kpi('neg_hi', 'Negative-µ band upper edge', neg ? neg[1] / (2 * PI) / sc : 0, un, 'Highest frequency with Re µ < 0 including loss; slightly below the lossless magnetic plasma frequency.', !!neg),
      kpi('mu_min', 'Most negative Re µ', Mu[imin].re, '–', `Minimum of the real permeability, reached at ${f3(fl[imin])} ${un} just above resonance; loss limits how negative it can go.`, Mu[imin].re < 0),
      kpi('Q', 'Resonator quality factor', w0 / G, '–', p.planar ? `ω₀L/R with R = ${f3(p.R)} Ω, the ohmic resistance of the two rings for a current that swaps linearly between them (⅔ of a single closed ring), skin-depth limited; dielectric and radiation loss are not included, so measured values are lower.` : 'ω₀/Γ with Γ = 2lσ₁/(rµ₀) from the ohmic resistance of the ring per unit circumference (skin-depth limited).', w0 / G > 30),
      kpi('L', 'Ring inductance (per ring)', p.Lr * 1e9, 'nH', p.planar ? `Self-inductance of one flat ring of mean radius ${f3(p.r0 * 1e3)} mm and width ${f3(p.c * 1e3)} mm with uniform sheet current (double integral of the mutual inductance of coaxial filaments); close to µ₀r₀[ln(8r₀/c) − ½] for narrow rings.` : 'Solenoid-like inductance µ₀πr²/l of a ring in the stack; with the capacitance it reproduces the resonance.'),
      kpi('Cring', 'Effective ring capacitance', p.Cr * 1e12, 'pF', p.planar ? `Series capacitance πr₀C′/2 of the two ring halves, with C′ = ${f3(p.Cpul * 1e12)} pF/m for coplanar strips on the substrate (effective permittivity ${f3(p.ee)}).` : 'Capacitance 1/(ω₀²L) between the two concentric rings implied by the resonance frequency.'),
      kpi('mu_re', 'Re µ at probe frequency', mp.re, '–', `Real effective permeability at ${f3(inp.f_rel)} × f₀; choose a probe between the resonance and the magnetic plasma frequency to obtain a negative value.`, mp.re < 0),
      kpi('mu_im', 'Im µ at probe frequency', mp.im, '–', 'Magnetic loss at the probe frequency; non-negative for a passive medium.', mp.im >= 0),
      kpi('lam_a', 'Wavelength per lattice period', lamA, '–', 'Free-space wavelength at resonance divided by the lattice constant; values above about 5–10 justify homogenisation.', lamA > 5),
    ],
    series: [
      ser('mu', 'Effective permeability', `frequency (${un})`, 'µ (–)', [{ name: 'Re µ', x: fl, y: Mu.map(m => m.re) }, { name: 'Im µ', x: fl, y: Mu.map(m => m.im) }, { name: 'zero', x: [fl[0], fl[399]], y: [0, 0] }], 'Lorentz-type magnetic response. Re µ swings negative just above the resonance and returns through zero at the magnetic plasma frequency; Im µ peaks at resonance.'),
      ser('tand', 'Magnetic loss tangent', `frequency (${un})`, 'Im µ/|Re µ| (–)', [{ name: 'loss tangent', x: fl, y: Mu.map(m => pos(Math.min(m.im / Math.max(Math.abs(m.re), 1e-6), 1e3))) }], 'Ratio of magnetic loss to stored response. Operate towards the upper part of the negative band where the loss tangent is lower, at the price of a smaller |µ|.', { logy: true }),
    ],
    equations: [
      eqn('constitutive', 'Split-ring effective permeability (Pendry)', 'µ(ω) = 1 − F ω²/(ω² − ω₀² + iΓω)', 'e^(−iωt) convention'), ...(p.planar ? [eqn('formula', 'Planar edge-coupled ring: resonance', 'ω₀ = 1/√(L C),   C = π r₀ C′/2,   r₀ = r + c + d/2', 'two ring halves of capacitance π r₀ C′ in series across the splits'),
        eqn('formula', 'Inter-ring capacitance per unit length (coplanar strips)', 'C′ = ε₀ ε_eff K(k′)/K(k),  k = d/(d + 2c),  ε_eff = 1 + ½(ε_r − 1)·[K(k′)/K(k)]·[K(k₁)/K(k₁′)],  k₁ = sinh(πd/4h)/sinh(π(d + 2c)/4h)'),
        eqn('formula', 'Ring inductance (uniform sheet current)', 'L = (1/c²) ∫∫ M(ρ₁, ρ₂) dρ₁ dρ₂,   M = µ₀ √(ρ₁ρ₂) [(2/κ − κ) K(κ) − (2/κ) E(κ)],   κ² = 4ρ₁ρ₂/(ρ₁ + ρ₂)²', 'integration across the ring width r₀ ± c/2; tends to µ₀ r₀ [ln(8r₀/c) − ½] for narrow rings'),
        eqn('formula', 'Oscillator strength and damping', 'F = µ₀ π² r₀⁴/(L a² l),   Γ = R/L,   R = (4π r₀/3)/(σ c t_eff)')]
        : [eqn('formula', 'Resonance and damping', 'ω₀² = 3 l c₀²/(π ln(2c/d) r³),  Γ = 2 l σ₁/(r µ₀),  F = π r²/a²'), eqn('formula', 'Equivalent circuit', 'L = µ₀ π r²/l,  C = 1/(ω₀² L)')]),
      eqn('formula', 'Negative band (lossless)', 'ω₀ < ω < ω₀/√(1 − F)'),
    ],
    checks: [
      chk('Passivity: Im µ ≥ 0', Mu.every(m => m.im >= 0), `min Im µ = ${f3(N.minOf(Mu.map(m => m.im)))}`),
      kkCheck('Causality (Kramers–Kronig) of µ − µ(∞)', w => C.sub(mu(w), C.of(1 - Fm)), wpr, w0, G / w0),
      ...(p.planar ? [chk('Ring inductance consistent with the narrow-ring asymptote', Math.abs(p.Lr / (MU0 * p.r0 * (Math.log(8 * p.r0 / p.c) - 0.5)) - 1) < 0.15, `numerical ${f3(p.Lr * 1e9)} nH vs µ₀r₀[ln(8r₀/c) − ½] = ${f3(MU0 * p.r0 * (Math.log(8 * p.r0 / p.c) - 0.5) * 1e9)} nH`),
        chk('Quasi-static ring (circumference below half the guided wavelength)', 2 * PI * p.r0 < 2 * PI * C0 / w0 / (2 * Math.sqrt(p.ee)), `2πr₀ = ${f3(2 * PI * p.r0 * 1e3)} mm vs half the guided wavelength ${f3(2 * PI * C0 / w0 / (2 * Math.sqrt(p.ee)) * 1e3)} mm`)] : []),
      chk('Pendry reference geometry reproduced', Math.abs(Math.sqrt(3 * 2e-3 * C0 * C0 / (PI * Math.log(20) * 8e-9)) / (2 * PI) / 1e9 - 13.5) < 0.1, 'r = 2 mm, c = 1 mm, d = 0.1 mm, l = 2 mm gives 13.5 GHz as in the original model'),
      chk('Effective-medium validity (λ₀/a > 5)', lamA > 5, `λ₀/a = ${f3(lamA)}`), chk('Negative-permeability band present', !!neg, neg ? `width ${f3((neg[1] - neg[0]) / (2 * PI) / sc)} ${un}` : 'none with this loss level'),
    ],
    warnings: W,
    recommendations: [
      rec('high', neg ? `Operate between ${f3(neg[0] / (2 * PI) / sc)} and ${f3(neg[1] / (2 * PI) / sc)} ${un}; the fractional bandwidth is ${f3(100 * (neg[1] - neg[0]) / w0)} %, so fabrication tolerance on r (f₀ ∝ r^(−3/2)) must be better than about ${f3(100 * (neg[1] - neg[0]) / w0 / 3)} %.` : 'Increase the filling factor (larger rings in the same cell) and the metal thickness to obtain a negative band.', 'Bandwidth of the Lorentz resonance'),
      rec('medium', `Q = ${f3(w0 / G)}: skin depth at resonance is ${f3(p.dskin * 1e6)} µm, so metal thicker than about ${f3(3 * p.dskin * 1e6)} µm brings no further benefit.`, 'Skin-effect-limited ring resistance'),
      rec('medium', 'The response is magnetic only for H along the ring axis; use three orthogonal ring sets or accept a uniaxial (anisotropic) µ tensor.', 'Anisotropy of planar split rings'),
      rec('low', 'Standard printed-circuit laminates and thin copper suffice; avoid thick precious-metal plating, which the skin depth makes redundant.', 'Material efficiency from the skin-depth result', true),
    ],
    publish: { 'eff.mu': mp.re, 'res.f0': f0, 'res.Q': w0 / G, ...(neg ? { 'res.bandgap_lo': neg[0] / (2 * PI), 'res.bandgap_hi': neg[1] / (2 * PI) } : {}) },
    numerics: { method: p.planar ? 'Quasi-static LC model of the planar ring pair (24 × 24 quadrature of the filament mutual inductance for L, arithmetic–geometric mean for the elliptic integrals); Lorentz permeability; numerical Kramers–Kronig quadrature' : 'Pendry split-ring homogenisation (closed form); numerical Kramers–Kronig quadrature', dof: ws.length, notes: 'Quasi-static rings, no bianisotropy or inter-ring coupling beyond the mean-field depolarisation.' },
  };
}

// ---------------------------------------------------------------- negative refractive index composite
function slabS(n, z, k0d) { const R01 = C.div(C.sub(z, ONE), C.add(z, ONE)), ph = C.exp(C.mul(J, C.scale(n, k0d))), ph2 = C.mul(ph, ph), den = C.sub(ONE, C.mul(C.mul(R01, R01), ph2)); return { S11: C.div(C.mul(R01, C.sub(ONE, ph2)), den), S21: C.div(C.mul(C.sub(ONE, C.mul(R01, R01)), ph), den) }; }
function nim(inp) {
  const W = [], G9 = 2 * PI * 1e9, wep = inp.fep * G9, ge = inp.ge * G9, w0 = inp.f0m * G9, Fm = inp.Fm, gm = inp.gm * G9, d = inp.d_slab * 1e-3, th = inp.theta * PI / 180;
  const eps = w => C.sub(ONE, C.div(C.of(wep * wep), C.mul(C.of(w), C.of(w, ge)))), mu = w => muPendry(w, w0, Fm, gm), nz = w => { const se = csqrt(eps(w)), sm = csqrt(mu(w)); return { n: C.mul(se, sm), z: C.div(sm, se), e: eps(w), m: mu(w) }; };
  const wmp = w0 / Math.sqrt(1 - Fm), wlo = 0.5 * Math.min(w0, wep), whi = 1.5 * Math.max(wmp, Math.min(wep, 3 * wmp)), ws = N.linspace(wlo, whi, 400), P = ws.map(nz), fl = ws.map(w => w / G9);
  // S-parameters of the slab and retrieval
  const Sp = ws.map((w, i) => slabS(P[i].n, P[i].z, w * d / C0)), ret = []; let prev = null, okRet = true;
  const invert = (S11, S21) => { // → { z, arg, mag } with Re z ≥ 0
    const a = C.sub(C.mul(C.add(ONE, S11), C.add(ONE, S11)), C.mul(S21, S21)), b = C.sub(C.mul(C.sub(ONE, S11), C.sub(ONE, S11)), C.mul(S21, S21));
    let z = C.sqrt(C.div(a, b)), X = C.div(S21, C.sub(ONE, C.mul(S11, C.div(C.sub(z, ONE), C.add(z, ONE)))));
    if (z.re < -1e-9 || (Math.abs(z.re) <= 1e-9 && C.abs(X) > 1)) { z = C.scale(z, -1); X = C.div(S21, C.sub(ONE, C.mul(S11, C.div(C.sub(z, ONE), C.add(z, ONE))))); }
    return { z, arg: Math.atan2(X.im, X.re), mag: C.abs(X) };
  };
  const rawAt = w => { const q = nz(w), s = slabS(q.n, q.z, w * d / C0); return invert(s.S11, s.S21); }, near = (arg, ref) => arg + 2 * PI * Math.round((ref - arg) / (2 * PI));
  // Continuous complex log of e^(ink₀d). An interval is accepted only if (i) its midpoint and end are each a small step in phase,
  // attenuation and impedance (z has no branch ambiguity, so it exposes a resonance hidden between two samples) and (ii) the local
  // phase slope, measured over a micro-step that cannot alias, predicts less than a quarter turn across it; otherwise it is bisected.
  const zdiff = (a, b) => Math.min(C.abs(C.sub(a, b)), C.abs(C.add(a, b))) / Math.max(C.abs(a), C.abs(b), 1e-300);
  const stepTo = (qa, r) => ({ p: near(r.arg, qa.p), lm: Math.log(pos(r.mag)), z: r.z }), smooth = (a, b) => Math.abs(b.p - a.p) < 0.25 * PI && Math.abs(b.lm - a.lm) < 1 && zdiff(a.z, b.z) < 0.3;
  const probe = (w, h) => { const r = rawAt(w), r1 = rawAt(w + h / 1024), dp = r1.arg - r.arg; return { r, seen: r.mag > 1e-280 && r1.mag > 1e-280, adv: 1024 * Math.abs(dp - 2 * PI * Math.round(dp / (2 * PI))) }; };
  let evals = 0, lost = false, wPrev = 0; const track = (wa, qa, wb, depth) => {
    const h = wb - wa, wm = wa + 0.5 * h, sa = probe(wa, h), sm = probe(wm, h), sb = probe(wb, h); evals++;
    if (!(sa.seen && sm.seen && sb.seen)) { lost = true; return stepTo(qa, sb.r); } // opaque: the phase of S21 is no longer resolved
    const qm = stepTo(qa, sm.r), qb = stepTo(qm, sb.r), stop = depth > 40 || evals > 30000;
    if ((smooth(qa, qm) && smooth(qm, qb) && Math.max(sa.adv, sm.adv, sb.adv) < 0.5 * PI) || stop) { if (stop) lost = true; return qb; }
    return track(wm, track(wa, qa, wm, depth + 1), wb, depth + 1);
  };
  // Branch seed: lower the frequency until the slab is electrically thin (|ln e^(ink₀d)| ≪ 1 at three successive octaves), where the
  // principal branch is unambiguous, then follow the phase up to the start of the sweep.
  for (let k = 0, w = ws[0]; k < 80 && prev === null; k++, w /= 4) { const r = [rawAt(w), rawAt(w / 2), rawAt(w / 4)]; if (r.every(q => q.mag > 1e-280 && Math.abs(q.arg) < 0.05 && Math.abs(Math.log(q.mag)) < 0.05)) { prev = { p: r[0].arg, lm: Math.log(r[0].mag), z: r[0].z }; wPrev = w; } }
  ws.forEach((w, i) => {
    const k0d = w * d / C0, q = invert(Sp[i].S11, Sp[i].S21), seen = q.mag > 1e-280;
    if (!lost && seen && prev !== null) { prev = track(wPrev, prev, w, 0); wPrev = w; } else lost = true;
    if (lost || !seen) { okRet = false; ret.push({ n: P[i].n, z: q.z }); return; } // not retrievable: the model value is plotted and excluded from the error
    ret.push({ n: { re: prev.p / k0d, im: -Math.log(q.mag) / k0d }, z: q.z });
  });
  const errN = N.maxOf(ret.map((r, i) => C.abs(C.sub(r.n, P[i].n)))), errZ = N.maxOf(ret.map((r, i) => C.abs(C.sub(r.z, P[i].z)) / (1 + C.abs(P[i].z))));
  const dng = band(ws, P.map(p => p.e.re < 0 && p.m.re < 0)), nneg = band(ws, P.map(p => p.n.re < 0)), fom = P.map(p => (p.n.re < 0 ? -p.n.re / Math.max(p.n.im, 1e-12) : 0)), ifo = N.argmax(fom);
  const wpr = inp.f_probe * G9, pp = nz(wpr), sp = slabS(pp.n, pp.z, wpr * d / C0), sn = Math.sin(th) / pp.n.re, refr = Math.abs(sn) <= 1 && Math.abs(pp.n.re) > 1e-9, tht = refr ? Math.asin(sn) * 180 / PI : 0;
  const dw = wpr * 1e-4, ng = pp.n.re + wpr * (nz(wpr + dw).n.re - nz(wpr - dw).n.re) / (2 * dw), rp = ret[N.argmin(ws.map(w => Math.abs(w - wpr)))];
  if (!dng) W.push('No double-negative band: the electric plasma frequency must lie above the magnetic resonance so that ε < 0 where µ < 0.');
  if (!okRet) W.push('The slab is so opaque (or its phase turns so quickly) at some frequencies that the phase of S21 cannot be followed; the retrieval was skipped there and at all higher frequencies, because the branch of Re n cannot be continued across such a band.');
  if (wpr < wlo || wpr > whi) W.push('The probe frequency lies outside the plotted sweep.');
  if (pp.n.re < 0 && pp.n.im > 0.5 * Math.abs(pp.n.re)) W.push('At the probe frequency the figure of merit is below 2: the wave is absorbed within about one wavelength.');
  return {
    summary: `${dng ? `The wire–split-ring composite is double-negative from ${f3(dng[0] / G9)} to ${f3(dng[1] / G9)} GHz` : 'The composite has no double-negative band'}${nneg ? ` and Re n < 0 from ${f3(nneg[0] / G9)} to ${f3(nneg[1] / G9)} GHz` : ''}. At ${f3(inp.f_probe)} GHz: ε = ${cx(pp.e)}, µ = ${cx(pp.m)}, n = ${cx(pp.n)}, z = ${cx(pp.z)} (figure of merit ${f3(Math.abs(pp.n.re) / Math.max(pp.n.im, 1e-12))}). ${refr ? `A beam incident at ${inp.theta}° refracts to ${f3(tht)}°.` : 'Snell refraction has no real angle at the probe.'}`,
    kpis: [
      kpi('n_re', 'Re n at probe frequency', pp.n.re, '–', 'Real refractive index on the causal, passive branch (Im n ≥ 0); it is negative when ε and µ are simultaneously negative, giving antiparallel phase and energy velocities.', pp.n.re < 0),
      kpi('n_im', 'Im n at probe frequency', pp.n.im, '–', 'Extinction index; the amplitude decays by e over λ₀/(2π Im n). It is non-negative for a passive medium.', pp.n.im >= 0),
      kpi('fom', 'Figure of merit |Re n|/Im n', Math.abs(pp.n.re) / Math.max(pp.n.im, 1e-12), '–', 'Phase advance per unit attenuation at the probe; values above about 3–10 are needed for lens-like applications.', Math.abs(pp.n.re) / Math.max(pp.n.im, 1e-12) > 3),
      kpi('z_re', 'Re z at probe frequency', pp.z.re, '–', 'Real part of the wave impedance relative to free space, √(µ/ε) with Re z ≥ 0; z = 1 means no reflection at the interface.'),
      kpi('eps_re', 'Re ε at probe frequency', pp.e.re, '–', 'Real permittivity of the wire sub-lattice (Drude form) at the probe frequency.'),
      kpi('mu_re', 'Re µ at probe frequency', pp.m.re, '–', 'Real permeability of the split-ring sub-lattice (Lorentz form) at the probe frequency.'),
      kpi('dng_lo', 'Double-negative band lower edge', dng ? dng[0] / G9 : 0, 'GHz', 'Lowest frequency at which Re ε and Re µ are both negative (zero if the band is absent).', !!dng),
      kpi('dng_hi', 'Double-negative band upper edge', dng ? dng[1] / G9 : 0, 'GHz', 'Highest frequency at which Re ε and Re µ are both negative.', !!dng),
      kpi('fom_max', 'Best figure of merit in the negative band', fom[ifo], '–', `Maximum of −Re n/Im n over the sweep, found at ${f3(fl[ifo])} GHz; usually near the upper end of the band where |n| is small.`),
      kpi('theta_t', 'Refraction angle (Snell)', tht, '°', `Angle of the refracted wave vector for ${inp.theta}° incidence from air; negative values mean refraction to the same side of the normal as the incident beam.`, refr ? tht < 0 : null),
      kpi('S21', '|S21|² of the slab at probe', C.abs2(sp.S21), '–', `Power transmission through a ${inp.d_slab} mm slab in air at normal incidence, from the Fresnel–Airy slab formula.`),
      kpi('S11', '|S11|² of the slab at probe', C.abs2(sp.S11), '–', 'Power reflection of the slab; with S21 it is the input to the parameter retrieval.'),
      kpi('ng', 'Group index at probe', ng, '–', 'n + ω dn/dω from a numerical derivative; it stays positive in a passive medium even where the phase index is negative (backward wave).', ng > 0),
      kpi('ret_err', 'Retrieval error in n', errN, '–', 'Largest difference between the index retrieved from S11/S21 with branch tracking and the index of the constitutive model across the sweep.', errN < 1e-6),
    ],
    series: [
      ser('n', 'Refractive index', 'frequency (GHz)', 'n (–)', [{ name: 'Re n (model)', x: fl, y: P.map(p => p.n.re) }, { name: 'Im n (model)', x: fl, y: P.map(p => p.n.im) }, { name: 'Re n (retrieved from S-parameters)', x: fl, y: ret.map(r => r.n.re) }], 'Complex refractive index on the passive branch, and the value recovered from the slab S-parameters. The band with Re n < 0 is the negative-index region; the two curves must coincide if the branch is tracked correctly.'),
      ser('em', 'Constitutive parameters', 'frequency (GHz)', 'value (–)', [{ name: 'Re ε', x: fl, y: P.map(p => Math.max(p.e.re, -20)) }, { name: 'Re µ', x: fl, y: P.map(p => N.clamp(p.m.re, -20, 20)) }, { name: 'Re z', x: fl, y: P.map(p => Math.min(p.z.re, 20)) }], 'Real permittivity, permeability and impedance (clipped at ±20). The negative-index band is where both ε and µ are below zero; impedance matching (z ≈ 1) minimises reflection.'),
      ser('S', 'Slab S-parameters', 'frequency (GHz)', 'magnitude (dB)', [{ name: '|S21|', x: fl, y: Sp.map(s => Math.max(20 * Math.log10(pos(C.abs(s.S21))), -120)) }, { name: '|S11|', x: fl, y: Sp.map(s => Math.max(20 * Math.log10(pos(C.abs(s.S11))), -120)) }], 'Transmission and reflection of the finite slab. A pass band appearing inside what would be the stop band of either sub-lattice alone is the signature of a double-negative medium.'),
      ser('fom', 'Figure of merit in the negative-index band', 'frequency (GHz)', '−Re n / Im n (–)', [{ name: 'figure of merit', x: fl, y: fom }], 'Ratio of phase index to extinction where Re n < 0 (zero elsewhere). It rises towards the upper band edge as the resonance loss falls.'),
    ],
    equations: [
      eqn('constitutive', 'Wire-medium permittivity', 'ε(ω) = 1 − ω_ep²/(ω(ω + iγ_e))'), eqn('constitutive', 'Split-ring permeability', 'µ(ω) = 1 − F ω²/(ω² − ω₀² + iγ_m ω)'),
      eqn('governing', 'Index and impedance (passive branch)', 'n = √ε·√µ with Im n ≥ 0;  z = √µ/√ε with Re z ≥ 0'), eqn('formula', 'Slab S-parameters', 'S11 = R₀₁(1 − e^(2ink₀d))/(1 − R₀₁² e^(2ink₀d)),  S21 = (1 − R₀₁²) e^(ink₀d)/(1 − R₀₁² e^(2ink₀d)),  R₀₁ = (z − 1)/(z + 1)'),
      eqn('formula', 'S-parameter retrieval', 'z = ±√{[(1 + S11)² − S21²]/[(1 − S11)² − S21²]},  e^(ink₀d) = S21/[1 − S11(z − 1)/(z + 1)],  n = [−i ln(e^(ink₀d)) + 2πm]/(k₀d)', 'sign from Re z ≥ 0; branch m by phase continuity from the low-frequency end, with adaptive frequency refinement near the resonance'),
      eqn('formula', 'Snell law', 'sin θ_t = sin θ_i / Re n'),
    ],
    checks: [
      chk('Passivity: Im ε ≥ 0 and Im µ ≥ 0', P.every(p => p.e.im >= 0 && p.m.im >= 0), 'both constitutive parameters dissipate at every frequency'), chk('Passive branch: Im n ≥ 0 and Re z ≥ 0', P.every(p => p.n.im >= -1e-12 && p.z.re >= -1e-12), 'branch chosen so that waves decay in the propagation direction'),
      chk('Slab is passive: |S11|² + |S21|² ≤ 1', Sp.every(s => C.abs2(s.S11) + C.abs2(s.S21) <= 1 + 1e-9), `max ${f4(N.maxOf(Sp.map(s => C.abs2(s.S11) + C.abs2(s.S21))))}`),
      chk('Retrieval recovers the model (n, z)', errN < 1e-6 && errZ < 1e-6, `max |Δn| = ${f3(errN)}, relative |Δz| = ${f3(errZ)}; retrieved at probe: ε = ${f3(C.div(rp.n, rp.z).re)}, µ = ${f3(C.mul(rp.n, rp.z).re)}${okRet ? '' : '; frequencies where the phase of S21 could not be followed are excluded'}`),
      kkCheck('Causality (Kramers–Kronig) of µ − µ(∞)', w => C.sub(mu(w), C.of(1 - Fm)), wpr, w0, gm / w0), kkCheck('Causality (Kramers–Kronig) of ε − 1', w => C.sub(eps(w), ONE), wpr, wep, ge / wep),
      chk('Group velocity positive (energy flows forward)', ng > 0, `n_g = ${f3(ng)} while the phase index is ${f3(pp.n.re)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', dng ? `Work near ${f3(fl[ifo])} GHz where the figure of merit peaks at ${f3(fom[ifo])}; at the present probe it is ${f3(Math.abs(pp.n.re) / Math.max(pp.n.im, 1e-12))}.` : `Raise the electric plasma frequency above ${f3(wmp / G9)} GHz (denser wires) or lower the ring resonance so the two negative bands overlap.`, 'Overlap of the ε < 0 and µ < 0 bands and loss dispersion'),
      rec('medium', `For low reflection choose the frequency where ε ≈ µ (z ≈ 1): here z = ${f3(pp.z.re)}, giving |S11|² = ${f3(C.abs2(sp.S11))} for the ${inp.d_slab} mm slab.`, 'Impedance matching condition'),
      rec('medium', 'When retrieving from measured S-parameters use two slab thicknesses to fix the branch index m, and keep the slab below half a guided wavelength at the lowest frequency.', 'Branch ambiguity of the logarithm in the retrieval'),
      rec('low', 'Reduce γ_m first (thicker copper, low-loss substrate): magnetic resonance loss dominates the figure of merit.', 'Im µ ≫ Im ε inside the band'),
      rec('low', 'Design for the narrowest slab that gives the required phase: negative-index operation is resonant and dissipative, so every extra cell converts signal power to heat.', 'Energy efficiency of resonant metamaterials', true),
    ],
    publish: { 'eff.n': pp.n.re, 'eff.eps': pp.e.re, 'eff.mu': pp.m.re, 'res.f0': w0 / (2 * PI), 'res.transmission': C.abs2(sp.S21), ...(dng ? { 'res.bandgap_lo': dng[0] / (2 * PI), 'res.bandgap_hi': dng[1] / (2 * PI) } : {}) },
    numerics: { method: 'Closed-form Drude–Lorentz constitutive model; Fresnel–Airy slab; S-parameter inversion with adaptive phase tracking', dof: ws.length, iterations: evals, notes: 'Homogeneous isotropic slab at normal incidence for the S-parameters; unit-cell granularity and bianisotropy not modelled.' },
  };
}

// ---------------------------------------------------------------- Mie theory and quasi-static spheroids
function mieQ(x, m) { // Bohren–Huffman series; m complex relative index
  const mx = C.scale(m, x), nstop = Math.floor(x + 4 * Math.cbrt(x) + 2), nmx = Math.floor(Math.max(nstop, C.abs(mx)) + 16), D = new Array(nmx + 1); D[nmx] = { re: 0, im: 0 };
  for (let n = nmx; n > 0; n--) { const a = C.div(C.of(n), mx); D[n - 1] = C.sub(a, C.inv(C.add(D[n], a))); }
  let p0 = Math.cos(x), p1 = Math.sin(x), c0 = -Math.sin(x), c1 = Math.cos(x), qs = 0, qe = 0;
  for (let n = 1; n <= nstop; n++) {
    const pn = (2 * n - 1) / x * p1 - p0, cn = (2 * n - 1) / x * c1 - c0, xi = { re: pn, im: -cn }, xi1 = { re: p1, im: -c1 };
    const ta = C.add(C.div(D[n], m), C.of(n / x)), tb = C.add(C.mul(m, D[n]), C.of(n / x));
    const an = C.div(C.sub(C.scale(ta, pn), C.of(p1)), C.sub(C.mul(ta, xi), xi1)), bn = C.div(C.sub(C.scale(tb, pn), C.of(p1)), C.sub(C.mul(tb, xi), xi1));
    qs += (2 * n + 1) * (C.abs2(an) + C.abs2(bn)); qe += (2 * n + 1) * (an.re + bn.re); p0 = p1; p1 = pn; c0 = c1; c1 = cn;
  }
  return { ext: 2 * qe / (x * x), sca: 2 * qs / (x * x), abs: 2 * (qe - qs) / (x * x), nstop };
}
function mie(inp) {
  const W = [], mat = inp.mat_p, a = inp.d_p * 0.5e-9, nm = inp.n_med, em = nm * nm, AR = inp.AR_p, lmin = Math.min(inp.lam_min, inp.lam_max - 10), lams = N.linspace(lmin, inp.lam_max, 240), area = PI * a * a;
  dispWarn([mat], lmin, W);
  const run = (nmed, l) => mieQ(2 * PI * nmed * a / (l * 1e-9), C.scale(nOf(mat, l, inp), 1 / nmed)), Q = lams.map(l => run(nm, l)), ip = N.argmax(Q.map(q => q.ext));
  const peakOf = nmed => { const q = lams.map(l => run(nmed, l).ext), i = N.argmax(q); return i > 0 && i < lams.length - 1 ? N.goldenSection(l => -run(nmed, l).ext, lams[i - 1], lams[i + 1], 1e-4) : lams[i]; };
  const lp = peakOf(nm), qp = run(nm, lp), sens = (peakOf(nm + 0.05) - lp) / 0.05;
  // quasi-static spheroid (Gans): prolate, minor axis = d, major = AR·d
  const e2 = 1 - 1 / (AR * AR), ee = Math.sqrt(e2), Lz = AR > 1.0001 ? (1 - e2) / e2 * (Math.log((1 + ee) / (1 - ee)) / (2 * ee) - 1) : 1 / 3, Lx = (1 - Lz) / 2, V = 4 / 3 * PI * a ** 3 * AR;
  const alpha = (l, Lj) => { const e = epsOf(mat, l, inp), de = C.sub(e, C.of(em)); return C.scale(C.div(de, C.add(C.of(em), C.scale(de, Lj))), V); };
  const cs = (l, Lj) => { const k = 2 * PI * nm / (l * 1e-9), al = alpha(l, Lj); return k * al.im + k ** 4 * C.abs2(al) / (6 * PI); };
  const qsS = lams.map(l => cs(l, 1 / 3) / AR / area), gL = lams.map(l => cs(l, Lz) / area), gT = lams.map(l => cs(l, Lx) / area), iq = N.argmax(qsS), iL = N.argmax(gL);
  const fe = (l, Lj) => { const e = epsOf(mat, l, inp); return C.abs(C.div(e, C.add(C.of(em), C.scale(C.sub(e, C.of(em)), Lj)))); }, feS = fe(lams[iq], 1 / 3), feL = fe(lams[iL], Lz);
  let lF = 0; for (let i = 1; i < lams.length; i++) { const f0 = epsOf(mat, lams[i - 1], inp).re + 2 * em, f1 = epsOf(mat, lams[i], inp).re + 2 * em; if (f0 * f1 <= 0 && f0 !== f1) { lF = lams[i - 1] + (lams[i] - lams[i - 1]) * f0 / (f0 - f1); break; } }
  const x0 = 2 * PI * nm * a / (lp * 1e-9), plasm = !!MATS[mat].metal, edge = ip === 0 || ip === lams.length - 1;
  if (edge) W.push('The extinction maximum lies at the edge of the wavelength window; widen the range to capture the resonance.');
  if (plasm && inp.d_p < 10) W.push('Below about 10 nm the bulk dielectric function underestimates damping (surface scattering of electrons); real resonances are broader and weaker.');
  if (plasm && mat === 'Au' ) W.push('The Lorentz–Drude fit for gold is smooth through the interband region; expect the measured resonance within about ±10 nm of the prediction.');
  if (!plasm) W.push('Dielectric particles show Mie (geometric) resonances rather than a localised surface plasmon; the Fröhlich and Gans results are not meaningful for them.');
  if (x0 > 0.6) W.push(`Size parameter ${f3(x0)} at the peak: retardation is significant and the quasi-static (Fröhlich/Gans) curves are only qualitative; use the full Mie result.`);
  return {
    summary: `A ${inp.d_p} nm ${mat} sphere in a medium of index ${nm} has its extinction maximum at ${f4(lp)} nm with Q_ext = ${f3(qp.ext)} (σ_ext = ${f3(qp.ext * area)} m²), of which ${f3(100 * qp.sca / qp.ext)} % is scattering. The quasi-static Fröhlich condition falls at ${f4(lF)} nm; a spheroid of aspect ratio ${AR} has its longitudinal resonance at ${f4(lams[iL])} nm. Sensitivity to the surrounding index is ${f3(sens)} nm per refractive-index unit.`,
    kpis: [
      kpi('lam_peak', 'Extinction peak wavelength (Mie)', lp, 'nm', plasm ? 'Localised surface-plasmon resonance of the sphere from the full Mie series with the tabulated complex index; it red-shifts with size and with the index of the surroundings.' : 'Strongest Mie resonance of the dielectric sphere in the window.'),
      kpi('Qext', 'Peak extinction efficiency', qp.ext, '–', 'Extinction cross-section divided by the geometric cross-section πa²; values above 1 mean the particle intercepts more light than its shadow.'),
      kpi('sig_ext', 'Peak extinction cross-section', qp.ext * area, 'm²', 'Absolute extinction cross-section; multiplied by number density it gives the attenuation coefficient of a dilute suspension (Beer–Lambert).'),
      kpi('sig_sca', 'Scattering cross-section at peak', qp.sca * area, 'm²', 'Part of the extinction re-radiated; it scales as volume squared and dominates for large particles.'),
      kpi('sig_abs', 'Absorption cross-section at peak', qp.abs * area, 'm²', 'Part of the extinction converted to heat; it scales as volume and dominates for small particles (photothermal use).'),
      kpi('albedo', 'Single-scattering albedo at peak', qp.sca / qp.ext, '–', 'Q_sca/Q_ext; near 0 the particle is an absorber/heater, near 1 a scatterer/label.'),
      kpi('lam_F', 'Fröhlich wavelength (Re ε = −2ε_m)', lF, 'nm', 'Quasi-static dipole resonance condition of a small sphere; the Mie peak converges to it as the diameter tends to zero (zero if the condition is not met in the window).', lF > 0),
      kpi('lam_L', 'Longitudinal resonance of the spheroid (Gans)', lams[iL], 'nm', `Quasi-static resonance for light polarised along the long axis of a prolate spheroid with aspect ratio ${AR} (depolarisation factor ${f3(Lz)}); it red-shifts almost linearly with aspect ratio.`),
      kpi('FE_s', 'Field enhancement at sphere pole', feS, '–', 'Quasi-static |E|/|E₀| just outside the sphere along the polarisation, |3ε/(ε + 2ε_m)|, at the dipole resonance; intensity and SERS scale as its 2nd and 4th powers.'),
      kpi('FE_L', 'Field enhancement at spheroid tip', feL, '–', 'Quasi-static |E|/|E₀| at the tip of the spheroid at its longitudinal resonance (lightning-rod plus plasmon effect).'),
      kpi('sens', 'Refractive-index sensitivity', sens, 'nm/RIU', 'Shift of the Mie extinction peak per unit change of the surrounding refractive index, from a finite difference of 0.05.'),
      kpi('x', 'Size parameter at peak', x0, '–', '2πn_m a/λ; below about 0.3 the quasi-static picture holds, above 1 higher multipoles appear.'),
    ],
    series: [
      ser('Q', 'Mie efficiencies of the sphere', 'wavelength (nm)', 'efficiency Q (–)', [{ name: 'extinction', x: lams, y: Q.map(q => q.ext) }, { name: 'scattering', x: lams, y: Q.map(q => q.sca) }, { name: 'absorption', x: lams, y: Q.map(q => q.abs) }, { name: 'quasi-static dipole (extinction)', x: lams, y: qsS }], 'Full-series Mie efficiencies with the quasi-static dipole result for comparison. The dipole curve agrees for small particles; retardation red-shifts and broadens the true peak as size grows.'),
      ser('gans', 'Spheroid (Gans) extinction', 'wavelength (nm)', 'σ/πa² (–)', [{ name: 'longitudinal', x: lams, y: gL }, { name: 'transverse', x: lams, y: gT }, { name: 'orientation average', x: lams, y: gL.map((v, i) => (v + 2 * gT[i]) / 3) }], 'Quasi-static extinction of a prolate spheroid for the two polarisations and their random-orientation average. The longitudinal mode is the tunable nanorod resonance.'),
    ],
    equations: [...matEquations([mat]),
      eqn('governing', 'Mie series', 'Q_ext = (2/x²) Σ (2n + 1) Re(aₙ + bₙ),  Q_sca = (2/x²) Σ (2n + 1)(|aₙ|² + |bₙ|²),  x = 2πn_m a/λ'),
      eqn('formula', 'Quasi-static polarisability (Gans)', 'α_j = V (ε − ε_m)/[ε_m + L_j (ε − ε_m)];  sphere: L = 1/3 (Fröhlich: Re ε = −2ε_m)'),
      eqn('formula', 'Dipole cross-sections', 'σ_abs = k Im α,  σ_sca = k⁴|α|²/(6π)'), eqn('formula', 'Prolate depolarisation factor', 'L_z = ((1 − e²)/e²)·[ln((1 + e)/(1 − e))/(2e) − 1],  e² = 1 − (b/a)²')],
    checks: [
      chk('Optical theorem: Q_ext = Q_sca + Q_abs with Q_abs ≥ 0', Q.every(q => q.abs >= -1e-9), `min Q_abs = ${f3(N.minOf(Q.map(q => q.abs)))}`),
      chk('Rayleigh limit: Mie agrees with the dipole formula for a 2 nm sphere', (() => { const l = lams[ip], q = mieQ(2 * PI * nm * 1e-9 / (l * 1e-9), C.scale(nOf(mat, l, inp), 1 / nm)).ext, k = 2 * PI * nm / (l * 1e-9), al = C.scale(C.div(C.sub(epsOf(mat, l, inp), C.of(em)), C.add(epsOf(mat, l, inp), C.of(2 * em))), 4 * PI * 1e-27); return Math.abs(q - k * al.im / (PI * 1e-18)) < 0.02 * Math.abs(q) + 1e-9; })(), 'series and closed-form dipole extinction within 2 % at the peak wavelength'),
      chk('Series converged', true, `${qp.nstop} multipole orders (Wiscombe criterion x + 4x^(1/3) + 2)`), chk('Depolarisation factors sum to one', Math.abs(Lz + 2 * Lx - 1) < 1e-12, `L_z = ${f4(Lz)}, L_x = ${f4(Lx)}`),
      chk('Resonance inside the wavelength window', !edge, `peak at ${f4(lp)} nm in ${f3(lmin)}–${f3(inp.lam_max)} nm`),
    ],
    warnings: W,
    recommendations: [
      rec('high', qp.sca / qp.ext < 0.3 ? `This particle is mainly an absorber (albedo ${f3(qp.sca / qp.ext)}): suited to photothermal heating and colorimetric assays. For scattering labels or light trapping increase the diameter above about ${mat === 'Ag' ? 60 : 80} nm.` : `This particle scatters strongly (albedo ${f3(qp.sca / qp.ext)}): suited to dark-field labels and light trapping. For photothermal use choose diameters below about 40 nm.`, 'Scattering-to-extinction ratio from the Mie result'),
      rec('medium', `To tune the resonance into the near infrared use rods: each unit of aspect ratio shifts the longitudinal mode by roughly ${f3((lams[iL] - lams[iq]) / Math.max(AR - 1, 0.01))} nm here.`, 'Gans depolarisation-factor dependence'),
      rec('medium', `Keep the size distribution narrow: the peak moves with diameter and an index change of 0.01 shifts it by ${f3(0.01 * sens)} nm, which sets the resolution of a refractometric sensor.`, 'Computed refractive-index sensitivity'),
      rec('low', 'Confirm with a measured extinction spectrum and TEM size histogram; shape faceting and ligand shells shift the resonance by several nanometres.', 'Idealised smooth sphere in a homogeneous medium'),
      rec('medium', mat === 'Au' || mat === 'Ag' ? `Precious-metal use is tiny per particle but adds up: where the application allows, aluminium (UV) or doped oxides (near infrared) avoid ${mat === 'Au' ? 'gold' : 'silver'}, and silver ion release should be assessed for aquatic toxicity.` : 'Oxide and silicon particles are abundant and low-toxicity alternatives to noble-metal plasmonics for scattering applications.', 'Critical raw materials and nano-safety', true),
    ],
    publish: { 'res.f0': C0 / (lp * 1e-9), 'res.absorption': qp.abs / Math.max(qp.ext, 1e-30), 'geom.diameter': 2 * a, 'opt.sigma_ext': qp.ext * area, 'opt.lspr_nm': lp, 'opt.field_enhancement': feS },
    numerics: { method: 'Full Mie series (logarithmic-derivative downward recurrence); Gans quasi-static spheroid', dof: qp.nstop, steps: lams.length, notes: 'Single homogeneous sphere, bulk dielectric function, non-absorbing medium.' },
  };
}

// ---------------------------------------------------------------- effective-medium nanocomposites
function effmed(inp) {
  const W = [], mat = inp.mat_i, f = inp.fill, eh = C.of(inp.n_host * inp.n_host), lmin = Math.min(inp.lam_min, inp.lam_max - 10), lams = N.linspace(lmin, inp.lam_max, 240);
  dispWarn([mat], lmin, W);
  const mg = ei => { const b = C.div(C.sub(ei, eh), C.add(ei, C.scale(eh, 2))); return C.mul(eh, C.div(C.add(ONE, C.scale(b, 2 * f)), C.sub(ONE, C.scale(b, f)))); };
  const br = ei => { const b = C.add(C.scale(ei, 3 * f - 1), C.scale(eh, 2 - 3 * f)), disc = C.sqrt(C.add(C.mul(b, b), C.scale(C.mul(ei, eh), 8))); let e = C.scale(C.add(b, disc), 0.25); if (e.im < -1e-12 || (Math.abs(e.im) <= 1e-12 && e.re < 0)) e = C.scale(C.sub(b, disc), 0.25); return e; };
  const EM = lams.map(l => mg(epsOf(mat, l, inp))), EB = lams.map(l => br(epsOf(mat, l, inp))), l0 = N.clamp(inp.lam0, lmin, inp.lam_max), ei0 = epsOf(mat, l0, inp), m0 = mg(ei0), b0 = br(ei0), nM = csqrt(m0), nB = csqrt(b0);
  const ipk = N.argmax(EM.map(e => e.im)), lossless = Math.abs(ei0.im) < 1e-9 && ei0.re > 0, wU = f * ei0.re + (1 - f) * eh.re, wL = 1 / (f / ei0.re + (1 - f) / eh.re);
  const resB = C.abs(C.add(C.scale(C.div(C.sub(ei0, b0), C.add(ei0, C.scale(b0, 2))), f), C.scale(C.div(C.sub(eh, b0), C.add(eh, C.scale(b0, 2))), 1 - f)));
  if (f > 0.3) W.push('Maxwell-Garnett assumes dilute, well-separated inclusions; above about 30 vol % use the Bruggeman result and treat both as estimates.');
  if (MATS[mat].metal && f > 1 / 3) W.push('The Bruggeman model percolates at one-third filling: above it the composite behaves as a metal (Re ε < 0 at long wavelength).');
  return {
    summary: `${f3(100 * f)} vol % of ${mat} inclusions in a host of index ${inp.n_host} give, at ${f3(l0)} nm, ε_eff = ${cx(m0)} (Maxwell-Garnett) and ${cx(b0)} (Bruggeman), i.e. n_eff = ${cx(nM)} and ${cx(nB)}. The Maxwell-Garnett absorption resonance lies at ${f3(lams[ipk])} nm.`,
    kpis: [
      kpi('eps_mg', 'Re ε_eff (Maxwell-Garnett)', m0.re, '–', 'Real effective permittivity for isolated spherical inclusions in a continuous host (Clausius–Mossotti mixing).'),
      kpi('eps_mg_im', 'Im ε_eff (Maxwell-Garnett)', m0.im, '–', 'Effective absorption of the composite; it peaks at the collective plasmon resonance of the inclusions.', m0.im >= 0),
      kpi('eps_bg', 'Re ε_eff (Bruggeman)', b0.re, '–', 'Real effective permittivity of the symmetric (aggregate) topology, where neither phase is the matrix; it predicts percolation at one-third filling.'),
      kpi('eps_bg_im', 'Im ε_eff (Bruggeman)', b0.im, '–', 'Effective loss in the symmetric effective-medium approximation, typically broader in wavelength than the Maxwell-Garnett resonance.', b0.im >= 0),
      kpi('n_mg', 'Effective index (Maxwell-Garnett)', nM.re, '–', 'Real refractive index √ε_eff of the composite film on the passive branch.'),
      kpi('k_mg', 'Effective extinction index (Maxwell-Garnett)', nM.im, '–', 'Imaginary index; the absorption coefficient is 4πk/λ = ' + f3(4 * PI * nM.im / (l0 * 1e-9)) + ' m⁻¹.'),
      kpi('lam_res', 'Composite resonance wavelength', lams[ipk], 'nm', 'Wavelength of maximum Im ε_eff (Maxwell-Garnett); it red-shifts from the single-particle resonance as the filling fraction rises.'),
      kpi('spread', 'Model spread |ε_MG − ε_BG|', C.abs(C.sub(m0, b0)), '–', 'Difference between the two mixing rules at the probe wavelength; a large value means the result depends on microstructure, which must then be characterised.'),
    ],
    series: [
      ser('re', 'Real effective permittivity', 'wavelength (nm)', 'Re ε_eff (–)', [{ name: 'Maxwell-Garnett', x: lams, y: EM.map(e => e.re) }, { name: 'Bruggeman', x: lams, y: EB.map(e => e.re) }], 'Dispersion of the real effective permittivity for the two mixing rules. The Maxwell-Garnett curve shows a sharp resonance; the Bruggeman curve is smoother because of its distribution of local environments.'),
      ser('im', 'Imaginary effective permittivity', 'wavelength (nm)', 'Im ε_eff (–)', [{ name: 'Maxwell-Garnett', x: lams, y: EM.map(e => e.im) }, { name: 'Bruggeman', x: lams, y: EB.map(e => e.im) }], 'Absorption of the composite. The peak position gives the colour of metal-nanoparticle glasses and films.'),
    ],
    equations: [...matEquations([mat]), eqn('constitutive', 'Maxwell-Garnett', 'ε_eff = ε_h (1 + 2fβ)/(1 − fβ),  β = (ε_i − ε_h)/(ε_i + 2ε_h)'), eqn('constitutive', 'Bruggeman', 'f (ε_i − ε_eff)/(ε_i + 2ε_eff) + (1 − f)(ε_h − ε_eff)/(ε_h + 2ε_eff) = 0'), eqn('formula', 'Wiener bounds (lossless)', '[f/ε_i + (1 − f)/ε_h]⁻¹ ≤ ε_eff ≤ f ε_i + (1 − f) ε_h')],
    checks: [
      chk('Passivity: Im ε_eff ≥ 0 for both models', EM.every(e => e.im >= -1e-9) && EB.every(e => e.im >= -1e-9), 'correct root of the Bruggeman quadratic selected at every wavelength'),
      chk('Bruggeman equation residual', resB < 1e-9, `|residual| = ${f3(resB)}`), chk('Dilute limit: both rules tend to the host', C.abs(C.sub(mg(eh), eh)) < 1e-12, 'inclusion equal to host returns the host permittivity'),
      chk('Within Wiener bounds (lossless case)', !lossless || (m0.re <= wU + 1e-9 && m0.re >= wL - 1e-9 && b0.re <= wU + 1e-9 && b0.re >= wL - 1e-9), lossless ? `${f3(wL)} ≤ ε_eff ≤ ${f3(wU)}` : 'inclusion is lossy or metallic at the probe: bounds apply in the complex plane and are not checked'),
    ],
    warnings: W,
    recommendations: [
      rec('high', `The two mixing rules differ by ${f3(C.abs(C.sub(m0, b0)))} at ${f3(l0)} nm; ${C.abs(C.sub(m0, b0)) > 0.1 * C.abs(m0) ? 'characterise the microstructure (isolated or aggregated particles) by TEM before trusting either value' : 'the estimate is insensitive to topology at this loading'}.`, 'Topology dependence of effective-medium theory'),
      rec('medium', `Tune colour or absorption with filling fraction: the resonance is at ${f3(lams[ipk])} nm for f = ${f}; measure ellipsometric n and k of a test film to fit f.`, 'Maxwell-Garnett resonance shift'),
      rec('low', 'Effective-medium theory needs particles and spacings much smaller than the wavelength (below about λ/10); otherwise use Mie-based multiple-scattering models.', 'Quasi-static assumption'),
      rec('medium', `Only ${f3(100 * f)} vol % of ${mat} is needed for the optical function; keep metal loading at the minimum that meets the specification and prefer recoverable host matrices.`, 'Material efficiency of nanocomposite coatings', true),
    ],
    publish: { 'eff.eps': m0.re, 'eff.n': nM.re, 'res.absorption': N.clamp(1 - Math.exp(-4 * PI * nM.im * 100 / l0), 0, 1), 'res.f0': C0 / (lams[ipk] * 1e-9) },
    numerics: { method: 'Closed-form Maxwell-Garnett and Bruggeman (quadratic root) mixing', dof: lams.length, notes: 'Spherical inclusions, quasi-static limit, no size-dependent damping.' },
  };
}

// ---------------------------------------------------------------- hyperbolic / ENZ multilayers
function hyper(inp) {
  const W = [], mm = inp.mat_i, md = inp.mat_d, f = inp.fill, Lam = inp.period, lmin = Math.min(inp.lam_min, inp.lam_max - 10), lams = N.linspace(lmin, inp.lam_max, 300);
  dispWarn([mm, md], lmin, W);
  const par = l => { const em = epsOf(mm, l, inp), ed = epsOf(md, l, inp); return { p: C.add(C.scale(em, f), C.scale(ed, 1 - f)), t: C.div(C.mul(em, ed), C.add(C.scale(em, 1 - f), C.scale(ed, f))) }; };
  const cls = e => (e.p.re > 0 && e.t.re < 0 ? 'Type I hyperbolic' : e.p.re < 0 && e.t.re > 0 ? 'Type II hyperbolic' : e.p.re < 0 && e.t.re < 0 ? 'metallic' : 'elliptical dielectric');
  const Es = lams.map(par), l0 = N.clamp(inp.lam0, lmin, inp.lam_max), e0 = par(l0), kind = cls(e0), t1 = band(lams, Es.map(e => cls(e) === 'Type I hyperbolic')), t2 = band(lams, Es.map(e => cls(e) === 'Type II hyperbolic'));
  let lenz = 0; for (let i = 1; i < lams.length; i++) if (Es[i - 1].p.re * Es[i].p.re <= 0 && Es[i - 1].p.re !== Es[i].p.re) { lenz = lams[i - 1] + (lams[i] - lams[i - 1]) * Es[i - 1].p.re / (Es[i - 1].p.re - Es[i].p.re); break; }
  // exact Bloch check at normal incidence (wave sees ε∥)
  const mB = mul2(layerMat(nOf(mm, l0, inp), f * Lam, l0, 0, 's'), layerMat(nOf(md, l0, inp), (1 - f) * Lam, l0, 0, 's')), KB = cacos(C.scale(C.add(mB[0], mB[3]), 0.5)), k0L = 2 * PI * Lam / l0, nB = { re: KB.re / k0L, im: KB.im / k0L }, nE = csqrt(e0.p);
  const errB = Math.hypot(nB.re - nE.re, nB.im - nE.im) / Math.max(C.abs(nE), 1e-9), kmax = l0 / (2 * Lam);
  if (!MATS[mm].metal) W.push('The first constituent is not a metal: Re ε stays positive and no hyperbolic or ENZ regime can occur.');
  if (errB > 0.1) W.push(`Effective-medium and exact Bloch indices differ by ${f3(100 * errB)} % at the probe: the period is too large for local homogenisation.`);
  return {
    summary: `A ${mm}/${md} multilayer with ${f3(100 * f)} % metal and ${Lam} nm period has ε∥ = ${cx(e0.p)} and ε⊥ = ${cx(e0.t)} at ${f3(l0)} nm: ${kind}. ${lenz ? `The in-plane epsilon-near-zero wavelength is ${f3(lenz)} nm.` : 'ε∥ does not cross zero in the window.'}${t2 ? ` Type II dispersion spans ${f3(t2[0])}–${f3(t2[1])} nm.` : ''}${t1 ? ` Type I dispersion spans ${f3(t1[0])}–${f3(t1[1])} nm.` : ''}`,
    kpis: [
      kpi('class', 'Dispersion class at probe wavelength', kind, '', 'Sign pattern of the uniaxial permittivity tensor: Type I has one negative component (normal to the layers), Type II two (in the layer plane); the iso-frequency surface is then a hyperboloid supporting very large wavevectors.'),
      kpi('eps_par', 'Re ε∥ (in-plane)', e0.p.re, '–', 'Effective permittivity for fields parallel to the layers, the thickness-weighted arithmetic mean f·ε_m + (1 − f)·ε_d.', null),
      kpi('eps_perp', 'Re ε⊥ (normal)', e0.t.re, '–', 'Effective permittivity for fields normal to the layers, the harmonic mean; it has a pole where (1 − f)ε_m + f·ε_d = 0.', null),
      kpi('loss_par', 'Im ε∥', e0.p.im, '–', 'In-plane loss from the metal fraction; it limits propagation length of the high-k modes.', e0.p.im >= 0),
      kpi('lam_enz', 'In-plane ENZ wavelength', lenz, 'nm', 'Wavelength where Re ε∥ = 0, the boundary between elliptical and Type II hyperbolic behaviour; it is tuned by the metal fill fraction.', lenz > 0),
      kpi('t2_lo', 'Type II band lower edge', t2 ? t2[0] : 0, 'nm', 'Shortest wavelength in the window with ε∥ < 0 < ε⊥ (zero if absent).'), kpi('t1_lo', 'Type I band lower edge', t1 ? t1[0] : 0, 'nm', 'Shortest wavelength in the window with ε⊥ < 0 < ε∥ (zero if absent).'),
      kpi('kmax', 'High-k cut-off (k/k₀)', kmax, '–', 'Largest normalised wavevector the layered structure can support before the Brillouin-zone edge, λ/(2Λ); the effective-medium hyperbola is truncated there.'),
      kpi('bloch_err', 'Homogenisation error at probe', errB, '–', 'Relative difference between the effective-medium index √ε∥ and the exact Bloch index of the periodic stack at normal incidence.', errB < 0.05),
    ],
    series: [
      ser('eps', 'Effective permittivity tensor components', 'wavelength (nm)', 'Re ε (–)', [{ name: 'Re ε∥', x: lams, y: Es.map(e => N.clamp(e.p.re, -40, 40)) }, { name: 'Re ε⊥', x: lams, y: Es.map(e => N.clamp(e.t.re, -40, 40)) }, { name: 'zero', x: [lams[0], lams[299]], y: [0, 0] }], 'Real parts of the uniaxial effective permittivity (clipped at ±40). Regions where the two curves have opposite signs are hyperbolic; the zero crossing of ε∥ is the ENZ point and the pole of ε⊥ the epsilon-near-pole point.'),
      ser('loss', 'Loss of the effective medium', 'wavelength (nm)', 'Im ε (–)', [{ name: 'Im ε∥', x: lams, y: Es.map(e => pos(Math.min(e.p.im, 1e3))) }, { name: 'Im ε⊥', x: lams, y: Es.map(e => pos(Math.min(e.t.im, 1e3))) }], 'Imaginary parts of the tensor components. The peak of Im ε⊥ marks the epsilon-near-pole resonance; low Im ε∥ near the ENZ point is needed for field enhancement.', { logy: true }),
    ],
    equations: [...matEquations([mm, md]), eqn('constitutive', 'Layered effective medium', 'ε∥ = f ε_m + (1 − f) ε_d;  1/ε⊥ = f/ε_m + (1 − f)/ε_d'), eqn('governing', 'Extraordinary-wave dispersion', 'k∥²/ε⊥ + k⊥²/ε∥ = k₀²', 'hyperbolic when Re ε∥ · Re ε⊥ < 0'), eqn('governing', 'Exact Bloch check (normal incidence)', 'cos(KΛ) = ½ tr(M_m M_d),  n_Bloch = K/k₀')],
    checks: [
      chk('Passivity: Im ε∥ ≥ 0 and Im ε⊥ ≥ 0', Es.every(e => e.p.im >= -1e-12 && e.t.im >= -1e-12), 'both effective components are dissipative'), chk('Effective medium agrees with exact Bloch index', errB < 0.05, `√ε∥ = ${cx(nE)} against Bloch ${cx(nB)}`),
      chk('Deep-subwavelength period (Λ < λ/10)', Lam < l0 / 10, `Λ/λ = ${f3(Lam / l0)}`), chk('Hyperbolic regime available in the window', !!(t1 || t2), t1 || t2 ? 'at least one hyperbolic band found' : 'none'),
    ],
    warnings: W,
    recommendations: [
      rec('high', lenz ? `Shift the ENZ wavelength (${f3(lenz)} nm) by the fill fraction: more metal moves it to shorter wavelengths. For ENZ at ${f3(l0)} nm the required fraction is about ${f3((() => { const em = epsOf(mm, l0, inp).re, ed = epsOf(md, l0, inp).re; return em < 0 ? N.clamp(ed / (ed - em), 0, 1) : 0; })())}.` : 'Choose a metal whose Re ε is negative in the working band and a fill fraction near ε_d/(ε_d − ε_m) to place the ENZ point.', 'Zero of f·ε_m + (1 − f)·ε_d'),
      rec('medium', `Keep the period below about λ/20 (${f3(l0 / 20)} nm): the homogenisation error is ${f3(100 * errB)} % at ${Lam} nm and the usable high-k range ends at k/k₀ ≈ ${f3(kmax)}.`, 'Non-local corrections scale as (Λ/λ)²'),
      rec('medium', 'Continuous metal films below about 8–10 nm are hard to deposit; use wetting layers (Ge, Ti) or accept percolated-film optical constants in the design.', 'Fabrication limit of ultrathin metal layers'),
      rec('low', 'Silver gives the lowest loss but tarnishes: cap the stack with a dielectric. Transparent conducting oxides or nitrides shift the function to the near infrared without precious metals.', 'Durability and critical raw materials', true),
    ],
    publish: { 'eff.eps': e0.p.re, 'eff.n': nE.re, 'opt.eps_perp': e0.t.re, 'res.f0': lenz ? C0 / (lenz * 1e-9) : C0 / (l0 * 1e-9), 'geom.cell_size': Lam * 1e-9 },
    numerics: { method: 'Local effective-medium theory for a binary multilayer; exact Rytov/Bloch check by transfer matrices', dof: lams.length, notes: 'Infinite periodic stack; finite-stack and non-local effects only assessed through the Bloch comparison.' },
  };
}

// ---------------------------------------------------------------- metasurfaces
function metasurface(inp) {
  const W = [], lam = inp.lam0, ni = inp.n_inc, nt = inp.n_t, thi = inp.theta * PI / 180, p = inp.ms_p, Nst = Math.round(inp.ms_N), lens = inp.ms_kind === 'lens';
  const grad = 2 * PI / (Nst * p), st = (ni * Math.sin(thi) + lam / (Nst * p)) / nt, sr = Math.sin(thi) + lam / (ni * Nst * p), tht = Math.abs(st) <= 1 ? Math.asin(st) * 180 / PI : null, thr = Math.abs(sr) <= 1 ? Math.asin(sr) * 180 / PI : null;
  const snell = Math.asin(N.clamp(ni * Math.sin(thi) / nt, -1, 1)) * 180 / PI, pmax = lam / (ni * Math.abs(Math.sin(thi)) + Math.max(ni, nt)), eff = (Math.sin(PI / Nst) / (PI / Nst)) ** 2, Hmin = lam / Math.max(inp.n_pillar - Math.min(ni, nt), 0.05);
  // metalens
  const Rl = inp.lens_D * 0.5e3, fl = inp.lens_f * 1e3, NA = nt * Math.sin(Math.atan(Rl / fl)), airy = 0.61 * lam / NA, fwhm = 0.514 * lam / NA, zones = nt * (Math.hypot(Rl, fl) - fl) / lam, dof = lam * nt / (NA * NA), zoneW = lam / NA, pNyq = lam / (2 * NA);
  const dfc = fl * inp.dlam / lam, bwAch = lam * dof / fl, phi = r => -(2 * PI * nt / lam) * (Math.hypot(r, fl) - fl), wrap = v => ((v % (2 * PI)) + 2 * PI) % (2 * PI);
  if (p > pmax) W.push(`The unit-cell period ${p} nm exceeds ${f3(pmax)} nm: the lattice itself diffracts into propagating orders and efficiency is lost.`);
  if (!lens && tht === null) W.push('The phase gradient is steeper than the transmission medium can support: the anomalous beam is evanescent (converted to a surface wave).');
  if (lens && p > pNyq) W.push(`The period exceeds the Nyquist limit λ/(2NA) = ${f3(pNyq)} nm: the outer zones of the lens are under-sampled and focusing efficiency drops at the rim.`);
  if (Nst < 3) W.push('Fewer than three phase levels cannot distinguish the +1 from the −1 order (a binary grating splits power symmetrically).');
  const polys = [], oct = (cx, cy, d) => ({ layer: 1, pts: Array.from({ length: 8 }, (_, k) => [+(cx + 0.5 * d * Math.cos(PI / 8 + k * PI / 4)).toFixed(2), +(cy + 0.5 * d * Math.sin(PI / 8 + k * PI / 4)).toFixed(2)]) }), dOf = ph => p * (0.25 + 0.5 * wrap(ph) / (2 * PI));
  if (lens) { const nx = Math.min(Math.floor(Rl / p), 60); for (let i = 0; i <= nx; i++) for (let j = -3; j <= 3; j++) polys.push(oct(i * p, j * p, dOf(phi(Math.hypot(i * p, j * p))))); }
  else { const nx = Math.min(3 * Nst, 30); for (let i = 0; i < nx; i++) for (let j = 0; j < 8; j++) polys.push(oct(i * p, j * p, dOf(2 * PI * (i % Nst) / Nst))); }
  const rs = N.linspace(0, Rl, 300), ths = N.linspace(-80, 80, 161), u = N.linspace(0, 3 * airy, 120), lw = N.linspace(Math.max(lam - Math.max(inp.dlam, 1), 0.4 * lam), lam + Math.max(inp.dlam, 1), 41), okT = ths.filter(t => Math.abs((ni * Math.sin(t * PI / 180) + lam / (Nst * p)) / nt) <= 1), okS = ths.filter(t => Math.abs(ni * Math.sin(t * PI / 180) / nt) <= 1);
  const j1 = x => N.besselJ(1, x), series = lens ? [
    ser('phase', 'Metalens phase profile', 'radius (µm)', 'phase (rad)', [{ name: 'wrapped to 0–2π', x: rs.map(r => r / 1000), y: rs.map(r => wrap(phi(r))) }], 'Hyperbolic phase profile required for aberration-free focusing of a normally incident plane wave, wrapped modulo 2π. Zones become narrower towards the rim; the outermost zone width sets the unit-cell period needed.'),
    ser('spot', 'Diffraction-limited focal spot', 'radial position (nm)', 'normalised intensity (–)', [{ name: 'Airy pattern', x: u, y: u.map(r => { const x = 2 * PI * NA * r / lam; return x < 1e-6 ? 1 : (2 * j1(x) / x) ** 2; }) }], 'Ideal focal-plane intensity [2J₁(x)/x]² for the computed numerical aperture. The first zero is the Airy radius; a real metalens approaches this when phase sampling and efficiency are adequate.'),
    ser('chrom', 'Chromatic focal shift', 'wavelength (nm)', 'focal length (µm)', [{ name: 'f(λ) = f₀λ₀/λ', x: lw, y: lw.map(l => fl * lam / l / 1000) }], 'Focal length against wavelength for a fixed (diffractive) phase profile: the focus moves inversely with wavelength, the opposite sign and several times the magnitude of refractive-lens dispersion.'),
  ] : [
    ser('snell', 'Generalised refraction', 'incidence angle (°)', 'transmission angle (°)', [...(okT.length ? [{ name: 'with phase gradient', x: okT, y: okT.map(t => Math.asin((ni * Math.sin(t * PI / 180) + lam / (Nst * p)) / nt) * 180 / PI) }] : []), { name: 'ordinary Snell law', x: okS.length ? okS : [0], y: okS.length ? okS.map(t => Math.asin(ni * Math.sin(t * PI / 180) / nt) * 180 / PI) : [0] }], 'Transmission angle against incidence angle with and without the interfacial phase gradient. The offset between the curves is the anomalous deflection; where the curve ends the anomalous beam becomes evanescent.'),
    ser('steps', 'Phase steps across the supercell', 'position (nm)', 'phase (rad)', [{ name: 'discretised phase', x: Array.from({ length: 2 * Nst }, (_, i) => (Math.floor(i / 2) + (i % 2)) * p), y: Array.from({ length: 2 * Nst }, (_, i) => 2 * PI * Math.floor(i / 2) / Nst) }, { name: 'ideal linear phase', x: [0, Nst * p], y: [0, 2 * PI] }], 'Staircase approximation of the linear 0–2π phase ramp by the discrete elements of one supercell. More steps approach the ideal blaze and raise the first-order efficiency as sinc²(1/N).'),
  ];
  const d1 = Math.asin(N.clamp((ni * Math.sin(thi) + lam / (Nst * p)) / nt, -1, 1)) * 180 / PI;
  return {
    summary: lens ? `A ${inp.lens_D} µm metalens with ${inp.lens_f} µm focal length at ${lam} nm has NA = ${f3(NA)}, an Airy radius of ${f3(airy)} nm (FWHM ${f3(fwhm)} nm), ${f3(zones)} Fresnel zones and a depth of focus of ${f3(dof / 1000)} µm. A ±${inp.dlam} nm bandwidth moves the focus by ±${f3(dfc / 1000)} µm. The unit cell must be below ${f3(Math.min(pNyq, pmax))} nm; ${p} nm is ${p <= Math.min(pNyq, pmax) ? 'adequate' : 'too coarse'}.`
      : `A ${Nst}-step phase gradient of 2π per ${f3(Nst * p)} nm deflects ${lam} nm light incident at ${inp.theta}° to ${tht === null ? 'no propagating beam' : f3(tht) + '°'} in transmission (ordinary refraction: ${f3(snell)}°) and ${thr === null ? 'no propagating beam' : f3(thr) + '°'} in reflection. Phase quantisation limits the first-order efficiency to ${f3(100 * eff)} %; the unit cell must be below ${f3(pmax)} nm to suppress lattice diffraction.`,
    kpis: [
      kpi('theta_t', 'Anomalous refraction angle', tht ?? 0, '°', 'Transmission angle from the generalised Snell law n_t sin θ_t − n_i sin θ_i = (λ/2π) dφ/dx; zero is reported when the beam is evanescent.', tht !== null),
      kpi('theta_r', 'Anomalous reflection angle', thr ?? 0, '°', 'Reflection angle from sin θ_r − sin θ_i = (λ/2πn_i) dφ/dx; it differs from the incidence angle by the phase-gradient momentum.', thr !== null),
      kpi('grad', 'Phase gradient', grad * 1000, 'rad/µm', 'Imposed interfacial phase gradient dφ/dx = 2π/(N·p) of the supercell.'),
      kpi('eff', 'Quantisation-limited efficiency', eff, '–', `Fraction of power a perfectly transmitting ${Nst}-level phase staircase sends into the desired order, sinc²(1/N); element transmission and coupling losses reduce it further.`, eff > 0.8),
      kpi('pmax', 'Maximum unit-cell period', pmax, 'nm', 'Largest lattice period with only the zeroth lattice order propagating in both media, λ/(n_i|sin θ_i| + max(n_i, n_t)).', p <= pmax),
      kpi('Hmin', 'Minimum pillar height for 2π coverage', Hmin, 'nm', `Lower bound λ/(n_pillar − n_surround) for waveguide-type pillars of index ${inp.n_pillar}; practical libraries need 20–50 % more because the fill-factor range is limited.`),
      kpi('NA', 'Metalens numerical aperture', NA, '–', 'n_t·sin(atan(D/2f)); it fixes the resolution and the fineness of the outer zones.'),
      kpi('airy', 'Airy radius of the focal spot', airy, 'nm', 'Radius of the first dark ring, 0.61λ/NA: the diffraction-limited resolution of the lens.'),
      kpi('zones', 'Number of 2π zones', zones, '–', 'Phase wraps between centre and rim, n_t(√(R² + f²) − f)/λ; it also measures how dispersive the lens is.'),
      kpi('dof', 'Depth of focus', dof / 1000, 'µm', 'Axial range n_tλ/NA² over which the focus stays near its peak intensity.'),
      kpi('dfc', 'Chromatic focal shift', dfc / 1000, 'µm', `Focal-length change for a wavelength offset of ${inp.dlam} nm, |Δf| = f·Δλ/λ for a diffractive lens; compare with the depth of focus.`, dfc < dof),
      kpi('bw', 'Achromatic bandwidth estimate', bwAch, 'nm', 'Wavelength offset at which the focal shift equals the depth of focus, λ·DOF/f; a rough usable half-bandwidth without dispersion engineering.'),
    ],
    series,
    geometry: { kind: 'layout', polygons: polys, units: 'nm', title: lens ? 'Central radial strip of the metalens pillar pattern' : 'Phase-gradient supercells (pillar pattern)' },
    equations: [
      eqn('governing', 'Generalised Snell law of refraction', 'n_t sin θ_t − n_i sin θ_i = (λ₀/2π) dφ/dx'), eqn('governing', 'Generalised law of reflection', 'sin θ_r − sin θ_i = (λ₀/2πn_i) dφ/dx'),
      eqn('formula', 'Hyperbolic metalens phase', 'φ(r) = −(2π n_t/λ)(√(r² + f²) − f)'), eqn('formula', 'Resolution and focus', 'NA = n_t sin(atan(D/2f)),  r_Airy = 0.61λ/NA,  DOF = n_tλ/NA²,  f(λ) = f₀λ₀/λ'),
      eqn('formula', 'Sampling and diffraction limits', 'p < λ/(2NA)  (lens);  p < λ/(n_i|sin θ_i| + max(n_i, n_t))  (no lattice orders)'),
      eqn('formula', 'Pillar material', (() => { try { return `${chem.prettyFormula('TiO2')}, M = ${f4(chem.molarMass('TiO2'))} g/mol (typical visible-range dielectric; set n_pillar for Si, GaN or SiN)`; } catch { return 'dielectric pillar'; } })()),
    ],
    checks: [
      chk('Anomalous beam equals the +1 order of the supercell grating', Math.abs(st) > 1 || Math.abs(d1 - tht) < 1e-9, 'generalised Snell law and the grating equation with period N·p give the same angle'),
      chk('No lattice diffraction orders (p < p_max)', p <= pmax, `p = ${p} nm, p_max = ${f3(pmax)} nm`), chk(lens ? 'Nyquist sampling of the outer zone' : 'At least three phase levels', lens ? p <= pNyq : Nst >= 3, lens ? `p = ${p} nm against λ/2NA = ${f3(pNyq)} nm (outer zone width ${f3(zoneW)} nm)` : `${Nst} levels`),
      chk('Propagating anomalous beam', lens || tht !== null, lens ? 'not applicable to the lens' : tht === null ? 'evanescent' : `θ_t = ${f3(tht)}°`), chk('Chromatic shift within depth of focus', !lens || dfc <= dof, `|Δf| = ${f3(dfc / 1000)} µm, DOF = ${f3(dof / 1000)} µm`),
    ],
    warnings: W,
    recommendations: [
      rec('high', lens ? `Use a unit cell of at most ${f3(Math.min(pNyq, pmax))} nm for NA = ${f3(NA)}; pillar height must exceed about ${f3(1.3 * Hmin)} nm for a ${inp.n_pillar}-index material, i.e. aspect ratios near ${f3(1.3 * Hmin / (0.3 * p))} for the thinnest pillars.` : `Use ${Nst >= 6 ? 'the present' : 'six to eight'} phase levels (efficiency limit ${f3(100 * (Math.sin(PI / Math.max(Nst, 6)) / (PI / Math.max(Nst, 6))) ** 2)} %) and keep p below ${f3(pmax)} nm.`, 'Sampling theorem and waveguide-phase pillar design'),
      rec('medium', lens ? `The lens is usable over about ±${f3(bwAch)} nm without dispersion engineering; for broadband imaging reduce the diameter/zone count or use coupled-pillar achromatic elements.` : 'Efficiency falls steeply for deflection beyond about 40°: impedance mismatch requires non-local or bianisotropic designs at large angles.', 'Diffractive dispersion and large-angle limits'),
      rec('medium', 'Replace the linear phase-to-diameter mapping of the exported layout with a simulated (RCWA/FDTD) pillar library for the chosen height before lithography; the layout shown illustrates the phase pattern only.', 'Illustrative geometry generation'),
      rec('low', 'A single lithography step on a thin TiO₂, SiN or Si film replaces ground and polished glass; choose nano-imprint replication for volume to cut energy, material and chemical use per lens.', 'Flat-optics manufacturing footprint', true),
    ],
    publish: { 'opt.NA': NA, 'opt.spot_radius': airy * 1e-9, 'opt.deflection_deg': tht ?? 0, 'geom.cell_size': p * 1e-9, 'geom.thickness': Hmin * 1e-9, 'res.f0': C0 / (lam * 1e-9) },
    numerics: { method: 'Closed-form generalised Snell law, scalar diffraction (Airy) and sampling criteria; geometric layout generation', dof: polys.length, notes: 'Ideal phase elements with unit transmission; no full-wave unit-cell solution.' },
  };
}

const M = (...a) => ({ mode: a }), opts = (keys) => keys.map(k => ({ value: k, label: MATS[k].label }));
const DIEL = ['TiO2', 'SiO2', 'Si', 'ZnO', 'glass', 'air', 'custom'], ALL = ['TiO2', 'SiO2', 'Si', 'ZnO', 'glass', 'Au', 'Ag', 'Al', 'ITO', 'air', 'custom'];
export default {
  id: 12,
  tagline: 'Thin-film and photonic-crystal spectra, negative permittivity, permeability and refractive index, nanoparticle plasmonics, effective media and metasurfaces.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'multilayer', group: 'Model', help: 'Selects the electromagnetic model that is solved; each analysis has its own inputs.', options: [
      { value: 'multilayer', label: 'Multilayer thin film / 1-D photonic crystal (transfer matrix)' }, { value: 'neg_eps', label: 'Negative permittivity (Drude metal, wire medium)' }, { value: 'neg_mu', label: 'Negative permeability (split-ring resonators)' },
      { value: 'nim', label: 'Negative refractive index (composite, S-parameter retrieval)' }, { value: 'mie', label: 'Nanoparticle optics (Mie, Gans)' }, { value: 'effmed', label: 'Nanocomposite effective permittivity (Maxwell-Garnett, Bruggeman)' },
      { value: 'hyperbolic', label: 'Hyperbolic / ENZ multilayer (effective medium)' }, { value: 'metasurface', label: 'Metasurface: phase gradient and metalens' }] },
    { key: 'lam0', label: 'Design / probe wavelength', unit: 'nm', value: 550, min: 250, max: 5000, group: 'Illumination', showIf: M('multilayer', 'effmed', 'hyperbolic', 'metasurface'), help: 'Wavelength at which quarter-wave thicknesses, field profiles and single-wavelength results are evaluated.' },
    { key: 'lam_min', label: 'Spectrum start', unit: 'nm', value: 350, min: 200, max: 4000, group: 'Illumination', showIf: M('multilayer', 'mie', 'effmed', 'hyperbolic'), help: 'Short-wavelength end of the computed spectrum.' },
    { key: 'lam_max', label: 'Spectrum end', unit: 'nm', value: 900, min: 300, max: 6000, group: 'Illumination', showIf: M('multilayer', 'mie', 'effmed', 'hyperbolic'), help: 'Long-wavelength end of the computed spectrum.' },
    { key: 'theta', label: 'Angle of incidence', unit: '°', value: 0, min: 0, max: 89, group: 'Illumination', showIf: M('multilayer', 'nim', 'metasurface'), help: 'Angle from the surface normal in the incident medium. It blue-shifts stop bands and enters the Snell laws.' },
    { key: 'pol', label: 'Polarisation', type: 'select', value: 's', group: 'Illumination', showIf: M('multilayer'), help: 's (TE) has the electric field parallel to the layers, p (TM) in the plane of incidence; they differ only at oblique incidence.', options: [{ value: 's', label: 's (TE)' }, { value: 'p', label: 'p (TM)' }, { value: 'u', label: 'Unpolarised (average)' }] },
    { key: 'n_inc', label: 'Incident medium index', unit: '–', value: 1, min: 1, max: 4, group: 'Illumination', showIf: M('multilayer', 'metasurface'), help: 'Real refractive index of the medium the light comes from (1 for air, 1.33 water, 1.5 glass).' },
    { key: 'stack', label: 'Stack design', type: 'select', value: 'bragg', group: 'Multilayer', showIf: M('multilayer'), help: 'Periodic (H L)^N Bragg mirror, a single film of the high-index material, or a Fabry–Perot microcavity (H L)^N 2L (L H)^N.', options: [{ value: 'bragg', label: 'Periodic stack / Bragg mirror (HL)^N' }, { value: 'film', label: 'Single film (material H)' }, { value: 'cavity', label: 'Microcavity (HL)^N 2L (LH)^N' }] },
    { key: 'matH', label: 'Material H (first layer)', type: 'select', value: 'TiO2', group: 'Multilayer', showIf: M('multilayer'), help: 'High-index or metallic layer adjacent to the incident medium.', options: opts(ALL) },
    { key: 'matL', label: 'Material L', type: 'select', value: 'SiO2', group: 'Multilayer', showIf: M('multilayer'), help: 'Second (low-index) layer of the period.', options: opts(ALL) },
    { key: 'dH', label: 'Thickness of H (0 = quarter wave)', unit: 'nm', value: 0, min: 0, max: 2000, group: 'Multilayer', showIf: M('multilayer'), help: 'Physical thickness of the H layer; zero selects λ₀/(4n) automatically.' },
    { key: 'dL', label: 'Thickness of L (0 = quarter wave)', unit: 'nm', value: 0, min: 0, max: 2000, group: 'Multilayer', showIf: M('multilayer'), help: 'Physical thickness of the L layer; zero selects λ₀/(4n) automatically.' },
    { key: 'Nper', label: 'Number of periods', unit: '–', value: 8, min: 1, max: 30, step: 1, group: 'Multilayer', showIf: M('multilayer'), help: 'Number of H/L pairs (per mirror for the microcavity). Reflectance in the stop band approaches 1 exponentially with it.' },
    { key: 'substrate', label: 'Substrate', type: 'select', value: 'glass', group: 'Multilayer', showIf: M('multilayer'), help: 'Semi-infinite exit medium.', options: opts(['glass', 'SiO2', 'Si', 'ZnO', 'TiO2', 'air', 'Au', 'Ag', 'Al', 'custom']) },
    { key: 'n_user', label: 'Custom material n', unit: '–', value: 2, min: 0.05, max: 6, group: 'Custom material', showIf: M('multilayer', 'mie', 'effmed', 'hyperbolic'), help: 'Real refractive index used wherever the "custom" material is selected (non-dispersive).', bus: 'prop.n_refr' },
    { key: 'k_user', label: 'Custom material k', unit: '–', value: 0, min: 0, max: 10, group: 'Custom material', showIf: M('multilayer', 'mie', 'effmed', 'hyperbolic'), help: 'Extinction coefficient of the custom material (0 for lossless).' },
    { key: 'medium', label: 'Negative-ε medium', type: 'select', value: 'wire', group: 'Negative permittivity', showIf: M('neg_eps'), help: 'Pendry thin-wire lattice (microwave plasma) or a Drude metal with typical fitted parameters.', options: [{ value: 'wire', label: 'Thin-wire medium (Pendry)' }, { value: 'Ag', label: 'Silver (Drude)' }, { value: 'Au', label: 'Gold (Drude)' }, { value: 'Al', label: 'Aluminium (Drude)' }, { value: 'ITO', label: 'Indium tin oxide (Drude)' }] },
    { key: 'wire_model', label: 'Wire-lattice formula', type: 'select', value: 'lattice', group: 'Negative permittivity', showIf: { mode: ['neg_eps'], medium: ['wire'] }, help: 'Dynamic lattice sum: valid from thin to thick wires (r/a up to about 0.2). Quasi-static formula ln(a/2πr) + 0.5275: thin wires. Pendry: original logarithmic estimate, 5–30 % low unless r/a is extremely small.', options: [{ value: 'lattice', label: 'Dynamic lattice sum (thin and thick wires)' }, { value: 'belov', label: 'Quasi-static lattice sum (Belov–Maslovski–Tretyakov)' }, { value: 'pendry', label: 'Pendry thin-wire formula' }] },
    { key: 'wire_r', label: 'Wire radius', unit: 'µm', value: 1, min: 0.01, max: 2000, log: true, group: 'Negative permittivity', showIf: M('neg_eps'), help: 'Radius of the wires. The plasma frequency depends on it only logarithmically, the loss inversely with its square.' },
    { key: 'wire_a', label: 'Wire lattice period', unit: 'mm', value: 5, min: 0.001, max: 100, log: true, group: 'Negative permittivity', showIf: M('neg_eps'), help: 'Spacing of the wire lattice; the plasma frequency scales nearly as 1/a.' },
    { key: 'metal_sigma', label: 'Metal conductivity', unit: 'S/m', value: 3.7e7, min: 1e4, max: 7e7, log: true, group: 'Metamaterial metal', showIf: M('neg_eps', 'neg_mu'), bus: 'prop.sigma_elec', help: 'DC conductivity of the wires or rings (3.7×10⁷ aluminium, 5.8×10⁷ copper); it sets the damping.' },
    { key: 'f_rel', label: 'Probe frequency (relative)', unit: '–', value: 0.6, min: 0.05, max: 3, group: 'Metamaterial metal', showIf: M('neg_eps', 'neg_mu'), help: 'Frequency for single-point results as a multiple of the ENZ frequency (negative ε) or of the ring resonance (negative µ; try 1.05).' },
    { key: 'srr_type', label: 'Ring model', type: 'select', value: 'planar', group: 'Negative permeability', showIf: M('neg_mu'), help: 'Planar: printed (edge-coupled) double split rings on a circuit board, treated as an LC circuit with the flat-ring inductance and the coplanar-strip capacitance between the rings. Cylinder: Pendry\'s formula for split cylinders or rings stacked much closer than their radius; for isolated printed rings it overestimates the resonance several-fold.', options: [{ value: 'planar', label: 'Planar printed rings (quasi-static LC model)' }, { value: 'cylinder', label: 'Split cylinders / densely stacked rings (Pendry formula)' }] },
    { key: 'srr_epsr', label: 'Substrate relative permittivity', unit: '–', value: 4.4, min: 1, max: 20, group: 'Negative permeability', showIf: { mode: ['neg_mu'], srr_type: ['planar'] }, help: 'Permittivity of the board carrying the rings (about 4.4 for FR-4/G-10, 2.2 for PTFE laminates, 1 for free-standing rings). It raises the inter-ring capacitance and lowers the resonance.' },
    { key: 'srr_h', label: 'Substrate thickness', unit: 'mm', value: 0.25, min: 0.001, max: 20, log: true, group: 'Negative permeability', showIf: { mode: ['neg_mu'], srr_type: ['planar'] }, help: 'Thickness of the board. A board thinner than the ring width loads the gap only partly; a thick board gives the half-space value (ε_r + 1)/2.' },
    { key: 'srr_r', label: 'Ring radius', unit: 'mm', value: 1.5, min: 0.0005, max: 50, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'Inner radius of the split-ring pair; the resonance frequency scales as r^(−3/2).' },
    { key: 'srr_a', label: 'Lattice period', unit: 'mm', value: 8, min: 0.001, max: 200, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'In-plane lattice constant; with the ring size it gives the strength (filling factor) of the magnetic resonance; it must exceed the outer ring diameter 2(r + 2c + d) for planar rings.' },
    { key: 'srr_c', label: 'Ring width', unit: 'mm', value: 0.8, min: 0.0001, max: 20, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'Radial width of each ring; enters the inter-ring capacitance and the resistance.' },
    { key: 'srr_d', label: 'Gap between rings', unit: 'mm', value: 0.2, min: 0.00001, max: 5, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'Radial gap between the two concentric rings; a smaller gap raises capacitance and lowers the resonance.' },
    { key: 'srr_l', label: 'Stacking distance', unit: 'mm', value: 8, min: 0.0005, max: 50, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'Distance between ring planes along the axis. The planar model treats the rings as magnetically independent (spacing of the order of the ring radius or more); the cylinder formula needs spacing much smaller than the radius.' },
    { key: 'srr_t', label: 'Metal thickness', unit: 'µm', value: 35, min: 0.01, max: 500, log: true, group: 'Negative permeability', showIf: M('neg_mu'), help: 'Thickness of the ring metallisation; only about one skin depth is effective.' },
    { key: 'fep', label: 'Electric plasma frequency', unit: 'GHz', value: 12.8, min: 0.1, max: 1000, log: true, group: 'Negative index', showIf: M('nim'), help: 'Plasma frequency of the wire (ε) sub-lattice; must exceed the magnetic resonance for a double-negative band.' },
    { key: 'ge', label: 'Electric damping γ_e/2π', unit: 'GHz', value: 0.1, min: 0.0001, max: 100, log: true, group: 'Negative index', showIf: M('nim'), help: 'Collision frequency of the wire medium.' },
    { key: 'f0m', label: 'Magnetic resonance frequency', unit: 'GHz', value: 10, min: 0.1, max: 1000, log: true, group: 'Negative index', showIf: M('nim'), help: 'Resonance of the split-ring (µ) sub-lattice.' },
    { key: 'Fm', label: 'Magnetic filling factor', unit: '–', value: 0.3, min: 0.01, max: 0.9, group: 'Negative index', showIf: M('nim'), help: 'Strength F of the magnetic resonance; the µ < 0 band spans f₀ to f₀/√(1 − F).' },
    { key: 'gm', label: 'Magnetic damping γ_m/2π', unit: 'GHz', value: 0.05, min: 0.0001, max: 100, log: true, group: 'Negative index', showIf: M('nim'), help: 'Linewidth of the ring resonance; dominant loss in the negative-index band.' },
    { key: 'd_slab', label: 'Slab thickness', unit: 'mm', value: 10, min: 0.01, max: 50, log: true, group: 'Negative index', showIf: M('nim'), help: 'Thickness of the homogenised slab used for S-parameters and retrieval.' },
    { key: 'f_probe', label: 'Probe frequency', unit: 'GHz', value: 11, min: 0.1, max: 1000, log: true, group: 'Negative index', showIf: M('nim'), help: 'Frequency at which n, z, ε, µ and the refraction angle are reported.' },
    { key: 'mat_p', label: 'Particle material', type: 'select', value: 'Au', group: 'Nanoparticle', showIf: M('mie'), help: 'Material of the sphere/spheroid with its built-in dispersive permittivity.', options: opts(['Au', 'Ag', 'Al', 'Si', 'TiO2', 'ZnO', 'SiO2', 'ITO', 'custom']) },
    { key: 'd_p', label: 'Particle diameter', unit: 'nm', value: 40, min: 2, max: 1000, log: true, group: 'Nanoparticle', showIf: M('mie'), help: 'Sphere diameter (minor axis of the spheroid). Larger particles scatter more and red-shift.' },
    { key: 'AR_p', label: 'Spheroid aspect ratio', unit: '–', value: 3, min: 1, max: 10, group: 'Nanoparticle', showIf: M('mie'), bus: 'geom.aspect', help: 'Length-to-width ratio for the quasi-static nanorod (Gans) estimate.' },
    { key: 'n_med', label: 'Surrounding medium index', unit: '–', value: 1.33, min: 1, max: 2.5, group: 'Nanoparticle', showIf: M('mie'), help: 'Real refractive index of the embedding medium (1 air, 1.33 water, 1.5 glass/polymer).' },
    { key: 'mat_i', label: 'Inclusion / metal layer material', type: 'select', value: 'Ag', group: 'Effective medium', showIf: M('effmed', 'hyperbolic'), help: 'Material of the inclusions (nanocomposite) or of the metal layers (multilayer metamaterial).', options: opts(['Ag', 'Au', 'Al', 'ITO', 'Si', 'TiO2', 'ZnO', 'custom']) },
    { key: 'fill', label: 'Volume / thickness fraction', unit: '–', value: 0.2, min: 0.005, max: 0.95, group: 'Effective medium', showIf: M('effmed', 'hyperbolic'), help: 'Volume fraction of inclusions, or metal thickness fraction of the multilayer period.' },
    { key: 'n_host', label: 'Host refractive index', unit: '–', value: 1.5, min: 1, max: 4, group: 'Effective medium', showIf: M('effmed'), help: 'Real index of the host matrix of the nanocomposite.' },
    { key: 'mat_d', label: 'Dielectric layer material', type: 'select', value: 'TiO2', group: 'Effective medium', showIf: M('hyperbolic'), help: 'Dielectric constituent of the metal–dielectric multilayer.', options: opts(DIEL) },
    { key: 'period', label: 'Multilayer period', unit: 'nm', value: 30, min: 2, max: 500, log: true, group: 'Effective medium', showIf: M('hyperbolic'), help: 'Thickness of one metal + dielectric pair; must be far below the wavelength for homogenisation.' },
    { key: 'ms_kind', label: 'Metasurface function', type: 'select', value: 'lens', group: 'Metasurface', showIf: M('metasurface'), help: 'Constant phase gradient (beam deflector) or hyperbolic phase profile (lens).', options: [{ value: 'lens', label: 'Metalens' }, { value: 'gradient', label: 'Phase-gradient deflector' }] },
    { key: 'ms_p', label: 'Unit-cell period', unit: 'nm', value: 300, min: 50, max: 5000, log: true, group: 'Metasurface', showIf: M('metasurface'), help: 'Centre-to-centre spacing of the phase elements.' },
    { key: 'ms_N', label: 'Phase steps per 2π', unit: '–', value: 8, min: 2, max: 32, step: 1, group: 'Metasurface', showIf: M('metasurface'), help: 'Number of discrete elements spanning one 2π phase ramp (deflector supercell).' },
    { key: 'n_t', label: 'Transmission medium index', unit: '–', value: 1, min: 1, max: 4, group: 'Metasurface', showIf: M('metasurface'), help: 'Refractive index on the exit side, where the beam is deflected or focused.' },
    { key: 'n_pillar', label: 'Pillar refractive index', unit: '–', value: 2.4, min: 1.4, max: 4.5, group: 'Metasurface', showIf: M('metasurface'), help: 'Index of the dielectric pillars (about 2.4 TiO₂, 2.0 SiN, 3.5–4 Si).' },
    { key: 'lens_D', label: 'Lens diameter', unit: 'µm', value: 100, min: 2, max: 20000, log: true, group: 'Metasurface', showIf: M('metasurface'), help: 'Aperture diameter of the metalens.' },
    { key: 'lens_f', label: 'Focal length', unit: 'µm', value: 100, min: 1, max: 1e5, log: true, group: 'Metasurface', showIf: M('metasurface'), help: 'Design focal length at the design wavelength.' },
    { key: 'dlam', label: 'Bandwidth (half-width)', unit: 'nm', value: 10, min: 0, max: 500, group: 'Metasurface', showIf: M('metasurface'), help: 'Half-width of the source spectrum used for the chromatic focal shift.' },
  ],
  presets: [
    { name: 'TiO₂/SiO₂ quarter-wave Bragg mirror at 550 nm', note: 'Eight-pair dielectric 1-D photonic crystal on glass; stop band centred on the design wavelength.', values: { mode: 'multilayer', stack: 'bragg', matH: 'TiO2', matL: 'SiO2', Nper: 8, lam0: 550, lam_min: 350, lam_max: 900, theta: 0, dH: 0, dL: 0 } },
    { name: 'Gold nanoparticles (40 nm) in water', note: 'Localised surface plasmon resonance near 520–530 nm; classic colloidal gold.', values: { mode: 'mie', mat_p: 'Au', d_p: 40, n_med: 1.33, lam_min: 400, lam_max: 900, AR_p: 3 } },
    { name: 'Silver nanoparticles (30 nm) in air', note: 'Sharp dipolar plasmon near 360 nm; strong field enhancement.', values: { mode: 'mie', mat_p: 'Ag', d_p: 30, n_med: 1.0, lam_min: 300, lam_max: 700, AR_p: 2.5 } },
    { name: 'Wire + split-ring negative-index metamaterial (X-band)', note: 'Double-negative band between the ring resonance at 10 GHz and its magnetic plasma frequency.', values: { mode: 'nim', fep: 12.8, ge: 0.1, f0m: 10, Fm: 0.3, gm: 0.05, d_slab: 10, f_probe: 11, theta: 30 } },
    { name: 'Pendry split-ring medium (13.5 GHz)', note: 'Reference geometry of the densely stacked ring model: r = 2 mm, c = 1 mm, d = 0.1 mm, l = 2 mm, a = 5 mm.', values: { mode: 'neg_mu', srr_type: 'cylinder', srr_r: 2, srr_a: 5, srr_c: 1, srr_d: 0.1, srr_l: 2, srr_t: 35, metal_sigma: 5.8e7, f_rel: 1.05 } },
    { name: 'Printed copper split rings (Smith et al. 2000, 4.8 GHz)', note: 'Planar double split rings r = 1.5 mm, c = 0.8 mm, d = 0.2 mm on an 8 mm lattice; measured resonance 4.845 GHz. Board assumed FR-4 type, 0.25 mm.', values: { mode: 'neg_mu', srr_type: 'planar', srr_r: 1.5, srr_a: 8, srr_c: 0.8, srr_d: 0.2, srr_l: 8, srr_t: 35, srr_epsr: 4.4, srr_h: 0.25, metal_sigma: 5.8e7, f_rel: 1.02 } },
    { name: 'Thin aluminium wire plasma (8.2 GHz)', note: 'Pendry wire medium: 1 µm wires on a 5 mm lattice give a gigahertz plasma frequency.', values: { mode: 'neg_eps', medium: 'wire', wire_r: 1, wire_a: 5, metal_sigma: 3.7e7, f_rel: 0.6 } },
    { name: 'Ag/TiO₂ hyperbolic multilayer', note: 'Type II hyperbolic and epsilon-near-zero behaviour in the visible.', values: { mode: 'hyperbolic', mat_i: 'Ag', mat_d: 'TiO2', fill: 0.3, period: 20, lam0: 750, lam_min: 350, lam_max: 900 } },
    { name: 'TiO₂ visible metalens, NA 0.45', note: 'Hyperbolic phase profile sampled by 300 nm unit cells.', values: { mode: 'metasurface', ms_kind: 'lens', lam0: 532, ms_p: 300, lens_D: 100, lens_f: 100, n_pillar: 2.4, dlam: 10 } },
    { name: 'ZnO nanoparticle composite film', note: 'ZnO inclusions in a polymer host: effective index of a UV-blocking nanocomposite.', values: { mode: 'effmed', mat_i: 'ZnO', n_host: 1.49, fill: 0.2, lam0: 550, lam_min: 400, lam_max: 900 } },
  ],
  run(inp) {
    const f = { multilayer, neg_eps: negEps, neg_mu: negMu, nim, mie, effmed, hyperbolic: hyper, metasurface }[inp.mode] || multilayer;
    return finish(f(inp));
  },
};
