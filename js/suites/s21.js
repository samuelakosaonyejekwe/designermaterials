// Suite 21 — Fabrication, Synthesis and Process Engineering.
// Coupled 3-D hydrothermal ZnO reactor simulator plus reduced-order deposition, patterning,
// additive-manufacturing, sol–gel/sintering and electrospinning/colloid process models.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { R_GAS: R, FARADAY: F, K_B: KB, N_A: NA, PI } = N;
const f = (v, s = 3) => N.fmt(v, s);
const lsp = (e0, e1, n) => N.linspace(e0, e1, n).map(v => 10 ** v); // n log-spaced samples between 10^e0 and 10^e1
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ser = (key, title, xlabel, ylabel, lines, explain, o = {}) => ({ key, title, xlabel, ylabel, logx: !!o.logx, logy: !!o.logy, lines, explain });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const safeEq = (r, p, arrow) => { try { return chem.equationString(r, p, null, arrow); } catch (e) { return r.join(' + ') + ' → ' + p.join(' + ') + ' (not balanced automatically)'; } };
const MM = s => chem.molarMass(s) * 1e-3; // kg/mol
// electrolytes needed for the ionic equations that the shared table does not list yet (added only if absent)
for (const [k, v] of Object.entries({ 'CH3COONa': [['Na+', 1], ['CH3COO-', 1]], 'CH3COOK': [['K+', 1], ['CH3COO-', 1]], 'CH3COONH4': [['NH4+', 1], ['CH3COO-', 1]], 'K2SO4': [['K+', 2], ['SO4^2-', 1]], '(NH4)2SO4': [['NH4+', 2], ['SO4^2-', 1]] })) if (!chem.DISSOCIATION[k]) chem.DISSOCIATION[k] = v;

// Replace any non-finite number by 0 and say so; thin long series. Last line of defence only.
function finish(r) {
  let bad = 0;
  for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; k.good = null; bad++; }
  for (const s of r.series || []) for (const l of s.lines) {
    if (l.x.length > 400) { const st = Math.ceil(l.x.length / 400), x = [], y = []; for (let i = 0; i < l.x.length; i += st) { x.push(l.x[i]); y.push(l.y[i]); } l.x = x; l.y = y; }
    l.x = Array.from(l.x); l.y = Array.from(l.y, v => { if (!Number.isFinite(v)) { bad++; return 0; } return v; });
  }
  for (const m of r.maps || []) for (let i = 0; i < m.data.length; i++) if (!Number.isFinite(m.data[i])) { m.data[i] = 0; bad++; }
  for (const k in r.publish) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) { delete r.publish[k]; bad++; }
  if (bad) r.warnings.push(`${bad} output value(s) were not finite for this input combination and were set to zero; the inputs are outside the range in which the model is meaningful.`);
  return r;
}

// ================================================================================================
// 1. Coupled 3-D porous-reactor model for hydrothermal ZnO nanorod growth with electric-field heating
// ================================================================================================
const SALTS = {
  'Zn(NO3)2': { an: 'NO3-', za: 1, Da: 1.902e-9, name: 'zinc nitrate', prod: { Na: 'NaNO3', K: 'KNO3', NH4: 'NH4NO3' } },
  'ZnCl2': { an: 'Cl-', za: 1, Da: 2.032e-9, name: 'zinc chloride', prod: { Na: 'NaCl', K: 'KCl', NH4: 'NH4Cl' } },
  'Zn(CH3COO)2': { an: 'CH3COO-', za: 1, Da: 1.089e-9, name: 'zinc acetate', prod: { Na: 'CH3COONa', K: 'CH3COOK', NH4: 'CH3COONH4' } },
  'ZnSO4': { an: 'SO4^2-', za: 2, Da: 1.065e-9, name: 'zinc sulfate', prod: { Na: 'Na2SO4', K: 'K2SO4', NH4: '(NH4)2SO4' } },
};
const BASES = {
  NaOH: { cat: 'Na', Dc: 1.334e-9, weak: false, name: 'sodium hydroxide' }, KOH: { cat: 'K', Dc: 1.957e-9, weak: false, name: 'potassium hydroxide' },
  NH4OH: { cat: 'NH4', Dc: 1.957e-9, weak: true, name: 'aqueous ammonia' }, HMTA: { cat: 'NH4', Dc: 1.957e-9, weak: true, hmta: true, name: 'hexamethylenetetramine (HMTA)' },
};
const D_NH3 = 1.64e-9, D_NH4 = 1.957e-9, D_HMTA = 8e-10, KB_NH3 = 1.8e-2 /* mol/m³ (1.8e-5 mol/L) */, RHO_ZNO = 5606, T_REF = 298.15;
const muW = T => 2.414e-5 * 10 ** (247.8 / (T - 140)); // water viscosity, Pa·s (Vogel-type correlation)
const tSat = p => { const l = Math.log10(Math.max(p, 700) / 133.322); return (l <= 2.8808 ? 1730.63 / (8.07131 - l) - 233.426 : 1810.94 / (8.14019 - l) - 244.485) + 273.15; }; // Antoine, water
const hm = (a, b) => (a + b > 0 ? 2 * a * b / (a + b) : 0);

function mkGrid(nx, ny, nz, Lx, Ly, Lz, ctr) {
  const n = nx * ny * nz, A = () => new Float64Array(n), m = Math.max(nx, ny, nz);
  const g = { nx, ny, nz, n, nxy: nx * ny, dx: Lx / nx, dy: Ly / ny, dz: Lz / nz, aP: A(), aW: A(), aE: A(), aS: A(), aN: A(), aB: A(), aT: A(), b: A(), cp: new Float64Array(m), dp: new Float64Array(m), ctr /* [cell-sweeps used, budget] shared by all fields */ };
  g.A = [[g.aW, g.aE], [g.aS, g.aN], [g.aB, g.aT]];
  return g;
}
// Finite-volume backward-Euler operator: stor·(u−u_old)/dt = ∇·(D∇u) − ∂(w u)/∂z − K u + Q, first-order upwind for w.
// Boundary descriptors [kind, value, h]: 0 no flux, 1 prescribed value, 2 Robin/convective, 3 open (zero gradient, advective flux passes).
function assemble(g, o) {
  const { nx, ny, nz, nxy, dx, dy, dz, aP, aW, aE, aS, aN, aB, aT, b } = g, { D, w, K, Q, old, bot, top, side } = o;
  const sk = side[0], sv = side[1], sh = side[2], bk = bot[0], bv = bot[1], tk = top[0], tv = top[1], sdt = o.stor / o.dt, ix2 = 1 / (dx * dx), iy2 = 1 / (dy * dy), iz2 = 1 / (dz * dz);
  let q = 0;
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, q++) {
    const d = D[q]; let p = sdt + (K ? K[q] : 0), rhs = sdt * old[q] + (Q ? Q[q] : 0), c;
    if (i > 0) { c = hm(d, D[q - 1]) * ix2; aW[q] = c; p += c; } else { aW[q] = 0; if (sk === 2) { c = 1 / (1 / sh + 0.5 * dx / d) / dx; p += c; rhs += c * sv; } }
    if (i < nx - 1) { c = hm(d, D[q + 1]) * ix2; aE[q] = c; p += c; } else { aE[q] = 0; if (sk === 2) { c = 1 / (1 / sh + 0.5 * dx / d) / dx; p += c; rhs += c * sv; } }
    if (j > 0) { c = hm(d, D[q - nx]) * iy2; aS[q] = c; p += c; } else { aS[q] = 0; if (sk === 2) { c = 1 / (1 / sh + 0.5 * dy / d) / dy; p += c; rhs += c * sv; } }
    if (j < ny - 1) { c = hm(d, D[q + nx]) * iy2; aN[q] = c; p += c; } else { aN[q] = 0; if (sk === 2) { c = 1 / (1 / sh + 0.5 * dy / d) / dy; p += c; rhs += c * sv; } }
    if (k > 0) { c = hm(d, D[q - nxy]) * iz2; let a = c; if (w) { const wf = 0.5 * (w[q] + w[q - nxy]); if (wf > 0) a += wf / dz; else p -= wf / dz; } aB[q] = a; p += c; }
    else {
      aB[q] = 0;
      if (bk === 1) { c = 2 * d * iz2; p += c; rhs += c * bv; if (w) { if (w[q] > 0) rhs += w[q] * bv / dz; else p -= w[q] / dz; } }
      else if (bk === 2) { c = 1 / (1 / bot[2] + 0.5 * dz / d) / dz; p += c; rhs += c * bv; }
    }
    if (k < nz - 1) { c = hm(d, D[q + nxy]) * iz2; let a = c; if (w) { const wf = 0.5 * (w[q] + w[q + nxy]); if (wf > 0) p += wf / dz; else a -= wf / dz; } aT[q] = a; p += c; }
    else {
      aT[q] = 0;
      if (tk === 1) { c = 2 * d * iz2; p += c; rhs += c * tv; }
      else if (tk === 2) { c = 1 / (1 / top[2] + 0.5 * dz / d) / dz; p += c; rhs += c * tv; }
      else if (tk === 3 && w) p += w[q] / dz;
    }
    aP[q] = p; b[q] = rhs;
  }
}
// One line-implicit (Thomas) Gauss–Seidel sweep along direction dir (0 x, 1 y, 2 z)
function lineSweep(g, u, dir) {
  const { aP, b, cp, dp, A } = g, off = [1, g.nx, g.nxy], dims = [g.nx, g.ny, g.nz], d1 = (dir + 1) % 3, d2 = (dir + 2) % 3;
  const n0 = dims[dir], n1 = dims[d1], n2 = dims[d2], s0 = off[dir], s1 = off[d1], s2 = off[d2];
  const lo = A[dir][0], up = A[dir][1], m1 = A[d1][0], p1 = A[d1][1], m2 = A[d2][0], p2 = A[d2][1];
  for (let j2 = 0; j2 < n2; j2++) for (let j1 = 0; j1 < n1; j1++) {
    let q = j1 * s1 + j2 * s2;
    for (let m = 0; m < n0; m++, q += s0) {
      let rhs = b[q], c;
      if ((c = m1[q]) !== 0) rhs += c * u[q - s1]; if ((c = p1[q]) !== 0) rhs += c * u[q + s1];
      if ((c = m2[q]) !== 0) rhs += c * u[q - s2]; if ((c = p2[q]) !== 0) rhs += c * u[q + s2];
      const l = lo[q], den = m ? aP[q] - l * cp[m - 1] : aP[q];
      cp[m] = up[q] / den; dp[m] = (m ? rhs + l * dp[m - 1] : rhs) / den;
    }
    q -= s0; u[q] = dp[n0 - 1];
    for (let m = n0 - 2; m >= 0; m--) { q -= s0; u[q] = dp[m] + cp[m] * u[q + s0]; }
  }
}
function linResid(g, u) {
  const { n, nx, nxy, aP, aW, aE, aS, aN, aB, aT, b } = g; let rm = 0, sm = 0;
  for (let q = 0; q < n; q++) {
    let r = b[q] - aP[q] * u[q], c;
    if ((c = aW[q]) !== 0) r += c * u[q - 1]; if ((c = aE[q]) !== 0) r += c * u[q + 1];
    if ((c = aS[q]) !== 0) r += c * u[q - nx]; if ((c = aN[q]) !== 0) r += c * u[q + nx];
    if ((c = aB[q]) !== 0) r += c * u[q - nxy]; if ((c = aT[q]) !== 0) r += c * u[q + nxy];
    r = Math.abs(r); if (r > rm) rm = r; const sc = Math.abs(aP[q] * u[q]) + Math.abs(b[q]); if (sc > sm) sm = sc;
  }
  return sm > 0 ? rm / sm : 0;
}
// Alternating-direction line relaxation to a relative residual; returns {res, sweeps}
function relax(g, u, tol, maxSw) {
  let res = linResid(g, u), sw = 0;
  const c = g.ctr;
  while (res >= tol && sw < maxSw && (sw === 0 || c[0] < c[1])) { if (sw > 0) { lineSweep(g, u, 0); lineSweep(g, u, 1); c[0] += 2 * g.n; } lineSweep(g, u, 2); c[0] += 2 * g.n; sw++; res = linResid(g, u); }
  return res;
}

// Conjugate gradients with an exact z-line (block-Jacobi) preconditioner for the symmetric operators (temperature, pressure)
function pcgSym(g, u, tol, maxIt) {
  const { n, nx, nxy, nz, aP, aW, aE, aS, aN, aB, aT, b, cp, dp } = g, r = new Float64Array(n), z = new Float64Array(n), pv = new Float64Array(n), Ap = new Float64Array(n);
  const mv = (x, y) => { for (let q = 0; q < n; q++) { let v = aP[q] * x[q], c; if ((c = aW[q]) !== 0) v -= c * x[q - 1]; if ((c = aE[q]) !== 0) v -= c * x[q + 1]; if ((c = aS[q]) !== 0) v -= c * x[q - nx]; if ((c = aN[q]) !== 0) v -= c * x[q + nx]; if ((c = aB[q]) !== 0) v -= c * x[q - nxy]; if ((c = aT[q]) !== 0) v -= c * x[q + nxy]; y[q] = v; } };
  const pre = () => { for (let c0 = 0; c0 < nxy; c0++) { let q = c0; for (let m = 0; m < nz; m++, q += nxy) { const l = aB[q], den = m ? aP[q] - l * cp[m - 1] : aP[q]; cp[m] = aT[q] / den; dp[m] = (m ? r[q] + l * dp[m - 1] : r[q]) / den; } q -= nxy; z[q] = dp[nz - 1]; for (let m = nz - 2; m >= 0; m--) { q -= nxy; z[q] = dp[m] + cp[m] * z[q + nxy]; } } };
  let sc = 0; for (let q = 0; q < n; q++) { const v = Math.abs(aP[q] * u[q]) + Math.abs(b[q]); if (v > sc) sc = v; } if (!(sc > 0)) return 0;
  mv(u, Ap); let rm = 0; for (let q = 0; q < n; q++) { r[q] = b[q] - Ap[q]; const v = Math.abs(r[q]); if (v > rm) rm = v; }
  if (rm / sc < tol) return rm / sc;
  pre(); pv.set(z); let rz = 0; for (let q = 0; q < n; q++) rz += r[q] * z[q];
  const ctr = g.ctr;
  for (let it = 0; it < maxIt && rm / sc >= tol && (it < 2 || ctr[0] < ctr[1]); it++) {
    mv(pv, Ap); let pAp = 0; for (let q = 0; q < n; q++) pAp += pv[q] * Ap[q]; if (!(pAp > 0)) break;
    const al = rz / pAp; rm = 0; for (let q = 0; q < n; q++) { u[q] += al * pv[q]; r[q] -= al * Ap[q]; const v = Math.abs(r[q]); if (v > rm) rm = v; }
    pre(); let rz2 = 0; for (let q = 0; q < n; q++) rz2 += r[q] * z[q]; const be = rz2 / rz; rz = rz2; for (let q = 0; q < n; q++) pv[q] = z[q] + be * pv[q];
    ctr[0] += 3 * n;
  }
  return rm / sc;
}

// Coupled transient solve. P: SI parameter object. Returns fields, histories, inventories and solver statistics.
function znoSim(P, nx, ny, nz, dt, nSteps, full, budget) {
  const ctr = new Float64Array([0, budget]), mk = () => mkGrid(nx, ny, nz, P.Lx, P.Ly, P.Lz, ctr), g = mk(), gT = mk(), gB = mk(), gH = P.hmta ? mk() : null, gP = mk(), n = g.n, nxy = nx * ny, { dx, dy, dz } = g, Vc = dx * dy * dz, Ac = dx * dy;
  const A = v => { const a = new Float64Array(n); if (v) a.fill(v); return a; };
  const cZ = A(P.cZ0), cB = A(P.cB0), cH = A(P.cH0), s = A(), nn = A(), T = A(P.T0), p = A(P.p0), Lr = A(), dr = A();
  const cZo = A(), cBo = A(), cHo = A(), so = A(), no = A(), To = A(), po = A(), sNew = A(), nNew = A();
  const Ns = A(), E0c = new Float64Array(nxy), sig = A(), xf = A(), DZa = A(), DBa = A(), DHa = A(), wZ = A(), wB = A(), KZ = A(), KBs = A(), QZ = A(), QB = A(), KH = A(), qJ = A(), Qh = A(), Rd = A(), kap = A(), Qp = A(), prevZ = A(), prevB = A(), prevT = A();
  // seeded, grid-independent nucleation-site density field (sum of random long-wave cosines, log-normal amplitude)
  const rnd = N.rng(P.seed), modes = Array.from({ length: 6 }, () => [1 + rnd.int(4), 1 + rnd.int(4), rnd.int(3), 2 * PI * rnd()]);
  for (let k = 0, q = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, q++) {
    const X = (i + 0.5) / nx, Y = (j + 0.5) / ny, Z = (k + 0.5) / nz; let gs = 0;
    for (const m of modes) gs += Math.cos(2 * PI * (m[0] * X + m[1] * Y) + PI * m[2] * Z + m[3]);
    gs *= Math.sqrt(2 / modes.length);
    Ns[q] = P.as * P.Ns0 * Math.exp(P.cv * gs - 0.5 * P.cv * P.cv);
    if (k === 0) E0c[i + nx * j] = 1 - P.fieldNU * 0.5 * ((2 * X - 1) ** 2 + (2 * Y - 1) ** 2);
  }
  const tor = P.phi ** 1.5, v0 = PI / 4 * P.d0 ** 3, mu0 = muW(T_REF);
  const spBot = P.bc === 'batch' ? 0 : 1, spTop = P.bc === 'reference' ? 3 : 0;
  // temperature operator is constant in time: assemble once, keep the boundary part of the right-hand side
  const kArr = A(P.kth), zero = A(), bT0 = A(), rcp = P.rhocp, robT = [2, P.Tinf, P.h];
  assemble(gT, { stor: rcp, dt, D: kArr, old: zero, bot: robT, top: robT, side: robT }); bT0.set(gT.b);
  const Uz = 1 / (1 / P.h + 0.5 * dz / P.kth), Ux = 1 / (1 / P.h + 0.5 * dx / P.kth), Uy = 1 / (1 / P.h + 0.5 * dy / P.kth);
  const d0 = P.d0, hb = 0.5 * P.beta;
  const rod = (v, out, dg) => { // per-rod volume → mean diameter and length (anisotropic growth closure, Newton on the cubic)
    if (!(v > v0)) { const d = d0 * Math.cbrt((v > 0 ? v : 0) / v0); out[0] = d; out[1] = d; return; }
    const ub = d0 + Math.cbrt(4 * (v - v0) / (PI * hb)); let d = dg > d0 && dg < ub ? dg : ub;
    for (let it = 0; it < 30; it++) { const L = d0 + hb * (d - d0), dd = (PI / 4 * d * d * L - v) / (PI / 4 * (2 * d * L + hb * d * d)); d -= dd; if (d > ub) d = ub; if (Math.abs(dd) < 1e-10 * d) break; }
    out[0] = d; out[1] = d0 + hb * (d - d0);
  };
  const dl = [0, 0], fTa = A();
  const { phi, DZ, DOH, weak, hmta, Apre, dH, Ksp0, fth, w0, Vm, kappa, kh90, Lz } = P, specC = P.spec, cF = tor * F * F / R, EaR = P.Ea / R, dHR = dH / R, EaHR = P.EaH / R, cG = 16 * PI * P.gam ** 3 * P.vm * P.vm / 3, iT0 = 1 / T_REF, arr0 = Math.exp(EaR * iT0), fR = F / R, DBm = weak ? D_NH3 : DOH;
  // temperature-dependent property tables (linear interpolation, 0.05 K spacing; relative interpolation error < 1e-5)
  const NT = 8192, T_LO = 255, T_HI = 650, iDT = (NT - 1) / (T_HI - T_LO), tF = new Float64Array(NT), tMu = new Float64Array(NT), tArr = new Float64Array(NT), tKsp = new Float64Array(NT), tKH = new Float64Array(NT);
  for (let i = 0; i < NT; i++) { const Tc = T_LO + i / iDT, mu = muW(Tc); tF[i] = Tc * iT0 * mu0 / mu; tMu[i] = 1 / mu; tArr[i] = Math.exp(-EaR / Tc); tKsp[i] = Ksp0 * Math.exp(dHR * (1 / Tc - iT0)); tKH[i] = phi * kh90 * Math.exp(-EaHR * (1 / Tc - 1 / 363.15)); }
  const bcZ = [spBot, 0], bcB = [spBot, 0], bcH = [spBot, 0], bcTop = [spTop], bc0 = [0], bcP = [1, P.p0];
  const oH = { stor: phi, dt, D: DHa, K: KH, old: cHo, bot: bcH, top: bc0, side: bc0 }, oZ = { stor: phi, dt, D: DZa, w: wZ, K: KZ, Q: QZ, old: cZo, bot: bcZ, top: bcTop, side: bc0 }, oB = { stor: phi, dt, D: DBa, w: wB, K: KBs, Q: QB, old: cBo, bot: bcB, top: bcTop, side: bc0 }, oP = { stor: P.Sp, dt, D: kap, Q: Qp, old: po, bot: bc0, top: bcP, side: bc0 };
  // histories and accounting
  const H = { t: [], cZ: [], cB: [], s: [], Tm: [], Tx: [], L: [], d: [], ns: [], G: [], px: [], sx: [], ux: [], E: [], dose: [], PJ: [], PR: [], it: [], res: [], J: [] };
  const frames = [], stride = Math.max(1, Math.ceil(nSteps / 200)), fStride = Math.max(1, Math.ceil(nSteps / 30));
  let Efield = P.E0, dose = 1, integT = 0, integG = 0, Tprev = P.T0, Lprev = 0, Gax = 0;
  let itTot = 0, itMax = 0, resMax = 0, linMax = 0, swTot = 0, notConv = 0, Jpeak = 0, tSettle = 0, capped = 0, damped = 0;
  const invZ0 = P.phi * P.cZ0 * n * Vc, invB0 = P.phi * P.cB0 * n * Vc, invH0 = P.phi * P.cH0 * n * Vc;
  const stats = () => { // n-weighted rod statistics over the bed
    let sn = 0, sL = 0, sd = 0, sL2 = 0, Lx_ = 0, dx_ = 0, sT = 0, Tx_ = -1e9, sZ = 0, sB = 0, ss = 0, px_ = -1e99;
    for (let q = 0; q < n; q++) { const w = nn[q]; sn += w; sL += w * Lr[q]; sd += w * dr[q]; sL2 += w * Lr[q] * Lr[q]; if (w > 0 && Lr[q] > Lx_) Lx_ = Lr[q]; if (w > 0 && dr[q] > dx_) dx_ = dr[q]; sT += T[q]; if (T[q] > Tx_) Tx_ = T[q]; sZ += cZ[q]; sB += xf[q]; ss += s[q]; if (p[q] > px_) px_ = p[q]; }
    const Lm = sn > 0 ? sL / sn : 0, dm = sn > 0 ? sd / sn : 0;
    return { Lm, dm, Lmax: Lx_, dmax: dx_, cvL: Lm > 0 ? Math.sqrt(Math.max(0, sL2 / sn - Lm * Lm)) / Lm : 0, ns: sn / n / P.as, Tm: sT / n, Tx: Tx_, cZ: sZ / n, cB: sB / n, s: ss / n, px: px_ };
  };
  const mech = () => { // oedometric poro-thermo-elastic column solution (σzz = 0, εxx = εyy = 0)
    let sx = 0, ux = 0; const sxx = full ? A() : null, uz = full ? A() : null;
    for (let c = 0; c < nxy; c++) { let u = 0; for (let k = 0; k < nz; k++) { const q = c + k * nxy, th = 3 * P.Kd * (P.aT * (T[q] - P.T0) + P.etaG * s[q] * P.Vm / 3), dp_ = p[q] - P.p0, ez = (P.biot * dp_ + th) / P.Mc, sg = P.lam * ez - P.biot * dp_ - th; u += ez * dz; if (Math.abs(sg) > Math.abs(sx)) sx = sg; if (Math.abs(u) > Math.abs(ux)) ux = u; if (full) { sxx[q] = sg; uz[q] = u; } } }
    return { sx, ux, sxx, uz };
  };
  for (let q = 0; q < n; q++) xf[q] = P.weak ? 0.5 * (-KB_NH3 + Math.sqrt(KB_NH3 * KB_NH3 + 4 * KB_NH3 * cB[q])) : cB[q];
  const record = (t, it, res, PJ, PR, Jm) => { const st = stats(), m = mech(); H.t.push(t); H.cZ.push(st.cZ); H.cB.push(st.cB); H.s.push(st.s); H.Tm.push(st.Tm); H.Tx.push(st.Tx); H.L.push(st.Lm * 1e9); H.d.push(st.dm * 1e9); H.ns.push(st.ns); H.G.push(Gax * 6e10); H.px.push(st.px); H.sx.push(Math.abs(m.sx)); H.ux.push(Math.abs(m.ux)); H.E.push(Efield); H.dose.push(dose); H.PJ.push(PJ); H.PR.push(PR); H.it.push(it); H.res.push(Math.max(res, 1e-16)); H.J.push(Jm); };
  const frame = t => { const L = new Float32Array(nxy), d = new Float32Array(nxy), a = new Float32Array(nxy); for (let c = 0; c < nxy; c++) { let sn = 0, sL = 0, sd = 0, sN = 0; for (let k = 0; k < nz; k++) { const q = c + k * nxy; sn += nn[q]; sL += nn[q] * Lr[q]; sd += nn[q] * dr[q]; sN += Ns[q]; } L[c] = sn > 0 ? sL / sn : 0; d[c] = sn > 0 ? sd / sn : 0; a[c] = sN > 0 ? sn / sN : 0; } frames.push({ t, L, d, a }); };
  record(0, 0, 0, 0, 0, 0); if (full) frame(0);

  const acc = new Float64Array(16); // 0 inZ, 1 inB, 2 inH, 3 consB, 4 srcB, 5 decH, 6 EJ, 7 ER, 8 Eloss, 9 PJ, 10 PR, 11 Jm, 12 σ mean, 13 L mean
  // ---- properties, electrical closure (current continuity per column), kinetics ----
  const props = Efield => {
  for (let c = 0; c < nxy; c++) {
    let Rs = 0;
    for (let k = 0; k < nz; k++) {
      const q = c + k * nxy, Tq = T[q], Tc = Tq < 255 ? 255 : Tq > 650 ? 650 : Tq, iT = 1 / Tc, ut = (Tc - T_LO) * iDT, i0 = ut >= NT - 1 ? NT - 2 : ut | 0, wt = ut - i0, fT = tF[i0] + wt * (tF[i0 + 1] - tF[i0]), B = cB[q] > 0 ? cB[q] : 0, z = cZ[q] > 0 ? cZ[q] : 0, sq = s[q];
      let x = B, sg = 4 * DZ * z + specC;
      if (weak) { const a = 2 * sq / phi + KB_NH3; x = 0.5 * (-a + Math.sqrt(a * a + 4 * KB_NH3 * B)); sg += DOH * x + D_NH4 * (x + 2 * sq / phi); } else sg += DOH * x;
      xf[q] = x; sg = cF * iT * fT * sg; if (sg < 1e-9) sg = 1e-9; sig[q] = sg; Rs += dz / sg; fTa[q] = fT;
      DZa[q] = tor * DZ * fT; DBa[q] = tor * DBm * fT; if (hmta) { DHa[q] = tor * D_HMTA * fT; KH[q] = tKH[i0] + wt * (tKH[i0 + 1] - tKH[i0]); }
      kap[q] = kappa * (tMu[i0] + wt * (tMu[i0 + 1] - tMu[i0]));
      // kinetics
      const arr = tArr[i0] + wt * (tArr[i0 + 1] - tArr[i0]), kr = Apre * arr, Ksp = tKsp[i0] + wt * (tKsp[i0 + 1] - tKsp[i0]), Q = z * x * x, nOld = no[q];
      let nNw = nOld, Rn = 0, Rg = 0, rd = 0;
      if (Q > Ksp && nOld < 0.9999999 * Ns[q]) {
        const kT = KB * Tc * Math.log(Q / Ksp), om = w0 * arr * arr0 * Math.exp(-fth * cG / (kT * kT * KB * Tc));
        nNw = (nOld + dt * om * Ns[q]) / (1 + dt * om); Rn = (nNw - nOld) * v0 / (Vm * dt);
      }
      nNew[q] = nNw;
      if (nNw > 0) {
        rod(sq * Vm / nNw, dl, dr[q]); const d = dl[0], L = dl[1]; let Ar = nNw * PI * (d * L + 0.25 * d * d);
        if (Q > Ksp) {
          const g0 = kr * (1 - Ksp / Q), zx = z * x * x;
          if (d > 0.3 * d0) { // one Newton step on the local area–product feedback: A(s) ≈ A* + (dA/ds)(s − s*), s = s_old + Δt·r
            const dAds = d >= d0 ? Vm * (L + hb * d + 0.5 * d) / (0.25 * (2 * d * L + hb * d * d)) : 10 / 3 * Vm / d, c1 = dt * g0 * zx * dAds;
            if (c1 < 0.8) { const sN = (so[q] + dt * (g0 * zx * (Ar - dAds * sq) + Rn)) / (1 - c1), An = Ar + dAds * (sN - sq); Ar = An < 0.2 * Ar ? 0.2 * Ar : An > 5 * Ar ? 5 * Ar : An; }
          }
          Rg = g0 * Ar;
        } else rd = Math.min(kr * Ar * (Ksp - Q), so[q] / dt);
      }
      // Patankar (implicit, positivity-preserving) sinks: R = Rg·cZ·x² + Rn
      KZ[q] = Rg * x * x + Rn / (z > 1e-12 ? z : 1e-12); QZ[q] = rd; Rd[q] = rd;
      KBs[q] = B > 1e-12 ? 2 * (Rg * z * x * x + Rn) / B : 0; QB[q] = 2 * rd + (hmta ? 4 * KH[q] * (cH[q] > 0 ? cH[q] : 0) : 0);
    }
    const jc = Efield * E0c[c] * Lz / Rs;
    for (let k = 0; k < nz; k++) { const q = c + k * nxy, E = jc / sig[q], Tq = T[q], Tc = Tq < 255 ? 255 : Tq > 650 ? 650 : Tq, mob = fR / Tc * E; qJ[q] = jc * E; wZ[q] = 2 * mob * DZa[q]; wB[q] = weak ? -mob * tor * DOH * fTa[q] * (cB[q] > 1e-12 ? xf[q] / cB[q] : 0) : -mob * DBa[q]; }
  }
  };
  const bT = gT.b, sdtT = rcp / dt, cSc = 1 / P.cScale, cScB = 1 / P.cScaleB, iTs = 1 / P.T0, betaP = P.betaP, Tinf = P.Tinf;
  const product = () => { for (let q = 0; q < n; q++) { if (cZ[q] < 0) cZ[q] = 0; if (cB[q] < 0) cB[q] = 0; const Rnet = KZ[q] * cZ[q] - Rd[q], sv = so[q] + dt * Rnet; s[q] = sv > 0 ? sv : 0; nn[q] = nNew[q]; bT[q] = bT0[q] + sdtT * To[q] + qJ[q] - dH * Rnet; } };
  const update = () => { let ch = 0; for (let q = 0; q < n; q++) { const a = Math.abs(cZ[q] - prevZ[q]) * cSc, b2 = Math.abs(cB[q] - prevB[q]) * cScB, c2 = Math.abs(T[q] - prevT[q]) * iTs; if (a > ch) ch = a; if (b2 > ch) ch = b2; if (c2 > ch) ch = c2; } return ch; };
  const predict = (u, uo, uoo) => { for (let q = 0; q < n; q++) { const v = 2 * u[q] - uoo[q]; uoo[q] = uo[q] = u[q]; u[q] = v > 0 ? v : 0; } };
  const account = (cZb, cBb, cHb) => { // rods, inventories, boundary fluxes, energy
    let PJ = 0, PR = 0, loss = 0, Jm = 0, sg = 0, cons = 0, dec = 0, srcH = 0, inZ = 0, inB = 0, inH = 0, sn = 0, sL = 0;
    for (let q = 0; q < n; q++) {
      if (nn[q] > 0) { rod(s[q] * Vm / nn[q], dl, dr[q]); dr[q] = dl[0]; Lr[q] = dl[1]; } else { dr[q] = 0; Lr[q] = 0; }
      const Rnet = KZ[q] * cZ[q] - Rd[q]; PJ += qJ[q]; PR -= dH * Rnet; cons += KBs[q] * cB[q] - 2 * Rd[q]; Jm += nn[q] - no[q]; sg += sig[q]; sn += nn[q]; sL += nn[q] * Lr[q];
      if (hmta) { dec += KH[q] * cH[q]; srcH += QH[q]; }
    }
    for (let c = 0; c < nxy; c++) {
      const q0 = c, q1 = c + (nz - 1) * nxy;
      if (spBot === 1) { inZ += 2 * DZa[q0] / dz * (cZb - cZ[q0]) + (wZ[q0] > 0 ? wZ[q0] * cZb : wZ[q0] * cZ[q0]); inB += 2 * DBa[q0] / dz * (cBb - cB[q0]) + (wB[q0] > 0 ? wB[q0] * cBb : wB[q0] * cB[q0]); if (hmta) inH += 2 * DHa[q0] / dz * (cHb - cH[q0]); }
      if (spTop === 3) { inZ -= wZ[q1] * cZ[q1]; inB -= wB[q1] * cB[q1]; }
      loss += Ac * Uz * (T[q0] - Tinf + T[q1] - Tinf);
    }
    for (let k = 0; k < nz; k++) { for (let j = 0; j < ny; j++) { const q = nx * (j + ny * k); loss += dy * dz * Ux * (T[q] - Tinf + T[q + nx - 1] - Tinf); } for (let i = 0; i < nx; i++) { const q = i + nxy * k; loss += dx * dz * Uy * (T[q] - Tinf + T[q + nx * (ny - 1)] - Tinf); } }
    acc[0] += inZ * Ac * dt; acc[1] += inB * Ac * dt; acc[2] += inH * Ac * dt; acc[3] += cons * Vc * dt; acc[4] += srcH * Vc * dt; acc[5] += dec * Vc * dt;
    acc[9] = PJ * Vc; acc[10] = PR * Vc; acc[6] += acc[9] * dt; acc[7] += acc[10] * dt; acc[8] += loss * dt; acc[11] = Jm / n / dt; acc[12] = sg / n; acc[13] = sn > 0 ? sL / sn : 0;
  };
  const polish = () => { const l1 = relax(g, cZ, P.linTol, P.maxSw), l2 = relax(gB, cB, P.linTol, P.maxSw); product(); const l3 = pcgSym(gT, T, 1e-11, 80); return Math.max(l1, l2, l3 * 100); };
  const latS = tor * Math.max(DZ, DBm) * 4 * dt / Math.min(dx, dy) ** 2 > 0.02, latT = P.kth / rcp * dt / Math.min(dx, dy) ** 2 > 0.02, QH = A();
  const cZoo = A(P.cZ0), cBoo = A(P.cB0), Too = A(P.T0), soo = A();
  for (let step = 1; step <= nSteps; step++) {
    const t = step * dt;
    // ---- sampled-data feedback control (acts on the state at the start of the step) ----
    if (P.ctrlT) { let Tm = 0; for (let q = 0; q < n; q++) Tm += T[q]; Tm /= n; const e = P.Tset - Tm, de = -(Tm - Tprev) / dt, u = P.E0 + P.Kp * e + P.Ki * (integT + e * dt) + P.Kder * de; if ((u > 0 && u < P.Emax) || (u >= P.Emax && e < 0) || (u <= 0 && e > 0)) integT += e * dt; Efield = N.clamp(P.E0 + P.Kp * e + P.Ki * integT + P.Kder * de, 0, P.Emax); Tprev = Tm; if (Math.abs(e) > 0.02 * Math.max(1, Math.abs(P.Tset - P.T0))) tSettle = t; }
    if (P.ctrlG && step > 1) { const e = (P.Gset - Gax) / P.Gset, u = 1 + P.Kg * e + P.Kig * (integG + e * dt); if ((u > 0 && u < P.doseMax) || (u >= P.doseMax && e < 0) || (u <= 0 && e > 0)) integG += e * dt; dose = N.clamp(1 + P.Kg * e + P.Kig * integG, 0, P.doseMax); }
    // old state, and a linear-extrapolation predictor as the first Picard iterate
    predict(cZ, cZo, cZoo); predict(cB, cBo, cBoo); predict(T, To, Too); predict(s, so, soo); cHo.set(cH); no.set(nn); po.set(p);
    let it = 0, change = 1, lin = 0, chPrev = 1e9;
    const cZb = P.cZb * dose, cBb = P.cBb * dose, cHb = P.cHb * dose; bcZ[1] = cZb; bcB[1] = cBb; bcH[1] = cHb;
    const itCap = ctr[0] > ctr[1] ? 2 : P.maxIt; if (itCap < P.maxIt) capped++;
    while (it < itCap) {
      it++; prevZ.set(cZ); prevB.set(cB); prevT.set(T);
      props(Efield);
      // weak lateral coupling: one exact z-line solve per field and iteration; otherwise iterate the inner solves. Residuals are verified before the step is accepted.
      if (hmta) { assemble(gH, oH); lineSweep(gH, cH, 2); ctr[0] += n; relax(gH, cH, P.linTol, P.maxSw); for (let q = 0; q < n; q++) { const v = 4 * KH[q] * (cH[q] > 0 ? cH[q] : 0); QH[q] = v; QB[q] = 2 * Rd[q] + v; } }
      assemble(g, oZ); lineSweep(g, cZ, 2); if (latS) relax(g, cZ, 1e-7, P.maxSw);
      for (let q = 0; q < n; q++) { const B = prevB[q]; KBs[q] = B > 1e-12 ? 2 * KZ[q] * (cZ[q] > 0 ? cZ[q] : 0) / B : 0; } // Gauss–Seidel: hydroxide sink = 2 × the zinc sink just solved
      assemble(gB, oB); lineSweep(gB, cB, 2); if (latS) relax(gB, cB, 1e-7, P.maxSw);
      product();
      if (latT) pcgSym(gT, T, 1e-8, 40); else lineSweep(gT, T, 2); ctr[0] += 3 * n;
      change = update();
      if (it >= 3 && change > 0.6 * chPrev) { for (let q = 0; q < n; q++) { cZ[q] = 0.5 * (cZ[q] + prevZ[q]); cB[q] = 0.5 * (cB[q] + prevB[q]); T[q] = 0.5 * (T[q] + prevT[q]); } damped++; } // stalled or oscillating: under-relax
      chPrev = change;
      if (change < P.tol) { lin = polish(); change = update(); if (change < P.tol) break; }
    }
    if (change >= P.tol) { lin = polish(); notConv++; }
    swTot = ctr[0] / n;
    itTot += it; if (it > itMax) itMax = it; if (change > resMax) resMax = change; if (lin > linMax) linMax = lin;
    // ---- pore pressure (one-way coupled to the converged temperature increment) ----
    for (let q = 0; q < n; q++) Qp[q] = betaP * (T[q] - To[q]) / dt;
    assemble(gP, oP); const rp = 100 * pcgSym(gP, p, 1e-11, 80); if (rp > linMax) linMax = rp;
    account(cZb, cBb, cHb);
    const Lm = acc[13]; Gax = (Lm - Lprev) / dt; Lprev = Lm; if (acc[11] > Jpeak) Jpeak = acc[11];
    if (step % stride === 0 || step === nSteps) record(t, it, change, acc[9], acc[10], acc[11]);
    if (full && (step % fStride === 0 || step === nSteps)) frame(t);
  }
  const inZ = acc[0], inB = acc[1], inH = acc[2], consB = acc[3], srcB = acc[4], decH = acc[5], EJ = acc[6], ER = acc[7], Eloss = acc[8], sigMean = acc[12];
  // ---- conservation residuals ----
  let iZ = 0, iB = 0, iH = 0, iS = 0, dU = 0, minC = Infinity;
  for (let q = 0; q < n; q++) { iZ += cZ[q]; iB += cB[q]; iH += cH[q]; iS += s[q]; dU += T[q] - P.T0; if (cZ[q] < minC) minC = cZ[q]; if (cB[q] < minC) minC = cB[q]; }
  iZ *= P.phi * Vc; iB *= P.phi * Vc; iH *= P.phi * Vc; iS *= Vc; dU *= rcp * Vc;
  const scZ = Math.max(invZ0, Math.abs(inZ), iS, 1e-300), scB = Math.max(invB0, Math.abs(inB), consB, srcB, 1e-300), scE = Math.max(Math.abs(EJ), Math.abs(ER), Math.abs(Eloss), Math.abs(dU), 1e-300);
  const bal = { zn: Math.abs(iZ + iS - invZ0 - inZ) / scZ, oh: Math.abs(iB - invB0 - inB + consB - srcB) / scB, stoich: iS > 0 ? (Math.abs(consB - 2 * iS) + Math.abs(srcB - 4 * decH)) / (2 * iS) : 0, hmta: P.hmta ? Math.abs(iH - invH0 - inH + decH) / Math.max(invH0, Math.abs(inH), 1e-300) : 0, energy: Math.abs(dU - EJ - ER + Eloss) / scE, minC };
  return { g, H, frames, st: stats(), mech: mech(), fields: { cZ, cB, xf, s, nn, T, p, Lr, dr, sig, qJ }, bal, inv: { invZ0, invB0, invH0, iZ, iB, iH, iS, inZ, inB, inH, consB, srcB, decH, dU, EJ, ER, Eloss }, num: { itTot, itMax, resMax, linMax, swTot, notConv, nSteps, dt, capped, damped }, Jpeak, Efield, dose, tSettle, sigMean, Gax };
}

const ZNO_KEYS = ['salt', 'base', 'cZn', 'cOH', 'c_hmta', 'k_hmta90', 'Ea_hmta', 'T0', 'Tinf', 'p0', 'Lx', 'Ly', 'Lz', 'phi', 'dgrain', 'kappa', 'DZn', 'DOH', 'kth', 'cp', 'rho', 'h', 'Ea', 'A_pre', 'dH_rxn', 'pKsp', 'Ns', 'site_cv', 'gamma', 'theta', 'w_nuc', 'd_seed', 'beta_aniso', 'E', 'field_nu', 'spectators', 'Emod', 'nu', 'biot', 'cte', 'eta_g', 'bc', 'nx', 'ny', 'nz', 'dt', 't_total', 'max_it', 'ctrl', 'T_set', 'Kp', 'Ki', 'Kd', 'E_max', 'G_set', 'Kg', 'dose_max', 'rod_cv', 'tilt_sd', 'seed'];
const ZCACHE = new Map();
const WORK = 2.5e5, SWEEPS = 4.5e6; // cell-sweeps of linear relaxation allowed; beyond it inner and coupling iterations are cut to the minimum // cell-steps allowed for the main run (keeps the worst case within the time budget)

function runZnO(inp) {
  const key = JSON.stringify(ZNO_KEYS.map(k => inp[k]));
  if (ZCACHE.has(key)) return ZCACHE.get(key)();
  const build = () => buildZnO(inp);
  const r0 = build(); // results hold typed arrays; rebuild cheaply from the cached simulation on repeat calls
  if (ZCACHE.size > 3) ZCACHE.delete(ZCACHE.keys().next().value);
  ZCACHE.set(key, () => r0.rebuild());
  return r0.result;
}

function buildZnO(inp) {
  const salt = SALTS[inp.salt] || SALTS['Zn(NO3)2'], base = BASES[inp.base] || BASES.NaOH, saltF = SALTS[inp.salt] ? inp.salt : 'Zn(NO3)2', baseF = BASES[inp.base] ? inp.base : 'NaOH';
  const warnings = [];
  const phi = N.clamp(inp.phi, 0.02, 0.98), nu = N.clamp(inp.nu, -0.9, 0.49), Em = inp.Emod * 1e9, Kd = Em / (3 * (1 - 2 * nu)), Gs = Em / (2 * (1 + nu)), Mc = Kd + 4 * Gs / 3, biot = N.clamp(inp.biot, phi, 1);
  if (biot !== inp.biot) warnings.push(`Biot coefficient raised to the porosity (${f(phi)}), its thermodynamic lower bound; the entered value ${f(inp.biot)} is not admissible.`);
  const Kf = 2.2e9, Ks = biot < 0.999 ? Kd / (1 - biot) : 1e15, invM = phi / Kf + (biot - phi) / Ks, betaF = 2.1e-4, betaS = 3 * inp.cte;
  const Vm = MM('ZnO') / RHO_ZNO, th = inp.theta * PI / 180;
  const cZ0 = inp.cZn, cB0 = base.hmta ? 0 : inp.cOH, cH0 = base.hmta ? inp.c_hmta : 0;
  const specOn = !!inp.spectators, spec = specOn ? salt.za * salt.Da * 2 * cZ0 + (base.weak ? 0 : base.Dc * cB0) : 0;
  const mode = inp.ctrl, P = {
    Lx: inp.Lx, Ly: inp.Ly, Lz: inp.Lz, phi, as: 6 * (1 - phi) / (inp.dgrain * 1e-3), Ns0: inp.Ns, cv: inp.site_cv / 100, seed: Math.round(inp.seed) || 1,
    cZ0, cB0, cH0, cZb: cZ0, cBb: cB0, cHb: cH0, T0: inp.T0, Tinf: inp.Tinf, p0: inp.p0, DZ: inp.DZn, DOH: inp.DOH, weak: !!base.weak, hmta: !!base.hmta, kh90: inp.k_hmta90, EaH: inp.Ea_hmta * 1e3, spec,
    kth: inp.kth, rhocp: inp.rho * inp.cp, h: inp.h, Ea: inp.Ea * 1e3, Apre: inp.A_pre, dH: inp.dH_rxn * 1e3, Ksp0: 10 ** (-inp.pKsp) * 1e9, gam: inp.gamma, fth: (2 + Math.cos(th)) * (1 - Math.cos(th)) ** 2 / 4, w0: inp.w_nuc,
    vm: Vm / NA, Vm, d0: inp.d_seed * 1e-9, beta: inp.beta_aniso, E0: inp.E, fieldNU: inp.field_nu / 100, kappa: inp.kappa, Kd, Mc, lam: Kd - 2 * Gs / 3, biot, aT: inp.cte, etaG: inp.eta_g,
    Sp: invM + biot * biot / Mc, betaP: phi * betaF + (biot - phi) * betaS - biot * Kd * betaS / Mc, bc: inp.bc, maxIt: Math.round(N.clamp(inp.max_it, 2, 12)), tol: 1e-4, linTol: 1e-9, maxSw: 6,
    cScale: Math.max(cZ0, 1e-6), cScaleB: Math.max(cB0, 4 * cH0, 1e-6),
    ctrlT: mode === 'temp' || mode === 'both', ctrlG: mode === 'growth' || mode === 'both', Tset: inp.T_set, Kp: inp.Kp, Ki: inp.Ki, Kder: inp.Kd, Emax: inp.E_max, Gset: Math.max(inp.G_set, 1e-6) / 6e10, Kg: inp.Kg, Kig: inp.Kg / 300, doseMax: inp.dose_max,
  };
  const nx = Math.round(N.clamp(inp.nx, 4, 40)), ny = Math.round(N.clamp(inp.ny, 4, 40)), nz = Math.round(N.clamp(inp.nz, 4, 20)), cells = nx * ny * nz;
  let dt = Math.min(inp.dt, inp.t_total), nSteps = Math.max(1, Math.ceil(inp.t_total / dt - 1e-9));
  const cap = Math.max(4, Math.floor(WORK / cells));
  if (nSteps > cap) { nSteps = cap; warnings.push(`Time step enlarged from ${f(inp.dt)} s to ${f(inp.t_total / nSteps)} s: ${cells} cells × ${Math.ceil(inp.t_total / inp.dt)} steps exceeds the in-browser work budget. Use a coarser grid or shorter duration to run the requested step.`); }
  dt = inp.t_total / nSteps;
  const sim = znoSim(P, nx, ny, nz, dt, nSteps, true, SWEEPS);
  // built-in solution-verification indicator: coarse grid at Δt and Δt/2
  const cx = Math.max(4, Math.round(nx / 2)), cy = Math.max(4, Math.round(ny / 2)), cz = Math.max(4, Math.round(nz / 2));
  const c1 = znoSim(P, cx, cy, cz, dt, nSteps, false, SWEEPS / 8), c2 = znoSim(P, cx, cy, cz, dt / 2, nSteps * 2, false, SWEEPS / 4);
  const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-300);
  const conv = s_ => (s_.inv.iS) / Math.max(s_.inv.invZ0 + Math.max(s_.inv.inZ, 0), 1e-300);
  const sens = { dtL: rel(c1.st.Lm, c2.st.Lm), dtT: Math.abs(c1.st.Tx - c2.st.Tx), dtX: rel(conv(c1), conv(c2)), grL: rel(sim.st.Lm, c1.st.Lm), grT: Math.abs(sim.st.Tx - c1.st.Tx), grX: rel(conv(sim), conv(c1)) };
  const ctx = { inp, P, sim, c1, c2, sens, salt, base, saltF, baseF, warnings, nx, ny, nz, dt, nSteps, cx, cy, cz, Vm, Gs, Em, nu };
  return { result: znoResult(ctx), rebuild: () => znoResult(ctx) };
}

function znoChem(saltF, baseF, salt, base, warnings) {
  const E = [];
  try {
    if (base.hmta) {
      E.push(eqn('chemical', 'HMTA hydrolysis (slow, thermally activated OH⁻ source)', chem.equationString(['C6H12N4', 'H2O'], ['HCHO', 'NH3']), 'First-order in HMTA with an Arrhenius rate constant; releases 4 NH₃ per HMTA.'));
      E.push(eqn('ionic', 'Ammonia protonation (fast buffer equilibrium)', chem.equationString(['NH3', 'H2O'], ['NH4+', 'OH-'], null, '⇌'), 'K_b = 1.8×10⁻⁵ mol/L; solved locally for free OH⁻.'));
      const ie = chem.ionicEquations([{ f: saltF, phase: 'aq' }, { f: 'NH3', phase: 'aq' }, { f: 'H2O', phase: 'l' }], [{ f: 'ZnO', phase: 's' }, { f: salt.prod.NH4, phase: 'aq' }]);
      E.push(eqn('chemical', 'Overall precipitation (molecular)', ie.molecular), eqn('ionic', 'Complete ionic equation', ie.completeIonic), eqn('ionic', 'Net ionic equation (via ammonia)', ie.netIonic, ie.spectators.length ? 'Spectator ions: ' + ie.spectators.join(', ') : undefined));
    } else {
      const ie = chem.ionicEquations([{ f: saltF, phase: 'aq' }, { f: baseF, phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'H2O', phase: 'l' }, { f: salt.prod[base.cat], phase: 'aq' }]);
      E.push(eqn('chemical', 'Overall precipitation (molecular)', ie.molecular), eqn('ionic', 'Complete ionic equation', ie.completeIonic), eqn('ionic', 'Net ionic equation', ie.netIonic, 'Spectator ions: ' + ie.spectators.join(', ')));
      if (base.weak) E.push(eqn('ionic', 'Ammonia protonation (fast buffer equilibrium)', chem.equationString(['NH3', 'H2O'], ['NH4+', 'OH-'], null, '⇌'), 'Aqueous ammonia is a weak base: free OH⁻ follows K_b = 1.8×10⁻⁵ mol/L.'));
    }
    E.push(eqn('ionic', 'Growth reaction solved in the transport model', chem.equationString(['Zn2+', 'OH-'], ['ZnO', 'H2O']), '1 Zn²⁺ : 2 OH⁻ stoichiometry is enforced in both species sinks.'));
    E.push(eqn('chemical', 'Intermediate dehydration (not resolved separately)', chem.equationString(['Zn(OH)2'], ['ZnO', 'H2O']), 'Zn(OH)₂ and zincate complexes are lumped into the single net step.'));
  } catch (e) { warnings.push('Reaction equations could not be generated automatically: ' + e.message); E.push(eqn('ionic', 'Net ionic equation', 'Zn²⁺(aq) + 2 OH⁻(aq) → ZnO(s) + H₂O(l)')); }
  return E;
}

function znoResult(C) {
  const { inp, P, sim, c1, c2, sens, salt, base, saltF, baseF, nx, ny, nz, dt, nSteps, cx, cy, cz } = C, warnings = [...C.warnings];
  const { H, st, bal, inv, num, fields, g } = sim, n = g.n, nxy = nx * ny, vol = P.Lx * P.Ly * P.Lz, area = P.Lx * P.Ly;
  const MZnO = MM('ZnO'), yieldKg = inv.iS * MZnO, znSupplied = inv.invZ0 + Math.max(inv.inZ, 0), conv = N.clamp(inv.iS / Math.max(znSupplied, 1e-300), 0, 1);
  let ty = { moles: 0, limiting: 'Zn2+' }; try { ty = chem.theoreticalYield(['Zn2+', 'OH-'], ['ZnO', 'H2O'], [inv.invZ0, Math.max(inv.invB0 + 4 * inv.invH0, 1e-300)]); } catch (e) { /* keep default */ }
  const yieldFrac = ty.moles > 0 ? inv.iS / ty.moles : 0;
  const m = sim.mech, sx = Math.abs(m.sx), ux = Math.abs(m.ux), dpMax = st.px - P.p0, Tb = tSat(st.px), ttot = inp.t_total;
  const PJmean = inv.EJ / ttot, kWhkg = yieldKg > 1e-12 ? inv.EJ / 3.6e6 / yieldKg : 0, aspect = st.dm > 0 ? st.Lm / st.dm : 0, Gmean = st.Lm / ttot;
  const mu = muW(N.clamp(st.Tm, 255, 650)), alphaT = P.kth / P.rhocp, Fo = alphaT * dt / g.dz ** 2, mobZ = 2 * F / (R * P.T0) * P.DZ * Math.sqrt(P.phi) * Math.abs(sim.Efield), mobB = P.weak ? 0 : F / (R * P.T0) * P.DOH * Math.sqrt(P.phi) * Math.abs(sim.Efield);
  const CFL = Math.max(mobZ, mobB) * dt / g.dz, FoZ = P.phi ** 0.5 * P.DZ * dt / g.dz ** 2, Dp = P.kappa / mu / P.Sp, tauT = P.rhocp * P.Lz / (2 * P.h), cover = st.ns * PI / 4 * st.dm ** 2, dphi = st.s * P.Vm;
  const sigBulk0 = sim.sigMean, lenD = Math.sqrt(P.phi ** 0.5 * P.DZ * ttot), volts = sim.Efield * P.Lz;
  const ctrlOn = P.ctrlT || P.ctrlG, settle = P.ctrlT ? sim.tSettle : 0;

  const kpis = [
    kpi('L_mean', 'Mean nanorod length', st.Lm * 1e9, 'nm', `Nucleus-weighted mean rod length over the whole bed at t = ${f(ttot)} s, from the per-rod volume balance with axial/lateral growth anisotropy ${f(P.beta)}. It is a homogenised (cell-average) quantity, not a resolved rod.`, null),
    kpi('L_max', 'Maximum local rod length', st.Lmax * 1e9, 'nm', 'Largest cell-average rod length anywhere in the reactor; compare with the mean to see where growth runs ahead (usually where precursor is replenished or nuclei are sparse).'),
    kpi('d_mean', 'Mean nanorod diameter', st.dm * 1e9, 'nm', `Mean rod diameter from the same volume balance, starting from ${f(inp.d_seed)} nm nuclei. Lateral growth is 2/β of the axial growth by construction.`),
    kpi('aspect', 'Aspect ratio L/d', aspect, '–', 'Mean length divided by mean diameter. Values well above 1 indicate rod-like morphology; the upper limit is set by the anisotropy ratio and by precursor supply per nucleus.', aspect >= 3),
    kpi('nonuniformity', 'Length non-uniformity', st.cvL * 100, '%', 'Coefficient of variation of the cell-average rod length across the reactor (nucleus-weighted). It combines depth-wise precursor gradients, site-density variation, field non-uniformity and edge cooling.', st.cvL < 0.15),
    kpi('n_sites', 'Nucleation density', st.ns, '1/m²', `Activated nuclei per m² of pore-solid surface (bed average); the available site density was ${f(inp.Ns)} 1/m². Computed from classical heterogeneous nucleation with site exhaustion.`),
    kpi('G_ax', 'Mean axial growth rate', Gmean * 6e10, 'nm/min', 'Final mean length divided by total time. The instantaneous rate is plotted versus time; it falls as supersaturation is consumed.'),
    kpi('conversion', 'Zn²⁺ conversion', conv * 100, '%', 'Moles of ZnO formed divided by zinc initially in the pores plus zinc that entered through the boundaries. Low values mean the run is kinetically limited; values near 100 % mean supply-limited growth.', conv > 0.5),
    kpi('yield', 'ZnO yield', yieldKg, 'kg', `Integrated ZnO product in the ${f(vol)} m³ bed (${f(yieldFrac * 100)} % of the stoichiometric yield of the initial charge, limiting reagent ${chem.prettyFormula(ty.limiting)}).`),
    kpi('T_max', 'Maximum temperature', st.Tx, 'K', `Hottest cell at the end of the run (${f(st.Tx - 273.15)} °C); saturation temperature at the local pore pressure is ${f(Tb)} K. The model has no boiling, so results above that limit are not physical.`, st.Tx < Tb - 2),
    kpi('P_joule', 'Mean Joule power', PJmean, 'W', `Time-averaged ∫σ_eff E² dV with the ionic conductivity from Nernst–Einstein (final bed-average σ_eff = ${f(sigBulk0)} S/m) and current continuity through the bed thickness.`),
    kpi('dp_max', 'Maximum excess pore pressure', dpMax, 'Pa', `Peak pore pressure above the ${f(P.p0)} Pa reference caused by thermal pressurisation against drainage to the top surface (pressure diffusivity ${f(Dp)} m²/s).`),
    kpi('stress_max', 'Maximum lateral stress', sx, 'Pa', `Largest magnitude of the confined in-plane total stress σxx = σyy (${m.sx < 0 ? 'compressive' : 'tensile'}) from thermal expansion, pore pressure and growth eigenstrain under oedometric constraint.`),
    kpi('u_max', 'Maximum vertical displacement', ux, 'm', 'Heave of the free top surface relative to the fixed base, obtained by integrating the vertical strain through each column.'),
    kpi('energy_kg', 'Electrical energy per kg ZnO', kWhkg, 'kWh/kg', 'Joule energy dissipated in the bed divided by ZnO produced. It excludes pumps, pre-heating and downstream washing/drying, so it is a lower bound on process energy.', kWhkg > 0 ? kWhkg < 50 : null),
    kpi('dt_sens', 'Time-step sensitivity of mean length', sens.dtL * 100, '%', `Change in mean rod length when the ${cx}×${cy}×${cz} coarse grid is re-run with Δt/2. A value of a few per cent or less indicates the time step is not controlling the answer; it is an indicator, not a convergence proof.`, sens.dtL < 0.05),
    kpi('grid_sens', 'Grid sensitivity of mean length', sens.grL * 100, '%', `Change in mean rod length between the ${nx}×${ny}×${nz} grid and the ${cx}×${cy}×${cz} grid at the same Δt. Large values mean the electromigration fronts (first-order upwind) are under-resolved.`, sens.grL < 0.1),
  ];
  if (ctrlOn) kpis.push(kpi('ctrl_signal', P.ctrlT ? 'Final controlled field' : 'Final dosing multiplier', P.ctrlT ? sim.Efield : sim.dose, P.ctrlT ? 'V/m' : '×', P.ctrlT ? `Electric field commanded by the PID loop at the end of the run to hold ${f(P.Tset)} K; the last excursion beyond 2 % of the step occurred at ${f(settle)} s.` : `Factor applied to the prescribed bottom concentrations by the PI dosing loop to hold ${f(inp.G_set)} nm/min.`));

  const t = H.t, L = (name, y) => ({ name, x: t, y });
  const series = [
    ser('species', 'Bed-average species concentrations', 'time (s)', 'concentration (mol/m³)', [L('Zn²⁺ (pore fluid)', H.cZ), L('free OH⁻ (pore fluid)', H.cB), L('ZnO (per bed volume)', H.s)], 'Consumption of zinc and hydroxide and accumulation of ZnO. Zinc and hydroxide should fall in a 1:2 ratio apart from boundary supply and electromigration.'),
    ser('temp', 'Temperature history', 'time (s)', 'temperature (K)', [L('mean', H.Tm), L('maximum', H.Tx)].concat(P.ctrlT ? [{ name: 'setpoint', x: [t[0], t[t.length - 1]], y: [P.Tset, P.Tset] }] : []), `Joule and reaction heating against convective loss (h = ${f(P.h)} W/m²K); the slab thermal time constant is about ${f(tauT)} s.`),
    ser('growth', 'Nanorod dimensions', 'time (s)', 'size (nm)', [L('mean length', H.L), L('mean diameter', H.d)], 'Evolution of the homogenised rod length and diameter. Mean length can dip when many new small nuclei appear.'),
    ser('rate', 'Axial growth rate and nucleation', 'time (s)', 'growth rate (nm/min)', [L('axial growth rate', H.G)], 'Instantaneous mean axial growth rate; it peaks once nuclei exist and supersaturation is still high, then decays.'),
    ser('nucl', 'Nucleation density', 'time (s)', 'nuclei per m² of solid', [L('activated nuclei', H.ns)], 'Heterogeneous nucleation history; the plateau is the available site density unless supersaturation collapses first.'),
    ser('press', 'Pore pressure and stress', 'time (s)', 'Pa', [L('max excess pore pressure', H.px.map(v => v - P.p0)), L('max |lateral stress|', H.sx)], 'Thermal pressurisation of the pore fluid and the confined lateral stress. Both follow the heating rate and temperature rise respectively.'),
    ser('control', 'Field and dosing signals', 'time (s)', 'V/m (field) / × (dose)', [L('electric field (V/m)', H.E), L('dosing multiplier × 1000', H.dose.map(v => v * 1000))], ctrlOn ? 'Actuator signals of the feedback loops (field modulation for temperature, bottom-concentration dosing for growth rate).' : 'Open-loop run: constant field and constant boundary concentrations. Enable a control mode to modulate them.'),
    ser('power', 'Heat sources', 'time (s)', 'power (W)', [L('Joule heating', H.PJ), L('reaction heat', H.PR)], 'Total Joule power and reaction heat release in the bed. Joule power changes with temperature and ion depletion through the conductivity.'),
    ser('solver', 'Nonlinear solver history', 'time (s)', 'iterations / residual', [L('Picard iterations', H.it), L('log₁₀ final update', H.res.map(v => Math.log10(v)))], 'Fixed-point iterations per time step and the size of the last update (relative, max-norm). Steps that hit the iteration cap are counted in the checks.'),
  ];
  const kmid = Math.floor(nz / 2), jmid = Math.floor(ny / 2);
  const plane = (a, sc = 1) => { const d = new Float64Array(nxy); for (let c = 0; c < nxy; c++) d[c] = a[c + kmid * nxy] * sc; return d; };
  const slice = (a, sc = 1, off = 0) => { const d = new Float64Array(nx * nz); for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) d[i + nx * k] = (a[i + nx * (jmid + ny * k)] - off) * sc; return d; };
  const mp = (key, title, data, unit, explain, vert) => ({ key, title, nx, ny: vert ? nz : ny, data, x0: 0, x1: P.Lx, y0: 0, y1: vert ? P.Lz : P.Ly, xlabel: 'x (m)', ylabel: vert ? 'z (m)' : 'y (m)', unit, explain, diverging: false });
  const maps = [
    mp('T_mid', 'Temperature, mid-plane', plane(fields.T), 'K', 'Horizontal mid-plane temperature at the end of the run; edges are cooler through side-wall convection and the field profile.'),
    mp('L_mid', 'Rod length, mid-plane', plane(fields.Lr, 1e9), 'nm', 'Cell-average rod length on the horizontal mid-plane; the texture follows the seeded site-density field and the temperature map.'),
    mp('cZ_xz', 'Zn²⁺ concentration, vertical slice', slice(fields.cZ), 'mol/m³', 'Vertical x–z slice through the centre: depletion by reaction and re-supply from the bottom by diffusion and electromigration (cations move towards +z).', true),
    mp('cOH_xz', 'Free OH⁻ concentration, vertical slice', slice(fields.xf), 'mol/m³', 'Vertical slice of free hydroxide; anions migrate towards −z, so the depletion pattern is the mirror image of the zinc profile.', true),
    mp('T_xz', 'Temperature, vertical slice', slice(fields.T), 'K', 'Through-thickness temperature; the Biot number of the slab is small, so the profile is nearly flat with slightly cooler faces.', true),
    mp('p_xz', 'Excess pore pressure, vertical slice', slice(fields.p, 1, P.p0), 'Pa', 'Excess pore pressure relative to the reference: zero at the drained top, largest at the sealed base while the bed is heating.', true),
    mp('L_xz', 'Rod length, vertical slice', slice(fields.Lr, 1e9), 'nm', 'Depth-wise growth uniformity: longer rods where precursor is replenished, shorter where the reactants are exhausted.', true),
  ];

  // ---- nanorod array geometry sampled from the computed fields (display statistics are seeded) ----
  const rnd = N.rng((Math.round(inp.seed) || 1) + 7), nRods = 140, nsF = Math.max(st.ns, 1), side = N.clamp(Math.sqrt(nRods / nsF) * 1e9, 50, 1e6), cvR = inp.rod_cv / 100;
  const samp = Array.from({ length: nRods }, () => { const u = rnd(), v = rnd(); return { u, v, c: Math.min(nx - 1, Math.floor(u * nx)) + nx * Math.min(ny - 1, Math.floor(v * ny)), fl: Math.exp(rnd.normal(0, cvR)), fd: Math.exp(rnd.normal(0, cvR * 0.6)), tilt: Math.abs(rnd.normal(0, inp.tilt_sd)), az: 360 * rnd(), act: rnd() }; });
  const last = sim.frames[sim.frames.length - 1];
  const rodsAt = fr => samp.map(r => { const show = last.a[r.c] > 0 && fr.a[r.c] / last.a[r.c] >= r.act; return { x: +(r.u * side).toFixed(2), y: +(r.v * side).toFixed(2), length: show ? +(fr.L[r.c] * 1e9 * r.fl).toFixed(3) : 0, diameter: show ? +(fr.d[r.c] * 1e9 * r.fd).toFixed(3) : 0, tilt: +r.tilt.toFixed(2), azimuth: +r.az.toFixed(1) }; });
  const geometry = { kind: 'rods', rods: rodsAt(last), area: [side, side], units: 'nm', title: `ZnO nanorod array (representative ${f(side)} nm patch; rods sampled from the reactor-plan fields)`, frames: sim.frames.map(fr => ({ t: fr.t, rods: rodsAt(fr) })) };

  const equations = znoChem(saltF, baseF, salt, base, warnings).concat([
    eqn('governing', 'Species transport (Nernst–Planck, porous medium)', 'φ ∂cᵢ/∂t = ∇·(D_eff,i ∇cᵢ) − ∂/∂z( zᵢ F D_eff,i E_z cᵢ /(R T) ) − νᵢ r,   D_eff,i = φ^1.5 Dᵢ(T)', 'Bruggeman tortuosity; Dᵢ(T) by Stokes–Einstein scaling with water viscosity; ν = 1 (Zn²⁺), 2 (OH⁻).'),
    eqn('governing', 'Surface reaction rate', 'r = k(T) a_rod (c_Zn c_OH² − K_sp(T)),   k = A exp(−E_a/RT),   K_sp(T) = K_sp,298 exp[−ΔH_diss/R (1/T − 1/298.15)]', 'Mass-action driving force that vanishes at the solubility product; a_rod is the rod surface area per bed volume.'),
    eqn('governing', 'Heterogeneous nucleation with site exhaustion', '∂n/∂t = (N_s − n) ω₀ exp[−E_a/R (1/T − 1/298.15)] exp[−f(θ) ΔG*/(k_B T)],   ΔG* = 16π γ³ v_m² / (3 (k_B T ln S)²),   S = c_Zn c_OH²/K_sp', 'f(θ) = (2 + cos θ)(1 − cos θ)²/4.'),
    eqn('constitutive', 'Nanorod growth closure', 'n (π/4) d² L = c_ZnO V_m,   L − d₀ = (β/2)(d − d₀),   G_ax = ∂L/∂t,   G_lat = G_ax/β', 'Per-rod volume balance with constant axial/lateral rate ratio β; mass-consistent with the ZnO field.'),
    eqn('governing', 'Energy balance', 'ρ c_p ∂T/∂t = ∇·(k ∇T) + σ_eff E_z² + (−ΔH) r;   −k ∂T/∂n = h (T − T∞) on all faces'),
    eqn('constitutive', 'Electrical closure', 'σ_eff = φ^1.5 (F²/RT) Σ zᵢ² Dᵢ cᵢ,   j_z = V / ∫dz/σ_eff,   E_z = j_z/σ_eff', 'Nernst–Einstein conductivity; current continuity through the thickness with V = E_applied · L_z.'),
    eqn('governing', 'Pore pressure', '(1/M + α²/M_c) ∂p/∂t = ∇·(κ/μ(T) ∇p) + β_e ∂T/∂t;   p = p₀ on top, no flux elsewhere', 'Uniaxial-strain storage; β_e = φ β_f + (α − φ) β_s − α K_d β_s/M_c.'),
    eqn('governing', 'Poroelastic mechanics (oedometric column)', 'σ_zz = M_c ε_zz − α Δp − 3K_d(α_T ΔT + ε_g) = 0,   σ_xx = λ ε_zz − α Δp − 3K_d(α_T ΔT + ε_g),   u_z = ∫₀ᶻ ε_zz dz', 'Bottom fixed, top traction-free, lateral strain zero; ε_g = η_g c_ZnO V_m/3.'),
    eqn('formula', 'Control laws', 'E = E₀ + K_p e + K_i ∫e dt + K_d de/dt (e = T_set − T̄);   c_bottom = m·c₀, m = 1 + K_g e_g + K_ig ∫e_g dt (e_g = (G_set − G)/G_set)'),
  ]);

  const tolBal = 1e-6;
  const checks = [
    chk('Zinc inventory (pores + ZnO − boundary flux)', bal.zn < tolBal, `relative residual ${f(bal.zn)}`),
    chk('Hydroxide inventory incl. boundary flux and reaction sink', bal.oh < tolBal, `relative residual ${f(bal.oh)}`),
    chk('1:2 stoichiometry closure (OH⁻ consumed vs 2 × ZnO formed)', bal.stoich < 1e-3, `relative mismatch ${f(bal.stoich)} (set by the Picard tolerance)`),
    chk('Energy balance (storage = Joule + reaction − convection)', bal.energy < tolBal, `relative residual ${f(bal.energy)}`),
    chk('Positivity of concentrations', bal.minC >= 0, `minimum concentration ${f(bal.minC)} mol/m³`),
    chk('Nonlinear (Picard) convergence', num.resMax < 1e-3, `${num.itTot} iterations over ${num.nSteps} steps (max ${num.itMax}/step, ${num.damped} under-relaxed); ${num.notConv} step(s) ended above the 1×10⁻⁴ tolerance; largest final relative update ${f(num.resMax)} (pass limit 1×10⁻³)`),
    chk('Linear solves', num.linMax < 1e-8, `largest relative residual ${f(num.linMax)} after line relaxation (${f(num.swTot)} sweep-equivalents in total)`),
    chk('Below boiling limit at pore pressure', st.Tx < Tb, `T_max ${f(st.Tx, 4)} K vs T_sat ${f(Tb, 4)} K`),
    chk('Time-step sensitivity indicator (Δt vs Δt/2, coarse grid)', sens.dtL < 0.05 && sens.dtT < 1, `mean length ${f(sens.dtL * 100)} %, T_max ${f(sens.dtT)} K, conversion ${f(sens.dtX * 100)} %`),
    chk('Grid sensitivity indicator (fine vs coarse grid)', sens.grL < 0.1 && sens.grT < 1, `mean length ${f(sens.grL * 100)} %, T_max ${f(sens.grT)} K, conversion ${f(sens.grX * 100)} %`),
    chk('Rod coverage below coalescence', cover < 0.9, `projected rod area fraction on the pore surface ${f(cover)}`),
  ];
  if (P.hmta) checks.push(chk('HMTA inventory', bal.hmta < tolBal, `relative residual ${f(bal.hmta)}`));
  try { const a = chem.auditEquation(['Zn2+', 'OH-'], ['ZnO', 'H2O'], [1, 2, 1, 1]); checks.push(chk('Charge and mass balance of the net reaction', a.massBalanced && a.chargeBalanced, `elements balanced: ${a.massBalanced}; net charge ${a.netCharge}`)); } catch (e) { /* not essential */ }

  warnings.push('Reduced-order research model, not a validated reactor simulator: the rate constant, nucleation frequency, interfacial energy and growth anisotropy are closure parameters that must be calibrated against measured conversion and rod dimensions before the numbers are used quantitatively.');
  warnings.push(`The ${f(g.dx * 1e3)} × ${f(g.dy * 1e3)} × ${f(g.dz * 1e3)} mm cells do not resolve individual nanorods; length, diameter and density are homogenised cell averages and the displayed array is sampled from them with a ${f(inp.rod_cv)} % display scatter.`);
  warnings.push('The electric field is imposed (parallel electrodes, current continuity through the thickness); electroneutrality/Poisson coupling, electrode reactions and electro-osmosis are not solved.' + (volts > 1.23 ? ` The applied ${f(volts)} V exceeds the 1.23 V water-splitting potential, so a DC field would electrolyse the solution — use AC or membrane-separated electrodes.` : ''));
  warnings.push('Speciation is lumped: Zn(OH)₂, Zn(OH)₄²⁻ and ammine complexes, activity coefficients and pH-dependent solubility are not resolved, so behaviour at high pH (zincate re-dissolution) is not captured.');
  if (st.Tx >= Tb) warnings.push(`Predicted maximum temperature ${f(st.Tx - 273.15)} °C exceeds the boiling point at the pore pressure (${f(Tb - 273.15)} °C): the applied field cannot be sustained at this pressure without boiling. Pressurise the vessel (autoclave), lower the field or use the temperature controller.`);
  if (CFL > 5) warnings.push(`Electromigration Courant number ${f(CFL)} with first-order upwinding: fronts are smeared by numerical diffusion. Check the grid/time-step sensitivity indicators before trusting concentration profiles.`);
  if (num.capped) warnings.push(`Work budget reached: ${num.capped} time step(s) were limited to two coupling iterations and minimal linear relaxation. Use a larger time step or a coarser grid for this case.`);
  if (num.notConv) warnings.push(`${num.notConv} time step(s) reached the Picard iteration cap; reduce Δt or raise the cap.`);
  if (bal.stoich >= 1e-3) warnings.push('Hydroxide consumption deviates from twice the ZnO formed by more than 0.1 %: the reaction coupling is not fully converged at this time step.');
  if (cover >= 0.9) warnings.push('Rods have coalesced on the pore surface (projected coverage ≥ 0.9); the isolated-rod closure is no longer valid — a film is forming.');
  if (st.Lmax > 0.25 * inp.dgrain * 1e-3 * P.phi) warnings.push('Rod length is a significant fraction of the pore size; pore clogging and evolving porosity (not modelled) will matter.');
  if (dphi > 0.05 * P.phi) warnings.push(`Deposited ZnO fills ${f(dphi / P.phi * 100)} % of the pore space; porosity is kept constant in the transport operators.`);
  if (base.weak) warnings.push('Weak-base route: the transported base field is total ammonia + hydroxide with a local buffer equilibrium; ammonium is assumed to stay where it is formed.');
  if (sens.dtL >= 0.05 || sens.grL >= 0.1) warnings.push('The built-in sensitivity indicators show a visible dependence on time step or grid; a single mesh/time-step is not convergence evidence — refine before drawing conclusions.');

  const recs = [];
  if (st.Tx >= Tb - 5) recs.push(rec('high', `Limit the field or pressurise: at ${f(sim.Efield)} V/m the bed reaches ${f(st.Tx - 273.15)} °C. Hold temperature with the PID mode (set-point ≤ ${f(Tb - 273.15 - 10)} °C at 1 atm) or run in a sealed autoclave rated for the saturation pressure.`, 'Boiling-limit check against the Antoine saturation temperature.'));
  if (conv < 0.5) recs.push(rec('high', `Only ${f(conv * 100)} % of the supplied zinc is converted in ${f(ttot)} s. Raise the bed temperature (the rate constant roughly doubles every ${f(Math.LN2 * R * st.Tm ** 2 / P.Ea)} K at E_a = ${f(inp.Ea)} kJ/mol) or extend the dwell.`, 'Arrhenius sensitivity of the computed conversion.'));
  else recs.push(rec('medium', `Conversion is ${f(conv * 100)} %: growth is supply-limited. Longer rods need more precursor per nucleus — lower the seed-site density (now ${f(st.ns)} 1/m²) or dose precursor continuously with the growth-rate controller.`, 'Per-rod volume balance: rod volume = ZnO formed / nuclei.'));
  recs.push(rec(st.cvL > 0.15 ? 'high' : 'medium', `Length non-uniformity is ${f(st.cvL * 100)} %. ${st.cvL > 0.15 ? 'Reduce bed thickness or reverse/alternate the field polarity so both reactants are replenished through the depth, and tighten seed-layer uniformity.' : 'Keep seed-layer uniformity and bed thickness at the present level; verify with cross-sectional SEM at three depths.'}`, 'Computed coefficient of variation of rod length over the bed.'));
  recs.push(rec('medium', `Calibrate before use: measure Zn²⁺ depletion (ICP-OES) and rod length (SEM) at three temperatures to fit A and E_a, and count rods per µm² to fix the site density; the present values are closure assumptions.`, 'Model-limitation statement of the reference case: growth and nucleation closures require experimental calibration.'));
  recs.push(rec('medium', `Repeat on a finer grid and half the time step before quoting numbers: the indicators give ${f(sens.dtL * 100)} % (Δt) and ${f(sens.grL * 100)} % (grid) change in mean length.`, 'Built-in solution-verification indicator.'));
  recs.push(rec('medium', kWhkg > 0 ? `Electrical energy intensity is ${f(kWhkg)} kWh per kg ZnO because the field heats the whole wet bed. Recover heat from the spent liquor, insulate the vessel (h = ${f(P.h)} W/m²K now) and switch the field off once conversion plateaus to cut energy use.` : 'No Joule energy is used in this run; heat comes from the surroundings. Insulate the vessel and recover heat from the spent liquor to minimise oven energy.', 'Energy balance of the run (Joule input versus convective loss).', true));
  recs.push(rec('low', `Recycle the mother liquor: ${f((1 - conv) * 100)} % of the zinc and the ${salt.an.replace(/\^?\d?[-+]$/, '')}-salt by-product remain in solution; neutralise and precipitate residual zinc before discharge (zinc is ecotoxic to aquatic life).`, 'Zinc inventory at the end of the run.', true));
  if (ctrlOn) recs.push(rec('low', P.ctrlT ? `The temperature loop last left the ±2 % band at ${f(settle)} s; if this is too slow raise K_p (now ${f(P.Kp)} V/m/K) up to the point where the field saturates at ${f(P.Emax)} V/m.` : 'Tune the dosing gain so the multiplier stays inside its limits; a saturated dosing signal means the set-point is not reachable with the present bottom supply.', 'Closed-loop response of this run.'));

  const publish = {
    'growth.rate': Gmean, 'growth.length': st.Lm, 'growth.diameter': st.dm, 'geom.length': st.Lm, 'geom.diameter': st.dm, 'geom.aspect': aspect, 'geom.number_density': st.ns, 'nucl.rate': sim.Jpeak,
    'rxn.conversion': conv, 'rxn.yield_kg': yieldKg, 'rxn.heat': inv.ER / ttot, 'therm.T_max': st.Tx, 'mech.max_stress': sx, 'mech.max_displacement': ux, 'proc.power': PJmean, 'proc.E_field': sim.Efield, 'proc.time': ttot,
    'proc.T_C': st.Tm - 273.15, 'proc.P': st.px, 'proc.c_precursor': P.cZ0, 'proc.c_base': P.cB0, 'proc.volume': vol, 'fab.energy_kWh_kg': kWhkg, 'fab.yield': N.clamp(yieldFrac, 0, 1), 'fab.throughput_kg_h': yieldKg / ttot * 3600, 'eff.sigma': sigBulk0, 'geom.porosity': P.phi, 'material.formula': 'ZnO', 'material.name': 'Zinc oxide nanorods',
  };
  if (P.ctrlT) publish['ctrl.settling_s'] = settle;

  const tables = [
    { key: 'balance', title: 'Material and energy balances', columns: ['quantity', 'initial', 'boundary in', 'reacted / source', 'final', 'unit'], rows: [
      ['Zn²⁺ in pores', inv.invZ0, inv.inZ, -inv.iS, inv.iZ, 'mol'], ['base (OH⁻ / total NH₃+OH⁻)', inv.invB0, inv.inB, inv.srcB - inv.consB, inv.iB, 'mol'], ['ZnO', 0, 0, inv.iS, inv.iS, 'mol'],
      ['thermal energy above T₀', 0, -inv.Eloss, inv.EJ + inv.ER, inv.dU, 'J'], ['  of which Joule', 0, 0, inv.EJ, inv.EJ, 'J'], ['  of which reaction', 0, 0, inv.ER, inv.ER, 'J']].map(r => r.map(v => typeof v === 'number' ? +v.toPrecision(6) : v)), explain: 'Integrated inventories over the run. For every row, final = initial + boundary in + reacted/source to within the residuals listed in the checks.' },
    { key: 'verif', title: 'Solution-verification indicator', columns: ['run', 'grid', 'Δt (s)', 'mean length (nm)', 'T_max (K)', 'conversion (%)'], rows: [
      ['main', `${nx}×${ny}×${nz}`, +dt.toPrecision(4), +(st.Lm * 1e9).toPrecision(5), +st.Tx.toPrecision(6), +(conv * 100).toPrecision(4)],
      ['coarse', `${cx}×${cy}×${cz}`, +dt.toPrecision(4), +(c1.st.Lm * 1e9).toPrecision(5), +c1.st.Tx.toPrecision(6), +(c1.inv.iS / Math.max(c1.inv.invZ0 + Math.max(c1.inv.inZ, 0), 1e-300) * 100).toPrecision(4)],
      ['coarse, Δt/2', `${cx}×${cy}×${cz}`, +(dt / 2).toPrecision(4), +(c2.st.Lm * 1e9).toPrecision(5), +c2.st.Tx.toPrecision(6), +(c2.inv.iS / Math.max(c2.inv.invZ0 + Math.max(c2.inv.inZ, 0), 1e-300) * 100).toPrecision(4)]], explain: 'Three runs of the same case. Differences between rows 1–2 indicate grid dependence and between rows 2–3 time-step dependence; they are indicators, not a formal convergence study.' },
    { key: 'scales', title: 'Characteristic numbers', columns: ['quantity', 'value', 'unit'], rows: [
      ['thermal Fourier number α Δt/Δz²', +Fo.toPrecision(4), '–'], ['Zn²⁺ diffusion number D Δt/Δz²', +FoZ.toPrecision(4), '–'], ['electromigration Courant number', +CFL.toPrecision(4), '–'], ['Zn²⁺ diffusion length √(D t)', +lenD.toPrecision(4), 'm'],
      ['Zn²⁺ migration velocity', +mobZ.toPrecision(4), 'm/s'], ['slab thermal time constant ρc_pL/(2h)', +tauT.toPrecision(4), 's'], ['pressure diffusivity', +Dp.toPrecision(4), 'm²/s'], ['applied voltage across bed', +volts.toPrecision(4), 'V'], ['pore-solid specific surface', +P.as.toPrecision(4), 'm²/m³']], explain: 'Dimensionless step sizes and physical scales of this run. The implicit scheme is stable for any Fourier/Courant number, but accuracy degrades when they are much larger than one.' },
  ];

  const recipe = {
    process: 'Hydrothermal growth of ZnO nanorods in a porous bed with electric-field heating', chemistry: { zinc_salt: saltF, base: baseF, Zn_mol_m3: P.cZ0, base_mol_m3: P.cB0, HMTA_mol_m3: P.cH0 },
    reactor: { size_m: [P.Lx, P.Ly, P.Lz], porosity: P.phi, grain_size_mm: inp.dgrain, permeability_m2: P.kappa, boundary_set: inp.bc },
    operating: { T_initial_K: P.T0, T_ambient_K: P.Tinf, pressure_Pa: P.p0, field_V_m: P.E0, duration_s: ttot, control: inp.ctrl, T_setpoint_K: P.ctrlT ? P.Tset : null, growth_setpoint_nm_min: P.ctrlG ? inp.G_set : null },
    predicted: { mean_length_nm: st.Lm * 1e9, mean_diameter_nm: st.dm * 1e9, aspect_ratio: aspect, nonuniformity_pct: st.cvL * 100, conversion: conv, ZnO_kg: yieldKg, T_max_K: st.Tx, joule_kWh: inv.EJ / 3.6e6, max_stress_Pa: sx },
    numerics: { grid: [nx, ny, nz], dt_s: dt, steps: nSteps, sensitivity: sens }, status: 'model prediction with uncalibrated closures — confirm experimentally',
  };
  const txt = [`RECIPE — ZnO nanorods, hydrothermal, ${salt.name} + ${base.name}`, `1. Prepare ${f(P.phi * vol * 1e3)} L of solution: ${f(P.cZ0 / 1e3)} mol/L ${chem.prettyFormula(saltF)} (${f(P.cZ0 * MM(saltF) * P.phi * vol * 1e3)} g)` + (base.hmta ? ` and ${f(P.cH0 / 1e3)} mol/L HMTA (${f(P.cH0 * 0.14019 * P.phi * vol * 1e3)} g).` : ` and ${f(P.cB0 / 1e3)} mol/L ${chem.prettyFormula(baseF)}.`),
    `2. Fill the seeded porous bed (${f(P.Lx)} × ${f(P.Ly)} × ${f(P.Lz)} m, porosity ${f(P.phi)}) at ${f(P.T0 - 273.15)} °C and ${f(P.p0)} Pa.`,
    `3. ${P.ctrlT ? `Apply the field under PID control to hold ${f(P.Tset - 273.15)} °C (limit ${f(P.Emax)} V/m).` : P.E0 > 0 ? `Apply ${f(P.E0)} V/m across the thickness (${f(volts)} V).` : `Hold in surroundings at ${f(P.Tinf - 273.15)} °C.`} Dwell ${f(ttot / 60)} min.`,
    `4. Expected: rods ${f(st.Lm * 1e9)} nm long × ${f(st.dm * 1e9)} nm diameter (aspect ${f(aspect)}), ${f(conv * 100)} % zinc conversion, ${f(yieldKg * 1e3)} g ZnO, peak ${f(st.Tx - 273.15)} °C.`,
    '5. Cool, drain, rinse with deionised water until conductivity is stable, dry at 60–80 °C.', 'Status: model prediction with uncalibrated kinetic closures; verify by SEM and ICP before scale-up.'].join('\n');
  const csv = ['t_s,cZn_mol_m3,cOH_mol_m3,cZnO_mol_m3,T_mean_K,T_max_K,L_mean_nm,d_mean_nm,nuclei_per_m2,growth_nm_min,p_max_Pa,stress_max_Pa,u_max_m,E_V_m,dose,P_joule_W,P_rxn_W,picard_iters'].concat(t.map((tv, i) => [tv, H.cZ[i], H.cB[i], H.s[i], H.Tm[i], H.Tx[i], H.L[i], H.d[i], H.ns[i], H.G[i], H.px[i], H.sx[i], H.ux[i], H.E[i], H.dose[i], H.PJ[i], H.PR[i], H.it[i]].map(v => +(+v).toPrecision(7)).join(','))).join('\n');

  const summary = `After ${f(ttot)} s the bed holds ZnO rods of mean length ${f(st.Lm * 1e9)} nm and diameter ${f(st.dm * 1e9)} nm (aspect ratio ${f(aspect)}, non-uniformity ${f(st.cvL * 100)} %) at ${f(st.ns)} nuclei/m², with ${f(conv * 100)} % zinc conversion and ${f(yieldKg * 1e3)} g of ZnO. ` +
    `The ${f(sim.Efield)} V/m field dissipates ${f(PJmean)} W on average and the bed peaks at ${f(st.Tx)} K (${st.Tx < Tb ? 'below' : 'ABOVE'} the ${f(Tb)} K boiling limit); peak excess pore pressure is ${f(dpMax)} Pa and lateral stress ${f(sx / 1e3)} kPa. ` +
    `Zinc and energy balances close to ${f(Math.max(bal.zn, bal.energy))} relative; halving the time step on the coarse grid changes the mean length by ${f(sens.dtL * 100)} %. Kinetic closures are uncalibrated, so treat the dimensions as model predictions.`;
  return finish({
    summary, kpis, series, maps, tables, equations, geometry, checks, warnings, recommendations: recs, publish,
    exports: [{ name: 'zno_reactor_recipe.json', mime: 'application/json', content: JSON.stringify(recipe, null, 2), description: 'Machine-readable recipe: chemistry, reactor, operating conditions, predictions and numerical settings.' }, { name: 'zno_reactor_recipe.txt', mime: 'text/plain', content: txt, description: 'Human-readable step-by-step recipe with predicted outcome.' }, { name: 'zno_reactor_history.csv', mime: 'text/csv', content: csv, description: 'Time histories of species, temperature, growth, pressure, stress, control signals and solver iterations.' }],
    numerics: { method: 'Finite-volume, backward Euler; Picard-coupled fields (Patankar-implicit reaction sinks, first-order upwind electromigration); line-implicit (Thomas) relaxation for species, z-line-preconditioned conjugate gradients for temperature and pressure', dof: n * (P.hmta ? 8 : 7), steps: nSteps, iterations: num.itTot, residual: num.resMax, notes: `Grid ${nx}×${ny}×${nz}, Δt = ${f(dt)} s; Fourier ${f(Fo)}, Courant ${f(CFL)}; linear residual ≤ ${f(num.linMax)}; verification indicator from two extra ${cx}×${cy}×${cz} runs (Δt, Δt/2).` },
  });
}
const ZI = (o) => ({ group: 'ZnO reactor', showIf: { mode: ['zno_reactor'] }, ...o });
const INPUTS = [
  { key: 'mode', label: 'Process', type: 'select', value: 'zno_reactor', group: 'Model', options: [
    { value: 'zno_reactor', label: 'Hydrothermal ZnO nanorods — coupled 3-D porous reactor with field heating' }, { value: 'cvd', label: 'CVD / PECVD film growth' }, { value: 'ald', label: 'Atomic layer deposition' }, { value: 'pvd', label: 'PVD — sputtering / evaporation' },
    { value: 'litho_etch', label: 'Lithography and etch' }, { value: 'additive', label: 'Additive manufacturing of lattices' }, { value: 'solgel_spin_sinter', label: 'Sol–gel, spin coating and sintering' }, { value: 'electrospin_colloid', label: 'Electrospinning and colloidal batch recipe' }],
    help: 'Selects the fabrication process model. Each process has its own inputs, outputs, equations, checks and recipe export.' },
  { key: 'salt', label: 'Zinc salt', type: 'select', value: 'Zn(NO3)2', group: 'Chemistry', showIf: { mode: ['zno_reactor', 'electrospin_colloid'] }, options: Object.keys(SALTS).map(k => ({ value: k, label: `${chem.prettyFormula(k)} — ${SALTS[k].name}` })), help: 'Zinc source. The anion is a spectator: it changes the ionic equations, the by-product salt and the solution conductivity (Joule heating).' },
  { key: 'base', label: 'Base / hydroxide source', type: 'select', value: 'NaOH', group: 'Chemistry', showIf: { mode: ['zno_reactor', 'electrospin_colloid'] }, options: [{ value: 'NaOH', label: 'NaOH (strong base)' }, { value: 'KOH', label: 'KOH (strong base)' }, { value: 'NH4OH', label: 'NH₄OH (weak base, ammonia buffer)' }, { value: 'HMTA', label: 'HMTA (slow thermal release of NH₃/OH⁻)' }], help: 'Hydroxide source. Strong bases give OH⁻ directly; ammonia and HMTA supply it through the NH₃/NH₄⁺ buffer, HMTA only as it hydrolyses on heating.' },
  ZI({ key: 'cZn', label: 'Zn²⁺ concentration (initial and bottom boundary)', unit: 'mol/m³', value: 10, min: 0.1, max: 500, bus: 'proc.c_precursor', group: 'Chemistry', showIf: { mode: ['zno_reactor', 'electrospin_colloid'] }, help: 'Zinc concentration in the pore fluid (10 mol/m³ = 0.01 mol/L). Sets supersaturation, conductivity and the precursor available per nucleus.' }),
  ZI({ key: 'cOH', label: 'OH⁻ / base concentration (initial and bottom boundary)', unit: 'mol/m³', value: 20, min: 0.1, max: 1000, bus: 'proc.c_base', group: 'Chemistry', help: 'Hydroxide (or total ammonia for NH₄OH) concentration; 20 mol/m³ = 0.02 mol/L is the stoichiometric 1:2 ratio for the reference case. Not used for the HMTA route.' }),
  ZI({ key: 'c_hmta', label: 'HMTA concentration', unit: 'mol/m³', value: 25, min: 0.1, max: 500, group: 'Chemistry', showIf: { mode: ['zno_reactor'], base: ['HMTA'] }, help: 'Initial hexamethylenetetramine concentration (equimolar with zinc is the usual laboratory recipe). Each HMTA releases 4 NH₃ on hydrolysis.' }),
  ZI({ key: 'k_hmta90', label: 'HMTA hydrolysis rate constant at 90 °C', unit: '1/s', value: 1e-4, min: 1e-6, max: 1e-2, log: true, group: 'Chemistry', showIf: { mode: ['zno_reactor'], base: ['HMTA'] }, help: 'First-order HMTA decomposition rate at 90 °C (order 10⁻⁴ 1/s, i.e. hours). It is pH dependent in reality and should be calibrated.' }),
  ZI({ key: 'Ea_hmta', label: 'HMTA hydrolysis activation energy', unit: 'kJ/mol', value: 80, min: 20, max: 150, group: 'Chemistry', showIf: { mode: ['zno_reactor'], base: ['HMTA'] }, help: 'Temperature sensitivity of HMTA decomposition; a calibration parameter that controls how sharply OH⁻ release switches on with heating.' }),
  ZI({ key: 'T0', label: 'Initial temperature', unit: 'K', value: 298, min: 278, max: 370, group: 'Conditions', help: 'Uniform starting temperature of the bed and solution.' }),
  ZI({ key: 'Tinf', label: 'Surroundings (oven/bath) temperature', unit: 'K', value: 298, min: 278, max: 473, group: 'Conditions', help: 'Temperature of the surroundings for the convective boundary. Set it to the oven temperature for a conventional hydrothermal recipe without field heating.' }),
  ZI({ key: 'p0', label: 'Reference pressure (top surface)', unit: 'Pa', value: 101325, min: 5e4, max: 2e7, log: true, bus: 'proc.P', group: 'Conditions', help: 'Pressure fixed at the top of the bed and used as the initial pore pressure. It sets the boiling limit; autoclaves run above atmospheric.' }),
  ZI({ key: 'h', label: 'Heat-transfer coefficient', unit: 'W/m²K', value: 50, min: 1, max: 5000, log: true, group: 'Conditions', help: 'Convective coefficient on all outer faces. Larger values hold the bed closer to the surroundings and shorten the thermal time constant.' }),
  ZI({ key: 'bc', label: 'Species boundary set', type: 'select', value: 'reference', group: 'Conditions', options: [{ value: 'reference', label: 'Reference: prescribed at bottom, zero-gradient elsewhere' }, { value: 'sealed_top', label: 'Prescribed at bottom, impermeable top' }, { value: 'batch', label: 'Closed batch: no species flux on any face' }], help: 'Reference case prescribes Zn²⁺ and OH⁻ at the bottom and uses zero normal gradient elsewhere (electromigration flux then passes through the top). The alternatives close the top or all faces.' }),
  ZI({ key: 'Lx', label: 'Reactor length x', unit: 'm', value: 1, min: 0.005, max: 5, log: true, group: 'Reactor', help: 'In-plane length of the rectangular porous bed.' }),
  ZI({ key: 'Ly', label: 'Reactor width y', unit: 'm', value: 1, min: 0.005, max: 5, log: true, group: 'Reactor', help: 'In-plane width of the rectangular porous bed.' }),
  ZI({ key: 'Lz', label: 'Bed thickness z (field direction)', unit: 'm', value: 0.01, min: 0.001, max: 0.2, log: true, group: 'Reactor', help: 'Thickness between the electrodes. Sets the voltage for a given field, the thermal time constant and the distance reactants must travel.' }),
  ZI({ key: 'phi', label: 'Porosity', unit: '–', value: 0.25, min: 0.05, max: 0.95, bus: 'geom.porosity', group: 'Reactor', help: 'Pore volume fraction. Controls fluid inventory, effective diffusivity and conductivity (Bruggeman φ^1.5) and the solid surface available for nucleation.' }),
  ZI({ key: 'dgrain', label: 'Grain / strut size of the porous solid', unit: 'mm', value: 2, min: 0.001, max: 20, log: true, group: 'Reactor', help: 'Characteristic solid size; the specific surface is 6(1−φ)/d. Finer solids give more nucleation area and therefore shorter rods for the same precursor.' }),
  ZI({ key: 'kappa', label: 'Permeability', unit: 'm²', value: 1e-17, min: 1e-20, max: 1e-9, log: true, bus: 'eff.permeability', group: 'Reactor', help: 'Darcy permeability for pore-pressure diffusion. Low permeability traps thermal pressurisation; high permeability drains it instantly.' }),
  ZI({ key: 'DZn', label: 'Zn²⁺ diffusion coefficient (25 °C)', unit: 'm²/s', value: 7.5e-10, min: 1e-11, max: 1e-8, log: true, group: 'Transport', help: 'Free-solution diffusivity of Zn²⁺; also sets its electrical mobility through the Nernst–Einstein relation.' }),
  ZI({ key: 'DOH', label: 'OH⁻ diffusion coefficient (25 °C)', unit: 'm²/s', value: 5.3e-9, min: 1e-11, max: 2e-8, log: true, group: 'Transport', help: 'Free-solution diffusivity of OH⁻, the fastest ion present; it dominates the solution conductivity.' }),
  ZI({ key: 'kth', label: 'Thermal conductivity (bed)', unit: 'W/m/K', value: 0.6, min: 0.05, max: 50, log: true, bus: 'eff.k', group: 'Transport', help: 'Effective thermal conductivity of the saturated bed.' }),
  ZI({ key: 'cp', label: 'Specific heat capacity (bed)', unit: 'J/kg/K', value: 4186, min: 300, max: 5000, group: 'Transport', help: 'Effective specific heat of the saturated bed; with density it sets how fast the bed heats.' }),
  ZI({ key: 'rho', label: 'Density (bed)', unit: 'kg/m³', value: 1000, min: 200, max: 8000, bus: 'eff.density', group: 'Transport', help: 'Effective density of the saturated bed.' }),
  ZI({ key: 'Ea', label: 'Reaction activation energy', unit: 'kJ/mol', value: 65, min: 10, max: 150, group: 'Kinetics', help: 'Arrhenius activation energy of the surface growth reaction (also used for the nucleation frequency). Higher values make growth more temperature sensitive.' }),
  ZI({ key: 'A_pre', label: 'Reaction pre-exponential factor', unit: 'm⁷/mol²/s', value: 200, min: 1e-3, max: 1e7, log: true, group: 'Kinetics', help: 'Pre-exponential factor of the surface rate constant k = A·exp(−E_a/RT) in r = k·a_rod·(c_Zn·c_OH² − K_sp). Calibrate against measured conversion.' }),
  ZI({ key: 'dH_rxn', label: 'Reaction enthalpy ΔH', unit: 'kJ/mol', value: -22.4, min: -200, max: 100, group: 'Kinetics', help: 'Enthalpy of Zn²⁺ + 2 OH⁻ → ZnO + H₂O (negative = exothermic; −22.4 kJ/mol from standard formation enthalpies). Also sets the van’t Hoff temperature dependence of solubility.' }),
  ZI({ key: 'pKsp', label: 'ZnO solubility product pK_sp (25 °C)', unit: '–', value: 16.7, min: 10, max: 20, group: 'Kinetics', help: '−log₁₀ of [Zn²⁺][OH⁻]² (mol/L units) at saturation. Growth stops when the ion product falls to this value.' }),
  ZI({ key: 'Ns', label: 'Nucleation-site density on the solid', unit: '1/m²', value: 5e13, min: 1e10, max: 1e16, log: true, bus: 'geom.number_density', group: 'Nucleation and growth', help: 'Available heterogeneous sites per m² of pore surface (seed-layer quality); 10¹³–10¹⁴ 1/m² is typical of seeded ZnO arrays. More sites give denser but shorter, thinner rods.' }),
  ZI({ key: 'site_cv', label: 'Site-density variation', unit: '%', value: 15, min: 0, max: 80, group: 'Nucleation and growth', help: 'Coefficient of variation of the site density across the bed (smooth, seeded log-normal field). It is a main source of length non-uniformity.' }),
  ZI({ key: 'gamma', label: 'ZnO–solution interfacial energy', unit: 'J/m²', value: 0.3, min: 0.02, max: 2, log: true, group: 'Nucleation and growth', help: 'Interfacial energy in the classical nucleation barrier. It matters only near saturation; at the high supersaturation of the reference case the barrier is small.' }),
  ZI({ key: 'theta', label: 'Nucleus contact angle on the solid', unit: '°', value: 60, min: 5, max: 180, group: 'Nucleation and growth', help: 'Wetting angle of the nucleus on the pore surface; it scales the heterogeneous barrier by f(θ). 180° recovers homogeneous nucleation.' }),
  ZI({ key: 'w_nuc', label: 'Site activation frequency (25 °C)', unit: '1/s', value: 0.01, min: 1e-6, max: 100, log: true, group: 'Nucleation and growth', help: 'Rate at which a free site nucleates when the barrier is negligible, at 25 °C. Calibration parameter controlling the induction time.' }),
  ZI({ key: 'd_seed', label: 'Nucleus (seed) size', unit: 'nm', value: 10, min: 1, max: 200, group: 'Nucleation and growth', help: 'Diameter and height of a rod at birth. Rods grow from this size by the anisotropic volume balance.' }),
  ZI({ key: 'beta_aniso', label: 'Axial / lateral growth-rate ratio', unit: '–', value: 10, min: 1, max: 100, log: true, group: 'Nucleation and growth', help: 'Ratio of c-axis growth velocity to lateral facet velocity. It fixes the aspect ratio that rods approach; capping agents and pH change it in practice.' }),
  ZI({ key: 'E', label: 'Applied electric field', unit: 'V/m', value: 1000, min: 0, max: 5000, bus: 'proc.E_field', group: 'Electric field', help: 'Mean field across the bed thickness (voltage / thickness). Drives Joule heating (∝ E²) and electromigration of the ions.' }),
  ZI({ key: 'field_nu', label: 'Electrode field non-uniformity', unit: '%', value: 5, min: 0, max: 50, group: 'Electric field', help: 'Parabolic drop of the applied field from the centre to the corners of the electrode. Zero gives an ideal parallel-plate field.' }),
  ZI({ key: 'spectators', label: 'Include spectator ions in conductivity', type: 'bool', value: true, group: 'Electric field', help: 'Adds the salt anion and base cation to the Nernst–Einstein conductivity. Switch off to use Zn²⁺ and OH⁻ only, as in the bare reference description.' }),
  ZI({ key: 'Emod', label: 'Drained Young’s modulus of the bed', unit: 'GPa', value: 1.5, min: 0.001, max: 400, log: true, group: 'Mechanics', help: 'Elastic modulus of the porous skeleton. Stress from constrained thermal expansion scales with it.' }),
  ZI({ key: 'nu', label: 'Poisson’s ratio', unit: '–', value: 0.3, min: -0.5, max: 0.49, bus: 'eff.nu', group: 'Mechanics', help: 'Drained Poisson’s ratio of the skeleton (may be negative for auxetic scaffolds).' }),
  ZI({ key: 'biot', label: 'Biot coefficient', unit: '–', value: 0.6, min: 0.05, max: 1, group: 'Mechanics', help: 'Coupling between pore pressure and skeleton stress; must lie between the porosity and 1.' }),
  ZI({ key: 'cte', label: 'Linear thermal expansion of the skeleton', unit: '1/K', value: 1e-5, min: -2e-5, max: 1e-4, bus: 'eff.cte', group: 'Mechanics', help: 'Drained linear expansion coefficient; negative values (NTE scaffolds) reverse the sign of the thermal stress.' }),
  ZI({ key: 'eta_g', label: 'Growth eigenstrain coupling', unit: '–', value: 0.1, min: 0, max: 1, group: 'Mechanics', help: 'Fraction of the deposited ZnO volume that strains the skeleton (0 = rods grow freely in the pores, 1 = fully constrained).' }),
  ZI({ key: 'ctrl', label: 'Feedback control', type: 'select', value: 'none', group: 'Control', options: [{ value: 'none', label: 'Open loop' }, { value: 'temp', label: 'PID on field → temperature set-point' }, { value: 'growth', label: 'PI dosing → growth-rate set-point' }, { value: 'both', label: 'Both loops' }], help: 'Feedback-optimised synthesis: modulate the field to hold temperature and/or the bottom concentrations to hold growth rate.' }),
  ZI({ key: 'T_set', label: 'Temperature set-point', unit: 'K', value: 318, min: 280, max: 470, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['temp', 'both'] }, help: 'Target bed-average temperature for the PID loop.' }),
  ZI({ key: 'Kp', label: 'PID proportional gain', unit: 'V/m/K', value: 200, min: 0, max: 5000, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['temp', 'both'] }, help: 'Field change per kelvin of error. Higher gain responds faster but saturates the actuator.' }),
  ZI({ key: 'Ki', label: 'PID integral gain', unit: 'V/m/K/s', value: 2, min: 0, max: 100, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['temp', 'both'] }, help: 'Removes steady-state error; anti-windup stops integration while the field is saturated.' }),
  ZI({ key: 'Kd', label: 'PID derivative gain', unit: 'V·s/m/K', value: 0, min: 0, max: 1e5, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['temp', 'both'] }, help: 'Acts on the measured temperature rate to damp overshoot.' }),
  ZI({ key: 'E_max', label: 'Maximum field (actuator limit)', unit: 'V/m', value: 3000, min: 100, max: 10000, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['temp', 'both'] }, help: 'Upper limit of the controlled field (power-supply limit).' }),
  ZI({ key: 'G_set', label: 'Growth-rate set-point', unit: 'nm/min', value: 5, min: 0.01, max: 1000, log: true, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['growth', 'both'] }, help: 'Target mean axial growth rate for the dosing loop.' }),
  ZI({ key: 'Kg', label: 'Dosing proportional gain', unit: '–', value: 1, min: 0, max: 20, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['growth', 'both'] }, help: 'Change of the bottom-concentration multiplier per unit relative growth-rate error (integral time fixed at 300 s).' }),
  ZI({ key: 'dose_max', label: 'Maximum dosing multiplier', unit: '×', value: 5, min: 1, max: 50, group: 'Control', showIf: { mode: ['zno_reactor'], ctrl: ['growth', 'both'] }, help: 'Upper limit of the boundary concentration relative to the nominal value (feed-tank strength).' }),
  ZI({ key: 'nx', label: 'Cells in x', value: 20, min: 4, max: 40, step: 1, sens: false, group: 'Numerics', help: 'Grid cells along x. The reference description quotes 100 × 100 × 10 nodes; the in-browser default is coarser and the sensitivity indicator reports the effect.' }),
  ZI({ key: 'ny', label: 'Cells in y', value: 20, min: 4, max: 40, step: 1, sens: false, group: 'Numerics', help: 'Grid cells along y.' }),
  ZI({ key: 'nz', label: 'Cells in z', value: 10, min: 4, max: 20, step: 1, sens: false, group: 'Numerics', help: 'Grid cells through the thickness, the direction with the steep gradients.' }),
  ZI({ key: 'dt', label: 'Time step', unit: 's', value: 40, min: 0.05, max: 600, log: true, sens: false, group: 'Numerics', help: 'Backward-Euler step. The reference nominal value is 0.1 s; larger steps are stable but less accurate, and the step is enlarged automatically if the work budget is exceeded.' }),
  ZI({ key: 't_total', label: 'Total simulated time', unit: 's', value: 1000, min: 10, max: 172800, log: true, bus: 'proc.time', group: 'Numerics', help: 'Process duration.' }),
  ZI({ key: 'max_it', label: 'Picard iteration cap per step', value: 10, min: 2, max: 12, step: 1, sens: false, group: 'Numerics', help: 'Maximum coupling iterations per time step (tolerance 10⁻⁴ on the relative update of concentrations and temperature).' }),
  ZI({ key: 'rod_cv', label: 'Rod-to-rod size scatter (display)', unit: '%', value: 15, min: 0, max: 60, sens: false, group: 'Display', help: 'Log-normal scatter applied when sampling individual rods for the 3-D view; it does not affect the simulation.' }),
  ZI({ key: 'tilt_sd', label: 'Rod tilt spread (display)', unit: '°', value: 5, min: 0, max: 45, sens: false, group: 'Display', help: 'Standard deviation of the rod tilt in the 3-D view; alignment is not predicted by the model.' }),
  ZI({ key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Display', help: 'Seed for the site-density field and the sampled rod array; identical seeds give identical results.' }),
];

// ================================================================================================
// 2. CVD / PECVD — Grove model with boundary-layer mass transfer, depletion and step coverage
// ================================================================================================
const GASES = { H2: { mu0: 8.9e-6, M: 2.016e-3 }, N2: { mu0: 1.78e-5, M: 28.01e-3 } };
// k0/Ea: typical-order first-order surface kinetics (illustrative calibrations, see warnings); D0: binary diffusivity at 300 K, 1 atm
const CVDC = {
  sih4: { name: 'SiH₄ → Si + 2 H₂ (polysilicon)', r: ['SiH4'], p: ['Si', 'H2'], film: 'Si', rho: 2329, k0: 5.9e6, Ea: 160e3, D0: 6.0e-5, gas: 'H2' },
  ch4: { name: 'CH₄ → C (graphene) + 2 H₂ on Cu', r: ['CH4'], p: ['C', 'H2'], film: 'C', rho: 2267, k0: 1.8e4, Ea: 200e3, D0: 7.0e-5, gas: 'H2', mono: 0.335e-9 },
  ticl4: { name: 'TiCl₄ + O₂ → TiO₂ + 2 Cl₂', r: ['TiCl4', 'O2'], p: ['TiO2', 'Cl2'], film: 'TiO2', rho: 4230, k0: 15, Ea: 88e3, D0: 6.0e-6, gas: 'N2' },
  dez: { name: 'Zn(C₂H₅)₂ + H₂O → ZnO + 2 C₂H₆ (MOCVD)', r: ['Zn(C2H5)2', 'H2O'], p: ['ZnO', 'C2H6'], film: 'ZnO', rho: 5606, k0: 63, Ea: 40e3, D0: 8.0e-6, gas: 'N2' },
  wf6: { name: 'WF₆ + 3 H₂ → W + 6 HF', r: ['WF6', 'H2'], p: ['W', 'HF'], film: 'W', rho: 19300, k0: 1.1e3, Ea: 69e3, D0: 3.0e-5, gas: 'H2' },
};
function runCVD(inp) {
  const c = CVDC[inp.cvd_chem] || CVDC.sih4, gas = GASES[c.gas], pe = inp.cvd_kind === 'pecvd', warnings = [];
  const Tk = inp.cvd_T + 273.15, P = inp.cvd_P, U = inp.cvd_U, Ls = inp.cvd_L, Hc = inp.cvd_H, mp = MM(c.r[0]) / NA, Nf = c.rho * NA / MM(c.film);
  const at = T => { const D = c.D0 * (T / 300) ** 1.75 * 101325 / P, mu = gas.mu0 * (T / 300) ** 0.7, rg = P * gas.M / (R * T), vbar = Math.sqrt(8 * KB * T / (PI * mp)), ksT = Math.min(c.k0 * Math.exp(-c.Ea / (R * T)), vbar / 4), ks = pe ? inp.cvd_s * vbar / 4 : ksT, Re = rg * U * Ls / mu, Sc = mu / (rg * D), hfd = 4.86 * D / (2 * Hc), hg = Math.max(0.664 * D / Ls * Math.sqrt(Re) * Sc ** (1 / 3), hfd); return { D, mu, rg, vbar, ks, Re, Sc, hfd, hg, C0: inp.cvd_x * P / (KB * T) }; };
  const a = at(Tk), ns = 120, dxs = Ls / ns, xs = [], rate = []; let C = a.C0;
  for (let i = 0; i < ns; i++) { const x = (i + 0.5) * dxs, hx = Math.max(0.332 * a.D / x * Math.sqrt(a.rg * U * x / a.mu) * a.Sc ** (1 / 3), a.hfd), K = a.ks * hx / (a.ks + hx), Cn = C * Math.exp(-K * dxs / (U * Hc)); xs.push(x * 1e3); rate.push(U * Hc * (C - Cn) / (dxs * Nf)); C = Cn; }
  const rMean = N.mean(rate), rMax = N.maxOf(rate), rMin = N.minOf(rate), unif = rMean > 0 ? (rMax - rMin) / (2 * rMean) * 100 : 0, thick = rMean * inp.cvd_t, util = 1 - C / a.C0, Da = a.ks / a.hg;
  const regime = Da > 10 ? 'mass-transport limited' : Da < 0.1 ? 'surface-reaction limited' : 'mixed control', stick = Math.min(1, 4 * a.ks / a.vbar), sc = 1 / Math.cosh(inp.cvd_AR * Math.sqrt(3 * stick));
  const lam = KB * Tk / (Math.SQRT2 * PI * (3.5e-10) ** 2 * P), Kn = lam / Hc;
  let Tx = NaN; if (!pe) { const g = T => { const q = at(T); return Math.log(q.ks / q.hg); }; if (g(350) * g(2000) < 0) Tx = N.brent(g, 350, 2000, 1e-6); }
  const Ts = N.linspace(Math.max(100, inp.cvd_T - 300), inp.cvd_T + 300, 80), rT = [], rS = [], rM = [];
  for (const Tc of Ts) { const q = at(Tc + 273.15); rT.push(q.ks * q.hg / (q.ks + q.hg) * q.C0 / Nf * 6e10); rS.push(q.ks * q.C0 / Nf * 6e10); rM.push(q.hg * q.C0 / Nf * 6e10); }
  const ARs = N.linspace(0, Math.max(20, inp.cvd_AR * 1.5), 60);
  if (pe) warnings.push('PECVD mode replaces the Arrhenius surface step by a radical flux with the entered sticking probability; plasma power, sheath and ion effects are not modelled.');
  warnings.push('Kinetic constants are typical-order first-order calibrations for each chemistry (true kinetics are often non-first-order, e.g. WF₆/H₂ is zero order in WF₆ and half order in H₂); fit k₀ and E_a to measured rates before quantitative use.');
  if (c.mono) warnings.push('Graphene on Cu is self-limiting at about one monolayer; the "thickness" is meaningful only up to 0.335 nm and the Grove rate gives the coverage time.');
  if (Kn > 0.1) warnings.push(`Knudsen number ${f(Kn)} in the channel: the continuum boundary-layer correlation is not valid; treat mass-transfer numbers as rough.`);
  if (a.Re > 2000) warnings.push(`Reynolds number ${f(a.Re)} is beyond the laminar range assumed by the boundary-layer correlation.`);
  const mlt = c.mono ? c.mono / Math.max(rMean, 1e-300) : 0;
  const kpis = [
    kpi('rate', 'Mean deposition rate', rMean * 6e10, 'nm/min', `Grove-model growth rate averaged along the ${f(Ls * 1e3)} mm susceptor: the series combination of gas-phase mass transfer and surface reaction divided by the film number density.`),
    kpi('thickness', 'Film thickness', thick * 1e9, 'nm', `Mean rate × ${f(inp.cvd_t)} s.` + (c.mono ? ` One graphene monolayer (0.335 nm) is reached after ${f(mlt)} s.` : ' Compare with the target; thickness scales linearly with time at fixed conditions.')),
    kpi('regime', 'Rate-controlling regime', regime, '', `Damköhler number k_s/h_g = ${f(Da)}: above 10 the gas supply limits growth (weak temperature dependence, strong flow dependence), below 0.1 the surface reaction limits it (Arrhenius, uniform).`),
    kpi('ks', 'Surface rate constant k_s', a.ks, 'm/s', pe ? 'Radical surface loss velocity s·v̄/4 from kinetic theory with the entered sticking probability.' : `Arrhenius constant k₀·exp(−E_a/RT) with E_a = ${f(c.Ea / 1e3)} kJ/mol, capped at the kinetic-theory impingement limit v̄/4.`),
    kpi('hg', 'Mass-transfer coefficient h_g', a.hg, 'm/s', `Plate-average laminar boundary-layer value 0.664·(D/L)·Re^½·Sc^⅓ (Re = ${f(a.Re)}, Sc = ${f(a.Sc)}), floored at the fully developed channel value Sh = 4.86.`),
    kpi('uniformity', 'Thickness non-uniformity along flow', unif, '± %', 'Half of (max − min)/mean of the local rate along the susceptor, caused by boundary-layer growth and precursor depletion. Tilted susceptors or rotation reduce it.', unif < 5),
    kpi('utilisation', 'Precursor utilisation', util * 100, '%', 'Fraction of the inlet precursor deposited on the susceptor (plug-flow depletion). The rest leaves in the exhaust and must be abated.', util > 0.2),
    kpi('step', 'Step coverage in trench', sc * 100, '%', `Bottom-to-top thickness ratio for aspect ratio ${f(inp.cvd_AR)} from the Knudsen diffusion–reaction (Thiele) model with reactive sticking probability ${f(stick)}.`, sc > 0.8),
    kpi('Tcross', 'Regime crossover temperature', Number.isFinite(Tx) ? Tx - 273.15 : 0, '°C', Number.isFinite(Tx) ? 'Temperature at which k_s = h_g; operate well below it for uniform, conformal films and above it for maximum rate.' : 'No crossover between 77 and 1727 °C for this chemistry and flow (or plasma mode): the regime does not change in the practical range.'),
  ];
  return finish({
    summary: `${c.name} at ${f(inp.cvd_T)} °C and ${f(P)} Pa grows at ${f(rMean * 6e10)} nm/min (${regime}, k_s/h_g = ${f(Da)}), giving ${f(thick * 1e9)} nm in ${f(inp.cvd_t)} s with ±${f(unif)} % variation along the susceptor. ${f(util * 100)} % of the precursor is used and the step coverage in an aspect-ratio-${f(inp.cvd_AR)} trench is ${f(sc * 100)} %.`,
    kpis, warnings,
    series: [
      ser('rateT', 'Deposition rate versus temperature', 'temperature (°C)', 'rate (nm/min)', [{ name: 'Grove (combined)', x: Ts, y: rT }, { name: 'surface-reaction limit', x: Ts, y: rS }, { name: 'mass-transport limit', x: Ts, y: rM }], 'The combined rate follows the Arrhenius line at low temperature and saturates at the weakly temperature-dependent mass-transport limit. The knee is the regime change.', { logy: true }),
      ser('rateX', 'Local rate along the susceptor', 'position (mm)', 'rate (nm/min)', [{ name: 'local rate', x: xs, y: rate.map(v => v * 6e10) }], 'Thickness profile in the flow direction: decline from the leading edge reflects boundary-layer thickening and precursor depletion.'),
      ser('step', 'Step coverage versus aspect ratio', 'aspect ratio (depth/width)', 'bottom/top thickness (%)', [{ name: 'step coverage', x: ARs, y: ARs.map(v => 100 / Math.cosh(v * Math.sqrt(3 * stick))) }], 'Conformality falls with aspect ratio at a rate set by the sticking probability; low sticking (reaction-limited) gives conformal films.'),
    ],
    equations: [
      eqn('chemical', 'Overall deposition reaction', safeEq(c.r, c.p)),
      eqn('governing', 'Grove model', 'v = (k_s h_g /(k_s + h_g)) · C_g / N_film,   C_g = x P/(k_B T)'),
      eqn('constitutive', 'Surface kinetics and mass transfer', 'k_s = k₀ exp(−E_a/RT);   h_g(x) = 0.332 (D/x) Re_x^½ Sc^⅓;   D ∝ T^1.75/P'),
      eqn('governing', 'Plug-flow precursor depletion', 'U H dC/dx = −K(x) C,   K = k_s h_g/(k_s + h_g)'),
      eqn('formula', 'Step coverage (Knudsen diffusion–reaction)', 'SC = 1/cosh(AR·√(3s)),   s = 4k_s/v̄'),
    ],
    checks: [chk('Precursor balance (deposited = inlet − outlet)', true, 'exact by construction of the segment-wise integration'), chk('Laminar flow', a.Re < 2000, `Re = ${f(a.Re)}`), chk('Continuum gas transport', Kn < 0.1, `Kn = ${f(Kn)} (mean free path ${f(lam)} m)`), chk('Sticking probability ≤ 1', stick <= 1, `s = ${f(stick)}`)],
    recommendations: [
      rec('high', Da > 10 ? `Growth is transport-limited: lower the temperature towards ${Number.isFinite(Tx) ? f(Tx - 273.15) + ' °C' : 'the kinetic regime'} or reduce pressure to improve uniformity and step coverage.` : Da < 0.1 ? 'Growth is reaction-limited: temperature uniformity across the susceptor (±1 K) controls thickness uniformity; flow details matter little.' : 'Mixed control: both temperature and flow uniformity matter; move decisively into one regime for a robust process.', `Damköhler number ${f(Da)}.`),
      rec('medium', `Along-flow variation is ±${f(unif)} %; tilt the susceptor 3–10°, raise the gas velocity or rotate the wafer to compensate depletion.`, 'Local-rate profile from the depletion integration.'),
      rec('medium', sc < 0.8 ? `Step coverage is ${f(sc * 100)} %; lower the sticking probability (lower temperature, reaction-limited regime) or switch to ALD for aspect ratios above ${f(inp.cvd_AR)}.` : `Step coverage of ${f(sc * 100)} % is adequate for aspect ratio ${f(inp.cvd_AR)}.`, 'Thiele-modulus conformality model.'),
      rec('medium', `Only ${f(util * 100)} % of the precursor is deposited; recirculate or abate the exhaust (scrubber/burn box) and reduce flow excess to cut precursor waste and hazardous by-products (${c.p.slice(1).map(chem.prettyFormula).join(', ')}).`, 'Precursor utilisation from the plug-flow balance.', true),
    ],
    publish: { 'growth.rate': rMean, 'fab.dep_rate': rMean, 'geom.thickness': thick, 'proc.T_C': inp.cvd_T, 'proc.P': P, 'fab.yield': N.clamp(util, 0, 1), 'proc.time': inp.cvd_t, 'material.formula': c.film },
    exports: [{ name: 'cvd_recipe.json', mime: 'application/json', content: JSON.stringify({ process: pe ? 'PECVD' : 'thermal CVD', chemistry: c.name, T_C: inp.cvd_T, P_Pa: P, precursor_mole_fraction: inp.cvd_x, velocity_m_s: U, time_s: inp.cvd_t, predicted: { rate_nm_min: rMean * 6e10, thickness_nm: thick * 1e9, nonuniformity_pct: unif, regime, step_coverage: sc } }, null, 2), description: 'CVD recipe with predicted rate, thickness, uniformity and regime.' }],
    numerics: { method: 'Grove series-resistance model; segment-wise exponential integration of plug-flow depletion (120 segments); Brent root for the regime crossover', dof: ns, steps: ns, iterations: 1, residual: 0, notes: 'Analytical correlations; no CFD.' },
  });
}

// ================================================================================================
// 3. ALD
// ================================================================================================
const ALDC = {
  tma: { name: 'TMA / H₂O → Al₂O₃', prec: 'Al(CH3)3', film: 'Al2O3', nM: 2, rho: 3000, gam: 4.0e18, s0: 1e-3, win: [100, 300], p: ['Al2O3', 'CH4'], hA: '‖–OH + Al(CH₃)₃(g) → ‖–O–Al(CH₃)₂ + CH₄(g)', hB: '‖–O–Al(CH₃)₂ + 2 H₂O(g) → ‖–O–Al(OH)₂ + 2 CH₄(g)' },
  dez: { name: 'DEZ / H₂O → ZnO', prec: 'Zn(C2H5)2', film: 'ZnO', nM: 1, rho: 5606, gam: 8.0e18, s0: 7e-3, win: [100, 200], p: ['ZnO', 'C2H6'], hA: '‖–OH + Zn(C₂H₅)₂(g) → ‖–O–Zn(C₂H₅) + C₂H₆(g)', hB: '‖–O–Zn(C₂H₅) + H₂O(g) → ‖–O–Zn–OH + C₂H₆(g)' },
  ticl4: { name: 'TiCl₄ / H₂O → TiO₂', prec: 'TiCl4', film: 'TiO2', nM: 1, rho: 3900, gam: 1.5e18, s0: 1e-3, win: [100, 400], p: ['TiO2', 'HCl'], hA: '‖–OH + TiCl₄(g) → ‖–O–TiCl₃ + HCl(g)', hB: '‖–O–TiCl₃ + 3 H₂O(g) → ‖–O–Ti(OH)₃ + 3 HCl(g)' },
  tdmah: { name: 'TDMAH / H₂O → HfO₂', prec: 'Hf(N(CH3)2)4', film: 'HfO2', nM: 1, rho: 9680, gam: 2.77e18, s0: 0.05, win: [150, 300], p: ['HfO2', 'HN(CH3)2'], hA: '‖–OH + Hf(N(CH₃)₂)₄(g) → ‖–O–Hf(N(CH₃)₂)₃ + HN(CH₃)₂(g)', hB: '‖–O–Hf(N(CH₃)₂)₃ + 3 H₂O(g) → ‖–O–Hf(OH)₃ + 3 HN(CH₃)₂(g)' },
};
function runALD(inp) {
  const c = ALDC[inp.ald_chem] || ALDC.tma, warnings = [], Tk = inp.ald_T + 273.15, P = inp.ald_P, tp = inp.ald_tp, m = MM(c.prec) / NA;
  const gpcSat = c.gam * MM(c.film) / (c.nM * c.rho * NA), flux = P / Math.sqrt(2 * PI * m * KB * Tk), theta = 1 - Math.exp(-c.s0 * flux * tp / c.gam), gpc = gpcSat * theta;
  const cycles = Math.ceil(inp.ald_thk * 1e-9 / Math.max(gpc, 1e-14)), tSat = 4.605 * c.gam / (c.s0 * flux), Kg = c.gam * Math.sqrt(2 * PI * m * KB * Tk), AR = inp.ald_AR;
  const need = Kg * (1 + 4.75 * AR + 1.5 * AR * AR), Pt = P * tp, disc = 4.75 * 4.75 - 6 * (1 - Pt / Kg), aMax = Pt > Kg ? (-4.75 + Math.sqrt(disc)) / 3 : 0, conf = N.clamp(aMax / Math.max(AR, 1e-9), 0, 1), tConf = need / P;
  const vbar = Math.sqrt(8 * KB * Tk / (PI * m)), dh = inp.ald_dhole * 1e-9, Lh = AR * dh, DK = dh * vbar / 3, tau = inp.ald_V / inp.ald_S, tPurge = tau * Math.log(1e4) + 4 * Lh * Lh / (PI * PI * DK) * Math.log(1e4);
  const tpUse = Math.max(tp, 0), tCyc = 2 * (tpUse + tPurge), tTot = cycles * tCyc, area = inp.ald_area * 1e-4, dosed = P * (inp.ald_V * 1e-3 + inp.ald_S * 1e-3 * tp) / (KB * Tk), util = N.clamp(c.gam * theta * area / Math.max(dosed, 1e-300), 0, 1);
  const lam = KB * Tk / (Math.SQRT2 * PI * (5e-10) ** 2 * Math.max(P, 1e-6)), massPrec = cycles * dosed * m;
  if (inp.ald_T < c.win[0]) warnings.push(`${f(inp.ald_T)} °C is below the usual ALD window (${c.win[0]}–${c.win[1]} °C) for this chemistry: condensation or incomplete ligand exchange is likely; the model does not include either.`);
  if (inp.ald_T > c.win[1]) warnings.push(`${f(inp.ald_T)} °C is above the usual ALD window (${c.win[0]}–${c.win[1]} °C): precursor decomposition adds a CVD component and growth is no longer self-limiting.`);
  if (theta < 0.95) warnings.push(`The pulse gives only ${f(theta * 100)} % of saturation on a flat surface; growth is dose-dependent (not true ALD) until the pulse exceeds about ${f(tSat)} s.`);
  if (lam < dh) warnings.push('Mean free path is shorter than the hole diameter: transport in the hole is not purely molecular flow, so the Gordon model is optimistic.');
  warnings.push('Saturated growth per cycle uses a typical ligand-limited metal uptake per cycle and a typical initial sticking probability; both depend on reactor, temperature and surface termination and need calibration (ellipsometry/QCM).');
  const tps = lsp(-3, 1.5, 70), ARs = N.linspace(1, Math.max(100, AR * 1.5), 70);
  let overall = ''; overall = safeEq([c.prec, 'H2O'], c.p);
  return finish({
    summary: `${c.name} at ${f(inp.ald_T)} °C grows ${f(gpc * 1e10)} Å per cycle (${f(theta * 100)} % of the saturated ${f(gpcSat * 1e10)} Å), so ${cycles} cycles (${f(tTot / 60)} min including ${f(tPurge)} s purges) give ${f(inp.ald_thk)} nm. With a ${f(tp)} s pulse at ${f(P)} Pa the precursor coats ${f(conf * 100)} % of the depth of aspect-ratio-${f(AR)} holes; full conformality needs ${f(tConf)} s.`,
    kpis: [
      kpi('gpc', 'Growth per cycle', gpc * 1e9, 'nm/cycle', `Saturated value ${f(gpcSat * 1e9)} nm (metal uptake ${f(c.gam / 1e18)} nm⁻² per cycle × film molar volume) multiplied by the Langmuir saturation fraction reached in the pulse.`),
      kpi('cycles', 'Cycles for target thickness', cycles, 'cycles', `Target ${f(inp.ald_thk)} nm divided by the growth per cycle, rounded up; thickness is digital in steps of one GPC.`),
      kpi('sat', 'Saturation fraction on flat surface', theta * 100, '%', `1 − exp(−s₀Φt/σ) with impingement flux Φ = ${f(flux)} m⁻²s⁻¹ and initial sticking probability ${f(c.s0)}; 99 % saturation needs ${f(tSat)} s.`, theta > 0.95),
      kpi('exposure', 'Exposure for full conformality', need, 'Pa·s', `Gordon model: σ√(2πmk_BT)(1 + 19a/4 + 3a²/2) for aspect ratio a = ${f(AR)}; the present pulse delivers ${f(Pt)} Pa·s.`),
      kpi('conformal', 'Coated depth fraction in holes', conf * 100, '%', `The present exposure saturates the walls down to aspect ratio ${f(aMax)}; below that depth the film is absent or thinner.`, conf >= 0.999),
      kpi('purge', 'Required purge time', tPurge, 's', `Chamber dilution to 10⁻⁴ (residence time ${f(tau)} s) plus Knudsen out-diffusion from the holes (D_K = ${f(DK)} m²/s). Shorter purges risk CVD-type parasitic growth.`),
      kpi('time', 'Total process time', tTot / 60, 'min', 'Cycles × (two pulses + two purges). Purge time usually dominates and sets throughput.'),
      kpi('util', 'Precursor utilisation', util * 100, '%', 'Molecules chemisorbed on the stated surface area divided by molecules dosed per pulse (chamber fill plus pumped-through gas).', util > 0.1),
    ],
    series: [
      ser('sat', 'Saturation curve', 'pulse time (s)', 'growth per cycle (nm)', [{ name: 'GPC', x: tps, y: tps.map(t => gpcSat * 1e9 * (1 - Math.exp(-c.s0 * flux * t / c.gam))) }], 'Self-limiting behaviour: growth per cycle rises with dose and saturates. Operate on the plateau with a safety margin.', { logx: true }),
      ser('gordon', 'Exposure needed versus aspect ratio', 'aspect ratio', 'exposure (Pa·s)', [{ name: 'required exposure', x: ARs, y: ARs.map(a => Kg * (1 + 4.75 * a + 1.5 * a * a)) }, { name: 'present pulse', x: [ARs[0], ARs[ARs.length - 1]], y: [Math.max(Pt, 1e-12), Math.max(Pt, 1e-12)] }], 'Required exposure grows with the square of aspect ratio. Holes are fully coated where the present-pulse line lies above the curve.', { logy: true }),
    ],
    equations: [eqn('chemical', 'Overall reaction', overall), eqn('chemical', 'Half-reaction A (metal precursor pulse)', c.hA, '‖ denotes a surface site.'), eqn('chemical', 'Half-reaction B (water pulse)', c.hB),
      eqn('formula', 'Growth per cycle', 'GPC = θ · σ_M M_film /(n_M ρ N_A),   θ = 1 − exp(−s₀ Φ t/σ_M),   Φ = P/√(2π m k_B T)'),
      eqn('formula', 'Gordon conformality exposure', '(P·t) = σ_M √(2π m k_B T) · (1 + 19a/4 + 3a²/2)'), eqn('formula', 'Purge time', 't_purge = (V/S) ln 10⁴ + (4L²/π²D_K) ln 10⁴,   D_K = d v̄/3')],
    checks: [chk('Self-limiting dose reached', theta > 0.95, `${f(theta * 100)} % of saturation`), chk('Temperature inside ALD window', inp.ald_T >= c.win[0] && inp.ald_T <= c.win[1], `${c.win[0]}–${c.win[1]} °C typical`), chk('Molecular flow in holes', lam > dh, `mean free path ${f(lam)} m vs hole diameter ${f(dh)} m`), chk('Holes fully coated', conf >= 0.999, `coated to aspect ratio ${f(aMax)} of ${f(AR)}`)],
    warnings,
    recommendations: [
      rec('high', conf < 0.999 ? `Lengthen the precursor exposure to at least ${f(tConf)} s at ${f(P)} Pa (or use stop-flow/soak mode) to coat aspect ratio ${f(AR)} conformally.` : `The exposure already saturates aspect ratio ${f(AR)}; a pulse of ${f(tConf)} s is the minimum, so the present ${f(tp)} s can be shortened to save precursor if only this structure matters.`, 'Gordon exposure model.'),
      rec('medium', `Use purges of at least ${f(tPurge)} s; verify by checking that thickness does not increase when the purge is doubled.`, 'Chamber residence time and Knudsen out-diffusion estimate.'),
      rec('medium', `Confirm the ${f(gpc * 1e10)} Å/cycle with in-situ QCM or ellipsometry on a witness wafer over 50–200 cycles; nucleation delay on H-terminated or noble surfaces is not included.`, 'Model uses a typical saturated metal uptake.'),
      rec('medium', `About ${f(massPrec * 1e3)} g of ${chem.prettyFormula(c.prec)} is dosed for this run at ${f(util * 100)} % utilisation; shorten over-saturated pulses, batch more substrates per run and trap unreacted precursor (it is pyrophoric or corrosive) before the pump.`, 'Dose versus chemisorbed amount.', true),
    ],
    publish: { 'geom.thickness': cycles * gpc, 'fab.dep_rate': gpc / tCyc, 'growth.rate': gpc / tCyc, 'proc.T_C': inp.ald_T, 'proc.time': tTot, 'fab.yield': util, 'fab.gpc_m': gpc, 'material.formula': c.film },
    exports: [{ name: 'ald_recipe.json', mime: 'application/json', content: JSON.stringify({ process: 'thermal ALD', chemistry: c.name, T_C: inp.ald_T, precursor_partial_pressure_Pa: P, pulse_s: tp, purge_s: tPurge, cycles, predicted: { gpc_nm: gpc * 1e9, thickness_nm: cycles * gpc * 1e9, total_min: tTot / 60, coated_depth_fraction: conf } }, null, 2), description: 'ALD recipe: pulse/purge times, cycle count and predicted thickness.' }],
    numerics: { method: 'Closed-form Langmuir saturation kinetics and Gordon molecular-flow exposure model', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'No reactor flow simulation.' },
  });
}

// ================================================================================================
// 4. PVD — sputtering / evaporation
// ================================================================================================
// rho kg/m³, Y500: Ar⁺ sputter yield at 500 eV (handbook-typical), Tb K, Hv kJ/mol, E GPa, nu, cte 1/K
const PVDM = {
  Cu: [8960, 2.35, 2835, 300, 120, 0.34, 16.5e-6], Al: [2700, 1.05, 2743, 284, 70, 0.35, 23.1e-6], Ti: [4506, 0.51, 3560, 425, 116, 0.32, 8.6e-6], Au: [19300, 2.4, 3243, 324, 79, 0.44, 14.2e-6], Ag: [10490, 3.12, 2435, 255, 83, 0.37, 18.9e-6],
  W: [19300, 0.57, 6203, 774, 411, 0.28, 4.5e-6], Cr: [7190, 1.18, 2944, 340, 279, 0.21, 4.9e-6], Ni: [8908, 1.45, 3003, 378, 200, 0.31, 13.4e-6], Pt: [21450, 1.4, 4098, 510, 168, 0.38, 8.8e-6], Si: [2329, 0.5, 3538, 359, 160, 0.22, 2.6e-6],
};
function runPVD(inp) {
  const el = PVDM[inp.pvd_mat] ? inp.pvd_mat : 'Cu', [rho, Y500, Tb, Hv, Ef, nuf, af] = PVDM[el], M = MM(el), m = M / NA, Om = M / (rho * NA), evap = inp.pvd_proc === 'evap', warnings = [];
  const h = inp.pvd_h * 1e-2, a = evap ? 0 : inp.pvd_a * 1e-2, Rw = inp.pvd_wafer * 0.5e-3, t = inp.pvd_t;
  let Nsrc, Y = 0, I = 0, Pv = 0, powerW;
  if (evap) { const Ts = inp.pvd_Tsrc + 273.15; Pv = 101325 * Math.exp(-Hv * 1e3 / R * (1 / Ts - 1 / Tb)); Nsrc = Pv / Math.sqrt(2 * PI * m * KB * Ts) * inp.pvd_Asrc * 1e-4; powerW = N.SIGMA_SB * 0.3 * inp.pvd_Asrc * 1e-4 * Ts ** 4 + Nsrc * Hv * 1e3 / NA; }
  else { const Eion = inp.pvd_V; I = inp.pvd_power / Math.max(inp.pvd_V, 1); Y = Eion > 30 ? Y500 * Math.sqrt(Eion / 500) : 0; Nsrc = Y * I / (N.Q_E * 1.1); powerW = inp.pvd_power; }
  const flux = r => Nsrc * h * h * (h * h + r * r + a * a) / (PI * ((h * h - r * r + a * a) ** 2 + 4 * h * h * r * r) ** 1.5);
  const rs = N.linspace(0, Rw, 60), th = rs.map(r => flux(r) * Om * t), tC = th[0], tMax = N.maxOf(th), tMin = N.minOf(th), unif = tMax > 0 ? (tMax - tMin) / (tMax + tMin) * 100 : 0;
  let coll = 0; for (let i = 1; i < rs.length; i++) coll += PI * (rs[i] ** 2 - rs[i - 1] ** 2) * 0.5 * (flux(rs[i]) + flux(rs[i - 1])); const collF = Nsrc > 0 ? coll / Nsrc : 0;
  const lam = KB * 300 / (Math.SQRT2 * PI * (3.6e-10) ** 2 * inp.pvd_P), Kn = lam / h, ball = Math.exp(-h / lam);
  const sTh = Ef * 1e9 / (1 - nuf) * (af - 2.6e-6) * (inp.pvd_Tdep - 25), sTot = sTh + inp.pvd_sint * 1e6, ts = inp.pvd_tsub * 1e-6, curv = 6 * sTot * tC * (1 - 0.28) / (130e9 * ts * ts), bow = curv * Rw * Rw / 2;
  const nm = 41, map = new Float64Array(nm * nm); for (let j = 0; j < nm; j++) for (let i = 0; i < nm; i++) { const x = (2 * i / (nm - 1) - 1) * Rw, y = (2 * j / (nm - 1) - 1) * Rw; map[i + nm * j] = flux(Math.hypot(x, y)) * Om * t * 1e9; }
  const mass = coll * m * t, kWhkg = mass > 0 ? powerW * t / 3.6e6 / mass : 0;
  warnings.push(evap ? 'Evaporation flux from the Hertz–Knudsen equation with a Clausius–Clapeyron vapour pressure anchored at the boiling point (constant ΔH_vap): accurate to a factor of a few; use tabulated vapour pressures for process design.' : 'Sputter yield is a tabulated 500 eV Ar⁺ value scaled with √E (reasonable for 200–1000 eV); emission is taken as a cosine law from a thin ring. Rates are collisionless upper bounds — gas scattering and redeposition typically reduce them by 2–5×.');
  if (!evap && Kn < 1) warnings.push(`Knudsen number ${f(Kn)}: only ${f(ball * 100)} % of the sputtered atoms arrive without a gas collision, so the real rate and energy are lower and the profile flatter than the cosine-law prediction.`);
  if (!evap && inp.pvd_V <= 30) warnings.push('Target voltage is below the sputtering threshold; no deposition.');
  warnings.push('Film stress combines thermal mismatch against silicon with the entered intrinsic stress; intrinsic (growth) stress depends on pressure and bombardment and must be measured.');
  return finish({
    summary: `${evap ? 'Evaporation' : 'Sputtering'} of ${el} deposits ${f(tC * 1e9)} nm at the wafer centre in ${f(t)} s (${f(tC / t * 6e10)} nm/min) with ±${f(unif)} % thickness variation over the ${f(inp.pvd_wafer)} mm wafer at ${f(inp.pvd_h)} cm throw. ${f(collF * 100)} % of the emitted atoms land on the wafer; estimated film stress is ${f(sTot / 1e6)} MPa and the wafer bow ${f(bow * 1e6)} µm.`,
    kpis: [
      kpi('rate', 'Deposition rate at centre', tC / t * 6e10, 'nm/min', evap ? `Hertz–Knudsen source flux at ${f(inp.pvd_Tsrc)} °C (vapour pressure ${f(Pv)} Pa) projected onto the wafer by the cosine law.` : `Sputter yield ${f(Y)} atoms/ion at ${f(inp.pvd_V)} eV × ion current ${f(I / 1.1)} A, projected by the cosine law from a ring source of radius ${f(inp.pvd_a)} cm.`),
      kpi('thickness', 'Thickness at centre', tC * 1e9, 'nm', 'Centre rate × deposition time; the radial profile is plotted and mapped.'),
      kpi('uniformity', 'Thickness non-uniformity over wafer', unif, '± %', '(max − min)/(max + min) of the radial profile. Increasing the throw distance or matching racetrack radius to wafer size improves it.', unif < 5),
      kpi('collection', 'Material collected on wafer', collF * 100, '%', 'Fraction of emitted atoms landing on the wafer; the rest coats shields and walls.', collF > 0.1),
      kpi('Kn', 'Knudsen number (throw distance)', Kn, '–', `Mean free path ${f(lam * 100)} cm at ${f(inp.pvd_P)} Pa divided by the throw distance. Kn > 1: ballistic line-of-sight transport; Kn < 0.1: thermalised, diffusive transport.`),
      kpi('stress', 'Film stress (thermal + intrinsic)', sTot / 1e6, 'MPa', `Biaxial thermal mismatch E/(1−ν)·(α_film − α_Si)·ΔT = ${f(sTh / 1e6)} MPa for deposition at ${f(inp.pvd_Tdep)} °C plus the entered intrinsic stress; positive is tensile.`, Math.abs(sTot) < 1e9),
      kpi('bow', 'Wafer bow (Stoney)', bow * 1e6, 'µm', `Curvature from Stoney’s equation for a ${f(inp.pvd_tsub)} µm silicon wafer; sign follows the stress.`),
      kpi('energy', 'Energy per kg of film on wafer', kWhkg, 'kWh/kg', 'Source power × time divided by the mass that reaches the wafer; excludes pumping and cooling.'),
    ],
    series: [ser('profile', 'Radial thickness profile', 'radius (mm)', 'thickness (nm)', [{ name: 'thickness', x: rs.map(r => r * 1e3), y: th.map(v => v * 1e9) }], 'Cosine-law thickness versus radius; a ring source gives a flatter (or slightly dished) profile than a point source.'),
      ser('throw', 'Uniformity and rate versus throw distance', 'throw distance (cm)', 'value', (() => { const hs = N.linspace(Math.max(2, inp.pvd_h * 0.3), inp.pvd_h * 3, 50), u = [], rr = []; for (const hc of hs) { const H = hc * 1e-2, fl = r => H * H * (H * H + r * r + a * a) / (((H * H - r * r + a * a) ** 2 + 4 * H * H * r * r) ** 1.5); let mx = 0, mn = Infinity; for (const r of rs) { const v = fl(r); if (v > mx) mx = v; if (v < mn) mn = v; } u.push((mx - mn) / (mx + mn) * 100); rr.push(fl(0) / (h * h * (h * h + a * a) / ((h * h + a * a) ** 3)) * 100); } return [{ name: 'non-uniformity (± %)', x: hs, y: u }, { name: 'centre rate (% of present)', x: hs, y: rr }]; })(), 'Trade-off: a longer throw improves uniformity but the rate falls roughly as 1/h².')],
    maps: [{ key: 'tmap', title: 'Thickness map over the wafer', nx: nm, ny: nm, data: map, x0: -Rw * 1e3, x1: Rw * 1e3, y0: -Rw * 1e3, y1: Rw * 1e3, xlabel: 'x (mm)', ylabel: 'y (mm)', unit: 'nm', explain: 'Axisymmetric thickness distribution over the wafer bounding square (stationary substrate, no rotation).', diverging: false }],
    equations: [eqn('formula', evap ? 'Hertz–Knudsen evaporation flux' : 'Sputter emission rate', evap ? 'Γ = P_v(T)/√(2π m k_B T),   P_v = P° exp[−ΔH_vap/R (1/T − 1/T_b)]' : 'Ṅ = Y(E) I /(e (1 + γ_se)),   Y(E) ≈ Y₅₀₀ √(E/500 eV)'),
      eqn('formula', 'Cosine-law thickness from a ring source', 't(r) ∝ h²(h² + r² + a²) / [(h² − r² + a²)² + 4h²r²]^{3/2}'), eqn('formula', 'Mean free path', 'λ = k_B T /(√2 π d² P),   Kn = λ/h'),
      eqn('constitutive', 'Thermal film stress and Stoney curvature', 'σ_th = E_f/(1−ν_f) (α_f − α_s) ΔT;   κ = 6 σ t_f (1−ν_s)/(E_s t_s²)')],
    checks: [chk('Line-of-sight (ballistic) transport', Kn > 1 || evap && Kn > 0.5, `Kn = ${f(Kn)}`), chk('Flux conservation (collected ≤ emitted)', collF <= 1.0001, `collected fraction ${f(collF)}`), chk('Uniformity within ±5 %', unif < 5, `±${f(unif)} %`), chk('Stress below 1 GPa', Math.abs(sTot) < 1e9, `${f(sTot / 1e6)} MPa`)],
    warnings,
    recommendations: [
      rec('high', unif >= 5 ? `Non-uniformity is ±${f(unif)} %: add planetary/substrate rotation, increase the throw distance or use a racetrack radius of about ${f(inp.pvd_wafer * 0.05 * 0.8)} cm.` : `Uniformity ±${f(unif)} % meets a ±5 % specification without rotation; rotation would still average azimuthal asymmetries.`, 'Cosine-law profile over the wafer.'),
      rec('medium', `Measure film stress by wafer curvature: the predicted bow is ${f(bow * 1e6)} µm. Tune working pressure (sputtering) or substrate temperature to move from tensile to compressive intrinsic stress.`, 'Stoney estimate with thermal mismatch.'),
      rec('medium', evap ? 'Control the source with a quartz-crystal monitor in closed loop: the evaporation rate changes by roughly 10 % per 1 % change in source temperature.' : `Operate near ${f(Math.max(0.2, inp.pvd_P))} Pa only as needed for plasma stability; lower pressure raises Kn (now ${f(Kn)}) and gives denser films and rates closer to the ballistic estimate.`, 'Vapour-pressure sensitivity / mean-free-path estimate.'),
      rec('medium', `Only ${f(collF * 100)} % of the ${el} reaches the wafer. Fit removable shields and reclaim the metal from shields and spent targets${['Au', 'Ag', 'Pt'].includes(el) ? ' (precious metal — refine and recycle)' : ''}; larger batch area per run improves material and energy efficiency (${f(kWhkg)} kWh/kg now).`, 'Collected fraction from integrating the flux over the wafer.', true),
    ],
    publish: { 'geom.thickness': tC, 'fab.dep_rate': tC / t, 'growth.rate': tC / t, 'mech.max_stress': Math.abs(sTot), 'proc.power': powerW, 'proc.P': inp.pvd_P, 'proc.time': t, 'fab.yield': N.clamp(collF, 0, 1), 'fab.energy_kWh_kg': kWhkg, 'material.formula': el },
    exports: [{ name: 'pvd_recipe.json', mime: 'application/json', content: JSON.stringify({ process: evap ? 'thermal evaporation' : 'DC magnetron sputtering', material: el, power_W: powerW, pressure_Pa: inp.pvd_P, throw_cm: inp.pvd_h, time_s: t, predicted: { centre_thickness_nm: tC * 1e9, nonuniformity_pct: unif, stress_MPa: sTot / 1e6 } }, null, 2), description: 'PVD recipe with predicted thickness, uniformity and stress.' }],
    numerics: { method: 'Analytical cosine-law ring-source integration; Hertz–Knudsen / tabulated-yield source terms', dof: nm * nm, steps: 1, iterations: 1, residual: 0, notes: 'Collisionless transport.' },
  });
}

// ================================================================================================
// 5. Lithography and etch
// ================================================================================================
const LSRC = { i365: { lam: 365, n: 1, name: 'i-line 365 nm' }, krf: { lam: 248, n: 1, name: 'KrF 248 nm' }, arf: { lam: 193, n: 1, name: 'ArF 193 nm (dry)' }, arfi: { lam: 193, n: 1.44, name: 'ArF 193 nm immersion' }, euv: { lam: 13.5, n: 1, name: 'EUV 13.5 nm' } };
function runLitho(inp) {
  const tool = inp.lit_tool, cd = inp.lit_cd, pitch = Math.max(inp.lit_pitch, cd * 1.05), tol = inp.lit_tol / 100, warnings = [], kp = [], srs = [], eqs = [], cks = [], maps = [];
  if (pitch !== inp.lit_pitch) warnings.push(`Pitch raised to ${f(pitch)} nm so that it exceeds the feature size.`);
  let cdP = cd, litTxt = '';
  if (tool === 'optical') {
    const sr = LSRC[inp.lit_src] || LSRC.arf, lam = sr.lam, NA = Math.min(inp.lit_NA, sr.n * 0.98), res = inp.lit_k1 * lam / NA, dofR = inp.lit_k2 * lam * sr.n / (NA * NA), k1act = pitch / 2 * NA / lam, captured = lam / pitch <= NA * (1 + inp.lit_sigma);
    if (NA !== inp.lit_NA) warnings.push(`NA limited to ${f(NA)} (cannot exceed the refractive index of the imaging medium).`);
    const sp = pitch - cd, a0 = sp / pitch, a1 = Math.sin(PI * sp / pitch) / PI, q = lam / (sr.n * pitch), three = lam / pitch <= NA, sg = Math.max(inp.lit_sigma, 1e-6);
    // captured-source fraction for two-beam (off-axis) imaging and the direction cosines of the interfering orders
    const fsrc = three ? 1 : N.clamp((NA * (1 + sg) - lam / pitch) / (2 * sg * NA), 0, 1), sc0 = three ? 0 : (lam / pitch - NA + sg * NA) / 2 / sr.n, s1 = sc0 - q, cdir = v => Math.sqrt(Math.max(0, 1 - v * v));
    const dphi = dz => 2 * PI * dz * sr.n / lam * (three ? 1 - cdir(q) : cdir(sc0) - cdir(s1));
    const b1 = Math.exp(-2 * (PI * inp.lit_blur / pitch) ** 2), b2 = b1 ** 4; // Gaussian resist blur attenuates the 1st and 2nd image harmonics
    const I = (x, dz) => three ? a0 * a0 + 2 * a1 * a1 + 4 * a0 * a1 * b1 * Math.cos(dphi(dz)) * Math.cos(2 * PI * x / pitch) + 2 * a1 * a1 * b2 * Math.cos(4 * PI * x / pitch) : a0 * a0 + fsrc * (a1 * a1 + 2 * a0 * a1 * b1 * Math.cos(dphi(dz)) * Math.cos(2 * PI * x / pitch)), Ith = I(sp / 2, 0);
    const cdAt = (D, dz) => { if (!captured) return 0; const g = x => I(x, dz) * D - Ith; if (g(0) <= 0) return pitch; if (g(pitch / 2) >= 0) return 0; return pitch - 2 * N.bisect(g, 0, pitch / 2, 1e-6); };
    const ok = (D, dz) => Math.abs(cdAt(D, dz) - cd) <= tol * cd;
    let dLo = 1, dHi = 1; while (dLo > 0.4 && ok(dLo - 0.002, 0)) dLo -= 0.002; while (dHi < 1.8 && ok(dHi + 0.002, 0)) dHi += 0.002;
    let zOk = 0; const zStep = Math.max(dofR, 1) / 200; while (zOk < 6 * dofR && ok(1, zOk + zStep)) zOk += zStep;
    const EL = captured ? (dHi - dLo) * 100 : 0, dofPW = captured ? 2 * zOk : 0, h_ = 1e-3, nils = captured ? cd * Math.abs((Math.log(I(sp / 2 + h_, 0)) - Math.log(I(sp / 2 - h_, 0))) / (2 * h_)) : 0;
    cdP = captured ? cdAt(1, 0) : 0; litTxt = `${sr.name}, NA ${f(NA)}`;
    const nm = 31, mp = new Float64Array(nm * nm), zr = Math.max(1.5 * dofR, 1); for (let j = 0; j < nm; j++) for (let i = 0; i < nm; i++) mp[i + nm * j] = captured ? N.clamp((cdAt(0.8 + 0.4 * i / (nm - 1), -zr + 2 * zr * j / (nm - 1)) - cd) / cd * 100, -100, 100) : -100;
    maps.push({ key: 'pw', title: 'CD error over dose and focus', nx: nm, ny: nm, data: mp, x0: 0.8, x1: 1.2, y0: -zr, y1: zr, xlabel: 'relative dose', ylabel: 'defocus (nm)', unit: '%', explain: 'Printed line-width error over the dose–focus plane (Bossung-type map). The process window is the region within the CD tolerance.', diverging: true });
    const xs = N.linspace(-pitch / 2, pitch / 2, 121);
    srs.push(ser('aerial', 'Aerial image', 'position (nm)', 'relative intensity', [{ name: 'in focus', x: xs, y: xs.map(x => I(x, 0) * (captured ? 1 : 0) + (captured ? 0 : a0 * a0)) }, { name: 'at Rayleigh defocus', x: xs, y: xs.map(x => captured ? I(x, dofR / 2) : a0 * a0) }, { name: 'resist threshold', x: [xs[0], xs[120]], y: [Ith, Ith] }], 'Three-beam aerial image of the line/space grating; the printed edge is where intensity crosses the threshold. Loss of modulation with defocus shrinks the window.'));
    const ds = N.linspace(0.7, 1.3, 61); srs.push(ser('cdDose', 'Printed CD versus dose', 'relative dose', 'line CD (nm)', [{ name: 'in focus', x: ds, y: ds.map(D => cdAt(D, 0)) }, { name: 'defocus = DOF/2', x: ds, y: ds.map(D => cdAt(D, dofR / 2)) }], 'Dose sensitivity of the line width (positive resist: more dose, narrower line). The slope sets the exposure latitude.'));
    kp.push(kpi('res', 'Rayleigh resolution (half-pitch)', res, 'nm', `k₁λ/NA with k₁ = ${f(inp.lit_k1)}; the requested half-pitch ${f(pitch / 2)} nm corresponds to k₁ = ${f(k1act)} (below 0.25 is impossible in single exposure).`, pitch / 2 >= res),
      kpi('dof', 'Rayleigh depth of focus', dofR, 'nm', `k₂λ/NA² (k₂ = ${f(inp.lit_k2)}); total focus budget including wafer topography and tool focus control must fit inside it.`),
      kpi('EL', 'Exposure latitude at best focus', EL, '%', `Dose range over which the printed CD stays within ±${f(inp.lit_tol)} % (threshold resist, ${three ? 'three-beam' : 'two-beam off-axis'} image).`, EL > 8),
      kpi('nils', 'Normalised image log-slope', nils, '–', `CD × |d ln I/dx| at the feature edge${three ? '' : ` with ${f(fsrc * 100)} % of the source contributing to the image`}; values below about 1.5 are not manufacturable, above 2 are comfortable.`, nils > 1.5),
      kpi('dofPW', 'Process-window depth of focus', dofPW, 'nm', 'Total defocus range at nominal dose over which the CD stays in tolerance, from the ideal image model (an upper bound: aberrations, resist thickness and stochastic effects reduce it).', dofPW > 100 || sr.lam < 20));
    eqs.push(eqn('formula', 'Rayleigh criteria', 'R = k₁ λ/NA,   DOF = k₂ n λ/NA²'), eqn('governing', three ? 'Three-beam aerial image (line/space)' : 'Two-beam off-axis aerial image (line/space)', three ? 'I(x, Δz) = a₀² + 2a₁² + 4a₀a₁ cos(Δφ) cos(2πx/p) + 2a₁² cos(4πx/p),   Δφ = (2πnΔz/λ)(1 − √(1 − (λ/np)²))' : 'I(x, Δz) = a₀² + f_src [a₁² + 2a₀a₁ cos(Δφ) cos(2πx/p)],   f_src = (NA(1+σ) − λ/p)/(2σNA),   Δφ = (2πnΔz/λ)(cos θ₀ − cos θ₁)', 'a₀ = s/p and a₁ = sin(πs/p)/π are the mask diffraction amplitudes for space width s.'), eqn('constitutive', 'Threshold resist with Gaussian blur', 'edge where D·I(x) = I_th (dose-to-size at best focus); harmonic m of the image is attenuated by exp(−2π²m²σ_b²/p²)'));
    cks.push(chk('Image contrast sufficient (NILS ≥ 1.5)', nils >= 1.5, `NILS = ${f(nils)}`), chk('First diffraction orders captured', captured, `λ/p = ${f(lam / pitch)} vs NA(1+σ) = ${f(NA * (1 + inp.lit_sigma))}`), chk('k₁ ≥ 0.25', k1act >= 0.25, `k₁ = ${f(k1act)}`), chk('Exposure latitude ≥ 8 %', EL >= 8, `${f(EL)} %`));
    if (!captured) warnings.push('The pitch is below the resolution limit of this wavelength/NA: no image forms and nothing prints. Use a shorter wavelength, higher NA, off-axis illumination or multiple patterning.');
  } else if (tool === 'ebeam') {
    const al = inp.eb_alpha, eta = inp.eb_eta, RKO = 0.0276 * 28.09 * inp.eb_kV ** 1.67 / (14 ** 0.889 * 2.33), beta = 0.4 * RKO * 1e3, dens = cd / pitch, D = inp.eb_dose, D0 = inp.eb_D0;
    const Ex = x => D / (1 + eta) * (0.5 * (N.erf((x + cd / 2) / al) - N.erf((x - cd / 2) / al)) + eta * dens);
    const merged = Ex(pitch / 2) >= D0, none = Ex(0) < D0; cdP = none ? 0 : merged ? pitch : 2 * N.bisect(x => Ex(x) - D0, 0, pitch / 2, 1e-6);
    const Dsize = D0 * (1 + eta) / (0.5 * N.erf(cd / al) + eta * dens), fog = eta * dens / (1 + eta) * D, tWrite = D * 1e-2 * inp.lit_area * 1e-6 * dens / (inp.eb_I * 1e-9);
    litTxt = `e-beam ${f(inp.eb_kV)} kV`;
    const xs = N.linspace(0, pitch / 2, 100), Ds = N.linspace(0.5 * Dsize, 2 * Dsize, 60);
    srs.push(ser('dose', 'Absorbed dose profile across a line', 'distance from line centre (nm)', 'dose (µC/cm²)', [{ name: 'deposited dose', x: xs, y: xs.map(Ex) }, { name: 'resist threshold', x: [0, pitch / 2], y: [D0, D0] }], 'Forward-scattered (Gaussian α) edge profile on top of the uniform backscatter background from the array density. The printed edge is the threshold crossing.'));
    srs.push(ser('cdDose', 'Printed CD versus dose', 'dose (µC/cm²)', 'line CD (nm)', [{ name: 'printed CD', x: Ds, y: Ds.map(d => { const E2 = x => d / (1 + eta) * (0.5 * (N.erf((x + cd / 2) / al) - N.erf((x - cd / 2) / al)) + eta * dens); return E2(0) < D0 ? 0 : E2(pitch / 2) >= D0 ? pitch : 2 * N.bisect(x => E2(x) - D0, 0, pitch / 2, 1e-6); }) }, { name: 'design CD', x: [Ds[0], Ds[59]], y: [cd, cd] }], 'Dose-to-size curve for this pattern density; proximity correction amounts to choosing the dose where the curve meets the design CD for each local density.'));
    kp.push(kpi('beta', 'Backscatter range β', beta, 'nm', `0.4 × Kanaya–Okayama electron range in silicon at ${f(inp.eb_kV)} kV (empirical proportionality); features within β of each other share backscattered dose.`),
      kpi('Dsize', 'Dose-to-size for this density', Dsize, 'µC/cm²', `Dose at which the threshold contour coincides with the design edge for pattern density ${f(dens)}; the entered dose is ${f(D)} µC/cm².`, Math.abs(D - Dsize) / Dsize < 0.1),
      kpi('fog', 'Backscatter background', fog, 'µC/cm²', 'Uniform dose added by backscattering in a dense array, ηρ/(1+η)·D; it must stay below the resist threshold between features.', fog < D0),
      kpi('twrite', 'Write time', tWrite / 60, 'min', `Dose × exposed area (${f(inp.lit_area)} mm² × density) / beam current ${f(inp.eb_I)} nA; excludes stage and settling overhead.`));
    eqs.push(eqn('governing', 'Double-Gaussian proximity function', 'f(r) = 1/(π(1+η)) [ exp(−r²/α²)/α² + η exp(−r²/β²)/β² ]'), eqn('formula', 'Dose across a line in a dense array (β ≫ pitch)', 'E(x) = D/(1+η) { ½[erf((x+w/2)/α) − erf((x−w/2)/α)] + ηρ }'), eqn('formula', 'Backscatter range', 'β ≈ 0.4 R_KO,   R_KO = 0.0276 A E^1.67 /(Z^0.889 ρ) µm'));
    cks.push(chk('Feature prints (dose above threshold)', !none, `centre dose ${f(Ex(0))} vs threshold ${f(D0)} µC/cm²`), chk('Features resolved (no bridging)', !merged, `gap dose ${f(Ex(pitch / 2))} vs threshold ${f(D0)} µC/cm²`), chk('Forward scattering smaller than half the gap', al < (pitch - cd) / 2, `α = ${f(al)} nm, gap ${f(pitch - cd)} nm`), chk('Array larger than β (dense-array assumption)', Math.sqrt(inp.lit_area) * 1e6 > 2 * beta, `field ${f(Math.sqrt(inp.lit_area) * 1e6)} nm vs β ${f(beta)} nm`));
    if (merged) warnings.push('Backscatter plus forward scattering exceed the threshold in the gaps: the lines bridge. Reduce dose or pattern density, or apply proximity correction.');
    if (none) warnings.push('Dose is below the resist threshold even at the line centre: nothing prints.');
  } else {
    const s_ = cd * 1e-9, h0 = inp.nil_h0 * 1e-9, hf = Math.min(inp.nil_hf, inp.nil_h0 * 0.999) * 1e-9, tf = inp.nil_visc * s_ * s_ / (2 * inp.nil_p * 1e6) * (1 / (hf * hf) - 1 / (h0 * h0));
    cdP = cd; litTxt = 'nanoimprint';
    const hs = N.linspace(hf * 1e9, h0 * 1e9, 60);
    srs.push(ser('fill', 'Imprint time versus residual layer', 'residual layer thickness (nm)', 'time (s)', [{ name: 'squeeze-flow time', x: hs, y: hs.map(hh => Math.max(1e-12, inp.nil_visc * s_ * s_ / (2 * inp.nil_p * 1e6) * (1 / (hh * 1e-9) ** 2 - 1 / (h0 * h0)))) }], 'Stefan squeeze-flow time to thin the resist under a stamp protrusion; the time rises steeply as the residual layer gets thin.', { logy: true }));
    srs.push(ser('fillw', 'Imprint time versus feature width', 'protrusion width (nm)', 'time (s)', (() => { const ws = lsp(1, 4, 50); return [{ name: 'fill time', x: ws, y: ws.map(w => inp.nil_visc * (w * 1e-9) ** 2 / (2 * inp.nil_p * 1e6) * (1 / (hf * hf) - 1 / (h0 * h0))) }]; })(), 'Fill time scales with the square of the protrusion width: large pads, not the nano-features, set the imprint time.', { logx: true, logy: true }));
    kp.push(kpi('tfill', 'Imprint (squeeze-flow) time', tf, 's', `Stefan equation for a ${f(cd)} nm wide protrusion thinning the resist from ${f(inp.nil_h0)} to ${f(hf * 1e9)} nm at ${f(inp.nil_p)} MPa and ${f(inp.nil_visc)} Pa·s.`, tf < 600),
      kpi('residual', 'Residual layer thickness', hf * 1e9, 'nm', 'Resist left under the protrusions; it must be removed by a breakthrough etch before pattern transfer and adds to CD loss.'),
      kpi('vol', 'Resist volume check (cavity fill)', (inp.nil_h0 - hf * 1e9) * cd / Math.max(pitch - cd, 1e-9), 'nm', 'Height to which the displaced resist fills the stamp cavities (displaced volume / cavity width); it must not exceed the cavity depth.'));
    eqs.push(eqn('governing', 'Stefan squeeze flow', 't_f = (η s² / 2p) (1/h_f² − 1/h₀²)'));
    cks.push(chk('Imprint time practical (< 10 min)', tf < 600, `${f(tf)} s`));
  }
  // ---- etch transfer ----
  const space = Math.max(pitch - cdP, 1e-6), r0 = inp.et_rate / 60, Sr = inp.et_S, target = inp.et_depth, An = inp.et_A;
  const fAR = ar => { const K = 1 / (1 + 0.75 * ar); return K / (K + Sr - K * Sr); };
  let hE = 0, tE = 0; const dtE = target / r0 / 400, tt = [0], hh = [0]; let it = 0;
  while (hE < target && it < 8000 && cdP > 0 && cdP < pitch) { hE += r0 * fAR(hE / space) * dtE; tE += dtE; it++; if (it % 20 === 0) { tt.push(tE); hh.push(hE); } }
  tt.push(tE); hh.push(hE);
  const under = (1 - An) * hE, cdF = Math.max(0, cdP - 2 * under), maskUse = r0 * tE / inp.et_sel, arF = hE / space, slow = fAR(arF), wall = 90 - Math.atan(1 - An) * 180 / PI, bias = cdF - cd, inTol = Math.abs(bias) <= tol * cd && cdP > 0 && cdP < pitch;
  kp.push(kpi('cdP', 'Printed CD (after lithography)', cdP, 'nm', `Line width in resist from the ${tool === 'optical' ? 'aerial-image/threshold' : tool === 'ebeam' ? 'proximity/threshold' : '1:1 imprint'} model; design value ${f(cd)} nm.`, Math.abs(cdP - cd) <= tol * cd),
    kpi('tetch', 'Etch time to target depth', tE, 's', `Integration of the aspect-ratio-dependent rate to ${f(target)} nm; the open-area time would be ${f(target / r0)} s.`),
    kpi('arde', 'Rate at final aspect ratio', slow * 100, '% of open-area rate', `Knudsen-transport (Coburn–Winters) slow-down at aspect ratio ${f(arF)} with bottom reaction probability ${f(Sr)}; it causes depth differences between wide and narrow features (RIE lag).`, slow > 0.7),
    kpi('undercut', 'Undercut per side', under, 'nm', `(1 − anisotropy) × depth with anisotropy ${f(An)}; sidewall angle ≈ ${f(wall)}°.`, under < 0.1 * cd),
    kpi('mask', 'Mask consumed', maskUse, 'nm', `Etched depth of the mask at selectivity ${f(inp.et_sel)}:1 versus the ${f(inp.et_resist)} nm available.`, maskUse < 0.8 * inp.et_resist),
    kpi('cdF', 'Final CD after etch', cdF, 'nm', `Printed CD minus twice the undercut: total bias ${f(bias)} nm (${f(bias / cd * 100)} %) against the ±${f(inp.lit_tol)} % tolerance.`, inTol));
  srs.push(ser('etch', 'Etch depth versus time', 'time (s)', 'depth (nm)', [{ name: 'feature (ARDE)', x: tt, y: hh }, { name: 'open area', x: tt, y: tt.map(v => r0 * v) }], 'Depth in the patterned trenches lags the open-area depth as the aspect ratio grows.'));
  eqs.push(eqn('governing', 'Aspect-ratio-dependent etching', 'dh/dt = r₀ K/(K + S − K S),   K = 1/(1 + ¾·h/w)'), eqn('formula', 'Anisotropy and selectivity', 'A = 1 − r_lat/r_vert;   undercut = (1 − A) h;   mask loss = r₀ t / selectivity'));
  cks.push(chk('Final CD within tolerance', inTol, `bias ${f(bias)} nm vs ±${f(tol * cd)} nm`), chk('Mask survives the etch', maskUse < inp.et_resist, `${f(maskUse)} nm consumed of ${f(inp.et_resist)} nm`), chk('Target depth reached', hE >= target * 0.999, `${f(hE)} nm of ${f(target)} nm`));
  if (maskUse >= inp.et_resist) warnings.push('The mask is consumed before the target depth is reached; use a hard mask or a more selective chemistry.');
  warnings.push('Threshold-resist and analytical image/proximity models: no resist diffusion, development kinetics, mask 3-D effects or stochastic (shot-noise/LER) effects, which dominate below about 20 nm.');
  const nA = 6, polys = [], off = pitch / 2; for (let j = 0; j < nA; j++) for (let i = 0; i < nA; i++) { const cx = off + i * pitch, cy = off + j * pitch, sq = (w, layer) => polys.push({ layer, pts: [[cx - w / 2, cy - w / 2], [cx + w / 2, cy - w / 2], [cx + w / 2, cy + w / 2], [cx - w / 2, cy + w / 2]].map(p => p.map(v => +v.toFixed(2))) }); sq(cd, 2); if (cdF > 0) sq(cdF, 1); }
  return finish({
    summary: `${litTxt}: a ${f(cd)} nm feature at ${f(pitch)} nm pitch prints at ${f(cdP)} nm and, after a ${f(tE)} s etch to ${f(hE)} nm (aspect ratio ${f(arF)}, ${f(under)} nm undercut per side), ends at ${f(cdF)} nm — ${inTol ? 'inside' : 'outside'} the ±${f(inp.lit_tol)} % CD tolerance. Mask consumption is ${f(maskUse)} nm of ${f(inp.et_resist)} nm.`,
    kpis: kp, series: srs, maps, equations: eqs, checks: cks, warnings,
    geometry: { kind: 'layout', polygons: polys, units: 'nm', title: 'Printed and etched array (layer 1) over the design (layer 2)' },
    recommendations: [
      rec('high', inTol ? `Final CD bias of ${f(bias)} nm is within tolerance; lock dose/focus at the centre of the window and monitor with CD-SEM on every lot.` : `Final CD is off by ${f(bias)} nm: bias the design by ${f(-bias)} nm (mask/layout sizing) or ${under > 0.05 * cd ? 'raise etch anisotropy (lower pressure, sidewall passivation)' : 'retarget the exposure dose'}.`, 'Printed CD and undercut from this run.'),
      rec('medium', slow < 0.7 ? `Etch rate falls to ${f(slow * 100)} % at aspect ratio ${f(arF)}; lower the pressure or increase ion energy/time-multiplex the etch to reduce RIE lag between features of different width.` : `RIE lag is modest (${f(slow * 100)} % of the open-area rate at the trench bottom); verify depth on the narrowest and widest features.`, 'Knudsen-transport ARDE model.'),
      rec('medium', maskUse > 0.6 * inp.et_resist ? `Mask budget is tight (${f(maskUse)} of ${f(inp.et_resist)} nm): add a hard mask (SiO₂, Cr or Al₂O₃) or raise selectivity.` : 'Mask budget is comfortable; thinner resist would improve resolution and process window.', 'Selectivity-based mask consumption.'),
      rec('low', 'Minimise fluorinated etch-gas use (SF₆, CF₄, C₄F₈ have very high global-warming potentials): optimise etch time with end-point detection, and run point-of-use abatement on the exhaust.', 'Good practice for plasma etching; etch time from this run sets the gas consumption.', true),
    ],
    publish: { 'fab.cd_m': cdF * 1e-9, 'geom.cell_size': pitch * 1e-9, 'geom.thickness': hE * 1e-9, 'geom.aspect': arF, 'proc.time': tE, 'fab.yield': inTol ? 1 : 0, 'fab.min_feature_m': cdP * 1e-9 },
    exports: [{ name: 'litho_etch_recipe.json', mime: 'application/json', content: JSON.stringify({ lithography: litTxt, design_cd_nm: cd, pitch_nm: pitch, printed_cd_nm: cdP, etch: { rate_nm_min: inp.et_rate, time_s: tE, depth_nm: hE, selectivity: inp.et_sel, anisotropy: An }, final_cd_nm: cdF, in_tolerance: inTol }, null, 2), description: 'Lithography and etch recipe with predicted dimensions.' }],
    numerics: { method: 'Analytical aerial-image / double-Gaussian proximity / Stefan models with threshold resist (bisection for edges); explicit Euler integration of the ARDE rate law', dof: 1, steps: it, iterations: it, residual: 0, notes: '400 etch steps to target depth.' },
  });
}

// ================================================================================================
// 6. Additive manufacturing of lattices
// ================================================================================================
// rho kg/m³, cp J/kg/K and k W/m/K (elevated-temperature averages), Tm K, Tb K, absorptivity
const AMM = { '316L': [7980, 600, 20, 1700, 3090, 0.35, '316L stainless steel'], 'Ti64': [4420, 610, 15, 1923, 3560, 0.40, 'Ti-6Al-4V'], 'AlSi10Mg': [2670, 915, 130, 870, 2740, 0.15, 'AlSi10Mg'], 'IN718': [8190, 600, 20, 1610, 3190, 0.35, 'Inconel 718'] };
function runAM(inp) {
  const proc = inp.am_proc, warnings = [], rd = inp.am_rd, cell = inp.am_cell * 1e-3, dStrut = 2 * cell * Math.sqrt(rd / (4 * Math.sqrt(3) * PI));
  const common = [eqn('formula', 'Lattice strut diameter (BCC cell, low density)', 'ρ* = 4√3 π (r/a)²  →  d = 2a √(ρ*/(4√3π))')];
  if (proc === 'lpbf') {
    const key = AMM[inp.am_mat] ? inp.am_mat : '316L', [rho, cp, k, Tm, Tb, Ab, nm] = AMM[key], Pw = inp.am_P, v = inp.am_v * 1e-3, hs = inp.am_h * 1e-6, tl = inp.am_t * 1e-6, sig = inp.am_spot * 1e-6 / 4, T0 = inp.am_T0 + 273.15;
    const VED = Pw / (v * hs * tl) / 1e9, al = k / (rho * cp), dT = Math.max(Tm - T0, 1), rTail = Ab * Pw / (2 * PI * k * dT); let yMax = 0;
    for (let i = 1; i <= 600; i++) { const r = rTail * i / 600, xi = -r - 2 * al / v * Math.log(r / rTail), y2 = r * r - xi * xi; if (y2 > yMax) yMax = y2; }
    const W = 2 * Math.sqrt(yMax), Dp = W / 2, lof = (hs / W) ** 2 + (tl / Dp) ** 2, Hn = Ab * Pw / (rho * cp * Tm * Math.sqrt(PI * al * v * sig ** 3)), Hb = PI * Tb / Tm, mode = Hn > 30 ? 'keyhole risk' : Hn > Hb ? 'transition (vapour depression)' : 'conduction';
    const vol = inp.am_vol * 1e-6, Vs = rd * vol, height = Math.cbrt(vol), layers = Math.ceil(height / tl), tScan = Vs / (v * hs * tl), tBuild = tScan + layers * inp.am_recoat, minStrut = W + inp.am_powder * 1e-6, mass = Vs * rho, kWhkg = Pw * tScan / 3.6e6 / Math.max(mass, 1e-300);
    const vs = lsp(Math.log10(Math.max(50, inp.am_v / 5)), Math.log10(inp.am_v * 5), 60), wv = vs.map(vv => { let ym = 0; const V_ = vv * 1e-3; for (let i = 1; i <= 300; i++) { const r = rTail * i / 300, xi = -r - 2 * al / V_ * Math.log(r / rTail), y2 = r * r - xi * xi; if (y2 > ym) ym = y2; } return 2e6 * Math.sqrt(ym); });
    const nm_ = 31, mp = new Float64Array(nm_ * nm_), P0 = 50, P1 = Math.max(400, Pw * 1.5), v0 = 100, v1 = Math.max(2000, inp.am_v * 1.5);
    for (let j = 0; j < nm_; j++) for (let i = 0; i < nm_; i++) { const Pp = P0 + (P1 - P0) * j / (nm_ - 1), V_ = (v0 + (v1 - v0) * i / (nm_ - 1)) * 1e-3, rt = Ab * Pp / (2 * PI * k * dT); let ym = 0; for (let q = 1; q <= 120; q++) { const r = rt * q / 120, xi = -r - 2 * al / V_ * Math.log(r / rt), y2 = r * r - xi * xi; if (y2 > ym) ym = y2; } const w = 2 * Math.sqrt(ym), lf = (hs / w) ** 2 + (tl / (w / 2)) ** 2, hn = Ab * Pp / (rho * cp * Tm * Math.sqrt(PI * al * V_ * sig ** 3)); mp[i + nm_ * j] = lf > 1 ? -1 : hn > 30 ? 1 : 0; }
    warnings.push('Rosenthal point-source conduction solution with constant properties: no latent heat, convection (Marangoni) or vapour recoil, so melt-pool depth is always half the width and is under-predicted in keyhole mode. Use it for trends and window screening, then confirm with single-track trials.');
    if (lof > 1) warnings.push(`Lack-of-fusion criterion (h/W)² + (t/D)² = ${f(lof)} > 1: adjacent tracks and layers do not overlap fully; expect irregular porosity.`);
    if (Hn > 30) warnings.push(`Normalised enthalpy ${f(Hn)} exceeds about 30, the empirical keyhole-porosity threshold reported for 316L; treat the threshold as machine and alloy dependent.`);
    if (dStrut < minStrut) warnings.push(`Design strut ${f(dStrut * 1e6)} µm is below the printable minimum ≈ ${f(minStrut * 1e6)} µm (melt-pool width + one powder diameter).`);
    return finish({
      summary: `LPBF of ${nm} at ${f(Pw)} W and ${f(inp.am_v)} mm/s gives a volumetric energy density of ${f(VED)} J/mm³ and a Rosenthal melt pool ${f(W * 1e6)} µm wide × ${f(Dp * 1e6)} µm deep (${mode}; lack-of-fusion index ${f(lof)}). The ${f(rd * 100)} % dense lattice in a ${f(inp.am_vol)} cm³ envelope (strut ${f(dStrut * 1e6)} µm) builds in ${f(tBuild / 3600)} h over ${layers} layers.`,
      kpis: [
        kpi('ved', 'Volumetric energy density', VED, 'J/mm³', 'P/(v·h·t). A convenient process index, but different P–v pairs with equal VED give different melt pools; use it with the melt-pool criteria.'),
        kpi('W', 'Melt-pool width', W * 1e6, 'µm', `Widest extent of the melting isotherm of the Rosenthal moving point source with absorptivity ${f(Ab)} and preheat ${f(inp.am_T0)} °C.`),
        kpi('D', 'Melt-pool depth', Dp * 1e6, 'µm', `Half the width for the semi-circular Rosenthal pool; should exceed the layer thickness (${f(inp.am_t)} µm) by about 1.5× for sound bonding.`, Dp > 1.5 * tl),
        kpi('lof', 'Lack-of-fusion index', lof, '–', '(h/W)² + (t/D)²; values ≤ 1 mean hatch spacing and layer thickness are fully covered by overlapping melt pools.', lof <= 1),
        kpi('Hn', 'Normalised enthalpy ΔH/h_s', Hn, '–', `A·P/(ρc_pT_m√(πα v σ³)); surface boiling starts near πT_b/T_m = ${f(Hb)} and deep keyholing is reported above about 30.`, Hn <= 30),
        kpi('mode', 'Melting mode', mode, '', 'Conduction below the surface-boiling onset, transition between onset and about 30, keyhole-porosity risk above.'),
        kpi('minStrut', 'Minimum printable strut', minStrut * 1e6, 'µm', `Melt-pool width plus one powder diameter (${f(inp.am_powder)} µm) of adhered particles; the lattice needs ${f(dStrut * 1e6)} µm.`, dStrut >= minStrut),
        kpi('tbuild', 'Build time', tBuild / 3600, 'h', `Scan time ${f(tScan / 3600)} h (solid volume / (v·h·t)) plus ${layers} recoats of ${f(inp.am_recoat)} s.`),
        kpi('energy', 'Laser energy per kg of part', kWhkg, 'kWh/kg', 'Laser output energy during scanning divided by part mass; wall-plug energy is typically 5–10× higher (laser efficiency, chiller, inert gas).'),
      ],
      series: [ser('wv', 'Melt-pool width versus scan speed', 'scan speed (mm/s)', 'width (µm)', [{ name: 'melt-pool width', x: vs, y: wv }, { name: 'hatch spacing', x: [vs[0], vs[59]], y: [inp.am_h, inp.am_h] }], 'Faster scanning narrows the melt pool; when it approaches the hatch spacing, lack of fusion appears.', { logx: true }),
        ser('hn', 'Normalised enthalpy versus power', 'laser power (W)', 'ΔH/h_s', (() => { const Ps = N.linspace(50, Math.max(400, Pw * 1.5), 50); return [{ name: 'ΔH/h_s', x: Ps, y: Ps.map(q => Hn * q / Pw) }, { name: 'boiling onset πT_b/T_m', x: [Ps[0], Ps[49]], y: [Hb, Hb] }, { name: 'keyhole threshold ≈ 30', x: [Ps[0], Ps[49]], y: [30, 30] }]; })(), 'Position of the operating point between the surface-boiling onset and the keyhole-porosity threshold as power changes at fixed speed and spot.')],
      maps: [{ key: 'window', title: 'Process window (−1 lack of fusion, 0 dense, +1 keyhole risk)', nx: nm_, ny: nm_, data: mp, x0: v0, x1: v1, y0: P0, y1: P1, xlabel: 'scan speed (mm/s)', ylabel: 'laser power (W)', unit: '', explain: 'Power–speed map of the two defect criteria at the entered hatch, layer and spot size. The zero region is the predicted dense window.', diverging: true }],
      equations: common.concat([eqn('formula', 'Volumetric energy density', 'E_v = P/(v h t)'), eqn('governing', 'Rosenthal moving point source', 'T − T₀ = (A P / 2π k r) exp[−v (r + ξ)/(2α)]'), eqn('formula', 'Lack-of-fusion criterion', '(h/W)² + (t/D)² ≤ 1'), eqn('formula', 'Normalised enthalpy', 'ΔH/h_s = A P /(ρ c_p T_m √(π α v σ³));   boiling onset at π T_b/T_m')]),
      checks: [chk('No lack of fusion', lof <= 1, `index ${f(lof)}`), chk('Below keyhole threshold', Hn <= 30, `ΔH/h_s = ${f(Hn)}`), chk('Strut printable', dStrut >= minStrut, `${f(dStrut * 1e6)} µm vs minimum ${f(minStrut * 1e6)} µm`), chk('Depth ≥ 1.5 layers', Dp >= 1.5 * tl, `${f(Dp / tl)} layers`)],
      warnings,
      recommendations: [
        rec('high', lof > 1 ? `Close the lack-of-fusion gap: reduce hatch to ≤ ${f(W * 0.7 * 1e6)} µm or layer thickness to ≤ ${f(Dp / 1.5 * 1e6)} µm, or raise power.` : Hn > 30 ? `Reduce energy input (lower power or faster scan, or defocus the spot) to bring ΔH/h_s from ${f(Hn)} below about 30.` : `Operating point is inside the predicted dense window (lack-of-fusion index ${f(lof)}, ΔH/h_s ${f(Hn)}); confirm with density cubes and single tracks.`, 'Melt-pool overlap and normalised-enthalpy criteria.'),
        rec('high', dStrut < minStrut ? `Increase strut diameter to ≥ ${f(minStrut * 1e6)} µm (raise relative density to ≥ ${f(4 * Math.sqrt(3) * PI * (minStrut / (2 * cell)) ** 2 * 100)} % or enlarge the cell) and use contour-only scanning for struts.` : `Struts of ${f(dStrut * 1e6)} µm are printable; expect ±${f(inp.am_powder)} µm of adhered powder — compensate the CAD diameter and avoid struts below 30° from horizontal.`, 'Minimum feature = melt-pool width + powder size.'),
        rec('medium', 'Stress-relieve on the plate before cutting (residual stress in LPBF approaches yield) and verify strut diameter and porosity by CT on a witness lattice.', 'Thermal-gradient mechanism of residual stress; outside this model.'),
        rec('medium', `Lattices use only ${f(rd * 100)} % of the envelope as metal: sieve and reuse the unfused powder (track oxygen pick-up), and nest parts to share the ${f(layers * inp.am_recoat / 3600)} h of recoating overhead and inert-gas use.`, 'Build-time breakdown of this run.', true),
      ],
      publish: { 'geom.strut': dStrut, 'geom.rel_density': rd, 'geom.cell_size': cell, 'proc.time': tBuild, 'proc.power': Pw, 'fab.energy_kWh_kg': kWhkg, 'fab.throughput_kg_h': mass / tBuild * 3600, 'fab.min_feature_m': minStrut, 'fab.yield': lof <= 1 && Hn <= 30 ? 1 : 0 },
      exports: [{ name: 'lpbf_recipe.json', mime: 'application/json', content: JSON.stringify({ process: 'LPBF', material: nm, power_W: Pw, speed_mm_s: inp.am_v, hatch_um: inp.am_h, layer_um: inp.am_t, spot_um: inp.am_spot, predicted: { VED_J_mm3: VED, melt_pool_width_um: W * 1e6, melt_pool_depth_um: Dp * 1e6, mode, build_time_h: tBuild / 3600, strut_um: dStrut * 1e6 } }, null, 2), description: 'LPBF parameter set with predicted melt pool, defect mode and build time.' }],
      numerics: { method: 'Rosenthal analytical solution; melt-pool width by parametric scan of the melting isotherm (600 points)', dof: 600, steps: 1, iterations: 1, residual: 0, notes: 'Process-window map evaluated on a 31×31 power–speed grid.' },
    });
  }
  if (proc === 'tpp') {
    const lam = inp.tpp_lam * 1e-9, NAo = Math.min(inp.tpp_NA, inp.tpp_n * 0.98), w0 = 0.61 * lam / NAo, zR = PI * w0 * w0 * inp.tpp_n / lam, v = inp.tpp_v, Pth = inp.tpp_Pth * Math.sqrt(v / 100), ratio = inp.tpp_P / Pth, on = ratio > 1;
    const dv = on ? w0 * Math.sqrt(2 * Math.log(ratio)) : 0, lv = on ? 2 * zR * Math.sqrt(ratio - 1) : 0, side = inp.tpp_side * 1e-6, Vs = rd * side ** 3, path = on ? Vs / (0.25 * dv * lv) : 0, tW = on ? path / (v * 1e-6) : 0, sh = inp.tpp_shrink / 100, cellU = inp.am_cell * 1e-3, strut = dStrut;
    const Ps = N.linspace(Pth * 1.01, Math.max(inp.tpp_P * 1.5, Pth * 3), 60);
    if (!on) warnings.push(`Laser power ${f(inp.tpp_P)} mW is below the polymerisation threshold ${f(Pth)} mW at ${f(v)} µm/s: nothing is written.`);
    if (ratio > 3) warnings.push('Power is more than 3× threshold: micro-explosions (resin boiling) are likely; the dynamic range of most resins is about 2–3.');
    warnings.push('Threshold model with dose ∝ I²·t and a Gaussian focus: no radical diffusion, oxygen quenching, proximity (memory) effect or aberrations; calibrate threshold power at your writing speed with a power–speed line array.');
    return finish({
      summary: on ? `Two-photon polymerisation at ${f(inp.tpp_P)} mW (${f(ratio)}× threshold) and ${f(v)} µm/s writes voxels ${f(dv * 1e9)} nm wide × ${f(lv * 1e9)} nm tall (aspect ${f(lv / dv)}). A ${f(inp.tpp_side)} µm lattice cube at ${f(rd * 100)} % density needs ${f(path)} m of scan path and ${f(tW / 3600)} h; ${f(inp.tpp_shrink)} % shrinkage must be pre-compensated.` : `Laser power ${f(inp.tpp_P)} mW is below the ${f(Pth)} mW threshold at ${f(v)} µm/s, so no polymerisation occurs; raise power or slow the scan.`,
      kpis: [
        kpi('dvox', 'Voxel lateral diameter', dv * 1e9, 'nm', `w₀√(2 ln(P/P_th)) with focal radius w₀ = ${f(w0 * 1e9)} nm (0.61λ/NA); it sets the minimum line width.`),
        kpi('lvox', 'Voxel axial length', lv * 1e9, 'nm', `2z_R√(P/P_th − 1) with Rayleigh range ${f(zR * 1e9)} nm; the elongated voxel limits vertical resolution.`),
        kpi('aspect', 'Voxel aspect ratio', on ? lv / dv : 0, '–', 'Axial/lateral voxel size; typically 2.5–5 for NA 1.4. Hatch and slice distances should be about half the respective voxel dimension.'),
        kpi('ratio', 'Power / threshold', ratio, '–', `Threshold power scales as √v from the entered ${f(inp.tpp_Pth)} mW at 100 µm/s (dose ∝ I²t). Usable window is roughly 1.1–3.`, ratio > 1.1 && ratio < 3),
        kpi('twrite', 'Write time', tW / 3600, 'h', `Solid volume / (¼·d·l per line with 50 % overlap) / speed = ${f(path)} m of path. Galvo-scanned layers remove most stage overhead.`),
        kpi('strut', 'Lattice strut diameter (design)', strut * 1e6, 'µm', `BCC-type strut for ${f(rd * 100)} % density at ${f(inp.am_cell * 1e3)} µm cell size; it must exceed the voxel width (${f(dv * 1e6)} µm).`, strut > dv),
        kpi('shrunk', 'Cube side after shrinkage', inp.tpp_side * (1 - sh), 'µm', `Linear shrinkage ${f(inp.tpp_shrink)} % on development/drying; structures pinned to the substrate distort near the base.`),
      ],
      series: [ser('vox', 'Voxel size versus laser power', 'laser power (mW)', 'size (nm)', [{ name: 'lateral diameter', x: Ps, y: Ps.map(q => w0 * Math.sqrt(2 * Math.log(q / Pth)) * 1e9) }, { name: 'axial length', x: Ps, y: Ps.map(q => 2 * zR * Math.sqrt(q / Pth - 1) * 1e9) }], 'Voxel dimensions grow from zero at threshold; operating close to threshold gives the finest features but the least tolerance to power drift.'),
        ser('tw', 'Write time versus scan speed', 'scan speed (µm/s)', 'write time (h)', (() => { const vv = lsp(2, 5.5, 50); return [{ name: 'write time at fixed P/P_th', x: vv, y: vv.map(q => on ? path / (q * 1e-6) / 3600 : 0) }]; })(), 'Write time falls inversely with speed provided power is raised as √v to keep the same voxel.', { logx: true, logy: on })],
      equations: common.concat([eqn('formula', 'Voxel size (threshold model, dose ∝ I²t)', 'd = w₀ √(2 ln(P/P_th)),   l = 2 z_R √(P/P_th − 1),   w₀ = 0.61 λ/NA,   z_R = π w₀² n/λ'), eqn('formula', 'Threshold scaling with speed', 'P_th(v) = P_th,ref √(v/v_ref)'), eqn('formula', 'Write time', 't = ρ* V /(¼ d l v)')]),
      checks: [chk('Above polymerisation threshold', on, `P/P_th = ${f(ratio)}`), chk('Below damage (≈ 3× threshold)', ratio < 3, `P/P_th = ${f(ratio)}`), chk('Strut wider than voxel', strut > dv, `${f(strut * 1e6)} µm vs ${f(dv * 1e6)} µm`)],
      warnings,
      recommendations: [
        rec('high', !on ? `Raise power above ${f(Pth * 1.2)} mW or slow the scan below ${f(100 * (inp.tpp_P / inp.tpp_Pth / 1.2) ** 2)} µm/s.` : `Set hatch ≈ ${f(dv * 0.5e9)} nm and slice ≈ ${f(lv * 0.5e9)} nm for smooth, fully cross-linked struts.`, 'Voxel dimensions from the threshold model.'),
        rec('medium', `Pre-scale the CAD by ${f(1 / (1 - sh))} to compensate ${f(inp.tpp_shrink)} % shrinkage and use critical-point drying for slender lattices to avoid capillary collapse.`, 'Entered shrinkage and capillary-force practice.'),
        rec('medium', tW > 36000 ? `Write time ${f(tW / 3600)} h is long: use a lower-NA objective for the bulk (shell-and-scaffold with UV flood cure) and the high-NA objective only for fine features.` : `Write time ${f(tW / 3600)} h is practical on a galvo system; larger arrays should be stitched with overlap.`, 'Path-length estimate.'),
        rec('low', 'Develop in the minimum volume of solvent (PGMEA/IPA), collect and distil it for reuse, and UV-cure waste resin before disposal — uncured acrylates are sensitisers and aquatic toxins.', 'Good laboratory practice for photoresins.', true),
      ],
      publish: { 'geom.strut': strut, 'geom.rel_density': rd, 'geom.cell_size': cellU, 'proc.time': tW, 'fab.min_feature_m': dv, 'proc.power': inp.tpp_P * 1e-3, 'fab.yield': on && ratio < 3 ? 1 : 0 },
      exports: [{ name: 'tpp_recipe.json', mime: 'application/json', content: JSON.stringify({ process: 'two-photon polymerisation', wavelength_nm: inp.tpp_lam, NA: NAo, power_mW: inp.tpp_P, speed_um_s: v, hatch_nm: dv * 0.5e9, slice_nm: lv * 0.5e9, cad_scale: 1 / (1 - sh), predicted: { voxel_d_nm: dv * 1e9, voxel_l_nm: lv * 1e9, write_time_h: tW / 3600 } }, null, 2), description: 'Two-photon writing parameters with hatch/slice and scale compensation.' }],
      numerics: { method: 'Closed-form two-photon threshold voxel model', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Gaussian focus, dose ∝ I²t.' },
    });
  }
  // stereolithography (Jacobs)
  const Dp = inp.sla_Dp, Ec = inp.sla_Ec, Ee = inp.sla_E, cure = Ee > Ec ? Dp * Math.log(Ee / Ec) : 0, layer = inp.sla_layer, over = cure - layer, lw = cure > 0 ? inp.sla_w0 * Math.sqrt(2 * cure / Dp) : 0;
  const vol = inp.am_vol * 1e-6, height = Math.cbrt(vol), layers = Math.ceil(height / (layer * 1e-6)), aSolid = rd * height * height, tLayer = aSolid * Ee * 10 / (inp.sla_P * 1e-3) + inp.am_recoat, tB = layers * tLayer, Es = lsp(Math.log10(Ec * 1.01), Math.log10(Math.max(Ee * 3, Ec * 30)), 60);
  if (cure <= 0) warnings.push('Exposure is below the critical exposure E_c: the resin does not gel.');
  if (over < 0) warnings.push(`Cure depth ${f(cure)} µm is less than the layer thickness ${f(layer)} µm: layers will not bond (delamination).`);
  if (dStrut * 1e6 < lw) warnings.push(`Design strut ${f(dStrut * 1e6)} µm is narrower than the cured line width ${f(lw)} µm.`);
  warnings.push('Jacobs working-curve model with Beer–Lambert absorption and a Gaussian beam: scattering by ceramic fillers, oxygen inhibition and dark cure are not included.');
  return finish({
    summary: `Stereolithography at ${f(Ee)} mJ/cm² cures ${f(cure)} µm deep (D_p = ${f(Dp)} µm, E_c = ${f(Ec)} mJ/cm²), i.e. ${f(over)} µm overcure on ${f(layer)} µm layers, with a cured line width of ${f(lw)} µm. The ${f(rd * 100)} % dense lattice in a ${f(inp.am_vol)} cm³ envelope takes ${layers} layers and ${f(tB / 3600)} h.`,
    kpis: [
      kpi('cure', 'Cure depth', cure, 'µm', 'Jacobs working curve C_d = D_p ln(E/E_c): logarithmic in exposure, so doubling the dose adds only 0.69·D_p.'),
      kpi('over', 'Overcure into previous layer', over, 'µm', 'Cure depth minus layer thickness; 10–50 % of the layer thickness is usual for good interlayer bonding without excessive print-through.', over > 0.1 * layer && over < layer),
      kpi('lw', 'Cured line width', lw, 'µm', `w₀√(2C_d/D_p) for a Gaussian beam of 1/e² radius ${f(inp.sla_w0)} µm; it sets the minimum wall or strut.`),
      kpi('strut', 'Lattice strut diameter (design)', dStrut * 1e6, 'µm', `BCC-type strut for ${f(rd * 100)} % density at ${f(inp.am_cell)} mm cells; compare with the cured line width.`, dStrut * 1e6 >= lw),
      kpi('layers', 'Number of layers', layers, '', `Envelope height ${f(height * 1e3)} mm divided by the layer thickness.`),
      kpi('tbuild', 'Build time', tB / 3600, 'h', `Per layer: exposed area × exposure / laser power (${f(inp.sla_P)} mW) plus ${f(inp.am_recoat)} s recoating.`),
    ],
    series: [ser('wc', 'Working curve', 'exposure (mJ/cm²)', 'cure depth (µm)', [{ name: 'cure depth', x: Es, y: Es.map(e => Dp * Math.log(e / Ec)) }, { name: 'layer thickness', x: [Es[0], Es[59]], y: [layer, layer] }], 'Semi-log working curve: slope D_p, intercept E_c. The operating point must lie above the layer-thickness line.', { logx: true }),
      ser('lwE', 'Line width versus exposure', 'exposure (mJ/cm²)', 'line width (µm)', [{ name: 'cured line width', x: Es, y: Es.map(e => inp.sla_w0 * Math.sqrt(2 * Math.log(e / Ec))) }], 'Over-exposure broadens the cured line and degrades lateral resolution.', { logx: true })],
    equations: common.concat([eqn('formula', 'Jacobs working curve', 'C_d = D_p ln(E/E_c)'), eqn('formula', 'Cured line width (Gaussian beam)', 'L_w = w₀ √(2 C_d/D_p)')]),
    checks: [chk('Resin gels (E > E_c)', cure > 0, `E/E_c = ${f(Ee / Ec)}`), chk('Layers bond (cure depth > layer)', over > 0, `overcure ${f(over)} µm`), chk('Strut resolvable', dStrut * 1e6 >= lw, `${f(dStrut * 1e6)} µm vs line width ${f(lw)} µm`)],
    warnings,
    recommendations: [
      rec('high', over <= 0 ? `Raise exposure to at least ${f(Ec * Math.exp(1.2 * layer / Dp))} mJ/cm² (20 % overcure) or reduce the layer to ${f(cure / 1.2)} µm.` : `Overcure is ${f(over / layer * 100)} % of the layer; target 10–50 % — ${over > 0.5 * layer ? 'reduce exposure to limit print-through on overhangs' : 'the setting is appropriate'}.`, 'Jacobs working curve.'),
      rec('medium', `Measure D_p and E_c for each resin batch with a windowpane test (cure depth at 5–6 exposures); pigment or ceramic loading changes both.`, 'The two resin constants fully determine the model.'),
      rec('medium', 'Post-cure under UV and heat to complete conversion, and compensate the 2–5 % linear shrinkage of acrylates in the CAD model.', 'Typical photopolymer behaviour; not computed here.'),
      rec('low', 'Filter and reuse uncured resin from the vat, wash parts in a closed IPA loop with distillation recovery, and fully cure all resin waste before disposal.', 'Good practice for vat photopolymerisation.', true),
    ],
    publish: { 'geom.strut': dStrut, 'geom.rel_density': rd, 'geom.cell_size': cell, 'proc.time': tB, 'fab.min_feature_m': lw * 1e-6, 'geom.thickness': layer * 1e-6, 'fab.yield': over > 0 ? 1 : 0 },
    exports: [{ name: 'sla_recipe.json', mime: 'application/json', content: JSON.stringify({ process: 'stereolithography', Dp_um: Dp, Ec_mJ_cm2: Ec, exposure_mJ_cm2: Ee, layer_um: layer, predicted: { cure_depth_um: cure, overcure_um: over, line_width_um: lw, build_time_h: tB / 3600 } }, null, 2), description: 'SLA exposure settings and predicted cure geometry.' }],
    numerics: { method: 'Closed-form Jacobs working-curve model', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Beer–Lambert absorption, Gaussian beam.' },
  });
}

// ================================================================================================
// 7. Sol–gel, spin coating, drying and sintering
// ================================================================================================
const ALK = { TEOS: { f: 'Si(OC2H5)4', ox: 'SiO2', alc: 'C2H5OH', hyd: 'Si(OH)4', rho: 2200, name: 'tetraethyl orthosilicate (TEOS)' }, TTIP: { f: 'Ti(OC3H7)4', ox: 'TiO2', alc: 'C3H7OH', hyd: 'Ti(OH)4', rho: 3900, name: 'titanium isopropoxide (TTIP)' } };
function runSolGel(inp) {
  const a = ALK[inp.sg_alk] || ALK.TEOS, warnings = [], Tg = inp.sg_T + 273.15, tGel = inp.sg_tref * Math.exp(inp.sg_Ea * 1e3 / R * (1 / Tg - 1 / T_REF));
  const rho0 = N.clamp(inp.sg_rho0, 0.2, 0.95), x0 = N.clamp(inp.sg_c * 1e3 * MM(a.ox) / (a.rho * rho0), 1e-4, 0.6), om = inp.sp_rpm * 2 * PI / 60, eta = inp.sp_visc * 1e-3, ev = inp.sp_evap * Math.sqrt(om), rhoS = 900;
  const hSpin = w => x0 * Math.cbrt(3 * eta * inp.sp_evap * Math.sqrt(w) / (2 * (1 - x0) * rhoS * w * w)), hGel = hSpin(om);
  const Pc = 2 * inp.dr_gamma / (inp.dr_r * 1e-9), hc = (inp.dr_KIc * 1e6 / Pc) ** 2 / 1.976, crack = hGel > hc;
  const Th = inp.sint_T + 273.15, Tr = inp.sint_Tref + 273.15, ramp = inp.sint_ramp / 60, tRamp = Math.max(0, (Th - T_REF) / ramp), tEnd = tRamp + inp.sint_hold * 60, G0 = inp.sint_G0 * 1e-9, Q = inp.sint_Q * 1e3, Qg = inp.sint_Qg * 1e3;
  const Tt = t => Math.min(Th, T_REF + ramp * t), nst = 400;
  const sol = N.rk4((t, y) => { const T = Tt(t), G3 = Math.max(y[1], G0 ** 3); return [inp.sint_k * (G0 ** 3 / G3) * Math.exp(-Q / R * (1 / T - 1 / Tr)) * (1 - Math.min(y[0], 1)), inp.sint_Kg * 1e-27 * Math.exp(-Qg / R * (1 / T - 1 / Tr)), Math.exp(-Q / (R * T)) / T]; }, [rho0, G0 ** 3, 0], 0, tEnd, nst);
  const yE = sol.y[sol.y.length - 1], rhoF = N.clamp(yE[0], rho0, 1), GF = Math.cbrt(yE[1]), theta = yE[2], hF = hGel * rho0 / rhoF, lin = 1 - Math.cbrt(rho0 / rhoF);
  const rpm = lsp(Math.log10(500), Math.log10(8000), 50), tm = sol.t.map(v => v / 60);
  let E3 = []; E3.push(eqn('chemical', 'Hydrolysis', safeEq([a.f, 'H2O'], [a.hyd, a.alc])), eqn('chemical', 'Condensation', safeEq([a.hyd], [a.ox, 'H2O'])), eqn('chemical', 'Overall sol–gel reaction', safeEq([a.f, 'H2O'], [a.ox, a.alc])));
  if (crack) warnings.push(`Dried gel film ${f(hGel * 1e9)} nm exceeds the critical cracking thickness ${f(hc * 1e9)} nm: expect channel (mud) cracks. Deposit thinner layers and build up by multiple coats.`);
  if (rhoF < 0.92) warnings.push(`Final relative density ${f(rhoF)}: pores remain open (closure occurs near 0.92), so the film is permeable.`);
  warnings.push('Gel time is an Arrhenius extrapolation of the entered reference value (pH, water ratio and catalyst are not modelled). Densification is a first-order pore-closure law with cubic grain-size dependence and parabolic-type grain growth; rate constants are user calibrations.');
  warnings.push('Meyerhofer spin-coating model assumes a Newtonian sol with constant viscosity and evaporation ∝ √ω; it does not capture edge beads or striations.');
  return finish({
    summary: `${a.name} sol (${f(inp.sg_c)} mol/L) gels in ${f(tGel)} h at ${f(inp.sg_T)} °C. Spinning at ${f(inp.sp_rpm)} rpm leaves a ${f(hGel * 1e9)} nm gel film (${crack ? 'ABOVE' : 'below'} the ${f(hc * 1e9)} nm cracking limit); sintering at ${f(inp.sint_T)} °C for ${f(inp.sint_hold)} min densifies it from ${f(rho0)} to ${f(rhoF)} relative density, giving a ${f(hF * 1e9)} nm film with ${f(GF * 1e9)} nm grains.`,
    kpis: [
      kpi('tgel', 'Gel time', tGel, 'h', `Arrhenius scaling of the reference gel time (${f(inp.sg_tref)} h at 25 °C) with apparent activation energy ${f(inp.sg_Ea)} kJ/mol; the sol must be coated well before this (pot life).`),
      kpi('hgel', 'As-spun dried gel thickness', hGel * 1e9, 'nm', `Meyerhofer: x₀[3ηe/(2(1−x₀)ρω²)]^⅓ with gel solids fraction x₀ = ${f(x0)} and evaporation rate ${f(ev)} m/s. Thickness scales as ω^(−½).`),
      kpi('hc', 'Critical cracking thickness', hc * 1e9, 'nm', `Channel-cracking criterion h_c = (K_Ic/σ)²/Z with Z = 1.976 and drying stress taken as the capillary pressure 2γ/r = ${f(Pc / 1e6)} MPa (upper bound).`, !crack),
      kpi('rhoF', 'Final relative density', rhoF, '–', `Integrated densification over the ${f(inp.sint_ramp)} K/min ramp and ${f(inp.sint_hold)} min hold; starting value ${f(rho0)}.`, rhoF > 0.95),
      kpi('hF', 'Sintered film thickness', hF * 1e9, 'nm', 'A film constrained by the substrate shrinks only through its thickness: h = h_gel ρ₀/ρ.'),
      kpi('shrink', 'Free linear shrinkage (bulk part)', lin * 100, '%', 'Isotropic linear shrinkage 1 − (ρ₀/ρ)^⅓ for an unconstrained body of the same material; use it to scale moulds.'),
      kpi('grain', 'Final grain size', GF * 1e9, 'nm', `Cubic grain-growth law G³ − G₀³ = ∫K dt from ${f(inp.sint_G0)} nm; coarsening slows densification through the G⁻³ factor.`),
      kpi('theta', 'Master-sintering-curve work log₁₀Θ', Math.log10(Math.max(theta, 1e-300)), 'log₁₀(s/K)', 'Θ = ∫(1/T)exp(−Q/RT)dt; schedules with equal Θ give the same density, so it can be used to trade temperature against time.'),
    ],
    series: [ser('sint', 'Densification and grain growth', 'time (min)', 'relative density (–) / grain size (µm)', [{ name: 'relative density', x: tm, y: sol.y.map(y => Math.min(y[0], 1)) }, { name: 'grain size (µm)', x: tm, y: sol.y.map(y => Math.cbrt(y[1]) * 1e6) }, { name: 'temperature / 1000 (°C)', x: tm, y: sol.t.map(t => (Tt(t) - 273.15) / 1000) }], 'Density rises during the ramp and hold while grains coarsen; densification stalls if grains grow before pores close.'),
      ser('spin', 'Film thickness versus spin speed', 'spin speed (rpm)', 'gel thickness (nm)', [{ name: 'Meyerhofer thickness', x: rpm, y: rpm.map(r => hSpin(r * 2 * PI / 60) * 1e9) }, { name: 'cracking limit', x: [rpm[0], rpm[49]], y: [hc * 1e9, hc * 1e9] }], 'Thickness falls as the inverse square root of spin speed; stay below the cracking limit line for single coats.', { logx: true, logy: true })],
    equations: E3.concat([eqn('formula', 'Gel time', 't_gel = t_ref exp[E_a/R (1/T − 1/T_ref)]'), eqn('formula', 'Meyerhofer spin coating', 'h = x₀ [3 η e /(2 (1 − x₀) ρ ω²)]^⅓,   e = C √ω'), eqn('formula', 'Critical cracking thickness', 'h_c = (K_Ic/σ)² / Z,   σ ≈ 2γ/r_pore,   Z = 1.976'),
      eqn('governing', 'Densification and grain growth', 'dρ/dt = k_ref (G₀/G)³ exp[−Q/R (1/T − 1/T_ref)] (1 − ρ);   dG³/dt = K_ref exp[−Q_g/R (1/T − 1/T_ref)]'), eqn('formula', 'Shrinkage', 'ΔL/L₀ = 1 − (ρ₀/ρ)^⅓;   constrained film: h = h₀ ρ₀/ρ')]),
    checks: [chk('Film below critical cracking thickness', !crack, `${f(hGel * 1e9)} nm vs ${f(hc * 1e9)} nm`), chk('Density bounded (ρ ≤ 1) and monotonic', rhoF <= 1 && rhoF >= rho0, `ρ: ${f(rho0)} → ${f(rhoF)}`), chk('Closed porosity reached (ρ ≥ 0.92)', rhoF >= 0.92, `ρ = ${f(rhoF)}`), chk('Solid mass conserved in the film', true, `h·ρ = ${f(hGel * rho0 * 1e9)} = ${f(hF * rhoF * 1e9)} nm-equivalent of dense oxide`)],
    warnings,
    recommendations: [
      rec('high', crack ? `Reduce single-coat thickness below ${f(hc * 1e9)} nm: spin at ≥ ${f(inp.sp_rpm * (hGel / hc) ** 2)} rpm or dilute the sol to ≤ ${f(inp.sg_c * hc / hGel)} mol/L, and build thickness with ${Math.ceil(hGel / hc)} or more coats with intermediate bakes.` : `Single coats of ${f(hGel * 1e9)} nm are below the cracking limit (${f(hc * 1e9)} nm); up to ${Math.max(1, Math.floor(hc / hGel))}× thicker coats remain safe.`, 'Channel-cracking criterion under capillary drying stress.'),
      rec('medium', rhoF < 0.95 ? `Density reaches only ${f(rhoF)}: raise the hold temperature by ${f(R * Th * Th / Q * Math.log(3))} K (≈ 3× faster densification) or extend the hold; avoid much higher temperature where grain growth (now ${f(GF * 1e9)} nm) takes over.` : `Densification is essentially complete (${f(rhoF)}); shorten the hold or lower the temperature to limit grain growth (${f(GF * 1e9)} nm) and energy use.`, 'Integrated densification / grain-growth kinetics.'),
      rec('medium', `Coat within about ${f(tGel * 0.3)} h of mixing (30 % of gel time) and control humidity during spinning; viscosity rises steeply near gelation.`, 'Gel-time estimate.'),
      rec('medium', `Use ethanol or isopropanol rather than methoxyethanol as solvent, recover the alcohol by-product (${chem.prettyFormula(a.alc)}) from the bake exhaust, and prefer rapid thermal or two-step sintering to cut furnace energy for the ${f((tEnd / 3600))} h cycle.`, 'Solvent hazard ranking and furnace cycle time of this run.', true),
    ],
    publish: { 'geom.thickness': hF, 'geom.porosity': 1 - rhoF, 'geom.rel_density': rhoF, 'proc.T_C': inp.sint_T, 'proc.time': tEnd, 'fab.shrinkage': lin, 'geom.grain_size': GF, 'material.formula': a.ox },
    exports: [{ name: 'solgel_recipe.json', mime: 'application/json', content: JSON.stringify({ precursor: a.name, concentration_mol_L: inp.sg_c, gel_time_h: tGel, spin_rpm: inp.sp_rpm, gel_thickness_nm: hGel * 1e9, critical_thickness_nm: hc * 1e9, sinter: { T_C: inp.sint_T, ramp_K_min: inp.sint_ramp, hold_min: inp.sint_hold }, predicted: { density: rhoF, thickness_nm: hF * 1e9, grain_nm: GF * 1e9 } }, null, 2), description: 'Sol–gel coating and sintering recipe with predicted film properties.' }],
    numerics: { method: 'Closed-form Meyerhofer and cracking relations; RK4 integration of densification, grain growth and master-curve work', dof: 3, steps: nst, iterations: nst, residual: 0, notes: 'Fixed-step RK4 over ramp + hold.' },
  });
}

// ================================================================================================
// 8. Electrospinning and colloidal batch recipe
// ================================================================================================
function runSpinColloid(inp) {
  const warnings = [], gam = inp.es_gamma, Hc = inp.es_H, Ln = inp.es_L, Rn = inp.es_R * 0.1, lg = Math.log(2 * Ln / Rn) - 1.5, Vc = lg > 0 ? Math.sqrt(4 * Hc * Hc / (Ln * Ln) * lg * 0.117 * PI * Rn * gam) : 0;
  const Qf = inp.es_Q * 1e-6 / 3600, I = inp.es_I * 1e-9, cv = inp.es_c / 100, chi = 100, ht = Math.cbrt(gam * 1e-3 * N.EPS0 * Qf * Qf / (I * I) * 2 / (PI * (2 * Math.log(chi) - 3))), df = ht * Math.sqrt(cv), prod = Qf * cv * inp.es_rho * 3600 * inp.es_needles, lenRate = Qf * cv / (PI / 4 * df * df);
  if (lg <= 0) warnings.push('Needle length is too short relative to its radius for the Taylor critical-voltage formula (ln(2L/R) ≤ 1.5).');
  if (inp.es_V < Vc) warnings.push(`Applied voltage ${f(inp.es_V)} kV is below the Taylor critical voltage ${f(Vc)} kV: no stable jet (dripping).`);
  if (inp.es_V > 2.5 * Vc && Vc > 0) warnings.push('Voltage is more than 2.5× the critical value: multiple jets and corona discharge are likely.');
  warnings.push('Terminal jet diameter is the Fridrikh-type balance of surface tension and surface-charge repulsion; it ignores viscoelastic stresses and solvent evaporation during whipping, so it is a lower bound — measured fibres are usually 1–5× thicker. Jet current must be measured.');
  // ---- colloidal precipitation batch ----
  const saltF = SALTS[inp.salt] ? inp.salt : 'Zn(NO3)2', salt = SALTS[saltF], baseF = BASES[inp.base] ? inp.base : 'NaOH', base = BASES[baseF], VL = inp.col_V, nZn = inp.cZn * VL * 1e-3, X = N.clamp(inp.col_conv, 0, 1);
  let coef = [1, 2, 1, 1, 2], prodSalt = salt.prod[base.cat], ie = null;
  try { ie = base.hmta ? chem.ionicEquations([{ f: saltF, phase: 'aq' }, { f: 'NH3', phase: 'aq' }, { f: 'H2O', phase: 'l' }], [{ f: 'ZnO', phase: 's' }, { f: prodSalt, phase: 'aq' }]) : chem.ionicEquations([{ f: saltF, phase: 'aq' }, { f: baseF, phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'H2O', phase: 'l' }, { f: prodSalt, phase: 'aq' }]); coef = ie.coefficients; } catch (e) { warnings.push('Could not balance the precipitation reaction automatically: ' + e.message); }
  const nSaltBy = coef[coef.length - 1] / coef[0] * nZn * X, nBase = coef[1] / coef[0] * nZn * (base.hmta ? 0.25 : 1), mZnO = nZn * X * MM('ZnO'), nW = Math.round(inp.col_wash), wr = inp.col_wratio, loss = inp.col_loss / 100;
  const resid = (1 / (1 + wr)) ** nW, mOut = mZnO * (1 - loss) ** nW, yld = nZn > 0 ? mOut / (nZn * MM('ZnO')) : 0, Vwash = nW * wr * inp.col_ret * mZnO, solv = mOut > 0 ? (VL + Vwash) / mOut : 0, mBy = nSaltBy * MM(prodSalt), mSaltIn = nZn * MM(saltF), mBaseIn = nBase * MM(base.hmta ? 'C6H12N4' : baseF);
  const nws = N.linspace(0, 8, 9);
  const E = [eqn('formula', 'Taylor critical voltage', 'V_c² = (4H²/L²) [ln(2L/R) − 3/2] (0.117 π R γ)   (kV, cm, dyn/cm)'), eqn('formula', 'Terminal jet diameter', 'h_t = [γ ε₀ (Q/I)² · 2 /(π (2 ln χ − 3))]^⅓,   d_fibre = h_t √c_polymer', 'χ ≈ 100 is the whipping-to-jet length-scale ratio.')];
  if (ie) E.push(eqn('chemical', 'Batch precipitation (molecular)', ie.molecular), eqn('ionic', 'Complete ionic equation', ie.completeIonic), eqn('ionic', 'Net ionic equation', ie.netIonic));
  if (base.hmta) E.push(eqn('chemical', 'HMTA hydrolysis', safeEq(['C6H12N4', 'H2O'], ['HCHO', 'NH3'])));
  E.push(eqn('formula', 'Dilution washing', 'residual impurity fraction = (1/(1 + w))ⁿ;   yield = X (1 − ℓ)ⁿ'));
  return finish({
    summary: `Electrospinning: the Taylor-cone onset is ${f(Vc)} kV (applied ${f(inp.es_V)} kV); the terminal jet is ${f(ht * 1e9)} nm and the dry fibre about ${f(df * 1e9)} nm, produced at ${f(prod * 1e3)} g/h from ${f(inp.es_needles)} needle(s). Colloidal batch: ${f(VL)} L at ${f(inp.cZn)} mol/m³ gives ${f(mOut * 1e3)} g of ZnO after ${nW} washes (overall yield ${f(yld * 100)} %), using ${f(solv)} L of water per kg and leaving ${f(resid * 100)} % of the dissolved salt.`,
    kpis: [
      kpi('Vc', 'Taylor critical voltage', Vc, 'kV', `Onset of the Taylor cone for a ${f(Ln)} cm needle of ${f(inp.es_R)} mm radius at ${f(Hc)} cm gap and γ = ${f(gam)} mN/m; operate about 1.2–2× above it.`, inp.es_V >= Vc && inp.es_V <= 2.5 * Vc),
      kpi('jet', 'Terminal jet diameter', ht * 1e9, 'nm', `Balance of surface tension and charge repulsion at flow ${f(inp.es_Q)} mL/h and current ${f(inp.es_I)} nA; it scales as (Q/I)^⅔.`),
      kpi('fibre', 'Dry fibre diameter', df * 1e9, 'nm', `Jet diameter × √(polymer volume fraction ${f(cv)}) after solvent loss; a lower bound (see warnings).`),
      kpi('prod', 'Fibre production rate', prod * 1e3, 'g/h', `Flow rate × polymer fraction × density × ${f(inp.es_needles)} needle(s); single needles give 0.01–1 g/h, hence multi-needle or needleless heads for scale-up.`),
      kpi('len', 'Fibre length rate per needle', lenRate / 1e3, 'km/s', 'Volume conservation: polymer flow divided by fibre cross-section; illustrates the extreme draw ratio of the whipping jet.'),
      kpi('mZnO', 'Batch product mass after washing', mOut * 1e3, 'g', `${f(nZn)} mol Zn × conversion ${f(X)} × M(ZnO), less ${f(inp.col_loss)} % solids loss per wash.`),
      kpi('yield', 'Overall batch yield', yld * 100, '%', 'Product recovered divided by the stoichiometric maximum from the zinc charged.', yld > 0.8),
      kpi('resid', 'Residual dissolved salt after washing', resid * 100, '%', `Dilution washing with ${f(wr)} volumes of water per volume of retained liquor, ${nW} times.`, resid < 0.01),
      kpi('solvent', 'Water use', solv, 'L/kg', 'Batch volume plus wash water per kg of product recovered; the main waste stream of precipitation routes.'),
      kpi('byprod', 'By-product salt', mBy * 1e3, 'g', `${chem.prettyFormula(prodSalt)} formed stoichiometrically; it leaves in the mother liquor and washings.`),
    ],
    series: [ser('wash', 'Washing efficiency', 'number of washes', 'fraction (%)', [{ name: 'residual salt', x: nws, y: nws.map(n => Math.max(1e-10, (1 / (1 + wr)) ** n * 100)) }, { name: 'product yield', x: nws, y: nws.map(n => X * (1 - loss) ** n * 100) }], 'Each wash lowers the dissolved impurity geometrically but costs product and water; choose the fewest washes that meet the purity target.', { logy: true }),
      ser('jetQ', 'Fibre diameter versus flow rate', 'flow rate (mL/h)', 'diameter (nm)', (() => { const qs = lsp(-1.5, 1.5, 50); return [{ name: 'dry fibre (constant current)', x: qs, y: qs.map(q => df * (q / inp.es_Q) ** (2 / 3) * 1e9) }]; })(), 'At constant jet current the fibre diameter scales as Q^⅔; in practice current rises with flow, weakening the dependence.', { logx: true, logy: true })],
    equations: E,
    checks: [chk('Voltage above Taylor onset', inp.es_V >= Vc, `${f(inp.es_V)} kV vs ${f(Vc)} kV`), chk('Zinc mass balance of the batch', true, `Zn in ${f(nZn)} mol = ZnO ${f(nZn * X)} mol + unreacted ${f(nZn * (1 - X))} mol`), chk('Purity target (< 1 % residual salt)', resid < 0.01, `${f(resid * 100)} %`)],
    warnings,
    recommendations: [
      rec('high', inp.es_V < Vc ? `Raise the voltage to about ${f(Vc * 1.3)} kV or shorten the gap to start a stable jet.` : `Operate at ${f(Vc * 1.2)}–${f(Vc * 2)} kV; to thin the fibres below ${f(df * 1e9)} nm lower the flow rate or raise solution conductivity (salt additive) to increase jet current.`, 'Taylor onset and (Q/I)^⅔ scaling.'),
      rec('medium', resid >= 0.01 ? `Add ${Math.max(1, Math.ceil(Math.log(0.01) / Math.log(1 / (1 + wr))) - nW)} more wash(es) or raise the wash ratio to reach < 1 % residual ${chem.prettyFormula(prodSalt)}.` : `Washing meets < 1 % residual salt; ${nW > 1 && (1 / (1 + wr)) ** (nW - 1) < 0.01 ? 'one wash fewer would still meet it and save water and product.' : 'do not reduce the number of washes.'}`, 'Geometric dilution-washing model.'),
      rec('medium', `Charge ${f(mSaltIn * 1e3)} g ${chem.prettyFormula(saltF)} and ${f(mBaseIn * 1e3)} g ${base.hmta ? 'HMTA' : chem.prettyFormula(baseF)} for the ${f(VL)} L batch; add base slowly under vigorous mixing to avoid local over-saturation and broad size distributions.`, 'Stoichiometry of the balanced reaction.'),
      rec('medium', `Water use is ${f(solv)} L/kg: use counter-current washing or diafiltration to cut it by half or more, recover the ${chem.prettyFormula(prodSalt)} by-product (${f(mBy * 1e3)} g per batch) and electrospin from benign solvents (water, ethanol, acetic acid) in place of DMF or chloroform.`, 'Solvent and by-product inventory of this run.', true),
    ],
    publish: { 'geom.diameter': df, 'fab.throughput_kg_h': prod, 'rxn.yield_kg': mOut, 'rxn.conversion': X, 'fab.yield': N.clamp(yld, 0, 1), 'proc.volume': VL * 1e-3, 'proc.c_precursor': inp.cZn, 'fab.solvent_L_kg': solv },
    exports: [{ name: 'electrospin_colloid_recipe.json', mime: 'application/json', content: JSON.stringify({ electrospinning: { voltage_kV: inp.es_V, critical_kV: Vc, gap_cm: Hc, flow_mL_h: inp.es_Q, predicted_fibre_nm: df * 1e9, production_g_h: prod * 1e3 }, colloid_batch: { volume_L: VL, zinc_salt: saltF, zinc_salt_g: mSaltIn * 1e3, base: baseF, base_g: mBaseIn * 1e3, washes: nW, product_g: mOut * 1e3, yield: yld, water_L_per_kg: solv, byproduct: prodSalt, byproduct_g: mBy * 1e3 } }, null, 2), description: 'Electrospinning settings and precipitation batch sheet (charges, washes, yield, by-product).' }],
    numerics: { method: 'Closed-form Taylor and terminal-jet scaling relations; stoichiometric batch balance with geometric washing', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Reaction balanced with the shared chemistry kernel.' },
  });
}

const MI = (mode, group) => o => ({ group, showIf: { mode: [mode] }, ...o });
{
  const C = MI('cvd', 'CVD'), A_ = MI('ald', 'ALD'), V = MI('pvd', 'PVD'), L = MI('litho_etch', 'Lithography'), ET = MI('litho_etch', 'Etch'), M = MI('additive', 'Additive manufacturing'), G = MI('solgel_spin_sinter', 'Sol–gel and spin coating'), SN = MI('solgel_spin_sinter', 'Sintering'), ES = MI('electrospin_colloid', 'Electrospinning'), CO = MI('electrospin_colloid', 'Colloidal batch');
  const opt = (o, lab) => Object.keys(o).map(k => ({ value: k, label: lab(o[k], k) }));
  INPUTS.push(
    C({ key: 'cvd_chem', label: 'Chemistry', type: 'select', value: 'sih4', options: opt(CVDC, v => v.name), help: 'Deposition reaction; sets film density, typical surface kinetics, gas diffusivity and carrier gas.' }),
    C({ key: 'cvd_kind', label: 'Activation', type: 'select', value: 'thermal', options: [{ value: 'thermal', label: 'Thermal CVD (Arrhenius surface reaction)' }, { value: 'pecvd', label: 'PECVD (radical flux with sticking probability)' }], help: 'Thermal CVD uses the Arrhenius surface rate; PECVD replaces it by plasma-generated radicals arriving with a sticking probability.' }),
    C({ key: 'cvd_T', label: 'Substrate temperature', unit: '°C', value: 620, min: 100, max: 1300, bus: 'proc.T_C', help: 'Controls the surface rate exponentially and the gas properties weakly.' }),
    C({ key: 'cvd_P', label: 'Total pressure', unit: 'Pa', value: 100, min: 1, max: 101325, log: true, help: 'Low pressure raises diffusivity (D ∝ 1/P) and pushes the process to reaction-limited, uniform growth.' }),
    C({ key: 'cvd_x', label: 'Precursor mole fraction', unit: '–', value: 0.3, min: 1e-4, max: 1, log: true, help: 'Mole fraction of the limiting precursor in the feed; growth rate is proportional to it in the first-order model.' }),
    C({ key: 'cvd_U', label: 'Gas velocity', unit: 'm/s', value: 1, min: 0.01, max: 30, log: true, help: 'Mean velocity over the susceptor; thins the boundary layer (h_g ∝ √U) and reduces depletion.' }),
    C({ key: 'cvd_L', label: 'Susceptor length', unit: 'm', value: 0.2, min: 0.01, max: 2, log: true, help: 'Length in the flow direction over which boundary-layer growth and depletion act.' }),
    C({ key: 'cvd_H', label: 'Channel height', unit: 'm', value: 0.05, min: 0.002, max: 0.5, log: true, help: 'Free height above the susceptor; sets the gas inventory per unit area and the fully developed mass-transfer floor.' }),
    C({ key: 'cvd_t', label: 'Deposition time', unit: 's', value: 1800, min: 1, max: 86400, log: true, help: 'Film thickness is rate × time.' }),
    C({ key: 'cvd_AR', label: 'Trench aspect ratio', unit: '–', value: 5, min: 0.1, max: 100, log: true, help: 'Depth/width of the feature for the step-coverage estimate.' }),
    C({ key: 'cvd_s', label: 'Radical sticking probability (PECVD)', unit: '–', value: 0.1, min: 1e-4, max: 1, log: true, showIf: { mode: ['cvd'], cvd_kind: ['pecvd'] }, help: 'Probability that an arriving radical reacts; high values give fast but non-conformal growth.' }),
    A_({ key: 'ald_chem', label: 'Chemistry', type: 'select', value: 'tma', options: opt(ALDC, v => v.name), help: 'ALD process; sets the saturated metal uptake per cycle, film density and typical temperature window.' }),
    A_({ key: 'ald_T', label: 'Substrate temperature', unit: '°C', value: 200, min: 25, max: 500, bus: 'proc.T_C', help: 'Must lie inside the ALD window; affects impingement flux and the Gordon exposure only weakly (√T).' }),
    A_({ key: 'ald_P', label: 'Precursor partial pressure', unit: 'Pa', value: 10, min: 0.01, max: 1000, log: true, help: 'Partial pressure during the pulse; exposure is pressure × time.' }),
    A_({ key: 'ald_tp', label: 'Pulse time', unit: 's', value: 0.5, min: 0.005, max: 300, log: true, help: 'Duration of each precursor pulse; must reach saturation on flat surfaces and deliver the Gordon exposure for deep holes.' }),
    A_({ key: 'ald_thk', label: 'Target thickness', unit: 'nm', value: 20, min: 0.5, max: 500, log: true, help: 'Film thickness required; sets the cycle count.' }),
    A_({ key: 'ald_AR', label: 'Hole aspect ratio', unit: '–', value: 50, min: 1, max: 2000, log: true, help: 'Depth/diameter of the holes or pores to coat (use half the value for through-holes open at both ends).' }),
    A_({ key: 'ald_dhole', label: 'Hole diameter', unit: 'nm', value: 100, min: 5, max: 1e5, log: true, help: 'Sets the Knudsen diffusivity for purge-out from the holes.' }),
    A_({ key: 'ald_V', label: 'Chamber volume', unit: 'L', value: 1, min: 0.05, max: 500, log: true, help: 'With pumping speed gives the gas residence time that controls purge duration and precursor dose.' }),
    A_({ key: 'ald_S', label: 'Effective pumping speed', unit: 'L/s', value: 10, min: 0.1, max: 2000, log: true, help: 'Volumetric pumping (or purge-flow) speed at the chamber.' }),
    A_({ key: 'ald_area', label: 'Total surface area to coat', unit: 'cm²', value: 2000, min: 1, max: 1e7, log: true, help: 'Substrate plus internal surface (high-aspect structures multiply the area); used for precursor utilisation.' }),
    V({ key: 'pvd_proc', label: 'Source', type: 'select', value: 'sputter', options: [{ value: 'sputter', label: 'DC magnetron sputtering' }, { value: 'evap', label: 'Thermal / e-beam evaporation' }], help: 'Sputtering uses the ion-yield source term; evaporation the Hertz–Knudsen flux.' }),
    V({ key: 'pvd_mat', label: 'Target / source material', type: 'select', value: 'Cu', options: Object.keys(PVDM).map(k => ({ value: k, label: k })), help: 'Sets sputter yield, vapour pressure, density and film elastic/thermal constants.' }),
    V({ key: 'pvd_power', label: 'Sputter power', unit: 'W', value: 300, min: 10, max: 20000, log: true, bus: 'proc.power', showIf: { mode: ['pvd'], pvd_proc: ['sputter'] }, help: 'Discharge power; ion current = power/voltage.' }),
    V({ key: 'pvd_V', label: 'Target voltage', unit: 'V', value: 400, min: 100, max: 1500, showIf: { mode: ['pvd'], pvd_proc: ['sputter'] }, help: 'Sets the ion energy and hence the sputter yield (∝ √E).' }),
    V({ key: 'pvd_Tsrc', label: 'Evaporation source temperature', unit: '°C', value: 1300, min: 300, max: 3500, showIf: { mode: ['pvd'], pvd_proc: ['evap'] }, help: 'Melt-surface temperature; vapour pressure and rate rise exponentially with it.' }),
    V({ key: 'pvd_Asrc', label: 'Evaporating area', unit: 'cm²', value: 1, min: 0.01, max: 100, log: true, showIf: { mode: ['pvd'], pvd_proc: ['evap'] }, help: 'Area of the molten pool that emits vapour.' }),
    V({ key: 'pvd_P', label: 'Working pressure', unit: 'Pa', value: 0.5, min: 1e-5, max: 20, log: true, help: 'Sets the mean free path; sputtering runs at 0.1–2 Pa, evaporation below 10⁻³ Pa.' }),
    V({ key: 'pvd_h', label: 'Throw distance', unit: 'cm', value: 10, min: 2, max: 100, log: true, help: 'Source-to-substrate distance: rate falls about as 1/h² while uniformity improves.' }),
    V({ key: 'pvd_a', label: 'Racetrack (ring-source) radius', unit: 'cm', value: 2.5, min: 0, max: 30, showIf: { mode: ['pvd'], pvd_proc: ['sputter'] }, help: 'Mean radius of the magnetron erosion groove; zero gives a point source.' }),
    V({ key: 'pvd_wafer', label: 'Wafer diameter', unit: 'mm', value: 100, min: 5, max: 450, help: 'Substrate over which thickness uniformity is evaluated.' }),
    V({ key: 'pvd_t', label: 'Deposition time', unit: 's', value: 600, min: 1, max: 86400, log: true, help: 'Thickness is rate × time.' }),
    V({ key: 'pvd_Tdep', label: 'Substrate temperature during deposition', unit: '°C', value: 100, min: 20, max: 800, help: 'Thermal-mismatch stress develops on cooling from this temperature to 25 °C.' }),
    V({ key: 'pvd_sint', label: 'Intrinsic (growth) stress', unit: 'MPa', value: 0, min: -3000, max: 3000, help: 'Measured or assumed growth stress (negative = compressive) added to the thermal stress.' }),
    V({ key: 'pvd_tsub', label: 'Silicon substrate thickness', unit: 'µm', value: 525, min: 50, max: 2000, help: 'Used in Stoney’s equation for wafer bow.' }),
    L({ key: 'lit_tool', label: 'Patterning method', type: 'select', value: 'optical', options: [{ value: 'optical', label: 'Optical / EUV projection lithography' }, { value: 'ebeam', label: 'Electron-beam lithography' }, { value: 'nil', label: 'Nanoimprint lithography' }], help: 'Selects the image-formation model used to predict the printed CD before etch transfer.' }),
    L({ key: 'lit_cd', label: 'Design feature size (CD)', unit: 'nm', value: 70, min: 5, max: 5000, log: true, help: 'Target line width of the periodic pattern (e.g. a metasurface element).' }),
    L({ key: 'lit_pitch', label: 'Pitch', unit: 'nm', value: 140, min: 10, max: 20000, log: true, help: 'Period of the array; must exceed the feature size.' }),
    L({ key: 'lit_tol', label: 'CD tolerance', unit: '± %', value: 10, min: 1, max: 30, help: 'Allowed deviation of the final dimension; defines the process window.' }),
    L({ key: 'lit_src', label: 'Exposure source', type: 'select', value: 'arfi', options: opt(LSRC, v => v.name), showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Exposure wavelength (and immersion medium).' }),
    L({ key: 'lit_NA', label: 'Numerical aperture', unit: '–', value: 1.35, min: 0.1, max: 1.4, showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Projection-lens NA; resolution improves as 1/NA, depth of focus degrades as 1/NA².' }),
    L({ key: 'lit_sigma', label: 'Partial coherence σ', unit: '–', value: 0.8, min: 0, max: 0.98, showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Illumination fill; extends the capture condition of the first diffraction order to NA(1+σ).' }),
    L({ key: 'lit_blur', label: 'Resist blur (acid diffusion length)', unit: 'nm', value: 5, min: 0, max: 100, showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Gaussian blur of the latent image by acid diffusion and secondary electrons; it lowers contrast for small pitches.' }),
    L({ key: 'lit_k1', label: 'Process factor k₁', unit: '–', value: 0.3, min: 0.25, max: 1, showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Rayleigh resolution factor of the process (0.25 theoretical limit, 0.3–0.4 aggressive, > 0.6 relaxed).' }),
    L({ key: 'lit_k2', label: 'Focus factor k₂', unit: '–', value: 0.5, min: 0.2, max: 1.5, showIf: { mode: ['litho_etch'], lit_tool: ['optical'] }, help: 'Rayleigh depth-of-focus factor.' }),
    L({ key: 'eb_kV', label: 'Beam energy', unit: 'kV', value: 100, min: 5, max: 150, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Higher energy spreads backscattered electrons over a larger range β (a lower, flatter background).' }),
    L({ key: 'eb_dose', label: 'Exposure dose', unit: 'µC/cm²', value: 600, min: 5, max: 20000, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Areal charge delivered to the exposed shapes.' }),
    L({ key: 'eb_D0', label: 'Resist threshold dose', unit: 'µC/cm²', value: 250, min: 1, max: 10000, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Absorbed-dose threshold at which the (negative-tone) resist remains after development.' }),
    L({ key: 'eb_alpha', label: 'Forward-scattering range α', unit: 'nm', value: 8, min: 1, max: 200, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Gaussian width of the beam plus forward scattering in the resist; sets edge sharpness.' }),
    L({ key: 'eb_eta', label: 'Backscatter ratio η', unit: '–', value: 0.7, min: 0, max: 2, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Energy deposited by backscattered relative to forward electrons (≈ 0.7 on silicon, higher on heavy substrates).' }),
    L({ key: 'eb_I', label: 'Beam current', unit: 'nA', value: 1, min: 0.01, max: 200, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Sets the write time.' }),
    L({ key: 'lit_area', label: 'Patterned field area', unit: 'mm²', value: 1, min: 1e-4, max: 1e4, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['ebeam'] }, help: 'Area of the array to be written.' }),
    L({ key: 'nil_visc', label: 'Resist viscosity at imprint', unit: 'Pa·s', value: 1e4, min: 1e-3, max: 1e7, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['nil'] }, help: 'Thermal NIL polymers: 10³–10⁵ Pa·s above T_g; UV-NIL resins: about 10⁻² Pa·s.' }),
    L({ key: 'nil_p', label: 'Imprint pressure', unit: 'MPa', value: 5, min: 0.01, max: 100, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['nil'] }, help: 'Applied pressure driving squeeze flow.' }),
    L({ key: 'nil_h0', label: 'Initial resist thickness', unit: 'nm', value: 200, min: 10, max: 5000, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['nil'] }, help: 'Spin-coated resist thickness before imprint.' }),
    L({ key: 'nil_hf', label: 'Residual layer target', unit: 'nm', value: 20, min: 1, max: 1000, log: true, showIf: { mode: ['litho_etch'], lit_tool: ['nil'] }, help: 'Resist thickness left under the stamp protrusions.' }),
    ET({ key: 'et_rate', label: 'Open-area etch rate', unit: 'nm/min', value: 200, min: 1, max: 20000, log: true, help: 'Vertical etch rate of the film on unpatterned areas.' }),
    ET({ key: 'et_depth', label: 'Target etch depth', unit: 'nm', value: 100, min: 5, max: 1e5, log: true, help: 'Depth to transfer into the film or substrate.' }),
    ET({ key: 'et_sel', label: 'Selectivity to mask', unit: ':1', value: 5, min: 0.2, max: 500, log: true, help: 'Film etch rate divided by mask etch rate.' }),
    ET({ key: 'et_A', label: 'Anisotropy', unit: '–', value: 0.99, min: 0, max: 1, help: '1 − lateral/vertical rate; 1 is perfectly vertical, 0 isotropic (wet-etch-like).' }),
    ET({ key: 'et_S', label: 'Bottom reaction probability', unit: '–', value: 0.3, min: 0.01, max: 1, log: true, help: 'Probability that a neutral reaching the trench bottom reacts; with the Knudsen conductance it sets the aspect-ratio slow-down.' }),
    ET({ key: 'et_resist', label: 'Mask (resist) thickness', unit: 'nm', value: 100, min: 5, max: 20000, log: true, help: 'Available mask thickness to survive the etch.' }),
    M({ key: 'am_proc', label: 'Process', type: 'select', value: 'lpbf', options: [{ value: 'lpbf', label: 'Laser powder-bed fusion (metal)' }, { value: 'tpp', label: 'Two-photon polymerisation' }, { value: 'sla', label: 'Stereolithography (Jacobs)' }], help: 'Additive process for the lattice; each has its own feature-size and defect model.' }),
    M({ key: 'am_rd', label: 'Lattice relative density', unit: '–', value: 0.2, min: 0.01, max: 0.6, bus: 'geom.rel_density', help: 'Solid volume fraction of the lattice; with the cell size it fixes the strut diameter to be printed.' }),
    M({ key: 'am_cell', label: 'Unit-cell size', unit: 'mm', value: 2, min: 0.001, max: 50, log: true, help: 'Lattice period (use micrometre-scale values such as 0.01 mm for two-photon lattices).' }),
    M({ key: 'am_vol', label: 'Part envelope volume', unit: 'cm³', value: 10, min: 0.01, max: 5000, log: true, showIf: { mode: ['additive'], am_proc: ['lpbf', 'sla'] }, help: 'Bounding volume of the lattice part (taken as a cube for the layer count).' }),
    M({ key: 'am_recoat', label: 'Recoat time per layer', unit: 's', value: 10, min: 0, max: 120, showIf: { mode: ['additive'], am_proc: ['lpbf', 'sla'] }, help: 'Powder spreading or resin recoating time added for every layer.' }),
    M({ key: 'am_mat', label: 'Alloy', type: 'select', value: 'Ti64', options: Object.keys(AMM).map(k => ({ value: k, label: AMM[k][6] })), showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Powder alloy; sets thermal properties, melting/boiling points and absorptivity.' }),
    M({ key: 'am_P', label: 'Laser power', unit: 'W', value: 200, min: 20, max: 1000, bus: 'proc.power', showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Laser output power.' }),
    M({ key: 'am_v', label: 'Scan speed', unit: 'mm/s', value: 1000, min: 50, max: 5000, log: true, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Beam travel speed.' }),
    M({ key: 'am_h', label: 'Hatch spacing', unit: 'µm', value: 100, min: 20, max: 300, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Distance between adjacent scan tracks.' }),
    M({ key: 'am_t', label: 'Layer thickness', unit: 'µm', value: 30, min: 10, max: 150, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Powder layer thickness.' }),
    M({ key: 'am_spot', label: 'Laser spot diameter (4σ)', unit: 'µm', value: 120, min: 30, max: 500, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Beam diameter at the bed; small spots concentrate energy and promote keyholing.' }),
    M({ key: 'am_T0', label: 'Bed preheat temperature', unit: '°C', value: 100, min: 20, max: 800, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Baseplate/bed temperature; raises the melt-pool size slightly and lowers residual stress.' }),
    M({ key: 'am_powder', label: 'Powder size D50', unit: 'µm', value: 30, min: 5, max: 100, showIf: { mode: ['additive'], am_proc: ['lpbf'] }, help: 'Median particle diameter; partially fused particles add about one diameter to strut size.' }),
    M({ key: 'tpp_lam', label: 'Laser wavelength', unit: 'nm', value: 780, min: 500, max: 1100, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Femtosecond laser wavelength.' }),
    M({ key: 'tpp_NA', label: 'Objective NA', unit: '–', value: 1.4, min: 0.3, max: 1.45, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Numerical aperture; sets focal spot and voxel aspect ratio.' }),
    M({ key: 'tpp_n', label: 'Resin refractive index', unit: '–', value: 1.5, min: 1.3, max: 1.7, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Refractive index of the resin at the laser wavelength.' }),
    M({ key: 'tpp_P', label: 'Average laser power at focus', unit: 'mW', value: 15, min: 0.5, max: 100, log: true, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Average power delivered to the sample.' }),
    M({ key: 'tpp_Pth', label: 'Threshold power at 100 µm/s', unit: 'mW', value: 1, min: 0.05, max: 50, log: true, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Measured polymerisation threshold at the reference speed; scales as √v.' }),
    M({ key: 'tpp_v', label: 'Scan speed', unit: 'µm/s', value: 5000, min: 10, max: 5e5, log: true, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Writing speed (galvo systems reach 10⁴–10⁵ µm/s).' }),
    M({ key: 'tpp_side', label: 'Lattice cube side', unit: 'µm', value: 100, min: 5, max: 5000, log: true, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Edge length of the lattice block to write.' }),
    M({ key: 'tpp_shrink', label: 'Linear shrinkage', unit: '%', value: 5, min: 0, max: 30, showIf: { mode: ['additive'], am_proc: ['tpp'] }, help: 'Shrinkage on development and drying, to be pre-compensated.' }),
    M({ key: 'sla_Dp', label: 'Penetration depth D_p', unit: 'µm', value: 150, min: 5, max: 1000, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Depth at which irradiance falls to 1/e in the resin.' }),
    M({ key: 'sla_Ec', label: 'Critical exposure E_c', unit: 'mJ/cm²', value: 10, min: 0.5, max: 500, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Exposure at which the resin gels.' }),
    M({ key: 'sla_E', label: 'Applied exposure', unit: 'mJ/cm²', value: 15, min: 0.5, max: 5000, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Energy per area delivered to each layer.' }),
    M({ key: 'sla_layer', label: 'Layer thickness', unit: 'µm', value: 50, min: 5, max: 500, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Build-platform step between layers.' }),
    M({ key: 'sla_w0', label: 'Beam radius (1/e²)', unit: 'µm', value: 70, min: 5, max: 500, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Laser spot radius (or pixel half-size for DLP).' }),
    M({ key: 'sla_P', label: 'Laser power', unit: 'mW', value: 100, min: 1, max: 5000, log: true, showIf: { mode: ['additive'], am_proc: ['sla'] }, help: 'Optical power at the vat; sets exposure time per layer.' }),
    G({ key: 'sg_alk', label: 'Alkoxide precursor', type: 'select', value: 'TEOS', options: opt(ALK, v => v.name), help: 'Metal alkoxide; sets the hydrolysis–condensation equations and the oxide formed.' }),
    G({ key: 'sg_c', label: 'Alkoxide concentration', unit: 'mol/L', value: 0.5, min: 0.01, max: 3, log: true, help: 'Sets the oxide (gel solids) volume fraction of the sol and hence film thickness.' }),
    G({ key: 'sg_T', label: 'Ageing temperature', unit: '°C', value: 25, min: 0, max: 90, help: 'Gelation accelerates with temperature following the apparent activation energy.' }),
    G({ key: 'sg_tref', label: 'Gel time at 25 °C', unit: 'h', value: 100, min: 0.01, max: 5000, log: true, help: 'Measured gel time of this sol at 25 °C (depends strongly on pH and water ratio).' }),
    G({ key: 'sg_Ea', label: 'Gelation activation energy', unit: 'kJ/mol', value: 60, min: 10, max: 120, help: 'Apparent activation energy of gelation (40–80 kJ/mol typical for silica sols).' }),
    G({ key: 'sg_rho0', label: 'Gel relative density', unit: '–', value: 0.5, min: 0.2, max: 0.9, help: 'Relative density of the dried gel film (xerogel) before sintering.' }),
    G({ key: 'sp_rpm', label: 'Spin speed', unit: 'rpm', value: 3000, min: 300, max: 10000, log: true, help: 'Final spin speed; thickness ∝ ω^(−½).' }),
    G({ key: 'sp_visc', label: 'Sol viscosity', unit: 'mPa·s', value: 3, min: 0.3, max: 500, log: true, help: 'Viscosity of the sol at coating; thickness ∝ η^(⅓).' }),
    G({ key: 'sp_evap', label: 'Evaporation coefficient C', unit: 'm/s per √(rad/s)', value: 3e-8, min: 1e-9, max: 1e-6, log: true, help: 'Solvent evaporation rate e = C√ω (about 3×10⁻⁸ for ethanol-based sols in still air); calibrate from one measured thickness.' }),
    G({ key: 'dr_r', label: 'Gel pore radius', unit: 'nm', value: 5, min: 0.5, max: 500, log: true, help: 'Sets the capillary drying pressure 2γ/r.' }),
    G({ key: 'dr_gamma', label: 'Pore-liquid surface tension', unit: 'N/m', value: 0.022, min: 0.005, max: 0.08, help: 'Surface tension of the pore liquid (ethanol 0.022, water 0.072 N/m).' }),
    G({ key: 'dr_KIc', label: 'Gel fracture toughness', unit: 'MPa·√m', value: 0.012, min: 0.001, max: 1, log: true, help: 'Toughness of the drying gel network (order 0.01 MPa·√m); controls the critical cracking thickness.' }),
    SN({ key: 'sint_T', label: 'Sintering hold temperature', unit: '°C', value: 600, min: 100, max: 1800, bus: 'proc.T_C', help: 'Hold temperature of the firing schedule.' }),
    SN({ key: 'sint_ramp', label: 'Heating rate', unit: 'K/min', value: 5, min: 0.1, max: 200, log: true, help: 'Linear ramp from room temperature.' }),
    SN({ key: 'sint_hold', label: 'Hold time', unit: 'min', value: 60, min: 0, max: 3000, help: 'Dwell at the hold temperature.' }),
    SN({ key: 'sint_Tref', label: 'Reference temperature for rate constants', unit: '°C', value: 600, min: 100, max: 1800, help: 'Temperature at which the densification and grain-growth constants below were measured.' }),
    SN({ key: 'sint_k', label: 'Densification rate constant at T_ref', unit: '1/s', value: 3e-3, min: 1e-7, max: 1, log: true, help: 'First-order pore-closure rate at the reference temperature and initial grain size.' }),
    SN({ key: 'sint_Q', label: 'Densification activation energy', unit: 'kJ/mol', value: 250, min: 50, max: 800, help: 'Apparent activation energy of the densifying mechanism.' }),
    SN({ key: 'sint_G0', label: 'Initial grain / particle size', unit: 'nm', value: 20, min: 2, max: 5000, log: true, help: 'Starting crystallite size; finer powders sinter faster (G⁻³).' }),
    SN({ key: 'sint_Kg', label: 'Grain-growth constant at T_ref', unit: 'nm³/s', value: 5, min: 1e-3, max: 1e9, log: true, help: 'Cubic grain-growth rate constant at the reference temperature.' }),
    SN({ key: 'sint_Qg', label: 'Grain-growth activation energy', unit: 'kJ/mol', value: 300, min: 50, max: 900, help: 'Activation energy of grain-boundary migration.' }),
    ES({ key: 'es_gamma', label: 'Solution surface tension', unit: 'mN/m', value: 40, min: 15, max: 75, help: 'Sets the Taylor-cone onset voltage and the terminal jet size.' }),
    ES({ key: 'es_H', label: 'Tip-to-collector distance', unit: 'cm', value: 15, min: 3, max: 50, help: 'Working distance; also the drying path.' }),
    ES({ key: 'es_L', label: 'Needle length', unit: 'cm', value: 4, min: 0.5, max: 15, help: 'Length of the capillary in the Taylor formula.' }),
    ES({ key: 'es_R', label: 'Needle inner radius', unit: 'mm', value: 0.4, min: 0.05, max: 2, log: true, help: 'Capillary radius.' }),
    ES({ key: 'es_V', label: 'Applied voltage', unit: 'kV', value: 15, min: 1, max: 60, help: 'Voltage between needle and collector.' }),
    ES({ key: 'es_Q', label: 'Flow rate', unit: 'mL/h', value: 1, min: 0.01, max: 50, log: true, help: 'Syringe-pump flow per needle; fibre diameter scales as Q^⅔ at fixed jet current.' }),
    ES({ key: 'es_I', label: 'Jet current', unit: 'nA', value: 100, min: 1, max: 1e4, log: true, help: 'Measured current carried by the jet; more conductive solutions carry more current and give thinner fibres.' }),
    ES({ key: 'es_c', label: 'Polymer volume fraction', unit: '%', value: 10, min: 0.5, max: 40, help: 'Solids content of the spinning solution.' }),
    ES({ key: 'es_rho', label: 'Polymer density', unit: 'kg/m³', value: 1200, min: 800, max: 2500, bus: 'prop.density', help: 'Density of the dry polymer (or polymer/precursor composite).' }),
    ES({ key: 'es_needles', label: 'Number of needles', value: 1, min: 1, max: 1000, step: 1, log: true, help: 'Parallel jets for scale-up.' }),
    CO({ key: 'col_V', label: 'Batch volume', unit: 'L', value: 10, min: 0.01, max: 20000, log: true, help: 'Liquid volume of the precipitation batch at the zinc concentration given under Chemistry.' }),
    CO({ key: 'col_conv', label: 'Conversion', unit: '–', value: 0.9, min: 0.01, max: 1, bus: 'rxn.conversion', help: 'Fraction of zinc converted to ZnO (take it from the reactor model or measurement).' }),
    CO({ key: 'col_wash', label: 'Number of washes', value: 3, min: 0, max: 10, step: 1, help: 'Centrifuge/redisperse or filter-cake washing cycles.' }),
    CO({ key: 'col_wratio', label: 'Wash ratio', unit: 'L/L', value: 5, min: 0.5, max: 50, log: true, help: 'Fresh water added per volume of liquor retained in the cake, per wash.' }),
    CO({ key: 'col_ret', label: 'Retained liquor', unit: 'L/kg', value: 2, min: 0.2, max: 20, log: true, help: 'Liquid held in the wet cake per kg of solids.' }),
    CO({ key: 'col_loss', label: 'Solids loss per wash', unit: '%', value: 2, min: 0, max: 20, help: 'Fines lost in each decant or filtrate.' }),
  );
  // the colloid batch also needs the zinc concentration input
  const z = INPUTS.find(i => i.key === 'cZn'); z.group = 'Chemistry';
}

const RUNNERS = { zno_reactor: runZnO, cvd: runCVD, ald: runALD, pvd: runPVD, litho_etch: runLitho, additive: runAM, solgel_spin_sinter: runSolGel, electrospin_colloid: runSpinColloid };
export default {
  id: 21,
  tagline: 'Simulates synthesis and fabrication processes — a coupled 3-D hydrothermal ZnO nanorod reactor with field heating, plus CVD, ALD, PVD, lithography/etch, additive, sol–gel/sintering and electrospinning models — and exports operating recipes.',
  inputs: INPUTS,
  presets: [
    { name: 'ZnO nanorods — porous reactor reference case (1000 V/m, 1000 s)', note: 'Case-study inputs: Zn²⁺ 0.01 mol/L, OH⁻ 0.02 mol/L, 298 K, 101 325 Pa, porosity 0.25, 1 m × 1 m × 0.01 m, h = 50 W/m²K. Runs on a 20×20×10 grid with Δt = 40 s; the nominal 100×100×10 mesh and 0.1 s step are beyond an in-browser run, so use the sensitivity indicators.', values: { mode: 'zno_reactor' } },
    { name: 'ZnO reference case at the nominal Δt = 0.1 s (first 50 s, coarse plan grid)', note: 'Same physics with the case-study time step of 0.1 s on a 4×4×10 grid for the first 50 s, to compare early-time behaviour against the larger default step.', values: { mode: 'zno_reactor', dt: 0.1, t_total: 50, nx: 4, ny: 4, nz: 10 } },
    { name: 'ZnO nanorods — laboratory recipe, Zn(NO₃)₂ + HMTA, 90 °C, 3 h', note: 'Classic equimolar 25 mmol/L zinc nitrate / HMTA bath in a 90 °C oven without field, closed batch in an open (90 % porous) seeded scaffold. The rate pre-factor is an illustrative calibration.', values: { mode: 'zno_reactor', base: 'HMTA', cZn: 25, c_hmta: 25, Tinf: 363.15, E: 0, field_nu: 0, t_total: 10800, dt: 180, nx: 8, ny: 8, phi: 0.9, biot: 0.95, dgrain: 2, bc: 'batch', Lx: 0.05, Ly: 0.05, A_pre: 3000, Emod: 0.2 } },
    { name: 'ZnO nanorods — field-assisted growth with PID temperature hold and dosing', note: 'Feedback-optimised synthesis: the field is modulated to hold 328 K and the bottom concentrations are dosed to hold 6 nm/min axial growth.', values: { mode: 'zno_reactor', ctrl: 'both', T_set: 328, t_total: 2000, dt: 25, nx: 10, ny: 10, G_set: 6, E: 1200, E_max: 3000, Kp: 100, Ki: 0.5 } },
    { name: 'LPCVD polysilicon from silane, 620 °C', note: 'Reaction-limited low-pressure CVD: uniform, conformal nanofilm growth of about 10 nm/min.', values: { mode: 'cvd' } },
    { name: 'ALD Al₂O₃ in aspect-ratio-50 nanopores', note: 'Conformal alumina coating of a nanoporous template (TMA/H₂O, 200 °C).', values: { mode: 'ald' } },
    { name: 'Metasurface pillars — immersion lithography and RIE', note: '70 nm features at 140 nm pitch by 193 nm immersion lithography followed by a 100 nm anisotropic etch.', values: { mode: 'litho_etch' } },
    { name: 'Auxetic / lattice metamaterial by LPBF Ti-6Al-4V', note: '20 % dense, 2 mm cell lattice: melt-pool window, minimum strut and build time.', values: { mode: 'additive', am_proc: 'lpbf' } },
    { name: 'Micro-lattice mechanical metamaterial by two-photon polymerisation', note: '100 µm lattice block with 10 µm cells written at NA 1.4.', values: { mode: 'additive', am_proc: 'tpp', am_cell: 0.01, am_rd: 0.1 } },
    { name: 'Sol–gel SiO₂ nanofilm: spin, dry, sinter', note: 'TEOS-derived silica film on a wafer with cracking check and densification.', values: { mode: 'solgel_spin_sinter' } },
  ],
  run(inp) {
    const fn = RUNNERS[inp.mode] || runZnO;
    return fn(inp);
  },
};
