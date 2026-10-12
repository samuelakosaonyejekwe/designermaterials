// Suite 5 — Metamaterial architecture and unit-cell design.
// Mechanical lattice cells (periodic beam-frame homogenisation with analytical cross-checks), TPMS architectures,
// resonant electromagnetic and acoustic cells, one-dimensional Bloch band design, bistable / origami / kirigami
// cells, negative-thermal-expansion and negative-compressibility cells, graded / hierarchical / finite arrays.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import * as implicit from '../geometry/implicit.js';
import * as layout from '../geometry/layout.js';

const F = N.fmt, PI = Math.PI, { sqrt, abs, min, max, log, exp, sin, cos, tan, round, floor, hypot, atan2 } = Math;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const fin = (v, d = 0) => (Number.isFinite(v) ? v : d);
const rad = a => a * PI / 180, deg = a => a * 180 / PI;
const unitOf = L => (L >= 1e-4 ? ['mm', 1e3] : L >= 1e-7 ? ['µm', 1e6] : ['nm', 1e9]);
const len = L => { const [u, s] = unitOf(L); return `${F(L * s, 4)} ${u}`; };
const freq = f => (f >= 1e12 ? `${F(f / 1e12, 4)} THz` : f >= 1e9 ? `${F(f / 1e9, 4)} GHz` : f >= 1e6 ? `${F(f / 1e6, 4)} MHz` : f >= 1e3 ? `${F(f / 1e3, 4)} kHz` : `${F(f, 4)} Hz`);

// Typical minimum feature sizes and overhang limits of fabrication routes (order-of-magnitude process capability).
const PROC = {
  tpp: { name: 'two-photon polymerisation', res: 2e-7, overhang: 0, planar: false },
  sla: { name: 'SLA / DLP vat photopolymerisation', res: 5e-5, overhang: 30, planar: false },
  lpbf: { name: 'laser powder-bed fusion', res: 1.5e-4, overhang: 45, planar: false },
  fff: { name: 'fused filament fabrication', res: 4e-4, overhang: 45, planar: false },
  ebl: { name: 'electron-beam lithography', res: 1e-8, overhang: 0, planar: true },
  photolith: { name: 'optical (contact / i-line) photolithography', res: 1e-6, overhang: 0, planar: true },
  pcb: { name: 'printed-circuit-board etching', res: 1e-4, overhang: 0, planar: true },
};

// ---------------------------------------------------------------- unit-cell registry
const OWN = {
  square2D: { name: 'Square grid', dim: 2, make: () => ({ nodes: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], edges: [[0, 1], [0, 2], [2, 3], [1, 3]], size: [1, 1, 0] }) },
  rhombic2D: { name: 'Wine-rack (rhombic) cell', dim: 2, make: p => { const t = tan(rad(p.theta_deg ?? 30)); return { nodes: [[0, t / 2, 0], [0.5, 0, 0], [1, t / 2, 0], [0.5, t, 0]], edges: [[0, 1], [1, 2], [2, 3], [3, 0]], size: [1, t, 0] }; } },
};
const cellDef = type => OWN[type] || implicit.LATTICES[type];
const CELL_OPTS = [
  ['reentrant2D', '2-D re-entrant honeycomb (auxetic)'], ['honeycomb2D', '2-D hexagonal honeycomb'], ['triangular2D', '2-D triangular'], ['kagome2D', '2-D kagome'], ['square2D', '2-D square grid'],
  ['chiral2D', '2-D tetrachiral (rings and ligaments)'], ['rotatingSquares2D', '2-D rotating squares (auxetic)'], ['starHoneycomb2D', '2-D four-point star (auxetic)'],
  ['octet', '3-D octet truss'], ['bcc', '3-D body-centred cubic'], ['fcc', '3-D face-centred cubic (corner–face struts)'], ['simpleCubic', '3-D simple cubic'], ['kelvin', '3-D Kelvin (truncated octahedron)'], ['diamond', '3-D diamond'], ['reentrant3D', '3-D re-entrant (auxetic)'],
].map(([value, label]) => ({ value, label }));
const ANGLE_CELLS = ['honeycomb2D', 'reentrant2D', 'reentrant3D', 'rotatingSquares2D'];
const cellParams = (type, theta, hl) => (type === 'rotatingSquares2D' ? { theta_deg: min(44, theta) } : { theta_deg: theta, h_over_l: hl });

// Periodic reduction: merge nodes that are images of one another and remove duplicated boundary struts.
// strip > 0 stacks that many cells along y and leaves the two y-faces free (finite-width array).
function reduceCell(uc, dim, a, strip = 0) {
  const size = uc.size || [1, 1, dim === 3 ? 1 : 0], per = [0, 1, 2].map(i => (size[i] || 0) * a), reps = [], map = new Map(), els = [], seen = new Set(), ny = strip || 1;
  const rep = q => {
    const f = [0, 0, 0], w = [0, 0, 0];
    for (let i = 0; i < dim; i++) { let u = q[i] / size[i]; if (!(strip && i === 1)) { u -= floor(u + 1e-6); if (u < 0) u = 0; } f[i] = u * per[i]; w[i] = round(u * 1e5); }
    const k = w.join(); let id = map.get(k); if (id === undefined) { id = reps.length; reps.push(f); map.set(k, id); } return id;
  };
  for (let j = 0; j < ny; j++) for (const e of uc.edges) {
    if (strip && e[3]) continue;
    const P = uc.nodes[e[0]], Q = uc.nodes[e[1]], A = [P[0], P[1] + j * size[1], P[2] || 0], B = [Q[0] + (e[2] || 0) * size[0], Q[1] + ((e[3] || 0) + j) * size[1], (Q[2] || 0) + (e[4] || 0) * (size[2] || 0)];
    let d = [0, 1, 2].map(i => (B[i] - A[i]) * a); const L = hypot(d[0], d[1], d[2]), eps = 1e-9 * a; if (L < eps) continue;
    let ia = rep(A), ib = rep(B);
    if (ia > ib || (ia === ib && (d[0] < -eps || (abs(d[0]) <= eps && (d[1] < -eps || (abs(d[1]) <= eps && d[2] < 0)))))) { [ia, ib] = [ib, ia]; d = d.map(x => -x); }
    const k = ia + ';' + ib + ';' + d.map(x => round(x / a * 1e4)).join(); if (seen.has(k)) continue; seen.add(k); els.push({ a: ia, b: ib, d, L });
  }
  return { reps, els, per, dim, vol: dim === 3 ? per[0] * ny * per[1] * per[2] : per[0] * ny * per[1] };
}
function tessellate(uc, dim, n, a, t, scale) {   // lattice geometry of an own cell, n cells per direction
  const size = uc.size || [1, 1, dim === 3 ? 1 : 0], nodes = [], edges = [], idx = new Map(), seen = new Set(), nz = dim === 3 ? n : 1;
  const add = p => { const k = p.map(x => round(x * 1e5)).join(); let i = idx.get(k); if (i === undefined) { i = nodes.length; nodes.push(p.map(x => x * a * scale)); idx.set(k, i); } return i; };
  for (let k = 0; k < nz; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) for (const e of uc.edges) {
    const o = [i * size[0], j * size[1], k * (size[2] || 0)], p = add(uc.nodes[e[0]].map((x, m) => (x || 0) + o[m])), q = add(uc.nodes[e[1]].map((x, m) => (x || 0) + o[m])), key = min(p, q) * 1e6 + max(p, q);
    if (p !== q && !seen.has(key)) { seen.add(key); edges.push([p, q]); }
  }
  return { kind: 'lattice', nodes, edges, radius: t * scale / 2, units: '' };
}

// ---------------------------------------------------------------- periodic beam-frame homogenisation
const VP = [[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]];   // Voigt order xx, yy, zz, yz, xz, xy (engineering shear)
function section(dim, t, m) {
  const G = m.Es / (2 * (1 + m.nus));
  if (dim === 3) { const r = t / 2, A = PI * r * r, I = PI * r ** 4 / 4; return { A, I, c: r, EA: m.Es * A, EI: m.Es * I, GJ: G * 2 * I, GAs: 0.9 * G * A }; }
  const A = t, I = t ** 3 / 12; return { A, I, c: t / 2, EA: m.Es * A, EI: m.Es * I, GJ: G * 4 * I, GAs: 5 / 6 * G * A };   // per unit depth
}
function localK(L, EA, EI, GJ, ph) {
  const k = Array.from({ length: 12 }, () => new Float64Array(12)), s = (i, j, v) => { k[i][j] = v; k[j][i] = v; };
  const a = EA / L, t = GJ / L, b = EI / (L ** 3 * (1 + ph)), b1 = 12 * b, b2 = 6 * L * b, b3 = (4 + ph) * L * L * b, b4 = (2 - ph) * L * L * b;
  s(0, 0, a); s(6, 6, a); s(0, 6, -a); s(3, 3, t); s(9, 9, t); s(3, 9, -t);
  s(1, 1, b1); s(1, 5, b2); s(1, 7, -b1); s(1, 11, b2); s(5, 5, b3); s(5, 7, -b2); s(5, 11, b4); s(7, 7, b1); s(7, 11, -b2); s(11, 11, b3);
  s(2, 2, b1); s(2, 4, -b2); s(2, 8, -b1); s(2, 10, -b2); s(4, 4, b3); s(4, 8, b2); s(4, 10, b4); s(8, 8, b1); s(8, 10, b2); s(10, 10, b3);
  return k;
}
function solveMulti(K, B) {   // symmetric positive (semi-)definite K, several right-hand sides; Gaussian elimination with pivot guard
  const n = K.length, m = n ? B[0].length : 0; let scale = 0, ok = true; for (let i = 0; i < n; i++) scale = max(scale, abs(K[i][i]));
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (abs(K[i][k]) > abs(K[p][k])) p = i;
    if (p !== k) { [K[p], K[k]] = [K[k], K[p]]; [B[p], B[k]] = [B[k], B[p]]; }
    if (abs(K[k][k]) < 1e-13 * scale) { K[k][k] = 1e-9 * scale || 1; ok = false; }
    const row = K[k], piv = row[k];
    for (let i = k + 1; i < n; i++) { const f = K[i][k] / piv; if (f === 0) continue; const ri = K[i]; for (let j = k; j < n; j++) ri[j] -= f * row[j]; for (let j = 0; j < m; j++) B[i][j] -= f * B[k][j]; }
  }
  for (let i = n - 1; i >= 0; i--) { const ri = K[i]; for (let j = 0; j < m; j++) { let s = B[i][j]; for (let q = i + 1; q < n; q++) s -= ri[q] * B[q][j]; B[i][j] = s / ri[i]; } }
  return ok;
}
// Timoshenko beam elements on the periodic graph; macroscopic strain enters through the lattice translation of each strut.
function frameSolve(red, sec, comps, bend = 1) {
  const { reps, els, dim } = red, dm = dim === 3 ? [0, 1, 2, 3, 4, 5] : [0, 1, 5], nd = dm.length, nc = comps.length, n = reps.length * nd - dim;
  const fid = (node, k) => { const g = node * nd + k; return g < dim ? -1 : g - dim; };
  const K = Array.from({ length: n }, () => new Float64Array(n)), Kue = Array.from({ length: n }, () => new Float64Array(nc)), Kee = Array.from({ length: nc }, () => new Float64Array(nc)), data = [];
  for (const e of els) {
    const L = e.L, ex = e.d.map(x => x / L); let ey, ez;
    if (dim === 2) { ey = [-ex[1], ex[0], 0]; ez = [0, 0, 1]; } else { const r = abs(ex[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0]; ey = [r[1] * ex[2] - r[2] * ex[1], r[2] * ex[0] - r[0] * ex[2], r[0] * ex[1] - r[1] * ex[0]]; const q = hypot(ey[0], ey[1], ey[2]); ey = ey.map(x => x / q); ez = [ex[1] * ey[2] - ex[2] * ey[1], ex[2] * ey[0] - ex[0] * ey[2], ex[0] * ey[1] - ex[1] * ey[0]]; }
    const R = [ex, ey, ez], kl = localK(L, sec.EA, sec.EI * bend, sec.GJ * bend, bend === 1 ? 12 * sec.EI / (sec.GAs * L * L) : 0);
    const Kg = Array.from({ length: 12 }, () => new Float64Array(12));
    for (let bi = 0; bi < 4; bi++) for (let bj = 0; bj < 4; bj++) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { let s = 0; for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) s += R[p][i] * kl[3 * bi + p][3 * bj + q] * R[q][j]; Kg[3 * bi + i][3 * bj + j] = s; }
    const D = [0, 1, 2].map(() => new Float64Array(nc));
    comps.forEach((I, c) => { const [p, q] = VP[I]; if (p === q) D[p][c] = e.d[p]; else { D[p][c] = e.d[q] / 2; D[q][c] = e.d[p] / 2; } });
    const KB = Array.from({ length: 12 }, (_, g) => { const r = new Float64Array(nc); for (let c = 0; c < nc; c++) r[c] = Kg[g][6] * D[0][c] + Kg[g][7] * D[1][c] + Kg[g][8] * D[2][c]; return r; });
    const gi = [], fi = []; for (let s = 0; s < 2; s++) for (let k = 0; k < nd; k++) { gi.push(s * 6 + dm[k]); fi.push(fid(s ? e.b : e.a, k)); }
    for (let p = 0; p < gi.length; p++) { if (fi[p] < 0) continue; const row = K[fi[p]], kr = Kg[gi[p]]; for (let q = 0; q < gi.length; q++) if (fi[q] >= 0) row[fi[q]] += kr[gi[q]]; for (let c = 0; c < nc; c++) Kue[fi[p]][c] += KB[gi[p]][c]; }
    for (let c = 0; c < nc; c++) for (let c2 = 0; c2 < nc; c2++) Kee[c][c2] += D[0][c] * KB[6][c2] + D[1][c] * KB[7][c2] + D[2][c] * KB[8][c2];
    data.push({ R, kl, D, gi, fi, L });
  }
  const U = Kue.map(r => Float64Array.from(r, x => -x)), ok = solveMulti(K, U), C = Array.from({ length: nc }, () => new Array(nc).fill(0));
  for (let c = 0; c < nc; c++) for (let c2 = 0; c2 < nc; c2++) { let s = Kee[c][c2]; for (let i = 0; i < n; i++) s += Kue[i][c] * U[i][c2]; C[c][c2] = s / red.vol; }
  return { C, U, data, ok, dof: n };
}
function memberForces(sol, e) {   // axial force and largest end moment of every strut for the macroscopic strain vector e
  return sol.data.map(d => {
    const q = new Float64Array(12);
    d.gi.forEach((g, p) => { if (d.fi[p] >= 0) { let s = 0; for (let c = 0; c < e.length; c++) s += sol.U[d.fi[p]][c] * e[c]; q[g] = s; } });
    for (let i = 0; i < 3; i++) for (let c = 0; c < e.length; c++) q[6 + i] += d.D[i][c] * e[c];
    const ql = new Float64Array(12); for (let b = 0; b < 4; b++) for (let i = 0; i < 3; i++) ql[3 * b + i] = d.R[i][0] * q[3 * b] + d.R[i][1] * q[3 * b + 1] + d.R[i][2] * q[3 * b + 2];
    const f = i => { let s = 0; for (let j = 0; j < 12; j++) s += d.kl[i][j] * ql[j]; return s; };
    return { N: f(6), M: max(hypot(f(4), f(5)), hypot(f(10), f(11))), L: d.L };
  });
}
function invSym(C) { const n = C.length, sc = max(...C.map((r, i) => abs(r[i]))) || 1, A = C.map((r, i) => r.map((x, j) => x / sc + (i === j ? 1e-12 : 0))); try { return N.inverse(A).map(r => r.map(x => x / sc)); } catch { return C.map((r, i) => r.map((x, j) => (i === j ? 1 / (x || sc) : 0))); } }
function unionFraction(red, t) {   // solid fraction with node overlap counted once: quasi-random sampling of the union of struts
  const { reps, els, per, dim } = red, r = t / 2, r2 = r * r, segs = [], sh = dim === 3 ? [-1, 0, 1] : [0];
  for (const e of els) for (const i of [-1, 0, 1]) for (const j of [-1, 0, 1]) for (const k of sh) {
    const p = [reps[e.a][0] + i * per[0], reps[e.a][1] + j * per[1], reps[e.a][2] + k * per[2]]; let keep = true;
    for (let m = 0; m < dim; m++) if (max(p[m], p[m] + e.d[m]) < -r || min(p[m], p[m] + e.d[m]) > per[m] + r) keep = false;
    if (keep) segs.push([p[0], p[1], p[2], e.d[0], e.d[1], e.d[2], e.L * e.L]);
  }
  // Kronecker (additive-recurrence) sequence: a regular grid aliases with axis-aligned struts and biases the area of thin members
  const n = dim === 3 ? 90000 : 120000, g = dim === 3 ? 1.2207440846057596 : 1.324717957244746, a1 = 1 / g, a2 = a1 / g, a3 = a2 / g; let inside = 0;
  // bin the segments (bounding box widened by r) on a coarse grid of the cell so that each sample point tests only nearby struts
  const G = dim === 3 ? 8 : 16, Gz = dim === 3 ? G : 1, bins = Array.from({ length: G * G * Gz }, () => []), bi = (v, m, g) => min(g - 1, max(0, floor(v / per[m] * g)));
  for (const s of segs) { const lo = [0, 1, 2].map(m => min(s[m], s[m] + s[3 + m]) - r), hi = [0, 1, 2].map(m => max(s[m], s[m] + s[3 + m]) + r);
    for (let kz = dim === 3 ? bi(lo[2], 2, Gz) : 0; kz <= (dim === 3 ? bi(hi[2], 2, Gz) : 0); kz++) for (let jy = bi(lo[1], 1, G); jy <= bi(hi[1], 1, G); jy++) for (let ix = bi(lo[0], 0, G); ix <= bi(hi[0], 0, G); ix++) bins[ix + G * (jy + G * kz)].push(s); }
  for (let q = 1; q <= n; q++) {
    const x = (0.5 + q * a1) % 1 * per[0], y = (0.5 + q * a2) % 1 * per[1], z = dim === 3 ? (0.5 + q * a3) % 1 * per[2] : 0;
    for (const s of bins[bi(x, 0, G) + G * (bi(y, 1, G) + G * (dim === 3 ? bi(z, 2, Gz) : 0))]) { const ax = x - s[0], ay = y - s[1], az = z - s[2], u = min(1, max(0, (ax * s[3] + ay * s[4] + az * s[5]) / s[6])), bx = ax - u * s[3], by = ay - u * s[4], bz = az - u * s[5]; if (bx * bx + by * by + bz * bz <= r2) { inside++; break; } }
  }
  return inside / n;
}
// Effective properties of one cell. Primary load axis: x for 2-D cells, z for 3-D cells.
function cellModel(type, prm, a, tIn, m, light = false) {
  const def = cellDef(type), dim = def.dim, uc = def.make(prm || {}), red = reduceCell(uc, dim, a), Lmin = N.minOf(red.els.map(e => e.L)), t = min(tIn, 0.45 * Lmin), sec = section(dim, t, m);
  const comps = dim === 3 ? [0, 1, 2, 3, 4, 5] : [0, 1, 5], sol = frameSolve(red, sec, comps), S = invSym(sol.C), ia = dim === 3 ? 2 : 0;
  const rho1 = min(1, N.sum(red.els.map(e => e.L)) * sec.A / red.vol);
  const p = dim === 3
    ? { E: 1 / S[2][2], Eo: 1 / S[0][0], nu: -S[0][2] / S[2][2], nuo: -S[1][0] / S[0][0], G: 1 / S[4][4], Go: 1 / S[5][5], K: 1 / (S[0][0] + S[1][1] + S[2][2] + 2 * (S[0][1] + S[0][2] + S[1][2])) }
    : { E: 1 / S[0][0], Eo: 1 / S[1][1], nu: -S[1][0] / S[0][0], nuo: -S[0][1] / S[1][1], G: 1 / S[2][2], Go: 1 / S[2][2], K: 1 / (S[0][0] + S[1][1] + 2 * S[0][1]) };
  const out = { type, dim, uc, red, sec, sol, S, C: sol.C, t, capped: t < tIn, rho1, Lmin, ia, comps, ...p, Z: 2 * red.els.length / red.reps.length, j: red.reps.length, b: red.els.length };
  if (light) return out;
  // member forces for unit uniaxial stress along the primary axis → first yield and Euler buckling of the struts
  const mf = memberForces(sol, S.map(r => r[ia])); let smax = 0, bmin = Infinity, Nmax = 0;
  for (const f of mf) { smax = max(smax, abs(f.N) / sec.A + f.M * sec.c / sec.I); Nmax = max(Nmax, abs(f.N)); if (f.N > 0) bmin = min(bmin, PI * PI * sec.EI / (f.L * f.L) / f.N); }
  out.yieldRatio = smax > 0 ? 1 / smax : 0;                      // σ*_y / σ_ys
  out.buckle = Number.isFinite(bmin) ? bmin : null;               // macroscopic compressive stress at first strut buckling (Pa)
  // share of the stiffness carried by axial action: repeat with the bending and torsion stiffness removed
  const tr = frameSolve(red, sec, comps, 1e-7), St = invSym(tr.C), Et = 1 / St[ia][ia];
  out.stretch = min(1, max(0, fin(Et / out.E)));
  return out;
}
const dirE = (S, dim, nvec) => { const m = dim === 3 ? [nvec[0] ** 2, nvec[1] ** 2, nvec[2] ** 2, nvec[1] * nvec[2], nvec[0] * nvec[2], nvec[0] * nvec[1]] : [nvec[0] ** 2, nvec[1] ** 2, nvec[0] * nvec[1]]; let s = 0; for (let i = 0; i < m.length; i++) for (let j = 0; j < m.length; j++) s += m[i] * S[i][j] * m[j]; return 1 / s; };
// Closed-form references (ratios to E_s). Hexagonal cells: Gibson–Ashby with axial and shear terms.
function hexGeom(c) { const n = c.uc.nodes, s = c.red.per[0], l = hypot(n[2][0] - n[1][0], n[2][1] - n[1][1]) * s, h = hypot(n[1][0] - n[0][0], n[1][1] - n[0][1]) * s; return { l, h, th: atan2(n[2][1] - n[1][1], n[2][0] - n[1][0]) }; }
function reference(c, m) {
  const t = c.t, rho = c.rho1;
  if (c.type === 'honeycomb2D' || c.type === 'reentrant2D') {
    const { l, h, th } = hexGeom(c), r = t / l, hl = h / l, sn = sin(th), cs = cos(th), k = 2.4 + 1.5 * m.nus, safe = abs(sn) > 0.02;
    const E1f = safe ? r ** 3 * cs / ((hl + sn) * sn * sn) : null, nuF = safe ? cs * cs / ((hl + sn) * sn) : null, ct2 = (cs / sn) ** 2;
    const E1 = r ** 3 * cs / ((hl + sn) * (sn * sn + (k * sn * sn + cs * cs) * r * r)), nu = cs * cs * sn / (hl + sn) * (1 + (k - 1) * r * r) / (sn * sn + (k * sn * sn + cs * cs) * r * r);
    const E2 = r ** 3 * (hl + sn) / (cs ** 3 * (1 + (k + (sn / cs) ** 2 + 2 * hl / (cs * cs)) * r * r)), G = r ** 3 * (hl + sn) / (hl * hl * (1 + 2 * hl) * cs);
    const nb = N.interp1([1, 1.5, 2], [0.686, 0.76, 0.806], hl);
    return { label: 'Gibson–Ashby honeycomb model (flexure with axial and shear terms)', E: E1, Eo: E2, nu, nuo: nu * E2 / E1, G, E1f, nuF, l, h, th, hl, r, ct2,
      rho: r * (hl + 2) / (2 * cs * (hl + sn)), plastic: r * r / (2 * (hl + sn) * max(0.02, abs(sn))), buckle2: nb * nb * PI * PI / 24 * t ** 3 / (l * h * h * cs) };
  }
  if (c.type === 'triangular2D' || c.type === 'kagome2D') return { label: `${c.type === 'kagome2D' ? 'Kagome' : 'Triangular'} lattice, pin-jointed stretch model`, E: rho / 3, Eo: rho / 3, nu: 1 / 3, nuo: 1 / 3, G: rho / 8, plastic: rho / 3 };
  if (c.type === 'square2D') return { label: 'Square grid: axial along the walls, wall bending in shear', E: rho / 2, Eo: rho / 2, nu: 0, nuo: 0, G: (t / c.Lmin) ** 3 / 2, plastic: rho / 2 };
  if (c.type === 'octet') return { label: 'Octet truss, pin-jointed (Deshpande–Fleck–Ashby)', E: rho / 9, Eo: rho / 9, nu: 1 / 3, nuo: 1 / 3, G: rho / 12 };
  if (c.type === 'simpleCubic') return { label: 'Simple cubic: axial along the struts', E: rho / 3, Eo: rho / 3, nu: 0, nuo: 0, G: null };
  if (c.type === 'rotatingSquares2D') return { label: 'Rotating rigid squares (Grima–Evans kinematics)', E: null, nu: -1, nuo: -1, G: null };
  return null;
}

// ---------------------------------------------------------------- inputs
const SH = (m, extra) => ({ mode: Array.isArray(m) ? m : [m], ...(extra || {}) });
const RES = (...t) => SH('resonant', { resType: t }), BIS = (...t) => SH('bistable', { bisType: t }), NTE = (...t) => SH('negcte', { nteType: t });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'lattice', group: 'Model', options: [
    { value: 'lattice', label: 'Mechanical lattice unit cell (stiffness, Poisson ratio, strength, manufacturability)' }, { value: 'tpms', label: 'TPMS architecture (gyroid, Schwarz, IWP, Neovius)' },
    { value: 'resonant', label: 'Resonant unit cell (SRR, wire medium, Helmholtz, mass-in-mass, membrane)' }, { value: 'band', label: 'Periodic-structure band design (1-D Bloch, Bragg)' },
    { value: 'bistable', label: 'Bistable, origami and kirigami cells' }, { value: 'negcte', label: 'Negative thermal expansion / negative compressibility cells' },
    { value: 'graded', label: 'Graded, hierarchical and finite arrays' }], help: 'Selects the unit-cell design model. Each analysis builds the cell geometry and publishes its effective properties to the shared bus.' },
  { key: 'cellType', label: 'Lattice cell', type: 'select', value: 'reentrant2D', group: 'Lattice', options: CELL_OPTS, showIf: SH(['lattice', 'graded']), help: 'Strut topology of the unit cell. Two-dimensional cells are prismatic (in-plane wall thickness, unit depth); three-dimensional cells use circular struts.' },
  { key: 'cell', label: 'Cell size (period along x) a', unit: 'm', value: 5e-3, min: 1e-7, max: 0.5, log: true, group: 'Lattice', bus: 'geom.cell_size', showIf: SH(['lattice', 'tpms', 'resonant', 'band', 'bistable', 'negcte', 'graded']), help: 'Lattice period. For curved-beam cells it is the beam span; for resonant and band cells it is the lattice constant that must stay well below the operating wavelength.' },
  { key: 'strut', label: 'Strut diameter / wall thickness t', unit: 'm', value: 2.5e-4, min: 2e-8, max: 0.05, log: true, group: 'Lattice', bus: 'geom.strut', showIf: SH(['lattice', 'bistable', 'negcte', 'graded']), help: 'Diameter of circular struts (3-D cells) or in-plane wall thickness (2-D cells, beams). It is limited to 45 % of the shortest strut length. Stiffness scales with a power of t/a set by the deformation mode.' },
  { key: 'theta', label: 'Cell angle θ (magnitude)', unit: '°', value: 30, min: 2, max: 75, group: 'Lattice', showIf: SH(['lattice', 'graded', 'bistable']), help: 'Inclination of the oblique ribs to the x-axis for honeycomb cells (re-entrant cells use −θ), the rotation of each square for rotating squares (limited to 44°), or the re-entrant angle of the 3-D cell. Not used by fixed-topology cells.' },
  { key: 'hl', label: 'Rib length ratio h/l', unit: '–', value: 2, min: 0.3, max: 4, group: 'Lattice', showIf: SH(['lattice', 'graded']), help: 'Length of the ribs along y (h) over the oblique rib length (l) for honeycomb and re-entrant cells; it is raised automatically if opposite re-entrant vertices would touch.' },
  { key: 'Es', label: 'Constituent Young’s modulus E_s', unit: 'Pa', value: 1.14e11, min: 1e5, max: 1.2e12, log: true, group: 'Constituent', bus: 'prop.E', showIf: SH(['lattice', 'tpms', 'band', 'bistable', 'negcte', 'graded']), help: 'Young’s modulus of the solid from which the cell is made (default: Ti-6Al-4V). Effective moduli scale linearly with it.' },
  { key: 'nus', label: 'Constituent Poisson’s ratio', unit: '–', value: 0.34, min: 0, max: 0.49, group: 'Constituent', bus: 'prop.nu', showIf: SH(['lattice', 'tpms', 'graded']), help: 'Poisson’s ratio of the solid; it enters through the shear and torsional stiffness of the struts.' },
  { key: 'ys', label: 'Constituent yield strength', unit: 'Pa', value: 8.8e8, min: 1e4, max: 1e10, log: true, group: 'Constituent', bus: 'prop.yield', showIf: SH(['lattice', 'bistable', 'graded']), help: 'Yield (or fracture) strength of the solid, used for the first-yield strength of the lattice and the strain limit of flexures.' },
  { key: 'rhos', label: 'Constituent density', unit: 'kg/m³', value: 4430, min: 10, max: 23000, group: 'Constituent', bus: 'prop.density', showIf: SH(['lattice', 'tpms', 'band', 'graded']), help: 'Density of the solid; the effective density is the relative density times this value.' },
  { key: 'nView', label: 'Cells per side in the view', unit: '–', value: 3, min: 1, max: 4, step: 1, group: 'Lattice', sens: false, showIf: SH(['lattice', 'negcte']), help: 'Number of unit cells tessellated in each direction for the viewer and the exported node/edge files. It does not affect the homogenised properties.' },
  { key: 'process', label: 'Fabrication process', type: 'select', value: 'lpbf', group: 'Manufacturability', options: Object.keys(PROC).map(k => ({ value: k, label: `${PROC[k].name} (≈ ${len(PROC[k].res)} features)` })), showIf: SH(['lattice', 'tpms', 'resonant', 'graded']), help: 'Sets the typical minimum feature size and overhang limit against which struts, walls, traces and gaps are checked.' },
  { key: 'tolPct', label: 'Feature tolerance (±)', unit: '%', value: 8, min: 0, max: 40, group: 'Manufacturability', showIf: SH('lattice'), help: 'As-built deviation of the strut diameter or wall thickness, used for the finite-difference tolerance sensitivity of the effective properties.' },
  { key: 'tpmsType', label: 'TPMS family', type: 'select', value: 'gyroid', group: 'TPMS', showIf: SH('tpms'), options: [{ value: 'gyroid', label: 'Gyroid' }, { value: 'schwarzP', label: 'Schwarz P (primitive)' }, { value: 'schwarzD', label: 'Schwarz D (diamond)' }, { value: 'iwp', label: 'I-WP' }, { value: 'neovius', label: 'Neovius' }], help: 'Triply periodic minimal surface defined by its leading-order level-set function.' },
  { key: 'tpmsMode', label: 'Solid definition', type: 'select', value: 'sheet', group: 'TPMS', showIf: SH('tpms'), options: [{ value: 'sheet', label: 'Sheet (double-sided wall |f| ≤ c)' }, { value: 'solid', label: 'Solid network (one labyrinth f ≤ c)' }], help: 'Sheet structures thicken the surface into a wall; network structures fill one of the two channels.' },
  { key: 'rd', label: 'Relative density', unit: '–', value: 0.3, min: 0.05, max: 0.8, group: 'TPMS', bus: 'geom.rel_density', showIf: SH(['tpms', 'graded']), help: 'Solid volume fraction. For TPMS the level-set threshold is found so that this value is met; in the hierarchical analysis it is the total relative density to be split between two levels.' },
  { key: 'rdTo', label: 'Relative density at the far end (graded)', unit: '–', value: 0.3, min: 0.05, max: 0.8, group: 'TPMS', showIf: SH('tpms'), help: 'Set different from the relative density to grade the structure linearly along z.' },
  { key: 'plC', label: 'Stiffness power-law coefficient C', unit: '–', value: 0.5, min: 0.05, max: 1.5, group: 'TPMS', showIf: SH('tpms'), help: 'Coefficient of the empirical relation E*/E_s = C·ρ̄ⁿ. Reported fits for TPMS lattices lie roughly between 0.2 and 1; calibrate on tested specimens.' },
  { key: 'plN', label: 'Stiffness power-law exponent n', unit: '–', value: 1.4, min: 1, max: 3, group: 'TPMS', showIf: SH('tpms'), help: 'Exponent of the empirical relation: near 1 for stretch-dominated sheet structures, near 2 for bending-dominated networks.' },
  { key: 'resType', label: 'Resonant cell', type: 'select', value: 'srrwire', group: 'Resonator', showIf: SH('resonant'), options: [{ value: 'srrwire', label: 'SRR + wire (double-negative)' }, { value: 'srr', label: 'Split-ring resonator (negative μ)' }, { value: 'wire', label: 'Thin-wire medium (negative ε)' }, { value: 'helmholtz', label: 'Helmholtz resonator (negative modulus)' }, { value: 'massinmass', label: 'Mass-in-mass (negative effective mass)' }, { value: 'membrane', label: 'Membrane-type (negative density)' }], help: 'Electromagnetic cells use lumped L–C estimates; acoustic and elastic cells use lumped mass–spring–compliance models.' },
  { key: 'srrSize', label: 'SRR outer side', unit: 'm', value: 3.5e-3, min: 5e-8, max: 0.4, log: true, group: 'Resonator', showIf: RES('srr', 'srrwire'), help: 'Outer side length of the square ring; limited to 95 % of the cell size. Larger rings lower the resonance.' },
  { key: 'srrW', label: 'SRR trace width', unit: 'm', value: 4e-4, min: 5e-9, max: 0.05, log: true, group: 'Resonator', showIf: RES('srr', 'srrwire'), help: 'Width of the metal trace; limited to 30 % of the ring side.' },
  { key: 'srrGap', label: 'SRR gap', unit: 'm', value: 3e-4, min: 5e-9, max: 0.05, log: true, group: 'Resonator', showIf: RES('srr', 'srrwire'), help: 'Width of the split. A smaller gap raises the capacitance and lowers the resonance.' },
  { key: 'srrT', label: 'Metal thickness', unit: 'm', value: 3.5e-5, min: 5e-9, max: 5e-3, log: true, group: 'Resonator', showIf: RES('srr', 'srrwire'), help: 'Thickness of the metal layer (35 µm is standard PCB copper).' },
  { key: 'srrCadd', label: 'Additional gap capacitance', unit: 'F', value: 5e-14, min: 0, max: 1e-10, group: 'Resonator', showIf: RES('srr', 'srrwire'), help: 'Inter-ring capacitance of a double SRR or a lumped capacitor across the gap; zero for a plain single ring. It is the main handle for bringing the resonance below the wire plasma frequency.' },
  { key: 'sigma', label: 'Metal conductivity', unit: 'S/m', value: 5.8e7, min: 1e4, max: 7e7, log: true, group: 'Resonator', bus: 'prop.sigma_elec', showIf: RES('srr', 'srrwire', 'wire'), help: 'Conductivity of the ring and wire metal (copper 5.8e7 S/m); it sets the ohmic loss and the quality factor.' },
  { key: 'epsH', label: 'Host relative permittivity', unit: '–', value: 1, min: 1, max: 15, group: 'Resonator', bus: 'prop.eps_r', showIf: RES('srr', 'srrwire', 'wire'), help: 'Effective permittivity of the medium around the metal (about (1 + ε_r)/2 for traces on a substrate). It lowers both the resonance and the plasma frequency by its square root.' },
  { key: 'wireR', label: 'Wire radius', unit: 'm', value: 1e-4, min: 5e-9, max: 0.05, log: true, group: 'Resonator', showIf: RES('wire', 'srrwire'), help: 'Radius of the continuous wires (one per cell, parallel to the electric field); limited to 20 % of the cell size.' },
  { key: 'fluid', label: 'Fluid', type: 'select', value: 'air', group: 'Resonator', showIf: RES('helmholtz', 'membrane'), options: [{ value: 'air', label: 'Air (343 m/s, 1.204 kg/m³)' }, { value: 'water', label: 'Water (1481 m/s, 998 kg/m³)' }], help: 'Host fluid of the acoustic metamaterial at 20 °C.' },
  { key: 'hrV', label: 'Helmholtz cavity volume', unit: 'm³', value: 5e-8, min: 1e-18, max: 1e-2, log: true, group: 'Resonator', showIf: RES('helmholtz'), help: 'Volume of the resonator cavity; the resonance scales as V^(−1/2).' },
  { key: 'hrRn', label: 'Neck radius', unit: 'm', value: 5e-4, min: 1e-7, max: 0.05, log: true, group: 'Resonator', showIf: RES('helmholtz'), help: 'Radius of the circular neck. End corrections of 0.85 radius per flanged end are added to the neck length.' },
  { key: 'hrLn', label: 'Neck length', unit: 'm', value: 1e-3, min: 1e-7, max: 0.1, log: true, group: 'Resonator', showIf: RES('helmholtz'), help: 'Physical length of the neck (wall thickness for a perforated plate).' },
  { key: 'mimM1', label: 'Outer (shell) mass m₁', unit: 'kg', value: 1e-3, min: 1e-15, max: 10, log: true, group: 'Resonator', showIf: RES('massinmass'), help: 'Mass of the outer shell of each cell, connected to its neighbours.' },
  { key: 'mimM2', label: 'Inner (resonator) mass m₂', unit: 'kg', value: 2e-3, min: 1e-15, max: 10, log: true, group: 'Resonator', showIf: RES('massinmass'), help: 'Mass of the internal resonator. A larger ratio m₂/m₁ widens the negative-mass band.' },
  { key: 'mimK1', label: 'Lattice spring k₁', unit: 'N/m', value: 1e5, min: 1e-3, max: 1e10, log: true, group: 'Resonator', showIf: RES('massinmass'), help: 'Stiffness of the springs connecting neighbouring shells.' },
  { key: 'mimK2', label: 'Internal spring k₂', unit: 'N/m', value: 5e4, min: 1e-3, max: 1e10, log: true, group: 'Resonator', showIf: RES('massinmass'), help: 'Stiffness of the spring between shell and internal mass; with m₂ it sets the local resonance.' },
  { key: 'memR', label: 'Membrane radius', unit: 'm', value: 2e-3, min: 1e-6, max: 0.2, log: true, group: 'Resonator', showIf: RES('membrane'), help: 'Radius of the clamped circular membrane (equal to the tube radius).' },
  { key: 'memT', label: 'Membrane tension', unit: 'N/m', value: 50, min: 1e-3, max: 1e4, log: true, group: 'Resonator', showIf: RES('membrane'), help: 'In-plane tension per unit length; the cut-off frequency scales with its square root.' },
  { key: 'memRhoS', label: 'Membrane areal density', unit: 'kg/m²', value: 0.1, min: 1e-4, max: 10, log: true, group: 'Resonator', showIf: RES('membrane'), help: 'Mass per unit area of the membrane.' },
  { key: 'waveType', label: 'Wave type', type: 'select', value: 'elastic', group: 'Band design', showIf: SH('band'), options: [{ value: 'elastic', label: 'Elastic / acoustic (longitudinal, normal incidence)' }, { value: 'optical', label: 'Optical (Bragg stack, normal incidence)' }], help: 'Elastic stacks use modulus and density of the two phases; optical stacks use the two refractive indices.' },
  { key: 'fill', label: 'Fill fraction of phase A', unit: '–', value: 0.5, min: 0.05, max: 0.95, group: 'Band design', showIf: SH('band'), help: 'Thickness of layer A divided by the period. The gap is widest near the quarter-wave condition (equal transit time or optical thickness in both layers).' },
  { key: 'E2', label: 'Second-constituent modulus', unit: 'Pa', value: 3e9, min: 1e4, max: 1.2e12, log: true, group: 'Constituent', showIf: SH(['band', 'negcte']), help: 'Modulus of phase B of the periodic stack (default: epoxy), or of the base member of the bi-material triangle.' },
  { key: 'rho2', label: 'Second-constituent density', unit: 'kg/m³', value: 1200, min: 1, max: 23000, group: 'Constituent', showIf: SH('band'), help: 'Density of phase B of the elastic stack.' },
  { key: 'nA', label: 'Refractive index of phase A', unit: '–', value: 2.4, min: 1, max: 5, group: 'Band design', bus: 'prop.n_refr', showIf: SH('band', { waveType: ['optical'] }), help: 'Refractive index of the high-index layer (2.4 is typical of TiO₂ in the visible).' },
  { key: 'nB', label: 'Refractive index of phase B', unit: '–', value: 1.46, min: 1, max: 5, group: 'Band design', showIf: SH('band', { waveType: ['optical'] }), help: 'Refractive index of the low-index layer (1.46 for SiO₂).' },
  { key: 'fTarget', label: 'Target operating frequency', unit: 'Hz', value: 2e5, min: 1, max: 1e16, log: true, group: 'Band design', bus: 'res.f0', showIf: SH('band'), help: 'Frequency at which the centre of the first band gap is wanted; the required period follows from the scale invariance of the dispersion relation.' },
  { key: 'bisType', label: 'Cell', type: 'select', value: 'beam', group: 'Mechanism', showIf: SH('bistable'), options: [{ value: 'beam', label: 'Curved (cosine) beam — bistability and negative stiffness' }, { value: 'miura', label: 'Miura-ori fold pattern' }, { value: 'kirigami', label: 'Kirigami (staggered parallel cuts)' }, { value: 'rotsq', label: 'Rotating squares kinematics' }], help: 'Mechanism-type cells whose response is set by geometry: snap-through beams, rigid-foldable origami, cut patterns and rotating units.' },
  { key: 'bh', label: 'Beam apex height h', unit: 'm', value: 7.5e-4, min: 1e-8, max: 0.1, log: true, group: 'Mechanism', showIf: BIS('beam'), help: 'Initial rise of the cosine-shaped beam. Bistability requires h/t above about 2.31.' },
  { key: 'depth', label: 'Out-of-plane depth b', unit: 'm', value: 5e-3, min: 1e-8, max: 0.5, log: true, group: 'Mechanism', showIf: BIS('beam', 'kirigami'), help: 'Depth of the beam (or sheet thickness direction for forces); forces scale linearly with it.' },
  { key: 'miA', label: 'Miura parallelogram side a', unit: 'm', value: 1e-2, min: 1e-6, max: 1, log: true, group: 'Mechanism', showIf: BIS('miura'), help: 'Side of the parallelogram facet along the zigzag creases.' },
  { key: 'miB', label: 'Miura parallelogram side b', unit: 'm', value: 1e-2, min: 1e-6, max: 1, log: true, group: 'Mechanism', showIf: BIS('miura'), help: 'Side of the parallelogram facet along the straight creases.' },
  { key: 'miGamma', label: 'Miura sector angle γ', unit: '°', value: 60, min: 15, max: 85, group: 'Mechanism', showIf: BIS('miura'), help: 'Acute angle of the parallelogram facet. It fixes the magnitude of the in-plane Poisson’s ratio.' },
  { key: 'miTheta', label: 'Miura fold angle', unit: '°', value: 45, min: 2, max: 88, group: 'Mechanism', showIf: BIS('miura'), help: 'Dihedral angle between the facets and the base plane: 0° is flat (unfolded), 90° is fully folded.' },
  { key: 'kiLc', label: 'Kirigami cut length', unit: 'm', value: 1e-2, min: 1e-6, max: 1, log: true, group: 'Mechanism', showIf: BIS('kirigami'), help: 'Length of each slit, perpendicular to the pulling direction.' },
  { key: 'kiD', label: 'Kirigami row spacing', unit: 'm', value: 1e-3, min: 1e-7, max: 0.5, log: true, group: 'Mechanism', showIf: BIS('kirigami'), help: 'Distance between successive rows of slits along the pulling direction (width of each ligament).' },
  { key: 'kiG', label: 'Kirigami hinge width', unit: 'm', value: 1e-3, min: 1e-7, max: 0.5, log: true, group: 'Mechanism', showIf: BIS('kirigami'), help: 'Uncut material left between collinear slits; limited to 80 % of the cut length.' },
  { key: 'nteType', label: 'Cell', type: 'select', value: 'triangle', group: 'Mechanism', showIf: SH('negcte'), options: [{ value: 'triangle', label: 'Bi-material triangle lattice (tunable / negative thermal expansion)' }, { value: 'winerack', label: 'Wine-rack cell (negative linear compressibility)' }], help: 'Pin-jointed bi-material triangle whose height expansion is set by geometry, or a hinged rhombic cell whose linear compressibility changes sign.' },
  { key: 'alphaB', label: 'CTE of the base member', unit: '1/K', value: 2.3e-5, min: -1e-5, max: 3e-4, group: 'Mechanism', showIf: NTE('triangle'), help: 'Thermal expansion coefficient of the base (high-expansion) member of the triangle; default aluminium.' },
  { key: 'alphaS', label: 'CTE of the inclined members', unit: '1/K', value: 8.6e-6, min: -1e-5, max: 3e-4, group: 'Mechanism', bus: 'prop.cte', showIf: NTE('triangle'), help: 'Thermal expansion coefficient of the two inclined (low-expansion) members; default titanium.' },
  { key: 'nteTheta', label: 'Base angle of the triangle', unit: '°', value: 45, min: 10, max: 80, group: 'Mechanism', showIf: NTE('triangle'), help: 'Angle between the base and the inclined members. Shallow triangles amplify the base expansion and drive the height expansion negative.' },
  { key: 'wrTheta', label: 'Wine-rack strut angle', unit: '°', value: 30, min: 5, max: 85, group: 'Mechanism', showIf: NTE('winerack'), help: 'Angle between the struts and the x-axis. Below 45° the cell lengthens along x under hydrostatic pressure.' },
  { key: 'gradVar', label: 'Graded variable', type: 'select', value: 'strut', group: 'Grading', showIf: SH('graded'), options: [{ value: 'strut', label: 'Strut size (relative density)' }, { value: 'angle', label: 'Cell angle θ' }], help: 'Quantity varied from cell to cell across the graded array.' },
  { key: 'gradLaw', label: 'Grading law', type: 'select', value: 'linear', group: 'Grading', showIf: SH('graded'), options: [{ value: 'linear', label: 'Linear' }, { value: 'power', label: 'Quadratic' }, { value: 'exp', label: 'Exponential' }], help: 'Shape of the variation between the first and the last cell.' },
  { key: 'gradRatio', label: 'Strut ratio last/first', unit: '–', value: 2, min: 0.3, max: 3, group: 'Grading', showIf: SH('graded'), help: 'Strut size of the last cell divided by that of the first when the strut is graded.' },
  { key: 'theta2', label: 'Cell angle of the last cell', unit: '°', value: 60, min: 2, max: 75, group: 'Grading', showIf: SH('graded'), help: 'End value of the cell angle when the angle is graded (the first cell uses the cell angle θ).' },
  { key: 'nGrad', label: 'Cells across the graded array', unit: '–', value: 10, min: 3, max: 24, step: 1, group: 'Grading', sens: false, showIf: SH('graded'), help: 'Number of cells over which the property is graded.' },
  { key: 'nAcross', label: 'Largest finite width studied', unit: 'cells', value: 6, min: 2, max: 8, step: 1, group: 'Grading', sens: false, showIf: SH('graded'), help: 'The edge-effect study solves strips of 1 to this many cells across with free edges.' },
];
function clean(inp) {
  const v = {};
  for (const i of inputs) {
    const t = i.type || 'number', x = inp ? inp[i.key] : undefined;
    if (t === 'number') { const q = Number(x); v[i.key] = Number.isFinite(q) ? min(i.max, max(i.min, q)) : i.value; }
    else if (t === 'select') v[i.key] = i.options.some(o => o.value === x) ? x : i.value;
    else v[i.key] = x === undefined || x === null ? i.value : String(x);
  }
  for (const k of ['nView', 'nGrad', 'nAcross']) v[k] = round(v[k]);
  return v;
}
const optLabel = (key, val) => inputs.find(i => i.key === key).options.find(o => o.value === val).label;
const printable = (feature, v) => { const p = PROC[v.process]; return { p, ok: feature >= p.res, margin: feature / p.res }; };
function latticeGeometry(type, prm, a, t, n) {
  const [u, s] = unitOf(a), def = cellDef(type), n3 = def.dim === 3 ? min(n, 3) : n;
  const g = OWN[type] ? tessellate(def.make(prm), def.dim, n3, a, t, s) : implicit.buildLattice({ type, cells: [n3, n3, def.dim === 3 ? n3 : 1], cellSize_mm: a * s, strut_mm: t * s, params: prm });
  return { geo: { kind: 'lattice', nodes: g.nodes, edges: g.edges, radius: g.radius, units: u }, note: g.note || null, scale: s, n: n3 };
}
const latticeExports = (geo, desc) => [
  { name: 'lattice_nodes.csv', mime: 'text/csv', content: `id,x_${geo.units},y_${geo.units},z_${geo.units}\n` + geo.nodes.map((p, i) => `${i},${+p[0].toPrecision(7)},${+p[1].toPrecision(7)},${+(p[2] || 0).toPrecision(7)}`).join('\n') + '\n', description: 'Node coordinates of the tessellated lattice.' },
  { name: 'lattice_edges.csv', mime: 'text/csv', content: 'id,node_i,node_j\n' + geo.edges.map((e, i) => `${i},${e[0]},${e[1]}`).join('\n') + '\n', description: 'Strut connectivity (node index pairs) of the tessellated lattice.' },
  { name: 'unit_cell.json', mime: 'application/json', content: JSON.stringify(desc, null, 1), description: 'Unit-cell descriptor: lattice vectors, periodic nodes and struts, section and homogenised properties.' },
];
const cellDescriptor = (c, extra) => ({ type: c.type, dimension: c.dim, lattice_periods_m: c.red.per.slice(0, c.dim), nodes_m: c.red.reps.map(p => p.slice(0, c.dim)), struts: c.red.els.map(e => ({ i: e.a, j: e.b, vector_m: e.d.slice(0, c.dim), length_m: e.L })), section: c.dim === 3 ? { shape: 'circular', diameter_m: c.t } : { shape: 'rectangular, unit depth', thickness_m: c.t }, stiffness_matrix_Pa: c.C, voigt_order: c.dim === 3 ? 'xx,yy,zz,yz,xz,xy' : 'xx,yy,xy', ...extra });

// ---------------------------------------------------------------- mechanical lattice cell
function runLattice(v) {
  const warnings = [], m = { Es: v.Es, nus: v.nus }, type = v.cellType, def = cellDef(type), dim = def.dim, prm = cellParams(type, v.theta, v.hl), a = v.cell;
  const c = cellModel(type, prm, a, v.strut, m), t = c.t, ref = reference(c, m), ax = dim === 3 ? 'z' : 'x', ay = dim === 3 ? 'x' : 'y';
  if (c.capped) warnings.push(`The strut size was limited to ${len(t)} (45 % of the shortest strut, ${len(c.Lmin)}); thicker members are not beams and the frame model does not apply.`);
  const rhoU = c.rho1 >= 0.03 ? unionFraction(c.red, t) : null, rho = rhoU === null ? c.rho1 : min(c.rho1, max(rhoU, 0.5 * c.rho1));   // below 3 % the overlap is within the sampling error
  const lg = latticeGeometry(type, prm, a, t, v.nView), geo = lg.geo; if (lg.note) warnings.push(lg.note + '.');
  const Mblock = geo.edges.length - dim * geo.nodes.length + (dim === 3 ? 6 : 3), Mcell = c.b - dim * c.j, slender = t / c.Lmin;
  // strength
  const sy = c.yieldRatio * v.ys, sb = c.buckle, govern = sb !== null && sb < sy ? 'elastic buckling of struts' : 'first yield of struts', strength = sb !== null ? min(sy, sb) : sy;
  // anisotropy
  const Cn = c.C, zener = dim === 3 ? 2 * Cn[3][3] / (Cn[0][0] - Cn[0][1]) : 2 * Cn[2][2] / (Cn[0][0] - Cn[0][1]);
  const phis = N.linspace(0, dim === 3 ? 90 : 180, dim === 3 ? 31 : 61), Edir = phis.map(p => { const q = rad(p); return dirE(c.S, dim, dim === 3 ? [sin(q) / Math.SQRT2, sin(q) / Math.SQRT2, cos(q)] : [cos(q), sin(q), 0]) / v.Es; });
  const anis = N.maxOf(Edir) / max(1e-300, N.minOf(Edir));
  // relative-density sweep (strut size varied) and scaling exponent
  const facs = N.logspace(0.3, 2.5, 9), sw = facs.map(f => cellModel(type, prm, a, t * f, m, true)).filter(q => !q.capped || q.t >= t), sr = sw.map(q => q.rho1), sE = sw.map(q => q.E / v.Es);
  const lowI = sr.map((r, i) => i).filter(i => sr[i] < 0.35), fitI = lowI.length >= 3 ? lowI : sr.map((r, i) => i), pl = fitI.length >= 2 ? N.linreg(fitI.map(i => log(sr[i])), fitI.map(i => log(max(1e-300, sE[i])))) : { slope: 1, intercept: 0 };
  const nexp = fin(pl.slope, 1), Cpl = exp(fin(pl.intercept, 0));
  // angle sweep for parametric cells
  const angleCell = ANGLE_CELLS.includes(type); let ang = null;
  if (angleCell) {
    const ths = N.linspace(type === 'rotatingSquares2D' ? 3 : 5, type === 'rotatingSquares2D' ? 43 : 70, 14), ms = ths.map(th => cellModel(type, cellParams(type, th, v.hl), a, t, m, true));
    ang = { ths, nu: ms.map(q => q.nu), nuo: ms.map(q => q.nuo), E: ms.map(q => q.E / v.Es), ref: ms.map(q => { const r = reference(q, m); return r && r.nuF !== undefined ? r.nuF : r ? r.nu : null; }) };
  }
  // tolerance sensitivity (central finite differences on the strut size and, for parametric cells, ±1° on the angle)
  const tol = v.tolPct / 100, cp = cellModel(type, prm, a, t * (1 + max(tol, 0.01)), m, true), cm = cellModel(type, prm, a, t * (1 - max(tol, 0.01)), m, true), sc = tol / max(tol, 0.01);
  const dE = (cp.E - cm.E) / (2 * c.E) * sc, dNu = (cp.nu - cm.nu) / 2 * sc, dRho = (cp.rho1 - cm.rho1) / (2 * c.rho1) * sc, sensE = (log(cp.E) - log(cm.E)) / (log(cp.t) - log(cm.t));
  const tolRows = [[`Strut size ±${F(v.tolPct, 3)} %`, +(dE * 100).toPrecision(4), +dNu.toPrecision(3), +(dRho * 100).toPrecision(4)]];
  if (angleCell) { const p1 = cellModel(type, cellParams(type, min(75, v.theta + 1), v.hl), a, t, m, true), m1 = cellModel(type, cellParams(type, max(2, v.theta - 1), v.hl), a, t, m, true); tolRows.push(['Cell angle ±1°', +((p1.E - m1.E) / (2 * c.E) * 100).toPrecision(4), +((p1.nu - m1.nu) / 2).toPrecision(3), +((p1.rho1 - m1.rho1) / (2 * c.rho1) * 100).toPrecision(4)]); }
  // manufacturability
  const pr = printable(t, v), P = pr.p; let low = 0, horiz = 0; const angles = new Set();
  if (dim === 3) for (const e of c.red.els) { const el = deg(Math.asin(min(1, abs(e.d[2]) / e.L))); angles.add(round(el)); if (el < P.overhang - 0.5) low++; if (el < 1) horiz++; }
  const overhangOK = dim === 2 || P.overhang === 0 || low === 0;
  if (P.planar && dim === 3) warnings.push(`${P.name} is a planar patterning process and cannot build a three-dimensional strut lattice directly; a layer-by-layer or direct-write route is needed.`);
  if (slender > 0.2) warnings.push(`Strut slenderness t/L = ${F(slender, 3)} exceeds 0.2: beam theory overestimates compliance at the nodes and the frame results are approximate; use solid finite elements for final values.`);
  if (!c.sol.ok) warnings.push('The periodic frame has a near-singular stiffness matrix (a mechanism); the corresponding moduli are numerically regularised.');
  // admissibility
  const Cs = Cn.map(r => r.map(x => x / v.Es)), ev = N.eigSym(Cs).values, recip = abs(c.nu / c.E - c.nuo / c.Eo) / max(1e-300, abs(c.nu / c.E), abs(c.nuo / c.Eo));
  const checks = [
    chk('Stiffness matrix positive definite (stable material)', ev[0] > 0, `smallest eigenvalue of C/E_s = ${F(ev[0], 3)}; Poisson ratios may be negative but the strain energy stays positive`),
    chk('Reciprocity of the compliance matrix', dim === 3 || recip < 1e-3 || abs(c.nu) < 1e-6, dim === 3 ? `ν_${ax}${ay}/E_${ax} symmetric by construction; anisotropy ratio ${F(anis, 3)}` : `ν_xy/E_x = ${F(c.nu / c.E, 3)} against ν_yx/E_y = ${F(c.nuo / c.Eo, 3)} 1/Pa`),
    chk('Slender-strut (beam) model valid, t/L ≤ 0.2', slender <= 0.2, `t/L = ${F(slender, 3)} for the shortest strut`),
    chk(`Minimum feature printable by ${P.name}`, pr.ok, `strut ${len(t)} against a typical minimum feature of ${len(P.res)} (ratio ${F(pr.margin, 3)})`),
    chk('Overhang / support requirement', overhangOK, dim === 2 ? 'prismatic 2-D cell: built along the extrusion axis, no overhanging struts' : P.overhang === 0 ? `${P.name}: struts are self-supporting at any angle; strut elevations ${[...angles].sort((p, q) => p - q).join('°, ')}°` : `${low} of ${c.b} struts per cell lie below ${P.overhang}° from the build plane (${horiz} horizontal); strut elevations ${[...angles].sort((p, q) => p - q).join('°, ')}° with z as build direction`),
  ];
  if (ref && ref.E) checks.push(chk('Frame homogenisation against the analytical model (modulus)', abs(c.E / v.Es / ref.E - 1) < 0.15 || slender > 0.2, `E_${ax}/E_s = ${F(c.E / v.Es, 4)} (frame) against ${F(ref.E, 4)} (${ref.label}); difference ${F((c.E / v.Es / ref.E - 1) * 100, 3)} %`));
  if (ref && ref.nu !== null && ref.nu !== undefined) checks.push(chk('Frame homogenisation against the analytical model (Poisson ratio)', abs(c.nu - ref.nu) < 0.05 + 0.15 * abs(ref.nu) || type === 'rotatingSquares2D', `ν = ${F(c.nu, 4)} (frame) against ${F(ref.nu, 4)} (${ref.label})`));
  if (type === 'chiral2D') warnings.push('The drawn and analysed cell is tetrachiral. The classical Prall–Lakes result ν = −1 applies to the hexachiral lattice with rigid rings; the tetrachiral cell is anisotropic and its Poisson ratio along the lattice axes is close to zero.');
  if (type === 'rotatingSquares2D') warnings.push('The squares are modelled as welded strut frames joined rigidly at their corners, so the frame result includes hinge bending; the ideal hinged mechanism gives ν = −1 exactly with a stiffness set only by the hinges.');
  if (c.nu < 0) warnings.push(`Negative Poisson's ratio (ν = ${F(c.nu, 3)}): the response is geometric, valid for small strains about this configuration, and changes as the cell angle evolves under load.`);
  const dom = c.stretch > 0.5 ? 'stretch-dominated' : 'bending-dominated', pretty = optLabel('cellType', type);
  const hexLine = ref && ref.nuF !== undefined && ref.nuF !== null ? ` The thin-wall Gibson–Ashby flexure value is ν₁₂ = ${F(ref.nuF, 4)} for θ = ${F(deg(ref.th), 3)}°, h/l = ${F(ref.hl, 3)}.` : '';
  const kpis = [
    kpi('rho', 'Relative density', rho, '–', `Solid fraction of the cell = ${F(rho, 4)}. The first-order strut sum ΣA·L/V gives ${F(c.rho1, 4)}; ${rhoU !== null ? `counting node overlaps once (quasi-random sampling of the union of the struts, about 1 % sampling error) gives ${F(rhoU, 4)}` : 'node overlap is negligible at this slenderness'}.${ref && ref.rho ? ` Closed form: ${F(ref.rho, 4)}.` : ''}`),
    kpi('Z', 'Nodal connectivity Z', c.Z, '–', `Average number of struts meeting at a node of the periodic cell (${c.b} struts, ${c.j} nodes) = ${F(c.Z, 3)}. Rigidity of similarly situated pin-jointed nodes needs Z ≥ ${dim === 3 ? 12 : 6} (${dim === 3 ? '3-D' : '2-D'}).`),
    kpi('maxwell', 'Maxwell number', Mblock, '–', `M = b − ${dim}j + ${dim === 3 ? 6 : 3} = ${Mblock} for the displayed ${lg.n}-cell block (per periodic cell b − ${dim}j = ${Mcell}). M < 0 means mechanisms when pin-jointed, so stiffness must come from strut bending; M ≥ 0 is necessary but not sufficient for stretch domination.`, Mblock >= 0),
    kpi('stretch', 'Axial (stretch) share of stiffness', c.stretch, '–', `Modulus of the same cell with strut bending and torsion removed, divided by the full value = ${F(c.stretch, 3)}: the cell is ${dom} along ${ax}. The fitted scaling exponent E ∝ ρ̄ⁿ is n = ${F(nexp, 3)}.`, c.stretch > 0.5),
    kpi('Erel', `Relative modulus E_${ax}/E_s`, c.E / v.Es, '–', `Periodic Timoshenko-frame homogenisation: E_${ax} = ${F(c.E, 4)} Pa, that is ${F(c.E / v.Es, 4)} of the solid.${ref && ref.E ? ` ${ref.label}: ${F(ref.E, 4)}.` : ''} In the other principal direction E_${ay}/E_s = ${F(c.Eo / v.Es, 4)}.`),
    kpi('nu', `Poisson's ratio ν_${ax}${ay}`, c.nu, '–', `Lateral contraction along ${ay} per unit extension along ${ax} = ${F(c.nu, 4)}${c.nu < 0 ? ' (auxetic: the cell widens when stretched)' : ''}.${hexLine} The conjugate ratio is ${F(c.nuo, 4)}.`),
    kpi('Grel', 'Relative shear modulus G/E_s', c.G / v.Es, '–', `Shear modulus ${dim === 3 ? 'G_xz' : 'G_xy'} = ${F(c.G, 4)} Pa from the homogenised compliance.${ref && ref.G ? ` Analytical: ${F(ref.G, 4)} E_s.` : ''} Bulk-type modulus (${dim === 3 ? 'volumetric' : 'areal'}): ${F(c.K, 4)} Pa.`),
    kpi('E', `Effective modulus E_${ax}`, c.E, 'Pa', `E_s × ${F(c.E / v.Es, 4)} with E_s = ${F(v.Es, 4)} Pa; published as eff.E for the mechanics suites.`),
    kpi('yield', 'First-yield strength', sy, 'Pa', `Uniaxial stress along ${ax} at which the most highly stressed strut reaches σ_ys (axial plus bending fibre stress from the frame solution) = ${F(c.yieldRatio, 4)} σ_ys. It is a lower bound to plastic collapse.${ref && ref.plastic ? ` Analytical plastic-collapse estimate: ${F(ref.plastic, 4)} σ_ys.` : ''}`),
    kpi('buckle', 'Elastic-buckling limit (compression)', sb === null ? 'no strut in compression' : sb, sb === null ? '' : 'Pa', sb === null ? `No strut carries compression under uniaxial compression along ${ax} in the linear solution, so Euler buckling of struts does not limit this load case.` : `Compressive stress along ${ax} at which the most heavily loaded strut reaches its Euler load π²E_sI/L² (pinned ends) = ${F(sb, 4)} Pa. Rigid nodes raise it by up to 4; sway-permitting cells such as honeycombs can be lower.${ref && ref.buckle2 ? ` Gibson–Ashby honeycomb buckling for loading along y: ${F(ref.buckle2 * v.Es, 4)} Pa.` : ''} Governing failure: ${govern}.`, sb === null ? null : sb > sy),
    kpi('anis', 'Directional stiffness ratio E_max/E_min', anis, '–', `Ratio of the largest to the smallest Young's modulus over ${dim === 3 ? 'directions between [001], [111] and [110]' : 'all in-plane directions'} = ${F(anis, 4)}; the Zener-type ratio 2C_s/(C₁₁ − C₁₂) is ${F(zener, 4)} (1 for an isotropic cell).`, anis < 1.5),
    kpi('feature', 'Feature size / process resolution', pr.margin, '–', `Strut size ${len(t)} divided by the typical minimum feature of ${P.name} (${len(P.res)}) = ${F(pr.margin, 3)}. Values below about 2 give large relative geometric error.`, pr.margin >= 2),
    kpi('tolE', 'Stiffness change at tolerance', dE * 100, '%', `Change of E_${ax} for a ±${F(v.tolPct, 3)} % deviation of the strut size (central finite difference) = ±${F(abs(dE) * 100, 3)} %; the logarithmic sensitivity ∂ln E/∂ln t is ${F(sensE, 3)}.`, abs(dE) < 0.2),
  ];
  const series = [
    { key: 'E_rho', title: 'Relative modulus vs relative density (strut size varied)', xlabel: 'relative density ρ̄ (–)', ylabel: `E_${ax}/E_s (–)`, logx: true, logy: true, lines: [{ name: 'frame homogenisation', x: sr, y: sE }, { name: `fit ${F(Cpl, 3)}·ρ̄^${F(nexp, 3)}`, x: sr, y: sr.map(r => Cpl * r ** nexp) }, ...(ref && ref.E && ['octet', 'triangular2D', 'kagome2D', 'square2D', 'simpleCubic'].includes(type) ? [{ name: 'analytical (pin-jointed)', x: sr, y: sr.map(r => r * ref.E / c.rho1) }] : [])], explain: 'Slope near 1 indicates stretch-dominated behaviour; slopes of 2 (3-D) or 3 (2-D) indicate bending-dominated behaviour, which is far more compliant at low density.' },
    { key: 'E_dir', title: 'Directional Young’s modulus', xlabel: dim === 3 ? 'angle from [001] towards [110] (°)' : 'in-plane direction from x (°)', ylabel: 'E(direction)/E_s (–)', logy: anis > 30, lines: [{ name: 'E/E_s', x: phis, y: Edir }], explain: `Young's modulus for uniaxial loading in different directions, from the homogenised compliance matrix. ${dim === 3 ? '[111] lies at 54.7°.' : ''} A flat curve means an isotropic cell.` },
  ];
  if (ang) series.push(
    { key: 'nu_theta', title: `Poisson's ratio vs cell angle`, xlabel: 'cell angle θ (°)', ylabel: `ν_${ax}${ay} (–)`, lines: [{ name: 'frame homogenisation', x: ang.ths, y: ang.nu }, ...(ang.ref.every(q => q !== null && Number.isFinite(q)) ? [{ name: type === 'rotatingSquares2D' ? 'rigid-unit kinematics' : 'Gibson–Ashby flexure model', x: ang.ths, y: ang.ref }] : []), { name: 'conjugate ratio', x: ang.ths, y: ang.nuo }], explain: 'The Poisson ratio is set by the cell angle: re-entrant angles give negative values whose magnitude grows as the ribs approach the loading direction. The finite wall thickness reduces the magnitude relative to the thin-wall formula.' },
    { key: 'E_theta', title: 'Relative modulus vs cell angle', xlabel: 'cell angle θ (°)', ylabel: `E_${ax}/E_s (–)`, logy: true, lines: [{ name: 'frame homogenisation', x: ang.ths, y: ang.E }], explain: 'Stiffness along the primary axis at fixed cell size and strut size. Large negative Poisson ratios are bought with low stiffness in the conjugate direction.' });
  const eqs = [
    eqn('governing', 'Periodic homogenisation of the beam frame', 'u(x) = ε̄·x + ũ(x), ũ periodic;  K_uu ũ = −K_uε ε̄;  C* = (K_εε − K_εu K_uu⁻¹ K_uε)/V_cell', 'Timoshenko beam elements on the periodic graph of the cell; the macroscopic strain enters through the lattice translation carried by each strut.'),
    eqn('constitutive', 'Effective constants from the compliance S* = C*⁻¹', 'E_i = 1/S_ii;  ν_ij = −S_ji/S_ii;  G_ij = 1/S_kk (shear);  1/E(n) = Σ S_IJ m_I m_J'),
    eqn('formula', 'Relative density', dim === 3 ? 'ρ̄ = Σ πr²L / V_cell (first order);  octet: ρ̄ = 6√2π(r/l)² − 54.6(r/l)³ with node overlap' : 'ρ̄ = Σ tL / A_cell (first order);  honeycomb: ρ̄ = (t/l)(h/l + 2) / (2cosθ (h/l + sinθ))'),
    eqn('formula', 'Maxwell rigidity criterion', dim === 3 ? 'M = b − 3j + 6;  rigidity of similarly situated nodes: Z ≥ 12' : 'M = b − 2j + 3;  rigidity of similarly situated nodes: Z ≥ 6'),
    eqn('formula', 'Strut failure', 'σ_max = |N|/A + |M|c/I = σ_ys (first yield);  P_cr = π²E_sI/L² (Euler, pinned)'),
  ];
  if (ref && ref.th !== undefined) eqs.push(eqn('constitutive', 'Gibson–Ashby honeycomb (flexure model)', 'ν₁₂ = cos²θ / ((h/l + sinθ) sinθ);  E₁/E_s = (t/l)³ cosθ / ((h/l + sinθ) sin²θ);  E₂/E_s = (t/l)³ (h/l + sinθ)/cos³θ;  G₁₂/E_s = (t/l)³ (h/l + sinθ) / ((h/l)²(1 + 2h/l) cosθ)', 'θ < 0 (re-entrant) makes ν₁₂ and ν₂₁ negative; reciprocity ν₁₂E₂ = ν₂₁E₁.'));
  if (type === 'octet') eqs.push(eqn('constitutive', 'Octet truss (pin-jointed)', 'E/E_s = ρ̄/9;  G/E_s = ρ̄/12;  ν = 1/3;  Zener ratio 2'));
  if (type === 'triangular2D' || type === 'kagome2D') eqs.push(eqn('constitutive', 'Stretch-dominated planar lattices', 'E/E_s = ρ̄/3;  G/E_s = ρ̄/8;  ν = 1/3 (isotropic in the plane)'));
  if (type === 'rotatingSquares2D') eqs.push(eqn('constitutive', 'Rotating rigid squares', 'X = 2s(cos(φ) + sin(φ));  ν = −1 for all φ'));
  const nuTxt = c.nu < 0 ? 'auxetic' : 'conventional';
  return {
    summary: `The ${pretty} cell of ${len(a)} period with ${len(t)} ${dim === 3 ? 'struts' : 'walls'} has a relative density of ${F(rho, 3)} and is ${dom} (Z = ${F(c.Z, 3)}, axial share ${F(c.stretch * 100, 3)} %). Homogenisation gives E_${ax} = ${F(c.E / v.Es, 3)} E_s = ${F(c.E, 3)} Pa, a ${nuTxt} Poisson's ratio ν_${ax}${ay} = ${F(c.nu, 3)} and G = ${F(c.G / v.Es, 3)} E_s. First yield occurs at ${F(sy, 3)} Pa${sb !== null ? ` and strut buckling at ${F(sb, 3)} Pa, so ${govern} governs` : ''}. The strut is ${F(pr.margin, 3)} times the typical minimum feature of ${P.name}.`,
    kpis, series,
    tables: [
      { key: 'tol', title: 'Tolerance sensitivity (finite differences)', columns: ['Deviation', `ΔE_${ax} (%)`, `Δν_${ax}${ay} (absolute)`, 'Δρ̄ (%)'], rows: tolRows, explain: 'Change of the key effective properties for the stated manufacturing deviation, from re-solving the periodic frame with perturbed geometry.' },
      { key: 'C', title: 'Homogenised stiffness matrix C*/E_s', columns: ['', ...(dim === 3 ? ['xx', 'yy', 'zz', 'yz', 'xz', 'xy'] : ['xx', 'yy', 'xy'])], rows: Cs.map((r, i) => [(dim === 3 ? ['xx', 'yy', 'zz', 'yz', 'xz', 'xy'] : ['xx', 'yy', 'xy'])[i], ...r.map(x => +(abs(x) < 1e-14 ? 0 : x).toPrecision(4))]), explain: 'Effective stiffness in Voigt notation with engineering shear strains, normalised by the constituent modulus. Off-diagonal normal terms of negative sign indicate auxetic coupling.' },
    ],
    equations: eqs, geometry: geo, checks, warnings,
    recommendations: [
      rec('high', c.stretch > 0.5 ? `The cell is stretch-dominated: stiffness scales nearly linearly with density (n = ${F(nexp, 3)}), so it is efficient for stiff, light structures; check node strength, which governs once struts carry mainly axial load.` : `The cell is bending-dominated (n = ${F(nexp, 3)}): at ρ̄ = ${F(rho, 3)} it reaches only ${F(c.E / v.Es, 2)} E_s. Use it where compliance, energy absorption or a tailored Poisson ratio is wanted; for stiffness per mass choose an octet, kagome or triangular topology.`, 'Axial share of stiffness and density-scaling exponent from the frame solution.'),
      rec('high', sb !== null && sb < sy ? `Elastic buckling (${F(sb, 3)} Pa) precedes yield (${F(sy, 3)} Pa): thicken the struts (buckling scales as (t/L)², one power higher than yield) or shorten them with a smaller cell.` : `Yield (${F(sy, 3)} Pa) governs over strut buckling; a stronger constituent raises the strength directly, and filleted nodes reduce the bending stress peaks there.`, 'Member forces of the homogenised frame under uniaxial load.'),
      rec(pr.ok && overhangOK ? 'low' : 'high', pr.ok ? (overhangOK ? `The ${len(t)} strut is ${F(pr.margin, 3)} times the minimum feature of ${P.name} and all struts are self-supporting; ${pr.margin < 3 ? 'expect strut-diameter scatter of the order of the tolerance entered and compensate in the CAD.' : 'no geometric compensation should be needed.'}` : `${low} struts per cell fall below the ${P.overhang}° overhang limit of ${P.name}: rotate the cell relative to the build direction (for example build along [111]) or choose a self-supporting topology to avoid internal supports that cannot be removed.`) : `The ${len(t)} strut is below the ${len(P.res)} minimum feature of ${P.name}: enlarge the cell at constant t/a (properties are scale-free) or change to a finer process.`, 'Process resolution and overhang rules.'),
      rec('medium', `Hold the strut size within ±${F(v.tolPct, 3)} %: the stiffness changes by ±${F(abs(dE) * 100, 3)} % (sensitivity ${F(sensE, 3)}). Measure as-built struts by CT or SEM and feed the measured size back into this model.`, 'Finite-difference tolerance sensitivity.'),
      rec('low', `At ρ̄ = ${F(rho, 3)} the cell uses ${F((1 - rho) * 100, 3)} % less material than the solid; architected lattices made additively generate little scrap, so specify the lowest density that meets stiffness and strength and recycle unfused powder or uncured resin.`, 'Relative density; material and energy efficiency of additive manufacture.', true),
    ],
    publish: { 'geom.cell_size': a, 'geom.strut': t, 'geom.rel_density': rho, 'eff.E': c.E, 'eff.nu': c.nu, 'eff.G': c.G, 'eff.K': c.K, 'eff.density': rho * v.rhos, 'eff.yield': strength, 'eff.anisotropy': anis, 'geom.connectivity': c.Z },
    exports: latticeExports(geo, cellDescriptor(c, { relative_density: rho, E_primary_Pa: c.E, E_secondary_Pa: c.Eo, nu_primary: c.nu, nu_conjugate: c.nuo, G_Pa: c.G, primary_axis: ax })),
    numerics: { method: 'Periodic Timoshenko beam-frame homogenisation (direct stiffness method, static condensation of the fluctuation field); quasi-random sampled union volume for node overlap', dof: c.sol.dof, steps: sw.length + (ang ? ang.ths.length : 0) + 4, iterations: 1, residual: recip, notes: `${c.j} periodic nodes and ${c.b} struts per cell. Linear, small-strain response; nodes are rigid joints of zero size.` },
  };
}

// ---------------------------------------------------------------- TPMS
const TPMS_AREA = { gyroid: 3.0915, schwarzP: 2.3451, schwarzD: 3.8378, iwp: 3.4646 };   // minimal-surface area per unit cell, in units of a²
function voxSurface(g) {   // interface faces and mean chord lengths of solid (1) and void (0) along the three axes
  const { nx, ny, nz, data } = g, st = [1, nx, nx * ny], dims = [nx, ny, nz]; let faces = 0, sSum = 0, sCnt = 0, pSum = 0, pCnt = 0;
  for (let a = 0; a < 3; a++) { const b = (a + 1) % 3, c = (a + 2) % 3;
    for (let jb = 0; jb < dims[b]; jb++) for (let jc = 0; jc < dims[c]; jc++) { let p = jb * st[b] + jc * st[c], run = 0, prev = -1, first = true;
      for (let i = 0; i < dims[a]; i++, p += st[a]) { const s = data[p] ? 1 : 0; if (prev >= 0 && s !== prev) { faces++; if (!first) { if (prev) { sSum += run; sCnt++; } else { pSum += run; pCnt++; } } first = false; run = 0; } run++; prev = s; } } }
  return { faces, solidChord: sCnt ? sSum / sCnt : 0, poreChord: pCnt ? pSum / pCnt : 0 };
}
function hashinShtrikman(rho, Es, nus) {   // upper bounds for a porous isotropic solid
  const Ks = Es / (3 * (1 - 2 * nus)), Gs = Es / (2 * (1 + nus)), K = 4 * Gs * Ks * rho / (4 * Gs + 3 * Ks * (1 - rho)), G = Gs * rho * (9 * Ks + 8 * Gs) / (9 * Ks + 8 * Gs + 6 * (1 - rho) * (Ks + 2 * Gs));
  return { K, G, E: 9 * K * G / (3 * K + G), nu: (3 * K - 2 * G) / (2 * (3 * K + G)) };
}
function runTPMS(v) {
  const warnings = [], a = v.cell, sheet = v.tpmsMode === 'sheet', graded = abs(v.rdTo - v.rd) > 1e-6, res = 24, [u, s] = unitOf(a);
  const g = implicit.tpmsVoxels({ type: v.tpmsType, cells: [2, 2, graded ? 3 : 2], resolution: res, relDensity: v.rd, mode: v.tpmsMode, grading: graded ? { axis: 'z', from: v.rd, to: v.rdTo } : undefined, cellSize_mm: a * s });
  const uni = graded ? implicit.tpmsVoxels({ type: v.tpmsType, cells: [2, 2, 2], resolution: res, relDensity: v.rd, mode: v.tpmsMode }) : g;
  const h = a / uni.resolution, st = voxSurface(uni), nTot = uni.nx * uni.ny * uni.nz, rho = uni.relDensity, SvFace = st.faces * h * h * 2 / 3 / (nTot * h ** 3);   // 1/m, first-order face count
  // surface from the two-point statistics of one cell sampled at twice the resolution (13 directions, extrapolated to zero lag)
  const fine = implicit.tpmsVoxels({ type: v.tpmsType, cells: [1, 1, 1], resolution: 2 * res, relDensity: v.rd, mode: v.tpmsMode }), Sv = max(1e-30, implicit.twoPointSurface(fine.data, fine.nx, fine.ny, fine.nz) * fine.resolution / a);
  const wall = sheet ? 2 * rho / Sv : 4 * rho / Sv, pore = 4 * (1 - rho) / Sv, wallChord = st.solidChord * h, poreChord = st.poreChord * h, ssaMass = Sv / (rho * v.rhos);
  // level set ↔ relative density from the sampled field of one cell
  const fn = implicit.TPMS[v.tpmsType], ns = 22, vals = new Float64Array(ns ** 3); let p = 0;
  for (let k = 0; k < ns; k++) for (let j = 0; j < ns; j++) for (let i = 0; i < ns; i++) { const f = fn((i + 0.137) / ns * 2 * PI, (j + 0.291) / ns * 2 * PI, (k + 0.413) / ns * 2 * PI); vals[p++] = sheet ? abs(f) : f; }
  vals.sort(); const rs = N.linspace(0.02, 0.98, 49), levels = rs.map(r => vals[min(vals.length - 1, round(r * vals.length))]), level = N.interp1(rs, levels, rho), fmax = vals[vals.length - 1];
  // stiffness: empirical power law, capped by the Hashin–Shtrikman upper bound
  const pw = r => v.plC * r ** v.plN, hsAt = r => hashinShtrikman(r, v.Es, v.nus), hs = hsAt(rho), Erel = min(pw(rho), hs.E / v.Es), E = Erel * v.Es, G = E / (2 * (1 + hs.nu));
  // graded profile along z from the voxel slices
  const zs = [], rz = [], Ez = []; for (let k = 0; k < g.nz; k++) { let c = 0; for (let q = k * g.nx * g.ny, e = q + g.nx * g.ny; q < e; q++) c += g.data[q]; const r = c / (g.nx * g.ny); zs.push((k + 0.5) * g.voxel); rz.push(r); Ez.push(min(pw(r), hsAt(max(1e-6, min(1, r))).E / v.Es)); }
  const Eser = N.sum(Ez.map(e => 1 / max(1e-12, e))), Eseries = Ez.length / Eser, pr = printable(wall, v), A0 = TPMS_AREA[v.tpmsType], SvA = Sv * a, expect = A0 ? (sheet ? 2 * A0 : A0) : null;
  const areaComparable = expect && (sheet ? rho <= 0.35 : abs(rho - 0.5) < 0.2);
  const checks = [
    chk('Achieved relative density matches the target', abs(rho - v.rd) < 0.01, `achieved ${F(rho, 4)} against target ${F(v.rd, 4)} at ${uni.resolution} voxels per cell`),
    chk('Empirical stiffness below the Hashin–Shtrikman upper bound', pw(rho) <= hs.E / v.Es * 1.0001, `C·ρ̄ⁿ = ${F(pw(rho), 4)} against the bound ${F(hs.E / v.Es, 4)} (the smaller value is used)`),
    chk(`Wall printable by ${pr.p.name}`, pr.ok, `mean wall ${len(wall)} against a typical minimum feature of ${len(pr.p.res)}`),
    chk('Wall resolved by the voxel grid (≥ 2 voxels)', wall / h >= 2, `wall / voxel = ${F(wall / h, 3)}; thinner walls are under-resolved and the surface estimate degrades`),
    chk('Both phases continuous (self-supporting, drainable)', rho > 0.03 && rho < 0.97, sheet ? 'sheet TPMS: one connected wall separating two interpenetrating channel networks' : 'network TPMS: one solid labyrinth and one channel labyrinth'),
  ];
  if (areaComparable) checks.push(chk('Surface area against the minimal-surface reference', abs(SvA / expect - 1) < 0.2, `S_v·a = ${F(SvA, 4)} (two-point estimator) against ${sheet ? '2 × ' : ''}${A0} = ${F(expect, 4)} for the ${sheet ? 'two faces of a thin sheet' : 'mid-level surface'}`));
  if (wall / h < 2) warnings.push('The wall is thinner than two voxels: surface area and wall thickness are under-resolved. Raise the relative density or treat the values as indicative.');
  if (pore < 3 * pr.p.res && !pr.p.planar) warnings.push(`Channels of ${len(pore)} are less than three times the process feature size: trapped powder or resin will be hard to remove.`);
  if (pr.p.planar) warnings.push(`${pr.p.name} is a planar process and cannot produce a three-dimensional TPMS directly.`);
  if (pw(rho) > hs.E / v.Es) warnings.push('The chosen power-law coefficient and exponent exceed the Hashin–Shtrikman upper bound at this density; the bound was used instead.');
  const name = optLabel('tpmsType', v.tpmsType);
  return {
    summary: `A ${sheet ? 'sheet' : 'network'} ${name} of ${len(a)} cell size at relative density ${F(rho, 3)} corresponds to a level-set threshold of ${F(level, 3)}, has ${F(Sv, 3)} m²/m³ of surface (${F(SvA, 3)} per cell in units of 1/a), a mean wall of ${len(wall)} and channels of ${len(pore)}. The empirical power law gives E* = ${F(Erel, 3)} E_s = ${F(E, 3)} Pa (Hashin–Shtrikman bound ${F(hs.E / v.Es, 3)} E_s).${graded ? ` Graded from ${F(v.rd, 3)} to ${F(v.rdTo, 3)} along z, the through-thickness (series) modulus is ${F(Eseries, 3)} E_s.` : ''}`,
    kpis: [
      kpi('rho', 'Relative density', rho, '–', `Solid voxel fraction achieved = ${F(rho, 4)}; the threshold is set by the quantile of the sampled level-set field so the target is met to within one voxel count.`),
      kpi('level', 'Level-set threshold c', level, '–', `${sheet ? 'Solid where |f| ≤ c' : 'Solid where f ≤ c'}: c = ${F(level, 4)} on a field whose ${sheet ? 'magnitude' : 'value'} reaches ${F(fmax, 3)}. The curve of c against density is given below for other targets.`),
      kpi('Sv', 'Surface area per volume', Sv, '1/m', `Solid–void interface area per unit bulk volume = ${F(Sv, 4)} m²/m³ (two-point estimator S_v = −4·dS₂/dr along 13 lattice directions on one cell at ${fine.resolution} voxels per edge, extrapolated to zero lag; the first-order voxel face count × 2/3 gives ${F(SvFace, 4)}), i.e. ${F(SvA, 4)}/a; it scales inversely with cell size. Per unit solid mass: ${F(ssaMass, 4)} m²/kg.`),
      kpi('wall', sheet ? 'Mean wall thickness' : 'Mean strut diameter', wall, 'm', `${sheet ? '2ρ̄/S_v (solid volume over half the surface)' : '4ρ̄/S_v (equivalent cylinder)'} = ${len(wall)}; the mean solid chord along the axes is ${len(wallChord)}.`, pr.ok),
      kpi('pore', 'Channel (pore) size', pore, 'm', `Hydraulic diameter of the void, 4(1 − ρ̄)/S_v = ${len(pore)}; mean void chord ${len(poreChord)}. It governs permeability, cell ingrowth and powder removal.`),
      kpi('Erel', 'Relative modulus E*/E_s', Erel, '–', `Empirical power law C·ρ̄ⁿ with C = ${F(v.plC, 3)}, n = ${F(v.plN, 3)} gives ${F(pw(rho), 4)}; the Hashin–Shtrikman upper bound is ${F(hs.E / v.Es, 4)}. The power law is a fit form to be calibrated on tests, not a prediction from first principles.`),
      kpi('E', 'Effective modulus', E, 'Pa', `E_s × ${F(Erel, 4)} = ${F(E, 4)} Pa; effective density ${F(rho * v.rhos, 4)} kg/m³, so the specific stiffness is ${F(E / (rho * v.rhos), 4)} m²/s².`),
      kpi('Egrad', 'Through-thickness modulus of the graded stack', Eseries, 'E_s', `Series (Reuss) average of the slice moduli along z = ${F(Eseries, 4)} E_s${graded ? `; slices range from ${F(N.minOf(Ez), 3)} to ${F(N.maxOf(Ez), 3)} E_s` : ' (uniform structure: equal to the cell value up to voxel noise)'}.`),
    ],
    series: [
      { key: 'level', title: 'Level-set threshold vs relative density', xlabel: 'relative density ρ̄ (–)', ylabel: 'threshold c (–)', lines: [{ name: `${name}, ${sheet ? 'sheet' : 'network'}`, x: rs, y: levels }, { name: 'this design', x: [rho], y: [level] }], explain: 'Mapping between the threshold of the implicit function and the solid fraction, obtained from the cumulative distribution of the sampled field. Use it to set the threshold in CAD or slicer software.' },
      { key: 'E', title: 'Stiffness scaling with relative density', xlabel: 'relative density ρ̄ (–)', ylabel: 'E*/E_s (–)', logx: true, logy: true, lines: [{ name: `empirical C·ρ̄ⁿ`, x: rs, y: rs.map(pw) }, { name: 'Hashin–Shtrikman upper bound', x: rs, y: rs.map(r => hsAt(r).E / v.Es) }, { name: 'this design', x: [rho], y: [Erel] }], explain: 'The empirical law must stay below the rigorous upper bound for an isotropic porous solid. Sheet TPMS structures approach the bound more closely than strut lattices of equal density.' },
      { key: 'grad', title: 'Density and stiffness profile along z', xlabel: `z (${u})`, ylabel: 'ρ̄ and E*/E_s (–)', lines: [{ name: 'relative density', x: zs, y: rz }, { name: 'E*/E_s', x: zs, y: Ez }], explain: 'Slice-by-slice solid fraction of the generated voxel model and the corresponding local modulus. Within one cell the slice density oscillates with the surface; a graded design adds a trend.' },
    ],
    tables: [{ key: 'ref', title: 'Minimal-surface reference areas', columns: ['Surface', 'Area per cubic cell (a²)'], rows: Object.keys(TPMS_AREA).map(k => [optLabel('tpmsType', k), TPMS_AREA[k]]), explain: 'Area of one conventional cubic cell of the exact minimal surface. A thin sheet structure has about twice this area per cell volume; the voxel estimate carries an error of roughly ten per cent.' }],
    equations: [
      eqn('formula', `Level-set function — ${name}`, { gyroid: 'f = sin x cos y + sin y cos z + sin z cos x', schwarzP: 'f = cos x + cos y + cos z', schwarzD: 'f = sin x sin y sin z + sin x cos y cos z + cos x sin y cos z + cos x cos y sin z', iwp: 'f = 2(cos x cos y + cos y cos z + cos z cos x) − (cos 2x + cos 2y + cos 2z)', neovius: 'f = 3(cos x + cos y + cos z) + 4 cos x cos y cos z' }[v.tpmsType], 'x, y, z scaled by 2π/a. Leading Fourier terms of the minimal surface.'),
      eqn('formula', 'Solid definition and density', sheet ? 'solid: |f| ≤ c;  ρ̄(c) = P(|f| ≤ c)' : 'solid: f ≤ c;  ρ̄(c) = P(f ≤ c)'),
      eqn('formula', 'Wall and channel size', sheet ? 't_w = 2ρ̄/S_v;  d_h = 4(1 − ρ̄)/S_v;  S_v = (2/3)·N_faces·h²/V' : 'd_s = 4ρ̄/S_v;  d_h = 4(1 − ρ̄)/S_v;  S_v = (2/3)·N_faces·h²/V'),
      eqn('constitutive', 'Empirical stiffness scaling (Gibson–Ashby form)', 'E*/E_s = C ρ̄ⁿ  (C, n fitted to tests)'),
      eqn('constitutive', 'Hashin–Shtrikman upper bound (porous isotropic solid)', 'K_U = 4G_sK_sρ̄ / (4G_s + 3K_s(1 − ρ̄));  G_U = G_sρ̄(9K_s + 8G_s) / (9K_s + 8G_s + 6(1 − ρ̄)(K_s + 2G_s));  E_U = 9K_UG_U/(3K_U + G_U)'),
    ],
    geometry: { kind: 'voxels', nx: g.nx, ny: g.ny, nz: g.nz, data: g.data, voxel: g.voxel, units: u },
    checks, warnings,
    recommendations: [
      rec('high', pr.ok ? `The ${len(wall)} wall is ${F(pr.margin, 3)} times the minimum feature of ${pr.p.name}; ${pr.margin < 2 ? 'increase the cell size or density for a safer margin, because thin TPMS walls print with pinholes.' : 'the design is printable with margin.'}` : `The ${len(wall)} wall is below the ${len(pr.p.res)} minimum feature of ${pr.p.name}: at this density the cell size must be at least ${len(a * pr.p.res / wall * 1.5)} for a 1.5 × margin.`, 'Wall thickness 2ρ̄/S_v against process resolution.'),
      rec('medium', `Calibrate C and n on at least three densities of printed specimens; the defaults (${F(v.plC, 3)}, ${F(v.plN, 3)}) are only a starting point and published fits vary with process and wall quality.`, 'Empirical nature of the power law.'),
      rec('medium', graded ? `The graded stack is governed in series by its lightest slices (E = ${F(N.minOf(Ez), 3)} E_s); avoid grading below the printable wall thickness at the low-density end.` : `Grade the density (set a different far-end value) to place material where stresses are highest; the smooth level-set field keeps the surface continuous across the grading.`, 'Slice-wise density and stiffness profile.'),
      rec('low', `TPMS structures are self-supporting in most powder-bed and vat processes, so no support material is wasted; the ${len(pore)} open channels allow unfused powder or resin to be drained and reused.`, 'Continuous channel network; support-free manufacture.', true),
    ],
    publish: { 'geom.cell_size': a, 'geom.rel_density': rho, 'geom.strut': wall, 'geom.porosity': 1 - rho, 'geom.pore_size': pore, 'geom.ssa': ssaMass, 'eff.E': E, 'eff.nu': hs.nu, 'eff.G': G, 'eff.K': min(hs.K, E / (3 * (1 - 2 * hs.nu))), 'eff.density': rho * v.rhos },
    exports: [
      { name: `tpms_${v.tpmsType}_${g.nx}x${g.ny}x${g.nz}_uint8.raw`, mime: 'application/octet-stream', bytes: new Uint8Array(g.data), description: `Voxel model (1 = solid), x fastest, voxel ${F(g.voxel, 4)} ${u}.` },
      { name: 'tpms_design.json', mime: 'application/json', content: JSON.stringify({ type: v.tpmsType, mode: v.tpmsMode, cell_size_m: a, relative_density: rho, level_set_threshold: level, surface_per_volume_per_m: Sv, wall_m: wall, pore_m: pore, E_rel: Erel, grading: graded ? { axis: 'z', from: v.rd, to: v.rdTo } : null, dims: [g.nx, g.ny, g.nz], voxel_m: g.voxel / s }, null, 1), description: 'Design parameters and descriptors of the TPMS structure.' },
    ],
    numerics: { method: 'Level-set voxelisation with quantile thresholding; two-point surface-area estimator (13 lattice directions, extrapolated to zero lag) on one cell at twice the resolution; empirical power law bounded by Hashin–Shtrikman', dof: g.nx * g.ny * g.nz, steps: 1, iterations: 1, residual: abs(rho - v.rd), notes: `${uni.resolution} voxels per cell for the model, ${fine.resolution} for the surface estimate; surface-area error about 1–2 % while the wall spans at least four fine voxels, larger for thinner walls.` },
  };
}

// ---------------------------------------------------------------- resonant unit cells
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cdiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1] || 1e-300; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
const csqrtUp = z => { const r = hypot(z[0], z[1]); let re = sqrt(max(0, (r + z[0]) / 2)), im = sqrt(max(0, (r - z[0]) / 2)); if (z[1] < 0) re = -re; return im < 0 || (im === 0 && re < 0) ? [-re, -im] : [re, im]; };   // root with Im ≥ 0 (passive medium, e^(−iωt))
const rectP = (x0, y0, x1, y1, layer = 1) => ({ layer, pts: [[x0, y0], [x1, y0], [x1, y1], [x0, y1]] });
const negBand = (fs, ys) => { let lo = null, hi = null; for (let i = 0; i < fs.length; i++) if (ys[i] < 0) { if (lo === null) lo = fs[i]; hi = fs[i]; } return lo === null ? null : [lo, hi]; };
const FLUIDS = { air: { name: 'air', c: 343, rho: 1.204, mu: 1.81e-5 }, water: { name: 'water', c: 1481, rho: 998, mu: 1.0e-3 } };

function resEM(v, warnings) {
  const a = v.cell, type = v.resType, hasR = type !== 'wire', hasW = type !== 'srr', eh = v.epsH, c0 = N.C0, mu0 = N.MU0, e0 = N.EPS0;
  const size = min(v.srrSize, 0.95 * a), w = min(v.srrW, 0.3 * size), g = min(v.srrGap, 0.5 * size), tm = v.srrT, r = min(v.wireR, 0.2 * a);
  if (hasR && (size < v.srrSize || w < v.srrW || g < v.srrGap)) warnings.push(`Ring dimensions were limited to fit the cell: side ${len(size)}, trace ${len(w)}, gap ${len(g)}.`);
  if (hasW && r < v.wireR) warnings.push(`Wire radius limited to 20 % of the cell size (${len(r)}).`);
  // split ring: self-inductance of a square loop of rectangular conductor, gap + fringing + added capacitance
  const lm = size - w, req = 0.2235 * (w + tm), Lraw = 2 * mu0 * lm / PI * (log(lm / req) - 0.774), L = max(Lraw, 0.2 * mu0 * lm), Cgap = e0 * eh * w * tm / g, Csurf = 2 * e0 * eh * (w + tm) / PI * max(0.1, log(4 * lm / g)), C = Cgap + Csurf + v.srrCadd;
  if (hasR && Lraw < 0.2 * mu0 * lm) warnings.push('The conductor is too thick relative to the ring for the thin-loop inductance formula; a floor value was used and the resonance is only indicative.');
  const f0 = 1 / (2 * PI * sqrt(L * C)), w0 = 2 * PI * f0, skin = sqrt(2 / (w0 * mu0 * v.sigma)), Rr = (4 * lm - g) / (v.sigma * w * min(tm, 2 * skin)), Q = w0 * L / Rr, Fm = min(0.9, mu0 * lm ** 4 / (L * a ** 3)), fmp = f0 / sqrt(1 - Fm);
  // wire medium
  const lnP = log(a / r), lnB = log(a / (2 * PI * r)) + 0.5275, fpP = c0 / (a * sqrt(eh)) / sqrt(2 * PI * lnP), fp = c0 / (a * sqrt(eh)) / sqrt(2 * PI * lnB), gE = e0 * a * a * (2 * PI * fp) ** 2 / (PI * r * r * v.sigma) / (2 * PI);   // collision frequency in Hz
  const muAt = f => (hasR ? (() => { const q = cdiv([Fm * f * f, 0], [f * f - f0 * f0, f * f0 / Q]); return [1 - q[0], -q[1]]; })() : [1, 0]);
  const epAt = f => (hasW ? (() => { const q = cdiv([fp * fp, 0], [f * f, f * gE]); return [eh * (1 - q[0]), -eh * q[1]]; })() : [eh, 0]);
  const fref = hasR ? f0 : fp, fs = N.linspace(0.3 * fref, (hasR && hasW ? max(2.2 * f0, 1.3 * fp) : 2.2 * fref), 360), MU = fs.map(muAt), EP = fs.map(epAt), NN = fs.map((f, i) => csqrtUp(cmul(MU[i], EP[i])));
  const bMu = hasR ? negBand(fs, MU.map(z => z[0])) : null, bEp = hasW ? [0, fp] : null, bDng = hasR && hasW ? negBand(fs, MU.map((z, i) => (z[0] < 0 && EP[i][0] < 0 ? -1 : 1))) : null;
  const fe = bDng ? sqrt(bDng[0] * bDng[1]) : type === 'wire' ? 0.7 * fp : bMu ? sqrt(bMu[0] * bMu[1]) : 1.02 * f0, mu = muAt(fe), ep = epAt(fe), nn = csqrtUp(cmul(mu, ep)), aLam = a * fe * sqrt(eh) / c0, fom = abs(nn[0]) / max(1e-12, nn[1]);
  if (hasR && !bMu) warnings.push(`The ring loss is too high (Q = ${F(Q, 3)}) for the real part of μ to become negative: raise the conductivity or thickness, or increase the filling factor.`);
  if (type === 'srrwire' && !bDng) warnings.push(`No double-negative band: the magnetic resonance band (${freq(f0)}–${freq(fmp)}) must lie below the wire plasma frequency (${freq(fp)}) and μ must turn negative. Lower the ring resonance (more gap capacitance, larger ring) or thin the wires less.`);
  // layout
  const polys = []; if (hasR) for (const q of layout.resonator({ type: 'srr', params: { size, width: w, gap: g, shape: 'square' } })) polys.push(q);
  if (hasW) polys.push(rectP(-r, -a / 2, r, a / 2, 2));
  const band = type === 'srr' ? bMu : type === 'wire' ? bEp : bDng, bandName = type === 'srr' ? 'negative-μ band' : type === 'wire' ? 'negative-ε band' : 'double-negative band';
  const kpis = [];
  if (hasR) kpis.push(
    kpi('f0', 'Magnetic resonance f₀', f0, 'Hz', `1/(2π√(LC)) with L = ${F(L * 1e9, 4)} nH and C = ${F(C * 1e15, 4)} fF (gap ${F(Cgap * 1e15, 3)}, fringing ${F(Csurf * 1e15, 3)}, added ${F(v.srrCadd * 1e15, 3)} fF) = ${freq(f0)}. Lumped estimates of this kind are typically within 10–20 % of full-wave results.`),
    kpi('Q', 'Quality factor (ohmic)', Q, '–', `ω₀L/R with the ring resistance ${F(Rr, 3)} Ω limited by the skin depth ${len(skin)} = ${F(Q, 4)}. Radiation and dielectric loss are not included, so real values are lower.`, Q > 20),
    kpi('F', 'Magnetic filling factor F', Fm, '–', `μ₀A²/(L·V_cell) with loop area ${F(lm * lm * 1e6, 3)} mm² and cell volume a³ = ${F(Fm, 4)}; it sets the strength and width of the magnetic response.`),
    kpi('fmp', 'Magnetic plasma frequency', fmp, 'Hz', `f₀/√(1 − F) = ${freq(fmp)}: upper edge of the lossless negative-μ band, giving a fractional bandwidth of ${F((fmp - f0) / f0 * 100, 3)} %.`));
  if (hasW) kpis.push(kpi('fp', 'Wire-medium plasma frequency', fp, 'Hz', `Effective plasma frequency of the wire array, (c/a)/√(2π(ln(a/(2πr)) + 0.5275)) = ${freq(fp)} (the simpler ln(a/r) form gives ${freq(fpP)}). Below it Re ε < 0 for fields along the wires. Collision frequency ${freq(gE)}.`), kpi('fp_pendry', 'Wire-medium plasma frequency (leading-order logarithm)', fpP, 'Hz', `Pendry's original estimate (c/a)/√(2π ln(a/r)) = ${freq(fpP)}. It keeps only the leading logarithm of the wire inductance; the refined lattice sum used for the main value replaces ln(a/r) by ln(a/(2πr)) + 0.5275 and is ${F(100 * (fp / fpP - 1), 3)} % higher here. Quoted for comparison with values calculated from the original formula.`));
  kpis.push(
    kpi('band_lo', `Lower edge of the ${bandName}`, band ? band[0] : 0, 'Hz', band ? `From the complex effective parameters including loss: the ${bandName} starts at ${freq(band[0])}.` : `No ${bandName} exists for this geometry; zero is reported.`, !!band),
    kpi('band_hi', `Upper edge of the ${bandName}`, band ? band[1] : 0, 'Hz', band ? `The ${bandName} ends at ${freq(band[1])}; relative width ${F((band[1] - band[0]) / max(1e-300, (band[1] + band[0]) / 2) * 100, 3)} %.` : `No ${bandName} exists for this geometry; zero is reported.`, !!band),
    kpi('mu', 'Re μ_eff at the evaluation frequency', mu[0], '–', `Effective relative permeability at ${freq(fe)} (centre of the band): ${F(mu[0], 4)} + ${F(mu[1], 3)}i. A positive imaginary part means absorption (passive medium).`),
    kpi('eps', 'Re ε_eff at the evaluation frequency', ep[0], '–', `Effective relative permittivity at ${freq(fe)}: ${F(ep[0], 4)} + ${F(ep[1], 3)}i${hasW ? ' (Drude response of the wire array)' : ' (host only: no wires)'}.`),
    kpi('n', 'Re n_eff at the evaluation frequency', nn[0], '–', `n = √(εμ) on the branch with Im n ≥ 0: ${F(nn[0], 4)} + ${F(nn[1], 3)}i. ${nn[0] < 0 ? 'The real part is negative: phase and energy velocities are antiparallel.' : mu[0] * ep[0] < 0 ? 'ε and μ have opposite signs, so the wave is evanescent (stop band).' : 'Both parameters are positive.'} Figure of merit |Re n|/Im n = ${F(min(fom, 1e6), 3)}.`, type === 'srrwire' ? nn[0] < 0 : null),
    kpi('aLam', 'Cell size / wavelength a/λ', aLam, '–', `Lattice constant over the host wavelength at ${freq(fe)} = ${F(aLam, 3)}. Homogenisation is reasonable below about 0.1 and fails above about 0.25, where the cell diffracts.`, aLam < 0.25));
  const series = [
    { key: 'mu', title: 'Effective permeability and permittivity', xlabel: 'frequency (GHz)', ylabel: 'Re μ, Re ε (–)', lines: [{ name: 'Re μ', x: fs.map(f => f / 1e9), y: MU.map(z => max(-20, min(20, z[0]))) }, { name: 'Im μ', x: fs.map(f => f / 1e9), y: MU.map(z => min(20, z[1])) }, { name: 'Re ε', x: fs.map(f => f / 1e9), y: EP.map(z => max(-20, min(20, z[0]))) }], explain: 'Lorentz magnetic response of the ring and Drude electric response of the wires (clipped to ±20). The band where both real parts are negative supports backward-wave propagation.' },
    { key: 'n', title: 'Effective refractive index', xlabel: 'frequency (GHz)', ylabel: 'n (–)', lines: [{ name: 'Re n', x: fs.map(f => f / 1e9), y: NN.map(z => max(-10, min(10, z[0]))) }, { name: 'Im n', x: fs.map(f => f / 1e9), y: NN.map(z => min(10, z[1])) }], explain: 'Complex index on the passive branch (Im n ≥ 0). A negative real part with a small imaginary part marks the usable negative-index band; a large imaginary part alone marks a stop band.' },
  ];
  const eqs = [];
  if (hasR) eqs.push(eqn('formula', 'Split-ring L–C estimate', 'L ≈ (2μ₀l_m/π)[ln(l_m/ρ) − 0.774], ρ = 0.2235(w + t);  C = ε₀ε_h w t/g + (2ε₀ε_h(w + t)/π) ln(4l_m/g) + C_add;  f₀ = 1/(2π√(LC))', 'l_m: mean side of the square loop. The logarithmic fringing term is an estimate of the surface capacitance of the gap.'),
    eqn('constitutive', 'Effective permeability (Lorentz)', 'μ_eff(ω) = 1 − Fω² / (ω² − ω₀² + iγω),  F = μ₀A²/(L V_cell),  γ = R/L;  Re μ < 0 for ω₀ < ω < ω₀/√(1 − F)'));
  if (hasW) eqs.push(eqn('constitutive', 'Thin-wire medium (Drude)', 'ε_eff(ω) = ε_h[1 − ω_p² / (ω(ω + iγ_e))],  ω_p² = 2πc² / (ε_h a² [ln(a/(2πr)) + 0.5275]),  γ_e = ε₀a²ω_p²/(πr²σ)', 'Electric field parallel to the wires; spatial dispersion is neglected.'));
  eqs.push(eqn('governing', 'Refractive index of the homogenised medium', 'n = √(ε_eff μ_eff) with Im n ≥ 0;  Re n < 0 when Re ε < 0 and Re μ < 0'), eqn('formula', 'Homogenisation validity', 'a/λ = a f √ε_h / c ≲ 0.1–0.25'));
  const checks = [
    chk('Passivity: Im μ ≥ 0 and Im ε ≥ 0', MU.every(z => z[1] >= -1e-12) && EP.every(z => z[1] >= -1e-12), 'imaginary parts of the Lorentz and Drude responses are non-negative at all frequencies (e^(−iωt) convention)'),
    chk('Homogenisation validity a/λ < 0.25', aLam < 0.25, `a/λ = ${F(aLam, 3)} at ${freq(fe)}`),
    chk(`A ${bandName} exists`, !!band, band ? `${freq(band[0])} to ${freq(band[1])}` : 'none found over the scanned frequency range'),
  ];
  if (hasR) checks.push(chk('Thin-loop inductance formula applicable', lm / req > 5, `l_m/ρ = ${F(lm / req, 3)} (should exceed about 5)`));
  const cu = chem.prettyFormula('Cu');
  const recs = [
    rec('high', type === 'wire' ? `Operate below ${freq(fp)} for negative ε; ε reaches −1 near ${freq(fp / sqrt(2))}. To move the plasma frequency, change the lattice constant (f_p ∝ 1/a) — the wire radius acts only logarithmically.` : band ? `Design for operation near ${freq(fe)} inside the ${bandName} (${freq(band[0])}–${freq(band[1])}); the relative bandwidth of ${F((band[1] - band[0]) / fe * 100, 3)} % is intrinsic to the resonant mechanism and widens only with the filling factor F = ${F(Fm, 3)}.` : `Create the ${bandName} first: bring f₀ = ${freq(f0)} below ${hasW ? `the plasma frequency ${freq(fp)}` : 'the target band'} by adding gap capacitance or enlarging the ring, and raise Q.`, 'Lorentz/Drude effective-parameter model.'),
    rec('medium', aLam < 0.1 ? `a/λ = ${F(aLam, 3)}: the cell is deeply sub-wavelength and effective parameters are meaningful; verify with a full-wave S-parameter retrieval.` : `a/λ = ${F(aLam, 3)} is at the limit of homogenisation: lower the resonance (larger capacitance) at fixed cell size, and confirm by full-wave simulation because spatial dispersion will distort the retrieved parameters.`, 'Sub-wavelength criterion a/λ ≲ 0.1–0.25.'),
    rec('medium', hasR ? `Loss is set by the ring: Q = ${F(Q, 3)} (ohmic). Keep the metal thicker than two skin depths (${len(2 * skin)}) and use a low-loss substrate; dielectric and radiation loss will reduce the figure of merit below ${F(min(fom, 1e6), 3)}.` : `Wire loss is small (collision frequency ${freq(gE)}); the practical limit is the finite wire length and the junction quality between cells.`, 'Skin-depth resistance of the conductors.'),
    rec('low', `Printed ${cu} on a laminate needs only ${len(tm)} of metal: pattern by subtractive etching with etchant regeneration, or by additive printing, to minimise copper waste, and prefer halogen-free laminates for end-of-life recycling.`, 'Metal usage of planar resonant cells; circularity of PCB materials.', true),
  ];
  const pub = { 'geom.cell_size': a, 'eff.mu': mu[0], 'eff.eps': ep[0], 'eff.n': nn[0] };
  if (hasR) { pub['res.f0'] = f0; pub['res.Q'] = Q; } else pub['res.f0'] = fp;
  if (type === 'srr' && bMu) { pub['res.bandgap_lo'] = bMu[0]; pub['res.bandgap_hi'] = bMu[1]; }
  if (type === 'wire') { pub['res.bandgap_lo'] = 0; pub['res.bandgap_hi'] = fp; }
  if (type === 'srrwire') { if (bDng) { pub['res.passband_lo'] = bDng[0]; pub['res.passband_hi'] = bDng[1]; } pub['res.bandgap_lo'] = bDng ? bDng[1] : 0; pub['res.bandgap_hi'] = max(fp, bDng ? bDng[1] : 0); }
  const what = type === 'srr' ? 'split-ring' : type === 'wire' ? 'thin-wire' : 'split-ring + wire';
  return { polys, pitch: a, kpis, series, eqs, checks, recs, pub,
    summary: `The ${what} cell of ${len(a)} period ${hasR ? `resonates magnetically at ${freq(f0)} (L = ${F(L * 1e9, 3)} nH, C = ${F(C * 1e15, 3)} fF, Q = ${F(Q, 3)})` : ''}${hasR && hasW ? ' and ' : ''}${hasW ? `has a wire plasma frequency of ${freq(fp)}` : ''}. ${band ? `The ${bandName} spans ${freq(band[0])} to ${freq(band[1])}` : `No ${bandName} is obtained`}; at ${freq(fe)} the effective parameters are μ = ${F(mu[0], 3)}, ε = ${F(ep[0], 3)} and n = ${F(nn[0], 3)} + ${F(nn[1], 2)}i, with a/λ = ${F(aLam, 3)}.`,
    method: 'Lumped L–C (Lorentz) and Drude effective-medium models evaluated on a frequency grid', dof: fs.length };
}
function resAcoustic(v, warnings) {
  const a = v.cell, type = v.resType;
  if (type === 'helmholtz') {
    const fl = FLUIDS[v.fluid], V = v.hrV, rn = v.hrRn, Ln = v.hrLn, Sn = PI * rn * rn, Leff = Ln + 1.7 * rn, f0 = fl.c / (2 * PI) * sqrt(Sn / (V * Leff)), w0 = 2 * PI * f0, dv = sqrt(2 * fl.mu / (fl.rho * w0));
    const Q = max(0.5, Leff * rn / (dv * (Ln + 2 * rn))), Fg = min(5, V / a ** 3), fhi = f0 * sqrt(1 + Fg), B0 = fl.rho * fl.c ** 2, side = Math.cbrt(V), lam = fl.c / f0;
    const invB = f => { const q = cdiv([Fg * f0 * f0, 0], [f * f - f0 * f0, f * f0 / Q]); return [(1 - q[0]) / B0, -q[1] / B0]; };
    const fs = N.linspace(0.4 * f0, max(1.6 * fhi, 2 * f0), 320), IB = fs.map(invB), reB = IB.map(z => cdiv([1, 0], z)[0]), band = negBand(fs, IB.map(z => z[0])), fe = band ? sqrt(band[0] * band[1]) : 1.05 * f0, Be = cdiv([1, 0], invB(fe));
    const kq = fs.map((f, i) => { const k = csqrtUp(cmul([fl.rho * (2 * PI * f) ** 2, 0], IB[i])); return k[1] * a * 8.686; });   // attenuation per cell, dB
    if (rn < 3 * dv) warnings.push(`The neck radius (${len(rn)}) is less than three viscous boundary layers (${len(dv)}): the resonator is heavily damped and the lumped inertial model overestimates Q.`);
    if (V > a ** 3) warnings.push('The cavity is larger than the host cell volume a³: the cavities must sit outside the duct as side branches.');
    const cs = side, wall = 0.08 * cs;
    return { polys: layout.resonator({ type: 'helmholtzNeck2D', params: { cavityWidth: cs, cavityHeight: cs, wall, neckWidth: 2 * rn, neckLength: Ln } }), pitch: max(a, 1.5 * (cs + 2 * wall)),
      kpis: [
        kpi('f0', 'Helmholtz resonance f₀', f0, 'Hz', `(c/2π)√(S/(V·L_eff)) with neck area ${F(Sn * 1e6, 3)} mm² and effective length L + 1.7r = ${len(Leff)} (two flanged end corrections) = ${freq(f0)} in ${fl.name}.`),
        kpi('Q', 'Quality factor (viscous)', Q, '–', `L_eff·r/(δ_v(L + 2r)) with the viscous boundary layer δ_v = ${len(dv)} = ${F(Q, 4)}; thermal and radiation losses lower it further.`, Q > 10),
        kpi('F', 'Geometric factor F', Fg, '–', `Cavity volume over host volume per cell, V/a³ = ${F(Fg, 4)}; it sets the width of the negative-modulus band.`),
        kpi('band_lo', 'Negative-modulus band, lower edge', band ? band[0] : 0, 'Hz', band ? `Re(1/B_eff) < 0 from ${freq(band[0])}: the medium expands when compressed dynamically, so sound is evanescent.` : 'Damping prevents the real part of the inverse modulus from turning negative; zero is reported.', !!band),
        kpi('band_hi', 'Negative-modulus band, upper edge', band ? band[1] : 0, 'Hz', band ? `The band ends at ${freq(band[1])}; the lossless estimate is f₀√(1 + F) = ${freq(fhi)}.` : 'No negative-modulus band exists; zero is reported.', !!band),
        kpi('B', 'Re B_eff at mid-band', Be[0], 'Pa', `Effective bulk modulus at ${freq(fe)} = ${F(Be[0], 4)} Pa against ${F(B0, 4)} Pa for the host ${fl.name}. A negative value is dynamic (out-of-phase response of the resonator) and does not violate static stability.`),
        kpi('att', 'Peak attenuation per cell', N.maxOf(kq), 'dB', `Largest Im(k)·a × 8.686 over the scan = ${F(N.maxOf(kq), 4)} dB per cell from the effective-medium wavenumber k = ω√(ρ₀/B_eff).`),
        kpi('aLam', 'Cell size / wavelength a/λ', a / lam, '–', `a/λ at resonance = ${F(a / lam, 3)}; the cavity itself spans ${F(side / lam, 3)} λ. Both must be well below 0.25 for the lumped description.`, a / lam < 0.25),
      ],
      series: [
        { key: 'B', title: 'Effective bulk modulus', xlabel: 'frequency (Hz)', ylabel: 'Re B_eff / B₀ (–)', lines: [{ name: 'Re B_eff/B₀', x: fs, y: reB.map(b => max(-10, min(10, b / B0))) }], explain: 'Dynamic bulk modulus of the fluid loaded with resonators (clipped to ±10). It is negative just above the resonance, where the cavity breathes out of phase with the applied pressure.' },
        { key: 'att', title: 'Attenuation per unit cell', xlabel: 'frequency (Hz)', ylabel: 'attenuation (dB per cell)', lines: [{ name: 'Im(k)·a', x: fs, y: kq }], explain: 'Decay of a plane wave across one cell from the effective-medium wavenumber. The peak marks the stop band produced by the negative modulus.' },
      ],
      eqs: [eqn('formula', 'Helmholtz resonance with end corrections', 'f₀ = (c/2π)√(S_n / (V L_eff)),  L_eff = L + 2 × 0.85 r_n'), eqn('constitutive', 'Effective bulk modulus of the resonator-loaded fluid', '1/B_eff(ω) = (1/B₀)[1 − Fω₀² / (ω² − ω₀² + iΓω)],  F = V/V_cell,  Γ = ω₀/Q;  Re(1/B_eff) < 0 for ω₀ < ω < ω₀√(1 + F)'), eqn('governing', 'Wave propagation', 'k = ω√(ρ₀/B_eff);  attenuation per cell = Im(k)·a')],
      checks: [chk('Sub-wavelength cell and cavity', a / lam < 0.25 && side / lam < 0.25, `a/λ = ${F(a / lam, 3)}, cavity/λ = ${F(side / lam, 3)}`), chk('Passivity: Im(1/B_eff) ≥ 0', IB.every(z => z[1] >= -1e-30), 'the imaginary part of the compliance is non-negative at all frequencies (energy is dissipated, not generated)'), chk('Neck wider than three viscous boundary layers', rn >= 3 * dv, `r_n/δ_v = ${F(rn / dv, 3)}`), chk('A negative-modulus band exists', !!band, band ? `${freq(band[0])} to ${freq(band[1])}` : 'none')],
      recs: [
        rec('high', band ? `Use the band ${freq(band[0])}–${freq(band[1])} for sound blocking; widen it by raising V/a³ (now ${F(Fg, 3)}) — bandwidth grows as √(1 + F).` : 'Increase the neck radius or cavity volume ratio: the resonance is over-damped and no negative-modulus band forms.', 'Lumped effective-modulus model.'),
        rec('medium', `Tune the resonance with the neck: f₀ ∝ r_n/√(L_eff·V). A ±5 % error in neck radius shifts f₀ by about ±${F(5 * (1 - 0.85 * rn / Leff), 2)} %, so specify the neck tolerance first.`, 'Sensitivity of the Helmholtz formula.'),
        rec('medium', `Stagger several resonators with f₀ spaced by about the single-cell bandwidth (${F((fhi / f0 - 1) * 100, 3)} %) to cover a broader band.`, 'Additivity of resonator compliances at low filling.'),
        rec('low', 'Resonant liners achieve low-frequency attenuation with a fraction of the thickness and mass of porous absorbers, saving material; make them from a single recyclable polymer or metal to avoid mixed-material waste.', 'Sub-wavelength attenuation per unit mass; design for recycling.', true),
      ],
      pub: { 'geom.cell_size': a, 'res.f0': f0, 'res.Q': Q, 'eff.K': Be[0], 'eff.rho_dyn': fl.rho, ...(band ? { 'res.bandgap_lo': band[0], 'res.bandgap_hi': band[1] } : {}) },
      summary: `The Helmholtz resonator (cavity ${F(V * 1e9, 3)} mm³, neck radius ${len(rn)}, effective neck length ${len(Leff)}) resonates at ${freq(f0)} in ${fl.name} with a viscous Q of ${F(Q, 3)}. Arrayed with a ${len(a)} period (F = ${F(Fg, 3)}) it gives ${band ? `a negative effective bulk modulus from ${freq(band[0])} to ${freq(band[1])}` : 'no negative-modulus band because of damping'}, with a peak attenuation of ${F(N.maxOf(kq), 3)} dB per cell; a/λ = ${F(a / lam, 3)}.`,
      method: 'Lumped acoustic mass–compliance resonator with boundary-layer damping; effective-medium wavenumber', dof: fs.length };
  }
  if (type === 'massinmass') {
    const m1 = v.mimM1, m2 = v.mimM2, k1 = v.mimK1, k2 = v.mimK2, f2 = sqrt(k2 / m2) / (2 * PI), fhi = f2 * sqrt(1 + m2 / m1), eta = 0.02;
    const meff = f => { const q = cdiv([m2 * f2 * f2, 0], [f2 * f2 - f * f, -eta * f2 * f2]); return [m1 + q[0], q[1]]; };   // internal spring k₂(1 − iη) in the e^(−iωt) convention ⇒ Im m_eff ≥ 0
    const branch = K => { const w2 = (2 * PI * f2) ** 2, b = (m1 + m2) * w2 + K, d = sqrt(max(0, b * b - 4 * m1 * K * w2)); return [sqrt(max(0, (b - d) / (2 * m1))) / (2 * PI), sqrt((b + d) / (2 * m1)) / (2 * PI)]; };
    const qs = N.linspace(0, 1, 60), br = qs.map(q => branch(2 * k1 * (1 - cos(PI * q)))), flo = br[br.length - 1][0], fs = N.linspace(0.3 * f2, 1.6 * fhi, 320), ME = fs.map(meff), band = negBand(fs, ME.map(z => z[0]));
    const fe = sqrt(f2 * fhi), me = meff(fe), rhoS = (m1 + m2) / a ** 3, cL = a * sqrt(k1 / (m1 + m2)), aLam = a * flo / cL;
    const s = 0.8 * a, wf = 0.1 * a, si = 0.4 * a;
    return { polys: [rectP(-s / 2, -s / 2, s / 2, -s / 2 + wf), rectP(-s / 2, s / 2 - wf, s / 2, s / 2), rectP(-s / 2, -s / 2 + wf, -s / 2 + wf, s / 2 - wf), rectP(s / 2 - wf, -s / 2 + wf, s / 2, s / 2 - wf), rectP(-si / 2, -si / 2, si / 2, si / 2, 2)], pitch: a,
      kpis: [
        kpi('f0', 'Local resonance f₂', f2, 'Hz', `√(k₂/m₂)/2π = ${freq(f2)}: natural frequency of the internal mass on its spring with the shell held.`),
        kpi('band_lo', 'Band-gap lower edge', flo, 'Hz', `Top of the acoustic branch at the zone boundary (qa = π) from the two-mass dispersion relation = ${freq(flo)}.`),
        kpi('band_hi', 'Band-gap upper edge', fhi, 'Hz', `Bottom of the optical branch, f₂√(1 + m₂/m₁) = ${freq(fhi)}; between f₂ and this frequency the effective mass is negative.`),
        kpi('gap', 'Gap-to-mid-gap ratio', (fhi - flo) / ((fhi + flo) / 2), '–', `Relative width of the locally resonant gap = ${F((fhi - flo) / ((fhi + flo) / 2) * 100, 3)} %; it grows with the mass ratio m₂/m₁ = ${F(m2 / m1, 3)}.`),
        kpi('meff', 'Effective mass at the evaluation frequency', me[0], 'kg', `m_eff = m₁ + m₂ω₂²/(ω₂² − ω²) at ${freq(fe)} = ${F(me[0], 4)} kg against the static mass ${F(m1 + m2, 4)} kg. Negative means the cell accelerates against the applied force because the hidden resonator moves out of phase.`, me[0] < 0),
        kpi('rho', 'Effective dynamic density', me[0] / a ** 3, 'kg/m³', `m_eff/a³ = ${F(me[0] / a ** 3, 4)} kg/m³ (static ${F(rhoS, 4)} kg/m³); published as eff.rho_dyn.`),
        kpi('aLam', 'Cell size / wavelength at the gap', aLam, '–', `a·f/c_L with the long-wave speed a√(k₁/(m₁ + m₂)) = ${F(cL, 4)} m/s: ${F(aLam, 3)}. A locally resonant gap lies well below the Bragg condition a/λ = 0.5.`, aLam < 0.5),
      ],
      series: [
        { key: 'disp', title: 'Dispersion of the mass-in-mass chain', xlabel: 'reduced wavenumber qa/π (–)', ylabel: 'frequency (Hz)', lines: [{ name: 'acoustic branch', x: qs, y: br.map(b => b[0]) }, { name: 'optical branch', x: qs, y: br.map(b => b[1]) }], explain: 'Bloch dispersion of the infinite chain. No real wavenumber exists between the top of the acoustic branch and the bottom of the optical branch: this is the band gap.' },
        { key: 'meff', title: 'Effective mass vs frequency', xlabel: 'frequency (Hz)', ylabel: 'Re m_eff / (m₁ + m₂) (–)', lines: [{ name: 'Re m_eff', x: fs, y: ME.map(z => max(-10, min(10, z[0] / (m1 + m2)))) }], explain: 'Dynamic mass of one cell seen from outside (clipped to ±10, 2 % loss factor). It diverges at the local resonance and is negative up to f₂√(1 + m₂/m₁).' },
      ],
      eqs: [eqn('constitutive', 'Effective dynamic mass', 'm_eff(ω) = m₁ + m₂ω₂² / (ω₂² − ω²),  ω₂² = k₂/m₂;  m_eff < 0 for ω₂ < ω < ω₂√(1 + m₂/m₁)'), eqn('governing', 'Bloch dispersion of the chain', 'm_eff(ω) ω² = 2k₁(1 − cos qa)  ⇔  m₁ω⁴ − [(m₁ + m₂)ω₂² + K]ω² + Kω₂² = 0,  K = 2k₁(1 − cos qa)'), eqn('governing', 'Equations of motion', 'm₁ü₁ = k₁(u_{j+1} − 2u_j + u_{j−1}) + k₂(u₂ − u₁);  m₂ü₂ = −k₂(u₂ − u₁)')],
      checks: [chk('Gap edges ordered (f_lower < f₂ < f_upper)', flo < f2 * 1.0000001 && f2 < fhi, `${freq(flo)} < ${freq(f2)} < ${freq(fhi)}`), chk('Negative effective mass inside the gap', me[0] < 0, `m_eff = ${F(me[0], 3)} kg at ${freq(fe)}`), chk('Static limit recovered', abs(meff(1e-6 * f2)[0] / (m1 + m2) - 1) < 1e-3, `m_eff(ω → 0) = ${F(meff(1e-6 * f2)[0], 4)} kg = m₁ + m₂`), chk('Passivity: Im m_eff ≥ 0', ME.every(z => z[1] >= 0), 'positive loss factor in the internal spring gives a non-negative imaginary mass')],
      recs: [
        rec('high', `Place the working frequency between ${freq(flo)} and ${freq(fhi)}; to widen the gap raise m₂/m₁ (now ${F(m2 / m1, 3)}) — the upper edge moves as √(1 + m₂/m₁).`, 'Band edges of the two-mass chain.'),
        rec('medium', `Tune f₂ with the internal spring: a ±10 % stiffness tolerance shifts the gap by about ±5 %. Soft polymer coatings used as k₂ are temperature dependent, so check the gap over the service temperature range.`, 'f₂ = √(k₂/m₂)/2π.'),
        rec('medium', 'Add moderate damping in the internal spring to broaden the attenuation peak at the cost of depth; very high damping removes the negative-mass effect.', 'Complex effective mass with loss factor.'),
        rec('low', `Local resonance gives low-frequency isolation with a ${len(a)} cell instead of a Bragg structure many times larger, saving mass; avoid lead for the internal mass and use steel or tungsten-filled recyclable polymer.`, 'Sub-wavelength band gap; substitution of hazardous materials.', true),
      ],
      pub: { 'geom.cell_size': a, 'res.f0': f2, 'res.bandgap_lo': flo, 'res.bandgap_hi': fhi, 'eff.rho_dyn': me[0] / a ** 3, 'eff.density': rhoS },
      summary: `The mass-in-mass cell (m₁ = ${F(m1 * 1e3, 3)} g, m₂ = ${F(m2 * 1e3, 3)} g, k₁ = ${F(k1, 3)} N/m, k₂ = ${F(k2, 3)} N/m) has a local resonance at ${freq(f2)} and a band gap from ${freq(flo)} to ${freq(fhi)}. Between ${freq(f2)} and ${freq(fhi)} the effective mass is negative: ${F(me[0] * 1e3, 3)} g at ${freq(fe)}, i.e. a dynamic density of ${F(me[0] / a ** 3, 3)} kg/m³ against ${F(rhoS, 3)} kg/m³ statically.`,
      method: 'Closed-form two-degree-of-freedom Bloch dispersion; complex effective mass with 2 % loss factor', dof: 2 };
  }
  const fl = FLUIDS[v.fluid], R = min(v.memR, 0.5 * a), T = v.memT, rs = v.memRhoS, f1 = 2.405 / (2 * PI * R) * sqrt(T / rs), rp = fl.rho + rs / a, fc = sqrt(8 * T / (rp * R * R * a)) / (2 * PI), fe = 0.7 * fc;
  if (R < v.memR) warnings.push(`Membrane radius limited to half the cell size (${len(R)}).`);
  const fs = N.linspace(0.2 * fc, 2 * fc, 240), rho = fs.map(f => rp * (1 - fc * fc / (f * f))), re = rp * (1 - fc * fc / (fe * fe)), B0 = fl.rho * fl.c ** 2, att = fs.map(f => (f < fc ? 2 * PI * f * sqrt(-rp * (1 - fc * fc / (f * f)) / B0) * a * 8.686 : 0)), aLam = a * fc / fl.c;
  return { polys: layout.resonator({ type: 'nanodisc', params: { diameter: 2 * R, segments: 40 } }), pitch: max(a, 2.4 * R),
    kpis: [
      kpi('fc', 'Cut-off frequency f_c', fc, 'Hz', `√(8T/(ρ′R²a))/2π with the average density ρ′ = ρ₀ + ρ_s/a = ${F(rp, 4)} kg/m³ = ${freq(fc)}. Below it the effective density is negative and sound is evanescent.`),
      kpi('f1', 'First membrane resonance', f1, 'Hz', `(2.405/2πR)√(T/ρ_s) = ${freq(f1)} for the clamped circular membrane in vacuum; the lumped stiffness model holds below this frequency.`),
      kpi('rho', 'Effective density at 0.7 f_c', re, 'kg/m³', `ρ_eff = ρ′(1 − f_c²/f²) at ${freq(fe)} = ${F(re, 4)} kg/m³ (host ${fl.name}: ${F(fl.rho, 4)} kg/m³). The membrane tension acts as a restoring spring on the fluid column, giving a Drude-like negative density.`, re < 0),
      kpi('band_hi', 'Upper edge of the negative-density band', fc, 'Hz', `The negative-density (stop) band extends from zero frequency to f_c = ${freq(fc)}.`),
      kpi('att', 'Attenuation per cell at 0.7 f_c', 2 * PI * fe * sqrt(max(0, -re) / B0) * a * 8.686, 'dB', `Im(k)·a × 8.686 with k = ω√(ρ_eff/B₀) = ${F(2 * PI * fe * sqrt(max(0, -re) / B0) * a * 8.686, 4)} dB per membrane.`),
      kpi('km', 'Membrane stiffness (uniform pressure)', 8 * PI * T, 'N/m', `Force per mean displacement of a clamped circular membrane under uniform pressure, 8πT = ${F(8 * PI * T, 4)} N/m.`),
      kpi('aLam', 'Cell size / wavelength at f_c', aLam, '–', `a·f_c/c = ${F(aLam, 3)}; the description as an effective fluid needs this well below 0.25.`, aLam < 0.25),
    ],
    series: [
      { key: 'rho', title: 'Effective density vs frequency', xlabel: 'frequency (Hz)', ylabel: 'ρ_eff / ρ′ (–)', lines: [{ name: 'ρ_eff/ρ′', x: fs, y: rho.map(q => max(-25, q / rp)) }], explain: 'Drude-type dynamic density of a tube loaded with tensioned membranes: negative below the cut-off, rising towards the average density at high frequency.' },
      { key: 'att', title: 'Attenuation per cell', xlabel: 'frequency (Hz)', ylabel: 'attenuation (dB per cell)', lines: [{ name: 'Im(k)·a', x: fs, y: att }], explain: 'Evanescent decay per membrane below the cut-off frequency; above it the wave propagates and the lossless model gives no attenuation.' },
    ],
    eqs: [eqn('constitutive', 'Effective density of a membrane-loaded tube', 'ρ_eff(ω) = ρ′(1 − ω_c²/ω²),  ρ′ = ρ₀ + ρ_s/a,  ω_c² = k_m/(ρ′ A a),  k_m = 8πT'), eqn('formula', 'Clamped circular membrane', 'f₁ = (2.405/2πR)√(T/ρ_s);  w(r) = p(R² − r²)/(4T)'), eqn('governing', 'Wave propagation', 'k = ω√(ρ_eff/B₀): evanescent for ρ_eff < 0')],
    checks: [chk('Negative density below the cut-off', re < 0, `ρ_eff = ${F(re, 4)} kg/m³ at ${freq(fe)}`), chk('Cut-off below the first membrane resonance (lumped stiffness valid)', fc < 1.5 * f1, `f_c = ${freq(fc)}, f₁ = ${freq(f1)}`), chk('Sub-wavelength cell', aLam < 0.25, `a/λ = ${F(aLam, 3)} at f_c`), chk('High-frequency limit recovers the average density', abs(rho[rho.length - 1] / rp - 0.75) < 1e-6, `ρ_eff(2f_c) = ${F(rho[rho.length - 1] / rp, 4)} ρ′ = 1 − 1/4`)],
    recs: [
      rec('high', `Use the structure as a low-frequency barrier below ${freq(fc)}; raise the cut-off with higher tension or a smaller radius (f_c ∝ √T/R).`, 'Drude-type effective density.'),
      rec('medium', `Tension relaxes by creep in polymer films: a 20 % loss of tension lowers f_c by about 10 %. Use a low-creep film (polyimide, PET) and a rigid frame, and re-measure f_c after thermal cycling.`, 'f_c ∝ √T.'),
      rec('medium', 'Adding a small central mass creates a second, tunable antiresonance with a narrow total-reflection peak; model it with the full membrane modes before relying on it.', 'Limit of the single-stiffness lumped model.'),
      rec('low', `A ${F(rs * 1e3, 3)} g/m² membrane blocks low-frequency sound that would otherwise need a wall obeying the mass law at many kg/m², a large saving in material and transported mass.`, 'Negative-density barrier against the mass law; resource efficiency.', true),
    ],
    pub: { 'geom.cell_size': a, 'res.f0': fc, 'res.bandgap_lo': 0, 'res.bandgap_hi': fc, 'eff.rho_dyn': re },
    summary: `A tube in ${fl.name} loaded with clamped membranes (radius ${len(R)}, tension ${F(T, 3)} N/m, ${F(rs * 1e3, 3)} g/m²) every ${len(a)} behaves as a fluid with Drude-type density: negative below the cut-off ${freq(fc)} (first membrane resonance ${freq(f1)}). At ${freq(fe)} the effective density is ${F(re, 3)} kg/m³ and each cell attenuates by ${F(2 * PI * fe * sqrt(max(0, -re) / B0) * a * 8.686, 3)} dB.`,
    method: 'Lumped membrane stiffness and fluid inertia (closed form)', dof: 1 };
}
function runResonant(v) {
  const warnings = [], em = ['srr', 'wire', 'srrwire'].includes(v.resType), r = em ? resEM(v, warnings) : resAcoustic(v, warnings);
  const big = r.pitch >= 1e-4, s = big ? 1e6 : 1e9, unit = big ? 'µm' : 'nm';
  const lay = layout.array(r.polys.map(q => ({ layer: q.layer, pts: q.pts.map(p => [p[0] * s, p[1] * s]) })), { nx: 4, ny: 4, px: r.pitch * s, py: r.pitch * s }), st = layout.layoutStats(lay);
  const feat = min(st.minFeature ?? Infinity, st.minGap ?? Infinity), featM = Number.isFinite(feat) ? feat / s : r.pitch, pr = printable(featM, v);
  if (lay.capped) warnings.push(lay.capped + '.');
  r.checks.push(chk(`Smallest feature or gap resolvable by ${pr.p.name}`, pr.ok, `minimum feature ${st.minFeature === null ? 'n/a' : len(st.minFeature / s)}, minimum gap ${st.minGap === null ? 'n/a' : len(st.minGap / s)} against ${len(pr.p.res)}; fill factor ${F(st.fillFactor * 100, 3)} %`));
  r.kpis.push(kpi('feature', 'Smallest layout feature / process resolution', pr.margin, '–', `Smallest line width, slot or clearance in the 4 × 4 array (${len(featM)}) divided by the typical minimum feature of ${pr.p.name} (${len(pr.p.res)}) = ${F(pr.margin, 3)}.`, pr.ok));
  if (!pr.ok) warnings.push(`The smallest feature (${len(featM)}) is below the resolution of ${pr.p.name} (${len(pr.p.res)}).`);
  return {
    summary: r.summary, kpis: r.kpis, series: r.series,
    tables: [{ key: 'layout', title: 'Layout statistics (4 × 4 array)', columns: ['Quantity', 'Value'], rows: [['Polygons', st.nPolygons], ['Period', len(r.pitch)], ['Minimum feature', st.minFeature === null ? 'n/a' : len(st.minFeature / s)], ['Minimum gap', st.minGap === null ? 'n/a' : len(st.minGap / s)], ['Fill factor (%)', +(st.fillFactor * 100).toPrecision(4)], ['Process', pr.p.name]], explain: 'Geometric statistics of the planar pattern from the layout engine: smallest line or slot within a polygon, smallest clearance between polygons, and the patterned area fraction.' }],
    equations: r.eqs, geometry: { kind: 'layout', polygons: lay.polygons, units: unit }, checks: r.checks, warnings, recommendations: r.recs, publish: r.pub,
    exports: [{ name: 'resonator_array.json', mime: 'application/json', content: JSON.stringify({ units: unit, period: [r.pitch * s, r.pitch * s], cells: [4, 4], unit_cell: r.polys.map(q => ({ layer: q.layer, pts: q.pts.map(p => [+(p[0] * s).toPrecision(7), +(p[1] * s).toPrecision(7)]) })) }, null, 1), description: 'Unit-cell polygons and array definition for layout tools.' },
      { name: 'response.csv', mime: 'text/csv', content: [r.series[0].xlabel + ',' + r.series[0].lines.map(l => l.name).join(','), ...Array.from(r.series[0].lines[0].x, (x, i) => [x, ...r.series[0].lines.map(l => l.y[i])].map(q => (+q).toPrecision(6)).join(','))].join('\n') + '\n', description: 'Frequency response of the effective parameter shown in the first chart.' }],
    numerics: { method: r.method, dof: r.dof, steps: 1, iterations: 1, residual: 0, notes: 'Lumped-element effective-medium estimate; confirm with full-wave simulation and parameter retrieval.' },
  };
}

// ---------------------------------------------------------------- periodic-structure band design
function blochStack(c1, c2, zr, d1, d2) {
  const tB = d1 / c1 + d2 / c2, fB = 1 / (2 * tB), rhs = f => cos(2 * PI * f * d1 / c1) * cos(2 * PI * f * d2 / c2) - 0.5 * (zr + 1 / zr) * sin(2 * PI * f * d1 / c1) * sin(2 * PI * f * d2 / c2);
  return { fB, rhs };
}
function firstGap(st, n = 500) {
  const fmax = 1.9 * st.fB; let lo = null, hi = null, prev = 0, pv = st.rhs(0);
  for (let i = 1; i <= n; i++) { const f = fmax * i / n, val = st.rhs(f), out = abs(val) > 1, pout = abs(pv) > 1;
    if (out && !pout && lo === null) lo = N.bisect(x => abs(st.rhs(x)) - 1, prev, f, 1e-9 * st.fB);
    if (!out && pout && lo !== null) { hi = N.bisect(x => 1 - abs(st.rhs(x)), prev, f, 1e-9 * st.fB); break; }
    prev = f; pv = val; }
  return lo !== null && hi !== null ? [lo, hi] : null;
}
function runBand(v) {
  const warnings = [], a = v.cell, opt = v.waveType === 'optical', d1 = v.fill * a, d2 = a - d1;
  const c1 = opt ? N.C0 / v.nA : sqrt(v.Es / v.rhos), c2 = opt ? N.C0 / v.nB : sqrt(v.E2 / v.rho2), Z1 = opt ? v.nA : v.rhos * c1, Z2 = opt ? v.nB : v.rho2 * c2, zr = Z1 / Z2;
  const st = blochStack(c1, c2, zr, d1, d2), gap = firstGap(st), fB = st.fB, lo = gap ? gap[0] : fB, hi = gap ? gap[1] : fB, mid = (lo + hi) / 2, ratio = (hi - lo) / mid;
  if (!gap) warnings.push('No band gap was resolved: the impedance contrast is too small (or the fill fraction too extreme). The Bragg frequency is reported for both edges.');
  const fqw = c1 / (c1 + c2), refl = abs(Z1 - Z2) / (Z1 + Z2), qwRatio = 4 / PI * Math.asin(refl), stq = blochStack(c1, c2, zr, fqw * a, (1 - fqw) * a), gq = firstGap(stq), qwNum = gq ? (gq[1] - gq[0]) / ((gq[1] + gq[0]) / 2) : 0;
  // dispersion and attenuation
  const fmax = max(3.2 * fB, 1.2 * hi), fs = N.linspace(0, fmax, 400), R = fs.map(st.rhs), bx = [], by = [], att = R.map(r => (abs(r) > 1 ? Math.acosh(min(1e300, abs(r))) * 8.686 : 0));
  fs.forEach((f, i) => { if (abs(R[i]) <= 1) { bx.push(Math.acos(R[i]) / PI); by.push(f); } });
  const kap = gap ? Math.acosh(abs(st.rhs(mid))) : 0, n40 = kap > 0 ? Math.ceil(40 / (8.686 * kap)) : 0;
  // sweeps
  const fills = N.linspace(0.05, 0.95, 37), gf = fills.map(ff => { const g = firstGap(blochStack(c1, c2, zr, ff * a, (1 - ff) * a), 260); return g ? (g[1] - g[0]) / ((g[1] + g[0]) / 2) : 0; });
  const zrs = N.logspace(1.05, 30, 40), gz = zrs.map(z => 4 / PI * Math.asin((z - 1) / (z + 1)));
  const aReq = a * mid / v.fTarget, ceff = a / (d1 / c1 + d2 / c2), [u, s] = unitOf(a);
  const bragg = [['1-D stack / simple cubic (100)', 1], ['Square (11) / simple cubic (110)', 1 / Math.SQRT2], ['Hexagonal (10)', sqrt(3) / 2], ['Face-centred cubic (111)', 1 / sqrt(3)], ['Body-centred cubic (110)', 1 / Math.SQRT2]].map(([n, k]) => [n, +(k * a * s).toPrecision(4), +(ceff / (2 * k * a)).toPrecision(5)]);
  // voxel picture: four periods along x
  const ppc = 12, nx = 4 * ppc, data = new Uint8Array(nx * 8 * 8), nA = max(1, min(ppc - 1, round(v.fill * ppc))); for (let k = 0; k < 8; k++) for (let j = 0; j < 8; j++) for (let i = 0; i < nx; i++) data[i + nx * (j + 8 * k)] = i % ppc < nA ? 1 : 2;
  const fA = opt ? chem.prettyFormula('TiO2') : null, fBm = opt ? chem.prettyFormula('SiO2') : null;
  const Ereuss = 1 / (v.fill / v.Es + (1 - v.fill) / v.E2), rhoAv = v.fill * v.rhos + (1 - v.fill) * v.rho2, nEff = sqrt(v.fill * v.nA ** 2 + (1 - v.fill) * v.nB ** 2);
  const lamMid = (opt ? N.C0 : ceff) / mid;
  return {
    summary: `The ${opt ? 'optical' : 'elastic'} two-phase stack of period ${len(a)} (fill ${F(v.fill, 3)}, ${opt ? `indices ${F(v.nA, 3)} / ${F(v.nB, 3)}` : `wave speeds ${F(c1, 4)} / ${F(c2, 4)} m/s, impedance ratio ${F(zr, 3)}`}) has its first band gap from ${freq(lo)} to ${freq(hi)}, centred at ${freq(mid)} with a gap-to-mid-gap ratio of ${F(ratio * 100, 3)} %. The mid-gap attenuation is ${F(kap * 8.686, 3)} dB per cell. To centre the gap at ${freq(v.fTarget)} the period must be ${len(aReq)}.`,
    kpis: [
      kpi('lo', 'Band-gap lower edge', lo, 'Hz', `Lower edge of the first Bloch stop band from |cos qa| = 1 = ${freq(lo)}. Below it waves propagate with a real Bloch wavenumber.`),
      kpi('hi', 'Band-gap upper edge', hi, 'Hz', `Upper edge of the first stop band = ${freq(hi)}.`),
      kpi('mid', 'Mid-gap frequency', mid, 'Hz', `Centre of the gap = ${freq(mid)}; the Bragg estimate 1/(2(d₁/c₁ + d₂/c₂)) is ${freq(fB)}. The ${opt ? 'vacuum' : 'effective'} wavelength there is ${len(lamMid)}.`),
      kpi('ratio', 'Gap-to-mid-gap ratio', ratio, '–', `(f_hi − f_lo)/f_mid = ${F(ratio * 100, 4)} %. The quarter-wave stack of the same materials (fill ${F(fqw, 3)}) gives ${F(qwNum * 100, 4)} % numerically and (4/π)asin|Z₁ − Z₂|/(Z₁ + Z₂) = ${F(qwRatio * 100, 4)} % analytically.`, ratio > 0.1),
      kpi('att', 'Mid-gap attenuation', kap * 8.686, 'dB/cell', `Im(qa) = acosh|cos qa| at mid-gap = ${F(kap, 4)} Np per cell (amplitude); about ${n40 || 'many'} cells give 40 dB.`),
      kpi('aReq', 'Period for the target frequency', aReq, 'm', `The dispersion relation is scale invariant, so a_required = a·f_mid/f_target = ${len(aReq)} for a mid-gap at ${freq(v.fTarget)} (layers ${len(aReq * v.fill)} and ${len(aReq * (1 - v.fill))}).`),
      kpi('zr', opt ? 'Index contrast n_A/n_B' : 'Impedance contrast Z_A/Z_B', zr, '–', `${opt ? 'Refractive-index' : 'Acoustic-impedance (ρc)'} ratio of the two phases = ${F(zr, 4)}; the interface amplitude reflection coefficient is ${F(refl, 4)}. The gap width is governed by this contrast.`),
      kpi('fqw', 'Quarter-wave fill fraction', fqw, '–', `Fill fraction of phase A giving equal transit time (or optical thickness) in both layers, c₁/(c₁ + c₂) = ${F(fqw, 4)}; it maximises the first gap.`),
    ],
    series: [
      { key: 'bands', title: 'Bloch dispersion (pass bands)', xlabel: 'reduced wavenumber qa/π (–)', ylabel: 'frequency (Hz)', lines: [{ name: 'propagating bands', x: bx.length ? bx : [0], y: by.length ? by : [0] }], explain: 'Real Bloch wavenumber against frequency in the reduced zone. Frequency ranges with no points are band gaps, where the wave decays from cell to cell.' },
      { key: 'att', title: 'Attenuation in the gaps', xlabel: 'frequency (Hz)', ylabel: 'attenuation (dB per cell)', lines: [{ name: 'Im(qa)', x: fs, y: att }], explain: 'Evanescent decay per unit cell inside the stop bands, 8.686·acosh|cos qa|; zero in the pass bands of the lossless stack.' },
      { key: 'fill', title: 'Gap-to-mid-gap ratio vs fill fraction', xlabel: 'fill fraction of phase A (–)', ylabel: 'gap/mid-gap (–)', lines: [{ name: 'first gap', x: fills, y: gf }, { name: 'this design', x: [v.fill], y: [ratio] }], explain: 'Relative width of the first gap as the layer thickness ratio changes at fixed materials. The maximum lies at the quarter-wave condition.' },
      { key: 'contrast', title: 'Quarter-wave gap vs contrast', xlabel: opt ? 'index ratio (–)' : 'impedance ratio (–)', ylabel: 'gap/mid-gap (–)', logx: true, lines: [{ name: '(4/π)·asin((Z−1)/(Z+1))', x: zrs, y: gz }, { name: 'this material pair (quarter-wave)', x: [max(zr, 1 / zr)], y: [qwNum] }], explain: 'Widest achievable first gap for a given contrast. High-contrast pairs open wide gaps with few cells; low contrast needs many periods.' },
    ],
    tables: [{ key: 'bragg', title: 'Bragg condition for two- and three-dimensional lattices of the same constant', columns: ['Lattice plane family', `Plane spacing d (${u})`, 'Bragg frequency c_eff/(2d) (Hz)'], rows: bragg, explain: 'First-order Bragg frequencies for normal incidence on the densest plane families using the transit-time average speed of this material pair. They locate where directional gaps open in 2-D and 3-D periodic structures.' }],
    equations: [
      eqn('governing', 'Bloch dispersion of a two-layer period', `cos(qa) = cos(k₁d₁)cos(k₂d₂) − ½(Z₁/Z₂ + Z₂/Z₁) sin(k₁d₁) sin(k₂d₂),  kᵢ = ω/cᵢ,  ${opt ? 'Zᵢ ∝ nᵢ, cᵢ = c/nᵢ' : 'Zᵢ = ρᵢcᵢ, cᵢ = √(Eᵢ/ρᵢ)'}`),
      eqn('formula', 'Band gap and attenuation', '|cos qa| > 1 ⇒ stop band;  Im(qa) = acosh|cos qa|'),
      eqn('formula', 'Quarter-wave stack', 'd₁/c₁ = d₂/c₂;  Δf/f_mid = (4/π) asin(|Z₁ − Z₂|/(Z₁ + Z₂))'),
      eqn('formula', 'Bragg condition and scaling', 'f_B = c_eff/(2d);  a_required = a·f_mid/f_target'),
      ...(opt ? [eqn('formula', 'Example material pair', `${fA} (n ≈ 2.4) / ${fBm} (n ≈ 1.46) quarter-wave mirror`)] : []),
    ],
    geometry: { kind: 'voxels', nx, ny: 8, nz: 8, data, voxel: a * s / ppc, units: u },
    checks: [
      chk('Gap edges bracket the Bragg frequency', !gap || (lo <= fB * 1.02 && hi >= fB * 0.98), `f_lo = ${freq(lo)}, f_B = ${freq(fB)}, f_hi = ${freq(hi)}`),
      chk('Quarter-wave gap matches the analytical formula', abs(qwNum - qwRatio) < 0.01 + 0.02 * qwRatio, `numerical ${F(qwNum, 4)} against analytical ${F(qwRatio, 4)}`),
      chk('Long-wave limit: Bloch speed equals the effective-medium speed', (() => { const f = 1e-3 * fB, q = Math.acos(max(-1, min(1, st.rhs(f)))) / a, cb = 2 * PI * f / q, ce = opt ? N.C0 / nEff : sqrt(Ereuss / rhoAv); return abs(cb / ce - 1) < 0.01; })(), opt ? `effective index √(f n_A² + (1 − f) n_B²) = ${F(nEff, 4)}` : `Reuss modulus ${F(Ereuss, 4)} Pa and mean density ${F(rhoAv, 4)} kg/m³`),
      chk('A band gap exists', !!gap, gap ? `relative width ${F(ratio * 100, 3)} %` : 'contrast too low'),
    ],
    warnings,
    recommendations: [
      rec('high', `Set the period to ${len(aReq)} (layers ${len(aReq * v.fill)} / ${len(aReq * (1 - v.fill))}) to centre the gap at ${freq(v.fTarget)}; check that both layers are within the thickness control of the deposition or bonding process.`, 'Scale invariance of the Bloch dispersion.'),
      rec('medium', abs(v.fill - fqw) > 0.08 ? `Move the fill fraction from ${F(v.fill, 3)} towards the quarter-wave value ${F(fqw, 3)}: the gap would widen from ${F(ratio * 100, 3)} % to ${F(qwNum * 100, 3)} %.` : `The fill fraction is close to the quarter-wave optimum (${F(fqw, 3)}); further widening needs a higher contrast pair.`, 'Gap width against fill fraction.'),
      rec('medium', `Use at least ${n40 || 'several tens of'} periods for 40 dB rejection at mid-gap; fewer periods leave a leaky gap and finite-size ripples at the band edges.`, 'Evanescent decay per cell at mid-gap.'),
      rec('low', `Because attenuation is geometric (${F(kap * 8.686, 3)} dB per cell without absorption), the stack needs no lossy or heavy filler: choose the lowest-impact pair with sufficient contrast and avoid more periods than the rejection target requires.`, 'Interference-based rejection; material efficiency.', true),
    ],
    publish: { 'geom.cell_size': a, 'res.bandgap_lo': lo, 'res.bandgap_hi': hi, 'res.f0': mid, 'res.gap_ratio': ratio, 'geom.cell_size_required': aReq, ...(opt ? { 'eff.n': nEff, 'eff.eps': nEff * nEff } : { 'eff.E': Ereuss, 'eff.density': rhoAv }) },
    exports: [{ name: 'dispersion.csv', mime: 'text/csv', content: 'frequency_Hz,cos_qa,attenuation_dB_per_cell\n' + fs.map((f, i) => `${f.toPrecision(7)},${max(-1e6, min(1e6, R[i])).toPrecision(6)},${min(1e6, att[i]).toPrecision(5)}`).join('\n') + '\n', description: 'Bloch dispersion function and attenuation against frequency.' },
      { name: 'stack_design.json', mime: 'application/json', content: JSON.stringify({ wave: v.waveType, period_m: a, layers_m: [d1, d2], speeds_m_s: [c1, c2], contrast: zr, gap_Hz: [lo, hi], midgap_Hz: mid, gap_to_midgap: ratio, period_for_target_m: aReq }, null, 1), description: 'Layer thicknesses and band-gap summary.' }],
    numerics: { method: 'Transfer-matrix Bloch condition for a two-layer period; frequency scan with bisection refinement of the band edges', dof: 2, steps: 400 + 37 * 260, iterations: 60, residual: abs(qwNum - qwRatio), notes: 'Normal incidence, lossless and non-dispersive phases; one polarisation.' },
  };
}

// ---------------------------------------------------------------- bistable, origami and kirigami cells
function runBistable(v) {
  const warnings = [], type = v.bisType;
  if (type === 'beam') {
    const l = v.cell, t = min(v.strut, 0.2 * l), h = v.bh, b = v.depth, Q = h / t, I = b * t ** 3 / 12, p4 = PI ** 4, sc = v.Es * I * h / l ** 3;
    const F1 = d => 1.5 * p4 * Q * Q * d * ((d - 1.5) ** 2 - (0.25 - 4 / (3 * Q * Q))), F3 = d => 8 * p4 - 6 * p4 * d;
    const ds = N.linspace(0, 2, 401), g = ds.map(d => F1(d) - F3(d)), cross = []; for (let i = 1; i < g.length; i++) if (g[i - 1] * g[i] < 0) cross.push(i);
    const da = cross.length >= 2 ? ds[cross[0]] : null, db = cross.length >= 2 ? ds[cross[cross.length - 1]] : null, Fd = d => (da !== null && d >= da && d <= db ? F3(d) : F1(d)), Fs = ds.map(Fd);
    const iTop = N.argmax(Fs.slice(0, 201)), fTop = Fs[iTop], rest = Fs.slice(iTop), iBot = iTop + N.argmin(rest), fBot = Fs[iBot], bistable = fBot < 0;
    let n1 = null, n2 = null; for (let i = 1; i < Fs.length; i++) if (Fs[i] < Fs[i - 1]) { if (n1 === null) n1 = ds[i - 1]; n2 = ds[i]; }
    const kNeg = n1 !== null ? (Fd((n1 + n2) / 2 + 0.01) - Fd((n1 + n2) / 2 - 0.01)) / 0.02 * v.Es * I / l ** 3 : 0;
    let dEnd = 2; if (bistable) for (let i = Fs.length - 1; i > iBot; i--) if (Fs[i - 1] < 0 && Fs[i] >= 0) { dEnd = ds[i]; break; }
    let dMid = null; if (bistable) for (let i = iTop; i <= iBot; i++) if (Fs[i] <= 0) { dMid = ds[i]; break; }
    const Ein = N.trapz(ds.filter(d => d <= (dMid ?? 2)), Fs.filter((f, i) => ds[i] <= (dMid ?? 2))) * sc * h, Eback = bistable ? -N.trapz(ds.filter(d => d >= dMid && d <= dEnd), Fs.filter((f, i) => ds[i] >= dMid && ds[i] <= dEnd)) * sc * h : 0;
    const strain = 2 * PI * PI * t * h / (l * l) + 4 * PI * PI * t * t / (3 * l * l), allow = v.ys / v.Es;
    if (t < v.strut) warnings.push(`Beam thickness limited to 20 % of the span (${len(t)}).`);
    if (Q > 6) warnings.push('h/t above about 6: the beam snaps through higher modes and the two-mode force curve is approximate.');
    if (strain > allow) warnings.push(`The peak strain (${F(strain * 100, 3)} %) exceeds the elastic limit σ_ys/E_s = ${F(allow * 100, 3)} %: the beam yields before snapping. Increase the span-to-thickness ratio.`);
    const [u, s] = unitOf(l), nodes = [], edges = [], seg = 24, gapY = max(3 * t, 0.15 * h);
    for (let c = 0; c < 3; c++) for (const off of [0, gapY]) { const b0 = nodes.length; for (let i = 0; i <= seg; i++) { const x = i / seg; nodes.push([(c * 1.15 + x) * l * s, (h / 2 * (1 - cos(2 * PI * x)) + off) * s, 0]); if (i) edges.push([b0 + i - 1, b0 + i]); } if (off) edges.push([b0 - seg - 1 + seg / 2, b0 + seg / 2]); }
    return {
      summary: `The cosine-shaped beam (span ${len(l)}, thickness ${len(t)}, rise ${len(h)}, depth ${len(b)}) has h/t = ${F(Q, 3)} and is ${bistable ? 'bistable' : 'monostable (below the criterion h/t ≥ 2.31)'}. It snaps through at ${F(fTop * sc, 3)} N${bistable ? ` and back at ${F(fBot * sc, 3)} N, resting in a second stable position ${len(dEnd * h)} from the first` : ''}. The tangent stiffness is negative between ${F(n1 ?? 0, 3)} and ${F(n2 ?? 0, 3)} of the rise (${F(kNeg, 3)} N/m at the centre of that range). Peak strain ${F(strain * 100, 3)} % against an elastic limit of ${F(allow * 100, 3)} %.`,
      kpis: [
        kpi('Q', 'Geometry constant h/t', Q, '–', `Rise over thickness = ${F(Q, 4)}. A centrally clamped double beam is bistable for h/t ≥ 2.31 (Qiu–Lang–Slocum criterion); larger values deepen the second well.`, Q >= 2.31),
        kpi('ftop', 'Snap-through force', fTop * sc, 'N', `Largest force on the forward path, ${F(fTop, 4)}·E I h/l³ = ${F(fTop * sc, 4)} N (for large h/t it tends to about 740 E I h/l³).`),
        kpi('fbot', 'Snap-back (minimum) force', fBot * sc, 'N', `Minimum force on the path, ${F(fBot, 4)}·E I h/l³ = ${F(fBot * sc, 4)} N. A negative value means the beam holds its second position without external force.`, fBot < 0),
        kpi('kneg', 'Negative tangent stiffness', kNeg, 'N/m', `dF/dd in the centre of the unstable range = ${F(kNeg, 4)} N/m (−6π⁴EI/l³ when the second mode is suppressed). Under displacement control this branch is stable and can cancel a positive spring.`, kNeg < 0),
        kpi('range', 'Negative-stiffness travel', ((n2 ?? 0) - (n1 ?? 0)) * h, 'm', `The force falls with displacement from ${len((n1 ?? 0) * h)} to ${len((n2 ?? 0) * h)} of centre travel, i.e. over ${F(((n2 ?? 0) - (n1 ?? 0)) * 100, 3)} % of the rise.`),
        kpi('dend', 'Second stable position', dEnd * h, 'm', bistable ? `Centre displacement of the second zero-force state = ${F(dEnd, 4)} h = ${len(dEnd * h)}.` : `The beam is monostable; the full stroke 2h = ${len(2 * h)} is reported.`),
        kpi('barrier', 'Energy barrier (forward)', Ein, 'J', `Work needed to reach the unstable equilibrium, ∫F dd = ${F(Ein, 4)} J; the barrier for the return snap is ${F(Eback, 4)} J.`),
        kpi('strain', 'Peak strain', strain, '–', `ε_max ≈ 2π²th/l² + 4π²t²/(3l²) = ${F(strain * 100, 4)} % (bending plus axial compression during snap-through) against σ_ys/E_s = ${F(allow * 100, 4)} %.`, strain <= allow),
      ],
      series: [
        { key: 'fd', title: 'Force–displacement curve', xlabel: 'centre displacement d/h (–)', ylabel: 'force (N)', lines: [{ name: 'force', x: ds, y: Fs.map(f => f * sc) }, { name: 'first-mode solution', x: ds, y: ds.map(d => F1(d) * sc) }], explain: 'Transverse force at mid-span against centre displacement. The falling part is the negative-stiffness branch; a crossing of zero force with positive slope after it is the second stable state.' },
        { key: 'k', title: 'Tangent stiffness', xlabel: 'centre displacement d/h (–)', ylabel: 'stiffness dF/dd (N/m)', lines: [{ name: 'tangent stiffness', x: ds.slice(1), y: ds.slice(1).map((d, i) => (Fs[i + 1] - Fs[i]) / (ds[i + 1] - ds[i]) * v.Es * I / l ** 3) }], explain: 'Slope of the force–displacement curve. Negative values mark the range usable for quasi-zero-stiffness isolators or for energy trapping.' },
      ],
      tables: [], equations: [
        eqn('formula', 'Beam shape and normalisation', 'w(x) = (h/2)(1 − cos(2πx/l));  Q = h/t;  F = f l³/(E I h);  Δ = d/h'),
        eqn('constitutive', 'First-mode force (Qiu–Lang–Slocum)', 'F₁ = (3π⁴Q²/2)·Δ·(Δ − 3/2 + √(1/4 − 4/(3Q²)))·(Δ − 3/2 − √(1/4 − 4/(3Q²)))'),
        eqn('constitutive', 'Second mode suppressed (double beam clamped at the centre)', 'F₃ = 8π⁴ − 6π⁴Δ between the intersections with F₁;  k = −6π⁴EI/l³'),
        eqn('formula', 'Bistability criterion and strain', 'bistable ⇔ Q ≥ 2.31;  ε_max ≈ 2π²th/l² + 4π²t²/(3l²)'),
      ],
      geometry: { kind: 'lattice', nodes, edges, radius: t * s / 2, units: u },
      checks: [chk('Bistability criterion h/t ≥ 2.31', Q >= 2.31, `h/t = ${F(Q, 4)}; minimum force ${F(fBot, 4)} E I h/l³`), chk('Elastic strain limit respected', strain <= allow, `ε_max = ${F(strain * 100, 3)} % against ${F(allow * 100, 3)} %`), chk('Shallow-beam assumption (h ≪ l, t ≪ l)', h / l < 0.15 && t / l < 0.05, `h/l = ${F(h / l, 3)}, t/l = ${F(t / l, 3)}`), chk('Energy consistency: forward barrier exceeds the return barrier', !bistable || Ein >= Eback, `forward ${F(Ein, 3)} J, return ${F(Eback, 3)} J (the as-fabricated state is the deeper well)`)],
      warnings,
      recommendations: [
        rec('high', bistable ? `Bistable with margin (h/t = ${F(Q, 3)}); hold fabrication so that h/t stays above about 2.5 — thickening by ${F((1 - 2.31 / Q) * 100, 3)} % would remove the second state.` : `Raise the rise or thin the beam until h/t ≥ 2.5 (now ${F(Q, 3)}) to obtain a second stable state; below 2.31 the beam only softens.`, 'Qiu–Lang–Slocum criterion.'),
        rec('high', strain > allow ? `Reduce the strain below ${F(allow * 100, 3)} %: lengthen the span to at least ${len(l * sqrt(strain / allow))} at the same t and h, or use a material with higher σ_ys/E_s.` : `The strain margin is ${F(allow / strain, 3)}; fatigue-critical designs should keep it above about 2.`, 'ε_max ≈ 2π²th/l².'),
        rec('medium', `Pair the negative branch (${F(kNeg, 3)} N/m) with a positive spring of equal magnitude for a quasi-zero-stiffness isolator, and always clamp two parallel beams at the centre so that the asymmetric mode cannot twist the shuttle.`, 'Negative tangent stiffness on the F₃ branch.'),
        rec('low', 'Snap-through cells store and dissipate energy without plastic deformation, so impact absorbers built from them are reusable rather than single-use crushable foams, reducing waste.', 'Elastic, recoverable energy trapping.', true),
      ],
      publish: { 'geom.cell_size': l, 'geom.strut': t, 'mech.buckling_load': fTop * sc, 'mech.snap_back_force': fBot * sc, 'eff.k_negative': kNeg, 'mech.max_strain': strain, 'mech.safety_factor': allow / strain },
      exports: [{ name: 'force_displacement.csv', mime: 'text/csv', content: 'displacement_m,force_N\n' + ds.map((d, i) => `${(d * h).toPrecision(6)},${(Fs[i] * sc).toPrecision(6)}`).join('\n') + '\n', description: 'Force–displacement curve of the curved beam.' }],
      numerics: { method: 'Closed-form two-mode solution of the pre-curved clamped beam (Qiu–Lang–Slocum), evaluated on 401 displacement points', dof: 2, steps: 401, iterations: 1, residual: 0, notes: 'Shallow, slender, linear-elastic beam; second (asymmetric) buckling mode suppressed.' },
    };
  }
  if (type === 'miura') {
    const A = v.miA, B = v.miB, gm = rad(v.miGamma), th = rad(v.miTheta), tg = tan(gm);
    const dims = t => { const k = sqrt(1 + cos(t) ** 2 * tg * tg); return { H: A * sin(t) * sin(gm), S: B * cos(t) * tg / k, L: A * sqrt(1 - sin(t) ** 2 * sin(gm) ** 2), V: B / k }; };
    const d = dims(th), nuSL = -(cos(th) ** 2) * tg * tg, e = 1e-5, dp = dims(th + e), dm = dims(th - e), nuNum = -(log(dp.L) - log(dm.L)) / (log(dp.S) - log(dm.S)), nuSH = -(log(dp.H) - log(dm.H)) / (log(dp.S) - log(dm.S));
    const ts = N.linspace(2, 88, 44), sheetT = min(v.strut, 0.05 * min(A, B)), rhoRel = min(1, A * B * sin(gm) * sheetT / (max(1e-300, d.S * d.L * d.H))), flat = { S: B * sin(gm), L: A };
    const [u, s] = unitOf(max(A, B)), nodes = [], edges = [], nI = 5, nJ = 5, id = (i, j) => i * nJ + j;
    for (let i = 0; i < nI; i++) for (let j = 0; j < nJ; j++) nodes.push([i * d.S * s, (j * d.L + (i % 2 ? d.V : 0)) * s, (j % 2 ? d.H : 0) * s]);
    for (let i = 0; i < nI; i++) for (let j = 0; j < nJ; j++) { if (i + 1 < nI) edges.push([id(i, j), id(i + 1, j)]); if (j + 1 < nJ) edges.push([id(i, j), id(i, j + 1)]); }
    return {
      summary: `The Miura-ori cell (a = ${len(A)}, b = ${len(B)}, γ = ${F(v.miGamma, 3)}°) folded to ${F(v.miTheta, 3)}° has a unit of ${len(2 * d.S)} × ${len(2 * d.L)} × ${len(d.H)}. Its in-plane Poisson's ratio is ν_SL = −cos²θ·tan²γ = ${F(nuSL, 4)}: the sheet contracts in both in-plane directions as it folds. The out-of-plane ratio ν_SH is ${F(nuSH, 4)} (positive), and the folded core has a relative density of ${F(rhoRel, 3)} for a ${len(sheetT)} sheet.`,
      kpis: [
        kpi('nuSL', "In-plane Poisson's ratio ν_SL", nuSL, '–', `−(dL/L)/(dS/S) = −cos²θ·tan²γ = ${F(nuSL, 4)}; the finite-difference value from the fold kinematics is ${F(nuNum, 4)}. It is always negative and purely geometric (rigid facets, one degree of freedom).`, nuSL < 0),
        kpi('nuLS', "Conjugate ratio ν_LS", 1 / nuSL, '–', `1/ν_SL = ${F(1 / nuSL, 4)}: for a one-degree-of-freedom mechanism the two in-plane ratios are reciprocal.`),
        kpi('nuSH', "Out-of-plane Poisson's ratio ν_SH", nuSH, '–', `−(dH/H)/(dS/S) = ${F(nuSH, 4)}: the height grows as the sheet is compressed along S, a positive ratio of large magnitude near the flat state.`),
        kpi('S', 'Unit-cell width 2S', 2 * d.S, 'm', `2b·cosθ·tanγ/√(1 + cos²θ·tan²γ) = ${len(2 * d.S)}; ${F(d.S / flat.S * 100, 3)} % of the flat width.`),
        kpi('L', 'Unit-cell length 2L', 2 * d.L, 'm', `2a√(1 − sin²θ·sin²γ) = ${len(2 * d.L)}; ${F(d.L / flat.L * 100, 3)} % of the flat length.`),
        kpi('H', 'Folded height H', d.H, 'm', `a·sinθ·sinγ = ${len(d.H)}.`),
        kpi('area', 'Projected area / flat area', d.S * d.L / (flat.S * flat.L), '–', `Deployment ratio of the folded sheet = ${F(d.S * d.L / (flat.S * flat.L), 4)}; its inverse is the stowage ratio.`),
        kpi('rho', 'Relative density of the folded core', rhoRel, '–', `Sheet volume over the enclosing volume, a·b·sinγ·t/(S·L·H) = ${F(rhoRel, 4)} for a sheet thickness of ${len(sheetT)}.`),
      ],
      series: [
        { key: 'nu', title: "In-plane Poisson's ratio vs fold angle", xlabel: 'fold angle θ (°)', ylabel: 'ν_SL (–)', lines: [{ name: `γ = ${F(v.miGamma, 3)}°`, x: ts, y: ts.map(q => -(cos(rad(q)) ** 2) * tg * tg) }, { name: 'γ = 45°', x: ts, y: ts.map(q => -(cos(rad(q)) ** 2)) }], explain: 'The ratio is most negative near the flat state (−tan²γ) and tends to zero when fully folded. Larger sector angles give stronger auxetic coupling.' },
        { key: 'dims', title: 'Cell dimensions vs fold angle', xlabel: 'fold angle θ (°)', ylabel: `dimension (${u})`, lines: [{ name: 'S', x: ts, y: ts.map(q => dims(rad(q)).S * s) }, { name: 'L', x: ts, y: ts.map(q => dims(rad(q)).L * s) }, { name: 'H', x: ts, y: ts.map(q => dims(rad(q)).H * s) }], explain: 'Both in-plane dimensions shrink together as the pattern folds while the height grows: the signature of a negative in-plane and positive out-of-plane Poisson ratio.' },
      ],
      tables: [], equations: [
        eqn('formula', 'Miura-ori fold kinematics (rigid facets)', 'H = a sinθ sinγ;  S = b cosθ tanγ/√(1 + cos²θ tan²γ);  L = a√(1 − sin²θ sin²γ);  V = b/√(1 + cos²θ tan²γ)'),
        eqn('constitutive', "In-plane Poisson's ratio", 'ν_SL = −(dL/L)/(dS/S) = −cos²θ·tan²γ < 0;  ν_LS = 1/ν_SL'),
        eqn('formula', 'Out-of-plane ratio', 'ν_SH = −(dH/H)/(dS/S) = cot²θ·(1 + cos²θ tan²γ) > 0'),
      ],
      geometry: { kind: 'lattice', nodes, edges, radius: max(sheetT * s / 2, 0.004 * max(A, B) * s), units: u },
      checks: [chk('Closed form matches the fold kinematics', abs(nuNum - nuSL) < 1e-4 * max(1, abs(nuSL)), `finite difference ${F(nuNum, 6)} against −cos²θ tan²γ = ${F(nuSL, 6)}`), chk('Facet edge lengths preserved (rigid folding)', abs(hypot(d.S, d.V) / B - 1) < 1e-9 && abs(hypot(d.L, d.H) / A - 1) < 1e-9, `√(S² + V²)/b = ${F(hypot(d.S, d.V) / B, 6)}, √(L² + H²)/a = ${F(hypot(d.L, d.H) / A, 6)}`), chk('Negative in-plane Poisson ratio', nuSL < 0, `ν_SL = ${F(nuSL, 4)}`), chk('Thin-sheet assumption', sheetT < 0.05 * min(A, B) + 1e-15, `sheet thickness ${len(sheetT)} against facet side ${len(min(A, B))}`)],
      warnings,
      recommendations: [
        rec('high', `Choose γ for the coupling needed: ν_SL ranges from −tan²γ = ${F(-tg * tg, 3)} (flat) to 0 (folded). At the current fold it is ${F(nuSL, 3)}; a sector angle of 45° would give ${F(-(cos(th) ** 2), 3)}.`, 'ν_SL = −cos²θ tan²γ.'),
        rec('medium', 'Real sheets have finite thickness and crease stiffness: use thickness-accommodating hinge offsets or perforated creases, and expect facet bending to add a second, softer deformation mode.', 'Rigid-facet, zero-thickness kinematic idealisation.'),
        rec('medium', `As a sandwich core at θ = ${F(v.miTheta, 3)}° the relative density is ${F(rhoRel, 3)}; open channels along the straight creases allow ventilation and drainage, unlike closed honeycombs.`, 'Fold geometry and core density.'),
        rec('low', `Folding from flat sheet stows to ${F(d.S * d.L / (flat.S * flat.L) * 100, 3)} % of the deployed area and produces no cutting scrap, reducing transport volume and material waste.`, 'Deployment ratio; zero-offcut forming.', true),
      ],
      publish: { 'eff.nu': nuSL, 'geom.cell_size': 2 * d.S, 'geom.thickness': d.H, 'geom.rel_density': rhoRel, 'eff.nu_out_of_plane': nuSH },
      exports: [{ name: 'miura_vertices.csv', mime: 'text/csv', content: `x_${u},y_${u},z_${u}\n` + nodes.map(p => p.map(q => +q.toPrecision(6)).join(',')).join('\n') + '\n', description: 'Vertex coordinates of the folded 4 × 4 facet patch.' }],
      numerics: { method: 'Closed-form rigid-origami kinematics; central finite differences for verification', dof: 1, steps: 44, iterations: 1, residual: abs(nuNum - nuSL), notes: 'One-degree-of-freedom rigid folding; facets and creases carry no elastic energy.' },
    };
  }
  if (type === 'kirigami') {
    const Lc = v.kiLc, d = v.kiD, gH = min(v.kiG, 0.8 * Lc), t = v.strut, Lb = (Lc - gH) / 2, P = Lc + gH, emax = sqrt(1 + (Lb / d) ** 2) - 1, rot = deg(Math.atan(Lb / d));
    const Ein = 2 * (d / Lb) ** 3 * d / P, Eout = Ein * min(1, (t / d) ** 2), ecrit = min(emax, (t / Lb) ** 2 * PI * PI / 12), cutFrac = Lc / P, strainMat = 3 * t / (2 * Lb) * (rot * PI / 180);
    if (gH < v.kiG) warnings.push(`Hinge width limited to 80 % of the cut length (${len(gH)}).`);
    if (t > d) warnings.push('The sheet is thicker than the ligament width: the ligaments bend in-plane rather than buckling out of plane, and the high-stretch regime is not reached.');
    const big = P >= 1e-4, s = big ? 1e6 : 1e9, unit = big ? 'µm' : 'nm', polys = [], sw = 0.06 * d;
    for (let r = 0; r < 8; r++) for (let c = -1; c < 4; c++) { const x0 = (c + (r % 2 ? 0.5 : 0)) * P + gH / 2, y = (r + 0.5) * d; const xa = max(0, x0), xb = min(4 * P, x0 + Lc); if (xb - xa > 0.02 * Lc) polys.push({ layer: 1, pts: [[xa * s, (y - sw / 2) * s], [xb * s, (y - sw / 2) * s], [xb * s, (y + sw / 2) * s], [xa * s, (y + sw / 2) * s]] }); }
    const ratios = N.logspace(0.5, 30, 40);
    return {
      summary: `The staggered-slit kirigami sheet (cuts ${len(Lc)}, row spacing ${len(d)}, hinges ${len(gH)}, thickness ${len(t)}) can extend geometrically by up to ${F(emax * 100, 3)} % as the ${len(Lb)} ligaments rotate by ${F(rot, 3)}°. Its initial in-plane modulus is about ${F(Ein, 2)} E_s; after the ligaments buckle out of plane (beyond a strain of about ${F(ecrit * 100, 2)} %) the tangent modulus falls to roughly ${F(Eout, 2)} E_s.`,
      kpis: [
        kpi('emax', 'Geometric stretchability', emax, '–', `Rigid-ligament limit √(1 + (L_b/d)²) − 1 = ${F(emax * 100, 4)} %: the strain at which the ligaments are fully aligned with the pull. Beyond it the material itself must stretch.`),
        kpi('Lb', 'Ligament (beam) length', Lb, 'm', `(L_c − hinge)/2 = ${len(Lb)}: the free length of each strip between staggered hinges.`),
        kpi('Ein', 'Initial in-plane modulus / E_s', Ein, '–', `Fixed–guided beam estimate 2(d/L_b)³·d/P = ${F(Ein, 4)}: before buckling the ligaments bend in the plane of the sheet.`),
        kpi('Eout', 'Post-buckling modulus / E_s', Eout, '–', `In-plane value times (t/d)² = ${F(Eout, 4)}: once the ligaments twist and bend out of plane the sheet is softer by the ratio of the bending stiffnesses.`),
        kpi('ecrit', 'Buckling onset strain (estimate)', ecrit, '–', `Order-of-magnitude strain for out-of-plane buckling, (π²/12)(t/L_b)² = ${F(ecrit * 100, 4)} %; it marks the end of the stiff initial regime.`),
        kpi('rot', 'Ligament rotation at full extension', rot, '°', `atan(L_b/d) = ${F(rot, 4)}°.`),
        kpi('strain', 'Hinge bending strain at full extension', strainMat, '–', `Fibre strain of a ligament bent through the rotation over its length, about 3tφ/(2L_b) = ${F(strainMat * 100, 4)} %, against σ_ys/E_s = ${F(v.ys / v.Es * 100, 3)} %. Stress concentrates at the slit tips and is higher there.`, strainMat < v.ys / v.Es),
        kpi('cut', 'Cut fraction along a row', cutFrac, '–', `L_c/(L_c + hinge) = ${F(cutFrac, 4)}: the share of each row that is slit.`),
      ],
      series: [
        { key: 'emax', title: 'Stretchability vs slit geometry', xlabel: 'ligament length / row spacing L_b/d (–)', ylabel: 'geometric strain limit (–)', logx: true, lines: [{ name: '√(1 + (L_b/d)²) − 1', x: ratios, y: ratios.map(q => sqrt(1 + q * q) - 1) }, { name: 'this design', x: [Lb / d], y: [emax] }], explain: 'The extension limit depends only on the ratio of ligament length to row spacing: long slits closely spaced give very large stretch from an inextensible sheet.' },
        { key: 'ss', title: 'Idealised stress–strain response', xlabel: 'strain (–)', ylabel: 'stress / E_s (–)', lines: [{ name: 'bilinear estimate', x: [0, ecrit, emax], y: [0, Ein * ecrit, Ein * ecrit + Eout * (emax - ecrit)] }], explain: 'Stiff in-plane bending up to the buckling strain, then a long soft plateau while ligaments rotate out of plane, ending at the geometric limit where the response stiffens sharply.' },
      ],
      tables: [], equations: [
        eqn('formula', 'Geometric extension limit (rigid ligaments)', 'L_b = (L_c − g)/2;  ε_max = √(1 + (L_b/d)²) − 1;  φ_max = atan(L_b/d)'),
        eqn('constitutive', 'Beam estimates of the two regimes', 'E_in/E_s ≈ 2(d/L_b)³ d/(L_c + g);  E_out ≈ E_in (t/d)²;  ε_c ~ (π²/12)(t/L_b)²'),
      ],
      geometry: { kind: 'layout', polygons: polys, units: unit },
      checks: [chk('Out-of-plane regime available (t < d)', t < d, `t/d = ${F(t / d, 3)}`), chk('Hinge strain below the elastic limit', strainMat < v.ys / v.Es, `${F(strainMat * 100, 3)} % against ${F(v.ys / v.Es * 100, 3)} %`), chk('Slender ligaments (L_b ≥ 3d) for the beam model', Lb >= 3 * d, `L_b/d = ${F(Lb / d, 3)}`), chk('Extension limit positive', emax > 0, `ε_max = ${F(emax * 100, 3)} %`)],
      warnings,
      recommendations: [
        rec('high', `For more stretch raise L_b/d (now ${F(Lb / d, 3)}): doubling it would give about ${F((sqrt(1 + (2 * Lb / d) ** 2) - 1) * 100, 3)} %. Keep hinges at least as wide as the row spacing so they do not tear.`, 'Geometric extension limit.'),
        rec('high', strainMat >= v.ys / v.Es ? `Hinge strain (${F(strainMat * 100, 3)} %) exceeds the elastic limit: thin the sheet or lengthen the ligaments, and round the slit ends to cut the stress concentration.` : 'Round or key-hole the slit ends: the nominal hinge strain is acceptable but slit tips concentrate stress and initiate fatigue cracks.', 'Bending strain of rotating ligaments.'),
        rec('medium', `Operate in the soft plateau between ${F(ecrit * 100, 2)} % and about ${F(emax * 80, 3)} % strain for stretchable electronics; beyond that the ligaments carry membrane tension and conductors on them will crack.`, 'Bilinear regime estimate.'),
        rec('low', 'Slitting removes no material: kirigami obtains large compliance from a stiff, thin sheet with zero offcut, so thin recyclable mono-material films can replace elastomers.', 'Cut-only patterning; material substitution.', true),
      ],
      publish: { 'geom.cell_size': P, 'geom.strut': t, 'eff.E': Ein * v.Es, 'mech.max_strain': emax, 'eff.E_post_buckling': Eout * v.Es },
      exports: [{ name: 'kirigami_cuts.json', mime: 'application/json', content: JSON.stringify({ units: unit, cut_length: Lc * s, row_spacing: d * s, hinge_width: gH * s, slits: polys.map(p => p.pts) }, null, 1), description: 'Slit rectangles of the staggered cut pattern.' }],
      numerics: { method: 'Closed-form geometric and beam estimates', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Order-of-magnitude model: ligaments as fixed–guided beams, rigid hinges, no slit-tip mechanics.' },
    };
  }
  // rotating squares: rigid-unit kinematics
  const a = v.cell, phi = rad(min(44, v.theta)), sq = a / (2 * (cos(phi) + sin(phi))), w = min(v.strut, 0.2 * sq), Kh = v.Es * w * w / 12, Ehinge = 8 * Kh / (sq * sq * max(1e-6, 1 - sin(2 * phi))), rho = 1 / (1 + sin(2 * phi));
  const ps = N.linspace(0, 45, 46), X = p => 2 * sq * (cos(rad(p)) + sin(rad(p))), e = 1e-4, nuNum = -((log(X(deg(phi) + e)) - log(X(deg(phi) - e))) / (log(X(deg(phi) + e)) - log(X(deg(phi) - e)))), remain = sqrt(2) / (cos(phi) + sin(phi)) - 1;
  const lg = latticeGeometry('rotatingSquares2D', { theta_deg: deg(phi) }, a, min(w, 0.05 * a), 3);
  return {
    summary: `Rigid squares of side ${len(sq)} rotated by ${F(deg(phi), 3)}° tile a ${len(a)} cell and deform by counter-rotation only, so both in-plane dimensions change by the same factor: ν = −1 at every angle. The solid fraction is ${F(rho, 3)}; ${F(remain * 100, 3)} % of further isotropic expansion remains before the fully open state at 45°. With ${len(w)} flexure hinges the hinge-limited modulus is about ${F(Ehinge, 3)} Pa.`,
    kpis: [
      kpi('nu', "Poisson's ratio", -1, '–', `Both cell dimensions equal 2s(cosφ + sinφ), so the strains are equal and ν = −1 exactly for rigid squares with ideal hinges (finite-difference check: ${F(nuNum, 4)}). This is the isotropic lower bound in two dimensions.`, true),
      kpi('s', 'Square side', sq, 'm', `a/(2(cosφ + sinφ)) = ${len(sq)} for the ${len(a)} cell at φ = ${F(deg(phi), 3)}°.`),
      kpi('rho', 'Solid area fraction', rho, '–', `4s²/a² = 1/(1 + sin 2φ) = ${F(rho, 4)}; it falls from 1 (closed) to 0.5 (fully open at 45°).`),
      kpi('remain', 'Remaining linear expansion', remain, '–', `√2/(cosφ + sinφ) − 1 = ${F(remain * 100, 4)} % before the squares lock in the fully open configuration.`),
      kpi('areal', 'Maximum areal expansion from closed', 1, '–', 'The cell area doubles between the closed (φ = 0) and fully open (φ = 45°) states: (cosφ + sinφ)² = 1 + sin 2φ ranges from 1 to 2.'),
      kpi('E', 'Hinge-limited modulus', Ehinge, 'Pa', `8K_h/(s²(1 − sin 2φ)) per unit depth with a flexure-hinge stiffness K_h = E_s w²/12 (hinge length equal to its width ${len(w)}) = ${F(Ehinge, 4)} Pa. It rises without bound as the fully open state is approached.`),
      kpi('Kh', 'Hinge rotational stiffness per unit depth', Kh, 'N', `E_s w³/(12 l_h) with l_h = w: ${F(Kh, 4)} N·m per metre depth per radian.`),
    ],
    series: [
      { key: 'X', title: 'Cell size vs square rotation', xlabel: 'rotation of each square φ (°)', ylabel: 'cell size / closed size (–)', lines: [{ name: 'X/X₀ = Y/Y₀', x: ps, y: ps.map(p => cos(rad(p)) + sin(rad(p))) }], explain: 'Both in-plane dimensions follow the same curve, which is why the Poisson ratio is −1. The structure expands by 41 % in each direction between closed and fully open.' },
      { key: 'E', title: 'Hinge-limited modulus vs rotation', xlabel: 'rotation of each square φ (°)', ylabel: 'E (Pa)', logy: true, lines: [{ name: 'E', x: ps.slice(0, 45), y: ps.slice(0, 45).map(p => 8 * Kh / (sq * sq * max(1e-6, 1 - sin(2 * rad(p))))) }], explain: 'With elastic hinges the stiffness depends on the opening angle and diverges at the locked, fully open state (square side held fixed at the current value).' },
    ],
    tables: [], equations: [
      eqn('formula', 'Rotating-squares kinematics', 'X = Y = 2s(cos φ + sin φ);  ε_x = ε_y  ⇒  ν = −1;  A/A₀ = 1 + sin 2φ'),
      eqn('constitutive', 'Hinge-limited stiffness (Grima–Evans type)', 'E = 8K_h / (s²(1 − sin 2φ)) per unit depth;  K_h = E_s w³/(12 l_h)'),
    ],
    geometry: lg.geo,
    checks: [chk("Kinematic Poisson's ratio equals −1", abs(nuNum + 1) < 1e-9, `finite difference of ln X against ln Y: ${F(nuNum, 6)}`), chk('Thermodynamic admissibility in 2-D (−1 ≤ ν ≤ 1)', true, 'ν = −1 is the limiting isotropic value: the shear modulus is finite while the areal modulus is hinge-limited'), chk('Slender hinges (w ≤ 0.2 s)', w <= 0.2 * sq + 1e-15, `w/s = ${F(w / sq, 3)}`)],
    warnings: v.strut > 0.2 * sq ? [...warnings, `Hinge width limited to 20 % of the square side (${len(w)}).`] : warnings,
    recommendations: [
      rec('high', `Use thin flexure hinges (w/s = ${F(w / sq, 3)}) to stay close to ν = −1; thick hinges make the squares deform and move the ratio towards zero, as the frame analysis of the welded cell shows.`, 'Rigid-unit assumption.'),
      rec('medium', `Fabricate near φ = ${F(deg(phi), 3)}° with ${F(remain * 100, 3)} % expansion remaining; operating close to 45° gives a rapidly stiffening, locking response useful as an end stop.`, 'X ∝ cos φ + sin φ.'),
      rec('medium', 'Hinges see the full rotation strain: check the hinge fibre strain w·Δφ/(2l_h) against the elastic limit for the working stroke and consider living hinges in a fatigue-resistant polymer.', 'Flexure-hinge strain.'),
      rec('low', 'Perforated-sheet versions are cut from a single sheet with small slits, giving auxetic behaviour without assembly and with negligible scrap.', 'Cut-pattern realisation of rotating units; low waste.', true),
    ],
    publish: { 'eff.nu': -1, 'eff.E': Ehinge, 'geom.cell_size': a, 'geom.strut': w, 'geom.rel_density': rho },
    exports: latticeExports(lg.geo, { type: 'rotatingSquares2D', cell_size_m: a, square_side_m: sq, rotation_deg: deg(phi), nu: -1 }),
    numerics: { method: 'Closed-form rigid-unit kinematics; finite-difference verification', dof: 1, steps: 46, iterations: 1, residual: abs(nuNum + 1), notes: 'Rigid squares with ideal or flexure hinges.' },
  };
}

// ---------------------------------------------------------------- negative thermal expansion / negative compressibility
function runNegCTE(v) {
  const warnings = [], a = v.cell;
  if (v.nteType === 'triangle') {
    const th = rad(v.nteTheta), ab = v.alphaB, as = v.alphaS, cte = q => (as - ab * cos(q) ** 2) / sin(q) ** 2, ah = cte(th), ratio = as / ab, th0 = ab !== 0 && ratio > 0 && ratio < 1 ? deg(Math.acos(sqrt(ratio))) : null;
    const b = a, Ls = b / (2 * cos(th)), h = b / 2 * tan(th), t = min(v.strut, 0.2 * b), comp = Ls / (2 * v.Es * t * sin(th) ** 2) + b / (4 * v.E2 * t * tan(th) ** 2), kh = 1 / comp, Eh = kh * h / b;
    const dT = 100, dh = ah * h * dT, sens = deg(1) * 0; const dAdth = (cte(th + rad(0.5)) - cte(th - rad(0.5)));   // per degree
    const ths = N.linspace(10, 80, 71), uc = { nodes: [[0, 0, 0], [0.5, 0, 0], [1, 0, 0], [0.5, tan(th) / 2, 0]], edges: [[0, 1], [1, 2], [0, 3], [2, 3]], size: [1, tan(th) / 2, 0] }, [u, s] = unitOf(a), geo = tessellate(uc, 2, v.nView, a, t, s); geo.units = u;
    if (t < v.strut) warnings.push(`Member thickness limited to 20 % of the base length (${len(t)}).`);
    if (ab <= as) warnings.push('The base member does not expand more than the inclined members: the height expansion cannot be reduced below the constituent values. Put the high-expansion material in the base.');
    const al = chem.prettyFormula('Al'), ti = chem.prettyFormula('Ti');
    return {
      summary: `A pin-jointed triangle with a base of CTE ${F(ab * 1e6, 3)} ppm/K and inclined members of ${F(as * 1e6, 3)} ppm/K at a base angle of ${F(v.nteTheta, 3)}° has a height expansion coefficient of ${F(ah * 1e6, 3)} ppm/K${ah < 0 ? ' (negative: the apex moves down on heating)' : ''}. ${th0 !== null ? `Zero expansion is obtained at a base angle of ${F(th0, 3)}°; shallower triangles are negative.` : 'No zero-expansion angle exists for this pair.'} The stiffness of one triangle against an apex load is ${F(kh, 3)} N/m per unit depth.`,
      kpis: [
        kpi('cte', 'Effective CTE along the height', ah, '1/K', `α_h = (α_s − α_b cos²θ)/sin²θ = ${F(ah * 1e6, 4)} ppm/K. The base expansion lowers the apex geometrically; no thermal stress arises because the pin-jointed triangle is statically determinate.`, abs(ah) < abs(as)),
        kpi('th0', 'Zero-expansion base angle', th0 === null ? 'none' : th0, th0 === null ? '' : '°', th0 === null ? 'cos²θ₀ = α_s/α_b has no solution between 0 and 1 for these coefficients, so the height CTE cannot be tuned to zero.' : `cos²θ₀ = α_s/α_b gives θ₀ = ${F(th0, 4)}°; below it α_h is negative, above it positive and approaching α_s.`),
        kpi('ctex', 'CTE along the base', ab, '1/K', `Along the base the cell simply follows the base member: ${F(ab * 1e6, 4)} ppm/K. A lattice with zero or negative expansion in all in-plane directions needs a three-fold arrangement in which the bases are internal.`),
        kpi('careal', 'Areal CTE of the triangle', ab + ah, '1/K', `Sum of the base and height expansion coefficients, α_b + α_h = ${F((ab + ah) * 1e6, 4)} ppm/K; it gives the in-plane area change of the triangular cell per kelvin, and is zero for an area-preserving design.`),
        kpi('sens', 'Angle sensitivity of the CTE', dAdth, '1/K per °', `Change of α_h per degree of base angle = ${F(dAdth * 1e6, 4)} ppm/K per degree (central difference): the tolerance on the angle sets how closely zero expansion can be hit.`),
        kpi('kh', 'Apex stiffness per unit depth', kh, 'N/m²', `Inverse of L_s/(2E_s t sin²θ) + b/(4E_b t tan²θ) = ${F(kh, 4)} N/m per metre depth; the equivalent modulus along the height is ${F(Eh, 4)} Pa. Shallow, strongly negative designs are compliant.`),
        kpi('dh', 'Height change for ΔT = 100 K', dh, 'm', `Change in triangle height for a 100 K temperature rise, α_h·h·ΔT with h = ${len(h)}: ${F(dh * 1e6, 4)} µm. A negative value means the cell contracts along its height on heating.`),
      ],
      series: [
        { key: 'cte', title: 'Height CTE vs base angle', xlabel: 'base angle θ (°)', ylabel: 'α_h (ppm/K)', lines: [{ name: 'α_h', x: ths, y: ths.map(q => max(-500, min(500, cte(rad(q)) * 1e6))) }, { name: 'inclined members α_s', x: [10, 80], y: [as * 1e6, as * 1e6] }, { name: 'this design', x: [v.nteTheta], y: [max(-500, min(500, ah * 1e6))] }], explain: 'The height expansion (clipped to ±500 ppm/K) crosses zero at the design angle and becomes strongly negative for shallow triangles, at the price of stiffness.' },
        { key: 'k', title: 'Apex stiffness vs base angle', xlabel: 'base angle θ (°)', ylabel: 'stiffness per depth (N/m²)', logy: true, lines: [{ name: 'k_h', x: ths, y: ths.map(q => { const r = rad(q); return 1 / (b / (2 * cos(r)) / (2 * v.Es * t * sin(r) ** 2) + b / (4 * v.E2 * t * tan(r) ** 2)); }) }], explain: 'Stiffness of one triangle against a load at the apex. The trade-off between low or negative thermal expansion and stiffness is intrinsic to this mechanism.' },
      ],
      tables: [], equations: [
        eqn('formula', 'Geometry of the triangle', 'h² = L_s² − (b/2)²;  cos θ = b/(2L_s)'),
        eqn('constitutive', 'Effective thermal expansion along the height', 'α_h = (α_s L_s² − α_b (b/2)²)/h² = (α_s − α_b cos²θ)/sin²θ;  α_h = 0 at cos²θ₀ = α_s/α_b'),
        eqn('formula', 'Apex compliance (pin-jointed truss, unit depth)', '1/k_h = L_s/(2E_s t sin²θ) + b/(4E_b t tan²θ)'),
        eqn('formula', 'Example material pair', `${al} base (α ≈ 23 ppm/K) with ${ti} inclined members (α ≈ 8.6 ppm/K)`),
      ],
      geometry: geo,
      checks: [chk('High-expansion material in the base (α_b > α_s)', ab > as, `α_b = ${F(ab * 1e6, 3)}, α_s = ${F(as * 1e6, 3)} ppm/K`), chk('Thermodynamic consistency: constituents have conventional expansion', true, 'the negative effective value is geometric; the volume of each solid member still increases on heating'), chk('Limits recovered', abs(cte(rad(89.999)) - as) < 1e-6 * max(1e-9, abs(as)) + 1e-12, `α_h → α_s = ${F(as * 1e6, 3)} ppm/K as θ → 90°`), chk('Slender members (pin-joint idealisation)', t / Ls < 0.1, `t/L_s = ${F(t / Ls, 3)}; bonded or welded joints add bending restraint that pulls α_h back towards the constituent values`)],
      warnings,
      recommendations: [
        rec('high', th0 !== null ? `For zero expansion set the base angle to ${F(th0, 3)}° and hold it within ±${F(min(5, abs(1e-6 / (dAdth || 1e-12))), 2)}° for a residual below 1 ppm/K (sensitivity ${F(dAdth * 1e6, 3)} ppm/K per degree).` : 'Choose a pair with α_base > α_inclined > 0 so that a zero-expansion angle exists.', 'α_h(θ) and its slope.'),
        rec('medium', `Use real pin or flexure joints: rigid joints add bending restraint. With t/L_s = ${F(t / Ls, 3)} the idealisation is ${t / Ls < 0.1 ? 'reasonable' : 'poor'}; verify with a bi-material frame analysis including joint stiffness.`, 'Statically determinate truss assumption.'),
        rec('medium', `Check the stiffness penalty: k_h = ${F(kh, 3)} N/m² here. Steeper triangles are stiffer but less negative; pick the shallowest angle that still meets the stiffness requirement.`, 'Apex compliance of the truss.'),
        rec('low', 'A geometric zero-expansion lattice can replace Invar or glass-ceramics, avoiding nickel-rich alloys; keep the two metals separable at end of life (mechanical joints rather than brazing) for recycling.', 'Substitution of critical raw materials; design for disassembly.', true),
      ],
      publish: { 'eff.cte': ah, 'geom.cell_size': a, 'geom.strut': t, 'eff.cte_base': ab, 'eff.E': Eh },
      exports: latticeExports(geo, { type: 'bimaterial_triangle', base_length_m: b, base_angle_deg: v.nteTheta, cte_base_per_K: ab, cte_inclined_per_K: as, cte_height_per_K: ah, zero_cte_angle_deg: th0 }),
      numerics: { method: 'Closed-form kinematics of a pin-jointed bi-material triangle; unit-load method for the compliance', dof: 2, steps: 71, iterations: 1, residual: 0, notes: 'Small temperature change, temperature-independent coefficients, ideal pin joints.' },
    };
  }
  // wine-rack: hinged rhombic cell, negative linear compressibility
  const m = { Es: v.Es, nus: 0.3 }, prm = { theta_deg: v.wrTheta }, c = cellModel('rhombic2D', prm, a, v.strut, m), t = c.t, th = rad(v.wrTheta), l = a / (2 * cos(th));
  const an = (q, tt) => { const ll = a / (2 * cos(q)), r = ll / tt; return { bx: (r * sin(2 * q) - r ** 3 * cos(2 * q) * tan(q)) / v.Es, by: (r * sin(2 * q) + r ** 3 * cos(2 * q) / tan(q)) / v.Es }; };
  const bx = c.S[0][0] + c.S[0][1], by = c.S[1][1] + c.S[1][0], bA = bx + by, A0 = an(th, t), ths = N.linspace(6, 84, 40), fr = N.linspace(10, 80, 8).map(q => { const cc = cellModel('rhombic2D', { theta_deg: q }, a, v.strut, m, true); return [q, cc.S[0][0] + cc.S[0][1], cc.S[1][1] + cc.S[1][0]]; });
  const ev = N.eigSym(c.C.map(r => r.map(x => x / v.Es))).values, lg = latticeGeometry('rhombic2D', prm, a, t, v.nView), neg = bx < 0 ? 'x' : by < 0 ? 'y' : null;
  if (c.capped) warnings.push(`Strut thickness limited to ${len(t)} (45 % of the strut length).`);
  if (!neg) warnings.push('Neither linear compressibility is negative at this angle and slenderness: at exactly 45° the hinging term vanishes, and thick struts add positive axial compressibility.');
  return {
    summary: `The wine-rack cell with ${len(l)} struts at ${F(v.wrTheta, 3)}° to x and ${len(t)} walls has linear compressibilities β_x = ${F(bx, 3)} 1/Pa and β_y = ${F(by, 3)} 1/Pa under in-plane hydrostatic pressure${neg ? `: it lengthens along ${neg} when squeezed (negative linear compressibility)` : ''}. The areal compressibility β_x + β_y = ${F(bA, 3)} 1/Pa stays positive, as stability requires. Homogenised Poisson's ratio ν_xy = ${F(c.nu, 3)}.`,
    kpis: [
      kpi('bx', 'Linear compressibility along x', bx, '1/Pa', `β_x = S₁₁ + S₁₂ from the homogenised frame = ${F(bx, 4)} 1/Pa; the hinging-plus-axial closed form gives ${F(A0.bx, 4)} 1/Pa. Negative for strut angles below 45°: pressure closes the rhombus along y and pushes it out along x.`, bx < 0),
      kpi('by', 'Linear compressibility along y', by, '1/Pa', `β_y = S₂₂ + S₁₂ = ${F(by, 4)} 1/Pa (closed form ${F(A0.by, 4)} 1/Pa).`),
      kpi('bA', 'Areal compressibility', bA, '1/Pa', `β_x + β_y = ${F(bA, 4)} 1/Pa: positive, so the area still shrinks under pressure even though one dimension grows. The areal bulk modulus is ${F(1 / bA, 4)} Pa.`, bA > 0),
      kpi('ratio', 'Negative-to-positive ratio', min(bx, by) / max(bx, by), '–', `Smaller over larger linear compressibility = ${F(min(bx, by) / max(bx, by), 4)}; it approaches −1 for slender struts, where the mechanism is nearly area-preserving.`),
      kpi('nu', "Poisson's ratio ν_xy", c.nu, '–', `From the same compliance matrix: ${F(c.nu, 4)} (conjugate ${F(c.nuo, 4)}). Values above 1 in one direction accompany negative linear compressibility in an anisotropic cell.`),
      kpi('Ex', 'Modulus along x', c.E, 'Pa', `E_x = ${F(c.E, 4)} Pa and E_y = ${F(c.Eo, 4)} Pa: the cell is strongly anisotropic and soft, because it deforms by strut bending at the joints.`),
      kpi('rho', 'Relative density', c.rho1, '–', `ΣtL/A = ${F(c.rho1, 4)} for the ${len(t)} walls.`),
    ],
    series: [
      { key: 'beta', title: 'Linear compressibility vs strut angle', xlabel: 'strut angle to x (°)', ylabel: 'compressibility × E_s (–)', lines: [{ name: 'β_x (closed form)', x: ths, y: ths.map(q => an(rad(q), t).bx * v.Es) }, { name: 'β_y (closed form)', x: ths, y: ths.map(q => an(rad(q), t).by * v.Es) }, { name: 'β_x (frame)', x: fr.map(r => r[0]), y: fr.map(r => r[1] * v.Es) }, { name: 'β_y (frame)', x: fr.map(r => r[0]), y: fr.map(r => r[2] * v.Es) }], explain: 'The two linear compressibilities have opposite signs except near 45°, and swap roles across it. The frame solution includes shear deformation of the struts and follows the closed form closely for slender members.' },
      { key: 'area', title: 'Areal compressibility vs strut angle', xlabel: 'strut angle to x (°)', ylabel: '(β_x + β_y) × E_s (–)', logy: true, lines: [{ name: 'β_A (closed form)', x: ths, y: ths.map(q => { const r = an(rad(q), t); return max(1e-12, (r.bx + r.by) * v.Es); }) }], explain: 'The sum of the linear compressibilities is positive at every angle: negative linear compressibility never implies negative volume (area) compressibility.' },
    ],
    tables: [], equations: [
      eqn('formula', 'Wine-rack geometry', 'X = 2l cos θ;  Y = 2l sin θ;  A = 2l² sin 2θ'),
      eqn('constitutive', 'Linear compressibilities (joint bending + strut axial strain, unit depth)', 'β_x = (1/E_s)[(l/t) sin 2θ − (l/t)³ cos 2θ tan θ];  β_y = (1/E_s)[(l/t) sin 2θ + (l/t)³ cos 2θ cot θ]', 'Rigid joints: each strut bends in double curvature, equivalent hinge stiffness K = 3E_sI/l.'),
      eqn('formula', 'Areal compressibility (always positive)', 'β_A = β_x + β_y = (1/E_s)[2(l/t) sin 2θ + 2(l/t)³ cos² 2θ / sin 2θ] > 0'),
      eqn('governing', 'From the homogenised compliance', 'β_i = −∂ε_i/∂p = Σ_j S_ij  (σ = −p·I)'),
    ],
    geometry: lg.geo,
    checks: [chk('Areal (volumetric) compressibility positive', bA > 0, `β_x + β_y = ${F(bA, 3)} 1/Pa`), chk('Stiffness matrix positive definite', ev[0] > 0, `smallest eigenvalue of C/E_s = ${F(ev[0], 3)}`), chk('Frame solution agrees with the closed form', abs(bx - A0.bx) < 0.1 * max(abs(bx), abs(A0.bx)) + 0.02 * abs(by), `β_x: frame ${F(bx, 3)}, closed form ${F(A0.bx, 3)} 1/Pa`), chk('Negative linear compressibility present', !!neg, neg ? `along ${neg}` : 'none at this angle')],
    warnings,
    recommendations: [
      rec('high', neg ? `The cell expands along ${neg} under pressure by ${F(abs(min(bx, by)) * 1e9, 3)} per GPa; use angles well away from 45° (now ${F(v.wrTheta, 3)}°) for a strong effect, and note the response reverses for the other axis.` : 'Move the strut angle away from 45° or use more slender struts to obtain negative linear compressibility.', 'Sign of cos 2θ in the hinging term.'),
      rec('medium', `Slender struts strengthen the effect as (l/t)³ = ${F((l / t) ** 3, 3)} but make the cell soft (areal modulus ${F(1 / bA, 3)} Pa); set l/t from the pressure range so that strains stay small.`, 'Closed-form compressibilities.'),
      rec('medium', 'Pressure must act on the cell boundary, not inside the pores: seal the faces or use a penetrating-fluid-free configuration, otherwise the open framework sees no net load.', 'Definition of drained hydrostatic loading.'),
      rec('low', 'Negative-linear-compressibility elements make pressure sensors and actuators with large displacement from little material and no rare or toxic piezoelectric ceramics.', 'Geometric amplification; substitution of lead-containing materials.', true),
    ],
    publish: { 'eff.compressibility_x': bx, 'eff.compressibility_y': by, 'eff.K': 1 / bA, 'eff.E': c.E, 'eff.nu': c.nu, 'eff.G': c.G, 'geom.cell_size': a, 'geom.strut': t, 'geom.rel_density': c.rho1 },
    exports: latticeExports(lg.geo, cellDescriptor(c, { beta_x_per_Pa: bx, beta_y_per_Pa: by, strut_angle_deg: v.wrTheta })),
    numerics: { method: 'Periodic Timoshenko beam-frame homogenisation of the rhombic cell; closed-form hinge model for comparison', dof: c.sol.dof, steps: 9, iterations: 1, residual: abs(bx - A0.bx) / max(1e-300, abs(bx)), notes: 'Linear response about the stated angle; pressure applied to the cell boundary (drained framework).' },
  };
}

// ---------------------------------------------------------------- graded, hierarchical and finite arrays
function runGraded(v) {
  const warnings = [], m = { Es: v.Es, nus: v.nus }, type = v.cellType, def = cellDef(type), dim = def.dim, a = v.cell, n = v.nGrad, angleCell = ANGLE_CELLS.includes(type), byAngle = v.gradVar === 'angle' && angleCell;
  if (v.gradVar === 'angle' && !angleCell) warnings.push('This cell has no angle parameter: the strut size was graded instead.');
  const law = s => (v.gradLaw === 'linear' ? s : v.gradLaw === 'power' ? s * s : (exp(2 * s) - 1) / (exp(2) - 1));
  const base = cellModel(type, cellParams(type, v.theta, v.hl), a, v.strut, m), t0 = base.t, ax = dim === 3 ? 'z' : 'x';
  const cells = Array.from({ length: n }, (_, i) => { const w = law(i / (n - 1)), th = byAngle ? v.theta + (v.theta2 - v.theta) * w : v.theta, tt = byAngle ? t0 : t0 * (1 + (v.gradRatio - 1) * w); return { th, tt, c: cellModel(type, cellParams(type, th, v.hl), a, tt, m, true) }; });
  const idx = cells.map((_, i) => i + 1), Er = cells.map(q => q.c.E / v.Es), nus = cells.map(q => q.c.nu), rhos = cells.map(q => q.c.rho1), Eseries = n / N.sum(Er.map(e => 1 / max(1e-300, e))), Epar = N.mean(Er), rhoMean = N.mean(rhos);
  if (cells.some(q => q.c.capped)) warnings.push('Some cells reached the strut-size limit (45 % of the shortest strut); the grading saturates there.');
  // hierarchical (order 2) scaling from the fitted power law of this cell
  const fs = N.logspace(0.4, 1.6, 5).map(f => cellModel(type, cellParams(type, v.theta, v.hl), a, t0 * f, m, true)), pl = N.linreg(fs.map(q => log(q.rho1)), fs.map(q => log(max(1e-300, q.E / v.Es)))), nE = fin(pl.slope, 1), CE = exp(fin(pl.intercept, 0)), mB = dim === 3 ? 2 : 3;
  const kb = base.buckle !== null ? base.buckle / (v.Es * base.rho1 ** mB) : null, rT = v.rd, r1 = sqrt(rT), E1 = min(1, CE * rT ** nE), E2 = min(1, CE * r1 ** nE) ** 2, b1 = kb !== null ? kb * rT ** mB : null, b2 = kb !== null ? kb * min(1, CE * r1 ** nE) * r1 ** mB : null;
  // finite-width strips with free edges
  const sec = base.sec, Einf = dim === 3 ? base.Eo : base.E, widths = [], Efin = [];
  for (let w = 1; w <= v.nAcross; w++) { const red = reduceCell(base.uc, dim, a, w); if (red.reps.length * (dim === 3 ? 6 : 3) > 430) break; const sol = frameSolve(red, sec, dim === 3 ? [0, 2] : [0]), S = invSym(sol.C); widths.push(w); Efin.push(fin(1 / S[0][0] / Einf, 0)); }
  if (widths.length < v.nAcross) warnings.push(`The finite-width study stopped at ${widths.length} cells across to keep the frame model small.`);
  const edgeLast = Efin[Efin.length - 1], n95 = widths.find((w, i) => abs(Efin[i] - 1) < 0.05) ?? null;
  // periodic-boundary compatibility of a 2-cell supercell, from the generated geometry (faces one block period apart: re-entrant cells protrude beyond the period, so the bounding box is not the periodic face)
  const sup = latticeGeometry(type, cellParams(type, v.theta, v.hl), a, t0, 2), G = sup.geo, compat = [];
  for (let axn = 0; axn < dim; axn++) { const lo = N.minOf(G.nodes.map(p => p[axn])), hi = lo + sup.n * base.red.per[axn] * sup.scale, ext = hi - lo, tolr = 1e-4 * ext, keyOf = p => p.map((x, k) => (k === axn ? 0 : round(x / tolr / 10))).join(), A = new Set(G.nodes.filter(p => abs(p[axn] - lo) < tolr).map(keyOf)), B = G.nodes.filter(p => abs(p[axn] - hi) < tolr).map(keyOf); const matched = B.filter(k => A.has(k)).length; compat.push(['xyz'[axn], A.size, B.length, matched, A.size === B.length && matched === A.size ? 'compatible' : 'mismatch']); }
  const compatOK = compat.every(r => r[4] === 'compatible');
  // geometry: graded lattice along x
  const [u, s] = unitOf(a), nx = min(n, 8), gl = OWN[type] ? tessellate(def.make(cellParams(type, v.theta, v.hl)), dim, 3, a, t0, s) : implicit.buildLattice({ type, cells: [nx, 2, dim === 3 ? 2 : 1], cellSize_mm: a * s, strut_mm: t0 * s, params: cellParams(type, v.theta, v.hl), grading: byAngle ? undefined : { axis: 'x', from: 1, to: v.gradRatio } });
  const geo = { kind: 'lattice', nodes: gl.nodes, edges: gl.edges, radius: gl.radius, units: u, ...(gl.radii ? { radii: gl.radii } : {}) };
  if (byAngle) warnings.push('The viewer shows the uniform lattice at the first-cell angle: angle-graded tessellations are not conforming without transition cells, which is itself a design constraint.');
  const pr = printable(N.minOf(cells.map(q => q.c.t)), v), pretty = optLabel('cellType', type);
  return {
    summary: `Across ${n} ${pretty} cells graded ${byAngle ? `in angle from ${F(v.theta, 3)}° to ${F(v.theta2, 3)}°` : `in strut size by a factor ${F(v.gradRatio, 3)}`} (${v.gradLaw}), the modulus E_${ax}/E_s runs from ${F(Er[0], 3)} to ${F(Er[n - 1], 3)} and Poisson's ratio from ${F(nus[0], 3)} to ${F(nus[n - 1], 3)}; loaded along the grading the stack behaves as ${F(Eseries, 3)} E_s (series) and across it as ${F(Epar, 3)} E_s (parallel). A strip ${widths[widths.length - 1]} cells wide with free edges reaches ${F(edgeLast * 100, 3)} % of the infinite-array stiffness. At a total relative density of ${F(rT, 3)}, a two-level hierarchy changes the stiffness by a factor ${F(E2 / E1, 3)}${b1 ? ` and the strut-buckling strength by ${F(b2 / b1, 3)}` : ''}.`,
    kpis: [
      kpi('Eser', 'Series modulus along the grading / E_s', Eseries, '–', `Harmonic mean of the cell moduli (equal stress in every cell) = ${F(Eseries, 4)}; the softest cells dominate.`),
      kpi('Epar', 'Parallel modulus across the grading / E_s', Epar, '–', `Arithmetic mean of the cell moduli (equal strain) = ${F(Epar, 4)}; the ratio parallel/series ${F(Epar / Eseries, 3)} measures how strongly graded the stack is.`),
      kpi('Erange', 'Modulus ratio last/first cell', Er[n - 1] / Er[0], '–', `E_${ax} of the last cell over the first = ${F(Er[n - 1] / Er[0], 4)}, each from its own periodic frame homogenisation (local periodicity assumed).`),
      kpi('nuRange', "Poisson's ratio span", nus[n - 1] - nus[0], '–', `ν changes from ${F(nus[0], 4)} to ${F(nus[n - 1], 4)} across the array; a gradient in ν produces shape-morphing (synclastic or anticlastic curvature) under uniform load.`),
      kpi('rho', 'Mean relative density', rhoMean, '–', `Average of the first-order cell densities = ${F(rhoMean, 4)} (from ${F(N.minOf(rhos), 3)} to ${F(N.maxOf(rhos), 3)}).`),
      kpi('edge', 'Finite-width stiffness ratio', edgeLast, '–', `Uniaxial modulus of a strip ${widths[widths.length - 1]} cells wide with free lateral edges, relative to the infinite array = ${F(edgeLast, 4)}.${n95 ? ` Within 5 % of the bulk value from ${n95} cells across.` : ' Not yet within 5 % of the bulk value over the widths studied.'}`, abs(edgeLast - 1) < 0.1),
      kpi('hierE', 'Order-2 / order-1 stiffness at equal density', E2 / E1, '–', `With E/E_s = ${F(CE, 3)}·ρ̄^${F(nE, 3)} fitted to this cell, a two-level lattice with ρ̄ split equally (√ρ̄ at each level) has E₂/E_s = ${F(E2, 3)} against E₁/E_s = ${F(E1, 3)} at ρ̄ = ${F(rT, 3)}.`),
      kpi('hierB', 'Order-2 / order-1 buckling strength', b1 ? b2 / b1 : 'not applicable', b1 ? '–' : '', b1 ? `Strut Euler buckling scales as ρ̄^${mB} for this cell type; with lattice struts the ratio is C·ρ̄^((n−${mB})/2) = ${F(b2 / b1, 4)}. Hierarchy pays off at low density, where it exceeds 1.` : 'No strut is in compression for the primary load case, so no buckling comparison is made.', b1 ? b2 / b1 > 1 : null),
      kpi('feature', 'Thinnest strut / process resolution', pr.margin, '–', `Smallest strut in the graded array ${len(N.minOf(cells.map(q => q.c.t)))} against ${len(pr.p.res)} for ${pr.p.name}.`, pr.ok),
    ],
    series: [
      { key: 'profile', title: 'Property profile across the graded array', xlabel: 'cell index', ylabel: 'ρ̄, E/E_s, ν (–)', lines: [{ name: `E_${ax}/E_s`, x: idx, y: Er }, { name: 'relative density', x: idx, y: rhos }, { name: 'Poisson ratio', x: idx, y: nus }], explain: 'Each cell is homogenised separately with its own strut size or angle. Grading is valid when properties change slowly from cell to cell, so that each cell still sees nearly periodic surroundings.' },
      { key: 'edge', title: 'Finite-array edge effect', xlabel: 'cells across the width', ylabel: 'E_strip / E_infinite (–)', lines: [{ name: 'free-edge strip (frame solution)', x: widths, y: Efin }, { name: 'infinite array', x: [widths[0], widths[widths.length - 1]], y: [1, 1] }], explain: 'Uniaxial stiffness of strips periodic along the load direction with stress-free lateral edges. Boundary cells are less constrained, so narrow specimens deviate from the bulk value; test coupons should be wide enough to sit on the plateau.' },
      { key: 'hier', title: 'Hierarchical scaling', xlabel: 'total relative density ρ̄ (–)', ylabel: 'E/E_s (–)', logx: true, logy: true, lines: [{ name: 'order 1', x: N.logspace(0.002, 0.5, 30), y: N.logspace(0.002, 0.5, 30).map(r => min(1, CE * r ** nE)) }, { name: 'order 2', x: N.logspace(0.002, 0.5, 30), y: N.logspace(0.002, 0.5, 30).map(r => min(1, CE * sqrt(r) ** nE) ** 2) }], explain: 'Stiffness of a simple lattice and of a lattice whose struts are themselves lattices, at the same total density. Hierarchy costs stiffness for stretch-dominated cells (factor C) but raises buckling strength at low density.' },
    ],
    tables: [
      { key: 'compat', title: 'Supercell periodic-boundary compatibility (2-cell block)', columns: ['Axis', 'Nodes on min face', 'Nodes on max face', 'Matched pairs', 'Status'], rows: compat, explain: 'Nodes on opposite faces of the generated block are paired by their in-face coordinates. Every node must have a partner for periodic boundary conditions (or for joining blocks) to be applied node to node.' },
      { key: 'hier', title: 'Order-1 against order-2 lattice at equal total density', columns: ['Quantity', 'Order 1', 'Order 2'], rows: [['Relative density per level', +rT.toPrecision(4), +r1.toPrecision(4)], ['E/E_s', +E1.toPrecision(4), +E2.toPrecision(4)], ['Buckling strength / E_s', b1 ? +b1.toPrecision(4) : 'n/a', b2 ? +b2.toPrecision(4) : 'n/a']], explain: 'Scaling estimate from the fitted modulus law and the Euler-buckling law of this cell; the second-order lattice treats each strut as a homogenised lattice beam (separation of scales assumed).' },
    ],
    equations: [
      eqn('formula', 'Grading law', v.gradLaw === 'linear' ? 'p(i) = p₁ + (p_N − p₁)·s,  s = (i − 1)/(N − 1)' : v.gradLaw === 'power' ? 'p(i) = p₁ + (p_N − p₁)·s²' : 'p(i) = p₁ + (p_N − p₁)·(e^(2s) − 1)/(e² − 1)'),
      eqn('constitutive', 'Stack averages', 'series (along grading): 1/E = (1/N)Σ 1/Eᵢ;  parallel (across): E = (1/N)Σ Eᵢ'),
      eqn('constitutive', 'Hierarchical lattice (two levels, equal split)', 'ρ̄ = ρ̄₁ρ̄₂;  E₂/E_s = (Cρ̄₁ⁿ)(Cρ̄₂ⁿ) = C²ρ̄ⁿ;  σ_b,2/σ_b,1 = C ρ̄^((n−m)/2)'),
      eqn('governing', 'Finite strip', 'periodic along the load, free lateral edges: E_strip = (K_εε − K_εu K_uu⁻¹ K_uε)/V with ũ unconstrained across the width'),
    ],
    geometry: geo,
    checks: [
      chk('Periodic-boundary compatibility of the supercell', compatOK, compat.map(r => `${r[0]}: ${r[3]}/${r[1]} pairs`).join('; ')),
      chk('Finite-width stiffness converges towards the bulk value', Efin.length < 2 || abs(Efin[Efin.length - 1] - 1) <= abs(Efin[0] - 1) + 1e-6, `E_strip/E_∞ from ${F(Efin[0], 4)} (1 cell) to ${F(edgeLast, 4)} (${widths[widths.length - 1]} cells)`),
      chk('Slow grading (scale separation)', N.maxOf(Er.slice(1).map((e, i) => abs(e / Er[i] - 1))) < 0.5, `largest cell-to-cell modulus change ${F(N.maxOf(Er.slice(1).map((e, i) => abs(e / Er[i] - 1))) * 100, 3)} %`),
      chk(`Thinnest strut printable by ${pr.p.name}`, pr.ok, `${len(N.minOf(cells.map(q => q.c.t)))} against ${len(pr.p.res)}`),
      chk('Series modulus below parallel modulus', Eseries <= Epar * (1 + 1e-9), `${F(Eseries, 4)} ≤ ${F(Epar, 4)} E_s (harmonic ≤ arithmetic mean)`),
    ],
    warnings,
    recommendations: [
      rec('high', `Use at least ${n95 ?? `${v.nAcross + 1} or more`} cells across any load-bearing section or test coupon so that edge effects stay below 5 %; narrower sections must be analysed as finite structures, not with homogenised properties.`, 'Free-edge strip solutions.'),
      rec('medium', `Limit the cell-to-cell property step (now up to ${F(N.maxOf(Er.slice(1).map((e, i) => abs(e / Er[i] - 1))) * 100, 3)} %) by using more cells or a gentler law; abrupt steps concentrate stress at the interface between cells.`, 'Scale-separation requirement of graded homogenisation.'),
      rec('medium', b1 && b2 / b1 > 1 ? `At ρ̄ = ${F(rT, 3)} a two-level hierarchy raises buckling strength by ${F(b2 / b1, 3)} × for a stiffness factor of ${F(E2 / E1, 3)}: worthwhile where buckling governs and the finer struts (about one fifth of the macro strut) remain printable.` : `At ρ̄ = ${F(rT, 3)} hierarchy does not pay for this cell (stiffness × ${F(E2 / E1, 3)}); it becomes attractive only at much lower density where strut buckling governs.`, 'Hierarchical scaling laws.'),
      rec('low', `Grading places material only where needed: this array averages ρ̄ = ${F(rhoMean, 3)} instead of the ${F(N.maxOf(rhos), 3)} of its densest cell applied everywhere, a ${F((1 - rhoMean / N.maxOf(rhos)) * 100, 3)} % saving in feedstock and build energy.`, 'Mean against maximum relative density; material efficiency.', true),
    ],
    publish: { 'geom.cell_size': a, 'geom.strut': t0, 'geom.rel_density': rhoMean, 'eff.E': Eseries * v.Es, 'eff.nu': N.mean(nus), 'eff.G': base.G, 'eff.density': rhoMean * v.rhos, 'eff.edge_factor': edgeLast },
    exports: [...latticeExports(geo, cellDescriptor(base, { grading: { variable: byAngle ? 'angle_deg' : 'strut_m', law: v.gradLaw, values: cells.map(q => (byAngle ? q.th : q.c.t)) }, E_over_Es: Er, nu: nus, relative_density: rhos })), { name: 'graded_profile.csv', mime: 'text/csv', content: 'cell,strut_m,angle_deg,relative_density,E_over_Es,nu\n' + cells.map((q, i) => [i + 1, q.c.t, q.th, rhos[i], Er[i], nus[i]].map(x => +(+x).toPrecision(6)).join(',')).join('\n') + '\n', description: 'Cell-by-cell geometry and homogenised properties of the graded array.' }],
    numerics: { method: 'Cell-wise periodic beam-frame homogenisation; free-edge strip frames for the finite-width study; power-law fit for hierarchical scaling', dof: base.sol.dof, steps: n + widths.length + 6, iterations: 1, residual: 0, notes: 'Local periodicity assumed for graded cells; boundary struts of the strips keep their full section.' },
  };
}

const RUN = { lattice: runLattice, tpms: runTPMS, resonant: runResonant, band: runBand, bistable: runBistable, negcte: runNegCTE, graded: runGraded };

export default {
  id: 5,
  tagline: 'Unit-cell design of mechanical lattices, TPMS, resonant, band-gap, bistable, origami, negative-expansion and negative-compressibility metamaterials with homogenised properties and manufacturability checks.',
  inputs,
  presets: [
    { name: 'Re-entrant auxetic lattice, θ = −30°, h/l = 2', note: 'Gibson–Ashby reference case with ν₁₂ = −1 in the thin-wall limit; Ti-6Al-4V by laser powder-bed fusion.', values: { mode: 'lattice', cellType: 'reentrant2D', theta: 30, hl: 2, cell: 5e-3, strut: 2.5e-4, process: 'lpbf' } },
    { name: 'Regular hexagonal honeycomb (ν = 1 check)', note: 'θ = 30°, h/l = 1: isotropic in-plane with ν = 1 in the flexure model.', values: { mode: 'lattice', cellType: 'honeycomb2D', theta: 30, hl: 1, cell: 5e-3, strut: 1.5e-4, process: 'lpbf' } },
    { name: 'Octet nanolattice by two-photon polymerisation', note: 'Stretch-dominated polymer nanolattice: 10 µm cell, 0.5 µm struts (E/E_s ≈ ρ̄/9).', values: { mode: 'lattice', cellType: 'octet', cell: 1e-5, strut: 5e-7, Es: 3e9, nus: 0.35, ys: 7e7, rhos: 1200, process: 'tpp', nView: 2 } },
    { name: 'SRR + wire negative-index cell (X-band)', note: 'Printed copper split ring with added gap capacitance and a continuous wire: double-negative band near 8 GHz.', values: { mode: 'resonant', resType: 'srrwire', cell: 5e-3, srrSize: 3.5e-3, srrW: 4e-4, srrGap: 3e-4, srrT: 3.5e-5, srrCadd: 5e-14, wireR: 1e-4, epsH: 1, process: 'pcb' } },
    { name: 'Gyroid sheet lattice, ρ̄ = 0.3', note: 'Self-supporting TPMS for powder-bed fusion with 4 mm cells.', values: { mode: 'tpms', tpmsType: 'gyroid', tpmsMode: 'sheet', rd: 0.3, rdTo: 0.3, cell: 4e-3, process: 'lpbf' } },
    { name: 'TiO₂ / SiO₂ quarter-wave Bragg mirror at 550 nm', note: 'Optical one-dimensional photonic crystal: period 151.5 nm, fill 0.378.', values: { mode: 'band', waveType: 'optical', cell: 1.515e-7, fill: 0.378, nA: 2.4, nB: 1.46, fTarget: 5.45e14 } },
    { name: 'Bistable polymer snap-through beam', note: 'Double cosine beam, h/t = 3, span 20 mm.', values: { mode: 'bistable', bisType: 'beam', cell: 2e-2, strut: 2e-4, bh: 6e-4, depth: 5e-3, Es: 2.5e9, ys: 6e7 } },
    { name: 'Al / Ti bi-material zero-expansion triangle', note: 'Base angle chosen near the zero-CTE point of the aluminium–titanium pair.', values: { mode: 'negcte', nteType: 'triangle', alphaB: 2.3e-5, alphaS: 8.6e-6, nteTheta: 52.3, cell: 2e-2, strut: 1e-3 } },
    { name: 'Helmholtz-resonator acoustic cell in air', note: 'Negative effective bulk modulus just above 5 kHz.', values: { mode: 'resonant', resType: 'helmholtz', cell: 5e-3, hrV: 5e-8, hrRn: 5e-4, hrLn: 1e-3, fluid: 'air', process: 'sla' } },
  ],
  run(inp) { const v = clean(inp); return RUN[v.mode](v); },
};
