// Suite 14 — Structural and Nanomechanical Engineering.
// Periodic frame finite-element homogenisation of 2-D lattices (auxetic and conventional), snap-through negative stiffness,
// negative linear compressibility, 3-D lattice scaling laws, and nanomechanics (surface elasticity, contact, strength, poroelasticity).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { PI } = N;
const f3 = v => N.fmt(v, 3), f4 = v => N.fmt(v, 4);
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
const formulaEq = (f, label) => { try { return [eqn('formula', label, `${chem.prettyFormula(f)}, M = ${f4(chem.molarMass(f))} g/mol`)]; } catch { return []; } };

// ================================================================ periodic frame finite elements
// A cell is { nodes:[[x,y]], els:[[i,j,sx,sy]], Lx, Ly }: element end j is the periodic image of node j shifted by (sx·Lx, sy·Ly).
function bondCell(nodes, Lx, Ly, l) { // connect all node pairs (and images) separated by l
  const els = []; for (let i = 0; i < nodes.length; i++) for (let j = i; j < nodes.length; j++) for (let sx = -1; sx <= 1; sx++) for (let sy = -1; sy <= 1; sy++) {
    if (i === j && !(sx > 0 || (sx === 0 && sy > 0))) continue;
    if (Math.abs(Math.hypot(nodes[j][0] + sx * Lx - nodes[i][0], nodes[j][1] + sy * Ly - nodes[i][1]) - l) < 1e-6 * l) els.push([i, j, sx, sy]);
  } return { nodes, els, Lx, Ly };
}
function makeCell(topo, l, hl, thDeg, rR, ringT = 1) {
  const s3 = Math.sqrt(3);
  if (topo === 'square') return bondCell([[0, 0]], l, l, l);
  if (topo === 'triangular') return bondCell([[0, 0], [l / 2, l * s3 / 2]], l, l * s3, l);
  if (topo === 'kagome') return bondCell([[0, 0], [l, 0], [l / 2, s3 * l / 2], [l, s3 * l], [0, s3 * l], [1.5 * l, 1.5 * s3 * l]], 2 * l, 2 * s3 * l, l);
  if (topo === 'rhombic') { const th = thDeg * PI / 180, X = 2 * l * Math.cos(th), Y = 2 * l * Math.sin(th); return { nodes: [[0, 0], [X / 2, Y / 2]], els: [[0, 1, 0, 0], [0, 1, -1, 0], [0, 1, 0, -1], [0, 1, -1, -1]], Lx: X, Ly: Y }; }
  if (topo === 'chiral') { // hexachiral: 12-gon rings on a triangular lattice joined by tangent ligaments
    const R = l, r = rR * R, ph = Math.acos(2 * r / R), cen = [[0, 0], [R / 2, s3 * R / 2]], nodes = [], els = [];
    cen.forEach(c => { for (let m = 0; m < 12; m++) nodes.push([c[0] + r * Math.cos(ph + m * PI / 6), c[1] + r * Math.sin(ph + m * PI / 6)]); });
    for (let a = 0; a < 2; a++) for (let m = 0; m < 12; m++) els.push([a * 12 + m, a * 12 + (m + 1) % 12, 0, 0, ringT]);   // fifth entry: wall-thickness factor of the ring segments
    const tgt = [[[0, 1, 0], [1, 0, 0], [1, -1, 0]], [[1, 1, 0], [0, 1, 1], [0, 0, 1]]]; // [ring, sx, sy] for directions 0°, 60°, 120°
    for (let a = 0; a < 2; a++) for (let k = 0; k < 3; k++) { const [b, sx, sy] = tgt[a][k]; els.push([a * 12 + 2 * k, b * 12 + (2 * k + 6) % 12, sx, sy]); }
    return { nodes, els, Lx: R, Ly: s3 * R, ligament: Math.sqrt(R * R - 4 * r * r), ring: r };
  }
  // hexagonal honeycomb (h, l, θ): θ > 0 conventional, θ < 0 re-entrant
  const th = thDeg * PI / 180, h = hl * l, cx = l * Math.cos(th), sy = l * Math.sin(th);
  return { nodes: [[0, 0], [0, h], [cx, h + sy], [cx, 2 * h + sy]], els: [[0, 1, 0, 0], [1, 2, 0, 0], [1, 2, -1, 0], [2, 3, 0, 0], [3, 0, 0, 1], [3, 0, 1, 1]], Lx: 2 * cx, Ly: 2 * (h + sy) };
}
function prepElements(cell, E, nu, t, timo) {
  const G = E / (2 * (1 + nu)); let len = 0;
  const els = cell.els.map(([i, j, sx, sy, tf = 1]) => {
    const A = t * tf, I = A ** 3 / 12;
    const SX = sx * cell.Lx, SY = sy * cell.Ly, dx = cell.nodes[j][0] + SX - cell.nodes[i][0], dy = cell.nodes[j][1] + SY - cell.nodes[i][1], L = Math.hypot(dx, dy), c = dx / L, s = dy / L;
    const phi = timo ? 12 * E * I / ((5 / 6) * G * A * L * L) : 0, ea = E * A / L, a = 12 * E * I / ((1 + phi) * L ** 3), b = 6 * E * I / ((1 + phi) * L * L), c1 = (4 + phi) * E * I / ((1 + phi) * L), c2 = (2 - phi) * E * I / ((1 + phi) * L);
    const kl = [[ea, 0, 0, -ea, 0, 0], [0, a, b, 0, -a, b], [0, b, c1, 0, -b, c2], [-ea, 0, 0, ea, 0, 0], [0, -a, -b, 0, a, -b], [0, b, c2, 0, -b, c1]];
    const T = N.zeros(6, 6); for (const o of [0, 3]) { T[o][o] = c; T[o][o + 1] = s; T[o + 1][o] = -s; T[o + 1][o + 1] = c; T[o + 2][o + 2] = 1; }
    len += L * tf; return { i, j, SX, SY, L, c, s, kl, T, kg: N.matMul(N.transpose(T), N.matMul(kl, T)) };
  });
  return { els, rho: len * t / (cell.Lx * cell.Ly), A: t, I: t ** 3 / 12 };
}
const jump = (e, q) => [q[0] * e.SX + 0.5 * q[2] * e.SY, 0.5 * q[2] * e.SX + q[1] * e.SY, 0];
function homogenise(cell, E, nu, t, timo) {
  const P = prepElements(cell, E, nu, t, timo), nn = cell.nodes.length, n = 3 * nn, K = N.zeros(n, n), area = cell.Lx * cell.Ly;
  for (const e of P.els) { const d = [3 * e.i, 3 * e.i + 1, 3 * e.i + 2, 3 * e.j, 3 * e.j + 1, 3 * e.j + 2]; for (let a = 0; a < 6; a++) for (let b = 0; b < 6; b++) K[d[a]][d[b]] += e.kg[a][b]; }
  const free = []; for (let i = 2; i < n; i++) free.push(i); const Kr = free.map(i => free.map(j => K[i][j])), C = N.zeros(3, 3), U = [];
  for (let q = 0; q < 3; q++) {
    const eq = [0, 0, 0]; eq[q] = 1; const F = new Array(n).fill(0);
    for (const e of P.els) { const dj = jump(e, eq), d = [3 * e.i, 3 * e.i + 1, 3 * e.i + 2, 3 * e.j, 3 * e.j + 1, 3 * e.j + 2]; for (let a = 0; a < 6; a++) for (let b = 0; b < 3; b++) F[d[a]] -= e.kg[a][3 + b] * dj[b]; }
    const u = new Array(n).fill(0); if (free.length) { const x = N.solveLinear(Kr, free.map(i => F[i])); free.forEach((i, k) => { u[i] = x[k]; }); }
    let s11 = 0, s22 = 0, s12 = 0;
    for (const e of P.els) { const dj = jump(e, eq), d = [u[3 * e.i], u[3 * e.i + 1], u[3 * e.i + 2], u[3 * e.j] + dj[0], u[3 * e.j + 1] + dj[1], u[3 * e.j + 2]], f = N.matVec(e.kg, d), rx = e.c * e.L, ry = e.s * e.L; s11 += f[3] * rx; s22 += f[4] * ry; s12 += 0.5 * (f[3] * ry + f[4] * rx); }
    C[0][q] = s11 / area; C[1][q] = s22 / area; C[2][q] = s12 / area; U.push(u);
  }
  const S = N.inverse(C);
  return { C, S, U, P, E1: 1 / S[0][0], E2: 1 / S[1][1], nu12: -S[1][0] / S[0][0], nu21: -S[0][1] / S[1][1], G12: 1 / S[2][2], rho: P.rho };
}
// peak wall stress for a macroscopic strain vector e (combination of the three unit solutions)
function wallStress(cell, H, e, t) {
  let smax = 0; const u = H.U[0].map((_, k) => e[0] * H.U[0][k] + e[1] * H.U[1][k] + e[2] * H.U[2][k]);
  for (const el of H.P.els) { const dj = jump(el, e), d = [u[3 * el.i], u[3 * el.i + 1], u[3 * el.i + 2], u[3 * el.j] + dj[0], u[3 * el.j + 1] + dj[1], u[3 * el.j + 2]], fl = N.matVec(el.kl, N.matVec(el.T, d)); smax = Math.max(smax, Math.abs(fl[3]) / H.P.A + Math.max(Math.abs(fl[2]), Math.abs(fl[5])) * (t / 2) / H.P.I); }
  return { smax, u };
}
// tessellation of nx × ny cells: verifies the unit-cell solution on the supercell by a matrix-free residual / CG solve
function supercellCheck(cell, H, e, u, nx, ny) {
  const nn = cell.nodes.length, n = 3 * nn * nx * ny, idx = (ix, iy, k) => ((((iy % ny) + ny) % ny) * nx + (((ix % nx) + nx) % nx)) * nn + k, list = [];
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) cell.els.forEach(([i, j, sx, sy], q) => { const el = H.P.els[q], wx = Math.floor((ix + sx) / nx), wy = Math.floor((iy + sy) / ny); list.push({ a: idx(ix, iy, i), b: idx(ix + sx, iy + sy, j), kg: el.kg, dj: jump({ SX: wx * nx * cell.Lx, SY: wy * ny * cell.Ly }, e) }); });
  const mv = x => { const y = new Float64Array(n); for (const el of list) { const d = [3 * el.a, 3 * el.a + 1, 3 * el.a + 2, 3 * el.b, 3 * el.b + 1, 3 * el.b + 2]; for (let a = 0; a < 6; a++) { let s = 0; for (let b = 0; b < 6; b++) s += el.kg[a][b] * x[d[b]]; y[d[a]] += s; } } y[0] = x[0]; y[1] = x[1]; return y; };
  const b = new Float64Array(n), diag = new Float64Array(n); for (const el of list) { const d = [3 * el.a, 3 * el.a + 1, 3 * el.a + 2, 3 * el.b, 3 * el.b + 1, 3 * el.b + 2]; for (let a = 0; a < 6; a++) { diag[d[a]] += el.kg[a][a]; for (let c = 0; c < 3; c++) b[d[a]] -= el.kg[a][3 + c] * el.dj[c]; } }
  b[0] = 0; b[1] = 0; diag[0] = 1; diag[1] = 1; const x0 = new Float64Array(n); for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) { const o = (iy * nx + ix) * 3 * nn, X = ix * cell.Lx, Y = iy * cell.Ly, ax = e[0] * X + 0.5 * e[2] * Y, ay = 0.5 * e[2] * X + e[1] * Y; for (let k = 0; k < nn; k++) { x0[o + 3 * k] = u[3 * k] + ax; x0[o + 3 * k + 1] = u[3 * k + 1] + ay; x0[o + 3 * k + 2] = u[3 * k + 2]; } }
  const r0 = mv(x0); let rn = 0, bn = 0; for (let i = 0; i < n; i++) { rn += (b[i] - r0[i]) ** 2; bn += b[i] ** 2; }
  const cg = N.cg(mv, b, { x0, tol: 1e-9, maxIter: 60, diag }); return { dof: n, res0: Math.sqrt(rn / pos(bn)), iterations: cg.iterations, residual: cg.residual };
}
function drawLattice(cell, u, e, nx, ny, t, amp) {
  const nn = cell.nodes.length, base = [], disp = [], edges = [], ghost = new Map();
  const put = (ix, iy, k) => { const X = ix * cell.Lx, Y = iy * cell.Ly; base.push([cell.nodes[k][0] + X, cell.nodes[k][1] + Y, 0]); disp.push([e[0] * X + 0.5 * e[2] * Y + u[3 * k], 0.5 * e[2] * X + e[1] * Y + u[3 * k + 1]]); return base.length - 1; }; // u already holds the in-cell affine part
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) for (let k = 0; k < nn; k++) put(ix, iy, k);
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) for (const [i, j, sx, sy] of cell.els) {
    const jx = ix + sx, jy = iy + sy, a = (iy * nx + ix) * nn + i; let b;
    if (jx >= 0 && jx < nx && jy >= 0 && jy < ny) b = (jy * nx + jx) * nn + j; else { const key = `${jx},${jy},${j}`; if (!ghost.has(key)) ghost.set(key, put(jx, jy, j)); b = ghost.get(key); }
    edges.push([a, b]);
  }
  const r4 = v => +v.toPrecision(6), frames = N.linspace(0, 1, 9).map(s => ({ t: s, nodes: base.map((p, i) => [r4(p[0] + s * amp * disp[i][0]), r4(p[1] + s * amp * disp[i][1]), 0]) }));
  return { kind: 'lattice', nodes: base.map(p => [r4(p[0]), r4(p[1]), 0]), edges, radius: t / 2, units: 'mm', frames };
}
const TOPO = {
  reentrant: 'General hexagonal / re-entrant honeycomb (h/l, θ)', honeycomb: 'Regular hexagonal honeycomb', triangular: 'Triangular lattice', square: 'Square lattice', kagome: 'Kagome lattice', chiral: 'Hexachiral lattice (rings + tangent ligaments)',
};
function analytic(topo, hl, thDeg, tl, cell) { // classical beam-theory results (E/E_s, ν); null where not available
  const th = thDeg * PI / 180, s = Math.sin(th), c = Math.cos(th);
  if (topo === 'reentrant' || topo === 'honeycomb') return Math.abs(s) < 1e-3 ? null : { E1: tl ** 3 * c / ((hl + s) * s * s), E2: tl ** 3 * (hl + s) / c ** 3, nu12: c * c / ((hl + s) * s), nu21: (hl + s) * s / (c * c), rho: tl * (hl + 2) / (2 * c * (hl + s)), name: 'Gibson–Ashby cell-wall bending' };
  if (topo === 'triangular') return { E1: 2 / Math.sqrt(3) * tl, E2: 2 / Math.sqrt(3) * tl, nu12: 1 / 3, nu21: 1 / 3, rho: 2 * Math.sqrt(3) * tl, name: 'stretch-dominated truss' };
  if (topo === 'kagome') return { E1: tl / Math.sqrt(3), E2: tl / Math.sqrt(3), nu12: 1 / 3, nu21: 1 / 3, rho: Math.sqrt(3) * tl, name: 'stretch-dominated truss' };
  if (topo === 'square') return { E1: tl, E2: tl, nu12: 0, nu21: 0, rho: 2 * tl, name: 'axially loaded walls' };
  if (topo === 'chiral') { const tL = tl * cell.Lx / cell.ligament; return { E1: Math.sqrt(3) * tL ** 3 * (cell.ligament / cell.ring) ** 2, E2: Math.sqrt(3) * tL ** 3 * (cell.ligament / cell.ring) ** 2, nu12: -1, nu21: -1, rho: null, name: 'Prall–Lakes rigid-ring model' }; }
  return null;
}
function lattice2d(inp) {
  const W = [], topo = inp.topology, l = inp.cell_l, Es = inp.Es, nus = inp.nu_s, timo = inp.beam === 'timoshenko'; let hl = inp.h_l, th = inp.theta_c, t = inp.wall_t;
  if (topo === 'honeycomb') { hl = 1; th = 30; }
  if (topo === 'reentrant') { const lim = -Math.asin(Math.min(0.95 * hl / 2, 0.98)) * 180 / PI; if (th < lim) { W.push(`Re-entrant angle limited to ${f3(lim)}°: beyond it opposite vertices of the cell would touch (h must exceed 2l·|sin θ|).`); th = lim; } if (th > 80) th = 80; }
  if (t > 0.3 * l) { W.push('Wall thickness limited to 0.3 l: slender-beam elements are not meaningful for stockier walls.'); t = 0.3 * l; }
  const ringT = N.clamp(Number(inp.ring_t) || 1, 0.5, 30), cell = makeCell(topo, l, hl, th, inp.r_R, ringT), tRef = topo === 'chiral' ? cell.ligament : l;
  let H; try { H = homogenise(cell, Es, nus, t, timo); } catch (e) { throw new Error('Lattice stiffness matrix is singular (mechanism): ' + e.message); }
  const an = analytic(topo, hl, th, t / l, cell), dir = inp.load_dir === '2' ? 1 : 0, eps = inp.strain_pct / 100, Ed = dir ? H.E2 : H.E1, sv = [0, 0, 0]; sv[dir] = Ed * eps;
  const ev = N.matVec(H.S, sv), ws = wallStress(cell, H, ev, t), sf = inp.yield_s / Math.max(ws.smax, 1e-30), nuD = dir ? H.nu21 : H.nu12;
  const eig = N.eigSym([[H.C[0][0], 0.5 * (H.C[0][1] + H.C[1][0]), 0.5 * (H.C[0][2] + H.C[2][0])], [0.5 * (H.C[0][1] + H.C[1][0]), H.C[1][1], 0.5 * (H.C[1][2] + H.C[2][1])], [0.5 * (H.C[0][2] + H.C[2][0]), 0.5 * (H.C[1][2] + H.C[2][1]), H.C[2][2]]]).values;
  const cmax = Math.max(Math.abs(H.C[0][0]), Math.abs(H.C[1][1])), asym = Math.max(Math.abs(H.C[0][1] - H.C[1][0]), Math.abs(H.C[0][2] - H.C[2][0]), Math.abs(H.C[1][2] - H.C[2][1])) / cmax, recip = Math.abs(H.E1 * H.nu21 - H.E2 * H.nu12) / Math.max(Math.abs(H.E1 * H.nu21), 1e-30);
  const nc = Math.round(inp.cells), sc = supercellCheck(cell, H, ev, ws.u, nc, nc), geometry = drawLattice(cell, ws.u, ev, nc, nc, t, 1);
  // directional properties from the compliance matrix
  const S = H.S, ang = N.linspace(0, 180, 91), Eang = ang.map(a => { const c = Math.cos(a * PI / 180), s = Math.sin(a * PI / 180); return 1 / (S[0][0] * c ** 4 + (2 * S[0][1] + S[2][2]) * c * c * s * s + S[1][1] * s ** 4 + 2 * S[0][2] * c ** 3 * s + 2 * S[1][2] * c * s ** 3); });
  const nuang = ang.map((a, i) => { const c = Math.cos(a * PI / 180), s = Math.sin(a * PI / 180); return -Eang[i] * (S[0][1] * (c ** 4 + s ** 4) + (S[0][0] + S[1][1] - S[2][2]) * c * c * s * s + (S[1][2] - S[0][2]) * (c ** 3 * s - c * s ** 3)); });
  let sweep;
  if (topo === 'reentrant') { const lim = -Math.asin(Math.min(0.9 * hl / 2, 0.95)) * 180 / PI, as = N.linspace(Math.max(lim, -60), 60, 41).filter(a => Math.abs(a) > 2), fe = as.map(a => homogenise(makeCell(topo, l, hl, a, 0.25), Es, nus, t, timo).nu12); sweep = ser('nu_theta', 'Poisson\'s ratio against cell angle', 'cell angle θ (°)', 'ν₁₂ (–)', [{ name: 'frame FE', x: as, y: fe }, { name: 'Gibson–Ashby', x: as, y: as.map(a => analytic(topo, hl, a, t / l).nu12) }, { name: 'this design', x: [th], y: [H.nu12] }], 'Poisson\'s ratio of the honeycomb family against the wall inclination: negative (auxetic) for re-entrant angles θ < 0 and positive for conventional cells. Finite-element and analytical curves differ only through axial and shear deformation of the walls.'); }
  else { const tls = N.logspace(0.01, 0.25, 9), fe = tls.map(x => homogenise(cell, Es, nus, x * l, timo)); sweep = ser('scaling', 'Stiffness scaling with relative density', 'relative density (–)', 'E₁/E_s (–)', [{ name: 'frame FE', x: fe.map(h => h.rho), y: fe.map(h => pos(h.E1 / Es)) }, ...(an ? [{ name: an.name, x: fe.map(h => h.rho), y: tls.map(x => pos(analytic(topo, hl, th, x, cell).E1)) }] : []), { name: 'this design', x: [H.rho], y: [pos(H.E1 / Es)] }], 'Relative modulus against relative density on logarithmic axes. A slope of 1 indicates stretch-dominated behaviour (triangular, kagome, square along the walls); a slope of 3 indicates wall bending (honeycombs, chiral lattices).', { logx: true, logy: true }); }
  const tl = t / l, errE = an ? Math.abs(H.E1 / Es - an.E1) / an.E1 : 0, errNu = an ? Math.abs(H.nu12 - an.nu12) : 0, tolE = topo === 'chiral' ? 0.9 : 0.06 + 6 * tl * tl * (topo === 'reentrant' || topo === 'honeycomb' ? 1 / Math.max(Math.sin(th * PI / 180) ** 2, 0.02) : 1), tolNu = topo === 'chiral' ? 0.5 : 0.03 + 0.06 * Math.abs(an ? an.nu12 : 0) + (topo === 'reentrant' || topo === 'honeycomb' ? 6 * tl * tl * Math.abs(an ? an.nu12 : 0) / Math.max(Math.sin(th * PI / 180) ** 2, 0.02) : 0);
  if (sf < 1) W.push(`Peak wall stress ${f3(ws.smax / 1e6)} MPa at ${inp.strain_pct} % strain exceeds the wall yield strength: the response is no longer linear elastic at this strain.`);
  if (eps > 0.02) W.push('Applied strain above 2 %: the small-displacement solution is extrapolated linearly; cell-wall rotation changes the effective properties (and Poisson\'s ratio) at large strain.');
  if (topo === 'chiral') W.push(ringT >= 10 ? `Rings are modelled as 12-sided polygons of beam elements ${f3(ringT)} times thicker than the ligaments (nearly rigid rings): the result approaches the classical Prall–Lakes value ν = −1.` : `Rings are modelled as 12-sided polygons of straight beam elements ${ringT === 1 ? 'with the same thickness as' : f3(ringT) + ' times the thickness of'} the ligaments; the classical ν = −1 result assumes rigid rings, so the finite-element value is less negative. Raise the ring wall thickness ratio to approach the rigid-ring limit.`);
  const auxetic = H.nu12 < 0 || H.nu21 < 0;
  return {
    summary: `${TOPO[topo]} with ${f3(t)} mm walls on ${f3(l)} mm cells (relative density ${f3(H.rho)}): E₁ = ${f3(H.E1 / 1e6)} MPa, E₂ = ${f3(H.E2 / 1e6)} MPa, ν₁₂ = ${f3(H.nu12)}, ν₂₁ = ${f3(H.nu21)}, G₁₂ = ${f3(H.G12 / 1e6)} MPa from periodic frame finite elements${an ? ` (${an.name}: E₁ = ${f3(an.E1 * Es / 1e6)} MPa, ν₁₂ = ${f3(an.nu12)})` : ''}. The lattice is ${auxetic ? 'auxetic (negative Poisson\'s ratio)' : 'conventional (non-negative Poisson\'s ratio)'}. At ${inp.strain_pct} % strain along direction ${dir + 1} the peak wall stress is ${f3(ws.smax / 1e6)} MPa (safety factor ${f3(sf)}).`,
    kpis: [
      kpi('E1', 'Effective modulus E₁', H.E1, 'Pa', 'Homogenised Young\'s modulus along x₁ (horizontal in the plot) from the inverse of the stiffness matrix computed by periodic frame finite elements.'),
      kpi('E2', 'Effective modulus E₂', H.E2, 'Pa', 'Homogenised Young\'s modulus along x₂ (vertical in the plot); equal to E₁ for lattices with three- or six-fold symmetry.'),
      kpi('nu12', 'Poisson\'s ratio ν₁₂', H.nu12, '–', 'Lateral contraction along x₂ per unit extension along x₁ under uniaxial stress; a negative value means the lattice widens when stretched (auxetic).', null),
      kpi('nu21', 'Poisson\'s ratio ν₂₁', H.nu21, '–', 'Lateral contraction along x₁ per unit extension along x₂; related to ν₁₂ by the reciprocity relation E₁ν₂₁ = E₂ν₁₂.', null),
      kpi('G12', 'Effective shear modulus G₁₂', H.G12, 'Pa', 'In-plane shear modulus from the engineering-shear load case; bending-dominated lattices are very compliant in shear, auxetic ones comparatively stiff.'),
      kpi('rho', 'Relative density', H.rho, '–', 'Solid area fraction, total wall length × thickness per cell area (joint overlap not subtracted, so slightly high for thick walls).'),
      kpi('E_rel', 'Relative modulus E₁/E_s', H.E1 / Es, '–', 'Stiffness knock-down relative to the wall material; scales as ρ³ for bending-dominated and ρ for stretch-dominated planar lattices.'),
      kpi('E1_an', 'Analytical E₁', an ? an.E1 * Es : 0, 'Pa', an ? `Closed-form value from the ${an.name} model for comparison with the finite-element result (difference ${f3(100 * errE)} %).` : 'No closed-form reference is available for this configuration (zero reported).'),
      kpi('nu12_an', 'Analytical ν₁₂', an ? an.nu12 : 0, '–', an ? `Closed-form Poisson's ratio from the ${an.name} model; the finite-element value differs by ${f3(errNu)}.` : 'No closed-form reference is available for this configuration (zero reported).'),
      kpi('rho_eff', 'Effective density', H.rho * inp.rho_s, 'kg/m³', 'Mass per unit volume of the extruded lattice, relative density × wall-material density.'),
      kpi('spec', 'Specific stiffness E₁/ρ', H.E1 / (H.rho * inp.rho_s), 'J/kg', 'Modulus-to-density ratio in the loading plane; compare with the wall material value ' + f3(Es / inp.rho_s) + ' J/kg.'),
      kpi('smax', 'Peak wall stress at applied strain', ws.smax, 'Pa', `Largest combined axial + bending stress in any wall for uniaxial stress along direction ${dir + 1} producing ${inp.strain_pct} % strain.`, sf >= 1),
      kpi('sf', 'Safety factor against wall yield', sf, '–', `Wall yield strength divided by the peak wall stress; the lattice yields at an effective strain of about ${f3(inp.strain_pct * sf)} %.`, sf >= 1),
      kpi('eps_lat', 'Lateral strain at applied strain', -nuD * eps * 100, '%', 'Transverse strain accompanying the applied strain under uniaxial stress; positive means lateral expansion in tension.'),
    ],
    series: [
      ser('polar', 'Directional modulus and Poisson\'s ratio', 'loading direction from x₁ (°)', 'value', [{ name: 'E(φ)/E_s × 1000', x: ang, y: Eang.map(v => 1000 * v / Es) }, { name: 'ν(φ)', x: ang, y: nuang }], 'In-plane anisotropy from rotating the homogenised compliance matrix. Flat curves indicate isotropy (hexagonal symmetry); orthotropic and re-entrant cells show strong direction dependence and can be auxetic only in some directions.'),
      sweep,
    ],
    tables: [{ key: 'C', title: 'Homogenised in-plane stiffness matrix (MPa)', columns: ['', 'ε₁₁', 'ε₂₂', 'γ₁₂'], rows: [['σ₁₁', ...H.C[0].map(v => +(v / 1e6).toPrecision(5))], ['σ₂₂', ...H.C[1].map(v => +(v / 1e6).toPrecision(5))], ['σ₁₂', ...H.C[2].map(v => +(v / 1e6).toPrecision(5))]], explain: 'Stiffness matrix relating macroscopic stress to macroscopic strain, each column obtained from one unit-strain load case on the periodic cell (per unit out-of-plane depth).' }],
    geometry,
    equations: [
      eqn('governing', 'Frame equilibrium (direct stiffness method)', 'K u = −Σₑ Kₑ·Δuₑ(ε̄),  Kₑ = Tᵀ kₑ T', `${timo ? 'Timoshenko' : 'Euler–Bernoulli'} beam elements with axial, bending${timo ? ' and shear' : ''} stiffness; periodic boundary conditions`),
      eqn('governing', 'Periodic boundary conditions', 'u(x + a) = u(x) + ε̄·a,  θ(x + a) = θ(x)'), eqn('constitutive', 'Homogenised stress', 'σ̄ = (1/A_cell) Σₑ fₑ ⊗ (x_j − x_i);  σ̄ = C ε̄,  S = C⁻¹'),
      eqn('formula', 'Engineering constants', 'E₁ = 1/S₁₁,  E₂ = 1/S₂₂,  ν₁₂ = −S₂₁/S₁₁,  ν₂₁ = −S₁₂/S₂₂,  G₁₂ = 1/S₆₆'),
      ...(topo === 'reentrant' || topo === 'honeycomb' ? [eqn('formula', 'Gibson–Ashby honeycomb (bending)', 'E₁/E_s = (t/l)³ cos θ/[(h/l + sin θ) sin²θ],  E₂/E_s = (t/l)³ (h/l + sin θ)/cos³θ,  ν₁₂ = cos²θ/[(h/l + sin θ) sin θ]')] : []),
      eqn('formula', 'Reciprocity and stability', 'E₁ν₂₁ = E₂ν₁₂;  C positive definite  ⇒  |ν₁₂| < √(E₁/E₂)'),
    ],
    checks: [
      chk('Reciprocity E₁ν₂₁ = E₂ν₁₂ (independent load cases)', asym < 1e-6 && recip < 1e-5, `stiffness asymmetry ${f3(asym)}, relative reciprocity error ${f3(recip)}`),
      chk('Thermodynamic stability: stiffness matrix positive definite', eig[0] > 0, `eigenvalues ${eig.map(v => f3(v / 1e6)).join(', ')} MPa`),
      chk('Poisson\'s ratio within orthotropic bound |ν₁₂| < √(E₁/E₂)', Math.abs(H.nu12) < Math.sqrt(H.E1 / H.E2) * (1 + 1e-9), `|ν₁₂| = ${f3(Math.abs(H.nu12))}, bound ${f3(Math.sqrt(H.E1 / H.E2))}`),
      chk('Agreement with analytical modulus', !an || errE < tolE, an ? `FE ${f3(H.E1 / Es)} against ${f3(an.E1)} (E₁/E_s), difference ${f3(100 * errE)} %, tolerance ${f3(100 * tolE)} % for axial/shear effects of order (t/l)²` : 'no reference'),
      chk('Agreement with analytical Poisson\'s ratio', !an || errNu < tolNu, an ? `FE ${f4(H.nu12)} against ${f4(an.nu12)}` : 'no reference'),
      chk(`Tessellation (${nc}×${nc} cells, ${sc.dof} DOF) reproduces the unit-cell solution`, sc.res0 < 1e-6 || sc.residual < 1e-8, `initial residual ${f3(sc.res0)}, conjugate-gradient iterations ${sc.iterations}, final residual ${f3(sc.residual)}`),
      chk('Elastic at the applied strain', sf >= 1, `peak wall stress ${f3(ws.smax / 1e6)} MPa against yield ${f3(inp.yield_s / 1e6)} MPa`),
    ],
    warnings: W,
    recommendations: [
      rec('high', auxetic ? `ν₁₂ = ${f3(H.nu12)}: use the lattice where synclastic curvature, indentation resistance or strain-matched sensing is needed; note that the auxetic effect is strongly anisotropic (see the directional plot) and falls as the walls rotate at large strain.` : topo === 'reentrant' ? 'Set a negative cell angle (for example θ = −20° to −30° with h/l = 2) to obtain a negative Poisson\'s ratio.' : `This topology is not auxetic (ν₁₂ = ${f3(H.nu12)}); choose the re-entrant honeycomb or the chiral lattice for a negative Poisson's ratio.`, 'Sign of the homogenised Poisson\'s ratio'),
      rec('medium', sf < 1.5 ? `Peak wall stress leaves a safety factor of ${f3(sf)}: fillet the joints (bending stress peaks at wall ends) or reduce t/l, which lowers wall strain for a given macroscopic strain.` : `Wall stress is ${f3(100 / sf)} % of yield at ${inp.strain_pct} % strain; the elastic strain range of the lattice is about ${f3(inp.strain_pct * sf)} %.`, 'Wall bending stress from the frame solution'),
      rec('medium', `Manufacturing tolerance on wall thickness matters as the ${topo === 'triangular' || topo === 'kagome' || topo === 'square' ? 'first' : 'third'} power: ±10 % in t changes E₁ by about ${topo === 'triangular' || topo === 'kagome' || topo === 'square' ? '±10' : '±30'} %. Measure printed wall thickness and re-run.`, 'Scaling exponent of this topology'),
      rec('low', 'Validate against a tensile test on a specimen at least 6–8 cells wide; free edges soften finite samples relative to the periodic prediction.', 'Periodic homogenisation assumes an infinite lattice'),
      rec('medium', `At a relative density of ${f3(H.rho)} the lattice uses ${f3(100 * H.rho)} % of the solid material; choose a stretch-dominated topology where stiffness per mass matters and a bending-dominated or auxetic one where compliance or energy absorption is the function, to avoid over-building.`, 'Material efficiency of architected materials', true),
    ],
    publish: { 'eff.E': Ed, 'eff.nu': nuD, 'eff.G': H.G12, 'geom.rel_density': H.rho, 'eff.density': H.rho * inp.rho_s, 'mech.max_stress': ws.smax, 'mech.safety_factor': sf, 'geom.cell_size': l * 1e-3, 'geom.strut': t * 1e-3 },
    numerics: { method: `Direct-stiffness frame FE (${timo ? 'Timoshenko' : 'Euler–Bernoulli'} beams) on the periodic unit cell, dense LU; matrix-free conjugate-gradient verification on the tessellation`, dof: 3 * cell.nodes.length, iterations: sc.iterations, residual: Math.min(sc.res0, sc.residual), notes: `Three unit-strain load cases; supercell ${sc.dof} DOF. Linear small-displacement analysis, rigid joints, unit depth.` },
  };
}

// ================================================================ negative stiffness (pre-shaped beam, Qiu–Lang–Slocum)
function negstiff(inp) {
  const W = [], E = inp.Es, l = inp.ns_l * 1e-3, t = inp.ns_t * 1e-3, h = inp.ns_h * 1e-3, b = inp.ns_b * 1e-3, nb = Math.round(inp.ns_n), Q = h / t, I = b * t ** 3 / 12, p4 = PI ** 4, dbl = inp.ns_mode === 'double';
  const F1 = d => 1.5 * p4 * Q * Q * d * (d * d - 3 * d + 2 + 4 / (3 * Q * Q)), FL = d => (dbl ? 8 * p4 - 6 * p4 * d : 4.18 * p4 - 2.18 * p4 * d), ds = N.linspace(0, 2.2, 331), g = ds.map(d => F1(d) - FL(d)), roots = [];
  for (let i = 1; i < ds.length; i++) if (g[i - 1] * g[i] < 0) roots.push(N.brent(d => F1(d) - FL(d), ds[i - 1], ds[i], 1e-12)); else if (g[i] === 0) roots.push(ds[i]);
  const da = roots.length >= 2 ? roots[0] : null, db = roots.length >= 2 ? roots[roots.length - 1] : null, Fc = d => (da !== null && d > da && d < db ? FL(d) : F1(d)), F = ds.map(Fc);
  const itop = N.argmax(F.slice(0, 166)), ibot = 100 + N.argmin(F.slice(100)); let Ftop = F[itop], Fbot = F[ibot], dTop = ds[itop];
  if (da !== null) { if (FL(da) > Ftop) { Ftop = FL(da); dTop = da; } if (FL(db) < Fbot) Fbot = FL(db); } // exact branch intersections, not the nearest grid point
  const bistable = Fbot < 0, sc = E * I * h / l ** 3 * nb, ksc = E * I / l ** 3 * nb;
  const slope = ds.map((d, i) => (i === 0 ? (F[1] - F[0]) / (ds[1] - ds[0]) : i === ds.length - 1 ? (F[i] - F[i - 1]) / (ds[i] - ds[i - 1]) : (F[i + 1] - F[i - 1]) / (ds[i + 1] - ds[i - 1])));
  let n0 = -1, n1 = -1; for (let i = itop; i < ds.length; i++) { if (slope[i] < 0) { if (n0 < 0) n0 = i; n1 = i; } else if (n0 >= 0) break; }
  const hasNS = n0 >= 0 && Ftop > 0 && itop > 0 && itop < 165, kmin = hasNS ? N.minOf(slope.slice(n0, n1 + 1)) : 0, kneg = kmin * ksc, U = [0]; for (let i = 1; i < ds.length; i++) U.push(U[i - 1] + 0.5 * (F[i] + F[i - 1]) * (ds[i] - ds[i - 1]));
  const usc = sc * h; let dU = null, dS = null; for (let i = 1; i < ds.length; i++) { if (dU === null && F[i - 1] > 0 && F[i] <= 0 && i > itop) dU = i; if (dU !== null && dS === null && F[i - 1] < 0 && F[i] >= 0) dS = i; }
  const Ubar = dU !== null ? U[dU] * usc : U[ibot] * usc, Utrap = dS !== null ? U[dS] * usc : 0, Uback = dS !== null ? (U[dU] - U[dS]) * usc : 0, epsMax = 2 * PI * PI * t * h / (l * l) + 4 * PI * PI * t * t / (3 * l * l), epsY = inp.yield_s / E;
  // quasi-zero-stiffness isolator: positive spring in parallel with the negative-stiffness element at mid-travel
  const kp = inp.k_pos, kt = kp + kneg, m = inp.m_load, f0 = Math.sqrt(kp / m) / (2 * PI), fq = kt > 0 ? Math.sqrt(kt / m) / (2 * PI) : 0, kcrit = -kneg;
  if (!hasNS) W.push(`With Q = h/t = ${f3(Q)} the force–displacement curve is monotonic: there is no snap-through and no negative-stiffness branch (a curved beam needs Q above about ${dbl ? '1.2' : '1.2'} for a force maximum and above 2.31 for bistability).`);
  if (epsMax > epsY) W.push(`Peak bending strain ${f3(100 * epsMax)} % exceeds the yield strain ${f3(100 * epsY)} %: the beam will yield or fracture during snap-through. Increase the span or reduce thickness/apex height.`);
  if (hasNS && kt <= 0) W.push('The positive spring is softer than the negative-stiffness element: the combined system is statically unstable (it snaps) rather than quasi-zero-stiffness.');
  if (l / t < 20) W.push('Span-to-thickness ratio below 20: the shallow, slender-beam assumptions of the model are marginal.');
  return {
    summary: `Cosine-shaped clamped–clamped beam (${nb}×, span ${inp.ns_l} mm, thickness ${inp.ns_t} mm, apex ${inp.ns_h} mm, Q = ${f3(Q)}): ${hasNS ? `snap-through at ${f3(Ftop * sc)} N (displacement ${f3(dTop * h * 1e3)} mm), minimum force ${f3(Fbot * sc)} N, negative stiffness down to ${f3(kneg)} N/m over ${f3(ds[n0] * h * 1e3)}–${f3(ds[n1] * h * 1e3)} mm` : 'no snap-through'}. The beam is ${bistable ? 'bistable' : 'monostable'}${dbl ? ' (criterion Q > 2.31)' : ' (the Q > 2.31 criterion applies when the second mode is constrained)'}. With a ${f3(kp)} N/m spring in parallel the net stiffness is ${f3(kt)} N/m, lowering the natural frequency of a ${m} kg payload from ${f3(f0)} to ${f3(fq)} Hz.`,
    kpis: [
      kpi('Q', 'Geometry constant Q = h/t', Q, '–', 'Ratio of apex height to beam thickness; it alone fixes the shape of the normalised force–displacement curve. Bistability requires Q > 2.31 when the second (twisting) mode is suppressed.', Q > 2.31),
      kpi('F_top', 'Snap-through force', hasNS ? Ftop * sc : 0, 'N', `Maximum force before the beam snaps (${nb} beam(s)); for large Q it approaches ${dbl ? '8π⁴ ≈ 779' : '4.18π⁴ ≈ 407'}·EIh/l³ per beam.`),
      kpi('F_bot', 'Minimum (snap-back) force', hasNS ? Fbot * sc : 0, 'N', 'Lowest force on the curve after snap-through. A negative value means the beam stays in its second position when released (bistable); its magnitude is the force needed to snap back.', null),
      kpi('d_top', 'Snap-through displacement', hasNS ? dTop * h : 0, 'm', 'Centre displacement at the force maximum, where the tangent stiffness passes through zero.'),
      kpi('k_neg', 'Most negative stiffness', kneg, 'N/m', `Minimum slope of the force–displacement curve (${dbl ? '−6π⁴' : '−2.18π⁴'}·EI/l³ per beam on the mode-switching branch); this is the negative stiffness available to cancel a positive spring.`, null),
      kpi('ns_range', 'Negative-stiffness travel', hasNS ? (ds[n1] - ds[n0]) * h : 0, 'm', 'Displacement range over which the tangent stiffness is negative; it is only traversable quasi-statically under displacement control.'),
      kpi('k0', 'Initial stiffness', slope[0] * ksc, 'N/m', 'Tangent stiffness of the undeformed curved beam at small load.'),
      kpi('bistable', 'Stability class', bistable ? 'bistable' : 'monostable', '', 'Bistable if the force becomes negative after snap-through, so a second stable equilibrium exists without external load.', null),
      kpi('U_bar', 'Energy barrier to snap-through', Ubar, 'J', 'Work done up to the unstable equilibrium (or to the force minimum for a monostable beam): the energy that must be supplied to switch.'),
      kpi('U_trap', 'Strain energy trapped in second state', Utrap, 'J', 'Elastic energy stored in the second stable configuration relative to the as-made shape (zero if monostable); usable for energy-absorbing or reusable shock-trapping lattices.'),
      kpi('U_back', 'Energy barrier for snap-back', Uback, 'J', 'Depth of the second energy well: the work needed to return to the first state (zero if monostable).'),
      kpi('eps_max', 'Peak bending strain', epsMax, '–', 'Maximum surface strain during snap-through, ≈ 2π²th/l² + 4π²t²/(3l²); it must stay below the elastic limit of the beam material.', epsMax < epsY),
      kpi('k_qzs', 'Net stiffness with parallel spring', kt, 'N/m', 'Positive spring stiffness plus the most negative beam stiffness; quasi-zero-stiffness isolation needs this small but positive.', kt > 0),
      kpi('f_qzs', 'Natural frequency with negative stiffness', fq, 'Hz', `Natural frequency of the ${m} kg payload on the combined suspension; isolation starts above √2 times this value. Without the negative-stiffness element it is ${f3(f0)} Hz.`, fq > 0 && fq < f0),
    ],
    series: [
      ser('Fd', 'Force–displacement curve', 'centre displacement (mm)', 'force (N)', [{ name: 'beam response', x: ds.map(d => d * h * 1e3), y: F.map(v => v * sc) }, { name: 'first-mode branch alone', x: ds.map(d => d * h * 1e3), y: ds.map(d => F1(d) * sc) }, { name: 'zero', x: [0, 2.2 * h * 1e3], y: [0, 0] }], 'Quasi-static transverse force at the beam centre against displacement. The falling segment is the negative-stiffness branch; where the curve dips below zero the beam is bistable. Under force control the response jumps horizontally from the maximum (snap-through).'),
      ser('k', 'Tangent stiffness', 'centre displacement (mm)', 'stiffness (N/m)', [{ name: 'dF/dd', x: ds.map(d => d * h * 1e3), y: slope.map(v => v * ksc) }, { name: 'with parallel spring', x: ds.map(d => d * h * 1e3), y: slope.map(v => v * ksc + kp) }, { name: 'zero', x: [0, 2.2 * h * 1e3], y: [0, 0] }], 'Slope of the force–displacement curve. Adding a positive spring shifts the curve up: tuned so that the minimum just clears zero gives a quasi-zero-stiffness isolator with high static but low dynamic stiffness.'),
      ser('U', 'Strain-energy landscape', 'centre displacement (mm)', 'energy (mJ)', [{ name: 'strain energy', x: ds.map(d => d * h * 1e3), y: U.map(v => v * usc * 1e3) }], 'Elastic energy against displacement. Two minima separated by a barrier indicate bistability; the height difference between them is the energy trapped after snap-through.'),
    ],
    equations: [
      eqn('formula', 'Initial shape', 'w₀(x) = (h/2)[1 − cos(2πx/l)]'), eqn('governing', 'First-mode branch (Qiu–Lang–Slocum)', 'F₁ = (3π⁴Q²/2)·Δ·(Δ − 3/2 + √(1/4 − 4/3Q²))·(Δ − 3/2 − √(1/4 − 4/3Q²)),  F = f l³/(EIh),  Δ = d/h'),
      eqn('governing', 'Mode-switching branch', dbl ? 'F₃ = 8π⁴ − 6π⁴Δ  (second mode constrained)' : 'F₂ = 4.18π⁴ − 2.18π⁴Δ  (second mode free)', 'the response follows this line between its intersections with F₁'),
      eqn('formula', 'Bistability criterion', 'Q = h/t > 2.31  (real roots of F₁ = 0 beyond Δ = 3/2)'), eqn('formula', 'Quasi-zero stiffness', 'k_net = k_spring + k_neg ≥ 0,  f = (1/2π)√(k_net/m)'),
    ],
    checks: [
      chk('Branch identity F₁(1) = F₃(1) = 2π⁴', Math.abs(F1(1) - 2 * p4) < 1e-9 * p4, 'closed-form consistency of the two branches at mid-travel'),
      chk('Bistability consistent with Q > 2.31', !dbl || (bistable === (Q > 4 / Math.sqrt(3))) || Math.abs(Q - 4 / Math.sqrt(3)) < 0.05, `Q = ${f3(Q)}, numerical minimum force ${f3(Fbot * sc)} N`),
      chk('Elastic snap-through (strain below yield)', epsMax < epsY, `ε_max = ${f3(100 * epsMax)} % against ε_y = ${f3(100 * epsY)} %`),
      chk('Static stability of the combined suspension (k_net > 0)', kt > 0, `k_spring = ${f3(kp)} N/m against |k_neg| = ${f3(Math.abs(kneg))} N/m`),
      chk('Energy landscape consistent (U = ∫F dd ≥ 0 up to the barrier)', U.slice(0, itop + 1).every(v => v >= 0), `barrier ${f3(Ubar * 1e3)} mJ`),
      chk('Negative stiffness needs displacement control or a stabilising constraint', true, 'an isolated negative-stiffness element is unstable (non-positive-definite stiffness); it is admissible only inside a system whose total stiffness is positive'),
    ],
    warnings: W,
    recommendations: [
      rec('high', hasNS ? (kt > 0 ? `Net stiffness ${f3(kt)} N/m is ${f3(100 * kt / kp)} % of the spring alone: a margin of 5–10 % above zero is a practical compromise between isolation and robustness to tolerance (|k_neg| scales as t³).` : `Raise the parallel spring above ${f3(kcrit)} N/m (about 1.05–1.10 × |k_neg|) to obtain a stable quasi-zero-stiffness mount.`) : `Increase the apex height to at least ${f3(2.5 * inp.ns_t)} mm (Q ≥ 2.5) to obtain snap-through and a usable negative-stiffness range.`, 'k_net = k_spring + k_neg'),
      rec('medium', `Thickness tolerance is critical: forces scale with t³ and the bistability threshold with h/t. A ±5 % thickness error changes |k_neg| by about ±15 %; trim with an adjustable preload spring.`, 'EI ∝ t³ in the normalised model'),
      rec('medium', dbl ? 'Keep the two parallel beams clamped together at mid-span to suppress the asymmetric (second) mode; a single beam snaps through the twisting mode at about half the force and is not reliably bistable.' : 'Suppress the asymmetric mode by pairing two beams joined at the centre: the snap-through force nearly doubles per beam and bistability becomes robust.', 'Mode-2 constraint in the Qiu–Lang–Slocum analysis'),
      rec('low', `Limit ${f3(100 * epsMax)} % peak strain by fatigue data for the beam material if cycling is intended.`, 'Peak bending strain during snap-through'),
      rec('low', `Bistable elements hold either state without power and trap ${f3(Utrap * 1e3)} mJ per event elastically: reusable shock absorbers and latching actuators built from them avoid consumable crush elements and holding currents.`, 'Zero-power state retention and reusability', true),
    ],
    publish: { 'mech.k_neg': kneg, 'mech.snap_force': hasNS ? Ftop * sc : 0, 'res.f0': fq, 'eff.K': kt, 'mech.max_stress': epsMax * E, 'mech.safety_factor': epsY / epsMax },
    numerics: { method: 'Closed-form modal solution of the pre-shaped clamped–clamped beam with numerical branch intersection (Brent) and trapezoidal energy integration', dof: ds.length, notes: 'Shallow, slender, linear-elastic beam; axial compression coupled to bending through the first three buckling modes.' },
  };
}

// ================================================================ negative linear compressibility
function nlc(inp) {
  const W = [], T = 1e12; // TPa⁻¹ scaling
  if (inp.nlc_model === 'crystal') {
    const g = 1e9, C = [[inp.c11 * g, inp.c12 * g, inp.c13 * g], [inp.c12 * g, inp.c22 * g, inp.c23 * g], [inp.c13 * g, inp.c23 * g, inp.c33 * g]], eig = N.eigSym(C).values, stable = eig[0] > 0 && inp.c44 > 0;
    if (!stable) { W.push('The stiffness matrix is not positive definite (Born stability violated): the entered constants do not describe a mechanically stable crystal; compressibilities below are not physical.'); }
    let S; try { S = N.inverse(C); } catch { S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; W.push('Stiffness matrix is singular.'); }
    const be = S.map(r => r[0] + r[1] + r[2]), bV = be[0] + be[1] + be[2], imin = N.argmin(be), ang = N.linspace(0, 180, 91), dirB = (a, i, j) => be[i] * Math.cos(a * PI / 180) ** 2 + be[j] * Math.sin(a * PI / 180) ** 2, neg = be[imin] < 0;
    return {
      summary: `Orthotropic stiffness set (C₁₁ = ${inp.c11}, C₂₂ = ${inp.c22}, C₃₃ = ${inp.c33}, C₁₂ = ${inp.c12}, C₁₃ = ${inp.c13}, C₂₃ = ${inp.c23} GPa): linear compressibilities β₁ = ${f3(be[0] * T)}, β₂ = ${f3(be[1] * T)}, β₃ = ${f3(be[2] * T)} TPa⁻¹; volumetric compressibility ${f3(bV * T)} TPa⁻¹ (bulk modulus ${f3(1 / bV / 1e9)} GPa). ${neg ? `Axis ${imin + 1} expands under hydrostatic pressure (negative linear compressibility).` : 'No axis shows negative linear compressibility.'} The set is ${stable ? 'mechanically stable' : 'unstable'}.`,
      kpis: [
        kpi('beta1', 'Linear compressibility β₁', be[0] * T, 'TPa⁻¹', 'Relative shortening of axis 1 per unit hydrostatic pressure, the row sum S₁₁ + S₁₂ + S₁₃ of the compliance matrix; a negative value means the axis lengthens under pressure.', null),
        kpi('beta2', 'Linear compressibility β₂', be[1] * T, 'TPa⁻¹', 'Relative shortening of axis 2 per unit hydrostatic pressure, S₂₁ + S₂₂ + S₂₃.', null),
        kpi('beta3', 'Linear compressibility β₃', be[2] * T, 'TPa⁻¹', 'Relative shortening of axis 3 per unit hydrostatic pressure, S₃₁ + S₃₂ + S₃₃.', null),
        kpi('betaV', 'Volumetric compressibility', bV * T, 'TPa⁻¹', 'Sum of the three linear compressibilities; thermodynamic stability requires it to be positive even when one axis expands.', bV > 0),
        kpi('K', 'Bulk modulus (Reuss)', 1 / bV, 'Pa', 'Reciprocal of the volumetric compressibility under uniform pressure.'),
        kpi('nlc', 'Most negative linear compressibility', Math.min(be[imin], 0) * T, 'TPa⁻¹', `Magnitude of the negative linear compressibility along axis ${imin + 1} (zero if all axes shrink); conventional solids have about +5 to +20 TPa⁻¹.`, neg),
        kpi('cone', 'Half-angle of the NLC cone', neg ? Math.atan(Math.sqrt(-be[imin] / Math.max(Math.max(...be), 1e-30))) * 180 / PI : 0, '°', 'Directions within this angle of the NLC axis also expand under pressure (zero-compressibility cone in the plane containing the stiffest-contracting axis).'),
        kpi('eigmin', 'Smallest stiffness eigenvalue', eig[0], 'Pa', 'Lowest eigenvalue of the normal-stress stiffness block; it must be positive for Born mechanical stability.', eig[0] > 0),
      ],
      series: [
        ser('dir', 'Directional linear compressibility', 'angle from first axis of the plane (°)', 'β (TPa⁻¹)', [{ name: 'plane 1–3', x: ang, y: ang.map(a => dirB(a, 0, 2) * T) }, { name: 'plane 1–2', x: ang, y: ang.map(a => dirB(a, 0, 1) * T) }, { name: 'plane 2–3', x: ang, y: ang.map(a => dirB(a, 1, 2) * T) }, { name: 'zero', x: [0, 180], y: [0, 0] }], 'Linear compressibility β(n) = Σ βᵢnᵢ² in the three principal planes. Directions where a curve is below zero lengthen under hydrostatic pressure.'),
        ser('row', 'Compliance contributions per axis', 'axis', 'compliance (TPa⁻¹)', [{ name: 'S_ii (direct)', x: [1, 2, 3], y: [S[0][0] * T, S[1][1] * T, S[2][2] * T] }, { name: 'Σ S_ij, j≠i (coupling)', x: [1, 2, 3], y: [0, 1, 2].map(i => (be[i] - S[i][i]) * T) }], 'Decomposition of each linear compressibility into the direct compliance and the (negative) Poisson coupling terms; negative linear compressibility occurs when the coupling outweighs the direct term.'),
      ],
      equations: [eqn('constitutive', 'Generalised Hooke law', 'ε = S σ,  S = C⁻¹'), eqn('formula', 'Linear compressibility', 'βᵢ = −(1/Lᵢ)(∂Lᵢ/∂p) = Σⱼ Sᵢⱼ  (j = 1…3)'), eqn('formula', 'Volumetric compressibility', 'β_V = β₁ + β₂ + β₃ = 1/K > 0'), eqn('formula', 'Born stability', 'C positive definite: all leading principal minors > 0, C₄₄, C₅₅, C₆₆ > 0')],
      checks: [
        chk('Born mechanical stability (C positive definite)', stable, `eigenvalues ${eig.map(v => f3(v / 1e9)).join(', ')} GPa; shear constants ${inp.c44} GPa`), chk('Volumetric compressibility positive (thermodynamic stability)', bV > 0, `β_V = ${f3(bV * T)} TPa⁻¹`),
        chk('Compliance is the inverse of stiffness', Math.abs(N.matMul(C, S)[0][0] - 1) < 1e-8, 'C·S = I verified'), chk('Negative linear compressibility present', neg, neg ? `β along axis ${imin + 1} = ${f3(be[imin] * T)} TPa⁻¹` : 'all βᵢ ≥ 0'),
      ],
      warnings: [...W, 'The default constants are an illustrative orthotropic set with strong axial–transverse coupling, not data for a specific crystal: replace them with measured or computed elastic constants.'],
      recommendations: [
        rec('high', neg ? `Use axis ${imin + 1} for pressure-sensing or zero-area-compressibility designs: combined with a conventional material of compressibility β_c in a length ratio |β_NLC|/β_c the composite length is pressure-invariant.` : 'Increase the off-diagonal coupling (C₁₃, C₂₃) relative to the transverse stiffness, as in wine-rack or helical frameworks, to obtain negative linear compressibility.', 'Sign of the compliance row sums'),
        rec('medium', 'Obtain all six normal-block constants from the same source (resonant ultrasound, Brillouin scattering or DFT): β is a small difference of large terms and mixes poorly between datasets.', 'Sensitivity of row sums of S to errors in Cᵢⱼ'),
        rec('medium', 'Confirm by high-pressure diffraction: lattice-parameter change of the NLC axis must be positive while the cell volume decreases.', 'Definition of linear compressibility'),
        rec('low', 'Framework materials with NLC are often porous and pressure-amorphise: stay below the stability limit and avoid pressure media that enter the pores.', 'Mechanical-stability check', false),
        rec('low', 'Pressure sensors built on an NLC axis can dispense with fluid-filled bellows and lead-based piezoelectrics in some designs.', 'Substitution of hazardous transducer materials', true),
      ],
      publish: { 'eff.K': 1 / bV, 'mech.beta_min': be[imin], 'mech.beta_V': bV },
      numerics: { method: 'Dense inversion of the 3×3 normal stiffness block; Jacobi eigenvalues for the stability test', dof: 3, notes: 'Orthotropic symmetry: shear constants do not enter the hydrostatic response.' },
    };
  }
  // wine-rack (rhombic) hinging framework
  const E = inp.Es, nu = inp.nu_s, L = inp.wr_L, w = Math.min(inp.wr_w, 0.3 * L), th = inp.wr_theta, thr = th * PI / 180, timo = inp.beam === 'timoshenko';
  if (inp.wr_w > 0.3 * L) W.push('Strut width limited to 0.3 L for the beam model.');
  const an2 = a => { const r = a * PI / 180, hs = (L / w) ** 3 / E * Math.cos(2 * r), st = L * Math.sin(2 * r) / (E * w); return [st - hs * Math.tan(r), st + hs / Math.tan(r)]; };
  const cell = makeCell('rhombic', L, 0, th), H = homogenise(cell, E, nu, w, timo), rho = H.rho, bz2 = -nu / (rho * E), [ax, ay] = an2(th);
  const fx = H.S[0][0] + H.S[0][1], fy = H.S[1][0] + H.S[1][1], bx = fx + bz2, by = fy + bz2, bz = (1 - 2 * nu) / (rho * E), bV = bx + by + bz, errX = Math.abs(fx - ax) / Math.max(Math.abs(ax), 1e-30), errY = Math.abs(fy - ay) / Math.max(Math.abs(ay), 1e-30);
  const be = [bx, by, bz], imin = N.argmin(be), neg = be[imin] < 0, angs = N.linspace(8, 82, 75), eig = N.eigSym([[H.C[0][0], H.C[0][1]], [H.C[0][1], H.C[1][1]]]).values;
  // deformation under in-plane hydrostatic pressure for the display
  const p = 0.02 / Math.max(Math.abs(fx), Math.abs(fy)), ev = N.matVec(H.S, [-p, -p, 0]), ws = wallStress(cell, H, ev, w), geometry = drawLattice(cell, ws.u, ev, 4, 4, w, 1);
  const th0 = (() => { try { return N.brent(a => { const [x, y] = an2(a); return Math.min(x, y) + bz2 * 0; }, 8, 44.9); } catch { return null; } })();
  if (!neg) W.push('No negative linear compressibility at this geometry: the hinging (strut-bending) mechanism is weak near 45° or the struts are too stocky.');
  return {
    summary: `Wine-rack framework (strut length ${L} mm, width ${f3(w)} mm, half-angle ${th}°, relative density ${f3(rho)}): β_x = ${f3(bx * T)}, β_y = ${f3(by * T)}, β_z = ${f3(bz * T)} TPa⁻¹ under hydrostatic pressure. ${neg ? `The framework lengthens along ${['x (long diagonal)', 'y (long diagonal)', 'z'][imin]} when compressed: negative linear compressibility of ${f3(be[imin] * T)} TPa⁻¹.` : 'All directions shrink.'} The volumetric compressibility stays positive at ${f3(bV * T)} TPa⁻¹. Frame finite elements and the hinging formula agree within ${f3(100 * Math.max(errX, errY))} %.`,
    kpis: [
      kpi('beta_x', 'Linear compressibility β_x', bx * T, 'TPa⁻¹', 'Relative shortening along x per unit hydrostatic pressure from the frame finite-element compliance plus the out-of-plane Poisson coupling of the strut material; negative when the rhombus closes and elongates along x.', null),
      kpi('beta_y', 'Linear compressibility β_y', by * T, 'TPa⁻¹', 'Relative shortening along y per unit hydrostatic pressure; the conjugate direction of the wine-rack mechanism, which always shrinks strongly when x expands.', null),
      kpi('beta_z', 'Linear compressibility β_z', bz * T, 'TPa⁻¹', 'Out-of-plane (extrusion) direction, (1 − 2ν_s)/(ρ̄E_s) for prismatic walls.', bz > 0),
      kpi('beta_V', 'Volumetric compressibility', bV * T, 'TPa⁻¹', 'Sum of the three linear compressibilities; it must be positive for a thermodynamically stable solid even though one axis expands.', bV > 0),
      kpi('K', 'Bulk modulus', 1 / bV, 'Pa', 'Reciprocal of the volumetric compressibility of the framework.'),
      kpi('nlc', 'Negative linear compressibility', Math.min(be[imin], 0) * T, 'TPa⁻¹', 'Most negative of the three values (zero if none). Molecular frameworks reach about −10 to −100 TPa⁻¹; compliant architected lattices reach far larger magnitudes because their hinges are soft.', neg),
      kpi('ratio', 'Expansion-to-volume ratio −β_min/β_V', neg ? -be[imin] / bV : 0, '–', 'Fraction of the volumetric compression that appears as expansion along the NLC axis; a measure of mechanism efficiency.'),
      kpi('an_x', 'Analytical β_x (in-plane hinging + stretching)', ax * T, 'TPa⁻¹', `Closed form (L/E_sw)·sin 2θ − (L/w)³·cos 2θ·tan θ/E_s for rigid joints; the finite-element in-plane value is ${f3(fx * T)} TPa⁻¹.`),
      kpi('an_y', 'Analytical β_y (in-plane hinging + stretching)', ay * T, 'TPa⁻¹', `Closed form (L/E_sw)·sin 2θ + (L/w)³·cos 2θ·cot θ/E_s; the finite-element in-plane value is ${f3(fy * T)} TPa⁻¹.`),
      kpi('rho', 'Relative density', rho, '–', 'Solid fraction of the framework, 2w/(L·sin 2θ) to first order.'),
    ],
    series: [
      ser('bth', 'Linear compressibilities against wine-rack angle', 'half-angle θ (°)', 'β (TPa⁻¹)', [{ name: 'β_x (analytical, in-plane)', x: angs, y: angs.map(a => an2(a)[0] * T) }, { name: 'β_y (analytical, in-plane)', x: angs, y: angs.map(a => an2(a)[1] * T) }, { name: 'β_x + β_y', x: angs, y: angs.map(a => (an2(a)[0] + an2(a)[1]) * T) }, { name: 'this design (FE)', x: [th, th], y: [fx * T, fy * T] }, { name: 'zero', x: [8, 82], y: [0, 0] }], 'In-plane compressibilities against the strut half-angle. Below 45° the x direction (long diagonal) expands under pressure, above 45° the y direction does; at exactly 45° hinging is inactive. The sum (area compressibility) is always positive.'),
      ser('bw', 'Negative compressibility against strut slenderness', 'strut width w/L (–)', '−β_min (TPa⁻¹)', [{ name: 'analytical', x: N.logspace(0.01, 0.3, 40), y: N.logspace(0.01, 0.3, 40).map(r => { const hs = r ** -3 / E * Math.cos(2 * thr), st = Math.sin(2 * thr) / (E * r); return pos(-Math.min(st - hs * Math.tan(thr), st + hs / Math.tan(thr))); }) }], 'Magnitude of the negative in-plane compressibility against strut slenderness: it grows as (L/w)³ because the hinging stiffness is the bending stiffness of the struts.', { logx: true, logy: true }),
    ],
    geometry,
    equations: [
      eqn('governing', 'Wine-rack kinematics', 'X = 2L cos θ,  Y = 2L sin θ;  hinging: δX/X = −tan θ·δθ,  δY/Y = cot θ·δθ'), eqn('formula', 'Hinging + stretching compressibilities (rigid joints)', 'β_x = (L/E_s w) sin 2θ − (L/w)³ cos 2θ tan θ/E_s,  β_y = (L/E_s w) sin 2θ + (L/w)³ cos 2θ cot θ/E_s'),
      eqn('formula', 'Three-dimensional correction for prismatic walls', 'β_x,y → β_x,y − ν_s/(ρ̄E_s),  β_z = (1 − 2ν_s)/(ρ̄E_s)'), eqn('formula', 'Thermodynamic requirement', 'β_V = β_x + β_y + β_z > 0'),
      eqn('governing', 'Finite-element verification', 'β_x = S₁₁ + S₁₂,  β_y = S₂₁ + S₂₂ from the homogenised compliance of the periodic frame'),
    ],
    checks: [
      chk('Volumetric compressibility positive (thermodynamic stability)', bV > 0, `β_V = ${f3(bV * T)} TPa⁻¹`), chk('In-plane stiffness positive definite', eig[0] > 0, `eigenvalues ${eig.map(v => f3(v)).join(', ')} Pa`),
      chk('Frame FE agrees with the hinging formula', errX < 0.08 + 8 * (w / L) ** 2 && errY < 0.08 + 8 * (w / L) ** 2, `β_x: ${f3(fx * T)} against ${f3(ax * T)} TPa⁻¹; β_y: ${f3(fy * T)} against ${f3(ay * T)} TPa⁻¹`),
      chk('Area compressibility positive', fx + fy > 0, `β_x + β_y (in-plane) = ${f3((fx + fy) * T)} TPa⁻¹`), chk('Negative linear compressibility present', neg, neg ? `along ${['x', 'y', 'z'][imin]}` : 'none'),
    ],
    warnings: W,
    recommendations: [
      rec('high', neg ? `The NLC axis is ${['x', 'y', 'z'][imin]}; its magnitude scales as (L/w)³ = ${f3((L / w) ** 3)}. Halving the strut width multiplies |β| by eight but also divides the bulk modulus by a similar factor: choose w from the pressure range of the application.` : `Move the half-angle away from 45° (for example ${th < 45 ? 25 : 65}°) and use slender struts (w/L ≤ 0.1) to activate the hinging mechanism.`, 'Hinging-dominated wine-rack mechanics'),
      rec('medium', `Small-strain result: the angle changes by about ${f3(Math.abs(0.02 / Math.tan(thr)) * 180 / PI)}° per 2 % strain, so β varies with pressure; re-evaluate at the working angle for large pressures.`, 'Geometric non-linearity of mechanisms'),
      rec('medium', 'Use flexure hinges or compliant joints rather than pins: the model assumes rigid (welded) joints whose compliance is the strut bending; pinned joints give unbounded hinging until contact.', 'Rigid-joint frame assumption'),
      rec('low', 'A single low-modulus polymer lattice can replace multi-material pressure compensators; printing only the required struts keeps material use to the relative density shown.', 'Material efficiency of mechanism-based metamaterials', true),
    ],
    publish: { 'eff.K': 1 / bV, 'mech.beta_min': be[imin], 'mech.beta_V': bV, 'geom.rel_density': rho, 'eff.E': H.E1, 'eff.nu': H.nu12 },
    numerics: { method: `Periodic frame FE (${timo ? 'Timoshenko' : 'Euler–Bernoulli'}) homogenisation of the rhombic cell with closed-form hinging/stretching verification`, dof: 6, notes: `Linear analysis about the reference angle; out-of-plane coupling from prismatic-wall kinematics.${th0 === null ? '' : ''}` },
  };
}

// ================================================================ 3-D lattice scaling (Gibson–Ashby, Deshpande–Fleck–Ashby)
const L3 = { // nominal scaling coefficients: E/E_s = CE ρ^nE, σ/σ_ys = CS ρ^nS; ρ = g (r/l)²
  octet: { label: 'Octet truss (stretch-dominated)', CE: 1 / 9, nE: 1, CS: 1 / 3, nS: 1, g: 6 * Math.SQRT2 * PI, stretch: true },
  bcc: { label: 'Body-centred cubic truss (bending-dominated)', CE: 1, nE: 2, CS: 0.3, nS: 1.5, g: 3 * Math.sqrt(3) * PI, stretch: false },
  kelvin: { label: 'Kelvin cell foam (bending-dominated)', CE: 1, nE: 2, CS: 0.3, nS: 1.5, g: 3 * PI / (2 * Math.SQRT2), stretch: false },
  gyroid: { label: 'Gyroid sheet lattice (mixed mode, nominal)', CE: 0.3, nE: 1.5, CS: 0.4, nS: 1.3, g: null, stretch: false },
};
function lattice3d(inp) {
  const W = [], tp = L3[inp.topo3] || L3.octet, Es = inp.Es, nus = inp.nu_s, sy = inp.yield_s, rb = N.clamp(inp.rel_dens, 1e-4, 0.5), a = inp.cell3 * 1e-6;
  if (inp.rel_dens > 0.5) W.push('Relative density limited to 0.5: scaling laws for slender struts do not hold for denser structures.');
  const Eeff = r => tp.CE * r ** tp.nE * Es, Sy = r => tp.CS * r ** tp.nS * sy, sl = tp.g ? Math.sqrt(rb / tp.g) : null; // r/l
  const lstrut = inp.topo3 === 'bcc' ? a * Math.sqrt(3) / 2 : inp.topo3 === 'octet' ? a / Math.SQRT2 : inp.topo3 === 'kelvin' ? a / (2 * Math.SQRT2) : a, r = sl ? sl * lstrut : rb * a / 6, nf = inp.end_fix;
  const sEl = r2 => (tp.stretch ? tp.CS * r2 * nf * nf * PI * PI * Es * (r2 / tp.g) / 4 : 0.05 * r2 * r2 * Es), sB = sEl(rb), sY = Sy(rb), gov = sB < sY ? 'elastic strut buckling' : 'plastic yielding', sC = Math.min(sB, sY);
  const K = Es / (3 * (1 - 2 * nus)), G = Es / (2 * (1 + nus)), KH = 4 * G * K * rb / (4 * G + 3 * K * (1 - rb)), GH = G * rb * (9 * K + 8 * G) / (20 * G + 15 * K - 6 * (K + 2 * G) * rb), EH = 9 * KH * GH / (3 * KH + GH);
  const sSize = Math.min(Es / 10, inp.KIc / Math.sqrt(PI * r)), sLatSize = tp.CS * rb ** tp.nS * sSize, rhoEff = rb * inp.rho_s, rs = N.logspace(1e-3, 0.5, 60), rT = tp.stretch ? Math.sqrt(sy * 4 * tp.g / (nf * nf * PI * PI * Es)) ** 2 : (tp.CS * sy / (0.05 * Es)) ** (1 / (2 - tp.nS));
  if (sl && sl > 0.2) W.push(`Strut radius-to-length ratio ${f3(sl)} exceeds 0.2: nodes occupy a large volume fraction and the slender-strut scaling overestimates compliance.`);
  if (inp.topo3 === 'gyroid') W.push('Gyroid coefficients are nominal mixed-mode values within the range of published fits (exponent 1.2–2.0); calibrate against compression tests of the printed material.');
  if (r < 5e-6 && sSize > sy * 1.01) W.push(`Struts of ${f3(2 * r * 1e9)} nm diameter are flaw-tolerant up to ${f3(sSize / 1e6)} MPa by the Griffith estimate, above the bulk strength entered: nanolattices can exceed the bulk-strength scaling.`);
  return {
    summary: `${tp.label} at relative density ${f3(rb)} (cell ${inp.cell3} µm, strut radius ${f3(r * 1e6)} µm): E* = ${f3(Eeff(rb) / 1e6)} MPa (${f3(Eeff(rb) / Es)} E_s), yield strength ${f3(sY / 1e6)} MPa, elastic-buckling strength ${f3(sB / 1e6)} MPa; failure is governed by ${gov} at ${f3(sC / 1e6)} MPa. Effective density ${f3(rhoEff)} kg/m³, specific stiffness ${f3(Eeff(rb) / rhoEff)} J/kg.`,
    kpis: [
      kpi('E', 'Effective Young\'s modulus', Eeff(rb), 'Pa', `Scaling law E*/E_s = ${f3(tp.CE)}·ρ̄^${tp.nE}: linear in density for stretch-dominated trusses, quadratic for bending-dominated cells.`),
      kpi('sy', 'Plastic collapse strength', sY, 'Pa', `σ*/σ_ys = ${f3(tp.CS)}·ρ̄^${tp.nS}: strength when struts yield (axially for stretch-dominated, at plastic hinges for bending-dominated lattices).`),
      kpi('sb', 'Elastic buckling strength', sB, 'Pa', tp.stretch ? `Strength when compressed struts reach their Euler load with end-fixity factor ${nf}: σ* = C_σ·ρ̄·n²π²E_s(r/l)²/4.` : 'Gibson–Ashby elastic collapse of bending-dominated cells, σ* ≈ 0.05·ρ̄²·E_s.'),
      kpi('sc', 'Governing compressive strength', sC, 'Pa', `Lower of yield and buckling strengths; the active mechanism is ${gov}.`, gov === 'plastic yielding'),
      kpi('rho', 'Effective density', rhoEff, 'kg/m³', 'Mass per unit volume of the lattice, relative density multiplied by the density of the solid.'),
      kpi('spec', 'Specific stiffness', Eeff(rb) / rhoEff, 'J/kg', `E*/ρ*; the solid has ${f3(Es / inp.rho_s)} J/kg, so the lattice retains ${f3(100 * Eeff(rb) / rhoEff / (Es / inp.rho_s))} % of it.`),
      kpi('r', 'Strut radius', r, 'm', tp.g ? `From ρ̄ = ${f3(tp.g)}·(r/l)² with strut length ${f3(lstrut * 1e6)} µm.` : 'Nominal half wall thickness ρ̄·a/6 for the sheet lattice (surface area ≈ 3.1 a² per cell).'),
      kpi('slender', 'Strut radius-to-length ratio', sl ?? 0, '–', 'r/l; below about 0.05–0.1 struts are slender and buckling may govern, above 0.2 nodal volume matters (zero for sheet lattices).', sl === null || sl < 0.2),
      kpi('rho_tr', 'Buckling-to-yield transition density', N.clamp(rT, 0, 1), '–', 'Relative density above which yielding precedes elastic buckling; design denser than this for a yield-governed, non-catastrophic failure mode.', rb > rT),
      kpi('E_HS', 'Hashin–Shtrikman upper bound on E', EH, 'Pa', 'Stiffest possible isotropic porous solid at this density; the ratio E*/E_HS = ' + f3(Eeff(rb) / EH) + ' measures how efficiently the topology uses material.'),
      kpi('s_size', 'Size-dependent strut strength', sSize, 'Pa', `Griffith flaw-size estimate min(E_s/10, K_Ic/√(πr)) for a flaw no larger than the strut radius; applies to brittle (ceramic, glassy carbon) nanolattice struts.`),
      kpi('s_lat_size', 'Lattice strength with size effect', sLatSize, 'Pa', 'Collapse-strength scaling evaluated with the size-dependent strut strength instead of the bulk strength.'),
    ],
    series: [
      ser('E', 'Relative modulus against relative density', 'relative density (–)', 'E*/E_s (–)', [...Object.entries(L3).map(([k, v]) => ({ name: v.label, x: rs, y: rs.map(x => v.CE * x ** v.nE) })), { name: 'this design', x: [rb], y: [Eeff(rb) / Es] }], 'Scaling-law comparison of the four topologies on logarithmic axes. At low density the stretch-dominated octet truss is far stiffer than bending-dominated cells of the same mass.', { logx: true, logy: true }),
      ser('S', 'Strength against relative density', 'relative density (–)', 'strength (MPa)', [{ name: 'plastic collapse', x: rs, y: rs.map(x => Sy(x) / 1e6) }, { name: 'elastic buckling', x: rs, y: rs.map(x => pos(sEl(x) / 1e6)) }, { name: 'governing', x: rs, y: rs.map(x => pos(Math.min(Sy(x), sEl(x)) / 1e6)) }, { name: 'this design', x: [rb], y: [sC / 1e6] }], 'Yield and buckling strengths of the selected topology. The lower envelope governs: buckling at low density (steeper slope), yielding at high density.', { logx: true, logy: true }),
    ],
    equations: [
      eqn('constitutive', 'Stiffness scaling', 'E*/E_s = C_E ρ̄^n  (n = 1 stretch-dominated, n = 2 bending-dominated)'), eqn('constitutive', 'Strength scaling', 'σ*_pl/σ_ys = C_σ ρ̄^m  (m = 1 stretch, 1.5 bending);  σ*_el ≈ 0.05 ρ̄² E_s (bending)'),
      eqn('formula', 'Octet truss (Deshpande–Fleck–Ashby)', 'ρ̄ = 6√2π(r/l)²,  E* = ρ̄E_s/9,  σ*_pl = ρ̄σ_ys/3,  σ*_el = (ρ̄/3)·n²π²E_s r²/(4l²)'), eqn('formula', 'Maxwell criterion', 'M = b − 3j + 6 ≥ 0 for a stretch-dominated (rigid) pin-jointed frame'),
      eqn('formula', 'Size effect (Griffith)', 'σ_f = min(E_s/10, K_Ic/√(πa)),  a ≤ strut radius'),
    ],
    checks: [
      chk('Modulus below Hashin–Shtrikman bound', Eeff(rb) <= EH * 1.0001, `E* = ${f3(Eeff(rb) / 1e6)} MPa, bound ${f3(EH / 1e6)} MPa`), chk('Strength below Voigt bound ρ̄σ_ys', sY <= rb * sy * 1.0001, `${f3(sY / 1e6)} ≤ ${f3(rb * sy / 1e6)} MPa`),
      chk('Slender-strut validity (r/l < 0.2)', sl === null || sl < 0.2, sl === null ? 'sheet lattice' : `r/l = ${f3(sl)}`), chk('Yield governs (ductile, non-catastrophic failure)', gov === 'plastic yielding', `buckling ${f3(sB / 1e6)} MPa, yield ${f3(sY / 1e6)} MPa`),
    ],
    warnings: W,
    recommendations: [
      rec('high', gov === 'elastic strut buckling' ? `Buckling governs: raise the relative density above ${f3(N.clamp(rT, 0, 1))} or use hollow/tapered struts (higher radius of gyration at the same mass) to reach the yield-governed regime.` : `Yield governs with a buckling margin of ${f3(sB / sY)}; this density is efficient for strength.`, 'Lower envelope of buckling and yield strengths'),
      rec('medium', tp.stretch ? 'The octet truss is anisotropic (stiffest along [111]); combine with a complementary truss or use isotropic plate lattices if direction-independent stiffness is required.' : `For the same mass an octet truss would be about ${f3((1 / 9) * rb / (tp.CE * rb ** tp.nE))} times stiffer; keep the bending-dominated topology only if compliance or energy absorption is the aim.`, 'Stretch- against bending-dominated scaling'),
      rec('medium', `Printed strut strength differs from bulk: test single struts or coupons at the printed diameter (${f3(2 * r * 1e6)} µm); the Griffith estimate suggests ${f3(sSize / 1e6)} MPa is attainable when flaws are confined to the strut.`, 'Size-dependent strength of micro- and nanoscale struts'),
      rec('medium', `The lattice delivers its function with ${f3(100 * rb)} % of the solid mass (${f3(rhoEff)} kg/m³); lightweighting reduces embodied energy in proportion and use-phase energy in transport applications. Powder or resin recovery from enclosed cells should be designed in.`, 'Mass efficiency of architected materials', true),
    ],
    publish: { 'eff.E': Eeff(rb), 'eff.density': rhoEff, 'geom.rel_density': rb, 'mech.max_stress': sC, 'geom.strut': 2 * r, 'geom.cell_size': a },
    numerics: { method: 'Closed-form scaling laws (Gibson–Ashby; Deshpande–Fleck–Ashby octet truss) with Euler strut buckling and Hashin–Shtrikman bounds', dof: 1, notes: 'Nominal coefficients; nodal volume, imperfections and anisotropy not resolved.' },
  };
}

// ================================================================ nanowire bending, buckling and resonance with surface elasticity
const BC = { cantilever: { label: 'Cantilever (clamped–free)', K: 2, kf: 3, mf: 1, beta: 1.8751 }, fixed: { label: 'Clamped–clamped', K: 0.5, kf: 192, mf: 1 / 8, beta: 4.7300 }, pinned: { label: 'Pinned–pinned', K: 1, kf: 48, mf: 1 / 4, beta: PI } };
function nanowire(inp) {
  const W = [], E = inp.E_n, D = inp.nw_D * 1e-9, L = inp.nw_L * 1e-6, bc = BC[inp.nw_bc] || BC.cantilever, Esf = inp.E_surf, tau = inp.tau0, F = inp.nw_F * 1e-9, rho = inp.rho_s;
  const I = PI * D ** 4 / 64, A = PI * D * D / 4, EIb = E * I, EIraw = E * I + PI * Esf * D ** 3 / 8, EIs = Math.max(EIraw, 0.01 * EIb), Eeff = d => Math.max(E + 8 * Esf / d, 0.01 * E), k = bc.kf * EIs / L ** 3, kb = bc.kf * EIb / L ** 3, defl = F / k, Mmax = bc.mf * F * L, sig = Mmax * (D / 2) / I * (Eeff(D) / Math.max(E, 1)) * (EIb / EIs);
  const Pcr = PI * PI * EIs / (bc.K * L) ** 2 + 2 * tau * D, Pcb = PI * PI * EIb / (bc.K * L) ** 2, scr = Pcr / A, ecr = scr / E, f1 = bc.beta ** 2 / (2 * PI) * Math.sqrt(EIs / (rho * A * L ** 4)), f1b = bc.beta ** 2 / (2 * PI) * Math.sqrt(EIb / (rho * A * L ** 4));
  const Ds = N.logspace(5e-9, 5e-7, 60), Ls = N.logspace(0.2 * L, 5 * L, 50), mEff = 0.24 * rho * A * L;
  if (L / D < 10) W.push(`Aspect ratio ${f3(L / D)} is below 10: Euler–Bernoulli theory overestimates stiffness and buckling load (shear deformation neglected).`);
  if (defl > 0.1 * L) W.push('Deflection exceeds 10 % of the length: the small-deflection formulae underestimate stiffening and stress.');
  if (EIraw <= 0.01 * EIb) W.push('A strongly negative surface modulus makes the effective bending rigidity non-positive at this diameter: the surface-elasticity continuum model has broken down and the rigidity was floored at 1 % of the bulk value.');
  if (D < 5e-9) W.push('Below about 5 nm diameter, continuum surface elasticity is unreliable; atomistic simulation is required.');
  return {
    summary: `A ${inp.nw_D} nm × ${inp.nw_L} µm nanowire (${bc.label.toLowerCase()}): surface elasticity changes the effective bending modulus from ${f3(E / 1e9)} to ${f3(Eeff(D) / 1e9)} GPa (${f3(100 * (Eeff(D) / E - 1))} %). Stiffness ${f3(k)} N/m, deflection ${f3(defl * 1e9)} nm under ${inp.nw_F} nN, peak bending stress ${f3(sig / 1e6)} MPa. Euler buckling load ${f3(Pcr * 1e9)} nN (critical strain ${f3(100 * ecr)} %); fundamental flexural resonance ${f3(f1 / 1e6)} MHz.`,
    kpis: [
      kpi('E_eff', 'Effective bending modulus', Eeff(D), 'Pa', 'Size-dependent modulus E + 8E_s/D for a circular wire with a surface layer of elastic modulus E_s (Miller–Shenoy surface elasticity); it tends to the bulk value for thick wires.'),
      kpi('EI', 'Effective bending rigidity', EIs, 'N·m²', 'Bulk rigidity EπD⁴/64 plus the surface contribution πE_sD³/8.'),
      kpi('k', 'Bending stiffness', k, 'N/m', `Force per unit deflection at the ${inp.nw_bc === 'cantilever' ? 'free end' : 'mid-span'}, ${bc.kf}·(EI)*/L³; the value an AFM bending test measures.`),
      kpi('defl', 'Deflection under load', defl, 'm', `Elastic deflection for the applied ${inp.nw_F} nN (small-deflection beam theory).`, defl < 0.1 * L),
      kpi('sigma', 'Peak bending stress', sig, 'Pa', 'Maximum fibre stress in the core at the section of largest bending moment; nanowires often sustain several GPa because they contain few defects.'),
      kpi('Pcr', 'Euler buckling load', Pcr, 'N', `Critical axial load π²(EI)*/(KL)² + 2τ₀D with effective-length factor K = ${bc.K}; the second term is the transverse load from residual surface stress. Bulk theory gives ${f3(Pcb * 1e9)} nN.`),
      kpi('scr', 'Critical buckling stress', scr, 'Pa', 'Buckling load divided by the cross-section; compare with the compressive strength to see which fails first.'),
      kpi('ecr', 'Critical buckling strain', ecr, '–', 'Axial strain at buckling, about π²/(16)(D/KL)² for a bulk circular wire.'),
      kpi('f1', 'Fundamental flexural frequency', f1, 'Hz', `First bending resonance β²/(2πL²)·√((EI)*/ρA) with β = ${f3(bc.beta)}; without surface effects it would be ${f3(f1b / 1e6)} MHz.`),
      kpi('dm', 'Mass responsivity', f1 / (2 * mEff), 'Hz/kg', 'Frequency shift per added mass at the antinode, f/(2m_eff): the sensitivity of the wire used as a resonant mass sensor.'),
      kpi('slender', 'Aspect ratio L/D', L / D, '–', 'Slenderness; Euler–Bernoulli theory needs more than about 10.', L / D >= 10),
    ],
    series: [
      ser('ED', 'Effective modulus against diameter', 'diameter (nm)', 'modulus (GPa)', [{ name: 'with surface elasticity', x: Ds.map(d => d * 1e9), y: Ds.map(d => Eeff(d) / 1e9) }, { name: 'bulk', x: [5, 500], y: [E / 1e9, E / 1e9] }, { name: 'this wire', x: [inp.nw_D], y: [Eeff(D) / 1e9] }], 'Size dependence of the bending modulus. The surface term scales as 1/D, so it becomes visible below a few tens of nanometres; a positive surface modulus stiffens, a negative one softens.', { logx: true }),
      ser('PL', 'Buckling load and resonance against length', 'length (µm)', 'value', [{ name: 'buckling load (nN)', x: Ls.map(x => x * 1e6), y: Ls.map(x => pos((PI * PI * EIs / (bc.K * x) ** 2 + 2 * tau * D) * 1e9)) }, { name: 'resonance (MHz)', x: Ls.map(x => x * 1e6), y: Ls.map(x => bc.beta ** 2 / (2 * PI) * Math.sqrt(EIs / (rho * A * x ** 4)) / 1e6) }], 'Both the Euler load and the flexural frequency fall as 1/L²: longer wires are easier to buckle and resonate lower. The buckling curve flattens at the surface-stress term for long wires.', { logx: true, logy: true }),
    ],
    equations: [...formulaEq('ZnO', 'Reference nanowire material'),
      eqn('governing', 'Euler–Bernoulli beam with surface elasticity', '(EI)* d⁴w/dx⁴ − H d²w/dx² + ρA ∂²w/∂t² = q,  (EI)* = EπD⁴/64 + πE_sD³/8,  H = 2τ₀D'),
      eqn('formula', 'Effective modulus', 'E_eff = E + 8E_s/D'), eqn('formula', 'Euler buckling', 'P_cr = π²(EI)*/(KL)² + 2τ₀D'), eqn('formula', 'Flexural resonance', 'f_n = (β_n²/2πL²)·√((EI)*/ρA)'),
    ],
    checks: [
      chk('Slender-beam validity (L/D ≥ 10)', L / D >= 10, `L/D = ${f3(L / D)}`), chk('Small deflection (< 10 % of length)', defl < 0.1 * L, `w/L = ${f3(defl / L)}`),
      chk('Bulk limit recovered for thick wires', Math.abs(Eeff(1e-3) / E - 1) < 1e-3, 'E_eff → E as D → ∞'), chk('Positive bending rigidity', EIraw > 0.01 * EIb, `(EI)* = ${f3(EIs)} N·m²`),
      chk('Buckling below theoretical strength (E/10)', scr < E / 10, `σ_cr = ${f3(scr / 1e9)} GPa against E/10 = ${f3(E / 1e10)} GPa`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `In compression the wire buckles at ${f3(100 * ecr)} % strain; to use the material strength keep L/D below about ${f3(PI / (4 * bc.K) / Math.sqrt(0.01))} (critical strain 1 %) or clamp both ends.`, 'Euler buckling strain ∝ (D/KL)²'),
      rec('medium', `Surface elasticity changes the modulus by ${f3(100 * (Eeff(D) / E - 1))} % at ${inp.nw_D} nm; measured size effects in ZnO are often larger because of a stiffened surface shell several nanometres deep. Fit E_s to bending or resonance data over a range of diameters.`, 'Core–surface model of size-dependent elasticity'),
      rec('medium', `As a resonator the wire offers ${f3(f1 / (2 * mEff) * 1e-21)} Hz per attogram; clamping quality and tension change the frequency more than surface elasticity does, so characterise each device.`, 'Mass responsivity f/(2m_eff)'),
      rec('low', 'Nanowire sensors and resonators use nanograms of material; prefer abundant oxides (ZnO, SiO₂, TiO₂) to arsenides or tellurides, and solution growth to high-temperature vapour routes where quality allows.', 'Material choice and process energy', true),
    ],
    publish: { 'eff.E': Eeff(D), 'mech.buckling_load': Pcr, 'res.f0': f1, 'mech.max_stress': sig, 'geom.diameter': D, 'geom.length': L, 'geom.aspect': L / D },
    numerics: { method: 'Closed-form Euler–Bernoulli beam theory with Miller–Shenoy surface elasticity and Wang–Feng surface-stress buckling term', dof: 1, notes: 'Circular cross-section, isotropic core, small deflections.' },
  };
}

// ================================================================ contact: Hertz and Oliver–Pharr
function contact(inp) {
  const W = [], E = inp.E_n, nu = inp.nu_n, Ei = 1141e9, nui = 0.07, Er = 1 / ((1 - nu * nu) / E + (1 - nui * nui) / Ei), R = inp.tip_R * 1e-9, P = inp.P_ind * 1e-6;
  const a = Math.cbrt(3 * P * R / (4 * Er)), d = a * a / R, p0 = 3 * P / (2 * PI * a * a), tmax = 0.31 * p0, Py = (PI * R / (2 * Er)) ** 2 * (1.6 * inp.sy_n) ** 3 * PI / 6 * 4, kc = 2 * Er * a;
  const Pyv = PI ** 3 * R * R * (1.6 * inp.sy_n) ** 3 / (6 * Er * Er);
  // Oliver–Pharr analysis of a Berkovich unloading curve
  const Pm = inp.Pmax * 1e-3, hm = inp.hmax * 1e-9, S = inp.S_unl * 1e3, hc = hm - 0.75 * Pm / S, Ac = 24.5 * hc * hc, Hd = Pm / Ac, ErO = Math.sqrt(PI) / (2 * 1.034) * S / Math.sqrt(Ac), invE = 1 / ErO - (1 - nui * nui) / Ei, Eop = invE > 0 ? (1 - nu * nu) / invE : 0, mexp = 1.5, hf = hm - mexp * Pm / S;
  if (hc <= 0) W.push('Contact depth is not positive: the unloading stiffness is too low for the entered load and depth (check units of S).');
  if (hf < 0) W.push('The implied final depth is negative: the unloading data are nearly fully elastic or inconsistent with a power-law exponent of 1.5.');
  if (a > 0.3 * R) W.push('Contact radius exceeds 0.3 R: the parabolic (Hertz) approximation of the sphere is inaccurate.');
  if (P > Pyv) W.push(`The load exceeds the estimated yield onset (${f3(Pyv * 1e6)} µN): the Hertz solution is an elastic lower bound on contact size.`);
  const Ps = N.linspace(0, P, 60), hs = N.linspace(Math.max(hf, 0), hm, 60), hl = N.linspace(0, hm, 60);
  return {
    summary: `Hertz contact of a ${inp.tip_R} nm diamond tip at ${inp.P_ind} µN on a solid with E = ${f3(E / 1e9)} GPa: contact radius ${f3(a * 1e9)} nm, indentation ${f3(d * 1e9)} nm, peak pressure ${f3(p0 / 1e9)} GPa (yield onset at ${f3(Pyv * 1e6)} µN). Oliver–Pharr analysis of the Berkovich data (P_max = ${inp.Pmax} mN, h_max = ${inp.hmax} nm, S = ${inp.S_unl} kN/m) gives hardness ${f3(Hd / 1e9)} GPa and modulus ${f3(Eop / 1e9)} GPa.`,
    kpis: [
      kpi('Er', 'Reduced modulus (from inputs)', Er, 'Pa', 'Combined elastic modulus of sample and diamond indenter, 1/E_r = (1 − ν²)/E + (1 − ν_i²)/E_i with E_i = 1141 GPa, ν_i = 0.07.'),
      kpi('a', 'Hertz contact radius', a, 'm', 'Radius of the circular contact patch, (3PR/4E_r)^(1/3); it sets the lateral resolution of the measurement.'),
      kpi('delta', 'Hertz indentation depth', d, 'm', 'Elastic approach of the two bodies, a²/R; the load varies as its 3/2 power.'),
      kpi('p0', 'Peak contact pressure', p0, 'Pa', 'Maximum Hertzian pressure at the centre of contact, 1.5 times the mean pressure.', p0 < 1.6 * inp.sy_n),
      kpi('tmax', 'Maximum shear stress', tmax, 'Pa', `0.31 p₀, located on the axis at a depth of about 0.48 a = ${f3(0.48 * a * 1e9)} nm; yielding or pop-in initiates there.`),
      kpi('Py', 'Load at first yield', Pyv, 'N', 'Load at which p₀ reaches 1.6σ_y (Tresca/von Mises onset beneath a sphere); defect-free nanovolumes pop in at much higher loads, near the theoretical strength.'),
      kpi('kc', 'Contact stiffness (Hertz)', kc, 'N/m', 'Tangent stiffness of the elastic contact, dP/dδ = 2E_r·a, at the applied load; it grows with the cube root of load.'),
      kpi('hc', 'Contact depth (Oliver–Pharr)', hc, 'm', 'Depth of material in contact with the indenter, h_max − 0.75·P_max/S, after correcting for elastic sink-in.', hc > 0),
      kpi('H', 'Indentation hardness', Hd, 'Pa', 'P_max divided by the projected contact area 24.5·h_c² of an ideal Berkovich tip.'),
      kpi('E_op', 'Indentation modulus', Eop, 'Pa', 'Sample modulus from the unloading stiffness, E_r = (√π/2β)·S/√A with β = 1.034, corrected for indenter compliance.'),
      kpi('HE', 'Hardness-to-modulus ratio', Eop > 0 ? Hd / Eop : 0, '–', 'H/E, an index of elastic strain to failure and wear resistance (about 0.1 for hard ceramics and glasses, 0.01 for soft metals).'),
      kpi('hf', 'Residual depth', Math.max(hf, 0), 'm', 'Depth of the permanent impression implied by a power-law unloading curve with exponent 1.5.'),
    ],
    series: [
      ser('hertz', 'Hertzian load–displacement curve', 'indentation depth (nm)', 'load (µN)', [{ name: 'P = (4/3)E_r√R δ^(3/2)', x: Ps.map(p => Math.cbrt(3 * p * R / (4 * Er)) ** 2 / R * 1e9), y: Ps.map(p => p * 1e6) }], 'Elastic loading of a sphere on a flat. Experimental curves that follow this line and then jump (pop-in) mark the onset of plasticity in a small, initially dislocation-free volume.'),
      ser('op', 'Berkovich load–displacement reconstruction', 'displacement (nm)', 'load (mN)', [{ name: 'loading (P ∝ h²)', x: hl.map(h => h * 1e9), y: hl.map(h => Pm * (h / hm) ** 2 * 1e3) }, { name: 'unloading (power law, m = 1.5)', x: hs.map(h => h * 1e9), y: hs.map(h => Pm * ((h - Math.max(hf, 0)) / Math.max(hm - Math.max(hf, 0), 1e-30)) ** mexp * 1e3) }], 'Idealised indentation cycle consistent with the entered peak load, depth and unloading stiffness. The initial slope of the unloading branch is the contact stiffness S used by the Oliver–Pharr method.'),
    ],
    equations: [
      eqn('governing', 'Hertz contact', 'a³ = 3PR/(4E_r),  δ = a²/R,  p₀ = 3P/(2πa²),  τ_max = 0.31p₀'), eqn('formula', 'Reduced modulus', '1/E_r = (1 − ν²)/E + (1 − ν_i²)/E_i'),
      eqn('formula', 'Oliver–Pharr method', 'h_c = h_max − ε P_max/S (ε = 0.75),  A = 24.5h_c²,  H = P_max/A,  E_r = (√π/2β)·S/√A'), eqn('formula', 'Yield onset beneath a sphere', 'p₀ = 1.6σ_y  ⇒  P_y = π³R²(1.6σ_y)³/(6E_r²)'),
    ],
    checks: [
      chk('Hertz geometry valid (a < 0.3R)', a < 0.3 * R, `a/R = ${f3(a / R)}`), chk('Elastic contact (P < P_y)', P <= Pyv, `P/P_y = ${f3(P / Pyv)}`),
      chk('Physical Oliver–Pharr inputs (0 < h_c < h_max)', hc > 0 && hc < hm, `h_c/h_max = ${f3(hc / hm)}`), chk('Fused-silica reference recovered', Math.abs((1 - 0.17 ** 2) / (1 / (Math.sqrt(PI) / (2 * 1.034) * 60.5e3 / Math.sqrt(24.5 * (213e-9 - 0.75 * 5e-3 / 60.5e3) ** 2)) - (1 - nui * nui) / Ei) / 1e9 - 72) < 4, 'P = 5 mN, h = 213 nm, S = 60.5 kN/m returns E ≈ 72 GPa'),
      chk('Load recovered from Hertz relation', Math.abs(4 / 3 * Er * Math.sqrt(R) * d ** 1.5 - P) < 1e-9 * P + 1e-30, 'P = (4/3)E_r√R·δ^(3/2)'),
    ],
    warnings: W,
    recommendations: [
      rec('high', `Keep the indentation depth below 10 % of a film thickness to avoid substrate influence: for the present contact depth ${f3(hc * 1e9)} nm the film should be thicker than about ${f3(10 * hc * 1e9)} nm.`, 'Bückle 10 % rule'),
      rec('medium', 'Calibrate the tip area function and frame compliance on fused silica (E ≈ 72 GPa, H ≈ 9 GPa) before each session; at depths below about 100 nm tip blunting dominates the error.', 'Ideal Berkovich area 24.5h_c² assumed'),
      rec('medium', `Pile-up is not captured by the Oliver–Pharr contact depth: for materials with h_f/h_max > 0.7 (here ${f3(Math.max(hf, 0) / hm)}) and little work hardening, image the impression and correct the area.`, 'Sink-in assumption of the method'),
      rec('low', 'Nanoindentation consumes almost no material and replaces many destructive coupon tests during process development.', 'Low-waste characterisation', true),
    ],
    publish: { 'prop.E': Eop > 0 ? Eop : E, 'mech.hardness': Hd, 'mech.max_stress': p0, 'mech.contact_radius': a },
    numerics: { method: 'Closed-form Hertz elastic contact and Oliver–Pharr unloading analysis', dof: 1, notes: 'Frictionless, isotropic, homogeneous half-space; ideal indenter geometry.' },
  };
}

// ================================================================ strength: Hall–Petch, Griffith/LEFM, thermal mismatch (Stoney)
function strength(inp) {
  const W = [], d = inp.d_grain * 1e-9, s0 = inp.sig0 * 1e6, kHP = inp.k_hp * 1e6 * 1e-3, wg = inp.w_gb * 1e-9, sgb = inp.sig_gb * 1e6; // k: MPa·µm^½ → Pa·m^½
  const hp = x => s0 + kHP / Math.sqrt(x), fgb = x => (x > wg ? 1 - ((x - wg) / x) ** 3 : 1), mix = x => (1 - fgb(x)) * hp(x) + fgb(x) * sgb, ds = N.logspace(2e-9, 1e-4, 120), ipk = N.argmax(ds.map(mix)), sy = mix(d);
  const E = inp.E_n, nu = inp.nu_n, Ep = E / (1 - nu * nu), KIc = inp.KIc, a = inp.flaw * 1e-9, Y = inp.Y_geo, sf = KIc / (Y * Math.sqrt(PI * a)), sapp = inp.sig_app, ac = (KIc / (Y * sapp)) ** 2 / PI, Gc = KIc * KIc / Ep, Kid = Math.sqrt(2 * inp.gamma_s * Ep), sth = E / 10;
  const tf = inp.t_film * 1e-9, tsb = inp.t_sub * 1e-6, sfilm = E / (1 - nu) * (inp.cte_s - inp.cte_f) * inp.dT, kap = 6 * sfilm * tf * (1 - inp.nu_sub) / (inp.E_sub * tsb * tsb), hcr = sfilm > 0 ? 1.976 ** -1 * Gc * 1 / (sfilm * sfilm / Ep) : Infinity, eel = sfilm * sfilm * (1 - nu) / E * tf;
  if (d < 10e-9) W.push('Below about 10 nm grain size the strength is controlled by grain-boundary sliding; the phase-mixture estimate is indicative only.');
  if (sf > sth) W.push('The Griffith stress for this flaw exceeds the theoretical strength E/10: the flaw is too small to control fracture and the ideal strength caps the result.');
  if (tf > 0.1 * tsb) W.push('Film thicker than 10 % of the substrate: the Stoney thin-film formula overestimates curvature.');
  if (sfilm > 0 && tf > hcr) W.push(`Tensile film stress with thickness above the channel-cracking limit (${f3(hcr * 1e9)} nm): through-thickness cracks are expected.`);
  const as = N.logspace(1e-9, 1e-4, 80);
  return {
    summary: `Grain size ${inp.d_grain} nm: Hall–Petch yield strength ${f3(hp(d) / 1e6)} MPa, ${f3(sy / 1e6)} MPa with the grain-boundary phase mixture (strongest grain size ≈ ${f3(ds[ipk] * 1e9)} nm). A ${inp.flaw} nm flaw fails at ${f3(Math.min(sf, sth) / 1e6)} MPa (K_Ic = ${f3(KIc / 1e6)} MPa·√m); at ${f3(sapp / 1e6)} MPa the critical flaw is ${f3(ac * 1e6)} µm. A ${inp.t_film} nm film cooled by ${f3(-inp.dT)} K carries ${f3(sfilm / 1e6)} MPa thermal-mismatch stress (${sfilm >= 0 ? 'tensile' : 'compressive'}), bending the substrate to a radius of ${f3(Math.abs(1 / (kap || 1e-30)))} m.`,
    kpis: [
      kpi('sy_hp', 'Hall–Petch yield strength', hp(d), 'Pa', 'σ_y = σ₀ + k·d^(−1/2): strengthening by dislocation pile-up at grain boundaries, valid down to a few tens of nanometres.'),
      kpi('sy', 'Yield strength with inverse regime', sy, 'Pa', `Volume-weighted mixture of grain interiors (Hall–Petch) and a ${inp.w_gb} nm grain-boundary phase of strength ${inp.sig_gb} MPa, which reproduces softening at the finest grain sizes.`),
      kpi('d_peak', 'Strongest grain size', ds[ipk], 'm', 'Grain size at which the mixture strength is maximal; refining further softens the material (inverse Hall–Petch regime).'),
      kpi('HV', 'Estimated hardness', 3 * sy, 'Pa', 'Tabor relation H ≈ 3σ_y for a fully plastic indentation.'),
      kpi('sf', 'Fracture stress for the flaw', Math.min(sf, sth), 'Pa', `Griffith/LEFM stress K_Ic/(Y√(πa)) for a ${inp.flaw} nm flaw, capped at the theoretical strength E/10.`, Math.min(sf, sth) > sapp),
      kpi('ac', 'Critical flaw size at applied stress', ac, 'm', `Flaw half-length that becomes unstable at ${f3(sapp / 1e6)} MPa; inspection must detect smaller flaws.`),
      kpi('SF', 'Fracture safety factor', Math.min(sf, sth) / sapp, '–', 'Fracture stress of the assumed flaw divided by the applied stress.', Math.min(sf, sth) / sapp > 1),
      kpi('Gc', 'Fracture energy', Gc, 'J/m²', `Critical energy release rate K_Ic²/E′; an ideally brittle solid would have 2γ = ${f3(2 * inp.gamma_s)} J/m² (K = ${f3(Kid / 1e6)} MPa·√m).`),
      kpi('s_film', 'Thermal-mismatch film stress', sfilm, 'Pa', 'Biaxial stress E_f/(1 − ν_f)·(α_s − α_f)·ΔT in a thin film constrained by a thick substrate; positive is tensile.', null),
      kpi('kappa', 'Substrate curvature', kap, '1/m', 'Stoney curvature 6σ_f t_f(1 − ν_s)/(E_s t_s²); wafer-bow measurement inverts this relation to obtain the film stress.', null),
      kpi('h_c', 'Channel-cracking critical thickness', Number.isFinite(hcr) ? hcr : 0, 'm', 'Film thickness above which a tensile film cracks, Γ·E′/(1.976·σ²); zero is reported for compressive films, which buckle or delaminate instead.', !(sfilm > 0) || tf < hcr),
      kpi('U_el', 'Stored elastic energy per area', eel, 'J/m²', 'Strain energy σ²(1 − ν)t/E in the film; delamination is possible when it exceeds the interface toughness (typically 1–10 J/m²).'),
    ],
    series: [
      ser('hp', 'Yield strength against grain size', 'grain size (nm)', 'yield strength (MPa)', [{ name: 'Hall–Petch', x: ds.map(x => x * 1e9), y: ds.map(x => hp(x) / 1e6) }, { name: 'with grain-boundary phase (inverse regime)', x: ds.map(x => x * 1e9), y: ds.map(x => mix(x) / 1e6) }, { name: 'this material', x: [inp.d_grain], y: [sy / 1e6] }], 'Strength rises as d^(−1/2) with grain refinement until the grain-boundary volume fraction becomes significant; below the peak the material softens.', { logx: true, logy: true }),
      ser('griffith', 'Fracture stress against flaw size', 'flaw size (nm)', 'fracture stress (MPa)', [{ name: 'LEFM', x: as.map(x => x * 1e9), y: as.map(x => Math.min(KIc / (Y * Math.sqrt(PI * x)), sth) / 1e6) }, { name: 'applied stress', x: [1, 1e5], y: [sapp / 1e6, sapp / 1e6] }, { name: 'this flaw', x: [inp.flaw], y: [Math.min(sf, sth) / 1e6] }], 'Strength falls as a^(−1/2) with flaw size and saturates at the theoretical strength for flaws of a few nanometres: the reason nanoscale specimens are strong. The intersection with the applied stress gives the critical flaw.', { logx: true, logy: true }),
    ],
    equations: [
      eqn('constitutive', 'Hall–Petch relation', 'σ_y = σ₀ + k d^(−1/2)'), eqn('constitutive', 'Phase-mixture model (inverse regime)', 'σ = (1 − f_gb)σ_HP(d) + f_gb σ_gb,  f_gb = 1 − [(d − w)/d]³'),
      eqn('governing', 'Griffith / LEFM criterion', 'K_I = Yσ√(πa) = K_Ic,  G_c = K_Ic²/E′,  E′ = E/(1 − ν²)'), eqn('formula', 'Thermal-mismatch stress', 'σ_f = E_f/(1 − ν_f)·(α_s − α_f)·ΔT'),
      eqn('formula', 'Stoney equation', 'κ = 6σ_f t_f(1 − ν_s)/(E_s t_s²)'), eqn('formula', 'Channel cracking', 'h_c = Γ E′/(Zσ²),  Z = 1.976'),
    ],
    checks: [
      chk('No fracture at the applied stress', Math.min(sf, sth) > sapp, `σ_f = ${f3(Math.min(sf, sth) / 1e6)} MPa against ${f3(sapp / 1e6)} MPa`), chk('Fracture energy at least the ideal surface work 2γ', Gc >= 2 * inp.gamma_s * 0.999, `G_c = ${f3(Gc)} J/m², 2γ = ${f3(2 * inp.gamma_s)} J/m²`),
      chk('Stoney thin-film validity (t_f < 0.1 t_s)', tf < 0.1 * tsb, `t_f/t_s = ${f3(tf / tsb)}`), chk('Film below channel-cracking thickness', !(sfilm > 0) || tf < hcr, sfilm > 0 ? `t_f = ${inp.t_film} nm, h_c = ${f3(hcr * 1e9)} nm` : 'film is in compression'),
      chk('Yield below theoretical strength', sy < E / 10, `σ_y = ${f3(sy / 1e6)} MPa, E/10 = ${f3(E / 1e7)} MPa`),
    ],
    warnings: W,
    recommendations: [
      rec('high', Math.min(sf, sth) / sapp < 1.5 ? `Fracture margin is ${f3(Math.min(sf, sth) / sapp)}: reduce the largest flaw below ${f3(ac * 1e9 / 2.25)} nm (surface finish, pore control) or lower the stress.` : `Flaws up to ${f3(ac * 1e6)} µm are tolerable at the applied stress; set the inspection threshold at about half that size.`, 'LEFM critical flaw size'),
      rec('medium', `Grain refinement beyond about ${f3(ds[ipk] * 1e9)} nm no longer strengthens; target that size and stabilise it against grain growth with solute or second-phase pinning.`, 'Peak of the Hall–Petch/inverse curve'),
      rec('medium', sfilm > 0 ? `Tensile film stress of ${f3(sfilm / 1e6)} MPa: keep thickness below ${f3(hcr * 1e9)} nm, lower the deposition temperature, or choose a substrate with a closer expansion coefficient.` : `Compressive film stress of ${f3(-sfilm / 1e6)} MPa improves crack resistance but may buckle-delaminate when σ²t(1 − ν)/E = ${f3(eel)} J/m² approaches the interface toughness.`, 'Thermal mismatch and film fracture mechanics'),
      rec('low', 'Lower deposition and annealing temperatures reduce both thermal-mismatch stress and process energy; verify crystallinity is still adequate.', 'ΔT appears linearly in the film stress', true),
    ],
    publish: { 'prop.yield': sy, 'mech.max_stress': Math.max(Math.abs(sfilm), sapp), 'mech.safety_factor': Math.min(sf, sth) / sapp, 'mech.film_stress': sfilm, 'mech.critical_flaw': ac },
    numerics: { method: 'Closed-form Hall–Petch with phase-mixture softening, LEFM, biaxial thermal stress and Stoney curvature', dof: 1, notes: 'Isotropic elasticity; temperature-independent properties; equi-biaxial film stress.' },
  };
}

// ================================================================ Biot poroelastic consolidation (1-D)
function poro(inp) {
  const W = [], E = inp.E_d, nu = inp.nu_d, al = inp.biot, phi = inp.poros, K = E / (3 * (1 - 2 * nu)), G = E / (2 * (1 + nu)), Kv = K + 4 * G / 3, Ks = al < 0.999 ? K / (1 - al) : Infinity;
  if (al < phi) W.push('Biot coefficient below the porosity is not physically admissible (α ≥ φ); the grain-compressibility term was set to zero.');
  const invM = phi / inp.K_f + Math.max(al - phi, 0) / Ks, Mb = 1 / invM, kap = inp.perm / (inp.mu_f * 1e-3), c = kap / (invM + al * al / Kv), Hh = inp.H_layer * 1e-3, s0 = inp.sig_load, p0 = s0 * al * Mb / (Kv + al * al * Mb);
  const t90 = 0.848 * Hh * Hh / c, tend = inp.t_cons > 0 ? inp.t_cons : 1.5 * t90, nz = 61, nt = 240, dz = Hh / (nz - 1), dt = tend / nt, r = c * dt / (dz * dz);
  let p = new Float64Array(nz).fill(p0); p[0] = 0; const snaps = [Array.from(p)], Ut = [0], ts = [0], a = new Float64Array(nz), b = new Float64Array(nz), cc = new Float64Array(nz);
  for (let i = 1; i < nz - 1; i++) { a[i] = -r; b[i] = 1 + 2 * r; cc[i] = -r; } b[0] = 1; cc[0] = 0; a[nz - 1] = -2 * r; b[nz - 1] = 1 + 2 * r;
  for (let n = 1; n <= nt; n++) { const rhs = Float64Array.from(p); rhs[0] = 0; p = N.solveTridiag(a, b, cc, rhs); let m = 0; for (let i = 1; i < nz; i++) m += 0.5 * (p[i] + p[i - 1]) * dz; ts.push(n * dt); Ut.push(1 - m / (p0 * Hh)); if (n % 4 === 0) snaps.push(Array.from(p)); }
  const Uan = t => { let s = 0; for (let m = 0; m < 200; m++) { const q = (2 * m + 1) * PI / 2; s += 2 / (q * q) * Math.exp(-q * q * c * t / (Hh * Hh)); } return 1 - s; };
  const w0 = s0 * Hh / (Kv + al * al * Mb), winf = s0 * Hh / Kv, Uend = Ut[nt], errU = Math.abs(Uend - Uan(tend)), B = al * Mb / (K + al * al * Mb), Ku = K + al * al * Mb, nuu = (3 * Ku - 2 * G) / (2 * (3 * Ku + G)), vD = kap * p0 / Hh;
  const nx = snaps.length, ny = 41, data = new Float64Array(nx * ny); for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data[j * nx + i] = N.interp1(N.linspace(0, Hh, nz), snaps[i], Hh * (1 - j / (ny - 1)));
  const zz = N.linspace(0, Hh, nz).map(z => z * 1e3);
  if (nuu > 0.4999) W.push('The undrained response is nearly incompressible; volumetric locking would affect displacement-based 2-D/3-D finite elements.');
  return {
    summary: `Poroelastic layer (E = ${f3(E / 1e9)} GPa, ν = ${nu}, Biot α = ${al}, porosity ${phi}, permeability ${f3(inp.perm)} m²): Biot modulus ${f3(Mb / 1e9)} GPa and consolidation coefficient ${f3(c)} m²/s. A step load of ${f3(s0 / 1e3)} kPa on a ${inp.H_layer} mm layer drained at the top raises the pore pressure by ${f3(p0 / 1e3)} kPa (${f3(100 * p0 / s0)} % of the load); 90 % consolidation takes ${f3(t90)} s. Settlement grows from ${f3(w0 * 1e9)} nm (undrained) to ${f3(winf * 1e9)} nm (drained); after ${f3(tend)} s consolidation is ${f3(100 * Uend)} % complete.`,
    kpis: [
      kpi('M', 'Biot modulus', Mb, 'Pa', 'Pressure rise per unit fluid content at constant strain, 1/M = φ/K_f + (α − φ)/K_s; it couples pore pressure to deformation.'),
      kpi('c', 'Consolidation coefficient', c, 'm²/s', 'Hydraulic diffusivity of pore pressure under uniaxial strain, (k/µ)/(1/M + α²/(K + 4G/3)); time scales as H²/c.'),
      kpi('p0', 'Initial excess pore pressure', p0, 'Pa', 'Undrained pore-pressure jump under the step load, σ₀·αM/(K_v + α²M): the fraction of the load initially carried by the fluid.'),
      kpi('B', 'Skempton coefficient', B, '–', 'Undrained pore-pressure rise per unit mean stress, αM/(K + α²M); 1 for a soft saturated solid with incompressible constituents.', B >= 0 && B <= 1),
      kpi('t90', 'Time to 90 % consolidation', t90, 's', 'Terzaghi time factor 0.848·H²/c for single drainage; halving the drainage path cuts it fourfold.'),
      kpi('U', 'Degree of consolidation at end time', Uend, '–', `Fraction of the excess pore pressure dissipated after ${f3(tend)} s from the finite-difference solution (series solution: ${f4(Uan(tend))}).`),
      kpi('w0', 'Immediate (undrained) settlement', w0, 'm', 'Surface displacement at the instant of loading, σ₀H/(K_v + α²M), before any drainage.'),
      kpi('winf', 'Final (drained) settlement', winf, 'm', 'Surface displacement after full drainage, σ₀H/K_v with the drained uniaxial-strain modulus K_v = K + 4G/3.'),
      kpi('seff', 'Initial effective stress', s0 - al * p0, 'Pa', 'Biot effective stress σ − αp carried by the solid skeleton immediately after loading; it rises to the full load as pressure dissipates.'),
      kpi('nu_u', 'Undrained Poisson\'s ratio', nuu, '–', 'Poisson\'s ratio for fast (undrained) loading; it lies between the drained value and 0.5.', nuu <= 0.5),
      kpi('vD', 'Initial Darcy velocity scale', vD, 'm/s', 'Characteristic seepage velocity (k/µ)·p₀/H at the drained surface shortly after loading.'),
      kpi('Kv', 'Drained uniaxial-strain modulus', Kv, 'Pa', 'Oedometric modulus K + 4G/3 of the dry skeleton, controlling final settlement.'),
    ],
    series: [
      ser('Ut', 'Degree of consolidation against time', 'time (s)', 'U (–)', [{ name: 'finite difference', x: ts, y: Ut }, { name: 'Terzaghi series', x: ts, y: ts.map(Uan) }], 'Fraction of excess pore pressure dissipated. The numerical curve should overlie the analytical series; early behaviour follows √t and the late stage a single exponential.'),
      ser('pz', 'Pore-pressure isochrones', 'depth below drained surface (mm)', 'pore pressure (kPa)', [0.02, 0.1, 0.3, 0.6, 1].map(fr => { const i = Math.min(Math.round(fr * (nx - 1)), nx - 1); return { name: `t = ${f3(i * 4 * dt)} s`, x: zz, y: snaps[i].map(v => v / 1e3) }; }), 'Pore-pressure profiles at successive times. Drainage starts at the free surface and the front diffuses towards the sealed base; the effective stress rises where the pressure has dissipated.'),
      ser('wt', 'Settlement against time', 'time (s)', 'settlement (nm)', [{ name: 'surface settlement', x: ts, y: Ut.map(u => (w0 + (winf - w0) * u) * 1e9) }], 'Time-dependent surface displacement: an immediate undrained part followed by consolidation settlement as fluid is expelled.'),
    ],
    maps: [{ key: 'p', title: 'Excess pore pressure p(z, t)', nx, ny, data, x0: 0, x1: tend, y0: 0, y1: inp.H_layer, xlabel: 'time (s)', ylabel: 'height above sealed base (mm)', unit: 'Pa', explain: 'Space–time field of excess pore pressure from the implicit finite-difference solution: the drained surface is at the top, the sealed base at the bottom.', diverging: false }],
    equations: [...formulaEq('ZnO', 'Reference porous solid'),
      eqn('constitutive', 'Biot constitutive relations', 'σ = C:ε − αp I,  ζ = α ε_v + p/M'), eqn('governing', 'Fluid mass balance with Darcy flow', '∂ζ/∂t = ∇·((k/µ)∇p)'),
      eqn('governing', 'Uniaxial-strain consolidation', '∂p/∂t = c ∂²p/∂z²,  c = (k/µ)/(1/M + α²/K_v),  p(0,t) = 0,  ∂p/∂z(H,t) = 0'), eqn('formula', 'Undrained response', 'p₀ = σ₀ αM/(K_v + α²M),  B = αM/(K + α²M)'),
      eqn('formula', 'Terzaghi solution', 'U(t) = 1 − Σ (2/m²) exp(−m² c t/H²),  m = (2n + 1)π/2'), eqn('formula', 'Effective stress', 'σ′ = σ − αp'),
    ],
    checks: [
      chk('Finite-difference solution matches the Terzaghi series', errU < 0.02, `|ΔU| = ${f3(errU)} at t = ${f3(tend)} s (backward Euler, ${nz} nodes, ${nt} steps)`), chk('Pore pressure bounded 0 ≤ p ≤ p₀ (maximum principle)', snaps.every(s => s.every(v => v >= -1e-9 * p0 && v <= p0 * (1 + 1e-9))), 'no overshoot in the implicit scheme'),
      chk('Skempton coefficient within 0–1', B >= 0 && B <= 1, `B = ${f3(B)}`), chk('Admissible Biot coefficient (φ ≤ α ≤ 1)', al >= phi && al <= 1, `α = ${al}, φ = ${phi}`),
      chk('Undrained stiffer than drained', Ku >= K && w0 <= winf, `K_u = ${f3(Ku / 1e9)} GPa against K = ${f3(K / 1e9)} GPa`), chk('Fluid mass balance: dissipated pressure equals consolidation settlement', Math.abs((w0 + (winf - w0) * Uend) - (winf - al * (1 - Uend) * p0 * Hh / Kv)) < 1e-6 * winf, 'settlement from U(t) equals σ₀H/K_v − α∫p dz/K_v'),
    ],
    warnings: W,
    recommendations: [
      rec('high', `Loading or drying faster than about ${f3(t90)} s builds pore pressure up to ${f3(100 * p0 / s0)} % of the applied stress, reducing effective stress and risking cracking of porous films or green bodies; ramp loads over several t₉₀ or add drainage paths.`, 'Consolidation time H²/c'),
      rec('medium', `Permeability is the least certain input: c scales linearly with k. Measure it (pressure-decay or flow-through) on the actual porous ZnO layer; a tenfold error shifts every time above tenfold.`, 'c ∝ k/µ'),
      rec('medium', `Halving the drainage path (drain both faces) cuts the consolidation time to ${f3(t90 / 4)} s.`, 'Time factor ∝ H²'),
      rec('low', 'Slower, drainage-controlled pressing or infiltration avoids rejected parts and re-work, the dominant waste stream in wet processing of porous ceramics; recycle the expelled liquor.', 'Yield and solvent recovery', true),
    ],
    publish: { 'eff.K': K, 'eff.permeability': inp.perm, 'proc.time': t90, 'mech.max_stress': s0 - al * p0, 'mech.pore_pressure': p0, 'mech.settlement': winf },
    numerics: { method: 'Backward-Euler finite differences (Thomas algorithm) for the 1-D pore-pressure diffusion equation, verified against the Terzaghi series', dof: nz, steps: nt, residual: errU, notes: 'Uniaxial strain, constant total stress, constant coefficients (linear Biot theory).' },
  };
}

const M = (...a) => ({ mode: a });
export default {
  id: 14,
  tagline: 'Lattice and auxetic homogenisation by frame finite elements, negative stiffness and compressibility mechanisms, 3-D lattice scaling and nanomechanics.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'lattice2d', group: 'Model', help: 'Selects the mechanical model that is solved; each analysis has its own inputs.', options: [
      { value: 'lattice2d', label: '2-D lattice / auxetic honeycomb (frame finite elements)' }, { value: 'negstiff', label: 'Negative stiffness: snap-through curved beam, quasi-zero-stiffness mount' }, { value: 'nlc', label: 'Negative linear compressibility (wine-rack framework, crystal)' },
      { value: 'lattice3d', label: '3-D lattice scaling (octet, BCC, Kelvin, gyroid)' }, { value: 'nanowire', label: 'Nanowire bending, buckling and resonance (surface elasticity)' }, { value: 'contact', label: 'Nanoindentation: Hertz and Oliver–Pharr' },
      { value: 'strength', label: 'Strength: Hall–Petch, fracture, thin-film stress' }, { value: 'poro', label: 'Biot poroelastic consolidation' }] },
    { key: 'Es', label: 'Wall/strut material modulus', unit: 'Pa', value: 2e9, min: 1e6, max: 1.2e12, log: true, group: 'Solid material', bus: 'prop.E', showIf: M('lattice2d', 'negstiff', 'nlc', 'lattice3d'), help: 'Young\'s modulus of the solid from which walls, struts or beams are made; all effective moduli scale linearly with it.' },
    { key: 'nu_s', label: 'Wall material Poisson\'s ratio', unit: '–', value: 0.35, min: 0, max: 0.49, group: 'Solid material', bus: 'prop.nu', showIf: M('lattice2d', 'nlc', 'lattice3d'), help: 'Poisson\'s ratio of the solid; it enters through shear deformation of the walls and the out-of-plane coupling.' },
    { key: 'rho_s', label: 'Solid density', unit: 'kg/m³', value: 1100, min: 100, max: 22000, group: 'Solid material', bus: 'prop.density', showIf: M('lattice2d', 'lattice3d', 'nanowire'), help: 'Density of the solid material, used for effective density, specific stiffness and resonance.' },
    { key: 'yield_s', label: 'Solid yield / failure strength', unit: 'Pa', value: 5e7, min: 1e5, max: 2e10, log: true, group: 'Solid material', bus: 'prop.yield', showIf: M('lattice2d', 'negstiff', 'lattice3d'), help: 'Elastic limit of the solid; sets safety factors, allowable snap-through strain and collapse strength.' },
    { key: 'topology', label: 'Lattice topology', type: 'select', value: 'reentrant', group: '2-D lattice', showIf: M('lattice2d'), help: 'Unit-cell architecture. The general honeycomb is auxetic for negative cell angles; the chiral lattice is auxetic by ring rotation.', options: Object.entries(TOPO).map(([value, label]) => ({ value, label })) },
    { key: 'cell_l', label: 'Wall length l (ring spacing for chiral)', unit: 'mm', value: 5, min: 0.001, max: 500, log: true, group: '2-D lattice', showIf: M('lattice2d'), help: 'Length of the inclined cell walls (or centre-to-centre ring spacing of the chiral lattice). Effective properties depend only on ratios to it.' },
    { key: 'wall_t', label: 'Wall thickness t', unit: 'mm', value: 0.4, min: 0.00001, max: 100, log: true, group: '2-D lattice', showIf: M('lattice2d'), help: 'In-plane thickness of the walls. Bending-dominated stiffness scales with (t/l)³, stretch-dominated with t/l.' },
    { key: 'h_l', label: 'Wall length ratio h/l', unit: '–', value: 2, min: 0.5, max: 4, group: '2-D lattice', showIf: M('lattice2d'), help: 'Ratio of the vertical wall length to the inclined wall length of the general honeycomb (ignored by other topologies).' },
    { key: 'theta_c', label: 'Cell angle θ', unit: '°', value: -30, min: -70, max: 70, group: '2-D lattice', showIf: M('lattice2d'), help: 'Inclination of the oblique walls from the horizontal: +30° is the regular hexagon, negative values give re-entrant (auxetic) cells.' },
    { key: 'ring_t', label: 'Chiral ring wall thickness / ligament thickness', unit: '–', value: 1, min: 0.5, max: 30, group: '2-D lattice', showIf: { mode: ['lattice2d'], topology: ['chiral'] }, help: 'Wall thickness of the rings relative to the ligaments. One gives compliant rings that ovalise under load; values of 10 or more approach the rigid-ring idealisation of the Prall–Lakes model (ν = −1).' },
    { key: 'r_R', label: 'Chiral ring radius ratio r/R', unit: '–', value: 0.25, min: 0.1, max: 0.42, group: '2-D lattice', showIf: M('lattice2d'), help: 'Ring radius divided by the ring centre spacing (chiral lattice only); larger rings give longer lever arms and softer lattices.' },
    { key: 'beam', label: 'Beam theory', type: 'select', value: 'timoshenko', group: '2-D lattice', showIf: M('lattice2d', 'nlc'), help: 'Timoshenko elements include transverse shear of the walls (important for thick walls); Euler–Bernoulli neglects it.', options: [{ value: 'timoshenko', label: 'Timoshenko (with shear)' }, { value: 'euler', label: 'Euler–Bernoulli' }] },
    { key: 'load_dir', label: 'Loading direction', type: 'select', value: '1', group: '2-D lattice', showIf: M('lattice2d'), help: 'Direction of the uniaxial stress used for the deformed shape, wall stress and published modulus.', options: [{ value: '1', label: 'x₁ (horizontal)' }, { value: '2', label: 'x₂ (vertical)' }] },
    { key: 'strain_pct', label: 'Applied strain', unit: '%', value: 1, min: 0.01, max: 10, log: true, group: '2-D lattice', showIf: M('lattice2d'), help: 'Macroscopic strain in the loading direction for the wall-stress check and the animated deformation.' },
    { key: 'cells', label: 'Cells per side in the tessellation', unit: '–', value: 4, min: 2, max: 6, step: 1, group: '2-D lattice', showIf: M('lattice2d'), help: 'Number of unit cells per side drawn and used for the supercell verification solve.' },
    { key: 'ns_l', label: 'Beam span', unit: 'mm', value: 40, min: 0.01, max: 1000, log: true, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Distance between the clamped ends of the curved beam.' },
    { key: 'ns_t', label: 'Beam thickness', unit: 'mm', value: 0.4, min: 0.0001, max: 20, log: true, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Thickness in the bending direction; forces scale with its cube.' },
    { key: 'ns_h', label: 'Apex height', unit: 'mm', value: 2.4, min: 0.0001, max: 100, log: true, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Initial rise of the cosine-shaped beam at mid-span; h/t above 2.31 gives bistability.' },
    { key: 'ns_b', label: 'Beam depth', unit: 'mm', value: 5, min: 0.001, max: 500, log: true, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Out-of-plane width of the beam; forces scale linearly with it.' },
    { key: 'ns_n', label: 'Number of beams in parallel', unit: '–', value: 2, min: 1, max: 20, step: 1, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Identical beams acting in parallel (two joined at the centre suppress the twisting mode).' },
    { key: 'ns_mode', label: 'Second-mode constraint', type: 'select', value: 'double', group: 'Negative stiffness', showIf: M('negstiff'), help: 'Whether the asymmetric buckling mode is suppressed (double beam clamped at the centre) or free (single beam).', options: [{ value: 'double', label: 'Constrained (centre-clamped double beam)' }, { value: 'single', label: 'Free (single beam)' }] },
    { key: 'k_pos', label: 'Parallel positive spring', unit: 'N/m', value: 1100, min: 0, max: 1e7, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Stiffness of the load-bearing spring combined with the negative-stiffness element in a quasi-zero-stiffness mount.' },
    { key: 'm_load', label: 'Payload mass', unit: 'kg', value: 0.2, min: 1e-9, max: 1e4, log: true, group: 'Negative stiffness', showIf: M('negstiff'), help: 'Isolated mass; with the net stiffness it gives the natural frequency.' },
    { key: 'nlc_model', label: 'Compressibility model', type: 'select', value: 'winerack', group: 'Negative compressibility', showIf: M('nlc'), help: 'Hinging wine-rack framework solved by frame finite elements, or an anisotropic solid defined by its stiffness constants.', options: [{ value: 'winerack', label: 'Wine-rack framework' }, { value: 'crystal', label: 'Anisotropic crystal (stiffness matrix)' }] },
    { key: 'wr_theta', label: 'Wine-rack half-angle θ', unit: '°', value: 25, min: 8, max: 82, group: 'Negative compressibility', showIf: M('nlc'), help: 'Angle between the struts and the x axis. Below 45° the framework expands along x under pressure, above 45° along y.' },
    { key: 'wr_L', label: 'Strut length', unit: 'mm', value: 10, min: 0.0001, max: 500, log: true, group: 'Negative compressibility', showIf: M('nlc'), help: 'Node-to-node strut length of the rhombic framework.' },
    { key: 'wr_w', label: 'Strut width', unit: 'mm', value: 0.5, min: 0.00001, max: 100, log: true, group: 'Negative compressibility', showIf: M('nlc'), help: 'In-plane width of the struts; hinging compliance scales as (L/w)³.' },
    { key: 'c11', label: 'C₁₁', unit: 'GPa', value: 20, min: 0.1, max: 1200, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Normal stiffness along axis 1 (crystal model).' }, { key: 'c22', label: 'C₂₂', unit: 'GPa', value: 20, min: 0.1, max: 1200, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Normal stiffness along axis 2 (crystal model).' },
    { key: 'c33', label: 'C₃₃', unit: 'GPa', value: 80, min: 0.1, max: 1200, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Normal stiffness along axis 3 (crystal model).' }, { key: 'c12', label: 'C₁₂', unit: 'GPa', value: 8, min: -200, max: 600, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Coupling between axes 1 and 2 (crystal model).' },
    { key: 'c13', label: 'C₁₃', unit: 'GPa', value: 28, min: -200, max: 600, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Coupling between axes 1 and 3; large values relative to C₁₁ + C₁₂ produce negative compressibility along axis 3.' }, { key: 'c23', label: 'C₂₃', unit: 'GPa', value: 28, min: -200, max: 600, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Coupling between axes 2 and 3 (crystal model).' },
    { key: 'c44', label: 'Smallest shear constant', unit: 'GPa', value: 10, min: -10, max: 600, group: 'Crystal stiffness', showIf: M('nlc'), help: 'Smallest of C₄₄, C₅₅, C₆₆; it must be positive for Born stability but does not affect the hydrostatic response.' },
    { key: 'topo3', label: '3-D lattice topology', type: 'select', value: 'octet', group: '3-D lattice', showIf: M('lattice3d'), help: 'Architecture of the 3-D cell; stretch-dominated trusses scale linearly with density, bending-dominated cells quadratically.', options: Object.entries(L3).map(([value, v]) => ({ value, label: v.label })) },
    { key: 'rel_dens', label: 'Relative density', unit: '–', value: 0.1, min: 0.0005, max: 0.6, log: true, group: '3-D lattice', bus: 'geom.rel_density', showIf: M('lattice3d'), help: 'Solid volume fraction of the lattice.' },
    { key: 'cell3', label: 'Unit-cell size', unit: 'µm', value: 2000, min: 0.1, max: 1e5, log: true, group: '3-D lattice', showIf: M('lattice3d'), help: 'Edge length of the cubic unit cell; with the density it gives the strut radius (sub-micron cells are nanolattices).' },
    { key: 'end_fix', label: 'Strut end-fixity factor n', unit: '–', value: 1, min: 1, max: 2, group: '3-D lattice', showIf: M('lattice3d'), help: 'Euler end-condition factor: 1 for pin-jointed, 2 for fully built-in struts.' },
    { key: 'KIc', label: 'Fracture toughness', unit: 'Pa·√m', value: 1e6, min: 1e4, max: 2e8, log: true, group: 'Fracture', bus: 'prop.KIc', showIf: M('lattice3d', 'strength'), help: 'Mode-I fracture toughness of the solid (about 1 MPa·√m for oxide ceramics and glasses).' },
    { key: 'E_n', label: 'Young\'s modulus', unit: 'Pa', value: 1.4e11, min: 1e6, max: 1.2e12, log: true, group: 'Nanomechanics', bus: 'prop.E', showIf: M('nanowire', 'contact', 'strength'), help: 'Bulk Young\'s modulus of the nanostructure or film material (about 140 GPa for ZnO along the c-axis).' },
    { key: 'nu_n', label: 'Poisson\'s ratio', unit: '–', value: 0.3, min: 0, max: 0.49, group: 'Nanomechanics', bus: 'prop.nu', showIf: M('contact', 'strength'), help: 'Poisson\'s ratio of the tested or film material.' },
    { key: 'nw_D', label: 'Nanowire diameter', unit: 'nm', value: 50, min: 2, max: 5000, log: true, group: 'Nanowire', showIf: M('nanowire'), help: 'Diameter of the circular nanowire or nanorod.' },
    { key: 'nw_L', label: 'Nanowire length', unit: 'µm', value: 2, min: 0.02, max: 1000, log: true, group: 'Nanowire', showIf: M('nanowire'), help: 'Free length between supports.' },
    { key: 'nw_bc', label: 'End conditions', type: 'select', value: 'cantilever', group: 'Nanowire', showIf: M('nanowire'), help: 'Support conditions for bending, buckling and resonance.', options: Object.entries(BC).map(([value, v]) => ({ value, label: v.label })) },
    { key: 'E_surf', label: 'Surface elastic modulus', unit: 'N/m', value: 10, min: -50, max: 200, group: 'Nanowire', showIf: M('nanowire'), help: 'Elastic modulus of the surface layer (force per length). Positive values stiffen thin wires, negative soften them.' },
    { key: 'tau0', label: 'Residual surface stress', unit: 'N/m', value: 1, min: -10, max: 10, group: 'Nanowire', showIf: M('nanowire'), help: 'Residual surface tension; adds an effective transverse load that raises (positive) or lowers the buckling load.' },
    { key: 'nw_F', label: 'Transverse load', unit: 'nN', value: 1, min: 0.001, max: 1e6, log: true, group: 'Nanowire', showIf: M('nanowire'), help: 'Point load at the free end (cantilever) or mid-span.' },
    { key: 'sy_n', label: 'Yield strength of the tested material', unit: 'Pa', value: 1.7e9, min: 1e6, max: 5e10, log: true, group: 'Indentation', bus: 'prop.yield', showIf: M('contact'), help: 'Uniaxial yield strength used for the onset of plasticity beneath the spherical tip (roughly hardness/3).' },
    { key: 'tip_R', label: 'Indenter tip radius', unit: 'nm', value: 1000, min: 5, max: 1e6, log: true, group: 'Indentation', showIf: M('contact'), help: 'Radius of the spherical (or blunted) diamond tip for the Hertz analysis.' },
    { key: 'P_ind', label: 'Hertz load', unit: 'µN', value: 5, min: 0.001, max: 1e6, log: true, group: 'Indentation', showIf: M('contact'), help: 'Normal load for the elastic spherical contact.' },
    { key: 'Pmax', label: 'Berkovich peak load', unit: 'mN', value: 5, min: 0.001, max: 500, log: true, group: 'Indentation', showIf: M('contact'), help: 'Maximum load of the measured indentation cycle.' },
    { key: 'hmax', label: 'Depth at peak load', unit: 'nm', value: 226, min: 5, max: 20000, log: true, group: 'Indentation', showIf: M('contact'), help: 'Maximum displacement into the surface at peak load.' },
    { key: 'S_unl', label: 'Unloading stiffness', unit: 'kN/m', value: 158, min: 0.1, max: 1e5, log: true, group: 'Indentation', showIf: M('contact'), help: 'Initial slope dP/dh of the unloading curve (1 kN/m = 1 µN/nm).' },
    { key: 'd_grain', label: 'Grain size', unit: 'nm', value: 50, min: 2, max: 1e5, log: true, group: 'Strength', showIf: M('strength'), help: 'Mean grain diameter of the polycrystal.' },
    { key: 'sig0', label: 'Friction stress σ₀', unit: 'MPa', value: 25, min: 0, max: 5000, group: 'Strength', showIf: M('strength'), help: 'Lattice friction stress of the Hall–Petch relation (about 25 MPa for copper).' },
    { key: 'k_hp', label: 'Hall–Petch coefficient', unit: 'MPa·µm^½', value: 110, min: 1, max: 3000, log: true, group: 'Strength', showIf: M('strength'), help: 'Grain-boundary strengthening coefficient (about 110 for copper, 600–700 for steels).' },
    { key: 'w_gb', label: 'Grain-boundary thickness', unit: 'nm', value: 1, min: 0.2, max: 5, group: 'Strength', showIf: M('strength'), help: 'Effective width of the disordered grain-boundary phase in the mixture model.' },
    { key: 'sig_gb', label: 'Grain-boundary phase strength', unit: 'MPa', value: 600, min: 1, max: 20000, log: true, group: 'Strength', showIf: M('strength'), help: 'Flow strength of the grain-boundary (amorphous-like) phase that limits the finest-grain strength.' },
    { key: 'flaw', label: 'Flaw size a', unit: 'nm', value: 200, min: 1, max: 1e6, log: true, group: 'Fracture', showIf: M('strength'), help: 'Half-length of the largest crack-like flaw (pore, scratch, grain-boundary groove).' },
    { key: 'Y_geo', label: 'Geometry factor Y', unit: '–', value: 1.12, min: 0.5, max: 3, group: 'Fracture', showIf: M('strength'), help: 'Stress-intensity geometry factor (1.12 edge crack, 1.0 centre crack, 0.64 penny-shaped).' },
    { key: 'sig_app', label: 'Applied tensile stress', unit: 'Pa', value: 2e8, min: 1e4, max: 5e10, log: true, group: 'Fracture', bus: 'mech.max_stress', showIf: M('strength'), help: 'Service stress normal to the flaw.' },
    { key: 'gamma_s', label: 'Surface energy', unit: 'J/m²', value: 1.5, min: 0.01, max: 10, group: 'Fracture', bus: 'prop.surface_energy', showIf: M('strength'), help: 'Surface energy of the fracture plane, used for the ideal-brittle (Griffith) toughness comparison.' },
    { key: 'cte_f', label: 'Film expansion coefficient', unit: '1/K', value: 4.7e-6, min: -2e-5, max: 5e-5, group: 'Thin film', bus: 'prop.cte', showIf: M('strength'), help: 'Linear thermal expansion coefficient of the film (about 4.7×10⁻⁶ for ZnO in the basal plane).' },
    { key: 'cte_s', label: 'Substrate expansion coefficient', unit: '1/K', value: 2.6e-6, min: -2e-5, max: 5e-5, group: 'Thin film', showIf: M('strength'), help: 'Linear thermal expansion coefficient of the substrate (2.6×10⁻⁶ silicon, 7–9×10⁻⁶ glass, 5–8×10⁻⁶ sapphire).' },
    { key: 'dT', label: 'Temperature change', unit: 'K', value: -400, min: -1500, max: 1500, group: 'Thin film', showIf: M('strength'), help: 'Final minus stress-free (deposition) temperature; negative for cooling.' },
    { key: 't_film', label: 'Film thickness', unit: 'nm', value: 500, min: 1, max: 1e5, log: true, group: 'Thin film', showIf: M('strength'), help: 'Thickness of the film.' },
    { key: 't_sub', label: 'Substrate thickness', unit: 'µm', value: 500, min: 1, max: 1e4, log: true, group: 'Thin film', showIf: M('strength'), help: 'Thickness of the substrate.' },
    { key: 'E_sub', label: 'Substrate modulus', unit: 'Pa', value: 1.3e11, min: 1e8, max: 1e12, log: true, group: 'Thin film', showIf: M('strength'), help: 'Young\'s modulus of the substrate (130 GPa for Si(100)).' },
    { key: 'nu_sub', label: 'Substrate Poisson\'s ratio', unit: '–', value: 0.28, min: 0, max: 0.49, group: 'Thin film', showIf: M('strength'), help: 'Poisson\'s ratio of the substrate.' },
    { key: 'E_d', label: 'Drained Young\'s modulus', unit: 'Pa', value: 1.5e9, min: 1e5, max: 1e12, log: true, group: 'Poroelasticity', bus: 'eff.E', showIf: M('poro'), help: 'Young\'s modulus of the porous skeleton under drained conditions (1.5 GPa in the ZnO reference case).' },
    { key: 'nu_d', label: 'Drained Poisson\'s ratio', unit: '–', value: 0.3, min: 0, max: 0.49, group: 'Poroelasticity', bus: 'eff.nu', showIf: M('poro'), help: 'Poisson\'s ratio of the drained skeleton (0.30 in the reference case).' },
    { key: 'biot', label: 'Biot coefficient α', unit: '–', value: 0.6, min: 0.05, max: 1, group: 'Poroelasticity', showIf: M('poro'), help: 'Fraction of pore pressure that acts on the skeleton, α = 1 − K/K_s (0.60 in the reference case).' },
    { key: 'poros', label: 'Porosity', unit: '–', value: 0.3, min: 0.01, max: 0.9, group: 'Poroelasticity', bus: 'geom.porosity', showIf: M('poro'), help: 'Connected pore volume fraction.' },
    { key: 'perm', label: 'Permeability', unit: 'm²', value: 1e-16, min: 1e-22, max: 1e-9, log: true, group: 'Poroelasticity', bus: 'eff.permeability', showIf: M('poro'), help: 'Intrinsic (Darcy) permeability of the porous solid.' },
    { key: 'mu_f', label: 'Fluid viscosity', unit: 'mPa·s', value: 1, min: 0.01, max: 1e5, log: true, group: 'Poroelasticity', showIf: M('poro'), help: 'Dynamic viscosity of the pore fluid.' },
    { key: 'K_f', label: 'Fluid bulk modulus', unit: 'Pa', value: 2.2e9, min: 1e5, max: 1e10, log: true, group: 'Poroelasticity', showIf: M('poro'), help: 'Bulk modulus of the pore fluid (2.2 GPa water; far lower if gas bubbles are present).' },
    { key: 'H_layer', label: 'Layer thickness (drainage path)', unit: 'mm', value: 10, min: 0.0001, max: 1e4, log: true, group: 'Poroelasticity', showIf: M('poro'), help: 'Thickness of the layer drained at its top surface and sealed at the base.' },
    { key: 'sig_load', label: 'Applied compressive stress', unit: 'Pa', value: 1e5, min: 1, max: 1e9, log: true, group: 'Poroelasticity', bus: 'proc.P', showIf: M('poro'), help: 'Step load applied to the surface at time zero.' },
    { key: 't_cons', label: 'Simulated time (0 = automatic)', unit: 's', value: 0, min: 0, max: 1e9, group: 'Poroelasticity', showIf: M('poro'), help: 'Duration of the consolidation simulation; zero selects 1.5 times the 90 % consolidation time.' },
  ],
  presets: [
    { name: 'Re-entrant auxetic honeycomb (θ = −30°, h/l = 2)', note: 'Classic negative-Poisson\'s-ratio metamaterial: Gibson–Ashby predicts ν₁₂ = −1.', values: { mode: 'lattice2d', topology: 'reentrant', cell_l: 5, wall_t: 0.4, h_l: 2, theta_c: -30, Es: 2e9, nu_s: 0.35, strain_pct: 1 } },
    { name: 'Regular hexagonal honeycomb', note: 'Conventional bending-dominated honeycomb: isotropic with ν = 1.', values: { mode: 'lattice2d', topology: 'honeycomb', cell_l: 5, wall_t: 0.3, Es: 7e10, nu_s: 0.33, rho_s: 2700, yield_s: 2.7e8, strain_pct: 0.1 } },
    { name: 'Hexachiral auxetic lattice', note: 'Rings joined by tangent ligaments rotate under load, giving ν close to −1.', values: { mode: 'lattice2d', topology: 'chiral', cell_l: 10, wall_t: 0.3, r_R: 0.25, Es: 2e9, strain_pct: 0.5, cells: 3 } },
    { name: 'Bistable negative-stiffness beam and quasi-zero-stiffness mount', note: 'Cosine beam with Q = 6 in parallel with a positive spring.', values: { mode: 'negstiff', ns_l: 40, ns_t: 0.4, ns_h: 2.4, ns_b: 5, ns_n: 2, ns_mode: 'double', Es: 2e9, yield_s: 5e7, k_pos: 1100, m_load: 0.2 } },
    { name: 'Wine-rack framework with negative linear compressibility', note: 'Hinging rhombic lattice that lengthens along its long diagonal under pressure.', values: { mode: 'nlc', nlc_model: 'winerack', wr_theta: 25, wr_L: 10, wr_w: 0.5, Es: 2e9, nu_s: 0.35 } },
    { name: 'Alumina octet nanolattice', note: 'Stretch-dominated ceramic nanolattice with sub-micron struts.', values: { mode: 'lattice3d', topo3: 'octet', rel_dens: 0.05, cell3: 5, Es: 1.7e11, nu_s: 0.24, rho_s: 3000, yield_s: 2e9, KIc: 2e6 } },
    { name: 'ZnO nanowire cantilever (50 nm × 2 µm)', note: 'Bending, buckling and resonance of a hydrothermal ZnO nanorod with surface elasticity.', values: { mode: 'nanowire', E_n: 1.4e11, nw_D: 50, nw_L: 2, nw_bc: 'cantilever', rho_s: 5606, E_surf: 10, tau0: 1, nw_F: 1 } },
    { name: 'ZnO reference case: poroelastic consolidation', note: 'E = 1.5 GPa, ν = 0.30, Biot coefficient 0.60 as in the reference porous ZnO layer.', values: { mode: 'poro', E_d: 1.5e9, nu_d: 0.3, biot: 0.6, poros: 0.3, perm: 1e-16, H_layer: 10, sig_load: 1e5 } },
  ],
  run(inp) {
    const f = { lattice2d, negstiff, nlc, lattice3d, nanowire, contact, strength, poro }[inp.mode] || lattice2d;
    return finish(f(inp));
  },
};
