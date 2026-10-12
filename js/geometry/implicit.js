// Implicit / lattice geometry engine: triply periodic minimal surfaces, strut lattices (incl. auxetic and chiral
// cells), signed-distance CSG, iso-surface extraction and stochastic porous media. Pure ES module, no dependencies
// beyond the shared numeric kernel. Voxel index convention: i = x + nx*(y + ny*z).
import { rng } from '../core/num.js';

const { sin, cos, sqrt, abs, min, max, PI, floor, round, hypot } = Math;
const MAX_VOX = 2.2e6;
const clampi = (v, a, b) => min(b, max(a, round(Number(v) || 0)));

// ---------------------------------------------------------------- TPMS
// Level-set approximations (leading Fourier terms) of the minimal surfaces; period 2π in x, y, z.
export const TPMS = {
  gyroid: (x, y, z) => sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x),
  schwarzP: (x, y, z) => cos(x) + cos(y) + cos(z),
  schwarzD: (x, y, z) => sin(x) * sin(y) * sin(z) + sin(x) * cos(y) * cos(z) + cos(x) * sin(y) * cos(z) + cos(x) * cos(y) * sin(z),
  neovius: (x, y, z) => 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z),
  iwp: (x, y, z) => 2 * (cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x)) - (cos(2 * x) + cos(2 * y) + cos(2 * z)),
  fischerKochS: (x, y, z) => cos(2 * x) * sin(y) * cos(z) + cos(x) * cos(2 * y) * sin(z) + sin(x) * cos(y) * cos(2 * z),
  lidinoid: (x, y, z) => 0.5 * (sin(2 * x) * cos(y) * sin(z) + sin(2 * y) * cos(z) * sin(x) + sin(2 * z) * cos(x) * sin(y)) - 0.5 * (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)) + 0.15,
  frd: (x, y, z) => 4 * cos(x) * cos(y) * cos(z) - (cos(2 * x) * cos(2 * y) + cos(2 * y) * cos(2 * z) + cos(2 * z) * cos(2 * x)),
};

// Voxelised TPMS. mode 'sheet': solid where |f| <= level (double-sided wall); mode 'solid': solid where f <= level
// (one of the two labyrinths filled). Give `relDensity` (level found by quantile of the sampled field, so the target
// is met to within one voxel count) or `thickness` (dimensionless level-set width: sheet |f| <= thickness/2,
// solid f <= thickness). `grading` varies the target linearly along an axis.
export function tpmsVoxels({ type = 'gyroid', cells = [2, 2, 2], resolution = 24, thickness, relDensity, mode = 'sheet', grading, cellSize_mm = 1 } = {}) {
  const f = TPMS[type];
  if (!f) throw new Error(`Unknown TPMS type "${type}". Available: ${Object.keys(TPMS).join(', ')}`);
  const c = [0, 1, 2].map(i => clampi(cells[i] ?? 1, 1, 64));
  let n = clampi(resolution, 6, 128), capped = null;
  while (n > 6 && n * n * n * c[0] * c[1] * c[2] > MAX_VOX) n--;
  if (n !== clampi(resolution, 6, 128)) capped = `resolution reduced to ${n} voxels per cell to keep the grid below ${MAX_VOX} voxels`;
  const sheet = mode !== 'solid', N1 = n * n * n, u = new Float32Array(N1);
  for (let k = 0, p = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++, p++) {
    const v = f((i + 0.137) / n * 2 * PI, (j + 0.291) / n * 2 * PI, (k + 0.413) / n * 2 * PI); // unequal sub-voxel offsets avoid symmetry ties at the level set
    u[p] = sheet ? abs(v) : v;
  }
  const sorted = Float32Array.from(u).sort();
  const levelFor = rho => rho <= 0 ? -Infinity : rho >= 1 ? Infinity : sorted[min(N1 - 1, max(0, round(rho * N1) - 1))];
  const byDensity = relDensity !== undefined && relDensity !== null || thickness === undefined || thickness === null;
  const target0 = byDensity ? min(1, max(0, Number(relDensity ?? 0.3))) : Number(thickness);
  const toLevel = t => byDensity ? levelFor(min(1, max(0, t))) : (sheet ? t / 2 : t);
  const nx = n * c[0], ny = n * c[1], nz = n * c[2], data = new Uint8Array(nx * ny * nz);
  const ax = grading ? ({ x: 0, y: 1, z: 2 })[grading.axis] ?? 2 : -1, dims = [nx, ny, nz];
  const lev = ax < 0 ? null : Float32Array.from({ length: dims[ax] }, (_, i) => toLevel(grading.from + (grading.to - grading.from) * (dims[ax] > 1 ? i / (dims[ax] - 1) : 0)));
  const L0 = toLevel(target0);
  let count = 0;
  for (let k = 0, p = 0; k < nz; k++) for (let j = 0; j < ny; j++) {
    const row = n * ((j % n) + n * (k % n));
    for (let i = 0; i < nx; i++, p++) {
      const L = ax < 0 ? L0 : lev[ax === 0 ? i : ax === 1 ? j : k];
      if (u[row + (i % n)] <= L) { data[p] = 1; count++; }
    }
  }
  const out = { kind: 'voxels', nx, ny, nz, data, voxel: cellSize_mm / n, units: 'mm', type, mode: sheet ? 'sheet' : 'solid', relDensity: count / data.length, level: ax < 0 ? L0 : null, cells: c, resolution: n };
  if (grading) out.grading = { axis: 'xyz'[ax], from: grading.from, to: grading.to, of: byDensity ? 'relDensity' : 'thickness' };
  if (capped) out.capped = capped;
  return out;
}

// ---------------------------------------------------------------- strut lattices
const C8 = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]];
const F6 = [[.5, .5, 0], [.5, .5, 1], [.5, 0, .5], [.5, 1, .5], [0, .5, .5], [1, .5, .5]];
const dist = (a, b) => hypot(a[0] - b[0], a[1] - b[1], (a[2] || 0) - (b[2] || 0));
function byDist(nodes, d, keep = () => true) {
  const e = [];
  for (let i = 0; i < nodes.length; i++) for (let j = i + 1; j < nodes.length; j++) if (abs(dist(nodes[i], nodes[j]) - d) < 1e-6 * max(1, d) && keep(i, j)) e.push([i, j]);
  return e;
}
const deg = a => a * PI / 180;
// Gibson–Ashby hexagonal cell (vertical ribs h, inclined ribs l at angle θ to the horizontal; θ < 0 gives the re-entrant, auxetic bow-tie).
function hexCell(p, dflt) {
  let th = deg(p.theta_deg ?? dflt.theta), hl = p.h_over_l ?? dflt.hl, note;
  th = min(deg(80), max(deg(-80), th));
  if (hl + 2 * sin(th) < 0.05) { hl = 0.05 - 2 * sin(th); note = `h/l raised to ${hl.toFixed(3)} so that opposite re-entrant vertices do not touch`; }
  const l = 1, h = hl, c = l * cos(th), s = l * sin(th), W = 2 * c, H = 2 * (h + s), k = 1 / W;
  const nodes = [[0, 0], [0, h], [c, h + s], [c, 2 * h + s], [0, H], [W, 0], [W, h], [W, H]].map(q => [q[0] * k, q[1] * k, 0]);
  const nu12 = abs(sin(th)) < 1e-9 ? null : cos(th) ** 2 / ((hl + sin(th)) * sin(th));
  return { nodes, edges: [[0, 1], [1, 2], [2, 3], [3, 4], [2, 6], [5, 6], [3, 7]], size: [1, H * k, 0], note,
    analytic: { model: 'Gibson–Ashby flexure model, thin ribs (t/l ≪ 1), loading along x', nu12, l_mm_per_cell: k, relDensityFormula: '(t/l)(h/l + 2) / (2 cosθ (h/l + sinθ))', relDensityPerT_over_l: (hl + 2) / (2 * cos(th) * (hl + sin(th))) } };
}
function sitesCell(sites, W, H, d) {
  const nodes = sites.filter(q => q[0] > -1e-9 && q[0] < W + 1e-9 && q[1] > -1e-9 && q[1] < H + 1e-9).map(q => [q[0] / W, q[1] / W, 0]);
  return { nodes, edges: byDist(nodes, d / W), size: [1, H / W, 0] };
}
const cell = (name, dim, make, info = {}) => { const d = make({}); return { name, dim, ...info, nodes: d.nodes, edges: d.edges, size: d.size || [1, 1, dim === 3 ? 1 : 0], make }; };

export const LATTICES = {
  simpleCubic: cell('Simple cubic', 3, () => ({ nodes: C8, edges: byDist(C8, 1) })),
  bcc: cell('Body-centred cubic', 3, () => { const n = [...C8, [.5, .5, .5]]; return { nodes: n, edges: byDist(n, sqrt(3) / 2) }; }),
  fcc: cell('Face-centred cubic', 3, () => { const n = [...C8, ...F6]; return { nodes: n, edges: byDist(n, sqrt(2) / 2, i => i < 8) }; }),
  octet: cell('Octet truss', 3, () => { const n = [...C8, ...F6]; return { nodes: n, edges: byDist(n, sqrt(2) / 2) }; }),
  kelvin: cell('Kelvin (truncated octahedron)', 3, () => {
    const n = [];
    for (const [a, b, c] of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) for (const s1 of [-1, 1]) for (const s2 of [-1, 2]) { const v = [0, 0, 0]; v[a] = 0; v[b] = s1; v[c] = s2 < 0 ? -2 : 2; n.push(v.map(q => q / 4 + 0.5)); }
    return { nodes: n, edges: byDist(n, sqrt(2) / 4) };
  }),
  diamond: cell('Diamond', 3, () => { const n = [...C8, ...F6, [.25, .25, .25], [.75, .75, .25], [.75, .25, .75], [.25, .75, .75]]; return { nodes: n, edges: byDist(n, sqrt(3) / 4) }; }),
  reentrant3D: cell('3-D re-entrant (auxetic)', 3, p => {
    let th = deg(p.theta_deg ?? -20), hl = p.h_over_l ?? 2, note;
    th = min(deg(-2), max(deg(-75), -abs(th)));
    if (hl + 2 * sin(th) < 0.05) { hl = 0.05 - 2 * sin(th); note = `h/l raised to ${hl.toFixed(3)} to avoid self-contact`; }
    const c = cos(th), s = sin(th), W = 2 * c, H = 2 * (hl + s), k = 1 / W, m = W / 2;
    const mid = [[0, m], [W, m], [m, 0], [m, W]], nodes = [[m, m, hl + s], [m, m, 2 * hl + s]], edges = [[0, 1]];
    for (const [x, y] of mid) { const b = nodes.length; nodes.push([x, y, 0], [x, y, hl], [x, y, H]); edges.push([b, b + 1], [b + 1, 0], [b + 2, 1]); }
    return { nodes: nodes.map(q => q.map(v => v * k)), edges, size: [1, 1, H * k], note };
  }, { auxetic: true }),
  honeycomb2D: cell('Hexagonal honeycomb', 2, p => hexCell(p, { theta: 30, hl: 1 })),
  reentrant2D: cell('Re-entrant honeycomb (auxetic)', 2, p => hexCell({ ...p, theta_deg: -abs(p.theta_deg ?? 20) }, { theta: -20, hl: 2 }), { auxetic: true }),
  chiral2D: cell('Tetrachiral (rings with tangent ligaments)', 2, p => {
    const r = min(0.4, max(0.03, p.r_over_L ?? 0.18)), a = Math.acos(2 * r), seg = 12, nodes = [], edges = [];
    for (let i = 0; i < seg; i++) { const t = a + i * 2 * PI / seg; nodes.push([.5 + r * cos(t), .5 + r * sin(t), 0]); edges.push([i, (i + 1) % seg]); }
    [[1, .5], [.5, 1], [0, .5], [.5, 0]].forEach((q, i) => { nodes.push([q[0], q[1], 0]); edges.push([i * 3, seg + i]); });
    return { nodes, edges, analytic: { model: 'Prall–Lakes chiral honeycomb (ideal, rigid rings): ν → −1', nu12: -1 } };
  }, { auxetic: true }),
  kagome2D: cell('Kagome', 2, () => {
    const s = [], h = sqrt(3);
    for (let i = -2; i < 4; i++) for (let j = -1; j < 4; j++) { const x = 2 * i + j, y = h * j; s.push([x, y], [x + 1, y], [x + .5, y + h / 2]); }
    return sitesCell(s, 2, 2 * h, 1);
  }),
  triangular2D: cell('Triangular', 2, () => { const h = sqrt(3); return { nodes: [[0, 0, 0], [1, 0, 0], [.5, h / 2, 0], [0, h, 0], [1, h, 0]], edges: [[0, 1], [3, 4], [0, 2], [1, 2], [2, 3], [2, 4], [2, 2, 1, 0, 0]], size: [1, h, 0] }; }),
  rotatingSquares2D: cell('Rotating squares (auxetic)', 2, p => {
    const phi = deg(min(44, max(1, abs(p.theta_deg ?? 20)))), pch = 0.5, s = pch / (cos(phi) + sin(phi)), r = s / sqrt(2), nodes = [], edges = [];
    [[.25, .25, 1], [.75, .25, -1], [.25, .75, -1], [.75, .75, 1]].forEach(([cx, cy, sg]) => {
      const b = nodes.length;
      for (let q = 0; q < 4; q++) { const t = PI / 4 + q * PI / 2 + sg * phi; nodes.push([cx + r * cos(t), cy + r * sin(t), 0]); edges.push([b + q, b + (q + 1) % 4]); }
      edges.push([b, b + 2]);
    });
    return { nodes, edges, analytic: { model: 'Grima–Evans rigid rotating squares: ν = −1 at every opening angle', nu12: -1, squareSide_per_cell: s } };
  }, { auxetic: true }),
  starHoneycomb2D: cell('Four-point star honeycomb (auxetic)', 2, p => {
    const R = 0.34, ri = R * min(0.9, max(0.15, p.inner ?? 0.42)), nodes = [], edges = [];
    for (let q = 0; q < 8; q++) { const t = q * PI / 4, rr = q % 2 ? ri : R; nodes.push([.5 + rr * cos(t), .5 + rr * sin(t), 0]); edges.push([q, (q + 1) % 8]); }
    [[1, .5], [.5, 1], [0, .5], [.5, 0]].forEach((q, i) => { nodes.push([q[0], q[1], 0]); edges.push([2 * i, 8 + i]); });
    return { nodes, edges };
  }, { auxetic: true }),
};

// Tessellate a unit cell. Returns a `lattice` geometry plus first-order relative density (strut volume / cell volume,
// node overlap ignored, so it overestimates above ~0.2), nodal connectivity and the Maxwell rigidity number.
export function buildLattice({ type = 'octet', cells = [2, 2, 2], cellSize_mm = 10, strut_mm = 1, params = {}, grading } = {}) {
  const L = LATTICES[type];
  if (!L) throw new Error(`Unknown lattice type "${type}". Available: ${Object.keys(LATTICES).join(', ')}`);
  const uc = L.make(params || {}), dim = L.dim, size = uc.size || [1, 1, dim === 3 ? 1 : 0];
  const n = [clampi(cells[0] ?? 1, 1, 200), clampi(cells[1] ?? 1, 1, 200), dim === 3 ? clampi(cells[2] ?? 1, 1, 200) : 1];
  let capped = null;
  while (n[0] * n[1] * n[2] * uc.edges.length > 60000) { const k = n.indexOf(max(...n)); n[k]--; capped = 'cell count reduced to keep the lattice below 60 000 struts'; }
  const s = Number(cellSize_mm) > 0 ? Number(cellSize_mm) : 10, t = Number(strut_mm) > 0 ? Number(strut_mm) : s / 10;
  const nodes = [], edges = [], idx = new Map(), seen = new Set();
  const add = p => { const k = (round(p[0] / s * 4096) + 1048576) * 2097152 + round(p[1] / s * 4096) + 1048576, kz = round(p[2] / s * 4096); let m = idx.get(k); if (!m) idx.set(k, m = new Map()); let i = m.get(kz); if (i === undefined) { i = nodes.length; nodes.push(p); m.set(kz, i); } return i; }; // nodes closer than 1/4096 of a cell are merged
  const pos = (q, c) => [(q[0] + c[0] * size[0]) * s, (q[1] + c[1] * size[1]) * s, ((q[2] || 0) + c[2] * (size[2] || 0)) * s];
  for (let k = 0; k < n[2]; k++) for (let j = 0; j < n[1]; j++) for (let i = 0; i < n[0]; i++) for (const e of uc.edges) {
    const c2 = [i + (e[2] || 0), j + (e[3] || 0), k + (e[4] || 0)];
    if (c2[0] < 0 || c2[1] < 0 || c2[2] < 0 || c2[0] >= n[0] || c2[1] >= n[1] || c2[2] >= n[2]) continue;
    const a = add(pos(uc.nodes[e[0]], [i, j, k])), b = add(pos(uc.nodes[e[1]], c2));
    if (a === b) continue;
    const key = a < b ? a * 1e7 + b : b * 1e7 + a;
    if (!seen.has(key)) { seen.add(key); edges.push([a, b]); }
  }
  // periodic (per unit cell) bookkeeping: struts and nodes that are translates of one another by a cell vector count once
  const wrap = (v, a) => { if (!(size[a] > 0)) return 0; let u = (v / size[a]) % 1; if (u < 0) u += 1; u = round(u * 1e6) / 1e6; return u >= 1 ? 0 : u; };
  const eSeen = new Set(), nSeen = new Set(); let bEff = 0, lenEff = 0;
  for (const e of uc.edges) {
    const A = uc.nodes[e[0]], B = uc.nodes[e[1]], B2 = [B[0] + (e[2] || 0) * size[0], B[1] + (e[3] || 0) * size[1], (B[2] || 0) + (e[4] || 0) * (size[2] || 0)];
    let d = [B2[0] - A[0], B2[1] - A[1], B2[2] - (A[2] || 0)]; const L = hypot(d[0], d[1], d[2]); if (!(L > 0)) continue;
    const sg = (abs(d[0]) > 1e-9 ? d[0] : abs(d[1]) > 1e-9 ? d[1] : d[2]) < 0 ? -1 : 1; d = d.map(q => round(sg * q / L * 1e5));
    const key = [0, 1, 2].map(a => wrap(((A[a] || 0) + B2[a]) / 2, a)).join(',') + '|' + d.join(',');
    nSeen.add([0, 1, 2].map(a => wrap(A[a] || 0, a)).join(',')); nSeen.add([0, 1, 2].map(a => wrap(B2[a], a)).join(','));
    if (!eSeen.has(key)) { eSeen.add(key); bEff++; lenEff += L * s; }
  }
  const jEff = nSeen.size || 1;
  const r = t / 2, vol = dim === 3 ? size[0] * size[1] * size[2] * s ** 3 : size[0] * size[1] * s * s;
  const relDensity = min(1, dim === 3 ? PI * r * r * lenEff / vol : t * lenEff / vol);
  const j = nodes.length, b = edges.length, Zp = 2 * bEff / jEff;
  const stretch = dim === 3 ? Zp >= 12 - 1e-9 : Zp >= 6 - 1e-9 || type === 'kagome2D';
  const out = {
    kind: 'lattice', nodes, edges, radius: r, units: 'mm', type, dim, cells: n, cellSize_mm: s, strut_mm: t,
    period: [size[0] * s, size[1] * s, (size[2] || 0) * s],
    relDensity, relDensityNote: dim === 3 ? 'Σ πr²L / V per periodic cell; node overlap neglected (overestimates when ρ* ≳ 0.2)' : 'Σ tL / A per periodic cell with in-plane wall thickness t = strut_mm; node overlap neglected',
    nodeConnectivity: 2 * b / max(1, j), nodeConnectivityPeriodic: Zp,
    maxwell: dim === 3 ? b - 3 * j + 6 : b - 2 * j + 3, maxwellFormula: dim === 3 ? 'M = b − 3j + 6' : 'M = b − 2j + 3',
    dominated: stretch ? 'stretch' : 'bending',
    dominatedBasis: type === 'kagome2D' ? 'periodic kagome (Z = 4) carries macroscopic loads by stretching although it retains finite mechanisms' : `periodic connectivity Z = ${Zp.toFixed(2)} against the rigidity threshold ${dim === 3 ? 12 : 6} for similarly situated nodes (Deshpande–Ashby–Fleck criterion)`,
    auxetic: !!L.auxetic,
  };
  if (uc.analytic) out.analytic = uc.analytic;
  if (uc.note) out.note = uc.note;
  if (capped) out.capped = capped;
  if (grading && (grading.from !== undefined)) {
    const a = ({ x: 0, y: 1, z: 2 })[grading.axis] ?? (dim === 3 ? 2 : 1), ext = max(1e-12, n[a] * size[a] * s);
    out.radii = edges.map(e => { const m = (nodes[e[0]][a] + nodes[e[1]][a]) / 2 / ext; return r * (grading.from + (grading.to - grading.from) * m); });
    out.grading = { axis: 'xyz'[a], from: grading.from, to: grading.to, of: 'strut radius multiplier' };
  }
  return out;
}

// ---------------------------------------------------------------- signed-distance CSG
const v3 = c => c || [0, 0, 0];
const h2 = (a, b) => sqrt(a * a + b * b), h3 = (a, b, c) => sqrt(a * a + b * b + c * c); // Math.hypot is slow in hot loops
export function sdf(f, x, y, z) { return f(x, y, z); }
function rotMat([rx = 0, ry = 0, rz = 0] = []) { // R = Rz·Ry·Rx, degrees
  const [a, b, c] = [rx, ry, rz].map(deg), ca = cos(a), sa = sin(a), cb = cos(b), sb = sin(b), cc = cos(c), sc = sin(c);
  return [[cc * cb, cc * sb * sa - sc * ca, cc * sb * ca + sc * sa], [sc * cb, sc * sb * sa + cc * ca, sc * sb * ca - cc * sa], [-sb, cb * sa, cb * ca]];
}
Object.assign(sdf, {
  sphere: (r = 1, c) => { const [a, b, d] = v3(c); return (x, y, z) => h3(x - a, y - b, z - d) - r; },
  box: (h = [1, 1, 1], c) => { const [a, b, d] = v3(c); return (x, y, z) => { const qx = abs(x - a) - h[0], qy = abs(y - b) - h[1], qz = abs(z - d) - h[2]; return h3(max(qx, 0), max(qy, 0), max(qz, 0)) + min(max(qx, qy, qz), 0); }; },
  cylinder: (r = 1, h = 2, c) => { const [a, b, d] = v3(c); return (x, y, z) => { const q = h2(x - a, y - b) - r, w = abs(z - d) - h / 2; return h2(max(q, 0), max(w, 0)) + min(max(q, w), 0); }; },
  cone: (r1 = 1, r2 = 0, h = 2, c) => { const [a, b, d] = v3(c), k = h / hypot(h, r1 - r2); return (x, y, z) => { const t = min(1, max(0, (z - d) / h + 0.5)); return max(abs(z - d) - h / 2, (h2(x - a, y - b) - (r1 + (r2 - r1) * t)) * k); }; },
  capsule: (p, q, r = 1) => { const d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]], dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] || 1e-30; return (x, y, z) => { const t = min(1, max(0, ((x - p[0]) * d[0] + (y - p[1]) * d[1] + (z - p[2]) * d[2]) / dd)); return h3(x - p[0] - t * d[0], y - p[1] - t * d[1], z - p[2] - t * d[2]) - r; }; },
  torus: (R = 2, r = 0.5, c) => { const [a, b, d] = v3(c); return (x, y, z) => h2(h2(x - a, y - b) - R, z - d) - r; },
  union: (...f) => (x, y, z) => { let m = Infinity; for (const g of f) { const v = g(x, y, z); if (v < m) m = v; } return m; },
  intersect: (...f) => (x, y, z) => { let m = -Infinity; for (const g of f) { const v = g(x, y, z); if (v > m) m = v; } return m; },
  subtract: (a, ...f) => (x, y, z) => { let m = a(x, y, z); for (const g of f) { const v = -g(x, y, z); if (v > m) m = v; } return m; },
  smoothUnion: (a, b, k = 0.1) => (x, y, z) => { const p = a(x, y, z), q = b(x, y, z), h = max(0, min(1, 0.5 + 0.5 * (q - p) / k)); return q + (p - q) * h - k * h * (1 - h); },
  shell: (f, t = 0.1) => (x, y, z) => abs(f(x, y, z)) - t / 2,
  repeat: (f, p) => { const w = (v, q) => q > 0 ? v - q * round(v / q) : v; return (x, y, z) => f(w(x, p[0]), w(y, p[1]), w(z, p[2])); },
  translate: (f, t) => (x, y, z) => f(x - t[0], y - t[1], z - t[2]),
  scale: (f, s) => { const q = typeof s === 'number' ? [s, s, s] : s, m = min(abs(q[0]), abs(q[1]), abs(q[2])) || 1; return (x, y, z) => f(x / q[0], y / q[1], z / q[2]) * m; },
  rotate: (f, r) => { const R = Array.isArray(r) && Array.isArray(r[0]) ? r : rotMat(r); return (x, y, z) => f(R[0][0] * x + R[1][0] * y + R[2][0] * z, R[0][1] * x + R[1][1] * y + R[2][1] * z, R[0][2] * x + R[1][2] * y + R[2][2] * z); },
  rotationMatrix: rotMat,
});

// Sample an SDF on a regular grid (cell centres). resolution = voxels along the longest box edge (capped).
export function sdfToVoxels(sdfFn, bbox = [[-1, -1, -1], [1, 1, 1]], resolution = 48) {
  const [lo, hi] = bbox, ext = [0, 1, 2].map(i => max(1e-12, hi[i] - lo[i]));
  let n = clampi(resolution, 4, 256), capped = null, d;
  for (; ; n--) { d = [0, 1, 2].map(i => max(1, Math.ceil(ext[i] / (max(...ext) / n) - 1e-9))); if (d[0] * d[1] * d[2] <= MAX_VOX || n <= 4) break; }
  if (n !== clampi(resolution, 4, 256)) capped = `resolution reduced to ${n} to keep the grid below ${MAX_VOX} voxels`;
  const h = max(...ext) / n, [nx, ny, nz] = d, data = new Uint8Array(nx * ny * nz), field = new Float32Array(nx * ny * nz);
  let solid = 0;
  for (let k = 0, p = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, p++) {
    const v = sdfFn(lo[0] + (i + .5) * h, lo[1] + (j + .5) * h, lo[2] + (k + .5) * h);
    field[p] = v; if (v < 0) { data[p] = 1; solid++; }
  }
  const out = { kind: 'voxels', nx, ny, nz, data, voxel: h, units: 'mm', origin: [...lo], field, fieldKind: 'sdf', solidFraction: solid / data.length };
  if (capped) out.capped = capped;
  return out;
}

// Iso-surface extraction. Uses the marching-tetrahedra variant of marching cubes (each cube split into six
// tetrahedra sharing the main diagonal): unambiguous, no lookup tables, and watertight; it produces roughly twice the
// triangles of table-driven marching cubes. Input: a voxels geometry (solid = data > iso, default iso 0.5), or one
// carrying a signed-distance `field` (fieldKind 'sdf', or no `data`; solid = field < iso, default 0). The grid is padded by one layer of
// "outside", so surfaces are closed at the domain boundary.
export function marchingCubes(g, iso) {
  const { nx, ny, nz } = g, h = g.voxel || 1, o = g.origin || [0, 0, 0];
  const useField = !!g.field && (g.fieldKind === 'sdf' || !g.data);
  const X = nx + 2, Y = ny + 2, Z = nz + 2, XY = X * Y, NT = X * Y * Z, v = new Float32Array(NT);
  const lv = iso ?? (useField ? 0 : 0.5), src = useField ? g.field : g.data;
  if (!src || src.length < nx * ny * nz) throw new Error('marchingCubes: grid data missing or shorter than nx·ny·nz');
  for (let k = 0; k < Z; k++) for (let j = 0; j < Y; j++) for (let i = 0; i < X; i++) {
    const ci = min(nx - 1, max(0, i - 1)), cj = min(ny - 1, max(0, j - 1)), ck = min(nz - 1, max(0, k - 1));
    const raw = src[ci + nx * (cj + ny * ck)], val = useField ? lv - raw : (raw > lv ? 0.5 : -0.5);
    const pad = i === 0 || j === 0 || k === 0 || i === X - 1 || j === Y - 1 || k === Z - 1;
    v[i + X * (j + Y * k)] = pad ? (val > 0 ? -val : val || -1e-6) : val;
  }
  // vertices and faces are accumulated in typed arrays; edge → vertex lookup uses a flat table when it fits in memory
  const table = NT <= 1.2e6 ? new Int32Array(NT * 8) : null, vmap = table ? null : new Map();
  let vx = new Float64Array(3 * 8192), nv = 0, fc = new Int32Array(3 * 16384), nf = 0;
  const off = [0, 1, X, X + 1, XY, XY + 1, XY + X, XY + X + 1], TETS = [0, 1, 3, 7, 0, 1, 5, 7, 0, 2, 3, 7, 0, 2, 6, 7, 0, 4, 5, 7, 0, 4, 6, 7];
  let bi = 0, bj = 0, bk = 0, base = 0;
  // tetrahedron corners form a bitwise-monotone chain, so an edge is identified by its lower corner and the differing bits
  const vert = (p, q) => {
    const lo = (p & q) === p ? p : q, key = (base + off[lo]) * 8 + (p ^ q); let id = table ? table[key] - 1 : (vmap.get(key) ?? -1);
    if (id < 0) {
      const a = v[base + off[p]], t = a / (a - v[base + off[q]]), px = p & 1, py = p >> 1 & 1, pz = p >> 2, qx = q & 1, qy = q >> 1 & 1, qz = q >> 2;
      id = nv++; if (table) table[key] = id + 1; else vmap.set(key, id);
      if (3 * nv > vx.length) { const n = new Float64Array(vx.length * 2); n.set(vx); vx = n; }
      vx[3 * id] = o[0] + (bi + px + t * (qx - px) - 0.5) * h; vx[3 * id + 1] = o[1] + (bj + py + t * (qy - py) - 0.5) * h; vx[3 * id + 2] = o[2] + (bk + pz + t * (qz - pz) - 0.5) * h;
    }
    return id;
  };
  let d0 = 0, d1 = 0, d2 = 0;
  const tri = (a, b, c) => {
    const ax = vx[3 * a], ay = vx[3 * a + 1], az = vx[3 * a + 2], ux = vx[3 * b] - ax, uy = vx[3 * b + 1] - ay, uz = vx[3 * b + 2] - az, wx = vx[3 * c] - ax, wy = vx[3 * c + 1] - ay, wz = vx[3 * c + 2] - az;
    if (3 * nf + 3 > fc.length) { const n = new Int32Array(fc.length * 2); n.set(fc); fc = n; }
    fc[3 * nf] = a; if ((uy * wz - uz * wy) * d0 + (uz * wx - ux * wz) * d1 + (ux * wy - uy * wx) * d2 >= 0) { fc[3 * nf + 1] = b; fc[3 * nf + 2] = c; } else { fc[3 * nf + 1] = c; fc[3 * nf + 2] = b; } nf++;
  };
  const ins = new Int32Array(4), outs = new Int32Array(4);
  for (bk = 0; bk < Z - 1; bk++) for (bj = 0; bj < Y - 1; bj++) for (bi = 0; bi < X - 1; bi++) {
    base = bi + X * (bj + Y * bk); let m = 0;
    if (v[base] > 0) m |= 1; if (v[base + 1] > 0) m |= 2; if (v[base + X] > 0) m |= 4; if (v[base + X + 1] > 0) m |= 8; if (v[base + XY] > 0) m |= 16; if (v[base + XY + 1] > 0) m |= 32; if (v[base + XY + X] > 0) m |= 64; if (v[base + XY + X + 1] > 0) m |= 128;
    if (m === 0 || m === 255) continue;
    for (let ti = 0; ti < 24; ti += 4) {
      let ni = 0, no = 0;
      for (let q = 0; q < 4; q++) { const c = TETS[ti + q]; if (m >> c & 1) ins[ni++] = c; else outs[no++] = c; }
      if (ni === 0 || ni === 4) continue;
      d0 = d1 = d2 = 0;
      for (let q = 0; q < no; q++) { const c = outs[q]; d0 += (c & 1) / no; d1 += (c >> 1 & 1) / no; d2 += (c >> 2) / no; }
      for (let q = 0; q < ni; q++) { const c = ins[q]; d0 -= (c & 1) / ni; d1 -= (c >> 1 & 1) / ni; d2 -= (c >> 2) / ni; }
      if (ni === 1) tri(vert(ins[0], outs[0]), vert(ins[0], outs[1]), vert(ins[0], outs[2]));
      else if (ni === 3) tri(vert(outs[0], ins[0]), vert(outs[0], ins[1]), vert(outs[0], ins[2]));
      else { const p0 = vert(ins[0], outs[0]), p1 = vert(ins[0], outs[1]), p2 = vert(ins[1], outs[1]), p3 = vert(ins[1], outs[0]); tri(p0, p1, p2); tri(p0, p2, p3); }
    }
  }
  const vertices = new Array(nv), faces = new Array(nf);
  for (let i = 0; i < nv; i++) vertices[i] = [vx[3 * i], vx[3 * i + 1], vx[3 * i + 2]];
  for (let i = 0; i < nf; i++) faces[i] = [fc[3 * i], fc[3 * i + 1], fc[3 * i + 2]];
  return { kind: 'mesh', vertices, faces, units: g.units || 'mm', method: 'marching tetrahedra (6 per cube)' };
}

// Struts as closed prisms with `segments` sides and nodes as icosahedra. The result is a set of overlapping closed
// shells (not a Boolean union): slicers and voxelisers that union shells accept it; FE meshing needs a union first.
const ICO = (() => { const t = (1 + sqrt(5)) / 2, v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(p => p.map(q => q / hypot(1, t))); return { v, f: [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]] }; })();
export function latticeToMesh(lattice, { segments = 6, nodeSpheres = true, maxStruts = 20000 } = {}) {
  const seg = clampi(segments, 3, 24), V = [], F = [], N = lattice.nodes, r0 = lattice.radius || 0.05;
  const E = lattice.edges.length > maxStruts ? lattice.edges.slice(0, maxStruts) : lattice.edges;
  const rmax = new Float64Array(N.length);
  E.forEach((e, ei) => {
    const A = N[e[0]], B = N[e[1]], r = lattice.radii ? lattice.radii[ei] : r0;
    const d = [B[0] - A[0], B[1] - A[1], (B[2] || 0) - (A[2] || 0)], L = hypot(...d);
    if (!(L > 0)) return;
    rmax[e[0]] = max(rmax[e[0]], r); rmax[e[1]] = max(rmax[e[1]], r);
    const w = d.map(q => q / L), ref = abs(w[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    let u = [w[1] * ref[2] - w[2] * ref[1], w[2] * ref[0] - w[0] * ref[2], w[0] * ref[1] - w[1] * ref[0]]; const ul = hypot(...u); u = u.map(q => q / ul);
    const t = [w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0]], b = V.length;
    for (let s = 0; s < seg; s++) { const a = 2 * PI * s / seg, c = r * cos(a), sn = r * sin(a), q = [c * u[0] + sn * t[0], c * u[1] + sn * t[1], c * u[2] + sn * t[2]]; V.push([A[0] + q[0], A[1] + q[1], (A[2] || 0) + q[2]], [B[0] + q[0], B[1] + q[1], (B[2] || 0) + q[2]]); }
    for (let s = 0; s < seg; s++) { const a0 = b + 2 * s, a1 = b + 2 * ((s + 1) % seg); F.push([a0, a1, a1 + 1], [a0, a1 + 1, a0 + 1]); }
    for (let s = 1; s < seg - 1; s++) { F.push([b, b + 2 * (s + 1), b + 2 * s], [b + 1, b + 2 * s + 1, b + 2 * (s + 1) + 1]); }
  });
  if (nodeSpheres) N.forEach((p, i) => { if (!rmax[i]) return; const b = V.length, r = rmax[i] * 1.05; for (const q of ICO.v) V.push([p[0] + r * q[0], p[1] + r * q[1], (p[2] || 0) + r * q[2]]); for (const f of ICO.f) F.push([b + f[0], b + f[1], b + f[2]]); });
  const out = { kind: 'mesh', vertices: V, faces: F, units: lattice.units || 'mm', shells: 'overlapping closed shells (struts and node spheres), not a Boolean union' };
  if (E.length < lattice.edges.length) out.truncated = `only the first ${maxStruts} of ${lattice.edges.length} struts were meshed`;
  return out;
}

// ---------------------------------------------------------------- voxel analysis
// Breadth-first flood through voxels satisfying `inPhase`, started from the low face of `axis` (6-connectivity).
function flood(mask, want, nx, ny, nz, axis) {
  const N = nx * ny * nz, distv = new Int32Array(N).fill(-1), q = new Int32Array(N), nxy = nx * ny, last = [nx, ny, nz][axis] - 1;
  let head = 0, tail = 0, reached = 0, sum = 0;
  const seedStep = axis === 0 ? nx : 1, seedEnd = axis === 2 ? nxy : N;
  if (axis === 1) { for (let z = 0; z < nz; z++) for (let x = 0; x < nx; x++) { const p = x + nxy * z; if (mask[p] === want) { distv[p] = 0; q[tail++] = p; } } }
  else for (let p = 0; p < seedEnd; p += seedStep) if (mask[p] === want) { distv[p] = 0; q[tail++] = p; }
  while (head < tail) {
    const p = q[head++], d = distv[p] + 1, x = p % nx, yz = (p / nx) | 0, y = yz % ny, z = (yz / ny) | 0;
    if ((axis === 0 ? x : axis === 1 ? y : z) === last) { reached++; sum += d - 1; }
    let n;
    if (x > 0 && distv[n = p - 1] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
    if (x < nx - 1 && distv[n = p + 1] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
    if (y > 0 && distv[n = p - nx] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
    if (y < ny - 1 && distv[n = p + nx] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
    if (z > 0 && distv[n = p - nxy] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
    if (z < nz - 1 && distv[n = p + nxy] < 0 && mask[n] === want) { distv[n] = d; q[tail++] = n; }
  }
  return { percolates: reached > 0, tortuosity: reached > 0 && last > 0 ? (sum / reached) / last : null, connected: tail };
}

// Phase fractions, interfacial area, percolation, geometric tortuosity and chord-length (pore-size) statistics.
// Interfacial area per unit volume of a two-phase voxel field from its two-point (phase-change) statistics:
// S_v = −4·dS₂/dr at r = 0 (Debye; Berryman & Blair), i.e. S_v = 2·⟨C′(0)⟩ with C(r) the probability that two points a
// distance r apart lie in different phases, averaged over directions. Thirteen lattice directions are used (3 axes, 6 face
// diagonals, 4 space diagonals) with the Voronoi-cell weights of the Cauchy–Crofton discretisation on a cubic grid
// (Ohser & Mücklich), and C′(0) is extrapolated to zero lag from one and two lattice steps, C′(0) ≈ [4C(d) − C(2d)]/(2d)
// (Richardson). The first-order, axes-only form of the same estimator is the voxel face count × 2/3; it reads low where
// necks and gaps approach the voxel size and is biased for surfaces that are not isotropically oriented.
// mask: 0/1 per voxel (index = x + nx·(y + ny·z)). Returns S_v in units of 1/voxel.
const CROFTON_W = [0.091555782, 0.073961256, 0.070391280];   // weight per direction: axis, face diagonal, space diagonal (3·w₁ + 6·w₂ + 4·w₃ = 1)
const CROFTON_D = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1], [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]];
export function twoPointSurface(mask, nx, ny, nz) {
  let sv = 0, wsum = 0;
  for (const [dx, dy, dz] of CROFTON_D) {
    const nn = abs(dx) + abs(dy) + abs(dz), w = CROFTON_W[nn - 1], len = sqrt(nn), c = [0, 0];
    for (let lag = 1; lag <= 2; lag++) {
      const ox = dx * lag, oy = dy * lag, oz = dz * lag, off = ox + nx * (oy + ny * oz);
      const x0 = max(0, -ox), x1 = nx - max(0, ox), y0 = max(0, -oy), y1 = ny - max(0, oy), z0 = max(0, -oz), z1 = nz - max(0, oz);
      if (x1 <= x0 || y1 <= y0 || z1 <= z0) { c[lag - 1] = -1; continue; }
      let f = 0;
      for (let z = z0; z < z1; z++) for (let y = y0; y < y1; y++) { let p = x0 + nx * (y + ny * z); for (let x = x0; x < x1; x++, p++) if (mask[p] !== mask[p + off]) f++; }
      c[lag - 1] = f / ((x1 - x0) * (y1 - y0) * (z1 - z0));
    }
    if (c[0] < 0) continue;
    sv += w * (4 * c[0] - (c[1] < 0 ? 2 * c[0] : c[1])) / len; wsum += w;
  }
  return wsum ? sv / wsum : 0;
}
// Phase fractions, interfacial area, percolation, geometric tortuosity and chord-length (pore-size) statistics.
export function voxelStats(vox, { porePhase = 0, surfaceOnly = false } = {}) {
  const { nx, ny, nz, data } = vox, N = nx * ny * nz, h = vox.voxel || 1, counts = new Uint32Array(256);
  for (let p = 0; p < N; p++) counts[data[p]]++;
  const fractions = {}; counts.forEach((c, i) => { if (c) fractions[i] = c / N; });
  const out = { nVoxels: N, dims: [nx, ny, nz], voxel: h, units: vox.units || '', fractions, porePhase, porosity: counts[porePhase] / N, solidFraction: 1 - counts[porePhase] / N };
  if (N > 6e6) { out.note = 'grid above 6 million voxels: surface, percolation and chord statistics skipped'; return out; }
  const mask = new Uint8Array(N); for (let p = 0; p < N; p++) mask[p] = data[p] === porePhase ? 1 : 0;
  let facesN = 0; const chords = {}, st = [1, nx, nx * ny], dims = [nx, ny, nz];
  for (let a = 0; a < 3; a++) {
    const b = (a + 1) % 3, c = (a + 2) % 3, hist = new Uint32Array(dims[a] + 1); let sum = 0, cnt = 0;
    for (let jb = 0; jb < dims[b]; jb++) for (let jc = 0; jc < dims[c]; jc++) {
      let p = jb * st[b] + jc * st[c], run = 0, prev = -1;
      for (let i = 0; i < dims[a]; i++, p += st[a]) {
        const pore = mask[p];
        if (prev >= 0 && pore !== prev) facesN++;
        if (pore) run++; else if (run) { hist[run]++; sum += run; cnt++; run = 0; }
        prev = pore;
      }
      if (run) { hist[run]++; sum += run; cnt++; }
    }
    const hg = []; for (let i = 1; i < hist.length && hg.length < 64; i++) if (hist[i]) hg.push([i * h, hist[i]]);
    chords['xyz'[a]] = { mean: cnt ? sum / cnt * h : 0, count: cnt, histogram: hg };
  }
  out.interfaceFaces = facesN;
  out.surfaceAreaFaceCount = facesN * h * h * 2 / 3;
  out.surfaceArea = twoPointSurface(mask, nx, ny, nz) / h * (N * h ** 3);
  out.surfaceAreaNote = 'S_v = −4·dS₂/dr at zero lag from the two-point phase-change statistics along 13 lattice directions (Cauchy–Crofton weights, Richardson extrapolation over one and two lattice steps); surfaceAreaFaceCount is the first-order, axes-only value, voxel face count × h² × 2/3';
  out.ssa_per_volume = out.surfaceArea / (N * h ** 3);
  out.chords = chords;
  out.meanPoreChord = (chords.x.mean + chords.y.mean + chords.z.mean) / 3;
  out.poreSizeNote = 'mean pore chord length (line intercepts along x, y, z); for comparison, stereology gives 4φ/S_v';
  out.hydraulicDiameter = out.ssa_per_volume > 0 ? 4 * out.porosity / out.ssa_per_volume : null;
  if (surfaceOnly) return out;   // fractions, surface and chords only (no percolation search)
  out.percolates = {}; out.solidPercolates = {}; out.tortuosityGeometric = {};
  for (let a = 0; a < 3; a++) {
    const fp = flood(mask, 1, nx, ny, nz, a), fs = flood(mask, 0, nx, ny, nz, a), k = 'xyz'[a];
    out.percolates[k] = fp.percolates; out.solidPercolates[k] = fs.percolates; out.tortuosityGeometric[k] = fp.tortuosity;
  }
  out.tortuosityNote = 'geometric tortuosity = mean shortest 6-connected path between opposite faces / straight distance (city-block metric: an upper estimate)';
  return out;
}

// Blind-ant random walk in the pore space with periodic wrapping: tortuosity factor τ = D₀/D_eff = t / ⟨r²⟩ at late time.
function walkTortuosity(data, nx, ny, nz, rand, porePhase = 0) {
  const N = nx * ny * nz, pores = [];
  for (let p = 0; p < N && pores.length < 200000; p++) if (data[p] === porePhase) pores.push(p);
  if (pores.length < 8) return null;
  const W = 200, steps = max(400, max(nx, ny, nz) ** 2 >> 1), half = steps >> 1; let m1 = 0, m2 = 0;
  for (let w = 0; w < W; w++) {
    const p0 = pores[rand.int(pores.length)]; let x = p0 % nx, y = floor(p0 / nx) % ny, z = floor(p0 / (nx * ny)), ux = 0, uy = 0, uz = 0;
    for (let s = 1; s <= steps; s++) {
      const d = (rand() * 6) | 0, dx = d === 0 ? 1 : d === 1 ? -1 : 0, dy = d === 2 ? 1 : d === 3 ? -1 : 0, dz = d === 4 ? 1 : d === 5 ? -1 : 0;
      let X = x + dx, Y = y + dy, Z2 = z + dz;
      if (X < 0) X = nx - 1; else if (X === nx) X = 0; if (Y < 0) Y = ny - 1; else if (Y === ny) Y = 0; if (Z2 < 0) Z2 = nz - 1; else if (Z2 === nz) Z2 = 0;
      if (data[X + nx * (Y + ny * Z2)] === porePhase) { x = X; y = Y; z = Z2; ux += dx; uy += dy; uz += dz; }
      if (s === half) m1 += ux * ux + uy * uy + uz * uz;
    }
    m2 += ux * ux + uy * uy + uz * uz;
  }
  const slope = (m2 - m1) / W / (steps - half);
  return slope > 1e-4 ? min(1e4, 1 / slope) : null;
}

// Stochastic porous media. data: 1 = solid, 0 = pore. featureSize is in voxels (sphere / fibre diameter,
// correlation length, or foam cell size).
// light: true returns the voxel field with porosity and specific surface only (no percolation or tortuosity), for ensemble statistics.
export function porousMedium({ type = 'randomSpheres', size = [48, 48, 48], porosity = 0.5, featureSize = 8, seed = 1, voxel_um = 1, open = true, light = false } = {}) {
  let [nx, ny, nz] = [0, 1, 2].map(i => clampi(size[i] ?? 48, 4, 256)), capped = null;
  while (nx * ny * nz > 1.1e6) { nx = max(4, floor(nx * 0.9)); ny = max(4, floor(ny * 0.9)); nz = max(4, floor(nz * 0.9)); capped = 'grid reduced to keep the volume below 1.1 million voxels'; }
  const N = nx * ny * nz, data = new Uint8Array(N), rand = rng(seed ?? 1), phi = min(0.98, max(0.02, Number(porosity) || 0.5)), fs = min(max(nx, ny, nz), max(2, Number(featureSize) || 8));
  const quantileFill = (fld, solidBelow) => { // solid where field below the (1-φ) quantile, or above the φ quantile
    const s = Float32Array.from(fld).sort(), cut = solidBelow ? s[min(N - 1, max(0, round((1 - phi) * N) - 1))] : s[min(N - 1, round(phi * N))];
    for (let p = 0; p < N; p++) data[p] = solidBelow ? (fld[p] <= cut ? 1 : 0) : (fld[p] >= cut ? 1 : 0);
  };
  const stamp = (cx, cy, cz, r, wrap, dry = false) => {
    // voxels whose centre lies within r of the (sub-voxel) sphere centre: the centre is not snapped to the grid, so the digitised
    // volume and projected area of the spheres are unbiased on average
    const r2 = r * r; let added = 0;
    for (let zz = Math.ceil(cz - r); zz <= cz + r; zz++) for (let yy = Math.ceil(cy - r); yy <= cy + r; yy++) for (let xx = Math.ceil(cx - r); xx <= cx + r; xx++) {
      if ((xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2 > r2) continue;
      let x = xx, y = yy, z = zz;
      if (wrap) { x = (x % nx + nx) % nx; y = (y % ny + ny) % ny; z = (z % nz + nz) % nz; } else if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
      const p = x + nx * (y + ny * z); if (!data[p]) { if (!dry) data[p] = 1; added++; }
    }
    return added;
  };
  let model;
  if (type === 'randomSpheres') {
    model = 'Boolean model of overlapping spheres (periodic); porosity target φ = exp(−n·V_sphere), spheres added until reached';
    const r = fs / 2; let solid = 0, guard = 0;
    // the last sphere is placed only if it brings the solid fraction closer to the target than leaving it out
    while (solid < (1 - phi) * N && guard++ < 200000) { const cx = rand() * nx, cy = rand() * ny, cz = rand() * nz, add = stamp(cx, cy, cz, r, true, true); if (solid + add - (1 - phi) * N > (1 - phi) * N - solid) break; solid += stamp(cx, cy, cz, r, true); }
  } else if (type === 'fibres') {
    model = 'randomly oriented straight overlapping cylinders (isotropic directions), added until the solid fraction is reached';
    const r = max(0.75, fs / 2), Lmax = hypot(nx, ny, nz); let solid = 0, guard = 0;
    while (solid < (1 - phi) * N && guard++ < 20000) {
      const p = [rand() * nx, rand() * ny, rand() * nz], cz = 2 * rand() - 1, a = 2 * PI * rand(), sr = sqrt(1 - cz * cz), d = [sr * cos(a), sr * sin(a), cz];
      for (let t = -Lmax; t <= Lmax; t += max(0.5, r * 0.6)) { const x = p[0] + t * d[0], y = p[1] + t * d[1], z = p[2] + t * d[2]; if (x < -r || y < -r || z < -r || x > nx + r || y > ny + r || z > nz + r) continue; solid += stamp(x, y, z, r, false); }
    }
  } else if (type === 'voronoiFoam') {
    // one jittered seed per bin of edge ≈ featureSize: the three nearest seeds are searched in the 27 surrounding bins
    const bx = max(1, round(nx / fs)), by = max(1, round(ny / fs)), bz = max(1, round(nz / fs)), ns = bx * by * bz, local = bx >= 3 && by >= 3 && bz >= 3;
    model = `periodic Voronoi tessellation of ${ns} jittered-grid seeds; ${open ? 'open-cell: struts where d₃ − d₁ is small (Plateau borders)' : 'closed-cell: walls where d₂ − d₁ is small'}; threshold set by quantile to meet the porosity`;
    const sx = new Float32Array(ns), sy = new Float32Array(ns), sz = new Float32Array(ns), fld = new Float32Array(N);
    for (let k = 0, i = 0; k < bz; k++) for (let j = 0; j < by; j++) for (let ii = 0; ii < bx; ii++, i++) { sx[i] = (ii + rand()) * nx / bx; sy[i] = (j + rand()) * ny / by; sz[i] = (k + rand()) * nz / bz; }
    if (!local && ns * N > 3e7) capped = (capped ? capped + '; ' : '') + 'coarse foam on a large grid: generation is slow because fewer than three cells fit along an axis';
    const cand = new Int32Array(local ? 27 : ns); if (!local) for (let i = 0; i < ns; i++) cand[i] = i;
    for (let z = 0, p = 0; z < nz; z++) for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++, p++) {
      if (local && (x === 0 || (x * bx / nx | 0) !== ((x - 1) * bx / nx | 0))) {
        const ix = x * bx / nx | 0, iy = y * by / ny | 0, iz = z * bz / nz | 0; let c = 0;
        for (let dk = -1; dk <= 1; dk++) for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) cand[c++] = (ix + di + bx) % bx + bx * ((iy + dj + by) % by + by * ((iz + dk + bz) % bz));
      }
      let d1 = 1e30, d2 = 1e30, d3 = 1e30;
      for (let c = 0; c < cand.length; c++) {
        const i = cand[c]; let dx = abs(x - sx[i]), dy = abs(y - sy[i]), dz = abs(z - sz[i]);
        if (dx > nx / 2) dx = nx - dx; if (dy > ny / 2) dy = ny - dy; if (dz > nz / 2) dz = nz - dz;
        const d = dx * dx + dy * dy + dz * dz;
        if (d < d1) { d3 = d2; d2 = d1; d1 = d; } else if (d < d2) { d3 = d2; d2 = d; } else if (d < d3) d3 = d;
      }
      fld[p] = sqrt(open ? d3 : d2) - sqrt(d1);
    }
    quantileFill(fld, true);
  } else if (type === 'gaussianField') {
    model = 'level-cut smoothed random field (white noise, three periodic box-filter passes ≈ Gaussian kernel); threshold set by quantile to meet the porosity';
    let a = new Float32Array(N), b = new Float32Array(N); const R = max(1, round(fs / 4)), st = [1, nx, nx * ny], dims = [nx, ny, nz];
    for (let p = 0; p < N; p++) a[p] = rand() - 0.5;
    for (let pass = 0; pass < 3; pass++) for (let ax = 0; ax < 3; ax++) {
      const n = dims[ax], s = st[ax], o1 = dims[(ax + 1) % 3], o2 = dims[(ax + 2) % 3], s1 = st[(ax + 1) % 3], s2 = st[(ax + 2) % 3], r = min(R, (n - 1) >> 1);
      for (let i1 = 0; i1 < o1; i1++) for (let i2 = 0; i2 < o2; i2++) {
        const base = i1 * s1 + i2 * s2; let acc = 0;
        for (let q = -r; q <= r; q++) acc += a[base + ((q + n) % n) * s];
        for (let i = 0; i < n; i++) { b[base + i * s] = acc; acc += a[base + ((i + r + 1) % n) * s] - a[base + ((i - r + n) % n) * s]; }
      }
      [a, b] = [b, a];
    }
    quantileFill(a, false);
  } else throw new Error(`Unknown porous-medium type "${type}". Use randomSpheres, gaussianField, fibres or voronoiFoam.`);
  const vox = { kind: 'voxels', nx, ny, nz, data, voxel: voxel_um, units: 'µm', type, model, seed: seed ?? 1 };
  if (light) { const q = voxelStats(vox, { surfaceOnly: true }); return Object.assign(vox, { porosity: q.porosity, targetPorosity: phi, ssa: q.ssa_per_volume, light: true }); }
  const s = voxelStats(vox);
  Object.assign(vox, { porosity: s.porosity, targetPorosity: phi, ssa: s.ssa_per_volume, ssaUnits: '1/µm (interfacial area per bulk volume, two-point estimator extrapolated to zero lag)', ssaFaceCount: s.surfaceAreaFaceCount / (N * voxel_um ** 3), percolates: s.percolates, solidPercolates: s.solidPercolates, tortuosityGeometric: s.tortuosityGeometric, meanPoreChord: s.meanPoreChord, hydraulicDiameter: s.hydraulicDiameter });
  vox.tortuosity = (s.percolates.x || s.percolates.y || s.percolates.z) ? walkTortuosity(data, nx, ny, nz, rand) : null;
  vox.tortuosityNote = 'tortuosity factor τ = D₀/D_eff from a pore-space lattice random walk with periodic wrapping (200 walkers; statistical scatter of roughly ±15 %); null when the pore space does not percolate';
  if (capped) vox.capped = capped;
  return vox;
}
