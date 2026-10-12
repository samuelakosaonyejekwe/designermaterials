// Geometry engine 1 — atoms, molecules and crystals.
// Crystal library (conventional cells), builders (supercell, nanoparticle, nanosheet, nanotube, core–shell, slab, defects),
// analysis (bonds, coordination, density, RDF, d-spacings, kinematic powder XRD) and text writers/parsers for the
// common atomistic file formats. Lengths in Å unless a name says _nm. Pure functions, no DOM / network / fs.
import { N_A, rng } from '../core/num.js';
import { ELEMENTS, SYMBOLS } from '../core/chem.js';

const D2R = Math.PI / 180, BOHR = 0.529177210903;
export const MAX_ATOMS = 60000;            // default cap for generated clusters (keeps phones responsive)

// ───────────────────────────── small vector / cell helpers ─────────────────────────────
const det3 = m => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
function inv3(m) {
  const d = det3(m); if (!(Math.abs(d) > 1e-12) || !Number.isFinite(d)) throw new Error('Degenerate cell: lattice vectors are coplanar or not finite (zero volume)');
  const [a, b, c] = m;
  return [[(b[1] * c[2] - b[2] * c[1]) / d, (a[2] * c[1] - a[1] * c[2]) / d, (a[1] * b[2] - a[2] * b[1]) / d],
  [(b[2] * c[0] - b[0] * c[2]) / d, (a[0] * c[2] - a[2] * c[0]) / d, (a[2] * b[0] - a[0] * b[2]) / d],
  [(b[0] * c[1] - b[1] * c[0]) / d, (a[1] * c[0] - a[0] * c[1]) / d, (a[0] * b[1] - a[1] * b[0]) / d]];
}
const rowMul = (v, M) => [v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0], v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1], v[0] * M[0][2] + v[1] * M[1][2] + v[2] * M[2][2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = a => Math.hypot(a[0], a[1], a[2]);
const wrap1 = x => { let w = x - Math.floor(x); if (w > 1 - 1e-9) w = 0; return w; };
const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a; };
const r6 = v => { const x = Math.round(v * 1e6) / 1e6; return x === 0 ? 0 : x; };
const r4 = v => { const x = Math.round(v * 1e4) / 1e4; return x === 0 ? 0 : x; };
const massOf = el => ELEMENTS[el]?.mass ?? 0;
const radOf = el => ELEMENTS[el]?.radius ?? 1.5;
const enOf = el => ELEMENTS[el]?.en ?? 0;

/** Lattice vectors (rows a, b, c) in the standard orientation: a ∥ x, b in the xy-plane. */
export function latticeVectors({ a, b = a, c = a, alpha = 90, beta = 90, gamma = 90 }) {
  if (![a, b, c, alpha, beta, gamma].every(v => Number.isFinite(v) && v > 0)) throw new Error('latticeVectors: a, b, c and the three angles must be positive finite numbers');
  const cz = v => Math.abs(v) < 1e-12 ? 0 : v;
  const ca = cz(Math.cos(alpha * D2R)), cb = cz(Math.cos(beta * D2R)), cg = cz(Math.cos(gamma * D2R)), sg = Math.sin(gamma * D2R);
  const cy = (ca - cb * cg) / sg, z2 = 1 - cb * cb - cy * cy;
  if (!(z2 > 1e-12)) throw new Error(`latticeVectors: angles α=${alpha}, β=${beta}, γ=${gamma} do not form a valid cell`);
  return [[a, 0, 0], [cz(b * cg), b * sg, 0], [cz(c * cb), cz(c * cy), c * Math.sqrt(z2)]];
}
/** Cell parameters from lattice vectors. */
export function cellParams(cell) {
  const [A, B, C] = cell, a = len(A), b = len(B), c = len(C), ang = (u, v, lu, lv) => Math.acos(clamp(dot(u, v) / (lu * lv), -1, 1)) / D2R;
  return { a, b, c, alpha: ang(B, C, b, c), beta: ang(A, C, a, c), gamma: ang(A, B, a, b), volume: Math.abs(det3(cell)) };
}

// ───────────────────────────── crystal library ─────────────────────────────
const P = [[0, 0, 0]], F = [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], I = [[0, 0, 0], [.5, .5, .5]], RH = [[0, 0, 0], [2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]];
function xtal(name, formula, system, spacegroup, lat, sites, centring = P, notes = '') {
  const [a, b = a, c = a, alpha = 90, beta = 90, gamma = 90] = lat, basis = [];
  for (const [el, pts] of sites) for (const p of pts) for (const t of centring) {
    const f = [wrap1(p[0] + t[0]), wrap1(p[1] + t[1]), wrap1(p[2] + t[2])];
    if (!basis.some(q => q.el === el && [0, 1, 2].every(k => { const d = Math.abs(q.f[k] - f[k]); return Math.min(d, 1 - d) < 1e-5; }))) basis.push({ el, f });
  }
  return { name, formula, system, spacegroup, a, b, c, alpha, beta, gamma, basis, notes };
}
const T3 = 1 / 3, T23 = 2 / 3, Q = .25, H = .5;
const wurtzite = (A, B, u) => [[A, [[T3, T23, 0], [T23, T3, H]]], [B, [[T3, T23, u], [T23, T3, H + u]]]];
const zincblende = (A, B) => [[A, [[0, 0, 0]]], [B, [[Q, Q, Q]]]];
const diamond = A => [[A, [[0, 0, 0], [Q, Q, Q]]]];
const fccM = (el, name, a) => xtal(name, el, 'cubic', 'Fm-3m (225)', [a], [[el, [[0, 0, 0]]]], F, 'Face-centred cubic metal, room-temperature lattice parameter.');
const HEX = (a, c) => [a, a, c, 90, 90, 120];
const spinelO = x => [[x, x, x], [x, -x, -x], [-x, x, -x], [-x, -x, x], [Q - x, Q - x, Q - x], [Q - x, Q + x, Q + x], [Q + x, Q - x, Q + x], [Q + x, Q + x, Q - x]];

export const CRYSTALS = {
  ZnO_wurtzite: xtal('Zinc oxide (wurtzite)', 'ZnO', 'hexagonal', 'P6_3mc (186)', HEX(3.2495, 5.2069), wurtzite('Zn', 'O', 0.3819), P, 'Internal parameter u = 0.3819; polar along c, nanorods grow along [0001].'),
  ZnO_zincblende: xtal('Zinc oxide (zinc blende)', 'ZnO', 'cubic', 'F-43m (216)', [4.62], zincblende('Zn', 'O'), F, 'Metastable cubic polymorph stabilised on cubic substrates; a ≈ 4.62 Å (thin-film value).'),
  TiO2_rutile: xtal('Titanium dioxide (rutile)', 'TiO2', 'tetragonal', 'P4_2/mnm (136)', [4.5937, 4.5937, 2.9587], [['Ti', [[0, 0, 0], [H, H, H]]], ['O', [[0.3048, 0.3048, 0], [-0.3048, -0.3048, 0], [H + 0.3048, H - 0.3048, H], [H - 0.3048, H + 0.3048, H]]]], P, 'Oxygen parameter u = 0.3048.'),
  TiO2_anatase: xtal('Titanium dioxide (anatase)', 'TiO2', 'tetragonal', 'I4_1/amd (141)', [3.7845, 3.7845, 9.5143], [['Ti', [[0, 0, 0], [0, H, Q]]], ['O', [[0, 0, 0.2081], [0, H, 0.4581], [H, 0, 0.5419], [H, H, 0.2919]]]], I, 'Origin choice 1, oxygen parameter z = 0.2081.'),
  Si: xtal('Silicon', 'Si', 'cubic', 'Fd-3m (227)', [5.4310], diamond('Si'), F, 'Diamond cubic.'),
  Ge: xtal('Germanium', 'Ge', 'cubic', 'Fd-3m (227)', [5.6579], diamond('Ge'), F, 'Diamond cubic.'),
  C_diamond: xtal('Diamond', 'C', 'cubic', 'Fd-3m (227)', [3.5669], diamond('C'), F, 'Diamond cubic carbon.'),
  graphite: xtal('Graphite (2H, AB stacking)', 'C', 'hexagonal', 'P6_3/mmc (194)', HEX(2.464, 6.711), [['C', [[0, 0, Q], [0, 0, .75], [T3, T23, Q], [T23, T3, .75]]]], P, 'Bernal AB stacking, interlayer spacing c/2 = 3.356 Å.'),
  graphene: xtal('Graphene (monolayer, 20 Å vacuum)', 'C', 'hexagonal', 'P6/mmm (191)', HEX(2.46, 20), [['C', [[T3, T23, H], [T23, T3, H]]]], P, 'Single layer in a periodic box; c = 20 Å is vacuum, not a lattice parameter.'),
  hBN: xtal('Hexagonal boron nitride', 'BN', 'hexagonal', 'P6_3/mmc (194)', HEX(2.5040, 6.6612), [['B', [[T3, T23, Q], [T23, T3, .75]]], ['N', [[T23, T3, Q], [T3, T23, .75]]]], P, "AA' stacking: B sits above N in adjacent layers."),
  MoS2_2H: xtal('Molybdenum disulfide (2H)', 'MoS2', 'hexagonal', 'P6_3/mmc (194)', HEX(3.1604, 12.295), [['Mo', [[T3, T23, Q], [T23, T3, .75]]], ['S', [[T3, T23, 0.621], [T23, T3, 0.121], [T23, T3, 0.379], [T3, T23, 0.879]]]], P, 'Trigonal-prismatic Mo; sulfur parameter z = 0.621.'),
  GaAs: xtal('Gallium arsenide', 'GaAs', 'cubic', 'F-43m (216)', [5.6533], zincblende('Ga', 'As'), F, 'Zinc blende.'),
  GaN_wurtzite: xtal('Gallium nitride (wurtzite)', 'GaN', 'hexagonal', 'P6_3mc (186)', HEX(3.189, 5.185), wurtzite('Ga', 'N', 0.377), P, 'Internal parameter u = 0.377.'),
  CdSe_wurtzite: xtal('Cadmium selenide (wurtzite)', 'CdSe', 'hexagonal', 'P6_3mc (186)', HEX(4.2999, 7.0109), wurtzite('Cd', 'Se', 0.3767), P, 'Internal parameter u ≈ 0.377; common quantum-dot phase.'),
  CdS_zincblende: xtal('Cadmium sulfide (zinc blende)', 'CdS', 'cubic', 'F-43m (216)', [5.832], zincblende('Cd', 'S'), F, 'Hawleyite.'),
  ZnS_zincblende: xtal('Zinc sulfide (sphalerite)', 'ZnS', 'cubic', 'F-43m (216)', [5.4093], zincblende('Zn', 'S'), F, 'Sphalerite.'),
  Au: fccM('Au', 'Gold', 4.0782), Ag: fccM('Ag', 'Silver', 4.0853), Cu: fccM('Cu', 'Copper', 3.6149), Pt: fccM('Pt', 'Platinum', 3.9242),
  Pd: fccM('Pd', 'Palladium', 3.8907), Ni: fccM('Ni', 'Nickel', 3.5240), Al: fccM('Al', 'Aluminium', 4.0495),
  Fe_bcc: xtal('Iron (α, bcc)', 'Fe', 'cubic', 'Im-3m (229)', [2.8665], [['Fe', [[0, 0, 0]]]], I, 'Body-centred cubic ferrite.'),
  W: xtal('Tungsten', 'W', 'cubic', 'Im-3m (229)', [3.1652], [['W', [[0, 0, 0]]]], I, 'Body-centred cubic.'),
  NaCl: xtal('Sodium chloride (rock salt)', 'NaCl', 'cubic', 'Fm-3m (225)', [5.6402], [['Na', [[0, 0, 0]]], ['Cl', [[H, H, H]]]], F, 'Rock salt.'),
  MgO: xtal('Magnesium oxide (periclase)', 'MgO', 'cubic', 'Fm-3m (225)', [4.2112], [['Mg', [[0, 0, 0]]], ['O', [[H, H, H]]]], F, 'Rock salt.'),
  CeO2: xtal('Cerium dioxide (fluorite)', 'CeO2', 'cubic', 'Fm-3m (225)', [5.4113], [['Ce', [[0, 0, 0]]], ['O', [[Q, Q, Q], [.75, .75, .75]]]], F, 'Fluorite; oxygen fills all tetrahedral sites.'),
  SrTiO3: xtal('Strontium titanate (cubic perovskite)', 'SrTiO3', 'cubic', 'Pm-3m (221)', [3.905], [['Sr', [[0, 0, 0]]], ['Ti', [[H, H, H]]], ['O', [[H, H, 0], [H, 0, H], [0, H, H]]]], P, 'Ideal cubic perovskite.'),
  BaTiO3_cubic: xtal('Barium titanate (cubic, paraelectric)', 'BaTiO3', 'cubic', 'Pm-3m (221)', [4.009], [['Ba', [[0, 0, 0]]], ['Ti', [[H, H, H]]], ['O', [[H, H, 0], [H, 0, H], [0, H, H]]]], P, 'Paraelectric phase above the Curie point (≈ 393 K); tetragonal and ferroelectric at room temperature.'),
  CsCl: xtal('Caesium chloride', 'CsCl', 'cubic', 'Pm-3m (221)', [4.123], [['Cs', [[0, 0, 0]]], ['Cl', [[H, H, H]]]], P, 'B2 structure.'),
  Fe3O4: xtal('Magnetite (inverse spinel)', 'Fe3O4', 'cubic', 'Fd-3m (227)', [8.3967], [['Fe', [[0, 0, 0], [Q, Q, Q], [.625, .625, .625], [.625, .875, .875], [.875, .625, .875], [.875, .875, .625]]], ['O', spinelO(0.3799)]], F, 'Origin choice 1; 8 tetrahedral + 16 octahedral Fe, oxygen parameter x = 0.3799 (u = 0.2549 in origin choice 2). Fe²⁺/Fe³⁺ ordering on the octahedral site is not distinguished.'),
  SiO2_cristobalite: xtal('β-Cristobalite (ideal)', 'SiO2', 'cubic', 'Fd-3m (227)', [7.16], [['Si', [[0, 0, 0], [Q, Q, Q]]], ['O', [[.125, .125, .125], [.125, .375, .375], [.375, .125, .375], [.375, .375, .125]]]], F, 'Idealised high-temperature structure with linear Si–O–Si; the real phase is dynamically disordered, so the Si–O distance (1.55 Å) is shorter than the true bond.'),
  Al2O3_corundum: xtal('α-Alumina (corundum)', 'Al2O3', 'trigonal', 'R-3c (167)', HEX(4.7589, 12.991), [['Al', [[0, 0, 0.3522], [0, 0, 0.1478], [0, 0, -0.3522], [0, 0, 0.8522]]], ['O', [[0.3064, 0, Q], [0, 0.3064, Q], [-0.3064, -0.3064, Q], [-0.3064, 0, .75], [0, -0.3064, .75], [0.3064, 0.3064, .75]]]], RH, 'Hexagonal setting, 6 formula units; Al z = 0.3522, O x = 0.3064.'),
};

function resolveDef(keyOrDef) {
  if (typeof keyOrDef === 'string') {
    const k = CRYSTALS[keyOrDef] ? keyOrDef : Object.keys(CRYSTALS).find(q => q.toLowerCase() === keyOrDef.toLowerCase());
    if (!k) throw new Error(`Unknown crystal "${keyOrDef}". Available: ${Object.keys(CRYSTALS).join(', ')}`);
    return CRYSTALS[k];
  }
  if (!keyOrDef || typeof keyOrDef !== 'object' || !Array.isArray(keyOrDef.basis) || !keyOrDef.basis.length) throw new Error('Crystal definition needs lattice parameters and a non-empty basis:[{el, f:[u,v,w]}]');
  return keyOrDef;
}
/** Conventional unit cell of a library crystal (key) or a definition {a,b,c,alpha,beta,gamma,basis}. */
export function buildCell(keyOrDef) {
  const d = resolveDef(keyOrDef), cell = latticeVectors(d);
  const atoms = d.basis.map((s, i) => {
    if (!s || !Array.isArray(s.f) || s.f.length !== 3 || !s.f.every(Number.isFinite)) throw new Error(`Crystal definition: basis entry ${i + 1} needs el and f:[u,v,w]`);
    const f = s.f.map(wrap1), r = rowMul(f, cell); const a = { el: s.el, x: r6(r[0]), y: r6(r[1]), z: r6(r[2]), f };
    if (s.occ !== undefined && s.occ < 1) a.occ = s.occ; return a;
  });
  return { kind: 'atoms', atoms, cell, title: d.name || d.formula || 'crystal', formula: d.formula || formulaOf(composition(atoms)), system: d.system || null, spacegroup: d.spacegroup || null };
}
function asPeriodic(x, who) {
  if (x && typeof x === 'object' && x.kind === 'atoms') {
    if (!x.cell) throw new Error(`${who}: needs a periodic structure (a crystal key, a definition, or atoms with a cell)`);
    if (!x.atoms?.length) throw new Error(`${who}: structure has no atoms`);
    return x;
  }
  return buildCell(x);
}
const titleOf = x => typeof x === 'string' ? (CRYSTALS[x]?.formula || x) : (x?.formula || x?.title || x?.name || 'crystal');

/** Replicate a periodic structure nx × ny × nz. */
export function supercell(struct, reps = [1, 1, 1]) {
  if (!struct?.cell) throw new Error('supercell: structure has no cell');
  const n = [0, 1, 2].map(k => Math.max(1, Math.round(reps[k] ?? 1))), M = struct.cell, Mi = inv3(M), atoms = [];
  if (struct.atoms.length * n[0] * n[1] * n[2] > 2e6) throw new Error('supercell: more than 2 million atoms requested; reduce the repeat counts');
  const fr = struct.atoms.map(a => rowMul([a.x, a.y, a.z], Mi));
  for (let i = 0; i < n[0]; i++) for (let j = 0; j < n[1]; j++) for (let k = 0; k < n[2]; k++) struct.atoms.forEach((a, q) => {
    const g = [fr[q][0] + i, fr[q][1] + j, fr[q][2] + k], r = rowMul(g, M);
    atoms.push({ ...a, x: r6(r[0]), y: r6(r[1]), z: r6(r[2]), f: [g[0] / n[0], g[1] / n[1], g[2] / n[2]] });
  });
  return { ...struct, kind: 'atoms', atoms, cell: M.map((v, k) => v.map(c => c * n[k])), title: `${struct.title || 'structure'} ${n.join('×')} supercell`, repeats: n };
}

// ───────────────────────────── neighbour search ─────────────────────────────
// Calls cb(i, j, d) for every ordered pair (centre i, neighbour j or a periodic image of j) with 0 < d ≤ rcut.
function neighbourScan(struct, rcut, cb, stride = 1) {
  const at = struct.atoms, n = at.length; if (!n) return;
  const px = new Float64Array(n), py = new Float64Array(n), pz = new Float64Array(n);
  let ex, ey, ez, ei;
  if (struct.cell) {
    const M = struct.cell, Mi = inv3(M), del = [0, 1, 2].map(k => rcut * Math.hypot(Mi[0][k], Mi[1][k], Mi[2][k])), nm = del.map(Math.ceil);
    const X = [], Y = [], Z = [], J = [];
    for (let i = 0; i < n; i++) {
      const f = rowMul([at[i].x, at[i].y, at[i].z], Mi).map(wrap1), r = rowMul(f, M); px[i] = r[0]; py[i] = r[1]; pz[i] = r[2];
      for (let a = -nm[0]; a <= nm[0]; a++) { const g0 = f[0] + a; if (g0 < -del[0] || g0 > 1 + del[0]) continue;
        for (let b = -nm[1]; b <= nm[1]; b++) { const g1 = f[1] + b; if (g1 < -del[1] || g1 > 1 + del[1]) continue;
          for (let c = -nm[2]; c <= nm[2]; c++) { const g2 = f[2] + c; if (g2 < -del[2] || g2 > 1 + del[2]) continue;
            X.push(g0 * M[0][0] + g1 * M[1][0] + g2 * M[2][0]); Y.push(g0 * M[0][1] + g1 * M[1][1] + g2 * M[2][1]); Z.push(g0 * M[0][2] + g1 * M[1][2] + g2 * M[2][2]); J.push(i);
          } } }
    }
    ex = X; ey = Y; ez = Z; ei = J;
  } else { for (let i = 0; i < n; i++) { px[i] = at[i].x; py[i] = at[i].y; pz[i] = at[i].z; } ex = px; ey = py; ez = pz; ei = null; }
  const m = ex.length; let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
  for (let q = 0; q < m; q++) { if (ex[q] < x0) x0 = ex[q]; if (ex[q] > x1) x1 = ex[q]; if (ey[q] < y0) y0 = ey[q]; if (ey[q] > y1) y1 = ey[q]; if (ez[q] < z0) z0 = ez[q]; if (ez[q] > z1) z1 = ez[q]; }
  const nb = e => clamp(Math.floor(e / rcut), 1, 160), nx = nb(x1 - x0), ny = nb(y1 - y0), nz = nb(z1 - z0);
  const sx = nx / Math.max(x1 - x0, 1e-9), sy = ny / Math.max(y1 - y0, 1e-9), sz = nz / Math.max(z1 - z0, 1e-9);
  const bin = (v, v0, s, nn) => Math.min(nn - 1, Math.max(0, Math.floor((v - v0) * s)));
  const head = new Int32Array(nx * ny * nz).fill(-1), next = new Int32Array(m);
  for (let q = 0; q < m; q++) { const b = bin(ex[q], x0, sx, nx) + nx * (bin(ey[q], y0, sy, ny) + ny * bin(ez[q], z0, sz, nz)); next[q] = head[b]; head[b] = q; }
  const r2 = rcut * rcut;
  for (let i = 0; i < n; i += stride) {
    const bx = bin(px[i], x0, sx, nx), by = bin(py[i], y0, sy, ny), bz = bin(pz[i], z0, sz, nz);
    for (let c = Math.max(0, bz - 1); c <= Math.min(nz - 1, bz + 1); c++) for (let b = Math.max(0, by - 1); b <= Math.min(ny - 1, by + 1); b++) for (let a = Math.max(0, bx - 1); a <= Math.min(nx - 1, bx + 1); a++) {
      for (let q = head[a + nx * (b + ny * c)]; q >= 0; q = next[q]) {
        const dx = ex[q] - px[i], dy = ey[q] - py[i], dz = ez[q] - pz[i], d2 = dx * dx + dy * dy + dz * dz;
        if (d2 <= r2 && d2 > 1e-8) cb(i, ei ? ei[q] : q, Math.sqrt(d2));
      }
    }
  }
}

// Bond perception (heuristic, documented in analyse()).
function bondGraph(struct, wantList = false) {
  const at = struct.atoms, n = at.length, rad = at.map(a => radOf(a.el)), en = at.map(a => enOf(a.el));
  let rmax = 0; for (const r of rad) if (r > rmax) rmax = r;
  const rcut = Math.min(2.5 * rmax, 7), nbJ = Array.from({ length: n }, () => []), nbD = Array.from({ length: n }, () => []);
  let nn = Infinity;
  neighbourScan(struct, rcut, (i, j, d) => { if (d < nn) nn = d; if (d <= 1.25 * (rad[i] + rad[j])) { nbJ[i].push(j); nbD[i].push(d); } });
  // shortest contact to each neighbouring element, and to any clearly more electronegative neighbour
  const dEl = new Array(n), dAn = new Float64Array(n).fill(Infinity);
  for (let i = 0; i < n; i++) { const m = {}; for (let q = 0; q < nbJ[i].length; q++) { const j = nbJ[i][q], d = nbD[i][q], e = at[j].el; if (!(m[e] <= d)) m[e] = d; if (en[j] >= en[i] + 0.3 && d < dAn[i]) dAn[i] = d; } dEl[i] = m; }
  const cn = new Int32Array(n), pairs = {}, list = [];
  for (let i = 0; i < n; i++) for (let q = 0; q < nbJ[i].length; q++) {
    const j = nbJ[i][q], d = nbD[i][q];
    if (!(d <= 1.14 * dEl[i][at[j].el] || d <= 1.14 * dEl[j][at[i].el])) continue;                 // first shell only
    if (en[i] < 2.2 && en[j] < 2.2 && dAn[i] < d / 1.1 && dAn[j] < d / 1.1) continue;              // cation–cation contact bridged by anions
    cn[i]++;
    const A = at[i].el, B = at[j].el, key = enOf(A) < enOf(B) || (enOf(A) === enOf(B) && A <= B) ? `${A}-${B}` : `${B}-${A}`;
    const p = pairs[key] || (pairs[key] = { sum: 0, count: 0, min: Infinity, max: 0 }); p.sum += d; p.count++; if (d < p.min) p.min = d; if (d > p.max) p.max = d;
    if (wantList && i < j) list.push([i, j, d]);
  }
  return { cn, pairs, nn: Number.isFinite(nn) ? nn : null, list };
}

function composition(atoms) { const c = {}; for (const a of atoms) c[a.el] = (c[a.el] || 0) + 1; return c; }
function formulaOf(comp) {
  const els = Object.keys(comp); if (!els.length) return '';
  const counts = els.map(e => comp[e]); let g = counts.every(Number.isInteger) ? counts.reduce(gcd) : 1; if (!g) g = 1;
  if (comp.C && comp.H) els.sort((a, b) => (a === 'C' ? -2 : a === 'H' ? -1 : 0) - (b === 'C' ? -2 : b === 'H' ? -1 : 0) || (a < b ? -1 : 1));   // Hill order
  else els.sort((a, b) => enOf(a) - enOf(b) || (a < b ? -1 : 1));                                                                          // electropositive first
  return els.map(e => { const k = comp[e] / g; return e + (k === 1 ? '' : k); }).join('');
}

/**
 * Structure analysis. Bonds use a geometric heuristic: a contact is a bond when d ≤ 1.25 (r_cov,i + r_cov,j), it belongs to the
 * first shell of that element pair (d ≤ 1.14 × the shortest such contact of either atom), and it is not a cation–cation contact
 * bridged by closer anions. Suitable for crystals, clusters and simple molecules; it is not a chemical-perception engine.
 */
export function analyse(struct) {
  if (!struct?.atoms?.length) throw new Error('analyse: structure has no atoms');
  const at = struct.atoms, n = at.length, comp = composition(at), bg = bondGraph(struct);
  let M = 0, cx = 0, cy = 0, cz = 0;
  for (const a of at) { const m = massOf(a.el) || 1; M += m; cx += m * a.x; cy += m * a.y; cz += m * a.z; }
  cx /= M; cy /= M; cz /= M; let rg = 0; for (const a of at) rg += (massOf(a.el) || 1) * ((a.x - cx) ** 2 + (a.y - cy) ** 2 + (a.z - cz) ** 2);
  rg = Math.sqrt(rg / M);
  let volume = null, density = null, mass = 0; for (const a of at) mass += massOf(a.el) * (a.occ ?? 1);
  if (struct.cell) { volume = Math.abs(det3(struct.cell)); density = volume > 0 ? mass * 1e-3 / N_A / (volume * 1e-30) : null; }
  let nnA = bg.nn;
  if (nnA === null && n > 1 && n <= 3000 && !struct.cell) { let m = Infinity; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d = Math.hypot(at[i].x - at[j].x, at[i].y - at[j].y, at[i].z - at[j].z); if (d < m) m = d; } nnA = m; }
  const cnSum = {}, cnMax = {}; at.forEach((a, i) => { cnSum[a.el] = (cnSum[a.el] || 0) + bg.cn[i]; if (!(cnMax[a.el] >= bg.cn[i])) cnMax[a.el] = bg.cn[i]; });
  const coordination = {}; for (const e in comp) coordination[e] = cnSum[e] / comp[e];
  // Under-coordination: against the bulk coordination of the parent site when known (carved particles), otherwise against the
  // highest coordination found for that element. A periodic cell without a vacuum gap (> 5 Å) has no surface.
  const ref = Array.isArray(struct.bulkCN) ? struct.bulkCN : null; let under = 0, hasSurface = true;
  if (struct.cell && !ref) { const Mi = inv3(struct.cell); hasSurface = false; for (let k = 0; k < 3 && !hasSurface; k++) { const h = 1 / Math.hypot(Mi[0][k], Mi[1][k], Mi[2][k]), f = at.map(a => wrap1(a.x * Mi[0][k] + a.y * Mi[1][k] + a.z * Mi[2][k])).sort((p, q) => p - q); let g = f[0] + 1 - f[f.length - 1]; for (let i = 1; i < f.length; i++) if (f[i] - f[i - 1] > g) g = f[i] - f[i - 1]; if (g * h > 5) hasSurface = true; } }
  if (hasSurface) at.forEach((a, i) => { const full = ref && a.b !== undefined && ref[a.b] !== undefined ? ref[a.b] : cnMax[a.el]; if (bg.cn[i] < full) under++; });
  const bonds = Object.keys(bg.pairs).map(k => { const p = bg.pairs[k]; return { pair: k, length_A: p.sum / p.count, count: p.count / 2, min_A: p.min, max_A: p.max }; }).sort((a, b) => a.length_A - b.length_A);
  return { nAtoms: n, composition: comp, formula: formulaOf(comp), density_kg_m3: density, volume_A3: volume, nearestNeighbour_A: nnA, bonds, coordination, surfaceFraction: under / n, radiusOfGyration_A: rg, molarMass_g_mol: mass, periodic: !!struct.cell };
}

/** Radial distribution function g(r). Periodic structures use the cell density; finite clusters use an equivalent-sphere density estimate. */
export function rdf(struct, { rmax_A = 10, dr_A = 0.05 } = {}) {
  if (!struct?.atoms?.length) throw new Error('rdf: structure has no atoms');
  const n = struct.atoms.length, rmax = clamp(rmax_A > 0 ? rmax_A : 10, 1, 30); let dr = dr_A > 0 ? dr_A : 0.05; if (rmax / dr > 3000) dr = rmax / 3000;
  const nb = Math.max(1, Math.round(rmax / dr)), counts = new Float64Array(nb);
  let rho, normalisation;
  if (struct.cell) { rho = n / Math.abs(det3(struct.cell)); normalisation = 'periodic cell density'; }
  else {
    let cx = 0, cy = 0, cz = 0; for (const a of struct.atoms) { cx += a.x; cy += a.y; cz += a.z; } cx /= n; cy /= n; cz /= n;
    let s = 0; for (const a of struct.atoms) s += (a.x - cx) ** 2 + (a.y - cy) ** 2 + (a.z - cz) ** 2;
    const R = Math.max(Math.sqrt(5 / 3 * s / n), 1); rho = n / (4 / 3 * Math.PI * R ** 3);
    normalisation = 'finite cluster: density of the equivalent sphere (R = √(5/3)·Rg); g(r) decays with r because of the missing neighbours at the surface';
  }
  const stride = Math.max(1, Math.ceil(n * rho * 27 * rmax ** 3 / 6e7));
  let centres = 0; for (let i = 0; i < n; i += stride) centres++;
  neighbourScan(struct, rmax, (i, j, d) => { const b = Math.floor(d / dr); if (b < nb) counts[b]++; }, stride);
  const r = [], g = [], nPer = [], cn = []; let run = 0;
  for (let b = 0; b < nb; b++) { const r1 = b * dr, r2 = r1 + dr, c = counts[b] / centres; run += c; r.push(r6(r1 + dr / 2)); nPer.push(c); cn.push(run); g.push(c / (rho * 4 / 3 * Math.PI * (r2 ** 3 - r1 ** 3))); }
  return { r, g, n: nPer, cn, density_A3: rho, dr_A: dr, normalisation, sampledCentres: centres };
}

// ───────────────────────────── builders ─────────────────────────────
function carve(s, inside, rb) {
  const M = s.cell, Mi = inv3(M), nm = [0, 1, 2].map(k => Math.ceil(rb * Math.hypot(Mi[0][k], Mi[1][k], Mi[2][k])) + 1);
  const base = s.atoms.map(a => { const f = rowMul([a.x, a.y, a.z], Mi).map(wrap1); return rowMul(f, M); });
  const diag = len(M[0]) + len(M[1]) + len(M[2]), lim = (rb + diag) ** 2, out = [];
  for (let i = -nm[0]; i <= nm[0]; i++) for (let j = -nm[1]; j <= nm[1]; j++) for (let k = -nm[2]; k <= nm[2]; k++) {
    const ox = i * M[0][0] + j * M[1][0] + k * M[2][0], oy = i * M[0][1] + j * M[1][1] + k * M[2][1], oz = i * M[0][2] + j * M[1][2] + k * M[2][2];
    if (ox * ox + oy * oy + oz * oz > lim) continue;
    for (let q = 0; q < base.length; q++) { const x = ox + base[q][0], y = oy + base[q][1], z = oz + base[q][2]; if (inside(x, y, z)) out.push({ el: s.atoms[q].el, x: r4(x), y: r4(y), z: r4(z), b: q }); }
  }
  return out;
}
const SHAPES = ['sphere', 'cube', 'rod', 'hexrod', 'octahedron', 'wulff-like'];
/**
 * Carve a free nanoparticle from the bulk lattice (unrelaxed, unreconstructed surfaces; centred on a lattice point).
 * diameter_nm: sphere diameter / cube edge / rod diameter / hexagonal prism width across flats ({10-10} facets) /
 * octahedron vertex-to-vertex distance / {100} facet-to-facet distance of the truncated octahedron ('wulff-like',
 * {111}:{100} distance ratio √3/2 from the fcc broken-bond estimate; hexagonal crystals give a prism of height = width).
 * Rods lie along c (z); length_nm defaults to 3 × diameter. If the atom count would exceed maxAtoms the particle is scaled
 * down isotropically and truncated:true is returned.
 */
export function nanoparticle(keyOrDef, o = {}) {
  const s = asPeriodic(keyOrDef, 'nanoparticle'), shape = o.shape || 'sphere';
  if (!SHAPES.includes(shape)) throw new Error(`nanoparticle: unknown shape "${shape}". Use one of ${SHAPES.join(', ')}`);
  let d = (o.diameter_nm > 0 ? +o.diameter_nm : 3) * 10, L = (o.length_nm > 0 ? +o.length_nm : 3 * d / 10) * 10;
  const hexSys = /hex|trig/i.test(s.system || '') || Math.abs(cellParams(s.cell).gamma - 120) < 0.5;
  const maxAtoms = clamp(Math.round(o.maxAtoms ?? MAX_ATOMS), 50, 400000), rho = s.atoms.length / Math.abs(det3(s.cell)), S3 = Math.sqrt(3);
  if (shape === 'wulff-like' && hexSys) L = d;
  const vol = () => shape === 'sphere' ? Math.PI * d ** 3 / 6 : shape === 'cube' ? d ** 3 : shape === 'rod' ? Math.PI * d * d * L / 4 : shape === 'hexrod' || (shape === 'wulff-like' && hexSys) ? S3 / 2 * d * d * L : shape === 'octahedron' ? d ** 3 / 6 : 0.75 * d ** 3;
  let truncated = false; const est = rho * vol();
  if (est > maxAtoms) { const k = Math.cbrt(maxAtoms / est) * 0.985; d *= k; L *= k; truncated = true; }
  const R = d / 2, hz = L / 2, e = 1e-6, c30 = S3 / 2; let inside, rb;
  const hexIn = (x, y, z) => Math.abs(z) <= hz + e && Math.abs(y) <= R + e && Math.abs(c30 * x + y / 2) <= R + e && Math.abs(c30 * x - y / 2) <= R + e;
  if (shape === 'sphere') { inside = (x, y, z) => x * x + y * y + z * z <= R * R + e; rb = R; }
  else if (shape === 'cube') { inside = (x, y, z) => Math.abs(x) <= R + e && Math.abs(y) <= R + e && Math.abs(z) <= R + e; rb = R * S3; }
  else if (shape === 'rod') { inside = (x, y, z) => x * x + y * y <= R * R + e && Math.abs(z) <= hz + e; rb = Math.hypot(R, hz); }
  else if (shape === 'hexrod' || (shape === 'wulff-like' && hexSys)) { inside = hexIn; rb = Math.hypot(d / S3, hz); }
  else if (shape === 'octahedron') { inside = (x, y, z) => Math.abs(x) + Math.abs(y) + Math.abs(z) <= R + e; rb = R; }
  else { const h111 = S3 * (S3 / 2) * R; inside = (x, y, z) => Math.abs(x) <= R + e && Math.abs(y) <= R + e && Math.abs(z) <= R + e && Math.abs(x) + Math.abs(y) + Math.abs(z) <= h111 + e; rb = R * S3; }
  const atoms = carve(s, inside, rb + 0.01);
  let bulkCN; try { bulkCN = Array.from(bondGraph(s).cn); } catch { bulkCN = undefined; }
  const rods = shape === 'rod' || shape === 'hexrod';
  return { kind: 'atoms', atoms, cell: null, title: `${titleOf(keyOrDef)} ${shape} nanoparticle, ${(d / 10).toFixed(2)} nm${rods ? ` × ${(L / 10).toFixed(2)} nm` : ''}`, shape, diameter_nm: d / 10, length_nm: rods || (shape === 'wulff-like' && hexSys) ? L / 10 : null, nAtoms: atoms.length, truncated, bulkCN, source: titleOf(keyOrDef) };
}

/** Spherical core–shell particle: core_nm = core diameter, shell_nm = shell thickness. Lattices are not strained to match; shell atoms closer than 0.8 × the covalent-radius sum to a core atom are removed. */
export function coreShell(coreKey, shellKey, o = {}) {
  const sc = asPeriodic(coreKey, 'coreShell'), ss = asPeriodic(shellKey, 'coreShell');
  let rc = (o.core_nm > 0 ? +o.core_nm : 2) * 5, ts = (o.shell_nm > 0 ? +o.shell_nm : 0.5) * 10;
  const maxAtoms = clamp(Math.round(o.maxAtoms ?? MAX_ATOMS), 50, 400000), rhoC = sc.atoms.length / Math.abs(det3(sc.cell)), rhoS = ss.atoms.length / Math.abs(det3(ss.cell));
  const est = 4 / 3 * Math.PI * (rhoC * rc ** 3 + rhoS * ((rc + ts) ** 3 - rc ** 3)); let truncated = false;
  if (est > maxAtoms) { const k = Math.cbrt(maxAtoms / est) * 0.985; rc *= k; ts *= k; truncated = true; }
  const ro = rc + ts, core = carve(sc, (x, y, z) => x * x + y * y + z * z <= rc * rc, rc).map(a => ({ el: a.el, x: a.x, y: a.y, z: a.z, region: 'core' }));
  const shell = carve(ss, (x, y, z) => { const r2 = x * x + y * y + z * z; return r2 > rc * rc && r2 <= ro * ro; }, ro).map(a => ({ el: a.el, x: a.x, y: a.y, z: a.z, region: 'shell' }));
  const all = core.concat(shell), nC = core.length, drop = new Uint8Array(all.length);
  neighbourScan({ atoms: all, cell: null }, 2.6, (i, j, d) => { if (i >= nC && j < nC && d < 0.8 * (radOf(all[i].el) + radOf(all[j].el))) drop[i] = 1; });
  const atoms = all.filter((a, i) => !drop[i]), pc = cellParams(sc.cell), ps = cellParams(ss.cell);
  const sameType = sc.atoms.length === ss.atoms.length && Math.abs(pc.gamma - ps.gamma) < 0.5;
  return { kind: 'atoms', atoms, cell: null, title: `${titleOf(coreKey)}@${titleOf(shellKey)} core–shell, core ${(rc / 5).toFixed(2)} nm, shell ${(ts / 10).toFixed(2)} nm`, core_nm: rc / 5, shell_nm: ts / 10, diameter_nm: ro / 5, nCore: nC, nShell: atoms.length - nC, removedAtInterface: all.length - atoms.length, latticeMismatch_pct: sameType ? (ps.a - pc.a) / pc.a * 100 : null, truncated };
}

/** Random point defects (seeded): vacancies first, then substitution of a fraction of the remaining `from` atoms by `to`. */
export function addDefects(struct, o = {}) {
  if (!struct?.atoms?.length) throw new Error('addDefects: structure has no atoms');
  const rand = rng(o.seed ?? 1), n = struct.atoms.length, idx = Array.from({ length: n }, (_, i) => i);
  const pick = (arr, k) => { for (let i = 0; i < k; i++) { const j = i + rand.int(arr.length - i); [arr[i], arr[j]] = [arr[j], arr[i]]; } return arr.slice(0, k); };
  const nV = Math.min(n - 1, Math.round(clamp(+o.vacancyFraction || 0, 0, 1) * n)), gone = new Uint8Array(n); for (const i of pick(idx, nV)) gone[i] = 1;
  const atoms = struct.atoms.filter((a, i) => !gone[i]).map(a => ({ ...a })); let nS = 0;
  const sub = o.substitution;
  if (sub && sub.from && sub.to && sub.fraction > 0) {
    if (!ELEMENTS[sub.to]) throw new Error(`addDefects: unknown substituent element "${sub.to}"`);
    const cand = []; atoms.forEach((a, i) => { if (a.el === sub.from) cand.push(i); });
    nS = Math.round(clamp(+sub.fraction, 0, 1) * cand.length); for (const i of pick(cand, nS)) { atoms[i].el = sub.to; atoms[i].substituted = true; }
  }
  const out = { ...struct, atoms, title: `${struct.title || 'structure'} with defects`, defects: { vacancies: nV, substitutions: nS, vacancyFraction: nV / n, seed: o.seed ?? 1 } };
  delete out.formula; return out;
}

/** Single-wall nanotube by rolling a honeycomb sheet along Ch = n·a1 + m·a2. element 'C' (default) or 'BN'. */
export function nanotube(o = {}) {
  const n = Math.round(o.n ?? 10), m = Math.round(o.m ?? 0), bond = o.bond_A > 0 ? +o.bond_A : 1.42, element = o.element || 'C';
  if (!(n >= 0 && m >= 0 && n + m > 0)) throw new Error('nanotube: chiral indices need n ≥ 0, m ≥ 0 and n + m > 0');
  const els = element === 'BN' ? ['B', 'N'] : [element, element];
  if (!ELEMENTS[els[0]] || !ELEMENTS[els[1]]) throw new Error(`nanotube: unknown element "${element}"`);
  const a = Math.sqrt(3) * bond, a1 = [a, 0], a2 = [a / 2, a * Math.sqrt(3) / 2], N2 = n * n + n * m + m * m, C = a * Math.sqrt(N2), dR = gcd(2 * n + m, 2 * m + n);
  const T = Math.sqrt(3) * C / dR, perCell = 4 * N2 / dR, diameter = C / Math.PI, Rr = diameter / 2;
  let L = (o.length_nm > 0 ? +o.length_nm : 2) * 10, periodic = o.periodic !== false && T <= 1.5 * L, truncated = false;
  const maxAtoms = clamp(Math.round(o.maxAtoms ?? MAX_ATOMS), 50, 400000), perA = perCell / T;
  if (perA * L > maxAtoms) { L = maxAtoms / perA; truncated = true; if (T > L) periodic = false; }
  if (periodic) L = Math.max(1, Math.round(L / T - (truncated ? 0.5 : 0))) * T;
  const ch = [(n * a1[0] + m * a2[0]) / C, (n * a1[1] + m * a2[1]) / C], th = [-ch[1], ch[0]];      // unit vectors along Ch and the tube axis
  const corners = [[0, 0], [ch[0] * C, ch[1] * C], [th[0] * L, th[1] * L], [ch[0] * C + th[0] * L, ch[1] * C + th[1] * L]];
  const toIJ = p => { const j = p[1] / a2[1]; return [(p[0] - j * a2[0]) / a, j]; }, ij = corners.map(toIJ);
  const i0 = Math.floor(Math.min(...ij.map(v => v[0]))) - 2, i1 = Math.ceil(Math.max(...ij.map(v => v[0]))) + 2, j0 = Math.floor(Math.min(...ij.map(v => v[1]))) - 2, j1 = Math.ceil(Math.max(...ij.map(v => v[1]))) + 2;
  const basis = [[0, 0], [(a1[0] + a2[0]) / 3, (a1[1] + a2[1]) / 3]], atoms = [], eps = 1e-7, off = periodic ? (diameter + 20) / 2 : 0;
  for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) for (let b = 0; b < 2; b++) {
    const x = i * a1[0] + j * a2[0] + basis[b][0], y = i * a1[1] + j * a2[1] + basis[b][1], s = (x * ch[0] + y * ch[1]) / C, t = x * th[0] + y * th[1];
    if (s < -eps || s >= 1 - eps || t < -eps || t >= L - eps) continue;
    const phi = 2 * Math.PI * s; atoms.push({ el: els[b], x: r4(Rr * Math.cos(phi) + off), y: r4(Rr * Math.sin(phi) + off), z: r4(t) });
  }
  return { kind: 'atoms', atoms, cell: periodic ? [[diameter + 20, 0, 0], [0, diameter + 20, 0], [0, 0, L]] : null, title: `(${n},${m}) ${element === 'C' ? 'carbon' : element} nanotube, ${(L / 10).toFixed(2)} nm`, n, m, diameter_nm: diameter / 10, chiralAngle_deg: Math.atan2(Math.sqrt(3) * m, 2 * n + m) / D2R, metallic: element === 'C' && (n - m) % 3 === 0, length_nm: L / 10, translation_A: T, atomsPerUnitCell: perCell, periodic, truncated };
}

/**
 * Surface slab for Miller plane (hkl) of the conventional cell ([h,k,i,l] accepted for hexagonal). The in-plane cell is the
 * primitive 2-D cell of that lattice plane; `layers` counts atomic planes (planes closer than 0.25 Å are one layer), or
 * van der Waals layers when the stacking contains gaps ≥ 2.4 Å next to smaller ones (graphite, MoS2, hBN). Unrelaxed
 * bulk termination; polar or non-stoichiometric terminations are not corrected. Output: a ∥ x, surface normal ∥ z.
 */
export function slab(keyOrDef, o = {}) {
  const s = asPeriodic(keyOrDef, 'slab'); let hkl = (o.hkl || [0, 0, 1]).map(Number);
  if (hkl.length === 4) hkl = [hkl[0], hkl[1], hkl[3]];
  if (hkl.length !== 3 || hkl.some(v => !Number.isInteger(v)) || !hkl.some(v => v) || hkl.some(v => Math.abs(v) > 8)) throw new Error('slab: hkl must be three (or four, hexagonal) integers with |index| ≤ 8, not all zero');
  const g = hkl.reduce(gcd); hkl = hkl.map(v => v / g);
  const [h, k, l] = hkl, M = s.cell, R = Math.max(1, ...hkl.map(Math.abs)) + 2, vacuum = Math.max(0, o.vacuum_A ?? 15), inP = [], out = [];
  for (let u = -R; u <= R; u++) for (let v = -R; v <= R; v++) for (let w = -R; w <= R; w++) {
    const p = h * u + k * v + l * w; if (p !== 0 && p !== 1) continue; if (!u && !v && !w) continue;
    const c = rowMul([u, v, w], M); (p ? out : inP).push({ n: [u, v, w], c, L: len(c) });
  }
  const nice = (A, B) => (Math.abs(A.L - B.L) > 1e-6 ? A.L - B.L : 0) || A.n.filter(v => v).length - B.n.filter(v => v).length || (B.n[0] + B.n[1] + B.n[2]) - (A.n[0] + A.n[1] + A.n[2]) || B.n[0] - A.n[0] || B.n[1] - A.n[1] || B.n[2] - A.n[2];
  inP.sort(nice); out.sort(nice);
  const U = inP[0]; let V = null, best = Infinity;
  for (const c of inP) { const ar = len(cross(U.c, c.c)); if (ar > 1e-6 && ar < best - 1e-6) { best = ar; V = c; } }
  if (!V || !out.length) throw new Error('slab: could not construct an in-plane cell for this plane');
  const W = out[0]; let u = U.c, v = V.c; const w = W.c;
  if (dot(cross(u, v), w) < 0) v = v.map(c => -c);
  const nrm = cross(u, v), area = len(nrm), nh = nrm.map(c => c / area), dz = dot(w, nh), Si = inv3([u, v, w]);
  let pts = s.atoms.map(a => { const f = rowMul([a.x, a.y, a.z], Si).map(wrap1); return { el: a.el, r: rowMul(f, [u, v, w]), z: f[2] * dz }; }).sort((A, B) => A.z - B.z);
  let gi = 0, gmax = -1; for (let i = 0; i < pts.length; i++) { const gap = i + 1 < pts.length ? pts[i + 1].z - pts[i].z : pts[0].z + dz - pts[i].z; if (gap > gmax + 1e-9) { gmax = gap; gi = i; } }
  pts = pts.slice(gi + 1).concat(pts.slice(0, gi + 1).map(p => ({ el: p.el, r: [p.r[0] + w[0], p.r[1] + w[1], p.r[2] + w[2]], z: p.z + dz })));
  let planes = []; for (const p of pts) { const last = planes[planes.length - 1]; if (last && p.z - last[last.length - 1].z < 0.25) last.push(p); else planes.push([p]); }
  const top = pl => pl[pl.length - 1].z, gaps = planes.map((pl, i) => i + 1 < planes.length ? planes[i + 1][0].z - top(pl) : planes[0][0].z + dz - top(pl));
  let layerMode = 'atomic planes';
  if (Math.max(...gaps) >= 2.4 && Math.min(...gaps) < 2.4) { const merged = [[...planes[0]]]; for (let i = 1; i < planes.length; i++) { if (gaps[i - 1] >= 2.4) merged.push([...planes[i]]); else merged[merged.length - 1].push(...planes[i]); } planes = merged; layerMode = 'van der Waals layers'; }
  const nL = planes.length, per = s.atoms.length / nL; let layers = Math.max(1, Math.round(o.layers ?? 6)), truncated = false;
  if (layers * per > MAX_ATOMS) { layers = Math.max(1, Math.floor(MAX_ATOMS / per)); truncated = true; }
  const ex = u.map(c => c / len(u)), ey = cross(nh, ex), lu = len(u), v2 = [dot(v, ex), dot(v, ey)], raw = []; let zmin = Infinity, zmax = -Infinity;
  for (let q = 0; q < layers; q++) { const p = Math.floor(q / nL); for (const a of planes[q % nL]) { const r = [a.r[0] + p * w[0], a.r[1] + p * w[1], a.r[2] + p * w[2]], Z = dot(r, nh), Y = dot(r, ey), X = dot(r, ex), fv = wrap1(Y / v2[1]), fu = wrap1((X - (Y / v2[1]) * v2[0]) / lu); raw.push({ el: a.el, fu, fv, Z }); if (Z < zmin) zmin = Z; if (Z > zmax) zmax = Z; } }
  const thick = zmax - zmin, cz = thick + vacuum || 1;
  const atoms = raw.map(a => ({ el: a.el, x: r6(a.fu * lu + a.fv * v2[0]), y: r6(a.fv * v2[1]), z: r6(a.Z - zmin + vacuum / 2), f: [a.fu, a.fv, (a.Z - zmin + vacuum / 2) / cz] }));
  return { kind: 'atoms', atoms, cell: [[lu, 0, 0], [r6(v2[0]), v2[1], 0], [0, 0, cz]], title: `${titleOf(keyOrDef)} (${hkl.join(' ')}) slab, ${layers} layers`, hkl, layers, layersPerPeriod: nL, layerMode, d_hkl_A: dz, thickness_A: thick, vacuum_A: vacuum, area_A2: area, truncated };
}

/** Free-standing nanosheet: (001) slab of `layers` layers repeated nx × ny in-plane, with vacuum along z. */
export function nanosheet(keyOrDef, o = {}) {
  const base = slab(keyOrDef, { hkl: [0, 0, 1], layers: o.layers ?? 1, vacuum_A: o.vacuum_A ?? 15 });
  let nx = Math.max(1, Math.round(o.nx ?? 4)), ny = Math.max(1, Math.round(o.ny ?? 4)), truncated = base.truncated;
  while (nx * ny * base.atoms.length > MAX_ATOMS && (nx > 1 || ny > 1)) { if (nx >= ny) nx--; else ny--; truncated = true; }
  const sc = supercell(base, [nx, ny, 1]), p = cellParams(sc.cell), warnings = [];
  if (base.layers > 1 && base.layerMode === 'atomic planes' && base.d_hkl_A > 8) warnings.push('The source cell already contains vacuum along c; stacked layers are separated by that vacuum. Use the bulk layered crystal (e.g. graphite) for multilayers.');
  return { ...sc, title: `${titleOf(keyOrDef)} nanosheet ${nx}×${ny}, ${base.layers} layer${base.layers > 1 ? 's' : ''}`, nx, ny, size_nm: [p.a / 10, p.b / 10], truncated, warnings };
}

// ───────────────────────────── diffraction ─────────────────────────────
// Cromer–Mann coefficients f(s) = Σ a_i exp(−b_i s²) + c, s = sinθ/λ (Å⁻¹), neutral atoms: a1 b1 a2 b2 a3 b3 a4 b4 c.
const CM = {};
('H 0.489918 20.6593 0.262003 7.74039 0.196767 49.5519 0.049879 2.20159 0.001305|Li 1.1282 3.9546 0.7508 1.0524 0.6175 85.3905 0.4653 168.261 0.0377|B 2.0545 23.2185 1.3326 1.021 1.0979 60.3498 0.7068 0.1403 -0.1932|C 2.31 20.8439 1.02 10.2075 1.5886 0.5687 0.865 51.6512 0.2156|N 12.2126 0.0057 3.1322 9.8933 2.0125 28.9975 1.1663 0.5826 -11.529|O 3.0485 13.2771 2.2868 5.7011 1.5463 0.3239 0.867 32.9089 0.2508|F 3.5392 10.2825 2.6412 4.2944 1.517 0.2615 1.0243 26.1476 0.2776|Na 4.7626 3.285 3.1736 8.8422 1.2674 0.3136 1.1128 129.424 0.676|Mg 5.4204 2.8275 2.1735 79.2611 1.2269 0.3808 2.3073 7.1937 0.8584|Al 6.4202 3.0387 1.9002 0.7426 1.5936 31.5472 1.9646 85.0886 1.1151|Si 6.2915 2.4386 3.0353 32.3337 1.9891 0.6785 1.541 81.6937 1.1407|P 6.4345 1.9067 4.1791 27.157 1.78 0.526 1.4908 68.1645 1.1149|S 6.9053 1.4679 5.2034 22.2151 1.4379 0.2536 1.5863 56.172 0.8669|Cl 11.4604 0.0104 7.1964 1.1662 6.2556 18.5194 1.6455 47.7784 -9.5574|K 8.2186 12.7949 7.4398 0.7748 1.0519 213.187 0.8659 41.6841 1.4228|Ca 8.6266 10.4421 7.3873 0.6599 1.5899 85.7484 1.0211 178.437 1.3751|Ti 9.7595 7.8508 7.3558 0.5 1.6991 35.6338 1.9021 116.105 1.2807|V 10.2971 6.8657 7.3511 0.4385 2.0703 26.8938 2.0571 102.478 1.2199|Cr 10.6406 6.1038 7.3537 0.392 3.324 20.2626 1.4922 98.7399 1.1832|Mn 11.2819 5.3409 7.3573 0.3432 3.0193 17.8674 2.2441 83.7543 1.0896|Fe 11.7695 4.7611 7.3573 0.3072 3.5222 15.3535 2.3045 76.8805 1.0369|Co 12.2841 4.2791 7.3409 0.2784 4.0034 13.5359 2.3488 71.1692 1.0118|Ni 12.8376 3.8785 7.292 0.2565 4.4438 12.1763 2.38 66.3421 1.0341|Cu 13.338 3.5828 7.1676 0.247 5.6158 11.3966 1.6735 64.8126 1.191|Zn 14.0743 3.2655 7.0318 0.2333 5.1652 10.3163 2.41 58.7097 1.3041|Ga 15.2354 3.0669 6.7006 0.2412 4.3591 10.7805 2.9623 61.4135 1.7189|Ge 16.0816 2.8509 6.3747 0.2516 3.7068 11.4468 3.683 54.7625 2.1313|As 16.6723 2.6345 6.0701 0.2647 3.4313 12.9479 4.2779 47.7972 2.531|Se 17.0006 2.4098 5.8196 0.2726 3.9731 15.2372 4.3543 43.8163 2.8409|Sr 17.5663 1.5564 9.8184 14.0988 5.422 0.1664 2.6694 132.376 2.5064|Zr 17.8765 1.27618 10.948 11.916 5.41732 0.117622 3.65721 87.6627 2.06929|Mo 3.7025 0.2772 17.2356 1.0958 12.8876 11.004 3.7429 61.6584 4.3875|Pd 19.3319 0.698655 15.5017 7.98929 5.29537 25.2052 0.605844 76.8986 5.26593|Ag 19.2808 0.6446 16.6885 7.4726 4.8045 24.6605 1.0463 99.8156 5.179|Cd 19.2214 0.5946 17.6444 6.9089 4.461 24.7008 1.6029 87.4825 5.0694|In 19.1624 0.5476 18.5596 6.3776 4.2948 25.8499 2.0396 92.8029 4.9391|Sn 19.1889 5.8303 19.1005 0.5031 4.4585 26.8909 2.4663 83.9571 4.7821|Cs 20.3892 3.569 19.1062 0.3107 10.662 24.3879 1.4953 213.904 3.3352|Ba 20.3361 3.216 19.297 0.2756 10.888 20.2073 2.6959 167.202 2.7731|Ce 21.1671 2.81219 19.7695 0.226836 11.8513 17.6083 3.33049 127.113 1.86264|W 29.0818 1.72029 15.43 9.2259 14.4327 0.321703 5.11982 57.056 9.8875|Pt 27.0059 1.51293 17.7639 8.81174 15.7131 0.424593 5.7837 38.6103 11.6883|Au 16.8819 0.4611 18.5913 8.6216 25.5582 1.4826 5.86 36.3956 12.0658|Pb 31.0617 0.6902 13.0637 2.3576 18.442 8.618 5.9696 47.2579 13.4118')
  .split('|').forEach(r => { const t = r.split(' '); CM[t[0]] = t.slice(1).map(Number); });
/**
 * X-ray atomic scattering factor of a neutral atom at s = sinθ/λ (Å⁻¹). Tabulated elements use the four-Gaussian
 * Cromer–Mann form; all others use the Thomas–Fermi atom with Molière screening, f = Z Σ α_i μ_i²/(q² + μ_i²),
 * q = 4πs, μ_i = β_i Z^(1/3)/0.4685 Å (accurate to roughly 5–15 %). No anomalous dispersion.
 */
export function scatteringFactor(el, s) {
  const c = CM[el], s2 = s * s;
  if (c) return c[0] * Math.exp(-c[1] * s2) + c[2] * Math.exp(-c[3] * s2) + c[4] * Math.exp(-c[5] * s2) + c[6] * Math.exp(-c[7] * s2) + c[8];
  const Z = ELEMENTS[el]?.Z || 0; if (!Z) return 0;
  const k = Math.cbrt(Z) / 0.4685, q2 = (4 * Math.PI * s) ** 2, t = (al, be) => { const m2 = (be * k) ** 2; return al * m2 / (q2 + m2); };
  return Z * (t(0.1, 6) + t(0.55, 1.2) + t(0.35, 0.3));
}

const hklNice = (A, B) => (B.filter(v => v >= 0).length - A.filter(v => v >= 0).length) || B[0] - A[0] || B[1] - A[1] || B[2] - A[2];
// All reflections with d ≥ dmin, merged by d-spacing. F² summed over the merged set (so multiplicity is included).
function reflectionGroups(s, dmin, B = 0, maxIndex = Infinity) {
  const M = s.cell, Mi = inv3(M), V = Math.abs(det3(M)), at = s.atoms, n = at.length; let truncated = false;
  const work = d => 2 / 3 * Math.PI * V / d ** 3 * n;
  if (work(dmin) > 4e7) { dmin = Math.cbrt(2 / 3 * Math.PI * V * n / 4e7); truncated = true; }
  const as = [0, 1, 2].map(k => [Mi[0][k], Mi[1][k], Mi[2][k]]), hm = [0, 1, 2].map(k => Math.min(maxIndex, Math.floor(len(M[k]) / dmin + 1e-9)));
  const fr = at.map(a => rowMul([a.x, a.y, a.z], Mi)), els = [...new Set(at.map(a => a.el))], ei = at.map(a => els.indexOf(a.el)), occ = at.map(a => a.occ ?? 1);
  const g2max = 1 / (dmin * dmin) * (1 + 1e-9), refl = [], fE = new Float64Array(els.length), TWO_PI = 2 * Math.PI;
  for (let h = 0; h <= hm[0]; h++) for (let k = -hm[1]; k <= hm[1]; k++) for (let l = -hm[2]; l <= hm[2]; l++) {
    if (h === 0 && (k < 0 || (k === 0 && l <= 0))) continue;                                         // half-space (Friedel pairs counted ×2)
    const gx = h * as[0][0] + k * as[1][0] + l * as[2][0], gy = h * as[0][1] + k * as[1][1] + l * as[2][1], gz = h * as[0][2] + k * as[1][2] + l * as[2][2], g2 = gx * gx + gy * gy + gz * gz;
    if (g2 > g2max) continue;
    const sv = Math.sqrt(g2) / 2, dw = Math.exp(-B * sv * sv); let f0 = 0;
    for (let e = 0; e < els.length; e++) { fE[e] = scatteringFactor(els[e], sv) * dw; }
    let A = 0, Bs = 0; for (let q = 0; q < n; q++) { const ph = TWO_PI * (h * fr[q][0] + k * fr[q][1] + l * fr[q][2]), w = fE[ei[q]] * occ[q]; A += w * Math.cos(ph); Bs += w * Math.sin(ph); f0 += w; }
    const F2 = A * A + Bs * Bs; refl.push({ hkl: [h, k, l], g2, F2, on: F2 > 1e-6 * f0 * f0 });
  }
  refl.sort((a, b) => a.g2 - b.g2);
  const groups = [];
  for (const r of refl) {
    let g = groups[groups.length - 1];
    if (!g || r.g2 - g.g2 > 2e-6 * r.g2) { g = { g2: r.g2, d: 1 / Math.sqrt(r.g2), F2: 0, mult: 0, multLattice: 0, on: [], all: [] }; groups.push(g); }
    g.multLattice += 2; g.all.push(r.hkl); if (r.on) { g.F2 += 2 * r.F2; g.mult += 2; g.on.push(r.hkl); }
  }
  for (const g of groups) { const c = (g.on.length ? g.on : g.all).flatMap(v => [v, v.map(x => (x === 0 ? 0 : -x))]).sort(hklNice); g.hkl = c[0]; g.allowed = g.on.length > 0; delete g.on; delete g.all; }
  return { groups, truncated, dmin, V };
}

/** Lattice-plane spacings d ≥ dmin_A. Systematically absent reflections (|F|² = 0 for the given basis) are omitted unless includeForbidden is set. */
export function dSpacings(def, { dmin_A = 0.8, maxIndex, includeForbidden = false } = {}) {
  const s = asPeriodic(def, 'dSpacings'), { groups } = reflectionGroups(s, Math.max(0.3, dmin_A > 0 ? dmin_A : 0.8), 0, maxIndex > 0 ? Math.round(maxIndex) : Infinity);
  return groups.filter(g => includeForbidden || g.allowed).sort((a, b) => b.d - a.d).map(g => ({ hkl: g.hkl, d_A: g.d, multiplicity: g.allowed ? g.mult : g.multLattice, allowed: g.allowed }));
}

/**
 * Kinematic powder X-ray diffraction pattern of a periodic structure. I(hkl) = m |F|² LP with F = Σ f_j(s) exp(−B s²)
 * exp(2πi h·r_j), LP = (1 + cos²2θ)/(sin²θ cosθ). Peak shape: Thompson–Cox–Hastings pseudo-Voigt with Lorentzian width
 * from the Scherrer relation (K = 0.9) and Gaussian width from microstrain (4 ε tanθ) plus an instrumental term.
 * Single wavelength (no Kα2), no absorption, preferred orientation or anomalous dispersion.
 */
export function powderXRD(x, o = {}) {
  const s = asPeriodic(x, 'powderXRD'), lam = o.wavelength_A > 0 ? +o.wavelength_A : 1.5406;
  const tmin = clamp(o.twoThetaMin ?? 10, 0.5, 177), tmax = clamp(o.twoThetaMax ?? 90, tmin + 1, 179);
  let step = o.step > 0 ? +o.step : 0.02, npts = Math.floor((tmax - tmin) / step + 1e-9) + 1; if (npts > 40001) { npts = 40001; step = (tmax - tmin) / 40000; }
  const Lc = o.crystallite_nm > 0 ? o.crystallite_nm * 10 : Infinity, strain = Math.max(0, +o.strain || 0), inst = Math.max(0, o.instrument_deg ?? 0.08), B = Math.max(0, o.Biso_A2 ?? 0.5);
  const tExt = Math.min(179.5, tmax + 3), { groups, truncated } = reflectionGroups(s, lam / (2 * Math.sin(tExt / 2 * D2R)), B);
  const twoTheta = new Array(npts), y = new Float64Array(npts); for (let i = 0; i < npts; i++) twoTheta[i] = r6(tmin + i * step);
  const pk = [], LN2 = Math.LN2;
  for (const g of groups) {
    if (!g.allowed) continue; const sn = lam / (2 * g.d); if (sn >= 1) continue;
    const th = Math.asin(sn), tt = 2 * th / D2R; if (tt < tmin - 3 || tt > tExt) continue;
    const ct = Math.cos(th), lp = (1 + Math.cos(2 * th) ** 2) / (sn * sn * ct), Iint = g.F2 * lp;
    let gL = Number.isFinite(Lc) ? 0.9 * lam / (Lc * ct) / D2R : 0, gG = Math.hypot(inst, 4 * strain * Math.tan(th) / D2R);
    if (gL < 2 * step && gG < 2 * step) gG = 2 * step;                                              // keep peaks resolvable on the sampling grid
    const G = (gG ** 5 + 2.69269 * gG ** 4 * gL + 2.42843 * gG ** 3 * gL * gL + 4.47163 * gG * gG * gL ** 3 + 0.07842 * gG * gL ** 4 + gL ** 5) ** 0.2, q = gL / G, eta = clamp(1.36603 * q - 0.47719 * q * q + 0.11116 * q ** 3, 0, 1);
    pk.push({ hkl: g.hkl, tt, d: g.d, Iint, G, eta, mult: g.mult });
  }
  for (const p of pk) {
    const win = Math.max(40 * p.G * (p.eta > 0.02 ? 1 : 0.1), 5 * p.G), i0 = Math.max(0, Math.ceil((p.tt - win - tmin) / step)), i1 = Math.min(npts - 1, Math.floor((p.tt + win - tmin) / step));
    const cL = p.eta * 2 / (Math.PI * p.G), cG = (1 - p.eta) * 2 / p.G * Math.sqrt(LN2 / Math.PI), k = 4 / (p.G * p.G);
    for (let i = i0; i <= i1; i++) { const dx = tmin + i * step - p.tt, u = k * dx * dx; y[i] += p.Iint * (cL / (1 + u) + cG * Math.exp(-LN2 * u)); }
  }
  let ymax = 0; for (let i = 0; i < npts; i++) if (y[i] > ymax) ymax = y[i];
  const intensity = new Array(npts); for (let i = 0; i < npts; i++) intensity[i] = ymax > 0 ? Math.round(y[i] / ymax * 1e5) / 1e3 : 0;
  const inR = pk.filter(p => p.tt >= tmin && p.tt <= tmax), imax = inR.reduce((m, p) => Math.max(m, p.Iint), 0);
  const peaks = inR.map(p => ({ hkl: p.hkl, twoTheta: p.tt, d_A: p.d, I_rel: imax > 0 ? p.Iint / imax * 100 : 0, multiplicity: p.mult, fwhm_deg: p.G })).filter(p => p.I_rel >= 0.01).sort((a, b) => a.twoTheta - b.twoTheta);
  return { twoTheta, intensity, peaks, wavelength_A: lam, step, truncated, model: 'Kinematic powder diffraction: Cromer–Mann (or Thomas–Fermi/Molière) scattering factors, isotropic Debye–Waller factor, Lorentz-polarisation, Scherrer + microstrain pseudo-Voigt broadening' };
}

// ───────────────────────────── writers ─────────────────────────────
const fx = (v, n) => { const t = (+v).toFixed(n); return /^-0\.?0*$/.test(t) ? t.slice(1) : t; };
const pad = (v, n) => String(v).padStart(n);
function needAtoms(struct, who) { if (!struct || !Array.isArray(struct.atoms) || !struct.atoms.length) throw new Error(`${who}: structure has no atoms`); return struct; }
// Periodic structures pass through; finite ones are placed in an orthogonal box with padding.
function boxed(struct, padA = 10) {
  if (struct.cell) return { cell: struct.cell, atoms: struct.atoms, isBox: false };
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const a of struct.atoms) { const r = [a.x, a.y, a.z]; for (let k = 0; k < 3; k++) { if (r[k] < lo[k]) lo[k] = r[k]; if (r[k] > hi[k]) hi[k] = r[k]; } }
  const L = [0, 1, 2].map(k => hi[k] - lo[k] + 2 * padA);
  return { cell: [[L[0], 0, 0], [0, L[1], 0], [0, 0, L[2]]], atoms: struct.atoms.map(a => ({ ...a, x: a.x - lo[0] + padA, y: a.y - lo[1] + padA, z: a.z - lo[2] + padA })), isBox: true };
}
// Re-express a periodic structure with a ∥ x and b in the xy-plane (needed by formats that store only cell parameters or a restricted triclinic box).
function standardise(cell, atoms) {
  const p = cellParams(cell), M2 = latticeVectors(p), Mi = inv3(cell);
  return { p, cell: M2, atoms: atoms.map(a => { const f = rowMul([a.x, a.y, a.z], Mi), r = rowMul(f, M2); return { ...a, x: r[0], y: r[1], z: r[2], f }; }) };
}
const oneLine = t => String(t ?? '').replace(/[\r\n]+/g, ' ').trim();

export function toCIF(struct, o = {}) {
  needAtoms(struct, 'toCIF'); const b = boxed(struct, o.padding_A ?? 10), st = standardise(b.cell, b.atoms), p = st.p, comp = composition(struct.atoms), cnt = {};
  const name = (oneLine(o.title || struct.title || 'structure').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '') || 'structure').slice(0, 60);
  const L = [`data_${name}`, `_chemical_formula_sum '${Object.keys(comp).sort().map(e => e + (comp[e] === 1 ? '' : comp[e])).join(' ')}'`,
    `_cell_length_a ${fx(p.a, 6)}`, `_cell_length_b ${fx(p.b, 6)}`, `_cell_length_c ${fx(p.c, 6)}`, `_cell_angle_alpha ${fx(p.alpha, 5)}`, `_cell_angle_beta ${fx(p.beta, 5)}`, `_cell_angle_gamma ${fx(p.gamma, 5)}`, `_cell_volume ${fx(p.volume, 5)}`,
    "_symmetry_space_group_name_H-M 'P 1'", '_symmetry_Int_Tables_number 1', 'loop_', '_symmetry_equiv_pos_as_xyz', "'x, y, z'", 'loop_', '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z', '_atom_site_occupancy'];
  if (b.isBox) L.splice(1, 0, '# finite structure written in a padded orthogonal box (P 1)');
  for (const a of st.atoms) { cnt[a.el] = (cnt[a.el] || 0) + 1; L.push(`${a.el}${cnt[a.el]} ${a.el} ${fx(a.f[0], 6)} ${fx(a.f[1], 6)} ${fx(a.f[2], 6)} ${fx(a.occ ?? 1, 3)}`); }
  return L.join('\n') + '\n';
}
export function toXYZ(struct, o = {}) {
  needAtoms(struct, 'toXYZ');
  return [struct.atoms.length, oneLine(o.title || struct.title || ''), ...struct.atoms.map(a => `${a.el.padEnd(2)} ${pad(fx(a.x, 6), 14)} ${pad(fx(a.y, 6), 14)} ${pad(fx(a.z, 6), 14)}`)].join('\n') + '\n';
}
export function toExtXYZ(struct, o = {}) {
  needAtoms(struct, 'toExtXYZ'); const c = struct.cell, t = oneLine(o.title || struct.title || '').replace(/"/g, "'");
  const head = (c ? `Lattice="${c.flat().map(v => fx(v, 8)).join(' ')}" ` : '') + `Properties=species:S:1:pos:R:3 pbc="${c ? 'T T T' : 'F F F'}"` + (t ? ` comment="${t}"` : '');
  return [struct.atoms.length, head, ...struct.atoms.map(a => `${a.el.padEnd(2)} ${pad(fx(a.x, 8), 16)} ${pad(fx(a.y, 8), 16)} ${pad(fx(a.z, 8), 16)}`)].join('\n') + '\n';
}
export function toPOSCAR(struct, o = {}) {
  needAtoms(struct, 'toPOSCAR'); const b = boxed(struct, o.padding_A ?? 10), Mi = inv3(b.cell), els = [...new Set(b.atoms.map(a => a.el))];
  const L = [oneLine(o.title || struct.title || 'structure') || 'structure', '1.0', ...b.cell.map(v => v.map(c => pad(fx(c, 10), 18)).join(' ')), els.join(' '), els.map(e => b.atoms.filter(a => a.el === e).length).join(' '), 'Direct'];
  for (const e of els) for (const a of b.atoms) if (a.el === e) { const f = rowMul([a.x, a.y, a.z], Mi); L.push(f.map(v => pad(fx(v, 10), 16)).join(' ')); }
  return L.join('\n') + '\n';
}
export function toXSF(struct, o = {}) {
  needAtoms(struct, 'toXSF'); const row = a => `${pad(ELEMENTS[a.el]?.Z ?? 0, 3)} ${pad(fx(a.x, 8), 16)} ${pad(fx(a.y, 8), 16)} ${pad(fx(a.z, 8), 16)}`, L = [`# ${oneLine(o.title || struct.title || 'structure')}`];
  if (struct.cell) L.push('CRYSTAL', 'PRIMVEC', ...struct.cell.map(v => v.map(c => pad(fx(c, 8), 16)).join(' ')), 'PRIMCOORD', `${struct.atoms.length} 1`); else L.push('ATOMS');
  return L.concat(struct.atoms.map(row)).join('\n') + '\n';
}
export function toPDB(struct, o = {}) {
  needAtoms(struct, 'toPDB'); const L = [`TITLE     ${oneLine(o.title || struct.title || 'structure').slice(0, 60)}`]; let atoms = struct.atoms;
  if (struct.cell) { const st = standardise(struct.cell, struct.atoms), p = st.p; atoms = st.atoms; L.push(`CRYST1${pad(fx(p.a, 3), 9)}${pad(fx(p.b, 3), 9)}${pad(fx(p.c, 3), 9)}${pad(fx(p.alpha, 2), 7)}${pad(fx(p.beta, 2), 7)}${pad(fx(p.gamma, 2), 7)} P 1           1`); }
  atoms.forEach((a, i) => {
    if (Math.abs(a.x) >= 1e4 || Math.abs(a.y) >= 1e4 || Math.abs(a.z) >= 1e4) throw new Error('toPDB: coordinates exceed the PDB fixed-column range (±9999.999 Å)');
    const nm = (a.el.length === 1 ? ' ' + a.el : a.el.toUpperCase()).padEnd(4);
    L.push(`HETATM${pad((i + 1) % 100000, 5)} ${nm} UNL A   1    ${pad(fx(a.x, 3), 8)}${pad(fx(a.y, 3), 8)}${pad(fx(a.z, 3), 8)}${pad(fx(a.occ ?? 1, 2), 6)}  0.00          ${pad(a.el, 2)}`);
  });
  L.push('END'); return L.join('\n') + '\n';
}
/** LAMMPS data file. style: 'atomic' (default) | 'charge' | 'full'. The cell is rotated to the LAMMPS convention (a ∥ x, b in xy). */
export function toLAMMPSData(struct, o = {}) {
  needAtoms(struct, 'toLAMMPSData'); const style = o.style || 'atomic';
  if (!['atomic', 'charge', 'full'].includes(style)) throw new Error(`toLAMMPSData: unsupported atom style "${style}" (use atomic, charge or full)`);
  const b = boxed(struct, o.padding_A ?? 10), st = standardise(b.cell, b.atoms), M = st.cell, els = [...new Set(st.atoms.map(a => a.el))], tri = Math.abs(M[1][0]) + Math.abs(M[2][0]) + Math.abs(M[2][1]) > 1e-9;
  const L = [`# ${oneLine(o.title || struct.title || 'structure')}`, '', `${st.atoms.length} atoms`, `${els.length} atom types`, '', `0.0 ${fx(M[0][0], 8)} xlo xhi`, `0.0 ${fx(M[1][1], 8)} ylo yhi`, `0.0 ${fx(M[2][2], 8)} zlo zhi`];
  if (tri) L.push(`${fx(M[1][0], 8)} ${fx(M[2][0], 8)} ${fx(M[2][1], 8)} xy xz yz`);
  L.push('', 'Masses', '', ...els.map((e, i) => `${i + 1} ${massOf(e) || 1} # ${e}`), '', `Atoms # ${style}`, '');
  st.atoms.forEach((a, i) => { const t = els.indexOf(a.el) + 1, xyz = `${fx(a.x, 8)} ${fx(a.y, 8)} ${fx(a.z, 8)}`, q = fx(a.q ?? 0, 6); L.push(style === 'atomic' ? `${i + 1} ${t} ${xyz}` : style === 'charge' ? `${i + 1} ${t} ${q} ${xyz}` : `${i + 1} ${a.mol ?? 1} ${t} ${q} ${xyz}`); });
  return L.join('\n') + '\n';
}
/** MDL molfile (V2000). Bonds come from struct.bonds ([[i,j,order]]) or are perceived geometrically (all written as single bonds). */
export function toMOL(struct, o = {}) {
  needAtoms(struct, 'toMOL'); const n = struct.atoms.length;
  if (n > 999) throw new Error(`toMOL: V2000 molfiles hold at most 999 atoms (structure has ${n}); use XYZ, PDB or CIF instead`);
  let bonds = Array.isArray(struct.bonds) ? struct.bonds : bondGraph({ atoms: struct.atoms, cell: null }, true).list.map(b => [b[0], b[1], 1]);
  if (bonds.length > 999) bonds = bonds.slice(0, 999);
  const L = [oneLine(o.title || struct.title || 'structure').slice(0, 80), '  DesignerMaterials          3D', '', `${pad(n, 3)}${pad(bonds.length, 3)}  0  0  0  0  0  0  0  0999 V2000`];
  for (const a of struct.atoms) L.push(`${pad(fx(a.x, 4), 10)}${pad(fx(a.y, 4), 10)}${pad(fx(a.z, 4), 10)} ${a.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  for (const b of bonds) L.push(`${pad(b[0] + 1, 3)}${pad(b[1] + 1, 3)}${pad(Math.round(b[2] || 1), 3)}  0`);
  L.push('M  END'); return L.join('\n') + '\n';
}

// ───────────────────────────── parser helpers ─────────────────────────────
const splitLines = t => { if (typeof t !== 'string') throw new Error('Parser input must be text (string)'); return t.replace(/^﻿/, '').split(/\r\n|\r|\n/); };
const NUM_RE = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
function num(tok, fmt, what, ln) {
  let t = String(tok ?? '').trim().replace(/[dD](?=[+-]?\d)/, 'e');
  if (!NUM_RE.test(t)) throw new Error(`${fmt}: expected a number for ${what}${ln ? ` on line ${ln}` : ''}, got "${String(tok ?? '').slice(0, 30)}"`);
  return parseFloat(t);
}
const isNum = t => NUM_RE.test(String(t ?? '').trim().replace(/[dD](?=[+-]?\d)/, 'e'));
const isInt = t => /^[+-]?\d+$/.test(String(t ?? '').trim());
const toks = l => l.trim().split(/\s+/).filter(Boolean);
/** Element symbol from a label such as "Zn1", "O2-", "FE", "30" (atomic number); 'X' when unrecognised. */
function normEl(raw) {
  const s = String(raw ?? '').trim(); if (/^\d+$/.test(s)) return SYMBOLS[+s - 1] || 'X';
  const m = s.match(/^[A-Za-z]{1,2}/); if (!m) return 'X';
  const two = m[0][0].toUpperCase() + (m[0][1] || '').toLowerCase(); if (ELEMENTS[two]) return two;
  const one = two[0]; return one === 'D' || one === 'T' ? 'H' : ELEMENTS[one] ? one : 'X';
}
const BIO = new Set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HOH SOL WAT TIP3 DA DC DG DT A C G U'.split(' '));
function elFromAtomName(name, res = '') {
  const al = String(name).replace(/[^A-Za-z]/g, ''); if (!al) return 'X';
  const one = al[0].toUpperCase(), two = one + (al[1] || '').toLowerCase();
  if (BIO.has(res.trim().toUpperCase())) return ELEMENTS[one] ? one : normEl(al);
  if (al.length >= 2 && ELEMENTS[two] && !(one === 'H' && !['He', 'Hg', 'Hf'].includes(two))) return two;
  return ELEMENTS[one] ? one : normEl(al);
}
function finish(fmt, atoms, cell, title, extra = {}) {
  if (!atoms.length) throw new Error(`${fmt}: no atoms found`);
  for (let i = 0; i < atoms.length; i++) if (!Number.isFinite(atoms[i].x) || !Number.isFinite(atoms[i].y) || !Number.isFinite(atoms[i].z)) throw new Error(`${fmt}: atom ${i + 1} has a non-finite coordinate`);
  if (cell) { if (!cell.flat().every(Number.isFinite)) throw new Error(`${fmt}: lattice vectors are not finite numbers`); if (!(Math.abs(det3(cell)) > 1e-9)) throw new Error(`${fmt}: lattice vectors are coplanar (zero cell volume)`); }
  const warnings = extra.warnings || []; const nx = atoms.filter(a => a.el === 'X').length;
  if (nx) warnings.push(`${nx} atom(s) have an unrecognised element and are labelled X`);
  return { kind: 'atoms', atoms, cell: cell || null, title: oneLine(title) || fmt, format: fmt, ...extra, warnings };
}
const frac = t => { const m = String(t).trim().match(/^([+-]?[\d.]+(?:[eE][+-]?\d+)?)\/([\d.]+)$/); return m ? parseFloat(m[1]) / parseFloat(m[2]) : NaN; };

// ───────────────────────────── CIF ─────────────────────────────
function cifTokens(text) {
  const L = splitLines(text), out = [];
  for (let i = 0; i < L.length; i++) {
    let line = L[i];
    if (line[0] === ';') {
      const buf = [line.slice(1)]; const start = i + 1; i++;
      while (i < L.length && L[i][0] !== ';') buf.push(L[i++]);
      if (i >= L.length) throw new Error(`CIF: text field opened with ";" on line ${start} is never closed`);
      out.push({ v: buf.join('\n').trim(), q: true }); line = L[i].slice(1);
    }
    const n = line.length; let p = 0;
    while (p < n) {
      const ch = line[p];
      if (ch === ' ' || ch === '\t') { p++; continue; }
      if (ch === '#') break;
      if (ch === "'" || ch === '"') {
        let e = p + 1; while (e < n && !(line[e] === ch && (e + 1 >= n || line[e + 1] === ' ' || line[e + 1] === '\t'))) e++;
        if (e >= n) throw new Error(`CIF: unterminated quoted string on line ${i + 1}`);
        out.push({ v: line.slice(p + 1, e), q: true }); p = e + 1;
      } else { let e = p; while (e < n && line[e] !== ' ' && line[e] !== '\t') e++; out.push({ v: line.slice(p, e), q: false }); p = e; }
    }
  }
  return out;
}
const cifNum = v => { const t = String(v).replace(/\(\d*\)/g, '').trim(); if (t === '?' || t === '.' || t === '') return NaN; if (isNum(t)) return parseFloat(t); const f = frac(t); return f; };
function parseSymop(op) {
  const parts = String(op).toLowerCase().replace(/\s+/g, '').split(',');
  if (parts.length !== 3) throw new Error(`CIF: cannot parse symmetry operator "${op}" (expected three comma-separated terms)`);
  return parts.map(p => {
    const row = [0, 0, 0, 0], re = /([+-]?)(?:(\d*\.?\d+)(?:\/(\d*\.?\d+))?)?\*?([xyz])?/gy; let pos = 0, m;
    while (pos < p.length) {
      re.lastIndex = pos; m = re.exec(p);
      if (!m || !m[0].length || (m[2] === undefined && !m[4])) throw new Error(`CIF: cannot parse symmetry operator "${op}"`);
      const c = (m[1] === '-' ? -1 : 1) * (m[2] !== undefined ? parseFloat(m[2]) / (m[3] !== undefined ? parseFloat(m[3]) : 1) : 1);
      row[m[4] ? 'xyz'.indexOf(m[4]) : 3] += c; pos = re.lastIndex;
    }
    return row;
  });
}
/** CIF (also the common mmCIF _atom_site.* tags). Symmetry operators are applied and duplicate positions removed. */
export function parseCIF(text) {
  const tk = cifTokens(text), blocks = []; let cur = null;
  const tag = t => t.toLowerCase().replace(/\./g, '_');
  for (let i = 0; i < tk.length;) {
    const t = tk[i], lv = t.q ? '' : t.v.toLowerCase();
    if (lv.startsWith('data_')) { cur = { name: t.v.slice(5), items: {}, loops: [] }; blocks.push(cur); i++; continue; }
    if (!cur) { cur = { name: '', items: {}, loops: [] }; blocks.push(cur); }
    if (lv === 'loop_') {
      const tags = []; i++;
      while (i < tk.length && !tk[i].q && tk[i].v[0] === '_') tags.push(tag(tk[i++].v));
      const vals = [];
      while (i < tk.length) { const u = tk[i]; if (!u.q) { const w = u.v.toLowerCase(); if (u.v[0] === '_' || w === 'loop_' || w.startsWith('data_') || w.startsWith('save_') || w === 'stop_') break; } vals.push(u.v); i++; }
      if (tags.length) cur.loops.push({ tags, vals });
      continue;
    }
    if (!t.q && t.v[0] === '_') { const nx = tk[i + 1]; if (nx && (nx.q || (nx.v[0] !== '_' && nx.v.toLowerCase() !== 'loop_'))) { cur.items[tag(t.v)] = nx.v; i += 2; } else i++; continue; }
    i++;
  }
  const has = (lp, ...names) => names.some(n => lp.tags.includes(n));
  let blk = null, loop = null;
  for (const b of blocks) { const lp = b.loops.find(l => has(l, '_atom_site_fract_x', '_atom_site_cartn_x')); if (lp) { blk = b; loop = lp; break; } }
  if (!blk) throw new Error('CIF: no _atom_site loop with fractional or Cartesian coordinates found (this may be a reflection-data, dictionary or layout file rather than a crystal structure)');
  const nt = loop.tags.length;
  if (loop.vals.length % nt !== 0) throw new Error(`CIF: _atom_site loop has ${loop.vals.length} values, which is not a multiple of its ${nt} columns (a row is incomplete)`);
  const col = n => loop.tags.indexOf(n), useFrac = col('_atom_site_fract_x') >= 0, pre = useFrac ? '_atom_site_fract_' : '_atom_site_cartn_';
  const cx = col(pre + 'x'), cy = col(pre + 'y'), cz = col(pre + 'z'), cT = col('_atom_site_type_symbol'), cL = col('_atom_site_label') >= 0 ? col('_atom_site_label') : col('_atom_site_label_atom_id'), cO = col('_atom_site_occupancy'), cM = col('_atom_site_pdbx_pdb_model_num');
  if (cy < 0 || cz < 0) throw new Error('CIF: _atom_site loop lacks y or z coordinate columns');
  const it = blk.items, warnings = []; let cell = null, params = null;
  const cellTag = k => it['_cell_' + k];
  if (cellTag('length_a') !== undefined) {
    const v = ['length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma'].map((k, q) => { const raw = cellTag(k); if (raw === undefined) { if (q >= 3) return 90; throw new Error(`CIF: missing _cell_${k}`); } const x = cifNum(raw); if (!Number.isFinite(x) || x <= 0) throw new Error(`CIF: _cell_${k} is not a positive number ("${raw}")`); return x; });
    params = { a: v[0], b: v[1], c: v[2], alpha: v[3], beta: v[4], gamma: v[5] }; cell = latticeVectors(params);
  }
  if (useFrac && !cell) throw new Error('CIF: fractional coordinates given but the _cell_length_* / _cell_angle_* parameters are missing');
  let ops = null;
  for (const b of [blk]) for (const lp of b.loops) for (const nm of ['_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz']) { const c = lp.tags.indexOf(nm); if (c >= 0 && !ops) { ops = []; for (let r = 0; r + lp.tags.length <= lp.vals.length; r += lp.tags.length) ops.push(parseSymop(lp.vals[r + c])); } }
  if (!ops) { const single = it['_space_group_symop_operation_xyz'] || it['_symmetry_equiv_pos_as_xyz']; if (single) ops = [parseSymop(single)]; }
  const sg = it['_space_group_name_h-m_alt'] || it['_symmetry_space_group_name_h-m'] || it['_space_group_name_hall'] || null;
  const sgNo = cifNum(it['_space_group_it_number'] ?? it['_symmetry_int_tables_number'] ?? '');
  if (!ops || !useFrac) { if (useFrac && ((sg && !/^p\s*1$/i.test(sg.trim())) || sgNo > 1)) warnings.push(`Space group ${sg || sgNo} is declared but no symmetry operators are listed; only the atoms given in the file were used (the cell may be incomplete).`); ops = [[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]]; }
  const atoms = [], Mc = cell; let partial = false; const model0 = cM >= 0 ? loop.vals[cM] : null;
  for (let r = 0, row = 1; r < loop.vals.length; r += nt, row++) {
    if (cM >= 0 && loop.vals[r + cM] !== model0) continue;                                           // first model only
    const v = [cifNum(loop.vals[r + cx]), cifNum(loop.vals[r + cy]), cifNum(loop.vals[r + cz])];
    const label = cL >= 0 ? loop.vals[r + cL] : '', rawEl = cT >= 0 ? loop.vals[r + cT] : label;
    if (!v.every(Number.isFinite)) { if ([cx, cy, cz].some(c => /^[?.]$/.test(loop.vals[r + c]))) { warnings.push(`Atom site ${label || row} has undefined coordinates and was skipped.`); continue; } throw new Error(`CIF: non-numeric coordinate for atom site ${label || row}: "${loop.vals[r + cx]} ${loop.vals[r + cy]} ${loop.vals[r + cz]}"`); }
    const el = normEl(rawEl), oc = cO >= 0 ? cifNum(loop.vals[r + cO]) : 1, occ = Number.isFinite(oc) ? oc : 1; if (occ < 0.999) partial = true;
    if (!useFrac) { const a = { el, x: v[0], y: v[1], z: v[2] }; if (occ < 0.999) a.occ = occ; atoms.push(a); continue; }
    const orbit = [];
    for (const op of ops) {
      const f = op.map(rw => wrap1(rw[0] * v[0] + rw[1] * v[1] + rw[2] * v[2] + rw[3]));
      if (orbit.some(g => { const d = [0, 1, 2].map(k => { const t = Math.abs(g[k] - f[k]); return Math.min(t, 1 - t); }), c = rowMul(d, Mc); return Math.hypot(c[0], c[1], c[2]) < 0.05; })) continue;
      orbit.push(f);
    }
    for (const f of orbit) { const c = rowMul(f, Mc), a = { el, x: r6(c[0]), y: r6(c[1]), z: r6(c[2]), f }; if (label) a.label = label; if (occ < 0.999) a.occ = occ; atoms.push(a); }
  }
  if (partial) warnings.push('Partial site occupancies are present; every listed site is kept as an atom with its occupancy in `occ` (density and diffraction use it, visualisation shows all sites).');
  return finish('CIF', atoms, cell, it['_chemical_name_common'] || it['_chemical_name_mineral'] || it['_chemical_formula_sum'] || blk.name, { spacegroup: sg, nSymops: ops.length, warnings });
}

// ───────────────────────────── XYZ / extended XYZ ─────────────────────────────
export function parseXYZ(text) {
  const L = splitLines(text); let i = 0; while (i < L.length && !L[i].trim()) i++;
  if (i >= L.length) throw new Error('XYZ: file is empty');
  const frames = []; let first = null;
  const readAtoms = (from, n, sp, ps) => { const atoms = [];
    for (let q = 0; q < n; q++) { const ln = from + q, t = toks(L[ln] ?? ''); if (t.length < ps + 3) throw new Error(`XYZ: line ${ln + 1} should hold an element and three coordinates, got "${(L[ln] ?? '').trim().slice(0, 40)}"`);
      atoms.push({ el: normEl(t[sp]), x: num(t[ps], 'XYZ', 'x', ln + 1), y: num(t[ps + 1], 'XYZ', 'y', ln + 1), z: num(t[ps + 2], 'XYZ', 'z', ln + 1) }); }
    return atoms; };
  if (!isInt(L[i].trim())) {                                                                         // headerless coordinate list
    const rows = L.slice(i).map((l, k) => [l, k + i]).filter(([l]) => l.trim() && !/^\s*[#!]/.test(l));
    const atoms = rows.map(([l, k]) => { const t = toks(l); if (t.length < 4 || !isNum(t[1])) throw new Error(`XYZ: line ${k + 1} is neither an atom count nor an "element x y z" record: "${l.trim().slice(0, 40)}"`); return { el: normEl(t[0]), x: num(t[1], 'XYZ', 'x', k + 1), y: num(t[2], 'XYZ', 'y', k + 1), z: num(t[3], 'XYZ', 'z', k + 1) }; });
    return finish('XYZ', atoms, null, 'XYZ coordinates', { warnings: ['No atom-count header; every line was read as an atom.'] });
  }
  while (i < L.length) {
    while (i < L.length && !L[i].trim()) i++; if (i >= L.length) break;
    if (!isInt(L[i].trim())) { if (first) break; throw new Error(`XYZ: expected an atom count on line ${i + 1}`); }
    const n = parseInt(L[i], 10), comment = L[i + 1] ?? ''; if (n <= 0) throw new Error(`XYZ: atom count on line ${i + 1} must be positive`);
    const avail = L.slice(i + 2, i + 2 + n).filter(l => l.trim()).length; if (avail < n) throw new Error(`XYZ: header on line ${i + 1} declares ${n} atoms but only ${avail} coordinate lines follow`);
    let cell = null, sp = 0, ps = 1; const lm = comment.match(/Lattice\s*=\s*["'{]([^"'}]+)["'}]/i);
    if (lm) { const v = toks(lm[1]).map(t => num(t, 'extended XYZ', 'Lattice', i + 2)); if (v.length !== 9) throw new Error(`extended XYZ: Lattice needs 9 numbers, found ${v.length}`); cell = [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)]; }
    const pm = comment.match(/Properties\s*=\s*(\S+)/i);
    if (pm) { const p = pm[1].replace(/["']/g, '').split(':'); let c = 0, fs = -1, fp = -1; for (let q = 0; q + 2 < p.length; q += 3) { const nm = p[q].toLowerCase(), w = parseInt(p[q + 2], 10) || 1; if (fs < 0 && (nm === 'species' || nm === 'element' || nm === 'z')) fs = c; if (fp < 0 && (nm === 'pos' || nm === 'positions')) fp = c; c += w; } if (fp >= 0) { ps = fp; sp = fs >= 0 ? fs : 0; } }
    const atoms = readAtoms(i + 2, n, sp, ps), tm = comment.match(/comment\s*=\s*"([^"]*)"/i);
    const fr = { atoms, cell, title: tm ? tm[1] : (lm || pm ? '' : comment.trim()) }; if (!first) first = fr; if (frames.length < 500) frames.push({ t: frames.length, atoms, cell });
    i += 2 + n;
  }
  const extra = frames.length > 1 ? { frames, nFrames: frames.length } : {};
  return finish(first.cell ? 'extended XYZ' : 'XYZ', first.atoms, first.cell, first.title || 'XYZ structure', extra);
}

// ───────────────────────────── VASP POSCAR / CONTCAR ─────────────────────────────
export function parsePOSCAR(text) {
  const L = splitLines(text); if (L.length < 8) throw new Error('POSCAR: file is too short (needs comment, scale, three lattice vectors, counts, mode and coordinates)');
  const title = L[0].trim(), st = toks(L[1]), warnings = [];
  if (st.length !== 1 && st.length !== 3) throw new Error(`POSCAR: line 2 must hold the scaling factor (one number, or three), got "${L[1].trim()}"`);
  const sc = st.map(t => num(t, 'POSCAR', 'scaling factor', 2));
  let cell = [2, 3, 4].map(k => { const t = toks(L[k]); if (t.length < 3) throw new Error(`POSCAR: line ${k + 1} must hold a lattice vector (three numbers)`); return t.slice(0, 3).map(v => num(v, 'POSCAR', 'lattice vector', k + 1)); });
  let scale;
  if (sc.length === 3) { cell = cell.map(v => v.map((c, k) => c * sc[k])); scale = [1, 1, 1]; }
  else if (sc[0] < 0) { const s = Math.cbrt(-sc[0] / Math.abs(det3(cell))); cell = cell.map(v => v.map(c => c * s)); scale = [s, s, s]; }
  else { if (sc[0] === 0) throw new Error('POSCAR: scaling factor is zero'); cell = cell.map(v => v.map(c => c * sc[0])); scale = [sc[0], sc[0], sc[0]]; }
  let ln = 5, els = null, t = toks(L[ln].split('!')[0]);
  if (!t.length) throw new Error('POSCAR: line 6 (element symbols or atom counts) is empty');
  if (!t.every(isInt)) { els = t.map(s => normEl(s.split(/[_/]/)[0])); ln++; t = toks((L[ln] ?? '').split('!')[0]); }
  if (!t.length || !t.every(isInt)) throw new Error(`POSCAR: line ${ln + 1} must list the number of atoms per species, got "${(L[ln] ?? '').trim()}"`);
  const counts = t.map(Number); ln++;
  if (els && els.length !== counts.length) throw new Error(`POSCAR: ${els.length} element symbols but ${counts.length} atom counts`);
  if (!els) { const c = toks(title).map(s => s.replace(/[^A-Za-z]/g, '')).filter(s => ELEMENTS[s]); if (c.length === counts.length) { els = c; warnings.push('VASP 4 format without an element line: species were taken from the comment line.'); } else { els = counts.map(() => 'X'); warnings.push('VASP 4 format without an element line and no species in the comment line: elements are unknown (X).'); } }
  if (/^\s*[sS]/.test(L[ln] ?? '')) ln++;
  const mode = (L[ln] ?? '').trim(); if (!mode) throw new Error(`POSCAR: missing coordinate mode line (Direct or Cartesian) on line ${ln + 1}`);
  const cart = /^[cCkK]/.test(mode); ln++;
  const atoms = []; let selective = false;
  counts.forEach((cnt, si) => { for (let q = 0; q < cnt; q++, ln++) { const tt = toks(L[ln] ?? ''); if (tt.length < 3) throw new Error(`POSCAR: expected ${counts.reduce((a, b) => a + b, 0)} coordinate lines; line ${ln + 1} is missing or incomplete`);
    const v = tt.slice(0, 3).map(s => num(s, 'POSCAR', 'coordinate', ln + 1)), r = cart ? v.map((c, k) => c * scale[k]) : rowMul(v, cell), a = { el: els[si], x: r6(r[0]), y: r6(r[1]), z: r6(r[2]) };
    if (tt.length >= 6 && /^[TF]/i.test(tt[3])) { a.fixed = [tt[3], tt[4], tt[5]].map(s => /^F/i.test(s)); selective = true; } atoms.push(a); } });
  return finish('POSCAR', atoms, cell, title, { coordinateMode: cart ? 'Cartesian' : 'Direct', selectiveDynamics: selective, warnings });
}

// ───────────────────────────── XSF ─────────────────────────────
export function parseXSF(text) {
  const L = splitLines(text).map((l, i) => [l.replace(/#.*$/, '').trim(), i + 1]).filter(([l]) => l);
  let cell = null, conv = null, atoms = null, nFrames = 0, periodic = 'MOLECULE';
  const atomRow = ([l, ln]) => { const t = toks(l); if (t.length < 4) throw new Error(`XSF: line ${ln} should hold a species and three coordinates`); return { el: normEl(t[0]), x: num(t[1], 'XSF', 'x', ln), y: num(t[2], 'XSF', 'y', ln), z: num(t[3], 'XSF', 'z', ln) }; };
  for (let i = 0; i < L.length;) {
    const t = toks(L[i][0]), kw = t[0].toUpperCase();
    if (['CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE'].includes(kw)) { periodic = kw; i++; }
    else if (kw === 'PRIMVEC' || kw === 'CONVVEC') {
      if (i + 3 >= L.length) throw new Error(`XSF: ${kw} on line ${L[i][1]} is not followed by three vectors`);
      const v = [1, 2, 3].map(k => { const tt = toks(L[i + k][0]); if (tt.length < 3) throw new Error(`XSF: line ${L[i + k][1]} should hold a lattice vector`); return tt.slice(0, 3).map(s => num(s, 'XSF', 'lattice vector', L[i + k][1])); });
      if (kw === 'PRIMVEC') { if (!cell) cell = v; } else if (!conv) conv = v; i += 4;
    } else if (kw === 'PRIMCOORD') {
      const h = toks(L[i + 1]?.[0] ?? ''); if (!h.length || !isInt(h[0])) throw new Error(`XSF: PRIMCOORD on line ${L[i][1]} must be followed by "<number of atoms> 1"`);
      const n = +h[0]; if (i + 2 + n > L.length) throw new Error(`XSF: PRIMCOORD declares ${n} atoms but the file ends early`);
      if (!atoms) atoms = L.slice(i + 2, i + 2 + n).map(atomRow); nFrames++; i += 2 + n;
    } else if (kw === 'ATOMS') {
      const rows = []; i++; while (i < L.length) { const tt = toks(L[i][0]); if (tt.length >= 4 && isNum(tt[1]) && isNum(tt[2]) && isNum(tt[3]) && (isInt(tt[0]) || /^[A-Za-z]{1,2}\d*$/.test(tt[0]))) { rows.push(L[i]); i++; } else break; }
      if (!atoms) atoms = rows.map(atomRow); nFrames++;
    } else i++;
  }
  if (!atoms) throw new Error('XSF: no PRIMCOORD or ATOMS block found');
  return finish('XSF', atoms, cell || conv, 'XSF structure', { periodicity: periodic, nFrames });
}

// ───────────────────────────── PDB ─────────────────────────────
export function parsePDB(text) {
  const L = splitLines(text), atoms = [], warnings = []; let cell = null, title = '', models = 0, done = false;
  for (let i = 0; i < L.length; i++) {
    const l = L[i], rec = l.slice(0, 6).trim().toUpperCase();
    if (rec === 'TITLE' && !title) title = l.slice(10).trim();
    else if (rec === 'HEADER' && !title) title = l.slice(10, 50).trim();
    else if (rec === 'CRYST1') {
      let v = [[6, 15], [15, 24], [24, 33], [33, 40], [40, 47], [47, 54]].map(([a, b]) => l.slice(a, b).trim());
      if (!v.every(isNum)) v = toks(l).slice(1, 7); if (v.length < 6 || !v.every(isNum)) throw new Error(`PDB: malformed CRYST1 record on line ${i + 1}`);
      const p = v.map(Number); if (p.slice(0, 3).every(x => x > 1.001) && p.every(x => x > 0)) cell = latticeVectors({ a: p[0], b: p[1], c: p[2], alpha: p[3], beta: p[4], gamma: p[5] }); else warnings.push('CRYST1 holds a dummy unit cell (1 Å); treated as non-periodic.');
    } else if (rec === 'MODEL') models++;
    else if (rec === 'ENDMDL') done = true;
    else if ((rec === 'ATOM' || rec === 'HETATM') && !done) {
      let x = l.slice(30, 38), y = l.slice(38, 46), z = l.slice(46, 54);
      if (![x, y, z].every(isNum)) { const t = toks(l), k = t.findIndex((s, q) => q >= 3 && /\./.test(s) && isNum(s) && isNum(t[q + 1]) && isNum(t[q + 2])); if (k < 0) throw new Error(`PDB: cannot read coordinates on line ${i + 1}: "${l.trim().slice(0, 60)}"`); [x, y, z] = t.slice(k, k + 3); }
      const name = l.slice(12, 16), res = l.slice(17, 20); let el = l.slice(76, 78).trim();
      el = el && ELEMENTS[normEl(el)] && /^[A-Za-z]{1,2}$/.test(el) ? normEl(el) : (name[0] === ' ' || /\d/.test(name[0]) || rec === 'ATOM' && BIO.has(res.trim().toUpperCase())) ? normEl(name.replace(/[^A-Za-z]/g, '').slice(0, 1)) : elFromAtomName(name, res);
      const a = { el, x: +x, y: +y, z: +z, name: name.trim() }; if (res.trim()) a.res = res.trim(); atoms.push(a);
    }
  }
  if (!atoms.length) throw new Error('PDB: no ATOM or HETATM records found');
  if (models > 1) warnings.push(`File holds ${models} models; only the first was read.`);
  return finish('PDB', atoms, cell, title || 'PDB structure', { nModels: Math.max(1, models), warnings });
}

// ───────────────────────────── MDL MOL / SDF ─────────────────────────────
export function parseMOL(text) {
  const L = splitLines(text); if (L.some(l => /V3000/.test(l))) throw new Error('MOL: V3000 connection tables are not supported; export as V2000 (≤ 999 atoms) or as XYZ/PDB');
  let ci = L.findIndex(l => /V2000/.test(l)); if (ci < 0) ci = 3;
  const cl = L[ci]; if (cl === undefined) throw new Error('MOL: file is too short; the counts line (line 4) is missing');
  let na = parseInt(cl.slice(0, 3), 10), nb = parseInt(cl.slice(3, 6), 10); const ct = toks(cl);
  if (!(na > 0) || !Number.isFinite(nb) || !isInt(cl.slice(0, 3)) ) { na = parseInt(ct[0], 10); nb = parseInt(ct[1], 10); }
  if (!(na > 0) || !(nb >= 0)) throw new Error(`MOL: counts line (line ${ci + 1}) is malformed: "${cl.trim().slice(0, 40)}"`);
  if (L.length < ci + 1 + na) throw new Error(`MOL: counts line declares ${na} atoms but the file has only ${Math.max(0, L.length - ci - 1)} more lines`);
  const atoms = [], bonds = [];
  for (let q = 0; q < na; q++) { const ln = ci + 1 + q, l = L[ln]; let x = l.slice(0, 10), y = l.slice(10, 20), z = l.slice(20, 30), s = l.slice(31, 34).trim();
    if (![x, y, z].every(isNum) || !s) { const t = toks(l); if (t.length < 4) throw new Error(`MOL: atom line ${ln + 1} is malformed: "${l.trim().slice(0, 40)}"`); [x, y, z, s] = t; }
    atoms.push({ el: normEl(s), x: num(x, 'MOL', 'x', ln + 1), y: num(y, 'MOL', 'y', ln + 1), z: num(z, 'MOL', 'z', ln + 1) }); }
  for (let q = 0; q < nb; q++) { const ln = ci + 1 + na + q, l = L[ln]; if (l === undefined) throw new Error(`MOL: counts line declares ${nb} bonds but the file ends at line ${ln}`);
    let v = [l.slice(0, 3), l.slice(3, 6), l.slice(6, 9)]; if (!v.every(isInt)) v = toks(l).slice(0, 3); if (v.length < 3 || !v.every(isInt)) throw new Error(`MOL: bond line ${ln + 1} is malformed: "${l.trim().slice(0, 40)}"`);
    const a = +v[0] - 1, b = +v[1] - 1; if (a < 0 || b < 0 || a >= na || b >= na) throw new Error(`MOL: bond on line ${ln + 1} refers to an atom outside 1…${na}`); bonds.push([a, b, +v[2]]); }
  for (const l of L.slice(ci + 1 + na + nb)) { if (/^M\s+CHG/.test(l)) { const t = toks(l); for (let q = 3; q + 1 < t.length; q += 2) { const a = atoms[+t[q] - 1]; if (a) a.charge = +t[q + 1]; } } if (/^M\s+END/.test(l) || l.startsWith('$$$$')) break; }
  return finish('MOL', atoms, null, ci >= 3 ? L[ci - 3] : '', { bonds });
}
export function parseSDF(text) {
  const recs = String(text).replace(/^﻿/, '').split(/^\$\$\$\$.*$/m).filter(r => r.trim());
  if (!recs.length) throw new Error('SDF: file holds no molecule records');
  const first = recs[0].replace(/^(\r\n|\r|\n)/, ''), s = parseMOL(first), properties = {}, re = />\s*<([^>]+)>[^\n]*\r?\n([\s\S]*?)(?=\r?\n\s*\r?\n|\r?\n>|$)/g; let m;
  while ((m = re.exec(first))) properties[m[1].trim()] = m[2].trim();
  const warnings = s.warnings; if (recs.length > 1) warnings.push(`File holds ${recs.length} molecules; only the first was read.`);
  return { ...s, format: 'SDF', nMolecules: recs.length, properties, warnings };
}

// ───────────────────────────── Tripos MOL2 ─────────────────────────────
export function parseMOL2(text) {
  const L = splitLines(text), atoms = [], bonds = [], ids = {}; let sec = '', title = '', cell = null, mols = 0, nameNext = false;
  const ORD = { 1: 1, 2: 2, 3: 3, am: 1, ar: 1.5, du: 1, un: 1, nc: 0 };
  for (let i = 0; i < L.length; i++) {
    const l = L[i].trim(); if (!l || l[0] === '#') continue;
    if (l.startsWith('@<TRIPOS>')) { sec = l.slice(9).toUpperCase(); if (sec === 'MOLECULE') { mols++; nameNext = true; } continue; }
    if (mols > 1) continue;
    if (sec === 'MOLECULE' && nameNext) { title = l; nameNext = false; }
    else if (sec === 'ATOM') { const t = toks(l); if (t.length < 6) throw new Error(`MOL2: ATOM record on line ${i + 1} needs id, name, x, y, z and type`);
      const ty = t[5].split('.')[0], tyEl = /^[A-Za-z]{1,2}$/.test(ty) && ELEMENTS[ty[0].toUpperCase() + ty.slice(1).toLowerCase()] ? ty[0].toUpperCase() + ty.slice(1).toLowerCase() : null;
      const a = { el: tyEl || normEl(t[1]), x: num(t[2], 'MOL2', 'x', i + 1), y: num(t[3], 'MOL2', 'y', i + 1), z: num(t[4], 'MOL2', 'z', i + 1), name: t[1], type: t[5] }; if (t[8] !== undefined && isNum(t[8])) a.q = +t[8]; ids[t[0]] = atoms.length; atoms.push(a); }
    else if (sec === 'BOND') { const t = toks(l); if (t.length >= 4 && ids[t[1]] !== undefined && ids[t[2]] !== undefined) bonds.push([ids[t[1]], ids[t[2]], ORD[t[3].toLowerCase()] ?? 1]); }
    else if (sec === 'CRYSIN') { const t = toks(l); if (t.length >= 6 && t.slice(0, 6).every(isNum)) { const p = t.slice(0, 6).map(Number); cell = latticeVectors({ a: p[0], b: p[1], c: p[2], alpha: p[3], beta: p[4], gamma: p[5] }); } }
  }
  if (!atoms.length) throw new Error('MOL2: no @<TRIPOS>ATOM section with atoms found');
  const warnings = mols > 1 ? [`File holds ${mols} molecules; only the first was read.`] : [];
  return finish('MOL2', atoms, cell, title, { bonds, warnings });
}

// ───────────────────────────── LAMMPS data ─────────────────────────────
function elFromMass(m) { let best = 'X', bd = Infinity; for (const s of SYMBOLS) { const d = Math.abs(ELEMENTS[s].mass - m); if (d < bd) { bd = d; best = s; } } return bd <= Math.max(0.06, 0.012 * m) ? best : 'X'; }
/** LAMMPS data file, atom styles atomic / charge / molecular / full (from the "Atoms # style" hint, opts.style, or the column count). opts.elements maps type → element, e.g. ['Zn','O']. */
export function parseLAMMPSData(text, o = {}) {
  const raw = splitLines(text), title = (raw[0] || '').replace(/^#\s*/, '').trim(), warnings = [];
  const box = {}, tilt = [0, 0, 0], sections = {}; let nat = null, hint = '';
  let cur = null;
  for (let i = 1; i < raw.length; i++) {
    const full = raw[i], l = full.replace(/#.*$/, '').trim(); if (!l) continue;
    if (/^[A-Za-z]/.test(l)) { cur = { name: l.toLowerCase(), rows: [] }; sections[cur.name] = cur; if (cur.name === 'atoms') { const h = full.match(/#\s*(\w+)/); hint = h ? h[1].toLowerCase() : ''; } continue; }
    if (cur) { cur.rows.push([l, i + 1, (full.match(/#\s*(.*)$/) || [])[1] || '']); continue; }
    let m;
    if ((m = l.match(/^(\d+)\s+atoms\b/i))) nat = +m[1];
    else if ((m = l.match(/^(\S+)\s+(\S+)\s+([xyz])lo\s+[xyz]hi/i))) box[m[3].toLowerCase()] = [num(m[1], 'LAMMPS data', 'box bound', i + 1), num(m[2], 'LAMMPS data', 'box bound', i + 1)];
    else if ((m = l.match(/^(\S+)\s+(\S+)\s+(\S+)\s+xy\s+xz\s+yz/i))) { tilt[0] = num(m[1], 'LAMMPS data', 'tilt', i + 1); tilt[1] = num(m[2], 'LAMMPS data', 'tilt', i + 1); tilt[2] = num(m[3], 'LAMMPS data', 'tilt', i + 1); }
  }
  const A = sections.atoms; if (!A || !A.rows.length) throw new Error('LAMMPS data: no "Atoms" section found');
  for (const k of 'xyz') if (!box[k]) throw new Error(`LAMMPS data: header lacks the "${k}lo ${k}hi" box bounds`);
  const typeEl = {};
  if (sections.masses) for (const [l, ln, cm] of sections.masses.rows) { const t = toks(l); if (t.length < 2) throw new Error(`LAMMPS data: malformed Masses entry on line ${ln}`); const lab = toks(cm)[0]; typeEl[t[0]] = lab && ELEMENTS[normEl(lab)] && /^[A-Za-z]{1,2}\d*$/.test(lab) ? normEl(lab) : elFromMass(num(t[1], 'LAMMPS data', 'mass', ln)); }
  if (Array.isArray(o.elements)) o.elements.forEach((e, i) => { typeEl[i + 1] = normEl(e); }); else if (o.elements && typeof o.elements === 'object') for (const k in o.elements) typeEl[k] = normEl(o.elements[k]);
  let style = (o.style || hint || '').toLowerCase(); const nc = toks(A.rows[0][0]).length;
  if (!['atomic', 'charge', 'molecular', 'bond', 'angle', 'full'].includes(style)) {
    if (nc === 5 || nc === 8) style = 'atomic'; else if (nc === 7 || nc === 10) style = 'full';
    else if (nc === 6 || nc === 9) { style = A.rows.every(([l]) => { const t = toks(l); return isInt(t[1]) && isInt(t[2]); }) && !A.rows.every(([l]) => toks(l)[2] === '0') ? 'molecular' : 'charge'; warnings.push(`Atom style not stated; "${style}" assumed from the column layout.`); }
    else throw new Error(`LAMMPS data: cannot infer the atom style from ${nc} columns; state it as "Atoms # atomic|charge|molecular|full"`);
  }
  const lay = { atomic: [1, -1, -1, 2], charge: [1, 2, -1, 3], molecular: [2, -1, 1, 3], bond: [2, -1, 1, 3], angle: [2, -1, 1, 3], full: [2, 3, 1, 4] }[style];
  const org = [box.x[0], box.y[0], box.z[0]], rows = [];
  for (const [l, ln] of A.rows) { const t = toks(l); if (t.length < lay[3] + 3) throw new Error(`LAMMPS data: Atoms line ${ln} has ${t.length} columns; atom style "${style}" needs at least ${lay[3] + 3}`);
    const a = { el: typeEl[t[lay[0]]] || 'X', x: num(t[lay[3]], 'LAMMPS data', 'x', ln) - org[0], y: num(t[lay[3] + 1], 'LAMMPS data', 'y', ln) - org[1], z: num(t[lay[3] + 2], 'LAMMPS data', 'z', ln) - org[2], type: +t[lay[0]], id: +t[0] };
    if (lay[1] >= 0) a.q = num(t[lay[1]], 'LAMMPS data', 'charge', ln); if (lay[2] >= 0) a.mol = +t[lay[2]]; rows.push(a); }
  if (nat !== null && rows.length !== nat) throw new Error(`LAMMPS data: header declares ${nat} atoms but the Atoms section holds ${rows.length}`);
  rows.sort((a, b) => a.id - b.id);
  if (rows.some(a => a.el === 'X')) warnings.push('Some atom types have no Masses entry or an unrecognised mass; pass {elements:[…]} to name them.');
  const idIx = {}; rows.forEach((a, i) => { idIx[a.id] = i; }); const bonds = [];
  if (sections.bonds) for (const [l] of sections.bonds.rows) { const t = toks(l); if (t.length >= 4 && idIx[t[2]] !== undefined && idIx[t[3]] !== undefined) bonds.push([idIx[t[2]], idIx[t[3]], 1]); }
  const cell = [[box.x[1] - box.x[0], 0, 0], [tilt[0], box.y[1] - box.y[0], 0], [tilt[1], tilt[2], box.z[1] - box.z[0]]], extra = { atomStyle: style, origin: org, warnings }; if (bonds.length) extra.bonds = bonds;
  return finish('LAMMPS data', rows, cell, title, extra);
}

// ───────────────────────────── LAMMPS dump (text) ─────────────────────────────
/** LAMMPS text dump. opts.frames 'first' (default) | 'all'; opts.elements maps type → element. Returns the first frame plus frames:[{t, atoms, cell}]. */
export function parseLAMMPSDump(text, o = {}) {
  const L = splitLines(text), frames = [], all = o.frames === 'all', maxFrames = o.maxFrames > 0 ? o.maxFrames : 200, warnings = []; let i = 0, total = 0, t = 0, n = null, cell = null, org = [0, 0, 0];
  const mapEl = ty => Array.isArray(o.elements) ? normEl(o.elements[ty - 1] ?? '') : o.elements && o.elements[ty] ? normEl(o.elements[ty]) : 'X';
  while (i < L.length) {
    const l = L[i].trim(); if (!l.startsWith('ITEM:')) { i++; continue; }
    const item = l.slice(5).trim();
    if (item.startsWith('TIMESTEP')) { t = num(toks(L[i + 1] ?? '')[0], 'LAMMPS dump', 'timestep', i + 2); i += 2; }
    else if (item.startsWith('NUMBER OF ATOMS')) { n = num(toks(L[i + 1] ?? '')[0], 'LAMMPS dump', 'number of atoms', i + 2); i += 2; }
    else if (item.startsWith('BOX BOUNDS')) {
      const tri = /\bxy\b/.test(item), b = [1, 2, 3].map(k => { const tt = toks(L[i + k] ?? ''); if (tt.length < (tri ? 3 : 2)) throw new Error(`LAMMPS dump: malformed BOX BOUNDS on line ${i + k + 1}`); return tt.slice(0, 3).map(s => num(s, 'LAMMPS dump', 'box bound', i + k + 1)); });
      let [xlo, xhi, xy = 0] = b[0], [ylo, yhi, xz = 0] = b[1]; const [zlo, zhi, yz = 0] = b[2];
      if (tri) { xlo -= Math.min(0, xy, xz, xy + xz); xhi -= Math.max(0, xy, xz, xy + xz); ylo -= Math.min(0, yz); yhi -= Math.max(0, yz); }
      cell = [[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]]; org = [xlo, ylo, zlo]; i += 4;
    } else if (item.startsWith('ATOMS')) {
      if (n === null) throw new Error(`LAMMPS dump: "ITEM: ATOMS" on line ${i + 1} appears before "ITEM: NUMBER OF ATOMS"`);
      const cols = toks(item).slice(1), ix = nm => cols.indexOf(nm), cId = ix('id'), cTy = ix('type'), cEl = ix('element');
      const pick = ax => { for (const [sfx, sc] of [['', 0], ['u', 0], ['s', 1], ['su', 1]]) { const c = ix(ax + sfx); if (c >= 0) return [c, sc]; } return [-1, 0]; };
      const [cx, sx] = pick('x'), [cy, sy] = pick('y'), [cz, sz] = pick('z');
      if (cx < 0 || cy < 0 || cz < 0) throw new Error(`LAMMPS dump: the ATOMS columns (${cols.join(' ')}) hold no x/y/z, xs/ys/zs or xu/yu/zu coordinates`);
      if ((sx || sy || sz) && !cell) throw new Error('LAMMPS dump: scaled coordinates need a BOX BOUNDS item');
      if (L.length - i - 1 < n) throw new Error(`LAMMPS dump: frame at timestep ${t} declares ${n} atoms but the file ends early`);
      total++;
      if (frames.length < (all ? maxFrames : 1)) {
        const atoms = [];
        for (let q = 0; q < n; q++) { const ln = i + 1 + q, tt = toks(L[ln] ?? ''); if (tt.length < cols.length) throw new Error(`LAMMPS dump: atom line ${ln + 1} has ${tt.length} values, expected ${cols.length}`);
          let x = num(tt[cx], 'LAMMPS dump', 'x', ln + 1), y = num(tt[cy], 'LAMMPS dump', 'y', ln + 1), z = num(tt[cz], 'LAMMPS dump', 'z', ln + 1);
          if (sx || sy || sz) { const r = rowMul([x, y, z], cell); x = r[0]; y = r[1]; z = r[2]; } else { x -= org[0]; y -= org[1]; z -= org[2]; }
          const ty = cTy >= 0 ? +tt[cTy] : 1, a = { el: cEl >= 0 ? normEl(tt[cEl]) : mapEl(ty), x, y, z, type: ty }; if (cId >= 0) a.id = +tt[cId]; atoms.push(a); }
        if (cId >= 0) atoms.sort((a, b) => a.id - b.id);
        frames.push({ t, atoms, cell });
      }
      i += 1 + n; if (!all) break;
    } else i++;
  }
  if (!frames.length) throw new Error('LAMMPS dump: no "ITEM: ATOMS" frame found');
  if (total > frames.length) warnings.push(`Trajectory holds ${total} frames; ${frames.length} were read.`);
  if (frames[0].atoms.some(a => a.el === 'X')) warnings.push('Dump has no element column; pass {elements:[…]} (type → element) to name the atom types.');
  return finish('LAMMPS dump', frames[0].atoms, frames[0].cell, `LAMMPS dump, timestep ${frames[0].t}`, { frames, nFrames: frames.length, warnings });
}

// ───────────────────────────── GROMACS GRO ─────────────────────────────
export function parseGRO(text) {
  const L = splitLines(text); while (L.length && !L[L.length - 1].trim()) L.pop();
  if (L.length < 3) throw new Error('GRO: file is too short (title, atom count, atoms, box)');
  if (!isInt(L[1].trim())) throw new Error(`GRO: line 2 must hold the number of atoms, got "${L[1].trim().slice(0, 30)}"`);
  const n = +L[1].trim(); if (L.length < n + 3) throw new Error(`GRO: header declares ${n} atoms but the file has only ${L.length - 3} atom lines before the box line`);
  const atoms = [];
  for (let q = 0; q < n; q++) { const ln = 2 + q, l = L[ln], res = l.slice(5, 10).trim(), name = l.slice(10, 15).trim(), t = toks(l.slice(20));
    if (t.length < 3) throw new Error(`GRO: atom line ${ln + 1} lacks three coordinates`);
    atoms.push({ el: elFromAtomName(name, res), x: r6(num(t[0], 'GRO', 'x', ln + 1) * 10), y: r6(num(t[1], 'GRO', 'y', ln + 1) * 10), z: r6(num(t[2], 'GRO', 'z', ln + 1) * 10), name, res }); }
  const b = toks(L[n + 2]).map(s => num(s, 'GRO', 'box', n + 3)); let cell = null;
  if (b.length !== 3 && b.length !== 9) throw new Error(`GRO: box line must hold 3 or 9 numbers, found ${b.length}`);
  if (b[0] > 0 && b[1] > 0 && b[2] > 0) cell = (b.length === 9 ? [[b[0], b[3], b[4]], [b[5], b[1], b[6]], [b[7], b[8], b[2]]] : [[b[0], 0, 0], [0, b[1], 0], [0, 0, b[2]]]).map(v => v.map(c => r6(c * 10)));
  return finish('GRO', atoms, cell, L[0].trim());
}

// ───────────────────────────── CML ─────────────────────────────
export function parseCML(text) {
  const s = String(text); if (!/<(\w+:)?(cml|molecule|atomArray|atom)\b/.test(s)) throw new Error('CML: no <cml>, <molecule> or <atom> elements found');
  const attrs = t => { const o = {}, re = /([\w:.-]+)\s*=\s*(?:"([^"]*)"|'([^']*)')/g; let m; while ((m = re.exec(t))) o[m[1].replace(/^\w+:/, '')] = m[2] ?? m[3]; return o; };
  const P6 = {}; const cr = s.match(/<(?:\w+:)?crystal\b[\s\S]*?<\/(?:\w+:)?crystal>/);
  if (cr) { const re = /<(?:\w+:)?scalar\b([^>]*)>([^<]*)</g; let m; while ((m = re.exec(cr[0]))) { const a = attrs(m[1]), key = ((a.title || a.dictRef || a.id || '').match(/(alpha|beta|gamma|[abc])$/i) || [])[1]; if (key) P6[key.toLowerCase()] = cifNum(m[2]); }
    const re2 = /<(?:\w+:)?cellParameter\b([^>]*)>([^<]*)</g; while ((m = re2.exec(cr[0]))) { const a = attrs(m[1]), v = toks(m[2]).map(cifNum); if (a.type === 'length') [P6.a, P6.b, P6.c] = v; else if (a.type === 'angle') [P6.alpha, P6.beta, P6.gamma] = v; } }
  let cell = null; if (['a', 'b', 'c'].every(k => P6[k] > 0)) cell = latticeVectors({ a: P6.a, b: P6.b, c: P6.c, alpha: P6.alpha || 90, beta: P6.beta || 90, gamma: P6.gamma || 90 });
  const atoms = [], ids = {}, add = (a, where) => {
    let r; const el = normEl(a.elementType ?? '');
    if (a.x3 !== undefined) r = [a.x3, a.y3, a.z3].map(v => num(v, 'CML', 'x3/y3/z3', 0));
    else if (a.xFract !== undefined) { if (!cell) throw new Error('CML: fractional coordinates given without a <crystal> cell'); r = rowMul([a.xFract, a.yFract, a.zFract].map(v => num(String(v).replace(/\(\d*\)/, ''), 'CML', 'fractional coordinate', 0)), cell); }
    else if (a.x2 !== undefined) r = [num(a.x2, 'CML', 'x2', 0), num(a.y2, 'CML', 'y2', 0), 0];
    else throw new Error(`CML: atom ${a.id || where} has no x3/y3/z3, xFract/yFract/zFract or x2/y2 coordinates`);
    if (a.id) ids[a.id] = atoms.length; const at = { el, x: r[0], y: r[1], z: r[2] }; if (a.formalCharge) at.charge = +a.formalCharge; atoms.push(at);
  };
  let m; const reA = /<(?:\w+:)?atom\b([^>]*)>/g; while ((m = reA.exec(s))) add(attrs(m[1]), atoms.length + 1);
  if (!atoms.length) { const reArr = /<(?:\w+:)?atomArray\b([^>]*)>/g; while ((m = reArr.exec(s))) { const a = attrs(m[1]); if (!a.elementType) continue; const col = k => a[k] !== undefined ? toks(a[k]) : null, E = col('elementType'), I = col('atomID'), keys = ['x3', 'y3', 'z3', 'xFract', 'yFract', 'zFract', 'x2', 'y2'], C = {}; for (const k of keys) C[k] = col(k);
      E.forEach((e, q) => { const o = { elementType: e, id: I ? I[q] : undefined }; for (const k of keys) if (C[k]) o[k] = C[k][q]; add(o, q + 1); }); } }
  if (!atoms.length) throw new Error('CML: no atoms with an elementType attribute found');
  const bonds = [], reB = /<(?:\w+:)?bond\b([^>]*)>/g, ORD = { S: 1, D: 2, T: 3, A: 1.5 };
  while ((m = reB.exec(s))) { const a = attrs(m[1]), r = toks(a.atomRefs2 || ''); if (r.length === 2 && ids[r[0]] !== undefined && ids[r[1]] !== undefined) bonds.push([ids[r[0]], ids[r[1]], ORD[a.order] ?? (parseFloat(a.order) || 1)]); }
  const mt = s.match(/<(?:\w+:)?molecule\b([^>]*)>/), ma = mt ? attrs(mt[1]) : {}, nm = s.match(/<(?:\w+:)?name\b[^>]*>([^<]*)</);
  return finish('CML', atoms, cell, (nm && nm[1]) || ma.title || ma.id || 'CML structure', { bonds });
}

// ───────────────────────────── Gaussian input ─────────────────────────────
export function parseGaussianInput(text) {
  const L = splitLines(text); let i = 0; while (i < L.length && !L[i].trim()) i++;
  while (i < L.length && (L[i].trim().startsWith('%') || L[i].trim().startsWith('!'))) i++;
  if (i >= L.length || !L[i].trim().startsWith('#')) throw new Error('Gaussian input: no route section (a line starting with "#") found after the Link 0 commands');
  let route = ''; while (i < L.length && L[i].trim()) route += ' ' + L[i++].trim(); route = route.trim();
  if (/geom\s*=?\s*\(?\s*(allcheck|check)/i.test(route)) throw new Error('Gaussian input: the geometry is read from a checkpoint file (geom=check); the input holds no coordinates');
  while (i < L.length && !L[i].trim()) i++; let title = ''; while (i < L.length && L[i].trim()) title += ' ' + L[i++].trim();
  while (i < L.length && !L[i].trim()) i++;
  const cm = (L[i] ?? '').trim().split(/[\s,]+/).filter(Boolean);
  if (cm.length < 2 || !cm.every(isInt)) throw new Error(`Gaussian input: expected the "charge multiplicity" line on line ${i + 1}, got "${(L[i] ?? '').trim().slice(0, 40)}"`);
  i++; const atoms = [], tv = [];
  for (; i < L.length && L[i].trim(); i++) {
    const t = L[i].trim().split(/[\s,]+/).filter(Boolean), sym = t[0].split(/[(\-]/)[0], nums = t.slice(1);
    const k = nums.length >= 4 && /^(0|-1)$/.test(nums[0]) && isNum(nums[3]) ? 1 : 0, nNum = nums.filter(isNum).length;      // optional freeze code before x y z
    if (nums.length < k + 3 || ![0, 1, 2].every(q => isNum(nums[k + q])) || (k === 0 && isInt(nums[0]) && nNum % 2 === 0))
      throw new Error(`Gaussian input: line ${i + 1} is not a Cartesian atom record ("${L[i].trim().slice(0, 40)}"). Z-matrix (internal coordinate) geometries are not supported; supply Cartesian coordinates`);
    const r = [0, 1, 2].map(q => parseFloat(nums[k + q]));
    if (/^tv$/i.test(sym)) { tv.push(r); continue; } if (/^(x|bq|q)$/i.test(sym) || sym === '0' || sym === '-1') continue;
    atoms.push({ el: normEl(sym), x: r[0], y: r[1], z: r[2] });
  }
  if (!atoms.length) throw new Error('Gaussian input: no Cartesian atom records found after the charge/multiplicity line');
  const warnings = tv.length && tv.length < 3 ? [`${tv.length} translation vector(s) (1-D/2-D periodicity) found; a 3-D cell is needed, so the structure is treated as non-periodic.`] : [];
  return finish('Gaussian input', atoms, tv.length === 3 ? tv : null, title.trim(), { route, charge: +cm[0], multiplicity: +cm[1], translationVectors: tv, warnings });
}

// ───────────────────────────── Quantum ESPRESSO pw.x input ─────────────────────────────
export function parseQEInput(text) {
  const L = splitLines(text).map(l => l.replace(/!.*$/, '').replace(/^\s*#.*$/, '')), body = L.join('\n');
  const nl = re => { const m = body.match(re); return m ? parseFloat(m[1].replace(/[dD]/, 'e')) : null; };
  const ibrav = nl(/\bibrav\s*=\s*(-?\d+)/i), nat = nl(/\bnat\s*=\s*(\d+)/i), A = nl(/(?:^|[\s,])a\s*=\s*([-+\d.eEdD]+)/im), Bq = nl(/(?:^|[\s,])b\s*=\s*([-+\d.eEdD]+)/im), Cq = nl(/(?:^|[\s,])c\s*=\s*([-+\d.eEdD]+)/im);
  const cd = {}; { const re = /celldm\s*\(\s*(\d)\s*\)\s*=\s*([-+\d.eEdD]+)/gi; let m; while ((m = re.exec(body))) cd[m[1]] = parseFloat(m[2].replace(/[dD]/, 'e')); }
  const alat = cd[1] ? cd[1] * BOHR : A ? A : null;
  const val = (t, what, ln) => { if (isNum(t)) return num(t, 'QE input', what, ln); const f = frac(t); if (Number.isFinite(f)) return f; throw new Error(`QE input: expected a number for ${what} on line ${ln}, got "${t}"`); };
  const CARD = /^\s*(CELL_PARAMETERS|ATOMIC_POSITIONS|ATOMIC_SPECIES|K_POINTS|CONSTRAINTS|OCCUPATIONS|ATOMIC_FORCES|ATOMIC_VELOCITIES|HUBBARD|SOLVENTS|ADDITIONAL_K_POINTS)\b\s*[{(]?\s*(\w+)?/i;
  let cell = null, posUnit = null, rawPos = null, hadCard = false;
  for (let i = 0; i < L.length; i++) {
    const m = L[i].match(CARD); if (!m) continue; const kw = m[1].toUpperCase(), unit = (m[2] || '').toLowerCase();
    if (kw === 'CELL_PARAMETERS') {
      const v = [1, 2, 3].map(k => { const t = toks(L[i + k] ?? ''); if (t.length < 3) throw new Error(`QE input: CELL_PARAMETERS needs three vectors; line ${i + k + 1} is incomplete`); return t.slice(0, 3).map(s => val(s, 'cell vector', i + k + 1)); });
      const u = unit || (alat ? 'alat' : 'bohr'); let f;
      if (u === 'angstrom') f = 1; else if (u === 'bohr') f = BOHR; else if (u === 'alat') { if (!alat) throw new Error('QE input: CELL_PARAMETERS in alat units but neither celldm(1) nor A is set'); f = alat; } else throw new Error(`QE input: unknown CELL_PARAMETERS unit "${u}"`);
      cell = v.map(r => r.map(c => c * f)); i += 3;
    } else if (kw === 'ATOMIC_POSITIONS') {
      hadCard = true; posUnit = unit || 'alat'; rawPos = [];
      for (let q = i + 1; q < L.length; q++) { if (!L[q].trim()) { if (rawPos.length) break; continue; } if (CARD.test(L[q]) || /^\s*[&/]/.test(L[q])) break; const t = toks(L[q]); if (t.length < 4) throw new Error(`QE input: ATOMIC_POSITIONS line ${q + 1} needs a species label and three coordinates`); rawPos.push({ el: normEl(t[0]), label: t[0], v: t.slice(1, 4).map(s => val(s, 'atomic position', q + 1)) }); i = q; }
    }
  }
  if (!hadCard || !rawPos.length) throw new Error('QE input: no ATOMIC_POSITIONS card with atoms found');
  if (!cell && ibrav !== null && ibrav !== 0) {
    if (!alat) throw new Error('QE input: ibrav ≠ 0 needs celldm(1) or A to define the lattice');
    const ca = cd[3] ?? (Cq && A ? Cq / A : null), ba = cd[2] ?? (Bq && A ? Bq / A : null), a = alat;
    if (ibrav === 1) cell = [[a, 0, 0], [0, a, 0], [0, 0, a]];
    else if (ibrav === 2) cell = [[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]];
    else if (ibrav === 3) cell = [[a / 2, a / 2, a / 2], [-a / 2, a / 2, a / 2], [-a / 2, -a / 2, a / 2]];
    else if (ibrav === 4 && ca) cell = [[a, 0, 0], [-a / 2, a * Math.sqrt(3) / 2, 0], [0, 0, a * ca]];
    else if (ibrav === 6 && ca) cell = [[a, 0, 0], [0, a, 0], [0, 0, a * ca]];
    else if (ibrav === 8 && ca && ba) cell = [[a, 0, 0], [0, a * ba, 0], [0, 0, a * ca]];
    else throw new Error(`QE input: ibrav = ${ibrav} without a CELL_PARAMETERS card is not supported here (supported: 1, 2, 3, 4, 6, 8 with the needed celldm/A,B,C); re-export with ibrav = 0 and explicit CELL_PARAMETERS`);
  }
  if (posUnit === 'crystal_sg') throw new Error('QE input: ATOMIC_POSITIONS crystal_sg (space-group Wyckoff input) is not supported; use crystal or angstrom positions');
  if (!['alat', 'bohr', 'angstrom', 'crystal'].includes(posUnit)) throw new Error(`QE input: unknown ATOMIC_POSITIONS unit "${posUnit}"`);
  if (posUnit === 'crystal' && !cell) throw new Error('QE input: crystal (fractional) positions given but no cell could be determined');
  if (posUnit === 'alat' && !alat) throw new Error('QE input: positions in alat units but neither celldm(1) nor A is set');
  if (nat !== null && nat !== rawPos.length) throw new Error(`QE input: nat = ${nat} but ${rawPos.length} atomic positions were found`);
  const atoms = rawPos.map(p => { const r = posUnit === 'crystal' ? rowMul(p.v, cell) : p.v.map(c => c * (posUnit === 'angstrom' ? 1 : posUnit === 'bohr' ? BOHR : alat)); return { el: p.el, x: r6(r[0]), y: r6(r[1]), z: r6(r[2]), label: p.label }; });
  return finish('QE input', atoms, cell, 'Quantum ESPRESSO input', { ibrav, positionUnit: posUnit });
}

// ───────────────────────────── CASTEP .cell ─────────────────────────────
export function parseCASTEPCell(text) {
  const body = splitLines(text).map(l => l.replace(/[#!;].*$/, '')).join('\n'), blocks = {}, re = /%block\s+(\w+)\s*\n([\s\S]*?)%endblock\s+\w+/gi; let m;
  while ((m = re.exec(body))) blocks[m[1].toLowerCase()] = m[2].split('\n').map(l => l.trim()).filter(Boolean);
  const UNIT = { ang: 1, angstrom: 1, bohr: BOHR, a0: BOHR, nm: 10, m: 1e10, cm: 1e8 };
  const unitOf = rows => { if (rows.length && toks(rows[0]).length === 1 && !isNum(rows[0])) { const u = rows[0].toLowerCase(); if (!(u in UNIT)) throw new Error(`CASTEP cell: unknown length unit "${rows[0]}"`); rows.shift(); return UNIT[u]; } return 1; };
  let cell = null;
  if (blocks.lattice_cart) { const rows = [...blocks.lattice_cart], f = unitOf(rows); if (rows.length < 3) throw new Error('CASTEP cell: LATTICE_CART needs three lattice vectors'); cell = rows.slice(0, 3).map(r => { const t = toks(r); if (t.length < 3) throw new Error('CASTEP cell: a LATTICE_CART vector has fewer than three components'); return t.slice(0, 3).map(s => num(s, 'CASTEP cell', 'lattice vector', 0) * f); }); }
  else if (blocks.lattice_abc) { const rows = [...blocks.lattice_abc], f = unitOf(rows); if (rows.length < 2) throw new Error('CASTEP cell: LATTICE_ABC needs "a b c" and "alpha beta gamma" lines'); const l = toks(rows[0]).map(s => num(s, 'CASTEP cell', 'a b c', 0) * f), an = toks(rows[1]).map(s => num(s, 'CASTEP cell', 'angles', 0)); if (l.length < 3 || an.length < 3) throw new Error('CASTEP cell: LATTICE_ABC lines need three numbers each'); cell = latticeVectors({ a: l[0], b: l[1], c: l[2], alpha: an[0], beta: an[1], gamma: an[2] }); }
  else throw new Error('CASTEP cell: no %BLOCK LATTICE_CART or %BLOCK LATTICE_ABC found');
  const fracB = blocks.positions_frac, absB = blocks.positions_abs;
  if (!fracB && !absB) throw new Error('CASTEP cell: no %BLOCK POSITIONS_FRAC or %BLOCK POSITIONS_ABS found');
  const rows = [...(fracB || absB)], f = fracB ? 1 : unitOf(rows);
  const atoms = rows.map((r, k) => { const t = toks(r); if (t.length < 4) throw new Error(`CASTEP cell: position entry ${k + 1} needs a species and three coordinates ("${r.slice(0, 40)}")`); const v = t.slice(1, 4).map(s => num(s, 'CASTEP cell', 'position', 0)), c = fracB ? rowMul(v, cell) : v.map(x => x * f); return { el: normEl(t[0].split(':')[0]), x: r6(c[0]), y: r6(c[1]), z: r6(c[2]) }; });
  return finish('CASTEP cell', atoms, cell, 'CASTEP cell', { positions: fracB ? 'fractional' : 'absolute' });
}
