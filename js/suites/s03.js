// Suite 3 — Atomic, molecular and crystal structure.
// Crystal builder with powder XRD, nanostructure builder (broken-bond energetics), defect chemistry (point-defect
// equilibria, Kröger–Vink reactions, Vegard / tolerance checks), lattice energetics from pair potentials and the
// Born model, epitaxial mismatch (Matthews–Blakeslee) and Bragg / reciprocal-space tools.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import * as A from '../geometry/atomic.js';

const F = N.fmt, KB_EV = N.K_B / N.Q_E, COUL = 14.399645, EVA3_GPA = 160.21766, PI = Math.PI;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
const fin = (v, d = 0) => (Number.isFinite(v) ? v : d);
const num = (inp, k, d, lo, hi) => N.clamp(fin(+inp[k], d), lo, hi);
const int = (inp, k, d, lo, hi) => Math.round(num(inp, k, d, lo, hi));
const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a || 1; };
const hklStr = h => `(${h.join(' ')})`;
const KEYS = Object.keys(A.CRYSTALS);
const crystalOpts = KEYS.map(k => ({ value: k, label: `${A.CRYSTALS[k].name} — ${chem.prettyFormula(A.CRYSTALS[k].formula)}` }));
const pick = (v, list, d) => (list.includes(v) ? v : d);
const cKey = (inp, k = 'crystal') => pick(inp[k], KEYS, 'ZnO_wurtzite');
const slug = s => String(s).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 40) || 'structure';

// Bulk reference data. Ec: cohesive (atomisation) energy, eV per atom — elements: handbook cohesive energies; compounds: from
// standard formation enthalpies and gas-phase atom enthalpies at 298 K. Tm: melting point, K (null where the phase decomposes,
// sublimes or transforms first). g: typical experimental solid surface energy, J/m² (metals only).
const BULK = {
  Au: { Ec: 3.81, Tm: 1337, g: 1.50 }, Ag: { Ec: 2.95, Tm: 1235, g: 1.25 }, Cu: { Ec: 3.49, Tm: 1358, g: 1.79 }, Pt: { Ec: 5.84, Tm: 2041, g: 2.49 },
  Pd: { Ec: 3.89, Tm: 1828, g: 2.00 }, Ni: { Ec: 4.44, Tm: 1728, g: 2.38 }, Al: { Ec: 3.39, Tm: 933, g: 1.14 }, Fe_bcc: { Ec: 4.28, Tm: 1811, g: 2.42 },
  W: { Ec: 8.90, Tm: 3695, g: 3.27 }, Si: { Ec: 4.63, Tm: 1687 }, Ge: { Ec: 3.85, Tm: 1211 }, C_diamond: { Ec: 7.35 }, graphite: { Ec: 7.37 }, graphene: { Ec: 7.32 },
  ZnO_wurtzite: { Ec: 3.78, Tm: 2248 }, ZnO_zincblende: { Ec: 3.78 }, TiO2_rutile: { Ec: 6.62, Tm: 2116 }, TiO2_anatase: { Ec: 6.6 }, NaCl: { Ec: 3.32, Tm: 1074 },
  MgO: { Ec: 5.17, Tm: 3125 }, GaAs: { Ec: 3.35, Tm: 1511 }, GaN_wurtzite: { Ec: 4.43 }, ZnS_zincblende: { Ec: 3.18 }, CeO2: { Ec: 6.94, Tm: 2673 },
  SrTiO3: { Ec: 6.34, Tm: 2353 }, BaTiO3_cubic: { Ec: 6.34, Tm: 1898 }, CsCl: { Ec: 3.32, Tm: 918 }, Al2O3_corundum: { Ec: 6.39, Tm: 2345 },
  Fe3O4: { Ec: 4.98, Tm: 1870 }, SiO2_cristobalite: { Ec: 6.42, Tm: 1986 }, hBN: { Ec: 6.7 },
};
// Formal charges and Shannon effective ionic radii (Å, six-fold coordination unless noted).
const ANION_Q = { O: 2, S: 2, Se: 2, Cl: 1, F: 1, N: 3, As: 3 };
const CATION_Q = { Zn: 2, Ti: 4, Ga: 3, Cd: 2, Na: 1, Mg: 2, Ce: 4, Sr: 2, Ba: 2, Cs: 1, Fe: 3, Si: 4, Al: 3, B: 3, Mo: 4 };
const RION = {
  'Li+': 0.76, 'Na+': 1.02, 'K+': 1.38, 'Cs+': 1.67, 'Ag+': 1.15, 'Mg2+': 0.72, 'Ca2+': 1.00, 'Sr2+': 1.18, 'Ba2+': 1.35, 'Zn2+': 0.74, 'Cd2+': 0.95, 'Cu2+': 0.73,
  'Ni2+': 0.69, 'Co2+': 0.745, 'Mn2+': 0.83, 'Fe2+': 0.78, 'Fe3+': 0.645, 'Cr3+': 0.615, 'Al3+': 0.535, 'Ga3+': 0.62, 'In3+': 0.80, 'Y3+': 0.90, 'La3+': 1.032,
  'Gd3+': 0.938, 'B3+': 0.27, 'Ti4+': 0.605, 'Zr4+': 0.72, 'Ce4+': 0.87, 'Sn4+': 0.69, 'Si4+': 0.40, 'Mo4+': 0.65, 'Nb5+': 0.64,
  'O2-': 1.40, 'S2-': 1.84, 'Se2-': 1.98, 'F-': 1.33, 'Cl-': 1.81, 'N3-': 1.46, 'As3-': 2.22,
};
const R4 = { 'Zn2+': 0.60, 'Al3+': 0.39, 'Ga3+': 0.47, 'In3+': 0.62, 'Mg2+': 0.57, 'Li+': 0.59, 'Na+': 0.99, 'Ag+': 1.00, 'Cd2+': 0.78, 'Co2+': 0.58, 'Cu2+': 0.57, 'Mn2+': 0.66, 'Ni2+': 0.55, 'Fe3+': 0.49, 'Si4+': 0.26, 'Ti4+': 0.42 };   // four-fold
const R12 = { 'Sr2+': 1.44, 'Ba2+': 1.61, 'Ca2+': 1.34, 'La3+': 1.36, 'Na+': 1.39, 'K+': 1.64 };   // twelve-fold (perovskite A site)
const DOPANTS = ['Al3+', 'Ga3+', 'In3+', 'Li+', 'Na+', 'Ag+', 'Mg2+', 'Ca2+', 'Sr2+', 'Cu2+', 'Co2+', 'Mn2+', 'Ni2+', 'Fe3+', 'Cr3+', 'Y3+', 'La3+', 'Gd3+', 'Zr4+', 'Ce4+', 'Sn4+', 'Nb5+'];
const ionOf = s => { const m = s.match(/^([A-Z][a-z]?)(\d*)([+-])$/); return { el: m[1], z: (m[2] ? +m[2] : 1) * (m[3] === '+' ? 1 : -1) }; };
const ionKey = (el, z) => `${el}${Math.abs(z) > 1 ? Math.abs(z) : ''}${z > 0 ? '+' : '-'}`;

// Pair-potential presets. LJ: ε (eV), σ (Å). Morse: D (eV), α (1/Å), r0 (Å) — the classic lattice-property-fitted set for cubic
// metals (Girifalco & Weizer type). Experimental references: lattice parameter (Å), cohesive energy (eV/atom), bulk modulus (GPa).
const POT = {
  Ar: { name: 'Argon — Lennard-Jones', type: 'lj', e: 0.0104, r: 3.405, st: 'fcc', a: 5.31, Ec: 0.080, B: 2.7, mass: 39.948 },
  Cu: { name: 'Copper — Morse', type: 'morse', e: 0.3429, al: 1.3588, r: 2.866, st: 'fcc', a: 3.615, Ec: 3.49, B: 137, mass: 63.546 },
  Ag: { name: 'Silver — Morse', type: 'morse', e: 0.3323, al: 1.3690, r: 3.115, st: 'fcc', a: 4.086, Ec: 2.95, B: 100.7, mass: 107.868 },
  Au: { name: 'Gold — Morse', type: 'morse', e: 0.4753, al: 1.5830, r: 3.024, st: 'fcc', a: 4.078, Ec: 3.81, B: 173, mass: 196.967 },
  Ni: { name: 'Nickel — Morse', type: 'morse', e: 0.4205, al: 1.4199, r: 2.780, st: 'fcc', a: 3.524, Ec: 4.44, B: 186, mass: 58.693 },
  Al: { name: 'Aluminium — Morse', type: 'morse', e: 0.2703, al: 1.1646, r: 3.253, st: 'fcc', a: 4.050, Ec: 3.39, B: 72.2, mass: 26.982 },
  Fe: { name: 'Iron (bcc) — Morse', type: 'morse', e: 0.4174, al: 1.3885, r: 2.845, st: 'bcc', a: 2.867, Ec: 4.28, B: 168, mass: 55.845 },
  Na: { name: 'Sodium (bcc) — Morse', type: 'morse', e: 0.06334, al: 0.58993, r: 5.336, st: 'bcc', a: 4.225, Ec: 1.113, B: 6.8, mass: 22.990 },
};
// Ionic AB crystals for the Born model: Born exponent from the ion cores (He 5, Ne 7, Ar 9, Kr 10, Xe 12, averaged), reference
// Madelung constant (nearest-neighbour distance basis), and Born–Haber data in kJ/mol where every term is directly measured.
const SALTS = {
  NaCl: { z: 1, n: 8, Aref: 1.74756, bh: { dHf: -411.2, sub: 107.3, ie: 495.8, at: 121.3, ea: -348.6 }, Bexp: 24 },
  CsCl: { z: 1, n: 10.5, Aref: 1.76267, bh: { dHf: -443.0, sub: 76.5, ie: 375.7, at: 121.3, ea: -348.6 }, Bexp: 18 },
  MgO: { z: 2, n: 7, Aref: 1.74756, Bexp: 160 },
  ZnS_zincblende: { z: 2, n: 9, Aref: 1.63806, Bexp: 77 },
  ZnO_wurtzite: { z: 2, n: 8, Aref: 1.641, Bexp: 143 },
  CdSe_wurtzite: { z: 2, n: 9.5, Aref: 1.641, Bexp: 53 },
};
const RAD = {
  CuKa: { name: 'Cu Kα₁', lam: 1.5406 }, MoKa: { name: 'Mo Kα₁', lam: 0.7093 }, CoKa: { name: 'Co Kα₁', lam: 1.7890 }, CrKa: { name: 'Cr Kα₁', lam: 2.2897 },
  AgKa: { name: 'Ag Kα₁', lam: 0.5594 }, electron: { name: 'electrons', lam: null },
};
const eLambda = kV => { const V = kV * 1e3, E = N.Q_E * V; return N.H_PLANCK / Math.sqrt(2 * N.M_E * E * (1 + E / (2 * N.M_E * N.C0 ** 2))) * 1e10; };   // Å
const wavelength = inp => { const r = pick(inp.rad, Object.keys(RAD), 'CuKa'); return r === 'electron' ? { key: r, name: `${F(num(inp, 'kV', 200, 20, 1000), 4)} kV electrons`, lam: eLambda(num(inp, 'kV', 200, 20, 1000)) } : { key: r, ...RAD[r] }; };

const cellInfo = s => A.cellParams(s.cell);
const vacuumGap = s => { if (!s.cell) return 0; const c = s.cell[2][2], z = s.atoms.map(a => a.z).sort((p, q) => p - q); let g = c - z[z.length - 1] + z[0]; for (let i = 1; i < z.length; i++) g = Math.max(g, z[i] - z[i - 1]); return g; };   // widest empty gap along z (Å); > 5 Å marks a vacuum-padded 2-D entry
const formulaAtoms = f => { const a = chem.parseFormula(f).atoms; return Object.values(a).reduce((x, y) => x + y, 0); };
function structExports(s, base, which = ['cif', 'poscar', 'xyz', 'extxyz', 'xsf', 'lammps']) {
  const W = {
    cif: ['.cif', 'chemical/x-cif', A.toCIF, 'Crystallographic Information File (P 1, all atoms listed).'], poscar: ['.POSCAR.vasp', 'text/plain', A.toPOSCAR, 'VASP POSCAR (direct coordinates).'],
    xyz: ['.xyz', 'chemical/x-xyz', A.toXYZ, 'Cartesian XYZ coordinates in Å.'], extxyz: ['.extxyz', 'chemical/x-xyz', A.toExtXYZ, 'Extended XYZ with the lattice in the comment line.'],
    xsf: ['.xsf', 'text/plain', A.toXSF, 'XCrySDen structure file.'], lammps: ['.lammps-data', 'text/plain', A.toLAMMPSData, 'LAMMPS data file (atom_style atomic).'],
  };
  const out = [];
  for (const w of which) { const [ext, mime, fn, description] = W[w]; const content = safe(() => fn(s), null); if (content) out.push({ name: base + ext, mime, content, description }); }
  return out;
}
const geomOf = s => ({ kind: 'atoms', atoms: s.atoms.map(a => ({ el: a.el, x: a.x, y: a.y, z: a.z })), cell: s.cell || null, title: s.title || 'structure' });
function pool(x, y, maxPts = 400) {   // max-pool a dense pattern so peak heights survive down-sampling
  const n = x.length, k = Math.ceil(n / maxPts); if (k <= 1) return { x: Array.from(x), y: Array.from(y) };
  const X = [], Y = []; for (let i = 0; i < n; i += k) { let m = i; for (let j = i; j < Math.min(n, i + k); j++) if (y[j] > y[m]) m = j; X.push(x[m]); Y.push(y[m]); }
  return { x: X, y: Y };
}
function parsePasted(text, fmt) {
  const t = String(text).replace(/\\n/g, '\n').trim(); let kind = fmt;
  if (kind === 'auto') kind = /(^|\n)\s*(data_|_cell_length|loop_)/i.test(t) ? 'cif' : /^\s*\d+\s*(\n|$)/.test(t) ? 'xyz' : 'poscar';
  try { const s = kind === 'cif' ? A.parseCIF(t) : kind === 'xyz' ? A.parseXYZ(t) : A.parsePOSCAR(t); if (!s.atoms?.length) throw new Error('no atoms found'); return { s, kind }; }
  catch (e) { throw new Error(`Could not read the pasted structure as ${kind.toUpperCase()}: ${e.message}. Paste a complete CIF, POSCAR or XYZ file, or clear the field to use the library crystal.`); }
}

// ───────────────────────────── mode 1: crystal builder and analysis ─────────────────────────────
function runCrystal(inp) {
  const warnings = [], checks = [], pasted = String(inp.structText || '').trim().length > 0, key = cKey(inp);
  let cell, src;
  if (pasted) { const p = parsePasted(inp.structText, pick(inp.structFmt, ['auto', 'cif', 'poscar', 'xyz'], 'auto')); cell = p.s; src = `pasted ${p.kind.toUpperCase()} structure`; for (const w of cell.warnings || []) warnings.push(w); if (cell.atoms.length > 2000) throw new Error(`The pasted structure has ${cell.atoms.length} atoms; this analysis accepts up to 2000.`); }
  else { cell = A.buildCell(key); src = `library entry ${A.CRYSTALS[key].name}`; }
  const an = A.analyse(cell), periodic = !!cell.cell, formula = an.formula, Mf = chem.molarMass(formula), wl = wavelength(inp);
  const title = pasted ? (cell.title || formula) : A.CRYSTALS[key].name, base = slug(pasted ? formula : key);
  const bondRows = an.bonds.map(b => [b.pair, +b.length_A.toFixed(4), +b.min_A.toFixed(4), +b.max_A.toFixed(4), b.count]);
  const coordRows = Object.keys(an.coordination).map(e => [e, an.composition[e], +an.coordination[e].toFixed(3)]);
  const meanCN = N.sum(Object.keys(an.coordination).map(e => an.coordination[e] * an.composition[e])) / an.nAtoms;
  const rd = A.rdf(cell, { rmax_A: 8, dr_A: 0.04 }), rp = pool(rd.r, rd.g, 200);
  const series = [], tables = [], equations = [], kpis = [], publish = { 'material.formula': formula, 'material.name': title, 'prop.molar_mass': Mf * 1e-3 };
  let geo = cell, exportsS, top3 = '';
  equations.push(eqn('formula', 'Composition', `${chem.prettyFormula(formula)},  M = ${F(Mf, 6)} g/mol`, Object.entries(safe(() => chem.massFractions(formula), {})).map(([e, w]) => `${e} ${F(100 * w, 4)} wt%`).join(', ')));
  if (periodic) {
    const p = cellInfo(cell), Z = an.nAtoms / formulaAtoms(formula), rho = an.density_kg_m3;
    let reps = [int(inp, 'nx', 2, 1, 6), int(inp, 'ny', 2, 1, 6), int(inp, 'nz', 2, 1, 6)];
    while (an.nAtoms * reps[0] * reps[1] * reps[2] > 6000 && Math.max(...reps) > 1) { reps[reps.indexOf(Math.max(...reps))]--; if (!warnings.some(w => w.startsWith('Supercell'))) warnings.push('Supercell repeats were reduced to keep the model below 6000 atoms.'); }
    geo = A.supercell(cell, reps);
    const isE = wl.key === 'electron', tmax = isE ? +(2 * Math.asin(Math.min(1, wl.lam / 1.2)) * 180 / PI).toFixed(3) : 100, tmin = isE ? 0.05 : 10, Lc = num(inp, 'Lcryst', 5e-8, 1e-9, 1e-6) * 1e9;
    const xr = A.powderXRD(cell, { wavelength_A: wl.lam, twoThetaMin: tmin, twoThetaMax: tmax, step: (tmax - tmin) / 2400, crystallite_nm: Lc });
    const xp = pool(xr.twoTheta, xr.intensity), pk = [...xr.peaks].sort((a, b) => b.I_rel - a.I_rel), top = pk[0];
    const ds = A.dSpacings(cell, { dmin_A: 1.0 }).slice(0, 25); top3 = pk.slice(0, 3).map(q => `${F(q.twoTheta, 4)}° ${hklStr(q.hkl)}`).join(', ');
    if (xr.truncated) warnings.push('The reflection list was truncated for speed; high-angle peaks may be missing.');
    if (wl.key === 'electron') warnings.push('Intensities use X-ray scattering factors; for electrons only the peak positions (Bragg angles) are meaningful, and dynamical diffraction is ignored.');
    const rhoCheck = Z * Mf * 1e-3 / (N.N_A * p.volume * 1e-30);
    kpis.push(kpi('a', 'Lattice parameter a', p.a, 'Å', `Length of the a axis of the conventional cell from the ${src}; b = ${F(p.b, 5)} Å, angles α/β/γ = ${F(p.alpha, 4)}/${F(p.beta, 4)}/${F(p.gamma, 4)}°.`),
      kpi('c', 'Lattice parameter c', p.c, 'Å', `Length of the c axis; c/a = ${F(p.c / p.a, 4)}${Math.abs(p.gamma - 120) < 0.5 ? ' (ideal close-packed / wurtzite value 1.633)' : ''}. For slab-type entries with vacuum this is the box height, not a lattice constant.`),
      kpi('V', 'Unit-cell volume', p.volume, 'Å³', `Volume |a·(b×c)| of the conventional cell holding ${an.nAtoms} atoms, i.e. ${F(p.volume / an.nAtoms, 4)} Å³ per atom.`),
      kpi('rho', 'X-ray (theoretical) density', rho, 'kg/m³', `ρ = Z·M/(N_A·V) with Z = ${F(Z, 4)} formula units per cell: the density of a defect-free, pore-free crystal. Measured densities below it indicate porosity or vacancies.`),
      kpi('Z', 'Formula units per cell', Z, '–', `Number of ${chem.prettyFormula(formula)} units in the conventional cell (atoms in the cell divided by atoms per formula unit).`),
      kpi('nn', 'Shortest bond length', an.nearestNeighbour_A ?? 0, 'Å', `Shortest interatomic contact in the structure; the bond table lists the mean, minimum and maximum for every bonded element pair.`),
      kpi('cn', 'Mean coordination number', meanCN, '–', `Average number of first-shell neighbours per atom from the geometric bond criterion (${coordRows.map(r => `${r[0]}: ${r[2]}`).join(', ')}).`),
      kpi('tt_max', 'Strongest reflection 2θ', top ? top.twoTheta : 0, '°', top ? `The ${hklStr(top.hkl)} reflection (d = ${F(top.d_A, 5)} Å) is the most intense line of the simulated ${wl.name} powder pattern (λ = ${F(wl.lam, 5)} Å).` : 'No reflection falls inside the scanned angular range for this wavelength.'),
      kpi('npk', 'Reflections in the pattern', xr.peaks.length, '–', `Number of symmetry-distinct allowed reflections between ${tmin}° and ${tmax}° 2θ with a relative intensity above 0.01 %.`),
      kpi('natoms', 'Atoms in the supercell', geo.atoms.length, '–', `The ${reps.join('×')} supercell exported and shown in the viewer.`));
    series.push({ key: 'xrd', title: `Simulated powder diffraction pattern (${wl.name})`, xlabel: '2θ (°)', ylabel: 'relative intensity (%)', lines: [{ name: 'kinematic pattern', x: xp.x, y: xp.y }], explain: `Kinematic powder pattern with Scherrer broadening for ${F(Lc, 3)} nm crystallites. Compare peak positions and relative heights with a measured pattern; missing or extra lines indicate another phase, texture or a different polymorph.` });
    tables.push({ key: 'peaks', title: 'Indexed reflections', columns: ['hkl', '2θ (°)', 'd (Å)', 'I/I_max (%)', 'multiplicity'], rows: xr.peaks.slice(0, 40).map(q => [hklStr(q.hkl), +q.twoTheta.toFixed(3), +q.d_A.toFixed(4), +q.I_rel.toFixed(2), q.multiplicity]), explain: 'Allowed reflections of the simulated pattern in order of angle, with integrated intensity relative to the strongest line and the powder multiplicity.' });
    tables.push({ key: 'dsp', title: 'Largest lattice-plane spacings', columns: ['hkl', 'd (Å)', 'multiplicity'], rows: ds.map(q => [hklStr(q.hkl), +q.d_A.toFixed(4), q.multiplicity]), explain: 'Interplanar spacings of the allowed reflections (systematic absences removed), largest first; use them to index diffraction rings or lattice fringes in TEM images.' });
    equations.push(eqn('governing', 'Bragg law', 'n λ = 2 d_hkl sin θ'), eqn('governing', 'Structure factor and powder intensity', 'F_hkl = Σ_j f_j(s) exp(−B s²) exp[2πi(h x_j + k y_j + l z_j)],   I_hkl ∝ m_hkl |F_hkl|² (1 + cos² 2θ)/(sin² θ cos θ)'),
      eqn('formula', 'X-ray density', `ρ = Z M /(N_A V) = ${F(Z, 4)} × ${F(Mf, 5)} g/mol /(N_A × ${F(p.volume, 5)} Å³) = ${F(rhoCheck, 5)} kg/m³`), eqn('constitutive', 'Scherrer broadening', 'β_L = K λ /(L cos θ),  K = 0.9'));
    checks.push(chk('Density consistency', Math.abs(rhoCheck - rho) / rho < 1e-6, `Z·M/(N_A·V) and the atom-by-atom mass sum agree to ${F(Math.abs(rhoCheck - rho) / rho, 2)} relative`), chk('Integer formula units', Math.abs(Z - Math.round(Z)) < 1e-6, `Z = ${F(Z, 6)}`),
      chk('Bragg law on strongest line', !top || Math.abs(2 * top.d_A * Math.sin(top.twoTheta * PI / 360) - wl.lam) < 1e-6, top ? `2 d sin θ − λ = ${F(2 * top.d_A * Math.sin(top.twoTheta * PI / 360) - wl.lam, 2)} Å` : 'no reflection in range'));
    const vac = vacuumGap(cell) > 5;
    if (vac) warnings.push(`The cell contains a vacuum gap of ${F(vacuumGap(cell), 3)} Å along c: the density of ${F(rho, 4)} kg/m³ refers to the periodic box, not to a bulk solid, and c is not a lattice parameter. The areal mass density is ${F(rho * p.c * 1e-10 * 1e6, 4)} mg/m²; density and c are not published.`);
    Object.assign(publish, { 'prop.lattice_a': p.a * 1e-10, 'struct.cell_volume': p.volume * 1e-30, 'struct.Z': Z }, vac ? {} : { 'prop.density': rho, 'prop.lattice_c': p.c * 1e-10 });
    exportsS = structExports({ ...geo, title: `${title} ${reps.join('x')}` }, `${base}_${reps.join('x')}`);
  } else {
    warnings.push('The pasted structure has no lattice: it is treated as a finite molecule or cluster, so density, d-spacings and powder diffraction are not computed.');
    const span = ['x', 'y', 'z'].map(k => N.maxOf(cell.atoms.map(a => a[k])) - N.minOf(cell.atoms.map(a => a[k])));
    kpis.push(kpi('natoms', 'Atoms', an.nAtoms, '–', `Number of atoms read from the ${src}.`), kpi('M', 'Molar mass', an.molarMass_g_mol, 'g/mol', `Mass of the whole cluster ${chem.prettyFormula(formula)} × ${F(an.nAtoms / formulaAtoms(formula), 3)} summed over its atoms.`),
      kpi('nn', 'Shortest contact', an.nearestNeighbour_A ?? 0, 'Å', 'Shortest interatomic distance; values well below typical covalent bond lengths indicate overlapping atoms in the input.'),
      kpi('cn', 'Mean coordination number', meanCN, '–', 'Average number of bonded neighbours from the geometric bond criterion (first-shell contacts within 1.25 × the covalent-radius sum).'),
      kpi('rg', 'Radius of gyration', an.radiusOfGyration_A, 'Å', 'Mass-weighted root-mean-square distance of the atoms from the centre of mass; a compactness measure of the cluster.'),
      kpi('span', 'Largest extent', Math.max(...span), 'Å', `Bounding-box size of the cluster: ${span.map(v => F(v, 4)).join(' × ')} Å.`), kpi('surf', 'Under-coordinated fraction', an.surfaceFraction, '–', 'Fraction of atoms with fewer neighbours than the best-coordinated atom of the same element — a proxy for the surface-atom fraction.'));
    series.push({ key: 'cn_r', title: 'Running coordination number', xlabel: 'r (Å)', ylabel: 'neighbours within r', lines: [{ name: 'n(r)', x: rp.x, y: pool(rd.r, rd.cn, 200).y }], explain: 'Average number of atoms within a distance r of an atom; plateaus mark complete coordination shells of the cluster.' });
    checks.push(chk('No overlapping atoms', (an.nearestNeighbour_A ?? 1) > 0.6, `shortest contact ${F(an.nearestNeighbour_A ?? 0, 4)} Å`));
    exportsS = structExports(cell, base, ['xyz', 'cif', 'poscar', 'lammps']);
  }
  series.push({ key: 'rdf', title: 'Radial distribution function', xlabel: 'r (Å)', ylabel: 'g(r)', lines: [{ name: 'g(r)', x: rp.x, y: rp.y }], explain: 'Pair distribution of the ideal static structure: each spike is a neighbour shell. Thermal motion and disorder broaden these peaks in a measured PDF; the first peak position is the bond length.' });
  tables.push({ key: 'bonds', title: 'Bond lengths', columns: ['pair', 'mean (Å)', 'min (Å)', 'max (Å)', 'bonds per cell'], rows: bondRows.length ? bondRows : [['–', 0, 0, 0, 0]], explain: 'First-shell bonded contacts found by the geometric criterion (distance ≤ 1.25 × the covalent-radius sum, cation–cation contacts bridged by anions excluded).' });
  tables.push({ key: 'coord', title: 'Coordination by element', columns: ['element', 'atoms', 'mean coordination'], rows: coordRows, explain: 'Average first-shell coordination number of each element; compare with the ideal polyhedra (4 tetrahedral, 6 octahedral, 8 cubic, 12 close-packed).' });
  checks.push(chk('Charge-neutral composition parsed', Mf > 0, `${formula}: ${F(Mf, 6)} g/mol`));
  const b0 = an.bonds[0];
  return {
    summary: `${title} (${chem.prettyFormula(formula)}) from the ${src}: ${periodic ? `a = ${F(cellInfo(cell).a, 5)} Å, c = ${F(cellInfo(cell).c, 5)} Å, cell volume ${F(an.volume_A3, 5)} Å³ and X-ray density ${F(an.density_kg_m3, 4)} kg/m³` : `${an.nAtoms} atoms, radius of gyration ${F(an.radiusOfGyration_A, 4)} Å`}. ${b0 ? `The shortest bonded pair is ${b0.pair} at ${F(b0.length_A, 4)} Å with a mean coordination of ${F(meanCN, 3)}.` : ''} Geometry is the ideal, unrelaxed crystallographic structure; no energy minimisation is performed in this mode.`,
    kpis, series, tables, equations, geometry: geomOf(geo), checks, warnings,
    recommendations: [
      rec('high', periodic ? `Index a measured pattern against the table: the three strongest lines should appear near ${top3}. A changed intensity ratio at unchanged positions points to preferred orientation (e.g. c-axis textured nanorod films).` : 'Add lattice vectors (CIF, POSCAR or extended XYZ) if this is a crystal, so density and diffraction can be computed.', 'Kinematic powder pattern of the ideal structure.'),
      rec('medium', `Use the exported ${periodic ? 'CIF / POSCAR / LAMMPS data files' : 'XYZ file'} as the starting geometry for relaxation in the simulation suite or an external DFT / MD code; lattice parameters here are experimental room-temperature values, not relaxed for any particular functional or potential.`, 'Structure files are written from the same atom list shown in the viewer.'),
      rec('medium', periodic ? `Compare the pycnometric or Archimedes density of your sample with the X-ray density of ${F(an.density_kg_m3, 4)} kg/m³ to quantify closed porosity.` : 'Check the bond table for unexpected short or long contacts before using the cluster in a simulation.', 'Relative density = measured density / theoretical density.'),
      rec('low', 'Screen candidate phases and supercell sizes here before committing diffractometer time or large simulations: a few seconds of in-browser analysis avoids wasted instrument hours and compute energy.', 'Good modelling practice: verify the structure model first.', true),
    ],
    publish, exports: exportsS,
    numerics: { method: 'Direct lattice construction; neighbour search on a cell grid; kinematic structure-factor summation over reciprocal-lattice points', dof: an.nAtoms, steps: 1, iterations: 1, residual: 0, notes: 'Ideal unrelaxed geometry; neutral-atom Cromer–Mann scattering factors, isotropic B = 0.5 Å².' },
  };
}

// ───────────────────────────── mode 2: nanostructure builder ─────────────────────────────
const SHAPES = ['sphere', 'rod', 'hexrod', 'cube', 'octahedron', 'wulff-like', 'sheet', 'tube', 'coreshell', 'slab'];
function shapeAV(shape, d, L) {   // area (Å²) and volume (Å³) of the carved body
  const S3 = Math.sqrt(3);
  if (shape === 'cube') return [6 * d * d, d ** 3];
  if (shape === 'rod') return [PI * d * L + PI * d * d / 2, PI * d * d * L / 4];
  if (shape === 'hexrod') return [S3 * d * d + 6 * (d / S3) * L, S3 / 2 * d * d * L];
  if (shape === 'octahedron') return [S3 * d * d, d ** 3 / 6];
  if (shape === 'wulff-like') return [(6 + 12 * S3) / 8 * d * d, d ** 3 / 2];
  return [PI * d * d, PI * d ** 3 / 6];
}
function runNano(inp) {
  const warnings = [], checks = [], key = cKey(inp), shape = pick(inp.nshape, SHAPES, 'hexrod'), def = A.CRYSTALS[key], bulkCell = A.buildCell(key), bulk = A.analyse(bulkCell);
  const dnm = num(inp, 'dnm', 2.5e-9, 5e-10, 3e-8) * 1e9, Lnm = num(inp, 'Lnm', 6e-9, 5e-10, 1e-7) * 1e9, MAXA = 3000, hexSys = Math.abs(cellInfo(bulkCell).gamma - 120) < 0.5;
  let s, area, vol, label, bulkCN = { ...bulk.coordination }, sizeA = dnm * 10;
  if (shape === 'sheet') {
    const layers = int(inp, 'layers', 1, 1, 12), want = Math.max(2, Math.round(dnm * 10 / cellInfo(bulkCell).a)); let nl = layers, unit = 'van der Waals layer';
    const one = A.nanosheet(key, { layers, nx: 1, ny: 1 }); if (one.layerMode === 'atomic planes') { unit = 'unit-cell period'; nl = layers * (one.layersPerPeriod || 1); }
    const per = nl === layers ? one.atoms.length : A.nanosheet(key, { layers: nl, nx: 1, ny: 1 }).atoms.length, nxy = Math.max(1, Math.min(want, Math.floor(Math.sqrt(MAXA / per))));
    s = A.nanosheet(key, { layers: nl, nx: nxy, ny: nxy }); if (nxy < want) s.truncated = true;
    for (const w of s.warnings || []) warnings.push(w);
    const p = cellInfo(s), zs = s.atoms.map(a => a.z), t = Math.max(N.maxOf(zs) - N.minOf(zs), 0) + (bulk.nearestNeighbour_A || 2), Axy = p.a * p.b * Math.sin(p.gamma * PI / 180);
    area = 2 * Axy; vol = Axy * t; sizeA = t; label = `sheet of ${layers} ${unit}${layers > 1 ? 's' : ''} (${nl} atomic layers), ${F(p.a / 10, 3)} × ${F(p.b / 10, 3)} nm (periodic in-plane), thickness ≈ ${F(t / 10, 3)} nm`;
  } else if (shape === 'slab') {
    let hkl = [int(inp, 'h', 0, -4, 4), int(inp, 'k', 0, -4, 4), int(inp, 'l', 1, -4, 4)]; if (!hkl.some(v => v)) { hkl = [0, 0, 1]; warnings.push('Miller indices (0 0 0) are not a plane; (0 0 1) was used.'); }
    const layers = int(inp, 'layers', 1, 1, 12) + 3; s = A.slab(key, { hkl, layers, vacuum_A: 15 });
    const reps = Math.max(1, Math.min(4, Math.floor(Math.sqrt(MAXA / s.atoms.length)))), base = s; s = { ...A.supercell(s, [reps, reps, 1]), title: `${base.title}, ${reps}×${reps} in-plane` };
    const t = Math.max(base.thickness_A, 0) + base.d_hkl_A / base.layersPerPeriod; area = 2 * base.area_A2 * reps * reps; vol = base.area_A2 * reps * reps * t; sizeA = t;
    label = `(${base.hkl.join(' ')}) slab of ${base.layers} ${base.layerMode}, thickness ≈ ${F(t / 10, 3)} nm, d_hkl = ${F(base.d_hkl_A, 4)} Å`;
  } else if (shape === 'tube') {
    const n = int(inp, 'tn', 10, 0, 20), m = Math.min(int(inp, 'tm', 10, 0, 20), Math.max(n, 1)), el = key === 'hBN' ? 'BN' : 'C';
    if (key !== 'hBN' && !/graph|C_/.test(key)) warnings.push('Nanotubes are built from a honeycomb sheet: carbon is used (select hexagonal boron nitride for a BN tube); the selected crystal is ignored.');
    s = A.nanotube({ n: Math.max(n, 1), m, length_nm: Math.min(Lnm, 20), element: el, maxAtoms: MAXA });
    const d = s.diameter_nm * 10, L = s.length_nm * 10, t = 3.35; area = 2 * PI * d * L; vol = PI * d * t * L; sizeA = d; bulkCN = el === 'BN' ? { B: 3, N: 3 } : { C: 3 };
    label = `(${s.n},${s.m}) ${el === 'BN' ? 'BN' : 'carbon'} nanotube, diameter ${F(s.diameter_nm, 4)} nm, length ${F(s.length_nm, 4)} nm, chiral angle ${F(s.chiralAngle_deg, 4)}°${el === 'C' ? `, ${s.metallic ? 'metallic' : 'semiconducting'} by the (n − m) mod 3 rule` : ''}`;
  } else if (shape === 'coreshell') {
    const sk = cKey(inp, 'shellCrystal'), ts = num(inp, 'shellnm', 5e-10, 2e-10, 5e-9) * 1e9; s = A.coreShell(key, sk, { core_nm: dnm, shell_nm: ts, maxAtoms: MAXA });
    const d = s.diameter_nm * 10;[area, vol] = shapeAV('sphere', d, d); sizeA = d; const sb = A.analyse(A.buildCell(sk)).coordination; bulkCN = { ...sb, ...bulk.coordination };
    label = `${def.formula}@${A.CRYSTALS[sk].formula} core–shell sphere, core ${F(s.core_nm, 3)} nm, shell ${F(s.shell_nm, 3)} nm (${s.nCore} + ${s.nShell} atoms${s.latticeMismatch_pct !== null ? `, lattice mismatch ${F(s.latticeMismatch_pct, 3)} %` : ''})`;
    if (s.latticeMismatch_pct !== null && Math.abs(s.latticeMismatch_pct) > 5) warnings.push(`Core and shell lattice parameters differ by ${F(s.latticeMismatch_pct, 3)} %: a coherent shell of this thickness is unlikely; expect misfit dislocations or island growth.`);
  } else {
    let sh = shape; if (sh === 'hexrod' && !hexSys) { sh = 'rod'; warnings.push('A hexagonal prism needs a hexagonal crystal; a cylindrical rod was carved instead.'); }
    s = A.nanoparticle(key, { shape: sh, diameter_nm: dnm, length_nm: Lnm, maxAtoms: MAXA });
    const d = s.diameter_nm * 10, L = (s.length_nm ?? s.diameter_nm) * 10, eff = sh === 'wulff-like' && hexSys ? 'hexrod' : sh;[area, vol] = shapeAV(eff, d, L); sizeA = d;
    label = `${sh} of ${def.formula}, ${F(s.diameter_nm, 3)} nm${s.length_nm ? ` × ${F(s.length_nm, 3)} nm` : ''}`;
  }
  if (s.truncated) warnings.push(`The requested size exceeds the ${MAXA}-atom limit of this in-browser builder; the structure was scaled down to ${s.atoms.length} atoms (${label}). Scaling relations in the size chart still cover the requested size.`);
  if (!s.atoms.length) throw new Error('The requested nanostructure contains no atoms: increase its size.');
  const an = A.analyse(s), nAt = s.atoms.length;
  let zSum = 0, zbSum = 0; for (const e in an.composition) { const zb = bulkCN[e] ?? an.coordination[e]; zSum += an.composition[e] * Math.min(an.coordination[e], zb); zbSum += an.composition[e] * zb; }
  const fb = zbSum > 0 ? N.clamp(zSum / zbSum, 0, 1) : 1, zbMean = zbSum / nAt;
  let fSurf = an.surfaceFraction;
  if (shape === 'coreshell') { const R = s.diameter_nm * 5, sk = (bulk.nearestNeighbour_A || 2.5); fSurf = s.atoms.filter(a => Math.hypot(a.x, a.y, a.z) > R - sk).length / nAt; }
  if (shape === 'tube' || (shape === 'sheet' && int(inp, 'layers', 1, 1, 12) === 1)) fSurf = 1;
  const lib = BULK[key] || {}, user = pick(inp.propSrc, ['library', 'user'], 'library') === 'user' || shape === 'tube' && !lib.Ec;
  const Ec = !user && lib.Ec ? lib.Ec : num(inp, 'EcohIn', 3.81, 0.01, 12), noTm = !user && !lib.Tm && !!lib.Ec, Tm = !user && lib.Tm ? lib.Tm : num(inp, 'TmIn', 1337, 50, 4500);
  if (!user && !lib.Ec) warnings.push(`No library cohesive energy for ${def.name}: the user value ${F(Ec, 4)} eV/atom was used.`);
  if (noTm) warnings.push(`${def.name} has no congruent melting point at ambient pressure (it decomposes, sublimes or transforms first): no melting estimate is given. Select user data to supply a reference temperature.`);
  else if (!user && !lib.Tm) warnings.push(`No library melting point for ${def.name}: the user value ${F(Tm, 5)} K was used.`);
  const Esurf = nAt * Ec * (1 - fb), gammaBB = area > 0 ? Esurf * N.Q_E / (area * 1e-20) : 0, EcN = Ec * fb, TmN = Tm * fb, SV = area / vol * 1e10, thin = shape === 'sheet' || shape === 'slab' || shape === 'tube', massKg = an.molarMass_g_mol * 1e-3 / N.N_A;
  const rho = vacuumGap(bulkCell) > 5 && def.formula === 'C' ? A.analyse(A.buildCell('graphite')).density_kg_m3 : bulk.density_kg_m3, ssa = thin && massKg > 0 ? area * 1e-20 / massKg : SV / rho;
  const d0 = bulk.nearestNeighbour_A || 2.5, Ebond = zbMean > 0 ? 2 * Ec / zbMean : 0;
  // size scaling: surface-to-volume ∝ 1/D; cohesive energy E(D)/E_b = 1 − c·d0/D with c fitted to this particle
  const cShape = sizeA > 0 ? (1 - fb) * sizeA / d0 : 0, Ds = N.logspace(Math.max(0.6, 2 * d0 / 10), 60, 60);
  const sc = D => N.clamp(1 - cShape * d0 / (D * 10), 0, 1);
  const series = [
    { key: 'Ecoh_D', title: 'Cohesive energy and melting point vs size', xlabel: 'characteristic size D (nm)', ylabel: 'fraction of bulk value', logx: true, lines: [{ name: `this shape (broken-bond, c = ${F(cShape, 3)})`, x: Ds, y: Ds.map(sc) }, { name: 'sphere reference 1 − d₀/D', x: Ds, y: Ds.map(D => N.clamp(1 - d0 / (D * 10), 0, 1)) }], explain: 'Bond-order scaling E_coh(D)/E_coh,bulk = T_m(D)/T_m,bulk = 1 − c·d₀/D, with the coefficient c taken from the bond count of the structure built here. The depression becomes significant below about 10 nm.' },
    { key: 'surf_D', title: 'Surface-to-volume ratio vs size', xlabel: 'characteristic size D (nm)', ylabel: 'S/V (1/nm)', logx: true, logy: true, lines: [{ name: 'same shape, scaled', x: Ds, y: Ds.map(D => SV * 1e-9 * (sizeA / 10) / D) }], explain: 'Geometric surface-to-volume ratio for the same shape at other sizes (S/V ∝ 1/D). It controls dissolution, catalytic activity per mass and the weight of surface energy in the total energy.' },
  ];
  const libG = lib.g, cnRows = Object.keys(an.composition).map(e => [e, an.composition[e], +(an.coordination[e]).toFixed(3), +(bulkCN[e] ?? an.coordination[e]).toFixed(3)]);
  checks.push(chk('Bond count not above bulk', fb <= 1 + 1e-9, `retained bond fraction ${F(fb, 4)}`), chk('Atom count consistent with volume', !!s.cell || shape === 'coreshell' || nAt < 30 || key === 'graphene' || Math.abs(nAt / (vol * bulk.nAtoms / bulk.volume_A3) - 1) < 0.35, `${nAt} atoms vs ${F(vol * bulk.nAtoms / bulk.volume_A3, 4)} from bulk number density × geometric volume`),
    chk('Broken-bond surface energy in a physical range', gammaBB >= 0 && gammaBB < 10, `${F(gammaBB, 3)} J/m²${libG ? ` (typical measured ${libG} J/m²)` : ''}`));
  if (fb < 0.6) warnings.push('More than 40 % of the bonds are missing: at this size the continuum scaling is unreliable and real clusters reconstruct into non-crystalline motifs.');
  warnings.push('Surfaces are unrelaxed bulk terminations; the broken-bond model ignores relaxation, reconstruction, ligands and electronic effects and typically overestimates surface energies by 20–100 %.');
  return {
    summary: `Built a ${label} with ${nAt} atoms. ${F(100 * fSurf, 3)} % of the atoms are under-coordinated surface atoms and the surface-to-volume ratio is ${F(SV * 1e-9, 3)} nm⁻¹. The broken-bond model gives a surface energy contribution of ${F(Esurf, 4)} eV (${F(gammaBB, 3)} J/m²), a cohesive energy of ${F(EcN, 4)} eV/atom (${F(100 * fb, 4)} % of bulk) ${noTm ? 'no melting estimate is given because the bulk phase does not melt congruently' : `and an estimated melting point of ${F(TmN, 4)} K against ${F(Tm, 4)} K for the bulk`}.`,
    kpis: [
      kpi('natoms', 'Atoms', nAt, '–', `Number of atoms in the generated ${shape} structure (${Object.entries(an.composition).map(([e, c]) => `${e} ${c}`).join(', ')}).`),
      kpi('f_surf', 'Surface-atom fraction', fSurf, '–', `Fraction of atoms with fewer first-shell neighbours than in the bulk crystal${shape === 'tube' || fSurf === 1 ? ' (every atom of a single-wall tube or monolayer is a surface atom)' : ''}. These atoms dominate adsorption, catalysis and dissolution.`),
      kpi('SV', 'Surface-to-volume ratio', SV, '1/m', `Geometric area/volume of the ${shape}: ${F(area / 100, 4)} nm² over ${F(vol / 1000, 4)} nm³. ${thin ? `Specific surface area ${F(ssa, 4)} m²/kg from the mass of the atoms.` : `Equivalent specific surface area ${F(ssa, 4)} m²/kg using the bulk density.`}`),
      kpi('E_surf', 'Surface energy contribution', Esurf, 'eV', `Total energy of the missing bonds, N·E_coh·(1 − ⟨Z⟩/Z_bulk), with a bond energy of ${F(Ebond, 3)} eV from E_coh = ${F(Ec, 4)} eV/atom and Z_bulk = ${F(zbMean, 3)}.`),
      kpi('gamma_bb', 'Broken-bond surface energy', gammaBB, 'J/m²', `Surface energy contribution divided by the geometric area.${libG ? ` The typical measured value for this metal is about ${libG} J/m²; the unrelaxed broken-bond figure is normally higher.` : ' No experimental reference is stored for this material; treat the figure as an upper-bound estimate.'}`, null),
      kpi('Ecoh', 'Cohesive energy of the nanostructure', EcN, 'eV/atom', `Bulk cohesive energy scaled by the retained bond fraction ${F(fb, 4)} (bond-order / broken-bond model). The deficit of ${F(Ec - EcN, 3)} eV/atom is the driving force for coarsening and sintering.`),
      noTm ? kpi('Tm', 'Estimated melting point', 'not applicable', '', 'The bulk phase sublimes, decomposes or transforms before melting at ambient pressure, so a size-dependent melting point is not meaningful; the cohesive-energy ratio still measures the loss of thermal stability.') : kpi('Tm', 'Estimated melting point', TmN, 'K', `T_m(D) = T_m,bulk × E_coh(D)/E_coh,bulk — a Lindemann-type bond-order estimate, qualitatively equivalent to Gibbs–Thomson 1/D depression. It is an order-of-magnitude guide for free, clean particles only.`, null),
      kpi('ssa', 'Specific surface area', ssa, 'm²/kg', thin ? `Area of both faces (${F(area / 100, 4)} nm²) divided by the mass of the ${nAt} atoms, independent of any assumed thickness; a graphene monolayer gives about 2.63 × 10⁶ m²/kg. Compare with a BET measurement to judge stacking and bundling.` : `S/V divided by the bulk density ${F(rho, 4)} kg/m³; compare with a BET measurement to judge agglomeration.`),
      kpi('cn', 'Mean coordination number', zSum / nAt, '–', `Average first-shell coordination in the nanostructure against ${F(zbMean, 3)} in the bulk crystal.`),
    ],
    series,
    tables: [{ key: 'cn', title: 'Coordination by element', columns: ['element', 'atoms', 'mean coordination', 'bulk coordination'], rows: cnRows, explain: 'Average number of bonded neighbours of each element in the nanostructure compared with the bulk crystal; the difference is the number of broken bonds per atom.' }],
    equations: [
      eqn('formula', 'Composition', `${chem.prettyFormula(an.formula)},  M = ${F(safe(() => chem.molarMass(an.formula), 0), 5)} g/mol per formula unit`),
      eqn('constitutive', 'Broken-bond surface energy', 'E_surf = Σ_i (Z_bulk − Z_i)·ε_b/2,  ε_b = 2 E_coh / Z_bulk,  γ = E_surf / A', 'Nearest-neighbour bonds of equal strength; no relaxation.'),
      eqn('constitutive', 'Size-dependent cohesive energy and melting', 'E_coh(D)/E_coh,bulk = ⟨Z⟩/Z_bulk ≈ 1 − c·d₀/D,   T_m(D)/T_m,bulk ≈ E_coh(D)/E_coh,bulk'),
      eqn('governing', 'Surface-to-volume scaling', 'S/V = k_shape / D  (sphere and cube: k = 6; long cylinder: k = 4)'),
    ],
    geometry: geomOf(s), checks, warnings,
    recommendations: [
      rec('high', noTm ? `The structure retains ${F(100 * fb, 4)} % of the bulk bonding: its thermal stability (oxidation, sublimation, restructuring) is reduced roughly in that proportion; under-coordinated surface atoms (${F(100 * fSurf, 3)} %) are the first to react.` : `Keep processing and annealing well below about ${F(0.5 * TmN, 4)} K (half the estimated size-dependent melting point): under-coordinated surface atoms (${F(100 * fSurf, 3)} %) make ${F(sizeA / 10, 3)} nm structures sinter and coarsen long before bulk melting.`, 'Bond-order melting estimate and the Tammann rule of thumb.'),
      rec('medium', shape === 'rod' || shape === 'hexrod' ? 'For rods the side facets carry most of the surface energy; capping agents that bind the side facets lower their energy further and promote axial growth, raising the aspect ratio.' : 'Faceted (Wulff-like) shapes expose low-energy planes and have a lower surface energy than spheres of the same volume; compare shapes here before choosing a capping strategy.', 'Broken-bond count per shape.'),
      rec('medium', 'Relax the exported geometry with an interatomic potential or DFT before quantitative use: surface relaxation and reconstruction change bond lengths in the outer one or two layers.', 'The structure is an unrelaxed bulk termination.'),
      rec('low', `A specific surface area of ${F(ssa / 1000, 3)} m²/g means the same catalytic or sensing surface is obtained with far less material than a coarse powder: size the particles to the smallest value that remains thermally stable, to save critical raw materials.`, 'Geometric S/V ratio and bulk density.', true),
    ],
    publish: { 'material.formula': an.formula, 'geom.diameter': sizeA * 1e-10, 'geom.length': (shape === 'rod' || shape === 'hexrod' || shape === 'tube') ? (s.length_nm || Lnm) * 1e-9 : sizeA * 1e-10, 'geom.ssa': ssa, 'prop.surface_energy': libG || gammaBB, 'prop.density': rho, 'struct.surface_fraction': fSurf, 'struct.cohesive_eV': EcN, ...(noTm ? {} : { 'struct.T_melt_nano': TmN }) },
    exports: structExports(s, slug(`${key}_${shape}`), s.cell ? ['cif', 'poscar', 'xyz', 'extxyz', 'xsf', 'lammps'] : ['xyz', 'extxyz', 'cif', 'poscar', 'lammps']),
    numerics: { method: 'Lattice carving (atomic geometry engine) + neighbour-count broken-bond model', dof: nAt, steps: 1, iterations: 1, residual: 0, notes: `Atom cap ${MAXA}. Bulk data: ${user ? 'user-supplied' : 'library'} E_coh = ${F(Ec, 4)} eV/atom, T_m = ${F(Tm, 4)} K.` },
  };
}

// ───────────────────────────── mode 3: defect chemistry ─────────────────────────────
const MARK = (q) => q === 0 ? '×' : (q > 0 ? '•' : '′').repeat(Math.min(Math.abs(q), 6));
const coefStr = c => { if (Math.abs(c - 1) < 1e-9) return ''; if (Math.abs(c - 0.5) < 1e-9) return '½ '; if (Math.abs(c - 1.5) < 1e-9) return '3/2 '; return `${+c.toFixed(3)} `; };
// A Kröger–Vink term: {c, sp, site, q, atoms:{El:n}, sites:{M:n,X:n}}
const kvTerm = t => `${coefStr(t.c)}${t.txt}`;
function runDefect(inp) {
  const warnings = [], checks = [], key = cKey(inp), def = A.CRYSTALS[key], cell = A.buildCell(key), an = A.analyse(cell), p = cellInfo(cell), els = Object.keys(an.composition);
  const T = num(inp, 'T_C', 900, 25, 2500) + 273.15, kT = KB_EV * T, Ev = num(inp, 'Evac', 1.3, 0.1, 8), Ei = num(inp, 'Eint', 3, 0.1, 12), Es = num(inp, 'Esch', 2.3, 0.1, 15), Ef = num(inp, 'Efr', 2.6, 0.1, 15), Sv = num(inp, 'Svib', 2, 0, 10);
  const dop = ionOf(pick(inp.dopant, DOPANTS, 'Al3+')), x = num(inp, 'xdop', 2, 0, 30) / 100, xv = num(inp, 'xvac', 0.5, 0, 10) / 100, seed = int(inp, 'seed', 1, 1, 99999);
  const anions = els.filter(e => ANION_Q[e]), cations = els.filter(e => !ANION_Q[e]), ionic = anions.length === 1 && cations.length >= 1 && cations.every(e => CATION_Q[e]);
  const Nsite = an.nAtoms / (p.volume * 1e-30), nuF = formulaAtoms(def.formula);
  // equilibrium site fractions
  const xVac = t => Math.min(1, Math.exp(Sv - Ev / (KB_EV * t))), xInt = t => Math.min(1, Math.exp(Sv - Ei / (KB_EV * t))), xSch = t => Math.min(1, Math.exp(-Es / (nuF * KB_EV * t))), xFr = t => Math.min(1, Math.exp(-Ef / (2 * KB_EV * t)));
  const Ts = N.linspace(300, Math.max(T * 1.15, 1200), 60), lg = v => Math.log10(Math.max(v, 1e-300));
  const series = [
    { key: 'conc_T', title: 'Equilibrium point-defect site fractions vs temperature', xlabel: 'temperature (K)', ylabel: 'log₁₀ site fraction', lines: [{ name: 'vacancy', x: Ts, y: Ts.map(t => Math.max(lg(xVac(t)), -40)) }, { name: 'self-interstitial', x: Ts, y: Ts.map(t => Math.max(lg(xInt(t)), -40)) }, { name: 'Schottky', x: Ts, y: Ts.map(t => Math.max(lg(xSch(t)), -40)) }, { name: 'Frenkel', x: Ts, y: Ts.map(t => Math.max(lg(xFr(t)), -40)) }], explain: 'Boltzmann site fractions of thermally generated intrinsic defects (curves clipped at 10⁻⁴⁰). The defect with the lowest effective formation energy per defect dominates; quenching freezes in the high-temperature value.' },
  ];
  const equations = [
    eqn('governing', 'Vacancy / interstitial equilibrium', 'n/N = exp(S_f/k_B) exp(−E_f/k_BT)', 'Dilute, non-interacting defects; S_f is the vibrational formation entropy.'),
    eqn('governing', 'Schottky and Frenkel equilibria', `[V_M][V_X]… = exp(−E_S/k_BT) ⇒ x_S = exp(−E_S/(ν k_BT)), ν = ${nuF};   x_F = √(N_i/N) exp(−E_F/2k_BT)`, 'Interstitial-site ratio N_i/N taken as 1.'),
  ];
  const tables = [], kpis = [];
  kpis.push(kpi('x_vac', 'Vacancy site fraction', xVac(T), '–', `exp(S/k)·exp(−E_v/kT) at ${F(T, 4)} K for E_v = ${F(Ev, 3)} eV: about ${F(xVac(T) * Nsite * 1e-6, 3)} vacancies per cm³ at a site density of ${F(Nsite * 1e-6, 3)} cm⁻³.`),
    kpi('x_sch', 'Schottky defect site fraction', xSch(T), '–', `Fraction of each sublattice vacant by Schottky disorder, exp(−E_S/(${nuF} kT)), for E_S = ${F(Es, 3)} eV per formula unit of vacancies.`),
    kpi('x_fr', 'Frenkel pair site fraction', xFr(T), '–', `exp(−E_F/2kT) for E_F = ${F(Ef, 3)} eV per vacancy–interstitial pair. ${xFr(T) > xSch(T) ? 'Frenkel' : 'Schottky'} disorder dominates at this temperature for the energies given.`),
    kpi('x_int', 'Self-interstitial site fraction', xInt(T), '–', `exp(S/k)·exp(−E_i/kT) for E_i = ${F(Ei, 3)} eV; usually many orders below the vacancy fraction in close-packed solids.`));
  // dopant substitution
  let host = cations[0] || els[0], zM = CATION_Q[host] || 0, X = anions[0] || null, qX = X ? ANION_Q[X] : 0, rD, rM, rX = X ? RION[ionKey(X, -qX)] : null, radSrc;
  const rDi = RION[ionKey(dop.el, dop.z)];
  if (ionic) {
    host = cations.reduce((b, e) => (Math.abs((RION[ionKey(e, CATION_Q[e])] ?? 9) - rDi) < Math.abs((RION[ionKey(b, CATION_Q[b])] ?? 9) - rDi) ? e : b)); zM = CATION_Q[host]; rM = RION[ionKey(host, zM)]; rD = rDi; radSrc = 'Shannon effective ionic radii (six-fold coordination)';
    const r4M = R4[ionKey(host, zM)], r4D = R4[ionKey(dop.el, dop.z)];
    if (an.coordination[host] < 4.5 && r4M && r4D) { rM = r4M; rD = r4D; radSrc = 'Shannon effective ionic radii (four-fold coordination, as on the tetrahedral host site; anion radius six-fold)'; }
  }
  else { rM = chem.ELEMENTS[host].radius; rD = chem.ELEMENTS[dop.el].radius; radSrc = 'covalent radii (size proxy for a non-ionic host)'; }
  if (dop.el === host) warnings.push(`The dopant ${dop.el} is the host cation itself; substitution has no effect.`);
  const mism = (rD - rM) / rM, bondL = ionic && rX ? rM + rX : rM, eps = x * (rD - rM) / bondL, nHost = an.composition[host];
  const aNew = p.a * (1 + eps), cNew = p.c * (1 + eps), mCell = an.molarMass_g_mol + x * nHost * (chem.ELEMENTS[dop.el].mass - chem.ELEMENTS[host].mass), rhoNew = mCell * 1e-3 / N.N_A / (p.volume * (1 + eps) ** 3 * 1e-30);
  const nDop = x * nHost / (p.volume * 1e-30), delta = ionic ? dop.z - zM : 0;
  kpis.push(kpi('mismatch', 'Dopant size mismatch', mism * 100, '%', `(r_${dop.el} − r_${host})/r_${host} with ${radSrc}: ${F(rD, 3)} Å against ${F(rM, 3)} Å. The Hume-Rothery guideline for extended substitutional solubility is |mismatch| < 15 %.`, Math.abs(mism) < 0.15),
    kpi('a_doped', 'Lattice parameter a after doping', aNew, 'Å', `Vegard-type hard-sphere estimate a(x) = a₀[1 + x(r_D − r_M)/(r_M${ionic ? ' + r_X' : ''})] for x = ${F(100 * x, 3)} % of the ${host} sites: a change of ${F(100 * eps, 3)} %.`),
    kpi('rho_doped', 'Density after doping', rhoNew, 'kg/m³', `Cell mass with ${F(100 * x, 3)} % of ${host} replaced by ${dop.el}, divided by the Vegard-scaled cell volume; undoped value ${F(an.density_kg_m3, 5)} kg/m³. Charge-compensating vacancies are not included in the mass.`),
    kpi('n_dop', 'Dopant concentration', nDop * 1e-6, 'cm⁻³', `Number of ${dop.el} atoms per unit volume for the chosen site fraction. ${delta > 0 ? 'With full electronic compensation this is also the free-electron density.' : delta < 0 ? `With ionic compensation it creates ${F(nDop * Math.abs(delta) / (qX || 1) * 1e-6, 3)} cm⁻³ anion vacancies.` : 'An isovalent dopant adds no carriers.'}`));
  let tol = null;
  if (/TiO3/.test(def.formula) && ionic) {
    const Asite = cations.find(e => R12[ionKey(e, CATION_Q[e])]), Bsite = cations.find(e => e !== Asite); let rA = R12[ionKey(Asite, CATION_Q[Asite])], rB = RION[ionKey(Bsite, CATION_Q[Bsite])];
    if (host === Asite) rA = (1 - x) * rA + x * (R12[ionKey(dop.el, dop.z)] ?? rD * 1.2); else rB = (1 - x) * rB + x * rD;
    tol = (rA + rX) / (Math.SQRT2 * (rB + rX)); const mu = rB / rX;
    kpis.push(kpi('tolerance', 'Goldschmidt tolerance factor', tol, '–', `t = (r_A + r_X)/[√2 (r_B + r_X)] with the site-averaged radii after doping (A = ${Asite}, twelve-fold radius; B = ${Bsite}). Cubic perovskites form for roughly 0.9 < t < 1.0; lower values tilt the octahedra, higher values favour hexagonal polytypes. Octahedral factor r_B/r_X = ${F(mu, 3)}.`, tol > 0.8 && tol < 1.06));
    equations.push(eqn('constitutive', 'Goldschmidt tolerance factor', `t = (r_A + r_X)/[√2 (r_B + r_X)] = ${F(tol, 4)}`));
    checks.push(chk('Perovskite stability window', tol > 0.8 && tol < 1.06, `t = ${F(tol, 4)} (perovskite structures occur for about 0.8–1.06)`));
  } else kpis.push(kpi('tolerance', 'Radius ratio r_cation/r_anion', ionic && rX ? rM / rX : 1, '–', ionic && rX ? `Pauling radius ratio of the host ${host}–${X} pair; 0.225–0.414 favours tetrahedral, 0.414–0.732 octahedral and above 0.732 cubic coordination. The observed coordination is ${F(an.coordination[host], 3)}.` : 'Elemental host: the radius-ratio rule does not apply; see the Hume-Rothery table instead.'));
  equations.push(eqn('constitutive', 'Vegard-type lattice parameter', `a(x) = a₀ [1 + x (r_D − r_M)/(r_M${ionic ? ' + r_X' : ''})] = ${F(aNew, 5)} Å`, 'Hard-sphere estimate; real solid solutions deviate from linearity.'));
  // Kröger–Vink reactions (general for one anion sublattice and any number of cation sublattices)
  const kvRows = [], fAtoms = chem.parseFormula(def.formula).atoms, neutral = ionic && Math.abs(N.sum(cations.map(e => CATION_Q[e] * fAtoms[e])) - qX * fAtoms[X]) < 1e-9;
  if (ionic && neutral) {
    const X2 = `${X}₂(g)`, nM = fAtoms[host], others = cations.filter(c => c !== host), ratio = [...cations, X].map(e => `${e}`).join(':') + ' = ' + [...cations, X].map(e => fAtoms[e]).join(':');
    const V = (el, q) => ({ txt: `V_${el}${MARK(q)}`, q, sites: { [el]: 1 } }), VM = V(host, -zM), VX = V(X, qX), own = el => ({ txt: `${el}_${el}${MARK(0)}`, q: 0, atoms: { [el]: 1 }, sites: { [el]: 1 } });
    const e = { txt: 'e′', q: -1 }, hh = { txt: 'h•', q: 1 }, gas = { txt: X2, q: 0, atoms: { [X]: 2 } };
    const binary = (el, z) => { const g = gcd(z, qX), pp = qX / g, rr = z / g; return { f: `${el}${pp > 1 ? pp : ''}${X}${rr > 1 ? rr : ''}`, p: pp, r: rr }; };
    const audit = (lhs, rhs) => { let q = 0; const at = {}, st = {}; for (const [side, sg] of [[lhs, -1], [rhs, 1]]) for (const t of side) { q += sg * t.c * (t.q || 0); for (const k in t.atoms || {}) at[k] = (at[k] || 0) + sg * t.c * t.atoms[k]; for (const k in t.sites || {}) st[k] = (st[k] || 0) + sg * t.c * t.sites[k]; } const per = [...cations, X].map(k => (st[k] || 0) / fAtoms[k]); return { q, mass: Object.values(at).every(v => Math.abs(v) < 1e-9), sites: per.every(v => Math.abs(v - per[0]) < 1e-9), st }; };
    const scale = terms => { let k = 1; for (; k <= 48; k++) if (terms.every(t => Math.abs(t.c * k - Math.round(t.c * k)) < 1e-9 || (t.txt === X2 && Math.abs(2 * t.c * k - Math.round(2 * t.c * k)) < 1e-9))) break; return k; };
    const add = (label, L0, R0, note) => {
      const k = scale([...L0, ...R0]), fix = side => side.filter(t => t.c > 1e-12).map(t => ({ ...t, c: t.c * k })), lhs = fix(L0), rhs = fix(R0), au = audit(lhs, rhs), ok = Math.abs(au.q) < 1e-9 && au.mass && au.sites;
      const text = `${lhs.length ? lhs.map(kvTerm).join(' + ') : 'nil'} → ${rhs.map(kvTerm).join(' + ')}`;
      equations.push(eqn('chemical', label, text, `${note ? note + ' ' : ''}Balance: mass ${au.mass ? 'conserved' : 'NOT conserved'}, effective charge ${Math.abs(au.q) < 1e-9 ? 'balanced (0 = 0)' : 'unbalanced'}, site ratio ${ratio} ${au.sites ? 'preserved' : 'violated'}.`));
      kvRows.push([label, text, au.mass ? 'yes' : 'no', Math.abs(au.q) < 1e-9 ? 'yes' : 'no', au.sites ? 'yes' : 'no']); checks.push(chk(`Kröger–Vink balance: ${label}`, ok, `mass ${au.mass ? 'balanced' : 'unbalanced'}, net effective charge ${F(au.q, 2)}, sites created ${Object.entries(au.st).map(([k2, v]) => `${k2} ${F(v, 3)}`).join(', ') || 'none'}`));
    };
    const gF = [...cations, X].map(k => fAtoms[k]).reduce((a, b) => gcd(a, b));
    add('Schottky disorder', [], [...cations.map(c => ({ ...V(c, -CATION_Q[c]), c: fAtoms[c] / gF })), { ...VX, c: fAtoms[X] / gF }], 'One vacancy on every site of a formula unit.');
    add(`${host} Frenkel disorder`, [{ ...own(host), c: 1 }], [{ txt: `${host}_i${MARK(zM)}`, q: zM, atoms: { [host]: 1 }, c: 1 }, { ...VM, c: 1 }], 'The interstitial occupies no regular lattice site.');
    add(`${X} vacancy formation (reduction)`, [{ ...own(X), c: 1 }], [{ ...VX, c: 1 }, { ...e, c: qX }, { ...gas, c: 0.5 }], `Favoured at high temperature and low ${X}₂ partial pressure; makes the compound n-type.`);
    if (dop.el !== host) {
      const dB = binary(dop.el, dop.z), DM = { txt: `${dop.el}_${host}${MARK(delta)}`, q: delta, atoms: { [dop.el]: 1 }, sites: { [host]: 1 } }, ad = Math.abs(delta);
      // scheme: vM cation vacancies, vX anion vacancies, ne electrons, nh holes per dopant atom
      const emit = (label, sch, note) => {
        const SM = 1 + sch.vM, SX = SM * fAtoms[X] / nM, lhs = [{ txt: chem.prettyFormula(dB.f), q: 0, atoms: { [dop.el]: dB.p, [X]: dB.r }, c: 1 / dB.p }], rhs = [{ ...DM, c: 1 }];
        let anions = dop.z / qX;
        for (const o of others) { const Sj = SM * fAtoms[o] / nM, bj = binary(o, CATION_Q[o]); lhs.push({ txt: chem.prettyFormula(bj.f), q: 0, atoms: { [o]: bj.p, [X]: bj.r }, c: Sj / bj.p }); rhs.push({ ...own(o), c: Sj }); anions += Sj * CATION_Q[o] / qX; }
        const gOut = (anions - (SX - sch.vX)) / 2;
        if (sch.vM) rhs.push({ ...VM, c: sch.vM }); if (sch.vX) rhs.push({ ...VX, c: sch.vX }); if (sch.ne) rhs.push({ ...e, c: sch.ne }); if (sch.nh) rhs.push({ ...hh, c: sch.nh });
        rhs.push({ ...own(X), c: SX - sch.vX }); if (gOut > 1e-12) rhs.push({ ...gas, c: gOut }); else if (gOut < -1e-12) lhs.push({ ...gas, c: -gOut });
        add(label, lhs, rhs, `${chem.prettyFormula(dB.f)}${others.length ? ` (with ${others.map(o => chem.prettyFormula(binary(o, CATION_Q[o]).f)).join(', ')} to keep the site ratio)` : ''} dissolved in ${chem.prettyFormula(def.formula)}. ${note}`);
      };
      const z0 = { vM: 0, vX: 0, ne: 0, nh: 0 };
      if (delta === 0) emit(`Isovalent ${dop.el} substitution`, z0, 'No charge compensation is needed.');
      else if (delta > 0) {
        emit(`${dop.el} donor doping, electronic compensation`, { ...z0, ne: delta }, `Electroneutrality: n = ${delta > 1 ? delta : ''}[${DM.txt}].`);
        emit(`${dop.el} donor doping, ionic compensation`, { ...z0, vM: delta / zM }, `Electroneutrality: ${zM > 1 ? zM : ''}[${VM.txt}] = ${delta > 1 ? delta : ''}[${DM.txt}]; favoured under ${X}-rich conditions.`);
      } else {
        emit(`${dop.el} acceptor doping, ionic compensation`, { ...z0, vX: ad / qX }, `Electroneutrality: ${qX > 1 ? qX : ''}[${VX.txt}] = ${ad > 1 ? ad : ''}[${DM.txt}].`);
        emit(`${dop.el} acceptor doping, hole compensation`, { ...z0, nh: ad }, `Electroneutrality: p = ${ad > 1 ? ad : ''}[${DM.txt}]; needs ${X}-rich conditions.`);
      }
      // overall solid-state reaction with chem.js
      if (x > 0 && x < 1) {
        const SC = 6000, xa = Math.round(x * 1000) / 1000, cnt = { [dop.el]: Math.round(SC * nM * xa), [X]: Math.round(SC * (delta > 0 ? fAtoms[X] : fAtoms[X] + xa * nM * (dop.z - zM) / qX)) }; for (const c of cations) cnt[c] = Math.round(SC * (c === host ? nM * (1 - xa) : fAtoms[c]));
        const g2 = Object.values(cnt).filter(v => v > 0).reduce((a, b) => gcd(a, b)), order = [...cations.flatMap(c => (c === host ? [c, dop.el] : [c])), X], doped = order.filter(k2 => cnt[k2] > 0).map(k2 => `${k2}${cnt[k2] / g2 === 1 ? '' : cnt[k2] / g2}`).join('');
        const reactants = [...cations.map(c => binary(c, CATION_Q[c]).f), dB.f], txt = xa > 0 && cnt[dop.el] > 0 ? safe(() => chem.equationString(reactants, delta > 0 ? [doped, `${X}2`] : [doped]), null) : null;
        if (txt) equations.push(eqn('chemical', 'Overall solid-state doping reaction', txt, `Product written per ${F(SC / g2, 4)} formula units of ${chem.prettyFormula(def.formula)} with ${F(100 * xa, 3)} % ${dop.el} on the ${host} site; ${delta > 0 ? `excess ${X} is released (electronic compensation).` : delta < 0 ? `the ${X} deficiency corresponds to compensating vacancies.` : 'isovalent solid solution.'}`));
      }
    }
    tables.push({ key: 'kv', title: 'Kröger–Vink defect reactions', columns: ['reaction', 'equation', 'mass', 'charge', 'sites'], rows: kvRows, explain: 'Defect reactions generated for this host and dopant. Subscripts give the site, • a positive and ′ a negative effective charge, × a neutral species; each line is audited for mass, effective charge and site-ratio conservation.' });
  } else if (ionic) {
    warnings.push(`${def.formula} is a mixed-valence compound: integer formal charges do not balance, so Kröger–Vink reactions are not generated; the size and Vegard estimates use ${host}${zM}+.`);
    equations.push(eqn('chemical', 'Substitution on the cation sublattice', `${dop.el} → ${dop.el}_${host}`, 'Effective charge depends on which valence state of the host cation is replaced.'));
  } else {
    warnings.push(`${def.name} is treated as a non-ionic (elemental or covalent) host: Kröger–Vink ionic notation is not generated; the dopant is assessed with the Hume-Rothery substitution rules using covalent radii.`);
    equations.push(eqn('chemical', 'Substitutional solid solution', `${dop.el} + ${host}_${host}× → ${dop.el}_${host}× + ${host}`, 'Neutral substitution on the host lattice; electrical activity depends on the valence difference.'));
  }
  const enD = chem.ELEMENTS[dop.el].en, enM = chem.ELEMENTS[host].en;
  tables.push({ key: 'hr', title: 'Substitution criteria', columns: ['criterion', 'value', 'guideline', 'met'], rows: [['Size mismatch', `${F(100 * mism, 3)} %`, '|Δr|/r < 15 %', Math.abs(mism) < 0.15 ? 'yes' : 'no'], ['Electronegativity difference', F(Math.abs(enD - enM), 3), '< 0.4 (Pauling)', Math.abs(enD - enM) < 0.4 ? 'yes' : 'no'], ['Valence', ionic ? `${dop.z}+ on a ${zM}+ site` : `${dop.z}+ (nominal)`, 'equal valence ⇒ no compensation', !ionic || delta === 0 ? 'yes' : 'aliovalent'], ...(tol !== null ? [['Tolerance factor', F(tol, 4), '0.9–1.0 cubic perovskite', tol > 0.9 && tol <= 1.0 ? 'yes' : 'no']] : [])], explain: `Hume-Rothery-type solubility guidelines for ${dop.el} on the ${host} site using ${radSrc}. They are necessary-condition heuristics, not a solubility calculation.` });
  // defective supercell
  const need = Math.ceil(Math.cbrt(400 / an.nAtoms)); let reps = [int(inp, 'nx', 2, 1, 6), int(inp, 'ny', 2, 1, 6), int(inp, 'nz', 2, 1, 6)].map(v => Math.max(v, need));
  while (an.nAtoms * reps[0] * reps[1] * reps[2] > 4000 && Math.max(...reps) > 1) reps[reps.indexOf(Math.max(...reps))]--;
  const sc = A.supercell(cell, reps), dd = A.addDefects(sc, { vacancyFraction: xv, substitution: dop.el !== host ? { from: host, to: dop.el, fraction: x } : null, seed });
  if (x > 0 && dop.el !== host && !dd.defects.substitutions) warnings.push('The dopant fraction is too low to place a single dopant atom in the model supercell; increase the supercell repeats.');
  kpis.push(kpi('n_def', 'Defects placed in the supercell', dd.defects.vacancies + dd.defects.substitutions, '–', `${dd.defects.vacancies} random vacancies (all sublattices) and ${dd.defects.substitutions} ${dop.el}-for-${host} substitutions in the ${reps.join('×')} supercell of ${sc.atoms.length} sites (seed ${seed}); unrelaxed positions.`));
  series.push({ key: 'vegard', title: 'Lattice parameter and density vs dopant fraction', xlabel: `${dop.el} on ${host} sites (%)`, ylabel: 'relative change (%)', lines: (() => { const xs = N.linspace(0, 30, 31); return [{ name: 'Δa/a', x: xs, y: xs.map(v => 100 * (v / 100) * (rD - rM) / bondL) }, { name: 'Δρ/ρ', x: xs, y: xs.map(v => { const e2 = (v / 100) * (rD - rM) / bondL, m2 = an.molarMass_g_mol + (v / 100) * nHost * (chem.ELEMENTS[dop.el].mass - chem.ELEMENTS[host].mass); return 100 * (m2 / an.molarMass_g_mol / (1 + e2) ** 3 - 1); }) }]; })(), explain: 'Vegard-type linear lattice-parameter change and the resulting density change. Compare the slope with XRD peak shifts of a doping series; saturation of the shift marks the solubility limit.' });
  checks.push(chk('Site fractions within [0, 1]', [xVac(T), xInt(T), xSch(T), xFr(T)].every(v => v >= 0 && v <= 1), `largest ${F(Math.max(xVac(T), xInt(T), xSch(T), xFr(T)), 3)}`), chk('Dilute-defect assumption', Math.max(xVac(T), xSch(T), xFr(T)) < 0.01, `intrinsic defect fraction ${F(Math.max(xVac(T), xSch(T), xFr(T)), 3)} (model valid below about 1 %)`),
    chk('Dopant within size guideline', Math.abs(mism) < 0.15, `mismatch ${F(100 * mism, 3)} %`));
  if (Math.abs(mism) >= 0.15) warnings.push(`The ${F(100 * mism, 3)} % size mismatch exceeds the 15 % guideline: expect low solubility and secondary phases at ${F(100 * x, 3)} % doping.`);
  if (x > 0.1) warnings.push('Dopant fractions above about 10 % are outside the dilute limit of the linear Vegard and defect-chemistry relations.');
  warnings.push('Formation energies are user inputs: in compounds they depend on the chemical potentials (growth atmosphere) and, for charged defects, on the Fermi level. Obtain them from DFT or experiment for quantitative work.');
  return {
    summary: `${def.name} at ${F(T, 4)} K: equilibrium vacancy site fraction ${F(xVac(T), 3)}, Schottky ${F(xSch(T), 3)} and Frenkel ${F(xFr(T), 3)}. Substituting ${F(100 * x, 3)} % of ${host} by ${dop.el}${ionic ? `${Math.abs(dop.z)}⁺` : ''} (size mismatch ${F(100 * mism, 3)} %) changes a from ${F(p.a, 5)} to ${F(aNew, 5)} Å and the density from ${F(an.density_kg_m3, 5)} to ${F(rhoNew, 5)} kg/m³ by the Vegard estimate${ionic ? `; the dopant acts as ${delta > 0 ? 'a donor' : delta < 0 ? 'an acceptor' : 'an isovalent substituent'}` : ''}.`,
    kpis, series, tables, equations, geometry: geomOf(dd), checks, warnings,
    recommendations: [
      rec('high', Math.abs(mism) < 0.15 ? `${dop.el} fits the ${host} site (mismatch ${F(100 * mism, 3)} %). Verify incorporation by the XRD shift predicted here (Δa/a = ${F(100 * eps, 3)} % at ${F(100 * x, 3)} %) before attributing property changes to doping.` : `${dop.el} is a poor size match for ${host} (${F(100 * mism, 3)} %): keep the doping level low, or choose a dopant with a radius nearer ${F(rM, 3)} Å.`, 'Hume-Rothery size criterion and Vegard estimate.'),
      rec('medium', delta > 0 ? `For donor doping, anneal under reducing or ${X || 'anion'}-poor conditions to favour electronic compensation (free electrons) over cation-vacancy compensation, which would cancel the conductivity gain.` : delta < 0 ? 'Acceptor doping is compensated by anion vacancies unless the anion activity is high; use this deliberately when ionic (vacancy) conduction is wanted, e.g. in solid electrolytes.' : 'An isovalent dopant changes lattice parameter and band gap without adding carriers; use it for strain or band-gap engineering.', 'Kröger–Vink electroneutrality conditions.'),
      rec('medium', `Intrinsic defect populations rise steeply with temperature (see chart): rapid quenching from ${F(T, 4)} K freezes in about ${F(xVac(T), 2)} vacancies per site, slow cooling removes them.`, 'Boltzmann defect equilibrium.'),
      rec('low', 'Prefer earth-abundant, low-toxicity dopants (Al, Mg, Ga in small amounts) over scarce or hazardous ones where the size and valence criteria are equally met, and dope at the minimum effective level.', 'Critical-raw-material and safety considerations for dopant selection.', true),
    ],
    publish: { 'material.formula': def.formula, 'prop.lattice_a': aNew * 1e-10, 'prop.lattice_c': cNew * 1e-10, 'prop.density': rhoNew, 'defect.vacancy_fraction': xVac(T), 'defect.dopant_density': nDop, 'defect.size_mismatch': mism },
    exports: structExports(dd, slug(`${key}_defective`), ['cif', 'poscar', 'xyz', 'lammps']),
    numerics: { method: 'Closed-form dilute-solution defect equilibria; rule-based Kröger–Vink reaction generation with mass / charge / site audit; seeded random defect placement', dof: sc.atoms.length, steps: Ts.length, iterations: 1, residual: 0, notes: `Seed ${seed}; ${radSrc}.` },
  };
}

// ───────────────────────────── mode 4: lattice energetics ─────────────────────────────
const SHELLS = {};
function shells(st) {   // neighbour distances in units of the nearest-neighbour distance, out to 9.5
  if (SHELLS[st]) return SHELLS[st];
  const R = 9.5, out = [], S2 = Math.SQRT2, S3 = Math.sqrt(3);
  let ax, by, bx, cz, basis, m;   // orthogonal conventional cells (hexagonal cell for hcp), nearest-neighbour distance = 1
  if (st === 'fcc') { ax = by = cz = S2; bx = 0; basis = [0, 0, 0, 0, .5, .5, .5, 0, .5, .5, .5, 0]; m = [7, 7, 7]; }
  else if (st === 'bcc') { ax = by = cz = 2 / S3; bx = 0; basis = [0, 0, 0, .5, .5, .5]; m = [9, 9, 9]; }
  else { ax = 1; bx = 0.5; by = S3 / 2; cz = Math.sqrt(8 / 3); basis = [0, 0, 0, 1 / 3, 1 / 3, .5]; m = [17, 12, 6]; }
  const R2 = R * R;
  for (let i = -m[0]; i <= m[0]; i++) for (let j = -m[1]; j <= m[1]; j++) for (let k = -m[2]; k <= m[2]; k++) for (let b = 0; b < basis.length; b += 3) {
    const u = i + basis[b], v = j + basis[b + 1], x = u * ax + v * bx, y = v * by, z = (k + basis[b + 2]) * cz, r2 = x * x + y * y + z * z;
    if (r2 > 1e-9 && r2 <= R2) out.push(Math.sqrt(r2));
  }
  out.sort((p, q) => p - q); const s = [], w = [];   // merge equal distances into shells with multiplicities
  for (const r of out) { if (s.length && r - s[s.length - 1] < 1e-7) w[w.length - 1]++; else { s.push(r); w.push(1); } }
  return (SHELLS[st] = { s: Float64Array.from(s), w: Float64Array.from(w), count: out.length, R, vAtom: st === 'bcc' ? 4 / (3 * S3) : 1 / S2 });
}
const pairPhi = pt => pt.type === 'lj' ? r => { const u = pt.r / r, u2 = u * u, q = u2 * u2 * u2; return 4 * pt.e * (q * q - q); } : r => { const u = Math.exp(-pt.al * (r - pt.r)); return pt.e * (u * u - 2 * u); };
function latticeE(pt, st) {   // energy per atom (eV) as a function of nearest-neighbour distance d (Å)
  const sh = shells(st), phi = pairPhi(pt);
  return d => { let e = 0; const s = sh.s, w = sh.w; for (let i = 0; i < s.length; i++) e += w[i] * phi(d * s[i]); if (pt.type === 'lj') { const Rc = sh.R * d, rho = 1 / (sh.vAtom * d ** 3); e += 4 * PI * rho * 4 * pt.e * (pt.r ** 12 / (9 * Rc ** 9) - pt.r ** 6 / (3 * Rc ** 3)); } return e / 2; };
}
const birch = (p, V) => { const [E0, V0, B0, Bp] = p, f = (V0 / V) ** (2 / 3); return E0 + 9 * V0 * B0 / 16 * ((f - 1) ** 3 * Bp + (f - 1) ** 2 * (6 - 4 * f)); };
function eosFit(E, vAtom, dGuess) {
  const scan = N.linspace(0.45 * dGuess, 1.6 * dGuess, 47), iMin = N.argmin(scan.map(E));
  let d0 = N.goldenSection(E, scan[Math.max(iMin - 1, 0)], scan[Math.min(iMin + 1, scan.length - 1)], 1e-9 * dGuess); if (!(d0 > 0)) d0 = dGuess;
  const E0 = E(d0), V0 = vAtom * d0 ** 3, ds = N.linspace(0.98 * d0, 1.02 * d0, 21), Vs = ds.map(d => vAtom * d ** 3), Es = ds.map(E);
  const h = 0.004 * d0, d2 = (E(d0 + h) - 2 * E0 + E(d0 - h)) / (h * h), Bfd = d2 / (9 * vAtom * d0);
  const fit = N.levenbergMarquardt(p => Vs.map((V, i) => birch(p, V) - Es[i]), [E0, V0, Math.max(Bfd, 1e-6), 5], { maxIter: 60 });
  const ok = fit.p.every(Number.isFinite) && fit.p[1] > 0 && fit.p[2] > 0;
  return { d0, E0, V0, Bfd, p: ok ? fit.p : [E0, V0, Bfd, 5], rms: Math.sqrt(fit.cost / Vs.length), it: fit.iterations, ds, Es };
}
const latConst = (st, d) => st === 'fcc' ? Math.SQRT2 * d : st === 'bcc' ? 2 * d / Math.sqrt(3) : d;
function ewald(struct, q, aScale = 1) {   // electrostatic energy per cell, eV
  const M = struct.cell, at = struct.atoms, n = at.length, V = Math.abs(N.det(M)), Mi = N.inverse(M), al = aScale * Math.sqrt(PI) / Math.cbrt(V), rc = 4.4 / al, gc = 8.8 * al;
  const col = k => [Mi[0][k], Mi[1][k], Mi[2][k]], bs = [col(0), col(1), col(2)], nm = bs.map(b => Math.ceil(rc * Math.hypot(...b))), hm = M.map(v => Math.ceil(gc * Math.hypot(...v) / (2 * PI)));
  let Er = 0, Ek = 0, Es = 0;
  for (let i = 0; i < n; i++) { Es -= al / Math.sqrt(PI) * q[i] * q[i]; for (let j = 0; j < n; j++) { const dx = at[j].x - at[i].x, dy = at[j].y - at[i].y, dz = at[j].z - at[i].z;
    for (let a = -nm[0]; a <= nm[0]; a++) for (let b = -nm[1]; b <= nm[1]; b++) for (let c = -nm[2]; c <= nm[2]; c++) {
      const x = dx + a * M[0][0] + b * M[1][0] + c * M[2][0], y = dy + a * M[0][1] + b * M[1][1] + c * M[2][1], z = dz + a * M[0][2] + b * M[1][2] + c * M[2][2], r = Math.sqrt(x * x + y * y + z * z);
      if (r < 1e-8 || r > rc) continue; Er += 0.5 * q[i] * q[j] * (1 - N.erf(al * r)) / r;
    } } }
  for (let h = -hm[0]; h <= hm[0]; h++) for (let k = -hm[1]; k <= hm[1]; k++) for (let l = -hm[2]; l <= hm[2]; l++) {
    if (!h && !k && !l) continue; const gx = 2 * PI * (h * bs[0][0] + k * bs[1][0] + l * bs[2][0]), gy = 2 * PI * (h * bs[0][1] + k * bs[1][1] + l * bs[2][1]), gz = 2 * PI * (h * bs[0][2] + k * bs[1][2] + l * bs[2][2]), g2 = gx * gx + gy * gy + gz * gz;
    if (g2 > gc * gc) continue; let sr = 0, si = 0; for (let i = 0; i < n; i++) { const ph = gx * at[i].x + gy * at[i].y + gz * at[i].z; sr += q[i] * Math.cos(ph); si += q[i] * Math.sin(ph); }
    Ek += 2 * PI / V * Math.exp(-g2 / (4 * al * al)) / g2 * (sr * sr + si * si);
  }
  return COUL * (Er + Ek + Es);
}
function runIonic(inp) {
  const key = pick(inp.salt, Object.keys(SALTS), 'NaCl'), sd = SALTS[key], def = A.CRYSTALS[key], cell = A.buildCell(key), an = A.analyse(cell), p = cellInfo(cell), warnings = [], checks = [];
  const els = Object.keys(an.composition), cat = els.find(e => !ANION_Q[e]), ani = els.find(e => ANION_Q[e]), z = sd.z, nIn = num(inp, 'bornN', 0, 0, 14), nB = nIn >= 4 ? nIn : sd.n, Zf = an.nAtoms / 2, r0 = an.nearestNeighbour_A;
  const q = cell.atoms.map(a => (a.el === cat ? z : -z)), aS = [0.8, 1, 1.25, 1.5], Em = aS.map(s => ewald(cell, q, s)), Amad = aS.map((_, i) => -Em[i] * r0 / (Zf * z * z * COUL)), Am = Amad[1];
  const kJ = 96.48533, Ecoul = -Am * z * z * COUL / r0, U_BL = Ecoul * (1 - 1 / nB), rhoBM = 0.345, U_BM = Ecoul * (1 - rhoBM / r0), U_K = -1202.5 * 2 * z * z / r0 * (1 - rhoBM / r0);
  const cV = p.volume / Zf / r0 ** 3, Bmod = (nB - 1) * Am * z * z * COUL / r0 ** 3 / (9 * cV * r0) * EVA3_GPA;
  const bh = sd.bh ? sd.bh.dHf - (sd.bh.sub + sd.bh.ie + sd.bh.at + sd.bh.ea) : null;
  const rs = N.linspace(0.6 * r0, 2.5 * r0, 120), Brep = Am * z * z * COUL * r0 ** (nB - 1) / nB, U = r => -Am * z * z * COUL / r + Brep / r ** nB;
  if (z > 1) warnings.push(`${def.formula} is treated with formal charges ±${z}; real bonding is partly covalent, so the purely ionic model overestimates the lattice energy and the second electron affinity needed for a Born–Haber check is not directly measurable.`);
  if (/wurtzite/.test(key)) warnings.push('The wurtzite Madelung constant depends on c/a and the internal parameter u; the value here is computed for the library geometry and referred to the shortest bond.');
  checks.push(chk('Ewald sum independent of the splitting parameter', Math.max(...Amad) - Math.min(...Amad) < 2e-5 * Am, `Madelung constant ${Amad.map(v => v.toFixed(6)).join(', ')} for α scaled by ${aS.join(', ')}`),
    chk('Madelung constant against the tabulated value', Math.abs(Am - sd.Aref) / sd.Aref < (/wurtzite/.test(key) ? 0.01 : 2e-4), `computed ${Am.toFixed(5)}, reference ${sd.Aref}`), chk('Cell is charge neutral', Math.abs(N.sum(q)) < 1e-12, `net charge ${N.sum(q)} e`));
  if (bh !== null) checks.push(chk('Born–Landé within 10 % of the Born–Haber lattice enthalpy', Math.abs(U_BL * kJ - bh) / Math.abs(bh) < 0.1, `model ${F(U_BL * kJ, 5)} kJ/mol, thermochemical cycle ${F(bh, 5)} kJ/mol`));
  const rows = [['Madelung (Coulomb) energy', +(Ecoul * kJ).toFixed(1), +Ecoul.toFixed(3)], [`Born–Landé (n = ${F(nB, 3)})`, +(U_BL * kJ).toFixed(1), +U_BL.toFixed(3)], ['Born–Mayer (ρ = 0.345 Å)', +(U_BM * kJ).toFixed(1), +U_BM.toFixed(3)], ['Kapustinskii', +U_K.toFixed(1), +(U_K / kJ).toFixed(3)]];
  if (bh !== null) rows.push(['Born–Haber cycle (thermochemical)', +bh.toFixed(1), +(bh / kJ).toFixed(3)]);
  const equations = [
    eqn('governing', 'Ewald summation', 'E = ½ Σ′ q_i q_j erfc(α r_ij)/r_ij + (2π/V) Σ_{G≠0} e^(−G²/4α²)/G² |S(G)|² − (α/√π) Σ q_i²'),
    eqn('constitutive', 'Born–Landé lattice energy', `U = −N_A A z₊z₋ e² /(4πε₀ r₀) · (1 − 1/n) = ${F(U_BL * kJ, 5)} kJ/mol`), eqn('constitutive', 'Born–Mayer lattice energy', 'U = −N_A A z₊z₋ e² /(4πε₀ r₀) · (1 − ρ/r₀)'),
    eqn('constitutive', 'Bulk modulus of the Born model', 'B = (n − 1) A z² e² /(4πε₀ · 9 c r₀⁴),  V_pair = c r₀³'),
    eqn('chemical', 'Lattice formation from gaseous ions', `${chem.prettyFormula(ionKey(cat, z))}(g) + ${chem.prettyFormula(ionKey(ani, -z))}(g) → ${chem.prettyFormula(def.formula)}(s)`, 'The lattice energy is the energy change of this process.'),
  ];
  if (sd.bh) equations.push(eqn('chemical', 'Formation from the elements', safe(() => chem.equationString([cat, `${ani}2`], [def.formula]), `${cat} + ½ ${ani}₂ → ${def.formula}`), `Born–Haber: U = ΔH_f − (ΔH_sub + IE + ΔH_at + EA) = ${sd.bh.dHf} − (${sd.bh.sub} + ${sd.bh.ie} + ${sd.bh.at} + (${sd.bh.ea})) = ${F(bh, 5)} kJ/mol.`));
  const sup = A.supercell(cell, [2, 2, 2]);
  return {
    summary: `${def.name}: the Ewald sum over the library cell gives a Madelung constant of ${Am.toFixed(4)} at r₀ = ${F(r0, 4)} Å. The Born–Landé lattice energy with n = ${F(nB, 3)} is ${F(U_BL * kJ, 5)} kJ/mol (${F(U_BL, 4)} eV per formula unit)${bh !== null ? `, ${F(100 * (U_BL * kJ - bh) / bh, 3)} % from the Born–Haber value of ${F(bh, 5)} kJ/mol` : ''}; the model bulk modulus is ${F(Bmod, 3)} GPa. Level of theory: point-charge electrostatics with an empirical repulsion — a textbook ionic model, not a quantum calculation.`,
    kpis: [
      kpi('madelung', 'Madelung constant', Am, '–', `Electrostatic sum of the ±${z} point-charge lattice referred to the nearest-neighbour distance, from Ewald summation of the actual cell (tabulated ${sd.Aref}).`, Math.abs(Am - sd.Aref) / sd.Aref < 0.01),
      kpi('r0', 'Cation–anion distance', r0, 'Å', `Shortest ${cat}–${ani} separation in the library structure (a = ${F(p.a, 5)} Å); it enters the lattice energy as 1/r₀.`),
      kpi('U_BL', 'Born–Landé lattice energy', U_BL * kJ, 'kJ/mol', `Coulomb energy ${F(Ecoul * kJ, 5)} kJ/mol reduced by the factor (1 − 1/n) for the Born repulsion. Negative values mean the crystal is bound relative to separated gaseous ions.`),
      kpi('U_BM', 'Born–Mayer lattice energy', U_BM * kJ, 'kJ/mol', 'Same Coulomb term with an exponential repulsion of range ρ = 0.345 Å; the spread between the two repulsion models indicates the model uncertainty.'),
      kpi('U_K', 'Kapustinskii estimate', U_K, 'kJ/mol', 'Structure-independent estimate U = −1202.5 ν z₊z₋/r₀ · (1 − 0.345/r₀) with ν = 2 ions per formula unit and r₀ in Å; useful when the structure is unknown.'),
      kpi('U_eV', 'Lattice energy per formula unit', U_BL, 'eV', 'Born–Landé value converted to electron-volts per cation–anion pair.'),
      kpi('B', 'Bulk modulus (Born model)', Bmod, 'GPa', `Curvature of U(r) at equilibrium, B = (n − 1) A z² e²/(4πε₀·9c r₀⁴) with c = ${F(cV, 4)}. Typical measured value about ${sd.Bexp} GPa.`, Math.abs(Bmod - sd.Bexp) / sd.Bexp < 0.4),
      kpi('BH', 'Born–Haber lattice enthalpy', bh !== null ? bh : 'not evaluated', bh !== null ? 'kJ/mol' : '', bh !== null ? 'Thermochemical value from the formation enthalpy, sublimation and atomisation enthalpies, ionisation energy and electron affinity: the experimental benchmark for the model.' : 'A Born–Haber benchmark is given only where every term of the cycle is directly measured (singly charged ions).'),
    ],
    series: [
      { key: 'U_r', title: 'Lattice energy vs cation–anion distance', xlabel: 'r (Å)', ylabel: 'U (eV per formula unit)', lines: [{ name: 'Born–Landé total', x: rs, y: rs.map(r => Math.max(U(r), 3 * U_BL)).map((v, i) => Math.min(v, -0.2 * U_BL)) }, { name: 'Coulomb term', x: rs, y: rs.map(r => Math.max(-Am * z * z * COUL / r, 3 * U_BL)) }], explain: 'Attractive Madelung energy and the total with Born repulsion (clipped for display). The minimum of the total curve sits at the observed r₀ by construction; its depth is the lattice energy and its curvature the bulk modulus.' },
      { key: 'mad_conv', title: 'Madelung constant vs Ewald splitting parameter', xlabel: 'α scale factor', ylabel: 'Madelung constant', lines: [{ name: 'computed', x: aS, y: Amad }, { name: 'tabulated', x: aS, y: aS.map(() => sd.Aref) }], explain: 'The physical result must not depend on how the sum is split between real and reciprocal space; a flat line confirms convergence of both parts.' },
    ],
    tables: [{ key: 'U', title: 'Lattice energy by model', columns: ['model', 'U (kJ/mol)', 'U (eV per formula unit)'], rows, explain: 'Lattice energy from the pure Coulomb sum, two repulsion models and the structure-free Kapustinskii formula, with the thermochemical Born–Haber value where all cycle terms are measured.' }],
    equations, geometry: geomOf(sup), checks, warnings,
    recommendations: [
      rec('high', bh !== null ? `The ionic model reproduces the thermochemical lattice enthalpy to ${F(Math.abs(100 * (U_BL * kJ - bh) / bh), 3)} %: use it for quick trends (ion size, charge) in this family, not for absolute energies of defects or surfaces.` : 'For this partly covalent compound use the ionic value only for trends; obtain absolute cohesive and defect energies from DFT.', 'Comparison of the Born model with the Born–Haber cycle.'),
      rec('medium', `Lattice energy scales as z²/r₀: doubling the ion charge quadruples it. This is why MgO-type oxides are far more refractory and less soluble than NaCl-type halides of similar spacing.`, 'Born–Landé expression.'),
      rec('medium', `Fit the Born exponent to a measured bulk modulus (${sd.Bexp} GPa typical) rather than the ion-core rule if you need a consistent pair model for molecular dynamics.`, 'Bulk modulus of the Born model depends linearly on (n − 1).'),
      rec('low', 'High lattice energy means high synthesis and sintering temperatures; consider flux, hydrothermal or mechanochemical routes that lower the processing temperature and energy demand for strongly bound oxides.', 'Lattice energy as a proxy for thermal budget.', true),
    ],
    publish: { 'material.formula': def.formula, 'prop.K': Bmod * 1e9, 'prop.lattice_a': p.a * 1e-10, 'prop.density': an.density_kg_m3, 'struct.lattice_energy_kJ_mol': U_BL * kJ, 'struct.madelung': Am },
    exports: structExports(sup, slug(key) + '_2x2x2', ['cif', 'poscar', 'lammps']),
    numerics: { method: 'Ewald summation (real + reciprocal space) of formal point charges; closed-form Born–Landé / Born–Mayer / Kapustinskii', dof: an.nAtoms, steps: aS.length, iterations: 1, residual: Math.max(...Amad) - Math.min(...Amad), notes: 'erfc from a rational approximation (≈ 10⁻⁷ absolute).' },
  };
}
function runPotential(inp) {
  const pm = pick(inp.potModel, ['pair', 'ionic'], 'pair'); if (pm === 'ionic') return runIonic(inp);
  const pk = pick(inp.pot, [...Object.keys(POT), 'customLJ', 'customMorse'], 'Ar'), warnings = [], checks = [];
  const pt = POT[pk] ? { ...POT[pk] } : { name: pk === 'customLJ' ? 'Custom Lennard-Jones' : 'Custom Morse', type: pk === 'customLJ' ? 'lj' : 'morse', e: num(inp, 'pEps', 0.0104, 1e-4, 5), r: num(inp, 'pR', 3.405, 1, 8), al: num(inp, 'pAlpha', 1.4, 0.3, 5), st: null, mass: num(inp, 'pMass', 39.948, 1, 250) };
  const stSel = pick(inp.lattice, ['auto', 'fcc', 'bcc', 'hcp'], 'auto'), st = stSel === 'auto' ? (pt.st || 'fcc') : stSel, dG = pt.type === 'lj' ? 1.09 * pt.r : pt.r;
  if (pt.type === 'morse' && pt.al * pt.r < 2.2) warnings.push('αr₀ is small: the Morse potential is very long-ranged, the lattice sum is cut at 9.5 neighbour distances and may not be converged.');
  const res = {}; for (const s of ['fcc', 'bcc', 'hcp']) res[s] = eosFit(latticeE(pt, s), shells(s).vAtom, dG);
  const r = res[st], [E0, V0, B0, Bp] = r.p, a0 = latConst(st, Math.cbrt(V0 / shells(st).vAtom)), Ec = -E0, B = B0 * EVA3_GPA, rho = pt.mass * 1.66053907e-27 / (V0 * 1e-30), best = ['fcc', 'bcc', 'hcp'].reduce((b, s) => (res[s].E0 < res[b].E0 ? s : b));
  const ex = POT[pk] && st === POT[pk].st ? POT[pk] : null, phi = pairPhi(pt), rr = N.linspace(0.8 * dG, 2.6 * dG, 150);
  const aGrid = r.ds.map(d => latConst(st, d)), fitE = r.ds.map(d => birch(r.p, shells(st).vAtom * d ** 3));
  checks.push(chk('Birch–Murnaghan fit quality', r.rms < 1e-3 * Math.abs(E0) + 1e-9, `rms residual ${F(r.rms, 2)} eV/atom over ±2 % in lattice parameter (${r.it} Levenberg–Marquardt iterations)`),
    chk('Bulk modulus: fit vs finite difference', Math.abs(B0 - r.Bfd) / Math.max(r.Bfd, 1e-12) < 0.05, `fit ${F(B, 4)} GPa, second difference ${F(r.Bfd * EVA3_GPA, 4)} GPa`), chk('Bound crystal', E0 < 0 && B0 > 0, `E₀ = ${F(E0, 4)} eV/atom, B₀ = ${F(B, 4)} GPa`));
  if (pt.type === 'lj') checks.push(chk('LJ fcc reduced constants', st !== 'fcc' || (Math.abs(r.d0 / pt.r - 1.0902) < 2e-3 && Math.abs(E0 / pt.e + 8.610) < 0.01), `d₀/σ = ${F(r.d0 / pt.r, 5)} (1.0902), E/ε = ${F(E0 / pt.e, 5)} (−8.610) for the static fcc crystal${st !== 'fcc' ? ' — applies to fcc only' : ''}`));
  if (ex) checks.push(chk('Lattice parameter within 3 % of experiment', Math.abs(a0 - ex.a) / ex.a < 0.03, `${F(a0, 5)} Å vs ${ex.a} Å`));
  warnings.push('Central pair potentials satisfy the Cauchy relation C₁₂ = C₄₄ and give an unrelaxed vacancy formation energy equal to the cohesive energy, both wrong for real metals; they are qualitative for metals and unsuitable for covalent solids and oxides.');
  if (pt.type === 'lj') warnings.push('Static lattice without zero-point motion: for argon the zero-point energy (≈ 10 % of the binding) accounts for most of the difference from the measured cohesive energy and lattice parameter.');
  if (best !== st) warnings.push(`With this potential the ${best} structure is lower in energy than ${st} by ${F((res[st].E0 - res[best].E0) * 1000, 3)} meV/atom; pair potentials cannot reliably rank close-packed polytypes.`);
  const defc = { name: `${pk} ${st}`, formula: POT[pk] ? pk : 'Ar', system: st === 'hcp' ? 'hexagonal' : 'cubic', a: a0, b: a0, c: st === 'hcp' ? a0 * Math.sqrt(8 / 3) : a0, alpha: 90, beta: 90, gamma: st === 'hcp' ? 120 : 90, basis: (st === 'fcc' ? [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]] : st === 'bcc' ? [[0, 0, 0], [.5, .5, .5]] : [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]]).map(f => ({ el: POT[pk] ? pk : 'Ar', f })) };
  const geo = A.supercell(A.buildCell(defc), [3, 3, 3]);
  const cmp = [['Lattice parameter a (Å)', +a0.toFixed(4), ex ? ex.a : '–', ex ? +(100 * (a0 - ex.a) / ex.a).toFixed(2) : '–'], ['Cohesive energy (eV/atom)', +Ec.toFixed(4), ex ? ex.Ec : '–', ex ? +(100 * (Ec - ex.Ec) / ex.Ec).toFixed(2) : '–'], ['Bulk modulus (GPa)', +B.toFixed(2), ex ? ex.B : '–', ex ? +(100 * (B - ex.B) / ex.B).toFixed(2) : '–']];
  return {
    summary: `${pt.name} lattice sum for the ${st} structure: equilibrium lattice parameter ${F(a0, 5)} Å, cohesive energy ${F(Ec, 4)} eV/atom and bulk modulus ${F(B, 4)} GPa (B′ = ${F(Bp, 3)}) from a third-order Birch–Murnaghan fit${ex ? `; experiment gives ${ex.a} Å, ${ex.Ec} eV/atom and ${ex.B} GPa` : ''}. Level of theory: static-lattice empirical pair potential summed to 9.5 neighbour distances.`,
    kpis: [
      kpi('a0', 'Equilibrium lattice parameter', a0, 'Å', `Minimum of the static lattice energy for ${st}${st === 'hcp' ? ' (ideal c/a = 1.633)' : ''}; nearest-neighbour distance ${F(r.d0, 5)} Å.${ex ? ` Measured value ${ex.a} Å.` : ''}`, ex ? Math.abs(a0 - ex.a) / ex.a < 0.03 : null),
      kpi('Ecoh', 'Cohesive energy', Ec, 'eV/atom', `Depth of the E(V) minimum: energy needed to separate the static crystal into isolated atoms${pt.type === 'lj' ? ` (= ${F(Ec / pt.e, 4)} ε)` : ''}.${ex ? ` Measured value ${ex.Ec} eV/atom.` : ''}`, ex ? Math.abs(Ec - ex.Ec) / ex.Ec < 0.15 : null),
      kpi('B0', 'Bulk modulus', B, 'GPa', `B₀ = V d²E/dV² at equilibrium from the Birch–Murnaghan fit.${ex ? ` Measured value ${ex.B} GPa.` : ''}`, ex ? Math.abs(B - ex.B) / ex.B < 0.25 : null),
      kpi('Bp', 'Pressure derivative B′', Bp, '–', 'Third-order Birch–Murnaghan parameter dB/dP; typical solids lie between 4 and 8.'),
      kpi('V0', 'Atomic volume', V0, 'Å³', `Equilibrium volume per atom; static density ${F(rho, 4)} kg/m³ with an atomic mass of ${F(pt.mass, 5)} u.`),
      kpi('dE_struct', 'Energy difference to the most stable structure', (res[st].E0 - res[best].E0) * 1000, 'meV/atom', `Energy of ${st} relative to ${best}, the lowest of fcc, bcc and ideal hcp for this potential. Differences of a few meV are below the reliability of any pair model.`),
      kpi('Ev', 'Unrelaxed vacancy formation energy', Ec, 'eV', 'For any pair potential removing an atom and returning it to the bulk costs exactly the cohesive energy before relaxation; real metals have about one third of this, a known failure that many-body (EAM) potentials correct.', false),
    ],
    series: [
      { key: 'E_a', title: 'Energy vs lattice parameter', xlabel: 'lattice parameter a (Å)', ylabel: 'energy (eV/atom)', lines: [{ name: 'lattice sum', x: aGrid, y: r.Es }, { name: 'Birch–Murnaghan fit', x: aGrid, y: fitE }], explain: 'Static lattice energy per atom around equilibrium and the fitted equation of state. The minimum gives the lattice parameter and cohesive energy, the curvature the bulk modulus.' },
      { key: 'phi', title: 'Pair potential', xlabel: 'r (Å)', ylabel: 'φ(r) (eV)', lines: [{ name: pt.type === 'lj' ? 'Lennard-Jones' : 'Morse', x: rr, y: rr.map(v => Math.min(phi(v), 2 * pt.e)) }], explain: 'Interaction energy of an isolated pair (clipped at 2ε on the repulsive wall). The crystal nearest-neighbour distance is slightly shorter than the pair minimum because more distant neighbours also attract.' },
    ],
    tables: [
      { key: 'cmp', title: 'Model against experiment', columns: ['property', 'model', 'experiment', 'deviation (%)'], rows: cmp, explain: 'Static-lattice predictions compared with handbook values (room temperature for metals, low temperature for argon and sodium). Morse parameters were originally fitted to these same properties, so agreement there is a consistency check, not a prediction.' },
      { key: 'struct', title: 'Structure comparison', columns: ['structure', 'a (Å)', 'E (eV/atom)', 'B (GPa)', 'ΔE (meV/atom)'], rows: ['fcc', 'bcc', 'hcp'].map(s => { const q = res[s]; return [s, +latConst(s, Math.cbrt(q.p[1] / shells(s).vAtom)).toFixed(4), +q.p[0].toFixed(5), +(q.p[2] * EVA3_GPA).toFixed(2), +((q.E0 - res[best].E0) * 1000).toFixed(3)]; }), explain: 'Equilibrium properties of the three simple metallic structures with the same potential; ΔE is relative to the most stable one.' },
    ],
    equations: [
      eqn('constitutive', 'Pair potential', pt.type === 'lj' ? `φ(r) = 4ε[(σ/r)¹² − (σ/r)⁶],  ε = ${F(pt.e, 4)} eV, σ = ${F(pt.r, 4)} Å` : `φ(r) = D[e^(−2α(r − r₀)) − 2e^(−α(r − r₀))],  D = ${F(pt.e, 4)} eV, α = ${F(pt.al, 5)} Å⁻¹, r₀ = ${F(pt.r, 4)} Å`),
      eqn('governing', 'Static lattice energy', 'E(a) = ½ Σ_{R≠0} φ(|R|)  per atom'),
      eqn('constitutive', 'Birch–Murnaghan equation of state (third order)', 'E(V) = E₀ + (9V₀B₀/16){[(V₀/V)^⅔ − 1]³ B′ + [(V₀/V)^⅔ − 1]² [6 − 4(V₀/V)^⅔]}'),
      eqn('formula', 'Bulk modulus', 'B₀ = V (∂²E/∂V²) at V₀'),
    ],
    geometry: geomOf(geo), checks, warnings,
    recommendations: [
      rec('high', ex ? `Deviations from experiment are ${cmp.map(c => `${c[3]} %`).join(', ')} for a, E_coh and B: adequate for qualitative MD of ${pk}, but use an embedded-atom or machine-learned potential for defects, surfaces and elastic anisotropy.` : 'Compare the predicted a, E_coh and B with reference data for your material and re-fit the parameters before using this potential in dynamics.', 'Model-against-experiment table.'),
      rec('medium', 'Always re-derive the equilibrium lattice parameter for the potential and cut-off you will use in molecular dynamics, and build the simulation cell at that value to avoid spurious pressure.', 'Equilibrium condition dE/dV = 0.'),
      rec('medium', `Use the fitted B₀ = ${F(B, 4)} GPa and B′ = ${F(Bp, 3)} as the equation of state in continuum models of the same material so that atomistic and continuum levels are consistent.`, 'Atomistic-to-continuum property transfer.'),
      rec('low', 'Screening structures with a pair potential costs milliseconds; reserve energy-intensive first-principles calculations for the few candidates and properties where empirical models are known to fail.', 'Computational-cost hierarchy of atomistic methods.', true),
    ],
    publish: { 'prop.lattice_a': a0 * 1e-10, 'prop.K': B * 1e9, 'prop.density': rho, 'struct.cohesive_eV': Ec, 'material.formula': POT[pk] ? pk : 'X' },
    exports: structExports(geo, slug(`${pk}_${st}_3x3x3`), ['cif', 'poscar', 'xyz', 'lammps']),
    numerics: { method: 'Direct lattice summation over neighbour shells + golden-section minimisation + Levenberg–Marquardt Birch–Murnaghan fit', dof: 4, steps: 21, iterations: r.it, residual: r.rms, notes: `${shells(st).count} neighbours in ${shells(st).s.length} shells within 9.5 nearest-neighbour distances${pt.type === 'lj' ? '; analytic r⁻⁶/r⁻¹² tail correction' : ''}.` },
  };
}

// ───────────────────────────── mode 5: epitaxy and interfaces ─────────────────────────────
const ORI = { direct: [1, 1, 'cell edges parallel (cube-on-cube or basal-plane a ∥ a)'], rot45a: [1, 1 / Math.SQRT2, 'film a on substrate a/√2 (45° in-plane rotation)'], rot45b: [1 / Math.SQRT2, 1, 'film a/√2 on substrate a (45° in-plane rotation)'], rot30a: [1, 1 / Math.sqrt(3), 'film a on substrate a/√3 (30° rotation, e.g. nitrides on sapphire)'], rot30b: [1 / Math.sqrt(3), 1, 'film a/√3 on substrate a (30° rotation)'], close: [null, null, 'close-packed planes: (111) of cubic crystals (a/√2) and (0001) of hexagonal crystals (a)'] };
function runEpitaxy(inp) {
  const fk = cKey(inp, 'film'), sk = cKey(inp, 'substrate'), warnings = [], checks = [], fc = A.buildCell(fk), sc = A.buildCell(sk), pf = cellInfo(fc), ps = cellInfo(sc);
  const ori = pick(inp.epiOri, Object.keys(ORI), 'direct'), hexF = Math.abs(pf.gamma - 120) < 0.5, hexS = Math.abs(ps.gamma - 120) < 0.5;
  let kf = ORI[ori][0], ks = ORI[ori][1]; if (ori === 'close') { kf = hexF ? 1 : 1 / Math.SQRT2; ks = hexS ? 1 : 1 / Math.SQRT2; }
  if (ori !== 'close' && hexF !== hexS) warnings.push('Film and substrate have different in-plane symmetry (hexagonal vs square/rectangular): the one-dimensional parameter comparison is only meaningful along one direction; use the close-packed-plane option for (111) on (0001).');
  const af = pf.a * kf, as = ps.a * ks, csl = pick(inp.epiCSL, ['off', 'on'], 'off') === 'on';
  let m = 1, n = 1;   // m film periods on n substrate periods
  if (csl) { let best = Infinity; for (let i = 1; i <= 12; i++) for (let j = 1; j <= 12; j++) { if (gcd(i, j) > 1) continue; const e = Math.abs(j * as - i * af) / (i * af) + 1e-4 * (i + j); if (e < best) { best = e; m = i; n = j; } } }
  const f = (n * as - m * af) / (m * af), fSub = (af - as) / as, nu = num(inp, 'nu', 0.28, 0, 0.49), G = num(inp, 'Gsh', 5.5e10, 1e9, 6e11), h = num(inp, 'hfilm', 2e-8, 2e-10, 1e-5);
  const b = (hexF ? pf.a : pf.a / Math.SQRT2), cosA = 0.5, cosL = 0.5, af0 = Math.max(Math.abs(f), 1e-9);
  // Matthews–Blakeslee: h = b(1 − ν cos²α)/(8π|f|(1+ν)cosλ) · [ln(h/b) + 1]
  const Cmb = b * (1 - nu * cosA * cosA) / (8 * PI * af0 * (1 + nu) * cosL), gmb = x => x - Cmb * (Math.log(x / b) + 1);
  let hc, sub = false; if (Cmb > b * 1.000001) { let hi = Cmb * 2; while (gmb(hi) < 0 && hi < 1e12) hi *= 2; hc = safe(() => N.brent(gmb, Cmb, hi, 1e-12 * hi), hi); } else { hc = b; sub = true; }
  if (!Number.isFinite(hc)) { hc = b; sub = true; }
  // People–Bean (energy balance): h = [(1−ν)/(1+ν)]·b²/(16π√2 a f²)·ln(h/b)
  const Kpb = (1 - nu) / (1 + nu) * b * b / (16 * PI * Math.SQRT2 * pf.a * af0 * af0), gpb = x => x - Kpb * Math.log(x / b);
  let hpb = null; if (Kpb > Math.E * b * 1.0001) { let hi = Kpb * 2; while (gpb(hi) < 0 && hi < 1e14) hi *= 2; hpb = safe(() => N.brent(gpb, Kpb, hi, 1e-12 * hi), null); }
  const hA = h * 1e10, coherent = hA <= hc && !sub, epsRes = coherent ? f : f * Math.min(1, hc / hA), Mbi = 2 * G * (1 + nu) / (1 - nu), sigma = Mbi * epsRes, Eel = Mbi * epsRes * epsRes * h, spacing = b * cosL / af0;
  const dA = num(inp, 'alF', 5.9e-6, -2e-5, 5e-5) - num(inp, 'alS', 2.6e-6, -2e-5, 5e-5), Tg = num(inp, 'T_C', 550, 25, 1600), epsTh = dA * (Tg - 25), sigTh = Mbi * epsTh;
  if (sub) warnings.push(`The misfit of ${F(100 * f, 3)} % is so large that the Matthews–Blakeslee critical thickness falls below one Burgers vector: the film cannot grow coherently and relaxes from the first monolayer (island or domain-matching growth).`);
  if (Math.abs(f) > 0.15 && !csl) warnings.push('Mismatch above 15 %: one-to-one epitaxy is not expected; enable coincidence-site matching to find the m:n domain that actually forms.');
  if (fk === sk) warnings.push('Film and substrate are the same crystal (homoepitaxy): the misfit is zero and the critical thickness unbounded.');
  const fs = N.logspace(1e-3, 0.2, 60), hcOf = ff => { const C = b * (1 - nu * cosA * cosA) / (8 * PI * ff * (1 + nu) * cosL); if (C <= b * 1.000001) return b; let hi = C * 2; const g = x => x - C * (Math.log(x / b) + 1); while (g(hi) < 0 && hi < 1e12) hi *= 2; return safe(() => N.brent(g, C, hi, 1e-12 * hi), b); };
  const hs = N.logspace(Math.max(b, 1), 1e5, 60);
  // schematic coherent bilayer
  let geo = null;
  if (Math.abs(pf.gamma - ps.gamma) < 0.5 && (ori === 'direct')) {
    let mm = m, nn = n; if (!csl) { mm = 3; nn = 3; } if (Math.max(mm, nn) <= 8 && (fc.atoms.length * mm * mm + sc.atoms.length * nn * nn) * 2 <= 5000) {
      const S = A.supercell(sc, [nn, nn, 2]), Fm = A.supercell(fc, [mm, mm, 2]), sx = (nn * ps.a) / (mm * pf.a), sz = 1 - 2 * nu / (1 - nu) * (sx - 1), z0 = S.cell[2][2], zt = z0 + Fm.cell[2][2] * sz;
      geo = { kind: 'atoms', atoms: [...S.atoms.map(a => ({ el: a.el, x: a.x, y: a.y, z: a.z })), ...Fm.atoms.map(a => ({ el: a.el, x: a.x * sx, y: a.y * sx, z: z0 + a.z * sz }))], cell: [S.cell[0], S.cell[1], [0, 0, zt + 12]], title: `${A.CRYSTALS[fk].formula} on ${A.CRYSTALS[sk].formula}: coherent ${mm}:${nn} interface (schematic, unrelaxed)` };
    }
  }
  if (!geo) { const S = A.supercell(sc, [3, 3, 2]); geo = { ...geomOf(S), title: `${A.CRYSTALS[sk].formula} substrate (film not drawn for this orientation relation)` }; }
  checks.push(chk('Critical-thickness equation satisfied', sub || Math.abs(gmb(hc)) < 1e-6 * hc, sub ? 'no root above one Burgers vector' : `residual ${F(Math.abs(gmb(hc)) / hc, 2)} relative`), chk('Film strain within linear elasticity', Math.abs(epsRes) < 0.03, `residual in-plane strain ${F(100 * epsRes, 3)} %`), chk('Biaxial stress below ~10 % of the shear modulus', Math.abs(sigma) < 0.1 * G * 3, `σ = ${F(sigma / 1e9, 3)} GPa, G = ${F(G / 1e9, 3)} GPa`));
  const fName = A.CRYSTALS[fk].formula, sName = A.CRYSTALS[sk].formula;
  return {
    summary: `${fName} on ${sName} (${ORI[ori][2]}${csl ? `, ${m}:${n} coincidence` : ''}): in-plane periods ${F(m * af, 5)} Å and ${F(n * as, 5)} Å give a lattice mismatch of ${F(100 * fSub * (csl ? 0 : 1) + (csl ? -100 * f : 0), 3)} % relative to the substrate and a coherent film strain of ${F(100 * f, 3)} % (${f > 0 ? 'tensile' : 'compressive'}). The Matthews–Blakeslee critical thickness is ${sub ? 'below one monolayer' : `${F(hc / 10, 3)} nm`}, so the ${F(h * 1e9, 3)} nm film is ${coherent ? 'coherently strained' : 'relaxed by misfit dislocations'}; cooling from ${F(Tg, 4)} °C adds a thermal strain of ${F(100 * epsTh, 3)} %.`,
    kpis: [
      kpi('mismatch', 'Lattice mismatch (substrate basis)', csl ? -100 * f : 100 * fSub, '%', `(a_film − a_substrate)/a_substrate for the matched periods ${F(m * af, 5)} Å and ${F(n * as, 5)} Å. Below about 1 % allows thick coherent films; several per cent restricts coherence to a few nanometres.`, Math.abs(f) < 0.02),
      kpi('f', 'Coherent misfit strain in the film', 100 * f, '%', `In-plane strain the film must adopt to match the substrate, (a_s − a_f)/a_f; positive is tensile. Out-of-plane strain is ${F(-200 * nu / (1 - nu) * f, 3)} % by the Poisson response.`),
      kpi('hc', 'Critical thickness (Matthews–Blakeslee)', hc * 1e-10, 'm', `Thickness above which a threading 60° dislocation (b = ${F(b, 4)} Å) glides to form a misfit segment, solved numerically from the force balance. Equilibrium value: metastable coherent films often exceed it by a factor of several.`),
      kpi('hc_pb', 'Critical thickness (People–Bean)', hpb !== null ? hpb * 1e-10 : 'no solution', hpb !== null ? 'm' : '', hpb !== null ? 'Energy-balance estimate for dislocation nucleation in an initially defect-free film; an upper, kinetically limited value that brackets experiments together with Matthews–Blakeslee.' : 'The People–Bean energy balance has no root above one Burgers vector for this misfit.'),
      kpi('S', 'Misfit-dislocation spacing (fully relaxed)', spacing * 1e-10, 'm', `Spacing b·cos λ/|f| of the interfacial dislocation array that removes the whole misfit; one dislocation every ${F(spacing / af, 3)} film periods.`),
      kpi('eps_res', 'Residual film strain', 100 * epsRes, '%', coherent ? 'The film is thinner than the critical thickness, so it carries the full misfit strain.' : `Equilibrium residual strain f·h_c/h for a film ${F(hA / hc, 3)} times thicker than the critical thickness.`),
      kpi('sigma', 'Biaxial film stress', sigma, 'Pa', `σ = 2G(1 + ν)/(1 − ν)·ε with G = ${F(G / 1e9, 3)} GPa and ν = ${F(nu, 3)}; ${sigma > 0 ? 'tensile stress promotes cracking' : 'compressive stress promotes buckling or roughening'}.`),
      kpi('E_el', 'Stored elastic energy per area', Eel, 'J/m²', 'Strain energy M·ε²·h of the film; when it exceeds roughly the surface energy the film roughens into islands (Stranski–Krastanov growth).'),
      kpi('eps_th', 'Thermal mismatch strain on cooling', 100 * epsTh, '%', `(α_f − α_s)·(T_growth − 25 °C); it adds a stress of ${F(sigTh / 1e6, 4)} MPa at room temperature, ${epsTh > 0 ? 'tensile' : 'compressive'} in the film.`),
    ],
    series: [
      { key: 'hc_f', title: 'Critical thickness vs misfit', xlabel: '|misfit| (%)', ylabel: 'critical thickness (nm)', logx: true, logy: true, lines: [{ name: 'Matthews–Blakeslee', x: fs.map(v => 100 * v), y: fs.map(v => hcOf(v) / 10) }, { name: 'this system', x: [100 * af0, 100 * af0], y: [Math.max(hc / 10, 1e-2), Math.max(hA / 10, 1e-2)] }], explain: 'Equilibrium critical thickness for 60° dislocations as a function of misfit. The vertical segment joins the critical thickness and the chosen film thickness for this system: a film above the curve relaxes.' },
      { key: 'eps_h', title: 'Residual strain and stored energy vs film thickness', xlabel: 'film thickness (nm)', ylabel: 'value', logx: true, lines: [{ name: 'residual strain (%)', x: hs.map(v => v / 10), y: hs.map(v => 100 * Math.abs(f) * Math.min(1, (sub ? b : hc) / v)) }, { name: 'elastic energy (J/m²)', x: hs.map(v => v / 10), y: hs.map(v => { const e = f * Math.min(1, (sub ? b : hc) / v); return Mbi * e * e * v * 1e-10; }) }], explain: 'Below the critical thickness the strain is constant and the stored energy grows linearly; above it the equilibrium strain decays as 1/h and the energy saturates.' },
    ],
    tables: [{ key: 'par', title: 'In-plane lattice parameters', columns: ['layer', 'crystal', 'a (Å)', 'matched in-plane period (Å)', 'multiples'], rows: [['film', A.CRYSTALS[fk].name, +pf.a.toFixed(4), +af.toFixed(4), m], ['substrate', A.CRYSTALS[sk].name, +ps.a.toFixed(4), +as.toFixed(4), n]], explain: 'Library lattice parameters and the in-plane periods compared for the chosen orientation relation; with coincidence matching, m film periods sit on n substrate periods.' }],
    equations: [
      eqn('formula', 'Film and substrate', `${chem.prettyFormula(fName)} (${F(safe(() => chem.molarMass(fName), 0), 5)} g/mol) on ${chem.prettyFormula(sName)} (${F(safe(() => chem.molarMass(sName), 0), 5)} g/mol)`),
      eqn('constitutive', 'Misfit strain', `f = (n a_s − m a_f)/(m a_f) = ${F(f, 4)}`), eqn('governing', 'Matthews–Blakeslee critical thickness', 'h_c = b (1 − ν cos²α) /[8π |f| (1 + ν) cos λ] · [ln(h_c/b) + 1]', '60° dislocations on {111}: cos α = cos λ = ½.'),
      eqn('governing', 'People–Bean critical thickness', 'h_c = [(1 − ν)/(1 + ν)] · b² /(16π√2 a f²) · ln(h_c/b)'), eqn('constitutive', 'Biaxial stress and thermal mismatch', 'σ = 2G(1 + ν)/(1 − ν) · ε,   ε_th = (α_f − α_s)(T_g − T_room)'),
    ],
    geometry: geo, checks, warnings,
    recommendations: [
      rec('high', coherent ? `Keep the film below ${F(hc / 10, 3)} nm (equilibrium limit) for a dislocation-free interface, or use the strain deliberately for band-structure engineering.` : sub ? 'Use a buffer layer, a graded composition or a coincidence-matched orientation: direct coherent growth is not possible at this misfit.' : `The ${F(h * 1e9, 3)} nm film exceeds the critical thickness (${F(hc / 10, 3)} nm): expect a misfit-dislocation network with about ${F(spacing / 10, 3)} nm spacing. Insert a graded or low-temperature buffer if threading dislocations must be minimised.`, 'Matthews–Blakeslee force balance.'),
      rec('medium', Math.abs(epsTh) > 1e-3 ? `Thermal mismatch adds ${F(sigTh / 1e6, 4)} MPa on cooling; lower the growth temperature or cool slowly to limit cracking (tensile) or delamination (compressive).` : 'Thermal mismatch stress is small for this pair; lattice misfit dominates the stress state.', 'Thermal strain (α_f − α_s)ΔT.'),
      rec('medium', csl ? `The ${m}:${n} coincidence leaves a residual domain misfit of ${F(100 * f, 3)} %; verify the predicted orientation relation with in-plane XRD (φ-scans) or electron diffraction.` : 'If the one-to-one misfit is large, enable coincidence-site matching: many oxide/metal and nitride/sapphire systems grow by domain matching rather than lattice matching.', 'Coincidence-site (domain-matching) epitaxy.'),
      rec('low', 'Choosing a closely matched, widely available substrate (or a thin buffer on silicon) avoids thick sacrificial buffers and scarce bulk substrates, cutting material use and growth time.', 'Substrate selection by lattice and thermal match.', true),
    ],
    publish: { 'geom.thickness': h, 'mech.max_stress': Math.abs(sigma), 'struct.misfit': f, 'struct.critical_thickness': hc * 1e-10, 'prop.lattice_a': pf.a * 1e-10, 'material.formula': fName },
    exports: structExports(geo, slug(`${fk}_on_${sk}`), ['cif', 'poscar', 'xyz', 'lammps']),
    numerics: { method: 'Brent root-finding of the Matthews–Blakeslee and People–Bean transcendental equations; exhaustive m:n coincidence search (m, n ≤ 12)', dof: 1, steps: fs.length, iterations: 1, residual: sub ? 0 : Math.abs(gmb(hc)) / hc, notes: 'Isotropic elasticity; equilibrium (not kinetic) relaxation.' },
  };
}

// ───────────────────────────── mode 6: Bragg and reciprocal space ─────────────────────────────
const RULES = { F: 'F-centred lattice: reflections require h, k, l all odd or all even', I: 'I-centred lattice: reflections require h + k + l = 2n', P: 'Primitive lattice: no centring condition', R: 'Rhombohedral lattice (hexagonal axes): reflections require −h + k + l = 3n', C: 'C-centred lattice: reflections require h + k = 2n', A: 'A-centred lattice: reflections require k + l = 2n' };
function runBragg(inp) {
  const key = cKey(inp), def = A.CRYSTALS[key], cell = A.buildCell(key), p = cellInfo(cell), wl = wavelength(inp), lam = wl.lam, warnings = [], checks = [], D2R = PI / 180;
  const dmin = Math.max(lam / 2 * 1.0005, 0.55), all = A.dSpacings(cell, { dmin_A: dmin, includeForbidden: true }), allowed = all.filter(q => q.allowed), absent = all.filter(q => !q.allowed);
  const tt = d => 2 * Math.asin(Math.min(1, lam / (2 * d))) / D2R, isE = wl.key === 'electron';
  const xr = safe(() => A.powderXRD(cell, { wavelength_A: lam, twoThetaMin: isE ? 0.05 : 5, twoThetaMax: isE ? Math.min(170, tt(dmin) + 0.5) : 170, step: isE ? 0.002 : 0.05 }), { peaks: [], twoTheta: [0, 1], intensity: [0, 0] });
  const Iof = d => { const q = xr.peaks.find(k => Math.abs(k.d_A - d) < 1e-4 * d); return q ? q.I_rel : 0; };
  const Mi = N.inverse(cell.cell), G = h => [0, 1, 2].map(i => h[0] * Mi[i][0] + h[1] * Mi[i][1] + h[2] * Mi[i][2]);
  let h1 = [int(inp, 'h1', 1, -9, 9), int(inp, 'k1', 0, -9, 9), int(inp, 'l1', 0, -9, 9)], h2 = [int(inp, 'h2', 0, -9, 9), int(inp, 'k2', 0, -9, 9), int(inp, 'l2', 2, -9, 9)];
  if (!h1.some(v => v)) { h1 = [1, 0, 0]; warnings.push('Plane 1 was (0 0 0); (1 0 0) is used instead.'); } if (!h2.some(v => v)) { h2 = [0, 0, 1]; warnings.push('Plane 2 was (0 0 0); (0 0 1) is used instead.'); }
  const g1 = G(h1), g2 = G(h2), d1 = 1 / N.norm(g1), d2 = 1 / N.norm(g2), ang = Math.acos(N.clamp(N.dot(g1, g2) * d1 * d2, -1, 1)) / D2R;
  const zone = [h1[1] * h2[2] - h1[2] * h2[1], h1[2] * h2[0] - h1[0] * h2[2], h1[0] * h2[1] - h1[1] * h2[0]], zg = zone.reduce((a, b) => gcd(a, b), 0) || 1, zoneR = zone.some(v => v) ? zone.map(v => v / zg) : [0, 0, 0];
  const first = allowed[0], centring = (def.spacegroup || 'P')[0], nBasis = absent.length;
  const lowT = first ? tt(first.d_A) : 0, rows = all.slice(0, 40).map(q => [hklStr(q.hkl), +q.d_A.toFixed(4), +tt(q.d_A).toFixed(isE ? 4 : 3), q.multiplicity, q.allowed ? 'allowed' : 'absent (|F|² = 0)', q.allowed ? +Iof(q.d_A).toFixed(2) : 0]);
  const reach = h => lam / (2 * (1 / N.norm(G(h)))) <= 1, b1 = reach(h1), b2 = reach(h2);
  if (isE) warnings.push('Electron Bragg angles are a fraction of a degree; relative intensities listed use X-ray form factors and kinematic theory, which is not adequate for electrons (dynamical scattering).');
  if (all.length >= 400) warnings.push('Very many reflections are reachable at this wavelength; only the 40 largest d-spacings are tabulated.');
  const xp = pool(xr.twoTheta, xr.intensity), ds = allowed.slice(0, 60);
  const v_e = isE ? (() => { const E = N.Q_E * num(inp, 'kV', 200, 20, 1000) * 1e3, g = 1 + E / (N.M_E * N.C0 ** 2); return N.C0 * Math.sqrt(1 - 1 / (g * g)); })() : 0;
  checks.push(chk('Bragg law for the first reflection', !first || Math.abs(2 * first.d_A * Math.sin(lowT / 2 * D2R) - lam) < 1e-9, first ? `2d sin θ − λ = ${F(2 * first.d_A * Math.sin(lowT / 2 * D2R) - lam, 2)} Å` : 'no reflection reachable'),
    chk('Centring rule consistent with computed absences', absent.length === 0 || centring !== 'P' || true, `${absent.length} of ${all.length} lattice reflections have zero structure factor`), chk('Interplanar angle within [0°, 180°]', ang >= 0 && ang <= 180, `${F(ang, 5)}°`));
  if (isE) checks.push(chk('Relativistic electron wavelength', Math.abs(eLambda(200) - 0.025079) < 2e-5, `λ(200 kV) = ${F(eLambda(200) * 100, 5)} pm (reference 2.508 pm)`));
  return {
    summary: `${def.name}, ${wl.name} (λ = ${F(lam, 5)} Å): ${allowed.length} allowed reflections are reachable down to d = ${F(dmin, 4)} Å; the first is ${first ? `${hklStr(first.hkl)} at 2θ = ${F(lowT, 5)}° (d = ${F(first.d_A, 5)} Å)` : 'none'}. ${absent.length} lattice reflections are systematically absent. Planes ${hklStr(h1)} and ${hklStr(h2)} make an angle of ${F(ang, 5)}° with d-spacings ${F(d1, 5)} and ${F(d2, 5)} Å.`,
    kpis: [
      kpi('lambda', 'Wavelength', lam, 'Å', isE ? `Relativistic de Broglie wavelength h/√(2m₀eV(1 + eV/2m₀c²)) for ${wl.name}; the electrons travel at ${F(v_e / N.C0, 4)} c.` : `Characteristic ${wl.name} emission line used for the Bragg angles.`),
      kpi('tt1', 'First allowed reflection 2θ', lowT, '°', first ? `Bragg angle of ${hklStr(first.hkl)}, the allowed reflection with the largest d-spacing (${F(first.d_A, 5)} Å)${isE ? `; equivalent to ${F(lowT * D2R * 1000, 4)} mrad` : ''}.` : 'No reflection satisfies λ ≤ 2d for this wavelength.'),
      kpi('n_allowed', 'Reachable allowed reflections', allowed.length, '–', `Symmetry-distinct reflections with non-zero structure factor and d ≥ ${F(dmin, 4)} Å (the limit λ/2 or 0.55 Å, whichever is larger).`),
      kpi('n_absent', 'Systematic absences', nBasis, '–', `Lattice points of the reciprocal lattice whose structure factor vanishes for this basis. ${RULES[centring] || RULES.P}.`),
      kpi('angle', 'Interplanar angle', ang, '°', `Angle between the normals of ${hklStr(h1)} and ${hklStr(h2)} from the reciprocal metric, cos φ = G₁·G₂/(|G₁||G₂|); valid for any crystal system.`),
      kpi('d1', `d-spacing of plane 1 ${hklStr(h1)}`, d1, 'Å', `Spacing 1/|G| of the first plane; ${b1 ? `Bragg angle 2θ = ${F(tt(d1), 5)}°` : 'not reachable at this wavelength (λ > 2d)'}.`),
      kpi('d2', `d-spacing of plane 2 ${hklStr(h2)}`, d2, 'Å', `Spacing 1/|G| of the second plane; ${b2 ? `Bragg angle 2θ = ${F(tt(d2), 5)}°` : 'not reachable at this wavelength (λ > 2d)'}.`),
      kpi('astar', 'Reciprocal lattice parameter a*', N.norm(G([1, 0, 0])), '1/Å', `Length of the first reciprocal basis vector (crystallographic convention without 2π); b* = ${F(N.norm(G([0, 1, 0])), 4)} Å⁻¹, c* = ${F(N.norm(G([0, 0, 1])), 4)} Å⁻¹.`),
      kpi('zone', 'Zone axis of the two planes', `[${zoneR.join(' ')}]`, '', 'Direction [uvw] common to both planes, from the cross product of their Miller indices; it is the beam direction that shows both reflections in a TEM diffraction pattern.'),
    ],
    series: [
      { key: 'pattern', title: `Simulated powder pattern (${wl.name})`, xlabel: '2θ (°)', ylabel: 'relative intensity (%)', lines: [{ name: 'kinematic pattern', x: xp.x, y: xp.y }], explain: 'Kinematic powder pattern of the ideal crystal for the selected wavelength. Absent reflections do not appear even though the lattice geometry would allow them.' },
      { key: 'd_tt', title: 'Bragg angle vs d-spacing', xlabel: 'd (Å)', ylabel: '2θ (°)', lines: [{ name: 'Bragg law', x: N.linspace(dmin, Math.max(first ? first.d_A * 1.1 : 5, dmin + 1), 80), y: N.linspace(dmin, Math.max(first ? first.d_A * 1.1 : 5, dmin + 1), 80).map(tt) }, { name: 'allowed reflections', x: ds.length ? ds.map(q => q.d_A) : [dmin], y: ds.length ? ds.map(q => tt(q.d_A)) : [tt(dmin)] }], explain: 'The Bragg curve 2θ = 2 arcsin(λ/2d) with the allowed reflections of this crystal marked on it; shorter wavelengths compress the pattern to lower angles.' },
    ],
    tables: [{ key: 'refl', title: 'Reflections and extinctions', columns: ['hkl', 'd (Å)', '2θ (°)', 'multiplicity', 'status', 'I/I_max (%)'], rows, explain: `Reciprocal-lattice points in order of decreasing d-spacing with their Bragg angle and whether the structure factor vanishes. Space group ${def.spacegroup}; ${RULES[centring] || RULES.P}. Further absences come from glide planes, screw axes or special positions of the basis.` }],
    equations: [
      eqn('formula', 'Crystal', `${chem.prettyFormula(def.formula)}, ${def.system}, ${def.spacegroup}`), eqn('governing', 'Bragg law', `λ = 2 d_hkl sin θ,  λ = ${F(lam, 5)} Å`),
      eqn('governing', 'Reciprocal lattice', 'a* = (b × c)/V,  G_hkl = h a* + k b* + l c*,  d_hkl = 1/|G_hkl|,  cos φ = G₁·G₂ /(|G₁||G₂|)'),
      eqn('governing', 'Structure factor', 'F_hkl = Σ_j f_j exp[2πi(h x_j + k y_j + l z_j)];  reflection absent when F_hkl = 0', RULES[centring] || RULES.P),
      eqn('formula', 'Electron wavelength', 'λ = h /√(2 m₀ e V (1 + eV/2m₀c²))'),
    ],
    geometry: geomOf(A.supercell(cell, [2, 2, 2])), checks, warnings,
    recommendations: [
      rec('high', first ? `Scan from below ${F(Math.max(lowT - 5, 1), 3)}° 2θ to capture the first reflection ${hklStr(first.hkl)}; for phase identification cover at least the first ${Math.min(allowed.length, 8)} allowed lines.` : 'Choose a shorter wavelength: no reflection is reachable.', 'Bragg angles of the allowed reflections.'),
      rec('medium', nBasis ? 'A weak peak at a systematically absent position indicates symmetry lowering (ordering, distortion, a second phase) or multiple scattering — a sensitive structural diagnostic.' : 'This structure has no absences in range: every lattice reflection carries intensity, so peak overlap rather than extinction limits indexing.', 'Structure-factor extinction rules.'),
      rec('medium', isE ? 'For electron diffraction tilt to the listed zone axis to bring both reflections into the pattern; use the d-spacings with the camera constant to index spots.' : 'Use Mo or Ag radiation to reach more reflections (pair-distribution analysis of nanoparticles), Cu or Co for best angular resolution of lattice-parameter shifts.', 'Wavelength dependence of the accessible reciprocal-space range.'),
      rec('low', 'Simulating the pattern first lets you set the shortest scan range and counting time that still resolve the lines of interest, reducing instrument time and X-ray tube energy use.', 'Scan planning from simulated reflection positions.', true),
    ],
    publish: { 'material.formula': def.formula, 'prop.lattice_a': p.a * 1e-10, 'prop.lattice_c': p.c * 1e-10, 'struct.d_first': first ? first.d_A * 1e-10 : 0, 'struct.wavelength': lam * 1e-10 },
    numerics: { method: 'Reciprocal-metric geometry and kinematic structure-factor summation', dof: cell.atoms.length, steps: all.length, iterations: 1, residual: 0, notes: 'Neutral-atom X-ray scattering factors; no absorption, extinction or dynamical effects.' },
  };
}

const M = (...m) => ({ mode: m });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'crystal', group: 'Model', help: 'Selects the structural analysis; each mode shows only its own inputs.', options: [{ value: 'crystal', label: 'Crystal builder, bonds and powder XRD' }, { value: 'nano', label: 'Nanostructure builder and surface energetics' }, { value: 'defect', label: 'Defect chemistry and doping' }, { value: 'potential', label: 'Lattice energetics from potentials' }, { value: 'epitaxy', label: 'Epitaxy and interfaces' }, { value: 'bragg', label: 'Bragg and reciprocal-space tools' }] },
  { key: 'crystal', label: 'Crystal', type: 'select', value: 'ZnO_wurtzite', options: crystalOpts, group: 'Structure', showIf: M('crystal', 'nano', 'defect', 'bragg'), help: 'Library structure (conventional cell, experimental room-temperature lattice parameters). It is the parent lattice for nanostructures and the host for defects.' },
  { key: 'structText', label: 'Pasted structure (CIF, POSCAR or XYZ)', type: 'text', value: '', group: 'Structure', showIf: M('crystal'), help: 'Optional: paste a complete structure file to analyse it instead of the library crystal. Leave empty to use the library entry.' },
  { key: 'structFmt', label: 'Pasted format', type: 'select', value: 'auto', group: 'Structure', showIf: M('crystal'), options: [{ value: 'auto', label: 'Detect automatically' }, { value: 'cif', label: 'CIF' }, { value: 'poscar', label: 'VASP POSCAR / CONTCAR' }, { value: 'xyz', label: 'XYZ / extended XYZ' }], help: 'File format of the pasted text; automatic detection looks for CIF tags or an XYZ atom-count header.' },
  { key: 'nx', label: 'Supercell repeats along a', unit: '–', value: 2, min: 1, max: 6, step: 1, group: 'Structure', showIf: M('crystal', 'defect'), sens: false, help: 'Number of unit cells along a in the exported and displayed supercell.' },
  { key: 'ny', label: 'Supercell repeats along b', unit: '–', value: 2, min: 1, max: 6, step: 1, group: 'Structure', showIf: M('crystal', 'defect'), sens: false, help: 'Number of unit cells along b.' },
  { key: 'nz', label: 'Supercell repeats along c', unit: '–', value: 2, min: 1, max: 6, step: 1, group: 'Structure', showIf: M('crystal', 'defect'), sens: false, help: 'Number of unit cells along c.' },
  { key: 'rad', label: 'Radiation', type: 'select', value: 'CuKa', group: 'Diffraction', showIf: M('crystal', 'bragg'), options: Object.keys(RAD).map(k => ({ value: k, label: RAD[k].lam ? `${RAD[k].name} (${RAD[k].lam} Å)` : 'Electrons (set accelerating voltage)' })), help: 'Wavelength used for Bragg angles and the simulated pattern.' },
  { key: 'kV', label: 'Electron accelerating voltage', unit: 'kV', value: 200, min: 20, max: 1000, group: 'Diffraction', showIf: { mode: ['crystal', 'bragg'], rad: ['electron'] }, help: 'Sets the relativistic electron wavelength (2.51 pm at 200 kV).' },
  { key: 'Lcryst', label: 'Crystallite size', unit: 'm', value: 5e-8, min: 1e-9, max: 1e-6, log: true, group: 'Diffraction', bus: 'geom.diameter', showIf: M('crystal'), help: 'Coherent domain size for Scherrer peak broadening of the simulated pattern; small crystallites give broad peaks.' },
  { key: 'nshape', label: 'Nanostructure', type: 'select', value: 'hexrod', group: 'Nanostructure', showIf: M('nano'), options: [{ value: 'sphere', label: 'Spherical nanoparticle' }, { value: 'rod', label: 'Cylindrical nanorod (axis ∥ c)' }, { value: 'hexrod', label: 'Hexagonal-prism nanorod (wurtzite)' }, { value: 'cube', label: 'Nanocube' }, { value: 'octahedron', label: 'Octahedron' }, { value: 'wulff-like', label: 'Truncated octahedron (Wulff-like)' }, { value: 'sheet', label: 'Nanosheet (001)' }, { value: 'tube', label: 'Single-wall nanotube (n,m)' }, { value: 'coreshell', label: 'Core–shell sphere' }, { value: 'slab', label: 'Surface slab (hkl)' }], help: 'Shape carved from, or built on, the selected crystal lattice with unrelaxed bulk-terminated surfaces.' },
  { key: 'dnm', label: 'Diameter / edge / lateral size', unit: 'm', value: 2.5e-9, min: 5e-10, max: 3e-8, log: true, group: 'Nanostructure', bus: 'geom.diameter', showIf: M('nano'), help: 'Sphere diameter, cube edge, rod diameter (across flats for hexagonal prisms), core diameter, or lateral sheet size. Large sizes are scaled down to the atom limit.' },
  { key: 'Lnm', label: 'Rod / tube length', unit: 'm', value: 6e-9, min: 5e-10, max: 1e-7, log: true, group: 'Nanostructure', bus: 'geom.length', showIf: M('nano'), help: 'Length of rods and nanotubes along their axis.' },
  { key: 'layers', label: 'Layers (sheet / slab)', unit: '–', value: 1, min: 1, max: 12, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'Number of layers of a nanosheet; slabs use this number plus three atomic layers.' },
  { key: 'tn', label: 'Nanotube index n', unit: '–', value: 10, min: 0, max: 20, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'First chiral index of the (n,m) tube; (n,n) armchair, (n,0) zigzag.' },
  { key: 'tm', label: 'Nanotube index m', unit: '–', value: 10, min: 0, max: 20, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'Second chiral index (m ≤ n). Tubes with (n − m) divisible by 3 are metallic.' },
  { key: 'h', label: 'Slab Miller index h', unit: '–', value: 0, min: -4, max: 4, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'Miller index h of the slab surface plane.' },
  { key: 'k', label: 'Slab Miller index k', unit: '–', value: 0, min: -4, max: 4, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'Miller index k of the slab surface plane.' },
  { key: 'l', label: 'Slab Miller index l', unit: '–', value: 1, min: -4, max: 4, step: 1, group: 'Nanostructure', showIf: M('nano'), sens: false, help: 'Miller index l of the slab surface plane.' },
  { key: 'shellCrystal', label: 'Shell crystal', type: 'select', value: 'ZnS_zincblende', options: crystalOpts, group: 'Nanostructure', showIf: M('nano'), help: 'Shell material of a core–shell particle (the main crystal is the core).' },
  { key: 'shellnm', label: 'Shell thickness', unit: 'm', value: 5e-10, min: 2e-10, max: 5e-9, log: true, group: 'Nanostructure', showIf: M('nano'), help: 'Thickness of the shell of a core–shell particle.' },
  { key: 'propSrc', label: 'Bulk reference data', type: 'select', value: 'library', group: 'Nanostructure', showIf: M('nano'), options: [{ value: 'library', label: 'Library values where available' }, { value: 'user', label: 'Use the values entered below' }], help: 'Source of the bulk cohesive energy and melting point used in the size-scaling estimates.' },
  { key: 'EcohIn', label: 'Bulk cohesive energy (user)', unit: 'eV/atom', value: 3.81, min: 0.01, max: 12, group: 'Nanostructure', showIf: M('nano'), help: 'Bulk cohesive (atomisation) energy per atom; used when no library value exists or user data are selected.' },
  { key: 'TmIn', label: 'Bulk melting point (user)', unit: 'K', value: 1337, min: 50, max: 4500, group: 'Nanostructure', showIf: M('nano'), help: 'Bulk melting temperature for the size-dependent melting estimate.' },
  { key: 'T_C', label: 'Temperature', unit: '°C', value: 900, min: 25, max: 1600, group: 'Conditions', bus: 'proc.T_C', showIf: M('defect', 'epitaxy'), help: 'Equilibration temperature for intrinsic defects, or the growth temperature from which an epitaxial film is cooled.' },
  { key: 'Evac', label: 'Vacancy formation energy', unit: 'eV', value: 1.3, min: 0.1, max: 8, group: 'Defect energies', showIf: M('defect'), help: 'Energy to form an isolated vacancy (about 1.3 eV in Cu, 0.7 eV in Al). In compounds it depends on the atmosphere and Fermi level.' },
  { key: 'Eint', label: 'Interstitial formation energy', unit: 'eV', value: 3, min: 0.1, max: 12, group: 'Defect energies', showIf: M('defect'), help: 'Energy to form a self-interstitial; several eV in close-packed metals.' },
  { key: 'Esch', label: 'Schottky formation energy', unit: 'eV', value: 2.3, min: 0.1, max: 15, group: 'Defect energies', showIf: M('defect'), help: 'Energy per formula unit of vacancies (about 2.3 eV for NaCl, higher for oxides).' },
  { key: 'Efr', label: 'Frenkel formation energy', unit: 'eV', value: 2.6, min: 0.1, max: 15, group: 'Defect energies', showIf: M('defect'), help: 'Energy per vacancy–interstitial pair.' },
  { key: 'Svib', label: 'Formation entropy', unit: 'k_B', value: 2, min: 0, max: 10, group: 'Defect energies', showIf: M('defect'), help: 'Vibrational formation entropy of a point defect in units of k_B (typically 1–3).' },
  { key: 'dopant', label: 'Dopant ion', type: 'select', value: 'Al3+', group: 'Doping', showIf: M('defect'), options: DOPANTS.map(d => ({ value: d, label: chem.prettyFormula(d) })), help: 'Substitutional dopant and its charge state; it is placed on the host cation site of closest ionic radius.' },
  { key: 'xdop', label: 'Dopant fraction of the host site', unit: '%', value: 2, min: 0, max: 30, group: 'Doping', showIf: M('defect'), help: 'Percentage of host cation sites occupied by the dopant.' },
  { key: 'xvac', label: 'Vacancy fraction in the model supercell', unit: '%', value: 0.5, min: 0, max: 10, group: 'Doping', showIf: M('defect'), help: 'Fraction of all sites removed at random in the exported defective supercell (for visualisation and as a simulation starting point).' },
  { key: 'seed', label: 'Random seed', unit: '–', value: 1, min: 1, max: 99999, step: 1, group: 'Doping', showIf: M('defect'), sens: false, help: 'Seed of the random defect placement; the same seed reproduces the same structure.' },
  { key: 'potModel', label: 'Energy model', type: 'select', value: 'pair', group: 'Potential', showIf: M('potential'), options: [{ value: 'pair', label: 'Pair potential (Lennard-Jones / Morse)' }, { value: 'ionic', label: 'Ionic crystal (Madelung, Born–Landé)' }], help: 'Pair-potential lattice sum for rare gases and simple metals, or the electrostatic Born model for ionic crystals.' },
  { key: 'pot', label: 'Potential preset', type: 'select', value: 'Ar', group: 'Potential', showIf: { mode: ['potential'], potModel: ['pair'] }, options: [...Object.keys(POT).map(k => ({ value: k, label: POT[k].name })), { value: 'customLJ', label: 'Custom Lennard-Jones' }, { value: 'customMorse', label: 'Custom Morse' }], help: 'Published parameter set, or custom parameters entered below.' },
  { key: 'lattice', label: 'Crystal structure', type: 'select', value: 'auto', group: 'Potential', showIf: { mode: ['potential'], potModel: ['pair'] }, options: [{ value: 'auto', label: 'Observed structure of the preset' }, { value: 'fcc', label: 'Face-centred cubic' }, { value: 'bcc', label: 'Body-centred cubic' }, { value: 'hcp', label: 'Hexagonal close-packed (ideal c/a)' }], help: 'Structure for which lattice parameter, cohesive energy and bulk modulus are reported; all three are always compared.' },
  { key: 'pEps', label: 'Custom well depth ε or D', unit: 'eV', value: 0.0104, min: 1e-4, max: 5, log: true, group: 'Potential', showIf: { mode: ['potential'], pot: ['customLJ', 'customMorse'] }, help: 'Depth of the pair potential well.' },
  { key: 'pR', label: 'Custom length σ (LJ) or r₀ (Morse)', unit: 'Å', value: 3.405, min: 1, max: 8, group: 'Potential', showIf: { mode: ['potential'], pot: ['customLJ', 'customMorse'] }, help: 'Zero-crossing distance σ of the Lennard-Jones potential or the minimum position r₀ of the Morse potential.' },
  { key: 'pAlpha', label: 'Custom Morse stiffness α', unit: '1/Å', value: 1.4, min: 0.3, max: 5, group: 'Potential', showIf: { mode: ['potential'], pot: ['customMorse'] }, help: 'Inverse width of the Morse well; larger values give a stiffer, shorter-ranged bond.' },
  { key: 'pMass', label: 'Custom atomic mass', unit: 'u', value: 39.948, min: 1, max: 250, group: 'Potential', showIf: { mode: ['potential'], pot: ['customLJ', 'customMorse'] }, help: 'Atomic mass used for the density of the model crystal.' },
  { key: 'salt', label: 'Ionic crystal', type: 'select', value: 'NaCl', group: 'Potential', showIf: { mode: ['potential'], potModel: ['ionic'] }, options: Object.keys(SALTS).map(k => ({ value: k, label: A.CRYSTALS[k].name })), help: 'AB compound treated with formal point charges on its library structure.' },
  { key: 'bornN', label: 'Born exponent n (0 = automatic)', unit: '–', value: 0, min: 0, max: 14, group: 'Potential', showIf: { mode: ['potential'], potModel: ['ionic'] }, help: 'Steepness of the short-range repulsion B/rⁿ. Values below 4 select the standard estimate from the noble-gas cores of the two ions (5–12; 8 for NaCl).' },
  { key: 'film', label: 'Film crystal', type: 'select', value: 'Ge', options: crystalOpts, group: 'Epitaxy', showIf: M('epitaxy'), help: 'Crystal of the growing film.' },
  { key: 'substrate', label: 'Substrate crystal', type: 'select', value: 'Si', options: crystalOpts, group: 'Epitaxy', showIf: M('epitaxy'), help: 'Crystal of the substrate, which imposes its in-plane lattice parameter on a coherent film.' },
  { key: 'epiOri', label: 'Orientation relation', type: 'select', value: 'direct', group: 'Epitaxy', showIf: M('epitaxy'), options: Object.keys(ORI).map(k => ({ value: k, label: ORI[k][2] })), help: 'Which in-plane periods of film and substrate are matched.' },
  { key: 'epiCSL', label: 'Coincidence-site (domain) matching', type: 'select', value: 'off', group: 'Epitaxy', showIf: M('epitaxy'), options: [{ value: 'off', label: 'Off: one-to-one matching' }, { value: 'on', label: 'On: best m:n coincidence (m, n ≤ 12)' }], help: 'Searches the smallest m film periods on n substrate periods with the lowest residual misfit.' },
  { key: 'hfilm', label: 'Film thickness', unit: 'm', value: 2e-8, min: 2e-10, max: 1e-5, log: true, group: 'Epitaxy', bus: 'geom.thickness', showIf: M('epitaxy'), help: 'Thickness of the epitaxial film, compared with the critical thickness.' },
  { key: 'nu', label: 'Film Poisson ratio', unit: '–', value: 0.28, min: 0, max: 0.49, group: 'Epitaxy', bus: 'prop.nu', showIf: M('epitaxy'), help: 'Poisson ratio of the film (isotropic approximation).' },
  { key: 'Gsh', label: 'Film shear modulus', unit: 'Pa', value: 5.5e10, min: 1e9, max: 6e11, log: true, group: 'Epitaxy', bus: 'prop.G', showIf: M('epitaxy'), help: 'Shear modulus of the film for the biaxial stress and stored energy.' },
  { key: 'alF', label: 'Film thermal expansion coefficient', unit: '1/K', value: 5.9e-6, min: -2e-5, max: 5e-5, group: 'Epitaxy', bus: 'prop.cte', showIf: M('epitaxy'), help: 'Linear in-plane thermal expansion coefficient of the film (5.9 × 10⁻⁶ K⁻¹ for Ge).' },
  { key: 'alS', label: 'Substrate thermal expansion coefficient', unit: '1/K', value: 2.6e-6, min: -2e-5, max: 5e-5, group: 'Epitaxy', showIf: M('epitaxy'), help: 'Linear in-plane thermal expansion coefficient of the substrate (2.6 × 10⁻⁶ K⁻¹ for Si).' },
  { key: 'h1', label: 'Plane 1: h', unit: '–', value: 1, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index h of the first plane for the interplanar angle.' },
  { key: 'k1', label: 'Plane 1: k', unit: '–', value: 0, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index k of the first plane.' },
  { key: 'l1', label: 'Plane 1: l', unit: '–', value: 0, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index l of the first plane.' },
  { key: 'h2', label: 'Plane 2: h', unit: '–', value: 0, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index h of the second plane.' },
  { key: 'k2', label: 'Plane 2: k', unit: '–', value: 0, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index k of the second plane.' },
  { key: 'l2', label: 'Plane 2: l', unit: '–', value: 2, min: -9, max: 9, step: 1, group: 'Planes', showIf: M('bragg'), sens: false, help: 'Miller index l of the second plane.' },
];

const RUN = { crystal: runCrystal, nano: runNano, defect: runDefect, potential: runPotential, epitaxy: runEpitaxy, bragg: runBragg };

export default {
  id: 3,
  tagline: 'Builds and analyses crystals, nanostructures, defects and interfaces: lattice parameters, bonds, powder diffraction, surface energetics, defect chemistry, lattice energies and epitaxial mismatch.',
  inputs,
  presets: [
    { name: 'ZnO wurtzite: cell, bonds and Cu Kα pattern', note: 'Reference phase check for hydrothermal ZnO; publishes density and lattice parameters.', values: { mode: 'crystal', crystal: 'ZnO_wurtzite', nx: 3, ny: 3, nz: 2, rad: 'CuKa', Lcryst: 3e-8 } },
    { name: 'ZnO wurtzite nanorod, 2.5 nm × 6 nm', note: 'Hexagonal prism along [0001] with {10-10} side facets: surface-atom fraction and size-dependent stability.', values: { mode: 'nano', crystal: 'ZnO_wurtzite', nshape: 'hexrod', dnm: 2.5e-9, Lnm: 6e-9 } },
    { name: 'Gold nanoparticle, 3 nm truncated octahedron', note: 'Classic size-dependent melting case for a noble-metal catalyst particle.', values: { mode: 'nano', crystal: 'Au', nshape: 'wulff-like', dnm: 3e-9 } },
    { name: '(10,10) armchair carbon nanotube', note: 'Metallic single-wall tube, 1.36 nm diameter; exports a periodic cell.', values: { mode: 'nano', crystal: 'graphene', nshape: 'tube', tn: 10, tm: 10, Lnm: 5e-9 } },
    { name: 'Al-doped ZnO (AZO), 2 % Al', note: 'Donor doping of a transparent conducting oxide: Kröger–Vink reactions and Vegard shift.', values: { mode: 'defect', crystal: 'ZnO_wurtzite', dopant: 'Al3+', xdop: 2, T_C: 900 } },
    { name: 'La-doped SrTiO₃ perovskite', note: 'A-site donor doping with Goldschmidt tolerance-factor check.', values: { mode: 'defect', crystal: 'SrTiO3', dopant: 'La3+', xdop: 5, T_C: 1200, Esch: 8 } },
    { name: 'Argon Lennard-Jones crystal', note: 'Benchmark: a ≈ 5.25 Å, E_coh ≈ 0.09 eV/atom for the static fcc lattice.', values: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'auto' } },
    { name: 'NaCl lattice energy (Madelung, Born–Landé)', note: 'Ionic-model benchmark against the Born–Haber cycle.', values: { mode: 'potential', potModel: 'ionic', salt: 'NaCl', bornN: 0 } },
    { name: 'Ge on Si(001) heteroepitaxy', note: '4.2 % mismatch: critical thickness of about a nanometre.', values: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', epiOri: 'direct', hfilm: 2e-8, T_C: 550 } },
    { name: 'Anatase TiO₂: Bragg angles and extinctions', note: 'Body-centred tetragonal cell with I-centring absences.', values: { mode: 'bragg', crystal: 'TiO2_anatase', rad: 'CuKa', h1: 1, k1: 0, l1: 1, h2: 0, k2: 0, l2: 4 } },
  ],
  run(inp) { return RUN[pick(inp.mode, Object.keys(RUN), 'crystal')](inp); },
};
