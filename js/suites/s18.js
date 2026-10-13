// Suite 18 — Multiscale modelling and homogenisation.
// Analytical micromechanics bounds/estimates, numerical periodic RVE homogenisation (finite volumes + Q4 finite elements, CG),
// laminate (Backus) averaging, lattice-to-continuum scaling laws and single-crystal → polycrystal (Voigt–Reuss–Hill) transfer.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => note ? { kind, label, text, note } : { kind, label, text };
const rec = (priority, text, basis, sustainability) => sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis };
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const line = (name, x, y) => { const X = [], Y = [], st = Math.max(1, Math.ceil(x.length / 400)); for (let i = 0; i < x.length; i += st) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) { X.push(x[i]); Y.push(y[i]); } if (!X.length) { X.push(0); Y.push(0); } return { name, x: X, y: Y }; };
const GPa = v => F(v / 1e9, 4);

const MODES = [
  { value: 'bounds', label: 'Two-phase composite: bounds and mean-field estimates' },
  { value: 'rve', label: 'Numerical RVE homogenisation (periodic cell problems)' },
  { value: 'laminate', label: 'Nanolaminate / multilayer homogenisation' },
  { value: 'lattice', label: 'Lattice-to-continuum effective constants' },
  { value: 'atomistic', label: 'Single-crystal constants → polycrystal (Voigt–Reuss–Hill)' },
];
const CONST = ['bounds', 'rve', 'laminate'], CL = ['bounds', 'rve', 'laminate', 'lattice'];
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'bounds', options: MODES, group: 'Model', help: 'Selects the scale-bridging method. Each analysis uses the inputs shown for it.' },
  { key: 'E_m', label: 'Matrix / phase A / strut Young\'s modulus', unit: 'Pa', value: 3e9, min: 1e6, max: 1.2e12, log: true, group: 'Phase A (matrix)', bus: 'prop.E', showIf: { mode: CL }, help: 'Young\'s modulus of the continuous phase, of laminate layer A, or of the solid strut material of a lattice.' },
  { key: 'nu_m', label: 'Matrix Poisson\'s ratio', value: 0.35, min: -0.9, max: 0.49, step: 0.01, group: 'Phase A (matrix)', bus: 'prop.nu', showIf: { mode: CL }, help: 'Poisson\'s ratio of phase A. Negative (auxetic) constituents are admissible down to −1.' },
  { key: 'k_m', label: 'Matrix thermal conductivity', unit: 'W/m/K', value: 0.2, min: 0.01, max: 5000, log: true, group: 'Phase A (matrix)', bus: 'prop.k_thermal', showIf: { mode: CONST }, help: 'Thermal conductivity of phase A. The same formulas apply to permittivity, permeability and diffusivity.' },
  { key: 'sig_m', label: 'Matrix electrical conductivity', unit: 'S/m', value: 1e-12, min: 1e-16, max: 1e8, log: true, group: 'Phase A (matrix)', bus: 'prop.sigma_elec', showIf: { mode: ['bounds', 'laminate'] }, help: 'Electrical conductivity of phase A; with a conducting filler it governs the percolation transition.' },
  { key: 'cte_m', label: 'Matrix thermal expansion', unit: '1/K', value: 6e-5, min: -5e-5, max: 3e-4, step: 1e-6, group: 'Phase A (matrix)', bus: 'prop.cte', showIf: { mode: CONST }, help: 'Linear coefficient of thermal expansion of phase A.' },
  { key: 'rho_m', label: 'Matrix / strut density', unit: 'kg/m³', value: 1200, min: 10, max: 23000, group: 'Phase A (matrix)', bus: 'prop.density', showIf: { mode: CL }, help: 'Mass density of phase A or of the lattice strut material.' },
  { key: 'E_i', label: 'Inclusion / phase B Young\'s modulus', unit: 'Pa', value: 1.4e11, min: 1e6, max: 1.2e12, log: true, group: 'Phase B (inclusion)', showIf: { mode: CONST }, help: 'Young\'s modulus of the dispersed phase or of laminate layer B.' },
  { key: 'nu_i', label: 'Inclusion Poisson\'s ratio', value: 0.35, min: -0.9, max: 0.49, step: 0.01, group: 'Phase B (inclusion)', showIf: { mode: CONST }, help: 'Poisson\'s ratio of phase B.' },
  { key: 'k_i', label: 'Inclusion thermal conductivity', unit: 'W/m/K', value: 50, min: 0.01, max: 5000, log: true, group: 'Phase B (inclusion)', showIf: { mode: CONST }, help: 'Thermal conductivity of phase B.' },
  { key: 'sig_i', label: 'Inclusion electrical conductivity', unit: 'S/m', value: 1e-2, min: 1e-16, max: 1e8, log: true, group: 'Phase B (inclusion)', showIf: { mode: ['bounds', 'laminate'] }, help: 'Electrical conductivity of phase B.' },
  { key: 'cte_i', label: 'Inclusion thermal expansion', unit: '1/K', value: 4.3e-6, min: -5e-5, max: 3e-4, step: 1e-6, group: 'Phase B (inclusion)', showIf: { mode: CONST }, help: 'Linear thermal expansion of phase B. Negative values (e.g. ZrW₂O₈ ≈ −8.7×10⁻⁶ 1/K) allow zero or negative composite expansion.' },
  { key: 'rho_i', label: 'Inclusion density', unit: 'kg/m³', value: 5606, min: 10, max: 23000, group: 'Phase B (inclusion)', showIf: { mode: CONST }, help: 'Mass density of phase B.' },
  { key: 'vf', label: 'Volume fraction of phase B', value: 0.2, min: 0.005, max: 0.95, step: 0.005, group: 'Microstructure', showIf: { mode: CONST }, help: 'Volume fraction of inclusions, or thickness fraction of layer B in a laminate.' },
  { key: 'shape', label: 'Inclusion shape', type: 'select', value: 'sphere', options: [{ value: 'sphere', label: 'Spheres / equiaxed particles' }, { value: 'fibre', label: 'Fibres, nanowires, nanotubes' }, { value: 'platelet', label: 'Platelets, nanosheets' }], group: 'Microstructure', showIf: { mode: ['bounds'] }, help: 'Shape class for the Halpin–Tsai stiffness estimate and for the percolation threshold.' },
  { key: 'est', label: 'Particulate stiffness estimate', type: 'select', value: 'dem', options: [{ value: 'dem', label: 'Differential scheme (default; valid to high loading)' }, { value: 'gsc', label: 'Generalised self-consistent (Christensen–Lo)' }, { value: 'lielens', label: 'Lielens interpolation' }, { value: 'mt', label: 'Mori–Tanaka (dilute to moderate loading)' }, { value: 'sc', label: 'Self-consistent (interpenetrating phases)' }], group: 'Microstructure', showIf: { mode: ['bounds'], shape: ['sphere'] }, help: 'Model for the reported modulus of a sphere-filled composite. All agree below about 15 vol %. Mori–Tanaka equals the Hashin–Shtrikman lower bound for stiff particles and is 15–30 % low at 40–50 vol %; the differential scheme follows measured particulate-composite moduli to about 40 vol %; the Lielens interpolation is an empirical blend that stays close to data up to 50 vol %.' },
  { key: 'aspect', label: 'Aspect ratio', value: 10, min: 1, max: 1e4, log: true, group: 'Microstructure', bus: 'geom.aspect', showIf: { mode: ['bounds'] }, help: 'Length/diameter of fibres or width/thickness of platelets. High aspect ratios stiffen efficiently and percolate at low loading.' },
  { key: 'orient', label: 'Orientation', type: 'select', value: 'random3d', options: [{ value: 'aligned', label: 'Aligned (longitudinal property)' }, { value: 'random2d', label: 'Random in plane' }, { value: 'random3d', label: 'Random in 3-D' }], group: 'Microstructure', showIf: { mode: ['bounds'] }, help: 'Orientation distribution used to combine longitudinal and transverse Halpin–Tsai moduli.' },
  { key: 'poros', label: 'Porosity of the matrix', value: 0.1, min: 0, max: 0.9, step: 0.01, group: 'Microstructure', bus: 'geom.porosity', showIf: { mode: ['bounds'] }, help: 'Pore volume fraction for the porous-solid estimates (Mori–Tanaka voids, differential scheme, Gibson–Ashby).' },
  // rve
  { key: 'micro', label: 'Microstructure generator', type: 'select', value: 'discs', options: [{ value: 'discs', label: 'Random discs (fibre cross-sections / particles)' }, { value: 'voronoi', label: 'Random two-phase Voronoi grains' }, { value: 'checker', label: 'Checkerboard (50:50 benchmark)' }], group: 'RVE', showIf: { mode: ['rve'] }, help: 'Seeded periodic two-phase 2-D microstructure. The checkerboard has the exact conductivity √(k₁k₂) and serves as a benchmark.' },
  { key: 'nres', label: 'RVE resolution', type: 'select', value: '48', options: [{ value: '32', label: '32 × 32' }, { value: '48', label: '48 × 48' }, { value: '64', label: '64 × 64' }], group: 'RVE', showIf: { mode: ['rve'] }, help: 'Grid cells per side for the conductivity cell problem. The elastic problem is solved on at most 32 × 32 elements.' },
  { key: 'n_incl', label: 'Number of inclusions / grains', value: 24, min: 4, max: 80, step: 1, group: 'RVE', showIf: { mode: ['rve'] }, help: 'Number of discs or Voronoi grains in the cell. More features make the cell more representative but resolve each one more coarsely.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'RVE', showIf: { mode: ['rve'] }, help: 'Seed of the microstructure realisation.' },
  // laminate
  { key: 'period_nm', label: 'Bilayer period', unit: 'nm', value: 20, min: 0.5, max: 1e5, log: true, group: 'Laminate', showIf: { mode: ['laminate'] }, help: 'Thickness of one A + B bilayer. Short periods increase the number of interfaces and their thermal resistance.' },
  { key: 'n_bi', label: 'Number of bilayers', value: 20, min: 1, max: 2000, step: 1, group: 'Laminate', showIf: { mode: ['laminate'] }, help: 'Number of repeated bilayers in the stack; few layers break scale separation in bending.' },
  { key: 'R_int', label: 'Interface thermal resistance', unit: 'm²K/W', value: 2e-9, min: 0, max: 1e-6, group: 'Laminate', showIf: { mode: ['laminate'] }, help: 'Kapitza (thermal boundary) resistance per interface, typically 10⁻⁹–10⁻⁸ m²K/W for dissimilar solids.' },
  { key: 'dT', label: 'Temperature change', unit: 'K', value: 100, min: -500, max: 1500, group: 'Laminate', showIf: { mode: ['laminate'] }, help: 'Uniform temperature excursion from the stress-free state for thermal-mismatch stresses.' },
  // lattice
  { key: 'topo', label: 'Lattice topology', type: 'select', value: 'octet', options: [{ value: 'octet', label: 'Octet truss (stretch-dominated)' }, { value: 'sc', label: 'Simple cubic' }, { value: 'bcc', label: 'Body-centred cubic' }, { value: 'kelvin', label: 'Kelvin cell (bending-dominated)' }, { value: 'reentrant', label: 'Re-entrant 2-D honeycomb (auxetic)' }], group: 'Lattice', showIf: { mode: ['lattice'] }, help: 'Strut architecture. It decides whether stiffness scales linearly (stretching) or quadratically (bending) with relative density.' },
  { key: 'rel', label: 'Relative density', value: 0.1, min: 0.002, max: 0.5, log: true, group: 'Lattice', bus: 'geom.rel_density', showIf: { mode: ['lattice'] }, help: 'Solid volume fraction of the lattice. Slender-strut scaling laws hold below about 0.3.' },
  { key: 'sy_s', label: 'Strut material yield strength', unit: 'Pa', value: 7e7, min: 1e5, max: 2e10, log: true, group: 'Lattice', bus: 'prop.yield', showIf: { mode: ['lattice'] }, help: 'Yield or fracture strength of the solid strut material.' },
  { key: 'cell', label: 'Unit-cell size', unit: 'm', value: 2e-6, min: 1e-8, max: 0.1, log: true, group: 'Lattice', bus: 'geom.cell_size', showIf: { mode: ['lattice'] }, help: 'Lattice period. It sets the strut dimensions and the frequency limit of the effective-medium description.' },
  { key: 'L_struct', label: 'Smallest structural dimension', unit: 'm', value: 1e-4, min: 1e-7, max: 10, log: true, group: 'Lattice', showIf: { mode: ['lattice'] }, help: 'Smallest dimension of the component made of the lattice; at least about ten cells are needed for bulk effective properties.' },
  { key: 'f_dyn', label: 'Operating frequency', unit: 'Hz', value: 1e6, min: 1, max: 1e11, log: true, group: 'Lattice', showIf: { mode: ['lattice'] }, help: 'Highest frequency of interest for dynamic use; the wavelength must stay long compared with the cell.' },
  { key: 'theta', label: 'Re-entrant angle θ', unit: '°', value: -20, min: -60, max: 45, step: 1, group: 'Lattice', showIf: { mode: ['lattice'] }, help: 'Inclination of the oblique cell walls (Gibson–Ashby convention): negative angles give a re-entrant, auxetic cell; +30° is the regular honeycomb.' },
  { key: 'h_l', label: 'Wall length ratio h/l', value: 2, min: 0.5, max: 4, step: 0.1, group: 'Lattice', showIf: { mode: ['lattice'] }, help: 'Ratio of vertical to oblique wall length of the 2-D honeycomb cell.' },
  // atomistic
  { key: 'src', label: 'Source of elastic constants', type: 'select', value: 'cij', options: [{ value: 'cij', label: 'Given single-crystal constants Cij' }, { value: 'lj', label: 'Lennard-Jones pair potential (fcc lattice sums)' }], group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Either enter measured/DFT elastic constants or derive them from a pair potential by static lattice sums.' },
  { key: 'sym', label: 'Crystal symmetry', type: 'select', value: 'cubic', options: [{ value: 'cubic', label: 'Cubic (C11, C12, C44)' }, { value: 'hex', label: 'Hexagonal (C11, C12, C13, C33, C44)' }, { value: 'tet', label: 'Tetragonal (adds C66)' }], group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Symmetry class deciding which independent constants are used.' },
  { key: 'C11', label: 'C11', unit: 'GPa', value: 165.7, min: 1, max: 1200, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Longitudinal stiffness along a crystal axis.' },
  { key: 'C12', label: 'C12', unit: 'GPa', value: 63.9, min: -200, max: 800, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Transverse coupling stiffness.' },
  { key: 'C44', label: 'C44', unit: 'GPa', value: 79.6, min: 0.5, max: 700, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Shear stiffness on the cube / prism planes.' },
  { key: 'C13', label: 'C13', unit: 'GPa', value: 105.1, min: -200, max: 800, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Coupling between basal plane and c axis (hexagonal and tetragonal only).' },
  { key: 'C33', label: 'C33', unit: 'GPa', value: 210.9, min: 1, max: 1500, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Stiffness along the c axis (hexagonal and tetragonal only).' },
  { key: 'C66', label: 'C66', unit: 'GPa', value: 192, min: 0.5, max: 700, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'In-plane shear stiffness (tetragonal only; hexagonal uses (C11 − C12)/2).' },
  { key: 'rho_x', label: 'Crystal density', unit: 'kg/m³', value: 2329, min: 100, max: 23000, group: 'Crystal', bus: 'prop.density', showIf: { mode: ['atomistic'] }, help: 'Mass density for wave speeds and Debye temperature (ignored for the pair-potential source, which computes it).' },
  { key: 'formula', label: 'Chemical formula', type: 'text', value: 'Si', group: 'Crystal', bus: 'material.formula', showIf: { mode: ['atomistic'] }, help: 'Formula unit, used for molar mass and atoms per formula in the Debye temperature.' },
  { key: 'eps_lj', label: 'Lennard-Jones well depth ε', unit: 'eV', value: 0.0104, min: 1e-4, max: 2, log: true, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Pair-potential well depth (argon ≈ 0.0104 eV).' },
  { key: 'sig_lj', label: 'Lennard-Jones length σ', unit: 'Å', value: 3.4, min: 1.5, max: 6, step: 0.01, group: 'Crystal', showIf: { mode: ['atomistic'] }, help: 'Pair-potential zero-crossing distance (argon ≈ 3.40 Å).' },
];
function read(inp) {
  const p = {};
  for (const d of inputs) {
    let v = inp ? inp[d.key] : undefined;
    if (d.type === 'select') p[d.key] = d.options.some(o => o.value === v) ? v : d.value;
    else if (d.type === 'text') p[d.key] = typeof v === 'string' && v.trim() ? v.trim() : d.value;
    else { v = Number(v); p[d.key] = N.clamp(Number.isFinite(v) ? v : d.value, d.min, d.max); }
  }
  return p;
}
function finish(r) {
  r.warnings = r.warnings || [];
  for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; r.warnings.push(`KPI "${k.label}" could not be evaluated for these inputs and is reported as 0.`); }
  for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
  if (!Object.keys(r.publish).length) r.publish['eff.status'] = 'not evaluated';
  return r;
}

// ---------- isotropic elasticity helpers ----------
const KG = (E, nu) => ({ K: E / (3 * (1 - 2 * nu)), G: E / (2 * (1 + nu)) });
const Eof = (K, G) => 9 * K * G / (3 * K + G), nuof = (K, G) => (3 * K - 2 * G) / (2 * (3 * K + G));
const zeta = (K, G) => G * (9 * K + 8 * G) / (6 * (K + 2 * G));
const Kb = (f, K1, K2, G0) => 1 / ((1 - f) / (K1 + 4 * G0 / 3) + f / (K2 + 4 * G0 / 3)) - 4 * G0 / 3;
const Gb = (f, G1, G2, z0) => 1 / ((1 - f) / (G1 + z0) + f / (G2 + z0)) - z0;
function mech(f, m, i) { // m, i: {K,G}
  const Gmin = Math.min(m.G, i.G), Gmax = Math.max(m.G, i.G), Kmin = Math.min(m.K, i.K), Kmax = Math.max(m.K, i.K);
  const hsKlo = Kb(f, m.K, i.K, Gmin), hsKhi = Kb(f, m.K, i.K, Gmax), hsGlo = Gb(f, m.G, i.G, zeta(Kmin, Gmin)), hsGhi = Gb(f, m.G, i.G, zeta(Kmax, Gmax));
  const KV = (1 - f) * m.K + f * i.K, GV = (1 - f) * m.G + f * i.G, KR = 1 / ((1 - f) / m.K + f / i.K), GR = 1 / ((1 - f) / m.G + f / i.G);
  const mtK = Kb(f, m.K, i.K, m.G), mtG = Gb(f, m.G, i.G, zeta(m.K, m.G));
  let K = 0.5 * (hsKlo + hsKhi), G = 0.5 * (hsGlo + hsGhi), it = 0, d = 1;
  for (; it < 300 && d > 1e-12; it++) { const Kn = Kb(f, m.K, i.K, G), Gn = Gb(f, m.G, i.G, zeta(K, G)); d = Math.abs(Kn - K) / K + Math.abs(Gn - G) / G; K = 0.5 * (K + Kn); G = 0.5 * (G + Gn); }
  return { hsKlo, hsKhi, hsGlo, hsGhi, KV, GV, KR, GR, mtK, mtG, scK: K, scG: G, scIter: it, scRes: d };
}
function dem(f, m, i) { // differential effective medium, spherical inclusions added to the current effective medium
  const n = 120, r = N.rk4((x, y) => { const K = y[0], G = y[1], z = zeta(K, G), P = (K + 4 * G / 3) / (i.K + 4 * G / 3), Q = (G + z) / (i.G + z); return [(i.K - K) * P / (1 - x), (i.G - G) * Q / (1 - x)]; }, [m.K, m.G], 0, f, n);
  const y = r.y[r.y.length - 1]; return { K: y[0], G: y[1] };
}
// Generalised self-consistent (three-phase) model of Christensen & Lo for spheres: composite-sphere bulk modulus (= Mori–Tanaka K) and the
// positive root of A x² + B x + C = 0 for x = G/G_m (coefficients in the corrected form of the published erratum).
function gscKG(c, m, i) {
  const K = Kb(c, m.K, i.K, m.G), nm = nuof(m.K, m.G), ni = nuof(i.K, i.G), r = i.G / m.G;
  const e1 = (r - 1) * (49 - 50 * ni * nm) + 35 * r * (ni - 2 * nm) + 35 * (2 * ni - nm), e2 = 5 * ni * (r - 8) + 7 * (r + 4), e3 = r * (8 - 10 * nm) + (7 - 5 * nm);
  const c103 = Math.pow(c, 10 / 3), c73 = Math.pow(c, 7 / 3), c53 = Math.pow(c, 5 / 3), q = 63 * (r - 1) * e2 + 2 * e1 * e3;
  const A = 8 * (r - 1) * (4 - 5 * nm) * e1 * c103 - 2 * q * c73 + 252 * (r - 1) * e2 * c53 - 50 * (r - 1) * (7 - 12 * nm + 8 * nm * nm) * e2 * c + 4 * (7 - 10 * nm) * e2 * e3;
  const B = -4 * (r - 1) * (1 - 5 * nm) * e1 * c103 + 4 * q * c73 - 504 * (r - 1) * e2 * c53 + 150 * (r - 1) * (3 - nm) * nm * e2 * c + 3 * (15 * nm - 7) * e2 * e3;
  const C = 4 * (r - 1) * (5 * nm - 7) * e1 * c103 - 2 * q * c73 + 252 * (r - 1) * e2 * c53 + 25 * (r - 1) * (nm * nm - 7) * e2 * c - (7 + 5 * nm) * e2 * e3;
  const disc = B * B - 4 * A * C, roots = disc >= 0 && A !== 0 ? [(-B + Math.sqrt(disc)) / (2 * A), (-B - Math.sqrt(disc)) / (2 * A)].filter(x => x > 0 && Number.isFinite(x)) : [];
  const mtG = Gb(c, m.G, i.G, zeta(m.K, m.G)), x = roots.length ? roots.reduce((a, b) => (Math.abs(b * m.G - mtG) < Math.abs(a * m.G - mtG) ? b : a)) : NaN;
  return { K, G: Number.isFinite(x) ? x * m.G : mtG, ok: Number.isFinite(x) };
}
// Lielens interpolation: compliance-weighted mean of the Mori–Tanaka estimate and its phase-inverted counterpart with weight (f + f²)/2
function lielensKG(f, m, i) { const w = 0.5 * (f + f * f), K1 = Kb(f, m.K, i.K, m.G), G1 = Gb(f, m.G, i.G, zeta(m.K, m.G)), K2 = Kb(f, m.K, i.K, i.G), G2 = Gb(f, m.G, i.G, zeta(i.K, i.G)); return { K: 1 / ((1 - w) / K1 + w / K2), G: 1 / ((1 - w) / G1 + w / G2) }; }
const EST = { dem: 'differential scheme', gsc: 'generalised self-consistent (Christensen–Lo) model', lielens: 'Lielens interpolation', mt: 'Mori–Tanaka estimate', sc: 'self-consistent scheme' };
function estimate(est, f, m, i, M) { M = M || mech(f, m, i); return est === 'mt' ? { K: M.mtK, G: M.mtG } : est === 'sc' ? { K: M.scK, G: M.scG } : est === 'gsc' ? gscKG(f, m, i) : est === 'lielens' ? lielensKG(f, m, i) : dem(f, m, i); }
// transport (conductivity, permittivity, diffusivity): dimension d
const mg = (km, ki, f, d = 3) => Math.abs(ki - km) < 1e-300 ? km : km + f / (1 / (ki - km) + (1 - f) / (d * km));
const brug = (km, ki, f) => { const b = (3 * f - 1) * ki + (2 - 3 * f) * km; return 0.25 * (b + Math.sqrt(b * b + 8 * ki * km)); };
const halpin = (Em, Ei, f, z) => { const eta = (Ei / Em - 1) / (Ei / Em + z); return Em * (1 + z * eta * f) / (1 - eta * f); };
const levin = (f, am, ai, Km, Ki, Ks) => Math.abs(1 / Ki - 1 / Km) < 1e-9 / Km ? (1 - f) * am + f * ai : (1 - f) * am + f * ai + (ai - am) / (1 / Ki - 1 / Km) * (1 / Ks - ((1 - f) / Km + f / Ki));

function bounds(p) {
  const f = p.vf, m = KG(p.E_m, p.nu_m), i = KG(p.E_i, p.nu_i), M = mech(f, m, i), warnings = [];
  const Ehs = [Eof(M.hsKlo, M.hsGlo), Eof(M.hsKhi, M.hsGhi)], EV = Eof(M.KV, M.GV), ER = Eof(M.KR, M.GR), Emt = Eof(M.mtK, M.mtG), Esc = Eof(M.scK, M.scG), D = dem(f, m, i), Edem = Eof(D.K, D.G);
  const AR = p.aspect, zL = p.shape === 'sphere' ? 2 : 2 * AR, zT = 2;
  const htAt = ff => { const EL = halpin(p.E_m, p.E_i, ff, zL), ET = halpin(p.E_m, p.E_i, ff, zT); if (p.shape === 'sphere') return ET; if (p.orient === 'aligned') return EL; if (p.orient === 'random2d') return 0.375 * EL + 0.625 * ET; return p.shape === 'platelet' ? 0.49 * EL + 0.51 * ET : 0.2 * EL + 0.8 * ET; };
  // Particulate (sphere) composites: the design estimate is selectable; the default is the differential scheme, which stays close to measured
  // moduli of glass-bead/epoxy up to about 40 vol % (Smith, J. Res. NBS 80A, 1976: −1, −1, −3 % at 15, 30, 40 vol % and −14 % at 50 vol %),
  // where the Mori–Tanaka estimate (= Hashin–Shtrikman lower bound) is 3, 8, 15 and 30 % low.
  const Eht = htAt(f), useHT = p.shape !== 'sphere', est = useHT ? 'mt' : p.est, S = estimate(est, f, m, i, M), G0 = gscKG(f, m, i), L0 = lielensKG(f, m, i), Egsc = Eof(G0.K, G0.G), Eli = Eof(L0.K, L0.G);
  if (est === 'gsc' && !G0.ok) warnings.push('The quadratic of the generalised self-consistent model has no positive root for this combination; the Mori–Tanaka shear modulus was used instead.');
  const Eeff = useHT ? Eht : Eof(S.K, S.G), nuEff = nuof(S.K, S.G), Geff = useHT ? Eeff / (2 * (1 + nuEff)) : S.G, Keff = useHT ? Eeff / (3 * (1 - 2 * nuEff)) : S.K, estName = useHT ? 'Halpin–Tsai' : EST[est];
  // transport
  const lo = Math.min(p.k_m, p.k_i), hi = Math.max(p.k_m, p.k_i), flo = p.k_m <= p.k_i ? f : 1 - f; // fraction of the better conductor = phase "hi"
  const khsLo = mg(lo, hi, flo), khsHi = mg(hi, lo, 1 - flo), kmg = mg(p.k_m, p.k_i, f), kbr = brug(p.k_m, p.k_i, f), kW = [1 / ((1 - f) / p.k_m + f / p.k_i), (1 - f) * p.k_m + f * p.k_i];
  const fc = p.shape === 'sphere' ? 1 / 3 : Math.min(1 / 3, (p.shape === 'fibre' ? 0.7 : 1.15) / AR);
  const sLo = Math.min(p.sig_m, p.sig_i), sHi = Math.max(p.sig_m, p.sig_i), sW = [1 / ((1 - f) / p.sig_m + f / p.sig_i), (1 - f) * p.sig_m + f * p.sig_i];
  const percAt = ff => { if (p.sig_i <= p.sig_m) return mg(p.sig_m, p.sig_i, ff); const above = ff > fc ? p.sig_i * Math.pow((ff - fc) / (1 - fc), 2) : 0, below = ff < fc ? p.sig_m * Math.pow(fc / (fc - ff), 0.73) : Infinity; const lw = 1 / ((1 - ff) / p.sig_m + ff / p.sig_i), up = (1 - ff) * p.sig_m + ff * p.sig_i; return N.clamp(ff > fc ? Math.max(above, p.sig_m) : below, lw, up); };
  const sigEff = percAt(f), percolated = p.sig_i > p.sig_m && f > fc;
  // thermal expansion
  const cteAt = (ff, which) => { const Mm = mech(ff, m, i); const Ks = which === 'lo' ? Mm.hsKlo : which === 'hi' ? Mm.hsKhi : estimate(est, ff, m, i, Mm).K; return levin(ff, p.cte_m, p.cte_i, m.K, i.K, Ks); };
  const turner = ff => ((1 - ff) * p.cte_m * m.K + ff * p.cte_i * i.K) / ((1 - ff) * m.K + ff * i.K);
  const kerner = ff => (1 - ff) * p.cte_m + ff * p.cte_i + ff * (1 - ff) * (p.cte_i - p.cte_m) * (i.K - m.K) / ((1 - ff) * m.K + ff * i.K + 3 * i.K * m.K / (4 * m.G));
  const cteMT = cteAt(f), cteB = [cteAt(f, 'lo'), cteAt(f, 'hi')].sort((a, b) => a - b), cteT = turner(f), cteK = kerner(f);
  let fZero = 0; if (p.cte_m * p.cte_i < 0) { try { fZero = N.brent(ff => cteAt(ff), 1e-4, 0.9999); } catch { fZero = 0; } }
  const rho = (1 - f) * p.rho_m + f * p.rho_i;
  // porous matrix
  const po = p.poros, void0 = { K: m.K * 1e-9, G: m.G * 1e-9 }, Mp = mech(po, m, void0), Dp = po > 0 ? dem(Math.min(po, 0.95), m, void0) : m, EporMT = Eof(Mp.mtK, Mp.mtG), EporDEM = Eof(Dp.K, Dp.G), EporGA = p.E_m * (1 - po) ** 2, kPor = p.k_m * 2 * (1 - po) / (2 + po);
  if (!useHT && est === 'mt' && f > 0.25 && i.G > 3 * m.G) warnings.push(`Mori–Tanaka coincides with the Hashin–Shtrikman lower bound for stiff spheres and falls increasingly below measured moduli as the loading rises (about 8 % at 30 vol %, 15 % at 40 vol % and 30 % at 50 vol % for glass beads in epoxy); the differential scheme gives ${GPa(Edem)} GPa here.`);
  if (!useHT && f > 0.45 && i.G > 3 * m.G) warnings.push(`At ${F(100 * f, 3)} vol % stiff particles approach random close packing (about 64 vol % for equal spheres): particle contacts stiffen the composite beyond what any matrix-continuous mean-field estimate gives. Measured glass-bead/epoxy moduli at 50 vol % lie 14 % above the differential scheme and 6 % above the Lielens interpolation (${GPa(Eli)} GPa here).`);
  if (!useHT && est === 'sc' && f > 0.3 && f < 0.7) warnings.push('The symmetric self-consistent scheme treats the filler as partly contiguous and over-predicts the modulus of dispersed stiff particles above about 30 vol % (by 14 % at 30 vol % and 26 % at 40 vol % for glass beads in epoxy).');
  if (useHT && f * AR > 10 && p.orient !== 'aligned') warnings.push('High loading of high-aspect-ratio filler: randomly oriented fibres/platelets cannot pack beyond roughly 5/AR without alignment or agglomeration, so the dispersion assumed here may not be achievable.');
  if (Math.max(p.nu_m, p.nu_i) > 0.48) warnings.push('A nearly incompressible phase makes bulk-modulus-based estimates (and thermal expansion) very sensitive to Poisson\'s ratio.');
  if (useHT) warnings.push('Halpin–Tsai is a semi-empirical interpolation; its orientation averaging weights are approximations and the result may fall outside the isotropic Hashin–Shtrikman bounds when the filler is aligned.');
  const fs = N.linspace(0, 0.95, 39), sweep = fs.map(ff => mech(Math.max(ff, 1e-9), m, i));
  const within = (v, a, b, tol = 1e-6) => v >= Math.min(a, b) * (1 - tol) - 1e-300 && v <= Math.max(a, b) * (1 + tol) + 1e-300;
  const negCte = cteMT < 0;
  return {
    summary: `Two-phase composite with ${F(100 * f, 3)} vol % ${p.shape} inclusions: effective E = ${GPa(Eeff)} GPa (${estName}; Hashin–Shtrikman bounds ${GPa(Ehs[0])}–${GPa(Ehs[1])} GPa), thermal conductivity ${F(kmg, 4)} W/m/K (Maxwell-Garnett; Bruggeman ${F(kbr, 4)}), electrical conductivity ${F(sigEff, 3)} S/m (${p.sig_i <= p.sig_m ? 'filler less conducting than the matrix, no percolation effect' : (percolated ? 'above' : 'below') + ' the percolation threshold of ' + F(100 * fc, 3) + ' vol %'}), and thermal expansion ${F(cteMT * 1e6, 4)}×10⁻⁶ 1/K${negCte ? ' — a negative-expansion composite' : fZero > 0 ? `; zero expansion would need ${F(100 * fZero, 3)} vol % filler` : ''}.`,
    kpis: [
      kpi('E_eff', 'Effective Young\'s modulus', Eeff, 'Pa', `${useHT ? `Halpin–Tsai estimate for ${p.shape}s of aspect ratio ${F(AR)} (${p.orient})` : `Estimate for spherical inclusions by the ${estName}`}; ${F(Eeff / p.E_m, 4)} times the matrix modulus.${useHT ? '' : ` Other estimates: Mori–Tanaka ${GPa(Emt)}, generalised self-consistent ${GPa(Egsc)}, differential ${GPa(Edem)}, Lielens ${GPa(Eli)}, self-consistent ${GPa(Esc)} GPa.`}`, null),
      kpi('E_hs_lo', 'Hashin–Shtrikman lower bound on E', Ehs[0], 'Pa', 'Tightest lower bound on the modulus of any isotropic mixture of the two phases at this volume fraction (from the K and G bounds in Walpole form).', null),
      kpi('E_hs_hi', 'Hashin–Shtrikman upper bound on E', Ehs[1], 'Pa', 'Tightest upper bound for an isotropic mixture; reached by a microstructure in which the stiff phase is the continuous one.', null),
      kpi('E_sc', 'Self-consistent estimate of E', Esc, 'Pa', `Self-consistent (effective-medium) scheme for spheres, solved by fixed-point iteration in ${M.scIter} steps; it treats both phases symmetrically and suits interpenetrating microstructures.`, null),
      kpi('nu_eff', 'Effective Poisson\'s ratio', nuEff, '–', `From the bulk and shear moduli of the ${useHT ? 'Mori–Tanaka estimate' : estName}, ν = (3K − 2G)/(2(3K + G)).`, nuEff > -1 && nuEff < 0.5),
      kpi('G_eff', 'Effective shear modulus', Geff, 'Pa', `Shear modulus consistent with the reported E and ν (${useHT ? 'derived from the Halpin–Tsai modulus' : estName}).`, null),
      kpi('K_eff', 'Effective bulk modulus', Keff, 'Pa', 'Bulk modulus consistent with the reported E and ν; it also enters the thermal-expansion estimate through Levin\'s relation.', null),
      kpi('k_eff', 'Effective thermal conductivity', kmg, 'W/m/K', `Maxwell-Garnett (Maxwell–Eucken) estimate for dispersed inclusions; the symmetric Bruggeman value is ${F(kbr, 4)} W/m/K and the Hashin–Shtrikman bounds are ${F(khsLo, 4)}–${F(khsHi, 4)} W/m/K.`, null),
      kpi('sigma_eff', 'Effective electrical conductivity', sigEff, 'S/m', p.sig_i <= p.sig_m ? 'Maxwell-Garnett estimate for a poorly conducting filler dispersed in the conducting matrix; no percolation transition occurs in this case.' : `Percolation scaling σ ∝ (f − f_c)^2 above and (f_c − f)^−0.73 below the threshold, clipped to the Wiener bounds. With f_c = ${F(fc, 3)} the composite is ${percolated ? 'conducting' : 'not percolated'}.`, null),
      kpi('f_c', 'Percolation threshold', fc, '–', p.shape === 'sphere' ? 'Bruggeman effective-medium threshold of 1/3 for equiaxed particles (random packings percolate between about 0.16 and 0.33).' : `Excluded-volume estimate ${p.shape === 'fibre' ? '0.7/AR for slender rods' : '1.15/AR for thin discs'}, valid for randomly oriented, well-dispersed filler.`, null),
      kpi('cte_eff', 'Effective thermal expansion', cteMT, '1/K', `Levin relation with the bulk modulus of the ${useHT ? 'Mori–Tanaka estimate' : estName}. Schapery bounds: ${F(cteB[0] * 1e6, 4)} to ${F(cteB[1] * 1e6, 4)} ×10⁻⁶ 1/K; Turner ${F(cteT * 1e6, 4)}, Kerner ${F(cteK * 1e6, 4)} ×10⁻⁶ 1/K.`, null),
      kpi('f_zero', 'Filler fraction for zero expansion', fZero, '–', fZero > 0 ? 'Volume fraction at which the Levin expansion (with the selected bulk-modulus estimate) changes sign, possible because the two phases expand in opposite senses.' : 'No zero crossing exists because both phases have thermal expansion of the same sign; reported as 0.', null),
      kpi('rho_eff', 'Effective density', rho, 'kg/m³', 'Rule of mixtures on volume fractions, exact for any microstructure without porosity.', null),
      kpi('E_porous', 'Porous-matrix modulus', EporDEM, 'Pa', `Differential-scheme estimate for the matrix with ${F(100 * po, 3)} % spherical pores; Mori–Tanaka gives ${GPa(EporMT)} GPa and the Gibson–Ashby open-cell law ${GPa(EporGA)} GPa.`, null),
    ],
    series: [
      { key: 'E_f', title: 'Young\'s modulus versus volume fraction', xlabel: 'volume fraction of phase B', ylabel: 'E (GPa)', logy: p.E_i / p.E_m > 30 || p.E_m / p.E_i > 30, lines: [line('Voigt', fs, sweep.map(s => Eof(s.KV, s.GV) / 1e9)), line('Reuss', fs, sweep.map(s => Eof(s.KR, s.GR) / 1e9)), line('HS upper', fs, sweep.map(s => Eof(s.hsKhi, s.hsGhi) / 1e9)), line('HS lower', fs, sweep.map(s => Eof(s.hsKlo, s.hsGlo) / 1e9)), line('Mori–Tanaka', fs, sweep.map(s => Eof(s.mtK, s.mtG) / 1e9)), line('Self-consistent', fs, sweep.map(s => Eof(s.scK, s.scG) / 1e9)), line('Differential scheme', fs, fs.map(ff => { const d = dem(Math.max(ff, 1e-9), m, i); return Eof(d.K, d.G) / 1e9; })), line('Generalised self-consistent', fs, fs.map(ff => { const d = gscKG(Math.max(ff, 1e-9), m, i); return Eof(d.K, d.G) / 1e9; })), line('Halpin–Tsai', fs, fs.map(ff => htAt(ff) / 1e9))], explain: 'For stiff particles the estimates fan out above about 25 vol %: Mori–Tanaka follows the lower bound, the generalised self-consistent and differential schemes lie progressively higher and the symmetric self-consistent scheme highest. All estimates must lie between the Voigt and Reuss limits; isotropic ones also between the Hashin–Shtrikman bounds. The gap between the bounds shows how much microstructure matters at each loading.' },
      { key: 'k_f', title: 'Conductivity versus volume fraction', xlabel: 'volume fraction of phase B', ylabel: 'thermal conductivity (W/m/K)', logy: hi / lo > 30, lines: [line('Maxwell-Garnett', fs, fs.map(ff => mg(p.k_m, p.k_i, ff))), line('Bruggeman', fs, fs.map(ff => brug(p.k_m, p.k_i, ff))), line('Wiener parallel', fs, fs.map(ff => (1 - ff) * p.k_m + ff * p.k_i)), line('Wiener series', fs, fs.map(ff => 1 / ((1 - ff) / p.k_m + ff / p.k_i)))], explain: 'Maxwell-Garnett keeps the filler isolated at all loadings, Bruggeman percolates at one third. The same curves describe permittivity and diffusivity.' },
      { key: 'cte_f', title: 'Thermal expansion versus volume fraction', xlabel: 'volume fraction of phase B', ylabel: 'CTE (10⁻⁶ 1/K)', lines: [line('Levin + selected estimate', fs, fs.map(ff => cteAt(Math.max(ff, 1e-9)) * 1e6)), line('Schapery (HS lower K)', fs, fs.map(ff => cteAt(Math.max(ff, 1e-9), 'lo') * 1e6)), line('Schapery (HS upper K)', fs, fs.map(ff => cteAt(Math.max(ff, 1e-9), 'hi') * 1e6)), line('Turner', fs, fs.map(ff => turner(ff) * 1e6)), line('Kerner', fs, fs.map(ff => kerner(ff) * 1e6))], explain: 'A zero crossing appears when a negative-expansion filler is combined with a normal matrix. Stiff fillers pull the curve below the rule of mixtures.' },
      { key: 'sig_f', title: 'Electrical conductivity versus volume fraction', xlabel: 'volume fraction of phase B', ylabel: 'σ (S/m)', logy: true, lines: [line('percolation scaling', fs, fs.map(ff => Math.max(percAt(Math.min(Math.max(ff, 1e-6), 0.999)), 1e-300)))], explain: 'Insulator–conductor transition at the percolation threshold. The exponents (t = 2, s = 0.73) are the universal 3-D values; real composites often show larger t.' },
    ],
    tables: [...(useHT ? [] : [{ key: 'est', title: 'Estimates of the particulate-composite modulus', columns: ['Model', 'E (GPa)', 'ν', 'E/E_matrix'], rows: [['Mori–Tanaka (Hashin–Shtrikman bound of the matrix-continuous kind)', M.mtK, M.mtG], ['Generalised self-consistent (Christensen–Lo)', G0.K, G0.G], ['Differential scheme', D.K, D.G], ['Lielens interpolation', L0.K, L0.G], ['Self-consistent (symmetric)', M.scK, M.scG]].map(([nm, K, G]) => [nm, GPa(Eof(K, G)), F(nuof(K, G), 4), F(Eof(K, G) / p.E_m, 4)]), explain: 'Mean-field estimates for spheres at this loading. They agree at low loading; their spread at high loading is the model uncertainty. Against measured glass-bead/epoxy moduli the differential scheme is within 3 % up to 40 vol %.' }]), { key: 'porous', title: 'Porous-solid estimates for the matrix', columns: ['Model', 'E (GPa)', 'E/E_solid'], rows: [['Mori–Tanaka, spherical voids', GPa(EporMT), F(EporMT / p.E_m, 4)], ['Differential scheme, spherical voids', GPa(EporDEM), F(EporDEM / p.E_m, 4)], ['Gibson–Ashby open cell (1 − p)²', GPa(EporGA), F(EporGA / p.E_m, 4)], ['Maxwell conductivity factor k/k_solid', '—', F(kPor / p.k_m, 4)]], explain: 'Stiffness of the matrix containing the stated porosity by three standard models; the spread indicates the model uncertainty for the pore morphology.' }],
    equations: [
      eqn('formula', 'Hashin–Shtrikman bounds (Walpole form)', 'K* = [Σ fᵢ/(Kᵢ + 4G₀/3)]⁻¹ − 4G₀/3;  G* = [Σ fᵢ/(Gᵢ + ζ₀)]⁻¹ − ζ₀;  ζ = G(9K + 8G)/(6(K + 2G))', 'G₀, ζ₀ from the softest (lower) or stiffest (upper) phase'),
      eqn('formula', 'Mori–Tanaka / self-consistent', 'same expressions with (G₀, ζ₀) of the matrix (Mori–Tanaka) or of the effective medium (self-consistent)'),
      eqn('formula', 'Halpin–Tsai', 'E/E_m = (1 + ζηf)/(1 − ηf),  η = (E_f/E_m − 1)/(E_f/E_m + ζ),  ζ_L = 2·AR, ζ_T = 2'),
      eqn('formula', 'Maxwell-Garnett and Bruggeman', 'k_MG = k_m + f/[1/(k_i − k_m) + (1 − f)/(3k_m)];  Σ fᵢ(kᵢ − k*)/(kᵢ + 2k*) = 0'),
      eqn('formula', 'Levin relation (thermal expansion)', 'α* = ⟨α⟩ + (α_i − α_m)/(1/K_i − 1/K_m)·(1/K* − ⟨1/K⟩)'),
      eqn('formula', 'Differential scheme', 'dK/df = (K_i − K)P/(1 − f),  dG/df = (G_i − G)Q/(1 − f),  P = (K + 4G/3)/(K_i + 4G/3),  Q = (G + ζ)/(G_i + ζ)', 'inclusions added in infinitesimal increments to the current effective medium, starting from the matrix'),
      eqn('formula', 'Generalised self-consistent model (spheres)', 'K* = K_m + f/[1/(K_i − K_m) + (1 − f)/(K_m + 4G_m/3)];  A(G*/G_m)² + B(G*/G_m) + C = 0', 'composite sphere embedded in the effective medium; A, B, C are the Christensen–Lo polynomials in f^(1/3), G_i/G_m, ν_m and ν_i'),
      eqn('formula', 'Lielens interpolation', '1/K* = (1 − w)/K_MT + w/K_MT,inverse,  w = (f + f²)/2  (same for G)', 'interpolates between Mori–Tanaka with the matrix and with the filler as the continuous phase'),
    ],
    checks: [
      chk('Hashin–Shtrikman bounds ordered within Voigt–Reuss', M.KR <= M.hsKlo * (1 + 1e-9) && M.hsKlo <= M.hsKhi * (1 + 1e-9) && M.hsKhi <= M.KV * (1 + 1e-9) && M.GR <= M.hsGlo * (1 + 1e-9) && M.hsGlo <= M.hsGhi * (1 + 1e-9) && M.hsGhi <= M.GV * (1 + 1e-9), `K: ${GPa(M.KR)} ≤ ${GPa(M.hsKlo)} ≤ ${GPa(M.hsKhi)} ≤ ${GPa(M.KV)} GPa`),
      chk('Mori–Tanaka and self-consistent within HS bounds', within(M.mtK, M.hsKlo, M.hsKhi) && within(M.mtG, M.hsGlo, M.hsGhi) && within(M.scK, M.hsKlo, M.hsKhi, 1e-4) && within(M.scG, M.hsGlo, M.hsGhi, 1e-4), `G: MT ${GPa(M.mtG)}, SC ${GPa(M.scG)} in [${GPa(M.hsGlo)}, ${GPa(M.hsGhi)}] GPa`),
      chk('Design modulus within Voigt–Reuss limits', within(Eeff, ER, EV, 1e-6), `E = ${GPa(Eeff)} GPa in [${GPa(ER)}, ${GPa(EV)}] GPa`),
      ...(useHT ? [] : [chk('Selected estimate within the Hashin–Shtrikman bounds', within(S.K, M.hsKlo, M.hsKhi, 1e-3) && within(S.G, M.hsGlo, M.hsGhi, 1e-3), `${estName}: K = ${GPa(S.K)}, G = ${GPa(S.G)} GPa; bounds K [${GPa(M.hsKlo)}, ${GPa(M.hsKhi)}], G [${GPa(M.hsGlo)}, ${GPa(M.hsGhi)}] GPa`)]),
      chk('Conductivity estimates within HS bounds', within(kmg, khsLo, khsHi) && within(kbr, khsLo, khsHi), `[${F(khsLo, 4)}, ${F(khsHi, 4)}] W/m/K`),
      chk('Self-consistent iteration converged', M.scRes < 1e-8, `residual ${M.scRes.toExponential(2)} after ${M.scIter} iterations`),
      chk('Thermodynamic admissibility of ν', nuEff > -1 && nuEff < 0.5, `ν = ${F(nuEff, 4)}`),
    ],
    warnings,
    recommendations: [
      rec('high', percolated ? `The filler network is percolated (f = ${F(f, 3)} > f_c = ${F(fc, 3)}): conductivity is controlled by contact resistance between particles, so improve dispersion and contacts rather than adding more filler.` : `For an electrically conducting composite the loading must exceed about ${F(100 * fc, 3)} vol %; a higher aspect ratio lowers this threshold in proportion to 1/AR.`, 'Excluded-volume percolation threshold.'),
      rec('high', `Measure the modulus of one composition and compare with the bounds ${GPa(Ehs[0])}–${GPa(Ehs[1])} GPa: a value near the lower bound indicates isolated particles, near the upper bound a contiguous stiff network; agglomeration and weak interfaces push results below the Mori–Tanaka line. The estimates span ${GPa(Math.min(Emt, Egsc, Edem, Eli, Esc))}–${GPa(Math.max(Emt, Egsc, Edem, Eli, Esc))} GPa at this loading.`, 'Position between Hashin–Shtrikman bounds identifies the microstructure class.'),
      rec('medium', fZero > 0 ? `Target ${F(100 * fZero, 3)} vol % filler for near-zero thermal expansion and verify by dilatometry; thermal-mismatch stresses between the phases scale with (α_m − α_i)ΔT and can cause microcracking on cycling.` : 'To tailor thermal expansion further, choose a filler with lower (or negative) expansion and higher bulk modulus: stiff fillers constrain the matrix more than the rule of mixtures suggests.', 'Levin relation and Schapery bounds.'),
      rec('medium', `Use high-aspect-ratio filler to reach the target property at lower loading (here ${F(100 * f, 3)} vol %): less filler means lower embodied energy and cost, lower viscosity during processing and easier recycling of the matrix.`, 'Halpin–Tsai and percolation scaling with aspect ratio.', true),
    ],
    publish: { 'eff.E': Eeff, 'eff.nu': nuEff, 'eff.G': Geff, 'eff.K': Keff, 'eff.k': kmg, 'eff.sigma': sigEff, 'eff.cte': cteMT, 'eff.density': rho },
    numerics: { method: 'Closed-form bounds and mean-field estimates (Mori–Tanaka, Christensen–Lo quadratic, Lielens); fixed-point self-consistent iteration; RK4 integration of the differential scheme', dof: 2, steps: 120, iterations: M.scIter, residual: M.scRes, notes: 'Isotropic phases, perfect interfaces, spherical-inclusion Eshelby factors for Mori–Tanaka, self-consistent and differential schemes.' },
  };
}

// ---------- microstructure generation (periodic unit cell) ----------
function microGen(type, n, nIncl, vf, seed) {
  const rnd = N.rng(seed), ph = new Uint8Array(n * n), per = d => { d = Math.abs(d); return d > 0.5 ? 1 - d : d; };
  if (type === 'checker') { const m = [2, 4, 8, 16].reduce((b, v) => n % v === 0 && Math.abs(v - Math.sqrt(nIncl)) < Math.abs(b - Math.sqrt(nIncl)) ? v : b, 2); for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) ph[i + n * j] = (Math.floor(i * m / n) + Math.floor(j * m / n)) % 2; }
  else if (type === 'voronoi') {
    const sx = [], sy = [], sp = []; for (let k = 0; k < nIncl; k++) { sx.push(rnd()); sy.push(rnd()); sp.push(0); }
    const nB = N.clamp(Math.round(vf * nIncl), 1, nIncl - 1); { const idx = sx.map((_, k) => k); for (let k = idx.length - 1; k > 0; k--) { const j = rnd.int(k + 1);[idx[k], idx[j]] = [idx[j], idx[k]]; } for (let k = 0; k < nB; k++) sp[idx[k]] = 1; }
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const x = (i + 0.5) / n, y = (j + 0.5) / n; let best = 9, bk = 0; for (let k = 0; k < nIncl; k++) { const d = per(x - sx[k]) ** 2 + per(y - sy[k]) ** 2; if (d < best) { best = d; bk = k; } } ph[i + n * j] = sp[bk]; }
  } else {
    const r = Math.sqrt(vf / (nIncl * N.PI)), cx = [], cy = [];
    for (let k = 0; k < nIncl; k++) { let x = rnd(), y = rnd(); for (let t = 0; t < 150; t++) { let ok = true; for (let q = 0; q < cx.length && ok; q++) if (per(x - cx[q]) ** 2 + per(y - cy[q]) ** 2 < (2.05 * r) ** 2) ok = false; if (ok) break; x = rnd(); y = rnd(); } cx.push(x); cy.push(y); }
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const x = (i + 0.5) / n, y = (j + 0.5) / n; for (let k = 0; k < nIncl; k++) if (per(x - cx[k]) ** 2 + per(y - cy[k]) ** 2 < r * r) { ph[i + n * j] = 1; break; } }
  }
  let s = 0; for (let i = 0; i < ph.length; i++) s += ph[i];
  return { ph, f: s / ph.length };
}
// periodic conductivity cell problem, finite volumes with harmonic face conductances
function solveCond(ph, n, k0, k1, maxIter, wantField, tol = 1e-6, nLoad = 2) {
  const nn = n * n, ke = new Float64Array(nn), kn = new Float64Array(nn), dg = new Float64Array(nn), E = new Int32Array(nn), W = new Int32Array(nn), No = new Int32Array(nn), S = new Int32Array(nn), kc = i => ph[i] ? k1 : k0;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const P = i + n * j; E[P] = (i + 1) % n + n * j; W[P] = (i + n - 1) % n + n * j; No[P] = i + n * ((j + 1) % n); S[P] = i + n * ((j + n - 1) % n); }
  for (let P = 0; P < nn; P++) { const a = kc(P), b = kc(E[P]), c = kc(No[P]); ke[P] = 2 * a * b / (a + b); kn[P] = 2 * a * c / (a + c); }
  for (let P = 0; P < nn; P++) dg[P] = ke[P] + ke[W[P]] + kn[P] + kn[S[P]];
  const y = new Float64Array(nn), mv = x => { for (let P = 0; P < nn; P++) y[P] = dg[P] * x[P] - ke[P] * x[E[P]] - ke[W[P]] * x[W[P]] - kn[P] * x[No[P]] - kn[S[P]] * x[S[P]]; return y; };
  const K = [[0, 0], [0, 0]], Wen = [0, 0]; let iters = 0, res = 0, conv = true, field = null;
  for (let a = 0; a < nLoad; a++) {
    const b = new Float64Array(nn); for (let P = 0; P < nn; P++) b[P] = a === 0 ? ke[P] - ke[W[P]] : kn[P] - kn[S[P]];
    const sol = N.cg(mv, b, { tol, maxIter, diag: dg }), x = sol.x; iters += sol.iterations; res = Math.max(res, sol.residual); conv = conv && sol.converged;
    let qx = 0, qy = 0, en = 0; for (let P = 0; P < nn; P++) { const gx = (a === 0 ? 1 : 0) + x[E[P]] - x[P], gy = (a === 1 ? 1 : 0) + x[No[P]] - x[P]; qx += ke[P] * gx; qy += kn[P] * gy; en += ke[P] * gx * gx + kn[P] * gy * gy; }
    K[0][a] = qx / nn; K[1][a] = qy / nn; Wen[a] = en / nn;
    if (a === 0 && wantField) { field = new Float64Array(nn); for (let P = 0; P < nn; P++) { const fx = 0.5 * (ke[P] * (1 + x[E[P]] - x[P]) + ke[W[P]] * (1 + x[P] - x[W[P]])), fy = 0.5 * (kn[P] * (x[No[P]] - x[P]) + kn[S[P]] * (x[P] - x[S[P]])); field[P] = Math.hypot(fx, fy); } }
  }
  return { K, Wen, iters, res, conv, field };
}
// Q4 plane-strain element on the unit square
function q4(lam, mu) {
  const C = [[lam + 2 * mu, lam, 0], [lam, lam + 2 * mu, 0], [0, 0, mu]], Ke = new Float64Array(64), Fe = [new Float64Array(8), new Float64Array(8), new Float64Array(8)], g = [0.5 - 0.5 / Math.sqrt(3), 0.5 + 0.5 / Math.sqrt(3)];
  const Bat = (xi, et) => { const dx = [-(1 - et), 1 - et, et, -et], dy = [-(1 - xi), -xi, xi, 1 - xi], B = [new Float64Array(8), new Float64Array(8), new Float64Array(8)]; for (let a = 0; a < 4; a++) { B[0][2 * a] = dx[a]; B[1][2 * a + 1] = dy[a]; B[2][2 * a] = dy[a]; B[2][2 * a + 1] = dx[a]; } return B; };
  for (const xi of g) for (const et of g) { const B = Bat(xi, et), CB = C.map(r => { const o = new Float64Array(8); for (let c = 0; c < 8; c++) o[c] = r[0] * B[0][c] + r[1] * B[1][c] + r[2] * B[2][c]; return o; }); for (let r = 0; r < 8; r++) { for (let c = 0; c < 8; c++) Ke[r * 8 + c] += 0.25 * (B[0][r] * CB[0][c] + B[1][r] * CB[1][c] + B[2][r] * CB[2][c]); for (let l = 0; l < 3; l++) Fe[l][r] += 0.25 * (B[0][r] * C[0][l] + B[1][r] * C[1][l] + B[2][r] * C[2][l]); } }
  const Bc = Bat(0.5, 0.5), CBc = C.map(r => { const o = new Float64Array(8); for (let c = 0; c < 8; c++) o[c] = r[0] * Bc[0][c] + r[1] * Bc[1][c] + r[2] * Bc[2][c]; return o; });
  return { C, Ke, Fe, CBc };
}
function solveElastic(ph, n, mat0, mat1, maxIter) {
  const ne = n * n, nd = 2 * ne, el = [q4(mat0.lam, mat0.mu), q4(mat1.lam, mat1.mu)], dof = new Int32Array(ne * 8), dg = new Float64Array(nd);
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const e = i + n * j, ip = (i + 1) % n, jp = (j + 1) % n, nodes = [i + n * j, ip + n * j, ip + n * jp, i + n * jp]; for (let a = 0; a < 4; a++) { dof[e * 8 + 2 * a] = 2 * nodes[a]; dof[e * 8 + 2 * a + 1] = 2 * nodes[a] + 1; } const Ke = el[ph[e]].Ke; for (let a = 0; a < 8; a++) dg[dof[e * 8 + a]] += Ke[a * 9]; }
  const y = new Float64Array(nd), ue = new Float64Array(8);
  const mv = x => { y.fill(0); for (let e = 0; e < ne; e++) { const o = e * 8, Ke = el[ph[e]].Ke; for (let a = 0; a < 8; a++) ue[a] = x[dof[o + a]]; for (let r = 0; r < 8; r++) { let s = 0; const q = r * 8; for (let c = 0; c < 8; c++) s += Ke[q + c] * ue[c]; y[dof[o + r]] += s; } } return y; };
  const Cav = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; for (let e = 0; e < ne; e++) for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) Cav[r][c] += el[ph[e]].C[r][c] / ne;
  const Ceff = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], Wen = [0, 0, 0]; let iters = 0, res = 0, conv = true;
  for (let l = 0; l < 3; l++) {
    const f = new Float64Array(nd); for (let e = 0; e < ne; e++) { const Fe = el[ph[e]].Fe[l]; for (let a = 0; a < 8; a++) f[dof[e * 8 + a]] += Fe[a]; }
    const b = f.map(v => -v), sol = N.cg(mv, b, { tol: 1e-4, maxIter, diag: dg }), u = sol.x; iters += sol.iterations; res = Math.max(res, sol.residual); conv = conv && sol.converged;
    for (let r = 0; r < 3; r++) { let s = Cav[r][l]; for (let e = 0; e < ne; e++) { const CB = el[ph[e]].CBc[r], o = e * 8; let t = 0; for (let c = 0; c < 8; c++) t += CB[c] * u[dof[o + c]]; s += t / ne; } Ceff[r][l] = s; }
    const Ku = mv(u); let uKu = 0, uf = 0; for (let q = 0; q < nd; q++) { uKu += u[q] * Ku[q]; uf += u[q] * f[q]; }
    Wen[l] = Cav[l][l] + (2 * uf + uKu) / ne;
  }
  return { Ceff, Wen, iters, res, conv, Cav };
}

function rve(p) {
  const n = +p.nres, nIncl = Math.round(p.n_incl), warnings = [], gen = (seed, nn, ni) => microGen(p.micro, nn, ni, p.vf, seed);
  const main = gen(p.seed, n, nIncl), f = main.f, k0 = p.k_m, k1 = p.k_i, contrast = Math.max(k0, k1) / Math.min(k0, k1);
  const sol = solveCond(main.ph, n, k0, k1, 1500, true), K = sol.K, kxy = 0.5 * (K[0][1] + K[1][0]), tr = K[0][0] + K[1][1], dd = Math.sqrt(((K[0][0] - K[1][1]) / 2) ** 2 + kxy * kxy), lam1 = tr / 2 + dd, lam2 = tr / 2 - dd, kIso = tr / 2, anis = lam1 / lam2;
  const lo = Math.min(k0, k1), hi = Math.max(k0, k1), fHi = k1 >= k0 ? f : 1 - f, fLo = 1 - fHi, degenerate = hi / lo < 1 + 1e-9 || fHi < 1e-9 || fLo < 1e-9;
  const hs2 = [mg(lo, hi, fHi, 2), mg(hi, lo, fLo, 2)], wien = [1 / (fLo / lo + fHi / hi), fLo * lo + fHi * hi], tol = 0.02;
  const inWiener = lam2 >= wien[0] * (1 - tol) && lam1 <= wien[1] * (1 + tol);
  let traceOK = true, trDetail = 'single phase or equal conductivities: bounds coincide';
  if (!degenerate) { const l1 = Math.max(lam1, lo * (1 + 1e-9)), l2 = Math.max(lam2, lo * (1 + 1e-9)), h1 = Math.min(lam1, hi * (1 - 1e-9)), h2 = Math.min(lam2, hi * (1 - 1e-9)); const L = 1 / (l1 - lo) + 1 / (l2 - lo), RL = (2 / (hi - lo) + fLo / lo) / fHi, U = 1 / (hi - h1) + 1 / (hi - h2), RU = (2 / (hi - lo) - fHi / hi) / fLo; traceOK = L <= RL * (1 + 2 * tol) && U <= RU * (1 + 2 * tol) && kIso >= hs2[0] * (1 - tol) && kIso <= hs2[1] * (1 + tol); trDetail = `mean eigenvalue ${F(kIso, 4)} in [${F(hs2[0], 4)}, ${F(hs2[1], 4)}] W/m/K; trace bounds ${F(L / RL, 3)} and ${F(U / RU, 3)} of their limits`; }
  const hmErr = Math.max(Math.abs(sol.Wen[0] - K[0][0]) / K[0][0], Math.abs(sol.Wen[1] - K[1][1]) / K[1][1]), symErr = Math.abs(K[0][1] - K[1][0]) / kIso;
  // statistics: three realisations at full size, three at a quarter of the area
  const full = [K[0][0]], small = [], nS = n / 2, niS = Math.max(1, Math.round(nIncl / 4));
  for (let s = 1; s < 3; s++) { const g = gen(p.seed + s, n, nIncl), r = solveCond(g.ph, n, k0, k1, 500, false, 1e-4, 1); full.push(r.K[0][0]); }
  for (let s = 0; s < 3; s++) { const g = gen(p.seed + 10 + s, nS, niS), r = solveCond(g.ph, nS, k0, k1, 500, false, 1e-4, 1); small.push(r.K[0][0]); }
  const mFull = N.mean(full), cov = N.std(full) / mFull, mSmall = N.mean(small), sizeDiff = Math.abs(mSmall - mFull) / mFull, covS = N.std(small) / mSmall;
  // elasticity on at most 32 × 32 elements
  const ne = n > 32 ? n / 2 : n, gE = ne === n ? main : gen(p.seed, ne, nIncl), fE = gE.f, lame = (E, nu) => ({ lam: E * nu / ((1 + nu) * (1 - 2 * nu)), mu: E / (2 * (1 + nu)) }), m0 = lame(p.E_m, p.nu_m), m1 = lame(p.E_i, p.nu_i);
  const el = solveElastic(gE.ph, ne, m0, m1, 400), C = el.Ceff, kap = (C[0][0] + C[1][1] + 2 * C[0][1]) / 4, muA = C[2][2], muB = (C[0][0] + C[1][1] - 2 * C[0][1]) / 4, muE = 0.5 * (muA + muB), aE = muA / muB, lamE = kap - muE, Eeq = muE * (3 * lamE + 2 * muE) / (lamE + muE), nuEq = lamE / (2 * (lamE + muE));
  const hmEl = Math.max(...[0, 1, 2].map(l => Math.abs(el.Wen[l] - C[l][l]) / Math.abs(C[l][l]))), symEl = Math.max(Math.abs(C[0][1] - C[1][0]), Math.abs(C[0][2] - C[2][0]), Math.abs(C[1][2] - C[2][1])) / C[0][0];
  const k2 = m => m.lam + m.mu, ordered = (k2(m1) - k2(m0)) * (m1.mu - m0.mu) >= 0; let hsK = [0, 0], hsMu = [0, 0], elOK = true, elDetail = 'phases not well-ordered (κ and μ ranked oppositely): 2-D bounds not evaluated';
  if (ordered && fE > 1e-9 && fE < 1 - 1e-9) { const [s, t, fs, ft] = k2(m1) >= k2(m0) ? [m0, m1, 1 - fE, fE] : [m1, m0, fE, 1 - fE], ks = k2(s), kt = k2(t), same = Math.abs(kt - ks) < 1e-9 * ks, sameM = Math.abs(t.mu - s.mu) < 1e-9 * s.mu; hsK = same ? [ks, ks] : [ks + ft / (1 / (kt - ks) + fs / (ks + s.mu)), kt + fs / (1 / (ks - kt) + ft / (kt + t.mu))]; hsMu = sameM ? [s.mu, s.mu] : [s.mu + ft / (1 / (t.mu - s.mu) + fs * (ks + 2 * s.mu) / (2 * s.mu * (ks + s.mu))), t.mu + fs / (1 / (s.mu - t.mu) + ft * (kt + 2 * t.mu) / (2 * t.mu * (kt + t.mu)))]; elOK = kap >= hsK[0] * (1 - 0.03) && kap <= hsK[1] * (1 + 0.03); elDetail = `κ = ${GPa(kap)} GPa in [${GPa(hsK[0])}, ${GPa(hsK[1])}] GPa; μ = ${GPa(muE)} GPa, 2-D HS range [${GPa(hsMu[0])}, ${GPa(hsMu[1])}] GPa`; }
  const qMean = Math.hypot(K[0][0], K[1][0]), loc = N.maxOf(sol.field) / qMean, rho = (1 - f) * p.rho_m + f * p.rho_i;
  if (!sol.conv) warnings.push(`Conductivity CG solver stopped at residual ${sol.res.toExponential(1)} (contrast ${F(contrast, 3)}): results are approximate; reduce the property contrast or resolution.`);
  if (!el.conv) warnings.push(`Elastic CG solver stopped at residual ${el.res.toExponential(1)}: the effective stiffness is approximate for this stiffness contrast.`);
  if (p.micro === 'checker') warnings.push(`The checkerboard is fixed at 50 % volume fraction; the requested fraction is ignored. Conduction through the corner contacts is under-resolved by the grid: the numerical value is ${F(100 * Math.abs(kIso - Math.sqrt(k0 * k1)) / Math.sqrt(k0 * k1), 3)} % from the exact √(k₁k₂), and the error grows with property contrast.`);
  if (Math.abs(f - p.vf) > 0.03 && p.micro !== 'checker') warnings.push(`Realised volume fraction ${F(f, 3)} differs from the requested ${F(p.vf, 3)} (pixelation, overlap or few grains); bounds use the realised value.`);
  if (ne !== n) warnings.push(`Elastic cell problem solved on ${ne}×${ne} elements (same geometry, coarser raster, phase fraction ${F(fE, 3)}) to keep the run time short.`);
  warnings.push('2-D plane-strain cell: represents columnar/fibre cross-sections. Equiaxed 3-D particles conduct and stiffen differently — compare with the 3-D mean-field estimates.');
  const map = new Float64Array(n * n); for (let i = 0; i < n * n; i++) map[i] = main.ph[i];
  const checkerExact = Math.sqrt(k0 * k1);
  return {
    summary: `Periodic ${n}×${n} RVE (${p.micro}, phase-B fraction ${F(f, 3)}): effective conductivity tensor has principal values ${F(lam1, 4)} and ${F(lam2, 4)} W/m/K (anisotropy ${F(anis, 4)}), ${traceOK ? 'inside' : 'OUTSIDE'} the 2-D Hashin–Shtrikman bounds ${F(hs2[0], 4)}–${F(hs2[1], 4)} W/m/K. In-plane stiffness: C11 = ${GPa(C[0][0])}, C22 = ${GPa(C[1][1])}, C12 = ${GPa(C[0][1])}, C66 = ${GPa(C[2][2])} GPa (equivalent E ≈ ${GPa(Eeq)} GPa, ν ≈ ${F(nuEq, 3)}). Realisation scatter ${F(100 * cov, 3)} % and RVE-size change ${F(100 * sizeDiff, 3)} %.`,
    kpis: [
      kpi('k_iso', 'Effective conductivity (mean)', kIso, 'W/m/K', `Half the trace of the homogenised conductivity tensor from two periodic cell problems (macroscopic gradients along x and y), solved by finite volumes and conjugate gradients in ${sol.iters} iterations.`, traceOK),
      kpi('k_xx', 'k_xx', K[0][0], 'W/m/K', 'Effective conductivity for a macroscopic gradient along x: volume-averaged flux per unit applied gradient.', null),
      kpi('k_yy', 'k_yy', K[1][1], 'W/m/K', 'Effective conductivity for a macroscopic gradient along y. Differences from k_xx reveal anisotropy of this realisation.', null),
      kpi('anis_k', 'Conductivity anisotropy', anis, '–', 'Ratio of the principal conductivities. For a statistically isotropic microstructure it tends to 1 as the RVE grows; values above about 1.1 indicate too few inclusions.', anis < 1.15),
      kpi('C11', 'Effective C11 (plane strain)', C[0][0], 'Pa', 'In-plane stiffness component from three periodic finite-element cell problems (bilinear quadrilaterals, unit macroscopic strains).', null),
      kpi('C12', 'Effective C12', C[0][1], 'Pa', 'Coupling stiffness between the two in-plane normal strains.', null),
      kpi('C66', 'Effective in-plane shear C66', C[2][2], 'Pa', 'Shear stiffness for engineering shear strain in the cell axes.', null),
      kpi('E_eq', 'Equivalent isotropic E', Eeq, 'Pa', 'Young\'s modulus of the isotropic material having the same plane-strain bulk and mean shear moduli as the RVE (λ = κ − μ).', null),
      kpi('nu_eq', 'Equivalent Poisson\'s ratio', nuEq, '–', 'Poisson\'s ratio of the same equivalent isotropic material, ν = λ/(2(λ + μ)).', nuEq > -1 && nuEq < 0.5),
      kpi('anis_el', 'Elastic anisotropy 2C66/(C11+C22−2C12)·2', aE, '–', 'Ratio of the two in-plane shear moduli (axes shear over 45° shear); 1 for an in-plane isotropic response.', aE > 0.87 && aE < 1.15),
      kpi('cov', 'Realisation scatter of k', cov, '–', 'Coefficient of variation of k_xx over three seeded realisations at full RVE size: a direct measure of statistical representativeness.', cov < 0.05),
      kpi('size_err', 'RVE-size sensitivity', sizeDiff, '–', `Relative change of the realisation-averaged k_xx between cells of ${nS}×${nS} (${niS} features) and ${n}×${n} (${nIncl} features); small values indicate a converged RVE.`, sizeDiff < 0.05),
      kpi('hill', 'Hill–Mandel energy mismatch', Math.max(hmErr, hmEl), '–', 'Largest relative difference between the volume-averaged microscopic energy and the macroscopic energy computed from the averaged flux or stress. It should vanish to solver tolerance.', Math.max(hmErr, hmEl) < 1e-4),
      kpi('loc', 'Flux localisation factor', loc, '–', 'Maximum local flux magnitude divided by the mean flux for the x-gradient load case: the microscale amplification relevant to hot spots and breakdown.', null),
    ],
    series: [
      { key: 'conv', title: 'RVE-size and realisation convergence', xlabel: 'RVE side (grid cells)', ylabel: 'effective conductivity k_xx (W/m/K)', lines: [line('realisations', [...small.map(() => nS), ...full.map(() => n)], [...small, ...full]), line('mean', [nS, n], [mSmall, mFull]), line('2-D HS lower', [nS, n], [hs2[0], hs2[0]]), line('2-D HS upper', [nS, n], [hs2[1], hs2[1]])], explain: 'Each point is one seeded realisation. The scatter shrinks and the mean stabilises as the cell becomes representative; all points should lie between the Hashin–Shtrikman bounds.' },
      { key: 'cg', title: 'Estimates compared', xlabel: 'estimate index', ylabel: 'conductivity (W/m/K)', lines: [line('Wiener series, HS lower, RVE mean, HS upper, Wiener parallel', [1, 2, 3, 4, 5], [wien[0], hs2[0], kIso, hs2[1], wien[1]])], explain: 'Ordered comparison of the rigorous bounds with the numerical RVE result (index 3). The RVE value must sit between indices 2 and 4 for an isotropic microstructure.' },
    ],
    maps: [
      { key: 'phase', title: 'RVE microstructure', nx: n, ny: n, data: map, x0: 0, x1: 1, y0: 0, y1: 1, xlabel: 'x / L', ylabel: 'y / L', unit: 'phase (1 = B)', explain: 'Seeded periodic two-phase microstructure used for the cell problems; phase B is the inclusion phase.' },
      { key: 'flux', title: 'Local flux magnitude (gradient along x)', nx: n, ny: n, data: sol.field, x0: 0, x1: 1, y0: 0, y1: 1, xlabel: 'x / L', ylabel: 'y / L', unit: 'W/m² per unit gradient', explain: 'Magnitude of the local heat flux for a unit macroscopic temperature gradient along x. Bright paths are the preferred conduction channels; necks between inclusions concentrate flux.' },
    ],
    tables: [
      { key: 'ktensor', title: 'Effective conductivity tensor', columns: ['', 'x', 'y'], rows: [['x', F(K[0][0], 5), F(K[0][1], 5)], ['y', F(K[1][0], 5), F(K[1][1], 5)]], explain: 'Homogenised conductivity tensor in W/m/K. Off-diagonal terms should be equal (symmetry) and small for an isotropic microstructure.' },
      { key: 'ctensor', title: 'Effective plane-strain stiffness (GPa)', columns: ['', 'ε11', 'ε22', 'γ12'], rows: [['σ11', GPa(C[0][0]), GPa(C[0][1]), GPa(C[0][2])], ['σ22', GPa(C[1][0]), GPa(C[1][1]), GPa(C[1][2])], ['σ12', GPa(C[2][0]), GPa(C[2][1]), GPa(C[2][2])]], explain: 'Homogenised in-plane stiffness matrix in Voigt notation (engineering shear strain); it can be used directly as an anisotropic material card.' },
    ],
    equations: [
      eqn('governing', 'Conductivity cell problem', '∇·[k(x)(ē + ∇χ)] = 0 in the cell, χ periodic;  K*ē = ⟨k(ē + ∇χ)⟩'),
      eqn('governing', 'Elastic cell problem', '∇·[C(x) : (ε̄ + ε(ũ))] = 0, ũ periodic;  C* : ε̄ = ⟨σ⟩'),
      eqn('formula', 'Hill–Mandel condition', '⟨σ : ε⟩ = ⟨σ⟩ : ⟨ε⟩;   ⟨q·g⟩ = ⟨q⟩·⟨g⟩'),
      eqn('formula', '2-D Hashin–Shtrikman bounds (conductivity)', 'k₁ + f₂/[1/(k₂ − k₁) + f₁/(2k₁)] ≤ k* ≤ k₂ + f₁/[1/(k₁ − k₂) + f₂/(2k₂)]'),
      eqn('formula', 'Checkerboard benchmark (Dykhne)', 'k* = √(k₁k₂)', `= ${F(checkerExact, 4)} W/m/K for these phases`),
    ],
    geometry: { kind: 'voxels', nx: n, ny: n, nz: 1, data: main.ph, voxel: 1 / n, units: 'cell' },
    checks: [
      chk('Conductivity within Hashin–Shtrikman bounds', traceOK, trDetail),
      chk('Principal conductivities within Wiener bounds', inWiener, `[${F(lam2, 4)}, ${F(lam1, 4)}] in [${F(wien[0], 4)}, ${F(wien[1], 4)}] W/m/K`),
      chk('Hill–Mandel energy consistency (conduction)', hmErr < 1e-4, `relative mismatch ${hmErr.toExponential(2)}`),
      chk('Hill–Mandel energy consistency (elasticity)', hmEl < 1e-4, `relative mismatch ${hmEl.toExponential(2)}`),
      chk('Effective tensors symmetric', symErr < 1e-4 && symEl < 1e-4, `conductivity ${symErr.toExponential(1)}, stiffness ${symEl.toExponential(1)}`),
      chk('Plane-strain bulk modulus within 2-D HS bounds', elOK, elDetail),
      chk('Linear solvers converged', sol.conv && el.conv, `residuals ${sol.res.toExponential(1)} (conduction), ${el.res.toExponential(1)} (elasticity)`),
      chk('RVE statistically representative', cov < 0.05 && sizeDiff < 0.05, `scatter ${F(100 * cov, 3)} %, size change ${F(100 * sizeDiff, 3)} %`),
      ...(p.micro === 'checker' ? [chk('Checkerboard benchmark', Math.abs(kIso - checkerExact) / checkerExact < 0.15, `numerical ${F(kIso, 4)} vs exact ${F(checkerExact, 4)} W/m/K (corner singularities converge slowly)`)] : []),
    ],
    warnings,
    recommendations: [
      rec('high', cov > 0.05 || sizeDiff > 0.05 ? `Increase the number of inclusions (now ${nIncl}) or average more realisations: scatter ${F(100 * cov, 3)} % and size sensitivity ${F(100 * sizeDiff, 3)} % exceed the 5 % guideline for a representative volume.` : 'The cell is representative to within 5 %; use the mean and scatter as the property and its uncertainty in component models.', 'RVE convergence indicators.'),
      rec('high', 'Replace the synthetic microstructure with a segmented micrograph or tomography slice of the real material before quoting design values; mean-field shapes rarely match measured morphology.', 'Image-based homogenisation practice.'),
      rec('medium', loc > 3 ? `Local flux reaches ${F(loc, 3)} times the mean: check hot spots, electrical breakdown or interface failure at inclusion necks.` : 'Localisation is moderate; microscale failure is unlikely to initiate far below the homogenised limit.', `Localisation factor ${F(loc, 3)}.`),
      rec('medium', 'Use homogenised properties in component simulations instead of resolving every inclusion: the computational cost and energy use of design studies fall by orders of magnitude, with the error quantified by the scatter above.', 'Scale separation and computational homogenisation.', true),
    ],
    publish: { 'eff.k': kIso, 'eff.E': Eeq, 'eff.nu': nuEq, 'eff.G': muE, 'eff.K': kap - muE / 3, 'eff.density': rho, 'cred.rve_cov': cov },
    numerics: { method: 'Periodic cell problems: cell-centred finite volumes (conduction) and Q4 plane-strain finite elements (elasticity), Jacobi-preconditioned conjugate gradients', dof: n * n + 2 * ne * ne, steps: 5, iterations: sol.iters + el.iters, residual: Math.max(sol.res, el.res), notes: `Conduction ${n}×${n}, elasticity ${ne}×${ne}; 5 additional single-direction conduction solves for statistics.` },
  };
}

// ---------- laminate ----------
function laminate(p) {
  const fB = p.vf, fA = 1 - fB, per = p.period_nm * 1e-9, nb = Math.round(p.n_bi), H = nb * per, warnings = [];
  const lm = (E, nu) => ({ lam: E * nu / ((1 + nu) * (1 - 2 * nu)), mu: E / (2 * (1 + nu)) }), A = lm(p.E_m, p.nu_m), B = lm(p.E_i, p.nu_i), av = (a, b) => fA * a + fB * b;
  const backus = () => { const M = x => x.lam + 2 * x.mu, C33 = 1 / av(1 / M(A), 1 / M(B)), r = av(A.lam / M(A), B.lam / M(B)), C13 = r * C33, C11 = av(4 * A.mu * (A.lam + A.mu) / M(A), 4 * B.mu * (B.lam + B.mu) / M(B)) + r * r * C33, C44 = 1 / av(1 / A.mu, 1 / B.mu), C66 = av(A.mu, B.mu); return { C11, C12: C11 - 2 * C66, C13, C33, C44, C66 }; };
  const c = backus(), S = N.inverse([[c.C11, c.C12, c.C13], [c.C12, c.C11, c.C13], [c.C13, c.C13, c.C33]]), Ein = 1 / S[0][0], Eout = 1 / S[2][2], nu12 = -S[0][1] / S[0][0], nu13 = -S[0][2] / S[0][0];
  const kIn = av(p.k_m, p.k_i), nInt = Math.max(0, 2 * nb - 1), kOutIdeal = 1 / av(1 / p.k_m, 1 / p.k_i), kOut = H / (H / kOutIdeal + nInt * p.R_int), sIn = av(p.sig_m, p.sig_i), sOut = 1 / av(1 / p.sig_m, 1 / p.sig_i);
  const bA = p.E_m / (1 - p.nu_m), bB = p.E_i / (1 - p.nu_i), aIn = (fA * bA * p.cte_m + fB * bB * p.cte_i) / (fA * bA + fB * bB);
  const aOut = fA * (p.cte_m + 2 * p.nu_m / (1 - p.nu_m) * (p.cte_m - aIn)) + fB * (p.cte_i + 2 * p.nu_i / (1 - p.nu_i) * (p.cte_i - aIn));
  const sA = bA * (aIn - p.cte_m) * p.dT, sB = bB * (aIn - p.cte_i) * p.dT, rho = av(p.rho_m, p.rho_i);
  // discrete-stack bending stiffness versus homogenised plate
  const QA = p.E_m / (1 - p.nu_m ** 2), QB = p.E_i / (1 - p.nu_i ** 2); let Am = 0, Bm = 0, Dm = 0, z = -H / 2;
  for (let k = 0; k < nb; k++) for (const [Q, t] of [[QA, fA * per], [QB, fB * per]]) { const z1 = z + t; Am += Q * t; Bm += Q * (z1 * z1 - z * z) / 2; Dm += Q * (z1 ** 3 - z ** 3) / 3; z = z1; }
  const Dred = Dm - Bm * Bm / Am, Dhom = (Am / H) * H ** 3 / 12, dRatio = Dred / Dhom;
  const kapRatio = kOut / kOutIdeal;
  if (per < 2e-9) warnings.push('Bilayer period below about 2 nm: layers are only a few atomic planes thick, so bulk constituent properties and sharp interfaces are questionable.');
  if (nb < 5) warnings.push(`Only ${nb} bilayer(s): the stack is not scale-separated in bending (stiffness ratio ${F(dRatio, 3)}); use discrete laminate theory rather than homogenised properties.`);
  if (Math.max(Math.abs(sA), Math.abs(sB)) > 1e9) warnings.push('Thermal-mismatch stress exceeds 1 GPa: yielding, cracking or delamination is likely before this temperature change is reached.');
  const fs = N.linspace(0.01, 0.99, 50), at = ff => laminate._at(p, ff), periods = N.logspace(1e-9, 1e-6, 40);
  const nRows = 48, lay = new Uint8Array(16 * nRows), show = Math.min(nb, 6); for (let r = 0; r < nRows; r++) { const u = (r / nRows * show) % 1; for (let x = 0; x < 16; x++) lay[x + 16 * r] = u < fA ? 0 : 1; }
  return {
    summary: `${nb}-bilayer laminate of period ${F(p.period_nm)} nm (${F(100 * fB, 3)} % B, total ${F(H * 1e6, 3)} µm): in-plane E = ${GPa(Ein)} GPa, through-thickness E = ${GPa(Eout)} GPa; thermal conductivity ${F(kIn, 4)} W/m/K in-plane and ${F(kOut, 4)} W/m/K cross-plane (interfaces reduce it to ${F(100 * kapRatio, 3)} % of the ideal series value); expansion ${F(aIn * 1e6, 4)} (in-plane) and ${F(aOut * 1e6, 4)} ×10⁻⁶ 1/K (through-thickness).`,
    kpis: [
      kpi('E_in', 'In-plane Young\'s modulus', Ein, 'Pa', 'From the transversely isotropic stiffness obtained by Backus (exact layered-medium) averaging of the two isotropic layers; loading parallel to the layers.', null),
      kpi('E_out', 'Through-thickness Young\'s modulus', Eout, 'Pa', 'Modulus for loading normal to the layers, dominated by the softer layer (series coupling with lateral constraint).', null),
      kpi('nu12', 'In-plane Poisson\'s ratio', nu12, '–', 'Lateral in-plane contraction per unit in-plane extension, from the homogenised compliance.', nu12 > -1 && nu12 < 1),
      kpi('nu13', 'Through-thickness Poisson\'s ratio', nu13, '–', 'Thickness contraction per unit in-plane extension; it can differ strongly from the constituent values in high-contrast stacks.', null),
      kpi('G_out', 'Interlaminar shear modulus C44', c.C44, 'Pa', 'Shear stiffness for sliding of layers over each other (harmonic mean of the layer shear moduli).', null),
      kpi('k_in', 'In-plane thermal conductivity', kIn, 'W/m/K', 'Layers conduct in parallel: thickness-weighted arithmetic mean.', null),
      kpi('k_out', 'Cross-plane thermal conductivity', kOut, 'W/m/K', `Layers and ${nInt} interfaces in series, each interface adding ${F(p.R_int, 3)} m²K/W. Without interface resistance the value would be ${F(kOutIdeal, 4)} W/m/K.`, null),
      kpi('k_anis', 'Thermal anisotropy k∥/k⊥', kIn / kOut, '–', 'Ratio of in-plane to cross-plane conductivity; large values make the laminate a heat spreader in-plane and a barrier across.', null),
      kpi('cte_in', 'In-plane thermal expansion', aIn, '1/K', 'Biaxial-modulus-weighted mean of the layer expansions (layers constrained to equal in-plane strain).', null),
      kpi('cte_out', 'Through-thickness thermal expansion', aOut, '1/K', 'Free thickness expansion including the Poisson response of each layer to its in-plane mismatch stress.', null),
      kpi('s_A', 'Mismatch stress in layer A', sA, 'Pa', `Equi-biaxial stress in layer A after a temperature change of ${F(p.dT)} K, σ = E/(1 − ν)·(α∥ − α_A)ΔT.`, Math.abs(sA) < 1e9),
      kpi('s_B', 'Mismatch stress in layer B', sB, 'Pa', 'Equi-biaxial stress in layer B for the same temperature change; the layer forces balance (f_Aσ_A + f_Bσ_B = 0).', Math.abs(sB) < 1e9),
      kpi('D_ratio', 'Bending stiffness: discrete / homogenised', dRatio, '–', 'Reduced bending stiffness D − B²/A of the actual stacking sequence divided by that of a homogeneous plate with the same in-plane stiffness; it tends to 1 for many layers.', Math.abs(dRatio - 1) < 0.05),
      kpi('rho', 'Density', rho, 'kg/m³', 'Thickness-weighted mean density of the stack.', null),
    ],
    series: [
      { key: 'E_f', title: 'Moduli versus thickness fraction of B', xlabel: 'thickness fraction of B', ylabel: 'modulus (GPa)', lines: [line('E in-plane', fs, fs.map(ff => at(ff).Ein / 1e9)), line('E through-thickness', fs, fs.map(ff => at(ff).Eout / 1e9))], explain: 'In-plane stiffness follows a near-linear (parallel) mixture; the through-thickness modulus stays close to the soft layer until the stiff fraction is high.' },
      { key: 'k_period', title: 'Cross-plane conductivity versus bilayer period', xlabel: 'bilayer period (nm)', ylabel: 'k cross-plane (W/m/K)', logx: true, lines: [line('with interface resistance', periods.map(v => v * 1e9), periods.map(L => L / (L / kOutIdeal + 2 * p.R_int))), line('ideal series limit', periods.map(v => v * 1e9), periods.map(() => kOutIdeal))], explain: 'Interface (Kapitza) resistance dominates when the period is short: nanolaminates can conduct far less across the layers than either constituent.' },
      { key: 'a_f', title: 'Thermal expansion versus thickness fraction of B', xlabel: 'thickness fraction of B', ylabel: 'CTE (10⁻⁶ 1/K)', lines: [line('in-plane', fs, fs.map(ff => at(ff).aIn * 1e6)), line('through-thickness', fs, fs.map(ff => at(ff).aOut * 1e6))], explain: 'In-plane expansion is pulled towards the stiffer layer; the through-thickness value can exceed both constituents because constrained layers bulge outward.' },
    ],
    equations: [
      eqn('formula', 'Backus averaging (isotropic layers)', 'C33 = ⟨1/(λ+2μ)⟩⁻¹;  C13 = ⟨λ/(λ+2μ)⟩C33;  C11 = ⟨4μ(λ+μ)/(λ+2μ)⟩ + ⟨λ/(λ+2μ)⟩²C33;  C44 = ⟨1/μ⟩⁻¹;  C66 = ⟨μ⟩'),
      eqn('formula', 'Conductivity', 'k∥ = Σ fᵢkᵢ;   k⊥ = H / (Σ tᵢ/kᵢ + N_int·R_int)'),
      eqn('formula', 'Thermal expansion', 'α∥ = Σ fᵢEᵢ′αᵢ / Σ fᵢEᵢ′,  Eᵢ′ = Eᵢ/(1 − νᵢ);   α⊥ = Σ fᵢ[αᵢ + 2νᵢ/(1 − νᵢ)·(αᵢ − α∥)]'),
      eqn('constitutive', 'Classical laminate theory', '[N; M] = [A B; B D][ε⁰; κ],  A, B, D = Σ Q̄ₖ (zₖ − zₖ₋₁), ½(zₖ² − zₖ₋₁²), ⅓(zₖ³ − zₖ₋₁³)'),
    ],
    geometry: { kind: 'voxels', nx: 16, ny: nRows, nz: 1, data: lay, voxel: show * p.period_nm / nRows, units: 'nm' },
    checks: [
      chk('Homogenised stiffness positive definite', c.C44 > 0 && c.C66 > 0 && c.C33 > 0 && (c.C11 + c.C12) * c.C33 > 2 * c.C13 ** 2, `C11 = ${GPa(c.C11)}, C33 = ${GPa(c.C33)}, C13 = ${GPa(c.C13)} GPa`),
      chk('In-plane modulus between Reuss and Voigt', Ein <= Math.max(QA, QB) && Ein >= 1 / av(1 / p.E_m, 1 / p.E_i) * (1 - 1e-9), `E∥ = ${GPa(Ein)} GPa`),
      chk('Mismatch forces balance', Math.abs(fA * sA + fB * sB) <= 1e-6 * (Math.abs(sA) + Math.abs(sB) + 1), `f_Aσ_A + f_Bσ_B = ${F(fA * sA + fB * sB, 3)} Pa`),
      chk('Conductivity ordering k⊥ ≤ k∥', kOut <= kIn * (1 + 1e-12), `${F(kOut, 4)} ≤ ${F(kIn, 4)} W/m/K`),
      chk('Scale separation in bending', Math.abs(dRatio - 1) < 0.05, `D ratio ${F(dRatio, 4)} with ${nb} bilayers`),
    ],
    warnings,
    recommendations: [
      rec('high', kapRatio < 0.5 ? `Interfaces carry most of the cross-plane thermal resistance (${F(100 * (1 - kapRatio), 3)} %): exploit this for thermal barriers, or increase the period if heat must be removed through the stack.` : 'Interface resistance is minor at this period; shorten the period below about R_int·k to obtain a thermal-barrier effect.', 'Series thermal resistance with Kapitza interfaces.'),
      rec('high', Math.max(Math.abs(sA), Math.abs(sB)) > 3e8 ? `Mismatch stresses of ${F(sA / 1e6, 3)} and ${F(sB / 1e6, 3)} MPa arise over ${F(p.dT)} K: check against layer strength and interface toughness, or adjust the thickness ratio to shift stress into the tougher layer.` : 'Thermal-mismatch stresses are moderate for this temperature excursion.', 'Biaxial mismatch stress of constrained layers.'),
      rec('medium', nb < 5 ? 'Add layers or use a symmetric stacking sequence to remove bending–stretching coupling.' : 'The stack is well scale-separated: the homogenised transversely isotropic constants can be used as a material card for component analysis.', `Bending stiffness ratio ${F(dRatio, 3)}.`),
      rec('low', 'Choose deposition thicknesses that achieve the property target with the fewest layers: each additional layer adds vacuum-process time and energy roughly in proportion.', 'Number of deposition cycles scales with layer count.', true),
    ],
    publish: { 'eff.E': Ein, 'eff.nu': nu12, 'eff.G': c.C66, 'eff.k': kOut, 'eff.sigma': sIn, 'eff.cte': aIn, 'eff.density': rho, 'geom.thickness': H, 'eff.sigma_cross': sOut },
    numerics: { method: 'Backus layered-medium averaging; classical laminate A-B-D integration over the discrete stack', dof: 6, steps: 2 * nb, iterations: 1, residual: 0, notes: 'Isotropic layers, perfect mechanical bonding, bulk layer properties; interface thermal resistance lumped per interface.' },
  };
}
laminate._at = (p, fB) => {
  const fA = 1 - fB, lm = (E, nu) => ({ lam: E * nu / ((1 + nu) * (1 - 2 * nu)), mu: E / (2 * (1 + nu)) }), A = lm(p.E_m, p.nu_m), B = lm(p.E_i, p.nu_i), av = (a, b) => fA * a + fB * b, M = x => x.lam + 2 * x.mu;
  const C33 = 1 / av(1 / M(A), 1 / M(B)), r = av(A.lam / M(A), B.lam / M(B)), C13 = r * C33, C11 = av(4 * A.mu * (A.lam + A.mu) / M(A), 4 * B.mu * (B.lam + B.mu) / M(B)) + r * r * C33, C66 = av(A.mu, B.mu), C12 = C11 - 2 * C66;
  const S = N.inverse([[C11, C12, C13], [C12, C11, C13], [C13, C13, C33]]), bA = p.E_m / (1 - p.nu_m), bB = p.E_i / (1 - p.nu_i), aIn = (fA * bA * p.cte_m + fB * bB * p.cte_i) / (fA * bA + fB * bB);
  return { Ein: 1 / S[0][0], Eout: 1 / S[2][2], aIn, aOut: fA * (p.cte_m + 2 * p.nu_m / (1 - p.nu_m) * (p.cte_m - aIn)) + fB * (p.cte_i + 2 * p.nu_i / (1 - p.nu_i) * (p.cte_i - aIn)) };
};

// ---------- lattice-to-continuum ----------
const TOPO = {
  octet: { name: 'octet truss', Ltot: 24 / Math.SQRT2, l: 1 / Math.SQRT2, stretch: true, C: r => ({ C11: r / 6, C12: r / 12, C44: r / 12 }) },
  sc: { name: 'simple cubic', Ltot: 3, l: 1, stretch: true, C: r => ({ C11: r / 3, C12: 0, C44: r * r / 18 }) },
  bcc: { name: 'body-centred cubic', Ltot: 4 * Math.sqrt(3), l: Math.sqrt(3) / 2, stretch: false, C: r => { const K = r / 9, E = r * r, Cp = 3 * K * E / (9 * K - E); return { C11: K + 4 * Cp / 3, C12: K - 2 * Cp / 3, C44: r / 9 }; } },
  kelvin: { name: 'Kelvin cell', Ltot: 24 / (2 * Math.SQRT2), l: 1 / (2 * Math.SQRT2), stretch: false, C: r => { const K = r / 9, E = r * r, Cp = 3 * K * E / (9 * K - E), G = 3 * E / 8; return { C11: K + 4 * Cp / 3, C12: K - 2 * Cp / 3, C44: G }; } }, // hydrostatic load is carried by strut stretching (K = ρ̄E_s/9); E[100] and G follow the Gibson–Ashby bending law
};
function cubicFrom(c) { const Cp = (c.C11 - c.C12) / 2, K = (c.C11 + 2 * c.C12) / 3, S11 = (c.C11 + c.C12) / ((c.C11 - c.C12) * (c.C11 + 2 * c.C12)), S12 = -c.C12 / ((c.C11 - c.C12) * (c.C11 + 2 * c.C12)), S44 = 1 / c.C44, GV = (2 * Cp + 3 * c.C44) / 5, GR = 5 * Cp * c.C44 / (2 * c.C44 + 3 * Cp), G = 0.5 * (GV + GR); return { K, Cp, S11, S12, S44, E100: 1 / S11, E111: 1 / (S11 - 2 * (S11 - S12 - S44 / 2) / 3), nu100: -S12 / S11, A: c.C44 / Cp, G, E: Eof(K, G), nu: nuof(K, G), Edir: J => 1 / (S11 - 2 * (S11 - S12 - S44 / 2) * J) }; }
function latticeGeom(topo, a) { // one or a few unit cells as a strut network, coordinates in units of a
  let pts = [], target;
  const add = q => { if (!pts.some(o => Math.abs(o[0] - q[0]) + Math.abs(o[1] - q[1]) + Math.abs(o[2] - q[2]) < 1e-6)) pts.push(q); };
  if (topo === 'kelvin') { for (const perm of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) for (const s1 of [-1, 1]) for (const s2 of [-1, 1]) { const v = [0, 0, 0]; v[perm[1]] = s1 * 0.25; v[perm[2]] = s2 * 0.5; add(v.map(x => x + 0.5)); } target = Math.SQRT2 / 4; }
  else { const m = 2; for (let i = 0; i <= m; i++) for (let j = 0; j <= m; j++) for (let k = 0; k <= m; k++) { add([i, j, k]); if (topo === 'bcc' && i < m && j < m && k < m) add([i + 0.5, j + 0.5, k + 0.5]); if (topo === 'octet') { if (i < m && j < m) add([i + 0.5, j + 0.5, k]); if (i < m && k < m) add([i + 0.5, j, k + 0.5]); if (j < m && k < m) add([i, j + 0.5, k + 0.5]); } } target = TOPO[topo].l; }
  const edges = []; for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) if (Math.abs(Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1], pts[i][2] - pts[j][2]) - target) < 1e-6) edges.push([i, j]);
  return { nodes: pts.map(q => q.map(v => v * a)), edges };
}
function lattice(p) {
  const warnings = [], rel = p.rel, Es = p.E_m, a = p.cell, rhoEff = rel * p.rho_m;
  if (p.topo === 'reentrant') {
    let th = p.theta * N.PI / 180; if (Math.abs(th) < 2 * N.PI / 180) { th = (th < 0 ? -1 : 1) * 2 * N.PI / 180; warnings.push('Angle magnitude raised to 2°: at θ = 0 the bending model gives zero transverse strain coupling and the formulas degenerate.'); }
    let hl = p.h_l; const s = Math.sin(th), c = Math.cos(th); if (hl + 2 * s < 0.2) { hl = 0.2 - 2 * s; warnings.push(`h/l increased to ${F(hl, 3)} so that opposite re-entrant vertices do not touch (h/l > 2|sin θ|).`); }
    const g = hl + s, tl = Math.min(0.5, rel * 2 * c * g / (hl + 2)), l = a / (2 * c), t = tl * l;
    const E1 = Es * tl ** 3 * c / (g * s * s), E2 = Es * tl ** 3 * g / c ** 3, nu12 = c * c / (g * s), nu21 = g * s / (c * c), G12 = Es * tl ** 3 * g / (hl * hl * (1 + 2 * hl) * c);
    const stab = nu12 * nu21 < 1 + 1e-9, cT = Math.sqrt(G12 / rhoEff), lamT = cT / p.f_dyn, sep = lamT / a, nCells = p.L_struct / a, fStrut = 22.37 / (2 * N.PI * l * l) * (t / Math.sqrt(12)) * Math.sqrt(Es / p.rho_m);
    const valid = sep >= 10 && nCells >= 10 && p.f_dyn < 0.3 * fStrut, aux = nu12 < 0;
    if (tl > 0.25) warnings.push(`Wall thickness ratio t/l = ${F(tl, 3)} exceeds 0.25: slender-beam (bending-only) honeycomb theory over-estimates compliance; include axial and shear deformation.`);
    const ths = N.linspace(-60, 45, 106).filter(v => Math.abs(v) >= 2 && hl + 2 * Math.sin(v * N.PI / 180) > 0.2), nuOf = v => { const r = v * N.PI / 180; return Math.cos(r) ** 2 / ((hl + Math.sin(r)) * Math.sin(r)); }, E1of = v => { const r = v * N.PI / 180, tt = Math.min(0.5, rel * 2 * Math.cos(r) * (hl + Math.sin(r)) / (hl + 2)); return tt ** 3 * Math.cos(r) / ((hl + Math.sin(r)) * Math.sin(r) ** 2); };
    const nodes = [], edges = [], id = new Map(), nd = (x, y) => { const k = Math.round(x * 1e6) + '_' + Math.round(y * 1e6); if (!id.has(k)) { id.set(k, nodes.length); nodes.push([x * 1e3, y * 1e3, 0]); } return id.get(k); }, h = hl * l, cw = 2 * l * c, chh = 2 * (h + l * s);
    for (let j = 0; j < 3; j++) for (let i = 0; i < 4; i++) { const x = i * cw, y = j * chh, A0 = nd(x, y), B0 = nd(x, y + h), C0 = nd(x + l * c, y + h + l * s), D0 = nd(x + l * c, y + 2 * h + l * s); edges.push([A0, B0], [B0, C0], [C0, D0]); if (i > 0) edges.push([B0, nd(x - l * c, y + h + l * s)]); if (j < 2) { edges.push([D0, nd(x, y + chh)]); if (i < 3) edges.push([D0, nd(x + cw, y + chh)]); } }
    return {
      summary: `${aux ? 'Re-entrant (auxetic)' : 'Conventional'} 2-D honeycomb at θ = ${F(th * 180 / N.PI, 3)}°, h/l = ${F(hl, 3)}, relative density ${F(rel, 3)} (t/l = ${F(tl, 3)}): E₁ = ${F(E1 / 1e6, 4)} MPa, E₂ = ${F(E2 / 1e6, 4)} MPa, ν₁₂ = ${F(nu12, 4)}, ν₂₁ = ${F(nu21, 4)}, G₁₂ = ${F(G12 / 1e6, 4)} MPa. Homogenisation is ${valid ? 'valid' : 'NOT valid'} at ${F(p.f_dyn, 3)} Hz for a ${F(nCells, 3)}-cell-wide part (wavelength/cell = ${F(sep, 3)}).`,
      kpis: [
        kpi('E1', 'Effective modulus E₁', E1, 'Pa', 'Gibson–Ashby cell-wall bending model, loading along x₁ (perpendicular to the vertical walls); scales with (t/l)³.', null),
        kpi('E2', 'Effective modulus E₂', E2, 'Pa', 'Same model for loading along x₂ (parallel to the vertical walls).', null),
        kpi('nu12', 'Poisson\'s ratio ν₁₂', nu12, '–', `cos²θ/((h/l + sin θ)·sin θ): negative for re-entrant angles, so the lattice widens when stretched. It is purely geometric in the slender-wall limit.`, null),
        kpi('nu21', 'Poisson\'s ratio ν₂₁', nu21, '–', 'Reciprocal coupling (h/l + sin θ)·sin θ/cos²θ; reciprocity requires E₁ν₂₁ = E₂ν₁₂.', null),
        kpi('G12', 'In-plane shear modulus', G12, 'Pa', 'Gibson–Ashby shear modulus; auxetic honeycombs remain shear-compliant even when the magnitude of ν is large.', null),
        kpi('rho_eff', 'Effective density', rhoEff, 'kg/m³', 'Relative density times the solid density.', null),
        kpi('t_wall', 'Wall thickness', t, 'm', `From the relative density relation ρ̄ = (t/l)(h/l + 2)/(2cos θ(h/l + sin θ)) with oblique wall length ${F(l, 3)} m.`, tl <= 0.25),
        kpi('sep', 'Shear wavelength / cell size', sep, '–', `Shear wavelength ${F(lamT, 3)} m at ${F(p.f_dyn, 3)} Hz divided by the cell width; at least 10 is required for a quasi-static effective medium.`, sep >= 10),
        kpi('n_cells', 'Cells across the part', nCells, '–', 'Smallest structural dimension divided by the cell size; fewer than about 10 cells give edge-dominated, size-dependent stiffness.', nCells >= 10),
        kpi('f_strut', 'First wall bending resonance', fStrut, 'Hz', 'Clamped–clamped flexural resonance of an oblique wall; near it the lattice shows local-resonance dispersion that static homogenisation cannot represent.', p.f_dyn < 0.3 * fStrut),
      ],
      series: [
        { key: 'nu_th', title: 'Poisson\'s ratio versus cell angle', xlabel: 'θ (°)', ylabel: 'ν₁₂', lines: [line('ν₁₂(θ)', ths, ths.map(v => N.clamp(nuOf(v), -20, 20)))], explain: 'The sign of Poisson\'s ratio follows the sign of θ: re-entrant cells (θ < 0) are auxetic. Values are clipped at ±20 near θ = 0, where the model diverges.' },
        { key: 'E_th', title: 'Stiffness versus cell angle at fixed relative density', xlabel: 'θ (°)', ylabel: 'E₁ / E_s', logy: true, lines: [line('E₁/E_s', ths, ths.map(E1of))], explain: 'Relative modulus along x₁ at constant relative density. Auxetic behaviour is obtained at the price of low stiffness per unit mass because the walls bend.' },
      ],
      equations: [
        eqn('formula', 'Gibson–Ashby honeycomb moduli', 'E₁/E_s = (t/l)³ cos θ / [(h/l + sin θ) sin²θ];  E₂/E_s = (t/l)³ (h/l + sin θ)/cos³θ'),
        eqn('formula', 'Poisson\'s ratios', 'ν₁₂ = cos²θ / [(h/l + sin θ) sin θ];  ν₂₁ = (h/l + sin θ) sin θ / cos²θ'),
        eqn('formula', 'Shear modulus and relative density', 'G₁₂/E_s = (t/l)³ (h/l + sin θ) / [(h/l)²(1 + 2h/l) cos θ];  ρ̄ = (t/l)(h/l + 2)/[2 cos θ (h/l + sin θ)]'),
        eqn('formula', 'Scale separation', 'λ = c_T/f ≫ a;  L_struct/a ≳ 10;  f ≪ f_wall'),
      ],
      geometry: { kind: 'lattice', nodes, edges, radius: t * 1e3 / 2, units: 'mm' },
      checks: [
        chk('Reciprocity E₁ν₂₁ = E₂ν₁₂', Math.abs(E1 * nu21 - E2 * nu12) <= 1e-9 * Math.abs(E2 * nu12), `${F(E1 * nu21, 4)} vs ${F(E2 * nu12, 4)} Pa`),
        chk('Orthotropic stability ν₁₂ν₂₁ < 1', stab, `ν₁₂ν₂₁ = ${F(nu12 * nu21, 4)} (equals 1 in the inextensible-wall limit: marginal)`),
        chk('Slender walls', tl <= 0.25, `t/l = ${F(tl, 3)}`),
        chk('Scale separation (dynamic)', sep >= 10 && p.f_dyn < 0.3 * fStrut, `λ/a = ${F(sep, 3)}, f/f_wall = ${F(p.f_dyn / fStrut, 3)}`),
        chk('Scale separation (structural)', nCells >= 10, `${F(nCells, 3)} cells across`),
      ],
      warnings: [...warnings, ...(valid ? [] : ['Homogenisation validity limits are exceeded: use a discrete lattice model (finite-element frame or Bloch analysis) instead of effective constants.']), 'The bending-only model gives ν₁₂ν₂₁ = 1 exactly; including wall stretching and shear lowers the magnitude of ν and restores a strictly positive-definite stiffness.'],
      recommendations: [
        rec('high', aux ? `ν₁₂ = ${F(nu12, 3)} is geometric: verify it on a panel at least 10 cells wide, because clamped edges suppress the auxetic expansion locally.` : 'Use a negative wall angle (for example −20°) if synclastic curvature, shear resistance or indentation resistance through auxetic behaviour is the goal.', 'Gibson–Ashby honeycomb kinematics.'),
        rec('medium', valid ? 'Scale separation holds: the effective constants can be passed to a plate or panel analysis.' : 'Reduce the cell size or the operating frequency, or model the lattice discretely; effective constants are not meaningful under the present conditions.', `λ/a = ${F(sep, 3)}, cells across = ${F(nCells, 3)}.`),
        rec('medium', 'Check wall yielding and elastic buckling of the oblique walls under compression before relying on the linear constants at strains above a few tenths of a per cent.', 'Bending-dominated cells collapse at low strain.'),
        rec('low', `At ${F(100 * rel, 3)} % relative density the lattice uses ${F(100 * (1 - rel), 3)} % less material than a solid panel of the same envelope; choose the lowest density that meets stiffness and strength targets.`, 'Mass efficiency of architected materials.', true),
      ],
      publish: { 'eff.E': E2, 'eff.nu': nu21, 'eff.G': G12, 'eff.density': rhoEff, 'geom.strut': t, 'geom.rel_density': rel, 'geom.cell_size': a },
      numerics: { method: 'Closed-form Gibson–Ashby cellular-solid relations (cell-wall bending)', dof: 3, steps: 1, iterations: 1, residual: 0, notes: 'Slender Euler–Bernoulli walls, rigid joints, small strain; axial and shear wall deformation neglected.' },
    };
  }
  const T = TOPO[p.topo], cn = T.C(rel), c = { C11: cn.C11 * Es, C12: cn.C12 * Es, C44: cn.C44 * Es }, q = cubicFrom(c);
  const r = Math.sqrt(rel / (N.PI * T.Ltot)) * a, l = T.l * a, slender = l / (2 * r);
  const sigY = T.stretch ? 0.3 * rel * p.sy_s : 0.3 * Math.pow(rel, 1.5) * p.sy_s, sigB = T.stretch ? rel / 3 * N.PI ** 2 * Es * (r / l) ** 2 / 4 : 0.05 * rel * rel * Es, strength = Math.min(sigY, sigB);
  const cL = Math.sqrt(c.C11 / rhoEff), cT = Math.sqrt(c.C44 / rhoEff), lamT = cT / p.f_dyn, sep = lamT / a, nCells = p.L_struct / a, fStrut = 22.37 / (2 * N.PI * l * l) * (r / 2) * Math.sqrt(Es / p.rho_m), fBragg = cT / (2 * a);
  const valid = sep >= 10 && nCells >= 10 && p.f_dyn < 0.3 * fStrut, born = c.C11 - c.C12 > 0 && c.C11 + 2 * c.C12 > 0 && c.C44 > 0;
  if (rel > 0.3) warnings.push('Relative density above 0.3: nodes occupy a significant volume and slender-strut scaling laws over-estimate compliance.');
  if (p.topo === 'bcc' || p.topo === 'kelvin') warnings.push('Bending-dominated topology: the prefactor of the ρ̄² law (taken as 1, Gibson–Ashby) depends on strut cross-section and node geometry; calibrate on a unit-cell finite-element model.');
  if (p.topo === 'sc') warnings.push('Simple-cubic lattice is strongly anisotropic: stiff along the struts but shear and off-axis loading are carried by bending (C44 ∝ ρ̄²).');
  const rels = N.logspace(0.002, 0.5, 40), ths = N.linspace(0, 90, 46), g = latticeGeom(p.topo, a * 1e3);
  return {
    summary: `${T.name} lattice at relative density ${F(rel, 3)} (strut radius ${F(r * 1e9, 4)} nm, length ${F(l * 1e9, 4)} nm): C11 = ${F(c.C11 / 1e6, 4)}, C12 = ${F(c.C12 / 1e6, 4)}, C44 = ${F(c.C44 / 1e6, 4)} MPa; E[100] = ${F(q.E100 / 1e6, 4)} MPa, E[111] = ${F(q.E111 / 1e6, 4)} MPa, Zener ratio ${F(q.A, 3)}; strength ≈ ${F(strength / 1e6, 3)} MPa (${sigB < sigY ? 'strut buckling' : 'yield'}-limited). Homogenisation is ${valid ? 'valid' : 'NOT valid'} at ${F(p.f_dyn, 3)} Hz (wavelength/cell = ${F(sep, 3)}, ${F(nCells, 3)} cells across).`,
    kpis: [
      kpi('E100', 'Effective modulus along [100]', q.E100, 'Pa', `${T.stretch ? 'Stretch-dominated affine strut model (stiffness ∝ ρ̄)' : 'Bending-dominated Gibson–Ashby scaling (stiffness ∝ ρ̄²)'}: E/E_s = ${F(q.E100 / Es, 3)} at ρ̄ = ${F(rel, 3)}.`, null),
      kpi('E111', 'Effective modulus along [111]', q.E111, 'Pa', 'Directional modulus along the cube diagonal from the cubic compliance; its ratio to E[100] expresses the elastic anisotropy.', null),
      kpi('E_iso', 'Polycrystal-average modulus (VRH)', q.E, 'Pa', 'Voigt–Reuss–Hill average of the cubic constants: the modulus of a multi-domain or randomly oriented lattice.', null),
      kpi('nu', 'Poisson\'s ratio (VRH)', q.nu, '–', `Isotropic average; the single-crystal value along [100] is ${F(q.nu100, 3)}.`, q.nu > -1 && q.nu < 0.5),
      kpi('G', 'Shear modulus (VRH)', q.G, 'Pa', 'Voigt–Reuss–Hill shear modulus from C44 and (C11 − C12)/2.', null),
      kpi('K', 'Bulk modulus', q.K, 'Pa', '(C11 + 2C12)/3; hydrostatic loading is carried by strut stretching in all the cubic topologies considered.', null),
      kpi('zener', 'Zener anisotropy ratio', q.A, '–', '2C44/(C11 − C12): 1 for isotropy, 2 for the ideal octet truss; values far from 1 mean stiffness depends strongly on loading direction.', q.A > 0.5 && q.A < 2.5),
      kpi('strength', 'Effective strength', strength, 'Pa', `Minimum of plastic collapse (${F(sigY / 1e6, 3)} MPa, Ashby scaling 0.3ρ̄^${T.stretch ? '1' : '1.5'}σ_ys) and elastic strut buckling (${F(sigB / 1e6, 3)} MPa).`, null),
      kpi('rho_eff', 'Effective density', rhoEff, 'kg/m³', 'Relative density times the strut-material density.', null),
      kpi('r_strut', 'Strut radius', r, 'm', `Circular strut radius giving the stated relative density, from ρ̄ = πr²·L_tot/a³ with ${F(T.Ltot, 4)} strut lengths per unit cell; slenderness l/2r = ${F(slender, 3)}.`, slender > 2.5),
      kpi('sep', 'Shear wavelength / cell size', sep, '–', `Shear wave speed ${F(cT, 4)} m/s gives a wavelength of ${F(lamT, 3)} m at ${F(p.f_dyn, 3)} Hz; a ratio of at least 10 is required for quasi-static effective constants.`, sep >= 10),
      kpi('f_limit', 'Homogenisation frequency limit', Math.min(0.1 * cT / a, 0.3 * fStrut), 'Hz', `Lower of λ = 10a (${F(0.1 * cT / a, 3)} Hz) and 30 % of the first strut bending resonance (${F(fStrut, 3)} Hz). The first Bragg gap is expected near ${F(fBragg, 3)} Hz.`, p.f_dyn < Math.min(0.1 * cT / a, 0.3 * fStrut)),
      kpi('n_cells', 'Cells across the part', nCells, '–', 'Smallest structural dimension divided by the cell size; below about 10 the response is size dependent (boundary layers of incompletely constrained struts).', nCells >= 10),
    ],
    series: [
      { key: 'E_rho', title: 'Stiffness scaling with relative density', xlabel: 'relative density', ylabel: 'E[100] / E_s', logx: true, logy: true, lines: Object.entries(TOPO).map(([k, t]) => line(t.name, rels, rels.map(rr => { const cc = t.C(rr); return cubicFrom(cc).E100; }))), explain: 'Stretch-dominated lattices (slope 1) are far stiffer than bending-dominated ones (slope 2) at low density; the selected topology is one of the lines.' },
      { key: 'E_dir', title: 'Directional Young\'s modulus', xlabel: 'angle from [001] towards [110] (°)', ylabel: 'E (MPa)', lines: [line(T.name, ths, ths.map(t => { const s = Math.sin(t * N.PI / 180), cz = Math.cos(t * N.PI / 180), n1 = s / Math.SQRT2; return q.Edir(n1 * n1 * n1 * n1 + 2 * n1 * n1 * cz * cz) / 1e6; }))], explain: 'Modulus for uniaxial loading in the (1−10) plane, passing through [001] (0°), [111] (54.7°) and [110] (90°). A flat curve means an isotropic lattice.' },
    ],
    tables: [{ key: 'cij', title: 'Effective cubic stiffness constants', columns: ['Constant', 'Value (MPa)', 'Normalised by ρ̄E_s'], rows: [['C11', F(c.C11 / 1e6, 5), F(cn.C11 / rel, 4)], ['C12', F(c.C12 / 1e6, 5), F(cn.C12 / rel, 4)], ['C44', F(c.C44 / 1e6, 5), F(cn.C44 / rel, 4)]], explain: 'Homogenised stiffness constants in the cubic lattice axes, ready to be used as an anisotropic continuum material card.' }],
    equations: [
      eqn('constitutive', 'Affine strut (pin-jointed) stiffness', 'C_ijkl = E_s Σ_s ρ̄_s n_i n_j n_k n_l   (stretch-dominated topologies)'),
      eqn('formula', 'Octet truss', 'C11 = ρ̄E_s/6,  C12 = C44 = ρ̄E_s/12  ⇒  E[100] = ρ̄E_s/9,  A = 2'),
      eqn('formula', 'Gibson–Ashby bending scaling', 'E/E_s ≈ ρ̄²,  σ_pl/σ_ys ≈ 0.3ρ̄^1.5;   stretch: σ_pl/σ_ys ≈ 0.3ρ̄'),
      eqn('formula', 'Scale separation', 'λ = c_T/f ≥ 10a;  f ≤ 0.3 f_strut,  f_strut = (22.37/2πl²)√(E_sI/ρ_sA)'),
    ],
    geometry: { kind: 'lattice', nodes: g.nodes, edges: g.edges, radius: r * 1e3, units: 'mm' },
    checks: [
      chk('Born stability of effective constants', born, `C11 − C12 = ${F((c.C11 - c.C12) / 1e6, 4)} MPa, C44 = ${F(c.C44 / 1e6, 4)} MPa`),
      chk('Slender-strut assumption', rel <= 0.3 && slender > 2.5, `ρ̄ = ${F(rel, 3)}, l/2r = ${F(slender, 3)}`),
      chk('Scale separation (dynamic)', sep >= 10 && p.f_dyn < 0.3 * fStrut, `λ/a = ${F(sep, 3)}, f/f_strut = ${F(p.f_dyn / fStrut, 3)}`),
      chk('Scale separation (structural)', nCells >= 10, `${F(nCells, 3)} cells across the smallest dimension`),
      chk('Effective modulus below Voigt bound', q.E100 <= rel * Es * (1 + 1e-9), `E[100]/(ρ̄E_s) = ${F(q.E100 / (rel * Es), 3)}`),
    ],
    warnings: [...warnings, ...(valid ? [] : ['Homogenisation validity limits are exceeded for the stated frequency or part size: use discrete lattice or Bloch-wave analysis rather than effective constants.'])],
    recommendations: [
      rec('high', valid ? 'Scale separation is satisfied: export the effective constants and density as a continuum material card for component-level analysis.' : `Effective constants are not valid here. Keep the frequency below ${F(Math.min(0.1 * cT / a, 0.3 * fStrut), 3)} Hz and at least 10 cells across the part, or model the lattice explicitly.`, `λ/a = ${F(sep, 3)}, cells across = ${F(nCells, 3)}.`),
      rec('high', T.stretch ? 'Stretch-dominated topology: stiffness and strength scale linearly with density, the efficient choice for load-bearing parts. Check strut buckling at low density.' : 'Bending-dominated topology: good for energy absorption and compliance but mass-inefficient for stiffness; switch to an octet truss if stiffness per mass matters.', 'Maxwell rigidity criterion and Gibson–Ashby scaling.'),
      rec('medium', sigB < sigY ? `Failure is governed by elastic strut buckling (${F(sigB / 1e6, 3)} MPa): raise the relative density or use hollow or hierarchical struts to move into the yield-limited regime.` : `Failure is yield-limited (${F(sigY / 1e6, 3)} MPa); nanoscale struts may exceed the bulk strength used here through size effects — measure it.`, 'Competing collapse mechanisms.'),
      rec('low', `The lattice reaches E[100] = ${F(q.E100 / 1e6, 3)} MPa with ${F(100 * rel, 3)} % of the solid mass; minimising relative density for the required stiffness directly reduces raw-material and print-energy demand.`, 'Mass efficiency of architected materials.', true),
    ],
    publish: { 'eff.E': q.E100, 'eff.nu': q.nu100, 'eff.G': c.C44, 'eff.K': q.K, 'eff.density': rhoEff, 'geom.strut': 2 * r, 'geom.rel_density': rel, 'geom.cell_size': a, 'prop.c_sound': cL },
    numerics: { method: 'Closed-form affine-strut stiffness tensors and Gibson–Ashby scaling laws; cubic compliance algebra', dof: 3, steps: 1, iterations: 1, residual: 0, notes: 'Pin-jointed affine deformation for stretching; unit prefactor for bending-dominated cells; node volume neglected.' },
  };
}

// ---------- single crystal → polycrystal ----------
function ljConstants(epsJ, sig) { // static fcc lattice sums for the 12-6 potential
  const sh = []; for (let i = -9; i <= 9; i++) for (let j = -9; j <= 9; j++) for (let k = -9; k <= 9; k++) { if ((i + j + k) & 1 || (i === 0 && j === 0 && k === 0)) continue; const r2 = (i * i + j * j + k * k) / 2; if (r2 <= 36) sh.push([i / Math.SQRT2, j / Math.SQRT2, k / Math.SQRT2, r2]); } // in units of nearest-neighbour distance d
  const U = d => { let s = 0; for (const v of sh) { const x = sig * sig / (v[3] * d * d), x3 = x * x * x; s += 4 * epsJ * (x3 * x3 - x3); } return s / 2; };
  const d = N.goldenSection(U, 0.95 * sig, 1.3 * sig, 1e-9 * sig), Om = d ** 3 / Math.SQRT2; let c11 = 0, c12 = 0, pv = 0;
  for (const v of sh) { const r = Math.sqrt(v[3]) * d, x = (sig / r) ** 6, d1 = 4 * epsJ * (-12 * x * x + 6 * x) / r, d2 = 4 * epsJ * (156 * x * x - 42 * x) / (r * r), w = (d2 - d1 / r) / (r * r), X = v[0] * d, Y = v[1] * d; c11 += w * X ** 4; c12 += w * X * X * Y * Y; pv += d1 * r; }
  return { d, a: Math.SQRT2 * d, Om, C11: c11 / (2 * Om), C12: c12 / (2 * Om), C44: c12 / (2 * Om), Ecoh: -U(d), pressure: -pv / (6 * Om) };
}
function atomistic(p) {
  const warnings = []; let formula = p.formula, pf;
  try { pf = chem.parseFormula(formula); if (!Object.keys(pf.atoms).length) throw new Error('empty formula'); } catch (e) { warnings.push(`Formula "${formula}" could not be parsed (${e.message}); silicon is used for molar mass.`); formula = 'Si'; pf = chem.parseFormula('Si'); }
  const Mmol = chem.molarMass(formula) * 1e-3, nAt = N.sum(Object.values(pf.atoms));
  let sym = p.sym, C11 = p.C11 * 1e9, C12 = p.C12 * 1e9, C44 = p.C44 * 1e9, C13 = p.C13 * 1e9, C33 = p.C33 * 1e9, C66 = p.C66 * 1e9, rho = p.rho_x, lj = null;
  if (p.src === 'lj') { lj = ljConstants(p.eps_lj * N.Q_E, p.sig_lj * 1e-10); sym = 'cubic'; C11 = lj.C11; C12 = lj.C12; C44 = lj.C44; rho = Mmol / nAt / (N.N_A * lj.Om) * 1; warnings.push('Lennard-Jones constants are static-lattice (0 K, no zero-point motion) values for an fcc crystal; central pair forces impose the Cauchy relation C12 = C44, which metals and covalent crystals violate.'); }
  const C = N.zeros(6, 6);
  if (sym === 'cubic') { for (let i = 0; i < 3; i++) { C[i][i] = C11; C[i + 3][i + 3] = C44; for (let j = 0; j < 3; j++) if (i !== j) C[i][j] = C12; } }
  else { if (sym === 'hex') C66 = (C11 - C12) / 2; C[0][0] = C[1][1] = C11; C[2][2] = C33; C[0][1] = C[1][0] = C12; C[0][2] = C[2][0] = C[1][2] = C[2][1] = C13; C[3][3] = C[4][4] = C44; C[5][5] = C66; }
  const ev = N.eigSym(C).values, stable = ev[0] > 0;
  let S; try { S = N.inverse(C); } catch { S = null; }
  const KV = (C[0][0] + C[1][1] + C[2][2] + 2 * (C[0][1] + C[0][2] + C[1][2])) / 9, GV = (C[0][0] + C[1][1] + C[2][2] - (C[0][1] + C[0][2] + C[1][2]) + 3 * (C[3][3] + C[4][4] + C[5][5])) / 15;
  let KR = KV, GR = GV;
  if (S && stable) { KR = 1 / (S[0][0] + S[1][1] + S[2][2] + 2 * (S[0][1] + S[0][2] + S[1][2])); GR = 15 / (4 * (S[0][0] + S[1][1] + S[2][2]) - 4 * (S[0][1] + S[0][2] + S[1][2]) + 3 * (S[3][3] + S[4][4] + S[5][5])); }
  else warnings.push('The stiffness matrix is not positive definite (Born-unstable): Reuss averages cannot be formed and the Voigt values are reported instead. Check the entered constants.');
  const K = 0.5 * (KV + KR), G = 0.5 * (GV + GR), E = Eof(K, G), nu = nuof(K, G), AU = stable ? 5 * GV / GR + KV / KR - 6 : 0;
  const Gp = Math.max(G, 1e3), Kp = Math.max(K, 1e3), vl = Math.sqrt((Kp + 4 * Gp / 3) / rho), vt = Math.sqrt(Gp / rho), vm = Math.pow((2 / vt ** 3 + 1 / vl ** 3) / 3, -1 / 3), thD = N.H_PLANCK / N.K_B * Math.cbrt(3 * nAt * N.N_A * rho / (4 * N.PI * Mmol)) * vm;
  const Edir = n => { if (!S) return E; const [a, b, c] = n; return 1 / (S[0][0] * a ** 4 + S[1][1] * b ** 4 + S[2][2] * c ** 4 + (2 * S[0][1] + S[5][5]) * a * a * b * b + (2 * S[0][2] + S[4][4]) * a * a * c * c + (2 * S[1][2] + S[3][3]) * b * b * c * c); };
  const ths = N.linspace(0, 90, 46), dirOf = t => { const s = Math.sin(t * N.PI / 180), c = Math.cos(t * N.PI / 180); return sym === 'cubic' ? [s / Math.SQRT2, s / Math.SQRT2, c] : [s, 0, c]; }, Eang = ths.map(t => Edir(dirOf(t)) / 1e9);
  const Emax = N.maxOf(Eang), Emin = N.minOf(Eang), zen = sym === 'cubic' ? 2 * C44 / (C11 - C12) : C[3][3] / C[5][5];
  const Ts = N.linspace(5, Math.max(600, 1.5 * thD), 80), cv = T => { const xm = thD / T, xs = N.linspace(1e-6, xm, 120), ys = xs.map(x => { const ex = Math.exp(Math.min(x, 500)); return x ** 4 * ex / (ex - 1) ** 2; }); return 9 * nAt * N.R_GAS * (T / thD) ** 3 * N.trapz(xs, ys); };
  const born = sym === 'cubic' ? `C11 − C12 = ${GPa(C11 - C12)} > 0, C11 + 2C12 = ${GPa(C11 + 2 * C12)} > 0, C44 = ${GPa(C44)} > 0` : `C44 = ${GPa(C44)} > 0, C11 − |C12| = ${GPa(C11 - Math.abs(C12))} > 0, (C11 + C12)C33 − 2C13² = ${F(((C11 + C12) * C33 - 2 * C13 * C13) / 1e18, 4)} GPa² > 0${sym === 'tet' ? `, C66 = ${GPa(C66)} > 0` : ''}`;
  const pretty = chem.prettyFormula(formula), pugh = K / Gp;
  if (AU > 1.5) warnings.push(`Strong elastic anisotropy (A^U = ${F(AU, 3)}): Voigt and Reuss bounds are far apart and textured or nanostructured material may deviate markedly from the Hill average.`);
  return {
    summary: `${pretty} (${sym === 'cubic' ? 'cubic' : sym === 'hex' ? 'hexagonal' : 'tetragonal'}${lj ? ', Lennard-Jones fcc' : ''}): Voigt–Reuss–Hill polycrystal averages E = ${GPa(E)} GPa, G = ${GPa(G)} GPa, K = ${GPa(K)} GPa, ν = ${F(nu, 3)}; longitudinal and shear wave speeds ${F(vl, 4)} and ${F(vt, 4)} m/s; Debye temperature ${F(thD, 4)} K. The crystal is ${stable ? 'mechanically stable (Born criteria met)' : 'mechanically UNSTABLE for the entered constants'} with universal anisotropy index ${F(AU, 3)}.`,
    kpis: [
      kpi('E', 'Polycrystal Young\'s modulus (VRH)', E, 'Pa', `Hill average of the Voigt (uniform strain) and Reuss (uniform stress) polycrystal moduli: E = 9KG/(3K + G). Single-crystal values range from ${F(Emin, 4)} to ${F(Emax, 4)} GPa along the sampled directions.`, stable),
      kpi('G', 'Shear modulus (VRH)', G, 'Pa', `Mean of Voigt ${GPa(GV)} GPa and Reuss ${GPa(GR)} GPa; the gap between them measures the uncertainty from grain interaction.`, stable),
      kpi('K', 'Bulk modulus (VRH)', K, 'Pa', `Mean of Voigt ${GPa(KV)} GPa and Reuss ${GPa(KR)} GPa (identical for cubic crystals).`, stable),
      kpi('nu', 'Poisson\'s ratio', nu, '–', 'From the Hill-averaged bulk and shear moduli, ν = (3K − 2G)/(2(3K + G)).', nu > -1 && nu < 0.5),
      kpi('v_l', 'Longitudinal wave speed', vl, 'm/s', 'Compressional wave speed of the isotropic polycrystalline aggregate, √((K + 4G/3)/ρ).', null),
      kpi('v_t', 'Shear wave speed', vt, 'm/s', 'Shear (transverse) wave speed of the isotropic polycrystalline aggregate, √(G/ρ).', null),
      kpi('theta_D', 'Debye temperature', thD, 'K', `θ_D = (h/k_B)·(3nN_Aρ/4πM)^(1/3)·v_m with mean sound speed v_m = ${F(vm, 4)} m/s, ${nAt} atom(s) per formula unit and M = ${F(Mmol * 1e3, 5)} g/mol.`, null),
      kpi('A_U', 'Universal anisotropy index', AU, '–', 'A^U = 5G_V/G_R + K_V/K_R − 6: zero for an elastically isotropic crystal, increasing with anisotropy for any symmetry.', AU < 1.5),
      kpi('zener', sym === 'cubic' ? 'Zener ratio 2C44/(C11 − C12)' : 'Shear anisotropy C44/C66', zen, '–', 'Ratio of the two principal shear stiffnesses; 1 means the shear response is isotropic.', null),
      kpi('pugh', 'Pugh ratio K/G', pugh, '–', 'Empirical ductility indicator: values above about 1.75 suggest ductile behaviour, below it brittle behaviour.', null),
      kpi('C_min', 'Smallest stiffness eigenvalue', ev[0], 'Pa', 'Lowest eigenvalue of the 6×6 stiffness matrix; it must be positive for mechanical (Born) stability.', stable),
      kpi('rho', 'Density used', rho, 'kg/m³', lj ? `Computed from the equilibrium fcc lattice parameter a = ${F(lj.a * 1e10, 4)} Å of the pair potential.` : 'Entered crystal density, used for wave speeds and the Debye temperature.', null),
    ],
    series: [
      { key: 'E_dir', title: 'Directional Young\'s modulus of the single crystal', xlabel: sym === 'cubic' ? 'angle from [001] towards [110] (°)' : 'angle from the c axis towards the a axis (°)', ylabel: 'E (GPa)', lines: [line('single crystal', ths, Eang), line('VRH polycrystal', ths, ths.map(() => E / 1e9))], explain: 'Young\'s modulus for uniaxial stress along directions in a symmetry plane, compared with the isotropic polycrystal average. Nanowires and textured films take the single-crystal value of their growth axis.' },
      { key: 'cv', title: 'Debye heat capacity', xlabel: 'temperature (K)', ylabel: 'C_v (J/mol/K)', lines: [line('Debye model', Ts, Ts.map(cv)), line('Dulong–Petit limit', Ts, Ts.map(() => 3 * nAt * N.R_GAS))], explain: 'Lattice heat capacity per mole of formula units from the Debye model with the elastically derived Debye temperature; it approaches 3nR at high temperature.' },
    ],
    tables: [{ key: 'vrh', title: 'Voigt, Reuss and Hill averages (GPa)', columns: ['Quantity', 'Voigt', 'Reuss', 'Hill'], rows: [['K', GPa(KV), GPa(KR), GPa(K)], ['G', GPa(GV), GPa(GR), GPa(G)], ['E', GPa(Eof(KV, GV)), GPa(Eof(KR, GR)), GPa(E)], ['ν', F(nuof(KV, GV), 4), F(nuof(KR, GR), 4), F(nu, 4)]], explain: 'Upper (Voigt) and lower (Reuss) polycrystal bounds and their arithmetic (Hill) mean, the usual estimate for an untextured aggregate.' }],
    equations: [
      eqn('formula', 'Formula unit', `${pretty}:  M = ${F(Mmol * 1e3, 5)} g/mol, ${nAt} atom(s) per formula unit`),
      eqn('formula', 'Voigt averages', '9K_V = (C11+C22+C33) + 2(C12+C13+C23);  15G_V = (C11+C22+C33) − (C12+C13+C23) + 3(C44+C55+C66)'),
      eqn('formula', 'Reuss averages', '1/K_R = (S11+S22+S33) + 2(S12+S13+S23);  15/G_R = 4(S11+S22+S33) − 4(S12+S13+S23) + 3(S44+S55+S66)'),
      eqn('formula', 'Hill average and isotropic constants', 'K = ½(K_V + K_R), G = ½(G_V + G_R);  E = 9KG/(3K+G), ν = (3K−2G)/(2(3K+G))'),
      eqn('formula', 'Debye temperature', 'θ_D = (h/k_B)·[3nN_Aρ/(4πM)]^(1/3)·v_m,  v_m = [⅓(2/v_t³ + 1/v_l³)]^(−1/3)'),
      eqn('formula', 'Born stability', born),
      ...(lj ? [eqn('constitutive', 'Pair-potential lattice sums', 'φ(r) = 4ε[(σ/r)¹² − (σ/r)⁶];  C_ijkl = (1/2Ω) Σ_R (φ″ − φ′/r) R_iR_jR_kR_l / r²', `d = ${F(lj.d / (p.sig_lj * 1e-10), 5)}σ, cohesive energy ${F(lj.Ecoh / (p.eps_lj * N.Q_E), 4)}ε per atom`)] : []),
    ],
    checks: [
      chk('Born mechanical stability', stable, born),
      chk('Reuss ≤ Hill ≤ Voigt', KR <= KV * (1 + 1e-9) && GR <= GV * (1 + 1e-9), `G: ${GPa(GR)} ≤ ${GPa(G)} ≤ ${GPa(GV)} GPa`),
      chk('Poisson\'s ratio admissible', nu > -1 && nu < 0.5, `ν = ${F(nu, 4)}`),
      ...(lj ? [chk('Zero-pressure equilibrium of the lattice sum', Math.abs(lj.pressure) < 1e-4 * lj.C11, `residual pressure ${F(lj.pressure, 3)} Pa`), chk('Cauchy relation C12 = C44', Math.abs(lj.C12 - lj.C44) < 1e-9 * lj.C11, 'holds identically for central pair forces')] : []),
    ],
    warnings,
    recommendations: [
      rec('high', AU > 0.5 ? `Anisotropy is significant (E from ${F(Emin, 4)} to ${F(Emax, 4)} GPa): for nanowires, nanorods and textured films use the modulus along the actual growth axis, not the polycrystal average.` : 'The crystal is nearly isotropic: the Hill average can be used for single-crystalline nanostructures as well.', `Universal anisotropy index ${F(AU, 3)}.`),
      rec('medium', 'Pass the Hill-averaged E, ν and density to the continuum suites; keep the Voigt–Reuss gap as the model uncertainty for untextured polycrystals.', 'Voigt and Reuss values are rigorous bounds for a random aggregate.'),
      rec('medium', 'Elastic constants from DFT are typically within 5–10 % of experiment; validate one constant (for example by nanoindentation or resonant ultrasound) before design use, and expect surface-stress effects below about 20 nm.', 'Typical accuracy of first-principles elastic constants and size effects.'),
      rec('low', 'Deriving continuum properties from a small set of atomistic constants avoids repeating expensive simulations for each geometry, reducing computing energy for design exploration.', 'Hierarchical property transfer.', true),
    ],
    publish: { 'prop.E': E, 'prop.nu': nu, 'prop.G': G, 'prop.K': K, 'prop.c_sound': vl, 'prop.density': rho, 'prop.debye_T': thD, 'material.formula': formula, 'prop.molar_mass': Mmol },
    numerics: { method: lj ? 'Static fcc lattice sums of the 12-6 potential (golden-section energy minimisation) + Voigt–Reuss–Hill averaging' : 'Voigt–Reuss–Hill averaging from the 6×6 stiffness matrix (direct inversion, Jacobi eigenvalues)', dof: 36, steps: 1, iterations: 1, residual: lj ? Math.abs(lj.pressure) : 0, notes: 'Athermal elastic constants; no texture, grain-boundary or surface elasticity.' },
  };
}

export default {
  id: 18,
  tagline: 'Bridges scales: composite bounds and mean-field estimates, numerical RVE homogenisation, laminates, lattice effective constants and single-crystal to polycrystal transfer.',
  inputs,
  presets: [
    { name: 'ZnO nanoparticles in epoxy (20 vol %)', note: 'Reference nanocomposite: stiff, thermally conductive oxide particles in a polymer matrix.', values: { mode: 'bounds', shape: 'sphere', vf: 0.2 } },
    { name: 'Graphene nanoplatelet / polymer nanocomposite', note: 'High-aspect-ratio platelets: strong stiffening and a percolation threshold near 0.1 vol %.', values: { mode: 'bounds', shape: 'platelet', aspect: 1000, orient: 'random3d', vf: 0.02, E_m: 3e9, nu_m: 0.35, k_m: 0.2, sig_m: 1e-12, cte_m: 6e-5, rho_m: 1200, E_i: 1e12, nu_i: 0.19, k_i: 3000, sig_i: 1e6, cte_i: -1e-6, rho_i: 2200 } },
    { name: 'ZrW₂O₈–Cu zero-expansion composite', note: 'Negative-thermal-expansion filler in copper: composite expansion can be tuned through zero.', values: { mode: 'bounds', shape: 'sphere', vf: 0.6, E_m: 1.2e11, nu_m: 0.34, k_m: 400, sig_m: 5.8e7, cte_m: 1.7e-5, rho_m: 8960, E_i: 8.8e10, nu_i: 0.3, k_i: 0.8, sig_i: 1e-10, cte_i: -8.7e-6, rho_i: 5080 } },
    { name: 'Random-disc RVE, oxide fibres in polymer', note: 'Numerical homogenisation of a fibre cross-section with conductivity contrast 250 and stiffness contrast 47.', values: { mode: 'rve', micro: 'discs', vf: 0.3, n_incl: 24, nres: '48' } },
    { name: 'W / Al₂O₃ nanolaminate thermal barrier', note: 'Metal–oxide nanolaminate: interface resistance pushes cross-plane conductivity below both constituents.', values: { mode: 'laminate', E_m: 4.1e11, nu_m: 0.28, k_m: 170, sig_m: 1.8e7, cte_m: 4.5e-6, rho_m: 19300, E_i: 3.7e11, nu_i: 0.22, k_i: 30, sig_i: 1e-12, cte_i: 8e-6, rho_i: 3950, vf: 0.5, period_nm: 6, n_bi: 50, R_int: 3e-9, dT: 500 } },
    { name: 'Octet nanolattice (polymer, 10 % density)', note: 'Two-photon-printed stretch-dominated lattice metamaterial with a 2 µm cell.', values: { mode: 'lattice', topo: 'octet', rel: 0.1, E_m: 3e9, nu_m: 0.4, rho_m: 1200, sy_s: 7e7, cell: 2e-6, L_struct: 1e-4, f_dyn: 1e6 } },
    { name: 'Re-entrant auxetic honeycomb', note: 'Negative Poisson\'s ratio from re-entrant cell geometry (θ = −20°).', values: { mode: 'lattice', topo: 'reentrant', theta: -20, h_l: 2, rel: 0.08, E_m: 7e10, nu_m: 0.33, rho_m: 2700, sy_s: 2.5e8, cell: 5e-3, L_struct: 0.1, f_dyn: 500 } },
    { name: 'Silicon single crystal', note: 'Cubic Si: C11 = 165.7, C12 = 63.9, C44 = 79.6 GPa.', values: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 165.7, C12: 63.9, C44: 79.6, rho_x: 2329, formula: 'Si' } },
    { name: 'ZnO wurtzite', note: 'Hexagonal ZnO: C11 = 209.7, C12 = 121.1, C13 = 105.1, C33 = 210.9, C44 = 42.5 GPa.', values: { mode: 'atomistic', src: 'cij', sym: 'hex', C11: 209.7, C12: 121.1, C13: 105.1, C33: 210.9, C44: 42.5, rho_x: 5606, formula: 'ZnO' } },
    { name: 'Gold', note: 'Cubic Au, strongly anisotropic: C11 = 192.9, C12 = 163.8, C44 = 41.5 GPa.', values: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 192.9, C12: 163.8, C44: 41.5, rho_x: 19300, formula: 'Au' } },
    { name: 'TiO₂ rutile', note: 'Tetragonal rutile: C11 = 269, C12 = 177, C13 = 146, C33 = 480, C44 = 124, C66 = 192 GPa.', values: { mode: 'atomistic', src: 'cij', sym: 'tet', C11: 269, C12: 177, C13: 146, C33: 480, C44: 124, C66: 192, rho_x: 4250, formula: 'TiO2' } },
  ],
  run(inp) {
    const p = read(inp);
    const r = p.mode === 'rve' ? rve(p) : p.mode === 'laminate' ? laminate(p) : p.mode === 'lattice' ? lattice(p) : p.mode === 'atomistic' ? atomistic(p) : bounds(p);
    return finish(r);
  },
};
