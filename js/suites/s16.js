// Suite 16 — Microstructure, phase-field and defect evolution.
// Reduced-order engineering models: 2-D spectral Cahn–Hilliard, Monte-Carlo Potts grain growth, Kobayashi phase-field
// solidification, analytical grain-growth / coarsening / sintering / JMAK kinetics and point/line/planar defect equilibria.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => note ? { kind, label, text, note } : { kind, label, text };
const rec = (priority, text, basis, sustainability) => sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis };
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const line = (name, x, y) => { const X = [], Y = [], st = Math.max(1, Math.ceil(x.length / 400)); for (let i = 0; i < x.length; i += st) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) { X.push(x[i]); Y.push(y[i]); } if (!X.length) { X.push(0); Y.push(0); } return { name, x: X, y: Y }; };
const KB_EV = N.K_B / N.Q_E;

const MODES = [
  { value: 'spinodal', label: 'Spinodal decomposition (Cahn–Hilliard, 2-D spectral)' },
  { value: 'grain', label: 'Grain growth with particle pinning (Monte-Carlo Potts)' },
  { value: 'solid', label: 'Dendritic crystal growth (Kobayashi phase field)' },
  { value: 'kinetics', label: 'Analytical kinetics: grain growth, LSW, sintering, JMAK' },
  { value: 'defects', label: 'Defect equilibria: vacancies, Kröger–Vink, dislocations, segregation' },
];
const OXIDES = { ZnO: 5606, TiO2: 4230, CeO2: 7220, Al2O3: 3980, ZrO2: 5680, MgO: 3580 };

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'kinetics', options: MODES, group: 'Model', help: 'Selects the microstructure or defect model. Each analysis uses its own inputs shown below.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Model', showIf: { mode: ['spinodal', 'grain', 'solid'] }, help: 'Seed of the initial fluctuation field or Monte-Carlo sequence; the same seed reproduces the same microstructure.' },
  { key: 'dx_nm', label: 'Grid spacing', unit: 'nm', value: 1, min: 0.1, max: 1e4, log: true, group: 'Model', showIf: { mode: ['spinodal', 'grain'] }, help: 'Physical size of one grid cell. It converts simulated lengths (domain size, grain size) to nanometres.' },
  // spinodal
  { key: 'c0', label: 'Mean composition (atomic fraction of B)', value: 0.5, min: 0.05, max: 0.95, step: 0.01, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Average composition of the alloy or blend. Inside the spinodal (about 0.21–0.79 for this free energy) the mixture decomposes spontaneously; outside it is metastable.' },
  { key: 'noise', label: 'Initial fluctuation amplitude', value: 0.02, min: 1e-4, max: 0.2, log: true, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Amplitude of the seeded random composition fluctuation that triggers decomposition.' },
  { key: 'kappa', label: 'Gradient-energy coefficient κ/(f₀Δx²)', value: 1.5, min: 0.5, max: 6, step: 0.1, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Dimensionless gradient-energy coefficient. Larger values widen the interface and coarsen the initial spinodal wavelength.' },
  { key: 'f0', label: 'Free-energy density scale f₀', unit: 'J/m³', value: 1e8, min: 1e6, max: 1e10, log: true, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Height scale of the double-well free energy; with κ and the grid spacing it fixes the interfacial energy.' },
  { key: 'Dint', label: 'Interdiffusion coefficient magnitude', unit: 'm²/s', value: 1e-18, min: 1e-24, max: 1e-9, log: true, group: 'Cahn–Hilliard', bus: 'prop.D_diff', showIf: { mode: ['spinodal'] }, help: 'Sets the physical time unit Δx²/D of the dimensionless simulation.' },
  { key: 'ch_steps', label: 'Time steps', value: 160, min: 40, max: 400, step: 10, sens: false, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Number of semi-implicit spectral steps. More steps follow coarsening further.' },
  { key: 'ch_dt', label: 'Dimensionless time step', value: 1, min: 0.05, max: 2, step: 0.05, sens: false, group: 'Cahn–Hilliard', showIf: { mode: ['spinodal'] }, help: 'Time step in units of Δx²/D. The stabilised scheme is energy-stable, but accuracy degrades for large steps.' },
  // grain growth
  { key: 'Qstates', label: 'Number of orientations Q', value: 48, min: 8, max: 200, step: 1, group: 'Potts grain growth', showIf: { mode: ['grain'] }, help: 'Number of distinct grain orientations. Large Q suppresses artificial coalescence of like-oriented grains.' },
  { key: 'mcs', label: 'Monte-Carlo steps', value: 80, min: 20, max: 500, step: 10, sens: false, group: 'Potts grain growth', showIf: { mode: ['grain'] }, help: 'Number of Monte-Carlo steps (one attempted re-orientation per site per step) at zero simulation temperature.' },
  { key: 'fp', label: 'Pinning-particle fraction', value: 0.01, min: 0, max: 0.1, step: 0.002, group: 'Potts grain growth', showIf: { mode: ['grain', 'kinetics'] }, help: 'Volume (area) fraction of inert second-phase particles. They pin boundaries and set the Zener limiting grain size.' },
  // solidification
  { key: 'latent', label: 'Dimensionless latent heat K', value: 1.6, min: 0.8, max: 2.4, step: 0.05, group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Kobayashi latent-heat parameter; the dimensionless undercooling is 1/K. Larger K slows growth and promotes side-branching.' },
  { key: 'aniso', label: 'Interface anisotropy strength δ', value: 0.02, min: 0, max: 0.05, step: 0.005, group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Strength of the orientation dependence of the interface energy; zero gives a compact, seaweed-like crystal.' },
  { key: 'fold', label: 'Crystal symmetry', type: 'select', value: '6', options: [{ value: '4', label: 'Four-fold (cubic)' }, { value: '6', label: 'Six-fold (hexagonal, e.g. ZnO basal plane, ice)' }], group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Rotational symmetry of the anisotropy, which fixes the number of primary dendrite arms.' },
  { key: 'sol_steps', label: 'Time steps', value: 300, min: 100, max: 500, step: 20, sens: false, group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Number of explicit steps of 2×10⁻⁴ dimensionless time units.' },
  { key: 'dx_um', label: 'Grid spacing', unit: 'µm', value: 1, min: 0.01, max: 1000, log: true, group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Physical size of one grid cell of the dendrite simulation. The model length unit is Δx/0.03, so this choice fixes the diffusion length scale.' },
  { key: 'alpha_th', label: 'Thermal (or solute) diffusivity', unit: 'm²/s', value: 1e-9, min: 1e-12, max: 1e-3, log: true, group: 'Solidification', showIf: { mode: ['solid'] }, help: 'Diffusivity of the field that limits growth; converts the dimensionless tip velocity to m/s.' },
  // kinetics + shared thermal inputs
  { key: 'T_C', label: 'Temperature', unit: '°C', value: 500, min: 25, max: 1800, group: 'Conditions', bus: 'proc.T_C', showIf: { mode: ['grain', 'kinetics', 'defects'] }, help: 'Annealing, sintering or equilibration temperature. All rate constants and defect concentrations are Arrhenius functions of it.' },
  { key: 'time_h', label: 'Hold time', unit: 'h', value: 1, min: 0.001, max: 1000, log: true, group: 'Conditions', showIf: { mode: ['kinetics'] }, help: 'Isothermal hold time over which grain growth, coarsening, sintering and transformation are evaluated.' },
  { key: 'd0_nm', label: 'Initial grain size', unit: 'nm', value: 20, min: 1, max: 1e5, log: true, group: 'Grain growth', showIf: { mode: ['kinetics'] }, help: 'Starting mean grain diameter (e.g. as-synthesised nanocrystallite size).' },
  { key: 'M0', label: 'Boundary mobility pre-factor M₀', unit: 'm⁴/(J·s)', value: 1e-5, min: 1e-12, max: 1e2, log: true, group: 'Grain growth', showIf: { mode: ['grain', 'kinetics'] }, help: 'Pre-exponential grain-boundary mobility; M = M₀·exp(−Q/RT). Must be calibrated on time-resolved grain-size data.' },
  { key: 'Qgb', label: 'Boundary-migration activation energy', unit: 'kJ/mol', value: 200, min: 20, max: 700, group: 'Grain growth', showIf: { mode: ['grain', 'kinetics'] }, help: 'Activation energy of grain-boundary migration.' },
  { key: 'gamma_gb', label: 'Grain-boundary energy', unit: 'J/m²', value: 0.5, min: 0.05, max: 3, group: 'Grain growth', bus: 'prop.surface_energy', showIf: { mode: ['grain', 'kinetics'] }, help: 'Driving force of curvature-driven growth; also used as surface energy in the sintering model.' },
  { key: 'rp_nm', label: 'Pinning-particle radius', unit: 'nm', value: 10, min: 0.5, max: 5000, log: true, group: 'Grain growth', showIf: { mode: ['kinetics'] }, help: 'Radius of second-phase particles in the Zener–Smith limit D = 4r/(3f).' },
  { key: 'D0', label: 'Diffusivity pre-factor D₀', unit: 'm²/s', value: 1e-5, min: 1e-10, max: 1, log: true, group: 'Diffusion', showIf: { mode: ['kinetics'] }, help: 'Pre-exponential factor of the diffusivity that controls coarsening and sintering.' },
  { key: 'Qd', label: 'Diffusion activation energy', unit: 'kJ/mol', value: 180, min: 20, max: 700, group: 'Diffusion', showIf: { mode: ['kinetics'] }, help: 'Activation energy of the rate-controlling diffusion path.' },
  { key: 'c_sol', label: 'Equilibrium solubility', unit: 'mol/m³', value: 100, min: 1e-3, max: 1e5, log: true, group: 'Coarsening', showIf: { mode: ['kinetics'] }, help: 'Solubility of the precipitate species in the matrix (flat-interface value) for LSW coarsening.' },
  { key: 'Vm', label: 'Molar volume', unit: 'cm³/mol', value: 14.5, min: 3, max: 200, group: 'Coarsening', showIf: { mode: ['kinetics'] }, help: 'Molar volume of the coarsening or sintering phase (ZnO ≈ 14.5 cm³/mol).' },
  { key: 'gamma_int', label: 'Precipitate–matrix interface energy', unit: 'J/m²', value: 0.3, min: 0.01, max: 3, group: 'Coarsening', showIf: { mode: ['kinetics'] }, help: 'Interfacial energy driving Ostwald ripening.' },
  { key: 'r0_nm', label: 'Initial precipitate / particle radius', unit: 'nm', value: 5, min: 0.5, max: 1e4, log: true, group: 'Coarsening', showIf: { mode: ['kinetics'] }, help: 'Initial mean precipitate radius for LSW coarsening.' },
  { key: 'sinter', label: 'Sintering mechanism', type: 'select', value: 'gb', options: [{ value: 'vol', label: 'Volume (lattice) diffusion' }, { value: 'gb', label: 'Grain-boundary diffusion' }, { value: 'surf', label: 'Surface diffusion (no densification)' }], group: 'Sintering', showIf: { mode: ['kinetics'] }, help: 'Dominant initial-stage neck-growth mechanism. The diffusivity above is interpreted as that of the selected path.' },
  { key: 'dp_nm', label: 'Powder particle diameter', unit: 'nm', value: 200, min: 5, max: 1e5, log: true, group: 'Sintering', showIf: { mode: ['kinetics'] }, help: 'Particle diameter of the powder compact; neck growth is strongly size dependent (Herring scaling).' },
  { key: 'rho_green', label: 'Green relative density', value: 0.55, min: 0.3, max: 0.74, step: 0.01, group: 'Sintering', showIf: { mode: ['kinetics'] }, help: 'Relative density of the compact before sintering.' },
  { key: 'n_av', label: 'Avrami exponent', value: 2.5, min: 0.5, max: 4, step: 0.1, group: 'Transformation', showIf: { mode: ['kinetics'] }, help: 'JMAK exponent: about 3–4 for 3-D growth with nucleation, lower for site-saturated or diffusion-controlled growth.' },
  { key: 'k0_av', label: 'JMAK rate pre-factor', unit: '1/s', value: 1e8, min: 1, max: 1e16, log: true, group: 'Transformation', showIf: { mode: ['kinetics'] }, help: 'Pre-exponential factor of the JMAK rate constant k = k₀·exp(−Q/RT).' },
  { key: 'Q_av', label: 'Transformation activation energy', unit: 'kJ/mol', value: 180, min: 20, max: 700, group: 'Transformation', showIf: { mode: ['kinetics'] }, help: 'Effective activation energy of recrystallisation or phase transformation.' },
  // defects
  { key: 'oxide', label: 'Oxide', type: 'select', value: 'ZnO', options: Object.keys(OXIDES).map(o => ({ value: o, label: chem.prettyFormula(o) })), group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Host oxide for the Kröger–Vink defect reactions; sets site stoichiometry and site density.' },
  { key: 'Evf', label: 'Vacancy formation energy (metal)', unit: 'eV', value: 1.0, min: 0.3, max: 4, step: 0.05, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Formation energy of a lattice vacancy in the elemental metal or generic crystal (Cu ≈ 1.2, Al ≈ 0.7, Au ≈ 0.95 eV).' },
  { key: 'Svf', label: 'Vacancy formation entropy', unit: 'k_B', value: 1.5, min: 0, max: 6, step: 0.1, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Vibrational formation entropy in units of the Boltzmann constant.' },
  { key: 'dH_red', label: 'Oxide reduction enthalpy', unit: 'eV', value: 4.0, min: 1, max: 8, step: 0.1, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Enthalpy of O_O → V_O•• + 2e′ + ½O₂. Material specific: calibrate on conductivity or thermogravimetry versus pO₂ data.' },
  { key: 'dS_red', label: 'Oxide reduction entropy', unit: 'k_B', value: 10, min: 0, max: 25, step: 0.5, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Entropy of the reduction reaction (dominated by the released O₂ gas).' },
  { key: 'pO2', label: 'Oxygen partial pressure', unit: 'atm', value: 0.21, min: 1e-30, max: 100, log: true, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Oxygen activity of the annealing atmosphere; reducing atmospheres raise the oxygen-vacancy and electron concentrations.' },
  { key: 'x_donor', label: 'Donor dopant site fraction', value: 0, min: 0, max: 0.05, step: 0.001, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Fully ionised donor concentration per cation site (e.g. Al in ZnO). It fixes the electron density in the extrinsic regime.' },
  { key: 'dH_S', label: 'Schottky formation enthalpy', unit: 'eV', value: 5.0, min: 1, max: 12, step: 0.1, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Enthalpy to form one stoichiometric set of cation and anion vacancies.' },
  { key: 'dH_F', label: 'Cation Frenkel formation enthalpy', unit: 'eV', value: 4.0, min: 1, max: 12, step: 0.1, group: 'Point defects', showIf: { mode: ['defects'] }, help: 'Enthalpy to form one cation interstitial–vacancy pair.' },
  { key: 'rho_d', label: 'Dislocation density', unit: '1/m²', value: 1e14, min: 1e8, max: 1e16, log: true, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Initial total dislocation density (annealed metals ≈ 10¹⁰–10¹², heavily worked ≈ 10¹⁵–10¹⁶ m⁻²).' },
  { key: 'G_shear', label: 'Shear modulus', unit: 'Pa', value: 4.8e10, min: 1e9, max: 5e11, log: true, group: 'Dislocations', bus: 'prop.G', showIf: { mode: ['defects'] }, help: 'Shear modulus used in the Taylor relation and line-energy estimate.' },
  { key: 'b_nm', label: 'Burgers vector', unit: 'nm', value: 0.256, min: 0.1, max: 1, step: 0.001, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Magnitude of the Burgers vector (Cu 0.256 nm, Al 0.286 nm, Fe 0.248 nm).' },
  { key: 'alpha_T', label: 'Taylor constant α', value: 0.3, min: 0.1, max: 1, step: 0.01, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Dislocation interaction strength in σ = σ₀ + MαGb√ρ (typically 0.2–0.5).' },
  { key: 'sigma0', label: 'Friction (lattice + solute) stress', unit: 'MPa', value: 30, min: 0, max: 2000, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Stress contribution not due to forest dislocations.' },
  { key: 'k1', label: 'Storage coefficient k₁', unit: '1/m', value: 2e8, min: 1e6, max: 5e9, log: true, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Kocks–Mecking athermal storage coefficient in dρ/dε = M(k₁√ρ − k₂ρ).' },
  { key: 'k2', label: 'Dynamic-recovery coefficient k₂', value: 8, min: 0.5, max: 100, step: 0.5, group: 'Dislocations', showIf: { mode: ['defects'] }, help: 'Kocks–Mecking dynamic-recovery coefficient; increases with temperature and decreases with strain rate.' },
  { key: 'x_bulk', label: 'Bulk solute fraction', value: 0.001, min: 1e-6, max: 0.2, log: true, group: 'Segregation', showIf: { mode: ['defects'] }, help: 'Atomic fraction of segregating solute in the grain interior.' },
  { key: 'dG_seg', label: 'Segregation free energy', unit: 'kJ/mol', value: 40, min: 0, max: 150, group: 'Segregation', showIf: { mode: ['defects'] }, help: 'Magnitude of the free-energy gain when a solute atom moves from the bulk to a grain-boundary site (McLean).' },
  { key: 'dgrain_nm', label: 'Grain size for boundary inventory', unit: 'nm', value: 1000, min: 2, max: 1e6, log: true, group: 'Segregation', showIf: { mode: ['defects'] }, help: 'Grain diameter used to estimate the fraction of atoms residing in grain boundaries.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
function read(inp) {
  const p = {};
  for (const d of inputs) {
    let v = inp ? inp[d.key] : undefined;
    if (d.type === 'select') p[d.key] = d.options.some(o => o.value === v) ? v : d.value;
    else { v = Number(v); p[d.key] = N.clamp(Number.isFinite(v) ? v : d.value, d.min, d.max); }
  }
  return p;
}
function finish(r) {
  r.warnings = r.warnings || [];
  for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; r.warnings.push(`KPI "${k.label}" could not be evaluated for these inputs and is reported as 0.`); }
  for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
  if (!Object.keys(r.publish).length) r.publish['micro.status'] = 'not evaluated';
  return r;
}

// ---------- 2-D real FFT built from the kernel's 1-D transform ----------
// Only the half spectrum (columns 0..n/2) is stored: hr/hi are (n rows) × (n/2+1 columns). Two real rows are packed into
// one complex transform, which halves the number of kernel calls.
function rfft2(src, hr, hi, n, zr, zi, cr, ci) {
  const m = n / 2 + 1;
  for (let r = 0; r < n; r += 2) {
    for (let c = 0; c < n; c++) { zr[c] = src[r * n + c]; zi[c] = src[(r + 1) * n + c]; }
    N.fft(zr, zi, false);
    for (let c = 0; c < m; c++) { const d = (n - c) % n; hr[r * m + c] = 0.5 * (zr[c] + zr[d]); hi[r * m + c] = 0.5 * (zi[c] - zi[d]); hr[(r + 1) * m + c] = 0.5 * (zi[c] + zi[d]); hi[(r + 1) * m + c] = -0.5 * (zr[c] - zr[d]); }
  }
  for (let c = 0; c < m; c++) {
    for (let r = 0; r < n; r++) { cr[r] = hr[r * m + c]; ci[r] = hi[r * m + c]; }
    N.fft(cr, ci, false);
    for (let r = 0; r < n; r++) { hr[r * m + c] = cr[r]; hi[r * m + c] = ci[r]; }
  }
}
function irfft2(hr, hi, dst, n, zr, zi, cr, ci, wr, wi) {
  const m = n / 2 + 1;
  for (let c = 0; c < m; c++) {
    for (let r = 0; r < n; r++) { cr[r] = hr[r * m + c]; ci[r] = hi[r * m + c]; }
    N.fft(cr, ci, true);
    for (let r = 0; r < n; r++) { wr[r * m + c] = cr[r]; wi[r * m + c] = ci[r]; }
  }
  for (let r = 0; r < n; r += 2) {
    for (let c = 0; c < n; c++) { const k = c < m ? c : n - c, sg = c < m ? 1 : -1, ar = wr[r * m + k], ai = sg * wi[r * m + k], br = wr[(r + 1) * m + k], bi = sg * wi[(r + 1) * m + k]; zr[c] = ar - bi; zi[c] = ai + br; }
    N.fft(zr, zi, true);
    for (let c = 0; c < n; c++) { dst[r * n + c] = zr[c]; dst[(r + 1) * n + c] = zi[c]; }
  }
}

// ---------- Cahn–Hilliard ----------
function spinodal(p) {
  const n = 64, nn = n * n, rnd = N.rng(p.seed), phi0 = 2 * p.c0 - 1, kap = p.kappa, dt = p.ch_dt, steps = Math.round(p.ch_steps), a = 2;
  const dx = p.dx_nm * 1e-9, tau = dx * dx / p.Dint, warnings = [];
  const m = n / 2 + 1, nh = n * m, phi = new Float64Array(nn), g = new Float64Array(nn), hr = new Float64Array(nh), hi = new Float64Array(nh), gr = new Float64Array(nh), gi = new Float64Array(nh), wr = new Float64Array(nh), wi = new Float64Array(nh);
  const zr = new Float64Array(n), zi = new Float64Array(n), cr = new Float64Array(n), ci = new Float64Array(n);
  let m0 = 0; for (let i = 0; i < nn; i++) { phi[i] = p.noise * (2 * rnd() - 1); m0 += phi[i]; } m0 /= nn; for (let i = 0; i < nn; i++) phi[i] += phi0 - m0;
  const k2 = new Float64Array(nh), wgt = new Float64Array(nh), kq = j => 2 * N.PI * (j < n / 2 ? j : j - n) / n;
  for (let r = 0; r < n; r++) for (let c = 0; c < m; c++) { k2[r * m + c] = kq(r) ** 2 + kq(c) ** 2; wgt[r * m + c] = c === 0 || c === n / 2 ? 1 : 2; }
  rfft2(phi, hr, hi, n, zr, zi, cr, ci);
  const nS = 40, every = Math.max(1, Math.floor(steps / nS)), T = [], L = [], FE = [], MASS = [], frames = [], frameEvery = Math.max(1, Math.floor(steps / 10));
  const sample = t => {
    let fb = 0, s1 = 0, s0 = 0, fg = 0, mean = 0;
    for (let i = 0; i < nn; i++) { const q = phi[i] * phi[i] - 1; fb += 0.25 * q * q; mean += phi[i]; }
    for (let i = 1; i < nh; i++) { const S = wgt[i] * (hr[i] * hr[i] + hi[i] * hi[i]); fg += k2[i] * S; s1 += Math.sqrt(k2[i]) * S; s0 += S; }
    T.push(t); FE.push(fb / nn + 0.5 * kap * fg / (nn * nn)); MASS.push(mean / nn); L.push(s1 > 0 ? 2 * N.PI * s0 / s1 : n);
  };
  const snap = () => { const d = new Uint8Array(nn); for (let i = 0; i < nn; i++) d[i] = phi[i] > 0 ? 1 : 0; return d; };
  sample(0); frames.push({ t: 0, data: snap() });
  for (let s = 1; s <= steps; s++) {
    for (let i = 0; i < nn; i++) g[i] = phi[i] * phi[i] * phi[i] - phi[i];
    rfft2(g, gr, gi, n, zr, zi, cr, ci);
    for (let i = 0; i < nh; i++) { const q = k2[i], den = 1 + dt * (a * q + kap * q * q), num = 1 + dt * a * q; hr[i] = (num * hr[i] - dt * q * gr[i]) / den; hi[i] = (num * hi[i] - dt * q * gi[i]) / den; }
    irfft2(hr, hi, phi, n, zr, zi, cr, ci, wr, wi);
    if (s % every === 0 || s === steps) sample(s * dt);
    if (s % frameEvery === 0 && frames.length < 12) frames.push({ t: s * dt * tau, data: snap() });
  }
  let inc = 0; for (let i = 1; i < FE.length; i++) if (FE[i] > FE[i - 1] + 1e-10 * Math.abs(FE[0])) inc++;
  const drift = Math.abs(MASS[MASS.length - 1] - MASS[0]);
  const spin = 1 - 3 * phi0 * phi0, inside = spin > 0;
  const lamM = inside ? 2 * N.PI * Math.sqrt(2 * kap / spin) : NaN, wMax = inside ? spin * spin / (4 * kap) : 0;
  const half = Math.floor(T.length / 2), lt = [], ll = []; for (let i = Math.max(1, half); i < T.length; i++) { lt.push(Math.log(T[i])); ll.push(Math.log(L[i])); }
  const fit = lt.length > 2 ? N.linreg(lt, ll) : { slope: 0, r2: 0 }, nExp = Number.isFinite(fit.slope) ? fit.slope : 0;
  const Lf = L[L.length - 1], fracB = Array.from(phi).filter(v => v > 0).length / nn;
  const amp = (N.maxOf(phi) - N.minOf(phi)) / 2, decomposed = amp > 0.5;
  const wInt = 3.1 * Math.sqrt(kap), sigma = (2 * Math.SQRT2 / 3) * Math.sqrt(kap) * p.f0 * dx;
  if (!inside) warnings.push(`Mean composition ${F(p.c0)} lies outside the spinodal (0.211–0.789 for this double-well free energy): the mixture is metastable and decomposes by nucleation, which this noise-free deterministic model does not capture.`);
  if (Lf > n / 3) warnings.push('The domain size exceeds one third of the periodic box: late-time coarsening statistics are finite-size limited.');
  if (dt > 1) warnings.push('Large time step: the stabilised scheme remains energy-stable but under-predicts the early-stage growth rate.');
  const map = new Float64Array(nn); for (let i = 0; i < nn; i++) map[i] = 0.5 * (phi[i] + 1);
  return {
    summary: decomposed ? `Spinodal decomposition of a ${F(p.c0)} mixture on a 64×64 periodic grid: after ${F(steps * dt * tau)} s (${steps} steps) the bicontinuous/droplet structure has a characteristic domain size of ${F(Lf * p.dx_nm)} nm and coarsens as L ∝ t^${F(nExp, 3)}. The free energy fell monotonically by ${F(100 * (1 - FE[FE.length - 1] / FE[0]), 3)} % and mean composition was conserved to ${F(drift, 2)}.`
      : `The ${F(p.c0)} mixture did not decompose within ${steps} steps (composition amplitude ${F(amp, 3)}): ${inside ? 'the run is still in the linear incubation stage — increase the number of steps or the fluctuation amplitude' : 'the composition is outside the spinodal, so only nucleation and growth could separate it'}. The free-energy and mass checks below still apply.`,
    kpis: [
      kpi('L_domain', 'Characteristic domain size', Lf * dx, 'm', `Inverse first moment of the circularly averaged structure factor, L = 2π/⟨k⟩, at the end of the run (${F(Lf)} grid cells). It measures the mean spacing between like domains.`, Lf <= n / 3),
      kpi('n_coarsen', 'Coarsening exponent', nExp, '–', `Slope of log L against log t over the second half of the run (R² = ${F(fit.r2, 3)}). Diffusion-limited (LSW-type) coarsening gives 1/3; lower values indicate the run is still in the early stage.`, nExp > 0.2 && nExp < 0.45),
      kpi('lambda_max', 'Fastest-growing wavelength (linear theory)', inside ? lamM * dx : 0, 'm', 'Cahn linear stability result λ = 2π·√(2κ/(1−3φ₀²)) for this free energy: the wavelength that dominates the early spinodal pattern. Zero means the composition is outside the spinodal.', inside),
      kpi('t_incub', 'Early-stage amplification time', inside ? tau / wMax : 0, 's', 'Reciprocal of the maximum linear growth rate ω = (1−3φ₀²)²/(4κ) in units of Δx²/D: the e-folding time of the dominant composition wave.', null),
      kpi('frac_B', 'B-rich phase fraction', fracB, '–', `Area fraction with composition above 0.5. The lever rule for fully separated phases predicts ${F(p.c0, 3)}.`, null),
      kpi('sigma_int', 'Interfacial energy', sigma, 'J/m²', 'Equilibrium flat-interface energy σ = (2√2/3)·√κ·f₀·Δx of the quartic double-well model with the chosen κ, f₀ and grid spacing.', sigma > 0.005 && sigma < 3),
      kpi('dF', 'Free-energy reduction', 1 - FE[FE.length - 1] / FE[0], '–', 'Relative decrease of the total (bulk + gradient) free energy between the initial and final states. It must be non-negative for a thermodynamically consistent scheme.', FE[FE.length - 1] <= FE[0]),
      kpi('mass_err', 'Mean-composition drift', drift, '–', 'Absolute change of the spatially averaged order parameter over the run; the conserved Cahn–Hilliard dynamics must keep it at round-off level.', drift < 1e-9),
    ],
    series: [
      { key: 'L_t', title: 'Domain size versus time', xlabel: 'time (s)', ylabel: 'domain size (nm)', logx: true, logy: true, lines: [line('L(t)', T.slice(1).map(t => t * tau), L.slice(1).map(v => v * p.dx_nm))], explain: 'Log–log coarsening curve from the structure-factor first moment. A straight late-time segment of slope ≈ 1/3 indicates diffusion-limited coarsening.' },
      { key: 'F_t', title: 'Free energy versus time', xlabel: 'time (s)', ylabel: 'free energy per cell (f₀)', lines: [line('F(t)', T.map(t => t * tau), FE)], explain: 'Total free energy per grid cell in units of f₀. It must decrease monotonically; the steep drop marks phase separation, the slow tail is interface-area reduction.' },
    ],
    maps: [{ key: 'comp', title: 'Final composition field', nx: n, ny: n, data: map, x0: 0, x1: n * p.dx_nm, y0: 0, y1: n * p.dx_nm, xlabel: 'x (nm)', ylabel: 'y (nm)', unit: 'atomic fraction B', explain: 'Composition after the last step. Interconnected patterns form near 50:50 composition, isolated droplets for off-symmetric mixtures.' }],
    equations: [
      eqn('governing', 'Cahn–Hilliard equation', '∂φ/∂t = ∇·[M ∇(∂f/∂φ − κ∇²φ)],  φ = 2c − 1'),
      eqn('constitutive', 'Double-well free energy', 'F[φ] = ∫ [f₀·¼(φ² − 1)² + ½κ|∇φ|²] dV'),
      eqn('formula', 'Linear stability (Cahn)', 'ω(k) = M k²[(1 − 3φ₀²)f₀ − κk²],  k_max² = (1 − 3φ₀²)f₀/(2κ)'),
      eqn('formula', 'Semi-implicit spectral update', 'φ̂ⁿ⁺¹ = [(1 + aΔt k²)φ̂ⁿ − Δt k² FFT(φ³ − φ)ⁿ] / [1 + Δt(a k² + κk⁴)],  a = 2'),
    ],
    geometry: { kind: 'voxels', nx: n, ny: n, nz: 1, data: snap(), voxel: p.dx_nm, units: 'nm', frames },
    checks: [
      chk('Free energy decays monotonically', inc === 0, `${inc} of ${FE.length - 1} sampled intervals showed an increase`),
      chk('Mass (mean composition) conserved', drift < 1e-9, `drift ${drift.toExponential(2)} in order-parameter units`),
      chk('Interface resolved by the grid', wInt >= 3, `10–90 % interface width ≈ ${F(wInt, 3)} cells (≥ 3 needed)`),
      chk('Composition inside the spinodal', inside, `1 − 3φ₀² = ${F(spin, 3)}`),
      chk('Domain size small compared with the box', Lf <= n / 3, `L/box = ${F(Lf / n, 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', `Calibrate κ and f₀ so that the model interfacial energy (${F(sigma, 3)} J/m²) and interface width (${F(wInt * p.dx_nm, 3)} nm) match measured or CALPHAD/atomistic values before quoting absolute length scales.`, 'Interfacial energy and width are the two quantities the double-well parameters fix.'),
      rec('high', `Compare the predicted dominant wavelength (${F(inside ? lamM * p.dx_nm : 0, 3)} nm) with small-angle scattering peak positions from quenched samples to identify the mobility and the effective diffusivity.`, 'Cahn linear theory relates the scattering peak to κ and the spinodal driving force.'),
      rec('medium', nExp < 0.25 ? 'Run more steps or a larger time step to reach the scaling regime before using the coarsening exponent.' : 'Use the fitted coarsening law to choose the anneal time for a target domain size instead of repeated trial anneals.', `Fitted exponent ${F(nExp, 3)}.`),
      rec('medium', 'Use the lowest annealing temperature and shortest hold that reach the target domain size: coarsening follows t^(1/3), so doubling the size costs eight times the furnace time and energy.', 'Diffusion-limited coarsening scaling.', true),
    ],
    publish: { 'micro.domain_size': Lf * dx, 'prop.surface_energy': sigma, 'micro.phase_fraction': fracB },
    numerics: { method: 'Fourier-spectral Cahn–Hilliard, linearly stabilised semi-implicit time stepping', dof: nn, steps, iterations: steps, residual: drift, notes: `64×64 periodic grid, Δt = ${dt} Δx²/D, constant mobility, 2-D only.` },
  };
}

// ---------- Monte-Carlo Potts grain growth ----------
function grains(p) {
  const n = 64, nn = n * n, rnd = N.rng(p.seed), Q = Math.round(p.Qstates), steps = Math.round(p.mcs), warnings = [];
  const s = new Int16Array(nn); for (let i = 0; i < nn; i++) s[i] = rnd.int(Q);
  const nPart = Math.round(p.fp * nn); for (let k = 0, guard = 0; k < nPart && guard < 20 * nn; guard++) { const i = rnd.int(nn); if (s[i] >= 0) { s[i] = -1; k++; } }
  const nb = new Int32Array(nn * 8);
  for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) { const xm = (x + n - 1) % n, xp = (x + 1) % n, ym = (y + n - 1) % n, yp = (y + 1) % n, o = (y * n + x) * 8; nb[o] = y * n + xm; nb[o + 1] = y * n + xp; nb[o + 2] = ym * n + x; nb[o + 3] = yp * n + x; nb[o + 4] = ym * n + xm; nb[o + 5] = ym * n + xp; nb[o + 6] = yp * n + xm; nb[o + 7] = yp * n + xp; }
  const lab = new Int32Array(nn), stack = new Int32Array(nn);
  const count = () => { // connected components (4-neighbour, periodic)
    lab.fill(-1); const areas = [];
    for (let i = 0; i < nn; i++) { if (lab[i] >= 0 || s[i] < 0) continue; let top = 0, a = 0; stack[top++] = i; lab[i] = areas.length; while (top) { const c = stack[--top]; a++; for (let k = 0; k < 4; k++) { const j = nb[c * 8 + k]; if (lab[j] < 0 && s[j] === s[c]) { lab[j] = areas.length; stack[top++] = j; } } } areas.push(a); }
    return areas;
  };
  const T = [], D = [], NG = [], frames = [], sampleAt = new Set(N.logspace(1, steps, 26).map(Math.round)), frameAt = new Set(N.linspace(0, steps, 9).map(Math.round));
  const snap = () => { const d = new Uint8Array(nn); for (let i = 0; i < nn; i++) d[i] = s[i] < 0 ? 0 : 1 + (s[i] % 254); return d; };
  const rec1 = t => { const a = count(), mA = N.mean(a); T.push(t); NG.push(a.length); D.push(Math.sqrt(4 * mA / N.PI)); return a; };
  frames.push({ t: 0, data: snap() });
  let areas = [];
  for (let t = 1; t <= steps; t++) {
    for (let a = 0; a < nn; a++) {
      const i = rnd.int(nn); if (s[i] < 0) continue; const o = i * 8, cand = s[nb[o + rnd.int(8)]]; if (cand < 0 || cand === s[i]) continue;
      let dE = 0; for (let k = 0; k < 8; k++) { const v = s[nb[o + k]]; dE += (v !== cand ? 1 : 0) - (v !== s[i] ? 1 : 0); }
      if (dE <= 0) s[i] = cand;
    }
    if (sampleAt.has(t) || t === steps) areas = rec1(t);
    if (frameAt.has(t) && frames.length < 12) frames.push({ t, data: snap() });
  }
  const lt = [], ld = []; for (let i = 0; i < T.length; i++) if (T[i] >= 8 && NG[i] >= 8) { lt.push(Math.log(T[i])); ld.push(Math.log(D[i])); }
  const fit = lt.length > 3 ? N.linreg(lt, ld) : { slope: NaN, r2: 0 }, slope = Number.isFinite(fit.slope) ? fit.slope : 0, nGrowth = slope > 0.02 ? 1 / slope : 50;
  const lateT = [], lateD = []; for (let i = Math.floor(T.length * 0.7); i < T.length; i++) { lateT.push(Math.log(T[i])); lateD.push(Math.log(D[i])); }
  const late = lateT.length > 2 ? N.linreg(lateT, lateD).slope : slope, stagn = p.fp > 0 && late < 0.15;
  const dEq = areas.map(a => Math.sqrt(4 * a / N.PI)), dMean = N.mean(dEq), norm = dEq.map(v => v / dMean), h = N.histogram(norm, 12, 0, 3), sdLn = N.std(dEq.filter(v => v > 0).map(Math.log));
  const dxm = p.dx_nm * 1e-9, Df = D[D.length - 1], rPart = 1 / Math.sqrt(N.PI); // single-cell particle, equivalent radius in cells
  const zener3 = p.fp > 0 ? 4 * rPart / (3 * p.fp) : Infinity, zener2 = p.fp > 0 ? 2 * rPart / Math.sqrt(p.fp) : Infinity; // 3-D Zener–Smith diameter; 2-D spacing-controlled estimate
  const TK = p.T_C + 273.15, M = p.M0 * Math.exp(-p.Qgb * 1e3 / (N.R_GAS * TK)), tPhys = ((Df * dxm) ** 2 - (D[0] * dxm) ** 2) / (4 * M * p.gamma_gb);
  if (NG[NG.length - 1] < 20) warnings.push(`Only ${NG[NG.length - 1]} grains remain: grain-size statistics are poor and the periodic box influences growth. Use fewer Monte-Carlo steps.`);
  if (Q < 24) warnings.push('Small Q allows frequent coalescence of like-oriented grains, which accelerates apparent growth.');
  warnings.push('Monte-Carlo time is not physical time; the equivalent anneal time quoted uses the parabolic law with the supplied mobility and boundary energy.');
  const map = new Float64Array(nn); for (let i = 0; i < nn; i++) map[i] = s[i];
  return {
    summary: `Zero-temperature Potts simulation (64×64, Q = ${Q}, ${steps} MCS): the mean grain diameter grew from ${F(D[0] * p.dx_nm, 3)} to ${F(Df * p.dx_nm, 3)} nm with ${NG[NG.length - 1]} grains remaining and a growth exponent n = ${F(nGrowth, 3)} (ideal curvature-driven growth: 2). ${p.fp > 0 ? `With ${F(100 * p.fp, 3)} % pinning particles growth has ${stagn ? 'stagnated' : 'not yet stagnated'}; the 2-D pinning estimate is ${F(zener2 * p.dx_nm, 3)} nm.` : 'No pinning particles are present, so growth continues without limit.'}`,
    kpis: [
      kpi('d_mean', 'Final mean grain diameter', Df * dxm, 'm', `Equivalent-circle diameter of the mean grain area after ${steps} Monte-Carlo steps (${F(Df, 3)} grid cells), from a periodic connected-component analysis.`, null),
      kpi('n_grains', 'Grains remaining', NG[NG.length - 1], '–', 'Number of connected grains in the periodic box at the end of the run. Below about 20 the statistics and the growth law become unreliable.', NG[NG.length - 1] >= 20),
      kpi('n_growth', 'Grain-growth exponent n', nGrowth, '–', `From the log–log slope of mean diameter against Monte-Carlo time (d ∝ t^(1/n), R² = ${F(fit.r2, 3)}). Ideal normal grain growth gives n = 2; pinning and lattice anisotropy raise it.`, nGrowth > 1.6 && nGrowth < 3.2),
      kpi('zener_2d', 'Pinned grain size, 2-D estimate', Number.isFinite(zener2) ? zener2 * dxm : 0, 'm', 'Two-dimensional pinning estimate D ≈ 2r/√f: in a plane every boundary intersects particles, so stagnation occurs near the particle spacing. Zero means no particles.', null),
      kpi('zener_3d', 'Zener–Smith limit (3-D)', Number.isFinite(zener3) ? zener3 * dxm : 0, 'm', 'Classical three-dimensional limiting grain diameter D = 4r/(3f) for randomly distributed spherical particles of the same radius and fraction. Zero means no particles.', null),
      kpi('sd_ln', 'Width of grain-size distribution', sdLn, '–', 'Standard deviation of ln(d). Self-similar normal grain growth gives about 0.3–0.45; larger values indicate abnormal growth or pinning.', sdLn < 0.6),
      kpi('t_equiv', 'Equivalent anneal time', tPhys, 's', `Time for the same diameter change under the parabolic law d² − d₀² = 4Mγt with M = ${F(M, 3)} m⁴/(J·s) at ${F(p.T_C)} °C.`, null),
      kpi('late_slope', 'Late-time growth slope', late, '–', 'd ln d / d ln t over the last 30 % of samples. Values near zero indicate stagnation by particle pinning.', null),
    ],
    series: [
      { key: 'd_t', title: 'Mean grain size versus Monte-Carlo time', xlabel: 'time (MCS)', ylabel: 'mean grain diameter (nm)', logx: true, logy: true, lines: [line('simulation', T, D.map(v => v * p.dx_nm)), line('fit d ∝ t^(1/n)', T, T.map(t => Math.exp((fit.intercept ?? Math.log(D[0])) + slope * Math.log(t)) * p.dx_nm))], explain: 'Log–log growth curve with the fitted power law. A slope of ½ corresponds to parabolic growth; flattening indicates Zener pinning.' },
      { key: 'gsd', title: 'Grain-size distribution', xlabel: 'd / ⟨d⟩', ylabel: 'number of grains', lines: [line('final state', h.x, h.y)], explain: 'Histogram of equivalent grain diameters normalised by the mean. Normal growth produces a single self-similar peak near 1 with a cut-off below about 2.5.' },
    ],
    maps: [{ key: 'grains', title: 'Grain orientation map', nx: n, ny: n, data: map, x0: 0, x1: n * p.dx_nm, y0: 0, y1: n * p.dx_nm, xlabel: 'x (nm)', ylabel: 'y (nm)', unit: 'orientation id (−1 = particle)', explain: 'Final orientation index per cell; each uniformly coloured region is one grain and isolated cells of value −1 are pinning particles.' }],
    equations: [
      eqn('constitutive', 'Potts boundary energy', 'E = ½ J Σᵢ Σ_{j∈nn(i)} (1 − δ(sᵢ, sⱼ)),  8 neighbours'),
      eqn('governing', 'Monte-Carlo acceptance (T = 0)', 'sᵢ → s_neighbour accepted if ΔE ≤ 0'),
      eqn('governing', 'Curvature-driven boundary migration', 'v = M γ κ  ⇒  d² − d₀² = 4 M γ t'),
      eqn('formula', 'Zener–Smith limit', 'D_Z = 4 r / (3 f)   (3-D);   D ≈ 2 r / √f   (2-D section)'),
    ],
    geometry: { kind: 'voxels', nx: n, ny: n, nz: 1, data: snap(), voxel: p.dx_nm, units: 'nm', frames },
    checks: [
      chk('Growth exponent near parabolic', nGrowth > 1.6 && nGrowth < 3.2 || stagn, `n = ${F(nGrowth, 3)}${stagn ? ' (pinned)' : ''}`),
      chk('Enough grains for statistics', NG[NG.length - 1] >= 20, `${NG[NG.length - 1]} grains`),
      chk('Mean grain size never decreases', D.every((v, i) => i === 0 || v >= D[i - 1] * 0.98), 'monotonic within 2 % sampling noise'),
      chk('Particles conserved', Array.from(s).filter(v => v < 0).length === nPart, `${nPart} inert sites`),
    ],
    warnings,
    recommendations: [
      rec('high', p.fp > 0 ? `To hold the grain size below ${F(zener3 * p.dx_nm, 3)} nm in a bulk part, keep particle radius and fraction on the Zener line r/f = ${F(0.75 * zener3 * p.dx_nm, 3)} nm; finer particles at the same fraction pin more strongly.` : 'Add a fine, stable second phase (about 1–2 vol %) if the nanocrystalline grain size must survive annealing: without pinning, growth is unbounded.', 'Zener–Smith pinning.'),
      rec('medium', 'Calibrate the Monte-Carlo time scale by matching one measured grain size after a known anneal, then reuse the fitted mobility at other temperatures through its Arrhenius law.', 'Monte-Carlo steps carry no intrinsic time unit.'),
      rec('medium', 'Validate the simulated grain-size distribution against EBSD or TEM grain maps (at least 300 grains) rather than against the mean size only.', 'Distribution width is a sensitive indicator of abnormal growth.'),
      rec('low', 'Lower the annealing temperature rather than shortening the hold when grain growth must be limited: the mobility falls exponentially, which also reduces furnace energy use.', `Arrhenius mobility with Q = ${F(p.Qgb)} kJ/mol.`, true),
    ],
    publish: { 'micro.grain_size': Df * dxm, 'micro.growth_exponent': nGrowth },
    numerics: { method: 'Monte-Carlo Potts model (zero temperature, neighbour-orientation proposals, 8-neighbour energy)', dof: nn, steps, iterations: steps * nn, residual: 0, notes: 'Periodic 64×64 lattice; grain statistics by 4-connected component labelling.' },
  };
}

// ---------- Kobayashi phase-field solidification ----------
function solid(p) {
  const n = 96, nn = n * n, h = 0.03, dt = 2e-4, tauK = 3e-4, eb = 0.01, al = 0.9, gam = 10, jf = +p.fold, del = p.aniso, Kl = p.latent, steps = Math.round(p.sol_steps), rnd = N.rng(p.seed), warnings = [];
  const ph = new Float64Array(nn), a1 = new Float64Array(nn), a2 = new Float64Array(nn), e2 = new Float64Array(nn).fill(eb * eb), dph = new Float64Array(nn), cT = dt / (h * h); let T = new Float64Array(nn), lapT = new Float64Array(nn);
  const c = n / 2; for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) if ((x - c + 0.5) ** 2 + (y - c + 0.5) ** 2 < 16) ph[y * n + x] = 1;
  const heat0 = N.sum(T) - Kl * N.sum(ph);
  const tip = () => { const r = Math.floor(c) * n; for (let x = n - 2; x >= c; x--) if (ph[r + x] >= 0.5) { const d = ph[r + x] - ph[r + x + 1]; return (x + (d > 1e-12 ? (ph[r + x] - 0.5) / d : 0) - c + 0.5); } return 0; };
  const tt = [], tipX = [], fs = [], frames = [], every = Math.max(1, Math.floor(steps / 40)), fevery = Math.max(1, Math.floor(steps / 10));
  const snap = () => { const d = new Uint8Array(nn); for (let i = 0; i < nn; i++) d[i] = ph[i] > 0.5 ? 1 : 0; return d; };
  tt.push(0); tipX.push(tip()); fs.push(N.sum(ph) / nn); frames.push({ t: 0, data: snap() });
  let x0 = Math.floor(c) - 8, x1 = Math.floor(c) + 8, y0 = x0, y1 = x1; // active window around the diffuse interface
  for (let s = 1; s <= steps; s++) {
    x0 = Math.max(0, x0); y0 = Math.max(0, y0); x1 = Math.min(n - 1, x1); y1 = Math.min(n - 1, y1);
    for (let y = y0; y <= y1; y++) { const ym = y > 0 ? y - 1 : 0, yp = y < n - 1 ? y + 1 : n - 1; for (let x = x0; x <= x1; x++) { const xm = x > 0 ? x - 1 : 0, xp = x < n - 1 ? x + 1 : n - 1, i = y * n + x, px = (ph[y * n + xp] - ph[y * n + xm]) / (2 * h), py = (ph[yp * n + x] - ph[ym * n + x]) / (2 * h); const g2 = px * px + py * py; if (g2 < 1e-12) { a1[i] = 0; a2[i] = 0; e2[i] = eb * eb; continue; } const inv = 1 / Math.sqrt(g2), cr = px * inv, ci = py * inv; let pr = 1, pi = 0; for (let k = 0; k < jf; k++) { const t = pr * cr - pi * ci; pi = pr * ci + pi * cr; pr = t; } const e = eb * (1 + del * pr), ed = -eb * jf * del * pi; a1[i] = e * ed * py; a2[i] = e * ed * px; e2[i] = e * e; } }
    let nx0 = n, nx1 = -1, ny0 = n, ny1 = -1;
    for (let y = y0; y <= y1; y++) { const ym = y > 0 ? y - 1 : 0, yp = y < n - 1 ? y + 1 : n - 1; for (let x = x0; x <= x1; x++) { const xm = x > 0 ? x - 1 : 0, xp = x < n - 1 ? x + 1 : n - 1, i = y * n + x, il = y * n + xm, ir = y * n + xp, id = ym * n + x, iu = yp * n + x, q = ph[i]; if ((q === 0 || q === 1) && ph[il] === q && ph[ir] === q && ph[id] === q && ph[iu] === q && a1[il] === 0 && a1[ir] === 0 && a2[id] === 0 && a2[iu] === 0) { dph[i] = 0; continue; } const t1 = -(a1[ir] - a1[il]) / (2 * h), t2 = (a2[iu] - a2[id]) / (2 * h); const t3 = ((e2[ir] + e2[i]) * (ph[ir] - q) - (e2[i] + e2[il]) * (q - ph[il]) + (e2[iu] + e2[i]) * (ph[iu] - q) - (e2[i] + e2[id]) * (q - ph[id])) / (2 * h * h); const w = q * (1 - q), m = w > 1e-12 || w < -1e-12 ? al / N.PI * Math.atan(gam * (1 - T[i])) : 0; const d = dt / tauK * (t1 + t2 + t3 + w * (q - 0.5 + m + (w > 1e-4 ? 0.02 * (rnd() - 0.5) : 0))); dph[i] = d; if (d > 1e-9 || d < -1e-9) { if (x < nx0) nx0 = x; if (x > nx1) nx1 = x; if (y < ny0) ny0 = y; if (y > ny1) ny1 = y; } } }
    for (let y = 1; y < n - 1; y++) { const r = y * n; for (let i = r + 1; i < r + n - 1; i++) { const v = T[i] + cT * (T[i - 1] + T[i + 1] + T[i - n] + T[i + n] - 4 * T[i]); lapT[i] = v < 1e-30 && v > -1e-30 ? 0 : v; } }
    for (let q = 0; q < n; q++) for (const [x, y] of [[q, 0], [q, n - 1], [0, q], [n - 1, q]]) { const xm = x > 0 ? x - 1 : 0, xp = x < n - 1 ? x + 1 : n - 1, ym = y > 0 ? y - 1 : 0, yp = y < n - 1 ? y + 1 : n - 1, i = y * n + x; lapT[i] = T[i] + cT * (T[y * n + xm] + T[y * n + xp] + T[ym * n + x] + T[yp * n + x] - 4 * T[i]); }
    { const sw = T; T = lapT; lapT = sw; }
    for (let y = y0; y <= y1; y++) for (let x = x0; x <= x1; x++) { const i = y * n + x, d = dph[i]; if (d === 0) continue; let q = ph[i] + d; if (q < 1e-30 && q > -1e-30) q = 0; else if (q > 1 - 1e-13 && q < 1 + 1e-13) q = 1; T[i] += Kl * (q - ph[i]); ph[i] = q; dph[i] = 0; }
    if (nx1 >= 0) { x0 = Math.min(x0, nx0 - 3); x1 = Math.max(x1, nx1 + 3); y0 = Math.min(y0, ny0 - 3); y1 = Math.max(y1, ny1 + 3); }
    if (s % every === 0 || s === steps) { tt.push(s * dt); tipX.push(tip()); fs.push(N.sum(ph) / nn); }
    if (s % fevery === 0 && frames.length < 12) frames.push({ t: s * dt, data: snap() });
  }
  const heat1 = N.sum(T) - Kl * N.sum(ph), heatErr = Math.abs(heat1 - heat0) / Math.max(1e-12, Kl * N.sum(ph));
  const hf = Math.floor(tt.length / 2), reg = N.linreg(tt.slice(hf), tipX.slice(hf).map(v => v * h)), vTip = Number.isFinite(reg.slope) ? reg.slope : 0;
  let area = 0, per = 0; for (let y = 1; y < n - 1; y++) for (let x = 1; x < n - 1; x++) { const i = y * n + x; area += ph[i]; per += Math.hypot(ph[i + 1] - ph[i - 1], ph[i + n] - ph[i - n]) / 2; }
  const compact = area > 0 ? per * per / (4 * N.PI * area) : 1;
  const rAlong = tipX[tipX.length - 1], ang = N.PI / jf; let rOff = 0; for (let r = 1; r < n / 2 - 1; r += 0.25) { const x = Math.round(c - 0.5 + r * Math.cos(ang)), y = Math.round(c - 0.5 + r * Math.sin(ang)); if (ph[y * n + x] >= 0.5) rOff = r; }
  const armRatio = rOff > 0 ? rAlong / rOff : 1;
  const dxm = p.dx_um * 1e-6, L0 = dxm / h, t0 = L0 * L0 / p.alpha_th, vPhys = vTip * L0 / t0, hitWall = rAlong > n / 2 - 6;
  const Tmin = N.minOf(T), Tmax = N.maxOf(T), bounded = N.minOf(ph) > -0.05 && N.maxOf(ph) < 1.05;
  if (hitWall) warnings.push('The dendrite tip has reached the domain boundary: tip velocity and solid fraction are affected by the zero-flux walls. Reduce the number of steps.');
  if (del * (jf * jf - 1) > 1) warnings.push('Anisotropy exceeds the interface-stiffness limit δ(j²−1) < 1: missing orientations (facets) would appear and the smooth-interface model is ill-posed there.');
  warnings.push('Kobayashi model parameters are dimensionless and qualitative; the physical tip velocity follows only from the chosen length and diffusivity scaling.');
  const mp = new Float64Array(nn), mT = new Float64Array(nn); mp.set(ph); mT.set(T);
  const ext = n * p.dx_um;
  return {
    summary: `Kobayashi phase-field growth of a ${jf}-fold crystal into an undercooled melt (dimensionless undercooling ${F(1 / Kl, 3)}): after ${steps} steps the solid fraction is ${F(fs[fs.length - 1], 3)} and the primary tip advances at ${F(vTip, 3)} (dimensionless), i.e. ${F(vPhys, 3)} m/s with the chosen scaling. The shape factor P²/4πA = ${F(compact, 3)} indicates a ${compact > 2.5 ? 'dendritic, branched' : compact > 1.4 ? 'faceted/star-shaped' : 'compact'} morphology.`,
    kpis: [
      kpi('v_tip', 'Tip velocity', vPhys, 'm/s', `Slope of primary-arm tip position against time over the second half of the run (dimensionless ${F(vTip, 3)}), scaled with length L₀ = Δx/0.03 and time L₀²/α.`, null),
      kpi('r_tip', 'Primary arm length', rAlong * dxm, 'm', 'Distance from the seed centre to the φ = 0.5 contour along the preferred growth direction at the end of the run.', !hitWall),
      kpi('f_solid', 'Solid fraction', fs[fs.length - 1], '–', 'Area fraction with order parameter φ > 0 (integral of φ) at the end of the run; the adiabatic limit for complete recalescence is the dimensionless undercooling 1/K.', null),
      kpi('shape', 'Shape factor P²/(4πA)', compact, '–', 'Perimeter squared over 4π times area of the solid: 1 for a circle, larger for branched dendrites. Perimeter is evaluated as the integral of |∇φ|.', null),
      kpi('arm_ratio', 'Arm anisotropy ratio', armRatio, '–', `Ratio of crystal radius along a preferred direction to that half-way between two arms (at ${F(180 / jf, 3)}°). Values well above 1 show developed primary arms.`, armRatio > 1.1),
      kpi('T_int', 'Maximum melt temperature', Tmax, '–', 'Highest dimensionless temperature (0 = initial undercooled melt, 1 = equilibrium melting point) produced by latent-heat release; values near 1 show recalescence at the interface. Because the model has no re-melting constraint the solid interior can exceed 1 when K > 1.', Tmax <= Kl + 0.05),
      kpi('heat_err', 'Enthalpy balance error', heatErr, '–', 'Relative drift of the conserved enthalpy ∫(T − Kφ)dA with zero-flux boundaries; it measures the discretisation consistency of the coupled update.', heatErr < 1e-6),
    ],
    series: [
      { key: 'tip_t', title: 'Tip position versus time', xlabel: 'time (s)', ylabel: 'tip position (µm)', lines: [line('primary arm', tt.map(t => t * t0), tipX.map(v => v * p.dx_um))], explain: 'Position of the primary dendrite tip. After the initial transient the curve becomes linear: steady-state tip velocity selected by undercooling and anisotropy.' },
      { key: 'fs_t', title: 'Solid fraction versus time', xlabel: 'time (s)', ylabel: 'solid fraction', lines: [line('solid fraction', tt.map(t => t * t0), fs)], explain: 'Growth of the transformed fraction. Its curvature reflects the increasing interface length of the branching crystal.' },
    ],
    maps: [
      { key: 'phi', title: 'Crystal morphology (order parameter)', nx: n, ny: n, data: mp, x0: 0, x1: ext, y0: 0, y1: ext, xlabel: 'x (µm)', ylabel: 'y (µm)', unit: 'φ (1 = solid)', explain: 'Phase-field order parameter at the end of the run: the diffuse interface separates crystal (1) from melt or solution (0).' },
      { key: 'temp', title: 'Temperature (or supersaturation) field', nx: n, ny: n, data: mT, x0: 0, x1: ext, y0: 0, y1: ext, xlabel: 'x (µm)', ylabel: 'y (µm)', unit: 'dimensionless T', explain: 'Latent heat rejected at the interface warms the surroundings; steep gradients ahead of the tips drive fast growth, while warm pockets between arms slow it.' },
    ],
    equations: [
      eqn('governing', 'Anisotropic Allen–Cahn (Kobayashi)', 'τ ∂φ/∂t = −∂ₓ(ε ε′ ∂ᵧφ) + ∂ᵧ(ε ε′ ∂ₓφ) + ∇·(ε²∇φ) + φ(1−φ)(φ − ½ + m)'),
      eqn('constitutive', 'Anisotropy and driving force', `ε(θ) = ε̄[1 + δ cos(${jf}θ)], θ = angle of ∇φ,  m(T) = (α/π)·atan[γ(T_eq − T)]`),
      eqn('governing', 'Heat (or solute) balance', '∂T/∂t = ∇²T + K ∂φ/∂t'),
      eqn('formula', 'Interface-stiffness admissibility', 'ε + ε″ > 0  ⇔  δ(j² − 1) < 1'),
    ],
    geometry: { kind: 'voxels', nx: n, ny: n, nz: 1, data: snap(), voxel: p.dx_um, units: 'µm', frames },
    checks: [
      chk('Enthalpy conserved', heatErr < 1e-6, `relative drift ${heatErr.toExponential(2)}`),
      chk('Order parameter bounded', bounded, `φ ∈ [${F(N.minOf(ph), 3)}, ${F(N.maxOf(ph), 3)}]`),
      chk('Tip clear of the boundary', !hitWall, `tip at ${F(rAlong, 3)} of ${n / 2} cells`),
      chk('Interface stiffness positive', del * (jf * jf - 1) < 1, `δ(j²−1) = ${F(del * (jf * jf - 1), 3)}`),
      chk('Temperature within the adiabatic range', Tmin > -0.05 && Tmax < Kl + 0.05, `T ∈ [${F(Tmin, 3)}, ${F(Tmax, 3)}], adiabatic maximum K = ${F(Kl)}`),
    ],
    warnings,
    recommendations: [
      rec('high', 'Treat the morphology trends (arm number, branching against undercooling and anisotropy) as the result; calibrate the kinetic coefficient and capillary length against measured tip velocities before using absolute speeds.', 'The Kobayashi model is qualitative without thin-interface calibration.'),
      rec('medium', compact > 2.5 ? 'To obtain compact, facet-free crystals reduce the driving force (lower undercooling or supersaturation): the present conditions are in the dendritic regime.' : 'Raise the driving force if a high-surface-area dendritic morphology is wanted; the present conditions give a compact crystal.', `Shape factor ${F(compact, 3)}.`),
      rec('medium', 'Repeat with a refined grid (halved spacing) to confirm that arm spacing and tip speed are not set by grid anisotropy.', 'Interface-width and grid convergence requirement of phase-field models.'),
      rec('low', 'Operating at low supersaturation gives slower, more uniform growth with less secondary nucleation and less precursor lost to fines, which reduces chemical waste per usable crystal.', 'Morphological stability at low driving force.', true),
    ],
    publish: { 'growth.rate': Math.abs(vPhys), 'micro.solid_fraction': fs[fs.length - 1] },
    numerics: { method: 'Explicit finite differences (Kobayashi anisotropic phase field + heat equation)', dof: 2 * nn, steps, iterations: steps, residual: heatErr, notes: `96×96 grid, Δx = 0.03, Δt = 2×10⁻⁴, zero-flux boundaries, seeded interface noise.` },
  };
}

// ---------- analytical kinetics ----------
function kinetics(p) {
  const TK = p.T_C + 273.15, RT = N.R_GAS * TK, t = p.time_h * 3600, warnings = [];
  const M = p.M0 * Math.exp(-p.Qgb * 1e3 / RT), d0 = p.d0_nm * 1e-9, Kgg = 4 * M * p.gamma_gb;
  const dFree = Math.sqrt(d0 * d0 + Kgg * t), dZ = p.fp > 0 ? 4 * p.rp_nm * 1e-9 / (3 * p.fp) : Infinity;
  const ts = N.linspace(0, t, 120);
  // pinned growth: d(d²)/dt = 4Mγ(1 − d/d_Z), integrated in y = d²
  const pin = N.rk4((_, y) => [Number.isFinite(dZ) ? Kgg * Math.max(0, 1 - Math.sqrt(Math.max(y[0], 0)) / dZ) : Kgg], [d0 * d0], 0, t, 600);
  const dPin = N.linspace(0, t, 120).map(tv => Math.sqrt(Math.max(0, N.interp1(pin.t, pin.y.map(v => v[0]), tv)))), dEnd = dPin[dPin.length - 1];
  const D = p.D0 * Math.exp(-p.Qd * 1e3 / RT), Vm = p.Vm * 1e-6, Klsw = 8 * p.gamma_int * D * p.c_sol * Vm * Vm / (9 * RT), r0 = p.r0_nm * 1e-9, rEnd = Math.cbrt(r0 ** 3 + Klsw * t);
  const Om = Vm / N.N_A, kT = N.K_B * TK, Dp = p.dp_nm * 1e-9, gs = p.gamma_gb;
  const S = { vol: { n: 5, m: 3, B: 80 * D * gs * Om / kT, dens: true, name: 'volume diffusion' }, gb: { n: 6, m: 4, B: 20 * 5e-10 * D * gs * Om / kT, dens: true, name: 'grain-boundary diffusion (δ = 0.5 nm)' }, surf: { n: 7, m: 4, B: 56 * D * gs * Math.pow(Om, 4 / 3) / kT, dens: false, name: 'surface diffusion' } }[p.sinter];
  const neck = tv => Math.min(0.5, Math.pow(S.B * tv / Dp ** S.m, 1 / S.n)), xr = neck(t), tInit = Math.pow(0.3, S.n) * Dp ** S.m / S.B;
  const shrink = x => S.dens ? Math.min(0.12, (x / 2) ** 2) : 0, rho = x => Math.min(1, p.rho_green / (1 - shrink(x)) ** 3), rhoEnd = rho(xr);
  const kA = p.k0_av * Math.exp(-p.Q_av * 1e3 / RT), X = tv => 1 - Math.exp(-Math.pow(kA * tv, p.n_av)), Xend = X(t), t50 = Math.pow(Math.LN2, 1 / p.n_av) / kA;
  if (xr >= 0.3) warnings.push(`Neck ratio X/D reached ${F(xr, 3)}: the initial-stage sintering model is only valid below about 0.3 (reached after ${F(tInit, 3)} s); later densification needs intermediate/final-stage models.`);
  if (d0 >= dZ) warnings.push('Initial grain size already exceeds the Zener limit: the particles prevent any normal grain growth.');
  if (p.c_sol * Vm > 0.1) warnings.push('Solubility is not dilute: LSW theory assumes a dilute matrix and vanishing precipitate volume fraction.');
  const tl = N.logspace(Math.max(t * 1e-4, 1e-3), t, 100);
  return {
    summary: `At ${F(p.T_C)} °C for ${F(p.time_h)} h: grains grow from ${F(p.d0_nm)} to ${F(dEnd * 1e9, 3)} nm (${F(dFree * 1e9, 3)} nm without pinning; Zener limit ${Number.isFinite(dZ) ? F(dZ * 1e9, 3) + ' nm' : 'none'}), precipitates coarsen from ${F(p.r0_nm)} to ${F(rEnd * 1e9, 3)} nm radius (LSW), the sinter neck ratio reaches ${F(xr, 3)} by ${S.name} with relative density ${F(rhoEnd, 3)}, and the JMAK transformed fraction is ${F(Xend, 3)}.`,
    kpis: [
      kpi('d_final', 'Grain size after hold', dEnd, 'm', `Parabolic normal grain growth with Zener drag, d(d²)/dt = 4Mγ(1 − d/d_Z), using M = ${F(M, 3)} m⁴/(J·s). Without particles the result would be ${F(dFree * 1e9, 3)} nm.`, dEnd < 3 * d0),
      kpi('d_zener', 'Zener limiting grain size', Number.isFinite(dZ) ? dZ : 0, 'm', 'Zener–Smith limit D = 4r/(3f) for the given particle radius and fraction; grain growth stagnates when the mean size approaches it. Zero means no particles.', null),
      kpi('r_lsw', 'Precipitate radius after hold', rEnd, 'm', `LSW diffusion-limited coarsening r³ − r₀³ = Kt with K = 8γDc∞Vₘ²/(9RT) = ${F(Klsw, 3)} m³/s and D = ${F(D, 3)} m²/s.`, rEnd < 2 * r0),
      kpi('neck', 'Sinter neck ratio X/D', xr, '–', `Initial-stage neck growth (X/D)^${S.n} = B·t/D^${S.m} for ${S.name}. The model is valid up to about 0.3 and is capped at 0.5.`, xr < 0.3),
      kpi('rho_rel', 'Relative density after hold', rhoEnd, '–', S.dens ? 'From the initial-stage shrinkage ΔL/L₀ = (X/2D)² applied to the green density, ρ = ρ₀/(1 − ΔL/L₀)³.' : 'Surface diffusion grows necks without bringing particle centres together, so the density stays at the green value.', null),
      kpi('t_initial', 'Time to end of initial sintering stage', tInit, 's', 'Hold time at which the neck ratio reaches 0.3; beyond this pores form an interconnected network and intermediate-stage models apply.', null),
      kpi('X_jmak', 'Transformed fraction (JMAK)', Xend, '–', `Johnson–Mehl–Avrami–Kolmogorov fraction X = 1 − exp(−(kt)ⁿ) with k = ${F(kA, 3)} 1/s and n = ${F(p.n_av)}.`, null),
      kpi('t50', 'Half-transformation time', t50, 's', 'Time at which half of the material has recrystallised or transformed, t₅₀ = (ln 2)^(1/n)/k.', null),
    ],
    series: [
      { key: 'gg', title: 'Grain size versus time', xlabel: 'time (h)', ylabel: 'grain diameter (nm)', lines: [line('with Zener pinning', ts.map(v => v / 3600), dPin.map(v => v * 1e9)), line('unpinned parabolic law', ts.map(v => v / 3600), ts.map(v => Math.sqrt(d0 * d0 + Kgg * v) * 1e9))], explain: 'Mean grain diameter during the isothermal hold. The pinned curve saturates at the Zener limit while the unpinned parabolic law grows without bound.' },
      { key: 'lsw', title: 'Precipitate coarsening (LSW)', xlabel: 'time (h)', ylabel: 'mean radius (nm)', lines: [line('r(t)', ts.map(v => v / 3600), ts.map(v => Math.cbrt(r0 ** 3 + Klsw * v) * 1e9))], explain: 'Mean precipitate radius from the cubic coarsening law. Number density falls as 1/r³ at constant volume fraction, which reduces Orowan strengthening.' },
      { key: 'sint', title: 'Sintering: neck ratio and density', xlabel: 'time (s)', ylabel: 'X/D and relative density', logx: true, lines: [line('neck ratio X/D', tl, tl.map(neck)), line('relative density', tl, tl.map(v => rho(neck(v))))], explain: 'Initial-stage neck growth and the resulting densification. The curves are capped where the model loses validity (X/D ≈ 0.3–0.5).' },
      { key: 'jmak', title: 'Transformation kinetics (JMAK)', xlabel: 'time (s)', ylabel: 'transformed fraction', logx: true, lines: [line('X(t)', tl, tl.map(X))], explain: 'Sigmoidal transformed fraction on a logarithmic time axis; the Avrami exponent sets the steepness and the rate constant the position.' },
    ],
    tables: [{ key: 'rates', title: 'Rate constants at temperature', columns: ['Process', 'Law', 'Rate constant', 'Unit'], rows: [['Grain growth', 'd² − d₀² = K t', F(Kgg, 4), 'm²/s'], ['LSW coarsening', 'r³ − r₀³ = K t', F(Klsw, 4), 'm³/s'], ['Sintering (' + S.name + ')', `(X/D)^${S.n} = B t / D^${S.m}`, F(S.B, 4), `m^${S.m}/s`], ['JMAK', 'X = 1 − exp(−(kt)ⁿ)', F(kA, 4), '1/s'], ['Diffusivity', 'D = D₀ exp(−Q/RT)', F(D, 4), 'm²/s']], explain: 'Arrhenius-evaluated rate constants of each kinetic law at the hold temperature, for use in hand calculations or other suites.' }],
    equations: [
      eqn('governing', 'Normal grain growth with Zener drag', 'dd/dt = 2Mγ(1/d − 1/d_Z),  M = M₀ exp(−Q/RT),  d_Z = 4r/(3f)'),
      eqn('governing', 'LSW coarsening', 'r³ − r₀³ = [8 γ D c∞ Vₘ² / (9RT)] t'),
      eqn('formula', 'Initial-stage sintering', `(X/D)^${S.n} = B t / D^${S.m};  ΔL/L₀ = (X/2D)²;  ρ = ρ₀/(1 − ΔL/L₀)³`),
      eqn('formula', 'JMAK transformation', 'X(t) = 1 − exp[−(k t)ⁿ],  k = k₀ exp(−Q/RT)'),
    ],
    checks: [
      chk('Grain size does not exceed Zener limit', !Number.isFinite(dZ) || dEnd <= Math.max(d0, dZ) * 1.001, Number.isFinite(dZ) ? `d = ${F(dEnd * 1e9, 3)} nm, d_Z = ${F(dZ * 1e9, 3)} nm` : 'no particles'),
      chk('Initial-stage sintering model valid', xr < 0.3, `X/D = ${F(xr, 3)}`),
      chk('Relative density admissible', rhoEnd <= 1 && rhoEnd >= p.rho_green, `ρ = ${F(rhoEnd, 3)}`),
      chk('Transformed fraction within [0, 1]', Xend >= 0 && Xend <= 1, `X = ${F(Xend, 4)}`),
      chk('Dilute-solution assumption (LSW)', p.c_sol * Vm <= 0.1, `c∞Vₘ = ${F(p.c_sol * Vm, 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', dEnd > 2 * d0 ? `Grain size more than doubles during this hold. To retain the ${F(p.d0_nm)} nm structure lower the temperature, shorten the hold, or add pinning particles with r/f below ${F(0.75 * d0 * 1e9, 3)} nm.` : 'The nanostructure is retained under these conditions; confirm with XRD line broadening or TEM after the anneal.', `d/d₀ = ${F(dEnd / d0, 3)}.`),
      rec('high', 'Calibrate M₀ and Q on at least three isothermal grain-size series; mobility of nanocrystalline oxides varies by orders of magnitude with impurity content.', 'Arrhenius parameters dominate every result of this analysis.'),
      rec('medium', xr >= 0.3 ? 'Neck growth leaves the initial stage within the hold: use a two-step or field-assisted sintering schedule to densify without the accompanying grain growth.' : `Densification is still in the initial stage; ${F(tInit, 3)} s would be needed to complete it at this temperature.`, 'Initial-stage neck model.'),
      rec('medium', `Herring scaling: halving the particle size shortens the time for the same neck ratio by a factor ${F(2 ** S.m, 3)}, allowing a lower sintering temperature and a substantial saving in furnace energy.`, `Exponent m = ${S.m} of the selected mechanism.`, true),
    ],
    publish: { 'micro.grain_size': dEnd, 'geom.porosity': 1 - rhoEnd, 'geom.rel_density': rhoEnd, 'prop.D_diff': D, 'rxn.conversion': Xend, 'micro.precipitate_radius': rEnd },
    numerics: { method: 'Closed-form kinetic laws; RK4 integration of the pinned grain-growth ODE', dof: 1, steps: 600, iterations: 600, residual: 0, notes: 'Isothermal hold; no heating/cooling ramps.' },
  };
}

// ---------- defect equilibria ----------
function defects(p) {
  const TK = p.T_C + 273.15, kT = KB_EV * TK, warnings = [];
  const cv = Math.exp(p.Svf) * Math.exp(-p.Evf / kT);
  const ox = p.oxide, pf = chem.parseFormula(ox), metal = Object.keys(pf.atoms).find(e => e !== 'O'), a = pf.atoms[metal], b = pf.atoms.O, z = Math.round(2 * b / a);
  const Mox = chem.molarMass(ox) * 1e-3, nFU = OXIDES[ox] * N.N_A / Mox, nO = nFU * b, nCat = nFU * a;
  const Kred = Math.exp(p.dS_red) * Math.exp(-p.dH_red / kT);
  // electroneutrality n = 2[V_O] + [D], mass action [V_O] n² pO2^½ = K   (site fractions per formula unit)
  const dn = p.x_donor * a;
  const solve = pO2 => { const rhs = Kred / Math.sqrt(pO2); const f = lv => { const v = Math.exp(lv); return Math.log(v) + 2 * Math.log(2 * v + dn) - Math.log(rhs); }; let lv; try { lv = N.brent(f, -400, 80); } catch { lv = f(0) > 0 ? -400 : 80; } const v = Math.exp(lv); return { v, n: 2 * v + dn }; };
  const here = solve(p.pO2), vO = here.v / b, nE = here.n * nFU;
  const KS = Math.exp(-p.dH_S / kT), xS = Math.pow(KS / (a ** a * b ** b), 1 / (a + b)), xF = Math.exp(-p.dH_F / (2 * kT));
  const lp = N.linspace(-30, 2, 81), bro = lp.map(l => solve(10 ** l));
  const i0 = N.clamp(Math.round((Math.log10(p.pO2) + 30) / 0.4), 1, 79), slope = (Math.log10(bro[i0 + 1].v) - Math.log10(bro[i0 - 1].v)) / 0.8;
  if (here.v / b > 0.02) warnings.push(`Oxygen-vacancy site fraction ${F(here.v / b, 3)} exceeds the dilute limit (~1 %): defect interactions and phase decomposition are neglected by ideal mass action${here.v / b >= 1 ? '; a fraction of 1 or more is unphysical and means the oxide is not stable under these conditions' : ''}.`);
  // dislocations
  const bb = p.b_nm * 1e-9, Mt = 3.06, tay = r => p.sigma0 * 1e6 + Mt * p.alpha_T * p.G_shear * bb * Math.sqrt(r), sig = tay(p.rho_d), rhoSat = (p.k1 / p.k2) ** 2;
  const eps = N.linspace(0, 0.3, 80), rhoE = eps.map(e => (p.k1 / p.k2 + (Math.sqrt(p.rho_d) - p.k1 / p.k2) * Math.exp(-Mt * p.k2 * e / 2)) ** 2), sigE = rhoE.map(tay);
  const Est = 0.5 * p.G_shear * bb * bb * p.rho_d, spacing = 1 / Math.sqrt(p.rho_d);
  // McLean segregation
  const mcl = T => { const q = p.x_bulk / (1 - p.x_bulk) * Math.exp(p.dG_seg * 1e3 / (N.R_GAS * T)); return q / (1 + q); };
  const xgb = mcl(TK), beta = xgb / p.x_bulk, fgb = Math.min(1, 3 * 0.5 / p.dgrain_nm), need = fgb * xgb, depl = need > p.x_bulk;
  if (depl) warnings.push('The grain-boundary area is large enough to deplete the bulk solute: the McLean isotherm with fixed bulk composition over-estimates the coverage for this grain size.');
  const Ts = N.linspace(300, 1800, 76);
  let decomp = '', audit = { massBalanced: false };
  try { const bal = chem.balance([ox], [metal, 'O2']); decomp = chem.equationString([ox], [metal, 'O2'], bal.coefficients); audit = chem.auditEquation([ox], [metal, 'O2'], bal.coefficients); } catch (e) { warnings.push('Could not balance the oxide decomposition reaction: ' + e.message); }
  const P = chem.prettyFormula(ox), prime = '′'.repeat(z), dot = '•'.repeat(z);
  return {
    summary: `${P} at ${F(p.T_C)} °C and pO₂ = ${F(p.pO2, 3)} atm: oxygen-vacancy site fraction ${F(vO, 3)} with ${F(nE * 1e-6, 3)} electrons/cm³ (Brouwer slope ${F(slope, 3)}); thermal vacancy fraction in the reference metal ${F(cv, 3)}. A dislocation density of ${F(p.rho_d, 3)} m⁻² gives a Taylor flow stress of ${F(sig / 1e6, 4)} MPa, and McLean segregation enriches the boundaries ${F(beta, 3)}-fold to ${F(100 * xgb, 3)} at.%.`,
    kpis: [
      kpi('c_vac', 'Equilibrium vacancy fraction (metal)', cv, '–', `c_v = exp(S_f/k)·exp(−E_f/kT) with E_f = ${F(p.Evf)} eV. Near the melting point of metals this reaches about 10⁻⁴–10⁻³.`, null),
      kpi('VO', 'Oxygen-vacancy site fraction', vO, '–', `From the reduction equilibrium [V_O••]·n²·pO₂^½ = K_red with electroneutrality n = 2[V_O••] + [D•]; K_red = ${F(Kred, 3)}. Expressed per oxygen site of ${P}.`, vO < 0.02),
      kpi('n_e', 'Electron concentration', nE, '1/m³', `Free-electron density implied by the defect equilibrium (${F(nE * 1e-6, 3)} cm⁻³) using ${F(nFU, 3)} formula units per m³. It controls n-type conductivity.`, null),
      kpi('brouwer', 'Brouwer slope d log[V_O]/d log pO₂', slope, '–', 'Local slope of the Brouwer diagram: −1/6 in the intrinsic reduction regime (n = 2[V_O••]) and −1/2 where donors fix the electron density.', null),
      kpi('schottky', 'Schottky defect fraction', xS, '–', `Extent x of nil ⇌ ${a} V_${metal} + ${b} V_O from (ax)^a(bx)^b = exp(−ΔH_S/kT), entropy neglected.`, null),
      kpi('frenkel', 'Cation Frenkel defect fraction', xF, '–', 'Site fraction of interstitial–vacancy pairs, exp(−ΔH_F/2kT). The larger of the Schottky and Frenkel values identifies the dominant intrinsic disorder.', null),
      kpi('sigma_T', 'Taylor flow stress', sig, 'Pa', `σ = σ₀ + MαGb√ρ with Taylor factor M = 3.06, α = ${F(p.alpha_T)}; the dislocation term contributes ${F((sig - p.sigma0 * 1e6) / 1e6, 4)} MPa.`, null),
      kpi('rho_sat', 'Saturation dislocation density', rhoSat, '1/m²', 'Kocks–Mecking steady state (k₁/k₂)² where storage balances dynamic recovery; the flow stress saturates at the corresponding Taylor value.', null),
      kpi('E_stored', 'Stored dislocation energy', Est, 'J/m³', 'Line-energy estimate ½Gb²ρ: the driving pressure available for recovery and recrystallisation.', null),
      kpi('x_gb', 'Grain-boundary solute fraction', xgb, '–', `Langmuir–McLean isotherm with ΔG_seg = ${F(p.dG_seg)} kJ/mol: boundary sites are enriched by a factor ${F(beta, 3)} over the bulk.`, !depl),
    ],
    series: [
      { key: 'brouwer', title: 'Brouwer diagram', xlabel: 'log₁₀ pO₂ (atm)', ylabel: 'log₁₀ concentration (per formula unit)', lines: [line('[V_O••]', lp, bro.map(s => Math.log10(s.v))), line('n (electrons)', lp, bro.map(s => Math.log10(s.n)))], explain: 'Defect concentrations against oxygen partial pressure at the chosen temperature. A change of slope marks the transition from donor-controlled to reduction-controlled behaviour.' },
      { key: 'flow', title: 'Work hardening (Kocks–Mecking + Taylor)', xlabel: 'plastic strain', ylabel: 'flow stress (MPa)', lines: [line('σ(ε)', eps, sigE.map(v => v / 1e6))], explain: 'Flow stress as dislocation density evolves with strain. The curve saturates when storage and dynamic recovery balance.' },
      { key: 'seg', title: 'Grain-boundary segregation versus temperature', xlabel: 'temperature (K)', ylabel: 'boundary solute fraction', lines: [line('McLean isotherm', Ts, Ts.map(mcl))], explain: 'Equilibrium boundary coverage falls with temperature as entropy favours the bulk; low-temperature anneals therefore maximise segregation (if kinetics permit).' },
    ],
    tables: [{ key: 'defects', title: 'Intrinsic defect inventory', columns: ['Defect', 'Site fraction', 'Concentration (1/cm³)'], rows: [['V_O•• (reduction)', F(vO, 3), F(vO * nO * 1e-6, 3)], [`V_${metal} (Schottky)`, F(xS, 3), F(xS * nCat * 1e-6, 3)], [`${metal}_i (Frenkel)`, F(xF, 3), F(xF * nCat * 1e-6, 3)], ['electrons', F(here.n, 3), F(nE * 1e-6, 3)], ['donor D•', F(p.x_donor, 3), F(p.x_donor * nCat * 1e-6, 3)]], explain: 'Site fractions and volumetric concentrations of each defect species for the selected oxide, based on its handbook density and formula mass.' }],
    equations: [
      ...(decomp ? [eqn('chemical', 'Oxide decomposition (overall reduction)', decomp + '(g)', audit.massBalanced ? 'balanced by element' : 'audit failed')] : []),
      eqn('ionic', 'Reduction (Kröger–Vink)', 'O_O× ⇌ V_O•• + 2 e′ + ½ O₂(g)', 'K_red = [V_O••]·n²·pO₂^½'),
      eqn('ionic', 'Schottky disorder', `nil ⇌ ${a > 1 ? a + ' ' : ''}V_${metal}${prime} + ${b > 1 ? b + ' ' : ''}V_O••`),
      eqn('ionic', 'Cation Frenkel disorder', `${metal}_${metal}× ⇌ ${metal}_i${dot} + V_${metal}${prime}`),
      eqn('formula', 'Electroneutrality', 'n = 2[V_O••] + [D•]'),
      eqn('constitutive', 'Taylor hardening and Kocks–Mecking evolution', 'σ = σ₀ + M α G b √ρ;  dρ/dε = M(k₁√ρ − k₂ρ)'),
      eqn('formula', 'McLean segregation isotherm', 'X_gb/(1 − X_gb) = X_b/(1 − X_b)·exp(ΔG_seg/RT)'),
      eqn('formula', 'Thermal vacancies', 'c_v = exp(S_f/k_B)·exp(−E_f/k_BT)'),
    ],
    checks: [
      chk('Electroneutrality satisfied', Math.abs(here.n - 2 * here.v - dn) <= 1e-9 * here.n, `n − 2[V_O] − [D] = ${(here.n - 2 * here.v - dn).toExponential(2)}`),
      chk('Mass-action residual', Math.abs(Math.log(here.v * here.n ** 2 * Math.sqrt(p.pO2) / Kred)) < 1e-6, 'ln([V_O]n²pO₂^½/K) below 10⁻⁶'),
      chk('Dilute defect limit', vO < 0.02 && xS < 0.02 && xF < 0.02, `largest site fraction ${F(Math.max(vO, xS, xF), 3)}`),
      chk('Decomposition equation balanced', !!audit.massBalanced, decomp || 'not available'),
      chk('Solute inventory sufficient for segregation', !depl, `boundary demand ${F(need, 3)} vs bulk ${F(p.x_bulk, 3)}`),
      chk('Dislocation spacing exceeds core size', spacing > 5 * bb, `mean spacing ${F(spacing * 1e9, 3)} nm`),
    ],
    warnings,
    recommendations: [
      rec('high', `Fit ΔH_red and ΔS_red to conductivity or thermogravimetric data over at least four decades of pO₂ before using absolute defect densities; the −1/6 or −1/2 slope (here ${F(slope, 3)}) identifies the electroneutrality regime.`, 'Brouwer analysis.'),
      rec('medium', nE > 1e24 ? 'The electron density is high: anneal in oxygen or reduce donor content if a resistive or lightly doped oxide is required.' : 'Anneal in a reducing atmosphere or add a donor if higher n-type conductivity is required; the electron density scales as pO₂^(−1/6) without donors.', `n = ${F(nE * 1e-6, 3)} cm⁻³.`),
      rec('medium', p.rho_d > 0.5 * rhoSat ? 'Dislocation density is near saturation: further cold work adds little strength; use grain refinement or precipitates instead.' : `Cold work can still raise strength towards the saturation value of about ${F(tay(rhoSat) / 1e6, 4)} MPa.`, 'Kocks–Mecking saturation density.'),
      rec('medium', beta > 50 ? 'Strong segregation: exploit it to stabilise nanocrystalline grain size with a small solute addition, which avoids energy-intensive re-processing of coarsened material.' : 'Segregation is weak at this temperature; choose a solute with larger size or chemical mismatch if thermodynamic grain-size stabilisation is intended, so that lower-temperature, lower-energy processing remains possible.', `Enrichment factor ${F(beta, 3)}.`, true),
    ],
    publish: { 'prop.yield': sig, 'micro.dislocation_density': p.rho_d, 'micro.vacancy_fraction': vO, 'micro.electron_density': nE, 'material.formula': ox },
    numerics: { method: 'Mass-action defect equilibria with Brent root finding on electroneutrality; closed-form Kocks–Mecking and McLean relations', dof: 2, steps: 81, iterations: 81, residual: Math.abs(here.n - 2 * here.v - dn), notes: 'Ideal dilute solution; fully ionised defects; band-gap excitation and holes neglected.' },
  };
}

export default {
  id: 16,
  tagline: 'Simulates microstructure evolution — spinodal decomposition, grain growth, dendritic growth, coarsening, sintering — and equilibrium defect populations.',
  inputs,
  presets: [
    { name: 'Spinodal nanostructure, 50:50 alloy', note: 'Symmetric composition gives a bicontinuous structure, the classic route to nanoporous metals after dealloying.', values: { mode: 'spinodal', c0: 0.5, dx_nm: 1, kappa: 1.5, ch_steps: 160 } },
    { name: 'Off-symmetric blend, droplet morphology', note: 'A 35:65 composition inside the spinodal separates into droplets that coarsen by Ostwald ripening.', values: { mode: 'spinodal', c0: 0.35, dx_nm: 2, ch_steps: 200, noise: 0.05 } },
    { name: 'Nanocrystalline ZnO film with pinning particles', note: 'Grain growth of a nanocrystalline oxide film stabilised by 2 % second-phase particles.', values: { mode: 'grain', fp: 0.02, dx_nm: 2, mcs: 120, T_C: 500 } },
    { name: 'Hexagonal dendrite (six-fold)', note: 'Six-fold anisotropy typical of basal-plane growth of hexagonal crystals at high driving force.', values: { mode: 'solid', fold: '6', latent: 1.8, aniso: 0.02, dx_um: 0.5, alpha_th: 1e-9 } },
    { name: 'ZnO nanopowder anneal at 500 °C', note: 'Grain growth, coarsening and initial-stage sintering of 20 nm crystallites during calcination.', values: { mode: 'kinetics', T_C: 500, time_h: 2, d0_nm: 20, dp_nm: 40, sinter: 'gb' } },
    { name: 'ZnO oxygen vacancies, reducing anneal', note: 'Kröger–Vink equilibrium of ZnO at 700 °C in a low-oxygen atmosphere, with work-hardening and segregation estimates.', values: { mode: 'defects', oxide: 'ZnO', T_C: 700, pO2: 1e-10, x_donor: 0 } },
    { name: 'Lattice strut metal: dislocations and segregation', note: 'Copper-like strut material of a nanolattice metamaterial: Taylor hardening and boundary segregation in 50 nm grains.', values: { mode: 'defects', oxide: 'Al2O3', T_C: 300, rho_d: 5e14, G_shear: 4.8e10, dgrain_nm: 50 } },
  ],
  run(inp) {
    const p = read(inp);
    const r = p.mode === 'grain' ? grains(p) : p.mode === 'solid' ? solid(p) : p.mode === 'kinetics' ? kinetics(p) : p.mode === 'defects' ? defects(p) : spinodal(p);
    return finish(r);
  },
};
