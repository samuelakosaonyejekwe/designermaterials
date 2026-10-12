// Reference benchmarks for Suite 8 — surface, interface and colloidal engineering.
// Expected values are literals from standard texts or are evaluated here from closed-form
// expressions with CODATA constants; nothing below calls the suite.

const kB = 1.380649e-23, qe = 1.602176634e-19, NA = 6.02214076e23, E0 = 8.8541878128e-12, Rg = 8.314462618, g0 = 9.80665, PI = Math.PI;
const T25 = 298.15, epsW = 78.5, muW = 8.9e-4, kT = kB * T25;
const rad = d => d * PI / 180, deg = r => r * 180 / PI;

// Debye length (m) of a symmetric z:z electrolyte, concentration in mol/L
const debye = (c, z = 1, eps = epsW, T = T25) => Math.sqrt(eps * E0 * kB * T / (2 * z * z * qe * qe * NA * c * 1000));
// Grahame equation (C/m²)
const grahame = (psi, c, z = 1) => Math.sqrt(8 * epsW * E0 * kT * c * 1000 * NA) * Math.sinh(z * qe * psi / (2 * kT));

// DLVO barrier (kT) of two equal spheres from the closed-form potentials, located on a fine grid
const dlvoBarrier = ({ a, A, psi, c }) => {
  const k = 1 / debye(c); let best = -Infinity;
  for (let i = 0; i <= 200000; i++) {
    const h = 0.158e-9 * Math.pow(10, 3 * i / 200000);
    const vdw = -A / 6 * (2 * a * a / (h * (4 * a + h)) + 2 * a * a / (2 * a + h) ** 2 + Math.log(h * (4 * a + h) / (2 * a + h) ** 2));
    const edl = 2 * PI * epsW * E0 * a * psi * psi * Math.log(1 + Math.exp(-k * h));
    best = Math.max(best, (vdw + edl) / kT);
  }
  return best;
};

// exact synthetic isotherms (no noise) for parameter-recovery tests
const XS = [0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
const table = f => ['p_rel,q_mmol_g', ...XS.map(x => `${x},${f(x).toPrecision(12)}`)].join('\n');
const QM_STP = 1000 / 22413.97;                       // 1 cm³(STP)/g of gas in mmol/g
const betIso = table(x => QM_STP * 80 * x / ((1 - x) * (1 + 79 * x)));
const langIso = table(x => 2 * 5 * x / (1 + 5 * x));
const freuIso = table(x => 1.5 * Math.pow(x, 0.4));

// nanorod-array geometry helpers: number density giving a chosen hydraulic radius r_h = 2(1 − f)/(Nπd)
const nForRh = (d, rh) => 2 / (PI * d * (rh + d / 2));
const pitchGap = (N, d) => 1 / Math.sqrt(N) - d;

const stokesD = d => kT / (3 * PI * muW * d);

export default [
  // ------------------------------------------------------------------ double layer
  {
    id: 's08-debye-length-0p1M', suite: 8, title: 'Debye length of a 0.1 mol/L 1:1 electrolyte in water at 25 °C', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.1, zval: 1, T: 25, epsr: 78.5 }, kpi: 'debye',
    expected: debye(0.1) * 1e9, tol: 0.002, kind: 'analytical',
    reference: 'Debye–Hückel theory: κ⁻¹ = [εε₀kT/(2z²e²N_A c)]^½ (Israelachvili, Intermolecular and Surface Forces)',
    note: 'Evaluates to 0.962 nm with ε_r = 78.5.',
  },
  {
    id: 's08-debye-length-handbook-1M', suite: 8, title: 'Debye length rule 0.304 nm/√c at 1 mol/L (1:1, 25 °C)', mode: 'dlvo',
    inputs: { mode: 'dlvo', cel: 1, zval: 1, T: 25, epsr: 78.5 }, kpi: 'debye',
    expected: 0.304, tol: 0.01, kind: 'reference-data',
    reference: 'Israelachvili, Intermolecular and Surface Forces — κ⁻¹ = 0.304/√[NaCl] nm for 1:1 electrolytes at 25 °C',
    note: 'The tabulated coefficient is rounded to three figures and depends slightly on the permittivity adopted for water.',
  },
  {
    id: 's08-debye-length-2-2', suite: 8, title: 'Debye length of a 0.01 mol/L 2:2 electrolyte', mode: 'dlvo',
    inputs: { mode: 'dlvo', cel: 0.01, zval: 2, T: 25, epsr: 78.5 }, kpi: 'debye',
    expected: 1.52, tol: 0.01, kind: 'reference-data',
    reference: 'Israelachvili, Intermolecular and Surface Forces — κ⁻¹ = 0.152/√c nm for 2:2 electrolytes at 25 °C',
    note: 'Rounded handbook coefficient.',
  },
  {
    id: 's08-grahame-surface-charge', suite: 8, title: 'Grahame equation: surface charge at −100 mV in 0.1 mol/L 1:1 electrolyte', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.1, zval: 1, psi0: -100, dStern: 0, T: 25, epsr: 78.5 }, kpi: 'sigma',
    expected: grahame(-0.1, 0.1) * 1e3, tol: 0.002, kind: 'analytical',
    reference: 'Grahame equation σ = (8εε₀kT c N_A)^½ sinh(zeψ₀/2kT) (Israelachvili; Hunter, Foundations of Colloid Science)',
    note: 'No Stern layer, so the diffuse-layer potential equals the surface potential. About −127 mC/m², in line with σ = 0.117 sinh(ψ₀/51.4 mV)√c.',
  },
  {
    id: 's08-diffuse-layer-capacitance', suite: 8, title: 'Gouy–Chapman diffuse-layer capacitance at 0.1 mol/L near the point of zero charge', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.1, zval: 1, psi0: 1, dStern: 0, T: 25, epsr: 78.5 }, kpi: 'cap',
    expected: 228 * Math.sqrt(0.1) * Math.cosh(19.5 * 0.001), tol: 0.01, kind: 'reference-data',
    reference: 'Bard & Faulkner, Electrochemical Methods — C_d = 228 z c^½ cosh(19.5 z φ₀) µF/cm² for dilute aqueous solutions at 25 °C',
    note: 'Rounded textbook coefficients; gives 72.1 µF/cm².',
  },
  {
    id: 's08-stern-series-capacitance', suite: 8, title: 'Gouy–Chapman–Stern capacitance: Stern and diffuse layers in series (low potential)', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.1, zval: 1, psi0: 1, dStern: 0.3, epsS: 10, T: 25, epsr: 78.5 }, kpi: 'cap',
    expected: 100 / (debye(0.1) / (epsW * E0) + 0.3e-9 / (10 * E0)), tol: 0.002, kind: 'analytical',
    reference: 'Stern model: 1/C = d/(ε_S ε₀) + 1/(εε₀κ) in the Debye–Hückel limit (Bard & Faulkner, Electrochemical Methods)',
    note: 'At 1 mV the cosh factor of the diffuse capacitance differs from one by less than 0.01 %.',
  },
  {
    id: 's08-stern-potential-divider', suite: 8, title: 'Potential at the outer Helmholtz plane from the capacitive divider (low potential)', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.1, zval: 1, psi0: 1, dStern: 0.3, epsS: 10, T: 25, epsr: 78.5 }, kpi: 'psi_d',
    expected: 1 / (1 + (epsW * E0 / debye(0.1)) * 0.3e-9 / (10 * E0)), tol: 0.002, kind: 'analytical',
    reference: 'Stern model in the Debye–Hückel limit: ψ_d = ψ₀ C_Stern/(C_Stern + C_diffuse)',
    note: 'Linearised Poisson–Boltzmann; the nonlinear correction at 1 mV is negligible.',
  },
  {
    id: 's08-boltzmann-counterion-enrichment', suite: 8, title: 'Boltzmann enrichment: tenfold counter-ion concentration at ψ = −(kT/e) ln 10', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.01, zval: 1, psi0: -(kT / qe) * Math.LN10 * 1e3, dStern: 0, T: 25 }, kpi: 'c_wall',
    expected: 0.1, tol: 1e-4, kind: 'analytical',
    reference: 'Boltzmann distribution of ions in the diffuse layer, c = c∞ exp(−zeψ/kT)',
    note: 'The potential of 59.16 mV is the Nernstian decade at 25 °C.',
  },
  {
    id: 's08-gouy-chapman-solver', suite: 8, title: 'Nonlinear Poisson–Boltzmann finite-difference profile against the Gouy–Chapman solution', mode: 'edl',
    inputs: { mode: 'edl', cel: 0.01, zval: 1, psi0: -100, dStern: 0, T: 25 }, kpi: 'gc_err',
    expected: 0, abs: 2e-3, kind: 'analytical',
    reference: 'Gouy–Chapman solution tanh(y/4) = tanh(y₀/4) e^(−κx) of the planar Poisson–Boltzmann equation',
    note: 'Maximum deviation of the numerical profile, as a fraction of the wall potential (zeψ/kT ≈ 3.9, strongly nonlinear).',
  },
  {
    id: 's08-smoluchowski-mobility', suite: 8, title: 'Electrophoretic mobility in the thin-double-layer (Smoluchowski) limit', mode: 'edl',
    inputs: { mode: 'edl', dp: 1e-4, cel: 0.1, zval: 1, psi0: -25, dStern: 0, T: 25, epsr: 78.5, mu: 8.9e-4 }, kpi: 'mob',
    expected: epsW * E0 * (-0.025) / muW * 1e8, tol: 0.002, kind: 'analytical',
    reference: 'Smoluchowski equation μ_e = εε₀ζ/η for κa ≫ 1 (Hunter, Zeta Potential in Colloid Science)',
    note: 'κa ≈ 5×10⁴, where the Henry function is within 0.01 % of 3/2.',
  },
  {
    id: 's08-huckel-mobility', suite: 8, title: 'Electrophoretic mobility in the thick-double-layer (Hückel) limit', mode: 'edl',
    inputs: { mode: 'edl', dp: 1e-9, cel: 1e-6, zval: 1, psi0: -25, dStern: 0, T: 25, epsr: 78.5, mu: 8.9e-4 }, kpi: 'mob',
    expected: 2 * epsW * E0 * (-0.025) / (3 * muW) * 1e8, tol: 0.002, kind: 'analytical',
    reference: 'Hückel equation μ_e = 2εε₀ζ/(3η) for κa ≪ 1 (Hunter, Zeta Potential in Colloid Science)',
    note: 'κa ≈ 0.002.',
  },
  // ------------------------------------------------------------------ DLVO
  {
    id: 's08-dlvo-barrier', suite: 8, title: 'DLVO energy barrier of 100 nm spheres (A = 20 zJ, −30 mV, 0.01 mol/L)', mode: 'dlvo',
    inputs: { mode: 'dlvo', dp: 1e-7, AH: 20, psi0: -30, cel: 0.01, zval: 1, T: 25, epsr: 78.5 }, kpi: 'barrier',
    expected: dlvoBarrier({ a: 5e-8, A: 20e-21, psi: -0.03, c: 0.01 }), tol: 0.005, kind: 'analytical',
    reference: 'DLVO theory: Hamaker sphere–sphere attraction plus Hogg–Healy–Fuerstenau constant-potential repulsion (Hunter, Foundations of Colloid Science)',
    note: 'Reference maximum located independently on a 200 000-point grid of the closed-form potentials; the suite uses a 400-point logarithmic grid.',
  },
  {
    id: 's08-schulze-hardy-low-potential', suite: 8, title: 'Schulze–Hardy exponent in the low-potential limit (ccc ∝ z⁻²)', mode: 'dlvo',
    inputs: { mode: 'dlvo', dp: 1e-7, AH: 20, psi0: -30, cel: 0.01, T: 25 }, kpi: 'sh',
    expected: -2, abs: 0.02, kind: 'analytical',
    reference: 'DLVO theory at fixed low surface potential: the barrier depends on electrolyte only through κ ∝ z√c, so ccc ∝ z⁻² (Hunter, Foundations of Colloid Science)',
    note: 'The z⁻⁶ law applies only at high potential, where the nonlinear Gouy–Chapman term saturates; the suite uses the linearised repulsion.',
  },
  {
    id: 's08-smoluchowski-rapid-half-time', suite: 8, title: 'Smoluchowski rapid-coagulation half-time 3μ/(4kT n₀)', mode: 'dlvo',
    inputs: { mode: 'dlvo', dp: 1e-7, AH: 0.1, psi0: 0, phi: 1e-4, T: 25, mu: 8.9e-4 }, kpi: 't_half',
    expected: 3 * muW / (4 * kT * (1e-4 / (PI / 6 * 1e-21))), tol: 0.01, kind: 'analytical',
    reference: 'Smoluchowski theory of diffusion-limited (rapid) coagulation, t_½ = 3μ/(4kT n₀)',
    note: 'Uncharged spheres with the smallest admissible Hamaker constant (0.1 zJ); the residual attraction makes the Fuchs factor about 0.5 % below one.',
  },
  // ------------------------------------------------------------------ adsorption
  {
    id: 's08-bet-nitrogen-area', suite: 8, title: 'BET area for a monolayer of 1 cm³(STP)/g of nitrogen', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: betIso, sigN: 0.162 }, kpi: 'ssa',
    expected: 4.35, tol: 0.005, kind: 'reference-data',
    reference: 'Gregg & Sing, Adsorption, Surface Area and Porosity — 4.35 m² per cm³(STP) of N₂ monolayer with σ = 0.162 nm²',
    note: 'Exact synthetic BET isotherm (q_m = 1 cm³(STP)/g, C = 80) generated in this file.',
  },
  {
    id: 's08-bet-monolayer-recovery', suite: 8, title: 'BET fit recovers the monolayer capacity of an exact BET isotherm', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: betIso }, kpi: 'qm',
    expected: QM_STP, tol: 1e-3, kind: 'analytical',
    reference: 'Brunauer–Emmett–Teller isotherm q = q_m C x/[(1 − x)(1 + (C − 1)x)]',
    note: 'Noise-free data, so slope and intercept of the BET transform are exact.',
  },
  {
    id: 's08-bet-constant-recovery', suite: 8, title: 'BET fit recovers the constant C of an exact BET isotherm', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: betIso }, kpi: 'C',
    expected: 80, tol: 1e-3, kind: 'analytical',
    reference: 'Brunauer–Emmett–Teller isotherm, C = 1 + slope/intercept of the linear transform',
  },
  {
    id: 's08-langmuir-recovery', suite: 8, title: 'Langmuir fit recovers the capacity of an exact Langmuir isotherm', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: langIso }, kpi: 'qmL',
    expected: 2, tol: 1e-3, kind: 'analytical',
    reference: 'Langmuir isotherm q = q_m K x/(1 + K x); synthetic data with q_m = 2 mmol/g and K = 5',
  },
  {
    id: 's08-langmuir-half-coverage', suite: 8, title: 'Langmuir coverage of one half at x = 1/K', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: langIso, xq: 0.2 }, kpi: 'theta',
    expected: 0.5, tol: 1e-3, kind: 'analytical',
    reference: 'Langmuir isotherm: θ = Kx/(1 + Kx) = ½ at x = 1/K',
  },
  {
    id: 's08-freundlich-exponent-recovery', suite: 8, title: 'Freundlich fit recovers the exponent of an exact power-law isotherm', mode: 'adsorb',
    inputs: { mode: 'adsorb', iso: freuIso }, kpi: 'nF',
    expected: 0.4, tol: 1e-3, kind: 'analytical',
    reference: 'Freundlich isotherm q = K_F x^(1/n); synthetic data with K_F = 1.5 and 1/n = 0.4',
  },
  // ------------------------------------------------------------------ wetting
  {
    id: 's08-wenzel-angle', suite: 8, title: 'Wenzel contact angle for roughness r = 1 + π/10 and θ_Y = 60°', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 60, dRod: 1e-7, LRod: 1e-6, Nrod: 1e12 }, kpi: 'thW',
    expected: deg(Math.acos((1 + PI * 0.1) * Math.cos(rad(60)))), tol: 1e-4, kind: 'analytical',
    reference: 'Wenzel relation cos θ* = r cos θ_Y (de Gennes, Brochard-Wyart & Quéré, Capillarity and Wetting Phenomena)',
    note: 'Cylindrical posts: r = 1 + NπdL.',
  },
  {
    id: 's08-cassie-baxter-angle', suite: 8, title: 'Cassie–Baxter angle on a 10 % solid fraction with θ_Y = 110°', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 110, dRod: 1e-7, LRod: 1e-6, Nrod: 0.4 / (PI * 1e-14) }, kpi: 'thCB',
    expected: deg(Math.acos(0.1 * (1 + Math.cos(rad(110))) - 1)), tol: 1e-4, kind: 'analytical',
    reference: 'Cassie–Baxter relation cos θ* = f(1 + cos θ_Y) − 1 (de Gennes, Brochard-Wyart & Quéré, Capillarity and Wetting Phenomena)',
    note: 'About 159°: a superhydrophobic composite state.',
  },
  {
    id: 's08-hemiwicking-critical-angle', suite: 8, title: 'Critical contact angle for hemiwicking in a post array', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 60, dRod: 1e-7, LRod: 1e-6, Nrod: 1e12 }, kpi: 'thC',
    expected: deg(Math.acos((1 - PI * 0.0025) / (1 + PI * 0.1 - PI * 0.0025))), tol: 1e-4, kind: 'analytical',
    reference: 'Hemiwicking criterion cos θ_c = (1 − f)/(r − f) (de Gennes, Brochard-Wyart & Quéré, Capillarity and Wetting Phenomena)',
  },
  {
    id: 's08-young-laplace-capillary-pressure', suite: 8, title: 'Young–Laplace pressure 2γ/r for water in a 1 µm hydraulic radius', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 0, gLV: 0.072, dRod: 1e-6, LRod: 1e-6, Nrod: nForRh(1e-6, 1e-6) }, kpi: 'Pcap',
    expected: 0.144, tol: 1e-4, kind: 'analytical',
    reference: 'Young–Laplace equation ΔP = 2γ cos θ/r',
    note: 'Post density chosen so that the hydraulic radius 2(1 − f)/(Nπd) of the inter-rod space is exactly 1 µm.',
  },
  {
    id: 's08-kelvin-capillary-condensation', suite: 8, title: 'Kelvin equation: condensation humidity of water in a 10 nm radius meniscus at 25 °C', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 0, gLV: 0.072, Vml: 1.8e-5, T: 25, dRod: 2e-8, LRod: 1e-7, Nrod: nForRh(2e-8, 1e-8) }, kpi: 'RHc',
    expected: 100 * Math.exp(-2 * 0.072 * 1.8e-5 / (1e-8 * Rg * T25)), tol: 1e-3, kind: 'analytical',
    reference: 'Kelvin equation ln(p/p₀) = −2γV_m/(rRT) for a concave meniscus (Israelachvili, Intermolecular and Surface Forces)',
    note: 'About 90 % relative humidity; γV_m/RT ≈ 0.52 nm for water.',
  },
  {
    id: 's08-young-dupre-work-of-adhesion', suite: 8, title: 'Young–Dupré work of adhesion of water at θ_Y = 60°', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 60, gLV: 0.072 }, kpi: 'Wad',
    expected: 108, tol: 1e-4, kind: 'analytical',
    reference: 'Young–Dupré equation W = γ(1 + cos θ)',
  },
  {
    id: 's08-elastocapillary-critical-length', suite: 8, title: 'Critical rod length for capillary collapse from cantilever theory', mode: 'wetting',
    inputs: { mode: 'wetting', thetaY: 50, gLV: 0.072, dRod: 8e-8, LRod: 1.5e-6, Nrod: 4e13, Emod: 1.4e11 }, kpi: 'Lcrit',
    expected: (() => { const d = 8e-8, gap = pitchGap(4e13, d), q = 2 * 0.072 * Math.cos(rad(50)) / gap * d, EI = 1.4e11 * PI * d ** 4 / 64; return Math.pow((gap / 2) * 8 * EI / q, 0.25) * 1e6; })(),
    tol: 1e-3, kind: 'analytical',
    reference: 'Euler–Bernoulli cantilever under uniform load, tip deflection qL⁴/(8EI), with the Laplace pressure 2γ cos θ/g acting on the rod width (Timoshenko, Strength of Materials)',
    note: 'Collapse criterion: tip deflection equal to half the gap. A reduced-order worst-case estimate, benchmarked here only for internal consistency with beam theory.',
  },
  // ------------------------------------------------------------------ transport
  {
    id: 's08-stokes-einstein-100nm', suite: 8, title: 'Stokes–Einstein diffusivity of a 100 nm sphere in water at 25 °C', mode: 'transport',
    inputs: { mode: 'transport', dp: 1e-7, mu: 8.9e-4, T: 25 }, kpi: 'D',
    expected: stokesD(1e-7), tol: 1e-4, kind: 'analytical',
    reference: 'Stokes–Einstein relation D = kT/(3πμd)',
    note: '4.91×10⁻¹² m²/s.',
  },
  {
    id: 's08-hydrodynamic-diameter', suite: 8, title: 'Hydrodynamic diameter from a measured diffusivity of 2.0×10⁻¹¹ m²/s', mode: 'transport',
    inputs: { mode: 'transport', Dmeas: 2e-11, mu: 8.9e-4, T: 25 }, kpi: 'dH',
    expected: kT / (3 * PI * muW * 2e-11) * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Stokes–Einstein relation inverted, d_H = kT/(3πμD), as used in dynamic light scattering',
    note: '24.5 nm.',
  },
  {
    id: 's08-stokes-settling-1um', suite: 8, title: 'Stokes settling velocity of a 1 µm ZnO sphere in water', mode: 'transport',
    inputs: { mode: 'transport', dp: 1e-6, rhoP: 5606, rhoF: 997, mu: 8.9e-4, RCF: 1 }, kpi: 'v',
    expected: (5606 - 997) * g0 * 1e-12 / (18 * muW) * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Stokes law v = Δρ g d²/(18μ)',
    note: '2.82 µm/s.',
  },
  {
    id: 's08-centrifugal-settling', suite: 8, title: 'Stokes settling of a 50 nm gold sphere at 10 000 g', mode: 'transport',
    inputs: { mode: 'transport', dp: 5e-8, rhoP: 19300, rhoF: 997, mu: 8.9e-4, RCF: 1e4 }, kpi: 'v',
    expected: (19300 - 997) * g0 * 1e4 * 25e-16 / (18 * muW) * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Stokes law with centrifugal acceleration, v = Δρ (RCF·g) d²/(18μ)',
  },
  {
    id: 's08-gravitational-length', suite: 8, title: 'Sedimentation (barometric) length kT/(Δm g) of a 100 nm ZnO sphere', mode: 'transport',
    inputs: { mode: 'transport', dp: 1e-7, rhoP: 5606, rhoF: 997, T: 25, RCF: 1 }, kpi: 'lg',
    expected: kT / (PI / 6 * 1e-21 * (5606 - 997) * g0), tol: 1e-4, kind: 'analytical',
    reference: 'Perrin sedimentation–diffusion equilibrium c(z) = c(0) exp(−Δm g z/kT)',
  },
  {
    id: 's08-peclet-unity-crossover', suite: 8, title: 'Péclet number of one at the crossover diameter (12kT/πΔρg)^¼', mode: 'transport',
    inputs: { mode: 'transport', dp: Math.pow(12 * kT / (PI * (5606 - 997) * g0), 0.25), rhoP: 5606, rhoF: 997, mu: 8.9e-4, T: 25, RCF: 1 }, kpi: 'Pe',
    expected: 1, tol: 1e-4, kind: 'analytical',
    reference: 'Stokes settling and Stokes–Einstein diffusion combined: Pe = v a/D = πΔρ g d⁴/(12kT)',
  },
];
