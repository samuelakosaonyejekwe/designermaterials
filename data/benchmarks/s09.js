// Reference benchmarks for Suite 9 — fluid dynamics and porous-media transport.
// Expected values are literals from standard texts or are evaluated here from closed-form
// expressions with standard constants; nothing below calls the suite.

const PI = Math.PI, kB = 1.380649e-23, g0 = 9.80665;
const rhoW = 997, muW = 8.9e-4;                        // water at 25 °C
const rhoA = 1.205, muA = 1.81e-5, lamA = 6.6e-8;      // air at 20 °C, 101.3 kPa

// packed bed (Ergun) pressure drop in Pa
const ergun = ({ U, e, d, L, rho = rhoW, mu = muW }) => (150 * mu * (1 - e) ** 2 * U / (e ** 3 * d * d) + 1.75 * rho * (1 - e) * U * U / (e ** 3 * d)) * L;

// channel helpers
const area = d => PI * d * d / 4;
const qForRe = (Re, d, rho = rhoW, mu = muW) => Re * mu / (rho * d) * area(d);
const dpFromF = (f, Re, d, L, rho = rhoW, mu = muW) => { const U = Re * mu / (rho * d); return f * (L / d) * rho * U * U / 2; };
// Taylor–Aris case used below: default 1 mm tube, slow diffusing species replaced by a fast one so that Dτ/a² ≫ 1
const TA = (() => { const d = 1e-3, L = 5, Q = 1e-8, D = 1e-8, U = Q / area(d), Dax = D + U * U * d * d / (192 * D), Pe = U * L / Dax, s2 = 2 / Pe - 2 / (Pe * Pe) * (1 - Math.exp(-Pe)); return { d, L, Q, D, Dax, Pe, s2 }; })();
// laminar segregated-flow case: default tube, k chosen for kτ = 1
const TAU_DEF = 5 * area(1e-3) / 1e-7;
const E1_HALF = 0.5597736;                             // exponential integral E₁(0.5)

// stirred tank helpers
const RUSH = { D: 0.1, T: 0.3, N: 5 };                 // Rushton turbine, 300 rpm
const P_RUSH = 5.0 * rhoW * RUSH.N ** 3 * RUSH.D ** 5;

// particle helpers
const stokes = (d, rhoP, rho, mu) => (rhoP - rho) * g0 * d * d / (18 * mu);

export default [
  // ------------------------------------------------------------------ packed bed
  {
    id: 's09-kozeny-carman-permeability', suite: 9, title: 'Kozeny–Carman permeability of 100 µm grains at porosity 0.25', mode: 'bed',
    inputs: { mode: 'bed', por: 0.25, dg: 1e-4 }, kpi: 'K',
    expected: 0.25 ** 3 * 1e-8 / (180 * 0.75 ** 2), tol: 1e-4, kind: 'analytical',
    reference: 'Kozeny–Carman equation K = ε³d²/(180(1 − ε)²) (Bear, Dynamics of Fluids in Porous Media)',
    note: '1.54×10⁻¹² m², about 1.6 darcy.',
  },
  {
    id: 's09-ergun-viscous-regime', suite: 9, title: 'Ergun pressure drop in creeping flow (Blake–Kozeny term dominant)', mode: 'bed',
    inputs: { mode: 'bed', por: 0.25, dg: 1e-4, Lbed: 0.1, Abed: 1e-3, Q: 1e-7, rho: 997, mu: 8.9e-4 }, kpi: 'dP',
    expected: ergun({ U: 1e-4, e: 0.25, d: 1e-4, L: 0.1 }) / 1e3, tol: 1e-4, kind: 'analytical',
    reference: 'Ergun equation ΔP/L = 150μ(1 − ε)²U/(ε³d²) + 1.75ρ(1 − ε)U²/(ε³d) (Bird, Stewart & Lightfoot, Transport Phenomena)',
    note: 'Particle Reynolds number 0.011; 4.81 kPa over 0.1 m.',
  },
  {
    id: 's09-ergun-inertial-regime', suite: 9, title: 'Ergun pressure drop with dominant inertial (Burke–Plummer) term', mode: 'bed',
    inputs: { mode: 'bed', por: 0.4, dg: 5e-3, Lbed: 1, Abed: 1e-3, Q: 5e-5, rho: 997, mu: 8.9e-4 }, kpi: 'dP',
    expected: ergun({ U: 0.05, e: 0.4, d: 5e-3, L: 1 }) / 1e3, tol: 1e-4, kind: 'analytical',
    reference: 'Ergun equation (Bird, Stewart & Lightfoot, Transport Phenomena)',
    note: 'Water at 0.05 m/s through 5 mm spheres, particle Reynolds number about 280: the inertial term supplies about 85 % of the 9.7 kPa.',
  },
  {
    id: 's09-bruggeman-diffusivity', suite: 9, title: 'Bruggeman effective diffusivity D_eff = ε^1.5 D', mode: 'bed',
    inputs: { mode: 'bed', por: 0.25, Dmol: 7.5e-10 }, kpi: 'Deff',
    expected: Math.pow(0.25, 1.5) * 7.5e-10, tol: 1e-4, kind: 'analytical',
    reference: 'Bruggeman relation for the effective diffusivity of a porous medium, D_eff = ε^(3/2) D',
  },
  {
    id: 's09-bed-plug-flow-conversion', suite: 9, title: 'First-order conversion of a long packed bed approaches plug flow, 1 − exp(−kτ)', mode: 'bed',
    inputs: { mode: 'bed', por: 0.25, dg: 1e-4, Lbed: 0.1, Abed: 1e-3, Q: 1e-7, k1: 1e-3 }, kpi: 'X',
    expected: 1 - Math.exp(-0.25), tol: 2e-3, kind: 'analytical',
    reference: 'Plug-flow reactor design equation for a first-order reaction, X = 1 − exp(−Da) (Levenspiel, Chemical Reaction Engineering)',
    note: 'Residence time εLA/Q = 250 s, Da = 0.25. The suite uses the axial-dispersion model at a Péclet number near 2000, which lowers the conversion by about 0.01 %.',
  },
  // ------------------------------------------------------------------ 2-D Darcy flow
  {
    id: 's09-darcy2d-uniform-flow-rate', suite: 9, title: 'Uniform-permeability 2-D Darcy solution reproduces the 1-D law Q = KAΔP/(μL)', mode: 'darcy2d',
    inputs: { mode: 'darcy2d', sigK: 0, Kg: 1e-12, Lx: 1, Ly: 1, thick: 0.01, dP: 1e5, mu: 8.9e-4 }, kpi: 'Q',
    expected: 1e-12 * (1 * 0.01) * 1e5 / (8.9e-4 * 1), tol: 1e-6, kind: 'analytical',
    reference: 'Darcy\'s law for one-dimensional flow through a homogeneous slab (Bear, Dynamics of Fluids in Porous Media)',
    note: 'A linear pressure field is represented exactly by the finite-volume scheme; the residual is the conjugate-gradient tolerance.',
  },
  {
    id: 's09-darcy2d-uniform-keff', suite: 9, title: 'Effective permeability of a homogeneous field equals the cell permeability (fine grid, elongated slab)', mode: 'darcy2d',
    inputs: { mode: 'darcy2d', sigK: 0, Kg: 3e-13, Lx: 2, Ly: 0.5, grid: 'fine' }, kpi: 'Krat',
    expected: 1, tol: 1e-6, kind: 'conservation',
    reference: 'Darcy\'s law: K_eff = QμL/(AΔP) for a homogeneous medium',
  },
  {
    id: 's09-darcy2d-uniform-outflow', suite: 9, title: 'Outlet flow-uniformity index of one for a homogeneous field', mode: 'darcy2d',
    inputs: { mode: 'darcy2d', sigK: 0 }, kpi: 'unif',
    expected: 1, tol: 1e-6, kind: 'conservation',
    reference: 'Symmetry of one-dimensional Darcy flow: identical flux through every outlet cell',
  },
  {
    id: 's09-darcy2d-plug-flow-conversion', suite: 9, title: 'Reactive transport in uniform Darcy flow against plug-flow conversion', mode: 'darcy2d',
    inputs: { mode: 'darcy2d', sigK: 0, grid: 'fine', Kg: 1e-12, por: 0.25, Lx: 1, Ly: 1, dP: 1e5, mu: 8.9e-4, k1: 5e-4, Dmol: 7.5e-10 }, kpi: 'X',
    expected: 1 - Math.exp(-5e-4 * 0.25 * 1 / (1e-12 * 1e5 / (8.9e-4 * 1))), tol: 0.01, kind: 'analytical',
    reference: 'Plug-flow reactor with first-order consumption, X = 1 − exp(−kεL/U) (Levenspiel, Chemical Reaction Engineering)',
    note: 'Da = 1.11 and Pe ≈ 10⁶. The first-order upwind scheme on 72 cells behaves as 72 tanks in series, which lowers the conversion by about 0.4 %.',
  },
  {
    id: 's09-darcy2d-species-conservation', suite: 9, title: 'Species balance of the heterogeneous reactive-transport solution', mode: 'darcy2d',
    inputs: { mode: 'darcy2d' }, kpi: 'mb',
    expected: 0, abs: 1e-6, kind: 'conservation',
    reference: 'Steady species conservation: inflow = outflow + amount reacted',
    note: 'Default log-normal permeability field (σ_lnK = 1).',
  },
  {
    id: 's09-darcy2d-lognormal-geometric-mean', suite: 9, title: 'Ensemble-mean effective permeability of a 2-D isotropic log-normal field against the geometric mean', mode: 'darcy2d',
    inputs: { mode: 'darcy2d', sigK: 1, corr: 3, seed: 1, grid: 'medium', nens: 16 }, kpi: 'Krat_ens',
    expected: 1, tol: 0.03, kind: 'analytical',
    reference: 'Matheron\'s exact result for isotropic two-dimensional log-normal media: K_eff equals the geometric mean (Dagan, Flow and Transport in Porous Formations)',
    note: 'Mean over 16 independent seeded 48 × 32 realisations of a statistically isotropic field (σ_lnK = 1, correlation length about 1.4 cells of the finer spacing); the ensemble standard error is 1.0 % (KPI Krat_se), so the tolerance is three standard errors. Single realisations scatter by about 4 %. The harmonic face averaging lowers the mean by roughly 1 % at this correlation length.',
  },
  // ------------------------------------------------------------------ channel flow
  {
    id: 's09-hagen-poiseuille', suite: 9, title: 'Hagen–Poiseuille pressure drop of water in a 1 mm tube', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: 1e-3, Lch: 5, Q: 1e-7, mu: 8.9e-4, rho: 997 }, kpi: 'dP',
    expected: 128 * 8.9e-4 * 5 * 1e-7 / (PI * 1e-12) / 1e3, tol: 1e-4, kind: 'analytical',
    reference: 'Hagen–Poiseuille law ΔP = 128μLQ/(πd⁴)',
    note: 'Re = 143; 18.13 kPa.',
  },
  {
    id: 's09-laminar-wall-shear', suite: 9, title: 'Wall shear stress of Poiseuille flow, 32μQ/(πd³)', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: 1e-3, Lch: 5, Q: 1e-7, mu: 8.9e-4, rho: 997 }, kpi: 'tw',
    expected: 32 * 8.9e-4 * 1e-7 / (PI * 1e-9), tol: 1e-4, kind: 'analytical',
    reference: 'Poiseuille flow: τ_w = 8μU/d = 32μQ/(πd³)',
  },
  {
    id: 's09-moody-smooth-pipe', suite: 9, title: 'Turbulent smooth-pipe friction at Re = 10⁵ (f = 0.0180)', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: 0.05, Lch: 10, rough: 0, Q: qForRe(1e5, 0.05), mu: 8.9e-4, rho: 997 }, kpi: 'dP',
    expected: dpFromF(0.0180, 1e5, 0.05, 10) / 1e3, tol: 0.01, kind: 'reference-data',
    reference: 'Moody chart / Prandtl–von Kármán smooth-pipe law 1/√f = 2 log₁₀(Re√f) − 0.8, giving f = 0.0180 at Re = 10⁵ (White, Fluid Mechanics)',
    note: 'The friction factor is quoted to three significant figures.',
  },
  {
    id: 's09-fully-rough-pipe', suite: 9, title: 'Fully rough turbulent friction, ε/d = 0.01 (von Kármán)', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: 0.5, Lch: 100, rough: 0.005, Q: qForRe(1e7, 0.5), mu: 8.9e-4, rho: 997 }, kpi: 'dP',
    expected: dpFromF(Math.pow(2 * Math.log10(3.7 / 0.01), -2), 1e7, 0.5, 100) / 1e3, tol: 0.005, kind: 'analytical',
    reference: 'von Kármán fully rough law 1/√f = 2 log₁₀(3.7d/ε) (White, Fluid Mechanics)',
    note: 'f = 0.0379. At Re = 10⁷ the Colebrook equation used by the suite still carries a viscous contribution of about 0.2 %.',
  },
  {
    id: 's09-taylor-aris-dispersion', suite: 9, title: 'Taylor–Aris axial dispersion coefficient D + U²d²/(192D)', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: TA.d, Lch: TA.L, Q: TA.Q, Dmol: TA.D }, kpi: 'Dax',
    expected: TA.Dax, tol: 1e-4, kind: 'analytical',
    reference: 'Taylor–Aris dispersion in laminar tube flow, D_ax = D + a²U²/(48D)',
  },
  {
    id: 's09-tanks-in-series-variance', suite: 9, title: 'Equivalent tanks in series from the closed-vessel dispersion variance, N = 1/σ_θ²', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: TA.d, Lch: TA.L, Q: TA.Q, Dmol: TA.D }, kpi: 'Nt',
    expected: 1 / TA.s2, tol: 1e-4, kind: 'analytical',
    reference: 'Tanks-in-series model σ_θ² = 1/N and closed-vessel dispersion model σ_θ² = 2/Pe − 2(1 − e^(−Pe))/Pe² (Levenspiel, Chemical Reaction Engineering)',
    note: 'Péclet number from the Taylor–Aris coefficient evaluated independently in this file.',
  },
  {
    id: 's09-laminar-flow-reactor-conversion', suite: 9, title: 'First-order conversion in segregated laminar tube flow at kτ = 1', mode: 'channel',
    inputs: { mode: 'channel', geo: 'pipe', dh: 1e-3, Lch: 5, Q: 1e-7, Dmol: 7.5e-10, k1: 1 / TAU_DEF }, kpi: 'X',
    expected: 1 - (0.25 * E1_HALF + 0.5 * Math.exp(-0.5)), tol: 2e-3, kind: 'analytical',
    reference: 'Laminar-flow reactor with negligible diffusion: 1 − X = y²E₁(y) + (1 − y)e^(−y), y = kτ/2 (Levenspiel, Chemical Reaction Engineering)',
    note: 'X = 0.557, between the stirred-tank (0.500) and plug-flow (0.632) values. Dτ/a² = 0.12, so the suite selects the segregated-flow model; it integrates the residence-time distribution by Simpson\'s rule.',
  },
  // ------------------------------------------------------------------ stirred vessel
  {
    id: 's09-rushton-power', suite: 9, title: 'Shaft power of a Rushton turbine in turbulent flow (N_p = 5)', mode: 'stirred',
    inputs: { mode: 'stirred', imp: 'rushton', Dimp: RUSH.D, Tt: RUSH.T, Nrpm: RUSH.N * 60, rho: 997, mu: 8.9e-4 }, kpi: 'P',
    expected: P_RUSH, tol: 0.01, kind: 'reference-data',
    reference: 'P = N_p ρN³D⁵ with the turbulent power number of about 5 for a six-blade Rushton turbine in a baffled tank (Paul, Atiemo-Obeng & Kresta, Handbook of Industrial Mixing)',
    note: 'Re = 5.6×10⁴. Measured power numbers range from about 4.8 to 5.5 with blade thickness and scale; the comparison is against the nominal value.',
  },
  {
    id: 's09-zwietering-just-suspended', suite: 9, title: 'Zwietering just-suspended speed for 20 µm ZnO agglomerates', mode: 'stirred',
    inputs: { mode: 'stirred', imp: 'pbt', Dimp: 0.05, Tt: 0.15, dpS: 2e-5, rhoS: 5606, Xs: 1, Szw: 6, rho: 997, mu: 8.9e-4 }, kpi: 'Njs',
    expected: 60 * 6 * Math.pow(muW / rhoW, 0.1) * Math.pow(2e-5, 0.2) * Math.pow(g0 * (5606 - 997) / rhoW, 0.45) * Math.pow(1, 0.13) / Math.pow(0.05, 0.85), tol: 1e-4, kind: 'analytical',
    reference: 'Zwietering correlation N_js = S ν^0.1 d_p^0.2 (gΔρ/ρ_L)^0.45 X^0.13 / D^0.85 (Paul, Atiemo-Obeng & Kresta, Handbook of Industrial Mixing)',
    note: 'An empirical correlation: the benchmark checks its evaluation (exponents and units), not its accuracy for fine particles.',
  },
  {
    id: 's09-grenville-mixing-time', suite: 9, title: 'Turbulent mixing time from the Grenville–Nienow correlation', mode: 'stirred',
    inputs: { mode: 'stirred', imp: 'rushton', Dimp: RUSH.D, Tt: RUSH.T, Nrpm: RUSH.N * 60, rho: 997, mu: 8.9e-4 }, kpi: 'tmix',
    expected: 5.2 * RUSH.T ** 2 / (Math.cbrt(5.0) * RUSH.N * RUSH.D ** 2), tol: 0.005, kind: 'analytical',
    reference: 'Grenville–Nienow turbulent blending correlation N_p^⅓ Re Fo = 5.2 for H = T (Paul, Atiemo-Obeng & Kresta, Handbook of Industrial Mixing)',
    note: '5.5 s for a 0.3 m tank at 300 rpm.',
  },
  {
    id: 's09-kolmogorov-scale', suite: 9, title: 'Kolmogorov length scale (ν³/ε)^¼ from the mean dissipation', mode: 'stirred',
    inputs: { mode: 'stirred', imp: 'rushton', Dimp: RUSH.D, Tt: RUSH.T, Nrpm: RUSH.N * 60, rho: 997, mu: 8.9e-4 }, kpi: 'eta',
    expected: Math.pow((muW / rhoW) ** 3 / (P_RUSH / (rhoW * PI * RUSH.T ** 3 / 4)), 0.25) * 1e6, tol: 0.005, kind: 'analytical',
    reference: 'Kolmogorov microscale η = (ν³/ε)^¼ with ε = P/(ρV) (Pope, Turbulent Flows)',
    note: 'Liquid height equal to tank diameter; about 39 µm.',
  },
  {
    id: 's09-impeller-tip-speed', suite: 9, title: 'Impeller tip speed πND', mode: 'stirred',
    inputs: { mode: 'stirred', Dimp: 0.1, Tt: 0.3, Nrpm: 300 }, kpi: 'tip',
    expected: PI * 5 * 0.1, tol: 1e-6, kind: 'analytical',
    reference: 'Kinematics of a rotating impeller, v_tip = πND',
  },
  // ------------------------------------------------------------------ particle–fluid
  {
    id: 's09-stokes-terminal-velocity', suite: 9, title: 'Stokes terminal velocity of a 1 µm ZnO sphere in water', mode: 'particle',
    inputs: { mode: 'particle', dpS: 1e-6, rhoS: 5606, rho: 997, mu: 8.9e-4, lam: 0, phiS: 0 }, kpi: 'vt',
    expected: stokes(1e-6, 5606, rhoW, muW), tol: 1e-3, kind: 'analytical',
    reference: 'Stokes law v_t = Δρ g d²/(18μ)',
    note: 'Re ≈ 3×10⁻⁶, where the Schiller–Naumann correction used by the suite is below 0.01 %.',
  },
  {
    id: 's09-richardson-zaki-hindered', suite: 9, title: 'Richardson–Zaki hindered settling at 10 vol% (n = 4.65)', mode: 'particle',
    inputs: { mode: 'particle', dpS: 1e-6, rhoS: 5606, rho: 997, mu: 8.9e-4, lam: 0, phiS: 0.1 }, kpi: 'vh',
    expected: stokes(1e-6, 5606, rhoW, muW) * Math.pow(0.9, 4.65), tol: 1e-3, kind: 'analytical',
    reference: 'Richardson–Zaki equation v = v_t(1 − φ)^n with n = 4.65 in the Stokes regime',
    note: 'The hindrance factor is 0.613.',
  },
  {
    id: 's09-aerosol-settling-1um', suite: 9, title: 'Settling velocity of a 1 µm unit-density sphere in air at 20 °C', mode: 'particle',
    inputs: { mode: 'particle', dpS: 1e-6, rhoS: 1000, rho: rhoA, mu: muA, Tp: 20, lam: lamA, phiS: 0 }, kpi: 'vt',
    expected: 3.5e-5, tol: 0.03, kind: 'reference-data',
    reference: 'Hinds, Aerosol Technology — tabulated properties of airborne unit-density spheres at 293 K and 101 kPa (terminal settling velocity of a 1 µm particle about 3.5×10⁻⁵ m/s)',
    note: 'Includes the Cunningham slip correction (about 1.16). Handbook value quoted to two figures; slip-correction constants differ slightly between sources.',
  },
  {
    id: 's09-aerosol-diffusivity-100nm', suite: 9, title: 'Brownian diffusivity of a 0.1 µm particle in air at 20 °C', mode: 'particle',
    inputs: { mode: 'particle', dpS: 1e-7, rhoS: 1000, rho: rhoA, mu: muA, Tp: 20, lam: lamA, phiS: 0 }, kpi: 'D',
    expected: 6.9e-10, tol: 0.03, kind: 'reference-data',
    reference: 'Hinds, Aerosol Technology — tabulated properties of airborne particles at 293 K and 101 kPa (diffusion coefficient of a 0.1 µm particle about 6.9×10⁻¹⁰ m²/s)',
    note: 'Stokes–Einstein with a Cunningham slip correction of about 2.9. Handbook value quoted to two figures.',
  },
  {
    id: 's09-stokes-number', suite: 9, title: 'Stokes number ρ_p d²U/(18μL) of a 10 µm particle', mode: 'particle',
    inputs: { mode: 'particle', dpS: 1e-5, rhoS: 2000, rho: 997, mu: 8.9e-4, Uc: 1, Lc: 1e-4, lam: 0 }, kpi: 'Stk',
    expected: 2000 * 1e-10 * 1 / (18 * 8.9e-4 * 1e-4), tol: 1e-4, kind: 'analytical',
    reference: 'Particle relaxation time τ_p = ρ_p d²/(18μ) and Stokes number Stk = τ_p U/L',
  },
  {
    id: 's09-soret-coefficient', suite: 9, title: 'Soret coefficient D_T/D of a 20 µm particle in water', mode: 'particle',
    inputs: { mode: 'particle', dpS: 2e-5, mu: 8.9e-4, Tp: 25, DT: 1e-11, lam: 0 }, kpi: 'ST',
    expected: 1e-11 / (kB * 298.15 / (3 * PI * 8.9e-4 * 2e-5)), tol: 1e-4, kind: 'analytical',
    reference: 'Definition S_T = D_T/D with the Stokes–Einstein diffusivity D = kT/(3πμd)',
  },
];
