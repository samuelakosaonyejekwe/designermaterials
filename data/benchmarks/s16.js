// Suite 16 — reference benchmarks for microstructure, phase-field and defect evolution.
// Expected values are literals or closed forms evaluated here with CODATA constants; none of them calls the suite.
const KB = 1.380649e-23, QE = 1.602176634e-19, NA = 6.02214076e23, R = 8.314462618, KB_EV = KB / QE;
const K = c => c + 273.15;

// --- Cahn–Hilliard, quartic double well f = f0·¼(φ²−1)², gradient coefficient κ (in f0·Δx²) ---
const chK = 2, chF0 = 2e8, chDx = 2e-9, chC0 = 0.4, chPhi0 = 2 * chC0 - 1;
const sigmaCH = (2 * Math.SQRT2 / 3) * Math.sqrt(chK) * chF0 * chDx;                    // σ = ∫√(2κf) dφ
const lambdaCH = 2 * Math.PI * Math.sqrt(2 * chK / (1 - 3 * chPhi0 * chPhi0)) * chDx;   // k_max² = −f″/(2κ)

// --- analytical kinetics ---
const kin = { mode: 'kinetics', T_C: 700, time_h: 10, d0_nm: 20, M0: 1e-5, Qgb: 200, gamma_gb: 0.5, D0: 1e-5, Qd: 180, c_sol: 100, Vm: 14.5, gamma_int: 0.3, r0_nm: 5, dp_nm: 2000, rho_green: 0.55, n_av: 2.5, k0_av: 1e8, Q_av: 180 };
const kT = K(kin.T_C), kt = kin.time_h * 3600, RT = R * kT;
const Mgb = kin.M0 * Math.exp(-kin.Qgb * 1e3 / RT);
const dParabolic = Math.sqrt((kin.d0_nm * 1e-9) ** 2 + 4 * Mgb * kin.gamma_gb * kt);
const Dk = kin.D0 * Math.exp(-kin.Qd * 1e3 / RT), VmSI = kin.Vm * 1e-6;
const rLSW = Math.cbrt((kin.r0_nm * 1e-9) ** 3 + 8 * kin.gamma_int * Dk * kin.c_sol * VmSI * VmSI / (9 * RT) * kt);
const sin = { T_C: 600, time_h: 1, dp_nm: 2000 }, sT = K(sin.T_C), Ds = kin.D0 * Math.exp(-kin.Qd * 1e3 / (R * sT));
const neckVol = Math.pow(80 * Ds * kin.gamma_gb * (VmSI / NA) / (KB * sT) * sin.time_h * 3600 / (sin.dp_nm * 1e-9) ** 3, 1 / 5);
const rhoVol = kin.rho_green / (1 - (neckVol / 2) ** 2) ** 3;
const kAv = kin.k0_av * Math.exp(-kin.Q_av * 1e3 / RT);
const jm = { T_C: 500, time_h: 4 }, kAv500 = kin.k0_av * Math.exp(-kin.Q_av * 1e3 / (R * K(jm.T_C)));
const xJmak = 1 - Math.exp(-Math.pow(kAv500 * jm.time_h * 3600, kin.n_av));

// --- defect equilibria ---
const cu = { T_C: 1000, Evf: 1.28, Svf: 0 }, cVac = Math.exp(-cu.Evf / (KB_EV * K(cu.T_C)));
const red = { T_C: 900, dH_red: 4, dS_red: 10, pO2: 1e-5 }, Kred = Math.exp(red.dS_red) * Math.exp(-red.dH_red / (KB_EV * K(red.T_C)));
const vIntr = Math.cbrt(Kred / (4 * Math.sqrt(red.pO2)));          // n = 2[V_O], [V_O]·n²·pO2^½ = K  (per formula unit)
const nZnO = 5606 * NA / ((65.38 + 15.999) * 1e-3);                // formula units per m³ from handbook density and molar mass
const sch = { T_C: 1500, dH_S: 5 }, xSch = Math.pow(Math.exp(-sch.dH_S / (KB_EV * K(sch.T_C))) / (2 ** 2 * 3 ** 3), 1 / 5);
const tay = { G_shear: 4.8e10, b_nm: 0.256, rho_d: 1e14, alpha_T: 0.3, sigma0: 0 }, sigTaylor = 3.06 * tay.alpha_T * tay.G_shear * tay.b_nm * 1e-9 * Math.sqrt(tay.rho_d);
const seg = { T_C: 600, x_bulk: 0.002, dG_seg: 50 }, qSeg = seg.x_bulk / (1 - seg.x_bulk) * Math.exp(seg.dG_seg * 1e3 / (R * K(seg.T_C)));

export default [
  // ---------------- spinodal decomposition ----------------
  {
    id: 's16-ch-mass-conservation', suite: 16, title: 'Cahn–Hilliard: conservation of mean composition', mode: 'spinodal',
    inputs: { mode: 'spinodal', c0: 0.5, ch_steps: 200 }, kpi: 'mass_err', expected: 0, abs: 1e-10, kind: 'conservation',
    reference: 'Cahn & Hilliard, free energy of a non-uniform system — conserved (model B) dynamics ∂φ/∂t = ∇·(M∇μ)',
    note: 'The k = 0 Fourier mode is untouched by the spectral update, so the drift must stay at round-off level.',
  },
  {
    id: 's16-ch-interface-energy', suite: 16, title: 'Cahn–Hilliard: flat-interface energy of the quartic double well', mode: 'spinodal',
    inputs: { mode: 'spinodal', kappa: chK, f0: chF0, dx_nm: chDx * 1e9, ch_steps: 40 }, kpi: 'sigma_int', expected: sigmaCH, tol: 1e-6, kind: 'analytical',
    reference: 'Cahn & Hilliard — σ = ∫√(2κ·Δf) dφ = (2√2/3)·√(κ f₀) for f = f₀·¼(φ²−1)²',
  },
  {
    id: 's16-ch-fastest-wavelength', suite: 16, title: 'Spinodal decomposition: fastest-growing wavelength (Cahn linear theory)', mode: 'spinodal',
    inputs: { mode: 'spinodal', c0: chC0, kappa: chK, dx_nm: chDx * 1e9, ch_steps: 40 }, kpi: 'lambda_max', expected: lambdaCH, tol: 1e-6, kind: 'analytical',
    reference: 'Cahn, On spinodal decomposition — λ_max = 2π·√(−2κ/f″), f″ = f₀(3φ₀² − 1)',
  },
  {
    id: 's16-ch-symmetric-phase-fraction', suite: 16, title: 'Spinodal decomposition of a 50:50 mixture: equal phase fractions', mode: 'spinodal',
    inputs: { mode: 'spinodal', c0: 0.5 }, kpi: 'frac_B', expected: 0.5, tol: 0.03, kind: 'analytical',
    reference: 'Lever rule with the symmetric binodal φ = ±1 of the quartic double well',
    note: 'Measured as the area with c > 0.5 on a 64×64 grid after 120 steps; the residual is the sampling noise of one finite periodic box.',
  },
  // ---------------- Potts grain growth ----------------
  {
    id: 's16-potts-parabolic-growth', suite: 16, title: 'Normal grain growth without particles: parabolic law d² − d₀² ∝ t (exponent ½)', mode: 'grain',
    inputs: { mode: 'grain', fp: 0, mcs: 200 }, kpi: 'm_parabolic', expected: 0.5, tol: 0.05, kind: 'analytical',
    reference: 'Burke & Turnbull parabolic grain-growth law for curvature-driven boundary migration (d² − d₀² = K t)',
    note: 'Offset-corrected exponent m in ⟨A⟩ − A_ref ∝ (t − t_ref)^(2m), fitted to the mean grain area pooled over twelve seeded 64×64 zero-temperature Potts lattices from the first sample with at least 25 cells per grain. The reference law itself contains the initial size d₀, so the comparison is like for like. Jackknife standard error about 0.014 (KPI m_se); other seeds give 0.47–0.50. The uncorrected late-time slope d ln d/d ln t (KPI late_slope) reads 0.43–0.48 at 200 steps because d₀ is not negligible, and is not benchmarked.',
  },
  {
    id: 's16-potts-zener-limit', suite: 16, title: 'Zener–Smith limiting grain size for single-cell particles', mode: 'grain',
    inputs: { mode: 'grain', fp: 0.02, dx_nm: 5, mcs: 40 }, kpi: 'zener_3d', expected: 4 * (5e-9 / Math.sqrt(Math.PI)) / (3 * 0.02), tol: 1e-6, kind: 'analytical',
    reference: 'Zener (communicated by C. S. Smith) — D_lim = 4r/(3f)',
    note: 'One grid cell is treated as a particle of equal area, r = Δx/√π.',
  },
  // ---------------- dendritic solidification ----------------
  {
    id: 's16-kobayashi-enthalpy', suite: 16, title: 'Phase-field solidification: enthalpy conservation with zero-flux walls', mode: 'solid',
    inputs: { mode: 'solid' }, kpi: 'heat_err', expected: 0, abs: 1e-9, kind: 'conservation',
    reference: 'Kobayashi, Modeling and numerical simulations of dendritic crystal growth — ∂T/∂t = ∇²T + K ∂φ/∂t conserves ∫(T − Kφ) dA',
  },
  // ---------------- analytical kinetics ----------------
  {
    id: 's16-parabolic-grain-growth', suite: 16, title: 'Unpinned grain growth at 700 °C: d² − d₀² = 4Mγt', mode: 'kinetics',
    inputs: { ...kin, fp: 0 }, kpi: 'd_final', expected: dParabolic, tol: 0.002, kind: 'analytical',
    reference: 'Burke & Turnbull parabolic law with Arrhenius mobility, dd/dt = 2Mγ/d',
  },
  {
    id: 's16-zener-stagnation', suite: 16, title: 'Pinned grain growth saturates at the Zener limit', mode: 'kinetics',
    inputs: { ...kin, T_C: 1000, time_h: 1000, fp: 0.05, rp_nm: 10 }, kpi: 'd_final', expected: 4 * 10e-9 / (3 * 0.05), tol: 0.01, kind: 'analytical',
    reference: 'Zener–Smith pinning — stagnation of dd/dt = 2Mγ(1/d − 1/d_Z) at d_Z = 4r/(3f)',
    note: 'Long hold at high temperature so that the asymptote is reached.',
  },
  {
    id: 's16-lsw-cube-law', suite: 16, title: 'Ostwald ripening: LSW cube law', mode: 'kinetics',
    inputs: { ...kin }, kpi: 'r_lsw', expected: rLSW, tol: 1e-6, kind: 'analytical',
    reference: 'Lifshitz & Slyozov; Wagner — r³ − r₀³ = 8γDc∞Vₘ²t/(9RT)',
  },
  {
    id: 's16-sinter-neck-volume-diffusion', suite: 16, title: 'Initial-stage sintering neck ratio, volume diffusion', mode: 'kinetics',
    inputs: { ...kin, ...sin, sinter: 'vol' }, kpi: 'neck', expected: neckVol, tol: 1e-6, kind: 'analytical',
    reference: 'German, Sintering Theory and Practice — (X/D)⁵ = 80 D_v γ Ω t/(k_B T D³)',
  },
  {
    id: 's16-sinter-density', suite: 16, title: 'Initial-stage shrinkage and density, volume diffusion', mode: 'kinetics',
    inputs: { ...kin, ...sin, sinter: 'vol' }, kpi: 'rho_rel', expected: rhoVol, tol: 1e-6, kind: 'analytical',
    reference: 'German, Sintering Theory and Practice — ΔL/L₀ = (X/2D)², ρ = ρ₀/(1 − ΔL/L₀)³',
  },
  {
    id: 's16-jmak-fraction', suite: 16, title: 'JMAK transformed fraction after 4 h at 500 °C', mode: 'kinetics',
    inputs: { ...kin, ...jm }, kpi: 'X_jmak', expected: xJmak, tol: 1e-6, kind: 'analytical',
    reference: 'Johnson–Mehl–Avrami–Kolmogorov — X = 1 − exp[−(kt)ⁿ], k = k₀·exp(−Q/RT)',
  },
  {
    id: 's16-jmak-half-time', suite: 16, title: 'JMAK half-transformation time', mode: 'kinetics',
    inputs: { ...kin }, kpi: 't50', expected: Math.pow(Math.LN2, 1 / kin.n_av) / kAv, tol: 1e-6, kind: 'analytical',
    reference: 'Johnson–Mehl–Avrami–Kolmogorov — t₅₀ = (ln 2)^(1/n)/k',
  },
  // ---------------- defect equilibria ----------------
  {
    id: 's16-vacancy-fraction', suite: 16, title: 'Equilibrium vacancy fraction exp(−E_f/k_BT), E_f = 1.28 eV at 1000 °C', mode: 'defects',
    inputs: { mode: 'defects', ...cu }, kpi: 'c_vac', expected: cVac, tol: 1e-6, kind: 'analytical',
    reference: 'Equilibrium thermodynamics of point defects (e.g. Porter & Easterling, Phase Transformations in Metals and Alloys)',
    note: 'Formation entropy set to zero so that the bare Boltzmann factor is tested (copper-like formation energy).',
  },
  {
    id: 's16-oxygen-vacancy-intrinsic', suite: 16, title: 'ZnO reduction equilibrium: [V_O••] = (K_red/4)^⅓·pO₂^(−1/6)', mode: 'defects',
    inputs: { mode: 'defects', oxide: 'ZnO', ...red, x_donor: 0 }, kpi: 'VO', expected: vIntr, tol: 1e-5, kind: 'analytical',
    reference: 'Kröger–Vink mass action with electroneutrality n = 2[V_O••] (Kröger, The Chemistry of Imperfect Crystals)',
  },
  {
    id: 's16-brouwer-intrinsic-slope', suite: 16, title: 'Brouwer slope in the intrinsic reduction regime', mode: 'defects',
    inputs: { mode: 'defects', oxide: 'TiO2', ...red, x_donor: 0 }, kpi: 'brouwer', expected: -1 / 6, abs: 1e-4, kind: 'analytical',
    reference: 'Brouwer approximation — [V_O••] ∝ pO₂^(−1/6) when n = 2[V_O••]',
  },
  {
    id: 's16-brouwer-donor-slope', suite: 16, title: 'Brouwer slope in the donor-controlled regime', mode: 'defects',
    inputs: { mode: 'defects', oxide: 'ZnO', T_C: 600, pO2: 1, x_donor: 0.01 }, kpi: 'brouwer', expected: -0.5, abs: 1e-4, kind: 'analytical',
    reference: 'Brouwer approximation — [V_O••] ∝ pO₂^(−1/2) when donors fix n = [D•]',
  },
  {
    id: 's16-donor-electron-density', suite: 16, title: 'ZnO with 1 % fully ionised donors: electron density', mode: 'defects',
    inputs: { mode: 'defects', oxide: 'ZnO', T_C: 600, pO2: 1, x_donor: 0.01 }, kpi: 'n_e', expected: 0.01 * nZnO, tol: 0.005, kind: 'reference-data',
    reference: 'Cation site density of wurtzite ZnO from handbook density 5.606 g/cm³ and molar mass 81.38 g/mol (CRC Handbook of Chemistry and Physics)',
    note: 'n = [D•] in the extrinsic regime; the vacancy contribution is negligible at 600 °C in 1 atm oxygen.',
  },
  {
    id: 's16-schottky-alumina', suite: 16, title: 'Schottky disorder in Al₂O₃: (2x)²(3x)³ = exp(−ΔH_S/k_BT)', mode: 'defects',
    inputs: { mode: 'defects', oxide: 'Al2O3', ...sch }, kpi: 'schottky', expected: xSch, tol: 1e-6, kind: 'analytical',
    reference: 'Mass action for nil ⇌ 2 V_Al‴ + 3 V_O•• (Kingery, Bowen & Uhlmann, Introduction to Ceramics)',
  },
  {
    id: 's16-taylor-hardening', suite: 16, title: 'Taylor flow stress of copper at ρ = 10¹⁴ m⁻²', mode: 'defects',
    inputs: { mode: 'defects', ...tay }, kpi: 'sigma_T', expected: sigTaylor, tol: 1e-6, kind: 'analytical',
    reference: 'Taylor relation σ = MαGb√ρ with Taylor factor M = 3.06 for untextured fcc polycrystals',
  },
  {
    id: 's16-kocks-mecking-saturation', suite: 16, title: 'Kocks–Mecking saturation dislocation density', mode: 'defects',
    inputs: { mode: 'defects', k1: 3e8, k2: 6 }, kpi: 'rho_sat', expected: (3e8 / 6) ** 2, tol: 1e-9, kind: 'analytical',
    reference: 'Kocks & Mecking — steady state of dρ/dε = M(k₁√ρ − k₂ρ)',
  },
  {
    id: 's16-mclean-segregation', suite: 16, title: 'Langmuir–McLean grain-boundary segregation isotherm', mode: 'defects',
    inputs: { mode: 'defects', ...seg }, kpi: 'x_gb', expected: qSeg / (1 + qSeg), tol: 1e-6, kind: 'analytical',
    reference: 'McLean, Grain Boundaries in Metals — X_gb/(1 − X_gb) = X_b/(1 − X_b)·exp(ΔG/RT)',
  },
  // ------------------------------------------------------------------ published experiments
  {
    id: 's16-exp-simmons-balluffi-al-vacancies-melting', suite: 16, title: 'Aluminium at the melting point: measured equilibrium vacancy concentration (differential dilatometry)', mode: 'defects',
    inputs: { mode: 'defects', T_C: 660, Evf: 0.76, Svf: 2.4 }, kpi: 'c_vac', expected: 9.4e-4, tol: 0.15, kind: 'experiment',
    measured: { value: 9.4e-4, uncertainty: null, unit: 'site fraction' },
    source: { authors: 'Simmons, R. O.; Balluffi, R. W.', title: 'Measurements of Equilibrium Vacancy Concentrations in Aluminum', journal: 'Phys. Rev. 117, 52', year: 1960, doi: '10.1103/PhysRev.117.52', url: 'https://doi.org/10.1103/PhysRev.117.52' },
    conditions: '99.995 % aluminium bar; macroscopic length change ΔL/L (filar micrometer microscopes) and X-ray lattice-parameter change Δa/a measured at identical temperatures from 229 to 656 °C; added site fraction 3(ΔL/L − Δa/a) at the melting point (660 °C)',
    match: ['T_C', 'Evf', 'Svf'],
    reference: 'Simmons & Balluffi, Phys. Rev. 117, 52 (1960)',
    note: 'Abstract of the source: the added site fraction "could be represented by 3(ΔL/L − Δa/a) = exp(2.4)·exp(−0.76 eV/kT)" and "at the melting point (660 °C) the equilibrium concentration of additional sites is 9.4×10⁻⁴"; the authors conclude that the defects are predominantly single vacancies (divacancy share at melting estimated below about 15 %). The formation energy and entropy entered are the rounded values of that representation, so this checks the Arrhenius form against the stated end-point value; the suite gives 8.66×10⁻⁴ (−7.8 %). Rounding the formation energy to 0.01 eV alone changes the result by 13 % at this temperature, hence the 15 % tolerance. The source states no uncertainty for the value in the abstract. The numbers were read from the abstract in the open bibliographic record of the article (OpenAlex); the full text was not accessible.',
  },
  {
    id: 's16-exp-simmons-balluffi-au-vacancies-melting', suite: 16, title: 'Gold just below the melting temperature: measured equilibrium vacancy concentration', mode: 'defects',
    inputs: { mode: 'defects', T_C: 1063, Evf: 0.94, Svf: 1.0 }, kpi: 'c_vac', expected: 7.2e-4, tol: 0.15, kind: 'experiment',
    measured: { value: 7.2e-4, uncertainty: 0.6e-4, unit: 'site fraction' },
    source: { authors: 'Simmons, R. O.; Balluffi, R. W.', title: 'Measurement of Equilibrium Concentrations of Lattice Vacancies in Gold', journal: 'Phys. Rev. 125, 862', year: 1962, doi: '10.1103/PhysRev.125.862', url: 'https://doi.org/10.1103/PhysRev.125.862' },
    conditions: '99.999 % gold bar; length expansion and X-ray lattice-parameter expansion measured from 15 to 1057 °C; net added site fraction 3(ΔL/L − Δa/a) just below the melting temperature, obtained by an extrapolation of 6 °C',
    match: ['T_C', 'Evf', 'Svf'],
    reference: 'Simmons & Balluffi, Phys. Rev. 125, 862 (1962)',
    note: 'Abstract of the source: ΔN/N "just below the melting temperature, is (7.2 ± 0.6)×10⁻⁴ as obtained by an extrapolation of the data of only 6 °C" and "can be described by exp(1.0)·exp(−0.94 eV/kT)". The temperature entered, 1063 °C, is the upper end of the measured range (1057 °C) plus the stated 6 °C extrapolation. The suite gives 7.74×10⁻⁴ (+7.5 %), inside the stated experimental uncertainty of ±8 %; the tolerance of 15 % adds the rounding of the published energy (0.01 eV corresponds to 9 % at this temperature). The added sites include a minority of divacancies (the source: more than 80 % single vacancies for plausible binding energies), which the single-vacancy formula does not separate. The numbers were read from the abstract in the open bibliographic record of the article (OpenAlex); the full text was not accessible.',
  },
];
