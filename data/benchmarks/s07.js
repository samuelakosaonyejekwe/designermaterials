// Reference benchmarks for Suite 7 — chemical reaction, nucleation and growth.
// Expected values are closed-form results evaluated here from standard constants, exact invariants or
// widely tabulated data; none of them is taken from the suite's own output.
const R = 8.314462618, KB = 1.380649e-23, HP = 6.62607015e-34, FARADAY = 96485.33212, PI = Math.PI;

// --- Arrhenius / Eyring: exact synthetic data sets (the fit must return the parameters used to build them)
const A0 = 1.5e6, EA0 = 65e3;
const temps = [40, 50, 60, 70, 80, 90, 100];
const arrTable = 'T_C,k\n' + temps.map(t => `${t},${(A0 * Math.exp(-EA0 / (R * (t + 273.15)))).toExponential(12)}`).join('\n');
const DH0 = 60e3, DS0 = -50;
const eyrTable = 'T_C,k\n' + temps.map(t => { const T = t + 273.15; return `${t},${((KB * T / HP) * Math.exp(DS0 / R) * Math.exp(-DH0 / (R * T))).toExponential(12)}`; }).join('\n');
// Gas-phase decomposition of N2O5 (first order), the standard textbook Arrhenius data set
const n2o5Table = 'T_C,k\n25,3.46e-5\n35,1.35e-4\n45,4.98e-4\n55,1.50e-3\n65,4.87e-3';

// --- JMAK
const jm = { T: 673.15, n: 2.5, k0: 5e8, Ea: 150e3 };
const kJ = jm.k0 * Math.exp(-jm.Ea / (R * jm.T));
const tJ = 600;                                                     // s (10 min hold)
const kUnder = (() => { const T = 900, Teq = 1200, Bn = 300; return 1e9 * Math.exp(-200e3 / (R * T)) * Math.exp(-Bn * Teq * Teq / (T * (Teq - T) ** 2)); })();

// --- coarsening
const zno = { Vm: 81.379e-3 / 5606 };                               // m³/mol (ZnO 81.38 g/mol, 5606 kg/m³)
const rip = { T: 363.15, g: 0.25, c: 0.002, D: 1e-9, r0: 5e-9, t: 10800 };
const K_LSW = 8 * rip.g * rip.D * rip.c * zno.Vm ** 2 / (9 * R * rip.T);        // m³/s
const sm = { T: 298.15, mu: 0.89e-3, N0: 1e18, W: 1 };
const kSm = 4 * KB * sm.T / (3 * sm.mu * sm.W), tauSm = 1 / (kSm * sm.N0);

// --- electromigration through a stagnant film (free Zn²⁺, z = 2)
const Pe = 2 * FARADAY * 100 * 200e-6 / (R * 298.15);

export default [
  // ------------------------------------------------------------------ Arrhenius fitting
  {
    id: 's07-arrhenius-recover-Ea', suite: 7, title: 'Arrhenius fit recovers the activation energy of exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: arrTable, T: 90 }, kpi: 'Ea', expected: EA0 / 1000, tol: 1e-6, kind: 'analytical',
    reference: 'Arrhenius equation k = A exp(−E_a/RT); noise-free data generated with E_a = 65 kJ/mol, A = 1.5×10⁶'
  },
  {
    id: 's07-arrhenius-recover-lnA', suite: 7, title: 'Arrhenius fit recovers the pre-exponential factor of exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: arrTable, T: 90 }, kpi: 'lnA', expected: Math.log(A0), tol: 1e-6, kind: 'analytical',
    reference: 'Arrhenius equation; intercept of ln k against 1/T equals ln A'
  },
  {
    id: 's07-arrhenius-nonlinear-fit', suite: 7, title: 'Nonlinear (Levenberg–Marquardt) fit recovers E_a of exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: arrTable, T: 90 }, kpi: 'Ea_nl', expected: EA0 / 1000, tol: 1e-4, kind: 'analytical',
    reference: 'Arrhenius equation; least-squares minimum of noise-free data is the generating parameter set'
  },
  {
    id: 's07-arrhenius-prediction', suite: 7, title: 'Predicted rate constant at 75 °C from exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: arrTable, T: 75 }, kpi: 'k_pred', expected: A0 * Math.exp(-EA0 / (R * 348.15)), tol: 1e-6, kind: 'analytical',
    reference: 'Arrhenius equation evaluated at 348.15 K'
  },
  {
    id: 's07-eyring-recover-dH', suite: 7, title: 'Eyring fit recovers the activation enthalpy of exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: eyrTable, T: 70 }, kpi: 'dH', expected: DH0 / 1000, tol: 1e-6, kind: 'analytical',
    reference: 'Eyring equation k = (k_BT/h) exp(ΔS‡/R) exp(−ΔH‡/RT); data generated with ΔH‡ = 60 kJ/mol'
  },
  {
    id: 's07-eyring-recover-dS', suite: 7, title: 'Eyring fit recovers the activation entropy of exact data', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: eyrTable, T: 70 }, kpi: 'dS', expected: DS0, tol: 1e-5, kind: 'analytical',
    reference: 'Eyring equation; data generated with ΔS‡ = −50 J/mol/K (CODATA k_B and h)'
  },
  {
    id: 's07-n2o5-activation-energy', suite: 7, title: 'Activation energy of N₂O₅ decomposition from tabulated rate constants', mode: 'arrfit',
    inputs: { mode: 'arrfit', data: n2o5Table, T: 45 }, kpi: 'Ea', expected: 103, tol: 0.02, kind: 'reference-data',
    reference: 'First-order decomposition of N₂O₅, 25–65 °C: textbook Arrhenius data set with E_a ≈ 103 kJ/mol (Atkins, Physical Chemistry)',
    note: 'Measured rate constants quoted to three figures; the tolerance covers the fitting scatter.'
  },
  // ------------------------------------------------------------------ JMAK
  {
    id: 's07-jmak-rate-constant', suite: 7, title: 'JMAK: Arrhenius rate constant at 400 °C', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius' }, kpi: 'k', expected: kJ, tol: 1e-5, kind: 'analytical',
    reference: 'k = k₀ exp(−E_a/RT)'
  },
  {
    id: 's07-jmak-half-time', suite: 7, title: 'JMAK: half-transformation time (ln 2)^(1/n)/k', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius' }, kpi: 't50', expected: Math.pow(Math.LN2, 1 / jm.n) / kJ, tol: 1e-5, kind: 'analytical',
    reference: 'Johnson–Mehl–Avrami–Kolmogorov equation X = 1 − exp[−(kt)ⁿ] (Christian, The Theory of Transformations in Metals and Alloys)'
  },
  {
    id: 's07-jmak-fraction', suite: 7, title: 'JMAK: transformed fraction after a 10 min hold', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius', tmin: 10 }, kpi: 'X', expected: 1 - Math.exp(-Math.pow(kJ * tJ, jm.n)), tol: 1e-5, kind: 'analytical',
    reference: 'Johnson–Mehl–Avrami–Kolmogorov equation X = 1 − exp[−(kt)ⁿ]'
  },
  {
    id: 's07-jmak-t99', suite: 7, title: 'JMAK: time to 99 % transformation', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius' }, kpi: 't99', expected: Math.pow(Math.log(100), 1 / jm.n) / kJ, tol: 1e-5, kind: 'analytical',
    reference: 'Inversion of the JMAK equation, t = [−ln(1 − X)]^(1/n)/k'
  },
  {
    id: 's07-jmak-peak-rate-time', suite: 7, title: 'JMAK: time of maximum transformation rate', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius' }, kpi: 'tpeak', expected: Math.pow((jm.n - 1) / jm.n, 1 / jm.n) / kJ, tol: 1e-5, kind: 'analytical',
    reference: 'd²X/dt² = 0 for the JMAK equation gives (kt)ⁿ = (n − 1)/n'
  },
  {
    id: 's07-jmak-temperature-for-half', suite: 7, title: 'JMAK: temperature giving 50 % transformation in 60 min', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 400, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius', tmin: 60 }, kpi: 'Tnose',
    expected: jm.Ea / (R * Math.log(jm.k0 * 3600 / Math.pow(Math.LN2, 1 / jm.n))) - 273.15, tol: 1e-5, kind: 'analytical',
    reference: 'Arrhenius inversion T = E_a/[R ln(k₀ t/(ln 2)^(1/n))]'
  },
  {
    id: 's07-jmak-undercooling-rate', suite: 7, title: 'JMAK: rate constant with a nucleation-barrier (undercooling) factor', mode: 'jmak',
    inputs: { mode: 'jmak', Tj: 900 - 273.15, k0: 1e9, EaJ: 200, drive: 'undercool', Teq: 1200, Bn: 300 }, kpi: 'k', expected: kUnder, tol: 1e-5, kind: 'analytical',
    reference: 'Classical C-curve kinetics: k = k₀ exp(−E_a/RT) exp[−B T_eq²/(T ΔT²)] (Porter & Easterling, Phase Transformations in Metals and Alloys)'
  },
  // ------------------------------------------------------------------ ripening and aggregation
  {
    id: 's07-lsw-rate-constant', suite: 7, title: 'LSW ripening rate constant of ZnO in water at 90 °C', mode: 'ripen',
    inputs: { mode: 'ripen', chem: 'zn_nitrate', T: 90, gam: 0.25, cInf: 0.002, Dm: 1e-9 }, kpi: 'K_lsw', expected: K_LSW * 1e27, tol: 1e-3, kind: 'analytical',
    reference: 'Lifshitz–Slyozov–Wagner theory, K = 8γ D c∞ V_m²/(9RT)',
    note: 'Molar volume from M = 81.38 g/mol and ρ = 5606 kg/m³.'
  },
  {
    id: 's07-lsw-cube-law', suite: 7, title: 'LSW cube law: mean radius after 3 h', mode: 'ripen',
    inputs: { mode: 'ripen', chem: 'zn_nitrate', T: 90, gam: 0.25, cInf: 0.002, Dm: 1e-9, r0: 5, tmin: 180 }, kpi: 'r_end',
    expected: Math.cbrt(rip.r0 ** 3 + K_LSW * rip.t) * 1e9, tol: 1e-3, kind: 'analytical',
    reference: 'Lifshitz–Slyozov–Wagner theory, r̄³ − r̄₀³ = K t'
  },
  {
    id: 's07-lsw-doubling-time', suite: 7, title: 'LSW: time to double the mean radius, 7r₀³/K', mode: 'ripen',
    inputs: { mode: 'ripen', chem: 'zn_nitrate', T: 90, gam: 0.25, cInf: 0.002, Dm: 1e-9, r0: 5 }, kpi: 't2_rip', expected: 7 * rip.r0 ** 3 / K_LSW, tol: 1e-3, kind: 'analytical',
    reference: 'Lifshitz–Slyozov–Wagner cube law with r̄ = 2r̄₀'
  },
  {
    id: 's07-smoluchowski-rate-constant', suite: 7, title: 'Smoluchowski rapid-coagulation rate constant in water at 25 °C', mode: 'ripen',
    inputs: { mode: 'ripen', T: 25, mu: 0.89e-3, W: 1 }, kpi: 'k_smol', expected: kSm, tol: 1e-4, kind: 'analytical',
    reference: 'Smoluchowski theory of rapid Brownian coagulation, dN/dt = −kN² with k = 4k_BT/(3μ) ≈ 6×10⁻¹⁸ m³/s in water (Hiemenz & Rajagopalan, Principles of Colloid and Surface Chemistry)'
  },
  {
    id: 's07-smoluchowski-half-time', suite: 7, title: 'Smoluchowski coagulation half-time 3μ/(4k_BT N₀)', mode: 'ripen',
    inputs: { mode: 'ripen', T: 25, mu: 0.89e-3, W: 1, N0: 1e18 }, kpi: 't_half', expected: 3 * sm.mu / (4 * KB * sm.T * sm.N0), tol: 1e-4, kind: 'analytical',
    reference: 'Smoluchowski theory of rapid coagulation: t½ = 3μ/(4k_BT N₀), about 2×10¹¹/N₀ s with N₀ per cm³ for water'
  },
  {
    id: 's07-smoluchowski-number-decay', suite: 7, title: 'Smoluchowski: particle number after ten half-times', mode: 'ripen',
    inputs: { mode: 'ripen', T: 25, mu: 0.89e-3, W: 1, N0: 1e16, tmin: 10 * (100 * tauSm) / 60 }, kpi: 'N_rel', expected: 1 / 11, tol: 1e-4, kind: 'analytical',
    reference: 'Smoluchowski constant-kernel solution N/N₀ = 1/(1 + t/t½)',
    note: 'N₀ = 10¹⁶ m⁻³, for which the half-time is 16.2 s and the hold of ten half-times lies inside the input range.'
  },
  {
    id: 's07-fractal-aggregate-radius', suite: 7, title: 'Fractal aggregate radius r₀ n^(1/d_f)', mode: 'ripen',
    inputs: { mode: 'ripen', T: 25, mu: 0.89e-3, W: 1, N0: 1e18, tmin: 99 * tauSm / 60, r0: 10, df: 2 }, kpi: 'r_agg', expected: 100, tol: 1e-3, kind: 'analytical',
    reference: 'Fractal scaling of aggregates, R = r₀ n^(1/d_f): n = 100 primary particles with d_f = 2 gives R = 10 r₀'
  },
  // ------------------------------------------------------------------ LaMer burst nucleation
  {
    id: 's07-lamer-yield-zno', suite: 7, title: 'LaMer batch: stoichiometric ZnO yield at complete conversion', mode: 'lamer',
    inputs: { mode: 'lamer', chem: 'zn_nitrate', cP: 10, Vb: 5e-5, tmin: 180 }, kpi: 'yield', expected: 10 * 5e-5 * 81.379 * 1000, tol: 0.002, kind: 'conservation',
    reference: 'Mass balance: 0.5 mmol of Zn gives 0.5 mmol ZnO (M = 81.38 g/mol) = 40.69 mg',
    note: 'The residual dissolved monomer (solubility 0.002 mol/m³ plus Gibbs–Thomson excess) accounts for the small shortfall.'
  },
  {
    id: 's07-lamer-yield-gold', suite: 7, title: 'LaMer batch: stoichiometric gold yield (Turkevich-type recipe)', mode: 'lamer',
    inputs: { mode: 'lamer', chem: 'au', T: 100, tmin: 15, cP: 0.25, kf: 0.02, gam: 0.45, cInf: 1e-6, Dm: 1.5e-9, ks: 0.02, logA: 34, Vb: 1e-4 }, kpi: 'yield',
    expected: 0.25 * 1e-4 * 196.967 * 1000, tol: 0.002, kind: 'conservation',
    reference: 'Mass balance: 25 µmol of Au(III) reduced to 25 µmol Au (M = 196.97 g/mol) = 4.924 mg'
  },
  // ------------------------------------------------------------------ ZnO nanorod bath
  {
    id: 's07-zno-hmta-first-order', suite: 7, title: 'HMTA hydrolysis: first-order conversion in an isothermal bath', mode: 'zno',
    inputs: { mode: 'zno', T: 90, T0: 90, tmin: 180, kh: 3e-5, cHMTA: 25 }, kpi: 'hmta', expected: 1 - Math.exp(-3e-5 * 10800), tol: 0.005, kind: 'analytical',
    reference: 'First-order kinetics, conversion = 1 − exp(−kt)',
    note: 'The suite integrates with backward Euler; the tolerance covers its first-order time-discretisation error.'
  },
  {
    id: 's07-zno-salt-solution-ph', suite: 7, title: 'pH of a 0.025 mol/L zinc salt solution from the first hydrolysis step', mode: 'zno',
    inputs: { mode: 'zno', T: 25, T0: 25, cZn: 25, cHMTA: 0, cOH: 0, tmin: 10 }, kpi: 'pH_end', expected: 0.5 * (8.96 - Math.log10(0.025)), abs: 0.01, kind: 'analytical',
    reference: 'Weak-acid approximation [H⁺] = √(*β₁ C) with log *β₁ = −8.96 for Zn²⁺ (Baes & Mesmer, The Hydrolysis of Cations)'
  },
  {
    id: 's07-zno-no-base-no-growth', suite: 7, title: 'No hydroxide source: bath stays undersaturated and no zinc is consumed', mode: 'zno',
    inputs: { mode: 'zno', T: 90, cZn: 25, cHMTA: 0, cOH: 0, tmin: 180 }, kpi: 'conv', expected: 0, abs: 1e-6, kind: 'conservation',
    reference: 'Thermodynamic constraint: without supersaturation with respect to ZnO there is no growth or precipitation'
  },
  {
    id: 's07-zno-electromigration-factor', suite: 7, title: 'Nernst–Planck enhancement of film mass transfer for Zn²⁺', mode: 'zno',
    inputs: { mode: 'zno', T: 25, T0: 25, cZn: 25, cHMTA: 0, cOH: 0, tmin: 10, E: 100, delta: 200 }, kpi: 'enh', expected: Pe / (1 - Math.exp(-Pe)), tol: 0.002, kind: 'analytical',
    reference: 'Steady 1-D Nernst–Planck flux through a stagnant film in a uniform field: J/J_diffusion = Pe/(1 − e^(−Pe)), Pe = zFEδ/RT (Newman, Electrochemical Systems)',
    note: 'Acidic zinc salt solution, where zinc is essentially all free Zn²⁺ (z = 2); the suite uses the species-averaged charge.'
  },
  {
    id: 's07-zno-nanorod-length-3h', suite: 7, title: 'ZnO nanorods from 25 mmol/L zinc nitrate + HMTA, 90 °C, 3 h: rod length', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 25, cHMTA: 25, cOH: 0, T: 90, tmin: 180 }, kpi: 'L_final', expected: 1500, tol: 0.35, kind: 'reference-data',
    reference: 'Seeded hydrothermal growth of ZnO nanorod arrays (Vayssieres; Greene et al.): lengths of roughly 1–2 µm after about 3 h at 90 °C in equimolar 25 mmol/L baths',
    note: 'Order-of-magnitude check only: reported lengths scatter between about 1 and 2 µm with seed layer and vessel, and the surface rate constant of this reduced model is a calibration parameter. The tolerance spans that literature range.'
  },
];
