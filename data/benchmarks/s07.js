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
    id: 's07-zno-nanorod-length-3h', suite: 7, title: 'ZnO nanorods from 25 mmol/L zinc nitrate + HMTA, 90 °C, 3 h: rod length (calibration series)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 25, cHMTA: 25, cOH: 0, T: 90, tmin: 180, Vb: 1e-4 }, kpi: 'L_final', expected: 800, tol: 0.15, kind: 'calibration',
    measured: { value: 800, uncertainty: 59, unit: 'nm' },
    source: { authors: 'Mosalagae, K.; Murape, D. M.; Lepodise, L. M.', title: 'Effects of growth conditions on properties of CBD synthesized ZnO nanorods grown on ultrasonic spray pyrolysis deposited ZnO seed layers', journal: 'Heliyon 6 (7), e04458', year: 2020, doi: '10.1016/j.heliyon.2020.e04458', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7369615/' },
    conditions: 'Glass substrate with a ZnO seed layer from ultrasonic spray pyrolysis, suspended vertically in 100 mL of 25 mM Zn(NO₃)₂·6H₂O + 25 mM HMTA at 90 °C for 3 h',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Mosalagae, Murape & Lepodise, Heliyon 6 (2020) e04458 — Tables 1–3 (rod length against growth time, temperature and concentration)',
    note: 'Calibration case, not validation. Growth law: rods and bulk ZnO particles share the driving force V_m k (c − c_eq), where c_eq is the ZnO solubility at the pH set by the HMTA/ammonia buffer; the number of bulk particles is proportional to the zinc charged, so the rod length follows the cube root of the precipitated fraction rather than the concentration. Calibration set: the ten lengths of Tables 1–3 of this source. Fitted by least absolute deviation of ln(length): axial rate constant 1.3×10⁻⁶ m/s at 90 °C and 4.4×10¹² bulk particles per mole of zinc; the HMTA hydrolysis constant 2×10⁻⁵ 1/s was scanned (0.6–4.5×10⁻⁵, particle yield refitted) and is only weakly determined. Residuals (model − measured): Table 1, 25 mM at 90 °C after 1.5, 2.0, 2.5, 3.0 h (400 ± 27, 690 ± 45, 710 ± 47, 800 ± 59 nm): +62, −1, 0, −8 %; Table 2, 25 mM for 2 h at 60, 70, 80 °C (400, 657, 671 nm): 0, −22, −12 %; Table 3, 90 °C for 2 h at 10, 15, 30 mM (519, 521, 425 nm): +25, +28, +62 %. Mean absolute deviation of ln(length) 0.19. Two calibration points are not reproduced: the 1.5 h length (the model has the precipitation burst complete by then; the heating history of the bath is not reported) and the 30 mM length, which the source itself reports as shorter than at 10–25 mM. The tolerance of this entry covers the residual of the 3 h point only. Substrate area, seed density and heat-up time are not given in the paper and stay at the suite defaults.'
  },
  {
    id: 's07-exp-zno-nanowire-length-70c-12h', suite: 7, title: 'ZnO nanowires from 25 mmol/L zinc nitrate + HMT at 70 °C after 12 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 25, cHMTA: 25, cOH: 0, T: 70, tmin: 720, Vb: 2e-5 }, kpi: 'L_final', expected: 900, tol: 0.3, kind: 'experiment',
    measured: { value: 900, uncertainty: null, unit: 'nm' },
    source: { authors: 'Mason, A. D.; Roberts, T. F.; Conley Jr., J. F.; Price, D. T.; Allman, D. D. J.; McGuire, M. S.', title: 'Investigation of Growth Parameter Influence on Hydrothermally Grown ZnO Nanowires Using a Research Grade Microwave', journal: 'Author manuscript hosted by ON Semiconductor (a two-page contribution with the same title and authors is in the 2009 International Semiconductor Device Research Symposium, doi 10.1109/ISDRS.2009.5378342)', year: 2009, doi: null, url: 'https://www.onsemi.com/site/pdf/SSE_MicrowaveZnO.pdf' },
    conditions: '90 nm ALD ZnO seed layer on Si, 20 mL of 25 mM Zn(NO₃)₂ + 25 mM HMT (1:1) in a sealed 35 mL vial, 70 °C held by a temperature-controlled microwave reactor, slow stirring, 12 h',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Mason, Roberts, Conley, Price, Allman & McGuire, hydrothermal ZnO nanowires grown in a research-grade microwave reactor (author manuscript, ON Semiconductor / Oregon State University)',
    note: 'Value from the text of section 3.1.1 ("At an extremely long growth time of 12 hours, the 70°C / 25 mM NWs grew to 900 nm"; the same section states that the maximum length is limited to about 1 µm). Model 669 nm (−26 %). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). The model is low because at 70 °C its HMTA hydrolysis (activation energy 90 kJ/mol, assumed, not fitted) has released only 14 % of the base after 12 h. Bath volume 20 mL from the experimental section. Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series.'
  },
  {
    id: 's07-exp-zno-nanowire-length-85c-1h', suite: 7, title: 'ZnO nanowires from 25 mmol/L zinc nitrate + HMT at 85 °C after 1 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 25, cHMTA: 25, cOH: 0, T: 85, tmin: 60, Vb: 2e-5 }, kpi: 'L_final', expected: 696, tol: 0.3, kind: 'experiment',
    measured: { value: 696, uncertainty: null, unit: 'nm' },
    source: { authors: 'Mason, A. D.; Roberts, T. F.; Conley Jr., J. F.; Price, D. T.; Allman, D. D. J.; McGuire, M. S.', title: 'Investigation of Growth Parameter Influence on Hydrothermally Grown ZnO Nanowires Using a Research Grade Microwave', journal: 'Author manuscript hosted by ON Semiconductor (a two-page contribution with the same title and authors is in the 2009 International Semiconductor Device Research Symposium, doi 10.1109/ISDRS.2009.5378342)', year: 2009, doi: null, url: 'https://www.onsemi.com/site/pdf/SSE_MicrowaveZnO.pdf' },
    conditions: '90 nm ALD ZnO seed layer on Si, 20 mL of 25 mM Zn(NO₃)₂ + 25 mM HMT (1:1) in a sealed 35 mL vial, 85 °C held by a temperature-controlled microwave reactor, slow stirring, 60 min',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Mason, Roberts, Conley, Price, Allman & McGuire, hydrothermal ZnO nanowires grown in a research-grade microwave reactor (author manuscript, ON Semiconductor / Oregon State University)',
    note: '60-minute growth at 85 °C and 25 mM as retrieved from the source (696 nm). Model 574 nm (−17 %); the previous growth law, paced by HMTA hydrolysis alone, gave 60–70 % too little here. The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). The fast start comes from the retrograde solubility of ZnO in the HMTA-buffered bath, not from base release. The microwave reactor heats faster than the default 600 s time constant, which is kept. Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series.'
  },
  {
    id: 's07-exp-zno-nanorod-length-0p1M-80c-4h', suite: 7, title: 'ZnO nanorods from 0.1 mol/L zinc nitrate + HMT at 80 °C after 4 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 100, cHMTA: 100, cOH: 0, T: 80, tmin: 240, Vb: 1e-4 }, kpi: 'L_final', expected: 770, tol: 0.3, kind: 'experiment',
    measured: { value: 770, uncertainty: null, unit: 'nm' },
    source: { authors: 'Ridhuan, N. S.; Abdul Razak, K.; Lockman, Z.; Abdul Aziz, A.', title: 'Structural and Morphology of ZnO Nanorods Synthesized Using ZnO Seeded Growth Hydrothermal Method and Its Properties as UV Sensing', journal: 'PLoS ONE 7 (11), e50405', year: 2012, doi: '10.1371/journal.pone.0050405', url: 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050405' },
    conditions: 'About 200 nm RF-sputtered ZnO seed layer on Si(100), annealed at 400 °C; equimolar 0.1 M Zn(NO₃)₂ + HMT (solution pH 7) in a screw-capped bottle placed in an oven preheated to 80 °C; 4 h; lengths from FESEM cross-sections (averages, no uncertainty stated)',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Ridhuan, Abdul Razak, Lockman & Abdul Aziz, PLoS ONE 7, e50405 (2012), section 3.4 (growth-duration series)',
    note: 'Section 3.4 of the source: "average lengths of 770 nm, 865 nm, 948 nm, and 1.1 µm for growth durations of 4, 8, 16, and 24 h, respectively". Model 699 nm (−9 %). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). A four-fold higher concentration than the calibration baths: the old growth law (rate proportional to concentration) gave several micrometres here. The bath volume is not stated (0.1 L assumed; it has no influence because the nanorods take up under 0.1 % of the zinc). Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series. From the same source, 0.05 M for 4 h gave 321 nm, which the model does NOT reproduce (680 nm, +112 %): it predicts almost no difference between 0.05 and 0.1 M, whereas this laboratory found a factor 2.4; that point is not offered as a validation entry.'
  },
  {
    id: 's07-exp-zno-nanorod-length-0p1M-80c-8h', suite: 7, title: 'ZnO nanorods from 0.1 mol/L zinc nitrate + HMT at 80 °C after 8 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 100, cHMTA: 100, cOH: 0, T: 80, tmin: 480, Vb: 1e-4 }, kpi: 'L_final', expected: 865, tol: 0.3, kind: 'experiment',
    measured: { value: 865, uncertainty: null, unit: 'nm' },
    source: { authors: 'Ridhuan, N. S.; Abdul Razak, K.; Lockman, Z.; Abdul Aziz, A.', title: 'Structural and Morphology of ZnO Nanorods Synthesized Using ZnO Seeded Growth Hydrothermal Method and Its Properties as UV Sensing', journal: 'PLoS ONE 7 (11), e50405', year: 2012, doi: '10.1371/journal.pone.0050405', url: 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050405' },
    conditions: 'About 200 nm RF-sputtered ZnO seed layer on Si(100), annealed at 400 °C; equimolar 0.1 M Zn(NO₃)₂ + HMT (solution pH 7) in a screw-capped bottle placed in an oven preheated to 80 °C; 8 h; lengths from FESEM cross-sections (averages, no uncertainty stated)',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Ridhuan, Abdul Razak, Lockman & Abdul Aziz, PLoS ONE 7, e50405 (2012), section 3.4 (growth-duration series)',
    note: 'Section 3.4 of the source: "average lengths of 770 nm, 865 nm, 948 nm, and 1.1 µm for growth durations of 4, 8, 16, and 24 h, respectively". Model 779 nm (−10 %). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). A four-fold higher concentration than the calibration baths: the old growth law (rate proportional to concentration) gave several micrometres here. The bath volume is not stated (0.1 L assumed; it has no influence because the nanorods take up under 0.1 % of the zinc). Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series. From the same source, 0.05 M for 4 h gave 321 nm, which the model does NOT reproduce (680 nm, +112 %): it predicts almost no difference between 0.05 and 0.1 M, whereas this laboratory found a factor 2.4; that point is not offered as a validation entry.'
  },
  {
    id: 's07-exp-zno-nanorod-length-0p1M-80c-16h', suite: 7, title: 'ZnO nanorods from 0.1 mol/L zinc nitrate + HMT at 80 °C after 16 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 100, cHMTA: 100, cOH: 0, T: 80, tmin: 960, Vb: 1e-4 }, kpi: 'L_final', expected: 948, tol: 0.3, kind: 'experiment',
    measured: { value: 948, uncertainty: null, unit: 'nm' },
    source: { authors: 'Ridhuan, N. S.; Abdul Razak, K.; Lockman, Z.; Abdul Aziz, A.', title: 'Structural and Morphology of ZnO Nanorods Synthesized Using ZnO Seeded Growth Hydrothermal Method and Its Properties as UV Sensing', journal: 'PLoS ONE 7 (11), e50405', year: 2012, doi: '10.1371/journal.pone.0050405', url: 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050405' },
    conditions: 'About 200 nm RF-sputtered ZnO seed layer on Si(100), annealed at 400 °C; equimolar 0.1 M Zn(NO₃)₂ + HMT (solution pH 7) in a screw-capped bottle placed in an oven preheated to 80 °C; 16 h; lengths from FESEM cross-sections (averages, no uncertainty stated)',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Ridhuan, Abdul Razak, Lockman & Abdul Aziz, PLoS ONE 7, e50405 (2012), section 3.4 (growth-duration series)',
    note: 'Section 3.4 of the source: "average lengths of 770 nm, 865 nm, 948 nm, and 1.1 µm for growth durations of 4, 8, 16, and 24 h, respectively". Model 882 nm (−7 %). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). A four-fold higher concentration than the calibration baths: the old growth law (rate proportional to concentration) gave several micrometres here. The bath volume is not stated (0.1 L assumed; it has no influence because the nanorods take up under 0.1 % of the zinc). Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series. From the same source, 0.05 M for 4 h gave 321 nm, which the model does NOT reproduce (680 nm, +112 %): it predicts almost no difference between 0.05 and 0.1 M, whereas this laboratory found a factor 2.4; that point is not offered as a validation entry.'
  },
  {
    id: 's07-exp-zno-nanorod-length-0p1M-80c-24h', suite: 7, title: 'ZnO nanorods from 0.1 mol/L zinc nitrate + HMT at 80 °C after 24 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 100, cHMTA: 100, cOH: 0, T: 80, tmin: 1440, Vb: 1e-4 }, kpi: 'L_final', expected: 1100, tol: 0.3, kind: 'experiment',
    measured: { value: 1100, uncertainty: null, unit: 'nm' },
    source: { authors: 'Ridhuan, N. S.; Abdul Razak, K.; Lockman, Z.; Abdul Aziz, A.', title: 'Structural and Morphology of ZnO Nanorods Synthesized Using ZnO Seeded Growth Hydrothermal Method and Its Properties as UV Sensing', journal: 'PLoS ONE 7 (11), e50405', year: 2012, doi: '10.1371/journal.pone.0050405', url: 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0050405' },
    conditions: 'About 200 nm RF-sputtered ZnO seed layer on Si(100), annealed at 400 °C; equimolar 0.1 M Zn(NO₃)₂ + HMT (solution pH 7) in a screw-capped bottle placed in an oven preheated to 80 °C; 24 h; lengths from FESEM cross-sections (averages, no uncertainty stated)',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Ridhuan, Abdul Razak, Lockman & Abdul Aziz, PLoS ONE 7, e50405 (2012), section 3.4 (growth-duration series)',
    note: 'Section 3.4 of the source: "average lengths of 770 nm, 865 nm, 948 nm, and 1.1 µm for growth durations of 4, 8, 16, and 24 h, respectively". Model 927 nm (−16 %). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). A four-fold higher concentration than the calibration baths: the old growth law (rate proportional to concentration) gave several micrometres here. The bath volume is not stated (0.1 L assumed; it has no influence because the nanorods take up under 0.1 % of the zinc). Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series. From the same source, 0.05 M for 4 h gave 321 nm, which the model does NOT reproduce (680 nm, +112 %): it predicts almost no difference between 0.05 and 0.1 M, whereas this laboratory found a factor 2.4; that point is not offered as a validation entry.'
  },
  {
    id: 's07-exp-zno-nanorod-length-0p1M-90c-10h', suite: 7, title: 'ZnO nanorods from 0.1 mol/L zinc nitrate + HMT at 90 °C after 10 h: length (independent of the calibration)', mode: 'zno',
    inputs: { mode: 'zno', chem: 'zn_nitrate', cZn: 100, cHMTA: 100, cOH: 0, T: 90, tmin: 600, Vb: 1e-4 }, kpi: 'L_final', expected: 1200, tol: 0.3, kind: 'experiment',
    measured: { value: 1200, uncertainty: null, unit: 'nm' },
    source: { authors: 'Chae, K.-W.; Zhang, Q.; Kim, J. S.; Jeong, Y.-H.; Cao, G.', title: 'Low-temperature solution growth of ZnO nanotube arrays', journal: 'Beilstein J. Nanotechnol. 1, 128–134', year: 2010, doi: '10.3762/bjnano.1.15', url: 'https://www.beilstein-journals.org/bjnano/articles/1/15' },
    conditions: 'ZnO nanocrystallite seeds (20–50 nm) electrophoretically deposited on ITO glass and annealed at 500 °C; equimolar 0.1 M zinc nitrate + hexamethylenetetramine in a sealed vial heated in a programmable furnace to 90 °C and held for 10 h; length from SEM ("approximately 1.2 µm")',
    match: ['T', 'cZn', 'cHMTA', 'tmin'],
    reference: 'Chae, Zhang, Kim, Jeong & Cao, Beilstein J. Nanotechnol. 1, 128 (2010)',
    note: 'Text of the source: "The length of the ZnO nanorods observed by SEM was approximately 1.2 µm" (diameter about 200 nm). Model 921 nm (−23 %): the bath is exhausted and the length is at the model\'s ceiling, which is set by the calibrated particle yield (about 0.93 µm). The measured diameter of about 200 nm is not reproduced (80 nm with the default facet ratio and seed density; the electrophoretic seeds of this source are much sparser than the default). The suite grows rods and competing bulk ZnO particles under the same driving force, V_m k (c − c_eq), with c_eq the ZnO solubility at the buffered bath pH; its three rate parameters were calibrated on the Heliyon tables only (see s07-zno-nanorod-length-3h). Bath volume not stated (0.1 L assumed, without influence). Seed density, seeded area, boundary-layer thickness and heat-up time constant are not given in the source and stay at the suite defaults. Tolerance 30 %: lengths reported by different laboratories for nominally similar baths differ by this much (seed layer, vessel, heating), the sources give no uncertainty for these points, and the model is a reduced one calibrated on a single laboratory\'s series. A second sample of the same source, 0.1 M at 60 °C for 24 h, is described as rods of about 400 nm length that are wider than long; the model gives 655 nm and that point is not offered as a validation entry.'
  },
];
