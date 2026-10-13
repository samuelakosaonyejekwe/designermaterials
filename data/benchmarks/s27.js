// Reference benchmarks for Suite 27 — safety, standards, environment and device integration.
// Expected values are handbook data or closed-form results evaluated here from standard constants;
// none of them is taken from the suite's own output.
const R = 8.314462618, KB = 1.380649e-23, FARADAY = 96485.33212, PATM = 101325, G0 = 9.80665;

// --- occupational exposure: well-mixed room, G = 1 mg/min, no capture, Q = (ach/60)·V
const EXPO = { mode: 'exposure', G: 1, lev: 0, Vroom: 100, ach: 6, recirc: 0, kdep: 0, Vnf: 1, beta: 5 };
const Q_ROOM = 6 / 60 * 100;                                          // m³/min
// one-box TWA with a slowly ventilated room (λ = 0.01 1/min), 120 min task in a 480 min shift
const LAM_S = 0.6 / 60, Q_S = LAM_S * 100, T_ON = 120, E_ON = 1 - Math.exp(-LAM_S * T_ON);
const twaBox = (1 / Q_S) * (T_ON - E_ON / LAM_S + E_ON * (1 - Math.exp(-LAM_S * (480 - T_ON))) / LAM_S) / 480;
// two-zone model, source on for the whole shift: c(t) = (I − e^{At})·c_ss, so ∫c dt = T·c_ss + A⁻¹·(I − e^{AT})·c_ss.
// The slowest eigenvalue is about Q/V = 0.1 1/min, so e^{AT} (< 1e-20 at T = 480 min) is dropped.
const twaNear = (() => {
  const Vn = 1, Vf = 99, b = 5, a11 = -b / Vn, a12 = b / Vn, a21 = b / Vf, a22 = -(b + Q_ROOM) / Vf, det = a11 * a22 - a12 * a21;
  const cF = 1 / Q_ROOM, cN = cF + 1 / b;
  return (480 * cN + (a22 * cN - a12 * cF) / det) / 480;
})();

// --- aerosol mechanics in air at 20 °C and 1 atm
const T20 = 293.15, MU_AIR = 1.813e-5, MFP = 6.64e-8;                 // handbook viscosity (Pa·s) and mean free path (m)
const slip = d => { const kn = 2 * MFP / d; return 1 + kn * (1.257 + 0.4 * Math.exp(-1.1 / kn)); };
const dAl = 0.02, lnD = Math.log(dAl);                                // 20 nm in µm
const dfAl = (0.0155 / dAl) * (Math.exp(-0.416 * (lnD + 2.84) ** 2) + 19.11 * Math.exp(-0.482 * (lnD - 1.362) ** 2));

// --- water on the saturation line
// Wagner–Pruss auxiliary vapour-pressure equation (an independent formulation from the IF97 region-4 equation)
const psatWP = T => { const Tc = 647.096, th = 1 - T / Tc; return 22.064e6 * Math.exp(Tc / T * (-7.85951783 * th + 1.84408259 * th ** 1.5 - 11.7866497 * th ** 3 + 22.6807411 * th ** 3.5 - 15.9618719 * th ** 4 + 1.80122502 * th ** 7.5)); };
// steam-table densities (kg/m³): liquid at 20 °C, saturated liquid and vapour at 200 °C
const RHO20 = 998.2, RHOL200 = 864.7, RHOV200 = 7.86;
const FILL = 0.7, PHI200 = (FILL * RHO20 - RHOV200) / (RHOL200 - RHOV200);
const P200 = psatWP(473.15) + PATM * (473.15 / T20) * (1 - FILL) / (1 - PHI200);   // vapour + trapped air, Pa absolute
const PG200 = P200 - PATM;
const K2 = (33 / 25) ** 2;                                            // default shell: 50 mm bore, 8 mm wall

// --- hazards
const V_H2 = 1 / (2 * FARADAY) * R * T20 / PATM * 6e7;                // mL/min of hydrogen per ampere at 20 °C, 1 atm

// --- device integration (suite default film stack)
const SIG_DEF = 1.4e11 / (1 - 0.3) * (2.6e-6 - 4.7e-6) * (25 - 300);  // Pa
const BOW_DEF = 6 * SIG_DEF * 5e-7 * (1 - 0.28) / (1.3e11 * 5.25e-4 ** 2) * 0.1 ** 2 / 8;
const TJ_DEF = 25 + 2 * (5e-7 / (30 * 1e-4) + 5.25e-4 / (150 * 1e-4) + 50e-6 / (3 * 1e-4) + 5);

// --- sources of the published-experiment entries
const SRC_KIM = { authors: 'Kim, Mulholland, Kukuck, Pui', title: 'Slip Correction Measurements of Certified PSL Nanoparticles Using a Nanometer Differential Mobility Analyzer (Nano-DMA) for Knudsen Number From 0.5 to 83', journal: 'J. Res. Natl Inst. Stand. Technol.', year: 2005, doi: '10.6028/jres.110.005', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4849562/' };
const SRC_OSBORNE = { authors: 'Osborne, Stimson, Fiock, Ginnings', title: 'The pressure of saturated water vapor in the range 100° to 374° C', journal: 'Bur. Stand. J. Res.', year: 1933, doi: '10.6028/jres.010.014', url: 'https://nvlpubs.nist.gov/nistpubs/jres/10/jresv10n2p155_A2b.pdf' };
const EV = rows => 'event,probability\n' + rows.map(r => r.join(',')).join('\n');

export default [
  // ---------- occupational exposure ----------
  {
    id: 's27-one-box-steady-state', suite: 27, title: 'Well-mixed room: steady-state concentration C = G/Q', mode: 'exposure',
    inputs: { ...EXPO }, kpi: 'c_ff', expected: 1 / Q_ROOM, tol: 1e-6, kind: 'analytical',
    reference: 'Well-mixed (one-box) room model of occupational hygiene: C_ss = G/Q',
    note: 'G = 1 mg/min into 100 m³ at 6 air changes per hour (Q = 10 m³/min); no capture, recirculation or deposition.'
  },
  {
    id: 's27-near-field-steady-state', suite: 27, title: 'Near-field/far-field model: steady near-field concentration', mode: 'exposure',
    inputs: { ...EXPO }, kpi: 'c_nf', expected: 1 / Q_ROOM + 1 / 5, tol: 1e-6, kind: 'analytical',
    reference: 'Two-zone (near-field/far-field) exposure model: C_N,ss = G/Q + G/β',
    note: 'Inter-zone air flow β = 5 m³/min.'
  },
  {
    id: 's27-one-box-twa-transient', suite: 27, title: 'One-box model: 8-h TWA with build-up and decay', mode: 'exposure',
    inputs: { ...EXPO, ach: 0.6, tTask: T_ON }, kpi: 'twa_box', expected: twaBox, tol: 1e-3, kind: 'analytical',
    reference: 'Analytical solution of the first-order well-mixed room balance V·dC/dt = G − Q·C, integrated over the shift',
    note: 'Slow ventilation (0.6 air changes per hour) so that the room neither reaches steady state during the 120 min task nor clears by the end of the shift.'
  },
  {
    id: 's27-two-zone-twa-full-shift', suite: 27, title: 'Near-field/far-field model: 8-h TWA with the source on all shift', mode: 'exposure',
    inputs: { ...EXPO, tTask: 480 }, kpi: 'twa', expected: twaNear, tol: 5e-3, kind: 'analytical',
    reference: 'Matrix-exponential solution of the linear two-zone exposure model',
    note: 'Checks the Crank–Nicolson time integration (0.5 min step) against the exact integral of the two-exponential solution.'
  },
  {
    id: 's27-clearance-time', suite: 27, title: 'Time to clear 99 % of an airborne contaminant', mode: 'exposure',
    inputs: { ...EXPO }, kpi: 't_clear', expected: Math.log(100) / 0.1, tol: 1e-6, kind: 'analytical',
    reference: 'Exponential dilution of a well-mixed room: t = ln(C₀/C)/λ',
    note: '6 air changes per hour, i.e. λ = 0.1 per minute.'
  },

  // ---------- aerosol behaviour ----------
  {
    id: 's27-cunningham-50nm', suite: 27, title: 'Cunningham slip correction of a 50 nm particle in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 5e-8, Tair: 20 }, kpi: 'cc', expected: slip(5e-8), tol: 2e-3, kind: 'analytical',
    reference: 'Cunningham slip correction with the Davies coefficients (1.257, 0.4, 1.1) and a mean free path of air of 66.4 nm at 20 °C',
    note: 'Evaluates to about 5.0, the value tabulated in aerosol textbooks for 0.05 µm.'
  },
  {
    id: 's27-stokes-einstein-50nm', suite: 27, title: 'Stokes–Einstein diffusivity of a 50 nm particle in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 5e-8, Tair: 20 }, kpi: 'D', expected: KB * T20 * slip(5e-8) / (3 * Math.PI * MU_AIR * 5e-8), tol: 5e-3, kind: 'analytical',
    reference: 'Stokes–Einstein relation with slip correction, D = k_B·T·C_c/(3πμd); air viscosity 1.813×10⁻⁵ Pa·s at 20 °C',
    note: 'The suite obtains viscosity from the Sutherland formula; the reference uses the handbook value.'
  },
  {
    id: 's27-settling-1um-unit-density', suite: 27, title: 'Settling velocity of a 1 µm unit-density sphere in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 1e-6, rhoP: 1000, Tair: 20 }, kpi: 'vs', expected: 3.5e-5, tol: 0.02, kind: 'reference-data',
    reference: 'Hinds, Aerosol Technology — tabulated terminal settling velocity of unit-density spheres at 20 °C',
    note: 'Tabulated to two significant figures.'
  },
  {
    id: 's27-coagulation-continuum-limit', suite: 27, title: 'Brownian coagulation coefficient: continuum limit', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 1e-5, rhoP: 1000, Tair: 20 }, kpi: 'K', expected: 8 * KB * T20 * slip(1e-5) / (3 * MU_AIR), tol: 0.01, kind: 'analytical',
    reference: 'Smoluchowski continuum coagulation coefficient of equal spheres, K = 8·k_B·T·C_c/(3μ) (convention dN/dt = −½KN²)',
    note: 'At 10 µm the Fuchs transition correction used by the suite is still about 0.6 % below the pure continuum limit.'
  },
  {
    id: 's27-coagulation-free-molecular-limit', suite: 27, title: 'Brownian coagulation coefficient: free-molecular limit', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 1e-9, rhoP: 1000, Tair: 20 }, kpi: 'K', expected: 4 * Math.sqrt(6 * KB * T20 * 1e-9 / 1000), tol: 0.01, kind: 'analytical',
    reference: 'Kinetic-theory collision rate of equal spheres, K = √2·π·d²·c̄ = 4·√(6·k_B·T·d/ρ) (convention dN/dt = −½KN²)',
    note: 'At 1 nm (Knudsen number above 100) the Fuchs interpolation must reduce to the kinetic-theory value.'
  },
  {
    id: 's27-alveolar-deposition-20nm', suite: 27, title: 'Alveolar deposition fraction of 20 nm particles', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 2e-8, Tair: 20 }, kpi: 'df_al', expected: dfAl, tol: 1e-3, kind: 'analytical',
    reference: 'Hinds, Aerosol Technology — algebraic fit of the ICRP Publication 66 deposition model (adult, light exercise, nose breathing)',
    note: 'Evaluates to about 0.48, the alveolar deposition maximum near 20 nm; the reference is the published fit itself, so this verifies its transcription.'
  },

  // ---------- hydrothermal autoclave ----------
  {
    id: 's27-water-psat-100C', suite: 27, title: 'Water saturation pressure at 100 °C', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 100 }, kpi: 'psat', expected: 101.42e3, tol: 5e-4, kind: 'reference-data',
    reference: 'Steam tables (IAPWS formulations on ITS-90): p_sat(100 °C) = 101.42 kPa'
  },
  {
    id: 's27-water-psat-200C', suite: 27, title: 'Water saturation pressure at 200 °C', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 200 }, kpi: 'psat', expected: psatWP(473.15), tol: 5e-4, kind: 'reference-data',
    reference: 'Wagner–Pruss auxiliary vapour-pressure equation of the IAPWS scientific formulation (steam tables list 1.5549 MPa)',
    note: 'The suite uses the IAPWS-IF97 industrial equation, which gives 1.5547 MPa; the two IAPWS formulations differ by about 0.02 % here.'
  },
  {
    id: 's27-water-psat-if97-500K', suite: 27, title: 'IAPWS-IF97 saturation pressure: verification value at 500 K', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 226.85 }, kpi: 'psat', expected: 2.63889776e6, tol: 1e-7, kind: 'numerical-reference',
    reference: 'IAPWS Industrial Formulation 1997, region 4 — computer-program verification value p_s(500 K) = 2.63889776 MPa'
  },
  {
    id: 's27-critical-fill-200C', suite: 27, title: 'Fill fraction that makes a sealed vessel liquid-full at 200 °C', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 200 }, kpi: 'fcrit', expected: RHOL200 / RHO20, tol: 2e-3, kind: 'reference-data',
    reference: 'Steam tables: saturated liquid density 864.7 kg/m³ at 200 °C and 998.2 kg/m³ at 20 °C',
    note: 'Mass conservation of the liquid charge: f_crit = ρ′(T)/ρ(20 °C).'
  },
  {
    id: 's27-autoclave-total-pressure', suite: 27, title: 'Sealed autoclave at 200 °C, 70 % fill: vapour plus trapped air', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 200, fillPct: 70 }, kpi: 'P', expected: P200, tol: 5e-3, kind: 'analytical',
    reference: 'Lever rule for the two-phase charge with steam-table densities, Dalton and ideal-gas laws for the air sealed in at 20 °C and 1 atm',
    note: 'Air solubility and vessel expansion neglected in both the suite and the reference.'
  },
  {
    id: 's27-lame-hoop-stress', suite: 27, title: 'Lamé hoop stress at the bore of a thick-walled cylinder', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 200, fillPct: 70, IDmm: 50, wall: 8 }, kpi: 's_hoop', expected: PG200 * (K2 + 1) / (K2 - 1), tol: 5e-3, kind: 'analytical',
    reference: 'Lamé solution for an internally pressurised thick cylinder: σ_θ(rᵢ) = p·(rₒ² + rᵢ²)/(rₒ² − rᵢ²)',
    note: 'Gauge pressure from the steam-table reference above; the tolerance covers the property data, the stress formula itself is exact.'
  },
  {
    id: 's27-thin-wall-limit', suite: 27, title: 'Thin-walled limit of the hoop stress: p·r/t', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 200, fillPct: 70, IDmm: 1000, wall: 5 }, kpi: 's_hoop', expected: PG200 * 502.5 / 5, tol: 5e-3, kind: 'analytical',
    reference: 'Thin-walled pressure-vessel (membrane) formula σ_θ = p·r/t with the mean radius',
    note: 'Wall/radius = 0.01, where the Lamé bore stress and the membrane stress differ by less than 0.01 %.'
  },
  {
    id: 's27-thick-cylinder-first-yield', suite: 27, title: 'Allowable pressure of a closed-end thick cylinder (von Mises first yield)', mode: 'autoclave',
    inputs: { mode: 'autoclave', IDmm: 50, wall: 8, yieldT: 130, designSF: 1.5 }, kpi: 'mawp', expected: 130e6 * (K2 - 1) / (Math.sqrt(3) * K2) / 1.5, tol: 1e-6, kind: 'analytical',
    reference: 'Lamé stresses with the von Mises criterion at the bore of a closed-end cylinder: p_y = σ_y·(K² − 1)/(√3·K²)'
  },

  // ---------- flammable solvent and electrical screening ----------
  {
    id: 's27-ethanol-normal-boiling-point', suite: 27, title: 'Ethanol vapour pressure at its normal boiling point', mode: 'hazards',
    inputs: { mode: 'hazards', solvent: 'ethanol', Tsolv: 78.37 }, kpi: 'psat', expected: PATM, tol: 0.01, kind: 'reference-data',
    reference: 'Normal boiling point of ethanol, 78.37 °C (CRC Handbook of Chemistry and Physics)',
    note: 'Tests the Antoine constants at the upper end of their range.'
  },
  {
    id: 's27-ethanol-vapour-pressure-25C', suite: 27, title: 'Ethanol vapour pressure at 25 °C', mode: 'hazards',
    inputs: { mode: 'hazards', solvent: 'ethanol', Tsolv: 25 }, kpi: 'psat', expected: 7870, tol: 0.02, kind: 'reference-data',
    reference: 'Handbook vapour pressure of ethanol at 25 °C, 7.87 kPa (59 mmHg)'
  },
  {
    id: 's27-jones-rule-ethanol', suite: 27, title: 'Lower flammable limit of ethanol from stoichiometry (Jones rule)', mode: 'hazards',
    inputs: { mode: 'hazards', solvent: 'ethanol' }, kpi: 'lel_jones', expected: 0.55 * 100 / (1 + 3 / 0.21), tol: 5e-3, kind: 'analytical',
    reference: 'Jones rule LFL ≈ 0.55·C_st with C₂H₅OH + 3 O₂ → 2 CO₂ + 3 H₂O in air of 21 vol % oxygen',
    note: 'The rule itself is an empirical estimate (tabulated limit 3.3 vol %); the benchmark checks the stoichiometric concentration.'
  },
  {
    id: 's27-ohmic-cell-current', suite: 27, title: 'Ohmic current of a parallel-plate electrolyte cell', mode: 'hazards',
    inputs: { mode: 'hazards', Uapp: 25, gap: 10, Ael: 4, sigEl: 1 }, kpi: 'I', expected: 1 * 4e-4 * 25 / 0.01, tol: 1e-6, kind: 'analytical',
    reference: 'Ohm\'s law for a uniform conductor: I = σ·A·U/g'
  },
  {
    id: 's27-faraday-hydrogen', suite: 27, title: 'Electrolytic hydrogen evolution at 1 A', mode: 'hazards',
    inputs: { mode: 'hazards', Uapp: 25, gap: 10, Ael: 4, sigEl: 1, acdc: 'dc' }, kpi: 'h2', expected: V_H2, tol: 1e-4, kind: 'analytical',
    reference: 'Faraday\'s law of electrolysis (two electrons per H₂) with the ideal-gas molar volume at 20 °C and 1 atm',
    note: 'About 7.48 mL/min per ampere at 100 % current efficiency.'
  },

  // ---------- risk analysis ----------
  {
    id: 's27-rpn-product', suite: 27, title: 'FMEA risk priority number S × O × D', mode: 'risk',
    inputs: { mode: 'risk', fmeaSet: 'table', fmea: 'failure_mode,severity,occurrence,detection\nSeal leak,7,5,6\nLabel missing,2,2,2' }, kpi: 'rpn_max', expected: 210, abs: 1e-9, kind: 'analytical',
    reference: 'Failure mode and effects analysis: RPN = severity × occurrence × detection'
  },
  {
    id: 's27-fault-tree-and-gate', suite: 27, title: 'Fault tree: AND gate of two independent events', mode: 'risk',
    inputs: { mode: 'risk', tree: 'a*b', events: EV([['a', 0.02], ['b', 0.05]]) }, kpi: 'p_top', expected: 0.02 * 0.05, tol: 1e-9, kind: 'analytical',
    reference: 'Probability of the intersection of independent events: P = Π pᵢ'
  },
  {
    id: 's27-fault-tree-or-gate', suite: 27, title: 'Fault tree: OR gate of three independent events', mode: 'risk',
    inputs: { mode: 'risk', tree: 'a + b + c', events: EV([['a', 0.1], ['b', 0.2], ['c', 0.3]]) }, kpi: 'p_top', expected: 1 - 0.9 * 0.8 * 0.7, tol: 1e-9, kind: 'analytical',
    reference: 'Probability of the union of independent events: P = 1 − Π(1 − pᵢ)'
  },
  {
    id: 's27-fault-tree-two-out-of-three', suite: 27, title: 'Fault tree with repeated events: two-out-of-three voting', mode: 'risk',
    inputs: { mode: 'risk', tree: 'a*b + b*c + a*c', events: EV([['a', 0.1], ['b', 0.1], ['c', 0.1]]) }, kpi: 'p_top', expected: 3 * 0.1 ** 2 - 2 * 0.1 ** 3, tol: 1e-9, kind: 'analytical',
    reference: 'Reliability of a two-out-of-three system of identical independent components: P = 3p² − 2p³',
    note: 'Each event appears in two minimal cut sets, so gate-by-gate multiplication would be wrong; the exact value requires state enumeration.'
  },

  // ---------- environmental fate and effluent ----------
  {
    id: 's27-pec-dilution', suite: 27, title: 'Predicted environmental concentration by river dilution', mode: 'envfate',
    inputs: { mode: 'envfate', mProd: 1000, wNP: 1, relFrac: 10, wwtp: 90, Qriver: 10 }, kpi: 'pec', expected: 1000e3 * 0.01 * 0.1 * 0.1 / (10 * 365.25 * 86400) * 1e6, tol: 1e-3, kind: 'conservation',
    reference: 'Steady mass balance of a fully mixed river: PEC = emission/flow',
    note: '100 kg/yr reaching 10 m³/s; 1 kg/m³ = 10⁶ µg/L.'
  },
  {
    id: 's27-neutralisation-acid-demand', suite: 27, title: 'Sulfuric acid needed to neutralise 100 mol of NaOH', mode: 'envfate',
    inputs: { mode: 'envfate', wasteType: 'naoh_h2so4', wasteVol: 1000, wasteConc: 0.1 }, kpi: 'm_neut', expected: 50 * 98.08e-3, tol: 1e-3, kind: 'analytical',
    reference: 'Stoichiometry of 2 NaOH + H₂SO₄ → Na₂SO₄ + 2 H₂O with M(H₂SO₄) = 98.08 g/mol'
  },
  {
    id: 's27-neutralisation-salt-formed', suite: 27, title: 'Sodium sulfate formed from 100 mol of NaOH', mode: 'envfate',
    inputs: { mode: 'envfate', wasteType: 'naoh_h2so4', wasteVol: 1000, wasteConc: 0.1 }, kpi: 'm_salt', expected: 50 * 142.04e-3, tol: 1e-3, kind: 'analytical',
    reference: 'Stoichiometry of 2 NaOH + H₂SO₄ → Na₂SO₄ + 2 H₂O with M(Na₂SO₄) = 142.04 g/mol'
  },
  {
    id: 's27-zno-solubility-ph7', suite: 27, title: 'Equilibrium solubility of ZnO at pH 7', mode: 'envfate',
    inputs: { mode: 'envfate', nm: 'ZnO', pHw: 7 }, kpi: 'c_eq', expected: Math.pow(10, 11.2 - 14) * 65.38 * 1e3, tol: 1e-3, kind: 'analytical',
    reference: 'Mass action for ZnO(s) + 2 H⁺ ⇌ Zn²⁺ + H₂O with the handbook constant log K ≈ 11.2 at 25 °C',
    note: 'Free Zn²⁺ only (no hydroxo complexes, activity coefficients of 1), as documented by the suite; the benchmark checks the equilibrium expression and unit conversion, not the constant.'
  },
  {
    id: 's27-stokes-settling-in-water', suite: 27, title: 'Stokes settling of a 100 nm silver particle in water', mode: 'envfate',
    inputs: { mode: 'envfate', nm: 'Ag', dp: 1e-7 }, kpi: 'v_set', expected: (10490 - 998.2) * G0 * 1e-7 ** 2 / (18 * 1.002e-3), tol: 5e-3, kind: 'analytical',
    reference: 'Stokes law v = (ρ_p − ρ_w)·g·d²/(18μ) with water at 20 °C (998.2 kg/m³, 1.002 mPa·s) and silver 10 490 kg/m³'
  },

  // ---------- device integration ----------
  {
    id: 's27-thermal-mismatch-stress', suite: 27, title: 'Thermal-mismatch stress of a ZnO film on silicon', mode: 'device',
    inputs: { mode: 'device' }, kpi: 'sigma', expected: SIG_DEF, tol: 1e-6, kind: 'analytical',
    reference: 'Equi-biaxial thin-film thermal stress σ = E_f/(1 − ν_f)·(α_s − α_f)·ΔT',
    note: 'Cooling by 275 K with α_f > α_s gives 115.5 MPa tension.'
  },
  {
    id: 's27-stoney-bow', suite: 27, title: 'Stoney wafer bow under a stressed film', mode: 'device',
    inputs: { mode: 'device' }, kpi: 'bow', expected: BOW_DEF, tol: 1e-6, kind: 'analytical',
    reference: 'Stoney equation κ = 6·σ_f·t_f·(1 − ν_s)/(E_s·t_s²) with spherical-cap bow δ = κ·D²/8',
    note: '500 nm film on a 525 µm, 100 mm silicon wafer.'
  },
  {
    id: 's27-auxetic-nte-film-stress', suite: 27, title: 'Auxetic film with negative thermal expansion: mismatch stress', mode: 'device',
    inputs: { mode: 'device', Ef: 1e10, nuf: -0.5, ctef: -5e-6, ctes: 2.6e-6, Tdep: 200, Tuse: 25 }, kpi: 'sigma', expected: 1e10 / 1.5 * (2.6e-6 + 5e-6) * (25 - 200), tol: 1e-6, kind: 'analytical',
    reference: 'Equi-biaxial thin-film thermal stress σ = E_f/(1 − ν_f)·(α_s − α_f)·ΔT',
    note: 'ν_f = −0.5 lowers the biaxial modulus to E/1.5, and the film expands on cooling, so the stress is compressive (negative).'
  },
  {
    id: 's27-sheet-resistance', suite: 27, title: 'Sheet resistance of a 100 nm transparent conductor', mode: 'device',
    inputs: { mode: 'device', sigf: 1e6, tf: 1e-7 }, kpi: 'Rsheet', expected: 10, tol: 1e-9, kind: 'analytical',
    reference: 'Sheet resistance R_□ = 1/(σ·t)'
  },
  {
    id: 's27-junction-temperature', suite: 27, title: 'Junction temperature of a series thermal stack', mode: 'device',
    inputs: { mode: 'device' }, kpi: 'Tj', expected: TJ_DEF, tol: 1e-6, kind: 'analytical',
    reference: 'One-dimensional steady conduction: T_j = T_amb + P·(Σ tᵢ/(kᵢ·A) + R_sa)'
  },
  {
    id: 's27-moisture-ingress-half-saturation', suite: 27, title: 'Fickian moisture ingress: time to 50 % saturation at a sealed face', mode: 'device',
    inputs: { mode: 'device', Lenc: 1000, Dm: 1e-13, cThr: 50 }, kpi: 't_moist', expected: 4 / Math.PI ** 2 * Math.log(8 / Math.PI) * 1e-3 ** 2 / 1e-13, tol: 1e-3, kind: 'analytical',
    reference: 'Crank, The Mathematics of Diffusion — plane sheet with one face at constant concentration and one impermeable face',
    note: 'Leading-term solution t = (4/π²)·ln(8/π)·L²/D; the neglected higher terms change the result by about 0.02 %.'
  },
  // ---------- published experiments ----------
  {
    id: 's27-exp-slip-correction-psl-101nm', suite: 27, title: 'Experiment: slip correction factor of 100.7 nm polystyrene-latex spheres in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 100.7e-9, Tair: 22.35, rhoP: 1050 }, kpi: 'cc', expected: 2.897, tol: 0.03, kind: 'experiment',
    measured: { value: 2.897, uncertainty: 0.058, unit: '–' },
    source: SRC_KIM,
    conditions: 'Certified 100.7 nm polystyrene-latex spheres (standard reference material), electrical mobility measured with a nano-differential mobility analyser in air at 98.80 kPa and 295.5 K (Knudsen number 1.367)',
    match: ['dp', 'Tair'],
    reference: 'Kim, Mulholland, Kukuck & Pui, J. Res. Natl Inst. Stand. Technol. 110, 31 (2005)',
    note: 'Table 7a of the source, first row (98.80 kPa, 295.5 K). Expanded (95 %) relative uncertainty about 2 % (abstract). The suite evaluates the slip correction at 101.325 kPa (pressure is not an input); the 2.5 % lower pressure of the measurement raises the measured value by about 1.4 %, which is included in the tolerance.'
  },
  {
    id: 's27-exp-slip-correction-psl-269nm', suite: 27, title: 'Experiment: slip correction factor of 269 nm polystyrene-latex spheres in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 269e-9, Tair: 23.05, rhoP: 1050 }, kpi: 'cc', expected: 1.631, tol: 0.035, kind: 'experiment',
    measured: { value: 1.631, uncertainty: 0.041, unit: '–' },
    source: SRC_KIM,
    conditions: 'Certified 269 nm polystyrene-latex spheres (standard reference material), electrical mobility measured with a differential mobility analyser in air at 98.50 kPa and 296.2 K (Knudsen number 0.515)',
    match: ['dp', 'Tair'],
    reference: 'Kim, Mulholland, Kukuck & Pui, J. Res. Natl Inst. Stand. Technol. 110, 31 (2005)',
    note: 'Table 7b of the source, first row (a repeat at the same pressure gave 1.642). Expanded (95 %) relative uncertainty about 2.5 % (abstract). The suite evaluates at 101.325 kPa; the pressure difference accounts for about 1 %.'
  },
  {
    id: 's27-exp-slip-correction-psl-20nm', suite: 27, title: 'Experiment: slip correction factor of 19.90 nm polystyrene-latex spheres in air', mode: 'aerosol',
    inputs: { mode: 'aerosol', dp: 19.9e-9, Tair: 23.05, rhoP: 1050 }, kpi: 'cc', expected: 12.042, tol: 0.045, kind: 'experiment',
    measured: { value: 12.042, uncertainty: 0.36, unit: '–' },
    source: SRC_KIM,
    conditions: '19.90 nm polystyrene-latex spheres generated by electrospray, electrical mobility measured with a nano-differential mobility analyser in air at 98.30 kPa and 296.2 K (Knudsen number 6.973)',
    match: ['dp', 'Tair'],
    reference: 'Kim, Mulholland, Kukuck & Pui, J. Res. Natl Inst. Stand. Technol. 110, 31 (2005)',
    note: 'Table 7c of the source, first row (98.30 kPa, 296.2 K). Expanded (95 %) relative uncertainty about 3 % (abstract). At this size the slip correction is nearly proportional to the mean free path, so the 3 % lower pressure of the measurement raises the measured value by almost 3 % relative to the suite, which evaluates at 101.325 kPa; the tolerance covers this offset.'
  },
  {
    id: 's27-exp-water-vapour-pressure-180C', suite: 27, title: 'Experiment: measured pressure of saturated water vapour at 180 °C (static method)', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 180 }, kpi: 'psat', expected: 9.8958 * PATM, tol: 2e-3, kind: 'experiment',
    measured: { value: 9.8958, uncertainty: null, unit: 'standard atm' },
    source: SRC_OSBORNE,
    conditions: 'Pure water in a closed calorimeter, liquid and vapour in equilibrium at 180 °C (International Temperature Scale of 1927); pressure measured with a piston gauge, mean of the repeat determinations reduced to the even temperature',
    match: ['Tac'],
    reference: 'Osborne, Stimson, Fiock & Ginnings, Bur. Stand. J. Res. 10, 155 (1933)',
    note: 'Table 1 of the source, "Mean value of pressure at 180°" (9.8958 atm = 1.00269 MPa). No overall uncertainty is quoted for this entry; the nine individual determinations reduced to 180 °C lie between 9.8920 and 9.8968 atm. The 1927 temperature scale differs from ITS-90 by a few hundredths of a kelvin here, equivalent to about 0.1 % in pressure, which sets the tolerance.'
  },
  {
    id: 's27-exp-water-vapour-pressure-250C', suite: 27, title: 'Experiment: measured pressure of saturated water vapour at 250 °C (static method)', mode: 'autoclave',
    inputs: { mode: 'autoclave', Tac: 250 }, kpi: 'psat', expected: 39.2566 * PATM, tol: 2e-3, kind: 'experiment',
    measured: { value: 39.2566, uncertainty: null, unit: 'standard atm' },
    source: SRC_OSBORNE,
    conditions: 'Pure water in a closed calorimeter, liquid and vapour in equilibrium at 250 °C (International Temperature Scale of 1927); pressure measured with a piston gauge, mean of the repeat determinations reduced to the even temperature',
    match: ['Tac'],
    reference: 'Osborne, Stimson, Fiock & Ginnings, Bur. Stand. J. Res. 10, 155 (1933)',
    note: 'Table 1 of the source, "Mean value of pressure at 250°" (39.2566 atm = 3.9777 MPa). Tolerance as for the 180 °C point (temperature-scale difference of the order of 0.1 % in pressure).'
  },
];
