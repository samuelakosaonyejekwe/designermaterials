// Reference benchmarks for Suite 6 — thermodynamics, chemical equilibrium and phase behaviour.
// Expected values are handbook data or closed-form results evaluated here from standard constants;
// none of them is taken from the suite's own output.
const R = 8.314462618, KB = 1.380649e-23, NA = 6.02214076e23, QE = 1.602176634e-19, LN10 = Math.LN10;

// --- aqueous ZnO: independent mass-action evaluation with the Baes–Mesmer constants at 25 °C, I = 0
const ZN = { logKs: 11.14, logB: [-8.96, -16.9, -28.4, -41.2] };
const znoSol = pH => { let s = Math.pow(10, ZN.logKs - 2 * pH); ZN.logB.forEach((b, i) => { s += Math.pow(10, ZN.logKs - 2 * pH + b + (i + 1) * pH); }); return s; };   // mol/L
let znMin = Infinity, znMinPH = 0;
for (let p = 8; p <= 12; p += 1e-4) { const s = znoSol(p); if (s < znMin) { znMin = s; znMinPH = p; } }

// --- Davies activity coefficient with the classical Debye–Hückel constant for water at 25 °C
const A_DH = 0.509;
const davies = (z, I) => Math.pow(10, -A_DH * z * z * (Math.sqrt(I) / (1 + Math.sqrt(I)) - 0.3 * I));

// --- reaction equilibrium
const dG350 = -50e3 - 350 * -100;                                   // J/mol
const T25 = 298.15, Tk = 600, dCp = 80;
const dGkirch = (-40e3 + dCp * (Tk - T25)) - Tk * (-60 + dCp * Math.log(Tk / T25));

// --- regular solution
const OMEGA = 18e3;
const T_x25 = OMEGA / (2 * R * Math.log(3));                        // temperature at which the binodal is x = 0.25 / 0.75

// --- nanoscale thermodynamics
const au = { g: 0.27, Vm: 196.97e-3 / 19300, Tm: 1337.3, Hf: 12.55e3 };
const zno = { g: 0.25, Vm: 81.38e-3 / 5606 }, Tn = 298.15, S = 50;
const lcapZnO = 2 * zno.g * zno.Vm / (R * Tn);
const vZnO = zno.Vm / NA;
const barZnO = 16 * Math.PI * zno.g ** 3 * vZnO ** 2 / (3 * (KB * Tn * Math.log(S)) ** 2) / (KB * Tn);

export default [
  // ------------------------------------------------------------------ vapour pressure
  {
    id: 's06-water-psat-200C', suite: 6, title: 'Water saturation pressure at 200 °C', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 200 }, kpi: 'Psat', expected: 1.5547, tol: 0.01, kind: 'reference-data',
    reference: 'IAPWS-IF97 steam tables: saturation pressure of water at 200 °C = 1.5547 MPa',
    note: 'The suite uses a two-range Antoine correlation, accurate to a fraction of a per cent over 100–374 °C.'
  },
  {
    id: 's06-water-psat-100C', suite: 6, title: 'Water saturation pressure at the normal boiling point', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 100 }, kpi: 'Psat', expected: 0.101325, tol: 0.005, kind: 'reference-data',
    reference: 'Definition of the normal boiling point: 1 atm = 101 325 Pa (steam tables give 101.42 kPa at 100 °C on ITS-90)'
  },
  {
    id: 's06-water-psat-25C', suite: 6, title: 'Water saturation pressure at 25 °C', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 25 }, kpi: 'Psat', expected: 0.0031698, tol: 0.01, kind: 'reference-data',
    reference: 'Steam tables: saturation pressure of water at 25 °C = 3.1698 kPa'
  },
  {
    id: 's06-water-hvap-100C', suite: 6, title: 'Enthalpy of vaporisation of water at 100 °C', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 100 }, kpi: 'Hvap', expected: 40.65, tol: 0.03, kind: 'reference-data',
    reference: 'Steam tables: latent heat of water at 100 °C = 2257 kJ/kg = 40.65 kJ/mol',
    note: 'The suite applies Clausius–Clapeyron with an ideal vapour to the Antoine slope, which overestimates by about 2 % at 1 atm.'
  },
  {
    id: 's06-water-max-fill-200C', suite: 6, title: 'Water autoclave: initial fill that becomes liquid-full at 200 °C', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 200, fill: 60 }, kpi: 'phiCrit', expected: 864.66 / 997.05, tol: 0.005, kind: 'reference-data',
    reference: 'Steam tables: saturated-liquid density of water 864.66 kg/m³ at 200 °C and 997.05 kg/m³ at 25 °C'
  },
  {
    id: 's06-ethanol-boiling-point', suite: 6, title: 'Normal boiling point of ethanol', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'ethanol', T: 60, Pamb: 101325 }, kpi: 'Tb', expected: 78.3, abs: 0.6, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: normal boiling point of ethanol 78.3 °C'
  },
  {
    id: 's06-methanol-boiling-point', suite: 6, title: 'Normal boiling point of methanol', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'methanol', T: 50, Pamb: 101325 }, kpi: 'Tb', expected: 64.7, abs: 0.6, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: normal boiling point of methanol 64.7 °C'
  },
  {
    id: 's06-toluene-boiling-point', suite: 6, title: 'Normal boiling point of toluene', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'toluene', T: 80, Pamb: 101325 }, kpi: 'Tb', expected: 110.6, abs: 0.6, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: normal boiling point of toluene 110.6 °C'
  },
  {
    id: 's06-dmf-boiling-point', suite: 6, title: 'Normal boiling point of N,N-dimethylformamide', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'dmf', T: 100, Pamb: 101325 }, kpi: 'Tb', expected: 153, abs: 1.0, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: normal boiling point of DMF 153 °C'
  },
  {
    id: 's06-eg-boiling-point', suite: 6, title: 'Normal boiling point of ethylene glycol', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'eg', T: 150, Pamb: 101325 }, kpi: 'Tb', expected: 197.3, abs: 1.0, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: normal boiling point of ethylene glycol 197.3 °C'
  },
  // ------------------------------------------------------------------ aqueous speciation
  {
    id: 's06-zno-minimum-solubility', suite: 6, title: 'ZnO (zincite): minimum of the amphoteric solubility curve', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Zn', T: 25, I: 0, pHsrc: 'given', pH: 9.8 }, kpi: 'S_min', expected: znMin * 1000, tol: 0.01, kind: 'analytical',
    reference: 'Baes & Mesmer, The Hydrolysis of Cations — mononuclear Zn(II) hydrolysis constants and ZnO solubility constant; mass-action sum evaluated independently (≈ 1.8×10⁻⁶ mol/L)',
    note: 'Infinite dilution (I = 0) so that activity corrections vanish; the suite scans pH in steps of 0.05.'
  },
  {
    id: 's06-zno-ph-of-minimum', suite: 6, title: 'ZnO (zincite): pH of minimum solubility', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Zn', T: 25, I: 0, pHsrc: 'given', pH: 9.8 }, kpi: 'pH_min', expected: znMinPH, abs: 0.05, kind: 'analytical',
    reference: 'Baes & Mesmer, The Hydrolysis of Cations — minimum of Σ[Zn(OH)ₙ] at saturation, located independently',
    note: 'Absolute tolerance equals the pH step of the suite scan.'
  },
  {
    id: 's06-zno-solubility-ph7', suite: 6, title: 'ZnO solubility in the acid branch (pH 7)', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Zn', T: 25, I: 0, pHsrc: 'given', pH: 7 }, kpi: 'S_eq', expected: znoSol(7) * 1000, tol: 0.005, kind: 'analytical',
    reference: 'Baes & Mesmer, The Hydrolysis of Cations — log *Ks0 = 11.14 for ZnO + 2H⁺ ⇌ Zn²⁺ + H₂O; S ≈ 10^(11.14 − 2pH) plus hydroxo species'
  },
  {
    id: 's06-brucite-saturation-index', suite: 6, title: 'Brucite: saturation index of 0.1 mol/L Mg²⁺ at pH 10', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Mg', T: 25, I: 0, pHsrc: 'given', pH: 10, cM: 100 }, kpi: 'SI',
    expected: Math.log10(0.1) - (16.84 - 20 + Math.log10(1 + Math.pow(10, -11.44 + 10))), abs: 0.005, kind: 'analytical',
    reference: 'Baes & Mesmer, The Hydrolysis of Cations — Mg(OH)₂ + 2H⁺ ⇌ Mg²⁺ + 2H₂O, log *Ks0 = 16.84 (Ksp ≈ 10⁻¹¹·¹⁶); SI = log(C_T/S_eq)'
  },
  {
    id: 's06-davies-divalent', suite: 6, title: 'Davies activity coefficient of a divalent ion at I = 0.05 mol/L', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Zn', T: 25, I: 0.05 }, kpi: 'gammaM', expected: davies(2, 0.05), tol: 0.01, kind: 'analytical',
    reference: 'Davies equation, log γ = −A z²(√I/(1+√I) − 0.3 I) with A = 0.509 for water at 25 °C (Stumm & Morgan, Aquatic Chemistry)',
    note: 'The suite evaluates A from its own permittivity correlation (A ≈ 0.512), a 0.4 % difference in γ.'
  },
  {
    id: 's06-davies-trivalent', suite: 6, title: 'Davies activity coefficient of a trivalent ion at I = 0.1 mol/L', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Fe', T: 25, I: 0.1, pH: 3 }, kpi: 'gammaM', expected: davies(3, 0.1), tol: 0.02, kind: 'analytical',
    reference: 'Davies equation with A = 0.509 for water at 25 °C (Stumm & Morgan, Aquatic Chemistry)',
    note: 'z² = 9 amplifies the small difference in the Debye–Hückel constant.'
  },
  {
    id: 's06-strong-base-ph', suite: 6, title: 'pH of 0.01 mol/L strong base from the charge balance', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Mg', T: 25, I: 0, pHsrc: 'base', cBase: 10, cM: 1e-4 }, kpi: 'pH', expected: 12, abs: 0.005, kind: 'analytical',
    reference: 'pH = pKw + log[OH⁻] for a fully dissociated base at 25 °C (pKw = 14.00); metal at trace level (10⁻⁷ mol/L)'
  },
  {
    id: 's06-pkw-60C', suite: 6, title: 'Ion product of water at 60 °C', mode: 'speciation',
    inputs: { mode: 'speciation', metal: 'Zn', T: 60 }, kpi: 'pKw', expected: 13.02, abs: 0.08, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics: pKw of water = 13.02 at 60 °C',
    note: 'The suite uses a constant-enthalpy van\'t Hoff shift (55.8 kJ/mol), which drifts by a few hundredths of a unit per 35 K.'
  },
  // ------------------------------------------------------------------ reaction equilibrium
  {
    id: 's06-vant-hoff-logK', suite: 6, title: 'Equilibrium constant from ΔG° = ΔH° − TΔS°', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: 350, dH: -50, dS: -100, dCp: 0 }, kpi: 'log10K', expected: -dG350 / (R * 350 * LN10), tol: 1e-4, kind: 'analytical',
    reference: 'ΔG° = −RT ln K (Atkins, Physical Chemistry); ΔG°(350 K) = −15 kJ/mol'
  },
  {
    id: 's06-inversion-temperature', suite: 6, title: 'Inversion temperature ΔH°/ΔS°', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: 350, dH: -50, dS: -100, dCp: 0 }, kpi: 'Tinv', expected: 500, tol: 1e-6, kind: 'analytical',
    reference: 'ΔG° = 0 at T = ΔH°/ΔS° (Atkins, Physical Chemistry)'
  },
  {
    id: 's06-vant-hoff-slope', suite: 6, title: 'van\'t Hoff slope d ln K/dT', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: 350, dH: -50, dS: -100, dCp: 0 }, kpi: 'dlnK', expected: -50e3 / (R * 350 * 350), tol: 1e-4, kind: 'analytical',
    reference: 'van\'t Hoff equation d ln K/dT = ΔH°/(RT²) (Atkins, Physical Chemistry)'
  },
  {
    id: 's06-kirchhoff-dG', suite: 6, title: 'ΔG° at 600 K with a heat-capacity (Kirchhoff) correction', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: Tk, dH: -40, dS: -60, dCp }, kpi: 'dG', expected: dGkirch / 1000, tol: 1e-4, kind: 'analytical',
    reference: 'Kirchhoff\'s law: ΔH°(T) = ΔH°₂₉₈ + ΔCp(T − 298.15), ΔS°(T) = ΔS°₂₉₈ + ΔCp ln(T/298.15) (Atkins, Physical Chemistry)'
  },
  {
    id: 's06-conversion-K2-equimolar', suite: 6, title: 'A + B ⇌ C with K = 2 and 1 mol/L equimolar feed: conversion 0.5', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: 400, dH: 0, dS: R * Math.LN2, dCp: 0, bst: 1, ratioB: 1, cM: 1000 }, kpi: 'X', expected: 0.5, tol: 1e-4, kind: 'analytical',
    reference: 'Mass action K = X/[(1 − X)²c_A0/c°]; K = 2, c_A0 = c° gives X = 1/2 exactly'
  },
  {
    id: 's06-conversion-isomerisation', suite: 6, title: 'Isomerisation A ⇌ C: X = K/(1 + K)', mode: 'reaction',
    inputs: { mode: 'reaction', Tr: 400, dH: 0, dS: R * Math.log(3), dCp: 0, bst: 0, ratioB: 0 }, kpi: 'X', expected: 0.75, tol: 1e-4, kind: 'analytical',
    reference: 'Mass action for A ⇌ C with K = 3 (ΔS° = R ln 3, ΔH° = 0): X = K/(1 + K) = 0.75'
  },
  // ------------------------------------------------------------------ regular solution
  {
    id: 's06-regular-critical-temperature', suite: 6, title: 'Regular solution: consolute temperature Ω/2R', mode: 'regular',
    inputs: { mode: 'regular', omega: 18, Tr: 900, x0: 0.4 }, kpi: 'Tc', expected: OMEGA / (2 * R), tol: 1e-5, kind: 'analytical',
    reference: 'Regular-solution model, T_c = Ω/2R (Porter & Easterling, Phase Transformations in Metals and Alloys)'
  },
  {
    id: 's06-regular-binodal-quarter', suite: 6, title: 'Regular solution: binodal at x = 0.25 when Ω/RT = 2 ln 3', mode: 'regular',
    inputs: { mode: 'regular', omega: 18, Tr: T_x25, x0: 0.5 }, kpi: 'xa', expected: 0.25, tol: 1e-4, kind: 'analytical',
    reference: 'Common-tangent condition of the symmetric regular solution, ln[x/(1 − x)] = (Ω/RT)(2x − 1), solved exactly by x = 1/4 at Ω/RT = 2 ln 3'
  },
  {
    id: 's06-regular-spinodal', suite: 6, title: 'Regular solution: spinodal composition at T = 0.75 T_c', mode: 'regular',
    inputs: { mode: 'regular', omega: 18, Tr: 0.75 * OMEGA / (2 * R), x0: 0.5 }, kpi: 's1', expected: 0.25, tol: 1e-5, kind: 'analytical',
    reference: 'Spinodal of the regular solution, 4x(1 − x) = T/T_c, i.e. x = 1/4 at T = 0.75 T_c (Porter & Easterling)'
  },
  {
    id: 's06-regular-lever-rule', suite: 6, title: 'Regular solution: lever rule inside the miscibility gap', mode: 'regular',
    inputs: { mode: 'regular', omega: 18, Tr: T_x25, x0: 0.4 }, pub: 'thermo.phase_fraction', expected: (0.4 - 0.25) / 0.5, tol: 1e-3, kind: 'conservation',
    reference: 'Lever rule (mass balance) between binodal compositions 0.25 and 0.75'
  },
  {
    id: 's06-regular-activity-coefficient', suite: 6, title: 'Regular solution: activity coefficient γ_B = exp[Ω(1 − x)²/RT]', mode: 'regular',
    inputs: { mode: 'regular', omega: 10, Tr: 1000, x0: 0.2 }, kpi: 'gammaB', expected: Math.exp(10e3 * 0.64 / (R * 1000)), tol: 1e-5, kind: 'analytical',
    reference: 'Regular-solution activity coefficient, RT ln γ_B = Ω(1 − x_B)² (Porter & Easterling)'
  },
  // ------------------------------------------------------------------ nanoscale thermodynamics
  {
    id: 's06-gibbs-thomson-au-5nm', suite: 6, title: 'Gibbs–Thomson melting point of a 5 nm gold particle', mode: 'nano',
    inputs: { mode: 'nano', gamma: 0.27, Mw: 196.97, rho: 19300, dnm: 5, Tm: 1337.3, Hf: 12.55 }, kpi: 'Tm_d',
    expected: au.Tm * (1 - 4 * au.g * au.Vm / (au.Hf * 5e-9)), tol: 1e-4, kind: 'analytical',
    reference: 'Gibbs–Thomson relation T_m(d) = T_m[1 − 4γV_m/(ΔH_f d)] with handbook Au data (T_m 1337.3 K, ΔH_f 12.55 kJ/mol, ρ 19 300 kg/m³)',
    note: 'About 1100 K, consistent in magnitude with the measured melting-point depression of nanometre gold particles (Buffat & Borel).'
  },
  {
    id: 's06-capillary-length-zno', suite: 6, title: 'Capillary length 2γV_m/RT of ZnO at 25 °C', mode: 'nano',
    inputs: { mode: 'nano', T: 25, gamma: 0.25, Mw: 81.38, rho: 5606 }, kpi: 'lcap', expected: lcapZnO * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Ostwald–Freundlich (Kelvin) relation; capillary length 2γV_m/RT'
  },
  {
    id: 's06-ostwald-freundlich', suite: 6, title: 'Ostwald–Freundlich solubility enhancement of a 10 nm ZnO particle', mode: 'nano',
    inputs: { mode: 'nano', T: 25, gamma: 0.25, Mw: 81.38, rho: 5606, dnm: 10 }, kpi: 'sol_enh', expected: Math.exp(lcapZnO / 5e-9), tol: 1e-4, kind: 'analytical',
    reference: 'Ostwald–Freundlich equation c(r)/c∞ = exp[2γV_m/(rRT)]'
  },
  {
    id: 's06-critical-radius', suite: 6, title: 'Classical nucleation theory: critical radius at S = 50', mode: 'nano',
    inputs: { mode: 'nano', T: 25, gamma: 0.25, Mw: 81.38, rho: 5606, S }, kpi: 'r_crit', expected: lcapZnO / Math.log(S) * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Classical nucleation theory, r* = 2γV_m/(RT ln S) (Mullin, Crystallization)'
  },
  {
    id: 's06-nucleation-barrier', suite: 6, title: 'Classical nucleation theory: homogeneous barrier at S = 50', mode: 'nano',
    inputs: { mode: 'nano', T: 25, gamma: 0.25, Mw: 81.38, rho: 5606, S }, kpi: 'barrier', expected: barZnO, tol: 1e-4, kind: 'analytical',
    reference: 'Classical nucleation theory, ΔG* = 16πγ³v²/[3(k_BT ln S)²] (Mullin, Crystallization)'
  },
  {
    id: 's06-heterogeneous-factor', suite: 6, title: 'Heterogeneous nucleation: barrier ratio at a 60° contact angle', mode: 'nano',
    inputs: { mode: 'nano', T: 25, gamma: 0.25, Mw: 81.38, rho: 5606, S, theta: 60 }, kpi: 'barrier_het', expected: barZnO * 0.15625, tol: 1e-4, kind: 'analytical',
    reference: 'Volmer wetting factor f(θ) = (2 + cos θ)(1 − cos θ)²/4 = 5/32 at θ = 60°'
  },
  {
    id: 's06-driving-force-eV', suite: 6, title: 'Driving force k_BT ln S per formula unit', mode: 'nano',
    inputs: { mode: 'nano', T: 25, S }, kpi: 'dmu', expected: KB * Tn * Math.log(S) / QE, tol: 1e-4, kind: 'analytical',
    reference: 'Δμ = k_BT ln S with CODATA constants (k_BT = 25.69 meV at 298.15 K)'
  },
  // ------------------------------------------------------------------ published reference measurements (NIST Chemistry WebBook)
  {
    id: 's06-exp-water-psat-150C', suite: 6, title: 'Water saturation pressure at 150 °C (NIST WebBook, IAPWS-95)', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 150 }, kpi: 'Psat', expected: 0.47616, tol: 0.015, kind: 'experiment',
    measured: { value: 0.47616, uncertainty: 0.00012, unit: 'MPa' },
    source: { authors: 'Wagner, W.; Pruss, A.', title: 'The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use', journal: 'J. Phys. Chem. Ref. Data', year: 2002, doi: '10.1063/1.1461829', url: 'https://webbook.nist.gov/chemistry/fluid/' },
    conditions: 'Pure water on the liquid–vapour saturation line at 150 °C',
    match: ['solvent', 'T'],
    reference: 'Wagner & Pruss, J. Phys. Chem. Ref. Data 31, 387–535 (2002), as served by the NIST Chemistry WebBook (Thermophysical Properties of Fluid Systems)',
    note: 'Saturation-properties table of the NIST Chemistry WebBook fluid database (temperature increments, 150.00 °C: 0.47616 MPa). This is critically evaluated reference data: an equation of state fitted to the body of vapour-pressure measurements, with a stated uncertainty of 0.025 % in vapour pressure. The suite uses a two-range Antoine correlation; the tolerance is its documented accuracy.'
  },
  {
    id: 's06-exp-water-psat-250C', suite: 6, title: 'Water saturation pressure at 250 °C (NIST WebBook, IAPWS-95)', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'water', T: 250 }, kpi: 'Psat', expected: 3.9762, tol: 0.015, kind: 'experiment',
    measured: { value: 3.9762, uncertainty: 0.001, unit: 'MPa' },
    source: { authors: 'Wagner, W.; Pruss, A.', title: 'The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use', journal: 'J. Phys. Chem. Ref. Data', year: 2002, doi: '10.1063/1.1461829', url: 'https://webbook.nist.gov/chemistry/fluid/' },
    conditions: 'Pure water on the liquid–vapour saturation line at 250 °C (typical upper range of hydrothermal autoclaves)',
    match: ['solvent', 'T'],
    reference: 'Wagner & Pruss, J. Phys. Chem. Ref. Data 31, 387–535 (2002), as served by the NIST Chemistry WebBook (Thermophysical Properties of Fluid Systems)',
    note: 'Saturation-properties table of the NIST Chemistry WebBook fluid database (250.00 °C: 3.9762 MPa); critically evaluated reference data with a stated vapour-pressure uncertainty of 0.025 %. The Antoine correlation of the suite is about 0.8 % high here.'
  },
  {
    id: 's06-exp-methanol-psat-150C', suite: 6, title: 'Methanol saturation pressure at 150 °C (NIST WebBook)', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'methanol', T: 150 }, kpi: 'Psat', expected: 1.3933, tol: 0.015, kind: 'experiment',
    measured: { value: 1.3933, uncertainty: null, unit: 'MPa' },
    source: { authors: 'de Reuck, K. M.; Craven, R. J. B.', title: 'Methanol, International Thermodynamic Tables of the Fluid State – 12', journal: 'IUPAC, Blackwell Scientific Publications, London', year: 1993, doi: null, url: 'https://webbook.nist.gov/chemistry/fluid/' },
    conditions: 'Pure methanol on the liquid–vapour saturation line at 150 °C (solvothermal conditions)',
    match: ['solvent', 'T'],
    reference: 'de Reuck & Craven, International Thermodynamic Tables of the Fluid State – 12: Methanol (IUPAC, 1993), as served by the NIST Chemistry WebBook',
    note: 'Saturation-properties table of the NIST Chemistry WebBook fluid database (150.00 °C: 1.3933 MPa); critically evaluated reference data from an equation of state fitted to measurements. The WebBook states the equation-of-state uncertainty as generally 0.1 % in density and gives no separate figure for vapour pressure.'
  },
  {
    id: 's06-exp-toluene-psat-150C', suite: 6, title: 'Toluene saturation pressure at 150 °C (NIST WebBook)', mode: 'vapour',
    inputs: { mode: 'vapour', solvent: 'toluene', T: 150 }, kpi: 'Psat', expected: 0.27536, tol: 0.015, kind: 'experiment',
    measured: { value: 0.27536, uncertainty: 0.00014, unit: 'MPa' },
    source: { authors: 'Lemmon, E. W.; Span, R.', title: 'Short Fundamental Equations of State for 20 Industrial Fluids', journal: 'J. Chem. Eng. Data', year: 2006, doi: '10.1021/je050186n', url: 'https://webbook.nist.gov/chemistry/fluid/' },
    conditions: 'Pure toluene on the liquid–vapour saturation line at 150 °C',
    match: ['solvent', 'T'],
    reference: 'Lemmon & Span, J. Chem. Eng. Data 51, 785–850 (2006), as served by the NIST Chemistry WebBook (Thermophysical Properties of Fluid Systems)',
    note: 'Saturation-properties table of the NIST Chemistry WebBook fluid database (150.00 °C: 0.27536 MPa); critically evaluated reference data with a stated vapour-pressure uncertainty of 0.05 % between 305 and 425 K.'
  },
];
