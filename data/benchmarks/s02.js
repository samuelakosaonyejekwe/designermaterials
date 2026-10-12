// Reference benchmarks for suite 2 — Chemical Composition and Material Selection.
// Expected values are evaluated here from standard atomic masses, tabulated radii and closed-form relations.

// IUPAC abridged standard atomic masses (g/mol)
const A = { H: 1.008, N: 14.007, O: 15.999, Na: 22.990, Al: 26.982, Si: 28.085, Zn: 65.38, Ge: 72.630, Se: 78.971, Cd: 112.414 };
const M_ZnO = A.Zn + A.O, M_NO3 = A.N + 3 * A.O, M_H2O = 2 * A.H + A.O, M_ZnNit = A.Zn + 2 * M_NO3, M_ZnNit6 = M_ZnNit + 6 * M_H2O, M_NaOH = A.Na + A.O + A.H;
const N_A = 6.02214076e23, R = 8.314462618, FARADAY = 96485.33212;
// Shannon effective ionic radii (Å): A site twelve-fold, B site and oxide ion six-fold
const tGold = (rA, rB, rX = 1.40) => (rA + rX) / (Math.SQRT2 * (rB + rX));
// Goldschmidt twelve-coordinate metallic radii (Å) of the Cantor-alloy elements
const R12 = [1.251, 1.249, 1.241, 1.350, 1.246], rBar = R12.reduce((a, b) => a + b, 0) / 5, deltaCantor = 100 * Math.sqrt(R12.reduce((s, r) => s + 0.2 * (1 - r / rBar) ** 2, 0));
const RXN_NAOH = 'Zn(NO3)2(aq) + NaOH(aq) -> ZnO(s) + NaNO3(aq) + H2O(l)', RXN_DEC = 'Zn(NO3)2·6H2O(s) -> ZnO(s) + NO2(g) + O2(g) + H2O(g)';
const xAl = 0.02, M_AZO = M_ZnO + xAl * (A.Al - A.Zn), aSiGe = 0.5 * 5.431 + 0.5 * 5.658;

export default [
  {
    id: 's02-molar-mass-zno', suite: 2, title: 'Molar mass of zinc oxide', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO' }, kpi: 'M', expected: 81.38, tol: 2e-4, kind: 'reference-data',
    reference: 'IUPAC standard atomic masses: Zn 65.38, O 15.999 g/mol',
  },
  {
    id: 's02-molar-mass-zinc-nitrate-hexahydrate', suite: 2, title: 'Molar mass of zinc nitrate hexahydrate (hydrate notation)', mode: 'formula',
    inputs: { mode: 'formula', formula: 'Zn(NO3)2·6H2O' }, kpi: 'M', expected: M_ZnNit6, tol: 2e-4, kind: 'reference-data',
    reference: 'IUPAC standard atomic masses; Zn(NO₃)₂·6H₂O is tabulated as 297.5 g/mol in reagent catalogues',
  },
  {
    id: 's02-zno-yield-from-nitrate-hexahydrate', suite: 2, title: 'ZnO mass obtained per mass of Zn(NO₃)₂·6H₂O (thermal decomposition)', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_DEC, targetMass: 10, yieldPct: 100, excess: 0, purity: 100 }, kpi: 'AE', expected: M_ZnO / M_ZnNit6, tol: 2e-4, kind: 'analytical',
    reference: 'Stoichiometry of 2 Zn(NO₃)₂·6H₂O → 2 ZnO + 4 NO₂ + O₂ + 12 H₂O with IUPAC atomic masses',
    note: 'One formula unit of precursor yields one of ZnO, so the atom economy equals the ZnO-to-precursor mass ratio 0.2736.',
  },
  {
    id: 's02-pmi-nitrate-hexahydrate', suite: 2, title: 'Precursor mass per gram of ZnO at full yield (process mass intensity)', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_DEC, targetMass: 10, yieldPct: 100, excess: 0, purity: 100 }, kpi: 'pmi', expected: M_ZnNit6 / M_ZnO, tol: 2e-4, kind: 'analytical',
    reference: 'Stoichiometric mass ratio M(Zn(NO₃)₂·6H₂O)/M(ZnO) = 3.655',
  },
  {
    id: 's02-atom-economy-naoh-route', suite: 2, title: 'Atom economy of Zn(NO₃)₂ + 2 NaOH → ZnO + 2 NaNO₃ + H₂O', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_NAOH, targetMass: 10, yieldPct: 100, excess: 0 }, kpi: 'AE', expected: M_ZnO / (M_ZnNit + 2 * M_NaOH), tol: 2e-4, kind: 'analytical',
    reference: 'Atom economy (Trost): product formula mass over the sum of reactant formula masses of the balanced equation',
  },
  {
    id: 's02-e-factor-stoichiometric', suite: 2, title: 'Reagent E-factor at full yield without excess equals 1/AE − 1', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_NAOH, targetMass: 10, yieldPct: 100, excess: 0 }, kpi: 'Ef', expected: (M_ZnNit + 2 * M_NaOH) / M_ZnO - 1, tol: 2e-4, kind: 'conservation',
    reference: 'Mass balance of the balanced equation: waste = reactant mass − product mass (Sheldon E-factor)',
  },
  {
    id: 's02-theoretical-yield', suite: 2, title: 'Theoretical yield needed for 10 g isolated product at 85 % yield', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_NAOH, targetMass: 10, yieldPct: 85, excess: 10 }, kpi: 'm_th', expected: 10 / 0.85, tol: 1e-9, kind: 'conservation',
    reference: 'Definition of percentage yield: isolated mass = yield × theoretical mass',
  },
  {
    id: 's02-moles-of-product', suite: 2, title: 'Amount of ZnO corresponding to 10 g at full yield', mode: 'formula',
    inputs: { mode: 'formula', formula: 'ZnO', rxn: RXN_NAOH, targetMass: 10, yieldPct: 100, excess: 0 }, kpi: 'n_prod', expected: 10 / M_ZnO, tol: 2e-4, kind: 'analytical',
    reference: 'n = m/M with M(ZnO) = 81.38 g/mol',
  },
  {
    id: 's02-faraday-charge-zinc', suite: 2, title: 'Charge to deposit one mole of zinc (two-electron half-cell)', mode: 'formula',
    inputs: { mode: 'formula', formula: 'Zn', rxn: 'Zn2+(aq) + e- -> Zn(s)', targetMass: A.Zn, yieldPct: 100, excess: 0 }, kpi: 'Q', expected: 2 * FARADAY, tol: 2e-4, kind: 'analytical',
    reference: 'Faraday’s law of electrolysis: Q = zFn with F = 96 485 C/mol',
  },
  {
    id: 's02-goldschmidt-srtio3', suite: 2, title: 'Goldschmidt tolerance factor of SrTiO₃', mode: 'doping',
    inputs: { mode: 'doping', abx: 'SrTiO3' }, kpi: 't', expected: 1.00, tol: 0.005, kind: 'reference-data',
    reference: 'Goldschmidt tolerance factor with Shannon radii: SrTiO₃ is the textbook ideal cubic perovskite, t ≈ 1.00',
  },
  {
    id: 's02-goldschmidt-batio3', suite: 2, title: 'Goldschmidt tolerance factor of BaTiO₃', mode: 'doping',
    inputs: { mode: 'doping', abx: 'BaTiO3' }, kpi: 't', expected: tGold(1.61, 0.605), tol: 1e-6, kind: 'analytical',
    reference: 'Goldschmidt t = (r_A + r_O)/(√2 (r_B + r_O)) with Shannon radii Ba²⁺(XII) 1.61, Ti⁴⁺(VI) 0.605, O²⁻ 1.40 Å; t ≈ 1.06',
  },
  {
    id: 's02-goldschmidt-catio3', suite: 2, title: 'Goldschmidt tolerance factor of CaTiO₃ (tilted perovskite)', mode: 'doping',
    inputs: { mode: 'doping', abx: 'CaTiO3' }, kpi: 't', expected: tGold(1.34, 0.605), tol: 1e-6, kind: 'analytical',
    reference: 'Goldschmidt t with Shannon radii Ca²⁺(XII) 1.34, Ti⁴⁺(VI) 0.605, O²⁻ 1.40 Å; t ≈ 0.97',
  },
  {
    id: 's02-doping-atomic-percent', suite: 2, title: '2 % Al on the Zn site of ZnO expressed as at % of all atoms', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'ZnO', hostSite: 'Zn', dopant: 'Al', dopConc: 2, dopUnit: 'site' }, kpi: 'at', expected: 1.0, tol: 1e-9, kind: 'analytical',
    reference: 'Composition arithmetic: Zn₀.₉₈Al₀.₀₂O contains 0.02 Al in 2 atoms per formula unit',
  },
  {
    id: 's02-doping-weight-percent', suite: 2, title: '2 % Al on the Zn site of ZnO expressed as wt % Al', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'ZnO', hostSite: 'Zn', dopant: 'Al', dopConc: 2, dopUnit: 'site' }, kpi: 'wt', expected: 100 * xAl * A.Al / M_AZO, tol: 2e-4, kind: 'analytical',
    reference: 'w = x·A(Al)/M(Zn₁₋ₓAlₓO) with IUPAC atomic masses',
  },
  {
    id: 's02-doping-number-density', suite: 2, title: 'Al atoms per cm³ for 2 % site doping of ZnO (no lattice change)', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'ZnO', hostSite: 'Zn', dopant: 'Al', dopConc: 2, dopUnit: 'site', hostRho: 5606, aHost: 3.25, aEnd: 3.25 }, kpi: 'ccm', expected: xAl * 5606 * N_A / (M_ZnO * 1e-3) * 1e-6, tol: 2e-4, kind: 'analytical',
    reference: 'Zn-site density of ZnO, ρN_A/M = 4.15 × 10²² cm⁻³ (5606 kg/m³, 81.38 g/mol), times the site fraction',
  },
  {
    id: 's02-doping-cm3-round-trip', suite: 2, title: 'Dopant density entered in cm⁻³ is returned unchanged', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'ZnO', hostSite: 'Zn', dopant: 'Al', dopConc: 1e21, dopUnit: 'cm3' }, kpi: 'ccm', expected: 1e21, tol: 1e-6, kind: 'conservation',
    reference: 'Self-consistency of concentration units: converting cm⁻³ → site fraction → cm⁻³ must be the identity',
  },
  {
    id: 's02-vegard-sige-lattice', suite: 2, title: 'Vegard lattice parameter of Si₀.₅Ge₀.₅', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'Si', hostSite: 'Si', dopant: 'Ge', dopConc: 50, dopUnit: 'site', hostVal: 4, dopVal: 4, hostRho: 2329, aHost: 5.431, aEnd: 5.658 }, kpi: 'a', expected: aSiGe, tol: 1e-9, kind: 'analytical',
    reference: 'Vegard’s law a(x) = (1 − x)a_Si + x·a_Ge with a_Si = 5.431 Å and a_Ge = 5.658 Å',
  },
  {
    id: 's02-vegard-sige-density', suite: 2, title: 'Density of Si₀.₅Ge₀.₅ from Vegard scaling', mode: 'doping',
    inputs: { mode: 'doping', hostF: 'Si', hostSite: 'Si', dopant: 'Ge', dopConc: 50, dopUnit: 'site', hostVal: 4, dopVal: 4, hostRho: 2329, aHost: 5.431, aEnd: 5.658 }, kpi: 'rho', expected: 8 * (0.5 * A.Si + 0.5 * A.Ge) * 1e-3 / (N_A * (aSiGe * 1e-10) ** 3), tol: 0.002, kind: 'analytical',
    reference: 'X-ray density of the diamond-cubic alloy, ρ = 8·M̄/(N_A·a³), with the Vegard lattice parameter',
    note: 'The measured density of Si₀.₅Ge₀.₅ (about 3950 kg/m³) lies 0.7 % higher because the real lattice parameter bows below Vegard’s line.',
  },
  {
    id: 's02-entropy-equimolar-five', suite: 2, title: 'Configurational entropy of an equimolar five-component alloy', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Co:1, Cr:1, Fe:1, Mn:1, Ni:1' }, kpi: 'S', expected: Math.log(5), tol: 1e-9, kind: 'analytical',
    reference: 'Ideal mixing entropy ΔS = −RΣxᵢ ln xᵢ = R ln 5 = 1.609 R (Boltzmann; high-entropy-alloy definition)',
  },
  {
    id: 's02-entropy-non-equimolar', suite: 2, title: 'Configurational entropy of a 3:1 binary', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Cu:3, Au:1' }, kpi: 'S', expected: -(0.75 * Math.log(0.75) + 0.25 * Math.log(0.25)), tol: 1e-9, kind: 'analytical',
    reference: 'Ideal mixing entropy ΔS = −RΣxᵢ ln xᵢ',
  },
  {
    id: 's02-regular-solution-enthalpy', suite: 2, title: 'Regular-solution mixing enthalpy of an equimolar quinary with Ω = −10 kJ/mol', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Co:1, Cr:1, Fe:1, Mn:1, Ni:1', omega: -10 }, kpi: 'H', expected: -10 * (1 - 5 * 0.04) / 2, tol: 1e-9, kind: 'analytical',
    reference: 'Regular solution: ΔH = ΩΣ_{i<j} xᵢxⱼ = Ω(1 − Σxᵢ²)/2',
  },
  {
    id: 's02-ideal-gibbs-energy', suite: 2, title: 'Ideal Gibbs energy of mixing of the equimolar quinary at 1273 K', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Co:1, Cr:1, Fe:1, Mn:1, Ni:1', Talloy: 1273 }, kpi: 'Gid', expected: -1273 * R * Math.log(5) / 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'ΔG_ideal = −TΔS_mix = −RT ln 5 with R = 8.314 J/mol/K',
  },
  {
    id: 's02-cantor-vec', suite: 2, title: 'Valence electron concentration of CoCrFeMnNi', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Co:1, Cr:1, Fe:1, Mn:1, Ni:1' }, kpi: 'VEC', expected: (9 + 6 + 8 + 7 + 10) / 5, tol: 1e-9, kind: 'reference-data',
    reference: 'Valence electron counts (s + d electrons) Co 9, Cr 6, Fe 8, Mn 7, Ni 10; VEC = 8.0 for the Cantor alloy (fcc)',
  },
  {
    id: 's02-cantor-size-mismatch', suite: 2, title: 'Atomic size mismatch δ of CoCrFeMnNi', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Co:1, Cr:1, Fe:1, Mn:1, Ni:1' }, kpi: 'delta', expected: deltaCantor, tol: 0.005, kind: 'reference-data',
    reference: 'δ = 100·√(Σxᵢ(1 − rᵢ/r̄)²) with Goldschmidt twelve-coordinate metallic radii (metals reference tables); about 3.3 % for the Cantor alloy',
  },
  {
    id: 's02-regular-solution-critical-temperature', suite: 2, title: 'Critical temperature of a symmetric regular solution, Ω = +20 kJ/mol', mode: 'alloy',
    inputs: { mode: 'alloy', alloyText: 'Cu:1, Ni:1', omega: 20 }, kpi: 'Tc', expected: 20000 / (2 * R), tol: 1e-6, kind: 'analytical',
    reference: 'Regular-solution model: miscibility gap closes at T_c = Ω/(2R) (thermodynamics of solutions textbooks)',
  },
  {
    id: 's02-mixture-identity-density', suite: 2, title: 'Mixture of copper with copper returns the density of copper', mode: 'mixture',
    inputs: { mode: 'mixture', m1: 'copper', m2: 'copper', m3: 'none', vf1: 0.37 }, kpi: 'rho', expected: 8960, tol: 0.005, kind: 'reference-data',
    reference: 'Density of copper at room temperature, 8960 kg/m³ (metals handbooks)',
  },
  {
    id: 's02-mixture-identity-conductivity', suite: 2, title: 'Maxwell–Eucken estimate for identical phases returns the conductivity of copper', mode: 'mixture',
    inputs: { mode: 'mixture', m1: 'copper', m2: 'copper', m3: 'none', vf1: 0.37 }, kpi: 'kME', expected: 400, tol: 0.01, kind: 'reference-data',
    reference: 'Thermal conductivity of pure copper at 300 K, about 400 W/m/K (heat-transfer handbooks)',
  },
  {
    id: 's02-porosity-from-density', suite: 2, title: 'Porosity of a copper compact measured at 8064 kg/m³', mode: 'mixture',
    inputs: { mode: 'mixture', m1: 'copper', m2: 'copper', m3: 'none', vf1: 0.5, rhoMeas: 8064 }, kpi: 'por', expected: 1 - 8064 / 8960, abs: 0.003, kind: 'reference-data',
    reference: 'φ = 1 − ρ_measured/ρ_theoretical with the handbook density of copper, 8960 kg/m³',
  },
  {
    id: 's02-maxwell-eucken-conducting-filler', suite: 2, title: 'Epoxy with 30 vol % alumina: Maxwell–Eucken conductivity', mode: 'mixture',
    inputs: { mode: 'mixture', m1: 'epoxy', m2: 'alumina', m3: 'none', vf1: 0.7 }, kpi: 'kME', expected: 0.2 * (1 + 2 * 0.3) / (1 - 0.3), tol: 0.03, kind: 'analytical',
    reference: 'Maxwell (Maxwell–Eucken) model in the limit of highly conducting spheres: k/k_m = (1 + 2φ)/(1 − φ), with k_m = 0.2 W/m/K for cured epoxy',
    note: 'The finite conductivity of alumina (about 30 W/m/K, 150 × the matrix) lowers the exact Maxwell value by 1.6 % relative to this limit.',
  },
  {
    id: 's02-turner-zero-expansion', suite: 2, title: 'Zero-expansion composition of ZrW₂O₈ in copper (modulus-weighted Turner rule)', mode: 'mixture',
    inputs: { mode: 'mixture', m1: 'zrw2o8', m2: 'copper', m3: 'none', vf1: 0.5 }, kpi: 'vzero', expected: 16.5 * 120 / (16.5 * 120 + 8.7 * 88), tol: 0.05, kind: 'reference-data',
    reference: 'Turner model α_eff = ΣαᵢEᵢVᵢ/ΣEᵢVᵢ with handbook values: ZrW₂O₈ α ≈ −8.7 × 10⁻⁶ /K, E ≈ 88 GPa; copper α ≈ 16.5 × 10⁻⁶ /K, E ≈ 120 GPa',
    note: 'Young’s moduli replace the bulk moduli of the original Turner formula because no Poisson ratio is tabulated for ZrW₂O₈; with bulk moduli (73 and 140 GPa) the zero-expansion fraction would be about 0.78.',
  },
  {
    id: 's02-subst-cadmium-mass-fraction', suite: 2, title: 'Mass fraction of the hazard-flagged element in CdSe', mode: 'subst',
    inputs: { mode: 'subst', matId: 'cdse-qd' }, kpi: 'w_haz', expected: A.Cd / (A.Cd + A.Se), tol: 2e-4, kind: 'analytical',
    reference: 'Mass fraction of cadmium in CdSe from IUPAC atomic masses (58.7 wt %)',
  },
  {
    id: 's02-subst-mismatch-strain-zno-on-si', suite: 2, title: 'Thermal-mismatch strain of ZnO on silicon for ΔT = 100 K', mode: 'subst',
    inputs: { mode: 'subst', matId: 'zno-nanoparticles', subId: 'silicon', dT: 100 }, kpi: 'eps', expected: (2.6e-6 - 4.3e-6) * 100, tol: 0.1, kind: 'reference-data',
    reference: 'ε = (α_substrate − α_film)ΔT with handbook expansion coefficients: Si 2.6 × 10⁻⁶ /K, polycrystalline ZnO about 4.3 × 10⁻⁶ /K',
    note: 'Tabulated mean expansion coefficients of ZnO range from 4.0 to 4.8 × 10⁻⁶ /K (a- and c-axis values differ), hence the 10 % tolerance.',
  },
  {
    id: 's02-props-copper-diffusivity', suite: 2, title: 'Thermal diffusivity of copper', mode: 'props',
    inputs: { mode: 'props', matId: 'copper' }, kpi: 'alpha', expected: 1.17e-4, tol: 0.02, kind: 'reference-data',
    reference: 'Thermal diffusivity of pure copper at 300 K, 117 × 10⁻⁶ m²/s (heat-transfer textbook property tables)',
  },
  {
    id: 's02-props-copper-shear-modulus', suite: 2, title: 'Shear modulus of copper from E and ν', mode: 'props',
    inputs: { mode: 'props', matId: 'copper' }, kpi: 'G', expected: 4.6e10, tol: 0.05, kind: 'reference-data',
    reference: 'Shear modulus of polycrystalline copper, 44–48 GPa (metals handbooks)',
  },
  {
    id: 's02-props-silicon-bar-speed', suite: 2, title: 'Bar wave speed of silicon', mode: 'props',
    inputs: { mode: 'props', matId: 'silicon' }, kpi: 'cbar', expected: Math.sqrt(160e9 / 2329), tol: 0.03, kind: 'reference-data',
    reference: 'c = √(E/ρ) with the polycrystalline-average modulus of silicon (160 GPa) and its density 2329 kg/m³',
  },
];
