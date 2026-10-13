// Suite 11 — Electrical, electronic and electrochemical engineering: independent reference benchmarks.
// Expected values are handbook literals or closed-form expressions evaluated here with CODATA constants;
// nothing is taken from the suite's own output.

const R = 8.314462618, F = 96485.33212, KB = 1.380649e-23, QE = 1.602176634e-19, EPS0 = 8.8541878128e-12;
const HBAR = 1.054571817e-34, H = 6.62607015e-34, ME = 9.1093837015e-31, NA = 6.02214076e23;
const T25 = 298.15, T300C = 26.85;                         // 300 K expressed in the suite's °C input
const VT300 = KB * 300 / QE;                               // thermal voltage at 300 K, V
const H2_2M0 = HBAR * HBAR / (2 * ME) / QE * 1e18;         // ħ²/2m₀ = 0.0381 eV·nm²
const E2_4PE0 = QE / (4 * Math.PI * EPS0) * 1e9;           // e²/4πε₀ = 1.44 eV·nm
const A0 = 4 * Math.PI * EPS0 * HBAR * HBAR / (ME * QE * QE) * 1e9;   // Bohr radius, nm

// --- electrodeposition: copper in the pure Tafel regime at exactly 40.8 A/m²
const J_CU = 40.8, J0_CU = 1e-3, ETA_CU = Math.log(J_CU / J0_CU) * R * T25 / (0.5 * F);
const M_CU = 63.546e-3, RHO_CU = 8960;                      // kg/mol, kg/m³

// --- Debye length of a 1:1 electrolyte, λ = √(εRT/2F²I)
const debyeNm = (I_molm3, epsr, T) => Math.sqrt(epsr * EPS0 * R * T / (2 * F * F * I_molm3)) * 1e9;

// --- semiconductor closed forms
const dos = (m, T) => 2 * Math.pow(2 * Math.PI * m * ME * KB * T / (H * H), 1.5);      // effective density of states, m⁻³
const EPS_SI = 11.7 * EPS0;
const VBI_SI = VT300 * Math.log(1e18 * 1e16 / (1.0e10 * 1.0e10));                         // with the accepted nᵢ = 1.0×10¹⁰ cm⁻³
const NEFF_SI = 1e24 * 1e22 / (1e24 + 1e22);
const W_SI = Math.sqrt(2 * EPS_SI * VBI_SI / (QE * NEFF_SI)) * 1e9;                       // nm

// --- first bound state of a finite quantum well with a mass discontinuity (BenDaniel–Duke matching):
// k·tan(kw/2) = (m_w/m_b)·κ,  k = √(2m_wE)/ħ,  κ = √(2m_b(V₀ − E))/ħ
function wellLevel(V0, w, mw, mb) {
  const f = E => { const k = Math.sqrt(mw * E / H2_2M0), ka = Math.sqrt(mb * (V0 - E) / H2_2M0); return k * Math.tan(k * w / 2) - (mw / mb) * ka; };
  let lo = 1e-9, hi = Math.min(V0, H2_2M0 * Math.PI * Math.PI / (mw * w * w)) * (1 - 1e-12);
  for (let i = 0; i < 200; i++) { const mid = 0.5 * (lo + hi); if (f(mid) > 0) hi = mid; else lo = mid; }
  return 0.5 * (lo + hi) * 1000;                             // meV
}

// --- percolation: excluded-volume threshold of spherocylinders (used only as an input to the linear effective-medium check)
const PHIC_SPHERE = 1 - Math.exp(-2.8 / 8);

// --- Esaki diode, pure tunnelling term I = I_p (V/V_p) exp(1 − V/V_p)
const IP_E = 5e-3, VP_E = 0.06, RN_E = Math.E * VP_E / IP_E, RS_E = 1, C_E = 5e-12;

export default [
  // ------------------------------------------------------------------ electrodeposition
  {
    id: 's11-deposit-nernst-silver',
    suite: 11,
    title: 'Nernst potential of Ag⁺/Ag at 0.01 mol/L, 25 °C',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Ag', c_ion: 0.01, T: 25 },
    kpi: 'E_eq',
    expected: 0.7996 + (R * T25 / F) * Math.log(0.01),
    tol: 0.002,
    kind: 'reference-data',
    reference: 'Standard electrode potential E°(Ag⁺/Ag) = +0.7996 V (CRC Handbook of Chemistry and Physics) with the Nernst equation',
    note: '0.7996 − 2 × 0.05916 = 0.6813 V vs SHE; activities taken as concentrations.'
  },
  {
    id: 's11-deposit-nernst-copper',
    suite: 11,
    title: 'Nernst potential of Cu²⁺/Cu at 0.05 mol/L, 25 °C',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Cu', c_ion: 0.05, T: 25 },
    kpi: 'E_eq',
    expected: 0.3419 + (R * T25 / (2 * F)) * Math.log(0.05),
    tol: 0.01,
    kind: 'reference-data',
    reference: 'Standard electrode potential E°(Cu²⁺/Cu) = +0.3419 V (CRC Handbook of Chemistry and Physics) with the Nernst equation',
    note: 'Tabulated E° values for copper range from 0.337 to 0.342 V; the suite uses 0.340 V, a difference of 2 mV.'
  },
  {
    id: 's11-deposit-tafel-slope',
    suite: 11,
    title: 'Cathodic Tafel slope 2.303RT/(αF) for α = 0.5 at 25 °C',
    mode: 'deposit',
    inputs: { mode: 'deposit', alpha_c: 0.5, T: 25 },
    kpi: 'tafel',
    expected: Math.LN10 * R * T25 / (0.5 * F) * 1000,
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Tafel equation (Bard & Faulkner, Electrochemical Methods)',
    note: '118.3 mV per decade.'
  },
  {
    id: 's11-deposit-limiting-current',
    suite: 11,
    title: 'Diffusion-limited current density nFDc/δ (Cu²⁺, 0.05 mol/L, δ = 100 µm)',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Cu', c_ion: 0.05, D_ion: 7e-10, delta: 100 },
    kpi: 'j_lim',
    expected: 2 * F * 7e-10 * 50 / 100e-6,
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Nernst diffusion-layer model of the limiting current (Bard & Faulkner, Electrochemical Methods)'
  },
  {
    id: 's11-deposit-mass-transport-plateau',
    suite: 11,
    title: 'At high overpotential the deposition current saturates at nFDc/δ',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Cu', c_ion: 0.05, D_ion: 7e-10, delta: 100, eta: 1.0, j0: 10, alpha_c: 0.5 },
    kpi: 'j',
    expected: 2 * F * 7e-10 * 50 / 100e-6,
    tol: 1e-5,
    kind: 'analytical',
    reference: 'Nernst diffusion-layer model of the limiting current (Bard & Faulkner, Electrochemical Methods)'
  },
  {
    id: 's11-deposit-tafel-current',
    suite: 11,
    title: 'Tafel regime: j = j₀·exp(αFη/RT) far from the mass-transport limit',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Cu', T: 25, j0: J0_CU, alpha_c: 0.5, eta: ETA_CU, c_ion: 2, D_ion: 1e-8, delta: 1, pH: 4 },
    kpi: 'j',
    expected: J_CU,
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Butler–Volmer equation in the high-overpotential (Tafel) limit (Bard & Faulkner, Electrochemical Methods)',
    note: 'The overpotential is chosen here as (RT/αF)·ln(40.8/j₀) = 0.546 V; the limiting current is five orders of magnitude higher (3.9×10⁶ A/m²), so the mass-transport correction is about 10⁻⁵.'
  },
  {
    id: 's11-deposit-faraday-copper-rate',
    suite: 11,
    title: 'Faraday\'s law: copper at 40.8 A/m² grows at jM/(nFρ) = 1.50 nm/s',
    mode: 'deposit',
    inputs: { mode: 'deposit', couple: 'Cu', T: 25, j0: J0_CU, alpha_c: 0.5, eta: ETA_CU, c_ion: 2, D_ion: 1e-8, delta: 1, pH: 4 },
    kpi: 'rate',
    expected: J_CU * M_CU / (2 * F * RHO_CU) * 1e9,
    tol: 1e-3,
    kind: 'analytical',
    reference: 'Faraday\'s laws of electrolysis with M(Cu) = 63.546 g/mol and ρ(Cu) = 8960 kg/m³ (CRC Handbook of Chemistry and Physics)'
  },
  // ------------------------------------------------------------------ ionic transport
  {
    id: 's11-ionic-nernst-einstein-zno-reference',
    suite: 11,
    title: 'Nernst–Einstein conductivity of the ZnO reference electrolyte (0.01 M Zn²⁺, 0.02 M OH⁻, 25 °C)',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'zno', D_src: 'custom', cond_model: 'ne', c_cat: 10, c_an: 20, D_cat: 7.5e-10, D_an: 5.3e-9, T: 25 },
    kpi: 'sigma',
    expected: F * F / (R * T25) * (4 * 7.5e-10 * 10 + 1 * 5.3e-9 * 20),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Nernst–Einstein relation σ = (F²/RT)Σzᵢ²Dᵢcᵢ (Newman, Electrochemical Systems)',
    note: '0.511 S/m.'
  },
  {
    id: 's11-ionic-kcl-limiting-molar-conductivity',
    suite: 11,
    title: 'KCl at 1 mol/m³ from the tabulated ion diffusivities: limiting molar conductivity 149.85 S·cm²/mol',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'kcl', D_src: 'table', cond_model: 'ne', c_cat: 1, c_an: 1, T: 25 },
    kpi: 'sigma',
    expected: 149.85e-4 * 1,
    tol: 0.005,
    kind: 'reference-data',
    reference: 'Limiting ionic conductivities λ°(K⁺) = 73.5 and λ°(Cl⁻) = 76.35 S·cm²/mol at 25 °C (CRC Handbook of Chemistry and Physics; Robinson & Stokes, Electrolyte Solutions)',
    note: 'Infinite-dilution (Nernst–Einstein) option; the suite loads its tabulated ionic diffusion coefficients of K⁺ and Cl⁻ (1.957 and 2.032 ×10⁻⁹ m²/s) when KCl is selected.'
  },
  {
    id: 's11-ionic-hcl-transference-number',
    suite: 11,
    title: 'Cation transference number of dilute HCl at 25 °C',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'hcl', D_src: 'table', cond_model: 'ne', c_cat: 1, c_an: 1, T: 25 },
    kpi: 't_cat',
    expected: 0.821,
    tol: 0.005,
    kind: 'reference-data',
    reference: 'Limiting transference number of H⁺ in HCl, t₊ = 0.821 at 25 °C (Robinson & Stokes, Electrolyte Solutions)'
  },
  {
    id: 's11-ionic-debye-length-0p1M',
    suite: 11,
    title: 'Debye length of a 0.1 mol/L 1:1 electrolyte in water at 25 °C',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'nacl', c_cat: 100, c_an: 100, eps_r: 78.5, T: 25 },
    kpi: 'debye',
    expected: 0.962,
    tol: 0.003,
    kind: 'reference-data',
    reference: 'Debye length λ_D = 0.304 nm/√I for a 1:1 aqueous electrolyte at 25 °C (Israelachvili, Intermolecular and Surface Forces)'
  },
  {
    id: 's11-ionic-debye-length-closed-form',
    suite: 11,
    title: 'Debye length √(εRT/2F²I) for 0.001 mol/L KCl',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'kcl', c_cat: 1, c_an: 1, eps_r: 78.5, T: 25 },
    kpi: 'debye',
    expected: debyeNm(1, 78.5, T25),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Debye–Hückel theory (Israelachvili, Intermolecular and Surface Forces)',
    note: '9.62 nm.'
  },
  {
    id: 's11-ionic-electromigration-flux-enhancement',
    suite: 11,
    title: 'Steady Nernst–Planck flux to a sink: enhancement Pe/(1 − e^(−Pe)) at Pe = zFEL/RT = 1',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'nacl', c_cat: 10, c_an: 10, T: 25, Lgap: 0.001, Efield: R * T25 / F / 1e-6 },
    kpi: 'enh',
    expected: 1 / (1 - Math.exp(-1)),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Steady one-dimensional Nernst–Planck equation in a uniform field (Newman, Electrochemical Systems)',
    note: 'Monovalent cation, 1 µm path and a field of RT/F per micrometre (25.7 kV/m).'
  },
  // ------------------------------------------------------------------ impedance
  {
    id: 's11-eis-randles-characteristic-frequency',
    suite: 11,
    title: 'Randles circuit with an ideal capacitor: apex frequency 1/(2πR_ctC_dl)',
    mode: 'eis',
    inputs: { mode: 'eis', Rct: 100, Cdl: 20, n_cpe: 1, sigma_w: 0 },
    kpi: 'f_c',
    expected: 1 / (2 * Math.PI * 100 * 20e-6),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Randles equivalent circuit (Bard & Faulkner, Electrochemical Methods)',
    note: '79.58 Hz.'
  },
  {
    id: 's11-eis-randles-dc-limit',
    suite: 11,
    title: 'Randles circuit without diffusion: |Z| → R_s + R_ct at low frequency',
    mode: 'eis',
    inputs: { mode: 'eis', Rs: 10, Rct: 100, Cdl: 20, n_cpe: 1, sigma_w: 0, fmin: 1e-4 },
    kpi: 'Z_low',
    expected: 110,
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Randles equivalent circuit (Bard & Faulkner, Electrochemical Methods)'
  },
  {
    id: 's11-eis-randles-maximum-phase',
    suite: 11,
    title: 'R_s + (R_ct ∥ C): most capacitive phase −atan[R_ct/(2√(R_s(R_s + R_ct)))]',
    mode: 'eis',
    inputs: { mode: 'eis', Rs: 10, Rct: 100, Cdl: 20, n_cpe: 1, sigma_w: 0, fmin: 0.01, fmax: 1e5 },
    kpi: 'phase_min',
    expected: -Math.atan(100 / (2 * Math.sqrt(10 * 110))) * 180 / Math.PI,
    tol: 1e-3,
    kind: 'analytical',
    reference: 'Phase extremum of a resistor in series with a parallel RC element (elementary circuit analysis; Orazem & Tribollet, Electrochemical Impedance Spectroscopy)',
    note: '−56.4°. The suite samples 220 logarithmically spaced frequencies; the extremum is flat, so the sampling error is far below the tolerance.'
  },
  {
    id: 's11-eis-exchange-current-from-rct',
    suite: 11,
    title: 'Exchange current density from the charge-transfer resistance, j₀ = RT/(F R_ct A)',
    mode: 'eis',
    inputs: { mode: 'eis', Rct: 100, area: 1, T: 25 },
    kpi: 'j0',
    expected: R * T25 / (F * 100 * 1e-4),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Linearised Butler–Volmer equation for a one-electron reaction (Bard & Faulkner, Electrochemical Methods)'
  },
  // ------------------------------------------------------------------ semiconductors
  {
    id: 's11-semi-silicon-intrinsic-density',
    suite: 11,
    title: 'Intrinsic carrier density of silicon at 300 K',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'Si', T: T300C },
    kpi: 'ni',
    expected: 1.0e10,
    tol: 0.03,
    kind: 'reference-data',
    reference: 'Accepted intrinsic carrier density of silicon at 300 K, about 1×10¹⁰ cm⁻³ (Sze & Ng, Physics of Semiconductor Devices; Sproul & Green, intrinsic carrier concentration of silicon)',
    note: 'The suite uses the tabulated effective densities of states of silicon (2.8×10¹⁹ and 2.65×10¹⁹ cm⁻³ at 300 K, scaled as T^1.5) with the Varshni gap E_g(T) = 1.17 − 4.73×10⁻⁴ T²/(T + 636) eV, which gives 0.98×10¹⁰ cm⁻³. Published determinations lie between about 0.97 and 1.08×10¹⁰ cm⁻³.'
  },
  {
    id: 's11-semi-gaas-fermi-level',
    suite: 11,
    title: 'n-GaAs, 10¹⁶ cm⁻³ at 300 K: E_c − E_F = kT·ln(N_c/N_D) with N_c from m* = 0.067 m₀',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'GaAs', Nd: 1e16, T: T300C },
    kpi: 'EcEf',
    expected: VT300 * Math.log(dos(0.067, 300) / 1e22),
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Boltzmann statistics with the parabolic-band effective density of states N_c = 2(2πm*kT/h²)^(3/2) (Sze & Ng, Physics of Semiconductor Devices)',
    note: 'N_c = 4.35×10¹⁷ cm⁻³, E_c − E_F = 0.0975 eV.'
  },
  {
    id: 's11-semi-silicon-pn-built-in-potential',
    suite: 11,
    title: 'Silicon abrupt p–n junction, N_A = 10¹⁸, N_D = 10¹⁶ cm⁻³ at 300 K: built-in potential',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'Si', contact: 'pn', Na: 1e18, Nd: 1e16, Vbias: 0, T: T300C },
    kpi: 'Vbi',
    expected: VBI_SI,
    tol: 0.02,
    kind: 'analytical',
    reference: 'V_bi = (kT/q)·ln(N_A N_D/nᵢ²) with nᵢ = 1.0×10¹⁰ cm⁻³ (Sze & Ng, Physics of Semiconductor Devices)',
    note: '0.833 V. The reference uses the accepted nᵢ; the suite value differs by about 2 %, which the logarithm reduces to under 0.2 %.'
  },
  {
    id: 's11-semi-silicon-pn-depletion-width',
    suite: 11,
    title: 'Silicon abrupt p–n junction, N_A = 10¹⁸, N_D = 10¹⁶ cm⁻³ at 300 K: zero-bias depletion width',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'Si', contact: 'pn', Na: 1e18, Nd: 1e16, Vbias: 0, T: T300C },
    kpi: 'W',
    expected: W_SI,
    tol: 0.01,
    kind: 'analytical',
    reference: 'Abrupt-junction depletion approximation W = √(2ε_sV_bi(N_A + N_D)/(qN_AN_D)), ε_r(Si) = 11.7 (Sze & Ng, Physics of Semiconductor Devices)',
    note: 'About 0.33 µm. Uses V_bi from the accepted nᵢ = 1.0×10¹⁰ cm⁻³; the square root halves the small difference in V_bi.'
  },
  {
    id: 's11-semi-silicon-debye-length',
    suite: 11,
    title: 'Extrinsic Debye length of n-Si at 10¹⁶ cm⁻³, 300 K',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'Si', Nd: 1e16, T: T300C },
    kpi: 'LD',
    expected: Math.sqrt(EPS_SI * KB * 300 / (QE * QE * 1e22)) * 1e9,
    tol: 1e-4,
    kind: 'analytical',
    reference: 'L_D = √(ε_skT/q²N_D) (Sze & Ng, Physics of Semiconductor Devices)',
    note: '40.9 nm.'
  },
  {
    id: 's11-semi-brus-cdse-quantum-dot',
    suite: 11,
    title: 'Brus band gap of a CdSe quantum dot of radius 2.5 nm',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'CdSe', R_qd: 2.5 },
    kpi: 'Eg_qd',
    expected: 1.74 + H2_2M0 * Math.PI * Math.PI / (2.5 * 2.5) * (1 / 0.13 + 1 / 0.45) - 1.786 * E2_4PE0 / (10.2 * 2.5),
    tol: 1e-3,
    kind: 'analytical',
    reference: 'Brus effective-mass model, E = E_g + (ħ²π²/2R²)(1/m_e + 1/m_h) − 1.786e²/(4πεR), with CdSe handbook parameters E_g = 1.74 eV, m_e = 0.13, m_h = 0.45, ε_r = 10.2',
    note: '2.24 eV (554 nm); evaluated here from ħ²/2m₀ = 0.0381 eV·nm² and e²/4πε₀ = 1.44 eV·nm.'
  },
  {
    id: 's11-semi-cdse-exciton-bohr-radius',
    suite: 11,
    title: 'Exciton Bohr radius of CdSe, a_B = ε_r(m₀/µ)a₀',
    mode: 'semi',
    inputs: { mode: 'semi', semi_mat: 'CdSe' },
    kpi: 'aB',
    expected: 10.2 * (1 / 0.13 + 1 / 0.45) * A0,
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Hydrogenic (Wannier–Mott) exciton model (Kittel, Introduction to Solid State Physics)',
    note: '5.35 nm, consistent with the 5–6 nm usually quoted for CdSe.'
  },
  // ------------------------------------------------------------------ percolation
  {
    id: 's11-percolation-overlapping-spheres',
    suite: 11,
    title: 'Continuum percolation threshold of spheres (aspect ratio 1)',
    mode: 'percolation',
    inputs: { mode: 'percolation', AR: 1 },
    kpi: 'phic',
    expected: 28.95,
    tol: 0.03,
    kind: 'numerical-reference',
    reference: 'Monte-Carlo continuum percolation threshold of overlapping spheres, φ_c ≈ 0.2895 (Lorenz & Ziff; Torquato, Random Heterogeneous Materials)',
    note: 'The suite uses the excluded-volume invariant ⟨V_ex⟩N_c = 2.8 for spheres, which gives 29.5 vol %.'
  },
  {
    id: 's11-percolation-slender-rod-limit',
    suite: 11,
    title: 'Slender randomly oriented rods (aspect ratio 1000): φ_c ≈ 0.7/AR',
    mode: 'percolation',
    inputs: { mode: 'percolation', AR: 1000 },
    kpi: 'phic',
    expected: 0.07,
    tol: 0.01,
    kind: 'analytical',
    reference: 'Excluded-volume theory of percolation for sticks, ⟨V_ex⟩N_c ≈ 1.4 in the slender-rod limit (Balberg, Anderson, Alexander & Wagner, Physical Review B)',
    note: 'KPI in vol %. End-cap terms of the spherocylinder change the limit by a few tenths of a per cent at this aspect ratio.'
  },
  {
    id: 's11-percolation-linear-effective-medium',
    suite: 11,
    title: 'General effective-medium equation with s = t = 1 and an insulating matrix: σ = σ_f(φ − φ_c)/(1 − φ_c)',
    mode: 'percolation',
    inputs: { mode: 'percolation', AR: 1, phi: 50, sig_f: 1e4, sig_m: 1e-16, t_exp: 1, s_exp: 1 },
    kpi: 'sigma',
    expected: 1e4 * (0.5 - PHIC_SPHERE) / (1 - PHIC_SPHERE),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'McLachlan general effective-medium equation, which reduces to the linear Bruggeman-type law for unit exponents and zero matrix conductivity',
    note: 'φ_c = 1 − exp(−2.8/8) is the excluded-volume threshold for spheres, evaluated here in closed form.'
  },
  // ------------------------------------------------------------------ resonant tunnelling diode
  {
    id: 's11-rtd-resonant-level-gaas-algaas',
    suite: 11,
    title: 'Double-barrier resonance of a 5 nm GaAs well with thick 0.3 eV AlGaAs barriers: finite-well ground state',
    mode: 'rtd',
    inputs: { mode: 'rtd', V0: 0.3, b_nm: 6, w_nm: 5, m_w: 0.067, m_b: 0.092 },
    kpi: 'E1',
    expected: wellLevel(0.3, 5, 0.067, 0.092),
    tol: 0.005,
    kind: 'analytical',
    reference: 'Bound state of a finite square well with BenDaniel–Duke boundary conditions, k·tan(kw/2) = (m_w/m_b)κ (Davies, The Physics of Low-Dimensional Semiconductors)',
    note: 'With 6 nm barriers the quasi-bound level of the open structure differs from the bound state of the isolated well only by terms of order exp(−2κb) ≈ 10⁻⁴. The transcendental equation is solved here by bisection.'
  },
  {
    id: 's11-rtd-resonant-level-deep-wide-well',
    suite: 11,
    title: 'Double-barrier resonance of a 10 nm well with 1 eV barriers and equal masses',
    mode: 'rtd',
    inputs: { mode: 'rtd', V0: 1.0, b_nm: 4, w_nm: 10, m_w: 0.067, m_b: 0.067 },
    kpi: 'E1',
    expected: wellLevel(1.0, 10, 0.067, 0.067),
    tol: 0.005,
    kind: 'analytical',
    reference: 'Bound state of a finite square well (Davies, The Physics of Low-Dimensional Semiconductors)',
    note: 'About 42 meV, against 56 meV for an infinitely deep well of the same width.'
  },
  // ------------------------------------------------------------------ Esaki tunnel diode
  {
    id: 's11-esaki-negative-resistance',
    suite: 11,
    title: 'Tunnel-diode term I_p(V/V_p)exp(1 − V/V_p): minimum negative resistance −e·V_p/I_p at V = 2V_p',
    mode: 'esaki',
    inputs: { mode: 'esaki', Ip: 5, Vp: 0.06, Iv: 1e-4, Vv: 1.0, Is: 1e-16, Vmax: 0.6, T: 25 },
    kpi: 'Rn',
    expected: -RN_E,
    tol: 0.002,
    kind: 'analytical',
    reference: 'Empirical tunnel-diode current law (Sze & Ng, Physics of Semiconductor Devices), differentiated analytically',
    note: '−32.6 Ω. Excess and diffusion currents are made negligible; the suite differentiates the sampled curve numerically (2 mV steps).'
  },
  {
    id: 's11-esaki-peak-current',
    suite: 11,
    title: 'Tunnel diode: peak of the characteristic at (V_p, I_p)',
    mode: 'esaki',
    inputs: { mode: 'esaki', Ip: 5, Vp: 0.06, Iv: 1e-4, Vv: 1.0, Is: 1e-16, Vmax: 0.6, T: 25 },
    kpi: 'Ip',
    expected: 5,
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Empirical tunnel-diode current law (Sze & Ng, Physics of Semiconductor Devices)'
  },
  {
    id: 's11-esaki-resistive-cutoff-frequency',
    suite: 11,
    title: 'Tunnel diode resistive cut-off frequency √(|R_n|/R_s − 1)/(2π|R_n|C)',
    mode: 'esaki',
    inputs: { mode: 'esaki', Ip: 5, Vp: 0.06, Iv: 1e-4, Vv: 1.0, Is: 1e-16, Vmax: 0.6, T: 25, Rs_c: RS_E, Cj: 5 },
    kpi: 'fr',
    expected: Math.sqrt(RN_E / RS_E - 1) / (2 * Math.PI * RN_E * C_E) / 1e9,
    tol: 0.003,
    kind: 'analytical',
    reference: 'Small-signal equivalent circuit of the tunnel diode (Sze & Ng, Physics of Semiconductor Devices)',
    note: 'With |R_n| = e·V_p/I_p = 32.6 Ω, R_s = 1 Ω and C = 5 pF: 5.49 GHz.'
  },
  // ------------------------------------------------------------------ published experiments
  {
    id: 's11-exp-kcl-001-molal-conductivity',
    suite: 11,
    title: 'KCl primary conductivity standard, 0.01 mol/kg at 25 °C: measured electrolytic conductivity',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'kcl', D_src: 'table', cond_model: 'dho', c_cat: 9.968, c_an: 9.968, T: 25 },
    kpi: 'sigma',
    expected: 0.140823,
    tol: 0.02,
    kind: 'experiment',
    measured: { value: 1408.23, uncertainty: 0.42, unit: 'µS/cm' },
    source: { authors: 'Jameel, Wu, Pratt', title: 'Primary Standards and Standard Reference Materials for Electrolytic Conductivity', journal: 'NIST Special Publication 260-142', year: 2000, doi: null, url: 'https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication260-142.pdf' },
    conditions: 'Aqueous KCl, molality 0.01 mol/kg, 25 °C (ITS-90); absolute conductivity measurement; value excludes the conductivity of the solvent; stated uncertainty is the expanded uncertainty 2u_c',
    match: ['system', 'c_cat', 'c_an', 'T'],
    reference: 'Jameel, Wu & Pratt, NIST Special Publication 260-142 (2000), Table 3',
    note: 'Value from Table 3 (molality-based primary standards, 25 °C row). 0.01 mol/kg corresponds to 9.968 mol/m³ with a solution density of 997.5 kg/m³. The suite loads the tabulated limiting diffusivities of K⁺ and Cl⁻ (Λ° = 149.9 S·cm²/mol) and applies the Debye–Hückel–Onsager relaxation and electrophoretic terms with the ion-size factor 1/(1 + κa) at the default a = 0.35 nm, which was not adjusted to this datum; the model lies within 0.1 % of the measurement. The tolerance of 2 % covers the handbook uncertainty of the limiting conductivities, the water viscosity correlation and the choice of a (0.3–0.45 nm moves the result by ±0.2 %). The infinite-dilution (Nernst–Einstein) option is 6 % high at this concentration.'
  },
  {
    id: 's11-exp-kcl-01-molal-conductivity',
    suite: 11,
    title: 'KCl primary conductivity standard, 0.1 mol/kg at 25 °C: measured electrolytic conductivity',
    mode: 'ionic',
    inputs: { mode: 'ionic', system: 'kcl', D_src: 'table', cond_model: 'dho', c_cat: 99.43, c_an: 99.43, T: 25 },
    kpi: 'sigma',
    expected: 1.28246,
    tol: 0.05,
    kind: 'experiment',
    measured: { value: 12824.6, uncertainty: null, unit: 'µS/cm' },
    source: { authors: 'Wu, Y. C.; Koch, W. F.; Pratt, K. W.', title: 'Proposed New Electrolytic Conductivity Primary Standards for KCl Solutions', journal: 'J. Res. Natl. Inst. Stand. Technol. 96 (2), 191–201', year: 1991, doi: '10.6028/jres.096.008', url: 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930028/' },
    conditions: 'Aqueous KCl, molality 0.1 mol/kg, 25 °C; absolute (a.c. and d.c.) conductivity determination in a cell of known geometry; overall uncertainty stated as 0.03 % over 0–50 °C',
    match: ['system', 'c_cat', 'c_an', 'T'],
    reference: 'Wu, Koch & Pratt, J. Res. NIST 96, 191 (1991); the same value is carried in NIST Special Publication 260-142',
    note: 'Measured 0.0128246 S/cm. 0.1 mol/kg corresponds to 99.43 mol/m³ (solution density 1001.7 kg/m³, handbook value, not from the source). Same model and default ion-size parameter as the 0.01 mol/kg entry; the model is about 1 % low. At this concentration the finite-ion-size factor matters (a = 0.3–0.45 nm changes the result by −1 to +1 %; the limiting Onsager law without it is 7 % low), hence the 5 % tolerance; the theory is not offered above about 0.1 mol/L.'
  },
];
