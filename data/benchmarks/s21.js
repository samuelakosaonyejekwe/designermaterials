// Reference benchmarks for Suite 21 — Fabrication, Synthesis and Process Engineering.
// Every expected value is a literal handbook figure or is evaluated here from a closed-form expression with standard
// constants; nothing is taken from the suite's output.

const R = 8.314462618, F = 96485.33212, KB = 1.380649e-23, NA = 6.02214076e23, EPS0 = 8.8541878128e-12, PI = Math.PI;
// standard atomic weights (g/mol)
const AW = { H: 1.008, C: 12.011, N: 14.007, O: 15.999, Na: 22.990, Al: 26.982, Si: 28.085, Cu: 63.546, Zn: 65.38 };

// ---------------------------------------------------------------------------------------------------------------
// ZnO porous reactor — a nearly inert, weakly heated slab (low field, almost adiabatic faces, negligible kinetics)
// ---------------------------------------------------------------------------------------------------------------
const ZJ = { mode: 'zno_reactor', T0: 298.15, Tinf: 298.15, E: 100, h: 1, A_pre: 1e-3, field_nu: 0 };
const phi = 0.25, cZn = 10, cOH = 20, Lz = 0.01, rhoCp = 1000 * 4186, tRun = 1000;
// Nernst–Einstein conductivity with Bruggeman tortuosity; limiting ionic diffusivities at 25 °C: Zn²⁺ 7.5e-10 (suite default),
// OH⁻ 5.3e-9 (suite default), NO₃⁻ 1.902e-9 and Na⁺ 1.334e-9 m²/s (handbook values). Nitrate 2·cZn, sodium cOH by electroneutrality.
const sigma0 = phi ** 1.5 * F * F / (R * 298.15) * (4 * 7.5e-10 * cZn + 5.3e-9 * cOH + 1.902e-9 * 2 * cZn + 1.334e-9 * cOH);
const qJoule = sigma0 * ZJ.E ** 2;                      // W/m³
const tauSlab = rhoCp * Lz / (2 * ZJ.h);                // lumped slab cooled on both faces
const dTlumped = qJoule / rhoCp * tauSlab * (1 - Math.exp(-tRun / tauSlab));

// ---------------------------------------------------------------------------------------------------------------
// CVD — Grove model with a laminar boundary layer and plug-flow depletion (closed form of the stated model)
// ---------------------------------------------------------------------------------------------------------------
function grove({ T, P, x, U, L, H, ks }) {
  // property model stated by the suite for SiH₄ in H₂: D = 6e-5 (T/300)^1.75 (101325/P) m²/s, µ = 8.9e-6 (T/300)^0.7 Pa·s
  const D = 6.0e-5 * (T / 300) ** 1.75 * 101325 / P, mu = 8.9e-6 * (T / 300) ** 0.7, rho = P * 2.016e-3 / (R * T);
  const Sc = mu / (rho * D), ReL = rho * U * L / mu, C0 = x * P / (KB * T), Nfilm = 2329 * NA / (AW.Si * 1e-3);
  const hg = 0.664 * D / L * Math.sqrt(ReL) * Sc ** (1 / 3);      // plate-average laminar mass-transfer coefficient
  const c = 0.332 * D * Math.sqrt(rho * U / mu) * Sc ** (1 / 3);  // local h(x) = c / √x
  // ∫₀ᴸ ks·h/(ks + h) dx = 2c[√L − (c/ks) ln(1 + ks√L/c)]; plug flow: C(L) = C0 exp(−∫/(U H)); mean rate = U H (C0 − C(L))/(L N)
  const I = 2 * c * (Math.sqrt(L) - c / ks * Math.log(1 + ks * Math.sqrt(L) / c));
  return { hg, C0, Nfilm, mean: U * H * C0 * (1 - Math.exp(-I / (U * H))) / (L * Nfilm), hfd: 4.86 * D / (2 * H), hL: c / Math.sqrt(L) };
}
const vbar = (T, Mkg) => Math.sqrt(8 * R * T / (PI * Mkg));
const M_SIH4 = (AW.Si + 4 * AW.H) * 1e-3;
const ksSi = T => 5.9e6 * Math.exp(-160e3 / (R * T));            // suite's stated first-order Arrhenius calibration for silane
const T_LP = 620 + 273.15, lp = grove({ T: T_LP, P: 100, x: 0.3, U: 1, L: 0.2, H: 0.05, ks: ksSi(T_LP) });
const T_AP = 1200 + 273.15;
const apMix = grove({ T: T_AP, P: 101325, x: 0.01, U: 10, L: 0.2, H: 0.05, ks: ksSi(T_AP) });
const apMT = grove({ T: T_AP, P: 101325, x: 0.01, U: 10, L: 0.2, H: 0.5, ks: vbar(T_AP, M_SIH4) / 4 });

// ---------------------------------------------------------------------------------------------------------------
// ALD — impingement flux, Langmuir saturation and the Gordon exposure for holes (TMA, 4 Al nm⁻² per cycle)
// ---------------------------------------------------------------------------------------------------------------
const mTMA = (AW.Al + 3 * AW.C + 9 * AW.H) * 1e-3 / NA, T_ALD = 473.15, rootALD = Math.sqrt(2 * PI * mTMA * KB * T_ALD);
const aldSat = 1 - Math.exp(-1e-3 * (0.01 / rootALD) * 0.005 / 4e18);
const gordon = a => 4e18 * rootALD * (1 + 19 * a / 4 + 1.5 * a * a);

// ---------------------------------------------------------------------------------------------------------------
// PVD — thermal evaporation of copper from a small planar source (Hertz–Knudsen + cosine law)
// ---------------------------------------------------------------------------------------------------------------
const mCu = AW.Cu * 1e-3 / NA, Ts = 1300 + 273.15, hThrow = 0.10, Rw = 0.05;
const PvCu = 101325 * Math.exp(-300e3 / R * (1 / Ts - 1 / 2835));      // Clausius–Clapeyron anchored at T_b = 2835 K, ΔH_vap = 300 kJ/mol
const nDot = PvCu / Math.sqrt(2 * PI * mCu * KB * Ts) * 1e-4;          // atoms/s from 1 cm²
const rateCu = nDot / (PI * hThrow ** 2) * (AW.Cu * 1e-3 / (8960 * NA)); // m/s at the wafer centre
const edge = (1 + (Rw / hThrow) ** 2) ** -2;                            // t(R)/t(0) = cos⁴θ for a small planar source facing a parallel wafer
const stoneyBow = 3 * 200e6 * (rateCu * 600) * Rw ** 2 / (180.5e9 * (525e-6) ** 2); // δ = κR²/2, κ = 6σt_f/(M_s t_s²), M_s = 180.5 GPa for (100) Si

// ---------------------------------------------------------------------------------------------------------------
// Additive — Rosenthal moving point source: widest extent of the melting isotherm (evaluated numerically from the closed form)
// ---------------------------------------------------------------------------------------------------------------
function rosenthalWidth({ A, P, v, k, rho, cp, Tm, T0 }) {
  const al = k / (rho * cp), rMax = A * P / (2 * PI * k * (Tm - T0)); // isotherm radius directly behind the source
  let yMax = 0;
  for (let i = 1; i <= 200000; i++) { const r = rMax * i / 200000, xi = -r + 2 * al / v * Math.log(rMax / r); if (Math.abs(xi) <= r) { const y = Math.sqrt(r * r - xi * xi); if (y > yMax) yMax = y; } }
  return 2 * yMax;
}
const wTi64 = rosenthalWidth({ A: 0.4, P: 200, v: 1, k: 15, rho: 4420, cp: 610, Tm: 1923, T0: 373.15 }); // suite's Ti-6Al-4V property set

// ---------------------------------------------------------------------------------------------------------------
// Sol–gel / spin coating (Meyerhofer with evaporation rate e = C√ω)
// ---------------------------------------------------------------------------------------------------------------
const x0 = 0.5e3 * (AW.Si + 2 * AW.O) * 1e-3 / (2200 * 0.5); // solids (gel) volume fraction: 0.5 mol/L silica at 2200 kg/m³, gel relative density 0.5
const meyerhofer = rpm => { const w = rpm * 2 * PI / 60; return x0 * Math.cbrt(3 * 3e-3 * 3e-8 * Math.sqrt(w) / (2 * (1 - x0) * 900 * w * w)); };

// ---------------------------------------------------------------------------------------------------------------
// Electrospinning and colloid batch
// ---------------------------------------------------------------------------------------------------------------
const taylorKV = Math.sqrt(4 * 15 ** 2 / 4 ** 2 * (Math.log(2 * 4 / 0.04) - 1.5) * 0.117 * PI * 0.04 * 40); // kV with cm and dyn/cm
const Qes = 1e-6 / 3600, Ies = 100e-9;
const jet = Math.cbrt(0.040 * EPS0 * (Qes / Ies) ** 2 * 2 / (PI * (2 * Math.log(100) - 3)));
const M_ZNO = AW.Zn + AW.O, M_NANO3 = AW.Na + AW.N + 3 * AW.O;

export default [
  // ---- ZnO reactor -------------------------------------------------------------------------------------------
  { id: 's21-zno-nernst-einstein-conductivity', suite: 21, title: 'ZnO reactor: Nernst–Einstein conductivity of the pore electrolyte', mode: 'zno_reactor', inputs: ZJ,
    pub: 'eff.sigma', expected: sigma0, tol: 0.01, kind: 'analytical', reference: 'Newman & Thomas-Alyea, Electrochemical Systems — dilute-solution conductivity σ = (F²/RT) Σ zᵢ² Dᵢ cᵢ; Bruggeman correction φ^1.5',
    note: '10 mol/m³ Zn(NO₃)₂ + 20 mol/m³ NaOH at 25 °C in a bed of porosity 0.25. The run warms by 0.22 K, which raises the diffusivities by about 0.5 % (Stokes–Einstein), hence the 1 % tolerance.' },
  { id: 's21-zno-joule-power', suite: 21, title: 'ZnO reactor: Joule power σE²·V of a uniform bed', mode: 'zno_reactor', inputs: ZJ,
    kpi: 'P_joule', expected: qJoule * 1 * 1 * Lz, tol: 0.01, kind: 'analytical', reference: 'Joule’s law, volumetric heating q = σ E²',
    note: 'Uniform 100 V/m field across a 1 m × 1 m × 0.01 m bed; the time average includes the slight conductivity rise with temperature.' },
  { id: 's21-zno-lumped-temperature-rise', suite: 21, title: 'ZnO reactor: lumped-capacitance temperature rise under Joule heating', mode: 'zno_reactor', inputs: ZJ,
    kpi: 'T_max', expected: 298.15 + dTlumped, abs: 0.005, kind: 'analytical', reference: 'Incropera & DeWitt, Fundamentals of Heat and Mass Transfer — lumped-capacitance transient with uniform generation',
    note: 'ΔT = (σE²/ρc_p)·τ·(1 − e^(−t/τ)), τ = ρc_pL/(2h); Biot number 0.008. The absolute tolerance is 2 % of the 0.22 K rise; the hottest cell lies slightly above the slab mean.' },
  { id: 's21-zno-zinc-conservation', suite: 21, title: 'ZnO reactor reference case: zinc inventory closes', mode: 'zno_reactor', inputs: { mode: 'zno_reactor' },
    kpi: 'bal_zn', expected: 0, abs: 1e-6, kind: 'conservation', reference: 'Patankar, Numerical Heat Transfer and Fluid Flow — conservative finite-volume discretisation',
    note: 'Pore zinc + ZnO formed − initial inventory − integrated boundary flux, relative to the largest term.' },
  { id: 's21-zno-hydroxide-conservation', suite: 21, title: 'ZnO reactor reference case: hydroxide inventory closes', mode: 'zno_reactor', inputs: { mode: 'zno_reactor' },
    kpi: 'bal_oh', expected: 0, abs: 1e-6, kind: 'conservation', reference: 'Patankar, Numerical Heat Transfer and Fluid Flow — conservative finite-volume discretisation' },
  { id: 's21-zno-energy-conservation', suite: 21, title: 'ZnO reactor reference case: energy balance closes', mode: 'zno_reactor', inputs: { mode: 'zno_reactor' },
    kpi: 'bal_energy', expected: 0, abs: 1e-6, kind: 'conservation', reference: 'First law of thermodynamics: stored heat = Joule heat + reaction heat − convective loss' },
  { id: 's21-zno-stoichiometry', suite: 21, title: 'ZnO reactor reference case: two hydroxide ions consumed per ZnO formed', mode: 'zno_reactor', inputs: { mode: 'zno_reactor' },
    kpi: 'oh_per_zno', expected: 2, tol: 1e-3, kind: 'conservation', reference: 'Stoichiometry of Zn²⁺ + 2 OH⁻ → ZnO + H₂O',
    note: 'The two species sinks are coupled by Picard iteration, so the ratio closes to the coupling tolerance rather than to round-off.' },
  { id: 's21-zno-diffusion-length', suite: 21, title: 'ZnO reactor: Zn²⁺ diffusion length √(D_pore t)', mode: 'zno_reactor', inputs: { mode: 'zno_reactor' },
    kpi: 'len_diff', expected: Math.sqrt(Math.sqrt(0.25) * 7.5e-10 * 1000), tol: 1e-6, kind: 'analytical', reference: 'Crank, The Mathematics of Diffusion — diffusion length; Bruggeman relation D_eff = φ^1.5 D',
    note: 'Pore diffusivity D_eff/φ = φ^0.5·D = 3.75×10⁻¹⁰ m²/s over 1000 s gives 0.61 mm.' },

  // ---- CVD ---------------------------------------------------------------------------------------------------
  { id: 's21-cvd-reaction-limited', suite: 21, title: 'LPCVD polysilicon: surface-reaction limit of the Grove model, v → k_s C_g/N', mode: 'cvd', inputs: { mode: 'cvd' },
    kpi: 'rate', expected: ksSi(T_LP) * lp.C0 / lp.Nfilm * 6e10, tol: 0.01, kind: 'analytical', reference: 'Grove, Physics and Technology of Semiconductor Devices — CVD growth-rate model',
    note: 'k_s/h_g ≈ 10⁻⁴ at 620 °C and 100 Pa. The susceptor-mean rate lies 0.5 % below the limit because 1 % of the silane is depleted along the flow. k_s uses the suite’s stated Arrhenius calibration.' },
  { id: 's21-cvd-mass-transport-limited', suite: 21, title: 'Atmospheric CVD: mass-transport limit of the Grove model, v → h_g C_g/N', mode: 'cvd',
    inputs: { mode: 'cvd', cvd_kind: 'pecvd', cvd_s: 1, cvd_T: 1200, cvd_P: 101325, cvd_x: 0.01, cvd_U: 10, cvd_L: 0.2, cvd_H: 0.5 },
    kpi: 'rate', expected: apMT.hg * apMT.C0 / apMT.Nfilm * 6e10, tol: 0.01, kind: 'analytical', reference: 'Grove, Physics and Technology of Semiconductor Devices; Incropera & DeWitt — laminar flat-plate average Sh = 0.664 Re^½ Sc^⅓',
    note: 'Unit sticking probability makes k_s = v̄/4 ≫ h_g (Damköhler number ≈ 2000); Re ≈ 1200. Residual deviation is precursor depletion (0.3 %). Gas properties follow the suite’s stated D ∝ T^1.75/P and µ ∝ T^0.7 laws.' },
  { id: 's21-cvd-mixed-control-depletion', suite: 21, title: 'Atmospheric silicon CVD: mixed control with plug-flow depletion (closed-form integral)', mode: 'cvd',
    inputs: { mode: 'cvd', cvd_T: 1200, cvd_P: 101325, cvd_x: 0.01, cvd_U: 10, cvd_L: 0.2, cvd_H: 0.05 },
    kpi: 'rate', expected: apMix.mean * 6e10, tol: 0.005, kind: 'analytical', reference: 'Grove, Physics and Technology of Semiconductor Devices — series resistance k_s h_g/(k_s + h_g) integrated along a plug-flow channel',
    note: 'Mean rate = U·H·C₀(1 − exp(−∫K dx/UH))/(L·N) with ∫K dx = 2c[√L − (c/k_s) ln(1 + k_s√L/c)], h(x) = c/√x. Tests the suite’s marching integration.' },
  { id: 's21-cvd-step-coverage', suite: 21, title: 'CVD step coverage in a trench: Thiele-type profile 1/cosh(AR·√(3s))', mode: 'cvd',
    inputs: { mode: 'cvd', cvd_kind: 'pecvd', cvd_s: 0.1, cvd_AR: 5 },
    kpi: 'step', expected: 100 / Math.cosh(5 * Math.sqrt(0.3)), tol: 1e-4, kind: 'analytical', reference: 'Knudsen diffusion–reaction (Thiele modulus) model of deposition in a trench',
    note: 'Bottom-to-top thickness ratio for reactive sticking probability 0.1 and aspect ratio 5.' },

  // ---- ALD ---------------------------------------------------------------------------------------------------
  { id: 's21-ald-alumina-gpc', suite: 21, title: 'ALD Al₂O₃ from TMA/H₂O: growth per cycle', mode: 'ald', inputs: { mode: 'ald' },
    kpi: 'gpc', expected: 0.11, tol: 0.05, kind: 'reference-data', reference: 'Puurunen, Surface chemistry of atomic layer deposition: a case study for the trimethylaluminum/water process (J. Appl. Phys.) — about 1.1 Å per cycle near 200 °C',
    note: 'Measured values scatter between about 0.09 and 0.12 nm per cycle with temperature and reactor; the suite uses 4 Al nm⁻² per cycle and a 3000 kg/m³ film.' },
  { id: 's21-ald-langmuir-saturation', suite: 21, title: 'ALD: Langmuir saturation fraction of an under-dosed pulse', mode: 'ald', inputs: { mode: 'ald', ald_P: 0.01, ald_tp: 0.005 },
    kpi: 'sat', expected: aldSat * 100, tol: 0.005, kind: 'analytical', reference: 'Kinetic theory of gases — impingement flux Φ = P/√(2πmk_BT); first-order Langmuir adsorption θ = 1 − exp(−s₀Φt/σ)',
    note: 'TMA at 200 °C, 0.01 Pa for 5 ms, s₀ = 10⁻³, 4×10¹⁸ sites/m².' },
  { id: 's21-ald-gordon-exposure', suite: 21, title: 'ALD: exposure for conformal coating of aspect-ratio-50 holes', mode: 'ald', inputs: { mode: 'ald', ald_AR: 50 },
    kpi: 'exposure', expected: gordon(50), tol: 0.002, kind: 'analytical', reference: 'Gordon, Hausmann, Kim & Shepard, A kinetic model for step coverage by atomic layer deposition in narrow holes or trenches (Chem. Vap. Deposition)',
    note: 'P·t = S√(2πmk_BT)(1 + 19a/4 + 3a²/2) with S = 4×10¹⁸ m⁻² for TMA at 200 °C.' },

  // ---- PVD ---------------------------------------------------------------------------------------------------
  { id: 's21-pvd-evaporation-cos4-uniformity', suite: 21, title: 'Evaporation from a small planar source: cos⁴ thickness fall-off at the wafer edge', mode: 'pvd', inputs: { mode: 'pvd', pvd_proc: 'evap', pvd_h: 10, pvd_wafer: 100 },
    kpi: 'uniformity', expected: (1 - edge) / (1 + edge) * 100, tol: 1e-4, kind: 'analytical', reference: 'Ohring, Materials Science of Thin Films — surface-source thickness distribution t/t₀ = [1 + (r/h)²]⁻²',
    note: 'r/h = 0.5 gives t(edge)/t(centre) = 0.64, i.e. ±21.95 %.' },
  { id: 's21-pvd-evaporation-collection', suite: 21, title: 'Evaporation: fraction of a cosine source collected by a coaxial disc', mode: 'pvd', inputs: { mode: 'pvd', pvd_proc: 'evap', pvd_h: 10, pvd_wafer: 100 },
    kpi: 'collection', expected: Rw ** 2 / (Rw ** 2 + hThrow ** 2) * 100, tol: 0.001, kind: 'analytical', reference: 'Lambert cosine law — view factor from a small surface element to a coaxial parallel disc, R²/(R² + h²)' },
  { id: 's21-pvd-hertz-knudsen-rate', suite: 21, title: 'Evaporation of copper at 1300 °C: Hertz–Knudsen centre rate', mode: 'pvd', inputs: { mode: 'pvd', pvd_proc: 'evap', pvd_mat: 'Cu', pvd_Tsrc: 1300, pvd_Asrc: 1, pvd_h: 10 },
    kpi: 'rate', expected: rateCu * 6e10, tol: 0.005, kind: 'analytical', reference: 'Ohring, Materials Science of Thin Films — Hertz–Knudsen equation and cosine-law deposition',
    note: 'Vapour pressure from Clausius–Clapeyron anchored at the boiling point (2835 K, 300 kJ/mol), as the suite states; that estimate is itself only good to a factor of a few against tabulated vapour pressures.' },
  { id: 's21-pvd-thermal-mismatch-stress', suite: 21, title: 'Copper film on silicon deposited at 100 °C: thermal-mismatch stress', mode: 'pvd', inputs: { mode: 'pvd', pvd_mat: 'Cu', pvd_Tdep: 100, pvd_sint: 0 },
    kpi: 'stress', expected: 120e9 / (1 - 0.34) * (16.5e-6 - 2.6e-6) * 75 / 1e6, tol: 0.002, kind: 'analytical', reference: 'Freund & Suresh, Thin Film Materials — biaxial mismatch stress E_f/(1−ν_f)·Δα·ΔT',
    note: 'Handbook-typical Cu: E = 120 GPa, ν = 0.34, α = 16.5×10⁻⁶ K⁻¹; Si α = 2.6×10⁻⁶ K⁻¹. Tensile on cooling to 25 °C.' },
  { id: 's21-pvd-stoney-bow', suite: 21, title: 'Wafer bow from a 200 MPa film stress (Stoney)', mode: 'pvd', inputs: { mode: 'pvd', pvd_proc: 'evap', pvd_mat: 'Cu', pvd_Tdep: 25, pvd_sint: 200, pvd_tsub: 525, pvd_wafer: 100, pvd_t: 600 },
    kpi: 'bow', expected: stoneyBow * 1e6, tol: 0.01, kind: 'analytical', reference: 'Stoney equation; Freund & Suresh, Thin Film Materials — (100) silicon biaxial modulus 180.5 GPa',
    note: 'δ = 3σt_fR²/(M_s t_s²) with the film thickness from the Hertz–Knudsen benchmark above; 525 µm wafer, 50 mm radius.' },
  { id: 's21-pvd-knudsen-number', suite: 21, title: 'Sputtering gas: mean free path of argon at 0.5 Pa over a 10 cm throw', mode: 'pvd', inputs: { mode: 'pvd', pvd_P: 0.5, pvd_h: 10 },
    kpi: 'Kn', expected: KB * 300 / (Math.SQRT2 * PI * (3.6e-10) ** 2 * 0.5) / 0.10, tol: 0.002, kind: 'analytical', reference: 'Kinetic theory of gases — mean free path λ = k_BT/(√2 π d² P)',
    note: 'Hard-sphere diameter 0.36 nm at 300 K gives λ = 1.44 cm.' },

  // ---- Lithography and etch ----------------------------------------------------------------------------------
  { id: 's21-litho-rayleigh-resolution', suite: 21, title: 'Dry ArF lithography: Rayleigh resolution k₁λ/NA', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'optical', lit_src: 'arf', lit_NA: 0.93, lit_k1: 0.3 },
    kpi: 'res', expected: 0.3 * 193 / 0.93, tol: 1e-4, kind: 'analytical', reference: 'Mack, Fundamental Principles of Optical Lithography — Rayleigh criteria' },
  { id: 's21-litho-rayleigh-dof', suite: 21, title: 'Dry ArF lithography: Rayleigh depth of focus k₂λ/NA²', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'optical', lit_src: 'arf', lit_NA: 0.93, lit_k2: 0.5 },
    kpi: 'dof', expected: 0.5 * 193 / 0.93 ** 2, tol: 1e-4, kind: 'analytical', reference: 'Mack, Fundamental Principles of Optical Lithography — Rayleigh criteria' },
  { id: 's21-litho-immersion-resolution', suite: 21, title: 'ArF immersion lithography at NA 1.35: Rayleigh resolution', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'optical', lit_src: 'arfi', lit_NA: 1.35, lit_k1: 0.3 },
    kpi: 'res', expected: 0.3 * 193 / 1.35, tol: 1e-4, kind: 'analytical', reference: 'Mack, Fundamental Principles of Optical Lithography — Rayleigh criteria' },
  { id: 's21-ebeam-write-time', suite: 21, title: 'Electron-beam lithography: write time = dose × exposed area / current', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'ebeam', eb_dose: 600, lit_area: 1, eb_I: 1, lit_cd: 70, lit_pitch: 140 },
    kpi: 'twrite', expected: 600e-6 * 0.01 * 0.5 / 1e-9 / 60, tol: 1e-4, kind: 'analytical', reference: 'Charge conservation: t = D·A/I',
    note: '600 µC/cm² over 1 mm² at 50 % pattern density with 1 nA: 3000 s.' },
  { id: 's21-nil-stefan-time', suite: 21, title: 'Nanoimprint: Stefan squeeze-flow time under a line protrusion', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'nil', lit_cd: 70, nil_visc: 1e4, nil_p: 5, nil_h0: 200, nil_hf: 20 },
    kpi: 'tfill', expected: 1e4 * (70e-9) ** 2 / (2 * 5e6) * (1 / (20e-9) ** 2 - 1 / (200e-9) ** 2), tol: 0.001, kind: 'analytical', reference: 'Stefan equation for squeeze flow of a Newtonian film under a long strip, t = (ηs²/2p)(1/h_f² − 1/h₀²)' },
  { id: 's21-etch-arde-time', suite: 21, title: 'Reactive-ion etch with aspect-ratio-dependent rate: time to depth', mode: 'litho_etch', inputs: { mode: 'litho_etch', lit_tool: 'nil', lit_cd: 70, et_rate: 200, et_depth: 100, et_S: 0.3 },
    kpi: 'tetch', expected: (100 + 0.375 * 0.3 * 100 ** 2 / 70) / (200 / 60), tol: 0.005, kind: 'analytical', reference: 'Coburn & Winters, Conductance considerations in the reactive ion etching of high aspect ratio features (Appl. Phys. Lett.)',
    note: 'dh/dt = r₀K/(K + S − KS), K = 1/(1 + ¾h/w) integrates to t = (H + ⅜·S·H²/w)/r₀ for constant width; the suite also narrows the feature slightly by undercut.' },

  // ---- Additive manufacturing --------------------------------------------------------------------------------
  { id: 's21-lpbf-energy-density', suite: 21, title: 'Laser powder-bed fusion: volumetric energy density P/(v·h·t)', mode: 'additive', inputs: { mode: 'additive', am_proc: 'lpbf', am_P: 200, am_v: 1000, am_h: 100, am_t: 30 },
    kpi: 'ved', expected: 200 / (1000 * 0.1 * 0.03), tol: 1e-6, kind: 'analytical', reference: 'Definition of volumetric energy density in laser powder-bed fusion' },
  { id: 's21-lpbf-rosenthal-width', suite: 21, title: 'Ti-6Al-4V melt-pool width from the Rosenthal moving point source', mode: 'additive', inputs: { mode: 'additive', am_proc: 'lpbf', am_mat: 'Ti64', am_P: 200, am_v: 1000, am_T0: 100 },
    kpi: 'W', expected: wTi64 * 1e6, tol: 0.01, kind: 'analytical', reference: 'Rosenthal, The theory of moving sources of heat and its application to metal treatments (Trans. ASME)',
    note: 'Widest extent of the 1923 K isotherm of T − T₀ = AP/(2πkr)·exp(−v(r+ξ)/2α), found here by a dense scan; absorptivity 0.4 and the suite’s Ti-6Al-4V properties (k = 15 W/m/K, ρ = 4420 kg/m³, c_p = 610 J/kg/K).' },
  { id: 's21-sla-jacobs-cure-depth', suite: 21, title: 'Stereolithography: Jacobs cure depth D_p ln(E/E_c)', mode: 'additive', inputs: { mode: 'additive', am_proc: 'sla', sla_Dp: 150, sla_E: 15, sla_Ec: 10 },
    kpi: 'cure', expected: 150 * Math.log(1.5), tol: 1e-5, kind: 'analytical', reference: 'Jacobs, Rapid Prototyping & Manufacturing: Fundamentals of StereoLithography — working curve' },
  { id: 's21-sla-line-width', suite: 21, title: 'Stereolithography: cured line width of a Gaussian beam', mode: 'additive', inputs: { mode: 'additive', am_proc: 'sla', sla_Dp: 150, sla_E: 15, sla_Ec: 10, sla_w0: 70 },
    kpi: 'lw', expected: 70 * Math.sqrt(2 * Math.log(1.5)), tol: 1e-5, kind: 'analytical', reference: 'Jacobs, Rapid Prototyping & Manufacturing: Fundamentals of StereoLithography — L_w = 2w₀√(C_d/2D_p)' },
  { id: 's21-tpp-voxel-width', suite: 21, title: 'Two-photon polymerisation: lateral voxel size above threshold', mode: 'additive', inputs: { mode: 'additive', am_proc: 'tpp', tpp_lam: 780, tpp_NA: 1.4, tpp_P: 15, tpp_Pth: 1, tpp_v: 5000 },
    kpi: 'dvox', expected: 0.61 * 780 / 1.4 * Math.sqrt(2 * Math.log(15 / Math.sqrt(50))), tol: 0.001, kind: 'analytical', reference: 'Threshold model of two-photon polymerisation with a Gaussian focus, d = w₀√(2 ln(P/P_th)); Airy radius 0.61λ/NA',
    note: 'Threshold power scaled as √v from 1 mW at 100 µm/s to 7.07 mW at 5 mm/s (two-photon dose ∝ I²t), as the suite states.' },
  { id: 's21-lattice-bcc-strut', suite: 21, title: 'BCC strut lattice: strut diameter for 20 % relative density', mode: 'additive', inputs: { mode: 'additive', am_proc: 'sla', am_rd: 0.2, am_cell: 2 },
    kpi: 'strut', expected: 2 * 2000 * Math.sqrt(0.2 / (4 * Math.sqrt(3) * PI)), tol: 1e-4, kind: 'analytical', reference: 'Gibson & Ashby, Cellular Solids — slender-strut relative density; BCC cell with four body diagonals, ρ* = 4√3π(r/a)²',
    note: 'First-order strut-volume estimate without node-overlap correction.' },

  // ---- Sol–gel, spin coating, drying -------------------------------------------------------------------------
  { id: 's21-solgel-gel-time-arrhenius', suite: 21, title: 'Sol–gel: Arrhenius scaling of gel time from 25 °C to 50 °C', mode: 'solgel_spin_sinter', inputs: { mode: 'solgel_spin_sinter', sg_T: 50, sg_tref: 100, sg_Ea: 60 },
    kpi: 'tgel', expected: 100 * Math.exp(60e3 / R * (1 / 323.15 - 1 / 298.15)), tol: 1e-4, kind: 'analytical', reference: 'Brinker & Scherer, Sol–Gel Science — Arrhenius temperature dependence of gel time' },
  { id: 's21-spin-meyerhofer-3000rpm', suite: 21, title: 'Spin coating at 3000 rpm: Meyerhofer film thickness', mode: 'solgel_spin_sinter', inputs: { mode: 'solgel_spin_sinter', sp_rpm: 3000 },
    kpi: 'hgel', expected: meyerhofer(3000) * 1e9, tol: 0.002, kind: 'analytical', reference: 'Meyerhofer, Characteristics of resist films produced by spinning (J. Appl. Phys.)',
    note: 'h = x₀[3ηe/(2(1−x₀)ρω²)]^⅓ with e = C√ω; silica molar mass from standard atomic weights.' },
  { id: 's21-spin-meyerhofer-9000-scaling', suite: 21, title: 'Spin coating at three times the speed: thickness falls by √3 (h ∝ ω^−½)', mode: 'solgel_spin_sinter', inputs: { mode: 'solgel_spin_sinter', sp_rpm: 9000 },
    kpi: 'hgel', expected: meyerhofer(3000) * 1e9 / Math.sqrt(3), tol: 0.002, kind: 'analytical', reference: 'Meyerhofer, Characteristics of resist films produced by spinning (J. Appl. Phys.)',
    note: 'With the evaporation rate proportional to √ω the Meyerhofer thickness scales as ω^−½, so tripling the speed divides the thickness by √3. (The exponent −⅔ applies only to a speed-independent evaporation rate.)' },
  { id: 's21-drying-critical-cracking-thickness', suite: 21, title: 'Drying gel film: critical cracking thickness', mode: 'solgel_spin_sinter', inputs: { mode: 'solgel_spin_sinter', dr_r: 5, dr_gamma: 0.022, dr_KIc: 0.012 },
    kpi: 'hc', expected: (0.012e6 / (2 * 0.022 / 5e-9)) ** 2 / 1.976 * 1e9, tol: 1e-4, kind: 'analytical', reference: 'Hutchinson & Suo, Mixed mode cracking in layered materials — channelling crack, Z = 1.976; Young–Laplace capillary pressure 2γ/r',
    note: 'Drying stress taken as the full capillary pressure (8.8 MPa), an upper bound.' },

  // ---- Electrospinning and colloid batch ---------------------------------------------------------------------
  { id: 's21-electrospin-taylor-voltage', suite: 21, title: 'Electrospinning: critical voltage for a Taylor cone on a needle', mode: 'electrospin_colloid', inputs: { mode: 'electrospin_colloid', es_H: 15, es_L: 4, es_R: 0.4, es_gamma: 40 },
    kpi: 'Vc', expected: taylorKV, tol: 1e-4, kind: 'analytical', reference: 'Taylor, Electrically driven jets (Proc. R. Soc. A) — V_c² = (4H²/L²)(ln(2L/R) − 3/2)(0.117πRγ)',
    note: 'H = 15 cm, L = 4 cm, R = 0.04 cm, γ = 40 dyn/cm gives 11.2 kV.' },
  { id: 's21-electrospin-terminal-jet', suite: 21, title: 'Electrospinning: terminal jet diameter from charge repulsion versus surface tension', mode: 'electrospin_colloid', inputs: { mode: 'electrospin_colloid', es_gamma: 40, es_Q: 1, es_I: 100 },
    kpi: 'jet', expected: jet * 1e9, tol: 0.001, kind: 'analytical', reference: 'Fridrikh, Yu, Brenner & Rutledge, Controlling the fiber diameter during electrospinning (Phys. Rev. Lett.)',
    note: 'h_t = [γε₀(Q/I)²·2/(π(2 ln χ − 3))]^⅓ with χ = 100, the value the suite states.' },
  { id: 's21-colloid-zno-mass', suite: 21, title: 'Colloid batch: ZnO recovered from 0.1 mol zinc at 90 % conversion after three washes', mode: 'electrospin_colloid', inputs: { mode: 'electrospin_colloid', cZn: 10, col_V: 10, col_conv: 0.9, col_wash: 3, col_loss: 2 },
    kpi: 'mZnO', expected: 0.1 * 0.9 * M_ZNO * 0.98 ** 3, tol: 0.001, kind: 'analytical', reference: 'Stoichiometry of Zn²⁺ + 2 OH⁻ → ZnO + H₂O with standard atomic weights' },
  { id: 's21-colloid-byproduct-salt', suite: 21, title: 'Colloid batch: sodium nitrate by-product', mode: 'electrospin_colloid', inputs: { mode: 'electrospin_colloid', salt: 'Zn(NO3)2', base: 'NaOH', cZn: 10, col_V: 10, col_conv: 0.9 },
    kpi: 'byprod', expected: 2 * 0.1 * 0.9 * M_NANO3, tol: 0.001, kind: 'analytical', reference: 'Stoichiometry of Zn(NO₃)₂ + 2 NaOH → ZnO + H₂O + 2 NaNO₃' },
  { id: 's21-colloid-wash-dilution', suite: 21, title: 'Colloid batch: residual salt after three dilution washes', mode: 'electrospin_colloid', inputs: { mode: 'electrospin_colloid', col_wash: 3, col_wratio: 5 },
    kpi: 'resid', expected: 100 / 6 ** 3, tol: 1e-5, kind: 'analytical', reference: 'Repeated-dilution washing: residual fraction (1/(1 + w))ⁿ' },
];
