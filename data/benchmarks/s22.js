// Reference benchmarks for Suite 22 — material behaviour, damage and lifetime prediction.
// Expected values are handbook figures or closed-form results evaluated here from standard constants;
// none of them is taken from the suite's own output.
const R = 8.314462618, KB_EV = 8.617333262e-5, KB = 1.380649e-23, FARADAY = 96485.33212, PI = Math.PI;
const YEAR_S = 365.25 * 86400, YEAR_H = 8766;                       // Julian year
const K = c => c + 273.15;

// --- fatigue: Basquin with σ_f′ = 1000 MPa and b = −0.1 gives round numbers
const FAT = { fat_Kt: 1, fat_q: 1, fat_sf: 1000, fat_b: -0.1, fat_su: 950 };
const basquin = (sa, sf, b) => 0.5 * Math.pow(sa / sf, 1 / b);      // cycles
const N500 = basquin(500, 1000, -0.1);                               // = 2⁹ = 512
const N250 = basquin(250, 1000, -0.1);                               // = 2¹⁹ = 524288

// --- Paris law, ΔK in MPa√m, a in m, da/dN in m/cycle
const parisLife = (C, m, Y, ds, a0, ac) => m === 2
  ? Math.log(ac / a0) / (C * Y * Y * ds * ds * PI)
  : 2 * (Math.pow(a0, 1 - m / 2) - Math.pow(ac, 1 - m / 2)) / ((m - 2) * C * Math.pow(Y * ds * Math.sqrt(PI), m));
const acOf = (KIc, Y, smax) => Math.pow(KIc / (Y * smax), 2) / PI;  // m
const AC50 = acOf(50, 1, 100);
const FRAC = { mode: 'fracture', fr_ds: 1e8, fr_R: 0, fr_Y: 1, fr_KIc: 5e7, fr_a0: 1000, fr_adet: 2000, fr_C: 1e-11, fr_dKth: 0.05 };

// --- Deal–Grove: Arrhenius fits of the rate constants as tabulated in VLSI process texts (µm, h, eV)
const dgThickness = (A, B, t, tau) => 0.5 * A * (Math.sqrt(1 + 4 * B * (t + tau) / (A * A)) - 1);   // µm
const kT1000 = KB_EV * K(1000);
const B_DRY = 772 * Math.exp(-1.23 / kT1000), BA_DRY_100 = 3.71e6 * Math.exp(-2.0 / kT1000), A_DRY_100 = B_DRY / BA_DRY_100;
const X_DRY_100 = dgThickness(A_DRY_100, B_DRY, 1, 0);
const TAU_25 = (0.025 * 0.025 + A_DRY_100 * 0.025) / B_DRY;
const X_DRY_100_25 = dgThickness(A_DRY_100, B_DRY, 1, TAU_25);
// Deal and Grove's own tabulated constants (independent of the Arrhenius fits)
const DG_DRY_1000 = { A: 0.165, B: 0.0117, tau: 0.37 }, DG_WET_1100 = { A: 0.11, B: 0.510, tau: 0 };
const XI_DRY_111 = 0.5 * (-DG_DRY_1000.A + Math.sqrt(DG_DRY_1000.A ** 2 + 4 * DG_DRY_1000.B * DG_DRY_1000.tau));   // µm, initial oxide equivalent to τ
const VM_RATIO = (28.086 / 2.329) / (60.084 / 2.20);                // molar volume of Si / amorphous SiO₂

// --- corrosion (ASTM G102 form: rate = K·i·EW/ρ with K = 3.27e-3 mm·g/(µA·cm·yr))
const g102 = (EW, rho) => 3.27e-3 * 1 * EW / rho;                   // mm/yr at 1 µA/cm²

// --- dissolution of ZnO under sink conditions
const M_ZNO = (65.38 + 15.999) * 1e-3, RHO_ZNO = 5606, K_DISS = 2e-8, R0 = 15e-9;
const T_SINK_H = R0 * RHO_ZNO / (K_DISS * M_ZNO) / 3600;

// --- LSW coarsening: time for the mean particle volume to double
const LSW = { gamma: 0.3, c: 1e-7, Vm: 14.5e-6, D: 1e-9, T: K(60), r0: 10e-9 };
const K_LSW = 8 * LSW.gamma * LSW.c * LSW.Vm ** 2 * LSW.D / (9 * R * LSW.T);
const T_DOUBLE_H = LSW.r0 ** 3 / K_LSW / 3600;

// --- thin-film mechanics
const E_PS = 150e9 / (1 - 0.25 * 0.25);                              // plane-strain modulus, Pa
const FILM = { mode: 'degradation', deg_mech: 'delamination', dl_E: 150, dl_nu: 0.25, dl_h: 1, dl_dT: 0 };

// --- rotating-squares kinematics
const S = deg => Math.cos(deg * PI / 360) + Math.sin(deg * PI / 360);

// --- exact Weibull sample (β = 2, η = 1000 h) placed on the median-rank (Benard) plotting positions
const WB = { beta: 2, eta: 1000, n: 12 };
const WB_TABLE = 'time_to_failure_h\n' + Array.from({ length: WB.n }, (_, i) => (WB.eta * Math.pow(-Math.log(1 - (i + 0.7) / (WB.n + 0.4)), 1 / WB.beta)).toPrecision(14)).join('\n');
const DEG_EXP = 'time_h,performance\n' + Array.from({ length: 11 }, (_, i) => `${500 * i},${(100 * Math.exp(-1e-4 * 500 * i)).toPrecision(14)}`).join('\n');
const DEG_LIN = 'time_h,performance\n' + Array.from({ length: 9 }, (_, i) => `${500 * i},${100 - 0.01 * 500 * i}`).join('\n');
const SQRT_PI_2 = Math.sqrt(PI) / 2;                                 // Γ(3/2)

// --- published experiment: endurance lives (millions of revolutions) of 23 deep-groove ball bearings, transcribed from the
// specimen worksheet reproduced in appendix A of Lieblein & Zelen, J. Res. Natl. Bur. Stand. 57, 273 (1956), p. 286
// (lot of 25, the two bearings marked "Omitted" on the worksheet excluded as in the source).
const LZ_TABLE = 'time_to_failure_h\n' + [17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 51.96, 54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 105.12, 105.84, 127.92, 128.04, 173.40].join('\n');
const LZ_SOURCE = { authors: 'Lieblein, J.; Zelen, M.', title: 'Statistical Investigation of the Fatigue Life of Deep-Groove Ball Bearings', journal: 'J. Res. Natl. Bur. Stand. 57 (5), 273–316, Research Paper 2719', year: 1956, doi: '10.6028/jres.057.033', url: 'https://nvlpubs.nist.gov/nistpubs/jres/057/jresv57n5p273_A1b.pdf' };
const LZ_COND = 'Endurance test of a lot of deep-groove ball bearings (bore 20 mm, outside diameter 42 mm, nine 1/4 in balls, one row) at 2000 r.p.m. with jet-oil lubrication; lives in millions of revolutions for 23 bearings (two non-fatigue failures omitted by the laboratory, three bearings marked discontinued and counted at their recorded lives)';

export default [
  // ------------------------------------------------------------------ fatigue
  { id: 's22-basquin-life', suite: 22, title: 'Basquin life at half the fatigue strength coefficient', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'none', fat_sa: 5e8, fat_sm: 0 }, kpi: 'Nf', expected: 512, tol: 1e-6, kind: 'analytical',
    reference: 'Basquin relation σ_a = σ_f′(2N)^b (Dowling, Mechanical Behavior of Materials)',
    note: 'σ_a/σ_f′ = 0.5 and b = −0.1 give 2N = 2¹⁰ exactly.' },
  { id: 's22-goodman-amplitude', suite: 22, title: 'Goodman equivalent fully reversed amplitude', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'goodman', fat_sa: 2e8, fat_sm: 475 }, kpi: 'sar', expected: 400, tol: 1e-6, kind: 'analytical',
    reference: 'Goodman mean-stress relation σ_ar = σ_a/(1 − σ_m/σ_u)', note: 'Mean stress at half the ultimate strength doubles the equivalent amplitude.' },
  { id: 's22-swt-amplitude', suite: 22, title: 'Smith–Watson–Topper equivalent amplitude', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'swt', fat_sa: 2e8, fat_sm: 200 }, kpi: 'sar', expected: Math.sqrt(400 * 200), tol: 1e-6, kind: 'analytical',
    reference: 'Smith–Watson–Topper parameter σ_ar = √(σ_max σ_a)' },
  { id: 's22-miner-two-level', suite: 22, title: 'Palmgren–Miner: two-level block spectrum, damage sum of one at failure', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'none', fat_Dcrit: 1, fat_spectrum: 'stress_amplitude_MPa,mean_stress_MPa,cycles_per_block\n500,0,64\n250,0,65536' },
    kpi: 'blocks', expected: 1 / (64 / N500 + 65536 / N250), tol: 1e-6, kind: 'conservation',
    reference: 'Palmgren–Miner linear damage rule Σ nᵢ/Nᵢ = 1', note: 'Each level contributes 1/8 per block, so failure occurs after exactly 4 blocks.' },
  { id: 's22-weibull-median', suite: 22, title: 'Weibull life scatter: failure probability at the median life', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'none', fat_sa: 5e8, fat_sm: 0, fat_beta: 3, fat_Nt: 512 }, kpi: 'Pf', expected: 50, tol: 1e-6, kind: 'analytical',
    reference: 'Two-parameter Weibull distribution, F(median) = 0.5' },
  { id: 's22-weibull-b10-fatigue', suite: 22, title: 'Weibull B10 life from the median life', mode: 'fatigue',
    inputs: { mode: 'fatigue', ...FAT, fat_mean: 'none', fat_sa: 5e8, fat_sm: 0, fat_beta: 3 }, kpi: 'B10', expected: N500 * Math.pow(-Math.log(0.9) / Math.LN2, 1 / 3), tol: 1e-6, kind: 'analytical',
    reference: 'Two-parameter Weibull distribution, B10 = η(−ln 0.9)^(1/β) with η = N_median/(ln 2)^(1/β)' },
  { id: 's22-transition-life', suite: 22, title: 'Strain–life transition life (elastic = plastic strain amplitude)', mode: 'fatigue',
    inputs: { mode: 'fatigue', fat_E: 2e11, fat_sf: 1000, fat_b: -0.1, fat_ef: 0.5, fat_c: -0.6 }, kpi: 'Nt_trans', expected: 0.5 * Math.pow(0.5 * 2e11 / 1e9, 1 / (-0.1 + 0.6)), tol: 1e-6, kind: 'analytical',
    reference: 'Coffin–Manson and Basquin strain–life lines, 2N_t = (ε_f′E/σ_f′)^(1/(b−c))', note: 'ε_f′E/σ_f′ = 100 and b − c = 0.5 give 2N_t = 10⁴.' },

  // ------------------------------------------------------------------ fracture
  { id: 's22-critical-crack', suite: 22, title: 'Critical crack size from the fracture toughness', mode: 'fracture',
    inputs: { ...FRAC, fr_ds: 2e8, fr_KIc: 6e7 }, kpi: 'ac', expected: 0.09 / PI * 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'Linear-elastic fracture mechanics, K = Yσ√(πa) = K_Ic (Anderson, Fracture Mechanics)', note: 'K_Ic = 60 MPa√m, σ = 200 MPa, Y = 1: a_c = 0.09/π m.' },
  { id: 's22-paris-m3', suite: 22, title: 'Paris-law life, m = 3, closed-form integration', mode: 'fracture',
    inputs: { ...FRAC, fr_m: 3 }, kpi: 'Nf', expected: parisLife(1e-11, 3, 1, 100, 1e-3, AC50), tol: 1e-6, kind: 'analytical',
    reference: 'Paris–Erdogan law integrated at constant geometry factor (Anderson, Fracture Mechanics)', note: 'Δσ = 100 MPa, a₀ = 1 mm, a_c = 79.6 mm.' },
  { id: 's22-paris-m2', suite: 22, title: 'Paris-law life, m = 2 (logarithmic case)', mode: 'fracture',
    inputs: { ...FRAC, fr_m: 2 }, kpi: 'Nf', expected: parisLife(1e-11, 2, 1, 100, 1e-3, AC50), tol: 1e-6, kind: 'analytical',
    reference: 'Paris–Erdogan law, N = ln(a_c/a₀)/(C Y² Δσ² π) for m = 2' },
  { id: 's22-paris-m4-edge', suite: 22, title: 'Paris-law life, m = 4, edge crack with R = 0.5', mode: 'fracture',
    inputs: { ...FRAC, fr_m: 4, fr_Y: 1.12, fr_R: 0.5, fr_C: 1e-12, fr_a0: 500 }, kpi: 'Nf', expected: parisLife(1e-12, 4, 1.12, 100, 5e-4, acOf(50, 1.12, 200)), tol: 1e-6, kind: 'analytical',
    reference: 'Paris–Erdogan law integrated to a_c = (K_Ic/(Yσ_max))²/π with σ_max = Δσ/(1 − R)' },
  { id: 's22-threshold-flaw', suite: 22, title: 'Largest non-propagating flaw at the growth threshold', mode: 'fracture',
    inputs: { ...FRAC, fr_ds: 1.5e8, fr_dKth: 3 }, kpi: 'ath', expected: Math.pow(3 / 150, 2) / PI * 1e6, tol: 1e-6, kind: 'analytical',
    reference: 'ΔK_th = YΔσ√(πa) solved for a' },
  { id: 's22-griffith-strength', suite: 22, title: 'Griffith strength of a brittle solid with a 1 µm flaw', mode: 'fracture',
    inputs: { mode: 'fracture', fr_E: 7e10, fr_gamma: 1, fr_size: 2000 }, kpi: 'sGriffith', expected: Math.sqrt(2 * 7e10 * 1 / (PI * 1e-6)) / 1e9, tol: 1e-6, kind: 'analytical',
    reference: 'Griffith criterion σ_f = √(2Eγ/(πa))', note: 'Glass-like constants (E = 70 GPa, γ = 1 J/m²): about 0.21 GPa, far below the E/10 cap.' },

  // ------------------------------------------------------------------ creep and accelerated ageing
  { id: 's22-norton-rate', suite: 22, title: 'Norton power-law creep rate', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', cr_sigma: 100, cr_T: 600, cr_A: 0.002, cr_n: 5, cr_Q: 300 }, kpi: 'rate', expected: 0.002 * 1e10 * Math.exp(-3e5 / (R * K(600))), tol: 1e-4, kind: 'analytical',
    reference: 'Norton law ε̇ = Aσⁿexp(−Q/RT) with the CODATA gas constant' },
  { id: 's22-larson-miller', suite: 22, title: 'Larson–Miller rupture life at LMP = 25 000 and 1000 K', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', cr_sigma: 1, cr_T: 726.85, lm_a0: 25, lm_a1: -5.5, lm_C: 20 }, kpi: 'tr', expected: 1e5 / YEAR_H, tol: 1e-6, kind: 'analytical',
    reference: 'Larson–Miller parameter LMP = T(C + log₁₀ t_r), T in K and t_r in h', note: 'LMP/T − C = 5, i.e. 10⁵ h = 11.41 years.' },
  { id: 's22-arrhenius-af', suite: 22, title: 'Arrhenius acceleration factor, 0.7 eV, 40 °C to 85 °C', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', acc_Ea: 0.7, acc_Tuse: 40, acc_Ttest: 85 }, kpi: 'afT', expected: Math.exp(0.7 / KB_EV * (1 / K(40) - 1 / K(85))), tol: 1e-4, kind: 'analytical',
    reference: 'Arrhenius acceleration model with the CODATA Boltzmann constant (JEDEC-type reliability practice)', note: 'About 26.' },
  { id: 's22-peck-af', suite: 22, title: 'Peck temperature–humidity acceleration, 85 °C / 85 % RH against 40 °C / 50 % RH', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', acc_Ea: 0.7, acc_Tuse: 40, acc_Ttest: 85, acc_RHuse: 50, acc_RHtest: 85, acc_nRH: 2.7 }, kpi: 'afP', expected: Math.pow(85 / 50, 2.7) * Math.exp(0.7 / KB_EV * (1 / K(40) - 1 / K(85))), tol: 1e-4, kind: 'analytical',
    reference: 'Peck model AF = (RH_test/RH_use)ⁿ·exp[(E_a/k)(1/T_use − 1/T_test)]' },
  { id: 's22-coffin-manson-cycles', suite: 22, title: 'Coffin–Manson thermal-cycling equivalence', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', cr_years: 10, acc_dTuse: 50, acc_dTtest: 200, acc_m: 2, acc_fuse: 1 }, kpi: 'cyctest', expected: 10 * 365.25 / 16, tol: 1e-6, kind: 'analytical',
    reference: 'Coffin–Manson thermal-cycling acceleration AF = (ΔT_test/ΔT_use)^m', note: 'A fourfold temperature swing with m = 2 gives AF = 16.' },
  { id: 's22-norris-landzberg', suite: 22, title: 'Norris–Landzberg acceleration, −40/+125 °C test against a 40 ± 15 °C field cycle', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', acc_Tuse: 40, acc_dTuse: 30, acc_Ttest: 42.5, acc_dTtest: 165, acc_fuse: 2, acc_ftest: 48 }, kpi: 'afNL',
    expected: Math.pow(165 / 30, 1.9) * Math.pow(2 / 48, 1 / 3) * Math.exp(1414 * (1 / K(55) - 1 / K(125))), tol: 1e-4, kind: 'analytical',
    reference: 'Norris–Landzberg model AF = (ΔT_t/ΔT_u)^1.9 (f_u/f_t)^(1/3) exp[1414(1/T_max,u − 1/T_max,t)]',
    note: 'The suite takes the entered use and test temperatures as cycle means, so T_max = T + ΔT/2 (55 °C and 125 °C here).' },
  { id: 's22-success-run', suite: 22, title: 'Zero-failure sample size for 90 % reliability at 90 % confidence', mode: 'creep_ageing',
    inputs: { mode: 'creep_ageing', acc_rel: 0.9, acc_conf: 0.9 }, kpi: 'nsamp', expected: 22, abs: 0, kind: 'analytical',
    reference: 'Success-run theorem n = ln(1 − C)/ln R', note: 'ln 0.1/ln 0.9 = 21.85, rounded up to 22 samples.' },

  // ------------------------------------------------------------------ degradation: oxidation
  { id: 's22-deal-grove-dry-100', suite: 22, title: 'Deal–Grove dry oxidation of (100) Si, 1000 °C, 1 h, no initial oxide', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'oxidation', ox_T: 1000, ox_t: 1, ox_amb: 'dry', ox_orient: '100', ox_x0: 0 }, kpi: 'xox', expected: X_DRY_100 * 1e3, tol: 0.005, kind: 'analytical',
    reference: 'Deal–Grove linear–parabolic law with the Arrhenius rate constants tabulated in Jaeger, Introduction to Microelectronic Fabrication',
    note: 'About 38.5 nm. The linear–parabolic law alone underestimates thin dry oxides; texts add a 25 nm fictitious initial oxide (next benchmark).' },
  { id: 's22-deal-grove-dry-100-xi', suite: 22, title: 'Deal–Grove dry oxidation of (100) Si, 1000 °C, 1 h, with the 25 nm initial-oxide correction', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'oxidation', ox_T: 1000, ox_t: 1, ox_amb: 'dry', ox_orient: '100', ox_x0: 25 }, kpi: 'xox', expected: X_DRY_100_25 * 1e3, tol: 0.005, kind: 'analytical',
    reference: 'Deal–Grove law with τ = (x_i² + A x_i)/B for x_i = 25 nm (Jaeger, Introduction to Microelectronic Fabrication)', note: 'About 58 nm, the figure read from dry-oxidation charts.' },
  { id: 's22-deal-grove-table-dry-111', suite: 22, title: 'Dry oxidation of (111) Si, 1000 °C, 1 h, against Deal and Grove’s tabulated constants', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'oxidation', ox_T: 1000, ox_t: 1, ox_amb: 'dry', ox_orient: '111', ox_x0: XI_DRY_111 * 1e3 }, kpi: 'xox',
    expected: dgThickness(DG_DRY_1000.A, DG_DRY_1000.B, 1, DG_DRY_1000.tau) * 1e3, tol: 0.05, kind: 'reference-data',
    reference: 'Deal and Grove, General relationship for the thermal oxidation of silicon — measured constants for dry O₂ at 1000 °C (A = 0.165 µm, B = 0.0117 µm²/h, τ = 0.37 h)',
    note: 'The suite uses Arrhenius fits of the rate constants, which differ from the individual tabulated values by up to about 10 %; the thickness agrees more closely.' },
  { id: 's22-deal-grove-table-wet-111', suite: 22, title: 'Wet oxidation of (111) Si, 1100 °C, 1 h, against Deal and Grove’s tabulated constants', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'oxidation', ox_T: 1100, ox_t: 1, ox_amb: 'wet', ox_orient: '111', ox_x0: 0 }, kpi: 'xox',
    expected: dgThickness(DG_WET_1100.A, DG_WET_1100.B, 1, 0) * 1e3, tol: 0.05, kind: 'reference-data',
    reference: 'Deal and Grove, General relationship for the thermal oxidation of silicon — measured constants for wet O₂ at 1100 °C (A = 0.11 µm, B = 0.510 µm²/h)',
    note: 'Arrhenius-fitted constants in the suite against the tabulated values at one temperature.' },
  { id: 's22-silicon-consumed', suite: 22, title: 'Silicon consumed per unit oxide grown', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'oxidation', ox_T: 1000, ox_t: 1, ox_amb: 'dry', ox_orient: '100', ox_x0: 0 }, kpi: 'si', expected: VM_RATIO * X_DRY_100 * 1e3, tol: 0.01, kind: 'reference-data',
    reference: 'Molar volumes of Si (28.086 g/mol, 2.329 g/cm³) and thermal SiO₂ (60.084 g/mol, 2.20 g/cm³)', note: 'The ratio is 0.4415; the suite uses the rounded textbook value 0.44.' },

  // ------------------------------------------------------------------ degradation: corrosion
  { id: 's22-faraday-iron', suite: 22, title: 'Faraday corrosion rate of iron at 1 µA/cm²', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'corrosion', co_metal: 'Fe', co_i: 0.01 }, kpi: 'rate', expected: 0.0116, tol: 0.01, kind: 'reference-data',
    reference: 'ASTM G102 conversion of corrosion current density to penetration rate (Fe, n = 2)', note: 'Handbook value quoted to three figures.' },
  { id: 's22-faraday-iron-mass-loss', suite: 22, title: 'Mass-loss rate of iron at 1 µA/cm²', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'corrosion', co_metal: 'Fe', co_i: 0.01 }, kpi: 'ml', expected: 0.01 * 55.845 / (2 * FARADAY) * 86400, tol: 1e-4, kind: 'analytical',
    reference: 'Faraday’s law, ṁ = iM/(nF)', note: '0.250 g/m²/day.' },
  { id: 's22-faraday-zinc', suite: 22, title: 'Faraday corrosion rate of zinc at 1 µA/cm²', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'corrosion', co_metal: 'Zn', co_i: 0.01 }, kpi: 'rate', expected: g102(65.38 / 2, 7.14), tol: 0.005, kind: 'reference-data',
    reference: 'ASTM G102 constant K = 3.27×10⁻³ mm·g/(µA·cm·yr) with the equivalent weight and density of zinc', note: 'About 0.0150 mm/yr.' },
  { id: 's22-faraday-aluminium', suite: 22, title: 'Faraday corrosion rate of aluminium at 10 µA/cm² (three-electron dissolution)', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'corrosion', co_metal: 'Al', co_i: 0.1 }, kpi: 'rate', expected: 10 * g102(26.982 / 3, 2.70), tol: 0.005, kind: 'reference-data',
    reference: 'ASTM G102 constant with the equivalent weight (M/3) and density of aluminium' },

  // ------------------------------------------------------------------ degradation: dissolution, ripening
  { id: 's22-dissolution-sink-lifetime', suite: 22, title: 'Shrinking-particle lifetime of 30 nm ZnO under sink conditions', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'dissolution', di_mat: 'ZnO', di_medium: 'water', di_d: 3e-8, di_k: 2e-8 }, kpi: 'tsink', expected: T_SINK_H, tol: 0.002, kind: 'analytical',
    reference: 'Surface-reaction-controlled shrinking particle, t = r₀ρ/(kM) (Levenspiel, Chemical Reaction Engineering)', note: 'ZnO density 5606 kg/m³ and standard atomic weights.' },
  { id: 's22-dissolution-cubic-law', suite: 22, title: 'Shrinking-particle conversion at half the sink lifetime', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'dissolution', di_mat: 'ZnO', di_medium: 'water', di_d: 3e-8, di_k: 2e-8, di_conc: 0.001, di_cs: 1e4, di_t: T_SINK_H / 2 }, kpi: 'frac', expected: 87.5, tol: 0.003, kind: 'analytical',
    reference: 'Shrinking-particle model with reaction control, X = 1 − (1 − t/τ)³ (Levenspiel, Chemical Reaction Engineering)',
    note: 'Very dilute loading and high solubility remove the saturation term (c/c_s below 10⁻⁸); the tolerance covers the density used to set the exposure time.' },
  { id: 's22-lsw-volume-doubling', suite: 22, title: 'LSW coarsening: diameter when the mean particle volume has doubled', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'ripening', rp_d: 20, rp_T: 60, rp_gamma: 0.3, rp_cinf: 1e-7, rp_D: 1e-9, rp_Vm: 14.5, rp_t: T_DOUBLE_H }, kpi: 'd', expected: 20 * Math.cbrt(2), tol: 1e-4, kind: 'analytical',
    reference: 'Lifshitz–Slyozov–Wagner theory, r̄³ − r̄₀³ = 8γc∞V_m²Dt/(9RT)', note: 'The ageing time is r̄₀³/K, evaluated here with the CODATA gas constant.' },
  { id: 's22-smoluchowski-half-time', suite: 22, title: 'Smoluchowski rapid-coagulation half-time in water at 20 °C', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'ripening', rp_T: 20, rp_n0: 1e18, rp_W: 1 }, kpi: 'tauA', expected: 3 * 1.002e-3 / (4 * KB * K(20) * 1e18) / 3600, tol: 0.01, kind: 'reference-data',
    reference: 'Smoluchowski t½ = 3μ/(4k_BT n₀) with the handbook viscosity of water at 20 °C (1.002 mPa·s)', note: 'About 0.19 s for 10¹⁸ particles/m³; the suite uses a viscosity correlation.' },

  // ------------------------------------------------------------------ degradation: delamination
  { id: 's22-film-debond-energy', suite: 22, title: 'Steady-state energy release rate for film debonding', mode: 'degradation',
    inputs: { ...FILM, dl_sigma: 1000 }, kpi: 'G', expected: 1e18 * 1e-6 / (2 * E_PS), tol: 1e-6, kind: 'analytical',
    reference: 'Hutchinson and Suo, Mixed mode cracking in layered materials — G = σ²h(1 − ν²)/(2E)', note: '1 GPa in a 1 µm film of E = 150 GPa, ν = 0.25: 3.125 J/m².' },
  { id: 's22-film-critical-thickness', suite: 22, title: 'Critical film thickness for delamination', mode: 'degradation',
    inputs: { ...FILM, dl_sigma: 1000, dl_Gi: 5 }, kpi: 'hc', expected: 2 * E_PS * 5 / 1e18 * 1e6, tol: 1e-6, kind: 'analytical',
    reference: 'Hutchinson and Suo, Mixed mode cracking in layered materials — h_c = 2E′Γ_i/σ²' },
  { id: 's22-film-channel-cracking', suite: 22, title: 'Critical thickness for channel cracking of a tensile film', mode: 'degradation',
    inputs: { ...FILM, dl_sigma: 1000, dl_Gf: 10 }, kpi: 'hcc', expected: 10 * E_PS / (1.976 * 1e18) * 1e6, tol: 1e-6, kind: 'analytical',
    reference: 'Beuth, Cracking of thin bonded films in residual tension — G = 1.976σ²h/E′ without elastic mismatch' },
  { id: 's22-film-buckling-stress', suite: 22, title: 'Buckling stress of a clamped debonded strip', mode: 'degradation',
    inputs: { ...FILM, dl_sigma: -1000, dl_b: 20 }, kpi: 'sb', expected: PI * PI / 12 * E_PS * Math.pow(1 / 20, 2) / 1e6, tol: 1e-6, kind: 'analytical',
    reference: 'Hutchinson and Suo, Mixed mode cracking in layered materials — σ_c = (π²/12)E′(h/b)²' },
  { id: 's22-film-thermal-mismatch', suite: 22, title: 'Thermal mismatch stress of a film cooled on a lower-expansion substrate', mode: 'degradation',
    inputs: { ...FILM, dl_sigma: 0, dl_dalpha: 5e-6, dl_dT: -100 }, kpi: 'sigma', expected: 150e9 / 0.75 * 5e-6 * 100 / 1e6, tol: 1e-6, kind: 'analytical',
    reference: 'Biaxial thermal mismatch stress σ = [E/(1 − ν)](α_substrate − α_film)ΔT (Freund and Suresh, Thin Film Materials)',
    note: 'A film that contracts more than its substrate on cooling is left in tension: +100 MPa.' },

  // ------------------------------------------------------------------ degradation: electromigration, dielectric breakdown
  { id: 's22-black-current-scaling', suite: 22, title: 'Black’s equation: lifetime ratio at half the test current density', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'electromigration', em_j: 1.25e10, em_jref: 2.5e10, em_n: 2, em_T: 300, em_Tref: 300, em_tref: 300 }, kpi: 't50', expected: 4 * 300 / YEAR_H, tol: 1e-6, kind: 'analytical',
    reference: 'Black’s equation MTTF = A j⁻ⁿ exp(E_a/kT)', note: 'Halving the current density with n = 2 quadruples the median life.' },
  { id: 's22-black-thermal', suite: 22, title: 'Black’s equation: thermal acceleration between 300 °C test and 105 °C use', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'electromigration', em_T: 105, em_Tref: 300, em_Ea: 0.9 }, kpi: 'afT', expected: Math.exp(0.9 / KB_EV * (1 / K(105) - 1 / K(300))), tol: 1e-4, kind: 'analytical',
    reference: 'Black’s equation with the CODATA Boltzmann constant' },
  { id: 's22-lognormal-early-failure', suite: 22, title: 'Lognormal time to 0.1 % cumulative failures', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'electromigration', em_j: 2.5e10, em_jref: 2.5e10, em_T: 300, em_Tref: 300, em_tref: 300, em_sigma: 0.5 }, kpi: 't01', expected: 300 * Math.exp(0.5 * -3.090232) / YEAR_H, tol: 1e-4, kind: 'analytical',
    reference: 'Lognormal distribution, t_p = t₅₀exp(σz_p) with the standard normal quantile z₀.₀₀₁ = −3.0902' },
  { id: 's22-tddb-field', suite: 22, title: 'TDDB E-model: characteristic life 1 MV/cm below the test field', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'tddb', td_V: 9, td_tox: 10, td_T: 125, td_Tref: 125, td_Eref: 10, td_gamma: 2.5, td_tref: 100, td_area: 1 }, kpi: 't63', expected: 100 * Math.exp(2.5) / YEAR_H, tol: 1e-6, kind: 'analytical',
    reference: 'Thermochemical E-model t_BD ∝ exp(−γE_ox) (McPherson)' },
  { id: 's22-tddb-weibull-area', suite: 22, title: 'TDDB: weakest-link area scaling and Weibull early-failure time', mode: 'degradation',
    inputs: { mode: 'degradation', deg_mech: 'tddb', td_V: 10, td_tox: 10, td_T: 125, td_Tref: 125, td_Eref: 10, td_tref: 100, td_area: 100, td_beta: 2 }, kpi: 'tp', expected: 100 * Math.pow(1 / 100, 1 / 2) * Math.sqrt(-Math.log(1 - 1e-4)) / YEAR_H, tol: 1e-6, kind: 'analytical',
    reference: 'Weibull statistics: t₆₃ ∝ A^(−1/β), t_p = t₆₃[−ln(1 − p)]^(1/β)' },

  // ------------------------------------------------------------------ functional lifetime of metamaterials
  { id: 's22-detuning-thickness-loss', suite: 22, title: 'Resonance detuning by surface recession alone: time to leave a ±3 % band', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'detuning', det_band: 3, det_t0: 100, det_rate: 20, det_dE: 0, det_dcrit: 0 }, kpi: 'life', expected: 0.03 * 100e-6 / (2 * 20e-9), tol: 1e-6, kind: 'analytical',
    reference: 'Bending resonator, f ∝ thickness: f/f₀ = 1 − 2rt/t₀', note: '20 nm/yr per face on a 100 µm flexure: 75 years.' },
  { id: 's22-bistable-snap-force', suite: 22, title: 'Snap-through force of a cosine-shaped clamped beam', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'bistable', bi_l: 10, bi_t: 150, bi_h: 600, bi_w: 5000, bi_E: 2 }, kpi: 'Fsnap', expected: 740 * 2e9 * (5e-3 * Math.pow(150e-6, 3) / 12) * 600e-6 / 1e-6 * 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'Qiu, Lang and Slocum, A curved-beam bistable mechanism — F_top ≈ 740EIh/l³', note: 'Asymptotic formula for h/t well above 2.31; the benchmark verifies its evaluation, not the approximation.' },
  { id: 's22-switching-drive-scaling', suite: 22, title: 'Switching endurance: inverse-power scaling with drive amplitude', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'switching', sw_tech: 'pcm', sw_drive: 2, sw_T: 25, sw_scale: 1 }, kpi: 'eta', expected: 1e9 / 16, tol: 1e-6, kind: 'analytical',
    reference: 'Inverse-power-law life–stress model η ∝ V^(−p) (Nelson, Accelerated Testing)',
    note: 'Checks the scaling only: the baseline of 10⁹ cycles and exponent 4 for phase-change elements are the suite’s stated typical values, not reference data.' },
  { id: 's22-auxetic-kinematic-range', suite: 22, title: 'Rotating squares: tensile strain available from a 30° opening angle', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'auxetic', ax_theta: 30 }, kpi: 'range', expected: (2 / Math.sqrt(3) - 1) * 100, tol: 1e-6, kind: 'analytical',
    reference: 'Grima and Evans, Auxetic behavior from rotating squares — L ∝ cos(θ/2) + sin(θ/2)', note: '√2/(cos 15° + sin 15°) − 1 = 2/√3 − 1 = 15.47 %.' },
  { id: 's22-auxetic-hinge-rotation', suite: 22, title: 'Rotating squares: hinge rotation from 30° to 60°', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'auxetic', ax_theta: 30, ax_strain: (S(60) / S(30) - 1) * 100 }, kpi: 'dtheta', expected: 30, tol: 1e-6, kind: 'analytical',
    reference: 'Grima and Evans, Auxetic behavior from rotating squares — strain ε = L(θ)/L(θ₀) − 1 inverted for θ' },
  { id: 's22-auxetic-poisson', suite: 22, title: 'Rotating rigid squares: Poisson’s ratio', mode: 'functional',
    inputs: { mode: 'functional', fun_mech: 'auxetic' }, pub: 'eff.nu', expected: -1, abs: 1e-12, kind: 'analytical',
    reference: 'Grima and Evans, Auxetic behavior from rotating squares — ν = −1 for rigid units' },

  // ------------------------------------------------------------------ reliability synthesis
  { id: 's22-weibull-shape-recovery', suite: 22, title: 'Median-rank regression recovers the shape of an exact Weibull sample', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: WB_TABLE }, kpi: 'beta', expected: WB.beta, tol: 1e-6, kind: 'analytical',
    reference: 'Weibull probability plotting with Benard’s median-rank approximation (Abernethy, The New Weibull Handbook)', note: 'Twelve failure times placed exactly on the β = 2, η = 1000 h line.' },
  { id: 's22-weibull-b10', suite: 22, title: 'Weibull B10 life, β = 2, η = 1000 h', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: WB_TABLE }, kpi: 'B10', expected: WB.eta * Math.sqrt(-Math.log(0.9)), tol: 1e-5, kind: 'analytical',
    reference: 'Two-parameter Weibull distribution, B10 = η(−ln 0.9)^(1/β)', note: '324.6 h.' },
  { id: 's22-weibull-mttf', suite: 22, title: 'Weibull mean life, β = 2 (Rayleigh)', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: WB_TABLE }, kpi: 'mttf', expected: WB.eta * SQRT_PI_2, tol: 1e-5, kind: 'analytical',
    reference: 'Two-parameter Weibull distribution, MTTF = ηΓ(1 + 1/β) with Γ(3/2) = √π/2' },
  { id: 's22-weibull-target-reliability', suite: 22, title: 'Age at 95 % reliability, β = 2, η = 1000 h', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: WB_TABLE, rel_Rt: 0.95 }, kpi: 'tR', expected: WB.eta * Math.sqrt(-Math.log(0.95)), tol: 1e-5, kind: 'analytical',
    reference: 'Two-parameter Weibull distribution, t_R = η(−ln R)^(1/β)' },
  { id: 's22-rul-exponential', suite: 22, title: 'Remaining useful life from an exact exponential degradation path', mode: 'reliability',
    inputs: { mode: 'reliability', rel_deg: DEG_EXP, rel_trend: 'exp', rel_thr: 50 }, kpi: 'RUL', expected: Math.LN2 / 1e-4 - 5000, tol: 1e-6, kind: 'analytical',
    reference: 'Exponential decay y = y₀e^(−λt): threshold y₀/2 is reached at ln 2/λ', note: 'λ = 10⁻⁴ 1/h observed to 5000 h: 1931.5 h remain.' },
  { id: 's22-rul-linear', suite: 22, title: 'Remaining useful life from an exact linear degradation path', mode: 'reliability',
    inputs: { mode: 'reliability', rel_deg: DEG_LIN, rel_trend: 'lin', rel_thr: 40 }, kpi: 'RULlo', expected: 2000, tol: 1e-6, kind: 'analytical',
    reference: 'Linear trend extrapolation; with zero residual scatter the lower prediction bound equals the point estimate' },

  // ------------------------------------------------------------------ published experiments
  { id: 's22-exp-lieblein-zelen-bearing-weibull-slope', suite: 22, title: 'Experiment: Weibull slope of the endurance lives of 23 deep-groove ball bearings', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: LZ_TABLE }, kpi: 'beta', expected: 2.23, tol: 0.05, kind: 'experiment',
    measured: { value: 2.23, uncertainty: null, unit: '–' },
    source: LZ_SOURCE, conditions: LZ_COND, match: ['rel_fail'],
    reference: 'Lieblein & Zelen, J. Res. Natl. Bur. Stand. 57, 273 (1956), appendix A, specimen worksheet',
    note: 'The 23 endurance lives are transcribed from the specimen worksheet on p. 286 (they are in millions of revolutions; the suite labels the time axis in hours, which does not affect the shape). The reference value is the "Slope of Curve 2.23" entered on the same worksheet by the testing laboratory from its plot on Weibull-function coordinate paper. Median-rank regression in the suite gives 2.18 (−2.2 %). Tolerance 5 %: the laboratory value is a graphical fit, and the three bearings marked "Disc." are treated as failures at their recorded lives (the suite has no censored-data method). The source notes that the worksheet estimates were not part of the data it analysed.' },
  { id: 's22-exp-lieblein-zelen-bearing-b10', suite: 22, title: 'Experiment: B10 life of 23 deep-groove ball bearings', mode: 'reliability',
    inputs: { mode: 'reliability', rel_fail: LZ_TABLE }, kpi: 'B10', expected: 29, tol: 0.1, kind: 'experiment',
    measured: { value: 29, uncertainty: null, unit: '10⁶ revolutions' },
    source: LZ_SOURCE, conditions: LZ_COND, match: ['rel_fail'],
    reference: 'Lieblein & Zelen, J. Res. Natl. Bur. Stand. 57, 273 (1956), appendix A, specimen worksheet',
    note: 'Reference value "B-10 29." as entered on the specimen worksheet (p. 286) by the testing laboratory from its Weibull plot, quoted to two figures; data as in the previous entry, in millions of revolutions. The suite reports the B10 life of the better-fitting of its Weibull and lognormal regressions — here the lognormal — and obtains 31.2 (+7.5 %); its Weibull parameters (β = 2.18, η = 81.6) correspond to 29.1. Tolerance 10 % for the choice of distribution, the graphical reference value and the treatment of the three discontinued bearings as failures.' },
];
