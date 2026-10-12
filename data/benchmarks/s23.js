// Reference benchmarks for Suite 23 — Experimental characterisation and metrology.
// Every data table below is generated in this file from a closed-form forward model with known parameters
// (no noise, no call into the suite); the reduction must return the generating parameters.

const PI = Math.PI, D2R = PI / 180, LN2 = Math.LN2;
const K_B = 1.380649e-23, N_A = 6.02214076e23, Q_E = 1.602176634e-19, MU0 = 4e-7 * PI;
const csv = (head, rows) => head.join(',') + '\n' + rows.map(r => r.join(',')).join('\n');
const lin = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
const logs = (a, b, n) => lin(Math.log10(a), Math.log10(b), n).map(v => 10 ** v);

// ---------- X-ray diffraction: fcc copper, Cu Kα₁ ----------
const LAM = 1.5406, A_CU = 3.6149, K_SCH = 0.9;
const CU_HKL = [[1, 1, 1, 100], [2, 0, 0, 46], [2, 2, 0, 20], [3, 1, 1, 17]];          // relative integrated intensities of Cu powder
const cuTheta = ([h, k, l]) => Math.asin(LAM * Math.sqrt(h * h + k * k + l * l) / (2 * A_CU));
const gauss = (x, x0, w) => Math.exp(-4 * LN2 * ((x - x0) / w) ** 2), lorentz = (x, x0, w) => 1 / (1 + 4 * ((x - x0) / w) ** 2);
function cuPattern(widthDeg, shape, areaNorm) {                                          // widthDeg(θ) → FWHM in degrees 2θ
  return csv(['two_theta_deg', 'intensity'], lin(35, 93, 2901).map(tt => {
    let y = 40;
    for (const r of CU_HKL) { const th = cuTheta(r), w = widthDeg(th); y += 4000 * r[3] / 100 * (0.3 / w) / areaNorm * shape(tt, 2 * th / D2R, w); }
    return [tt.toFixed(2), y.toFixed(3)];
  }));
}
const D_SIZE = 30, INST = 0.08;                                                          // nm, degrees
const XRD_SIZE = cuPattern(th => Math.hypot(K_SCH * LAM / (D_SIZE * 10 * Math.cos(th)) / D2R, INST), gauss, 1);
const D_WH = 40, EPS_WH = 2e-3;                                                          // Lorentzian breadths add linearly
const XRD_WH = cuPattern(th => (K_SCH * LAM / (D_WH * 10 * Math.cos(th)) + 4 * EPS_WH * Math.tan(th)) / D2R, lorentz, 1);
const RHO_CU = 4 * 63.546e-3 / (N_A * (A_CU * 1e-10) ** 3);                              // 4 atoms per fcc cell

// ---------- dynamic light scattering: monodisperse 100 nm spheres in water at 25 °C ----------
const D_SE = K_B * 298.15 / (3 * PI * 0.89e-3 * 100e-9);                                 // Stokes–Einstein, 4.91e-12 m²/s
const Q_DLS = 4 * PI * 1.33 / 633e-9 * Math.sin(173 / 2 * D2R);
const DLS = csv(['tau_s', 'g2'], logs(1e-6, 1e-1, 80).map(t => [t.toExponential(6), (1 + 0.9 * Math.exp(-2 * D_SE * Q_DLS * Q_DLS * t)).toFixed(8)]));

const G_POLY = K_B * 298.15 / (3 * PI * 0.89e-3 * 50e-9) * Q_DLS * Q_DLS;                 // 50 nm z-average, PDI 0.1 (second cumulant µ₂ = 0.1·Γ²)
const DLS_POLY = csv(['tau_s', 'g2'], logs(1e-6, 1e-1, 80).map(t => { const g1 = Math.exp(-G_POLY * t) * (1 + 0.1 * G_POLY * G_POLY * t * t / 2); return [t.toExponential(6), (1 + 0.9 * g1 * g1).toFixed(8)]; }));

// ---------- counted diameters ----------
const TEM = csv(['diameter_nm'], [[10], [20], [30], [40], [50]]);

// ---------- BET isotherm (infinite-layer BET equation) ----------
const VM = 20, C_BET = 100;
const BET = csv(['p_rel', 'v_ads_cm3_g'], [0.01, 0.02, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7].map(x => [x, (VM * C_BET * x / ((1 - x) * (1 - x + C_BET * x))).toFixed(6)]));

const V_LAST = VM * C_BET * 0.7 / (0.3 * (0.3 + C_BET * 0.7));                           // uptake at the highest point, p/p₀ = 0.7
// non-porous solid whose uptake is a uniform film of Harkins–Jura thickness on 50 m²/g
const TPLOT = csv(['p_rel', 'v_ads_cm3_g'], lin(0.05, 0.9, 35).map(x => [x.toFixed(4), (50 / 1.547 * 0.1 * Math.sqrt(13.99 / (0.034 - Math.log10(x)))).toFixed(6)]));

// ---------- Tauc: ideal parabolic-band absorbers ----------
const tauc = (Eg, r) => csv(['wavelength_nm', 'absorbance'], lin(300, 520, 221).map(l => { const E = 1239.842 / l; return [l, (0.8 * Math.pow(Math.max(0, E - Eg), r) / E).toFixed(8)]; }));
const UV_DIRECT = tauc(3.30, 0.5), UV_INDIRECT = tauc(2.60, 2);

// ---------- LSPR: Lorentzian plasmon band on a flat interband background ----------
const LSPR_BG = 0.30, LSPR_PK = 0.30, LSPR_L0 = 520, LSPR_W = 45;                        // half-width 45 nm
const lsprA = l => LSPR_BG + LSPR_PK / (1 + ((l - LSPR_L0) / LSPR_W) ** 2);
const LSPR = csv(['wavelength_nm', 'absorbance'], lin(400, 800, 401).map(l => [l, lsprA(l).toFixed(8)]));

// ---------- Raman: defective graphitic carbon on strained silicon ----------
const RAMAN = csv(['shift_cm1', 'intensity'], lin(100, 3000, 2901).map(s => [s, (50 + 1000 * lorentz(s, 518.7, 5) + 300 * lorentz(s, 1350, 40) + 600 * lorentz(s, 1582, 22) + 450 * lorentz(s, 2690, 50)).toFixed(4)]));

// ---------- Oliver–Pharr: fused silica, ideal Berkovich tip ----------
const OP = (() => {
  const E = 72e9, nu = 0.17, H = 9e9, Ei = 1141e9, nui = 0.07, P = 10e-3, m = 1.3;
  const Er = 1 / ((1 - nu * nu) / E + (1 - nui * nui) / Ei), A = P / H, hc = Math.sqrt(A / 24.5), S = 2 * 1.034 / Math.sqrt(PI) * Er * Math.sqrt(A), hmax = hc + 0.75 * P / S, hf = hmax - m * P / S, rows = [];
  for (let i = 1; i <= 50; i++) { const h = hmax * i / 50; rows.push([(h * 1e9).toFixed(5), (10 * (h / hmax) ** 2).toFixed(7)]); }
  for (let i = 1; i <= 60; i++) { const h = hmax - (hmax - hf) * i / 61; rows.push([(h * 1e9).toFixed(5), (10 * ((h - hf) / (hmax - hf)) ** m).toFixed(7)]); }
  return { table: csv(['depth_nm', 'load_mN'], rows), S, hc: hc * 1e9 };
})();

// ---------- Hertz: rigid sphere R = 1 µm on a 5 MPa, ν = 0.4 half-space ----------
const HERTZ = csv(['indentation_nm', 'force_nN'], lin(-20, 80, 51).map(d => [d, (4 / 3 * 5e6 / (1 - 0.16) * Math.sqrt(1e-6) * Math.pow(Math.max(0, d * 1e-9), 1.5) * 1e9).toFixed(6)]));

// ---------- electrical ----------
const FPP = csv(['current_mA', 'voltage_mV'], lin(-1, 1, 11).map(i => [i.toFixed(2), (2 * i).toFixed(6)]));            // V/I = 2 Ω
const HALL_N = 1e25, HALL_T = 1e-6, HALL_I = 1e-3;
const HALL = csv(['B_T', 'V_hall_mV'], lin(-1, 1, 11).map(b => [b.toFixed(2), (-HALL_I * b / (Q_E * HALL_N * HALL_T) * 1e3 + 0.05).toFixed(8)]));
// Randles circuit with an ideal capacitor and a Warburg element, complex arithmetic written out
const EIS = csv(['f_Hz', 'Z_re_ohm', 'Z_im_ohm'], logs(1e5, 1e-2, 60).map(f => {
  const w = 2 * PI * f, Rs = 10, Rct = 100, Cdl = 2e-5, sw = 15, fr = Rct + sw / Math.sqrt(w), fi = -sw / Math.sqrt(w), d = fr * fr + fi * fi;
  const yr = fr / d, yi = -fi / d + w * Cdl, D = yr * yr + yi * yi;
  return [f.toExponential(6), (Rs + yr / D).toFixed(8), (-yi / D).toFixed(8)];
}));

// ---------- laser flash: adiabatic Parker rear-face solution, α = 1.0e-5 m²/s through 2 mm ----------
const FLASH = csv(['time_ms', 'signal'], lin(0, 400, 401).map(t => { const w = PI * PI * 1e-5 * t * 1e-3 / 4e-6; let s = 1; for (let k = 1; k < 80; k++) s += 2 * (k % 2 ? -1 : 1) * Math.exp(-k * k * w); return [t, (w < 0.04 ? 0 : s).toFixed(8)]; }));

// ---------- tensile: Ramberg–Osgood curves ε = σ/E + 0.002(σ/σ_y)^n (σ_y is the 0.2 % proof stress by construction) ----------
const tens = (E, sy, n, smax, nu) => csv(['strain', 'stress_MPa', 'transverse_strain'], lin(0, smax, 120).map(s => { const e = s / E + 0.002 * (s / sy) ** n; return [e.toExponential(8), s.toFixed(6), (-nu * e).toExponential(8)]; }));
const TENS_STEEL = tens(200000, 400, 15, 440, 0.30), TENS_AUX = tens(50, 1.5, 10, 1.7, -0.6);

// ---------- VSM: two tanh branches shifted by ±H_c ----------
const VSM_MS = 400, VSM_HC = 20, VSM_H0 = 30;
const VSM = csv(['H_kA_m', 'M_kA_m'], [...lin(1000, -1000, 401).map(H => [H, (VSM_MS * Math.tanh((H + VSM_HC) / VSM_H0)).toFixed(8)]), ...lin(-1000, 1000, 401).map(H => [H, (VSM_MS * Math.tanh((H - VSM_HC) / VSM_H0)).toFixed(8)])]);

// ---------- metrology ----------
const GUM_FPP = 'quantity,value,u_standard\nV_mV,1,0.01\nI_mA,1,0.02\nt_nm,100,2';       // 1 %, 2 %, 2 % → 3 % combined
const GUM_STOKES = `quantity,value,u_standard\nT_K,298.15,0\neta_mPas,0.89,0\nD_m2_s,${D_SE},0`;
const REP = csv(['value'], [[9.8], [10.0], [10.2], [10.1], [9.9]]);

export default [
  // ---- XRD ----
  { id: 's23-xrd-bragg-cu111', suite: 23, title: 'XRD: Bragg d-spacing of Cu (111) from a simulated pattern', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_SIZE }, kpi: 'd_main', expected: A_CU / Math.sqrt(3), tol: 0.001, kind: 'analytical', reference: "Bragg's law, d = a/√(h²+k²+l²) for a cubic cell; Cullity, Elements of X-Ray Diffraction", note: 'Gaussian peaks placed at the Bragg angles of fcc Cu (a = 3.6149 Å, Cu Kα₁).' },
  { id: 's23-xrd-lattice-cu', suite: 23, title: 'XRD: refined lattice parameter of copper', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_SIZE }, kpi: 'a', expected: 3.6149, tol: 0.001, kind: 'reference-data', reference: 'Lattice parameter of fcc copper at room temperature, 3.615 Å (Cullity, Elements of X-Ray Diffraction, appendix of crystal structures)', note: 'Requires the phase ranking to select Cu from the built-in library and index four reflections.' },
  { id: 's23-xrd-scherrer-30nm', suite: 23, title: 'XRD: Scherrer size returned for a pattern built for 30 nm crystallites', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_SIZE, xrd_inst: INST, xrd_K: K_SCH, xrd_corr: 'quadrature' }, kpi: 'size', expected: D_SIZE, tol: 0.03, kind: 'analytical', reference: 'Scherrer equation D = Kλ/(β cos θ) with Gaussian quadrature removal of the instrumental breadth; Klug & Alexander, X-Ray Diffraction Procedures', note: 'Peak FWHM generated as √(β² + b²) with b = 0.08°; tolerance covers the pseudo-Voigt fit on a 0.02° grid.' },
  { id: 's23-xrd-wh-strain', suite: 23, title: 'XRD: Williamson–Hall microstrain of a pattern with known strain broadening', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_WH, xrd_inst: 0, xrd_K: K_SCH }, kpi: 'strain', expected: EPS_WH, tol: 0.02, kind: 'analytical', reference: 'Williamson–Hall relation β cos θ = Kλ/D + 4ε sin θ; Williamson & Hall, Acta Metallurgica', note: 'Lorentzian peaks whose size and strain breadths add linearly (D = 40 nm, ε = 0.002); Lorentzian tails partly absorbed by the background clipping.' },
  { id: 's23-xrd-wh-size', suite: 23, title: 'XRD: Williamson–Hall size of the same pattern', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_WH, xrd_inst: 0, xrd_K: K_SCH }, kpi: 'wh_size', expected: D_WH, tol: 0.02, kind: 'analytical', reference: 'Williamson–Hall relation β cos θ = Kλ/D + 4ε sin θ; Williamson & Hall, Acta Metallurgica', note: 'The intercept is an extrapolation over four peaks, which amplifies width errors; with noise-free Lorentzian input the pseudo-Voigt fit returns the widths to better than 0.1 %.' },
  { id: 's23-xrd-density-cu', suite: 23, title: 'XRD: X-ray density of copper', mode: 'xrd', inputs: { mode: 'xrd', xrd_table: XRD_SIZE }, kpi: 'rho_x', expected: RHO_CU, tol: 0.003, kind: 'analytical', reference: 'ρ = Z·M/(N_A·a³) with Z = 4, M = 63.546 g/mol, a = 3.6149 Å (≈ 8.93 g/cm³; handbook bulk density 8.96 g/cm³)' },
  { id: 's23-xrd-zno-c', suite: 23, title: 'XRD: c lattice parameter of wurtzite ZnO from the default noisy pattern', mode: 'xrd', inputs: { mode: 'xrd' }, kpi: 'c', expected: 5.2066, tol: 0.002, kind: 'reference-data', reference: 'Wurtzite ZnO, a = 3.2498 Å, c = 5.2066 Å (ICDD powder diffraction file for zincite)', note: 'Default demonstration pattern with counting noise and background; phase must be identified as wurtzite ZnO first.' },
  // ---- particle sizing ----
  { id: 's23-dls-stokes-einstein-D', suite: 23, title: 'DLS: diffusion coefficient of 100 nm spheres in water at 25 °C', mode: 'dls', inputs: { mode: 'dls', dls_table: DLS, dls_T: 25, dls_visc: 0.89, dls_n: 1.33, dls_lambda: 633, dls_angle: 173 }, kpi: 'D', expected: 4.91e-12, tol: 0.01, kind: 'analytical', reference: 'Stokes–Einstein relation D = k_B·T/(3π·η·d) with η = 0.89 mPa·s (water, 25 °C)', note: 'Single-exponential g₂ generated from Γ = D·q², q = 4πn·sin(θ/2)/λ.' },
  { id: 's23-dls-diameter', suite: 23, title: 'DLS: z-average diameter of a monodisperse 100 nm sol', mode: 'dls', inputs: { mode: 'dls', dls_table: DLS }, kpi: 'd_z', expected: 100, tol: 0.01, kind: 'analytical', reference: 'Cumulant analysis of a single-exponential correlation function (Koppel); Stokes–Einstein relation' },
  { id: 's23-dls-pdi-zero', suite: 23, title: 'DLS: polydispersity index of a monodisperse sol', mode: 'dls', inputs: { mode: 'dls', dls_table: DLS }, kpi: 'pdi', expected: 0, abs: 0.01, kind: 'analytical', reference: 'Second cumulant of a single exponential is zero (Koppel cumulant expansion)' },
  { id: 's23-dls-pdi-cumulant', suite: 23, title: 'DLS: polydispersity index of a second-order cumulant correlation function (PDI = 0.1)', mode: 'dls', inputs: { mode: 'dls', dls_table: DLS_POLY }, kpi: 'pdi', expected: 0.1, tol: 0.02, kind: 'analytical', reference: 'Cumulant expansion g₁ = exp(−Γτ)·(1 + µ₂τ²/2), PDI = µ₂/Γ² (Frisken, Applied Optics)' },
  { id: 's23-counts-required-n', suite: 23, title: 'Counting statistics: particles needed for ±2 % on the mean (CV = 0.527)', mode: 'counts', inputs: { mode: 'counts', tem_table: TEM, tem_target: 2 }, kpi: 'need', expected: Math.ceil((1.96 * Math.sqrt(250) / 30 / 0.02) ** 2), tol: 1e-9, kind: 'analytical', reference: 'Normal-theory sample size n = (z₀.₉₇₅·CV/ε)², z = 1.96' },
  { id: 's23-counts-sauter', suite: 23, title: 'Counting statistics: Sauter mean D[3,2] of {10, 20, 30, 40, 50} nm', mode: 'counts', inputs: { mode: 'counts', tem_table: TEM }, kpi: 'd32', expected: 225000 / 5500, tol: 1e-6, kind: 'analytical', reference: 'Moment-ratio definition D[3,2] = Σd³/Σd² (Allen, Particle Size Measurement)' },
  { id: 's23-counts-d43', suite: 23, title: 'Counting statistics: volume-weighted mean D[4,3] of the same set', mode: 'counts', inputs: { mode: 'counts', tem_table: TEM }, kpi: 'd43', expected: 9790000 / 225000, tol: 1e-6, kind: 'analytical', reference: 'Moment-ratio definition D[4,3] = Σd⁴/Σd³ (Allen, Particle Size Measurement)' },
  { id: 's23-counts-sd', suite: 23, title: 'Counting statistics: sample standard deviation of the same set', mode: 'counts', inputs: { mode: 'counts', tem_table: TEM }, kpi: 'sd', expected: Math.sqrt(250), tol: 1e-6, kind: 'analytical', reference: 'Sample standard deviation with n − 1 divisor' },
  // ---- gas adsorption ----
  { id: 's23-bet-ssa', suite: 23, title: 'BET: specific surface area from an exact BET isotherm (v_m = 20 cm³/g)', mode: 'bet', inputs: { mode: 'bet', bet_table: BET, bet_sigma: 0.162 }, kpi: 'ssa', expected: VM * 4.353, tol: 0.005, kind: 'analytical', reference: 'Brunauer–Emmett–Teller equation; 1 cm³(STP) of N₂ monolayer covers 4.35 m² at 0.162 nm² per molecule (Gregg & Sing, Adsorption, Surface Area and Porosity)' },
  { id: 's23-bet-C', suite: 23, title: 'BET: C constant recovered from the same isotherm', mode: 'bet', inputs: { mode: 'bet', bet_table: BET }, kpi: 'C', expected: C_BET, tol: 0.01, kind: 'analytical', reference: 'Linear BET transform: C = slope/intercept + 1 (Gregg & Sing, Adsorption, Surface Area and Porosity)' },
  { id: 's23-bet-diameter', suite: 23, title: 'BET: equivalent sphere diameter d = 6/(ρ·S)', mode: 'bet', inputs: { mode: 'bet', bet_table: BET, rho: 4000 }, kpi: 'd_eq', expected: 6 / (4000 * VM * 4.353 * 1e3) * 1e9, tol: 0.005, kind: 'analytical', reference: 'Surface-to-volume ratio of a sphere, S = 6/(ρ·d)' },
  { id: 's23-bet-gurvich-pore-volume', suite: 23, title: 'Gas adsorption: Gurvich total pore volume from the last isotherm point', mode: 'bet', inputs: { mode: 'bet', bet_table: BET }, kpi: 'v_pore', expected: V_LAST * 1.547e-3, tol: 0.002, kind: 'analytical', reference: 'Gurvich rule with the liquid-nitrogen conversion 1 cm³(STP) = 1.547×10⁻³ cm³ liquid (Gregg & Sing, Adsorption, Surface Area and Porosity)' },
  { id: 's23-bet-tplot-area', suite: 23, title: 't-plot: external area of a non-porous solid carrying a Harkins–Jura film (50 m²/g)', mode: 'bet', inputs: { mode: 'bet', bet_table: TPLOT }, kpi: 's_ext', expected: 50, tol: 0.005, kind: 'analytical', reference: 'de Boer t-plot: S = 15.47·V/t (m²/g, V in cm³ STP/g, t in Å) with the Harkins–Jura thickness t = [13.99/(0.034 − log₁₀(p/p₀))]^½ Å (Gregg & Sing; Lowell et al., Characterization of Porous Solids and Powders)' },
  // ---- optical ----
  { id: 's23-tauc-direct', suite: 23, title: 'Tauc: direct gap of an ideal (hν − E_g)^½ absorber, E_g = 3.30 eV', mode: 'tauc', inputs: { mode: 'tauc', uv_table: UV_DIRECT, tauc_type: 'direct', tauc_base: 'subtract' }, kpi: 'Eg', expected: 3.30, tol: 0.005, kind: 'analytical', reference: 'Tauc relation (αhν)² ∝ hν − E_g for direct allowed transitions; Pankove, Optical Processes in Semiconductors' },
  { id: 's23-tauc-indirect', suite: 23, title: 'Tauc: indirect gap of an ideal (hν − E_g)² absorber, E_g = 2.60 eV', mode: 'tauc', inputs: { mode: 'tauc', uv_table: UV_INDIRECT, tauc_type: 'indirect', tauc_base: 'subtract' }, kpi: 'Eg', expected: 2.60, tol: 0.005, kind: 'analytical', reference: 'Tauc relation (αhν)^½ ∝ hν − E_g for indirect allowed transitions; Pankove, Optical Processes in Semiconductors' },
  { id: 's23-tauc-edge-wavelength', suite: 23, title: 'Tauc: absorption-edge wavelength λ = hc/E_g', mode: 'tauc', inputs: { mode: 'tauc', uv_table: UV_DIRECT, tauc_type: 'direct' }, kpi: 'lambda_g', expected: 1239.842 / 3.30, tol: 0.005, kind: 'analytical', reference: 'Photon energy–wavelength relation, hc = 1239.84 eV·nm' },
  { id: 's23-lspr-beer-lambert', suite: 23, title: 'LSPR: Beer–Lambert particle concentration', mode: 'lspr', inputs: { mode: 'lspr', lspr_table: LSPR, lspr_metal: 'Au', lspr_eps: 1e9, lspr_path: 0.5 }, kpi: 'conc', expected: (LSPR_BG + LSPR_PK) / (1e9 * 0.5), tol: 0.01, kind: 'analytical', reference: 'Beer–Lambert law A = ε·c·l', note: 'Peak absorbance 0.60 of a Lorentzian band; the five-point smoothing lowers the located maximum by well under 1 %.' },
  { id: 's23-lspr-haiss-size', suite: 23, title: 'LSPR: gold-sphere diameter from the A_spr/A_450 correlation', mode: 'lspr', inputs: { mode: 'lspr', lspr_table: LSPR, lspr_metal: 'Au' }, kpi: 'd_est', expected: Math.exp(3.00 * lsprA(520) / lsprA(450) - 2.20), tol: 0.02, kind: 'analytical', reference: 'Haiss, Thanh, Aveyard & Fernig, Analytical Chemistry: d = exp(B₁·A_spr/A_450 − B₂), B₁ = 3.00, B₂ = 2.20', note: 'Checks the implementation of the published correlation; the correlation itself scatters by 10–20 % against microscopy.' },
  // ---- Raman ----
  { id: 's23-raman-id-ig', suite: 23, title: 'Raman: I_D/I_G of Lorentzian bands with height ratio 0.5', mode: 'raman', inputs: { mode: 'raman', raman_table: RAMAN }, kpi: 'idig', expected: 0.5, tol: 0.03, kind: 'analytical', reference: 'Ratio of fitted band heights of the D (1350 cm⁻¹) and G (1582 cm⁻¹) bands', note: 'Baseline clipping under the broad D band costs a few percent of its height.' },
  { id: 's23-raman-la', suite: 23, title: 'Raman: graphitic in-plane crystallite size at 532 nm', mode: 'raman', inputs: { mode: 'raman', raman_table: RAMAN, raman_laser: 532 }, kpi: 'La', expected: 2.4e-10 * 532 ** 4 / 0.5, tol: 0.03, kind: 'analytical', reference: 'Cançado et al., Applied Physics Letters: L_a (nm) = 2.4×10⁻¹⁰·λ⁴·(I_D/I_G)⁻¹', note: 'Same baseline effect as the I_D/I_G benchmark.' },
  { id: 's23-raman-si-stress', suite: 23, title: 'Raman: stress from a −2 cm⁻¹ shift of the silicon phonon', mode: 'raman', inputs: { mode: 'raman', raman_table: RAMAN, raman_stress: 500 }, kpi: 'si_stress', expected: 1000, tol: 0.02, kind: 'analytical', reference: 'Uniaxial stress in (100) silicon ≈ −500 MPa per cm⁻¹ shift of the 520.7 cm⁻¹ phonon (De Wolf, Semiconductor Science and Technology)', note: 'Band placed at 518.7 cm⁻¹; tensile stress is positive.' },
  // ---- indentation ----
  { id: 's23-op-fused-silica-E', suite: 23, title: 'Oliver–Pharr: modulus of fused silica from a synthetic Berkovich curve', mode: 'op', inputs: { mode: 'op', op_table: OP.table, op_nu: 0.17, op_beta: 1.034, op_eps: 0.75, op_C1: 0, op_film: 0 }, kpi: 'E', expected: 72, tol: 0.01, kind: 'reference-data', reference: 'Fused-silica calibration standard, E ≈ 72 GPa (Oliver & Pharr, Journal of Materials Research)', note: 'Curve constructed from E = 72 GPa, H = 9 GPa with a diamond tip (1141 GPa, ν = 0.07) and unloading exponent 1.3.' },
  { id: 's23-op-fused-silica-H', suite: 23, title: 'Oliver–Pharr: hardness of fused silica from the same curve', mode: 'op', inputs: { mode: 'op', op_table: OP.table, op_nu: 0.17, op_beta: 1.034, op_eps: 0.75, op_C1: 0, op_film: 0 }, kpi: 'H', expected: 9, tol: 0.01, kind: 'reference-data', reference: 'Fused-silica calibration standard, H ≈ 9 GPa (Oliver & Pharr, Journal of Materials Research)' },
  { id: 's23-op-stiffness', suite: 23, title: 'Oliver–Pharr: unloading contact stiffness S = (2β/√π)·E_r·√A', mode: 'op', inputs: { mode: 'op', op_table: OP.table }, kpi: 'S', expected: OP.S, tol: 0.01, kind: 'analytical', reference: 'Sneddon stiffness relation with the Berkovich factor β = 1.034 (Oliver & Pharr, Journal of Materials Research)' },
  { id: 's23-hertz-modulus', suite: 23, title: 'Hertz: modulus of a 5 MPa half-space under a 1 µm rigid sphere', mode: 'hertz', inputs: { mode: 'hertz', hz_table: HERTZ, hz_R: 1000, hz_nu: 0.4 }, kpi: 'E', expected: 5e6, tol: 0.005, kind: 'analytical', reference: 'Hertz contact F = (4/3)·E/(1 − ν²)·√R·δ^(3/2); Johnson, Contact Mechanics' },
  { id: 's23-hertz-contact-radius', suite: 23, title: 'Hertz: contact radius a = √(R·δ) at 80 nm indentation', mode: 'hertz', inputs: { mode: 'hertz', hz_table: HERTZ, hz_R: 1000, hz_nu: 0.4 }, kpi: 'a', expected: Math.sqrt(1000 * 80), tol: 0.005, kind: 'analytical', reference: 'Hertz contact geometry; Johnson, Contact Mechanics' },
  { id: 's23-hertz-pressure', suite: 23, title: 'Hertz: peak contact pressure p₀ = (2/π)·E_r·√(δ/R)', mode: 'hertz', inputs: { mode: 'hertz', hz_table: HERTZ, hz_R: 1000, hz_nu: 0.4 }, kpi: 'p0', expected: 2 / PI * 5e6 / 0.84 * Math.sqrt(0.08), tol: 0.01, kind: 'analytical', reference: 'Hertz contact p₀ = 3F/(2πa²) = 2E_r·a/(πR); Johnson, Contact Mechanics' },
  // ---- electrical ----
  { id: 's23-fpp-sheet-factor', suite: 23, title: 'Four-point probe: infinite thin sheet, R_s = (π/ln 2)·V/I', mode: 'fpp', inputs: { mode: 'fpp', fpp_table: FPP, fpp_s: 0.1, fpp_d: 300, el_thick: 100 }, kpi: 'Rs', expected: 2 * 4.532, tol: 0.001, kind: 'analytical', reference: 'Collinear four-point-probe factor π/ln 2 = 4.532 (Smits; Schroder, Semiconductor Material and Device Characterization)', note: 'V/I = 2 Ω; sample 3000 probe spacings across so the lateral correction is negligible.' },
  { id: 's23-fpp-resistivity', suite: 23, title: 'Four-point probe: resistivity of a 100 nm film, ρ = R_s·t', mode: 'fpp', inputs: { mode: 'fpp', fpp_table: FPP, fpp_s: 0.1, fpp_d: 300, el_thick: 100 }, kpi: 'rho', expected: PI / LN2 * 2 * 100e-9, tol: 0.001, kind: 'analytical', reference: 'Thin-sheet four-point-probe relation ρ = (π/ln 2)·(V/I)·t (Schroder, Semiconductor Material and Device Characterization)' },
  { id: 's23-fpp-finite-disc', suite: 23, title: 'Four-point probe: correction for a circular sample 10 probe spacings across', mode: 'fpp', inputs: { mode: 'fpp', fpp_table: FPP, fpp_s: 1, fpp_d: 10, el_thick: 100 }, kpi: 'Rs', expected: 2 * PI / (LN2 + Math.log(103 / 97)), tol: 0.001, kind: 'analytical', reference: 'Smits correction for probes centred on a circular sample: R_s = (V/I)·π/{ln 2 + ln[((d/s)² + 3)/((d/s)² − 3)]} (≈ 4.17·V/I at d/s = 10)' },
  { id: 's23-fpp-bulk-limit', suite: 23, title: 'Four-point probe: thick-sample limit ρ = 2π·s·V/I', mode: 'fpp', inputs: { mode: 'fpp', fpp_table: FPP, fpp_s: 0.1, fpp_d: 300, el_thick: 1e6 }, kpi: 'rho', expected: 2 * PI * 1e-4 * 2, tol: 0.001, kind: 'analytical', reference: 'Semi-infinite four-point-probe relation ρ = 2π·s·V/I (Valdes; Schroder, Semiconductor Material and Device Characterization)', note: 'Thickness 10 probe spacings: the sinh thickness correction differs from the bulk limit by about 10⁻⁵.' },
  { id: 's23-hall-carrier-density', suite: 23, title: 'Hall effect: carrier density n = 1/(q·|R_H|)', mode: 'hall', inputs: { mode: 'hall', hall_table: HALL, el_thick: 1000, hall_I: 1 }, kpi: 'n', expected: HALL_N, tol: 0.001, kind: 'analytical', reference: 'Single-carrier Hall effect V_H = I·B/(q·n·t); Kittel, Introduction to Solid State Physics', note: 'n-type film, 10¹⁹ cm⁻³, 1 µm thick, 1 mA, with a constant misalignment offset that the slope must reject.' },
  { id: 's23-hall-coefficient', suite: 23, title: 'Hall effect: Hall coefficient R_H = −1/(q·n)', mode: 'hall', inputs: { mode: 'hall', hall_table: HALL, el_thick: 1000, hall_I: 1 }, kpi: 'RH', expected: -1 / (Q_E * HALL_N), tol: 0.001, kind: 'analytical', reference: 'Single-carrier Hall coefficient; Kittel, Introduction to Solid State Physics' },
  { id: 's23-vdp-symmetric', suite: 23, title: 'van der Pauw: symmetric sample, R_s = π·R/ln 2', mode: 'hall', inputs: { mode: 'hall', hall_table: HALL, hall_RA: 10, hall_RB: 10 }, kpi: 'Rs', expected: PI * 10 / LN2, tol: 1e-5, kind: 'analytical', reference: 'van der Pauw equation exp(−πR_A/R_s) + exp(−πR_B/R_s) = 1 with R_A = R_B' },
  { id: 's23-vdp-golden-ratio', suite: 23, title: 'van der Pauw: R_A = 2·R_B has the closed-form solution R_s = π·R_B/ln φ', mode: 'hall', inputs: { mode: 'hall', hall_table: HALL, hall_RA: 20, hall_RB: 10 }, kpi: 'Rs', expected: PI * 10 / Math.log((1 + Math.sqrt(5)) / 2), tol: 1e-5, kind: 'analytical', reference: 'van der Pauw equation: with x = exp(−πR_B/R_s), x² + x = 1 gives x = 1/φ (golden ratio)' },
  { id: 's23-hall-mobility', suite: 23, title: 'Hall mobility µ = |R_H|/(R_s·t)', mode: 'hall', inputs: { mode: 'hall', hall_table: HALL, el_thick: 1000, hall_I: 1, hall_RA: 10, hall_RB: 10 }, kpi: 'mu', expected: 1 / (Q_E * HALL_N) / (PI * 10 / LN2 * 1e-6), tol: 0.001, kind: 'analytical', reference: 'Drude relation σ = q·n·µ with the van der Pauw resistivity' },
  { id: 's23-eis-rct', suite: 23, title: 'Impedance: charge-transfer resistance recovered from an exact Randles spectrum', mode: 'eis', inputs: { mode: 'eis', eis_table: EIS }, kpi: 'Rct', expected: 100, tol: 0.01, kind: 'analytical', reference: 'Randles equivalent circuit Z = R_s + [jωC + 1/(R_ct + σ(1 − j)/√ω)]⁻¹; Bard & Faulkner, Electrochemical Methods', note: 'Generating values R_s = 10 Ω, R_ct = 100 Ω, C_dl = 20 µF, σ = 15 Ω·s^-½, 60 frequencies from 100 kHz to 10 mHz.' },
  { id: 's23-eis-rs', suite: 23, title: 'Impedance: series resistance of the same spectrum', mode: 'eis', inputs: { mode: 'eis', eis_table: EIS }, kpi: 'Rs', expected: 10, tol: 0.01, kind: 'analytical', reference: 'Randles equivalent circuit; Bard & Faulkner, Electrochemical Methods' },
  { id: 's23-eis-cdl', suite: 23, title: 'Impedance: double-layer capacitance of the same spectrum', mode: 'eis', inputs: { mode: 'eis', eis_table: EIS }, kpi: 'Cdl', expected: 2e-5, tol: 0.02, kind: 'analytical', reference: 'Randles equivalent circuit; for α = 1 the Brug effective capacitance equals the capacitor value', note: 'The fitted CPE exponent must reach its upper bound of 1.' },
  { id: 's23-eis-warburg', suite: 23, title: 'Impedance: Warburg coefficient of the same spectrum', mode: 'eis', inputs: { mode: 'eis', eis_table: EIS }, kpi: 'sigma_w', expected: 15, tol: 0.01, kind: 'analytical', reference: 'Semi-infinite Warburg impedance Z_W = σ(1 − j)/√ω; Bard & Faulkner, Electrochemical Methods' },
  // ---- thermal ----
  { id: 's23-flash-parker', suite: 23, title: 'Laser flash: diffusivity returned for an exact adiabatic Parker curve', mode: 'flash', inputs: { mode: 'flash', flash_table: FLASH, fl_L: 2 }, kpi: 'alpha', expected: 1e-5, tol: 0.01, kind: 'analytical', reference: 'Parker, Jenkins, Butler & Abbott, Journal of Applied Physics: α = 0.1388·L²/t½', note: 'Rear-face series solution for α = 1.0×10⁻⁵ m²/s through 2 mm sampled every 1 ms; five-point smoothing of the record.' },
  { id: 's23-flash-half-time', suite: 23, title: 'Laser flash: half-rise time t½ = 0.1388·L²/α', mode: 'flash', inputs: { mode: 'flash', flash_table: FLASH, fl_L: 2 }, kpi: 't_half', expected: 0.1388 * 4e-6 / 1e-5 * 1e3, tol: 0.01, kind: 'analytical', reference: 'Parker, Jenkins, Butler & Abbott, Journal of Applied Physics' },
  { id: 's23-flash-conductivity', suite: 23, title: 'Laser flash: conductivity k = α·ρ·c_p', mode: 'flash', inputs: { mode: 'flash', flash_table: FLASH, fl_L: 2, rho: 2700, fl_cp: 900 }, kpi: 'k', expected: 1e-5 * 2700 * 900, tol: 0.01, kind: 'analytical', reference: 'Definition of thermal diffusivity α = k/(ρ·c_p)' },
  // ---- mechanical ----
  { id: 's23-tensile-modulus', suite: 23, title: "Tensile: Young's modulus of a Ramberg–Osgood curve with E = 200 GPa", mode: 'tensile', inputs: { mode: 'tensile', tens_table: TENS_STEEL }, kpi: 'E', expected: 200000, tol: 0.005, kind: 'analytical', reference: 'Ramberg–Osgood relation ε = σ/E + 0.002·(σ/σ_y)^n', note: 'n = 15: the plastic term is below 10⁻⁵ of the elastic strain over the fitted segment.' },
  { id: 's23-tensile-proof-stress', suite: 23, title: 'Tensile: 0.2 % proof stress of the same curve', mode: 'tensile', inputs: { mode: 'tensile', tens_table: TENS_STEEL, tens_offset: 0.2 }, kpi: 'sy', expected: 400, tol: 0.005, kind: 'analytical', reference: 'Ramberg–Osgood relation: the plastic strain equals 0.2 % exactly at σ = σ_y' },
  { id: 's23-tensile-auxetic-nu', suite: 23, title: "Tensile: negative Poisson's ratio of an auxetic coupon (ν = −0.6)", mode: 'tensile', inputs: { mode: 'tensile', tens_table: TENS_AUX }, kpi: 'nu', expected: -0.6, tol: 0.002, kind: 'analytical', reference: "Definition ν = −ε_transverse/ε_axial; Lakes, Science (foams with negative Poisson's ratio)" },
  { id: 's23-tensile-resilience', suite: 23, title: 'Tensile: modulus of resilience σ_y²/(2E)', mode: 'tensile', inputs: { mode: 'tensile', tens_table: TENS_STEEL }, kpi: 'resil', expected: 400 * 400 / (2 * 200000), tol: 0.01, kind: 'analytical', reference: 'Elastic strain energy density at yield (Dieter, Mechanical Metallurgy)' },
  // ---- magnetic ----
  { id: 's23-vsm-coercivity', suite: 23, title: 'VSM: coercive field of a loop built with H_c = 20 kA/m', mode: 'vsm', inputs: { mode: 'vsm', vsm_table: VSM }, kpi: 'Hc', expected: VSM_HC, tol: 0.002, kind: 'analytical', reference: 'Zero crossing of M(H) = M_s·tanh((H ± H_c)/H₀)' },
  { id: 's23-vsm-remanence', suite: 23, title: 'VSM: remanence M_r = M_s·tanh(H_c/H₀)', mode: 'vsm', inputs: { mode: 'vsm', vsm_table: VSM }, kpi: 'Mr', expected: VSM_MS * Math.tanh(VSM_HC / VSM_H0), tol: 0.002, kind: 'analytical', reference: 'Zero-field value of M(H) = M_s·tanh((H ± H_c)/H₀)' },
  { id: 's23-vsm-saturation', suite: 23, title: 'VSM: saturation magnetisation of a fully saturated loop', mode: 'vsm', inputs: { mode: 'vsm', vsm_table: VSM }, kpi: 'Ms', expected: VSM_MS, tol: 0.002, kind: 'analytical', reference: 'High-field limit of M(H) = M_s·tanh((H ± H_c)/H₀)' },
  { id: 's23-vsm-loss', suite: 23, title: 'VSM: hysteresis loss per cycle W = 4·µ₀·H_c·M_s', mode: 'vsm', inputs: { mode: 'vsm', vsm_table: VSM }, kpi: 'loss', expected: 4 * MU0 * VSM_HC * 1e3 * VSM_MS * 1e3, tol: 0.005, kind: 'analytical', reference: 'Loop area µ₀·∮H dM of two identical saturating branches displaced by ±H_c (Cullity & Graham, Introduction to Magnetic Materials)', note: 'Exact when the maximum field is far beyond saturation (here 1000 kA/m against H₀ = 30 kA/m).' },
  { id: 's23-vsm-specific-magnetisation', suite: 23, title: 'VSM: specific magnetisation σ_s = M_s/ρ', mode: 'vsm', inputs: { mode: 'vsm', vsm_table: VSM, rho: 5000 }, kpi: 'sigma_s', expected: 80, tol: 0.002, kind: 'analytical', reference: 'Definition of mass magnetisation: 400 kA/m ÷ 5000 kg/m³ = 80 A·m²/kg' },
  // ---- metrology ----
  { id: 's23-gum-product-quotient', suite: 23, title: 'GUM: relative uncertainties of a product/quotient add in quadrature (1 %, 2 %, 2 % → 3 %)', mode: 'gum', inputs: { mode: 'gum', gum_red: 'fpp', gum_table: GUM_FPP }, kpi: 'uc', expected: 0.03 * PI / LN2 * 100e-9, tol: 0.001, kind: 'analytical', reference: 'JCGM 100 (GUM), law of propagation of uncertainty for uncorrelated inputs applied to ρ = (π/ln 2)·(V/I)·t' },
  { id: 's23-gum-expanded', suite: 23, title: 'GUM: expanded uncertainty U = k·u_c with k = 2.5', mode: 'gum', inputs: { mode: 'gum', gum_red: 'fpp', gum_table: GUM_FPP, gum_k: 2.5 }, kpi: 'U', expected: 2.5 * 0.03 * PI / LN2 * 100e-9, tol: 0.001, kind: 'analytical', reference: 'JCGM 100 (GUM), expanded uncertainty' },
  { id: 's23-gum-stokes-einstein', suite: 23, title: 'GUM model value: Stokes–Einstein diameter for D = 4.91×10⁻¹² m²/s in water at 25 °C', mode: 'gum', inputs: { mode: 'gum', gum_red: 'stokes', gum_table: GUM_STOKES }, kpi: 'y', expected: 100, tol: 1e-6, kind: 'analytical', reference: 'Stokes–Einstein relation d = k_B·T/(3π·η·D)' },
  { id: 's23-gum-variance-share', suite: 23, title: 'GUM: variance share of each of two equal dominant inputs (4/9)', mode: 'gum', inputs: { mode: 'gum', gum_red: 'fpp', gum_table: GUM_FPP }, kpi: 'dom_share', expected: 4 / 9, tol: 0.001, kind: 'analytical', reference: 'JCGM 100 (GUM), uncertainty budget: share = (cᵢuᵢ)²/u_c²' },
  { id: 's23-repeat-type-a', suite: 23, title: 'Repeatability: standard uncertainty of the mean of {9.8, 10.0, 10.2, 10.1, 9.9}', mode: 'repeat', inputs: { mode: 'repeat', rep_table: REP }, kpi: 'se', expected: Math.sqrt(0.1 / 4 / 5), tol: 1e-6, kind: 'analytical', reference: 'JCGM 100 (GUM), type-A evaluation u = s/√n' },
  { id: 's23-repeat-en', suite: 23, title: 'Calibration check: E_n number', mode: 'repeat', inputs: { mode: 'repeat', rep_table: REP, rep_ref: 9.9, rep_Uref: 0.1, rep_uB: 0 }, kpi: 'En', expected: 0.1 / Math.sqrt(4 * 0.005 + 0.01), tol: 1e-6, kind: 'analytical', reference: 'E_n = (x − X_ref)/√(U_lab² + U_ref²) with k = 2 expanded uncertainties (ISO 13528 proficiency-testing statistics)' },
  { id: 's23-repeat-expanded', suite: 23, title: 'Repeatability: expanded uncertainty U = 2·√(s²/n + u_B²)', mode: 'repeat', inputs: { mode: 'repeat', rep_table: REP, rep_uB: 0.05 }, kpi: 'U', expected: 2 * Math.sqrt(0.005 + 0.0025), tol: 1e-6, kind: 'analytical', reference: 'JCGM 100 (GUM), combination of type-A and type-B standard uncertainties, coverage factor 2' },
  { id: 's23-repeat-grubbs', suite: 23, title: 'Grubbs statistic G = max|x − x̄|/s of the same set', mode: 'repeat', inputs: { mode: 'repeat', rep_table: REP }, kpi: 'G', expected: 0.2 / Math.sqrt(0.025), tol: 1e-6, kind: 'analytical', reference: 'Grubbs test for a single outlier' },
];
