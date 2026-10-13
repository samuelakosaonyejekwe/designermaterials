// Reference benchmarks for Suite 24 — Verification, Calibration, Validation and Uncertainty.
// Every expected value is a literal or is computed here from a closed-form expression; nothing is taken from the suite.

// ---------- Richardson extrapolation / grid convergence index ----------
// Three-grid example of Roache: normalised solutions on fine / medium / coarse grids, refinement ratio 2, Fs = 1.25.
const RF = [0.97050, 0.96854, 0.96178], RR = 2, RFS = 1.25;
const R_P = Math.log((RF[2] - RF[1]) / (RF[1] - RF[0])) / Math.log(RR);                 // 1.786
const R_FEXT = RF[0] + (RF[0] - RF[1]) / (Math.pow(RR, R_P) - 1);                       // 0.9713
const R_GCI12 = 100 * RFS * Math.abs((RF[0] - RF[1]) / RF[0]) / (Math.pow(RR, R_P) - 1);   // 0.103 %
const R_GCI23 = 100 * RFS * Math.abs((RF[1] - RF[2]) / RF[1]) / (Math.pow(RR, R_P) - 1);   // 0.356 %
const ROACHE = { mode: 'gci', gh1: 1, gh2: 2, gh3: 4, gf1: RF[0], gf2: RF[1], gf3: RF[2], gFs: RFS, gOrder: 2 };
const ROACHE_REF = 'Roache, Verification and Validation in Computational Science and Engineering — grid convergence index; three-grid example as used in the NASA Glenn (NPARC Alliance) tutorial on spatial grid convergence';
// exact power-law error models f(h) = f0 + C·h^p sampled on three grids
const pw = (f0, C, p, h) => f0 + C * Math.pow(h, p);
const NONUNI = { mode: 'gci', gh1: 1, gh2: 1.5, gh3: 3, gf1: pw(2, 0.3, 2, 1), gf2: pw(2, 0.3, 2, 1.5), gf3: pw(2, 0.3, 2, 3), gOrder: 2 };
const THIRD = { mode: 'gci', gh1: 0.1, gh2: 0.15, gh3: 0.225, gf1: pw(5, -40, 3, 0.1), gf2: pw(5, -40, 3, 0.15), gf3: pw(5, -40, 3, 0.225), gOrder: 3, gFs: 1.25 };
const CELIK_REF = 'Richardson extrapolation with the fixed-point observed-order procedure for unequal refinement ratios (Celik et al., ASME Journal of Fluids Engineering procedure for estimation of discretisation uncertainty)';

// ---------- heat-equation verification problem: von Neumann closed forms of the discrete scheme ----------
// u_t = α·u_xx on [0, L], u = sin(πx/L) initially. Second-order central differences on nx intervals damp this mode at the
// rate λ_h = (4α/h²)·sin²(πh/2L); a θ-scheme multiplies its amplitude by a fixed factor every step.
const FO = 0.1, A_EXACT = Math.exp(-FO * Math.PI ** 2);
const zT = nx => FO * 4 * nx * nx * Math.sin(Math.PI / (2 * nx)) ** 2;                 // λ_h·t_end
const aCN = (nx, nt) => { const z = zT(nx) / nt; return Math.pow((1 - z / 2) / (1 + z / 2), nt); };
const aBE = (nx, nt) => Math.pow(1 + zT(nx) / nt, -nt);
const ERR_FINE = Math.abs(aCN(40, 4000) - A_EXACT) / A_EXACT * Math.sqrt(40 / (2 * 39));   // RMS over the 39 interior nodes
const QOI_ERR = 100 * Math.abs(aCN(40, 4000) ** 2 - A_EXACT ** 2) / A_EXACT ** 2;          // ∫u² dx = A²·L/2 on the grid
const ERR_DT_CN = Math.abs(aCN(20, 64) - aCN(20, 1024)) / A_EXACT;
const ERR_DT_BE = Math.abs(aBE(20, 64) - aCN(20, 1024)) / A_EXACT;
const MMS = { mode: 'mms', vcase: 'decay', scheme: 'cn', alpha: 1e-6, Lx: 0.01, Fo: FO, nx0: 10, nt0: 16, odeN: 8 };
const VN_REF = 'Von Neumann (Fourier) analysis of the θ-scheme with second-order central differences for the heat equation — e.g. Morton & Mayers, Numerical Solution of Partial Differential Equations';
const ORDER_NOTE = 'The formal order is reached only in the limit of vanishing step size; the absolute tolerance covers the pre-asymptotic deviation of the three-level study.';

// ---------- calibration: noiseless data generated from known parameters ----------
const RG = 8.314462618;
const table = (xs, f) => 'x,y\n' + xs.map(x => `${x},${f(x)}`).join('\n');
const FAMS = [
  ['arrhenius', 'Arrhenius law', [300, 320, 340, 360, 380, 400], [20, 50], ['ln A', 'activation energy (kJ/mol)'], (x, p) => Math.exp(p[0] - p[1] * 1e3 / (RG * x))],
  ['power', 'power law', [1, 2, 4, 8, 16], [2.5, 1.7], ['prefactor a', 'exponent n'], (x, p) => p[0] * Math.pow(x, p[1])],
  ['exp', 'exponential decay', [0, 1, 2, 3, 4, 5], [3, -0.4], ['amplitude a', 'rate b'], (x, p) => p[0] * Math.exp(p[1] * x)],
  ['linear', 'straight line', [0, 1, 2, 3, 4], [1.5, -0.75], ['intercept a', 'slope b'], (x, p) => p[0] + p[1] * x],
  ['langmuir', 'Langmuir isotherm', [2, 5, 10, 20, 50, 100, 200], [120, 0.08], ['capacity q_max', 'affinity K'], (x, p) => p[0] * p[1] * x / (1 + p[1] * x)],
  ['jmak', 'JMAK transformation', [10, 20, 30, 40, 50, 60, 80, 100], [0.02, 2.5], ['rate constant k', 'Avrami exponent n'], (x, p) => 1 - Math.exp(-Math.pow(p[0] * x, p[1]))],
  ['parabolic', 'parabolic growth law', [1, 4, 9, 16, 25], [4], ['rate constant k'], (x, p) => Math.sqrt(p[0] * x)],
  ['difflim', 'diffusion-limited growth', [0, 5, 10, 20, 30, 50], [5, 3], ['initial radius r₀', 'rate constant k'], (x, p) => Math.sqrt(p[0] * p[0] + p[1] * x)],
  ['firstorder', 'first-order kinetics', [1, 2, 4, 6, 10, 15], [10, 0.3], ['plateau y∞', 'rate constant k'], (x, p) => p[0] * (1 - Math.exp(-p[1] * x))],
];
const recovery = FAMS.flatMap(([key, name, xs, p, pn, f]) => p.map((v, j) => ({
  id: `s24-calib-recover-${key}-p${j}`, suite: 24, title: `Noiseless ${name}: recovery of the ${pn[j]}`, mode: 'calib',
  inputs: { mode: 'calib', cmodel: key, cdata: table(xs, x => f(x, p)), bayes: false, rxn: '' },
  kpi: 'p' + j, expected: v, tol: 1e-5, kind: 'analytical',
  reference: 'Identity test of nonlinear least squares: data generated from the model without noise have the generating parameters as the unique minimiser (χ² = 0) — e.g. Bates & Watts, Nonlinear Regression Analysis and Its Applications',
  note: 'The tolerance covers termination of the Levenberg–Marquardt iteration and finite-difference Jacobians only.',
})));

// ordinary least squares of a straight line: textbook closed forms
const LX = [1, 2, 3, 4, 5, 6], LY = [3.1, 4.9, 7.2, 8.8, 11.1, 13.2], LN = LX.length;
const sum = a => a.reduce((s, v) => s + v, 0), LSX = sum(LX), LSXX = sum(LX.map(v => v * v)), LMX = LSX / LN, LMY = sum(LY) / LN;
const LSxx = sum(LX.map(v => (v - LMX) ** 2)), L_B = sum(LX.map((v, i) => (v - LMX) * (LY[i] - LMY))) / LSxx, L_A = LMY - L_B * LMX;
const L_SSE = sum(LX.map((v, i) => (LY[i] - L_A - L_B * v) ** 2)), L_S2 = L_SSE / (LN - 2);
const L_SEB = Math.sqrt(L_S2 / LSxx), L_SEA = Math.sqrt(L_S2 * LSXX / (LN * LSxx)), L_RHO = LSX / Math.sqrt(LN * LSXX);
const OLS = { mode: 'calib', cmodel: 'linear', cdata: 'x,y\n' + LX.map((v, i) => `${v},${LY[i]}`).join('\n'), bayes: false, rxn: '' };
const OLS_REF = 'Simple linear regression: least-squares estimators, residual variance and parameter covariance — e.g. Draper & Smith, Applied Regression Analysis';

// ---------- validation metrics on a hand-computed table ----------
const VX = [1, 2, 3, 4, 5], VM = [10, 20, 30, 40, 50], VP = [11, 19.5, 32, 38, 51], UM = 0.6, UP = 0.8, VK = 1.5;
const VE = VP.map((v, i) => v - VM[i]), VUC = Math.hypot(UM, UP), VMEAN = sum(VM) / VM.length;
const VALID = { mode: 'valid', predSrc: 'table', kCov: VK, xApp: 7, vdata: 'x,measured,u_measured,predicted,u_predicted\n' + VX.map((x, i) => [x, VM[i], UM, VP[i], UP].join(',')).join('\n') };
const VAL_REF = 'Definitions of validation metrics (RMSE, bias, coefficient of determination, normalised error Eₙ with combined standard uncertainty) — ASME V&V 20 validation-uncertainty approach and the GUM for combining uncertainties';

// ---------- uncertainty propagation ----------
// NPV model with volume and discount rate fixed: NPV = −CAPEX + (price − unit cost)·Q·a, a = annuity factor ⇒ linear and additive
const ANN = (1 - Math.pow(1.1, -10)) / 0.1, NQ = 20000, NPV0 = -2e6 + (40 - 22) * NQ * ANN;
const NV = [(0.15 * 2e6) ** 2, (NQ * ANN * 0.10 * 40) ** 2, (NQ * ANN * 0.10 * 22) ** 2], NSD = Math.sqrt(sum(NV)), NSAMP = 20000;
const NPV = { mode: 'uq', umodel: 'npv', uqOverride: 'volume=20000,0; r=0.1,0', dist: 'normal', sampler: 'lhs', nSamp: NSAMP, rho12: 0, cvScale: 1, limRatio: 1.2, limSide: 'above' };
const LIN_REF = 'Exact moments and variance decomposition of a linear additive function of independent inputs: Var(Y) = Σcᵢ²σᵢ², Sᵢ = S_Tᵢ = cᵢ²σᵢ²/Var(Y) — Saltelli et al., Global Sensitivity Analysis: The Primer';
// standard normal distribution function (series for erf, adequate for |x| < 3)
const phi = x => { const t = x / Math.SQRT2; let s = t, term = t; for (let n = 1; n < 60; n++) { term *= -t * t / n; s += term / (2 * n + 1); } return 0.5 * (1 + 2 / Math.sqrt(Math.PI) * s); };
const NPF = 1 - phi((1.2 * NPV0 - NPV0) / NSD);
// product of independent lognormal inputs (Arrhenius growth length with Eₐ and T fixed): L = A·k(T)·t
const L0 = 4e7 * Math.exp(-60e3 / (RG * 363)) * 10800, LCV2 = (1 + 0.10 ** 2) * (1 + 0.02 ** 2) - 1;
const LOGN = { mode: 'uq', umodel: 'arrh', uqOverride: 'Ea=60,0; T=363,0', dist: 'lognormal', sampler: 'lhs', nSamp: NSAMP, rho12: 0, cvScale: 1 };
const LOGN_REF = 'Product of independent lognormal variables is lognormal: E[Y] = ΠE[Xᵢ], 1 + CV_Y² = Π(1 + CVᵢ²), median = E[Y]/√(1 + CV_Y²) — e.g. Ang & Tang, Probability Concepts in Engineering';
// DLVO barrier of two equal spheres at the nominal inputs, maximised on a fine separation grid
const DLVO0 = (() => {
  const kB = 1.380649e-23, e = 1.602176634e-19, NA = 6.02214076e23, eps = 78.5 * 8.8541878128e-12, T = 298.15, AH = 1e-19, psi = 0.050, Rp = 2.5e-8, I = 1e-3;
  const kap = Math.sqrt(2 * NA * e * e * I * 1e3 / (eps * kB * T)); let mx = -Infinity;
  for (let i = 0; i <= 200000; i++) { const h = 2e-10 * Math.pow(250, i / 200000), V = -AH * Rp / (12 * h) + 2 * Math.PI * eps * Rp * psi * psi * Math.log(1 + Math.exp(-kap * h)); if (V > mx) mx = V; }
  return mx / (kB * T);
})();

// ---------- published measured data sets with certified least-squares results (NIST Statistical Reference Datasets, StRD) ----------
// Observed data transcribed from the NIST/ITL StRD data files as [x, y] pairs in file order (the files list y first);
// the certified values quoted in the entries are those printed in the same files. https://www.itl.nist.gov/div898/strd/
const strd = rows => 'x,y\n' + rows.map(r => r.join(',')).join('\n');
// Norris.dat — J. Norris (NIST), calibration of ozone monitors; 36 observations, model y = B0 + B1·x
const STRD_NORRIS = strd([[0.2, 0.1], [337.4, 338.8], [118.2, 118.1], [884.6, 888.0], [10.1, 9.2], [226.5, 228.1], [666.3, 668.5], [996.3, 998.5], [448.6, 449.1], [777.0, 778.9], [558.2, 559.2], [0.4, 0.3],
  [0.6, 0.1], [775.5, 778.1], [666.9, 668.8], [338.0, 339.3], [447.5, 448.9], [11.6, 10.8], [556.0, 557.7], [228.1, 228.3], [995.8, 998.0], [887.6, 888.8], [120.2, 119.6], [0.3, 0.3],
  [0.3, 0.6], [556.8, 557.6], [339.1, 339.3], [887.2, 888.0], [999.0, 998.5], [779.0, 778.9], [11.1, 10.2], [118.3, 117.6], [229.2, 228.9], [669.1, 668.4], [448.9, 449.2], [0.5, 0.2]]);
// DanWood.dat — energy radiated by a carbon-filament lamp (y, per cm² per second) against absolute filament temperature (x, in 1000 K); 6 observations, model y = b1·x^b2
const STRD_DANWOOD = strd([[1.309, 2.138], [1.471, 3.421], [1.490, 3.597], [1.565, 4.340], [1.611, 4.882], [1.680, 5.660]]);
// Misra1a.dat — D. Misra (NIST, 1978), dental research in monomolecular adsorption: volume (y) against pressure (x); 14 observations, model y = b1·(1 − exp(−b2·x))
const STRD_MISRA1A = strd([[77.6, 10.07], [114.9, 14.73], [141.1, 17.94], [190.8, 23.93], [239.9, 29.61], [289.0, 35.18], [332.8, 40.02], [378.4, 44.82], [434.8, 50.76], [477.3, 55.05], [536.8, 61.01], [593.1, 66.40], [689.1, 75.47], [760.0, 81.78]]);
// BoxBOD.dat — biochemical oxygen demand (y, mg/L) against incubation time (x, days); 6 observations, model y = b1·(1 − exp(−b2·x))
const STRD_BOXBOD = strd([[1, 109], [2, 149], [3, 149], [5, 191], [7, 213], [10, 224]]);
const strdSrc = (name, title) => ({ authors: 'National Institute of Standards and Technology, Information Technology Laboratory', title: `Statistical Reference Datasets (StRD): ${title}`, journal: 'NIST Standard Reference Database 140 (data content last updated 2003)', year: 2003, doi: '10.18434/T43G6C', url: `https://www.itl.nist.gov/div898/strd/${name}` });
const STRD_NOTE = 'What is tested is the data reduction (unweighted least squares on the published observations), not a physical model; no sigma column is given, as in the source. The Metropolis sampling is switched off because it does not enter the compared estimate. NIST certifies the values to 11 significant digits; the tolerance of 10⁻⁶ allows for the termination criterion of the Levenberg–Marquardt iteration.';

// ---------- porous zirconia: measured tensile modulus against porosity (Pabst, Gregorová, Tichá & Týnová, Ceramics − Silikáty 48 (4) 165–174, 2004) ----------
// Table 2 of the paper, transcribed exactly as [porosity φ, tensile modulus E in GPa]; one specimen per pair. Z-P: potato starch, Z-C: corn starch.
const PABST_ZP = strd([[0.132, 136], [0.191, 113], [0.249, 96], [0.315, 69], [0.424, 38]]);
const PABST_ZC = strd([[0.073, 167], [0.158, 136], [0.200, 124], [0.298, 85], [0.415, 52]]);
const PABST_SRC = { authors: 'Pabst, W.; Gregorová, E.; Tichá, G.; Týnová, E.', title: 'Effective elastic properties of alumina-zirconia composite ceramics — Part 4. Tensile modulus of porous alumina and zirconia', journal: 'Ceramics − Silikáty 48 (4), 165–174', year: 2004, doi: null, url: 'https://www.irsm.cas.cz/materialy/cs_content/2004/Pabst_CS_2004_0003.pdf' };
const PABST_COND = 'Tetragonal zirconia (3 mol.% Y₂O₃, TZ-3YE) shaped by starch consolidation casting and fired at 1400 °C; cylindrical rods about 4 mm in diameter and 50–80 mm long; total porosity from the Archimedes bulk density (theoretical density 6.1 g/cm³), tensile modulus by the resonant frequency technique';
const PABST_NOTE = 'The x column of the calibration table is the porosity φ and the y column the modulus in GPa; the exponential family y = a·exp(b·x) is the Spriggs relation E = E₀·exp(−Bφ) with a = E₀ and b = −B (the paper\'s Model 2; its equation is typeset as an image in the retrieved PDF, and the agreement of both parameters to four figures confirms the form). Unweighted fit, as in the source (one specimen per point, no uncertainties given). The paper states its fit parameters to four significant figures; tolerance 0.1 %. Metropolis sampling switched off (it does not enter the estimate). What is tested is the reduction of published measurements, not the physical adequacy of the Spriggs relation, which the paper itself judges too rigid.';
export default [
  // ----- solution verification -----
  { id: 's24-gci-roache-order', suite: 24, title: 'Roache three-grid example: observed order of accuracy', mode: 'gci', inputs: ROACHE, kpi: 'order', expected: R_P, tol: 1e-6, kind: 'numerical-reference', reference: ROACHE_REF,
    note: 'p = ln[(f₃ − f₂)/(f₂ − f₁)]/ln r evaluated here from the published solutions 0.97050 / 0.96854 / 0.96178; the published value is 1.786.' },
  { id: 's24-gci-roache-extrapolate', suite: 24, title: 'Roache three-grid example: Richardson-extrapolated value', mode: 'gci', inputs: ROACHE, kpi: 'f_ext', expected: R_FEXT, tol: 1e-6, kind: 'numerical-reference', reference: ROACHE_REF,
    note: 'f_ext = f₁ + (f₁ − f₂)/(rᵖ − 1); the published value is 0.97130.' },
  { id: 's24-gci-roache-gci-fine', suite: 24, title: 'Roache three-grid example: fine-grid GCI', mode: 'gci', inputs: ROACHE, kpi: 'gci_fine', expected: R_GCI12, tol: 1e-6, kind: 'numerical-reference', reference: ROACHE_REF,
    note: 'GCI₁₂ = Fₛ·|(f₁ − f₂)/f₁|/(rᵖ − 1) with Fₛ = 1.25; the published value is 0.103 %.' },
  { id: 's24-gci-roache-gci-medium', suite: 24, title: 'Roache three-grid example: medium-grid GCI', mode: 'gci', inputs: ROACHE, kpi: 'gci_med', expected: R_GCI23, tol: 1e-6, kind: 'numerical-reference', reference: ROACHE_REF,
    note: 'GCI₂₃ = Fₛ·|(f₂ − f₃)/f₂|/(rᵖ − 1); the published value is 0.356 %.' },
  { id: 's24-gci-roache-asymptotic', suite: 24, title: 'Roache three-grid example: asymptotic-range ratio', mode: 'gci', inputs: ROACHE, kpi: 'asy', expected: RF[0] / RF[1], tol: 1e-6, kind: 'analytical', reference: ROACHE_REF,
    note: 'With the observed order from the same three solutions, GCI₂₃/(rᵖ·GCI₁₂) reduces exactly to f₁/f₂ (1.002 for this example).' },
  { id: 's24-gci-nonuniform-order', suite: 24, title: 'Exact second-order error model on unequal refinement ratios (1.5 and 2): observed order', mode: 'gci', inputs: NONUNI, kpi: 'order', expected: 2, tol: 1e-6, kind: 'analytical', reference: CELIK_REF,
    note: 'Solutions sampled from f(h) = 2 + 0.3·h² at h = 1, 1.5, 3; the fixed-point iteration must return p = 2 exactly.' },
  { id: 's24-gci-nonuniform-extrapolate', suite: 24, title: 'Exact second-order error model on unequal refinement ratios: extrapolated value', mode: 'gci', inputs: NONUNI, kpi: 'f_ext', expected: 2, tol: 1e-6, kind: 'analytical', reference: CELIK_REF,
    note: 'Richardson extrapolation of f(h) = 2 + 0.3·h² recovers the zero-grid-size value 2.' },
  { id: 's24-gci-third-order', suite: 24, title: 'Exact third-order error model, refinement ratio 1.5: observed order', mode: 'gci', inputs: THIRD, kpi: 'order', expected: 3, tol: 1e-6, kind: 'analytical', reference: CELIK_REF,
    note: 'Solutions sampled from f(h) = 5 − 40·h³ at h = 0.1, 0.15, 0.225.' },
  { id: 's24-gci-third-order-gci', suite: 24, title: 'Exact third-order error model: fine-grid GCI equals 1.25 × the true relative error', mode: 'gci', inputs: THIRD, kpi: 'gci_fine', expected: 100 * 1.25 * 40 * 0.1 ** 3 / pw(5, -40, 3, 0.1), tol: 1e-6, kind: 'analytical', reference: CELIK_REF,
    note: 'For an exact power-law error the Richardson error estimate |f₁ − f₂|/(rᵖ − 1) equals the true error C·h₁ᵖ, so GCI = Fₛ·C·h₁ᵖ/|f₁|.' },

  // ----- code verification -----
  { id: 's24-mms-space-order', suite: 24, title: 'Heat equation, central differences: observed spatial order', mode: 'mms', inputs: MMS, kpi: 'p_space', expected: 2, abs: 0.05, kind: 'analytical', reference: VN_REF, note: ORDER_NOTE },
  { id: 's24-mms-cn-order', suite: 24, title: 'Heat equation, Crank–Nicolson: observed temporal order', mode: 'mms', inputs: MMS, kpi: 'p_time', expected: 2, abs: 0.05, kind: 'analytical', reference: VN_REF, note: ORDER_NOTE },
  { id: 's24-mms-be-order', suite: 24, title: 'Heat equation, backward Euler: observed temporal order', mode: 'mms', inputs: { ...MMS, scheme: 'be' }, kpi: 'p_time', expected: 1, abs: 0.05, kind: 'analytical', reference: VN_REF, note: ORDER_NOTE },
  { id: 's24-mms-rk4-order', suite: 24, title: 'Logistic equation, classical Runge–Kutta: observed order', mode: 'mms', inputs: MMS, kpi: 'p_rk4', expected: 4, abs: 0.1, kind: 'analytical', reference: 'Classical fourth-order Runge–Kutta method, global error O(h⁴) — Hairer, Nørsett & Wanner, Solving Ordinary Differential Equations I; exact logistic solution y = 1/(1 + 9·e^(−t))', note: ORDER_NOTE },
  { id: 's24-mms-bdf1-order', suite: 24, title: 'Logistic equation, implicit Euler: observed order', mode: 'mms', inputs: MMS, kpi: 'p_bdf1', expected: 1, abs: 0.05, kind: 'analytical', reference: 'Implicit (backward) Euler method, global error O(h) — Hairer, Nørsett & Wanner, Solving Ordinary Differential Equations I', note: ORDER_NOTE },
  { id: 's24-mms-manufactured-space-order', suite: 24, title: 'Manufactured solution with time-periodic source: observed spatial order', mode: 'mms', inputs: { ...MMS, vcase: 'mms' }, kpi: 'p_space', expected: 2, abs: 0.05, kind: 'analytical', reference: 'Method of manufactured solutions for order-of-accuracy testing — Roache, Verification and Validation in Computational Science and Engineering; Oberkampf & Roy, Verification and Validation in Scientific Computing', note: ORDER_NOTE },
  { id: 's24-mms-manufactured-be-order', suite: 24, title: 'Manufactured solution with time-periodic source: backward-Euler temporal order', mode: 'mms', inputs: { ...MMS, vcase: 'mms', scheme: 'be' }, kpi: 'p_time', expected: 1, abs: 0.1, kind: 'analytical', reference: 'Method of manufactured solutions for order-of-accuracy testing — Roache, Verification and Validation in Computational Science and Engineering; Oberkampf & Roy, Verification and Validation in Scientific Computing', note: ORDER_NOTE },
  { id: 's24-mms-discretisation-error', suite: 24, title: 'Heat equation, 40 intervals: discretisation error of the decaying sine mode', mode: 'mms', inputs: MMS, kpi: 'err_fine', expected: ERR_FINE, tol: 1e-4, kind: 'analytical', reference: VN_REF,
    note: 'Closed form: the discrete amplitude is [(1 − z/2)/(1 + z/2)]^4000 with z = λ_h·Δt, λ_h = (4α/h²)·sin²(πh/2L); the error is taken as the RMS over the interior nodes relative to exp(−π²·Fo).' },
  { id: 's24-mms-integral-error', suite: 24, title: 'Heat equation, 40 intervals: true error of ∫u² dx', mode: 'mms', inputs: MMS, kpi: 'qoi_err', expected: QOI_ERR, tol: 1e-4, kind: 'analytical', reference: VN_REF,
    note: 'The rectangle sum of sin² over a full half-wave is exact, so the error of ∫u² dx is that of the squared discrete amplitude.' },
  { id: 's24-mms-cn-time-error', suite: 24, title: 'Crank–Nicolson, 64 steps: temporal error against the 1024-step reference', mode: 'mms', inputs: MMS, kpi: 'err_dt', expected: ERR_DT_CN, tol: 1e-4, kind: 'analytical', reference: VN_REF,
    note: 'Difference of the Crank–Nicolson amplification products for 64 and 1024 steps on the 20-interval grid, relative to the exact amplitude.' },
  { id: 's24-mms-be-time-error', suite: 24, title: 'Backward Euler, 64 steps: temporal error against the Crank–Nicolson reference', mode: 'mms', inputs: { ...MMS, scheme: 'be' }, kpi: 'err_dt', expected: ERR_DT_BE, tol: 1e-4, kind: 'analytical', reference: VN_REF,
    note: 'Backward-Euler amplification (1 + z)^(−64) against the 1024-step Crank–Nicolson product on the 20-interval grid, relative to the exact amplitude.' },

  // ----- calibration -----
  ...recovery,
  { id: 's24-calib-ols-slope', suite: 24, title: 'Straight-line fit of six points: least-squares slope', mode: 'calib', inputs: OLS, kpi: 'p1', expected: L_B, tol: 1e-6, kind: 'analytical', reference: OLS_REF, note: 'b = Sxy/Sxx for x = 1…6, y = 3.1, 4.9, 7.2, 8.8, 11.1, 13.2 (unweighted).' },
  { id: 's24-calib-ols-intercept', suite: 24, title: 'Straight-line fit of six points: least-squares intercept', mode: 'calib', inputs: OLS, kpi: 'p0', expected: L_A, tol: 1e-6, kind: 'analytical', reference: OLS_REF, note: 'a = ȳ − b·x̄.' },
  { id: 's24-calib-ols-rmse', suite: 24, title: 'Straight-line fit of six points: root-mean-square residual', mode: 'calib', inputs: OLS, kpi: 'rmse', expected: Math.sqrt(L_SSE / LN), tol: 1e-5, kind: 'analytical', reference: OLS_REF, note: 'RMSE = √(SSE/n) of the least-squares line.' },
  { id: 's24-calib-ols-residual-variance', suite: 24, title: 'Straight-line fit of six points: residual variance s² = SSE/(n − 2)', mode: 'calib', inputs: OLS, kpi: 'chi2red', expected: L_S2, tol: 1e-5, kind: 'analytical', reference: OLS_REF, note: 'Without a sigma column the fit is unweighted, so the reduced χ² equals the residual variance.' },
  { id: 's24-calib-ols-standard-error', suite: 24, title: 'Straight-line fit of six points: largest relative standard error', mode: 'calib', inputs: OLS, kpi: 'worst_rel', expected: 100 * Math.max(L_SEA / Math.abs(L_A), L_SEB / Math.abs(L_B)), tol: 1e-4, kind: 'analytical', reference: OLS_REF,
    note: 'se(b) = s/√Sxx, se(a) = s·√(Σx²/(n·Sxx)); the tolerance covers the finite-difference sensitivities.' },
  { id: 's24-calib-ols-correlation', suite: 24, title: 'Straight-line fit of six points: intercept–slope correlation', mode: 'calib', inputs: OLS, kpi: 'maxcorr', expected: L_RHO, tol: 1e-5, kind: 'analytical', reference: OLS_REF, note: '|ρ(a, b)| = Σx/√(n·Σx²) = 21/√546.' },
  { id: 's24-calib-ols-collinearity', suite: 24, title: 'Straight-line fit of six points: collinearity index', mode: 'calib', inputs: OLS, kpi: 'gamma', expected: 1 / Math.sqrt(1 - L_RHO), tol: 1e-5, kind: 'analytical', reference: 'Collinearity index γ = 1/√λ_min of the unit-column sensitivity matrix (Brun, Reichert & Künsch, practical identifiability analysis of large environmental simulation models)',
    note: 'For two parameters the normalised sensitivity Gram matrix is [[1, c], [c, 1]] with c = Σx/√(n·Σx²), so λ_min = 1 − c.' },
  { id: 's24-calib-ols-condition', suite: 24, title: 'Straight-line fit of six points: condition number of the scaled Fisher information', mode: 'calib', inputs: OLS, kpi: 'cond', expected: (1 + L_RHO) / (1 - L_RHO), tol: 1e-5, kind: 'analytical', reference: OLS_REF, note: 'Eigenvalues of [[1, c], [c, 1]] are 1 ± c.' },

  // ----- validation -----
  { id: 's24-valid-nrmse', suite: 24, title: 'Five-point validation table: normalised RMSE', mode: 'valid', inputs: VALID, kpi: 'nrmse', expected: 100 * Math.sqrt(sum(VE.map(v => v * v)) / VE.length) / VMEAN, tol: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'Errors 1, −0.5, 2, −2, 1 on measurements 10…50: RMSE = √2.05, mean measured value 30.' },
  { id: 's24-valid-bias', suite: 24, title: 'Five-point validation table: bias', mode: 'valid', inputs: VALID, kpi: 'bias', expected: 0.3, tol: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'Mean signed error (1 − 0.5 + 2 − 2 + 1)/5.' },
  { id: 's24-valid-mae', suite: 24, title: 'Five-point validation table: mean absolute error', mode: 'valid', inputs: VALID, kpi: 'mae', expected: 1.3, tol: 1e-9, kind: 'analytical', reference: VAL_REF },
  { id: 's24-valid-r2', suite: 24, title: 'Five-point validation table: coefficient of determination', mode: 'valid', inputs: VALID, kpi: 'r2', expected: 1 - 10.25 / 1000, tol: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'R² = 1 − SS_res/SS_tot = 1 − 10.25/1000.' },
  { id: 's24-valid-within', suite: 24, title: 'Five-point validation table: fraction within k·u_c', mode: 'valid', inputs: VALID, kpi: 'within', expected: 60, abs: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'u_c = √(0.6² + 0.8²) = 1 and k = 1.5: three of the five errors lie within ±1.5.' },
  { id: 's24-valid-chi2', suite: 24, title: 'Five-point validation table: reduced χ²', mode: 'valid', inputs: VALID, kpi: 'chi2', expected: sum(VE.map(v => (v / VUC) ** 2)) / VE.length, tol: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'Mean squared error in units of the combined standard uncertainty: 10.25/5.' },
  { id: 's24-valid-extrapolation', suite: 24, title: 'Five-point validation table: extrapolation of the application point', mode: 'valid', inputs: VALID, kpi: 'extrap', expected: 50, tol: 1e-9, kind: 'analytical', reference: VAL_REF, note: 'Application at x = 7 against a validated range 1…5: (7 − 5)/(5 − 1).' },

  // ----- uncertainty propagation and global sensitivity -----
  { id: 's24-uq-npv-annuity', suite: 24, title: 'Net present value at nominal inputs: annuity closed form', mode: 'uq', inputs: { ...NPV, uqOverride: '', cvScale: 0, nSamp: 500 }, kpi: 'mean', expected: NPV0, tol: 1e-9, kind: 'analytical', reference: 'Present value of an ordinary annuity, a = [1 − (1 + r)⁻ⁿ]/r (6.1446 for 10 % and 10 years) — e.g. Brealey, Myers & Allen, Principles of Corporate Finance',
    note: 'All input uncertainty switched off, so every sample equals the nominal value −2·10⁶ + 18·20000·a.' },
  { id: 's24-uq-linear-mean', suite: 24, title: 'Linear additive model: sampled mean', mode: 'uq', inputs: NPV, kpi: 'mean', expected: NPV0, abs: 3 * NSD / Math.sqrt(NSAMP), kind: 'analytical', reference: LIN_REF,
    note: 'NPV with volume and discount rate fixed is linear in CAPEX, price and unit cost, so E[Y] is the nominal value. The absolute tolerance is three standard errors σ/√n of simple random sampling (n = 20000); Latin-hypercube sampling is expected to do better.' },
  { id: 's24-uq-linear-std', suite: 24, title: 'Linear additive model: sampled standard deviation', mode: 'uq', inputs: NPV, kpi: 'std', expected: NSD, tol: 3 / Math.sqrt(2 * NSAMP), kind: 'analytical', reference: LIN_REF,
    note: 'σ_Y² = σ_CAPEX² + (Q·a)²·(σ_price² + σ_cost²). The tolerance is three standard errors 1/√(2n) of a sample standard deviation of normal data.' },
  { id: 's24-uq-linear-pce-std', suite: 24, title: 'Linear additive model: polynomial-chaos standard deviation', mode: 'uq', inputs: NPV, kpi: 'pce_sd', expected: NSD, tol: 1e-5, kind: 'analytical', reference: LIN_REF,
    note: 'A linear function of normal inputs lies in the span of the degree-1 Hermite basis, so the regression surrogate reproduces the variance without sampling error.' },
  { id: 's24-uq-linear-sobol-top', suite: 24, title: 'Linear additive model: total Sobol index of the dominant input (price)', mode: 'uq', inputs: NPV, kpi: 's_top', expected: NV[1] / sum(NV), tol: 0.03, kind: 'analytical', reference: LIN_REF,
    note: 'Variance share of the price term. The Jansen estimator on 24 000 base samples (the suite caps the Saltelli scheme at 120 000 model evaluations) has a relative standard error of about √(2/N) ≈ 0.9 %; the tolerance is three standard errors of this statistical estimate.' },
  { id: 's24-uq-linear-sobol-sum', suite: 24, title: 'Linear additive model: first-order Sobol indices sum to one', mode: 'uq', inputs: NPV, kpi: 'sum_s1', expected: 1, abs: 0.06, kind: 'conservation', reference: LIN_REF,
    note: 'Without interactions the first-order indices exhaust the variance. The sum of the three Saltelli estimates on 6000 base samples scatters by about 0.02 (one standard deviation) between seeds; the absolute tolerance is three such standard deviations.' },
  { id: 's24-uq-linear-exceedance', suite: 24, title: 'Linear additive model: probability of exceeding 1.2 × nominal', mode: 'uq', inputs: NPV, kpi: 'pfail', expected: NPF, abs: 3 * Math.sqrt(0.25 / NSAMP), kind: 'analytical', reference: LIN_REF,
    note: 'Y is normal, so P(Y > 1.2·μ) = 1 − Φ(0.2·μ/σ). The absolute tolerance is three binomial standard errors for n = 20000.' },
  { id: 's24-uq-lognormal-mean', suite: 24, title: 'Product of lognormal inputs: mean', mode: 'uq', inputs: LOGN, kpi: 'mean', expected: L0, tol: 3 * Math.sqrt(LCV2 / NSAMP), kind: 'analytical', reference: LOGN_REF,
    note: 'Arrhenius growth length with Eₐ and T fixed is the product A·k·t of independent lognormal inputs. Tolerance: three standard errors CV/√n of simple random sampling.' },
  { id: 's24-uq-lognormal-cv', suite: 24, title: 'Product of lognormal inputs: coefficient of variation', mode: 'uq', inputs: LOGN, kpi: 'cv', expected: 100 * Math.sqrt(LCV2), tol: 3 / Math.sqrt(2 * NSAMP), kind: 'analytical', reference: LOGN_REF,
    note: 'CV = √[(1 + 0.10²)(1 + 0.02²) − 1] = 10.20 %. Tolerance: three standard errors of a sample standard deviation.' },
  { id: 's24-uq-lognormal-median', suite: 24, title: 'Product of lognormal inputs: median', mode: 'uq', inputs: LOGN, kpi: 'p50', expected: L0 / Math.sqrt(1 + LCV2), tol: 3 * 1.2533 * Math.sqrt(LCV2 / NSAMP), kind: 'analytical', reference: LOGN_REF,
    note: 'Median of a lognormal variable = mean/√(1 + CV²). Tolerance: three standard errors 1.2533·CV/√n of a sample median.' },
  { id: 's24-uq-gibson-ashby-nominal', suite: 24, title: 'Bending-dominated lattice at 20 % relative density: nominal modulus', mode: 'uq', inputs: { mode: 'uq', umodel: 'gibson', cvScale: 0, nSamp: 500, uqOverride: '' }, kpi: 'mean', expected: 2.0e9 * 0.2 ** 2, tol: 1e-9, kind: 'analytical', reference: 'Gibson & Ashby, Cellular Solids — open-cell bending-dominated scaling E*/Eₛ = C·(ρ*/ρₛ)² with C ≈ 1',
    note: 'All input uncertainty switched off: the propagated mean must equal the deterministic model value 80 MPa.' },
  { id: 's24-uq-dlvo-barrier-nominal', suite: 24, title: 'DLVO energy barrier of two 50 nm spheres (ψ = 50 mV, 1 mM 1:1 electrolyte, A_H = 10⁻¹⁹ J)', mode: 'uq', inputs: { mode: 'uq', umodel: 'dlvo', cvScale: 0, nSamp: 500, uqOverride: '', dist: 'normal' }, kpi: 'mean', expected: DLVO0, tol: 0.01, kind: 'analytical', reference: 'DLVO theory with the Derjaguin approximation for equal spheres: V = −A_H·R/(12h) + 2π·ε·R·ψ²·ln(1 + e^(−κh)) — Israelachvili, Intermolecular and Surface Forces; Hunter, Foundations of Colloid Science',
    note: 'Reference maximised here on a 200 000-point separation grid with CODATA constants (Debye length 9.6 nm at 1 mM, 25 °C); the suite searches a 70-point logarithmic grid, hence the 1 % tolerance.' },

  // ----- credibility assessment -----
  { id: 's24-cred-no-validation', suite: 24, title: 'Gate sequence without independent validation data stops at gate 5', mode: 'cred', inputs: { mode: 'cred', nVal: 0, g7uq: true, g7sens: true, g8use: true, g8rev: true }, kpi: 'level', expected: 5, abs: 0, kind: 'analytical', reference: 'Sequential verification–calibration–validation logic: a later gate cannot be credited while an earlier one fails — Oberkampf & Roy, Verification and Validation in Scientific Computing',
    note: 'Rule check, exact by definition: evidence for gates 7 and 8 is declared but gate 6 (independent validation) has no dataset.' },
  { id: 's24-cred-all-gates', suite: 24, title: 'Gate sequence with complete evidence reaches gate 8', mode: 'cred', inputs: { mode: 'cred', nVal: 2, g7uq: true, g7sens: true, g8use: true, g8rev: true }, kpi: 'level', expected: 8, abs: 0, kind: 'analytical', reference: 'Sequential verification–calibration–validation logic — Oberkampf & Roy, Verification and Validation in Scientific Computing',
    note: 'Rule check, exact by definition: every gate has sufficient declared evidence within the default acceptance limits.' },
  { id: 's24-cred-gci-blocks', suite: 24, title: 'A grid convergence index above its limit stops the sequence at gate 3', mode: 'cred', inputs: { mode: 'cred', nVal: 2, gciPct: 5, gciTolC: 2, g7uq: true, g7sens: true, g8use: true, g8rev: true }, kpi: 'level', expected: 3, abs: 0, kind: 'analytical', reference: 'Sequential verification–calibration–validation logic — Oberkampf & Roy, Verification and Validation in Scientific Computing',
    note: 'Rule check, exact by definition: solution verification (gate 4) fails, so calibration and validation evidence is not credited.' },
  // ------------------------------------------------------------------ published measured data sets (NIST StRD)
  { id: 's24-exp-nist-norris-ozone-slope', suite: 24, title: 'NIST StRD Norris (calibration of ozone monitors, 36 observations): slope of the calibration line', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'linear', cdata: STRD_NORRIS, bayes: false }, kpi: 'p1', expected: 1.00211681802045, tol: 1e-6, kind: 'experiment',
    measured: { value: 1.00211681802045, uncertainty: 0.429796848199937e-3, unit: '–' },
    source: strdSrc('lls/data/Norris.shtml', 'Norris — linear least squares regression, lower level of difficulty'),
    conditions: 'Observed data from the calibration of ozone monitors (J. Norris, NIST): 36 pairs of readings between 0.1 and 999; straight line y = B0 + B1·x fitted by ordinary least squares',
    match: ['cmodel'],
    reference: 'NIST Statistical Reference Datasets, Norris (J. Norris, NIST, Calibration of Ozone Monitors) — certified regression statistics',
    note: 'Certified value B1 = 1.00211681802045 (standard deviation of the estimate 0.429796848199937E-03) from the "Certified Regression Statistics" block of Norris.dat; the 36 observations are the "Data" block of the same file. The suite also returns the certified intercept B0 = −0.262323073774029 to 4 parts in 10⁷. ' + STRD_NOTE },
  { id: 's24-exp-nist-norris-ozone-residual-variance', suite: 24, title: 'NIST StRD Norris (calibration of ozone monitors): residual mean square of the straight-line fit', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'linear', cdata: STRD_NORRIS, bayes: false }, kpi: 'chi2red', expected: 0.782864662630069, tol: 1e-6, kind: 'experiment',
    measured: { value: 0.782864662630069, uncertainty: null, unit: '(reading units)²' },
    source: strdSrc('lls/data/Norris.shtml', 'Norris — linear least squares regression, lower level of difficulty'),
    conditions: 'Observed data from the calibration of ozone monitors (J. Norris, NIST): 36 observations, 34 residual degrees of freedom',
    match: ['cmodel'],
    reference: 'NIST Statistical Reference Datasets, Norris — certified analysis-of-variance table',
    note: 'Residual mean square 0.782864662630069 (residual sum of squares 26.6173985294224 on 34 degrees of freedom) from the "Certified Analysis of Variance Table" of Norris.dat; its square root is the certified residual standard deviation 0.884796396144373. Without a sigma column the suite\'s reduced χ² is the residual sum of squares divided by the degrees of freedom, i.e. this mean square. ' + STRD_NOTE },
  { id: 's24-exp-nist-danwood-radiation-exponent', suite: 24, title: 'NIST StRD DanWood (energy radiated by a carbon-filament lamp against temperature): power-law exponent', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'power', cdata: STRD_DANWOOD, bayes: false }, kpi: 'p1', expected: 3.8604055871, tol: 1e-6, kind: 'experiment',
    measured: { value: 3.8604055871, uncertainty: 5.1726610913e-2, unit: '–' },
    source: strdSrc('nls/data/daniel_wood.shtml', 'DanWood — nonlinear least squares regression, lower level of difficulty'),
    conditions: 'Energy radiated from a carbon-filament lamp per cm² per second against the absolute temperature of the filament in units of 1000 K (1.309–1.680), 6 observations; data described in Daniel & Wood, Fitting Equations to Data (1980), originally from Keeping, Introduction to Statistical Inference (1962); model y = b1·x^b2',
    match: ['cmodel'],
    reference: 'NIST Statistical Reference Datasets, DanWood (Daniel & Wood, Fitting Equations to Data, 2nd ed., 1980) — certified values',
    note: 'Certified b2 = 3.8604055871E+00 (standard deviation 5.1726610913E-02) from the "Certified Values" block of DanWood.dat; the suite also returns the certified b1 = 7.6886226176E-01 to 1 part in 10⁸. The exponent lies 2.7 standard deviations below the Stefan–Boltzmann value of 4. ' + STRD_NOTE },
  { id: 's24-exp-nist-misra1a-adsorption-rate', suite: 24, title: 'NIST StRD Misra1a (monomolecular adsorption, volume against pressure): saturation-law rate constant', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'firstorder', cdata: STRD_MISRA1A, bayes: false }, kpi: 'p1', expected: 5.5015643181e-4, tol: 1e-6, kind: 'experiment',
    measured: { value: 5.5015643181e-4, uncertainty: 7.2668688436e-6, unit: '1/(pressure unit)' },
    source: strdSrc('nls/data/misra1a.shtml', 'Misra1a — nonlinear least squares regression, lower level of difficulty'),
    conditions: 'NIST study on dental research in monomolecular adsorption (D. Misra, 1978): adsorbed volume against pressure, 14 observations, pressures 77.6–760.0; model y = b1·(1 − exp(−b2·x))',
    match: ['cmodel'],
    reference: 'NIST Statistical Reference Datasets, Misra1a (D. Misra, NIST, Dental Research Monomolecular Adsorption Study, 1978) — certified values',
    note: 'Certified b2 = 5.5015643181E-04 (standard deviation 7.2668688436E-06) from the "Certified Values" block of Misra1a.dat; the suite also returns the certified b1 = 2.3894212918E+02 to 2 parts in 10⁹ and the certified residual standard deviation 1.0187876330E-01. The data cover only the nearly linear start of the saturation curve, so b1 and b2 are strongly correlated (the suite reports a collinearity index of about 29), which is what makes the set a useful test. The source states no units. ' + STRD_NOTE },
  { id: 's24-exp-nist-boxbod-ultimate-demand', suite: 24, title: 'NIST StRD BoxBOD (biochemical oxygen demand against incubation time): ultimate demand', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'firstorder', cdata: STRD_BOXBOD, bayes: false }, kpi: 'p0', expected: 213.80940889, tol: 1e-6, kind: 'experiment',
    measured: { value: 213.80940889, uncertainty: 12.354515176, unit: 'mg/L' },
    source: strdSrc('nls/data/boxbod.shtml', 'BoxBOD — nonlinear least squares regression, higher level of difficulty'),
    conditions: 'Biochemical oxygen demand in mg/L after incubation for 1, 2, 3, 5, 7 and 10 days (6 observations), data of Box, Hunter & Hunter, Statistics for Experimenters (1978); first-order model y = b1·(1 − exp(−b2·x))',
    match: ['cmodel'],
    reference: 'NIST Statistical Reference Datasets, BoxBOD (Box, Hunter & Hunter, Statistics for Experimenters, 1978) — certified values',
    note: 'Certified b1 = 2.1380940889E+02 (standard deviation 1.2354515176E+01) from the "Certified Values" block of BoxBOD.dat; the suite also returns the certified rate constant b2 = 5.4723748542E-01 per day to 2 parts in 10⁸ and the certified residual standard deviation 1.7088072423E+01 mg/L. NIST classes this set as of higher difficulty (six scattered points; the suite reports R² = 0.88). ' + STRD_NOTE },
  // ------------------------------------------------------------------ published measured data set with published fit (porous zirconia)
  { id: 's24-exp-pabst-zirconia-potato-spriggs-e0', suite: 24, title: 'Porous zirconia (potato-starch series): dense-material modulus E₀ from the Spriggs fit of measured modulus against porosity', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'exp', cdata: PABST_ZP, bayes: false }, kpi: 'p0', expected: 233.9, tol: 1e-3, kind: 'experiment',
    measured: { value: 233.9, uncertainty: null, unit: 'GPa' },
    source: PABST_SRC, conditions: PABST_COND + '; potato starch as pore former, porosity 0.132–0.424 (5 specimens)',
    match: ['cmodel'],
    reference: 'Pabst, Gregorová, Tichá & Týnová, Ceramics − Silikáty 48, 165 (2004), Tables 2 and 7',
    note: 'Measured data: Table 2 (Z-P columns). Published fit: Table 7, Model 2 — E₀ = 233.9 GPa, B = 3.890, correlation coefficient 0.990 (the paper remarks "E0 too high" against the 210 GPa of dense zirconia). The suite gives E₀ = 233.85 GPa and B = 3.8904 (R² = 0.980, i.e. r = 0.990). ' + PABST_NOTE },
  { id: 's24-exp-pabst-zirconia-corn-spriggs-b', suite: 24, title: 'Porous zirconia (corn-starch series): Spriggs porosity coefficient from measured modulus against porosity', mode: 'calib',
    inputs: { mode: 'calib', cmodel: 'exp', cdata: PABST_ZC, bayes: false }, kpi: 'p1', expected: -3.118, tol: 1e-3, kind: 'experiment',
    measured: { value: 3.118, uncertainty: null, unit: '–' },
    source: PABST_SRC, conditions: PABST_COND + '; corn starch as pore former, porosity 0.073–0.415 (5 specimens)',
    match: ['cmodel'],
    reference: 'Pabst, Gregorová, Tichá & Týnová, Ceramics − Silikáty 48, 165 (2004), Tables 2 and 9',
    note: 'Measured data: Table 2 (Z-C columns). Published fit: Table 9, Model 2 — E₀ = 216.5 GPa, B = 3.118, correlation coefficient 0.990. The suite parameter is b = −B, hence the expected value −3.118; the suite gives b = −3.1178 and E₀ = 216.49 GPa. ' + PABST_NOTE },
];
