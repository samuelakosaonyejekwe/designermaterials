// Suite 24 — Verification, Calibration, Validation and Uncertainty.
// Richardson/GCI solution verification, code verification by exact/manufactured solutions, weighted nonlinear
// least-squares and Bayesian calibration, identifiability, independent validation metrics, uncertainty propagation,
// global sensitivity (Sobol, Morris, polynomial chaos) and the eight-gate credibility assessment.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS;
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ex = v => Math.exp(N.clamp(v, -700, 700));
const pct = v => F(100 * v, 3) + ' %';

// ---------- table handling ----------
function parseTable(text, minCols, minRows, label) {
  const lines = String(text ?? '').split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#') && !l.startsWith('//'));
  const rows = []; let head = null;
  for (const l of lines) {
    const c = l.split(/[,;\t]/).map(s => s.trim()), v = c.map(s => s === '' ? NaN : Number(s));
    if (v.slice(0, minCols).every(Number.isFinite) && c.length >= minCols) rows.push(v); else if (!head && !rows.length) head = c;
  }
  if (rows.length < minRows) throw new Error(`${label}: at least ${minRows} numeric rows with ${minCols} columns are needed (found ${rows.length}). Paste CSV/TSV text with one header line.`);
  return { head: head || [], rows };
}
const csv = (head, rows) => head.join(',') + '\n' + rows.map(r => r.map(v => +(+v).toPrecision(5)).join(',')).join('\n');

// ---------- calibration model families ----------
const pos = a => a.filter(v => v > 0);
const FAM = {
  arrhenius: { label: 'Arrhenius  y = A·exp(−Eₐ/(R·x)), x = T in K', names: ['ln A', 'Eₐ (kJ/mol)'], eq: 'y = A·exp(−Eₐ/(R·T))', f: (x, p) => ex(p[0] - p[1] * 1e3 / (R * Math.max(x, 1e-9))), b: [[-200, 300], [0, 1000]],
    guess: (x, y) => { const i = y.map((v, k) => k).filter(k => y[k] > 0 && x[k] > 0); if (i.length < 2) return [0, 50]; const l = N.linreg(i.map(k => 1 / x[k]), i.map(k => Math.log(y[k]))); return [l.intercept, -l.slope * R / 1e3]; }, derived: p => [['A = exp(ln A)', ex(p[0])]] },
  power: { label: 'Power law  y = a·xⁿ', names: ['a', 'n'], eq: 'y = a·xⁿ', f: (x, p) => p[0] * Math.pow(Math.max(x, 1e-300), N.clamp(p[1], -50, 50)), b: [[-1e30, 1e30], [-20, 20]],
    guess: (x, y) => { const i = y.map((v, k) => k).filter(k => y[k] > 0 && x[k] > 0); if (i.length < 2) return [N.mean(y) || 1, 1]; const l = N.linreg(i.map(k => Math.log(x[k])), i.map(k => Math.log(y[k]))); return [Math.exp(l.intercept), l.slope]; } },
  exp: { label: 'Exponential decay/growth  y = a·exp(b·x)', names: ['a', 'b'], eq: 'y = a·exp(b·x)  (b < 0 decay, b > 0 growth)', f: (x, p) => p[0] * ex(p[1] * x), b: [[-1e30, 1e30], [-1e6, 1e6]],
    guess: (x, y) => { const i = y.map((v, k) => k).filter(k => y[k] > 0); if (i.length < 2) return [N.mean(y) || 1, 0]; const l = N.linreg(i.map(k => x[k]), i.map(k => Math.log(y[k]))); return [Math.exp(l.intercept), l.slope]; } },
  linear: { label: 'Linear  y = a + b·x', names: ['a', 'b'], eq: 'y = a + b·x', f: (x, p) => p[0] + p[1] * x, b: [[-1e30, 1e30], [-1e30, 1e30]], guess: (x, y) => { const l = N.linreg(x, y); return [l.intercept, l.slope]; } },
  langmuir: { label: 'Langmuir isotherm  y = q_max·K·x/(1 + K·x)', names: ['q_max', 'K'], eq: 'q = q_max·K·c/(1 + K·c)', f: (x, p) => p[0] * p[1] * x / Math.max(1e-12, 1 + p[1] * x), b: [[0, 1e30], [0, 1e30]],
    guess: (x, y) => [1.2 * Math.max(N.maxOf(y), 1e-12), 1 / Math.max(N.quantile(x.map(Math.abs), 0.5), 1e-12)] },
  jmak: { label: 'JMAK transformation  y = 1 − exp(−(k·x)ⁿ)', names: ['k', 'n'], eq: 'X(t) = 1 − exp(−(k·t)ⁿ)', f: (x, p) => 1 - ex(-Math.pow(Math.max(p[0] * x, 0), N.clamp(p[1], 0.05, 20))), b: [[0, 1e30], [0.2, 8]],
    guess: (x, y) => { const i = y.map((v, k) => k).filter(k => y[k] > 0.01 && y[k] < 0.99 && x[k] > 0); if (i.length < 2) return [1 / Math.max(N.quantile(x.map(Math.abs), 0.5), 1e-12), 2]; const l = N.linreg(i.map(k => Math.log(x[k])), i.map(k => Math.log(-Math.log(1 - y[k])))); const n = N.clamp(l.slope, 0.3, 6); return [Math.exp(l.intercept / n), n]; } },
  parabolic: { label: 'Parabolic growth  L = √(k·t)', names: ['k'], eq: 'L = √(k·t)', f: (x, p) => Math.sqrt(Math.max(0, p[0] * x)), b: [[0, 1e30]], guess: (x, y) => { const i = x.map((v, k) => k).filter(k => x[k] > 0); return [i.length ? N.mean(i.map(k => y[k] * y[k] / x[k])) : 1]; } },
  difflim: { label: 'Diffusion-limited growth  r = √(r₀² + k·t)', names: ['r₀', 'k'], eq: 'r² = r₀² + 2·D·Ω·(c∞ − c_eq)·t = r₀² + k·t', f: (x, p) => Math.sqrt(Math.max(0, p[0] * p[0] + p[1] * x)), b: [[0, 1e30], [0, 1e30]],
    guess: (x, y) => { const l = N.linreg(x, y.map(v => v * v)); return [Math.sqrt(Math.max(l.intercept, 1e-12 * Math.abs(l.intercept) + 1e-30)), Math.max(l.slope, 1e-30)]; } },
  firstorder: { label: 'First-order kinetics  y = y∞·(1 − exp(−k·x))', names: ['y∞', 'k'], eq: 'dy/dt = k·(y∞ − y)  ⇒  y = y∞·(1 − exp(−k·t))', f: (x, p) => p[0] * (1 - ex(-p[1] * x)), b: [[-1e30, 1e30], [0, 1e30]],
    guess: (x, y) => [1.1 * (N.maxOf(y) || 1), 2 / Math.max(N.maxOf(x.map(Math.abs)), 1e-12)] },
};
// Student-t 97.5 % quantile: Simpson quadrature of the density (substitution x = √ν·tan θ) inverted by a bracketed root search;
// the Cornish–Fisher expansion in 1/ν is kept for ν > 200, where it is accurate to better than 10⁻⁸
function t975(nu) {
  const v = Math.max(1, nu), z = 1.959964, z3 = z ** 3, z5 = z ** 5, z7 = z ** 7;
  if (v > 200) return z + (z3 + z) / (4 * v) + (5 * z5 + 16 * z3 + 3 * z) / (96 * v * v) + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * v ** 3);
  const c = Math.exp(N.gammaLn((v + 1) / 2) - N.gammaLn(v / 2)) / Math.sqrt(Math.PI);
  const area = t => { const b = Math.atan(t / Math.sqrt(v)), n = 200, h = b / n; let s = 1 + Math.pow(Math.cos(b), v - 1); for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * Math.pow(Math.cos(i * h), v - 1); return c * s * h / 3; };
  return N.brent(t => area(t) - 0.475, 1.9, 13, 1e-10);
}

function calibrate(X, Y, S, fam, scaleCov, fixed) {
  const m = X.length, g0 = fam.guess(X, Y).map((v, j) => N.clamp(Number.isFinite(v) ? v : 1, fam.b[j][0], fam.b[j][1]));
  const sc = g0.map(v => Math.abs(v) > 1e-300 ? Math.abs(v) : 1), n = g0.length;
  const res = p => X.map((x, i) => { const v = (fam.f(x, p) - Y[i]) / S[i]; return Number.isFinite(v) ? v : 1e150; });
  const bq = fam.b.map((b, j) => [b[0] / sc[j], b[1] / sc[j]]);
  const lm = N.levenbergMarquardt(q => res(q.map((v, j) => v * sc[j])), g0.map((v, j) => v / sc[j]), { maxIter: 200, tol: 1e-13, bounds: bq });
  const p = lm.p.map((v, j) => v * sc[j]), chi2 = lm.cost, dof = Math.max(1, m - n);
  // sensitivity matrix (central differences) and Fisher information
  const J = N.zeros(m, n);
  for (let j = 0; j < n; j++) { const h = 1e-5 * Math.max(Math.abs(p[j]), 1e-3 * sc[j]), a = [...p], b = [...p]; a[j] += h; b[j] -= h; for (let i = 0; i < m; i++) { const d = (fam.f(X[i], a) - fam.f(X[i], b)) / (2 * h); J[i][j] = Number.isFinite(d) ? d : 0; } }
  const FI = N.zeros(n, n); for (let i = 0; i < m; i++) for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) FI[a][b] += J[i][a] * J[i][b] / (S[i] * S[i]);
  const D = FI.map((r, j) => r[j] > 0 ? 1 / Math.sqrt(r[j]) : 1), Fs = FI.map((r, a) => r.map((v, b) => v * D[a] * D[b]));
  const eg = N.eigSym(Fs), lmax = Math.max(N.maxOf(eg.values), 1e-300), lmin = Math.max(N.minOf(eg.values), 0);
  const inv = N.zeros(n, n); for (let k = 0; k < n; k++) { const l = Math.max(eg.values[k], 1e-13 * lmax), v = eg.vectors[k]; for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) inv[a][b] += v[a] * v[b] / l; }
  const s2 = scaleCov ? chi2 / dof : 1, cov0 = inv.map((r, a) => r.map((v, b) => v * D[a] * D[b])), cov = cov0.map(r => r.map(v => v * s2));
  const se = cov.map((r, j) => Math.min(Math.sqrt(Math.abs(r[j])), 1e6 * Math.max(Math.abs(p[j]), sc[j])));
  const corr = cov.map((r, a) => r.map((v, b) => a === b ? 1 : N.clamp(v / Math.sqrt(Math.abs(cov[a][a] * cov[b][b]) || 1e-300), -1, 1)));
  const identifiable = lmin > 1e-10 * lmax;
  return { p, se, cov, cov0, corr, chi2, dof, s2, eig: eg.values, cond: lmax / Math.max(lmin, 1e-16 * lmax), gamma: 1 / Math.sqrt(Math.max(lmin, 1e-16 * lmax)), hist: lm.history, iters: lm.iterations, sc, res, n, m, identifiable, pred: X.map(x => fam.f(x, p)), fixed };
}
// re-optimise all parameters except j (profile likelihood)
function profile(X, Y, S, fam, fit, j, grid) {
  const n = fit.n, others = [...Array(n).keys()].filter(k => k !== j);
  return grid.map(v => {
    const full = q => { const p = [...fit.p]; p[j] = v; others.forEach((k, i) => p[k] = q[i] * fit.sc[k]); return p; };
    if (!others.length) { const r = fit.res(full([])); return N.dot(r, r); }
    const lm = N.levenbergMarquardt(q => fit.res(full(q)), others.map(k => fit.p[k] / fit.sc[k]), { maxIter: 40, tol: 1e-11, bounds: others.map(k => [fam.b[k][0] / fit.sc[k], fam.b[k][1] / fit.sc[k]]) });
    return lm.cost;
  });
}

// ---------- default synthetic datasets (forward model + seeded noise) ----------
const TRUE_AR = [Math.log(2.4e9), 60];   // ln A (nm/min), Ea (kJ/mol): ZnO nanorod axial growth-rate demonstration
function synthRows(fn, xs, rel, seed) { const r = N.rng(seed); return xs.map(x => { const y0 = fn(x), s = rel * Math.abs(y0); return [x, y0 + r.normal(0, s), s]; }); }
const CAL_DEF = csv(['T_K', 'growth_rate_nm_min', 'sigma'], synthRows(x => FAM.arrhenius.f(x, TRUE_AR), [333, 338, 343, 348, 353, 358, 363, 368], 0.04, 19));
const VAL_DEF = (() => { const r = N.rng(77); return csv(['T_K', 'measured', 'u_measured', 'predicted', 'u_predicted'], [336, 341, 346, 351, 356, 361, 366, 371].map(x => { const t = FAM.arrhenius.f(x, [Math.log(2.4e9) + 0.33, 61]), s = 0.04 * t, pr = FAM.arrhenius.f(x, TRUE_AR); return [x, t + r.normal(0, s), s, pr, 0.03 * pr]; })); })();
const JMAK_DEF = csv(['time_min', 'fraction_transformed', 'sigma'], synthRows(x => FAM.jmak.f(x, [0.02, 2.5]), [10, 20, 30, 40, 50, 60, 75, 90, 110, 130], 0.0, 5).map((r, i) => { const g = N.rng(500 + i); return [r[0], N.clamp(r[1] + g.normal(0, 0.015), 0.001, 0.999), 0.015]; }));
const LANG_DEF = csv(['c_mg_L', 'q_mg_g', 'sigma'], synthRows(x => FAM.langmuir.f(x, [120, 0.08]), [2, 5, 10, 20, 35, 50, 80, 120, 180], 0.05, 9));

// ---------- uncertainty-propagation models ----------
const HGRID = Array.from({ length: 70 }, (_, i) => 2e-10 * Math.pow(250, i / 69));   // surface separations 0.2–50 nm
const UQM = {
  arrh: { label: 'Arrhenius nanorod growth length  L = A·exp(−Eₐ/RT)·t', out: 'Nanorod length', unit: 'nm', eq: 'L = A·exp(−Eₐ/(R·T))·t',
    vars: [['A', 'nm/s', 4e7, 0.10], ['Ea', 'kJ/mol', 60, 0.02], ['T', 'K', 363, 0.005], ['t', 's', 10800, 0.02]], f: v => v[0] * ex(-v[1] * 1e3 / (R * Math.max(v[2], 1))) * v[3] },
  gibson: { label: 'Gibson–Ashby lattice stiffness  E* = C·Eₛ·ρ̄ⁿ', out: 'Effective lattice modulus', unit: 'Pa', eq: 'E* = C·Eₛ·(ρ*/ρₛ)ⁿ   (n ≈ 2 bending-dominated, n ≈ 1 stretching-dominated)',
    vars: [['C', '–', 1.0, 0.10], ['Es', 'Pa', 2.0e9, 0.05], ['rho_rel', '–', 0.20, 0.05], ['n', '–', 2.0, 0.03]], f: v => v[0] * v[1] * Math.pow(Math.max(v[2], 1e-9), v[3]) },
  srr: { label: 'Split-ring resonator frequency  f₀ = 1/(2π√(LC))', out: 'SRR resonance frequency', unit: 'Hz', eq: 'f₀ = 1/(2π√(L·C)),  L = μ₀R(ln(8R/(w+t)) − ½),  C = ε₀εᵣ·w·t/g + ε₀·½(1+εᵣ)·(2(w+t)/π)·ln(4R/g)',
    vars: [['R', 'm', 3e-3, 0.01], ['w', 'm', 5e-4, 0.03], ['t', 'm', 3.5e-5, 0.05], ['g', 'm', 3e-4, 0.05], ['eps_r', '–', 4.4, 0.03]],
    f: v => { const Rr = Math.max(v[0], 1e-9), w = Math.max(v[1], 1e-12), t = Math.max(v[2], 1e-12), g = Math.max(v[3], 1e-12), L = N.MU0 * Rr * Math.max(Math.log(8 * Rr / (w + t)) - 0.5, 0.05), C = N.EPS0 * v[4] * w * t / g + N.EPS0 * 0.5 * (1 + v[4]) * (2 * (w + t) / Math.PI) * Math.max(Math.log(4 * Rr / g), 0.05); return 1 / (2 * Math.PI * Math.sqrt(Math.max(L * C, 1e-300))); } },
  dlvo: { label: 'DLVO energy barrier of a colloid (equal spheres)', out: 'DLVO barrier height', unit: 'kT', eq: 'V(h) = −A_H·R/(12h) + 2π·ε·R·ψ²·ln(1 + e^(−κh)),  κ = √(2·N_A·e²·I/(ε·k_B·T)),  barrier = max_h V(h)',
    vars: [['A_H', 'J', 1e-19, 0.20], ['psi', 'V', 0.050, 0.10], ['R', 'm', 2.5e-8, 0.10], ['I', 'mol/L', 1e-3, 0.20]],
    f: v => { const T = 298.15, e = 78.5 * N.EPS0, kT = N.K_B * T, kap = Math.sqrt(2 * N.N_A * N.Q_E ** 2 * Math.max(v[3], 1e-12) * 1e3 / (e * kT)); let mx = -Infinity; for (let i = 0; i < 70; i++) { const h = HGRID[i], V = -v[0] * v[2] / (12 * h) + 2 * Math.PI * e * v[2] * v[1] * v[1] * Math.log(1 + Math.exp(-kap * h)); if (V > mx) mx = V; } return mx / kT; } },
  npv: { label: 'Net present value of a 10-year production project', out: 'Net present value', unit: 'USD', eq: 'NPV = −CAPEX + Σₜ₌₁¹⁰ (p − c)·Q/(1 + r)ᵗ',
    vars: [['CAPEX', 'USD', 2e6, 0.15], ['price', 'USD/kg', 40, 0.10], ['volume', 'kg/yr', 20000, 0.15], ['unit_cost', 'USD/kg', 22, 0.10], ['r', '–', 0.10, 0.15]],
    f: v => { let s = -v[0]; const r = Math.max(v[4], -0.9); for (let t = 1; t <= 10; t++) s += (v[1] - v[3]) * v[2] / Math.pow(1 + r, t); return s; } },
};

// ---------- gate definitions (Part VI of the CVVUQ framework) ----------
const GATES = ['Input data quality', 'Mathematical and model consistency', 'Code verification', 'Solution verification', 'Parameter calibration', 'Independent experimental validation', 'Predictive robustness', 'Application-specific qualification'];
const CLASSES = ['Not verified', 'Numerically verified for the specified problem class', 'Calibrated against identified datasets', 'Experimentally validated over a documented domain', 'Qualified for a particular engineering application following the required review'];

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'calib', group: 'Model', options: [
    { value: 'calib', label: 'Calibration and identifiability (least squares + Bayesian)' }, { value: 'valid', label: 'Independent validation against measurements' },
    { value: 'gci', label: 'Solution verification: Richardson extrapolation and GCI' }, { value: 'mms', label: 'Code verification: exact/manufactured solutions (self-run)' },
    { value: 'uq', label: 'Uncertainty propagation and global sensitivity' }, { value: 'cred', label: 'Credibility assessment: eight qualification gates' }],
    help: 'Selects which credibility activity is carried out. Calibration, verification and validation are separate activities: agreement with the data used for fitting is not validation.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 99999, step: 1, sens: false, group: 'Model', help: 'Seed of the pseudo-random generator used for Monte Carlo, Latin hypercube and MCMC sampling. The same seed reproduces the same result exactly.' },
  // calibration
  { key: 'cmodel', label: 'Model family', type: 'select', value: 'arrhenius', group: 'Calibration', showIf: { mode: ['calib', 'valid'] }, options: Object.entries(FAM).map(([value, m]) => ({ value, label: m.label })), help: 'Mathematical model fitted to the calibration table (x = first column, y = second column). Choose the family that represents the governing mechanism, not the one that merely gives the smallest residual.' },
  { key: 'cdata', label: 'Calibration data (x, y, sigma)', type: 'table', columns: ['x', 'y', 'sigma'], value: CAL_DEF, group: 'Calibration', showIf: { mode: ['calib', 'valid'] }, help: 'Observations used for parameter estimation: x, y and (optional) one-standard-deviation measurement uncertainty of y. Without a sigma column the fit is unweighted and the noise level is estimated from the residuals.' },
  { key: 'scaleCov', label: 'Scale covariance by reduced χ²', type: 'bool', value: true, group: 'Calibration', showIf: { mode: ['calib', 'valid'] }, help: 'When on, parameter covariance is multiplied by χ²/(m − n) so that mis-stated measurement uncertainties do not produce over-confident parameter errors. Switch off only when the sigmas are trusted absolutely.' },
  { key: 'bayes', label: 'Bayesian calibration (Metropolis MCMC)', type: 'bool', value: true, group: 'Calibration', showIf: { mode: ['calib'] }, help: 'Adds a random-walk Metropolis sampling of the posterior with uniform priors inside the physical parameter bounds, giving credible intervals and a posterior predictive band.' },
  { key: 'nMc', label: 'MCMC samples', value: 6000, min: 500, max: 30000, step: 500, sens: false, group: 'Calibration', showIf: { mode: ['calib'] }, help: 'Length of the Markov chain; the first quarter is discarded as burn-in. Longer chains reduce Monte Carlo error of the credible intervals.' },
  { key: 'rxn', label: 'Reaction being calibrated (optional)', type: 'text', value: 'Zn(NO3)2 + NaOH -> ZnO + NaNO3 + H2O', group: 'Calibration', showIf: { mode: ['calib'] }, help: 'Optional unbalanced reaction "A + B -> C + D" whose kinetics the data describe. It is balanced and reported with the fit so the calibrated record names the chemistry.' },
  // validation
  { key: 'vdata', label: 'Validation data (x, measured, u_measured, predicted, u_predicted)', type: 'table', columns: ['x', 'measured', 'u_measured', 'predicted', 'u_predicted'], value: VAL_DEF, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Independent measurements not used in calibration, with one-sigma uncertainties. The predicted columns are optional when the prediction is taken from the calibrated model.' },
  { key: 'predSrc', label: 'Source of predictions', type: 'select', value: 'model', group: 'Validation', showIf: { mode: ['valid'] }, options: [{ value: 'model', label: 'Model calibrated on the calibration table' }, { value: 'table', label: 'Predicted columns of the validation table' }], help: 'Either predict at the validation conditions with the model fitted to the calibration table (prediction uncertainty from the parameter covariance), or compare externally computed predictions supplied in the table.' },
  { key: 'kCov', label: 'Coverage factor k', value: 2, min: 1, max: 3, step: 0.1, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Multiplier of the combined standard uncertainty used for the agreement test |error| ≤ k·u_c. k = 2 corresponds to about 95 % for normal errors.' },
  { key: 'accRmse', label: 'Acceptance: normalised RMSE', unit: '%', value: 10, min: 0.1, max: 100, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Largest root-mean-square error, as a percentage of the mean measured magnitude, that the intended application tolerates. This is an application decision, not a universal threshold.' },
  { key: 'accBias', label: 'Acceptance: |bias|', unit: '%', value: 5, min: 0.1, max: 100, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Largest systematic deviation (mean predicted minus measured, relative to mean measured magnitude) tolerated by the application.' },
  { key: 'accCov', label: 'Acceptance: fraction within k·u_c', unit: '%', value: 90, min: 10, max: 100, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Minimum percentage of validation points whose error must lie inside the combined measurement-plus-prediction uncertainty band.' },
  { key: 'xApp', label: 'Application condition x', value: 380, min: -1e6, max: 1e6, group: 'Validation', showIf: { mode: ['valid'] }, help: 'Value of the independent variable at which the model is to be used. It is compared with the range covered by validation data to quantify extrapolation.' },
  // GCI
  { key: 'gh1', label: 'Grid size h₁ (fine)', value: 1, min: 1e-6, max: 1e3, log: true, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Representative cell or time-step size of the finest solution, in any consistent unit (only ratios matter).' },
  { key: 'gh2', label: 'Grid size h₂ (medium)', value: 2, min: 1e-6, max: 1e3, log: true, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Representative cell or time-step size of the medium solution. A refinement ratio above 1.3 is recommended.' },
  { key: 'gh3', label: 'Grid size h₃ (coarse)', value: 4, min: 1e-6, max: 1e3, log: true, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Representative cell or time-step size of the coarsest solution.' },
  { key: 'gf1', label: 'Result f₁ on fine grid', value: 0.97050, min: -1e9, max: 1e9, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Quantity of interest (peak stress, resonance frequency, conversion…) computed on the finest grid.' },
  { key: 'gf2', label: 'Result f₂ on medium grid', value: 0.96854, min: -1e9, max: 1e9, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Same quantity of interest computed on the medium grid.' },
  { key: 'gf3', label: 'Result f₃ on coarse grid', value: 0.96178, min: -1e9, max: 1e9, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Same quantity of interest computed on the coarsest grid.' },
  { key: 'gFs', label: 'Safety factor Fₛ', value: 1.25, min: 1, max: 3, step: 0.05, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Roache safety factor: 1.25 when three grids give an observed order, 3 when only two grids and an assumed order are available.' },
  { key: 'gOrder', label: 'Formal order of the scheme', value: 2, min: 1, max: 6, step: 1, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Theoretical order of accuracy of the discretisation; the observed order is compared with it and the conservative GCI caps the order at this value.' },
  { key: 'gTol', label: 'Acceptance: fine-grid GCI', unit: '%', value: 1, min: 0.001, max: 50, log: true, group: 'Solution verification', showIf: { mode: ['gci'] }, help: 'Largest numerical uncertainty of the fine-grid result that the application tolerates.' },
  // MMS
  { key: 'vcase', label: 'Verification problem', type: 'select', value: 'decay', group: 'Code verification', showIf: { mode: ['mms'] }, options: [{ value: 'decay', label: 'Exact solution: decaying sine mode of the heat equation' }, { value: 'mms', label: 'Manufactured solution with time-periodic source' }], help: 'Transient 1-D heat/diffusion problem with homogeneous Dirichlet boundaries. The exact case has no source; the manufactured case adds the source term that makes a chosen function an exact solution.' },
  { key: 'scheme', label: 'Time integration tested', type: 'select', value: 'cn', group: 'Code verification', showIf: { mode: ['mms'] }, options: [{ value: 'cn', label: 'Crank–Nicolson (formal order 2)' }, { value: 'be', label: 'Backward Euler (formal order 1)' }], help: 'Implicit time-stepping scheme whose temporal order is measured. Both are solved with the tridiagonal kernel solver.' },
  { key: 'alpha', label: 'Diffusivity α', unit: 'm²/s', value: 1e-6, min: 1e-12, max: 1e-2, log: true, bus: 'prop.D_diff', group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Thermal or mass diffusivity of the benchmark medium. Only the Fourier number αt/L² controls the numerical behaviour.' },
  { key: 'Lx', label: 'Domain length L', unit: 'm', value: 0.01, min: 1e-9, max: 10, log: true, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Length of the 1-D benchmark domain.' },
  { key: 'Fo', label: 'Final Fourier number αt/L²', value: 0.1, min: 0.005, max: 0.5, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Dimensionless end time of the benchmark. Larger values need more time steps to stay in the asymptotic range.' },
  { key: 'nx0', label: 'Coarsest grid intervals', value: 10, min: 6, max: 40, step: 2, sens: false, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Number of intervals of the coarsest spatial grid; the study uses this grid and two successive doublings.' },
  { key: 'nt0', label: 'Coarsest number of time steps', value: 16, min: 4, max: 64, step: 2, sens: false, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Number of time steps of the coarsest temporal resolution; the study uses this and two successive doublings.' },
  { key: 'odeN', label: 'ODE test: coarsest step count', value: 8, min: 4, max: 64, step: 2, sens: false, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'Number of steps used for the logistic-equation order check of the RK4 and implicit-Euler integrators (also doubled twice; implicit Euler uses four times as many steps).' },
  { key: 'ordTol', label: 'Acceptance: |observed − formal order|', value: 0.25, min: 0.02, max: 1, group: 'Code verification', showIf: { mode: ['mms'] }, help: 'A verification test passes when the observed order of accuracy lies within this distance of the formal order.' },
  // UQ
  { key: 'umodel', label: 'Model propagated', type: 'select', value: 'arrh', group: 'Uncertainty', showIf: { mode: ['uq'] }, options: Object.entries(UQM).map(([value, m]) => ({ value, label: m.label })), help: 'Built-in reduced-order model whose output uncertainty is computed from uncertain inputs. The input table of the result lists the nominal values and coefficients of variation in use.' },
  { key: 'uqOverride', label: 'Input overrides', type: 'text', value: '', group: 'Uncertainty', showIf: { mode: ['uq'] }, help: 'Optional overrides "name=mean,cv%" separated by semicolons, e.g. "T=368,0.3; Ea=55,4". Names are those listed in the input table of the result; unlisted inputs keep their defaults.' },
  { key: 'cvScale', label: 'Uncertainty multiplier', value: 1, min: 0, max: 4, step: 0.1, group: 'Uncertainty', showIf: { mode: ['uq'] }, help: 'Multiplies every input coefficient of variation; 0 removes all uncertainty and 2 doubles it. Useful for asking how much tighter the inputs must be controlled.' },
  { key: 'dist', label: 'Input distribution', type: 'select', value: 'normal', group: 'Uncertainty', showIf: { mode: ['uq'] }, options: [{ value: 'normal', label: 'Normal' }, { value: 'lognormal', label: 'Lognormal (same mean and CV)' }, { value: 'uniform', label: 'Uniform (same mean and standard deviation)' }], help: 'Marginal distribution applied to all uncertain inputs. Lognormal is appropriate for strictly positive, multiplicative quantities such as rate constants.' },
  { key: 'rho12', label: 'Correlation of inputs 1 and 2', value: 0, min: -0.95, max: 0.95, step: 0.05, group: 'Uncertainty', showIf: { mode: ['uq'] }, help: 'Gaussian-copula correlation between the first two inputs (for example the compensating A and Eₐ of an Arrhenius fit). Sobol indices are always evaluated for independent inputs.' },
  { key: 'sampler', label: 'Sampling', type: 'select', value: 'lhs', group: 'Uncertainty', showIf: { mode: ['uq'] }, options: [{ value: 'lhs', label: 'Latin hypercube' }, { value: 'mc', label: 'Simple Monte Carlo' }], help: 'Latin hypercube stratifies each input and usually converges faster than simple random sampling for the mean.' },
  { key: 'nSamp', label: 'Samples', value: 2000, min: 100, max: 20000, step: 100, sens: false, group: 'Uncertainty', showIf: { mode: ['uq'] }, help: 'Number of model evaluations for propagation. The standard error of the mean falls as 1/√n.' },
  { key: 'limRatio', label: 'Limit as multiple of nominal output', value: 1.2, min: 0, max: 5, step: 0.05, group: 'Uncertainty', showIf: { mode: ['uq'] }, help: 'Specification limit expressed relative to the nominal (mean-input) output; 0 gives the limit "zero", useful for the probability of a negative NPV.' },
  { key: 'limSide', label: 'Failure side', type: 'select', value: 'above', group: 'Uncertainty', showIf: { mode: ['uq'] }, options: [{ value: 'above', label: 'Output above the limit' }, { value: 'below', label: 'Output below the limit' }], help: 'Whether exceeding or falling short of the limit counts as failure for the probability estimate.' },
  // credibility
  { key: 'g1', label: 'Gate 1 evidence: input data records', type: 'select', value: 'complete', group: 'Credibility', showIf: { mode: ['cred'] }, options: [{ value: 'complete', label: 'Units, geometry, material records, traceability and uncertainty documented' }, { value: 'partial', label: 'Partly documented' }, { value: 'none', label: 'Not documented' }], help: 'Quality of the input data: unit consistency, geometry integrity, material records, experimental traceability and stated uncertainties.' },
  { key: 'g2', label: 'Gate 2: model consistency checks passed', type: 'bool', value: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Governing equations, constitutive models, sign conventions and physical admissibility have been reviewed and the automatic checks pass.' },
  { key: 'consRes', label: 'Conservation residual (relative)', value: 1e-8, min: 1e-16, max: 1, log: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest relative imbalance of mass, charge or energy reported by the solver for the calculation being qualified.' },
  { key: 'consTol', label: 'Acceptance: conservation residual', value: 1e-4, min: 1e-12, max: 1, log: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest conservation imbalance acceptable for the intended use.' },
  { key: 'g3', label: 'Gate 3: benchmark/MMS tests executed', type: 'bool', value: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Analytical-benchmark, manufactured-solution and regression tests exist for the solver version in use.' },
  { key: 'ordDev', label: 'Observed-order deviation', value: 0.05, min: 0, max: 3, group: 'Credibility', bus: 'cred.order_dev', showIf: { mode: ['cred'] }, help: 'Largest |observed − formal| order of accuracy found in the code-verification tests.' },
  { key: 'ordTolC', label: 'Acceptance: order deviation', value: 0.25, min: 0.02, max: 1, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest acceptable deviation between observed and formal order of accuracy.' },
  { key: 'gciPct', label: 'Fine-grid GCI', unit: '%', value: 0.5, min: 0, max: 100, group: 'Credibility', bus: 'cred.gci_pct', showIf: { mode: ['cred'] }, help: 'Numerical uncertainty of the quantity of interest from a three-grid study (or time-step study, whichever is larger).' },
  { key: 'gciTolC', label: 'Acceptance: GCI', unit: '%', value: 2, min: 0.001, max: 50, log: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest numerical uncertainty acceptable for the intended use.' },
  { key: 'g4conv', label: 'Gate 4: iterative convergence demonstrated', type: 'bool', value: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Residual and solution-update criteria were met and false convergence (stagnation, loose tolerances) was excluded.' },
  { key: 'g5', label: 'Gate 5: parameters calibrated on identified data', type: 'bool', value: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Uncertain parameters were estimated from identified, archived experimental or reference datasets with a reproducible fitting record.' },
  { key: 'relSe', label: 'Largest parameter relative standard error', unit: '%', value: 12, min: 0, max: 1000, group: 'Credibility', bus: 'cred.max_rel_se_pct', showIf: { mode: ['cred'] }, help: 'Worst relative standard error among calibrated parameters; large values signal weak practical identifiability.' },
  { key: 'relSeTol', label: 'Acceptance: parameter relative error', unit: '%', value: 30, min: 1, max: 500, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest acceptable relative standard error of a calibrated parameter.' },
  { key: 'nVal', label: 'Independent validation datasets', value: 0, min: 0, max: 50, step: 1, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Number of datasets, not used in calibration, against which predictions have been compared.' },
  { key: 'indep', label: 'Validation data independent of calibration', type: 'bool', value: true, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Confirms that none of the validation measurements were used for parameter fitting or model selection.' },
  { key: 'valRmse', label: 'Validation normalised RMSE', unit: '%', value: 8, min: 0, max: 1000, group: 'Credibility', bus: 'cred.rmse_pct', showIf: { mode: ['cred'] }, help: 'Normalised root-mean-square prediction error on the independent validation data.' },
  { key: 'accRmseC', label: 'Acceptance: validation RMSE', unit: '%', value: 10, min: 0.1, max: 100, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest validation error acceptable for the intended engineering use.' },
  { key: 'g7uq', label: 'Gate 7: uncertainty propagated to predictions', type: 'bool', value: false, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Input, parameter and numerical uncertainties have been propagated to prediction intervals.' },
  { key: 'g7sens', label: 'Gate 7: global sensitivity analysis done', type: 'bool', value: false, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'A global sensitivity ranking (Sobol, Morris or equivalent) identifies the inputs that control the prediction.' },
  { key: 'extrap', label: 'Extrapolation beyond validated domain', unit: '%', value: 0, min: 0, max: 500, group: 'Credibility', bus: 'cred.extrap_pct', showIf: { mode: ['cred'] }, help: 'Distance of the application condition outside the validated range, as a percentage of that range (0 = interpolation).' },
  { key: 'extrapMax', label: 'Acceptance: extrapolation', unit: '%', value: 10, min: 0, max: 200, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'Largest extrapolation the application owner accepts with documented justification.' },
  { key: 'g8use', label: 'Gate 8: intended use and limits documented', type: 'bool', value: false, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'A credibility record states the intended engineering task, consequence of error, applicability limits and known deficiencies.' },
  { key: 'g8rev', label: 'Gate 8: independent review/approval obtained', type: 'bool', value: false, group: 'Credibility', showIf: { mode: ['cred'] }, help: 'The qualification record has been reviewed and approved by the responsible authority where this is required.' },
];

// ---------- mode: Richardson / GCI ----------
function runGci(inp) {
  const w = [], pts = [[inp.gh1, inp.gf1], [inp.gh2, inp.gf2], [inp.gh3, inp.gf3]].sort((a, b) => a[0] - b[0]);
  const [h1, h2, h3] = pts.map(p => p[0]), [f1, f2, f3] = pts.map(p => p[1]);
  if (!(inp.gh1 < inp.gh2 && inp.gh2 < inp.gh3)) w.push('Grid sizes were not in the order fine < medium < coarse; the three solutions were re-ordered by grid size.');
  let r21 = h2 / h1, r32 = h3 / h2;
  if (r21 < 1.1 || r32 < 1.1) { w.push('Refinement ratio below 1.1: differences between solutions are dominated by round-off and iteration error; ratios were floored at 1.1 for the estimate.'); r21 = Math.max(r21, 1.1); r32 = Math.max(r32, 1.1); }
  const e21 = f2 - f1, e32 = f3 - f2, scl = Math.max(Math.abs(f1), Math.abs(f2), Math.abs(f3), 1e-300), degenerate = Math.abs(e21) < 1e-13 * scl || Math.abs(e32) < 1e-13 * scl;
  const osc = !degenerate && e32 / e21 < 0, s = osc ? -1 : 1; let p = inp.gOrder, it = 0;
  if (!degenerate) { const lr = Math.log(Math.abs(e32 / e21)); p = N.clamp(Math.abs(lr) / Math.log(r21), 0.05, 20); for (; it < 100; it++) { const q = Math.log((Math.pow(r21, p) - s) / (Math.pow(r32, p) - s)), pn = N.clamp(Math.abs(lr + q) / Math.log(r21), 0.05, 20); if (Math.abs(pn - p) < 1e-12) { p = pn; break; } p = pn; } }
  else w.push('At least two of the three solutions are identical to the precision supplied: the observed order cannot be estimated and the formal order is used.');
  const pc = N.clamp(p, 0.5, inp.gOrder), rp = Math.pow(r21, p), fext = (rp * f1 - f2) / (rp - 1), ea21 = Math.abs(e21) / Math.max(Math.abs(f1), 1e-300), eext = Math.abs(fext - f1) / Math.max(Math.abs(fext), 1e-300);
  const g21 = inp.gFs * ea21 / (rp - 1), g32 = inp.gFs * Math.abs(e32) / Math.max(Math.abs(f2), 1e-300) / (Math.pow(r32, p) - 1), gcons = inp.gFs * ea21 / (Math.pow(r21, pc) - 1);
  const asy = g21 > 0 ? g32 / (rp * g21) : 1, inAsy = Math.abs(asy - 1) < 0.1 && !osc, ordOk = Math.abs(p - inp.gOrder) <= 0.5;
  if (osc) w.push('Oscillatory convergence (successive differences change sign): Richardson extrapolation is unreliable; add a fourth grid or refine further.');
  if (!degenerate && !ordOk) w.push(`Observed order ${F(p, 3)} differs from the formal order ${inp.gOrder} by more than 0.5: the grids are probably outside the asymptotic range or another error source (iteration, boundary treatment, singularity) dominates.`);
  const C = (f3 - fext) / Math.pow(h3, p), hs = N.linspace(0, h3 * 1.05, 60), ok = g21 * 100 <= inp.gTol;
  // kernel cross-check for uniform refinement
  let kern = null; if (!degenerate && Math.abs(r21 - r32) < 1e-9 * r21) kern = N.gci(f1, f2, f3, r21, inp.gFs);
  return {
    summary: `Three-grid Richardson analysis gives an observed order of accuracy p = ${F(p, 4)} (formal ${inp.gOrder}) and an extrapolated value ${F(fext, 6)}. The fine-grid GCI is ${pct(g21)} (medium grid ${pct(g32)}); the asymptotic-range ratio is ${F(asy, 4)}${inAsy ? ', so the solutions are in the asymptotic range' : ', so asymptotic convergence is not demonstrated'}. Against the stated tolerance of ${inp.gTol} % the fine-grid numerical uncertainty is ${ok ? 'acceptable' : 'not acceptable'}.`,
    kpis: [
      kp('order', 'Observed order of accuracy p', p, '–', `Apparent convergence rate from the three solutions (Celik fixed-point procedure for refinement ratios ${F(r21, 3)} and ${F(r32, 3)}); it should approach the formal order ${inp.gOrder} of the scheme.`, ordOk && !osc),
      kp('f_ext', 'Richardson-extrapolated value', fext, '', 'Estimate of the zero-grid-size solution f_ext = (r²¹ᵖ·f₁ − f₂)/(r²¹ᵖ − 1); it is an estimate, not an exact value, and is only meaningful in the asymptotic range.', null),
      kp('gci_fine', 'Fine-grid GCI', 100 * g21, '%', `Grid convergence index Fₛ·|ε₂₁|/(rᵖ − 1) with Fₛ = ${inp.gFs}: a banded estimate of the discretisation uncertainty of the fine-grid result.`, ok),
      kp('gci_med', 'Medium-grid GCI', 100 * g32, '%', 'Grid convergence index of the medium-grid solution; shows whether the cheaper grid would already meet the tolerance.', 100 * g32 <= inp.gTol),
      kp('gci_cons', 'Conservative fine-grid GCI', 100 * gcons, '%', `GCI recomputed with the order capped at the formal order (p = ${F(pc, 3)}), which avoids optimistic error bands when the observed order is super-convergent.`, 100 * gcons <= inp.gTol),
      kp('asy', 'Asymptotic-range ratio', asy, '–', 'GCI₃₂/(rᵖ·GCI₂₁); a value within about 10 % of 1 indicates that the three grids lie in the asymptotic range of convergence.', inAsy),
      kp('ea21', 'Relative change fine–medium', 100 * ea21, '%', 'Approximate relative error |f₁ − f₂|/|f₁| between the two finest solutions, before applying the convergence-rate correction.', null),
      kp('eext', 'Extrapolated relative error', 100 * eext, '%', 'Estimated relative discretisation error of the fine-grid result with respect to the Richardson-extrapolated value.', null),
    ],
    series: [
      { key: 'conv', title: 'Quantity of interest versus grid size', xlabel: 'grid size h', ylabel: 'f', lines: [{ name: 'computed solutions', x: [h1, h2, h3], y: [f1, f2, f3] }, { name: 'f_ext + C·hᵖ', x: hs, y: hs.map(h => fext + C * Math.pow(h, p)) }], explain: 'The three computed solutions and the power-law error model fitted through them. The intercept at h = 0 is the Richardson-extrapolated value; a smooth monotone approach indicates asymptotic convergence.' },
      { key: 'err', title: 'Estimated discretisation error (log–log)', xlabel: 'grid size h', ylabel: '|f − f_ext|', logx: true, logy: true, lines: [{ name: 'estimated error', x: [h1, h2, h3], y: [f1, f2, f3].map(v => Math.max(Math.abs(v - fext), 1e-300)) }], explain: 'Estimated error of each solution against grid size on logarithmic axes. In the asymptotic range the points fall on a straight line whose slope is the observed order p.' },
    ],
    tables: [{ key: 'grids', title: 'Grid-convergence study', columns: ['Grid', 'h', 'f', 'Refinement ratio to next', 'GCI (%)'], rows: [['fine', h1, f1, +r21.toPrecision(5), +(100 * g21).toPrecision(4)], ['medium', h2, f2, +r32.toPrecision(5), +(100 * g32).toPrecision(4)], ['coarse', h3, f3, '–', '–']], explain: 'Input solutions ordered from fine to coarse with the refinement ratios used and the grid convergence index attached to the fine and medium solutions.' }],
    equations: [
      { kind: 'formula', label: 'Observed order (non-uniform refinement)', text: 'p = |ln|ε₃₂/ε₂₁| + q(p)| / ln r₂₁,   q(p) = ln[(r₂₁ᵖ − s)/(r₃₂ᵖ − s)],   s = sign(ε₃₂/ε₂₁)' },
      { kind: 'formula', label: 'Richardson extrapolation', text: 'f_ext = (r₂₁ᵖ·f₁ − f₂)/(r₂₁ᵖ − 1)' },
      { kind: 'formula', label: 'Grid convergence index', text: 'GCI₂₁ = Fₛ·|(f₁ − f₂)/f₁| / (r₂₁ᵖ − 1);   asymptotic range if GCI₃₂/(r₂₁ᵖ·GCI₂₁) ≈ 1' },
      { kind: 'governing', label: 'Truncation-error model', text: 'f(h) = f_exact + C·hᵖ + O(hᵖ⁺¹)' },
    ],
    checks: [
      { name: 'Monotone convergence', pass: !osc && !degenerate, detail: degenerate ? 'solutions identical: no information' : `ε₃₂/ε₂₁ = ${F(e32 / e21, 4)}` },
      { name: 'Asymptotic range', pass: inAsy, detail: `GCI₃₂/(rᵖ·GCI₂₁) = ${F(asy, 4)} (target 1 ± 0.1)` },
      { name: 'Observed order consistent with formal order', pass: ordOk, detail: `p = ${F(p, 4)} vs formal ${inp.gOrder}` },
      { name: 'Fine-grid GCI within acceptance tolerance', pass: ok, detail: `${pct(g21)} vs ${inp.gTol} %` },
      ...(kern ? [{ name: 'Kernel cross-check (uniform refinement)', pass: Math.abs(kern.order - p) < 1e-6 * Math.max(1, p), detail: `kernel gci(): p = ${F(kern.order, 6)}, f_ext = ${F(kern.extrapolated, 6)}` }] : []),
    ],
    warnings: w,
    recommendations: [
      { priority: ok ? 'low' : 'high', text: ok ? `Report the fine-grid result as ${F(f1, 6)} ± ${pct(g21)} (numerical uncertainty) and keep the grid triple in the verification record.` : `Refine further: to reach ${inp.gTol} % the grid size must fall to about ${F(h1 * Math.pow(inp.gTol / Math.max(100 * g21, 1e-300), 1 / p), 3)} (h₁ scaled by (tolerance/GCI)^(1/p)).`, basis: 'GCI scales as hᵖ in the asymptotic range.' },
      { priority: inAsy ? 'low' : 'high', text: inAsy ? 'The asymptotic check is satisfied; the same refinement level can be reused for neighbouring design points of the same problem class.' : 'Add a fourth, finer grid and repeat: the asymptotic-range ratio is not close to 1, so the reported GCI is not yet a defensible error band.', basis: 'Roache asymptotic-range criterion.' },
      { priority: 'medium', text: 'Repeat the study for the time step and for the iterative tolerance; the total numerical uncertainty is the combination of spatial, temporal and iterative contributions.', basis: 'Solution verification must cover discretisation, iterative and coupling errors.' },
      { priority: 'medium', text: 100 * g32 <= inp.gTol ? `The medium grid already meets the tolerance (GCI ${pct(g32)}): use it for parameter sweeps and optimisation to cut computing time and energy by roughly the cell-count ratio.` : 'Use the coarsest grid that meets the tolerance for design sweeps; over-resolved runs spend computing energy without reducing the dominant uncertainty.', basis: 'Computational cost grows as h⁻ᵈ while the error falls only as hᵖ.', sustainability: true },
    ],
    publish: { 'cred.gci': g21, 'cred.gci_pct': 100 * g21, 'cred.order': p, 'cred.order_dev': Math.abs(p - inp.gOrder) },
    numerics: { method: 'Richardson extrapolation with fixed-point iteration for the observed order (non-uniform refinement ratios)', dof: 3, steps: 1, iterations: it + 1, residual: Math.abs(asy - 1), notes: 'Requires monotone convergence and solutions in the asymptotic range; iterative error should be at least two orders below the discretisation error.' },
  };
}

// ---------- mode: code verification by exact / manufactured solution ----------
function heatSolve(nx, nt, th, a, L, T, mms) {
  const h = L / nx, dt = T / nt, r = a * dt / (h * h), n = nx - 1, k = (mms ? 2 : 1) * Math.PI / L, w = 2 * Math.PI / T;
  const amp = t => mms ? 1 + 0.5 * Math.sin(w * t) : Math.exp(-a * k * k * t), src = t => mms ? 0.5 * w * Math.cos(w * t) + a * k * k * (1 + 0.5 * Math.sin(w * t)) : 0;
  const sx = Float64Array.from({ length: n }, (_, i) => Math.sin(k * (i + 1) * h)); let u = sx.map(v => v * amp(0));
  const lo = new Float64Array(n).fill(-th * r), di = new Float64Array(n).fill(1 + 2 * th * r), up = new Float64Array(n).fill(-th * r), rhs = new Float64Array(n);
  for (let s = 0; s < nt; s++) {
    const t0 = s * dt, q = dt * (th * src(t0 + dt) + (1 - th) * src(t0));
    for (let i = 0; i < n; i++) rhs[i] = u[i] + (1 - th) * r * ((i > 0 ? u[i - 1] : 0) - 2 * u[i] + (i < n - 1 ? u[i + 1] : 0)) + q * sx[i];
    u = N.solveTridiag(lo, di, up, rhs);
  }
  const A = amp(T); let e2 = 0, em = 0, q2 = 0; for (let i = 0; i < n; i++) { const d = u[i] - A * sx[i]; e2 += d * d; em = Math.max(em, Math.abs(d)); q2 += u[i] * u[i] * h; }
  return { u, l2: Math.sqrt(e2 / n), max: em, qoi: q2, exactQoi: A * A * L / 2, amp: A };
}
const ordOf = (e1, e2, r = 2) => e1 > 0 && e2 > 0 ? Math.log(e1 / e2) / Math.log(r) : 0;
function runMms(inp) {
  const w = [], mms = inp.vcase === 'mms', th = inp.scheme === 'be' ? 1 : 0.5, pT = th === 1 ? 1 : 2, a = inp.alpha, L = inp.Lx, T = inp.Fo * L * L / a;
  const nx0 = Math.max(4, Math.round(inp.nx0)), nt0 = Math.max(2, Math.round(inp.nt0)), nxs = [nx0, 2 * nx0, 4 * nx0], nts = [nt0, 2 * nt0, 4 * nt0];
  // spatial study: Crank–Nicolson with a time step small enough that temporal error is negligible
  const sp = nxs.map(nx => heatSolve(nx, 4000, 0.5, a, L, T, mms)), es = sp.map(s => s.l2 / Math.max(Math.abs(s.amp), 1e-300)), pS = ordOf(es[1], es[2]), pS0 = ordOf(es[0], es[1]);
  // temporal study on a fixed grid against a reference with 64× more Crank–Nicolson steps (spatial error cancels)
  const nxt = 2 * nx0, ref = heatSolve(nxt, 64 * nt0, 0.5, a, L, T, mms), tm = nts.map(nt => heatSolve(nxt, nt, th, a, L, T, mms));
  const et = tm.map(s => { let m = 0; for (let i = 0; i < s.u.length; i++) m = Math.max(m, Math.abs(s.u[i] - ref.u[i])); return m / Math.max(Math.abs(s.amp), 1e-300); }), pTo = ordOf(et[1], et[2]), pT0 = ordOf(et[0], et[1]);
  // GCI of the integral quantity ∫u² dx on the three grids (kernel routine)
  const q = sp.map(s => s.qoi), dq1 = q[1] - q[2], dq2 = q[0] - q[1]; let g = { order: 2, extrapolated: q[2], gciFine: 0 };
  if (Math.abs(dq1) > 1e-15 * Math.abs(q[2]) && Math.abs(dq2) > 1e-15 * Math.abs(q[2]) && q[2] !== 0) g = N.gci(q[2], q[1], q[0], 2, 1.25);
  const gciF = Number.isFinite(g.gciFine) ? Math.abs(g.gciFine) : 0, qErr = Math.abs(q[2] - sp[2].exactQoi) / Math.max(sp[2].exactQoi, 1e-300);
  // ODE integrators: logistic equation y' = y(1 − y), y(0) = 0.1, t ∈ [0, 4]
  const y0 = 0.1, yT = 1 / (1 + (1 / y0 - 1) * Math.exp(-4)), fo = (t, y) => [y[0] * (1 - y[0])], n0 = Math.max(2, Math.round(inp.odeN)), ns = [n0, 2 * n0, 4 * n0];
  const eR = ns.map(n => Math.abs(N.rk4(fo, [y0], 0, 4, n).y.at(-1)[0] - yT)), eB = ns.map(n => Math.abs(N.bdf1(fo, [y0], 0, 4, 4 * n, { tol: 1e-13 }).y.at(-1)[0] - yT)), pR = ordOf(eR[1], eR[2]), pB = ordOf(eB[1], eB[2]);
  const tol = inp.ordTol, okS = Math.abs(pS - 2) <= tol, okT = Math.abs(pTo - pT) <= tol, okR = Math.abs(pR - 4) <= Math.max(tol, 0.3), okB = Math.abs(pB - 1) <= tol, all = okS && okT && okR && okB;
  if (!okT) w.push(`Temporal order ${F(pTo, 3)} deviates from the formal order ${pT}: with ${nt0} coarsest steps the product of time step and decay rate is too large for the asymptotic range; increase the number of time steps.`);
  if (!okS) w.push(`Spatial order ${F(pS, 3)} deviates from 2: increase the coarsest grid resolution.`);
  if (eR[2] < 1e-12) w.push('RK4 error on the finest ODE step is close to round-off; the observed order of RK4 is then contaminated by floating-point noise.');
  const hs = nxs.map(nx => L / nx), dts = nts.map(nt => T / nt), pos1 = v => Math.max(v, 1e-300), dev = Math.max(Math.abs(pS - 2), Math.abs(pTo - pT), Math.abs(pR - 4), Math.abs(pB - 1));
  return {
    summary: `Code verification on the transient heat equation (${mms ? 'manufactured solution with source' : 'exact decaying mode'}, Fo = ${F(inp.Fo, 3)}): the observed spatial order is ${F(pS, 4)} (formal 2) and the observed temporal order of ${th === 1 ? 'backward Euler' : 'Crank–Nicolson'} is ${F(pTo, 4)} (formal ${pT}). The ODE integrators give orders ${F(pR, 4)} (RK4, formal 4) and ${F(pB, 4)} (implicit Euler, formal 1). ${all ? 'All tests lie within' : 'Not all tests lie within'} the acceptance band of ±${tol}; this verifies the implementation for this problem class only, not the physics of any application.`,
    kpis: [
      kp('p_space', 'Observed spatial order', pS, '–', 'Slope of the relative L2 error against grid size between the two finest grids (second-order central differences should give 2). Errors are measured against the analytical solution.', okS),
      kp('p_time', 'Observed temporal order', pTo, '–', `Slope of the error against time step between the two smallest steps for the ${th === 1 ? 'backward-Euler' : 'Crank–Nicolson'} scheme, measured against a reference on the same grid so that the spatial error cancels.`, okT),
      kp('p_rk4', 'Observed order of RK4', pR, '–', 'Convergence rate of the explicit fourth-order Runge–Kutta kernel on the logistic equation with known analytical solution; the formal value is 4.', okR),
      kp('p_bdf1', 'Observed order of implicit Euler', pB, '–', 'Convergence rate of the stiff implicit (backward-Euler/BDF1 with Newton iteration) kernel on the logistic equation; the formal value is 1.', okB),
      kp('err_fine', 'Relative L2 error on finest grid', es[2], '–', `Discretisation error of the ${nxs[2]}-interval solution relative to the exact amplitude; it is a true error because the exact solution is known.`, null),
      kp('err_dt', 'Temporal error at smallest step', et[2], '–', `Maximum-norm difference to the temporal reference with ${nts[2]} time steps, relative to the solution amplitude.`, null),
      kp('gci', 'GCI of ∫u² dx (finest grid)', 100 * gciF, '%', 'Grid convergence index of the integral quantity computed from the three spatial grids with the kernel Richardson routine; it should bound the true error of that quantity.', gciF >= 0.8 * qErr),
      kp('qoi_err', 'True error of ∫u² dx', 100 * qErr, '%', 'Actual relative error of the integral quantity on the finest grid against its analytical value, used to test whether the GCI is a conservative estimate.', null),
    ],
    series: [
      { key: 'sp', title: 'Spatial convergence (log–log)', xlabel: 'grid size h (m)', ylabel: 'relative L2 error', logx: true, logy: true, lines: [{ name: 'computed error', x: hs, y: es.map(pos1) }, { name: 'slope 2 reference', x: hs, y: hs.map(h => pos1(es[2] * (h / hs[2]) ** 2)) }], explain: 'True discretisation error against grid size for three grids. The computed points should run parallel to the second-order reference line; a shallower slope reveals an implementation or boundary-treatment defect.' },
      { key: 'tm', title: 'Temporal convergence (log–log)', xlabel: 'time step Δt (s)', ylabel: 'relative max error', logx: true, logy: true, lines: [{ name: 'computed error', x: dts, y: et.map(pos1) }, { name: `slope ${pT} reference`, x: dts, y: dts.map(d => pos1(et[2] * (d / dts[2]) ** pT)) }], explain: 'Time-integration error against step size on a fixed grid. The slope equals the order of the time scheme when the steps are in the asymptotic range.' },
      { key: 'ode', title: 'ODE integrator convergence (log–log)', xlabel: 'step size h', ylabel: 'absolute error at t = 4', logx: true, logy: true, lines: [{ name: 'RK4', x: ns.map(n => 4 / n), y: eR.map(pos1) }, { name: 'implicit Euler (BDF1)', x: ns.map(n => 1 / n), y: eB.map(pos1) }], explain: 'Global error of the two kernel ODE integrators on the logistic equation. RK4 should fall sixteen-fold and implicit Euler two-fold per halving of the step.' },
    ],
    tables: [{ key: 'conv', title: 'Convergence record', columns: ['Study', 'Level 1', 'Level 2', 'Level 3', 'Error 1', 'Error 2', 'Error 3', 'Order 1→2', 'Order 2→3', 'Formal'], rows: [
      ['Space (intervals)', ...nxs, ...es.map(v => F(v, 4)), F(pS0, 4), F(pS, 4), 2], [`Time (${th === 1 ? 'BE' : 'CN'} steps)`, ...nts, ...et.map(v => F(v, 4)), F(pT0, 4), F(pTo, 4), pT],
      ['RK4 (steps)', ...ns, ...eR.map(v => F(v, 4)), F(ordOf(eR[0], eR[1]), 4), F(pR, 4), 4], ['BDF1 (steps)', ...ns.map(n => 4 * n), ...eB.map(v => F(v, 4)), F(ordOf(eB[0], eB[1]), 4), F(pB, 4), 1]], explain: 'Errors and observed orders at three refinement levels for each verification study. Agreement of both order columns shows that the finest levels are in the asymptotic range.' }],
    equations: [
      { kind: 'governing', label: 'Benchmark equation', text: mms ? '∂u/∂t = α·∂²u/∂x² + S(x,t),  u(0,t) = u(L,t) = 0' : '∂u/∂t = α·∂²u/∂x²,  u(0,t) = u(L,t) = 0' },
      { kind: 'formula', label: mms ? 'Manufactured solution and source' : 'Exact solution', text: mms ? 'u_m = [1 + ½·sin(2πt/T)]·sin(2πx/L);   S = ∂u_m/∂t − α·∂²u_m/∂x²' : 'u = exp(−α·π²·t/L²)·sin(πx/L)' },
      { kind: 'formula', label: 'θ-scheme (θ = 1 backward Euler, θ = ½ Crank–Nicolson)', text: '(uⁿ⁺¹ − uⁿ)/Δt = α·[θ·δₓ²uⁿ⁺¹ + (1 − θ)·δₓ²uⁿ] + θ·Sⁿ⁺¹ + (1 − θ)·Sⁿ' },
      { kind: 'formula', label: 'Observed order of accuracy', text: 'p = ln(e_coarse/e_fine)/ln(r),  r = 2' },
      { kind: 'governing', label: 'ODE benchmark', text: 'dy/dt = y·(1 − y),  y(0) = 0.1  ⇒  y = 1/(1 + 9·e^(−t))' },
    ],
    checks: [
      { name: 'Spatial order within tolerance', pass: okS, detail: `observed ${F(pS, 4)} vs formal 2 (±${tol})` },
      { name: 'Temporal order within tolerance', pass: okT, detail: `observed ${F(pTo, 4)} vs formal ${pT} (±${tol})` },
      { name: 'RK4 order within tolerance', pass: okR, detail: `observed ${F(pR, 4)} vs formal 4` },
      { name: 'Implicit Euler order within tolerance', pass: okB, detail: `observed ${F(pB, 4)} vs formal 1` },
      { name: 'Asymptotic range (consecutive orders agree)', pass: Math.abs(pS - pS0) < 0.3 && Math.abs(pTo - pT0) < 0.5, detail: `space ${F(pS0, 3)} → ${F(pS, 3)}; time ${F(pT0, 3)} → ${F(pTo, 3)}` },
      { name: 'GCI bounds the true error of the integral quantity', pass: gciF >= 0.8 * qErr, detail: `GCI ${pct(gciF)} vs true error ${pct(qErr)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: all ? 'low' : 'high', text: all ? 'Archive this run as a regression test: re-execute it whenever the solver kernels, the time integration or the boundary treatment change.' : 'Do not rely on results of the affected scheme until the failed order test is understood; first repeat with finer coarsest levels to exclude pre-asymptotic behaviour.', basis: 'Code verification by order-of-accuracy testing is the most sensitive implementation check.' },
      { priority: 'medium', text: `For transient production runs choose Crank–Nicolson when smooth solutions are expected (error ${F(et[2], 2)} at ${nts[2]} steps here) and backward Euler for stiff start-up transients, where its damping avoids oscillation.`, basis: 'Measured temporal error levels and stability properties of the θ-scheme.' },
      { priority: 'medium', text: 'Extend the manufactured-solution test to each additional physics term (reaction, variable coefficients, flux boundaries) before trusting coupled multiphysics results; verify each coupling separately.', basis: 'Gate 3 of the qualification sequence requires tests per solver and per coupling.' },
      { priority: 'low', text: `Second-order convergence means halving h quarters the error at about ${th === 1 ? 'four' : 'four'} times the work in 1-D transient problems; select the resolution from the required error rather than by habit to avoid wasted computing energy.`, basis: 'Observed convergence rates and cost scaling.', sustainability: true },
    ],
    publish: { 'cred.order': pS, 'cred.order_dev': dev, 'cred.gci': gciF, 'cred.gci_pct': 100 * gciF, 'cred.code_verified': all ? 1 : 0 },
    numerics: { method: 'Second-order central finite differences + θ-scheme (tridiagonal Thomas solver); RK4 and implicit Euler with Newton iteration for the ODE test', dof: nxs[2] - 1, steps: 4000, iterations: 3 * 4000 + 7 * nt0 + 64 * nt0, residual: es[2], notes: 'Spatial study uses Crank–Nicolson with 4000 steps; temporal study is referenced to a 64× finer Crank–Nicolson solution on the same grid.' },
  };
}

// ---------- shared: read calibration table ----------
function readCal(inp, w) {
  const t = parseTable(inp.cdata, 2, 3, 'Calibration table'), X = t.rows.map(r => r[0]), Y = t.rows.map(r => r[1]);
  const hasS = t.rows.every(r => Number.isFinite(r[2]) && r[2] > 0), S = hasS ? t.rows.map(r => r[2]) : X.map(() => 1);
  if (!hasS) w.push('No positive sigma column in the calibration table: the fit is unweighted and the noise level is estimated from the residuals.');
  if (inp.cdata === CAL_DEF) w.push('Calibration table: synthetic demonstration data — replace with measurements.');
  return { X, Y, S, hasS, head: t.head };
}

// ---------- mode: calibration ----------
function runCalib(inp) {
  const w = [], fam = FAM[inp.cmodel] || FAM.arrhenius, { X, Y, S, hasS, head } = readCal(inp, w), m = X.length;
  const fit = calibrate(X, Y, S, fam, inp.scaleCov || !hasS), n = fit.n, p = fit.p, tq = t975(fit.dof), red = fit.chi2 / fit.dof;
  if (m <= n) w.push(`Only ${m} observations for ${n} parameters: no degrees of freedom remain, so uncertainties are not meaningful.`);
  const rel = p.map((v, j) => fit.se[j] / Math.max(Math.abs(v), 1e-300)), worst = N.argmax(rel), rm = N.rmse(fit.pred, Y), r2 = N.r2(fit.pred, Y), ym = Math.max(N.mean(Y.map(Math.abs)), 1e-300);
  if (!fit.identifiable) w.push('The Fisher information matrix is numerically singular: at least one parameter combination cannot be determined from these data.');
  const maxCorr = n > 1 ? N.maxOf(fit.corr.flatMap((r, a) => r.filter((_, b) => b > a).map(Math.abs))) : 0;
  if (maxCorr > 0.95) w.push(`Parameters are strongly correlated (|ρ| up to ${F(maxCorr, 4)}): individual values are poorly determined although the fitted curve is well defined.`);
  if (hasS && !inp.scaleCov && (red > 3 || red < 0.2)) w.push(`Reduced χ² = ${F(red, 3)} is far from 1: the stated measurement uncertainties are inconsistent with the scatter about this model.`);
  // profile likelihood for the least-determined parameter
  const span = Math.max(4 * fit.se[worst], 1e-3 * Math.max(Math.abs(p[worst]), fit.sc[worst])), cap = Math.max(Math.abs(p[worst]), fit.sc[worst]);
  const half = Math.min(span, 3 * cap), pg = N.linspace(Math.max(p[worst] - half, fam.b[worst][0]), Math.min(p[worst] + half, fam.b[worst][1]), 17);
  const prof = profile(X, Y, S, fam, fit, worst, pg).map(c => Math.max(0, (c - fit.chi2) / Math.max(fit.s2, 1e-300))), cross = side => { const idx = side < 0 ? [...pg.keys()].filter(i => pg[i] < p[worst]).reverse() : [...pg.keys()].filter(i => pg[i] > p[worst]); let prev = [p[worst], 0]; for (const i of idx) { if (prof[i] >= 3.84) return prev[0] + (pg[i] - prev[0]) * (3.84 - prev[1]) / Math.max(prof[i] - prev[1], 1e-300); prev = [pg[i], prof[i]]; } return null; };
  const pl = cross(-1), ph = cross(1), profOk = pl !== null && ph !== null;
  if (!profOk) w.push(`Profile likelihood of ${fam.names[worst]} does not rise above the 95 % threshold on ${pl === null && ph === null ? 'either side' : pl === null ? 'the lower side' : 'the upper side'} within the scanned range: the parameter is practically non-identifiable in that direction.`);
  // D-optimal next experiment (determinant lemma: gain = 1 + gᵀ·F⁻¹·g/σ²)
  const xmin = N.minOf(X), xmax = N.maxOf(X), dx = (xmax - xmin) || Math.abs(xmax) || 1, relS = N.mean(S.map((s, i) => s / Math.max(Math.abs(Y[i]), 1e-300)));
  const cand = N.linspace(xmin > 0 ? Math.max(xmin - 0.3 * dx, xmin * 0.5) : xmin - 0.3 * dx, xmax + 0.3 * dx, 41), grad = x => p.map((v, j) => { const h = 1e-5 * Math.max(Math.abs(v), 1e-3 * fit.sc[j]), a = [...p], b = [...p]; a[j] += h; b[j] -= h; const d = (fam.f(x, a) - fam.f(x, b)) / (2 * h); return Number.isFinite(d) ? d : 0; });
  const gain = cand.map(x => { const g = grad(x), sg = hasS ? Math.max(relS * Math.abs(fam.f(x, p)), 1e-300) : 1, q = N.dot(g, N.matVec(fit.cov0, g)) / (sg * sg); return Number.isFinite(q) ? 1 + Math.max(q, 0) : 1; }), xBest = cand[N.argmax(gain)];
  // fitted curve and prediction band from the covariance (delta method)
  const xs = N.linspace(xmin - 0.05 * dx, xmax + 0.05 * dx, 80), yf = xs.map(x => fam.f(x, p)), ub = xs.map(x => { const g = grad(x); return tq * Math.sqrt(Math.max(0, N.dot(g, N.matVec(fit.cov, g)))); });
  const series = [
    { key: 'fit', title: 'Calibration data and fitted model', xlabel: head[0] || 'x', ylabel: head[1] || 'y', lines: [{ name: 'observations', x: X, y: Y }, { name: 'fitted model', x: xs, y: yf }, { name: '95 % confidence band (lower)', x: xs, y: yf.map((v, i) => v - ub[i]) }, { name: '95 % confidence band (upper)', x: xs, y: yf.map((v, i) => v + ub[i]) }], explain: 'Observations with the least-squares model and the 95 % confidence band of the fitted mean (delta method from the parameter covariance). The band widens rapidly outside the range of the data.' },
    { key: 'res', title: 'Weighted residuals', xlabel: head[0] || 'x', ylabel: '(model − data)/σ', lines: [{ name: 'weighted residual', x: X, y: X.map((x, i) => (fit.pred[i] - Y[i]) / S[i]) }, { name: 'zero', x: [xmin, xmax], y: [0, 0] }], explain: 'Residuals divided by their measurement uncertainty. A random scatter about zero within about ±2 supports the model; a trend or curvature indicates model-form error that no parameter value can remove.' },
    { key: 'prof', title: `Profile likelihood of ${fam.names[worst]}`, xlabel: fam.names[worst], ylabel: 'Δχ²', lines: [{ name: 'profile Δχ²', x: pg, y: prof }, { name: '95 % threshold (3.84)', x: [pg[0], pg.at(-1)], y: [3.84, 3.84] }], explain: 'Increase of the fit cost when the least-determined parameter is fixed and all others are re-optimised. Crossing the threshold on both sides gives a finite likelihood-based interval; a flat branch means the data cannot bound the parameter.' },
    { key: 'lm', title: 'Levenberg–Marquardt convergence history', xlabel: 'accepted iteration', ylabel: 'χ²', logy: true, lines: [{ name: 'weighted sum of squares', x: fit.hist.map((_, i) => i), y: fit.hist.map(v => Math.max(v, 1e-300)) }], explain: 'Decrease of the weighted sum of squared residuals over the accepted optimiser iterations. A plateau reached well before the iteration limit indicates a converged fit.' },
  ];
  // Bayesian calibration: random-walk Metropolis with uniform priors on the parameter bounds
  let post = null;
  if (inp.bayes) {
    const rnd = N.rng(Math.round(inp.seed) + 4242), nMc = Math.round(N.clamp(inp.nMc, 500, 30000)), burn = Math.floor(nMc / 4), s2 = Math.max(fit.s2, 1e-300);
    const lp = q => { for (let j = 0; j < n; j++) if (q[j] < fam.b[j][0] || q[j] > fam.b[j][1]) return -Infinity; const r = fit.res(q); return -0.5 * N.dot(r, r) / s2; };
    const cq = fit.cov.map((r, a) => r.map((v, b) => v * 2.38 * 2.38 / n)); let Lc; try { Lc = N.cholesky(cq.map((r, a) => r.map((v, b) => a === b ? v * (1 + 1e-9) + 1e-300 : v))); if (!Lc.flat().every(Number.isFinite)) throw 0; } catch { Lc = N.zeros(n, n).map((r, a) => r.map((_, b) => a === b ? 2.38 / Math.sqrt(n) * Math.max(fit.se[a], 1e-3 * fit.sc[a]) : 0)); }
    let cur = [...p], lc = lp(cur), acc = 0; const ch = [];
    for (let i = 0; i < nMc; i++) { const z = cur.map(() => rnd.normal()), prop = cur.map((v, a) => { let s = v; for (let b = 0; b <= a; b++) s += Lc[a][b] * z[b]; return s; }), l = lp(prop); if (Math.log(Math.max(rnd(), 1e-300)) < l - lc) { cur = prop; lc = l; acc++; } if (i >= burn) ch.push(cur); }
    const cols = p.map((_, j) => ch.map(c => c[j])), halfN = Math.floor(ch.length / 2);
    const rhat = cols.map(c => { const a = c.slice(0, halfN), b = c.slice(halfN, 2 * halfN), W = (N.variance(a) + N.variance(b)) / 2, B = halfN * N.variance([N.mean(a), N.mean(b)]); return W > 0 ? Math.sqrt(((halfN - 1) / halfN * W + B / halfN) / W) : 1; });
    const thin = Math.max(1, Math.floor(ch.length / 300)), sub = ch.filter((_, i) => i % thin === 0), sbar = hasS ? Math.sqrt(N.mean(S.map(s => s * s)) * (inp.scaleCov ? Math.max(red, 1e-300) : 1)) : Math.sqrt(red);
    const lo = [], hi = [], plo = [], phi = []; xs.forEach(x => { const v = sub.map(c => fam.f(x, c)), vp = v.map(u => u + sbar * rnd.normal()); lo.push(N.quantile(v, 0.025)); hi.push(N.quantile(v, 0.975)); plo.push(N.quantile(vp, 0.025)); phi.push(N.quantile(vp, 0.975)); });
    const inside = X.filter((x, i) => { const a = N.interp1(xs, plo, x), b = N.interp1(xs, phi, x); return Y[i] >= a && Y[i] <= b; }).length / m;
    post = { acc: acc / nMc, med: cols.map(c => N.quantile(c, 0.5)), lo: cols.map(c => N.quantile(c, 0.025)), hi: cols.map(c => N.quantile(c, 0.975)), rhat, inside, n: ch.length };
    if (post.acc < 0.08 || post.acc > 0.7) w.push(`MCMC acceptance rate ${pct(post.acc)} is outside the efficient range (about 10–50 %): the posterior is poorly scaled or bounded by the priors; treat credible intervals with caution.`);
    if (N.maxOf(rhat) > 1.1) w.push(`Split-chain R̂ = ${F(N.maxOf(rhat), 3)} exceeds 1.1: the chain has not converged; increase the number of MCMC samples.`);
    series.push({ key: 'ppc', title: 'Posterior predictive band', xlabel: head[0] || 'x', ylabel: head[1] || 'y', lines: [{ name: 'observations', x: X, y: Y }, { name: '95 % credible band of the mean (lower)', x: xs, y: lo }, { name: '95 % credible band of the mean (upper)', x: xs, y: hi }, { name: '95 % predictive band (lower)', x: xs, y: plo }, { name: '95 % predictive band (upper)', x: xs, y: phi }], explain: 'Bayesian credible band of the model mean and the wider predictive band that also contains measurement noise. About 95 % of the observations should fall inside the predictive band if the error model is adequate.' });
    cols.slice(0, 3).forEach((c, j) => { const h = N.histogram(c, 24); series.push({ key: 'post' + j, title: `Posterior distribution of ${fam.names[j]}`, xlabel: fam.names[j], ylabel: 'samples', lines: [{ name: 'posterior histogram', x: h.x, y: h.y }], explain: 'Histogram of the Markov-chain samples after burn-in. A narrow, single-peaked distribution inside the prior bounds indicates a well-informed parameter; a flat or truncated one indicates prior dominance.' }); });
  }
  // optional chemistry record
  const eqs = [{ kind: 'constitutive', label: 'Calibrated model', text: fam.eq }, { kind: 'formula', label: 'Weighted least squares (maximum likelihood for Gaussian errors)', text: 'χ²(θ) = Σᵢ [(yᵢ − f(xᵢ;θ))/σᵢ]²  →  min' }, { kind: 'formula', label: 'Fisher information and covariance', text: 'F = JᵀWJ,  Cov(θ) ≈ s²·F⁻¹,  s² = χ²/(m − n),  J = ∂f/∂θ' }, { kind: 'formula', label: 'Collinearity index', text: 'γ = 1/√λ_min(S̃ᵀS̃), S̃ = sensitivity matrix with unit-length columns' }];
  if (inp.bayes) eqs.push({ kind: 'formula', label: 'Bayesian posterior', text: 'p(θ|y) ∝ exp(−χ²(θ)/2s²)·p(θ),  p(θ) uniform on the physical bounds' });
  const rx = String(inp.rxn || '').split(/->|→|=/);
  if (rx.length === 2 && rx[0].trim() && rx[1].trim()) { try { const re = rx[0].split('+').map(s => s.trim()).filter(Boolean), pr = rx[1].split('+').map(s => s.trim()).filter(Boolean), b = chem.balance(re, pr); eqs.push({ kind: 'chemical', label: 'Reaction whose kinetics are calibrated', text: chem.equationString(re, pr, b.coefficients), note: 'Balanced from the reaction text supplied; the calibrated parameters belong to this reaction under the stated conditions only.' }); } catch (e) { w.push(`Reaction text could not be balanced (${e.message}).`); } }
  const tbl = [{ key: 'par', title: 'Calibrated parameters', columns: ['Parameter', 'Estimate', 'Std. error', 'Rel. error (%)', '95 % CI low', '95 % CI high', ...(post ? ['Posterior median', '95 % CrI low', '95 % CrI high', 'R̂'] : [])], rows: p.map((v, j) => [fam.names[j], F(v, 5), F(fit.se[j], 3), F(100 * rel[j], 3), F(v - tq * fit.se[j], 5), F(v + tq * fit.se[j], 5), ...(post ? [F(post.med[j], 5), F(post.lo[j], 5), F(post.hi[j], 5), F(post.rhat[j], 4)] : [])]).concat((fam.derived ? fam.derived(p) : []).map(d => [d[0], F(d[1], 5), '–', '–', '–', '–', ...(post ? ['–', '–', '–', '–'] : [])])), explain: 'Least-squares estimates with standard errors and t-based 95 % confidence intervals; with Bayesian calibration the posterior median and 95 % credible interval are added. Disagreement between the two interval types signals non-Gaussian or prior-limited parameters.' },
    { key: 'corr', title: 'Parameter correlation matrix', columns: ['', ...fam.names], rows: fit.corr.map((r, a) => [fam.names[a], ...r.map(v => +v.toFixed(5))]), explain: 'Correlation coefficients of the parameter estimates. Values near ±1 mean that two parameters compensate each other and cannot be estimated independently from these data.' },
    { key: 'fim', title: 'Identifiability diagnostics', columns: ['Diagnostic', 'Value'], rows: [...fit.eig.map((v, k) => [`Eigenvalue ${k + 1} of the scaled Fisher information`, F(v, 4)]), ['Condition number', F(fit.cond, 4)], ['Collinearity index γ', F(fit.gamma, 4)], [`Profile-likelihood 95 % interval of ${fam.names[worst]}`, `${pl === null ? 'open' : F(pl, 5)} … ${ph === null ? 'open' : F(ph, 5)}`], ['Most informative next measurement (D-optimal) at x', F(xBest, 5)]], explain: 'Eigenvalues of the Fisher information scaled to unit diagonal: a smallest eigenvalue near zero (large condition number, γ above roughly 10–15) marks a practically non-identifiable parameter combination. The D-optimal point maximises the gain in det(F).' }];
  const okId = fit.identifiable && fit.gamma < 15 && profOk;
  return {
    summary: `Weighted nonlinear least squares of the ${fam.label.split('  ')[0]} model on ${m} observations gives ${p.map((v, j) => `${fam.names[j]} = ${F(v, 4)} ± ${F(fit.se[j], 2)}`).join(', ')} (reduced χ² = ${F(red, 3)}, R² = ${F(r2, 4)}). The collinearity index is ${F(fit.gamma, 3)} and the largest parameter correlation ${F(maxCorr, 4)}, so the parameters are ${okId ? 'practically identifiable' : 'only weakly identifiable'} from these data.${post ? ` Metropolis sampling (acceptance ${pct(post.acc)}) puts ${pct(post.inside)} of the observations inside the 95 % predictive band.` : ''} This is a calibration result; it does not establish predictive validity.`,
    kpis: [
      ...p.slice(0, 3).map((v, j) => kp('p' + j, fam.names[j], v, '', `Least-squares estimate of ${fam.names[j]} with standard error ${F(fit.se[j], 3)} (${F(100 * rel[j], 3)} %); 95 % confidence interval ${F(v - tq * fit.se[j], 4)} to ${F(v + tq * fit.se[j], 4)}.`, rel[j] < 0.3)),
      kp('chi2red', 'Reduced χ²', red, '–', 'Weighted sum of squared residuals divided by the degrees of freedom. A value near 1 means the scatter is consistent with the stated measurement uncertainty; much larger values indicate model-form error or under-estimated sigmas.', hasS ? red < 3 && red > 0.1 : null),
      kp('r2', 'Coefficient of determination R²', r2, '–', 'Fraction of the variance of the observations reproduced by the fitted model. It measures goodness of fit to the calibration data only and says nothing about independent predictions.', r2 > 0.9),
      kp('rmse', 'Calibration RMSE', rm, head[1] || '', `Root-mean-square residual of the fit, equal to ${F(100 * rm / ym, 3)} % of the mean observed magnitude.`, null),
      kp('gamma', 'Collinearity index γ', fit.gamma, '–', 'Brun-type index of near-linear dependence between parameter sensitivities; values above roughly 10–15 mean a change in one parameter can be almost fully compensated by the others.', fit.gamma < 15),
      kp('cond', 'Condition number of Fisher information', fit.cond, '–', 'Ratio of the largest to smallest eigenvalue of the scaled Fisher information matrix; large values quantify how unevenly the data inform the parameter combinations.', fit.cond < 1e4),
      kp('maxcorr', 'Largest parameter correlation |ρ|', maxCorr, '–', 'Largest off-diagonal entry of the parameter correlation matrix. Close to 1 the two parameters are determined only as a combination.', maxCorr < 0.95),
      kp('worst_rel', 'Largest relative standard error', 100 * rel[worst], '%', `Relative standard error of the least-determined parameter (${fam.names[worst]}); this is the figure passed to the calibration gate of the credibility assessment.`, rel[worst] < 0.3),
      kp('x_next', 'Most informative next measurement', xBest, head[0] || '', 'Condition at which one additional measurement gives the largest increase of the determinant of the Fisher information (local D-optimal design within the scanned range).', null),
      ...(post ? [kp('acc', 'MCMC acceptance rate', 100 * post.acc, '%', 'Fraction of proposed Metropolis moves accepted. Rates of roughly 10–50 % indicate an efficiently scaled random walk.', post.acc > 0.08 && post.acc < 0.7), kp('ppc', 'Observations inside 95 % predictive band', 100 * post.inside, '%', 'Posterior predictive check on the calibration data: with an adequate error model about 95 % of the observations lie inside the band.', post.inside >= 0.8)] : []),
    ],
    series, tables: tbl, equations: eqs,
    checks: [
      { name: 'Degrees of freedom available', pass: m > n, detail: `${m} observations, ${n} parameters` },
      { name: 'Fisher information non-singular', pass: fit.identifiable, detail: `scaled eigenvalues ${fit.eig.map(v => F(v, 3)).join(', ')}` },
      { name: 'Practical identifiability (collinearity index < 15)', pass: fit.gamma < 15, detail: `γ = ${F(fit.gamma, 4)}` },
      { name: 'Profile likelihood bounded on both sides', pass: profOk, detail: `${fam.names[worst]}: ${pl === null ? 'open' : F(pl, 4)} … ${ph === null ? 'open' : F(ph, 4)}` },
      { name: 'Residuals consistent with stated uncertainty', pass: !hasS || (red < 3 && red > 0.1), detail: hasS ? `reduced χ² = ${F(red, 3)}` : 'unweighted fit: noise estimated from residuals' },
      { name: 'Parameters inside physical bounds', pass: p.every((v, j) => v > fam.b[j][0] && v < fam.b[j][1]), detail: p.map((v, j) => `${fam.names[j]} = ${F(v, 4)}`).join('; ') },
      ...(post ? [{ name: 'MCMC chain convergence (split R̂ < 1.1)', pass: N.maxOf(post.rhat) < 1.1, detail: `R̂ max = ${F(N.maxOf(post.rhat), 4)} over ${post.n} samples` }] : []),
    ],
    warnings: w,
    recommendations: [
      { priority: okId ? 'medium' : 'high', text: `Measure next at ${head[0] || 'x'} ≈ ${F(xBest, 4)}: of the conditions scanned this gives the largest gain in parameter information (det F increases by a factor ${F(N.maxOf(gain), 3)}).`, basis: 'Local D-optimal experimental design from the Fisher information matrix.' },
      { priority: maxCorr > 0.95 ? 'high' : 'low', text: maxCorr > 0.95 ? (inp.cmodel === 'arrhenius' ? 'Widen the temperature range or re-parameterise the Arrhenius law about a reference temperature (k = k_ref·exp[−Eₐ/R·(1/T − 1/T_ref)]) to decouple the pre-exponential factor from the activation energy.' : 'Add observations under conditions where the correlated parameters act differently (different range of x, an additional observable) or fix one of them from an independent measurement.') : 'Parameter correlations are moderate; keep the covariance matrix with the calibrated record and propagate it, not only the standard errors.', basis: `Largest parameter correlation |ρ| = ${F(maxCorr, 4)}.` },
      { priority: 'high', text: 'Reserve independent measurements, not used here, and run the validation analysis before using the calibrated model for prediction; agreement with the fitting data is not validation.', basis: 'Calibration and validation datasets must be separate (CVVUQ framework).' },
      { priority: 'medium', text: `Plan the remaining experiments with the D-optimal ranking instead of a uniform grid: the same parameter precision is reached with fewer syntheses, reducing precursor, solvent and energy consumption of the calibration campaign.`, basis: 'Information gain per experiment from the Fisher information.', sustainability: true },
    ],
    publish: { 'cred.cal_rmse': rm, 'cred.cal_rmse_pct': 100 * rm / ym, 'cred.max_rel_se_pct': 100 * rel[worst], 'cred.collinearity': fit.gamma, 'cred.chi2_red': red },
    exports: [{ name: 'calibration_record.json', mime: 'application/json', content: JSON.stringify({ model: fam.eq, parameters: fam.names.map((nm, j) => ({ name: nm, estimate: p[j], stderr: fit.se[j], ci95: [p[j] - tq * fit.se[j], p[j] + tq * fit.se[j]] })), covariance: fit.cov, correlation: fit.corr, chi2: fit.chi2, dof: fit.dof, observations: m, seed: inp.seed, weighted: hasS, covarianceScaledByReducedChi2: !!(inp.scaleCov || !hasS) }, null, 1), description: 'Reproducible calibration record: model, parameter estimates, covariance, fit statistics and settings.' }],
    numerics: { method: 'Levenberg–Marquardt weighted nonlinear least squares; central-difference sensitivities; Jacobi eigen-decomposition of the Fisher information' + (inp.bayes ? '; random-walk Metropolis MCMC' : ''), dof: n, steps: post ? post.n : 0, iterations: fit.iters, residual: fit.chi2, notes: 'Parameters are scaled by their initial guesses during optimisation. Confidence intervals are linearised (Wald) intervals; the profile likelihood and MCMC test their adequacy.' },
  };
}

// ---------- mode: validation ----------
function runValid(inp) {
  const w = [], t = parseTable(inp.vdata, 2, 3, 'Validation table'), X = t.rows.map(r => r[0]), M = t.rows.map(r => r[1]), n = X.length;
  const um = t.rows.map(r => Number.isFinite(r[2]) && r[2] >= 0 ? r[2] : 0); if (um.every(v => v === 0)) w.push('No measurement uncertainties supplied: uncertainty-based agreement metrics use prediction uncertainty only and are optimistic about the data.');
  if (inp.vdata === VAL_DEF) w.push('Validation table: synthetic demonstration data — replace with measurements.');
  const hasPred = t.rows.every(r => Number.isFinite(r[3])); let P, up, src = inp.predSrc, cal = null, reuse = 0, calN = 0;
  try { const cw = []; cal = readCal(inp, cw); calN = cal.X.length; reuse = X.filter((x, i) => cal.X.some((cx, k) => Math.abs(cx - x) <= 1e-9 * Math.max(1, Math.abs(x)) && Math.abs(cal.Y[k] - M[i]) <= 1e-9 * Math.max(1e-300, Math.abs(M[i])))).length; if (src === 'model') w.push(...cw); } catch (e) { if (src === 'model') { if (!hasPred) throw e; src = 'table'; w.push('Calibration table unusable (' + e.message + '); predictions taken from the validation table.'); } }
  if (src === 'table' && !hasPred) { src = 'model'; w.push('The validation table has no complete "predicted" column: predictions were computed with the model calibrated on the calibration table.'); }
  const fam = FAM[inp.cmodel] || FAM.arrhenius;
  if (src === 'model') { const fit = calibrate(cal.X, cal.Y, cal.S, fam, inp.scaleCov || !cal.hasS); P = X.map(x => fam.f(x, fit.p)); up = X.map(x => { const g = fit.p.map((v, j) => { const h = 1e-5 * Math.max(Math.abs(v), 1e-3 * fit.sc[j]), a = [...fit.p], b = [...fit.p]; a[j] += h; b[j] -= h; const d = (fam.f(x, a) - fam.f(x, b)) / (2 * h); return Number.isFinite(d) ? d : 0; }); return Math.sqrt(Math.max(0, N.dot(g, N.matVec(fit.cov, g)))); }); }
  else { P = t.rows.map(r => r[3]); up = t.rows.map(r => Number.isFinite(r[4]) && r[4] >= 0 ? r[4] : 0); }
  P = P.map(v => Number.isFinite(v) ? v : 0); up = up.map(v => Number.isFinite(v) ? Math.min(v, 1e150) : 0);
  const e = P.map((v, i) => v - M[i]), mm = Math.max(N.mean(M.map(Math.abs)), 1e-300), rm = N.rmse(P, M), ma = N.mae(P, M), bias = N.mean(e), r2 = N.r2(P, M), emax = N.maxOf(e.map(Math.abs));
  const uc = um.map((v, i) => Math.hypot(v, up[i])), hasU = uc.some(v => v > 0), k = inp.kCov, En = e.map((v, i) => uc[i] > 0 ? Math.abs(v) / (k * uc[i]) : (v === 0 ? 0 : 99));
  const within = En.filter(v => v <= 1).length / n, chi = hasU ? N.mean(e.map((v, i) => uc[i] > 0 ? (v / uc[i]) ** 2 : 0)) : 0;
  const cover = up.some(v => v > 0) ? e.filter((v, i) => Math.abs(v) <= 1.96 * Math.hypot(up[i], um[i])).length / n : 0;
  const se = N.std(e) / Math.sqrt(n), tstat = se > 0 ? bias / se : 0, biasSig = Math.abs(tstat) > t975(n - 1);
  const nr = 100 * rm / mm, nb = 100 * Math.abs(bias) / mm, okR = nr <= inp.accRmse, okB = nb <= inp.accBias, okC = !hasU || 100 * within >= inp.accCov, identical = reuse === n && n > 0, indep = reuse === 0;
  const meets = okR && okB && okC && !identical, x0 = N.minOf(X), x1 = N.maxOf(X), rng = (x1 - x0) || 1, ext = inp.xApp < x0 ? (x0 - inp.xApp) / rng : inp.xApp > x1 ? (inp.xApp - x1) / rng : 0;
  if (identical) w.push('The validation table is identical to the calibration table: this is a goodness-of-fit check, not independent validation, and cannot support a validated classification.');
  else if (reuse > 0) w.push(`${reuse} of ${n} validation points also appear in the calibration table: remove them, since reused data cannot count as independent evidence.`);
  if (ext > 0) w.push(`The application condition x = ${F(inp.xApp, 5)} lies ${pct(ext)} of the validated range outside it: this is extrapolation, for which wider uncertainty alone does not demonstrate validity.`);
  if (n < 5) w.push(`Only ${n} validation points: all metrics have large sampling uncertainty.`);
  const lim = [Math.min(N.minOf(M), N.minOf(P)), Math.max(N.maxOf(M), N.maxOf(P))], hz = N.histogram(e.map((v, i) => uc[i] > 0 ? v / uc[i] : 0), Math.min(12, Math.max(4, Math.round(Math.sqrt(n)) + 2)));
  return {
    summary: `Comparison of ${n} ${indep ? 'independent ' : ''}validation points with predictions from ${src === 'model' ? 'the model calibrated on the calibration table' : 'the table'}: RMSE ${F(rm, 4)} (${F(nr, 3)} % of the mean measured value), bias ${F(bias, 3)} (${F(nb, 3)} %), R² = ${F(r2, 4)}; ${pct(within)} of the points agree within ${k}·u_c and the reduced χ² is ${F(chi, 3)}. Against the stated acceptance criteria (RMSE ≤ ${inp.accRmse} %, |bias| ≤ ${inp.accBias} %, coverage ≥ ${inp.accCov} %) the prediction ${meets ? 'meets' : 'does not meet'} the requirement over x = ${F(x0, 4)}…${F(x1, 4)}. The statement applies to these criteria and this domain only.`,
    kpis: [
      kp('rmse', 'RMSE', rm, t.head[1] || '', `Root-mean-square prediction error over the validation set, ${F(nr, 3)} % of the mean measured magnitude; compared with the acceptance value of ${inp.accRmse} %.`, okR),
      kp('nrmse', 'Normalised RMSE', nr, '%', 'RMSE divided by the mean absolute measured value. This is the figure passed to the validation gate of the credibility assessment.', okR),
      kp('mae', 'Mean absolute error', ma, t.head[1] || '', 'Average magnitude of the prediction error; less sensitive to single outliers than the RMSE, so a large RMSE/MAE ratio points to a few large deviations.', null),
      kp('bias', 'Bias (predicted − measured)', bias, t.head[1] || '', `Mean signed error, ${F(nb, 3)} % of the mean measured value; t = ${F(tstat, 3)}, so the bias is ${biasSig ? 'statistically significant at the 5 % level' : 'not statistically distinguishable from zero'}.`, okB),
      kp('r2', 'R² of prediction', r2, '–', 'Coefficient of determination of the predictions against the independent measurements (1 − SS_res/SS_tot); it can be negative when the model is worse than the mean of the data.', r2 > 0.8),
      kp('within', `Fraction within ${k}·u_c`, 100 * within, '%', 'Share of validation points whose error is smaller than the expanded combined uncertainty of measurement and prediction (normalised error Eₙ ≤ 1).', okC),
      kp('chi2', 'Reduced χ² of validation', chi, '–', 'Mean squared error in units of the combined standard uncertainty. Values near 1 mean errors are explained by stated uncertainties; values well above 1 reveal unexplained model discrepancy.', hasU ? chi < 3 : null),
      kp('cover', 'Coverage of 95 % prediction interval', 100 * cover, '%', 'Observed fraction of measurements inside the 95 % interval built from prediction and measurement uncertainty; it should be close to the nominal 95 % if the uncertainty statement is reliable.', up.some(v => v > 0) ? cover >= 0.8 : null),
      kp('emax', 'Largest absolute error', emax, t.head[1] || '', 'Worst single deviation in the validation set; relevant when the consequence of error is governed by the extreme rather than the average.', null),
      kp('extrap', 'Extrapolation of application point', 100 * ext, '%', 'Distance of the stated application condition outside the range of the validation data, as a percentage of that range; 0 means the model would be used inside the validated domain.', ext === 0),
      kp('verdict', 'Assessment against acceptance criteria', identical ? 'not independent' : meets ? 'criteria met' : 'criteria not met', '', 'Outcome relative to the user-defined acceptance criteria and validated domain; it is not an absolute statement of model validity and does not transfer to other conditions.', meets),
    ],
    series: [
      { key: 'parity', title: 'Predicted versus measured', xlabel: 'measured', ylabel: 'predicted', lines: [{ name: 'validation points', x: M, y: P }, { name: 'perfect agreement', x: lim, y: lim }], explain: 'Parity plot: points on the diagonal are predicted exactly. A systematic offset or slope different from one reveals bias or a wrong sensitivity of the model.' },
      { key: 'resx', title: 'Prediction error across the validation domain', xlabel: t.head[0] || 'x', ylabel: 'predicted − measured', lines: [{ name: 'error', x: X, y: e }, { name: `+${k}·u_c`, x: X, y: uc.map(v => k * v) }, { name: `−${k}·u_c`, x: X, y: uc.map(v => -k * v) }], explain: 'Signed error at each validation condition with the expanded combined uncertainty band. A trend with x shows where in the domain the model degrades and therefore where its applicability ends.' },
      { key: 'hist', title: 'Distribution of normalised residuals', xlabel: '(predicted − measured)/u_c', ylabel: 'count', lines: [{ name: 'normalised residuals', x: hz.x, y: hz.y }], explain: 'Histogram of the errors divided by their combined standard uncertainty. For a reliable uncertainty statement it is centred on zero with a spread of about one.' },
    ],
    tables: [{ key: 'pts', title: 'Point-by-point comparison', columns: [t.head[0] || 'x', 'Measured', 'u_meas', 'Predicted', 'u_pred', 'Error', 'Eₙ'], rows: X.map((x, i) => [x, F(M[i], 5), F(um[i], 3), F(P[i], 5), F(up[i], 3), F(e[i], 3), F(En[i], 3)]), explain: 'Measured and predicted values with uncertainties and the normalised error Eₙ = |error|/(k·u_c); points with Eₙ above 1 disagree by more than the uncertainties explain.' },
      { key: 'acc', title: 'Acceptance criteria', columns: ['Criterion', 'Result', 'Limit', 'Status'], rows: [['Normalised RMSE (%)', F(nr, 4), inp.accRmse, okR ? 'met' : 'not met'], ['|Bias| (%)', F(nb, 4), inp.accBias, okB ? 'met' : 'not met'], [`Fraction within ${k}·u_c (%)`, F(100 * within, 4), inp.accCov, hasU ? (okC ? 'met' : 'not met') : 'no uncertainties'], ['Independence from calibration data', `${reuse} of ${n} points reused`, 0, indep ? 'met' : 'not met'], ['Application inside validated domain (%)', F(100 * ext, 3), 0, ext === 0 ? 'met' : 'extrapolation']], explain: 'Each validation metric against the limit defined for the intended application. Acceptance limits are engineering decisions that depend on measurement capability and consequence of error.' }],
    equations: [
      { kind: 'formula', label: 'Error metrics', text: 'RMSE = √(Σeᵢ²/n),  MAE = Σ|eᵢ|/n,  bias = Σeᵢ/n,  eᵢ = ŷᵢ − yᵢ' },
      { kind: 'formula', label: 'Combined uncertainty and normalised error', text: 'u_c² = u_meas² + u_pred²,   Eₙ = |e|/(k·u_c),   χ²_red = (1/n)·Σ(eᵢ/u_c,ᵢ)²' },
      { kind: 'formula', label: 'Coefficient of determination', text: 'R² = 1 − Σ(yᵢ − ŷᵢ)²/Σ(yᵢ − ȳ)²' },
      ...(src === 'model' ? [{ kind: 'constitutive', label: 'Model used for prediction', text: fam.eq, note: 'Parameters from the calibration table; prediction uncertainty u_pred² = gᵀ·Cov(θ)·g.' }] : []),
    ],
    checks: [
      { name: 'Validation data independent of calibration data', pass: indep, detail: `${reuse} of ${n} validation points found in the ${calN}-row calibration table` },
      { name: 'Normalised RMSE within acceptance', pass: okR, detail: `${F(nr, 4)} % vs ${inp.accRmse} %` },
      { name: 'Bias within acceptance', pass: okB, detail: `${F(nb, 4)} % vs ${inp.accBias} %` },
      { name: 'Agreement within combined uncertainty', pass: okC, detail: hasU ? `${pct(within)} within ${k}·u_c vs ${inp.accCov} %` : 'no uncertainties available' },
      { name: 'Application point inside validated domain', pass: ext === 0, detail: `x_app = ${F(inp.xApp, 5)}, validated ${F(x0, 5)} … ${F(x1, 5)}` },
      { name: 'Sample size adequate for the metrics', pass: n >= 5, detail: `${n} validation points` },
    ],
    warnings: w,
    recommendations: [
      { priority: meets ? 'medium' : 'high', text: meets ? `Record the model as meeting the stated criteria for x = ${F(x0, 4)}…${F(x1, 4)} with this parameter set and model version; any change of parameters or code requires the comparison to be repeated.` : `Do not use the model for decisions at this accuracy: ${!okR ? 'the RMSE exceeds the limit' : !okB ? 'the bias exceeds the limit' : identical ? 'the data are not independent' : 'too few points agree within uncertainty'}. Inspect the error-versus-x plot to decide between recalibration and a change of model form.`, basis: 'Outcome of the acceptance-criteria table.' },
      { priority: biasSig ? 'high' : 'low', text: biasSig ? `The bias of ${F(bias, 3)} is systematic (t = ${F(tstat, 3)}): look for a missing mechanism or an offset in the measurement chain before adding a discrepancy term.` : 'No statistically significant bias is detected; the remaining error is consistent with random scatter at this sample size.', basis: 'One-sample t-test on the signed errors.' },
      { priority: ext > 0 ? 'high' : 'medium', text: ext > 0 ? `Extend the validation experiments to cover x = ${F(inp.xApp, 5)} before using the model there; the present evidence ends at ${F(inp.xApp < x0 ? x0 : x1, 5)}.` : 'Add validation points near the edges of the intended operating window so that the documented domain fully contains the application.', basis: 'Validation-domain and extrapolation assessment.' },
      { priority: 'medium', text: 'Use the validated model to screen process conditions before running further syntheses, and place new experiments only where the error band is widest; this limits consumption of precursors, solvent and energy in the validation campaign.', basis: 'Error-versus-condition distribution of this comparison.', sustainability: true },
    ],
    publish: { 'cred.rmse': rm / mm, 'cred.rmse_pct': nr, 'cred.bias_pct': nb, 'cred.val_within': within, 'cred.extrap_pct': 100 * ext, 'cred.independent': indep ? 1 : 0 },
    numerics: { method: 'Direct evaluation of validation metrics' + (src === 'model' ? ' with Levenberg–Marquardt calibration and delta-method prediction uncertainty' : ''), dof: n, steps: 1, iterations: 1, residual: rm, notes: 'Errors are assumed independent between validation points; correlated measurement errors would reduce the effective sample size.' },
  };
}

// ---------- mode: uncertainty propagation and sensitivity ----------
function runUq(inp) {
  const w = [], mdl = UQM[inp.umodel] || UQM.arrh, d = mdl.vars.length, seed = Math.round(inp.seed), ov = {};
  for (const part of String(inp.uqOverride || '').split(';')) { const mt = part.match(/^\s*([A-Za-z_]\w*)\s*=\s*([-+0-9.eE]+)\s*(?:,\s*([-+0-9.eE]+))?/); if (!part.trim()) continue; if (!mt || !mdl.vars.some(v => v[0] === mt[1]) || !Number.isFinite(+mt[2])) { w.push(`Override "${part.trim()}" ignored: use name=mean,cv% with one of ${mdl.vars.map(v => v[0]).join(', ')}.`); continue; } ov[mt[1]] = [+mt[2], mt[3] !== undefined && Number.isFinite(+mt[3]) ? Math.abs(+mt[3]) / 100 : null]; }
  const V = mdl.vars.map(([name, unit, m0, cv0]) => { const m = ov[name] ? ov[name][0] : m0, cv = Math.min((ov[name] && ov[name][1] !== null ? ov[name][1] : cv0) * inp.cvScale, 1.5), s = Math.abs(m) * cv, dist = inp.dist === 'lognormal' && m <= 0 ? 'normal' : inp.dist, sg = Math.sqrt(Math.log(1 + cv * cv)); return { name, unit, m, cv, s, dist, sg, mu: m > 0 ? Math.log(m) - sg * sg / 2 : 0 }; });
  const tr = (z, v) => v.dist === 'lognormal' ? Math.exp(v.mu + v.sg * z) : v.dist === 'uniform' ? v.m + v.s * Math.sqrt(3) * (2 * N.normCdf(z) - 1) : (v.m > 0 ? Math.max(v.m + v.s * z, 1e-6 * v.m) : v.m + v.s * z);
  const fz = z => { const y = mdl.f(z.map((q, i) => tr(q, V[i]))); return Number.isFinite(y) ? y : 0; }, nom = fz(V.map(v => v.dist === 'lognormal' ? v.sg / 2 : 0)), y0 = mdl.f(V.map(v => v.m));
  const n = Math.round(N.clamp(inp.nSamp, 100, 20000)), rnd = N.rng(seed), rho = N.clamp(inp.rho12, -0.95, 0.95);
  const Z = inp.sampler === 'lhs' ? N.lhs(n, d, rnd).map(r => r.map(u => N.normInv(N.clamp(u, 1e-12, 1 - 1e-12)))) : Array.from({ length: n }, () => V.map(() => rnd.normal()));
  const corr = z => { if (!rho) return z; const c = [...z]; c[1] = rho * z[0] + Math.sqrt(1 - rho * rho) * z[1]; return c; }, Yv = Z.map(z => fz(corr(z)));
  const mean = N.mean(Yv), sd = N.std(Yv), q = [0.025, 0.5, 0.975].map(p => N.quantile(Yv, p)), limit = inp.limRatio * y0, above = inp.limSide !== 'below';
  const pf = Yv.filter(v => above ? v > limit : v < limit).length / n, sePf = Math.sqrt(Math.max(pf * (1 - pf), 1 / n) / n);
  // running mean and standard deviation
  const cn = [...new Set(N.logspace(Math.log10(20), Math.log10(n), 40).map(Math.round))].filter(k => k >= 2 && k <= n), cm = [], cs = []; { let s = 0, s2 = 0, j = 0; for (let i = 0; i < n; i++) { s += Yv[i]; s2 += Yv[i] * Yv[i]; if (i + 1 === cn[j]) { const mu = s / (i + 1); cm.push(mu); cs.push(Math.sqrt(Math.max(0, (s2 - (i + 1) * mu * mu) / i))); j++; } } }
  // Saltelli scheme for first-order (Saltelli 2010) and total (Jansen) Sobol indices: independent inputs
  const Ns = N.clamp(2 * n, 2000, 6000), r2 = N.rng(seed + 101), A = Array.from({ length: Ns }, () => V.map(() => r2.normal())), B = Array.from({ length: Ns }, () => V.map(() => r2.normal())), fA = A.map(fz), fB = B.map(fz);
  const all = fA.concat(fB), VY = N.variance(all), f0 = N.mean(all), S1 = [], ST = [];
  for (let i = 0; i < d; i++) { let s1 = 0, st = 0; for (let k = 0; k < Ns; k++) { const ab = [...A[k]]; ab[i] = B[k][i]; const fab = fz(ab); s1 += (fB[k] - f0) * (fab - fA[k]); st += (fA[k] - fab) ** 2; } S1.push(VY > 0 ? s1 / Ns / VY : 0); ST.push(VY > 0 ? st / (2 * Ns) / VY : 0); }
  const flat = !(VY > 1e-24 * Math.max(mean * mean, 1e-300)); if (flat) { w.push('Output variance is zero (no input uncertainty): sensitivity indices cannot be defined and are reported as zero.'); S1.fill(0); ST.fill(0); }
  // Morris elementary effects (r trajectories, 4 levels, Δ = 2/3) on the 2.5–97.5 % quantile range
  const r3 = N.rng(seed + 202), rT = 30, EE = V.map(() => []), u2z = u => N.normInv(0.025 + 0.95 * u);
  for (let tI = 0; tI < rT; tI++) { let u = V.map(() => r3.int(2) / 3), fy = fz(u.map(u2z)); const order = V.map((_, i) => i); for (let i = d - 1; i > 0; i--) { const k = r3.int(i + 1);[order[i], order[k]] = [order[k], order[i]]; } for (const i of order) { const un = [...u]; un[i] += 2 / 3; const fn = fz(un.map(u2z)); EE[i].push((fn - fy) / (2 / 3)); u = un; fy = fn; } }
  const muS = EE.map(a => N.mean(a.map(Math.abs))), sgE = EE.map(a => N.std(a));
  // degree-2 Hermite polynomial chaos by least-squares regression on the standard-normal germ
  const basis = z => { const b = [1]; for (let i = 0; i < d; i++) b.push(z[i]); for (let i = 0; i < d; i++) b.push(z[i] * z[i] - 1); for (let i = 0; i < d; i++) for (let j = i + 1; j < d; j++) b.push(z[i] * z[j]); return b; };
  const nb = 1 + 2 * d + d * (d - 1) / 2, nTr = Math.min(Ns, 1200), sclY = Math.max(Math.sqrt(VY), Math.abs(N.mean(all)), 1e-300); let c = new Array(nb).fill(0), pceOk = true;
  try { c = N.lstsq(A.slice(0, nTr).map(basis), fA.slice(0, nTr).map(v => v / sclY), 1e-10).map(v => v * sclY); if (!c.every(Number.isFinite)) throw 0; } catch { pceOk = false; c = new Array(nb).fill(0); c[0] = N.mean(fA); w.push('Polynomial-chaos regression was singular; surrogate reduced to the sample mean.'); }
  const pMean = c[0]; let pVar = 0; const pS1 = V.map(() => 0), pST = V.map(() => 0);
  for (let i = 0; i < d; i++) { const v = c[1 + i] ** 2 + 2 * c[1 + d + i] ** 2; pS1[i] += v; pST[i] += v; pVar += v; } { let k = 1 + 2 * d; for (let i = 0; i < d; i++) for (let j = i + 1; j < d; j++) { const v = c[k++] ** 2; pST[i] += v; pST[j] += v; pVar += v; } }
  const pR2 = N.r2(B.slice(0, 600).map(z => N.dot(basis(z), c)), fB.slice(0, 600)), pSd = Math.sqrt(pVar), sumS = N.sum(S1), rank = V.map((_, i) => i).sort((a, b) => ST[b] - ST[a]), top = rank[0];
  const mcI = N.mean(fA), sdI = Math.sqrt(N.variance(fA));
  if (rho) w.push(`Inputs 1 and 2 are correlated (ρ = ${rho}) in the propagation; the Sobol, Morris and polynomial-chaos results are computed for independent inputs and do not include the correlation.`);
  if (V.some(v => v.dist === 'normal' && v.cv > 0.3 && v.m > 0)) w.push('A normal distribution with a coefficient of variation above 30 % reaches non-physical negative values; samples are truncated at 10⁻⁶ of the mean. Prefer the lognormal distribution.');
  if (pR2 < 0.9 && !flat) w.push(`The quadratic polynomial-chaos surrogate reproduces only R² = ${F(pR2, 3)} of independent samples: the response is strongly nonlinear over the input range and surrogate-based indices are approximate.`);
  const h = N.histogram(Yv, 30), cvOut = Math.abs(mean) > 0 ? sd / Math.abs(mean) : 0, sob = (v) => +N.clamp(v, -9, 9).toFixed(4);
  return {
    summary: `${inp.sampler === 'lhs' ? 'Latin-hypercube' : 'Monte Carlo'} propagation of ${d} uncertain inputs through the ${mdl.out.toLowerCase()} model (${n} samples) gives a mean of ${F(mean, 4)} ${mdl.unit} with standard deviation ${F(sd, 3)} ${mdl.unit} (CV ${pct(cvOut)}); the 95 % interval is ${F(q[0], 4)} to ${F(q[2], 4)} ${mdl.unit}. The probability that the output is ${above ? 'above' : 'below'} ${F(limit, 4)} ${mdl.unit} is ${pct(pf)} ± ${pct(sePf)}. ${flat ? 'No input uncertainty is active.' : `${V[top].name} dominates the variance (total Sobol index ${F(ST[top], 3)}); the first-order indices sum to ${F(sumS, 3)}.`}`,
    kpis: [
      kp('mean', `Mean ${mdl.out.toLowerCase()}`, mean, mdl.unit, `Sample mean of the model output over ${n} realisations; the nominal (mean-input) value is ${F(y0, 4)} ${mdl.unit}, and a difference between the two reflects nonlinearity of the model.`, null),
      kp('std', 'Standard deviation', sd, mdl.unit, 'Sample standard deviation of the output caused by the stated input uncertainties; it does not include model-form or numerical error.', null),
      kp('cv', 'Coefficient of variation', 100 * cvOut, '%', 'Output standard deviation relative to the mean: the relative prediction uncertainty that follows from the input scatter alone.', cvOut < 0.2),
      kp('p025', '2.5th percentile', q[0], mdl.unit, 'Lower bound of the central 95 % interval of the output distribution, from the empirical quantiles of the samples.', null),
      kp('p50', 'Median', q[1], mdl.unit, 'Median of the output distribution; for skewed (for example lognormal-like) outputs it differs from the mean.', null),
      kp('p975', '97.5th percentile', q[2], mdl.unit, 'Upper bound of the central 95 % interval of the output distribution.', null),
      kp('pfail', `Probability output ${above ? '>' : '<'} limit`, pf, '–', `Fraction of samples ${above ? 'above' : 'below'} the limit ${F(limit, 4)} ${mdl.unit}, with binomial standard error ${F(sePf, 2)}; probabilities below about 10/n are not resolved by this sample size.`, pf < 0.05),
      kp('s_top', `Total Sobol index of ${V[top].name}`, ST[top], '–', `Share of the output variance that involves ${V[top].name}, including its interactions (Jansen estimator on ${Ns} base samples). Reducing the uncertainty of this input gives the largest reduction of prediction uncertainty.`, null),
      kp('sum_s1', 'Sum of first-order Sobol indices', sumS, '–', 'For an additive model the first-order indices sum to 1; a smaller sum is the share of variance caused by interactions between inputs.', sumS <= 1.1),
      kp('pce_r2', 'Surrogate R² on independent samples', pR2, '–', 'Predictive coefficient of determination of the degree-2 Hermite polynomial-chaos surrogate on samples not used to fit it; near 1 the surrogate may replace the model for uncertainty propagation.', pR2 > 0.9),
      kp('pce_sd', 'Surrogate standard deviation', pSd, mdl.unit, `Output standard deviation obtained analytically from the polynomial-chaos coefficients (independent inputs), to be compared with the sampling value ${F(sdI, 4)} ${mdl.unit} as a cross-check.`, sdI > 0 ? Math.abs(pSd - sdI) / sdI < 0.15 : null),
      kp('se_mean', 'Standard error of the mean', sd / Math.sqrt(n), mdl.unit, 'Sampling uncertainty of the estimated mean; it falls with the square root of the number of samples and shows whether the sample size is sufficient.', null),
    ],
    series: [
      { key: 'hist', title: `Distribution of ${mdl.out.toLowerCase()}`, xlabel: `${mdl.out} (${mdl.unit})`, ylabel: 'samples', lines: [{ name: 'output histogram', x: h.x, y: h.y }, { name: 'limit', x: [limit, limit], y: [0, N.maxOf(h.y)] }], explain: 'Histogram of the propagated output with the specification limit. Skewness shows nonlinear response; the area beyond the limit is the estimated probability of failure.' },
      { key: 'conv', title: 'Convergence of mean and standard deviation', xlabel: 'number of samples', ylabel: mdl.unit, logx: true, lines: [{ name: 'running mean', x: cn.slice(0, cm.length), y: cm }, { name: 'running standard deviation', x: cn.slice(0, cs.length), y: cs }], explain: 'Running estimates of the output mean and standard deviation against sample size. Both curves must flatten before the statistics are quoted; residual drift indicates too few samples.' },
      { key: 'sobol', title: 'Global sensitivity indices', xlabel: 'input number', ylabel: 'Sobol index', lines: [{ name: 'first-order S₁ (Saltelli)', x: V.map((_, i) => i + 1), y: S1.map(sob) }, { name: 'total S_T (Jansen)', x: V.map((_, i) => i + 1), y: ST.map(sob) }, { name: 'first-order from polynomial chaos', x: V.map((_, i) => i + 1), y: pS1.map(v => pVar > 0 ? v / pVar : 0) }], explain: 'Variance-based sensitivity of the output to each input (numbered as in the input table). The gap between total and first-order index is the interaction contribution; agreement with the surrogate values cross-checks both estimates.' },
      { key: 'morris', title: 'Morris screening (μ* versus σ)', xlabel: 'μ* (mean absolute elementary effect)', ylabel: 'σ of elementary effects', lines: [{ name: 'inputs', x: muS, y: sgE }], explain: 'Each point is one input: large μ* marks an influential input, large σ relative to μ* marks nonlinearity or interaction. Points near the origin may be fixed at nominal values.' },
    ],
    tables: [{ key: 'in', title: 'Inputs and sensitivity ranking', columns: ['#', 'Input', 'Unit', 'Mean', 'CV (%)', 'Distribution', 'S₁', 'S_T', 'S₁ (PCE)', 'S_T (PCE)', 'Morris μ*', 'Morris σ'], rows: V.map((v, i) => [i + 1, v.name, v.unit, F(v.m, 5), F(100 * v.cv, 3), v.dist, sob(S1[i]), sob(ST[i]), sob(pVar > 0 ? pS1[i] / pVar : 0), sob(pVar > 0 ? pST[i] / pVar : 0), F(muS[i], 3), F(sgE[i], 3)]), explain: 'Uncertain inputs with their distributions and the sensitivity measures from three independent methods. Inputs with small total index and small μ* contribute little to prediction uncertainty.' },
      { key: 'x', title: 'Cross-check of propagation methods', columns: ['Quantity', 'Sampling (independent inputs)', 'Polynomial chaos', `Propagation as specified (ρ = ${rho})`], rows: [['Mean', F(mcI, 5), F(pMean, 5), F(mean, 5)], ['Standard deviation', F(sdI, 5), F(pSd, 5), F(sd, 5)]], explain: 'Mean and standard deviation from direct sampling and from the analytical moments of the polynomial-chaos expansion. Close agreement indicates that both the sample size and the surrogate order are adequate.' }],
    equations: [
      { kind: 'constitutive', label: 'Model propagated', text: mdl.eq },
      { kind: 'formula', label: 'Monte Carlo estimators', text: 'μ̂ = (1/n)·Σyᵢ,  σ̂² = Σ(yᵢ − μ̂)²/(n − 1),  P̂_f = (1/n)·Σ 1[yᵢ beyond limit],  SE(μ̂) = σ̂/√n' },
      { kind: 'formula', label: 'Sobol indices (Saltelli sampling)', text: 'Sᵢ = (1/N)·Σ f(B)ₖ·[f(A_B⁽ⁱ⁾)ₖ − f(A)ₖ]/V(Y),   S_Tᵢ = (1/2N)·Σ [f(A)ₖ − f(A_B⁽ⁱ⁾)ₖ]²/V(Y)' },
      { kind: 'formula', label: 'Polynomial chaos expansion (degree 2, Hermite)', text: 'Y ≈ c₀ + Σcᵢξᵢ + Σcᵢᵢ(ξᵢ² − 1) + Σ_{i<j} cᵢⱼξᵢξⱼ;   E[Y] = c₀,   V[Y] = Σcᵢ² + 2Σcᵢᵢ² + Σcᵢⱼ²' },
      { kind: 'formula', label: 'Morris elementary effect', text: 'EEᵢ = [f(u + Δeᵢ) − f(u)]/Δ,  μ*ᵢ = mean|EEᵢ|,  Δ = 2/3 on a 4-level grid' },
    ],
    checks: [
      { name: 'First-order Sobol indices sum ≤ 1', pass: sumS <= 1.1, detail: `ΣSᵢ = ${F(sumS, 4)} (estimator noise allows a small excess)` },
      { name: 'Total index ≥ first-order index', pass: ST.every((v, i) => v >= S1[i] - 0.08), detail: V.map((v, i) => `${v.name}: ${F(S1[i], 2)}/${F(ST[i], 2)}`).join('; ') },
      { name: 'Statistical convergence of the mean', pass: sd / Math.sqrt(n) <= 0.01 * Math.max(Math.abs(mean), sd), detail: `SE(mean) = ${F(sd / Math.sqrt(n), 3)} ${mdl.unit}` },
      { name: 'Surrogate cross-check of the standard deviation', pass: flat || Math.abs(pSd - sdI) <= 0.15 * sdI, detail: `sampling ${F(sdI, 4)} vs polynomial chaos ${F(pSd, 4)} ${mdl.unit}` },
      { name: 'Failure probability resolved by the sample size', pass: pf === 0 || pf === 1 || Math.min(pf, 1 - pf) * n >= 10, detail: `${Math.round(pf * n)} of ${n} samples beyond the limit` },
    ],
    warnings: w,
    recommendations: [
      { priority: 'high', text: flat ? 'Enter realistic input uncertainties: with zero scatter the propagation only reproduces the nominal value.' : `Tighten the control or measurement of ${V[top].name} first: it carries ${pct(N.clamp(ST[top], 0, 1))} of the output variance, so halving its uncertainty removes up to ${pct(N.clamp(0.75 * S1[top], 0, 1))} of the variance.`, basis: 'Variance decomposition by Sobol indices.' },
      { priority: 'medium', text: `Inputs ${rank.slice(-Math.max(1, Math.floor(d / 2))).map(i => V[i].name).join(', ')} have the smallest total indices; they can be fixed at nominal values in calibration and optimisation to reduce the dimension of the problem.`, basis: 'Factor fixing from total Sobol indices and Morris μ*.' },
      { priority: pf > 0.05 ? 'high' : 'low', text: pf > 0.05 ? `The probability of ${above ? 'exceeding' : 'falling below'} the limit is ${pct(pf)}: shift the nominal design by about ${F(Math.abs(N.quantile(Yv, above ? 0.95 : 0.05) - limit), 3)} ${mdl.unit} or reduce the dominant input scatter to bring it under 5 %.` : `The estimated probability beyond the limit is ${pct(pf)}; confirm small probabilities with a larger sample before using them as a reliability figure.`, basis: 'Empirical exceedance probability and its binomial standard error.' },
      { priority: 'medium', text: pR2 > 0.95 ? `Use the quadratic surrogate (R² = ${F(pR2, 4)}) for further propagation and optimisation loops: it replaces thousands of model runs and cuts the computing energy of robust-design studies.` : 'Before replacing the model by a surrogate in optimisation loops, raise its order or narrow the input ranges; an inaccurate surrogate wastes experimental and computing resources on false optima.', basis: 'Independent-sample accuracy of the polynomial-chaos surrogate.', sustainability: true },
    ],
    publish: { 'cred.uq_cv': cvOut, 'cred.uq_mean': mean, 'cred.uq_std': sd, 'cred.p_exceed': pf, 'cred.top_sobol': N.clamp(ST[top], 0, 1) },
    numerics: { method: `${inp.sampler === 'lhs' ? 'Latin hypercube' : 'Monte Carlo'} sampling with Gaussian copula; Saltelli/Jansen Sobol estimators; Morris one-at-a-time trajectories; least-squares Hermite polynomial chaos`, dof: d, steps: n, iterations: n + Ns * (d + 2) + rT * (d + 1), residual: sd / Math.sqrt(n), notes: `Sobol base sample ${Ns}; polynomial chaos with ${nb} terms fitted on ${nTr} samples${pceOk ? '' : ' (failed)'}; Morris with ${rT} trajectories.` },
  };
}

// ---------- mode: credibility assessment ----------
function runCred(inp) {
  const ev = [
    { ok: inp.g1 === 'complete', miss: inp.g1 === 'complete' ? [] : ['Complete the input data record: units, geometry integrity, material records, traceability and measurement uncertainty.'], num: null },
    { ok: !!inp.g2 && inp.consRes <= inp.consTol, miss: [...(inp.g2 ? [] : ['Review governing equations, constitutive models and admissibility checks.']), ...(inp.consRes <= inp.consTol ? [] : [`Reduce the conservation residual from ${F(inp.consRes, 2)} to below ${F(inp.consTol, 2)}.`])], num: ['Conservation residual', inp.consRes, inp.consTol] },
    { ok: !!inp.g3 && inp.ordDev <= inp.ordTolC, miss: [...(inp.g3 ? [] : ['Run analytical-benchmark or manufactured-solution tests for the solver version in use.']), ...(inp.ordDev <= inp.ordTolC ? [] : [`Resolve the order-of-accuracy deviation of ${F(inp.ordDev, 3)} (limit ${inp.ordTolC}).`])], num: ['Order deviation', inp.ordDev, inp.ordTolC] },
    { ok: !!inp.g4conv && inp.gciPct <= inp.gciTolC, miss: [...(inp.g4conv ? [] : ['Demonstrate iterative convergence (residual and update criteria, no stagnation).']), ...(inp.gciPct <= inp.gciTolC ? [] : [`Refine mesh/time step: GCI ${F(inp.gciPct, 3)} % exceeds ${F(inp.gciTolC, 3)} %.`])], num: ['GCI (%)', inp.gciPct, inp.gciTolC] },
    { ok: !!inp.g5 && inp.relSe <= inp.relSeTol, miss: [...(inp.g5 ? [] : ['Calibrate uncertain parameters against identified, archived datasets.']), ...(inp.relSe <= inp.relSeTol ? [] : [`Improve identifiability: parameter relative error ${F(inp.relSe, 3)} % exceeds ${F(inp.relSeTol, 3)} %.`])], num: ['Parameter relative error (%)', inp.relSe, inp.relSeTol] },
    { ok: !!inp.indep && inp.nVal >= 1 && inp.valRmse <= inp.accRmseC, miss: [...(inp.nVal >= 1 ? [] : ['Obtain at least one validation dataset that was not used for calibration.']), ...(inp.indep ? [] : ['Separate validation data from calibration data; reused data are not independent evidence.']), ...(inp.valRmse <= inp.accRmseC ? [] : [`Validation RMSE ${F(inp.valRmse, 3)} % exceeds the acceptance value ${F(inp.accRmseC, 3)} %.`])], num: ['Validation RMSE (%)', inp.valRmse, inp.accRmseC] },
    { ok: !!inp.g7uq && !!inp.g7sens && inp.extrap <= inp.extrapMax, miss: [...(inp.g7uq ? [] : ['Propagate input, parameter and numerical uncertainty to prediction intervals.']), ...(inp.g7sens ? [] : ['Carry out a global sensitivity analysis of the prediction.']), ...(inp.extrap <= inp.extrapMax ? [] : [`Application lies ${F(inp.extrap, 3)} % outside the validated domain (limit ${F(inp.extrapMax, 3)} %): extend validation.`])], num: ['Extrapolation (%)', inp.extrap, inp.extrapMax] },
    { ok: !!inp.g8use && !!inp.g8rev, miss: [...(inp.g8use ? [] : ['Document intended use, consequence of error, applicability limits and known deficiencies.']), ...(inp.g8rev ? [] : ['Obtain the required independent review and approval of the qualification record.'])], num: null },
  ];
  let level = 0; while (level < 8 && ev[level].ok) level++;
  const status = ev.map((g, i) => g.ok ? (i < level ? 'pass' : 'blocked by earlier gate') : (g.miss.length && i === 0 && inp.g1 === 'partial' ? 'incomplete' : 'fail'));
  const cls = level < 4 ? 0 : level === 4 ? 1 : level === 5 ? 2 : level < 8 ? 3 : 4, ownOk = ev.filter(g => g.ok).length, missing = ev.flatMap((g, i) => g.miss.map(m => [i + 1, GATES[i], m]));
  const nums = ev.map((g, i) => g.num ? [i + 1, ...g.num, g.num[2] > 0 ? g.num[1] / g.num[2] : (g.num[1] > 0 ? 99 : 0)] : null).filter(Boolean), first = ev.findIndex(g => !g.ok);
  const w = [];
  if (ev.some((g, i) => g.ok && i >= level)) w.push('Evidence exists for later gates, but they cannot be marked as passed while an earlier gate fails: later evidence rests on the earlier steps.');
  if (level < 6) w.push('Results of this model remain unvalidated predictions: independent experimental validation (gate 6) has not been passed.');
  if (inp.extrap > 0) w.push(`The intended application lies ${F(inp.extrap, 3)} % outside the validated domain; increased uncertainty alone is not evidence that the model remains physically valid there.`);
  const prio = i => i === first ? 'high' : i < 6 ? 'medium' : 'low', recs = missing.slice(0, 5).map(([g, name, m]) => ({ priority: prio(g - 1), text: `Gate ${g} (${name}): ${m}`, basis: g - 1 === first ? 'First failing gate: it blocks every later gate.' : 'Evidence still missing in the qualification sequence.' }));
  if (!recs.length) recs.push({ priority: 'low', text: 'All eight gates are supported by the declared evidence. Keep the record under version control and repeat gates 3–6 after any change of model, parameters, solver or software version.', basis: 'Validation claims are tied to specific model versions and parameter sets.' });
  recs.push({ priority: 'medium', text: 'Record the credibility class with every result exported from the platform so that downstream users see the evidence state, not only the numbers.', basis: 'Credibility labels are evidence-backed states, not marketing claims.' });
  recs.push({ priority: 'medium', text: level >= 6 ? 'Use the qualified model to replace part of the physical trial programme within its documented domain; each avoided synthesis or build saves precursors, solvent and energy.' : 'Close the open gates before committing to scale-up trials: decisions based on an unqualified model risk failed batches and the material, solvent and energy they consume.', basis: 'Resource cost of experiments versus the credibility level reached.', sustainability: true });
  recs.push({ priority: 'low', text: 'Define acceptance limits (GCI, validation RMSE, extrapolation) from the consequence of error of the specific application rather than adopting the default values.', basis: 'A universal percentage threshold is inappropriate across different phenomena.' });
  return {
    summary: `With the declared evidence ${level} consecutive qualification gate${level === 1 ? ' is' : 's are'} passed (${ownOk} of 8 have sufficient evidence of their own). The resulting credibility class is "${CLASSES[cls]}". ${first >= 0 ? `The first failing gate is gate ${first + 1} (${GATES[first]}); ${missing.length} item${missing.length === 1 ? ' of evidence is' : 's of evidence are'} still missing.` : 'No evidence item is missing.'} The classification reflects the entries made here and must be backed by the underlying records.`,
    kpis: [
      kp('class', 'Credibility class', CLASSES[cls], '', 'Evidence state of the model according to the five-level classification of the CVVUQ framework; it follows from the highest gate reached without any earlier gate failing.', cls >= 3),
      kp('level', 'Consecutive gates passed', level, 'of 8', 'Number of qualification gates passed in sequence starting from gate 1. A later gate cannot be credited while an earlier one fails.', level >= 6),
      kp('own', 'Gates with sufficient own evidence', ownOk, 'of 8', 'Number of gates whose individual evidence meets the stated criteria, irrespective of order; the difference to the consecutive count shows evidence that is currently blocked.', null),
      kp('missing', 'Missing evidence items', missing.length, '', 'Count of individual evidence items that are absent or outside their acceptance limits; they are listed in the missing-evidence table.', missing.length === 0),
      kp('gci_margin', 'GCI relative to acceptance', inp.gciTolC > 0 ? inp.gciPct / inp.gciTolC : 0, '–', 'Ratio of the numerical uncertainty (grid convergence index) to its acceptance limit; values below 1 satisfy the solution-verification criterion.', inp.gciPct <= inp.gciTolC),
      kp('rmse_margin', 'Validation RMSE relative to acceptance', inp.accRmseC > 0 ? inp.valRmse / inp.accRmseC : 0, '–', 'Ratio of the validation error to the acceptance limit set for the intended use; below 1 the validation criterion is met provided independent data exist.', inp.valRmse <= inp.accRmseC && inp.nVal >= 1),
      kp('nval', 'Independent validation datasets', inp.nVal, '', 'Number of datasets not used for calibration against which predictions were compared; without at least one the model cannot be classed as experimentally validated.', inp.nVal >= 1),
      kp('extrap', 'Extrapolation beyond validated domain', inp.extrap, '%', 'Distance of the application from the validated range relative to that range. Predictions outside the domain are flagged regardless of their uncertainty band.', inp.extrap <= inp.extrapMax),
    ],
    series: [
      { key: 'gates', title: 'Gate status', xlabel: 'gate number', ylabel: 'status (1 pass, 0.5 blocked, 0 fail)', lines: [{ name: 'status', x: ev.map((_, i) => i + 1), y: ev.map((g, i) => g.ok ? (i < level ? 1 : 0.5) : 0) }, { name: 'own evidence sufficient', x: ev.map((_, i) => i + 1), y: ev.map(g => g.ok ? 1 : 0) }], explain: 'Status of the eight qualification gates in sequence. A gate drawn at 0.5 has adequate evidence of its own but is blocked because an earlier gate fails.' },
      { key: 'margin', title: 'Quantitative evidence relative to acceptance limits', xlabel: 'gate number', ylabel: 'value / limit', lines: [{ name: 'value/limit', x: nums.map(r => r[0]), y: nums.map(r => Math.min(r[4], 99)) }, { name: 'acceptance (1.0)', x: [nums[0][0], nums.at(-1)[0]], y: [1, 1] }], explain: 'Each numerical criterion divided by its acceptance limit. Points above the line at 1 fail their gate; the distance below the line is the available margin.' },
    ],
    tables: [
      { key: 'gates', title: 'Qualification gates', columns: ['Gate', 'Name', 'Own evidence', 'Status', 'Key figure', 'Value', 'Limit'], rows: ev.map((g, i) => [i + 1, GATES[i], g.ok ? 'sufficient' : 'insufficient', status[i], g.num ? g.num[0] : '–', g.num ? F(g.num[1], 4) : '–', g.num ? F(g.num[2], 4) : '–']), explain: 'Evidence and resulting status of each gate of the standard qualification sequence, with the quantitative figure that decides it where one exists.' },
      { key: 'miss', title: 'Evidence still missing', columns: ['Gate', 'Name', 'Missing evidence'], rows: missing.length ? missing : [['–', '–', 'none']], explain: 'Itemised list of the evidence that must be produced, in gate order, before the next credibility class can be claimed.' },
      { key: 'cls', title: 'Credibility classes', columns: ['Class', 'Requires gates', 'Reached'], rows: CLASSES.map((c, i) => [c, ['–', '1–4', '1–5', '1–6', '1–8'][i], i === cls ? 'current' : i < cls ? 'yes' : 'no']), explain: 'The five evidence-backed classes and the gates each requires. The current class is the highest whose gates are all passed in sequence.' },
    ],
    equations: [
      { kind: 'formula', label: 'Gate sequence rule', text: 'gate k passed ⇔ evidence(k) sufficient ∧ gates 1…k−1 passed;   level = max{k : gates 1…k passed}' },
      { kind: 'formula', label: 'Classification', text: 'level < 4 → not verified;  4 → numerically verified;  5 → calibrated;  6–7 → experimentally validated over a documented domain;  8 → qualified for the application' },
      { kind: 'formula', label: 'Quantitative criteria', text: 'GCI ≤ GCI_tol,  |p_obs − p_formal| ≤ tol,  SE(θ)/|θ| ≤ tol,  RMSE_val ≤ RMSE_acc,  extrapolation ≤ limit' },
    ],
    checks: ev.map((g, i) => ({ name: `Gate ${i + 1}: ${GATES[i]}`, pass: g.ok && i < level, detail: g.ok ? (i < level ? 'evidence sufficient and all earlier gates passed' : 'evidence sufficient but blocked by an earlier gate') : g.miss.join(' ') })),
    warnings: w, recommendations: recs,
    publish: { 'cred.class': CLASSES[cls], 'cred.level': level, 'cred.gates_passed': level, 'cred.gci': inp.gciPct / 100, 'cred.rmse': inp.valRmse / 100 },
    exports: [{ name: 'credibility_record.json', mime: 'application/json', content: JSON.stringify({ class: CLASSES[cls], consecutiveGatesPassed: level, gates: ev.map((g, i) => ({ gate: i + 1, name: GATES[i], ownEvidenceSufficient: g.ok, status: status[i], keyFigure: g.num ? { name: g.num[0], value: g.num[1], limit: g.num[2] } : null, missing: g.miss })) }, null, 1), description: 'Machine-readable credibility record with gate status, key figures and missing evidence.' }],
    numerics: { method: 'Rule-based sequential gate evaluation of declared evidence', dof: 8, steps: 1, iterations: 1, residual: 0, notes: 'The assessment evaluates the evidence entered; it does not itself verify that the underlying records exist.' },
  };
}

export default {
  id: 24,
  tagline: 'Quantifies numerical error, calibrates and validates models against data, propagates uncertainty and assigns an evidence-based credibility class.',
  inputs,
  presets: [
    { name: 'ZnO nanorod growth rate: Arrhenius calibration', note: 'Axial growth rate of hydrothermal ZnO nanorods versus temperature (60–95 °C): shows the classic ln A–Eₐ correlation over a narrow temperature window.', values: { mode: 'calib', cmodel: 'arrhenius', cdata: CAL_DEF, rxn: 'Zn(NO3)2 + NaOH -> ZnO + NaNO3 + H2O' } },
    { name: 'Au nanoparticle colloid: DLVO barrier uncertainty', note: 'Propagates uncertainty of Hamaker constant (2.5×10⁻¹⁹ J, typical of gold across water), surface potential, radius and ionic strength of a 50 nm gold colloid to the aggregation barrier.', values: { mode: 'uq', umodel: 'dlvo', uqOverride: 'A_H=2.5e-19,20', dist: 'lognormal', limRatio: 0.5, limSide: 'below', nSamp: 2000 } },
    { name: 'Auxetic lattice FE model: three-mesh GCI', note: 'Effective modulus (MPa) of a re-entrant auxetic lattice from three beam-element meshes with refinement ratio 2; checks mesh independence before homogenisation.', values: { mode: 'gci', gh1: 0.25, gh2: 0.5, gh3: 1, gf1: 41.82, gf2: 42.31, gf3: 44.05, gOrder: 2, gTol: 2 } },
    { name: 'Metasurface split-ring resonance: tolerance analysis', note: 'Fabrication tolerances of a printed split-ring resonator propagated to its resonance frequency with Sobol ranking of the dimensions.', values: { mode: 'uq', umodel: 'srr', dist: 'normal', limRatio: 1.03, limSide: 'above', nSamp: 3000 } },
    { name: 'TiO₂ anatase crystallisation: JMAK calibration', note: 'Transformed fraction versus annealing time fitted with the JMAK law; rate constant and Avrami exponent with Bayesian credible intervals.', values: { mode: 'calib', cmodel: 'jmak', cdata: JMAK_DEF, rxn: 'Ti(OH)4 -> TiO2 + H2O' } },
    { name: 'Dye adsorption on nanoporous carbon: Langmuir fit', note: 'Adsorption isotherm of a nanoporous sorbent; capacity and affinity constant with identifiability diagnostics.', values: { mode: 'calib', cmodel: 'langmuir', cdata: LANG_DEF, rxn: '' } },
    { name: 'Credibility record of a growth model', note: 'Verified and calibrated model without independent validation: illustrates why the class stops at "calibrated".', values: { mode: 'cred', nVal: 0, g7uq: true, g7sens: true } },
  ],
  run(inp, ctx) {
    const r = ({ gci: runGci, mms: runMms, calib: runCalib, valid: runValid, uq: runUq, cred: runCred }[inp.mode] || runCalib)(inp, ctx);
    // guard against non-finite values reaching the interface
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    r.warnings.push('Reduced-order statistical and numerical procedures: results quantify evidence for the stated data, model and criteria and do not by themselves qualify a model for an application.');
    return r;
  },
};
