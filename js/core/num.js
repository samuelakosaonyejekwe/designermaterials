// Shared numerical kernel: linear algebra, ODE integration, root finding, optimisation,
// statistics, sampling and spectral tools used by every engineering suite.

export const PI = Math.PI;
export const K_B = 1.380649e-23, N_A = 6.02214076e23, R_GAS = 8.314462618, Q_E = 1.602176634e-19,
  EPS0 = 8.8541878128e-12, MU0 = 1.25663706212e-6, C0 = 299792458, H_PLANCK = 6.62607015e-34,
  HBAR = 1.054571817e-34, FARADAY = 96485.33212, SIGMA_SB = 5.670374419e-8, M_E = 9.1093837015e-31, G0 = 9.80665;

// ---------- basics ----------
export const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => n === 1 ? a : a + (b - a) * i / (n - 1));
export const logspace = (a, b, n) => linspace(Math.log10(a), Math.log10(b), n).map(v => 10 ** v);
export const sum = a => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
export const mean = a => a.length ? sum(a) / a.length : NaN;
export const variance = a => { const m = mean(a); let s = 0; for (const v of a) s += (v - m) ** 2; return a.length > 1 ? s / (a.length - 1) : 0; };
export const std = a => Math.sqrt(variance(a));
export const minOf = a => { let m = Infinity; for (const v of a) if (v < m) m = v; return m; };
export const maxOf = a => { let m = -Infinity; for (const v of a) if (v > m) m = v; return m; };
export const argmax = a => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[k]) k = i; return k; };
export const argmin = a => { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] < a[k]) k = i; return k; };
export function quantile(a, q) { const s = [...a].sort((x, y) => x - y); const p = (s.length - 1) * q, i = Math.floor(p); return s[i] + (s[Math.min(i + 1, s.length - 1)] - s[i]) * (p - i); }
export function histogram(a, bins = 20, lo = minOf(a), hi = maxOf(a)) {
  if (hi === lo) hi = lo + 1; const w = (hi - lo) / bins, c = new Array(bins).fill(0);
  for (const v of a) c[clamp(Math.floor((v - lo) / w), 0, bins - 1)]++;
  return { x: c.map((_, i) => lo + (i + 0.5) * w), y: c, width: w };
}
export function interp1(xs, ys, x) {
  if (x <= xs[0]) return ys[0]; const n = xs.length; if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}
export function trapz(xs, ys) { let s = 0; for (let i = 1; i < xs.length; i++) s += 0.5 * (ys[i] + ys[i - 1]) * (xs[i] - xs[i - 1]); return s; }
export function fmt(v, sig = 4) {
  if (v === null || v === undefined || Number.isNaN(v)) return '—';
  if (typeof v !== 'number') return String(v);
  if (!Number.isFinite(v)) return v > 0 ? '∞' : '−∞';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) { const [m, e] = v.toExponential(sig - 1).split('e'); return `${+m}×10${sup(String(+e))}`; }
  return String(+v.toPrecision(sig));
}
const SUP = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '-': '⁻', '+': '⁺' };
export const sup = s => [...s].map(c => SUP[c] ?? c).join('');

// ---------- random ----------
export function rng(seed = 12345) {
  let a = seed >>> 0;
  const u = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  u.normal = (m = 0, s = 1) => { let x = 0, y = 0; while (x === 0) x = u(); y = u(); return m + s * Math.sqrt(-2 * Math.log(x)) * Math.cos(2 * PI * y); };
  u.lognormal = (mu, sg) => Math.exp(u.normal(mu, sg));
  u.int = n => Math.floor(u() * n);
  return u;
}
export function lhs(n, d, rand = rng(1)) {
  const out = Array.from({ length: n }, () => new Array(d));
  for (let j = 0; j < d; j++) { const p = Array.from({ length: n }, (_, i) => i); for (let i = n - 1; i > 0; i--) { const k = rand.int(i + 1);[p[i], p[k]] = [p[k], p[i]]; } for (let i = 0; i < n; i++) out[i][j] = (p[i] + rand()) / n; }
  return out;
}

// ---------- special functions ----------
export function erf(x) { const s = Math.sign(x); x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x); const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return s * y; }
export const normCdf = (x, m = 0, s = 1) => 0.5 * (1 + erf((x - m) / (s * Math.SQRT2)));
export function normInv(p) { // Acklam
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239], b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572], c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783], d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  if (p <= 0) return -Infinity; if (p >= 1) return Infinity; let q, r;
  if (p < 0.02425) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p > 1 - 0.02425) { q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  q = p - 0.5; r = q * q; return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
export function gammaLn(x) { const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.001208650973866179, -0.000005395239384953]; let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t); let s = 1.000000000190015; for (const c of g) s += c / ++y; return -t + Math.log(2.5066282746310005 * s / x); }
export const gamma = x => Math.exp(gammaLn(x));
export function besselJ(n, x) { // integer order, integral representation (Simpson)
  const N = 200; let s = 0; for (let i = 0; i <= N; i++) { const t = PI * i / N, w = i === 0 || i === N ? 1 : i % 2 ? 4 : 2; s += w * Math.cos(n * t - x * Math.sin(t)); } return s / (3 * N);
}

// ---------- complex ----------
export const C = {
  of: (re, im = 0) => ({ re, im }),
  add: (a, b) => ({ re: a.re + b.re, im: a.im + b.im }),
  sub: (a, b) => ({ re: a.re - b.re, im: a.im - b.im }),
  mul: (a, b) => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }),
  div: (a, b) => { const d = b.re * b.re + b.im * b.im; return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d }; },
  scale: (a, s) => ({ re: a.re * s, im: a.im * s }),
  abs: a => Math.hypot(a.re, a.im), abs2: a => a.re * a.re + a.im * a.im, arg: a => Math.atan2(a.im, a.re),
  conj: a => ({ re: a.re, im: -a.im }),
  exp: a => { const e = Math.exp(a.re); return { re: e * Math.cos(a.im), im: e * Math.sin(a.im) }; },
  sqrt: a => { const r = Math.hypot(a.re, a.im); const re = Math.sqrt((r + a.re) / 2), im = Math.sign(a.im || 1) * Math.sqrt((r - a.re) / 2); return { re, im }; },
  cos: a => ({ re: Math.cos(a.re) * Math.cosh(a.im), im: -Math.sin(a.re) * Math.sinh(a.im) }),
  sin: a => ({ re: Math.sin(a.re) * Math.cosh(a.im), im: Math.cos(a.re) * Math.sinh(a.im) }),
  inv: a => C.div({ re: 1, im: 0 }, a),
};

// ---------- linear algebra (dense, arrays of rows) ----------
export const zeros = (n, m) => m === undefined ? new Array(n).fill(0) : Array.from({ length: n }, () => new Array(m).fill(0));
export const eye = n => { const I = zeros(n, n); for (let i = 0; i < n; i++) I[i][i] = 1; return I; };
export const matVec = (A, x) => A.map(r => { let s = 0; for (let j = 0; j < r.length; j++) s += r[j] * x[j]; return s; });
export const matMul = (A, B) => A.map(r => B[0].map((_, j) => { let s = 0; for (let k = 0; k < B.length; k++) s += r[k] * B[k][j]; return s; }));
export const transpose = A => A[0].map((_, j) => A.map(r => r[j]));
export const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
export const norm = a => Math.sqrt(dot(a, a));
export function solveLinear(Ain, bin) {
  const n = bin.length, A = Ain.map(r => [...r]), b = [...bin];
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i;
    if (Math.abs(A[p][k]) < 1e-300) throw new Error('Singular matrix');
    if (p !== k) { [A[p], A[k]] = [A[k], A[p]];[b[p], b[k]] = [b[k], b[p]]; }
    for (let i = k + 1; i < n; i++) { const f = A[i][k] / A[k][k]; if (f === 0) continue; for (let j = k; j < n; j++) A[i][j] -= f * A[k][j]; b[i] -= f * b[k]; }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = s / A[i][i]; }
  return x;
}
export function inverse(A) { const n = A.length, I = eye(n); return transpose(I.map(e => solveLinear(A, e))); }
export function det(Ain) { const n = Ain.length, A = Ain.map(r => [...r]); let d = 1; for (let k = 0; k < n; k++) { let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i; if (A[p][k] === 0) return 0; if (p !== k) { [A[p], A[k]] = [A[k], A[p]]; d = -d; } d *= A[k][k]; for (let i = k + 1; i < n; i++) { const f = A[i][k] / A[k][k]; for (let j = k; j < n; j++) A[i][j] -= f * A[k][j]; } } return d; }
// Thomas algorithm: a sub, b diag, c super, d rhs
export function solveTridiag(a, b, c, d) {
  const n = d.length, cp = new Float64Array(n), dp = new Float64Array(n), x = new Float64Array(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) { const m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (d[i] - a[i] * dp[i - 1]) / m; }
  x[n - 1] = dp[n - 1]; for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}
// Conjugate gradient with Jacobi preconditioner; matvec(x) -> Ax, diag optional
export function cg(matvec, b, { x0, tol = 1e-8, maxIter = 2000, diag } = {}) {
  const n = b.length, x = x0 ? Float64Array.from(x0) : new Float64Array(n); let r = new Float64Array(n); const Ax = matvec(x);
  for (let i = 0; i < n; i++) r[i] = b[i] - Ax[i];
  const z = new Float64Array(n), p = new Float64Array(n); const pre = () => { for (let i = 0; i < n; i++) z[i] = diag ? r[i] / diag[i] : r[i]; };
  pre(); p.set(z); let rz = dot(r, z); const b0 = norm(b) || 1; let it = 0, res = norm(r) / b0; const hist = [res];
  while (res > tol && it < maxIter) {
    const Ap = matvec(p), al = rz / dot(p, Ap);
    for (let i = 0; i < n; i++) { x[i] += al * p[i]; r[i] -= al * Ap[i]; }
    pre(); const rz2 = dot(r, z), be = rz2 / rz; rz = rz2; for (let i = 0; i < n; i++) p[i] = z[i] + be * p[i];
    res = norm(r) / b0; hist.push(res); it++;
  }
  return { x, iterations: it, residual: res, history: hist, converged: res <= tol };
}
// Symmetric eigenproblem (cyclic Jacobi). Returns ascending values and column vectors.
export function eigSym(Ain, maxSweep = 80) {
  const n = Ain.length, A = Ain.map(r => [...r]), V = eye(n);
  for (let sw = 0; sw < maxSweep; sw++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] ** 2; if (off < 1e-22) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(A[p][q]) < 1e-300) continue;
      const th = (A[q][q] - A[p][p]) / (2 * A[p][q]), t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
      for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
    }
  }
  const idx = A.map((_, i) => i).sort((i, j) => A[i][i] - A[j][j]);
  return { values: idx.map(i => A[i][i]), vectors: idx.map(i => V.map(r => r[i])) };
}
export function cholesky(A) { const n = A.length, L = zeros(n, n); for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) { let s = A[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k]; if (i === j) { if (s <= 0) throw new Error('Matrix not positive definite'); L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j]; } return L; }
// Generalised symmetric-definite K v = w M v
export function eigGen(K, M) {
  const n = K.length, L = cholesky(M), Li = inverse(L), A = matMul(matMul(Li, K), transpose(Li));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) A[i][j] = A[j][i] = 0.5 * (A[i][j] + A[j][i]);
  const e = eigSym(A), LiT = transpose(Li); return { values: e.values, vectors: e.vectors.map(v => matVec(LiT, v)) };
}
// Hermitian eigenvalues for complex H = Hr + i Hi (via real 2n x 2n embedding); returns n ascending values
export function eigHermitian(Hr, Hi) {
  const n = Hr.length, B = zeros(2 * n, 2 * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { B[i][j] = Hr[i][j]; B[i + n][j + n] = Hr[i][j]; B[i][j + n] = -Hi[i][j]; B[i + n][j] = Hi[i][j]; }
  const v = eigSym(B).values; return v.filter((_, i) => i % 2 === 0);
}
// Least squares via normal equations with Tikhonov option
export function lstsq(A, b, lambda = 0) { const At = transpose(A), N = matMul(At, A); for (let i = 0; i < N.length; i++) N[i][i] += lambda; return solveLinear(N, matVec(At, b)); }
export function polyfit(x, y, deg) { return lstsq(x.map(v => Array.from({ length: deg + 1 }, (_, k) => v ** k)), y); }
export const polyval = (c, x) => c.reduceRight((s, ck) => s * x + ck, 0);

// ---------- ODE ----------
export function rk4(f, y0, t0, t1, n) {
  const h = (t1 - t0) / n, t = [t0], y = [[...y0]]; let yc = [...y0], tc = t0; const m = y0.length;
  for (let i = 0; i < n; i++) {
    const k1 = f(tc, yc), k2 = f(tc + h / 2, yc.map((v, j) => v + h / 2 * k1[j])), k3 = f(tc + h / 2, yc.map((v, j) => v + h / 2 * k2[j])), k4 = f(tc + h, yc.map((v, j) => v + h * k3[j]));
    yc = yc.map((v, j) => v + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); tc += h; t.push(tc); y.push([...yc]);
    if (m && !yc.every(Number.isFinite)) break;
  }
  return { t, y };
}
// Adaptive Dormand–Prince 5(4)
export function rk45(f, y0, t0, t1, { rtol = 1e-6, atol = 1e-9, hmax = Infinity, maxSteps = 20000 } = {}) {
  const c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1], a = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9], [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729], [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656], [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]], b5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0], b4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];
  const n = y0.length; let t = t0, y = [...y0], h = Math.min((t1 - t0) / 100, hmax); const T = [t], Y = [[...y]]; let steps = 0, rejected = 0;
  while (t < t1 && steps < maxSteps) {
    if (t + h > t1) h = t1 - t; const k = [];
    for (let s = 0; s < 7; s++) { const ys = y.map((v, j) => { let acc = v; for (let q = 0; q < s; q++) acc += h * a[s][q] * k[q][j]; return acc; }); k.push(f(t + c[s] * h, ys)); }
    let err = 0; const yn = new Array(n);
    for (let j = 0; j < n; j++) { let s5 = 0, s4 = 0; for (let s = 0; s < 7; s++) { s5 += b5[s] * k[s][j]; s4 += b4[s] * k[s][j]; } yn[j] = y[j] + h * s5; const sc = atol + rtol * Math.max(Math.abs(y[j]), Math.abs(yn[j])); err = Math.max(err, Math.abs(h * (s5 - s4)) / sc); }
    if (err <= 1 || h < 1e-14 * Math.max(1, Math.abs(t))) { t += h; y = yn; T.push(t); Y.push([...y]); } else rejected++;
    h = Math.min(hmax, h * clamp(0.9 * Math.pow(Math.max(err, 1e-10), -0.2), 0.2, 5)); steps++;
  }
  return { t: T, y: Y, steps, rejected, completed: t >= t1 };
}
// Implicit backward Euler with Newton iteration and finite-difference Jacobian (stiff systems)
export function bdf1(f, y0, t0, t1, n, { tol = 1e-9, maxNewton = 25, nonNegative = false } = {}) {
  const h = (t1 - t0) / n, m = y0.length, T = [t0], Y = [[...y0]]; let y = [...y0], t = t0, newtonTotal = 0;
  for (let s = 0; s < n; s++) {
    let yn = [...y]; const tn = t + h;
    for (let it = 0; it < maxNewton; it++) {
      const fy = f(tn, yn), G = yn.map((v, i) => v - y[i] - h * fy[i]), J = zeros(m, m);
      for (let j = 0; j < m; j++) { const e = Math.max(1e-8, Math.abs(yn[j]) * 1e-7), yp = [...yn]; yp[j] += e; const fp = f(tn, yp); for (let i = 0; i < m; i++) J[i][j] = (i === j ? 1 : 0) - h * (fp[i] - fy[i]) / e; }
      const dx = solveLinear(J, G.map(v => -v)); let dn = 0;
      for (let i = 0; i < m; i++) { yn[i] += dx[i]; if (nonNegative && yn[i] < 0) yn[i] = 0; dn = Math.max(dn, Math.abs(dx[i]) / (Math.abs(yn[i]) + 1e-12)); }
      newtonTotal++; if (dn < tol) break;
    }
    y = yn; t = tn; T.push(t); Y.push([...y]);
  }
  return { t: T, y: Y, newtonIterations: newtonTotal };
}

// ---------- roots ----------
export function brent(f, a, b, tol = 1e-12, maxIter = 200) {
  let fa = f(a), fb = f(b); if (fa * fb > 0) throw new Error('Root not bracketed');
  if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a];[fa, fb] = [fb, fa]; } let c = a, fc = fa, d = 0, flag = true;
  for (let i = 0; i < maxIter && fb !== 0 && Math.abs(b - a) > tol; i++) {
    let s = fa !== fc && fb !== fc ? a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb)) : b - fb * (b - a) / (fb - fa);
    const cond = (s - (3 * a + b) / 4) * (s - b) > 0 || (flag && Math.abs(s - b) >= Math.abs(b - c) / 2) || (!flag && Math.abs(s - b) >= Math.abs(c - d) / 2);
    if (cond) { s = (a + b) / 2; flag = true; } else flag = false;
    const fs = f(s); d = c; c = b; fc = fb; if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
    if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a];[fa, fb] = [fb, fa]; }
  }
  return b;
}
export function bisect(f, a, b, tol = 1e-10) { let fa = f(a); for (let i = 0; i < 200 && b - a > tol; i++) { const m = 0.5 * (a + b), fm = f(m); if (fa * fm <= 0) b = m; else { a = m; fa = fm; } } return 0.5 * (a + b); }
export function newtonSystem(F, x0, { tol = 1e-10, maxIter = 60, damping = 1 } = {}) {
  let x = [...x0]; const n = x.length, hist = [];
  for (let it = 0; it < maxIter; it++) {
    const f = F(x), r = norm(f); hist.push(r); if (r < tol) return { x, iterations: it, residual: r, history: hist, converged: true };
    const J = zeros(n, n); for (let j = 0; j < n; j++) { const e = Math.max(1e-8, Math.abs(x[j]) * 1e-7), xp = [...x]; xp[j] += e; const fp = F(xp); for (let i = 0; i < n; i++) J[i][j] = (fp[i] - f[i]) / e; }
    let dx; try { dx = solveLinear(J, f.map(v => -v)); } catch { break; }
    let lam = damping; for (let ls = 0; ls < 12; ls++) { const xt = x.map((v, i) => v + lam * dx[i]); if (norm(F(xt)) < r || ls === 11) { x = xt; break; } lam /= 2; }
  }
  const r = norm(F(x)); return { x, iterations: maxIter, residual: r, history: hist, converged: r < tol };
}

// ---------- optimisation ----------
export function goldenSection(f, a, b, tol = 1e-8) { const g = (Math.sqrt(5) - 1) / 2; let c = b - g * (b - a), d = a + g * (b - a); while (Math.abs(b - a) > tol) { if (f(c) < f(d)) b = d; else a = c; c = b - g * (b - a); d = a + g * (b - a); } return (a + b) / 2; }
export function nelderMead(f, x0, { maxIter = 600, tol = 1e-10, step = 0.1, bounds } = {}) {
  const n = x0.length, proj = x => bounds ? x.map((v, i) => clamp(v, bounds[i][0], bounds[i][1])) : x, F = x => f(proj(x));
  let S = [x0, ...x0.map((_, i) => x0.map((v, j) => j === i ? v + (Math.abs(v) > 1e-12 ? step * Math.abs(v) : step) : v))].map(x => ({ x, f: F(x) })); const hist = []; let it = 0;
  for (; it < maxIter; it++) {
    S.sort((a, b) => a.f - b.f); hist.push(S[0].f); if (Math.abs(S[n].f - S[0].f) <= tol * (Math.abs(S[0].f) + tol)) break;
    const c = zeros(n); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += S[i].x[j] / n;
    const pt = k => c.map((v, j) => v + k * (S[n].x[j] - v)); const xr = pt(-1), fr = F(xr);
    if (fr < S[0].f) { const xe = pt(-2), fe = F(xe); S[n] = fe < fr ? { x: xe, f: fe } : { x: xr, f: fr }; }
    else if (fr < S[n - 1].f) S[n] = { x: xr, f: fr };
    else { const xc = pt(fr < S[n].f ? -0.5 : 0.5), fc = F(xc); if (fc < Math.min(fr, S[n].f)) S[n] = { x: xc, f: fc }; else for (let i = 1; i <= n; i++) { S[i].x = S[i].x.map((v, j) => S[0].x[j] + 0.5 * (v - S[0].x[j])); S[i].f = F(S[i].x); } }
  }
  S.sort((a, b) => a.f - b.f); return { x: proj(S[0].x), f: S[0].f, iterations: it, history: hist };
}
export function differentialEvolution(f, bounds, { pop = 24, gens = 80, F = 0.7, CR = 0.9, seed = 7 } = {}) {
  const r = rng(seed), d = bounds.length, P = Array.from({ length: pop }, () => bounds.map(([a, b]) => a + r() * (b - a))), fit = P.map(f), hist = [];
  for (let g = 0; g < gens; g++) {
    for (let i = 0; i < pop; i++) {
      let a, b, c; do a = r.int(pop); while (a === i); do b = r.int(pop); while (b === i || b === a); do c = r.int(pop); while (c === i || c === a || c === b);
      const jr = r.int(d), u = P[i].map((v, j) => (r() < CR || j === jr) ? clamp(P[a][j] + F * (P[b][j] - P[c][j]), bounds[j][0], bounds[j][1]) : v), fu = f(u);
      if (fu <= fit[i]) { P[i] = u; fit[i] = fu; }
    }
    hist.push(minOf(fit));
  }
  const k = argmin(fit); return { x: P[k], f: fit[k], history: hist, population: P, fitness: fit };
}
export function particleSwarm(f, bounds, { swarm = 24, iters = 80, w = 0.72, c1 = 1.49, c2 = 1.49, seed = 3 } = {}) {
  const r = rng(seed), d = bounds.length, X = Array.from({ length: swarm }, () => bounds.map(([a, b]) => a + r() * (b - a))), V = X.map(x => x.map(() => 0)), pb = X.map(x => [...x]), pf = X.map(f); let g = argmin(pf), gx = [...pb[g]], gf = pf[g]; const hist = [];
  for (let it = 0; it < iters; it++) { for (let i = 0; i < swarm; i++) { for (let j = 0; j < d; j++) { V[i][j] = w * V[i][j] + c1 * r() * (pb[i][j] - X[i][j]) + c2 * r() * (gx[j] - X[i][j]); X[i][j] = clamp(X[i][j] + V[i][j], bounds[j][0], bounds[j][1]); } const fx = f(X[i]); if (fx < pf[i]) { pf[i] = fx; pb[i] = [...X[i]]; if (fx < gf) { gf = fx; gx = [...X[i]]; } } } hist.push(gf); }
  return { x: gx, f: gf, history: hist };
}
export const dominates = (a, b) => a.every((v, i) => v <= b[i]) && a.some((v, i) => v < b[i]);
export function paretoFront(objs) { return objs.map((o, i) => i).filter(i => !objs.some((p, j) => j !== i && dominates(p, objs[i]))); }
// Compact NSGA-II (minimisation of all objectives)
export function nsga2(fobj, bounds, { pop = 40, gens = 40, seed = 11, pm } = {}) {
  const r = rng(seed), d = bounds.length; pm = pm ?? 1 / d;
  let P = Array.from({ length: pop }, () => bounds.map(([a, b]) => a + r() * (b - a))).map(x => ({ x, o: fobj(x) }));
  const rank = S => { const fronts = []; let rem = S.map((_, i) => i); while (rem.length) { const fr = rem.filter(i => !rem.some(j => j !== i && dominates(S[j].o, S[i].o))); fronts.push(fr); rem = rem.filter(i => !fr.includes(i)); } fronts.forEach((fr, k) => { const m = S[0].o.length; fr.forEach(i => { S[i].rank = k; S[i].cd = 0; }); for (let q = 0; q < m; q++) { const s = [...fr].sort((a, b) => S[a].o[q] - S[b].o[q]), span = S[s.at(-1)].o[q] - S[s[0]].o[q] || 1; S[s[0]].cd = S[s.at(-1)].cd = Infinity; for (let i = 1; i < s.length - 1; i++) S[s[i]].cd += (S[s[i + 1]].o[q] - S[s[i - 1]].o[q]) / span; } }); };
  rank(P); const pick = () => { const a = P[r.int(pop)], b = P[r.int(pop)]; return a.rank < b.rank || (a.rank === b.rank && a.cd > b.cd) ? a : b; };
  for (let g = 0; g < gens; g++) {
    const Q = []; while (Q.length < pop) { const a = pick().x, b = pick().x, u = r(); const x = a.map((v, j) => { let c = r() < 0.9 ? u * v + (1 - u) * b[j] : v; if (r() < pm) c += r.normal(0, 0.1 * (bounds[j][1] - bounds[j][0])); return clamp(c, bounds[j][0], bounds[j][1]); }); Q.push({ x, o: fobj(x) }); }
    const U = P.concat(Q); rank(U); U.sort((a, b) => a.rank - b.rank || b.cd - a.cd); P = U.slice(0, pop);
  }
  rank(P); return { population: P, front: P.filter(p => p.rank === 0) };
}
// Levenberg–Marquardt: minimise sum (model(p)[i]-y[i])^2 ; resid(p) returns residual vector
export function levenbergMarquardt(resid, p0, { maxIter = 100, tol = 1e-10, bounds } = {}) {
  let p = [...p0], lam = 1e-3; const proj = q => bounds ? q.map((v, i) => clamp(v, bounds[i][0], bounds[i][1])) : q; let r = resid(p), cost = dot(r, r); const hist = [cost]; let J = null, it = 0;
  for (; it < maxIter; it++) {
    const m = r.length, n = p.length; J = zeros(m, n);
    for (let j = 0; j < n; j++) { const e = Math.max(1e-8, Math.abs(p[j]) * 1e-6), q = [...p]; q[j] += e; const rq = resid(q); for (let i = 0; i < m; i++) J[i][j] = (rq[i] - r[i]) / e; }
    const Jt = transpose(J), H = matMul(Jt, J), g = matVec(Jt, r); let improved = false;
    for (let k = 0; k < 12; k++) {
      const A = H.map((row, i) => row.map((v, j) => i === j ? v * (1 + lam) + 1e-300 : v)); let dp; try { dp = solveLinear(A, g.map(v => -v)); } catch { lam *= 10; continue; }
      const q = proj(p.map((v, i) => v + dp[i])), rq = resid(q), cq = dot(rq, rq);
      if (cq < cost) { const rel = (cost - cq) / (cost + 1e-300); p = q; r = rq; cost = cq; lam = Math.max(lam / 10, 1e-12); improved = true; hist.push(cost); if (rel < tol) k = 99; break; } lam *= 10;
    }
    if (!improved || (hist.length > 1 && Math.abs(hist.at(-2) - hist.at(-1)) <= tol * (cost + tol))) break;
  }
  // covariance estimate
  let cov = null, se = null; try { const m = r.length, n = p.length, s2 = cost / Math.max(1, m - n); cov = inverse(matMul(transpose(J), J)).map(row => row.map(v => v * s2)); se = cov.map((row, i) => Math.sqrt(Math.abs(row[i]))); } catch { /* singular: parameters not identifiable */ }
  return { p, cost, iterations: it, history: hist, covariance: cov, stderr: se, residuals: r };
}

// ---------- statistics ----------
export function linreg(x, y) { const n = x.length, mx = mean(x), my = mean(y); let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; } const slope = sxy / sxx, intercept = my - slope * mx; return { slope, intercept, r2: syy ? sxy * sxy / (sxx * syy) : 1 }; }
export const rmse = (a, b) => Math.sqrt(mean(a.map((v, i) => (v - b[i]) ** 2)));
export const mae = (a, b) => mean(a.map((v, i) => Math.abs(v - b[i])));
export function r2(pred, obs) { const m = mean(obs); let ss = 0, st = 0; for (let i = 0; i < obs.length; i++) { ss += (obs[i] - pred[i]) ** 2; st += (obs[i] - m) ** 2; } return st ? 1 - ss / st : 1; }
export function correlation(a, b) { const ma = mean(a), mb = mean(b); let s = 0, sa = 0, sb = 0; for (let i = 0; i < a.length; i++) { s += (a[i] - ma) * (b[i] - mb); sa += (a[i] - ma) ** 2; sb += (b[i] - mb) ** 2; } return sa && sb ? s / Math.sqrt(sa * sb) : 0; }
export function weibullFit(data) { // median-rank regression
  const s = [...data].filter(v => v > 0).sort((a, b) => a - b), n = s.length, x = s.map(Math.log), y = s.map((_, i) => Math.log(-Math.log(1 - (i + 0.7) / (n + 0.4)))); const f = linreg(x, y); return { shape: f.slope, scale: Math.exp(-f.intercept / f.slope), r2: f.r2 };
}
// Richardson extrapolation / grid convergence index for three solutions on grids with refinement ratio r (fine->coarse f1,f2,f3)
export function gci(f1, f2, f3, r = 2, Fs = 1.25) { const p = Math.log(Math.abs((f3 - f2) / (f2 - f1))) / Math.log(r), fx = f1 + (f1 - f2) / (r ** p - 1), g = Fs * Math.abs((f2 - f1) / f1) / (r ** p - 1); return { order: p, extrapolated: fx, gciFine: g }; }

// ---------- spectral ----------
export function fft(re, im, inverse = false) {
  const n = re.length; if (n & (n - 1)) throw new Error('FFT length must be a power of two');
  for (let i = 1, j = 0; i < n; i++) { let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit; if (i < j) { [re[i], re[j]] = [re[j], re[i]];[im[i], im[j]] = [im[j], im[i]]; } }
  for (let len = 2; len <= n; len <<= 1) { const ang = 2 * PI / len * (inverse ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang); for (let i = 0; i < n; i += len) { let cr = 1, ci = 0; for (let k = 0; k < len / 2; k++) { const a = i + k, b = a + len / 2, tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr; re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti; const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t; } } }
  if (inverse) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  return { re, im };
}
export function findPeaks(y, { minProminence = 0, minDistance = 1 } = {}) {
  const pk = []; for (let i = 1; i < y.length - 1; i++) if (y[i] > y[i - 1] && y[i] >= y[i + 1]) { let l = i, r = i, lm = y[i], rm = y[i]; while (l > 0 && y[l - 1] <= y[i]) { l--; lm = Math.min(lm, y[l]); } while (r < y.length - 1 && y[r + 1] <= y[i]) { r++; rm = Math.min(rm, y[r]); } const prom = y[i] - Math.max(lm, rm); if (prom >= minProminence) pk.push({ i, y: y[i], prominence: prom }); }
  pk.sort((a, b) => b.y - a.y); const keep = []; for (const p of pk) if (keep.every(q => Math.abs(q.i - p.i) >= minDistance)) keep.push(p); return keep.sort((a, b) => a.i - b.i);
}
