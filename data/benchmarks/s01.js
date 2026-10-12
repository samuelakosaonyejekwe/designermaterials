// Reference benchmarks for suite 1 — Material Requirements and Functional Design.
// Expected values are closed-form results or handbook data evaluated here; none is taken from the suite.

// standard normal distribution function (Abramowitz & Stegun 26.2.17, |error| < 7.5e-8)
const Phi = x => { const t = 1 / (1 + 0.2316419 * Math.abs(x)), d = Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI), p = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))); return x >= 0 ? 1 - p : p; };

const STACK4 = 'contributor,nominal_mm,tolerance_mm,sensitivity\nA,10,0.1,1\nB,10,0.1,1\nC,10,0.1,1\nD,10,0.1,1';
const STACK345 = 'contributor,nominal_mm,tolerance_mm,sensitivity\nA,20,0.03,1\nB,5,0.02,2';   // sensitivity-weighted tolerances 0.03 and 0.04

export default [
  {
    id: 's01-stack-worst-case', suite: 1, title: 'Tolerance stack-up: worst-case sum of four equal contributors', mode: 'alloc',
    inputs: { mode: 'alloc', stack: STACK4, stackLimit: 0.2, kSigma: 3 }, kpi: 'wc', expected: 0.4, tol: 1e-9, kind: 'analytical',
    reference: 'Worst-case tolerance analysis: arithmetic sum of the contributor tolerances (dimensional-management textbooks)',
  },
  {
    id: 's01-stack-rss-equal', suite: 1, title: 'Tolerance stack-up: root-sum-square of four equal contributors', mode: 'alloc',
    inputs: { mode: 'alloc', stack: STACK4, stackLimit: 0.2, kSigma: 3 }, kpi: 'rss', expected: 0.1 * Math.sqrt(4), tol: 1e-9, kind: 'analytical',
    reference: 'Statistical (RSS) tolerance analysis: T = √(Σ tᵢ²) for independent, centred contributors',
    note: 'n equal contributors give T_rss = t·√n, i.e. half the worst case for n = 4.',
  },
  {
    id: 's01-stack-rss-sensitivity', suite: 1, title: 'RSS stack with sensitivities: 3-4-5 combination', mode: 'alloc',
    inputs: { mode: 'alloc', stack: STACK345, stackLimit: 0.1, kSigma: 3 }, kpi: 'rss', expected: 0.05, tol: 1e-9, kind: 'analytical',
    reference: 'First-order variance propagation: T = √(Σ (∂y/∂xᵢ · tᵢ)²)',
    note: 'Contributors 1 × 0.03 mm and 2 × 0.02 mm combine to exactly 0.05 mm.',
  },
  {
    id: 's01-stack-cp', suite: 1, title: 'Assembly capability index when the limit equals three stack standard deviations', mode: 'alloc',
    inputs: { mode: 'alloc', stack: STACK4, stackLimit: 0.2, kSigma: 3 }, kpi: 'cp', expected: 1, tol: 1e-9, kind: 'analytical',
    reference: 'Process capability index C_p = T/(3σ) (statistical process control textbooks)',
  },
  {
    id: 's01-stack-monte-carlo-1sigma', suite: 1, title: 'Monte Carlo reject fraction for a limit of one stack standard deviation', mode: 'alloc',
    inputs: { mode: 'alloc', stack: STACK4, stackLimit: 0.2 / 3, kSigma: 3, seed: 1 }, kpi: 'p_out', expected: 2 * (1 - Phi(1)), tol: 0.06, kind: 'analytical',
    reference: 'Normal distribution: P(|X| > σ) = 0.3173',
    note: 'Seeded sample of 4000 assemblies: the binomial standard error is 2.3 % of the expected value, so the tolerance is about 2.5 standard errors.',
  },
  {
    id: 's01-margin-of-safety', suite: 1, title: 'Margin of safety with derating and required factor of safety', mode: 'alloc',
    inputs: { mode: 'alloc', cap: 2.5e8, dem: 9e7, derate: 0.8, fosReq: 1.5 }, kpi: 'mos', expected: 0.8 * 2.5e8 / (1.5 * 9e7) - 1, tol: 1e-9, kind: 'analytical',
    reference: 'Margin of safety MoS = allowable/(FoS × applied) − 1 (aerospace structural design practice)',
  },
  {
    id: 's01-reliability-index', suite: 1, title: 'First-order second-moment reliability index', mode: 'alloc',
    inputs: { mode: 'alloc', cap: 2e8, dem: 1e8, covCap: 15, covDem: 40 }, kpi: 'beta', expected: 2, tol: 1e-9, kind: 'analytical',
    reference: 'Cornell reliability index β = (μ_R − μ_S)/√(σ_R² + σ_S²) (structural reliability textbooks)',
    note: 'σ_R = 30 MPa and σ_S = 40 MPa combine to 50 MPa; the 100 MPa mean margin gives β = 2 exactly.',
  },
  {
    id: 's01-failure-probability', suite: 1, title: 'Probability of failure for β = 2', mode: 'alloc',
    inputs: { mode: 'alloc', cap: 2e8, dem: 1e8, covCap: 15, covDem: 40 }, kpi: 'pf', expected: 0.0227501, tol: 1e-3, kind: 'analytical',
    reference: 'Standard normal table: Φ(−2) = 0.02275',
  },
  {
    id: 's01-series-reliability', suite: 1, title: 'Series system of ten items of reliability 0.99', mode: 'alloc',
    inputs: { mode: 'alloc', nSub: 10, Ritem: 0.99, relW: '' }, kpi: 'Rach', expected: Math.pow(0.99, 10), tol: 1e-9, kind: 'analytical',
    reference: 'Series reliability block diagram: R_sys = Π Rᵢ (reliability engineering textbooks)',
  },
  {
    id: 's01-reliability-allocation', suite: 1, title: 'Equal apportionment of a system reliability of 0.90 to five series items', mode: 'alloc',
    inputs: { mode: 'alloc', Rsys: 0.9, nSub: 5, relW: '' }, kpi: 'Ri', expected: Math.pow(0.9, 1 / 5), tol: 1e-9, kind: 'analytical',
    reference: 'Equal-apportionment reliability allocation: Rᵢ = R_sys^(1/n)',
  },
  {
    id: 's01-allocated-failure-rate', suite: 1, title: 'Constant failure rate equivalent to the allocated reliability', mode: 'alloc',
    inputs: { mode: 'alloc', Rsys: 0.9, nSub: 5, relW: '', tMission: 1000 }, kpi: 'lambda', expected: -Math.log(0.9) / 5 / 1000, tol: 1e-9, kind: 'analytical',
    reference: 'Exponential life model: R = exp(−λt)',
  },
  {
    id: 's01-pi-count-reactor', suite: 1, title: 'Buckingham Pi: seven reactor variables in M, L, T', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'reactor' }, kpi: 'npi', expected: 4, abs: 0, kind: 'analytical',
    reference: 'Buckingham Pi theorem: number of independent groups = n − rank of the dimensional matrix',
  },
  {
    id: 's01-pi-count-lattice-rank-deficient', suite: 1, title: 'Buckingham Pi: rank-deficient dimensional matrix (lattice mechanics)', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'lattice' }, kpi: 'npi', expected: 5, abs: 0, kind: 'analytical',
    reference: 'Buckingham Pi theorem: number of independent groups = n − rank of the dimensional matrix',
    note: 'Seven variables with dimensions of stress or length only: three base dimensions occur but the rank is 2, so there are 5 groups, not 4.',
  },
  {
    id: 's01-pi-count-heat', suite: 1, title: 'Buckingham Pi: six heat-up variables in M, L, T, Θ', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'heat' }, kpi: 'npi', expected: 2, abs: 0, kind: 'analytical',
    reference: 'Buckingham Pi theorem; the two groups are the Biot and Fourier numbers',
  },
  {
    id: 's01-reynolds-number', suite: 1, title: 'Pore Reynolds number of water at 1 mm/s in 100 µm pores', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'reactor', pU: 1e-3, pd: 1e-4, pRho: 998, pMu: 1e-3 }, kpi: 'pi1', expected: 998 * 1e-3 * 1e-4 / 1e-3, tol: 1e-9, kind: 'analytical',
    reference: 'Definition Re = ρUd/μ (fluid mechanics textbooks)',
  },
  {
    id: 's01-schmidt-number', suite: 1, title: 'Schmidt number of a dissolved ion in water', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'reactor', pRho: 998, pMu: 1e-3, pD: 1e-9 }, kpi: 'pi2', expected: 1e-3 / (998 * 1e-9), tol: 1e-9, kind: 'analytical',
    reference: 'Definition Sc = μ/(ρD); about 10³ for ions in water (transport-phenomena textbooks)',
  },
  {
    id: 's01-fourier-number', suite: 1, title: 'Fourier number of a ZnO body (L = 1 mm) after 10 s', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'heat', pKth: 50, pLh: 1e-3, pRhoS: 5606, pCp: 494, pT: 10 }, kpi: 'pi2', expected: 50 * 10 / (5606 * 494 * 1e-6), tol: 1e-9, kind: 'analytical',
    reference: 'Definition Fo = kt/(ρc_pL²) (heat-transfer textbooks)',
  },
  {
    id: 's01-frequency-ratio', suite: 1, title: 'Frequency ratio of a 10 g, 4 kN/m resonator driven at 120 Hz', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'resonator', pM: 0.01, pKs: 4000, pF: 120 }, kpi: 'pi1', expected: 2 * Math.PI * 120 * Math.sqrt(0.01 / 4000), tol: 1e-9, kind: 'analytical',
    reference: 'Single-degree-of-freedom oscillator: r = ω/ω₀, ω₀ = √(k/m) (vibration textbooks)',
  },
  {
    id: 's01-magnification-factor', suite: 1, title: 'Dynamic magnification of an undamped oscillator at r = √2', mode: 'pi',
    inputs: { mode: 'pi', piProb: 'resonator', pM: 0.01, pKs: 4000, pC: 0, pF: Math.sqrt(2) * Math.sqrt(4000 / 0.01) / (2 * Math.PI), pF0: 1 }, kpi: 'pi3', expected: 1, tol: 1e-9, kind: 'analytical',
    reference: 'Forced vibration: kX/F₀ = 1/|1 − r²| without damping, which equals 1 at r = √2 (vibration textbooks)',
  },
  {
    id: 's01-ahp-consistent-cr', suite: 1, title: 'AHP: consistency ratio of a perfectly consistent matrix', mode: 'ahp',
    inputs: { mode: 'ahp', ahpNames: 'a, b, c', ahpMatrix: '1, 5/3, 5/2; 3/5, 1, 3/2; 2/5, 2/3, 1' }, kpi: 'CR', expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Saaty, The Analytic Hierarchy Process: λ_max = n and CR = 0 for a consistent reciprocal matrix',
    note: 'Matrix built as aᵢⱼ = wᵢ/wⱼ from the weights 0.5, 0.3, 0.2.',
  },
  {
    id: 's01-ahp-consistent-weight', suite: 1, title: 'AHP: principal eigenvector recovers the generating weights', mode: 'ahp',
    inputs: { mode: 'ahp', ahpNames: 'a, b, c', ahpMatrix: '1, 5/3, 5/2; 3/5, 1, 3/2; 2/5, 2/3, 1' }, kpi: 'w_top', expected: 0.5, tol: 1e-9, kind: 'analytical',
    reference: 'Saaty, The Analytic Hierarchy Process: A·w = n·w for aᵢⱼ = wᵢ/wⱼ',
  },
  {
    id: 's01-ahp-3x3-eigenvalue', suite: 1, title: 'AHP: principal eigenvalue of a 3 × 3 reciprocal matrix (closed form)', mode: 'ahp',
    inputs: { mode: 'ahp', ahpNames: 'a, b, c', ahpMatrix: '1, 3, 5; 1/3, 1, 3; 1/5, 1/3, 1' }, kpi: 'lambda', expected: 1 + Math.cbrt(3 * 3 / 5) + Math.cbrt(5 / (3 * 3)), tol: 1e-9, kind: 'analytical',
    reference: 'Characteristic polynomial of a 3 × 3 positive reciprocal matrix: λ_max = 1 + (a₁₂a₂₃/a₁₃)^⅓ + (a₁₃/(a₁₂a₂₃))^⅓',
  },
  {
    id: 's01-ahp-3x3-cr', suite: 1, title: 'AHP: consistency ratio of the same 3 × 3 matrix with Saaty random index 0.58', mode: 'ahp',
    inputs: { mode: 'ahp', ahpNames: 'a, b, c', ahpMatrix: '1, 3, 5; 1/3, 1, 3; 1/5, 1/3, 1' }, kpi: 'CR', expected: (Math.cbrt(9 / 5) + Math.cbrt(5 / 9) - 2) / 2 / 0.58, tol: 1e-8, kind: 'analytical',
    reference: 'Saaty, The Analytic Hierarchy Process: CI = (λ_max − n)/(n − 1), CR = CI/RI with RI(3) = 0.58',
  },
  {
    id: 's01-rank-bandgap-match', suite: 1, title: 'Semiconductor screening: best band-gap match to 3.37 eV (ZnO) among nanomaterials', mode: 'rank',
    inputs: { mode: 'rank', func: 'semi', kind: 'nanomaterial', EgTarget: 3.37, rhoMax: 25000, Emin: 0, wPerf: 1, wRho: 0, wCost: 0, wCo2: 0, wTemp: 0 }, kpi: 'M1', expected: 1, tol: 0.02, kind: 'reference-data',
    reference: 'Room-temperature band gap of wurtzite ZnO, 3.37 eV (semiconductor data handbooks)',
    note: 'With all weight on the index the leader is the library entry closest to the target gap; the tolerance covers the 3.3–3.4 eV spread of tabulated ZnO values.',
  },
  {
    id: 's01-spec-zno-gap-window', suite: 1, title: 'Specification check: ZnO nanorods inside a 3.0–3.6 eV band-gap window', mode: 'spec',
    inputs: { mode: 'spec', func: 'semi', kind: 'nanomaterial', specId: 'zno-nanorods', EgTarget: 3.3, EgMin: 3.0, EgMax: 3.6, rhoMax: 25000, Emin: 0, Top: 150 }, kpi: 'n_fail', expected: 0, abs: 0, kind: 'reference-data',
    reference: 'Wurtzite ZnO: band gap 3.37 eV, melting point about 2248 K (semiconductor data handbooks)',
    note: 'A 150 °C service temperature is far below half the melting point, so no requirement may fail.',
  },
  {
    id: 's01-spec-requirement-count', suite: 1, title: 'Specification sheet: one functional requirement plus one per active limit', mode: 'spec',
    inputs: { mode: 'spec', func: 'semi', kind: 'nanomaterial', specId: 'zno-nanorods', EgMin: 3.0, EgMax: 3.6, rhoMax: 25000, Emin: 0, Top: 150 }, kpi: 'n_req', expected: 4, abs: 0, kind: 'conservation',
    reference: 'Requirements traceability: every stated limit maps to exactly one requirement',
    note: 'Active limits: lower and upper band gap and the service-temperature rule; density and modulus limits are switched off.',
  },
];
