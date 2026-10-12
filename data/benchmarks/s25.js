// Reference benchmarks for suite 25 — Inverse Design, Optimisation and Artificial Intelligence.
// Each inverse problem is posed so that its solution is known in closed form or from a published benchmark;
// the expected values are evaluated here and never taken from the suite.

const PI = Math.PI, R_GAS = 8.314462618, HBAR = 1.054571817e-34, M_E = 9.1093837015e-31, Q_E = 1.602176634e-19, EPS0 = 8.8541878128e-12;

// --- single quarter-wave layer and quarter-wave stack at the design wavelength (characteristic-admittance algebra)
const rQW = (n1, ns) => ((ns - n1 * n1) / (ns + n1 * n1)) ** 2;
const rStack = (nH, nL, ns, pairs) => { const Y = ns * (nH / nL) ** (2 * pairs); return ((1 - Y) / (1 + Y)) ** 2; };
// --- Brus equation for CdSe (E_g 1.74 eV, m_e 0.13, m_h 0.45, ε_r 10.6): closed-form root of the quadratic in 1/R
const mu = 1 / (1 / 0.13 + 1 / 0.45), Abrus = HBAR ** 2 * PI ** 2 / (2 * mu * M_E) / Q_E, Bbrus = 1.786 * Q_E / (4 * PI * EPS0 * 10.6), dE = 1239.84 / 580 - 1.74;
const Rbrus = (-Bbrus + Math.sqrt(Bbrus ** 2 + 4 * dE * Abrus)) / (2 * dE);
// --- Gans depolarisation factor of a prolate spheroid and its inverse
const gansL = ar => { const e = Math.sqrt(1 - 1 / (ar * ar)); return (1 - e * e) / (e * e) * (Math.log((1 + e) / (1 - e)) / (2 * e) - 1); };
const arOfL = L => { let lo = 1.001, hi = 30; for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; if (gansL(m) > L) lo = m; else hi = m; } return (lo + hi) / 2; };
// gold at 800 nm (1.550 eV), interpolated between the Johnson & Christy points 1.51 eV (n = 0.16, k = 5.083) and 1.64 eV (n = 0.14, k = 4.542)
const e151 = 0.16 ** 2 - 5.083 ** 2, e164 = 0.14 ** 2 - 4.542 ** 2, eAu800 = e151 + (1239.84 / 800 - 1.51) / 0.13 * (e164 - e151), epsWater = 1.33 ** 2, L800 = epsWater / (epsWater - eAu800);
// --- Arrhenius/depletion growth model at its reference point (90 °C, default kinetics): 40 mM, 3 h
const gRef = 3 * 3600 * (1 - Math.exp(-1)), dcRef = (40 - 0.5) / 25, Lref = 10 / 60 * dcRef * gRef, Dref = 20 + 0.4 / 60 * dcRef * gRef;
// --- smooth test function for the Gaussian-process surrogate
const yT = x => 10 + 5 * Math.sin(x / 2), TABLE = 'x,y\n' + [0, 1, 2, 3, 4, 5, 6, 7].map(x => `${x},${yT(x).toFixed(6)}`).join('\n');

const MIM = { mode: 'resonator', resType: 'mim', fLo: 100, fHi: 150, m1: 0.1, k1: 1e6, muMax: 3, k2Min: 10, k2Max: 1e7 };
const AUX = { mode: 'auxetic', nuT: -1, EratT: 1e-3, rhoMax: 0.3, cellL: 5, tMin: 0.05, seed: 1 };
const MBB = { mode: 'topo', tProb: 'mbb', nelx: 60, nely: 20, volfrac: 0.5, penal: 3, rminF: 1.5, tIter: 80 };

export default [
  {
    id: 's25-auxetic-target-poisson', suite: 25, title: 'Inverse honeycomb design reaches ν₁₂ = −1', mode: 'auxetic',
    inputs: AUX, kpi: 'nu12', expected: -1, tol: 0.002, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: ν₁₂ = cos²θ/((h/l + sin θ) sin θ) attains −1 inside the admissible re-entrant range (for example θ = −30°, h/l = 2), so an exact optimiser must return the target',
  },
  {
    id: 's25-auxetic-target-stiffness', suite: 25, title: 'Inverse honeycomb design reaches E₁/E_s = 10⁻³ simultaneously', mode: 'auxetic',
    inputs: AUX, kpi: 'E1', expected: 1e-3, tol: 0.002, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: E₁/E_s = (t/l)³ cos θ/((h/l + sin θ) sin²θ) is free to match any stiffness through t/l once θ and h/l fix ν₁₂',
  },
  {
    id: 's25-auxetic-reciprocal-poisson', suite: 25, title: 'Flexure-model identity ν₂₁ = 1/ν₁₂ at the optimum', mode: 'auxetic',
    inputs: AUX, kpi: 'nu21', expected: -1, tol: 0.002, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: in the cell-wall bending model ν₁₂ν₂₁ = 1 (reciprocity E₁ν₂₁ = E₂ν₁₂)',
  },
  {
    id: 's25-ar-single-layer-silicon', suite: 25, title: 'Single Al₂O₃ layer on silicon: minimum reflectance', mode: 'optical',
    inputs: { mode: 'optical', optGoal: 'ar', nL: 'Al2O3', nSub: 'silicon', nLay: 1, lam1: 1000, lam2: 1001, dMin: 10, dMax: 300 }, kpi: 'Rmean', expected: rQW(1.76, 3.48), tol: 0.005, kind: 'analytical',
    reference: 'Quarter-wave antireflection layer: R = ((n₀n_s − n₁²)/(n₀n_s + n₁²))², zero when n₁ = √(n₀n_s) (Macleod, Thin-Film Optical Filters; Born & Wolf)',
    note: 'n₁ = 1.76 is close to √3.48 = 1.87, leaving R = 0.34 %; the optimiser must find the quarter-wave thickness.',
  },
  {
    id: 's25-ar-mgf2-on-glass', suite: 25, title: 'Single MgF₂ layer on crown glass: residual reflectance 1.26 %', mode: 'optical',
    inputs: { mode: 'optical', optGoal: 'ar', nL: 'MgF2', nSub: 'glass', nLay: 1, lam1: 550, lam2: 551, dMin: 10, dMax: 300 }, kpi: 'Rmean', expected: rQW(1.38, 1.52), tol: 0.005, kind: 'analytical',
    reference: 'Quarter-wave MgF₂ (n = 1.38) on glass (n = 1.52): R = 1.26 % against 4.26 % uncoated (Hecht, Optics; Macleod, Thin-Film Optical Filters)',
  },
  {
    id: 's25-hr-quarter-wave-stack', suite: 25, title: 'Twelve-layer Ta₂O₅/SiO₂ mirror: peak reflectance of the quarter-wave stack', mode: 'optical',
    inputs: { mode: 'optical', optGoal: 'hr', nH: 'Ta2O5', nL: 'SiO2', nSub: 'silica', nLay: 12, lam1: 1064, lam2: 1065, dMin: 10, dMax: 400 }, kpi: 'Rmean', expected: rStack(2.1, 1.46, 1.46, 6), tol: 0.002, kind: 'analytical',
    reference: 'Quarter-wave stack air|(HL)ᴺ|substrate: Y = n_s(n_H/n_L)^(2N), R = ((1 − Y)/(1 + Y))², the maximum reflectance attainable at one wavelength (Macleod, Thin-Film Optical Filters)',
  },
  {
    id: 's25-optical-energy-conservation', suite: 25, title: 'Lossless multilayer: R + T = 1', mode: 'optical',
    inputs: { mode: 'optical', optGoal: 'ar', nH: 'TiO2', nL: 'SiO2', nSub: 'glass', nLay: 6, lam1: 450, lam2: 650 }, kpi: 'cons', expected: 0, abs: 1e-9, kind: 'conservation',
    reference: 'Energy conservation for non-absorbing films (Born & Wolf, Principles of Optics)',
  },
  {
    id: 's25-mass-in-mass-mass-ratio', suite: 25, title: 'Negative-mass band 100–150 Hz: required mass ratio', mode: 'resonator',
    inputs: MIM, kpi: 'mu', expected: (150 / 100) ** 2 - 1, tol: 1e-4, kind: 'analytical',
    reference: 'Huang, Sun & Huang, negative effective mass density in acoustic metamaterials: band ω₂ < ω < ω₂√(1 + m₂/m₁), hence m₂/m₁ = (f_hi/f_lo)² − 1',
  },
  {
    id: 's25-mass-in-mass-inner-stiffness', suite: 25, title: 'Negative-mass band 100–150 Hz: inner spring stiffness', mode: 'resonator',
    inputs: MIM, kpi: 'k2', expected: 0.1 * 1.25 * (2 * PI * 100) ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Local resonance ω₂ = √(k₂/m₂) at the lower band edge: k₂ = m₂(2πf_lo)²',
  },
  {
    id: 's25-split-ring-target', suite: 25, title: 'Split-ring inverse design meets the 10 GHz target', mode: 'resonator',
    inputs: { mode: 'resonator', resType: 'srr', fSrr: 1e10, epsSub: 4.4, tMet: 35, minFeat: 0.15 }, kpi: 'f0', expected: 1e10, tol: 0.01, kind: 'analytical',
    reference: 'L–C resonance f₀ = 1/(2π√(LC)) is monotonic in the ring radius, so a root exists inside the admissible radius range and must be returned',
  },
  {
    id: 's25-helmholtz-target', suite: 25, title: 'Helmholtz inverse design meets the 250 Hz target', mode: 'resonator',
    inputs: { mode: 'resonator', resType: 'helm', fHelm: 250, cAir: 343, rMinH: 2, VmaxH: 500 }, kpi: 'f0', expected: 250, tol: 0.01, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics: f₀ = (c/2π)√(S/(V·L′)) can be met exactly by the cavity volume V = S·c²/((2πf)²L′)',
  },
  {
    id: 's25-brus-cdse-diameter', suite: 25, title: 'CdSe quantum dot emitting at 580 nm: Brus diameter', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'qd', qdMat: 'CdSe', lamE: 580 }, kpi: 'D', expected: 2 * Rbrus * 1e9, tol: 0.001, kind: 'analytical',
    reference: 'Brus equation E = E_g + ħ²π²/(2μR²) − 1.786e²/(4πε₀ε_rR) with E_g = 1.74 eV, m_e = 0.13, m_h = 0.45, ε_r = 10.6, solved as a quadratic in 1/R with CODATA constants',
  },
  {
    id: 's25-cdse-bohr-radius', suite: 25, title: 'Exciton Bohr radius of CdSe', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'qd', qdMat: 'CdSe', lamE: 580 }, kpi: 'aB', expected: 5.6, tol: 0.02, kind: 'reference-data',
    reference: 'Exciton Bohr radius of bulk CdSe, about 5.6 nm (quantum-dot literature; a_B = ε_r·a₀·m₀/μ)',
  },
  {
    id: 's25-gold-sphere-plasmon', suite: 25, title: 'Small gold sphere in water: dipole plasmon wavelength', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 400, nMed: 1.33, dEq: 20 }, kpi: 'lam', expected: 520, tol: 0.02, kind: 'reference-data',
    reference: 'Localised surface-plasmon band of 10–20 nm gold colloids in water at about 520 nm (Fröhlich condition Re ε = −2ε_m; Bohren & Huffman; Kreibig & Vollmer)',
    note: 'A target below the sphere resonance returns the sphere itself; the Drude–Lorentz fit and the quasi-static limit place the peak within a few nanometres of the measured band.',
  },
  {
    id: 's25-gans-sphere-limit', suite: 25, title: 'Depolarisation factor of a sphere', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 400, nMed: 1.33 }, kpi: 'L', expected: 1 / 3, tol: 1e-9, kind: 'analytical',
    reference: 'Depolarisation factors of an ellipsoid sum to one; by symmetry L = 1/3 for a sphere (Bohren & Huffman, Absorption and Scattering of Light by Small Particles)',
  },
  {
    id: 's25-gans-needle-factor', suite: 25, title: 'Depolarisation factor of a prolate spheroid of aspect ratio 12', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 1800, nMed: 1.33 }, kpi: 'L', expected: gansL(12), tol: 1e-6, kind: 'analytical',
    reference: 'Gans theory: L = ((1 − e²)/e²)(ln((1 + e)/(1 − e))/(2e) − 1), e² = 1 − 1/AR² (Bohren & Huffman)',
    note: 'A target beyond the largest admissible aspect ratio returns the bound AR = 12.',
  },
  {
    id: 's25-gold-nanorod-800nm-factor', suite: 25, title: 'Gold nanorod resonant at 800 nm in water: depolarisation factor', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 800, nMed: 1.33, dEq: 20 }, kpi: 'L', expected: L800, tol: 0.02, kind: 'reference-data',
    reference: 'Gans resonance condition Re ε = −ε_m(1 − L)/L with the measured permittivity of gold at 1.55 eV, Re ε ≈ −24.2 (Johnson & Christy, optical constants of the noble metals)',
    note: 'The suite uses an analytical Drude–Lorentz fit and locates the absorption maximum rather than the Fröhlich point; both differences are below 1 % here.',
  },
  {
    id: 's25-gold-nanorod-800nm-aspect', suite: 25, title: 'Gold nanorod resonant at 800 nm in water: aspect ratio', mode: 'nano',
    inputs: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 800, nMed: 1.33, dEq: 20 }, kpi: 'AR', expected: arOfL(L800), tol: 0.03, kind: 'reference-data',
    reference: 'Gans theory for a prolate gold spheroid with the Johnson & Christy permittivity at 800 nm: aspect ratio about 4.3 (measured seed-grown rods need about 4 because of their spherocylinder shape)',
  },
  {
    id: 's25-recipe-reference-temperature', suite: 25, title: 'Recipe optimisation recovers the temperature that sets the aspect ratio', mode: 'recipe',
    inputs: { mode: 'recipe', Ltar: Lref, Dtar: Dref, errTol: 0.5, TmaxR: 120, cMaxR: 100, seed: 1 }, kpi: 'T', expected: 90, abs: 0.3, kind: 'analytical',
    reference: 'Arrhenius kinetics: L/(D − D₀) = (k_ax/k_rad)(T_ref)·exp(−(E_ax − E_rad)/R·(1/T − 1/T_ref)); targets generated at the 90 °C reference point can only be met at 90 °C',
    note: 'Targets are the closed-form model output at 90 °C, 40 mM, 3 h; 0.5 % quality tolerance corresponds to ±0.3 K through d ln(AR)/dT = (E_ax − E_rad)/(RT²).',
  },
  {
    id: 's25-recipe-aspect-ratio', suite: 25, title: 'Recipe optimisation reproduces the target length-to-diameter ratio', mode: 'recipe',
    inputs: { mode: 'recipe', Ltar: Lref, Dtar: Dref, errTol: 0.5, TmaxR: 120, cMaxR: 100, seed: 1 }, kpi: 'AR', expected: Lref / Dref, tol: 0.007, kind: 'analytical',
    reference: 'First-order growth with exponential precursor depletion: L = k_ax(c₀ − c_eq)τ(1 − e^(−t/τ)), D = D₀ + k_rad(c₀ − c_eq)τ(1 − e^(−t/τ))',
  },
  {
    id: 's25-mbb-initial-compliance', suite: 25, title: 'MBB half-beam 60 × 20: compliance of the uniform starting design', mode: 'topo',
    inputs: MBB, kpi: 'c0', expected: 1007.02, tol: 0.001, kind: 'numerical-reference',
    reference: 'Sigmund, a 99 line topology optimization code written in Matlab: top(60,20,0.5,3.0,1.5) reports an objective of 1007.02 in the first iteration',
  },
  {
    id: 's25-mbb-optimised-compliance', suite: 25, title: 'MBB half-beam 60 × 20: optimised compliance', mode: 'topo',
    inputs: MBB, kpi: 'c', expected: 203.3, tol: 0.005, kind: 'numerical-reference',
    reference: 'Sigmund, a 99 line topology optimization code written in Matlab: top(60,20,0.5,3.0,1.5) converges to a compliance of about 203.3',
    note: 'The suite limits the run to 77 iterations for this grid; the published run stops later, and the objective changes by less than 0.1 % in between.',
  },
  {
    id: 's25-mbb-volume-constraint', suite: 25, title: 'Topology optimisation holds the volume fraction', mode: 'topo',
    inputs: MBB, kpi: 'vol', expected: 0.5, abs: 1e-3, kind: 'conservation',
    reference: 'Optimality-criteria update with bisection on the Lagrange multiplier enforces Σρₑ = V_f·N (Bendsøe & Sigmund, Topology Optimization)',
  },
  {
    id: 's25-bayes-growth-temperature', suite: 25, title: 'Bayesian optimisation locates the growth temperature of a known optimum', mode: 'bayes',
    inputs: { mode: 'bayes', bModel: 'growth1d', Ltar: Lref, n0: 6, nBo: 15, seed: 1 }, kpi: 'xbest', expected: 90, abs: 0.2, kind: 'analytical',
    reference: 'Objective (L(T)/L* − 1)² with L* generated at 90 °C by the Arrhenius/depletion model has its unique zero at 90 °C',
    note: 'The acquisition grid has a spacing of 0.175 K, which bounds the attainable accuracy.',
  },
  {
    id: 's25-bayes-auxetic-angle', suite: 25, title: 'Bayesian optimisation finds the re-entrant angle giving ν₁₂ = −1 at h/l = 2', mode: 'bayes',
    inputs: { mode: 'bayes', bModel: 'lattice', nuT: -1, n0: 6, nBo: 20, seed: 1 }, kpi: 'xbest', expected: -30, abs: 0.35, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: cos²θ/((2 + sin θ) sin θ) = −1 has the single root θ = −30° in the re-entrant range',
    note: 'The acquisition grid has a spacing of 0.31°.',
  },
  {
    id: 's25-gp-interpolation-end-point', suite: 25, title: 'Gaussian process reproduces a training value (first point)', mode: 'bayes',
    inputs: { mode: 'bayes', bModel: 'table', bTable: TABLE, bGoal: 'max', xq: 0 }, kpi: 'q_mu', expected: yT(0), tol: 1e-3, kind: 'analytical',
    reference: 'Rasmussen & Williams, Gaussian Processes for Machine Learning: the posterior mean of a noise-free GP interpolates the training data',
    note: 'The marginal likelihood selects the smallest jitter (10⁻⁵ of the signal variance) for smooth data, so the interpolation is exact to about 10⁻⁴.',
  },
  {
    id: 's25-gp-interpolation-interior-point', suite: 25, title: 'Gaussian process reproduces a training value (interior point)', mode: 'bayes',
    inputs: { mode: 'bayes', bModel: 'table', bTable: TABLE, bGoal: 'max', xq: 3 / 7 }, kpi: 'q_mu', expected: yT(3), tol: 1e-3, kind: 'analytical',
    reference: 'Rasmussen & Williams, Gaussian Processes for Machine Learning: the posterior mean of a noise-free GP interpolates the training data',
  },
  {
    id: 's25-gp-table-incumbent', suite: 25, title: 'Table data: the incumbent is the largest tabulated response', mode: 'bayes',
    inputs: { mode: 'bayes', bModel: 'table', bTable: TABLE, bGoal: 'max', xq: 0.5 }, kpi: 'best', expected: Math.max(...[0, 1, 2, 3, 4, 5, 6, 7].map(x => +yT(x).toFixed(6))), tol: 1e-9, kind: 'conservation',
    reference: 'Definition of the incumbent in expected-improvement Bayesian optimisation (Jones, Schonlau & Welch, efficient global optimization)',
  },
  {
    id: 's25-knn-identical-composition-density', suite: 25, title: 'Nearest-neighbour screening returns the density of an identical composition (ZnO)', mode: 'knn',
    inputs: { mode: 'knn', kProp: 'density', kMat: 'zno-nanoparticles', kNN: 3 }, kpi: 'pred', expected: 5606, tol: 0.005, kind: 'reference-data',
    reference: 'Density of zinc oxide, 5606 kg/m³ (CRC Handbook of Chemistry and Physics); an inverse-distance-weighted neighbour model must reproduce a training entry of the same composition',
    note: 'The library holds ZnO in two morphologies; the selected entry is withheld and its twin lies at zero descriptor distance.',
  },
  {
    id: 's25-knn-identical-composition-bandgap', suite: 25, title: 'Nearest-neighbour screening returns the band gap of an identical composition (ZnO)', mode: 'knn',
    inputs: { mode: 'knn', kProp: 'bandgap_eV', kMat: 'zno-nanoparticles', kNN: 3 }, kpi: 'pred', expected: 3.37, tol: 0.005, kind: 'reference-data',
    reference: 'Room-temperature band gap of wurtzite ZnO, 3.37 eV (semiconductor data handbooks)',
  },
];
