// Reference benchmarks for suite 4 — Nanomaterial Geometry and Morphology.
// Expected values are closed-form geometry, distribution moments and published packing results evaluated here.

const PI = Math.PI, N_A = 6.02214076e23;
const RHO = { ZnO: 5606, Au: 19300, C: 2260 };   // handbook bulk densities, kg/m³
// log-normal moment ratios (Hatch–Choate): D[p,q] = d_med·exp((p + q)/2 · ln²σ_g)
const ln2 = g => Math.log(g) ** 2;
// gamma function (Lanczos, g = 7)
const gamma = z => { const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; z -= 1; let x = c[0]; for (let i = 1; i < 9; i++) x += c[i] / (z + i); const t = z + 7.5; return Math.sqrt(2 * PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x; };
// prolate spheroid, semi-axes a = 10 nm, c = 20 nm
const eP = Math.sqrt(1 - 0.25), areaProlate = 2 * PI * 1e-16 * (1 + 2 / eP * Math.asin(eP));
// thin circular rod and thin sheet
const rodSSA = (rho, d, L) => 4 / (rho * d) * (1 + d / (2 * L));
// Boolean model of overlapping spheres (diameter D, porosity φ): S_v = −(6/D)·φ·ln φ
const phiB = 0.25, DB = 2e-8, SvB = -6 / DB * phiB * Math.log(phiB);
// island growth: hemispherical caps, Avrami overlap
const NI = 1e15, r99 = Math.sqrt(Math.log(100) / (PI * NI)), covAt = h => 1 - Math.exp(-NI * PI * Math.cbrt(3 * h / (2 * PI * NI)) ** 2);
const ROD = { mode: 'rods', material: 'ZnO', arrangement: 'square', tilt: 0, cvL: 0, cvD: 0, d: 1e-7, L: 1e-6, nDens: 1e13 };
const POP = { mode: 'population', popSource: 'generate', material: 'Au', d: 1.5e-8, nSample: 20000, seed: 1 };

export default [
  {
    id: 's04-zno-sphere-ssa', suite: 4, title: 'Specific surface area of a 10 nm ZnO sphere', mode: 'single',
    inputs: { mode: 'single', material: 'ZnO', shape: 'sphere', d: 1e-8 }, kpi: 'ssa', expected: 6 / (RHO.ZnO * 1e-8) / 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'SSA = 6/(ρd) for dense spheres with ρ(ZnO) = 5606 kg/m³: 107 m²/g (BET-equivalent diameter relation)',
  },
  {
    id: 's04-au-cube-ssa', suite: 4, title: 'Specific surface area of a 20 nm gold cube', mode: 'single',
    inputs: { mode: 'single', material: 'Au', shape: 'cube', d: 2e-8 }, kpi: 'ssa', expected: 6 / (RHO.Au * 2e-8) / 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'Cube geometry: A/V = 6/d, the same as a sphere of diameter d; ρ(Au) = 19 300 kg/m³',
  },
  {
    id: 's04-cube-sphericity', suite: 4, title: 'Wadell sphericity of a cube', mode: 'single',
    inputs: { mode: 'single', material: 'Au', shape: 'cube', d: 2e-8 }, kpi: 'psi', expected: Math.cbrt(PI / 6), tol: 1e-9, kind: 'analytical',
    reference: 'Wadell sphericity ψ = π^(1/3)(6V)^(2/3)/A = (π/6)^(1/3) = 0.806 for a cube (particle-technology textbooks)',
  },
  {
    id: 's04-octahedron-sphericity', suite: 4, title: 'Wadell sphericity of a regular octahedron', mode: 'single',
    inputs: { mode: 'single', material: 'Au', shape: 'octahedron', d: 2e-8 }, kpi: 'psi', expected: Math.cbrt(PI) * Math.pow(2 * Math.SQRT2, 2 / 3) / (2 * Math.sqrt(3)), tol: 1e-9, kind: 'analytical',
    reference: 'Regular octahedron of edge a: V = (√2/3)a³, A = 2√3·a²; ψ = 0.846 (particle-technology textbooks)',
  },
  {
    id: 's04-prolate-spheroid-area', suite: 4, title: 'Surface area of a 20 nm × 40 nm prolate spheroid', mode: 'single',
    inputs: { mode: 'single', material: 'ZnO', shape: 'spheroid', d: 2e-8, L: 4e-8 }, kpi: 'A', expected: areaProlate, tol: 1e-9, kind: 'analytical',
    reference: 'Prolate spheroid: A = 2πa²[1 + (c/(a·e))·asin e], e² = 1 − a²/c² (21.48 a² for c = 2a; mathematical handbooks)',
  },
  {
    id: 's04-rod-ssa', suite: 4, title: 'Specific surface area of a 50 nm × 10 µm ZnO wire', mode: 'single',
    inputs: { mode: 'single', material: 'ZnO', shape: 'rod_circ', d: 5e-8, L: 1e-5 }, kpi: 'ssa', expected: rodSSA(RHO.ZnO, 5e-8, 1e-5) / 1e3, tol: 1e-9, kind: 'analytical',
    reference: 'Cylinder with end faces: A/V = (4/d)(1 + d/(2L)), tending to 4/(ρd) for long wires',
  },
  {
    id: 's04-rod-sauter-diameter', suite: 4, title: 'Sauter diameter of a long cylinder tends to 1.5 d', mode: 'single',
    inputs: { mode: 'single', material: 'ZnO', shape: 'rod_circ', d: 5e-8, L: 1e-5 }, kpi: 'd32', expected: 1.5 * 5e-8 / (1 + 5e-8 / 2e-5), tol: 1e-9, kind: 'analytical',
    reference: 'd₃₂ = 6V/A = 1.5d/(1 + d/(2L)) for a cylinder (particle-technology textbooks)',
  },
  {
    id: 's04-graphitic-sheet-ssa', suite: 4, title: 'Sheet of graphitic carbon one layer spacing thick: SSA of graphene', mode: 'single',
    inputs: { mode: 'single', material: 'C', shape: 'sheet', d: 1e-6, tw: 3.35e-10 }, kpi: 'ssa', expected: 2630, tol: 0.01, kind: 'reference-data',
    reference: 'Theoretical specific surface area of graphene, 2630 m²/g; equivalently 2/(ρt) with graphite density 2260 kg/m³ and layer spacing 0.335 nm',
  },
  {
    id: 's04-gold-atoms-5nm', suite: 4, title: 'Number of atoms in a 5 nm gold sphere', mode: 'single',
    inputs: { mode: 'single', material: 'Au', shape: 'sphere', d: 5e-9 }, kpi: 'nAtoms', expected: PI / 6 * 5e-9 ** 3 * RHO.Au * N_A / 0.196967, tol: 0.001, kind: 'analytical',
    reference: 'N = (π/6)d³·ρ·N_A/M with ρ = 19 300 kg/m³ and M = 196.967 g/mol: about 3860 atoms',
  },
  {
    id: 's04-ssa-tolerance-sensitivity', suite: 4, title: 'Worst-case SSA change of a sphere for ±5 % diameter tolerance', mode: 'single',
    inputs: { mode: 'single', material: 'ZnO', shape: 'sphere', d: 1e-8, tolPct: 5 }, kpi: 'ssaTol', expected: 5, tol: 1e-3, kind: 'analytical',
    reference: 'First-order error propagation: SSA ∝ 1/d, so |ΔSSA/SSA| = |Δd/d|',
  },
  {
    id: 's04-lognormal-number-mean', suite: 4, title: 'Log-normal population: number mean (Hatch–Choate)', mode: 'population',
    inputs: { ...POP, dist: 'lognormal', gsd: 1.3 }, kpi: 'd10m', expected: 1.5e-8 * Math.exp(0.5 * ln2(1.3)), tol: 0.01, kind: 'analytical',
    reference: 'Hatch–Choate relation D[1,0] = d_med·exp(0.5 ln²σ_g) (aerosol and particle-size textbooks)',
    note: 'Seeded sample of 20 000 particles; the tolerance covers the sampling error.',
  },
  {
    id: 's04-lognormal-sauter-mean', suite: 4, title: 'Log-normal population: Sauter mean (Hatch–Choate)', mode: 'population',
    inputs: { ...POP, dist: 'lognormal', gsd: 1.3 }, kpi: 'd32', expected: 1.5e-8 * Math.exp(2.5 * ln2(1.3)), tol: 0.01, kind: 'analytical',
    reference: 'Hatch–Choate relation D[3,2] = d_med·exp(2.5 ln²σ_g)',
    note: 'Seeded sample of 20 000 particles; the tolerance covers the sampling error.',
  },
  {
    id: 's04-lognormal-volume-mean', suite: 4, title: 'Log-normal population: volume-weighted mean (Hatch–Choate)', mode: 'population',
    inputs: { ...POP, dist: 'lognormal', gsd: 1.3 }, kpi: 'd43', expected: 1.5e-8 * Math.exp(3.5 * ln2(1.3)), tol: 0.015, kind: 'analytical',
    reference: 'Hatch–Choate relation D[4,3] = d_med·exp(3.5 ln²σ_g)',
    note: 'Fourth-moment statistics of 20 000 sampled particles carry about 0.5 % sampling error.',
  },
  {
    id: 's04-lognormal-population-ssa', suite: 4, title: 'Specific surface area of a log-normal gold sol', mode: 'population',
    inputs: { ...POP, dist: 'lognormal', gsd: 1.3 }, kpi: 'ssa', expected: 6 / (RHO.Au * 1.5e-8 * Math.exp(2.5 * ln2(1.3))) / 1e3, tol: 0.01, kind: 'analytical',
    reference: 'SSA = 6/(ρ·D[3,2]) with the Hatch–Choate Sauter mean',
  },
  {
    id: 's04-normal-sauter-mean', suite: 4, title: 'Normal population (CV 15 %): Sauter mean', mode: 'population',
    inputs: { ...POP, dist: 'normal', cvPct: 15 }, kpi: 'd32', expected: 1.5e-8 * (1 + 3 * 0.15 ** 2) / (1 + 0.15 ** 2), tol: 0.005, kind: 'analytical',
    reference: 'Moments of the normal distribution: E[d³]/E[d²] = μ(1 + 3c²)/(1 + c²), c = σ/μ',
  },
  {
    id: 's04-weibull-number-mean', suite: 4, title: 'Weibull population (shape 4): number mean from the median', mode: 'population',
    inputs: { ...POP, dist: 'weibull', wk: 4 }, kpi: 'd10m', expected: 1.5e-8 / Math.pow(Math.log(2), 0.25) * gamma(1.25), tol: 0.005, kind: 'analytical',
    reference: 'Weibull distribution: median = λ(ln 2)^(1/k), mean = λΓ(1 + 1/k)',
  },
  {
    id: 's04-rod-array-coverage', suite: 4, title: 'Areal coverage of vertical circular rods on a square pattern', mode: 'rods',
    inputs: { ...ROD, section: 'round' }, kpi: 'coverage', expected: 1e13 * PI * 1e-14 / 4, tol: 1e-6, kind: 'analytical',
    reference: 'Coverage = n·πd²/4 for 10 rods/µm² of 100 nm diameter',
  },
  {
    id: 's04-rod-array-area-enhancement', suite: 4, title: 'Surface-area enhancement of a vertical nanorod array', mode: 'rods',
    inputs: { ...ROD, section: 'round' }, kpi: 'sef', expected: 1 + 1e13 * PI * 1e-7 * 1e-6, tol: 1e-6, kind: 'analytical',
    reference: 'Roughness factor of a rod forest: 1 + n·πdL (side walls added to the projected area)',
  },
  {
    id: 's04-rod-array-hexagonal-section', suite: 4, title: 'Areal coverage of hexagonal-prism rods (d across flats)', mode: 'rods',
    inputs: { ...ROD, section: 'hex' }, kpi: 'coverage', expected: 1e13 * Math.sqrt(3) / 2 * 1e-14, tol: 1e-6, kind: 'analytical',
    reference: 'Regular hexagon of across-flats width d: area (√3/2)d²',
  },
  {
    id: 's04-rod-array-mass-loading', suite: 4, title: 'Mass loading of a ZnO nanorod array', mode: 'rods',
    inputs: { ...ROD, section: 'round' }, kpi: 'mass', expected: RHO.ZnO * 1e13 * PI * 1e-14 / 4 * 1e-6 * 100, tol: 1e-6, kind: 'analytical',
    reference: 'm″ = ρ·n·(πd²/4)·L, converted to mg/cm²',
  },
  {
    id: 's04-rod-array-ssa', suite: 4, title: 'Specific surface area of standing rods (side wall and top)', mode: 'rods',
    inputs: { ...ROD, section: 'round' }, kpi: 'ssa', expected: 4 / (RHO.ZnO * 1e-7) * (1 + 1e-7 / 4e-6) / 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'Cylinder attached at its base: (πdL + πd²/4)/(ρ·πd²L/4) = (4/(ρd))(1 + d/(4L))',
  },
  {
    id: 's04-rod-array-order-parameter', suite: 4, title: 'Orientation order parameter of perfectly vertical rods', mode: 'rods',
    inputs: { ...ROD, section: 'round' }, kpi: 'P2', expected: 1, tol: 1e-9, kind: 'analytical',
    reference: 'Hermans orientation factor ⟨(3cos²θ − 1)/2⟩ = 1 for perfect alignment',
  },
  {
    id: 's04-boolean-model-surface', suite: 4, title: 'Overlapping-sphere microstructure: specific surface of the Boolean model', mode: 'porous',
    inputs: { mode: 'porous', material: 'ZnO', pType: 'randomSpheres', phi: phiB, feat: DB, featVox: 12, pN: 48, seed: 1 }, kpi: 'ssa_ens', expected: SvB / (RHO.ZnO * (1 - phiB)) / 1e3, tol: 0.03, kind: 'analytical',
    reference: 'Fully penetrable spheres: S_v = −(6/D)·φ·ln φ (Torquato, Random Heterogeneous Materials)',
    note: 'Mean over twelve seeded 48³ volumes at the finest resolution the suite offers (12 voxels per sphere diameter). The surface is taken from the two-point function along 13 lattice directions extrapolated to zero lag (S_v = −4·dS₂/dr); the ensemble standard error is about 0.8 % (KPI ssa_se) and the truncation error about −1 %. At 9 voxels per diameter the truncation error grows to about −2 %, and a single volume scatters by about 3 %.',
  },
  {
    id: 's04-kozeny-carman-boolean', suite: 4, title: 'Kozeny–Carman permeability of the overlapping-sphere medium', mode: 'porous',
    inputs: { mode: 'porous', material: 'ZnO', pType: 'randomSpheres', phi: phiB, feat: DB, featVox: 12, pN: 48, seed: 1 }, kpi: 'k_ens', expected: phiB ** 3 / (5 * SvB * SvB), tol: 0.04, kind: 'analytical',
    reference: 'Kozeny–Carman equation k = φ³/(5·S_v²) with the Boolean-model specific surface (Carman; Torquato, Random Heterogeneous Materials)',
    note: 'k depends on 1/S_v², so the standard error (about 0.8 %) and truncation error (about −1 %) of the ensemble surface area appear doubled in the permeability.',
  },
  {
    id: 's04-rsa-jamming-coverage', suite: 4, title: 'Random sequential adsorption of discs: jamming coverage', mode: 'aggregate',
    inputs: { mode: 'aggregate', aggType: 'rsa', rsaL: 60, rsaTau: 400, seed: 1 }, kpi: 'thInf', expected: 0.547, tol: 0.02, kind: 'numerical-reference',
    reference: 'Jamming limit of random sequential adsorption of discs on a plane, θ_J = 0.547 (Feder; Hinrichsen, Feder & Jøssang)',
    note: 'Single seeded realisation on a 60 × 60 diameter periodic domain extrapolated with Feder’s law; realisation scatter is about ±0.005.',
  },
  {
    id: 's04-dla-fractal-dimension', suite: 4, title: 'Three-dimensional diffusion-limited aggregation: fractal dimension', mode: 'aggregate',
    inputs: { mode: 'aggregate', aggType: 'dla', nPart: 1500, stick: 1, seed: 1 }, kpi: 'Df', expected: 2.5, tol: 0.05, kind: 'numerical-reference',
    reference: 'Particle–cluster DLA in three dimensions, D_f ≈ 2.5 (Witten & Sander; Meakin)',
    note: 'Ensemble fit over twelve seeded on-lattice clusters of 1500 particles, window R_g ≥ 4 lattice spacings (about one decade in N). Other seeds give 2.40–2.50 with a standard error of about 0.04 (KPI Df_se); the short scaling range biases the value low by 1–2 %, and the asymptotic value needs far larger clusters.',
  },
  {
    id: 's04-film-areal-mass', suite: 4, title: 'Areal mass of a ZnO / Al₂O₃ / TiO₂ stack', mode: 'film',
    inputs: { mode: 'film', stack: 'layer,formula,thickness_nm,density_kg_m3\nSeed,ZnO,30,5606\nBarrier,Al2O3,10,3100\nFilm,TiO2,60,3890' }, kpi: 'massA', expected: (30 * 5606 + 10 * 3100 + 60 * 3890) * 1e-9 * 1e5, tol: 1e-9, kind: 'analytical',
    reference: 'm″ = Σρᵢtᵢ (quartz-crystal microbalance relation), in µg/cm²',
  },
  {
    id: 's04-film-closure-thickness', suite: 4, title: 'Island growth: nominal thickness for 99 % coverage', mode: 'film',
    inputs: { mode: 'film', islandN: NI }, kpi: 'h99', expected: 2 * PI * NI * r99 ** 3 / 3, tol: 1e-6, kind: 'analytical',
    reference: 'Hemispherical islands with random (Avrami–Kolmogorov) overlap: θ = 1 − exp(−Nπr²), h = N(2/3)πr³',
  },
  {
    id: 's04-film-island-coverage', suite: 4, title: 'Island growth: coverage of a 60 nm top layer at 10¹⁵ islands/m²', mode: 'film',
    inputs: { mode: 'film', islandN: NI }, kpi: 'cov', expected: covAt(6e-8), tol: 1e-6, kind: 'analytical',
    reference: 'Avrami–Kolmogorov extended-area relation for randomly nucleated hemispherical islands',
  },
  {
    id: 's04-film-roughness-ratio', suite: 4, title: 'Gaussian rough surface: R_a = R_q·√(2/π)', mode: 'film',
    inputs: { mode: 'film', rmsT: 2e-9, hurst: 0.7, xi: 5e-8, scan: 1e-6, seed: 1 }, kpi: 'ra', expected: 2e-9 * Math.sqrt(2 / PI), tol: 0.03, kind: 'analytical',
    reference: 'Normal height distribution: mean absolute deviation = σ√(2/π) (surface-metrology textbooks: R_a ≈ 0.80 R_q)',
    note: 'One seeded 64 × 64 realisation of a correlated surface; the sampled height distribution is only approximately normal.',
  },
  {
    id: 's04-step-coverage-line-of-sight', suite: 4, title: 'Line-of-sight flux at the bottom of an aspect-ratio-20 trench', mode: 'film',
    inputs: { mode: 'film', featAR: 20 }, kpi: 'scPVD', expected: 1 / Math.sqrt(1 + 4 * 400), tol: 1e-9, kind: 'analytical',
    reference: 'Cosine-law flux from a diffuse source through a slot: sin(arctan(w/2h)) = 1/√(1 + 4·AR²) (vacuum-deposition view-factor geometry)',
  },
  {
    id: 's04-step-coverage-knudsen', suite: 4, title: 'Diffusion–reaction step coverage in a via (Thiele modulus)', mode: 'film',
    inputs: { mode: 'film', featAR: 20, sticking: 0.01 }, kpi: 'scCVD', expected: 1 / Math.cosh(20 * Math.sqrt(0.03)), tol: 1e-9, kind: 'analytical',
    reference: 'First-order wall reaction with Knudsen diffusion in a cylindrical pore: n(H)/n(0) = 1/cosh Φ, Φ = AR·√(3s) (Thiele; chemical reaction engineering textbooks)',
  },
  {
    id: 's04-ald-exposure-multiplier', suite: 4, title: 'ALD saturation exposure of an aspect-ratio-20 hole', mode: 'film',
    inputs: { mode: 'film', featAR: 20 }, kpi: 'ald', expected: 1 + 19 / 4 * 20 + 1.5 * 400, tol: 1e-9, kind: 'analytical',
    reference: 'Gordon, Hausmann, Kim & Shepard kinetic model of conformal ALD in narrow holes: exposure ∝ 1 + (19/4)a + (3/2)a²',
  },
  // ------------------------------------------------------------------ published experiments
  {
    id: 's04-exp-rsa-random-parking-limit', suite: 4, title: 'Measured random-parking (jamming) limit of spheres adsorbed on a plane', mode: 'aggregate',
    inputs: { mode: 'aggregate', aggType: 'rsa', rsaL: 60, rsaTau: 400, seed: 1 }, kpi: 'thInf',
    expected: 0.55, tol: 0.03, kind: 'experiment',
    measured: { value: 0.55, uncertainty: 0.01, unit: '– (area fraction)' },
    source: { authors: 'Onoda, G. Y.; Liniger, E. G.', title: 'Experimental determination of the random-parking limit in two dimensions', journal: 'Phys. Rev. A', year: 1986, doi: '10.1103/PhysRevA.33.715', url: 'https://journals.aps.org/pra/abstract/10.1103/PhysRevA.33.715' },
    conditions: 'Spherical colloidal particles adsorbing irreversibly, without overlap, on a flat substrate until no further particle fits; covered area fraction measured',
    match: ['aggType'],
    reference: 'Onoda & Liniger, Phys. Rev. A 33, 715 (1986)',
    note: 'Value 0.55 ± 0.01 from the abstract. The suite runs a seeded random sequential adsorption of discs on the largest domain (60 × 60 diameters, 400 attempts per disc area) and extrapolates to the jamming limit with Feder\'s law; realisation scatter is about ±0.005, so the tolerance is the experimental uncertainty plus that scatter.',
  },
];
