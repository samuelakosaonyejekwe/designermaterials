// Reference benchmarks for Suite 17 — Quantum, atomistic and molecular simulation.
// Every expected value is a literal reference value or a closed form evaluated here from CODATA constants;
// nothing in this file calls the suite.
const PI = Math.PI;
const H = 6.62607015e-34, HBAR = H / (2 * PI), KB = 1.380649e-23, QE = 1.602176634e-19, U = 1.66053906660e-27, M0 = 9.1093837015e-31;
const KB_EV = KB / QE;
const HB2M = HBAR * HBAR / (2 * M0) / QE * 1e18;          // ħ²/2m₀ in eV nm²
const COUL = QE / (4 * PI * 8.8541878128e-12) * 1e9;      // e²/4πε₀ in eV nm

// Lennard-Jones argon (ε/k_B = 119.8 K, σ = 3.405 Å) and the fcc lattice sums A₁₂, A₆
const LJ_E = 119.8 * KB, LJ_S = 3.405e-10, A12 = 12.13188, A6 = 14.45392;
const ljNN = Math.pow(2 * A12 / A6, 1 / 6) * LJ_S;       // static nearest-neighbour distance, 1.0902 σ

// graphene / nanotubes
const ACC = 1.42e-10, G0 = 2.7, AG = Math.sqrt(3) * ACC;

// finite square well, equal masses: even-parity ground state from k tan(kW/2) = κ (bisection)
function wellE1(W, V0, m) {
  const f = E => { const k = Math.sqrt(E * m / HB2M), kap = Math.sqrt((V0 - E) * m / HB2M); return k * Math.tan(k * W / 2) - kap; };
  let lo = 1e-9, hi = Math.min(V0, HB2M * PI * PI / (m * W * W)) * (1 - 1e-9);
  for (let i = 0; i < 200; i++) { const mid = 0.5 * (lo + hi); if (f(mid) > 0) hi = mid; else lo = mid; }
  return 0.5 * (lo + hi);
}

// ideal monatomic gas reference pressure k_BT/Λ³ and Langmuir constant K = exp(E/k_BT) Λ³/k_BT
const langK = (E_eV, T, mass_u) => { const lam = H / Math.sqrt(2 * PI * mass_u * U * KB * T); return Math.exp(E_eV / (KB_EV * T)) * lam ** 3 / (KB * T); };
const K300 = langK(0.5, 300, 28);

// Müller–Brown stationary-point energies (surface units)
const MB_A = -146.6995, MB_B = -80.7678, MB_C = -108.1667, MB_S1 = -40.6648, MB_S2 = -72.2489;

const dia = { m1: 65.38 * U, m2: 15.999 * U, K: 50, a: 3.0e-10 };

export default [
  // ── molecular dynamics ──
  {
    id: 's17-md-nve-energy-conservation', suite: 17, title: 'Liquid argon, NVE: total energy is conserved', mode: 'md',
    inputs: { mode: 'md', mdSys: 'bulk', mdPot: 'Ar', mdEns: 'nve', ncell: 3, dens: 0.75, T: 94.4, dt: 10, nEq: 300, nProd: 1000, seed: 1 },
    kpi: 'drift', expected: 0, abs: 2e-4, kind: 'conservation',
    reference: 'Hamiltonian dynamics conserve total energy; velocity-Verlet integration (Allen & Tildesley, Computer Simulation of Liquids)',
    note: 'Relative drift of the total energy over 10 ps with a 10 fs step; a symplectic integrator with a continuous (shifted-force) potential should hold it far below 10⁻³.'
  },
  {
    id: 's17-md-lj-solid-nn-distance', suite: 17, title: 'Solid Lennard-Jones argon at 5 K: nearest-neighbour distance 1.09 σ', mode: 'md',
    inputs: { mode: 'md', mdSys: 'bulk', mdPot: 'Ar', mdEns: 'nve', ncell: 3, dens: 1, T: 5, dt: 5, nEq: 100, nProd: 400, seed: 1 },
    kpi: 'r_peak', expected: ljNN * 1e10, tol: 0.01, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — fcc Lennard-Jones lattice sums: R₀/σ = (2A₁₂/A₆)^(1/6) = 1.09',
    note: 'First peak of g(r); the histogram bin width (about 0.07 Å) limits the resolution to roughly 1 %.'
  },
  // ── lattice dynamics ──
  {
    id: 's17-phonon-lj-argon-lattice-constant', suite: 17, title: 'fcc Lennard-Jones argon: static lattice constant', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Ar' },
    pub: 'prop.lattice_a', expected: Math.SQRT2 * ljNN, tol: 0.002, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — a = √2 (2A₁₂/A₆)^(1/6) σ with A₁₂ = 12.13188, A₆ = 14.45392 (σ = 3.405 Å gives 5.25 Å)',
    note: 'The suite sums the potential to five nearest-neighbour distances rather than to infinity, which shifts the minimum by a few parts in 10⁴.'
  },
  {
    id: 's17-phonon-lj-argon-bulk-modulus', suite: 17, title: 'fcc Lennard-Jones argon: bulk modulus ≈ 75 ε/σ³', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Ar' },
    kpi: 'B', expected: 4 * LJ_E * A12 * Math.pow(A6 / A12, 2.5) / LJ_S ** 3, tol: 0.02, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — B = 4ε A₁₂ (A₆/A₁₂)^(5/2)/σ³ = 75 ε/σ³ for the fcc Lennard-Jones crystal',
    note: 'The suite obtains B from long-wavelength sound velocities of the dynamical matrix with a finite neighbour sum (five nearest-neighbour distances); the truncation removes about 1 % of the attractive lattice sum.'
  },
  {
    id: 's17-phonon-argon-debye-temperature', suite: 17, title: 'Solid argon: Debye temperature', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Ar' },
    kpi: 'thD', expected: 92, tol: 0.05, kind: 'reference-data',
    reference: 'Kittel, Introduction to Solid State Physics — table of Debye temperatures (Ar: 92 K)',
    note: 'Harmonic lattice dynamics of the static Lennard-Jones crystal; zero-point expansion and anharmonicity, which soften the real crystal by a few per cent, are not modelled.'
  },
  {
    id: 's17-phonon-dulong-petit', suite: 17, title: 'fcc crystal at T ≫ θ_D: Dulong–Petit heat capacity 3 k_B per atom', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Ar', T: 3000 },
    kpi: 'cv', expected: 3, tol: 0.002, kind: 'analytical',
    reference: 'Dulong–Petit law as the high-temperature limit of the harmonic (Einstein/Debye) heat capacity — Ashcroft & Mermin, Solid State Physics',
    note: 'At T ≈ 30 θ_D the leading quantum correction (θ_D/T)²/20 is below 10⁻⁴.'
  },
  {
    id: 's17-phonon-copper-lattice-constant', suite: 17, title: 'fcc copper from the Morse potential: lattice constant', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Cu' },
    pub: 'prop.lattice_a', expected: 3.615e-10, tol: 0.01, kind: 'reference-data',
    reference: 'Measured lattice parameter of copper, 3.615 Å (Kittel; CRC Handbook of Chemistry and Physics)',
    note: 'The Morse parameters were fitted to lattice properties with a longer-ranged neighbour sum than the five nearest-neighbour distances used here.'
  },
  {
    id: 's17-phonon-copper-bulk-modulus', suite: 17, title: 'fcc copper from the Morse potential: bulk modulus', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'fcc', mdPot: 'Cu' },
    kpi: 'B', expected: 137e9, tol: 0.05, kind: 'reference-data',
    reference: 'Bulk modulus of copper, 137 GPa (Kittel, Introduction to Solid State Physics — table of bulk moduli)',
    note: 'Pair-potential model fitted to compressibility; the truncated neighbour sum and the static (0 K) lattice account for a few per cent.'
  },
  {
    id: 's17-phonon-diatomic-optical-frequency', suite: 17, title: 'Diatomic chain: zone-centre optical frequency', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'diatomic', phM1: 65.38, phM2: 15.999, phK: 50, phA: 3.0 },
    kpi: 'nu_max', expected: Math.sqrt(2 * dia.K * (1 / dia.m1 + 1 / dia.m2)) / (2 * PI) * 1e-12, tol: 1e-4, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — two atoms per primitive basis: ω² = 2C(1/M₁ + 1/M₂) at K = 0'
  },
  {
    id: 's17-phonon-diatomic-gap', suite: 17, title: 'Diatomic chain: acoustic–optical gap at the zone boundary', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'diatomic', phM1: 65.38, phM2: 15.999, phK: 50, phA: 3.0 },
    kpi: 'gap', expected: (Math.sqrt(2 * dia.K / dia.m2) - Math.sqrt(2 * dia.K / dia.m1)) / (2 * PI) * 1e-12, tol: 1e-4, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — zone-boundary frequencies √(2C/M₁) and √(2C/M₂) of the diatomic chain'
  },
  {
    id: 's17-phonon-diatomic-sound-velocity', suite: 17, title: 'Diatomic chain: long-wavelength sound velocity', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'diatomic', phM1: 65.38, phM2: 15.999, phK: 50, phA: 3.0 },
    kpi: 'vs', expected: dia.a * Math.sqrt(dia.K / (2 * (dia.m1 + dia.m2))), tol: 1e-4, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — acoustic branch ω² = ½C K²a²/(M₁ + M₂) for small K'
  },
  {
    id: 's17-phonon-monatomic-cutoff', suite: 17, title: 'Monatomic chain: cut-off frequency 2√(K/m)/2π', mode: 'phonon',
    inputs: { mode: 'phonon', phModel: 'mono', phM1: 63.546, phK: 30, phA: 2.5 },
    kpi: 'nu_max', expected: Math.sqrt(30 / (63.546 * U)) / PI * 1e-12, tol: 1e-4, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — monatomic chain ω = √(4C/M) |sin(Ka/2)|'
  },
  // ── electronic structure ──
  {
    id: 's17-graphene-pi-bandwidth', suite: 17, title: 'Graphene nearest-neighbour tight binding: π band width 6γ₀', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'graphene', gamma0: 2.7, tbS: 0 },
    kpi: 'width', expected: 6 * G0, tol: 1e-6, kind: 'analytical',
    reference: 'Wallace band structure of graphite; Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — E(Γ) = ±3γ₀ for s = 0'
  },
  {
    id: 's17-graphene-dirac-point', suite: 17, title: 'Graphene: zero gap at the K point', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'graphene', gamma0: 2.7, tbS: 0 },
    kpi: 'gap', expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — f(K) = 0, the π and π* bands touch at K'
  },
  {
    id: 's17-graphene-fermi-velocity', suite: 17, title: 'Graphene: Fermi velocity 3γ₀a_cc/2ħ', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'graphene', gamma0: 2.7, tbS: 0 },
    kpi: 'vF', expected: 1.5 * G0 * QE * ACC / HBAR, tol: 1e-4, kind: 'analytical',
    reference: 'Castro Neto et al., The electronic properties of graphene (Reviews of Modern Physics) — v_F = 3ta/2ħ'
  },
  {
    id: 's17-graphene-van-hove', suite: 17, title: 'Graphene: van Hove saddle at M, E = γ₀', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'graphene', gamma0: 2.7, tbS: 0 },
    kpi: 'E_M', expected: G0, tol: 1e-6, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — |f(M)| = 1'
  },
  {
    id: 's17-cnt-armchair-metallic', suite: 17, title: '(10,10) armchair nanotube is metallic', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'cnt', tn: 10, tm: 10, gamma0: 2.7, tbS: 0 },
    kpi: 'gap', expected: 0, abs: 1e-3, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — (n − m) mod 3 = 0 rule; all armchair tubes are metallic'
  },
  {
    id: 's17-cnt-armchair-diameter', suite: 17, title: '(10,10) nanotube diameter a√(n² + nm + m²)/π', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'cnt', tn: 10, tm: 10 },
    kpi: 'd', expected: AG * Math.sqrt(300) / PI, tol: 1e-4, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — d_t = a√(n² + nm + m²)/π (1.36 nm for (10,10))'
  },
  {
    id: 's17-cnt-zigzag-gap-exact', suite: 17, title: '(10,0) zigzag nanotube: zone-folding gap', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'cnt', tn: 10, tm: 0, gamma0: 2.7, tbS: 0 },
    kpi: 'gap', expected: 2 * G0 * Math.abs(1 + 2 * Math.cos(7 * PI / 10)), tol: 0.002, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — zigzag sub-bands E_q(0) = ±γ₀|1 + 2cos(qπ/n)|',
    note: 'Exact nearest-neighbour zone-folding value (0.948 eV); the suite searches the minimum on a discrete axial k-grid.'
  },
  {
    id: 's17-cnt-zigzag-gap-rule', suite: 17, title: '(10,0) zigzag nanotube: gap against 2a_ccγ₀/d', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'cnt', tn: 10, tm: 0, gamma0: 2.7, tbS: 0 },
    kpi: 'gap', expected: 2 * ACC * G0 / (AG * 10 / PI), tol: 0.05, kind: 'analytical',
    reference: 'Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes — E_g = 2a_ccγ₀/d_t (linear-band approximation)',
    note: 'The rule comes from the linearised Dirac cone; trigonal warping of the full π band shifts the gap of this 0.78 nm tube by about 3 %.'
  },
  {
    id: 's17-ssh-gap', suite: 17, title: 'Dimerised (SSH) chain: band gap 2|t₁ − t₂|', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'ssh', t1: 2.4, t2: 3.0 },
    kpi: 'gap', expected: 1.2, tol: 1e-6, kind: 'analytical',
    reference: 'Su, Schrieffer & Heeger model of polyacetylene — E(k) = ±√(t₁² + t₂² + 2t₁t₂cos ka)'
  },
  {
    id: 's17-ssh-edge-states', suite: 17, title: 'SSH chain ending on weak bonds: two zero-energy edge states', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'ssh', t1: 2.4, t2: 3.0 },
    kpi: 'edge', expected: 2, abs: 0, kind: 'analytical',
    reference: 'Bulk–boundary correspondence of the Su–Schrieffer–Heeger model (winding number 1 for t₂ > t₁): one mid-gap state per end'
  },
  {
    id: 's17-brus-cdse-4nm', suite: 17, title: 'CdSe quantum dot, 4 nm: Brus equation', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'brus', qdMat: 'CdSe', qdD: 4e-9 },
    kpi: 'Eg_brus', expected: 1.74 + HB2M * PI * PI / 4 * (1 / 0.13 + 1 / 0.45) - 1.786 * COUL / (10.6 * 2), tol: 1e-3, kind: 'analytical',
    reference: 'Brus effective-mass model: E = E_g + (ħ²π²/2R²)(1/m_e + 1/m_h) − 1.786e²/(4πε₀ε_rR), with CdSe E_g = 1.74 eV, m_e = 0.13, m_h = 0.45, ε_r = 10.6'
  },
  {
    id: 's17-brus-cdse-bohr-radius', suite: 17, title: 'CdSe exciton Bohr radius', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'brus', qdMat: 'CdSe', qdD: 4e-9 },
    kpi: 'aB', expected: 5.6e-9, tol: 0.03, kind: 'reference-data',
    reference: 'Widely quoted bulk exciton Bohr radius of CdSe, about 5.6 nm (quantum-dot literature following Brus and Ekimov)',
    note: 'Hydrogenic estimate ε_r a₀ m₀/μ with handbook effective masses; quoted literature values range from 4.9 to 5.6 nm depending on the masses adopted.'
  },
  {
    id: 's17-well-infinite-level', suite: 17, title: 'Particle in a box: GaAs, 6 nm, ground level h²/8m*L²', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'well', wW: 6e-9, wMw: 0.067 },
    kpi: 'E1_inf', expected: H * H / (8 * 0.067 * M0 * 36e-18) / QE, tol: 1e-4, kind: 'analytical',
    reference: 'Infinite square well, E_n = n²h²/(8mL²) — any quantum mechanics text (e.g. Griffiths, Introduction to Quantum Mechanics)'
  },
  {
    id: 's17-well-finite-ground-state', suite: 17, title: 'Finite square well (0.3 eV, 6 nm, m* = 0.067): numerical ground state', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'well', wW: 6e-9, wV0: 0.3, wMw: 0.067, wMb: 0.067 },
    kpi: 'E1', expected: wellE1(6, 0.3, 0.067), tol: 0.01, kind: 'analytical',
    reference: 'Finite square well, even-parity matching condition k tan(kW/2) = κ — Griffiths, Introduction to Quantum Mechanics',
    note: 'The suite solves the Schrödinger equation by finite differences on about 100 grid points; second-order discretisation error is a few tenths of a per cent.'
  },
  {
    id: 's17-well-bound-state-count', suite: 17, title: 'Finite square well: number of bound states', mode: 'tb',
    inputs: { mode: 'tb', tbModel: 'well', wW: 6e-9, wV0: 0.3, wMw: 0.067, wMb: 0.067 },
    kpi: 'n_bound', expected: 1 + Math.floor(Math.sqrt(0.3 * 0.067 / HB2M) * 6 / PI), abs: 0, kind: 'analytical',
    reference: 'Finite square well: N = 1 + ⌊√(2mV₀)W/(πħ)⌋ bound states — Griffiths, Introduction to Quantum Mechanics'
  },
  // ── Monte Carlo ──
  {
    id: 's17-ising-onsager-tc-kelvin', suite: 17, title: '2-D Ising model: Onsager critical temperature for J = 10 meV', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'ising', Jmev: 10 },
    kpi: 'Tc_K', expected: 2 / Math.log(1 + Math.SQRT2) * 0.010 / KB_EV, tol: 1e-4, kind: 'analytical',
    reference: 'Onsager exact solution of the square-lattice Ising model: k_BT_c/J = 2/ln(1 + √2) = 2.269'
  },
  {
    id: 's17-ising-susceptibility-peak', suite: 17, title: '2-D Ising Monte Carlo, 16×16: susceptibility peak near T_c', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'ising', mcL: 16, mcSweeps: 600, mcH: 0, seed: 1 },
    kpi: 'Tc_chi', expected: 2 / Math.log(1 + Math.SQRT2), tol: 0.07, kind: 'analytical',
    reference: 'Onsager exact solution of the square-lattice Ising model: k_BT_c/J = 2/ln(1 + √2) = 2.269',
    note: 'The peak is located on a temperature grid of spacing 0.15 J/k_B and a 16×16 lattice shifts the susceptibility maximum above the infinite-lattice T_c by an amount ∝ 1/L; both effects are documented by the suite.'
  },
  {
    id: 's17-ising-spontaneous-magnetisation', suite: 17, title: '2-D Ising Monte Carlo: magnetisation at k_BT/J = 2', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'ising', mcL: 16, mcSweeps: 600, mcT: 2.0, mcH: 0, seed: 1 },
    kpi: 'm', expected: Math.pow(1 - Math.pow(Math.sinh(1), -4), 0.125), tol: 0.02, kind: 'analytical',
    reference: 'Onsager–Yang spontaneous magnetisation of the square-lattice Ising model: m = [1 − sinh⁻⁴(2J/k_BT)]^(1/8)',
    note: 'Finite 16×16 lattice and a few hundred correlated samples: statistical and finite-size error of about 1 %.'
  },
  {
    id: 's17-latticegas-langmuir-coverage', suite: 17, title: 'Lattice gas without interactions: Langmuir coverage', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'latticegas', wmev: 0, Eads: 0.5, Tads: 300, mAds: 28, Pgas: 1e3, mcL: 16, mcSweeps: 600, seed: 1 },
    kpi: 'theta', expected: K300 * 1e3 / (1 + K300 * 1e3), abs: 0.01, kind: 'analytical',
    reference: 'Langmuir isotherm θ = KP/(1 + KP) with K = exp(E_ads/k_BT) Λ³/k_BT for immobile adsorption from a monatomic ideal gas (Hill, An Introduction to Statistical Thermodynamics)',
    note: 'Grand-canonical Monte Carlo on 256 sites; the statistical error of the mean coverage is a few 10⁻³.'
  },
  {
    id: 's17-latticegas-half-coverage-pressure', suite: 17, title: 'Lattice gas without interactions: half-coverage pressure 1/K', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'latticegas', wmev: 0, Eads: 0.5, Tads: 300, mAds: 28, mcL: 16, mcSweeps: 600, seed: 1 },
    kpi: 'logP_half', expected: -Math.log10(K300), abs: 0.03, kind: 'analytical',
    reference: 'Langmuir isotherm: θ = ½ at P = 1/K (Hill, An Introduction to Statistical Thermodynamics)',
    note: 'Interpolated between simulated pressure points one unit of ln(KP) apart; the tolerance is in decades of pressure.'
  },
  {
    id: 's17-latticegas-critical-temperature', suite: 17, title: 'Square lattice gas: critical temperature for w = 20 meV', mode: 'mc',
    inputs: { mode: 'mc', mcModel: 'latticegas', wmev: 20 },
    kpi: 'Tc', expected: 0.020 / 4 * 2 / Math.log(1 + Math.SQRT2) / KB_EV, tol: 1e-4, kind: 'analytical',
    reference: 'Lattice-gas ↔ Ising mapping (J = w/4) with the Onsager critical temperature — Hill, Statistical Mechanics; Lee & Yang'
  },
  // ── kinetic Monte Carlo ──
  {
    id: 's17-kmc-adatom-diffusivity', suite: 17, title: 'Adatom terrace diffusion coefficient ν₀a²/4 · exp(−E_s/k_BT)', mode: 'kmc',
    inputs: { mode: 'kmc', T_C: 227, Es: 0.75, nu0: 1e13, aLat: 2.9e-10 },
    kpi: 'D', expected: 1e13 * Math.exp(-0.75 / (KB_EV * 500.15)) * 2.9e-10 ** 2 / 4, tol: 1e-4, kind: 'analytical',
    reference: 'Random-walk diffusion on a square lattice, D = Γa²/4 with Arrhenius hop rate Γ = ν₀exp(−E_s/k_BT)'
  },
  {
    id: 's17-kmc-random-deposition-roughness', suite: 17, title: 'Deposition without surface diffusion: Poisson roughness W = √θ', mode: 'kmc',
    inputs: { mode: 'kmc', T_C: 27, Es: 2.5, En: 0.25, flux: 1, kCov: 2, kL: 48, seed: 1 },
    kpi: 'W', expected: Math.SQRT2, tol: 0.05, kind: 'analytical',
    reference: 'Random deposition: column heights are Poisson distributed, W² = θ (Barabási & Stanley, Fractal Concepts in Surface Growth)',
    note: 'With a 2.5 eV barrier at 300 K no hop occurs; the sample standard deviation over 2304 columns carries a statistical error of about 1.6 %.'
  },
  {
    id: 's17-kmc-deposition-time', suite: 17, title: 'Deposition of 2 ML at 1 ML/s takes 2 s', mode: 'kmc',
    inputs: { mode: 'kmc', T_C: 27, Es: 2.5, flux: 1, kCov: 2, kL: 48, seed: 1 },
    kpi: 't', expected: 2, tol: 0.05, kind: 'conservation',
    reference: 'Residence-time (Bortz–Kalos–Lebowitz) algorithm: mean simulated time equals deposited coverage divided by flux',
    note: 'Sum of 4608 exponential waiting times: relative statistical error 1/√4608 ≈ 1.5 %.'
  },
  {
    id: 's17-kmc-island-density-exponent', suite: 17, title: 'Sub-monolayer nucleation: island density ∝ (D/F)^−1/3', mode: 'kmc',
    inputs: { mode: 'kmc', T_C: 227, Es: 0.75, En: 0.25, Ees: 0, flux: 1, kCov: 2, kL: 32, seed: 1 },
    kpi: 'chi', expected: 1 / 3, tol: 0.10, kind: 'numerical-reference',
    reference: 'Rate-equation nucleation theory (Venables): χ = i*/(i* + 2) = 1/3 for a critical nucleus of one atom',
    note: 'Six flux values on a 32×32 lattice (a handful of islands each) give counting noise of several per cent in the fitted slope; at E_n/k_BT ≈ 5.8 dimer dissociation is slow but not strictly absent.'
  },
  // ── reaction pathways ──
  {
    id: 's17-neb-muller-brown-main-barrier', suite: 17, title: 'Müller–Brown surface: barrier from the deep minimum (saddle −40.66)', mode: 'neb',
    inputs: { mode: 'neb', pes: 'mb', nebPath: 'AB', nImg: 15, pesE: 0.005 },
    kpi: 'Ea', expected: (MB_S1 - MB_A) * 0.005, tol: 0.002, kind: 'numerical-reference',
    reference: 'Müller & Brown model potential (Theoretica Chimica Acta): minimum −146.700 at (−0.558, 1.442), saddle −40.665 at (−0.822, 0.624)'
  },
  {
    id: 's17-neb-muller-brown-second-barrier', suite: 17, title: 'Müller–Brown surface: barrier from the intermediate (saddle −72.25)', mode: 'neb',
    inputs: { mode: 'neb', pes: 'mb', nebPath: 'BC', nImg: 15, pesE: 0.005 },
    kpi: 'Ea', expected: (MB_S2 - MB_B) * 0.005, tol: 0.005, kind: 'numerical-reference',
    reference: 'Müller & Brown model potential (Theoretica Chimica Acta): minimum −80.768 at (−0.050, 0.467), saddle −72.249 at (0.212, 0.293)'
  },
  {
    id: 's17-neb-muller-brown-reaction-energy', suite: 17, title: 'Müller–Brown surface: energy difference between minima A and C', mode: 'neb',
    inputs: { mode: 'neb', pes: 'mb', nebPath: 'AC', nImg: 15, pesE: 0.005 },
    kpi: 'dE', expected: (MB_C - MB_A) * 0.005, tol: 0.002, kind: 'numerical-reference',
    reference: 'Müller & Brown model potential (Theoretica Chimica Acta): minima −146.700 and −108.167'
  },
  {
    id: 's17-neb-cosine-hop-barrier', suite: 17, title: 'Adatom hop on a cosine corrugation: bridge-site barrier', mode: 'neb',
    inputs: { mode: 'neb', pes: 'hop', nImg: 15, pesE: 0.005 },
    kpi: 'Ea', expected: 100 * 0.005, tol: 0.002, kind: 'analytical',
    reference: 'V = 50[2 − cos 2πx − cos 2πy]: V(½, 0) − V(0, 0) = 100 surface units (closed form)'
  },
  {
    id: 's17-neb-cosine-hop-vineyard', suite: 17, title: 'Adatom hop on a cosine corrugation: Vineyard attempt frequency', mode: 'neb',
    inputs: { mode: 'neb', pes: 'hop', nImg: 15, pesE: 0.005, pesL: 1, nebMass: 63.546 },
    kpi: 'nu', expected: Math.sqrt(50 * 4 * PI * PI * 0.005 * QE / (63.546 * U * 1e-20)) / (2 * PI), tol: 0.002, kind: 'analytical',
    reference: 'Vineyard harmonic transition-state theory, ν* = Πν_min/Π′ν_saddle; for this separable surface ν* = (1/2π)√(V″/m) with V″ = 50(2π)² at the hollow site',
    note: 'The suite evaluates the Hessians by central finite differences.'
  },
  {
    id: 's17-neb-cosine-hop-diffusivity', suite: 17, title: 'Adatom hop: surface diffusion coefficient kℓ²/4 at 500 K', mode: 'neb',
    inputs: { mode: 'neb', pes: 'hop', nImg: 15, pesE: 0.005, pesL: 1, nebMass: 63.546, dim: 2, T_C: 226.85 },
    kpi: 'D', expected: Math.sqrt(50 * 4 * PI * PI * 0.005 * QE / (63.546 * U * 1e-20)) / (2 * PI) * Math.exp(-0.5 / (KB_EV * 500)) * 1e-20 / 4, tol: 0.01, kind: 'analytical',
    reference: 'Harmonic transition-state theory k = ν*exp(−E_a/k_BT) with random-walk diffusion D = kℓ²/(2d), d = 2',
    note: 'The exponential amplifies the residual NEB barrier error: 0.1 % in E_a is about 1 % in D at E_a/k_BT ≈ 11.6.'
  },
  // ── input decks ──
  {
    id: 's17-deck-zno-cell-volume', suite: 17, title: 'Wurtzite ZnO: primitive-cell volume (√3/2)a²c', mode: 'deck',
    inputs: { mode: 'deck', dkCrystal: 'ZnO_wurtzite', dkRep: 1 },
    kpi: 'vol', expected: Math.sqrt(3) / 2 * 3.2495 ** 2 * 5.2069, tol: 0.002, kind: 'reference-data',
    reference: 'Measured wurtzite ZnO lattice parameters a = 3.2495 Å, c = 5.2069 Å (standard powder-diffraction reference pattern for zincite)'
  },
  {
    id: 's17-deck-zno-kmesh', suite: 17, title: 'Wurtzite ZnO: k-points for 0.25 Å⁻¹ spacing (9×9×5)', mode: 'deck',
    inputs: { mode: 'deck', dkCrystal: 'ZnO_wurtzite', dkRep: 1, kspacing: 0.25 },
    kpi: 'nk', expected: Math.ceil(4 * PI / (Math.sqrt(3) * 3.2495) / 0.25) ** 2 * Math.ceil(2 * PI / 5.2069 / 0.25), abs: 0, kind: 'analytical',
    reference: 'Reciprocal lattice of a hexagonal cell: |b₁| = |b₂| = 4π/(√3 a), |b₃| = 2π/c; n_i = ⌈|b_i|/Δk⌉'
  },
  {
    id: 's17-deck-zno-plane-waves', suite: 17, title: 'Wurtzite ZnO: plane waves per band at 60 Ry', mode: 'deck',
    inputs: { mode: 'deck', dkCrystal: 'ZnO_wurtzite', dkRep: 1, ecut: 60 },
    kpi: 'npw', expected: Math.sqrt(3) / 2 * 3.2495 ** 2 * 5.2069 / 0.529177210903 ** 3 * Math.pow(60, 1.5) / (6 * PI * PI), tol: 0.003, kind: 'analytical',
    reference: 'Free-electron sphere count N_pw = Ω k_cut³/(6π²) with k_cut² = E_cut in Rydberg atomic units (Martin, Electronic Structure)'
  },
];
