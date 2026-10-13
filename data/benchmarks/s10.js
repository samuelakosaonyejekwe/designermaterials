// Suite 10 — Heat transfer and thermal engineering: independent reference benchmarks.
// Every expected value is a literal or is evaluated here from a closed-form textbook expression with
// standard constants; nothing is taken from the suite's own output.

const SB = 5.670374419e-8;                 // Stefan–Boltzmann constant, W/m²/K⁴ (exact SI value)

// --- transient slab: shared property set (water-like medium, the suite's reference slab)
const k = 0.6, rho = 1000, cp = 4186, Lz = 0.01, h = 50, sig = 0.25, Ef = 1000;
const qJ = sig * Ef * Ef;                  // Joule power density σE², W/m³
const tauSlab = rho * cp * Lz / (2 * h);   // lumped slab time constant ρc_pL/2h, s

// centre temperature of a slab (full thickness L) with both faces stepped to T_s:
// θ_c = (4/π) Σ (−1)ⁿ/(2n+1) exp[−(2n+1)² π² α t / L²]
const slabCentre = (Fo) => { let s = 0; for (let n = 0; n < 60; n++) s += (n % 2 ? -1 : 1) / (2 * n + 1) * Math.exp(-((2 * n + 1) ** 2) * Math.PI ** 2 * Fo); return 4 / Math.PI * s; };
const aS = 1e-6, LS = 0.02, tS = 60, FoS = aS * tS / (LS * LS);            // k = 1, ρc_p = 10⁶ → α = 10⁻⁶ m²/s

// --- thermal metamaterial: neutral coated cylinder (shell k_s on core k_0)
const k0c = 0.5, ksc = 20, ac = 0.006, cc = 0.012, pc = (ac / cc) ** 2;
const kNeutral = ksc * ((k0c + ksc) + (k0c - ksc) * pc) / ((k0c + ksc) - (k0c - ksc) * pc);
const bCl = 0.0095, cCl = 0.012, kbgCl = 2.3;

// --- thermal diode made of two bulk segments with k_A ∝ T and k_B constant (equal length, equal k at T_m)
const TmD = 300, dTD = 200, ThD = TmD + dTD / 2, TcD = TmD - dTD / 2, LD = 1e-3, kD = 1;
const TiF = -TmD + Math.sqrt(TmD * TmD + ThD * ThD + 2 * TmD * TcD), JfD = kD * (TiF - TcD) / LD;   // A hot
const TiR = -TmD + Math.sqrt(TmD * TmD + TcD * TcD + 2 * TmD * ThD), JrD = kD * (ThD - TiR) / LD;   // B hot

export default [
  // ------------------------------------------------------------------ transient slab
  {
    id: 's10-slab-adiabatic-joule-heating-rate',
    suite: 10,
    title: 'Insulated slab with uniform Joule heating: mean temperature T₀ + σE²t/(ρc_p)',
    mode: 'slab',
    inputs: { mode: 'slab', dim: '1d', bcA: 'insulated', bcB: 'insulated', alS: 0, qR: 0, sig, Ef, k, rho, cp, Lz, T0: 298, Tamb: 298, tEnd: 100 },
    kpi: 'Tmean',
    expected: 298 + qJ * 100 / (rho * cp),
    tol: 1e-4,
    kind: 'analytical',
    reference: 'First law for an adiabatic body with a uniform volumetric source (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)',
    note: 'Reference slab, σ = 0.25 S/m and E = 1000 V/m: σE²/(ρc_p) = 0.0597 K/s, a rise of 5.97 K in 100 s. The tolerance of 10⁻⁴ of the absolute temperature is 0.5 % of the rise.'
  },
  {
    id: 's10-slab-joule-power-density',
    suite: 10,
    title: 'Joule power density σE² of the reference slab',
    mode: 'slab',
    inputs: { mode: 'slab', sig: 0.25, Ef: 1000 },
    kpi: 'qJ',
    expected: 2.5e5,
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Joule\'s law in local form, q = σE²'
  },
  {
    id: 's10-slab-steady-volumetric-heating',
    suite: 10,
    title: 'Steady slab with uniform heating and convective faces: T_max = T∞ + qL/2h + qL²/8k',
    mode: 'slab',
    inputs: { mode: 'slab', dim: '1d', bcA: 'robin', bcB: 'robin', alS: 0, qR: 0, sig, Ef, k, rho, cp, Lz, h, T0: 298, Tamb: 298, tEnd: 30 * tauSlab },
    kpi: 'Tmax',
    expected: 298 + qJ * Lz / (2 * h) + qJ * Lz * Lz / (8 * k),
    tol: 2e-4,
    kind: 'analytical',
    reference: 'Plane wall with uniform heat generation and symmetric convection (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)',
    note: 'Run for 30 time constants. The rise above ambient is 30.2 K, so 2×10⁻⁴ of the absolute value is 0.2 % of the rise; the hottest cell centre lies half a cell from the mid-plane.'
  },
  {
    id: 's10-slab-2d-steady-reduces-to-1d',
    suite: 10,
    title: '2-D cross-section, full electrode coverage and a very wide slab: 1-D steady maximum recovered',
    mode: 'slab',
    inputs: { mode: 'slab', dim: '2d', hfrac: 1, Lxw: 50, bcA: 'robin', bcB: 'robin', alS: 0, qR: 0, sig, Ef, k, rho, cp, Lz, h, T0: 298, Tamb: 298, tEnd: 30 * tauSlab },
    kpi: 'Tmax',
    expected: 298 + qJ * Lz / (2 * h) + qJ * Lz * Lz / (8 * k),
    tol: 3e-4,
    kind: 'analytical',
    reference: 'Plane wall with uniform heat generation and symmetric convection (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)',
    note: 'Edge cooling only affects a margin a few thicknesses wide, so the interior of a 50 m wide, 10 mm thick slab is one-dimensional. 16 cells through the thickness.'
  },
  {
    id: 's10-slab-fourier-series-centre-temperature',
    suite: 10,
    title: 'Slab with both faces stepped to a fixed temperature: centre temperature from the Fourier series (Fo = 0.15)',
    mode: 'slab',
    inputs: { mode: 'slab', dim: '1d', bcA: 'fixed', bcB: 'fixed', Ef: 0, qR: 0, alS: 0, k: 1, rho: 1000, cp: 1000, Lz: LS, T0: 1000, Tamb: 200, tEnd: tS },
    kpi: 'Tmax',
    expected: 200 + 800 * slabCentre(FoS),
    tol: 0.005,
    kind: 'analytical',
    reference: 'Carslaw & Jaeger, Conduction of Heat in Solids — slab with prescribed surface temperature',
    note: 'θ_c = (4/π)Σ(−1)ⁿ/(2n+1)·exp[−(2n+1)²π²αt/L²], evaluated here. Backward Euler with 240 steps and 60 cells; the deviation is the first-order time-discretisation error.'
  },
  {
    id: 's10-slab-semi-infinite-heat-uptake',
    suite: 10,
    title: 'Semi-infinite solid with a stepped surface temperature: absorbed heat 2kΔT√(t/πα)',
    mode: 'slab',
    inputs: { mode: 'slab', dim: '1d', bcA: 'fixed', bcB: 'insulated', Ef: 0, qR: 0, alS: 0, k: 1, rho: 1000, cp: 1000, Lz: 0.05, T0: 300, Tamb: 1300, tEnd: 100 },
    kpi: 'Tmean',
    expected: 300 + 2 * 1 * 1000 * Math.sqrt(100 / (Math.PI * 1e-6)) / (1000 * 1000 * 0.05),
    tol: 0.01,
    kind: 'analytical',
    reference: 'Carslaw & Jaeger, Conduction of Heat in Solids — semi-infinite solid, error-function solution',
    note: 'Thermal penetration √(αt) = 10 mm in a 50 mm slab, so the far face is not yet felt (erfc(2.5) ≈ 4×10⁻⁴). The mean temperature rise equals the absorbed heat divided by ρc_pL. The step change at t = 0 is the hardest case for backward Euler, hence 1 %.'
  },
  {
    id: 's10-slab-energy-balance-closure',
    suite: 10,
    title: 'Energy conservation of the nonlinear reference run (stored = generated − lost)',
    mode: 'slab',
    inputs: { mode: 'slab' },
    kpi: 'eb',
    expected: 0,
    abs: 1e-8,
    kind: 'conservation',
    reference: 'First law of thermodynamics applied to the whole slab',
    note: 'Default reference slab with temperature-dependent Joule and Arrhenius reaction heating; relative residual of the finite-volume solution.'
  },
  // ------------------------------------------------------------------ lumped body, wall, fin
  {
    id: 's10-lumped-time-constant',
    suite: 10,
    title: 'Lumped-capacitance time constant ρc_pV/(hA) (bare vessel)',
    mode: 'lumped',
    inputs: { mode: 'lumped', Lw: 0, Lins: 0, h: 10, rho: 1000, cp: 4186, Vl: 1e-3, Al: 0.06 },
    kpi: 'tau',
    expected: 1000 * 4186 * 1e-3 / (10 * 0.06),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Lumped-capacitance method (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)'
  },
  {
    id: 's10-lumped-composite-wall-U',
    suite: 10,
    title: 'Overall coefficient of a plane composite wall: 1/U = L_w/k_w + L_ins/k_ins + 1/h',
    mode: 'lumped',
    inputs: { mode: 'lumped', Lw: 0.005, kw: 16, Lins: 0.02, kins: 0.04, h: 50 },
    kpi: 'U',
    expected: 1 / (0.005 / 16 + 0.02 / 0.04 + 1 / 50),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Series thermal resistances of a plane composite wall (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)'
  },
  {
    id: 's10-lumped-newton-cooling-time',
    suite: 10,
    title: 'Newtonian cooling from 363 K to 313 K in 298 K surroundings: t = τ ln[(T_s − T∞)/(T_c − T∞)]',
    mode: 'lumped',
    inputs: { mode: 'lumped', Lw: 0, Lins: 0, h: 10, rho: 1000, cp: 4186, Vl: 1e-3, Al: 0.06, Tsoak: 363, Tcool: 313, Tamb: 298 },
    kpi: 'tCool',
    expected: 1000 * 4186 * 1e-3 / (10 * 0.06) * Math.log((363 - 298) / (313 - 298)),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Newton\'s law of cooling, lumped-capacitance solution'
  },
  {
    id: 's10-lumped-pin-fin-efficiency',
    suite: 10,
    title: 'Pin fin with adiabatic tip: efficiency tanh(mL)/(mL), m = √(4h/kd)',
    mode: 'lumped',
    inputs: { mode: 'lumped', h: 50, dfin: 0.005, Lfin: 0.05, kfin: 200 },
    kpi: 'etaF',
    expected: Math.tanh(0.05 * Math.sqrt(4 * 50 / (200 * 0.005))) / (0.05 * Math.sqrt(4 * 50 / (200 * 0.005))),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Extended-surface (fin) theory, adiabatic-tip solution (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)'
  },
  // ------------------------------------------------------------------ radiative furnace
  {
    id: 's10-furnace-stefan-boltzmann-equilibrium',
    suite: 10,
    title: 'Black wafer between a black heater (1000 K) and a black wall (300 K): T⁴ = (T_h⁴ + T_w⁴)/2',
    mode: 'furnace',
    inputs: { mode: 'furnace', Th: 1000, Twall: 300, epsW: 1, epsH: 1, hg: 0 },
    kpi: 'Tss',
    expected: Math.pow((1000 ** 4 + 300 ** 4) / 2, 0.25),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Stefan–Boltzmann law; radiative equilibrium of a thin plate between two black surfaces',
    note: 'No convection. 842.6 K.'
  },
  {
    id: 's10-furnace-stefan-boltzmann-flux',
    suite: 10,
    title: 'Net radiative flux through the black wafer at equilibrium: σ(T_h⁴ − T_w⁴)/2',
    mode: 'furnace',
    inputs: { mode: 'furnace', Th: 1000, Twall: 300, epsW: 1, epsH: 1, hg: 0 },
    kpi: 'qin',
    expected: SB * (1000 ** 4 - 300 ** 4) / 2 / 1e3,
    tol: 1e-5,
    kind: 'analytical',
    reference: 'Stefan–Boltzmann law with σ = 5.670374×10⁻⁸ W/m²/K⁴; a black radiation shield halves the exchange',
    note: 'KPI in kW/m²: 28.12 kW/m².'
  },
  {
    id: 's10-furnace-grey-parallel-plates',
    suite: 10,
    title: 'Grey wafer, radiation only: T⁴ = (ε_eff T_h⁴ + ε_w T_w⁴)/(ε_eff + ε_w), ε_eff = 1/(1/ε_w + 1/ε_h − 1)',
    mode: 'furnace',
    inputs: { mode: 'furnace', Th: 1273, Twall: 300, epsW: 0.65, epsH: 0.85, hg: 0 },
    kpi: 'Tss',
    expected: (() => { const e = 1 / (1 / 0.65 + 1 / 0.85 - 1); return Math.pow((e * 1273 ** 4 + 0.65 * 300 ** 4) / (e + 0.65), 0.25); })(),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Grey-body exchange between infinite parallel plates and with a large enclosure (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)'
  },
  // ------------------------------------------------------------------ nanoscale conduction
  {
    id: 's10-nano-kinetic-mean-free-path',
    suite: 10,
    title: 'Grey phonon mean free path of silicon from kinetic theory, Λ = 3k/(C_v v)',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'grey', kb: 148, Cv: 1.66e6, vs: 6400 },
    kpi: 'Lam',
    expected: 3 * 148 / (1.66e6 * 6400) * 1e9,
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Kinetic theory of thermal conductivity, k = C v Λ/3 (Kittel, Introduction to Solid State Physics)',
    note: 'About 42 nm with the stated bulk values; the real spectrum-averaged value for silicon is several times longer, which the suite states.'
  },
  {
    id: 's10-nano-fuchs-sondheimer-thick-film',
    suite: 10,
    title: 'Fuchs–Sondheimer thin film, diffuse surfaces, thick-film limit: k/k_bulk = 1 − 3Λ/(8d)',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'grey', kb: 100, Cv: 1e6, vs: 3000, dch: 1e-6, spec: 0 },
    kpi: 'kFilm',
    expected: 100 * (1 - 3 / 80),
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Fuchs–Sondheimer theory of size effects in thin films (Ziman, Electrons and Phonons)',
    note: 'Λ = 100 nm, d = 10Λ. The terms neglected in the asymptote are of order e^(−10).'
  },
  {
    id: 's10-nano-maxwell-effective-conductivity',
    suite: 10,
    title: 'Maxwell (Maxwell-Garnett) conductivity of a dilute sphere dispersion, perfect interfaces',
    mode: 'nano',
    inputs: { mode: 'nano', kpar: 50, kmat: 0.6, phiP: 0.05 },
    kpi: 'kMax',
    expected: 0.6 * (50 + 2 * 0.6 + 2 * 0.05 * (50 - 0.6)) / (50 + 2 * 0.6 - 0.05 * (50 - 0.6)),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Maxwell, A Treatise on Electricity and Magnetism — effective conductivity of a dispersion of spheres'
  },
  {
    id: 's10-nano-kapitza-critical-radius',
    suite: 10,
    title: 'Interface resistance cancels the particle benefit at a = a_K/(1 − k_m/k_p): composite equals matrix',
    mode: 'nano',
    inputs: { mode: 'nano', RK: 1e-8, kmat: 0.6, kpar: 50, phiP: 0.2, aP: 1e-8 * 0.6 / (1 - 0.6 / 50) },
    kpi: 'kHJ',
    expected: 0.6,
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Hasselman & Johnson effective conductivity of composites with interfacial thermal barrier resistance; Kapitza radius a_K = R_K k_m',
    note: 'At the critical radius the particle with its interface is thermally equivalent to the matrix for any volume fraction.'
  },
  // ------------------------------------------------------------------ thermal metamaterials
  {
    id: 's10-meta-laminate-wiener-bounds',
    suite: 10,
    title: 'Cu/polymer laminate anisotropy: k_∥/k_⊥ from the arithmetic and harmonic means',
    mode: 'meta',
    inputs: { mode: 'meta', kA: 400, kBl: 0.15, fA: 0.5 },
    kpi: 'aniso',
    expected: (0.5 * 400 + 0.5 * 0.15) * (0.5 / 400 + 0.5 / 0.15),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Wiener bounds: parallel and series conduction in a layered medium'
  },
  {
    id: 's10-meta-bilayer-cloak-condition',
    suite: 10,
    title: 'Bilayer thermal cloak: with k₂ = k_b(c² + b²)/(c² − b²) and an insulating inner layer the external field is undisturbed',
    mode: 'meta',
    inputs: { mode: 'meta', kbg: kbgCl, kin: 0.001, kout: kbgCl * (cCl * cCl + bCl * bCl) / (cCl * cCl - bCl * bCl), kcore: 200, ra: 0.006, rb: bCl, rc: cCl },
    kpi: 'dist',
    expected: 0,
    abs: 1e-3,
    kind: 'analytical',
    reference: 'Exact bilayer thermal cloak design condition for a 2-D (cylindrical) cloak with a perfectly insulating inner layer (Han et al., Physical Review Letters)',
    note: 'The inner layer is 0.001 W/m/K rather than exactly zero (lower bound of the input), which leaves a residual dipole of order k_inner/k_outer; hence an absolute tolerance of 10⁻³ of the applied gradient.'
  },
  {
    id: 's10-meta-neutral-coated-cylinder',
    suite: 10,
    title: 'Neutral coated cylinder: background set to the Hashin coated-cylinder conductivity gives zero external disturbance',
    mode: 'meta',
    inputs: { mode: 'meta', kcore: k0c, kin: ksc, kout: ksc, ra: ac, rb: 0.009, rc: cc, kbg: kNeutral },
    kpi: 'dist',
    expected: 0,
    abs: 1e-9,
    kind: 'analytical',
    reference: 'Hashin coated-cylinder assemblage / neutral inclusion (Milton, The Theory of Composites)',
    note: 'Both shell layers are given the same conductivity so the shell is homogeneous; k_b = k_s[(k₀+k_s)+(k₀−k_s)(a/c)²]/[(k₀+k_s)−(k₀−k_s)(a/c)²] is evaluated here.'
  },
  {
    id: 's10-meta-bare-inclusion-dipole',
    suite: 10,
    title: 'Uncloaked cylinder: dipole strength at r = c is |(k_b − k₀)/(k_b + k₀)|(a/c)²',
    mode: 'meta',
    inputs: { mode: 'meta', kbg: 2.3, kcore: 200, ra: 0.006, rc: 0.012 },
    kpi: 'bare',
    expected: (200 - 2.3) / (200 + 2.3) * 0.25,
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Circular inclusion in a uniform gradient (2-D Laplace equation; Carslaw & Jaeger, Conduction of Heat in Solids)'
  },
  {
    id: 's10-meta-concentrator-factor',
    suite: 10,
    title: 'Radially anisotropic shell: concentration factor (c/a)^(1 − √(k_θ/k_r))',
    mode: 'meta',
    inputs: { mode: 'meta', kA: 400, kBl: 0.15, fA: 0.5, ra: 0.006, rb: 0.0095, rc: 0.012 },
    kpi: 'conc',
    expected: Math.pow(2, 1 - Math.sqrt(1 / ((0.5 * 400 + 0.5 * 0.15) * (0.5 / 400 + 0.5 / 0.15)))),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Power-law solution r^(±√(k_θ/k_r)) of the conduction equation in a cylindrically anisotropic shell',
    note: 'Radial-spoke laminate: k_r = k_∥, k_θ = k_⊥; impedance-matched exterior assumed, as stated by the suite.'
  },
  // ------------------------------------------------------------------ negative / zero thermal expansion
  {
    id: 's10-nte-turner-cte',
    suite: 10,
    title: 'Turner CTE of Cu / ZrW₂O₈ at 40 vol% filler',
    mode: 'nte',
    inputs: { mode: 'nte', aM: 17e-6, KM: 140, aN: -9e-6, KN: 72, Vf: 0.4 },
    kpi: 'aT',
    expected: (17 * 140 * 0.6 + (-9) * 72 * 0.4) / (140 * 0.6 + 72 * 0.4),
    tol: 1e-9,
    kind: 'analytical',
    reference: 'Turner model of composite thermal expansion (bulk-modulus-weighted mean)',
    note: 'KPI in ppm/K: 10.36 ppm/K.'
  },
  {
    id: 's10-nte-turner-zero-expansion-fraction',
    suite: 10,
    title: 'Filler fraction for zero expansion (Turner): V₀ = α_mK_m/(α_mK_m − α_fK_f)',
    mode: 'nte',
    inputs: { mode: 'nte', aM: 60e-6, KM: 5, aN: -9e-6, KN: 72 },
    kpi: 'V0',
    expected: 60 * 5 / (60 * 5 + 9 * 72),
    tol: 1e-8,
    kind: 'analytical',
    reference: 'Turner model of composite thermal expansion, solved for zero net expansion',
    note: 'Epoxy / ZrW₂O₈: 31.6 vol%.'
  },
  {
    id: 's10-nte-levin-hashin-shtrikman-cte',
    suite: 10,
    title: 'Schapery bound: Levin relation with the Hashin–Shtrikman bulk modulus (stiffer, higher-CTE matrix as reference phase)',
    mode: 'nte',
    inputs: { mode: 'nte', aM: 17e-6, KM: 140, GM: 48, aN: -9e-6, KN: 72, GN: 34, Vf: 0.4 },
    kpi: 'aU',
    expected: (() => { const KM = 140, GM = 48, KN = 72, v = 0.4, K = KM + v / (1 / (KN - KM) + 3 * (1 - v) / (3 * KM + 4 * GM)); return 17 + (-9 - 17) * (1 / K - 1 / KM) / (1 / KN - 1 / KM); })(),
    tol: 1e-8,
    kind: 'analytical',
    reference: 'Levin relation between effective thermal expansion and bulk modulus; Hashin–Shtrikman bounds; Schapery bounds on composite CTE',
    note: 'KPI in ppm/K. The Hashin–Shtrikman bound built on the copper matrix (the upper bulk-modulus bound) gives the upper CTE bound here.'
  },
  {
    id: 's10-nte-bimaterial-triangle-lattice',
    suite: 10,
    title: 'Pin-jointed bi-material triangle (Al base, Invar legs, 30° half angle): negative CTE of the height',
    mode: 'nte',
    inputs: { mode: 'nte', aHi: 23e-6, aLo: 1.2e-6, phiL: 30 },
    kpi: 'aLat',
    expected: (1.2 - 23 * 0.25) / 0.75,
    tol: 1e-8,
    kind: 'analytical',
    reference: 'Kinematics of a pin-jointed bi-material triangle, h² = l² − (b/2)² (Lakes; Steeves et al., low-expansion lattices)',
    note: 'α_h = (α_leg − α_base sin²φ)/cos²φ = −6.07 ppm/K: the expanding base pulls the apex in.'
  },
  // ------------------------------------------------------------------ thermal rectification
  {
    id: 's10-rect-linear-junction-series-conductance',
    suite: 10,
    title: 'Linear, symmetric junction: heat current ΔT/(L_A/k_A + 1/G_i + L_B/k_B) and no rectification',
    mode: 'rect',
    inputs: { mode: 'rect', nA: 0, nB: 0, sA: 0, sB: 0, del0: 0, kA0: 1, kB0: 2, LA: 1e-6, LB: 1e-6, Gi0: 2e5, Tm: 300, dTop: 100, dTmax: 200 },
    kpi: 'Jf',
    expected: 100 / (1e-6 / 1 + 1 / 2e5 + 1e-6 / 2),
    tol: 1e-6,
    kind: 'analytical',
    reference: 'Series thermal resistances (Fourier\'s law with a constant interface conductance)'
  },
  {
    id: 's10-rect-linear-junction-no-rectification',
    suite: 10,
    title: 'Temperature-independent junction: rectification ratio exactly one',
    mode: 'rect',
    inputs: { mode: 'rect', nA: 0, nB: 0, sA: 0, sB: 0, del0: 0, kA0: 1, kB0: 2, LA: 1e-6, LB: 1e-6, Gi0: 2e5, Tm: 300, dTop: 100, dTmax: 200 },
    kpi: 'R',
    expected: 1,
    tol: 1e-6,
    kind: 'conservation',
    reference: 'Reciprocity of linear heat conduction: rectification requires temperature-dependent (nonlinear) properties'
  },
  {
    id: 's10-rect-bulk-diode-rectification',
    suite: 10,
    title: 'Two-segment bulk thermal diode (k_A ∝ T, k_B constant): rectification ratio from the Kirchhoff-transform solution',
    mode: 'rect',
    inputs: { mode: 'rect', nA: 1, nB: 0, sA: 0, sB: 0, del0: 0, kA0: kD, kB0: kD, LA: LD, LB: LD, Gi0: 1e10, Tm: TmD, dTop: dTD, dTmax: dTD },
    kpi: 'R',
    expected: Math.max(JfD, JrD) / Math.min(JfD, JrD),
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Kirchhoff transformation for temperature-dependent conductivity (Carslaw & Jaeger, Conduction of Heat in Solids); two-segment thermal diode',
    note: 'Interface temperature from the quadratic flux-continuity condition, evaluated here. The interface conductance is set seven orders above the segment conductances so that it plays no part.'
  },
  {
    id: 's10-rect-bulk-diode-forward-current',
    suite: 10,
    title: 'Two-segment bulk thermal diode: forward heat current',
    mode: 'rect',
    inputs: { mode: 'rect', nA: 1, nB: 0, sA: 0, sB: 0, del0: 0, kA0: kD, kB0: kD, LA: LD, LB: LD, Gi0: 1e10, Tm: TmD, dTop: dTD, dTmax: dTD },
    kpi: 'Jf',
    expected: JfD,
    tol: 1e-4,
    kind: 'analytical',
    reference: 'Kirchhoff transformation for temperature-dependent conductivity (Carslaw & Jaeger, Conduction of Heat in Solids)',
    note: 'Segment A (hot side, k rising with T) and segment B of equal length; J = k_B(T_i − T_c)/L with T_i from the quadratic.'
  },
  // ------------------------------------------------------------------ published experiments
  {
    id: 's10-exp-nanofluid-al2o3-eg-3vol',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 3.1 vol% Al₂O₃ / ethylene glycol nanofluid at 303 K',
    mode: 'nano',
    inputs: { mode: 'nano', kpar: 36, kmat: 0.2463, phiP: 0.031, aP: 21.5e-9 },
    kpi: 'kMax',
    expected: 0.2652,
    tol: 0.04,
    kind: 'experiment',
    measured: { value: 0.2652, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Pastoriza-Gallego, Lugo, Legido, Piñeiro', title: 'Thermal conductivity and viscosity measurements of ethylene glycol-based Al2O3 nanofluids', journal: 'Nanoscale Res. Lett.', year: 2011, doi: '10.1186/1556-276X-6-221', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3211279/' },
    conditions: 'Al₂O₃ nanoparticles of 43 nm nominal size dispersed in ethylene glycol, volume fraction 0.031, transient hot-wire (needle probe) method, 303.15 K, atmospheric pressure',
    match: ['phiP', 'kmat', 'kpar'],
    reference: 'Pastoriza-Gallego et al., Nanoscale Res. Lett. 6, 221 (2011)',
    note: 'Table 2 of the source (303.15 K column); base-fluid conductivity 0.2463 W/m/K from Table 1 and particle conductivity 36 W/m/K as used by the authors. The source gives no uncertainty for individual points (its base-fluid data agree with the literature to 1.8 % on average). The Maxwell value is compared, so no interface resistance is needed; as the source also finds, Maxwell over-predicts the enhancement itself (9.4 % predicted against 7.7 % measured), which is a 1.6 % difference in the conductivity.'
  },
  {
    id: 's10-exp-nanofluid-al2o3-eg-5vol',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 4.8 vol% Al₂O₃ / ethylene glycol nanofluid at 303 K',
    mode: 'nano',
    inputs: { mode: 'nano', kpar: 36, kmat: 0.2463, phiP: 0.048, aP: 21.5e-9 },
    kpi: 'kMax',
    expected: 0.2773,
    tol: 0.04,
    kind: 'experiment',
    measured: { value: 0.2773, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Pastoriza-Gallego, Lugo, Legido, Piñeiro', title: 'Thermal conductivity and viscosity measurements of ethylene glycol-based Al2O3 nanofluids', journal: 'Nanoscale Res. Lett.', year: 2011, doi: '10.1186/1556-276X-6-221', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3211279/' },
    conditions: 'Al₂O₃ nanoparticles of 43 nm nominal size dispersed in ethylene glycol, volume fraction 0.048, transient hot-wire (needle probe) method, 303.15 K, atmospheric pressure',
    match: ['phiP', 'kmat', 'kpar'],
    reference: 'Pastoriza-Gallego et al., Nanoscale Res. Lett. 6, 221 (2011)',
    note: 'Table 2 of the source (303.15 K column); inputs as for the 3.1 vol% case. Maxwell over-predicts the enhancement (14.8 % against 12.6 % measured), a 2 % difference in the conductivity; the deviation grows with loading (3.7 % at 6.6 vol%).'
  },
  {
    id: 's10-exp-nte-zrw2o8-rubber-cte',
    suite: 10,
    title: 'Experiment: thermal expansion of HNBR rubber filled with 39.7 vol% ZrW₂O₈ (rubbery state)',
    mode: 'nte',
    inputs: { mode: 'nte', aM: 187e-6, KM: 2.2, GM: 0.03, aN: -9e-6, KN: 88 / (3 * (1 - 2 * 0.3)), GN: 88 / (2 * (1 + 0.3)), Vf: 0.397 },
    kpi: 'aU',
    expected: 0.54 * 187,
    tol: 0.12,
    kind: 'experiment',
    measured: { value: 54, uncertainty: null, unit: '% of the CTE of the unfilled rubber' },
    source: { authors: 'Akulichev, Alcock, Tiwari, Echtermeyer', title: 'Thermomechanical properties of zirconium tungstate / hydrogenated nitrile butadiene rubber (HNBR) composites for low-temperature applications', journal: 'J. Mater. Sci.', year: 2016, doi: '10.1007/s10853-016-0236-6', url: 'https://arxiv.org/abs/1710.05195' },
    conditions: 'Hydrogenated nitrile butadiene rubber filled with 39.7 vol% ZrW₂O₈ particles, linear thermal expansion measured between 298 K and 403 K (above the glass transition of the rubber)',
    match: ['Vf', 'aM', 'aN', 'KM', 'KN'],
    reference: 'Akulichev et al., J. Mater. Sci. (2016), author version arXiv:1710.05195',
    note: 'Section 3.2 of the source: the CTE of the most highly filled compound is 54 % of that of the unfilled rubber; with the matrix CTE of 187×10⁻⁶ K⁻¹ listed in the source (Table 2) this is 101 ppm/K, the expected value here. Constituent data from the same table: ZrW₂O₈ E = 88 GPa, ν = 0.3, α = −9×10⁻⁶ K⁻¹; rubber E = 6.6 MPa, ν = 0.4995 (K = 2.2 GPa, G = 2.2 MPa). The matrix shear modulus is entered at the lower bound of the input (0.03 GPa); the bound does not change at this level. The upper Schapery bound (equal to the Kerner estimate, which the source finds adequate) is compared. The Turner estimate fails for this soft matrix (−0.5 ppm/K) because it assumes equal strain in both phases. No uncertainty is quoted in the text (error bars are shown only graphically).'
  },
  {
    id: 's10-exp-si-nanowire-120nm',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 120 nm silicon nanowire at room temperature',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'spectral', kb: 148, vs: 6400, nat: 5.0e28, nbas: 2, Gam: 2.0e-4, Tn: 300, spec: 0, dch: 120e-9 },
    kpi: 'kWire',
    expected: 36,
    tol: 0.1,
    kind: 'experiment',
    measured: { value: 36, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Lee, J.; Lee, W.; Lim, J.; Yu, Y.; Kong, Q.; Urban, J. J.; Yang, P.', title: 'Thermal Transport in Silicon Nanowires at High Temperature up to 700 K', journal: 'Nano Lett. 16, 4133–4140', year: 2016, doi: '10.1021/acs.nanolett.6b00956', url: 'https://nanowires.berkeley.edu/wp-content/uploads/2017/01/Nano-Letters-2016-Lee.pdf' },
    conditions: 'Single-crystal silicon nanowires grown by the vapour–liquid–solid method (smooth surfaces), suspended between microfabricated heater and sensor membranes with Pt contacts; thermal contact resistance accounted for; room temperature',
    match: ['dch', 'Tn', 'kb'],
    reference: 'Lee, Lee, Lim, Yu, Kong, Urban & Yang, Nano Lett. 16, 4133 (2016)',
    note: 'Text of the source: "the thermal conductivity of VLS Si NWs with diameters of 120, 70, and 40 nm at room temperature is 36, 28, and 22 W m⁻¹ K⁻¹, respectively" (no uncertainty stated for the individual values). Suite inputs are bulk silicon data only (148 W/m/K, mean sound velocity 6400 m/s, 5.0×10²⁸ atoms/m³, two atoms per primitive cell, natural-isotope Γ = 2.0×10⁻⁴) with fully diffuse surfaces; nothing is fitted to nanowire data. Spectral Born–von Kármán model: 37.7 W/m/K (+5 %). The grey single-mean-free-path model used previously gives 110 W/m/K. Tolerance: model simplifications (isotropic single-branch dispersion, ω² Umklapp law) and the unstated measurement uncertainty.'
  },
  {
    id: 's10-exp-si-nanowire-70nm',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 70 nm silicon nanowire at room temperature',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'spectral', kb: 148, vs: 6400, nat: 5.0e28, nbas: 2, Gam: 2.0e-4, Tn: 300, spec: 0, dch: 70e-9 },
    kpi: 'kWire',
    expected: 28,
    tol: 0.1,
    kind: 'experiment',
    measured: { value: 28, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Lee, J.; Lee, W.; Lim, J.; Yu, Y.; Kong, Q.; Urban, J. J.; Yang, P.', title: 'Thermal Transport in Silicon Nanowires at High Temperature up to 700 K', journal: 'Nano Lett. 16, 4133–4140', year: 2016, doi: '10.1021/acs.nanolett.6b00956', url: 'https://nanowires.berkeley.edu/wp-content/uploads/2017/01/Nano-Letters-2016-Lee.pdf' },
    conditions: 'Single-crystal silicon nanowires grown by the vapour–liquid–solid method (smooth surfaces), suspended between microfabricated heater and sensor membranes with Pt contacts; thermal contact resistance accounted for; room temperature',
    match: ['dch', 'Tn', 'kb'],
    reference: 'Lee, Lee, Lim, Yu, Kong, Urban & Yang, Nano Lett. 16, 4133 (2016)',
    note: 'Text of the source: "the thermal conductivity of VLS Si NWs with diameters of 120, 70, and 40 nm at room temperature is 36, 28, and 22 W m⁻¹ K⁻¹, respectively" (no uncertainty stated for the individual values). Suite inputs are bulk silicon data only (148 W/m/K, mean sound velocity 6400 m/s, 5.0×10²⁸ atoms/m³, two atoms per primitive cell, natural-isotope Γ = 2.0×10⁻⁴) with fully diffuse surfaces; nothing is fitted to nanowire data. Spectral Born–von Kármán model: 27.6 W/m/K (−1 %). The grey single-mean-free-path model used previously gives 93 W/m/K. Tolerance: model simplifications (isotropic single-branch dispersion, ω² Umklapp law) and the unstated measurement uncertainty.'
  },
  {
    id: 's10-exp-si-nanowire-40nm',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 40 nm silicon nanowire at room temperature',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'spectral', kb: 148, vs: 6400, nat: 5.0e28, nbas: 2, Gam: 2.0e-4, Tn: 300, spec: 0, dch: 40e-9 },
    kpi: 'kWire',
    expected: 22,
    tol: 0.2,
    kind: 'experiment',
    measured: { value: 22, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Lee, J.; Lee, W.; Lim, J.; Yu, Y.; Kong, Q.; Urban, J. J.; Yang, P.', title: 'Thermal Transport in Silicon Nanowires at High Temperature up to 700 K', journal: 'Nano Lett. 16, 4133–4140', year: 2016, doi: '10.1021/acs.nanolett.6b00956', url: 'https://nanowires.berkeley.edu/wp-content/uploads/2017/01/Nano-Letters-2016-Lee.pdf' },
    conditions: 'Single-crystal silicon nanowires grown by the vapour–liquid–solid method (smooth surfaces), suspended between microfabricated heater and sensor membranes with Pt contacts; thermal contact resistance accounted for; room temperature',
    match: ['dch', 'Tn', 'kb'],
    reference: 'Lee, Lee, Lim, Yu, Kong, Urban & Yang, Nano Lett. 16, 4133 (2016)',
    note: 'Text of the source: "the thermal conductivity of VLS Si NWs with diameters of 120, 70, and 40 nm at room temperature is 36, 28, and 22 W m⁻¹ K⁻¹, respectively" (no uncertainty stated for the individual values). Suite inputs are bulk silicon data only (148 W/m/K, mean sound velocity 6400 m/s, 5.0×10²⁸ atoms/m³, two atoms per primitive cell, natural-isotope Γ = 2.0×10⁻⁴) with fully diffuse surfaces; nothing is fitted to nanowire data. Spectral Born–von Kármán model: 19.0 W/m/K (−14 %). The grey single-mean-free-path model used previously gives 72 W/m/K. Tolerance: model simplifications (isotropic single-branch dispersion, ω² Umklapp law) and the unstated measurement uncertainty.'
  },
  {
    id: 's10-exp-si-nanowire-115nm-li',
    suite: 10,
    title: 'Experiment: thermal conductivity of a 115 nm silicon nanowire at 300 K (first suspended-device measurements)',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'spectral', kb: 148, vs: 6400, nat: 5.0e28, nbas: 2, Gam: 2.0e-4, Tn: 300, spec: 0, dch: 115e-9 },
    kpi: 'kWire',
    expected: 40,
    tol: 0.15,
    kind: 'experiment',
    measured: { value: 40, uncertainty: null, unit: 'W/m/K' },
    source: { authors: 'Li, D.; Wu, Y.; Kim, P.; Shi, L.; Yang, P.; Majumdar, A.', title: 'Thermal conductivity of individual silicon nanowires', journal: 'Appl. Phys. Lett. 83, 2934–2936', year: 2003, doi: '10.1063/1.1616981', url: 'https://doi.org/10.1063/1.1616981' },
    conditions: 'Individual single-crystal intrinsic silicon nanowires (vapour–liquid–solid growth) of 22, 37, 56 and 115 nm diameter on a microfabricated suspended device, 20–320 K',
    match: ['dch', 'Tn', 'kb'],
    reference: 'Li, Wu, Kim, Shi, Yang & Majumdar, Appl. Phys. Lett. 83, 2934 (2003); 300 K values as summarised by Lee et al., Nano Lett. 16, 4133 (2016)',
    note: 'The original article gives its results as curves; the number used here is the summary in Lee et al. (2016): for diameters between 22 and 115 nm "the values at 300 K were ranging between 7 and 40 W m⁻¹ K⁻¹", so about 40 W/m/K belongs to the 115 nm wire (read to roughly ±10 %). Model 36.9 W/m/K (−8 %); inputs as for the other nanowire entries. The lower end of the same range, about 7 W/m/K for the 22 nm wire, is NOT reproduced: the model gives 12.2 W/m/K (+74 %). Diffuse-boundary models with bulk dispersion are known to over-predict that wire; the suite warns below 30 nm and the point is not offered as a validation entry.'
  },
  {
    id: 's10-exp-si-film-260nm-in-plane',
    suite: 10,
    title: 'Experiment: in-plane thermal conductivity of a 260 nm silicon-on-insulator device layer',
    mode: 'nano',
    inputs: { mode: 'nano', kmodel: 'spectral', kb: 148, vs: 6400, nat: 5.0e28, nbas: 2, Gam: 2.0e-4, Tn: 300, spec: 0, dch: 260e-9 },
    kpi: 'kFilm',
    expected: 90,
    tol: 0.2,
    kind: 'experiment',
    measured: { value: 90, uncertainty: 10, unit: 'W/m/K' },
    source: { authors: 'Aubain, M. S.; Bandaru, P. R.', title: 'In-plane Thermal Conductivity Determination in Silicon on Insulator (SOI) Structures Through Thermoreflectance Measurements', journal: 'MRS Online Proceedings Library 1267, DD12-01', year: 2010, doi: '10.1557/PROC-1267-DD12-01', url: 'https://doi.org/10.1557/PROC-1267-DD12-01' },
    conditions: '260 nm single-crystal silicon device layer of an SOI wafer; temperature field of microfabricated a.c. Joule heaters mapped by scanning thermoreflectance and fitted with a numerical heat-conduction model; the abstract does not state the temperature (taken as 300 K)',
    match: ['dch', 'kb'],
    reference: 'Aubain & Bandaru, MRS Online Proc. Libr. 1267, DD12-01 (2010)',
    note: 'Abstract of the source: the in-plane conductivity of the device layer "is reduced to 90±10 W/(m.K)". The spectral model with the Fuchs–Sondheimer factor applied to every mode and fully diffuse surfaces gives 74.2 W/m/K (−18 %, 1.6 times the quoted uncertainty); the grey model gives 139 W/m/K (+54 %). The model is on the low side for films: buried-oxide and native-oxide interfaces of SOI layers are partly specular for long-wavelength phonons (a specularity of 0.4 would reproduce the value, but none was fitted). Bond-and-etch-back SOI layers of 0.42 µm were reported about 40 % below bulk at room temperature (Asheghi et al., J. Heat Transfer, 1998), against 43 % from the model.'
  },
];
