// Reference benchmarks for suite 3 — Atomic, Molecular and Crystal Structure.
// Expected values are handbook data or closed-form results evaluated here from CODATA constants.

const N_A = 6.02214076e23, K_B_EV = 8.617333262e-5, H = 6.62607015e-34, M_E = 9.1093837015e-31, Q_E = 1.602176634e-19, C0 = 299792458;
const KE2 = 14.399645;                 // e²/(4πε₀) in eV·Å
const EV_KJ = 96.48533;                // 1 eV per particle in kJ/mol
const EVA3_GPA = 160.21766;            // 1 eV/Å³ in GPa
const deg = Math.PI / 180, CU_KA1 = 1.5406;
const dFromTwoTheta = tt => CU_KA1 / (2 * Math.sin(tt / 2 * deg));
const eLambda = kV => { const E = Q_E * kV * 1e3; return H / Math.sqrt(2 * M_E * E * (1 + E / (2 * M_E * C0 * C0))) * 1e10; };   // Å

// Lennard-Jones fcc lattice sums (Ashcroft & Mermin): A₁₂ = 12.13188, A₆ = 14.45392; hcp: 12.13229, 14.45489
const A12 = 12.13188, A6 = 14.45392, EPS_AR = 0.0104, SIG_AR = 3.405;
const ljD0 = Math.pow(2 * A12 / A6, 1 / 6), ljE = A6 * A6 / (2 * A12), ljB = 4 * A12 * Math.pow(A6 / A12, 2.5);
const ljEhcp = 14.45489 ** 2 / (2 * 12.13229);

// Matthews–Blakeslee critical thickness, 60° dislocations (cos α = cos λ = ½), solved by fixed-point iteration (upper root)
const aSi = 5.431, aGe = 5.658, fGe = (aSi - aGe) / aGe, bGe = aGe / Math.SQRT2;
function matthewsBlakeslee(f, b, nu) { const C = b * (1 - nu / 4) / (8 * Math.PI * Math.abs(f) * (1 + nu) * 0.5); let h = 10 * b; for (let i = 0; i < 500; i++) h = C * (Math.log(h / b) + 1); return h; }   // Å
function peopleBean(f, b, a, nu) { const K = (1 - nu) / (1 + nu) * b * b / (16 * Math.PI * Math.SQRT2 * a * f * f); let h = 100 * b; for (let i = 0; i < 500; i++) h = K * Math.log(h / b); return h; }   // Å

// Wurtzite ZnO, powder-diffraction reference cell: a = 3.2498 Å, c = 5.2066 Å, Z = 2
const aZnO = 3.2498, cZnO = 5.2066, vZnO = Math.sqrt(3) / 2 * aZnO * aZnO * cZnO, M_ZNO = 65.38 + 15.999;

export default [
  {
    id: 's03-zno-xray-density', suite: 3, title: 'X-ray density of wurtzite ZnO', mode: 'crystal',
    inputs: { mode: 'crystal', crystal: 'ZnO_wurtzite' }, kpi: 'rho', expected: 2 * M_ZNO * 1e-3 / (N_A * vZnO * 1e-30), tol: 0.002, kind: 'reference-data',
    reference: 'ρ = Z·M/(N_A·V) with the powder-diffraction reference cell of zincite, a = 3.2498 Å, c = 5.2066 Å, Z = 2 (5675 kg/m³; handbook bulk density 5606 kg/m³)',
  },
  {
    id: 's03-si-xray-density', suite: 3, title: 'X-ray density of silicon', mode: 'crystal',
    inputs: { mode: 'crystal', crystal: 'Si' }, kpi: 'rho', expected: 2329, tol: 0.001, kind: 'reference-data',
    reference: 'Density of single-crystal silicon at room temperature, 2329 kg/m³ (semiconductor data handbooks)',
  },
  {
    id: 's03-au-xray-density', suite: 3, title: 'X-ray density of gold', mode: 'crystal',
    inputs: { mode: 'crystal', crystal: 'Au' }, kpi: 'rho', expected: 19300, tol: 0.003, kind: 'reference-data',
    reference: 'Density of gold at room temperature, 19 300 kg/m³ (metals handbooks)',
  },
  {
    id: 's03-zno-strongest-line', suite: 3, title: 'Strongest Cu Kα reflection of ZnO: (101)', mode: 'crystal',
    inputs: { mode: 'crystal', crystal: 'ZnO_wurtzite', rad: 'CuKa' }, kpi: 'tt_max', expected: 36.25, tol: 0.001, kind: 'reference-data',
    reference: 'Powder-diffraction reference pattern of zincite: (100) 31.77°, (002) 34.42°, (101) 36.25° 2θ for Cu Kα₁, with (101) the 100 % line',
  },
  {
    id: 's03-zno-100-bragg-angle', suite: 3, title: 'Bragg angle of ZnO (100), Cu Kα₁', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'ZnO_wurtzite', rad: 'CuKa' }, kpi: 'tt1', expected: 31.77, tol: 0.001, kind: 'reference-data',
    reference: 'Powder-diffraction reference pattern of zincite: (100) at 31.77° 2θ for Cu Kα₁',
  },
  {
    id: 's03-zno-002-d-spacing', suite: 3, title: 'd-spacing of ZnO (002) from its tabulated Bragg angle 34.42°', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'ZnO_wurtzite', rad: 'CuKa', h1: 0, k1: 0, l1: 2, h2: 1, k2: 0, l2: 1 }, kpi: 'd1', expected: dFromTwoTheta(34.42), tol: 0.001, kind: 'reference-data',
    reference: 'Bragg law λ = 2d sin θ with the reference angle of zincite (002), 34.42° 2θ, λ(Cu Kα₁) = 1.5406 Å',
  },
  {
    id: 's03-zno-101-d-spacing', suite: 3, title: 'd-spacing of ZnO (101) from its tabulated Bragg angle 36.25°', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'ZnO_wurtzite', rad: 'CuKa', h1: 0, k1: 0, l1: 2, h2: 1, k2: 0, l2: 1 }, kpi: 'd2', expected: dFromTwoTheta(36.25), tol: 0.001, kind: 'reference-data',
    reference: 'Bragg law λ = 2d sin θ with the reference angle of zincite (101), 36.25° 2θ, λ(Cu Kα₁) = 1.5406 Å',
  },
  {
    id: 's03-zno-interplanar-angle', suite: 3, title: 'Angle between (002) and (101) of hexagonal ZnO', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'ZnO_wurtzite', h1: 0, k1: 0, l1: 2, h2: 1, k2: 0, l2: 1 }, kpi: 'angle', expected: Math.atan(2 * cZnO / (Math.sqrt(3) * aZnO)) / deg, tol: 0.001, kind: 'analytical',
    reference: 'Hexagonal reciprocal metric: tan φ = |a*|/|c*| = 2c/(√3·a) for (001) against (101) (crystallography textbooks)',
  },
  {
    id: 's03-si-111-bragg-angle', suite: 3, title: 'First allowed reflection of silicon: (111), Cu Kα₁', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'Si', rad: 'CuKa' }, kpi: 'tt1', expected: 28.44, tol: 0.001, kind: 'reference-data',
    reference: 'Silicon powder line-position standard: (111) at 28.44° 2θ for Cu Kα₁; (100) and (110) are absent for the diamond structure',
  },
  {
    id: 's03-cubic-interplanar-angle', suite: 3, title: 'Angle between (100) and (111) in a cubic crystal', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'Si', h1: 1, k1: 0, l1: 0, h2: 1, k2: 1, l2: 1 }, kpi: 'angle', expected: Math.acos(1 / Math.sqrt(3)) / deg, tol: 1e-8, kind: 'analytical',
    reference: 'cos φ = (h₁h₂ + k₁k₂ + l₁l₂)/(|h₁||h₂|) for cubic crystals: 54.74° (crystallography textbooks)',
  },
  {
    id: 's03-electron-wavelength-200kv', suite: 3, title: 'Relativistic electron wavelength at 200 kV', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'Si', rad: 'electron', kV: 200 }, kpi: 'lambda', expected: 0.02508, tol: 5e-4, kind: 'reference-data',
    reference: 'Electron wavelength table for transmission electron microscopy: 2.508 pm at 200 kV (Williams & Carter, Transmission Electron Microscopy)',
  },
  {
    id: 's03-electron-wavelength-100kv', suite: 3, title: 'Relativistic electron wavelength at 100 kV', mode: 'bragg',
    inputs: { mode: 'bragg', crystal: 'Si', rad: 'electron', kV: 100 }, kpi: 'lambda', expected: eLambda(100), tol: 1e-6, kind: 'analytical',
    reference: 'λ = h/√(2m₀eV(1 + eV/2m₀c²)) with CODATA constants; 3.701 pm at 100 kV',
  },
  {
    id: 's03-epitaxy-ge-si-mismatch', suite: 3, title: 'Lattice mismatch of Ge on Si', mode: 'epitaxy',
    inputs: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', epiOri: 'direct' }, kpi: 'mismatch', expected: 4.18, tol: 0.002, kind: 'reference-data',
    reference: 'Lattice parameters Si 5.431 Å and Ge 5.658 Å give the widely quoted 4.18 % mismatch (semiconductor heteroepitaxy textbooks)',
  },
  {
    id: 's03-epitaxy-matthews-blakeslee', suite: 3, title: 'Matthews–Blakeslee critical thickness of Ge on Si(001)', mode: 'epitaxy',
    inputs: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', epiOri: 'direct', nu: 0.28 }, kpi: 'hc', expected: matthewsBlakeslee(fGe, bGe, 0.28) * 1e-10, tol: 0.01, kind: 'analytical',
    reference: 'Matthews & Blakeslee force balance for a single strained layer, h_c = b(1 − ν cos²α)/(8π|f|(1 + ν)cos λ)·[ln(h_c/b) + 1], 60° dislocations',
    note: 'Independent fixed-point solution with a_Ge = 5.658 Å; about 1.2 nm, consistent with the few-monolayer coherent limit of Ge on Si.',
  },
  {
    id: 's03-epitaxy-people-bean', suite: 3, title: 'People–Bean critical thickness of Ge on Si(001)', mode: 'epitaxy',
    inputs: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', epiOri: 'direct', nu: 0.28 }, kpi: 'hc_pb', expected: peopleBean(fGe, bGe, aGe, 0.28) * 1e-10, tol: 0.01, kind: 'analytical',
    reference: 'People & Bean energy-balance criterion, h_c = [(1 − ν)/(1 + ν)]·b²/(16π√2·a·f²)·ln(h_c/b)',
  },
  {
    id: 's03-epitaxy-coherent-stress', suite: 3, title: 'Biaxial stress of a coherent Ge film on Si', mode: 'epitaxy',
    inputs: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', epiOri: 'direct', nu: 0.28, Gsh: 5.5e10, hfilm: 5e-10 }, kpi: 'sigma', expected: 2 * 5.5e10 * 1.28 / 0.72 * fGe, tol: 0.002, kind: 'analytical',
    reference: 'Thin-film biaxial stress σ = 2G(1 + ν)/(1 − ν)·ε = E/(1 − ν)·ε with the full misfit strain (thin-film mechanics textbooks)',
    note: 'A 0.5 nm film is below the critical thickness, so it carries the whole misfit strain of −4.0 %.',
  },
  {
    id: 's03-epitaxy-thermal-strain', suite: 3, title: 'Thermal-mismatch strain of Ge on Si cooled from 550 °C', mode: 'epitaxy',
    inputs: { mode: 'epitaxy', film: 'Ge', substrate: 'Si', T_C: 550, alF: 5.9e-6, alS: 2.6e-6 }, kpi: 'eps_th', expected: 100 * (5.9e-6 - 2.6e-6) * 525, tol: 1e-9, kind: 'analytical',
    reference: 'ε_th = (α_film − α_substrate)(T_growth − T_room) with α_Ge = 5.9 × 10⁻⁶ /K and α_Si = 2.6 × 10⁻⁶ /K',
  },
  {
    id: 's03-lj-fcc-cohesive-energy', suite: 3, title: 'Static Lennard-Jones fcc crystal: cohesive energy 8.61 ε (argon parameters)', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'Ecoh', expected: ljE * EPS_AR, tol: 0.001, kind: 'analytical',
    reference: 'Lennard-Jones lattice sums for fcc, A₁₂ = 12.13188 and A₆ = 14.45392: E = −εA₆²/(2A₁₂) = −8.610 ε (Ashcroft & Mermin, Solid State Physics; Kittel)',
  },
  {
    id: 's03-lj-fcc-lattice-parameter', suite: 3, title: 'Static Lennard-Jones fcc crystal: nearest-neighbour distance 1.09 σ', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'a0', expected: Math.SQRT2 * ljD0 * SIG_AR, tol: 0.0005, kind: 'analytical',
    reference: 'Lennard-Jones lattice sums for fcc: r_nn = (2A₁₂/A₆)^(1/6) σ = 1.0902 σ, a = √2·r_nn (Ashcroft & Mermin, Solid State Physics)',
  },
  {
    id: 's03-lj-fcc-bulk-modulus', suite: 3, title: 'Static Lennard-Jones fcc crystal: bulk modulus 75.2 ε/σ³', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'B0', expected: ljB * EPS_AR / SIG_AR ** 3 * EVA3_GPA, tol: 0.01, kind: 'analytical',
    reference: 'Lennard-Jones lattice sums for fcc: B = 4A₁₂(A₆/A₁₂)^(5/2) ε/σ³ = 75.2 ε/σ³ (Ashcroft & Mermin, Solid State Physics)',
    note: 'The suite obtains B₀ from a third-order Birch–Murnaghan fit over ±2 % in lattice parameter; the fit form differs slightly from the exact 12–6 equation of state.',
  },
  {
    id: 's03-lj-pressure-derivative', suite: 3, title: 'Lennard-Jones solid: pressure derivative of the bulk modulus', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'Bp', expected: 8, tol: 0.01, kind: 'analytical',
    reference: 'Mie (m–n) lattice equation of state: B′ = (m + n + 6)/3 = 8 for the 12–6 potential',
  },
  {
    id: 's03-lj-fcc-hcp-energy-difference', suite: 3, title: 'Lennard-Jones fcc–hcp energy difference', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'dE_struct', expected: (ljEhcp - ljE) * EPS_AR * 1e3, abs: 6e-4, kind: 'analytical',
    reference: 'Lennard-Jones lattice sums: hcp (A₁₂ = 12.13229, A₆ = 14.45489) lies 0.01 % below fcc (Ashcroft & Mermin; Kittel)',
    note: 'A difference of 9 µeV per atom between two sums of about 90 meV: the absolute tolerance corresponds to the six-figure precision of the tabulated lattice sums and the 9.5-neighbour-distance cut-off with continuum tail used by the suite.',
  },
  {
    id: 's03-madelung-nacl', suite: 3, title: 'Madelung constant of the rock-salt structure', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'NaCl' }, kpi: 'madelung', expected: 1.74756, tol: 1e-5, kind: 'analytical',
    reference: 'Madelung constant of NaCl referred to the nearest-neighbour distance, 1.74756 (Kittel, Introduction to Solid State Physics)',
  },
  {
    id: 's03-madelung-cscl', suite: 3, title: 'Madelung constant of the caesium-chloride structure', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'CsCl' }, kpi: 'madelung', expected: 1.76267, tol: 1e-5, kind: 'analytical',
    reference: 'Madelung constant of CsCl referred to the nearest-neighbour distance, 1.76267 (Kittel, Introduction to Solid State Physics)',
  },
  {
    id: 's03-madelung-zincblende', suite: 3, title: 'Madelung constant of the zinc-blende structure', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'ZnS_zincblende' }, kpi: 'madelung', expected: 1.6381, tol: 5e-5, kind: 'analytical',
    reference: 'Madelung constant of zinc blende referred to the nearest-neighbour distance, 1.6381 (Kittel, Introduction to Solid State Physics)',
  },
  {
    id: 's03-born-lande-nacl', suite: 3, title: 'Born–Landé lattice energy of NaCl (n = 8)', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'NaCl', bornN: 8 }, kpi: 'U_BL', expected: -1.74756 * KE2 * EV_KJ / (5.640 / 2) * (1 - 1 / 8), tol: 0.001, kind: 'analytical',
    reference: 'Born–Landé equation U = −N_A·A·e²/(4πε₀r₀)·(1 − 1/n) with r₀ = a/2 = 2.820 Å; about −753 kJ/mol (inorganic chemistry textbooks)',
  },
  {
    id: 's03-born-haber-nacl', suite: 3, title: 'Born–Haber lattice enthalpy of NaCl', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'NaCl' }, kpi: 'BH', expected: -787, tol: 0.003, kind: 'reference-data',
    reference: 'Thermochemical (Born–Haber) lattice enthalpy of NaCl, −787 kJ/mol (inorganic and physical chemistry textbooks)',
  },
  {
    id: 's03-born-bulk-modulus-nacl', suite: 3, title: 'Bulk modulus of NaCl from the Born model', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'ionic', salt: 'NaCl', bornN: 8 }, kpi: 'B', expected: 24.4, tol: 0.05, kind: 'reference-data',
    reference: 'Measured room-temperature bulk modulus of NaCl, about 24 GPa (Kittel, Introduction to Solid State Physics, alkali-halide table)',
    note: 'The Born exponent n = 8 is the ion-core rule value, not fitted to compressibility, so a few per cent deviation is expected.',
  },
  {
    id: 's03-schottky-nacl', suite: 3, title: 'Schottky defect fraction in NaCl at 800 °C (2.3 eV per pair)', mode: 'defect',
    inputs: { mode: 'defect', crystal: 'NaCl', T_C: 800, Esch: 2.3 }, kpi: 'x_sch', expected: Math.exp(-2.3 / (2 * K_B_EV * 1073.15)), tol: 0.001, kind: 'analytical',
    reference: 'Schottky equilibrium in an AB crystal: n/N = exp(−E_S/2k_BT) (Kittel, Introduction to Solid State Physics)',
  },
  {
    id: 's03-vacancy-fraction', suite: 3, title: 'Equilibrium vacancy fraction for E_v = 1.28 eV at 1000 °C without formation entropy', mode: 'defect',
    inputs: { mode: 'defect', crystal: 'Cu', T_C: 1000, Evac: 1.28, Svib: 0 }, kpi: 'x_vac', expected: Math.exp(-1.28 / (K_B_EV * 1273.15)), tol: 0.001, kind: 'analytical',
    reference: 'Boltzmann vacancy concentration n/N = exp(−E_v/k_BT) (Kittel; physical metallurgy textbooks)',
  },
  {
    id: 's03-goldschmidt-srtio3', suite: 3, title: 'Goldschmidt tolerance factor of undoped SrTiO₃', mode: 'defect',
    inputs: { mode: 'defect', crystal: 'SrTiO3', xdop: 0 }, kpi: 'tolerance', expected: 1.00, tol: 0.005, kind: 'reference-data',
    reference: 'Goldschmidt tolerance factor with Shannon radii: SrTiO₃ is the textbook ideal cubic perovskite, t ≈ 1.00',
  },
  {
    id: 's03-dopant-density-azo', suite: 3, title: 'Al atoms per cm³ for 2 % of the Zn sites in ZnO', mode: 'defect',
    inputs: { mode: 'defect', crystal: 'ZnO_wurtzite', dopant: 'Al3+', xdop: 2 }, kpi: 'n_dop', expected: 0.02 * 2 / (vZnO * 1e-24), tol: 0.002, kind: 'analytical',
    reference: 'Two Zn sites per wurtzite cell of volume (√3/2)a²c = 47.62 Å³: Zn-site density 4.20 × 10²² cm⁻³',
  },
  {
    id: 's03-sphere-surface-to-volume', suite: 3, title: 'Surface-to-volume ratio of a 3 nm sphere', mode: 'nano',
    inputs: { mode: 'nano', crystal: 'Au', nshape: 'sphere', dnm: 3e-9 }, kpi: 'SV', expected: 6 / 3e-9, tol: 1e-6, kind: 'analytical',
    reference: 'Sphere geometry: S/V = 6/d',
  },
  {
    id: 's03-graphene-specific-surface-area', suite: 3, title: 'Specific surface area of a graphene monolayer (both faces)', mode: 'nano',
    inputs: { mode: 'nano', crystal: 'graphene', nshape: 'sheet', layers: 1 }, kpi: 'ssa', expected: 2.63e6, tol: 0.005, kind: 'reference-data',
    reference: 'Theoretical specific surface area of graphene, 2630 m²/g (two faces of a honeycomb sheet with a = 2.46 Å)',
  },
  {
    id: 's03-nanotube-diameter', suite: 3, title: 'Diameter of a (10,10) carbon nanotube', mode: 'nano',
    inputs: { mode: 'nano', crystal: 'graphene', nshape: 'tube', tn: 10, tm: 10, Lnm: 5e-9 }, pub: 'geom.diameter', expected: 2.46e-10 * Math.sqrt(300) / Math.PI, tol: 0.002, kind: 'analytical',
    reference: 'd = a√(n² + nm + m²)/π with a = 2.46 Å: 1.356 nm for (10,10) (Saito, Dresselhaus & Dresselhaus, Physical Properties of Carbon Nanotubes)',
  },
  {
    id: 's03-nanotube-specific-surface-area', suite: 3, title: 'Specific surface area of a single-wall carbon nanotube (inner and outer wall)', mode: 'nano',
    inputs: { mode: 'nano', crystal: 'graphene', nshape: 'tube', tn: 10, tm: 10, Lnm: 5e-9 }, kpi: 'ssa', expected: 2.63e6, tol: 0.01, kind: 'reference-data',
    reference: 'A single-wall nanotube is a rolled graphene sheet: 2630 m²/g counting both wall faces (1315 m²/g for the outer face alone)',
  },
  // ------------------------------------------------------------------ published experiments
  {
    id: 's03-exp-argon-lattice-constant-0k', suite: 3, title: 'Solid argon: measured X-ray lattice constant extrapolated to 0 K against the Lennard-Jones crystal with zero-point motion', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'a0_zp', expected: 5.30017 * 1.00202, tol: 0.005, kind: 'experiment',
    measured: { value: 5.30017, uncertainty: 0.00008, unit: 'kxu' },
    source: { authors: 'Peterson, Batchelder, Simmons', title: 'Measurements of X-Ray Lattice Constant, Thermal Expansivity, and Isothermal Compressibility of Argon Crystals', journal: 'Phys. Rev.', year: 1966, doi: '10.1103/PhysRev.150.703', url: 'https://doi.org/10.1103/PhysRev.150.703' },
    conditions: 'Argon crystals of 99.998 % purity grown by directional solidification; rotating-camera X-ray method between 2.3 and 83.2 K, lattice constant extrapolated to 0 K',
    match: ['pot', 'potModel', 'lattice'],
    reference: 'Peterson, Batchelder & Simmons, Phys. Rev. 150, 703 (1966)',
    note: 'Value stated in the abstract in kX units; converted with 1 kxu = 1.00202 Å (5.3109 Å). The suite minimises the 12-6 lattice energy (gas-phase argon parameters ε = 0.0104 eV, σ = 3.405 Å) plus the quasi-harmonic zero-point energy (9/8)k_BΘ_D, with Θ_D from the lattice-sum Einstein frequency and no adjustable constant: 5.315 Å, +0.1 %. The static lattice alone gives 5.250 Å (−1.2 %), so the earlier shortfall was the missing quantum expansion. The tolerance covers the Debye-spectrum approximation and the pair-potential parameters.'
  },
  {
    id: 's03-exp-argon-bulk-modulus-4k', suite: 3, title: 'Solid argon: measured isothermal compressibility at 4.25 K against the Lennard-Jones crystal with zero-point motion', mode: 'potential',
    inputs: { mode: 'potential', potModel: 'pair', pot: 'Ar', lattice: 'fcc' }, kpi: 'B0_zp', expected: 1 / 3.75e-10 / 1e9, tol: 0.05, kind: 'experiment',
    measured: { value: 3.75e-11, uncertainty: null, unit: 'cm²/dyn (isothermal compressibility)' },
    source: { authors: 'Peterson, Batchelder, Simmons', title: 'Measurements of X-Ray Lattice Constant, Thermal Expansivity, and Isothermal Compressibility of Argon Crystals', journal: 'Phys. Rev.', year: 1966, doi: '10.1103/PhysRev.150.703', url: 'https://doi.org/10.1103/PhysRev.150.703' },
    conditions: 'Argon crystals of 99.998 % purity; isothermal compressibility from the X-ray lattice constant of crystals under pressure, value at 4.25 K',
    match: ['pot', 'potModel', 'lattice'],
    reference: 'Peterson, Batchelder & Simmons, Phys. Rev. 150, 703 (1966)',
    note: 'Compressibility 3.75×10⁻¹¹ cm²/dyn from the abstract, i.e. B = 2.67 GPa. With the zero-point term the model gives 2.66 GPa (−0.2 %); the static pair model gives 3.18 GPa (+19 %). The difference is therefore the quantum correction (lattice expansion of 1.2 % acting on a potential with B′ ≈ 8, plus the volume dependence of the zero-point energy), not a limitation of the pair potential itself. No uncertainty is quoted in the abstract; the 5 % tolerance covers the Debye-spectrum approximation of the zero-point energy and the neglect of three-body (Axilrod–Teller) forces. The model Debye temperature (79 K) is 14 % below the calorimetric 92 K, which bounds the accuracy of the correction.'
  },
];
