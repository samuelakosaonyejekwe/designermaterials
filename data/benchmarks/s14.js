// Reference benchmarks for Suite 14 — Structural and Nanomechanical Engineering.
// Every expected value is a literal handbook figure or is evaluated here from a closed-form expression; nothing is taken from the suite.
const PI = Math.PI, S3 = Math.sqrt(3), rad = d => d * PI / 180;

// ---- 2-D lattices (slender walls, t/l = 0.01, so that the bending / stretching closed forms apply)
const ES = 2e9, TL = 0.01, slender = { mode: 'lattice2d', cell_l: 5, wall_t: 0.05, Es: ES };
const gaNu12 = (hl, th) => Math.cos(th) ** 2 / ((hl + Math.sin(th)) * Math.sin(th));
const gaE1 = (hl, th) => TL ** 3 * Math.cos(th) / ((hl + Math.sin(th)) * Math.sin(th) ** 2);
const gaE2 = (hl, th) => TL ** 3 * (hl + Math.sin(th)) / Math.cos(th) ** 3;

// ---- curved-beam snap-through (Qiu, Lang and Slocum), second mode constrained
const NS = { E: 2e9, l: 40e-3, t: 0.4e-3, b: 5e-3, n: 2 }, nsEI = NS.E * NS.b * NS.t ** 3 / 12;
const nsS = Q => Math.sqrt(1 - 16 / (3 * Q * Q)); // branch intersections F₁ = F₃ at Δ = 1 ± s
const nsIn = { mode: 'negstiff', ns_l: 40, ns_t: 0.4, ns_b: 5, ns_n: 2, ns_mode: 'double', Es: 2e9 };

// ---- wine-rack framework: hinging + stretching of rigid-jointed struts, prismatic-wall out-of-plane coupling
const WR = { E: 2e9, nu: 0.35, L: 10, w: 0.1, th: rad(25) }, wrRho = 2 * WR.w / (WR.L * Math.sin(2 * WR.th));
const wrStretch = WR.L * Math.sin(2 * WR.th) / (WR.E * WR.w), wrHinge = (WR.L / WR.w) ** 3 * Math.cos(2 * WR.th) / WR.E;
const wrBx = wrStretch - wrHinge * Math.tan(WR.th) - WR.nu / (wrRho * WR.E), wrBy = wrStretch + wrHinge / Math.tan(WR.th) - WR.nu / (wrRho * WR.E), wrBz = (1 - 2 * WR.nu) / (wrRho * WR.E);
const wrIn = { mode: 'nlc', nlc_model: 'winerack', wr_theta: 25, wr_L: 10, wr_w: 0.1, Es: 2e9, nu_s: 0.35, beam: 'euler' };

// ---- nanowire (bulk limit: surface terms switched off)
const NW = { E: 140e9, D: 50e-9, L: 2e-6, rho: 5606 }, nwI = PI * NW.D ** 4 / 64, nwA = PI * NW.D ** 2 / 4;
const nwIn = { mode: 'nanowire', E_n: 140e9, nw_D: 50, nw_L: 2, rho_s: 5606, E_surf: 0, tau0: 0 };

// ---- contact: diamond indenter E = 1141 GPa, ν = 0.07 (the values used by Oliver and Pharr)
const erOf = (E, nu) => 1 / ((1 - nu * nu) / E + (1 - 0.07 ** 2) / 1141e9);
const HZ = { P: 5e-6, R: 1e-6, Er: erOf(140e9, 0.3) }, hzA = Math.cbrt(3 * HZ.P * HZ.R / (4 * HZ.Er));
// synthetic Berkovich unloading data generated from fused silica (E = 72 GPa, ν = 0.17, H = 9 GPa) at 5 mN
const FS = { E: 72e9, nu: 0.17, H: 9e9, P: 5e-3 }, fsA = FS.P / FS.H, fsHc = Math.sqrt(fsA / 24.5), fsS = 2 * 1.034 * erOf(FS.E, FS.nu) * Math.sqrt(fsA / PI), fsHmax = fsHc + 0.75 * FS.P / fsS;

// ---- poroelastic consolidation (default layer)
const PO = { E: 1.5e9, nu: 0.3, al: 0.6, phi: 0.3, k: 1e-16, mu: 1e-3, Kf: 2.2e9, H: 10e-3, s: 1e5 };
const poK = PO.E / (3 * (1 - 2 * PO.nu)), poG = PO.E / (2 * (1 + PO.nu)), poKv = poK + 4 * poG / 3, poInvM = PO.phi / PO.Kf + (PO.al - PO.phi) * (1 - PO.al) / poK;
const poC = (PO.k / PO.mu) / (poInvM + PO.al ** 2 / poKv), poIn = { mode: 'poro', E_d: 1.5e9, nu_d: 0.3, biot: 0.6, poros: 0.3, perm: 1e-16, mu_f: 1, K_f: 2.2e9, H_layer: 10, sig_load: 1e5 };

export default [
  // ------------------------------------------------------------------ 2-D lattices
  {
    id: 's14-regular-honeycomb-nu', suite: 14, title: 'Regular hexagonal honeycomb: in-plane Poisson ratio', mode: 'lattice2d',
    inputs: { ...slender, topology: 'honeycomb' }, kpi: 'nu12', expected: 1.0, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model for hexagonal honeycombs',
    note: 'Slender walls (t/l = 0.01); the reference neglects axial and shear deformation of the walls, which lower ν by a term of order (t/l)².',
  },
  {
    id: 's14-regular-honeycomb-modulus', suite: 14, title: 'Regular hexagonal honeycomb: E/E_s = 2.3094 (t/l)³', mode: 'lattice2d',
    inputs: { ...slender, topology: 'honeycomb' }, kpi: 'E_rel', expected: 4 / S3 * TL ** 3, tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — E*/E_s = (4/√3)(t/l)³ for regular hexagons',
    note: 'Frame finite elements include axial and shear wall compliance, an O((t/l)²) correction to the bending-only reference.',
  },
  {
    id: 's14-regular-honeycomb-shear', suite: 14, title: 'Regular hexagonal honeycomb: G/E_s = 0.5774 (t/l)³', mode: 'lattice2d',
    inputs: { ...slender, topology: 'honeycomb' }, kpi: 'G12', expected: ES * TL ** 3 / S3, tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — in-plane shear modulus of regular hexagons, G* = E*/(2(1 + ν*)) with ν* = 1',
  },
  {
    id: 's14-regular-honeycomb-density', suite: 14, title: 'Regular hexagonal honeycomb: relative density (2/√3) t/l', mode: 'lattice2d',
    inputs: { ...slender, topology: 'honeycomb' }, kpi: 'rho', expected: 2 / S3 * TL, tol: 1e-6, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — relative density of hexagonal honeycombs (thin-wall limit)',
  },
  {
    id: 's14-reentrant-nu-minus-one', suite: 14, title: 'Re-entrant honeycomb θ = −30°, h/l = 2: ν₁₂ = −1 (auxetic)', mode: 'lattice2d',
    inputs: { ...slender, topology: 'reentrant', h_l: 2, theta_c: -30 }, kpi: 'nu12', expected: -1.0, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — ν₁₂ = cos²θ/[(h/l + sin θ) sin θ] evaluated for a re-entrant cell',
    note: 'Slender walls (t/l = 0.01); bending-only reference.',
  },
  {
    id: 's14-reentrant-nu-theta-20', suite: 14, title: 'Re-entrant honeycomb θ = −20°, h/l = 2: ν₁₂ from the bending model', mode: 'lattice2d',
    inputs: { ...slender, topology: 'reentrant', h_l: 2, theta_c: -20 }, kpi: 'nu12', expected: gaNu12(2, rad(-20)), tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — ν₁₂ = cos²θ/[(h/l + sin θ) sin θ]',
  },
  {
    id: 's14-reentrant-modulus-e1', suite: 14, title: 'Re-entrant honeycomb θ = −30°, h/l = 2: E₁/E_s = 2.3094 (t/l)³', mode: 'lattice2d',
    inputs: { ...slender, topology: 'reentrant', h_l: 2, theta_c: -30 }, kpi: 'E1', expected: ES * gaE1(2, rad(-30)), tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — E₁/E_s = (t/l)³ cos θ/[(h/l + sin θ) sin²θ]',
  },
  {
    id: 's14-reentrant-modulus-e2', suite: 14, title: 'Re-entrant honeycomb θ = −20°, h/l = 2: E₂ from the bending model', mode: 'lattice2d',
    inputs: { ...slender, topology: 'reentrant', h_l: 2, theta_c: -20 }, kpi: 'E2', expected: ES * gaE2(2, rad(-20)), tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — E₂/E_s = (t/l)³ (h/l + sin θ)/cos³θ',
  },
  {
    id: 's14-triangular-modulus', suite: 14, title: 'Triangular lattice: E/E_s = 1.1547 t/l (stretch-dominated)', mode: 'lattice2d',
    inputs: { ...slender, topology: 'triangular' }, kpi: 'E_rel', expected: 2 / S3 * TL, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — fully triangulated (equilateral) lattice, E*/E_s = (2/√3) t/l',
    note: 'Pin-jointed reference; the rigid joints of the frame model add a bending contribution of order (t/l)².',
  },
  {
    id: 's14-triangular-nu', suite: 14, title: 'Triangular lattice: ν = 1/3', mode: 'lattice2d',
    inputs: { ...slender, topology: 'triangular' }, kpi: 'nu12', expected: 1 / 3, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — in-plane Poisson ratio of the equilateral triangular lattice',
  },
  {
    id: 's14-kagome-modulus', suite: 14, title: 'Kagome lattice: E/E_s = 0.5774 t/l', mode: 'lattice2d',
    inputs: { ...slender, topology: 'kagome' }, kpi: 'E_rel', expected: TL / S3, tol: 0.005, kind: 'analytical',
    reference: 'Hyun & Torquato; Wang & McDowell — in-plane stiffness of the kagome lattice, E*/E_s = ρ̄/3 with ρ̄ = √3 t/l',
  },
  {
    id: 's14-square-modulus', suite: 14, title: 'Square lattice loaded along the walls: E/E_s = t/l', mode: 'lattice2d',
    inputs: { ...slender, topology: 'square' }, kpi: 'E_rel', expected: TL, tol: 1e-6, kind: 'analytical',
    reference: 'Rule of mixtures for walls loaded axially (Gibson & Ashby, Cellular Solids — square cells)',
  },
  {
    id: 's14-square-nu-zero', suite: 14, title: 'Square lattice loaded along the walls: ν = 0', mode: 'lattice2d',
    inputs: { ...slender, topology: 'square' }, kpi: 'nu12', expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — square cells loaded along the walls have no lateral coupling',
  },
  {
    id: 's14-hexachiral-nu', suite: 14, title: 'Hexachiral lattice: ν → −1 (rigid-ring limit)', mode: 'lattice2d',
    inputs: { mode: 'lattice2d', topology: 'chiral', cell_l: 10, wall_t: 0.05, r_R: 0.1, Es: ES, beam: 'euler' }, kpi: 'nu12', expected: -1.0, tol: 0.06, kind: 'analytical',
    reference: 'Prall & Lakes, Properties of a chiral honeycomb with a Poisson\'s ratio of −1',
    note: 'The reference assumes rigid rings. The suite models each ring as a 12-sided polygon of beams with the ligament thickness, so the rings ovalise and ν is a few per cent less negative; the tolerance reflects that documented model difference.',
  },
  // ------------------------------------------------------------------ negative stiffness
  {
    id: 's14-snap-through-force', suite: 14, title: 'Cosine curved beam, Q = 6: snap-through force π⁴(2 + 6s)·EIh/l³', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 2.4 }, kpi: 'F_top', expected: NS.n * PI ** 4 * (2 + 6 * nsS(6)) * nsEI * 2.4e-3 / NS.l ** 3, tol: 0.002, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — intersection of the first-mode branch with F₃ = 8π⁴ − 6π⁴Δ',
    note: 's = √(1 − 16/(3Q²)); two beams in parallel.',
  },
  {
    id: 's14-snap-through-force-740', suite: 14, title: 'Cosine curved beam, Q = 6: peak force ≈ 740 EIh/l³', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 2.4 }, kpi: 'F_top', expected: NS.n * 740 * nsEI * 2.4e-3 / NS.l ** 3, tol: 0.02, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — design approximation f_top ≈ 740 EIh/l³',
    note: 'The quoted coefficient is a rounded design value; the exact coefficient depends weakly on Q (734 at Q = 6, 779 as Q → ∞).',
  },
  {
    id: 's14-snap-back-force', suite: 14, title: 'Cosine curved beam, Q = 6: minimum force π⁴(2 − 6s)·EIh/l³ (bistable, negative)', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 2.4 }, kpi: 'F_bot', expected: NS.n * PI ** 4 * (2 - 6 * nsS(6)) * nsEI * 2.4e-3 / NS.l ** 3, tol: 0.002, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — second intersection of F₁ with F₃',
  },
  {
    id: 's14-snap-through-displacement', suite: 14, title: 'Cosine curved beam, Q = 6: snap-through displacement h(1 − s)', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 2.4 }, kpi: 'd_top', expected: 2.4e-3 * (1 - nsS(6)), tol: 0.002, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — first intersection of F₁ with F₃',
  },
  {
    id: 's14-negative-stiffness-slope', suite: 14, title: 'Cosine curved beam: negative stiffness −6π⁴EI/l³ per beam', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 2.4 }, kpi: 'k_neg', expected: -NS.n * 6 * PI ** 4 * nsEI / NS.l ** 3, tol: 0.001, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — slope of the third-mode branch F₃ = 8π⁴ − 6π⁴Δ',
  },
  {
    id: 's14-bistability-threshold', suite: 14, title: 'Bistability threshold Q = 4/√3 ≈ 2.31: minimum force is zero', mode: 'negstiff',
    inputs: { ...nsIn, ns_h: 0.4 * 4 / S3 }, kpi: 'F_bot', expected: 0, abs: 1e-3, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, A curved-beam bistable mechanism — bistability requires Q = h/t > 2.31',
    note: 'Absolute tolerance 1 mN against a snap-through force of about 0.6 N for this beam.',
  },
  // ------------------------------------------------------------------ negative linear compressibility
  {
    id: 's14-winerack-negative-beta-x', suite: 14, title: 'Wine-rack framework θ = 25°: negative linear compressibility along x', mode: 'nlc',
    inputs: wrIn, kpi: 'beta_x', expected: wrBx * 1e12, tol: 0.01, kind: 'analytical',
    reference: 'Wine-rack hinging mechanism with strut stretching (Baughman et al.; Grima et al.) — β_x = (L/E_s w) sin 2θ − (L/w)³ cos 2θ tan θ/E_s − ν_s/(ρ̄E_s)',
    note: 'Slender struts (w/L = 0.01), Euler–Bernoulli elements; the closed form treats hinging and stretching as uncoupled.',
  },
  {
    id: 's14-winerack-beta-y', suite: 14, title: 'Wine-rack framework θ = 25°: conjugate (positive) compressibility along y', mode: 'nlc',
    inputs: wrIn, kpi: 'beta_y', expected: wrBy * 1e12, tol: 0.01, kind: 'analytical',
    reference: 'Wine-rack hinging mechanism with strut stretching — β_y = (L/E_s w) sin 2θ + (L/w)³ cos 2θ cot θ/E_s − ν_s/(ρ̄E_s)',
  },
  {
    id: 's14-winerack-beta-v-positive', suite: 14, title: 'Wine-rack framework: volumetric compressibility stays positive', mode: 'nlc',
    inputs: wrIn, kpi: 'beta_V', expected: (wrBx + wrBy + wrBz) * 1e12, tol: 0.01, kind: 'analytical',
    reference: 'Thermodynamic stability: β_V = β_x + β_y + β_z > 0 (Nye, Physical Properties of Crystals)',
    note: 'Expected value is the sum of the three closed-form linear compressibilities; it is positive although β_x is negative.',
  },
  {
    id: 's14-crystal-silicon-bulk-modulus', suite: 14, title: 'Cubic crystal (silicon constants): bulk modulus (C₁₁ + 2C₁₂)/3', mode: 'nlc',
    inputs: { mode: 'nlc', nlc_model: 'crystal', c11: 165.7, c22: 165.7, c33: 165.7, c12: 63.9, c13: 63.9, c23: 63.9, c44: 79.6 }, kpi: 'K', expected: 97.8e9, tol: 0.005, kind: 'reference-data',
    reference: 'Silicon elastic constants C₁₁ = 165.7, C₁₂ = 63.9 GPa (McSkimin); bulk modulus 97.8 GPa',
  },
  {
    id: 's14-crystal-negative-beta3', suite: 14, title: 'Transversely isotropic stiffness set: β₃ = (C₁₁ + C₁₂ − 2C₁₃)/[C₃₃(C₁₁ + C₁₂) − 2C₁₃²] < 0', mode: 'nlc',
    inputs: { mode: 'nlc', nlc_model: 'crystal', c11: 20, c22: 20, c33: 80, c12: 8, c13: 28, c23: 28, c44: 10 }, kpi: 'beta3', expected: (20 + 8 - 2 * 28) / (80 * (20 + 8) - 2 * 28 ** 2) * 1e3, tol: 1e-6, kind: 'analytical',
    reference: 'Nye, Physical Properties of Crystals — linear compressibility of hexagonal/tetragonal crystals from the stiffness constants',
    note: 'Illustrative constants with strong axial–transverse coupling; β₃ = −41.67 TPa⁻¹.',
  },
  // ------------------------------------------------------------------ 3-D lattice scaling
  {
    id: 's14-octet-modulus', suite: 14, title: 'Octet truss: E*/E_s = ρ̄/9', mode: 'lattice3d',
    inputs: { mode: 'lattice3d', topo3: 'octet', rel_dens: 0.1, Es: 2e9 }, kpi: 'E', expected: 0.1 * 2e9 / 9, tol: 1e-6, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, Effective properties of the octet-truss lattice material',
  },
  {
    id: 's14-octet-buckling-strength', suite: 14, title: 'Octet truss: elastic-buckling strength (ρ̄/3)·π²E_s(r/l)²/4, pin-jointed', mode: 'lattice3d',
    inputs: { mode: 'lattice3d', topo3: 'octet', rel_dens: 0.1, Es: 2e9, end_fix: 1 }, kpi: 'sb', expected: (0.1 / 3) * PI ** 2 * 2e9 * (0.1 / (6 * Math.SQRT2 * PI)) / 4, tol: 1e-6, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby — Euler buckling of solid circular struts with ρ̄ = 6√2π(r/l)²',
  },
  {
    id: 's14-octet-strut-radius', suite: 14, title: 'Octet truss: strut radius from ρ̄ = 6√2π(r/l)², l = a/√2', mode: 'lattice3d',
    inputs: { mode: 'lattice3d', topo3: 'octet', rel_dens: 0.05, cell3: 5 }, kpi: 'r', expected: Math.sqrt(0.05 / (6 * Math.SQRT2 * PI)) * 5e-6 / Math.SQRT2, tol: 1e-6, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby — first-order relative density of the octet truss',
  },
  {
    id: 's14-bending-foam-modulus', suite: 14, title: 'Bending-dominated open cell: E*/E_s = ρ̄²', mode: 'lattice3d',
    inputs: { mode: 'lattice3d', topo3: 'kelvin', rel_dens: 0.05, Es: 70e9 }, kpi: 'E', expected: 0.05 ** 2 * 70e9, tol: 1e-6, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — open-cell foam scaling E*/E_s ≈ ρ̄² (proportionality constant ≈ 1)',
  },
  // ------------------------------------------------------------------ nanowire
  {
    id: 's14-euler-buckling-pinned', suite: 14, title: 'Euler buckling of a pinned–pinned circular wire: P_cr = π²EI/L²', mode: 'nanowire',
    inputs: { ...nwIn, nw_bc: 'pinned' }, kpi: 'Pcr', expected: PI ** 2 * NW.E * nwI / NW.L ** 2, tol: 1e-6, kind: 'analytical',
    reference: 'Timoshenko & Gere, Theory of Elastic Stability — Euler column',
  },
  {
    id: 's14-euler-buckling-clamped', suite: 14, title: 'Euler buckling of a clamped–clamped wire: P_cr = 4π²EI/L²', mode: 'nanowire',
    inputs: { ...nwIn, nw_bc: 'fixed' }, kpi: 'Pcr', expected: 4 * PI ** 2 * NW.E * nwI / NW.L ** 2, tol: 1e-6, kind: 'analytical',
    reference: 'Timoshenko & Gere, Theory of Elastic Stability — effective length factor 0.5 for built-in ends',
  },
  {
    id: 's14-cantilever-stiffness', suite: 14, title: 'Cantilever tip stiffness 3EI/L³', mode: 'nanowire',
    inputs: { ...nwIn, nw_bc: 'cantilever' }, kpi: 'k', expected: 3 * NW.E * nwI / NW.L ** 3, tol: 1e-6, kind: 'analytical',
    reference: 'Euler–Bernoulli beam theory (Timoshenko, Strength of Materials)',
  },
  {
    id: 's14-clamped-beam-frequency', suite: 14, title: 'Clamped–clamped wire: f₁ = (4.730²/2π)√(EI/ρAL⁴)', mode: 'nanowire',
    inputs: { ...nwIn, nw_bc: 'fixed' }, kpi: 'f1', expected: 4.730041 ** 2 / (2 * PI) * Math.sqrt(NW.E * nwI / (NW.rho * nwA * NW.L ** 4)), tol: 1e-4, kind: 'analytical',
    reference: 'Blevins, Formulas for Natural Frequency and Mode Shape — clamped–clamped beam, λ₁ = 4.7300',
  },
  {
    id: 's14-surface-elasticity-modulus', suite: 14, title: 'Circular wire with surface elasticity: E_eff = E + 8E_s/D', mode: 'nanowire',
    inputs: { ...nwIn, nw_D: 20, E_surf: 10 }, kpi: 'E_eff', expected: 140e9 + 8 * 10 / 20e-9, tol: 1e-6, kind: 'analytical',
    reference: 'Miller & Shenoy, Size-dependent elastic properties of nanosized structural elements',
  },
  // ------------------------------------------------------------------ contact
  {
    id: 's14-hertz-contact-radius', suite: 14, title: 'Hertz sphere-on-flat: contact radius (3PR/4E_r)^(1/3)', mode: 'contact',
    inputs: { mode: 'contact', E_n: 140e9, nu_n: 0.3, tip_R: 1000, P_ind: 5 }, kpi: 'a', expected: hzA, tol: 1e-6, kind: 'analytical',
    reference: 'Johnson, Contact Mechanics — Hertz theory; diamond indenter E = 1141 GPa, ν = 0.07',
  },
  {
    id: 's14-hertz-peak-pressure', suite: 14, title: 'Hertz sphere-on-flat: peak pressure (6PE_r²/π³R²)^(1/3)', mode: 'contact',
    inputs: { mode: 'contact', E_n: 140e9, nu_n: 0.3, tip_R: 1000, P_ind: 5 }, kpi: 'p0', expected: Math.cbrt(6 * HZ.P * HZ.Er ** 2 / (PI ** 3 * HZ.R ** 2)), tol: 1e-6, kind: 'analytical',
    reference: 'Johnson, Contact Mechanics — Hertz theory',
  },
  {
    id: 's14-hertz-approach', suite: 14, title: 'Hertz sphere-on-flat: approach δ = (9P²/16RE_r²)^(1/3)', mode: 'contact',
    inputs: { mode: 'contact', E_n: 140e9, nu_n: 0.3, tip_R: 1000, P_ind: 5 }, kpi: 'delta', expected: Math.cbrt(9 * HZ.P ** 2 / (16 * HZ.R * HZ.Er ** 2)), tol: 1e-6, kind: 'analytical',
    reference: 'Johnson, Contact Mechanics — Hertz theory',
  },
  {
    id: 's14-oliver-pharr-silica-modulus', suite: 14, title: 'Oliver–Pharr analysis recovers E = 72 GPa from fused-silica unloading data', mode: 'contact',
    inputs: { mode: 'contact', nu_n: 0.17, Pmax: 5, hmax: fsHmax * 1e9, S_unl: fsS * 1e-3 }, kpi: 'E_op', expected: 72e9, tol: 0.002, kind: 'reference-data',
    reference: 'Oliver & Pharr, An improved technique for determining hardness and elastic modulus — fused silica calibration (E ≈ 72 GPa, H ≈ 9 GPa)',
    note: 'Peak depth and unloading stiffness are generated in this file from E = 72 GPa, ν = 0.17, H = 9 GPa for an ideal Berkovich tip (A = 24.5h_c², ε = 0.75, β = 1.034); the suite must invert them.',
  },
  {
    id: 's14-oliver-pharr-silica-hardness', suite: 14, title: 'Oliver–Pharr analysis recovers H = 9 GPa from fused-silica unloading data', mode: 'contact',
    inputs: { mode: 'contact', nu_n: 0.17, Pmax: 5, hmax: fsHmax * 1e9, S_unl: fsS * 1e-3 }, kpi: 'H', expected: 9e9, tol: 0.002, kind: 'reference-data',
    reference: 'Oliver & Pharr, An improved technique for determining hardness and elastic modulus — fused silica calibration',
  },
  // ------------------------------------------------------------------ strength
  {
    id: 's14-hall-petch-copper', suite: 14, title: 'Hall–Petch: copper, d = 1 µm, σ₀ = 25 MPa, k = 0.11 MPa·√m → 135 MPa', mode: 'strength',
    inputs: { mode: 'strength', d_grain: 1000, sig0: 25, k_hp: 110 }, kpi: 'sy_hp', expected: 135e6, tol: 1e-6, kind: 'analytical',
    reference: 'Hall–Petch relation σ_y = σ₀ + k d^(−1/2) with the textbook copper constants (Meyers & Chawla, Mechanical Behavior of Materials)',
  },
  {
    id: 's14-griffith-fracture-stress', suite: 14, title: 'LEFM fracture stress K_Ic/(Y√(πa))', mode: 'strength',
    inputs: { mode: 'strength', KIc: 1e6, flaw: 1000, Y_geo: 1, E_n: 140e9 }, kpi: 'sf', expected: 1e6 / Math.sqrt(PI * 1e-6), tol: 1e-6, kind: 'analytical',
    reference: 'Griffith / Irwin fracture criterion (Anderson, Fracture Mechanics)',
    note: 'A 1 µm flaw in a 1 MPa·√m solid fails at 564 MPa.',
  },
  {
    id: 's14-critical-flaw-size', suite: 14, title: 'Critical flaw size (K_Ic/Yσ)²/π', mode: 'strength',
    inputs: { mode: 'strength', KIc: 1e6, Y_geo: 1.12, sig_app: 2e8 }, kpi: 'ac', expected: (1e6 / (1.12 * 2e8)) ** 2 / PI, tol: 1e-6, kind: 'analytical',
    reference: 'Linear elastic fracture mechanics (Anderson, Fracture Mechanics)',
  },
  {
    id: 's14-fracture-energy', suite: 14, title: 'Irwin relation G_c = K_Ic²(1 − ν²)/E', mode: 'strength',
    inputs: { mode: 'strength', KIc: 1e6, E_n: 140e9, nu_n: 0.3 }, kpi: 'Gc', expected: 1e12 * (1 - 0.09) / 140e9, tol: 1e-6, kind: 'analytical',
    reference: 'Irwin relation for plane strain (Anderson, Fracture Mechanics)',
  },
  {
    id: 's14-thermal-mismatch-stress', suite: 14, title: 'Biaxial thermal-mismatch film stress E_f/(1 − ν_f)·(α_s − α_f)·ΔT', mode: 'strength',
    inputs: { mode: 'strength', E_n: 140e9, nu_n: 0.3, cte_f: 4.7e-6, cte_s: 2.6e-6, dT: -400 }, kpi: 's_film', expected: 140e9 / 0.7 * (2.6e-6 - 4.7e-6) * -400, tol: 1e-6, kind: 'analytical',
    reference: 'Freund & Suresh, Thin Film Materials — equi-biaxial mismatch stress',
  },
  {
    id: 's14-stoney-curvature', suite: 14, title: 'Stoney curvature 6σ_f t_f(1 − ν_s)/(E_s t_s²)', mode: 'strength',
    inputs: { mode: 'strength', E_n: 140e9, nu_n: 0.3, cte_f: 4.7e-6, cte_s: 2.6e-6, dT: -400, t_film: 500, t_sub: 500, E_sub: 130e9, nu_sub: 0.28 }, kpi: 'kappa',
    expected: 6 * 168e6 * 500e-9 * (1 - 0.28) / (130e9 * (500e-6) ** 2), tol: 1e-6, kind: 'analytical',
    reference: 'Stoney equation (Freund & Suresh, Thin Film Materials)',
  },
  // ------------------------------------------------------------------ poroelastic consolidation
  {
    id: 's14-terzaghi-t90', suite: 14, title: 'Terzaghi consolidation: t₉₀ = 0.848 H²/c_v', mode: 'poro',
    inputs: poIn, kpi: 't90', expected: 0.848 * PO.H ** 2 / poC, tol: 1e-4, kind: 'analytical',
    reference: 'Terzaghi one-dimensional consolidation, time factor T_v = 0.848 at U = 90 % (Terzaghi, Peck & Mesri, Soil Mechanics in Engineering Practice); Biot diffusivity from Wang, Theory of Linear Poroelasticity',
  },
  {
    id: 's14-terzaghi-u50', suite: 14, title: 'Terzaghi consolidation: U = 50 % at T_v = 0.197', mode: 'poro',
    inputs: { ...poIn, t_cons: 0.197 * PO.H ** 2 / poC }, kpi: 'U', expected: 0.5, tol: 0.01, kind: 'analytical',
    reference: 'Terzaghi one-dimensional consolidation, tabulated U(T_v) (Terzaghi, Peck & Mesri)',
    note: 'The suite integrates the diffusion equation by backward-Euler finite differences (61 nodes, 240 steps); the tabulated time factor is rounded to three figures.',
  },
  {
    id: 's14-terzaghi-u90', suite: 14, title: 'Terzaghi consolidation: U = 90 % at T_v = 0.848', mode: 'poro',
    inputs: { ...poIn, t_cons: 0.848 * PO.H ** 2 / poC }, kpi: 'U', expected: 0.9, tol: 0.005, kind: 'analytical',
    reference: 'Terzaghi one-dimensional consolidation, tabulated U(T_v) (Terzaghi, Peck & Mesri)',
  },
  {
    id: 's14-skempton-incompressible-grains', suite: 14, title: 'Skempton coefficient with incompressible grains: B = 1/(1 + φK/K_f)', mode: 'poro',
    inputs: { ...poIn, biot: 1 }, kpi: 'B', expected: 1 / (1 + PO.phi * poK / PO.Kf), tol: 1e-6, kind: 'analytical',
    reference: 'Skempton pore-pressure coefficient (Wang, Theory of Linear Poroelasticity)',
  },
  {
    id: 's14-drained-settlement', suite: 14, title: 'Final drained settlement σ₀H/(K + 4G/3)', mode: 'poro',
    inputs: poIn, kpi: 'winf', expected: PO.s * PO.H / poKv, tol: 1e-6, kind: 'analytical',
    reference: 'Oedometric (uniaxial-strain) compression of a linear elastic layer (Wang, Theory of Linear Poroelasticity)',
  },
];
