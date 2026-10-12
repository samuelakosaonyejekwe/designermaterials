// Reference benchmarks for suite 5 — Metamaterial Architecture and Unit-Cell Design.
// Expected values are closed-form results of the cited models, evaluated here; none is taken from the suite.

const PI = Math.PI, C0 = 299792458, MU0 = 4e-7 * PI * 1.00000000055, rad = a => a * PI / 180;
const A = 5e-3;                                              // lattice period used for the cell benchmarks, m

// --- Gibson–Ashby hexagonal honeycomb, flexure model (ratios to E_s); cell period a = 2·l·cos θ
const lHex = A / (2 * Math.cos(rad(30))), rHex = 5e-5 / lHex;
const gaE1 = (r, th, hl) => r ** 3 * Math.cos(th) / ((hl + Math.sin(th)) * Math.sin(th) ** 2);
const gaG = (r, th, hl) => r ** 3 * (hl + Math.sin(th)) / (hl * hl * (1 + 2 * hl) * Math.cos(th));
const gaRho = (r, th, hl) => r * (hl + 2) / (2 * Math.cos(th) * (hl + Math.sin(th)));
// --- octet truss: strut length l = a/√2, radius r
const octRho1 = (t, a) => 6 * Math.SQRT2 * PI * (t / 2 / (a / Math.SQRT2)) ** 2;
const octRho = (t, a) => octRho1(t, a) - 54.6 * (t / 2 / (a / Math.SQRT2)) ** 3;
// --- Hashin–Shtrikman upper bound of a porous isotropic solid
function hsE(rho, Es, nu) { const K = Es / (3 * (1 - 2 * nu)), G = Es / (2 * (1 + nu)), Ku = 4 * G * K * rho / (4 * G + 3 * K * (1 - rho)), Gu = G * rho * (9 * K + 8 * G) / (9 * K + 8 * G + 6 * (1 - rho) * (K + 2 * G)); return 9 * Ku * Gu / (3 * Ku + Gu); }
// --- gyroid level set by independent quadrature: threshold c with P(f ≤ c) = 0.30 on a 64³ mid-point grid
const gyroidLevel = (() => { const n = 64, v = new Float64Array(n ** 3); let p = 0; for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const x = (i + 0.5) / n * 2 * PI, y = (j + 0.5) / n * 2 * PI, z = (k + 0.5) / n * 2 * PI; v[p++] = Math.sin(x) * Math.cos(y) + Math.sin(y) * Math.cos(z) + Math.sin(z) * Math.cos(x); } v.sort(); return v[Math.round(0.3 * v.length)]; })();
// --- resonant cells
const fpPendry = (a, r) => C0 / (a * Math.sqrt(2 * PI * Math.log(a / r)));
const fpBelov = (a, r) => C0 / (a * Math.sqrt(2 * PI * (Math.log(a / (2 * PI * r)) + 0.5275)));
const lmS = 3.5e-3 - 4e-4, Lsrr = 2 * MU0 * lmS / PI * (Math.log(lmS / (0.2235 * (4e-4 + 3.5e-5))) - 0.774);
const fHelm = 343 / (2 * PI) * Math.sqrt(PI * (5e-4) ** 2 / (5e-8 * (1e-3 + 1.7 * 5e-4)));
const mim = { m1: 1e-3, m2: 2e-3, k1: 1e5, k2: 5e4 }, w2sq = mim.k2 / mim.m2, Kpi = 4 * mim.k1, bq = (mim.m1 + mim.m2) * w2sq + Kpi;
const mimLo = Math.sqrt((bq - Math.sqrt(bq * bq - 4 * mim.m1 * Kpi * w2sq)) / (2 * mim.m1)) / (2 * PI), mimF2 = Math.sqrt(w2sq) / (2 * PI);
// --- one-dimensional stacks
const nH = 2.4, nL = 1.46, fillQW = nL / (nH + nL), aQW = 550e-9 / 4 * (1 / nH + 1 / nL);
const c1 = Math.sqrt(1.14e11 / 4430), c2 = Math.sqrt(3e9 / 1200), Z1 = 4430 * c1, Z2 = 1200 * c2;
// --- curved-beam (Qiu–Lang–Slocum) preset: E = 2.5 GPa, b = 5 mm, l = 20 mm
const EIl3 = (t) => 2.5e9 * 5e-3 * t ** 3 / 12 / (2e-2) ** 3;
// --- wine-rack cell, θ = 30°, slender struts
const lW = A / (2 * Math.cos(rad(30))), rW = lW / 1e-4, betaX = (rW * Math.sin(rad(60)) - rW ** 3 * Math.cos(rad(60)) * Math.tan(rad(30))) / 1.14e11;

const HEX = { mode: 'lattice', cellType: 'honeycomb2D', theta: 30, hl: 1, cell: A, strut: 5e-5, Es: 1.14e11, nus: 0.34 };
const REE = { mode: 'lattice', cellType: 'reentrant2D', theta: 30, hl: 2, cell: A, strut: 5e-5, Es: 1.14e11, nus: 0.34 };
const OCT = { mode: 'lattice', cellType: 'octet', cell: 1e-5, strut: 1e-7, Es: 3e9, nus: 0.35, process: 'tpp' };
const BEAM = { mode: 'bistable', bisType: 'beam', cell: 2e-2, depth: 5e-3, Es: 2.5e9, ys: 6e7 };

export default [
  {
    id: 's05-regular-honeycomb-nu', suite: 5, title: 'Regular hexagonal honeycomb: in-plane Poisson ratio', mode: 'lattice',
    inputs: HEX, kpi: 'nu', expected: 1.0, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model: ν₁₂ = cos²θ/((h/l + sin θ) sin θ) = 1 for θ = 30°, h/l = 1',
    note: 'Slender walls (t/l = 0.017); axial and shear wall deformation, included in the frame model, lower ν by about 0.1 %.',
  },
  {
    id: 's05-regular-honeycomb-modulus', suite: 5, title: 'Regular hexagonal honeycomb: in-plane Young’s modulus', mode: 'lattice',
    inputs: HEX, kpi: 'Erel', expected: gaE1(rHex, rad(30), 1), tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: E₁/E_s = (t/l)³ cos θ/((h/l + sin θ) sin²θ) = (4/√3)(t/l)³ for regular hexagons',
    note: 'The frame model adds axial and shear compliance of the walls, about 0.2 % at t/l = 0.017.',
  },
  {
    id: 's05-regular-honeycomb-shear', suite: 5, title: 'Regular hexagonal honeycomb: in-plane shear modulus', mode: 'lattice',
    inputs: HEX, kpi: 'Grel', expected: gaG(rHex, rad(30), 1), tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: G₁₂/E_s = (t/l)³ (h/l + sin θ)/((h/l)²(1 + 2h/l) cos θ) = (1/√3)(t/l)³, consistent with isotropy G = E/(2(1 + ν))',
  },
  {
    id: 's05-regular-honeycomb-density', suite: 5, title: 'Regular hexagonal honeycomb: relative density', mode: 'lattice',
    inputs: HEX, kpi: 'rho', expected: gaRho(rHex, rad(30), 1), tol: 1e-6, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: ρ*/ρ_s = (t/l)(h/l + 2)/(2 cos θ (h/l + sin θ)) = (2/√3)(t/l)',
  },
  {
    id: 's05-reentrant-honeycomb-nu', suite: 5, title: 'Re-entrant honeycomb θ = −30°, h/l = 2: negative Poisson ratio', mode: 'lattice',
    inputs: REE, kpi: 'nu', expected: -1.0, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model with re-entrant angle: ν₁₂ = cos²θ/((h/l + sin θ) sin θ) = −1',
    note: 'Auxetic reference case; the finite wall thickness (t/l = 0.017) reduces |ν| by about 0.1 %.',
  },
  {
    id: 's05-reentrant-honeycomb-modulus', suite: 5, title: 'Re-entrant honeycomb θ = −30°, h/l = 2: Young’s modulus along x', mode: 'lattice',
    inputs: REE, kpi: 'Erel', expected: gaE1(rHex, rad(-30), 2), tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: E₁/E_s = (t/l)³ cos θ/((h/l + sin θ) sin²θ)',
  },
  {
    id: 's05-octet-modulus', suite: 5, title: 'Octet truss: relative modulus E/E_s = ρ̄/9', mode: 'lattice',
    inputs: OCT, kpi: 'Erel', expected: octRho1(1e-7, 1e-5) / 9, tol: 0.01, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, effective properties of the octet-truss lattice material: E/E_s = ρ̄/9 along the cube axes (pin-jointed)',
    note: 'Very slender struts (r/l = 0.007) so that the bending contribution of the rigid joints is below 0.1 %.',
  },
  {
    id: 's05-octet-shear-modulus', suite: 5, title: 'Octet truss: relative shear modulus G/E_s = ρ̄/12', mode: 'lattice',
    inputs: OCT, kpi: 'Grel', expected: octRho1(1e-7, 1e-5) / 12, tol: 0.01, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, effective properties of the octet-truss lattice material: G/E_s = ρ̄/12',
  },
  {
    id: 's05-octet-poisson', suite: 5, title: 'Octet truss: Poisson ratio 1/3', mode: 'lattice',
    inputs: OCT, kpi: 'nu', expected: 1 / 3, tol: 0.005, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, effective properties of the octet-truss lattice material: ν = 1/3 along the cube axes',
  },
  {
    id: 's05-octet-connectivity', suite: 5, title: 'Octet truss: nodal connectivity satisfies the rigidity criterion', mode: 'lattice',
    inputs: OCT, kpi: 'Z', expected: 12, abs: 0, kind: 'analytical',
    reference: 'Maxwell rigidity criterion for similarly situated nodes in three dimensions: Z = 12 (Deshpande, Ashby & Fleck, foam topology: bending versus stretching dominated architectures)',
  },
  {
    id: 's05-octet-stretch-dominated', suite: 5, title: 'Octet truss: stiffness carried by axial action', mode: 'lattice',
    inputs: OCT, kpi: 'stretch', expected: 1, tol: 0.01, kind: 'analytical',
    reference: 'Deshpande, Ashby & Fleck: a lattice satisfying Maxwell’s criterion is stretch-dominated, its stiffness is unchanged when the struts are pin-jointed',
  },
  {
    id: 's05-octet-density-node-overlap', suite: 5, title: 'Octet truss: relative density with node overlap', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'octet', cell: A, strut: 5e-4 }, kpi: 'rho', expected: octRho(5e-4, A), tol: 0.015, kind: 'analytical',
    reference: 'Octet-truss relative density with the higher-order nodal correction, ρ̄ = 6√2π(r/l)² − 54.6(r/l)³ (Dong, Deshpande & Wadley, mechanical response of octet-truss lattices)',
    note: 'The suite samples the union of the struts with 90 000 quasi-random points (about 1 % sampling error).',
  },
  {
    id: 's05-simple-cubic-density', suite: 5, title: 'Simple cubic lattice: exact solid fraction of three crossed cylinders', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'simpleCubic', cell: A, strut: 1e-3 }, kpi: 'rho', expected: 3 * PI * 0.1 ** 2 - 8 * Math.SQRT2 * 0.1 ** 3, tol: 0.015, kind: 'analytical',
    reference: 'Steinmetz solids: union of three mutually orthogonal cylinders of radius r through one point = 3πr²a − 8√2·r³ per cubic cell',
    note: 'Quasi-random sampling of the strut union, about 1 % sampling error.',
  },
  {
    id: 's05-square-grid-density', suite: 5, title: 'Square grid: exact area fraction with wall crossing counted once', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'square2D', cell: A, strut: 5e-4 }, kpi: 'rho', expected: 2 * 0.1 - 0.1 ** 2, tol: 0.01, kind: 'analytical',
    reference: 'Plane geometry: two orthogonal walls of thickness t per square cell of side a cover 2(t/a) − (t/a)²',
  },
  {
    id: 's05-triangular-modulus', suite: 5, title: 'Triangular lattice: E/E_s = ρ̄/3', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'triangular2D', cell: A, strut: 5e-5 }, kpi: 'Erel', expected: 2 * Math.sqrt(3) * 0.01 / 3, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids, and Fleck, Deshpande & Ashby: fully triangulated planar lattice, ρ̄ = 2√3·t/l, E = ρ̄E_s/3 (isotropic)',
  },
  {
    id: 's05-triangular-poisson', suite: 5, title: 'Triangular lattice: Poisson ratio 1/3', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'triangular2D', cell: A, strut: 5e-5 }, kpi: 'nu', expected: 1 / 3, tol: 0.005, kind: 'analytical',
    reference: 'Stretch-dominated isotropic planar lattices: ν = 1/3 (Gibson & Ashby, Cellular Solids)',
  },
  {
    id: 's05-kagome-modulus', suite: 5, title: 'Kagome lattice: E/E_s = ρ̄/3', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'kagome2D', cell: A, strut: 5e-5 }, kpi: 'Erel', expected: 2 * Math.sqrt(3) * 0.01 / 3, tol: 0.005, kind: 'analytical',
    reference: 'Kagome lattice with struts of length a/2: ρ̄ = √3·t/l and E = ρ̄E_s/3, the same as the triangular lattice at equal density (Hyun & Torquato; Fleck, Deshpande & Ashby)',
    note: 'Struts of half the cell period: t/l = 0.02, ρ̄ = √3 × 0.02.',
  },
  {
    id: 's05-square-grid-modulus', suite: 5, title: 'Square grid loaded along its walls: E/E_s = t/a', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'square2D', cell: A, strut: 5e-5 }, kpi: 'Erel', expected: 0.01, tol: 1e-6, kind: 'analytical',
    reference: 'Rule of mixtures for walls aligned with the load: one wall of thickness t per period a',
  },
  {
    id: 's05-square-grid-shear', suite: 5, title: 'Square grid in shear: wall bending, G/E_s = (t/l)³/2', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'square2D', cell: A, strut: 5e-5 }, kpi: 'Grel', expected: 0.01 ** 3 / 2, tol: 0.005, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids: square cells sheared along the cell walls deform by wall bending, G = E_s(t/l)³/2',
  },
  {
    id: 's05-rotating-squares-frame', suite: 5, title: 'Rotating-squares cell with thin connections: ν → −1', mode: 'lattice',
    inputs: { mode: 'lattice', cellType: 'rotatingSquares2D', theta: 30, cell: A, strut: 2e-5 }, kpi: 'nu', expected: -1, tol: 0.02, kind: 'analytical',
    reference: 'Grima & Evans, auxetic behaviour from rotating squares: ν = −1 for rigid squares on ideal hinges',
    note: 'The suite solves a welded frame, so the squares are not perfectly rigid; with t/a = 0.004 the frame value is within 1.5 % of the kinematic limit.',
  },
  {
    id: 's05-gyroid-level-set', suite: 5, title: 'Network gyroid at 30 % solid: level-set threshold', mode: 'tpms',
    inputs: { mode: 'tpms', tpmsType: 'gyroid', tpmsMode: 'solid', rd: 0.3, rdTo: 0.3 }, kpi: 'level', expected: gyroidLevel, tol: 0.01, kind: 'numerical-reference',
    reference: 'Volume fraction of the gyroid nodal approximation sin x cos y + sin y cos z + sin z cos x ≤ c, by independent mid-point quadrature on a 64³ grid (c ≈ −0.61; the linear rule ρ̄ ≈ 0.5 + c/3 gives −0.60)',
  },
  {
    id: 's05-schwarz-p-symmetric-level', suite: 5, title: 'Schwarz P network at 50 % solid: threshold zero by symmetry', mode: 'tpms',
    inputs: { mode: 'tpms', tpmsType: 'schwarzP', tpmsMode: 'solid', rd: 0.5, rdTo: 0.5 }, kpi: 'level', expected: 0, abs: 0.01, kind: 'analytical',
    reference: 'cos x + cos y + cos z changes sign under the translation (π, π, π), so the surface f = 0 divides space into two congruent labyrinths (Schwarz primitive surface)',
  },
  {
    id: 's05-gyroid-surface-area', suite: 5, title: 'Gyroid mid-surface area per unit cell', mode: 'tpms',
    inputs: { mode: 'tpms', tpmsType: 'gyroid', tpmsMode: 'solid', rd: 0.5, rdTo: 0.5, cell: 4e-3 }, kpi: 'Sv', expected: 3.0915 / 4e-3, tol: 0.02, kind: 'reference-data',
    reference: 'Area of the Schoen gyroid minimal surface per conventional cubic cell, 3.0915 a² (Schoen; tabulated in the minimal-surface literature)',
    note: 'The suite evaluates S_v = −4·dS₂/dr from the two-point statistics of one cell at 48 voxels per edge (13 lattice directions, extrapolated to zero lag). The trigonometric level set used for the gyroid differs from the exact minimal surface by well under 1 % in area.',
  },
  {
    id: 's05-tpms-power-law', suite: 5, title: 'TPMS stiffness from the empirical power law', mode: 'tpms',
    inputs: { mode: 'tpms', tpmsType: 'gyroid', tpmsMode: 'sheet', rd: 0.3, rdTo: 0.3, plC: 0.5, plN: 2 }, kpi: 'Erel', expected: 0.5 * 0.3 ** 2, tol: 0.005, kind: 'analytical',
    reference: 'Gibson–Ashby scaling form E*/E_s = C·ρ̄ⁿ with C = 0.5, n = 2',
    note: 'The achieved voxel density differs from the 0.30 target by less than 0.1 %.',
  },
  {
    id: 's05-tpms-hashin-shtrikman-cap', suite: 5, title: 'TPMS stiffness capped by the Hashin–Shtrikman upper bound', mode: 'tpms',
    inputs: { mode: 'tpms', tpmsType: 'gyroid', tpmsMode: 'sheet', rd: 0.3, rdTo: 0.3, plC: 1.5, plN: 1, Es: 1.14e11, nus: 0.34 }, kpi: 'Erel', expected: hsE(0.3, 1, 0.34), tol: 0.005, kind: 'analytical',
    reference: 'Hashin & Shtrikman upper bounds on the bulk and shear moduli of a two-phase isotropic material with one void phase',
    note: 'A power law above the bound (1.5·ρ̄) must be replaced by the bound.',
  },
  {
    id: 's05-wire-medium-pendry', suite: 5, title: 'Thin-wire medium: Pendry’s leading-order plasma frequency for his example (a = 5 mm, r = 1 µm)', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'wire', cell: 5e-3, wireR: 1e-6, epsH: 1 }, kpi: 'fp_pendry', expected: 8.2e9, tol: 0.005, kind: 'reference-data',
    reference: 'Pendry, Holden, Stewart & Youngs, extremely low frequency plasmons in metallic mesostructures: ω_p² = 2πc²/(a² ln(a/r)), about 8.2 GHz for 1 µm wires on a 5 mm lattice',
    note: 'Like-for-like comparison: the KPI fp_pendry is the leading-order formula quoted in the reference. The suite’s main value (KPI fp) uses the refined Belov–Tretyakov–Viitanen logarithm ln(a/(2πr)) + 0.5275, which lies 8.7 % higher for this geometry and is checked separately in s05-wire-medium-refined.',
  },
  {
    id: 's05-wire-medium-refined', suite: 5, title: 'Thin-wire medium: refined plasma frequency', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'wire', cell: 5e-3, wireR: 1e-4, epsH: 1 }, kpi: 'fp', expected: fpBelov(5e-3, 1e-4), tol: 1e-4, kind: 'analytical',
    reference: 'Belov, Tretyakov & Viitanen, dispersion and reflection properties of artificial media formed by regular lattices of ideally conducting wires: ω_p² = 2πc²/(a²[ln(a/(2πr)) + 0.5275])',
  },
  {
    id: 's05-srr-lc-resonance', suite: 5, title: 'Split ring loaded with a 100 pF capacitor: L–C resonance', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'srr', cell: 5e-3, srrSize: 3.5e-3, srrW: 4e-4, srrGap: 3e-4, srrT: 3.5e-5, srrCadd: 1e-10, epsH: 1 }, kpi: 'f0', expected: 1 / (2 * PI * Math.sqrt(Lsrr * 1e-10)), tol: 0.002, kind: 'analytical',
    reference: 'Self-inductance of a square loop of side s and equivalent round-wire radius ρ = 0.2235(w + t): L = (2μ₀s/π)[ln(s/ρ) − 0.774] (Grover, Inductance Calculations); f₀ = 1/(2π√(LC))',
    note: 'The 100 pF load dominates the gap and fringing capacitance (about 0.03 pF), which the reference neglects.',
  },
  {
    id: 's05-helmholtz-resonance', suite: 5, title: 'Helmholtz resonator in air: resonance frequency', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'helmholtz', fluid: 'air', hrV: 5e-8, hrRn: 5e-4, hrLn: 1e-3 }, kpi: 'f0', expected: fHelm, tol: 1e-4, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics: f₀ = (c/2π)√(S/(V·L′)), L′ = L + 2 × 0.85r for a flanged neck, c = 343 m/s',
  },
  {
    id: 's05-helmholtz-band-edge', suite: 5, title: 'Helmholtz-resonator array: upper edge of the negative-modulus band', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'helmholtz', fluid: 'air', cell: 5e-3, hrV: 5e-8, hrRn: 5e-4, hrLn: 1e-3 }, kpi: 'band_hi', expected: fHelm * Math.sqrt(1 + 5e-8 / 5e-3 ** 3), tol: 0.02, kind: 'analytical',
    reference: 'Fang, Xi, Xu, Ambati, Srituravanich, Sun & Zhang, ultrasonic metamaterials with negative modulus: 1/B_eff ∝ 1 − Fω₀²/(ω² − ω₀² + iΓω), negative for ω₀ < ω < ω₀√(1 + F) without loss',
    note: 'Viscous loss in the neck (Q ≈ 15) pulls the zero crossing of Re(1/B) about 1 % below the lossless edge.',
  },
  {
    id: 's05-mass-in-mass-resonance', suite: 5, title: 'Mass-in-mass cell: local resonance', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'massinmass', mimM1: mim.m1, mimM2: mim.m2, mimK1: mim.k1, mimK2: mim.k2 }, kpi: 'f0', expected: mimF2, tol: 1e-9, kind: 'analytical',
    reference: 'Huang, Sun & Huang, negative effective mass density in acoustic metamaterials: ω₂ = √(k₂/m₂)',
  },
  {
    id: 's05-mass-in-mass-upper-edge', suite: 5, title: 'Mass-in-mass cell: upper edge of the negative-mass band', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'massinmass', mimM1: mim.m1, mimM2: mim.m2, mimK1: mim.k1, mimK2: mim.k2 }, kpi: 'band_hi', expected: mimF2 * Math.sqrt(1 + mim.m2 / mim.m1), tol: 1e-9, kind: 'analytical',
    reference: 'Huang, Sun & Huang: m_eff = m₁ + m₂ω₂²/(ω₂² − ω²) is negative for ω₂ < ω < ω₂√(1 + m₂/m₁)',
  },
  {
    id: 's05-mass-in-mass-lower-edge', suite: 5, title: 'Mass-in-mass chain: top of the acoustic branch', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'massinmass', mimM1: mim.m1, mimM2: mim.m2, mimK1: mim.k1, mimK2: mim.k2 }, kpi: 'band_lo', expected: mimLo, tol: 1e-6, kind: 'analytical',
    reference: 'Bloch dispersion of the mass-in-mass chain at qa = π: m₁ω⁴ − [(m₁ + m₂)ω₂² + 4k₁]ω² + 4k₁ω₂² = 0 (Huang, Sun & Huang)',
  },
  {
    id: 's05-mass-in-mass-negative-mass', suite: 5, title: 'Mass-in-mass cell: negative effective mass at √(f₂·f_upper)', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'massinmass', mimM1: mim.m1, mimM2: mim.m2, mimK1: mim.k1, mimK2: mim.k2 }, kpi: 'meff', expected: mim.m1 + mim.m2 / (1 - Math.sqrt(1 + mim.m2 / mim.m1)), tol: 0.005, kind: 'analytical',
    reference: 'Huang, Sun & Huang: m_eff = m₁ + m₂/(1 − ω²/ω₂²), evaluated at ω²/ω₂² = √(1 + m₂/m₁)',
    note: 'The suite adds a 2 % loss factor to the internal spring, which changes the real part by about 0.1 %.',
  },
  {
    id: 's05-membrane-first-mode', suite: 5, title: 'Clamped circular membrane: fundamental frequency', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'membrane', memR: 2e-3, memT: 50, memRhoS: 0.1 }, kpi: 'f1', expected: 2.404826 / (2 * PI * 2e-3) * Math.sqrt(50 / 0.1), tol: 1e-4, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics: f₀₁ = (j₀₁/2πR)√(T/ρ_s), j₀₁ = 2.4048',
  },
  {
    id: 's05-membrane-static-stiffness', suite: 5, title: 'Clamped circular membrane: stiffness against uniform pressure', mode: 'resonant',
    inputs: { mode: 'resonant', resType: 'membrane', memT: 50 }, kpi: 'km', expected: 8 * PI * 50, tol: 1e-9, kind: 'analytical',
    reference: 'Membrane statics: w(r) = p(R² − r²)/(4T), mean deflection pR²/(8T), so force/mean deflection = 8πT',
  },
  {
    id: 's05-bragg-mirror-gap', suite: 5, title: 'TiO₂/SiO₂ quarter-wave stack: relative width of the first stop band', mode: 'band',
    inputs: { mode: 'band', waveType: 'optical', cell: aQW, fill: fillQW, nA: nH, nB: nL }, kpi: 'ratio', expected: 4 / PI * Math.asin((nH - nL) / (nH + nL)), tol: 0.001, kind: 'analytical',
    reference: 'Yeh, Optical Waves in Layered Media; Macleod, Thin-Film Optical Filters: Δf/f₀ = (4/π) asin((n_H − n_L)/(n_H + n_L)) for a quarter-wave stack',
  },
  {
    id: 's05-bragg-mirror-centre', suite: 5, title: 'Quarter-wave stack designed for 550 nm: mid-gap frequency', mode: 'band',
    inputs: { mode: 'band', waveType: 'optical', cell: aQW, fill: fillQW, nA: nH, nB: nL }, kpi: 'mid', expected: C0 / 550e-9, tol: 0.001, kind: 'analytical',
    reference: 'Quarter-wave condition n_H d_H = n_L d_L = λ₀/4 (Macleod, Thin-Film Optical Filters)',
  },
  {
    id: 's05-bragg-mirror-attenuation', suite: 5, title: 'Quarter-wave stack: mid-gap decay per period', mode: 'band',
    inputs: { mode: 'band', waveType: 'optical', cell: aQW, fill: fillQW, nA: nH, nB: nL }, kpi: 'att', expected: 8.685889638 * Math.log(nH / nL), tol: 0.001, kind: 'analytical',
    reference: 'Bloch condition at the quarter-wave point: cos qa = −½(n_H/n_L + n_L/n_H), hence Im(qa) = ln(n_H/n_L) (Yeh, Optical Waves in Layered Media)',
  },
  {
    id: 's05-phononic-stack-gap', suite: 5, title: 'Ti-6Al-4V/epoxy quarter-wave elastic stack: relative gap width', mode: 'band',
    inputs: { mode: 'band', waveType: 'elastic', fill: c1 / (c1 + c2), Es: 1.14e11, rhos: 4430, E2: 3e9, rho2: 1200 }, kpi: 'ratio', expected: 4 / PI * Math.asin((Z1 - Z2) / (Z1 + Z2)), tol: 0.001, kind: 'analytical',
    reference: 'Bloch (Rytov) dispersion of a two-layer elastic stack at equal transit times: Δf/f₀ = (4/π) asin(|Z₁ − Z₂|/(Z₁ + Z₂)), Z = ρc',
  },
  {
    id: 's05-quarter-wave-fill', suite: 5, title: 'Fill fraction for equal optical thickness', mode: 'band',
    inputs: { mode: 'band', waveType: 'optical', nA: nH, nB: nL }, kpi: 'fqw', expected: fillQW, tol: 1e-9, kind: 'analytical',
    reference: 'n_H d_H = n_L d_L ⇒ d_H/(d_H + d_L) = n_L/(n_H + n_L)',
  },
  {
    id: 's05-bistable-negative-stiffness', suite: 5, title: 'Curved double beam: negative stiffness on the constrained branch', mode: 'bistable',
    inputs: { ...BEAM, strut: 2e-4, bh: 6e-4 }, kpi: 'kneg', expected: -6 * PI ** 4 * EIl3(2e-4), tol: 0.001, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, a curved-beam bistable mechanism: with the second mode suppressed F = (8π⁴ − 6π⁴Δ)·EIh/l³, i.e. k = −6π⁴EI/l³',
  },
  {
    id: 's05-bistable-snap-force', suite: 5, title: 'Curved double beam (h/t = 10): snap-through force', mode: 'bistable',
    inputs: { ...BEAM, strut: 1e-4, bh: 1e-3 }, kpi: 'ftop', expected: 740 * EIl3(1e-4) * 1e-3, tol: 0.05, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, a curved-beam bistable mechanism: f_top ≈ 740·EIh/l³ for h/t ≥ 6',
    note: 'The quoted coefficient is an approximation of the exact branch intersection, which tends to 8π⁴ = 779 as h/t grows (762 at h/t = 10).',
  },
  {
    id: 's05-bistable-snap-back-force', suite: 5, title: 'Curved double beam (h/t = 10): snap-back force', mode: 'bistable',
    inputs: { ...BEAM, strut: 1e-4, bh: 1e-3 }, kpi: 'fbot', expected: -370 * EIl3(1e-4) * 1e-3, tol: 0.05, kind: 'analytical',
    reference: 'Qiu, Lang & Slocum, a curved-beam bistable mechanism: f_bot ≈ −370·EIh/l³ for h/t ≥ 6',
    note: 'Approximate coefficient; the exact lower intersection tends to −4π⁴ = −390.',
  },
  {
    id: 's05-miura-poisson', suite: 5, title: 'Miura-ori: in-plane Poisson ratio at γ = 60°, θ = 45°', mode: 'bistable',
    inputs: { mode: 'bistable', bisType: 'miura', miGamma: 60, miTheta: 45 }, kpi: 'nuSL', expected: -(Math.cos(rad(45)) ** 2) * Math.tan(rad(60)) ** 2, tol: 1e-9, kind: 'analytical',
    reference: 'Schenk & Guest, geometry of Miura-folded metamaterials: ν_SL = −cos²θ·tan²γ (always negative)',
  },
  {
    id: 's05-miura-height', suite: 5, title: 'Miura-ori: folded height', mode: 'bistable',
    inputs: { mode: 'bistable', bisType: 'miura', miA: 1e-2, miB: 1e-2, miGamma: 60, miTheta: 45 }, kpi: 'H', expected: 1e-2 * Math.sin(rad(45)) * Math.sin(rad(60)), tol: 1e-9, kind: 'analytical',
    reference: 'Schenk & Guest, geometry of Miura-folded metamaterials: H = a·sin θ·sin γ',
  },
  {
    id: 's05-kirigami-stretch-limit', suite: 5, title: 'Staggered-slit kirigami: geometric extension limit', mode: 'bistable',
    inputs: { mode: 'bistable', bisType: 'kirigami', kiLc: 1e-2, kiD: 1e-3, kiG: 1e-3 }, kpi: 'emax', expected: Math.sqrt(1 + 4.5 ** 2) - 1, tol: 1e-9, kind: 'analytical',
    reference: 'Rigid-ligament kinematics: a ligament of length L_b spanning rows d apart straightens to √(L_b² + d²)',
  },
  {
    id: 's05-rotating-squares-area-fraction', suite: 5, title: 'Rotating squares at 30°: solid area fraction', mode: 'bistable',
    inputs: { mode: 'bistable', bisType: 'rotsq', theta: 30 }, kpi: 'rho', expected: 1 / (1 + Math.sin(rad(60))), tol: 1e-9, kind: 'analytical',
    reference: 'Grima & Evans rotating-squares kinematics: cell side 2s(cos φ + sin φ), so 4s²/a² = 1/(1 + sin 2φ)',
  },
  {
    id: 's05-bimaterial-triangle-cte', suite: 5, title: 'Al/Ti bi-material triangle at 45°: negative height expansion', mode: 'negcte',
    inputs: { mode: 'negcte', nteType: 'triangle', alphaB: 2.3e-5, alphaS: 8.6e-6, nteTheta: 45 }, kpi: 'cte', expected: (8.6e-6 - 2.3e-5 * 0.5) / 0.5, tol: 1e-9, kind: 'analytical',
    reference: 'Pin-jointed bi-material triangle (Miller, Mackenzie, Smith & Evans; Steeves et al.): α_h = (α_s − α_b cos²θ)/sin²θ',
  },
  {
    id: 's05-bimaterial-triangle-zero-angle', suite: 5, title: 'Al/Ti bi-material triangle: zero-expansion base angle', mode: 'negcte',
    inputs: { mode: 'negcte', nteType: 'triangle', alphaB: 2.3e-5, alphaS: 8.6e-6 }, kpi: 'th0', expected: Math.acos(Math.sqrt(8.6 / 23)) / PI * 180, tol: 1e-9, kind: 'analytical',
    reference: 'α_h = 0 at cos²θ₀ = α_s/α_b; 52.3° for aluminium (23 × 10⁻⁶ /K) and titanium (8.6 × 10⁻⁶ /K)',
  },
  {
    id: 's05-wine-rack-negative-compressibility', suite: 5, title: 'Wine-rack cell at 30°: negative linear compressibility along x', mode: 'negcte',
    inputs: { mode: 'negcte', nteType: 'winerack', wrTheta: 30, cell: A, strut: 1e-4, Es: 1.14e11 }, kpi: 'bx', expected: betaX, tol: 0.02, kind: 'analytical',
    reference: 'Rigid-jointed rhombic framework under in-plane hydrostatic stress (strut bending in double curvature plus axial strain): β_x = (1/E_s)[(l/t) sin 2θ − (l/t)³ cos 2θ tan θ] (wine-rack mechanism of Baughman et al.; Grima et al.)',
    note: 'Derived here from beam theory; the frame solution also contains shear deformation of the struts (0.4 % at l/t = 29).',
  },
  {
    id: 's05-graded-uniform-limit', suite: 5, title: 'Graded array with a strut ratio of one reduces to the uniform honeycomb', mode: 'graded',
    inputs: { mode: 'graded', cellType: 'honeycomb2D', theta: 30, hl: 1, cell: A, strut: 5e-5, gradRatio: 1, gradVar: 'strut', nGrad: 6, nAcross: 3 }, kpi: 'Eser', expected: gaE1(rHex, rad(30), 1), tol: 0.005, kind: 'analytical',
    reference: 'Series (Reuss) average of identical cells equals the cell modulus; Gibson & Ashby regular honeycomb E/E_s = (4/√3)(t/l)³',
  },
];
