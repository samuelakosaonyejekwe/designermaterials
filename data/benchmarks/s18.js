// Suite 18 — reference benchmarks for multiscale modelling and homogenisation.
// Every expected value is a literal handbook number or is evaluated here from a textbook closed form; nothing is taken from the suite.

const KG = (E, nu) => ({ K: E / (3 * (1 - 2 * nu)), G: E / (2 * (1 + nu)) });
const Eof = (K, G) => 9 * K * G / (3 * K + G);
// Hashin–Shtrikman expressions in Hashin's original form with phase 1 as the reference ("coating") phase
const hsK = (f2, p1, p2) => p1.K + f2 / (1 / (p2.K - p1.K) + 3 * (1 - f2) / (3 * p1.K + 4 * p1.G));
const hsG = (f2, p1, p2) => p1.G + f2 / (1 / (p2.G - p1.G) + 6 * (1 - f2) * (p1.K + 2 * p1.G) / (5 * p1.G * (3 * p1.K + 4 * p1.G)));

// default two-phase system of the suite: epoxy-like matrix and ZnO-like inclusion
const epoxy = KG(3e9, 0.35), oxide = KG(1.4e11, 0.35), f = 0.2;
// stiff matrix with soft spherical inclusions (Mori–Tanaka then coincides with the Hashin–Shtrikman upper bound)
const stiffM = KG(7e10, 0.3), softI = KG(2e9, 0.4), fSoft = 0.3;
const mtG = (() => { const beta = 6 * (stiffM.K + 2 * stiffM.G) / (5 * (3 * stiffM.K + 4 * stiffM.G)); return stiffM.G + fSoft * (softI.G - stiffM.G) / (1 + (1 - fSoft) * beta * (softI.G - stiffM.G) / stiffM.G); })();
const mtK = (() => { const alpha = 3 * stiffM.K / (3 * stiffM.K + 4 * stiffM.G); return stiffM.K + fSoft * (softI.K - stiffM.K) / (1 + (1 - fSoft) * alpha * (softI.K - stiffM.K) / stiffM.K); })();
// Kerner thermal expansion of a particulate composite (default phases)
const kerner = (1 - f) * 6e-5 + f * 4.3e-6 + f * (1 - f) * (4.3e-6 - 6e-5) * (oxide.K - epoxy.K) / ((1 - f) * epoxy.K + f * oxide.K + 3 * oxide.K * epoxy.K / (4 * epoxy.G));

// RVE with equal shear moduli (Hill): μ = 3e9/(2·1.35) in both phases, phase B has ν = 0.2
const muEq = 3e9 / 2.7, lamA = 2 * muEq * 0.35 / (1 - 0.7), lamB = 2 * muEq * 0.2 / (1 - 0.4), E_B_eq = 2 * muEq * 1.2;
const hillC11 = 1 / (0.5 / (lamA + 2 * muEq) + 0.5 / (lamB + 2 * muEq));

// laminate cases
const lamEin = (() => { const fA = 0.6, fB = 0.4, EA = 7e10, EB = 2e11, nA = 0.35, nB = 0.2; return fA * EA + fB * EB + fA * fB * EA * EB * (nA - nB) ** 2 / (fA * EA * (1 - nB * nB) + fB * EB * (1 - nA * nA)); })();
const lamKout = (() => { const per = 10e-9, nb = 40, H = nb * per; return H / (H * (0.5 / 170 + 0.5 / 30) + (2 * nb - 1) * 3e-9); })();
const lamCte = (() => { const bA = 7e10 / (1 - 0.35), bB = 2e11 / (1 - 0.2); return (0.6 * bA * 2.3e-5 + 0.4 * bB * 1.2e-5) / (0.6 * bA + 0.4 * bB); })();

// silicon (cubic) and zinc oxide (hexagonal) Voigt–Reuss–Hill closed forms, constants in GPa
const si = (() => { const C11 = 165.7, C12 = 63.9, C44 = 79.6, K = (C11 + 2 * C12) / 3, GV = (C11 - C12 + 3 * C44) / 5, GR = 5 * (C11 - C12) * C44 / (4 * C44 + 3 * (C11 - C12)), G = 0.5 * (GV + GR); return { K, G, E: Eof(K, G), nu: (3 * K - 2 * G) / (2 * (3 * K + G)), vl: Math.sqrt((K + 4 * G / 3) * 1e9 / 2329) }; })();
const zno = (() => { const C11 = 209.7, C12 = 121.1, C13 = 105.1, C33 = 210.9, C44 = 42.5, C66 = (C11 - C12) / 2, M = C11 + C12 + 2 * C33 - 4 * C13, C2 = (C11 + C12) * C33 - 2 * C13 * C13; const KV = (2 * (C11 + C12) + C33 + 4 * C13) / 9, KR = C2 / M, GV = (M + 12 * C44 + 12 * C66) / 30, GR = 2.5 * C2 * C44 * C66 / (3 * KV * C44 * C66 + C2 * (C44 + C66)); return { K: 0.5 * (KV + KR), G: 0.5 * (GV + GR) }; })();
// Lennard-Jones argon, static fcc crystal: lattice sums A12, A6 of the fcc lattice
const lj = (() => { const kB = 1.380649e-23, qe = 1.602176634e-19, NA = 6.02214076e23, A12 = 12.13188, A6 = 14.45392, eps = 119.8 * kB, sig = 3.405e-10, d = sig * Math.pow(2 * A12 / A6, 1 / 6); return { epsEV: eps / qe, K: 4 * eps * A12 * Math.pow(A6 / A12, 2.5) / sig ** 3, rho: 39.948e-3 / NA / (d ** 3 / Math.SQRT2) }; })();

export default [
  // ---------- two-phase bounds and mean-field estimates ----------
  {
    id: 's18-hs-lower-E', suite: 18, title: 'Hashin–Shtrikman lower bound on Young\'s modulus (soft matrix, stiff spheres)', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', vf: f, E_m: 3e9, nu_m: 0.35, E_i: 1.4e11, nu_i: 0.35 }, kpi: 'E_hs_lo',
    expected: Eof(hsK(f, epoxy, oxide), hsG(f, epoxy, oxide)), tol: 1e-6, kind: 'analytical',
    reference: 'Hashin & Shtrikman, variational bounds for the elastic moduli of multiphase materials (J. Mech. Phys. Solids)',
    note: 'Bounds on K and G evaluated in Hashin\'s original form with the soft phase as reference; E = 9KG/(3K + G).'
  },
  {
    id: 's18-hs-upper-E', suite: 18, title: 'Hashin–Shtrikman upper bound on Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', vf: f, E_m: 3e9, nu_m: 0.35, E_i: 1.4e11, nu_i: 0.35 }, kpi: 'E_hs_hi',
    expected: Eof(hsK(1 - f, oxide, epoxy), hsG(1 - f, oxide, epoxy)), tol: 1e-6, kind: 'analytical',
    reference: 'Hashin & Shtrikman, variational bounds for the elastic moduli of multiphase materials (J. Mech. Phys. Solids)',
    note: 'Stiff phase taken as the reference phase.'
  },
  {
    id: 's18-mori-tanaka-G', suite: 18, title: 'Mori–Tanaka shear modulus, soft spheres in a stiff matrix (Eshelby tensor form)', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', est: 'mt', vf: fSoft, E_m: 7e10, nu_m: 0.3, E_i: 2e9, nu_i: 0.4 }, kpi: 'G_eff',
    expected: mtG, tol: 1e-6, kind: 'analytical',
    reference: 'Benveniste, re-examination of the Mori–Tanaka method (Mech. Mater.); Eshelby spherical-inclusion factor β = 6(K + 2G)/(5(3K + 4G))'
  },
  {
    id: 's18-mori-tanaka-K', suite: 18, title: 'Mori–Tanaka bulk modulus, soft spheres in a stiff matrix (Eshelby tensor form)', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', est: 'mt', vf: fSoft, E_m: 7e10, nu_m: 0.3, E_i: 2e9, nu_i: 0.4 }, kpi: 'K_eff',
    expected: mtK, tol: 1e-6, kind: 'analytical',
    reference: 'Benveniste, re-examination of the Mori–Tanaka method (Mech. Mater.); Eshelby spherical-inclusion factor α = 3K/(3K + 4G)'
  },
  {
    id: 's18-self-consistent-voids', suite: 18, title: 'Self-consistent scheme with spherical voids at ν = 0.2: E/E₀ = 1 − 2f', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', vf: 0.3, E_m: 1e12, nu_m: 0.2, E_i: 1e6, nu_i: 0.2 }, kpi: 'E_sc',
    expected: 1e12 * (1 - 2 * 0.3), tol: 1e-3, kind: 'analytical',
    reference: 'Budiansky / Hill self-consistent scheme; ν = 0.2 is its fixed point for porous solids (see Christensen, Mechanics of Composite Materials)',
    note: 'The inclusion is 10⁶ times softer than the matrix rather than a true void, which changes the result by parts per million.'
  },
  {
    id: 's18-dem-porous', suite: 18, title: 'Differential scheme, spherical pores at ν = 0.2: E/E₀ = (1 − p)²', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', E_m: 5e10, nu_m: 0.2, poros: 0.3 }, kpi: 'E_porous',
    expected: 5e10 * 0.7 ** 2, tol: 1e-3, kind: 'analytical',
    reference: 'Differential effective-medium scheme for porous solids (Norris; Zimmerman, Compressibility of Sandstones): ν = 0.2 fixed point',
    note: 'The suite integrates the scheme numerically (RK4, 120 steps).'
  },
  {
    id: 's18-maxwell-garnett-k', suite: 18, title: 'Maxwell (Maxwell-Garnett) conductivity of dispersed spheres', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', vf: 0.2, k_m: 0.2, k_i: 50 }, kpi: 'k_eff',
    expected: 0.2 * (50 + 2 * 0.2 + 2 * 0.2 * (50 - 0.2)) / (50 + 2 * 0.2 - 0.2 * (50 - 0.2)), tol: 1e-6, kind: 'analytical',
    reference: 'Maxwell, A Treatise on Electricity and Magnetism — conductivity of a dilute suspension of spheres'
  },
  {
    id: 's18-maxwell-insulating', suite: 18, title: 'Conducting matrix with insulating spheres: σ/σ_m = 2(1 − f)/(2 + f)', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', vf: 0.25, sig_m: 1e6, sig_i: 1e-16 }, kpi: 'sigma_eff',
    expected: 1e6 * 2 * 0.75 / 2.25, tol: 1e-6, kind: 'analytical',
    reference: 'Maxwell, A Treatise on Electricity and Magnetism — insulating-sphere limit (Maxwell–Eucken relation)'
  },
  {
    id: 's18-halpin-tsai-long-fibre', suite: 18, title: 'Halpin–Tsai, aligned fibres of very high aspect ratio: rule of mixtures', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'fibre', orient: 'aligned', aspect: 1e4, vf: 0.3, E_m: 3e9, E_i: 7e10 }, kpi: 'E_eff',
    expected: 0.7 * 3e9 + 0.3 * 7e10, tol: 0.005, kind: 'analytical',
    reference: 'Halpin & Kardos, The Halpin–Tsai equations: a review (Polym. Eng. Sci.) — continuous-fibre limit ζ → ∞',
    note: 'At aspect ratio 10⁴ (ζ = 2×10⁴) the Halpin–Tsai longitudinal modulus is within about 0.1 % of the Voigt mixture rule.'
  },
  {
    id: 's18-kerner-cte', suite: 18, title: 'Thermal expansion of a particulate composite (Kerner = Levin with the Hashin–Shtrikman bulk modulus)', mode: 'bounds',
    inputs: { mode: 'bounds', shape: 'sphere', est: 'mt', vf: f, E_m: 3e9, nu_m: 0.35, E_i: 1.4e11, nu_i: 0.35, cte_m: 6e-5, cte_i: 4.3e-6 }, kpi: 'cte_eff',
    expected: kerner, tol: 1e-6, kind: 'analytical',
    reference: 'Kerner, elastic and thermo-elastic properties of composite media (Proc. Phys. Soc.); Levin relation; Schapery bounds'
  },
  // ---------- numerical RVE ----------
  {
    id: 's18-rve-checkerboard-k', suite: 18, title: 'Numerical RVE, checkerboard: Dykhne result k* = √(k₁k₂)', mode: 'rve',
    inputs: { mode: 'rve', micro: 'checker', nres: '64', k_m: 1, k_i: 2 }, kpi: 'k_iso',
    expected: Math.SQRT2, tol: 0.01, kind: 'analytical',
    reference: 'Dykhne, conductivity of a two-dimensional two-phase system (Sov. Phys. JETP); Keller duality theorem',
    note: 'Finite-volume cell problem; the corner singularities of the checkerboard converge slowly, so moderate contrast (2:1) is used.'
  },
  {
    id: 's18-rve-hill-equal-shear-C11', suite: 18, title: 'Numerical RVE, phases of equal shear modulus: C11 = ⟨1/(λ + 2μ)⟩⁻¹ (Hill)', mode: 'rve',
    inputs: { mode: 'rve', micro: 'checker', nres: '32', E_m: 3e9, nu_m: 0.35, E_i: E_B_eq, nu_i: 0.2 }, kpi: 'C11',
    expected: hillC11, tol: 0.005, kind: 'analytical',
    reference: 'Hill, elastic properties of reinforced solids — exact result for composites with uniform shear modulus (J. Mech. Phys. Solids)',
    note: 'Exact for any microstructure when both phases share the same shear modulus; checkerboard used because its phase fraction is exactly 0.5.'
  },
  {
    id: 's18-rve-hill-equal-shear-C66', suite: 18, title: 'Numerical RVE, phases of equal shear modulus: effective shear modulus unchanged', mode: 'rve',
    inputs: { mode: 'rve', micro: 'checker', nres: '32', E_m: 3e9, nu_m: 0.35, E_i: E_B_eq, nu_i: 0.2 }, kpi: 'C66',
    expected: muEq, tol: 0.005, kind: 'analytical',
    reference: 'Hill, elastic properties of reinforced solids — exact result for composites with uniform shear modulus (J. Mech. Phys. Solids)'
  },
  // ---------- laminate ----------
  {
    id: 's18-laminate-E-in', suite: 18, title: 'Bilayer laminate, in-plane modulus with Poisson mismatch', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.4, E_m: 7e10, nu_m: 0.35, E_i: 2e11, nu_i: 0.2 }, kpi: 'E_in',
    expected: lamEin, tol: 1e-6, kind: 'analytical',
    reference: 'Classical laminate theory for bonded isotropic layers (Jones, Mechanics of Composite Materials): E∥ = ΣfE + f_Af_BE_AE_B(ν_A − ν_B)²/(f_AE_A(1 − ν_B²) + f_BE_B(1 − ν_A²))',
    note: 'Exceeds the plain rule of mixtures because the layers constrain each other\'s Poisson contraction.'
  },
  {
    id: 's18-laminate-E-out', suite: 18, title: 'Bilayer laminate, through-thickness modulus for ν = 0: Reuss (series) mean', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.25, E_m: 3e9, nu_m: 0, E_i: 1.4e11, nu_i: 0 }, kpi: 'E_out',
    expected: 1 / (0.75 / 3e9 + 0.25 / 1.4e11), tol: 1e-6, kind: 'analytical',
    reference: 'Backus averaging of a layered medium; Reuss iso-stress limit'
  },
  {
    id: 's18-laminate-G-out', suite: 18, title: 'Bilayer laminate, interlaminar shear modulus: harmonic mean', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.4, E_m: 7e10, nu_m: 0.35, E_i: 2e11, nu_i: 0.2 }, kpi: 'G_out',
    expected: 1 / (0.6 / (7e10 / 2.7) + 0.4 / (2e11 / 2.4)), tol: 1e-6, kind: 'analytical',
    reference: 'Backus, long-wave elastic anisotropy produced by horizontal layering (J. Geophys. Res.)'
  },
  {
    id: 's18-laminate-k-in', suite: 18, title: 'Layered medium, in-plane conductivity: arithmetic mean', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.5, k_m: 170, k_i: 30 }, kpi: 'k_in',
    expected: 100, tol: 1e-9, kind: 'analytical', reference: 'Wiener bounds — layers conducting in parallel'
  },
  {
    id: 's18-laminate-k-out', suite: 18, title: 'Layered medium, cross-plane conductivity with Kapitza interfaces (series resistances)', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.5, k_m: 170, k_i: 30, period_nm: 10, n_bi: 40, R_int: 3e-9 }, kpi: 'k_out',
    expected: lamKout, tol: 1e-6, kind: 'analytical',
    reference: 'Wiener series bound with thermal boundary resistance added in series (Cahill et al., nanoscale thermal transport review)',
    note: '40 bilayers have 79 internal interfaces.'
  },
  {
    id: 's18-laminate-cte-in', suite: 18, title: 'Bilayer laminate, in-plane thermal expansion (biaxial-modulus weighting)', mode: 'laminate',
    inputs: { mode: 'laminate', vf: 0.4, E_m: 7e10, nu_m: 0.35, cte_m: 2.3e-5, E_i: 2e11, nu_i: 0.2, cte_i: 1.2e-5 }, kpi: 'cte_in',
    expected: lamCte, tol: 1e-6, kind: 'analytical',
    reference: 'Equi-biaxial thermal mismatch of bonded layers (Freund & Suresh, Thin Film Materials)'
  },
  // ---------- lattice-to-continuum ----------
  {
    id: 's18-octet-E100', suite: 18, title: 'Octet truss: E[100]/E_s = ρ̄/9', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'octet', rel: 0.05, E_m: 1e11 }, kpi: 'E100',
    expected: 0.05 * 1e11 / 9, tol: 1e-6, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, Effective properties of the octet-truss lattice material (J. Mech. Phys. Solids)',
    note: 'Pin-jointed slender-strut limit.'
  },
  {
    id: 's18-octet-E111', suite: 18, title: 'Octet truss: E[111]/E_s = ρ̄/5', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'octet', rel: 0.05, E_m: 1e11 }, kpi: 'E111',
    expected: 0.05 * 1e11 / 5, tol: 1e-6, kind: 'analytical',
    reference: 'Deshpande, Fleck & Ashby, Effective properties of the octet-truss lattice material (J. Mech. Phys. Solids)'
  },
  {
    id: 's18-simple-cubic-E100', suite: 18, title: 'Simple-cubic lattice: E[100]/E_s = ρ̄/3', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'sc', rel: 0.06, E_m: 1e11 }, kpi: 'E100',
    expected: 0.06 * 1e11 / 3, tol: 1e-6, kind: 'analytical',
    reference: 'Affine pin-jointed strut model: one third of the struts are aligned with the load (Gibson & Ashby, Cellular Solids)'
  },
  {
    id: 's18-kelvin-K', suite: 18, title: 'Kelvin-cell open foam: bulk modulus K/E_s = ρ̄/9 (strut stretching under hydrostatic load)', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'kelvin', rel: 0.02, E_m: 1e11 }, kpi: 'K',
    expected: 0.02 * 1e11 / 9, tol: 1e-6, kind: 'analytical',
    reference: 'Zhu, Knott & Mills, Analysis of the elastic properties of open-cell foams with tetrakaidecahedral cells (J. Mech. Phys. Solids); Warren & Kraynik',
    note: 'Uniform dilatation is an affine deformation of the strut network, so the struts carry axial force only.'
  },
  {
    id: 's18-regular-honeycomb-E', suite: 18, title: 'Regular hexagonal honeycomb: E/E_s = (4/√3)(t/l)³ = 1.5ρ̄³', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'reentrant', theta: 30, h_l: 1, rel: 0.05, E_m: 7e10 }, kpi: 'E1',
    expected: 1.5 * 0.05 ** 3 * 7e10, tol: 1e-6, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model for hexagonal honeycombs',
    note: 'Slender walls: t/l = (√3/2)ρ̄; axial and shear wall deformation neglected in both suite and reference.'
  },
  {
    id: 's18-regular-honeycomb-G', suite: 18, title: 'Regular hexagonal honeycomb: G/E_s = (1/√3)(t/l)³ (isotropic, ν = 1)', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'reentrant', theta: 30, h_l: 1, rel: 0.05, E_m: 7e10 }, kpi: 'G12',
    expected: 0.375 * 0.05 ** 3 * 7e10, tol: 1e-6, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model for hexagonal honeycombs'
  },
  {
    id: 's18-reentrant-nu', suite: 18, title: 'Re-entrant honeycomb, θ = −30°, h/l = 2: ν₁₂ = −1 (auxetic)', mode: 'lattice',
    inputs: { mode: 'lattice', topo: 'reentrant', theta: -30, h_l: 2, rel: 0.05 }, kpi: 'nu12',
    expected: -1, tol: 1e-9, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — ν₁₂ = cos²θ/((h/l + sin θ) sin θ); Masters & Evans re-entrant honeycomb'
  },
  // ---------- single crystal → polycrystal ----------
  {
    id: 's18-si-vrh-E', suite: 18, title: 'Silicon: Voigt–Reuss–Hill Young\'s modulus from C11 = 165.7, C12 = 63.9, C44 = 79.6 GPa', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 165.7, C12: 63.9, C44: 79.6, rho_x: 2329, formula: 'Si' }, kpi: 'E',
    expected: si.E * 1e9, tol: 1e-6, kind: 'analytical',
    reference: 'Hill, The elastic behaviour of a crystalline aggregate (Proc. Phys. Soc.); cubic closed forms for G_V and G_R',
    note: 'Evaluates to about 163 GPa, the commonly quoted polycrystalline silicon modulus.'
  },
  {
    id: 's18-si-vrh-nu', suite: 18, title: 'Silicon: Voigt–Reuss–Hill Poisson ratio', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 165.7, C12: 63.9, C44: 79.6, rho_x: 2329, formula: 'Si' }, kpi: 'nu',
    expected: 0.22, tol: 0.02, kind: 'reference-data',
    reference: 'Hopcroft, Nix & Kenny, What is the Young\'s modulus of silicon? (J. Microelectromech. Syst.) — polycrystalline average ν ≈ 0.22',
    note: 'Literature value quoted to two figures.'
  },
  {
    id: 's18-si-longitudinal-speed', suite: 18, title: 'Silicon polycrystal: longitudinal wave speed √((K + 4G/3)/ρ)', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 165.7, C12: 63.9, C44: 79.6, rho_x: 2329, formula: 'Si' }, kpi: 'v_l',
    expected: si.vl, tol: 1e-6, kind: 'analytical',
    reference: 'Isotropic elastodynamics with Hill-averaged moduli (Auld, Acoustic Fields and Waves in Solids)'
  },
  {
    id: 's18-si-debye-temperature', suite: 18, title: 'Silicon: Debye temperature from elastic constants', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 165.7, C12: 63.9, C44: 79.6, rho_x: 2329, formula: 'Si' }, kpi: 'theta_D',
    expected: 645, tol: 0.03, kind: 'reference-data',
    reference: 'Kittel, Introduction to Solid State Physics — table of Debye temperatures (Si: 645 K)',
    note: 'Room-temperature elastic constants and the isotropic mean sound speed (Anderson method) versus the low-temperature calorimetric value.'
  },
  {
    id: 's18-zno-hex-K', suite: 18, title: 'Wurtzite ZnO: Hill bulk modulus from hexagonal closed forms', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'hex', C11: 209.7, C12: 121.1, C13: 105.1, C33: 210.9, C44: 42.5, rho_x: 5606, formula: 'ZnO' }, kpi: 'K',
    expected: zno.K * 1e9, tol: 1e-6, kind: 'analytical',
    reference: 'Voigt and Reuss averages for hexagonal crystals (Nye, Physical Properties of Crystals; Hill)',
    note: 'K_R = ((C11 + C12)C33 − 2C13²)/(C11 + C12 + 2C33 − 4C13).'
  },
  {
    id: 's18-zno-hex-G', suite: 18, title: 'Wurtzite ZnO: Hill shear modulus from hexagonal closed forms', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'hex', C11: 209.7, C12: 121.1, C13: 105.1, C33: 210.9, C44: 42.5, rho_x: 5606, formula: 'ZnO' }, kpi: 'G',
    expected: zno.G * 1e9, tol: 1e-6, kind: 'analytical',
    reference: 'Voigt and Reuss averages for hexagonal crystals (Nye, Physical Properties of Crystals; Hill)'
  },
  {
    id: 's18-isotropic-anisotropy-index', suite: 18, title: 'Isotropic stiffness (C44 = (C11 − C12)/2): universal anisotropy index is zero', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 200, C12: 100, C44: 50 }, kpi: 'A_U',
    expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Ranganathan & Ostoja-Starzewski, Universal elastic anisotropy index (Phys. Rev. Lett.)'
  },
  {
    id: 's18-lj-argon-bulk-modulus', suite: 18, title: 'Lennard-Jones argon (σ = 3.405 Å, ε/k_B = 119.8 K), static fcc: B = 4εA₁₂(A₆/A₁₂)^(5/2)/σ³ ≈ 75.2 ε/σ³', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'lj', eps_lj: lj.epsEV, sig_lj: 3.405, formula: 'Ar' }, kpi: 'K',
    expected: lj.K, tol: 0.01, kind: 'analytical',
    reference: 'Ashcroft & Mermin, Solid State Physics — cohesive energy and bulk modulus of the noble-gas solids (fcc lattice sums A₁₂ = 12.13, A₆ = 14.45)',
    note: 'The suite truncates the lattice sums at six nearest-neighbour distances, which omits about 0.2 % of the attractive sum.'
  },
  {
    id: 's18-lj-argon-density', suite: 18, title: 'Lennard-Jones argon, static fcc: density from the equilibrium spacing d = 1.09σ', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'lj', eps_lj: lj.epsEV, sig_lj: 3.405, formula: 'Ar' }, kpi: 'rho',
    expected: lj.rho, tol: 0.005, kind: 'analytical',
    reference: 'Ashcroft & Mermin, Solid State Physics — nearest-neighbour distance (2A₁₂/A₆)^(1/6) σ of the fcc Lennard-Jones crystal',
    note: 'Corresponds to a lattice parameter of about 5.25 Å (measured solid argon: about 5.3 Å including zero-point expansion).'
  },
  // ---- published experiments ----
  {
    id: 's18-exp-smith-glass-epoxy-e-015', suite: 18, title: 'Glass-bead/epoxy composite, 14.9 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.1494 }, kpi: 'E_eff', expected: 4.19e9, tol: 0.05, kind: 'experiment',
    measured: { value: 4.19, uncertainty: 0.14, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 2 (volume fraction 0.1494). Phase constants as given in the source: matrix E = 3.01 GPa and ν = 0.394 (section 5; the three unfilled castings of table 2 measured 2.99–3.03 GPa), glass E = 76 GPa (manufacturer, section 2) and ν = 0.23 (literature value quoted in section 2). Uncertainty is the probable error of a single observation. The suite default for spheres, the differential scheme, gives 4.14 GPa (−1.3 %); the tolerance is the specimen scatter (probable error 3 %) with a margin for the phase constants. Mori–Tanaka is 3 % low here.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-e-030', suite: 18, title: 'Glass-bead/epoxy composite, 30.3 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.3025 }, kpi: 'E_eff', expected: 6.01e9, tol: 0.05, kind: 'experiment',
    measured: { value: 6.01, uncertainty: 0.21, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 2 (volume fraction 0.3025). Phase constants as given in the source: matrix E = 3.01 GPa and ν = 0.394 (section 5; the three unfilled castings of table 2 measured 2.99–3.03 GPa), glass E = 76 GPa (manufacturer, section 2) and ν = 0.23 (literature value quoted in section 2). Uncertainty is the probable error of a single observation. The suite default for spheres, the differential scheme, gives 5.98 GPa (−0.5 %); the tolerance is the specimen scatter (probable error 3.5 %) with a margin for the phase constants. The source itself shows (figure 1) that the Hashin–Shtrikman lower bound, which the Mori–Tanaka estimate reproduces for spheres, falls increasingly below the data as the filler fraction rises (8 % here); that estimate was the suite default before and remains selectable.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-e-040', suite: 18, title: 'Glass-bead/epoxy composite, 39.6 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.396 }, kpi: 'E_eff', expected: 7.94e9, tol: 0.06, kind: 'experiment',
    measured: { value: 7.94, uncertainty: 0.12, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 2 (volume fraction 0.3960); phase constants as for the entries above. The differential scheme (suite default for spheres) gives 7.70 GPa (−3.0 %); the tolerance covers specimen scatter and the phase constants. For comparison at this loading: generalised self-consistent −7 %, Lielens +5 %, Mori–Tanaka −15 %, symmetric self-consistent +26 %.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-e-050', suite: 18, title: 'Glass-bead/epoxy composite, 49.6 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.496 }, kpi: 'E_eff', expected: 12.09e9, tol: 0.15, kind: 'experiment',
    measured: { value: 12.09, uncertainty: 0.09, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 2 (volume fraction 0.4960); phase constants as for the entries above. The differential scheme gives 10.36 GPa (−14.3 %): at this loading the aged composite is markedly stiffer than every matrix-continuous mean-field estimate (generalised self-consistent −20 %, Mori–Tanaka −30 %), which the suite states in a warning above 45 vol %; the tolerance is that documented model deviation, not measurement scatter. The freshly cast series of the same source at 49.8 vol % is reproduced within 5 % (next entries), so the discrepancy is specific to this point rather than systematic.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-e-050-lielens', suite: 18, title: 'Glass-bead/epoxy composite, 49.6 vol % filler: Lielens interpolation against the measured modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', est: 'lielens', poros: 0, vf: 0.496 }, kpi: 'E_eff', expected: 12.09e9, tol: 0.08, kind: 'experiment',
    measured: { value: 12.09, uncertainty: 0.09, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens',
    match: ['vf', 'E_m', 'E_i', 'est'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Same measurement as the previous entry, compared with the selectable Lielens interpolation (11.42 GPa, −5.6 %), the estimate the suite points to above 45 vol %. The interpolation weight (f + f²)/2 is the published one and was not adjusted.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-fresh-e-045', suite: 18, title: 'Glass-bead/epoxy composite tested fresh, 45.1 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 2.68e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.4513 }, kpi: 'E_eff', expected: 8.17e9, tol: 0.06, kind: 'experiment',
    measured: { value: 8.17, uncertainty: 0.23, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested 12–14 days after casting; mean of four specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 1 (volume fraction 0.4513, tested after 12 days). Matrix constants for the fresh series as given in section 5 of the source: E = 2.68 GPa (mean of the unfilled castings of table 1, which scatter between 2.48 and 3.21 GPa) and ν = 0.394. The differential scheme gives 8.18 GPa (+0.1 %); the tolerance reflects the specimen scatter (3 %) and the uncertainty of the fresh-matrix modulus.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-fresh-e-050', suite: 18, title: 'Glass-bead/epoxy composite tested fresh, 49.8 vol % filler: measured Young\'s modulus', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 2.68e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.4981 }, kpi: 'E_eff', expected: 9.95e9, tol: 0.1, kind: 'experiment',
    measured: { value: 9.95, uncertainty: 0.9, unit: 'GPa' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested 12–14 days after casting; mean of three specimens',
    match: ['vf', 'E_m', 'E_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 1 (volume fraction 0.4981, tested after 14 days; probable error 9 %). Matrix constants for the fresh series as given in section 5 of the source (E = 2.68 GPa, ν = 0.394). The differential scheme gives 9.47 GPa (−4.8 %); the tolerance is the quoted probable error of the measurement.'
  },
  {
    id: 's18-exp-smith-glass-epoxy-nu-030', suite: 18, title: 'Glass-bead/epoxy composite, 30.3 vol % filler: measured Poisson\'s ratio', mode: 'bounds',
    inputs: { mode: 'bounds', E_m: 3.01e9, nu_m: 0.394, E_i: 7.6e10, nu_i: 0.23, shape: 'sphere', poros: 0, vf: 0.3025 }, kpi: 'nu_eff', expected: 0.375, tol: 0.06, kind: 'experiment',
    measured: { value: 0.375, uncertainty: 0.017, unit: '–' },
    source: { authors: 'Smith, J. C.', title: 'Experimental values for the elastic constants of a particulate-filled glassy polymer', journal: 'J. Res. Natl. Bur. Stand. A', year: 1976, doi: '10.6028/jres.080A.008', url: 'https://nvlpubs.nist.gov/nistpubs/jres/80A/jresv80An1p45_A1b.pdf' },
    conditions: 'Glass spheres (1–30 µm, soda-lime crown glass, silane-coupled) in DGEBA/TETA epoxy; dumbbell tensile specimens at 0.05 cm/min, 24 °C, tested about 200 days after casting (stored at 54 % R.H.); mean of four specimens; transverse strain by LVDT extensometer simultaneously with the axial strain',
    match: ['vf', 'nu_m', 'nu_i'],
    reference: 'Smith, J. Res. Natl. Bur. Stand. 80A, 45 (1976)',
    note: 'Measured value from table 2 (volume fraction 0.3025). Phase constants as given in the source: matrix E = 3.01 GPa and ν = 0.394 (section 5; the three unfilled castings of table 2 measured 2.99–3.03 GPa), glass E = 76 GPa (manufacturer, section 2) and ν = 0.23 (literature value quoted in section 2). Uncertainty is the probable error of a single observation. The source notes considerable scatter in the measured Poisson\'s ratios; the tolerance is the quoted probable error (about 5 %). The differential scheme gives 0.360 (−4 %).'
  },
  {
    id: 's18-exp-copper-polycrystal-e', suite: 18, title: 'Copper: polycrystal Young\'s modulus from measured single-crystal constants', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 169.1, C12: 122.2, C44: 75.42, rho_x: 8932, formula: 'Cu' }, kpi: 'E', expected: 123.5e9, tol: 0.05, kind: 'experiment',
    measured: { value: 123.5, uncertainty: 0.7, unit: 'GPa' },
    source: { authors: 'Ledbetter, H. M. & Naimon, E. R.', title: 'Elastic properties of metals and alloys. II. Copper', journal: 'J. Phys. Chem. Ref. Data', year: 1974, doi: '10.1063/1.3253150', url: 'https://srd.nist.gov/jpcrdreprint/1.3253150.pdf' },
    conditions: 'Annealed polycrystalline copper at room temperature; mean of the published static and dynamic measurements compiled in the source (outliers beyond 1.5 σ excluded), uncertainty is the standard error of the mean',
    match: ['C11', 'C12', 'C44'],
    reference: 'Ledbetter & Naimon, J. Phys. Chem. Ref. Data 3, 897 (1974)',
    note: 'Measured polycrystal value from table 10 (E = 1.235 ± 0.007 × 10¹¹ N/m²). Inputs are the recommended ultrasonic single-crystal stiffnesses of the same table (C11 = 169.1, C12 = 122.2, C44 = 75.42 GPa), measured independently on single crystals; density 8932 kg/m³ from table 1. The source reports that the arithmetic Voigt–Reuss–Hill average (127.3 GPa) overestimates the measured polycrystal modulus of this strongly anisotropic metal (Zener ratio 3.2) by about 3 %, which the tolerance reflects.'
  },
  {
    id: 's18-exp-copper-polycrystal-vl', suite: 18, title: 'Copper: measured longitudinal ultrasonic velocity of the polycrystal', mode: 'atomistic',
    inputs: { mode: 'atomistic', src: 'cij', sym: 'cubic', C11: 169.1, C12: 122.2, C44: 75.42, rho_x: 8932, formula: 'Cu' }, kpi: 'v_l', expected: 4790, tol: 0.03, kind: 'experiment',
    measured: { value: 0.479e6, uncertainty: null, unit: 'cm/s' },
    source: { authors: 'Ledbetter, H. M. & Naimon, E. R.', title: 'Elastic properties of metals and alloys. II. Copper', journal: 'J. Phys. Chem. Ref. Data', year: 1974, doi: '10.1063/1.3253150', url: 'https://srd.nist.gov/jpcrdreprint/1.3253150.pdf' },
    conditions: 'Polycrystalline copper, ultrasonic pulse measurement at 2.5 MHz, room temperature (Haskins 1961, as compiled in the source)',
    match: ['C11', 'C12', 'C44', 'rho_x'],
    reference: 'Ledbetter & Naimon, J. Phys. Chem. Ref. Data 3, 897 (1974) — entry Haskins (1961)',
    note: 'Measured value from the comments column of table 9 (Haskins 1961: longitudinal velocity 0.479 × 10⁶ cm/s). Inputs are the recommended single-crystal stiffnesses of table 10 and the density of table 1; the specimen density of the velocity measurement is not given in the source. The velocity depends on the Hill-averaged K + 4G/3, so the 3 % shear-modulus bias of the average enters with a weight of about one sixth.'
  },
];
