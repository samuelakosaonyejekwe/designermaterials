// Reference benchmarks for Suite 19 — large-scale component and system analysis.
// Every expected value is a literal or is computed here from a closed-form expression with standard constants;
// nothing in this file calls the suite.
const PI = Math.PI, C0 = 299792458, deg = 180 / PI;

// ---- lattice panels (frame finite elements) ----
const Es = 7e10, cellW = 0.01;                       // member modulus (Pa), cell width (m)
const tThin = 1e-4;                                  // slender wall used for the unit-cell comparisons (m)
const lHex = cellW / (2 * Math.cos(PI / 6));         // regular hexagon wall length for a cell width of 10 mm
const thRe = -20 / deg, hlRe = 2;                    // re-entrant cell: θ = −20°, h/l = 2
const lTri = cellW;                                  // triangular lattice member length

// ---- Kirchhoff plate ----
const Dplate = (E, h, nu) => E * h ** 3 / (12 * (1 - nu * nu));
const steel = { E_p: 2e11, nu_p: 0.3, rho_p: 7850, ph: 0.01 }, Dsteel = Dplate(steel.E_p, steel.ph, steel.nu_p);
const aux = { E_p: 2e9, nu_p: -0.5, rho_p: 400, ph: 0.01 }, Daux = Dplate(aux.E_p, aux.ph, aux.nu_p);

// ---- resonator chain ----
const kS = 1e5, mLight = 0.01, mHeavy = 0.03;
const rM = Math.sqrt(mLight / mHeavy), cosMid = 3 - 2 * (rM + 1 / rM);      // cos(qa) of the diatomic chain at ω² = 2k/√(m₁m₂)
const kappaMid = Math.acosh(Math.abs(cosMid));

// ---- aperture ----
const fGHz10 = 1e10, lam = C0 / fGHz10, pitch = 0.015, dLam = pitch / lam;

// ---- reactor (cylinder with two flat ends; water charge; 316L wall) ----
const Vr = 0.05, HD = 1.5, Dv = Math.cbrt(4 * Vr / (PI * HD)), Av = PI * Dv * (HD * Dv) + PI * Dv * Dv / 2;
const kWall = 16, rhoWall = 8000, cpWall = 500, tWall = 0.008, kIns = 0.04;
const Rout = ins => tWall / (2 * kWall * Av) + ins / (kIns * Av) + 1 / (10 * Av);        // wall half + insulation + external film (h = 10 W/m²/K)
const Cliq = 1000 * 4186 * Vr, Cwall = rhoWall * cpWall * Av * tWall;
const RinFast = 1 / (1e4 * Av) + tWall / (2 * kWall * Av), tauLag = RinFast * Cliq * Cwall / (Cliq + Cwall);

// ---- reliability ----
const gamma15 = Math.sqrt(PI) / 2;                   // Γ(3/2)

export default [
  // ------------------------------------------------------------------ finite lattice panels
  {
    id: 's19-square-grid-modulus', suite: 19, title: 'Square grid in tension along its members: effective modulus', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'square', lcase: 'tension', ncx: 8, ncy: 5, cell: cellW, t_strut: 5e-4, depth: 0.01, E_s: Es, strain: 0.2 },
    kpi: 'eff_fe', expected: Es * 5e-4 / cellW, tol: 1e-4, kind: 'analytical',
    reference: 'Hand statics: parallel bars in uniaxial tension, E* = E_s·t/a (Gibson & Ashby, Cellular Solids — square cells loaded along the walls)',
    note: 'Exact for the frame model: the vertical members stretch uniformly and the horizontal ones stay unloaded.',
  },
  {
    id: 's19-square-grid-member-stress', suite: 19, title: 'Square grid in tension: member stress equals E_s·ε', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'square', lcase: 'tension', ncx: 8, ncy: 5, cell: cellW, t_strut: 5e-4, depth: 0.01, E_s: Es, strain: 0.2 },
    kpi: 's_max', expected: Es * 0.002, tol: 1e-4, kind: 'analytical',
    reference: 'Hooke\'s law for an axially loaded bar, σ = E·ε',
  },
  {
    id: 's19-square-grid-euler-margin', suite: 19, title: 'Square grid in compression: pinned-end Euler margin of the members', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'square', lcase: 'compression', ncx: 8, ncy: 5, cell: cellW, t_strut: 5e-4, depth: 0.01, E_s: Es, strain: 0.2 },
    kpi: 'buck', expected: PI ** 2 * (5e-4 / cellW) ** 2 / (12 * 0.002), tol: 1e-4, kind: 'analytical',
    reference: 'Euler column formula P_cr = π²EI/L² for a rectangular section, divided by the member force E·A·ε',
    note: 'P_cr/N = π²(t/L)²/(12ε) for a pinned–pinned strut of rectangular section.',
  },
  {
    id: 's19-triangular-lattice-modulus', suite: 19, title: 'Triangular (stretch-dominated) lattice: in-plane Young\'s modulus', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'triangular', lcase: 'tension', ncx: 20, ncy: 12, cell: cellW, t_strut: tThin, depth: 0.01, E_s: Es, strain: 0.01 },
    kpi: 'eff_fe', expected: 2 / Math.sqrt(3) * Es * tThin / lTri, tol: 0.02, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — equilateral-triangle cells: E* = (2/√3)·E_s·t/l (pin-jointed, axial deformation only)',
    note: '20 × 12 cells, t/l = 0.01. The reference ignores wall bending (a relative contribution of order (t/l)²) and the free side edges of the finite panel.',
  },
  {
    id: 's19-honeycomb-poisson', suite: 19, title: 'Regular hexagonal honeycomb panel: in-plane Poisson\'s ratio', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'honeycomb', lcase: 'tension', ncx: 20, ncy: 12, cell: cellW, t_strut: tThin, depth: 0.01, E_s: Es, strain: 0.01 },
    kpi: 'nu_fe', expected: 1.0, tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — cell-wall bending model for regular hexagonal honeycombs (ν* = 1)',
    note: 'Slender walls (t/l ≈ 0.017): axial wall deformation, neglected in the reference, lowers ν by about (t/l)².',
  },
  {
    id: 's19-reentrant-poisson', suite: 19, title: 'Re-entrant (auxetic) honeycomb panel: negative Poisson\'s ratio', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'reentrant', lcase: 'tension', th_re: 20, h_l: hlRe, ncx: 20, ncy: 12, cell: cellW, t_strut: tThin, depth: 0.01, E_s: Es, strain: 0.01 },
    kpi: 'nu_fe', expected: (hlRe + Math.sin(thRe)) * Math.sin(thRe) / Math.cos(thRe) ** 2, tol: 0.01, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — ν₂₁ = (h/l + sin θ)·sin θ / cos²θ with θ < 0 for the re-entrant cell',
    note: 'θ = −20°, h/l = 2 gives ν = −0.642. Bending-only reference; measured on the free edges in the middle half of a 20 × 12-cell panel.',
  },
  {
    id: 's19-honeycomb-modulus', suite: 19, title: 'Regular hexagonal honeycomb panel: in-plane Young\'s modulus', mode: 'lattice_fe',
    inputs: { mode: 'lattice_fe', ctype: 'honeycomb', lcase: 'tension', ncx: 20, ncy: 12, cell: cellW, t_strut: tThin, depth: 0.01, E_s: Es, strain: 0.01 },
    kpi: 'eff_fe', expected: 4 / Math.sqrt(3) * Es * (tThin / lHex) ** 3, tol: 0.07, kind: 'analytical',
    reference: 'Gibson & Ashby, Cellular Solids — regular hexagonal honeycomb, E* = (4/√3)·E_s·(t/l)³',
    note: 'The reference is the infinite lattice. A finite panel is softer because the wall columns on its two free zig-zag edges are only half supported; the deficit is close to one cell width in Nx (about 5 % at the maximum of 20 cells, 10 % at 10 cells) and is what the suite reports as its finite-size effect.',
  },
  // ------------------------------------------------------------------ homogenised plate
  {
    id: 's19-plate-square-deflection', suite: 19, title: 'Simply supported square plate under uniform pressure: centre deflection', mode: 'plate',
    inputs: { mode: 'plate', pa: 1, pb: 1, q_p: 1000, Nx_p: 0, ...steel },
    kpi: 'w_max', expected: 0.00406 * 1000 * 1 ** 4 / Dsteel, tol: 0.005, kind: 'analytical',
    reference: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells — Navier solution, w_max = 0.00406·q·a⁴/D',
    note: 'The tabulated coefficient is rounded to three significant figures.',
  },
  {
    id: 's19-plate-square-stress', suite: 19, title: 'Simply supported square plate (ν = 0.3): maximum bending stress', mode: 'plate',
    inputs: { mode: 'plate', pa: 1, pb: 1, q_p: 1000, Nx_p: 0, ...steel },
    kpi: 's_max', expected: 6 * 0.0479 * 1000 * 1 ** 2 / steel.ph ** 2, tol: 0.005, kind: 'analytical',
    reference: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells — centre moment M = 0.0479·q·a² for ν = 0.3; σ = 6M/h²',
    note: 'The tabulated coefficient is rounded to three significant figures.',
  },
  {
    id: 's19-plate-2to1-deflection', suite: 19, title: 'Simply supported 2:1 rectangular plate: centre deflection', mode: 'plate',
    inputs: { mode: 'plate', pa: 2, pb: 1, q_p: 1000, Nx_p: 0, ...steel },
    kpi: 'w_max', expected: 0.01013 * 1000 * 1 ** 4 / Dsteel, tol: 0.005, kind: 'analytical',
    reference: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells — w_max = 0.01013·q·b⁴/D for a/b = 2 (b the shorter side)',
  },
  {
    id: 's19-plate-square-frequency', suite: 19, title: 'Simply supported square plate: fundamental frequency', mode: 'plate',
    inputs: { mode: 'plate', pa: 1, pb: 1, q_p: 1000, Nx_p: 0, ...steel },
    kpi: 'f1', expected: (2 * PI ** 2 / 1 ** 2) * Math.sqrt(Dsteel / (steel.rho_p * steel.ph)) / (2 * PI), tol: 1e-4, kind: 'analytical',
    reference: 'Leissa, Vibration of Plates — simply supported rectangular plate, ω₁₁ = (2π²/a²)·√(D/ρh) for a square',
  },
  {
    id: 's19-plate-square-buckling', suite: 19, title: 'Simply supported square plate: uniaxial buckling load', mode: 'plate',
    inputs: { mode: 'plate', pa: 1, pb: 1, q_p: 1000, Nx_p: 1e4, ...steel },
    kpi: 'N_cr', expected: 4 * PI ** 2 * Dsteel / 1 ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Timoshenko & Gere, Theory of Elastic Stability — N_cr = 4π²D/b² (buckling coefficient k = 4)',
  },
  {
    id: 's19-plate-long-buckling', suite: 19, title: 'Simply supported 3:1 plate compressed along its length: buckling load', mode: 'plate',
    inputs: { mode: 'plate', pa: 3, pb: 1, q_p: 1000, Nx_p: 1e4, ...steel },
    kpi: 'N_cr', expected: 4 * PI ** 2 * Dsteel / 1 ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Timoshenko & Gere, Theory of Elastic Stability — k = 4 whenever a/b is an integer (three half-waves here)',
  },
  {
    id: 's19-plate-auxetic-deflection', suite: 19, title: 'Auxetic square plate (ν = −0.5): centre deflection', mode: 'plate',
    inputs: { mode: 'plate', pa: 0.5, pb: 0.5, q_p: 2000, Nx_p: 0, ...aux },
    kpi: 'w_max', expected: 0.00406 * 2000 * 0.5 ** 4 / Daux, tol: 0.005, kind: 'analytical',
    reference: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells — w_max = 0.00406·q·a⁴/D with D = Eh³/[12(1 − ν²)]',
    note: 'The deflection coefficient of the simply supported plate does not depend on ν; a negative ν enters only through the rigidity D, which is 1/(1 − ν²) = 1.333 times that of a ν = 0 plate.',
  },
  // ------------------------------------------------------------------ resonator chain
  {
    id: 's19-chain-gap-lower', suite: 19, title: 'Diatomic chain: lower band-gap edge (top of the acoustic branch)', mode: 'modal',
    inputs: { mode: 'modal', n_cells: 10, m1: mLight, m2: mHeavy, k_s: kS },
    kpi: 'gap_lo', expected: Math.sqrt(2 * kS / mHeavy) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — diatomic linear chain, ω² = 2C/M₁ and 2C/M₂ at the zone boundary',
  },
  {
    id: 's19-chain-gap-upper', suite: 19, title: 'Diatomic chain: upper band-gap edge (bottom of the optical branch)', mode: 'modal',
    inputs: { mode: 'modal', n_cells: 10, m1: mLight, m2: mHeavy, k_s: kS },
    kpi: 'gap_hi', expected: Math.sqrt(2 * kS / mLight) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Kittel, Introduction to Solid State Physics — diatomic linear chain, ω² = 2C/M₁ and 2C/M₂ at the zone boundary',
  },
  {
    id: 's19-chain-midgap-attenuation', suite: 19, title: 'Diatomic chain: evanescent decay per cell at mid-gap', mode: 'modal',
    inputs: { mode: 'modal', n_cells: 10, m1: mLight, m2: mHeavy, k_s: kS },
    kpi: 'att_cell', expected: 20 / Math.LN10 * kappaMid, tol: 1e-3, kind: 'analytical',
    reference: 'Brillouin, Wave Propagation in Periodic Structures — complex Bloch wavenumber of the diatomic lattice inside the stop band',
    note: 'At ω² = 2k/√(m₁m₂) the dispersion relation gives cos(qa) = 3 − 2(r + 1/r) with r = √(m₁/m₂); attenuation = (20/ln 10)·acosh|cos qa| dB per cell.',
  },
  {
    id: 's19-chain-monatomic-fundamental', suite: 19, title: 'Fixed–free uniform chain of 20 masses: fundamental frequency', mode: 'modal',
    inputs: { mode: 'modal', n_cells: 10, m1: mLight, m2: mLight, k_s: kS },
    kpi: 'f1', expected: 2 * Math.sqrt(kS / mLight) * Math.sin(PI / (2 * (2 * 20 + 1))) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Closed-form eigenvalues of the fixed–free spring–mass chain, ω_j = 2√(k/m)·sin[(2j − 1)π / (2(2N + 1))] (Meirovitch, Fundamentals of Vibrations)',
  },
  {
    id: 's19-chain-monatomic-highest', suite: 19, title: 'Fixed–free uniform chain of 20 masses: highest natural frequency', mode: 'modal',
    inputs: { mode: 'modal', n_cells: 10, m1: mLight, m2: mLight, k_s: kS },
    kpi: 'f_max', expected: 2 * Math.sqrt(kS / mLight) * Math.sin(39 * PI / (2 * (2 * 20 + 1))) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Closed-form eigenvalues of the fixed–free spring–mass chain, ω_j = 2√(k/m)·sin[(2j − 1)π / (2(2N + 1))] (Meirovitch, Fundamentals of Vibrations)',
  },
  // ------------------------------------------------------------------ finite aperture
  {
    id: 's19-array-uniform-beamwidth', suite: 19, title: 'Uniform 40-element broadside array: half-power beamwidth', mode: 'array',
    inputs: { mode: 'array', nx: 40, ny: 1, pitch, freq: fGHz10, steer: 0, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'hpbw', expected: 0.886 * lam / (40 * pitch) * deg, tol: 0.01, kind: 'analytical',
    reference: 'Balanis, Antenna Theory — uniform linear array, HPBW ≈ 0.886·λ/(N·d) rad at broadside',
    note: 'Large-N approximation of the reference; the suite reads the −3 dB points from a pattern sampled every 0.25°.',
  },
  {
    id: 's19-array-uniform-sidelobe', suite: 19, title: 'Uniform 40-element array: first sidelobe level', mode: 'array',
    inputs: { mode: 'array', nx: 40, ny: 1, pitch, freq: fGHz10, steer: 0, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'sll', expected: -13.26, abs: 0.1, kind: 'analytical',
    reference: 'Balanis, Antenna Theory — first sidelobe of a uniformly illuminated line source, −13.26 dB',
    note: 'The continuous-aperture value; a discrete 40-element array is about 0.02 dB higher and the 0.25° pattern sampling can miss the lobe peak by a few hundredths of a dB.',
  },
  {
    id: 's19-array-aperture-directivity', suite: 19, title: 'Error-free 16 × 16 aperture at broadside: directivity 4πA/λ²', mode: 'array',
    inputs: { mode: 'array', nx: 16, ny: 16, pitch, freq: fGHz10, steer: 0, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'gain', expected: 10 * Math.log10(4 * PI * 256 * dLam ** 2), abs: 0.01, kind: 'analytical',
    reference: 'Balanis, Antenna Theory — maximum directivity of a uniformly illuminated aperture, D = 4πA/λ²',
    note: 'Pitch 15 mm at 10 GHz (0.5003 λ): 29.06 dBi.',
  },
  {
    id: 's19-array-phase-gradient', suite: 19, title: 'Beam steering to 30°: progressive inter-element phase', mode: 'array',
    inputs: { mode: 'array', nx: 16, ny: 16, pitch, freq: fGHz10, steer: 30, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'grad', expected: -360 * dLam * 0.5, tol: 1e-5, kind: 'analytical',
    reference: 'Phased-array steering condition Δφ = −k·d·sin θ₀ (Mailloux, Phased Array Antenna Handbook)',
  },
  {
    id: 's19-array-steered-beam', suite: 19, title: 'Error-free array steered to 30°: realised beam direction', mode: 'array',
    inputs: { mode: 'array', nx: 16, ny: 16, pitch, freq: fGHz10, steer: 30, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'beam', expected: 30, abs: 0.25, kind: 'analytical',
    reference: 'Phased-array steering condition Δφ = −k·d·sin θ₀ (Mailloux, Phased Array Antenna Handbook)',
    note: 'Absolute tolerance equals the 0.25° sampling of the computed pattern.',
  },
  {
    id: 's19-array-ruze-loss', suite: 19, title: 'Gain loss for 15° rms random phase error', mode: 'array',
    inputs: { mode: 'array', nx: 16, ny: 16, pitch, freq: fGHz10, steer: 0, taper: 'uniform', phase_err: 15, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'loss_ruze', expected: 10 * Math.log10(Math.E) * (15 / deg) ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Ruze, antenna tolerance theory — G/G₀ = exp(−σ²) for rms phase error σ in radians',
  },
  {
    id: 's19-array-grating-lobe-scan', suite: 19, title: 'Grating-lobe-free scan limit for a 20 mm pitch at 10 GHz', mode: 'array',
    inputs: { mode: 'array', nx: 16, ny: 16, pitch: 0.02, freq: fGHz10, steer: 0, taper: 'uniform', phase_err: 0, amp_err: 0, bits: 0, fail: 0 },
    kpi: 'scan_max', expected: Math.asin(lam / 0.02 - 1) * deg, tol: 1e-5, kind: 'analytical',
    reference: 'Grating-lobe condition d/λ < 1/(1 + |sin θ₀|) (Balanis, Antenna Theory; Mailloux, Phased Array Antenna Handbook)',
  },
  // ------------------------------------------------------------------ scaled reactor
  {
    id: 's19-reactor-psat-200C', suite: 19, title: 'Autogenous pressure of a water-filled vessel at 200 °C', mode: 'reactor',
    inputs: { mode: 'reactor', T_set: 200, P_heat: 5e4 },
    kpi: 'P_sat', expected: 1.5547e6, tol: 0.01, kind: 'reference-data',
    reference: 'Steam tables (IAPWS-IF97): saturation pressure of water at 200 °C = 1.5547 MPa',
    note: 'The suite uses the Antoine correlation, which is accurate to a few tenths of a per cent here.',
  },
  {
    id: 's19-reactor-psat-100C', suite: 19, title: 'Saturation pressure of water at the normal boiling point', mode: 'reactor',
    inputs: { mode: 'reactor', T_set: 100 },
    kpi: 'P_sat', expected: 101325, tol: 0.01, kind: 'reference-data',
    reference: 'Definition of the normal boiling point of water (100 °C at 101.325 kPa; steam tables give 101.42 kPa at 100.00 °C)',
    note: 'Antoine correlation; 1 % covers the correlation error at the junction of its two temperature ranges.',
  },
  {
    id: 's19-reactor-holding-power', suite: 19, title: '50 L insulated vessel at 90 °C: steady holding power', mode: 'reactor',
    inputs: { mode: 'reactor', V_r: Vr, HD, T_set: 90, T_amb: 25, wall_mat: 'steel', wall_mm: 8, ins_mm: 50, h_ext: 10 },
    kpi: 'P_hold', expected: (90 - 25) / Rout(0.05), tol: 1e-4, kind: 'analytical',
    reference: 'Series thermal resistances Q = ΔT / (t_w/2k_w A + t_ins/k_ins A + 1/h_o A) (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer)',
    note: 'Closed cylinder H/D = 1.5 with both end caps (A = πDH + πD²/2), insulation conductivity 0.04 W/m/K, plane-wall resistances on the vessel surface area.',
  },
  {
    id: 's19-reactor-heatup-time', suite: 19, title: 'Heavily insulated, well-stirred 50 L vessel: heat-up time', mode: 'reactor',
    inputs: { mode: 'reactor', V_r: Vr, HD, T_set: 90, T_amb: 25, P_heat: 5000, wall_mat: 'steel', wall_mm: 8, ins_mm: 500, h_int: 1e4, h_ext: 10 },
    kpi: 't_heat', expected: (Cliq + Cwall) * (90 - 25) / 5000 + tauLag, tol: 0.005, kind: 'analytical',
    reference: 'First-law energy balance of a lumped system, t = (C_liquid + C_wall)·ΔT/P, plus the asymptotic wall-to-liquid lag R_i·C_f·C_w/(C_f + C_w) of the two-capacitance network',
    note: 'Water 1000 kg/m³ and 4186 J/kg/K, 316L wall 8000 kg/m³ and 500 J/kg/K. With 0.5 m of insulation the heat loss is below 0.1 % of the heater power, which is the only term the reference omits.',
  },
  {
    id: 's19-reactor-hoop-stress', suite: 19, title: 'Thin-walled vessel at 200 °C: hoop stress from the gauge pressure', mode: 'reactor',
    inputs: { mode: 'reactor', V_r: Vr, HD, T_set: 200, P_heat: 5e4, wall_mm: 8 },
    kpi: 'hoop', expected: (1.5547e6 - 101325) * Dv / (2 * tWall), tol: 0.01, kind: 'analytical',
    reference: 'Thin-walled cylinder σ_θ = p·D/(2t) with the steam-table saturation pressure of water at 200 °C',
  },
  // ------------------------------------------------------------------ reliability block diagram
  {
    id: 's19-reliability-series-weibull', suite: 19, title: '100 Weibull stages in series with an exponential common element', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 100, m_par: 1, beta_w: 2.5, eta_h: 2e5, beta_c: 1, eta_c: 5e5, t_miss: 8760 },
    kpi: 'R', expected: Math.exp(-100 * (8760 / 2e5) ** 2.5 - 8760 / 5e5), tol: 1e-6, kind: 'analytical',
    reference: 'Series-system reliability as the product of component reliabilities with R_c = exp[−(t/η)^β] (O\'Connor & Kleyner, Practical Reliability Engineering)',
  },
  {
    id: 's19-reliability-weakest-link', suite: 19, title: 'Weakest-link size effect: characteristic life of 100 units in series', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 100, m_par: 1, beta_w: 2.5, eta_h: 2e5 },
    kpi: 'eta_ser', expected: 2e5 * 100 ** (-1 / 2.5), tol: 1e-6, kind: 'analytical',
    reference: 'Weibull weakest-link scaling η_n = η·n^(−1/β)',
  },
  {
    id: 's19-reliability-series-b10', suite: 19, title: 'B10 life of 100 Weibull units in series', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 100, m_par: 1, beta_w: 2.5, eta_h: 2e5, beta_c: 1, eta_c: 1e9 },
    kpi: 'B10', expected: 2e5 * 100 ** (-1 / 2.5) * (-Math.log(0.9)) ** (1 / 2.5), tol: 1e-3, kind: 'analytical',
    reference: 'Weibull quantile t_p = η_n·[−ln(1 − p)]^(1/β) with the weakest-link scale η_n = η·n^(−1/β)',
    note: 'The common element is given its longest admissible life (10⁹ h); its residual contribution to the hazard is of order 10⁻⁴.',
  },
  {
    id: 's19-reliability-parallel-mttf', suite: 19, title: 'Two exponential components in active parallel: MTTF = 1.5/λ', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 1, m_par: 2, beta_w: 1, eta_h: 1e4, beta_c: 1, eta_c: 1e9, t_miss: 1000 },
    kpi: 'mttf', expected: 1.5e4, tol: 2e-3, kind: 'analytical',
    reference: 'Active-parallel redundancy of two identical exponential units, MTTF = 3/(2λ) (O\'Connor & Kleyner, Practical Reliability Engineering)',
    note: 'Trapezoidal integration of R(t) on 400 points, truncated where R falls to about 10⁻⁶.',
  },
  {
    id: 's19-reliability-parallel-R', suite: 19, title: 'Two exponential components in active parallel: reliability at t = η', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 1, m_par: 2, beta_w: 1, eta_h: 1e4, beta_c: 1, eta_c: 1e9, t_miss: 1e4 },
    kpi: 'R', expected: 1 - (1 - Math.exp(-1)) ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Parallel-system reliability R = 1 − (1 − R_c)^m (O\'Connor & Kleyner, Practical Reliability Engineering)',
    note: 'R = 2e⁻¹ − e⁻² = 0.6004; the 10⁹ h common element lowers it by 10⁻⁵ relative.',
  },
  {
    id: 's19-reliability-rayleigh-mttf', suite: 19, title: 'Single Weibull component with β = 2: mean life η·Γ(1 + 1/β)', mode: 'reliability',
    inputs: { mode: 'reliability', n_ser: 1, m_par: 1, beta_w: 2, eta_h: 1e4, beta_c: 1, eta_c: 1e9, t_miss: 1000 },
    kpi: 'mttf', expected: 1e4 * gamma15, tol: 2e-3, kind: 'analytical',
    reference: 'Mean of the Weibull distribution, MTTF = η·Γ(1 + 1/β), with Γ(3/2) = √π/2',
  },
];
