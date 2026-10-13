// Reference benchmarks for Suite 15 — Acoustic, Phononic and Elastic-Wave Engineering.
// Every expected value is a literal handbook figure or is evaluated here from a closed-form expression; nothing is taken from the suite.
const PI = Math.PI, KB = 1.380649e-23;
const AIR = { rho: 1.204, c: 343, mu: 1.81e-5 }; // air at 20 °C, 1 atm
const hz = w => w / (2 * PI);

// ---- mass-in-mass lattice
const MM = { m1: 0.05, m2: 0.1, k1: 1e5, k2: 4e4, a: 0.02 }, mmW0 = Math.sqrt(MM.k2 / MM.m2);
const mmB = (MM.m1 + MM.m2) * mmW0 ** 2 + 4 * MM.k1, mmD = Math.sqrt(mmB ** 2 - 16 * MM.m1 * MM.k1 * mmW0 ** 2); // band edges at qa = π
const mmIn = { mode: 'massinmass', m1: 0.05, m2: 0.1, k1: 1e5, k2: 4e4, a_cell: 20, eta: 0 };

// ---- Helmholtz-resonator / membrane duct in air
const HR = { Sd: PI * 0.015 ** 2, rn: 5e-3, ln: 5e-3, V: 30e-6, d: 0.04, T: 300, sig: 0.1 };
const hrW = AIR.c * Math.sqrt(PI * HR.rn ** 2 / (HR.V * (HR.ln + 1.7 * HR.rn))), hrF = HR.V / (HR.Sd * HR.d);
const hrRho = AIR.rho + HR.sig / HR.d, hrWc = Math.sqrt(8 * PI * HR.T / (hrRho * HR.Sd * HR.d));
const hrIn = { mode: 'helmholtz', fluid: 'air', duct_D: 30, neck_r: 5, neck_l: 5, cav_V: 30, d_cell: 40, memb_T: 300, memb_sig: 0.1, eta: 0 };
const hrKrel = x => 1 / (1 - hrF / (x * x - 1)), hrRrel = x => hrRho / AIR.rho * (1 - (hrWc / (x * hrW)) ** 2); // lossless, x = ω/ω_H

// ---- quarter-wave steel/epoxy stack (handbook-typical ρ, c_L): equal transit time 0.1 µs per layer
const ST = { rho: 7850, c: 5900 }, EP = { rho: 1180, c: 2540 }, zS = ST.rho * ST.c, zE = EP.rho * EP.c, qwR = (zS - zE) / (zS + zE), qwFB = 1 / (2 * (590e-6 / ST.c + 254e-6 / EP.c));
const qwIn = { mode: 'phononic', matA: 'steel', matB: 'epoxy', dA: 590, dB: 254, eta: 0, Ncell: 4 };
const tA = 1e-3 / ST.c, tB = 1e-3 / EP.c; // default 1 mm / 1 mm stack

// ---- resonator-loaded panel
const PN = { E: 70e9, rho: 2700, nu: 0.33 }, rc = AIR.rho * AIR.c;
const tlRes = (() => { const ms = 2.7, mu = 0.3, eta = 0.03, x = 2 * PI * 300 * ms / (2 * rc); return 10 * Math.log10((1 + x * mu / eta) ** 2 + (x * (1 + mu)) ** 2); })();

// ---- micro-perforated panel (Maa): independent evaluation of the absorption peak
const maa = (f, d, t, sg, D) => { const w = 2 * PI * f, k = d * Math.sqrt(w * AIR.rho / (4 * AIR.mu)), r = 32 * AIR.mu * t / (sg * AIR.rho * AIR.c * d * d) * (Math.sqrt(1 + k * k / 32) + Math.SQRT2 / 32 * k * d / t), x = w * t / (sg * AIR.c) * (1 + 1 / Math.sqrt(9 + k * k / 2) + 0.85 * d / t) - 1 / Math.tan(w * D / AIR.c); return 4 * r / ((1 + r) ** 2 + x * x); };
const maaPeak = (() => { const g = f => maa(f, 0.3e-3, 0.3e-3, 0.01, 0.03); let lo = 300, hi = 3000; for (let i = 0; i < 200; i++) { const a = lo + (hi - lo) / 3, b = hi - (hi - lo) / 3; if (g(a) < g(b)) lo = a; else hi = b; } return { f: 0.5 * (lo + hi), a: g(0.5 * (lo + hi)) }; })();

// ---- silicon (C₁₁ = 165.7, C₁₂ = 63.9, C₄₄ = 79.6 GPa, ρ = 2329 kg/m³)
const SI = { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9, rho: 2329 };
const siGV = (SI.c11 - SI.c12 + 3 * SI.c44) / 5, siGR = 5 * (SI.c11 - SI.c12) * SI.c44 / (4 * SI.c44 + 3 * (SI.c11 - SI.c12)), siK = (SI.c11 + 2 * SI.c12) / 3, siG = 0.5 * (siGV + siGR);

// ---- silicon nanobeam 10 µm × 300 nm × 100 nm
const NB = { E: 169e9, rho: 2329, L: 10e-6, w: 300e-9, t: 100e-9 }, nbA = NB.w * NB.t, nbI = NB.w * NB.t ** 3 / 12, nbR = Math.sqrt(NB.E * nbI / (NB.rho * nbA));
const nbIn = { mode: 'nanobeam', E_b: 169e9, rho_b: 2329, section: 'rect', L_b: 10, w_b: 300, t_b: 100, sig_r: 0, Es_b: 0, T_b: 25 };
const zW = 998 * 1482;

export default [
  // ------------------------------------------------------------------ mass-in-mass lattice (negative effective mass)
  {
    id: 's15-mim-local-resonance', suite: 15, title: 'Mass-in-mass cell: local resonance √(k₂/m₂)/2π', mode: 'massinmass',
    inputs: mmIn, kpi: 'f0', expected: hz(mmW0), tol: 1e-6, kind: 'analytical',
    reference: 'Huang, Sun & Huang, On the negative effective mass density in acoustic metamaterials',
  },
  {
    id: 's15-mim-gap-upper-edge', suite: 15, title: 'Mass-in-mass lattice: upper gap edge ω₂√(1 + m₂/m₁)', mode: 'massinmass',
    inputs: mmIn, kpi: 'gap_hi', expected: hz(mmW0 * Math.sqrt(1 + MM.m2 / MM.m1)), tol: 1e-6, kind: 'analytical',
    reference: 'Huang, Sun & Huang, On the negative effective mass density in acoustic metamaterials — m_eff = 0 at ω₂√(1 + m₂/m₁)',
  },
  {
    id: 's15-mim-gap-lower-edge', suite: 15, title: 'Mass-in-mass lattice: lower gap edge (acoustic branch at qa = π)', mode: 'massinmass',
    inputs: mmIn, kpi: 'gap_lo', expected: hz(Math.sqrt((mmB - mmD) / (2 * MM.m1))), tol: 1e-6, kind: 'analytical',
    reference: 'Dispersion relation m_eff(ω)ω² = 2k₁(1 − cos qa) of the mass-in-mass chain (Huang, Sun & Huang)',
  },
  {
    id: 's15-mim-negative-effective-mass', suite: 15, title: 'Mass-in-mass cell: m_eff = m₁ + m₂ω₂²/(ω₂² − ω²) is negative at 1.08 f₀', mode: 'massinmass',
    inputs: { ...mmIn, f_rel: 1.08 }, kpi: 'm_eff', expected: MM.m1 + MM.m2 / (1 - 1.08 ** 2), tol: 1e-6, kind: 'analytical',
    reference: 'Huang, Sun & Huang, On the negative effective mass density in acoustic metamaterials',
    note: 'Lossless case (η = 0).',
  },
  {
    id: 's15-mim-long-wave-speed', suite: 15, title: 'Mass-in-mass lattice: long-wave speed a√(k₁/(m₁ + m₂))', mode: 'massinmass',
    inputs: mmIn, kpi: 'c_lw', expected: MM.a * Math.sqrt(MM.k1 / (MM.m1 + MM.m2)), tol: 1e-6, kind: 'analytical',
    reference: 'Brillouin, Wave Propagation in Periodic Structures — monatomic chain with the total static mass',
  },
  // ------------------------------------------------------------------ Helmholtz resonators and membranes (negative modulus / density)
  {
    id: 's15-helmholtz-resonance', suite: 15, title: 'Helmholtz resonator: f₀ = (c/2π)√(S/(V·l_eff)), l_eff = l + 1.7r', mode: 'helmholtz',
    inputs: hrIn, kpi: 'fH', expected: hz(hrW), tol: 1e-6, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics — Helmholtz resonator with flanged end corrections (0.85r per end)',
  },
  {
    id: 's15-helmholtz-negative-modulus-edge', suite: 15, title: 'Helmholtz-resonator duct: negative-modulus band ends at f_H√(1 + F)', mode: 'helmholtz',
    inputs: { ...hrIn, dn_type: 'hr' }, kpi: 'negK_hi', expected: hz(hrW) * Math.sqrt(1 + hrF), tol: 0.01, kind: 'analytical',
    reference: 'Fang et al., Ultrasonic metamaterials with negative modulus — Lorentz-type effective compressibility',
    note: 'The suite locates the band on a 400-point frequency sweep (spacing about 0.5 % of this frequency) and includes viscous neck loss; the reference is lossless.',
  },
  {
    id: 's15-helmholtz-negative-modulus-value', suite: 15, title: 'Helmholtz-resonator duct: K_eff/K₀ = 1/(1 − F/(x² − 1)) at x = 1.2', mode: 'helmholtz',
    inputs: { ...hrIn, dn_type: 'hr', f_rel: 1.2 }, kpi: 'K_rel', expected: hrKrel(1.2), tol: 0.01, kind: 'analytical',
    reference: 'Fang et al., Ultrasonic metamaterials with negative modulus — Lorentz-type effective compressibility',
    note: 'Lossless reference; the suite adds thermoviscous neck damping (Q ≈ 57), which changes the real part slightly away from resonance.',
  },
  {
    id: 's15-membrane-cutoff', suite: 15, title: 'Membrane-loaded duct: negative-density cut-off ω_c = √(8πT/(ρ′Sd))', mode: 'helmholtz',
    inputs: hrIn, kpi: 'fc', expected: hz(hrWc), tol: 1e-6, kind: 'analytical',
    reference: 'Lee et al., Acoustic metamaterial with negative density — Drude-type density ρ′(1 − ω_c²/ω²) of a duct with clamped membranes',
  },
  {
    id: 's15-membrane-negative-density', suite: 15, title: 'Membrane-loaded duct: ρ_eff/ρ₀ = (ρ′/ρ₀)(1 − ω_c²/ω²) below cut-off', mode: 'helmholtz',
    inputs: { ...hrIn, f_rel: 1.2 }, kpi: 'rho_rel', expected: hrRrel(1.2), tol: 1e-3, kind: 'analytical',
    reference: 'Lee et al., Acoustic metamaterial with negative density',
    note: 'Evaluated at 1.2 f_H with the membrane loss switched off (the suite keeps a residual damping of 10⁻⁴ ω_c).',
  },
  {
    id: 's15-double-negative-phase-velocity', suite: 15, title: 'Double-negative duct: backward wave, v_p = −c₀√(K_rel/ρ_rel) at 1.2 f_H', mode: 'helmholtz',
    inputs: { ...hrIn, f_rel: 1.2 }, kpi: 'vp', expected: -AIR.c * Math.sqrt(hrKrel(1.2) / hrRrel(1.2)), tol: 0.01, kind: 'analytical',
    reference: 'Lee et al., Composite acoustic medium with simultaneously negative density and modulus; Veselago sign convention for double-negative media',
    note: 'Lossless reference against the suite value with neck damping.',
  },
  // ------------------------------------------------------------------ layered phononic crystal
  {
    id: 's15-phononic-bragg-frequency', suite: 15, title: 'Quarter-wave steel/epoxy stack: Bragg frequency 1/(2(d_A/c_A + d_B/c_B)) = 2.5 MHz', mode: 'phononic',
    inputs: qwIn, kpi: 'fB', expected: 2.5e6, tol: 1e-6, kind: 'analytical',
    reference: 'Rytov dispersion relation for a periodically layered medium (Brekhovskikh, Waves in Layered Media)',
  },
  {
    id: 's15-phononic-quarter-wave-gap-width', suite: 15, title: 'Quarter-wave stack: relative gap width (4/π)·asin(|Z_A − Z_B|/(Z_A + Z_B))', mode: 'phononic',
    inputs: qwIn, kpi: 'gap_w', expected: 4 / PI * Math.asin(qwR), tol: 1e-4, kind: 'analytical',
    reference: 'Quarter-wave stack stop-band width (Yeh, Optical Waves in Layered Media; identical for acoustic impedances)',
  },
  {
    id: 's15-phononic-quarter-wave-gap-edge', suite: 15, title: 'Quarter-wave stack: lower gap edge f_B(1 − (2/π)·asin r)', mode: 'phononic',
    inputs: qwIn, kpi: 'gap_lo', expected: qwFB * (1 - 2 / PI * Math.asin(qwR)), tol: 1e-4, kind: 'analytical',
    reference: 'Quarter-wave stack stop-band edges (Yeh, Optical Waves in Layered Media)',
  },
  {
    id: 's15-phononic-midgap-attenuation', suite: 15, title: 'Quarter-wave stack: mid-gap attenuation 8.686·ln(Z_A/Z_B) dB per period', mode: 'phononic',
    inputs: qwIn, kpi: 'att', expected: 20 * Math.log10(zS / zE), tol: 0.002, kind: 'analytical',
    reference: 'Bloch wavenumber at the centre of a quarter-wave stop band, cosh(Im qΛ) = ½(Z_A/Z_B + Z_B/Z_A) (Yeh, Optical Waves in Layered Media)',
    note: 'The suite samples the Bloch wavenumber at the sweep point nearest the gap centre.',
  },
  {
    id: 's15-phononic-interface-reflection', suite: 15, title: 'Steel/epoxy interface: power reflection ((Z_A − Z_B)/(Z_A + Z_B))²', mode: 'phononic',
    inputs: qwIn, kpi: 'R', expected: qwR ** 2, tol: 1e-6, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics — normal-incidence reflection between two media',
  },
  {
    id: 's15-phononic-effective-speed', suite: 15, title: 'Layered medium, long-wave limit: c = 1/√(⟨ρ⟩⟨1/ρc²⟩)', mode: 'phononic',
    inputs: { mode: 'phononic', matA: 'steel', matB: 'epoxy', dA: 1000, dB: 1000 }, kpi: 'c_eff',
    expected: 1 / Math.sqrt(0.5 * (ST.rho + EP.rho) * 0.5 * (1 / (ST.rho * ST.c ** 2) + 1 / (EP.rho * EP.c ** 2))), tol: 1e-6, kind: 'analytical',
    reference: 'Effective-medium (Backus / Rytov) limit of a periodically layered medium',
  },
  // ------------------------------------------------------------------ sound-insulating panel
  {
    id: 's15-mass-law-normal-incidence', suite: 15, title: 'Mass law, normal incidence: TL = 10·log₁₀[1 + (πfm″/ρ₀c₀)²]', mode: 'panel',
    inputs: { mode: 'panel', rho_p: 2700, h_p: 4, mu_r: 0.01, f_res: 1000 }, kpi: 'TL_mass', expected: 10 * Math.log10(1 + (PI * 1000 * 10.8 * 1.01 / rc) ** 2), tol: 1e-6, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics — transmission through a limp wall',
  },
  {
    id: 's15-mass-law-handbook-form', suite: 15, title: 'Mass law, handbook form: TL ≈ 20·log₁₀(f·m″) − 42.5 dB', mode: 'panel',
    inputs: { mode: 'panel', rho_p: 2700, h_p: 4, mu_r: 0.01, f_res: 1000 }, kpi: 'TL_mass', expected: 20 * Math.log10(1000 * 10.8 * 1.01) - 42.5, tol: 0.01, kind: 'reference-data',
    reference: 'Normal-incidence mass law as tabulated in noise-control handbooks (Bies & Hansen, Engineering Noise Control)',
    note: 'The handbook constant is rounded and assumes ρ₀c₀ ≈ 415 Pa·s/m; the suite uses 413 Pa·s/m (air at 20 °C).',
  },
  {
    id: 's15-coincidence-frequency', suite: 15, title: '1 mm aluminium plate: coincidence frequency (c₀²/2π)√(m″/D)', mode: 'panel',
    inputs: { mode: 'panel', E_p: 70e9, rho_p: 2700, nu_p: 0.33, h_p: 1 }, kpi: 'fc', expected: AIR.c ** 2 / (2 * PI) * Math.sqrt(2.7 / (PN.E * 1e-9 / (12 * (1 - PN.nu ** 2)))), tol: 1e-6, kind: 'analytical',
    reference: 'Cremer, Heckl & Petersson, Structure-Borne Sound — critical frequency of a thin plate (about 12 kHz·mm for aluminium)',
  },
  {
    id: 's15-panel-negative-effective-mass', suite: 15, title: 'Resonator-loaded panel: m″_eff = m″[1 + µ/(1 − (f/f_r)²)] is negative at 1.08 f_r', mode: 'panel',
    inputs: { mode: 'panel', rho_p: 2700, h_p: 1, mu_r: 0.3, f_res: 300, eta: 0, f_rel: 1.08 }, kpi: 'm_eff', expected: 2.7 * (1 + 0.3 / (1 - 1.08 ** 2)), tol: 1e-4, kind: 'analytical',
    reference: 'Locally resonant plate with spring–mass resonators (Xiao, Wen & Wen, Sound transmission loss of metamaterial-based thin plates)',
    note: 'The suite enforces a minimum resonator loss factor of 10⁻⁴.',
  },
  {
    id: 's15-panel-tl-at-resonance', suite: 15, title: 'Resonator-loaded panel: TL at the tuned frequency with hysteretic damping', mode: 'panel',
    inputs: { mode: 'panel', rho_p: 2700, h_p: 1, mu_r: 0.3, f_res: 300, eta: 0.03 }, kpi: 'TL_res', expected: tlRes, tol: 1e-4, kind: 'analytical',
    reference: 'Limp-wall transmission τ = |1 + iωm″_eff/(2ρ₀c₀)|⁻² (Kinsler et al.) with m″_eff = m″(1 + µ − iµ/η) at resonance',
  },
  // ------------------------------------------------------------------ micro-perforated absorber
  {
    id: 's15-mpp-peak-absorption', suite: 15, title: 'Micro-perforated panel (d = t = 0.3 mm, 1 %, 30 mm cavity): peak absorption', mode: 'mpp',
    inputs: { mode: 'mpp', d_hole: 0.3, t_mpp: 0.3, perf: 1, D_cav: 30 }, kpi: 'a_pk', expected: maaPeak.a, tol: 1e-3, kind: 'analytical',
    reference: 'Maa, Potential of microperforated panel absorber — approximate impedance formulae, α = 4r/[(1 + r)² + (ωm − cot(ωD/c₀))²]',
    note: 'Expected value from an independent maximisation of the Maa formulae in this file (air at 20 °C).',
  },
  {
    id: 's15-mpp-peak-frequency', suite: 15, title: 'Micro-perforated panel: frequency of peak absorption', mode: 'mpp',
    inputs: { mode: 'mpp', d_hole: 0.3, t_mpp: 0.3, perf: 1, D_cav: 30 }, kpi: 'f_pk', expected: maaPeak.f, tol: 1e-3, kind: 'analytical',
    reference: 'Maa, Potential of microperforated panel absorber',
  },
  {
    id: 's15-quarter-wave-tube', suite: 15, title: 'Quarter-wave tube: fundamental c₀/(4(L + 0.85r))', mode: 'mpp',
    inputs: { mode: 'mpp', L_qw: 1000, r_qw: 20, phi_qw: 20 }, kpi: 'f_qw', expected: AIR.c / (4 * (1 + 0.85 * 0.02)), tol: 0.02, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics — closed pipe with flanged open-end correction',
    note: 'The suite adds thermoviscous wall loss (lowers the peak by about 1 % for this wide tube) and reads the peak from a logarithmic sweep with 1.3 % spacing.',
  },
  // ------------------------------------------------------------------ Christoffel waves in crystals
  {
    id: 's15-si-100-longitudinal', suite: 15, title: 'Silicon [100]: longitudinal wave speed 8433 m/s', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: 90, azi_ang: 0 }, kpi: 'vL', expected: 8433, tol: 0.001, kind: 'reference-data',
    reference: 'McSkimin, Measurement of elastic constants of silicon; tabulated in Auld, Acoustic Fields and Waves in Solids',
  },
  {
    id: 's15-si-100-shear', suite: 15, title: 'Silicon [100]: shear wave speed 5843 m/s', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: 90, azi_ang: 0 }, kpi: 'vT1', expected: 5843, tol: 0.001, kind: 'reference-data',
    reference: 'McSkimin, Measurement of elastic constants of silicon; tabulated in Auld, Acoustic Fields and Waves in Solids',
  },
  {
    id: 's15-si-110-longitudinal', suite: 15, title: 'Silicon [110]: longitudinal wave speed 9134 m/s', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: 90, azi_ang: 45 }, kpi: 'vL', expected: 9134, tol: 0.001, kind: 'reference-data',
    reference: 'McSkimin, Measurement of elastic constants of silicon; tabulated in Auld, Acoustic Fields and Waves in Solids',
  },
  {
    id: 's15-si-110-fast-shear', suite: 15, title: 'Silicon [110]: fast shear ([001]-polarised) 5843 m/s', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: 90, azi_ang: 45 }, kpi: 'vT2', expected: 5843, tol: 0.001, kind: 'reference-data',
    reference: 'McSkimin, Measurement of elastic constants of silicon; tabulated in Auld, Acoustic Fields and Waves in Solids',
  },
  {
    id: 's15-si-110-slow-shear', suite: 15, title: 'Silicon [110]: slow shear ([1−10]-polarised) 4673 m/s', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: 90, azi_ang: 45 }, kpi: 'vT1', expected: 4673, tol: 0.001, kind: 'reference-data',
    reference: 'McSkimin, Measurement of elastic constants of silicon; tabulated in Auld, Acoustic Fields and Waves in Solids',
  },
  {
    id: 's15-si-111-longitudinal', suite: 15, title: 'Silicon [111]: v_L = √((C₁₁ + 2C₁₂ + 4C₄₄)/3ρ)', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: Math.acos(1 / Math.sqrt(3)) * 180 / PI, azi_ang: 45 }, kpi: 'vL', expected: Math.sqrt((SI.c11 + 2 * SI.c12 + 4 * SI.c44) / (3 * SI.rho)), tol: 1e-6, kind: 'analytical',
    reference: 'Christoffel equation along the cube diagonal of a cubic crystal (Auld, Acoustic Fields and Waves in Solids)',
  },
  {
    id: 's15-si-111-shear', suite: 15, title: 'Silicon [111]: degenerate shear speed √((C₁₁ − C₁₂ + C₄₄)/3ρ)', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si', pol_ang: Math.acos(1 / Math.sqrt(3)) * 180 / PI, azi_ang: 45 }, kpi: 'vT1', expected: Math.sqrt((SI.c11 - SI.c12 + SI.c44) / (3 * SI.rho)), tol: 1e-6, kind: 'analytical',
    reference: 'Christoffel equation along the cube diagonal of a cubic crystal (Auld, Acoustic Fields and Waves in Solids)',
  },
  {
    id: 's15-si-zener-ratio', suite: 15, title: 'Silicon: Zener anisotropy ratio 2C₄₄/(C₁₁ − C₁₂) = 1.56', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si' }, kpi: 'A', expected: 1.56, tol: 0.005, kind: 'reference-data',
    reference: 'Hirth & Lothe, Theory of Dislocations — anisotropy ratios of cubic crystals',
  },
  {
    id: 's15-cu-zener-ratio', suite: 15, title: 'Copper: Zener anisotropy ratio 3.21', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Cu' }, kpi: 'A', expected: 3.21, tol: 0.005, kind: 'reference-data',
    reference: 'Hirth & Lothe, Theory of Dislocations — anisotropy ratios of cubic crystals',
  },
  {
    id: 's15-si-hill-shear-modulus', suite: 15, title: 'Silicon polycrystal: Voigt–Reuss–Hill shear modulus', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si' }, kpi: 'G', expected: siG, tol: 1e-6, kind: 'analytical',
    reference: 'Hill, The elastic behaviour of a crystalline aggregate — G_V = (C₁₁ − C₁₂ + 3C₄₄)/5, G_R = 5(C₁₁ − C₁₂)C₄₄/[4C₄₄ + 3(C₁₁ − C₁₂)]',
  },
  {
    id: 's15-si-hill-young-modulus', suite: 15, title: 'Silicon polycrystal: Hill Young modulus 9KG/(3K + G) ≈ 163 GPa', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Si' }, kpi: 'E', expected: 9 * siK * siG / (3 * siK + siG), tol: 1e-6, kind: 'analytical',
    reference: 'Hill, The elastic behaviour of a crystalline aggregate',
  },
  {
    id: 's15-custom-cubic-diamond', suite: 15, title: 'Custom cubic crystal (diamond constants): [100] longitudinal speed √(C₁₁/ρ)', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'custom', cc11: 1079, cc12: 124, cc44: 578, rho_x: 3515, pol_ang: 90, azi_ang: 0 }, kpi: 'vL', expected: Math.sqrt(1079e9 / 3515), tol: 1e-6, kind: 'analytical',
    reference: 'Christoffel equation along a cube axis (Auld, Acoustic Fields and Waves in Solids); diamond constants after McSkimin & Andreatch',
  },
  // ------------------------------------------------------------------ nanobeam resonators and interfaces
  {
    id: 's15-clamped-beam-f1', suite: 15, title: 'Clamped–clamped silicon bridge: f₁ = (4.730²/2π)√(EI/ρAL⁴)', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'fixed' }, kpi: 'f1', expected: 4.730041 ** 2 / (2 * PI) * nbR / NB.L ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Blevins, Formulas for Natural Frequency and Mode Shape — clamped–clamped beam, λ₁ = 4.7300',
  },
  {
    id: 's15-clamped-beam-f2', suite: 15, title: 'Clamped–clamped silicon bridge: second mode, λ₂ = 7.8532', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'fixed' }, kpi: 'f2', expected: 7.853205 ** 2 / (2 * PI) * nbR / NB.L ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Blevins, Formulas for Natural Frequency and Mode Shape — clamped–clamped beam',
  },
  {
    id: 's15-cantilever-f1', suite: 15, title: 'Cantilever: f₁ = (1.8751²/2π)√(EI/ρAL⁴)', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'cantilever' }, kpi: 'f1', expected: 1.875104 ** 2 / (2 * PI) * nbR / NB.L ** 2, tol: 1e-4, kind: 'analytical',
    reference: 'Blevins, Formulas for Natural Frequency and Mode Shape — clamped–free beam, λ₁ = 1.8751',
  },
  {
    id: 's15-cantilever-stiffness', suite: 15, title: 'Cantilever: effective stiffness equals the static tip stiffness Ewt³/(4L³)', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'cantilever' }, kpi: 'k_eff', expected: NB.E * NB.w * NB.t ** 3 / (4 * NB.L ** 3), tol: 0.001, kind: 'analytical',
    reference: 'Euler–Bernoulli cantilever, k = 3EI/L³ (Sader et al., calibration of rectangular atomic force microscope cantilevers: m_eff ≈ 0.2427 m)',
  },
  {
    id: 's15-tensioned-pinned-beam', suite: 15, title: 'Pinned–pinned beam under 100 MPa tension: f₁ = (π/2L²)√(EI/ρA)·√(1 + σAL²/π²EI)', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'pinned', sig_r: 100 }, kpi: 'f1', expected: PI / (2 * NB.L ** 2) * nbR * Math.sqrt(1 + 100e6 * nbA * NB.L ** 2 / (PI ** 2 * NB.E * nbI)), tol: 1e-6, kind: 'analytical',
    reference: 'Exact solution for a simply supported beam with axial load (Timoshenko, Vibration Problems in Engineering)',
  },
  {
    id: 's15-thermomechanical-amplitude', suite: 15, title: 'Pinned–pinned beam at 25 °C: thermal amplitude √(k_BT/k_eff), k_eff = π⁴EI/(2L³)', mode: 'nanobeam',
    inputs: { ...nbIn, bc_b: 'pinned' }, kpi: 'x_th', expected: Math.sqrt(KB * 298.15 / (PI ** 4 * NB.E * nbI / (2 * NB.L ** 3))), tol: 1e-5, kind: 'analytical',
    reference: 'Equipartition theorem applied to the fundamental mode (modal mass m/2 for a sine mode shape)',
  },
  {
    id: 's15-steel-water-transmission', suite: 15, title: 'Steel–water interface: power transmission 4Z₁Z₂/(Z₁ + Z₂)² ≈ 12 %', mode: 'nanobeam',
    inputs: { mode: 'nanobeam', matA: 'steel', matB: 'water' }, kpi: 'T_ac', expected: 4 * zS * zW / (zS + zW) ** 2, tol: 1e-6, kind: 'analytical',
    reference: 'Kinsler, Frey, Coppens & Sanders, Fundamentals of Acoustics — normal-incidence transmission between two media (steel 7850 kg/m³, 5900 m/s; water 998 kg/m³, 1482 m/s)',
  },
  // ---- published experiments ----
  {
    id: 's15-exp-copper-polycrystal-e', suite: 15, title: 'Copper: measured polycrystal Young\'s modulus against the Hill average of the crystal constants', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Cu' }, kpi: 'E', expected: 123.5e9, tol: 0.05, kind: 'experiment',
    measured: { value: 123.5, uncertainty: 0.7, unit: 'GPa' },
    source: { authors: 'Ledbetter, H. M. & Naimon, E. R.', title: 'Elastic properties of metals and alloys. II. Copper', journal: 'J. Phys. Chem. Ref. Data', year: 1974, doi: '10.1063/1.3253150', url: 'https://srd.nist.gov/jpcrdreprint/1.3253150.pdf' },
    conditions: 'Annealed polycrystalline copper at room temperature; mean of the published static and dynamic measurements compiled in the source (outliers beyond 1.5 σ excluded), uncertainty is the standard error of the mean',
    match: ['xtal'],
    reference: 'Ledbetter & Naimon, J. Phys. Chem. Ref. Data 3, 897 (1974)',
    note: 'Measured polycrystal value from table 10 (E = 1.235 ± 0.007 × 10¹¹ N/m²). The suite computes the modulus from its built-in single-crystal constants of copper (168.4, 121.4, 75.4 GPa; the source recommends 169.1, 122.2, 75.42 GPa). The source reports that the arithmetic Voigt–Reuss–Hill average overestimates the measured polycrystal modulus of this strongly anisotropic metal by about 3 %, which the tolerance reflects.'
  },
  {
    id: 's15-exp-copper-polycrystal-g', suite: 15, title: 'Copper: measured polycrystal shear modulus against the Hill average of the crystal constants', mode: 'christoffel',
    inputs: { mode: 'christoffel', xtal: 'Cu' }, kpi: 'G', expected: 45.4e9, tol: 0.07, kind: 'experiment',
    measured: { value: 45.4, uncertainty: 1.2, unit: 'GPa' },
    source: { authors: 'Ledbetter, H. M. & Naimon, E. R.', title: 'Elastic properties of metals and alloys. II. Copper', journal: 'J. Phys. Chem. Ref. Data', year: 1974, doi: '10.1063/1.3253150', url: 'https://srd.nist.gov/jpcrdreprint/1.3253150.pdf' },
    conditions: 'Annealed polycrystalline copper at room temperature; mean of the published torsional and ultrasonic measurements compiled in the source, uncertainty is the standard error of the mean',
    match: ['xtal'],
    reference: 'Ledbetter & Naimon, J. Phys. Chem. Ref. Data 3, 897 (1974)',
    note: 'Measured polycrystal value from table 10 (G = 0.454 ± 0.012 × 10¹¹ N/m²). The suite uses its built-in single-crystal constants of copper. The source gives 47.3 GPa for the arithmetic Voigt–Reuss–Hill average and notes that it lies above the observed shear modulus; the tolerance is the 2.6 % standard error of the measurement plus that 4 % averaging bias.'
  },
];
