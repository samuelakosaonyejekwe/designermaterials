// Suite 12 — reference benchmarks (electromagnetic, photonic and optical engineering).
// Every expected value is either a closed-form result evaluated here from standard constants or a widely tabulated
// handbook value; none is taken from the suite's own output.
const PI = Math.PI, C0 = 299792458, EPS0 = 8.8541878128e-12, MU0 = 4e-7 * PI, DEG = PI / 180;
const sq = v => v * v;

// --- thin films and quarter-wave stacks
const fresnelS = (n, thDeg) => { const c = Math.cos(thDeg * DEG), s = Math.sin(thDeg * DEG), q = Math.sqrt(n * n - s * s); return sq((c - q) / (c + q)); };
const qwR = (n0, nH, nL, ns, N) => { const Y = (nH / nL) ** (2 * N) * ns; return sq((n0 - Y) / (n0 + Y)); };
const qwBw = (nH, nL) => (4 / PI) * Math.asin((nH - nL) / (nH + nL));
const N_BK7_D = 1.5168, N_SILICA_D = 1.4584, N_RUTILE_O = 2.613, N_BK7_NA = 1.5167; // handbook indices (587.6 nm / 589.3 nm)
const yMirror = (N_RUTILE_O / N_SILICA_D) ** 16 * N_BK7_NA, tMirror = 4 * yMirror / sq(1 + yMirror);

// --- Pendry thin-wire plasma (r = 1 µm, a = 5 mm, aluminium) and Drude form
const wireA = 5e-3, wireR = 1e-6, wireSig = 3.7e7;
const wpWire = Math.sqrt(2 * PI * C0 * C0 / (wireA * wireA * Math.log(wireA / wireR)));
const gWire = EPS0 * wireA * wireA * wpWire * wpWire / (PI * wireR * wireR * wireSig);
const wEnz = Math.sqrt(wpWire * wpWire - gWire * gWire);
const epsWire = w => 1 - wpWire * wpWire / (w * w + gWire * gWire);

// --- Pendry split-ring medium (r = 2 mm, c = 1 mm, d = 0.1 mm, l = 2 mm, a = 5 mm)
const sr = { r: 2e-3, c: 1e-3, d: 1e-4, l: 2e-3, a: 5e-3 };
const f0Srr = Math.sqrt(3 * sr.l * C0 * C0 / (PI * Math.log(2 * sr.c / sr.d) * sr.r ** 3)) / (2 * PI);
const Fsrr = PI * sr.r * sr.r / (sr.a * sr.a), fmpSrr = f0Srr / Math.sqrt(1 - Fsrr);
const muLossless = (x, F) => 1 - F * x * x / (x * x - 1); // x = f/f0

// --- negative-index composite (lossless limits)
const nimEps = (f, fep) => 1 - sq(fep / f), nimMu = (f, f0, F) => 1 - F * f * f / (f * f - f0 * f0);

// --- small-particle (Rayleigh) limits of Mie theory
const xOf = (dnm, lam, nmed) => PI * dnm * nmed / lam;
const qScaRayleigh = (x, m) => (8 / 3) * x ** 4 * sq((m * m - 1) / (m * m + 2));
const imClausius = (n, k) => { const a = n * n - k * k, b = 2 * n * k; return (b * (a + 2) - (a - 1) * b) / (sq(a + 2) + b * b); }; // Im[(m²−1)/(m²+2)]

// --- mixing rules
const hsLower = (ei, eh, f) => eh + f / (1 / (ei - eh) + (1 - f) / (3 * eh));
const bruggemanHalf = (ei, eh) => ((ei + eh) + Math.sqrt(sq(ei + eh) + 32 * ei * eh)) / 8;
// layered medium: metal ε = (0.05 + 3i)², dielectric ε = 1, metal fraction 0.2
const emRe = 0.05 * 0.05 - 9, emIm = 2 * 0.05 * 3, fm = 0.2;
const epsParRe = fm * emRe + (1 - fm);
const epsPerpRe = (() => { const dr = (1 - fm) * emRe + fm, di = (1 - fm) * emIm; return (emRe * dr + emIm * di) / (dr * dr + di * di); })();

const BRAGG = { mode: 'multilayer', stack: 'bragg', matH: 'custom', n_user: 2, k_user: 0, matL: 'air', substrate: 'air', Nper: 3, lam0: 550, lam_min: 350, lam_max: 900, theta: 0, n_inc: 1, dH: 0, dL: 0 };
const BARE = { mode: 'multilayer', stack: 'film', matH: 'air', substrate: 'custom', n_user: 1.5, k_user: 0, n_inc: 1, lam0: 550, dH: 0 };
const WIRE = { mode: 'neg_eps', medium: 'wire', wire_r: 1, wire_a: 5, metal_sigma: 3.7e7 };
const SRR = { mode: 'neg_mu', srr_r: 2, srr_a: 5, srr_c: 1, srr_d: 0.1, srr_l: 2, srr_t: 35, metal_sigma: 5.8e7 };
const MATCHED = { mode: 'nim', fep: 10 * Math.SQRT2, f0m: 10 * Math.sqrt(0.75), Fm: 0.5, f_probe: 10, ge: 0.0001, gm: 0.0001, d_slab: 10, theta: 30 };
const NIM_LOWLOSS = { mode: 'nim', fep: 12.8, f0m: 10, Fm: 0.3, f_probe: 11, ge: 0.0001, gm: 0.0001, d_slab: 10, theta: 0 };

export default [
  // ------------------------------------------------------------------ multilayers
  { id: 's12-fresnel-normal-glass', suite: 12, title: 'Fresnel reflectance of an air–glass interface at normal incidence', mode: 'multilayer',
    inputs: { ...BARE }, kpi: 'R0', expected: sq((1.5 - 1) / (1.5 + 1)), tol: 1e-6, kind: 'analytical',
    reference: 'Born & Wolf, Principles of Optics — Fresnel formulae, R = ((n − 1)/(n + 1))²',
    note: 'Substrate n = 1.5 with an air "film" so that the transfer-matrix solution reduces to a single interface.' },
  { id: 's12-fresnel-s-45deg', suite: 12, title: 'Fresnel s-polarised reflectance at 45° on n = 1.5', mode: 'multilayer',
    inputs: { ...BARE, theta: 45, pol: 's' }, kpi: 'R0', expected: fresnelS(1.5, 45), tol: 1e-6, kind: 'analytical',
    reference: 'Born & Wolf, Principles of Optics — Fresnel formulae for oblique incidence' },
  { id: 's12-brewster-angle', suite: 12, title: 'Brewster angle: p-polarised reflectance vanishes at atan(n)', mode: 'multilayer',
    inputs: { ...BARE, theta: Math.atan(1.5) / DEG, pol: 'p' }, kpi: 'R0', expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Hecht, Optics — Brewster angle tan θ_B = n₂/n₁' },
  { id: 's12-lossless-energy-balance', suite: 12, title: 'Lossless stack absorbs nothing (R + T = 1)', mode: 'multilayer',
    inputs: { ...BRAGG }, kpi: 'A0', expected: 0, abs: 1e-9, kind: 'conservation',
    reference: 'Energy conservation for non-absorbing dielectric multilayers (Macleod, Thin-Film Optical Filters)' },
  { id: 's12-quarter-wave-stack-reflectance', suite: 12, title: 'Quarter-wave stack: peak reflectance of three n = 2 / n = 1 pairs', mode: 'multilayer',
    inputs: { ...BRAGG }, kpi: 'R0', expected: qwR(1, 2, 1, 1, 3), tol: 1e-6, kind: 'analytical',
    reference: 'Macleod, Thin-Film Optical Filters — quarter-wave stack admittance Y = (n_H/n_L)^(2N) n_s',
    note: 'Non-dispersive custom index so that the closed form is exact: R = (63/65)².' },
  { id: 's12-quarter-wave-stopband-width', suite: 12, title: 'Quarter-wave stack: fractional stop-band width', mode: 'multilayer',
    inputs: { ...BRAGG }, kpi: 'bw', expected: qwBw(2, 1), tol: 1e-4, kind: 'analytical',
    reference: 'Yeh, Optical Waves in Layered Media — Δω/ω₀ = (4/π) asin((n_H − n_L)/(n_H + n_L))' },
  { id: 's12-quarter-wave-stopband-edge', suite: 12, title: 'Quarter-wave stack: short-wavelength band edge', mode: 'multilayer',
    inputs: { ...BRAGG }, kpi: 'gap_lo', expected: 550 / (1 + qwBw(2, 1) / 2), tol: 1e-4, kind: 'analytical',
    reference: 'Yeh, Optical Waves in Layered Media — band edges at ω₀(1 ± Δω/2ω₀) for a quarter-wave stack' },
  { id: 's12-ar-coating-bk7', suite: 12, title: 'Quarter-wave n = 1.38 antireflection layer on BK7 at 587.6 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'custom', n_user: 1.38, k_user: 0, substrate: 'glass', lam0: 587.56, dH: 0, theta: 0, n_inc: 1 }, kpi: 'R0',
    expected: sq((N_BK7_D - 1.38 * 1.38) / (N_BK7_D + 1.38 * 1.38)), tol: 0.003, kind: 'analytical',
    reference: 'Hecht, Optics — single quarter-wave layer R = ((n₀n_s − n₁²)/(n₀n_s + n₁²))², with the catalogue index n_d = 1.5168 of BK7',
    note: 'Checks the BK7 Sellmeier dispersion and the film solution together.' },
  { id: 's12-index-bk7-d-line', suite: 12, title: 'Refractive index of BK7 at the helium d line (587.56 nm)', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'glass', lam0: 587.56 }, pub: 'prop.n_refr', expected: N_BK7_D, tol: 2e-4, kind: 'reference-data',
    reference: 'Schott optical glass catalogue — N-BK7, n_d = 1.5168' },
  { id: 's12-index-fused-silica', suite: 12, title: 'Refractive index of fused silica at 589.3 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'SiO2', lam0: 589.3 }, pub: 'prop.n_refr', expected: N_SILICA_D, tol: 2e-4, kind: 'reference-data',
    reference: 'Malitson dispersion of fused silica; CRC Handbook of Chemistry and Physics, n_D = 1.4584' },
  { id: 's12-index-rutile', suite: 12, title: 'Ordinary refractive index of rutile TiO₂ at 589.3 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'TiO2', lam0: 589.3 }, pub: 'prop.n_refr', expected: N_RUTILE_O, tol: 0.002, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics — rutile, n_o = 2.613 at the sodium D line' },
  { id: 's12-index-silicon-1550', suite: 12, title: 'Refractive index of silicon at 1550 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'Si', lam0: 1550, lam_min: 1000, lam_max: 2000 }, pub: 'prop.n_refr', expected: 3.476, tol: 0.003, kind: 'reference-data',
    reference: 'Palik, Handbook of Optical Constants of Solids — crystalline silicon in the near infrared (n ≈ 3.476 at 1.55 µm, room temperature)' },
  { id: 's12-tio2-sio2-mirror-transmittance', suite: 12, title: 'Eight-pair TiO₂/SiO₂ quarter-wave mirror on BK7: residual transmittance at 589.3 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'bragg', matH: 'TiO2', matL: 'SiO2', substrate: 'glass', Nper: 8, lam0: 589.3, theta: 0, n_inc: 1, dH: 0, dL: 0 }, kpi: 'T0', expected: tMirror, tol: 0.02, kind: 'analytical',
    reference: 'Macleod, Thin-Film Optical Filters — T = 4Y/(1 + Y)², Y = (n_H/n_L)^16 n_s, with handbook indices 2.613 / 1.4584 / 1.5167',
    note: 'The transmittance scales as (n_L/n_H)^16, so an index error of 0.03 % changes it by about 1 %; hence the 2 % tolerance.' },
  { id: 's12-aluminium-mirror', suite: 12, title: 'Normal-incidence reflectance of an opaque aluminium film at 550 nm', mode: 'multilayer',
    inputs: { mode: 'multilayer', stack: 'film', matH: 'Al', dH: 200, substrate: 'air', lam0: 550, theta: 0, n_inc: 1 }, kpi: 'R0', expected: 0.92, tol: 0.02, kind: 'reference-data',
    reference: 'CRC Handbook of Chemistry and Physics / Palik — evaporated aluminium reflects about 92 % in the green',
    note: 'Lorentz–Drude fit of the measured optical constants; tabulated reflectances scatter by about 1 % with film preparation.' },

  // ------------------------------------------------------------------ negative permittivity
  { id: 's12-wire-plasma-pendry-value', suite: 12, title: 'Thin-wire medium: plasma frequency of 1 µm wires on a 5 mm lattice', mode: 'neg_eps',
    inputs: { ...WIRE }, kpi: 'fp', expected: 8.2, tol: 0.01, kind: 'reference-data',
    reference: 'Pendry, Holden, Stewart & Youngs, "Extremely low frequency plasmons in metallic mesostructures" — about 8.2 GHz for this geometry' },
  { id: 's12-wire-plasma-closed-form', suite: 12, title: 'Thin-wire medium: ω_p² = 2πc²/(a² ln(a/r))', mode: 'neg_eps',
    inputs: { ...WIRE }, kpi: 'fp', expected: wpWire / (2 * PI) / 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Pendry et al., thin-wire effective plasma frequency' },
  { id: 's12-wire-damping', suite: 12, title: 'Thin-wire medium: damping rate γ = ε₀a²ω_p²/(πr²σ)', mode: 'neg_eps',
    inputs: { ...WIRE }, kpi: 'gamma', expected: gWire / (2 * PI) / 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Pendry et al., thin-wire medium — γ ≈ 0.1 ω_p for aluminium wires' },
  { id: 's12-drude-enz-frequency', suite: 12, title: 'Drude medium: Re ε changes sign at √(ω_p² − γ²)', mode: 'neg_eps',
    inputs: { ...WIRE }, kpi: 'f_enz', expected: wEnz / (2 * PI) / 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Drude model, ε(ω) = 1 − ω_p²/(ω(ω + iγ)) (Ashcroft & Mermin, Solid State Physics)' },
  { id: 's12-drude-crossover-zero', suite: 12, title: 'Drude medium: Re ε = 0 at the epsilon-near-zero frequency', mode: 'neg_eps',
    inputs: { ...WIRE, f_rel: 1 }, kpi: 'eps_re', expected: 0, abs: 1e-9, kind: 'analytical',
    reference: 'Drude model (Ashcroft & Mermin, Solid State Physics)' },
  { id: 's12-drude-negative-eps', suite: 12, title: 'Drude medium: negative permittivity at 0.6 of the crossover frequency', mode: 'neg_eps',
    inputs: { ...WIRE, f_rel: 0.6 }, pub: 'eff.eps', expected: epsWire(0.6 * wEnz), tol: 1e-5, kind: 'analytical',
    reference: 'Drude model, Re ε = 1 − ω_p²/(ω² + γ²)' },
  { id: 's12-aluminium-volume-plasmon', suite: 12, title: 'Aluminium: bulk plasma (Re ε = 0) frequency', mode: 'neg_eps',
    inputs: { mode: 'neg_eps', medium: 'Al' }, kpi: 'f_enz', expected: 15.0 * 241.79893, tol: 0.03, kind: 'reference-data',
    reference: 'Kittel, Introduction to Solid State Physics — volume plasmon energy of aluminium, about 15 eV (electron energy loss)',
    note: 'KPI in THz (1 eV = 241.799 THz). The suite uses a single free-electron Drude term fitted to infrared data, which places the crossover within 2 % of the measured plasmon energy.' },

  // ------------------------------------------------------------------ negative permeability
  { id: 's12-srr-resonance', suite: 12, title: 'Split-ring medium: magnetic resonance frequency', mode: 'neg_mu',
    inputs: { ...SRR }, kpi: 'f0', expected: f0Srr / 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Pendry, Holden, Robbins & Stewart, "Magnetism from conductors and enhanced nonlinear phenomena" — ω₀² = 3lc²/(π ln(2c/d) r³)',
    note: 'Geometry of the original paper (r = 2 mm, c = 1 mm, d = 0.1 mm, l = 2 mm), which resonates at about 13.5 GHz.' },
  { id: 's12-srr-magnetic-plasma', suite: 12, title: 'Split-ring medium: magnetic plasma frequency ω₀/√(1 − F)', mode: 'neg_mu',
    inputs: { ...SRR }, kpi: 'fmp', expected: fmpSrr / 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Pendry et al., split-ring effective permeability with F = πr²/a²' },
  { id: 's12-srr-negative-band-lower', suite: 12, title: 'Split-ring medium: lower edge of the negative-µ band', mode: 'neg_mu',
    inputs: { ...SRR }, kpi: 'neg_lo', expected: f0Srr / 1e9, tol: 0.005, kind: 'analytical',
    reference: 'Pendry et al. — Re µ < 0 for ω₀ < ω < ω₀/√(1 − F)',
    note: 'The band edges are read from a 400-point frequency sweep, so they are resolved to about 0.3 % of ω₀.' },
  { id: 's12-srr-negative-band-upper', suite: 12, title: 'Split-ring medium: upper edge of the negative-µ band', mode: 'neg_mu',
    inputs: { ...SRR }, kpi: 'neg_hi', expected: fmpSrr / 1e9, tol: 0.005, kind: 'analytical',
    reference: 'Pendry et al. — Re µ < 0 for ω₀ < ω < ω₀/√(1 − F)',
    note: 'Edge read from a 400-point frequency sweep (resolution about 0.3 % of ω₀).' },
  { id: 's12-srr-mu-above-resonance', suite: 12, title: 'Split-ring medium: Re µ at 1.2 ω₀ (negative permeability)', mode: 'neg_mu',
    inputs: { ...SRR, f_rel: 1.2 }, kpi: 'mu_re', expected: muLossless(1.2, Fsrr), tol: 1e-3, kind: 'analytical',
    reference: 'Pendry et al. — µ = 1 − Fω²/(ω² − ω₀² + iΓω), evaluated in the low-loss limit',
    note: 'Copper rings give Q ≈ 1750, so the loss correction to Re µ at 1.2 ω₀ is of order 10⁻⁵.' },

  // ------------------------------------------------------------------ negative refractive index
  { id: 's12-nim-matched-index', suite: 12, title: 'Double-negative medium with ε = µ = −1: refractive index', mode: 'nim',
    inputs: { ...MATCHED }, kpi: 'n_re', expected: -1, tol: 1e-3, kind: 'analytical',
    reference: 'Veselago, "The electrodynamics of substances with simultaneously negative values of ε and µ" — n = −√(εµ)',
    note: 'Drude ε and Lorentz µ tuned so that both equal −1 at 10 GHz (f_ep = 10√2 GHz, f₀ = 10√0.75 GHz, F = 0.5) with negligible loss.' },
  { id: 's12-nim-matched-impedance', suite: 12, title: 'Double-negative medium with ε = µ = −1: wave impedance', mode: 'nim',
    inputs: { ...MATCHED }, kpi: 'z_re', expected: 1, tol: 1e-3, kind: 'analytical',
    reference: 'Veselago medium, z = √(µ/ε) = 1 (impedance matched to free space)' },
  { id: 's12-nim-negative-refraction', suite: 12, title: 'Negative refraction: Snell angle for n = −1 at 30° incidence', mode: 'nim',
    inputs: { ...MATCHED }, kpi: 'theta_t', expected: -30, tol: 1e-3, kind: 'analytical',
    reference: 'Veselago; Pendry, "Negative refraction makes a perfect lens" — sin θ_t = sin θ_i / n' },
  { id: 's12-nim-matched-slab-reflection', suite: 12, title: 'Matched n = −1 slab: no reflection', mode: 'nim',
    inputs: { ...MATCHED }, kpi: 'S11', expected: 0, abs: 1e-6, kind: 'analytical',
    reference: 'Fresnel–Airy slab formula with z = 1 (Pendry, "Negative refraction makes a perfect lens")' },
  { id: 's12-nim-index-closed-form', suite: 12, title: 'Wire + split-ring composite: n = −√(εµ) at 11 GHz', mode: 'nim',
    inputs: { ...NIM_LOWLOSS }, kpi: 'n_re', expected: -Math.sqrt(nimEps(11, 12.8) * nimMu(11, 10, 0.3)), tol: 1e-3, kind: 'analytical',
    reference: 'Smith, Padilla, Vier, Nemat-Nasser & Schultz, "Composite medium with simultaneously negative permeability and permittivity" — Drude ε and Lorentz µ' },
  { id: 's12-nim-band-lower', suite: 12, title: 'Double-negative band: lower edge at the magnetic resonance', mode: 'nim',
    inputs: { ...NIM_LOWLOSS }, kpi: 'dng_lo', expected: 10, tol: 0.005, kind: 'analytical',
    reference: 'Smith et al. — ε < 0 and µ < 0 for f₀ < f < min(f_ep, f₀/√(1 − F))',
    note: 'Edge read from a 400-point frequency sweep (resolution about 0.4 % here).' },
  { id: 's12-nim-band-upper', suite: 12, title: 'Double-negative band: upper edge at the magnetic plasma frequency', mode: 'nim',
    inputs: { ...NIM_LOWLOSS }, kpi: 'dng_hi', expected: 10 / Math.sqrt(1 - 0.3), tol: 0.005, kind: 'analytical',
    reference: 'Smith et al. — ε < 0 and µ < 0 for f₀ < f < min(f_ep, f₀/√(1 − F))',
    note: 'Edge read from a 400-point frequency sweep (resolution about 0.4 % here).' },
  { id: 's12-nim-retrieval', suite: 12, title: 'S-parameter retrieval returns the index of the constitutive model', mode: 'nim',
    inputs: { mode: 'nim', fep: 12.8, ge: 0.01, f0m: 10, Fm: 0.3, gm: 0.01, d_slab: 10, f_probe: 11 }, kpi: 'ret_err', expected: 0, abs: 1e-6, kind: 'conservation',
    reference: 'Smith, Schultz, Markoš & Soukoulis, "Determination of effective permittivity and permeability of metamaterials from reflection and transmission coefficients"',
    note: 'Round trip model → slab S-parameters → inversion with branch tracking, across a resonance of quality factor 1000.' },

  // ------------------------------------------------------------------ Mie theory and plasmonic particles
  { id: 's12-mie-rayleigh-scattering', suite: 12, title: 'Mie series in the Rayleigh limit: scattering efficiency of a 10 nm dielectric sphere', mode: 'mie',
    inputs: { mode: 'mie', mat_p: 'custom', n_user: 1.5, k_user: 0, d_p: 10, n_med: 1, lam_min: 600, lam_max: 900 }, kpi: 'Qext',
    expected: qScaRayleigh(xOf(10, 600, 1), 1.5), tol: 0.003, kind: 'analytical',
    reference: 'Bohren & Huffman, Absorption and Scattering of Light by Small Particles — Q_sca = (8/3)x⁴|(m² − 1)/(m² + 2)|²',
    note: 'Lossless sphere, so extinction equals scattering; the maximum of the 600–900 nm window is at 600 nm (x = 0.052).' },
  { id: 's12-mie-rayleigh-absorption', suite: 12, title: 'Mie series in the Rayleigh limit: absorption cross-section of a 5 nm absorbing sphere', mode: 'mie',
    inputs: { mode: 'mie', mat_p: 'custom', n_user: 1.5, k_user: 0.1, d_p: 5, n_med: 1, lam_min: 600, lam_max: 900 }, kpi: 'sig_abs',
    expected: 4 * xOf(5, 600, 1) * imClausius(1.5, 0.1) * PI * sq(2.5e-9), tol: 0.003, kind: 'analytical',
    reference: 'Bohren & Huffman — Q_abs = 4x Im[(m² − 1)/(m² + 2)]' },
  { id: 's12-gold-colloid-40nm', suite: 12, title: 'Gold nanospheres (40 nm) in water: extinction peak', mode: 'mie',
    inputs: { mode: 'mie', mat_p: 'Au', d_p: 40, n_med: 1.33, lam_min: 400, lam_max: 800 }, kpi: 'lam_peak', expected: 525, tol: 0.015, kind: 'reference-data',
    reference: 'Widely tabulated plasmon band of citrate gold colloids (about 520–530 nm for 20–40 nm spheres); Kreibig & Vollmer, Optical Properties of Metal Clusters',
    note: 'Bulk Lorentz–Drude dielectric function; measured peaks scatter by a few nm with size distribution and ligand shell.' },
  { id: 's12-silver-colloid-20nm', suite: 12, title: 'Silver nanospheres (20 nm) in water: extinction peak', mode: 'mie',
    inputs: { mode: 'mie', mat_p: 'Ag', d_p: 20, n_med: 1.33, lam_min: 300, lam_max: 700 }, kpi: 'lam_peak', expected: 400, tol: 0.04, kind: 'reference-data',
    reference: 'Widely tabulated plasmon band of small silver colloids in water (about 400 nm); Kreibig & Vollmer, Optical Properties of Metal Clusters',
    note: 'The Lorentz–Drude fit for silver is smooth near the interband edge and places the resonance about 10 nm to the red of measurements.' },
  { id: 's12-silver-froehlich-air', suite: 12, title: 'Silver sphere in air: Fröhlich wavelength (Re ε = −2)', mode: 'mie',
    inputs: { mode: 'mie', mat_p: 'Ag', d_p: 10, n_med: 1, lam_min: 300, lam_max: 700 }, kpi: 'lam_F', expected: 355, tol: 0.05, kind: 'reference-data',
    reference: 'Johnson & Christy optical constants of silver (Re ε = −2 near 3.5 eV); Bohren & Huffman, Fröhlich condition ε′ = −2ε_m',
    note: 'The analytic Lorentz–Drude fit used by the suite crosses −2 near 368 nm, 3–4 % to the red of the tabulated constants.' },

  // ------------------------------------------------------------------ effective media
  { id: 's12-maxwell-garnett-hashin-shtrikman', suite: 12, title: 'Maxwell-Garnett mixture equals the Hashin–Shtrikman lower bound', mode: 'effmed',
    inputs: { mode: 'effmed', mat_i: 'custom', n_user: 2.5, k_user: 0, n_host: 1.5, fill: 0.3, lam0: 550 }, kpi: 'eps_mg', expected: hsLower(6.25, 2.25, 0.3), tol: 1e-6, kind: 'analytical',
    reference: 'Hashin & Shtrikman bounds for two-phase media; Milton, The Theory of Composites' },
  { id: 's12-maxwell-garnett-dilute', suite: 12, title: 'Maxwell-Garnett dilute limit ε_h(1 + 3fβ)', mode: 'effmed',
    inputs: { mode: 'effmed', mat_i: 'custom', n_user: 2, k_user: 0, n_host: 1.5, fill: 0.005, lam0: 550 }, kpi: 'eps_mg', expected: 2.25 * (1 + 3 * 0.005 * (4 - 2.25) / (4 + 2 * 2.25)), tol: 1e-4, kind: 'analytical',
    reference: 'Clausius–Mossotti dilute limit (Landau & Lifshitz, Electrodynamics of Continuous Media)',
    note: 'First-order expansion in f; the neglected term is of order f²β² ≈ 10⁻⁶.' },
  { id: 's12-bruggeman-symmetric', suite: 12, title: 'Bruggeman mixture at equal volume fractions', mode: 'effmed',
    inputs: { mode: 'effmed', mat_i: 'custom', n_user: 2.5, k_user: 0, n_host: 1.5, fill: 0.5, lam0: 550 }, kpi: 'eps_bg', expected: bruggemanHalf(6.25, 2.25), tol: 1e-6, kind: 'analytical',
    reference: 'Bruggeman symmetric effective-medium approximation for spheres, closed-form root at f = 1/2' },
  { id: 's12-multilayer-eps-parallel', suite: 12, title: 'Metal–dielectric multilayer: in-plane permittivity (Type II hyperbolic)', mode: 'hyperbolic',
    inputs: { mode: 'hyperbolic', mat_i: 'custom', n_user: 0.05, k_user: 3, mat_d: 'air', fill: 0.2, period: 10, lam0: 550 }, kpi: 'eps_par', expected: epsParRe, tol: 1e-6, kind: 'analytical',
    reference: 'Rytov effective-medium theory of layered media — ε∥ = fε_m + (1 − f)ε_d',
    note: 'Negative in-plane and positive normal permittivity: the Type II hyperbolic regime.' },
  { id: 's12-multilayer-eps-perpendicular', suite: 12, title: 'Metal–dielectric multilayer: normal permittivity', mode: 'hyperbolic',
    inputs: { mode: 'hyperbolic', mat_i: 'custom', n_user: 0.05, k_user: 3, mat_d: 'air', fill: 0.2, period: 10, lam0: 550 }, kpi: 'eps_perp', expected: epsPerpRe, tol: 1e-6, kind: 'analytical',
    reference: 'Rytov effective-medium theory of layered media — 1/ε⊥ = f/ε_m + (1 − f)/ε_d' },

  // ------------------------------------------------------------------ metasurfaces
  { id: 's12-generalised-snell-refraction', suite: 12, title: 'Phase-gradient metasurface: anomalous refraction angle', mode: 'metasurface',
    inputs: { mode: 'metasurface', ms_kind: 'gradient', lam0: 600, ms_p: 300, ms_N: 8, theta: 0, n_inc: 1, n_t: 1.5 }, kpi: 'theta_t', expected: Math.asin(600 / (8 * 300) / 1.5) / DEG, tol: 1e-6, kind: 'analytical',
    reference: 'Yu et al., "Light propagation with phase discontinuities: generalized laws of reflection and refraction" — n_t sin θ_t − n_i sin θ_i = (λ/2π) dφ/dx' },
  { id: 's12-generalised-snell-reflection', suite: 12, title: 'Phase-gradient metasurface: anomalous reflection angle', mode: 'metasurface',
    inputs: { mode: 'metasurface', ms_kind: 'gradient', lam0: 600, ms_p: 300, ms_N: 8, theta: 0, n_inc: 1, n_t: 1.5 }, kpi: 'theta_r', expected: Math.asin(0.25) / DEG, tol: 1e-6, kind: 'analytical',
    reference: 'Yu et al., generalised law of reflection' },
  { id: 's12-four-level-efficiency', suite: 12, title: 'Four-level phase quantisation: first-order efficiency', mode: 'metasurface',
    inputs: { mode: 'metasurface', ms_kind: 'gradient', lam0: 600, ms_p: 300, ms_N: 4 }, kpi: 'eff', expected: sq(Math.sin(PI / 4) / (PI / 4)), tol: 1e-6, kind: 'analytical',
    reference: 'Goodman, Introduction to Fourier Optics — N-level phase grating efficiency sinc²(1/N) (81 % for four levels)' },
  { id: 's12-metalens-numerical-aperture', suite: 12, title: 'Metalens numerical aperture (D = 200 µm, f = 75 µm)', mode: 'metasurface',
    inputs: { mode: 'metasurface', ms_kind: 'lens', lam0: 632.8, lens_D: 200, lens_f: 75, n_t: 1 }, kpi: 'NA', expected: 0.8, tol: 1e-6, kind: 'analytical',
    reference: 'Born & Wolf, Principles of Optics — NA = n sin(atan(D/2f))' },
  { id: 's12-metalens-airy-radius', suite: 12, title: 'Metalens diffraction-limited spot: Airy radius 0.61λ/NA', mode: 'metasurface',
    inputs: { mode: 'metasurface', ms_kind: 'lens', lam0: 632.8, lens_D: 200, lens_f: 75, n_t: 1 }, kpi: 'airy', expected: 0.61 * 632.8 / 0.8, tol: 1e-6, kind: 'analytical',
    reference: 'Born & Wolf, Principles of Optics — Airy pattern of a circular aperture (first zero at 0.61λ/NA)' },
];
