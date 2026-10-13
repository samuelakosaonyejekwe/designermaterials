// Suite 13 — reference benchmarks (magnetic and spintronic engineering).
// Expected values are closed-form results evaluated here from CODATA constants and the stated inputs, or classic
// tabulated results (Stoner–Wohlfarth ensemble, Kittel 28 GHz/T); none is taken from the suite's own output.
const PI = Math.PI, KB = 1.380649e-23, MU0 = 4e-7 * PI, QE = 1.602176634e-19, HBAR = 1.054571817e-34;
const GAMMA = 1.76085963e11, GHZ_PER_T = 28.0249514; // electron gyromagnetic ratio (rad s⁻¹ T⁻¹) and γ/2π
const sq = v => v * v, sphereV = d => PI / 6 * d ** 3;

// --- reference particle: magnetite-like, 16 nm core, 3 nm shell, water at 25 °C
const Ms = 4.8e5, K = 1.35e4, A = 1.2e-11, RHO = 5180, T25 = 298.15, kT = KB * T25, V16 = sphereV(16e-9), TAU0 = 1e-9;
const sigma = K * V16 / kT, tauN = TAU0 * Math.exp(sigma), tauB = 3 * 1e-3 * sphereV(22e-9) / kT, tauE = tauN * tauB / (tauN + tauB);
const phi = 1e-3, chi0 = MU0 * phi * Ms * Ms * V16 / (3 * kT);
const Hac = 100, fac = 5e5, wt = 2 * PI * fac * tauE, Pvol = PI * MU0 * chi0 * fac * Hac * Hac * wt / (1 + wt * wt);
const xiDC = MU0 * Ms * V16 * 4e5 / kT;
const SPM = { mode: 'spm', material: 'custom', Ms, K, A_ex: 12, rho: RHO, T: 25, d_nm: 16, sigma_ln: 0, tau_m: TAU0 * Math.exp(25), coat: 3, visc: 1, phi_v: 0.1, H_ac: 0.1, f_ac: 500, H_dc: 400 };

// --- Stoner–Wohlfarth
const Hk = 2 * K / (MU0 * Ms), astroid = psiDeg => { const p = psiDeg * PI / 180; return (Math.cos(p) ** (2 / 3) + Math.sin(p) ** (2 / 3)) ** -1.5; };
const SW = { mode: 'sw', material: 'custom', Ms, K, T: -268, d_nm: 100, hmax: 1.5, psi: 30 };

// --- macrospin LLG: damped reversal in an antiparallel field B > B_K from the thermal tilt θ₀ = √(kT/2KV)
const Bk = 2 * K / Ms, T37 = 310.15, th0 = Math.sqrt(KB * T37 / (2 * K * V16));
const reversalIntegral = (b, a, u0) => -Math.log(1 - u0) / (2 * (b - a)) + Math.log(1 + u0) / (2 * (b + a)) + a / (b * b - a * a) * Math.log((b - a * u0) / b); // ∫₀^u0 du/((1−u²)(b−au))
const tReversal = alpha => (1 + alpha * alpha) / (alpha * GAMMA) * reversalIntegral(0.2, Bk, Math.cos(th0));
const LLG = { mode: 'llg', material: 'custom', Ms, K, d_nm: 16, T: 37, B_pulse: 200, th_H: 0, J_stt: 0 };

// --- ferromagnetic resonance (custom ferrite, anisotropy at its lower bound so that H_K = 0.33 kA/m ≪ H)
const Hb = 3e5, HkF = 2 * 100 / (MU0 * Ms), f0F = GHZ_PER_T * MU0 * (Hb + HkF), fmF = GHZ_PER_T * MU0 * Ms;
const FMR = { mode: 'fmr', material: 'custom', Ms, K: 100, H_bias: 300, alpha: 0.0001, f_fill: 1, f_rel: 1.1 };

// --- spintronics (40 nm junction, 1.3 nm free layer)
const area = PI * sq(40e-9) / 4, tf = 1.3e-9, MsF = 1e6, Keff = 1.3e-3 / tf - 0.5 * MU0 * MsF * MsF, Delta = Keff * area * tf / kT;
const SPIN = { mode: 'spin', T: 25, P1: 0.6, P2: 0.6, mtj_d: 40, t_free: 1.3, K_i: 1.3, Ms_free: 1e6, alpha: 0.01, RA: 10 };

export default [
  // ------------------------------------------------------------------ superparamagnetism and hyperthermia
  { id: 's13-blocking-temperature', suite: 13, title: 'Blocking temperature T_B = KV/(25 k_B)', mode: 'spm',
    inputs: { ...SPM }, kpi: 'TB', expected: K * V16 / (25 * KB), tol: 1e-4, kind: 'analytical',
    reference: 'Bean & Livingston, "Superparamagnetism"; Cullity & Graham, Introduction to Magnetic Materials — KV = 25 k_B T_B',
    note: 'Measurement time set to τ₀e²⁵ ≈ 72 s (τ₀ = 1 ns) so that ln(τ_m/τ₀) = 25, the classic laboratory criterion.' },
  { id: 's13-superparamagnetic-diameter', suite: 13, title: 'Superparamagnetic critical diameter (25 k_B T = KV)', mode: 'spm',
    inputs: { ...SPM }, kpi: 'd_sp', expected: Math.cbrt(6 * 25 * kT / (PI * K)) * 1e9, tol: 1e-4, kind: 'analytical',
    reference: 'Cullity & Graham, Introduction to Magnetic Materials — superparamagnetic limit' },
  { id: 's13-neel-relaxation-time', suite: 13, title: 'Néel relaxation time τ₀ exp(KV/k_BT)', mode: 'spm',
    inputs: { ...SPM }, kpi: 'tauN', expected: tauN, tol: 1e-3, kind: 'analytical',
    reference: 'Néel–Arrhenius law (Néel 1949; Coey, Magnetism and Magnetic Materials)',
    note: 'Attempt time τ₀ = 1 ns, as documented by the suite; KV/k_BT = 7.03 here.' },
  { id: 's13-brown-relaxation-time', suite: 13, title: 'Brownian relaxation time 3ηV_h/k_BT', mode: 'spm',
    inputs: { ...SPM }, kpi: 'tauB', expected: tauB, tol: 1e-4, kind: 'analytical',
    reference: 'Debye rotational diffusion; Rosensweig, "Heating magnetic fluid with alternating magnetic field"',
    note: 'Hydrodynamic diameter 22 nm (16 nm core + 3 nm shell) in a 1 mPa·s carrier.' },
  { id: 's13-langevin-initial-susceptibility', suite: 13, title: 'Langevin low-field limit: χ₀ = µ₀φM_s²V/(3k_BT)', mode: 'spm',
    inputs: { ...SPM }, kpi: 'chi', expected: chi0, tol: 1e-4, kind: 'analytical',
    reference: 'Langevin paramagnetism of monodisperse single-domain particles (Rosensweig, Ferrohydrodynamics)' },
  { id: 's13-langevin-high-field', suite: 13, title: 'Langevin high-field limit: M/M_s = 1 − 1/ξ', mode: 'spm',
    inputs: { ...SPM }, kpi: 'sat', expected: 1 - 1 / xiDC, tol: 1e-5, kind: 'analytical',
    reference: 'Langevin function L(ξ) = coth ξ − 1/ξ for ξ ≫ 1 (Cullity & Graham)',
    note: 'ξ = µ₀M_sVH/k_BT ≈ 126 at 400 kA/m, where coth ξ = 1 to machine precision.' },
  { id: 's13-rosensweig-sar', suite: 13, title: 'Linear-response heating: specific absorption rate of a monodisperse ferrofluid', mode: 'spm',
    inputs: { ...SPM }, kpi: 'SAR', expected: Pvol / (RHO * phi), tol: 2e-3, kind: 'analytical',
    reference: 'Rosensweig, "Heating magnetic fluid with alternating magnetic field" — P = πµ₀χ₀fH₀² ωτ/(1 + ω²τ²)',
    note: 'Field amplitude 0.1 kA/m (ξ = 0.03) so that the field-dependent correction to χ₀ is below 0.01 %.' },
  { id: 's13-domain-wall-width-iron', suite: 13, title: 'Bloch wall width π√(A/K) for iron-like constants', mode: 'spm',
    inputs: { ...SPM, Ms: 1.71e6, K: 4.8e4, A_ex: 21 }, kpi: 'dw', expected: PI * Math.sqrt(21e-12 / 4.8e4) * 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Landau–Lifshitz 180° wall, δ = π√(A/K) (Coey, Magnetism and Magnetic Materials)' },
  { id: 's13-exchange-length-iron', suite: 13, title: 'Exchange length √(2A/µ₀M_s²) for iron-like constants', mode: 'spm',
    inputs: { ...SPM, Ms: 1.71e6, K: 4.8e4, A_ex: 21 }, kpi: 'lex', expected: Math.sqrt(2 * 21e-12 / (MU0 * sq(1.71e6))) * 1e9, tol: 1e-5, kind: 'analytical',
    reference: 'Micromagnetic exchange length (Hubert & Schäfer, Magnetic Domains); about 3.4 nm for iron' },

  // ------------------------------------------------------------------ Stoner–Wohlfarth
  { id: 's13-anisotropy-field', suite: 13, title: 'Anisotropy field H_K = 2K/(µ₀M_s)', mode: 'sw',
    inputs: { ...SW }, kpi: 'Hk', expected: Hk, tol: 1e-6, kind: 'analytical',
    reference: 'Stoner & Wohlfarth, "A mechanism of magnetic hysteresis in heterogeneous alloys"' },
  { id: 's13-sw-aligned-coercivity', suite: 13, title: 'Stoner–Wohlfarth particle with the field along the easy axis: H_c = H_K', mode: 'sw',
    inputs: { ...SW, psi: 0 }, kpi: 'Hc1', expected: Hk, tol: 0.005, kind: 'analytical',
    reference: 'Stoner & Wohlfarth — rectangular loop for ψ = 0',
    note: 'The loop is traced with a field step of 0.015 H_K, which bounds the resolution of the switching field.' },
  { id: 's13-sw-aligned-remanence', suite: 13, title: 'Stoner–Wohlfarth particle with the field along the easy axis: M_r/M_s = 1', mode: 'sw',
    inputs: { ...SW, psi: 0 }, kpi: 'Mr1', expected: 1, tol: 1e-6, kind: 'analytical',
    reference: 'Stoner & Wohlfarth — rectangular loop for ψ = 0' },
  { id: 's13-sw-remanence-30deg', suite: 13, title: 'Stoner–Wohlfarth particle at ψ = 30°: remanence cos ψ', mode: 'sw',
    inputs: { ...SW, psi: 30 }, kpi: 'Mr1', expected: Math.cos(PI / 6), tol: 1e-4, kind: 'analytical',
    reference: 'Stoner & Wohlfarth — at zero field the moment lies on the easy axis' },
  { id: 's13-sw-astroid-30deg', suite: 13, title: 'Stoner–Wohlfarth switching field at ψ = 30° (astroid)', mode: 'sw',
    inputs: { ...SW, psi: 30 }, kpi: 'Hsw', expected: astroid(30) * Hk, tol: 0.01, kind: 'analytical',
    reference: 'Stoner–Wohlfarth astroid, h_sw = (cos^(2/3)ψ + sin^(2/3)ψ)^(−3/2)',
    note: 'Numerical switching field located on a field grid of step 0.015 H_K.' },
  { id: 's13-sw-random-remanence', suite: 13, title: 'Randomly oriented Stoner–Wohlfarth ensemble: M_r/M_s = 0.5', mode: 'sw',
    inputs: { ...SW }, kpi: 'Mr', expected: 0.5, tol: 1e-3, kind: 'analytical',
    reference: 'Stoner & Wohlfarth — isotropic ensemble of uniaxial particles' },
  { id: 's13-sw-random-coercivity', suite: 13, title: 'Randomly oriented Stoner–Wohlfarth ensemble: H_c = 0.479 H_K', mode: 'sw',
    inputs: { ...SW }, kpi: 'Hc', expected: 0.479 * Hk, tol: 0.01, kind: 'numerical-reference',
    reference: 'Stoner & Wohlfarth — coercivity of the isotropic ensemble, 0.479 H_K',
    note: 'The suite averages 24 orientations on a field grid of step 0.015 H_K; the quadrature and grid limit agreement to about 0.6 %.' },
  { id: 's13-energy-product-bound', suite: 13, title: 'Energy-product bound µ₀M_r²/4 of the isotropic ensemble', mode: 'sw',
    inputs: { ...SW }, kpi: 'BHmax', expected: MU0 * sq(Ms / 2) / 4, tol: 1e-3, kind: 'analytical',
    reference: 'Maximum energy product of an ideal square-loop magnet, (BH)max = µ₀M_r²/4 (Coey, Magnetism and Magnetic Materials)' },

  // ------------------------------------------------------------------ macrospin LLG
  { id: 's13-llg-anisotropy-resonance', suite: 13, title: 'Zero-field precession frequency γB_K/2π', mode: 'llg',
    inputs: { ...LLG, alpha: 0.1, t_pulse: 5, t_end: 5 }, kpi: 'f_res', expected: GHZ_PER_T * Bk, tol: 1e-5, kind: 'analytical',
    reference: 'Kittel, "On the theory of ferromagnetic resonance absorption" — uniaxial anisotropy field' },
  { id: 's13-llg-damping-time', suite: 13, title: 'Gilbert damping time 1/(αγB_K)', mode: 'llg',
    inputs: { ...LLG, alpha: 0.1, t_pulse: 5, t_end: 5 }, kpi: 'tau', expected: 1e12 / (0.1 * GAMMA * Bk), tol: 1e-5, kind: 'analytical',
    reference: 'Landau–Lifshitz–Gilbert equation, small-angle relaxation about the easy axis' },
  { id: 's13-llg-damped-reversal-time', suite: 13, title: 'LLG integration: damped reversal time in an antiparallel field (α = 0.1)', mode: 'llg',
    inputs: { ...LLG, alpha: 0.1, t_pulse: 5, t_end: 5 }, kpi: 't_sw', expected: tReversal(0.1) * 1e12, tol: 0.01, kind: 'analytical',
    reference: 'Kikuchi, "On the minimum of magnetization reversal time" — axially symmetric LLG, dθ/dt = αγ(B − B_K cos θ) sin θ/(1 + α²), integrated in closed form',
    note: 'Reversal from the thermal tilt θ₀ = √(k_BT/2KV) = 15.6° to θ = 90° in a 200 mT field; the reported time is the first integration step past m_z = 0.' },
  { id: 's13-llg-damped-reversal-low-damping', suite: 13, title: 'LLG integration: damped reversal time for α = 0.01 (about 50 precession periods)', mode: 'llg',
    inputs: { ...LLG, alpha: 0.01, t_pulse: 20, t_end: 20 }, kpi: 't_sw', expected: tReversal(0.01) * 1e12, tol: 0.01, kind: 'analytical',
    reference: 'Kikuchi — closed-form reversal time of the axially symmetric LLG equation',
    note: 'Tests the adaptive integrator over a long precessional trajectory; the time scales as (1 + α²)/α.' },
  { id: 's13-llg-static-switching-field', suite: 13, title: 'Static switching field for a pulse 30° from the reversed easy axis', mode: 'llg',
    inputs: { ...LLG, alpha: 0.1, th_H: 30, t_pulse: 5, t_end: 5 }, kpi: 'Bsw', expected: astroid(30) * Bk * 1e3, tol: 1e-5, kind: 'analytical',
    reference: 'Stoner–Wohlfarth astroid' },
  { id: 's13-stt-critical-current', suite: 13, title: 'Spin-transfer-torque threshold J_c0 = 2eαM_s t B_K/(ħP)', mode: 'llg',
    inputs: { ...LLG, alpha: 0.01, P_spin: 0.6, t_free: 1.3, t_pulse: 2, t_end: 5 }, kpi: 'Jc0', expected: 2 * QE * 0.01 * Ms * 1.3e-9 * Bk / (HBAR * 0.6), tol: 1e-4, kind: 'analytical',
    reference: 'Slonczewski, "Current-driven excitation of magnetic multilayers"; Sun, macrospin instability threshold for a uniaxial free layer' },

  // ------------------------------------------------------------------ FMR and negative permeability
  { id: 's13-kittel-sphere', suite: 13, title: 'Kittel resonance of a sphere: 28.02 GHz per tesla', mode: 'fmr',
    inputs: { ...FMR, shape: 'sphere' }, kpi: 'fK', expected: GHZ_PER_T * MU0 * Hb, tol: 0.003, kind: 'analytical',
    reference: 'Kittel, "On the theory of ferromagnetic resonance absorption" — f = γµ₀H/2π for a sphere (γ/2π = 28.025 GHz/T)',
    note: 'The suite adds the anisotropy field along the bias; at the lowest allowed K it is 0.11 % of the 300 kA/m bias.' },
  { id: 's13-kittel-film-in-plane', suite: 13, title: 'Kittel resonance of an in-plane magnetised film: γµ₀√(H(H + M_s))/2π', mode: 'fmr',
    inputs: { ...FMR, shape: 'film_ip' }, kpi: 'fK', expected: GHZ_PER_T * MU0 * Math.sqrt(Hb * (Hb + Ms)), tol: 0.003, kind: 'analytical',
    reference: 'Kittel formula for a tangentially magnetised thin film',
    note: 'Residual anisotropy field of 0.33 kA/m shifts the result by about 0.1 %.' },
  { id: 's13-kittel-film-out-of-plane', suite: 13, title: 'Kittel resonance of a perpendicularly magnetised film: γµ₀(H − M_s)/2π', mode: 'fmr',
    inputs: { ...FMR, shape: 'film_oop', H_bias: 800 }, kpi: 'fK', expected: GHZ_PER_T * MU0 * (8e5 - Ms), tol: 0.003, kind: 'analytical',
    reference: 'Kittel formula for a normally magnetised thin film',
    note: 'Residual anisotropy field of 0.33 kA/m shifts the result by about 0.1 %.' },
  { id: 's13-polder-negative-band-lower', suite: 13, title: 'Gyrotropic ferrite: lower edge of the negative-µ_eff band, √(f₀(f₀ + f_m))', mode: 'fmr',
    inputs: { ...FMR, shape: 'film_ip' }, kpi: 'neg_lo', expected: Math.sqrt(f0F * (f0F + fmF)), tol: 0.01, kind: 'analytical',
    reference: 'Polder permeability tensor; Pozar, Microwave Engineering — µ_eff = (µ² − κ²)/µ for propagation transverse to the bias',
    note: 'Band edges are read from a 400-point frequency sweep, which resolves them to about 0.3 % here.' },
  { id: 's13-polder-negative-band-upper', suite: 13, title: 'Gyrotropic ferrite: upper edge of the negative-µ_eff band, f₀ + f_m', mode: 'fmr',
    inputs: { ...FMR, shape: 'film_ip' }, kpi: 'neg_hi', expected: f0F + fmF, tol: 0.01, kind: 'analytical',
    reference: 'Polder permeability tensor; Pozar, Microwave Engineering — µ_eff changes sign at ω₀ + ω_m',
    note: 'Band edges are read from a 400-point frequency sweep, which resolves them to about 0.3 % here.' },
  { id: 's13-polder-negative-mu', suite: 13, title: 'Gyrotropic ferrite: negative µ_eff at 1.1 times the resonance frequency', mode: 'fmr',
    inputs: { ...FMR, shape: 'film_ip' }, kpi: 'mu_re',
    expected: (() => { const w2 = 1.21 * f0F * (f0F + fmF); return (sq(f0F + fmF) - w2) / (f0F * (f0F + fmF) - w2); })(), tol: 1e-3, kind: 'analytical',
    reference: 'Polder tensor in the lossless limit, µ_eff = ((ω₀ + ω_m)² − ω²)/(ω₀(ω₀ + ω_m) − ω²)',
    note: 'Gilbert damping 10⁻⁴, so the loss correction to Re µ_eff is negligible.' },

  // ------------------------------------------------------------------ spintronics
  { id: 's13-julliere-tmr', suite: 13, title: 'Julliere tunnel magnetoresistance for P₁ = P₂ = 0.5', mode: 'spin',
    inputs: { ...SPIN, P1: 0.5, P2: 0.5 }, kpi: 'TMR0', expected: 200 / 3, tol: 1e-6, kind: 'analytical',
    reference: 'Julliere, "Tunneling between ferromagnetic films" — TMR = 2P₁P₂/(1 − P₁P₂)' },
  { id: 's13-julliere-tmr-asymmetric', suite: 13, title: 'Julliere tunnel magnetoresistance for P₁ = 0.6, P₂ = 0.4', mode: 'spin',
    inputs: { ...SPIN, P1: 0.6, P2: 0.4 }, pub: 'mag.TMR', expected: 2 * 0.24 / (1 - 0.24), tol: 1e-6, kind: 'analytical',
    reference: 'Julliere model' },
  { id: 's13-mtj-parallel-resistance', suite: 13, title: 'Junction resistance from the resistance–area product', mode: 'spin',
    inputs: { ...SPIN }, kpi: 'RP', expected: 10e-12 / area, tol: 1e-6, kind: 'analytical',
    reference: 'R_P = RA/(πd²/4)' },
  { id: 's13-mtj-thermal-stability', suite: 13, title: 'Perpendicular free layer: thermal stability factor Δ = K_eff V/k_BT', mode: 'spin',
    inputs: { ...SPIN }, kpi: 'Delta', expected: Delta, tol: 1e-4, kind: 'analytical',
    reference: 'Ikeda et al., "A perpendicular-anisotropy CoFeB–MgO magnetic tunnel junction" — K_eff = K_i/t − µ₀M_s²/2' },
  { id: 's13-mtj-critical-current', suite: 13, title: 'Macrospin switching current I_c0 = (4eα/ħP) Δ k_BT', mode: 'spin',
    inputs: { ...SPIN }, kpi: 'Ic0', expected: 4 * QE * 0.01 * Delta * kT / (HBAR * 0.6) * 1e6, tol: 1e-4, kind: 'analytical',
    reference: 'Sun, "Spin-current interaction with a monodomain magnetic body" — perpendicular free layer with spin-torque efficiency P' },
  { id: 's13-cpp-gmr-two-current', suite: 13, title: 'CPP-GMR of a trilayer in the two-current limit: ΔR/R_P = β²/(1 − β²)', mode: 'spin',
    inputs: { ...SPIN, beta: 0.5, rho_F: 300, t_F: 50, rho_N: 0.5, t_N: 0.5, mfp: 200, tau_sf: 10000 }, pub: 'mag.GMR', expected: 0.25 / 0.75, tol: 1e-3, kind: 'analytical',
    reference: 'Valet & Fert, "Theory of the perpendicular magnetoresistance in magnetic multilayers" — series-resistor limit (layers thin against the spin-diffusion length)',
    note: 'Spacer resistance and spin relaxation made negligible (0.5 nm spacer, spin-diffusion length 32 µm).' },
  { id: 's13-spin-diffusion-length', suite: 13, title: 'Spin-diffusion length √(v_F λ τ_sf/3) of a copper spacer', mode: 'spin',
    inputs: { ...SPIN, mfp: 20, tau_sf: 10 }, pub: 'mag.l_sf', expected: Math.sqrt(1.57e6 * 20e-9 / 3 * 10e-12), tol: 1e-4, kind: 'analytical',
    reference: 'Diffusive spin transport, l_sf = √(Dτ_sf) with D = v_Fλ/3 and the Fermi velocity of copper 1.57×10⁶ m/s (Ashcroft & Mermin, Solid State Physics)' },
  // ---------------------------------------------------------------- published experiments
  { id: 's13-exp-tmr-co-al2o3-nife-4k', suite: 13, title: 'Experiment: tunnel magnetoresistance of a Co/Al₂O₃/Ni₈₀Fe₂₀ junction at 4.2 K against the Julliere value from measured spin polarisations', mode: 'spin',
    inputs: { mode: 'spin', P1: 0.35, P2: 0.45 }, kpi: 'TMR0', expected: 100 * 0.273 / (1 - 0.273), tol: 0.08, kind: 'experiment',
    measured: { value: 27.3, uncertainty: null, unit: '% (ΔR relative to the antiparallel-state resistance)' },
    source: { authors: 'Ziese', title: 'Extrinsic magnetotransport phenomena in ferromagnetic oxides', journal: 'Rep. Prog. Phys.', year: 2002, doi: '10.1088/0034-4885/65/2/202', url: 'https://arxiv.org/abs/cond-mat/0111263' },
    conditions: 'Co/Al₂O₃/Ni₈₀Fe₂₀ tunnel junction at 4.2 K and low bias (measurement of Moodera, Nowak and van de Veerdonk, Phys. Rev. Lett. 80, 2941, 1998); electrode spin polarisations 35 % (Co) and 45 % (NiFe) from superconducting tunnelling spectroscopy through Al₂O₃',
    match: ['P1', 'P2'],
    reference: 'Moodera, Nowak & van de Veerdonk, Phys. Rev. Lett. 80, 2941 (1998), as reported in Ziese, Rep. Prog. Phys. 65, 143 (2002)',
    note: 'Values from the section on ferromagnetic tunnelling junctions of the review (retrieved as arXiv:cond-mat/0111263): measured 20.2 %, 27.1 % and 27.3 % at 295 K, 77 K and 4.2 K; P_Co = 35 % and P_NiFe = 45 %; Julliere value 27.2 %. The review compares the measurement with 2PP′/(1 + PP′), i.e. ΔR over the antiparallel resistance; the suite reports ΔR/R_P = 2P₁P₂/(1 − P₁P₂), so the measured 27.3 % is converted to 0.273/(1 − 0.273) = 37.6 %. The primary paper was not retrieved; no uncertainty is given. Tolerance: a change of 0.01 in either polarisation moves the result by about 3 %.' },
  { id: 's13-exp-yig-sphere-kittel-mode', suite: 13, title: 'Experiment: uniform (Kittel) mode of a YIG sphere tuned to a 10.565 GHz cavity at about 370 mT', mode: 'fmr',
    inputs: { mode: 'fmr', material: 'YIG', shape: 'sphere', H_bias: 0.370 / MU0 / 1e3 }, kpi: 'fK', expected: 10.565, tol: 0.03, kind: 'experiment',
    measured: { value: 10.565, uncertainty: null, unit: 'GHz' },
    source: { authors: 'Tabuchi, Ishino, Ishikawa, Yamazaki, Usami, Nakamura', title: 'Hybridizing Ferromagnetic Magnons and Microwave Photons in the Quantum Limit', journal: 'Phys. Rev. Lett.', year: 2014, doi: '10.1103/PhysRevLett.113.083603', url: 'https://arxiv.org/abs/1405.1913' },
    conditions: 'Single-crystal YIG sphere, 0.5 mm diameter, in a copper microwave cavity (TE101 mode at 10.565 GHz); static field of approximately 370 mT along ⟨100⟩ brings the Kittel mode into resonance with the cavity; millikelvin temperature',
    match: ['H_bias', 'material', 'shape'],
    reference: 'Tabuchi et al., Phys. Rev. Lett. 113, 083603 (2014)',
    note: 'Cavity frequency and static field from the description of the experimental set-up in the text (arXiv:1405.1913). The field is stated only as "approximately 370 mT" (two to three significant figures, about ±1.5 %). The suite adds its room-temperature uniaxial anisotropy field (8.7 mT) to the applied field, whereas the experiment has the field along the cubic hard axis at low temperature; this sign difference is within the tolerance. For a sphere the demagnetising terms cancel, so the test is of the gyromagnetic ratio and the Kittel sphere relation; the temperature input does not enter this mode.' },
];
