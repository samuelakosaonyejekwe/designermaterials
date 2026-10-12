// Reference benchmarks for Suite 20 — Active tuning, programmability and control.
// Every expected value is computed here from a closed-form textbook result and standard constants; nothing is taken
// from the suite's own output.
const PI = Math.PI, EPS0 = 8.8541878128e-12, MU0 = 1.25663706212e-6, C0 = 299792458, Z0 = 376.730313668, R_GAS = 8.314462618;
const GAMMA_BAR = 28.025e9; // electron gyromagnetic ratio γ/2π, Hz/T (g ≈ 2.002)

// ---- PID: first-order-plus-dead-time plant K e^(−θs)/(τs + 1) built from a lumped heat balance ----
const RHO_CP_WATER = 1000 * 4186; // J/(m³·K)
const pidPlant = { V_reac: 1e-3, C_vessel: 1500, UA: 5, P_heater: 1000, dead: 200 };
const tauP = (RHO_CP_WATER * pidPlant.V_reac + pidPlant.C_vessel) / pidPlant.UA, Kp = pidPlant.P_heater / pidPlant.UA, thP = pidPlant.dead;
// saturated heat-up of the default 1 L reactor: T − T_amb = K(1 − e^(−t/τ)) at full power
const tauD = (RHO_CP_WATER * 1e-3 + 1500) / 2, KD = 400 / 2, stepD = 90 - 25;
const riseD = tauD * Math.log((KD - 0.1 * stepD) / (KD - 0.9 * stepD));

// ---- MPC: Arrhenius recipe temperature for half the reference growth rate ----
const mpcT = 1 / (1 / 363.15 + R_GAS * Math.log(2) / 60e3) - 273.15;

// ---- LQR: ẍ + ẋ/Q + x = u, weights q₁ (position), q₂ (velocity), r; symmetric root locus gives
// closed-loop polynomial s² + b s + c with c² = 1 + q₁/r and b² = 2c − 2 + 1/Q² + q₂/r ----
const lq = { f0: 1e4, Q: 200, q1: 15, q2: 1, r: 1 };
const cL = Math.sqrt(1 + lq.q1 / lq.r), bL = Math.sqrt(2 * cL - 2 + 1 / lq.Q ** 2 + lq.q2 / lq.r), w0L = 2 * PI * lq.f0;

// ---- varactor ----
const va = { L: 10e-9, Cg: 0.5e-12, Cj0: 0.5e-12, Vbi: 0.5, M: 0.5, V: 1.5, Vmax: 12, Rs: 2 };
const CjV = V => va.Cj0 / Math.pow(1 + V / va.Vbi, va.M), fV = V => 1 / (2 * PI * Math.sqrt(va.L * (va.Cg + CjV(V))));

// ---- piezo stack ----
const pz = { n: 200, t: 50e-6, A: 49e-6, Y: 60e9, k: 0.6, d33: 500e-12, V: 60, eps: 2000, m: 0.5, alpha: 3 };
const kscP = pz.Y * pz.A / (pz.n * pz.t);

// ---- liquid crystal ----
const lc = { K: 6.2e-12, de: 13, g1: 0.08, d: 4e-6, no: 1.52, ne: 1.74, lam: 1550e-9 };

// ---- MEMS parallel plate (negligible dielectric: t_d/ε_d = 10 pm) ----
const me = { k: 5, A: 4e4 * 1e-12, g: 3e-6 + 1e-9 / 100, f0: 2e4 };
const VpiM = Math.sqrt(8 * me.k * me.g ** 3 / (27 * EPS0 * me.A));
// voltage that places the stable equilibrium at x = g/6: V² = 2k x (g − x)²/(ε₀A)  ⇒  V/V_pi = √(25/32)
const VsixthM = VpiM * Math.sqrt(25 / 32);

// ---- ferrite (Polder tensor, wave ⟂ bias) ----
const fe = { Ms: 1.4e5, H: 1e5, f: 6e9 };
const f0F = GAMMA_BAR * MU0 * fe.H, fmF = GAMMA_BAR * MU0 * fe.Ms;
const muEffLossless = ((f0F + fmF) ** 2 - fe.f ** 2) / (f0F * (f0F + fmF) - fe.f ** 2);
const hF = fe.f / (GAMMA_BAR * MU0), HneedF = (-fe.Ms + Math.sqrt(fe.Ms ** 2 + 4 * hF * hF)) / 2;

// ---- von Mises truss ----
const vm = { kb: 5e4, w: 20e-3, h: 3e-3 };
const L0v = Math.hypot(vm.w, vm.h), uStar = Math.sqrt(Math.pow(L0v * vm.w * vm.w, 2 / 3) - vm.w * vm.w);
const FsnapV = 2 * vm.kb * uStar * (L0v / Math.cbrt(L0v * vm.w * vm.w) - 1);

// ---- quasi-zero-stiffness mount ----
const qz = { kv: 2000, kh: 1500, L0: 0.1, a: 0.0625, m: 5, zeta: 0.02, f: 2 };
const K0q = qz.kv - 2 * qz.kh * (qz.L0 / qz.a - 1), fnq = Math.sqrt(K0q / qz.m) / (2 * PI);
const zq = qz.zeta * Math.sqrt(qz.kv / K0q), rq = qz.f / fnq;
const Tq = Math.sqrt((1 + (2 * zq * rq) ** 2) / ((1 - rq * rq) ** 2 + (2 * zq * rq) ** 2));

// ---- phase-change cell ----
const gstTau = 5 * 6150 * 218 * (50e-12 * 200e-9) / 1e-5; // heated-volume factor × ρ c_p V / G_th

export default [
  // ----- PID -----
  {
    id: 's20-pid-ziegler-nichols-gain', suite: 20, title: 'Ziegler–Nichols open-loop (reaction-curve) PID gain', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'zn', ...pidPlant, noise: 0 },
    kpi: 'Kc', expected: 1.2 * tauP / (Kp * thP), tol: 0.002, kind: 'analytical',
    reference: 'Ziegler & Nichols open-loop tuning rule, K_c = 1.2 τ/(K θ) (Seborg, Edgar, Mellichamp & Doyle, Process Dynamics and Control)',
    note: 'Plant from the heat balance: τ = (ρc_pV + C_vessel)/UA with water at 1000 kg/m³ and 4186 J/(kg·K), K = P/UA.'
  },
  {
    id: 's20-pid-cohen-coon-gain', suite: 20, title: 'Cohen–Coon PID gain', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'cc', ...pidPlant, noise: 0 },
    kpi: 'Kc', expected: tauP / (Kp * thP) * (4 / 3 + thP / (4 * tauP)), tol: 0.002, kind: 'analytical',
    reference: 'Cohen–Coon tuning rule for a PID controller, K_c = (τ/Kθ)(4/3 + θ/4τ) (Seborg, Edgar, Mellichamp & Doyle, Process Dynamics and Control)'
  },
  {
    id: 's20-pid-simc-gain', suite: 20, title: 'SIMC PI gain with τ_c = θ', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'simc', lam: 1, ...pidPlant, noise: 0 },
    kpi: 'Kc', expected: tauP / (Kp * 2 * thP), tol: 0.002, kind: 'analytical',
    reference: 'Skogestad, Simple analytic rules for model reduction and PID controller tuning (SIMC): K_c = τ/[K(τ_c + θ)]'
  },
  {
    id: 's20-pid-simc-gain-margin', suite: 20, title: 'SIMC PI with τ_c = θ and T_i = τ: gain margin π', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'simc', lam: 1, ...pidPlant, Ts: 0.1, t_sim: 2000, noise: 0, dist: 0 },
    kpi: 'GM', expected: 20 * Math.log10(PI), tol: 0.003, kind: 'analytical',
    reference: 'Skogestad SIMC rule: for τ ≤ 8θ the loop is L = e^(−θs)/(2θs), giving GM = π (3.14)',
    note: 'τ = 1137 s < 8θ = 1600 s so the integral time cancels the plant pole exactly. The suite adds a half-sample delay (0.05 s on 200 s), worth 0.002 dB.'
  },
  {
    id: 's20-pid-simc-phase-margin', suite: 20, title: 'SIMC PI with τ_c = θ and T_i = τ: phase margin 61.4°', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'simc', lam: 1, ...pidPlant, Ts: 0.1, t_sim: 2000, noise: 0, dist: 0 },
    kpi: 'PM', expected: 90 - 0.5 * 180 / PI, tol: 0.003, kind: 'analytical',
    reference: 'Skogestad SIMC rule: L = e^(−θs)/(2θs) crosses unity at ω = 1/(2θ), PM = 90° − 28.6° = 61.4°'
  },
  {
    id: 's20-pid-simc-peak-sensitivity', suite: 20, title: 'SIMC PI with τ_c = θ: peak sensitivity M_s = 1.59', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'simc', lam: 1, ...pidPlant, Ts: 0.1, t_sim: 2000, noise: 0, dist: 0 },
    kpi: 'Ms', expected: 1.59, tol: 0.01, kind: 'numerical-reference',
    reference: 'Skogestad, Simple analytic rules for model reduction and PID controller tuning: robustness table for τ_c = θ (GM 3.14, PM 61.4°, M_s 1.59)',
    note: 'Reference quoted to three figures; the suite evaluates the maximum on a 500-point logarithmic frequency grid.'
  },
  {
    id: 's20-pid-pure-delay-gain-margin', suite: 20, title: 'Proportional loop around a pure delay: gain margin 1/(K K_c)', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'manual', V_reac: 1e-5, C_vessel: 0, UA: 1e5, P_heater: 5e6, dead: 100, Kc_man: 0.01, Ti_man: 1e6, Td_man: 0, Ts: 1, t_sim: 7200, noise: 0, dist: 0, T_set: 60, T_amb: 25 },
    kpi: 'GM', expected: 20 * Math.log10(2), tol: 0.005, kind: 'analytical',
    reference: 'Nyquist criterion for K K_c e^(−θs): |L| = K K_c at every frequency, so GM = 1/(K K_c) = 2 (6.02 dB) (Åström & Murray, Feedback Systems)',
    note: 'Plant time constant 0.4 ms is negligible against the 100 s delay; integral time 10⁶ s makes the controller effectively proportional.'
  },
  {
    id: 's20-pid-saturated-rise-time', suite: 20, title: 'Power-limited heat-up: 10–90 % rise time of the first-order plant', mode: 'pid',
    inputs: { mode: 'pid', tuning: 'simc', noise: 0, Ts: 1 },
    kpi: 't_rise', expected: riseD, tol: 0.005, kind: 'analytical',
    reference: 'First-order step response T − T_amb = K(1 − e^(−t/τ)) at full heater power: t_rise = τ ln[(K − 0.1ΔT)/(K − 0.9ΔT)]',
    note: 'Default 1 L reactor (K = 200 K, τ = 2843 s, 65 K step); the heater stays saturated until beyond 90 % of the step. Resolution limited to the 1 s sampling interval.'
  },
  // ----- MPC -----
  {
    id: 's20-mpc-arrhenius-recipe-temperature', suite: 20, title: 'Recipe temperature for half the reference growth rate (Arrhenius inversion)', mode: 'mpc',
    inputs: { mode: 'mpc', manip: 'temp', L_target: 1e-6, t_growth: 10000, G_ref: 2e-10, Ea: 60, T_lo: 60, T_hi: 97, bias: 0, deplete: 0, len_noise: 0, tau_act: 1 },
    kpi: 'z_mean', expected: mpcT, abs: 0.05, kind: 'analytical',
    reference: 'Arrhenius law: 1/T = 1/T_ref + (R/E_a) ln(G_ref/G) with T_ref = 363.15 K, E_a = 60 kJ/mol, G = G_ref/2 (77.8 °C)',
    note: 'No plant–model mismatch, depletion or noise, so the controller holds the nominal recipe; absolute tolerance in kelvin.'
  },
  {
    id: 's20-mpc-nominal-length', suite: 20, title: 'Nominal growth without mismatch reaches the target length', mode: 'mpc',
    inputs: { mode: 'mpc', manip: 'dose', L_target: 1e-6, t_growth: 10000, G_ref: 2e-10, bias: 0, deplete: 0, len_noise: 0, tau_act: 1 },
    kpi: 'L_final', expected: 1e-6, tol: 0.002, kind: 'conservation',
    reference: 'Integral of a constant growth rate: L = G·t, with first-order concentration dependence G = G_ref·c/c_ref and c/c_ref = L/(G_ref t) = 0.5'
  },
  // ----- LQR -----
  {
    id: 's20-lqr-closed-loop-frequency', suite: 20, title: 'LQR of an oscillator: closed-loop natural frequency ω₀(1 + q₁/r)^¼', mode: 'lqg',
    inputs: { mode: 'lqg', f0_r: lq.f0, Q_r: lq.Q, q_pos: lq.q1, q_vel: lq.q2, r_u: lq.r },
    kpi: 'f_cl', expected: lq.f0 * Math.sqrt(cL), tol: 1e-4, kind: 'analytical',
    reference: 'Symmetric root locus of the LQR problem (Kailath, Linear Systems; Anderson & Moore, Optimal Control): closed-loop c² = 1 + q₁/r for ẍ + ẋ/Q + x = u',
    note: 'q₁/r = 15 doubles the resonance frequency exactly.'
  },
  {
    id: 's20-lqr-closed-loop-damping', suite: 20, title: 'LQR of an oscillator: closed-loop damping ratio', mode: 'lqg',
    inputs: { mode: 'lqg', f0_r: lq.f0, Q_r: lq.Q, q_pos: lq.q1, q_vel: lq.q2, r_u: lq.r },
    kpi: 'zeta_cl', expected: bL / (2 * Math.sqrt(cL)), tol: 1e-4, kind: 'analytical',
    reference: 'Symmetric root locus of the LQR problem: b² = 2c − 2 + 1/Q² + q₂/r, ζ = b/(2√c) (Anderson & Moore, Optimal Control: Linear Quadratic Methods)'
  },
  {
    id: 's20-lqr-open-loop-ringdown', suite: 20, title: 'Free ring-down of a resonator to 2 % amplitude', mode: 'lqg',
    inputs: { mode: 'lqg', f0_r: lq.f0, Q_r: lq.Q },
    kpi: 'ts_open', expected: 2 * lq.Q * Math.log(50) / w0L, tol: 1e-4, kind: 'analytical',
    reference: 'Damped harmonic oscillator: amplitude envelope exp(−ω₀t/2Q), so t(2 %) = 2Q ln 50/ω₀'
  },
  // ----- varactor -----
  {
    id: 's20-varactor-lc-resonance', suite: 20, title: 'Varactor-tuned LC resonance at a bias of 3 V_bi (abrupt junction)', mode: 'varactor',
    inputs: { mode: 'varactor', L_srr: 10, C_gap: 0.5, Cj0: 0.5, Vbi: 0.5, Mj: 0.5, V_bias: 1.5, V_max: 12, R_s: 2 },
    kpi: 'f_bias', expected: 1 / (2 * PI * Math.sqrt(va.L * (va.Cg + va.Cj0 / 2))), tol: 1e-4, kind: 'analytical',
    reference: 'Thomson resonance formula f = 1/(2π√(LC)) with the abrupt-junction depletion capacitance C_j = C_j0/√(1 + V/V_bi) (Sze, Physics of Semiconductor Devices)',
    note: 'At V = 3 V_bi the junction capacitance is exactly C_j0/2, so C = 0.75 pF and f = 1.838 GHz.'
  },
  {
    id: 's20-varactor-tuning-range', suite: 20, title: 'Varactor tuning range between 0 V and 12 V', mode: 'varactor',
    inputs: { mode: 'varactor', L_srr: 10, C_gap: 0.5, Cj0: 0.5, Vbi: 0.5, Mj: 0.5, V_bias: 1.5, V_max: 12, R_s: 2 },
    kpi: 'range', expected: 100 * 2 * (fV(12) - fV(0)) / (fV(12) + fV(0)), tol: 1e-4, kind: 'analytical',
    reference: 'f ∝ 1/√C with C_j(12 V) = C_j0/5 for V_bi = 0.5 V: f_max/f_min = √(1.0/0.6)'
  },
  {
    id: 's20-varactor-quality-factor', suite: 20, title: 'Series-RLC quality factor √(L/C)/R', mode: 'varactor',
    inputs: { mode: 'varactor', L_srr: 10, C_gap: 0.5, Cj0: 0.5, Vbi: 0.5, Mj: 0.5, V_bias: 0, V_max: 12, R_s: 2 },
    kpi: 'Q', expected: Math.sqrt(10e-9 / 1e-12) / 2, tol: 1e-4, kind: 'analytical',
    reference: 'Series RLC circuit: Q = (1/R)√(L/C) = 100 Ω / 2 Ω = 50 (Pozar, Microwave Engineering)'
  },
  // ----- phase-change -----
  {
    id: 's20-pcm-shunt-sheet-transmission', suite: 20, title: 'Transmission of a resistive shunt sheet at series resonance', mode: 'pcm',
    inputs: { mode: 'pcm', pcm_mat: 'vo2', thick_nm: 5, R_ms: 100, X0: 300, f0_ms: 1e12 },
    kpi: 'T_off', expected: (200 / (200 + Z0)) ** 2, tol: 0.005, kind: 'analytical',
    reference: 'Transmission-line theory: a shunt sheet impedance Z_s between two free-space half-spaces transmits T = |2Z_s/(2Z_s + Z₀)|², Z₀ = 376.73 Ω (Pozar, Microwave Engineering)',
    note: 'At resonance Z_s = R = 100 Ω. The 5 nm insulating VO₂ film (sheet resistance 4 MΩ per square) across the gap adds about 0.02 Ω, hence the 0.5 % tolerance.'
  },
  {
    id: 's20-pcm-thermal-rc-rise', suite: 20, title: 'Thermal RC heating: 63.2 % of the steady rise is reached after one time constant', mode: 'pcm',
    inputs: { mode: 'pcm', pcm_mat: 'gst', area_um2: 50, thick_nm: 200, G_th: 1e-5, Cth_mult: 5, T_base: 25, P_mw: 1e3 * 1e-5 * 140 / (1 - Math.exp(-1)) },
    kpi: 't_on', expected: gstTau, tol: 0.002, kind: 'analytical',
    reference: 'First-order step response 1 − e^(−t/τ) = 0.632 at t = τ = C_th/G_th (Incropera & DeWitt, lumped-capacitance method)',
    note: 'Heater power chosen so that the switching temperature (165 °C in the model, 140 K above base) is 1 − 1/e of the steady rise. C_th uses the model\'s typical Ge₂Sb₂Te₅ film values ρ = 6150 kg/m³ and c_p = 218 J/(kg·K) times the heated-volume factor 5; no latent heat on crystallisation.'
  },
  // ----- piezo -----
  {
    id: 's20-piezo-free-stroke', suite: 20, title: 'Multilayer piezo stack free stroke n·d₃₃·V', mode: 'piezo',
    inputs: { mode: 'piezo', d33: 500, n_lay: 200, t_lay: 50, A_pz: 49, Y_pz: 60, k33: 0.6, epsr: 2000, V_pz: 60, m_pz: 0.5, alpha_sh: 3 },
    kpi: 'stroke', expected: pz.n * pz.d33 * pz.V, tol: 1e-6, kind: 'analytical',
    reference: 'Converse piezoelectric effect for n layers electrically in parallel: Δ = n d₃₃ V (6 µm here) (Jaffe, Cook & Jaffe, Piezoelectric Ceramics)'
  },
  {
    id: 's20-piezo-blocked-force', suite: 20, title: 'Piezo stack blocked force', mode: 'piezo',
    inputs: { mode: 'piezo', d33: 500, n_lay: 200, t_lay: 50, A_pz: 49, Y_pz: 60, k33: 0.6, epsr: 2000, V_pz: 60, m_pz: 0.5, alpha_sh: 3 },
    kpi: 'F_block', expected: pz.Y * pz.A * pz.d33 * pz.V / pz.t, tol: 1e-6, kind: 'analytical',
    reference: 'Blocked force F_b = Y A d₃₃ E with E = V/t (stiffness YA/L times free stroke): 1764 N'
  },
  {
    id: 's20-piezo-short-open-tuning', suite: 20, title: 'Short- to open-circuit resonance shift for k₃₃ = 0.6', mode: 'piezo',
    inputs: { mode: 'piezo', d33: 500, n_lay: 200, t_lay: 50, A_pz: 49, Y_pz: 60, k33: 0.6, epsr: 2000, V_pz: 60, m_pz: 0.5, alpha_sh: 3 },
    kpi: 'f_range', expected: 25, tol: 1e-6, kind: 'analytical',
    reference: 'Electromechanical coupling: k_oc/k_sc = 1/(1 − k₃₃²), so f_oc/f_sc = 1/√(1 − 0.36) = 1.25 (IEEE Standard on Piezoelectricity definitions)'
  },
  {
    id: 's20-piezo-capacitive-shunt-frequency', suite: 20, title: 'Resonance of a mass on a capacitively shunted stack (C_sh = 3 C_S)', mode: 'piezo',
    inputs: { mode: 'piezo', d33: 500, n_lay: 200, t_lay: 50, A_pz: 49, Y_pz: 60, k33: 0.6, epsr: 2000, V_pz: 60, m_pz: 0.5, alpha_sh: 3 },
    kpi: 'f_sh', expected: Math.sqrt(kscP * (1 + (0.36 / 0.64) / 4) / pz.m) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Capacitive shunt of a piezoelectric element: k = k_sc[1 + k²/(1 − k²)·C_S/(C_S + C_sh)] (Hagood & von Flotow shunt model; Preumont, Mechatronics)'
  },
  {
    id: 's20-piezo-stack-capacitance', suite: 20, title: 'Free capacitance of the multilayer stack', mode: 'piezo',
    inputs: { mode: 'piezo', d33: 500, n_lay: 200, t_lay: 50, A_pz: 49, Y_pz: 60, k33: 0.6, epsr: 2000, V_pz: 60, m_pz: 0.5, alpha_sh: 3 },
    kpi: 'C', expected: pz.n * EPS0 * pz.eps * pz.A / pz.t, tol: 1e-4, kind: 'analytical',
    reference: 'Parallel-plate capacitance of n layers in parallel, C = n ε₀ ε_r A/t (3.47 µF)'
  },
  // ----- liquid crystal -----
  {
    id: 's20-lc-freedericksz-threshold', suite: 20, title: 'Fréedericksz threshold voltage of a planar nematic cell', mode: 'lc',
    inputs: { mode: 'lc', K11: 6.2, d_eps: 13, gamma1: 0.08, d_lc: 4 },
    kpi: 'V_th', expected: PI * Math.sqrt(lc.K / (EPS0 * lc.de)), tol: 1e-4, kind: 'analytical',
    reference: 'Fréedericksz transition, V_th = π√(K₁₁/(ε₀Δε)) (de Gennes & Prost, The Physics of Liquid Crystals)',
    note: 'K₁₁ = 6.2 pN and Δε = 13 are of the order of the values of 5CB; V_th ≈ 0.73 V.'
  },
  {
    id: 's20-lc-relaxation-time', suite: 20, title: 'Nematic free-relaxation time γ₁d²/(K₁₁π²)', mode: 'lc',
    inputs: { mode: 'lc', K11: 6.2, d_eps: 13, gamma1: 0.08, d_lc: 4 },
    kpi: 't_off', expected: lc.g1 * lc.d ** 2 / (lc.K * PI * PI), tol: 1e-4, kind: 'analytical',
    reference: 'Director relaxation of a planar cell with strong anchoring (de Gennes & Prost, The Physics of Liquid Crystals; Yang & Wu, Fundamentals of Liquid Crystal Devices)'
  },
  {
    id: 's20-lc-saturated-phase', suite: 20, title: 'Saturated phase modulation 2πΔn·d/λ', mode: 'lc',
    inputs: { mode: 'lc', n_o: 1.52, n_e: 1.74, d_lc: 4, lam_nm: 1550, V_lc: 100 },
    kpi: 'phase', expected: 2 * PI * (lc.ne - lc.no) * lc.d / lc.lam, tol: 1e-3, kind: 'analytical',
    reference: 'Retardation of a uniaxial layer switched from n_e to n_o: Δφ = 2π(n_e − n_o)d/λ (Yariv & Yeh, Optical Waves in Crystals)',
    note: 'At 100 V rms the director is fully reoriented in the single-angle model.'
  },
  // ----- MEMS -----
  {
    id: 's20-mems-pull-in-voltage', suite: 20, title: 'Parallel-plate electrostatic pull-in voltage', mode: 'mems',
    inputs: { mode: 'mems', k_m: 5, A_m: 4e4, g0: 3, t_d: 1, eps_d: 100, V_m: 1, f_mech: 2e4 },
    kpi: 'V_pi', expected: VpiM, tol: 1e-4, kind: 'analytical',
    reference: 'Pull-in of a spring-suspended parallel plate, V_pi = √(8k g³/(27ε₀A)) at x = g/3 (Senturia, Microsystem Design)',
    note: 'Dielectric made electrically negligible (t_d/ε_d = 10 pm); the expected value includes it in g.'
  },
  {
    id: 's20-mems-deflection-sixth-gap', suite: 20, title: 'Static deflection g/6 at V = V_pi·√(25/32)', mode: 'mems',
    inputs: { mode: 'mems', k_m: 5, A_m: 4e4, g0: 3, t_d: 1, eps_d: 100, V_m: VsixthM, f_mech: 2e4 },
    kpi: 'x', expected: me.g / 6, tol: 1e-4, kind: 'analytical',
    reference: 'Exact root of k x = ε₀AV²/(2(g − x)²): x = g/6 requires V² = (25/108)·k g³/(ε₀A) = (25/32)V_pi² (Senturia, Microsystem Design)'
  },
  {
    id: 's20-mems-spring-softening', suite: 20, title: 'Electrostatic spring softening at x = g/6 (negative electrostatic stiffness)', mode: 'mems',
    inputs: { mode: 'mems', k_m: 5, A_m: 4e4, g0: 3, t_d: 1, eps_d: 100, V_m: VsixthM, f_mech: 2e4 },
    kpi: 'f_V', expected: me.f0 * Math.sqrt(0.6), tol: 1e-4, kind: 'analytical',
    reference: 'k_eff = k − ε₀AV²/(g − x)³ = k[1 − 2x/(g − x)]; at x = g/6 the electrostatic stiffness is −0.4k and f = f₀√0.6'
  },
  // ----- ferrite -----
  {
    id: 's20-ferrite-kittel-resonance', suite: 20, title: 'Ferromagnetic resonance of transversely biased YIG (Kittel)', mode: 'ferrite',
    inputs: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5 },
    kpi: 'f_res', expected: GAMMA_BAR * MU0 * Math.sqrt(fe.H * (fe.H + fe.Ms)), tol: 1e-4, kind: 'analytical',
    reference: 'Kittel formula f = (γμ₀/2π)√(H(H + M_s)) with γ/2π = 28.025 GHz/T (Pozar, Microwave Engineering; Lax & Button, Microwave Ferrites and Ferrimagnetics)'
  },
  {
    id: 's20-ferrite-anti-resonance', suite: 20, title: 'Anti-resonance (μ_eff = 0) of the transversely biased ferrite', mode: 'ferrite',
    inputs: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5 },
    kpi: 'f_ar', expected: GAMMA_BAR * MU0 * (fe.H + fe.Ms), tol: 1e-4, kind: 'analytical',
    reference: 'Zero of μ_eff = (μ² − κ²)/μ at ω = ω₀ + ω_m (Pozar, Microwave Engineering)'
  },
  {
    id: 's20-ferrite-negative-permeability', suite: 20, title: 'Negative effective permeability between resonance and anti-resonance', mode: 'ferrite',
    inputs: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5, f_op: 6e9, alpha_g: 1e-5 },
    kpi: 'mu_re', expected: muEffLossless, tol: 0.002, kind: 'analytical',
    reference: 'Lossless Polder tensor: μ_eff = [(ω₀ + ω_m)² − ω²]/[ω₀(ω₀ + ω_m) − ω²] (Pozar, Microwave Engineering)',
    note: 'Gilbert damping set to 10⁻⁵ so that loss does not alter the real part; μ_eff ≈ −5.7 at 6 GHz.'
  },
  {
    id: 's20-ferrite-static-permeability', suite: 20, title: 'Low-frequency limit μ_eff → 1 + M_s/H', mode: 'ferrite',
    inputs: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5, f_op: 1e7 },
    kpi: 'mu_re', expected: 1 + fe.Ms / fe.H, tol: 1e-4, kind: 'analytical',
    reference: 'Static susceptibility of a saturated ferrite transverse to the bias, χ = M_s/H'
  },
  {
    id: 's20-ferrite-bias-for-resonance', suite: 20, title: 'Bias field placing the Kittel resonance at 6 GHz', mode: 'ferrite',
    inputs: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5, f_op: 6e9 },
    kpi: 'H_need', expected: HneedF, tol: 1e-4, kind: 'analytical',
    reference: 'Inversion of the Kittel formula: H = ½[−M_s + √(M_s² + 4(ω/γμ₀)²)]'
  },
  // ----- bistable -----
  {
    id: 's20-bistable-snap-through-force', suite: 20, title: 'Von Mises truss: exact snap-through (limit) load', mode: 'bistable',
    inputs: { mode: 'bistable', k_bar: 5e4, w_b: 20, h_b: 3, k_v: 0 },
    kpi: 'F_snap', expected: FsnapV, tol: 1e-4, kind: 'analytical',
    reference: 'Two-bar (von Mises) truss with linear bars: F = 2k u(L₀/√(w² + u²) − 1) is extremal at w² + u² = (L₀w²)^(2/3) (Bažant & Cedolin, Stability of Structures)',
    note: 'Shallow-truss estimate 2k h³/(3√3 w²) = 1.299 N; the exact value is slightly lower.'
  },
  {
    id: 's20-bistable-energy-barrier', suite: 20, title: 'Von Mises truss: energy barrier equals the bar strain energy in the flat state', mode: 'bistable',
    inputs: { mode: 'bistable', k_bar: 5e4, w_b: 20, h_b: 3, k_v: 0 },
    kpi: 'E_bar', expected: vm.kb * (L0v - vm.w) ** 2, tol: 1e-3, kind: 'conservation',
    reference: 'Work–energy theorem: ∫₀^h F dx = 2·½k(L₀ − w)², the strain energy of both bars when the truss is flat',
    note: 'The suite integrates the force curve with a 200-point trapezoidal rule.'
  },
  {
    id: 's20-bistable-negative-stiffness', suite: 20, title: 'Von Mises truss: negative tangent stiffness in the flat configuration', mode: 'bistable',
    inputs: { mode: 'bistable', k_bar: 5e4, w_b: 20, h_b: 3, k_v: 0 },
    kpi: 'k_neg', expected: -2 * vm.kb * (L0v / vm.w - 1), tol: 1e-3, kind: 'analytical',
    reference: 'dF/dx at u = 0: compressed bars of force k(L₀ − w) give a geometric stiffness −2k(L₀ − w)/w (Bažant & Cedolin, Stability of Structures)'
  },
  // ----- QZS -----
  {
    id: 's20-qzs-net-stiffness', suite: 20, title: 'Three-spring quasi-zero-stiffness mount: net stiffness at equilibrium', mode: 'qzs',
    inputs: { mode: 'qzs', kq_v: 2000, kq_h: 1500, L0_q: 100, a_q: 62.5, m_q: 5, zeta_q: 0.02, f_exc: 2 },
    kpi: 'K0', expected: 200, tol: 1e-6, kind: 'analytical',
    reference: 'Carrella, Brennan & Waters three-spring QZS model: K₀ = k_v − 2k_h(L₀/a − 1) = 2000 − 3000·0.6 N/m'
  },
  {
    id: 's20-qzs-zero-stiffness-length', suite: 20, title: 'Compressed length giving exactly zero stiffness', mode: 'qzs',
    inputs: { mode: 'qzs', kq_v: 2000, kq_h: 1500, L0_q: 100, a_q: 62.5 },
    kpi: 'a_star', expected: 0.06, tol: 1e-6, kind: 'analytical',
    reference: 'K₀ = 0 ⇒ a* = 2k_hL₀/(k_v + 2k_h) = 60 mm (Carrella, Brennan & Waters QZS condition)'
  },
  {
    id: 's20-qzs-natural-frequency', suite: 20, title: 'Natural frequency reduced by √(K₀/k_v)', mode: 'qzs',
    inputs: { mode: 'qzs', kq_v: 2000, kq_h: 1500, L0_q: 100, a_q: 62.5, m_q: 5 },
    kpi: 'f_n', expected: Math.sqrt(200 / 5) / (2 * PI), tol: 1e-5, kind: 'analytical',
    reference: 'Single-degree-of-freedom oscillator f_n = √(K/m)/2π (Den Hartog, Mechanical Vibrations): 1.0066 Hz against 3.183 Hz for the linear mount'
  },
  {
    id: 's20-qzs-transmissibility', suite: 20, title: 'Transmissibility of the tuned mount at 2 Hz', mode: 'qzs',
    inputs: { mode: 'qzs', kq_v: 2000, kq_h: 1500, L0_q: 100, a_q: 62.5, m_q: 5, zeta_q: 0.02, f_exc: 2 },
    kpi: 'T_q', expected: Tq, tol: 1e-4, kind: 'analytical',
    reference: 'Viscously damped single-degree-of-freedom transmissibility T = √{[1 + (2ζr)²]/[(1 − r²)² + (2ζr)²]} (Den Hartog, Mechanical Vibrations)',
    note: 'The damper is unchanged, so ζ rises to ζ₀√(k_v/K₀) = 0.0632 when the stiffness is cut to 10 %.'
  },
];
