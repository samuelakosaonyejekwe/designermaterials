// Reference benchmarks for Suite 26 — Industrial Scale-Up, Automation and Quality.
// Every expected value is a literal or is computed here from a closed-form expression; nothing is taken from the suite.
const sum = a => a.reduce((s, v) => s + v, 0);

// ---------- ideal reactors ----------
const RG = 8.314462618, M_ZNO = (65.38 + 15.999) / 1e3;                         // kg/mol, standard atomic weights of Zn and O
const KIN1 = { mode: 'reactor', route: 'zno', order: 1, k0: 1e-3, Ea: 0, T: 90, c0: 50, X: 0.9, demand: 10, hours: 6000, isoYield: 92 };   // k = 10⁻³ s⁻¹ at any temperature
const KIN2 = { ...KIN1, order: 2, k0: 1e-5 };                                   // k·c₀ = 5·10⁻⁴ s⁻¹
const Q_FEED = 10e3 / M_ZNO / (50 * 0.9 * 0.92) / (6000 * 3600);                // m³/s of feed for 10 t/yr ZnO (1 mol ZnO per mol Zn salt)
const TAU_CSTR = 0.9 / (1e-3 * (1 - 0.9)), V_CSTR = TAU_CSTR * Q_FEED;          // V = F_A0·X/(k·C_A0·(1 − X))
const LEV_REF = 'Ideal-reactor design equations for nth-order kinetics at constant density — Levenspiel, Chemical Reaction Engineering; Fogler, Elements of Chemical Reaction Engineering';
// laboratory stirred tank of the default recipe: 1 L liquid at 75 % fill, H/D = 1.2, impeller d = 0.4·D, Po = 5, 400 rpm, ρ = 1000 kg/m³, μ = 0.5 mPa·s
const V_LAB = 1e-3, D_LAB = Math.cbrt(4 * (V_LAB / 0.75) / (Math.PI * 1.2)), d_LAB = 0.4 * D_LAB, N_LAB = 400 / 60;
const PV_LAB = 5 * 1000 * N_LAB ** 3 * d_LAB ** 5 / V_LAB, TIP_LAB = Math.PI * N_LAB * d_LAB;
const MIX = { mode: 'reactor', rtype: 'batch', route: 'zno', Vlab: 1, fill: 75, aspect: 1.2, dRatio: 0.4, Np: 5, Nlab: 400, rho: 1000, mu: 5e-4 };
const MIX_REF = 'Stirred-tank scale-up at geometric similarity with P = Po·ρ·N³·d⁵ — Paul, Atiemo-Obeng & Kresta (eds.), Handbook of Industrial Mixing; Nienow, turbulent blend-time correlation N·θ = 5.2·Po^(−1/3)·(T/D)²';

// ---------- statistical process control ----------
// six subgroups of five identical patterns: every range is 4, the grand mean 100
const PAT5 = [98, 99, 100, 101, 102], G5 = 6, D2_5 = 2.326, A2_5 = 0.577;
const spc = (groups) => 'subgroup,value\n' + groups.flatMap((g, i) => g.map(v => `${i + 1},${v}`)).join('\n');
const SPC5 = spc(Array.from({ length: G5 }, () => PAT5)), SW5 = 4 / D2_5, SO5 = Math.sqrt(G5 * 10 / (G5 * 5 - 1));
const SPC_REF = 'Shewhart control-chart constants (d₂ = 2.326, A₂ = 0.577 for subgroups of five; d₂ = 1.128, E₂ = 2.660 for moving ranges of two) and capability-index definitions — Montgomery, Introduction to Statistical Quality Control';
// individuals alternating 10, 12: every moving range is 2
const IMR = spc(Array.from({ length: 20 }, (_, i) => [i % 2 ? 12 : 10]));
// four subgroups of twelve with the same pattern: X̄–S chart, σ = s̄/c₄
const PAT12 = [94, 95, 96, 97, 98, 99, 101, 102, 103, 104, 105, 106], S12 = Math.sqrt(sum(PAT12.map(v => (v - 100) ** 2)) / 11);
const C4_12 = Math.sqrt(2 / 11) * 120 / (945 / 32 * Math.sqrt(Math.PI));        // c₄ = √(2/(n−1))·Γ(n/2)/Γ((n−1)/2) with Γ(6) = 120, Γ(11/2) = (945/32)·√π  → 0.9776

// ---------- two-level factorial designs ----------
const doe = (head, rows) => head.join(',') + '\n' + rows.map(r => r.join(',')).join('\n');
const y2 = (a, b) => 50 + 10 * a + 5 * b + 2.5 * a * b;                         // effects: A = 20, B = 10, AB = 5
const C2 = [[-1, -1], [1, -1], [-1, 1], [1, 1]];
const DOE2 = doe(['A', 'B', 'y'], C2.map(([a, b]) => [a, b, y2(a, b)]));
const DOE2R = doe(['A', 'B', 'y'], [...C2.flatMap(([a, b]) => [[a, b, y2(a, b) + 1], [a, b, y2(a, b) - 1]]), [0, 0, 54], [0, 0, 56]]);   // replicates ±1 about the model, two centre points averaging 55
const y3 = (a, b, c) => 100 + 10 * a + 5 * b + 3 * c + 2 * a * b + a * c + 0.5 * b * c + 0.25 * a * b * c;   // effects 20, 10, 6, 4, 2, 1, 0.5
const DOE3 = doe(['A', 'B', 'C', 'y'], Array.from({ length: 8 }, (_, i) => { const x = [i & 1, (i >> 1) & 1, (i >> 2) & 1].map(v => 2 * v - 1); return [...x, y3(...x)]; }));
const T975_1 = Math.tan(0.475 * Math.PI);                                       // Student-t 97.5 % quantile for 1 degree of freedom (Cauchy) = 12.706
const DOE_REF = 'Effects of two-level factorial designs (effect = ȳ₊ − ȳ₋) and their standard error 2s/√N — Box, Hunter & Hunter, Statistics for Experimenters; Montgomery, Design and Analysis of Experiments';
const LENTH_REF = "Lenth's pseudo standard error for unreplicated factorials, ME = t(0.975; m/3)·PSE — Lenth, Quick and Easy Analysis of Unreplicated Factorials (Technometrics); Student-t quantiles from standard tables";

// ---------- yield, throughput, OEE ----------
const AD = 0.5 * 1 * 0.4;                                                       // defect density × area × critical fraction = 0.2
const YLD = { mode: 'yield', D0: 0.5, area: 1, critFrac: 0.4, alphaC: 2 };
const YLD_REF = 'Defect-limited yield models (Poisson, Murphy, Seeds, negative binomial) — Murphy, Cost-size optima of monolithic integrated circuits (Proc. IEEE); Stapper, Armstrong & Saji, Integrated circuit yield statistics (Proc. IEEE)';
const STEP_HEAD = 'step,cycle_min,machines,yield_pct,availability_pct\n';
const ONE = { mode: 'yield', steps: STEP_HEAD + 'single step,10,1,95,90', perf: 92, D0: 1e-3, area: 1e-3, critFrac: 0.01, ymodel: 'poisson', units: 12000, hours: 6000, tBatch: 6, tClean: 1, tChange: 12, nPerCamp: 10, batchSize: 25 };
const Y_ONE = 0.95 * Math.exp(-1e-8);
const TWO = { mode: 'yield', steps: STEP_HEAD + 'first,10,1,90,100\nsecond,12,1,80,100', perf: 100, D0: 1e-3, area: 1e-3, critFrac: 0.01, ymodel: 'poisson', units: 12000, hours: 6000 };
const LINE_REF = 'Rolled throughput yield, takt time, bottleneck capacity and OEE = availability × performance × quality — Hopp & Spearman, Factory Physics; Nakajima, Introduction to TPM';

// ---------- automation ----------
const SIMC_REF = 'SIMC rules for a first-order-plus-dead-time process: Kc = τ/(K·(τc + θ)), Tᵢ = min(τ, 4·(τc + θ)) — Skogestad, Simple analytic rules for model reduction and PID controller tuning (Journal of Process Control)';
const TS_PM = 0.1, PM_REF = 90 - (0.5 + TS_PM / (4 * 60)) * 180 / Math.PI;    // loop 1/((τc + θ)·s)·e^(−θs) with τc = θ: ω_c = 1/(2θ), PM = 90° − 0.5 rad = 61.4°, less the half-sample hold delay
// crossed gauge study 2 parts × 2 operators × 2 repeats: value = 100 ± 10 (part) ± 2 (operator) ± 1 (interaction) ± 1 (repeat)
const MSA = 'part,operator,value\n' + [-1, 1].flatMap(p => [-1, 1].flatMap(o => [-1, 1].map(e => [p > 0 ? 2 : 1, o > 0 ? 2 : 1, 100 + 10 * p + 2 * o + p * o + e].join(',')))).join('\n');
// mean squares: MS_P = 800, MS_O = 32, MS_PO = 8, MS_E = 2  ⇒  σ²_repeat = 2, σ²_PO = (8 − 2)/2 = 3, σ²_oper = (32 − 8)/4 = 6, σ²_part = (800 − 8)/4 = 198
const V_GRR = 2 + 3 + 6, V_PART = 198;
const MSA_REF = 'Gauge R&R by two-way crossed random-effects ANOVA (expected mean squares) — Montgomery, Introduction to Statistical Quality Control; AIAG Measurement Systems Analysis manual (ndc = 1.41·σ_part/σ_GRR)';
// process signal: 20 baseline points alternating 9, 11 (mean 10, σ = MR̄/1.128 = 2/1.128), then a sustained shift of +2.5σ
const S0 = 2 / 1.128, SIG = spc([...Array.from({ length: 20 }, (_, i) => [i % 2 ? 11 : 9]), ...Array.from({ length: 20 }, () => [10 + 2.5 * S0])]);
const SHIFT = { mode: 'auto', spcData: SIG, nBase: 20, lamE: 0.2, LE: 3, kC: 0.5, hC: 5 };
const SHIFT_REF = 'Tabular CUSUM and EWMA control-chart recursions — Montgomery, Introduction to Statistical Quality Control';

// ---------- published measured data sets (NIST/SEMATECH e-Handbook of Statistical Methods) ----------
// Section 6.3.2.2, table "Example of moving range": flow rate of the first ten batches of a new process, transcribed in batch order
const NIST_FLOW = 'batch,flowrate\n' + [49.6, 47.6, 49.9, 51.3, 47.8, 51.2, 52.6, 52.4, 53.6, 52.1].map((v, i) => `${i + 1},${v}`).join('\n');
// Section 5.6.1.1, table "Data Used in the Analysis": eddy-current probe impedance (ohm) of a 2³ full factorial in number of wire turns (90/180),
// winding distance (0.38/1.14) and wire gauge (40/48), coded −1/+1, in the standard order of the source
const NIST_EDDY = 'turns,winding_distance,wire_gauge,impedance_ohm\n' + [[-1, -1, -1, 1.70], [1, -1, -1, 4.57], [-1, 1, -1, 0.55], [1, 1, -1, 3.39], [-1, -1, 1, 1.51], [1, -1, 1, 4.59], [-1, 1, 1, 0.67], [1, 1, 1, 4.29]].map(r => r.join(',')).join('\n');
const NIST_HB = { authors: 'National Institute of Standards and Technology; SEMATECH', journal: 'NIST/SEMATECH e-Handbook of Statistical Methods (web edition, pages undated; retrieved October 2026)', year: null, doi: null };
const NIST_FLOW_SRC = { ...NIST_HB, title: 'Individuals Control Charts (section 6.3.2.2)', url: 'https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc322.htm' };
const NIST_FLOW_COND = 'Flow rate of the first 10 batches of a new process, one measurement per batch (individuals / moving-range chart of 9 moving ranges, x̄ = 50.81, MR̄ = 1.8778)';

export default [
  // ----- reactor design -----
  { id: 's26-batch-first-order-time', suite: 26, title: 'Batch reactor, first-order kinetics: time for 90 % conversion', mode: 'reactor', inputs: { ...KIN1, rtype: 'batch' }, kpi: 't_rxn', expected: Math.log(10) / 1e-3, tol: 1e-9, kind: 'analytical', reference: LEV_REF, note: 'X = 1 − e^(−kt) ⇒ t = −ln(1 − X)/k = ln 10/10⁻³ s⁻¹ = 2302.6 s.' },
  { id: 's26-cstr-first-order-time', suite: 26, title: 'CSTR, first-order kinetics: residence time for 90 % conversion', mode: 'reactor', inputs: { ...KIN1, rtype: 'cstr' }, kpi: 't_rxn', expected: 9000, tol: 1e-9, kind: 'analytical', reference: LEV_REF, note: 'τ = X/(k·(1 − X)) = 0.9/(10⁻³·0.1).' },
  { id: 's26-cstr-first-order-volume', suite: 26, title: 'CSTR, first-order kinetics: volume for 10 t/yr of ZnO', mode: 'reactor', inputs: { ...KIN1, rtype: 'cstr' }, kpi: 'V_liq', expected: V_CSTR, tol: 1e-3, kind: 'analytical', reference: LEV_REF,
    note: 'V = F_A0·X/(k·C_A0·(1 − X)) with the feed fixed by 10 t/yr, 6000 h/yr, 92 % isolation yield and M(ZnO) = 81.379 g/mol; the tolerance covers the atomic-weight table.' },
  { id: 's26-cstr-feed-rate', suite: 26, title: 'Continuous reactor: feed rate from the product mass balance', mode: 'reactor', inputs: { ...KIN1, rtype: 'cstr' }, kpi: 'n_batch', expected: Q_FEED * 3.6e6, tol: 1e-3, kind: 'conservation', reference: LEV_REF,
    note: 'Q = ṁ_product/(M·C_A0·X·η) in L/h (1 mol ZnO per mol zinc salt); the tolerance covers the atomic-weight table.' },
  { id: 's26-pfr-second-order-time', suite: 26, title: 'Plug-flow reactor, second-order kinetics: residence time for 90 % conversion', mode: 'reactor', inputs: { ...KIN2, rtype: 'pfr' }, kpi: 't_rxn', expected: 0.9 / (1e-5 * 50 * 0.1), tol: 1e-9, kind: 'analytical', reference: LEV_REF, note: 'τ = X/(k·C_A0·(1 − X)) = 18 000 s.' },
  { id: 's26-cstr-second-order-time', suite: 26, title: 'CSTR, second-order kinetics: residence time for 90 % conversion', mode: 'reactor', inputs: { ...KIN2, rtype: 'cstr' }, kpi: 't_rxn', expected: 0.9 / (1e-5 * 50 * 0.1 ** 2), tol: 1e-9, kind: 'analytical', reference: LEV_REF, note: 'τ = X/(k·C_A0·(1 − X)²) = 180 000 s, ten times the plug-flow value.' },
  { id: 's26-batch-zero-order-time', suite: 26, title: 'Batch reactor, zero-order kinetics: time for 90 % conversion', mode: 'reactor', inputs: { ...KIN1, rtype: 'batch', order: 0, k0: 0.01 }, kpi: 't_rxn', expected: 50 * 0.9 / 0.01, tol: 1e-9, kind: 'analytical', reference: LEV_REF, note: 't = C_A0·X/k = 4500 s.' },
  { id: 's26-arrhenius-rate-constant', suite: 26, title: 'Arrhenius rate constant at 90 °C (k₀ = 1.4·10⁵ s⁻¹, Eₐ = 60 kJ/mol)', mode: 'reactor', inputs: { mode: 'reactor', order: 1, k0: 1.4e5, Ea: 60, T: 90 }, kpi: 'k', expected: 1.4e5 * Math.exp(-60e3 / (RG * 363.15)), tol: 1e-6, kind: 'analytical', reference: 'Arrhenius equation k = k₀·exp(−Eₐ/RT) with R = 8.314462618 J/(mol·K) — Atkins, Physical Chemistry' },
  { id: 's26-adiabatic-temperature-rise', suite: 26, title: 'Adiabatic temperature rise of the dilute aqueous recipe', mode: 'reactor', inputs: { mode: 'reactor', dHr: -30, c0: 50, X: 0.9, rho: 1000, cp: 4180 }, kpi: 'dT_ad', expected: 30e3 * 50 * 0.9 / (1000 * 4180), tol: 1e-9, kind: 'analytical', reference: 'Adiabatic energy balance ΔT_ad = (−ΔH_r)·C_A0·X/(ρ·cₚ) — Fogler, Elements of Chemical Reaction Engineering' },
  { id: 's26-heating-energy', suite: 26, title: 'Sensible heat per kilogram of ZnO for a 0.05 mol/L recipe heated from 20 to 90 °C', mode: 'reactor', inputs: { ...KIN1, rtype: 'batch', T0: 20, rho: 1000, cp: 4180 }, kpi: 'e_heat', expected: 1000 * 4180 * 70 / (50 * 0.9 * 0.92 * M_ZNO) / 3.6e6, tol: 1e-3, kind: 'conservation', reference: 'Sensible-heat balance Q = ρ·cₚ·ΔT per volume of liquid, divided by the product mass formed in that volume',
    note: 'About 24 kWh/kg: one kilogram of ZnO requires roughly 300 L of solution at this concentration. The tolerance covers the atomic-weight table.' },
  // ----- heat-transfer and mixing scale-up -----
  { id: 's26-area-volume-scaling', suite: 26, title: 'Geometric similarity: jacket area per volume scales as V^(−1/3)', mode: 'reactor', inputs: { ...KIN1, rtype: 'cstr', Vlab: 1 }, kpi: 'av_ratio', expected: Math.pow(V_CSTR / V_LAB, -1 / 3), tol: 1e-3, kind: 'analytical', reference: 'Surface-to-volume scaling of geometrically similar vessels, A/V ∝ V^(−1/3) — Zlokarnik, Scale-up in Chemical Engineering',
    note: 'Scale ratio from the closed-form CSTR volume (about 1240 times the 1 L laboratory reactor); the tolerance covers the atomic-weight table.' },
  { id: 's26-constant-power-per-volume', suite: 26, title: 'Constant P/V rule: power per volume at scale equals the laboratory value', mode: 'reactor', inputs: { ...MIX, mixRule: 'pv' }, kpi: 'PV', expected: PV_LAB, tol: 1e-6, kind: 'analytical', reference: MIX_REF,
    note: 'Laboratory value Po·ρ·N³·d⁵/V = 270 W/m³ for a 1 L tank (D = 112 mm, d = 45 mm, 400 rpm).' },
  { id: 's26-constant-tip-speed', suite: 26, title: 'Constant tip-speed rule: tip speed at scale equals the laboratory value', mode: 'reactor', inputs: { ...MIX, mixRule: 'tip' }, kpi: 'tip', expected: TIP_LAB, tol: 1e-6, kind: 'analytical', reference: MIX_REF, note: 'π·N·d = 0.94 m/s in the laboratory tank.' },
  { id: 's26-constant-blend-time', suite: 26, title: 'Constant mixing-time rule: turbulent blend time', mode: 'reactor', inputs: { ...MIX, mixRule: 'time' }, kpi: 'theta', expected: 5.2 / (Math.cbrt(5) * N_LAB * 0.4 ** 2), tol: 1e-6, kind: 'analytical', reference: MIX_REF, note: 'θ = 5.2/(Po^(1/3)·N·(d/T)²) = 2.85 s at 400 rpm, independent of scale when N is kept.' },
  { id: 's26-micromixing-time', suite: 26, title: 'Engulfment micromixing time at constant P/V', mode: 'reactor', inputs: { ...MIX, mixRule: 'pv' }, kpi: 'tau_E', expected: 17.24 * Math.sqrt((5e-4 / 1000) / (PV_LAB / 1000)), tol: 1e-4, kind: 'analytical', reference: 'Engulfment model of micromixing, τ_E = 17.24·√(ν/ε) — Bałdyga & Bourne, Turbulent Mixing and Chemical Reactions',
    note: 'Mean dissipation ε = P/(ρV) = 0.27 W/kg, ν = 5·10⁻⁷ m²/s ⇒ τ_E ≈ 23 ms.' },
  { id: 's26-tip-speed-at-constant-pv', suite: 26, title: 'Constant P/V rule: tip speed grows as V^(1/9)', mode: 'reactor', inputs: { ...KIN1, ...MIX, rtype: 'cstr', mixRule: 'pv' }, kpi: 'tip', expected: TIP_LAB * Math.pow(V_CSTR / V_LAB, 1 / 9), tol: 1e-3, kind: 'analytical', reference: MIX_REF,
    note: 'N ∝ d^(−2/3) at constant P/V, so π·N·d ∝ d^(1/3) ∝ V^(1/9); scale ratio from the closed-form CSTR volume.' },

  // ----- statistical process control -----
  { id: 's26-spc-sigma-within', suite: 26, title: 'X̄–R chart, subgroups of five: within-subgroup σ = R̄/d₂', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 94, target: 100 }, kpi: 'sigma_w', expected: SW5, tol: 1e-9, kind: 'reference-data', reference: SPC_REF, note: 'Every subgroup is 98…102, so R̄ = 4 and σ = 4/2.326.' },
  { id: 's26-spc-ucl', suite: 26, title: 'X̄–R chart, subgroups of five: upper control limit x̿ + A₂·R̄', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 94, target: 100 }, kpi: 'ucl', expected: 100 + A2_5 * 4, tol: 1e-4, kind: 'reference-data', reference: SPC_REF, note: 'A₂ = 0.577 is tabulated to three decimals (3/(d₂·√5) = 0.5768), hence the 0.01 % tolerance.' },
  { id: 's26-spc-cp', suite: 26, title: 'Capability index Cp = (USL − LSL)/(6σ)', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 94, target: 100 }, kpi: 'cp', expected: 12 / (6 * SW5), tol: 1e-9, kind: 'analytical', reference: SPC_REF },
  { id: 's26-spc-cpk', suite: 26, title: 'Capability index Cpk for an off-centre process', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 96, target: 100 }, kpi: 'cpk', expected: 4 / (3 * SW5), tol: 1e-9, kind: 'analytical', reference: SPC_REF, note: 'Cpk = min(USL − μ, μ − LSL)/(3σ) with μ = 100, limits 96 and 106.' },
  { id: 's26-spc-pp', suite: 26, title: 'Performance index Pp from the overall standard deviation', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 94, target: 100 }, kpi: 'pp', expected: 12 / (6 * SO5), tol: 1e-9, kind: 'analytical', reference: SPC_REF, note: 'Overall s = √(60/29) for thirty values with Σ(x − x̄)² = 60.' },
  { id: 's26-spc-cpm', suite: 26, title: 'Taguchi index Cpm with the target one unit off the mean', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 106, lsl: 94, target: 101 }, kpi: 'cpm', expected: 12 / (6 * Math.sqrt(SO5 ** 2 + 1)), tol: 1e-9, kind: 'analytical', reference: SPC_REF, note: 'Cpm = (USL − LSL)/(6·√(σ² + (μ − T)²)).' },
  { id: 's26-spc-ppm-three-sigma', suite: 26, title: 'Expected nonconforming fraction with limits at ±3σ', mode: 'spc', inputs: { mode: 'spc', spcData: SPC5, usl: 100 + 3 * SO5, lsl: 100 - 3 * SO5, target: 100 }, kpi: 'ppm', expected: 2699.8, tol: 1e-3, kind: 'reference-data', reference: 'Standard normal distribution: P(|Z| > 3) = 0.0026998 (2700 ppm) — standard normal tables, e.g. Abramowitz & Stegun, Handbook of Mathematical Functions',
    note: 'The suite evaluates Φ with a rational approximation of erf accurate to about 10⁻⁷, which is 10⁻⁴ of this tail probability.' },
  { id: 's26-spc-individuals-ucl', suite: 26, title: 'Individuals–moving-range chart: upper control limit x̄ + 2.660·MR̄', mode: 'spc', inputs: { mode: 'spc', spcData: IMR, usl: 20, lsl: 2, target: 11 }, kpi: 'ucl', expected: 11 + 2.660 * 2, tol: 1e-4, kind: 'reference-data', reference: SPC_REF, note: 'Values alternate 10, 12: mean 11, every moving range 2. E₂ = 3/1.128 is tabulated as 2.660.' },
  { id: 's26-spc-xbar-s-sigma', suite: 26, title: 'X̄–S chart, subgroups of twelve: within-subgroup σ = s̄/c₄', mode: 'spc', inputs: { mode: 'spc', spcData: spc(Array.from({ length: 4 }, () => PAT12)), usl: 115, lsl: 85, target: 100 }, kpi: 'sigma_w', expected: S12 / C4_12, tol: 1e-6, kind: 'analytical', reference: SPC_REF,
    note: 'c₄(12) = √(2/11)·Γ(6)/Γ(5.5) = 0.9776 evaluated exactly here; the tolerance covers the suite\'s log-gamma approximation.' },

  // ----- design of experiments -----
  { id: 's26-doe-main-effect', suite: 26, title: '2² factorial on an exact response: largest main effect', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2, goal: 'max' }, kpi: 'top', expected: 20, tol: 1e-8, kind: 'analytical', reference: DOE_REF, note: 'y = 50 + 10·x₁ + 5·x₂ + 2.5·x₁x₂: effect of A = ȳ₊ − ȳ₋ = 20.' },
  { id: 's26-doe-intercept', suite: 26, title: '2² factorial on an exact response: mean response', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2, goal: 'max' }, kpi: 'mean', expected: 50, tol: 1e-8, kind: 'analytical', reference: DOE_REF },
  { id: 's26-doe-optimum', suite: 26, title: '2² factorial on an exact response: predicted maximum at the (+, +) corner', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2, goal: 'max' }, kpi: 'yopt', expected: 67.5, tol: 1e-8, kind: 'analytical', reference: DOE_REF, note: 'A multilinear model takes its extremes at the vertices: 50 + 10 + 5 + 2.5.' },
  { id: 's26-doe-lenth-three-effects', suite: 26, title: "Unreplicated 2² design: Lenth margin of error", mode: 'doe', inputs: { mode: 'doe', doeData: DOE2, goal: 'max' }, kpi: 'me', expected: T975_1 * 15, tol: 1e-4, kind: 'analytical', reference: LENTH_REF,
    note: 'Effects 20, 10, 5: s₀ = 1.5·10 = 15, no effect exceeds 2.5·s₀, PSE = 15; m/3 = 1 degree of freedom, t = tan(0.475π) = 12.706.' },
  { id: 's26-doe-lenth-seven-effects', suite: 26, title: "Unreplicated 2³ design: Lenth margin of error", mode: 'doe', inputs: { mode: 'doe', doeData: DOE3, goal: 'max' }, kpi: 'me', expected: 3.764 * 4.5, tol: 1e-3, kind: 'reference-data', reference: LENTH_REF,
    note: 'Effects 20, 10, 6, 4, 2, 1, 0.5: s₀ = 1.5·4 = 6, effects below 15 have median 3, PSE = 4.5; t(0.975; 7/3) = 3.764 as tabulated by Lenth for seven contrasts.' },
  { id: 's26-doe-replicated-margin', suite: 26, title: 'Replicated 2² design: margin of error from the residual t-test', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2R, goal: 'max' }, kpi: 'me', expected: 2.776 * 1, tol: 5e-4, kind: 'reference-data', reference: DOE_REF,
    note: 'Replicates at ±1 about the model: SSE = 8 on 4 degrees of freedom, s = √2, standard error of an effect 2s/√8 = 1; t(0.975; 4) = 2.776 from standard tables (three decimals).' },
  { id: 's26-doe-residual-sd', suite: 26, title: 'Replicated 2² design: residual standard deviation', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2R, goal: 'max' }, kpi: 'sres', expected: Math.SQRT2, tol: 1e-8, kind: 'analytical', reference: DOE_REF },
  { id: 's26-doe-curvature', suite: 26, title: 'Replicated 2² design with centre points: curvature estimate', mode: 'doe', inputs: { mode: 'doe', doeData: DOE2R, goal: 'max' }, kpi: 'curv', expected: 5, tol: 1e-8, kind: 'analytical', reference: DOE_REF, note: 'Centre-point mean 55 minus factorial mean 50.' },

  // ----- yield, throughput, OEE -----
  { id: 's26-yield-poisson', suite: 26, title: 'Poisson yield at A·D₀ = 0.2', mode: 'yield', inputs: { ...YLD, ymodel: 'poisson' }, kpi: 'y_def', expected: Math.exp(-AD), tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'Y = e^(−0.2) = 0.8187.' },
  { id: 's26-yield-murphy', suite: 26, title: 'Murphy yield at A·D₀ = 0.2', mode: 'yield', inputs: { ...YLD, ymodel: 'murphy' }, kpi: 'y_def', expected: ((1 - Math.exp(-AD)) / AD) ** 2, tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'Y = [(1 − e^(−0.2))/0.2]² = 0.8215.' },
  { id: 's26-yield-seeds', suite: 26, title: 'Seeds yield at A·D₀ = 0.2', mode: 'yield', inputs: { ...YLD, ymodel: 'seeds' }, kpi: 'y_def', expected: 1 / 1.2, tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'Y = 1/(1 + 0.2) = 0.8333.' },
  { id: 's26-yield-negative-binomial', suite: 26, title: 'Negative-binomial yield at A·D₀ = 0.2, α = 2', mode: 'yield', inputs: { ...YLD, ymodel: 'negbin' }, kpi: 'y_def', expected: 1 / 1.1 ** 2, tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'Y = (1 + 0.2/2)^(−2) = 0.8264.' },
  { id: 's26-yield-density-target-poisson', suite: 26, title: 'Defect density giving 90 % Poisson yield', mode: 'yield', inputs: { ...YLD, ymodel: 'poisson' }, kpi: 'd_target', expected: -Math.log(0.9) / 0.4, tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'D = −ln 0.9/(A·f_crit) = 0.263 cm⁻² for 1 cm² with a critical-area fraction of 0.4.' },
  { id: 's26-yield-density-target-seeds', suite: 26, title: 'Defect density giving 90 % Seeds yield for a process already above 90 %', mode: 'yield', inputs: { mode: 'yield', D0: 0.001, area: 1, critFrac: 0.01, ymodel: 'seeds' }, kpi: 'd_target', expected: (1 / 0.9 - 1) / 0.01, tol: 1e-9, kind: 'analytical', reference: YLD_REF, note: 'A·D₀ = 1/0.9 − 1 = 0.111 ⇒ D = 11.1 cm⁻², independent of the present defect density.' },
  { id: 's26-rolled-throughput-yield', suite: 26, title: 'Rolled throughput yield of the seven-step patterning line', mode: 'yield', inputs: { mode: 'yield' }, kpi: 'rty', expected: 0.995 * 0.99 * 0.97 * 0.99 * 0.98 * 0.97 * 0.995, tol: 1e-9, kind: 'analytical', reference: LINE_REF, note: 'Product of the first-pass yields 99.5, 99, 97, 99, 98, 97 and 99.5 %.' },
  { id: 's26-takt-time', suite: 26, title: 'Takt time for 12 000 units in 6000 h', mode: 'yield', inputs: { mode: 'yield', units: 12000, hours: 6000 }, kpi: 'takt', expected: 30, tol: 1e-9, kind: 'analytical', reference: LINE_REF },
  { id: 's26-oee-product', suite: 26, title: 'Single-step line: OEE = availability × performance × quality', mode: 'yield', inputs: ONE, kpi: 'oee', expected: 0.90 * 0.92 * Y_ONE, tol: 1e-9, kind: 'analytical', reference: LINE_REF, note: '0.90 × 0.92 × 0.95 = 0.7866 (defect-limited yield 1 − 10⁻⁸).' },
  { id: 's26-line-capacity', suite: 26, title: 'Single-step line: good-unit capacity', mode: 'yield', inputs: ONE, kpi: 'cap', expected: 0.90 * 0.92 * 6 * Y_ONE, tol: 1e-9, kind: 'analytical', reference: LINE_REF, note: 'One machine at 10 min per unit: 6 /h × availability × performance × yield.' },
  { id: 's26-campaign-output', suite: 26, title: 'Campaign schedule: annual good output', mode: 'yield', inputs: ONE, kpi: 'camp', expected: 6000 / (10 * (6 + 1) + 12) * 10 * 25 * Y_ONE, tol: 1e-9, kind: 'analytical', reference: LINE_REF, note: 'Campaign of ten 6 h batches with 1 h cleaning each and a 12 h changeover: 82 h; 25 units per batch at 95 % yield.' },
  { id: 's26-bottleneck-utilisation', suite: 26, title: 'Two-step line: bottleneck utilisation with downstream yield loss', mode: 'yield', inputs: TWO, kpi: 'util', expected: 0.5 / Math.exp(-1e-8), tol: 1e-9, kind: 'analytical', reference: LINE_REF,
    note: 'Step 1: 6 /h × 0.9 × 0.8 = 4.32 good/h; step 2: 5 /h × 0.8 = 4.0 good/h (bottleneck); demand 2 good/h ⇒ 50 %.' },

  // ----- automation readiness -----
  { id: 's26-simc-gain', suite: 26, title: 'SIMC PI tuning: controller gain', mode: 'auto', inputs: { mode: 'auto', tauP: 600, deadT: 60, Kproc: 1.5, lam: 1 }, kpi: 'kc', expected: 600 / (1.5 * 120), tol: 1e-9, kind: 'analytical', reference: SIMC_REF, note: 'τ = 600 s, θ = 60 s, K = 1.5, τc = θ.' },
  { id: 's26-simc-integral-time', suite: 26, title: 'SIMC PI tuning: integral time for a lag-dominant process', mode: 'auto', inputs: { mode: 'auto', tauP: 600, deadT: 60, Kproc: 1.5, lam: 1 }, kpi: 'ti', expected: 480, tol: 1e-9, kind: 'analytical', reference: SIMC_REF, note: 'Tᵢ = min(600, 4·(60 + 60)) = 480 s.' },
  { id: 's26-simc-integral-time-tau', suite: 26, title: 'SIMC PI tuning: integral time equals τ when τ < 4·(τc + θ)', mode: 'auto', inputs: { mode: 'auto', tauP: 100, deadT: 60, Kproc: 1.5, lam: 1 }, kpi: 'ti', expected: 100, tol: 1e-9, kind: 'analytical', reference: SIMC_REF },
  { id: 's26-simc-phase-margin', suite: 26, title: 'SIMC PI tuning with τc = θ: phase margin 61.4°', mode: 'auto', inputs: { mode: 'auto', tauP: 100, deadT: 60, Kproc: 1.5, lam: 1, Ts: TS_PM }, kpi: 'pm', expected: PM_REF, tol: 1e-5, kind: 'analytical', reference: SIMC_REF,
    note: 'With Tᵢ = τ the loop is e^(−θs)/((τc + θ)·s): ω_c = 1/(2θ) and PM = 90° − 0.5 rad = 61.35°, reduced by ω_c·Tₛ/2 = 0.02° for the 0.1 s sample-and-hold.' },
  { id: 's26-nyquist-interval', suite: 26, title: 'Nyquist sampling interval for a 20 s time constant', mode: 'auto', inputs: { mode: 'auto', tauFast: 20 }, kpi: 'ts_nyq', expected: Math.PI * 20, tol: 1e-9, kind: 'analytical', reference: 'Sampling theorem f_s > 2·f_max with the first-order corner frequency f_max = 1/(2πτ) — Åström & Wittenmark, Computer-Controlled Systems' },
  { id: 's26-gauge-rr', suite: 26, title: 'Crossed gauge study with exact variance components: %GRR', mode: 'auto', inputs: { mode: 'auto', msa: MSA }, kpi: 'grr', expected: 100 * Math.sqrt(V_GRR / (V_GRR + V_PART)), tol: 1e-8, kind: 'analytical', reference: MSA_REF,
    note: 'Expected mean squares give σ²_repeat = 2, σ²_interaction = 3, σ²_operator = 6, σ²_part = 198: %GRR = √(11/209) = 22.9 %.' },
  { id: 's26-gauge-ndc', suite: 26, title: 'Crossed gauge study with exact variance components: number of distinct categories', mode: 'auto', inputs: { mode: 'auto', msa: MSA }, kpi: 'ndc', expected: 5, abs: 0, kind: 'analytical', reference: MSA_REF, note: 'ndc = ⌊1.41·√(198/11)⌋ = ⌊5.98⌋.' },
  { id: 's26-cusum-alarm', suite: 26, title: 'Tabular CUSUM (k = 0.5, h = 5): alarm three observations after a 2.5σ shift', mode: 'auto', inputs: SHIFT, kpi: 'cusum', expected: 23, abs: 0, kind: 'analytical', reference: SHIFT_REF,
    note: 'C⁺ is 0.064 at the end of the 20-point baseline and then grows by 2.5 − 0.5 = 2 per observation: 2.06, 4.06, 6.06 > 5 at observation 23.' },
  { id: 's26-ewma-alarm', suite: 26, title: 'EWMA (λ = 0.2, L = 3): alarm three observations after a 2.5σ shift', mode: 'auto', inputs: SHIFT, kpi: 'ewma', expected: 23, abs: 0, kind: 'analytical', reference: SHIFT_REF,
    note: 'In units of σ the EWMA deviation after the shift is 0.55, 0.94, 1.25 against the asymptotic limit L·√(λ/(2 − λ)) = 1.0, first exceeded at observation 23.' },
  { id: 's26-cusum-arl0', suite: 26, title: 'Two-sided tabular CUSUM (k = 0.5, h = 5): in-control average run length', mode: 'auto', inputs: { mode: 'auto', kC: 0.5, hC: 5 }, kpi: 'arl0', expected: 465, tol: 0.02, kind: 'numerical-reference', reference: 'Tabulated average run lengths of the tabular CUSUM (ARL₀ = 465 for k = ½, h = 5) — Montgomery, Introduction to Statistical Quality Control; Hawkins & Olwell, Cumulative Sum Charts and Charting for Quality Improvement',
    note: "The suite uses Siegmund's closed-form approximation (h replaced by h + 1.166), which is accurate to about 1–2 % against the tabulated Markov-chain value." },
  // ----- published experiments -----
  { id: 's26-exp-nist-flowrate-imr-ucl', suite: 26, title: 'NIST flow-rate batches: upper control limit of the individuals chart', mode: 'spc', inputs: { mode: 'spc', spcData: NIST_FLOW, usl: 56, lsl: 45, target: 50 }, kpi: 'ucl', expected: 55.8041, tol: 1e-4, kind: 'experiment',
    measured: { value: 55.8041, uncertainty: null, unit: 'flow-rate units of the source' }, source: NIST_FLOW_SRC, conditions: NIST_FLOW_COND, match: ['spcData'],
    reference: 'NIST/SEMATECH e-Handbook of Statistical Methods, section 6.3.2.2 (Individuals Control Charts)',
    note: 'Measured data and reduced value from the section "Example of moving range" / "Limits for the moving range chart": UCL = 50.81 + 3·1.8778/1.128 = 55.8041. The ten measurements are transcribed from the table of the source; the suite recomputes the mean, the mean moving range and the limit from them. The source rounds the mean moving range to 1.8778, hence the tolerance of 0.01 %. The specification limits and the target entered here are not given by the source and do not enter the control limits.' },
  { id: 's26-exp-nist-flowrate-imr-lcl', suite: 26, title: 'NIST flow-rate batches: lower control limit of the individuals chart', mode: 'spc', inputs: { mode: 'spc', spcData: NIST_FLOW, usl: 56, lsl: 45, target: 50 }, kpi: 'lcl', expected: 45.8159, tol: 1e-4, kind: 'experiment',
    measured: { value: 45.8159, uncertainty: null, unit: 'flow-rate units of the source' }, source: NIST_FLOW_SRC, conditions: NIST_FLOW_COND, match: ['spcData'],
    reference: 'NIST/SEMATECH e-Handbook of Statistical Methods, section 6.3.2.2 (Individuals Control Charts)',
    note: 'Same data set; LCL = 50.81 − 3·1.8778/1.128 = 45.8159 as stated under "Limits for the moving range chart". Tolerance 0.01 % for the rounding of the published intermediate values. The specification limits and the target entered here are not given by the source and do not enter the control limits.' },
  { id: 's26-exp-nist-eddy-current-turns-effect', suite: 26, title: 'NIST eddy-current probe 2³ factorial: effect of the number of wire turns on the measured impedance', mode: 'doe', inputs: { mode: 'doe', doeData: NIST_EDDY, goal: 'max' }, kpi: 'top', expected: 3.1025, tol: 1e-6, kind: 'experiment',
    measured: { value: 3.1025, uncertainty: null, unit: 'ohm' },
    source: { ...NIST_HB, title: 'Eddy Current Probe Sensitivity Case Study (sections 5.6.1.1 Background and Data, 5.6.1.3 Interaction Effects, 5.6.1.5 Estimate Main and Interaction Effects)', url: 'https://www.itl.nist.gov/div898/handbook/pri/section6/pri611.htm' },
    conditions: 'Eddy-current crack-detector probe (study by Capobianco, Splett and Iyer, NIST): 2³ full factorial, eight runs in random order, number of wire turns 90/180, winding distance 0.38/1.14, wire gauge 40/48; response: probe impedance (sensitivity) in ohm',
    match: ['doeData'],
    reference: 'NIST/SEMATECH e-Handbook of Statistical Methods, section 5.6.1 (Eddy Current Probe Sensitivity Case Study)',
    note: 'The eight measured impedances are transcribed from the table in section 5.6.1.1. Section 5.6.1.5 tabulates the least-squares coefficient of X1 (number of turns) as 1.55125 and states that the effects are twice the coefficients; section 5.6.1.3 names X1 the most important effect with magnitude 3.1025 (printed there with a minus sign, although the tabulated coefficient and the data are positive). The suite identifies the same factor as the largest effect and computes it from the data; agreement is exact to rounding because the design is orthogonal. The grand mean of the source (2.65875) is reproduced as well.' },
];
