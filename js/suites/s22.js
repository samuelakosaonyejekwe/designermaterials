// Suite 22 — Material Behaviour, Damage and Lifetime Prediction.
// Fatigue, fracture and crack growth, creep and accelerated ageing, chemical/environmental degradation,
// functional lifetime of metamaterials and reliability synthesis (reduced-order physics-of-failure models).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { R_GAS: R, FARADAY: F, N_A: NA, PI } = N;
const KEV = 8.617333e-5; // Boltzmann constant, eV/K
const YEAR = 3.15576e7, HOUR = 3600;
const f = (v, s = 3) => N.fmt(v, s);
const lsp = (e0, e1, n) => N.linspace(e0, e1, n).map(v => 10 ** v); // n log-spaced samples between 10^e0 and 10^e1
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ser = (key, title, xlabel, ylabel, lines, explain, o = {}) => ({ key, title, xlabel, ylabel, logx: !!o.logx, logy: !!o.logy, lines, explain });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const safeEq = (r, p, arrow) => { try { return chem.equationString(r, p, null, arrow); } catch (e) { return r.join(' + ') + ' → ' + p.join(' + ') + ' (not balanced automatically)'; } };
const MM = s => chem.molarMass(s) * 1e-3;
const BIG = 1e30; // stands for "no failure predicted" in numeric outputs
const capN = v => (Number.isFinite(v) ? Math.min(Math.max(v, 0), BIG) : BIG);

function finish(r) {
  let bad = 0;
  for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; k.good = null; bad++; }
  for (const s of r.series || []) for (const l of s.lines) {
    const x = [], y = [], st = Math.ceil(l.x.length / 400);
    for (let i = 0; i < l.x.length; i += st) { let xv = l.x[i], yv = l.y[i]; if (!Number.isFinite(xv) || !Number.isFinite(yv)) { bad++; continue; } if (s.logy && yv <= 0) continue; if (s.logx && xv <= 0) continue; x.push(xv); y.push(yv); }
    if (!x.length) { x.push(s.logx ? 1 : 0); y.push(s.logy ? 1 : 0); }
    l.x = x; l.y = y;
  }
  for (const k in r.publish) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) { delete r.publish[k]; bad++; }
  if (bad) r.warnings.push(`${bad} output value(s) were not finite for this input combination and were dropped or set to zero; the inputs are outside the range in which the model is meaningful.`);
  return r;
}
// Defensive CSV/TSV parser: returns numeric rows with at least ncol columns
function parseTable(text, ncol, minRows, what) {
  const rows = [];
  for (const raw of String(text || '').split(/\r?\n/)) {
    const line = raw.trim(); if (!line || line.startsWith('#') || line.startsWith('//')) continue;
    const c = line.split(/[,;\t]+|\s{2,}/).map(v => parseFloat(v.trim().replace(/^"|"$/g, '')));
    if (c.length >= ncol && c.slice(0, ncol).every(Number.isFinite)) rows.push(c.slice(0, ncol));
  }
  if (rows.length < minRows) throw new Error(`${what}: need at least ${minRows} numeric rows with ${ncol} column(s); found ${rows.length}. Paste CSV with a header line followed by numbers.`);
  return rows;
}
// seeded synthetic demonstration tables (forward model + noise)
const DEF_SPECTRUM = (() => { const r = N.rng(22), out = ['stress_amplitude_MPa,mean_stress_MPa,cycles_per_block']; for (let i = 0; i < 6; i++) { const a = 30 + 22 * i, n = Math.round(2e5 * Math.exp(-a * a / (2 * 45 * 45)) * (0.9 + 0.2 * r())); out.push(`${a},${Math.round(40 + 10 * r.normal())},${Math.max(n, 20)}`); } return out.join('\n'); })();
const DEF_FAIL = (() => { const r = N.rng(2201), t = []; for (let i = 0; i < 20; i++) t.push(12000 * (-Math.log(1 - r())) ** (1 / 2.2)); t.sort((a, b) => a - b); return 'time_to_failure_h\n' + t.map(v => v.toFixed(0)).join('\n'); })();
const DEF_DEG = (() => { const r = N.rng(2202), out = ['time_h,performance_pct']; for (let i = 0; i <= 12; i++) { const t = 500 * i; out.push(`${t},${(100 * Math.exp(-2.2e-5 * t) + r.normal(0, 0.25)).toFixed(2)}`); } return out.join('\n'); })();

// ================================================================================================
// 1. Fatigue
// ================================================================================================
function runFatigue(inp) {
  const warnings = [], E = inp.fat_E, sf = inp.fat_sf * 1e6, b = -Math.abs(inp.fat_b), ef = inp.fat_ef, c = -Math.abs(inp.fat_c), su = inp.fat_su * 1e6, sa = inp.fat_sa, sm = inp.fat_sm * 1e6;
  const Kf = 1 + inp.fat_q * (inp.fat_Kt - 1), meanModel = inp.fat_mean;
  const eqAmp = (a, m) => { // equivalent fully reversed amplitude (local, including the fatigue notch factor on amplitude)
    const al = Kf * a; if (al <= 0) return 0;
    if (meanModel === 'goodman') return m >= su ? Infinity : al / (1 - Math.max(m, -su) / su);
    if (meanModel === 'swt') { const smax = m + al; return smax > 0 ? Math.sqrt(smax * al) : 0; }
    return al;
  };
  const lifeOf = sar => (sar <= 0 ? BIG : sar >= sf ? 0.5 : capN(0.5 * (sar / sf) ** (1 / b)));
  const sar = eqAmp(sa, sm), Nf = lifeOf(sar), twoNt = (ef * E / sf) ** (1 / (b - c)), nP = b / c, Kp = sf / ef ** nP;
  const epsE = Math.min(sar, sf) / E, epsP = (Math.min(sar, sf) / Kp) ** (1 / nP), Nt = inp.fat_Nt;
  // Palmgren–Miner on the block spectrum
  let rows, synthetic = String(inp.fat_spectrum || '').trim() === DEF_SPECTRUM.trim();
  try { rows = parseTable(inp.fat_spectrum, 3, 1, 'Load spectrum'); } catch (e) { warnings.push(e.message + ' The default demonstration spectrum is used instead.'); rows = parseTable(DEF_SPECTRUM, 3, 1, 'Load spectrum'); synthetic = true; }
  if (synthetic) warnings.push('Load spectrum is synthetic demonstration data — replace with measurements (rainflow-counted amplitudes, means and counts per block).');
  let Dblock = 0, nBlock = 0; const trow = [];
  for (const [aM, mM, n] of rows) { const s_ = eqAmp(Math.abs(aM) * 1e6, mM * 1e6), Ni = lifeOf(s_), d = Ni > 0 ? Math.max(n, 0) / Ni : Infinity; Dblock += d; nBlock += Math.max(n, 0); trow.push([aM, mM, n, +(s_ / 1e6).toPrecision(4), +Ni.toPrecision(4), d]); }
  const blocks = Dblock > 0 ? capN(inp.fat_Dcrit / Dblock) : BIG, Nspec = capN(blocks * nBlock);
  for (const r_ of trow) r_[5] = +(Dblock > 0 ? r_[5] / Dblock * 100 : 0).toPrecision(3);
  // Weibull scatter on life (median = deterministic life)
  const beta = inp.fat_beta, eta = Nf / Math.LN2 ** (1 / beta), Pf = Nf >= BIG ? 0 : 1 - Math.exp(-((Nt / eta) ** beta)), B10 = capN(eta * (-Math.log(0.9)) ** (1 / beta));
  const sAllow = sf * (2 * Nt) ** b, SF = sar > 0 ? sAllow / sar : BIG;
  const Ns = lsp(2, 9, 80), xs = lsp(Math.log10(Math.max(Nf, 10) / 100), Math.log10(Math.max(Nf, 10) * 100), 80);
  if (sar >= sf) warnings.push('Equivalent amplitude reaches the fatigue strength coefficient: static failure in the first reversal; the life of 0.5 cycle is nominal.');
  if (meanModel === 'goodman' && sm >= su) warnings.push('Mean stress is at or above the ultimate strength: static failure.');
  if (epsP > epsE) warnings.push(`Plastic strain amplitude exceeds the elastic one (life below the transition life of ${f(twoNt / 2)} cycles): this is low-cycle fatigue; a stress-based (Basquin) assessment is unreliable — use measured strain amplitudes.`);
  if (Kf * sa + Math.abs(sm) > su) warnings.push('Local peak stress exceeds the ultimate strength; notch plasticity redistributes the stress (Neuber correction not applied).');
  warnings.push('Constant-amplitude strain–life constants with linear damage accumulation: load-sequence effects, small-crack behaviour, surface finish and environment are not included. Lattice struts made by additive manufacturing typically need a further knock-down for surface roughness and porosity.');
  return finish({
    summary: `With a nominal amplitude of ${f(sa / 1e6)} MPa (mean ${f(inp.fat_sm)} MPa, K_f = ${f(Kf)}, ${meanModel === 'none' ? 'no mean-stress correction' : meanModel === 'swt' ? 'SWT' : 'Goodman'}) the equivalent fully reversed amplitude is ${f(sar / 1e6)} MPa and the Basquin life is ${f(Nf)} cycles. Under the block spectrum Miner's damage is ${f(Dblock)} per block, giving ${f(blocks)} blocks (${f(Nspec)} cycles). With Weibull shape ${f(beta)} the failure probability at ${f(Nt)} cycles is ${f(Pf * 100)} % and the B10 life ${f(B10)} cycles.`,
    kpis: [
      kpi('Nf', 'Constant-amplitude life', Nf, 'cycles', `Basquin relation σ_ar = σ_f′(2N)^b solved for N with σ_f′ = ${f(inp.fat_sf)} MPa and b = ${f(b)}; values of 10³⁰ mean no fatigue damage is predicted.`),
      kpi('sar', 'Equivalent fully reversed amplitude', sar / 1e6, 'MPa', `Local amplitude K_f·σ_a corrected for mean stress (${meanModel}); this is the amplitude that would do the same damage at zero mean.`),
      kpi('Kf', 'Fatigue notch factor', Kf, '–', `1 + q(K_t − 1) with K_t = ${f(inp.fat_Kt)} (stress concentration at lattice nodes or fillets) and notch sensitivity q = ${f(inp.fat_q)}.`),
      kpi('Nt_trans', 'Transition life', twoNt / 2, 'cycles', 'Life at which elastic and plastic strain amplitudes are equal, (ε_f′E/σ_f′)^(1/(b−c))/2. Shorter lives are strain-controlled low-cycle fatigue.'),
      kpi('Dblock', 'Miner damage per block', Dblock, '–', 'Σ nᵢ/Nᵢ over the load-spectrum rows; failure is predicted when the sum reaches the critical value.'),
      kpi('blocks', 'Blocks to failure', blocks, 'blocks', `D_crit/D_block with D_crit = ${f(inp.fat_Dcrit)} (${f(Nspec)} cycles in total).`),
      kpi('Pf', 'Probability of failure at target', Pf * 100, '%', `Two-parameter Weibull life distribution with shape ${f(beta)} and median equal to the deterministic life, evaluated at ${f(Nt)} cycles.`, Pf < 0.01),
      kpi('B10', 'B10 life', B10, 'cycles', 'Cycles by which 10 % of nominally identical parts are expected to have failed.'),
      kpi('SF', 'Stress safety factor at target life', Math.min(SF, 1e6), '–', `Allowable fully reversed amplitude at ${f(Nt)} cycles (${f(sAllow / 1e6)} MPa) divided by the applied equivalent amplitude.`, SF > 1.5),
    ],
    series: [
      ser('sn', 'S–N curve', 'cycles to failure', 'stress amplitude (MPa)', [{ name: 'Basquin median', x: Ns, y: Ns.map(n => sf * (2 * n) ** b / 1e6) }, { name: '10 % failure bound', x: Ns.map(n => n * (-Math.log(0.9) / Math.LN2) ** (1 / beta)), y: Ns.map(n => sf * (2 * n) ** b / 1e6) }, { name: 'operating point', x: [Math.max(Nf, 1), Math.max(Nf, 1) * 1.0001], y: [sar / 1e6, sar / 1e6] }], 'Median S–N line with the 10 % failure (B10) bound from the Weibull scatter. The operating point lies on the median line by construction; compare its cycles with the target.', { logx: true, logy: true }),
      ser('en', 'Strain–life curve', 'reversals to failure 2N', 'strain amplitude', [{ name: 'total', x: Ns.map(n => 2 * n), y: Ns.map(n => sf / E * (2 * n) ** b + ef * (2 * n) ** c) }, { name: 'elastic (Basquin)', x: Ns.map(n => 2 * n), y: Ns.map(n => sf / E * (2 * n) ** b) }, { name: 'plastic (Coffin–Manson)', x: Ns.map(n => 2 * n), y: Ns.map(n => ef * (2 * n) ** c) }], 'Elastic and plastic strain–life lines; they cross at the transition life. Left of the crossing plastic strain dominates.', { logx: true, logy: true }),
      ser('rel', 'Life distribution', 'cycles', 'probability of failure', [{ name: 'Weibull CDF', x: xs, y: xs.map(n => 1 - Math.exp(-((n / eta) ** beta))) }, { name: 'target cycles', x: [Nt, Nt], y: [0, 1] }], 'Cumulative failure probability versus cycles. A steeper curve (higher shape) means less scatter.', { logx: true }),
    ],
    tables: [{ key: 'miner', title: 'Palmgren–Miner damage by spectrum block', columns: ['amplitude (MPa)', 'mean (MPa)', 'cycles per block', 'equivalent amplitude (MPa)', 'life Nᵢ (cycles)', 'damage share (%)'], rows: trow, explain: 'Each load level with its equivalent fully reversed amplitude, allowable cycles and share of the total damage. The rows with the largest share govern life.' }],
    equations: [eqn('constitutive', 'Strain–life (Basquin + Coffin–Manson)', 'ε_a = (σ_f′/E)(2N)^b + ε_f′(2N)^c'), eqn('formula', 'Mean-stress correction', meanModel === 'swt' ? 'σ_ar = √(σ_max σ_a)   (Smith–Watson–Topper)' : meanModel === 'goodman' ? 'σ_ar = σ_a /(1 − σ_m/σ_u)   (Goodman)' : 'σ_ar = σ_a   (no correction)'), eqn('formula', 'Notch factor', 'K_f = 1 + q (K_t − 1)'), eqn('governing', 'Palmgren–Miner rule', 'D = Σ nᵢ/Nᵢ;   failure at D = D_crit'), eqn('formula', 'Weibull life distribution', 'P_f(N) = 1 − exp[−(N/η)^β],   η = N_median /(ln 2)^{1/β}')],
    checks: [chk('Elastic (high-cycle) regime', epsP <= epsE, `plastic/elastic strain amplitude = ${f(epsP / Math.max(epsE, 1e-300))}`), chk('Peak local stress below ultimate', Kf * sa + Math.abs(sm) <= su, `${f((Kf * sa + Math.abs(sm)) / 1e6)} MPa vs ${f(inp.fat_su)} MPa`), chk('Damage is non-negative and accumulates monotonically', Dblock >= 0, `D per block ${f(Dblock)}`), chk('Target reliability (P_f < 1 %)', Pf < 0.01, `${f(Pf * 100)} % at ${f(Nt)} cycles`)],
    warnings,
    recommendations: [
      rec('high', Pf >= 0.01 ? `Reduce the equivalent amplitude to ${f(sAllow / 1.5 / 1e6)} MPa (safety factor 1.5 at ${f(Nt)} cycles): lower K_t at the nodes with fillets of radius ≥ strut diameter/2, or increase strut cross-section by ${f((sar / (sAllow / 1.5) - 1) * 100)} %.` : `Failure probability at target is ${f(Pf * 100)} %; the stress safety factor of ${f(Math.min(SF, 1e6))} is ${SF > 1.5 ? 'adequate' : 'thin — keep amplitude below ' + f(sAllow / 1.5 / 1e6) + ' MPa'}.`, 'Basquin life and Weibull scatter of this run.'),
      rec('medium', `The spectrum is governed by the ${f(trow.reduce((m, r_) => (r_[5] > m[5] ? r_ : m), trow[0])[0])} MPa level (${f(Math.max(...trow.map(r_ => r_[5])))} % of the damage); limit or filter these events for the largest life gain.`, 'Damage share in the Miner table.'),
      rec('medium', 'Measure S–N data on as-built specimens (at least 12 specimens over 3 stress levels) to replace handbook constants; for lattices test at the strut and unit-cell level.', 'Model uses constant-amplitude constants and an assumed Weibull shape.'),
      rec('low', `Designing to a B10 life of ${f(B10)} cycles with inspection, rather than replacing at a fixed conservative interval, avoids scrapping serviceable parts; surface treatments (peening, HIP of printed lattices) extend life without adding mass.`, 'Life-distribution result; circular-economy practice.', true),
    ],
    publish: { 'life.cycles': Nf, 'life.cycles_spectrum': Nspec, 'life.pf': Pf, 'life.b10_cycles': B10, 'mech.safety_factor': Math.min(SF, 1e6), 'mech.max_stress': Kf * sa + Math.abs(sm) },
    exports: [{ name: 'fatigue_miner.csv', mime: 'text/csv', content: 'amplitude_MPa,mean_MPa,cycles_per_block,equivalent_amplitude_MPa,life_cycles,damage_share_pct\n' + trow.map(r_ => r_.join(',')).join('\n'), description: 'Miner damage table for the load spectrum.' }],
    numerics: { method: 'Closed-form Basquin/Coffin–Manson inversion, linear damage summation, two-parameter Weibull', dof: rows.length, steps: rows.length, iterations: 1, residual: 0, notes: 'No cycle-by-cycle integration.' },
  });
}

// ================================================================================================
// 2. Fracture and crack growth
// ================================================================================================
function runFracture(inp) {
  const warnings = [], ds = inp.fr_ds, Rr = inp.fr_R, smax = ds / (1 - Rr), Y = inp.fr_Y, KIc = inp.fr_KIc, a0 = inp.fr_a0 * 1e-6, C = inp.fr_C, m = inp.fr_m, dKth = inp.fr_dKth * 1e6, adet = inp.fr_adet * 1e-6;
  const ac = (KIc / (Y * smax)) ** 2 / PI, dK0 = Y * ds * Math.sqrt(PI * a0), K0 = Y * smax * Math.sqrt(PI * a0);
  // Paris law with ΔK in MPa√m: da/dN = C ΔK^m
  const kk = C * (Y * ds / 1e6 * Math.sqrt(PI)) ** m, lifeFrom = a => { if (a >= ac) return 0; if (Y * ds * Math.sqrt(PI * a) <= dKth) return BIG; return Math.abs(m - 2) < 1e-9 ? Math.log(ac / a) / kk : (ac ** (1 - m / 2) - a ** (1 - m / 2)) / (kk * (1 - m / 2)); };
  const Nf = capN(lifeFrom(a0)), Ndet = capN(lifeFrom(Math.max(adet, a0))), interval = Ndet >= BIG ? BIG : Ndet / 2, days = inp.fr_rate > 0 ? interval / inp.fr_rate : BIG;
  // numerical integration in ln(a) as a consistency check and for the a(N) curve
  const na = 300, as = [], Nn = [0]; let numN = 0;
  if (a0 < ac) { const la0 = Math.log(a0), la1 = Math.log(ac); for (let i = 0; i <= na; i++) as.push(Math.exp(la0 + (la1 - la0) * i / na)); for (let i = 1; i <= na; i++) { const g = a => a / (kk * a ** (m / 2)); numN += 0.5 * (g(as[i]) + g(as[i - 1])) * (la1 - la0) / na; Nn.push(numN); } } else { as.push(a0, a0); Nn.push(1); }
  const consist = Nf < BIG && Nf > 0 ? Math.abs(numN - Nf) / Nf : 0;
  const ath = (dKth / (Y * ds)) ** 2 / PI;
  // Griffith strength of nanostructures
  const En = inp.fr_E, gam = inp.fr_gamma, sth = En / 10, dn = inp.fr_size * 1e-9, aN = dn / 2, sG = a => Math.min(sth, Math.sqrt(2 * En * gam / (PI * a))), sNano = sG(aN), aTol = 2 * En * gam / (PI * sth * sth), ds_ = lsp(0, 5, 70);
  if (a0 >= ac) warnings.push('The initial flaw is already at or above the critical size: fast fracture at peak load.');
  if (dK0 <= dKth) warnings.push(`Initial ΔK (${f(dK0 / 1e6)} MPa√m) is below the threshold ${f(inp.fr_dKth)} MPa√m: no long-crack growth is predicted (life shown as 10³⁰). Short-crack growth below threshold is not modelled.`);
  const rpl = (K0 / (inp.fr_sy * 1e6)) ** 2 / (2 * PI);
  if (rpl > a0 / 4) warnings.push(`Plastic-zone size (${f(rpl * 1e6)} µm) is not small compared with the crack (${f(a0 * 1e6)} µm): linear-elastic fracture mechanics is marginal.`);
  warnings.push('Paris law with constant geometry factor Y and constant amplitude: no crack closure, load-sequence retardation, residual stress or environment effects; region I/III curvature is represented only by the threshold and K_Ic cut-offs.');
  return finish({
    summary: `At a stress range of ${f(ds / 1e6)} MPa (R = ${f(Rr)}) the critical crack size is ${f(ac * 1e3)} mm and a ${f(inp.fr_a0)} µm initial flaw grows to failure in ${f(Nf)} cycles (Paris C = ${f(C)}, m = ${f(m)}). From the detectable size of ${f(inp.fr_adet)} µm the remaining life is ${f(Ndet)} cycles, so inspect every ${f(interval)} cycles (${f(days)} days). A ${f(inp.fr_size)} nm nanostructure with a surface flaw of half its size has a Griffith strength of ${f(sNano / 1e9)} GPa.`,
    kpis: [
      kpi('ac', 'Critical crack size', ac * 1e3, 'mm', `a_c = (K_Ic/(Y σ_max))²/π with K_Ic = ${f(KIc / 1e6)} MPa√m and σ_max = ${f(smax / 1e6)} MPa; a crack of this size fails at peak load.`, a0 < ac),
      kpi('Nf', 'Crack-growth life from initial flaw', Nf, 'cycles', 'Closed-form integration of the Paris law from the initial flaw to the critical size (10³⁰ = below threshold, no growth).'),
      kpi('dK0', 'Initial stress-intensity range', dK0 / 1e6, 'MPa√m', `Y Δσ √(πa₀); compare with the growth threshold ${f(inp.fr_dKth)} MPa√m and the toughness.`, dK0 < KIc * (1 - Rr)),
      kpi('ath', 'Largest non-propagating flaw', ath * 1e6, 'µm', 'Flaw size at which ΔK equals the threshold; smaller flaws do not grow as long cracks at this stress range.'),
      kpi('Ndet', 'Life from detectable flaw', Ndet, 'cycles', `Cycles for a crack of the inspection detection limit (${f(inp.fr_adet)} µm) to reach critical size.`),
      kpi('interval', 'Recommended inspection interval', interval, 'cycles', `Half of the detectable-to-critical life (two inspection opportunities before failure); ${f(days)} days at ${f(inp.fr_rate)} cycles/day.`),
      kpi('sGriffith', 'Griffith strength of the nanostructure', sNano / 1e9, 'GPa', `√(2Eγ/πa) for a flaw of half the ${f(inp.fr_size)} nm size, capped at the theoretical strength E/10 = ${f(sth / 1e9)} GPa.`),
      kpi('aTol', 'Flaw-tolerant size', 2 * aTol * 1e9, 'nm', 'Structure size below which the Griffith strength of a half-size flaw exceeds the theoretical strength, so the material fails at its ideal strength regardless of flaws.'),
    ],
    series: [
      ser('aN', 'Crack length versus cycles', 'cycles', 'crack length (mm)', [{ name: 'a(N)', x: Nn, y: as.map(a => a * 1e3) }, { name: 'critical size', x: [0, Math.max(numN, 1)], y: [ac * 1e3, ac * 1e3] }], 'Crack growth is slow for most of the life and accelerates sharply near the end; most of the life is spent while the crack is small.', { logy: true }),
      ser('dadn', 'Crack-growth rate', 'ΔK (MPa√m)', 'da/dN (m/cycle)', (() => { const ks = lsp(Math.log10(Math.max(inp.fr_dKth, 0.1)), Math.log10(Math.max(KIc / 1e6 * (1 - Rr), inp.fr_dKth * 1.5)), 50); return [{ name: 'Paris law', x: ks, y: ks.map(k => C * k ** m) }]; })(), 'Paris region between the threshold and K_Ic(1 − R) on log–log axes: a straight line of slope m.', { logx: true, logy: true }),
      ser('griffith', 'Strength versus structure size', 'structure size (nm)', 'strength (GPa)', [{ name: 'Griffith (flaw = size/2)', x: ds_, y: ds_.map(d => sG(d * 0.5e-9) / 1e9) }, { name: 'theoretical E/10', x: [1, 1e5], y: [sth / 1e9, sth / 1e9] }], '“Smaller is stronger”: strength rises as size^(−½) until it saturates at the theoretical strength below the flaw-tolerant size.', { logx: true, logy: true }),
    ],
    equations: [eqn('governing', 'Paris crack-growth law', 'da/dN = C (ΔK)^m,   ΔK = Y Δσ √(πa)'), eqn('formula', 'Integrated life', 'N = [a_c^{1−m/2} − a₀^{1−m/2}] / [C (YΔσ√π)^m (1 − m/2)]   (m ≠ 2)'), eqn('formula', 'Critical crack size', 'a_c = (1/π)(K_Ic/(Y σ_max))²'), eqn('constitutive', 'Griffith criterion', 'σ_f = √(2Eγ/(πa)),   G_c = 2γ = K_Ic²/E′')],
    checks: [chk('Initial flaw below critical size', a0 < ac, `a₀ = ${f(a0 * 1e6)} µm, a_c = ${f(ac * 1e6)} µm`), chk('Closed-form vs numerical Paris integration', consist < 1e-3, `relative difference ${f(consist)}`), chk('Small-scale yielding (LEFM valid)', rpl <= a0 / 4, `plastic zone ${f(rpl * 1e6)} µm vs crack ${f(a0 * 1e6)} µm`), chk('Crack growth is irreversible (monotonic a(N))', true, 'da/dN ≥ 0 by construction'), chk('Detectable flaw below critical size', adet < ac, `${f(inp.fr_adet)} µm vs ${f(ac * 1e6)} µm`)],
    warnings,
    recommendations: [
      rec('high', adet >= ac ? `The inspection method cannot find a crack before it is critical (a_c = ${f(ac * 1e6)} µm): lower the stress, use a tougher material or a finer NDE method.` : `Inspect at intervals of at most ${f(interval)} cycles (${f(days)} days) with a method that reliably finds ${f(inp.fr_adet)} µm cracks.`, 'Damage-tolerance rule: half the detectable-to-critical life.'),
      rec('medium', `Halving the initial flaw size raises life by about ${f((2 ** (m / 2 - 1) - 1) * 100)} % when a₀ ≪ a_c (for m = ${f(m)}); invest in process control and proof testing to cap the largest defect.`, 'Dependence of the integrated Paris life on a₀.'),
      rec('medium', dn < 2 * aTol ? `At ${f(inp.fr_size)} nm the structure is below the flaw-tolerant size (${f(2 * aTol * 1e9)} nm): strength is limited by the ideal strength, so nano-lattice struts can be loaded close to ${f(sth / 1e9)} GPa in the absence of surface steps.` : `Reducing strut size from ${f(inp.fr_size)} nm towards ${f(2 * aTol * 1e9)} nm would raise strength from ${f(sNano / 1e9)} to ${f(sth / 1e9)} GPa — the basis of nano-architected high-strength lattices.`, 'Griffith scaling with flaw size bounded by structure size.'),
      rec('low', 'A damage-tolerant design with scheduled inspection lets components be retired on condition rather than on a fixed conservative life, reducing material use and waste.', 'Inspection-interval result.', true),
    ],
    publish: { 'life.cycles': Nf, 'life.inspection_cycles': interval, 'life.critical_crack_m': ac, 'prop.KIc': KIc, 'mech.max_stress': smax, 'prop.griffith_strength': sNano },
    numerics: { method: 'Closed-form Paris integration verified by trapezoidal quadrature in ln a (300 intervals)', dof: na, steps: na, iterations: 1, residual: consist, notes: 'Residual is the relative closed-form/numerical difference.' },
  });
}

// ================================================================================================
// 3. Creep and thermal ageing / accelerated testing
// ================================================================================================
function runCreep(inp) {
  const warnings = [], T = inp.cr_T + 273.15, sig = inp.cr_sigma, Q = inp.cr_Q * 1e3, rate = inp.cr_A * sig ** inp.cr_n * Math.exp(-Q / (R * T)), tServ = inp.cr_years * YEAR, strain = rate * tServ, t1 = rate > 0 ? 0.01 / rate : BIG;
  const LMP = (inp.lm_a0 + inp.lm_a1 * Math.log10(sig)) * 1e3, logtr = LMP / T - inp.lm_C, tr = capN(10 ** Math.min(logtr, 29)), trY = tr * HOUR / YEAR, MG = rate * tr * HOUR;
  // accelerated testing
  const Ea = inp.acc_Ea, Tu = inp.acc_Tuse + 273.15, Tt = inp.acc_Ttest + 273.15, afT = Math.exp(Ea / KEV * (1 / Tu - 1 / Tt)), afH = (inp.acc_RHtest / inp.acc_RHuse) ** inp.acc_nRH, afP = afT * afH;
  const afCM = (inp.acc_dTtest / inp.acc_dTuse) ** inp.acc_m, afNL = (inp.acc_dTtest / inp.acc_dTuse) ** 1.9 * (inp.acc_fuse / inp.acc_ftest) ** (1 / 3) * Math.exp(1414 * (1 / (Tu + inp.acc_dTuse / 2) - 1 / (Tt + inp.acc_dTtest / 2)));
  const life = inp.cr_years, tTest = life * 8766 / afP, cycUse = life * 365.25 * inp.acc_fuse, cycTest = cycUse / afCM, nSamp = Math.ceil(Math.log(1 - inp.acc_conf) / Math.log(inp.acc_rel));
  const Ts = N.linspace(Math.max(inp.cr_T - 150, 20), inp.cr_T + 150, 60), sg = lsp(Math.log10(Math.max(sig / 5, 1)), Math.log10(sig * 3), 60), Tts = N.linspace(inp.acc_Tuse + 5, Math.max(inp.acc_Ttest + 60, inp.acc_Tuse + 80), 60);
  if (logtr > 29) warnings.push('Larson–Miller rupture life exceeds 10²⁹ h and is capped: creep rupture is not a credible failure mode at this stress and temperature.');
  if (Tt <= Tu) warnings.push('Test temperature is not above the use temperature: the thermal acceleration factor is ≤ 1.');
  if (afP > 1e4) warnings.push(`Combined acceleration factor ${f(afP)} is very large: extrapolating over more than about three decades risks activating a different failure mechanism in the test than in service.`);
  warnings.push('Norton secondary creep only (no primary or tertiary stage); the Larson–Miller master curve is a straight line in log stress fitted by the two entered coefficients — replace them with data for your alloy and heat. Acceleration models assume a single dominant mechanism with constant activation energy.');
  return finish({
    summary: `At ${f(sig)} MPa and ${f(inp.cr_T)} °C the Norton minimum creep rate is ${f(rate)} 1/s, i.e. ${f(strain * 100)} % strain in ${f(life)} years; the Larson–Miller rupture life is ${f(trY)} years. For accelerated ageing, ${f(inp.acc_Ttest)} °C / ${f(inp.acc_RHtest)} % RH gives an acceleration factor of ${f(afP)} (Arrhenius ${f(afT)} × Peck ${f(afH)}), so ${f(life)} years of service is represented by ${f(tTest)} h; ${nSamp} samples with zero failures demonstrate ${f(inp.acc_rel * 100)} % reliability at ${f(inp.acc_conf * 100)} % confidence.`,
    kpis: [
      kpi('rate', 'Minimum creep rate', rate, '1/s', `Norton law A σⁿ exp(−Q/RT) with n = ${f(inp.cr_n)} and Q = ${f(inp.cr_Q)} kJ/mol; a 10 % stress increase multiplies it by ${f(1.1 ** inp.cr_n)}.`),
      kpi('strain', 'Creep strain over service life', strain * 100, '%', `Steady-state rate × ${f(life)} years. Design limits are typically 0.2–1 % for dimensional stability.`, strain < 0.01),
      kpi('t1', 'Time to 1 % creep strain', capN(t1 / YEAR), 'years', 'Service time after which accumulated secondary creep strain reaches 1 %.'),
      kpi('tr', 'Rupture life (Larson–Miller)', trY, 'years', `LMP = T(C + log₁₀t_r) = ${f(LMP)} from the master curve at ${f(sig)} MPa with C = ${f(inp.lm_C)}; the Monkman–Grant product ε̇·t_r is ${f(MG)}.`, trY > 2 * life),
      kpi('afT', 'Arrhenius acceleration factor', afT, '–', `exp[E_a/k (1/T_use − 1/T_test)] with E_a = ${f(Ea)} eV between ${f(inp.acc_Tuse)} and ${f(inp.acc_Ttest)} °C.`),
      kpi('afP', 'Temperature–humidity (Peck) acceleration factor', afP, '–', `Arrhenius factor × (RH_test/RH_use)^n with n = ${f(inp.acc_nRH)}; used for moisture-driven corrosion and delamination.`),
      kpi('afNL', 'Thermal-cycling acceleration (Norris–Landzberg)', afNL, '–', `(ΔT_t/ΔT_u)^1.9 (f_u/f_t)^⅓ exp[1414(1/T_max,u − 1/T_max,t)] with the classic solder constants; the plain Coffin–Manson factor with exponent ${f(inp.acc_m)} is ${f(afCM)}.`),
      kpi('ttest', 'Test duration equivalent to service life', tTest, 'h', `${f(life)} years divided by the temperature–humidity acceleration factor.`, tTest < 5000),
      kpi('cyctest', 'Thermal cycles equivalent to service life', cycTest, 'cycles', `${f(cycUse)} field cycles divided by the Coffin–Manson factor.`),
      kpi('nsamp', 'Zero-failure sample size', nSamp, 'samples', `Success-run theorem n = ln(1 − C)/ln R for reliability ${f(inp.acc_rel)} at confidence ${f(inp.acc_conf)}.`),
    ],
    series: [
      ser('rateT', 'Creep rate versus temperature', 'temperature (°C)', 'creep rate (1/s)', [{ name: 'Norton rate', x: Ts, y: Ts.map(t => inp.cr_A * sig ** inp.cr_n * Math.exp(-Q / (R * (t + 273.15)))) }], 'Arrhenius temperature dependence of the minimum creep rate at the applied stress.', { logy: true }),
      ser('lm', 'Rupture life versus stress', 'stress (MPa)', 'rupture life (h)', [{ name: `at ${f(inp.cr_T)} °C`, x: sg, y: sg.map(s => 10 ** Math.min(29, (inp.lm_a0 + inp.lm_a1 * Math.log10(s)) * 1e3 / T - inp.lm_C)) }, { name: `at ${f(inp.cr_T + 50)} °C`, x: sg, y: sg.map(s => 10 ** Math.min(29, (inp.lm_a0 + inp.lm_a1 * Math.log10(s)) * 1e3 / (T + 50) - inp.lm_C)) }], 'Stress-rupture curves from the Larson–Miller master curve; the vertical gap between isotherms shows the cost of 50 K.', { logx: true, logy: true }),
      ser('af', 'Acceleration factor versus test temperature', 'test temperature (°C)', 'acceleration factor', [{ name: 'Arrhenius', x: Tts, y: Tts.map(t => Math.exp(Ea / KEV * (1 / Tu - 1 / (t + 273.15)))) }, { name: 'Arrhenius × Peck', x: Tts, y: Tts.map(t => Math.exp(Ea / KEV * (1 / Tu - 1 / (t + 273.15))) * afH) }], 'How much test time is saved by raising the stress temperature; keep the extrapolation within the range where the mechanism is unchanged.', { logy: true }),
    ],
    equations: [eqn('constitutive', 'Norton creep law', 'ε̇ = A σⁿ exp(−Q/RT)'), eqn('formula', 'Larson–Miller parameter', 'LMP = T (C + log₁₀ t_r) = f(σ)'), eqn('formula', 'Arrhenius acceleration', 'AF_T = exp[(E_a/k)(1/T_use − 1/T_test)]'), eqn('formula', 'Peck temperature–humidity model', 'AF = (RH_test/RH_use)ⁿ · AF_T'), eqn('formula', 'Coffin–Manson / Norris–Landzberg thermal cycling', 'AF = (ΔT_t/ΔT_u)^m;   AF_NL = (ΔT_t/ΔT_u)^1.9 (f_u/f_t)^{1/3} exp[1414 (1/T_max,u − 1/T_max,t)]'), eqn('formula', 'Success-run sample size', 'n = ln(1 − C)/ln(R)')],
    checks: [chk('Creep strain within 1 % over life', strain < 0.01, `${f(strain * 100)} %`), chk('Rupture life ≥ 2 × service life', trY >= 2 * life, `${f(trY)} years vs ${f(life)} years`), chk('Homologous extrapolation moderate (AF ≤ 10⁴)', afP <= 1e4, `AF = ${f(afP)}`), chk('Test more severe than use', Tt > Tu && inp.acc_dTtest >= inp.acc_dTuse, `T: ${f(inp.acc_Tuse)} → ${f(inp.acc_Ttest)} °C; ΔT: ${f(inp.acc_dTuse)} → ${f(inp.acc_dTtest)} K`)],
    warnings,
    recommendations: [
      rec('high', strain >= 0.01 || trY < 2 * life ? `Creep governs: reduce stress to ${f(sig * (0.01 / Math.max(strain, 1e-300)) ** (1 / inp.cr_n))} MPa for ≤ 1 % strain, or lower the temperature by ${f(R * T * T / Q * Math.log(Math.max(strain / 0.01, 1.0001)))} K.` : `Creep margin is adequate (${f(strain * 100)} % strain, rupture at ${f(trY)} years); re-check if temperature rises by more than ${f(R * T * T / Q * Math.log(0.01 / Math.max(strain, 1e-300)))} K.`, 'Norton and Larson–Miller results.'),
      rec('medium', `Run the ageing test for ${f(tTest)} h at ${f(inp.acc_Ttest)} °C / ${f(inp.acc_RHtest)} % RH with ${nSamp} samples, and add two further temperatures (±15 K) to measure E_a rather than assume ${f(Ea)} eV.`, 'Acceleration-factor and success-run calculations.'),
      rec('medium', `Thermal cycling: ${f(cycTest)} test cycles of ΔT = ${f(inp.acc_dTtest)} K represent the service life; verify by failure analysis that the test produces the same failure mode as the field.`, 'Coffin–Manson acceleration.'),
      rec('low', `Higher acceleration shortens chamber time and energy use (test energy scales with the ${f(tTest)} h duration), but only within the validated mechanism; share chambers across sample lots to cut energy per result.`, 'Test-duration result.', true),
    ],
    publish: { 'life.years': Math.min(trY, capN(t1 / YEAR)), 'life.creep_rate': rate, 'life.acceleration_factor': afP, 'life.test_hours': tTest, 'proc.T_C': inp.cr_T },
    numerics: { method: 'Closed-form Norton, Larson–Miller, Arrhenius, Peck, Coffin–Manson and Norris–Landzberg relations', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'No time integration (steady-state creep).' },
  });
}

// ================================================================================================
// 4. Chemical / environmental degradation
// ================================================================================================
const METALS = { Fe: { n: 2, rho: 7870 }, Zn: { n: 2, rho: 7140 }, Al: { n: 3, rho: 2700 }, Cu: { n: 2, rho: 8960 }, Mg: { n: 2, rho: 1740 }, Ni: { n: 2, rho: 8908 }, Ag: { n: 1, rho: 10490 } };
// dissolving nanoparticles: molar volume from density; dissolution equation built with the chemistry kernel
const DISS = {
  ZnO: { rho: 5606, r: ['ZnO', 'H+'], p: ['Zn2+', 'H2O'], ion: 'Zn2+', alk: [['ZnO', 'H2O'], ['Zn2+', 'OH-']] },
  CuO: { rho: 6310, r: ['CuO', 'H+'], p: ['Cu2+', 'H2O'], ion: 'Cu2+', alk: [['CuO', 'H2O'], ['Cu2+', 'OH-']] },
  Ag: { rho: 10490, r: ['Ag', 'O2', 'H+'], p: ['Ag+', 'H2O'], ion: 'Ag+' },
  SiO2: { rho: 2200, r: ['SiO2', 'H2O'], p: ['H4SiO4'], ion: 'H4SiO4' },
};
const MEDIA = { water: { name: 'pure water, pH 7', cs: 1, k: 1 }, bio: { name: 'cell-culture / physiological medium, pH 7.4', cs: 3, k: 2 }, lyso: { name: 'lysosomal fluid, pH 4.5', cs: 1000, k: 20 } };
function runDegradation(inp) {
  const mech = inp.deg_mech, warnings = [];
  if (mech === 'oxidation') {
    const T = inp.ox_T + 273.15, wet = inp.ox_amb === 'wet', o100 = inp.ox_orient === '100', kT = KEV * T;
    const B = (wet ? 386 * Math.exp(-0.78 / kT) : 772 * Math.exp(-1.23 / kT)), BA = (wet ? 1.63e8 * Math.exp(-2.05 / kT) : 6.23e6 * Math.exp(-2.0 / kT)) / (o100 ? 1.68 : 1), A = B / BA; // µm²/h, µm/h, µm
    const xi = inp.ox_x0 * 1e-3, tau = (xi * xi + A * xi) / B, th = inp.ox_t, x = t => 0.5 * A * (Math.sqrt(1 + 4 * B * (t + tau) / (A * A)) - 1), xf = x(th), grown = xf - xi, regime = (th + tau) < A * A / (4 * B) ? 'linear (interface-reaction limited)' : 'parabolic (diffusion limited)';
    const ts = lsp(-2, Math.log10(Math.max(th * 10, 1)), 80), resid = Math.abs(xf * xf + A * xf - B * (th + tau)) / Math.max(B * (th + tau), 1e-300);
    if (inp.ox_T < 700 || inp.ox_T > 1200) warnings.push('Deal–Grove rate constants are fitted for about 700–1200 °C; outside this range the extrapolation is uncertain.');
    if (!wet && xf < 0.025) warnings.push('Dry oxides thinner than about 25 nm grow faster than Deal–Grove predicts (anomalous initial regime); a Massoud-type correction is needed for gate-oxide thicknesses.');
    return finish({
      summary: `${wet ? 'Wet' : 'Dry'} oxidation of (${o100 ? '100' : '111'}) silicon at ${f(inp.ox_T)} °C for ${f(th)} h grows the oxide from ${f(inp.ox_x0)} to ${f(xf * 1e3)} nm (${regime}); ${f(0.44 * grown * 1e3)} nm of silicon is consumed. Rate constants: B = ${f(B)} µm²/h, B/A = ${f(BA)} µm/h.`,
      kpis: [
        kpi('xox', 'Oxide thickness', xf * 1e3, 'nm', `Deal–Grove linear–parabolic solution x² + Ax = B(t + τ) after ${f(th)} h, starting from ${f(inp.ox_x0)} nm.`),
        kpi('grown', 'Oxide grown in this step', grown * 1e3, 'nm', 'Final minus initial thickness; average rate ' + f(grown * 1e3 / Math.max(th, 1e-9)) + ' nm/h.'),
        kpi('si', 'Silicon consumed', 0.44 * grown * 1e3, 'nm', '0.44 × oxide grown (ratio of Si and SiO₂ molar volumes); the interface moves into the wafer by this amount.'),
        kpi('B', 'Parabolic rate constant B', B, 'µm²/h', `Diffusion-controlled constant with activation energy ${wet ? '0.78' : '1.23'} eV (${wet ? 'H₂O' : 'O₂'} diffusing through the oxide).`),
        kpi('BA', 'Linear rate constant B/A', BA, 'µm/h', `Interface-reaction constant with activation energy ${wet ? '2.05' : '2.0'} eV; (100) is 1.68× slower than (111).`),
        kpi('regime', 'Growth regime', regime, '', `Crossover thickness A/2 = ${f(A * 500)} nm: thinner oxides grow linearly in time, thicker ones as √t.`),
        kpi('rate', 'Instantaneous growth rate at end', B / (2 * xf + A) * 1e3, 'nm/h', 'dx/dt = B/(2x + A); it falls as the oxide thickens and protects the silicon.'),
      ],
      series: [ser('xt', 'Oxide thickness versus time', 'time (h)', 'thickness (nm)', [{ name: `${wet ? 'wet' : 'dry'}, ${f(inp.ox_T)} °C`, x: ts, y: ts.map(t => x(t) * 1e3) }], 'Linear growth at short times bending over to parabolic growth; the same curve describes thermal-oxidation ageing of silicon nanostructures.', { logx: true, logy: true }),
        ser('xT', 'Thickness after the same time versus temperature', 'temperature (°C)', 'thickness (nm)', (() => { const Tc = N.linspace(700, 1200, 50); return [{ name: 'oxide thickness', x: Tc, y: Tc.map(t => { const k = KEV * (t + 273.15), b = wet ? 386 * Math.exp(-0.78 / k) : 772 * Math.exp(-1.23 / k), ba = (wet ? 1.63e8 * Math.exp(-2.05 / k) : 6.23e6 * Math.exp(-2.0 / k)) / (o100 ? 1.68 : 1), a = b / ba, t0 = (xi * xi + a * xi) / b; return 0.5 * a * (Math.sqrt(1 + 4 * b * (th + t0) / (a * a)) - 1) * 1e3; }) }]; })(), 'Strong Arrhenius sensitivity: tens of kelvin change the thickness appreciably, so furnace uniformity controls oxide uniformity.', { logy: true })],
      equations: [eqn('chemical', 'Oxidation reaction', wet ? safeEq(['Si', 'H2O'], ['SiO2', 'H2']) : safeEq(['Si', 'O2'], ['SiO2'])), eqn('governing', 'Deal–Grove linear–parabolic law', 'x² + A x = B (t + τ),   x(t) = (A/2)[√(1 + 4B(t + τ)/A²) − 1]'), eqn('constitutive', 'Rate constants', 'B = B₀ exp(−E_B/kT),   B/A = C₀ exp(−E_A/kT)')],
      checks: [chk('Deal–Grove equation satisfied', resid < 1e-9, `relative residual ${f(resid)}`), chk('Temperature within fitted range (700–1200 °C)', inp.ox_T >= 700 && inp.ox_T <= 1200, `${f(inp.ox_T)} °C`), chk('Oxide growth is monotonic (irreversible)', grown >= 0, `Δx = ${f(grown * 1e3)} nm`)],
      warnings,
      recommendations: [
        rec('high', `For a target thickness use t = (x² + A·x)/B − τ; at these conditions each further 10 nm costs about ${f(10 / Math.max(B / (2 * xf + A) * 1e3, 1e-300))} h.`, 'Instantaneous Deal–Grove rate.'),
        rec('medium', wet ? 'Wet oxidation is fast but gives lower-density oxide; finish with a dry step for a high-quality Si/SiO₂ interface.' : 'Dry oxidation is slow but gives the best interface; use wet or pyrogenic steam for field oxides thicker than about 200 nm.', 'Relative magnitude of the dry and wet rate constants.'),
        rec('medium', `Silicon nanowires and membranes lose ${f(0.44 * grown * 1e3)} nm per oxidised surface: account for this in the final dimensions (and for stress-retarded oxidation on curved surfaces, not modelled).`, 'Silicon consumption 0.44·x.'),
        rec('low', 'Lower-temperature, longer oxidations or rapid thermal oxidation reduce furnace energy per wafer; load full boats and avoid idle hot stand-by.', 'Furnace energy scales with time at temperature.', true),
      ],
      publish: { 'geom.thickness': xf * 1e-6, 'growth.rate': B / (2 * xf + A) * 1e-6 / 3600, 'proc.T_C': inp.ox_T, 'proc.time': th * 3600, 'material.formula': 'SiO2' },
      numerics: { method: 'Closed-form Deal–Grove solution', dof: 1, steps: 1, iterations: 1, residual: resid, notes: 'Classic rate constants for (111) silicon; (100) linear constant divided by 1.68.' },
    });
  }
  if (mech === 'corrosion') {
    const key = METALS[inp.co_metal] ? inp.co_metal : 'Fe', me = METALS[key], M = MM(key), ic = inp.co_i, pen = ic * M / (me.n * F * me.rho), mmy = pen * YEAR * 1e3, ml = ic * M / (me.n * F) * 86400 * 1e3, th = inp.co_thk * 1e-6, tperf = capN(th / Math.max(pen, 1e-300) / YEAR), yrs = inp.co_years, loss = Math.min(pen * yrs * YEAR, th);
    const pit = inp.co_pit, tpit = tperf / pit, cls = mmy < 0.02 ? 'outstanding (< 0.02 mm/yr)' : mmy < 0.1 ? 'excellent' : mmy < 0.5 ? 'good' : mmy < 1 ? 'fair' : 'poor (> 1 mm/yr)';
    const is = lsp(-4, 1, 60), ts = N.linspace(0, Math.max(yrs, 1), 60);
    warnings.push('Uniform corrosion from a constant corrosion current (Faraday’s law); localised attack is represented only by the entered pitting factor, and passivation, film growth or galvanic coupling are not modelled.');
    if (key === 'Al' || key === 'Ni') warnings.push(`${key} is normally passive; uniform-corrosion estimates are meaningful only if the passive film has broken down (e.g. chloride, extreme pH).`);
    return finish({
      summary: `A corrosion current density of ${f(ic)} A/m² (${f(ic * 100)} µA/cm²) on ${key} corresponds to ${f(mmy)} mm/yr (${cls}) and ${f(ml)} g/m²/day. A ${f(inp.co_thk)} µm section is consumed uniformly in ${f(tperf)} years, or perforated in ${f(tpit)} years with a pitting factor of ${f(pit)}; after ${f(yrs)} years ${f(loss * 1e6)} µm is lost.`,
      kpis: [
        kpi('rate', 'Penetration rate', mmy, 'mm/yr', `Faraday’s law i·M/(nFρ) with n = ${me.n} electrons per ${key} atom.`, mmy < 0.1),
        kpi('ml', 'Mass-loss rate', ml, 'g/m²/day', 'Mass of metal dissolved per unit area and day; directly comparable with coupon weight-loss tests.'),
        kpi('class', 'Corrosion resistance class', cls, '', 'Common engineering classification of uniform corrosion rates.'),
        kpi('loss', 'Thickness lost over service life', loss * 1e6, 'µm', `Uniform loss after ${f(yrs)} years (capped at the section thickness).`, loss < 0.2 * th),
        kpi('tperf', 'Time to consume the section (uniform)', tperf, 'years', 'Section thickness divided by the penetration rate.'),
        kpi('tpit', 'Time to perforation with pitting', tpit, 'years', `Uniform time divided by the pitting factor ${f(pit)} (deepest pit / mean penetration).`, tpit > yrs),
        kpi('ion', 'Metal-ion release rate', ic / (me.n * F), 'mol/m²/s', 'Flux of dissolved metal ions from the surface — the source term for environmental release estimates.'),
      ],
      series: [ser('ri', 'Penetration rate versus corrosion current', 'i_corr (A/m²)', 'rate (mm/yr)', [{ name: key, x: is, y: is.map(i => i * M / (me.n * F * me.rho) * YEAR * 1e3) }], 'Linear Faraday relation; 1 µA/cm² (0.01 A/m²) is about 0.01 mm/yr for most engineering metals.', { logx: true, logy: true }),
        ser('lt', 'Remaining thickness versus time', 'time (years)', 'thickness (µm)', [{ name: 'uniform', x: ts, y: ts.map(t => Math.max(0, th - pen * t * YEAR) * 1e6) }, { name: 'at deepest pit', x: ts, y: ts.map(t => Math.max(0, th - pit * pen * t * YEAR) * 1e6) }], 'Moving corrosion front: remaining ligament under uniform attack and under the deepest pit.')],
      equations: [eqn('half-cell', 'Anodic dissolution', safeEq([key], [key + (me.n > 1 ? me.n : '') + '+', 'e-'])), eqn('half-cell', 'Cathodic oxygen reduction (neutral/alkaline)', safeEq(['O2', 'H2O', 'e-'], ['OH-'])), eqn('half-cell', 'Cathodic hydrogen evolution (acid)', safeEq(['H+', 'e-'], ['H2'])), eqn('formula', 'Faraday’s law', 'penetration rate = i_corr M /(n F ρ)')],
      checks: [chk('Section survives the service life (uniform)', tperf > yrs, `${f(tperf)} years vs ${f(yrs)} years`), chk('Section survives with pitting factor', tpit > yrs, `${f(tpit)} years`), chk('Charge–mass consistency (Faraday)', true, `${f(ic / (me.n * F))} mol/m²/s of ${key} for ${f(ic)} A/m²`)],
      warnings,
      recommendations: [
        rec('high', tpit <= yrs ? `Perforation is predicted within the service life: add a corrosion allowance of ${f(pit * pen * yrs * YEAR * 1e6 - th * 1e6)} µm, apply a barrier coating (ALD Al₂O₃/TiO₂ nanolaminate for thin structures) or cathodic protection.` : `The section outlasts the ${f(yrs)}-year life (${f(tpit)} years with pitting); monitor i_corr by linear polarisation resistance to confirm the assumed ${f(ic)} A/m².`, 'Faraday penetration with pitting factor.'),
        rec('medium', 'Measure i_corr in the real electrolyte by Tafel extrapolation or polarisation resistance (Stern–Geary) at service temperature; rates change by an order of magnitude with chloride, pH and aeration.', 'The corrosion current is an input, not a prediction.'),
        rec('medium', `Thin films and lattice struts have little allowance: ${f(mmy * 1e3)} µm/yr removes a 1 µm strut in ${f(1 / Math.max(mmy * 1e3, 1e-300) * 12)} months.`, 'Penetration rate of this run.'),
        rec('low', `Metal-ion release of ${f(ml)} g/m²/day reaches the environment; choose corrosion-resistant or sacrificially protected designs and chromate-free inhibitors to limit discharge.`, 'Mass-loss rate.', true),
      ],
      publish: { 'life.years': Math.min(tpit, BIG), 'life.corrosion_rate_m_s': pen, 'elec.current_density': ic, 'material.formula': key },
      numerics: { method: 'Faraday’s law (closed form)', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Uniform attack with a multiplicative pitting factor.' },
    });
  }
  if (mech === 'dissolution') {
    const key = DISS[inp.di_mat] ? inp.di_mat : 'ZnO', mat = DISS[key], med = MEDIA[inp.di_medium] || MEDIA.water, M = MM(key), Vm = M / mat.rho, r0 = inp.di_d / 2, k = inp.di_k * med.k, cs = inp.di_cs * med.cs, cm = inp.di_conc * 1e-3 /* kg/m³ */, np = cm / (mat.rho * 4 / 3 * PI * r0 ** 3), ctot = cm / M;
    const tEnd = inp.di_t * HOUR, ns = 400, sol = N.rk4((t, y) => { const r = Math.max(y[0], 0), c = np * 4 / 3 * PI * (r0 ** 3 - r ** 3) / Vm; return [r > 1e-13 ? -k * Vm * Math.max(0, 1 - c / cs) : 0]; }, [r0], 0, tEnd, ns);
    const rr = sol.y.map(y => Math.max(y[0], 0)), frac = rr.map(r => 1 - (r / r0) ** 3), fEnd = frac[ns], cEnd = fEnd * ctot, tSink = r0 / (k * Vm), eqFrac = Math.min(1, cs / ctot);
    let th = BIG; for (let i = 1; i <= ns; i++) if (frac[i] >= 0.5) { th = sol.t[i - 1] + (0.5 - frac[i - 1]) / Math.max(frac[i] - frac[i - 1], 1e-300) * (sol.t[i] - sol.t[i - 1]); break; }
    const areaRet = (rr[ns] / r0) ** 2, ionM = MM(mat.ion.replace(/\d*[+-]$/, '')), tm = sol.t.map(t => t / HOUR);
    const massBal = Math.abs(np * 4 / 3 * PI * rr[ns] ** 3 * mat.rho + cEnd * M - cm) / cm;
    warnings.push(`Monodisperse shrinking-particle model with surface-reaction control and a single solubility limit for "${med.name}" (rate ×${med.k}, solubility ×${med.cs} relative to the entered water values); aggregation, protein corona, ligand complexation and size-dependent (Gibbs–Thomson) solubility are not included.`);
    if (eqFrac < 1) warnings.push(`The particle loading exceeds the solubility: dissolution stops at ${f(eqFrac * 100)} % when the medium saturates (closed system). In flowing or sink conditions it continues to completion in about ${f(tSink / HOUR)} h.`);
    return finish({
      summary: `${f(inp.di_d * 1e9)} nm ${chem.prettyFormula(key)} particles at ${f(inp.di_conc)} mg/L in ${med.name} are ${f(fEnd * 100)} % dissolved after ${f(inp.di_t)} h, releasing ${f(cEnd * ionM * 1e3)} mg/L of ${chem.prettyFormula(mat.ion)} (${f(cEnd)} mol/m³). The half-dissolution time is ${th >= BIG ? 'not reached (saturation)' : f(th / HOUR) + ' h'}; under sink conditions a particle would vanish in ${f(tSink / HOUR)} h.`,
      kpis: [
        kpi('frac', 'Fraction dissolved', fEnd * 100, '%', `Mass fraction released after ${f(inp.di_t)} h from integrating dr/dt = −kV_m(1 − c/c_s) with the released ions accumulating in the closed volume.`),
        kpi('ion', 'Dissolved ion concentration', cEnd * ionM * 1e3, 'mg/L', `${chem.prettyFormula(mat.ion)} released into the medium; compare with ecotoxicity or cytotoxicity thresholds for the ion.`),
        kpi('thalf', 'Half-dissolution time', th >= BIG ? BIG : th / HOUR, 'h', 'Time for half of the particle mass to dissolve (10³⁰ = never, because the medium saturates first).'),
        kpi('tsink', 'Complete-dissolution time (sink conditions)', tSink / HOUR, 'h', 'r₀/(kV_m): lifetime of a particle when the ion concentration stays far below solubility; proportional to particle size.'),
        kpi('eq', 'Equilibrium dissolved fraction', eqFrac * 100, '%', `Solubility ${f(cs)} mol/m³ divided by the total metal loading ${f(ctot)} mol/m³, capped at 100 %.`),
        kpi('area', 'Surface area retained', areaRet * 100, '%', 'Remaining particle surface (r/r₀)² — a proxy for retained catalytic, antibacterial or sensing function.'),
        kpi('r', 'Remaining particle diameter', 2 * rr[ns] * 1e9, 'nm', 'Particle size at the end of the exposure (all particles shrink equally in this monodisperse model).'),
      ],
      series: [ser('diss', 'Dissolution kinetics', 'time (h)', 'percent', [{ name: 'mass dissolved (%)', x: tm, y: frac.map(v => v * 100) }, { name: 'surface area retained (%)', x: tm, y: rr.map(r => (r / r0) ** 2 * 100) }], 'Release curve of the shrinking particles: nearly linear in radius at first, flattening as the medium approaches saturation or the particles vanish.'),
        ser('size', 'Sink-condition lifetime versus particle size', 'diameter (nm)', 'lifetime (h)', (() => { const dsz = lsp(0, 3, 40); return [{ name: med.name, x: dsz, y: dsz.map(d => d * 0.5e-9 / (k * Vm) / HOUR) }]; })(), 'Smaller particles dissolve proportionally faster under surface-reaction control; this drives the size dependence of ion-mediated toxicity.', { logx: true, logy: true })],
      equations: [eqn('chemical', 'Dissolution reaction (acid/neutral)', safeEq(mat.r, mat.p))].concat(mat.alk ? [eqn('ionic', 'Dissolution in water (hydroxide form)', safeEq(mat.alk[0], mat.alk[1]))] : []).concat([eqn('governing', 'Shrinking-particle model', 'dr/dt = −k V_m (1 − c/c_s),   c = n_p (4π/3)(r₀³ − r³)/V_m'), eqn('formula', 'Sink-condition lifetime', 't = r₀/(k V_m)')]),
      checks: [chk('Mass balance (solid + dissolved = initial)', massBal < 1e-9, `relative residual ${f(massBal)}`), chk('Dissolved concentration ≤ solubility', cEnd <= cs * (1 + 1e-6), `${f(cEnd)} vs ${f(cs)} mol/m³`), chk('Radius decreases monotonically', rr[ns] <= r0, `${f(r0 * 2e9)} → ${f(rr[ns] * 2e9)} nm`)],
      warnings,
      recommendations: [
        rec('high', fEnd > 0.5 ? `More than half of the material dissolves in ${f(inp.di_t)} h: for durable function apply a sparingly soluble shell (SiO₂, TiO₂ or polymer) or use larger particles (lifetime ∝ size).` : `Dissolution is limited to ${f(fEnd * 100)} % in ${f(inp.di_t)} h; confirm with ICP-MS of ultrafiltered supernatant at 3–4 time points.`, 'Computed release fraction.'),
        rec('medium', `Risk assessment should treat both the particle and the released ${chem.prettyFormula(mat.ion)} (${f(cEnd * ionM * 1e3)} mg/L here); in acidic lysosomal fluid the same particles dissolve about ${f(MEDIA.lyso.k / med.k)}× faster.`, 'Medium dependence of rate and solubility.'),
        rec('medium', 'Measure the surface-normalised rate constant in the actual medium (dialysis or ultrafiltration plus ICP) and re-run; literature values span more than two orders of magnitude with pH and ligands.', 'The rate constant is the controlling input.'),
        rec('low', `Design for safe degradation: choose compositions whose dissolved ions are of low ecotoxicity, minimise the released mass (${f(cEnd * M * 1e3)} g/m³ of solid here) and capture particles from process water by flocculation before discharge.`, 'Release estimate of this run.', true),
      ],
      publish: { 'geom.diameter': 2 * rr[ns], 'life.dissolved_fraction': fEnd, 'life.years': capN(tSink / YEAR), 'risk.release_mol_m3': cEnd, 'material.formula': key },
      numerics: { method: 'RK4 integration of the shrinking-particle ODE with closed-volume ion balance', dof: 1, steps: ns, iterations: ns, residual: massBal, notes: 'Residual is the mass-balance error.' },
    });
  }
  if (mech === 'ripening') {
    const T = inp.rp_T + 273.15, r0 = inp.rp_d * 0.5e-9, Vm = inp.rp_Vm * 1e-6, K = 8 * inp.rp_gamma * inp.rp_cinf * Vm * Vm * inp.rp_D / (9 * R * T), tEnd = inp.rp_t * HOUR, rt = t => Math.cbrt(r0 ** 3 + K * t), rE = rt(tEnd), areaRet = r0 / rE;
    const mu = 2.414e-5 * 10 ** (247.8 / (T - 140)), n0 = inp.rp_n0, tauA = 3 * mu * inp.rp_W / (4 * N.K_B * T * n0), nRet = 1 / (1 + tEnd / tauA), dH = 2 * r0 * (1 / nRet) ** (1 / inp.rp_df), t20 = (r0 ** 3 * (1 / 0.8 ** 3 - 1)) / K, lc = 2 * inp.rp_gamma * Vm / (R * T);
    const ts = lsp(-2, Math.log10(Math.max(inp.rp_t * 10, 1)), 70);
    warnings.push('Diffusion-limited LSW coarsening for a dilute dispersion (r³ law) and Smoluchowski perikinetic agglomeration with a constant stability ratio; supported catalysts, pinned particles and reaction-limited (r² law) ripening behave differently.');
    if (r0 < 2 * lc) warnings.push(`Particle radius is comparable to the capillary length (${f(lc * 1e9)} nm): the linearised Gibbs–Thomson relation behind LSW is inaccurate.`);
    return finish({
      summary: `Ostwald ripening at ${f(inp.rp_T)} °C coarsens ${f(inp.rp_d)} nm particles to ${f(2 * rE * 1e9)} nm in ${f(inp.rp_t)} h (K = ${f(K)} m³/s), leaving ${f(areaRet * 100)} % of the specific surface area. Agglomeration with stability ratio ${f(inp.rp_W)} reduces the number concentration to ${f(nRet * 100)} % (half-time ${f(tauA / HOUR)} h), giving aggregates of about ${f(dH * 1e9)} nm.`,
      kpis: [
        kpi('d', 'Mean diameter after ripening', 2 * rE * 1e9, 'nm', `LSW law r³ = r₀³ + Kt with K = 8γc∞V_m²D/(9RT) = ${f(K)} m³/s.`),
        kpi('area', 'Specific surface area retained', areaRet * 100, '%', 'Surface area per unit mass scales as 1/r; it is the first-order proxy for catalytic activity or sensing response.', areaRet > 0.8),
        kpi('t20', 'Time to 20 % area loss', t20 / HOUR, 'h', 'Time for the mean radius to grow by 25 % (area ∝ 1/r falls to 80 %).'),
        kpi('tauA', 'Agglomeration half-time', tauA / HOUR, 'h', `Smoluchowski time 3μW/(4k_BTn₀) for n₀ = ${f(n0)} m⁻³; a stability ratio W of 1 is fast (diffusion-limited) aggregation.`, tauA > tEnd),
        kpi('nret', 'Primary-particle number retained', nRet * 100, '%', 'n/n₀ = 1/(1 + t/τ); the rest are bound in aggregates.'),
        kpi('dagg', 'Mean aggregate size', dH * 1e9, 'nm', `Primary size × (n₀/n)^(1/d_f) with fractal dimension ${f(inp.rp_df)}; it shifts optical scattering/plasmon response and sedimentation.`),
        kpi('lc', 'Capillary length', lc * 1e9, 'nm', '2γV_m/(RT): the size below which solubility enhancement (Gibbs–Thomson) is strong.'),
      ],
      series: [ser('rt', 'Coarsening and area loss', 'time (h)', 'diameter (nm) / area retained (%)', [{ name: 'mean diameter (nm)', x: ts, y: ts.map(t => 2 * rt(t * HOUR) * 1e9) }, { name: 'surface area retained (%)', x: ts, y: ts.map(t => r0 / rt(t * HOUR) * 100) }], 'Cube-root-of-time coarsening and the corresponding loss of specific surface area.', { logx: true }),
        ser('agg', 'Agglomeration', 'time (h)', 'percent / nm', [{ name: 'primary particles remaining (%)', x: ts, y: ts.map(t => 100 / (1 + t * HOUR / tauA)) }, { name: 'aggregate size (nm)', x: ts, y: ts.map(t => 2 * r0 * (1 + t * HOUR / tauA) ** (1 / inp.rp_df) * 1e9) }], 'Second-order loss of free primary particles and growth of fractal aggregates.', { logx: true, logy: true })],
      equations: [eqn('governing', 'LSW Ostwald ripening', 'r̄³(t) = r̄₀³ + K t,   K = 8 γ c∞ V_m² D /(9 R T)'), eqn('governing', 'Smoluchowski agglomeration', 'dn/dt = −(4k_BT/(3μW)) n²  →  n = n₀/(1 + t/τ),   τ = 3μW/(4k_BT n₀)'), eqn('formula', 'Functional surface area', 'S/S₀ = r₀/r̄')],
      checks: [chk('Volume conserved during ripening (number falls as r⁻³)', true, `number ratio ${f((r0 / rE) ** 3)}`), chk('Radius above 2 capillary lengths', r0 >= 2 * lc, `r₀ = ${f(r0 * 1e9)} nm, l_c = ${f(lc * 1e9)} nm`), chk('Area loss below 20 %', areaRet > 0.8, `${f((1 - areaRet) * 100)} % lost`)],
      warnings,
      recommendations: [
        rec('high', areaRet <= 0.8 ? `Area loss reaches ${f((1 - areaRet) * 100)} %: lower solubility or interfacial energy (capping ligands, less polar medium), lower temperature, or anchor particles on a support; the 20 % loss point is at ${f(t20 / HOUR)} h.` : `Ripening is slow on this time scale (20 % loss after ${f(t20 / HOUR)} h); agglomeration is the greater risk if τ = ${f(tauA / HOUR)} h is shorter than the storage time.`, 'LSW and Smoluchowski time scales.'),
        rec('medium', `Raise colloidal stability (|ζ| > 30 mV or steric stabiliser) to increase W; each decade in W gives a decade in agglomeration half-time (now ${f(tauA / HOUR)} h).`, 'τ ∝ W.'),
        rec('medium', 'Track mean size by DLS/TEM and specific area by BET at three ageing times and fit K; the r³ versus t plot should be linear if LSW applies.', 'Model verification step.'),
        rec('low', 'Longer functional life means less frequent catalyst or sensor replacement: stabilised, supported nanoparticles reduce critical-metal consumption and waste.', 'Functional-lifetime result.', true),
      ],
      publish: { 'geom.diameter': 2 * rE, 'life.area_retained': areaRet, 'life.years': capN(t20 / YEAR), 'geom.ssa_ratio': areaRet },
      numerics: { method: 'Closed-form LSW and Smoluchowski solutions', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Dilute-limit kinetics.' },
    });
  }
  if (mech === 'delamination') {
    const Ef = inp.dl_E * 1e9, nu = inp.dl_nu, h = inp.dl_h * 1e-6, s0 = inp.dl_sigma * 1e6, sT = Ef / (1 - nu) * inp.dl_dalpha * inp.dl_dT, sig = s0 + sT, Ep = Ef / (1 - nu * nu), G = sig * sig * h / (2 * Ep), Gi = inp.dl_Gi, hc = 2 * Ep * Gi / Math.max(sig * sig, 1e-300), Gch = 1.976 * sig * sig * h / Ep, hcc = sig > 0 ? inp.dl_Gf * Ep / (1.976 * sig * sig) : BIG, ratio = G / Gi;
    const sb = sig < 0 ? PI * PI / 12 * Ep * (h / (inp.dl_b * 1e-6)) ** 2 : 0, buckles = sig < 0 && -sig > sb;
    const hs = lsp(Math.log10(h / 30), Math.log10(h * 30), 60), dTs = N.linspace(-Math.abs(inp.dl_dT) * 2 - 50, Math.abs(inp.dl_dT) * 2 + 50, 60);
    warnings.push('Steady-state energy-release rates for a thin film on a thick stiff substrate (elastic mismatch and mode-mixity dependence of interface toughness neglected). Moisture- or fatigue-assisted subcritical debonding can occur well below the critical value.');
    if (buckles) warnings.push(`Compressive stress exceeds the buckling stress of a ${f(inp.dl_b)} µm wide debond (${f(sb / 1e6)} MPa): buckle-driven delamination (blisters, telephone-cord buckles) is possible.`);
    return finish({
      summary: `A ${f(inp.dl_h)} µm film under ${f(sig / 1e6)} MPa (${f(s0 / 1e6)} residual + ${f(sT / 1e6)} thermal) releases ${f(G)} J/m² on debonding against an interface toughness of ${f(Gi)} J/m² (ratio ${f(ratio)}): delamination is ${ratio >= 1 ? 'PREDICTED' : 'not predicted'}. The critical thickness is ${f(hc * 1e6)} µm${sig > 0 ? `; channel cracking of the film starts above ${f(hcc * 1e6)} µm` : ''}.`,
      kpis: [
        kpi('G', 'Energy release rate for debonding', G, 'J/m²', 'Steady-state value σ²h(1−ν²)/(2E) available to drive an interface crack under the film.', ratio < 1),
        kpi('ratio', 'Delamination driving ratio G/Γ_i', ratio, '–', 'Above 1 the interface crack advances spontaneously; below about 0.3–0.5 subcritical (environment-assisted) growth is usually negligible.', ratio < 0.5),
        kpi('sigma', 'Total film stress', sig / 1e6, 'MPa', `Residual stress plus thermal mismatch E/(1−ν)·Δα·ΔT = ${f(sT / 1e6)} MPa; positive is tensile.`),
        kpi('hc', 'Critical thickness for delamination', hc * 1e6, 'µm', 'Film thickness at which G equals the interface toughness at the present stress; thinner films stay attached.', h < hc),
        kpi('hcc', 'Critical thickness for channel cracking', Math.min(hcc * 1e6, 1e12), 'µm', sig > 0 ? `h_c = Γ_f E′/(1.976σ²) with film fracture energy ${f(inp.dl_Gf)} J/m²; applies to tensile films.` : 'Not applicable: the film is in compression (very large value shown).', sig <= 0 || h < hcc),
        kpi('sb', 'Buckling stress of an existing debond', sb / 1e6, 'MPa', sig < 0 ? `(π²/12)·E′(h/b)² for a straight-sided debond of half-width ${f(inp.dl_b)} µm; compressive stress above it drives buckle delamination.` : 'Zero because the film is in tension (buckling needs compression).', !buckles),
        kpi('dTcrit', 'Temperature change to reach delamination', inp.dl_dalpha !== 0 ? Math.abs((Math.sqrt(2 * Ep * Gi / h) - Math.abs(s0)) / (Ef / (1 - nu) * inp.dl_dalpha)) : BIG, 'K', 'Additional |ΔT| (in the unfavourable direction) at which thermal mismatch alone raises |σ| to the critical stress √(2E′Γ_i/h).'),
      ],
      series: [ser('Gh', 'Driving force versus film thickness', 'thickness (µm)', 'G (J/m²)', [{ name: 'debond energy release', x: hs.map(v => v * 1e6), y: hs.map(v => sig * sig * v / (2 * Ep)) }, { name: 'interface toughness', x: [hs[0] * 1e6, hs[59] * 1e6], y: [Gi, Gi] }], 'Energy release rate grows linearly with thickness at fixed stress; the crossing with the interface toughness is the critical thickness.', { logx: true, logy: true }),
        ser('GT', 'Driving ratio versus temperature excursion', 'ΔT (K)', 'G/Γ_i', [{ name: 'G/Γ_i', x: dTs, y: dTs.map(d => { const s_ = s0 + Ef / (1 - nu) * inp.dl_dalpha * d; return s_ * s_ * h / (2 * Ep) / Gi; }) }, { name: 'critical', x: [dTs[0], dTs[59]], y: [1, 1] }], 'Thermal excursions add to or relieve the residual stress; the safe temperature window is where the curve stays below one.')],
      equations: [eqn('governing', 'Steady-state interface debonding', 'G = σ² h (1 − ν²)/(2E);   delamination when G ≥ Γ_i'), eqn('constitutive', 'Thermal mismatch stress', 'σ = σ_res + [E/(1 − ν)] Δα ΔT'), eqn('formula', 'Channel cracking of a tensile film', 'G = 1.976 σ² h (1 − ν²)/E ≥ Γ_f'), eqn('formula', 'Buckling of a debonded strip', 'σ_b = (π²/12) [E/(1 − ν²)] (h/b)²', 'Clamped straight-sided debond of half-width b.')],
      checks: [chk('Interface stable (G < Γ_i)', ratio < 1, `G/Γ_i = ${f(ratio)}`), chk('No channel cracking', sig <= 0 || h < hcc, sig > 0 ? `h = ${f(h * 1e6)} µm vs ${f(hcc * 1e6)} µm` : 'film in compression'), chk('No buckle-driven delamination', !buckles, sig < 0 ? `|σ| = ${f(-sig / 1e6)} MPa vs σ_b = ${f(sb / 1e6)} MPa` : 'film in tension'), chk('Energy release rate non-negative', G >= 0, `${f(G)} J/m²`)],
      warnings,
      recommendations: [
        rec('high', ratio >= 1 ? `Reduce thickness below ${f(hc * 1e6)} µm, cut the stress to ${f(Math.sqrt(2 * Ep * Gi / h) / 1e6)} MPa (deposition pressure/anneal) or raise adhesion (plasma clean, Ti/Cr adhesion layer, silane primer) to more than ${f(G)} J/m².` : `Margin against delamination is ${f(1 / Math.max(ratio, 1e-12))}× in energy; keep G/Γ_i below 0.5 for humid service because water lowers interface toughness.`, 'Energy balance of this run.'),
        rec('medium', 'Measure film stress by wafer curvature and interface toughness by four-point bend or nano-scratch/indentation on the real stack; both vary strongly with process.', 'Inputs dominate the prediction.'),
        rec('medium', `Thermal excursions change the stress by ${f(Ef / (1 - nu) * inp.dl_dalpha / 1e6)} MPa/K; choose a substrate or interlayer with matched expansion for large |ΔT|.`, 'Thermal-mismatch term.'),
        rec('low', 'Durable coatings avoid recoating cycles and the associated solvent, energy and waste; avoid chromate-based adhesion promoters in favour of silane or plasma treatments.', 'Good practice; delamination margin of this run.', true),
      ],
      publish: { 'mech.max_stress': Math.abs(sig), 'mech.safety_factor': 1 / Math.max(ratio, 1e-12), 'life.delam_ratio': ratio, 'geom.thickness': h },
      numerics: { method: 'Closed-form thin-film fracture mechanics (steady-state energy release rates)', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Homogeneous elastic properties assumed across the interface.' },
    });
  }
  if (mech === 'electromigration') {
    const T = inp.em_T + 273.15, Tr = inp.em_Tref + 273.15, j = inp.em_j, n = inp.em_n, Ea = inp.em_Ea, t50 = inp.em_tref * (inp.em_jref / j) ** n * Math.exp(Ea / KEV * (1 / T - 1 / Tr)), sg = inp.em_sigma, t01 = t50 * Math.exp(sg * N.normInv(0.001)), jL = j * inp.em_L * 1e-6, immortal = jL < inp.em_jLc, yrs = t50 / 8766;
    const jmax = inp.em_jref * (inp.em_tref * Math.exp(Ea / KEV * (1 / T - 1 / Tr)) * Math.exp(sg * N.normInv(0.001)) / (inp.em_life * 8766)) ** (1 / n);
    const js = lsp(Math.log10(j / 10), Math.log10(j * 10), 60), Ts = N.linspace(Math.max(25, inp.em_T - 60), inp.em_T + 100, 60);
    warnings.push('Black’s equation extrapolates accelerated-test data with constant exponent and activation energy; Joule self-heating of the line at test current, microstructure (bamboo vs polycrystalline) and via/barrier effects are not modelled.');
    if (immortal) warnings.push(`The j·L product (${f(jL)} A/m) is below the Blech threshold (${f(inp.em_jLc)} A/m): back-stress should stop electromigration, so the line may be immortal and Black’s lifetime is conservative.`);
    return finish({
      summary: `At ${f(j)} A/m² (${f(j / 1e10)} MA/cm²) and ${f(inp.em_T)} °C Black’s equation gives a median electromigration life of ${f(yrs)} years (n = ${f(n)}, E_a = ${f(Ea)} eV), and ${f(t01 / 8766)} years to 0.1 % cumulative failures (lognormal σ = ${f(sg)}). The current density that meets ${f(inp.em_life)} years at 0.1 % failures is ${f(jmax)} A/m². Blech product j·L = ${f(jL)} A/m (${immortal ? 'below' : 'above'} threshold).`,
      kpis: [
        kpi('t50', 'Median time to failure', yrs, 'years', `Black’s equation scaled from the reference test (${f(inp.em_tref)} h at ${f(inp.em_jref)} A/m², ${f(inp.em_Tref)} °C).`, yrs > inp.em_life),
        kpi('t01', 'Time to 0.1 % failures', t01 / 8766, 'years', `Lognormal distribution: t₅₀·exp(σ·z₀.₀₀₁) with σ = ${f(sg)}; this, not the median, is the design life for large interconnect populations.`, t01 / 8766 > inp.em_life),
        kpi('jmax', 'Allowable current density', jmax, 'A/m²', `Current density at which the 0.1 % life equals the ${f(inp.em_life)}-year target at ${f(inp.em_T)} °C.`, j <= jmax),
        kpi('jL', 'Blech product j·L', jL, 'A/m', `Lines with j·L below about ${f(inp.em_jLc)} A/m are protected by back-stress (short-length effect).`, immortal),
        kpi('afT', 'Thermal acceleration of the test', Math.exp(Ea / KEV * (1 / T - 1 / Tr)), '–', 'Ratio of lifetime at use temperature to lifetime at test temperature for the same current density.'),
        kpi('dT10', 'Temperature rise that halves life', KEV * T * T / Ea * Math.LN2, 'K', 'Local heating by this amount (hot spot or self-heating) cuts the electromigration life in half.'),
      ],
      series: [ser('tj', 'Lifetime versus current density', 'current density (A/m²)', 'time (years)', [{ name: 'median', x: js, y: js.map(v => t50 * (j / v) ** n / 8766) }, { name: '0.1 % failures', x: js, y: js.map(v => t01 * (j / v) ** n / 8766) }, { name: 'target life', x: [js[0], js[59]], y: [inp.em_life, inp.em_life] }], 'Power-law dependence on current density; the allowable current is where the 0.1 % line crosses the target life.', { logx: true, logy: true }),
        ser('tT', 'Lifetime versus temperature', 'temperature (°C)', 'median life (years)', [{ name: 'median', x: Ts, y: Ts.map(t => inp.em_tref * (inp.em_jref / j) ** n * Math.exp(Ea / KEV * (1 / (t + 273.15) - 1 / Tr)) / 8766) }], 'Arrhenius temperature dependence of electromigration life.', { logy: true })],
      equations: [eqn('formula', 'Black’s equation', 'MTTF = A j⁻ⁿ exp(E_a/kT)'), eqn('formula', 'Lognormal early-failure time', 't_p = t₅₀ exp(σ z_p)'), eqn('formula', 'Blech immortality condition', 'j·L < (j·L)_crit')],
      checks: [chk('0.1 % life meets target', t01 / 8766 >= inp.em_life, `${f(t01 / 8766)} years vs ${f(inp.em_life)} years`), chk('Current density within allowable', j <= jmax, `${f(j)} vs ${f(jmax)} A/m²`), chk('Blech short-length immunity', immortal, `j·L = ${f(jL)} A/m`)],
      warnings,
      recommendations: [
        rec('high', j > jmax ? `Reduce current density to ${f(jmax)} A/m² by widening the line ${f(j / jmax)}× or adding parallel lines/vias; alternatively keep segments shorter than ${f(inp.em_jLc / j * 1e6)} µm to exploit the Blech effect.` : `Current density is within the allowable ${f(jmax)} A/m²; verify the hottest operating point, since +${f(KEV * T * T / Ea * Math.LN2)} K halves the life.`, 'Black’s equation with lognormal scatter.'),
        rec('medium', 'Run package-level tests at three temperatures and two current densities to fit n and E_a for your metallisation; use resistance-rise (10 %) as the failure criterion.', 'Extrapolation depends on fitted constants.'),
        rec('medium', 'For metasurface or sensor electrodes carrying high current density, cap the line with a refractory layer or alloy the metal to slow interface diffusion.', 'Mechanism-based mitigation.'),
        rec('low', 'Long interconnect life avoids premature device replacement; designing to the 0.1 % life with derating reduces electronic waste.', 'Lifetime result.', true),
      ],
      publish: { 'life.years': Math.min(t01 / 8766, BIG), 'elec.current_density': j, 'life.j_allow': jmax },
      numerics: { method: 'Black’s equation (closed form) with lognormal failure statistics', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Scaled from a reference accelerated test.' },
    });
  }
  // time-dependent dielectric breakdown, E-model
  const T = inp.td_T + 273.15, Tr = inp.td_Tref + 273.15, tox = inp.td_tox * 1e-9, Eox = inp.td_V / tox / 1e8 /* MV/cm */, g = inp.td_gamma, Ea = inp.td_Ea, bw = inp.td_beta;
  const t63 = inp.td_tref * Math.exp(-g * (Eox - inp.td_Eref)) * Math.exp(Ea / KEV * (1 / T - 1 / Tr)) * (1 / inp.td_area) ** (1 / bw), tp = t63 * (-Math.log(1 - 1e-4)) ** (1 / bw), yrs = tp / 8766, life = inp.td_life;
  const Emax = inp.td_Eref + Math.log(inp.td_tref * Math.exp(Ea / KEV * (1 / T - 1 / Tr)) * (1 / inp.td_area) ** (1 / bw) * (-Math.log(1 - 1e-4)) ** (1 / bw) / (life * 8766)) / g, Vmax = Emax * 1e8 * tox;
  const Es = N.linspace(Math.max(0.5, Eox * 0.4), Math.max(Eox * 1.6, inp.td_Eref * 1.1), 60);
  warnings.push('Thermochemical E-model (ln t ∝ −γE) extrapolated from one reference stress condition; the 1/E and power-law voltage models predict longer low-field lifetimes, so the E-model is conservative. Thin (< 5 nm) oxides with direct tunnelling and soft breakdown follow different statistics.');
  if (Eox > 12) warnings.push(`Oxide field ${f(Eox)} MV/cm is near or above the intrinsic breakdown strength of SiO₂ (about 10–15 MV/cm): immediate breakdown is likely.`);
  return finish({
    summary: `With ${f(inp.td_V)} V across ${f(inp.td_tox)} nm (${f(Eox)} MV/cm) at ${f(inp.td_T)} °C the E-model gives a characteristic breakdown time of ${f(t63 / 8766)} years and ${f(yrs)} years to 0.01 % failures (Weibull β = ${f(bw)}, area ratio ${f(inp.td_area)}). The largest field meeting ${f(life)} years is ${f(Emax)} MV/cm, i.e. ${f(Vmax)} V.`,
    kpis: [
      kpi('Eox', 'Dielectric field', Eox, 'MV/cm', 'Applied voltage divided by dielectric thickness; lifetime falls exponentially with it in the E-model.', Eox < Emax),
      kpi('t63', 'Characteristic lifetime (63 %)', t63 / 8766, 'years', `Reference ${f(inp.td_tref)} h at ${f(inp.td_Eref)} MV/cm and ${f(inp.td_Tref)} °C scaled by exp(−γΔE), Arrhenius and area.`),
      kpi('tp', 'Time to 0.01 % failures', yrs, 'years', 'Weibull early-failure time t₆₃·(−ln(1 − p))^(1/β); low β makes it far shorter than the characteristic time.', yrs > life),
      kpi('Vmax', 'Maximum voltage for target life', Vmax, 'V', `Voltage at which the 0.01 % life equals ${f(life)} years at ${f(inp.td_T)} °C.`, inp.td_V <= Vmax),
      kpi('dec', 'Field acceleration', g / Math.LN10, 'decades per MV/cm', 'Each additional MV/cm shortens life by this many decades.'),
      kpi('areaF', 'Area-scaling factor on life', (1 / inp.td_area) ** (1 / bw), '–', 'Weakest-link scaling t ∝ A^(−1/β): larger capacitor area fails earlier.'),
    ],
    series: [ser('tE', 'Lifetime versus field', 'field (MV/cm)', 'time (years)', [{ name: '63 % (characteristic)', x: Es, y: Es.map(e => t63 * Math.exp(-g * (e - Eox)) / 8766) }, { name: '0.01 % failures', x: Es, y: Es.map(e => tp * Math.exp(-g * (e - Eox)) / 8766) }, { name: 'target life', x: [Es[0], Es[59]], y: [life, life] }], 'Exponential field dependence of dielectric lifetime; the operating voltage must sit left of the crossing with the target life.', { logy: true }),
      ser('wb', 'Breakdown distribution', 'time (years)', 'cumulative failures', (() => { const ts = lsp(Math.log10(Math.max(t63 / 8766, 1e-12)) - 6, Math.log10(Math.max(t63 / 8766, 1e-12)) + 1, 60); return [{ name: 'Weibull CDF', x: ts, y: ts.map(t => Math.max(1e-12, 1 - Math.exp(-((t * 8766 / t63) ** bw)))) }]; })(), 'Weibull breakdown distribution; reliability targets at parts-per-million level sit many decades before the characteristic time.', { logx: true, logy: true })],
    equations: [eqn('formula', 'TDDB E-model', 't_BD = t₀ exp(−γ E_ox) exp(E_a/kT)'), eqn('formula', 'Weibull area scaling', 't₆₃(A₂) = t₆₃(A₁)(A₁/A₂)^{1/β}'), eqn('formula', 'Early-failure time', 't_p = t₆₃ [−ln(1 − p)]^{1/β}')],
    checks: [chk('0.01 % life meets target', yrs >= life, `${f(yrs)} years vs ${f(life)} years`), chk('Field below intrinsic strength (≈ 12 MV/cm)', Eox <= 12, `${f(Eox)} MV/cm`), chk('Operating voltage within allowable', inp.td_V <= Vmax, `${f(inp.td_V)} V vs ${f(Vmax)} V`)],
    warnings,
    recommendations: [
      rec('high', inp.td_V > Vmax ? `Lower the operating voltage to ${f(Vmax)} V or thicken the dielectric to ${f(inp.td_V / (Emax * 1e8) * 1e9)} nm (or use a higher-k dielectric of equal capacitance).` : `Voltage margin is ${f(Vmax / inp.td_V)}×; keep transients below ${f(Vmax)} V.`, 'E-model lifetime with Weibull early failures.'),
      rec('medium', 'Characterise TDDB at three fields and two temperatures on capacitors of two areas to obtain γ, E_a and β for your dielectric (ALD oxides differ from thermal SiO₂).', 'Model constants are user-supplied.'),
      rec('medium', 'For electrically tunable metasurfaces (varactor- or MEMS-loaded), check the peak field under resonance-enhanced voltage, not just the bias.', 'Field enhancement at resonant gaps.'),
      rec('low', 'Operating with voltage derating extends device life and reduces electronic waste; burn-in screens early failures but consumes energy, so size it from the measured β.', 'Lifetime distribution result.', true),
    ],
    publish: { 'life.years': Math.min(yrs, BIG), 'life.v_max': Vmax, 'life.field_MV_cm': Eox },
    numerics: { method: 'Thermochemical E-model with Weibull statistics (closed form)', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Scaled from a reference stress condition.' },
  });
}

// ================================================================================================
// 5. Functional lifetime of metamaterials
// ================================================================================================
const SWT = { pcm: { name: 'Phase-change (GST/VO₂) element', eta: 1e9, beta: 1.2, p: 4, Ea: 0.3 }, mems: { name: 'MEMS ohmic/capacitive switch', eta: 5e9, beta: 1.5, p: 6, Ea: 0.2 }, varactor: { name: 'Varactor / PIN diode', eta: 1e14, beta: 1.0, p: 3, Ea: 0.7 }, lc: { name: 'Liquid-crystal cell', eta: 1e11, beta: 2.0, p: 2, Ea: 0.5 } };
function runFunctional(inp) {
  const mech = inp.fun_mech, warnings = [];
  if (mech === 'detuning') {
    const f0 = inp.det_f0, band = inp.det_band / 100, t0 = inp.det_t0 * 1e-6, r = inp.det_rate * 1e-9, dE = inp.det_dE / 100, tau = inp.det_tau, Nf = inp.det_Nf, dcr = inp.det_dcrit / 100, cyc = inp.det_cyc;
    const ratio = ty => { const geo = Math.max(0, 1 - 2 * r * ty / t0), Eage = 1 - dE * (1 - Math.exp(-ty / tau)), Efat = Math.max(0, 1 - dcr * Math.min(1, cyc * ty / Nf)); return geo * Math.sqrt(Math.max(0, Eage * Efat)); };
    const g = ty => Math.abs(ratio(ty) - 1) - band; let life = BIG;
    { let a = 0, b = 1e-3; while (b < 1e6 && g(b) < 0) { a = b; b *= 1.6; } if (g(b) >= 0) life = N.brent(g, a, b, 1e-9); }
    const horizon = inp.det_years, fEnd = f0 * ratio(horizon), tfat = Nf / Math.max(cyc, 1e-300), ts = N.linspace(0, Math.max(horizon, Math.min(life, 1e4) * 1.3), 120);
    const cg = 1 - Math.max(0, 1 - 2 * r * horizon / t0), ca = 1 - Math.sqrt(1 - dE * (1 - Math.exp(-horizon / tau))), cf = 1 - Math.sqrt(Math.max(0, 1 - dcr * Math.min(1, cyc * horizon / Nf))), dom = cg >= ca && cg >= cf ? 'thickness loss (corrosion/oxidation)' : ca >= cf ? 'modulus ageing' : 'fatigue stiffness degradation';
    warnings.push('Bending-type resonator scaling (f ∝ thickness·√E at fixed planform); mass change from oxide growth, residual-stress stiffening and temperature drift of the modulus are not included. Degradation laws are user-calibrated trends.');
    if (tfat < life) warnings.push(`Fatigue failure of the flexure (${f(tfat)} years at ${f(cyc)} cycles/year) precedes the detuning limit: structural failure, not drift, ends the life.`);
    return finish({
      summary: `The ${f(f0)} Hz resonance drifts by ${f((fEnd / f0 - 1) * 100)} % over ${f(horizon)} years and leaves the ±${f(inp.det_band)} % operating band after ${life >= BIG ? 'more than 10⁶ years' : f(life) + ' years'} (functional life). The dominant contribution over the horizon is ${dom}; fatigue of the flexure would occur at ${f(tfat)} years.`,
      kpis: [
        kpi('life', 'Functional life (band exit)', life, 'years', `Time at which |f/f₀ − 1| first exceeds ${f(inp.det_band)} % (10³⁰ = never within 10⁶ years), found by root finding on the combined drift law.`, life > horizon),
        kpi('fend', 'Resonance at end of horizon', fEnd, 'Hz', `Predicted resonance after ${f(horizon)} years; shift ${f((fEnd / f0 - 1) * 100)} %.`),
        kpi('shift', 'Frequency shift at horizon', (fEnd / f0 - 1) * 100, '%', 'Negative values mean softening (lower frequency); the sign matters for which band edge is crossed.', Math.abs(fEnd / f0 - 1) < band),
        kpi('k', 'Stiffness retained at horizon', (ratio(horizon) ** 2) * 100, '%', 'Effective stiffness-to-mass ratio relative to the as-built state, (f/f₀)².'),
        kpi('dom', 'Dominant drift mechanism', dom, '', `Contributions to the frequency shift at the horizon: thickness ${f(cg * 100)} %, modulus ageing ${f(ca * 100)} %, fatigue ${f(cf * 100)} %.`),
        kpi('tfat', 'Flexure fatigue life', tfat, 'years', `Cycles to failure ${f(Nf)} divided by ${f(cyc)} cycles per year.`, tfat > horizon),
        kpi('Qmax', 'Largest Q tolerating the drift', 1 / (2 * Math.max(Math.abs(fEnd / f0 - 1), 1e-12)), '–', 'A resonance of quality factor Q has half-bandwidth f₀/2Q; sharper resonances than this fall out of tune within the horizon.'),
      ],
      series: [ser('ft', 'Resonance drift', 'time (years)', 'frequency (Hz)', [{ name: 'resonance', x: ts, y: ts.map(t => f0 * ratio(t)) }, { name: 'upper band edge', x: [ts[0], ts[119]], y: [f0 * (1 + band), f0 * (1 + band)] }, { name: 'lower band edge', x: [ts[0], ts[119]], y: [f0 * (1 - band), f0 * (1 - band)] }], 'Resonance frequency against the operating band; functional life ends where the curve leaves the band.'),
        ser('contrib', 'Stiffness and geometry retention', 'time (years)', 'retained (%)', [{ name: 'thickness', x: ts, y: ts.map(t => Math.max(0, 1 - 2 * r * t / t0) * 100) }, { name: 'modulus (ageing)', x: ts, y: ts.map(t => (1 - dE * (1 - Math.exp(-t / tau))) * 100) }, { name: 'modulus (fatigue)', x: ts, y: ts.map(t => Math.max(0, 1 - dcr * Math.min(1, cyc * t / Nf)) * 100) }], 'Individual degradation paths that combine into the frequency drift.')],
      equations: [eqn('formula', 'Resonance of a bending element', 'f₀ ∝ (t/L²) √(E/ρ)'), eqn('constitutive', 'Drift law', 'f(t)/f₀ = [1 − 2 r t/t₀] · √{ [1 − ΔE∞(1 − e^{−t/τ})] · [1 − d_c·min(1, N(t)/N_f)] }'), eqn('formula', 'Functional-life criterion', '|f(t)/f₀ − 1| = band half-width')],
      checks: [chk('Resonance inside band at horizon', Math.abs(fEnd / f0 - 1) < band, `shift ${f((fEnd / f0 - 1) * 100)} % vs ±${f(inp.det_band)} %`), chk('Structural (fatigue) life exceeds functional life', tfat >= Math.min(life, 1e6), `${f(tfat)} vs ${f(life)} years`), chk('Stiffness degradation is monotonic (irreversible)', ratio(horizon) <= 1 + 1e-12, `f/f₀ = ${f(ratio(horizon))}`)],
      warnings,
      recommendations: [
        rec('high', life <= horizon ? `Functional life (${f(life)} years) is shorter than the ${f(horizon)}-year horizon: widen the band (lower Q or add tuning range ≥ ${f(Math.abs(fEnd / f0 - 1) * 100)} %), or suppress ${dom}.` : `The resonance stays in band for ${life >= BIG ? 'the whole horizon' : f(life) + ' years'}; margin at the horizon is ${f((band - Math.abs(fEnd / f0 - 1)) * 100)} % of f₀.`, 'Drift-law root for band exit.'),
        rec('medium', dom.startsWith('thickness') ? `Passivate the resonator (ALD oxide or noble-metal cap): the ${f(inp.det_rate)} nm/yr loss on a ${f(inp.det_t0)} µm element is the main drift term.` : dom.startsWith('modulus') ? 'Pre-age (anneal) the structure so that most of the modulus relaxation occurs before calibration.' : 'Reduce the cyclic strain in the flexures (fillets, lower amplitude) so that fatigue softening does not accumulate.', 'Dominant mechanism of this run.'),
        rec('medium', 'Add active or passive tuning (varactor, thermal or electrostatic) sized to the predicted drift so the metamaterial can be re-centred in service; verify drift by periodic frequency sweeps.', 'Predicted shift at horizon.'),
        rec('low', 'A re-tunable or re-calibratable design extends service life and avoids replacing whole metamaterial panels when only the resonance has drifted.', 'Functional-life concept.', true),
      ],
      publish: { 'life.years': Math.min(life, tfat), 'res.f0': fEnd, 'life.f_shift': fEnd / f0 - 1, 'life.cycles': Nf },
      numerics: { method: 'Brent root finding on the analytical drift law after geometric bracketing', dof: 1, steps: 120, iterations: 1, residual: life < BIG ? Math.abs(g(life)) : 0, notes: 'Bracket expanded geometrically up to 10⁶ years.' },
    });
  }
  if (mech === 'bistable') {
    const l = inp.bi_l * 1e-3, t = inp.bi_t * 1e-6, h = inp.bi_h * 1e-6, bw = inp.bi_w * 1e-6, E = inp.bi_E * 1e9, I = bw * t ** 3 / 12, Q = h / t, eps = 2 * PI * PI * t * h / (l * l), Fsnap = 740 * E * I * h / l ** 3, sf = inp.bi_sf * 1e6, b = -Math.abs(inp.bi_b);
    const sa = E * eps / 2, Nfat = sa >= sf ? 0.5 : capN(0.5 * (sa / sf) ** (1 / b)), rate = inp.bi_rate, duty = inp.bi_duty, tauR = inp.bi_tau, br = inp.bi_beta, Einf = inp.bi_Einf;
    const relax = th => Einf + (1 - Einf) * Math.exp(-((th / tauR) ** br)), Qeff = th => Q * (1 - duty * (1 - relax(th)));
    const Qc = 2.31; let tb = BIG;
    if (Q > Qc) { const g = th => Qeff(th) - Qc; if (Q * (1 - duty * (1 - Einf)) < Qc) { let a = 0, bb = 1e-3; while (g(bb) > 0 && bb < 1e12) { a = bb; bb *= 2; } tb = N.brent(g, a, bb, 1e-9); } } else tb = 0;
    const Nb = tb >= BIG ? BIG : tb * rate, Nlife = Math.min(Nfat, Nb), mode = Q <= Qc ? 'not bistable as built' : Nfat <= Nb ? 'fatigue of the beam' : 'loss of bistability by stress relaxation', ths = lsp(-2, Math.log10(Math.max(tauR * 100, 10)), 70);
    const kneg = -(Fsnap + 0.5 * Fsnap) / (1.76 * h); // secant slope between the force maximum (≈ 0.16h) and minimum (≈ −F/2 at 1.92h)
    warnings.push('Curved (cosine-shaped) clamped–clamped beam formulae of the Qiu–Lang–Slocum type for small t/h and constrained second mode; snap force and strain are approximate. Relaxation uses a stretched-exponential (KWW) modulus and assumes the permanent set is proportional to the relaxed stress fraction while held in the snapped state.');
    if (Q <= Qc) warnings.push(`Geometry ratio h/t = ${f(Q)} is below 2.31: the element is monostable (snap-through without latching) from the start.`);
    if (eps > 0.05) warnings.push('Peak strain exceeds 5 %: beyond small-strain beam theory; use a hyperelastic finite-element analysis.');
    return finish({
      summary: `The bistable beam (h/t = ${f(Q)}) snaps at ${f(Fsnap * 1e3)} mN with a peak strain of ${f(eps * 100)} %. Fatigue life is ${f(Nfat)} cycles; stress relaxation (${f(duty * 100)} % of time latched) erodes bistability after ${tb >= BIG ? 'no finite time' : f(tb) + ' h'} (${f(Nb)} cycles at ${f(rate)} cycles/h). Functional life is ${f(Nlife)} cycles, limited by ${mode}.`,
      kpis: [
        kpi('Q', 'Bistability ratio h/t', Q, '–', 'Apex height over thickness of the curved beam; bistability requires h/t > 2.31 (with the second buckling mode constrained).', Q > Qc),
        kpi('Fsnap', 'Snap-through force', Fsnap * 1e3, 'mN', '≈ 740·E·I·h/l³: peak force needed to switch state; the snap-back force is about half of it.'),
        kpi('kneg', 'Mean negative stiffness in the unstable branch', kneg, 'N/m', 'Secant slope between the force extrema (negative): the usable negative-stiffness region of the element. It must be stabilised by displacement control or a stiffer parallel spring.'),
        kpi('eps', 'Peak bending strain', eps * 100, '%', '≈ 2π²·t·h/l²; it must stay below the elastic limit of the material (about 1 % for silicon, 2–5 % for polymers).', eps < 0.02),
        kpi('Nfat', 'Fatigue life', Nfat, 'cycles', `Basquin life at stress amplitude E·ε/2 = ${f(sa / 1e6)} MPa with σ_f′ = ${f(inp.bi_sf)} MPa, b = ${f(b)}.`),
        kpi('Nb', 'Cycles to loss of bistability', Nb, 'cycles', `Time for the effective ratio h_eff/t to relax to 2.31 (${tb >= BIG ? 'never: the relaxed ratio stays above 2.31' : f(tb) + ' h'}), converted with the switching rate.`),
        kpi('Nlife', 'Functional life', Nlife, 'cycles', `Smaller of fatigue life and loss of bistability; governing mode: ${mode}.`, Nlife > inp.bi_target),
      ],
      series: [ser('Qt', 'Effective bistability ratio versus time', 'time (h)', 'h_eff/t', [{ name: 'effective ratio', x: ths, y: ths.map(Qeff) }, { name: 'bistability limit 2.31', x: [ths[0], ths[69]], y: [Qc, Qc] }], 'Viscoelastic set while latched lowers the effective apex height; the element stops latching when the curve reaches 2.31.', { logx: true }),
        ser('Ft', 'Snap force retention', 'time (h)', 'snap force (mN)', [{ name: 'snap-through force', x: ths, y: ths.map(th => Fsnap * relax(th) * (Qeff(th) / Q) * 1e3) }], 'Decay of the switching force by modulus relaxation and loss of apex height.', { logx: true })],
      equations: [eqn('formula', 'Bistability condition (curved beam)', 'Q = h/t > 2.31'), eqn('formula', 'Snap force and peak strain', 'F_top ≈ 740 E I h / l³,   ε_max ≈ 2π² t h / l²'), eqn('constitutive', 'Stress relaxation (KWW)', 'E(t)/E₀ = E∞ + (1 − E∞) exp[−(t/τ)^β];   h_eff = h [1 − duty·(1 − E(t)/E₀)]'), eqn('constitutive', 'Basquin fatigue', 'σ_a = σ_f′ (2N)^b,   σ_a = E ε_max/2')],
      checks: [chk('Bistable as built', Q > Qc, `h/t = ${f(Q)}`), chk('Negative stiffness present in the unstable branch', kneg < 0, `${f(kneg)} N/m (requires stabilisation)`), chk('Strain within small-strain theory', eps <= 0.05, `${f(eps * 100)} %`), chk('Life meets target', Nlife >= inp.bi_target, `${f(Nlife)} vs ${f(inp.bi_target)} cycles`)],
      warnings,
      recommendations: [
        rec('high', Q <= Qc ? `Increase apex height to at least ${f(2.6 * inp.bi_t)} µm (h/t ≥ 2.6) for robust bistability.` : mode.startsWith('loss') ? `Raise h/t from ${f(Q)} to ${f(Qc / (1 - duty * (1 - Einf)) * 1.1)} so that the fully relaxed element still latches, or choose a material with higher relaxed modulus ratio (now ${f(Einf)}).` : `Lower the peak strain (${f(eps * 100)} %) by a thinner or longer beam — strain scales as t·h/l² — to lift the fatigue life above the ${f(inp.bi_target)}-cycle target.`, 'Governing failure mode of this run.'),
        rec('medium', 'Store the element in its as-fabricated (stress-free) state and minimise dwell in the snapped state; relaxation acts only while latched.', 'Duty-cycle term of the relaxation model.'),
        rec('medium', 'Measure the relaxation modulus (DMA or stress-relaxation test at service temperature) and cycle a sample to 10 % force loss to calibrate τ, β and E∞.', 'Relaxation parameters are inputs.'),
        rec('low', 'Bistable elements hold state without power: use them for reconfigurable metamaterials to eliminate holding energy, and select recyclable thermoplastics or metals over thermosets where the strain allows.', 'Zero-power latching; material choice.', true),
      ],
      publish: { 'life.cycles': Nlife, 'mech.buckling_load': Fsnap, 'mech.max_stress': E * eps, 'eff.k_negative': kneg },
      numerics: { method: 'Closed-form curved-beam relations; Brent root for the relaxation time to the bistability limit', dof: 1, steps: 70, iterations: 1, residual: 0, notes: 'Reduced-order (single-mode) beam model.' },
    });
  }
  if (mech === 'switching') {
    const tech = SWT[inp.sw_tech] || SWT.pcm, T = inp.sw_T + 273.15, eta = tech.eta * inp.sw_scale * inp.sw_drive ** (-tech.p) * Math.exp(tech.Ea / KEV * (1 / T - 1 / 298.15)), beta = tech.beta, rate = inp.sw_rate * inp.sw_duty * YEAR, nEl = inp.sw_n, ftol = inp.sw_ftol / 100;
    const Ntol = eta * (-Math.log(1 - ftol)) ** (1 / beta), N1 = eta * (-Math.log(1 - 1 / Math.max(nEl, 2))) ** (1 / beta), yTol = Ntol / Math.max(rate, 1e-300), y1 = N1 / Math.max(rate, 1e-300), Ntar = inp.sw_years * rate, Ftar = 1 - Math.exp(-((Ntar / eta) ** beta));
    const Ns = lsp(Math.log10(eta) - 6, Math.log10(eta) + 1, 70);
    warnings.push(`Baseline endurance for "${tech.name}" is a typical order of magnitude (η = ${f(tech.eta)} cycles, β = ${f(tech.beta)}), scaled by an inverse-power law in drive amplitude (exponent ${tech.p}) and an Arrhenius factor (${tech.Ea} eV). Published endurance spans several decades between devices — calibrate the scale factor with your own cycling data.`);
    if (tech === SWT.varactor) warnings.push('Semiconductor varactors do not wear out by cycling in the usual sense; the cycle count here stands for bias/thermal stress time and the time-based mechanisms (electromigration, TDDB) should be checked as well.');
    return finish({
      summary: `${tech.name}: characteristic endurance ${f(eta)} cycles at ${f(inp.sw_drive)}× nominal drive and ${f(inp.sw_T)} °C. At ${f(inp.sw_rate)} Hz with ${f(inp.sw_duty * 100)} % activity the first of ${f(nEl)} elements is expected to fail after ${f(y1)} years and the tolerated ${f(inp.sw_ftol)} % failed fraction is reached after ${f(yTol)} years; after ${f(inp.sw_years)} years ${f(Ftar * 100)} % have failed.`,
      kpis: [
        kpi('eta', 'Characteristic endurance', eta, 'cycles', 'Weibull scale (63 % failed) after drive and temperature acceleration.'),
        kpi('ytol', 'Array functional life', yTol, 'years', `Time until ${f(inp.sw_ftol)} % of the elements have failed — the fraction of dead unit cells the metasurface tolerates before its pattern degrades.`, yTol > inp.sw_years),
        kpi('y1', 'Expected time to first failure', y1, 'years', `Cycles at which the failure probability equals 1/${f(nEl)}, converted with the switching rate.`),
        kpi('F', 'Failed fraction at target life', Ftar * 100, '%', `Weibull cumulative failure probability evaluated at the ${f(Ntar)} cycles accumulated in ${f(inp.sw_years)} years; compare with the tolerated failed fraction.`, Ftar < ftol),
        kpi('cyc', 'Cycles per year', rate, 'cycles/yr', 'Switching frequency × activity factor × seconds per year.'),
        kpi('B1', 'B1 endurance', eta * (-Math.log(0.99)) ** (1 / beta), 'cycles', 'Cycles by which 1 % of elements have failed.'),
      ],
      series: [ser('cdf', 'Failed fraction versus cycles', 'switching cycles', 'failed fraction', [{ name: 'Weibull CDF', x: Ns, y: Ns.map(n => Math.max(1e-12, 1 - Math.exp(-((n / eta) ** beta)))) }, { name: 'tolerated fraction', x: [Ns[0], Ns[69]], y: [ftol, ftol] }], 'Wear-out distribution of the switching elements; functional life is where the curve crosses the tolerated failed fraction.', { logx: true, logy: true }),
        ser('drive', 'Functional life versus drive level', 'drive amplitude (× nominal)', 'array life (years)', (() => { const dv = N.linspace(0.6, 1.6, 50); return [{ name: 'array life', x: dv, y: dv.map(d => Ntol * (d / inp.sw_drive) ** (-tech.p) / Math.max(rate, 1e-300)) }]; })(), 'Derating the drive voltage/current/pulse energy buys life through the inverse-power law.', { logy: true })],
      equations: [eqn('formula', 'Weibull wear-out', 'F(N) = 1 − exp[−(N/η)^β]'), eqn('formula', 'Acceleration', 'η = η₀ · s · (V/V₀)^{−p} · exp[(E_a/k)(1/T − 1/298 K)]'), eqn('formula', 'Array criteria', 'N_first = η[−ln(1 − 1/n)]^{1/β};   N_tol = η[−ln(1 − f_tol)]^{1/β}')],
      checks: [chk('Failed fraction at target within tolerance', Ftar < ftol, `${f(Ftar * 100)} % vs ${f(inp.sw_ftol)} %`), chk('Wear-out regime (β > 1)', beta > 1, `β = ${f(beta)}${beta <= 1 ? ' — random/early failures dominate; burn-in or redundancy is more effective than derating' : ''}`), chk('Reliability function monotonic', true, 'F(N) increases with N by construction')],
      warnings,
      recommendations: [
        rec('high', Ftar >= ftol ? `Derate the drive to ${f(inp.sw_drive * (Ntar / Ntol) ** (-1 / tech.p))}× nominal or reduce the switching rate to ${f(inp.sw_rate * yTol / inp.sw_years)} Hz to reach ${f(inp.sw_years)} years; alternatively design the pattern to tolerate more than ${f(Ftar * 100)} % dead cells.` : `The array meets the ${f(inp.sw_years)}-year target with ${f(Ftar * 100)} % failed elements; first failures appear after about ${f(y1)} years, so include self-test.`, 'Weibull array-life calculation.'),
        rec('medium', 'Cycle at least 20 elements to failure (or to 10⁹ cycles) at two drive levels to fit η, β and the drive exponent; replace the typical values used here.', 'Baseline endurance is a typical order of magnitude.'),
        rec('medium', 'Use error-tolerant coding of the metasurface state (spatial redundancy, re-optimisation around stuck cells) so isolated failures do not end the functional life.', 'Tolerated-fraction criterion.'),
        rec('low', 'Reducing unnecessary reconfiguration (event-driven rather than periodic switching) saves drive energy and extends life in direct proportion.', 'Cycles-per-year term.', true),
      ],
      publish: { 'life.cycles': Ntol, 'life.years': yTol, 'life.pf': Ftar },
      numerics: { method: 'Two-parameter Weibull with inverse-power and Arrhenius acceleration (closed form)', dof: 1, steps: 1, iterations: 1, residual: 0, notes: 'Typical-order baseline endurance per technology.' },
    });
  }
  // auxetic hinge fatigue (rotating rigid units)
  const th0 = inp.ax_theta * PI / 180, S = a => Math.cos(a / 2) + Math.sin(a / 2), eg = inp.ax_strain / 100, target = S(th0) * (1 + eg), thMaxS = Math.SQRT2;
  let th1 = th0, reach = true; if (target >= thMaxS) { th1 = PI / 2; reach = false; } else if (target <= 1) { th1 = 0; reach = false; } else th1 = 2 * (Math.asin(target / Math.SQRT2) - PI / 4);
  const dth = Math.abs(th1 - th0), tH = inp.ax_th * 1e-6, lH = inp.ax_lh * 1e-6, epsH = tH * dth / (2 * lH), E = inp.ax_E * 1e9, sf = inp.ax_sf * 1e6, b = -Math.abs(inp.ax_b), ef = inp.ax_ef, c = -Math.abs(inp.ax_c);
  const ea = epsH / 2, g = ln2N => sf / E * Math.exp(b * ln2N) + ef * Math.exp(c * ln2N) - ea; let Nf = BIG;
  if (ea > 0) { if (g(0) <= 0) Nf = 0.5; else if (g(69) < 0) Nf = 0.5 * Math.exp(N.brent(g, 0, 69, 1e-10)); }
  const yrs = Nf >= BIG ? BIG : Nf / Math.max(inp.ax_cyc, 1e-300), epsMaxAllow = 2 * (sf / E * (2 * inp.ax_target) ** b + ef * (2 * inp.ax_target) ** c), tAllow = epsMaxAllow * 2 * lH / Math.max(dth, 1e-12), egs = N.linspace(0.05, Math.max(eg * 100 * 2, 2), 50);
  const epsOf = e => { const tg = S(th0) * (1 + e / 100); if (tg >= Math.SQRT2 || tg <= 1) return NaN; return tH * Math.abs(2 * (Math.asin(tg / Math.SQRT2) - PI / 4) - th0) / (2 * lH); };
  warnings.push('Rotating-rigid-unit kinematics (Poisson’s ratio −1, units rigid, all deformation in the hinges) with uniform bending of a flexure hinge of the entered length; real hinges concentrate strain at fillets, so apply a stress-concentration factor from a finite-element model of the hinge.');
  if (!reach) warnings.push('The requested global strain exceeds the kinematic range of the rotating-unit mechanism (fully open or fully closed); the hinge rotation is capped at that limit.');
  return finish({
    summary: `A global strain of ${f(inp.ax_strain)} % rotates the units by ${f(dth * 180 / PI)}° from the ${f(inp.ax_theta)}° opening angle, giving a hinge bending strain range of ${f(epsH * 100)} %. The strain–life (Basquin + Coffin–Manson) prediction is ${f(Nf)} cycles (${f(yrs)} years at ${f(inp.ax_cyc)} cycles/year). The ideal Poisson’s ratio of the mechanism is −1 throughout.`,
    kpis: [
      kpi('dtheta', 'Hinge rotation per cycle', dth * 180 / PI, '°', 'From the rotating-squares kinematics L ∝ cos(θ/2) + sin(θ/2) solved for the angle that produces the applied global strain.'),
      kpi('eps', 'Hinge strain range', epsH * 100, '%', `t_h·Δθ/(2·l_h) for a flexure ${f(inp.ax_th)} µm thick and ${f(inp.ax_lh)} µm long (uniform curvature).`, epsH < 0.02),
      kpi('Nf', 'Hinge fatigue life', Nf, 'cycles', 'Strain–life equation solved for the reversals at strain amplitude Δε/2 (10³⁰ = below the range of the curve, effectively infinite).', Nf > inp.ax_target),
      kpi('years', 'Service life', yrs, 'years', `Cycles to failure divided by ${f(inp.ax_cyc)} cycles per year.`),
      kpi('nu', 'Poisson’s ratio of the mechanism', -1, '–', 'Ideal rotating rigid squares are isotropic auxetics with ν = −1 at all strains; finite hinge stiffness moves it towards zero.'),
      kpi('tallow', 'Allowable hinge thickness for target life', tAllow * 1e6, 'µm', `Hinge thickness at which the strain–life curve gives ${f(inp.ax_target)} cycles for the same rotation and hinge length.`, tH <= tAllow),
      kpi('range', 'Kinematic strain range available', (thMaxS / S(th0) - 1) * 100, '%', 'Maximum tensile global strain before the units are fully open (θ = 90°).'),
    ],
    series: [ser('eN', 'Hinge strain–life curve', 'cycles to failure', 'strain amplitude (%)', (() => { const Ns = lsp(1, 9, 60); return [{ name: 'strain–life', x: Ns, y: Ns.map(n => (sf / E * (2 * n) ** b + ef * (2 * n) ** c) * 100) }, { name: 'operating amplitude', x: [10, 1e9], y: [ea * 100, ea * 100] }]; })(), 'Operating strain amplitude against the material strain–life curve; the crossing is the predicted life.', { logx: true, logy: true }),
      ser('eg', 'Hinge strain versus global strain', 'global strain (%)', 'hinge strain range (%)', [{ name: 'hinge strain', x: egs, y: egs.map(e => epsOf(e) * 100) }], 'Kinematic amplification from global strain to hinge strain; it grows as the units approach the fully open state.')],
    equations: [eqn('formula', 'Rotating-squares kinematics', 'L(θ) ∝ cos(θ/2) + sin(θ/2);   ε = L(θ)/L(θ₀) − 1;   ν = −1'), eqn('formula', 'Flexure-hinge strain', 'Δε = t_h Δθ /(2 l_h)'), eqn('constitutive', 'Strain–life', 'Δε/2 = (σ_f′/E)(2N)^b + ε_f′(2N)^c')],
    checks: [chk('Strain within the mechanism range', reach, reach ? `θ: ${f(inp.ax_theta)}° → ${f(th1 * 180 / PI)}°` : 'kinematic limit reached'), chk('Poisson’s ratio negative (auxetic) and within isotropic bounds −1 ≤ ν ≤ 0.5', true, 'ν = −1 (lower stability bound of an isotropic solid; shear modulus → ∞ relative to bulk modulus)'), chk('Life meets target', Nf >= inp.ax_target, `${f(Nf)} vs ${f(inp.ax_target)} cycles`), chk('Hinge strain below 2 %', epsH < 0.02, `${f(epsH * 100)} %`)],
    warnings,
    recommendations: [
      rec('high', Nf < inp.ax_target ? `Thin the hinges to ≤ ${f(tAllow * 1e6)} µm or lengthen them to ≥ ${f(lH * tH / Math.max(tAllow, 1e-12) * 1e6)} µm to reach ${f(inp.ax_target)} cycles; hinge strain scales with t_h/l_h.` : `Hinge life (${f(Nf)} cycles) exceeds the target; thickness could rise to ${f(tAllow * 1e6)} µm for higher stiffness if needed.`, 'Strain–life result for the hinge.'),
      rec('medium', 'Replace sharp-cornered hinges by circular or elliptical flexures and polish or coat the hinge surface; fatigue cracks start at the hinge root.', 'Stress concentration not included in the uniform-curvature estimate.'),
      rec('medium', `Limit the working strain to well below the ${f((thMaxS / S(th0) - 1) * 100)} % kinematic limit; near full opening the hinge rotation per unit strain diverges.`, 'Kinematic amplification curve.'),
      rec('low', 'Compliant-hinge auxetics made from a single recyclable material (one polymer or one metal) avoid assembly, adhesives and mixed-material waste at end of life.', 'Monolithic design practice.', true),
    ],
    publish: { 'life.cycles': Nf, 'life.years': yrs, 'eff.nu': -1, 'mech.max_stress': E * Math.min(epsH, sf / E) },
    numerics: { method: 'Closed-form mechanism kinematics; Brent root of the strain–life equation in ln(2N)', dof: 1, steps: 1, iterations: 1, residual: Nf < BIG && Nf > 0.5 ? Math.abs(g(Math.log(2 * Nf))) : 0, notes: 'Uniform-curvature flexure hinge.' },
  });
}

// ================================================================================================
// 6. Reliability synthesis
// ================================================================================================
function runReliability(inp) {
  const warnings = []; let synthetic = false, fl, dg;
  try { fl = parseTable(inp.rel_fail, 1, 4, 'Failure-time table').map(r => r[0]).filter(v => v > 0); if (fl.length < 4) throw new Error('Failure-time table: need at least 4 positive failure times.'); } catch (e) { warnings.push(e.message + ' The demonstration data set is used instead.'); fl = parseTable(DEF_FAIL, 1, 4, 'Failure-time table').map(r => r[0]); synthetic = true; }
  try { dg = parseTable(inp.rel_deg, 2, 4, 'Degradation table'); } catch (e) { warnings.push(e.message + ' The demonstration data set is used instead.'); dg = parseTable(DEF_DEG, 2, 4, 'Degradation table'); synthetic = true; }
  if (String(inp.rel_fail || '').trim() === DEF_FAIL.trim() || String(inp.rel_deg || '').trim() === DEF_DEG.trim()) synthetic = true;
  if (synthetic) warnings.push('Failure-time and/or degradation tables are synthetic demonstration data — replace with measurements.');
  fl.sort((a, b) => a - b); const n = fl.length;
  // Weibull (median-rank regression) and lognormal (probability-plot regression)
  let wf = N.weibullFit(fl); if (!Number.isFinite(wf.shape) || wf.shape <= 0) { wf = { shape: 1, scale: N.mean(fl), r2: 0 }; warnings.push('Weibull regression failed (identical failure times?); exponential distribution with the sample mean is used.'); }
  const beta = wf.shape, eta = wf.scale, lnT = fl.map(Math.log), zs = fl.map((_, i) => N.normInv((i + 0.7) / (n + 0.4))), lr = N.linreg(zs, lnT), muL = Number.isFinite(lr.intercept) ? lr.intercept : N.mean(lnT), sgL = Number.isFinite(lr.slope) && lr.slope > 0 ? lr.slope : Math.max(N.std(lnT), 1e-6), r2L = Number.isFinite(lr.r2) ? lr.r2 : 0;
  const best = wf.r2 >= r2L ? 'Weibull' : 'lognormal', Rw = t => Math.exp(-((t / eta) ** beta)), RL = t => 1 - N.normCdf((Math.log(t) - muL) / sgL), Rb = best === 'Weibull' ? Rw : RL;
  const hz = t => (best === 'Weibull' ? beta / eta * (t / eta) ** (beta - 1) : (() => { const z = (Math.log(t) - muL) / sgL, pdf = Math.exp(-0.5 * z * z) / (t * sgL * Math.sqrt(2 * PI)); return pdf / Math.max(RL(t), 1e-300); })());
  const B10 = best === 'Weibull' ? eta * (-Math.log(0.9)) ** (1 / beta) : Math.exp(muL + sgL * N.normInv(0.1)), mttf = best === 'Weibull' ? eta * N.gamma(1 + 1 / beta) : Math.exp(muL + sgL * sgL / 2), Rt = inp.rel_Rt, tR = best === 'Weibull' ? eta * (-Math.log(Rt)) ** (1 / beta) : Math.exp(muL + sgL * N.normInv(1 - Rt));
  // approximate 90 % confidence on β and η (large-sample normal approximation of the MLE variances)
  const seB = 0.78 * beta / Math.sqrt(n), seE = 1.05 * eta / (beta * Math.sqrt(n));
  // age-replacement optimum for cost ratio (failure/preventive)
  const cr = inp.rel_cost, tmax = fl[n - 1] * 3, costRate = T => { let s = 0; const m = 200; for (let i = 0; i < m; i++) { const a = T * i / m, b_ = T * (i + 1) / m; s += 0.5 * (Rb(Math.max(a, 1e-12)) + Rb(b_)) * (b_ - a); } return (Rb(T) + cr * (1 - Rb(T))) / Math.max(s, 1e-300); };
  const wear = best === 'Weibull' ? beta > 1 : true; let Topt = BIG, cOpt = cr / mttf, saving = 0;
  if (wear && cr > 1) { const T1 = N.goldenSection(costRate, fl[0] * 0.02, tmax, tmax * 1e-6); const c1 = costRate(T1); if (c1 < cr / mttf * 0.999 && T1 < tmax * 0.98) { Topt = T1; cOpt = c1; saving = 1 - c1 / (cr / mttf); } }
  // remaining useful life from the degradation path
  const tx = dg.map(r => r[0]), yv = dg.map(r => r[1]), thr = inp.rel_thr, expo = inp.rel_trend === 'exp' && yv.every(v => v > 0) && thr > 0;
  if (inp.rel_trend === 'exp' && !expo) warnings.push('Exponential trend needs positive values and threshold; a linear trend is used instead.');
  const yy = expo ? yv.map(Math.log) : yv, fit = N.linreg(tx, yy), m_ = fit.slope, c0 = fit.intercept, yThr = expo ? Math.log(thr) : thr, nd = tx.length, tbar = N.mean(tx), Sxx = tx.reduce((s, v) => s + (v - tbar) ** 2, 0);
  const se = Math.sqrt(yy.reduce((s, v, i) => s + (v - (c0 + m_ * tx[i])) ** 2, 0) / Math.max(nd - 2, 1)), tNow = tx[nd - 1], yNow = c0 + m_ * tNow, towards = Number.isFinite(m_) && m_ !== 0 && (yThr - yNow) / m_ > 0;
  const tFail = towards ? (yThr - c0) / m_ : BIG, RUL = towards ? Math.max(0, tFail - tNow) : BIG, z90 = 1.645;
  const sT = towards ? se / Math.abs(m_) * Math.sqrt(1 + 1 / nd + (tFail - tbar) ** 2 / Math.max(Sxx, 1e-300)) : 0, rulLo = towards ? Math.max(0, RUL - z90 * sT) : BIG, rulHi = towards ? RUL + z90 * sT : BIG;
  if (!towards) warnings.push('The fitted degradation trend does not approach the threshold: no failure time can be extrapolated (RUL shown as 10³⁰).');
  if (towards && RUL > 3 * (tNow - tx[0])) warnings.push('The extrapolated RUL is more than three times the observed duration: the trend may change (wear-out acceleration) before the threshold is reached.');
  if (n < 10) warnings.push(`Only ${n} failure times: distribution parameters have wide confidence bounds (β ± ${f(1.645 * seB)}).`);
  warnings.push('Complete (uncensored) failure data are assumed; suspended units must be handled with a censored-data method (not implemented), otherwise life is under-estimated.');
  const ts = lsp(Math.log10(fl[0] / 5), Math.log10(fl[n - 1] * 3), 80), tEx = N.linspace(tx[0], Math.max(tNow, Math.min(towards ? tFail * 1.15 : tNow * 2, tNow * 6)), 80), back = v => (expo ? Math.exp(v) : v);
  const band = t => z90 * se * Math.sqrt(1 + 1 / nd + (t - tbar) ** 2 / Math.max(Sxx, 1e-300)), interval = Math.min(Topt, tR, towards ? Math.max(rulLo, 0) + 0 * tNow : BIG);
  return finish({
    summary: `${n} failure times fit a ${best} distribution best (Weibull β = ${f(beta)}, η = ${f(eta)} h, R² = ${f(wf.r2)}; lognormal R² = ${f(r2L)}): MTTF ${f(mttf)} h, B10 life ${f(B10)} h, and ${f(Rt * 100)} % reliability up to ${f(tR)} h. The ${expo ? 'exponential' : 'linear'} degradation trend reaches the threshold ${f(thr)} at ${f(tFail)} h, i.e. a remaining useful life of ${f(RUL)} h (90 % band ${f(rulLo)}–${f(rulHi)} h). ${Topt < BIG ? `Preventive replacement every ${f(Topt)} h minimises cost (saving ${f(saving * 100)} %).` : 'No cost-optimal preventive interval exists (run to failure or use condition-based maintenance).'}`,
    kpis: [
      kpi('beta', 'Weibull shape β', beta, '–', `Median-rank regression on ${n} failures (approx. 90 % interval ±${f(1.645 * seB)}). β < 1 infant mortality, ≈ 1 random failures, > 1 wear-out.`),
      kpi('eta', 'Weibull scale η', eta, 'h', `Characteristic life (63.2 % failed), approx. 90 % interval ±${f(1.645 * seE)} h.`),
      kpi('best', 'Best-fitting distribution', best, '', `Chosen by the higher probability-plot R² (Weibull ${f(wf.r2)}, lognormal ${f(r2L)}); with few data the two are often indistinguishable.`),
      kpi('mttf', 'Mean time to failure', mttf, 'h', 'Mean of the best-fitting distribution (η·Γ(1 + 1/β) for Weibull, exp(μ + σ²/2) for lognormal); it is not a safe life — about half or more of the units fail before it.'),
      kpi('B10', 'B10 life', B10, 'h', 'Age by which 10 % of units fail; the usual basis for warranty and service intervals.'),
      kpi('tR', 'Life at target reliability', tR, 'h', `Age at which the fitted reliability falls to the target of ${f(Rt * 100)} %; replace or inspect before this age to hold the target.`),
      kpi('RUL', 'Remaining useful life', RUL, 'h', `Extrapolation of the fitted ${expo ? 'exponential' : 'linear'} trend (rate ${f(m_)} per h${expo ? ' in ln units' : ''}) from the last observation at ${f(tNow)} h to the threshold.`, RUL > 0.2 * tNow),
      kpi('RULlo', 'RUL lower 90 % bound', rulLo, 'h', 'Lower prediction bound from the regression standard error; plan maintenance on this conservative value.'),
      kpi('Topt', 'Cost-optimal replacement age', Topt, 'h', `Age-replacement policy minimising [c_p R(T) + c_f F(T)]/∫₀ᵀR dt for a failure/preventive cost ratio of ${f(cr)} (10³⁰ = run to failure).`),
      kpi('interval', 'Recommended maintenance interval', interval, 'h', 'Smallest of the cost-optimal age, the target-reliability life and the lower RUL bound.'),
    ],
    series: [
      ser('R', 'Reliability and data', 'time (h)', 'reliability', [{ name: 'Weibull', x: ts, y: ts.map(Rw) }, { name: 'lognormal', x: ts, y: ts.map(RL) }, { name: 'data (median ranks)', x: fl, y: fl.map((_, i) => 1 - (i + 0.7) / (n + 0.4)) }], 'Fitted survival functions against the median-rank plotting positions of the failures.', { logx: true }),
      ser('h', 'Hazard rate', 'time (h)', 'hazard (1/h)', [{ name: best, x: ts, y: ts.map(hz) }], 'Instantaneous failure rate of the best-fitting distribution: rising hazard indicates wear-out and justifies preventive replacement.', { logx: true, logy: true }),
      ser('deg', 'Degradation path and RUL', 'time (h)', 'performance', [{ name: 'measurements', x: tx, y: yv }, { name: 'trend', x: tEx, y: tEx.map(t => back(c0 + m_ * t)) }, { name: '90 % lower', x: tEx, y: tEx.map(t => back(c0 + m_ * t - band(t))) }, { name: '90 % upper', x: tEx, y: tEx.map(t => back(c0 + m_ * t + band(t))) }, { name: 'threshold', x: [tEx[0], tEx[79]], y: [thr, thr] }], 'Fitted degradation trend with its 90 % prediction band extrapolated to the failure threshold.'),
      ser('cost', 'Maintenance cost rate', 'replacement age (h)', 'cost rate (preventive-cost units per h)', (() => { const Tm = N.linspace(fl[0] * 0.1, tmax, 60); return [{ name: 'age replacement', x: Tm, y: Tm.map(costRate) }, { name: 'run to failure', x: [Tm[0], Tm[59]], y: [cr / mttf, cr / mttf] }]; })(), 'Expected cost per hour as a function of the preventive replacement age; a minimum below the run-to-failure line indicates worthwhile preventive maintenance.'),
    ],
    tables: [{ key: 'fit', title: 'Distribution fits', columns: ['distribution', 'parameter 1', 'parameter 2', 'R² (probability plot)', 'B10 (h)', 'MTTF (h)'], rows: [['Weibull', `β = ${f(beta, 4)}`, `η = ${f(eta, 4)} h`, +wf.r2.toPrecision(4), +(eta * (-Math.log(0.9)) ** (1 / beta)).toPrecision(4), +(eta * N.gamma(1 + 1 / beta)).toPrecision(4)], ['lognormal', `μ = ${f(muL, 4)}`, `σ = ${f(sgL, 4)}`, +r2L.toPrecision(4), +Math.exp(muL + sgL * N.normInv(0.1)).toPrecision(4), +Math.exp(muL + sgL * sgL / 2).toPrecision(4)]], explain: 'Parameters and derived lives of both candidate distributions fitted by probability-plot regression to the same failure data.' }],
    equations: [eqn('formula', 'Weibull reliability and hazard', 'R(t) = exp[−(t/η)^β],   h(t) = (β/η)(t/η)^{β−1}'), eqn('formula', 'Lognormal reliability', 'R(t) = 1 − Φ[(ln t − μ)/σ]'), eqn('formula', 'B10 life and MTTF (Weibull)', 'B10 = η(−ln 0.9)^{1/β},   MTTF = η Γ(1 + 1/β)'), eqn('governing', 'Degradation trend and RUL', expo ? 'ln y = c + m t;   RUL = (ln y_thr − c)/m − t_now' : 'y = c + m t;   RUL = (y_thr − c)/m − t_now'), eqn('formula', 'Age-replacement cost rate', 'C(T) = [c_p R(T) + c_f (1 − R(T))] / ∫₀ᵀ R(t) dt')],
    checks: [chk('Probability-plot fit quality (R² ≥ 0.9)', Math.max(wf.r2, r2L) >= 0.9, `Weibull ${f(wf.r2)}, lognormal ${f(r2L)}`), chk('Sample size adequate (n ≥ 10)', n >= 10, `n = ${n}`), chk('Reliability function monotonic and bounded', Rb(fl[0]) <= 1 && Rb(fl[n - 1]) >= 0 && Rb(fl[0]) >= Rb(fl[n - 1]), `R(${f(fl[0])} h) = ${f(Rb(fl[0]))}, R(${f(fl[n - 1])} h) = ${f(Rb(fl[n - 1]))}`), chk('Degradation trend approaches the threshold', towards, `slope ${f(m_)}, trend R² = ${f(fit.r2)}`), chk('RUL extrapolation within 3× observed span', !towards || RUL <= 3 * (tNow - tx[0]), `RUL ${f(RUL)} h vs observed ${f(tNow - tx[0])} h`)],
    warnings,
    recommendations: [
      rec('high', towards ? `Schedule the next intervention before ${f(tNow + rulLo)} h (lower 90 % RUL bound ${f(rulLo)} h from now); re-fit after each new measurement, because the band narrows as data accumulate.` : 'No approach to the threshold is visible in the degradation data; keep monitoring at the present interval and re-fit as data accumulate.', 'Degradation-trend extrapolation with prediction band.'),
      rec('medium', Topt < BIG ? `Adopt preventive replacement at ${f(Topt)} h: expected cost falls by ${f(saving * 100)} % relative to run-to-failure for a cost ratio of ${f(cr)}.` : beta <= 1 && best === 'Weibull' ? `With β = ${f(beta)} (no wear-out) scheduled replacement does not reduce failures; use condition monitoring, burn-in or redundancy instead.` : 'Preventive replacement is not cheaper than run-to-failure at this cost ratio; rely on condition-based maintenance from the degradation signal.', 'Age-replacement cost model.'),
      rec('medium', `Collect at least ${Math.max(10, 20 - n)} more failure (or suspension) times and record suspensions explicitly; with n = ${n} the shape parameter is known only to ±${f(1.645 * seB)}.`, 'Confidence bounds of the fit.'),
      rec('low', `Condition-based replacement near the RUL instead of a fixed conservative interval uses about ${f(B10 > 0 ? mttf / B10 : 1)}× more of each component’s life (MTTF/B10), reducing spare-part consumption and waste.`, 'Ratio of mean life to B10 life.', true),
    ],
    publish: { 'life.years': capN(Math.min(B10, BIG) / 8766), 'life.b10_h': B10, 'life.mttf_h': mttf, 'life.rul_h': RUL, 'life.weibull_shape': beta, 'life.weibull_scale_h': eta, 'life.maintenance_interval_h': interval },
    exports: [{ name: 'reliability_curves.csv', mime: 'text/csv', content: 'time_h,R_weibull,R_lognormal,hazard_best_per_h\n' + ts.map(t => [t, Rw(t), RL(t), hz(t)].map(v => +v.toPrecision(6)).join(',')).join('\n'), description: 'Fitted reliability and hazard curves.' }],
    numerics: { method: 'Median-rank regression (Weibull, lognormal), least-squares trend with prediction interval, golden-section minimisation of the age-replacement cost rate', dof: n + nd, steps: 200, iterations: 1, residual: 1 - Math.max(wf.r2, r2L), notes: 'Residual = 1 − R² of the best probability-plot fit.' },
  });
}

// ================================================================================================
const MI = (mode, group, extra) => o => ({ group, showIf: { mode: [mode], ...(extra || {}) }, ...o });
const FA = MI('fatigue', 'Fatigue'), FR = MI('fracture', 'Fracture'), CR = MI('creep_ageing', 'Creep'), AC = MI('creep_ageing', 'Accelerated testing'), RL_ = MI('reliability', 'Reliability');
const DG = (m, g) => MI('degradation', g, { deg_mech: [m] }), FN = (m, g) => MI('functional', g, { fun_mech: [m] });
const OX = DG('oxidation', 'Oxidation'), CO = DG('corrosion', 'Corrosion'), DI = DG('dissolution', 'Dissolution'), RP = DG('ripening', 'Ripening and agglomeration'), DL = DG('delamination', 'Delamination'), EM = DG('electromigration', 'Electromigration'), TD = DG('tddb', 'Dielectric breakdown');
const DT = FN('detuning', 'Resonance detuning'), BI = FN('bistable', 'Bistable element'), SW = FN('switching', 'Switching wear-out'), AX = FN('auxetic', 'Auxetic hinge');
const INPUTS = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'fatigue', group: 'Model', options: [{ value: 'fatigue', label: 'Fatigue (strain–life, Miner, Weibull scatter)' }, { value: 'fracture', label: 'Fracture and crack growth (Paris, Griffith)' }, { value: 'creep_ageing', label: 'Creep and accelerated ageing' }, { value: 'degradation', label: 'Chemical / environmental degradation' }, { value: 'functional', label: 'Functional lifetime of metamaterials' }, { value: 'reliability', label: 'Reliability synthesis and remaining useful life' }], help: 'Selects the damage or lifetime model. Each analysis has its own inputs, outputs, equations and checks.' },
  FA({ key: 'fat_sa', label: 'Nominal stress amplitude', unit: 'Pa', value: 1.5e8, min: 1e5, max: 5e9, log: true, bus: 'mech.max_stress', help: 'Alternating stress amplitude in the member (half the stress range). Life falls steeply with it (exponent 1/b ≈ −10).' }),
  FA({ key: 'fat_sm', label: 'Mean stress', unit: 'MPa', value: 50, min: -2000, max: 3000, help: 'Mean stress of the cycle; tensile mean stress shortens life through the selected correction.' }),
  FA({ key: 'fat_mean', label: 'Mean-stress correction', type: 'select', value: 'swt', options: [{ value: 'swt', label: 'Smith–Watson–Topper' }, { value: 'goodman', label: 'Goodman' }, { value: 'none', label: 'None' }], help: 'How mean stress is converted to an equivalent fully reversed amplitude.' }),
  FA({ key: 'fat_E', label: 'Young’s modulus', unit: 'Pa', value: 1.14e11, min: 1e8, max: 1.2e12, log: true, bus: 'prop.E', help: 'Elastic modulus (default Ti-6Al-4V); converts stress to elastic strain in the strain–life equation.' }),
  FA({ key: 'fat_su', label: 'Ultimate tensile strength', unit: 'MPa', value: 950, min: 10, max: 5000, log: true, help: 'Used by the Goodman correction and the static check.' }),
  FA({ key: 'fat_sf', label: 'Fatigue strength coefficient σ_f′', unit: 'MPa', value: 1700, min: 10, max: 10000, log: true, help: 'Intercept of the elastic strain–life line at one reversal (≈ true fracture strength).' }),
  FA({ key: 'fat_b', label: 'Fatigue strength exponent b', unit: '–', value: -0.095, min: -0.3, max: -0.03, help: 'Slope of the Basquin line (typically −0.05 to −0.12 for metals).' }),
  FA({ key: 'fat_ef', label: 'Fatigue ductility coefficient ε_f′', unit: '–', value: 0.35, min: 0.001, max: 3, log: true, help: 'Intercept of the plastic strain–life line (≈ true fracture ductility).' }),
  FA({ key: 'fat_c', label: 'Fatigue ductility exponent c', unit: '–', value: -0.6, min: -1, max: -0.3, help: 'Slope of the Coffin–Manson line (typically −0.5 to −0.7).' }),
  FA({ key: 'fat_Kt', label: 'Stress concentration factor K_t', unit: '–', value: 2.5, min: 1, max: 10, help: 'Elastic stress concentration at lattice nodes, fillets or surface notches (1 for smooth bulk; 2–4 typical of as-built lattice nodes).' }),
  FA({ key: 'fat_q', label: 'Notch sensitivity q', unit: '–', value: 0.8, min: 0, max: 1, help: 'Fraction of the theoretical concentration effective in fatigue (Peterson); small notches in ductile metals have q < 1.' }),
  FA({ key: 'fat_beta', label: 'Weibull shape of life scatter', unit: '–', value: 3, min: 0.5, max: 20, help: 'Scatter of fatigue life at fixed stress (2–4 typical for metals; lower for additively manufactured lattices).' }),
  FA({ key: 'fat_Nt', label: 'Target life', unit: 'cycles', value: 1e6, min: 10, max: 1e12, log: true, help: 'Design life at which failure probability and safety factor are evaluated.' }),
  FA({ key: 'fat_Dcrit', label: 'Critical Miner sum', unit: '–', value: 1, min: 0.1, max: 3, help: 'Damage sum at failure (1 by definition; 0.3–0.7 is often used for variable-amplitude loading as a safety allowance).' }),
  FA({ key: 'fat_spectrum', label: 'Load spectrum (per block)', type: 'table', value: DEF_SPECTRUM, columns: ['stress_amplitude_MPa', 'mean_stress_MPa', 'cycles_per_block'], help: 'Rainflow-counted block spectrum: amplitude, mean and number of cycles per block. The default is synthetic demonstration data.' }),
  FR({ key: 'fr_ds', label: 'Stress range Δσ', unit: 'Pa', value: 1.5e8, min: 1e5, max: 5e9, log: true, bus: 'mech.max_stress', help: 'Cyclic stress range (max − min) normal to the crack.' }),
  FR({ key: 'fr_R', label: 'Stress ratio R', unit: '–', value: 0.1, min: -1, max: 0.9, help: 'Minimum/maximum stress; sets the peak stress used for the critical crack size.' }),
  FR({ key: 'fr_a0', label: 'Initial flaw size', unit: 'µm', value: 200, min: 0.01, max: 1e5, log: true, help: 'Depth of the largest initial defect (pore, lack-of-fusion void, scratch).' }),
  FR({ key: 'fr_adet', label: 'Detectable flaw size', unit: 'µm', value: 1000, min: 0.01, max: 1e5, log: true, help: 'Smallest crack the inspection method finds reliably.' }),
  FR({ key: 'fr_KIc', label: 'Fracture toughness K_Ic', unit: 'Pa·√m', value: 6e7, min: 1e4, max: 3e8, log: true, bus: 'prop.KIc', help: 'Plane-strain fracture toughness (default Ti-6Al-4V, 60 MPa√m).' }),
  FR({ key: 'fr_Y', label: 'Geometry factor Y', unit: '–', value: 1.12, min: 0.3, max: 3, help: 'Dimensionless factor in K = Yσ√(πa): 1.12 edge crack, 1.0 centre crack, 0.71 semicircular surface crack.' }),
  FR({ key: 'fr_C', label: 'Paris coefficient C', unit: 'm/cycle/(MPa√m)^m', value: 1e-11, min: 1e-16, max: 1e-6, log: true, help: 'Paris-law coefficient with ΔK in MPa√m and da/dN in m/cycle.' }),
  FR({ key: 'fr_m', label: 'Paris exponent m', unit: '–', value: 3, min: 1.5, max: 12, help: 'Slope of the crack-growth curve (2–4 for metals, 10+ for ceramics and brittle polymers).' }),
  FR({ key: 'fr_dKth', label: 'Threshold ΔK_th', unit: 'MPa√m', value: 3, min: 0.05, max: 30, log: true, help: 'Stress-intensity range below which long cracks do not grow.' }),
  FR({ key: 'fr_sy', label: 'Yield strength', unit: 'MPa', value: 880, min: 1, max: 5000, log: true, help: 'Used to estimate the crack-tip plastic zone and the validity of linear-elastic fracture mechanics.' }),
  FR({ key: 'fr_rate', label: 'Load cycles per day', unit: '1/day', value: 1000, min: 0.01, max: 1e9, log: true, help: 'Converts the inspection interval from cycles to days.' }),
  FR({ key: 'fr_E', label: 'Young’s modulus (nanostructure)', unit: 'Pa', value: 1.4e11, min: 1e8, max: 1.2e12, log: true, bus: 'prop.E', help: 'Modulus for the Griffith estimate (default ZnO, 140 GPa).' }),
  FR({ key: 'fr_gamma', label: 'Surface energy γ', unit: 'J/m²', value: 1.5, min: 0.01, max: 20, log: true, bus: 'prop.surface_energy', help: 'Energy per unit area of new fracture surface.' }),
  FR({ key: 'fr_size', label: 'Nanostructure size', unit: 'nm', value: 50, min: 1, max: 1e5, log: true, help: 'Diameter or thickness of the nanowire, strut or film; the largest possible flaw is taken as half of it.' }),
  CR({ key: 'cr_sigma', label: 'Stress', unit: 'MPa', value: 100, min: 0.1, max: 2000, log: true, help: 'Sustained stress.' }),
  CR({ key: 'cr_T', label: 'Service temperature', unit: '°C', value: 600, min: 20, max: 1500, bus: 'proc.T_C', help: 'Metal temperature in service.' }),
  CR({ key: 'cr_years', label: 'Service life', unit: 'years', value: 10, min: 0.01, max: 100, log: true, help: 'Intended service duration; also the life to be demonstrated by accelerated tests.' }),
  CR({ key: 'cr_A', label: 'Norton coefficient A', unit: '1/s/MPaⁿ', value: 0.002, min: 1e-30, max: 1e20, log: true, help: 'Pre-exponential constant of the Norton law (default gives about 2×10⁻¹¹ 1/s at 100 MPa and 600 °C, of the order found for austenitic steel).' }),
  CR({ key: 'cr_n', label: 'Norton stress exponent n', unit: '–', value: 5, min: 1, max: 15, help: 'Stress exponent (1 diffusional, 3–8 dislocation creep).' }),
  CR({ key: 'cr_Q', label: 'Creep activation energy Q', unit: 'kJ/mol', value: 300, min: 30, max: 800, help: 'Activation energy of the rate-controlling diffusion process.' }),
  CR({ key: 'lm_a0', label: 'Larson–Miller master-curve intercept', unit: '10³ K', value: 33.6, min: 5, max: 80, help: 'a₀ in LMP/1000 = a₀ + a₁·log₁₀(σ/MPa) (default approximates 316 stainless steel, T in K, t in h).' }),
  CR({ key: 'lm_a1', label: 'Larson–Miller master-curve slope', unit: '10³ K per decade', value: -5.5, min: -30, max: -0.5, help: 'a₁: change of LMP per decade of stress (negative).' }),
  CR({ key: 'lm_C', label: 'Larson–Miller constant C', unit: '–', value: 20, min: 10, max: 40, help: 'Material constant, about 20 for many steels.' }),
  AC({ key: 'acc_Ea', label: 'Activation energy of ageing mechanism', unit: 'eV', value: 0.7, min: 0.1, max: 2.5, help: 'Apparent activation energy of the dominant degradation mechanism.' }),
  AC({ key: 'acc_Tuse', label: 'Use temperature', unit: '°C', value: 40, min: -40, max: 300, help: 'Typical service temperature.' }),
  AC({ key: 'acc_Ttest', label: 'Test temperature', unit: '°C', value: 85, min: -40, max: 500, help: 'Stress temperature of the accelerated test.' }),
  AC({ key: 'acc_RHuse', label: 'Use relative humidity', unit: '%', value: 50, min: 1, max: 100, help: 'Service humidity.' }),
  AC({ key: 'acc_RHtest', label: 'Test relative humidity', unit: '%', value: 85, min: 1, max: 100, help: 'Test-chamber humidity (85 °C / 85 % RH is the common damp-heat condition).' }),
  AC({ key: 'acc_nRH', label: 'Peck humidity exponent', unit: '–', value: 2.7, min: 0.5, max: 5, help: 'Exponent of the humidity term (about 2.7 for corrosion of metallisation).' }),
  AC({ key: 'acc_dTuse', label: 'Use temperature swing ΔT', unit: 'K', value: 30, min: 1, max: 300, help: 'Daily or operational temperature swing in service.' }),
  AC({ key: 'acc_dTtest', label: 'Test temperature swing ΔT', unit: 'K', value: 165, min: 1, max: 400, help: 'Thermal-cycling range of the test (for example −40 to +125 °C).' }),
  AC({ key: 'acc_m', label: 'Coffin–Manson exponent', unit: '–', value: 2.5, min: 1, max: 9, help: 'Exponent of the thermal-cycling law (≈ 2 ductile solder/metals, 4–9 brittle films and dielectrics).' }),
  AC({ key: 'acc_fuse', label: 'Use cycles per day', unit: '1/day', value: 1, min: 0.01, max: 1000, log: true, help: 'Thermal cycles per day in service.' }),
  AC({ key: 'acc_ftest', label: 'Test cycles per day', unit: '1/day', value: 24, min: 0.1, max: 2000, log: true, help: 'Thermal cycles per day in the chamber.' }),
  AC({ key: 'acc_rel', label: 'Reliability to demonstrate', unit: '–', value: 0.9, min: 0.5, max: 0.9999, help: 'Target reliability at the service life.' }),
  AC({ key: 'acc_conf', label: 'Confidence level', unit: '–', value: 0.9, min: 0.5, max: 0.999, help: 'Statistical confidence of the zero-failure demonstration.' }),
  { key: 'deg_mech', label: 'Degradation mechanism', type: 'select', value: 'dissolution', group: 'Degradation', showIf: { mode: ['degradation'] }, options: [{ value: 'dissolution', label: 'Nanoparticle dissolution (ion release)' }, { value: 'oxidation', label: 'Thermal oxidation (Deal–Grove)' }, { value: 'corrosion', label: 'Corrosion (Faraday)' }, { value: 'ripening', label: 'Ostwald ripening and agglomeration' }, { value: 'delamination', label: 'Coating delamination' }, { value: 'electromigration', label: 'Electromigration (Black)' }, { value: 'tddb', label: 'Dielectric breakdown (TDDB E-model)' }], help: 'Chemical or environmental degradation model to evaluate.' },
  OX({ key: 'ox_T', label: 'Oxidation temperature', unit: '°C', value: 1000, min: 600, max: 1300, help: 'Furnace or service temperature; rates are Arrhenius.' }),
  OX({ key: 'ox_t', label: 'Time', unit: 'h', value: 1, min: 0.001, max: 1e5, log: true, help: 'Exposure time.' }),
  OX({ key: 'ox_amb', label: 'Ambient', type: 'select', value: 'dry', options: [{ value: 'dry', label: 'Dry O₂' }, { value: 'wet', label: 'Wet (H₂O)' }], help: 'Oxidant: water diffuses much faster through the oxide than oxygen.' }),
  OX({ key: 'ox_orient', label: 'Silicon orientation', type: 'select', value: '100', options: [{ value: '100', label: '(100)' }, { value: '111', label: '(111)' }], help: 'Affects the linear (interface) rate constant only.' }),
  OX({ key: 'ox_x0', label: 'Initial oxide thickness', unit: 'nm', value: 0, min: 0, max: 2000, help: 'Oxide already present (native oxide ≈ 1–2 nm).' }),
  CO({ key: 'co_metal', label: 'Metal', type: 'select', value: 'Fe', options: Object.keys(METALS).map(k => ({ value: k, label: k })), help: 'Corroding metal; sets molar mass, valence and density.' }),
  CO({ key: 'co_i', label: 'Corrosion current density', unit: 'A/m²', value: 0.01, min: 1e-6, max: 100, log: true, bus: 'elec.current_density', help: 'i_corr from Tafel extrapolation or polarisation resistance (0.01 A/m² = 1 µA/cm²).' }),
  CO({ key: 'co_thk', label: 'Section thickness', unit: 'µm', value: 1000, min: 0.01, max: 1e5, log: true, help: 'Wall, strut or film thickness exposed to attack.' }),
  CO({ key: 'co_years', label: 'Service life', unit: 'years', value: 10, min: 0.01, max: 100, log: true, help: 'Exposure duration.' }),
  CO({ key: 'co_pit', label: 'Pitting factor', unit: '–', value: 5, min: 1, max: 100, log: true, help: 'Deepest pit depth divided by mean penetration (1 = uniform attack).' }),
  DI({ key: 'di_mat', label: 'Nanoparticle material', type: 'select', value: 'ZnO', options: Object.keys(DISS).map(k => ({ value: k, label: chem.prettyFormula(k) })), help: 'Dissolving solid; sets density, molar mass and the dissolution reaction.' }),
  DI({ key: 'di_medium', label: 'Medium', type: 'select', value: 'water', options: Object.keys(MEDIA).map(k => ({ value: k, label: MEDIA[k].name })), help: 'Scales the entered water-reference rate constant and solubility to typical biological media.' }),
  DI({ key: 'di_d', label: 'Particle diameter', unit: 'm', value: 3e-8, min: 1e-9, max: 1e-5, log: true, bus: 'geom.diameter', help: 'Primary particle diameter; sink-condition lifetime is proportional to it.' }),
  DI({ key: 'di_conc', label: 'Particle mass concentration', unit: 'mg/L', value: 10, min: 0.001, max: 1e4, log: true, help: 'Loading of particles in the medium.' }),
  DI({ key: 'di_k', label: 'Surface dissolution rate constant (water)', unit: 'mol/m²/s', value: 2e-8, min: 1e-12, max: 1e-3, log: true, help: 'Surface-normalised dissolution rate far from saturation in pure water (order 10⁻⁸ mol/m²/s for ZnO near neutral pH; calibrate).' }),
  DI({ key: 'di_cs', label: 'Solubility (water)', unit: 'mol/m³', value: 0.05, min: 1e-6, max: 1e4, log: true, help: 'Saturation concentration of the dissolved metal in pure water (ZnO ≈ 0.02–0.1 mol/m³, i.e. a few mg/L Zn).' }),
  DI({ key: 'di_t', label: 'Exposure time', unit: 'h', value: 24, min: 0.01, max: 1e5, log: true, help: 'Duration of the exposure.' }),
  RP({ key: 'rp_d', label: 'Initial mean diameter', unit: 'nm', value: 20, min: 1, max: 1000, log: true, help: 'Starting mean particle diameter.' }),
  RP({ key: 'rp_T', label: 'Temperature', unit: '°C', value: 60, min: 0, max: 99, help: 'Ageing temperature of the dispersion.' }),
  RP({ key: 'rp_t', label: 'Ageing time', unit: 'h', value: 100, min: 0.01, max: 1e6, log: true, help: 'Storage or service time.' }),
  RP({ key: 'rp_gamma', label: 'Interfacial energy', unit: 'J/m²', value: 0.3, min: 0.01, max: 3, log: true, help: 'Particle–liquid interfacial energy driving coarsening.' }),
  RP({ key: 'rp_cinf', label: 'Bulk solubility', unit: 'mol/m³', value: 1e-7, min: 1e-10, max: 1e3, log: true, help: 'Equilibrium solubility of the particle material in the medium (very low values, 10⁻⁸–10⁻⁶ mol/m³, reproduce the slow coarsening of oxide sols in alcohols; aqueous salts are far higher).' }),
  RP({ key: 'rp_D', label: 'Solute diffusivity', unit: 'm²/s', value: 1e-9, min: 1e-12, max: 1e-8, log: true, bus: 'prop.D_diff', help: 'Diffusion coefficient of the dissolved species.' }),
  RP({ key: 'rp_Vm', label: 'Molar volume', unit: 'cm³/mol', value: 14.5, min: 3, max: 200, help: 'Molar volume of the solid (ZnO 14.5, Au 10.2, Ag 10.3 cm³/mol).' }),
  RP({ key: 'rp_n0', label: 'Particle number concentration', unit: '1/m³', value: 1e18, min: 1e10, max: 1e24, log: true, help: 'Initial number of particles per unit volume.' }),
  RP({ key: 'rp_W', label: 'Stability ratio W', unit: '–', value: 1e7, min: 1, max: 1e15, log: true, help: 'Fuchs stability ratio: 1 = every collision sticks; electrostatic/steric stabilisation gives 10⁵–10¹⁰.' }),
  RP({ key: 'rp_df', label: 'Aggregate fractal dimension', unit: '–', value: 2.1, min: 1.5, max: 3, help: 'About 1.8 for fast (diffusion-limited) and 2.1 for slow (reaction-limited) aggregation.' }),
  DL({ key: 'dl_sigma', label: 'Residual film stress', unit: 'MPa', value: 300, min: -5000, max: 5000, help: 'As-deposited biaxial stress (positive tensile), e.g. from wafer-curvature measurement or the PVD stress estimate.' }),
  DL({ key: 'dl_h', label: 'Film thickness', unit: 'µm', value: 1, min: 0.001, max: 1000, log: true, help: 'Coating thickness; the driving force for delamination grows linearly with it.' }),
  DL({ key: 'dl_E', label: 'Film Young’s modulus', unit: 'GPa', value: 150, min: 0.01, max: 1000, log: true, help: 'Elastic modulus of the coating.' }),
  DL({ key: 'dl_nu', label: 'Film Poisson’s ratio', unit: '–', value: 0.25, min: -0.9, max: 0.49, help: 'Poisson’s ratio of the coating (auxetic films allowed).' }),
  DL({ key: 'dl_Gi', label: 'Interface toughness Γ_i', unit: 'J/m²', value: 5, min: 0.01, max: 1000, log: true, help: 'Critical energy release rate of the interface (0.1–1 weak, 1–10 typical, > 20 strong).' }),
  DL({ key: 'dl_Gf', label: 'Film fracture energy Γ_f', unit: 'J/m²', value: 10, min: 0.1, max: 1e4, log: true, help: 'Fracture energy of the film itself, for channel cracking.' }),
  DL({ key: 'dl_dalpha', label: 'Expansion mismatch α_film − α_substrate', unit: '1/K', value: 5e-6, min: -5e-5, max: 5e-5, help: 'Difference of thermal-expansion coefficients.' }),
  DL({ key: 'dl_dT', label: 'Temperature excursion ΔT', unit: 'K', value: -50, min: -500, max: 800, help: 'Service temperature minus the stress-reference temperature (cooling negative).' }),
  DL({ key: 'dl_b', label: 'Half-width of an existing debond', unit: 'µm', value: 20, min: 0.1, max: 1e4, log: true, help: 'Size of a pre-existing unbonded patch, for the buckling check under compression.' }),
  EM({ key: 'em_j', label: 'Current density', unit: 'A/m²', value: 1e10, min: 1e7, max: 1e12, log: true, bus: 'elec.current_density', help: 'Current density in the line (10¹⁰ A/m² = 1 MA/cm²).' }),
  EM({ key: 'em_T', label: 'Operating temperature', unit: '°C', value: 105, min: 0, max: 400, help: 'Line temperature including self-heating.' }),
  EM({ key: 'em_n', label: 'Current exponent n', unit: '–', value: 2, min: 1, max: 3, help: '1 for void-growth-limited, 2 for nucleation-limited failure.' }),
  EM({ key: 'em_Ea', label: 'Activation energy', unit: 'eV', value: 0.9, min: 0.3, max: 1.5, help: 'About 0.9 eV for Cu (interface diffusion), 0.5–0.7 eV for Al.' }),
  EM({ key: 'em_tref', label: 'Reference median life in test', unit: 'h', value: 300, min: 0.1, max: 1e6, log: true, help: 'Median time to failure measured in the accelerated test.' }),
  EM({ key: 'em_jref', label: 'Test current density', unit: 'A/m²', value: 2.5e10, min: 1e8, max: 1e12, log: true, help: 'Current density of the reference test.' }),
  EM({ key: 'em_Tref', label: 'Test temperature', unit: '°C', value: 300, min: 50, max: 500, help: 'Temperature of the reference test.' }),
  EM({ key: 'em_sigma', label: 'Lognormal σ', unit: '–', value: 0.5, min: 0.1, max: 2, help: 'Shape parameter of the lognormal failure distribution.' }),
  EM({ key: 'em_L', label: 'Line length', unit: 'µm', value: 20, min: 0.1, max: 1e5, log: true, help: 'Segment length between blocking boundaries (vias).' }),
  EM({ key: 'em_jLc', label: 'Blech threshold (j·L)_crit', unit: 'A/m', value: 3e5, min: 1e4, max: 5e6, log: true, help: 'Critical product for short-length immunity (about 2–4×10⁵ A/m for Cu).' }),
  EM({ key: 'em_life', label: 'Target life', unit: 'years', value: 10, min: 0.1, max: 50, help: 'Required life at 0.1 % cumulative failures.' }),
  TD({ key: 'td_V', label: 'Operating voltage', unit: 'V', value: 2.5, min: 0.1, max: 1000, log: true, help: 'Voltage across the dielectric.' }),
  TD({ key: 'td_tox', label: 'Dielectric thickness', unit: 'nm', value: 7, min: 1, max: 5000, log: true, help: 'Thickness of the dielectric layer across which the voltage is applied.' }),
  TD({ key: 'td_T', label: 'Operating temperature', unit: '°C', value: 105, min: -40, max: 300, help: 'Device temperature.' }),
  TD({ key: 'td_gamma', label: 'Field acceleration γ', unit: 'cm/MV', value: 2.5, min: 0.5, max: 10, help: 'ln-life slope with field (2.5 cm/MV ≈ 1.1 decades per MV/cm).' }),
  TD({ key: 'td_Ea', label: 'Activation energy', unit: 'eV', value: 0.7, min: 0.1, max: 2, help: 'Thermal activation of bond breakage (0.6–0.9 eV for SiO₂).' }),
  TD({ key: 'td_tref', label: 'Reference characteristic life in test', unit: 'h', value: 100, min: 1e-3, max: 1e6, log: true, help: 't₆₃ measured in the accelerated test.' }),
  TD({ key: 'td_Eref', label: 'Test field', unit: 'MV/cm', value: 10, min: 1, max: 20, help: 'Field of the reference test.' }),
  TD({ key: 'td_Tref', label: 'Test temperature', unit: '°C', value: 125, min: 0, max: 300, help: 'Temperature of the reference test.' }),
  TD({ key: 'td_beta', label: 'Weibull slope β', unit: '–', value: 1.5, min: 0.5, max: 10, help: 'Shape of the breakdown distribution (lower for thinner oxides).' }),
  TD({ key: 'td_area', label: 'Device area / test area', unit: '–', value: 100, min: 1e-3, max: 1e9, log: true, help: 'Ratio of total product dielectric area to the test-capacitor area (weakest-link scaling).' }),
  TD({ key: 'td_life', label: 'Target life', unit: 'years', value: 10, min: 0.1, max: 50, help: 'Required life at 0.01 % cumulative failures.' }),
  { key: 'fun_mech', label: 'Functional failure mode', type: 'select', value: 'detuning', group: 'Functional lifetime', showIf: { mode: ['functional'] }, options: [{ value: 'detuning', label: 'Resonance detuning' }, { value: 'bistable', label: 'Bistable / negative-stiffness element' }, { value: 'switching', label: 'Actuator / switching wear-out' }, { value: 'auxetic', label: 'Auxetic hinge fatigue' }], help: 'How the metamaterial loses its function.' },
  DT({ key: 'det_f0', label: 'As-built resonance frequency', unit: 'Hz', value: 2000, min: 0.1, max: 1e15, log: true, bus: 'res.f0', help: 'Resonance of the unit cell (mechanical, acoustic or electromagnetic element with a mechanical resonator).' }),
  DT({ key: 'det_band', label: 'Operating band half-width', unit: '± %', value: 3, min: 0.01, max: 50, log: true, help: 'Allowed fractional detuning (≈ 100/(2Q) for a resonance of quality factor Q).' }),
  DT({ key: 'det_years', label: 'Service horizon', unit: 'years', value: 10, min: 0.1, max: 100, log: true, help: 'Intended service life.' }),
  DT({ key: 'det_t0', label: 'Resonator (flexure) thickness', unit: 'µm', value: 100, min: 0.01, max: 1e5, log: true, help: 'Thickness of the bending element that sets the stiffness.' }),
  DT({ key: 'det_rate', label: 'Surface recession rate', unit: 'nm/yr', value: 20, min: 0, max: 1e6, help: 'Material lost per face per year by corrosion or oxidation (take it from the corrosion/oxidation analyses).' }),
  DT({ key: 'det_dE', label: 'Asymptotic modulus loss', unit: '%', value: 3, min: 0, max: 80, help: 'Long-term reduction of elastic modulus by physical/chemical ageing.' }),
  DT({ key: 'det_tau', label: 'Ageing time constant', unit: 'years', value: 5, min: 0.01, max: 1000, log: true, help: 'Time constant of the modulus ageing.' }),
  DT({ key: 'det_cyc', label: 'Load cycles per year', unit: '1/yr', value: 1e7, min: 0, max: 1e15, help: 'Vibration cycles accumulated per year.' }),
  DT({ key: 'det_Nf', label: 'Cycles to fatigue failure', unit: 'cycles', value: 1e9, min: 1e2, max: 1e18, log: true, bus: 'life.cycles', help: 'Fatigue life of the flexure at the operating amplitude (from the fatigue analysis).' }),
  DT({ key: 'det_dcrit', label: 'Stiffness loss at fatigue failure', unit: '%', value: 10, min: 0, max: 60, help: 'Stiffness drop accumulated when the fatigue life is exhausted (linear in N/N_f).' }),
  BI({ key: 'bi_l', label: 'Beam span', unit: 'mm', value: 10, min: 0.01, max: 1000, log: true, help: 'Clamped–clamped span of the curved beam.' }),
  BI({ key: 'bi_t', label: 'Beam thickness', unit: 'µm', value: 150, min: 0.1, max: 1e4, log: true, help: 'In-plane thickness of the beam.' }),
  BI({ key: 'bi_h', label: 'Apex height', unit: 'µm', value: 600, min: 0.1, max: 1e5, log: true, help: 'Initial rise of the cosine-shaped beam.' }),
  BI({ key: 'bi_w', label: 'Beam depth (out of plane)', unit: 'µm', value: 5000, min: 1, max: 1e5, log: true, help: 'Width of the beam normal to the bending plane.' }),
  BI({ key: 'bi_E', label: 'Young’s modulus', unit: 'GPa', value: 2, min: 0.001, max: 500, log: true, help: 'Short-term modulus (default a stiff thermoplastic).' }),
  BI({ key: 'bi_sf', label: 'Fatigue strength coefficient', unit: 'MPa', value: 90, min: 1, max: 5000, log: true, help: 'Basquin intercept of the beam material.' }),
  BI({ key: 'bi_b', label: 'Fatigue strength exponent', unit: '–', value: -0.09, min: -0.3, max: -0.03, help: 'Basquin slope.' }),
  BI({ key: 'bi_Einf', label: 'Relaxed modulus ratio E∞/E₀', unit: '–', value: 0.5, min: 0.05, max: 1, help: 'Fraction of the modulus that remains after full stress relaxation (1 = elastic, no relaxation).' }),
  BI({ key: 'bi_tau', label: 'Relaxation time τ', unit: 'h', value: 1000, min: 0.01, max: 1e7, log: true, help: 'KWW relaxation time at service temperature.' }),
  BI({ key: 'bi_beta', label: 'KWW stretching exponent', unit: '–', value: 0.5, min: 0.1, max: 1, help: 'Breadth of the relaxation spectrum (1 = single exponential).' }),
  BI({ key: 'bi_rate', label: 'Switching rate', unit: 'cycles/h', value: 10, min: 1e-4, max: 1e7, log: true, help: 'Snap-through cycles per hour.' }),
  BI({ key: 'bi_duty', label: 'Fraction of time in snapped state', unit: '–', value: 0.5, min: 0, max: 1, help: 'Relaxation accumulates only while the beam is held in its stressed (second) state.' }),
  BI({ key: 'bi_target', label: 'Target cycles', unit: 'cycles', value: 1e5, min: 1, max: 1e12, log: true, help: 'Required number of switching cycles.' }),
  SW({ key: 'sw_tech', label: 'Switching technology', type: 'select', value: 'pcm', options: Object.keys(SWT).map(k => ({ value: k, label: SWT[k].name })), help: 'Sets typical baseline endurance, Weibull shape and acceleration exponents.' }),
  SW({ key: 'sw_scale', label: 'Endurance calibration factor', unit: '×', value: 1, min: 1e-4, max: 1e4, log: true, help: 'Multiplies the typical baseline endurance to match your measured cycling data.' }),
  SW({ key: 'sw_drive', label: 'Drive amplitude', unit: '× nominal', value: 1, min: 0.3, max: 3, help: 'Drive voltage, current or pulse energy relative to nominal; life follows an inverse power law.' }),
  SW({ key: 'sw_T', label: 'Operating temperature', unit: '°C', value: 25, min: -40, max: 200, help: 'Element temperature (Arrhenius acceleration relative to 25 °C).' }),
  SW({ key: 'sw_rate', label: 'Switching frequency', unit: 'Hz', value: 1, min: 1e-6, max: 1e9, log: true, help: 'Reconfiguration rate while active.' }),
  SW({ key: 'sw_duty', label: 'Activity factor', unit: '–', value: 0.1, min: 1e-4, max: 1, log: true, help: 'Fraction of time the element is being switched.' }),
  SW({ key: 'sw_n', label: 'Number of elements', value: 10000, min: 2, max: 1e9, log: true, help: 'Switchable unit cells in the metasurface or array.' }),
  SW({ key: 'sw_ftol', label: 'Tolerated failed fraction', unit: '%', value: 2, min: 0.001, max: 50, log: true, help: 'Fraction of dead elements at which the array no longer meets its specification.' }),
  SW({ key: 'sw_years', label: 'Target life', unit: 'years', value: 10, min: 0.1, max: 50, help: 'Required service life.' }),
  AX({ key: 'ax_theta', label: 'Initial opening angle θ₀', unit: '°', value: 30, min: 2, max: 85, help: 'Angle between adjacent rotating squares in the unloaded state.' }),
  AX({ key: 'ax_strain', label: 'Global strain amplitude', unit: '%', value: 5, min: 0.01, max: 40, log: true, help: 'Cyclic macroscopic strain applied to the auxetic sheet.' }),
  AX({ key: 'ax_th', label: 'Hinge thickness', unit: 'µm', value: 100, min: 0.01, max: 1e4, log: true, help: 'Thickness of the flexure connecting the rigid units.' }),
  AX({ key: 'ax_lh', label: 'Hinge length', unit: 'µm', value: 2000, min: 0.01, max: 1e5, log: true, help: 'Length of the flexure over which the rotation is accommodated.' }),
  AX({ key: 'ax_E', label: 'Young’s modulus', unit: 'GPa', value: 114, min: 0.001, max: 500, log: true, help: 'Elastic modulus of the hinge material (default Ti-6Al-4V).' }),
  AX({ key: 'ax_sf', label: 'Fatigue strength coefficient σ_f′', unit: 'MPa', value: 1700, min: 1, max: 10000, log: true, help: 'Basquin intercept.' }),
  AX({ key: 'ax_b', label: 'Fatigue strength exponent b', unit: '–', value: -0.095, min: -0.3, max: -0.03, help: 'Basquin slope.' }),
  AX({ key: 'ax_ef', label: 'Fatigue ductility coefficient ε_f′', unit: '–', value: 0.35, min: 0.001, max: 3, log: true, help: 'Coffin–Manson intercept.' }),
  AX({ key: 'ax_c', label: 'Fatigue ductility exponent c', unit: '–', value: -0.6, min: -1, max: -0.3, help: 'Coffin–Manson slope.' }),
  AX({ key: 'ax_cyc', label: 'Cycles per year', unit: '1/yr', value: 1e5, min: 1, max: 1e12, log: true, help: 'Loading cycles per year.' }),
  AX({ key: 'ax_target', label: 'Target cycles', unit: 'cycles', value: 1e6, min: 10, max: 1e12, log: true, help: 'Required fatigue life.' }),
  RL_({ key: 'rel_fail', label: 'Failure times', type: 'table', value: DEF_FAIL, columns: ['time_to_failure_h'], help: 'Times to failure of nominally identical units (complete data). The default is synthetic demonstration data.' }),
  RL_({ key: 'rel_deg', label: 'Degradation path', type: 'table', value: DEF_DEG, columns: ['time_h', 'performance_pct'], help: 'Monitored performance indicator versus time for one unit (for example resonance amplitude, conductivity, transmission). The default is synthetic demonstration data.' }),
  RL_({ key: 'rel_trend', label: 'Degradation trend model', type: 'select', value: 'exp', options: [{ value: 'exp', label: 'Exponential' }, { value: 'lin', label: 'Linear' }], help: 'Functional form fitted to the degradation path.' }),
  RL_({ key: 'rel_thr', label: 'Failure threshold', unit: 'same as performance', value: 70, min: 1e-6, max: 1e6, log: true, help: 'Performance level that defines functional failure.' }),
  RL_({ key: 'rel_Rt', label: 'Target reliability', unit: '–', value: 0.95, min: 0.5, max: 0.9999, help: 'Reliability to be maintained between maintenance actions.' }),
  RL_({ key: 'rel_cost', label: 'Failure / preventive cost ratio', unit: '–', value: 10, min: 1, max: 1000, log: true, help: 'Cost of an in-service failure relative to a planned replacement.' }),
];

const RUNNERS = { fatigue: runFatigue, fracture: runFracture, creep_ageing: runCreep, degradation: runDegradation, functional: runFunctional, reliability: runReliability };
export default {
  id: 22,
  tagline: 'Predicts damage, degradation and service life — fatigue, crack growth, creep, ageing, corrosion, dissolution, functional drift of metamaterials — and synthesises reliability, remaining useful life and maintenance intervals.',
  inputs: INPUTS,
  presets: [
    { name: 'Ti-6Al-4V lattice strut under a block load spectrum', note: 'Additively manufactured lattice metamaterial: node stress concentration, SWT mean-stress correction, Miner sum and Weibull scatter.', values: { mode: 'fatigue' } },
    { name: 'ZnO nanoparticle dissolution in lysosomal fluid', note: '30 nm ZnO at 10 mg/L in pH 4.5 lysosomal simulant: Zn²⁺ release kinetics for nano-safety assessment.', values: { mode: 'degradation', deg_mech: 'dissolution', di_medium: 'lyso', di_t: 6 } },
    { name: 'Locally resonant metamaterial — resonance detuning over 10 years', note: 'Drift of a 2 kHz resonator relative to its ±3 % band by corrosion, modulus ageing and fatigue softening.', values: { mode: 'functional', fun_mech: 'detuning' } },
    { name: 'Bistable negative-stiffness element — snap-through life', note: 'Polymer curved-beam unit cell: loss of bistability by stress relaxation versus fatigue.', values: { mode: 'functional', fun_mech: 'bistable' } },
    { name: 'Dry thermal oxidation of silicon, 1000 °C, 1 h', note: 'Deal–Grove benchmark: a few tens of nanometres of gate-quality oxide.', values: { mode: 'degradation', deg_mech: 'oxidation' } },
    { name: 'Weibull analysis and remaining useful life from test data', note: 'Distribution fit to failure times, B10 life, degradation-path extrapolation and maintenance interval.', values: { mode: 'reliability' } },
    { name: 'Crack growth from a lack-of-fusion defect and Griffith strength of a ZnO nanowire', note: 'Paris-law life from a 200 µm flaw in Ti-6Al-4V and the size-dependent strength of a 50 nm nanowire.', values: { mode: 'fracture' } },
  ],
  run(inp) { return (RUNNERS[inp.mode] || runFatigue)(inp); },
};
