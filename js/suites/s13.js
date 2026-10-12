// Suite 13 — Magnetic and Spintronic Engineering.
// Stoner–Wohlfarth hysteresis, macrospin Landau–Lifshitz–Gilbert dynamics with spin-transfer torque, superparamagnetic
// relaxation and hyperthermia, ferromagnetic resonance with gyrotropic (negative) permeability, and spintronic figures of merit.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { PI, K_B, MU0, Q_E, HBAR } = N, GAMMA = 1.76085963e11; // rad s⁻¹ T⁻¹
const f3 = v => N.fmt(v, 3), f4 = v => N.fmt(v, 4);
const cx = z => `${f3(z.re)} ${z.im < 0 ? '−' : '+'} ${f3(Math.abs(z.im))}i`;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ser = (key, title, xlabel, ylabel, lines, explain, o = {}) => ({ key, title, xlabel, ylabel, logx: !!o.logx, logy: !!o.logy, lines, explain });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const pos = v => (v > 1e-300 ? v : 1e-300);
function finish(o) {
  o.warnings = o.warnings || []; o.tables = o.tables || []; let bad = 0;
  for (const k of o.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; bad++; }
  for (const s of o.series) for (const l of s.lines) {
    if (l.x.length > 400) { const st = (l.x.length - 1) / 399, ix = Array.from({ length: 400 }, (_, i) => Math.round(i * st)); l.x = ix.map(i => l.x[i]); l.y = ix.map(i => l.y[i]); }
    l.x = Array.from(l.x); l.y = Array.from(l.y, v => { if (!Number.isFinite(v)) { bad++; return 0; } return v; });
  }
  for (const k of Object.keys(o.publish)) if (typeof o.publish[k] === 'number' && !Number.isFinite(o.publish[k])) delete o.publish[k];
  if (bad) o.warnings.push('Some outputs were not computable for this input combination and are reported as 0 (non-finite values suppressed).');
  return o;
}

// Handbook-typical room-temperature values: Ms (A/m), |K| (J/m³), exchange stiffness A (J/m), density (kg/m³)
const MAG = {
  Fe3O4: { label: 'Magnetite Fe₃O₄', f: 'Fe3O4', Ms: 4.8e5, K: 1.35e4, A: 1.2e-11, rho: 5180, syn: [['FeCl2', 'FeCl3', 'NaOH'], ['Fe3O4', 'NaCl', 'H2O']] },
  gFe2O3: { label: 'Maghemite γ-Fe₂O₃', f: 'Fe2O3', Ms: 3.9e5, K: 4.7e3, A: 1.0e-11, rho: 4900, syn: [['FeCl3', 'NaOH'], ['Fe2O3', 'NaCl', 'H2O']] },
  Co: { label: 'Cobalt (hcp)', f: 'Co', Ms: 1.44e6, K: 4.5e5, A: 3.0e-11, rho: 8900 }, Fe: { label: 'Iron (bcc)', f: 'Fe', Ms: 1.71e6, K: 4.8e4, A: 2.1e-11, rho: 7874 },
  Ni: { label: 'Nickel (fcc)', f: 'Ni', Ms: 4.9e5, K: 5.7e3, A: 0.9e-11, rho: 8908 }, CoFe2O4: { label: 'Cobalt ferrite CoFe₂O₄', f: 'CoFe2O4', Ms: 4.25e5, K: 2.0e5, A: 1.0e-11, rho: 5290, syn: [['CoCl2', 'FeCl3', 'NaOH'], ['CoFe2O4', 'NaCl', 'H2O']] },
  FePt: { label: 'FePt (L1₀)', f: 'FePt', Ms: 1.14e6, K: 6.6e6, A: 1.0e-11, rho: 15200 }, YIG: { label: 'Yttrium iron garnet Y₃Fe₅O₁₂', f: 'Y3Fe5O12', Ms: 1.4e5, K: 610, A: 3.7e-12, rho: 5170 },
  Py: { label: 'Permalloy Ni₈₀Fe₂₀', f: 'Ni4Fe', Ms: 8.0e5, K: 500, A: 1.3e-11, rho: 8700 }, custom: { label: 'Custom (use the numeric inputs)', f: null },
};
function material(inp) { const m = MAG[inp.material] || MAG.Fe3O4; return m.f ? m : { ...m, Ms: inp.Ms, K: inp.K, A: inp.A_ex * 1e-12, rho: inp.rho }; }
function matEqs(inp, m, W) {
  const out = []; if (!m.f) return out;
  try { out.push(eqn('formula', 'Magnetic material', `${chem.prettyFormula(m.f)}, M = ${f4(chem.molarMass(m.f))} g/mol, M_s = ${f3(m.Ms)} A/m, K = ${f3(m.K)} J/m³, A = ${f3(m.A)} J/m`)); } catch { /* skip */ }
  if (m.syn) { try { const io = chem.ionicEquations(m.syn[0].map(f => ({ f, phase: 'aq' })), m.syn[1].map((f, i) => ({ f, phase: i === 0 ? 's' : f === 'H2O' ? 'l' : 'aq' }))); out.push(eqn('chemical', 'Co-precipitation synthesis (molecular)', io.molecular), eqn('ionic', 'Complete ionic equation', io.completeIonic), eqn('ionic', 'Net ionic equation', io.netIonic, 'spectators: ' + io.spectators.join(', '))); } catch (e) { W.push('Synthesis equation: ' + e.message); } }
  return out;
}
const sustain = (inp, m) => (['Co', 'CoFe2O4', 'FePt'].includes(inp.material) ? rec('medium', `${m.label} relies on cobalt or platinum, both critical raw materials with supply and toxicity concerns; check whether iron oxides or Fe–Ni alloys meet the specification, and plan recovery at end of life.`, 'Critical raw materials', true) : rec('low', 'Iron- and nickel-based magnetic materials are abundant; minimise the magnetic volume to the amount the function requires and avoid rare-earth or cobalt additions unless coercivity demands them.', 'Material efficiency and abundance', true));
const sphereV = d => PI / 6 * d ** 3;

// ---------------------------------------------------------------- Stoner–Wohlfarth
const NTH = 720;
function swLoop(psi, hmax, nH) { // follows the metastable minimum of e(θ) = sin²θ − 2h cos(θ − ψ)
  const dth = 2 * PI / NTH, s2 = new Float64Array(NTH), cp = new Float64Array(NTH); for (let i = 0; i < NTH; i++) { s2[i] = Math.sin(i * dth) ** 2; cp[i] = Math.cos(i * dth - psi); }
  const hs = [...N.linspace(hmax, -hmax, nH), ...N.linspace(-hmax, hmax, nH)], ms = []; let i = Math.round(psi / dth) % NTH;
  for (const h of hs) {
    const e = k => s2[k] - 2 * h * cp[k];
    for (let it = 0; it < NTH; it++) { const a = (i + 1) % NTH, b = (i + NTH - 1) % NTH, ea = e(a), eb = e(b), ec = e(i); if (ea < ec && ea <= eb) i = a; else if (eb < ec) i = b; else break; }
    const ea = e((i + 1) % NTH), eb = e((i + NTH - 1) % NTH), ec = e(i), den = eb - 2 * ec + ea, off = den > 1e-14 ? 0.5 * (eb - ea) / den : 0;
    ms.push(Math.cos((i + N.clamp(off, -0.5, 0.5)) * dth - psi));
  }
  return { h: hs, m: ms };
}
function loopStats(h, m, nH) { // descending branch = first nH points
  let hc = 0, mr = 0, area = 0;
  for (let i = 1; i < nH; i++) { if (m[i - 1] > 0 && m[i] <= 0) hc = Math.abs(h[i - 1] + (h[i] - h[i - 1]) * m[i - 1] / (m[i - 1] - m[i])); if (h[i - 1] > 0 && h[i] <= 0) mr = m[i - 1] + (m[i] - m[i - 1]) * h[i - 1] / (h[i - 1] - h[i]); }
  for (let i = 1; i < h.length; i++) area += 0.5 * (m[i] + m[i - 1]) * (h[i] - h[i - 1]);
  return { hc, mr, area: Math.abs(area) };
}
function sw(inp) {
  const m = material(inp), W = [], T = inp.T + 273.15, Hk = 2 * m.K / (MU0 * m.Ms), psi = N.clamp(inp.psi, 0, 90) * PI / 180, hmax = inp.hmax, nH = 201;
  const one = swLoop(psi, hmax, nH), st1 = loopStats(one.h, one.m, nH), nA = 24, ens = new Float64Array(2 * nH);
  for (let j = 0; j < nA; j++) { const l = swLoop(Math.acos((j + 0.5) / nA), hmax, nH); for (let i = 0; i < 2 * nH; i++) ens[i] += l.m[i] / nA; }
  const stE = loopStats(one.h, Array.from(ens), nH), hswAn = 1 / (Math.cos(psi) ** (2 / 3) + Math.sin(psi) ** (2 / 3)) ** 1.5;
  let jmp = 0, hsw = 0; for (let i = 1; i < nH; i++) { const d = one.m[i - 1] - one.m[i]; if (d > jmp) { jmp = d; hsw = Math.abs(one.h[i]); } }
  const V = sphereV(inp.d_nm * 1e-9), sig = m.K * V / (K_B * T), lnft = Math.log(Math.max(inp.tau_m, 1e-8) / 1e-9), red = sig > lnft ? 1 - (lnft / sig) ** 0.77 : 0, HcT = stE.hc * Hk * red;
  const loss = stE.area * 2 * m.K, Pm = loss * inp.f_ac * 1e3 / m.rho, dh = 2 * hmax / (nH - 1);
  if (hmax < 1) W.push('The field sweep does not reach the anisotropy field: particles with easy axes near the field direction do not switch and the loop is a minor loop.');
  if (red === 0) W.push(`At ${inp.d_nm} nm and ${inp.T} °C the particles are superparamagnetic on the ${f3(inp.tau_m)} s time scale (KV/kT = ${f3(sig)}): the zero-temperature Stoner–Wohlfarth loop collapses and coercivity vanishes.`);
  if (inp.d_nm * 1e-9 > 72 * Math.sqrt(m.A * m.K) / (MU0 * m.Ms * m.Ms)) W.push('The particle exceeds the single-domain size estimate; domain walls or vortex states reduce the coercivity well below the coherent-rotation value.');
  const hd = one.h.slice(0, nH), ha = one.h.slice(nH);
  return {
    summary: `${m.label}: anisotropy field H_K = ${f3(Hk / 1e3)} kA/m (µ₀H_K = ${f3(MU0 * Hk * 1e3)} mT). A single particle with the field at ${inp.psi}° to its easy axis has coercivity ${f3(st1.hc)} H_K and remanence ${f3(st1.mr)} M_s; a randomly oriented ensemble gives H_c = ${f3(stE.hc)} H_K = ${f3(stE.hc * Hk / 1e3)} kA/m and M_r/M_s = ${f3(stE.mr)}. With thermal activation at ${inp.T} °C a ${inp.d_nm} nm particle retains ${f3(100 * red)} % of that coercivity. The hysteresis loss is ${f3(loss)} J/m³ per cycle.`,
    kpis: [
      kpi('Hk', 'Anisotropy field', Hk, 'A/m', 'H_K = 2K/(µ₀M_s): the field needed to saturate along the hard axis and the upper bound of coercivity for coherent rotation.'),
      kpi('Hc1', 'Single-particle coercivity', st1.hc * Hk, 'A/m', `Field at which the magnetisation component along the field passes through zero for the chosen angle ${inp.psi}°, from tracking the metastable energy minimum.`),
      kpi('Hsw', 'Single-particle switching field', hsw * Hk, 'A/m', `Field of the irreversible jump; the Stoner–Wohlfarth astroid gives ${f3(hswAn)} H_K at this angle.`),
      kpi('Mr1', 'Single-particle remanence ratio', st1.mr, '–', 'M_r/M_s at zero field for the chosen angle, equal to cos ψ for coherent rotation.'),
      kpi('Hc', 'Ensemble coercivity (0 K)', stE.hc * Hk, 'A/m', 'Coercivity of non-interacting, randomly oriented uniaxial particles; the classical result is 0.479 H_K.'),
      kpi('Mr', 'Ensemble remanence ratio', stE.mr, '–', 'M_r/M_s of the random ensemble; exactly 0.5 for uniaxial anisotropy, so higher measured values indicate texture or cubic anisotropy.'),
      kpi('HcT', 'Ensemble coercivity with thermal activation', HcT, 'A/m', `Sharrock-type reduction H_c[1 − (ln(τ_m/τ₀)·kT/KV)^0.77] for ${inp.d_nm} nm particles measured over ${f3(inp.tau_m)} s.`, red > 0),
      kpi('loss', 'Hysteresis loss per cycle', loss, 'J/m³', 'Loop area µ₀∮M dH of the ensemble at 0 K; it is the energy dissipated per unit magnetic volume in each field cycle.'),
      kpi('P', 'Specific loss power at the stated frequency', Pm, 'W/kg', `Loop area × frequency (${inp.f_ac} kHz) per kilogram of magnetic material; an upper bound because thermal activation narrows real loops.`),
      kpi('BHmax', 'Energy-product bound', MU0 * (stE.mr * m.Ms) ** 2 / 4, 'J/m³', 'Upper limit µ₀M_r²/4 of (BH)max for the isotropic ensemble if the coercivity is sufficient; relevant to permanent-magnet nanocomposites.'),
      kpi('sigma', 'Thermal stability ratio KV/kT', sig, '–', 'Anisotropy barrier over thermal energy; values above about 25 are blocked for laboratory times and above 40–60 for data storage.', sig > 25),
    ],
    series: [
      ser('loop', 'Hysteresis loops', 'reduced field H/H_K (–)', 'M/M_s (–)', [{ name: `single particle, ψ = ${inp.psi}° (descending)`, x: hd, y: one.m.slice(0, nH) }, { name: 'single particle (ascending)', x: ha, y: one.m.slice(nH) }, { name: 'random ensemble (descending)', x: hd, y: Array.from(ens.slice(0, nH)) }, { name: 'random ensemble (ascending)', x: ha, y: Array.from(ens.slice(nH)) }], 'Magnetisation component along the field against reduced field. The single-particle loop shows reversible rotation ending in an abrupt switch; averaging over random easy-axis directions gives the rounded ensemble loop.'),
      ser('astroid', 'Switching field against field angle', 'field angle ψ (°)', 'H_sw/H_K (–)', [{ name: 'Stoner–Wohlfarth astroid', x: N.linspace(0, 90, 91), y: N.linspace(0, 90, 91).map(a => 1 / (Math.cos(a * PI / 180) ** (2 / 3) + Math.sin(a * PI / 180) ** (2 / 3)) ** 1.5) }, { name: 'this run (numerical)', x: [inp.psi], y: [hsw] }], 'Angular dependence of the switching field for coherent rotation: H_K along the easy and hard axes, with a minimum of H_K/2 at 45°. The point is the jump located in the numerical loop.'),
    ],
    equations: [...matEqs(inp, m, W),
      eqn('governing', 'Stoner–Wohlfarth energy', 'E/V = K sin²θ − µ₀ M_s H cos(θ − ψ);  ∂E/∂θ = 0 with ∂²E/∂θ² > 0 (metastable branch)'), eqn('formula', 'Astroid switching field', 'H_sw = H_K / (cos^(2/3)ψ + sin^(2/3)ψ)^(3/2),  H_K = 2K/(µ₀M_s)'),
      eqn('formula', 'Thermal reduction of coercivity', 'H_c(T) = H_c(0)·[1 − (kT·ln(τ_m/τ₀)/KV)^0.77]'), eqn('formula', 'Hysteresis loss', 'w = µ₀ ∮ M dH'),
    ],
    checks: [
      chk('Switching field matches the astroid', Math.abs(hsw - hswAn) < 2.5 * dh || hmax < hswAn, `numerical ${f4(hsw)} H_K, analytical ${f4(hswAn)} H_K (field step ${f3(dh)})`),
      chk('Ensemble remanence equals 0.5', hmax < 1.2 || Math.abs(stE.mr - 0.5) < 0.01, `M_r/M_s = ${f4(stE.mr)}`), chk('Ensemble coercivity equals 0.479 H_K', hmax < 1.2 || Math.abs(stE.hc - 0.479) < 0.015, `H_c/H_K = ${f4(stE.hc)}`),
      chk('Magnetisation bounded |M| ≤ M_s', one.m.every(v => Math.abs(v) <= 1 + 1e-9), 'coherent rotation conserves the moment magnitude'), chk('Loop is dissipative (area ≥ 0)', stE.area >= 0, `reduced area ${f4(stE.area)}`),
      chk('Thermally stable on the measurement time scale', red > 0, `KV/kT = ${f3(sig)} against ln(τ_m/τ₀) = ${f3(lnft)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', red === 0 ? `Increase the diameter above about ${f3(Math.cbrt(6 * lnft * K_B * T / (PI * m.K)) * 1e9)} nm (or choose a higher-anisotropy material) if remanence is required; at present the particles are superparamagnetic.` : `Thermal activation removes ${f3(100 * (1 - red))} % of the intrinsic coercivity; for recording or permanent-magnet use raise KV/kT from ${f3(sig)} towards 60.`, 'Sharrock thermal-activation law'),
      rec('medium', `Aligning easy axes (field-assisted assembly or texture) raises M_r/M_s from 0.5 towards 1 and the coercivity from 0.48 to 1.0 H_K (${f3(Hk / 1e3)} kA/m).`, 'Angular dependence of the Stoner–Wohlfarth loop'),
      rec('medium', 'Measured coercivities far below this estimate indicate incoherent reversal (curling, domain walls) or interparticle dipolar coupling; check particle size against the single-domain limit and dilute the sample.', 'Brown\'s paradox: coherent rotation is an upper bound'),
      rec('low', `For low-loss soft-magnetic use choose low K: the loss of ${f3(loss)} J/m³ per cycle scales directly with anisotropy.`, 'Loop area ∝ K'),
      sustain(inp, m),
    ],
    publish: { 'mag.Hc': HcT > 0 ? HcT : stE.hc * Hk, 'prop.Ms': m.Ms, 'mag.Hk': Hk, 'mag.Mr_ratio': stE.mr, 'mag.loss_J_m3': loss, 'material.formula': m.f || 'custom' },
    numerics: { method: 'Energy minimisation on a 0.5° angular grid following the metastable branch, parabolic refinement; 24-direction ensemble average', dof: NTH, steps: 2 * nH * (nA + 1), notes: 'Non-interacting single-domain particles, uniaxial anisotropy, quasi-static field, 0 K loop with separate thermal correction.' },
  };
}

// ---------------------------------------------------------------- macrospin LLG (+ spin-transfer torque)
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
function llg(inp) {
  const m = material(inp), W = [], T = inp.T + 273.15, al = inp.alpha, Bk = 2 * m.K / m.Ms, Ba = inp.B_pulse * 1e-3, thH = inp.th_H * PI / 180, Bdir = [Math.sin(thH), 0, -Math.cos(thH)];
  const V = sphereV(inp.d_nm * 1e-9), sig = m.K * V / (K_B * T), th0 = N.clamp(Math.sqrt(1 / (2 * Math.max(sig, 1e-9))), 0.01, 0.5), tf = inp.t_free * 1e-9, Jc = inp.J_stt, Pz = inp.P_spin;
  const BJ = HBAR * Pz * Jc / (2 * Q_E * m.Ms * tf), Jc0 = 2 * Q_E * al * m.Ms * tf * Bk / (HBAR * Pz), p = [0, 0, -1];
  const fL = GAMMA * (Ba + Bk) / (2 * PI), tp = inp.t_pulse * 1e-9; let tend = inp.t_end * 1e-9; const tcap = 150 / Math.max(fL, 1e6);
  if (tend > tcap) { W.push(`Simulated time shortened from ${f3(inp.t_end)} ns to ${f3(tcap * 1e9)} ns (150 precession periods) to bound the integration cost.`); tend = tcap; }
  const rhs = on => (t, y) => {
    const B = [on * Ba * Bdir[0], on * Ba * Bdir[1], on * Ba * Bdir[2] + Bk * y[2]], mxB = cross(y, B), mmB = cross(y, mxB), mxp = cross(y, p), mmp = cross(y, mxp), c = -GAMMA / (1 + al * al);
    return [0, 1, 2].map(i => c * (mxB[i] + al * mmB[i] + BJ * (mmp[i] - al * mxp[i])));
  };
  const y0 = [Math.sin(th0), 0, Math.cos(th0)], t1 = Math.min(tp, tend), opt = { rtol: 1e-7, atol: 1e-9, hmax: 0.05 / Math.max(fL, 1e6), maxSteps: 12000 };
  const s1 = N.rk45(rhs(1), y0, 0, t1, opt); let ts = s1.t, ys = s1.y, steps = s1.steps, done = s1.completed;
  if (tend > t1) { const s2 = N.rk45(rhs(0), ys[ys.length - 1], t1, tend, opt); ts = ts.concat(s2.t.slice(1)); ys = ys.concat(s2.y.slice(1)); steps += s2.steps; done = done && s2.completed; }
  if (!done) W.push('The adaptive integrator reached its step limit before the end time; the trajectory shown is truncated.');
  const mz = ys.map(y => y[2]), mx = ys.map(y => y[0]), nrm = N.maxOf(ys.map(y => Math.abs(Math.hypot(y[0], y[1], y[2]) - 1)));
  let tsw = 0, tset = 0; for (let i = 1; i < ts.length; i++) { if (!tsw && mz[i - 1] > 0 && mz[i] <= 0) tsw = ts[i]; if (!tset && mz[i] < -0.9) tset = ts[i]; }
  const switched = mz[mz.length - 1] < 0, hsw = 1 / (Math.cos(thH) ** (2 / 3) + Math.sin(thH) ** (2 / 3)) ** 1.5, Bsw = hsw * Bk, en = (y, on) => -m.K * y[2] * y[2] - m.Ms * on * Ba * (y[0] * Bdir[0] + y[2] * Bdir[2]);
  let mono = true; if (Jc === 0) { const scale = Math.max(m.K, m.Ms * Ba, 1e-30); for (let i = 1; i < ts.length; i++) { const on = ts[i] <= t1 ? 1 : 0, onp = ts[i - 1] < t1 ? 1 : 0; if (on === onp && en(ys[i], on) - en(ys[i - 1], on) > 1e-5 * scale) { mono = false; break; } } }
  const fFMR = GAMMA * Bk / (2 * PI), tau = 1 / (al * GAMMA * Math.max(Bk, 1e-9));
  if (!switched) W.push(Jc === 0 ? `The moment did not reverse: the pulse of ${inp.B_pulse} mT is below the static switching field of ${f3(Bsw * 1e3)} mT at this angle, or the pulse is too short for the precession to carry it over the barrier.` : 'The moment did not reverse within the simulated time: field and spin-torque drive are below threshold.');
  return {
    summary: `${m.label} macrospin (µ₀H_K = ${f3(Bk * 1e3)} mT, α = ${al}) under a ${inp.B_pulse} mT, ${inp.t_pulse} ns pulse at ${inp.th_H}° from the reversed easy axis${Jc ? ` with ${f3(Jc)} A/m² spin-polarised current` : ''}: the moment ${switched ? `reverses, crossing m_z = 0 at ${f3(tsw * 1e12)} ps` : 'does not reverse'}. Zero-field resonance is ${f3(fFMR / 1e9)} GHz and the damping time ${f3(tau * 1e12)} ps. The critical spin-torque current density is ${f3(Jc0)} A/m².`,
    kpis: [
      kpi('switched', 'Reversal outcome', switched ? 'switched' : 'not switched', '', 'Whether the easy-axis magnetisation component is negative at the end of the simulated time.', switched),
      kpi('t_sw', 'Switching time (m_z = 0)', tsw * 1e12, 'ps', 'Time at which the easy-axis component first crosses zero (0 if it never does); precessional switching takes about half a precession period of the applied field.'),
      kpi('t_set', 'Settling time (m_z < −0.9)', tset * 1e12, 'ps', 'Time to come within about 25° of the reversed state (0 if not reached); ringing decays with the Gilbert time constant.'),
      kpi('f_res', 'Zero-field resonance frequency', fFMR / 1e9, 'GHz', 'Small-angle precession frequency γµ₀H_K/2π about the anisotropy axis (Kittel formula for a uniaxial macrospin).'),
      kpi('tau', 'Gilbert damping time', tau * 1e12, 'ps', 'Relaxation time 1/(αγµ₀H_K) of the precession amplitude; write pulses must be followed by a few of these before the next operation.'),
      kpi('Bsw', 'Static switching field', Bsw * 1e3, 'mT', `Stoner–Wohlfarth astroid field at ${inp.th_H}° from the easy axis; pulses below it cannot switch quasi-statically.`, Ba >= Bsw),
      kpi('Jc0', 'Critical spin-torque current density', Jc0, 'A/m²', 'Zero-temperature threshold 2eαM_s t·µ₀H_K/(ħP) for the Slonczewski torque to overcome damping in a perpendicular macrospin.', Jc === 0 ? null : Jc > Jc0),
      kpi('BJ', 'Spin-torque effective field', BJ * 1e3, 'mT', 'Slonczewski torque amplitude ħPJ/(2eM_s t) expressed as a field; it competes with the damping field αµ₀H_K = ' + f3(al * Bk * 1e3) + ' mT.'),
      kpi('theta0', 'Initial thermal tilt', th0 * 180 / PI, '°', `Root-mean-square equilibrium cone angle √(kT/2KV) for a ${inp.d_nm} nm particle, used as the initial condition; reversal cannot start from perfect alignment.`),
      kpi('mz_end', 'Final easy-axis component', mz[mz.length - 1], '–', 'm_z at the end of the simulation: +1 is the initial state, −1 the reversed state.'),
    ],
    series: [
      ser('mz', 'Magnetisation dynamics', 'time (ps)', 'component (–)', [{ name: 'm_z (easy axis)', x: ts.map(t => t * 1e12), y: mz }, { name: 'm_x', x: ts.map(t => t * 1e12), y: mx }], 'Easy-axis and transverse components of the unit magnetisation. Reversal appears as m_z passing through zero; the oscillation of m_x is the precession, decaying at the Gilbert rate.'),
      ser('traj', 'Trajectory projected on the x–y plane', 'm_x (–)', 'm_y (–)', [{ name: 'trajectory', x: mx, y: ys.map(y => y[1]) }], 'Projection of the magnetisation path on the plane normal to the easy axis: an expanding and then contracting spiral shows the moment crossing the equator and relaxing to the reversed pole.'),
    ],
    equations: [...matEqs(inp, m, W),
      eqn('governing', 'Landau–Lifshitz–Gilbert–Slonczewski equation', 'dm/dt = −γ m×B_eff + α m×dm/dt − γ B_J m×(m×p),  B_eff = B_app + µ₀H_K m_z ẑ'),
      eqn('governing', 'Explicit (Landau–Lifshitz) form integrated', '(1 + α²) dm/dt = −γ[m×B_eff + α m×(m×B_eff) + B_J(m×(m×p) − α m×p)]'),
      eqn('formula', 'Spin-torque field and threshold', 'B_J = ħPJ/(2eM_s t),  J_c0 = 2eαM_s t µ₀H_K/(ħP)'), eqn('formula', 'Resonance and damping', 'f = γµ₀H_K/2π,  τ = 1/(αγµ₀H_K)'),
    ],
    checks: [
      chk('Magnetisation norm conserved', nrm < 1e-3, `max ||m| − 1| = ${f3(nrm)}`), chk('Energy non-increasing in each constant-field interval (Gilbert dissipation)', mono, Jc === 0 ? 'no energy gain between accepted steps' : 'skipped: spin-transfer torque is non-conservative and can pump energy'),
      chk('Integration completed', done, `${steps} adaptive steps`), chk('Reversal achieved', switched, switched ? `m_z = 0 at ${f3(tsw * 1e12)} ps` : 'final m_z > 0'),
      chk('Barrier is thermally stable (KV/kT > 40)', sig > 40, `KV/kT = ${f3(sig)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', switched ? `Reversal completes in ${f3((tset || tsw) * 1e12)} ps; trimming the pulse to about ${f3(Math.max(tsw, 1e-13) * 1.2e12)} ps (just beyond the m_z = 0 crossing) avoids back-switching by ringing and saves drive energy.` : `Raise the pulse field above ${f3(Bsw * 1e3)} mT, tilt it towards 45° from the easy axis (astroid minimum H_K/2), or add spin-torque current above ${f3(Jc0)} A/m².`, 'Astroid and precessional switching condition'),
      rec('medium', `Damping α = ${al}: higher damping shortens ringing (τ = ${f3(tau * 1e12)} ps) but raises the spin-torque threshold in proportion; choose it for the dominant write mechanism.`, 'Opposite roles of α in field and spin-torque switching'),
      rec('medium', 'Use micromagnetic simulation for elements larger than a few exchange lengths: non-uniform reversal lowers the real switching field and current below the macrospin values.', 'Macrospin approximation'),
      rec('low', `Spin-torque writing at J_c0 = ${f3(Jc0)} A/m² dissipates Joule heat in the junction; reducing M_s·t and α lowers the write energy per bit.`, 'Write-energy scaling of spin-transfer torque', true),
    ],
    publish: { 'mag.switching_time': tsw, 'res.f0': fFMR, 'elec.current_density': Jc0, 'prop.Ms': m.Ms, 'mag.Hk': Bk / MU0 },
    numerics: { method: 'Adaptive Dormand–Prince 5(4) integration of the explicit LLG–Slonczewski equation', dof: 3, steps, notes: 'Single macrospin, uniaxial anisotropy, no thermal field during the trajectory (thermal tilt only as initial condition).' },
  };
}

// ---------------------------------------------------------------- superparamagnetism and hyperthermia
const lang = x => (Math.abs(x) < 1e-4 ? x / 3 : 1 / Math.tanh(x) - 1 / x);
function spm(inp) {
  const m = material(inp), W = [], T = inp.T + 273.15, kT = K_B * T, d0 = inp.d_nm * 1e-9, sg = inp.sigma_ln, tau0 = 1e-9, eta = inp.visc * 1e-3, coat = inp.coat * 1e-9, H0 = inp.H_ac * 1e3, f = inp.f_ac * 1e3, w = 2 * PI * f, phi = inp.phi_v / 100;
  const tN = d => tau0 * Math.exp(Math.min(m.K * sphereV(d) / kT, 600)), tB = d => 3 * eta * sphereV(d + 2 * coat) / kT, tE = d => { const a = tN(d), b = tB(d); return a * b / (a + b); };
  const lnm = Math.log(Math.max(inp.tau_m, 1e-8) / tau0), TB = m.K * sphereV(d0) / (K_B * lnm), dSP = Math.cbrt(6 * lnm * kT / (PI * m.K)), Dsd = 72 * Math.sqrt(m.A * m.K) / (MU0 * m.Ms * m.Ms), dw = PI * Math.sqrt(m.A / m.K), lex = Math.sqrt(2 * m.A / (MU0 * m.Ms * m.Ms));
  // log-normal size distribution (number-weighted), 41-node quadrature in ln d
  const zs = N.linspace(-3.5, 3.5, 41), ds = zs.map(z => d0 * Math.exp(sg * z)), wn = zs.map(z => Math.exp(-z * z / 2)), vt = N.sum(wn.map((v, i) => v * sphereV(ds[i]))), wv = wn.map((v, i) => v * sphereV(ds[i]) / vt);
  const Mof = (H, list = ds, wts = wv) => m.Ms * N.sum(list.map((d, i) => wts[i] * lang(MU0 * m.Ms * sphereV(d) * H / kT)));
  const chi2 = (d, ph, ff = f) => { const Vd = sphereV(d), xi = MU0 * m.Ms * H0 * Vd / kT, chi0 = (MU0 * ph * m.Ms * m.Ms * Vd / (3 * kT)) * (xi < 1e-4 ? 1 : 3 / xi * lang(xi)), wt = 2 * PI * ff * tE(d); return chi0 * wt / (1 + wt * wt); };
  const Pvol = (list = ds, wts = wv, ff = f) => PI * MU0 * ff * H0 * H0 * N.sum(list.map((d, i) => chi2(d, phi * wts[i], ff)));
  const P = Pvol(), SAR = P / (m.rho * Math.max(phi, 1e-12)), ILP = SAR / (H0 * H0 * f) * 1e9, dTdt = P / 4.18e6, xi0 = MU0 * m.Ms * H0 * sphereV(d0) / kT, chi_i = MU0 * phi * m.Ms * m.Ms * N.sum(ds.map((d, i) => wv[i] * sphereV(d))) / (3 * kT), Hf = H0 * f;
  const dd = N.linspace(3e-9, Math.max(40e-9, 2.5 * d0), 160), sarD = dd.map(d => Pvol([d], [1]) / (m.rho * Math.max(phi, 1e-12))), dopt = dd[N.argmax(sarD)], Hs = N.linspace(-inp.H_dc * 1e3, inp.H_dc * 1e3, 161);
  const sig0 = m.K * sphereV(d0) / kT, blocked = sig0 > lnm, sat = Mof(inp.H_dc * 1e3) / m.Ms;
  if (blocked) W.push(`The mean particle is blocked at ${inp.T} °C on a ${f3(inp.tau_m)} s time scale (T_B = ${f3(TB)} K): expect hysteresis, and the equilibrium Langevin curve is then an anhysteretic idealisation.`);
  if (xi0 > 1) W.push(`Field parameter ξ = ${f3(xi0)} exceeds 1: linear-response theory overestimates heating at this amplitude; dynamic-hysteresis models are needed.`);
  if (Hf > 5e9) W.push(`H·f = ${f3(Hf)} A/m/s exceeds the commonly used clinical tolerance of about 5×10⁹ A/m/s (eddy-current heating of tissue).`);
  if (d0 > Dsd) W.push('Mean diameter exceeds the single-domain estimate; the macrospin relaxation model does not apply to multi-domain particles.');
  if (phi > 0.05) W.push('Above about 5 vol % dipolar interactions between particles alter relaxation times and susceptibility; results are for the non-interacting limit.');
  return {
    summary: `${m.label} particles of ${inp.d_nm} nm (σ = ${sg}) at ${inp.T} °C: KV/kT = ${f3(sig0)}, Néel time ${f3(tN(d0))} s, Brown time ${f3(tB(d0))} s; blocking temperature ${f3(TB)} K, so the mean particle is ${blocked ? 'blocked' : 'superparamagnetic'} (critical diameter ${f3(dSP * 1e9)} nm). Single-domain limit ≈ ${f3(Dsd * 1e9)} nm. In a ${inp.H_ac} kA/m, ${inp.f_ac} kHz field the linear-response SAR is ${f3(SAR)} W/kg (ILP ${f3(ILP)} nH·m²/kg), heating a ${inp.phi_v} vol % fluid at ${f3(dTdt)} K/s.`,
    kpis: [
      kpi('d_sp', 'Superparamagnetic critical diameter', dSP * 1e9, 'nm', `Diameter below which the Néel time is shorter than the measurement time ${f3(inp.tau_m)} s at this temperature, from KV = kT·ln(τ_m/τ₀) (≈ 25 kT for 100 s).`, d0 < dSP),
      kpi('TB', 'Blocking temperature', TB, 'K', 'Temperature below which the mean particle appears blocked, KV/(k·ln(τ_m/τ₀)); it is the ZFC peak position for a narrow distribution.'),
      kpi('sigma', 'Barrier ratio KV/kT', sig0, '–', 'Anisotropy energy barrier over thermal energy for the median particle; relaxation time grows exponentially with it.'),
      kpi('tauN', 'Néel relaxation time', tN(d0), 's', 'Time for the moment to flip over the anisotropy barrier inside the particle, τ₀·exp(KV/kT) with τ₀ = 10⁻⁹ s.'),
      kpi('tauB', 'Brown relaxation time', tB(d0), 's', `Rotational diffusion time of the whole particle with a ${inp.coat} nm shell in a fluid of ${inp.visc} mPa·s, 3ηV_h/kT.`),
      kpi('tauE', 'Effective relaxation time', tE(d0), 's', 'Parallel combination τ_Nτ_B/(τ_N + τ_B); heating peaks when 2πf·τ = 1, i.e. at ' + f3(1 / (2 * PI * tE(d0))) + ' Hz for the median particle.'),
      kpi('Dsd', 'Single-domain critical diameter', Dsd * 1e9, 'nm', 'Estimate 72√(AK)/(µ₀M_s²) above which a domain wall lowers the energy of a sphere; valid for strong anisotropy, an order-of-magnitude guide for soft materials.', d0 < Dsd),
      kpi('dw', 'Domain-wall width', dw * 1e9, 'nm', 'Bloch wall width π√(A/K); particles much smaller than this cannot support a wall.'),
      kpi('lex', 'Exchange length', lex * 1e9, 'nm', 'Length √(2A/µ₀M_s²) over which exchange keeps the magnetisation uniform against magnetostatic energy; it sets the micromagnetic cell size.'),
      kpi('chi', 'Initial susceptibility of the fluid', chi_i, '–', 'Langevin initial susceptibility µ₀φM_s²⟨V⟩/(3kT) of the suspension (SI, dimensionless), volume-averaged over the size distribution.'),
      kpi('sat', 'Magnetisation at the DC field', sat, '–', `M/M_s of the particle ensemble at ${inp.H_dc} kA/m from the size-averaged Langevin function.`),
      kpi('SAR', 'Specific absorption rate', SAR, 'W/kg', 'Heating power per kilogram of magnetic material from Rosensweig linear-response theory, averaged over the log-normal size distribution.'),
      kpi('ILP', 'Intrinsic loss power', ILP, 'nH·m²/kg', 'SAR normalised by H²f, which allows comparison between laboratories; commercial ferrofluids give roughly 0.2–4.'),
      kpi('P', 'Volumetric heating power', P, 'W/m³', `Power dissipated per cubic metre of fluid at ${inp.phi_v} vol % loading, P = πµ₀χ″fH₀².`),
      kpi('dTdt', 'Initial heating rate', dTdt, 'K/s', 'Adiabatic temperature rise rate of the aqueous suspension (4.18 MJ/m³/K).'),
      kpi('d_opt', 'Optimal monodisperse diameter', dopt * 1e9, 'nm', `Core diameter that maximises SAR at ${inp.f_ac} kHz for this material and fluid; the optimum is sharp, so polydispersity costs heating power.`),
      kpi('Hf', 'Field–frequency product', Hf, 'A/m/s', 'H₀·f, compared with the physiological tolerance of roughly 5×10⁹ A/m/s (4.85×10⁸ in the stricter Atkinson–Brezovich criterion).', Hf <= 5e9),
    ],
    series: [
      ser('MH', 'Equilibrium magnetisation curve', 'field (kA/m)', 'M/M_s (–)', [{ name: 'log-normal ensemble', x: Hs.map(h => h / 1e3), y: Hs.map(h => Mof(h) / m.Ms) }, { name: 'monodisperse (median size)', x: Hs.map(h => h / 1e3), y: Hs.map(h => lang(MU0 * m.Ms * sphereV(d0) * h / kT)) }], 'Langevin magnetisation of the superparamagnetic ensemble. A size distribution steepens the low-field slope (large particles) and slows the approach to saturation (small particles); fitting both yields the magnetic size and its spread.'),
      ser('tau', 'Relaxation times against core diameter', 'core diameter (nm)', 'time (s)', [{ name: 'Néel', x: dd.map(d => d * 1e9), y: dd.map(d => pos(Math.min(tN(d), 1e30))) }, { name: 'Brown', x: dd.map(d => d * 1e9), y: dd.map(d => pos(tB(d))) }, { name: 'effective', x: dd.map(d => d * 1e9), y: dd.map(d => pos(tE(d))) }, { name: '1/(2πf)', x: [dd[0] * 1e9, dd[159] * 1e9], y: [1 / w, 1 / w] }], 'Néel time rises exponentially with volume while Brown time rises only linearly; the faster process controls relaxation. Maximum heating occurs where the effective time crosses 1/(2πf).', { logy: true }),
      ser('sar', 'SAR against core diameter (monodisperse)', 'core diameter (nm)', 'SAR (W/kg)', [{ name: 'linear-response SAR', x: dd.map(d => d * 1e9), y: sarD }, { name: 'this distribution', x: [inp.d_nm], y: [SAR] }], 'Predicted heating of monodisperse particles against core size at the chosen field. The narrow maximum corresponds to 2πf·τ_N ≈ 1; the marked point is the value for the actual log-normal distribution.'),
    ],
    equations: [...matEqs(inp, m, W),
      eqn('constitutive', 'Néel–Arrhenius and Brown relaxation', 'τ_N = τ₀ exp(KV/kT),  τ_B = 3ηV_h/kT,  1/τ = 1/τ_N + 1/τ_B'), eqn('formula', 'Blocking condition', 'KV = kT ln(τ_m/τ₀)  ⇒  T_B = KV/(k ln(τ_m/τ₀))'),
      eqn('constitutive', 'Langevin magnetisation', 'M = φ M_s ⟨L(ξ)⟩_V,  L(ξ) = coth ξ − 1/ξ,  ξ = µ₀M_sVH/kT'), eqn('governing', 'Linear-response heating (Rosensweig)', 'P = πµ₀χ″fH₀²,  χ″ = χ₀·ωτ/(1 + (ωτ)²),  χ₀ = (µ₀φM_s²V/3kT)·(3/ξ)(coth ξ − 1/ξ)'),
      eqn('formula', 'Characteristic lengths', 'D_sd ≈ 72√(AK)/(µ₀M_s²),  δ_w = π√(A/K),  l_ex = √(2A/µ₀M_s²)'),
    ],
    checks: [
      chk('Linear-response validity (ξ < 1)', xi0 < 1, `ξ = ${f3(xi0)} at the field amplitude`), chk('Field–frequency product within tolerance', Hf <= 5e9, `H₀f = ${f3(Hf)} A/m/s`),
      chk('Single-domain particle', d0 < Dsd, `d = ${inp.d_nm} nm against D_sd ≈ ${f3(Dsd * 1e9)} nm`), chk('Size-distribution weights normalised', Math.abs(N.sum(wv) - 1) < 1e-12, 'volume weights sum to 1'),
      chk('Langevin curve bounded and odd', Math.abs(sat) <= 1 && Math.abs(Mof(0)) < 1e-9 * m.Ms, `M(0) = 0, |M/M_s| ≤ 1 (${f3(sat)} at the DC field)`), chk('Heating power non-negative (passivity)', P >= 0, `χ″ ≥ 0 gives P = ${f3(P)} W/m³`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `For hyperthermia at ${inp.f_ac} kHz the optimum core size is about ${f3(dopt * 1e9)} nm; the present ${inp.d_nm} nm median with σ = ${sg} delivers ${f3(100 * SAR / Math.max(N.maxOf(sarD), 1e-30))} % of the monodisperse maximum. Narrowing the distribution matters more than raising the field.`, 'Sharp SAR maximum at 2πf·τ = 1'),
      rec('medium', blocked ? 'The particles are blocked: for colloidal stability and zero remanence (biomedical, ferrofluid use) reduce the core size below the superparamagnetic limit or choose a lower-anisotropy phase.' : `The particles are superparamagnetic (no remanence), which suppresses magnetic agglomeration; keep cores below ${f3(dSP * 1e9)} nm to stay in this regime.`, 'Blocking criterion'),
      rec('medium', `Verify core size magnetically: fit the measured M(H) with the Langevin ensemble and compare with TEM; a magnetically dead surface layer of about 1 nm is common in ferrites and lowers the effective M_s.`, 'Langevin granulometry'),
      rec('low', `Stay below H₀f ≈ 5×10⁹ A/m/s (now ${f3(Hf)}); raising frequency rather than amplitude is more effective only while 2πf·τ < 1.`, 'Tissue eddy-current limit and Debye response'),
      sustain(inp, m),
    ],
    publish: { 'mag.SAR': SAR, 'mag.TB': TB, 'prop.Ms': m.Ms, 'rxn.heat': P, 'geom.diameter': d0, 'mag.tau_eff': tE(d0), 'material.formula': m.f || 'custom' },
    numerics: { method: 'Closed-form Néel–Brown relaxation and Langevin statistics; 41-node quadrature over the log-normal size distribution; Rosensweig linear-response loss', dof: zs.length, steps: dd.length, notes: 'Non-interacting uniaxial particles, attempt time 1 ns, no temperature dependence of M_s or K.' },
  };
}

// ---------------------------------------------------------------- FMR and gyrotropic permeability
const SHAPES = { sphere: ['Sphere', [1 / 3, 1 / 3, 1 / 3]], film_ip: ['Thin film, magnetised in plane', [1, 0, 0]], film_oop: ['Thin film, magnetised out of plane', [0, 0, 1]], rod: ['Long rod / nanowire, magnetised along axis', [0.5, 0.5, 0]] };
function fmr(inp) {
  const m = material(inp), W = [], al = inp.alpha, H = inp.H_bias * 1e3, Hk = 2 * m.K / (MU0 * m.Ms), [, Nd] = SHAPES[inp.shape] || SHAPES.sphere, ff = inp.f_fill, C = N.C;
  const kit = h => { const a = h + (Nd[0] - Nd[2]) * m.Ms + Hk, b = h + (Nd[1] - Nd[2]) * m.Ms + Hk; return a > 0 && b > 0 ? GAMMA * MU0 * Math.sqrt(a * b) / (2 * PI) : 0; };
  let Hi = H - Nd[2] * m.Ms + Hk; const fK = kit(H);
  if (Hi <= 0) { W.push(`The bias field does not saturate the sample against its demagnetising field (internal field ${f3(Hi / 1e3)} kA/m ≤ 0): the uniform-precession formulae do not apply; a nominal internal field of 1 kA/m is used.`); Hi = 1e3; }
  const w0 = GAMMA * MU0 * Hi, wm = GAMMA * MU0 * m.Ms, wV = Math.sqrt(w0 * (w0 + wm)), wAR = w0 + wm;
  const tens = w => { const wc = { re: w0, im: -al * w }, den = C.sub(C.mul(wc, wc), C.of(w * w)), mu = C.add(C.of(1), C.div(C.scale(wc, wm), den)), ka = C.div(C.of(wm * w), den), muV = C.div(C.sub(C.mul(mu, mu), C.mul(ka, ka)), mu), muP = C.add(C.of(1), C.div(C.of(wm), { re: w0 - w, im: -al * w })); return { mu, ka, muV, muP, muC: { re: 1 + ff * (muV.re - 1), im: ff * muV.im } }; };
  const ws = N.linspace(0.3 * w0, Math.max(1.5 * wAR, 2 * w0), 400), Tn = ws.map(tens), fl = ws.map(w => w / (2 * PI * 1e9));
  const bnd = (arr) => { let lo = 0, hi = 0; arr.forEach((v, i) => { if (v < 0) { if (!lo) lo = ws[i]; hi = ws[i]; } }); return lo ? [lo, hi] : null; }, bV = bnd(Tn.map(t => t.muV.re)), bC = bnd(Tn.map(t => t.muC.re)), bP = bnd(Tn.map(t => t.muP.re));
  const wpr = inp.f_rel * wV, tp = tens(wpr), dH = 2 * al * wV / (GAMMA * MU0), imin = N.argmin(Tn.map(t => t.muV.re)), Hs = N.linspace(0, Math.max(2 * H, 1e4), 120);
  if (!bV) W.push('Damping is too large for the effective permeability to become negative above resonance.');
  if (ff < 1 && !bC) W.push(`At a magnetic fill fraction of ${ff} the composite permeability stays positive: dilution weakens the resonance below the level needed for µ < 0.`);
  return {
    summary: `${m.label}, ${SHAPES[inp.shape][0].toLowerCase()}, biased at ${inp.H_bias} kA/m: Kittel resonance at ${f3(fK / 1e9)} GHz (linewidth ${f3(dH / 1e3)} kA/m for α = ${al}). The intrinsic transverse (Voigt) permeability is negative from ${f3(wV / (2e9 * PI))} to ${f3(wAR / (2e9 * PI))} GHz in the lossless limit${bV ? ` (computed with loss: ${f3(bV[0] / (2e9 * PI))}–${f3(bV[1] / (2e9 * PI))} GHz, minimum ${f3(Tn[imin].muV.re)})` : ''}. At ${f3(wpr / (2e9 * PI))} GHz µ_eff = ${cx(tp.muV)}.`,
    kpis: [
      kpi('fK', 'Kittel resonance frequency', fK / 1e9, 'GHz', `Uniform-precession frequency of the ${SHAPES[inp.shape][0].toLowerCase()} including shape (demagnetising) and uniaxial anisotropy fields; the measured FMR absorption peak.`, fK > 0),
      kpi('f0', 'Intrinsic precession frequency', w0 / (2e9 * PI), 'GHz', 'γµ₀H_i/2π for the internal static field H_i = H − N_zM_s + H_K; pole of the Polder tensor.'),
      kpi('fm', 'Magnetisation frequency', wm / (2e9 * PI), 'GHz', 'γµ₀M_s/2π: sets the strength of the gyrotropic response and the width of the negative-permeability band.'),
      kpi('neg_lo', 'Negative-µ band lower edge', (bV ? bV[0] : wV) / (2e9 * PI), 'GHz', 'Onset of Re µ_eff < 0 for a wave propagating across the bias field (antiresonance √(ω₀(ω₀ + ω_m)) in the lossless limit).', !!bV),
      kpi('neg_hi', 'Negative-µ band upper edge', (bV ? bV[1] : wAR) / (2e9 * PI), 'GHz', 'Upper end of the negative band, ω₀ + ω_m in the lossless limit; above it the ferrite is again transparent.', !!bV),
      kpi('mu_min', 'Most negative Re µ_eff', Tn[imin].muV.re, '–', `Minimum of the real effective permeability with loss, reached at ${f3(fl[imin])} GHz.`, Tn[imin].muV.re < 0),
      kpi('mu_re', 'Re µ_eff at probe', tp.muV.re, '–', `Real effective (Voigt) permeability (µ² − κ²)/µ at ${f3(inp.f_rel)} × the lower band edge.`, tp.muV.re < 0),
      kpi('mu_im', 'Im µ_eff at probe', tp.muV.im, '–', 'Magnetic loss at the probe frequency; non-negative for a passive magnet.', tp.muV.im >= -1e-9),
      kpi('mu_c', 'Composite Re µ at probe', tp.muC.re, '–', `Permeability of a laminate with magnetic fill fraction ${ff} for RF field in the layer plane, 1 + f(µ_eff − 1).`, tp.muC.re < 0),
      kpi('kappa', 'Off-diagonal κ at probe', tp.ka.re, '–', 'Gyrotropic (non-reciprocal) element of the Polder tensor; it is the source of Faraday rotation, isolation and circulation.'),
      kpi('dH', 'FMR linewidth', dH, 'A/m', 'Full width 2αω/(γµ₀) of the absorption line at the band-edge frequency; measuring it gives the Gilbert damping.'),
      kpi('Q', 'Resonance quality factor', 1 / (2 * al), '–', 'Approximately 1/(2α): number of precession cycles over which the stored energy decays.', al < 0.02),
    ],
    series: [
      ser('mu', 'Effective permeability against frequency', 'frequency (GHz)', 'µ (–)', [{ name: 'Re µ_eff (Voigt)', x: fl, y: Tn.map(t => N.clamp(t.muV.re, -60, 60)) }, { name: 'Im µ_eff', x: fl, y: Tn.map(t => Math.min(t.muV.im, 60)) }, { name: 'Re µ₊ (resonant circular)', x: fl, y: Tn.map(t => N.clamp(t.muP.re, -60, 60)) }, { name: 'Re µ composite', x: fl, y: Tn.map(t => N.clamp(t.muC.re, -60, 60)) }, { name: 'zero', x: [fl[0], fl[399]], y: [0, 0] }], 'Gyromagnetic permeability (clipped at ±60). The transverse effective permeability is negative between the antiresonance and ω₀ + ω_m; the resonant circular polarisation is negative from ω₀ upward. Dilution in a composite shrinks the band.'),
      ser('kittel', 'Resonance frequency against bias field', 'bias field (kA/m)', 'frequency (GHz)', [{ name: 'Kittel mode', x: Hs.map(h => h / 1e3), y: Hs.map(h => kit(h) / 1e9) }, { name: 'operating point', x: [inp.H_bias], y: [fK / 1e9] }], 'Field tuning of the uniform mode for the chosen shape (zero where the sample is not saturated). The slope approaches γµ₀/2π ≈ 35 MHz per kA/m; the negative-permeability band moves with it.'),
    ],
    equations: [...matEqs(inp, m, W),
      eqn('formula', 'Kittel equation', 'ω = γµ₀ √{[H + (N_x − N_z)M_s + H_K]·[H + (N_y − N_z)M_s + H_K]}', `demagnetising factors (N_x, N_y, N_z) = (${Nd.map(v => f3(v)).join(', ')})`),
      eqn('constitutive', 'Polder permeability tensor', 'µ = 1 + ω_m(ω₀ − iαω)/[(ω₀ − iαω)² − ω²],  κ = ω_mω/[(ω₀ − iαω)² − ω²]', 'e^(−iωt) convention; ω₀ = γµ₀H_i, ω_m = γµ₀M_s'),
      eqn('constitutive', 'Effective permeabilities', 'µ_eff = (µ² − κ²)/µ  (wave ⟂ bias);  µ± = 1 + ω_m/(ω₀ − iαω ∓ ω)  (circular, wave ∥ bias)'),
      eqn('formula', 'Negative band (lossless)', '√(ω₀(ω₀ + ω_m)) < ω < ω₀ + ω_m'),
    ],
    checks: [
      chk('Passivity: Im µ_eff ≥ 0 and Im µ₊ ≥ 0', Tn.every(t => t.muV.im >= -1e-9 && t.muP.im >= -1e-9), `min Im µ_eff = ${f3(N.minOf(Tn.map(t => t.muV.im)))}`),
      chk('Sample saturated (internal field > 0)', H - Nd[2] * m.Ms + Hk > 0, `H_i = ${f3((H - Nd[2] * m.Ms + Hk) / 1e3)} kA/m`),
      chk('Lossless band edges reproduced', !bV || al > 0.02 || (Math.abs(bV[0] - wV) < 0.03 * wV && Math.abs(bV[1] - wAR) < 0.03 * wAR), bV ? `numerical ${f3(bV[0] / (2e9 * PI))}–${f3(bV[1] / (2e9 * PI))} GHz against analytical ${f3(wV / (2e9 * PI))}–${f3(wAR / (2e9 * PI))} GHz` : 'no negative band'),
      chk('Static limit µ(0) = 1 + M_s/H_i', Math.abs(tens(1e-3 * w0).mu.re - (1 + wm / w0)) < 1e-3 * (1 + wm / w0), `µ(0) = ${f4(1 + wm / w0)}`), chk('Negative-permeability band present', !!bV, bV ? `width ${f3((bV[1] - bV[0]) / (2e9 * PI))} GHz` : 'none'),
    ],
    warnings: W,
    recommendations: [
      rec('high', bV ? `Use ${f3(bV[0] / (2e9 * PI))}–${f3(bV[1] / (2e9 * PI))} GHz for negative-µ operation and tune it with the bias field (about ${f3(GAMMA * MU0 / (2 * PI) / 1e6 * 1e3)} MHz per kA/m); pair with a wire array or a negative-ε semiconductor for a tunable negative-index medium.` : 'Choose a lower-damping ferrite (YIG, α ≈ 10⁻⁴–10⁻³) or raise the bias to obtain a negative band.', 'Gyromagnetic negative permeability above FMR'),
      rec('medium', `Loss at the probe is Im µ = ${f3(tp.muV.im)}; operate towards the upper band edge where |µ| is smaller but the loss tangent is lowest.`, 'Lorentz-type dispersion of the Polder tensor'),
      rec('medium', `The ${SHAPES[inp.shape][0].toLowerCase()} needs ${f3(Math.max(Nd[2] * m.Ms - Hk, 0) / 1e3)} kA/m just to overcome its own demagnetising field; in-plane magnetised films or axially magnetised rods need the smallest bias magnets.`, 'Shape anisotropy in the Kittel equation'),
      rec('low', 'Self-biased hexaferrites or nanowire arrays with large shape anisotropy remove the permanent magnet, saving rare-earth material and volume.', 'Bias-magnet elimination', true),
    ],
    publish: { 'eff.mu': tp.muV.re, 'res.f0': fK, 'res.Q': 1 / (2 * al), 'res.bandgap_lo': (bV ? bV[0] : wV) / (2 * PI), 'res.bandgap_hi': (bV ? bV[1] : wAR) / (2 * PI), 'prop.Ms': m.Ms },
    numerics: { method: 'Closed-form Kittel and Polder (linearised LLG) permeability on a 400-point frequency sweep', dof: ws.length, notes: 'Uniform precession (no spin waves or magnetostatic modes), saturated single-domain state, small-signal response.' },
  };
}

// ---------------------------------------------------------------- spintronics
function spin(inp) {
  const W = [], T = inp.T + 273.15, kT = K_B * T, P1 = inp.P1, P2 = inp.P2, tmr0 = 2 * P1 * P2 / (1 - P1 * P2), tmrV = v => tmr0 / (1 + (v / inp.V_half) ** 2), Vr = inp.V_read, tmr = tmrV(Vr);
  const d = inp.mtj_d * 1e-9, tf = inp.t_free * 1e-9, area = PI * d * d / 4, V = area * tf, Ms = inp.Ms_free, Keff = inp.K_i * 1e-3 / tf - 0.5 * MU0 * Ms * Ms, Delta = Math.max(Keff, 0) * V / kT, ret = 1e-9 * Math.exp(Math.min(Delta, 600)), yrs = ret / 3.156e7;
  const Bk = 2 * Math.max(Keff, 0) / Ms, Jc0 = 2 * Q_E * inp.alpha * Ms * tf * Bk / (HBAR * Math.max(P1, 0.01)), Ic0 = Jc0 * area, RP = inp.RA * 1e-12 / area, RAP = RP * (1 + tmr), dV = Vr * tmr / (2 + tmr), Vc = Ic0 * RP, Ew = Ic0 * Ic0 * RP * 10e-9;
  // CPP-GMR two-current series-resistor model (per unit area, Ω·m²)
  const be = inp.beta, rhoS = inp.rho_F * 1e-8 / (1 - be * be), tF = inp.t_F * 1e-9, tN = inp.t_N * 1e-9, rU = 2 * rhoS * (1 - be) * tF, rD = 2 * rhoS * (1 + be) * tF, rN = 2 * inp.rho_N * 1e-8 * tN;
  const Ds = 1.57e6 * inp.mfp * 1e-9 / 3, lsf = Math.sqrt(Ds * inp.tau_sf * 1e-12), att = Math.exp(-tN / lsf), RPg = (2 * rU + rN) * (2 * rD + rN) / (2 * (rU + rD + rN)), RAPg0 = (rU + rD + rN) / 2, RAPg = RPg + (RAPg0 - RPg) * att, gmr = (RAPg - RPg) / RPg;
  if (Keff <= 0) W.push('Interfacial anisotropy does not overcome the demagnetising energy: the free layer is in-plane magnetised and the perpendicular stability factor is zero. Thin the layer or raise K_i.');
  if (Vr > inp.V_half) W.push('Read voltage exceeds V_half: more than half of the zero-bias TMR is lost to inelastic tunnelling.');
  if (Vc > 1) W.push(`The critical switching voltage ${f3(Vc)} V approaches the breakdown voltage of an MgO barrier (about 1–1.5 V): lower RA or the critical current.`);
  if (tN > lsf) W.push('The spacer is thicker than its spin-diffusion length: most of the spin accumulation relaxes before reaching the second ferromagnet.');
  const Vs = N.linspace(0, 1, 80), ds = N.linspace(10, 100, 91), tNs = N.logspace(0.5, 2000, 80);
  return {
    summary: `Julliere TMR for P = ${P1}/${P2} is ${f3(100 * tmr0)} % at zero bias and ${f3(100 * tmr)} % at ${Vr} V (read signal ${f3(dV * 1e3)} mV). A ${inp.mtj_d} nm, ${inp.t_free} nm perpendicular free layer has K_eff = ${f3(Keff)} J/m³, thermal stability Δ = ${f3(Delta)} (retention ${yrs > 1e30 ? '> 10³⁰' : f3(yrs)} years) and a critical current of ${f3(Ic0 * 1e6)} µA (${f3(Jc0)} A/m²). The CPP spin valve gives GMR = ${f3(100 * gmr)} % with a spacer spin-diffusion length of ${f3(lsf * 1e9)} nm.`,
    kpis: [
      kpi('TMR0', 'Zero-bias TMR (Julliere)', 100 * tmr0, '%', 'Optimistic-definition tunnelling magnetoresistance (R_AP − R_P)/R_P = 2P₁P₂/(1 − P₁P₂) from the spin polarisations of the two electrodes.'),
      kpi('TMR', 'TMR at read voltage', 100 * tmr, '%', `TMR reduced by the empirical bias dependence 1/(1 + (V/V_half)²) at ${Vr} V.`, tmr > 1),
      kpi('RP', 'Parallel-state resistance', RP, 'Ω', `RA product ${inp.RA} Ω·µm² divided by the junction area; it must match the selector/transistor impedance.`),
      kpi('RAP', 'Antiparallel-state resistance', RAP, 'Ω', 'Junction resistance with the two magnetisations antiparallel, R_P(1 + TMR), evaluated at the read voltage.'),
      kpi('dV', 'Read signal', dV * 1e3, 'mV', 'Voltage difference between the two states when read against a mid-point reference at the read voltage; sense amplifiers typically need more than about 20–50 mV.', dV > 0.03),
      kpi('Keff', 'Effective perpendicular anisotropy', Keff, 'J/m³', 'K_i/t − µ₀M_s²/2: interfacial anisotropy per thickness minus thin-film shape anisotropy; positive values give an out-of-plane easy axis.', Keff > 0),
      kpi('Delta', 'Thermal stability factor Δ', Delta, '–', 'K_eff·V/kT of the free layer treated as a single domain; about 60 is needed for ten-year retention of a large array.', Delta >= 60),
      kpi('ret', 'Retention time', Math.min(yrs, 1e30), 'years', 'Néel–Arrhenius lifetime τ₀·exp(Δ) of a single bit with τ₀ = 1 ns (capped at 10³⁰ years in the display).', yrs > 10),
      kpi('Ic0', 'Critical switching current', Ic0 * 1e6, 'µA', 'Zero-temperature spin-transfer-torque threshold (4eα/ħP)·Δ·kT for a perpendicular free layer.'),
      kpi('Jc0', 'Critical current density', Jc0, 'A/m²', 'Threshold current divided by junction area; around 10¹⁰ A/m² is typical of CoFeB/MgO junctions.', Jc0 < 5e10),
      kpi('eff', 'Switching efficiency Δ/I_c0', Ic0 > 0 ? Delta / (Ic0 * 1e6) : 0, '1/µA', 'Stability obtained per unit write current; it depends only on damping and polarisation in the macrospin picture (ħP/4eαkT).'),
      kpi('Ew', 'Write energy (10 ns at I_c0)', Ew * 1e15, 'fJ', 'I²R·t estimate for a 10 ns pulse at the threshold current in the parallel state; practical overdrive raises it several-fold.'),
      kpi('GMR', 'CPP-GMR ratio', 100 * gmr, '%', 'Resistance change of a ferromagnet/spacer/ferromagnet trilayer from the two-current series-resistor model, attenuated by spin relaxation in the spacer.'),
      kpi('lsf', 'Spin-diffusion length of the spacer', lsf * 1e9, 'nm', `√(D·τ_sf) with D = v_Fλ/3 for a mean free path of ${inp.mfp} nm and spin-flip time ${inp.tau_sf} ps; spin information decays over this length.`, tN < lsf),
    ],
    series: [
      ser('tmrV', 'TMR against bias voltage', 'bias (V)', 'TMR (%)', [{ name: 'TMR(V)', x: Vs, y: Vs.map(v => 100 * tmrV(v)) }, { name: 'read point', x: [Vr], y: [100 * tmr] }], 'Bias dependence of the tunnelling magnetoresistance. The signal voltage V·TMR/(2 + TMR) passes through a maximum near V_half, which guides the choice of read voltage.'),
      ser('delta', 'Thermal stability against junction diameter', 'diameter (nm)', 'Δ (–)', [{ name: 'Δ = K_eff V/kT', x: ds, y: ds.map(x => Math.max(Keff, 0) * PI * (x * 1e-9) ** 2 / 4 * tf / kT) }, { name: 'ten-year target', x: [10, 100], y: [60, 60] }], 'Macrospin stability factor against diameter. Real junctions above about 30–40 nm reverse by domain-wall nucleation, so Δ saturates instead of following the area law.'),
      ser('gmr', 'CPP-GMR against spacer thickness', 'spacer thickness (nm)', 'GMR (%)', [{ name: 'GMR', x: tNs, y: tNs.map(t => { const rn = 2 * inp.rho_N * 1e-8 * t * 1e-9, rp = (2 * rU + rn) * (2 * rD + rn) / (2 * (rU + rD + rn)); return 100 * ((rU + rD + rn) / 2 - rp) * Math.exp(-t * 1e-9 / lsf) / rp; }) }], 'Decay of the spin-valve signal with spacer thickness from added spin-independent resistance and spin relaxation over the spin-diffusion length.', { logx: true }),
    ],
    equations: [
      eqn('formula', 'Julliere model', 'TMR = (R_AP − R_P)/R_P = 2P₁P₂/(1 − P₁P₂);  TMR(V) ≈ TMR₀/(1 + (V/V_half)²)'), eqn('formula', 'Thermal stability and retention', 'Δ = K_eff V/kT,  K_eff = K_i/t − µ₀M_s²/2,  τ = τ₀ e^Δ'),
      eqn('formula', 'Spin-transfer-torque threshold', 'I_c0 = (4eα/ħP)·Δ·kT'), eqn('constitutive', 'Two-current series-resistor model', 'r↑,↓ = 2ρ*_F(1 ∓ β)t_F,  r_N = 2ρ_N t_N;  R_P = (2r↑ + r_N)(2r↓ + r_N)/[2(r↑ + r↓ + r_N)],  R_AP = (r↑ + r↓ + r_N)/2'),
      eqn('governing', 'Spin diffusion equation', 'd²µ_s/dx² = µ_s/l_sf²,  l_sf = √(D τ_sf)'),
      eqn('formula', 'Typical stack', (() => { try { return `CoFeB / ${chem.prettyFormula('MgO')} (M = ${f4(chem.molarMass('MgO'))} g/mol) / CoFeB tunnel junction; Co/Cu/Co spin valve (${chem.prettyFormula('Cu')}, M = ${f4(chem.molarMass('Cu'))} g/mol spacer)`; } catch { return 'CoFeB/MgO/CoFeB'; } })()),
    ],
    checks: [
      chk('Polarisations physical (|P| < 1)', P1 < 1 && P2 < 1, `P₁P₂ = ${f3(P1 * P2)}`), chk('Ten-year retention (Δ ≥ 60)', Delta >= 60, `Δ = ${f3(Delta)}`), chk('Perpendicular easy axis (K_eff > 0)', Keff > 0, `K_eff = ${f3(Keff)} J/m³`),
      chk('Write voltage below barrier breakdown (≈ 1 V)', Vc < 1, `I_c0·R_P = ${f3(Vc)} V`), chk('R_AP > R_P in both models', RAP > RP && RAPg >= RPg, `GMR ${f3(100 * gmr)} %, TMR ${f3(100 * tmr)} %`),
      chk('Macrospin identity I_c0 = (4eα/ħP)ΔkT', Keff <= 0 || Math.abs(Ic0 - 4 * Q_E * inp.alpha * Delta * kT / (HBAR * Math.max(P1, 0.01))) < 1e-6 * Ic0, 'threshold from fields and from the stability factor agree'),
    ],
    warnings: W,
    recommendations: [
      rec('high', Delta < 60 ? `Δ = ${f3(Delta)} is below the retention target: increase the diameter to about ${f3(inp.mtj_d * Math.sqrt(60 / Math.max(Delta, 1e-6)))} nm, use a double CoFeB/MgO interface (doubling K_i), or thin the free layer.` : `Δ = ${f3(Delta)} meets ten-year retention; the write current of ${f3(Ic0 * 1e6)} µA scales with Δ, so do not over-design stability.`, 'Retention–write-current trade-off'),
      rec('medium', `Read near ${f3(Math.min(Vr, inp.V_half))} V: the signal is ${f3(dV * 1e3)} mV; keep the read current well below I_c0 (now ${f3(100 * (Vr / RP) / Math.max(Ic0, 1e-30))} % of it) to avoid read disturb.`, 'Bias dependence of TMR and read-disturb margin'),
      rec('medium', `Lower damping reduces write current proportionally (α = ${inp.alpha}); avoid heavy-metal capping layers adjacent to the free layer that add spin pumping.`, 'I_c0 ∝ α'),
      rec('low', `Non-volatile magnetic memory removes standby power: at ${f3(Ew * 1e15)} fJ per write the energy cost is paid only on writes, favouring write-infrequent workloads. Tantalum, cobalt and platinum-group layers are critical materials: keep seed and pinning layers thin.`, 'Energy and critical-material footprint of the stack', true),
    ],
    publish: { 'mag.TMR': tmr0, 'mag.delta': Delta, 'elec.current_density': Jc0, 'mag.GMR': gmr, 'mag.l_sf': lsf, 'life.years': Math.min(yrs, 1e30) },
    numerics: { method: 'Closed-form Julliere, Néel–Arrhenius, Slonczewski threshold and two-current resistor network', dof: 1, notes: 'Macrospin free layer, thin-film demagnetising factor N_z = 1, bias dependence of TMR empirical.' },
  };
}

const M = (...a) => ({ mode: a });
export default {
  id: 13,
  tagline: 'Hysteresis, magnetisation dynamics, superparamagnetism and hyperthermia, ferromagnetic resonance with negative permeability, and spintronic device figures.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'spm', group: 'Model', help: 'Selects the magnetic model that is solved; each analysis has its own inputs.', options: [
      { value: 'spm', label: 'Superparamagnetic nanoparticles and hyperthermia' }, { value: 'sw', label: 'Stoner–Wohlfarth hysteresis' }, { value: 'llg', label: 'Macrospin LLG switching (field / spin torque)' },
      { value: 'fmr', label: 'Ferromagnetic resonance and negative permeability' }, { value: 'spin', label: 'Spintronics: TMR, GMR, MTJ stability' }] },
    { key: 'material', label: 'Magnetic material', type: 'select', value: 'Fe3O4', group: 'Material', showIf: M('spm', 'sw', 'llg', 'fmr'), help: 'Preset saturation magnetisation, anisotropy constant, exchange stiffness and density (room-temperature handbook values); choose "custom" to use the numeric inputs.', options: Object.entries(MAG).map(([value, m]) => ({ value, label: m.label })) },
    { key: 'Ms', label: 'Saturation magnetisation (custom)', unit: 'A/m', value: 4.8e5, min: 1e4, max: 2e6, log: true, group: 'Material', bus: 'prop.Ms', showIf: M('spm', 'sw', 'llg', 'fmr'), help: 'Used only when the material is "custom". Sets moment per volume, demagnetising energy and the gyromagnetic bandwidth.' },
    { key: 'K', label: 'Uniaxial anisotropy constant (custom)', unit: 'J/m³', value: 1.35e4, min: 100, max: 1e7, log: true, group: 'Material', showIf: M('spm', 'sw', 'llg', 'fmr'), help: 'Used only when the material is "custom". Sets the energy barrier KV and the anisotropy field 2K/µ₀M_s.' },
    { key: 'A_ex', label: 'Exchange stiffness (custom)', unit: 'pJ/m', value: 12, min: 0.5, max: 50, group: 'Material', showIf: M('spm', 'sw'), help: 'Used only when the material is "custom". Enters wall width, exchange length and the single-domain size.' },
    { key: 'rho', label: 'Density (custom)', unit: 'kg/m³', value: 5180, min: 1000, max: 22000, group: 'Material', bus: 'prop.density', showIf: M('spm', 'sw'), help: 'Used only when the material is "custom". Converts volumetric loss to specific loss power.' },
    { key: 'T', label: 'Temperature', unit: '°C', value: 37, min: -268, max: 500, group: 'Conditions', bus: 'proc.T_C', help: 'Sets the thermal energy kT; the temperature dependence of M_s and K is not included.' },
    { key: 'd_nm', label: 'Particle core diameter', unit: 'nm', value: 16, min: 2, max: 500, log: true, group: 'Particle', showIf: M('spm', 'sw', 'llg'), help: 'Median magnetic core diameter (sphere). The barrier KV grows with its cube.' },
    { key: 'sigma_ln', label: 'Log-normal width σ', unit: '–', value: 0.15, min: 0, max: 0.6, group: 'Particle', showIf: M('spm'), help: 'Standard deviation of ln(diameter); 0.1 is nearly monodisperse, 0.3 typical of co-precipitation.' },
    { key: 'tau_m', label: 'Measurement time', unit: 's', value: 100, min: 1e-8, max: 3e9, log: true, group: 'Conditions', showIf: M('spm', 'sw'), help: 'Time window of the experiment (100 s magnetometry, 10⁻⁸ s Mössbauer, 3×10⁸ s storage).' },
    { key: 'coat', label: 'Hydrodynamic shell thickness', unit: 'nm', value: 3, min: 0, max: 100, group: 'Particle', showIf: M('spm'), help: 'Surfactant or polymer coating added to the core radius for Brownian rotation.' },
    { key: 'visc', label: 'Carrier viscosity', unit: 'mPa·s', value: 1, min: 0.2, max: 1e4, log: true, group: 'Conditions', showIf: M('spm'), help: 'Viscosity of the carrier fluid (0.7 water at 37 °C; very high values represent immobilised particles in tissue or gels).' },
    { key: 'phi_v', label: 'Particle volume fraction', unit: 'vol %', value: 0.1, min: 0.001, max: 20, log: true, group: 'Particle', showIf: M('spm'), help: 'Magnetic-core volume fraction of the suspension (0.1 vol % of magnetite ≈ 5 mg/mL).' },
    { key: 'H_ac', label: 'AC field amplitude', unit: 'kA/m', value: 3, min: 0.1, max: 100, log: true, group: 'Alternating field', showIf: M('spm'), help: 'Amplitude of the alternating field; linear-response heating scales with its square.' },
    { key: 'f_ac', label: 'AC field frequency', unit: 'kHz', value: 500, min: 1, max: 2000, log: true, group: 'Alternating field', showIf: M('spm', 'sw'), help: 'Frequency of the alternating field (100–500 kHz for hyperthermia).' },
    { key: 'H_dc', label: 'Maximum DC field for M(H)', unit: 'kA/m', value: 400, min: 1, max: 5000, log: true, group: 'Alternating field', showIf: M('spm'), help: 'Range of the plotted equilibrium magnetisation curve.' },
    { key: 'psi', label: 'Field angle to easy axis', unit: '°', value: 30, min: 0, max: 90, group: 'Hysteresis', showIf: M('sw'), help: 'Angle between the applied field and the particle easy axis for the single-particle loop.' },
    { key: 'hmax', label: 'Maximum field H_max/H_K', unit: '–', value: 1.5, min: 0.2, max: 3, group: 'Hysteresis', showIf: M('sw'), help: 'Sweep amplitude in units of the anisotropy field; values below 1 give minor loops.' },
    { key: 'alpha', label: 'Gilbert damping α', unit: '–', value: 0.01, min: 1e-4, max: 1, log: true, group: 'Dynamics', showIf: M('llg', 'fmr', 'spin'), help: 'Dimensionless damping (10⁻⁴ YIG, 0.005–0.01 CoFeB and permalloy, 0.05–0.1 ferrites and nanoparticles).' },
    { key: 'B_pulse', label: 'Pulse field µ₀H', unit: 'mT', value: 60, min: 0, max: 20000, group: 'Dynamics', showIf: M('llg'), help: 'Amplitude of the reversing field pulse.' },
    { key: 'th_H', label: 'Pulse angle from reversed easy axis', unit: '°', value: 30, min: 0, max: 90, group: 'Dynamics', showIf: M('llg'), help: 'Tilt of the pulse field away from the antiparallel direction; 45° minimises the switching field.' },
    { key: 't_pulse', label: 'Pulse duration', unit: 'ns', value: 2, min: 0.001, max: 100, log: true, group: 'Dynamics', showIf: M('llg'), help: 'Length of the field pulse; the field is zero afterwards.' },
    { key: 't_end', label: 'Simulated time', unit: 'ns', value: 30, min: 0.01, max: 200, log: true, group: 'Dynamics', showIf: M('llg'), help: 'Total integration time (automatically shortened to 150 precession periods).' },
    { key: 'J_stt', label: 'Spin-polarised current density', unit: 'A/m²', value: 0, min: 0, max: 1e13, group: 'Dynamics', bus: 'elec.current_density', showIf: M('llg'), help: 'Current density through the free layer with polariser antiparallel to the initial state (0 switches spin torque off).' },
    { key: 'P_spin', label: 'Spin polarisation of the current', unit: '–', value: 0.6, min: 0.05, max: 1, group: 'Dynamics', showIf: M('llg'), help: 'Spin-torque efficiency of the polarising layer.' },
    { key: 't_free', label: 'Free-layer thickness', unit: 'nm', value: 1.3, min: 0.5, max: 20, group: 'Dynamics', showIf: M('llg', 'spin'), help: 'Thickness of the switching layer; spin torque and interfacial anisotropy both scale as 1/t.' },
    { key: 'H_bias', label: 'Static bias field', unit: 'kA/m', value: 300, min: 1, max: 3000, log: true, group: 'Resonance', showIf: M('fmr'), help: 'Applied static field along z that saturates the sample and sets the resonance.' },
    { key: 'shape', label: 'Sample shape', type: 'select', value: 'sphere', group: 'Resonance', showIf: M('fmr'), help: 'Demagnetising geometry entering the Kittel equation.', options: Object.entries(SHAPES).map(([value, s]) => ({ value, label: s[0] })) },
    { key: 'f_fill', label: 'Magnetic fill fraction', unit: '–', value: 1, min: 0.01, max: 1, group: 'Resonance', showIf: M('fmr'), help: 'Volume fraction of magnetic material in a laminate metamaterial (1 = solid ferrite).' },
    { key: 'f_rel', label: 'Probe frequency (relative)', unit: '–', value: 1.1, min: 0.2, max: 3, group: 'Resonance', showIf: M('fmr'), help: 'Probe frequency as a multiple of the lower edge of the negative-permeability band.' },
    { key: 'P1', label: 'Spin polarisation, electrode 1', unit: '–', value: 0.6, min: 0.05, max: 0.99, group: 'Tunnel junction', showIf: M('spin'), help: 'Tunnelling spin polarisation (about 0.4–0.5 Co/Fe with AlOx, effectively 0.6–0.85 for CoFeB/MgO).' },
    { key: 'P2', label: 'Spin polarisation, electrode 2', unit: '–', value: 0.6, min: 0.05, max: 0.99, group: 'Tunnel junction', showIf: M('spin'), help: 'Tunnelling spin polarisation of the second electrode.' },
    { key: 'V_read', label: 'Read voltage', unit: 'V', value: 0.15, min: 0.01, max: 1, group: 'Tunnel junction', showIf: M('spin'), help: 'Bias applied during the read operation.' },
    { key: 'V_half', label: 'TMR half-voltage', unit: 'V', value: 0.5, min: 0.1, max: 2, group: 'Tunnel junction', showIf: M('spin'), help: 'Bias at which TMR falls to half its zero-bias value (0.4–0.7 V for MgO junctions).' },
    { key: 'RA', label: 'Resistance–area product', unit: 'Ω·µm²', value: 10, min: 0.5, max: 1e4, log: true, group: 'Tunnel junction', showIf: M('spin'), help: 'Barrier RA in the parallel state (5–20 Ω·µm² for STT-MRAM).' },
    { key: 'mtj_d', label: 'Junction diameter', unit: 'nm', value: 40, min: 8, max: 300, log: true, group: 'Tunnel junction', showIf: M('spin'), help: 'Diameter of the circular free layer.' },
    { key: 'K_i', label: 'Interfacial anisotropy', unit: 'mJ/m²', value: 1.3, min: 0.1, max: 5, group: 'Tunnel junction', showIf: M('spin'), help: 'Perpendicular interface anisotropy of CoFeB/MgO (about 1.3 mJ/m² per interface).' },
    { key: 'Ms_free', label: 'Free-layer saturation magnetisation', unit: 'A/m', value: 1.0e6, min: 2e5, max: 2e6, group: 'Tunnel junction', showIf: M('spin'), help: 'M_s of the free layer (about 1.0–1.3×10⁶ A/m for CoFeB).' },
    { key: 'rho_F', label: 'Ferromagnet resistivity', unit: 'µΩ·cm', value: 15, min: 1, max: 300, log: true, group: 'Spin valve', showIf: M('spin'), help: 'Resistivity of the ferromagnetic layers of the CPP spin valve.' },
    { key: 'beta', label: 'Bulk spin asymmetry β', unit: '–', value: 0.5, min: 0, max: 0.95, group: 'Spin valve', showIf: M('spin'), help: 'Spin asymmetry of the ferromagnet resistivity (about 0.4–0.7 for Co, CoFe, NiFe).' },
    { key: 't_F', label: 'Ferromagnet thickness', unit: 'nm', value: 3, min: 0.5, max: 50, group: 'Spin valve', showIf: M('spin'), help: 'Thickness of each ferromagnetic layer.' },
    { key: 'rho_N', label: 'Spacer resistivity', unit: 'µΩ·cm', value: 3, min: 0.5, max: 100, log: true, group: 'Spin valve', showIf: M('spin'), help: 'Resistivity of the non-magnetic spacer (thin-film Cu about 3 µΩ·cm).' },
    { key: 't_N', label: 'Spacer thickness', unit: 'nm', value: 3, min: 0.5, max: 2000, log: true, group: 'Spin valve', showIf: M('spin'), help: 'Thickness of the non-magnetic spacer.' },
    { key: 'mfp', label: 'Spacer mean free path', unit: 'nm', value: 20, min: 0.5, max: 200, log: true, group: 'Spin valve', showIf: M('spin'), help: 'Electron mean free path in the spacer; with the Fermi velocity it gives the diffusion constant.' },
    { key: 'tau_sf', label: 'Spin-flip time', unit: 'ps', value: 10, min: 0.01, max: 1e4, log: true, group: 'Spin valve', showIf: M('spin'), help: 'Spin relaxation time in the spacer (picoseconds in Cu at room temperature, far shorter in Pt).' },
  ],
  presets: [
    { name: 'Magnetite nanoparticles for hyperthermia (16 nm, 500 kHz)', note: 'Fe₃O₄ cores below the superparamagnetic limit heated in a low-amplitude field where linear-response theory holds.', values: { mode: 'spm', material: 'Fe3O4', d_nm: 16, sigma_ln: 0.15, T: 37, H_ac: 3, f_ac: 500, visc: 0.7, coat: 3, phi_v: 0.1 } },
    { name: 'Cobalt ferrite nanoparticles: blocked, hard magnetic', note: 'High anisotropy makes 12 nm CoFe₂O₄ blocked at room temperature.', values: { mode: 'spm', material: 'CoFe2O4', d_nm: 12, sigma_ln: 0.2, T: 25, H_ac: 10, f_ac: 300 } },
    { name: 'FePt recording grains: Stoner–Wohlfarth loop', note: 'L1₀ FePt grains of 6 nm are thermally stable with very high coercivity.', values: { mode: 'sw', material: 'FePt', d_nm: 6, psi: 30, hmax: 1.5, T: 25, tau_m: 3e8 } },
    { name: 'YIG ferrite metamaterial: negative permeability above FMR', note: 'Low-damping garnet biased at 300 kA/m; µ < 0 between antiresonance and ω₀ + ω_m.', values: { mode: 'fmr', material: 'YIG', H_bias: 300, shape: 'sphere', alpha: 0.001, f_fill: 1, f_rel: 1.1 } },
    { name: 'CoFeB/MgO STT-MRAM bit (40 nm)', note: 'Perpendicular magnetic tunnel junction: TMR, retention and switching current.', values: { mode: 'spin', P1: 0.65, P2: 0.65, mtj_d: 40, t_free: 1.3, K_i: 1.3, Ms_free: 1.0e6, alpha: 0.006, RA: 10, T: 25 } },
    { name: 'Precessional field switching of a cobalt nanoparticle', note: 'Macrospin LLG reversal by a tilted field pulse.', values: { mode: 'llg', material: 'Co', d_nm: 10, alpha: 0.02, B_pulse: 500, th_H: 30, t_pulse: 0.5, t_end: 1.5, T: 25 } },
  ],
  run(inp) {
    const f = { spm, sw, llg, fmr, spin }[inp.mode] || spm;
    return finish(f(inp));
  },
};
