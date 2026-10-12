// Suite 15 — Acoustic, Phononic and Elastic-Wave Engineering.
// Locally resonant (mass-in-mass) lattices with negative effective mass, Helmholtz-resonator and membrane ducts with negative
// bulk modulus / density, layered phononic crystals, sound-insulating and absorbing panels, Christoffel waves in crystals and
// flexural resonances of nanobeams.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { PI, C, K_B } = N;
const f3 = v => N.fmt(v, 3), f4 = v => N.fmt(v, 4);
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const ser = (key, title, xlabel, ylabel, lines, explain, o = {}) => ({ key, title, xlabel, ylabel, logx: !!o.logx, logy: !!o.logy, lines, explain });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const pos = v => (v > 1e-300 ? v : 1e-300), dB = v => 20 * Math.log10(pos(v));
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
const ONE = { re: 1, im: 0 }, J = { re: 0, im: 1 };
const up = z => (z.im < 0 || (z.im === 0 && z.re < 0) ? { re: -z.re, im: -z.im } : z), csqrt = z => up(C.sqrt(z));
const cacos = z => { if (!(C.abs(z) < 1e7)) return { re: Math.abs(Math.atan2(z.im, z.re)), im: Math.log(Math.min(2 * C.abs(z), 1e300)) }; const u = C.add(z, C.mul(J, C.sqrt(C.sub(ONE, C.mul(z, z))))); return { re: Math.abs(Math.atan2(u.im, u.re)), im: Math.abs(Math.log(pos(C.abs(u)))) }; };
const bands = (xs, flag) => { const out = []; let s = -1; for (let i = 0; i <= xs.length; i++) { const on = i < xs.length && flag[i]; if (on && s < 0) s = i; if (!on && s >= 0) { out.push([s, i - 1]); s = -1; } } return out; };
const formulaEq = (f, label) => { try { return [eqn('formula', label, `${chem.prettyFormula(f)}, M = ${f4(chem.molarMass(f))} g/mol`)]; } catch { return []; } };
function kkMismatch(chi, w, ws) { // Re χ(ω) = (2/π) P∫ ω′ Im χ/(ω′² − ω²) dω′ with singularity subtraction
  const n = 20000, a = w * chi(w).im; let s = 0;
  for (let i = 0; i < n; i++) { const t = (i + 0.5) * PI / (2 * n), wp = ws * Math.tan(t), d = wp * wp - w * w; if (Math.abs(d) < 1e-12 * w * w) continue; s += (wp * chi(wp).im - a) / d * ws / (Math.cos(t) ** 2); }
  return Math.abs((2 / PI) * s * PI / (2 * n) - chi(w).re) / Math.max(Math.abs(chi(w).re), 1e-12);
}
// homogeneous slab between identical half-spaces: n·k₀L = φ (complex phase), z = relative impedance
function slabT(phi, z) { const R = C.div(C.sub(z, ONE), C.add(z, ONE)), e = C.exp(C.mul(J, phi)), e2 = C.mul(e, e), den = C.sub(ONE, C.mul(C.mul(R, R), e2)); return { t: C.div(C.mul(C.sub(ONE, C.mul(R, R)), e), den), r: C.div(C.mul(R, C.sub(ONE, e2)), den) }; }

// ================================================================ negative effective mass: mass-in-mass lattice
function massInMass(inp) {
  const W = [], m1 = inp.m1, m2 = inp.m2, k1 = inp.k1, k2 = inp.k2, eta = inp.eta, a = inp.a_cell * 1e-3, Nc = Math.round(inp.Ncell), w0 = Math.sqrt(k2 / m2), f0 = w0 / (2 * PI), Mt = m1 + m2;
  const meff = (w, e = eta) => { const kc = { re: k2, im: k2 * e }; return C.add(C.of(m1), C.div(C.scale(kc, m2), C.sub(kc, C.of(m2 * w * w)))); };
  const qa = (w, e = eta) => cacos(C.sub(ONE, C.div(C.scale(meff(w, e), w * w), { re: 2 * k1, im: 2 * k1 * e })));
  // lossless band edges: m₁x² − (M ω₀² + 4k₁)x + 4k₁ω₀² = 0, x = ω²
  const bq = Mt * w0 * w0 + 4 * k1, disc = Math.sqrt(Math.max(bq * bq - 16 * m1 * k1 * w0 * w0, 0)), wLo = Math.sqrt((bq - disc) / (2 * m1)), wOpt = Math.sqrt((bq + disc) / (2 * m1)), wHi = w0 * Math.sqrt(1 + m2 / m1);
  const wmax = Math.max(1.25 * wOpt, 2.2 * wHi), ws = N.linspace(wmax / 400, wmax, 400), fl = ws.map(w => w / (2 * PI)), Q = ws.map(w => qa(w)), Q0 = ws.map(w => qa(w, 0)), Me = ws.map(w => meff(w));
  // finite chain: base excitation u₀ = 1, dynamic-stiffness (tridiagonal) solve for N cells, free far end
  const chain = w => { const kc = { re: k1, im: k1 * eta }, D = C.sub(C.scale(kc, 2), C.scale(meff(w), w * w)), Dl = C.sub(kc, C.scale(meff(w), w * w)), cp = [], dp = []; let den = Nc === 1 ? Dl : D; cp.push(C.div(C.scale(kc, -1), den)); dp.push(C.div(kc, den));
    for (let j = 1; j < Nc; j++) { const dj = j === Nc - 1 ? Dl : D; den = C.add(dj, C.mul(kc, cp[j - 1])); cp.push(C.div(C.scale(kc, -1), den)); dp.push(C.div(C.mul(kc, dp[j - 1]), den)); } return C.abs(dp[Nc - 1]); };
  const Tr = ws.map(chain), inGap = ws.map(w => w > wLo && w < wHi), iGap = ws.map((w, i) => i).filter(i => inGap[i]), iMin = iGap.length ? iGap[N.argmin(iGap.map(i => Tr[i]))] : 0, attMax = N.maxOf(Q.map(q => q.im)) * 8.686, iAtt = N.argmax(Q.map(q => q.im));
  const wpr = inp.f_rel * w0, mp = meff(wpr), qp = qa(wpr), cL = a * Math.sqrt(k1 / Mt), rhoDyn = mp.re / a ** 3, negMin = N.minOf(Me.map(m => m.re)), pred = iGap.length ? -8.686 * Nc * Q[iMin].im : 0;
  if (wLo >= w0) W.push('The Bragg cut-off of the host chain lies below the local resonance: the gap is then a merged Bragg/resonance gap and the negative-mass band lies entirely above the pass band.');
  if (eta > 0.2) W.push('Loss factor above 0.2 smears the resonance: the effective mass may no longer become negative.');
  if (2 * PI * cL / w0 < 6 * a) W.push('The host wavelength at the local resonance is shorter than about six cells: the effective-mass description is close to the Bragg regime.');
  return {
    summary: `The inner resonator (${f3(m2)} kg on ${f3(k2)} N/m) resonates at ${f3(f0)} Hz. The effective dynamic mass is negative from ${f3(f0)} to ${f3(wHi / (2 * PI))} Hz and the lattice has a band gap from ${f3(wLo / (2 * PI))} to ${f3(wHi / (2 * PI))} Hz (${f3(100 * (wHi - wLo) / (0.5 * (wHi + wLo)))} % relative width). Peak attenuation is ${f3(attMax)} dB per cell; a ${Nc}-cell chain transmits down to ${f3(dB(Tr[iMin]))} dB inside the gap. At ${f3(wpr / (2 * PI))} Hz m_eff = ${f3(mp.re)} kg against a static mass of ${f3(Mt)} kg.`,
    kpis: [
      kpi('f0', 'Local resonance frequency', f0, 'Hz', 'Natural frequency √(k₂/m₂)/2π of the inner mass on its spring; the effective mass diverges and changes sign here.'),
      kpi('neg_hi', 'Upper edge of negative-mass band', wHi / (2 * PI), 'Hz', 'Frequency f₀√(1 + m₂/m₁) where the effective mass returns through zero; between f₀ and this value the cell accelerates opposite to the applied force.'),
      kpi('gap_lo', 'Band-gap lower edge', wLo / (2 * PI), 'Hz', 'Top of the acoustic branch of the infinite lossless lattice (Bloch wavenumber π/a), from the quadratic m₁ω⁴ − (Mω₀² + 4k₁)ω² + 4k₁ω₀² = 0.'),
      kpi('gap_hi', 'Band-gap upper edge', wHi / (2 * PI), 'Hz', 'Bottom of the optical branch, where m_eff = 0; the gap between the edges exists for any lattice spacing because it comes from the local resonance, not from Bragg scattering.'),
      kpi('gap_w', 'Relative gap width', (wHi - wLo) / (0.5 * (wHi + wLo)), '–', 'Gap width divided by its centre frequency; it grows with the mass ratio m₂/m₁.'),
      kpi('m_eff', 'Effective mass at probe frequency', mp.re, 'kg', `Real part of m₁ + m₂ω₀²/(ω₀² − ω²) at ${f3(inp.f_rel)} × f₀ including damping; negative values lie in the band highlighted above.`, mp.re < 0),
      kpi('m_min', 'Most negative effective mass', negMin, 'kg', 'Minimum of Re m_eff over the sweep, reached just above the resonance; damping limits its magnitude.', negMin < 0),
      kpi('rho_dyn', 'Effective dynamic density at probe', rhoDyn, 'kg/m³', 'Effective mass divided by the volume a³ of a cubic cell of side equal to the lattice constant.', null),
      kpi('att', 'Peak attenuation per cell', attMax, 'dB/cell', `Largest imaginary Bloch wavenumber (8.686·Im(qa)), at ${f3(fl[iAtt])} Hz; a few cells give strong isolation near the resonance.`),
      kpi('T_min', 'Minimum transmissibility of the finite chain', Math.max(dB(Tr[iMin]), -300), 'dB', `Displacement of the last of ${Nc} cells relative to the driven base inside the gap, from the direct dynamic-stiffness solution; the infinite-lattice estimate −8.686·N·Im(qa) is ${f3(pred)} dB.`, Tr[iMin] < 1),
      kpi('c_lw', 'Long-wave speed', cL, 'm/s', 'Low-frequency phase velocity a√(k₁/(m₁ + m₂)): well below the gap the lattice behaves as an ordinary rod with the total static mass.'),
      kpi('f_opt', 'Upper cut-off frequency', wOpt / (2 * PI), 'Hz', 'Top of the optical branch; above it the discrete lattice passes no waves (Bragg-type cut-off).'),
    ],
    series: [
      ser('meff', 'Effective dynamic mass', 'frequency (Hz)', 'm_eff / (m₁ + m₂) (–)', [{ name: 'Re m_eff', x: fl, y: Me.map(m => N.clamp(m.re / Mt, -10, 10)) }, { name: 'Im m_eff', x: fl, y: Me.map(m => N.clamp(m.im / Mt, -10, 10)) }, { name: 'zero', x: [fl[0], fl[399]], y: [0, 0] }], 'Effective mass of one cell normalised by its static mass (clipped at ±10). It equals 1 at low frequency, diverges at the local resonance and is negative up to f₀√(1 + m₂/m₁).'),
      ser('disp', 'Bloch dispersion', 'q·a/π (–)', 'frequency (Hz)', [{ name: 'Re q (propagating)', x: Q.map(q => q.re / PI), y: fl }, { name: 'Im q (attenuating)', x: Q.map(q => q.im / PI), y: fl }, { name: 'Re q, lossless', x: Q0.map(q => q.re / PI), y: fl }], 'Real and imaginary Bloch wavenumber against frequency. The acoustic and optical branches are separated by the locally resonant gap, where Re q is pinned and Im q peaks at the resonance.'),
      ser('T', 'Transmissibility of the finite chain', 'frequency (Hz)', 'transmissibility (dB)', [{ name: `${Nc}-cell chain (dynamic stiffness)`, x: fl, y: Tr.map(v => Math.max(dB(v), -200)) }, { name: 'infinite-lattice estimate −8.686·N·Im(qa)', x: fl, y: Q.map(q => Math.max(-8.686 * Nc * q.im, -200)) }], 'Displacement transmissibility of the finite chain compared with the attenuation predicted from the infinite lattice. Pass bands show the N resonances of the finite chain; the trough coincides with the band gap.'),
    ],
    equations: [
      eqn('governing', 'Cell equations of motion', 'm₁ü₁ = k₁(u_(j+1) − 2u_j + u_(j−1)) + k₂(u₂ − u₁);  m₂ü₂ = −k₂(u₂ − u₁)'), eqn('constitutive', 'Effective dynamic mass', 'm_eff(ω) = m₁ + m₂ω₀²/(ω₀² − ω²),  ω₀² = k₂/m₂'),
      eqn('governing', 'Bloch dispersion', 'cos(qa) = 1 − m_eff(ω)ω²/(2k₁)'), eqn('formula', 'Negative-mass band and gap', 'ω₀ < ω < ω₀√(1 + m₂/m₁);  gap from the upper root of the acoustic branch to ω₀√(1 + m₂/m₁)'),
      eqn('governing', 'Finite chain (dynamic stiffness)', '[k₁(u_(j−1) + u_(j+1)) − (2k₁ − m_eff ω²)u_j] = 0,  u₀ = 1,  free end', 'structural damping k → k(1 + iη), e^(+iωt) convention'),
    ],
    checks: [
      chk('Static limit m_eff(0) = m₁ + m₂', Math.abs(meff(1e-6 * w0, 0).re - Mt) < 1e-6 * Mt, `m_eff(0) = ${f4(meff(1e-6 * w0, 0).re)} kg`), chk('Passivity: Im m_eff ≤ 0 (e^(+iωt), dissipative resonator)', Me.every(m => m.im <= 1e-12 * Mt), `min Im m_eff = ${f3(N.minOf(Me.map(m => m.im)))} kg`),
      chk('Negative effective mass present', negMin < 0, `min Re m_eff = ${f3(negMin)} kg`), chk('Finite chain attenuates inside the predicted gap', iGap.length > 0 && Tr[iMin] < 1, `minimum ${f3(dB(Tr[iMin]))} dB at ${f3(fl[iMin])} Hz; infinite-lattice estimate ${f3(pred)} dB (differences come from the free end and the drive point)`),
      chk('Lossless gap edges are band edges (|cos qa| = 1)', Math.abs(Math.abs(1 - meff(wLo, 0).re * wLo * wLo / (2 * k1)) - 1) < 1e-6, 'acoustic-branch top satisfies cos(qa) = −1'), chk('Low-frequency transmissibility tends to 1', Tr[0] > 0.9 && Tr[0] < 1.5, `|u_N/u₀| = ${f4(Tr[0])} at ${f3(fl[0])} Hz`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `To move the gap, tune the resonator: the lower part of the gap follows f₀ = ${f3(f0)} Hz and its width follows the mass ratio (now m₂/m₁ = ${f3(m2 / m1)}). Doubling the ratio widens the negative-mass band from ${f3(100 * (Math.sqrt(1 + m2 / m1) - 1))} % to ${f3(100 * (Math.sqrt(1 + 2 * m2 / m1) - 1))} % of f₀.`, 'ω₀√(1 + m₂/m₁) upper edge'),
      rec('medium', `Use moderate damping: η = ${eta} gives ${f3(attMax)} dB per cell at the peak; higher loss lowers the peak but broadens the useful attenuation band and suppresses the pass-band resonances of the finite chain.`, 'Trade-off between peak attenuation and bandwidth'),
      rec('medium', `About ${Math.max(2, Math.ceil(40 / Math.max(attMax, 1e-6)))} cells give 40 dB at the most attenuated frequency; fewer are needed than for a Bragg crystal, whose spacing would have to be half a wavelength (${f3(cL / (2 * f0))} m here).`, 'Sub-wavelength nature of locally resonant gaps'),
      rec('low', 'A negative effective mass exists only in steady harmonic response of a passive cell with hidden internal motion; it does not violate Newton\'s law and cannot be used below or above its narrow band.', 'Interpretation and validity'),
      rec('low', `Locally resonant isolators achieve low-frequency attenuation with a total added mass of ${f3(100 * m2 / Mt)} % of the cell instead of massive barriers, saving material; avoid lead cores in favour of steel or tungsten-polymer composites.`, 'Mass efficiency and hazardous-substance avoidance', true),
    ],
    publish: { 'res.f0': f0, 'res.bandgap_lo': wLo / (2 * PI), 'res.bandgap_hi': wHi / (2 * PI), 'eff.rho_dyn': rhoDyn, 'res.transmission': Math.min(Tr[iMin] ** 2, 1), 'geom.cell_size': a },
    numerics: { method: 'Closed-form Bloch dispersion with complex effective mass; complex tridiagonal (Thomas) dynamic-stiffness solve for the finite chain', dof: Nc, steps: ws.length, notes: 'Linear lumped lattice, hysteretic damping, inner mass condensed exactly.' },
  };
}

// ================================================================ negative bulk modulus / density: Helmholtz resonators and membranes in a duct
const FLUID = { air: { label: 'Air (20 °C)', rho: 1.204, c: 343, mu: 1.81e-5 }, water: { label: 'Water (20 °C)', rho: 998, c: 1482, mu: 1.0e-3 } };
function helmholtz(inp) {
  const W = [], fl = FLUID[inp.fluid] || FLUID.air, rho0 = fl.rho, c0 = fl.c, K0 = rho0 * c0 * c0, Sd = PI * (inp.duct_D * 0.5e-3) ** 2, rn = inp.neck_r * 1e-3, ln = inp.neck_l * 1e-3, Vc = inp.cav_V * 1e-6, d = inp.d_cell * 1e-3, Nc = Math.round(inp.Ncell);
  const Sn = PI * rn * rn, leff = ln + 1.7 * rn, wH = c0 * Math.sqrt(Sn / (Vc * leff)), fH = wH / (2 * PI), F = Vc / (Sd * d), Gam = Math.sqrt(2 * fl.mu * wH / rho0) * (ln + 2 * rn) / (rn * leff), hasK = inp.dn_type !== 'membrane', hasR = inp.dn_type !== 'hr';
  const rhoP = rho0 + inp.memb_sig / d, wc = Math.sqrt(8 * PI * inp.memb_T / (rhoP * Sd * d)), gm = Math.max(inp.eta, 1e-4) * wc;
  const invK = w => (hasK ? C.scale(C.sub(ONE, C.div(C.of(F * wH * wH), { re: w * w - wH * wH, im: Gam * w })), 1 / K0) : { re: 1 / K0, im: 0 });
  const rho = w => (hasR ? C.scale(C.sub(ONE, C.div(C.of(wc * wc), C.mul(C.of(w), { re: w, im: gm }))), rhoP) : { re: rho0, im: 0 });
  const prop = w => { const sr = csqrt(rho(w)), sk = csqrt(invK(w)), k = C.scale(C.mul(sr, sk), w), z = C.div(sr, sk), Kc = C.inv(invK(w)); return { k, z, K: Kc, rho: rho(w), vp: Math.abs(k.re) > 1e-12 ? w / k.re : 0, evan: k.im > Math.abs(k.re) }; };
  const wref = hasK ? wH : wc, wmax = 2.2 * Math.max(hasK ? wH * Math.sqrt(1 + F) : 0, hasR ? wc : 0, wref), ws = N.linspace(wmax / 400, wmax, 400), fr = ws.map(w => w / (2 * PI)), P = ws.map(prop);
  const negK = bands(ws, P.map(p => p.K.re < 0)), negR = bands(ws, P.map(p => p.rho.re < 0)), dng = bands(ws, P.map(p => p.K.re < 0 && p.rho.re < 0)), big = b => (b.length ? b.reduce((x, y) => (y[1] - y[0] > x[1] - x[0] ? y : x)) : null), bK = big(negK), bR = big(negR), bD = big(dng);
  const wpr = inp.f_rel * wref, pp = prop(wpr), Tn = ws.map((w, i) => C.abs(slabT(C.scale(P[i].k, Nc * d), C.scale(P[i].z, 1 / (rho0 * c0))).t)), tp = C.abs(slabT(C.scale(pp.k, Nc * d), C.scale(pp.z, 1 / (rho0 * c0))).t), lam = 2 * PI * c0 / wref;
  if (lam < 6 * d) W.push(`The free wavelength at resonance (${f3(lam * 1e3)} mm) is less than six periods: the homogenised description is unreliable and Bragg scattering interferes.`);
  if (rn > 0.5 * inp.duct_D * 0.5e-3 * 2) W.push('The neck is as wide as the duct; the lumped Helmholtz model assumes a small side branch.');
  if (hasK && !bK) W.push('Viscous loss in the neck is too high for Re K to become negative; widen the neck or increase the cavity fraction.');
  if (hasK && hasR && !bD) W.push('The negative-modulus and negative-density bands do not overlap: raise the membrane tension (higher cut-off) or lower the Helmholtz frequency to obtain a double-negative band.');
  const kk = hasK ? (Gam / wH < 3e-4 ? chk('Causality (Kramers–Kronig) of the Lorentz compressibility', true, 'resonance too sharp for the numerical quadrature; poles lie in the lower half-plane because Γ > 0') : (() => { const mm = kkMismatch(w => C.scale(C.div(C.of(F * wH * wH), { re: w * w - wH * wH, im: Gam * w }), -1), 1.3 * wH, wH); return chk('Causality (Kramers–Kronig) of the Lorentz compressibility', mm < 0.03, `Hilbert transform of Im(K₀/K) reproduces Re(K₀/K) − 1 within ${f3(100 * mm)} %`); })()) : chk('Causality of the Drude-type density', gm > 0, 'pole at ω = −iγ in the lower half-plane');
  const fb = b => (b ? `${f3(ws[b[0]] / (2 * PI))}–${f3(ws[b[1]] / (2 * PI))} Hz` : 'none');
  return {
    summary: `${hasK ? `Helmholtz resonators (neck ${inp.neck_r} mm × ${inp.neck_l} mm, cavity ${inp.cav_V} cm³) every ${inp.d_cell} mm resonate at ${f3(fH)} Hz and make the effective bulk modulus negative over ${fb(bK)}. ` : ''}${hasR ? `Membranes (tension ${inp.memb_T} N/m) make the effective density negative below ${f3(wc / (2 * PI))} Hz (${fb(bR)}). ` : ''}${hasK && hasR ? `Double-negative band: ${fb(bD)}. ` : ''}At ${f3(wpr / (2 * PI))} Hz: K_eff = ${f3(pp.K.re / K0)} K₀, ρ_eff = ${f3(pp.rho.re / rho0)} ρ₀, phase velocity ${f3(pp.vp)} m/s (${pp.evan ? 'evanescent: stop band' : pp.vp < 0 ? 'backward wave' : 'forward wave'}); ${Nc} cells transmit ${f3(dB(tp))} dB.`,
    kpis: [
      kpi('fH', 'Helmholtz resonance frequency', fH, 'Hz', 'c₀/2π·√(S_n/(V·l_eff)) with the neck length increased by 1.7 r for the two flanged ends; the bulk-modulus resonance.'),
      kpi('F', 'Cavity volume fraction', F, '–', 'Cavity volume per period divided by duct volume per period; the strength of the modulus resonance, giving K < 0 from f_H to f_H√(1 + F).'),
      kpi('negK_lo', 'Negative-modulus band lower edge', bK ? ws[bK[0]] / (2 * PI) : 0, 'Hz', 'Lowest frequency in the sweep with Re K_eff < 0 including viscous loss (zero if absent or not modelled).', !hasK || !!bK),
      kpi('negK_hi', 'Negative-modulus band upper edge', bK ? ws[bK[1]] / (2 * PI) : 0, 'Hz', `Highest frequency with Re K_eff < 0; the lossless value is f_H√(1 + F) = ${f3(fH * Math.sqrt(1 + F))} Hz.`, !hasK || !!bK),
      kpi('fc', 'Membrane cut-off frequency', hasR ? wc / (2 * PI) : 0, 'Hz', 'Below this frequency the membrane-loaded duct has negative effective density, ω_c² = 8πT/(ρ′S·d) (zero if no membranes).'),
      kpi('dng_lo', 'Double-negative band lower edge', bD ? ws[bD[0]] / (2 * PI) : 0, 'Hz', 'Lowest frequency where both Re K_eff and Re ρ_eff are negative: sound propagates again there, with negative phase velocity.', !(hasK && hasR) || !!bD),
      kpi('dng_hi', 'Double-negative band upper edge', bD ? ws[bD[1]] / (2 * PI) : 0, 'Hz', 'Highest frequency where both effective parameters are negative.', !(hasK && hasR) || !!bD),
      kpi('K_rel', 'Re K_eff/K₀ at probe', pp.K.re / K0, '–', `Effective bulk modulus relative to the fluid at ${f3(inp.f_rel)} × the reference resonance; negative means the medium expands when compressed dynamically.`, null),
      kpi('rho_rel', 'Re ρ_eff/ρ₀ at probe', pp.rho.re / rho0, '–', 'Effective density relative to the fluid; negative means acceleration opposes the pressure gradient.', null),
      kpi('vp', 'Phase velocity at probe', pp.vp, 'm/s', 'ω/Re k on the passive branch (Im k ≥ 0): negative in a double-negative band (phase travels towards the source while energy travels away), zero if evanescent.', null),
      kpi('att', 'Attenuation at probe', 8.686 * pp.k.im * d, 'dB/cell', 'Amplitude decay per period from Im k; large in single-negative bands (stop bands) and small in double-negative or double-positive bands.'),
      kpi('TL', 'Transmission loss of the array at probe', -dB(tp), 'dB', `Transmission loss of ${Nc} periods treated as a homogeneous slab in the plain duct.`),
      kpi('Q', 'Resonator quality factor', wH / Gam, '–', 'ω_H/Γ from thermoviscous boundary-layer loss in the neck; it sets how negative the modulus can become.', wH / Gam > 10),
    ],
    series: [
      ser('par', 'Effective parameters', 'frequency (Hz)', 'relative value (–)', [{ name: 'Re K_eff/K₀', x: fr, y: P.map(p => N.clamp(p.K.re / K0, -5, 5)) }, { name: 'Re ρ_eff/ρ₀', x: fr, y: P.map(p => N.clamp(p.rho.re / rho0, -5, 5)) }, { name: 'zero', x: [fr[0], fr[399]], y: [0, 0] }], 'Real effective bulk modulus and density relative to the fluid (clipped at ±5). Where only one is negative the wave is evanescent; where both are negative it propagates backwards.'),
      ser('k', 'Dispersion and attenuation', 'frequency (Hz)', 'k·d/π (–)', [{ name: 'Re k', x: fr, y: P.map(p => N.clamp(p.k.re * d / PI, -3, 3)) }, { name: 'Im k', x: fr, y: P.map(p => Math.min(p.k.im * d / PI, 3)) }], 'Real and imaginary wavenumber in units of π/d. Negative Re k marks the double-negative (backward-wave) band; large Im k marks single-negative stop bands. Values of |k·d/π| approaching 1 are outside the homogenisation limit.'),
      ser('T', 'Transmission through the array', 'frequency (Hz)', 'transmission (dB)', [{ name: `${Nc} periods`, x: fr, y: Tn.map(v => Math.max(dB(v), -200)) }], 'Pressure transmission of the finite array in the duct. Deep notches are single-negative bands; a restored pass band inside them is the double-negative region.'),
    ],
    equations: [
      eqn('formula', 'Helmholtz resonance', 'ω_H = c₀√(S_n/(V·l_eff)),  l_eff = l_n + 1.7r_n'), eqn('constitutive', 'Effective bulk modulus', '1/K_eff = (1/K₀)[1 − Fω_H²/(ω² − ω_H² + iΓω)],  F = V/(S_d·d)', 'e^(−iωt) convention'),
      eqn('constitutive', 'Effective density with membranes', 'ρ_eff = ρ′[1 − ω_c²/(ω(ω + iγ))],  ρ′ = ρ₀ + σ_m/d,  ω_c² = 8πT/(ρ′S_d·d)'), eqn('governing', 'Wavenumber and impedance (passive branch)', 'k = ω√ρ_eff·√(1/K_eff) with Im k ≥ 0;  Z = √ρ_eff/√(1/K_eff) with Re Z ≥ 0'),
      eqn('formula', 'Neck loss', 'Γ = √(2µω_H/ρ₀)·(l_n + 2r_n)/(r_n·l_eff)'),
    ],
    checks: [
      chk('Passivity: Im(1/K_eff) ≥ 0 and Im ρ_eff ≥ 0', ws.every(w => invK(w).im >= -1e-18 && rho(w).im >= -1e-12), 'both effective parameters dissipate at every frequency (e^(−iωt))'), chk('Passive branch: Im k ≥ 0, Re Z ≥ 0', P.every(p => p.k.im >= -1e-12 && p.z.re >= -1e-9), 'waves decay in the direction of energy flow'),
      chk('Transmission bounded |t| ≤ 1', Tn.every(v => v <= 1 + 1e-9), `max |t| = ${f4(N.maxOf(Tn))}`), kk,
      chk('Sub-wavelength period (λ > 6d)', lam >= 6 * d, `λ/d = ${f3(lam / d)} at the reference resonance`), chk('Negative band present', hasK ? !!bK : !!bR, hasK ? `K < 0 over ${fb(bK)}` : `ρ < 0 over ${fb(bR)}`),
      chk('Backward wave in the double-negative band', !(hasK && hasR && bD) || P[Math.round((bD[0] + bD[1]) / 2)].vp < 0, bD ? `phase velocity ${f3(P[Math.round((bD[0] + bD[1]) / 2)].vp)} m/s at mid-band` : 'no double-negative band to test'),
    ],
    warnings: W,
    recommendations: [
      rec('high', hasK ? `Tune the stop band with the cavity: f_H = ${f3(fH)} Hz scales as 1/√V and the bandwidth as √(1 + F) (F = ${f3(F)}). ${hasR && !bD ? `For a double-negative band the membrane cut-off (${f3(wc / (2 * PI))} Hz) must exceed f_H.` : ''}` : `The density is negative below ${f3(wc / (2 * PI))} Hz; raise the membrane tension or shorten the period to extend the band.`, 'Lorentz resonance of the compressibility / Drude form of the density'),
      rec('medium', `Neck loss gives Q = ${f3(wH / Gam)}: narrow necks (below about ${f3(10 * Math.sqrt(2 * fl.mu / (rho0 * wH)) * 1e3)} mm, ten viscous boundary layers) are heavily damped and suit absorbers; wide necks suit negative-modulus transmission lines.`, 'Viscous boundary-layer thickness at resonance'),
      rec('medium', 'Verify by a four-microphone impedance-tube measurement and retrieve K_eff and ρ_eff from the transfer matrix; evanescent coupling between neighbouring resonators shifts the band at small spacing.', 'Lumped-element homogenisation'),
      rec('low', 'Resonator liners attenuate low frequencies in a fraction of the thickness of porous absorbers and contain no fibres; they can be moulded from a single recyclable polymer.', 'Material use and fibre-free construction', true),
    ],
    publish: { 'res.f0': hasK ? fH : wc / (2 * PI), 'eff.K': pp.K.re, 'eff.rho_dyn': pp.rho.re, 'res.transmission': Math.min(tp * tp, 1), ...(bK ? { 'res.bandgap_lo': ws[bK[0]] / (2 * PI), 'res.bandgap_hi': ws[bK[1]] / (2 * PI) } : bR ? { 'res.bandgap_lo': 0, 'res.bandgap_hi': ws[bR[1]] / (2 * PI) } : {}) },
    numerics: { method: 'Lumped-element homogenisation (Lorentz compressibility, Drude-type density) with passive-branch square roots and Fresnel–Airy slab transmission', dof: ws.length, notes: 'Plane waves in the duct below its cut-on frequency; resonator interaction and radiation loss neglected.' },
  };
}

// ================================================================ 1-D phononic crystal (layered)
const SOLID = { // density kg/m³, longitudinal speed m/s (handbook-typical)
  steel: { label: 'Steel', rho: 7850, c: 5900, f: 'Fe' }, Al: { label: 'Aluminium', rho: 2700, c: 6320, f: 'Al' }, epoxy: { label: 'Epoxy', rho: 1180, c: 2540 }, PMMA: { label: 'PMMA', rho: 1190, c: 2700 }, W: { label: 'Tungsten', rho: 19300, c: 5200, f: 'W' },
  Si: { label: 'Silicon [100]', rho: 2329, c: 8433, f: 'Si' }, SiO2: { label: 'Fused silica', rho: 2200, c: 5970, f: 'SiO2' }, ZnO: { label: 'Zinc oxide (c-axis)', rho: 5606, c: 6133, f: 'ZnO' }, GaAs: { label: 'GaAs [100]', rho: 5317, c: 4730, f: 'GaAs' }, AlAs: { label: 'AlAs [100]', rho: 3760, c: 5650, f: 'AlAs' },
  water: { label: 'Water', rho: 998, c: 1482, f: 'H2O' }, Pb: { label: 'Lead', rho: 11340, c: 2160, f: 'Pb' },
};
function phononic(inp) {
  const W = [], A = SOLID[inp.matA] || SOLID.steel, B = SOLID[inp.matB] || SOLID.epoxy, dA = inp.dA * 1e-6, dB2 = inp.dB * 1e-6, Lam = dA + dB2, Nc = Math.round(inp.Ncell), eta = inp.eta, ZA = A.rho * A.c, ZB = B.rho * B.c, fB = 1 / (2 * (dA / A.c + dB2 / B.c));
  const lay = (w, m, d, e) => { const kd = C.scale({ re: 1, im: -e / 2 }, w * d / m.c), zr = C.scale({ re: 1, im: e / 2 }, m.rho * m.c / ZA), c = C.cos(kd), s = C.sin(kd); return [c, C.div(s, zr), C.scale(C.mul(zr, s), -1), c]; };
  const mul2 = (a, b) => [C.add(C.mul(a[0], b[0]), C.mul(a[1], b[2])), C.add(C.mul(a[0], b[1]), C.mul(a[1], b[3])), C.add(C.mul(a[2], b[0]), C.mul(a[3], b[2])), C.add(C.mul(a[2], b[1]), C.mul(a[3], b[3]))];
  const per = (w, e) => mul2(lay(w, B, dB2, e), lay(w, A, dA, e)), tr = w => { const m = per(w, 0); return 0.5 * (m[0].re + m[3].re); }, rytov = w => Math.cos(w * dA / A.c) * Math.cos(w * dB2 / B.c) - 0.5 * (ZA / ZB + ZB / ZA) * Math.sin(w * dA / A.c) * Math.sin(w * dB2 / B.c);
  const fmax = 3.3 * fB, fs = N.linspace(fmax / 400, fmax, 400), ws = fs.map(f => 2 * PI * f), trs = ws.map(tr), gaps = bands(fs, trs.map(v => Math.abs(v) > 1)).slice(0, 3).map(([i, j]) => { const e = (a, b) => (a < 0 || b < 0 || a >= 400 || b >= 400 ? fs[N.clamp(a, 0, 399)] : N.brent(f => Math.abs(tr(2 * PI * f)) - 1, fs[a], fs[b], fs[1] * 1e-9)); try { return [e(i - 1, i), e(j, j + 1)]; } catch { return [fs[i], fs[j]]; } });
  const Kb = ws.map(w => { const m = per(w, eta); return cacos(C.scale(C.add(m[0], m[3]), 0.5)); });
  const trans = w => { let m = [ONE, { re: 0, im: 0 }, { re: 0, im: 0 }, ONE]; const p = per(w, eta); for (let i = 0; i < Nc; i++) { m = mul2(p, m); const g = Math.max(C.abs(m[0]), C.abs(m[1]), C.abs(m[2]), C.abs(m[3])); if (g > 1e100) return 0; }
    const q = { re: 0, im: -1 }, q2 = C.mul(q, q), den = C.add(C.sub(C.sub(C.mul(q, m[0]), C.mul(q2, m[1])), m[2]), C.mul(q, m[3])), num = C.sub(C.sub(C.add(m[2], C.mul(q, m[3])), C.mul(q, m[0])), C.mul(q2, m[1])), r = C.div(num, den), t = C.add(C.mul(m[0], C.add(ONE, r)), C.mul(C.mul(m[1], q), C.sub(ONE, r))); return Math.min(C.abs(t), 1); };
  const T = ws.map(trans), g1 = gaps[0] || null, fc = g1 ? 0.5 * (g1[0] + g1[1]) : fB, ic = N.argmin(fs.map(f => Math.abs(f - fc))), attC = 8.686 * Kb[ic].im, qw = Math.abs(dA / A.c - dB2 / B.c) < 0.02 * (dA / A.c), bwAn = (4 / PI) * Math.asin(Math.abs(ZA - ZB) / (ZA + ZB)), Tmin = g1 ? N.minOf(T.filter((_, i) => fs[i] >= g1[0] && fs[i] <= g1[1])) : N.minOf(T);
  const Rint = ((ZA - ZB) / (ZA + ZB)) ** 2, ryErr = N.maxOf(ws.map((w, i) => Math.abs(rytov(w) - trs[i])));
  if (!g1) W.push('No band gap in the scanned range: the impedance contrast between the two layers is too small.');
  if (Math.abs(ZA - ZB) / (ZA + ZB) < 0.05) W.push('Impedance contrast below 5 %: gaps are very narrow and need many periods to become visible.');
  if (inp.matA === 'water' || inp.matB === 'water') W.push('A fluid layer carries only longitudinal waves; the model treats normal incidence, where this is exact, but oblique incidence would couple to shear waves in the solid.');
  return {
    summary: `${A.label}/${B.label} superlattice (${f3(inp.dA)}/${f3(inp.dB)} µm, period ${f3(Lam * 1e6)} µm): acoustic impedances ${f3(ZA / 1e6)} and ${f3(ZB / 1e6)} MRayl. ${g1 ? `The first Bragg gap spans ${f3(g1[0])}–${f3(g1[1])} Hz (centre ${f3(fc)} Hz, ${f3(100 * (g1[1] - g1[0]) / fc)} % wide)` : 'No gap was found'}; the Bragg estimate is ${f3(fB)} Hz. Mid-gap attenuation is ${f3(attC)} dB per period and ${Nc} periods transmit down to ${f3(Math.max(dB(Tmin), -300))} dB in the gap.`,
    kpis: [
      kpi('fB', 'Bragg frequency estimate', fB, 'Hz', 'Frequency at which the period is half an effective wavelength, 1/(2(d_A/c_A + d_B/c_B)); centre of the first gap.'),
      kpi('gap_lo', 'First gap lower edge', g1 ? g1[0] : 0, 'Hz', 'Lower band edge from the Bloch condition |½ tr M| = 1 of the lossless period transfer matrix (zero if no gap).', !!g1),
      kpi('gap_hi', 'First gap upper edge', g1 ? g1[1] : 0, 'Hz', 'Upper band edge of the first Bragg gap, refined by root finding on the Bloch condition (zero if no gap).', !!g1),
      kpi('gap_w', 'Relative gap width', g1 ? (g1[1] - g1[0]) / fc : 0, '–', `Gap width over centre frequency; a quarter-wave stack gives (4/π)·asin(|Z_A − Z_B|/(Z_A + Z_B)) = ${f3(bwAn)}.`),
      kpi('n_gaps', 'Gaps found below 3.3 f_B', gaps.length, '–', 'Number of stop bands in the scanned range (up to three are listed); even-order gaps close for quarter-wave stacks.'),
      kpi('att', 'Mid-gap attenuation', attC, 'dB/period', 'Imaginary Bloch wavenumber at the gap centre, 8.686·Im(qΛ): each period reduces the amplitude by this amount.'),
      kpi('Tmin', 'Minimum transmission in the gap', Math.max(dB(Tmin), -300), 'dB', `Amplitude transmission of ${Nc} periods between half-spaces of material A (values below −300 dB are reported as −300 dB).`, Tmin < 0.5),
      kpi('ZA', 'Impedance of layer A', ZA, 'Pa·s/m', 'Characteristic acoustic impedance ρc of layer A for longitudinal waves.'), kpi('ZB', 'Impedance of layer B', ZB, 'Pa·s/m', 'Characteristic acoustic impedance ρc of layer B for longitudinal waves.'),
      kpi('R', 'Single-interface power reflection', Rint, '–', '((Z_A − Z_B)/(Z_A + Z_B))²: the reflection at each interface whose coherent addition opens the gap.'),
      kpi('c_eff', 'Low-frequency effective speed', Lam / (dA / A.c + dB2 / B.c) * Math.sqrt(1 / (1 + (dA / A.c) * (dB2 / B.c) / (dA / A.c + dB2 / B.c) ** 2 * (ZA / ZB + ZB / ZA - 2))), 'm/s', 'Long-wavelength phase velocity of the layered medium, Λ/√(⟨ρ⟩⟨1/M⟩), below the time-average speed because of impedance mismatch.'),
    ],
    series: [
      ser('band', 'Bloch band structure', 'q·Λ/π (–)', 'frequency (Hz)', [{ name: 'Re q', x: Kb.map(k => k.re / PI), y: fs }, { name: 'Im q', x: Kb.map(k => k.im / PI), y: fs }], 'Folded dispersion of the infinite superlattice. Gaps open at the zone centre and edge, where Re q is pinned and Im q is non-zero.'),
      ser('T', 'Transmission of the finite crystal', 'frequency (Hz)', 'transmission (dB)', [{ name: `${Nc} periods`, x: fs, y: T.map(v => Math.max(dB(v), -200)) }], 'Amplitude transmission through the finite stack. Stop bands deepen linearly in dB with the number of periods; ripples in the pass bands are Fabry–Perot resonances of the whole stack.'),
    ],
    equations: [...(A.f ? formulaEq(A.f, `Layer A: ${A.label}`) : []), ...(B.f ? formulaEq(B.f, `Layer B: ${B.label}`) : []),
      eqn('governing', 'Elastodynamic wave equation (1-D)', 'ρ ∂²u/∂t² = ∂/∂x(ρc² ∂u/∂x)'), eqn('governing', 'Layer transfer matrix', '[u; σ]ₓ₊d = [[cos kd, sin kd/(ωZ)], [−ωZ sin kd, cos kd]]·[u; σ]ₓ,  Z = ρc'),
      eqn('governing', 'Bloch (Rytov) dispersion', 'cos(qΛ) = cos(k_A d_A) cos(k_B d_B) − ½(Z_A/Z_B + Z_B/Z_A) sin(k_A d_A) sin(k_B d_B)'), eqn('formula', 'Bragg condition', 'f_B = 1/[2(d_A/c_A + d_B/c_B)]'),
    ],
    checks: [
      chk('Transfer-matrix trace equals the Rytov formula', ryErr < 1e-9, `max difference ${f3(ryErr)}`), chk('Bragg frequency lies inside the first gap', !g1 || (fB > g1[0] && fB < g1[1]), g1 ? `f_B = ${f3(fB)} Hz within ${f3(g1[0])}–${f3(g1[1])} Hz (the gap is asymmetric about f_B for unequal layer transit times)` : 'no gap'),
      chk('Quarter-wave gap width reproduced', !qw || !g1 || Math.abs((g1[1] - g1[0]) / fc - bwAn) < 0.02 * bwAn + 1e-3, qw ? `numerical ${f4(g1 ? (g1[1] - g1[0]) / fc : 0)}, analytical ${f4(bwAn)}` : 'layers are not quarter-wave matched; check skipped'),
      chk('Passivity: |t| ≤ 1', T.every(v => v <= 1 + 1e-9), `max |t| = ${f4(N.maxOf(T))}`), chk('Band gap present', !!g1, g1 ? `${gaps.length} gap(s)` : 'none'),
    ],
    warnings: W,
    recommendations: [
      rec('high', g1 ? `For a gap centred at a target frequency f, scale both thicknesses by ${f3(fc)}/f; quarter-wave layers (d = c/4f: ${f3(A.c / (4 * fc) * 1e6)} and ${f3(B.c / (4 * fc) * 1e6)} µm here) give the widest first gap for this material pair.` : 'Choose materials with a larger impedance ratio (for example a dense metal against a polymer) to open a gap.', 'Bragg condition and quarter-wave optimum'),
      rec('medium', `${Math.max(1, Math.ceil(40 / Math.max(attC, 1e-6)))} periods are enough for 40 dB mid-gap rejection; beyond that, material damping (η = ${eta}) rather than reflection limits performance.`, 'Attenuation per period at mid-gap'),
      rec('medium', 'Thickness errors shift the gap in proportion: hold layer thickness to a few per cent, and at gigahertz–terahertz frequencies (nanometre layers) account for interface roughness, which scatters phonons diffusely.', 'Sensitivity of the Bragg frequency to period'),
      rec('low', 'Phononic mirrors confine vibration without added damping material or active control, and thin-film superlattices reduce thermal conductivity for thermoelectrics without scarce elements when oxide or silicon-based pairs are chosen.', 'Passive function and material choice', true),
    ],
    publish: { 'res.f0': fc, 'res.transmission': Math.min(Tmin * Tmin, 1), ...(g1 ? { 'res.bandgap_lo': g1[0], 'res.bandgap_hi': g1[1] } : {}), 'geom.cell_size': Lam },
    numerics: { method: 'Complex transfer-matrix method for longitudinal waves at normal incidence; Brent refinement of band edges', dof: 2, steps: fs.length, notes: 'Perfectly bonded isotropic layers, hysteretic damping, embedding medium equal to material A.' },
  };
}

// ================================================================ sound-insulating panels: mass law against locally resonant plate
const AIR = FLUID.air;
function panel(inp) {
  const W = [], rc = AIR.rho * AIR.c, h = inp.h_p * 1e-3, ms = inp.rho_p * h, D = inp.E_p * h ** 3 / (12 * (1 - inp.nu_p ** 2)), fcz = AIR.c ** 2 / (2 * PI) * Math.sqrt(ms / D), fr = inp.f_res, wr = 2 * PI * fr, mu = inp.mu_r, eta = Math.max(inp.eta, 1e-4);
  const meff = w => { const kc = { re: 1, im: eta }; return C.scale(C.add(ONE, C.scale(C.div(kc, C.sub(kc, C.of((w / wr) ** 2))), mu)), ms); };
  const TLm = (f, m = ms) => 10 * Math.log10(1 + (PI * f * m / rc) ** 2), TLr = f => { const w = 2 * PI * f, me = meff(w), zr = C.scale(C.mul(J, me), w / (2 * rc)); return 10 * Math.log10(C.abs2(C.add(ONE, zr))); };
  const fs = N.logspace(Math.max(fr / 20, 5), Math.max(fr * 20, 200), 400), tl = fs.map(TLr), tm = fs.map(f => TLm(f)), tmT = fs.map(f => TLm(f, ms * (1 + mu))), ipk = N.argmax(tl.map((v, i) => v - tmT[i])), fdip = fr * Math.sqrt(1 + mu), gain = tl[ipk] - tmT[ipk];
  const gb = bands(fs, tl.map((v, i) => v - tmT[i] > 3)), gw = gb.length ? gb.reduce((x, y) => (y[1] - y[0] > x[1] - x[0] ? y : x)) : null, wpr = inp.f_rel * wr, mp = meff(wpr), rhoDyn = mp.re / h, Me = fs.map(f => meff(2 * PI * f));
  if (fr > 0.5 * fcz) W.push(`The resonator frequency is above half the coincidence frequency (${f3(fcz)} Hz): plate bending stiffness and coincidence, not included in the limp-panel model, control the transmission there.`);
  if (fdip < fs[399] && TLr(fdip) < TLm(fdip, ms * (1 + mu)) - 3) W.push(`Above the resonance the effective mass passes through zero near ${f3(fdip)} Hz and the panel becomes more transparent than the mass law: place this dip away from the noise spectrum or add damping.`);
  return {
    summary: `A ${inp.h_p} mm panel (${f3(ms)} kg/m²) with resonators adding ${f3(100 * mu)} % mass tuned to ${f3(fr)} Hz: transmission loss at resonance is ${f3(TLr(fr))} dB against ${f3(TLm(fr, ms * (1 + mu)))} dB for a plain panel of the same total mass (gain ${f3(TLr(fr) - TLm(fr, ms * (1 + mu)))} dB). The effective surface mass is negative from ${f3(fr)} to ${f3(fdip)} Hz, with a transmission dip near ${f3(fdip)} Hz. Coincidence frequency of the bare plate: ${f3(fcz)} Hz.`,
    kpis: [
      kpi('ms', 'Surface mass of the plate', ms, 'kg/m²', 'Mass per unit area ρh of the host panel; the mass law gains 6 dB per doubling of it or of frequency.'),
      kpi('TL_res', 'Transmission loss at resonance', TLr(fr), 'dB', 'Normal-incidence sound transmission loss of the resonator-loaded limp panel at the tuned frequency, limited by resonator damping.'),
      kpi('TL_mass', 'Mass-law transmission loss (equal mass)', TLm(fr, ms * (1 + mu)), 'dB', '10·log₁₀[1 + (πf·m″/ρ₀c₀)²] for a homogeneous panel with the same total mass at the same frequency; field-incidence values are about 5 dB lower.'),
      kpi('gain', 'Peak gain over the mass law', gain, 'dB', `Largest excess of the metamaterial panel over the equal-mass panel, found at ${f3(fs[ipk])} Hz.`, gain > 3),
      kpi('bw', 'Bandwidth with more than 3 dB gain', gw ? Math.log2(fs[gw[1]] / fs[gw[0]]) : 0, 'octaves', 'Width of the frequency band in which the resonant panel beats the equal-mass law by at least 3 dB; it grows with mass ratio and damping.'),
      kpi('f_dip', 'Anti-resonance (transmission dip)', fdip, 'Hz', 'Frequency f_r√(1 + µ) at which the effective mass vanishes and the panel is acoustically transparent.'),
      kpi('m_eff', 'Effective surface mass at probe', mp.re, 'kg/m²', `Re m″_eff at ${f3(inp.f_rel)} × f_r; negative between the resonance and the anti-resonance.`, mp.re < 0),
      kpi('rho_dyn', 'Effective dynamic density at probe', rhoDyn, 'kg/m³', 'Effective surface mass divided by the plate thickness: the dynamic density of the equivalent homogeneous layer.', null),
      kpi('fc', 'Coincidence frequency', fcz, 'Hz', 'c₀²/2π·√(m″/D): above it bending waves radiate efficiently and the limp-panel mass law overestimates insulation.', fr < 0.5 * fcz),
      kpi('tau', 'Transmission coefficient at resonance', 10 ** (-TLr(fr) / 10), '–', 'Fraction of incident sound power transmitted at the tuned frequency.'),
    ],
    series: [
      ser('TL', 'Sound transmission loss', 'frequency (Hz)', 'transmission loss (dB)', [{ name: 'locally resonant panel', x: fs, y: tl }, { name: 'mass law, same total mass', x: fs, y: tmT }, { name: 'mass law, bare plate', x: fs, y: tm }], 'Normal-incidence transmission loss. The resonant panel rises above the 6 dB/octave mass law around the tuned frequency and dips below it at the anti-resonance; elsewhere it follows the mass law of its static mass.', { logx: true }),
      ser('m', 'Effective surface mass', 'frequency (Hz)', 'm″_eff / m″ (–)', [{ name: 'Re', x: fs, y: Me.map(m => N.clamp(m.re / ms, -20, 20)) }, { name: 'Im', x: fs, y: Me.map(m => N.clamp(m.im / ms, -20, 20)) }], 'Effective mass per unit area relative to the bare plate (clipped at ±20): large and positive below resonance, negative just above it, zero at the anti-resonance.', { logx: true }),
    ],
    equations: [
      eqn('formula', 'Mass law (normal incidence)', 'TL = 10 log₁₀[1 + (ω m″/2ρ₀c₀)²]'), eqn('constitutive', 'Effective surface mass with resonators', 'm″_eff = m″[1 + µ ω_r²(1 + iη)/(ω_r²(1 + iη) − ω²)]'),
      eqn('governing', 'Limp-panel transmission', 'τ = |1 + iω m″_eff/(2ρ₀c₀)|⁻²'), eqn('formula', 'Coincidence frequency', 'f_c = (c₀²/2π)√(m″/D),  D = Eh³/[12(1 − ν²)]'),
    ],
    checks: [
      chk('Low-frequency limit follows the total-mass law', Math.abs(tl[0] - tmT[0]) < 0.5, `${f3(tl[0])} dB against ${f3(tmT[0])} dB at ${f3(fs[0])} Hz`), chk('Passivity: Im m″_eff ≤ 0 and TL ≥ 0', Me.every(m => m.im <= 1e-12) && tl.every(v => v >= -1e-9), 'the panel absorbs or reflects but never amplifies'),
      chk('Negative effective mass present', N.minOf(Me.map(m => m.re)) < 0, `min Re m″_eff/m″ = ${f3(N.minOf(Me.map(m => m.re)) / ms)}`), chk('Tuned below coincidence (limp-panel model valid)', fr < 0.5 * fcz, `f_r/f_c = ${f3(fr / fcz)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `Tune resonators to the dominant tone: the gain of ${f3(gain)} dB is confined to about ${f3(gw ? Math.log2(fs[gw[1]] / fs[gw[0]]) : 0)} octaves. For broadband noise distribute several tunings across the band or raise the mass ratio above ${f3(mu)}.`, 'Narrow-band nature of local resonance'),
      rec('medium', `Keep the anti-resonance dip at ${f3(fdip)} Hz away from other noise components, or damp it: a loss factor of 0.05–0.1 fills the dip at a small cost in peak insulation.`, 'm_eff = 0 at f_r√(1 + µ)'),
      rec('medium', 'Resonator spacing must be below half the bending wavelength at the tuned frequency for the smeared-mass model to hold; verify with an impedance-tube or two-room measurement.', 'Homogenisation of discrete resonators'),
      rec('medium', `The same insulation at ${f3(fr)} Hz from mass alone would need about ${f3(10 ** ((TLr(fr) - TLm(fr, ms * (1 + mu))) / 20))} times the surface mass: resonant panels save material and, in vehicles, fuel.`, 'Equal-mass comparison with the mass law', true),
    ],
    publish: { 'res.f0': fr, 'res.transmission': 10 ** (-TLr(fr) / 10), 'eff.rho_dyn': rhoDyn, 'res.bandgap_lo': fr, 'res.bandgap_hi': fdip },
    numerics: { method: 'Closed-form limp-panel transmission with complex effective surface mass', dof: fs.length, notes: 'Normal incidence, infinite panel, no bending stiffness or coincidence effects.' },
  };
}

// ================================================================ absorbers: micro-perforated panel (Maa) and quarter-wave tubes
function mpp(inp) {
  const W = [], rho0 = AIR.rho, c0 = AIR.c, mu = AIR.mu, d = inp.d_hole * 1e-3, t = inp.t_mpp * 1e-3, sg = inp.perf / 100, Dc = inp.D_cav * 1e-3, Lq = inp.L_qw * 1e-3, rq = inp.r_qw * 1e-3, phq = inp.phi_qw / 100;
  const maa = f => { const w = 2 * PI * f, k = d * Math.sqrt(w * rho0 / (4 * mu)), r = 32 * mu * t / (sg * rho0 * c0 * d * d) * (Math.sqrt(1 + k * k / 32) + Math.SQRT2 / 32 * k * d / t), wm = w * t / (sg * c0) * (1 + 1 / Math.sqrt(9 + k * k / 2) + 0.85 * d / t), x = wm - 1 / Math.tan(w * Dc / c0); return { r, x, k, a: 4 * r / ((1 + r) ** 2 + x * x) }; };
  const qw = f => { const w = 2 * PI * f, dv = Math.sqrt(2 * mu / (rho0 * w)), dh = dv / Math.sqrt(0.71), b1 = (dv + 0.4 * dh) / (2 * rq), b2 = (dv - 0.4 * dh) / (2 * rq), kc = C.scale({ re: 1 + b1, im: -b1 }, w / c0), zc = { re: 1 + b2, im: -b2 }, kl = C.scale(kc, Lq + 0.85 * rq), z = C.scale(C.mul(C.mul({ re: 0, im: -1 }, zc), C.div(C.cos(kl), C.sin(kl))), 1 / phq); return 1 - C.abs2(C.div(C.sub(z, ONE), C.add(z, ONE))); };
  const fs = N.logspace(50, 10000, 400), A = fs.map(maa), al = A.map(a => a.a), aq = fs.map(qw), ip = N.argmax(al), fpk = ip > 0 && ip < 399 ? N.goldenSection(f => -maa(f).a, fs[ip - 1], fs[ip + 1], 1e-3) : fs[ip], pk = maa(fpk);
  let i1 = ip, i2 = ip; while (i1 > 0 && al[i1 - 1] >= 0.5) i1--; while (i2 < 399 && al[i2 + 1] >= 0.5) i2++; const bw = al[ip] >= 0.5 ? Math.log2(fs[i2] / fs[i1]) : 0, fq0 = c0 / (4 * (Lq + 0.85 * rq)), iq = (() => { let b = -1; fs.forEach((f, i) => { if (f > 0.6 * fq0 && f < 1.4 * fq0 && (b < 0 || aq[i] > aq[b])) b = i; }); return b < 0 ? N.argmax(aq) : b; })(), nrc = N.mean([250, 500, 1000, 2000].map(f => maa(f).a));
  if (d > 1e-3) W.push('Hole diameter above 1 mm: the panel is an ordinary perforate with little viscous resistance; Maa\'s micro-perforated theory needs sub-millimetre holes.');
  if (pk.k > 10 || pk.k < 1) W.push(`Perforate constant k = ${f3(pk.k)} at the peak is outside the range 1–10 in which Maa's approximation is accurate and wide-band absorption is obtained.`);
  if (sg > 0.1) W.push('Perforation above 10 %: hole interaction is neglected by the model.');
  return {
    summary: `Micro-perforated panel (holes ${inp.d_hole} mm, thickness ${inp.t_mpp} mm, perforation ${inp.perf} %, cavity ${inp.D_cav} mm): peak normal-incidence absorption ${f3(pk.a)} at ${f3(fpk)} Hz, with α ≥ 0.5 over ${f3(bw)} octaves (${f3(fs[i1])}–${f3(fs[i2])} Hz) and a relative resistance of ${f3(pk.r)} ρ₀c₀. A quarter-wave tube array (${inp.L_qw} mm, open area ${inp.phi_qw} %) has its fundamental peak at ${f3(fs[iq])} Hz with α = ${f3(aq[iq])}.`,
    kpis: [
      kpi('a_pk', 'Peak absorption coefficient', pk.a, '–', 'Maximum normal-incidence absorption of the micro-perforated panel in front of its air cavity (Maa model); 1 means all incident sound energy is dissipated in the holes.', pk.a > 0.8),
      kpi('f_pk', 'Frequency of peak absorption', fpk, 'Hz', 'Resonance of the hole mass against the cavity stiffness, where the reactance ωm − cot(ωD/c₀) vanishes.'),
      kpi('bw', 'Half-absorption bandwidth', bw, 'octaves', 'Width of the band with α ≥ 0.5 around the peak; micro-perforated absorbers reach 2–3 octaves when the resistance is near ρ₀c₀ and the mass reactance is small.', bw > 1.5),
      kpi('r', 'Relative acoustic resistance', pk.r, '–', 'Specific resistance of the perforate normalised by ρ₀c₀ at the peak; values near 1 give complete absorption at resonance.', pk.r > 0.5 && pk.r < 2),
      kpi('k', 'Perforate constant', pk.k, '–', 'Ratio of hole diameter to viscous boundary-layer thickness, d·√(ωρ₀/4µ); between about 1 and 10 the hole is viscosity-dominated.', pk.k >= 1 && pk.k <= 10),
      kpi('nrc', 'Mean absorption 250–2000 Hz', nrc, '–', 'Average of the normal-incidence absorption at 250, 500, 1000 and 2000 Hz (a noise-reduction-coefficient-like figure; diffuse-field values differ).'),
      kpi('f_qw', 'Quarter-wave tube peak frequency', fs[iq], 'Hz', `Frequency of the fundamental absorption peak of the tube array; the ideal value c₀/(4L_eff) is ${f3(fq0)} Hz, lowered slightly by wall losses.`),
      kpi('a_qw', 'Quarter-wave tube peak absorption', aq[iq], '–', 'Fundamental-peak absorption of the rigid-backed tube array with thermoviscous wall loss; it is maximal when the tube resistance divided by the open-area ratio equals ρ₀c₀.'),
      kpi('thick', 'Absorber depth over wavelength at peak', (Dc + t) * fpk / c0, '–', 'Total depth of the micro-perforated absorber relative to the wavelength at its peak; a quarter-wave porous layer would need 0.25.'),
    ],
    series: [
      ser('alpha', 'Absorption coefficient', 'frequency (Hz)', 'α (–)', [{ name: 'micro-perforated panel (Maa)', x: fs, y: al }, { name: 'quarter-wave tube array', x: fs, y: aq }], 'Normal-incidence absorption against frequency. The micro-perforated panel gives a broad resonance with repeats at higher cavity modes; the tube array gives sharp peaks at odd multiples of its quarter-wave frequency.', { logx: true }),
      ser('z', 'Impedance of the micro-perforated absorber', 'frequency (Hz)', 'z/ρ₀c₀ (–)', [{ name: 'resistance', x: fs, y: A.map(a => a.r) }, { name: 'reactance', x: fs, y: A.map(a => N.clamp(a.x, -20, 20)) }, { name: 'zero', x: [50, 10000], y: [0, 0] }], 'Normalised surface resistance and reactance (clipped at ±20). Absorption peaks where the reactance crosses zero; bandwidth is wide when the reactance slope there is small compared with 1 + r.', { logx: true }),
    ],
    equations: [
      eqn('constitutive', 'Maa impedance of a micro-perforated panel', 'r = (32µt/σρ₀c₀d²)[√(1 + k²/32) + (√2/32)k·d/t],  ωm = (ωt/σc₀)[1 + 1/√(9 + k²/2) + 0.85d/t],  k = d√(ωρ₀/4µ)'),
      eqn('governing', 'Absorption with a backing cavity', 'α = 4r/[(1 + r)² + (ωm − cot(ωD/c₀))²]'), eqn('governing', 'Quarter-wave tube array', 'z = −i Z_c cot(k_c L_eff)/φ,  α = 1 − |(z − 1)/(z + 1)|²,  f₁ = c₀/(4L_eff)', 'Kirchhoff wide-tube wavenumber and impedance with viscous and thermal boundary layers'),
    ],
    checks: [
      chk('Absorption bounded 0 ≤ α ≤ 1', al.every(v => v >= 0 && v <= 1 + 1e-12) && aq.every(v => v >= -1e-9 && v <= 1 + 1e-9), 'passive absorber'), chk('Maa validity range 1 ≤ k ≤ 10 at the peak', pk.k >= 1 && pk.k <= 10, `k = ${f3(pk.k)}`),
      chk('Sub-millimetre perforation', d <= 1e-3, `d = ${inp.d_hole} mm`), chk('Quarter-wave peak near c₀/4L', Math.abs(fs[iq] - fq0) < 0.15 * fq0 || aq[iq] < 0.05, `peak ${f3(fs[iq])} Hz against ${f3(fq0)} Hz`),
      chk('Peak inside the scanned range', ip > 0 && ip < 399, `${f3(fpk)} Hz within 50–10000 Hz`),
    ],
    warnings: W.filter(Boolean),
    recommendations: [
      rec('high', pk.r > 2 ? `Resistance ${f3(pk.r)} ρ₀c₀ is too high: raise the perforation ratio or the hole diameter to bring it towards 1.` : pk.r < 0.5 ? `Resistance ${f3(pk.r)} ρ₀c₀ is too low: use smaller holes (resistance ∝ 1/d²) or a lower perforation ratio.` : `Resistance ${f3(pk.r)} ρ₀c₀ is well matched; shift the band with the cavity depth (deeper cavity lowers the peak, now ${f3(fpk)} Hz).`, 'Impedance matching r ≈ 1 at resonance'),
      rec('medium', 'For more than three octaves stack two micro-perforated layers with different cavities, or partition the cavity into cells of several depths.', 'Single-resonance bandwidth limit'),
      rec('medium', 'Hole diameter tolerance matters as 1/d²: laser- or punch-made holes with burrs or taper change the resistance; measure in an impedance tube before production.', 'Maa resistance term'),
      rec('medium', 'Micro-perforated metal, glass or polymer sheets absorb without mineral fibres or foams: they are cleanable, non-shedding and fully recyclable, which suits hygienic and outdoor use.', 'Fibre-free absorber', true),
    ],
    publish: { 'res.absorption': pk.a, 'res.f0': fpk, 'res.bandgap_lo': fs[i1], 'res.bandgap_hi': fs[i2] },
    numerics: { method: 'Closed-form Maa impedance and Kirchhoff tube model on a logarithmic frequency grid; golden-section peak refinement', dof: fs.length, notes: 'Normal incidence, locally reacting absorber, linear (low-amplitude) regime.' },
  };
}

// ================================================================ Christoffel equation
const XT = { // stiffness in GPa, density kg/m³ (room-temperature handbook values)
  Si: { label: 'Silicon (cubic)', sym: 'cubic', c11: 165.7, c12: 63.9, c44: 79.6, rho: 2329, f: 'Si' }, GaAs: { label: 'Gallium arsenide (cubic)', sym: 'cubic', c11: 118.8, c12: 53.8, c44: 59.4, rho: 5317, f: 'GaAs' },
  C: { label: 'Diamond (cubic)', sym: 'cubic', c11: 1079, c12: 124, c44: 578, rho: 3515, f: 'C' }, Cu: { label: 'Copper (cubic)', sym: 'cubic', c11: 168.4, c12: 121.4, c44: 75.4, rho: 8960, f: 'Cu' },
  ZnO: { label: 'Zinc oxide (hexagonal)', sym: 'hex', c11: 209.7, c12: 121.1, c13: 105.1, c33: 210.9, c44: 42.5, rho: 5606, f: 'ZnO' }, GaN: { label: 'Gallium nitride (hexagonal)', sym: 'hex', c11: 390, c12: 145, c13: 106, c33: 398, c44: 105, rho: 6150, f: 'GaN' },
  custom: { label: 'Custom cubic (use the numeric inputs)', sym: 'cubic' },
};
function christoffel(inp) {
  const W = [], x = inp.xtal === 'custom' ? { ...XT.custom, c11: inp.cc11, c12: inp.cc12, c44: inp.cc44, rho: inp.rho_x } : (XT[inp.xtal] || XT.Si), g = 1e9, C6 = N.zeros(6, 6), hex = x.sym === 'hex';
  const c11 = x.c11 * g, c12 = x.c12 * g, c44 = x.c44 * g, c13 = hex ? x.c13 * g : c12, c33 = hex ? x.c33 * g : c11, c66 = hex ? (c11 - c12) / 2 : c44;
  C6[0][0] = C6[1][1] = c11; C6[2][2] = c33; C6[0][1] = C6[1][0] = c12; C6[0][2] = C6[2][0] = C6[1][2] = C6[2][1] = c13; C6[3][3] = C6[4][4] = c44; C6[5][5] = c66;
  const V = [[0, 5, 4], [5, 1, 3], [4, 3, 2]], gam = n => { const G = N.zeros(3, 3); for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { let s = 0; for (let j = 0; j < 3; j++) for (let l = 0; l < 3; l++) s += C6[V[i][j]][V[k][l]] * n[j] * n[l]; G[i][k] = s; } return G; };
  const solve = (th, ph) => { const n = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)], e = N.eigSym(gam(n)), v = e.values.map(l => Math.sqrt(Math.max(l, 0) / x.rho)), pol = e.vectors.map(u => Math.abs(u[0] * n[0] + u[1] * n[1] + u[2] * n[2])); return { v, pol, n }; };
  const th = inp.pol_ang * PI / 180, ph = inp.azi_ang * PI / 180, s0 = solve(th, ph), ths = N.linspace(0, 180, 181), sw = ths.map(a => solve(a * PI / 180, ph));
  // Voigt–Reuss–Hill isotropic averages
  let S6; const eig6 = N.eigSym(C6).values, stable = eig6[0] > 0; try { S6 = N.inverse(C6); } catch { S6 = N.zeros(6, 6); }
  const KV = (c11 + c11 + c33 + 2 * (c12 + 2 * c13)) / 9, GV = (c11 + c11 + c33 - (c12 + 2 * c13) + 3 * (2 * c44 + c66)) / 15, KR = 1 / (S6[0][0] + S6[1][1] + S6[2][2] + 2 * (S6[0][1] + S6[0][2] + S6[1][2])), GR = 15 / (4 * (S6[0][0] + S6[1][1] + S6[2][2]) - 4 * (S6[0][1] + S6[0][2] + S6[1][2]) + 3 * (S6[3][3] + S6[4][4] + S6[5][5]));
  const K = 0.5 * (KV + KR), G = 0.5 * (GV + GR), dn = Math.abs(3 * K + G) > 1e-6 ? 3 * K + G : 1e-6, E = 9 * K * G / dn, nu = (3 * K - 2 * G) / (2 * dn), vL = Math.sqrt(Math.max(K + 4 * G / 3, 0) / x.rho), vT = Math.sqrt(Math.max(G, 0) / x.rho), vD = vL > 0 && vT > 0 ? (1 / 3 * (1 / vL ** 3 + 2 / vT ** 3)) ** (-1 / 3) : 0, AZ = hex ? c44 / (c66 || 1e-30) : 2 * c44 / ((c11 - c12) || 1e-30), AU = 5 * GV / (GR || 1e-30) + KV / (KR || 1e-30) - 6;
  const vmaxL = N.maxOf(sw.map(s => s.v[2])), vminL = N.minOf(sw.map(s => s.v[2])), ref = hex ? Math.sqrt(c33 / x.rho) : Math.sqrt(c11 / x.rho), ax = solve(0, 0);
  if (!stable) W.push('The stiffness matrix is not positive definite: the Born stability conditions are violated and real wave speeds do not exist in some directions.');
  return {
    summary: `${x.label}: along the direction (polar ${inp.pol_ang}°, azimuth ${inp.azi_ang}°) the quasi-longitudinal speed is ${f3(s0.v[2])} m/s and the two quasi-shear speeds are ${f3(s0.v[1])} and ${f3(s0.v[0])} m/s. Elastic anisotropy factor ${f3(AZ)} (1 = isotropic). Polycrystal (Voigt–Reuss–Hill) averages: E = ${f3(E / 1e9)} GPa, G = ${f3(G / 1e9)} GPa, K = ${f3(K / 1e9)} GPa, ν = ${f3(nu)}, v_L = ${f3(vL)} m/s, v_T = ${f3(vT)} m/s.`,
    kpis: [
      kpi('vL', 'Quasi-longitudinal phase velocity', s0.v[2], 'm/s', 'Fastest root of the Christoffel equation in the chosen direction; its polarisation is within ' + f3(Math.acos(Math.min(s0.pol[2], 1)) * 180 / PI) + '° of the propagation direction.'),
      kpi('vT2', 'Fast quasi-shear velocity', s0.v[1], 'm/s', 'Second root: the faster of the two shear-like waves; the two shear speeds coincide along high-symmetry (acoustic) axes.'),
      kpi('vT1', 'Slow quasi-shear velocity', s0.v[0], 'm/s', 'Slowest root of the Christoffel equation; it dominates the low-temperature specific heat and thermal conductivity.'),
      kpi('Z', 'Longitudinal acoustic impedance', x.rho * s0.v[2], 'Pa·s/m', 'ρv for the quasi-longitudinal wave in this direction; governs transmission across interfaces and transducer matching.'),
      kpi('A', 'Elastic anisotropy factor', AZ, '–', hex ? 'Shear anisotropy C₄₄/C₆₆ of the hexagonal crystal (1 for isotropy in the prism planes).' : 'Zener ratio 2C₄₄/(C₁₁ − C₁₂): 1 for an isotropic solid, 1.56 for silicon, above 3 for copper.'),
      kpi('AU', 'Universal anisotropy index', AU, '–', '5G_V/G_R + K_V/K_R − 6: zero for an isotropic crystal, applicable to any symmetry.'),
      kpi('E', 'Polycrystalline Young\'s modulus (VRH)', E, 'Pa', 'Hill average of the Voigt (uniform strain) and Reuss (uniform stress) bounds for a random polycrystal.'),
      kpi('G', 'Polycrystalline shear modulus (VRH)', G, 'Pa', `Hill average shear modulus; the Voigt and Reuss bounds are ${f3(GV / 1e9)} and ${f3(GR / 1e9)} GPa.`),
      kpi('K', 'Bulk modulus (VRH)', K, 'Pa', 'Hill average bulk modulus (Voigt and Reuss coincide for cubic crystals).'),
      kpi('nu', 'Polycrystalline Poisson\'s ratio', nu, '–', 'From the Hill moduli, (3K − 2G)/(2(3K + G)).'),
      kpi('vD', 'Debye mean sound velocity', vD, 'm/s', 'Average [⅓(v_L⁻³ + 2v_T⁻³)]^(−1/3) of the isotropic speeds; with the atomic density it gives the Debye temperature.'),
      kpi('aniso_v', 'Longitudinal velocity anisotropy', (vmaxL - vminL) / Math.max(vminL, 1e-9), '–', 'Relative spread of the quasi-longitudinal speed over the scanned plane.'),
    ],
    series: [
      ser('v', 'Phase velocities against direction', 'polar angle from the z (c) axis (°)', 'phase velocity (m/s)', [{ name: 'quasi-longitudinal', x: ths, y: sw.map(s => s.v[2]) }, { name: 'fast quasi-shear', x: ths, y: sw.map(s => s.v[1]) }, { name: 'slow quasi-shear', x: ths, y: sw.map(s => s.v[0]) }], `Christoffel phase velocities in the plane containing the z axis at azimuth ${inp.azi_ang}°. Crossings of the shear branches mark acoustic axes; the variation of each branch measures the elastic anisotropy.`),
      ser('slow', 'Slowness curves', 'sin θ / v (µs/m)', 'cos θ / v (µs/m)', [{ name: 'quasi-longitudinal', x: sw.map(s => 1e6 * s.n[0] / Math.max(s.v[2], 1) / Math.max(Math.cos(ph), 1e-9) * Math.cos(ph)), y: sw.map(s => 1e6 * s.n[2] / Math.max(s.v[2], 1)) }, { name: 'slow quasi-shear', x: sw.map(s => 1e6 * Math.hypot(s.n[0], s.n[1]) / Math.max(s.v[0], 1)), y: sw.map(s => 1e6 * s.n[2] / Math.max(s.v[0], 1)) }], 'Slowness (inverse phase velocity) surfaces in the scanned plane. The energy-flow direction is normal to the slowness curve, so concave regions of the shear sheet produce phonon focusing.'),
    ],
    tables: [{ key: 'C', title: 'Stiffness matrix used (GPa, Voigt notation)', columns: ['', '1', '2', '3', '4', '5', '6'], rows: C6.map((r, i) => [String(i + 1), ...r.map(v => +(v / g).toPrecision(5))]), explain: 'Elastic stiffness constants in the crystal axes; hexagonal crystals have C₆₆ = (C₁₁ − C₁₂)/2 and the z axis along c.' }],
    equations: [...(x.f ? formulaEq(x.f, 'Crystal') : []),
      eqn('governing', 'Christoffel equation', '(Γ_ik − ρv²δ_ik) u_k = 0,  Γ_ik = C_ijkl n_j n_l'), eqn('constitutive', 'Generalised Hooke law', 'σ_ij = C_ijkl ε_kl'),
      eqn('formula', 'Born stability', hex ? 'C₁₁ > |C₁₂|,  C₃₃(C₁₁ + C₁₂) > 2C₁₃²,  C₄₄ > 0' : 'C₁₁ − C₁₂ > 0,  C₁₁ + 2C₁₂ > 0,  C₄₄ > 0'), eqn('formula', 'Voigt–Reuss–Hill averages', 'K = ½(K_V + K_R),  G = ½(G_V + G_R),  v_L = √((K + 4G/3)/ρ),  v_T = √(G/ρ)'),
    ],
    checks: [
      chk('Born mechanical stability', stable, `smallest stiffness eigenvalue ${f3(eig6[0] / g)} GPa`), chk('Axial longitudinal speed equals √(C/ρ)', Math.abs(ax.v[2] - ref) < 1e-6 * ref, `${f4(ax.v[2])} m/s along z against ${f4(ref)} m/s`),
      chk('Christoffel tensor trace invariant: ρΣv² = tr Γ', Math.abs(x.rho * N.sum(s0.v.map(v => v * v)) - (gam(s0.n)[0][0] + gam(s0.n)[1][1] + gam(s0.n)[2][2])) < 1e-6 * c11, 'sum of the three squared speeds matches the trace'),
      chk('Reuss ≤ Voigt bounds', GR <= GV * (1 + 1e-9) && KR <= KV * (1 + 1e-9), `G: ${f3(GR / g)} ≤ ${f3(GV / g)} GPa; K: ${f3(KR / g)} ≤ ${f3(KV / g)} GPa`), chk('Eigen-polarisations orthogonal (pure or quasi modes)', s0.pol[2] > 0.7, `longitudinal character ${f3(s0.pol[2])}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `Orient devices to the wave needed: in this plane the longitudinal speed varies by ${f3(100 * (vmaxL - vminL) / Math.max(vminL, 1e-9))} %, so resonator frequency and delay-line timing depend on crystal cut by the same amount.`, 'Direction dependence from the Christoffel solution'),
      rec('medium', `For polycrystalline or nanocrystalline films use the Hill averages (E = ${f3(E / 1e9)} GPa, ν = ${f3(nu)}); textured films (for example c-axis ZnO) should instead use the single-crystal value along the texture axis.`, 'Voigt–Reuss–Hill averaging'),
      rec('medium', 'Piezoelectric stiffening raises wave speeds in ZnO and GaN by a few per cent in some directions; include the electromechanical coupling for transducer and filter design.', 'Purely elastic Christoffel tensor used'),
      rec('low', 'Silicon, ZnO and AlN-family acoustic devices avoid lead-based piezoelectrics; choosing them where coupling suffices removes a regulated substance from the product.', 'Hazardous-substance substitution', true),
    ],
    publish: { 'prop.c_sound': vL, 'prop.E': E, 'prop.G': G, 'prop.K': K, 'prop.nu': nu, 'prop.density': x.rho, ...(x.f ? { 'material.formula': x.f } : {}) },
    numerics: { method: 'Jacobi eigen-solution of the 3×3 Christoffel tensor on a 1° direction sweep; dense inversion for Reuss averages', dof: 3, steps: ths.length, notes: 'Linear elasticity, no piezoelectric stiffening, no dispersion (long-wavelength limit).' },
  };
}

// ================================================================ nanobeam flexural resonances and impedance mismatch
const BCS = { cantilever: { label: 'Cantilever (clamped–free)', b: [1.8751, 4.6941, 7.8548], me: 0.2427 }, fixed: { label: 'Clamped–clamped (bridge)', b: [4.7300, 7.8532, 10.9956], me: 0.3965 }, pinned: { label: 'Pinned–pinned', b: [PI, 2 * PI, 3 * PI], me: 0.5 } };
function nanobeam(inp) {
  const W = [], E = inp.E_b, rho = inp.rho_b, L = inp.L_b * 1e-6, w = inp.w_b * 1e-9, t = inp.t_b * 1e-9, circ = inp.section === 'circular', bc = BCS[inp.bc_b] || BCS.fixed, T = inp.T_b + 273.15;
  const A = circ ? PI * t * t / 4 : w * t, I0 = circ ? PI * t ** 4 / 64 : w * t ** 3 / 12, EIraw = E * I0 + (circ ? PI * inp.Es_b * t ** 3 / 8 : inp.Es_b * (w * t * t / 2 + t ** 3 / 6)), EI = Math.max(EIraw, 0.01 * E * I0), m = rho * A * L, sig = inp.sig_r * 1e6;
  const tens = (n, l = L) => { if (inp.bc_b === 'cantilever') return 1; const g = inp.bc_b === 'pinned' ? 1 / ((n + 1) ** 2 * PI * PI) : [0.02458, 0.0121, 0.0072][n]; return Math.sqrt(Math.max(1 + g * sig * A * l * l / EI, 0)); };
  const fn = (n, l = L, ei = EI) => bc.b[n] ** 2 / (2 * PI * l * l) * Math.sqrt(ei / (rho * A)) * tens(n, l), f = [0, 1, 2].map(n => fn(n)), f0b = bc.b[0] ** 2 / (2 * PI * L * L) * Math.sqrt(E * I0 / (rho * A));
  const me = bc.me * m, keff = me * (2 * PI * f[0]) ** 2, resp = f[0] / (2 * me), xth = keff > 0 ? Math.sqrt(K_B * T / keff) : 0, Qf = inp.Q_b, lw = f[0] / Qf, dm = 2 * me / Qf * 0.01, Pcr = inp.bc_b === 'cantilever' ? 0 : (inp.bc_b === 'fixed' ? 4 : 1) * PI * PI * EI / (L * L), sbk = Pcr / A;
  const A1 = SOLID[inp.matA] || SOLID.Si, B1 = SOLID[inp.matB] || SOLID.SiO2, Z1 = A1.rho * A1.c, Z2 = B1.rho * B1.c, Tp = 4 * Z1 * Z2 / (Z1 + Z2) ** 2, Zm = Math.sqrt(Z1 * Z2), thc = A1.c < B1.c ? Math.asin(A1.c / B1.c) * 180 / PI : 90;
  if (L / t < 10) W.push(`Length-to-thickness ratio ${f3(L / t)} is below 10: Euler–Bernoulli theory overestimates the frequencies (shear and rotary inertia neglected), increasingly for higher modes.`);
  if (inp.bc_b !== 'cantilever' && sig < 0 && -sig >= sbk) W.push(`Compressive residual stress exceeds the Euler buckling stress (${f3(sbk / 1e6)} MPa): the beam is buckled and the straight-beam frequency is not meaningful.`);
  if (EIraw <= 0.01 * E * I0) W.push('Negative surface elasticity has driven the bending rigidity non-positive; it was floored at 1 % of the bulk value.');
  if (t < 10e-9) W.push('Below about 10 nm thickness, surface elasticity, surface stress and native oxide change the stiffness by tens of per cent; treat the result as an estimate.');
  const Ls = N.logspace(0.2 * L, 5 * L, 60), ts = N.logspace(0.2 * t, 5 * t, 60), xs = N.linspace(0, 1, 81);
  const shape = (n) => { const b = bc.b[n]; if (inp.bc_b === 'pinned') return xs.map(x => Math.sin(b * x)); const sgm = inp.bc_b === 'cantilever' ? (Math.cosh(b) + Math.cos(b)) / (Math.sinh(b) + Math.sin(b)) : (Math.cosh(b) - Math.cos(b)) / (Math.sinh(b) - Math.sin(b)), y = xs.map(x => Math.cosh(b * x) - Math.cos(b * x) - sgm * (Math.sinh(b * x) - Math.sin(b * x))), mx = N.maxOf(y.map(Math.abs)) || 1; return y.map(v => v / mx); };
  return {
    summary: `${circ ? `Nanowire of diameter ${inp.t_b} nm` : `Nanobeam ${inp.w_b} × ${inp.t_b} nm`}, length ${inp.L_b} µm, ${bc.label.toLowerCase()}: flexural resonances at ${f3(f[0] / 1e6)}, ${f3(f[1] / 1e6)} and ${f3(f[2] / 1e6)} MHz (bulk, stress-free fundamental ${f3(f0b / 1e6)} MHz). Effective mass ${f3(me * 1e18)} fg, stiffness ${f3(keff)} N/m, mass responsivity ${f3(resp * 1e-21)} Hz/ag, thermal amplitude ${f3(xth * 1e12)} pm at ${inp.T_b} °C. Across a ${A1.label}/${B1.label} interface ${f3(100 * Tp)} % of longitudinal acoustic power is transmitted.`,
    kpis: [
      kpi('f1', 'Fundamental flexural frequency', f[0], 'Hz', `First bending resonance β₁²/(2πL²)·√(EI/ρA) with β₁ = ${f3(bc.b[0])}, including surface elasticity and residual stress where applicable.`),
      kpi('f2', 'Second flexural frequency', f[1], 'Hz', `Second bending mode (β₂ = ${f3(bc.b[1])}); the ratio f₂/f₁ = ${f3(f[1] / f[0])} identifies the boundary condition and the tension state experimentally.`),
      kpi('f3', 'Third flexural frequency', f[2], 'Hz', `Third bending mode (β₃ = ${f3(bc.b[2])}).`.padEnd(60, ' ') + 'Higher modes are more sensitive to shear deformation.'),
      kpi('m_eff', 'Effective modal mass', me, 'kg', `${bc.me} × the beam mass for the fundamental mode referred to the point of maximum deflection.`),
      kpi('k_eff', 'Effective stiffness', keff, 'N/m', 'm_eff·(2πf₁)²: the spring constant of the equivalent lumped resonator.'),
      kpi('resp', 'Mass responsivity', resp, 'Hz/kg', 'Frequency shift per unit mass added at the antinode, f₁/(2m_eff); smaller beams are more sensitive.'),
      kpi('dm', 'Mass resolution at 1 % of linewidth', dm, 'kg', `Mass change that shifts the resonance by one hundredth of the linewidth f/Q for Q = ${Qf}; phase-locked detection typically reaches this level.`),
      kpi('lw', 'Resonance linewidth', lw, 'Hz', 'Full width f₁/Q of the resonance; Q is limited by clamping, surface and gas damping.'),
      kpi('x_th', 'Thermomechanical amplitude', xth, 'm', 'Root-mean-square thermal vibration √(kT/k_eff) of the fundamental mode; it is the noise floor of displacement detection and the natural calibration signal.'),
      kpi('surf', 'Stiffness change from surface elasticity', EI / (E * I0) - 1, '–', 'Relative change of bending rigidity from the surface elastic modulus; it scales as 1/thickness.'),
      kpi('T_ac', 'Acoustic power transmission at the interface', Tp, '–', `4Z₁Z₂/(Z₁ + Z₂)² for longitudinal waves at normal incidence between ${A1.label} and ${B1.label}; the remainder is reflected.`),
      kpi('Z_match', 'Matching-layer impedance', Zm, 'Pa·s/m', 'Impedance √(Z₁Z₂) of a quarter-wave layer that would give full transmission at its design frequency.'),
    ],
    series: [
      ser('fL', 'Resonance frequencies against length', 'length (µm)', 'frequency (MHz)', [0, 1, 2].map(n => ({ name: `mode ${n + 1}`, x: Ls.map(l => l * 1e6), y: Ls.map(l => pos(fn(n, l) / 1e6)) })), 'Flexural frequencies fall as 1/L² for a stress-free beam; tensile stress flattens the curve towards the 1/L string limit.', { logx: true, logy: true }),
      ser('shape', 'Mode shapes', 'position x/L (–)', 'normalised deflection (–)', [0, 1, 2].map(n => ({ name: `mode ${n + 1}`, x: xs, y: shape(n) })), 'Euler–Bernoulli eigenfunctions for the chosen supports. Mass sensing is most sensitive at antinodes and insensitive at nodes; strain sensing is most sensitive where curvature is greatest (the clamps).'),
      ser('ft', 'Fundamental frequency against thickness', 'thickness (nm)', 'frequency (MHz)', [{ name: 'with surface elasticity', x: ts.map(x => x * 1e9), y: ts.map(tt => { const a = circ ? PI * tt * tt / 4 : w * tt, i0 = circ ? PI * tt ** 4 / 64 : w * tt ** 3 / 12, ei = Math.max(E * i0 + (circ ? PI * inp.Es_b * tt ** 3 / 8 : inp.Es_b * (w * tt * tt / 2 + tt ** 3 / 6)), 0.01 * E * i0); return pos(bc.b[0] ** 2 / (2 * PI * L * L) * Math.sqrt(ei / (rho * a)) / 1e6); }) }, { name: 'bulk elasticity', x: ts.map(x => x * 1e9), y: ts.map(tt => pos(bc.b[0] ** 2 / (2 * PI * L * L) * Math.sqrt(E * (circ ? tt * tt / 16 : tt * tt / 12) / rho) / 1e6)) }], 'Frequency is proportional to thickness for a bulk-like beam; the surface term makes thin beams deviate from the straight line.', { logx: true, logy: true }),
    ],
    equations: [
      eqn('governing', 'Euler–Bernoulli beam equation', '(EI)* ∂⁴w/∂x⁴ − σA ∂²w/∂x² + ρA ∂²w/∂t² = 0'), eqn('formula', 'Flexural resonance', 'f_n = (β_n²/2πL²)·√((EI)*/ρA)·√(1 + γ_n σAL²/(EI)*)'),
      eqn('formula', 'Surface elasticity', circ ? '(EI)* = EπD⁴/64 + πE_sD³/8' : '(EI)* = Ewt³/12 + E_s(wt²/2 + t³/6)'), eqn('formula', 'Sensing figures', 'm_eff = c·ρAL,  ∂f/∂m = −f/(2m_eff),  ⟨x²⟩ = kT/k_eff'),
      eqn('formula', 'Acoustic mismatch', 'T = 4Z₁Z₂/(Z₁ + Z₂)²,  Z = ρc,  Z_match = √(Z₁Z₂)'),
    ],
    checks: [
      chk('Slender-beam validity (L/t ≥ 10)', L / t >= 10, `L/t = ${f3(L / t)}`), chk('Mode ratio matches the eigenvalue ratio', Math.abs(sig) > 0 || Math.abs(f[1] / f[0] - (bc.b[1] / bc.b[0]) ** 2) < 1e-9, `f₂/f₁ = ${f4(f[1] / f[0])} against (β₂/β₁)² = ${f4((bc.b[1] / bc.b[0]) ** 2)}`),
      chk('Not buckled by residual stress', inp.bc_b === 'cantilever' || sig >= 0 || -sig < sbk, inp.bc_b === 'cantilever' ? 'free end: no axial stress' : `σ_r = ${inp.sig_r} MPa, buckling at −${f3(sbk / 1e6)} MPa`),
      chk('Energy conservation at the interface (R + T = 1)', Math.abs(Tp + ((Z1 - Z2) / (Z1 + Z2)) ** 2 - 1) < 1e-12, `T = ${f4(Tp)}, R = ${f4(1 - Tp)}`), chk('Thermal amplitude far below thickness (linear regime)', xth < 0.1 * t, `x_th/t = ${f3(xth / t)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', `For higher mass sensitivity shorten the beam: responsivity scales as 1/L³·… here ${f3(resp * 1e-21)} Hz/ag; halving L raises f₁ fourfold and the responsivity eightfold, but also raises the required detection bandwidth.`, 'f ∝ t/L² and m_eff ∝ wtL'),
      rec('medium', inp.bc_b === 'cantilever' ? 'Cantilevers are immune to residual film stress; use them when the deposition stress is poorly controlled.' : `Residual stress of ${inp.sig_r} MPa changes f₁ by ${f3(100 * (tens(0) - 1))} %: measure wafer bow or f₂/f₁ (${f3(f[1] / f[0])} against ${f3((bc.b[1] / bc.b[0]) ** 2)} stress-free) to determine it.`, 'Tension stiffening of doubly clamped beams'),
      rec('medium', `Q = ${Qf}: operate in vacuum for nanobeams (gas damping dominates above about 1 mbar) and anneal or passivate surfaces, whose losses grow as the surface-to-volume ratio.`, 'Dissipation scaling in NEMS'),
      rec('low', `Only ${f3(100 * Tp)} % of acoustic power crosses the ${A1.label}/${B1.label} interface; a quarter-wave layer of impedance ${f3(Zm / 1e6)} MRayl improves coupling, whereas a large mismatch is useful to confine energy (higher Q) or reduce phonon heat flow.`, 'Acoustic impedance mismatch'),
      rec('low', 'NEMS resonant sensors operate at picowatt drive power and replace reagent-consuming assays in some detection tasks; silicon and nitride devices use no scarce elements.', 'Operating power and material footprint', true),
    ],
    publish: { 'res.f0': f[0], 'res.Q': Qf, 'geom.length': L, 'geom.thickness': t, 'res.transmission': Tp, 'mech.k_eff': keff },
    numerics: { method: 'Closed-form Euler–Bernoulli eigenfrequencies with surface-elasticity and tension corrections; acoustic-mismatch formula', dof: 3, notes: 'Uniform prismatic beam, ideal clamps, small amplitude; tension factors γ_n are standard approximations for clamped–clamped beams.' },
  };
}

const M = (...a) => ({ mode: a }), solidOpts = Object.entries(SOLID).map(([value, v]) => ({ value, label: v.label }));
export default {
  id: 15,
  tagline: 'Negative effective mass, bulk modulus and density in resonant acoustic media, phononic band gaps, sound insulation and absorption, crystal acoustics and nanobeam resonators.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'massinmass', group: 'Model', help: 'Selects the wave model that is solved; each analysis has its own inputs.', options: [
      { value: 'massinmass', label: 'Negative effective mass: mass-in-mass lattice' }, { value: 'helmholtz', label: 'Negative bulk modulus / density: Helmholtz resonators and membranes' }, { value: 'phononic', label: '1-D phononic crystal (layered)' },
      { value: 'panel', label: 'Sound insulation: mass law and locally resonant panel' }, { value: 'mpp', label: 'Sound absorption: micro-perforated panel, quarter-wave tubes' }, { value: 'christoffel', label: 'Elastic waves in crystals (Christoffel)' }, { value: 'nanobeam', label: 'Nanobeam flexural resonances and impedance mismatch' }] },
    { key: 'eta', label: 'Loss factor η', unit: '–', value: 0.02, min: 0, max: 0.5, group: 'Damping', showIf: M('massinmass', 'helmholtz', 'phononic', 'panel'), help: 'Structural (hysteretic) loss factor of springs, membranes or layers. It limits how negative the effective parameters become and smooths resonances.' },
    { key: 'Ncell', label: 'Number of cells / periods', unit: '–', value: 8, min: 1, max: 40, step: 1, group: 'Finite sample', showIf: M('massinmass', 'helmholtz', 'phononic'), help: 'Number of unit cells in the finite sample whose transmission is computed.' },
    { key: 'f_rel', label: 'Probe frequency (relative to resonance)', unit: '–', value: 1.08, min: 0.1, max: 5, group: 'Finite sample', showIf: M('massinmass', 'helmholtz', 'panel'), help: 'Frequency for single-point results as a multiple of the resonance frequency; values slightly above 1 lie in the negative band.' },
    { key: 'm1', label: 'Outer (host) mass m₁', unit: 'kg', value: 0.05, min: 1e-21, max: 1e3, log: true, group: 'Mass-in-mass lattice', showIf: M('massinmass'), help: 'Mass of the outer shell of each cell, connected to its neighbours.' },
    { key: 'm2', label: 'Inner (resonator) mass m₂', unit: 'kg', value: 0.1, min: 1e-21, max: 1e3, log: true, group: 'Mass-in-mass lattice', showIf: M('massinmass'), help: 'Mass of the internal resonator; the ratio m₂/m₁ sets the width of the negative-mass band.' },
    { key: 'k1', label: 'Lattice spring k₁', unit: 'N/m', value: 1e5, min: 1e-6, max: 1e12, log: true, group: 'Mass-in-mass lattice', showIf: M('massinmass'), help: 'Stiffness of the springs joining neighbouring cells; sets the host wave speed and Bragg cut-off.' },
    { key: 'k2', label: 'Resonator spring k₂', unit: 'N/m', value: 4e4, min: 1e-6, max: 1e12, log: true, group: 'Mass-in-mass lattice', showIf: M('massinmass'), help: 'Stiffness of the spring holding the inner mass; with m₂ it sets the local resonance.' },
    { key: 'a_cell', label: 'Lattice constant', unit: 'mm', value: 20, min: 1e-6, max: 1e4, log: true, group: 'Mass-in-mass lattice', showIf: M('massinmass'), help: 'Cell spacing; used for wave speed and effective density, not for the gap frequencies.' },
    { key: 'dn_type', label: 'Resonant elements', type: 'select', value: 'both', group: 'Resonant duct', showIf: M('helmholtz'), help: 'Helmholtz resonators give a negative bulk modulus, membranes a negative density; together they give a double-negative band.', options: [{ value: 'both', label: 'Helmholtz resonators + membranes (double negative)' }, { value: 'hr', label: 'Helmholtz resonators only (negative modulus)' }, { value: 'membrane', label: 'Membranes only (negative density)' }] },
    { key: 'fluid', label: 'Fluid', type: 'select', value: 'air', group: 'Resonant duct', showIf: M('helmholtz'), help: 'Fluid filling the duct and resonators.', options: Object.entries(FLUID).map(([value, v]) => ({ value, label: v.label })) },
    { key: 'duct_D', label: 'Duct diameter', unit: 'mm', value: 30, min: 0.1, max: 1000, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Diameter of the main duct; with the period it gives the duct volume per cell.' },
    { key: 'neck_r', label: 'Neck radius', unit: 'mm', value: 5, min: 0.01, max: 200, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Radius of the resonator neck; larger necks raise the resonance and lower the viscous loss.' },
    { key: 'neck_l', label: 'Neck length', unit: 'mm', value: 5, min: 0.01, max: 500, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Physical neck length before end corrections.' },
    { key: 'cav_V', label: 'Cavity volume', unit: 'cm³', value: 30, min: 1e-6, max: 1e5, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Volume of each resonator cavity; resonance frequency scales as V^(−1/2).' },
    { key: 'd_cell', label: 'Period along the duct', unit: 'mm', value: 40, min: 0.05, max: 2000, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Spacing of resonators/membranes; must be much shorter than the wavelength.' },
    { key: 'memb_T', label: 'Membrane tension', unit: 'N/m', value: 300, min: 0.1, max: 1e5, log: true, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Tension per unit length of each membrane; sets the cut-off below which the density is negative.' },
    { key: 'memb_sig', label: 'Membrane surface density', unit: 'kg/m²', value: 0.1, min: 0, max: 50, group: 'Resonant duct', showIf: M('helmholtz'), help: 'Mass per unit area of each membrane, added to the fluid inertia.' },
    { key: 'matA', label: 'Material A', type: 'select', value: 'steel', group: 'Layered media', showIf: M('phononic', 'nanobeam'), help: 'First layer of the phononic crystal (also the embedding medium), or the first medium of the interface in the nanobeam analysis.', options: solidOpts },
    { key: 'matB', label: 'Material B', type: 'select', value: 'epoxy', group: 'Layered media', showIf: M('phononic', 'nanobeam'), help: 'Second layer of the phononic crystal, or the second medium of the interface.', options: solidOpts },
    { key: 'dA', label: 'Thickness of layer A', unit: 'µm', value: 1000, min: 0.001, max: 1e6, log: true, group: 'Layered media', showIf: M('phononic'), help: 'Layer thickness (1000 µm = 1 mm; nanometre layers give gigahertz–terahertz gaps).' },
    { key: 'dB', label: 'Thickness of layer B', unit: 'µm', value: 1000, min: 0.001, max: 1e6, log: true, group: 'Layered media', showIf: M('phononic'), help: 'Thickness of the second layer.' },
    { key: 'E_p', label: 'Plate Young\'s modulus', unit: 'Pa', value: 7e10, min: 1e6, max: 1e12, log: true, group: 'Panel', bus: 'prop.E', showIf: M('panel'), help: 'Modulus of the host plate; enters only the coincidence frequency.' },
    { key: 'rho_p', label: 'Plate density', unit: 'kg/m³', value: 2700, min: 10, max: 22000, log: true, group: 'Panel', bus: 'prop.density', showIf: M('panel'), help: 'Density of the host plate; with thickness it gives the surface mass.' },
    { key: 'nu_p', label: 'Plate Poisson\'s ratio', unit: '–', value: 0.33, min: 0, max: 0.49, group: 'Panel', bus: 'prop.nu', showIf: M('panel'), help: 'Poisson\'s ratio of the plate (bending stiffness).' },
    { key: 'h_p', label: 'Plate thickness', unit: 'mm', value: 1, min: 0.001, max: 500, log: true, group: 'Panel', showIf: M('panel'), help: 'Thickness of the host plate.' },
    { key: 'mu_r', label: 'Resonator mass ratio', unit: '–', value: 0.3, min: 0.01, max: 5, log: true, group: 'Panel', showIf: M('panel'), help: 'Total resonator mass per area divided by the plate surface mass.' },
    { key: 'f_res', label: 'Resonator tuning frequency', unit: 'Hz', value: 300, min: 5, max: 20000, log: true, group: 'Panel', bus: 'res.f0', showIf: M('panel'), help: 'Natural frequency of the attached resonators.' },
    { key: 'd_hole', label: 'Hole diameter', unit: 'mm', value: 0.3, min: 0.02, max: 5, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Diameter of the micro-perforations; resistance scales as 1/d².' },
    { key: 't_mpp', label: 'Panel thickness', unit: 'mm', value: 0.3, min: 0.02, max: 20, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Thickness of the perforated sheet (hole length).' },
    { key: 'perf', label: 'Perforation ratio', unit: '%', value: 1, min: 0.05, max: 30, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Open-area fraction of the sheet (0.5–2 % for micro-perforated absorbers).' },
    { key: 'D_cav', label: 'Cavity depth', unit: 'mm', value: 30, min: 1, max: 2000, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Air gap between the sheet and the rigid wall; deeper cavities absorb lower frequencies.' },
    { key: 'L_qw', label: 'Quarter-wave tube length', unit: 'mm', value: 100, min: 1, max: 5000, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Length of the rigid-backed tubes of the comparison absorber.' },
    { key: 'r_qw', label: 'Quarter-wave tube radius', unit: 'mm', value: 2, min: 0.2, max: 100, log: true, group: 'Absorber', showIf: M('mpp'), help: 'Tube radius; narrow tubes add the wall loss needed for absorption.' },
    { key: 'phi_qw', label: 'Tube open-area ratio', unit: '%', value: 20, min: 0.5, max: 100, group: 'Absorber', showIf: M('mpp'), help: 'Fraction of the surface occupied by tube openings.' },
    { key: 'xtal', label: 'Crystal', type: 'select', value: 'ZnO', group: 'Crystal', showIf: M('christoffel'), help: 'Single crystal with handbook elastic constants; choose "custom" to enter cubic constants.', options: Object.entries(XT).map(([value, v]) => ({ value, label: v.label })) },
    { key: 'pol_ang', label: 'Polar angle from z (c) axis', unit: '°', value: 0, min: 0, max: 180, group: 'Crystal', showIf: M('christoffel'), help: 'Propagation direction measured from the crystal z axis ([001] cubic, c-axis hexagonal).' },
    { key: 'azi_ang', label: 'Azimuth from x axis', unit: '°', value: 0, min: 0, max: 360, group: 'Crystal', showIf: M('christoffel'), help: 'Azimuth of the propagation direction; also fixes the plane of the direction sweep (45° gives the (1−10) plane of cubic crystals).' },
    { key: 'cc11', label: 'Custom C₁₁', unit: 'GPa', value: 165.7, min: 1, max: 1500, group: 'Crystal', showIf: M('christoffel'), help: 'Used only for the custom cubic crystal.' },
    { key: 'cc12', label: 'Custom C₁₂', unit: 'GPa', value: 63.9, min: -100, max: 1000, group: 'Crystal', showIf: M('christoffel'), help: 'Used only for the custom cubic crystal.' },
    { key: 'cc44', label: 'Custom C₄₄', unit: 'GPa', value: 79.6, min: 0.1, max: 1000, group: 'Crystal', showIf: M('christoffel'), help: 'Used only for the custom cubic crystal.' },
    { key: 'rho_x', label: 'Custom density', unit: 'kg/m³', value: 2329, min: 100, max: 23000, group: 'Crystal', bus: 'prop.density', showIf: M('christoffel'), help: 'Used only for the custom cubic crystal.' },
    { key: 'E_b', label: 'Beam Young\'s modulus', unit: 'Pa', value: 1.69e11, min: 1e6, max: 1.2e12, log: true, group: 'Nanobeam', bus: 'prop.E', showIf: M('nanobeam'), help: 'Modulus along the beam axis (169 GPa for Si ⟨110⟩, about 140 GPa for ZnO c-axis nanowires).' },
    { key: 'rho_b', label: 'Beam density', unit: 'kg/m³', value: 2329, min: 100, max: 23000, group: 'Nanobeam', bus: 'prop.density', showIf: M('nanobeam'), help: 'Mass density of the beam material.' },
    { key: 'section', label: 'Cross-section', type: 'select', value: 'rect', group: 'Nanobeam', showIf: M('nanobeam'), help: 'Rectangular beam (width × thickness) or circular nanowire (diameter = thickness input).', options: [{ value: 'rect', label: 'Rectangular' }, { value: 'circular', label: 'Circular (nanowire)' }] },
    { key: 'L_b', label: 'Beam length', unit: 'µm', value: 10, min: 0.05, max: 5000, log: true, group: 'Nanobeam', showIf: M('nanobeam'), help: 'Free length between supports; frequency scales as 1/L².' },
    { key: 'w_b', label: 'Beam width', unit: 'nm', value: 300, min: 5, max: 1e6, log: true, group: 'Nanobeam', showIf: M('nanobeam'), help: 'In-plane width (rectangular section); it does not affect the out-of-plane frequency.' },
    { key: 't_b', label: 'Thickness / diameter', unit: 'nm', value: 100, min: 2, max: 1e6, log: true, group: 'Nanobeam', showIf: M('nanobeam'), help: 'Thickness in the vibration direction, or the nanowire diameter; frequency is proportional to it.' },
    { key: 'bc_b', label: 'Supports', type: 'select', value: 'fixed', group: 'Nanobeam', showIf: M('nanobeam'), help: 'Boundary conditions of the beam.', options: Object.entries(BCS).map(([value, v]) => ({ value, label: v.label })) },
    { key: 'sig_r', label: 'Residual axial stress', unit: 'MPa', value: 0, min: -2000, max: 5000, group: 'Nanobeam', showIf: M('nanobeam'), help: 'Tensile (positive) or compressive stress in a doubly supported beam; tension raises the frequencies.' },
    { key: 'Es_b', label: 'Surface elastic modulus', unit: 'N/m', value: 0, min: -50, max: 200, group: 'Nanobeam', showIf: M('nanobeam'), help: 'Surface elasticity (force per length) for the size-dependent bending rigidity; zero gives the bulk result.' },
    { key: 'Q_b', label: 'Quality factor', unit: '–', value: 5000, min: 1, max: 1e7, log: true, group: 'Nanobeam', showIf: M('nanobeam'), help: 'Measured or expected mechanical quality factor (10–100 in air, 10³–10⁵ in vacuum).' },
    { key: 'T_b', label: 'Temperature', unit: '°C', value: 25, min: -270, max: 800, group: 'Nanobeam', bus: 'proc.T_C', showIf: M('nanobeam'), help: 'Temperature for the thermomechanical noise amplitude.' },
  ],
  presets: [
    { name: 'Mass-in-mass metamaterial bar (100 Hz local resonance)', note: 'Negative effective mass between 100 and 174 Hz with a mass ratio of 2.', values: { mode: 'massinmass', m1: 0.05, m2: 0.1, k1: 1e5, k2: 4e4, a_cell: 20, eta: 0.02, Ncell: 8, f_rel: 1.08 } },
    { name: 'Helmholtz-resonator and membrane duct (double-negative acoustics)', note: 'Air-filled duct with side resonators (negative modulus) and membranes (negative density).', values: { mode: 'helmholtz', dn_type: 'both', fluid: 'air', duct_D: 30, neck_r: 5, neck_l: 5, cav_V: 30, d_cell: 40, memb_T: 300, memb_sig: 0.1, Ncell: 8, f_rel: 1.15 } },
    { name: 'Steel/epoxy phononic crystal (ultrasonic Bragg gap)', note: 'Millimetre layers give a wide stop band near 0.9 MHz.', values: { mode: 'phononic', matA: 'steel', matB: 'epoxy', dA: 1000, dB: 1000, Ncell: 6, eta: 0.01 } },
    { name: 'GaAs/AlAs nanoscale superlattice (sub-terahertz phonon mirror)', note: 'Nanometre periods fold acoustic phonons and open gigahertz–terahertz gaps.', values: { mode: 'phononic', matA: 'GaAs', matB: 'AlAs', dA: 0.004, dB: 0.004, Ncell: 20, eta: 0.001 } },
    { name: 'ZnO single crystal: wave speeds and polycrystal averages', note: 'Hexagonal ZnO elastic constants; reference case for the nanorod material.', values: { mode: 'christoffel', xtal: 'ZnO', pol_ang: 0, azi_ang: 0 } },
    { name: 'Silicon NEMS bridge resonator', note: 'Doubly clamped 10 µm × 300 nm × 100 nm silicon nanobeam.', values: { mode: 'nanobeam', E_b: 1.69e11, rho_b: 2329, section: 'rect', L_b: 10, w_b: 300, t_b: 100, bc_b: 'fixed', Q_b: 5000, matA: 'Si', matB: 'SiO2' } },
    { name: 'Locally resonant sound-insulating panel (300 Hz)', note: '1 mm aluminium sheet with 30 % added resonator mass.', values: { mode: 'panel', E_p: 7e10, rho_p: 2700, h_p: 1, mu_r: 0.3, f_res: 300, eta: 0.03, f_rel: 1.08 } },
    { name: 'Micro-perforated panel absorber', note: '0.3 mm holes, 1 % perforation, 30 mm cavity: fibre-free broadband absorber.', values: { mode: 'mpp', d_hole: 0.3, t_mpp: 0.3, perf: 1, D_cav: 30 } },
  ],
  run(inp) {
    const f = { massinmass: massInMass, helmholtz, phononic, panel, mpp, christoffel, nanobeam }[inp.mode] || massInMass;
    return finish(f(inp));
  },
};
