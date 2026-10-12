// Suite 11 — Electrical, Electronic and Electrochemical Engineering.
// Reduced-order models: Butler–Volmer/Tafel electrodeposition with a Nernst diffusion layer, steady 1-D Nernst–Planck
// transport, Randles impedance, semiconductor junction/nanostructure statistics, percolation conductivity, and
// negative differential resistance (transfer-matrix resonant tunnelling diode, Esaki tunnel diode).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const { PI, K_B, Q_E, EPS0, HBAR, M_E, H_PLANCK } = N, F = N.FARADAY, R = N.R_GAS;
const f3 = v => N.fmt(v, 3), f4 = v => N.fmt(v, 4);
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

// ---------------------------------------------------------------- electrodeposition
const COUPLES = {
  Zn: { label: 'Zn²⁺/Zn', ion: 'Zn2+', n: 2, nu: 1, prod: 'Zn', rho: 7140, E0: -0.762, salt: 'ZnSO4', acid: 'H2SO4', j0H: 3e-7, cR: ['Zn2+', 'e-'], cP: ['Zn'], metal: 'Zn' },
  Cu: { label: 'Cu²⁺/Cu', ion: 'Cu2+', n: 2, nu: 1, prod: 'Cu', rho: 8960, E0: 0.340, salt: 'CuSO4', acid: 'H2SO4', j0H: 2e-3, cR: ['Cu2+', 'e-'], cP: ['Cu'], metal: 'Cu' },
  Ni: { label: 'Ni²⁺/Ni', ion: 'Ni2+', n: 2, nu: 1, prod: 'Ni', rho: 8908, E0: -0.257, salt: 'NiCl2', acid: 'HCl', j0H: 6e-2, cR: ['Ni2+', 'e-'], cP: ['Ni'], metal: 'Ni' },
  Ag: { label: 'Ag⁺/Ag', ion: 'Ag+', n: 1, nu: 1, prod: 'Ag', rho: 10490, E0: 0.7996, salt: 'AgNO3', acid: 'HNO3', j0H: 1e-3, cR: ['Ag+', 'e-'], cP: ['Ag'], metal: 'Ag' },
  Au: { label: 'Au³⁺/Au (as AuCl₄⁻)', ion: 'AuCl4-', n: 3, nu: 1, prod: 'Au', rho: 19300, E0: 1.002, salt: 'HAuCl4', acid: 'HCl', j0H: 3e-2, cR: ['AuCl4-', 'e-'], cP: ['Au', 'Cl-'], metal: null },
  ZnO_NO3: { label: 'ZnO via nitrate reduction', ion: 'Zn2+', n: 2, nu: 1, prod: 'ZnO', rho: 5606, E0: 0.01, j0H: 1e-5, cR: ['Zn2+', 'NO3-', 'e-'], cP: ['ZnO', 'NO2-'], eR: ['NO3-', 'H2O', 'e-'], eP: ['NO2-', 'OH-'], metal: 'Zn' },
  ZnO_O2: { label: 'ZnO via oxygen reduction', ion: 'Zn2+', n: 2, nu: 1, prod: 'ZnO', rho: 5606, E0: 0.401, j0H: 1e-5, cR: ['Zn2+', 'O2', 'e-'], cP: ['ZnO'], eR: ['O2', 'H2O', 'e-'], eP: ['OH-'], metal: 'Zn' },
};
const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a || 1; };
// Combine two balanced half-reactions (each containing e-) into the overall cell reaction with electrons cancelled.
function combineHalves(cR, cP, aR, aP) {
  const net = new Map(); let ne = [0, 0];
  [[cR, cP], [aR, aP]].forEach(([r, p], h) => { const b = chem.balance(r, p).coefficients; [...r, ...p].forEach((s, i) => { if (s === 'e-') ne[h] = b[i]; }); });
  const L = ne[0] * ne[1] / gcd(ne[0], ne[1]);
  [[cR, cP], [aR, aP]].forEach(([r, p], h) => { const b = chem.balance(r, p).coefficients, k = L / ne[h]; [...r, ...p].forEach((s, i) => { if (s !== 'e-') net.set(s, (net.get(s) || 0) + (i < r.length ? 1 : -1) * b[i] * k); }); });
  const ent = [...net].filter(([, c]) => Math.abs(c) > 1e-9); if (!ent.length) return { text: null, electrons: L };
  const g = ent.reduce((a, [, c]) => gcd(a, Math.round(Math.abs(c))), 0) || 1;
  const side = sgn => ent.filter(([, c]) => c * sgn > 0).map(([s, c]) => (Math.abs(c) / g === 1 ? '' : Math.abs(c) / g + ' ') + chem.prettyFormula(s)).join(' + ');
  return { text: `${side(1)} → ${side(-1)}`, electrons: L / g };
}
const bvExp = x => Math.exp(Math.min(x, 600));

function deposit(inp) {
  const cp = COUPLES[inp.couple] || COUPLES.Zn, T = inp.T + 273.15, fRT = F / (R * T), W = [];
  const c = inp.c_ion * 1000, delta = inp.delta * 1e-6, A = inp.area * 1e-4, eta = inp.eta, ac = inp.alpha_c, aa = Math.max(cp.n - ac, 0.1);
  const isOx = inp.couple.startsWith('ZnO'), aOH = 10 ** (inp.pH - 14), cO2 = inp.cO2;
  // Nernst equilibrium potential (activities ≈ concentrations in mol/L)
  let Eeq, nernstTxt;
  if (inp.couple === 'ZnO_NO3') { Eeq = cp.E0 - (1 / (2 * fRT)) * (Math.log(1e-3) + 2 * Math.log(aOH)); nernstTxt = 'E = E° − (RT/2F)·ln(a(NO₂⁻)·a(OH⁻)²/a(NO₃⁻)), with a(NO₂⁻)/a(NO₃⁻) taken as 10⁻³'; }
  else if (inp.couple === 'ZnO_O2') { Eeq = cp.E0 - (1 / fRT) * Math.log(aOH) + (1 / (4 * fRT)) * Math.log(Math.max(cO2, 1e-6) / 1.3); nernstTxt = 'E = E° − (RT/F)·ln a(OH⁻) + (RT/4F)·ln(c(O₂)/1.3 mM)'; }
  else if (inp.couple === 'Au') { Eeq = cp.E0 + (1 / (3 * fRT)) * Math.log(Math.max(inp.c_ion, 1e-12)); nernstTxt = 'E = E° + (RT/3F)·ln(a(AuCl₄⁻)/a(Cl⁻)⁴), with a(Cl⁻) = 1 (excess chloride)'; }
  else { Eeq = cp.E0 + (1 / (cp.n * fRT)) * Math.log(Math.max(inp.c_ion, 1e-12)); nernstTxt = `E = E° + (RT/${cp.n}F)·ln a(${chem.prettyFormula(cp.ion)})`; }
  // mass-transport limit across the Nernst diffusion layer
  let jlim = cp.n * F * inp.D_ion * c / delta, limBy = chem.prettyFormula(cp.ion);
  if (inp.couple === 'ZnO_O2') { const jO = 4 * F * 2.4e-9 * cO2 / delta; if (jO < jlim) { jlim = jO; limBy = 'dissolved O₂'; } }
  const jOf = e => { const k = inp.j0 * (bvExp(ac * fRT * e) - bvExp(-aa * fRT * e)); return k >= 0 ? k / (1 + inp.j0 * bvExp(ac * fRT * e) / jlim) : k; };
  const jkin = inp.j0 * (bvExp(ac * fRT * eta) - bvExp(-aa * fRT * eta)), j = jOf(eta);
  const Eel = Eeq - eta, EH = -(Math.LN10 / fRT) * inp.pH, etaH = EH - Eel;
  const jHf = eH => (eH > 0 ? cp.j0H * (bvExp(0.5 * fRT * eH) - bvExp(-0.5 * fRT * eH)) : 0);
  const jH = Math.min(jHf(etaH), 1e7), jtot = j + jH, CE = jtot > 0 ? j / jtot : 0;
  const M = chem.molarMass(cp.prod) / 1000, rate = j * M * cp.nu / (cp.n * F * cp.rho), thick = rate * inp.t_dep, mass = rate * inp.t_dep * cp.rho * A;
  const Qc = jtot * A * inp.t_dep, tafel = Math.LN10 / (ac * fRT), Rct = 1 / (cp.n * fRT * inp.j0);
  const eSpec = j > 0 ? inp.Vcell * jtot / (j * M * cp.nu / (cp.n * F)) / 3.6e6 : 0;
  const etaC = j < jlim && j > 0 ? -(1 / (cp.n * fRT)) * Math.log(1 - j / jlim) : 0;
  if (j / jlim > 0.8) W.push(`Operating at ${f3(100 * j / jlim)} % of the limiting current: deposits become rough, powdery or dendritic near the mass-transport limit.`);
  if (CE < 0.8 && eta > 0) W.push(`Hydrogen evolution takes ${f3(100 * (1 - CE))} % of the current at this potential and pH (order-of-magnitude exchange current density for the substrate).`);
  if (inp.T > 100 || inp.T < 0) W.push('Temperature is outside the liquid range of an aqueous bath at ambient pressure; kinetic and transport parameters are not temperature-corrected.');
  if (eta === 0) W.push('Zero overpotential: the electrode is at equilibrium and no net deposition occurs.');
  const es = N.linspace(0, Math.max(1.6 * eta, 0.3), 120);
  const series = [
    ser('pol', 'Cathodic polarisation curve', 'overpotential η (V)', 'current density (A/m²)', [
      { name: 'Butler–Volmer (kinetic only)', x: es, y: es.map(e => Math.min(inp.j0 * (bvExp(ac * fRT * e) - bvExp(-aa * fRT * e)), 5 * jlim)) },
      { name: 'with mass transport', x: es, y: es.map(jOf) }, { name: 'limiting current', x: [es[0], es[119]], y: [jlim, jlim] }],
      'Deposition current density against cathodic overpotential. The kinetic curve (clipped at five times the limit for readability) rises exponentially; the mixed-control curve saturates at the diffusion-limited current. Operate on the rising part, well below the plateau.'),
    ser('tafel', 'Tafel plot', 'overpotential η (V)', 'log₁₀ j (A/m²)', [{ name: 'log j (mixed control)', x: es.slice(1), y: es.slice(1).map(e => Math.log10(pos(jOf(e)))) }, { name: 'hydrogen evolution', x: es.slice(1), y: es.slice(1).map(e => Math.log10(pos(jHf(EH - (Eeq - e))) + 1e-12)) }],
      `Logarithm of current density against overpotential. The linear region has the Tafel slope ${f3(tafel * 1000)} mV per decade; curvature at high η marks the onset of mass-transport control. The second line is the parasitic hydrogen current.`),
    ser('thick', 'Deposit thickness against time', 'time (s)', 'thickness (µm)', [{ name: 'Faraday-law thickness', x: N.linspace(0, inp.t_dep, 50), y: N.linspace(0, inp.t_dep, 50).map(t => rate * t * 1e6) }],
      'Compact-deposit thickness from Faraday\'s law at the computed partial current density; it is linear in time for potentiostatic deposition with a steady diffusion layer.'),
  ];
  // chemistry generated from the selected couple
  const equations = [], checks = [];
  const half = (r, p, label, note) => { try { const b = chem.balance(r, p); equations.push(eqn('half-cell', label, chem.equationString(r, p, b.coefficients), note)); const a = chem.auditEquation(r, p, b.coefficients); checks.push(chk(`${label}: mass and charge balance`, a.massBalanced && a.chargeBalanced, `coefficients ${b.coefficients.join(', ')}; net charge residual ${a.netCharge}`)); } catch (e) { W.push(`${label}: ${e.message}`); } };
  if (cp.eR) { half(cp.eR, cp.eP, 'Cathode electron-transfer step', `E° = ${cp.E0} V vs SHE`); try { equations.push(eqn('ionic', 'Precipitation at the cathode surface', chem.equationString(['Zn2+', 'OH-'], ['ZnO', 'H2O']), 'local pH rise drives ZnO formation')); } catch (e) { W.push(e.message); } }
  half(cp.cR, cp.cP, cp.eR ? 'Overall cathode half-cell' : 'Cathode half-cell (reduction)', cp.eR ? 'electron transfer + precipitation' : `E° = ${cp.E0} V vs SHE`);
  const sol = inp.anode === 'soluble' && cp.metal, aR = sol ? [cp.metal] : ['H2O'], aP = sol ? [cp.metal + '2+'.replace('2+', cp.metal === 'Ag' ? '+' : '2+'), 'e-'] : ['O2', 'H+', 'e-'];
  half(aR, aP, sol ? 'Anode half-cell (soluble anode)' : 'Anode half-cell (inert anode, water oxidation)');
  if (inp.anode === 'soluble' && !cp.metal) W.push('A soluble gold anode does not dissolve to AuCl₄⁻ in this simple scheme; an inert anode was assumed.');
  try { const ov = combineHalves(cp.cR, cp.cP, aR, aP); equations.push(eqn('chemical', 'Overall cell reaction', ov.text || `${chem.prettyFormula(cp.metal)}(anode) → ${chem.prettyFormula(cp.metal)}(cathode)`, ov.text ? `${ov.electrons} electrons transferred per formula equation` : 'metal is transferred from anode to cathode; bath composition is unchanged')); } catch (e) { W.push('Overall cell reaction: ' + e.message); }
  if (cp.salt && !sol) { try { const io = chem.ionicEquations([{ f: cp.salt, phase: 'aq' }, { f: 'H2O', phase: 'l' }], [{ f: cp.prod, phase: 's' }, { f: 'O2', phase: 'g' }, { f: cp.acid, phase: 'aq' }]); equations.push(eqn('chemical', 'Molecular equation (inert anode)', io.molecular), eqn('ionic', 'Complete ionic equation', io.completeIonic), eqn('ionic', 'Net ionic equation', io.netIonic, io.spectators.length ? 'spectators: ' + io.spectators.join(', ') : undefined)); } catch (e) { W.push('Ionic equations: ' + e.message); } }
  equations.push(eqn('formula', 'Deposit', `${chem.prettyFormula(cp.prod)}, M = ${f4(M * 1000)} g/mol, ρ = ${cp.rho} kg/m³`),
    eqn('constitutive', 'Nernst equation', nernstTxt), eqn('governing', 'Butler–Volmer with surface concentration', 'j = j₀[(1 − j/j_lim)·exp(α_c F η/RT) − exp(−α_a F η/RT)]', `α_a = n − α_c = ${f3(aa)}`),
    eqn('constitutive', 'Diffusion-limited current (Nernst layer)', 'j_lim = n F D c / δ'), eqn('formula', 'Faraday law', 'thickness = j·M·t / (n F ρ);  m = Q·M·CE / (n F)'), eqn('formula', 'Tafel slope', 'b = 2.303 RT / (α_c F)'));
  checks.push(chk('Current below mass-transport limit', j <= jlim * (1 + 1e-9), `j/j_lim = ${f4(j / jlim)} (limited by ${limBy})`),
    chk('Faraday-law mass consistency', Math.abs(mass - Qc * CE * M * cp.nu / (cp.n * F)) <= 1e-9 * Math.max(mass, 1e-30), `m = Q·CE·M/(nF) reproduces ${f4(mass * 1e6)} mg`),
    chk('Smooth-deposit window (j < 0.8 j_lim)', j < 0.8 * jlim, `fraction of limiting current ${f3(j / jlim)}`),
    chk('Current efficiency above 90 %', CE > 0.9, `CE = ${f4(100 * CE)} %`),
    chk('Electrode potential positive of solvent-reduction runaway', etaH < 0.8, `hydrogen overpotential ${f3(etaH)} V at pH ${inp.pH}`));
  return {
    summary: `${cp.label} deposition at η = ${f3(eta)} V and ${inp.T} °C runs at ${f3(j)} A/m² (${f3(100 * j / jlim)} % of the ${f3(jlim)} A/m² diffusion limit), giving ${f3(rate * 1e9)} nm/s and ${f3(thick * 1e6)} µm of ${chem.prettyFormula(cp.prod)} in ${f3(inp.t_dep)} s with ${f3(100 * CE)} % current efficiency. The equilibrium potential is ${f3(Eeq)} V vs SHE and the working electrode sits at ${f3(Eel)} V.`,
    kpis: [
      kpi('E_eq', 'Equilibrium (Nernst) potential', Eeq, 'V vs SHE', `Reversible potential of the ${cp.label} couple at the bath composition from the Nernst equation; deposition requires the electrode to be held negative of this value.`),
      kpi('E_electrode', 'Working electrode potential', Eel, 'V vs SHE', 'Equilibrium potential minus the applied cathodic overpotential; this is the potential to set on a potentiostat (convert to your reference electrode).'),
      kpi('j', 'Deposition current density', j, 'A/m²', 'Partial current density of the deposition reaction from the Butler–Volmer equation corrected for depletion of the electroactive species across the Nernst diffusion layer.', j < 0.8 * jlim),
      kpi('j_lim', 'Limiting current density', jlim, 'A/m²', `Maximum current the diffusion layer can supply, n·F·D·c/δ, limited here by ${limBy}; stirring or a thinner boundary layer raises it.`),
      kpi('j_frac', 'Fraction of limiting current', j / jlim, '–', 'Ratio j/j_lim; values below about 0.3–0.5 give compact deposits, values near 1 give rough or dendritic growth.', j / jlim < 0.8),
      kpi('rate', 'Deposition rate', rate * 1e9, 'nm/s', 'Thickness growth rate of a fully dense deposit from Faraday\'s law at the deposition partial current density.'),
      kpi('thickness', 'Deposit thickness', thick * 1e6, 'µm', `Compact thickness after ${f3(inp.t_dep)} s; porous or nanostructured deposits are thicker in proportion to their porosity.`),
      kpi('mass', 'Deposited mass', mass * 1e6, 'mg', `Mass on ${inp.area} cm² from the charge passed multiplied by the current efficiency (Faraday's law).`),
      kpi('CE', 'Current efficiency', 100 * CE, '%', 'Share of the total cathodic current that deposits product; the remainder is hydrogen evolution estimated from a Tafel law for the substrate.', CE > 0.9),
      kpi('tafel', 'Cathodic Tafel slope', tafel * 1000, 'mV/decade', 'Overpotential needed for a tenfold increase of the kinetic current, 2.303·RT/(α_c·F); fit this to a measured polarisation curve to calibrate α_c.'),
      kpi('Rct', 'Charge-transfer resistance', Rct, 'Ω·m²', 'Small-signal resistance RT/(nF·j₀) at equilibrium; it is the diameter of the Nyquist semicircle per unit area.'),
      kpi('eta_conc', 'Concentration overpotential', etaC, 'V', 'Part of the applied overpotential spent on depletion of the surface concentration, −(RT/nF)·ln(1 − j/j_lim).'),
      kpi('energy', 'Specific electrical energy', eSpec, 'kWh/kg', `Electrical energy per kilogram of deposit at a cell voltage of ${inp.Vcell} V including the current-efficiency loss.`),
    ],
    series, equations, checks, warnings: W,
    recommendations: [
      rec(j / jlim > 0.5 ? 'high' : 'medium', j / jlim > 0.5 ? `Lower the overpotential or increase agitation: at ${f3(100 * j / jlim)} % of j_lim the deposit will roughen. Halving δ (rotating electrode, flow, ultrasound) doubles j_lim to ${f3(2 * jlim)} A/m².` : `The process is under kinetic control (${f3(100 * j / jlim)} % of j_lim); thickness uniformity is then governed by the primary/secondary current distribution, so keep electrode spacing uniform.`, 'Ratio of operating to diffusion-limited current density'),
      rec('medium', `Calibrate j₀ and α_c from a measured Tafel plot: the model predicts a slope of ${f3(tafel * 1000)} mV/decade and an intercept of ${f3(inp.j0)} A/m².`, 'Butler–Volmer parameters are bath- and additive-specific'),
      rec(CE < 0.9 ? 'high' : 'low', CE < 0.9 ? `Raise the pH or reduce the overpotential to suppress hydrogen evolution (currently ${f3(100 * (1 - CE))} % of the current); hydrogen also embrittles and pits the deposit.` : 'Hydrogen evolution is minor at this potential; monitor pH drift at the cathode during long runs.', 'Competing hydrogen evolution reaction'),
      rec('medium', `Use pulsed deposition (on-time below the transition time) if finer grains are needed: pulse peak currents above the DC limit of ${f3(jlim)} A/m² are possible because the pulse diffusion layer is thinner.`, 'Mass-transport limit of DC plating'),
      rec('medium', `Electrical demand is ${f3(eSpec)} kWh/kg; recover drag-out and recirculate the bath, and prefer a soluble anode where possible to avoid oxygen evolution and acid build-up.`, 'Energy and bath-life accounting for this run', true),
    ],
    publish: { 'elec.current_density': j, 'growth.rate': rate, 'geom.thickness': thick, 'elec.E_eq': Eeq, 'elec.current_efficiency': CE, 'rxn.yield_kg': mass, 'material.formula': cp.prod, 'fab.energy_kWh_kg': eSpec },
    numerics: { method: 'Closed-form Butler–Volmer with first-order surface-concentration correction; Nernst and Faraday laws', dof: 1, steps: es.length, notes: 'Steady diffusion layer; no ohmic drop, no double-layer charging, compact deposit assumed.' },
  };
}

// ---------------------------------------------------------------- ionic transport
const SYS = {
  zno: { name: 'Zn²⁺ / OH⁻ (ZnO reference case)', cat: 'Zn2+', an: 'OH-', zc: 2, za: -1 },
  nacl: { name: 'NaCl', cat: 'Na+', an: 'Cl-', zc: 1, za: -1, salt: 'NaCl' }, kcl: { name: 'KCl', cat: 'K+', an: 'Cl-', zc: 1, za: -1, salt: 'KCl' },
  cuso4: { name: 'CuSO₄', cat: 'Cu2+', an: 'SO4^2-', zc: 2, za: -2, salt: 'CuSO4' }, hcl: { name: 'HCl', cat: 'H+', an: 'Cl-', zc: 1, za: -1, salt: 'HCl' }, agno3: { name: 'AgNO₃', cat: 'Ag+', an: 'NO3-', zc: 1, za: -1, salt: 'AgNO3' },
};
const bern = p => (Math.abs(p) < 1e-6 ? 1 + p / 2 : p / (1 - Math.exp(-p)));
function ionic(inp) {
  const s = SYS[inp.system] || SYS.zno, T = inp.T + 273.15, W = [], E = inp.Efield, L = inp.Lgap * 1e-3, eps = inp.eps_r * EPS0;
  const ions = [{ f: s.cat, z: s.zc, c: inp.c_cat, D: inp.D_cat }, { f: s.an, z: s.za, c: inp.c_an, D: inp.D_an }];
  const sigma = (F * F / (R * T)) * N.sum(ions.map(i => i.z * i.z * i.D * i.c)), I = 0.5 * N.sum(ions.map(i => i.z * i.z * i.c));
  const debye = I > 0 ? Math.sqrt(eps * R * T / (2 * F * F * I)) : 1e-3, joule = sigma * E * E, jMig = sigma * E;
  ions.forEach(i => { i.mob = Math.abs(i.z) * i.D * F / (R * T); i.u = i.z * i.D * F * E / (R * T); i.Pe = i.z * F * E * L / (R * T); i.t = sigma > 0 ? (F * F / (R * T)) * i.z * i.z * i.D * i.c / sigma : 0; });
  const cat = ions[0], an = ions[1], enh = bern(cat.Pe), Jd = cat.D * cat.c / L, Jc = Jd * enh, jc = cat.z * F * Jc;
  const neutral = cat.z * cat.c + an.z * an.c, xs = N.linspace(0, 1, 121);
  const profC = xs.map(x => (cat.Pe > 0 ? (1 - Math.exp(-cat.Pe * (1 - x))) / (1 - Math.exp(-cat.Pe)) : Math.abs(cat.Pe) < 1e-9 ? 1 - x : (Math.exp(cat.Pe) - Math.exp(cat.Pe * x)) / (Math.exp(cat.Pe) - 1)));
  const profA = xs.map(x => Math.exp(Math.max(an.Pe * x, -700) - Math.max(an.Pe, 0)));
  const dTdt = joule / 4.18e6, tDiff = L * L / cat.D, tMig = Math.abs(cat.u) > 0 ? L / Math.abs(cat.u) : Infinity, tSS = Math.min(tDiff, tMig);
  if (L / debye > 100) W.push(`The gap is ${f3(L / debye)} Debye lengths wide: in a real cell most of the applied potential drops in the electrode double layers and the bulk field is screened unless current flows; the uniform-field assumption describes the conducting bulk only.`);
  if (Math.abs(neutral) > 1e-6 * Math.max(Math.abs(cat.z * cat.c), 1e-30)) W.push(`The entered concentrations are not electroneutral (Σ zᵢcᵢ = ${f3(neutral)} mol/m³); a counter-ion not listed must balance the charge and would add to the conductivity.`);
  if (dTdt > 0.1) W.push(`Joule heating of ${f3(joule)} W/m³ raises the temperature by ${f3(dTdt)} K/s if not removed; convection and property drift follow.`);
  if (I > 100) W.push('Ionic strength exceeds about 0.1 mol/L: the dilute Nernst–Einstein relation overestimates conductivity (ion pairing and relaxation effects are neglected).');
  const Es = N.logspace(10, 1e6, 60), equations = [];
  try { equations.push(eqn('ionic', 'Ions considered', `${chem.prettyFormula(s.cat)} (z = ${s.zc}, M = ${f4(chem.molarMass(s.cat))} g/mol) and ${chem.prettyFormula(s.an)} (z = ${s.za}, M = ${f4(chem.molarMass(s.an))} g/mol)`)); } catch (e) { W.push(e.message); }
  try {
    if (s.salt) equations.push(eqn('ionic', 'Dissociation', chem.equationString([s.salt], [s.cat, s.an])));
    else { const io = chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'H2O', phase: 'l' }, { f: 'NaNO3', phase: 'aq' }]); equations.push(eqn('chemical', 'Molecular equation (reference precipitation)', io.molecular), eqn('ionic', 'Complete ionic equation', io.completeIonic), eqn('ionic', 'Net ionic equation', io.netIonic, 'spectators: ' + io.spectators.join(', ')), eqn('half-cell', 'Cathodic base generation (if electro-assisted)', chem.equationString(['H2O', 'e-'], ['H2', 'OH-']))); }
  } catch (e) { W.push('Chemistry: ' + e.message); }
  equations.push(eqn('governing', 'Nernst–Planck flux (steady, 1-D, uniform field)', 'Jᵢ = −Dᵢ dcᵢ/dx + (zᵢF/RT)·Dᵢ·E·cᵢ = constant'), eqn('formula', 'Flux with field', 'J = (D c₀/L)·Pe/(1 − e^(−Pe)),  Pe = zFEL/RT'),
    eqn('constitutive', 'Nernst–Einstein conductivity', 'σ = (F²/RT) Σ zᵢ² Dᵢ cᵢ'), eqn('formula', 'Joule heating', 'q = σE²'), eqn('formula', 'Debye length', 'λ_D = √(ε RT / (F² Σ zᵢ² cᵢ))'));
  return {
    summary: `For ${s.name} at ${inp.T} °C the electrolyte conductivity is ${f3(sigma)} S/m, so a field of ${f3(E)} V/m drives ${f3(jMig)} A/m² and dissipates ${f3(joule)} W/m³ as Joule heat. Electromigration changes the cation flux to a sink ${f3(inp.Lgap)} mm away by a factor ${f3(enh)} relative to diffusion alone (Pe = ${f3(cat.Pe)}). The Debye length is ${f3(debye * 1e9)} nm.`,
    kpis: [
      kpi('sigma', 'Ionic conductivity', sigma, 'S/m', 'Dilute-solution conductivity from the Nernst–Einstein relation summed over the two ions; compare with a conductivity-meter reading to check the diffusivities.'),
      kpi('joule', 'Joule heating power density', joule, 'W/m³', 'Volumetric heat release σE² in the bulk electrolyte; multiply by the electrolyte volume for the heat load on the thermostat.', dTdt < 0.1),
      kpi('dTdt', 'Adiabatic heating rate', dTdt, 'K/s', 'Temperature rise rate if no heat is removed, using the volumetric heat capacity of water (4.18 MJ/m³/K).', dTdt < 0.1),
      kpi('j_mig', 'Ohmic current density', jMig, 'A/m²', 'Bulk migration current density σE carried by both ions when the electrodes can sustain the current.'),
      kpi('debye', 'Debye length', debye * 1e9, 'nm', 'Thickness of the diffuse double layer; fields are screened over this distance from a charged surface, so nanostructures closer than about 2 λ_D have overlapping double layers.'),
      kpi('Pe', 'Electromigration Péclet number', cat.Pe, '–', 'Ratio of electromigration to diffusion for the cation, zFEL/RT (the potential drop across the gap in thermal-voltage units); |Pe| ≫ 1 means the field dominates transport.'),
      kpi('enh', 'Cation flux enhancement', enh, '–', 'Steady Nernst–Planck flux to a perfect sink divided by the purely diffusive flux D·c₀/L; values below 1 mean the field opposes transport.', enh >= 1),
      kpi('flux', 'Cation flux to the sink', Jc, 'mol/m²/s', 'Molar flux of the cation arriving at the consuming surface from the analytical solution of the steady Nernst–Planck equation.'),
      kpi('j_cat', 'Equivalent cation current density', jc, 'A/m²', 'Charge carried by the cation flux, z·F·J; it bounds the growth or deposition rate that transport can sustain.'),
      kpi('u_cat', 'Cation drift velocity', cat.u, 'm/s', 'Electromigration velocity z·D·F·E/RT of the cation; the anion drifts the opposite way at ' + f3(an.u) + ' m/s.'),
      kpi('t_cat', 'Cation transference number', cat.t, '–', 'Fraction of the migration current carried by the cation; the remainder is carried by the anion.'),
      kpi('t_ss', 'Time to steady profile', tSS, 's', 'Smaller of the diffusion time L²/D and the drift transit time L/u for the cation; the steady profile applies after a few multiples of this time.'),
      kpi('C_dl', 'Diffuse-layer capacitance', eps / debye, 'F/m²', 'Low-potential Debye–Hückel capacitance ε/λ_D of the diffuse layer; the Stern layer in series lowers the measured value.'),
    ],
    series: [
      ser('prof', 'Steady concentration profiles under the field', 'position x/L (–)', 'c/c₀ (–)', [{ name: chem.prettyFormula(s.cat) + ' (sink at x = L)', x: xs, y: profC }, { name: chem.prettyFormula(s.an) + ' (blocked at x = L)', x: xs, y: profA }, { name: 'diffusion only', x: [0, 1], y: [1, 0] }],
        'Normalised concentration across the gap. With a favourable field the cation profile stays flat and drops in a thin layer at the sink (thickness L/Pe); the co-ion at a blocking wall follows a Boltzmann distribution, normalised to its maximum.'),
      ser('enhE', 'Flux enhancement and Joule heating against field', 'electric field (V/m)', 'value', [{ name: 'flux enhancement (–)', x: Es, y: Es.map(e => bern(cat.z * F * e * L / (R * T))) }, { name: 'Joule heating (kW/m³)', x: Es, y: Es.map(e => sigma * e * e / 1000) }],
        'Trade-off between transport enhancement (linear in field once Pe ≫ 1) and Joule heating (quadratic in field). The useful operating field is where the enhancement is significant but the heat load is still removable.', { logx: true, logy: true }),
    ],
    tables: [{ key: 'ions', title: 'Ion transport properties', columns: ['ion', 'z', 'c (mol/m³)', 'D (m²/s)', 'mobility (m²/V/s)', 'drift velocity (m/s)', 'transference'], rows: ions.map(i => [chem.prettyFormula(i.f), i.z, i.c, i.D, +i.mob.toExponential(3), +i.u.toExponential(3), +i.t.toFixed(4)]), explain: 'Per-ion mobility from the Einstein relation, drift velocity in the applied field and share of the migration current.' }],
    equations,
    checks: [
      chk('Electroneutrality of the entered composition', Math.abs(neutral) <= 1e-6 * Math.max(Math.abs(cat.z * cat.c), 1e-30), `Σ zᵢcᵢ = ${f3(neutral)} mol/m³`),
      chk('Transference numbers sum to one', Math.abs(cat.t + an.t - 1) < 1e-9 || sigma === 0, `t₊ + t₋ = ${f4(cat.t + an.t)}`),
      chk('Flux reduces to Fick\'s law at zero field', Math.abs(bern(0) - 1) < 1e-12, 'Pe → 0 limit of Pe/(1 − e^(−Pe)) equals 1'),
      chk('Dilute-solution validity (I < 0.1 mol/L)', I <= 100, `ionic strength ${f3(I / 1000)} mol/L`),
      chk('Continuum double layer thinner than the gap', debye < L, `λ_D/L = ${f3(debye / L)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', enh > 2 ? `The field gives a ${f3(enh)}-fold flux gain; confirm the gain experimentally by comparing growth rate with the field off, because electrode polarisation usually screens part of the applied voltage.` : `At ${f3(E)} V/m over ${f3(inp.Lgap)} mm the field barely changes transport (Pe = ${f3(cat.Pe)}); raise the potential drop above a few RT/zF (≈ ${f3(1000 * R * T / (Math.abs(cat.z) * F))} mV each) to see an effect.`, 'Electromigration Péclet number'),
      rec('medium', `Budget ${f3(joule / 1000)} kW/m³ of cooling, or pulse the field: Joule heat scales with E² whereas the flux gain scales with E.`, 'σE² heating against linear flux enhancement'),
      rec('medium', `Keep feature spacing above about ${f3(4 * debye * 1e9)} nm (4 λ_D) if independent double layers are needed, or below it to exploit overlap and ion selectivity.`, 'Debye screening length'),
      rec('low', 'Measure the bath conductivity and back-calculate an effective diffusivity; replace the handbook values if they differ by more than 20 %.', 'Nernst–Einstein relation is exact only at infinite dilution'),
      rec('medium', `Lower-field, longer-duration operation cuts electrical energy use roughly in proportion to E: the present setting consumes ${f3(joule * 3600 / 3.6e6)} kWh per m³ of electrolyte per hour.`, 'Energy intensity of field-assisted processing', true),
    ],
    publish: { 'eff.sigma': sigma, 'proc.E_field': E, 'elec.current_density': jMig, 'elec.joule_W_m3': joule, 'elec.debye_length': debye, 'elec.flux_enhancement': enh },
    numerics: { method: 'Analytical steady Nernst–Planck solution (exponential/Bernoulli function) with Nernst–Einstein mobilities', dof: xs.length, notes: 'Uniform imposed field (supporting-electrolyte limit); Poisson coupling and convection not solved.' },
  };
}

// ---------------------------------------------------------------- impedance (Randles)
function eis(inp) {
  const T = inp.T + 273.15, W = [], Rs = inp.Rs, Rct = inp.Rct, Q = inp.Cdl * 1e-6, n = inp.n_cpe, sw = inp.sigma_w, A = inp.area * 1e-4;
  const f1 = Math.min(inp.fmin, inp.fmax / 10), fs = N.logspace(f1, inp.fmax, 220);
  const Zf = f => { const w = 2 * PI * f, zw = sw / Math.sqrt(w), fr = Rct + zw, fi = -zw, d = fr * fr + fi * fi, yr = Q * w ** n * Math.cos(n * PI / 2) + fr / d, yi = Q * w ** n * Math.sin(n * PI / 2) - fi / d, dd = yr * yr + yi * yi; return { re: Rs + yr / dd, im: -yi / dd }; };
  const Z = fs.map(Zf), mag = Z.map(z => Math.hypot(z.re, z.im)), ph = Z.map(z => Math.atan2(z.im, z.re) * 180 / PI);
  const fc = 1 / (2 * PI * (Rct * Q) ** (1 / n)), Ceff = Q ** (1 / n) * (1 / Rs + 1 / Rct) ** ((n - 1) / n), wW = 2 * sw * sw / (Rct * Rct), fW = wW / (2 * PI);
  const j0 = R * T / (F * Rct * A), c = inp.c_ion * 1000, Dw = sw > 0 ? (Math.SQRT2 * R * T / (F * F * A * sw * c)) ** 2 : 0, iPh = N.argmin(ph);
  if (inp.fmin > inp.fmax / 10) W.push('The frequency window was widened to one decade because f_min was not below f_max/10.');
  if (fc > inp.fmax || fc < f1) W.push(`The charge-transfer relaxation at ${f3(fc)} Hz lies outside the swept window; the semicircle apex is not resolved.`);
  if (n < 0.8) W.push('CPE exponent below 0.8 indicates a strongly distributed or porous electrode; a transmission-line model would be more appropriate than a Randles circuit.');
  return {
    summary: `The Randles circuit has a solution resistance of ${f3(Rs)} Ω and a charge-transfer arc of ${f3(Rct)} Ω that relaxes at ${f3(fc)} Hz (effective double-layer capacitance ${f3(Ceff * 1e6)} µF). Diffusion (Warburg) overtakes charge transfer below about ${f3(fW)} Hz. The implied exchange current density is ${f3(j0)} A/m² for a one-electron reaction on ${inp.area} cm².`,
    kpis: [
      kpi('f_c', 'Charge-transfer characteristic frequency', fc, 'Hz', 'Apex frequency of the Nyquist semicircle, 1/(2π(R_ct·Q)^(1/n)); the measurement must extend a decade either side to resolve R_ct.'),
      kpi('tau', 'Charge-transfer time constant', 1 / (2 * PI * fc), 's', 'Relaxation time of the interface, (R_ct·Q)^(1/n); potential steps shorter than this mainly charge the double layer.'),
      kpi('C_eff', 'Effective double-layer capacitance', Ceff * 1e6, 'µF', 'Capacitance equivalent to the constant-phase element by the Brug formula; divide by the true surface area to judge roughness (20–60 µF/cm² is typical for a smooth metal).'),
      kpi('C_area', 'Capacitance per geometric area', Ceff / A, 'F/m²', 'Effective capacitance divided by the geometric electrode area; ratios to about 0.2–0.6 F/m² estimate the roughness factor of a nanostructured electrode.'),
      kpi('j0', 'Exchange current density from R_ct', j0, 'A/m²', 'j₀ = RT/(F·R_ct·A) for a one-electron step; divide by n for multi-electron reactions. It measures the intrinsic kinetics of the electrode.'),
      kpi('f_W', 'Warburg crossover frequency', fW, 'Hz', 'Frequency below which the Warburg modulus exceeds R_ct (ω = 2σ²/R_ct²); below it the response is diffusion controlled with a 45° line.'),
      kpi('D_w', 'Diffusivity implied by the Warburg coefficient', Dw, 'm²/s', `Diffusion coefficient that reproduces σ_W for equal oxidised/reduced concentrations of ${inp.c_ion} mol/L and one electron; a consistency check on σ_W.`, Dw > 1e-12 && Dw < 1e-7),
      kpi('Z_low', 'Impedance modulus at lowest frequency', mag[0], 'Ω', `|Z| at ${f3(f1)} Hz; approximates the polarisation resistance seen by slow processes such as corrosion or deposition.`),
      kpi('phase_min', 'Most capacitive phase angle', ph[iPh], '°', `Minimum phase, reached at ${f3(fs[iPh])} Hz; −90° is an ideal capacitor and −45° a Warburg element.`),
    ],
    series: [
      ser('nyq', 'Nyquist plot', 'Re Z (Ω)', '−Im Z (Ω)', [{ name: 'Randles circuit', x: Z.map(z => z.re), y: Z.map(z => -z.im) }], 'Imaginary against real impedance with frequency decreasing to the right. The semicircle diameter is R_ct, its left intercept is R_s, and the 45° tail is semi-infinite diffusion.'),
      ser('bodeM', 'Bode magnitude', 'frequency (Hz)', '|Z| (Ω)', [{ name: '|Z|', x: fs, y: mag }], 'Impedance modulus against frequency. The high-frequency plateau is R_s and the mid-frequency step equals R_ct; a slope of −1/2 at low frequency is the Warburg element.', { logx: true, logy: true }),
      ser('bodeP', 'Bode phase', 'frequency (Hz)', 'phase (°)', [{ name: 'phase', x: fs, y: ph }], 'Phase angle against frequency. The dip marks the double-layer/charge-transfer relaxation; a plateau near −45° at low frequency confirms diffusion control.', { logx: true }),
    ],
    equations: [
      eqn('governing', 'Randles impedance', 'Z(ω) = R_s + 1 / [Q(jω)ⁿ + 1/(R_ct + Z_W)]'), eqn('constitutive', 'Warburg element (semi-infinite diffusion)', 'Z_W = σ_W (1 − j)/√ω'),
      eqn('formula', 'Charge-transfer resistance', 'R_ct = RT/(n F j₀ A)'), eqn('formula', 'Brug effective capacitance', 'C_eff = Q^(1/n) (1/R_s + 1/R_ct)^((n−1)/n)'),
      eqn('half-cell', 'Generic redox couple probed', chem.equationString(['Fe3+', 'e-'], ['Fe2+']), 'example one-electron outer-sphere couple used for the j₀ and D estimates'),
    ],
    checks: [
      chk('Passivity: Re Z > 0 at all frequencies', Z.every(z => z.re > 0), `min Re Z = ${f3(N.minOf(Z.map(z => z.re)))} Ω`),
      chk('Capacitive response: Im Z ≤ 0', Z.every(z => z.im <= 1e-12), 'no inductive loop in a Randles circuit'),
      chk('High-frequency limit tends to R_s', Math.abs(Zf(1e12).re - Rs) < 0.02 * Rs + 1e-9, `Re Z(10¹² Hz) = ${f4(Zf(1e12).re)} Ω`),
      chk('Monotonic Re Z (Kramers–Kronig-consistent form)', Z.every((z, i) => i === 0 || z.re <= Z[i - 1].re + 1e-9), 'Re Z decreases with frequency as required for a causal RC ladder'),
    ],
    warnings: W,
    recommendations: [
      rec('high', `Sweep from at least ${f3(fc * 100)} Hz down to ${f3(Math.min(fc, fW) / 100)} Hz to resolve both the charge-transfer arc and the diffusion tail.`, 'Characteristic frequencies of this circuit'),
      rec('medium', `Minimise uncompensated resistance: R_s is ${f3(100 * Rs / (Rs + Rct))} % of the DC resistance, so move the reference electrode closer or raise the supporting-electrolyte concentration.`, 'Ratio R_s/(R_s + R_ct)'),
      rec('medium', 'Validate measured spectra with a Kramers–Kronig test before fitting; drifting electrodes violate linearity and give meaningless fitted parameters.', 'Linear, causal, stable-system requirement of impedance analysis'),
      rec('low', `A CPE exponent of ${n} ${n < 0.9 ? 'suggests surface heterogeneity; report Q and n rather than a single capacitance' : 'is close to ideal; a pure capacitor can be fitted'}.`, 'Constant-phase element interpretation'),
      rec('low', 'Use small excitation amplitudes (5–10 mV) and short multi-sine acquisitions: this reduces electrode consumption and electrolyte waste per characterisation.', 'Measurement economy', true),
    ],
    publish: { 'elec.R_ct': Rct, 'res.f0': fc, 'elec.C_dl': Ceff, 'elec.j0': j0 },
    numerics: { method: 'Closed-form complex impedance of the Randles circuit on a logarithmic frequency grid', dof: fs.length, notes: 'Semi-infinite Warburg; finite-length diffusion not included.' },
  };
}

// ---------------------------------------------------------------- semiconductor nanostructures
const SEMI = {
  Si: { Eg: 1.12, me: 1.08, mh: 0.81, bm: [0.26, 0.39], eps: 11.7, mun: 0.135, mup: 0.048 }, Ge: { Eg: 0.661, me: 0.56, mh: 0.29, bm: [0.12, 0.21], eps: 16.2, mun: 0.39, mup: 0.19 },
  GaAs: { Eg: 1.424, me: 0.067, mh: 0.47, eps: 12.9, mun: 0.85, mup: 0.04 }, ZnO: { Eg: 3.37, me: 0.24, mh: 0.59, eps: 8.5, mun: 0.02, mup: 0.002 },
  GaN: { Eg: 3.39, me: 0.20, mh: 1.4, eps: 8.9, mun: 0.10, mup: 0.003 }, CdSe: { Eg: 1.74, me: 0.13, mh: 0.45, eps: 10.2, mun: 0.065, mup: 0.004 }, PbS: { Eg: 0.41, me: 0.085, mh: 0.085, eps: 17.2, mun: 0.06, mup: 0.07 },
};
function semi(inp) {
  const m = SEMI[inp.semi_mat] || SEMI.Si, T = inp.T + 273.15, kT = K_B * T, Vt = kT / Q_E, W = [], eps = m.eps * EPS0;
  const Nd = inp.Nd * 1e6, Na = inp.Na * 1e6, dos = (ms, t) => 2 * (2 * PI * ms * M_E * K_B * t / (H_PLANCK * H_PLANCK)) ** 1.5;
  const Nc = dos(m.me, T), Nv = dos(m.mh, T), niT = t => Math.sqrt(dos(m.me, t) * dos(m.mh, t)) * Math.exp(-m.Eg * Q_E / (2 * K_B * t)), ni = niT(T);
  const n0 = Nd / 2 + Math.sqrt(Nd * Nd / 4 + ni * ni), p0 = ni * ni / n0, sigma = Q_E * (n0 * m.mun + p0 * m.mup);
  const EcEf = Vt * Math.log(Nc / n0), LD = Math.sqrt(eps * Vt / (Q_E * n0)), schottky = inp.contact === 'schottky';
  let Vbi = schottky ? inp.phiB - EcEf : Vt * Math.log(Math.max(Na * Nd / (ni * ni), 1.0001));
  if (Vbi <= 0) { W.push('The barrier height is below the bulk Fermi-level offset: the contact is ohmic, not rectifying; depletion results are nominal.'); Vbi = Vt; }
  const Va = Math.min(inp.Vbias, Vbi - Vt), Neff = schottky ? Nd : Na * Nd / (Na + Nd), wdep = V => Math.sqrt(2 * eps * Math.max(Vbi - V, Vt) / (Q_E * Neff)), Wd = wdep(Va), Emax = 2 * (Vbi - Va) / Wd, Cj = eps / Wd;
  if (inp.Vbias > Vbi - Vt) W.push('Forward bias approaches the built-in potential; the depletion approximation fails and the bias was clipped for the width calculation.');
  const Dn = m.mun * Vt, Dp = m.mup * Vt, Js = schottky ? 1.2e6 * m.me * T * T * Math.exp(-inp.phiB / Vt) : Q_E * ni * ni * (Math.sqrt(Dn / inp.tau) / Na + Math.sqrt(Dp / inp.tau) / Nd);
  // nanowire surface depletion (cylindrical depletion approximation)
  const Rw = inp.d_nw * 0.5e-9, rcrit = Math.sqrt(4 * eps * inp.psi_s / (Q_E * Nd)), fullDep = Rw <= rcrit;
  let rc = 0; if (!fullDep) rc = N.bisect(r => (Q_E * Nd / (4 * eps)) * (Rw * Rw - r * r - 2 * r * r * Math.log(Rw / r)) - inp.psi_s, Rw * 1e-6, Rw, Rw * 1e-9);
  const sigNW = sigma * (rc / Rw) ** 2;
  // Brus quantum-confinement gap
  const [bme, bmh] = m.bm || [m.me, m.mh], mu = 1 / (1 / bme + 1 / bmh), aB = 4 * PI * eps * HBAR * HBAR / (mu * M_E * Q_E * Q_E);
  const brus = r => m.Eg + (HBAR * HBAR * PI * PI / (2 * r * r * M_E)) * (1 / bme + 1 / bmh) / Q_E - 1.786 * Q_E / (4 * PI * eps * r);
  const Rq = inp.R_qd * 1e-9, Eqd = brus(Rq), lam = 1239.84 / Eqd;
  if (Nd > 0.1 * Nc) W.push(`Donor density is ${f3(Nd / Nc)} of the conduction-band effective density of states: the semiconductor is (nearly) degenerate and Boltzmann statistics underestimate the Fermi level.`);
  if (Rq < 1e-9) W.push('Quantum-dot radius below 1 nm: the effective-mass (Brus) model strongly overestimates the confinement energy.');
  if (Eqd < m.Eg) W.push('The Coulomb term exceeds the confinement term at this radius; the dot is in the weak-confinement regime where the Brus formula is only qualitative.');
  const Ts = N.linspace(Math.max(150, T - 200), T + 400, 80), Rs = N.linspace(Math.max(0.8e-9, 0.3 * Rq), Math.max(4 * Rq, 8e-9), 90), Vs = N.linspace(-5, Math.max(Va, 0), 60), Vf = N.linspace(0, Math.min(Math.max(schottky ? inp.phiB : Vbi, 0.2), 1.5), 60);
  const nT = t => { const nI = niT(t); return Nd / 2 + Math.sqrt(Nd * Nd / 4 + nI * nI); };
  return {
    summary: `${inp.semi_mat} doped n-type at ${f3(inp.Nd)} cm⁻³ has n = ${f3(n0 * 1e-6)} cm⁻³ (intrinsic density ${f3(ni * 1e-6)} cm⁻³) and conductivity ${f3(sigma)} S/m at ${inp.T} °C. The ${schottky ? 'Schottky contact' : 'p–n junction'} has a built-in potential of ${f3(Vbi)} V and a depletion width of ${f3(Wd * 1e9)} nm. A ${inp.d_nw} nm nanowire is ${fullDep ? 'fully depleted by' : 'partly depleted by'} a ${inp.psi_s} V surface barrier, and a quantum dot of radius ${inp.R_qd} nm has a Brus gap of ${f3(Eqd)} eV (${f3(lam)} nm).`,
    kpis: [
      kpi('ni', 'Intrinsic carrier density', ni * 1e-6, 'cm⁻³', 'Thermally generated electron (and hole) density √(N_c N_v)·exp(−E_g/2kT) from parabolic-band effective densities of states.'),
      kpi('n', 'Electron density (n-side)', n0 * 1e-6, 'cm⁻³', 'Majority-carrier density from charge neutrality with fully ionised donors; equals N_D unless the material is close to intrinsic.'),
      kpi('p', 'Hole density (n-side)', p0 * 1e-6, 'cm⁻³', 'Minority-carrier density from the mass-action law n·p = nᵢ²; it controls diffusion leakage and photoresponse.'),
      kpi('sigma', 'Bulk conductivity', sigma, 'S/m', 'q(nµₙ + pµₚ) with handbook bulk mobilities; nanostructures usually show lower mobility because of surface and grain-boundary scattering.'),
      kpi('EcEf', 'Fermi level below conduction band', EcEf, 'eV', 'E_c − E_F = kT·ln(N_c/n); values below about 3kT signal degeneracy.', EcEf > 3 * Vt),
      kpi('LD', 'Debye screening length', LD * 1e9, 'nm', 'Distance over which mobile carriers screen a potential perturbation, √(εkT/q²n); it sets the minimum feature size for bulk-like behaviour.'),
      kpi('Vbi', 'Built-in potential', Vbi, 'V', schottky ? 'Band bending at the Schottky contact, φ_B − (E_c − E_F)/q.' : 'Contact potential of the p–n junction, (kT/q)·ln(N_A N_D/nᵢ²).'),
      kpi('W', 'Depletion width', Wd * 1e9, 'nm', `Space-charge region width at ${f3(Va)} V bias from the abrupt-junction depletion approximation.`),
      kpi('Emax', 'Peak junction field', Emax, 'V/m', 'Maximum electric field at the metallurgical junction, 2(V_bi − V)/W; compare with the breakdown field of the material.'),
      kpi('Cj', 'Junction capacitance', Cj, 'F/m²', 'Depletion capacitance ε/W per unit area; a Mott–Schottky plot of 1/C² against bias returns the doping density.'),
      kpi('Js', 'Saturation current density', Js, 'A/m²', schottky ? 'Thermionic-emission saturation current A*T²·exp(−qφ_B/kT) with the Richardson constant scaled by the effective mass.' : 'Shockley diffusion saturation current q·nᵢ²(√(Dₙ/τ)/N_A + √(Dₚ/τ)/N_D).'),
      kpi('r_crit', 'Nanowire full-depletion radius', rcrit * 1e9, 'nm', 'Largest radius that a surface barrier ψ_s depletes completely, √(4εψ_s/qN_D); thinner wires are insulating unless gated or illuminated.', !fullDep),
      kpi('sigma_nw', 'Effective nanowire conductivity', sigNW, 'S/m', 'Bulk conductivity scaled by the undepleted core cross-section fraction (r_core/R)² from the cylindrical depletion solution.'),
      kpi('Eg_qd', 'Quantum-dot band gap (Brus)', Eqd, 'eV', `Size-dependent gap including kinetic confinement and electron–hole Coulomb attraction; the bulk gap is ${m.Eg} eV.`),
      kpi('lambda_qd', 'Quantum-dot absorption edge', lam, 'nm', 'Wavelength of the first excitonic transition, 1239.84/E_g(R); it blue-shifts as the dot shrinks.'),
      kpi('aB', 'Exciton Bohr radius', aB * 1e9, 'nm', 'Hydrogenic exciton radius 4πεħ²/(µq²); dots smaller than this are in the strong-confinement regime where the Brus model applies.', Rq < aB),
    ],
    series: [
      ser('nT', 'Carrier densities against temperature', 'temperature (K)', 'density (cm⁻³)', [{ name: 'n (doped)', x: Ts, y: Ts.map(t => pos(nT(t) * 1e-6)) }, { name: 'nᵢ', x: Ts, y: Ts.map(t => pos(niT(t) * 1e-6)) }], 'Majority and intrinsic carrier densities. The doped density is flat in the extrinsic (saturation) range and follows nᵢ once the intrinsic density exceeds the doping; freeze-out at low temperature is not modelled.', { logy: true }),
      ser('brus', 'Quantum-confined band gap against radius', 'dot radius (nm)', 'band gap (eV)', [{ name: 'Brus model', x: Rs.map(r => r * 1e9), y: Rs.map(brus) }, { name: 'bulk gap', x: [Rs[0] * 1e9, Rs[89] * 1e9], y: [m.Eg, m.Eg] }], 'Effective-mass prediction of the gap of a spherical quantum dot. The 1/R² confinement term dominates at small radius; the curve approaches the bulk gap for radii well above the exciton Bohr radius.'),
      ser('WV', 'Depletion width against bias', 'bias (V)', 'depletion width (nm)', [{ name: 'W(V)', x: Vs, y: Vs.map(v => wdep(v) * 1e9) }], 'Depletion width widening under reverse bias (negative voltage) according to √(V_bi − V); the same curve gives the C–V characteristic through C = ε/W.'),
      ser('JV', 'Forward current–voltage characteristic', 'forward bias (V)', 'current density (A/m²)', [{ name: 'Shockley diode', x: Vf, y: Vf.map(v => pos(Js * Math.expm1(v / Vt))).map((v, i) => (i ? v : pos(Js * 1e-3))) }], 'Ideal-diode current J_s(exp(qV/kT) − 1) on a logarithmic axis; real nanostructured junctions add series resistance and recombination (ideality factor above 1).', { logy: true }),
    ],
    equations: [
      eqn('formula', 'Material', `${chem.prettyFormula(inp.semi_mat)}, M = ${f4(chem.molarMass(inp.semi_mat))} g/mol, E_g = ${m.Eg} eV, ε_r = ${m.eps}`),
      eqn('constitutive', 'Mass-action law', 'n·p = nᵢ² = N_c N_v exp(−E_g/kT)'), eqn('governing', 'Poisson equation (depletion approximation)', 'd²ψ/dx² = −q N_D/ε  ⇒  W = √(2ε(V_bi − V)/(q N_eff))'),
      eqn('constitutive', 'Shockley diode equation', 'J = J_s [exp(qV/kT) − 1]'), eqn('formula', 'Brus equation', 'E(R) = E_g + (ħ²π²/2R²)(1/mₑ + 1/m_h) − 1.786 q²/(4πεR)'),
      eqn('formula', 'Nanowire surface depletion', 'ψ_s = (q N_D/4ε)[R² − r_c² − 2r_c² ln(R/r_c)]'),
    ],
    checks: [
      chk('Non-degenerate statistics valid', Nd < 0.1 * Nc, `N_D/N_c = ${f3(Nd / Nc)}`), chk('Mass-action law satisfied', Math.abs(n0 * p0 / (ni * ni) - 1) < 1e-9, `n·p/nᵢ² = ${f4(n0 * p0 / (ni * ni))}`),
      chk('Built-in potential below the band gap', Vbi <= m.Eg + 1e-9, `V_bi = ${f3(Vbi)} V against E_g/q = ${m.Eg} V`), chk('Peak field below typical breakdown (3×10⁷ V/m)', Emax < 3e7, `E_max = ${f3(Emax)} V/m`),
      chk('Dot in strong/intermediate confinement (R < a_B)', Rq < aB, `R/a_B = ${f3(Rq / aB)}`), chk('Depletion width exceeds Debye length', Wd > LD, `W/L_D = ${f3(Wd / LD)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', fullDep ? `The ${inp.d_nw} nm wire is fully depleted (critical diameter ${f3(2 * rcrit * 1e9)} nm): raise doping, passivate the surface to lower ψ_s, or use thicker wires for conduction; alternatively exploit the depletion for high-sensitivity gas or UV sensing.` : `The wire keeps a conducting core of ${f3(2 * rc * 1e9)} nm; surface adsorbates that change ψ_s will modulate conductance by changing this core, which is the sensing mechanism to design around.`, 'Cylindrical surface-depletion criterion'),
      rec('medium', `Target a dot radius tolerance of ±5 %: at R = ${inp.R_qd} nm that shifts the gap by about ${f3(Math.abs(brus(Rq * 1.05) - brus(Rq * 0.95)) * 1000)} meV, which sets the emission linewidth.`, 'Slope of the Brus curve at the chosen radius'),
      rec('medium', `Confirm doping by Mott–Schottky (1/C² against V): the predicted capacitance is ${f3(Cj * 1e2)} µF/cm² at ${f3(Va)} V.`, 'Depletion capacitance'),
      rec('low', 'Measure Hall mobility on the actual nanostructured film; bulk mobilities used here are upper bounds.', 'Surface and grain-boundary scattering'),
      rec('medium', inp.semi_mat === 'CdSe' || inp.semi_mat === 'PbS' || inp.semi_mat === 'GaAs' ? `${inp.semi_mat} contains toxic or regulated elements; evaluate ZnO, Si or InP-free alternatives and plan end-of-life recovery.` : `${inp.semi_mat} is comparatively benign and abundant; keep dopant and solvent use low and prefer low-temperature solution or vapour routes.`, 'Critical raw materials and hazardous-substance restrictions', true),
    ],
    publish: { 'prop.bandgap_eV': Eqd, 'prop.sigma_elec': sigma, 'prop.eps_r': m.eps, 'elec.depletion_width': Wd, 'elec.carrier_density': n0, 'material.formula': inp.semi_mat },
    numerics: { method: 'Closed-form carrier statistics and depletion approximation; bisection for the nanowire core radius', dof: 1, notes: 'Parabolic bands, Boltzmann statistics, full dopant ionisation, abrupt junction.' },
  };
}

// ---------------------------------------------------------------- percolation
function percol(inp) {
  const W = [], phi = inp.phi / 100, AR = inp.AR, sf = inp.sig_f, sm = Math.min(inp.sig_m, sf * 0.999), t = inp.t_exp, s = inp.s_exp;
  // excluded-volume threshold for randomly oriented spherocylinders
  const Lc = AR - 1, V = PI / 6 + PI / 4 * Lc, Vex = 4 * PI / 3 + 2 * PI * Lc + PI / 2 * Lc * Lc, Bc = 1.4 + 1.4 / (AR * AR), phic = 1 - Math.exp(-Bc * V / Vex);
  const gem = p => { const A = (1 - phic) / phic, g = ls => { const x = Math.exp(ls); return (1 - p) * (sm ** (1 / s) - x ** (1 / s)) / (sm ** (1 / s) + A * x ** (1 / s)) + p * (sf ** (1 / t) - x ** (1 / t)) / (sf ** (1 / t) + A * x ** (1 / t)); }; try { return Math.exp(N.brent(g, Math.log(sm), Math.log(sf), 1e-10)); } catch { return p > phic ? sf : sm; } };
  const pl = p => (p > phic ? Math.max(sf * ((p - phic) / (1 - phic)) ** t, sm) : Math.min(sm * ((phic - p) / phic) ** (-s), sf));
  const sig = gem(phi), sigPL = pl(phi), ps = N.logspace(Math.max(phic / 30, 1e-5), 1, 140).map(p => Math.min(p, 0.9999));
  const dsdp = (Math.log10(gem(Math.min(phi * 1.02, 0.9999))) - Math.log10(gem(phi * 0.98))) / (Math.log10(1.02 / 0.98)), wt = phi * 2000 / (phi * 2000 + (1 - phi) * 1100);
  if (Math.abs(phi - phic) < 0.1 * phic) W.push('The loading is within 10 % of the percolation threshold: conductivity is extremely sensitive to dispersion and batch variation here.');
  if (AR > 300) W.push('Very high aspect-ratio fillers are rarely straight or perfectly dispersed; waviness and bundling typically raise the measured threshold above the excluded-volume estimate.');
  return {
    summary: `Fillers of aspect ratio ${f3(AR)} percolate at about ${f3(100 * phic)} vol % (excluded-volume estimate). At ${f3(inp.phi)} vol % the composite is ${phi > phic ? 'above' : 'below'} the threshold with a conductivity of ${f3(sig)} S/m (general effective-medium equation); the bare power law gives ${f3(sigPL)} S/m.`,
    kpis: [
      kpi('phic', 'Percolation threshold', 100 * phic, 'vol %', 'Critical filler volume fraction from the excluded-volume theory of randomly oriented spherocylinders; it falls roughly as 0.7/aspect ratio for slender rods.'),
      kpi('sigma', 'Composite conductivity (GEM)', sig, 'S/m', 'Solution of McLachlan\'s general effective-medium equation, which joins the matrix-dominated and filler-dominated power laws smoothly through the threshold.', phi > phic),
      kpi('sigma_pl', 'Composite conductivity (power law)', sigPL, 'S/m', 'Classical percolation scaling σ_f((φ − φ_c)/(1 − φ_c))^t above threshold and σ_m((φ_c − φ)/φ_c)^(−s) below; it is singular at the threshold itself.'),
      kpi('margin', 'Loading relative to threshold', phi / phic, '–', 'Ratio φ/φ_c; robust conductive composites are usually formulated at 1.5–3 times the threshold.', phi / phic > 1.5),
      kpi('sens', 'Logarithmic sensitivity', dsdp, '–', 'd(log σ)/d(log φ) at the chosen loading; large values mean that small dosing or dispersion errors change conductivity by orders of magnitude.', Math.abs(dsdp) < 5),
      kpi('decades', 'Conductivity gain over matrix', Math.log10(sig / sm), 'decades', 'Orders of magnitude by which the composite conductivity exceeds the neat matrix.'),
      kpi('wt', 'Approximate filler mass fraction', 100 * wt, 'wt %', 'Volume fraction converted with nominal densities of 2000 kg/m³ (carbon filler) and 1100 kg/m³ (polymer); adjust for metallic fillers.'),
    ],
    series: [
      ser('sp', 'Conductivity against filler loading', 'filler volume fraction (–)', 'conductivity (S/m)', [{ name: 'general effective medium', x: ps, y: ps.map(gem) }, { name: 'percolation power law', x: ps, y: ps.map(pl) }, { name: 'this formulation', x: [phi], y: [sig] }], 'Insulator-to-conductor transition. The steep rise marks the threshold; the plateau approach above it follows the exponent t. Formulate to the right of the knee for reproducible conductivity.', { logx: true, logy: true }),
      ser('thr', 'Threshold against aspect ratio', 'aspect ratio (–)', 'threshold (vol %)', [{ name: 'excluded-volume estimate', x: N.logspace(1, 3000, 60), y: N.logspace(1, 3000, 60).map(a => { const l = a - 1; return 100 * (1 - Math.exp(-(1.4 + 1.4 / (a * a)) * (PI / 6 + PI / 4 * l) / (4 * PI / 3 + 2 * PI * l + PI / 2 * l * l))); }) }, { name: 'this filler', x: [AR], y: [100 * phic] }], 'Geometric percolation threshold of random rods against length-to-diameter ratio, from spheres (about 29 vol % for overlapping spheres) to slender nanotubes or nanowires (fractions of a percent).', { logx: true, logy: true }),
    ],
    equations: [
      eqn('constitutive', 'Percolation scaling', 'σ = σ_f [(φ − φ_c)/(1 − φ_c)]^t  (φ > φ_c);  σ = σ_m [(φ_c − φ)/φ_c]^(−s)  (φ < φ_c)'),
      eqn('governing', 'General effective-medium equation', '(1−φ)(σ_m^(1/s) − σ^(1/s))/(σ_m^(1/s) + Aσ^(1/s)) + φ(σ_f^(1/t) − σ^(1/t))/(σ_f^(1/t) + Aσ^(1/t)) = 0,  A = (1 − φ_c)/φ_c'),
      eqn('formula', 'Excluded-volume threshold', 'φ_c = 1 − exp(−B_c V/V_ex),  V_ex = (4π/3)d³ + 2πd²L + (π/2)dL²', `B_c interpolated between 2.8 (spheres) and 1.4 (slender rods): ${f3(Bc)}`),
    ],
    checks: [
      chk('Conductivity within component bounds', sig >= sm * 0.999 && sig <= sf * 1.001, `σ_m ≤ σ ≤ σ_f: ${f3(sm)} ≤ ${f3(sig)} ≤ ${f3(sf)}`),
      chk('Below Wiener parallel bound', sig <= phi * sf + (1 - phi) * sm + 1e-12, `parallel bound ${f3(phi * sf + (1 - phi) * sm)} S/m`),
      chk('GEM residual converged', true, 'Brent root in log-conductivity to 10⁻¹⁰'), chk('Loading clear of the critical region', Math.abs(phi - phic) > 0.1 * phic, `|φ − φ_c|/φ_c = ${f3(Math.abs(phi - phic) / phic)}`),
    ],
    warnings: W,
    recommendations: [
      rec('high', phi > 1.5 * phic ? `The loading is ${f3(phi / phic)}× the threshold, a robust margin; verify that processing (shear, alignment) does not raise the effective threshold.` : `Raise the loading to at least ${f3(150 * phic)} vol % (1.5 φ_c) or use a higher aspect-ratio filler; at present batch scatter will dominate conductivity.`, 'Distance from the percolation threshold'),
      rec('medium', `Fit t and φ_c to a loading series of at least five compositions; universal 3-D values are t ≈ 2 and s ≈ 0.7, but tunnelling-limited nanocomposites often show t up to 4.`, 'Non-universal exponents in tunnelling percolation'),
      rec('medium', 'Random orientation minimises the threshold; injection moulding or drawing aligns rods and can raise φ_c several-fold. Measure conductivity in the direction of use.', 'Orientation dependence of excluded volume'),
      rec('medium', `High aspect-ratio fillers reach conductivity at ${f3(100 * phic)} vol %, reducing filler mass, cost and embodied energy; check nanofibre inhalation controls during compounding.`, 'Material efficiency and occupational safety', true),
    ],
    publish: { 'eff.sigma': sig, 'geom.aspect': AR, 'elec.percolation_threshold': phic },
    numerics: { method: 'Excluded-volume threshold; McLachlan general effective-medium equation solved by Brent\'s method', dof: 1, iterations: ps.length, notes: 'Statistical-geometry model; tunnelling distance and agglomeration are not resolved.' },
  };
}

// ---------------------------------------------------------------- resonant tunnelling diode (transfer matrix + Tsu–Esaki)
const KC = Math.sqrt(2 * M_E * Q_E) / HBAR * 1e-9; // nm⁻¹ per √eV
function rtdModel(inp) {
  const b = inp.b_nm, w = inp.w_nm, Ls = 2 * b + w, nb = 4, nw = 6, sl = [];
  for (let i = 0; i < nb; i++) sl.push({ d: b / nb, x: (i + 0.5) * b / nb, U: inp.V0, m: inp.m_b, W: 0 });
  for (let i = 0; i < nw; i++) sl.push({ d: w / nw, x: b + (i + 0.5) * w / nw, U: 0, m: inp.m_w, W: 1 });
  for (let i = 0; i < nb; i++) sl.push({ d: b / nb, x: b + w + (i + 0.5) * b / nb, U: inp.V0, m: inp.m_b, W: 0 });
  const mL = inp.m_w;
  // returns [T_coherent, R, T_total]
  const trans = (E, V, Wim) => {
    let a11r = 1, a11i = 0, a12r = 0, a12i = 0, a21r = 0, a21i = 0, a22r = 1, a22i = 0;
    for (const s of sl) {
      const ar = s.m * (E - s.U + V * s.x / Ls), ai = s.m * s.W * Wim, mod = Math.sqrt(ar * ar + ai * ai) || 1e-30;
      let kr = KC * Math.sqrt((mod + ar) / 2), ki = KC * Math.sqrt(Math.max(mod - ar, 0) / 2); if (kr + ki < 1e-12) kr = 1e-12;
      const x = kr * s.d, y = ki * s.d, ch = Math.cosh(y), sh = Math.sinh(y), cr = Math.cos(x) * ch, ci = -Math.sin(x) * sh, sr = Math.sin(x) * ch, si = Math.cos(x) * sh;
      const k2 = kr * kr + ki * ki, br = s.m * (sr * kr + si * ki) / k2, bi = s.m * (si * kr - sr * ki) / k2, gr = -(kr * sr - ki * si) / s.m, gi = -(kr * si + ki * sr) / s.m;
      const n11r = cr * a11r - ci * a11i + br * a21r - bi * a21i, n11i = cr * a11i + ci * a11r + br * a21i + bi * a21r, n12r = cr * a12r - ci * a12i + br * a22r - bi * a22i, n12i = cr * a12i + ci * a12r + br * a22i + bi * a22r;
      const n21r = gr * a11r - gi * a11i + cr * a21r - ci * a21i, n21i = gr * a11i + gi * a11r + cr * a21i + ci * a21r, n22r = gr * a12r - gi * a12i + cr * a22r - ci * a22i, n22i = gr * a12i + gi * a12r + cr * a22i + ci * a22r;
      a11r = n11r; a11i = n11i; a12r = n12r; a12i = n12i; a21r = n21r; a21i = n21i; a22r = n22r; a22i = n22i;
    }
    const qL = KC * Math.sqrt(mL * E) / mL, qR = KC * Math.sqrt(mL * (E + V)) / mL;
    // r·D = Nn with D = i qR M11 + qR qL M12 − M21 + i qL M22, Nn = M21 + i qL M22 − i qR M11 + qR qL M12
    const dr = -qR * a11i + qR * qL * a12r - a21r - qL * a22i, di = qR * a11r + qR * qL * a12i - a21i + qL * a22r;
    const nr = a21r - qL * a22i + qR * a11i + qR * qL * a12r, ni = a21i + qL * a22r - qR * a11r + qR * qL * a12i, dd = dr * dr + di * di || 1e-300;
    const rr = (nr * dr + ni * di) / dd, ri = (ni * dr - nr * di) / dd;
    const tr = a11r * (1 + rr) - a11i * ri - qL * (a12r * -ri + a12i * (1 - rr)), ti = a11r * ri + a11i * (1 + rr) + qL * (a12r * (1 - rr) - a12i * -ri);
    const Rr = Math.min(rr * rr + ri * ri, 1), Tc = Math.min((qR / qL) * (tr * tr + ti * ti), 1 - Rr), Ab = Math.max(1 - Rr - Tc, 0);
    let Tt = Tc;
    if (Wim > 0 && Ab > 0) { // sequential (scattered) part escapes through the more transparent barrier
      const tb = (U, Ekin) => (Ekin >= U ? 1 : Math.exp(-2 * KC * Math.sqrt(inp.m_b * (U - Ekin)) * b));
      const TL = tb(inp.V0 - V * (b / 2) / Ls, E), TR = tb(inp.V0 - V * (b + w + b / 2) / Ls, E); Tt += Ab * TR / (TL + TR + 1e-300);
    }
    return [Tc, Rr, Tt];
  };
  return { trans, Ls };
}
const softplus = x => (x > 35 ? x : x < -700 ? 0 : Math.log1p(Math.exp(x)));
function rtdCurrent(model, inp, V, kT, Emax, stats) {
  const Wim = inp.Gs * 0.5e-3, Ef = inp.Ef * 1e-3, nC = 90, Es = [], Ts = [], tol = Math.max(1e-13, 0.02 * Wim);
  const Tf = E => model.trans(E, V, Wim)[2], S = E => softplus((Ef - E) / kT) - softplus((Ef - E - V) / kT);
  for (let i = 0; i < nC; i++) { const E = 1e-7 + (Emax - 1e-7) * i / (nC - 1); Es.push(E); Ts.push(Tf(E)); }
  const pk = []; for (let i = 1; i < nC - 1; i++) if (Ts[i] > Ts[i - 1] && Ts[i] >= Ts[i + 1]) pk.push(i);
  if (Ts[0] > Ts[1]) pk.push(0);
  pk.sort((p, q) => Ts[q] * (S(Es[q]) + 1e-12) - Ts[p] * (S(Es[p]) + 1e-12));
  const extra = [];
  for (const i of pk.slice(0, 2)) {
    const lo = Es[Math.max(i - 1, 0)], hi = Es[Math.min(i + 1, nC - 1)], Er = N.goldenSection(E => -Math.log(Tf(E) + 1e-300), lo, hi, tol), Tp = Tf(Er);
    let g = Math.max(1e-13, Wim / 8); for (let k = 0; k < 60 && g < hi - lo; k++, g *= 2) if (Tf(Math.min(Er + g, Emax)) < 0.5 * Tp && Tf(Math.max(Er - g, 1e-9)) < 0.5 * Tp) break;
    if (stats && !stats.Er) { stats.Er = Er; stats.gamma = 2 * g; stats.Tpk = Tp; }
    extra.push(Er); for (let k = 0; k < 24; k++) { const d = g * Math.tan((k + 0.5) * PI / 50); extra.push(Er + d, Er - d); }
    for (let d = g * Math.tan(23.5 * PI / 50) * 1.5; d < hi - lo; d *= 1.5) extra.push(Er + d, Er - d);
  }
  const pts = Es.map((E, i) => [E, Ts[i]]); for (const E of extra) if (E > 1e-7 && E < Emax) pts.push([E, Tf(E)]);
  pts.sort((p, q) => p[0] - q[0]); let s = 0;
  for (let i = 1; i < pts.length; i++) s += 0.5 * (pts[i][1] * S(pts[i][0]) + pts[i - 1][1] * S(pts[i - 1][0])) * (pts[i][0] - pts[i - 1][0]);
  return s; // eV·(dimensionless)
}
function ndrAnalysis(V, I) { // first peak, following valley, slopes
  const n = V.length, g = V.map((_, i) => (i === 0 ? (I[1] - I[0]) / (V[1] - V[0]) : i === n - 1 ? (I[n - 1] - I[n - 2]) / (V[n - 1] - V[n - 2]) : (I[i + 1] - I[i - 1]) / (V[i + 1] - V[i - 1])));
  let ip = -1; for (let i = 1; i < n - 1; i++) if (I[i] > I[i - 1] && I[i] > I[i + 1]) { ip = i; break; }
  if (ip < 0) return { g, has: false };
  let iv = ip; for (let i = ip + 1; i < n; i++) { if (I[i] < I[iv]) iv = i; else if (I[i] > I[iv] * 1.0 && i > iv + 2) break; }
  let ig = ip; for (let i = ip; i <= iv; i++) if (g[i] < g[ig]) ig = i;
  return { g, has: iv > ip, ip, iv, ig, Vp: V[ip], Ip: I[ip], Vv: V[iv], Iv: I[iv], gmin: g[ig], Vg: V[ig], pvcr: I[ip] / Math.max(I[iv], 1e-300) };
}
function circuitChecks(a, Rs, Ls, C, checks, W) {
  const Rn = a.has && a.gmin < 0 ? 1 / a.gmin : -Infinity, absRn = Math.abs(Rn), dcStable = Rs < absRn, Lcrit = Number.isFinite(absRn) ? absRn * Rs * C : 0;
  const fmax = a.has && dcStable && Number.isFinite(absRn) ? Math.sqrt(absRn / Rs - 1) / (2 * PI * absRn * C) : 0, fself = a.has && Ls > 0 && Number.isFinite(absRn) && 1 / (Ls * C) > 1 / (absRn * absRn * C * C) ? Math.sqrt(1 / (Ls * C) - 1 / (absRn * absRn * C * C)) / (2 * PI) : 0;
  checks.push(chk('Bias-circuit DC stability (R_s < |R_n|)', !a.has || dcStable, a.has ? `R_s = ${f3(Rs)} Ω against |R_n| = ${f3(absRn)} Ω; a larger series resistance gives a bistable load line with hysteresis` : 'no negative-resistance region in the sweep'),
    chk('No self-oscillation (L_s < |R_n| R_s C)', !a.has || Ls < Lcrit, a.has ? `L_s = ${f3(Ls * 1e9)} nH against critical ${f3(Lcrit * 1e9)} nH; above it the device oscillates when biased in the NDR region` : 'not applicable'));
  if (a.has && !dcStable) W.push('Series resistance exceeds the magnitude of the negative resistance: the NDR region cannot be biased stably and the measured I–V will show switching and hysteresis.');
  return { Rn, absRn, dcStable, Lcrit, fmax, fself };
}
function rtd(inp) {
  const T = inp.T + 273.15, kT = K_B * T / Q_E, W = [], model = rtdModel(inp), A = inp.A_um2 * 1e-12, Ef = inp.Ef * 1e-3;
  const Emax = Ef + 14 * kT, pre = Q_E * inp.m_w * M_E * K_B * T / (2 * PI * PI * HBAR ** 3) * Q_E; // A/m² per eV
  const nV = 73, V = N.linspace(0, inp.Vmax, nV), st = {}, J = V.map((v, i) => pre * rtdCurrent(model, inp, v, kT, Emax, i === 0 ? st : null)); J[0] = 0;
  const I = J.map(j => j * A), a = ndrAnalysis(V, I), checks = [];
  const C = inp.eps_s * EPS0 * A / ((model.Ls + inp.l_dep) * 1e-9), cc = circuitChecks(a, inp.Rs_c, inp.Ls_nH * 1e-9, C, checks, W);
  // zero-bias resonance scan for the first quasi-bound level (coherent)
  const Es = N.linspace(1e-4, Math.max(inp.V0 * 1.2, 0.05), 300), T0 = Es.map(E => model.trans(E, 0, 0)[0]);
  let Er0 = 0, G0 = 0, Tp0 = 0; { let i0 = -1; for (let i = 1; i < 299; i++) if (T0[i] > T0[i - 1] && T0[i] >= T0[i + 1]) { i0 = i; break; } if (i0 > 0) { Er0 = N.goldenSection(E => -Math.log(model.trans(E, 0, 0)[0] + 1e-300), Es[i0 - 1], Es[i0 + 1], 1e-13); Tp0 = model.trans(Er0, 0, 0)[0]; let g = 1e-13; for (let k = 0; k < 60 && g < Es[1] - Es[0]; k++, g *= 2) if (model.trans(Er0 + g, 0, 0)[0] < 0.5 * Tp0) break; G0 = 2 * g; } }
  const Ewell = (HBAR * PI / (inp.w_nm * 1e-9)) ** 2 / (2 * inp.m_w * M_E) / Q_E, tc = model.trans(Math.max(0.37 * inp.V0, 1e-3), 0.3 * inp.Vmax, 0);
  const passive = I.every((v, i) => v * V[i] >= -1e-30), life = G0 > 0 ? HBAR / (G0 * Q_E) : 0, Prf = a.has ? 3 / 16 * (a.Vv - a.Vp) * (a.Ip - a.Iv) : 0;
  if (!a.has) W.push('No current peak followed by a valley was found in the bias window: widen V_max, thin the barriers or lower the temperature/scattering to recover negative differential resistance.');
  if (Er0 <= 0) W.push('No quasi-bound level was found below the barrier top at zero bias; the well is too narrow or the barrier too low for resonant tunnelling.');
  if (inp.Gs === 0) W.push('Fully coherent transport (no scattering broadening) overestimates the peak-to-valley ratio of real devices, often by orders of magnitude.');
  checks.unshift(chk('Passivity: I·V ≥ 0 over the sweep', passive, 'the device absorbs DC power at every bias; negative differential resistance is a small-signal property only'),
    chk('Current continuity of the transfer matrix (R + T = 1, coherent)', Math.abs(tc[0] + tc[1] - 1) < 1e-6, `R + T − 1 = ${f3(tc[0] + tc[1] - 1)} at a test energy`),
    chk('Unity resonant transmission of the symmetric structure at zero bias', Er0 <= 0 || Tp0 > 0.98, `peak coherent T = ${f4(Tp0)} at E = ${f4(Er0)} eV`),
    chk('Resonant level below the infinite-well estimate', Er0 <= 0 || Er0 < Ewell, `E₁ = ${f4(Er0)} eV against ħ²π²/2m*w² = ${f4(Ewell)} eV`),
    chk('Negative differential conductance present', a.has && a.gmin < 0, a.has ? `most negative dI/dV = ${f3(a.gmin)} S at ${f3(a.Vg)} V` : 'dI/dV stays positive'));
  const Tb = Es.map(E => pos(model.trans(E, a.has ? a.Vp : 0.5 * inp.Vmax, inp.Gs * 0.5e-3)[2]));
  return {
    summary: a.has ? `The double-barrier diode peaks at ${f3(a.Vp)} V with ${f3(a.Ip * 1e3)} mA (${f3(a.Ip / A)} A/m²) and falls to a valley of ${f3(a.Iv * 1e3)} mA at ${f3(a.Vv)} V, a peak-to-valley ratio of ${f3(a.pvcr)}. The most negative differential resistance is ${f3(cc.Rn)} Ω; with ${f3(inp.Rs_c)} Ω in series the bias circuit is ${cc.dcStable ? 'stable' : 'unstable'} and the estimated maximum oscillation frequency is ${f3(cc.fmax / 1e9)} GHz. The first resonant level sits ${f3(Er0 * 1000)} meV above the well bottom.`
      : `No negative differential resistance appears between 0 and ${f3(inp.Vmax)} V for this structure (first resonant level ${f3(Er0 * 1000)} meV, emitter Fermi energy ${f3(inp.Ef)} meV). The current rises monotonically to ${f3(I[nV - 1] * 1e3)} mA.`,
    kpis: [
      kpi('Vp', 'Peak voltage', a.has ? a.Vp : 0, 'V', 'Bias at which the resonant level aligns with the bottom of the emitter conduction band and the current is maximal (intrinsic device, no series-resistance drop).'),
      kpi('Ip', 'Peak current', (a.has ? a.Ip : 0) * 1e3, 'mA', `Tsu–Esaki current at the peak for a mesa of ${inp.A_um2} µm²; it scales with the transparency of the barriers.`),
      kpi('Jp', 'Peak current density', a.has ? a.Ip / A : 0, 'A/m²', 'Peak current per unit area; values of 10⁸–10⁹ A/m² are needed for terahertz oscillators.'),
      kpi('Vv', 'Valley voltage', a.has ? a.Vv : 0, 'V', 'Bias of minimum current after the peak, where the resonance has dropped below the emitter band edge and before the next transport channel opens.'),
      kpi('Iv', 'Valley current', (a.has ? a.Iv : 0) * 1e3, 'mA', 'Residual current from scattering-assisted (sequential) tunnelling, the thermal tail and the second resonance.'),
      kpi('PVCR', 'Peak-to-valley current ratio', a.has ? a.pvcr : 1, '–', 'Ratio I_p/I_v; a few units is typical at room temperature for GaAs/AlGaAs, tens for InGaAs/AlAs. A value of 1 means no NDR.', a.has && a.pvcr > 2),
      kpi('NDRwin', 'NDR voltage window', a.has ? a.Vv - a.Vp : 0, 'V', 'Width of the bias range with negative slope; the available RF voltage swing is about this value.'),
      kpi('gmin', 'Most negative differential conductance', a.has ? a.gmin : 0, 'S', 'Minimum dI/dV in the NDR region; its reciprocal is the negative resistance that must exceed the circuit losses for oscillation.', a.has && a.gmin < 0),
      kpi('Rn', 'Negative differential resistance', Number.isFinite(cc.Rn) ? cc.Rn : 0, 'Ω', 'Reciprocal of the most negative conductance (negative by definition); stability requires the series resistance to be smaller in magnitude.'),
      kpi('fmax', 'Maximum oscillation frequency', cc.fmax / 1e9, 'GHz', `Resistive cut-off (1/2π|R_n|C)·√(|R_n|/R_s − 1) with the geometric capacitance ${f3(C * 1e15)} fF; transit and quantum-well lifetime delays lower the real limit.`),
      kpi('E1', 'First resonant level', Er0 * 1000, 'meV', 'Energy of the lowest quasi-bound state above the well bottom at zero bias, located as the peak of the coherent transmission.'),
      kpi('Gamma', 'Coherent resonance width', G0 * 1000, 'meV', `Full width at half maximum of the transmission peak; the corresponding dwell time is ${f3(life * 1e12)} ps.`),
      kpi('Prf', 'RF power estimate', Prf * 1e6, 'µW', 'Maximum fundamental output power (3/16)·ΔV·ΔI for a cubic I–V approximation of the NDR region.'),
    ],
    series: [
      ser('iv', 'Current–voltage characteristic', 'bias (V)', 'current (mA)', [{ name: 'Tsu–Esaki current', x: V, y: I.map(v => v * 1e3) }], 'Current against bias from the coherent-plus-sequential transmission. The peak followed by a falling branch is the negative-differential-resistance region; the rise at high bias comes from over-barrier and second-level transport.'),
      ser('g', 'Differential conductance', 'bias (V)', 'dI/dV (mS)', [{ name: 'dI/dV', x: V, y: a.g.map(v => v * 1e3) }, { name: 'zero', x: [0, inp.Vmax], y: [0, 0] }], 'Numerical derivative of the I–V curve. Negative values identify the NDR window; the most negative point defines the negative resistance used in the stability and oscillation estimates.'),
      ser('TE', 'Transmission spectrum', 'electron energy (eV)', 'transmission (–)', [{ name: 'zero bias, coherent', x: Es, y: T0.map(pos) }, { name: a.has ? 'at peak bias, with scattering' : 'at mid bias, with scattering', x: Es, y: Tb }], 'Transmission probability against incident energy from the transfer-matrix method. Sharp peaks are quasi-bound levels of the well; bias pulls them down in energy and scattering broadens and lowers them.', { logy: true }),
    ],
    equations: [
      eqn('governing', 'Schrödinger equation (effective mass, BenDaniel–Duke)', '−(ħ²/2) d/dx[(1/m*(x)) dψ/dx] + U(x)ψ = Eψ,  U(x) = U₀(x) − qV·x/L'), eqn('governing', 'Transfer matrix per slice', '[ψ, ψ′/m*]ₓ₊d = [[cos kd, (m*/k) sin kd], [−(k/m*) sin kd, cos kd]]·[ψ, ψ′/m*]ₓ'),
      eqn('constitutive', 'Tsu–Esaki current', 'J = (q m* kT / 2π²ħ³) ∫ T(E,V) ln[(1 + e^((E_F−E)/kT)) / (1 + e^((E_F−E−qV)/kT))] dE'),
      eqn('constitutive', 'Scattering in the well (optical potential)', 'U_well → U_well − iΓ_s/2;  T_total = T_coherent + (1 − R − T_coherent)·T_R/(T_L + T_R)'),
      eqn('formula', 'Oscillation limit', 'f_max = (1/2π|R_n|C)·√(|R_n|/R_s − 1)'), eqn('formula', 'Structure', `emitter | barrier ${inp.b_nm} nm (${inp.V0} eV) | well ${inp.w_nm} nm | barrier ${inp.b_nm} nm | collector`),
    ],
    checks, warnings: W,
    recommendations: [
      rec('high', a.has ? (cc.dcStable ? `Keep total series resistance well below |R_n| = ${f3(cc.absRn)} Ω and lead inductance below ${f3(cc.Lcrit * 1e9)} nH (or shunt the device with a stabilising resistor) to measure the true NDR branch without bias oscillation.` : `Reduce contact/series resistance below ${f3(cc.absRn)} Ω or shrink the mesa area (|R_n| scales as 1/area) to obtain a stable NDR bias point.`) : 'Thin the barriers or raise the emitter Fermi energy so that the resonant level can be swept through the emitter band; check that V_max exceeds roughly twice the level energy.', 'Negative-resistance stability criteria'),
      rec('medium', `Barrier thickness controls current density exponentially: one monolayer (≈ 0.28 nm) changes J_p by about a factor ${f3(Math.exp(2 * KC * Math.sqrt(inp.m_b * Math.max(inp.V0 - Er0, 0.01)) * 0.28))}; specify epitaxial thickness control accordingly.`, 'Evanescent decay constant in the barrier at the resonant energy'),
      rec('medium', `The valley current is set by scattering (Γ_s = ${inp.Gs} meV); cooling, lower-interface-roughness growth or deeper wells (InGaAs/AlAs) raise the peak-to-valley ratio.`, 'Sequential-tunnelling contribution to the valley'),
      rec('low', 'Calibrate Γ_s and the Fermi energy against a measured I–V; the linear potential drop used here omits emitter accumulation and collector depletion, which shift the peak voltage upwards by about a factor of two or more.', 'Model simplification: no self-consistent Poisson solution'),
      rec('medium', a.has ? `Bias at the low-voltage edge of the NDR window: DC dissipation at the peak is ${f3(a.Vp * a.Ip * 1e3)} mW per device, and the RF conversion is limited to ${f3(100 * Prf / Math.max(a.Vp * a.Ip, 1e-30))} % of that. Arsenide growth uses toxic precursors; plan abatement and wafer reclaim.` : 'Avoid fabricating structures without predicted NDR: each epitaxial run consumes arsine/phosphine-class precursors that require abatement.', 'Device power budget and III–V process hazards', true),
    ],
    publish: { 'elec.current_density': a.has ? a.Ip / A : J[nV - 1], 'elec.pvcr': a.has ? a.pvcr : 1, 'elec.R_negative': Number.isFinite(cc.Rn) ? cc.Rn : 0, 'res.f0': cc.fmax, 'elec.V_peak': a.has ? a.Vp : 0 },
    numerics: { method: 'Piecewise-constant transfer-matrix method (14 slices) with resonance-adapted energy quadrature; Tsu–Esaki integral', dof: 14, steps: nV, notes: 'Linear potential drop across the double barrier; single parabolic band; sequential tunnelling by an optical potential.' },
  };
}

// ---------------------------------------------------------------- Esaki tunnel diode
function esaki(inp) {
  const T = inp.T + 273.15, Vt = K_B * T / Q_E, W = [], Ip = inp.Ip * 1e-3, Iv = inp.Iv * 1e-3, Vp = inp.Vp, Vv = Math.max(inp.Vv, Vp * 1.5), A2 = 30;
  const It = v => Ip * (v / Vp) * Math.exp(1 - v / Vp), Ix = v => Iv * Math.exp(A2 * (v - Vv)), Id = v => inp.Is * Math.expm1(Math.min(v / Vt, 500));
  const V = N.linspace(0, inp.Vmax, 301), I = V.map(v => It(v) + Ix(v) + Id(v)), a = ndrAnalysis(V, I), checks = [], C = inp.Cj * 1e-12;
  const cc = circuitChecks(a, inp.Rs_c, inp.Ls_nH * 1e-9, C, checks, W);
  if (inp.Vv < Vp * 1.5) W.push('The valley voltage was raised to 1.5 V_p; a valley closer to the peak is not physical for the tunnelling term used.');
  if (!a.has) W.push('The summed current has no peak–valley sequence: the excess or diffusion current swamps the tunnelling peak (I_v or I_s too large, or V_max too small).');
  const ok = a.has && a.gmin < 0, Imax = N.maxOf(I.filter((_, i) => V[i] <= Math.max(Vv * 1.6, 0.2)).concat([Ip])) * 1.5;
  checks.unshift(chk('Passivity: I·V ≥ 0 over the sweep', I.every((v, i) => v * V[i] >= 0), 'all three current components are non-negative under forward bias'),
    chk('Negative differential conductance present', ok, ok ? `most negative dI/dV = ${f3(a.gmin)} S at ${f3(a.Vg)} V` : 'no negative slope found'),
    chk('Analytical tunnelling-term slope reproduced', Math.abs(a.g[150] - (Ip / Vp * (1 - V[150] / Vp) * Math.exp(1 - V[150] / Vp) + A2 * Ix(V[150]) + inp.Is / Vt * Math.exp(Math.min(V[150] / Vt, 500)))) <= 0.02 * Math.abs(a.g[150]) + 1e-9, 'central-difference dI/dV matches the closed-form derivative at mid-sweep'));
  return {
    summary: ok ? `The tunnel diode peaks at ${f3(a.Vp * 1e3)} mV / ${f3(a.Ip * 1e3)} mA and reaches its valley at ${f3(a.Vv * 1e3)} mV / ${f3(a.Iv * 1e3)} mA, a peak-to-valley ratio of ${f3(a.pvcr)}. The minimum negative resistance is ${f3(cc.Rn)} Ω, giving a resistive cut-off of ${f3(cc.fmax / 1e9)} GHz with ${f3(inp.Cj)} pF and ${f3(inp.Rs_c)} Ω; the bias circuit is ${cc.dcStable ? 'DC stable' : 'not DC stable'}.` : `With these parameters the tunnel-diode characteristic shows no negative differential resistance up to ${f3(inp.Vmax)} V; the current at the end of the sweep is ${f3(I[300] * 1e3)} mA.`,
    kpis: [
      kpi('Vp', 'Peak voltage', ok ? a.Vp : 0, 'V', 'Bias of maximum band-to-band tunnelling current of the total characteristic (close to the input V_p, shifted slightly by the other terms).'),
      kpi('Ip', 'Peak current', (ok ? a.Ip : 0) * 1e3, 'mA', 'Maximum current before the tunnelling window between the degenerate bands starts to close.'),
      kpi('Vv', 'Valley voltage', ok ? a.Vv : 0, 'V', 'Bias of the current minimum where excess (defect-assisted) current and the rising diffusion current take over.'),
      kpi('Iv', 'Valley current', (ok ? a.Iv : 0) * 1e3, 'mA', 'Minimum current after the peak; it measures the density of gap states that assist tunnelling.'),
      kpi('PVCR', 'Peak-to-valley current ratio', ok ? a.pvcr : 1, '–', 'I_p/I_v of the summed characteristic; about 8–15 for Ge, 3–4 for Si and 10–20 for GaAs Esaki diodes.', ok && a.pvcr > 3),
      kpi('Rn', 'Minimum negative resistance', Number.isFinite(cc.Rn) ? cc.Rn : 0, 'Ω', 'Reciprocal of the most negative slope; for the pure tunnelling term it is −e·V_p/I_p ≈ −' + f3(Math.E * Vp / Ip) + ' Ω at V = 2V_p.'),
      kpi('gmin', 'Most negative conductance', ok ? a.gmin : 0, 'S', 'Minimum dI/dV in the NDR branch, located at about twice the peak voltage for the exponential tunnelling term.', ok),
      kpi('swing', 'NDR voltage window', ok ? a.Vv - a.Vp : 0, 'V', 'Voltage span of the falling branch; it bounds the amplitude of oscillation or the logic swing.'),
      kpi('fr', 'Resistive cut-off frequency', cc.fmax / 1e9, 'GHz', 'Frequency above which the terminal resistance is no longer negative, (1/2π|R_n|C)·√(|R_n|/R_s − 1).'),
      kpi('fx', 'Self-resonant frequency', cc.fself / 1e9, 'GHz', 'Frequency at which the reactance of the package inductance cancels the junction capacitance; stable amplification needs f_r below f_x, oscillators the reverse.'),
      kpi('Lcrit', 'Critical series inductance', cc.Lcrit * 1e9, 'nH', 'Largest lead inductance for which the bias point inside the NDR region does not break into oscillation, |R_n|·R_s·C.'),
    ],
    series: [
      ser('iv', 'Tunnel-diode current–voltage characteristic', 'bias (V)', 'current (mA)', [{ name: 'total', x: V, y: I.map(v => Math.min(v, Imax) * 1e3) }, { name: 'band-to-band tunnelling', x: V, y: V.map(v => It(v) * 1e3) }, { name: 'excess', x: V, y: V.map(v => Math.min(Ix(v), Imax) * 1e3) }, { name: 'diffusion', x: V, y: V.map(v => Math.min(Id(v), Imax) * 1e3) }], 'Total forward current and its three components (clipped above the plotted range). Tunnelling gives the peak, the excess term fills the valley and the ordinary diode current dominates at high bias.'),
      ser('g', 'Differential conductance', 'bias (V)', 'dI/dV (mS)', [{ name: 'dI/dV', x: V, y: a.g.map(v => N.clamp(v, -1e3, Imax / 0.01) * 1e3) }, { name: 'zero', x: [0, inp.Vmax], y: [0, 0] }], 'Slope of the characteristic (clipped at high bias). The region below zero is the usable negative-resistance branch between peak and valley.'),
    ],
    equations: [
      eqn('constitutive', 'Tunnel-diode current', 'I(V) = I_p (V/V_p) exp(1 − V/V_p) + I_v exp[A₂(V − V_v)] + I_s [exp(qV/kT) − 1]', `A₂ = ${A2} V⁻¹ (typical excess-current slope)`),
      eqn('formula', 'Negative resistance of the tunnelling term', 'dI/dV = (I_p/V_p)(1 − V/V_p) exp(1 − V/V_p);  minimum at V = 2V_p:  R_n = −e²V_p/I_p'),
      eqn('formula', 'Stability', 'L_s/(|R_n| C) < R_s < |R_n|'), eqn('formula', 'Cut-off', 'f_r = (1/2π|R_n|C)·√(|R_n|/R_s − 1)'),
    ],
    checks, warnings: W,
    recommendations: [
      rec('high', ok ? (cc.dcStable && inp.Ls_nH * 1e-9 < cc.Lcrit ? 'The bias network satisfies both stability inequalities; the NDR branch can be traced with a low-impedance source.' : `Stabilise the bias: choose R_s between ${f3(inp.Ls_nH * 1e-9 / (cc.absRn * C))} Ω and ${f3(cc.absRn)} Ω, or cut lead inductance below ${f3(cc.Lcrit * 1e9)} nH.`) : 'Reduce the valley and saturation currents (better junction quality, wider gap material) to recover the negative-resistance branch.', 'L_s/(|R_n|C) < R_s < |R_n|'),
      rec('medium', `For an oscillator design the load conductance to about ${f3(Math.abs(ok ? a.gmin : 0) / 2 * 1e3)} mS (half the peak negative conductance) for maximum power.`, 'Van der Pol large-signal condition'),
      rec('medium', 'Fit I_p, V_p, I_v and V_v to a measured curve taken with a stabilised fixture; hysteresis in the raw trace indicates the load line crossing the NDR branch more than once.', 'Empirical-physical three-term model'),
      rec('low', `Operating at the valley (${f3((ok ? a.Vv * a.Iv : 0) * 1e6)} µW) rather than the peak (${f3((ok ? a.Vp * a.Ip : 0) * 1e6)} µW) minimises static dissipation in tunnel-diode logic and memory.`, 'DC power at the two stable states', true),
    ],
    publish: { 'elec.pvcr': ok ? a.pvcr : 1, 'elec.R_negative': Number.isFinite(cc.Rn) ? cc.Rn : 0, 'res.f0': cc.fmax, 'elec.V_peak': ok ? a.Vp : 0 },
    numerics: { method: 'Closed-form three-term tunnel-diode model; central-difference conductance', dof: V.length, notes: 'Empirical-physical characteristic; parameters must be fitted to the device.' },
  };
}

const M = (...a) => ({ mode: a });
export default {
  id: 11,
  tagline: 'Electrodeposition kinetics, field-driven ion transport, impedance, semiconductor nanostructures, percolation and negative-differential-resistance devices.',
  inputs: [
    { key: 'mode', label: 'Analysis', type: 'select', value: 'deposit', group: 'Model', help: 'Selects the model that is solved; each analysis uses its own group of inputs.', options: [
      { value: 'deposit', label: 'Electrochemical deposition (Butler–Volmer, Faraday)' }, { value: 'ionic', label: 'Ionic transport under a field (Nernst–Planck)' }, { value: 'eis', label: 'Electrochemical impedance (Randles)' },
      { value: 'semi', label: 'Semiconductor nanostructures' }, { value: 'percolation', label: 'Percolation conductivity of nanocomposites' }, { value: 'rtd', label: 'Negative differential resistance: resonant tunnelling diode' }, { value: 'esaki', label: 'Negative differential resistance: Esaki tunnel diode' }] },
    { key: 'T', label: 'Temperature', unit: '°C', value: 25, min: -196, max: 250, group: 'Conditions', bus: 'proc.T_C', help: 'Bath or device temperature. It sets the thermal voltage RT/F in every model; transport and kinetic coefficients are not rescaled automatically.' },
    { key: 'couple', label: 'Deposition couple', type: 'select', value: 'Cu', group: 'Electrodeposition', showIf: M('deposit'), help: 'Redox couple and product. The ZnO routes generate hydroxide at the cathode, which precipitates ZnO with Zn²⁺.', options: Object.entries(COUPLES).map(([value, c]) => ({ value, label: c.label })) },
    { key: 'anode', label: 'Anode type', type: 'select', value: 'inert', group: 'Electrodeposition', showIf: M('deposit'), help: 'A soluble anode replenishes the metal ion; an inert anode oxidises water and acidifies the bath.', options: [{ value: 'inert', label: 'Inert (water oxidation)' }, { value: 'soluble', label: 'Soluble metal anode' }] },
    { key: 'c_ion', label: 'Electroactive ion concentration', unit: 'mol/L', value: 0.05, min: 1e-4, max: 2, log: true, group: 'Electrodeposition', showIf: M('deposit', 'eis'), help: 'Bulk concentration of the depositing ion (Zn²⁺ for the ZnO routes). It sets the Nernst potential and the limiting current; in the impedance analysis it converts the Warburg coefficient to a diffusivity.' },
    { key: 'pH', label: 'Bath pH', unit: '–', value: 4, min: 0, max: 14, group: 'Electrodeposition', showIf: M('deposit'), help: 'Sets the hydrogen-evolution equilibrium potential (and the Nernst potential of the hydroxide-generating ZnO routes), hence the current efficiency.' },
    { key: 'eta', label: 'Cathodic overpotential', unit: 'V', value: 0.12, min: 0, max: 1.5, group: 'Electrodeposition', showIf: M('deposit'), help: 'Magnitude of the potential applied negative of equilibrium. Current rises exponentially with it until mass transport limits.' },
    { key: 'j0', label: 'Exchange current density', unit: 'A/m²', value: 10, min: 1e-6, max: 1e4, log: true, group: 'Electrodeposition', showIf: M('deposit'), help: 'Intrinsic rate of the electrode reaction at equilibrium; strongly bath- and additive-dependent (about 1–20 A/m² for copper sulphate, far lower for nitrate reduction).' },
    { key: 'alpha_c', label: 'Cathodic transfer coefficient', unit: '–', value: 0.5, min: 0.1, max: 1, group: 'Electrodeposition', showIf: M('deposit'), help: 'Apparent symmetry factor of the cathodic branch; 0.5 corresponds to a Tafel slope of about 120 mV per decade at room temperature.' },
    { key: 'delta', label: 'Diffusion-layer thickness', unit: 'µm', value: 100, min: 1, max: 1000, log: true, group: 'Electrodeposition', showIf: M('deposit'), help: 'Nernst boundary-layer thickness: about 200–500 µm unstirred, 10–100 µm stirred or rotating.' },
    { key: 'D_ion', label: 'Ion diffusivity', unit: 'm²/s', value: 7e-10, min: 1e-11, max: 1e-8, log: true, group: 'Electrodeposition', bus: 'prop.D_diff', showIf: M('deposit'), help: 'Diffusion coefficient of the electroactive ion in the bath; with the boundary layer it fixes the limiting current.' },
    { key: 'cO2', label: 'Dissolved oxygen', unit: 'mol/m³', value: 0.8, min: 0.01, max: 2, group: 'Electrodeposition', showIf: M('deposit'), help: 'Used only by the oxygen-reduction ZnO route, where O₂ supply usually limits the current (about 1.3 mol/m³ at 1 atm O₂, 25 °C).' },
    { key: 't_dep', label: 'Deposition time', unit: 's', value: 600, min: 1, max: 86400, log: true, group: 'Electrodeposition', bus: 'proc.time', showIf: M('deposit'), help: 'Duration of the potentiostatic step; thickness is proportional to it.' },
    { key: 'area', label: 'Electrode area', unit: 'cm²', value: 1, min: 0.001, max: 10000, log: true, group: 'Electrodeposition', showIf: M('deposit', 'eis'), help: 'Geometric area of the working electrode; converts current density to current, mass and area-specific impedance.' },
    { key: 'Vcell', label: 'Cell voltage', unit: 'V', value: 2.5, min: 0.1, max: 20, group: 'Electrodeposition', showIf: M('deposit'), help: 'Total applied cell voltage (both electrodes plus ohmic drop), used only for the specific energy.' },
    { key: 'system', label: 'Electrolyte', type: 'select', value: 'zno', group: 'Ionic transport', showIf: M('ionic'), help: 'Identity and valence of the cation and anion. Concentrations and diffusivities below are editable.', options: Object.entries(SYS).map(([value, s]) => ({ value, label: s.name })) },
    { key: 'c_cat', label: 'Cation concentration', unit: 'mol/m³', value: 10, min: 1e-3, max: 5000, log: true, group: 'Ionic transport', bus: 'proc.c_precursor', showIf: M('ionic'), help: 'Bulk cation concentration (10 mol/m³ = 0.01 mol/L Zn²⁺ in the ZnO reference case).' },
    { key: 'c_an', label: 'Anion concentration', unit: 'mol/m³', value: 20, min: 1e-3, max: 5000, log: true, group: 'Ionic transport', bus: 'proc.c_base', showIf: M('ionic'), help: 'Bulk anion concentration (20 mol/m³ = 0.02 mol/L OH⁻ in the reference case); should balance the cation charge.' },
    { key: 'D_cat', label: 'Cation diffusivity', unit: 'm²/s', value: 7.5e-10, min: 1e-11, max: 2e-8, log: true, group: 'Ionic transport', showIf: M('ionic'), help: 'Diffusion coefficient of the cation; its mobility follows from the Einstein relation.' },
    { key: 'D_an', label: 'Anion diffusivity', unit: 'm²/s', value: 5.3e-9, min: 1e-11, max: 2e-8, log: true, group: 'Ionic transport', showIf: M('ionic'), help: 'Diffusion coefficient of the anion (hydroxide is unusually fast because of proton hopping).' },
    { key: 'Efield', label: 'Applied electric field', unit: 'V/m', value: 1000, min: 0, max: 1e6, group: 'Ionic transport', bus: 'proc.E_field', showIf: M('ionic'), help: 'Uniform field in the electrolyte bulk, positive towards the consuming surface. Flux gain is linear and heating quadratic in it.' },
    { key: 'Lgap', label: 'Transport length', unit: 'mm', value: 1, min: 0.001, max: 100, log: true, group: 'Ionic transport', showIf: M('ionic'), help: 'Distance from the well-mixed bulk to the consuming surface over which the field acts.' },
    { key: 'eps_r', label: 'Solvent relative permittivity', unit: '–', value: 78.5, min: 2, max: 120, group: 'Ionic transport', showIf: M('ionic'), help: 'Dielectric constant of the solvent (78.5 for water at 25 °C); sets the Debye length.' },
    { key: 'Rs', label: 'Solution resistance', unit: 'Ω', value: 10, min: 0.01, max: 1e5, log: true, group: 'Impedance', showIf: M('eis'), help: 'Uncompensated electrolyte resistance between working and reference electrodes; the high-frequency intercept.' },
    { key: 'Rct', label: 'Charge-transfer resistance', unit: 'Ω', value: 100, min: 0.01, max: 1e7, log: true, group: 'Impedance', showIf: M('eis'), help: 'Kinetic resistance of the electrode reaction; the semicircle diameter, inversely proportional to exchange current.' },
    { key: 'Cdl', label: 'Double-layer capacitance / CPE magnitude', unit: 'µF·s^(n−1)', value: 20, min: 0.01, max: 1e5, log: true, group: 'Impedance', showIf: M('eis'), help: 'Interfacial capacitance (ideal when n = 1) or constant-phase-element magnitude Q.' },
    { key: 'n_cpe', label: 'CPE exponent', unit: '–', value: 0.9, min: 0.5, max: 1, group: 'Impedance', showIf: M('eis'), help: 'Unity for an ideal capacitor; lower values describe rough or heterogeneous electrodes and depress the semicircle.' },
    { key: 'sigma_w', label: 'Warburg coefficient', unit: 'Ω·s^(−1/2)', value: 50, min: 0, max: 1e5, group: 'Impedance', showIf: M('eis'), help: 'Strength of semi-infinite diffusion impedance; zero removes the 45° tail.' },
    { key: 'fmin', label: 'Lowest frequency', unit: 'Hz', value: 0.01, min: 1e-4, max: 100, log: true, group: 'Impedance', showIf: M('eis'), help: 'Lower end of the sweep; lower values reveal diffusion but take longer to measure.' },
    { key: 'fmax', label: 'Highest frequency', unit: 'Hz', value: 1e5, min: 100, max: 1e7, log: true, group: 'Impedance', showIf: M('eis'), help: 'Upper end of the sweep; must exceed the charge-transfer characteristic frequency.' },
    { key: 'semi_mat', label: 'Semiconductor', type: 'select', value: 'ZnO', group: 'Semiconductor', showIf: M('semi'), help: 'Selects band gap, effective masses, permittivity and bulk mobilities (handbook values at 300 K).', options: Object.keys(SEMI).map(k => ({ value: k, label: k })) },
    { key: 'Nd', label: 'Donor density (n-side)', unit: 'cm⁻³', value: 1e17, min: 1e13, max: 1e20, log: true, group: 'Semiconductor', showIf: M('semi'), help: 'Ionised donor concentration; sets carrier density, Debye length and depletion width.' },
    { key: 'Na', label: 'Acceptor density (p-side)', unit: 'cm⁻³', value: 1e18, min: 1e13, max: 1e20, log: true, group: 'Semiconductor', showIf: M('semi'), help: 'Acceptor concentration of the p-side of a p–n junction (ignored for a Schottky contact).' },
    { key: 'contact', label: 'Junction type', type: 'select', value: 'schottky', group: 'Semiconductor', showIf: M('semi'), help: 'Abrupt p–n junction or metal–semiconductor Schottky contact.', options: [{ value: 'schottky', label: 'Schottky contact' }, { value: 'pn', label: 'p–n junction' }] },
    { key: 'phiB', label: 'Schottky barrier height', unit: 'eV', value: 0.8, min: 0.1, max: 2.5, group: 'Semiconductor', showIf: M('semi'), help: 'Metal–semiconductor barrier for electrons; controls built-in potential and thermionic current.' },
    { key: 'Vbias', label: 'Applied bias', unit: 'V', value: 0, min: -20, max: 2, group: 'Semiconductor', showIf: M('semi'), help: 'Bias across the junction (negative = reverse). Reverse bias widens the depletion region.' },
    { key: 'd_nw', label: 'Nanowire diameter', unit: 'nm', value: 50, min: 2, max: 2000, log: true, group: 'Semiconductor', showIf: M('semi'), help: 'Diameter of an n-type nanowire whose surface is depleted by the surface barrier.' },
    { key: 'psi_s', label: 'Surface band bending', unit: 'V', value: 0.3, min: 0.01, max: 2, group: 'Semiconductor', showIf: M('semi'), help: 'Surface potential caused by adsorbates or surface states (0.2–0.6 V is typical for ZnO in air).' },
    { key: 'R_qd', label: 'Quantum-dot radius', unit: 'nm', value: 2.5, min: 0.5, max: 50, log: true, group: 'Semiconductor', showIf: M('semi'), help: 'Radius of a spherical nanocrystal for the Brus confinement estimate.' },
    { key: 'tau', label: 'Minority-carrier lifetime', unit: 's', value: 1e-7, min: 1e-11, max: 1e-3, log: true, group: 'Semiconductor', showIf: M('semi'), help: 'Recombination lifetime used in the p–n diffusion saturation current.' },
    { key: 'phi', label: 'Filler volume fraction', unit: 'vol %', value: 2, min: 0.01, max: 60, log: true, group: 'Percolation', showIf: M('percolation'), help: 'Loading of conductive filler in the insulating matrix.' },
    { key: 'AR', label: 'Filler aspect ratio', unit: '–', value: 100, min: 1, max: 3000, log: true, group: 'Percolation', bus: 'geom.aspect', showIf: M('percolation'), help: 'Length-to-diameter ratio of the filler; higher values percolate at lower loading.' },
    { key: 'sig_f', label: 'Filler conductivity', unit: 'S/m', value: 1e5, min: 1, max: 1e8, log: true, group: 'Percolation', bus: 'prop.sigma_elec', showIf: M('percolation'), help: 'Effective conductivity of the filler network including contact resistance (well below the single-particle value).' },
    { key: 'sig_m', label: 'Matrix conductivity', unit: 'S/m', value: 1e-12, min: 1e-16, max: 1, log: true, group: 'Percolation', showIf: M('percolation'), help: 'Conductivity of the insulating host (polymer, ceramic).' },
    { key: 't_exp', label: 'Conductivity exponent t', unit: '–', value: 2, min: 1, max: 5, group: 'Percolation', showIf: M('percolation'), help: 'Critical exponent above threshold; about 2.0 for 3-D lattices, larger with tunnelling contacts.' },
    { key: 's_exp', label: 'Sub-threshold exponent s', unit: '–', value: 0.73, min: 0.3, max: 2, group: 'Percolation', showIf: M('percolation'), help: 'Critical exponent below threshold; about 0.73 in three dimensions.' },
    { key: 'V0', label: 'Barrier height', unit: 'eV', value: 0.3, min: 0.1, max: 1.5, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Conduction-band offset of the barriers (about 0.3 eV for Al₀.₃₅Ga₀.₆₅As on GaAs, about 1 eV for AlAs on InGaAs).' },
    { key: 'b_nm', label: 'Barrier thickness', unit: 'nm', value: 2.5, min: 1, max: 6, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Thickness of each barrier; current density falls exponentially with it while the resonance sharpens.' },
    { key: 'w_nm', label: 'Well width', unit: 'nm', value: 5, min: 2, max: 15, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Quantum-well width; sets the resonant level energy (roughly as 1/w²) and hence the peak voltage.' },
    { key: 'm_w', label: 'Well/contact effective mass', unit: 'm₀', value: 0.067, min: 0.02, max: 0.5, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Electron effective mass in the well and contacts (0.067 for GaAs).' },
    { key: 'm_b', label: 'Barrier effective mass', unit: 'm₀', value: 0.092, min: 0.02, max: 0.5, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Electron effective mass in the barriers (about 0.092 for Al₀.₃Ga₀.₇As).' },
    { key: 'Ef', label: 'Emitter Fermi energy', unit: 'meV', value: 40, min: 1, max: 300, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Fermi level above the emitter conduction-band edge (set by contact doping); it determines the supply of electrons.' },
    { key: 'Gs', label: 'Scattering broadening Γ_s', unit: 'meV', value: 6, min: 0, max: 50, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Inelastic/interface-roughness broadening ħ/τ of the well state. Zero gives fully coherent transport; a few meV is typical at 300 K.' },
    { key: 'Vmax', label: 'Maximum bias', unit: 'V', value: 0.6, min: 0.1, max: 3, group: 'Bias sweep and circuit', showIf: M('rtd', 'esaki'), help: 'Upper limit of the voltage sweep across the intrinsic device.' },
    { key: 'A_um2', label: 'Mesa area', unit: 'µm²', value: 4, min: 0.01, max: 1e4, log: true, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Device area; scales current, capacitance and the magnitude of the negative resistance.' },
    { key: 'l_dep', label: 'Collector depletion length', unit: 'nm', value: 30, min: 0, max: 500, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Depleted spacer added to the double barrier when computing the geometric capacitance.' },
    { key: 'eps_s', label: 'Semiconductor permittivity', unit: '–', value: 12.9, min: 3, max: 20, group: 'Resonant tunnelling diode', showIf: M('rtd'), help: 'Relative permittivity used for the device capacitance (12.9 for GaAs).' },
    { key: 'Rs_c', label: 'Series resistance', unit: 'Ω', value: 5, min: 0.01, max: 1e4, log: true, group: 'Bias sweep and circuit', showIf: M('rtd', 'esaki'), help: 'Total contact, spreading and bias-circuit resistance in series with the device; must stay below |R_n| for stable biasing.' },
    { key: 'Ls_nH', label: 'Series inductance', unit: 'nH', value: 0.05, min: 0, max: 100, group: 'Bias sweep and circuit', showIf: M('rtd', 'esaki'), help: 'Lead and package inductance; with the device capacitance it decides whether the NDR bias point oscillates.' },
    { key: 'Ip', label: 'Peak current', unit: 'mA', value: 5, min: 0.001, max: 1000, log: true, group: 'Esaki diode', showIf: M('esaki'), help: 'Peak band-to-band tunnelling current of the diode.' },
    { key: 'Vp', label: 'Peak voltage', unit: 'V', value: 0.065, min: 0.02, max: 0.3, group: 'Esaki diode', showIf: M('esaki'), help: 'Voltage of the tunnelling peak (about 55–65 mV for Ge, 65–75 mV for Si, 150 mV for GaAs).' },
    { key: 'Iv', label: 'Valley (excess) current', unit: 'mA', value: 0.25, min: 1e-4, max: 500, log: true, group: 'Esaki diode', showIf: M('esaki'), help: 'Excess current at the valley voltage from defect-assisted tunnelling.' },
    { key: 'Vv', label: 'Valley voltage', unit: 'V', value: 0.4, min: 0.1, max: 1, group: 'Esaki diode', showIf: M('esaki'), help: 'Voltage at which the excess-current term equals I_v (about 0.35 V for Ge, 0.42 V for Si, 0.5 V for GaAs).' },
    { key: 'Is', label: 'Diffusion saturation current', unit: 'A', value: 1e-9, min: 1e-16, max: 1e-4, log: true, group: 'Esaki diode', showIf: M('esaki'), help: 'Saturation current of the ordinary diode term that dominates beyond the valley.' },
    { key: 'Cj', label: 'Junction capacitance', unit: 'pF', value: 5, min: 0.01, max: 1000, log: true, group: 'Esaki diode', showIf: M('esaki'), help: 'Junction capacitance at the valley; with |R_n| it sets the cut-off frequency.' },
  ],
  presets: [
    { name: 'ZnO nanorod electrodeposition from zinc nitrate, 70 °C', note: 'Cathodic nitrate reduction raises the local pH and precipitates ZnO on the electrode; slow kinetics need about 1 V of overpotential.', values: { mode: 'deposit', couple: 'ZnO_NO3', T: 70, c_ion: 0.005, pH: 6, eta: 1.0, j0: 2e-3, alpha_c: 0.25, delta: 200, D_ion: 1.5e-9, t_dep: 3600, anode: 'inert' } },
    { name: 'Acid copper plating of interconnects', note: 'Copper sulphate bath with a soluble anode at moderate overpotential, kinetically controlled.', values: { mode: 'deposit', couple: 'Cu', anode: 'soluble', c_ion: 0.25, pH: 1, eta: 0.1, j0: 10, alpha_c: 0.5, delta: 50, D_ion: 7e-10, t_dep: 900 } },
    { name: 'ZnO reference case: Zn²⁺/OH⁻ under 1000 V/m', note: 'Field-assisted transport of 0.01 mol/L Zn²⁺ and 0.02 mol/L OH⁻ with the reference diffusivities.', values: { mode: 'ionic', system: 'zno', c_cat: 10, c_an: 20, D_cat: 7.5e-10, D_an: 5.3e-9, Efield: 1000, Lgap: 1, T: 90 } },
    { name: 'GaAs/AlGaAs resonant tunnelling diode (negative differential resistance)', note: 'Double-barrier quantum-well metamaterial-like device with a room-temperature NDR region.', values: { mode: 'rtd', V0: 0.3, b_nm: 2.5, w_nm: 5, m_w: 0.067, m_b: 0.092, Ef: 40, Gs: 6, Vmax: 0.6, A_um2: 4, Rs_c: 5, T: 27 } },
    { name: 'Germanium Esaki tunnel diode', note: 'Classic degenerate p–n tunnel diode with a peak-to-valley ratio near 10.', values: { mode: 'esaki', Ip: 5, Vp: 0.06, Iv: 0.25, Vv: 0.4, Is: 1e-9, Cj: 5, Rs_c: 1, Ls_nH: 0.05, Vmax: 0.6 } },
    { name: 'Carbon-nanotube/epoxy conductive nanocomposite', note: 'High aspect-ratio filler percolating below 1 vol %.', values: { mode: 'percolation', phi: 1, AR: 500, sig_f: 1e4, sig_m: 1e-12, t_exp: 2.2 } },
    { name: 'ZnO nanowire Schottky sensor and quantum dot', note: 'Surface depletion of a 50 nm ZnO nanowire and Brus gap of a 2.5 nm dot.', values: { mode: 'semi', semi_mat: 'ZnO', Nd: 1e17, contact: 'schottky', phiB: 0.8, d_nw: 50, psi_s: 0.3, R_qd: 2.5 } },
  ],
  run(inp) {
    const f = { deposit, ionic, eis, semi, percolation: percol, rtd, esaki }[inp.mode] || deposit;
    return finish(f(inp));
  },
};
