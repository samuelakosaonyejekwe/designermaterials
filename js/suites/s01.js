// Suite 1 — Material Requirements and Functional Design.
// Requirement screening of the materials library with Ashby performance indices, weighted-sum and TOPSIS ranking,
// Pareto sets and weight sensitivity; tolerance stack-up, safety margins and reliability allocation; Buckingham Pi
// dimensional analysis; analytic hierarchy process; material specification sheet with requirement traceability.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import { MATERIALS, findMaterial, findClass } from '../../data/materials.js';

const F = N.fmt;
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const num = v => typeof v === 'number' && Number.isFinite(v) ? v : null;
const pct = v => F(100 * v, 3) + ' %';
const SUITE = Object.fromEntries(['Material Requirements and Functional Design', 'Chemical Composition and Material Selection', 'Atomic, Molecular and Crystal Structure', 'Nanomaterial Geometry and Morphology', 'Metamaterial Architecture and Unit-Cell Design', 'Thermodynamics, Chemical Equilibrium and Phase Behaviour', 'Chemical Reaction, Nucleation and Growth', 'Surface, Interface and Colloidal Engineering', 'Fluid Dynamics and Porous-Media Transport', 'Heat Transfer and Thermal Engineering', 'Electrical, Electronic and Electrochemical Engineering', 'Electromagnetic, Photonic and Optical Engineering', 'Magnetic and Spintronic Engineering', 'Structural and Nanomechanical Engineering', 'Acoustic, Phononic and Elastic-Wave Engineering', 'Microstructure, Phase-Field and Defect Evolution', 'Quantum, Atomistic and Molecular Simulation', 'Multiscale Modelling and Homogenisation', 'Large-Scale Component and System Analysis', 'Active Tuning, Programmability and Control', 'Fabrication, Synthesis and Process Engineering', 'Material Behaviour, Damage and Lifetime Prediction', 'Experimental Characterisation and Metrology', 'Verification, Calibration, Validation and Uncertainty', 'Inverse Design, Optimisation and Artificial Intelligence', 'Industrial Scale-Up, Automation and Quality', 'Safety, Standards, Environment and Device Integration', 'Economics, Sustainability, Data and Lifecycle', 'Nanomaterials and Metamaterials Reverse Engineering and Reconstruction'].map((n, i) => [i + 1, n]));

// ---------- library property access (effective values for architected materials, never mixed with constituent values) ----------
const SHAPE = { '0D': 6, '1D': 4, '2D': 2, '3D': 6 };
function getter(m) {
  const p = m.props || {}, e = m.effective || null, arch = !!e && (num(e.density) !== null || num(e.rel_density) !== null);
  const rho = e && num(e.density) !== null ? e.density : e && num(e.rel_density) !== null && num(p.density) !== null ? e.rel_density * p.density : num(p.density);
  return { arch, rho, rhoBulk: num(p.density), E: arch ? num(e.E) : num(p.E), sy: arch ? (num(e.yield) ?? num(e.strength_Pa)) : num(p.yield), nu: e && num(e.nu) !== null ? e.nu : num(p.nu), cte: e && num(e.cte) !== null ? e.cte : num(p.cte),
    k: arch ? null : num(p.k_thermal), sig: arch ? null : num(p.sigma_elec), n: e && num(e.n) !== null ? e.n : num(p.n_refr), Eg: num(p.bandgap_eV), c: num(p.c_sound), Tm: num(p.melting_K), cost: num(p.cost_usd_kg), co2: num(p.co2_kg_kg), pvr: e ? num(e.pvr) : null, shape: SHAPE[m.dimensionality] || 6 };
}
const has = (...v) => v.every(x => x !== null);
const LOAD = {
  tie: [['E/ρ (stiff tie)', 'Pa·m³/kg', g => has(g.E, g.rho) && g.rho > 0 ? g.E / g.rho : null], ['σ_y/ρ (strong tie)', 'Pa·m³/kg', g => has(g.sy, g.rho) && g.rho > 0 ? g.sy / g.rho : null]],
  beam: [['E^(1/2)/ρ (stiff beam)', 'Pa^½·m³/kg', g => has(g.E, g.rho) && g.rho > 0 && g.E > 0 ? Math.sqrt(g.E) / g.rho : null], ['σ_y^(2/3)/ρ (strong beam)', 'Pa^⅔·m³/kg', g => has(g.sy, g.rho) && g.rho > 0 && g.sy > 0 ? Math.pow(g.sy, 2 / 3) / g.rho : null]],
  plate: [['E^(1/3)/ρ (stiff plate)', 'Pa^⅓·m³/kg', g => has(g.E, g.rho) && g.rho > 0 && g.E > 0 ? Math.cbrt(g.E) / g.rho : null], ['σ_y^(1/2)/ρ (strong plate)', 'Pa^½·m³/kg', g => has(g.sy, g.rho) && g.rho > 0 && g.sy > 0 ? Math.sqrt(g.sy) / g.rho : null]],
};
const NANO = ['metal-nanoparticle', 'oxide-nanoparticle', 'quantum-dot', 'magnetic-nanoparticle', 'fullerene', 'nanorod', 'nanowire', 'nanotube', 'nanofibre', 'graphene', 'tmd', 'mxene', 'hbn', 'nanosheet', 'nanoporous', 'nanocrystalline', 'core-shell', 'nanolattice'];
const FUNCS = {
  stiff: { label: 'Lightweight stiff structure', classes: null, idx: inp => LOAD[inp.loading] || LOAD.beam, eq: 'M₁ = E^(1/2)/ρ (beam), E/ρ (tie), E^(1/3)/ρ (plate);  M₂ = σ_y^(2/3)/ρ, σ_y/ρ, σ_y^(1/2)/ρ', next: [14, 18, 5] },
  auxetic: { label: 'Energy-absorbing / auxetic structure', classes: ['auxetic'], idx: () => [['−ν (auxetic strength)', '–', g => g.nu !== null && g.nu < 0 ? -g.nu : null], ['σ_y²/(2Eρ) (specific elastic energy)', 'J/kg', g => has(g.sy, g.E, g.rho) && g.E > 0 && g.rho > 0 ? g.sy * g.sy / (2 * g.E * g.rho) : null]], eq: 'M₁ = −ν (ν < 0);  M₂ = σ_y²/(2Eρ);  isotropic stability bound −1 < ν < 0.5', next: [5, 14, 25] },
  viso: { label: 'Vibration isolation (negative / quasi-zero stiffness)', classes: ['negative-stiffness'], idx: () => [['σ_y/E (elastic strain capacity of the flexure material)', '–', g => has(g.sy, g.E) && g.E > 0 ? g.sy / g.E : null], ['σ_y²/(Eρ) (specific spring energy)', 'J/kg', g => has(g.sy, g.E, g.rho) && g.E > 0 && g.rho > 0 ? g.sy * g.sy / (g.E * g.rho) : null]], eq: 'M₁ = σ_y/E (snap-through flexures stay elastic);  k_total = k₊ + k₋ ≥ 0 for static stability', next: [5, 15, 20] },
  thermal: { label: 'Thermal management', classes: null, idx: () => [['k (heat spreading)', 'W/m/K', g => g.k], ['k/α (thermal distortion resistance)', 'W/m', g => has(g.k, g.cte) && g.cte > 0 ? g.k / g.cte : null]], eq: 'M₁ = k;  M₂ = k/α (minimum thermal distortion at given heat flux)', next: [10, 22] },
  nte: { label: 'Zero / negative thermal expansion', classes: ['negative-thermal-expansion'], idx: () => [['−α (contraction on heating)', '1/K', g => g.cte === null ? null : -g.cte], ['k/|α| (distortion resistance)', 'W/m', g => has(g.k, g.cte) && g.cte !== 0 ? g.k / Math.abs(g.cte) : null]], eq: 'M₁ = −α;  composite target α_eff = Σ αᵢKᵢVᵢ / Σ KᵢVᵢ → 0 (Turner)', next: [10, 18, 5] },
  optical: { label: 'Optical / photonic (metasurface, negative index)', classes: ['negative-index', 'optical-metasurface', 'photonic-crystal', 'hyperbolic', 'chiral', 'epsilon-near-zero', 'tunable', 'nonlinear'], idx: () => [['|n| (index magnitude; effective where reported, else constituent)', '–', g => g.n === null ? null : Math.abs(g.n)], ['−n (negative-index strength)', '–', g => g.n !== null && g.n < 0 ? -g.n : null]], eq: 'M₁ = |n|;  phase coverage of a pillar of height H: Δφ = 2π(n − 1)H/λ;  n < 0 requires Re ε < 0 and Re μ < 0 in a passive medium', next: [12, 5, 25] },
  emi: { label: 'EM shielding / negative ε or μ', classes: ['negative-permittivity', 'negative-permeability', 'negative-index', 'epsilon-near-zero', 'metal-nanoparticle', 'nanowire', 'mxene', 'graphene', 'nanotube', 'constituent'], idx: () => [['σ (electrical conductivity)', 'S/m', g => g.sig], ['σ/ρ (specific conductivity)', 'S·m²/kg', g => has(g.sig, g.rho) && g.rho > 0 ? g.sig / g.rho : null]], eq: 'M₁ = σ;  skin depth δ = √(2/(ωμσ));  absorption loss A ≈ 8.69·t/δ dB;  Drude ε(ω) = 1 − ω_p²/(ω² + iγω)', next: [12, 11, 5] },
  acoustic: { label: 'Acoustic insulation (negative mass / modulus)', classes: ['negative-effective-mass', 'phononic-crystal', 'acoustic-metasurface', 'seismic', 'pentamode', 'bianisotropic-willis'], idx: () => [['Z = ρc (characteristic impedance of the constituent)', 'Pa·s/m', g => has(g.rhoBulk, g.c) ? g.rhoBulk * g.c : null], ['ρ (mass-law surface density per thickness)', 'kg/m³', g => g.rhoBulk]], eq: 'M₁ = Z = ρc;  mass law TL ≈ 20·log₁₀(f·m″) − 47 dB;  ρ_eff(ω) = ρ₀ + ρ_r·ω₀²/(ω₀² − ω²) < 0 just above ω₀', next: [15, 5, 25] },
  sensing: { label: 'Sensing / catalysis (high surface area nanomaterial)', classes: NANO, idx: inp => [['SSA = s/(ρ·d) (geometric specific surface area)', 'm²/kg', g => g.rhoBulk !== null && g.rhoBulk > 0 ? g.shape / (g.rhoBulk * Math.max(inp.dchar, 0.1) * 1e-9) : null], ['1/ρ (specific volume)', 'm³/kg', g => g.rhoBulk !== null && g.rhoBulk > 0 ? 1 / g.rhoBulk : null]], eq: 'M₁ = SSA = s/(ρ·d), s = 6 (sphere), 4 (rod), 2 (sheet of thickness d);  dense, smooth, monodisperse geometry', next: [4, 8, 7] },
  semi: { label: 'Semiconductor / optoelectronic nanomaterial', classes: null, idx: inp => [['Band-gap match 1 − |E_g − E_g*|/E_g*', '–', g => g.Eg !== null && g.Eg > 0 ? Math.max(0, 1 - Math.abs(g.Eg - inp.EgTarget) / inp.EgTarget) : null], ['λ_edge = hc/E_g', 'nm', g => g.Eg !== null && g.Eg > 0 ? 1239.84 / g.Eg : null]], eq: 'M₁ = 1 − |E_g − E_g*|/E_g*;  absorption edge λ = hc/E_g = 1239.84 nm·eV / E_g', next: [11, 17, 12] },
  ndr: { label: 'Electrical switching (negative differential resistance)', classes: ['negative-differential-resistance'], idx: () => [['Peak-to-valley current ratio', '–', g => g.pvr], ['E_g (constituent band gap)', 'eV', g => g.Eg]], eq: 'M₁ = I_peak/I_valley;  NDR region: dI/dV < 0;  stable biasing needs |R_load| < |R_NDR| (voltage-controlled N-type)', next: [11, 17, 20] },
  ndtr: { label: 'Thermal logic (negative differential thermal resistance)', classes: ['negative-differential-thermal-resistance'], idx: () => [['k (constituent thermal conductivity)', 'W/m/K', g => g.k], ['T_m (thermal stability)', 'K', g => g.Tm]], eq: 'NDTR: ∂J/∂(ΔT) < 0 over part of the operating range;  thermal amplification β = |∂J_C/∂J_G| > 1 needs NDTR', next: [10, 17, 16] },
};
const OBJ = { index: ['Performance index M₁', 'max', null], density: ['Density (kg/m³)', 'min', g => g.rho], E: ['Young’s modulus (Pa)', 'max', g => g.E], yield: ['Yield strength (Pa)', 'max', g => g.sy], k: ['Thermal conductivity (W/m/K)', 'max', g => g.k], cte: ['CTE (1/K)', 'min', g => g.cte], cost: ['Cost (USD/kg)', 'min', g => g.cost], co2: ['Embodied CO₂ (kg/kg)', 'min', g => g.co2], Tm: ['Melting point (K)', 'max', g => g.Tm] };

function screen(inp) {
  const fn = FUNCS[inp.func] || FUNCS.stiff, idx = fn.idx(inp);
  let pool = MATERIALS.filter(m => inp.kind === 'all' || m.kind === inp.kind || (m.kind === 'nanostructured-metamaterial' && inp.kind !== 'constituent'));
  if (inp.classOnly && fn.classes) pool = pool.filter(m => fn.classes.includes(m.classId));
  const cons = [];
  const add = (on, id, name, unit, get, sense, lim) => { if (on) cons.push({ id, name, unit, get, sense, lim }); };
  add(inp.rhoMax < 25000, 'C-RHO', 'Density', 'kg/m³', g => g.rho, 'max', inp.rhoMax);
  add(inp.Emin > 0, 'C-E', 'Young’s modulus', 'Pa', g => g.E, 'min', inp.Emin * 1e9);
  add(inp.yMin > 0, 'C-SY', 'Yield strength', 'Pa', g => g.sy, 'min', inp.yMin * 1e6);
  add(inp.nuMax < 0.5, 'C-NU', 'Poisson’s ratio', '–', g => g.nu, 'max', inp.nuMax);
  add(inp.kMin > 0, 'C-K', 'Thermal conductivity', 'W/m/K', g => g.k, 'min', inp.kMin);
  add(inp.sigMin > 0, 'C-SIG', 'Electrical conductivity', 'S/m', g => g.sig, 'min', inp.sigMin);
  add(inp.cteMax < 100, 'C-CTE', 'Thermal expansion coefficient', '1/K', g => g.cte, 'max', inp.cteMax * 1e-6);
  add(inp.EgMin > 0, 'C-EGLO', 'Band gap (lower)', 'eV', g => g.Eg, 'min', inp.EgMin);
  add(inp.EgMax < 10, 'C-EGHI', 'Band gap (upper)', 'eV', g => g.Eg, 'max', inp.EgMax);
  add(true, 'C-TOP', `Service limit ${F(inp.homol, 2)}·T_m`, 'K', g => g.Tm === null ? null : inp.homol * g.Tm, 'min', inp.Top + 273.15);
  add(inp.costMax < 1e6, 'C-COST', 'Material cost', 'USD/kg', g => g.cost, 'max', inp.costMax);
  add(inp.co2Max < 1e4, 'C-CO2', 'Embodied CO₂', 'kg/kg', g => g.co2, 'max', inp.co2Max);
  const margin = (c, v) => (c.sense === 'max' ? c.lim - v : v - c.lim) / Math.max(Math.abs(c.lim), 1e-12);
  const rows = pool.map(m => { const g = getter(m), st = cons.map(c => { const v = c.get(g); return v === null ? { v: null, s: 'n.d.', mg: null } : { v, mg: margin(c, v), s: margin(c, v) >= 0 ? 'pass' : 'FAIL' }; }); return { m, g, st, ok: st.every(s => s.s !== 'FAIL'), nd: st.filter(s => s.s === 'n.d.').length, M: idx.map(d => num(d[2](g))) }; });
  const surv = rows.filter(r => r.ok), elim = rows.filter(r => !r.ok);
  // criteria: performance index, low density, low cost, low CO2, temperature capability
  let w = [inp.wPerf, inp.wRho, inp.wCost, inp.wCo2, inp.wTemp].map(v => Math.max(0, v)), zeroW = false; if (N.sum(w) <= 0) { w = [1, 1, 1, 1, 1]; zeroW = true; }
  const crit = [['Performance index', 'max', r => r.M[0]], ['Density', 'min', r => r.g.rho], ['Cost', 'min', r => r.g.cost], ['Embodied CO₂', 'min', r => r.g.co2], ['Melting point', 'max', r => r.g.Tm]];
  const X = surv.map(r => crit.map(c => num(c[2](r))));
  const Sn = surv.map(() => new Array(crit.length).fill(null)), Tn = surv.map(() => new Array(crit.length).fill(null)), logged = [];
  crit.forEach((c, j) => {
    const v = X.map(x => x[j]).filter(x => x !== null); if (!v.length) { logged.push(false); return; }
    const lo = N.minOf(v), hi = N.maxOf(v), lg = lo > 0 && hi / lo > 100; logged.push(lg); const t = x => lg ? Math.log(x) : x, a = t(lo), b = t(hi), nrm = Math.sqrt(N.sum(v.map(x => x * x))) || 1;
    X.forEach((x, i) => { if (x[j] === null) return; const s = b > a ? (t(x[j]) - a) / (b - a) : 1; Sn[i][j] = c[1] === 'max' ? s : (b > a ? 1 - s : 1); Tn[i][j] = x[j] / nrm; });
  });
  const wsm = ww => Sn.map(s => { let a = 0, b = 0; s.forEach((v, j) => { if (v !== null && ww[j] > 0) { a += ww[j] * v; b += ww[j]; } }); return b > 0 ? a / b : null; });
  const cover = Sn.map(s => { let b = 0; s.forEach((v, j) => { if (v !== null) b += w[j]; }); return b / N.sum(w); });
  const order = sc => surv.map((_, i) => i).filter(i => sc[i] !== null).sort((a, b) => (X[b][0] !== null) - (X[a][0] !== null) || sc[b] - sc[a] || (surv[a].m.id < surv[b].m.id ? -1 : 1));
  const score = wsm(w), ord = order(score);
  // TOPSIS in the subspace of criteria with data
  const best = crit.map((c, j) => { const v = Tn.map(t => t[j]).filter(x => x !== null); return v.length ? [c[1] === 'max' ? N.maxOf(v) : N.minOf(v), c[1] === 'max' ? N.minOf(v) : N.maxOf(v)] : null; });
  const sw = N.sum(w), tops = Tn.map(t => { let dp = 0, dm = 0, k = 0; t.forEach((v, j) => { if (v === null || !(w[j] > 0)) return; const wj = w[j] / sw; dp += (wj * (v - best[j][0])) ** 2; dm += (wj * (v - best[j][1])) ** 2; k++; }); if (!k) return null; dp = Math.sqrt(dp); dm = Math.sqrt(dm); return dp + dm > 0 ? dm / (dp + dm) : 1; });
  const ordT = order(tops.map((v, i) => score[i] === null ? null : v));
  // weight sensitivity ±20 %
  const topN = Math.min(5, ord.length), scen = [];
  w.forEach((wk, k) => { if (!(wk > 0)) return; for (const s of [-0.2, 0.2]) { const ww = [...w]; ww[k] = wk * (1 + s); const o = order(wsm(ww)); let moved = 0; for (let i = 0; i < topN; i++) if (o[i] !== ord[i]) moved++; scen.push({ name: `${crit[k][0]} ${s > 0 ? '+' : '−'}20 %`, o, moved, top: o[0] }); } });
  return { fn, idx, pool, cons, rows, surv, elim, crit, w, zeroW, X, Sn, score, cover, ord, tops, ordT, scen, topN, logged };
}
function spearman(a, b) { const rk = v => { const o = v.map((x, i) => i).sort((i, j) => v[i] - v[j]), r = new Array(v.length); o.forEach((i, k) => r[i] = k); return r; }; return a.length > 2 ? N.correlation(rk(a), rk(b)) : 1; }
function chemEquations(m) {
  const eq = [];
  if (m.formula) { try { const M = chem.molarMass(m.formula), mf = chem.massFractions(m.formula); eq.push({ kind: 'formula', label: `Formula of ${m.name}`, text: `${chem.prettyFormula(m.formula)}:  M = ${F(M, 5)} g/mol;  mass fractions ${Object.entries(mf).map(([e, v]) => `${e} ${F(100 * v, 4)} %`).join(', ')}` }); } catch { /* non-stoichiometric label */ } }
  const s = (m.synthesis || []).find(x => x.reaction);
  if (s) { try { const fix = l => l.map(x => ({ f: x.f, phase: x.phase || '' })), io = chem.ionicEquations(fix(s.reaction.reactants), fix(s.reaction.products)); eq.push({ kind: 'chemical', label: `Library synthesis route: ${s.route}`, text: io.molecular, note: s.notes || undefined }); if (io.netIonic !== io.molecular && io.completeIonic !== io.molecular) { eq.push({ kind: 'ionic', label: 'Complete ionic equation', text: io.completeIonic }); eq.push({ kind: 'ionic', label: 'Net ionic equation', text: io.netIonic }); } } catch { /* reaction given as a scheme that is not a single balanced step */ } }
  return eq;
}
function publishTop(top, pub) {
  if (!top) return pub; const m = top.m, p = m.props || {}, e = m.effective || {};
  pub['material.name'] = m.name; if (m.formula) pub['material.formula'] = m.formula; pub['material.class'] = m.classId; pub['material.id'] = m.id;
  for (const k of ['density', 'E', 'nu', 'yield', 'KIc', 'k_thermal', 'cp', 'cte', 'sigma_elec', 'eps_r', 'n_refr', 'bandgap_eV', 'Ms', 'surface_energy', 'c_sound']) if (num(p[k]) !== null) pub['prop.' + k] = p[k];
  if (num(p.E) !== null && num(p.nu) !== null && p.nu > -1 && p.nu < 0.5) { pub['prop.G'] = p.E / (2 * (1 + p.nu)); pub['prop.K'] = p.E / (3 * (1 - 2 * p.nu)); }
  if (m.formula) { try { pub['prop.molar_mass'] = chem.molarMass(m.formula) / 1000; } catch { /* skip */ } }
  for (const [k, b] of [['E', 'eff.E'], ['nu', 'eff.nu'], ['cte', 'eff.cte'], ['n', 'eff.n'], ['eps', 'eff.eps'], ['mu', 'eff.mu'], ['rel_density', 'geom.rel_density'], ['f0_Hz', 'res.f0']]) if (num(e[k]) !== null) pub[b] = e[k];
  if (top.g.arch && top.g.rho !== null) pub['eff.density'] = top.g.rho;
  return pub;
}

// ---------- mode: requirements → candidate ranking ----------
function runRank(inp) {
  const S = screen(inp), { fn, idx, cons, surv, elim, ord, ordT, score, tops, cover, X, crit, w } = S, W = [];
  const top = ord.length ? surv[ord[0]] : null, name = r => r.m.name;
  if (S.zeroW) W.push('All weights were zero: equal weights were used for the ranking.');
  if (!S.pool.length) W.push('No library entry belongs to the selected kind and function-specific classes: widen the material kind or switch off the class restriction.');
  else if (!surv.length) W.push('Every candidate is eliminated by the hard constraints: the requirement set is infeasible within this library. See the elimination table for the governing constraints.');
  const ndTot = N.sum(surv.map(r => r.nd)), noIdx = surv.filter(r => r.M[0] === null).length;
  if (ndTot) W.push(`${ndTot} constraint evaluations on surviving candidates could not be made because the library holds no value ("n.d." in the compliance matrix). Those candidates are not screened on that constraint and must be checked by measurement.`);
  if (noIdx) W.push(`${noIdx} surviving candidate${noIdx === 1 ? ' has' : 's have'} no data for the primary performance index and ${noIdx === 1 ? 'is' : 'are'} ranked after the others on the remaining criteria only.`);
  if (surv.some(r => !r.g.arch && r.m.kind !== 'nanomaterial' && r.m.kind !== 'constituent')) W.push('For resonant and device-type metamaterials the tabulated density, modulus, conductivity and cost are those of the constituent material; the functional response (ε, μ, n, ρ_eff, NDR…) is a narrow-band effective property reported separately in the library.');
  W.push('Library values are typical handbook figures for bulk or representative material; nanoscale size effects, processing and supplier variation are not included.');
  // Pareto on two chosen objectives
  const ox = OBJ[inp.objX] || OBJ.density, oy = OBJ[inp.objY] || OBJ.index, gv = (o, r) => o[2] ? num(o[2](r.g)) : r.M[0];
  const pp = surv.map((r, i) => ({ i, x: gv(ox, r), y: gv(oy, r) })).filter(p => p.x !== null && p.y !== null);
  const pf = pp.length ? N.paretoFront(pp.map(p => [ox[1] === 'min' ? p.x : -p.x, oy[1] === 'min' ? p.y : -p.y])).map(k => pp[k]).sort((a, b) => a.x - b.x) : [];
  const lx = pp.length > 1 && pp.every(p => p.x > 0) && N.maxOf(pp.map(p => p.x)) / N.minOf(pp.map(p => p.x)) > 100, ly = pp.length > 1 && pp.every(p => p.y > 0) && N.maxOf(pp.map(p => p.y)) / N.minOf(pp.map(p => p.y)) > 100;
  // conflicts
  const conf = [];
  for (let a = 0; a < cons.length; a++) for (let b = a + 1; b < cons.length; b++) { let na = 0, nb = 0, nab = 0, nd = 0; for (const r of S.rows) { const sa = r.st[a].s, sb = r.st[b].s; if (sa === 'n.d.' || sb === 'n.d.') continue; nd++; if (sa === 'pass') na++; if (sb === 'pass') nb++; if (sa === 'pass' && sb === 'pass') nab++; } if (nd >= 2 && na && nb && nab < 0.34 * Math.min(na, nb)) conf.push([`${cons[a].name} ↔ ${cons[b].name}`, na, nb, nab, nab === 0 ? 'conflict: jointly infeasible in the library' : 'tension: few candidates meet both']); }
  const rankable = ord.filter(i => X[i][0] !== null);
  crit.forEach((c, j) => { if (j === 0 || !(w[j] > 0)) return; const k = rankable.filter(i => X[i][j] !== null); if (k.length < 4) return; const rho = spearman(k.map(i => X[i][0]), k.map(i => c[1] === 'max' ? X[i][j] : -X[i][j])); if (rho < -0.3) conf.push([`${crit[0][0]} ↔ ${c[0]} (objective trade-off)`, k.length, k.length, '–', `Spearman ρ = ${F(rho, 3)}: better index goes with worse ${c[0].toLowerCase()}`]); });
  const elimBy = cons.map((c, j) => S.rows.filter(r => r.st[j].s === 'FAIL').length), gov = cons.length && N.maxOf(elimBy) > 0 ? N.argmax(elimBy) : -1;
  const both = ord.filter(i => tops[i] !== null), agree = both.length > 2 ? spearman(both.map(i => score[i]), both.map(i => tops[i])) : 1;
  const topMg = top ? top.st.map(s => s.mg).filter(v => v !== null) : [], minMg = topMg.length ? N.minOf(topMg) : 0;
  const rev = S.scen.filter(s => s.moved > 0).length, topChange = S.scen.filter(s => s.top !== ord[0]).length;
  const show = ord.slice(0, 25), unitM = idx[0][1];
  const recs = [];
  if (top) {
    recs.push({ priority: 'high', text: `Carry ${name(top)} forward as the lead candidate (weighted score ${F(score[ord[0]], 3)}, ${idx[0][0]} = ${F(top.M[0], 4)} ${unitM}) and keep ${ord.slice(1, 3).map(i => name(surv[i])).join(' and ') || 'no further candidate'} as alternative${ord.length > 2 ? 's' : ''}; confirm the governing properties on the actual product form before design freeze.`, basis: 'Weighted-sum ranking of candidates that pass every hard constraint for which data exist.' });
    recs.push({ priority: topChange ? 'high' : 'low', text: topChange ? `The first rank changes in ${topChange} of ${S.scen.length} weight scenarios (±20 %): agree the weights with the stakeholders (e.g. by the analytic hierarchy process mode) before committing, because the selection is weight-sensitive.` : `The first rank is unchanged in all ${S.scen.length} ±20 % weight scenarios, so the choice does not hinge on the exact weights.`, basis: 'One-at-a-time weight perturbation with renormalisation.' });
    if (top.nd || cover[ord[0]] < 0.999) recs.push({ priority: 'medium', text: `Close the data gaps of ${name(top)}: ${top.nd} constraint${top.nd === 1 ? '' : 's'} not screened and ${pct(1 - cover[ord[0]])} of the criterion weight without data. Measure or source these properties (suite 23) rather than assuming compliance.`, basis: 'Missing library values are reported as not screened and excluded from the score.' });
    const hz = (top.m.hazards || [])[0];
    recs.push({ priority: 'medium', text: `${top.g.co2 !== null ? `Embodied CO₂ of the lead candidate is about ${F(top.g.co2, 3)} kg/kg` : 'No embodied-CO₂ figure is held for the lead candidate: obtain one before comparing routes'}${top.g.cost !== null ? ` and its cost about ${F(top.g.cost, 3)} USD/kg` : ''}. ${hz ? 'Handling note from the library: ' + hz + '.' : ''} Include cost and CO₂ as weighted criteria (suite 28) so that the selection does not optimise performance alone.`, basis: 'Library cost, embodied-carbon and hazard records of the top candidate.', sustainability: true });
  } else recs.push({ priority: 'high', text: gov >= 0 ? `Relax the governing constraint "${cons[gov].name}" (it eliminates ${elimBy[gov]} of ${S.pool.length} candidates) or widen the material kind/class; no library entry meets the present set.` : 'Widen the material kind or switch off the class restriction to obtain candidates.', basis: 'Constraint satisfaction over the library.' });
  if (conf.length) recs.push({ priority: 'medium', text: `Resolve the requirement conflict ${conf[0][0]}: ${conf[0][4]}. Decide which limit is negotiable, or consider an architected (lattice/composite) solution that decouples the two properties.`, basis: 'Pairwise joint-feasibility count over the library.' });
  recs.push({ priority: 'low', text: `Continue with ${fn.next.map(i => `suite ${i} (${SUITE[i]})`).join(', ')} to turn the selected candidate into a quantified design.`, basis: 'Function-specific analysis chain of the platform.' });
  if (recs.length < 4) recs.push({ priority: 'low', text: 'Record the weights, limits and library version with the ranking so that the selection can be reproduced and audited.', basis: 'Requirements traceability practice.' });
  if (!recs.some(r => r.sustainability)) recs.push({ priority: 'medium', text: 'When candidates become available, compare their embodied CO₂, critical-raw-material content and hazard records before performance is fine-tuned.', basis: 'Early-stage selection fixes most of the life-cycle impact.', sustainability: true });
  const series = [];
  if (pp.length) series.push({ key: 'pareto', title: `Selection chart: ${oy[0]} vs ${ox[0]}`, xlabel: ox[0], ylabel: oy[0] + (oy[2] ? '' : ` [${unitM}]`), logx: lx, logy: ly, lines: [{ name: 'feasible candidates', x: pp.map(p => p.x), y: pp.map(p => p.y), points: true }, { name: 'Pareto-optimal set', x: pf.map(p => p.x), y: pf.map(p => p.y), points: true }], explain: `Each point is a library candidate that passes the hard constraints and has data for both axes. The Pareto-optimal set (${ox[1] === 'min' ? 'low' : 'high'} ${ox[0]}, ${oy[1] === 'min' ? 'low' : 'high'} ${oy[0]}) contains the candidates that cannot be improved on one axis without losing on the other.` });
  if (show.length) series.push({ key: 'scores', title: 'Ranked candidate scores', xlabel: 'rank (weighted sum)', ylabel: 'score (0–1)', lines: [{ name: 'weighted sum', x: show.map((_, k) => k + 1), y: show.map(i => score[i]) }, { name: 'TOPSIS closeness', x: show.map((_, k) => k + 1), y: show.map(i => tops[i] ?? 0), dash: true }, { name: 'data coverage', x: show.map((_, k) => k + 1), y: show.map(i => cover[i]), dash: true }], explain: 'Weighted-sum score and TOPSIS closeness of the ranked candidates. Where the two curves disagree in order the ranking depends on the aggregation method; the coverage curve shows what fraction of the criterion weight is backed by library data.' });
  if (S.scen.length && S.topN) series.push({ key: 'wsens', title: 'Weight sensitivity of the leading candidates', xlabel: 'scenario number (see table)', ylabel: 'rank', lines: ord.slice(0, S.topN).map(i => ({ name: name(surv[i]).slice(0, 28), x: S.scen.map((_, k) => k + 1), y: S.scen.map(s => s.o.indexOf(i) + 1) })), explain: 'Rank of the baseline top candidates when each weight is raised or lowered by 20 % and the weights renormalised. Flat lines mean a robust ranking; crossings are rank reversals.' });
  if (cons.length) series.push({ key: 'elim', title: 'Candidates eliminated by each constraint', xlabel: 'constraint number (order of the compliance matrix)', ylabel: 'candidates', lines: [{ name: 'eliminated', x: cons.map((_, i) => i + 1), y: elimBy }, { name: 'no data (not screened)', x: cons.map((_, i) => i + 1), y: cons.map((_, j) => S.rows.filter(r => r.st[j].s === 'n.d.').length), dash: true }], explain: 'Number of candidates removed by each active hard constraint and the number that could not be screened on it for lack of library data. The tallest bar is the requirement that costs the most design freedom.' });
  if (series.length < 2) series.push({ key: 'status', title: 'Screening outcome', xlabel: 'stage (1 pool, 2 feasible, 3 eliminated)', ylabel: 'candidates', lines: [{ name: 'candidates', x: [1, 2, 3], y: [S.pool.length, surv.length, elim.length] }], explain: 'Number of library candidates before screening, after the hard constraints and removed by them. An empty feasible set means the requirement set cannot be met within this library.' });
  const eq = [
    { kind: 'governing', label: 'Constraint screening', text: 'feasible set = { i : gₖ(xᵢ) ≤ 0 for every constraint k with data };  missing data → not screened' },
    { kind: 'formula', label: `Performance indices: ${fn.label}`, text: fn.eq },
    { kind: 'formula', label: 'Weighted-sum score', text: 'Sᵢ = Σⱼ wⱼ·sᵢⱼ / Σⱼ wⱼ over criteria with data;  sᵢⱼ = (xᵢⱼ − min)/(max − min) (log scale when the range exceeds 100×), inverted for cost-type criteria' },
    { kind: 'formula', label: 'TOPSIS', text: 'rᵢⱼ = xᵢⱼ/√(Σᵢ xᵢⱼ²);  d± = √(Σⱼ (wⱼ(rᵢⱼ − rⱼ±))²);  Cᵢ = d⁻/(d⁺ + d⁻)' },
    { kind: 'formula', label: 'Design margin', text: 'm = (limit − value)/|limit| for upper limits,  m = (value − limit)/|limit| for lower limits;  m ≥ 0 complies' },
    { kind: 'formula', label: 'Pareto optimality', text: 'x* is Pareto-optimal ⇔ ∄ x : fₖ(x) ≤ fₖ(x*) ∀k and fₗ(x) < fₗ(x*) for some l' },
    ...(top ? chemEquations(top.m) : []),
  ];
  const csvRows = ord.map((i, k) => [k + 1, surv[i].m.id, '"' + name(surv[i]) + '"', surv[i].m.classId, F(score[i], 5), tops[i] === null ? '' : F(tops[i], 5), surv[i].M[0] === null ? '' : surv[i].M[0], F(cover[i], 3)].join(','));
  return {
    summary: top ? `Of ${S.pool.length} library candidates for "${fn.label}", ${surv.length} pass the ${cons.length} active hard constraint${cons.length === 1 ? '' : 's'} and ${elim.length} are eliminated${gov >= 0 ? ` (most by "${cons[gov].name}": ${elimBy[gov]})` : ''}. The leading candidate is ${name(top)} with weighted score ${F(score[ord[0]], 3)} and ${idx[0][0]} = ${F(top.M[0], 4)} ${unitM}; TOPSIS ranks ${ordT.length ? name(surv[ordT[0]]) : '–'} first. ${pf.length} candidate${pf.length === 1 ? ' is' : 's are'} Pareto-optimal on the chosen axes and the top rank changes in ${topChange} of ${S.scen.length} weight scenarios.` : `No candidate out of ${S.pool.length} library entries for "${fn.label}" satisfies the active hard constraints, so no ranking can be given. ${gov >= 0 ? `The governing constraint is "${cons[gov].name}", which eliminates ${elimBy[gov]} candidates.` : ''} Relax the limits or widen the candidate classes.`,
    kpis: [
      kp('top', 'Leading candidate', top ? name(top) : 'none', '', 'Library entry with the highest weighted-sum score among the candidates that pass all hard constraints for which data exist; it is a screening result, not a qualified selection.', top ? true : false),
      kp('n_pool', 'Candidates considered', S.pool.length, '', 'Number of library entries of the chosen kind (and, when the class restriction is on, of the classes belonging to the stated function) before any constraint is applied.'),
      kp('n_feasible', 'Feasible candidates', surv.length, '', 'Candidates not violating any hard constraint. Constraints with missing data are not applied, so some of these are only partly screened (see compliance matrix).', surv.length > 0),
      kp('n_elim', 'Eliminated candidates', elim.length, '', 'Candidates violating at least one hard constraint; the elimination table names the constraint, the library value and the limit for each.'),
      kp('score', 'Weighted-sum score of the leader', top ? score[ord[0]] : 0, '–', 'Normalised score between 0 (worst of the feasible set on every criterion) and 1 (best on every criterion), aggregated with the stated weights over the criteria that have data.', top ? score[ord[0]] > 0.6 : false),
      kp('M1', `${idx[0][0]} of the leader`, top && top.M[0] !== null ? top.M[0] : 0, unitM, 'Primary Ashby-type performance index for the stated function, evaluated from library properties; a higher value means less mass (or better function) for the same duty.'),
      kp('topsis', 'TOPSIS closeness of the leader', top && tops[ord[0]] !== null ? tops[ord[0]] : 0, '–', 'Relative closeness to the ideal solution (1 = ideal, 0 = anti-ideal) computed in the subspace of criteria for which the candidate has data.'),
      kp('agree', 'Rank agreement (weighted sum vs TOPSIS)', agree, '–', 'Spearman rank correlation between the two multi-criteria methods; values near 1 mean the ranking does not depend on the aggregation rule.', agree > 0.7),
      kp('margin', 'Smallest design margin of the leader', minMg, '–', 'Smallest normalised distance of the leading candidate from any active limit, (limit − value)/|limit| in the complying direction; a small value means a nearly active constraint.', top ? minMg > 0.1 : null),
      kp('pareto', 'Pareto-optimal candidates', pf.length, '', `Number of feasible candidates that are non-dominated on the two chosen objectives (${ox[0]} and ${oy[0]}).`),
      kp('reversals', 'Weight scenarios with rank reversal', rev, `of ${S.scen.length}`, `Number of ±20 % single-weight scenarios in which the order of the top ${S.topN} candidates changes; the top rank itself changes in ${topChange}.`, rev === 0),
    ],
    series,
    tables: [
      { key: 'rank', title: 'Ranked candidates', columns: ['Rank', 'Candidate', 'Class', 'Basis', 'Score', 'TOPSIS', idx[0][0], idx[1][0], 'Density (kg/m³)', 'Data coverage', 'TOPSIS rank'], rows: show.length ? show.map((i, k) => { const r = surv[i]; return [k + 1, name(r), (findClass(r.m.classId) || {}).name || r.m.classId, r.g.arch ? 'effective' : r.m.kind === 'nanomaterial' || r.m.kind === 'constituent' ? 'bulk' : 'constituent', F(score[i], 3), tops[i] === null ? 'n.d.' : F(tops[i], 3), r.M[0] === null ? 'n.d.' : F(r.M[0], 4), r.M[1] === null ? 'n.d.' : F(r.M[1], 4), r.g.rho === null ? 'n.d.' : F(r.g.rho, 4), pct(cover[i]), ordT.indexOf(i) + 1]; }) : [['–', 'none', '–', '–', '–', '–', '–', '–', '–', '–', '–']], explain: 'Feasible candidates in order of weighted-sum score (first 25). "Basis" states whether mechanical values are effective properties of the architecture, bulk values, or values of the constituent of a resonant/device metamaterial. "n.d." marks properties absent from the library.' },
      { key: 'compliance', title: 'Compliance matrix', columns: ['Candidate', ...cons.map(c => `${c.id} ${c.sense === 'max' ? '≤' : '≥'} ${F(c.lim, 4)} ${c.unit}`), 'Min margin'], rows: (show.length ? show.map(i => surv[i]) : elim.slice(0, 10)).map(r => { const mg = r.st.map(s => s.mg).filter(v => v !== null); return [name(r), ...r.st.map(s => s.v === null ? 'n.d.' : `${s.s} (${F(s.v, 3)})`), mg.length ? F(N.minOf(mg), 3) : 'n.d.']; }), explain: 'Status of each listed candidate against every active hard constraint with the library value in brackets. "n.d." means no data: the constraint was not applied and compliance is unknown.' },
      { key: 'elim', title: 'Eliminated candidates and governing constraints', columns: ['Candidate', 'Violated constraints (value vs limit)'], rows: elim.length ? elim.slice(0, 40).map(r => [name(r), r.st.map((s, j) => s.s === 'FAIL' ? `${cons[j].name}: ${F(s.v, 3)} ${cons[j].sense === 'max' ? '>' : '<'} ${F(cons[j].lim, 3)} ${cons[j].unit}` : null).filter(Boolean).join('; ')]) : [['none', '–']], explain: 'Every eliminated candidate with the constraint(s) that removed it (first 40). Use it to see which requirement costs the most candidates.' },
      { key: 'conflict', title: 'Conflicting requirements and trade-offs', columns: ['Pair', 'Meet first', 'Meet second', 'Meet both', 'Assessment'], rows: conf.length ? conf : [['none detected', '–', '–', '–', 'no pair of active constraints is jointly restrictive and no strong objective trade-off was found']], explain: 'Pairs of constraints that few or no library candidates satisfy together although each is satisfiable alone, and criteria whose ranks are negatively correlated with the performance index among feasible candidates.' },
      { key: 'pareto', title: 'Pareto-optimal set', columns: ['Candidate', ox[0], oy[0]], rows: pf.length ? pf.map(p => [name(surv[p.i]), F(p.x, 4), F(p.y, 4)]) : [['none', '–', '–']], explain: 'Non-dominated feasible candidates on the two selected objectives, ordered along the first axis. Any choice outside this set is worse on both objectives than some member of it.' },
      { key: 'wsens', title: 'Weight-sensitivity scenarios', columns: ['No.', 'Scenario', 'Top candidate', `Positions changed in top ${S.topN}`], rows: S.scen.length ? S.scen.map((s, k) => [k + 1, s.name, s.o.length ? name(surv[s.o[0]]) : '–', s.moved]) : [['–', 'no weighted criterion', '–', 0]], explain: 'Outcome of the ranking when one weight at a time is changed by ±20 % and all weights are renormalised. A change of the top candidate signals a decision that depends on the weighting.' },
    ],
    equations: eq,
    checks: [
      { name: 'Feasible set not empty', pass: surv.length > 0, detail: `${surv.length} of ${S.pool.length} candidates pass` },
      { name: 'Leader fully screened', pass: !!top && top.nd === 0, detail: top ? `${top.nd} of ${cons.length} constraints lack data for ${name(top)}` : 'no leader' },
      { name: 'Weights normalised', pass: Math.abs(N.sum(w.map(v => v / N.sum(w))) - 1) < 1e-12, detail: `w = [${w.map(v => F(v / N.sum(w), 3)).join(', ')}] for [${crit.map(c => c[0]).join(', ')}]` },
      { name: 'Scores within [0, 1]', pass: score.every(v => v === null || (v >= -1e-12 && v <= 1 + 1e-12)), detail: ord.length ? `range ${F(score[ord.at(-1)], 3)} to ${F(score[ord[0]], 3)}` : 'no scores' },
      { name: 'Ranking reproducible', pass: true, detail: 'deterministic screening; ties broken by library id' },
      { name: 'Top rank robust to ±20 % weights', pass: topChange === 0 && !!top, detail: `${topChange} of ${S.scen.length} scenarios change the first rank` },
      { name: 'Poisson’s ratio of the leader admissible', pass: !top || top.g.nu === null || (top.g.nu >= -1 && top.g.nu <= 0.5) || top.g.arch || !!top.m.effective, detail: top && top.g.nu !== null ? `ν = ${F(top.g.nu, 3)} (isotropic bound −1 < ν < 0.5; anisotropic lattices may lie outside)` : 'no ν value' },
    ],
    warnings: W, recommendations: recs,
    publish: publishTop(top, { 'opt.objective': top ? score[ord[0]] : 0, 'req.n_feasible': surv.length, 'req.function': fn.label, 'req.T_op_C': inp.Top }),
    exports: [{ name: 'candidate_ranking.csv', mime: 'text/csv', content: 'rank,id,name,class,weighted_score,topsis,index_M1,data_coverage\n' + csvRows.join('\n'), description: 'Full ranking of feasible candidates with scores, primary index and data coverage.' }],
    numerics: { method: 'Constraint screening + Ashby performance indices + weighted-sum and TOPSIS multi-criteria ranking; non-dominated sorting for the Pareto set', dof: S.pool.length, steps: cons.length, iterations: 1 + S.scen.length, residual: 0, notes: `${crit.filter((c, j) => S.logged[j]).map(c => c[0]).join(', ') || 'No criterion'} normalised on a log scale. Partial-data rule: criteria without data are left out and weights renormalised per candidate.` },
  };
}

// ---------- mode: requirement allocation and margins ----------
const STACK_DEF = 'contributor,nominal_mm,tolerance_mm,sensitivity\nlattice cell pitch x4,20.000,0.060,1\nstrut thickness,0.500,0.025,2\nnode fillet,0.150,0.020,2\nbase plate,2.000,0.030,1\nbond line,0.050,0.015,1\nthermal growth at ΔT,0.012,0.004,1';
function parseStack(text) {
  const out = [];
  for (const l of String(text ?? '').split(/\r?\n/).map(s => s.trim()).filter(s => s && !s.startsWith('#') && !s.startsWith('//'))) { const c = l.split(/[,;\t]/).map(s => s.trim()), nom = Number(c[1]), tol = Number(c[2]); if (c.length < 3 || !Number.isFinite(nom) || !Number.isFinite(tol)) continue; const s = Number(c[3]); out.push({ name: c[0] || `item ${out.length + 1}`, nom, tol: Math.abs(tol), sens: c.length > 3 && Number.isFinite(s) ? s : 1 }); }
  if (out.length < 2) throw new Error(`Tolerance table: at least 2 rows "name, nominal, tolerance[, sensitivity]" are needed (found ${out.length}). Paste CSV/TSV text with one header line.`);
  return out.slice(0, 200);
}
function runAlloc(inp) {
  const it = parseStack(inp.stack), W = [], ks = Math.max(inp.kSigma, 0.5), r = N.rng(inp.seed ?? 1);
  if (String(inp.stack).trim() === STACK_DEF) W.push('Illustrative demonstration data in the tolerance table — replace with the tolerances of your own drawing or process capability study.');
  const nom = N.sum(it.map(c => c.sens * c.nom)), wc = N.sum(it.map(c => Math.abs(c.sens) * c.tol)), varc = it.map(c => (c.sens * c.tol / ks) ** 2), sig = Math.sqrt(N.sum(varc)), rss = ks * sig, lim = inp.stackLimit;
  const nMc = Math.round(N.clamp(2e6 / it.length, 4000, 40000)), sm = new Array(nMc); /* up to 40 000 assemblies, fewer for long tolerance tables */ let out = 0; for (let k = 0; k < nMc; k++) { let d = 0; for (const c of it) d += c.sens * r.normal(0, c.tol / ks); sm[k] = d; if (Math.abs(d) > lim) out++; }
  const pAn = sig > 0 ? 2 * (1 - N.normCdf(lim / sig)) : 0, cp = sig > 0 ? lim / (3 * sig) : 99, h = N.histogram(sm, 30), dom = N.argmax(varc), vs = N.sum(varc) || 1;
  // strength margin and reliability index
  const cap = inp.cap, dem = Math.max(inp.dem, 1e-300), allow = inp.derate * cap, fos = cap / dem, fosD = allow / dem, mos = allow / (inp.fosReq * dem) - 1, sR = inp.covCap / 100 * cap, sS = inp.covDem / 100 * dem, sd = Math.hypot(sR, sS), beta = sd > 0 ? (cap - dem) / sd : (cap >= dem ? 99 : -99), pf = N.clamp(N.normCdf(-beta), 0, 1);
  // reliability allocation
  let wt = String(inp.relW ?? '').split(/[,;\s]+/).map(Number).filter(v => Number.isFinite(v) && v > 0); const nS = Math.round(inp.nSub); if (wt.length !== nS) { if (wt.length) W.push(`The complexity-weight list has ${wt.length} entries but ${nS} subsystems are declared: equal apportionment is used.`); wt = new Array(nS).fill(1); }
  const Rs = N.clamp(inp.Rsys, 1e-9, 1 - 1e-12), ws = N.sum(wt), Ri = wt.map(v => Math.pow(Rs, v / ws)), th = Math.max(inp.tMission, 1e-9), lam = Ri.map(v => -Math.log(v) / th), Rach = Math.pow(N.clamp(inp.Ritem, 0, 1), nS), prodR = Ri.reduce((a, b) => a * b, 1);
  if (wc > lim && rss <= lim) W.push('The stack meets the limit statistically (RSS) but not in the worst case: this is acceptable only if the contributors are independent, centred and in statistical control.');
  if (beta < 3 && beta > -90) W.push(`Reliability index β = ${F(beta, 3)} is below 3: the capacity and demand distributions overlap appreciably (normal, independent variables assumed).`);
  const nn = [...Array(20).keys()].map(i => i + 1);
  return {
    summary: `The stack of ${it.length} contributors has a nominal of ${F(nom, 5)} mm with ±${F(wc, 3)} mm worst-case and ±${F(rss, 3)} mm RSS (${F(ks, 2)}σ) variation against a limit of ±${F(lim, 3)} mm; ${pct(out / nMc)} of ${nMc} Monte Carlo assemblies fall outside. The factor of safety is ${F(fos, 3)} (${F(fosD, 3)} after derating to ${pct(inp.derate)}), giving a margin of safety of ${F(mos, 3)} against the required ${F(inp.fosReq, 3)} and a reliability index β = ${F(N.clamp(beta, -99, 99), 3)}. A series system of ${nS} subsystems with target ${F(Rs, 6)} requires per-item reliabilities of at least ${F(N.minOf(Ri), 6)}.`,
    kpis: [
      kp('wc', 'Worst-case stack variation', wc, '± mm', 'Arithmetic sum of |sensitivity|·tolerance of all contributors: the variation if every dimension sits at its unfavourable limit simultaneously.', wc <= lim),
      kp('rss', 'RSS stack variation', rss, '± mm', `Root-sum-square combination assuming independent, centred, normally distributed contributors whose tolerance equals ${F(ks, 2)} standard deviations.`, rss <= lim),
      kp('ratio', 'RSS / worst-case ratio', wc > 0 ? rss / wc : 0, '–', 'How much the statistical stack is smaller than the worst case; it tends to 1/√n for n equal contributors and to 1 when one contributor dominates.'),
      kp('p_out', 'Fraction outside the limit (Monte Carlo)', out / nMc, '–', `Share of ${nMc} simulated assemblies whose deviation exceeds the allowable limit; the normal-theory value is ${F(pAn, 3)}.`, out / nMc < 0.0027),
      kp('cp', 'Assembly capability index C_p', Math.min(cp, 99), '–', 'Allowable limit divided by three standard deviations of the stack; 1.33 or more is a common target for a capable assembly process.', cp >= 1.33),
      kp('fos', 'Factor of safety', fos, '–', 'Mean capacity divided by mean demand before derating; it must exceed the required factor of safety with the derating applied.', fosD >= inp.fosReq),
      kp('mos', 'Margin of safety', mos, '–', 'Derated capacity /(required factor of safety × demand) − 1; zero or more means the requirement is met with the stated factor.', mos >= 0),
      kp('beta', 'Reliability index β', N.clamp(beta, -99, 99), '–', 'Distance between mean capacity and mean demand in units of their combined standard deviation (first-order second-moment, normal variables).', beta >= 3),
      kp('pf', 'Probability of failure', pf, '–', 'P(capacity < demand) = Φ(−β) under the normal-independent assumption; it is sensitive to the distribution tails, so treat it as an indicator.', pf < 1e-3),
      kp('Ri', 'Tightest allocated item reliability', N.maxOf(Ri), '–', 'Highest per-subsystem reliability target produced by the apportionment; subsystems with small complexity weight receive the tightest target.'),
      kp('lambda', 'Largest allocated failure rate', N.maxOf(lam), '1/h', 'Constant failure rate that the most complex subsystem may have over the mission time, λ = −ln(Rᵢ)/t.'),
      kp('Rach', 'System reliability at the stated item reliability', Rach, '–', `Series reliability R_item^n if each of the ${nS} subsystems achieved the stated item reliability; compare with the target ${F(Rs, 6)}.`, Rach >= Rs),
    ],
    series: [
      { key: 'hist', title: 'Monte Carlo distribution of the stack deviation', xlabel: 'deviation from nominal (mm)', ylabel: 'count', lines: [{ name: 'assemblies', x: h.x, y: h.y }, { name: '+limit', x: [lim, lim], y: [0, N.maxOf(h.y)], dash: true }, { name: '−limit', x: [-lim, -lim], y: [0, N.maxOf(h.y)], dash: true }], explain: 'Histogram of the simulated assembly deviation with the allowable limits. The area outside the dashed lines is the predicted reject fraction.' },
      { key: 'contrib', title: 'Variance contribution of each contributor', xlabel: 'contributor number', ylabel: 'share of stack variance (%)', lines: [{ name: 'variance share', x: it.map((_, i) => i + 1), y: varc.map(v => 100 * v / vs) }, { name: 'worst-case share', x: it.map((_, i) => i + 1), y: it.map(c => wc > 0 ? 100 * Math.abs(c.sens) * c.tol / wc : 0), dash: true }], explain: 'Share of the total variance (statistical) and of the worst-case sum carried by each contributor. Tightening the largest bar is the efficient way to reduce the stack.' },
      { key: 'rel', title: 'Series reliability versus number of subsystems', xlabel: 'number of subsystems in series', ylabel: 'reliability', lines: [{ name: 'system at stated item reliability', x: nn, y: nn.map(n => Math.pow(N.clamp(inp.Ritem, 0, 1), n)) }, { name: 'required item reliability (equal split)', x: nn, y: nn.map(n => Math.pow(Rs, 1 / n)), dash: true }], explain: 'System reliability falls geometrically with the number of series items, and the per-item requirement rises towards 1. It shows why part count reduction is a reliability measure.' },
    ],
    tables: [
      { key: 'stack', title: 'Tolerance contributors', columns: ['No.', 'Contributor', 'Nominal (mm)', 'Tolerance (± mm)', 'Sensitivity', 'Variance share (%)'], rows: it.map((c, i) => [i + 1, c.name, c.nom, c.tol, c.sens, F(100 * varc[i] / vs, 3)]), explain: 'Parsed contributors with their sensitivity (∂assembly/∂dimension) and share of the statistical stack variance.' },
      { key: 'alloc', title: 'Reliability allocation', columns: ['Subsystem', 'Complexity weight', 'Allocated reliability', 'Allowed failure rate (1/h)', 'MTBF (h)'], rows: Ri.map((v, i) => [i + 1, wt[i], F(v, 7), F(lam[i], 4), lam[i] > 0 ? F(1 / lam[i], 4) : '∞']), explain: 'Apportionment of the system reliability target to series subsystems in proportion to their complexity weights, with the equivalent constant failure rate over the mission time.' },
    ],
    equations: [
      { kind: 'formula', label: 'Worst-case stack-up', text: 'T_wc = Σᵢ |∂y/∂xᵢ|·tᵢ' },
      { kind: 'formula', label: 'Root-sum-square stack-up', text: 'T_rss = k·√(Σᵢ (∂y/∂xᵢ · tᵢ/k)²);  C_p = T_allow/(3σ_y)' },
      { kind: 'formula', label: 'Factor and margin of safety', text: 'FoS = R/S;  MoS = D·R/(FoS_req·S) − 1   (D = derating factor)' },
      { kind: 'formula', label: 'Reliability index', text: 'β = (μ_R − μ_S)/√(σ_R² + σ_S²);  P_f = Φ(−β)' },
      { kind: 'formula', label: 'Series reliability allocation', text: 'R_sys = Πᵢ Rᵢ;  Rᵢ = R_sys^(wᵢ/Σw);  λᵢ = −ln(Rᵢ)/t' },
    ],
    checks: [
      { name: 'RSS ≤ worst case', pass: rss <= wc * (1 + 1e-9) || ks > 3, detail: `RSS ${F(rss, 4)} mm vs worst case ${F(wc, 4)} mm` },
      { name: 'Monte Carlo agrees with normal theory', pass: Math.abs(out / nMc - pAn) <= 3 * Math.sqrt(Math.max(pAn * (1 - pAn), 1e-6) / nMc) + 1e-3, detail: `sampled ${F(out / nMc, 3)} vs analytical ${F(pAn, 3)}` },
      { name: 'Allocated reliabilities reproduce the system target', pass: Math.abs(prodR - Rs) < 1e-9, detail: `Π Rᵢ = ${F(prodR, 8)} vs target ${F(Rs, 8)}` },
      { name: 'Margin of safety non-negative', pass: mos >= 0, detail: `MoS = ${F(mos, 4)} with FoS_req = ${F(inp.fosReq, 3)} and derating ${F(inp.derate, 3)}` },
      { name: 'Stack within limit (statistical)', pass: rss <= lim, detail: `±${F(rss, 4)} mm vs ±${F(lim, 4)} mm` },
    ],
    warnings: W,
    recommendations: [
      { priority: rss > lim ? 'high' : 'medium', text: `Tighten or redesign "${it[dom].name}" first: it carries ${F(100 * varc[dom] / vs, 3)} % of the stack variance. Halving its tolerance reduces the RSS stack to about ±${F(ks * Math.sqrt(Math.max(vs - 0.75 * varc[dom], 0)), 3)} mm.`, basis: 'Variance decomposition of the RSS stack.' },
      { priority: mos < 0 ? 'high' : 'low', text: mos < 0 ? `Raise the capacity to at least ${F(inp.fosReq * dem / Math.max(inp.derate, 1e-9), 4)} Pa or reduce the demand to ${F(allow / inp.fosReq, 4)} Pa to restore a non-negative margin.` : `The margin of safety is ${F(mos, 3)}; mass or cost could be removed until it approaches zero, provided the scatter data are trustworthy.`, basis: 'Margin of safety with derating and the required factor of safety.' },
      { priority: Rach < Rs ? 'high' : 'low', text: Rach < Rs ? `The stated item reliability gives only ${F(Rach, 5)} for ${nS} series items: raise item reliability to ${F(Math.pow(Rs, 1 / nS), 6)}, reduce the part count or add redundancy to the weakest item.` : 'The stated item reliability meets the system target; keep the allocation as the subsystem requirement in the specification.', basis: 'Series reliability block model.' },
      { priority: 'medium', text: 'Loosen tolerances whose variance share is below about 2 %: they do not influence the stack, and relaxed tolerances reduce scrap, rework and inspection energy.', basis: 'Contributors with negligible variance share.', sustainability: true },
    ],
    publish: { 'mech.safety_factor': fos, 'req.margin_of_safety': mos, 'req.stack_rss_mm': rss, 'req.stack_wc_mm': wc, 'req.beta': N.clamp(beta, -99, 99), 'req.R_item': N.maxOf(Ri), 'fab.yield': 1 - out / nMc },
    numerics: { method: 'Worst-case and RSS tolerance analysis with Monte Carlo confirmation; first-order second-moment reliability; series reliability apportionment', dof: it.length, steps: nMc, iterations: nMc, residual: Math.abs(out / nMc - pAn), notes: 'Normal, independent, centred contributors; seeded sampling.' },
  };
}

// ---------- mode: Buckingham Pi ----------
const frac = a => { for (const d of [1, 2, 3, 4, 6]) { const n = Math.round(a * d); if (Math.abs(a * d - n) < 1e-6) return d === 1 ? String(n) : `${n}/${d}`; } return F(a, 3); };
const pw = (s, a) => Math.abs(a) < 1e-9 ? '' : Math.abs(a - 1) < 1e-9 ? s : `${s}^${frac(a).includes('/') || a < 0 ? '(' + frac(a) + ')' : frac(a)}`;
const PIP = {
  reactor: { label: 'Porous reactor synthesis (flow, diffusion, first-order reaction)', vars: i => [['U', 'superficial velocity', 'm/s', [0, 1, -1, 0], i.pU], ['d', 'pore / particle size', 'm', [0, 1, 0, 0], i.pd], ['ρ', 'fluid density', 'kg/m³', [1, -3, 0, 0], i.pRho], ['μ', 'dynamic viscosity', 'Pa·s', [1, -1, -1, 0], i.pMu], ['L', 'bed length', 'm', [0, 1, 0, 0], i.pL], ['D', 'diffusivity', 'm²/s', [0, 2, -1, 0], i.pD], ['k', 'rate constant', '1/s', [0, 0, -1, 0], i.pK]],
    named: [['Re = ρUd/μ', { 'ρ': 1, U: 1, d: 1, 'μ': -1 }, v => v < 1 ? 'creeping (Darcy) flow in the pores' : v < 150 ? 'inertial laminar pore flow (Forchheimer correction needed)' : 'transitional to turbulent pore flow'], ['Sc = μ/(ρD)', { 'μ': 1, 'ρ': -1, D: -1 }, v => v > 100 ? 'liquid-like: momentum diffuses far faster than species' : 'gas-like: comparable momentum and species diffusion'], ['Pe = UL/D', { U: 1, L: 1, D: -1 }, v => v > 100 ? 'advection-dominated axial transport (plug-flow-like)' : v > 1 ? 'mixed advection and diffusion' : 'diffusion-dominated'], ['Da = kL/U', { k: 1, L: 1, U: -1 }, v => v > 10 ? 'reaction fast relative to residence time: high conversion, transport-limited' : v > 0.1 ? 'reaction and residence time comparable' : 'kinetically limited: low conversion per pass'], ['φ² = kd²/D', { k: 1, d: 2, D: -1 }, v => v > 1 ? 'pore-diffusion limitation (Thiele modulus > 1): concentration gradients inside pores' : 'uniform concentration inside pores']] },
  heat: { label: 'Heat-up of a particle, film or substrate', vars: i => [['k', 'thermal conductivity', 'W/m/K', [1, 1, -3, -1], i.pKth], ['L', 'characteristic length', 'm', [0, 1, 0, 0], i.pLh], ['ρ', 'density', 'kg/m³', [1, -3, 0, 0], i.pRhoS], ['c_p', 'specific heat', 'J/kg/K', [0, 2, -2, -1], i.pCp], ['h', 'heat-transfer coefficient', 'W/m²/K', [1, 0, -3, -1], i.pH], ['t', 'time', 's', [0, 0, 1, 0], i.pT]],
    named: [['Bi = hL/k', { h: 1, L: 1, k: -1 }, v => v < 0.1 ? 'lumped: the body is nearly isothermal while heating' : v < 10 ? 'internal and external resistances comparable: solve the conduction problem' : 'surface follows the fluid temperature; internal conduction controls'], ['Fo = kt/(ρc_pL²)', { k: 1, t: 1, 'ρ': -1, c_p: -1, L: -2 }, v => v > 1 ? 'thermal wave has crossed the body (near steady state when Bi is large)' : v > 0.2 ? 'one-term series solution valid' : 'early time: semi-infinite behaviour'], ['Bi·Fo = ht/(ρc_pL)', { h: 1, t: 1, 'ρ': -1, c_p: -1, L: -1 }, v => `lumped response θ = exp(−Bi·Fo) = ${F(Math.exp(-Math.min(v, 700)), 3)} of the initial temperature difference remains`]] },
  resonator: { label: 'Resonator dynamics (local resonator of a metamaterial)', vars: i => { const w = 2 * Math.PI * i.pF, r2 = i.pM * w * w / i.pKs, X = i.pF0 / (i.pKs * Math.sqrt((1 - r2) ** 2 + (i.pC * w / i.pKs) ** 2 + 1e-300)); return [['m', 'resonator mass', 'kg', [1, 0, 0, 0], i.pM], ['k', 'spring stiffness', 'N/m', [1, 0, -2, 0], i.pKs], ['X', 'steady-state amplitude (computed from the driven oscillator)', 'm', [0, 1, 0, 0], X], ['c', 'damping coefficient', 'N·s/m', [1, 0, -1, 0], i.pC], ['ω', 'angular drive frequency', 'rad/s', [0, 0, -1, 0], w], ['F₀', 'force amplitude', 'N', [1, 1, -2, 0], i.pF0]]; },
    named: [['r = ω√(m/k)', { 'ω': 1, m: 0.5, k: -0.5 }, v => v < 0.9 ? 'below resonance: resonator moves in phase (effective mass increased)' : v <= 1.1 ? 'near resonance: response limited by damping' : 'above resonance: resonator moves out of phase, the regime where a mass-in-mass cell shows negative effective mass'], ['ζ = c/(2√(km))', { c: 1, k: -0.5, m: -0.5 }, v => v < 0.05 ? `lightly damped, Q ≈ ${F(1 / (2 * Math.max(v, 1e-12)), 3)}: narrow, deep band` : v < 1 ? 'moderately damped: broadened, shallower resonance' : 'overdamped: no resonance'], ['kX/F₀', { k: 1, X: 1, 'F₀': -1 }, v => `dynamic magnification ${F(v, 3)} relative to the static deflection`]] },
  lattice: { label: 'Lattice mechanics (cell-wall bending)', vars: i => [['E_s', 'solid modulus', 'Pa', [1, -1, -2, 0], i.pEs], ['l', 'cell wall length', 'm', [0, 1, 0, 0], i.pLc], ['t', 'wall thickness', 'm', [0, 1, 0, 0], i.pTs], ['σ_ys', 'solid yield strength', 'Pa', [1, -1, -2, 0], i.pSy], ['σ', 'applied stress', 'Pa', [1, -1, -2, 0], i.pSig], ['W', 'specimen width', 'm', [0, 1, 0, 0], i.pW], ['E*', 'effective modulus', 'Pa', [1, -1, -2, 0], i.pEeff]],
    named: [['t/l', { t: 1, l: -1 }, v => v < 0.25 ? 'slender walls: beam (bending) theory applies' : 'stocky walls: shear and axial deformation matter, beam theory underestimates stiffness'], ['E*/E_s', { 'E*': 1, E_s: -1 }, v => `relative stiffness ${F(v, 3)}`], ['(E*/E_s)/(t/l)³', { 'E*': 1, E_s: -1, t: -3, l: 3 }, v => v > 0.3 && v < 30 ? 'of order one: consistent with bending-dominated scaling E* ∝ E_s(t/l)³' : 'far from order one: not bending-dominated, or inconsistent inputs'], ['σ/(σ_ys(t/l)²)', { 'σ': 1, σ_ys: -1, t: -2, l: 2 }, v => v < 0.3 ? 'well below plastic collapse of the walls (σ_pl ~ σ_ys(t/l)²)' : v < 1 ? 'approaching plastic collapse' : 'above the plastic-collapse scaling: wall yielding expected'], ['W/l', { W: 1, l: -1 }, v => v >= 8 ? 'enough cells across the specimen for size-independent effective properties' : 'few cells across: edge effects reduce the measured stiffness']] },
};
function rankOf(A) { const M = A.map(r => [...r]); let r = 0; const piv = []; for (let c = 0; c < M[0].length && r < M.length; c++) { let p = r; for (let i = r + 1; i < M.length; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i; if (Math.abs(M[p][c]) < 1e-9) continue;[M[r], M[p]] = [M[p], M[r]]; for (let i = r + 1; i < M.length; i++) { const f = M[i][c] / M[r][c]; for (let j = c; j < M[0].length; j++) M[i][j] -= f * M[r][j]; } piv.push(c); r++; } return { rank: r, piv }; }
function runPi(inp) {
  const pr = PIP[inp.piProb] || PIP.reactor, V = pr.vars(inp).map(v => ({ s: v[0], name: v[1], unit: v[2], d: v[3], v: Math.max(Math.abs(v[4]), 1e-300) })), n = V.length, W = [];
  const D = [0, 1, 2, 3].map(k => V.map(v => v.d[k])), { rank, piv } = rankOf(D), rep = piv, rest = V.map((_, j) => j).filter(j => !rep.includes(j));
  const used = [0, 1, 2, 3].filter(k => D[k].some(x => x !== 0)), groups = rest.map(j => {
    const A = used.map(k => rep.map(c => D[k][c])), b = used.map(k => -D[k][j]); let a; try { a = N.lstsq(A, b); } catch { a = rep.map(() => 0); }
    const res = N.maxOf(used.map((k, q) => Math.abs(N.dot(A[q], a) - b[q]))), val = Math.exp(Math.log(V[j].v) + N.sum(rep.map((c, q) => a[q] * Math.log(V[c].v))));
    return { expr: [V[j].s, ...rep.map((c, q) => pw(V[c].s, a[q]))].filter(Boolean).join('·'), val, res };
  });
  const byS = Object.fromEntries(V.map((v, j) => [v.s, j])), named = pr.named.map(([nm, ex, interp]) => { let lv = 0; const dim = [0, 0, 0, 0]; for (const s in ex) { const v = V[byS[s]]; lv += ex[s] * Math.log(v.v); for (let k = 0; k < 4; k++) dim[k] += ex[s] * v.d[k]; } const val = Math.exp(N.clamp(lv, -690, 690)); return { nm, val, res: N.maxOf(dim.map(Math.abs)), txt: interp(val) }; });
  const maxRes = Math.max(N.maxOf(groups.map(g => g.res).concat(0)), N.maxOf(named.map(g => g.res)));
  if (rank < used.length) W.push(`The rank of the dimensional matrix (${rank}) is lower than the number of base dimensions present (${used.length}): some dimensions only occur in a fixed combination, so the number of Π groups is n − rank, not n − number of dimensions.`);
  // problem-specific curve
  let curve;
  if (inp.piProb === 'heat') { const Bi = named[0].val, tau = inp.pRhoS * inp.pCp * inp.pLh / Math.max(inp.pH, 1e-300), t = N.linspace(0, 5 * tau, 80); curve = { key: 'curve', title: 'Lumped heat-up response', xlabel: 'Bi·Fo = t/τ', ylabel: 'θ = (T − T∞)/(T₀ − T∞)', lines: [{ name: 'θ = exp(−Bi·Fo)', x: t.map(v => v / tau), y: t.map(v => Math.exp(-v / tau)) }, { name: 'this case', x: [named[2].val], y: [Math.exp(-Math.min(named[2].val, 700))], points: true }], explain: `Dimensionless lumped-capacitance response; it is an accurate description only for Bi < 0.1 (here Bi = ${F(Bi, 3)}). The time constant is τ = ρc_pL/h = ${F(tau, 3)} s.` }; }
  else if (inp.piProb === 'resonator') { const z = named[1].val, r = N.linspace(0.05, 3, 160), meff = r.map(x => 1 / (1 - x * x + (Math.abs(1 - x * x) < 1e-3 ? 1e-3 : 0))); curve = { key: 'curve', title: 'Magnification and effective-mass factor versus frequency ratio', xlabel: 'r = ω/ω₀', ylabel: 'kX/F₀ and m_eff contribution', lines: [{ name: 'kX/F₀', x: r, y: r.map(x => 1 / Math.sqrt((1 - x * x) ** 2 + (2 * z * x) ** 2 + 1e-12)) }, { name: 'undamped 1/(1 − r²) (clipped ±10)', x: r, y: meff.map(v => N.clamp(v, -10, 10)), dash: true }, { name: 'this case', x: [named[0].val], y: [named[2].val], points: true }], explain: 'Dynamic magnification of the driven resonator and the undamped factor 1/(1 − r²) that multiplies the resonator mass in the effective mass of a mass-in-mass cell: it is negative for r > 1, which is the origin of the negative-mass band.' }; }
  else if (inp.piProb === 'lattice') { const x = N.logspace(0.005, 0.5, 60), tl = named[0].val; curve = { key: 'curve', title: 'Relative stiffness scaling', xlabel: 't/l', ylabel: 'E*/E_s', logx: true, logy: true, lines: [{ name: 'bending-dominated (t/l)³', x, y: x.map(v => v ** 3) }, { name: 'stretching-dominated (t/l)', x, y: x, dash: true }, { name: 'this case', x: [N.clamp(tl, 1e-6, 1e3)], y: [Math.max(named[1].val, 1e-300)], points: true }], explain: 'Dimensionless stiffness against wall slenderness with the two limiting scalings (unit prefactor). The position of the case point shows whether the stated effective modulus is consistent with a bending- or stretching-dominated architecture.' }; }
  else { const U = N.logspace(inp.pU / 100, inp.pU * 100, 60); curve = { key: 'curve', title: 'Damköhler and Reynolds numbers versus velocity', xlabel: 'superficial velocity (m/s)', ylabel: 'dimensionless group', logx: true, logy: true, lines: [{ name: 'Da = kL/U', x: U, y: U.map(u => inp.pK * inp.pL / u) }, { name: 'Re = ρUd/μ', x: U, y: U.map(u => inp.pRho * u * inp.pd / inp.pMu), dash: true }, { name: 'plug-flow conversion 1 − exp(−Da)', x: U, y: U.map(u => Math.max(1 - Math.exp(-inp.pK * inp.pL / u), 1e-300)) }], explain: 'How the two controlling groups move when the flow rate is changed at fixed geometry and chemistry. Conversion of a first-order reaction in ideal plug flow depends on Da alone, which is why scale-up at constant Da preserves conversion.' }; }
  const eq = [
    { kind: 'governing', label: 'Buckingham Pi theorem', text: `f(x₁,…,xₙ) = 0  ⇒  φ(Π₁,…,Π₍ₙ₋ᵣ₎) = 0;  n = ${n}, r = rank(D) = ${rank}  ⇒  ${n - rank} independent groups` },
    { kind: 'formula', label: 'Dimensional homogeneity', text: 'Πⱼ = xⱼ·Πᵢ xᵢ^aᵢ over the repeating variables with D_rep·a = −dⱼ (rows M, L, T, Θ)' },
    ...named.map(g => ({ kind: 'formula', label: 'Conventional group', text: `${g.nm} = ${F(g.val, 4)}` })),
  ];
  if (inp.piProb === 'reactor') { try { const io = chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }]); eq.push({ kind: 'chemical', label: 'Example reaction to which k refers (ZnO precipitation)', text: io.molecular }, { kind: 'ionic', label: 'Net ionic equation', text: io.netIonic, note: 'The rate constant k is the pseudo-first-order constant of the limiting species.' }); } catch { /* optional */ } }
  return {
    summary: `The ${n} variables of the "${pr.label}" problem give a dimensional matrix of rank ${rank}, so ${n - rank} independent dimensionless groups describe it. ${named.slice(0, 3).map(g => `${g.nm} = ${F(g.val, 3)} (${g.txt})`).join('; ')}. All groups were checked to be dimensionless to ${F(maxRes, 2)}.`,
    kpis: [
      kp('n', 'Number of variables', n, '', 'Count of dimensional variables entering the problem statement; adding or omitting a variable changes the number of groups.'),
      kp('rank', 'Rank of the dimensional matrix', rank, '', 'Number of independent base-dimension combinations spanned by the variables, obtained by Gaussian elimination of the matrix of exponents.'),
      kp('npi', 'Independent Π groups', n - rank, '', 'n − rank according to the Buckingham theorem: the number of dimensionless parameters that must be matched for complete similarity between model and full scale.'),
      ...named.map((g, i) => kp('pi' + (i + 1), g.nm, g.val, '–', `Conventional dimensionless group evaluated from the inputs. Regime reading: ${g.txt}.`)),
      kp('resid', 'Largest dimension residual', maxRes, '–', 'Largest remaining exponent of any base dimension over all groups; zero (to rounding) confirms that every group is dimensionless.', maxRes < 1e-8),
    ],
    series: [
      { key: 'groups', title: 'Magnitude of the dimensionless groups', xlabel: 'group number (conventional groups, in table order)', ylabel: 'log₁₀ Π', lines: [{ name: 'conventional groups', x: named.map((_, i) => i + 1), y: named.map(g => Math.log10(Math.max(g.val, 1e-300))) }, { name: 'Π = 1', x: [1, named.length], y: [0, 0], dash: true }], explain: 'Order of magnitude of each conventional group. Groups far from unity identify the effect that dominates (large) or can be neglected (small) in a reduced model.' },
      curve,
    ],
    tables: [
      { key: 'vars', title: 'Variables and dimensions', columns: ['Symbol', 'Quantity', 'Unit', 'M', 'L', 'T', 'Θ', 'Value (SI)', 'Role'], rows: V.map((v, j) => [v.s, v.name, v.unit, ...v.d, F(v.v, 4), rep.includes(j) ? 'repeating' : 'forms a Π group']), explain: 'Dimensional matrix of the problem (exponents of mass, length, time and temperature) and the values used. Repeating variables are the pivot columns of the elimination.' },
      { key: 'pi', title: 'Π groups from the Buckingham procedure', columns: ['Group', 'Expression', 'Value', 'Dimension residual'], rows: groups.map((g, i) => [`Π${i + 1}`, g.expr, F(g.val, 4), F(g.res, 2)]), explain: 'A complete, independent set of dimensionless groups built from each non-repeating variable and the repeating variables. Any other valid set is a product of powers of these.' },
      { key: 'named', title: 'Conventional groups and regime interpretation', columns: ['Group', 'Value', 'Interpretation'], rows: named.map(g => [g.nm, F(g.val, 4), g.txt]), explain: 'The same information expressed in the conventional named groups, with the regime each value indicates according to standard order-of-magnitude criteria.' },
    ],
    equations: eq,
    checks: [
      { name: 'All groups dimensionless', pass: maxRes < 1e-8, detail: `largest residual exponent ${F(maxRes, 2)}` },
      { name: 'Number of groups equals n − rank', pass: groups.length === n - rank, detail: `${groups.length} groups for n = ${n}, rank = ${rank}` },
      { name: 'Conventional groups span no more than n − rank dimensions', pass: named.length <= n - rank + 1, detail: `${named.length} conventional groups listed; at most ${n - rank} are independent` },
    ],
    warnings: [...W, 'Dimensional analysis identifies the governing groups and regimes; it does not give the functional relation between them, which needs a model or experiments.'],
    recommendations: [
      { priority: 'high', text: `Design experiments and scale-up around the ${n - rank} groups rather than the ${n} raw variables: this reduces a full-factorial study at 3 levels from ${F(3 ** n, 3)} to ${F(3 ** (n - rank), 3)} runs.`, basis: 'Buckingham Pi reduction of the parameter space.' },
      { priority: 'medium', text: `Keep ${named[0].nm.split(' = ')[0]} and ${named[1].nm.split(' = ')[0]} equal between laboratory and production scale; where both cannot be matched, hold the one that controls the rate-limiting step (${named[0].txt}).`, basis: 'Similarity requires equal dimensionless groups.' },
      { priority: 'medium', text: 'Report measured results in dimensionless form so that data from different sizes and materials collapse onto one curve and can be reused.', basis: 'Dimensional homogeneity.' },
      { priority: 'low', text: 'Fewer, group-based experiments mean less precursor, solvent and energy consumed in process development; plan the test matrix on the groups before ordering materials.', basis: 'Reduction of the experimental matrix by dimensional analysis.', sustainability: true },
    ],
    publish: Object.assign({ 'dim.rank': rank, 'dim.n_groups': n - rank }, Object.fromEntries(named.map(g => ['dim.' + g.nm.split(' = ')[0].replace(/[^A-Za-z]/g, '').slice(0, 8) || 'g', g.val]))),
    numerics: { method: 'Gaussian elimination of the dimensional matrix; exponents of each group by least squares on the repeating variables', dof: n, steps: n - rank, iterations: 1, residual: maxRes, notes: 'Exponents are reported as rational numbers with denominators up to 6.' },
  };
}

// ---------- mode: analytic hierarchy process ----------
const AHP_DEF = '1, 3, 5, 7; 1/3, 1, 3, 5; 1/5, 1/3, 1, 3; 1/7, 1/5, 1/3, 1', RI = [0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49];
function parseAhp(text) { const rows = String(text ?? '').split(/[;\n]+/).map(s => s.trim()).filter(Boolean).map(r => r.split(/[,\s\t]+/).filter(Boolean).map(t => { const m = t.split('/'); const v = m.length === 2 ? Number(m[0]) / Number(m[1]) : Number(t); return Number.isFinite(v) && v > 0 ? v : NaN; })); const n = rows.length; return n >= 2 && n <= 10 && rows.every(r => r.length === n && r.every(Number.isFinite)) ? rows : null; }
function runAhp(inp) {
  const W = []; let A = parseAhp(inp.ahpMatrix); if (!A) { A = parseAhp(AHP_DEF); W.push('The pairwise matrix could not be read as a square 2×2 … 10×10 matrix of positive numbers or fractions; the default example matrix was used instead.'); }
  const n = A.length, names = String(inp.ahpNames ?? '').split(/[,;]+/).map(s => s.trim()).filter(Boolean); while (names.length < n) names.push(`criterion ${names.length + 1}`);
  let recDev = 0; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) recDev = Math.max(recDev, Math.abs(A[i][j] * A[j][i] - 1));
  if (recDev > 0.02) W.push(`The matrix is not reciprocal (largest |aᵢⱼ·aⱼᵢ − 1| = ${F(recDev, 3)}); the upper triangle was used and the lower triangle replaced by reciprocals.`);
  if (recDev > 0.02) for (let i = 0; i < n; i++) { A[i][i] = 1; for (let j = 0; j < i; j++) A[i][j] = 1 / A[j][i]; }
  let w = new Array(n).fill(1 / n), lam = n, it = 0; const hist = [];
  for (; it < 200; it++) { const y = N.matVec(A, w), s = N.sum(y), w2 = y.map(v => v / s), d = N.maxOf(w2.map((v, i) => Math.abs(v - w[i]))); lam = s; w = w2; hist.push(Math.max(d, 1e-17)); if (d < 1e-13) break; }
  const Aw = N.matVec(A, w); lam = N.mean(Aw.map((v, i) => v / w[i]));
  const CI = n > 2 ? (lam - n) / (n - 1) : 0, CR = n > 2 ? CI / RI[n - 1] : 0, gm = A.map(r => Math.exp(N.mean(r.map(Math.log)))), gs = N.sum(gm), wg = gm.map(v => v / gs), dGm = N.maxOf(w.map((v, i) => Math.abs(v - wg[i])));
  let worst = [0, 1, 0]; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const e = Math.abs(Math.log(A[i][j] * w[j] / w[i])); if (e > worst[2]) worst = [i, j, e]; }
  const ok = CR < 0.1, top = N.argmax(w), ord = w.map((_, i) => i).sort((a, b) => w[b] - w[a]);
  if (!ok) W.push(`Consistency ratio ${F(CR, 3)} ≥ 0.1: the judgements are too inconsistent for the weights to be relied on. Revisit the comparison ${names[worst[0]]} vs ${names[worst[1]]} first.`);
  return {
    summary: `The ${n}×${n} pairwise-comparison matrix gives priority weights ${ord.map(i => `${names[i]} ${F(w[i], 3)}`).join(', ')}. The principal eigenvalue is λ_max = ${F(lam, 5)}, the consistency index ${F(CI, 3)} and the consistency ratio ${F(CR, 3)}, which is ${ok ? 'below' : 'above'} the 0.1 acceptance limit. Power iteration converged in ${it + 1} steps.`,
    kpis: [
      kp('top', 'Highest-priority criterion', names[top], '', `Criterion with the largest priority weight (${F(w[top], 3)}); it dominates a weighted-sum selection made with these weights.`),
      kp('w_top', 'Largest weight', w[top], '–', 'Component of the normalised principal eigenvector for the highest-priority criterion; all weights sum to one.'),
      kp('w_min', 'Smallest weight', N.minOf(w), '–', 'Weight of the least important criterion; criteria below a few per cent hardly influence the ranking and may be dropped.'),
      kp('lambda', 'Principal eigenvalue λ_max', lam, '–', `For a perfectly consistent matrix λ_max equals the matrix order n = ${n}; the excess measures inconsistency of the judgements.`),
      kp('CI', 'Consistency index', CI, '–', 'CI = (λ_max − n)/(n − 1): average of the remaining eigenvalues with opposite sign; zero for perfectly transitive judgements.'),
      kp('CR', 'Consistency ratio', CR, '–', `CI divided by the random index RI = ${RI[n - 1]} of a random reciprocal matrix of the same order; Saaty's acceptance limit is 0.1.`, ok),
      kp('iters', 'Power iterations', it + 1, '', 'Number of matrix–vector multiplications needed for the weight vector to converge to 10⁻¹³; convergence is geometric with the ratio of the two largest eigenvalue moduli.'),
      kp('dgm', 'Largest difference to geometric-mean weights', dGm, '–', 'Maximum deviation between eigenvector weights and row-geometric-mean weights; the two coincide for a consistent matrix, so a large value is another inconsistency indicator.', dGm < 0.02),
    ],
    series: [
      { key: 'weights', title: 'Priority weights', xlabel: 'criterion number', ylabel: 'weight', lines: [{ name: 'principal eigenvector', x: w.map((_, i) => i + 1), y: w }, { name: 'row geometric mean', x: w.map((_, i) => i + 1), y: wg, dash: true }], explain: 'Priority weights of the criteria from the eigenvector method and from the geometric-mean approximation. Close agreement indicates consistent judgements.' },
      { key: 'conv', title: 'Power-iteration convergence', xlabel: 'iteration', ylabel: 'largest weight change', logy: true, lines: [{ name: 'max |Δw|', x: hist.map((_, i) => i + 1), y: hist }], explain: 'Change of the weight vector per iteration on a logarithmic axis; the straight decline shows geometric convergence to the principal eigenvector.' },
    ],
    tables: [
      { key: 'w', title: 'Criteria weights', columns: ['Criterion', 'Eigenvector weight', 'Geometric-mean weight', 'Rank'], rows: w.map((v, i) => [names[i], F(v, 4), F(wg[i], 4), ord.indexOf(i) + 1]), explain: 'Weights to enter as the weighting factors of the ranking mode. They express the relative importance implied by the pairwise judgements.' },
      { key: 'A', title: 'Pairwise-comparison matrix as used', columns: ['', ...names.slice(0, n)], rows: A.map((r, i) => [names[i], ...r.map(v => F(v, 3))]), explain: 'Entry aᵢⱼ states how much more important criterion i is than criterion j on the 1–9 Saaty scale; the lower triangle holds the reciprocals.' },
      { key: 'incons', title: 'Judgement residuals', columns: ['Pair', 'Judgement aᵢⱼ', 'Implied wᵢ/wⱼ', 'Ratio'], rows: A.flatMap((r, i) => r.map((v, j) => j > i ? [`${names[i]} vs ${names[j]}`, F(v, 3), F(w[i] / w[j], 3), F(v * w[j] / w[i], 3)] : null).filter(Boolean)), explain: 'Each judgement compared with the ratio implied by the final weights; ratios far from 1 identify the judgements that cause inconsistency.' },
    ],
    equations: [
      { kind: 'governing', label: 'Eigenvector priority', text: 'A·w = λ_max·w,  Σ wᵢ = 1,  aᵢⱼ = 1/aⱼᵢ > 0' },
      { kind: 'formula', label: 'Power iteration', text: 'w⁽ᵏ⁺¹⁾ = A·w⁽ᵏ⁾ / ‖A·w⁽ᵏ⁾‖₁' },
      { kind: 'formula', label: 'Consistency', text: `CI = (λ_max − n)/(n − 1) = ${F(CI, 4)};  CR = CI/RI(n) = ${F(CR, 4)}  (accept if CR < 0.1)` },
    ],
    checks: [
      { name: 'Consistency ratio below 0.1', pass: ok, detail: `CR = ${F(CR, 4)}` },
      { name: 'Reciprocal matrix', pass: recDev <= 0.02, detail: `largest |aᵢⱼ·aⱼᵢ − 1| = ${F(recDev, 3)}` },
      { name: 'λ_max ≥ n', pass: lam >= n - 1e-9, detail: `λ_max = ${F(lam, 6)}, n = ${n}` },
      { name: 'Eigen-residual', pass: N.maxOf(Aw.map((v, i) => Math.abs(v - lam * w[i]))) < 1e-6, detail: `‖A·w − λw‖∞ = ${F(N.maxOf(Aw.map((v, i) => Math.abs(v - lam * w[i]))), 2)}` },
      { name: 'Weights sum to one', pass: Math.abs(N.sum(w) - 1) < 1e-12, detail: `Σw = ${F(N.sum(w), 8)}` },
    ],
    warnings: W,
    recommendations: [
      { priority: ok ? 'low' : 'high', text: ok ? `Use the weights in the ranking mode (performance ${F(w[0], 3)} …); the judgements are acceptably consistent.` : `Re-elicit the judgement "${names[worst[0]]} vs ${names[worst[1]]}" (entered ${F(A[worst[0]][worst[1]], 3)}, implied ${F(w[worst[0]] / w[worst[1]], 3)}) and recompute until CR < 0.1.`, basis: 'Saaty consistency test and judgement residuals.' },
      { priority: 'medium', text: `Criteria with weight below 0.05 (${w.filter(v => v < 0.05).length} here) can be handled as hard constraints instead of weighted objectives, which simplifies the decision.`, basis: 'Small eigenvector components.' },
      { priority: 'medium', text: 'Have at least two stakeholders fill the matrix independently and aggregate by the element-wise geometric mean; it preserves reciprocity and exposes disagreement.', basis: 'Group AHP practice.' },
      { priority: 'medium', text: 'Include cost, embodied CO₂ and critical-raw-material risk as explicit criteria in the comparison so that their weight is decided deliberately, not left at zero.', basis: 'Multi-criteria decisions omit what is not listed.', sustainability: true },
    ],
    publish: Object.assign({ 'req.ahp_CR': CR, 'req.ahp_lambda': lam, 'opt.objective': w[top] }, Object.fromEntries(w.slice(0, 6).map((v, i) => [`req.w${i + 1}`, v]))),
    exports: [{ name: 'ahp_weights.csv', mime: 'text/csv', content: 'criterion,weight\n' + w.map((v, i) => `"${names[i]}",${v}`).join('\n'), description: 'Criteria weights from the principal eigenvector.' }],
    numerics: { method: 'Power iteration for the principal eigenvector of the reciprocal comparison matrix', dof: n, steps: it + 1, iterations: it + 1, residual: hist.at(-1), notes: 'Random indices RI after Saaty for n ≤ 10.' },
  };
}

// ---------- mode: feasibility and specification sheet ----------
const VERIF = { 'C-RHO': ['helium pycnometry or Archimedes density', 2], 'C-E': ['tensile test or nanoindentation', 14], 'C-SY': ['tensile / compression test', 14], 'C-NU': ['digital image correlation of lateral strain', 14], 'C-K': ['laser-flash diffusivity with c_p and density', 10], 'C-SIG': ['four-point probe or van der Pauw measurement', 11], 'C-CTE': ['dilatometry over the service range', 10], 'C-EGLO': ['UV–vis absorption (Tauc analysis)', 11], 'C-EGHI': ['UV–vis absorption (Tauc analysis)', 11], 'C-TOP': ['thermal analysis (DSC/TGA) and creep or ageing test', 22], 'C-COST': ['supplier quotation at the planned volume', 28], 'C-CO2': ['life-cycle inventory of the chosen route', 28] };
function runSpec(inp) {
  const S = screen(inp), { fn, idx, cons, surv, ord, score } = S, W = [], forced = findMaterial(String(inp.specId ?? '').trim());
  let top = ord.length ? surv[ord[0]] : null, forcedNote = '';
  if (String(inp.specId ?? '').trim()) { if (forced) { top = S.rows.find(r => r.m.id === forced.id) || (() => { const g = getter(forced); return { m: forced, g, st: cons.map(c => { const v = c.get(g), mg = v === null ? null : (c.sense === 'max' ? c.lim - v : v - c.lim) / Math.max(Math.abs(c.lim), 1e-12); return v === null ? { v: null, s: 'n.d.', mg: null } : { v, mg, s: mg >= 0 ? 'pass' : 'FAIL' }; }), nd: 0, M: idx.map(d => num(d[2](g))) }; })(); forcedNote = ` (chosen by library id "${forced.id}")`; } else W.push(`Library id "${inp.specId}" was not found: the top-ranked candidate is used instead.`); }
  const proj = String(inp.project || 'Unnamed project').slice(0, 80);
  if (!top) {
    W.push('No candidate satisfies the hard constraints, so no specification can be assembled. The requirement list is still exported for review.');
  }
  const reqs = [{ id: 'REQ-001', text: `The material shall provide the function "${fn.label}".`, acc: `${idx[0][0]} reported and reviewed`, st: top ? (top.M[0] !== null ? 'pass' : 'n.d.') : 'open', val: top && top.M[0] !== null ? `${F(top.M[0], 4)} ${idx[0][1]}` : 'n.d.', ver: 'functional test on a representative specimen', suite: fn.next[0], mg: null },
    ...cons.map((c, j) => ({ id: `REQ-${String(j + 2).padStart(3, '0')}`, text: `${c.name} shall be ${c.sense === 'max' ? 'at most' : 'at least'} ${F(c.lim, 4)} ${c.unit}.`, acc: `${c.sense === 'max' ? '≤' : '≥'} ${F(c.lim, 4)} ${c.unit}`, st: top ? top.st[j].s : 'open', val: top && top.st[j].v !== null ? `${F(top.st[j].v, 4)} ${c.unit}` : 'n.d.', ver: VERIF[c.id][0], suite: VERIF[c.id][1], mg: top ? top.st[j].mg : null, cid: c.id }))];
  const nPass = reqs.filter(r => r.st === 'pass').length, nFail = reqs.filter(r => r.st === 'FAIL').length, nNd = reqs.length - nPass - nFail, feas = reqs.length ? nPass / reqs.length : 0, mgs = reqs.map(r => r.mg).filter(v => v !== null), minMg = mgs.length ? N.minOf(mgs) : 0;
  const m = top ? top.m : null, cls = m ? findClass(m.classId) : null, syn = m ? (m.synthesis || [])[0] : null, nextS = m ? [...new Set([...fn.next, ...(m.suites || []).slice(0, 6), 23, 24])].filter(i => SUITE[i]).slice(0, 8) : fn.next;
  if (nNd && top) W.push(`${nNd} requirement${nNd === 1 ? '' : 's'} cannot be verified from library data (status "n.d."): compliance is unknown, not assumed.`);
  if (cls && cls.caveats) W.push(`Class caveat (${cls.name}): ${cls.caveats}`);
  W.push('The specification is a preliminary screening document built from typical library values; it is not a qualified material specification until the acceptance tests have been carried out.');
  const spec = { project: proj, function: fn.label, operatingTemperature_C: inp.Top, material: m ? { id: m.id, name: m.name, formula: m.formula, class: m.classId, kind: m.kind, dimensionality: m.dimensionality } : null, requirements: reqs.map(r => ({ id: r.id, statement: r.text, acceptance: r.acc, libraryValue: r.val, status: r.st, margin: r.mg, verification: r.ver, analysisSuite: r.suite })), properties: m ? Object.fromEntries(Object.entries(m.props || {}).filter(([, v]) => num(v) !== null)) : {}, effective: m ? m.effective : null, synthesis: syn ? { route: syn.route, T_C: syn.T_C, P_bar: syn.P_bar, time_h: syn.time_h, precursors: syn.precursors } : null, hazards: m ? m.hazards : [], nextSuites: nextS.map(i => ({ id: i, name: SUITE[i] })), status: 'preliminary — not validated' };
  const txt = [`MATERIAL SPECIFICATION (preliminary) — ${proj}`, `Function: ${fn.label}`, `Material: ${m ? `${m.name}${m.formula ? ' [' + m.formula + ']' : ''} — ${cls ? cls.name : m.classId}` : 'no feasible candidate'}`, '', 'Requirements and acceptance criteria', ...reqs.map(r => `${r.id}  ${r.text}  Acceptance: ${r.acc}.  Library value: ${r.val}.  Status: ${r.st}.  Verification: ${r.ver}.`), '', syn ? `Reference synthesis / fabrication route: ${syn.route} (${syn.T_C ? syn.T_C.join('–') + ' °C' : 'temperature n.d.'})` : 'Reference route: n.d.', m && m.hazards && m.hazards.length ? 'Hazards: ' + m.hazards.join('; ') : 'Hazards: n.d.', 'Recommended next analyses: ' + nextS.map(i => `suite ${i} ${SUITE[i]}`).join('; '), 'Status: preliminary screening result; acceptance tests outstanding.'].join('\n');
  const pub = publishTop(top, { 'opt.objective': top && ord.length && top === surv[ord[0]] ? score[ord[0]] : feas, 'req.feasibility': feas, 'req.T_op_C': inp.Top, 'proc.T_C': syn && syn.T_C ? (syn.T_C[0] + syn.T_C[1]) / 2 : inp.Top });
  if (syn && syn.time_h) pub['proc.time'] = 3600 * (syn.time_h[0] + syn.time_h[1]) / 2; if (syn && syn.P_bar) pub['proc.P'] = 1e5 * (syn.P_bar[0] + syn.P_bar[1]) / 2;
  const show = ord.slice(0, 8);
  return {
    summary: m ? `Specification "${proj}" for ${fn.label.toLowerCase()} is assembled around ${m.name}${forcedNote}: ${nPass} of ${reqs.length} requirements pass on library data, ${nFail} fail and ${nNd} cannot be verified, a feasibility score of ${pct(feas)}. The smallest design margin is ${F(minMg, 3)}. ${syn ? `Reference route: ${syn.route}.` : ''} The sheet lists acceptance criteria, verification methods and the follow-on analyses for each requirement.` : `No library candidate satisfies the hard constraints for ${fn.label.toLowerCase()}, so specification "${proj}" contains the ${reqs.length} requirements with open status only. Relax the limits or widen the candidate classes in the ranking mode.`,
    kpis: [
      kp('material', 'Specified material', m ? m.name : 'none', '', 'Library entry on which the specification is built: the top-ranked feasible candidate unless a library id was entered.', !!m),
      kp('n_req', 'Requirements', reqs.length, '', 'Number of traceable requirements generated from the stated function and active limits; each carries an id, an acceptance criterion and a verification method.'),
      kp('n_pass', 'Requirements passed on library data', nPass, '', 'Requirements whose acceptance criterion is met by the typical library value; they still need confirmation by test on the actual product form.', nPass === reqs.length),
      kp('n_fail', 'Requirements failed', nFail, '', 'Requirements violated by the library value of the specified material; non-zero only when a material was forced by id.', nFail === 0),
      kp('n_nd', 'Requirements not verifiable', nNd, '', 'Requirements for which the library has no value; they are open items for measurement, not assumed compliances.', nNd === 0),
      kp('feas', 'Feasibility score', feas, '–', 'Fraction of requirements demonstrably met on present data. It rises when open items are closed by measurement and falls when a test fails.', feas > 0.8),
      kp('margin', 'Smallest design margin', minMg, '–', 'Smallest normalised margin over all verifiable limits; it identifies the requirement most at risk from property scatter.', m ? minMg > 0.1 : null),
      kp('Tsyn', 'Reference process temperature', syn && syn.T_C ? (syn.T_C[0] + syn.T_C[1]) / 2 : 0, '°C', 'Mid-point of the temperature window of the first library synthesis or fabrication route of the specified material (0 when none is recorded).'),
    ],
    series: [
      { key: 'margins', title: 'Design margin per requirement', xlabel: 'requirement number', ylabel: 'normalised margin', lines: [{ name: 'margin', x: reqs.map((_, i) => i + 1), y: reqs.map(r => r.mg === null ? 0 : N.clamp(r.mg, -5, 5)) }, { name: 'limit (0)', x: [1, Math.max(reqs.length, 2)], y: [0, 0], dash: true }], explain: 'Normalised margin of the specified material against each requirement (clipped to ±5; zero for requirements without data). Negative values are non-compliances; small positive ones need tight property control.' },
      { key: 'alts', title: 'Scores of alternative candidates', xlabel: 'rank', ylabel: 'weighted-sum score', lines: [{ name: 'feasible alternatives', x: show.length ? show.map((_, k) => k + 1) : [1], y: show.length ? show.map(i => score[i]) : [0] }], explain: 'Weighted-sum scores of the best feasible candidates, showing how far the alternatives lie behind the specified material and which could serve as second source.' },
    ],
    tables: [
      { key: 'req', title: 'Requirement traceability', columns: ['Id', 'Requirement', 'Acceptance criterion', 'Library value', 'Status', 'Verification method', 'Analysis suite'], rows: reqs.map(r => [r.id, r.text, r.acc, r.val, r.st, r.ver, `${r.suite} ${SUITE[r.suite] || ''}`]), explain: 'Each requirement with its identifier, acceptance criterion, present evidence and the method and suite by which it is to be verified. "n.d." marks open items.' },
      { key: 'props', title: 'Property sheet of the specified material', columns: ['Property', 'Value (SI)'], rows: m ? [...Object.entries(m.props || {}).filter(([, v]) => num(v) !== null).map(([k, v]) => [k, F(v, 4)]), ...Object.entries(m.effective || {}).filter(([, v]) => num(v) !== null).map(([k, v]) => ['effective ' + k, F(v, 4)])] : [['–', '–']], explain: 'Typical library values of the specified material in SI units; effective entries are narrow-band or architecture-specific values with the conditions given in the library note.' },
      { key: 'next', title: 'Recommended next suites and experiments', columns: ['Suite', 'Purpose'], rows: nextS.map(i => [`${i} ${SUITE[i]}`, i === 23 ? 'measure the properties marked n.d. and confirm the governing ones' : i === 24 ? 'quantify model credibility before relying on predictions' : fn.next.includes(i) ? 'function-specific design analysis' : 'material-specific analysis listed in the library']), explain: 'Follow-on analyses needed to turn the preliminary specification into a verified design, in the order suggested by the function and the material record.' },
    ],
    equations: [
      { kind: 'formula', label: 'Feasibility score', text: 'F = N_pass / N_requirements;  open items (no data) count as not yet demonstrated' },
      { kind: 'formula', label: 'Design margin', text: 'm = (limit − value)/|limit| (upper limit),  m = (value − limit)/|limit| (lower limit)' },
      { kind: 'formula', label: `Performance index: ${fn.label}`, text: fn.eq },
      ...(m ? chemEquations(m) : []),
    ],
    checks: [
      { name: 'Material assigned', pass: !!m, detail: m ? m.id : 'no feasible candidate' },
      { name: 'No failed requirement', pass: nFail === 0 && !!m, detail: `${nFail} failed` },
      { name: 'All requirements verifiable from data', pass: nNd === 0 && !!m, detail: `${nNd} without library data` },
      { name: 'Every requirement has an acceptance criterion and verification method', pass: reqs.every(r => r.acc && r.ver), detail: `${reqs.length} requirements traced` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: nNd || !m ? `Close the ${m ? nNd : reqs.length} open requirement${nNd === 1 ? '' : 's'} by measurement before design freeze${m ? `: ${reqs.filter(r => r.st === 'n.d.').slice(0, 3).map(r => `${r.id} by ${r.ver}`).join('; ')}` : ''}.` : 'All requirements are met on library data: confirm the two smallest margins by test on the actual product form and supplier.', basis: 'Requirement traceability table.' },
      { priority: 'medium', text: m && minMg < 0.1 && mgs.length ? `The tightest requirement has a margin of only ${F(minMg, 3)}: specify an incoming-inspection limit for that property or select the next candidate with a larger margin.` : 'Fix supplier, grade and product form in the specification; library values are typical and vary with processing.', basis: 'Design margins against property scatter.' },
      { priority: 'medium', text: `Proceed with ${nextS.slice(0, 3).map(i => `suite ${i} (${SUITE[i]})`).join(', ')} and feed the results back into this sheet.`, basis: 'Function- and material-specific analysis chain.' },
      { priority: 'medium', text: m ? `${(m.hazards || [])[0] ? 'Add the handling requirement "' + m.hazards[0] + '" to the specification. ' : ''}State embodied CO₂ and end-of-life route as explicit requirements (suite 28) so that sustainability is verified like any other property.` : 'When a candidate is found, add embodied CO₂, hazard and end-of-life requirements to the list.', basis: 'Library hazard record; life-cycle requirements.', sustainability: true },
    ],
    publish: pub,
    exports: [{ name: 'material_specification.json', mime: 'application/json', content: JSON.stringify(spec, null, 1), description: 'Machine-readable preliminary specification with requirement traceability.' }, { name: 'material_specification.txt', mime: 'text/plain', content: txt, description: 'Text specification sheet with acceptance criteria and next steps.' }],
    numerics: { method: 'Constraint screening and weighted-sum ranking followed by rule-based requirement traceability', dof: S.pool.length, steps: reqs.length, iterations: 1, residual: 0, notes: 'Requirement status is derived from library values only.' },
  };
}

const R = (mode) => ({ mode: [mode] }), RS = { mode: ['rank', 'spec'] };
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'rank', group: 'Model', options: [{ value: 'rank', label: 'Requirements → candidate ranking' }, { value: 'alloc', label: 'Requirement allocation, tolerances and margins' }, { value: 'pi', label: 'Dimensional analysis (Buckingham Pi)' }, { value: 'ahp', label: 'Analytic hierarchy process (weights)' }, { value: 'spec', label: 'Feasibility and specification sheet' }], help: 'Selects the requirements-engineering task. Ranking and specification screen the materials library; the other modes allocate tolerances and reliability, reduce variables to dimensionless groups, or derive criteria weights.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 99999, step: 1, sens: false, group: 'Model', help: 'Seed of the pseudo-random generator used for the Monte Carlo tolerance simulation. The same seed reproduces the same result exactly.' },
  { key: 'func', label: 'Application function', type: 'select', value: 'stiff', group: 'Function', showIf: RS, options: Object.entries(FUNCS).map(([value, f]) => ({ value, label: f.label })), help: 'What the material has to do. It selects the candidate classes and the performance index used for ranking.' },
  { key: 'kind', label: 'Material kind', type: 'select', value: 'all', group: 'Function', showIf: RS, options: [{ value: 'all', label: 'All library entries' }, { value: 'nanomaterial', label: 'Nanomaterials' }, { value: 'metamaterial', label: 'Metamaterials' }, { value: 'constituent', label: 'Constituent / bulk materials' }], help: 'Restricts the candidate pool to one kind of library entry. Nanostructured metamaterials are included under both nanomaterials and metamaterials.' },
  { key: 'classOnly', label: 'Restrict to function-specific classes', type: 'bool', value: true, group: 'Function', showIf: RS, help: 'When on, functions tied to a material class (auxetic, negative stiffness, negative thermal expansion, …) only consider members of those classes. Switch off to let every entry with data compete.' },
  { key: 'loading', label: 'Loading mode (structural index)', type: 'select', value: 'beam', group: 'Function', showIf: { mode: ['rank', 'spec'], func: ['stiff'] }, options: [{ value: 'tie', label: 'Tie in tension: E/ρ, σ_y/ρ' }, { value: 'beam', label: 'Beam in bending: E^(1/2)/ρ, σ_y^(2/3)/ρ' }, { value: 'plate', label: 'Plate in bending: E^(1/3)/ρ, σ_y^(1/2)/ρ' }], help: 'Loading geometry for the lightweight-structure function; it fixes the exponent of the Ashby index because the free section dimension differs.' },
  { key: 'dchar', label: 'Characteristic nanoscale size', unit: 'nm', value: 50, min: 1, max: 1000, log: true, group: 'Function', showIf: { mode: ['rank', 'spec'], func: ['sensing'] }, help: 'Particle diameter, rod diameter or sheet thickness used for the geometric specific surface area SSA = s/(ρd). Smaller size raises the surface area in inverse proportion.' },
  { key: 'EgTarget', label: 'Target band gap', unit: 'eV', value: 3.3, min: 0.1, max: 8, step: 0.05, group: 'Function', showIf: { mode: ['rank', 'spec'], func: ['semi'] }, bus: 'prop.bandgap_eV', help: 'Desired band gap for the optoelectronic function; candidates are scored by relative closeness to it. 3.3 eV corresponds to a 376 nm absorption edge.' },
  { key: 'rhoMax', label: 'Maximum density', unit: 'kg/m³', value: 5000, min: 1, max: 25000, log: true, group: 'Limits', showIf: RS, help: 'Upper limit on (effective) density. Set to the upper bound 25000 to switch the constraint off.' },
  { key: 'Emin', label: 'Minimum Young’s modulus', unit: 'GPa', value: 1, min: 0, max: 1200, group: 'Limits', showIf: RS, help: 'Lower limit on (effective) Young’s modulus. Zero switches the constraint off.' },
  { key: 'yMin', label: 'Minimum yield strength', unit: 'MPa', value: 0, min: 0, max: 5000, group: 'Limits', showIf: RS, help: 'Lower limit on yield or compressive strength. Zero switches the constraint off; many ceramics and nanomaterials have no library yield value and are then not screened.' },
  { key: 'nuMax', label: 'Maximum Poisson’s ratio', value: 0.5, min: -1, max: 0.5, step: 0.05, group: 'Limits', showIf: RS, help: 'Upper limit on Poisson’s ratio; enter a negative value to demand auxetic behaviour. 0.5 switches the constraint off.' },
  { key: 'kMin', label: 'Minimum thermal conductivity', unit: 'W/m/K', value: 0, min: 0, max: 3000, group: 'Limits', showIf: RS, help: 'Lower limit on thermal conductivity for heat-spreading duties. Zero switches the constraint off.' },
  { key: 'sigMin', label: 'Minimum electrical conductivity', unit: 'S/m', value: 0, min: 0, max: 1e8, group: 'Limits', showIf: RS, help: 'Lower limit on electrical conductivity for conductors, electrodes and shielding. Zero switches the constraint off; values for architected lattices are not held in the library and are then not screened.' },
  { key: 'cteMax', label: 'Maximum thermal expansion coefficient', unit: '10⁻⁶/K', value: 100, min: -20, max: 100, group: 'Limits', showIf: RS, help: 'Upper limit on the linear expansion coefficient; zero or a negative value demands zero or negative thermal expansion. 100 switches the constraint off.' },
  { key: 'EgMin', label: 'Minimum band gap', unit: 'eV', value: 0, min: 0, max: 10, step: 0.1, group: 'Limits', showIf: RS, help: 'Lower limit on the band gap (transparency or insulation requirement). Zero switches the constraint off.' },
  { key: 'EgMax', label: 'Maximum band gap', unit: 'eV', value: 10, min: 0, max: 10, step: 0.1, group: 'Limits', showIf: RS, help: 'Upper limit on the band gap (absorption requirement). 10 switches the constraint off.' },
  { key: 'Top', label: 'Operating temperature', unit: '°C', value: 100, min: -200, max: 2500, group: 'Limits', showIf: RS, bus: 'proc.T_C', help: 'Highest service temperature. Candidates are screened with the homologous-temperature rule T_op ≤ f·T_m; phase transitions, oxidation and glass transition are not covered by this rule.' },
  { key: 'homol', label: 'Allowed homologous temperature f', value: 0.5, min: 0.2, max: 0.95, step: 0.05, group: 'Limits', showIf: RS, help: 'Fraction of the absolute melting point up to which a material is admitted. About 0.3–0.5 is a usual creep-onset screening value for load-bearing use; higher values suit unloaded functional parts.' },
  { key: 'costMax', label: 'Maximum material cost', unit: 'USD/kg', value: 1e6, min: 0.1, max: 1e6, log: true, group: 'Limits', showIf: RS, help: 'Upper limit on the material price. The upper bound 10⁶ switches the constraint off; only part of the library holds cost data.' },
  { key: 'co2Max', label: 'Maximum embodied CO₂', unit: 'kg/kg', value: 1e4, min: 0.1, max: 1e4, log: true, group: 'Limits', showIf: RS, help: 'Upper limit on embodied carbon of primary production. The upper bound 10⁴ switches the constraint off.' },
  { key: 'wPerf', label: 'Weight: performance index', value: 0.5, min: 0, max: 1, step: 0.05, group: 'Weights', showIf: RS, help: 'Relative importance of the function-specific performance index. Weights are renormalised to sum to one.' },
  { key: 'wRho', label: 'Weight: low density', value: 0.2, min: 0, max: 1, step: 0.05, group: 'Weights', showIf: RS, help: 'Relative importance of low density (mass). Use zero if mass is already contained in the performance index.' },
  { key: 'wCost', label: 'Weight: low cost', value: 0.1, min: 0, max: 1, step: 0.05, group: 'Weights', showIf: RS, help: 'Relative importance of low material cost per kilogram. Candidates without cost data are scored on the remaining criteria.' },
  { key: 'wCo2', label: 'Weight: low embodied CO₂', value: 0.1, min: 0, max: 1, step: 0.05, group: 'Weights', showIf: RS, help: 'Relative importance of low embodied carbon. Candidates without data are scored on the remaining criteria.' },
  { key: 'wTemp', label: 'Weight: temperature capability', value: 0.1, min: 0, max: 1, step: 0.05, group: 'Weights', showIf: RS, help: 'Relative importance of a high melting point as a proxy for thermal stability.' },
  { key: 'objX', label: 'Pareto objective on x', type: 'select', value: 'density', group: 'Pareto', showIf: R('rank'), options: Object.entries(OBJ).map(([value, o]) => ({ value, label: `${o[0]} (${o[1]}imise)` })), help: 'First objective of the two-objective Pareto analysis and x axis of the selection chart.' },
  { key: 'objY', label: 'Pareto objective on y', type: 'select', value: 'index', group: 'Pareto', showIf: R('rank'), options: Object.entries(OBJ).map(([value, o]) => ({ value, label: `${o[0]} (${o[1]}imise)` })), help: 'Second objective of the Pareto analysis and y axis of the selection chart.' },
  { key: 'project', label: 'Project / product designation', type: 'text', value: 'Lightweight panel demonstrator', group: 'Specification', showIf: R('spec'), help: 'Free-text designation printed on the specification sheet and stored in the exported record.' },
  { key: 'specId', label: 'Force library id (optional)', type: 'text', value: '', group: 'Specification', showIf: R('spec'), help: 'Leave empty to specify the top-ranked candidate. Enter a library id (for example zno-nanorods or reentrant-honeycomb) to check a material of your own choice against the requirements.' },
  // allocation
  { key: 'stack', label: 'Tolerance contributors (name, nominal, tolerance, sensitivity)', type: 'table', columns: ['contributor', 'nominal_mm', 'tolerance_mm', 'sensitivity'], value: STACK_DEF, group: 'Tolerance stack-up', showIf: R('alloc'), help: 'One row per dimension in the chain: nominal value, symmetric ± tolerance and sensitivity (∂assembly/∂dimension, default 1; use the number of repeats for repeated features).' },
  { key: 'stackLimit', label: 'Allowable assembly deviation', unit: '± mm', value: 0.1, min: 0.001, max: 10, log: true, group: 'Tolerance stack-up', showIf: R('alloc'), help: 'Largest deviation of the assembled dimension that the function tolerates; compared with worst-case and statistical stacks.' },
  { key: 'kSigma', label: 'Tolerance = k standard deviations', value: 3, min: 1, max: 6, step: 0.5, group: 'Tolerance stack-up', showIf: R('alloc'), help: 'Process capability assumed for each contributor: its tolerance equals k standard deviations. k = 3 corresponds to C_p = 1.' },
  { key: 'cap', label: 'Capacity (strength)', unit: 'Pa', value: 2.5e8, min: 1e3, max: 1e11, log: true, group: 'Margins', showIf: R('alloc'), bus: 'prop.yield', help: 'Mean resistance of the material or component, for example yield strength or buckling stress.' },
  { key: 'dem', label: 'Demand (applied stress)', unit: 'Pa', value: 9e7, min: 1e3, max: 1e11, log: true, group: 'Margins', showIf: R('alloc'), bus: 'mech.max_stress', help: 'Mean applied stress or load effect in the same units as the capacity.' },
  { key: 'covCap', label: 'Capacity scatter (COV)', unit: '%', value: 7, min: 0, max: 50, group: 'Margins', showIf: R('alloc'), help: 'Coefficient of variation of the capacity. Wrought metals are typically a few per cent; ceramics, printed lattices and nanomaterials considerably more.' },
  { key: 'covDem', label: 'Demand scatter (COV)', unit: '%', value: 10, min: 0, max: 50, group: 'Margins', showIf: R('alloc'), help: 'Coefficient of variation of the demand (load uncertainty).' },
  { key: 'fosReq', label: 'Required factor of safety', value: 1.5, min: 1, max: 10, step: 0.1, group: 'Margins', showIf: R('alloc'), help: 'Factor of safety demanded by the design rule or customer; the margin of safety is computed against it.' },
  { key: 'derate', label: 'Derating factor', value: 0.8, min: 0.1, max: 1, step: 0.05, group: 'Margins', showIf: R('alloc'), help: 'Fraction of the rated capacity that may be used in service (temperature, ageing or environment allowance). 1 means no derating.' },
  { key: 'Rsys', label: 'System reliability target', value: 0.99, min: 0.5, max: 0.999999, group: 'Reliability', showIf: R('alloc'), help: 'Required probability that the whole series system survives the mission time.' },
  { key: 'nSub', label: 'Number of series subsystems', value: 5, min: 1, max: 50, step: 1, group: 'Reliability', showIf: R('alloc'), help: 'Number of items that must all work (series reliability block diagram).' },
  { key: 'relW', label: 'Complexity weights (optional)', type: 'text', value: '1, 2, 1, 3, 1', group: 'Reliability', showIf: R('alloc'), help: 'Comma-separated relative complexity or expected failure share of each subsystem; a larger weight receives a larger share of the allowed unreliability. Leave empty for equal apportionment.' },
  { key: 'tMission', label: 'Mission time', unit: 'h', value: 10000, min: 1, max: 1e6, log: true, group: 'Reliability', showIf: R('alloc'), help: 'Operating time over which the reliability target applies; used to convert item reliabilities into constant failure rates.' },
  { key: 'Ritem', label: 'Achieved item reliability', value: 0.995, min: 0.5, max: 1, group: 'Reliability', showIf: R('alloc'), help: 'Reliability each subsystem is currently believed to reach, to compare the resulting system reliability with the target.' },
  // Buckingham Pi
  { key: 'piProb', label: 'Dimensional-analysis problem', type: 'select', value: 'reactor', group: 'Dimensional analysis', showIf: R('pi'), options: Object.entries(PIP).map(([value, p]) => ({ value, label: p.label })), help: 'Built-in list of variables whose dimensional matrix is analysed. The values entered below are used to evaluate the resulting groups.' },
  { key: 'pU', label: 'Superficial velocity U', unit: 'm/s', value: 1e-3, min: 1e-7, max: 10, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, help: 'Volumetric flow divided by the reactor cross-section.' },
  { key: 'pL', label: 'Bed / reactor length L', unit: 'm', value: 0.1, min: 1e-4, max: 100, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, help: 'Length of the reacting zone in the flow direction.' },
  { key: 'pd', label: 'Pore / particle size d', unit: 'm', value: 1e-4, min: 1e-9, max: 0.1, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, bus: 'geom.diameter', help: 'Characteristic pore or particle dimension of the porous medium.' },
  { key: 'pRho', label: 'Fluid density ρ', unit: 'kg/m³', value: 998, min: 0.05, max: 2e4, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, help: 'Density of the reacting fluid (water at 20 °C by default).' },
  { key: 'pMu', label: 'Dynamic viscosity μ', unit: 'Pa·s', value: 1e-3, min: 1e-6, max: 100, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, help: 'Dynamic viscosity of the fluid.' },
  { key: 'pD', label: 'Diffusivity D', unit: 'm²/s', value: 1e-9, min: 1e-14, max: 1e-3, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, bus: 'prop.D_diff', help: 'Molecular diffusivity of the limiting species (about 10⁻⁹ m²/s for ions in water).' },
  { key: 'pK', label: 'First-order rate constant k', unit: '1/s', value: 1e-2, min: 1e-8, max: 1e4, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['reactor'] }, help: 'Pseudo-first-order rate constant of the limiting species.' },
  { key: 'pKth', label: 'Thermal conductivity k', unit: 'W/m/K', value: 50, min: 0.005, max: 3000, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, bus: 'prop.k_thermal', help: 'Thermal conductivity of the body being heated.' },
  { key: 'pLh', label: 'Characteristic length L (V/A)', unit: 'm', value: 1e-3, min: 1e-9, max: 1, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, help: 'Volume-to-surface ratio of the body (radius/3 for a sphere, half thickness for a plate).' },
  { key: 'pRhoS', label: 'Solid density ρ', unit: 'kg/m³', value: 5606, min: 1, max: 25000, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, bus: 'prop.density', help: 'Density of the body being heated (ZnO by default).' },
  { key: 'pCp', label: 'Specific heat c_p', unit: 'J/kg/K', value: 494, min: 50, max: 5000, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, bus: 'prop.cp', help: 'Specific heat capacity of the body.' },
  { key: 'pH', label: 'Heat-transfer coefficient h', unit: 'W/m²/K', value: 50, min: 0.1, max: 1e5, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, help: 'Convective coefficient at the surface: about 5–25 for still air, 50–1000 for forced liquids.' },
  { key: 'pT', label: 'Time t', unit: 's', value: 10, min: 1e-6, max: 1e6, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['heat'] }, bus: 'proc.time', help: 'Elapsed heating time at which the Fourier number is evaluated.' },
  { key: 'pM', label: 'Resonator mass m', unit: 'kg', value: 0.01, min: 1e-15, max: 1e3, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['resonator'] }, help: 'Moving mass of the local resonator.' },
  { key: 'pKs', label: 'Spring stiffness k', unit: 'N/m', value: 4000, min: 1e-6, max: 1e10, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['resonator'] }, help: 'Stiffness connecting the resonator mass to the host.' },
  { key: 'pC', label: 'Damping coefficient c', unit: 'N·s/m', value: 0.5, min: 0, max: 1e5, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['resonator'] }, help: 'Viscous damping of the resonator; it limits the peak response and the depth of the band.' },
  { key: 'pF', label: 'Drive frequency f', unit: 'Hz', value: 120, min: 1e-3, max: 1e12, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['resonator'] }, bus: 'res.f0', help: 'Frequency of the harmonic excitation.' },
  { key: 'pF0', label: 'Force amplitude F₀', unit: 'N', value: 1, min: 1e-12, max: 1e6, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['resonator'] }, help: 'Amplitude of the harmonic force acting on the resonator.' },
  { key: 'pEs', label: 'Solid modulus E_s', unit: 'Pa', value: 2e9, min: 1e5, max: 1.2e12, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, bus: 'prop.E', help: 'Young’s modulus of the cell-wall material.' },
  { key: 'pSy', label: 'Solid yield strength σ_ys', unit: 'Pa', value: 5e7, min: 1e4, max: 1e10, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, bus: 'prop.yield', help: 'Yield strength of the cell-wall material.' },
  { key: 'pTs', label: 'Wall thickness t', unit: 'm', value: 5e-4, min: 1e-9, max: 0.1, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, bus: 'geom.strut', help: 'Thickness of the cell walls or struts.' },
  { key: 'pLc', label: 'Cell wall length l', unit: 'm', value: 5e-3, min: 1e-8, max: 1, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, bus: 'geom.cell_size', help: 'Length of the inclined cell wall (unit-cell scale).' },
  { key: 'pSig', label: 'Applied stress σ', unit: 'Pa', value: 1e4, min: 1, max: 1e10, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, help: 'Macroscopic stress applied to the lattice.' },
  { key: 'pW', label: 'Specimen width W', unit: 'm', value: 0.05, min: 1e-6, max: 10, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, help: 'Overall width of the specimen or component, to count cells across it.' },
  { key: 'pEeff', label: 'Effective modulus E*', unit: 'Pa', value: 2e6, min: 1, max: 1e12, log: true, group: 'Dimensional analysis', showIf: { mode: ['pi'], piProb: ['lattice'] }, bus: 'eff.E', help: 'Measured or computed effective modulus of the lattice, to compare with the scaling laws.' },
  // AHP
  { key: 'ahpNames', label: 'Criteria names', type: 'text', value: 'performance index, low density, low cost, low CO₂', group: 'AHP', showIf: R('ahp'), help: 'Comma-separated names of the criteria, in the order of the rows of the comparison matrix.' },
  { key: 'ahpMatrix', label: 'Pairwise-comparison matrix', type: 'text', value: AHP_DEF, group: 'AHP', showIf: R('ahp'), help: 'Rows separated by semicolons or new lines, entries by commas; fractions such as 1/3 are allowed. Entry (i, j) states how much more important criterion i is than j on the 1–9 scale.' },
];

export default {
  id: 1,
  tagline: 'Turns functional requirements into screened and ranked material candidates, allocated tolerances and margins, dimensionless groups, criteria weights and a traceable preliminary specification.',
  inputs,
  presets: [
    { name: 'ZnO nanorod UV sensor: semiconductor selection', note: 'Wide-gap semiconductor nanomaterial for a UV photodetector: band gap 3.0–3.6 eV near 3.3 eV, moderate service temperature.', values: { mode: 'rank', func: 'semi', kind: 'nanomaterial', EgTarget: 3.3, EgMin: 3.0, EgMax: 3.6, rhoMax: 25000, Emin: 0, Top: 150, wPerf: 0.6, wRho: 0.1, wCost: 0.1, wCo2: 0.1, wTemp: 0.1, objX: 'density', objY: 'index' } },
    { name: 'Auxetic energy-absorbing lattice', note: 'Negative Poisson’s ratio below −0.5 for an impact-protection panel; candidates restricted to the auxetic class.', values: { mode: 'rank', func: 'auxetic', kind: 'all', nuMax: -0.5, rhoMax: 25000, Emin: 0, Top: 60, objX: 'density', objY: 'index' } },
    { name: 'Negative-index / metasurface optics', note: 'Optical metasurface and negative-index candidates ranked by index magnitude; no mechanical limits.', values: { mode: 'rank', func: 'optical', kind: 'all', rhoMax: 25000, Emin: 0, Top: 80, wPerf: 0.7, wRho: 0, wCost: 0.1, wCo2: 0.1, wTemp: 0.1, objX: 'Tm', objY: 'index' } },
    { name: 'Zero-expansion optical bench', note: 'Near-zero or negative thermal expansion below 1.5 × 10⁻⁶ /K at up to 100 °C.', values: { mode: 'rank', func: 'nte', cteMax: 1.5, rhoMax: 25000, Emin: 0, Top: 100, objX: 'cte', objY: 'Tm' } },
    { name: 'Lattice panel: tolerance stack and safety margin', note: 'Four-cell printed lattice on a base plate: stack-up against ±0.1 mm, strength margin with derating and reliability split over five subsystems.', values: { mode: 'alloc' } },
    { name: 'Hydrothermal flow reactor: dimensionless groups', note: 'Porous-bed synthesis reactor: Reynolds, Schmidt, Péclet, Damköhler and Thiele numbers from seven variables.', values: { mode: 'pi', piProb: 'reactor' } },
    { name: 'Mass-in-mass resonator similarity', note: 'Local resonator of an acoustic metamaterial driven just above resonance, the negative-effective-mass regime.', values: { mode: 'pi', piProb: 'resonator', pM: 0.01, pKs: 4000, pC: 0.2, pF: 110, pF0: 1 } },
    { name: 'Specification sheet: ZnO nanorod sensor', note: 'Preliminary specification with traceable requirements for hydrothermal ZnO nanorods.', values: { mode: 'spec', func: 'semi', kind: 'nanomaterial', EgTarget: 3.3, EgMin: 3.0, EgMax: 3.6, rhoMax: 25000, Emin: 0, Top: 150, project: 'ZnO nanorod UV sensor', specId: 'zno-nanorods' } },
  ],
  run(inp, ctx) {
    const r = ({ rank: runRank, alloc: runAlloc, pi: runPi, ahp: runAhp, spec: runSpec }[inp.mode] || runRank)(inp, ctx);
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    return r;
  },
};
