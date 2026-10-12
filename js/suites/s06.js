// Suite 6 — Thermodynamics, chemical equilibrium and phase behaviour.
// Aqueous hydrolysis speciation and oxide/hydroxide solubility, reaction equilibrium, regular-solution
// phase diagrams, solvent vapour pressure in sealed vessels and nanoscale (capillarity) thermodynamics.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS, LN10 = Math.LN10, T25 = 298.15;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };

// Hydrolysis constants log *β_n for M^z+ + n H2O ⇌ M(OH)_n^(z-n) + n H+ and solid dissolution log *Ks0 for
// solid + z H+ ⇌ M^z+ + (water), 25 °C, infinite dilution (typical critically-compiled values, Baes–Mesmer type).
// dHs: dissolution enthalpy (J/mol) from standard formation enthalpies; null = not applied.
const METALS = {
  Zn: { ion: 'Zn2+', z: 2, solid: 'ZnO', solidName: 'zincite, ZnO', logKs: 11.14, dHs: -89.3e3, logB: [-8.96, -16.9, -28.4, -41.2], salt: 'Zn(NO3)2', by: 'NaNO3', water: true },
  Cu: { ion: 'Cu2+', z: 2, solid: 'CuO', solidName: 'tenorite, CuO', logKs: 7.62, dHs: -63.7e3, logB: [-8.0, -17.3, -27.8, -39.6], salt: 'Cu(NO3)2', by: 'NaNO3', water: true },
  Fe: { ion: 'Fe3+', z: 3, solid: 'Fe(OH)3', solidName: 'ferrihydrite, Fe(OH)₃(am)', logKs: 4.89, dHs: null, logB: [-2.19, -5.67, -12.56, -21.6], salt: 'Fe(NO3)3', by: 'NaNO3', water: false },
  Al: { ion: 'Al3+', z: 3, solid: 'Al(OH)3', solidName: 'gibbsite, Al(OH)₃', logKs: 8.11, dHs: -95.4e3, logB: [-4.97, -9.3, -15.0, -23.0], salt: 'Al(NO3)3', by: 'NaNO3', water: false },
  Ni: { ion: 'Ni2+', z: 2, solid: 'Ni(OH)2', solidName: 'β-Ni(OH)₂', logKs: 10.8, dHs: -96.0e3, logB: [-9.86, -19.0, -30.0, -44.0], salt: 'Ni(NO3)2', by: 'NaNO3', water: false },
  Mg: { ion: 'Mg2+', z: 2, solid: 'Mg(OH)2', solidName: 'brucite, Mg(OH)₂', logKs: 16.84, dHs: -113.4e3, logB: [-11.44], salt: 'MgCl2', by: 'NaCl', water: false },
};
const DH_W = 55.8e3;   // enthalpy of water autoionisation, J/mol
const el = m => m.ion.replace(/\d*\+$/, '');
const complexFormula = (m, n) => { const q = m.z - n, core = n === 1 ? `${el(m)}OH` : `${el(m)}(OH)${n}`; if (q === 0) return core; return `[${core}]${Math.abs(q) > 1 ? Math.abs(q) : ''}${q > 0 ? '+' : '-'}`; };
const vantHoff = (logK, dH, T) => logK - (dH / (R * LN10)) * (1 / T - 1 / T25);
const epsWater = tC => 87.74 - 0.4008 * tC + 9.398e-4 * tC * tC - 1.41e-6 * tC ** 3;
const davies = (q, I, T) => { const A = 1.8246e6 / Math.pow(Math.max(epsWater(T - 273.15), 5) * T, 1.5), s = Math.sqrt(I); return -A * q * q * (s / (1 + s) - 0.3 * I); };   // log10 γ

function aqueousModel(m, T, I) {
  const pKw = -vantHoff(-14.0, DH_W, T);
  const logKs = m.dHs === null ? m.logKs : vantHoff(m.logKs, m.dHs, T);
  const logB = m.logB.map((b, i) => vantHoff(b, (i + 1) * DH_W, T));
  const lg = q => davies(q, I, T);
  // state at a given pH: species fractions, log10 equilibrium solubility (mol/L)
  const at = pH => {
    const lr = logB.map((b, i) => b + (i + 1) * pH + lg(m.z) - lg(m.z - (i + 1)));   // log([M(OH)n]/[M])
    const mx = Math.max(0, ...lr); let D = Math.pow(10, -mx); for (const v of lr) D += Math.pow(10, v - mx);
    const logD = mx + Math.log10(D), alpha = [Math.pow(10, -logD), ...lr.map(v => Math.pow(10, v - logD))];
    return { alpha, logS: logKs - m.z * pH - lg(m.z) + logD };
  };
  return { pKw, logKs, logB, lg, at };
}

const SOLV = {
  water: { name: 'water', f: 'H2O', ant: [[8.07131, 1730.63, 233.426, 1, 100], [8.14019, 1810.94, 244.485, 100, 374]], Tc: 647.096, Pc: 22.064e6, w: 0.3443, rho0: 997, kap: 4.5e-10 },
  ethanol: { name: 'ethanol', f: 'C2H5OH', ant: [[8.20417, 1642.89, 230.3, -57, 80], [7.68117, 1332.04, 199.2, 80, 240]], Tc: 514.0, Pc: 6.14e6, w: 0.645, rho0: 785, kap: 1.1e-9 },
  methanol: { name: 'methanol', f: 'CH3OH', ant: [[8.08097, 1582.271, 239.726, 15, 100]], Tc: 512.6, Pc: 8.08e6, w: 0.565, rho0: 787, kap: 1.2e-9 },
  eg: { name: 'ethylene glycol', f: 'C2H6O2', ant: [[8.09083, 2088.936, 203.454, 50, 200]], Tc: 720, Pc: 8.0e6, w: 0.507, rho0: 1110, kap: 3.7e-10 },
  // DMF: the tabulated low-temperature fit (30–90 °C) extrapolates to a boiling point 3 K too low, so the upper range is an
  // Antoine curve through the 90 °C value of that fit, the normal boiling point (153.0 °C) and the critical point.
  dmf: { name: 'N,N-dimethylformamide', f: 'C3H7NO', ant: [[6.92796, 1400.87, 196.43, 30, 90], [7.53108, 1907.83, 257.262, 90, 153]], Tc: 649.6, Pc: 4.42e6, w: 0.318, rho0: 944, kap: 6.5e-10 },
  toluene: { name: 'toluene', f: 'C7H8', ant: [[6.95464, 1344.8, 219.482, 6, 137]], Tc: 591.75, Pc: 4.108e6, w: 0.264, rho0: 862, kap: 9e-10 },
};
const antSet = (s, tC) => s.ant.find(a => tC <= a[4]) || s.ant[s.ant.length - 1];
const psat = (s, tC) => { const t = Math.min(tC, s.Tc - 273.15), a = antSet(s, t); return 133.322 * Math.pow(10, a[0] - a[1] / (a[2] + t)); };   // Pa
function rhoLiq(key, T) {   // saturated-liquid density, kg/m³
  const s = SOLV[key], Tt = Math.min(T, 0.995 * s.Tc);
  if (key === 'water') { const t = 1 - Tt / s.Tc; return 322 * (1 + 1.99274064 * t ** (1 / 3) + 1.09965342 * t ** (2 / 3) - 0.510839303 * t ** (5 / 3) - 1.75493479 * t ** (16 / 3) - 45.5170352 * t ** (43 / 3) - 6.74694450e5 * t ** (110 / 3)); }
  const Z = 0.29056 - 0.08775 * s.w, e = x => Math.pow(1 - x / s.Tc, 2 / 7);
  return s.rho0 * Math.pow(Z, e(T25) - e(Tt));
}

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'speciation', group: 'Model', options: [
    { value: 'speciation', label: 'Aqueous speciation and oxide/hydroxide solubility' }, { value: 'reaction', label: 'Reaction equilibrium (van\'t Hoff)' },
    { value: 'regular', label: 'Binary regular-solution phase diagram' }, { value: 'vapour', label: 'Solvent vapour pressure / autoclave' },
    { value: 'nano', label: 'Nanoscale thermodynamics (capillarity, nucleation)' }], help: 'Selects the thermodynamic model. Each analysis uses its own inputs and publishes its own results to the shared bus.' },
  { key: 'T', label: 'Temperature', unit: '°C', value: 25, min: 0, max: 400, group: 'Conditions', bus: 'proc.T_C', showIf: { mode: ['speciation', 'vapour', 'nano'] }, help: 'Solution or vessel temperature. Equilibrium constants are shifted from 25 °C with the van\'t Hoff equation; vapour pressure follows the Antoine equation.' },
  { key: 'metal', label: 'Metal cation', type: 'select', value: 'Zn', group: 'Chemistry', showIf: { mode: ['speciation'] }, options: [{ value: 'Zn', label: 'Zn²⁺ / ZnO (zincite)' }, { value: 'Cu', label: 'Cu²⁺ / CuO (tenorite)' }, { value: 'Fe', label: 'Fe³⁺ / Fe(OH)₃ (ferrihydrite)' }, { value: 'Al', label: 'Al³⁺ / Al(OH)₃ (gibbsite)' }, { value: 'Ni', label: 'Ni²⁺ / Ni(OH)₂' }, { value: 'Mg', label: 'Mg²⁺ / Mg(OH)₂ (brucite)' }], help: 'Metal whose hydrolysis species and least-soluble oxide or hydroxide are modelled.' },
  { key: 'cM', label: 'Total metal (or reactant A) concentration', unit: 'mol/m³', value: 10, min: 1e-4, max: 5000, log: true, group: 'Chemistry', bus: 'proc.c_precursor', showIf: { mode: ['speciation', 'reaction'] }, help: 'Total dissolved concentration before precipitation (10 mol/m³ = 0.01 mol/L). It sets the saturation index and the pH window in which the solid can precipitate.' },
  { key: 'pHsrc', label: 'pH definition', type: 'select', value: 'given', group: 'Chemistry', showIf: { mode: ['speciation'] }, options: [{ value: 'given', label: 'pH specified' }, { value: 'base', label: 'pH computed from strong base added' }], help: 'Either fix the pH, or let the charge balance compute it from the amount of strong base (e.g. NaOH) added to the metal-salt solution.' },
  { key: 'pH', label: 'Solution pH', unit: '', value: 10.5, min: 0, max: 14, step: 0.1, group: 'Chemistry', bus: 'proc.pH', showIf: { mode: ['speciation'] }, help: 'pH (activity scale) at which speciation, solubility and supersaturation are evaluated when the pH is specified.' },
  { key: 'cBase', label: 'Strong base added (OH⁻)', unit: 'mol/m³', value: 20, min: 0, max: 10000, group: 'Chemistry', bus: 'proc.c_base', showIf: { mode: ['speciation'] }, help: 'Hydroxide added as a fully dissociated base. Used when the pH is computed from the charge balance; the anion of the metal salt is taken as non-complexing.' },
  { key: 'I', label: 'Ionic strength', unit: 'mol/L', value: 0.05, min: 0, max: 1, step: 0.01, group: 'Chemistry', showIf: { mode: ['speciation'] }, help: 'Used in the Davies activity-coefficient equation. The Davies model is reliable up to about 0.5 mol/L.' },
  { key: 'Tr', label: 'Temperature', unit: 'K', value: 350, min: 200, max: 3000, group: 'Conditions', showIf: { mode: ['reaction', 'regular'] }, help: 'Absolute temperature for the reaction-equilibrium and regular-solution analyses.' },
  { key: 'dH', label: 'Standard reaction enthalpy ΔH°', unit: 'kJ/mol', value: -50, min: -1000, max: 1000, group: 'Reaction', showIf: { mode: ['reaction'] }, help: 'Enthalpy change at 298.15 K. Negative values (exothermic) make the equilibrium constant fall as temperature rises.' },
  { key: 'dS', label: 'Standard reaction entropy ΔS°', unit: 'J/mol/K', value: -100, min: -1000, max: 1000, group: 'Reaction', showIf: { mode: ['reaction'] }, help: 'Entropy change at 298.15 K (1 mol/L standard state for solutes). Together with ΔH° it fixes ΔG° and K at every temperature.' },
  { key: 'dCp', label: 'Reaction heat-capacity change ΔCp', unit: 'J/mol/K', value: 0, min: -500, max: 500, group: 'Reaction', showIf: { mode: ['reaction'] }, help: 'Kirchhoff correction of ΔH° and ΔS° with temperature. Leave at zero for the plain van\'t Hoff form.' },
  { key: 'bst', label: 'Stoichiometric coefficient b of B', unit: '', value: 1, min: 0, max: 4, step: 1, group: 'Reaction', showIf: { mode: ['reaction'] }, help: 'Reaction A + b B ⇌ C. Zero gives an isomerisation A ⇌ C.' },
  { key: 'ratioB', label: 'Feed ratio B/A', unit: 'mol/mol', value: 2, min: 0, max: 100, group: 'Reaction', showIf: { mode: ['reaction'] }, help: 'Initial moles of B per mole of A. Excess B pushes the equilibrium conversion of A upwards.' },
  { key: 'omega', label: 'Interaction parameter Ω', unit: 'kJ/mol', value: 18, min: -20, max: 60, group: 'Phase behaviour', showIf: { mode: ['regular'] }, help: 'Regular-solution exchange energy. Positive values favour unmixing with a critical temperature Ω/2R; negative values favour ordering and no miscibility gap.' },
  { key: 'x0', label: 'Overall mole fraction of B', unit: '', value: 0.3, min: 0.001, max: 0.999, step: 0.01, group: 'Phase behaviour', showIf: { mode: ['regular'] }, help: 'Overall composition used for the lever rule and for the stability classification.' },
  { key: 'solvent', label: 'Solvent', type: 'select', value: 'water', group: 'Vessel', showIf: { mode: ['vapour'] }, options: [{ value: 'water', label: 'Water' }, { value: 'ethanol', label: 'Ethanol' }, { value: 'methanol', label: 'Methanol' }, { value: 'eg', label: 'Ethylene glycol' }, { value: 'dmf', label: 'DMF' }, { value: 'toluene', label: 'Toluene' }], help: 'Solvent sealed in the vessel. Antoine constants, critical point and liquid expansion are solvent specific.' },
  { key: 'fill', label: 'Liquid fill fraction at 25 °C', unit: '%', value: 70, min: 5, max: 99, group: 'Vessel', showIf: { mode: ['vapour'] }, help: 'Fraction of the vessel volume occupied by liquid when sealed. Liquids expand on heating; a vessel that becomes liquid-full sees a very steep (hydraulic) pressure rise.' },
  { key: 'Prate', label: 'Vessel pressure rating', unit: 'MPa', value: 6, min: 0.1, max: 100, group: 'Vessel', sens: false, showIf: { mode: ['vapour'] }, help: 'Maximum allowable working pressure of the autoclave, used for the safety margin.' },
  { key: 'Pamb', label: 'Ambient / sealing pressure', unit: 'Pa', value: 101325, min: 1000, max: 1e6, group: 'Vessel', bus: 'proc.P', showIf: { mode: ['vapour'] }, help: 'Pressure at which the vessel is sealed; also used for the normal boiling point.' },
  { key: 'gamma', label: 'Solid–liquid interfacial energy', unit: 'J/m²', value: 0.25, min: 0.005, max: 3, group: 'Nanoscale', bus: 'prop.surface_energy', showIf: { mode: ['nano'] }, help: 'Interfacial free energy of the particle against its solution or melt. It enters the nucleation barrier to the third power.' },
  { key: 'Mw', label: 'Molar mass of the solid', unit: 'g/mol', value: 81.38, min: 10, max: 1000, group: 'Nanoscale', showIf: { mode: ['nano'] }, help: 'Molar mass per formula unit (81.38 g/mol for ZnO); with density it gives the molar volume.' },
  { key: 'rho', label: 'Solid density', unit: 'kg/m³', value: 5606, min: 500, max: 25000, group: 'Nanoscale', bus: 'prop.density', showIf: { mode: ['nano'] }, help: 'Bulk density of the solid phase.' },
  { key: 'dnm', label: 'Particle diameter', unit: 'nm', value: 10, min: 0.5, max: 1000, log: true, group: 'Nanoscale', showIf: { mode: ['nano'] }, help: 'Diameter at which the size-dependent melting point and solubility are evaluated.' },
  { key: 'S', label: 'Supersaturation ratio', unit: '', value: 50, min: 1.01, max: 1e8, log: true, group: 'Nanoscale', bus: 'proc.supersaturation', showIf: { mode: ['nano'] }, help: 'Ratio of the ion activity product (or concentration) to its equilibrium value. Drives nucleation through ln S.' },
  { key: 'Tm', label: 'Bulk melting temperature', unit: 'K', value: 2248, min: 200, max: 4500, group: 'Nanoscale', showIf: { mode: ['nano'] }, help: 'Melting point of the bulk solid for the Gibbs–Thomson relation.' },
  { key: 'Hf', label: 'Enthalpy of fusion', unit: 'kJ/mol', value: 70, min: 1, max: 500, group: 'Nanoscale', showIf: { mode: ['nano'] }, help: 'Molar latent heat of melting used in the Gibbs–Thomson melting-point depression.' },
  { key: 'theta', label: 'Contact angle of nucleus on substrate', unit: '°', value: 60, min: 1, max: 180, group: 'Nanoscale', showIf: { mode: ['nano'] }, help: 'Wetting angle of the nucleus on a seed or wall. Small angles lower the heterogeneous nucleation barrier strongly; 180° equals homogeneous nucleation.' },
  { key: 'logA', label: 'log₁₀ nucleation prefactor', unit: 'log(1/m³/s)', value: 33, min: 20, max: 42, group: 'Nanoscale', sens: false, showIf: { mode: ['nano'] }, help: 'Kinetic prefactor of classical nucleation theory. It is uncertain by several orders of magnitude; the rate is therefore an order-of-magnitude estimate.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
const num = (inp, k) => { const d = DEF[k]; let v = +inp[k]; if (!Number.isFinite(v)) v = d.value; return N.clamp(v, d.min, d.max); };
const sel = (inp, k) => (DEF[k].options.some(o => o.value === inp[k]) ? inp[k] : DEF[k].value);

// ---------------------------------------------------------------- speciation
function runSpeciation(inp) {
  const mk = sel(inp, 'metal'), m = METALS[mk], T = num(inp, 'T') + 273.15, I = num(inp, 'I'), CT = num(inp, 'cM') / 1000, cB = num(inp, 'cBase') / 1000;
  const src = sel(inp, 'pHsrc'), A = aqueousModel(m, T, I), warnings = [], checks = [];
  const g1 = Math.pow(10, A.lg(1));
  const charge = (pH, equilibrated) => {
    const s = A.at(pH), Cd = equilibrated ? Math.min(CT, Math.pow(10, Math.min(s.logS, 6))) : CT;
    let q = 0; s.alpha.forEach((a, n) => { q += (m.z - n) * a * Cd; });
    return q + Math.pow(10, -pH) / g1 - Math.pow(10, pH - A.pKw) / g1 + cB - m.z * CT;
  };
  const solvePH = eq => safe(() => N.brent(p => charge(p, eq), -1, 15.5, 1e-10), N.bisect(p => -charge(p, eq), -1, 15.5, 1e-9));
  let pH = num(inp, 'pH'), pHeq = pH;
  if (src === 'base') { pH = solvePH(false); pHeq = solvePH(true); }
  const st = A.at(pH), logS = st.logS, Seq = Math.pow(10, Math.min(logS, 6)), SI = Math.log10(CT) - logS, Sr = Math.pow(10, N.clamp(SI, -30, 30));
  const stEq = A.at(pHeq), SeqEq = Math.pow(10, Math.min(stEq.logS, 6)), fPrec = Math.max(0, 1 - SeqEq / CT);
  const names = [m.ion, ...m.logB.map((_, i) => complexFormula(m, i + 1))], pretty = names.map(chem.prettyFormula);
  const dom = N.argmax(st.alpha);
  // scans
  const pHs = N.linspace(0, 14, 281), sc = pHs.map(p => A.at(p)), logSs = sc.map(s => s.logS);
  const iMin = N.argmin(logSs), pHmin = pHs[iMin], Smin = Math.pow(10, logSs[iMin]);
  const lc = Math.log10(CT); let lo = null, hi = null;
  for (let i = 0; i < pHs.length; i++) if (logSs[i] < lc) { if (lo === null) lo = i > 0 ? pHs[i - 1] + (pHs[i] - pHs[i - 1]) * (logSs[i - 1] - lc) / (logSs[i - 1] - logSs[i]) : pHs[0]; hi = i < pHs.length - 1 && logSs[i + 1] >= lc ? pHs[i] + (pHs[i + 1] - pHs[i]) * (lc - logSs[i]) / (logSs[i + 1] - logSs[i]) : pHs[i]; }
  const dmu = R * T * Math.log(Sr);
  // equations
  const equations = [];
  equations.push(eqn('ionic', `Dissolution of ${m.solidName}`, safe(() => chem.equationString([m.solid, 'H+'], [m.ion, 'H2O'], null, '⇌'), `${m.solid} + ${m.z} H⁺ ⇌ ${m.ion} + water`), `log *Ks0 = ${F(A.logKs, 4)} at ${F(T - 273.15, 4)} °C (${m.logKs} at 25 °C)`));
  m.logB.forEach((b, i) => equations.push(eqn('ionic', `Hydrolysis step n = ${i + 1}`, safe(() => chem.equationString([m.ion, 'H2O'], [names[i + 1], 'H+'], null, '⇌'), `${m.ion} + ${i + 1} H₂O ⇌ ${names[i + 1]} + ${i + 1} H⁺`), `log *β${i + 1} = ${F(A.logB[i], 4)} at ${F(T - 273.15, 4)} °C (${b} at 25 °C)`)));
  equations.push(eqn('ionic', 'Water autoionisation', safe(() => chem.equationString(['H2O'], ['H+', 'OH-'], null, '⇌'), 'H₂O ⇌ H⁺ + OH⁻'), `pKw = ${F(A.pKw, 4)}`));
  const ie = safe(() => chem.ionicEquations([{ f: m.salt, phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: m.solid, phase: 's' }, { f: m.by, phase: 'aq' }, ...(m.water ? [{ f: 'H2O', phase: 'l' }] : [])]), null);
  if (ie) equations.push(eqn('chemical', 'Precipitation (molecular)', ie.molecular), eqn('ionic', 'Complete ionic equation', ie.completeIonic), eqn('ionic', 'Net ionic equation', ie.netIonic, `spectators: ${ie.spectators.join(', ')}`));
  equations.push(eqn('governing', 'Mass action with activities', '*βₙ = a(M(OH)ₙ) · a(H⁺)ⁿ / a(Mᶻ⁺);   *Ks0 = a(Mᶻ⁺) / a(H⁺)ᶻ'),
    eqn('governing', 'Solubility and saturation', 'S_eq(pH) = Σₙ [M(OH)ₙ] at saturation;   SI = log₁₀(IAP / Ks0);   S = IAP / Ks0;   Δμ = R T ln S'),
    eqn('constitutive', 'Davies activity coefficients', 'log γᵢ = −A zᵢ² ( √I / (1 + √I) − 0.3 I ),   A = 1.8246×10⁶ (εᵣ T)^(−3/2)'),
    eqn('constitutive', 'van\'t Hoff temperature correction', 'log K(T) = log K(298) − ΔH° / (R ln 10) · (1/T − 1/298.15)'));
  if (src === 'base') equations.push(eqn('governing', 'Charge balance for pH', 'Σₙ (z − n)[M(OH)ₙ] + [H⁺] − [OH⁻] + [Na⁺] − z·C_M(anion) = 0'));
  // checks and warnings
  const sumA = N.sum(st.alpha);
  checks.push(chk('Species fractions sum to one', Math.abs(sumA - 1) < 1e-9, `Σα − 1 = ${F(sumA - 1, 2)}`));
  checks.push(chk('Davies model validity (I ≤ 0.5 mol/L)', I <= 0.5, `ionic strength ${F(I, 3)} mol/L`));
  checks.push(chk('Temperature within van\'t Hoff extrapolation range (≤ 150 °C)', T <= 423.15, `T = ${F(T - 273.15, 4)} °C; constant ΔH° assumed`));
  if (src === 'base') checks.push(chk('Charge balance closed', Math.abs(charge(pH, false)) < 1e-7 * (CT + cB + 1e-7) + 1e-10, `residual ${F(charge(pH, false), 2)} mol/L (metastable solution)`));
  if (I > 0.5) warnings.push('Ionic strength above 0.5 mol/L: the Davies equation is outside its range; a Pitzer or SIT model is needed for quantitative work.');
  if (T > 423.15) warnings.push('Above 150 °C the constant-ΔH° van\'t Hoff extrapolation of the 25 °C constants becomes unreliable (heat-capacity terms neglected).');
  if (m.dHs === null) warnings.push('No dissolution enthalpy applied for the amorphous solid: its solubility constant is held at the 25 °C value.');
  warnings.push('Hydrolysis enthalpies are approximated as n × 55.8 kJ/mol (isocoulombic assumption); polynuclear species, carbonate, ammine and anion complexes are not included. The ionic strength is taken as given, not iterated.');
  const sat = SI > 0;
  const kpis = [
    kpi('pH', src === 'base' ? 'pH of the mixed (metastable) solution' : 'pH evaluated', pH, '', src === 'base' ? `Charge-balance pH with all ${el(m)} still dissolved as hydroxo species, i.e. before any solid forms. After precipitation to equilibrium the pH relaxes to ${F(pHeq, 4)}.` : 'pH on the activity scale at which the speciation and solubility are evaluated, taken from the input.'),
    kpi('S_eq', 'Equilibrium solubility at this pH', Seq * 1000, 'mol/m³', `Sum of all dissolved ${el(m)} hydroxo species in equilibrium with ${m.solidName}; ${F(Seq, 3)} mol/L. Concentrations above this value are supersaturated.`),
    kpi('SI', 'Saturation index', SI, '', 'SI = log₁₀(IAP/Ks0). Positive means the solid can precipitate, negative means it dissolves, zero is equilibrium.'),
    kpi('S', 'Supersaturation ratio', Sr, '', 'Ratio of the ion activity product to the solubility constant at the stated total concentration; the driving force passed to the nucleation and growth suite.'),
    kpi('dmu', 'Thermodynamic driving force Δμ', dmu / 1000, 'kJ/mol', 'Δμ = R T ln S per mole of metal. Positive values drive precipitation; a few kJ/mol is typical of controlled crystal growth, tens of kJ/mol of burst nucleation.'),
    kpi('dominant', 'Dominant dissolved species', `${pretty[dom]} (${F(100 * st.alpha[dom], 3)} %)`, '', 'Hydroxo complex carrying the largest share of dissolved metal at this pH; its charge controls adsorption on growing crystal faces.'),
    kpi('pH_min', 'pH of minimum solubility', pHmin, '', 'pH at which the amphoteric solubility curve passes through its minimum; precipitation yield is highest and dissolution losses lowest here.'),
    kpi('S_min', 'Minimum solubility', Smin * 1000, 'mol/m³', `Lowest total dissolved concentration reachable by pH adjustment alone, ${F(Smin, 3)} mol/L; it limits the attainable yield.`),
    kpi('win_lo', 'Precipitation window, lower pH', lo === null ? 'none' : lo, '', 'Lowest pH at which the stated total concentration exceeds the solubility of the solid. Below it the solution is undersaturated.'),
    kpi('win_hi', 'Precipitation window, upper pH', hi === null ? 'none' : hi, '', 'Highest pH at which the solid is still supersaturated; above it soluble hydroxo anions redissolve the solid (14 means it stays insoluble to the end of the scan).'),
    kpi('f_prec', 'Fraction precipitated at equilibrium', fPrec, '', `Share of the metal that leaves solution if the system reaches equilibrium with ${m.solidName} at pH ${F(pHeq, 4)}; the thermodynamic yield limit.`, fPrec > 0.9 ? true : null),
    kpi('gammaM', `Activity coefficient of ${pretty[0]}`, Math.pow(10, A.lg(m.z)), '', 'Davies activity coefficient of the free cation at the stated ionic strength; values well below one show that concentrations and activities differ appreciably.'),
    kpi('pKw', 'pKw at temperature', A.pKw, '', 'Ion product of water from the van\'t Hoff equation; neutral pH is pKw/2 and falls as temperature rises.'),
  ];
  const cl = v => N.clamp(v, -12, 2);
  const series = [
    { key: 'dist', title: `${el(m)} species distribution vs pH`, xlabel: 'pH', ylabel: 'fraction of dissolved metal', lines: names.map((_, n) => ({ name: pretty[n], x: pHs, y: sc.map(s => s.alpha[n]) })), explain: 'Fraction of dissolved metal present as each hydroxo complex. Crossing points lie at the stepwise hydrolysis pK values; the neutral complex dominates near the solubility minimum.' },
    { key: 'sol', title: `Solubility of ${m.solidName}`, xlabel: 'pH', ylabel: 'log₁₀ concentration (mol/L)', lines: [{ name: 'equilibrium solubility', x: pHs, y: logSs.map(cl) }, { name: 'total metal in solution', x: [0, 14], y: [cl(lc), cl(lc)] }, { name: 'operating point', x: [pH, pH], y: [cl(logS), cl(lc)] }], explain: 'Amphoteric solubility curve (clipped to the range 10⁻¹² to 10² mol/L). Where the horizontal total-concentration line lies above the curve the solid is supersaturated.' },
    { key: 'si', title: 'Saturation index vs pH', xlabel: 'pH', ylabel: 'SI', lines: [{ name: 'SI', x: pHs, y: logSs.map(v => N.clamp(lc - v, -15, 15)) }, { name: 'equilibrium', x: [0, 14], y: [0, 0] }], explain: 'Saturation index of the solid at the stated total concentration (clipped to ±15). The positive region is the precipitation window.' },
  ];
  const tables = [{ key: 'eq', title: 'Equilibria and constants', columns: ['equilibrium', 'log K at 25 °C', `log K at ${F(T - 273.15, 4)} °C`], rows: [[`${chem.prettyFormula(m.solid)} dissolution`, m.logKs, +A.logKs.toFixed(3)], ...m.logB.map((b, i) => [`→ ${pretty[i + 1]}`, b, +A.logB[i].toFixed(3)]), ['water autoionisation', -14, +(-A.pKw).toFixed(3)]], explain: 'Constants used in this run, written for reactions with H⁺ (cumulative hydrolysis constants and acid dissolution constant).' },
    { key: 'spec', title: 'Speciation at the operating pH', columns: ['species', 'fraction', 'concentration (mol/m³)', 'log γ'], rows: names.map((_, n) => [pretty[n], +st.alpha[n].toExponential(3), +(st.alpha[n] * CT * 1000).toExponential(3), +A.lg(m.z - n).toFixed(3)]), explain: 'Distribution of the total dissolved metal over hydroxo species before any solid forms, with the Davies activity coefficients used.' }];
  const recommendations = [
    rec(sat ? 'high' : 'medium', sat ? `The solution is supersaturated by a factor ${F(Sr, 3)} (SI ${F(SI, 3)}). ${Sr > 1000 ? 'At this level homogeneous burst nucleation in the bulk competes strongly with seeded growth; lower the supersaturation (slow base release, e.g. HMTA or urea hydrolysis, or a complexing agent) if oriented growth on seeds is wanted.' : 'This is a moderate driving force suited to seeded growth; hold it steady by controlled base addition.'}` : `The solution is undersaturated (SI ${F(SI, 3)}); no ${chem.prettyFormula(m.solid)} can form. Move the pH into the window ${lo === null ? '(none exists at this concentration — raise the concentration)' : `${F(lo, 3)}–${F(hi, 3)}`}.`, 'Saturation index and supersaturation ratio from the mass-action model.'),
    rec('medium', `Operate near pH ${F(pHmin, 3)} to maximise yield: the minimum solubility is ${F(Smin, 3)} mol/L, so at ${F(CT, 3)} mol/L total metal at most ${F(100 * Math.max(0, 1 - Smin / CT), 4)} % can be recovered as solid.`, 'Minimum of the amphoteric solubility curve.'),
    rec('medium', `Keep the pH below ${hi === null ? 'the redissolution limit' : F(hi, 3)} during washing and ageing; hydroxo anions (${pretty[pretty.length - 1]}) redissolve the product at high pH.`, 'Upper edge of the precipitation window.'),
    rec('low', `Residual dissolved ${el(m)} in the mother liquor is about ${F(SeqEq * 1000, 3)} mol/m³ at the final pH ${F(pHeq, 3)}. Adjusting the spent bath to the minimum-solubility pH before discharge recovers the metal as solid and lowers the metal load of the effluent.`, 'Equilibrium solubility at the end-point pH.', true),
    rec('low', 'Confirm the predicted window with a titration of the actual bath: ligands (ammonia, citrate, acetate, amines) and carbonate from air shift both solubility and speciation.', 'Model omits complexing ligands and polynuclear species.'),
  ];
  return {
    summary: `At pH ${F(pH, 4)} and ${F(T - 273.15, 4)} °C the solubility of ${m.solidName} is ${F(Seq, 3)} mol/L, so a total ${el(m)} concentration of ${F(CT, 3)} mol/L has a saturation index of ${F(SI, 3)} (supersaturation ratio ${F(Sr, 3)}). The dominant dissolved species is ${pretty[dom]}. The solid is ${lo === null ? 'undersaturated at every pH for this concentration' : `supersaturated between pH ${F(lo, 3)} and ${F(hi, 3)}`}, with minimum solubility ${F(Smin, 3)} mol/L at pH ${F(pHmin, 3)}.`,
    kpis, series, tables, equations, checks, warnings, recommendations,
    publish: { 'proc.supersaturation': Math.max(Sr, 1e-30), 'proc.pH': pH, 'thermo.solubility': Seq * 1000, 'thermo.SI': SI, 'thermo.driving_force': dmu, 'thermo.pH_min_solubility': pHmin },
    numerics: { method: 'Closed-form mass-action speciation with Davies activities; Brent root finding on the charge balance for pH', dof: names.length + 1, steps: pHs.length, iterations: src === 'base' ? 2 : 0, residual: Math.abs(sumA - 1), notes: 'Ionic strength fixed by the user; constants shifted from 25 °C by the van\'t Hoff equation.' },
  };
}

// ---------------------------------------------------------------- reaction equilibrium
function runReaction(inp) {
  const T = num(inp, 'Tr'), dH0 = num(inp, 'dH') * 1e3, dS0 = num(inp, 'dS'), dCp = num(inp, 'dCp'), b = Math.round(num(inp, 'bst')), th = num(inp, 'ratioB'), cA0 = num(inp, 'cM');
  const c0 = 1000, warnings = [], checks = [];
  const thermo = Tk => { const dH = dH0 + dCp * (Tk - T25), dS = dS0 + dCp * Math.log(Tk / T25), dG = dH - Tk * dS; return { dH, dS, dG, lnK: -dG / (R * Tk) }; };
  const conv = lnK => {
    const Xmax = b > 0 ? Math.min(1, th / b) : 1; if (Xmax <= 0) return 0;
    const f = X => lnK - (Math.log(X) - Math.log(1 - X) - b * Math.log(Math.max((th - b * X) * cA0 / c0, 1e-300)));
    const a = Xmax * 1e-13, z = Xmax * (1 - 1e-13);
    if (f(a) <= 0) return 0; if (f(z) >= 0) return Xmax;
    return safe(() => N.brent(f, a, z, 1e-14), N.bisect(X => -f(X), a, z, 1e-13));
  };
  const s = thermo(T), X = conv(s.lnK), log10K = s.lnK / LN10;
  const Tinv = dH0 * dS0 > 0 && dCp === 0 ? dH0 / dS0 : null;
  const Ts = N.linspace(Math.max(200, T - 150), Math.min(3000, T + 250), 81), th2 = Ts.map(thermo);
  const cA = cA0 * (1 - X), cB = cA0 * Math.max(0, th - b * X), cC = cA0 * X, Q = -s.dH * cA0 * X;
  const Kcheck = (cC / c0) / Math.max((cA / c0) * Math.pow(cB / c0, b), 1e-300), resK = X > 1e-9 && X < 1 - 1e-9 && (b === 0 || cB > 1e-12) ? Math.abs(Math.log(Kcheck) - s.lnK) : 0;
  checks.push(chk('Equilibrium constant reproduced by the computed composition', resK < 1e-5, `|ln Q − ln K| = ${F(resK, 2)}`));
  checks.push(chk('Conversion within stoichiometric bounds', X >= 0 && X <= (b > 0 ? Math.min(1, th / b) : 1) + 1e-12, `X = ${F(X, 4)}`));
  checks.push(chk('Gibbs–Helmholtz consistency', Math.abs((thermo(T + 0.5).lnK - thermo(T - 0.5).lnK) - s.dH / (R * T * T)) < 1e-3 * (Math.abs(s.dH / (R * T * T)) + 1e-9) + 1e-9, `d ln K/dT = ${F(s.dH / (R * T * T), 3)} 1/K (analytic)`));
  if (b > 0 && th < b) warnings.push(`B is the limiting reactant (feed ratio ${F(th, 3)} < b = ${b}); conversion of A cannot exceed ${F(th / b, 3)}.`);
  warnings.push('Ideal-solution model: activity coefficients are taken as one and the 1 mol/L standard state is used for all solutes.');
  const dir = s.dH < 0 ? 'exothermic' : 'endothermic';
  return {
    summary: `At ${F(T, 4)} K the reaction A + ${b} B ⇌ C has ΔG° = ${F(s.dG / 1000, 4)} kJ/mol and K = 10^${F(log10K, 4)}. With A at ${F(cA0, 3)} mol/m³ and a B/A feed ratio of ${F(th, 3)} the equilibrium conversion of A is ${F(100 * X, 4)} %. The reaction is ${dir}, so K ${s.dH < 0 ? 'falls' : 'rises'} with temperature.`,
    kpis: [
      kpi('dG', 'ΔG° at temperature', s.dG / 1000, 'kJ/mol', 'Standard Gibbs energy change ΔG° = ΔH° − TΔS° (with Kirchhoff correction if ΔCp ≠ 0). Negative values favour products at unit activities.', s.dG < 0),
      kpi('log10K', 'log₁₀ K', log10K, '', 'Equilibrium constant from ΔG° = −RT ln K, referred to 1 mol/L. Each unit of log K is one decade of product-to-reactant ratio.'),
      kpi('X', 'Equilibrium conversion of A', X, '', 'Fraction of A converted when the reaction quotient equals K; found by solving the mass-action equation with Brent\'s method. Kinetics may prevent this limit being reached.', X > 0.9 ? true : null),
      kpi('cC', 'Equilibrium product concentration', cC, 'mol/m³', 'Concentration of product C at equilibrium for the stated feed; the theoretical yield per unit volume.'),
      kpi('cB', 'Residual B concentration', cB, 'mol/m³', 'Unreacted co-reactant B at equilibrium; shows how much excess reagent ends up in the waste stream.'),
      kpi('Tinv', 'Inversion temperature (ΔG° = 0)', Tinv && Tinv > 0 ? Tinv : 'none', 'K', 'Temperature ΔH°/ΔS° at which K = 1. It exists only when ΔH° and ΔS° have the same sign (reported for ΔCp = 0).'),
      kpi('Q', 'Heat released at equilibrium', Q / 1e6, 'MJ/m³', 'Reaction heat per cubic metre of solution, −ΔH·c_A0·X. Positive means heat must be removed; negative means heat must be supplied.'),
      kpi('dlnK', 'Temperature sensitivity d ln K/dT', s.dH / (R * T * T), '1/K', 'van\'t Hoff slope ΔH°/(RT²): the fractional change of K per kelvin at the operating temperature.'),
    ],
    series: [
      { key: 'vh', title: 'van\'t Hoff plot', xlabel: '1000/T (1/K)', ylabel: 'ln K', lines: [{ name: 'ln K', x: Ts.map(t => 1000 / t), y: th2.map(v => v.lnK) }], explain: 'ln K against reciprocal temperature; the slope is −ΔH°/R and the intercept ΔS°/R. Curvature appears only when ΔCp ≠ 0.' },
      { key: 'xT', title: 'Equilibrium conversion vs temperature', xlabel: 'T (K)', ylabel: 'conversion of A', lines: [{ name: 'X_eq', x: Ts, y: th2.map(v => conv(v.lnK)) }], explain: 'Thermodynamic conversion limit across temperature for the stated feed. Exothermic reactions lose conversion on heating although they run faster.' },
      { key: 'gT', title: 'ΔG° vs temperature', xlabel: 'T (K)', ylabel: 'ΔG° (kJ/mol)', lines: [{ name: 'ΔG°', x: Ts, y: th2.map(v => v.dG / 1000) }, { name: 'zero', x: [Ts[0], Ts[Ts.length - 1]], y: [0, 0] }], explain: 'Standard Gibbs energy of reaction; where it crosses zero the equilibrium constant equals one.' },
    ],
    equations: [
      eqn('chemical', 'Reaction', `A + ${b} B ⇌ C`, 'generic stoichiometry set by the input b'),
      eqn('governing', 'Gibbs energy and equilibrium constant', 'ΔG°(T) = ΔH°(T) − T ΔS°(T) = −R T ln K'),
      eqn('governing', 'van\'t Hoff / Gibbs–Helmholtz', 'd ln K / dT = ΔH° / (R T²)'),
      eqn('constitutive', 'Kirchhoff correction', 'ΔH°(T) = ΔH°₂₉₈ + ΔCp (T − 298.15);   ΔS°(T) = ΔS°₂₉₈ + ΔCp ln(T/298.15)'),
      eqn('governing', 'Mass action', 'K = (c_C/c°) / [ (c_A/c°) (c_B/c°)ᵇ ],  c_A = c_A0(1 − X),  c_B = c_A0(θ − bX),  c_C = c_A0 X'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', s.dH < 0 ? `The reaction is exothermic: every 10 K rise multiplies K by ${F(Math.exp(10 * s.dH / (R * T * T)), 3)}. Run as cool as the kinetics allow, or use a falling temperature profile (hot start for rate, cool finish for conversion).` : `The reaction is endothermic: every 10 K rise multiplies K by ${F(Math.exp(10 * s.dH / (R * T * T)), 3)}. Raising temperature improves both rate and equilibrium conversion.`, 'van\'t Hoff slope at the operating temperature.'),
      rec('medium', X < 0.9 ? `Equilibrium limits conversion to ${F(100 * X, 3)} %. Raise the B/A ratio, remove product continuously (precipitation, stripping) or change temperature to shift the limit.` : `Equilibrium is favourable (${F(100 * X, 4)} % conversion); the process is likely rate limited, so check the kinetics suite before changing conditions.`, 'Computed equilibrium conversion.'),
      rec('medium', `Plan for a heat duty of ${F(Math.abs(Q) / 1e6, 3)} MJ per m³ of batch (${Q > 0 ? 'cooling' : 'heating'}).`, 'Reaction enthalpy times converted amount.'),
      rec('low', `Residual B is ${F(cB, 3)} mol/m³. Choose the smallest excess that still meets the conversion target; unreacted reagent is raw-material cost and effluent load.`, 'Equilibrium composition.', true),
    ],
    publish: { 'rxn.conversion': X, 'thermo.log10K': log10K, 'thermo.dG': s.dG },
    numerics: { method: 'Closed-form ΔG°(T) with Kirchhoff correction; Brent root finding on the logarithmic mass-action equation', dof: 1, steps: Ts.length, iterations: Ts.length + 1, residual: resK, notes: 'Ideal solution, 1 mol/L standard state.' },
  };
}

// ---------------------------------------------------------------- regular solution
function runRegular(inp) {
  const T = num(inp, 'Tr'), Om = num(inp, 'omega') * 1e3, x0 = num(inp, 'x0'), a = Om / (R * T), Tc = Om / (2 * R), warnings = [], checks = [];
  const gm = x => R * T * (x * Math.log(x) + (1 - x) * Math.log(1 - x)) + Om * x * (1 - x);
  const gap = Om > 0 && T < Tc;
  let xa = 0.5, xb = 0.5, s1 = 0.5, s2 = 0.5;
  if (gap) {
    const g = x => Math.log(x / (1 - x)) - a * (2 * x - 1);
    xa = g(1e-15) >= 0 ? Math.exp(-a) : safe(() => N.brent(g, 1e-15, 0.5 - 1e-9, 1e-15), 0.5); xb = 1 - xa;
    const d = Math.sqrt(Math.max(0, 1 - 2 / a)); s1 = 0.5 * (1 - d); s2 = 0.5 * (1 + d);
  }
  const state = !gap ? 'single phase (stable)' : (x0 <= xa || x0 >= xb) ? 'single phase (stable)' : (x0 > s1 && x0 < s2) ? 'unstable (spinodal decomposition)' : 'metastable (nucleation and growth)';
  const fB = gap && x0 > xa && x0 < xb ? (x0 - xa) / (xb - xa) : (x0 >= xb && gap ? 1 : 0);
  const d2 = R * T / (x0 * (1 - x0)) - 2 * Om;
  const xs = N.linspace(0.002, 0.998, 250);
  const lines = [];
  if (Om > 0) {
    lines.push({ name: 'binodal', x: xs, y: xs.map(x => Math.abs(x - 0.5) < 1e-9 ? Tc : Om * (1 - 2 * x) / (R * Math.log((1 - x) / x))) }, { name: 'spinodal', x: xs, y: xs.map(x => 2 * Om * x * (1 - x) / R) });
  } else lines.push({ name: 'no miscibility gap (Ω ≤ 0)', x: [0, 1], y: [0, 0] });
  lines.push({ name: 'state point', x: [x0, x0], y: [T * 0.98, T * 1.02] });
  const gline = [{ name: 'ΔG_mix', x: xs, y: xs.map(x => gm(x) / 1000) }];
  if (gap) gline.push({ name: 'common tangent', x: [xa, xb], y: [gm(Math.max(xa, 1e-300)) / 1000 || 0, gm(Math.min(xb, 1 - 1e-16)) / 1000 || 0] });
  const resB = gap && xa > 1e-14 ? Math.abs(Math.log(xa / (1 - xa)) - a * (2 * xa - 1)) : 0;
  checks.push(chk('Equal chemical potentials on the binodal', resB < 1e-8, `residual ${F(resB, 2)} (dimensionless)`));
  checks.push(chk('Thermodynamic stability at the overall composition (∂²G/∂x² > 0)', d2 > 0, `∂²G/∂x² = ${F(d2, 4)} J/mol`));
  checks.push(chk('Spinodal lies inside the binodal', !gap || (s1 >= xa - 1e-12 && s2 <= xb + 1e-12), gap ? `binodal ${F(xa, 4)}–${F(xb, 4)}, spinodal ${F(s1, 4)}–${F(s2, 4)}` : 'no two-phase region at this temperature'));
  if (Om <= 0) warnings.push('Ω ≤ 0: the regular-solution model predicts complete miscibility (ordering tendencies are not modelled).');
  warnings.push('The regular solution is symmetric with a temperature-independent Ω and ideal configurational entropy; real systems need sub-regular (Redlich–Kister) or CALPHAD descriptions.');
  const gA = Math.exp(N.clamp(a * x0 * x0, -50, 50)), gB = Math.exp(N.clamp(a * (1 - x0) ** 2, -50, 50));
  return {
    summary: `With Ω = ${F(Om / 1000, 4)} kJ/mol the critical temperature is ${Om > 0 ? F(Tc, 4) + ' K' : 'absent (no miscibility gap)'}. At ${F(T, 4)} K and x_B = ${F(x0, 3)} the mixture is ${state}${gap ? `; coexisting compositions are x_B = ${F(xa, 3)} and ${F(xb, 3)}, and the spinodal spans ${F(s1, 3)}–${F(s2, 3)}` : ''}. Free energy of mixing at the state point is ${F(gm(x0) / 1000, 4)} kJ/mol.`,
    kpis: [
      kpi('Tc', 'Critical (consolute) temperature', Om > 0 ? Tc : 0, 'K', 'T_c = Ω/2R. Above it the two components mix in all proportions; below it a miscibility gap opens. Zero is reported when Ω ≤ 0.'),
      kpi('state', 'Stability at the state point', state, '', 'Classification from the position of the overall composition relative to binodal and spinodal: stable, metastable (needs nucleation) or unstable (spontaneous spinodal decomposition).'),
      kpi('xa', 'Binodal composition, B-poor phase', xa, '', 'Equilibrium mole fraction of B in the B-poor phase from the common-tangent condition; the solubility limit of B in A.'),
      kpi('xb', 'Binodal composition, B-rich phase', xb, '', 'Equilibrium mole fraction of B in the B-rich phase; by symmetry of the model it equals 1 − x_a.'),
      kpi('s1', 'Spinodal composition (lower)', s1, '', 'Composition where ∂²G/∂x² = 0. Between the two spinodal points the mixture decomposes without a nucleation barrier.'),
      kpi('fB', 'Phase fraction of the B-rich phase', fB, '', 'Lever rule (x₀ − x_a)/(x_b − x_a) at the overall composition; zero or one outside the two-phase region.'),
      kpi('gmix', 'Free energy of mixing at x₀', gm(x0) / 1000, 'kJ/mol', 'ΔG_mix of the homogeneous mixture; inside the gap the two-phase state has a lower free energy given by the common tangent.'),
      kpi('gammaB', 'Activity coefficient of B at x₀', gB, '', `γ_B = exp(Ω(1−x)²/RT); values above one indicate positive deviation from ideality (γ_A = ${F(gA, 4)}).`),
    ],
    series: [
      { key: 'pd', title: 'Regular-solution phase diagram', xlabel: 'mole fraction of B', ylabel: 'T (K)', lines, explain: 'Binodal (coexistence) and spinodal (stability limit) curves. The region between them is metastable; inside the spinodal the mixture is unstable. The short vertical marker is the state point.' },
      { key: 'g', title: 'Free energy of mixing at the selected temperature', xlabel: 'mole fraction of B', ylabel: 'ΔG_mix (kJ/mol)', lines: gline, explain: 'A double-well curve signals a miscibility gap; the common tangent touches the curve at the two coexisting compositions.' },
    ],
    equations: [
      eqn('governing', 'Free energy of mixing', 'ΔG_mix = R T [x ln x + (1 − x) ln(1 − x)] + Ω x (1 − x)'),
      eqn('governing', 'Binodal (equal chemical potentials)', 'ln[x / (1 − x)] = (Ω / R T)(2x − 1)'),
      eqn('governing', 'Spinodal (stability limit)', '∂²ΔG_mix/∂x² = R T / [x(1 − x)] − 2Ω = 0'),
      eqn('formula', 'Critical temperature and lever rule', 'T_c = Ω / 2R;   f_β = (x₀ − x_α) / (x_β − x_α)'),
      eqn('constitutive', 'Activity coefficients', 'ln γ_A = Ω x_B² / R T;   ln γ_B = Ω x_A² / R T'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', gap ? (state.startsWith('unstable') ? `The composition lies inside the spinodal: it will decompose spontaneously into a fine interconnected two-phase structure. Use this deliberately for nanostructuring, or homogenise above ${F(Tc, 4)} K and quench to suppress it.` : state.startsWith('metastable') ? 'The composition is metastable: phase separation needs nucleation, so quenching rate and seeds control the precipitate density.' : `The composition is inside the single-phase field at this temperature; cooling below ${F(Math.abs(x0 - 0.5) < 1e-9 ? Tc : Om * (1 - 2 * x0) / (R * Math.log((1 - x0) / x0)), 4)} K would enter the gap.`) : 'No miscibility gap at this temperature: the solid or liquid solution is stable at every composition.', 'Position relative to binodal and spinodal.'),
      rec('medium', `Homogenisation requires holding above ${Om > 0 ? F(Tc, 4) + ' K' : 'no minimum temperature'}; ageing between spinodal and binodal gives coarser, nucleated precipitates.`, 'Critical temperature Ω/2R.'),
      rec('medium', 'Fit Ω to at least two measured solvus points before using the diagram for alloy or mixed-oxide design; an asymmetric measured solvus indicates that a sub-regular model is needed.', 'Symmetric single-parameter model.'),
      rec('low', 'Choosing the heat-treatment temperature just above the solvus instead of well above it saves furnace energy without changing the equilibrium state.', 'Solvus temperature from the binodal.', true),
    ],
    publish: { 'thermo.Tc': Om > 0 ? Tc : 0, 'thermo.binodal_lo': xa, 'thermo.binodal_hi': xb, 'thermo.phase_fraction': fB },
    numerics: { method: 'Brent root finding on the common-tangent condition; closed-form spinodal', dof: 1, steps: xs.length, iterations: 1, residual: resB, notes: 'Symmetric regular solution.' },
  };
}

// ---------------------------------------------------------------- vapour pressure / autoclave
function runVapour(inp) {
  const key = sel(inp, 'solvent'), s = SOLV[key], tC = num(inp, 'T'), T = tC + 273.15, phi0 = num(inp, 'fill') / 100, Prate = num(inp, 'Prate') * 1e6, Pamb = num(inp, 'Pamb');
  const warnings = [], checks = [], tcC = s.Tc - 273.15, superc = T >= s.Tc;
  const Ps = psat(s, tC), a = antSet(s, Math.min(tC, tcC));
  const tb = N.clamp(a[1] / (a[0] - Math.log10(Pamb / 133.322)) - a[2], -100, tcC);
  const tb2 = (() => { const a2 = antSet(s, tb); return N.clamp(a2[1] / (a2[0] - Math.log10(Pamb / 133.322)) - a2[2], -100, tcC); })();
  const phiAt = (Tk, p0 = phi0) => p0 * rhoLiq(key, T25) / rhoLiq(key, Tk);
  const phi = phiAt(T), full = phi >= 1;
  let Tfull = null; if (phiAt(0.995 * s.Tc) >= 1) Tfull = phiAt(T25) >= 1 ? T25 : N.bisect(Tk => phiAt(Tk) - 1, T25, 0.995 * s.Tc, 1e-4);
  const airAt = (Tk, ph) => Pamb * (Tk / T25) * Math.min((1 - phi0) / Math.max(1 - ph, 1e-9), 50);
  const Ptot = Tk => { const ph = phiAt(Tk), tc = Tk - 273.15; if (ph < 1) return psat(s, tc) + airAt(Tk, ph); return Math.min(psat(s, Tfull - 273.15) + airAt(Tfull, 0.98) + Math.log(ph) / s.kap, 1e9); };
  const P = Ptot(T), Hvap = R * T * T * a[1] * LN10 / (a[2] + Math.min(tC, tcC)) ** 2;
  const phiCrit = rhoLiq(key, T) / rhoLiq(key, T25);
  const Ts = N.linspace(25, Math.min(tcC, 400), 120);
  const inRange = tC >= s.ant[0][3] - 1e-9 && tC <= s.ant[s.ant.length - 1][4] + 1e-9;
  checks.push(chk('Temperature inside the Antoine correlation range', inRange, `${F(s.ant[0][3], 4)}–${F(s.ant[s.ant.length - 1][4], 4)} °C for ${s.name}; operating ${F(tC, 4)} °C`));
  checks.push(chk('Below the critical temperature', !superc, `T_c = ${F(tcC, 4)} °C, P_c = ${F(s.Pc / 1e6, 4)} MPa`));
  checks.push(chk('Vessel not liquid-full (no hydraulic lock)', !full, `liquid fill at temperature ${F(100 * Math.min(phi, 9.99), 4)} %`));
  checks.push(chk('Pressure below vessel rating', P <= Prate, `${F(P / 1e6, 4)} MPa vs rating ${F(Prate / 1e6, 4)} MPa`));
  if (!inRange) warnings.push('Temperature is outside the fitted range of the Antoine constants: the vapour pressure is an extrapolation.');
  if (superc) warnings.push(`Above the critical temperature of ${s.name} there is no vapour–liquid equilibrium: the value shown is the Antoine pressure at the critical temperature and is only a lower bound. The real pressure depends on the fill density and needs an equation of state.`);
  if (full) warnings.push(`Hydraulic lock: the liquid expands to fill the vessel at about ${F(Tfull - 273.15, 4)} °C. Beyond that the pressure follows the liquid compressibility (≈ ${F(1 / s.kap / 1e6 * 0.001, 3)} MPa per 0.1 % of thermal expansion) and can exceed any practical rating.`);
  else if (phi > 0.9) warnings.push(`Liquid occupies ${F(100 * phi, 3)} % of the vessel at temperature; the margin to hydraulic lock is small.`);
  if (key !== 'water') warnings.push('Liquid expansion uses the Rackett equation with the Yamada–Gunn compressibility estimate (accuracy of a few per cent); thermal decomposition of the solvent, which adds gas pressure, is not modelled.');
  if (key === 'dmf' && tC > 150) warnings.push('DMF decomposes slowly above about 150 °C to dimethylamine and CO, which raises the pressure beyond the vapour-pressure value.');
  const fills = [0.5, 0.7, 0.8, 0.9];
  return {
    summary: `At ${F(tC, 4)} °C the saturation pressure of ${s.name} is ${F(Ps / 1e6, 4)} MPa and the estimated total pressure in a vessel sealed ${F(100 * phi0, 3)} % full is ${F(P / 1e6, 4)} MPa (${F(100 * P / Prate, 3)} % of the ${F(Prate / 1e6, 3)} MPa rating). The liquid expands to occupy ${F(100 * Math.min(phi, 9.99), 4)} % of the vessel${full ? ' — the vessel is liquid-full and in hydraulic lock' : ''}.`,
    kpis: [
      kpi('Psat', 'Saturation vapour pressure', Ps / 1e6, 'MPa', `Antoine equation for ${s.name}; this is the autogenous pressure of the pure solvent while both liquid and vapour are present.`),
      kpi('Ptot', 'Estimated total vessel pressure', P / 1e6, 'MPa', 'Vapour pressure plus the trapped air compressed into the shrinking headspace (ideal gas, no dissolution), or the liquid-compressibility pressure once the vessel is liquid-full.', P <= Prate),
      kpi('margin', 'Pressure safety factor', Prate / P, '', 'Vessel rating divided by the estimated pressure. Values below about 1.5 leave too little margin for temperature overshoot and gas generation.', Prate / P >= 1.5),
      kpi('phiT', 'Liquid fill fraction at temperature', Math.min(phi, 9.99), '', 'Initial fill multiplied by the density ratio ρ(25 °C)/ρ(T) of the saturated liquid. One or more means the liquid has no room left to expand.', phi < 0.9),
      kpi('phiCrit', 'Maximum safe initial fill', phiCrit, '', 'Initial fill fraction that would make the vessel exactly liquid-full at the operating temperature; stay well below it.'),
      kpi('Tfull', 'Temperature of hydraulic lock', Tfull === null ? 'not reached below T_c' : Tfull - 273.15, '°C', 'Temperature at which the expanding liquid fills the whole vessel for the chosen initial fill.'),
      kpi('Tb', 'Boiling point at ambient pressure', tb2, '°C', 'Temperature at which the Antoine vapour pressure equals the ambient pressure; above it the solvent can only stay liquid in a sealed vessel.'),
      kpi('Hvap', 'Enthalpy of vaporisation (Clausius–Clapeyron)', Hvap / 1000, 'kJ/mol', 'ΔH_vap = R T² d ln P/dT from the slope of the Antoine curve (ideal vapour, negligible liquid volume); it overestimates near the critical point.'),
      kpi('Tr', 'Reduced temperature T/T_c', T / s.Tc, '', 'Closeness to the critical point. Above about 0.9 liquid density and vapour pressure change steeply and the simple correlations lose accuracy.', T / s.Tc < 0.9),
    ],
    series: [
      { key: 'pT', title: `Pressure of ${s.name} vs temperature`, xlabel: 'T (°C)', ylabel: 'pressure (MPa)', logy: true, lines: [{ name: 'saturation pressure', x: Ts, y: Ts.map(t => psat(s, t) / 1e6) }, { name: `total, ${F(100 * phi0, 3)} % fill`, x: Ts, y: Ts.map(t => Ptot(t + 273.15) / 1e6) }, { name: 'vessel rating', x: [Ts[0], Ts[Ts.length - 1]], y: [Prate / 1e6, Prate / 1e6] }], explain: 'Antoine saturation curve up to the critical temperature, and the estimated total pressure including trapped air. A sharp upturn of the total-pressure curve marks hydraulic lock.' },
      { key: 'fill', title: 'Liquid fill fraction vs temperature', xlabel: 'T (°C)', ylabel: 'liquid volume / vessel volume', lines: [...fills.map(f => ({ name: `${100 * f} % initial fill`, x: Ts, y: Ts.map(t => Math.min(phiAt(t + 273.15, f), 1.5)) })), { name: 'liquid-full', x: [Ts[0], Ts[Ts.length - 1]], y: [1, 1] }], explain: 'Thermal expansion of the saturated liquid for several initial fills (clipped at 1.5). A curve reaching one means the vessel becomes liquid-full at that temperature.' },
    ],
    equations: [
      eqn('formula', 'Solvent', chem.prettyFormula(s.f), `${s.name}, M = ${F(chem.molarMass(s.f), 5)} g/mol`),
      eqn('constitutive', 'Antoine equation', 'log₁₀ P_sat[mmHg] = A − B / (C + T[°C])', `A = ${a[0]}, B = ${a[1]}, C = ${a[2]}`),
      eqn('governing', 'Clausius–Clapeyron', 'd ln P_sat / dT = ΔH_vap / (R T²)'),
      eqn('constitutive', 'Saturated-liquid density', key === 'water' ? 'ρ/ρ_c = 1 + Σ bᵢ τ^(nᵢ),  τ = 1 − T/T_c  (saturated-water correlation)' : 'ρ(T)/ρ(T₀) = Z_RA^[(1 − T₀/T_c)^(2/7) − (1 − T/T_c)^(2/7)],  Z_RA = 0.29056 − 0.08775 ω  (Rackett, Yamada–Gunn)'),
      eqn('governing', 'Fill fraction and total pressure', 'φ(T) = φ₀ ρ(T₀)/ρ(T);   P = P_sat(T) + P₀ (T/T₀)(1 − φ₀)/(1 − φ(T));   liquid-full: ΔP = ln φ / κ_T'),
    ],
    checks, warnings,
    recommendations: [
      rec(full || P > Prate ? 'high' : 'medium', full ? `Reduce the initial fill below ${F(100 * phiCrit * 0.9, 3)} % (90 % of the liquid-full limit) or lower the temperature: the vessel is in hydraulic lock.` : P > Prate ? `Estimated pressure exceeds the rating. Lower the temperature, the fill or change to a higher-boiling solvent.` : `Keep the initial fill at or below ${F(100 * Math.min(phiCrit * 0.9, 0.8), 3)} % for this temperature; the common 80 % rule is not safe for every solvent and temperature.`, 'Liquid expansion and total-pressure estimate.'),
      rec('high', 'Fit a rupture disc or relief valve rated below the vessel limit and never rely on the vapour-pressure value alone: reactions that release gas (NH₃, CO₂, NOₓ from nitrates, decomposition products) add to it.', 'Autogenous pressure is a lower bound when gases are generated.'),
      rec('medium', `Cool the vessel below ${F(tb2, 3)} °C (the ambient boiling point) before opening, to avoid flashing of superheated solvent.`, 'Boiling point from the Antoine equation.'),
      rec('low', `Operating at the lowest temperature that gives the required phase lowers both pressure class and heating energy; for water, going from 200 to 180 °C lowers the vapour pressure by a factor of ${F(psat(SOLV.water, 200) / psat(SOLV.water, 180), 3)}. Recover and redistil organic solvents instead of discarding the mother liquor.`, 'Steep exponential dependence of vapour pressure on temperature.', true),
    ],
    publish: { 'proc.P': P, 'thermo.P_sat': Ps, 'thermo.fill_at_T': Math.min(phi, 9.99), 'proc.T_C': tC },
    numerics: { method: 'Antoine correlation, saturated-liquid density correlation, bisection for the liquid-full temperature', dof: 1, steps: Ts.length, iterations: 1, residual: 0, notes: 'Vapour mass neglected in the fill balance; trapped air treated as an ideal gas.' },
  };
}

// ---------------------------------------------------------------- nanoscale thermodynamics
function runNano(inp) {
  const T = num(inp, 'T') + 273.15, g = num(inp, 'gamma'), Vm = num(inp, 'Mw') * 1e-3 / num(inp, 'rho'), d = num(inp, 'dnm') * 1e-9, r = d / 2, S = num(inp, 'S');
  const Tm = num(inp, 'Tm'), Hf = num(inp, 'Hf') * 1e3, th = num(inp, 'theta') * Math.PI / 180, A = Math.pow(10, num(inp, 'logA')), warnings = [], checks = [];
  const v = Vm / N.N_A, kT = N.K_B * T, lnS = Math.log(S), lcap = 2 * g * Vm / (R * T);
  const tmOf = rr => Tm * (1 - 2 * g * Vm / (Hf * rr)), TmR = Math.max(0, tmOf(r));
  const solOf = rr => Math.exp(Math.min(lcap / rr, 230)), enh = solOf(r);
  const rc = lcap / lnS, nstar = (4 / 3) * Math.PI * rc ** 3 / v, dG = 16 * Math.PI * g ** 3 * v * v / (3 * (kT * lnS) ** 2), bar = dG / kT;
  const ct = Math.cos(th), fth = (2 + ct) * (1 - ct) ** 2 / 4, barHet = bar * fth;
  const J = A * Math.exp(-Math.min(bar, 700)), Jhet = A * Math.exp(-Math.min(barHet, 700));
  const rm = Math.cbrt(3 * v / (4 * Math.PI));
  checks.push(chk('Critical nucleus larger than one formula unit', nstar >= 1, `n* = ${F(nstar, 3)} formula units, r* = ${F(rc * 1e9, 3)} nm`));
  checks.push(chk('Continuum capillarity applicable (diameter ≥ 2 nm)', d >= 2e-9, `diameter ${F(d * 1e9, 3)} nm`));
  checks.push(chk('Gibbs–Thomson melting point positive', tmOf(r) > 0, `T_m(d) = ${F(tmOf(r), 4)} K`));
  checks.push(chk('Heterogeneous barrier does not exceed homogeneous barrier', barHet <= bar * (1 + 1e-12), `f(θ) = ${F(fth, 3)}`));
  if (nstar < 1) warnings.push('The critical nucleus is smaller than one formula unit: classical nucleation theory is not meaningful at this supersaturation (spinodal-like or non-classical nucleation regime).');
  if (d < 2e-9) warnings.push('Below about 2 nm the bulk interfacial energy and the continuum capillarity approximation are unreliable (Tolman-length corrections).');
  if (tmOf(r) <= 0) warnings.push('Gibbs–Thomson gives a non-physical melting point for this size; the linear relation is only valid for small relative depressions.');
  warnings.push('Classical nucleation theory rates are uncertain by many orders of magnitude because of the cubic dependence on interfacial energy; use them for trends and thresholds, not absolute rates.');
  const ds = N.logspace(1, 200, 80), Ss = N.logspace(Math.max(1.05, S / 50), S * 50, 80);
  const rr = N.linspace(rm * 0.2, Math.max(3 * rc, 4 * rm), 120), dmu = kT * lnS;
  return {
    summary: `For a ${F(d * 1e9, 3)} nm particle (γ = ${F(g, 3)} J/m², molar volume ${F(Vm * 1e6, 4)} cm³/mol) the Gibbs–Thomson melting point is ${F(TmR, 4)} K versus ${F(Tm, 4)} K in bulk, and the solubility is ${F(enh, 3)} times the flat-surface value. At supersaturation ${F(S, 3)} the critical nucleus radius is ${F(rc * 1e9, 3)} nm with a homogeneous barrier of ${F(bar, 3)} kT (${F(barHet, 3)} kT on the substrate).`,
    kpis: [
      kpi('Tm_d', 'Size-dependent melting point', TmR, 'K', `Gibbs–Thomson relation T_m(r) = T_m[1 − 2γV_m/(ΔH_f r)]; a depression of ${F(Tm - TmR, 3)} K relative to bulk. It sets sintering and annealing limits for nanoparticles.`),
      kpi('sol_enh', 'Solubility enhancement c(r)/c∞', enh, '', 'Ostwald–Freundlich factor exp(2γV_m/(rRT)). Small particles are more soluble than large ones, which drives Ostwald ripening.'),
      kpi('lcap', 'Capillary length 2γV_m/RT', lcap * 1e9, 'nm', 'Size scale below which curvature changes the solubility by more than a factor e; particles much larger than this behave as bulk.'),
      kpi('r_crit', 'Critical nucleus radius', rc * 1e9, 'nm', 'r* = 2γV_m/(RT ln S). Clusters smaller than this redissolve, larger ones grow.'),
      kpi('n_crit', 'Formula units in the critical nucleus', nstar, '', 'Volume of the critical nucleus divided by the molecular volume. Values of order one or less mean the classical picture has broken down.', nstar >= 1),
      kpi('barrier', 'Homogeneous nucleation barrier', bar, 'kT', 'ΔG*/kT = 16πγ³v²/(3(kT ln S)²)/kT. Below about 20–30 kT nucleation is fast (burst); above about 60–80 kT it is practically suppressed.'),
      kpi('barrier_het', 'Heterogeneous nucleation barrier', barHet, 'kT', `Homogeneous barrier multiplied by f(θ) = (2 + cos θ)(1 − cos θ)²/4 = ${F(fth, 3)}; explains why seeded substrates nucleate first.`),
      kpi('J', 'Homogeneous nucleation rate', J, '1/m³/s', 'J = A exp(−ΔG*/kT) with the stated prefactor. An order-of-magnitude estimate; zero means below floating-point range.'),
      kpi('dmu', 'Driving force per formula unit', dmu / N.Q_E, 'eV', 'Δμ = kT ln S, the free-energy gain per formula unit transferred from solution to crystal.'),
    ],
    series: [
      { key: 'tm', title: 'Melting point vs particle diameter', xlabel: 'diameter (nm)', ylabel: 'T_m (K)', logx: true, lines: [{ name: 'Gibbs–Thomson', x: ds, y: ds.map(x => Math.max(0, tmOf(x * 0.5e-9))) }, { name: 'bulk', x: [ds[0], ds[ds.length - 1]], y: [Tm, Tm] }], explain: 'Melting-point depression grows as 1/d; the curve is clipped at zero where the linear relation fails.' },
      { key: 'sol', title: 'Solubility enhancement vs particle diameter', xlabel: 'diameter (nm)', ylabel: 'c(r)/c∞', logx: true, logy: true, lines: [{ name: 'Ostwald–Freundlich', x: ds, y: ds.map(x => solOf(x * 0.5e-9)) }], explain: 'Equilibrium solubility of a curved particle relative to a flat surface. The steep rise below the capillary length is the origin of ripening.' },
      { key: 'dg', title: 'Free energy of cluster formation', xlabel: 'cluster radius (nm)', ylabel: 'ΔG / kT', lines: [{ name: 'ΔG(r)', x: rr.map(x => x * 1e9), y: rr.map(x => (4 * Math.PI * x * x * g - (4 / 3) * Math.PI * x ** 3 / v * dmu) / kT) }], explain: 'Surface term (∝ r²) against volume term (∝ r³). The maximum is the nucleation barrier at the critical radius.' },
      { key: 'bs', title: 'Nucleation barrier vs supersaturation', xlabel: 'supersaturation ratio', ylabel: 'ΔG*/kT', logx: true, logy: true, lines: [{ name: 'homogeneous', x: Ss, y: Ss.map(x => Math.max(1e-6, 16 * Math.PI * g ** 3 * v * v / (3 * (kT * Math.log(x)) ** 2) / kT)) }, { name: 'heterogeneous', x: Ss, y: Ss.map(x => Math.max(1e-6, fth * 16 * Math.PI * g ** 3 * v * v / (3 * (kT * Math.log(x)) ** 2) / kT)) }], explain: 'The barrier falls as 1/ln²S. The gap between the two curves is the window in which growth on seeds proceeds without bulk nucleation.' },
    ],
    equations: [
      eqn('governing', 'Gibbs–Thomson melting-point depression', 'T_m(r) = T_m,∞ [1 − 2 γ_sl V_m / (ΔH_f r)]'),
      eqn('governing', 'Ostwald–Freundlich (Kelvin) solubility', 'c(r) = c_∞ exp[2 γ V_m / (r R T)]'),
      eqn('governing', 'Classical nucleation theory', 'ΔG(r) = 4π r² γ − (4π r³ / 3v) kT ln S;   r* = 2 γ v / (kT ln S);   ΔG* = 16π γ³ v² / [3 (kT ln S)²]'),
      eqn('constitutive', 'Heterogeneous nucleation and rate', 'ΔG*_het = f(θ) ΔG*,  f(θ) = (2 + cos θ)(1 − cos θ)² / 4;   J = A exp(−ΔG*/kT)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', bar < 25 ? `The homogeneous barrier is only ${F(bar, 3)} kT: bulk nucleation will be fast. For seeded or epitaxial growth lower the supersaturation until the homogeneous barrier exceeds about 60 kT while the heterogeneous barrier (${F(barHet, 3)} kT) stays low.` : bar > 80 ? `The homogeneous barrier of ${F(bar, 3)} kT suppresses bulk nucleation; ${barHet < 40 ? 'growth will occur on seeds or walls only, which is the desired regime for oriented nanostructures.' : 'even heterogeneous nucleation is slow — add seeds or raise the supersaturation.'}` : `The homogeneous barrier (${F(bar, 3)} kT) is in the sensitive range where small changes of supersaturation or interfacial energy switch nucleation on and off; control the feed rate tightly.`, 'Classical nucleation barrier.'),
      rec('medium', `Particles below about ${F(2 * lcap * 1e9, 3)} nm diameter are markedly more soluble than bulk; expect ripening during ageing and protect small sizes with ligands or rapid quenching.`, 'Capillary length from the Ostwald–Freundlich relation.'),
      rec('medium', `Keep post-synthesis annealing below roughly ${F(0.5 * TmR, 4)} K (half the size-dependent melting point) to limit sintering and coarsening of ${F(d * 1e9, 3)} nm particles.`, 'Gibbs–Thomson melting point and the Tammann rule of thumb.'),
      rec('low', 'A well-wetted seed layer (small contact angle) lets growth proceed at lower supersaturation and temperature, cutting precursor waste to bulk precipitate and heating energy.', 'Heterogeneous barrier reduction factor f(θ).', true),
    ],
    publish: { 'nucl.rate': J, 'nucl.r_crit': rc, 'nucl.barrier_kT': bar, 'nucl.rate_het_prefactored': Jhet, 'thermo.capillary_length': lcap },
    numerics: { method: 'Closed-form capillarity relations (Gibbs–Thomson, Ostwald–Freundlich, classical nucleation theory)', dof: 1, steps: ds.length, iterations: 1, residual: 0, notes: 'Spherical particles, size-independent interfacial energy.' },
  };
}

const RUN = { speciation: runSpeciation, reaction: runReaction, regular: runRegular, vapour: runVapour, nano: runNano };

export default {
  id: 6,
  tagline: 'Aqueous speciation and oxide solubility, reaction equilibrium, regular-solution phase diagrams, autoclave vapour pressure and nanoscale capillarity thermodynamics.',
  inputs,
  presets: [
    { name: 'ZnO hydrothermal reference: 0.01 M Zn²⁺ + 0.02 M OH⁻ at 298 K', note: 'Reference ZnO precipitation case; the pH follows from the charge balance and the supersaturation is published for the growth suite.', values: { mode: 'speciation', metal: 'Zn', cM: 10, cBase: 20, pHsrc: 'base', T: 25, I: 0.05 } },
    { name: 'ZnO growth bath at 90 °C, pH 6.8', note: 'Near-neutral HMTA-type bath used for seeded nanorod growth.', values: { mode: 'speciation', metal: 'Zn', cM: 25, pHsrc: 'given', pH: 6.8, T: 90, I: 0.08 } },
    { name: 'Ferrihydrite precipitation from Fe(III) at pH 3', note: 'Strongly hydrolysing trivalent cation: precipitation starts in acid solution.', values: { mode: 'speciation', metal: 'Fe', cM: 1, pHsrc: 'given', pH: 3, T: 25, I: 0.1 } },
    { name: 'Water-filled autoclave, 200 °C, 80 % fill', note: 'Typical hydrothermal synthesis condition; checks pressure and liquid expansion.', values: { mode: 'vapour', solvent: 'water', T: 200, fill: 80, Prate: 6 } },
    { name: 'Gold nanoparticle, 5 nm: melting and nucleation', note: 'Classic size-dependent melting case.', values: { mode: 'nano', gamma: 0.27, Mw: 196.97, rho: 19300, dnm: 5, S: 1000, Tm: 1337.3, Hf: 12.55, T: 100, theta: 90 } },
    { name: 'Au–Ni-like miscibility gap (T_c ≈ 1083 K)', note: 'Regular-solution equivalent of a metallic miscibility gap used for spinodal nanostructuring.', values: { mode: 'regular', omega: 18, Tr: 900, x0: 0.4 } },
  ],
  run(inp) { return RUN[sel(inp, 'mode')](inp); },
};
