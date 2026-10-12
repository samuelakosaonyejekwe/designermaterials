// Suite 26 — Industrial Scale-Up, Automation and Quality.
// Reactor design equations and heat/mixing scale-up, statistical process control and capability, two-level factorial
// design of experiments, defect-limited yield / throughput / OEE, and automation readiness (sampling, PI loop, gauge R&R, EWMA/CUSUM).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS;
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const pct = v => F(100 * v, 3) + ' %';
// Student-t 97.5 % quantile: Simpson quadrature of the density (substitution x = √ν·tan θ) inverted by a bracketed root search;
// the Cornish–Fisher expansion in 1/ν is kept for ν > 200, where it is accurate to better than 10⁻⁸
function t975(nu) {
  const v = Math.max(1, nu), z = 1.959964, z3 = z ** 3, z5 = z ** 5, z7 = z ** 7;
  if (v > 200) return z + (z3 + z) / (4 * v) + (5 * z5 + 16 * z3 + 3 * z) / (96 * v * v) + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * v ** 3);
  const c = Math.exp(N.gammaLn((v + 1) / 2) - N.gammaLn(v / 2)) / Math.sqrt(Math.PI);
  const area = t => { const b = Math.atan(t / Math.sqrt(v)), n = 200, h = b / n; let s = 1 + Math.pow(Math.cos(b), v - 1); for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * Math.pow(Math.cos(i * h), v - 1); return c * s * h / 3; };
  return N.brent(t => area(t) - 0.475, 1.9, 13, 1e-10);
}
// generic delimited-text reader: returns header cells and rows of raw string cells
function readRows(text, minCols, minRows, label) {
  const lines = String(text ?? '').split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#') && !l.startsWith('//'));
  const cells = lines.map(l => l.split(/[,;\t]/).map(s => s.trim())).filter(c => c.length >= minCols);
  let head = []; if (cells.length && cells[0].slice(1).some(c => c === '' || !Number.isFinite(Number(c)))) head = cells.shift();
  if (cells.length < minRows) throw new Error(`${label}: at least ${minRows} data rows with ${minCols} columns are needed (found ${cells.length}). Paste CSV/TSV text with one header line.`);
  return { head, rows: cells };
}
const csv = (head, rows) => head.join(',') + '\n' + rows.map(r => r.join(',')).join('\n');

// ---------- chemistry of the production routes ----------
const ROUTES = {
  zno: { label: 'ZnO nanorods: Zn(NO₃)₂ + NaOH (hydrothermal)', re: [{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], pr: [{ f: 'ZnO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }] },
  ag: { label: 'Ag nanoparticles: AgNO₃ + NaBH₄ (borohydride reduction)', re: [{ f: 'AgNO3', phase: 'aq' }, { f: 'NaBH4', phase: 'aq' }, { f: 'H2O', phase: 'l' }], pr: [{ f: 'Ag', phase: 's' }, { f: 'H2', phase: 'g' }, { f: 'B(OH)3', phase: 'aq' }, { f: 'NaNO3', phase: 'aq' }] },
  au: { label: 'Au nanoparticles: HAuCl₄ + sodium citrate (Turkevich)', re: [{ f: 'HAuCl4', phase: 'aq' }, { f: 'Na3C6H5O7', phase: 'aq' }], pr: [{ f: 'Au', phase: 's' }, { f: 'Na2C5H4O5', phase: 'aq' }, { f: 'CO2', phase: 'g' }, { f: 'NaCl', phase: 'aq' }, { f: 'HCl', phase: 'aq' }] },
  tio2: { label: 'TiO₂ nanoparticles: Ti(OC₃H₇)₄ hydrolysis (sol–gel)', re: [{ f: 'Ti(OC3H7)4', phase: 'l' }, { f: 'H2O', phase: 'l' }], pr: [{ f: 'TiO2', phase: 's' }, { f: 'C3H8O', phase: 'l' }] },
};

// ---------- default synthetic tables (forward model + seeded noise) ----------
function spcSynth(mu, sd, shift, groups, n, seed, dec = 1) { const r = N.rng(seed), rows = []; for (let g = 1; g <= groups; g++) { const m = mu + (g > groups - 6 ? shift : 0) + r.normal(0, 0.25 * sd); for (let i = 0; i < n; i++) rows.push([g, (m + r.normal(0, sd)).toFixed(dec)]); } return csv(['subgroup', 'value'], rows); }
const SPC_DEF = spcSynth(1500, 38, 45, 25, 5, 26);        // nanorod length (nm), 25 batches × 5 rods, late mean shift
const SPC_LAT = spcSynth(250, 3.2, 0, 30, 4, 61);         // printed lattice strut diameter (µm)
const DOE_DEF = (() => { const r = N.rng(2026), rows = [], lv = [[70, 110], [10, 50], [1, 6]]; for (let rep = 0; rep < 2; rep++) for (let i = 0; i < 8; i++) { const x = [i & 1, (i >> 1) & 1, (i >> 2) & 1].map(b => 2 * b - 1); rows.push([...x.map((v, j) => lv[j][(v + 1) / 2]), (1200 + 250 * x[0] + 120 * x[1] + 180 * x[2] + 90 * x[0] * x[2] - 20 * x[0] * x[1] + r.normal(0, 25)).toFixed(0)]); } for (let c = 0; c < 3; c++) rows.push([90, 30, 3.5, (1140 + r.normal(0, 25)).toFixed(0)]); return csv(['T_C', 'c_mM', 't_h', 'length_nm'], rows); })();
const MSA_DEF = (() => { const r = N.rng(314), rows = [], ob = [0, 4, -3]; for (let p = 1; p <= 10; p++) { const tv = 1500 + r.normal(0, 40); for (let o = 1; o <= 3; o++) for (let k = 0; k < 2; k++) rows.push([p, o, (tv + ob[o - 1] + r.normal(0, 5)).toFixed(1)]); } return csv(['part', 'operator', 'value'], rows); })();
const STEPS_DEF = csv(['step', 'cycle_min', 'machines', 'yield_pct', 'availability_pct'], [['substrate clean', 5, 1, 99.5, 95], ['resist coat and bake', 6, 1, 99, 95], ['nano-imprint or e-beam exposure', 30, 2, 97, 85], ['develop', 8, 1, 99, 95], ['metal or dielectric deposition', 45, 2, 98, 88], ['lift-off or etch', 20, 1, 97, 92], ['optical inspection', 10, 1, 99.5, 97]]);
const STEPS_ZNO = csv(['step', 'cycle_min', 'machines', 'yield_pct', 'availability_pct'], [['precursor make-up', 40, 1, 99.5, 95], ['hydrothermal reactor batch', 300, 2, 97, 90], ['centrifuge and wash', 90, 1, 96, 92], ['drying', 240, 2, 99, 95], ['milling and packing', 45, 1, 99, 95]]);
const STEPS_AM = csv(['step', 'cycle_min', 'machines', 'yield_pct', 'availability_pct'], [['build preparation', 20, 1, 99.5, 97], ['photopolymer printing', 240, 4, 92, 85], ['solvent wash', 25, 1, 98, 95], ['UV and thermal post-cure', 90, 2, 98, 92], ['support removal', 30, 1, 96, 95], ['CT or optical inspection', 35, 1, 99, 95]]);

const SH = m => ({ mode: Array.isArray(m) ? m : [m] });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'reactor', group: 'Model', options: [
    { value: 'reactor', label: 'Reactor, heat-transfer and mixing scale-up' }, { value: 'spc', label: 'Statistical process control and capability' }, { value: 'doe', label: 'Design of experiments (2ᵏ factorial)' },
    { value: 'yield', label: 'Manufacturing yield, throughput and OEE' }, { value: 'auto', label: 'Automation readiness: sampling, control loop, gauge R&R, anomaly detection' }],
    help: 'Selects the scale-up or quality analysis. Each uses its own inputs and tables; results are engineering estimates that must be confirmed at pilot scale.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 99999, step: 1, sens: false, group: 'Model', help: 'Seed for randomised run orders and Latin-hypercube plans; the same seed reproduces the same plan.' },
  // reactor
  { key: 'rtype', label: 'Reactor type', type: 'select', value: 'batch', group: 'Reactor', showIf: SH('reactor'), options: [{ value: 'batch', label: 'Batch stirred tank' }, { value: 'cstr', label: 'Continuous stirred tank (CSTR)' }, { value: 'pfr', label: 'Plug-flow (tubular) reactor' }], help: 'Ideal reactor model used for sizing. Batch and plug flow need the same reaction time for the same conversion; a CSTR needs more volume for positive reaction orders.' },
  { key: 'route', label: 'Synthesis route', type: 'select', value: 'zno', group: 'Reactor', showIf: SH('reactor'), options: Object.entries(ROUTES).map(([value, r]) => ({ value, label: r.label })), help: 'Chemistry used for the stoichiometric mass balance. The first reactant listed is the limiting precursor to which concentration and conversion refer.' },
  { key: 'order', label: 'Reaction order n', value: 1, min: 0, max: 3, step: 0.5, group: 'Reactor', showIf: SH('reactor'), help: 'Apparent order of the rate law −r = k·cⁿ in the limiting precursor. Higher orders slow the approach to full conversion and penalise the CSTR.' },
  { key: 'k0', label: 'Pre-exponential factor k₀', unit: 's⁻¹·(mol/m³)^(1−n)', value: 1.4e5, min: 1e-3, max: 1e15, log: true, group: 'Reactor', showIf: SH('reactor'), help: 'Arrhenius pre-exponential factor of the rate constant, from laboratory kinetics (Suite 7 or the calibration in Suite 24).' },
  { key: 'Ea', label: 'Activation energy Eₐ', unit: 'kJ/mol', value: 60, min: 0, max: 300, group: 'Reactor', showIf: SH('reactor'), help: 'Arrhenius activation energy; with k₀ it sets how strongly the reaction time falls when the temperature is raised.' },
  { key: 'T', label: 'Reaction temperature', unit: '°C', value: 90, min: 20, max: 250, bus: 'proc.T_C', group: 'Reactor', showIf: SH('reactor'), help: 'Operating temperature of the reactor. Above the normal boiling point of the solvent a pressure vessel is required (see Suite 27).' },
  { key: 'c0', label: 'Precursor feed concentration', unit: 'mol/m³', value: 50, min: 0.01, max: 5000, log: true, bus: 'proc.c_precursor', group: 'Reactor', showIf: SH('reactor'), help: 'Initial or inlet concentration of the limiting precursor (50 mol/m³ = 0.05 mol/L). It sets the volumetric productivity.' },
  { key: 'X', label: 'Target conversion', value: 0.9, min: 0.05, max: 0.999, step: 0.01, bus: 'rxn.conversion', group: 'Reactor', showIf: SH('reactor'), help: 'Fraction of the limiting precursor converted in the reactor. The last few per cent cost disproportionate reaction time.' },
  { key: 'demand', label: 'Annual product demand', unit: 't/yr', value: 10, min: 0.001, max: 1e5, log: true, group: 'Reactor', showIf: SH('reactor'), help: 'Mass of isolated product required per year; the reactor volume and number of batches follow from it.' },
  { key: 'hours', label: 'Operating hours per year', unit: 'h/yr', value: 6000, min: 500, max: 8760, group: 'Production', showIf: SH(['reactor', 'yield']), help: 'Scheduled production time per year. 8000 h corresponds to continuous operation with planned shutdowns, 4000 h to two shifts on working days.' },
  { key: 'isoYield', label: 'Isolation yield', unit: '%', value: 92, min: 10, max: 100, group: 'Reactor', showIf: SH('reactor'), help: 'Fraction of the formed product recovered after separation, washing and drying.' },
  { key: 'tTurn', label: 'Batch turnaround time', unit: 'h', value: 2, min: 0, max: 48, group: 'Reactor', showIf: SH('reactor'), help: 'Non-reactive time per batch: charging, discharge, cooling and cleaning. It often dominates the cycle of fast reactions.' },
  { key: 'fill', label: 'Working fill fraction', unit: '%', value: 75, min: 20, max: 95, group: 'Reactor', showIf: SH('reactor'), help: 'Liquid volume as a percentage of vessel volume; the remainder is head space for thermal expansion, foaming and gas release.' },
  { key: 'Vlab', label: 'Laboratory reactor liquid volume', unit: 'L', value: 1, min: 0.01, max: 1000, log: true, group: 'Scale-up', showIf: SH('reactor'), help: 'Liquid volume of the laboratory or pilot reactor in which the recipe was developed; all scale-up ratios refer to it.' },
  { key: 'aspect', label: 'Vessel height/diameter ratio', value: 1.2, min: 0.5, max: 4, step: 0.1, group: 'Scale-up', showIf: SH('reactor'), help: 'Geometric ratio kept constant during scale-up (geometric similarity). Taller vessels offer more jacket area per volume.' },
  { key: 'U', label: 'Overall jacket heat-transfer coefficient', unit: 'W/m²/K', value: 400, min: 20, max: 3000, group: 'Scale-up', showIf: SH('reactor'), help: 'Overall coefficient between jacket fluid and reactor contents. Typical values are 150–700 W/m²/K for jacketed stirred tanks with aqueous contents.' },
  { key: 'dTj', label: 'Jacket temperature approach', unit: 'K', value: 20, min: 5, max: 100, group: 'Scale-up', showIf: SH('reactor'), help: 'Difference between jacket and target reaction temperature available for heating (and for removing reaction heat).' },
  { key: 'T0', label: 'Feed temperature', unit: '°C', value: 20, min: 5, max: 80, group: 'Scale-up', showIf: SH('reactor'), help: 'Temperature of the charged liquid before heating.' },
  { key: 'dHr', label: 'Reaction enthalpy', unit: 'kJ/mol', value: -30, min: -500, max: 200, group: 'Scale-up', showIf: SH('reactor'), help: 'Enthalpy of reaction per mole of limiting precursor converted (negative = exothermic). It determines the cooling duty and the adiabatic temperature rise.' },
  { key: 'rho', label: 'Liquid density', unit: 'kg/m³', value: 1000, min: 500, max: 2000, group: 'Scale-up', showIf: SH('reactor'), help: 'Density of the reaction mixture at operating temperature.' },
  { key: 'cp', label: 'Liquid heat capacity', unit: 'J/kg/K', value: 4180, min: 1000, max: 5000, group: 'Scale-up', showIf: SH('reactor'), help: 'Specific heat capacity of the reaction mixture.' },
  { key: 'mu', label: 'Liquid viscosity', unit: 'Pa·s', value: 5e-4, min: 1e-4, max: 1, log: true, group: 'Scale-up', showIf: SH('reactor'), help: 'Dynamic viscosity at operating temperature (water: 1.0 mPa·s at 20 °C, 0.31 mPa·s at 90 °C). It sets the Reynolds number and the micromixing time.' },
  { key: 'Nlab', label: 'Laboratory stirrer speed', unit: 'rpm', value: 400, min: 20, max: 3000, group: 'Mixing', showIf: SH('reactor'), help: 'Impeller speed used in the laboratory recipe.' },
  { key: 'Np', label: 'Impeller power number', value: 5, min: 0.2, max: 8, step: 0.1, group: 'Mixing', showIf: SH('reactor'), help: 'Turbulent power number: about 5 for a Rushton turbine, 1.3 for a pitched-blade turbine and 0.3 for a hydrofoil.' },
  { key: 'dRatio', label: 'Impeller/tank diameter ratio', value: 0.4, min: 0.2, max: 0.7, step: 0.05, group: 'Mixing', showIf: SH('reactor'), help: 'Ratio of impeller to tank diameter, kept constant on scale-up.' },
  { key: 'mixRule', label: 'Mixing scale-up rule', type: 'select', value: 'pv', group: 'Mixing', showIf: SH('reactor'), options: [{ value: 'pv', label: 'Constant power per volume' }, { value: 'tip', label: 'Constant tip speed' }, { value: 'time', label: 'Constant mixing time' }, { value: 're', label: 'Constant Reynolds number' }], help: 'Similarity criterion used to choose the stirrer speed at scale. All four are tabulated; the selected one feeds the key figures.' },
  { key: 'tInd', label: 'Nucleation induction time', unit: 's', value: 5, min: 1e-4, max: 1e4, log: true, group: 'Mixing', showIf: SH('reactor'), help: 'Characteristic time for nucleation after the reagents meet (from Suite 7). If micromixing is slower than this, nuclei form in unmixed, locally supersaturated fluid and the size distribution broadens.' },
  { key: 'Qunit', label: 'Flow per microreactor unit', unit: 'mL/min', value: 10, min: 0.01, max: 1000, log: true, group: 'Numbering-up', showIf: SH('reactor'), help: 'Design throughput of one continuous-flow microreactor channel or module.' },
  { key: 'dch', label: 'Microreactor channel diameter', unit: 'mm', value: 1, min: 0.05, max: 10, log: true, group: 'Numbering-up', showIf: SH('reactor'), help: 'Hydraulic diameter of the flow channel; smaller channels improve heat transfer and mixing at the price of pressure drop and fouling risk.' },
  // SPC
  { key: 'spcData', label: 'Subgroup measurements (subgroup, value)', type: 'table', columns: ['subgroup', 'value'], value: SPC_DEF, group: 'Quality data', showIf: SH(['spc', 'auto']), help: 'One row per measurement in time order: subgroup (batch or sample) identifier and measured value, e.g. nanorod length or film thickness. Subgroups of one measurement give an individuals/moving-range chart.' },
  { key: 'usl', label: 'Upper specification limit', value: 1650, min: -1e9, max: 1e9, group: 'Quality data', showIf: SH(['spc', 'auto']), help: 'Largest acceptable value of the quality characteristic, in the unit of the measurements.' },
  { key: 'lsl', label: 'Lower specification limit', value: 1350, min: -1e9, max: 1e9, group: 'Quality data', showIf: SH(['spc', 'auto']), help: 'Smallest acceptable value of the quality characteristic.' },
  { key: 'target', label: 'Target value', value: 1500, min: -1e9, max: 1e9, group: 'Quality data', showIf: SH('spc'), help: 'Nominal value of the characteristic, used for the Taguchi capability index Cpm.' },
  // DoE
  { key: 'doeData', label: 'Factorial runs (factor columns…, response)', type: 'table', columns: ['A', 'B', 'C', 'response'], value: DOE_DEF, group: 'Experiments', showIf: SH('doe'), help: 'Two-level factorial runs: two to four factor columns (actual or coded levels) followed by the response in the last column. Replicates and centre points are allowed.' },
  { key: 'goal', label: 'Optimisation goal', type: 'select', value: 'max', group: 'Experiments', showIf: SH('doe'), options: [{ value: 'max', label: 'Maximise the response' }, { value: 'min', label: 'Minimise the response' }, { value: 'target', label: 'Hit a target value' }], help: 'Objective used to locate the predicted optimum inside the experimental region.' },
  { key: 'yTarget', label: 'Response target', value: 1500, min: -1e9, max: 1e9, group: 'Experiments', showIf: SH('doe'), help: 'Desired response when the goal is to hit a target (same unit as the response column).' },
  { key: 'factors', label: 'Factors for new plans (name:low:high; …)', type: 'text', value: 'T_C:70:110; c_mM:10:50; t_h:1:6', group: 'Experiments', showIf: SH('doe'), help: 'Two to four factors with their low and high levels for the experimental plans offered as downloads (full factorial, central composite, Latin hypercube).' },
  { key: 'nCentre', label: 'Centre points in generated plans', value: 3, min: 0, max: 10, step: 1, sens: false, group: 'Experiments', showIf: SH('doe'), help: 'Number of centre-point runs added to the generated plans to estimate pure error and curvature.' },
  { key: 'nLhs', label: 'Latin-hypercube runs', value: 12, min: 4, max: 200, step: 1, sens: false, group: 'Experiments', showIf: SH('doe'), help: 'Number of runs of the space-filling Latin-hypercube plan.' },
  // yield
  { key: 'D0', label: 'Defect density', unit: '1/cm²', value: 0.5, min: 1e-3, max: 100, log: true, group: 'Yield', showIf: SH('yield'), help: 'Density of yield-killing defects (particles, pattern collapse, voids, missing struts) per unit area of product.' },
  { key: 'area', label: 'Product area per unit', unit: 'cm²', value: 1, min: 1e-3, max: 1e4, log: true, group: 'Yield', showIf: SH('yield'), help: 'Patterned or functional area of one device, die or panel.' },
  { key: 'critFrac', label: 'Critical-area fraction', value: 0.4, min: 0.01, max: 1, step: 0.01, group: 'Yield', showIf: SH('yield'), help: 'Fraction of the area in which a defect actually causes failure; dense metasurface patterns approach 1, sparse layouts are much lower.' },
  { key: 'alphaC', label: 'Defect clustering parameter α', value: 2, min: 0.3, max: 20, step: 0.1, group: 'Yield', showIf: SH('yield'), help: 'Cluster parameter of the negative-binomial model: small α means strongly clustered defects (higher yield than Poisson at equal density), large α approaches the Poisson limit.' },
  { key: 'ymodel', label: 'Yield model', type: 'select', value: 'negbin', group: 'Yield', showIf: SH('yield'), options: [{ value: 'poisson', label: 'Poisson (random defects)' }, { value: 'murphy', label: 'Murphy (triangular density distribution)' }, { value: 'seeds', label: 'Seeds (exponential density distribution)' }, { value: 'negbin', label: 'Negative binomial (clustered defects)' }], help: 'Defect-limited yield model used for the key figures; all four are compared in the charts.' },
  { key: 'steps', label: 'Process steps (step, cycle_min, machines, yield_pct, availability_pct)', type: 'table', columns: ['step', 'cycle_min', 'machines', 'yield_pct', 'availability_pct'], value: STEPS_DEF, group: 'Line', showIf: SH('yield'), help: 'One row per process step: name, cycle time per unit per machine in minutes, number of parallel machines, first-pass yield and availability in per cent.' },
  { key: 'units', label: 'Annual demand (good units)', unit: '1/yr', value: 12000, min: 1, max: 1e9, log: true, group: 'Line', showIf: SH('yield'), help: 'Number of good units (devices, wafers, panels or batches) required per year.' },
  { key: 'perf', label: 'Performance efficiency', unit: '%', value: 92, min: 10, max: 100, group: 'Line', showIf: SH('yield'), help: 'Actual running speed relative to ideal cycle time (minor stops and reduced speed), the performance factor of OEE.' },
  { key: 'tBatch', label: 'Campaign: batch cycle time', unit: 'h', value: 6, min: 0.1, max: 500, group: 'Campaign', showIf: SH('yield'), help: 'Duration of one production batch in the campaign schedule.' },
  { key: 'tClean', label: 'Campaign: cleaning per batch', unit: 'h', value: 1, min: 0, max: 100, group: 'Campaign', showIf: SH('yield'), help: 'Cleaning or reset time after each batch.' },
  { key: 'tChange', label: 'Campaign: changeover time', unit: 'h', value: 12, min: 0, max: 500, group: 'Campaign', showIf: SH('yield'), help: 'Product changeover (deep clean, requalification) at the end of each campaign.' },
  { key: 'nPerCamp', label: 'Batches per campaign', value: 10, min: 1, max: 500, step: 1, group: 'Campaign', showIf: SH('yield'), help: 'Number of consecutive batches of one product before a changeover.' },
  { key: 'batchSize', label: 'Batch size', unit: 'kg or units', value: 25, min: 0.001, max: 1e6, log: true, group: 'Campaign', showIf: SH('yield'), help: 'Product output of one batch before yield losses.' },
  // automation
  { key: 'tauP', label: 'Process time constant', unit: 's', value: 600, min: 0.01, max: 1e5, log: true, group: 'Control', showIf: SH('auto'), help: 'Dominant first-order time constant of the controlled variable (e.g. reactor temperature responding to heater power).' },
  { key: 'deadT', label: 'Dead time', unit: 's', value: 60, min: 0, max: 1e4, group: 'Control', showIf: SH('auto'), help: 'Transport and measurement delay between actuator change and first sensor response.' },
  { key: 'Kproc', label: 'Process gain', unit: 'output/input', value: 1.5, min: 0.001, max: 1000, log: true, group: 'Control', showIf: SH('auto'), help: 'Steady-state change of the controlled variable per unit change of the manipulated variable (e.g. K per % heater power).' },
  { key: 'tauFast', label: 'Fastest relevant phenomenon', unit: 's', value: 20, min: 1e-4, max: 1e5, log: true, group: 'Control', showIf: SH('auto'), help: 'Time constant of the fastest process event that must be observed (nucleation burst, exotherm onset, pressure spike). It sets the sensor sampling requirement.' },
  { key: 'lam', label: 'Closed-loop tuning ratio τc/θ', value: 1, min: 0.3, max: 10, step: 0.1, group: 'Control', showIf: SH('auto'), help: 'Desired closed-loop time constant relative to the dead time in the SIMC tuning rule: 1 is tight control, larger values are more robust and slower.' },
  { key: 'Ts', label: 'Controller sampling interval', unit: 's', value: 5, min: 1e-4, max: 1e4, log: true, group: 'Control', showIf: SH('auto'), help: 'Sampling and update interval of the digital controller and data acquisition.' },
  { key: 'chan', label: 'Sensor channels', value: 12, min: 1, max: 2000, step: 1, group: 'Control', showIf: SH('auto'), help: 'Number of logged sensor channels of the digital twin (temperature, pressure, pH, flow, power…).' },
  { key: 'msa', label: 'Gauge study (part, operator, value)', type: 'table', columns: ['part', 'operator', 'value'], value: MSA_DEF, group: 'Measurement system', showIf: SH('auto'), help: 'Crossed gauge repeatability and reproducibility study: every operator measures every part at least twice.' },
  { key: 'lamE', label: 'EWMA smoothing λ', value: 0.2, min: 0.02, max: 1, step: 0.01, group: 'Anomaly detection', showIf: SH('auto'), help: 'Weight of the newest observation in the exponentially weighted moving average; small values detect small sustained shifts.' },
  { key: 'LE', label: 'EWMA limit width L', value: 3, min: 1.5, max: 5, step: 0.1, group: 'Anomaly detection', showIf: SH('auto'), help: 'Width of the EWMA control limits in standard deviations of the statistic.' },
  { key: 'kC', label: 'CUSUM reference value k', unit: 'σ', value: 0.5, min: 0.1, max: 2, step: 0.05, group: 'Anomaly detection', showIf: SH('auto'), help: 'Slack of the tabular CUSUM in units of σ; half the shift to be detected is the usual choice.' },
  { key: 'hC', label: 'CUSUM decision interval h', unit: 'σ', value: 5, min: 1, max: 10, step: 0.1, group: 'Anomaly detection', showIf: SH('auto'), help: 'Alarm threshold of the CUSUM in units of σ; 4–5 gives a long in-control run length.' },
  { key: 'nBase', label: 'Baseline points', value: 50, min: 5, max: 5000, step: 1, sens: false, group: 'Anomaly detection', showIf: SH('auto'), help: 'Number of initial observations (assumed in control) used to estimate the baseline mean and standard deviation of the signal.' },
];

// ---------- mode: reactor / heat / mixing scale-up ----------
function runReactor(inp) {
  const w = [], rt = ROUTES[inp.route] || ROUTES.zno, re = rt.re.map(s => s.f), pr = rt.pr.map(s => s.f), bal = chem.balance(re, pr), ion = chem.ionicEquations(rt.re, rt.pr), aud = chem.auditEquation(re, pr, bal.coefficients);
  const Mp = chem.molarMass(pr[0]) / 1e3, nuP = bal.products[0] / bal.reactants[0], n = inp.order, T = inp.T + 273.15, k = inp.k0 * Math.exp(-inp.Ea * 1e3 / (R * T)), c0 = inp.c0, X = N.clamp(inp.X, 0.001, 0.9995);
  const kc = Math.max(k * Math.pow(c0, n - 1), 1e-300);                       // 1/s, characteristic first-order rate
  const tB = Xv => Math.abs(n - 1) < 1e-9 ? -Math.log(1 - Xv) / kc : (Math.pow(1 - Xv, 1 - n) - 1) / ((n - 1) * kc), tC = Xv => Xv / (kc * Math.pow(1 - Xv, n));
  const tRxn = inp.rtype === 'cstr' ? tC(X) : tB(X), eta = inp.isoYield / 100, fill = inp.rtype === 'pfr' ? 1 : inp.fill / 100, hrs = inp.hours * 3600;
  const molYr = inp.demand * 1e3 / Mp, Vyr = molYr / (c0 * X * nuP * eta);    // m³ of reaction liquid per year
  const rho = inp.rho, cp = inp.cp, mu = inp.mu, a = inp.aspect, Tj = inp.T + inp.dTj, dT = inp.T - inp.T0;
  const geom = Vl => { const Vv = Vl / (inp.fill / 100), D = Math.cbrt(4 * Vv / (Math.PI * a)), Hl = Vl / (Math.PI * D * D / 4), A = Math.PI * D * Hl + Math.PI * D * D / 4; return { Vv, D, A, AV: A / Vl, tHeat: dT > 0 ? rho * cp * Vl / (inp.U * A) * Math.log((Tj - inp.T0) / (Tj - inp.T)) : 0 }; };
  let Vl, Q, nBatch = 0, tCyc = 0, g;
  if (inp.rtype === 'batch') { Vl = inp.Vlab / 1e3; for (let i = 0; i < 60; i++) { g = geom(Vl); tCyc = tRxn + g.tHeat + inp.tTurn * 3600; nBatch = hrs / tCyc; const Vn = Vyr / nBatch; if (Math.abs(Vn - Vl) < 1e-10 * Vn) { Vl = Vn; break; } Vl = Vn; } g = geom(Vl); tCyc = tRxn + g.tHeat + inp.tTurn * 3600; nBatch = hrs / tCyc; Q = Vyr / hrs; }
  else { Q = Vyr / hrs; Vl = tRxn * Q; g = geom(Vl); }
  const Vves = Vl / fill, lab = geom(inp.Vlab / 1e3), scale = Vl / (inp.Vlab / 1e3), kgh = inp.demand * 1e3 / inp.hours;
  if (scale < 1) w.push(`The required reactor (${F(Vl * 1e3, 3)} L liquid) is smaller than the laboratory reactor: no scale-up is needed for this demand; ratios below 1 describe a scale-down.`);
  // numerical check of the design equation: integrate dX/dt = k·c0^(n−1)·(1−X)ⁿ
  const tChk = tB(X), sol = N.rk4((t, y) => [kc * Math.pow(Math.max(1 - y[0], 0), n)], [0], 0, tChk, 400), Xnum = sol.y.at(-1)[0];
  // heat
  const dHr = inp.dHr * 1e3, qHeatPeak = dT > 0 ? inp.U * g.A * (Tj - inp.T0) : 0, qCont = rho * cp * Q * Math.max(dT, 0), qRxn = inp.rtype === 'batch' ? Math.abs(dHr) * k * Math.pow(c0, n) * Vl : Math.abs(dHr) * Q * c0 * X;
  const qRem = inp.U * g.A * inp.dTj, heatRatio = dHr < 0 ? qRxn / Math.max(qRem, 1e-300) : 0, dTad = -dHr * c0 * X / (rho * cp), labRatio = dHr < 0 ? Math.abs(dHr) * k * Math.pow(c0, n) * (inp.Vlab / 1e3) / (inp.U * lab.A * inp.dTj) : 0;
  const kgBatch = inp.rtype === 'batch' ? inp.demand * 1e3 / nBatch : kgh * tRxn / 3600, eHeat = rho * cp * Math.max(dT, 0) / (c0 * X * nuP * eta * Mp) / 3.6e6;   // kWh per kg product (sensible heat of the liquid)
  if (heatRatio > 1) w.push(`Peak reaction heat (${F(qRxn / 1e3, 3)} kW) exceeds the jacket removal capacity (${F(qRem / 1e3, 3)} kW) at production scale: temperature control is lost unless feed is dosed, the concentration reduced or internal coils added.`);
  if (inp.T > 100) w.push('Operation above 100 °C with an aqueous medium requires a rated pressure vessel; check saturation pressure, fill fraction and relief in Suite 27.');
  // mixing
  const nu = mu / rho, dL = inp.dRatio * lab.D, dS = inp.dRatio * g.D, N1 = inp.Nlab / 60, VL = inp.Vlab / 1e3;
  const mix = (Ns, d, V) => { const P = inp.Np * rho * Ns ** 3 * d ** 5, eps = Math.max(P / (rho * V), 1e-300); return { N: Ns, P, PV: P / V, tip: Math.PI * Ns * d, Re: rho * Ns * d * d / mu, th: 5.2 / (Math.cbrt(inp.Np) * Ns * inp.dRatio ** 2), tauE: 17.24 * Math.sqrt(nu / eps), eta: Math.pow(nu ** 3 / eps, 0.25), shear: 11 * Ns }; };
  const rules = { pv: ['Constant P/V', r => Math.pow(r, 2 / 3)], tip: ['Constant tip speed', r => r], time: ['Constant mixing time', () => 1], re: ['Constant Reynolds number', r => r * r] }, m1 = mix(N1, dL, VL);
  const mAll = Object.fromEntries(Object.entries(rules).map(([key, [, fn]]) => [key, mix(N1 * fn(dL / dS), dS, Vl)])), ms = mAll[inp.mixRule] || mAll.pv;
  const DaMicro = ms.tauE / inp.tInd, DaLab = m1.tauE / inp.tInd, DaMacro = ms.th / (1 / kc), risk = DaMicro > 1 ? 'high' : DaMicro > 0.1 ? 'moderate' : 'low';
  if (m1.Re < 1e4) w.push(`Laboratory impeller Reynolds number ${F(m1.Re, 3)} is below 10⁴: the flow is transitional, the power number is not constant and the turbulent mixing correlations are only indicative.`);
  if (ms.Re < 1e4) w.push(`Impeller Reynolds number at scale (${F(ms.Re, 3)}) is below 10⁴ with the selected rule: turbulent correlations do not hold; choose another scale-up rule.`);
  if (DaMicro > 1) w.push(`Micromixing time at scale (${F(ms.tauE, 3)} s) exceeds the nucleation induction time (${F(inp.tInd, 3)} s): nucleation occurs in unmixed fluid, so a broader particle-size distribution than in the laboratory must be expected.`);
  // numbering-up
  const Qu = inp.Qunit * 1e-6 / 60, dch = inp.dch * 1e-3, tauP = tB(X), nUnits = Math.max(1, Math.ceil(Q / Qu)), Vu = tauP * Qu, Lch = Vu / (Math.PI * dch * dch / 4), ReCh = 4 * rho * Qu / (Math.PI * dch * mu), dP = 128 * mu * Lch * Qu / (Math.PI * dch ** 4), avCh = 4 / dch;
  if (ReCh > 2100) w.push(`Microreactor channel Reynolds number ${F(ReCh, 3)} is above the laminar limit: the Hagen–Poiseuille pressure drop underestimates the real value.`);
  if (dP > 2e6) w.push(`Channel pressure drop ${F(dP / 1e5, 3)} bar is impractical: shorten the residence time (higher temperature), enlarge the channel or reduce flow per unit.`);
  // series
  const Xs = N.linspace(0.05, 0.99, 60), ts = N.linspace(0, tChk * 1.3, 80), Vs = N.logspace(Math.log10(Math.min(VL, Vl)), Math.log10(Math.max(VL * 10, Vl * 1.5)), 40);
  const xOfT = t => { if (Math.abs(n - 1) < 1e-9) return 1 - Math.exp(-kc * t); const b = 1 + (n - 1) * kc * t; return n > 1 ? 1 - Math.pow(b, 1 / (1 - n)) : (b <= 0 ? 1 : 1 - Math.pow(b, 1 / (1 - n))); };
  const kAt = dTk => inp.k0 * Math.exp(-inp.Ea * 1e3 / (R * (T + dTk))) * Math.pow(c0, n - 1), xAt = (t, kk) => { if (Math.abs(n - 1) < 1e-9) return 1 - Math.exp(-kk * t); const b = 1 + (n - 1) * kk * t; return b <= 0 ? 1 : 1 - Math.pow(b, 1 / (1 - n)); };
  const massRows = [...rt.re.map((s, i) => [chem.prettyFormula(s.f), 'reactant', bal.reactants[i], F(bal.reactants[i] / bal.reactants[0] * X * c0 * Vl * chem.molarMass(s.f) / 1e3 * (inp.rtype === 'batch' ? 1 : 3600 / Math.max(tRxn, 1e-300)), 4)]), ...rt.pr.map((s, i) => [chem.prettyFormula(s.f), 'product', bal.products[i], F(bal.products[i] / bal.reactants[0] * X * c0 * Vl * chem.molarMass(s.f) / 1e3 * (inp.rtype === 'batch' ? 1 : 3600 / Math.max(tRxn, 1e-300)), 4)])];
  const mIn = N.sum(rt.re.map((s, i) => bal.reactants[i] * chem.molarMass(s.f))), mOut = N.sum(rt.pr.map((s, i) => bal.products[i] * chem.molarMass(s.f)));
  const typeName = { batch: 'batch stirred tank', cstr: 'CSTR', pfr: 'plug-flow reactor' }[inp.rtype] || 'batch stirred tank';
  return {
    summary: `To supply ${F(inp.demand, 3)} t/yr of ${chem.prettyFormula(pr[0])} at ${inp.T} °C (k = ${F(kc, 3)} s⁻¹) a ${typeName} needs ${F(Vl, 3)} m³ of reacting liquid (vessel ${F(Vves, 3)} m³, ${F(scale, 3)}× the laboratory volume) with a ${inp.rtype === 'batch' ? `reaction time of ${F(tRxn / 3600, 3)} h and ${F(nBatch, 3)} batches per year` : `residence time of ${F(tRxn / 60, 3)} min at ${F(Q * 3.6e6, 3)} L/h`}. Jacket area per volume falls from ${F(lab.AV, 3)} to ${F(g.AV, 3)} m⁻¹ and the heat-up time rises from ${F(lab.tHeat / 60, 3)} to ${F(g.tHeat / 60, 3)} min. With the rule "${(rules[inp.mixRule] || rules.pv)[0]}" the micromixing time is ${F(ms.tauE, 3)} s against an induction time of ${F(inp.tInd, 3)} s (Da = ${F(DaMicro, 3)}): ${risk} risk of particle-size broadening.`,
    kpis: [
      kp('k', 'Rate constant at temperature', k, 's⁻¹·(mol/m³)^(1−n)', `Arrhenius rate constant k = k₀·exp(−Eₐ/RT) at ${inp.T} °C; the equivalent first-order rate k·c₀ⁿ⁻¹ is ${F(kc, 3)} s⁻¹.`, null),
      kp('t_rxn', inp.rtype === 'batch' ? 'Batch reaction time' : 'Residence time', tRxn, 's', `Time required to reach ${pct(X)} conversion from the ideal ${typeName} design equation for order ${n}; equal to ${F(tRxn / 3600, 3)} h.`, null),
      kp('V_liq', 'Reacting liquid volume', Vl, 'm³', `Liquid hold-up needed to meet the annual demand at the stated conversion, isolation yield and operating hours; ${F(scale, 3)} times the laboratory volume.`, null),
      kp('V_vessel', 'Vessel volume', Vves, 'm³', `Liquid volume divided by the working fill fraction; diameter ${F(g.D, 3)} m at height/diameter ${a}.`, null),
      kp('n_batch', inp.rtype === 'batch' ? 'Batches per year' : 'Volumetric flow rate', inp.rtype === 'batch' ? nBatch : Q * 3.6e6, inp.rtype === 'batch' ? '1/yr' : 'L/h', inp.rtype === 'batch' ? `Number of batches that fit into ${inp.hours} h/yr with a cycle of ${F(tCyc / 3600, 3)} h (reaction + heat-up + turnaround); each batch yields ${F(kgBatch, 3)} kg of isolated product.` : 'Continuous feed rate of reaction liquid required by the annual demand at the stated conversion and isolation yield.', null),
      kp('kg_h', 'Product throughput', kgh, 'kg/h', 'Average isolated product mass per operating hour implied by the annual demand and the operating schedule.', null),
      kp('av_ratio', 'Jacket area/volume relative to laboratory', g.AV / lab.AV, '–', `Surface-to-volume ratio at scale (${F(g.AV, 3)} m⁻¹) divided by the laboratory value (${F(lab.AV, 3)} m⁻¹); it scales as V^(−1/3) and limits heating, cooling and heat-loss behaviour.`, g.AV / lab.AV > 0.2),
      kp('t_heat', 'Heat-up time at scale', g.tHeat, 's', `Time to heat the charge from ${inp.T0} to ${inp.T} °C with a jacket at ${F(Tj, 3)} °C (lumped energy balance); the laboratory reactor needs ${F(lab.tHeat, 3)} s.`, g.tHeat < Math.max(tRxn, 600)),
      kp('q_heat', inp.rtype === 'batch' ? 'Peak heating duty' : 'Feed pre-heating duty', inp.rtype === 'batch' ? qHeatPeak : qCont, 'W', inp.rtype === 'batch' ? 'Jacket duty U·A·(T_j − T₀) at the start of heat-up; the utility system must deliver this power.' : 'Sensible heat ρ·cₚ·Q·(T − T₀) required to bring the continuous feed to reaction temperature (before heat recovery).', null),
      kp('heat_ratio', 'Reaction heat / jacket removal capacity', heatRatio, '–', `Peak reaction heat release (${F(qRxn, 3)} W) divided by U·A·ΔT of the jacket; the laboratory value is ${F(labRatio, 3)}. Above 1 the vessel cannot hold its temperature.`, heatRatio < 0.5),
      kp('dT_ad', 'Adiabatic temperature rise', dTad, 'K', 'Temperature change of the reacting liquid if all reaction heat stayed in the liquid (−ΔH·c₀·X/ρcₚ); a screening figure for loss-of-cooling scenarios.', Math.abs(dTad) < 50),
      kp('PV', 'Power per volume at scale', ms.PV, 'W/m³', `Stirrer power input per liquid volume with the selected rule at ${F(ms.N * 60, 3)} rpm (laboratory: ${F(m1.PV, 3)} W/m³ at ${inp.Nlab} rpm).`, ms.PV < 5000),
      kp('tip', 'Impeller tip speed at scale', ms.tip, 'm/s', `Peripheral impeller speed π·N·d, a measure of maximum shear (laboratory: ${F(m1.tip, 3)} m/s); high values can break agglomerates or fragile nanostructures.`, ms.tip < 8),
      kp('tau_E', 'Micromixing time at scale', ms.tauE, 's', `Engulfment micromixing time 17.2·√(ν/ε) from the mean energy dissipation (laboratory: ${F(m1.tauE, 3)} s); the Kolmogorov length is ${F(ms.eta * 1e6, 3)} µm.`, DaMicro < 1),
      kp('Da_mix', 'Mixing Damköhler number', DaMicro, '–', `Micromixing time divided by the nucleation induction time (laboratory: ${F(DaLab, 3)}). Below about 0.1 nucleation is kinetically controlled; above 1 it is mixing-controlled and the size distribution broadens.`, DaMicro < 0.1 ? true : DaMicro > 1 ? false : null),
      kp('theta', 'Bulk blend time at scale', ms.th, 's', `95 % homogenisation time from the turbulent correlation Nθ = 5.2·Po^(−1/3)·(T/D)² (laboratory: ${F(m1.th, 3)} s); it is ${F(DaMacro, 3)} times the reaction time scale 1/(k·c₀ⁿ⁻¹).`, DaMacro < 0.1),
      kp('n_units', 'Microreactor units (numbering-up)', nUnits, '', `Number of parallel flow channels of ${inp.Qunit} mL/min needed for the same production; each needs ${F(Lch, 3)} m of ${inp.dch} mm channel and ${F(dP / 1e5, 3)} bar pressure drop, with ${F(avCh / g.AV, 3)} times the area per volume of the tank.`, nUnits < 500),
      kp('e_heat', 'Sensible heating energy', eHeat, 'kWh/kg', 'Energy to heat the reaction liquid to temperature per kilogram of isolated product, without heat recovery or losses; dilute recipes are dominated by heating solvent.', eHeat < 20),
    ],
    series: [
      { key: 'lev', title: 'Required time versus conversion', xlabel: 'conversion X', ylabel: 'time (s)', logy: true, lines: [{ name: 'batch / plug flow', x: Xs, y: Xs.map(tB) }, { name: 'CSTR', x: Xs, y: Xs.map(tC) }, { name: 'design point', x: [X], y: [tRxn] }], explain: 'Reaction or residence time needed for a given conversion in the three ideal reactors. The gap between the CSTR and the plug-flow curve widens with conversion and reaction order.' },
      { key: 'xt', title: 'Conversion history', xlabel: 'time (s)', ylabel: 'conversion X', lines: [{ name: `${inp.T} °C`, x: ts, y: ts.map(xOfT) }, { name: `${inp.T + 10} °C`, x: ts, y: ts.map(t => xAt(t, kAt(10))) }, { name: `${inp.T - 10} °C`, x: ts, y: ts.map(t => xAt(t, kAt(-10))) }], explain: 'Batch (or plug-flow) conversion against time at the operating temperature and 10 K either side. The spread shows how tightly temperature must be controlled to hold the cycle time.' },
      { key: 'heat', title: 'Heat-transfer scale-up', xlabel: 'liquid volume (m³)', ylabel: 'A/V (1/m) and heat-up time (min)', logx: true, logy: true, lines: [{ name: 'jacket area per volume (1/m)', x: Vs, y: Vs.map(v => geom(v).AV) }, { name: 'heat-up time (min)', x: Vs, y: Vs.map(v => Math.max(geom(v).tHeat / 60, 1e-6)) }], explain: 'Jacket area per unit volume falls as V^(−1/3) and the heat-up time rises as V^(1/3) for geometrically similar vessels. The same trend governs how fast reaction heat can be removed.' },
      { key: 'mix', title: 'Micromixing time under the four scale-up rules', xlabel: 'liquid volume (m³)', ylabel: 'micromixing time (s)', logx: true, logy: true, lines: [...Object.entries(rules).map(([key, [nm, fn]]) => ({ name: nm, x: Vs, y: Vs.map(v => { const gg = geom(v); return mix(N1 * fn(dL / (inp.dRatio * gg.D)), inp.dRatio * gg.D, v).tauE; }) })), { name: 'induction time', x: [Vs[0], Vs.at(-1)], y: [inp.tInd, inp.tInd] }], explain: 'Engulfment micromixing time against scale for each similarity rule. Where a curve rises above the induction time, nucleation becomes mixing-limited and laboratory particle sizes are not reproduced.' },
    ],
    tables: [
      { key: 'rules', title: 'Mixing scale-up rules compared', columns: ['Rule', 'Speed (rpm)', 'Power (W)', 'P/V (W/m³)', 'Tip speed (m/s)', 'Reynolds', 'Blend time (s)', 'Micromixing time (s)', 'Mean shear rate (1/s)', 'Da (micromix/induction)'], rows: [['Laboratory', F(N1 * 60, 4), F(m1.P, 3), F(m1.PV, 3), F(m1.tip, 3), F(m1.Re, 3), F(m1.th, 3), F(m1.tauE, 3), F(m1.shear, 3), F(DaLab, 3)], ...Object.entries(rules).map(([key, [nm]]) => { const m = mAll[key]; return [nm + (key === inp.mixRule ? ' (selected)' : ''), F(m.N * 60, 4), F(m.P, 3), F(m.PV, 3), F(m.tip, 3), F(m.Re, 3), F(m.th, 3), F(m.tauE, 3), F(m.shear, 3), F(m.tauE / inp.tInd, 3)]; })], explain: 'Stirrer speed and resulting hydrodynamic figures at production scale for each similarity rule. No single rule preserves power input, shear, blend time and Reynolds number simultaneously; the controlling mechanism decides which to keep.' },
      { key: 'mass', title: inp.rtype === 'batch' ? 'Stoichiometric mass balance per batch (kg)' : 'Stoichiometric mass balance per hour (kg/h)', columns: ['Species', 'Role', 'Coefficient', inp.rtype === 'batch' ? 'Mass converted or formed (kg)' : 'Mass rate (kg/h)'], rows: massRows, explain: 'Masses consumed and formed at the target conversion from the balanced equation, before isolation losses and without excess reagent. Solvent is not included.' },
      { key: 'num', title: 'Scale-up versus numbering-up', columns: ['Quantity', 'Stirred vessel at scale', 'Microreactor array'], rows: [['Units', 1, nUnits], ['Volume per unit (L)', F(Vl * 1e3, 4), F(Vu * 1e3, 4)], ['Area per volume (1/m)', F(g.AV, 4), F(avCh, 4)], ['Characteristic mixing time (s)', F(ms.tauE, 3), F(dch * dch / (4 * 1e-9), 3) + ' (diffusion across channel, D = 10⁻⁹ m²/s)'], ['Pressure drop (bar)', '–', F(dP / 1e5, 3)], ['Channel Reynolds number', '–', F(ReCh, 3)]], explain: 'Comparison of one large stirred vessel with parallel continuous-flow channels of laboratory dimensions. Numbering-up keeps heat-transfer and residence-time behaviour identical to the laboratory but needs uniform flow distribution and fouling control.' },
    ],
    equations: [
      { kind: 'chemical', label: 'Balanced overall reaction', text: ion.molecular },
      { kind: 'ionic', label: 'Complete ionic equation', text: ion.completeIonic },
      { kind: 'ionic', label: 'Net ionic equation', text: ion.netIonic, note: ion.spectators.length ? 'Spectator ions: ' + [...new Set(ion.spectators)].join(', ') : undefined },
      { kind: 'constitutive', label: 'Rate law', text: '−r = k·cⁿ,   k = k₀·exp(−Eₐ/(R·T))' },
      { kind: 'governing', label: 'Ideal reactor design equations', text: 'batch/PFR: t = c₀·∫₀ˣ dX/(−r);   CSTR: τ = c₀·X/(−r(X));   V = τ·Q' },
      { kind: 'governing', label: 'Jacket heat-up (lumped)', text: 'ρ·cₚ·V·dT/dt = U·A·(T_j − T)  ⇒  t_heat = (ρcₚV/UA)·ln[(T_j − T₀)/(T_j − T)]' },
      { kind: 'formula', label: 'Stirrer power and similarity rules', text: 'P = Po·ρ·N³·d⁵;  P/V = const ⇒ N₂ = N₁(d₁/d₂)^(2/3);  tip speed ⇒ N₂ = N₁d₁/d₂;  mixing time ⇒ N₂ = N₁;  Re ⇒ N₂ = N₁(d₁/d₂)²' },
      { kind: 'formula', label: 'Mixing times and Damköhler number', text: 'N·θ₉₅ = 5.2·Po^(−1/3)·(T/d)²;   τ_E = 17.2·√(ν/ε),  ε = P/(ρV);   Da = τ_E/t_ind' },
      { kind: 'formula', label: 'Laminar channel pressure drop', text: 'ΔP = 128·μ·L·Q/(π·d⁴)' },
    ],
    checks: [
      { name: 'Reaction equation balanced (mass and charge)', pass: aud.massBalanced && aud.chargeBalanced, detail: `Σ reactant mass ${F(mIn, 6)} g/mol-eq vs Σ product mass ${F(mOut, 6)} g/mol-eq` },
      { name: 'Design equation reproduced by numerical integration', pass: Math.abs(Xnum - X) < 2e-3, detail: `RK4 integration to t = ${F(tChk, 4)} s gives X = ${F(Xnum, 5)} (target ${F(X, 5)})` },
      { name: 'Turbulent mixing regime at scale (Re > 10⁴)', pass: ms.Re >= 1e4, detail: `Re = ${F(ms.Re, 3)} (laboratory ${F(m1.Re, 3)})` },
      { name: 'Reaction heat removable through the jacket', pass: heatRatio < 1, detail: dHr < 0 ? `heat release/removal capacity = ${F(heatRatio, 3)}` : 'endothermic or athermal reaction: no cooling demand' },
      { name: 'Micromixing faster than nucleation (Da < 1)', pass: DaMicro < 1, detail: `τ_E = ${F(ms.tauE, 3)} s, t_ind = ${F(inp.tInd, 3)} s` },
      { name: 'Laminar flow in microreactor channel', pass: ReCh < 2100, detail: `Re = ${F(ReCh, 3)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: DaMicro > 1 ? 'high' : 'medium', text: DaMicro > 1 ? `Feed the base or reducing agent through a dip pipe into the impeller discharge, or pre-mix the streams in a static or jet mixer before the vessel: the vessel-average micromixing time (${F(ms.tauE, 3)} s) is longer than the induction time (${F(inp.tInd, 3)} s).` : `Keep the reagent addition point in the impeller zone; with Da = ${F(DaMicro, 3)} the particle size is ${DaMicro < 0.1 ? 'kinetically controlled and should transfer to scale' : 'partly mixing-sensitive, so confirm the size distribution at an intermediate scale'}.`, basis: 'Mixing Damköhler number (engulfment micromixing time versus induction time).' },
      { priority: heatRatio > 0.5 || g.tHeat > tRxn ? 'high' : 'medium', text: g.tHeat > tRxn ? `Heat-up (${F(g.tHeat / 60, 3)} min) is longer than the reaction time, so much of the conversion occurs during the ramp at scale: pre-heat the feeds through an external exchanger or redefine the recipe on thermal history instead of hold time.` : `Thermal history at scale differs from the laboratory (heat-up ${F(g.tHeat / 60, 3)} min against ${F(lab.tHeat / 60, 3)} min); programme the laboratory ramp to mimic the plant when establishing the growth recipe.`, basis: 'Surface-to-volume ratio falls as V^(−1/3) on geometric scale-up.' },
      { priority: 'medium', text: inp.rtype === 'cstr' && n > 0 ? `A plug-flow or batch reactor would need ${pct(tB(X) / tC(X))} of the CSTR residence time for the same conversion; use tubular flow or a cascade of tanks if the size distribution tolerates it.` : nUnits <= 50 ? `Numbering-up is practical here (${nUnits} channels): it preserves laboratory mixing and heat transfer and avoids re-development of the recipe.` : `Numbering-up would need ${nUnits} parallel channels; a stirred vessel with an external rapid-mixing nucleation stage is the more practical route at this demand.`, basis: 'Ideal-reactor design equations and channel count for equal throughput.' },
      { priority: 'medium', text: `Heating the dilute reaction liquid costs ${F(eHeat, 3)} kWh per kg of product. Raising the precursor concentration, where colloidal stability allows, and recovering heat from the discharged batch into the next feed cut energy use and solvent volume almost proportionally.`, basis: 'Sensible-heat demand per kilogram of product from the mass and energy balance.', sustainability: true },
      { priority: 'low', text: 'Validate the scale-up at an intermediate pilot volume (about ten times the laboratory scale) before committing equipment: the correlations used are vessel-averaged and do not resolve local dissipation near the feed point.', basis: 'Validity limits of the stirred-tank correlations.' },
    ],
    publish: { 'proc.volume': Vl, 'proc.flow': Q, 'proc.time': tRxn, 'proc.power': ms.P, 'proc.T_C': inp.T, 'rxn.conversion': X, 'rxn.yield_kg': kgBatch, 'rxn.heat': qRxn, 'fab.throughput_kg_h': kgh, 'fab.energy_kWh_kg': eHeat, 'material.formula': pr[0] },
    numerics: { method: 'Closed-form ideal-reactor design equations with fixed-point iteration for the batch cycle; RK4 verification of the conversion integral', dof: 1, steps: 400, iterations: 60, residual: Math.abs(Xnum - X), notes: 'Constant density, isothermal reaction, single nth-order rate law; mixing from turbulent stirred-tank correlations with vessel-mean energy dissipation.' },
  };
}

// ---------- shared: grouped quality data ----------
const D2 = [0, 0, 1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078], A2 = [0, 0, 1.880, 1.023, 0.729, 0.577, 0.483, 0.419, 0.373, 0.337, 0.308], D3 = [0, 0, 0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223], D4 = [0, 0, 3.267, 2.574, 2.282, 2.114, 2.004, 1.924, 1.864, 1.816, 1.777];
const c4 = n => Math.sqrt(2 / (n - 1)) * Math.exp(N.gammaLn(n / 2) - N.gammaLn((n - 1) / 2));
function readGroups(inp, w) {
  const t = readRows(inp.spcData, 1, 8, 'Subgroup measurement table'), two = t.rows.every(r => r.length >= 2 && Number.isFinite(Number(r[1])));
  const pairs = two ? t.rows.map(r => [r[0], Number(r[1])]) : t.rows.map((r, i) => [String(i + 1), Number(r[0])]).filter(p => Number.isFinite(p[1]));
  if (pairs.length < 8) throw new Error('Subgroup measurement table: at least 8 numeric measurements are needed.');
  const order = [], map = new Map(); for (const [g, v] of pairs) { if (!map.has(g)) { map.set(g, []); order.push(g); } map.get(g).push(v); }
  let groups = order.map(g => map.get(g)); const sizes = groups.map(g => g.length), cnt = {}; sizes.forEach(s => cnt[s] = (cnt[s] || 0) + 1);
  let n = +Object.keys(cnt).sort((a, b) => cnt[b] - cnt[a] || b - a)[0];
  if (groups.length < 4 && n > 1) { groups = pairs.map(p => [p[1]]); n = 1; w.push('Fewer than four subgroups: every measurement is treated as an individual value.'); }
  else if (sizes.some(s => s !== n)) { const before = groups.length; groups = groups.filter(g => g.length >= n).map(g => g.slice(0, n)); w.push(`Subgroup sizes are unequal: ${before - groups.length} subgroup(s) smaller than the common size ${n} were omitted and larger ones truncated.`); }
  if ((inp.mode === 'spc' && inp.spcData === SPC_DEF) || (inp.mode === 'auto' && inp.spcData === SPC_DEF)) w.push('Subgroup measurement table: synthetic demonstration data — replace with measurements.');
  return { groups, n, all: groups.flat(), unit: t.head[1] || '' };
}

// ---------- mode: SPC ----------
function runSpc(inp) {
  const w = [], { groups, n, all } = readGroups(inp, w), m = groups.length, Nn = all.length, xb = groups.map(N.mean), gm = N.mean(all);
  let sw, cl2, ucl2, lcl2, disp, dispName, type;
  if (n === 1) { const mr = all.slice(1).map((v, i) => Math.abs(v - all[i])), mrb = N.mean(mr); sw = mrb / 1.128; disp = mr; cl2 = mrb; ucl2 = 3.267 * mrb; lcl2 = 0; dispName = 'moving range'; type = 'Individuals–moving range (I–MR)'; }
  else if (n <= 10) { const rg = groups.map(g => N.maxOf(g) - N.minOf(g)), rb = N.mean(rg); sw = rb / D2[n]; disp = rg; cl2 = rb; ucl2 = D4[n] * rb; lcl2 = D3[n] * rb; dispName = 'range'; type = 'X̄–R'; }
  else { const sg = groups.map(N.std), sb = N.mean(sg), c = c4(n), k3 = 3 * Math.sqrt(1 - c * c) / c; sw = sb / c; disp = sg; cl2 = sb; ucl2 = (1 + k3) * sb; lcl2 = Math.max(0, 1 - k3) * sb; dispName = 'standard deviation'; type = 'X̄–S'; }
  const sx = sw / Math.sqrt(n), ucl = gm + 3 * sx, lcl = gm - 3 * sx, so = N.std(all), pts = n === 1 ? all : xb, z = pts.map(v => sx > 0 ? (v - gm) / sx : 0);
  // Western Electric rules
  const viol = []; const run = (len, need, thr, name) => { for (let i = len - 1; i < z.length; i++) { const seg = z.slice(i - len + 1, i + 1); if (seg.filter(v => v > thr).length >= need || seg.filter(v => v < -thr).length >= need) viol.push([i + 1, name]); } };
  z.forEach((v, i) => { if (Math.abs(v) > 3) viol.push([i + 1, 'Rule 1: one point beyond 3σ']); }); run(3, 2, 2, 'Rule 2: two of three beyond 2σ (same side)'); run(5, 4, 1, 'Rule 3: four of five beyond 1σ (same side)'); run(8, 8, 0, 'Rule 4: eight consecutive on one side');
  const dispOut = disp.map((v, i) => v > ucl2 || v < lcl2 ? i + 1 : 0).filter(Boolean), r1 = viol.filter(v => v[1].startsWith('Rule 1')).length, inControl = viol.length === 0 && dispOut.length === 0;
  let usl = inp.usl, lsl = inp.lsl; if (usl < lsl) { [usl, lsl] = [lsl, usl]; w.push('Upper specification limit was below the lower limit; the two were exchanged.'); }
  const tol = usl - lsl, swS = Math.max(sw, 1e-12 * Math.max(Math.abs(gm), 1)), soS = Math.max(so, 1e-12 * Math.max(Math.abs(gm), 1));
  if (tol <= 0) w.push('Specification limits coincide: capability indices are zero by definition.');
  const cap = s => ({ cp: tol / (6 * s), cpk: Math.min(usl - gm, gm - lsl) / (3 * s), ppm: 1e6 * (N.normCdf(lsl, gm, s) + 1 - N.normCdf(usl, gm, s)) }), cw = cap(swS), co = cap(soS);
  const cpm = tol / (6 * Math.sqrt(soS * soS + (gm - inp.target) ** 2)), zst = 3 * cw.cpk, cpkLo = cw.cpk - 1.645 * Math.sqrt(1 / (9 * Nn) + cw.cpk ** 2 / (2 * (Nn - 1)));
  const m2 = N.mean(all.map(v => (v - gm) ** 2)) || 1e-300, sk = N.mean(all.map(v => (v - gm) ** 3)) / m2 ** 1.5, ku = N.mean(all.map(v => (v - gm) ** 4)) / (m2 * m2) - 3, jb = Nn / 6 * (sk * sk + ku * ku / 4), normal = jb < 5.99;
  const justified = inControl && normal && Nn >= 100;
  if (!inControl) w.push('The process is not in statistical control: capability indices describe a mixture of conditions and are not a prediction of future performance.');
  if (!normal) w.push(`Data depart from normality (Jarque–Bera ${F(jb, 3)} > 5.99): the ppm estimates from the normal model are unreliable; transform the data or use empirical percentiles.`);
  if (Nn < 100) w.push(`Only ${Nn} measurements: capability indices have wide confidence bounds (lower 95 % bound of Cpk ${F(cpkLo, 3)}); at least 100 measurements in 25 subgroups are customary.`);
  const idx = pts.map((_, i) => i + 1), h = N.histogram(all, 20), hy = N.maxOf(h.y), offc = cw.cp > 0 ? 1 - cw.cpk / cw.cp : 0;
  return {
    summary: `${type} chart of ${Nn} measurements in ${m} subgroups of ${n}: centre line ${F(gm, 5)}, control limits ${F(lcl, 5)} to ${F(ucl, 5)}, within-subgroup σ ${F(sw, 3)}. ${viol.length + dispOut.length ? `${viol.length} Western-Electric signal(s) on the location chart and ${dispOut.length} point(s) outside the ${dispName} limits show special-cause variation.` : 'No rule violations: the process shows only common-cause variation.'} Against the specification ${F(lsl, 5)}–${F(usl, 5)}: Cp = ${F(cw.cp, 3)}, Cpk = ${F(cw.cpk, 3)}, Ppk = ${F(co.cpk, 3)}, about ${F(co.ppm, 3)} ppm expected out of specification (long-term). ${justified ? 'The capability figures are statistically supported.' : 'The capability figures are indicative only (see warnings).'}`,
    kpis: [
      kp('mean', 'Process mean (centre line)', gm, '', `Grand average of all measurements; it lies ${F(Math.abs(gm - (usl + lsl) / 2), 3)} from the centre of the specification, which costs ${pct(offc)} of the potential capability.`, null),
      kp('sigma_w', 'Within-subgroup standard deviation', sw, '', `Short-term process standard deviation estimated from the average ${dispName} (${n === 1 ? 'MR̄/1.128' : n <= 10 ? 'R̄/d₂' : 's̄/c₄'}); it excludes drift between subgroups.`, null),
      kp('sigma_o', 'Overall standard deviation', so, '', 'Sample standard deviation of all measurements, including batch-to-batch shifts and drift; the ratio to the within-subgroup value shows how much variation comes from between-batch effects.', so < 1.3 * sw),
      kp('ucl', 'Upper control limit', ucl, '', 'Centre line plus three standard errors of the plotted statistic. Control limits describe the voice of the process and are unrelated to the specification limits.', null),
      kp('lcl', 'Lower control limit', lcl, '', 'Centre line minus three standard errors of the plotted statistic; points below it signal a special cause.', null),
      kp('viol', 'Western-Electric rule signals', viol.length, '', `Number of rule violations on the location chart (${r1} beyond 3σ, the rest run and zone rules). Any signal calls for investigation of an assignable cause before capability is quoted.`, viol.length === 0),
      kp('cp', 'Cp (potential capability)', cw.cp, '–', 'Specification width divided by six within-subgroup standard deviations: the capability the process would have if perfectly centred.', cw.cp >= 1.33),
      kp('cpk', 'Cpk (actual short-term capability)', cw.cpk, '–', `Distance from the mean to the nearer specification limit in units of three within-subgroup standard deviations; the lower 95 % confidence bound is ${F(cpkLo, 3)}.`, cw.cpk >= 1.33),
      kp('pp', 'Pp (overall performance)', co.cp, '–', 'Specification width divided by six overall standard deviations; lower than Cp when the process mean wanders between subgroups.', co.cp >= 1.33),
      kp('ppk', 'Ppk (overall performance, centring included)', co.cpk, '–', 'Long-term performance index using the overall standard deviation; the figure customers experience over many batches.', co.cpk >= 1.33),
      kp('cpm', 'Cpm (Taguchi index)', cpm, '–', `Capability relative to the target ${F(inp.target, 5)}: penalises deviation of the mean from target as well as spread.`, cpm >= 1),
      kp('ppm', 'Expected out-of-specification (long term)', co.ppm, 'ppm', `Parts per million outside the specification from the normal model with the overall standard deviation (short-term estimate: ${F(cw.ppm, 3)} ppm).`, co.ppm < 6210),
      kp('zst', 'Process sigma level', zst, 'σ', 'Distance from the mean to the nearer specification limit in short-term standard deviations (3·Cpk); the conventional "six sigma" scale assumes this short-term figure.', zst >= 4),
    ],
    series: [
      { key: 'xbar', title: n === 1 ? 'Individuals chart' : 'X̄ chart', xlabel: n === 1 ? 'observation' : 'subgroup', ylabel: 'value', lines: [{ name: n === 1 ? 'individual value' : 'subgroup mean', x: idx, y: pts }, { name: 'centre line', x: [1, idx.length], y: [gm, gm] }, { name: 'UCL', x: [1, idx.length], y: [ucl, ucl] }, { name: 'LCL', x: [1, idx.length], y: [lcl, lcl] }], explain: 'Location chart with three-sigma control limits. Points outside the limits, or non-random runs on one side of the centre line, indicate that an assignable cause has shifted the process.' },
      { key: 'disp', title: `${dispName[0].toUpperCase() + dispName.slice(1)} chart`, xlabel: n === 1 ? 'observation' : 'subgroup', ylabel: dispName, lines: [{ name: dispName, x: disp.map((_, i) => i + (n === 1 ? 2 : 1)), y: disp }, { name: 'centre line', x: [1, idx.length], y: [cl2, cl2] }, { name: 'UCL', x: [1, idx.length], y: [ucl2, ucl2] }, { name: 'LCL', x: [1, idx.length], y: [lcl2, lcl2] }], explain: 'Dispersion chart of the short-term variation. It must be in control before the location chart is interpreted, because the location limits are computed from this spread.' },
      { key: 'hist', title: 'Distribution against specification limits', xlabel: 'value', ylabel: 'count', lines: [{ name: 'measurements', x: h.x, y: h.y }, { name: 'LSL', x: [lsl, lsl], y: [0, hy] }, { name: 'USL', x: [usl, usl], y: [0, hy] }, { name: 'normal model (overall σ)', x: h.x, y: h.x.map(x => Nn * h.width * Math.exp(-0.5 * ((x - gm) / soS) ** 2) / (soS * Math.sqrt(2 * Math.PI))) }], explain: 'Histogram of all measurements with the specification limits and the fitted normal curve. The tail area beyond the limits is the predicted out-of-specification fraction.' },
    ],
    tables: [
      { key: 'viol', title: 'Rule violations', columns: [n === 1 ? 'Observation' : 'Subgroup', 'Signal'], rows: viol.length || dispOut.length ? [...viol.slice(0, 40), ...dispOut.map(i => [i, `${dispName} outside control limits`])] : [['–', 'none']], explain: 'Each Western-Electric signal with the subgroup at which it completes. Clusters of signals from one point onward indicate a sustained shift rather than isolated outliers.' },
      { key: 'cap', title: 'Capability summary', columns: ['Index', 'Within (short term)', 'Overall (long term)'], rows: [['Standard deviation', F(sw, 4), F(so, 4)], ['Cp / Pp', F(cw.cp, 4), F(co.cp, 4)], ['Cpk / Ppk', F(cw.cpk, 4), F(co.cpk, 4)], ['Expected ppm out of specification', F(cw.ppm, 4), F(co.ppm, 4)], ['Skewness / excess kurtosis', F(sk, 3), F(ku, 3)], ['Jarque–Bera statistic (5.99 = 5 % limit)', F(jb, 4), normal ? 'normality not rejected' : 'normality rejected']], explain: 'Short-term and long-term capability side by side. A large gap between the two columns means that removing between-batch variation, not reducing within-batch scatter, is the route to better capability.' },
    ],
    equations: [
      { kind: 'formula', label: 'Control limits', text: n === 1 ? 'x̄ ± 3·MR̄/1.128;   UCL_MR = 3.267·MR̄' : n <= 10 ? `x̿ ± A₂·R̄;   UCL_R = D₄·R̄,  LCL_R = D₃·R̄   (n = ${n}: A₂ = ${A2[n]}, D₃ = ${D3[n]}, D₄ = ${D4[n]}, d₂ = ${D2[n]})` : 'x̿ ± 3·s̄/(c₄·√n);   s̄·(1 ± 3·√(1 − c₄²)/c₄)' },
      { kind: 'formula', label: 'Capability indices', text: 'Cp = (USL − LSL)/(6σ_w),  Cpk = min(USL − μ, μ − LSL)/(3σ_w),  Pp, Ppk: same with overall σ,  Cpm = (USL − LSL)/(6·√(σ² + (μ − T)²))' },
      { kind: 'formula', label: 'Expected nonconforming fraction', text: 'p = Φ((LSL − μ)/σ) + 1 − Φ((USL − μ)/σ)' },
    ],
    checks: [
      { name: 'Location chart in control', pass: viol.length === 0, detail: `${viol.length} Western-Electric signal(s)` },
      { name: 'Dispersion chart in control', pass: dispOut.length === 0, detail: `${dispOut.length} point(s) outside limits` },
      { name: 'Normality not rejected (Jarque–Bera)', pass: normal, detail: `JB = ${F(jb, 4)} (limit 5.99)` },
      { name: 'Sample size sufficient for capability', pass: Nn >= 100 && m >= 25, detail: `${Nn} measurements in ${m} subgroups` },
      { name: 'Capability statistically justified', pass: justified, detail: justified ? 'in control, approximately normal, adequate sample' : 'one or more prerequisites not met' },
      { name: 'Short-term capability Cpk ≥ 1.33', pass: cw.cpk >= 1.33, detail: `Cpk = ${F(cw.cpk, 4)}, lower 95 % bound ${F(cpkLo, 4)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: viol.length ? 'high' : 'low', text: viol.length ? `Investigate the assignable cause starting at ${n === 1 ? 'observation' : 'subgroup'} ${viol[0][0]} (${viol[0][1]}): compare precursor lot, temperature trace, stirring and operator for that period before adjusting the recipe.` : 'Keep the chart running with the present limits and react only to rule violations; adjusting an in-control process increases its variation.', basis: 'Western-Electric run rules on the location chart.' },
      { priority: cw.cpk < 1.33 ? 'high' : 'low', text: cw.cpk < 1.33 ? (cw.cp >= 1.33 ? `Re-centre the process: shifting the mean by ${F((usl + lsl) / 2 - gm, 3)} raises Cpk from ${F(cw.cpk, 3)} to the potential value ${F(cw.cp, 3)} without reducing variation.` : `Reduce variation: σ must fall from ${F(sw, 3)} to ${F(tol / 8, 3)} for Cp = 1.33; use the factorial-design mode to find the factors that control the spread.`) : `Capability is adequate (Cpk ${F(cw.cpk, 3)}); sampling frequency can be reduced once 25 further in-control subgroups confirm it.`, basis: 'Cp and Cpk against the customary 1.33 requirement.' },
      { priority: so > 1.3 * sw ? 'medium' : 'low', text: so > 1.3 * sw ? `Overall σ (${F(so, 3)}) exceeds within-subgroup σ (${F(sw, 3)}): batch-to-batch effects dominate. Standardise precursor lots, ramp rates and ageing times between batches.` : 'Within- and between-subgroup variation are comparable: the remaining scatter is inherent to the process at the present settings.', basis: 'Ratio of overall to within-subgroup standard deviation.' },
      { priority: 'medium', text: `About ${F(co.ppm / 1e4, 3)} % of output is predicted outside specification; each rejected batch wastes its precursors, solvent and process energy. Detecting the shift at the first rule violation and holding the mean on target is the cheapest resource saving available.`, basis: 'Long-term expected nonconforming fraction.', sustainability: true },
    ],
    publish: { 'fab.yield': N.clamp(1 - co.ppm / 1e6, 0, 1), 'qual.cpk': cw.cpk, 'qual.ppk': co.cpk, 'qual.ppm': co.ppm, 'qual.sigma': so },
    numerics: { method: `Shewhart ${type} chart with tabulated control-chart constants; normal-model capability`, dof: Nn, steps: m, iterations: 1, residual: 0, notes: 'Capability indices assume a stable, normally distributed process; confidence bound of Cpk from the Bissell approximation.' },
  };
}

// ---------- mode: design of experiments ----------
function parseFactors(text) { const out = []; for (const p of String(text || '').split(';')) { const c = p.split(':').map(s => s.trim()); if (c.length >= 3 && c[0] && Number.isFinite(+c[1]) && Number.isFinite(+c[2]) && +c[1] !== +c[2]) out.push({ name: c[0].replace(/[,;\s]+/g, '_'), lo: Math.min(+c[1], +c[2]), hi: Math.max(+c[1], +c[2]) }); } return out.slice(0, 4); }
function runDoe(inp) {
  const w = [], t = readRows(inp.doeData, 3, 4, 'Factorial run table'), nc = Math.min(...t.rows.map(r => r.length)), rowsN = t.rows.map(r => [...r.slice(0, nc - 1), r[nc - 1]].map(Number)).filter(r => r.every(Number.isFinite));
  if (rowsN.length < 4) throw new Error('Factorial run table: at least 4 fully numeric runs are needed.');
  if (inp.doeData === DOE_DEF) w.push('Factorial run table: synthetic demonstration data — replace with measurements.');
  let k = nc - 1; if (k > 4) { w.push(`Only the first four of ${k} factor columns are analysed (2ᵏ analysis is limited to k ≤ 4).`); k = 4; }
  const names = Array.from({ length: k }, (_, j) => t.head[j] || 'ABCD'[j]), yName = t.head[nc - 1] || 'response', Y0 = rowsN.map(r => r[nc - 1]);
  const lo = names.map((_, j) => N.minOf(rowsN.map(r => r[j]))), hi = names.map((_, j) => N.maxOf(rowsN.map(r => r[j]))), act = names.map((_, j) => hi[j] > lo[j]);
  act.forEach((a, j) => { if (!a) w.push(`Factor ${names[j]} does not vary in the table and is excluded from the model.`); });
  const code = r => names.map((_, j) => act[j] ? (r[j] - (hi[j] + lo[j]) / 2) / ((hi[j] - lo[j]) / 2) : 0), Xc = rowsN.map(code), isF = Xc.map(x => x.every((v, j) => !act[j] || Math.abs(Math.abs(v) - 1) < 1e-6));
  const fx = Xc.filter((_, i) => isF[i]), fy = Y0.filter((_, i) => isF[i]), cy = Y0.filter((_, i) => !isF[i] && Xc[i].every(v => Math.abs(v) < 1e-6)), nf = fy.length, fac = names.map((_, j) => j).filter(j => act[j]);
  if (nf < 3 || !fac.length) throw new Error('Factorial run table: fewer than three runs lie at the low/high factor levels; a two-level design is required.');
  const distinct = new Set(fx.map(x => x.map(v => Math.sign(v)).join(''))).size;
  // model terms: all factor subsets up to the highest interaction order the design supports
  let subsets = []; const kk = fac.length; for (let mask = 1; mask < (1 << kk); mask++) subsets.push(fac.filter((_, b) => mask & (1 << b)));
  subsets.sort((a, b) => a.length - b.length || a[0] - b[0]); let maxOrd = kk; while (maxOrd > 1 && 1 + subsets.filter(s => s.length <= maxOrd).length > distinct) maxOrd--;
  if (maxOrd < kk) w.push(`The table contains ${distinct} distinct factorial points for ${kk} factors: interactions above order ${maxOrd} are aliased and were not estimated.`);
  subsets = subsets.filter(s => s.length <= maxOrd); if (1 + subsets.length > distinct) { subsets = subsets.slice(0, Math.max(1, distinct - 1)); w.push('The design is incomplete even for main effects; only the first estimable terms are fitted and the result is unreliable.'); }
  const term = (x, s) => s.reduce((p, j) => p * x[j], 1), row = x => [1, ...subsets.map(s => term(x, s))], ysc = Math.max(N.maxOf(fy.map(Math.abs)), 1e-300);
  let beta; try { beta = N.lstsq(fx.map(row), fy.map(v => v / ysc), 1e-10).map(v => v * ysc); } catch { beta = [N.mean(fy), ...subsets.map(() => 0)]; w.push('Least-squares system singular: effects could not be estimated.'); }
  const eff = subsets.map((_, i) => 2 * beta[i + 1]), tn = subsets.map(s => s.map(j => names[j]).join('×')), pred = fx.map(x => N.dot(row(x), beta)), sse = N.sum(fy.map((v, i) => (v - pred[i]) ** 2)), dof = nf - 1 - subsets.length;
  const sst = N.sum(fy.map(v => (v - N.mean(fy)) ** 2)), r2 = sst > 0 ? 1 - sse / sst : 1, mAbs = eff.map(Math.abs), med = a => N.quantile(a, 0.5);
  const s0 = 1.5 * med(mAbs), pse = 1.5 * med(mAbs.filter(v => v < 2.5 * s0).length ? mAbs.filter(v => v < 2.5 * s0) : mAbs);
  let sres = 0, seEff, me, method;
  if (dof >= 2) { sres = Math.sqrt(sse / dof); seEff = 2 * sres / Math.sqrt(nf); me = t975(dof) * seEff; method = `t-test with ${dof} residual degrees of freedom`; }
  else { seEff = pse; me = t975(Math.max(1, eff.length / 3)) * pse; method = "Lenth's pseudo standard error (unreplicated design)"; }
  const sig = eff.map(v => Math.abs(v) > me && me > 0), nSig = sig.filter(Boolean).length, ord = eff.map((_, i) => i).sort((a, b) => mAbs[b] - mAbs[a]), top = ord[0];
  if (!nSig) w.push('No effect exceeds the margin of error: within the tested ranges none of the factors has a demonstrable influence; the optimum is located with the full model and is not statistically supported.');
  const useT = nSig ? sig : sig.map(() => true), yhat = x => beta[0] + N.sum(subsets.map((s, i) => useT[i] ? beta[i + 1] * term(x, s) : 0));
  // optimum: multilinear model ⇒ extremes at vertices; target search on a grid
  const lev = inp.goal === 'target' ? N.linspace(-1, 1, 9) : [-1, 1]; let best = null, bestS = Infinity;
  const rec = (j, x) => { if (j === names.length) { const v = yhat(x), sc = inp.goal === 'max' ? -v : inp.goal === 'min' ? v : Math.abs(v - inp.yTarget); if (sc < bestS) { bestS = sc; best = [...x]; } return; } if (!act[j]) { x[j] = 0; rec(j + 1, x); return; } for (const l of lev) { x[j] = l; rec(j + 1, x); } }; rec(0, names.map(() => 0));
  const yOpt = yhat(best), xAct = best.map((v, j) => (hi[j] + lo[j]) / 2 + v * (hi[j] - lo[j]) / 2), curv = cy.length ? N.mean(cy) - beta[0] : 0, seCurv = cy.length && sres > 0 ? sres * Math.sqrt(1 / nf + 1 / cy.length) : 0, curvSig = cy.length > 0 && seCurv > 0 && Math.abs(curv) > t975(Math.max(dof, 1)) * seCurv;
  if (curvSig) w.push(`Centre points deviate from the factorial average by ${F(curv, 3)} (significant): the response is curved inside the region, so a two-level model cannot locate an interior optimum; augment to a central composite design.`);
  if (!cy.length) w.push('No centre points in the table: curvature cannot be tested and the linear-plus-interaction model is assumed.');
  // generated plans
  let fs = parseFactors(inp.factors); if (fs.length < 2) { fs = fac.map(j => ({ name: names[j], lo: lo[j], hi: hi[j] })); if (fs.length < 2) fs = [{ name: 'A', lo: -1, hi: 1 }, { name: 'B', lo: -1, hi: 1 }]; w.push('Factor definition text could not be read (need at least two "name:low:high" entries): plans were generated from the factor ranges of the table.'); }
  const kf = fs.length, rnd = N.rng(Math.round(inp.seed) + 26), actual = c => c.map((v, j) => +((fs[j].hi + fs[j].lo) / 2 + v * (fs[j].hi - fs[j].lo) / 2).toPrecision(6)), nC = Math.round(N.clamp(inp.nCentre, 0, 10));
  const corners = Array.from({ length: 1 << kf }, (_, i) => fs.map((_, j) => (i >> j) & 1 ? 1 : -1)), centre = Array.from({ length: nC }, () => fs.map(() => 0)), alpha = Math.pow(1 << kf, 0.25);
  const axial = fs.flatMap((_, j) => [-alpha, alpha].map(a => fs.map((_, q) => q === j ? a : 0))), shuffle = a => { const b = [...a]; for (let i = b.length - 1; i > 0; i--) { const q = rnd.int(i + 1);[b[i], b[q]] = [b[q], b[i]]; } return b; };
  const plan = (pts, nm) => csv(['run', ...fs.map(f => f.name), ...fs.map(f => f.name + '_coded'), 'point_type'], shuffle(pts).map((p, i) => [i + 1, ...actual(p.c), ...p.c.map(v => +v.toFixed(4)), p.t]));
  const ff = [...corners.map(c => ({ c, t: 'factorial' })), ...centre.map(c => ({ c, t: 'centre' }))], ccd = [...ff, ...axial.map(c => ({ c, t: 'axial' }))], nL = Math.round(N.clamp(inp.nLhs, 4, 200)), lh = N.lhs(nL, kf, rnd).map(u => ({ c: u.map(v => 2 * v - 1), t: 'lhs' }));
  const hq = ord.slice().reverse().map((i, r) => [N.normInv(0.5 + 0.5 * (r + 0.5) / eff.length), mAbs[i]]);
  const eqText = `${yName} = ${F(beta[0], 4)}` + subsets.map((s, i) => useT[i] ? ` ${beta[i + 1] < 0 ? '−' : '+'} ${F(Math.abs(beta[i + 1]), 4)}·${s.map(j => 'x_' + names[j]).join('·')}` : '').join('');
  return {
    summary: `Two-level factorial analysis of ${nf} factorial runs (${fac.length} factors${cy.length ? `, ${cy.length} centre points` : ''}): the largest effect is ${tn[top]} = ${F(eff[top], 4)} ${yName}; ${nSig} of ${eff.length} effects exceed the margin of error ${F(me, 3)} (${method}). The model explains R² = ${F(r2, 4)} of the variation. The predicted ${inp.goal === 'max' ? 'maximum' : inp.goal === 'min' ? 'minimum' : 'closest approach to the target'} inside the tested region is ${F(yOpt, 4)} at ${names.map((nm, j) => `${nm} = ${F(xAct[j], 4)}`).join(', ')}.`,
    kpis: [
      kp('top', `Largest effect (${tn[top]})`, eff[top], yName, `Change of the mean response when ${tn[top]} goes from its low to its high level (twice the regression coefficient in coded units).`, null),
      kp('nsig', 'Significant effects', nSig, `of ${eff.length}`, `Number of main and interaction effects larger than the margin of error, judged by ${method}.`, nSig > 0),
      kp('me', 'Margin of error of an effect', me, yName, 'Half-width of the 95 % interval for a single effect; effects smaller than this cannot be distinguished from experimental noise.', null),
      kp('r2', 'R² of the factorial model', r2, '–', 'Fraction of the response variation over the factorial runs explained by the fitted terms; with a saturated model it is 1 by construction and carries no information.', dof > 0 ? r2 > 0.8 : null),
      kp('sres', dof >= 2 ? 'Residual standard deviation' : 'Pseudo standard error (Lenth)', dof >= 2 ? sres : pse, yName, dof >= 2 ? `Run-to-run standard deviation estimated from ${dof} residual degrees of freedom (replicates and omitted interactions).` : 'Robust estimate of the standard error of an effect from the smaller effects, used when no replicates are available.', null),
      kp('yopt', 'Predicted optimum response', yOpt, yName, `Response predicted by the ${nSig ? 'significant-terms' : 'full'} model at the best combination of factor levels inside the tested region; it must be confirmed by a verification run.`, null),
      kp('curv', 'Curvature (centre − factorial mean)', curv, yName, cy.length ? `Difference between the mean of ${cy.length} centre points and the factorial average; ${curvSig ? 'significant, so quadratic terms are needed' : 'not significant at the 5 % level'}.` : 'No centre points were supplied, so curvature is not estimated and is reported as zero.', cy.length ? !curvSig : null),
      kp('mean', 'Mean response (intercept)', beta[0], yName, 'Average response over the factorial points, the intercept of the coded regression model.', null),
    ],
    series: [
      { key: 'half', title: 'Half-normal plot of effects', xlabel: 'half-normal quantile', ylabel: '|effect|', lines: [{ name: 'effects', x: hq.map(p => p[0]), y: hq.map(p => p[1]) }, { name: 'noise line (PSE)', x: [0, hq.at(-1)[0]], y: [0, pse * hq.at(-1)[0]] }, { name: 'margin of error', x: [0, hq.at(-1)[0]], y: [me, me] }], explain: 'Absolute effects against half-normal quantiles. Inactive effects fall on a straight line through the origin; effects well above the line and above the margin of error are real.' },
      { key: 'main', title: 'Main-effects plot', xlabel: 'coded level', ylabel: `mean ${yName}`, lines: fac.map(j => ({ name: names[j], x: [-1, 1], y: [-1, 1].map(l => N.mean(fy.filter((_, i) => Math.sign(fx[i][j]) === l))) })), explain: 'Average response at the low and high level of each factor. Steep lines mark influential factors; the plot ignores interactions, which must be read from the effects table.' },
      { key: 'par', title: 'Predicted versus observed', xlabel: `observed ${yName}`, ylabel: `predicted ${yName}`, lines: [{ name: 'factorial runs', x: fy, y: fx.map(yhat) }, { name: 'perfect agreement', x: [N.minOf(fy), N.maxOf(fy)], y: [N.minOf(fy), N.maxOf(fy)] }], explain: 'Model predictions (significant terms only) against the measured responses. Scatter about the diagonal reflects experimental noise and omitted terms.' },
    ],
    tables: [
      { key: 'eff', title: 'Estimated effects (ranked)', columns: ['Term', 'Effect', 'Coefficient (coded)', 'Std. error', '|Effect|/margin', 'Significant'], rows: ord.map(i => [tn[i], F(eff[i], 4), F(beta[i + 1], 4), F(seEff, 3), F(me > 0 ? mAbs[i] / me : 0, 3), sig[i] ? 'yes' : 'no']), explain: 'Main and interaction effects ranked by magnitude with their standard error. A significant interaction means the effect of one factor depends on the level of another, so the factors cannot be optimised one at a time.' },
      { key: 'opt', title: 'Predicted optimum inside the tested region', columns: ['Factor', 'Low', 'High', 'Coded setting', 'Actual setting'], rows: names.map((nm, j) => [nm, lo[j], hi[j], +best[j].toFixed(3), F(xAct[j], 5)]), explain: 'Factor settings that give the best predicted response. A setting on the boundary of the region indicates that the true optimum may lie outside it, in the direction of steepest ascent.' },
    ],
    equations: [
      { kind: 'constitutive', label: 'Fitted regression model (coded units, significant terms)', text: eqText, note: 'x = (level − centre)/half-range, so each x runs from −1 to +1.' },
      { kind: 'formula', label: 'Effect estimate', text: 'effect = ȳ(+) − ȳ(−) = 2·β;   β = (XᵀX)⁻¹Xᵀy' },
      { kind: 'formula', label: "Lenth's pseudo standard error", text: 's₀ = 1.5·median|eⱼ|,  PSE = 1.5·median{|eⱼ| : |eⱼ| < 2.5·s₀},  ME = t₀.₉₇₅,m/3·PSE' },
    ],
    checks: [
      { name: 'Complete two-level design', pass: distinct === (1 << kk), detail: `${distinct} of ${1 << kk} factorial points present` },
      { name: 'Balanced design (orthogonal effect estimates)', pass: fac.every(j => Math.abs(N.sum(fx.map(x => x[j]))) < 1e-9), detail: 'sum of coded levels per factor is zero when balanced' },
      { name: 'Independent error estimate available', pass: dof >= 2, detail: dof >= 2 ? `${dof} residual degrees of freedom` : 'unreplicated: pseudo standard error used' },
      { name: 'No significant curvature', pass: !curvSig, detail: cy.length ? `curvature ${F(curv, 3)} ± ${F(seCurv, 3)}` : 'no centre points: not tested' },
      { name: 'At least one effect distinguishable from noise', pass: nSig > 0, detail: `${nSig} significant; margin of error ${F(me, 3)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: 'high', text: `Run a confirmation experiment at ${names.map((nm, j) => `${nm} = ${F(xAct[j], 4)}`).join(', ')}; the model predicts ${F(yOpt, 4)} ${yName} (± about ${F(me, 3)}). Agreement confirms the model, disagreement points to curvature or an uncontrolled factor.`, basis: 'Predicted optimum of the fitted factorial model.' },
      { priority: curvSig ? 'high' : 'medium', text: curvSig ? 'Add the axial runs of the central composite plan (download below) to the existing factorial and centre points and fit a quadratic response surface before choosing the operating point.' : best.some((v, j) => act[j] && Math.abs(v) === 1) ? `The optimum lies on the edge of the region for ${names.filter((_, j) => act[j] && Math.abs(best[j]) === 1).join(', ')}: move the next design along the path of steepest ${inp.goal === 'min' ? 'descent' : 'ascent'} in proportion to the coefficients.` : 'The target is met inside the region; tighten the factor ranges around this point for a robustness study.', basis: 'Curvature test and location of the predicted optimum.' },
      { priority: 'medium', text: nSig < eff.length ? `Factors and interactions without significant effect (${ord.filter(i => !sig[i]).slice(0, 4).map(i => tn[i]).join(', ') || 'none'}) can be set to their most economical level: shortest time, lowest temperature and lowest precursor concentration reduce energy and material use at no cost in the response.` : 'All tested effects are active; none of the factors can be relaxed for cost without changing the response.', basis: 'Effects below the margin of error.', sustainability: true },
      { priority: 'low', text: `A full factorial for ${kf} factors needs ${1 << kf} runs (+${nC} centre points), the central composite ${ccd.length}, the Latin hypercube ${nL}; randomise the run order as given in the downloaded plans to protect against drift.`, basis: 'Generated experimental plans.' },
    ],
    publish: { 'opt.objective': yOpt, 'doe.r2': r2, 'doe.n_significant': nSig, 'doe.top_effect': eff[top] },
    exports: [
      { name: 'plan_full_factorial.csv', mime: 'text/csv', content: plan(ff), description: `2^${kf} full factorial with ${nC} centre points in randomised run order (actual and coded levels).` },
      { name: 'plan_central_composite.csv', mime: 'text/csv', content: plan(ccd), description: `Rotatable central composite design (α = ${F(alpha, 4)}); axial points lie outside the low/high range and must be checked for feasibility.` },
      { name: 'plan_latin_hypercube.csv', mime: 'text/csv', content: plan(lh), description: `Space-filling Latin-hypercube plan with ${nL} runs for surrogate modelling.` },
    ],
    numerics: { method: 'Least-squares estimation of two-level factorial effects; Lenth pseudo standard error or residual t-test; exhaustive vertex/grid search for the optimum', dof: subsets.length + 1, steps: nf, iterations: Math.pow(lev.length, fac.length), residual: Math.sqrt(sse / Math.max(nf, 1)), notes: 'Linear-plus-interaction model in coded units; valid only inside the tested factor ranges.' },
  };
}

// ---------- mode: yield, throughput, OEE ----------
function runYield(inp) {
  const w = [], t = readRows(inp.steps, 5, 1, 'Process step table'), S = t.rows.map((r, i) => ({ name: r[0] || `step ${i + 1}`, ct: Math.max(+r[1], 1e-6), mach: Math.max(Math.round(+r[2]) || 1, 1), y: N.clamp(+r[3] / 100, 1e-6, 1), av: N.clamp(+r[4] / 100, 1e-3, 1) })).filter(s => Number.isFinite(s.ct) && Number.isFinite(s.y) && Number.isFinite(s.av));
  if (!S.length) throw new Error('Process step table: no row with numeric cycle time, machines, yield and availability.');
  if ([STEPS_DEF, STEPS_ZNO, STEPS_AM].includes(inp.steps)) w.push('Process step table: illustrative demonstration data — replace with the cycle times, yields and availabilities of your line.');
  const AD = inp.D0 * inp.area * inp.critFrac, al = inp.alphaC, Ym = { poisson: x => Math.exp(-x), murphy: x => x < 1e-8 ? 1 : ((1 - Math.exp(-x)) / x) ** 2, seeds: x => 1 / (1 + x), negbin: x => Math.pow(1 + x / al, -al) }, Yd = (Ym[inp.ymodel] || Ym.negbin)(AD);
  const rty = S.reduce((p, s) => p * s.y, 1), yTot = rty * Yd, perf = inp.perf / 100, hrs = inp.hours;
  const down = S.map((_, i) => S.slice(i).reduce((p, s) => p * s.y, 1) * Yd), cap = S.map(s => s.mach * s.av * perf * 60 / s.ct), good = cap.map((c, i) => c * down[i]), bn = N.argmin(good), lineCap = good[bn];
  const need = inp.units / hrs, util = good.map(g => need / g), takt = hrs * 60 / inp.units, annual = lineCap * hrs, oee = S[bn].av * perf * yTot, starts = inp.units / yTot;
  const camp = inp.nPerCamp * (inp.tBatch + inp.tClean) + inp.tChange, nCamp = hrs / camp, batches = nCamp * inp.nPerCamp, campCap = batches * inp.batchSize * yTot, upFrac = inp.nPerCamp * inp.tBatch / camp;
  if (util[bn] > 1) w.push(`The bottleneck "${S[bn].name}" is loaded to ${pct(util[bn])}: the line cannot meet the demand of ${F(inp.units, 4)} good units per year (capacity ${F(annual, 4)}).`);
  if (AD > 3) w.push(`Mean number of killer defects per unit A·D = ${F(AD, 3)}: the yield models diverge strongly here and the result is dominated by the assumed defect clustering.`);
  const As = N.logspace(Math.log10(inp.area / 100), Math.log10(inp.area * 20), 50), Ds = N.logspace(Math.log10(inp.D0 / 50), Math.log10(inp.D0 * 20), 50), mods = [['Poisson', Ym.poisson], ['Murphy', Ym.murphy], ['Seeds', Ym.seeds], [`Negative binomial (α = ${al})`, Ym.negbin]];
  const worstY = N.argmin(S.map(s => s.y)), x90 = inp.ymodel === 'poisson' ? -Math.log(0.9) : inp.ymodel === 'seeds' ? 1 / 0.9 - 1 : inp.ymodel === 'murphy' ? N.brent(x => Ym.murphy(x) - 0.9, 1e-6, 2) : al * (Math.pow(0.9, -1 / al) - 1), dTarget = x90 / (inp.area * inp.critFrac);   // A·D₀ at which the selected model gives Y = 0.9, converted to a defect density
  return {
    summary: `Defect-limited yield is ${pct(Yd)} (${inp.ymodel} model, ${F(AD, 3)} killer defects per unit) and the rolled throughput yield of ${S.length} steps is ${pct(rty)}, giving ${pct(yTot)} overall; ${F(starts, 4)} starts are needed for ${F(inp.units, 4)} good units per year. The takt time is ${F(takt, 3)} min; the bottleneck is "${S[bn].name}" at ${pct(util[bn])} utilisation, limiting the line to ${F(lineCap, 3)} good units/h (${F(annual, 4)} per year). OEE at the bottleneck is ${pct(oee)}; the batch campaign schedule yields ${F(campCap, 4)} good kg or units per year.`,
    kpis: [
      kp('y_def', 'Defect-limited yield', Yd, '–', `Probability that a unit of ${inp.area} cm² (critical fraction ${inp.critFrac}) carries no killer defect at ${inp.D0} defects/cm², from the ${inp.ymodel} model.`, Yd > 0.8),
      kp('ad', 'Killer defects per unit (A·D₀)', AD, '–', 'Mean number of yield-killing defects falling on the critical area of one unit; the single parameter that controls every area-yield model.', AD < 0.3),
      kp('rty', 'Rolled throughput yield', rty, '–', `Product of the first-pass yields of all ${S.length} process steps: the probability that a unit passes every step without rework or scrap.`, rty > 0.85),
      kp('y_tot', 'Overall yield', yTot, '–', 'Rolled throughput yield multiplied by the defect-limited yield: good units out per unit started.', yTot > 0.7),
      kp('takt', 'Takt time', takt, 'min', `Available production time divided by demand: one good unit must leave the line every ${F(takt, 3)} min to meet ${F(inp.units, 4)} units per year.`, null),
      kp('util', `Bottleneck utilisation (${S[bn].name})`, util[bn], '–', 'Required rate divided by the effective capacity of the most constrained step, allowing for downstream yield losses, availability and performance; above about 0.85 queues grow rapidly, above 1 demand cannot be met.', util[bn] < 0.85),
      kp('cap', 'Line capacity (good units)', lineCap, '1/h', 'Good units per operating hour deliverable by the bottleneck step after downstream yield losses.', lineCap >= need),
      kp('annual', 'Annual capacity (good units)', annual, '1/yr', `Line capacity multiplied by ${inp.hours} operating hours per year, to be compared with the demand of ${F(inp.units, 4)}.`, annual >= inp.units),
      kp('oee', 'Overall equipment effectiveness', oee, '–', `Availability × performance × quality at the bottleneck (${F(S[bn].av, 3)} × ${F(perf, 3)} × ${F(yTot, 3)}); world-class discrete manufacturing is usually quoted around 0.85.`, oee > 0.6),
      kp('scrap', 'Scrap rate', 1 - yTot, '–', 'Fraction of started units lost to step failures and killer defects; each carries the material and energy already invested up to the point of loss.', 1 - yTot < 0.2),
      kp('camp', 'Campaign annual output', campCap, 'kg or units/yr', `Good output of the batch campaign schedule: ${F(batches, 4)} batches per year in ${F(nCamp, 3)} campaigns; ${pct(upFrac)} of scheduled time is productive.`, null),
      kp('d_target', 'Defect density for 90 % defect yield', dTarget, '1/cm²', 'Defect density at which the selected yield model predicts 90 % defect-limited yield for this unit area; the cleanliness target for the patterning environment.', inp.D0 <= dTarget),
    ],
    series: [
      { key: 'ya', title: 'Yield versus unit area', xlabel: 'unit area (cm²)', ylabel: 'defect-limited yield', logx: true, lines: [...mods.map(([nm, fn]) => ({ name: nm, x: As, y: As.map(a => fn(inp.D0 * a * inp.critFrac)) })), { name: 'design point', x: [inp.area], y: [Yd] }], explain: 'Defect-limited yield against device or panel area for the four classical models at the stated defect density. Large-area metasurfaces and lattice panels fall on the steep part of the curve, where clustering assumptions matter most.' },
      { key: 'yd', title: 'Yield versus defect density', xlabel: 'defect density (1/cm²)', ylabel: 'defect-limited yield', logx: true, lines: mods.map(([nm, fn]) => ({ name: nm, x: Ds, y: Ds.map(d => fn(d * inp.area * inp.critFrac)) })), explain: 'Sensitivity of yield to the cleanliness of the process for the present unit area. It shows what defect density must be reached for an economically acceptable yield.' },
      { key: 'util', title: 'Utilisation by process step', xlabel: 'step number', ylabel: 'utilisation', lines: [{ name: 'utilisation at demand', x: S.map((_, i) => i + 1), y: util }, { name: 'full load', x: [1, Math.max(S.length, 2)], y: [1, 1] }], explain: 'Load of each step at the demanded output, with yield losses downstream taken into account. The highest point is the bottleneck; steps far below it have spare capacity.' },
    ],
    tables: [
      { key: 'steps', title: 'Step capacity and loading', columns: ['#', 'Step', 'Cycle (min)', 'Machines', 'Yield (%)', 'Availability (%)', 'Capacity (units/h)', 'Good-unit capacity (1/h)', 'Utilisation'], rows: S.map((s, i) => [i + 1, s.name + (i === bn ? ' (bottleneck)' : ''), s.ct, s.mach, +(100 * s.y).toFixed(2), +(100 * s.av).toFixed(1), F(cap[i], 4), F(good[i], 4), F(util[i], 3)]), explain: 'Effective capacity of each step (machines × availability × performance ÷ cycle time) and the good-unit rate it can support after later yield losses. The smallest good-unit capacity sets the line output.' },
      { key: 'ym', title: 'Yield models compared at the design point', columns: ['Model', 'Formula', 'Yield'], rows: [['Poisson', 'exp(−A·D₀)', F(Ym.poisson(AD), 4)], ['Murphy', '[(1 − e^(−A·D₀))/(A·D₀)]²', F(Ym.murphy(AD), 4)], ['Seeds', '1/(1 + A·D₀)', F(Ym.seeds(AD), 4)], ['Negative binomial', '(1 + A·D₀/α)^(−α)', F(Ym.negbin(AD), 4)]], explain: 'The four area-yield models evaluated with the same defect density and critical area. Their spread is the model-form uncertainty of the yield forecast until clustering is measured.' },
    ],
    equations: [
      { kind: 'formula', label: 'Defect-limited yield models', text: 'Y_Poisson = e^(−AD₀);  Y_Murphy = [(1 − e^(−AD₀))/(AD₀)]²;  Y_Seeds = 1/(1 + AD₀);  Y_NB = (1 + AD₀/α)^(−α)' },
      { kind: 'formula', label: 'Rolled throughput yield', text: 'RTY = Πᵢ yᵢ;   Y_total = RTY·Y_defect;   starts = demand/Y_total' },
      { kind: 'formula', label: 'Takt time, capacity and OEE', text: 'takt = t_available/demand;  Cᵢ = mᵢ·aᵢ·p/t_c,ᵢ;  line capacity = minᵢ(Cᵢ·Πⱼ≥ᵢ yⱼ);  OEE = availability·performance·quality' },
      { kind: 'formula', label: 'Campaign capacity', text: 'T_campaign = n·(t_batch + t_clean) + t_changeover;  output = (t_available/T_campaign)·n·batch size·Y_total' },
    ],
    checks: [
      { name: 'Demand within line capacity', pass: annual >= inp.units, detail: `capacity ${F(annual, 4)} vs demand ${F(inp.units, 4)} good units/yr` },
      { name: 'Bottleneck utilisation below 85 %', pass: util[bn] < 0.85, detail: `${S[bn].name}: ${pct(util[bn])}` },
      { name: 'Yield probabilities within [0, 1]', pass: Yd >= 0 && Yd <= 1 && rty > 0 && rty <= 1, detail: `Y_defect = ${F(Yd, 4)}, RTY = ${F(rty, 4)}` },
      { name: 'Model ordering Poisson ≤ Murphy ≤ Seeds', pass: Ym.poisson(AD) <= Ym.murphy(AD) + 1e-12 && Ym.murphy(AD) <= Ym.seeds(AD) + 1e-12, detail: `${F(Ym.poisson(AD), 4)} ≤ ${F(Ym.murphy(AD), 4)} ≤ ${F(Ym.seeds(AD), 4)}` },
      { name: 'Unit flow balance (starts × yield = good units)', pass: Math.abs(starts * yTot - inp.units) < 1e-9 * inp.units, detail: `${F(starts, 5)} starts × ${F(yTot, 4)}` },
    ],
    warnings: w,
    recommendations: [
      { priority: util[bn] > 0.85 ? 'high' : 'medium', text: util[bn] > 0.85 ? `Relieve the bottleneck "${S[bn].name}": ${Math.ceil(S[bn].mach * util[bn] / 0.85)} machines instead of ${S[bn].mach}, or a cycle time of ${F(S[bn].ct * 0.85 / util[bn], 3)} min, brings utilisation to 85 %.` : `The line has ${pct(1 - util[bn])} head-room at the bottleneck; demand can grow to about ${F(annual * 0.85, 4)} units per year before the next investment.`, basis: 'Step utilisation with downstream yield losses.' },
      { priority: Yd < 0.8 ? 'high' : 'low', text: Yd < 0.8 ? `Reduce the killer-defect density from ${inp.D0} to ${F(dTarget, 3)} cm⁻² (cleaner environment, filtered resist or resin, in-line inspection) or tile the product into smaller units: yield falls exponentially with critical area.` : 'Defect-limited yield is not the main loss; concentrate on the step yields.', basis: 'Area-yield model at the design point.' },
      { priority: 'medium', text: `The lowest step yield is "${S[worstY].name}" (${pct(S[worstY].y)}); raising it by one percentage point increases good output by ${pct(0.01 / S[worstY].y)} at unchanged cost, and moving inspection ahead of the most expensive step avoids adding value to units that are already defective.`, basis: 'Rolled throughput yield decomposition.' },
      { priority: 'medium', text: `A scrap rate of ${pct(1 - yTot)} means ${F(starts - inp.units, 4)} units per year are produced and discarded, with their precursors, solvents, substrate and energy. Rework loops, recovery of noble-metal or substrate material from scrap and earlier defect detection reduce this waste directly.`, basis: 'Starts required versus good units delivered.', sustainability: true },
      { priority: 'low', text: `Productive time in the campaign schedule is ${pct(upFrac)}; lengthening campaigns from ${inp.nPerCamp} to ${inp.nPerCamp * 2} batches raises it to ${pct(2 * inp.nPerCamp * inp.tBatch / (2 * inp.nPerCamp * (inp.tBatch + inp.tClean) + inp.tChange))} at the cost of inventory.`, basis: 'Campaign cycle-time accounting.' },
    ],
    publish: { 'fab.yield': yTot, 'fab.throughput_units_h': lineCap, 'fab.oee': oee, 'fab.rty': rty, 'fab.takt_min': takt },
    numerics: { method: 'Closed-form area-yield models, serial-line capacity analysis and campaign time accounting', dof: S.length, steps: S.length, iterations: 1, residual: 0, notes: 'Deterministic capacity model without queueing variability; real lines need buffer analysis above about 85 % utilisation.' },
  };
}

// ---------- mode: automation readiness ----------
function runAuto(inp) {
  const w = [], tau = inp.tauP, th = inp.deadT, K = inp.Kproc, tc = Math.max(inp.lam * th, 0.05 * tau), Kc = tau / (K * (tc + th)), Ti = Math.min(tau, 4 * (tc + th)), Ts = inp.Ts;
  const fMax = 1 / (2 * Math.PI * inp.tauFast), TsNyq = 1 / (2 * fMax), TsRec = Math.min(inp.tauFast / 10, tau / 10, Math.max(th, tc) / 5 || Infinity);
  // loop frequency response: gain crossover and phase margin (continuous loop, then sample-and-hold delay Ts/2)
  const Lmag = om => Kc * K * Math.sqrt(1 + 1 / (Ti * om) ** 2) / Math.sqrt(1 + (tau * om) ** 2); let wc = 1 / tau;
  try { wc = Math.exp(N.brent(x => Math.log(Lmag(Math.exp(x))), Math.log(1e-6 / tau), Math.log(1e6 / Math.min(tau, tc + th + 1e-9)))); } catch { w.push('Gain crossover could not be bracketed; phase margin evaluated at ω = 1/τ.'); }
  const ph = om => -Math.PI / 2 + Math.atan(Ti * om) - Math.atan(tau * om) - th * om, pm = (Math.PI + ph(wc)) * 180 / Math.PI, pmS = pm - wc * Ts / 2 * 180 / Math.PI;
  // discrete closed-loop simulation: exact zero-order-hold step of the first-order lag, PI updated every Ts
  // (a sampling interval below the simulation step is treated as continuous control)
  const tEnd = 12 * (tc + th) + 4 * tau, dt0 = tEnd / 4000, fastS = Ts <= dt0, sub = fastS ? 1 : Math.ceil(Ts / dt0), dt = fastS ? dt0 : Ts / sub, TsC = fastS ? dt0 : Ts, nSt = Math.min(Math.ceil(tEnd / dt), 8000), a = Math.exp(-dt / tau), dly = Math.min(Math.round(th / dt), 8000), buf = new Float64Array(dly + 1);
  let y = 0, u = 0, integ = 0, over = 0, tset = 0, iae = 0, stable = true; const tt = [], yy = [], uu = [];
  for (let i = 0; i < nSt; i++) { if (i % sub === 0) { const e = 1 - y; integ += e * TsC / Ti; u = Kc * (e + integ); } const bi = i % (dly + 1), ud = buf[bi]; buf[bi] = u; y = a * y + K * (1 - a) * ud; if (!Number.isFinite(y) || Math.abs(y) > 1e6) { stable = false; y = N.clamp(y || 0, -1e6, 1e6); integ = 0; } const tNow = (i + 1) * dt; over = Math.max(over, y - 1); if (Math.abs(y - 1) > 0.02) tset = tNow; iae += Math.abs(1 - y) * dt; if (i % Math.ceil(nSt / 300) === 0) { tt.push(tNow); yy.push(y); uu.push(N.clamp(u * K, -1e6, 1e6)); } }
  const settled = stable && tset < 0.98 * nSt * dt; if (!settled) w.push(`The simulated loop does not settle within ${F(nSt * dt, 3)} s at a sampling interval of ${F(Ts, 3)} s: the sampling is too slow or the tuning too tight for this dead time.`);
  if (Ts > TsNyq) w.push(`Sampling interval ${F(Ts, 3)} s violates the Nyquist condition for the fastest phenomenon (τ = ${F(inp.tauFast, 3)} s needs Ts < ${F(TsNyq, 3)} s): fast events are aliased and cannot be reconstructed.`);
  else if (Ts > inp.tauFast / 10) w.push(`Sampling interval ${F(Ts, 3)} s satisfies Nyquist but not the ten-samples-per-time-constant rule (${F(inp.tauFast / 10, 3)} s): the shape of fast transients is poorly resolved.`);
  const bytes = inp.chan * 86400 / Ts * 8;
  // gauge R&R (crossed two-way ANOVA with interaction)
  const mt = readRows(inp.msa, 3, 8, 'Gauge study table'), obs = mt.rows.map(r => [r[0], r[1], Number(r[2])]).filter(r => Number.isFinite(r[2])); if (obs.length < 8) throw new Error('Gauge study table: at least 8 numeric measurements are needed.');
  if (inp.msa === MSA_DEF) w.push('Gauge study table: synthetic demonstration data — replace with measurements.');
  const parts = [...new Set(obs.map(r => r[0]))], ops = [...new Set(obs.map(r => r[1]))], cell = new Map(); obs.forEach(([p, o, v]) => { const key = p + '|' + o; if (!cell.has(key)) cell.set(key, []); cell.get(key).push(v); });
  const p = parts.length, o = ops.length, rr = Math.min(...parts.flatMap(pp => ops.map(oo => (cell.get(pp + '|' + oo) || []).length)));
  let grr = { ok: false, ev: 0, av: 0, pv: 0, tot: 0 };
  if (p >= 2 && rr >= 2) {
    const get = (pp, oo) => cell.get(pp + '|' + oo).slice(0, rr), allv = parts.flatMap(pp => ops.flatMap(oo => get(pp, oo))), g = N.mean(allv), pm_ = parts.map(pp => N.mean(ops.flatMap(oo => get(pp, oo)))), om = ops.map(oo => N.mean(parts.flatMap(pp => get(pp, oo))));
    const ssP = o * rr * N.sum(pm_.map(v => (v - g) ** 2)), ssO = p * rr * N.sum(om.map(v => (v - g) ** 2)); let ssE = 0, ssC = 0; parts.forEach((pp, i) => ops.forEach((oo, j) => { const c = get(pp, oo), cm = N.mean(c); ssC += rr * (cm - pm_[i] - om[j] + g) ** 2; ssE += N.sum(c.map(v => (v - cm) ** 2)); }));
    const msP = ssP / (p - 1), msO = o > 1 ? ssO / (o - 1) : 0, msPO = o > 1 ? ssC / ((p - 1) * (o - 1)) : 0, msE = ssE / (p * o * (rr - 1));
    const vE = msE, vPO = o > 1 ? Math.max(0, (msPO - msE) / rr) : 0, vO = o > 1 ? Math.max(0, (msO - msPO) / (p * rr)) : 0, vP = Math.max(0, (msP - (o > 1 ? msPO : msE)) / (o * rr));
    grr = { ok: true, ev: Math.sqrt(vE), av: Math.sqrt(vO + vPO), pv: Math.sqrt(vP), tot: Math.sqrt(vE + vO + vPO + vP), inter: Math.sqrt(vPO) };
    if (cell.size !== p * o || [...cell.values()].some(c => c.length !== rr)) w.push(`Gauge study is unbalanced: analysed with the first ${rr} replicate(s) of every part–operator cell.`);
    if (o < 2) w.push('Only one operator in the gauge study: reproducibility cannot be estimated and is taken as zero.');
  } else w.push('Gauge study needs at least two parts and two replicates per part–operator cell: gauge R&R is not evaluated and reported as zero.');
  const sg = Math.hypot(grr.ev, grr.av), pGrr = grr.tot > 0 ? sg / grr.tot : 0, tol = Math.abs(inp.usl - inp.lsl), pTol = tol > 0 ? 6 * sg / tol : 0, ndc = sg > 0 ? Math.max(1, Math.floor(1.41 * grr.pv / sg)) : 0, msaClass = !grr.ok ? 'not evaluated' : pGrr < 0.1 ? 'acceptable' : pGrr <= 0.3 ? 'marginal' : 'unacceptable';
  // anomaly detection on the process signal (individual values in time order)
  const { all } = readGroups(inp, w), nb = Math.min(Math.max(5, Math.round(inp.nBase)), Math.max(5, Math.floor(all.length * 0.6))), base = all.slice(0, nb), mu0 = N.mean(base), s0 = Math.max(N.mean(base.slice(1).map((v, i) => Math.abs(v - base[i]))) / 1.128, 1e-12 * Math.max(Math.abs(mu0), 1));
  const lamE = inp.lamE, zE = [], ue = [], le = []; let zz = mu0, ewAlarm = 0, cp = 0, cm = 0, cuAlarm = 0; const cpS = [], cmS = [];
  all.forEach((v, i) => { zz = lamE * v + (1 - lamE) * zz; const hw = inp.LE * s0 * Math.sqrt(lamE / (2 - lamE) * (1 - Math.pow(1 - lamE, 2 * (i + 1)))); zE.push(zz); ue.push(mu0 + hw); le.push(mu0 - hw); if (!ewAlarm && Math.abs(zz - mu0) > hw) ewAlarm = i + 1; const zs = (v - mu0) / s0; cp = Math.max(0, cp + zs - inp.kC); cm = Math.max(0, cm - zs - inp.kC); cpS.push(cp); cmS.push(cm); if (!cuAlarm && (cp > inp.hC || cm > inp.hC)) cuAlarm = i + 1; });
  const arl = d => { const D = d - inp.kC, b = inp.hC + 1.166; return Math.abs(D) < 1e-9 ? b * b : (Math.exp(-2 * D * b) + 2 * D * b - 1) / (2 * D * D); }, arl0 = 1 / (2 / arl(0)), arl1 = 1 / (1 / arl(1) + 1 / arl(-1));
  const ix = all.map((_, i) => i + 1), thin = (xs, ys) => { const st = Math.ceil(xs.length / 380); return [xs.filter((_, i) => i % st === 0), ys.filter((_, i) => i % st === 0)]; };
  const L = (name, ys) => { const [x, yv] = thin(ix, ys); return { name, x, y: yv }; };
  return {
    summary: `To resolve the fastest phenomenon (τ = ${F(inp.tauFast, 3)} s) the sensors must sample faster than every ${F(TsNyq, 3)} s (Nyquist) and preferably every ${F(TsRec, 3)} s; the chosen ${F(Ts, 3)} s is ${Ts <= TsRec ? 'adequate' : Ts <= TsNyq ? 'marginal' : 'too slow'}. SIMC tuning of the first-order-plus-dead-time loop gives Kc = ${F(Kc, 3)}, Tᵢ = ${F(Ti, 3)} s with a phase margin of ${F(pmS, 3)}° including the sampling delay; the simulated set-point step ${settled ? `settles in ${F(tset, 3)} s with ${pct(Math.max(over, 0))} overshoot` : 'does not settle'}. Gauge R&R is ${pct(pGrr)} of total variation (${msaClass}); on the process signal the EWMA ${ewAlarm ? `alarms at observation ${ewAlarm}` : 'raises no alarm'} and the CUSUM ${cuAlarm ? `at observation ${cuAlarm}` : 'raises no alarm'}.`,
    kpis: [
      kp('ts_nyq', 'Nyquist sampling interval', TsNyq, 's', `Largest sampling interval that still satisfies the sampling theorem for a signal bandwidth of 1/(2πτ) = ${F(fMax, 3)} Hz of the fastest phenomenon; slower sampling aliases it.`, Ts <= TsNyq),
      kp('ts_rec', 'Recommended sampling interval', TsRec, 's', 'Practical requirement of about ten samples per fastest time constant and per closed-loop time constant, whichever is smaller; anti-alias filters should cut off below half the sampling rate.', Ts <= TsRec),
      kp('kc', 'Controller gain Kc', Kc, 'input/output', `Proportional gain from the SIMC rule Kc = τ/(K·(τc + θ)) with τc = ${F(tc, 3)} s for the first-order-plus-dead-time model.`, null),
      kp('ti', 'Integral time Tᵢ', Ti, 's', 'Integral (reset) time from the SIMC rule min(τ, 4(τc + θ)); the smaller branch speeds up disturbance rejection for slow processes.', null),
      kp('pm', 'Phase margin with sampling delay', pmS, '°', `Phase margin of the loop at the gain crossover ${F(wc, 3)} rad/s, reduced by the half-sample hold delay (continuous value ${F(pm, 3)}°); margins below about 30° give oscillatory control.`, pmS >= 30),
      kp('settle', 'Set-point settling time (2 %)', settled ? tset : nSt * dt, 's', settled ? 'Time after which the simulated closed-loop step response stays within ±2 % of the set point with the chosen sampling interval.' : 'The response did not settle within the simulated horizon; the horizon length is reported and the loop must be retuned.', settled),
      kp('over', 'Overshoot', 100 * Math.max(over, 0), '%', 'Peak excess of the controlled variable over the set point in the simulated step; for temperature-sensitive nucleation it should stay within a few per cent.', over < 0.1),
      kp('data', 'Raw data volume', bytes / 1e6, 'MB/day', `Storage for ${inp.chan} channels of 8-byte values at the chosen sampling interval, before compression; it sizes the historian of the digital twin.`, null),
      kp('grr', 'Gauge R&R (% of total variation)', 100 * pGrr, '%', `Combined repeatability and reproducibility standard deviation relative to the total; below 10 % is customarily acceptable, 10–30 % marginal, above 30 % unacceptable. Relative to the tolerance it is ${pct(pTol)}.`, grr.ok ? pGrr < 0.3 : null),
      kp('ndc', 'Number of distinct categories', ndc, '', 'Number of non-overlapping groups of parts the measurement system can distinguish (1.41·σ_part/σ_GRR); five or more are needed for process control.', grr.ok ? ndc >= 5 : null),
      kp('ewma', 'First EWMA alarm', ewAlarm, 'observation', ewAlarm ? `Observation at which the exponentially weighted moving average first leaves its ±${inp.LE}σ limits (baseline from the first ${nb} points).` : 'The EWMA statistic stays inside its control limits for the whole record; zero denotes no alarm.', ewAlarm === 0),
      kp('cusum', 'First CUSUM alarm', cuAlarm, 'observation', cuAlarm ? `Observation at which the tabular CUSUM exceeds h = ${inp.hC}σ; the shift began roughly ${F(Math.max(cpS[cuAlarm - 1], cmS[cuAlarm - 1]) / Math.max(inp.kC, 0.1), 2)} or fewer observations earlier.` : 'The cumulative sums stay below the decision interval for the whole record; zero denotes no alarm.', cuAlarm === 0),
      kp('arl0', 'CUSUM in-control run length', arl0, 'observations', `Average number of observations between false alarms for the two-sided CUSUM (Siegmund approximation); a sustained 1σ shift is detected after about ${F(arl1, 3)} observations.`, arl0 > 200),
    ],
    series: [
      { key: 'step', title: 'Closed-loop set-point step', xlabel: 'time (s)', ylabel: 'controlled variable (set point = 1)', lines: [{ name: 'controlled variable', x: tt, y: yy }, { name: 'controller output × process gain', x: tt, y: uu }, { name: 'set point', x: [0, tt.at(-1)], y: [1, 1] }], explain: 'Simulated response of the sampled PI loop to a unit set-point change. Look for overshoot, oscillation and the size of the initial controller move, which must stay within the actuator range.' },
      { key: 'ewma', title: 'EWMA chart of the process signal', xlabel: 'observation', ylabel: 'value', lines: [L('EWMA statistic', zE), L('upper limit', ue), L('lower limit', le), L('raw signal', all)], explain: 'Exponentially weighted moving average with time-varying limits. It smooths noise and reveals small sustained shifts earlier than a Shewhart chart.' },
      { key: 'cusum', title: 'Tabular CUSUM', xlabel: 'observation', ylabel: 'cumulative sum (σ units)', lines: [L('C⁺ (upward shift)', cpS), L('C⁻ (downward shift)', cmS), { name: 'decision interval h', x: [1, all.length], y: [inp.hC, inp.hC] }], explain: 'Upper and lower cumulative sums of standardised deviations beyond the slack k. A steadily rising branch marks the onset of a shift; the alarm is raised when it crosses h.' },
    ],
    tables: [
      { key: 'msa', title: 'Gauge R&R variance components', columns: ['Source', 'Std. deviation', '% of total variation', '% of tolerance (6σ)'], rows: [['Repeatability (equipment)', grr.ev, grr.ev], ['Reproducibility (operator + interaction)', grr.av, grr.av], ['Gauge R&R', sg, sg], ['Part-to-part', grr.pv, grr.pv], ['Total', grr.tot, grr.tot]].map(([nm, s]) => [nm, F(s, 4), F(grr.tot > 0 ? 100 * s / grr.tot : 0, 3), F(tol > 0 ? 600 * s / tol : 0, 3)]), explain: `Variance components from the crossed ANOVA of ${p} parts × ${o} operators × ${rr} replicates. Repeatability belongs to the instrument, reproducibility to operators and procedure; percentages are ratios of standard deviations and do not add to 100.` },
      { key: 'samp', title: 'Sampling and data requirements', columns: ['Requirement', 'Value'], rows: [['Bandwidth of fastest phenomenon (Hz)', F(fMax, 4)], ['Nyquist minimum sampling rate (Hz)', F(2 * fMax, 4)], ['Ten-samples-per-time-constant rate (Hz)', F(10 / inp.tauFast, 4)], ['Chosen sampling rate (Hz)', F(1 / Ts, 4)], ['Anti-alias filter cut-off (Hz), ≤ 0.4 × sampling rate', F(0.4 / Ts, 4)], ['Closed-loop time constant τc (s)', F(tc, 4)], ['Samples per closed-loop time constant', F(tc / Ts, 4)]], explain: 'Sampling rates derived from the process time constants compared with the chosen rate. The digital twin can only reproduce dynamics that the sensors resolve.' },
    ],
    equations: [
      { kind: 'governing', label: 'Process model (first order plus dead time)', text: 'τ·dy/dt + y = K·u(t − θ)   ⇔   G(s) = K·e^(−θs)/(τs + 1)' },
      { kind: 'constitutive', label: 'PI controller and SIMC tuning', text: 'u = Kc·[e + (1/Tᵢ)∫e dt];   Kc = τ/(K(τc + θ)),  Tᵢ = min(τ, 4(τc + θ))' },
      { kind: 'formula', label: 'Sampling requirements', text: 'f_s > 2·f_max (Nyquist),  f_max ≈ 1/(2πτ_fast);  practical: T_s ≤ τ_fast/10;  hold delay ≈ T_s/2 ⇒ ΔPM = −ω_c·T_s/2' },
      { kind: 'formula', label: 'Gauge R&R (ANOVA)', text: 'σ²_repeat = MS_E;  σ²_oper = (MS_O − MS_PO)/(p·r);  σ²_PO = (MS_PO − MS_E)/r;  σ²_part = (MS_P − MS_PO)/(o·r);  %GRR = 100·σ_GRR/σ_total' },
      { kind: 'formula', label: 'EWMA and CUSUM', text: 'zᵢ = λxᵢ + (1 − λ)zᵢ₋₁,  limits μ₀ ± L·σ·√[λ/(2 − λ)·(1 − (1 − λ)²ⁱ)];   C⁺ᵢ = max(0, C⁺ᵢ₋₁ + (xᵢ − μ₀)/σ − k),  alarm if C > h' },
    ],
    checks: [
      { name: 'Nyquist condition for the fastest phenomenon', pass: Ts <= TsNyq, detail: `Ts = ${F(Ts, 3)} s vs limit ${F(TsNyq, 3)} s` },
      { name: 'Ten samples per fastest time constant', pass: Ts <= inp.tauFast / 10, detail: `Ts = ${F(Ts, 3)} s vs ${F(inp.tauFast / 10, 3)} s` },
      { name: 'Phase margin ≥ 30° including sampling delay', pass: pmS >= 30, detail: `${F(pmS, 4)}° (continuous ${F(pm, 4)}°)` },
      { name: 'Simulated loop stable and settled', pass: settled, detail: settled ? `settled in ${F(tset, 4)} s, overshoot ${pct(Math.max(over, 0))}` : 'no settling within the simulated horizon' },
      { name: 'Measurement system adequate (%GRR ≤ 30 %, ndc ≥ 5)', pass: grr.ok && pGrr <= 0.3 && ndc >= 5, detail: grr.ok ? `%GRR = ${pct(pGrr)}, ndc = ${ndc}` : 'not evaluated' },
      { name: 'Process signal free of alarms', pass: !ewAlarm && !cuAlarm, detail: `EWMA: ${ewAlarm || 'none'}; CUSUM: ${cuAlarm || 'none'}` },
    ],
    warnings: w,
    recommendations: [
      { priority: Ts > TsNyq ? 'high' : Ts > TsRec ? 'medium' : 'low', text: Ts > TsRec ? `Shorten the sampling interval to ${F(TsRec, 3)} s or less and place an analogue anti-alias filter below ${F(0.4 / TsRec, 3)} Hz; at ${F(Ts, 3)} s the fastest process events are ${Ts > TsNyq ? 'aliased' : 'under-resolved'}.` : `Sampling at ${F(Ts, 3)} s resolves the process; the rate can be relaxed to ${F(TsRec, 3)} s to cut data volume by ${pct(1 - Ts / TsRec)} without loss of information.`, basis: 'Nyquist criterion and ten-samples-per-time-constant rule.' },
      { priority: pmS < 30 || !settled ? 'high' : 'low', text: pmS < 30 || !settled ? `Detune the loop (τc/θ = ${F(Math.max(inp.lam * 2, 2), 2)}) or reduce dead time by moving the sensor closer to the reaction zone: the phase margin of ${F(pmS, 3)}° is too small for robust temperature control.` : `The loop is robust (phase margin ${F(pmS, 3)}°); record Kc = ${F(Kc, 3)} and Tᵢ = ${F(Ti, 3)} s as the starting tuning and verify with a plant step test.`, basis: 'SIMC tuning, frequency-response margins and step simulation.' },
      { priority: grr.ok && pGrr > 0.3 ? 'high' : 'medium', text: !grr.ok ? 'Carry out a crossed gauge study (10 parts × 3 operators × 2–3 repeats) before using the measurement for control charts or capability.' : pGrr > 0.1 ? `Improve the measurement system before tightening the process: ${grr.ev > grr.av ? `repeatability dominates (σ = ${F(grr.ev, 3)}), so improve the instrument, sampling statistics or image-analysis protocol` : `reproducibility dominates (σ = ${F(grr.av, 3)}), so standardise the operator procedure or automate the evaluation`}.` : 'The measurement system is adequate for process control; repeat the study after any change of instrument, software or procedure.', basis: 'ANOVA variance components of the gauge study.' },
      { priority: 'medium', text: `Use the EWMA/CUSUM alarms as an early-warning layer of the digital twin: detecting a 1σ drift after about ${F(arl1, 3)} observations instead of at final inspection prevents off-specification batches and the precursor, solvent and energy they would consume.`, basis: 'Average run length of the tabular CUSUM for a 1σ shift.', sustainability: true },
    ],
    publish: { 'ctrl.settling_s': settled ? tset : nSt * dt, 'ctrl.kc': Kc, 'ctrl.ti_s': Ti, 'ctrl.ts_s': TsRec, 'ctrl.phase_margin_deg': pmS, 'qual.grr_pct': 100 * pGrr },
    numerics: { method: 'Exact zero-order-hold discretisation of the FOPDT process with sampled PI control; Brent search for the gain crossover; two-way crossed ANOVA; EWMA and tabular CUSUM recursions', dof: 2, steps: nSt, iterations: nSt, residual: Math.abs(1 - y), notes: 'Linear single-loop model without actuator saturation or noise; ARL from the Siegmund approximation.' },
  };
}

export default {
  id: 26,
  tagline: 'Sizes reactors and scales heat transfer and mixing, monitors quality statistically, analyses factorial experiments and quantifies yield, throughput and automation readiness.',
  inputs,
  presets: [
    { name: 'ZnO nanorods, hydrothermal batch, 10 t/yr', note: 'Scale-up of a 1 L, 90 °C hydrothermal recipe to production: batch time, vessel size, heat-up and micromixing versus nucleation.', values: { mode: 'reactor', rtype: 'batch', route: 'zno', T: 90, c0: 50, demand: 10, X: 0.9 } },
    { name: 'Ag nanoparticle colloid, continuous flow', note: 'Fast borohydride reduction at room temperature in plug flow: mixing-limited nucleation and numbering-up of millifluidic reactors.', values: { mode: 'reactor', rtype: 'pfr', route: 'ag', order: 1, k0: 2e7, Ea: 45, T: 25, T0: 20, c0: 1, X: 0.99, demand: 0.05, tInd: 0.02, dHr: -150, Vlab: 0.1, Nlab: 800, Qunit: 20, mu: 9e-4, isoYield: 95 } },
    { name: 'Lithographic metasurface line: yield and bottleneck', note: 'Nano-imprint metasurface dies of 1 cm²: defect-limited yield, rolled throughput yield and takt/bottleneck analysis of the patterning line.', values: { mode: 'yield', steps: STEPS_DEF, D0: 0.5, area: 1, critFrac: 0.4, ymodel: 'negbin', units: 7500, hours: 4000 } },
    { name: '3-D printed auxetic lattice: strut-diameter SPC', note: 'Strut diameter (µm) of printed re-entrant lattices, four struts per build: control chart and capability against 250 ± 12 µm.', values: { mode: 'spc', spcData: SPC_LAT, usl: 262, lsl: 238, target: 250 } },
    { name: 'ZnO nanorod length: factorial experiment', note: 'Replicated 2³ factorial with centre points on temperature, zinc concentration and growth time; finds the active effects and the best corner.', values: { mode: 'doe', doeData: DOE_DEF, goal: 'target', yTarget: 1500 } },
    { name: 'ZnO powder plant: batch line throughput', note: 'Five-step batch line (make-up, reactor, centrifuge, drying, packing) with campaign schedule; units are batches.', values: { mode: 'yield', steps: STEPS_ZNO, D0: 0.001, area: 1, critFrac: 0.01, units: 1500, hours: 6000, batchSize: 8, tBatch: 5, tClean: 1.5 } },
    { name: 'Printed lattice panels: additive line', note: 'Photopolymer lattice panels of 100 cm² with missing-strut defects; print step is the bottleneck.', values: { mode: 'yield', steps: STEPS_AM, D0: 0.004, area: 100, critFrac: 0.6, ymodel: 'murphy', units: 2400, hours: 5000 } },
  ],
  run(inp, ctx) {
    const r = ({ reactor: runReactor, spc: runSpc, doe: runDoe, yield: runYield, auto: runAuto }[inp.mode] || runReactor)(inp, ctx);
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    r.warnings.push('Reduced-order engineering estimates (ideal reactors, vessel-averaged correlations, normal-distribution statistics): confirm at pilot scale and with plant data before committing equipment or releasing product.');
    return r;
  },
};
