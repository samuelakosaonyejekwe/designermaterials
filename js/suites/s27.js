// Suite 27 — Safety, Standards, Environment and Device Integration.
// Occupational exposure (one-box and near-field/far-field), aerosol behaviour and lung deposition, hydrothermal autoclave
// pressure safety, flammability / dust / electrical screening, FMEA and fault-tree analysis, environmental fate and
// effluent neutralisation, and device-integration budgets (film stress, thermal stack, contacts, moisture ingress).
// Modelling of this kind supports a safety case; it never by itself establishes regulatory compliance or safety.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS, PATM = 101325;
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const pct = v => F(100 * v, 3) + ' %';
const csv = (head, rows) => head.join(',') + '\n' + rows.map(r => r.join(',')).join('\n');
function readRows(text, minCols, minRows, label) {
  const lines = String(text ?? '').split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#') && !l.startsWith('//'));
  const cells = lines.map(l => l.split(/[,;\t]/).map(s => s.trim())).filter(c => c.length >= minCols);
  if (cells.length && cells[0].slice(1).some(c => c === '' || !Number.isFinite(Number(c)))) cells.shift();
  if (cells.length < minRows) throw new Error(`${label}: at least ${minRows} data rows with ${minCols} columns are needed (found ${cells.length}). Paste CSV/TSV text with one header line.`);
  return cells;
}
const NOTE = 'Modelling alone never establishes regulatory compliance or safety: these screening estimates must be confirmed by measurement, formal risk assessment and the applicable national regulations and standards.';

// ---------- reference data ----------
// Recommended exposure limits (8-h TWA) that are widely cited; verify against current national regulation before use.
const LIMITS = {
  user: { label: 'User-entered limit', v: null },
  cnt: { label: 'Carbon nanotubes / nanofibres — NIOSH REL 1 µg/m³ (elemental carbon)', v: 0.001 },
  tio2_uf: { label: 'Ultrafine (nanoscale) TiO₂ — NIOSH REL 0.3 mg/m³', v: 0.3 },
  tio2_fine: { label: 'Fine TiO₂ — NIOSH REL 2.4 mg/m³', v: 2.4 },
};
const APF = [[10, 'half-mask air-purifying respirator with high-efficiency particulate filters'], [25, 'loose-fitting powered air-purifying respirator'], [50, 'full-facepiece air-purifying respirator'], [1000, 'tight-fitting full-facepiece powered or supplied-air respirator'], [10000, 'pressure-demand self-contained breathing apparatus']];
const BAND = [[1, 1, 1, 2], [1, 1, 2, 3], [2, 2, 3, 4], [3, 3, 4, 4]];   // rows: severity low → very high; columns: probability extremely unlikely → probable
const BANDCTRL = ['', 'general ventilation', 'fume hood or local exhaust ventilation', 'containment (glove box, enclosed transfer)', 'seek specialist advice before work proceeds'];
// Antoine constants (log10 P[mmHg] = A − B/(C + T[°C])), flammability limits (vol %), flash point (°C), density (kg/L)
const SOLV = {
  ethanol: { label: 'Ethanol', f: 'C2H6O', A: 8.20417, B: 1642.89, C: 230.300, lel: 3.3, uel: 19, fp: 13, bp: 78.4, rho: 0.789 },
  methanol: { label: 'Methanol', f: 'CH4O', A: 8.08097, B: 1582.271, C: 239.726, lel: 6.0, uel: 36, fp: 11, bp: 64.7, rho: 0.792 },
  isopropanol: { label: 'Isopropanol', f: 'C3H8O', A: 8.11778, B: 1580.92, C: 219.61, lel: 2.0, uel: 12.7, fp: 12, bp: 82.3, rho: 0.786 },
  acetone: { label: 'Acetone', f: 'C3H6O', A: 7.11714, B: 1210.595, C: 229.664, lel: 2.5, uel: 12.8, fp: -20, bp: 56.1, rho: 0.791 },
  toluene: { label: 'Toluene', f: 'C7H8', A: 6.95464, B: 1344.8, C: 219.482, lel: 1.1, uel: 7.1, fp: 4, bp: 110.6, rho: 0.867 },
  glycol: { label: 'Ethylene glycol', f: 'C2H6O2', A: 8.09083, B: 2088.936, C: 203.454, lel: 3.2, uel: 15.3, fp: 111, bp: 197.3, rho: 1.113 },
};
// powders: oxide formed and its standard enthalpy of formation (kJ/mol oxide); null = not combustible
const DUST = {
  ZnO: { label: 'ZnO nanopowder (fully oxidised)', ox: null }, TiO2: { label: 'TiO₂ nanopowder (fully oxidised)', ox: null }, Ag: { label: 'Ag nanopowder (noble metal)', ox: null },
  Al: { label: 'Al nanopowder', ox: 'Al2O3', dH: -1675.7 }, Ti: { label: 'Ti nanopowder', ox: 'TiO2', dH: -944.0 }, Zn: { label: 'Zn nanopowder', ox: 'ZnO', dH: -350.5 }, Mg: { label: 'Mg nanopowder', ox: 'MgO', dH: -601.6 }, Fe: { label: 'Fe nanopowder', ox: 'Fe2O3', dH: -824.2 }, C: { label: 'Carbon nanomaterial (carbon black, CNT, graphene)', ox: 'CO2', dH: -393.5 },
};
const ENVM = {
  ZnO: { label: 'ZnO', f: 'ZnO', ion: 'Zn2+', metal: 'Zn', logK: 11.2, rho: 5610, re: ['ZnO', 'H+'], pr: ['Zn2+', 'H2O'] },
  CuO: { label: 'CuO', f: 'CuO', ion: 'Cu2+', metal: 'Cu', logK: 7.65, rho: 6310, re: ['CuO', 'H+'], pr: ['Cu2+', 'H2O'] },
  Ag: { label: 'Ag (oxidative dissolution)', f: 'Ag', ion: 'Ag+', metal: 'Ag', logK: null, rho: 10490, re: ['Ag', 'O2', 'H+'], pr: ['Ag+', 'H2O'] },
  TiO2: { label: 'TiO₂ (practically insoluble)', f: 'TiO2', ion: null, metal: 'Ti', logK: null, rho: 4230, re: null, pr: null },
};
const WASTE = {
  naoh_hcl: { label: 'Alkaline mother liquor (NaOH) neutralised with HCl', re: ['NaOH', 'HCl'], pr: ['NaCl', 'H2O'] },
  naoh_h2so4: { label: 'Alkaline mother liquor (NaOH) neutralised with H₂SO₄', re: ['NaOH', 'H2SO4'], pr: ['Na2SO4', 'H2O'] },
  hno3_naoh: { label: 'Nitric-acid waste neutralised with NaOH', re: ['HNO3', 'NaOH'], pr: ['NaNO3', 'H2O'] },
  hcl_naoh: { label: 'Hydrochloric-acid waste neutralised with NaOH', re: ['HCl', 'NaOH'], pr: ['NaCl', 'H2O'] },
};
const FM = (rows) => csv(['failure_mode', 'severity', 'occurrence', 'detection'], rows);
const FMEA = {
  hydro: FM([['Autoclave overfilled: hydraulic over-pressure on heating', 10, 3, 4], ['Seal or liner failure: hot caustic leak', 8, 4, 4], ['Vessel opened before fully cooled', 9, 3, 3], ['Oven thermostat failure: over-temperature', 8, 2, 5], ['Wrong precursor concentration: off-specification nanorods', 4, 5, 3], ['Dry nanopowder released during transfer or milling', 7, 5, 6], ['Caustic splash during solution make-up', 6, 4, 2], ['Nitrate-rich effluent discharged untreated', 6, 3, 5]]),
  colloid: FM([['Hydrogen evolution from borohydride in a closed vessel', 8, 4, 4], ['Nitrogen-oxide exposure during aqua-regia cleaning', 8, 3, 3], ['Silver nitrate or gold salt skin contact', 4, 5, 2], ['Uncontrolled aggregation: batch lost', 5, 5, 3], ['Nanoparticle aerosol from sonication or spraying', 7, 4, 6], ['Effluent with dissolved silver to drain', 7, 3, 5], ['Contaminated glassware: irreproducible particle size', 4, 6, 5]]),
  litho: FM([['Solvent vapour ignition at the spin coater', 9, 2, 4], ['Hydrofluoric-acid etchant skin contact', 10, 2, 3], ['Charging during e-beam exposure: pattern distortion', 5, 5, 4], ['Particle contamination: pattern defects', 5, 7, 4], ['Lift-off redeposition: shorted meta-atoms', 4, 6, 3], ['Alkaline developer skin exposure', 9, 2, 4], ['Substrate breakage in handling', 3, 4, 2], ['Plasma-etch exhaust abatement failure', 8, 2, 5]]),
  am: FM([['Skin sensitisation by uncured photopolymer resin', 6, 5, 3], ['Isopropanol wash vapour ignition', 8, 3, 4], ['Missing or fused lattice struts', 6, 6, 4], ['Warpage during post-cure: unit cell distortion', 5, 5, 3], ['Ultraviolet exposure of eyes or skin', 6, 3, 3], ['Metal powder dust ignition (powder-bed processes)', 9, 2, 5], ['Lattice damage during support removal', 4, 5, 3], ['Uncured resin waste disposed as general waste', 5, 4, 5]]),
};
const EV_DEF = csv(['event', 'probability'], [['overfill', 0.005], ['relief_fail', 0.02], ['thermostat_fail', 0.002], ['trip_fail', 0.05], ['corrosion', 0.001], ['inspection_miss', 0.1], ['early_opening', 0.005], ['interlock_fail', 0.01], ['procedure_violation', 0.02]]);
const TREE_DEF = 'overfill*relief_fail + thermostat_fail*trip_fail + corrosion*inspection_miss + early_opening*(interlock_fail + procedure_violation)';
const STANDARDS = [
  ['ISO/TS 80004 series', 'Nanotechnologies — vocabulary', 'Use the defined terms (nano-object, nanoparticle, nanofibre, agglomerate…) in all records.'],
  ['ISO/TR 13121', 'Nanomaterial risk evaluation', 'Follow a documented process for identifying, evaluating and communicating nanomaterial risks.'],
  ['ISO 10993 series', 'Biological evaluation of medical devices', 'Required when the material contacts the body as part of a medical device; biological testing cannot be replaced by modelling.'],
  ['ISO 14040 / ISO 14044', 'Life-cycle assessment principles and requirements', 'Environmental claims need an assessment that follows these principles (see Suite 28).'],
  ['ISO 9001', 'Quality management systems', 'Controlled procedures, records and change control for production and testing.'],
  ['IEC 62607 series', 'Nanomanufacturing — key control characteristics', 'Measure and report the key control characteristics defined for the material class.'],
  ['ASTM E2456', 'Standard terminology relating to nanotechnology', 'Terminology reference for specifications and test reports.'],
  ['EU REACH (nanoform requirements)', 'Registration and characterisation of nanoforms of substances', 'Nanoforms placed on the EU market need substance-specific characterisation and registration data.'],
];
const QUAL = [
  ['Temperature cycling', 'Thermal-mismatch stress, film cracking, delamination', 'Film stress and energy-release rate'],
  ['Damp-heat exposure (elevated temperature and humidity)', 'Moisture ingress, corrosion, optical drift', 'Fickian ingress time through the encapsulant'],
  ['High-temperature storage or operating life', 'Diffusion, grain growth, contact degradation', 'Junction temperature'],
  ['Thermal shock', 'Brittle fracture of films and substrates', 'Film stress and wafer bow'],
  ['Mechanical shock and vibration', 'Fracture of nanorods, lattice struts and bonds', 'Structural margins (Suite 14 and 19)'],
  ['Adhesion test (tape, scratch or pull-off)', 'Interface toughness of films and coatings', 'Energy-release rate against interface toughness'],
  ['Electrical continuity and contact-resistance monitoring', 'Contact degradation, electromigration', 'Resistance budget'],
  ['Ultraviolet or weathering exposure', 'Photo-degradation of polymers and encapsulants', 'Not modelled here'],
  ['Particle-release (abrasion) test', 'Release of nano-objects during use and wear', 'Release fraction in the environmental-fate mode'],
];

// ---------- property functions ----------
// IAPWS-IF97 saturation-pressure equation (region 4), T in K → Pa, valid 273.15–647.096 K
const NW = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849, 0.65017534844798e3];
function psatW(T) { T = N.clamp(T, 273.15, 647.096); const th = T + NW[8] / (T - NW[9]), A = th * th + NW[0] * th + NW[1], B = NW[2] * th * th + NW[3] * th + NW[4], C = NW[5] * th * th + NW[6] * th + NW[7]; return 1e6 * Math.pow(2 * C / (-B + Math.sqrt(B * B - 4 * A * C)), 4); }
// IAPWS auxiliary equations for saturated liquid and vapour density of water (kg/m³)
const TC = 647.096, RHOC = 322;
const rhoL = T => { const t = 1 - N.clamp(T, 273.15, TC) / TC; return RHOC * (1 + 1.99274064 * t ** (1 / 3) + 1.09965342 * t ** (2 / 3) - 0.510839303 * t ** (5 / 3) - 1.75493479 * t ** (16 / 3) - 45.5170352 * t ** (43 / 3) - 6.74694450e5 * t ** (110 / 3)); };
const rhoV = T => { const t = 1 - N.clamp(T, 273.15, TC) / TC; return RHOC * Math.exp(-2.03150240 * t ** (2 / 6) - 2.68302940 * t ** (4 / 6) - 5.38626492 * t ** (8 / 6) - 17.2991605 * t ** 3 - 44.7586581 * t ** (37 / 6) - 63.9201063 * t ** (71 / 6)); };
const kappaW = Tc => 4.5e-10 * Math.exp(0.67 * ((Tc - 50) / 150) ** 2);   // approximate isothermal compressibility of liquid water, 1/Pa (±20 %)
const antoine = (s, Tc) => 133.322 * Math.pow(10, s.A - s.B / (s.C + Tc));
const muAir = T => 1.716e-5 * Math.pow(T / 273.15, 1.5) * (273.15 + 110.4) / (T + 110.4);
const mfpAir = T => 6.64e-8 * T / 293.15;
const slip = (d, T) => { const kn = 2 * mfpAir(T) / d; return 1 + kn * (1.257 + 0.4 * Math.exp(-1.1 / kn)); };
// simplified ICRP 66 deposition fits (Hinds, Aerosol Technology, 2nd ed.): d in µm, adult, light exercise, nose breathing
function lung(dum) {
  const d = N.clamp(dum, 0.001, 100), l = Math.log(d), IF = 1 - 0.5 * (1 - 1 / (1 + 0.00076 * Math.pow(d, 2.8)));
  const ha = IF * (1 / (1 + Math.exp(6.84 + 1.183 * l)) + 1 / (1 + Math.exp(0.924 - 1.885 * l)));
  const tb = (0.00352 / d) * (Math.exp(-0.234 * (l + 3.40) ** 2) + 63.9 * Math.exp(-0.819 * (l - 1.61) ** 2));
  const al = (0.0155 / d) * (Math.exp(-0.416 * (l + 2.84) ** 2) + 19.11 * Math.exp(-0.482 * (l - 1.362) ** 2));
  const tot = IF * (0.0587 + 0.911 / (1 + Math.exp(4.77 + 1.485 * l)) + 0.943 / (1 + Math.exp(0.508 - 2.58 * l)));
  return { IF, ha, tb, al, tot };
}
// Fuchs coagulation coefficient for equal spheres (m³/s), convention dN/dt = −½·K·N²
function coagK(d, rho, T) { const mu = muAir(T), D = N.K_B * T * slip(d, T) / (3 * Math.PI * mu * d), m = rho * Math.PI * d ** 3 / 6, c = Math.sqrt(8 * N.K_B * T / (Math.PI * m)), l = 8 * D / (Math.PI * c), g = ((d + l) ** 3 - Math.pow(d * d + l * l, 1.5)) / (3 * d * l) - d; return 8 * Math.PI * D * d / (d / (d + Math.SQRT2 * g) + 4 * Math.SQRT2 * D / (c * d)); }

const SH = m => ({ mode: Array.isArray(m) ? m : [m] });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'exposure', group: 'Model', options: [
    { value: 'exposure', label: 'Occupational exposure to airborne nanoparticles' }, { value: 'aerosol', label: 'Aerosol behaviour and lung deposition' }, { value: 'autoclave', label: 'Hydrothermal autoclave pressure safety' },
    { value: 'hazards', label: 'Flammable solvent, dust and electrical screening' }, { value: 'risk', label: 'Risk analysis: FMEA, fault tree, risk matrix' }, { value: 'envfate', label: 'Environmental fate, release and effluent' }, { value: 'device', label: 'Device integration and standards checklist' }],
    help: 'Selects the safety, environmental or integration assessment. All results are screening estimates that support, but never replace, measurement and formal assessment.' },
  // shared particle description
  { key: 'dp', label: 'Particle diameter', unit: 'm', value: 5e-8, min: 1e-9, max: 1e-5, log: true, bus: 'geom.diameter', group: 'Particles', showIf: SH(['exposure', 'aerosol', 'envfate']), help: 'Diameter of the airborne or suspended particles (50 nm = 5×10⁻⁸ m). For agglomerates use the agglomerate size; for nanorods the equivalent mobility diameter.' },
  { key: 'rhoP', label: 'Particle density', unit: 'kg/m³', value: 5610, min: 500, max: 22000, bus: 'prop.density', group: 'Particles', showIf: SH(['exposure', 'aerosol']), help: 'Material density of the particles (ZnO 5610, TiO₂ 4230, Ag 10490, Au 19300, carbon about 2000 kg/m³). Used to convert mass to number concentration and for settling.' },
  // exposure
  { key: 'G', label: 'Emission rate', unit: 'mg/min', value: 0.5, min: 1e-4, max: 1000, log: true, group: 'Source', showIf: SH('exposure'), help: 'Mass of nanomaterial released to air per minute during the task (weighing, pouring, harvesting, milling). Measure with a dustiness test or take from task-based data.' },
  { key: 'tTask', label: 'Task duration', unit: 'min', value: 60, min: 1, max: 480, group: 'Source', showIf: SH('exposure'), help: 'Time during which the source emits and the worker stays next to it. Afterwards the worker is assumed to remain in the room for the rest of the 8-h shift.' },
  { key: 'lev', label: 'Local exhaust capture efficiency', unit: '%', value: 70, min: 0, max: 99.99, group: 'Controls', showIf: SH('exposure'), help: 'Fraction of the emission captured at source by a fume hood, enclosure or extraction arm before it reaches the room air.' },
  { key: 'Vroom', label: 'Room volume', unit: 'm³', value: 60, min: 5, max: 5000, log: true, group: 'Room', showIf: SH('exposure'), help: 'Air volume of the room in which the task is carried out.' },
  { key: 'ach', label: 'General ventilation', unit: '1/h', value: 6, min: 0.1, max: 60, group: 'Room', showIf: SH('exposure'), help: 'Fresh-air changes per hour of the room. Laboratories typically have 6–12, offices 1–3.' },
  { key: 'recirc', label: 'Recirculated filtered air', unit: '1/h', value: 2, min: 0, max: 60, group: 'Room', showIf: SH('exposure'), help: 'Room air changes per hour passed through a recirculating filter unit.' },
  { key: 'filtEff', label: 'Filter efficiency', unit: '%', value: 99.95, min: 0, max: 100, group: 'Room', showIf: SH('exposure'), help: 'Single-pass particle removal efficiency of the recirculation filter at the particle size of interest.' },
  { key: 'kdep', label: 'Deposition loss rate', unit: '1/h', value: 0.2, min: 0, max: 10, group: 'Room', showIf: SH('exposure'), help: 'First-order loss of particles to room surfaces by diffusion and settling. It is size dependent; 0.1–0.5 per hour is typical for 20–200 nm particles.' },
  { key: 'Vnf', label: 'Near-field volume', unit: 'm³', value: 1, min: 0.1, max: 20, group: 'Room', showIf: SH('exposure'), help: 'Volume of the zone around the source that contains the worker\'s breathing zone (about 1 m³ for bench work).' },
  { key: 'beta', label: 'Near-field/far-field air exchange', unit: 'm³/min', value: 5, min: 0.2, max: 30, group: 'Room', showIf: SH('exposure'), help: 'Air flow exchanged between the near field and the rest of the room; 3–10 m³/min corresponds to typical random room air speeds.' },
  { key: 'limitSrc', label: 'Exposure limit', type: 'select', value: 'user', group: 'Limit', showIf: SH('exposure'), options: Object.entries(LIMITS).map(([value, l]) => ({ value, label: l.label })), help: 'Limit used for comparison. Only a few widely cited recommended limits are built in; verify every limit against current national regulation and the supplier\'s safety data sheet.' },
  { key: 'oel', label: 'User exposure limit (8-h TWA)', unit: 'mg/m³', value: 0.3, min: 1e-5, max: 100, log: true, group: 'Limit', showIf: SH('exposure'), help: 'Occupational exposure limit or in-house benchmark for the material, applied when "user-entered limit" is selected.' },
  { key: 'sevBand', label: 'Hazard severity band', type: 'select', value: 'high', group: 'Control banding', showIf: SH('exposure'), options: [{ value: 'low', label: 'Low (well-characterised, low toxicity)' }, { value: 'medium', label: 'Medium' }, { value: 'high', label: 'High (or hazard unknown)' }, { value: 'vhigh', label: 'Very high (carcinogenic, mutagenic, fibre-like)' }], help: 'Severity band from the toxicological knowledge of the material and its parent substance. Unknown hazard is treated as high.' },
  { key: 'form', label: 'Physical form handled', type: 'select', value: 'dry', group: 'Control banding', showIf: SH('exposure'), options: [{ value: 'dry', label: 'Dry powder (dusty)' }, { value: 'damp', label: 'Damp powder, paste or granulate' }, { value: 'liquid', label: 'Liquid suspension' }, { value: 'solid', label: 'Embedded in a solid matrix' }], help: 'Physical form determines dustiness and therefore the probability of exposure.' },
  { key: 'mHandle', label: 'Amount handled per task', unit: 'mg', value: 5000, min: 0.1, max: 1e7, log: true, group: 'Control banding', showIf: SH('exposure'), help: 'Mass of nanomaterial handled in one task.' },
  { key: 'nWork', label: 'Workers with similar exposure', value: 2, min: 1, max: 100, step: 1, group: 'Control banding', showIf: SH('exposure'), help: 'Number of people carrying out the task.' },
  { key: 'freq', label: 'Tasks per week', value: 5, min: 0.1, max: 40, group: 'Control banding', showIf: SH('exposure'), help: 'How often the task is carried out.' },
  // aerosol
  { key: 'N0', label: 'Number concentration', unit: '1/cm³', value: 1e5, min: 100, max: 1e9, log: true, group: 'Aerosol', showIf: SH('aerosol'), help: 'Initial particle number concentration of the aerosol; background indoor air holds roughly 10³–10⁴ particles per cm³.' },
  { key: 'Tair', label: 'Air temperature', unit: '°C', value: 20, min: 0, max: 60, group: 'Aerosol', showIf: SH('aerosol'), help: 'Temperature of the workplace air; it changes viscosity and mean free path slightly.' },
  { key: 'Hroom', label: 'Settling height', unit: 'm', value: 2.5, min: 0.1, max: 20, group: 'Aerosol', showIf: SH('aerosol'), help: 'Height through which particles must settle to be removed (room height or breathing-zone height).' },
  { key: 'cAir', label: 'Mass concentration breathed', unit: 'mg/m³', value: 0.05, min: 1e-6, max: 100, log: true, bus: 'risk.exposure_mg_m3', group: 'Aerosol', showIf: SH('aerosol'), help: 'Time-weighted airborne mass concentration in the breathing zone, for the deposited-dose estimate.' },
  { key: 'vent', label: 'Breathing rate', unit: 'L/min', value: 25, min: 5, max: 100, group: 'Aerosol', showIf: SH('aerosol'), help: 'Minute ventilation of the worker: about 9 L/min at rest, 25 L/min for light work, 50 L/min for heavy work.' },
  // autoclave
  { key: 'Tac', label: 'Autoclave temperature', unit: '°C', value: 180, min: 100, max: 350, bus: 'proc.T_C', group: 'Autoclave', showIf: SH('autoclave'), help: 'Set temperature of the hydrothermal treatment. The saturation pressure of water rises steeply: 1.0 MPa at 180 °C, 1.55 MPa at 200 °C, 3.98 MPa at 250 °C.' },
  { key: 'fillPct', label: 'Fill fraction at 20 °C', unit: '%', value: 70, min: 5, max: 98, group: 'Autoclave', showIf: SH('autoclave'), help: 'Liquid volume as a percentage of the liner volume when loaded. Water expands on heating; an over-filled vessel becomes liquid-full and the pressure then rises almost without limit.' },
  { key: 'Vliner', label: 'Liner volume', unit: 'mL', value: 100, min: 5, max: 50000, log: true, group: 'Autoclave', showIf: SH('autoclave'), help: 'Internal volume of the liner or vessel.' },
  { key: 'IDmm', label: 'Shell inner diameter', unit: 'mm', value: 50, min: 5, max: 1000, group: 'Autoclave', showIf: SH('autoclave'), help: 'Inner diameter of the pressure-bearing steel shell (not of the polymer liner).' },
  { key: 'wall', label: 'Shell wall thickness', unit: 'mm', value: 8, min: 0.5, max: 200, group: 'Autoclave', showIf: SH('autoclave'), help: 'Wall thickness of the steel shell.' },
  { key: 'yieldT', label: 'Yield strength at temperature', unit: 'MPa', value: 130, min: 20, max: 1500, group: 'Autoclave', showIf: SH('autoclave'), help: 'Yield strength of the shell material at the operating temperature from the material certificate. Austenitic stainless steels lose roughly a quarter of their room-temperature yield strength by 200 °C.' },
  { key: 'designSF', label: 'Required safety factor on yield', value: 1.5, min: 1, max: 5, step: 0.1, group: 'Autoclave', showIf: SH('autoclave'), help: 'Minimum ratio of yield strength to equivalent stress demanded by the design basis.' },
  // hazards
  { key: 'solvent', label: 'Solvent', type: 'select', value: 'ethanol', group: 'Flammable solvent', showIf: SH('hazards'), options: Object.entries(SOLV).map(([value, s]) => ({ value, label: s.label })), help: 'Flammable solvent used in solvothermal synthesis, washing or resist processing. Flammability limits and flash points are handbook values to be verified against the safety data sheet.' },
  { key: 'Tsolv', label: 'Liquid temperature', unit: '°C', value: 25, min: -20, max: 120, group: 'Flammable solvent', showIf: SH('hazards'), help: 'Temperature of the exposed solvent surface. Above the flash point the vapour above the liquid is ignitable.' },
  { key: 'Vsp', label: 'Solvent volume evaporated', unit: 'mL', value: 100, min: 0.1, max: 1e5, log: true, group: 'Flammable solvent', showIf: SH('hazards'), help: 'Volume of solvent that can evaporate into the enclosure (spill, wash bath, vessel venting).' },
  { key: 'Venc', label: 'Enclosure volume', unit: 'm³', value: 0.5, min: 0.005, max: 2000, log: true, group: 'Enclosure', showIf: SH('hazards'), help: 'Free air volume of the oven, glove box, cabinet or room into which vapour, dust or gas is released.' },
  { key: 'achEnc', label: 'Enclosure ventilation', unit: '1/h', value: 10, min: 0, max: 600, group: 'Enclosure', showIf: SH('hazards'), help: 'Air changes per hour of the enclosure.' },
  { key: 'tEvap', label: 'Evaporation time', unit: 'min', value: 30, min: 0.1, max: 1440, log: true, group: 'Flammable solvent', showIf: SH('hazards'), help: 'Time over which the solvent volume evaporates; shorter times need more ventilation.' },
  { key: 'dustMat', label: 'Powder', type: 'select', value: 'ZnO', group: 'Dust', showIf: SH('hazards'), options: Object.entries(DUST).map(([value, s]) => ({ value, label: s.label })), help: 'Powder handled. Fully oxidised ceramics do not burn; metal and carbon nanopowders can form explosible dust clouds and may ignite at very low energies.' },
  { key: 'mDust', label: 'Powder mass that can be dispersed', unit: 'g', value: 50, min: 0.01, max: 1e6, log: true, group: 'Dust', showIf: SH('hazards'), help: 'Mass of powder that could be raised into a cloud inside the enclosure in an upset.' },
  { key: 'mec', label: 'Minimum explosible concentration', unit: 'g/m³', value: 60, min: 5, max: 1000, log: true, group: 'Dust', showIf: SH('hazards'), help: 'Lowest dust concentration that propagates a flame, from a standard dust-explosion test of the actual powder. Many combustible dusts lie in the range 30–125 g/m³; the default is only a placeholder for test data.' },
  { key: 'Uapp', label: 'Applied voltage', unit: 'V', value: 30, min: 0.1, max: 10000, log: true, group: 'Electrical', showIf: SH('hazards'), help: 'Voltage across the electrodes of a field-assisted reactor or electrodeposition cell.' },
  { key: 'acdc', label: 'Supply type', type: 'select', value: 'dc', group: 'Electrical', showIf: SH('hazards'), options: [{ value: 'dc', label: 'Direct current' }, { value: 'ac', label: 'Alternating current (50/60 Hz)' }], help: 'Supply type; conventional touch-voltage limits differ for AC and ripple-free DC, and only DC produces net electrolysis gas.' },
  { key: 'gap', label: 'Electrode gap', unit: 'mm', value: 10, min: 0.01, max: 500, log: true, group: 'Electrical', showIf: SH('hazards'), help: 'Distance between the electrodes.' },
  { key: 'Ael', label: 'Electrode area', unit: 'cm²', value: 4, min: 0.01, max: 1e4, log: true, group: 'Electrical', showIf: SH('hazards'), help: 'Wetted electrode area carrying the current.' },
  { key: 'sigEl', label: 'Electrolyte conductivity', unit: 'S/m', value: 1, min: 1e-5, max: 100, log: true, group: 'Electrical', showIf: SH('hazards'), help: 'Conductivity of the liquid between the electrodes (deionised water about 10⁻⁴, 0.05 mol/L salt solution about 0.5–1 S/m).' },
  { key: 'Vliq', label: 'Liquid volume', unit: 'mL', value: 100, min: 0.1, max: 1e6, log: true, group: 'Electrical', showIf: SH('hazards'), help: 'Volume of liquid heated by the current.' },
  // risk
  { key: 'fmeaSet', label: 'Failure-mode list', type: 'select', value: 'table', group: 'FMEA', showIf: SH('risk'), options: [{ value: 'table', label: 'Use the editable table below' }, { value: 'hydro', label: 'Built-in: hydrothermal oxide synthesis' }, { value: 'colloid', label: 'Built-in: noble-metal colloid synthesis' }, { value: 'litho', label: 'Built-in: lithographic metasurface fabrication' }, { value: 'am', label: 'Built-in: additive manufacturing of lattices' }], help: 'Either analyse the editable table or load a built-in illustrative list for the process type. Built-in ratings are starting points for a team review, not assessed values.' },
  { key: 'fmea', label: 'Failure modes (failure_mode, severity, occurrence, detection)', type: 'table', columns: ['failure_mode', 'severity', 'occurrence', 'detection'], value: FMEA.hydro, group: 'FMEA', showIf: SH('risk'), help: 'One row per failure mode with severity, occurrence and detection ratings from 1 (best) to 10 (worst), agreed by the assessment team.' },
  { key: 'rpnAct', label: 'RPN action threshold', value: 100, min: 10, max: 1000, group: 'FMEA', showIf: SH('risk'), help: 'Risk priority number at or above which corrective action is required. Any failure mode with severity 9 or 10 requires action regardless of RPN.' },
  { key: 'events', label: 'Basic events (event, probability per demand)', type: 'table', columns: ['event', 'probability'], value: EV_DEF, group: 'Fault tree', showIf: SH('risk'), help: 'Basic events of the fault tree with their probability of occurring on one demand (one batch or operation). Names must match those used in the tree expression.' },
  { key: 'tree', label: 'Top-event logic (+ = OR, * = AND)', type: 'text', value: TREE_DEF, group: 'Fault tree', showIf: SH('risk'), help: 'Boolean expression of the top event in terms of basic events, with + for OR, * for AND and parentheses for grouping.' },
  { key: 'nDem', label: 'Demands per year', unit: '1/yr', value: 50, min: 1, max: 1e6, log: true, group: 'Fault tree', showIf: SH('risk'), help: 'Number of batches or operations per year to which the per-demand probabilities apply.' },
  { key: 'consClass', label: 'Consequence class of the top event', value: 4, min: 1, max: 5, step: 1, group: 'Fault tree', showIf: SH('risk'), help: 'Severity category of the top event on a 1–5 scale (1 negligible, 3 lost-time injury or significant release, 5 fatality or major environmental damage).' },
  // environmental fate
  { key: 'nm', label: 'Nanomaterial', type: 'select', value: 'ZnO', group: 'Release', showIf: SH('envfate'), options: Object.entries(ENVM).map(([value, s]) => ({ value, label: s.label })), help: 'Material released. Dissolving oxides act mainly through their ions; persistent particles through the particulate form.' },
  { key: 'mProd', label: 'Product placed on the market', unit: 't/yr', value: 100, min: 0.001, max: 1e6, log: true, group: 'Release', showIf: SH('envfate'), help: 'Annual mass of the nano-enabled product (coating, composite, sunscreen, textile) used in the region considered.' },
  { key: 'wNP', label: 'Nanomaterial content', unit: 'wt %', value: 2, min: 0.001, max: 100, log: true, group: 'Release', showIf: SH('envfate'), help: 'Mass fraction of nanomaterial in the product.' },
  { key: 'relFrac', label: 'Release during use and wear', unit: '%/yr', value: 5, min: 0.001, max: 100, log: true, group: 'Release', showIf: SH('envfate'), help: 'Fraction of the contained nanomaterial released to wastewater per year by abrasion, washing or weathering.' },
  { key: 'wwtp', label: 'Wastewater treatment removal', unit: '%', value: 90, min: 0, max: 99.99, group: 'Release', showIf: SH('envfate'), help: 'Fraction of the released nanomaterial retained in sewage sludge; the remainder reaches surface water.' },
  { key: 'Qriver', label: 'Receiving water flow', unit: 'm³/s', value: 20, min: 0.01, max: 1e5, log: true, group: 'Water body', showIf: SH('envfate'), help: 'Mean flow of the river that dilutes the treated effluent.' },
  { key: 'pnec', label: 'Predicted no-effect concentration', unit: 'µg/L', value: 1, min: 1e-4, max: 1e5, log: true, group: 'Water body', showIf: SH('envfate'), help: 'No-effect concentration for the aquatic compartment from ecotoxicological data for this material; it must be supplied by the user from a qualified source.' },
  { key: 'pHw', label: 'Water pH', value: 7.8, min: 4, max: 10, step: 0.1, group: 'Water body', showIf: SH('envfate'), help: 'pH of the receiving water; it controls the solubility of oxides such as ZnO and CuO.' },
  { key: 'csUser', label: 'Solubility for non-oxide routes', unit: 'mg/L', value: 0.05, min: 1e-6, max: 1e4, log: true, group: 'Water body', showIf: SH('envfate'), help: 'Dissolved concentration limit used for materials without a built-in solubility constant (Ag, TiO₂).' },
  { key: 'kDiss', label: 'Surface dissolution rate constant', unit: 'mol/m²/s', value: 1e-7, min: 1e-12, max: 1e-3, log: true, group: 'Water body', showIf: SH('envfate'), help: 'Surface-reaction-controlled dissolution flux far from saturation, from dissolution experiments on the actual particles.' },
  { key: 'alphaAtt', label: 'Attachment efficiency α', value: 0.1, min: 1e-4, max: 1, log: true, group: 'Water body', showIf: SH('envfate'), help: 'Probability that a collision leads to attachment: near 1 in hard or saline water, 0.001–0.1 when natural organic matter stabilises the particles.' },
  { key: 'spm', label: 'Suspended particulate matter', unit: 'mg/L', value: 15, min: 0.1, max: 1000, log: true, group: 'Water body', showIf: SH('envfate'), help: 'Concentration of natural suspended solids to which nanoparticles attach (heteroaggregation).' },
  { key: 'depth', label: 'Water depth', unit: 'm', value: 2, min: 0.1, max: 100, group: 'Water body', showIf: SH('envfate'), help: 'Depth of the water column for the sedimentation time.' },
  { key: 'wasteType', label: 'Process effluent and neutraliser', type: 'select', value: 'naoh_hcl', group: 'Effluent', showIf: SH('envfate'), options: Object.entries(WASTE).map(([value, s]) => ({ value, label: s.label })), help: 'Acidic or alkaline process effluent and the reagent used to neutralise it before discharge.' },
  { key: 'wasteVol', label: 'Effluent volume', unit: 'L', value: 1000, min: 0.1, max: 1e7, log: true, group: 'Effluent', showIf: SH('envfate'), help: 'Volume of effluent to be neutralised (per batch or per day).' },
  { key: 'wasteConc', label: 'Effluent acid/base concentration', unit: 'mol/L', value: 0.1, min: 1e-5, max: 10, log: true, group: 'Effluent', showIf: SH('envfate'), help: 'Concentration of the acid or base in the effluent, from titration.' },
  // device
  { key: 'Ef', label: 'Film Young\'s modulus', unit: 'Pa', value: 1.4e11, min: 1e7, max: 1.2e12, log: true, bus: 'prop.E', group: 'Film on substrate', showIf: SH('device'), help: 'Elastic modulus of the functional film or nanostructured layer (ZnO about 140 GPa, TiO₂ about 230 GPa, Au 79 GPa, polymers 1–5 GPa).' },
  { key: 'nuf', label: 'Film Poisson\'s ratio', value: 0.3, min: -0.9, max: 0.49, step: 0.01, bus: 'prop.nu', group: 'Film on substrate', showIf: SH('device'), help: 'Poisson\'s ratio of the film; auxetic layers have negative values, which lowers the biaxial modulus E/(1 − ν).' },
  { key: 'ctef', label: 'Film thermal expansion', unit: '1/K', value: 4.7e-6, min: -5e-5, max: 3e-4, bus: 'prop.cte', group: 'Film on substrate', showIf: SH('device'), help: 'Linear thermal-expansion coefficient of the film; negative values (negative-thermal-expansion materials) are allowed.' },
  { key: 'tf', label: 'Film thickness', unit: 'm', value: 5e-7, min: 1e-9, max: 1e-3, log: true, bus: 'geom.thickness', group: 'Film on substrate', showIf: SH('device'), help: 'Thickness of the film or coating (500 nm = 5×10⁻⁷ m).' },
  { key: 'Es', label: 'Substrate Young\'s modulus', unit: 'Pa', value: 1.3e11, min: 1e8, max: 1e12, log: true, group: 'Film on substrate', showIf: SH('device'), help: 'Elastic modulus of the substrate (silicon ⟨100⟩ 130 GPa, fused silica 72 GPa, sapphire 400 GPa, polyimide 2.5 GPa).' },
  { key: 'nus', label: 'Substrate Poisson\'s ratio', value: 0.28, min: 0, max: 0.49, step: 0.01, group: 'Film on substrate', showIf: SH('device'), help: 'Poisson\'s ratio of the substrate.' },
  { key: 'ctes', label: 'Substrate thermal expansion', unit: '1/K', value: 2.6e-6, min: -1e-5, max: 3e-4, group: 'Film on substrate', showIf: SH('device'), help: 'Linear thermal-expansion coefficient of the substrate (silicon 2.6×10⁻⁶, fused silica 0.5×10⁻⁶, sapphire about 5×10⁻⁶ K⁻¹).' },
  { key: 'ts', label: 'Substrate thickness', unit: 'm', value: 5.25e-4, min: 1e-5, max: 0.02, log: true, group: 'Film on substrate', showIf: SH('device'), help: 'Thickness of the substrate (a 100 mm silicon wafer is 525 µm).' },
  { key: 'Dw', label: 'Substrate diameter', unit: 'm', value: 0.1, min: 1e-3, max: 0.5, group: 'Film on substrate', showIf: SH('device'), help: 'Diameter (or longest span) over which the bow is evaluated.' },
  { key: 'Tdep', label: 'Deposition (stress-free) temperature', unit: '°C', value: 300, min: -50, max: 1200, group: 'Film on substrate', showIf: SH('device'), help: 'Temperature at which the film was formed and is assumed free of thermal stress.' },
  { key: 'Tuse', label: 'Service temperature', unit: '°C', value: 25, min: -200, max: 600, group: 'Film on substrate', showIf: SH('device'), help: 'Temperature at which the stress is evaluated.' },
  { key: 'Gc', label: 'Film or interface toughness', unit: 'J/m²', value: 5, min: 0.05, max: 1000, log: true, group: 'Film on substrate', showIf: SH('device'), help: 'Critical energy-release rate for film cracking (tension) or delamination (compression); 1–10 J/m² is typical for ceramic films and weak interfaces.' },
  { key: 'Pdev', label: 'Device power dissipation', unit: 'W', value: 2, min: 1e-4, max: 5000, log: true, bus: 'proc.power', group: 'Thermal stack', showIf: SH('device'), help: 'Heat generated in the active layer that must flow through the stack to ambient.' },
  { key: 'Adev', label: 'Heat-flow area', unit: 'cm²', value: 1, min: 1e-4, max: 1e4, log: true, group: 'Thermal stack', showIf: SH('device'), help: 'Area through which heat flows (one-dimensional conduction, no heat spreading).' },
  { key: 'kf', label: 'Film thermal conductivity', unit: 'W/m/K', value: 30, min: 0.01, max: 2500, log: true, bus: 'prop.k_thermal', group: 'Thermal stack', showIf: SH('device'), help: 'Through-thickness thermal conductivity of the film; nanostructured and porous layers are far below bulk values.' },
  { key: 'ks', label: 'Substrate thermal conductivity', unit: 'W/m/K', value: 150, min: 0.05, max: 2500, log: true, group: 'Thermal stack', showIf: SH('device'), help: 'Thermal conductivity of the substrate (silicon about 150, fused silica 1.4, sapphire 35, polyimide 0.2 W/m/K).' },
  { key: 'tTim', label: 'Interface material thickness', unit: 'µm', value: 50, min: 0, max: 2000, group: 'Thermal stack', showIf: SH('device'), help: 'Bond-line thickness of the thermal interface material or adhesive.' },
  { key: 'kTim', label: 'Interface material conductivity', unit: 'W/m/K', value: 3, min: 0.05, max: 400, log: true, group: 'Thermal stack', showIf: SH('device'), help: 'Thermal conductivity of the interface material (greases and filled adhesives 1–5, solders 20–60 W/m/K).' },
  { key: 'Rsa', label: 'Heat-sink-to-ambient resistance', unit: 'K/W', value: 5, min: 0.01, max: 500, log: true, group: 'Thermal stack', showIf: SH('device'), help: 'Thermal resistance from the back of the stack to ambient air (heat sink, package or free convection).' },
  { key: 'Tamb', label: 'Ambient temperature', unit: '°C', value: 25, min: -55, max: 150, group: 'Thermal stack', showIf: SH('device'), help: 'Temperature of the surroundings.' },
  { key: 'Tjmax', label: 'Maximum junction temperature', unit: '°C', value: 125, min: 40, max: 600, group: 'Thermal stack', showIf: SH('device'), help: 'Highest allowed temperature of the active layer.' },
  { key: 'sigf', label: 'Film electrical conductivity', unit: 'S/m', value: 1e5, min: 1e-6, max: 7e7, log: true, bus: 'prop.sigma_elec', group: 'Electrical budget', showIf: SH('device'), help: 'Conductivity of the conducting film (doped ZnO or ITO 10⁵–10⁶, Au 4.1×10⁷, Ag 6.3×10⁷ S/m).' },
  { key: 'squares', label: 'Trace length/width (squares)', value: 10, min: 0.01, max: 1e5, log: true, group: 'Electrical budget', showIf: SH('device'), help: 'Number of squares of the conducting path: length divided by width.' },
  { key: 'rhoc', label: 'Specific contact resistivity', unit: 'Ω·cm²', value: 1e-5, min: 1e-9, max: 10, log: true, group: 'Electrical budget', showIf: SH('device'), help: 'Area-specific resistance of the metal–film contact, from transfer-length measurements.' },
  { key: 'Ac', label: 'Contact area', unit: 'mm²', value: 0.01, min: 1e-6, max: 1e4, log: true, group: 'Electrical budget', showIf: SH('device'), help: 'Area of one electrical contact (two contacts are counted).' },
  { key: 'Idev', label: 'Operating current', unit: 'A', value: 0.01, min: 1e-9, max: 1000, log: true, group: 'Electrical budget', showIf: SH('device'), help: 'Current carried by the trace and contacts.' },
  { key: 'Rmax', label: 'Resistance budget', unit: 'Ω', value: 500, min: 1e-4, max: 1e9, log: true, group: 'Electrical budget', showIf: SH('device'), help: 'Largest series resistance the circuit tolerates.' },
  { key: 'Lenc', label: 'Encapsulant thickness', unit: 'µm', value: 1000, min: 0.01, max: 1e4, log: true, group: 'Encapsulation', showIf: SH('device'), help: 'Thickness of the moisture-barrier layer between ambient and the sensitive surface.' },
  { key: 'Dm', label: 'Moisture diffusivity', unit: 'm²/s', value: 1e-13, min: 1e-18, max: 1e-9, log: true, bus: 'prop.D_diff', group: 'Encapsulation', showIf: SH('device'), help: 'Diffusion coefficient of water in the encapsulant at service temperature (epoxies and silicones 10⁻¹³–10⁻¹¹, inorganic barriers many orders lower).' },
  { key: 'cThr', label: 'Critical moisture level at the device', unit: '% of saturation', value: 50, min: 1, max: 99, group: 'Encapsulation', showIf: SH('device'), help: 'Moisture concentration at the protected surface, relative to saturation, at which degradation starts.' },
];

// ---------- mode: occupational exposure ----------
function runExposure(inp) {
  const w = [], Gp = inp.G * (1 - inp.lev / 100), V = inp.Vroom, Vn = Math.min(inp.Vnf, 0.5 * V), Vf = V - Vn, lam = (inp.ach + inp.recirc * inp.filtEff / 100 + inp.kdep) / 60, Qe = lam * V, b = inp.beta;
  if (Vn < inp.Vnf) w.push('Near-field volume was limited to half the room volume.');
  const dt = 0.5, nS = 960, tT = N.clamp(inp.tTask, dt, 480), Cn = [0], Cf = [0], Cb = [0], ts = [0]; let cn = 0, cf = 0, emitted = 0, removed = 0, twa2 = 0, twa1 = 0, peak = 0;
  const a11 = -b / Vn, a12 = b / Vn, a21 = b / Vf, a22 = -(b + Qe) / Vf, h = dt / 2, m11 = 1 - h * a11, m12 = -h * a12, m21 = -h * a21, m22 = 1 - h * a22, det = m11 * m22 - m12 * m21;
  const one = t => t <= tT ? Gp / Qe * (1 - Math.exp(-lam * t)) : Gp / Qe * (1 - Math.exp(-lam * tT)) * Math.exp(-lam * (t - tT));
  for (let i = 0; i < nS; i++) {
    const t0 = i * dt, on = Math.max(0, Math.min(t0 + dt, tT) - t0) / dt, g = Gp * on, r1 = (1 + h * a11) * cn + h * a12 * cf + dt * g / Vn, r2 = h * a21 * cn + (1 + h * a22) * cf;
    const n1 = (r1 * m22 - m12 * r2) / det, f1 = (m11 * r2 - m21 * r1) / det; emitted += g * dt; removed += Qe * 0.5 * (cf + f1) * dt;
    const w0 = t0 < tT ? cn : cf, w1 = t0 + dt <= tT ? n1 : f1; twa2 += 0.5 * (w0 + w1) * dt / 480; twa1 += 0.5 * (one(t0) + one(t0 + dt)) * dt / 480; cn = n1; cf = f1; peak = Math.max(peak, cn);
    if ((i + 1) % 4 === 0) { ts.push(t0 + dt); Cn.push(cn); Cf.push(cf); Cb.push(one(t0 + dt)); }
  }
  const inv = Vn * cn + Vf * cf, mbErr = Math.abs(emitted - removed - inv) / Math.max(emitted, 1e-300), ssF = Gp / Qe, ssN = ssF + Gp / b, mp = inp.rhoP * Math.PI * inp.dp ** 3 / 6, num = c => c * 1e-6 / mp / 1e6;   // mg/m³ → particles/cm³
  const L = LIMITS[inp.limitSrc] || LIMITS.user, oel = L.v ?? inp.oel, hq = twa2 / oel, apf = APF.find(a => a[0] >= hq), needLev = hq > 0.5 ? N.clamp(1 - (1 - inp.lev / 100) * 0.5 / hq, 0, 1) : inp.lev / 100, tClear = Math.log(100) / lam;
  // control banding (probability score patterned on the CB Nanotool factors)
  const sev = { low: 0, medium: 1, high: 2, vhigh: 3 }[inp.sevBand] ?? 2, sc = ({ dry: 30, damp: 15, liquid: 7.5, solid: 0 }[inp.form] ?? 30) + (inp.mHandle > 100 ? 25 : inp.mHandle >= 10 ? 12.5 : 6.25) + (inp.nWork > 15 ? 15 : inp.nWork >= 11 ? 10 : inp.nWork >= 6 ? 5 : 0) + (inp.freq >= 5 ? 15 : inp.freq >= 1 ? 10 : 5) + (tT > 240 ? 15 : tT >= 60 ? 10 : tT >= 30 ? 5 : 0);
  const pb = sc <= 25 ? 0 : sc <= 50 ? 1 : sc <= 75 ? 2 : 3, rl = BAND[sev][pb], pName = ['extremely unlikely', 'less likely', 'likely', 'probable'][pb], present = inp.lev >= 99 ? 3 : inp.lev >= 50 ? 2 : 1;
  if (hq > 1) w.push(`The predicted 8-h TWA (${F(twa2, 3)} mg/m³) exceeds the limit of ${F(oel, 3)} mg/m³ by a factor ${F(hq, 3)}: engineering controls must be improved; respiratory protection is only an interim measure.`);
  if (rl > present) w.push(`Control banding assigns risk level ${rl} (${BANDCTRL[rl]}), above the level of control implied by the present capture efficiency.`);
  if (inp.limitSrc !== 'user') w.push('The built-in limit is a widely cited recommended exposure limit: verify it against current national regulation before use.');
  if (num(peak) > 1e7) w.push('The number concentration quoted assumes airborne primary particles; at such concentrations coagulation rapidly reduces the count (see the aerosol mode).');
  return {
    summary: `With ${F(inp.G, 3)} mg/min emitted for ${F(tT, 3)} min, ${inp.lev} % capture at source and ${F(lam * 60, 3)} h⁻¹ effective air cleaning, the near-field concentration approaches ${F(ssN, 3)} mg/m³ (far field ${F(ssF, 3)} mg/m³), about ${F(num(peak), 3)} particles/cm³ of ${F(inp.dp * 1e9, 3)} nm at the peak. The 8-h time-weighted average for a worker at the source is ${F(twa2, 3)} mg/m³ (one-box model ${F(twa1, 3)} mg/m³), ${F(hq, 3)} times the limit of ${F(oel, 3)} mg/m³. Control banding gives risk level ${rl}: ${BANDCTRL[rl]}. These are model estimates to be confirmed by workplace monitoring.`,
    kpis: [
      kp('twa', '8-h TWA (near-field/far-field model)', twa2, 'mg/m³', 'Time-weighted average over 480 min for a worker in the near field during the task and in the far field for the rest of the shift, from the two-zone model.', hq < 0.5 ? true : hq > 1 ? false : null),
      kp('twa_box', '8-h TWA (one-box model)', twa1, 'mg/m³', `Same average from the well-mixed single-zone model; it ignores the concentration gradient near the source and here predicts ${F(twa1 / Math.max(twa2, 1e-300), 3)} times the two-zone value.`, null),
      kp('hq', 'Exposure / limit', hq, '–', `Ratio of the two-zone 8-h TWA to the limit of ${F(oel, 3)} mg/m³ (${L.label}). Good practice aims below 0.5; above 1 the limit is predicted to be exceeded.`, hq < 0.5 ? true : hq > 1 ? false : null),
      kp('c_nf', 'Near-field steady-state concentration', ssN, 'mg/m³', 'Concentration in the breathing zone if the task continued indefinitely: far-field value plus G′/β, the local increment that general ventilation cannot remove.', null),
      kp('c_ff', 'Far-field steady-state concentration', ssF, 'mg/m³', 'Well-mixed room concentration at steady state, G′/Q_eff, to which other occupants of the room are exposed.', null),
      kp('n_peak', 'Peak number concentration (near field)', num(peak), '1/cm³', `Peak mass concentration ${F(peak, 3)} mg/m³ converted with the mass of one ${F(inp.dp * 1e9, 3)} nm particle; an upper bound, because agglomerates carry more mass per particle.`, null),
      kp('apf', 'Assigned protection factor needed', Math.max(hq, 1), '–', hq > 1 ? `Factor by which a respirator must reduce the inhaled concentration; the smallest conventional class that covers it is APF ${apf ? apf[0] : '> 10000'}${apf ? ' (' + apf[1] + ')' : ''}. Verify respirator classes nationally.` : 'No respiratory protection factor is needed to meet the limit on these estimates (value 1 means none).', hq <= 1),
      kp('lev_need', 'Capture efficiency for TWA ≤ 50 % of limit', 100 * needLev, '%', 'Local exhaust capture efficiency required to bring the two-zone TWA to half the limit, all else unchanged; exposure scales with the uncaptured fraction.', needLev * 100 <= inp.lev + 1e-9),
      kp('lam', 'Effective air-cleaning rate', lam * 60, '1/h', `Sum of fresh-air changes (${inp.ach}), filtered recirculation (${F(inp.recirc * inp.filtEff / 100, 3)}) and surface deposition (${inp.kdep}) per hour.`, null),
      kp('t_clear', 'Time to clear 99 % after the task', tClear, 'min', 'Time for the room concentration to fall to 1 % of its value after emission stops (ln 100 divided by the cleaning rate): the waiting time before re-entry without protection.', tClear < 60),
      kp('band', 'Control band', `RL${rl}: ${BANDCTRL[rl]}`, '', `Risk level from the 4 × 4 control-banding matrix: severity band "${inp.sevBand}" and exposure probability "${pName}" (score ${sc} of 100). It prescribes the minimum engineering control independently of the concentration estimate.`, rl <= present),
    ],
    series: [
      { key: 'ct', title: 'Airborne mass concentration over the shift', xlabel: 'time (min)', ylabel: 'concentration (mg/m³)', lines: [{ name: 'near field (breathing zone at source)', x: ts, y: Cn }, { name: 'far field (room)', x: ts, y: Cf }, { name: 'one-box model', x: ts, y: Cb }, { name: 'exposure limit', x: [0, 480], y: [oel, oel] }], explain: 'Concentration histories of the two-zone and one-box models over an 8-h shift. The near field rises within minutes to well above the room average and collapses when the task ends; the room then clears at the ventilation rate.' },
      { key: 'ach', title: '8-h TWA versus general ventilation', xlabel: 'air changes per hour', ylabel: '8-h TWA (mg/m³)', logx: true, logy: true, lines: (() => { const xs = N.logspace(-1, Math.log10(60), 30), f = (a, levv) => { const l2 = (a + inp.recirc * inp.filtEff / 100 + inp.kdep) / 60, g2 = inp.G * (1 - levv), far = g2 / (l2 * V), tOn = tT, ff = (far * (tOn - (1 - Math.exp(-l2 * tOn)) / l2) + far * (1 - Math.exp(-l2 * tOn)) * (1 - Math.exp(-l2 * (480 - tOn))) / l2) / 480; return Math.max(ff + g2 / b * tOn / 480, 1e-300); }; return [{ name: 'present capture efficiency', x: xs, y: xs.map(a => f(a, inp.lev / 100)) }, { name: 'no local exhaust', x: xs, y: xs.map(a => f(a, 0)) }, { name: 'exposure limit', x: [0.1, 60], y: [oel, oel] }]; })(), explain: 'Shift-average exposure against room air-change rate, with and without capture at source (quasi-steady near-field increment). General ventilation lowers the room background but cannot remove the near-field increment; capture at source acts on both.' },
    ],
    tables: [
      { key: 'lim', title: 'Built-in recommended exposure limits', columns: ['Material', '8-h TWA limit (mg/m³)', 'Status'], rows: [['Carbon nanotubes and nanofibres (elemental carbon)', 0.001, 'NIOSH recommended exposure limit'], ['Ultrafine (nanoscale) titanium dioxide', 0.3, 'NIOSH recommended exposure limit'], ['Fine titanium dioxide', 2.4, 'NIOSH recommended exposure limit'], ['Any other nanomaterial', 'user entry', 'obtain from national regulation or the safety data sheet']], explain: 'Small set of widely cited recommended limits offered as defaults. Limits differ between countries and change over time, so each value must be verified against current national regulation.' },
      { key: 'cb', title: 'Control-banding matrix (risk level)', columns: ['Severity \\ Probability', 'Extremely unlikely', 'Less likely', 'Likely', 'Probable'], rows: ['Low', 'Medium', 'High', 'Very high'].map((nm, i) => [nm + (i === sev ? ' (this case)' : ''), ...BAND[i].map((v, j) => 'RL' + v + (i === sev && j === pb ? ' ◄' : ''))]), explain: 'Four-by-four matrix of the control-banding approach (CB Nanotool type). RL1: general ventilation; RL2: fume hood or local exhaust; RL3: containment; RL4: seek specialist advice. The probability score used here is a simplified version of the published factor list.' },
    ],
    equations: [
      { kind: 'governing', label: 'One-box (well-mixed) model', text: 'V·dC/dt = G·(1 − ε_LEV) − (Q + η_f·Q_r + k_dep·V)·C' },
      { kind: 'governing', label: 'Near-field / far-field model', text: 'V_N·dC_N/dt = G′ + β·(C_F − C_N);   V_F·dC_F/dt = β·(C_N − C_F) − Q_eff·C_F' },
      { kind: 'formula', label: 'Steady state and time-weighted average', text: 'C_F,ss = G′/Q_eff,  C_N,ss = C_F,ss + G′/β;   TWA = (1/480 min)·∫C_worker dt;   n = C/(ρ·π·d³/6)' },
    ],
    checks: [
      { name: 'Aerosol mass balance', pass: mbErr < 1e-6, detail: `emitted ${F(emitted, 4)} mg = removed ${F(removed, 4)} mg + airborne ${F(inv, 4)} mg (residual ${F(mbErr, 2)})` },
      { name: 'Predicted exposure below the limit', pass: hq <= 1, detail: `TWA/limit = ${F(hq, 3)}` },
      { name: 'Predicted exposure below half the limit', pass: hq <= 0.5, detail: 'margin for model and measurement uncertainty' },
      { name: 'Control level matches the control band', pass: rl <= present, detail: `band RL${rl}; present control equivalent to RL${present}` },
      { name: 'Two-zone and one-box models consistent in the far field', pass: Math.abs(Cf.at(-1) - Cb.at(-1)) <= 0.2 * Math.max(Cb.at(-1), 1e-12) + 1e-12, detail: `end of shift: ${F(Cf.at(-1), 3)} vs ${F(Cb.at(-1), 3)} mg/m³` },
    ],
    warnings: w,
    recommendations: [
      { priority: hq > 1 ? 'high' : hq > 0.5 ? 'medium' : 'low', text: hq > 0.5 ? `Raise capture at source from ${inp.lev} % to at least ${F(100 * needLev, 3)} % (ventilated enclosure or glove box instead of an open bench): this brings the predicted TWA to half the limit, which ${F(inp.ach * 2, 3)} air changes of general ventilation would not achieve.` : `Present controls keep the predicted TWA at ${pct(hq)} of the limit; maintain the exhaust performance (face-velocity checks) since exposure rises in proportion to the uncaptured fraction.`, basis: 'Near-field increment G′/β is unaffected by general ventilation.' },
      { priority: rl >= 3 ? 'high' : 'medium', text: `Apply the control band RL${rl} (${BANDCTRL[rl]}) as the minimum standard for this task${inp.form === 'dry' ? '; handling the material as a suspension or paste instead of dry powder would lower the probability band' : ''}.`, basis: `Control banding: severity "${inp.sevBand}", probability score ${sc}.` },
      { priority: 'medium', text: `Confirm the estimate with task-based monitoring (particle counter plus filter sampling) and wait about ${F(tClear, 3)} min after the task before removing protection or opening the enclosure.`, basis: 'Clearance time ln(100)/λ and the uncertainty of emission-rate estimates.' },
      { priority: 'medium', text: `Capture efficiency of ${inp.lev} % means ${F(inp.G * inp.lev / 100 * tT, 3)} mg per task goes to the exhaust filter: collect it on high-efficiency filters, dispose of them as nanomaterial waste and prefer wet transfer, which also prevents product loss to air.`, basis: 'Mass captured by local exhaust per task.', sustainability: true },
    ],
    publish: { 'risk.exposure_mg_m3': twa2, 'risk.index': hq, 'risk.band': rl, 'risk.number_conc_cm3': num(peak) },
    numerics: { method: 'Crank–Nicolson integration of the linear two-zone model (0.5 min step) with analytical one-box solution', dof: 2, steps: nS, iterations: nS, residual: mbErr, notes: 'Zones are perfectly mixed; coagulation and size change are not included; emission is constant during the task.' },
  };
}

// ---------- mode: aerosol behaviour ----------
function runAerosol(inp) {
  const w = [], T = inp.Tair + 273.15, d = inp.dp, mu = muAir(T), cc = slip(d, T), D = N.K_B * T * cc / (3 * Math.PI * mu * d), vs = inp.rhoP * d * d * 9.80665 * cc / (18 * mu), K = coagK(d, inp.rhoP, T), n0 = inp.N0 * 1e6;
  const tHalf = 2 / (K * n0), kn = 2 * mfpAir(T) / d, tauR = inp.rhoP * d * d * cc / (18 * mu), xd = Math.sqrt(2 * D), dep = lung(d * 1e6), tSet = inp.Hroom / vs;
  const ds = N.logspace(-9, -5, 80), vS = ds.map(x => inp.rhoP * x * x * 9.80665 * slip(x, T) / (18 * mu)), xD = ds.map(x => Math.sqrt(2 * N.K_B * T * slip(x, T) / (3 * Math.PI * mu * x)));
  let dx = 0; try { dx = Math.exp(N.brent(l => { const x = Math.exp(l); return Math.log(inp.rhoP * x * x * 9.80665 * slip(x, T) / (18 * mu)) - Math.log(Math.sqrt(2 * N.K_B * T * slip(x, T) / (3 * Math.PI * mu * x))); }, Math.log(1e-9), Math.log(1e-4))); } catch { dx = 0; }
  const air = inp.vent * 1e-3 * 480, dose = r => inp.cAir * air * r, lg = ds.map(x => lung(x * 1e6)), tt = N.linspace(0, 4 * 3600, 80), sumR = dep.ha + dep.tb + dep.al;
  if (d * 1e6 < 0.001 || d * 1e6 > 100) w.push('Particle size outside the 1 nm – 100 µm range of the deposition fits; fractions are evaluated at the range limit.');
  if (tHalf < 600) w.push(`Coagulation halves the number concentration in ${F(tHalf, 3)} s: the aerosol will not remain at this size and number; measured size distributions shift towards agglomerates.`);
  const mm = x => 1000 * x;
  return {
    summary: `A ${F(d * 1e9, 3)} nm particle (ρ = ${inp.rhoP} kg/m³) in air at ${inp.Tair} °C has a slip correction of ${F(cc, 3)}, a diffusion coefficient of ${F(D, 3)} m²/s and a settling velocity of ${F(vs, 3)} m/s: it needs ${F(tSet / 3600, 3)} h to settle ${inp.Hroom} m, so ${xd > vs ? 'diffusion dominates its motion and it stays airborne until removed by ventilation, coagulation or deposition' : 'gravitational settling dominates its motion'}. At ${F(inp.N0, 3)} cm⁻³ coagulation halves the number in ${F(tHalf / 60, 3)} min. Predicted deposition on inhalation is ${pct(dep.tot)} in total, of which ${pct(dep.al)} in the alveolar region (simplified ICRP fit).`,
    kpis: [
      kp('cc', 'Cunningham slip correction', cc, '–', `Correction of Stokes drag for particles comparable to the gas mean free path (Knudsen number ${F(kn, 3)}); values well above 1 mean the particle moves between gas molecules rather than through a continuum.`, null),
      kp('D', 'Diffusion coefficient', D, 'm²/s', 'Stokes–Einstein diffusion coefficient with slip correction; the root-mean-square Brownian displacement in one second is ' + F(xd * 1e6, 3) + ' µm.', null),
      kp('vs', 'Terminal settling velocity', vs, 'm/s', `Stokes settling velocity with slip correction (relaxation time ${F(tauR, 3)} s); in one second the particle falls ${F(vs * 1e6, 3)} µm.`, null),
      kp('t_settle', 'Settling time through the room height', tSet, 's', `Time to settle ${inp.Hroom} m in still air (${F(tSet / 3600, 3)} h). Values of many hours or days mean settling does not clean the air and deposits on floors are not the main sink.`, tSet < 3600),
      kp('K', 'Coagulation coefficient', K, 'm³/s', 'Brownian coagulation coefficient of equal spheres from the Fuchs interpolation between free-molecular and continuum regimes (convention dN/dt = −½KN²).', null),
      kp('t_half', 'Coagulation half-life', tHalf, 's', `Time for the number concentration to halve by self-coagulation at the initial ${F(inp.N0, 3)} cm⁻³, 2/(K·N₀); it scales inversely with concentration.`, tHalf > 3600 ? null : false),
      kp('d_cross', 'Diffusion/settling crossover diameter', dx, 'm', 'Diameter at which settling distance and Brownian displacement per second are equal for this density; smaller particles are transported mainly by diffusion.', null),
      kp('df_tot', 'Total deposition fraction', dep.tot, '–', 'Fraction of inhaled particles of this size deposited anywhere in the respiratory tract (simplified ICRP 66 fit for adults at light exercise, nose breathing).', null),
      kp('df_al', 'Alveolar deposition fraction', dep.al, '–', 'Fraction deposited in the gas-exchange region, from which insoluble nanoparticles clear slowly; it peaks near 20 nm.', dep.al < 0.2),
      kp('df_tb', 'Tracheobronchial deposition fraction', dep.tb, '–', 'Fraction deposited in the conducting airways, cleared mainly by mucociliary transport.', null),
      kp('df_ha', 'Head-airway deposition fraction', dep.ha, '–', `Fraction deposited in nose, mouth and throat, including the inhalable fraction ${F(dep.IF, 3)} for this size.`, null),
      kp('dose_al', 'Alveolar deposited mass per shift', dose(dep.al), 'mg', `Mass deposited in the alveolar region in 8 h at ${F(inp.cAir, 3)} mg/m³ and ${inp.vent} L/min (${F(air, 3)} m³ inhaled); total deposited mass is ${F(dose(dep.tot), 3)} mg.`, null),
    ],
    series: [
      { key: 'mot', title: 'Settling versus diffusion by particle size', xlabel: 'particle diameter (m)', ylabel: 'displacement in 1 s (m)', logx: true, logy: true, lines: [{ name: 'gravitational settling', x: ds, y: vS }, { name: 'Brownian displacement (rms)', x: ds, y: xD }, { name: 'this particle', x: [d, d], y: [Math.min(vs, xd), Math.max(vs, xd)] }], explain: 'Distance travelled in one second by settling and by Brownian motion. Below the crossover diameter diffusion to surfaces and airways dominates; above it particles settle and impact.' },
      { key: 'lung', title: 'Regional lung deposition versus size', xlabel: 'particle diameter (m)', ylabel: 'deposition fraction', logx: true, lines: [{ name: 'total', x: ds, y: lg.map(l => l.tot) }, { name: 'alveolar', x: ds, y: lg.map(l => l.al) }, { name: 'tracheobronchial', x: ds, y: lg.map(l => l.tb) }, { name: 'head airways', x: ds, y: lg.map(l => l.ha) }, { name: 'this particle (total)', x: [d], y: [dep.tot] }], explain: 'Deposition fractions of the simplified ICRP model across the size range. The minimum near 0.3 µm separates diffusion-driven deposition of nanoparticles from impaction and settling of coarse particles.' },
      { key: 'coag', title: 'Number concentration decay by coagulation', xlabel: 'time (s)', ylabel: 'number concentration (1/cm³)', logy: true, lines: [{ name: 'N(t), constant-size approximation', x: tt, y: tt.map(t => inp.N0 / (1 + K * n0 * t / 2)) }], explain: 'Decay of particle number by Brownian self-coagulation with a constant coagulation coefficient. Mass concentration is unchanged: the same mass ends up in fewer, larger agglomerates.' },
    ],
    tables: [{ key: 'size', title: 'Behaviour at reference sizes', columns: ['Diameter (nm)', 'Slip correction', 'Diffusion coeff. (m²/s)', 'Settling velocity (m/s)', 'Total deposition', 'Alveolar deposition'], rows: [5, 20, 50, 100, 300, 1000, 5000].map(nmD => { const x = nmD * 1e-9, c = slip(x, T), l = lung(nmD / 1000); return [nmD, F(c, 4), F(N.K_B * T * c / (3 * Math.PI * mu * x), 3), F(inp.rhoP * x * x * 9.80665 * c / (18 * mu), 3), F(l.tot, 3), F(l.al, 3)]; }), explain: 'Transport properties and deposition fractions for a ladder of particle sizes at the stated density and temperature, for comparison with the selected particle.' }],
    equations: [
      { kind: 'formula', label: 'Slip correction, diffusion and settling', text: 'C_c = 1 + Kn·(1.257 + 0.4·e^(−1.1/Kn)),  Kn = 2λ/d;   D = k_B·T·C_c/(3πμd);   v_s = ρ_p·d²·g·C_c/(18μ)' },
      { kind: 'governing', label: 'Monodisperse coagulation', text: 'dN/dt = −½·K·N²  ⇒  N(t) = N₀/(1 + K·N₀·t/2),  t½ = 2/(K·N₀)' },
      { kind: 'formula', label: 'Regional deposition (simplified ICRP 66 fit, d in µm)', text: 'DF_AL = (0.0155/d)·[exp(−0.416(ln d + 2.84)²) + 19.11·exp(−0.482(ln d − 1.362)²)];  DF_TB = (0.00352/d)·[exp(−0.234(ln d + 3.40)²) + 63.9·exp(−0.819(ln d − 1.61)²)]', note: 'Fits of the ICRP Publication 66 model as given by Hinds (Aerosol Technology, 2nd ed.) for adults, light exercise, nose breathing; unit-density spheres.' },
      { kind: 'formula', label: 'Deposited dose', text: 'm_dep = C·V̇·t·DF' },
    ],
    checks: [
      { name: 'Deposition fractions within [0, 1]', pass: [dep.tot, dep.al, dep.tb, dep.ha].every(v => v >= 0 && v <= 1), detail: `total ${F(dep.tot, 3)}, AL ${F(dep.al, 3)}, TB ${F(dep.tb, 3)}, HA ${F(dep.ha, 3)}` },
      { name: 'Regional fractions consistent with total', pass: Math.abs(sumR - dep.tot) < 0.08, detail: `HA + TB + AL = ${F(sumR, 3)} vs total ${F(dep.tot, 3)} (independent fits)` },
      { name: 'Stokes regime (particle Reynolds number < 0.1)', pass: 1.2 * vs * d / mu < 0.1, detail: `Re_p = ${F(1.2 * vs * d / mu, 2)}` },
      { name: 'Aerosol number stable over one hour', pass: tHalf > 3600, detail: `coagulation half-life ${F(tHalf, 3)} s` },
    ],
    warnings: w,
    recommendations: [
      { priority: dep.al > 0.2 ? 'high' : 'medium', text: dep.al > 0.2 ? `Particles of ${F(d * 1e9, 3)} nm deposit efficiently in the alveolar region (${pct(dep.al)}): control at source and use high-efficiency particulate filtration; settling and coarse-dust controls are ineffective at this size.` : `Alveolar deposition is ${pct(dep.al)} at this size; ${d > 1e-6 ? 'head-airway deposition dominates, so inhalable-dust controls and good housekeeping are the priority' : 'maintain source control and filtration'}.`, basis: 'Regional deposition fractions of the simplified ICRP model.' },
      { priority: 'medium', text: tSet > 3600 ? `Do not rely on settling to clear the air: the particles need ${F(tSet / 3600, 3)} h to fall ${inp.Hroom} m. Clearance must come from ventilation and filtration.` : `Particles settle within ${F(tSet / 60, 3)} min: wet-wipe or vacuum (high-efficiency filter) surfaces, since settled material is resuspended by dry sweeping.`, basis: 'Settling time through the room height.' },
      { priority: 'medium', text: `Report exposure with both mass and number metrics: at ${F(inp.N0, 3)} cm⁻³ the count halves in ${F(tHalf / 60, 3)} min by coagulation while the mass is conserved, so number readings depend strongly on where and when they are taken.`, basis: 'Coagulation half-life 2/(K·N₀).' },
      { priority: 'low', text: 'Fit recirculating high-efficiency filter units with low-pressure-drop media and run them on demand (particle-counter interlock): diffusion capture of nanoparticles is efficient even at reduced face velocity, which saves fan energy.', basis: 'Diffusion-dominated transport below the crossover diameter.', sustainability: true },
    ],
    publish: { 'risk.df_alveolar': dep.al, 'risk.df_total': dep.tot, 'risk.dose_alveolar_mg': dose(dep.al), 'risk.settling_m_s': vs },
    numerics: { method: 'Closed-form aerosol mechanics (Stokes–Cunningham, Fuchs coagulation) and algebraic ICRP deposition fits; Brent root search for the crossover diameter', dof: 1, steps: 80, iterations: 1, residual: Math.abs(sumR - dep.tot), notes: 'Spherical particles, still air at 1 atm; deposition fits are population averages and take the entered diameter as the governing diameter.' },
  };
}

// ---------- mode: hydrothermal autoclave ----------
function runAutoclave(inp) {
  const w = [], T = inp.Tac + 273.15, T0 = 293.15, f = inp.fillPct / 100, V = inp.Vliner * 1e-6, r20 = rhoL(T0), rl = rhoL(T), rv = rhoV(T), ps = psatW(T), ravg = f * r20, fCrit = rl / r20;
  const state = fl => { const ra = fl * r20; if (ra < rl * 0.9995) { const phi = N.clamp((ra - rv) / (rl - rv), 0, 1), pAir = PATM * (T / T0) * (1 - fl) / Math.max(1 - phi, 1e-6); return { phi, P: ps + Math.min(pAir, 5e8), full: false }; } return { phi: 1, P: ps + Math.log(ra / rl) / kappaW(inp.Tac) + PATM * (T / T0) * (1 - fl) / 5e-4, full: true }; };
  const st = state(f), P = st.P, Pg = P - PATM;
  let Tfull = TC; if (ravg > RHOC) { try { Tfull = N.brent(x => rhoL(x) - ravg, T0, TC - 0.01); } catch { Tfull = T0; } }
  const ri = inp.IDmm / 2e3, ro = ri + inp.wall / 1e3, Kr = ro / ri, k2 = Kr * Kr, sth = Pg * (k2 + 1) / (k2 - 1), sz = Pg / (k2 - 1), svm = Math.sqrt(3) * Pg * k2 / (k2 - 1), sthin = Pg * ri / (inp.wall / 1e3), Y = inp.yieldT * 1e6;
  const sf = svm > 0 ? Y / svm : 99, Pyield = Y * (k2 - 1) / (Math.sqrt(3) * k2), Pcoll = 2 / Math.sqrt(3) * Y * Math.log(Kr), mawp = Pyield / inp.designSF, Pbd = Pg / 0.8;
  const mL = st.phi * V * rl, Vg = (1 - st.phi) * V, Tb = 373.15, eL = T > Tb ? mL * 4400 * ((T - Tb) - Tb * Math.log(T / Tb)) : 0, eV = Pg * Vg / 0.33, E = eL + eV, tnt = E / 4.184e6 * 1e3;
  if (st.full) w.push(`At ${inp.Tac} °C the vessel is liquid-full (critical fill ${pct(fCrit)}): the pressure is set by the compressibility of water and the elasticity of the vessel, not by the vapour pressure. The value quoted is an order-of-magnitude estimate; the vessel must be expected to fail or vent.`);
  else if (f > 0.9 * fCrit) w.push(`Fill fraction ${inp.fillPct} % is within 10 % of the critical fill ${pct(fCrit)} at ${inp.Tac} °C: a small overshoot in temperature or filling makes the vessel liquid-full.`);
  if (sf < inp.designSF) w.push(`Equivalent stress at the bore (${F(svm / 1e6, 3)} MPa) leaves a safety factor of only ${F(sf, 3)} against yield at temperature, below the required ${inp.designSF}.`);
  if (inp.Tac > 220) w.push('Above about 200–220 °C fluoropolymer liners creep and lose sealing force; use a liner material and closure rated for the temperature.');
  if (inp.Tac > 300) w.push('Above 300 °C water approaches its critical point (374 °C, 22.06 MPa): density and pressure change very steeply with temperature and fill.');
  w.push('Only the cylindrical shell is assessed: closures, threads, seals and the liner must be verified separately against the manufacturer\'s rating, and dissolved salts or decomposing reagents (gas generation) can raise the pressure above that of pure water.');
  const Ts = N.linspace(20, 350, 67), fs = N.linspace(0.05, 0.98, 60), rr = N.linspace(ri, ro, 30), selfP = psatW(373.15), selfP2 = psatW(473.15);
  const Pt = Ts.map(tc => { const TT = tc + 273.15, rL2 = rhoL(TT), rV2 = rhoV(TT); if (ravg >= rL2 * 0.9995) return null; const phi = N.clamp((ravg - rV2) / (rL2 - rV2), 0, 1); return (psatW(TT) + PATM * (TT / T0) * (1 - f) / Math.max(1 - phi, 1e-6)) / 1e6; });
  return {
    summary: `At ${inp.Tac} °C water has a saturation pressure of ${F(ps / 1e6, 4)} MPa; with ${inp.fillPct} % fill the liquid occupies ${pct(st.phi)} of the vessel and the total pressure (vapour plus heated trapped air) is ${F(P / 1e6, 4)} MPa. ${st.full ? 'The vessel is liquid-full, which is a hydraulic over-pressure condition.' : `The vessel becomes liquid-full at ${pct(fCrit)} fill, or at ${F(Tfull - 273.15, 4)} °C for the present fill.`} The thick-wall equivalent stress at the bore is ${F(svm / 1e6, 3)} MPa, a safety factor of ${F(sf, 3)} against yield; the shell starts to yield at ${F(Pyield / 1e6, 3)} MPa. Stored energy is ${F(E / 1e3, 3)} kJ, about ${F(tnt, 3)} g TNT equivalent.`,
    kpis: [
      kp('psat', 'Water saturation pressure', ps, 'Pa', `Vapour pressure of pure water at ${inp.Tac} °C from the IAPWS-IF97 saturation equation (${F(ps / 1e5, 4)} bar absolute).`, null),
      kp('P', 'Total internal pressure', P, 'Pa', st.full ? 'Estimated pressure of the liquid-full vessel from the isothermal compressibility of water; indicative only and far above any normal operating value.' : `Saturation pressure plus the partial pressure of the air sealed in at 20 °C, compressed into the remaining ${pct(1 - st.phi)} gas space and heated (ideal gas, air solubility neglected).`, P < mawp),
      kp('phi', 'Liquid volume fraction at temperature', st.phi, '–', `Fraction of the vessel occupied by liquid at ${inp.Tac} °C: saturated liquid density falls from ${F(r20, 4)} to ${F(rl, 4)} kg/m³, so the charge expands by ${pct(r20 / rl - 1)}.`, st.phi < 0.95),
      kp('fcrit', 'Critical fill fraction', fCrit, '–', `Fill at 20 °C above which the vessel is completely full of liquid at ${inp.Tac} °C (ρ′(T)/ρ(20 °C)); the present fill is ${pct(f / fCrit)} of it.`, f < 0.9 * fCrit),
      kp('Tfull', 'Temperature at which the vessel is liquid-full', Tfull - 273.15, '°C', ravg > RHOC ? 'Temperature at which the expanding liquid fills the vessel for the present fill; beyond it the pressure rises by roughly 1 MPa per kelvin.' : 'The mean density is below the critical density, so the liquid level falls on heating and the vessel never becomes liquid-full (critical temperature shown).', Tfull - T > 20),
      kp('svm', 'Equivalent stress at the bore', svm, 'Pa', `Von Mises stress at the inner wall of the closed-end thick cylinder (Lamé): hoop ${F(sth / 1e6, 3)} MPa, axial ${F(sz / 1e6, 3)} MPa, radial −${F(Pg / 1e6, 3)} MPa; the thin-wall hoop formula gives ${F(sthin / 1e6, 3)} MPa.`, sf >= inp.designSF),
      kp('sf', 'Safety factor against yield', sf, '–', `Yield strength at temperature (${inp.yieldT} MPa) divided by the equivalent stress; the required value is ${inp.designSF}.`, sf >= inp.designSF),
      kp('mawp', 'Allowable working pressure of the shell', mawp, 'Pa', `Gauge pressure at which the bore reaches yield (${F(Pyield / 1e6, 3)} MPa) divided by the required safety factor; full plastic collapse of the wall is estimated at ${F(Pcoll / 1e6, 3)} MPa.`, mawp >= Pg),
      kp('pbd', 'Burst-disc setting', Pbd, 'Pa', `Disc rating that keeps the operating pressure at 80 % of the burst pressure (${F(Pbd / 1e6, 3)} MPa gauge); it must not exceed the allowable working pressure of ${F(mawp / 1e6, 3)} MPa.`, Pbd <= mawp),
      kp('E', 'Stored energy', E, 'J', `Energy released on sudden rupture: exergy of the superheated liquid (${F(eL, 3)} J) plus expansion of the vapour space (${F(eV, 3)} J, Brode).`, E < 1e4),
      kp('tnt', 'TNT equivalent', tnt, 'g', 'Stored energy divided by 4.184 MJ/kg; a screening figure for shielding and separation distance, not a blast calculation.', tnt < 5),
    ],
    series: [
      { key: 'pt', title: 'Pressure versus temperature', xlabel: 'temperature (°C)', ylabel: 'pressure (MPa)', logy: true, lines: [{ name: 'water saturation pressure', x: Ts, y: Ts.map(tc => psatW(tc + 273.15) / 1e6) }, { name: `total pressure at ${inp.fillPct} % fill`, x: Ts.filter((_, i) => Pt[i] !== null), y: Pt.filter(v => v !== null) }, { name: 'allowable working pressure', x: [20, 350], y: [Math.max(mawp + PATM, 1) / 1e6, Math.max(mawp + PATM, 1) / 1e6] }, { name: 'operating point', x: [inp.Tac], y: [P / 1e6] }], explain: 'Saturation pressure of water and total pressure in the sealed vessel against temperature. The total-pressure curve ends where the liquid fills the vessel; beyond that point pressure rises almost vertically.' },
      { key: 'pf', title: 'Pressure versus fill fraction', xlabel: 'fill fraction at 20 °C', ylabel: 'pressure (MPa)', logy: true, lines: [{ name: `pressure at ${inp.Tac} °C`, x: fs, y: fs.map(x => state(x).P / 1e6) }, { name: 'critical fill', x: [fCrit, fCrit], y: [ps / 1e6, Math.max(state(0.98).P, 2 * ps) / 1e6] }, { name: 'this vessel', x: [f], y: [P / 1e6] }], explain: 'Pressure at the operating temperature as a function of how full the vessel was loaded. Below the critical fill the pressure is close to the vapour pressure; above it the vessel is liquid-full and the pressure is hydraulic.' },
      { key: 'lame', title: 'Stress through the vessel wall', xlabel: 'radius (m)', ylabel: 'stress (MPa)', lines: [{ name: 'hoop stress', x: rr, y: rr.map(r => Pg / (k2 - 1) * (1 + ro * ro / (r * r)) / 1e6) }, { name: 'radial stress', x: rr, y: rr.map(r => Pg / (k2 - 1) * (1 - ro * ro / (r * r)) / 1e6) }, { name: 'axial stress', x: [ri, ro], y: [sz / 1e6, sz / 1e6] }], explain: 'Lamé stress distribution in the closed-end cylinder. Hoop stress is highest at the bore and the radial stress equals minus the pressure there, so yielding always starts at the inner surface.' },
    ],
    tables: [{ key: 'st', title: 'Water properties on the saturation line', columns: ['T (°C)', 'p_sat (MPa)', 'ρ liquid (kg/m³)', 'ρ vapour (kg/m³)', 'Critical fill (%)'], rows: [100, 120, 150, 180, 200, 220, 250, 300].map(tc => { const TT = tc + 273.15; return [tc, F(psatW(TT) / 1e6, 4), F(rhoL(TT), 4), F(rhoV(TT), 4), F(100 * rhoL(TT) / r20, 4)]; }), explain: 'Saturation pressure and saturated densities of water from the IAPWS correlations, with the fill fraction at 20 °C that makes a sealed vessel liquid-full at each temperature.' }],
    equations: [
      { kind: 'chemical', label: 'Phase equilibrium in the sealed vessel', text: `${chem.prettyFormula('H2O(l)')} ⇌ ${chem.prettyFormula('H2O(g)')}`, note: 'Pure water; dissolved electrolytes lower the vapour pressure slightly, gas-forming reactions raise the total pressure.' },
      { kind: 'constitutive', label: 'Saturation pressure (IAPWS-IF97, region 4)', text: 'p_sat = [2C/(−B + √(B² − 4AC))]⁴ MPa,  A, B, C quadratic in ϑ = T + n₉/(T − n₁₀)' },
      { kind: 'formula', label: 'Fill and trapped air', text: 'φ_L = (f·ρ₂₀ − ρ″)/(ρ′ − ρ″);  f_crit = ρ′(T)/ρ₂₀;  p_air = p₀·(T/T₀)·(1 − f)/(1 − φ_L);  liquid-full: Δp ≈ ln(ρ̄/ρ′)/κ_T' },
      { kind: 'governing', label: 'Thick-wall cylinder (Lamé), closed ends', text: 'σ_θ = p·(rᵢ²/(rₒ² − rᵢ²))·(1 + rₒ²/r²),  σ_r = p·(rᵢ²/(rₒ² − rᵢ²))·(1 − rₒ²/r²),  σ_z = p·rᵢ²/(rₒ² − rᵢ²);  σ_vM(rᵢ) = √3·p·K²/(K² − 1),  K = rₒ/rᵢ' },
      { kind: 'formula', label: 'Stored energy', text: 'E = m_L·c_p·[(T − T_b) − T_b·ln(T/T_b)] + (p − p₀)·V_g/(γ − 1);   m_TNT = E/4.184 MJ·kg⁻¹' },
    ],
    checks: [
      { name: 'Saturation correlation reproduces reference points', pass: Math.abs(selfP - 101418) < 150 && Math.abs(selfP2 - 1.5549e6) < 2000, detail: `p_sat(100 °C) = ${F(selfP / 1e3, 6)} kPa (101.42), p_sat(200 °C) = ${F(selfP2 / 1e6, 5)} MPa (1.5549)` },
      { name: 'Vessel not liquid-full at temperature', pass: !st.full, detail: `fill ${inp.fillPct} % vs critical ${pct(fCrit)}` },
      { name: 'Fill at least 10 % below critical fill', pass: f <= 0.9 * fCrit, detail: `margin ${pct(1 - f / fCrit)}` },
      { name: 'Safety factor against yield adequate', pass: sf >= inp.designSF, detail: `${F(sf, 3)} vs required ${inp.designSF}` },
      { name: 'Burst-disc setting below allowable working pressure', pass: Pbd <= mawp, detail: `${F(Pbd / 1e6, 3)} MPa vs ${F(mawp / 1e6, 3)} MPa (gauge)` },
      { name: 'Below the critical point of water', pass: T < TC - 20, detail: `T = ${F(T, 5)} K, T_c = 647.096 K` },
    ],
    warnings: w,
    recommendations: [
      { priority: st.full || f > 0.9 * fCrit ? 'high' : 'medium', text: st.full || f > 0.9 * fCrit ? `Reduce the fill to at most ${F(100 * Math.min(0.8, 0.85 * fCrit), 3)} % at 20 °C: the present ${inp.fillPct} % is ${st.full ? 'above' : 'too close to'} the critical fill of ${pct(fCrit)} at ${inp.Tac} °C.` : `Keep the fill at or below ${F(100 * Math.min(0.8, 0.85 * fCrit), 3)} % and mark the liner with the maximum level; the margin to liquid-full is ${F(Tfull - T, 3)} K in temperature.`, basis: 'Thermal expansion of saturated liquid water against the critical fill fraction.' },
      { priority: sf < inp.designSF ? 'high' : 'low', text: sf < inp.designSF ? `Do not operate at ${inp.Tac} °C: the shell needs a wall of at least ${F(1e3 * ri * (Math.sqrt(1 / (1 - Math.min(Math.sqrt(3) * Pg * inp.designSF / Y, 0.98))) - 1), 3)} mm for the required safety factor, or the temperature must be lowered until the pressure is below ${F(mawp / 1e6, 3)} MPa.` : `The shell has a safety factor of ${F(sf, 3)}; confirm that the closure and thread rating from the manufacturer also exceeds ${F(P / 1e6, 3)} MPa at ${inp.Tac} °C.`, basis: 'Lamé equivalent stress at the bore against yield strength at temperature.' },
      { priority: 'high', text: `Fit a burst disc or relief device rated near ${F(Pbd / 1e6, 3)} MPa (gauge), never above ${F(mawp / 1e6, 3)} MPa, heat the vessel in an oven with independent over-temperature cut-out and let it cool below 40 °C before opening; the contents store about ${F(tnt, 3)} g TNT equivalent.`, basis: 'Operating pressure, allowable working pressure and stored energy.' },
      { priority: 'low', text: `Lower-temperature routes save energy and reduce hazard together: at ${Math.max(inp.Tac - 40, 100)} °C the pressure falls to ${F(psatW(T - 40) / 1e6, 3)} MPa. Evaluate whether longer time or a mineraliser reaches the same product, and insulate the oven or use a heated block for small vessels.`, basis: 'Steep temperature dependence of the saturation pressure.', sustainability: true },
    ],
    publish: { 'proc.P': P, 'proc.T_C': inp.Tac, 'mech.max_stress': svm, 'mech.safety_factor': sf, 'risk.index': inp.designSF / Math.max(sf, 1e-9), 'risk.stored_energy_J': E },
    numerics: { method: 'IAPWS-IF97 saturation pressure and IAPWS saturated-density correlations; lever rule for the two-phase fill; Lamé thick-wall solution; Brent root search for the liquid-full temperature', dof: 1, steps: 67, iterations: 1, residual: Math.abs(selfP - 101418) / 101418, notes: 'Pure water, thermodynamic equilibrium, rigid vessel, ideal-gas air; closure and thread strength are not modelled.' },
  };
}

// ---------- mode: flammable solvent, dust, electrical ----------
function runHazards(inp) {
  const w = [], s = SOLV[inp.solvent] || SOLV.ethanol, Tk = inp.Tsolv + 273.15, ps = antoine(s, inp.Tsolv), sat = Math.min(100, 100 * ps / PATM), M = chem.molarMass(s.f) / 1e3, nSolv = inp.Vsp * 1e-3 * s.rho / M;
  const cb = chem.balance([s.f, 'O2'], ['CO2', 'H2O']), nuO2 = cb.reactants[1] / cb.reactants[0], cst = 100 / (1 + 4.76 * nuO2), lelJ = 0.55 * cst;
  let tLel = s.fp; try { tLel = N.brent(tc => 100 * antoine(s, tc) / PATM - s.lel, -70, 250); } catch { /* keep tabulated flash point */ }
  const volAll = Math.min(sat, 100 * nSolv * R * Tk / (PATM * inp.Venc)), fracLel = volAll / s.lel, qv = nSolv / (inp.tEvap * 60) * R * Tk / PATM, Qv = inp.achEnc * inp.Venc / 3600, volSS = 100 * qv / (Qv + qv), fracSS = volSS / s.lel, qNeed = qv * (100 / (0.25 * s.lel) - 1), achNeed = qNeed * 3600 / inp.Venc;
  const above = inp.Tsolv >= s.fp, inRange = sat >= s.lel && sat <= s.uel;
  if (inp.Tsolv > s.bp) w.push(`${s.label} is above its normal boiling point (${s.bp} °C): the Antoine fit is extrapolated and a pressure vessel is needed; any release flashes to vapour.`);
  if (above) w.push(`Liquid temperature ${inp.Tsolv} °C is at or above the flash point of ${s.label.toLowerCase()} (${s.fp} °C): the vapour above the liquid is ignitable (saturated concentration ${F(sat, 3)} vol %, flammable range ${s.lel}–${s.uel} vol %).`);
  if (fracLel > 0.25) w.push(`Complete evaporation of ${inp.Vsp} mL in ${inp.Venc} m³ gives ${F(volAll, 3)} vol %, i.e. ${pct(fracLel)} of the lower flammable limit: above the customary 25 % criterion.`);
  // dust
  const dm = DUST[inp.dustMat] || DUST.ZnO, cDust = inp.mDust / inp.Venc, o2 = 0.2095 * PATM / (R * 293.15); let dP = 0, burn = 0, dEq = null, comb = !!dm.ox;
  if (comb) { const db = chem.balance([inp.dustMat, 'O2'], [dm.ox]), nEl = cDust / chem.molarMass(inp.dustMat); burn = Math.min(nEl, o2 * db.reactants[0] / db.reactants[1]); dP = 0.4 * burn * (db.products[0] / db.reactants[0]) * Math.abs(dm.dH) * 1e3; dEq = chem.equationString([inp.dustMat, 'O2'], [dm.ox], db.coefficients); }
  const explosible = comb && cDust >= inp.mec;
  if (explosible) w.push(`${dm.label}: a dispersed cloud of ${F(cDust, 3)} g/m³ exceeds the minimum explosible concentration entered (${inp.mec} g/m³). Metal and carbon nanopowders can ignite at very low energies (static discharge) and some are pyrophoric; obtain dust-explosion test data for the actual powder.`);
  else if (comb) w.push(`${dm.label} is combustible; the cloud concentration of ${F(cDust, 3)} g/m³ is below the entered minimum explosible concentration, but deposits can still burn and local concentrations near a release are higher.`);
  // electrical
  const E = inp.Uapp / (inp.gap * 1e-3), Rc = inp.gap * 1e-3 / (inp.sigEl * inp.Ael * 1e-4), I = inp.Uapp / Rc, P = inp.Uapp * I, rate = P / (1000 * 4180 * inp.Vliq * 1e-6), tBoil = rate > 0 ? Math.max(0, 100 - 25) / rate : 1e12, uLim = inp.acdc === 'ac' ? 50 : 120, iBody = inp.Uapp / 1000 * 1e3;
  const nH2 = inp.acdc === 'dc' ? I / (2 * N.FARADAY) : 0, vH2 = nH2 * R * 293.15 / PATM, h2ss = 100 * vH2 / (Qv + vH2), tH2 = vH2 > 0 ? 0.04 * inp.Venc / vH2 : 1e12;
  if (inp.Uapp > uLim) w.push(`Applied voltage ${F(inp.Uapp, 3)} V exceeds the conventional extra-low-voltage touch limit (${uLim} V ${inp.acdc.toUpperCase()}): live parts must be insulated or enclosed with interlocks. The body current through 1 kΩ would be ${F(iBody, 3)} mA.`);
  if (tBoil < 600) w.push(`Joule heating of ${F(P, 3)} W raises ${inp.Vliq} mL by ${F(rate, 3)} K/s: the liquid would boil in ${F(tBoil, 3)} s without cooling or current limit.`);
  if (h2ss >= 1) w.push(`Electrolytic hydrogen reaches ${F(h2ss, 3)} vol % in the enclosure at steady state (lower flammable limit 4 vol %): provide ventilation or inerting.`);
  const Ts = N.linspace(-20, Math.min(s.bp, 120), 60), ratio = Math.max(fracLel / 0.25, fracSS / 0.25, explosible ? 1.5 : 0, inp.Uapp / uLim, h2ss / 1);
  return {
    summary: `${s.label} at ${inp.Tsolv} °C has a vapour pressure of ${F(ps / 1e3, 3)} kPa (saturated ${F(sat, 3)} vol %; flammable range ${s.lel}–${s.uel} vol %, flash point ${s.fp} °C). Evaporating ${inp.Vsp} mL into ${inp.Venc} m³ gives ${pct(fracLel)} of the lower flammable limit unventilated and ${pct(fracSS)} at ${inp.achEnc} air changes per hour. ${dm.label}: ${comb ? `combustible, ${F(cDust, 3)} g/m³ if dispersed, adiabatic over-pressure bound ${F(dP / 1e5, 3)} bar` : 'not combustible, no dust-explosion hazard from the powder itself'}. The cell draws ${F(I, 3)} A and ${F(P, 3)} W at ${F(inp.Uapp, 3)} V (${F(E, 3)} V/m), heating the liquid at ${F(rate * 60, 3)} K/min.`,
    kpis: [
      kp('psat', 'Solvent vapour pressure', ps, 'Pa', `Antoine vapour pressure of ${s.label.toLowerCase()} at ${inp.Tsolv} °C; the saturated vapour concentration above the liquid is ${F(sat, 3)} vol %.`, sat < s.lel),
      kp('t_lel', 'Temperature at which saturation reaches the LEL', tLel, '°C', `Liquid temperature at which the equilibrium vapour concentration equals the lower flammable limit of ${s.lel} vol %: a thermodynamic estimate of the flash point (tabulated closed-cup value ${s.fp} °C).`, inp.Tsolv < tLel - 5),
      kp('frac_lel', 'Unventilated release as fraction of LEL', fracLel, '–', `Vapour concentration (${F(volAll, 3)} vol %) after all ${inp.Vsp} mL evaporates into the enclosure, divided by the lower flammable limit; 0.25 is the customary upper bound for safe operation.`, fracLel <= 0.25),
      kp('frac_ss', 'Ventilated steady state as fraction of LEL', fracSS, '–', `Steady vapour concentration (${F(volSS, 3)} vol %) when the solvent evaporates over ${inp.tEvap} min into the ventilated enclosure, divided by the lower flammable limit.`, fracSS <= 0.25),
      kp('ach_need', 'Ventilation for 25 % of LEL', achNeed, '1/h', `Air changes per hour of the enclosure that dilute the evaporating vapour to a quarter of the lower flammable limit (${F(qNeed * 1e3, 3)} L/s of air).`, achNeed <= inp.achEnc),
      kp('lel_jones', 'LEL estimated from stoichiometry', lelJ, 'vol %', `Jones rule: 0.55 × stoichiometric concentration (${F(cst, 3)} vol %), from the balanced combustion equation; compared with the tabulated ${s.lel} vol % as a plausibility check.`, Math.abs(lelJ - s.lel) / s.lel < 0.35),
      kp('c_dust', 'Dispersed dust concentration', cDust, 'g/m³', `Concentration if ${inp.mDust} g of powder were suspended uniformly in the enclosure, compared with the minimum explosible concentration entered (${inp.mec} g/m³).`, !explosible),
      kp('dp_dust', 'Adiabatic dust over-pressure bound', dP, 'Pa', comb ? `Upper bound of the pressure rise for constant-volume, adiabatic combustion of the oxygen-limited dust (${F(burn, 3)} mol/m³ burnt), Δp = (γ − 1)·q; real explosions reach lower values.` : 'The selected powder is fully oxidised or noble and does not release combustion heat; the bound is zero.', dP < 1e4),
      kp('I', 'Cell current', I, 'A', `Ohmic current through the electrolyte of resistance ${F(Rc, 3)} Ω between plane electrodes (uniform field, electrode polarisation neglected).`, null),
      kp('P', 'Joule heating power', P, 'W', `Electrical power dissipated in the liquid (volumetric ${F(inp.sigEl * E * E, 3)} W/m³); it heats ${inp.Vliq} mL at ${F(rate, 3)} K/s adiabatically.`, tBoil > 600),
      kp('touch', 'Voltage / touch-voltage limit', inp.Uapp / uLim, '–', `Applied voltage relative to the conventional extra-low-voltage limit of ${uLim} V for ${inp.acdc === 'ac' ? 'alternating' : 'ripple-free direct'} current; above 1 protection against direct contact is mandatory.`, inp.Uapp <= uLim),
      kp('h2', 'Electrolytic hydrogen generation', vH2 * 6e7, 'mL/min', inp.acdc === 'dc' ? `Faraday-law hydrogen evolution at the cathode in aqueous electrolyte (upper bound, 100 % current efficiency); unventilated, the enclosure reaches 4 vol % after ${F(tH2 / 60, 3)} min.` : 'No net electrolysis gas is produced with alternating current at mains frequency in this estimate.', h2ss < 1),
    ],
    series: [
      { key: 'vp', title: 'Saturated vapour concentration versus temperature', xlabel: 'liquid temperature (°C)', ylabel: 'vapour concentration (vol %)', logy: true, lines: [{ name: 'saturated vapour', x: Ts, y: Ts.map(tc => Math.min(100, 100 * antoine(s, tc) / PATM)) }, { name: 'lower flammable limit', x: [Ts[0], Ts.at(-1)], y: [s.lel, s.lel] }, { name: 'upper flammable limit', x: [Ts[0], Ts.at(-1)], y: [s.uel, s.uel] }, { name: 'operating point', x: [N.clamp(inp.Tsolv, Ts[0], Ts.at(-1))], y: [Math.max(sat, 1e-6)] }], explain: 'Equilibrium vapour concentration above the liquid against temperature. Where the curve lies between the flammable limits the head space is ignitable; above the upper limit it becomes ignitable as soon as it is diluted with air.' },
      { key: 'vent', title: 'Vapour concentration versus enclosure ventilation', xlabel: 'air changes per hour', ylabel: 'fraction of LEL', logx: true, logy: true, lines: (() => { const xs = N.logspace(-1, Math.log10(600), 40); return [{ name: 'steady-state vapour', x: xs, y: xs.map(a => Math.max(100 * qv / (a * inp.Venc / 3600 + qv) / s.lel, 1e-12)) }, { name: '25 % of LEL', x: [0.1, 600], y: [0.25, 0.25] }]; })(), explain: 'Steady vapour concentration, as a fraction of the lower flammable limit, against the air-change rate of the enclosure. The crossing with the 25 % line gives the minimum ventilation.' },
      { key: 'joule', title: 'Adiabatic temperature rise by Joule heating', xlabel: 'time (s)', ylabel: 'temperature (°C)', lines: (() => { const te = Math.min(Math.max(tBoil * 1.2, 60), 36000), xs = N.linspace(0, te, 50); return [{ name: 'liquid temperature (no cooling)', x: xs, y: xs.map(t => Math.min(25 + rate * t, 100)) }, { name: 'boiling point of water', x: [0, te], y: [100, 100] }]; })(), explain: 'Temperature of the liquid if all electrical power stays in it. The slope shows how quickly a field-assisted reactor overheats when cooling or current limiting fails.' },
    ],
    tables: [{ key: 'solv', title: 'Solvent flammability data (verify against the safety data sheet)', columns: ['Solvent', 'Formula', 'Flash point (°C)', 'LEL (vol %)', 'UEL (vol %)', 'Boiling point (°C)', `Saturated vapour at ${inp.Tsolv} °C (vol %)`], rows: Object.values(SOLV).map(x => [x.label, chem.prettyFormula(x.f), x.fp, x.lel, x.uel, x.bp, F(Math.min(100, 100 * antoine(x, inp.Tsolv) / PATM), 3)]), explain: 'Handbook flammability data of common process solvents with the saturated vapour concentration at the entered temperature, for choosing a less hazardous solvent.' }],
    equations: [
      { kind: 'chemical', label: `Combustion of ${s.label.toLowerCase()}`, text: chem.equationString([s.f, 'O2'], ['CO2', 'H2O'], cb.coefficients), note: `Stoichiometric concentration in air ${F(cst, 3)} vol %.` },
      ...(dEq ? [{ kind: 'chemical', label: 'Dust combustion', text: dEq, note: `ΔH_f of the oxide ${dm.dH} kJ/mol.` }] : []),
      { kind: 'half-cell', label: 'Cathode (aqueous electrolyte)', text: chem.equationString(['H2O', 'e-'], ['H2', 'OH-']) },
      { kind: 'half-cell', label: 'Anode (aqueous electrolyte)', text: chem.equationString(['OH-'], ['O2', 'H2O', 'e-']) },
      { kind: 'constitutive', label: 'Antoine vapour pressure', text: `log₁₀ p[mmHg] = A − B/(C + T[°C]),  A = ${s.A}, B = ${s.B}, C = ${s.C}` },
      { kind: 'formula', label: 'Dilution ventilation', text: 'x_ss = q_v/(Q + q_v),  q_v = ṅ·R·T/p;   Q(25 % LEL) = q_v·(1/(0.25·LEL) − 1)' },
      { kind: 'formula', label: 'Ohmic cell and Joule heating', text: 'R = g/(σ·A),  I = U/R,  P = U·I = σ·E²·V_cell,  dT/dt = P/(ρ·c_p·V);   ṅ(H₂) = I/(2F)' },
      { kind: 'formula', label: 'Constant-volume combustion bound', text: 'Δp ≈ (γ − 1)·q,  q = heat released per unit volume (oxygen- or fuel-limited)' },
    ],
    checks: [
      { name: 'Liquid below its flash point', pass: !above, detail: `${inp.Tsolv} °C vs flash point ${s.fp} °C` },
      { name: 'Unventilated release below 25 % of LEL', pass: fracLel <= 0.25, detail: `${pct(fracLel)} of LEL` },
      { name: 'Ventilated steady state below 25 % of LEL', pass: fracSS <= 0.25, detail: `${pct(fracSS)} of LEL at ${inp.achEnc} h⁻¹` },
      { name: 'Combustion equation balanced', pass: chem.auditEquation([s.f, 'O2'], ['CO2', 'H2O'], cb.coefficients).massBalanced, detail: `ν(O₂) = ${F(nuO2, 4)} per mole of solvent` },
      { name: 'Dust cloud below minimum explosible concentration', pass: !explosible, detail: comb ? `${F(cDust, 3)} vs ${inp.mec} g/m³` : 'powder not combustible' },
      { name: 'Voltage within extra-low-voltage touch limit', pass: inp.Uapp <= uLim, detail: `${F(inp.Uapp, 4)} V vs ${uLim} V ${inp.acdc.toUpperCase()}` },
      { name: 'Hydrogen below 25 % of its lower flammable limit', pass: h2ss < 1, detail: `${F(h2ss, 3)} vol % at steady state (limit 1 vol %)` },
    ],
    warnings: w,
    recommendations: [
      { priority: fracLel > 0.25 || above ? 'high' : 'medium', text: fracLel > 0.25 ? `Limit the solvent inventory in the enclosure to ${F(inp.Vsp * 0.25 / fracLel, 3)} mL or provide at least ${F(achNeed, 3)} air changes per hour with flow monitoring; exclude ignition sources (hot plates, non-rated motors, static) from the zone.` : above ? 'The vapour above the liquid is ignitable although the diluted concentration is low: keep containers closed, bond and earth metal vessels and use equipment rated for the hazardous zone.' : `The inventory is within 25 % of the LEL even if fully evaporated; keep it that way by fixing the maximum charge at ${F(inp.Vsp * 0.25 / Math.max(fracLel, 1e-9), 3)} mL.`, basis: 'Vapour concentration against the lower flammable limit.' },
      { priority: explosible ? 'high' : 'low', text: comb ? `Handle ${dm.label.toLowerCase()} under inert gas or as a wetted paste, earth all equipment, avoid dust deposits and commission dust-explosion testing (minimum ignition energy, minimum explosible concentration, maximum pressure) before scaling up.` : `${dm.label} does not burn; the remaining powder hazards are inhalation exposure and, for mixtures with organic binders or solvents, the flammability of those components.`, basis: 'Combustibility of the powder and dispersed concentration.' },
      { priority: inp.Uapp > uLim || tBoil < 600 ? 'high' : 'medium', text: `Fit the field-assisted reactor with a current limit near ${F(I * 1.2, 3)} A, a temperature interlock and ${inp.Uapp > uLim ? 'interlocked covers over all live parts' : 'insulated terminals'}; ${nH2 > 0 ? `vent the cell head space, which receives ${F(vH2 * 6e7, 3)} mL/min of hydrogen.` : 'check electrode connections for heating.'}`, basis: 'Ohmic current, Joule heating rate and touch-voltage limit.' },
      { priority: 'medium', text: `Among the listed solvents, ${Object.values(SOLV).sort((a, b) => b.fp - a.fp)[0].label.toLowerCase()} has the highest flash point; where the chemistry allows, a high-flash-point or aqueous medium removes the ignition hazard at source and reduces solvent emissions. Condense and recover evaporated solvent rather than diluting it to exhaust.`, basis: 'Inherently safer design: substitution before ventilation.', sustainability: true },
    ],
    publish: { 'risk.index': ratio, 'risk.lel_fraction': fracLel, 'proc.power': P, 'proc.E_field': E, 'elec.current_A': I },
    numerics: { method: 'Antoine vapour pressure, ideal-gas dilution balances, chem-kernel stoichiometry, Ohmic cell model; Brent root search for the LEL temperature', dof: 1, steps: 60, iterations: 1, residual: 0, notes: 'Well-mixed enclosure, ideal gas, constant properties. Flammability limits apply at ambient temperature and pressure in air; they widen at elevated temperature.' },
  };
}

// ---------- mode: FMEA, fault tree, risk matrix ----------
function parseTree(src) {
  const tok = String(src).replace(/\bAND\b/gi, '*').replace(/\bOR\b/gi, '+').match(/[A-Za-z_]\w*|[+*()]/g) || []; let i = 0;
  const fac = () => { const t = tok[i++]; if (t === '(') { const e = expr(); if (tok[i++] !== ')') throw new Error('missing closing parenthesis'); return e; } if (!t || /^[+*()]$/.test(t)) throw new Error(`unexpected "${t ?? 'end'}"`); return { ev: t.toLowerCase() }; };
  const term = () => { const a = [fac()]; while (tok[i] === '*') { i++; a.push(fac()); } return a.length > 1 ? { op: 'and', a } : a[0]; };
  const expr = () => { const a = [term()]; while (tok[i] === '+') { i++; a.push(term()); } return a.length > 1 ? { op: 'or', a } : a[0]; };
  const e = expr(); if (i < tok.length) throw new Error(`unexpected "${tok[i]}"`); return e;
}
const evalTree = (n, st) => n.ev !== undefined ? !!st[n.ev] : n.op === 'and' ? n.a.every(c => evalTree(c, st)) : n.a.some(c => evalTree(c, st));
const gateP = (n, p) => n.ev !== undefined ? (p[n.ev] ?? 0) : n.op === 'and' ? n.a.reduce((s, c) => s * gateP(c, p), 1) : 1 - n.a.reduce((s, c) => s * (1 - gateP(c, p)), 1);
const leaves = (n, out = []) => { if (n.ev !== undefined) out.push(n.ev); else n.a.forEach(c => leaves(c, out)); return out; };
function runRisk(inp) {
  const w = [], src = inp.fmeaSet !== 'table' && FMEA[inp.fmeaSet] ? FMEA[inp.fmeaSet] : inp.fmea, cl = v => N.clamp(Math.round(Number(v)) || 1, 1, 10);
  const fm = readRows(src, 4, 1, 'Failure-mode table').map(r => ({ name: r[0] || 'unnamed', S: cl(r[1]), O: cl(r[2]), D: cl(r[3]) })).map(f => ({ ...f, rpn: f.S * f.O * f.D })).sort((a, b) => b.rpn - a.rpn || b.S - a.S);
  if (Object.values(FMEA).includes(src)) w.push('Failure-mode table: illustrative default ratings — replace with the ratings agreed by your assessment team.');
  const tot = N.sum(fm.map(f => f.rpn)), act = fm.filter(f => f.rpn >= inp.rpnAct || f.S >= 9), nS9 = fm.filter(f => f.S >= 9).length, mit = fm.map((f, i) => i < 3 ? f.S * Math.max(1, f.O - 2) * Math.min(f.D, 2) : f.rpn), cum = []; fm.reduce((s, f) => { cum.push(100 * (s + f.rpn) / tot); return s + f.rpn; }, 0);
  const grid = new Float64Array(100); fm.forEach(f => { grid[(f.S - 1) * 10 + (f.O - 1)] += 1; });
  // fault tree
  let tree; try { tree = parseTree(inp.tree); } catch (e) { tree = parseTree(TREE_DEF); w.push(`Fault-tree expression could not be read (${e.message}); the default tree was used.`); }
  const p = {}; let evRows = []; try { evRows = readRows(inp.events, 2, 1, 'Basic-event table'); } catch (e) { w.push(e.message + ' All basic-event probabilities taken as zero.'); }
  evRows.forEach(r => { const v = Number(r[1]); if (Number.isFinite(v)) p[String(r[0]).toLowerCase()] = N.clamp(v, 0, 1); });
  if (inp.events === EV_DEF) w.push('Basic-event table: illustrative probabilities — replace with plant-specific or generic reliability data.');
  const lv = leaves(tree), names = [...new Set(lv)], miss = names.filter(nm => !(nm in p)); if (miss.length) { w.push(`No probability given for basic event(s) ${miss.join(', ')}: taken as zero.`); miss.forEach(nm => p[nm] = 0); }
  const repeated = lv.length !== names.length, nE = names.length; let pTop, exact = nE <= 16;
  const enumP = fix => { let s = 0; for (let m = 0; m < (1 << nE); m++) { let pr = 1; const st = {}; for (let j = 0; j < nE; j++) { const on = fix && fix[0] === names[j] ? fix[1] : !!(m & (1 << j)); if (fix && fix[0] === names[j]) { if (!!(m & (1 << j)) !== fix[1]) { pr = 0; break; } } else pr *= on ? p[names[j]] : 1 - p[names[j]]; st[names[j]] = on; } if (pr > 0 && evalTree(tree, st)) s += pr; } return s; };
  if (exact) pTop = enumP(null); else { pTop = gateP(tree, p); if (repeated) w.push('More than 16 basic events with repeated events: the gate-by-gate result assumes independence of branches and is approximate.'); }
  const imp = names.map(nm => { let hi, lo; if (exact) { hi = enumP([nm, true]); lo = enumP([nm, false]); } else { hi = gateP(tree, { ...p, [nm]: 1 }); lo = gateP(tree, { ...p, [nm]: 0 }); } return { nm, p: p[nm], birn: hi - lo, rrw: lo > 0 ? pTop / lo : (pTop > 0 ? 1e6 : 1), fv: pTop > 0 ? (pTop - lo) / pTop : 0 }; }).sort((a, b) => b.fv - a.fv);
  const pGate = gateP(tree, p), fYr = inp.nDem * pTop, pYr = 1 - Math.pow(1 - pTop, inp.nDem), Lb = fYr < 1e-4 ? 1 : fYr < 1e-3 ? 2 : fYr < 1e-2 ? 3 : fYr < 1e-1 ? 4 : 5, Cc = N.clamp(Math.round(inp.consClass), 1, 5), score = Lb * Cc, cls = score <= 4 ? 'low' : score <= 9 ? 'medium (reduce as far as reasonably practicable)' : score <= 16 ? 'high' : 'intolerable';
  const ks = N.logspace(-1, 1, 21), scaleP = k => { const q = {}; names.forEach(nm => q[nm] = Math.min(1, p[nm] * k)); return Math.max(gateP(tree, q), 1e-300); }, top = fm[0], tI = imp[0] || { nm: '–', fv: 0, birn: 0 };
  if (!exact || !repeated) { /* gate formula exact when no event is repeated */ }
  return {
    summary: `FMEA of ${fm.length} failure modes: the highest risk priority number is ${top.rpn} ("${top.name}"), ${act.length} mode(s) require action (RPN ≥ ${inp.rpnAct} or severity ≥ 9) and the first three modes carry ${pct((cum[Math.min(2, cum.length - 1)]) / 100)} of the total RPN. The fault tree with ${nE} basic events gives a top-event probability of ${F(pTop, 3)} per demand, i.e. ${F(fYr, 3)} expected events per year at ${F(inp.nDem, 3)} demands; ${tI.nm} is the most important basic event (${pct(tI.fv)} of the risk). Likelihood band ${Lb} with consequence class ${Cc} gives a risk-matrix score of ${score} of 25: ${cls}.`,
    kpis: [
      kp('rpn_max', 'Highest risk priority number', top.rpn, 'of 1000', `Severity ${top.S} × occurrence ${top.O} × detection ${top.D} of the top-ranked failure mode "${top.name}". RPN ranks modes for attention; it is an ordinal product, not a probability.`, top.rpn < inp.rpnAct),
      kp('n_act', 'Failure modes requiring action', act.length, `of ${fm.length}`, `Modes with RPN at or above ${inp.rpnAct} or with severity 9–10 (${nS9} such modes), which require action whatever their RPN.`, act.length === 0),
      kp('rpn_tot', 'Total RPN', tot, '', `Sum of all risk priority numbers; reducing occurrence by two points and bringing detection to 2 for the three top modes would lower it to ${F(N.sum(mit), 4)}.`, null),
      kp('p_top', 'Top-event probability per demand', pTop, '–', exact ? `Exact probability of the top event by enumeration of all ${1 << nE} basic-event states, assuming independent basic events (gate-by-gate formula: ${F(pGate, 3)}).` : 'Probability of the top event from the gate-by-gate formulae assuming independent branches.', pTop < 1e-4),
      kp('f_year', 'Expected top events per year', fYr, '1/yr', `Per-demand probability multiplied by ${F(inp.nDem, 3)} demands per year; the probability of at least one event in a year is ${F(pYr, 3)}.`, fYr < 1e-3),
      kp('fv_top', `Risk share of "${tI.nm}"`, tI.fv, '–', 'Fussell–Vesely importance: fraction of the top-event probability that disappears if this basic event is made impossible. The largest value marks the most effective place for a safeguard.', null),
      kp('score', 'Risk-matrix score', score, 'of 25', `Likelihood band ${Lb} (from the annual frequency) × consequence class ${Cc} on a 5 × 5 matrix: ${cls}.`, score <= 4 ? true : score > 9 ? false : null),
      kp('class', 'Risk class', cls, '', 'Qualitative class of the top event: low (≤ 4), medium (5–9), high (10–16), intolerable (> 16). Acceptance criteria are set by the organisation and regulator, not by this calculation.', score <= 9),
    ],
    series: [
      { key: 'pareto', title: 'Pareto chart of risk priority numbers', xlabel: 'failure mode (ranked)', ylabel: 'RPN and cumulative share (%)', lines: [{ name: 'RPN', x: fm.map((_, i) => i + 1), y: fm.map(f => f.rpn) }, { name: 'cumulative share (%)', x: fm.map((_, i) => i + 1), y: cum }, { name: 'action threshold', x: [1, Math.max(fm.length, 2)], y: [inp.rpnAct, inp.rpnAct] }], explain: 'Failure modes ranked by RPN with the cumulative share of the total. The modes to the left of the knee in the cumulative curve deserve the first corrective actions.' },
      { key: 'imp', title: 'Importance of basic events', xlabel: 'basic event (ranked)', ylabel: 'importance', lines: [{ name: 'Fussell–Vesely (share of risk)', x: imp.map((_, i) => i + 1), y: imp.map(e => e.fv) }, { name: 'Birnbaum (sensitivity ∂P_top/∂p)', x: imp.map((_, i) => i + 1), y: imp.map(e => e.birn) }], explain: 'Two importance measures per basic event, ranked as in the importance table. Fussell–Vesely shows where the present risk comes from; Birnbaum shows where the top event is most sensitive to a change in reliability.' },
      { key: 'sens', title: 'Top-event probability versus basic-event reliability', xlabel: 'multiplier on all basic-event probabilities', ylabel: 'top-event probability per demand', logx: true, logy: true, lines: [{ name: 'top event', x: ks, y: ks.map(scaleP) }], explain: 'Effect of making every basic event more or less likely by a common factor. A slope of two on these logarithmic axes shows that the top event needs two simultaneous failures: the benefit of redundant safeguards.' },
    ],
    maps: [{ key: 'matrix', title: 'FMEA severity–occurrence matrix', nx: 10, ny: 10, data: grid, x0: 0.5, x1: 10.5, y0: 0.5, y1: 10.5, xlabel: 'occurrence rating', ylabel: 'severity rating', unit: 'failure modes', explain: 'Number of failure modes in each severity–occurrence cell. Entries towards the upper right combine serious consequences with frequent occurrence and must be reduced first.' }],
    tables: [
      { key: 'fmea', title: 'Failure modes ranked by RPN', columns: ['Rank', 'Failure mode', 'S', 'O', 'D', 'RPN', 'Criticality S×O', 'Action required'], rows: fm.map((f, i) => [i + 1, f.name, f.S, f.O, f.D, f.rpn, f.S * f.O, f.rpn >= inp.rpnAct || f.S >= 9 ? 'yes' : 'no']), explain: 'Failure modes sorted by risk priority number with the severity–occurrence criticality. High-severity modes are flagged even when a low occurrence or good detection keeps their RPN small.' },
      { key: 'imp', title: 'Basic-event importance', columns: ['Rank', 'Basic event', 'Probability', 'Fussell–Vesely', 'Birnbaum', 'Risk-reduction worth'], rows: imp.map((e, i) => [i + 1, e.nm, F(e.p, 3), F(e.fv, 3), F(e.birn, 3), F(Math.min(e.rrw, 1e6), 3)]), explain: 'Importance measures of each basic event. Risk-reduction worth is the factor by which the top-event probability falls if the event is eliminated.' },
    ],
    equations: [
      { kind: 'formula', label: 'Risk priority number', text: 'RPN = S × O × D   (each rated 1–10)' },
      { kind: 'formula', label: 'Fault-tree gates (independent events)', text: 'P(AND) = Π pᵢ;   P(OR) = 1 − Π(1 − pᵢ)' },
      { kind: 'formula', label: 'Top event for this tree', text: 'TOP = ' + String(inp.tree).replace(/\*/g, ' · ').replace(/\s+/g, ' ').slice(0, 300), note: '· = AND, + = OR' },
      { kind: 'formula', label: 'Importance and annual frequency', text: 'I_B = P(top | e) − P(top | ē);  I_FV = [P(top) − P(top | ē)]/P(top);  f = n_demand·P(top);  P(≥1 per year) = 1 − (1 − P)ⁿ' },
    ],
    checks: [
      { name: 'No failure mode above the action threshold', pass: act.length === 0, detail: `${act.length} mode(s) with RPN ≥ ${inp.rpnAct} or S ≥ 9` },
      { name: 'Gate formula agrees with exact enumeration', pass: !exact || repeated || Math.abs(pGate - pTop) <= 1e-9 * Math.max(pTop, 1e-300) + 1e-15, detail: repeated ? `repeated events: gate formula ${F(pGate, 4)} vs exact ${F(pTop, 4)}` : `${F(pGate, 6)} vs ${F(pTop, 6)}` },
      { name: 'All basic events have probabilities', pass: miss.length === 0, detail: miss.length ? `missing: ${miss.join(', ')}` : `${nE} events defined` },
      { name: 'Probabilities within [0, 1]', pass: pTop >= 0 && pTop <= 1, detail: `P(top) = ${F(pTop, 4)}` },
      { name: 'Top-event risk not high or intolerable', pass: score <= 9, detail: `score ${score} of 25` },
    ],
    warnings: w,
    recommendations: [
      { priority: act.length ? 'high' : 'low', text: act.length ? `Act first on "${top.name}" (RPN ${top.rpn}): ${top.D >= top.O ? 'detection is the weakest rating, so add monitoring or an interlock that reveals the fault before it has consequences' : 'occurrence is the weakest rating, so remove the cause by design (inherently safer condition, poka-yoke, automation)'}; severity can only be lowered by changing the process itself.` : 'No failure mode exceeds the action threshold; review the FMEA after any change of recipe, equipment or scale.', basis: 'RPN ranking and its weakest factor.' },
      { priority: score > 9 ? 'high' : 'medium', text: `Strengthen the safeguard against "${tI.nm}": eliminating it removes ${pct(tI.fv)} of the top-event probability; adding an independent protection layer with a failure probability of 0.1 on that branch lowers the annual frequency to about ${F(fYr * (1 - 0.9 * tI.fv), 3)} per year.`, basis: 'Fussell–Vesely importance from the fault tree.' },
      { priority: 'medium', text: `Review common-cause failures: the calculation assumes independent basic events, yet shared utilities, one operator or one maintenance error can defeat several safeguards at once, which would raise the top-event probability above ${F(pTop, 3)}.`, basis: 'Independence assumption of the fault-tree quantification.' },
      { priority: 'medium', text: 'Include environmental and waste-related failure modes (untreated effluent, powder loss to exhaust, solvent release) in the same FMEA and give them severity ratings that reflect their environmental consequence; prevention at source avoids both the release and the lost material.', basis: 'Completeness of the hazard identification.', sustainability: true },
    ],
    publish: { 'risk.index': score / 9, 'risk.rpn_max': top.rpn, 'risk.p_top': pTop, 'risk.f_year': fYr, 'risk.matrix_score': score },
    numerics: { method: exact ? 'Exact fault-tree quantification by enumeration of basic-event states; RPN ranking' : 'Gate-by-gate fault-tree quantification; RPN ranking', dof: nE, steps: exact ? 1 << nE : 1, iterations: (exact ? (1 << nE) : 1) * (2 * nE + 1), residual: Math.abs(pGate - pTop), notes: 'Basic events are assumed statistically independent; probabilities are per demand.' },
  };
}

// ---------- mode: environmental fate ----------
function runEnv(inp) {
  const w = [], m = ENVM[inp.nm] || ENVM.ZnO, d = inp.dp, T = 293.15, mu = 1.0e-3, rw = 998, rel = inp.mProd * 1e3 * inp.wNP / 100 * inp.relFrac / 100, toW = rel * (1 - inp.wwtp / 100), toS = rel - toW;
  const pec = toW / (inp.Qriver * 3.156e7), pecUg = pec * 1e6, rq = pecUg / inp.pnec, mp = m.rho * Math.PI * d ** 3 / 6, n0 = pec / mp;
  const kB = 8 * N.K_B * T / (3 * mu), tHomo = 2 / (inp.alphaAtt * kB * Math.max(n0, 1e-300)), dS = 2e-6, rS = 2500, nS = inp.spm * 1e-3 / (rS * Math.PI * dS ** 3 / 6), b12 = 2 * N.K_B * T / (3 * mu) * (d + dS) ** 2 / (d * dS), kHet = inp.alphaAtt * b12 * nS, tHet = Math.log(2) / kHet;
  const vP = (m.rho - rw) * 9.80665 * d * d / (18 * mu), vS = (rS - rw) * 9.80665 * dS * dS / (18 * mu), tSetP = inp.depth / vP, tSetS = inp.depth / vS;
  // dissolution
  const wMet = chem.massFractions(m.f)[m.metal], Mmet = chem.ELEMENTS[m.metal].mass, cEq = m.logK !== null ? Math.pow(10, m.logK - 2 * inp.pHw) * Mmet * 1e3 : inp.csUser, pecMet = pecUg * wMet / 1e3, under = pecMet < cEq, Vm = chem.molarMass(m.f) / 1e3 / m.rho;
  const sat = N.clamp(pecMet / cEq, 0, 1), tSurf = (d / 2) / (inp.kDiss * Vm * Math.max(1 - sat, 1e-6)), tDiff = m.rho * (d / 2) ** 2 / (2 * 7e-10 * Math.max(cEq / wMet * 1e-3 * (1 - sat), 1e-300));
  const dissolves = m.ion && under;
  if (rq >= 1) w.push(`Risk quotient ${F(rq, 3)} ≥ 1: the predicted concentration reaches the no-effect concentration entered; refine the assessment and reduce release.`);
  if (m.logK === null && m.ion) w.push('Silver dissolves oxidatively and is rapidly sulfidised in wastewater; the equilibrium solubility entered by the user stands in for these processes, which are not modelled.');
  if (!m.ion) w.push('TiO₂ is treated as persistent particulate matter: no dissolution, fate governed by heteroaggregation and sedimentation.');
  if (dissolves) w.push(`The released ${m.label} is undersaturated (${F(pecMet, 3)} mg/L as ${m.metal} against a solubility of ${F(cEq, 3)} mg/L at pH ${inp.pHw}): it is expected to dissolve, so the aquatic hazard is that of the ${chem.prettyFormula(m.ion)} ion and the no-effect concentration must refer to the dissolved metal.`);
  w.push('The no-effect concentration is a user input and release fractions are scenario assumptions; the risk quotient is a screening indicator, not an environmental risk assessment.');
  // effluent neutralisation
  const wt = WASTE[inp.wasteType] || WASTE.naoh_hcl, ion = chem.ionicEquations(wt.re.map(f => ({ f, phase: 'aq' })), [{ f: wt.pr[0], phase: 'aq' }, { f: 'H2O', phase: 'l' }]), co = ion.coefficients, nW = inp.wasteVol * inp.wasteConc, nN = nW * co[1] / co[0], mN = nN * chem.molarMass(wt.re[1]) / 1e3, mSalt = nW * co[2] / co[0] * chem.molarMass(wt.pr[0]) / 1e3, nH2O = nW * co[3] / co[0], dTn = nH2O * 57.3e3 / (inp.wasteVol * 4180), tds = mSalt * 1e3 / inp.wasteVol;
  const aud = chem.auditEquation(wt.re, wt.pr, co), mIn = nW * chem.molarMass(wt.re[0]) / 1e3 + mN, mOut = mSalt + nH2O * chem.molarMass('H2O') / 1e3;
  if (dTn > 20) w.push(`Neutralisation heats the effluent by about ${F(dTn, 3)} K (adiabatic): dose the neutraliser slowly with stirring and cooling.`);
  const Qs = N.logspace(Math.log10(inp.Qriver / 100), Math.log10(inp.Qriver * 100), 40), th = N.linspace(0, Math.min(8 * tHet, 3.156e7), 60), pHs = N.linspace(4, 10, 49), tr = N.linspace(0, Math.min(tSurf, 3.156e9), 40);
  const eqs = [];
  if (m.re) { const c = chem.balance(m.re, m.pr); eqs.push({ kind: 'ionic', label: `Dissolution of ${m.label}`, text: chem.equationString(m.re.map((f, i) => f + (i === 0 ? '(s)' : f === 'O2' ? '(aq)' : '(aq)')), m.pr.map(f => f + (f === 'H2O' ? '(l)' : '(aq)')), c.coefficients), note: m.logK !== null ? `log K = ${m.logK} (25 °C, handbook value): [${m.metal}²⁺] = 10^(log K − 2·pH)` : 'Oxidative dissolution; kinetically controlled.' }); }
  else eqs.push({ kind: 'formula', label: 'Material', text: `${chem.prettyFormula(m.f)}(s) — no significant dissolution at pH 4–10`, note: 'Persistent particulate form.' });
  eqs.push({ kind: 'chemical', label: 'Effluent neutralisation', text: ion.molecular }, { kind: 'ionic', label: 'Complete ionic equation', text: ion.completeIonic }, { kind: 'ionic', label: 'Net ionic equation', text: ion.netIonic, note: 'Spectator ions: ' + [...new Set(ion.spectators)].join(', ') },
    { kind: 'formula', label: 'Release and predicted environmental concentration', text: 'ṁ_rel = m_product·w_NP·f_release;   PEC = ṁ_rel·(1 − η_WWTP)/Q_river;   RQ = PEC/PNEC' },
    { kind: 'governing', label: 'Aggregation (Smoluchowski, Brownian)', text: 'homo: dN/dt = −α·(4k_BT/3μ)·N²;   hetero: dN/dt = −α·β₁₂·N_SPM·N,  β₁₂ = (2k_BT/3μ)·(d₁ + d₂)²/(d₁d₂)' },
    { kind: 'formula', label: 'Sedimentation and dissolution', text: 'v_s = (ρ_p − ρ_w)·g·d²/(18μ);   surface control: dr/dt = −k_s·V_m·(1 − c/c_eq);   diffusion limit: t = ρ·r₀²/(2·D·c_s)' });
  return {
    summary: `Use and wear release ${F(rel, 3)} kg/yr of ${m.label}; after ${inp.wwtp} % retention in sludge ${F(toW, 3)} kg/yr reach a river of ${F(inp.Qriver, 3)} m³/s, giving a predicted environmental concentration of ${F(pecUg, 3)} µg/L and a risk quotient of ${F(rq, 3)} against the entered no-effect concentration. In the water the particles attach to natural suspended matter with a half-time of ${F(tHet / 3600, 3)} h (α = ${inp.alphaAtt}) and then settle within ${F(tSetS / 86400, 3)} days; ${m.ion ? `equilibrium solubility is ${F(cEq, 3)} mg/L as ${m.metal}, so the release ${under ? 'dissolves' : 'remains particulate'}` : 'the material does not dissolve'}. Neutralising ${F(inp.wasteVol, 3)} L of ${inp.wasteConc} mol/L effluent needs ${F(mN, 3)} kg of ${chem.prettyFormula(wt.re[1])} and forms ${F(mSalt, 3)} kg of ${chem.prettyFormula(wt.pr[0])}.`,
    kpis: [
      kp('release', 'Nanomaterial released per year', rel, 'kg/yr', `Product mass × nanomaterial content × annual release fraction; ${F(toS, 3)} kg/yr of it are retained in sewage sludge.`, null),
      kp('to_water', 'Emission to surface water', toW, 'kg/yr', `Release remaining after ${inp.wwtp} % removal in wastewater treatment; the sludge pathway (soil, incineration) carries the rest and is not followed further here.`, null),
      kp('pec', 'Predicted environmental concentration', pecUg, 'µg/L', `Annual emission to water diluted in the mean river flow (complete mixing, no upstream load), equal to ${F(n0 / 1e6, 3)} primary particles per cm³.`, rq < 1),
      kp('rq', 'Risk quotient PEC/PNEC', rq, '–', `Predicted concentration divided by the no-effect concentration of ${inp.pnec} µg/L entered. Below 1 no risk is indicated at screening level; at or above 1 refinement or risk reduction is required.`, rq < 0.1 ? true : rq >= 1 ? false : null),
      kp('t_het', 'Heteroaggregation half-time', tHet, 's', `Half-time for attachment to natural suspended particles (2 µm, ${inp.spm} mg/L) with attachment efficiency ${inp.alphaAtt}; homoaggregation at the predicted concentration would take ${F(tHomo / 3.156e7, 3)} years and is irrelevant.`, null),
      kp('v_set', 'Settling velocity of primary particles', vP, 'm/s', `Stokes velocity of a free ${F(d * 1e9, 3)} nm particle: ${F(tSetP / 3.156e7, 3)} years to settle ${inp.depth} m, so free nanoparticles do not settle; after attachment the carrier particles settle in ${F(tSetS / 86400, 3)} days.`, null),
      kp('c_eq', 'Equilibrium solubility', cEq, `mg/L as ${m.metal}`, m.logK !== null ? `Dissolved ${m.metal} in equilibrium with ${m.label} at pH ${inp.pHw} from log K = ${m.logK} (activity corrections and complexation neglected); it falls a hundred-fold per pH unit.` : 'User-entered solubility limit for a material without a built-in solubility constant.', null),
      kp('t_diss', 'Dissolution time (surface-controlled)', m.ion ? tSurf : 0, 's', m.ion ? `Time for complete dissolution of one particle with the entered surface rate constant, ${F(tSurf / 3600, 3)} h; the diffusion-limited lower bound is ${F(tDiff, 3)} s.` : 'The material is treated as insoluble; no dissolution time applies and zero is shown.', null),
      kp('m_neut', 'Neutraliser demand', mN, 'kg', `Mass of ${chem.prettyFormula(wt.re[1])} (${F(nN, 4)} mol) required by the stoichiometry ${co[0]}:${co[1]} for ${F(nW, 4)} mol of ${chem.prettyFormula(wt.re[0])} in the effluent.`, null),
      kp('m_salt', 'Salt formed', mSalt, 'kg', `Mass of ${chem.prettyFormula(wt.pr[0])} produced, a dissolved-salt load of ${F(tds, 3)} g/L in the neutralised effluent that remains to be discharged or treated.`, tds < 10),
      kp('dT_neut', 'Adiabatic temperature rise on neutralisation', dTn, 'K', 'Heating of the effluent by the heat of neutralisation (57.3 kJ per mole of water formed, strong acid with strong base, dilute solution).', dTn < 20),
    ],
    series: [
      { key: 'dil', title: 'Predicted concentration versus river flow', xlabel: 'river flow (m³/s)', ylabel: 'concentration (µg/L)', logx: true, logy: true, lines: [{ name: 'PEC', x: Qs, y: Qs.map(q => toW / (q * 3.156e7) * 1e6) }, { name: 'PNEC (user)', x: [Qs[0], Qs.at(-1)], y: [inp.pnec, inp.pnec] }, { name: 'this scenario', x: [inp.Qriver], y: [pecUg] }], explain: 'Predicted environmental concentration against the diluting river flow. Small receiving waters reach the no-effect concentration at releases that are harmless in large rivers.' },
      { key: 'het', title: 'Free nanoparticles remaining in the water column', xlabel: 'time (s)', ylabel: 'fraction of free particles', lines: [{ name: 'heteroaggregation with suspended matter', x: th, y: th.map(t => Math.exp(-kHet * t)) }], explain: 'First-order loss of free nanoparticles by attachment to natural suspended particles. After attachment the fate follows that of the carrier particles, which settle to the sediment.' },
      { key: 'ph', title: 'Solubility versus pH', xlabel: 'pH', ylabel: `solubility (mg/L as ${m.metal})`, logy: true, lines: [{ name: 'equilibrium solubility', x: pHs, y: pHs.map(ph => m.logK !== null ? Math.max(Math.pow(10, m.logK - 2 * ph) * Mmet * 1e3, 1e-300) : inp.csUser) }, { name: 'released concentration', x: [4, 10], y: [Math.max(pecMet, 1e-300), Math.max(pecMet, 1e-300)] }], explain: 'Equilibrium dissolved-metal concentration against pH compared with the released concentration. Where the release lies below the curve the particles dissolve and the ionic form governs toxicity.' },
      { key: 'r', title: 'Particle radius during dissolution', xlabel: 'time (s)', ylabel: 'radius (nm)', lines: [{ name: 'surface-controlled shrinkage', x: tr, y: tr.map(t => m.ion ? Math.max(0, (d / 2 - inp.kDiss * Vm * (1 - sat) * t)) * 1e9 : d / 2 * 1e9) }], explain: 'Shrinking of a single particle under surface-reaction control far from saturation. Small particles disappear proportionally faster, shifting the size distribution during ageing.' },
    ],
    tables: [{ key: 'bal', title: 'Release and effluent balance', columns: ['Quantity', 'Value', 'Unit'], rows: [['Nanomaterial in products', F(inp.mProd * 1e3 * inp.wNP / 100, 4), 'kg/yr'], ['Released during use', F(rel, 4), 'kg/yr'], ['Retained in sewage sludge', F(toS, 4), 'kg/yr'], ['Emitted to surface water', F(toW, 4), 'kg/yr'], [`Effluent ${chem.prettyFormula(wt.re[0])}`, F(nW, 4), 'mol'], [`Neutraliser ${chem.prettyFormula(wt.re[1])}`, F(mN, 4), 'kg'], [`Salt ${chem.prettyFormula(wt.pr[0])}`, F(mSalt, 4), 'kg'], ['Water formed', F(nH2O * 0.018015, 4), 'kg']], explain: 'Mass flows of the release scenario and of the effluent neutralisation. Material retained in sludge is not eliminated: it is transferred to soil or to incineration residues.' }],
    equations: eqs,
    checks: [
      { name: 'Release mass balance', pass: Math.abs(rel - toW - toS) < 1e-9 * Math.max(rel, 1e-300), detail: `${F(rel, 4)} = ${F(toW, 4)} (water) + ${F(toS, 4)} (sludge) kg/yr` },
      { name: 'Risk quotient below 1', pass: rq < 1, detail: `PEC ${F(pecUg, 3)} µg/L vs PNEC ${inp.pnec} µg/L` },
      { name: 'Neutralisation equation balanced (mass and charge)', pass: aud.massBalanced && aud.chargeBalanced, detail: ion.molecular },
      { name: 'Neutralisation mass balance', pass: Math.abs(mIn - mOut) < 1e-6 * Math.max(mIn, 1e-300), detail: `reactants ${F(mIn, 6)} kg = products ${F(mOut, 6)} kg` },
      { name: 'Dissolved-salt load of the effluent below 10 g/L', pass: tds < 10, detail: `${F(tds, 3)} g/L` },
    ],
    warnings: w,
    recommendations: [
      { priority: rq >= 1 ? 'high' : rq > 0.1 ? 'medium' : 'low', text: rq > 0.1 ? `Reduce the release: bind the nanomaterial in the matrix (lower release fraction) or cut the content; the release fraction must fall from ${inp.relFrac} to ${F(inp.relFrac * 0.1 / rq, 3)} %/yr for a risk quotient of 0.1.` : `The screening risk quotient is ${F(rq, 3)}; document the scenario and the source of the no-effect concentration, and repeat the estimate for the smallest receiving water in the region.`, basis: 'PEC/PNEC risk characterisation.' },
      { priority: 'medium', text: dissolves ? `Because ${m.label} dissolves at pH ${inp.pHw}, assess the ${chem.prettyFormula(m.ion)} ion against the dissolved-metal quality standard and consider the background metal concentration of the river.` : 'Because the particles persist and attach to suspended matter, assess the sediment compartment and sludge-amended soil, where they accumulate, not only the water column.', basis: 'Solubility against released concentration and heteroaggregation half-time.' },
      { priority: 'medium', text: `Neutralise the effluent in a stirred, cooled tank with pH control (expected rise ${F(dTn, 3)} K) and check discharge limits for ${chem.prettyFormula(wt.pr[0])}; ${mSalt > 10 ? 'at this salt load consider recovering the base or acid by membrane processes instead of neutralising it.' : 'the salt load is modest.'}`, basis: 'Neutralisation stoichiometry and heat of neutralisation.' },
      { priority: 'medium', text: `Recycle wash water and mother liquor in the synthesis: each litre reused avoids ${F(mN / inp.wasteVol * 1e3, 3)} g of neutraliser and ${F(tds, 3)} g of salt discharge, and precipitating dissolved ${m.metal} as hydroxide before discharge returns the metal to the process.`, basis: 'Reagent and salt load per litre of effluent.', sustainability: true },
    ],
    publish: { 'risk.index': rq, 'env.pec_ug_L': pecUg, 'env.rq': rq, 'env.release_kg_yr': rel, 'env.salt_kg': mSalt },
    numerics: { method: 'Closed-form release/dilution balance, Smoluchowski aggregation rates, Stokes settling, solubility equilibrium and chem-kernel stoichiometry', dof: 1, steps: 60, iterations: 1, residual: Math.abs(mIn - mOut) / Math.max(mIn, 1e-300), notes: 'Steady, fully mixed river at 20 °C; spherical particles; activity coefficients of 1; suspended matter represented by 2 µm particles of density 2500 kg/m³.' },
  };
}

// ---------- mode: device integration ----------
function runDevice(inp) {
  const w = [], dT = inp.Tuse - inp.Tdep, Mf = inp.Ef / (1 - inp.nuf), sig = Mf * (inp.ctes - inp.ctef) * dT, kap = 6 * sig * inp.tf * (1 - inp.nus) / (inp.Es * inp.ts * inp.ts), bow = kap * inp.Dw * inp.Dw / 8, tens = sig > 0;
  const Z = tens ? 1.976 : 0.5, Ebar = inp.Ef / (1 - inp.nuf * inp.nuf), G = Z * sig * sig * inp.tf / Ebar, tCrit = sig !== 0 ? inp.Gc * Ebar / (Z * sig * sig) : 1, sfG = G > 0 ? inp.Gc / G : 1e6;
  if (inp.tf / inp.ts > 0.05) w.push(`Film thickness is ${pct(inp.tf / inp.ts)} of the substrate thickness: the Stoney thin-film formula overestimates curvature above about 5 %.`);
  if (Math.abs(bow) > inp.ts) w.push('Predicted bow exceeds the substrate thickness: the plate deforms nonlinearly and the Stoney estimate is only qualitative.');
  if (G > inp.Gc) w.push(`Energy-release rate ${F(G, 3)} J/m² exceeds the toughness ${inp.Gc} J/m²: ${tens ? 'channel cracking of the film' : 'buckle-driven delamination'} is predicted; the film must be thinner than ${F(tCrit * 1e9, 3)} nm or deposited closer to service temperature.`);
  if (inp.nuf < 0) w.push('Negative Poisson\'s ratio entered for the film: the biaxial modulus E/(1 − ν) is lower than for a conventional material, which reduces mismatch stress; confirm that the effective properties apply at film scale.');
  // thermal stack
  const A = inp.Adev * 1e-4, Rf = inp.tf / (inp.kf * A), Rs = inp.ts / (inp.ks * A), Rt = inp.tTim * 1e-6 / (inp.kTim * A), Rtot = Rf + Rs + Rt + inp.Rsa, Tj = inp.Tamb + inp.Pdev * Rtot, Pmax = Math.max(0, (inp.Tjmax - inp.Tamb) / Rtot), parts = [['Film', Rf], ['Substrate', Rs], ['Interface material', Rt], ['Heat sink to ambient', inp.Rsa]];
  if (Tj > inp.Tjmax) w.push(`Junction temperature ${F(Tj, 4)} °C exceeds the limit of ${inp.Tjmax} °C; the stack can carry only ${F(Pmax, 3)} W.`);
  // electrical
  const Rsh = 1 / (inp.sigf * inp.tf), Rtr = Rsh * inp.squares, Rc1 = inp.rhoc / (inp.Ac * 1e-2), Rel = Rtr + 2 * Rc1, Vd = inp.Idev * Rel, Pl = inp.Idev * inp.Idev * Rel, LT = Math.sqrt(inp.rhoc / Rsh) * 1e-2;
  if (Rel > inp.Rmax) w.push(`Series resistance ${F(Rel, 3)} Ω exceeds the budget of ${F(inp.Rmax, 3)} Ω (${Rtr > 2 * Rc1 ? 'trace' : 'contacts'} dominate).`);
  // moisture ingress: concentration at the protected (impermeable) face of a slab exposed on one side
  const L = inp.Lenc * 1e-6, cw = t => { let s = 0; for (let n = 0; n < 60; n++) { const k = (2 * n + 1), term = (n % 2 ? -1 : 1) / k * Math.exp(-inp.Dm * k * k * Math.PI * Math.PI * t / (4 * L * L)); s += term; if (Math.abs(term) < 1e-14) break; } return N.clamp(1 - 4 / Math.PI * s, 0, 1); }, tChar = L * L / inp.Dm, thr = inp.cThr / 100;
  let tThr = tChar; try { tThr = Math.exp(N.brent(x => cw(Math.exp(x)) - thr, Math.log(tChar * 1e-4), Math.log(tChar * 50))); } catch { /* keep characteristic time */ }
  const tLag = tChar / 6, yr = 3.156e7;
  if (tThr < yr) w.push(`Moisture reaches ${inp.cThr} % of saturation at the device after ${F(tThr / 86400, 3)} days: a polymer encapsulant of this thickness only delays ingress; an inorganic barrier layer is needed for multi-year life.`);
  const Tds = N.linspace(Math.min(inp.Tuse, inp.Tdep), Math.max(inp.Tuse, inp.Tdep) + 1, 40), tfs = N.logspace(Math.log10(inp.tf / 20), Math.log10(inp.tf * 20), 40), tms = N.logspace(Math.log10(tChar * 1e-3), Math.log10(tChar * 5), 60), cumT = [inp.Tamb]; [inp.Rsa, Rt, Rs, Rf].forEach(r => cumT.push(cumT.at(-1) + inp.Pdev * r));
  const ok = [G <= inp.Gc, Tj <= inp.Tjmax, Rel <= inp.Rmax, tThr >= yr], worst = Math.max(G / inp.Gc, (Tj - inp.Tamb) / Math.max(inp.Tjmax - inp.Tamb, 1e-9), Rel / inp.Rmax);
  return {
    summary: `Cooling from ${inp.Tdep} to ${inp.Tuse} °C puts the film under ${F(Math.abs(sig) / 1e6, 3)} MPa of ${tens ? 'tensile' : 'compressive'} thermal-mismatch stress and bows the ${F(inp.Dw * 1e3, 3)} mm substrate by ${F(Math.abs(bow) * 1e6, 3)} µm (Stoney); the energy-release rate is ${F(G, 3)} J/m² against a toughness of ${inp.Gc} J/m². The thermal stack of ${F(Rtot, 3)} K/W gives a junction temperature of ${F(Tj, 4)} °C at ${F(inp.Pdev, 3)} W (limit ${inp.Tjmax} °C). Sheet resistance is ${F(Rsh, 3)} Ω/sq and the series resistance ${F(Rel, 3)} Ω against a budget of ${F(inp.Rmax, 3)} Ω. Moisture reaches ${inp.cThr} % of saturation at the device after ${F(tThr / 86400, 3)} days. ${ok.filter(Boolean).length} of 4 integration budgets are met; none of this demonstrates compliance with a product standard.`,
    kpis: [
      kp('sigma', 'Film thermal-mismatch stress', sig, 'Pa', `Biaxial stress E_f/(1 − ν_f)·(α_s − α_f)·ΔT for ΔT = ${F(dT, 4)} K; positive is tensile (cracking risk), negative compressive (buckling and delamination risk).`, G <= inp.Gc),
      kp('bow', 'Substrate bow', bow, 'm', `Centre deflection κ·D²/8 over ${F(inp.Dw * 1e3, 3)} mm from the Stoney curvature (radius ${F(kap !== 0 ? Math.abs(1 / kap) : 1e12, 3)} m); ${tens ? 'the coated side becomes concave' : 'the coated side becomes convex'}. It affects lithography focus and bonding.`, Math.abs(bow) < 5e-5),
      kp('G', 'Energy-release rate', G, 'J/m²', `Driving force Z·σ²·t_f/Ē_f for ${tens ? 'channel cracking (Z = 1.976)' : 'delamination of the compressed film (Z = 0.5)'}; failure is predicted when it exceeds the toughness of ${inp.Gc} J/m².`, G <= inp.Gc),
      kp('t_crit', 'Critical film thickness', tCrit, 'm', 'Film thickness at which the energy-release rate equals the toughness at the present stress; thinner films store too little elastic energy to crack or delaminate.', inp.tf <= tCrit),
      kp('Rth', 'Total thermal resistance', Rtot, 'K/W', `One-dimensional series resistance of film, substrate, interface material and heat sink; the largest contribution is ${parts.slice().sort((a, b) => b[1] - a[1])[0][0].toLowerCase()} (${pct(parts.slice().sort((a, b) => b[1] - a[1])[0][1] / Rtot)}).`, null),
      kp('Tj', 'Junction temperature', Tj, '°C', `Ambient temperature plus power × total thermal resistance; the maximum power for the ${inp.Tjmax} °C limit is ${F(Pmax, 3)} W. Heat spreading is neglected, so small sources on thick substrates run cooler than predicted and thin poorly conducting layers hotter.`, Tj <= inp.Tjmax),
      kp('Rsheet', 'Sheet resistance', Rsh, 'Ω/sq', `1/(σ·t) of the conducting film; the trace of ${F(inp.squares, 3)} squares contributes ${F(Rtr, 3)} Ω.`, null),
      kp('Rel', 'Series electrical resistance', Rel, 'Ω', `Trace plus two contacts of ${F(Rc1, 3)} Ω each (ρ_c/A_c); at ${F(inp.Idev, 3)} A the voltage drop is ${F(Vd, 3)} V and the loss ${F(Pl, 3)} W. The transfer length of the contact is ${F(LT * 1e6, 3)} µm.`, Rel <= inp.Rmax),
      kp('t_moist', `Time to ${inp.cThr} % moisture saturation at the device`, tThr, 's', `Fickian diffusion through ${inp.Lenc} µm of encapsulant on an impermeable device surface (${F(tThr / yr, 3)} years); the diffusion time lag L²/6D is ${F(tLag / 86400, 3)} days.`, tThr >= yr),
      kp('budgets', 'Integration budgets met', ok.filter(Boolean).length, 'of 4', 'Count of the four screening criteria satisfied: film integrity, junction temperature, resistance budget and one-year moisture protection.', ok.every(Boolean)),
    ],
    series: [
      { key: 'sT', title: 'Film stress versus service temperature', xlabel: 'temperature (°C)', ylabel: 'film stress (MPa)', lines: [{ name: 'thermal-mismatch stress', x: Tds, y: Tds.map(t => Mf * (inp.ctes - inp.ctef) * (t - inp.Tdep) / 1e6) }, { name: 'zero stress', x: [Tds[0], Tds.at(-1)], y: [0, 0] }], explain: 'Biaxial film stress between service and deposition temperature. The stress vanishes at the deposition temperature and grows linearly with the temperature excursion; its sign follows the difference in expansion coefficients.' },
      { key: 'bt', title: 'Bow and energy-release rate versus film thickness', xlabel: 'film thickness (m)', ylabel: 'bow (µm) and G/Γ', logx: true, logy: true, lines: [{ name: '|bow| (µm)', x: tfs, y: tfs.map(t => Math.max(Math.abs(6 * sig * t * (1 - inp.nus) / (inp.Es * inp.ts * inp.ts) * inp.Dw * inp.Dw / 8) * 1e6, 1e-12)) }, { name: 'G/Γ', x: tfs, y: tfs.map(t => Math.max(Z * sig * sig * t / Ebar / inp.Gc, 1e-12)) }, { name: 'failure criterion G/Γ = 1', x: [tfs[0], tfs.at(-1)], y: [1, 1] }], explain: 'Both the wafer bow and the cracking driving force rise in proportion to film thickness at fixed stress. The crossing of G/Γ with 1 marks the critical thickness.' },
      { key: 'stack', title: 'Temperature through the thermal stack', xlabel: 'interface (0 ambient → 4 junction)', ylabel: 'temperature (°C)', lines: [{ name: 'temperature', x: [0, 1, 2, 3, 4], y: cumT }, { name: 'junction limit', x: [0, 4], y: [inp.Tjmax, inp.Tjmax] }], explain: 'Temperature at each interface from ambient through heat sink, interface material and substrate to the film. The steepest segment is the layer that limits heat removal.' },
      { key: 'moist', title: 'Moisture at the protected surface', xlabel: 'time (s)', ylabel: 'concentration / saturation', logx: true, lines: [{ name: 'moisture at the device surface', x: tms, y: tms.map(cw) }, { name: 'critical level', x: [tms[0], tms.at(-1)], y: [thr, thr] }], explain: 'Build-up of moisture at the device surface under the encapsulant after exposure to a humid environment. Diffusion delays but does not prevent saturation; only the time scale L²/D is at the designer\'s disposal.' },
    ],
    tables: [
      { key: 'rth', title: 'Thermal resistance budget', columns: ['Layer', 'Resistance (K/W)', 'Share (%)', 'Temperature rise (K)'], rows: parts.map(([nm, r]) => [nm, F(r, 4), F(100 * r / Rtot, 3), F(inp.Pdev * r, 3)]), explain: 'Contribution of each layer to the total thermal resistance and to the temperature rise at the stated power. The dominant layer is where a design change is most effective.' },
      { key: 'qual', title: 'Qualification test plan', columns: ['Test', 'Failure mechanisms addressed', 'Related model result'], rows: QUAL, explain: 'Generic qualification tests for a device containing nanomaterial or metamaterial layers, linked to the model quantity that anticipates each mechanism. Test conditions and durations must be taken from the product standard that applies.' },
      { key: 'std', title: 'Standards and regulation checklist', columns: ['Document', 'Subject', 'What it requires of this work', 'Status'], rows: STANDARDS.map(r => [...r, 'evidence required — not established by modelling']), explain: 'Standards and regulatory frameworks commonly relevant to nanomaterial products, described generically. Applicability and current editions must be checked for the specific product and market; a calculation cannot establish conformity.' },
    ],
    equations: [
      { kind: 'constitutive', label: 'Thermal-mismatch film stress', text: 'σ_f = [E_f/(1 − ν_f)]·(α_s − α_f)·(T − T_dep)' },
      { kind: 'formula', label: 'Stoney curvature and bow', text: 'κ = 6·σ_f·t_f·(1 − ν_s)/(E_s·t_s²);   δ = κ·D²/8' },
      { kind: 'formula', label: 'Film cracking / delamination criterion', text: 'G = Z·σ_f²·t_f·(1 − ν_f²)/E_f ≥ Γ   (Z = 1.976 channel crack, 0.5 delamination)' },
      { kind: 'governing', label: 'Series thermal resistance', text: 'T_j = T_amb + P·Σ(tᵢ/(kᵢ·A)) + P·R_sa' },
      { kind: 'formula', label: 'Sheet and contact resistance', text: 'R_□ = 1/(σ·t_f),  R = R_□·(L/W) + 2·ρ_c/A_c,  L_T = √(ρ_c/R_□)' },
      { kind: 'governing', label: 'Moisture diffusion through the encapsulant', text: '∂c/∂t = D·∂²c/∂x²;  c(L,t)/c_s = 1 − (4/π)·Σ [(−1)ⁿ/(2n+1)]·exp(−D(2n+1)²π²t/(4L²));  t_lag = L²/(6D)' },
    ],
    checks: [
      { name: 'Film integrity (G ≤ Γ)', pass: ok[0], detail: `G = ${F(G, 3)} J/m² vs Γ = ${inp.Gc} J/m²` },
      { name: 'Thin-film (Stoney) assumption valid', pass: inp.tf / inp.ts <= 0.05 && Math.abs(bow) <= inp.ts, detail: `t_f/t_s = ${F(inp.tf / inp.ts, 3)}, |bow|/t_s = ${F(Math.abs(bow) / inp.ts, 3)}` },
      { name: 'Junction temperature within limit', pass: ok[1], detail: `${F(Tj, 4)} °C vs ${inp.Tjmax} °C` },
      { name: 'Electrical resistance within budget', pass: ok[2], detail: `${F(Rel, 4)} Ω vs ${F(inp.Rmax, 4)} Ω` },
      { name: 'Moisture protection for at least one year', pass: ok[3], detail: `${F(tThr / yr, 3)} years to ${inp.cThr} % saturation` },
      { name: 'Thermal energy balance of the stack', pass: Math.abs(cumT.at(-1) - Tj) < 1e-9 * Math.max(Math.abs(Tj), 1), detail: `sum of layer temperature rises reproduces T_j = ${F(Tj, 5)} °C` },
    ],
    warnings: w,
    recommendations: [
      { priority: !ok[0] ? 'high' : 'medium', text: !ok[0] ? `Reduce film stress or thickness: keep the film below ${F(tCrit * 1e9, 3)} nm, lower the deposition temperature, or choose a substrate with expansion closer to ${F(inp.ctef * 1e6, 3)}×10⁻⁶ K⁻¹; a compliant or graded interlayer also lowers the driving force.` : `Film integrity has a margin of ${F(sfG, 3)} on energy-release rate; verify the toughness by an adhesion test, because interface toughness varies strongly with surface preparation.`, basis: 'Energy-release rate against film or interface toughness.' },
      { priority: !ok[1] ? 'high' : 'low', text: !ok[1] ? `Lower the thermal resistance of the ${parts.slice().sort((a, b) => b[1] - a[1])[0][0].toLowerCase()}, which carries ${pct(parts.slice().sort((a, b) => b[1] - a[1])[0][1] / Rtot)} of the total, or limit the power to ${F(Pmax, 3)} W.` : `Thermal margin is ${F(inp.Tjmax - Tj, 3)} K at ${F(inp.Pdev, 3)} W; the stack can dissipate up to ${F(Pmax, 3)} W.`, basis: 'Series thermal-resistance budget.' },
      { priority: !ok[3] ? 'high' : 'medium', text: !ok[3] ? `Add an inorganic barrier (oxide or nitride layer, glass lid) or increase the polymer thickness to ${F(Math.sqrt(yr / tThr) * inp.Lenc, 3)} µm for one year of protection; ingress time scales with thickness squared.` : 'Encapsulation delays moisture beyond one year at the stated diffusivity; confirm with a damp-heat test, since edge paths and interfaces usually leak faster than the bulk.', basis: 'Fickian ingress time through the encapsulant.' },
      { priority: 'medium', text: 'Design the stack for disassembly and material recovery: avoid permanent bonding of noble-metal or critical-element films to low-value substrates where a releasable interface meets the thermal and mechanical budgets, and document the layer materials in the product passport (Suite 28).', basis: 'End-of-life recoverability of the integrated layers.', sustainability: true },
      { priority: 'high', text: 'Treat the standards checklist as a list of evidence still to be produced: qualification testing, biological or environmental testing and conformity assessment by the responsible bodies are required before any compliance statement.', basis: 'Modelling alone never establishes regulatory compliance or safety.' },
    ],
    publish: { 'mech.max_stress': Math.abs(sig), 'mech.safety_factor': sfG, 'therm.T_max': Tj + 273.15, 'risk.index': worst, 'elec.resistance_ohm': Rel, 'life.years': tThr / yr },
    numerics: { method: 'Closed-form thin-film mechanics (Stoney), series thermal and electrical resistances, Fourier-series solution of Fickian diffusion with Brent root search', dof: 1, steps: 60, iterations: 60, residual: Math.abs(cumT.at(-1) - Tj), notes: 'Elastic, temperature-independent properties; one-dimensional heat flow without spreading; constant diffusivity and boundary humidity.' },
  };
}

export default {
  id: 27,
  tagline: 'Screens worker exposure, process hazards, environmental release and device-integration budgets of nanomaterial and metamaterial work, and lists the evidence a safety or compliance case still needs.',
  inputs,
  presets: [
    { name: 'ZnO nanorods: hydrothermal autoclave at 180 °C', note: 'Sealed 100 mL lined autoclave, 70 % fill: saturation pressure, margin to liquid-full, shell stress and stored energy.', values: { mode: 'autoclave', Tac: 180, fillPct: 70, Vliner: 100, IDmm: 50, wall: 8, yieldT: 130 } },
    { name: 'Ag nanoparticle powder: weighing and transfer', note: 'Dry silver nanopowder handled on a bench with partial capture; limit left to the user (enter the applicable national value).', values: { mode: 'exposure', G: 0.2, dp: 4e-8, rhoP: 10490, lev: 80, tTask: 30, limitSrc: 'user', oel: 0.01, sevBand: 'high', form: 'dry', mHandle: 2000 } },
    { name: 'Carbon nanotube handling against the NIOSH REL', note: 'Fibre-like nanomaterial with the 1 µg/m³ recommended limit: shows why containment rather than ventilation is required.', values: { mode: 'exposure', G: 0.02, dp: 1e-7, rhoP: 1800, lev: 95, tTask: 45, limitSrc: 'cnt', sevBand: 'vhigh', form: 'dry', mHandle: 500 } },
    { name: 'Au colloid synthesis: FMEA and fault tree', note: 'Citrate/borohydride colloid laboratory process with built-in failure-mode list.', values: { mode: 'risk', fmeaSet: 'colloid', nDem: 100, consClass: 3, tree: 'h2_buildup*(vent_blocked + cap_sealed)*ignition_source + acid_splash*ppe_missing', events: 'event,probability\nh2_buildup,0.05\nvent_blocked,0.02\ncap_sealed,0.03\nignition_source,0.01\nacid_splash,0.01\nppe_missing,0.05' } },
    { name: 'Lithographic metasurface on fused silica: integration', note: 'TiO₂ metasurface layer deposited at 250 °C on fused silica with a transparent conductor and polymer encapsulation.', values: { mode: 'device', Ef: 2.3e11, nuf: 0.27, ctef: 9e-6, tf: 6e-7, Es: 7.2e10, nus: 0.17, ctes: 5e-7, ts: 5e-4, Dw: 0.05, Tdep: 250, ks: 1.4, kf: 5, Pdev: 0.2, sigf: 5e5, Gc: 8 } },
    { name: '3-D printed auxetic lattice: isopropanol wash', note: 'Solvent-vapour flammability of the wash station for photopolymer lattices, with resin-free powder hazards switched to carbon filler.', values: { mode: 'hazards', solvent: 'isopropanol', Tsolv: 25, Vsp: 500, Venc: 0.3, achEnc: 20, dustMat: 'C', mDust: 5, Uapp: 24 } },
    { name: 'ZnO in a coating: environmental release', note: 'Wear release of nano-ZnO from an exterior coating to a medium-sized river, with neutralisation of alkaline mother liquor.', values: { mode: 'envfate', nm: 'ZnO', mProd: 200, wNP: 3, relFrac: 2, Qriver: 15, pnec: 1, pHw: 7.8 } },
  ],
  run(inp, ctx) {
    const r = ({ exposure: runExposure, aerosol: runAerosol, autoclave: runAutoclave, hazards: runHazards, risk: runRisk, envfate: runEnv, device: runDevice }[inp.mode] || runExposure)(inp, ctx);
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    r.warnings.push(NOTE);
    return r;
  },
};
