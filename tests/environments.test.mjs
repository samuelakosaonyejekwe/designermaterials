// Tests for the creation-environments knowledge base (data/environments.js).
// Usage: node -e "import('./tests/environments.test.mjs').then(m=>m.default()).then(r=>console.log(r))"
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { ENV_CONDITIONS, CONDITION_CATEGORIES, METAMATERIAL_ROUTES, recommendEnvironments } from '../data/environments.js';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

export default async function () {
  let pass = 0, fail = 0; const log = [];
  const check = (cond, msg) => { if (cond) pass++; else { fail++; log.push('✗ ' + msg); } };
  const isNum = v => typeof v === 'number' && Number.isFinite(v);
  const isStr = v => typeof v === 'string' && v.trim().length > 0;
  const strList = (a, min) => Array.isArray(a) && a.length >= min && a.every(isStr);

  const spec = JSON.parse(fs.readFileSync(path.join(root, 'data/spec.json'), 'utf8'));
  const envList = spec.environments;
  const KEYS = ['T_C', 'P_Pa', 'time', 'atmosphere', 'medium', 'energy', 'scale', 'chemical', 'thermal', 'transport', 'fields', 'mechanical', 'controls', 'outputsControlled', 'suites', 'hazards', 'sustainability', 'resolution_nm', 'rate'];

  // 1. Coverage and shape of every entry
  check(Object.keys(ENV_CONDITIONS).length === 146, `expected 146 entries, found ${Object.keys(ENV_CONDITIONS).length}`);
  check(envList.length === 146, 'spec.json should list 146 environments');
  for (let id = 1; id <= 146; id++) {
    const c = ENV_CONDITIONS[id];
    if (!c) { check(false, `environment ${id} missing`); continue; }
    for (const k of KEYS) check(k in c && c[k] !== undefined, `[${id}] key ${k} missing`);
    for (const k of ['T_C', 'P_Pa', 'time']) {
      const r = c[k];
      check(Array.isArray(r) && r.length === 2 && isNum(r[0]) && isNum(r[1]) && r[0] <= r[1], `[${id}] ${k} must be an ordered numeric pair: ${JSON.stringify(r)}`);
    }
    check(c.T_C[0] >= -273.15 && c.T_C[1] <= 25000, `[${id}] temperature outside physical range`);
    check(c.P_Pa[0] >= 0 && c.P_Pa[1] <= 1e12, `[${id}] pressure outside physical range`);
    check(c.time[0] > 0, `[${id}] duration must be positive`);
    for (const k of ['atmosphere', 'medium', 'energy', 'scale', 'sustainability']) check(isStr(c[k]), `[${id}] ${k} must be non-empty text`);
    for (const k of ['chemical', 'thermal', 'transport', 'mechanical', 'hazards']) check(strList(c[k], 1), `[${id}] ${k} must be a non-empty list of strings`);
    check(Array.isArray(c.fields) && c.fields.every(isStr), `[${id}] fields must be a list of strings (may be empty)`);
    check(strList(c.controls, 3), `[${id}] needs >= 3 controls`);
    check(strList(c.outputsControlled, 3), `[${id}] needs >= 3 controlled outputs`);
    check(Array.isArray(c.suites) && c.suites.length >= 3 && c.suites.every(s => Number.isInteger(s) && s >= 1 && s <= 29) && new Set(c.suites).size === c.suites.length, `[${id}] suites must be unique ids within 1–29`);
    const r = c.resolution_nm;
    check(r === null || (Array.isArray(r) && r.length === 2 && isNum(r[0]) && isNum(r[1]) && r[0] > 0 && r[0] <= r[1]), `[${id}] resolution_nm must be null or an ordered positive pair`);
    check(c.rate === null || isStr(c.rate), `[${id}] rate must be text or null`);
    check(c.sustainability.length <= 260, `[${id}] sustainability text too long for mobile`);
  }

  // 2. Spot checks against well-known process windows
  const E = ENV_CONDITIONS, within = (r, v) => r[0] <= v && v <= r[1];
  check(within(E[1].T_C, 90) && within(E[1].T_C, 200), 'hydrothermal window should include 90 and 200 °C');
  check(within(E[27].T_C, 1000), 'CVD window should include graphene growth at ~1000 °C');
  check(E[29].T_C[0] <= 100 && E[29].T_C[1] <= 500 && within(E[29].P_Pa, 100), 'ALD should be 50–400 °C at ~10–1000 Pa');
  check(E[33].P_Pa[0] >= 0.01 && E[33].P_Pa[1] <= 10, 'magnetron sputtering works at ~0.1–5 Pa');
  check(E[34].P_Pa[0] <= 1e-7 && E[34].P_Pa[1] <= 1e-3, 'MBE is an ultra-high-vacuum process');
  check(E[36].T_C[0] >= 1000 && E[36].T_C[1] >= 2000, 'flame synthesis is 1500–2500 °C');
  check(E[44].T_C[1] >= 1800, 'spark plasma sintering reaches ~2000 °C');
  check(E[63].resolution_nm && E[63].resolution_nm[0] <= 10 && E[63].resolution_nm[1] <= 30, 'e-beam lithography resolves 5–20 nm');
  check(E[67].resolution_nm && E[67].resolution_nm[0] >= 50 && E[67].resolution_nm[1] <= 500, 'two-photon polymerisation voxels are 100–300 nm');
  check(E[83].T_C[1] >= 1700 && /Ar|N₂/.test(E[83].atmosphere), 'LPBF melt pool exceeds alloy melting point under inert gas');
  check(E[126].P_Pa[0] >= 4e9 && E[126].P_Pa[1] <= 1e10 && within(E[126].T_C, 1500), 'HPHT is ~5–6 GPa at 1300–1600 °C');
  check(E[127].P_Pa[1] <= 1e-6, 'UHV must be below 1e-6 Pa');
  check(E[137].P_Pa[1] >= 1e11, 'diamond-anvil cells exceed 100 GPa');
  check(E[132].P_Pa[0] >= 7.3e6 && E[132].T_C[0] >= 31, 'supercritical CO₂ is above 31 °C and 7.38 MPa');
  check(E[141].P_Pa[0] >= 1e7, 'deep-sea vents sit under tens of MPa');
  for (const id of [61, 62, 64, 66, 67, 68, 86, 87, 88]) check(E[id].P_Pa[0] >= 5e4, `[${id}] should be an ambient-pressure environment`);
  for (const id of [30, 31, 32, 34, 60, 63]) check(E[id].P_Pa[1] < 100, `[${id}] should be a vacuum environment`);
  for (const id of [54, 69]) check(E[id].resolution_nm === null && /mask/.test(E[id].rate), `[${id}] etch tools transfer patterns: resolution is set by the mask`);
  for (const id of [60, 62, 63, 64, 65, 66, 67, 68, 83, 86, 87, 88]) check(Array.isArray(E[id].resolution_nm), `[${id}] patterning/printing environment needs resolution_nm`);
  for (const id of [1, 4, 6, 26, 36, 139, 146]) check(E[id].resolution_nm === null, `[${id}] non-patterning environment should have resolution_nm null`);

  // 3. Condition categories (specification areas 30–34)
  check(Array.isArray(CONDITION_CATEGORIES) && CONDITION_CATEGORIES.length === 5, 'five condition categories expected');
  check(CONDITION_CATEGORIES.map(c => c.key).join() === 'chemical,thermal,transport,fields,mechanical', 'condition category keys/order');
  const expectN = { chemical: 24, thermal: 25, transport: 26, fields: 24, mechanical: 26 };
  for (const c of CONDITION_CATEGORIES) {
    check(isStr(c.title) && strList(c.parameters, 20), `category ${c.key} needs a title and parameter list`);
    check(c.parameters.length === expectN[c.key], `category ${c.key}: expected ${expectN[c.key]} parameters, found ${c.parameters.length}`);
    check(new Set(c.parameters).size === c.parameters.length, `category ${c.key} has duplicate parameters`);
    check(c.key in E[1], `category key ${c.key} must exist on environment entries`);
  }
  const allParams = CONDITION_CATEGORIES.flatMap(c => c.parameters);
  for (const p of ['Supersaturation', 'pH', 'Vacuum levels', 'Reynolds number', 'Polarisation', 'Cleanroom contamination']) check(allParams.includes(p), `parameter "${p}" missing from categories`);

  // 4. Metamaterial routes (table J)
  const classes = (spec.metaEnv || []).map(m => m.cls);
  check(classes.length === 15, 'spec.json metaEnv should list 15 classes');
  check(Array.isArray(METAMATERIAL_ROUTES) && METAMATERIAL_ROUTES.length === 15, '15 metamaterial routes expected');
  for (const cls of classes) check(METAMATERIAL_ROUTES.some(r => r.cls === cls), `route missing for class "${cls}"`);
  for (const r of METAMATERIAL_ROUTES) {
    check(Array.isArray(r.envIds) && r.envIds.length >= 3 && r.envIds.every(i => Number.isInteger(i) && i >= 1 && i <= 138), `[${r.cls}] envIds must be engineered environment ids`);
    check(new Set(r.envIds).size === r.envIds.length, `[${r.cls}] duplicate envIds`);
    check(isStr(r.note), `[${r.cls}] note missing`);
    const cited = (r.note.match(/\b\d{2,3}\b/g) || []).map(Number).filter(n => n >= 1 && n <= 146 && E[n]);
    check(cited.every(n => r.envIds.includes(n) || n === 146), `[${r.cls}] note cites ids not in envIds`);
  }
  const route = c => METAMATERIAL_ROUTES.find(r => r.cls === c)?.envIds || [];
  check(route('Optical metamaterials').includes(63) && route('Optical metamaterials').includes(64), 'optical metamaterials need e-beam and nanoimprint lithography');
  check(route('Mechanical metamaterials').includes(67), 'mechanical metamaterials need two-photon polymerisation');
  check(route('Auxetic metamaterials').includes(83) && route('Auxetic metamaterials').includes(88), 'auxetics need metal AM and polymer printing');
  check(route('Metasurfaces').includes(64) && route('Metasurfaces').includes(69), 'metasurfaces need nanoimprint and etching');

  // 5. Recommender
  const top = (q, n, withList = true) => recommendEnvironments(q, withList ? envList : undefined).slice(0, n).map(r => r.id);
  const shape = recommendEnvironments({ kind: 'nanomaterial', dimensionality: '1D', materialFamily: 'oxide', scale: 'lab' }, envList);
  check(shape.length === 12, 'recommender returns the top 12');
  check(shape.every(r => Number.isInteger(r.id) && isNum(r.score) && Array.isArray(r.why) && r.why.length > 0 && r.why.every(isStr)), 'each hit needs id, finite score and explanatory why lines');
  check(shape.every((r, i) => i === 0 || shape[i - 1].score >= r.score), 'hits sorted by descending score');
  check(JSON.stringify(shape) === JSON.stringify(recommendEnvironments({ kind: 'nanomaterial', dimensionality: '1D', materialFamily: 'oxide', scale: 'lab' }, envList)), 'recommender is deterministic');
  for (const r of shape) {
    const sum = r.why.reduce((s, w) => s + Number(w.replace('−', '-').split(' ')[0]), 0);
    check(Math.abs(sum - r.score) < 1e-9, `[${r.id}] why lines must add up to the score (${sum} vs ${r.score})`);
  }

  const rods = top({ kind: 'nanomaterial', dimensionality: '1D', materialFamily: 'oxide', scale: 'lab' }, 8);
  log.push('oxide 1D nanorods (lab): ' + rods.join(', '));
  for (const id of [1, 15, 46]) check(rods.includes(id), `oxide 1-D nanorods at lab scale should rank environment ${id} in the top 8 (got ${rods.join(', ')})`);
  check(rods[0] === 1, `hydrothermal should lead for oxide nanorods (got ${rods[0]})`);
  const rodsNoList = top({ kind: 'nanomaterial', dimensionality: '1D', materialFamily: 'oxide', scale: 'lab' }, 8, false);
  check(rodsNoList.includes(1) && rodsNoList.includes(46), 'recommender works without envList');
  const rodsCheap = top({ kind: 'nanomaterial', dimensionality: '1D', materialFamily: 'oxide', scale: 'lab', T_max_C: 120, needs: ['low-cost', 'vacuum-free', 'alignment'] }, 5);
  log.push('oxide nanorods, ≤120 °C, low-cost, vacuum-free, aligned: ' + rodsCheap.join(', '));
  for (const id of [1, 15, 46]) check(rodsCheap.includes(id), `low-temperature aligned oxide nanorods should rank ${id} in the top 5 (got ${rodsCheap.join(', ')})`);

  const mech = top({ kind: 'metamaterial', dimensionality: '3D', materialFamily: 'polymer', featureSize_nm: 500, needs: ['3D-architecture'] }, 3);
  log.push('3-D polymer mechanical metamaterial, 500 nm: ' + mech.join(', '));
  check(mech[0] === 67, `two-photon polymerisation should lead for 500 nm polymer lattices (got ${mech.join(', ')})`);
  check(!mech.includes(88) && !mech.includes(86), 'filament and vat printers cannot resolve 500 nm and must not be in the top 3');

  const meta = top({ kind: 'metamaterial', dimensionality: '2D', featureSize_nm: 50, scale: 'lab', needs: ['periodic'] }, 4);
  log.push('optical metasurface prototype, 50 nm, lab: ' + meta.join(', '));
  check(meta.includes(63) && meta.includes(64), `optical metasurface with 50 nm features should rank e-beam and nanoimprint lithography in the top 4 (got ${meta.join(', ')})`);
  check(!meta.includes(68) && !meta.includes(91), 'direct laser writing and inkjet cannot resolve 50 nm');

  const aux = top({ kind: 'metamaterial', dimensionality: '3D', materialFamily: 'metal', featureSize_nm: 500000, scale: 'industrial', needs: ['3D-architecture'] }, 3);
  log.push('metal auxetic lattice, 0.5 mm struts, industrial: ' + aux.join(', '));
  check(aux.includes(83) && aux.includes(84), `metal lattices with 0.5 mm struts should rank powder-bed fusion (got ${aux.join(', ')})`);

  const gr = top({ kind: 'nanomaterial', dimensionality: '2D', materialFamily: 'carbon', needs: ['high-purity'] }, 3);
  log.push('2-D carbon, high purity: ' + gr.join(', '));
  check(gr.includes(27), `high-purity 2-D carbon should rank CVD (got ${gr.join(', ')})`);

  const epi = top({ kind: 'nanomaterial', dimensionality: '2D', materialFamily: 'semiconductor', needs: ['epitaxial', 'high-purity'] }, 3);
  log.push('epitaxial semiconductor films: ' + epi.join(', '));
  check(epi.includes(34), `epitaxial semiconductor films should rank MBE (got ${epi.join(', ')})`);

  const cold = recommendEnvironments({ kind: 'nanomaterial', dimensionality: '2D', materialFamily: 'oxide', T_max_C: 150, needs: ['conformal'] }, envList);
  log.push('conformal oxide film ≤150 °C: ' + cold.slice(0, 4).map(r => r.id).join(', '));
  check(cold.slice(0, 4).some(r => r.id === 29), 'conformal oxide films under 150 °C should rank ALD in the top 4');
  check(cold.every(r => E[r.id].T_C[0] <= 150), 'no top hit may need more than the thermal budget');

  for (const q of [{}, undefined, { kind: 'x', dimensionality: '9D', materialFamily: '?', featureSize_nm: -1, T_max_C: NaN, scale: 'z', needs: 'bad' }, { needs: [null, 3, 'vacuum-free'] }]) {
    let ok = true, r;
    try { r = recommendEnvironments(q, envList); ok = r.length === 12 && r.every(x => isNum(x.score)); } catch { ok = false; }
    check(ok, `recommender must tolerate input ${JSON.stringify(q)}`);
  }
  const any = recommendEnvironments({}, envList);
  check(any.every(r => r.id < 139), 'natural environments are never recommended by default');

  // 6. Size budget and house rules
  const src = fs.readFileSync(path.join(root, 'data/environments.js'), 'utf8');
  const bytes = Buffer.byteLength(src);
  log.push(`environments.js size: ${(bytes / 1024).toFixed(1)} kB`);
  check(bytes < 160 * 1024, `environments.js must stay under 160 kB (is ${(bytes / 1024).toFixed(1)} kB)`);
  check(!/^\s*import\s|\bimport\(|require\(|\bwindow\.\w|\bdocument\.\w|\bfetch\(|XMLHttpRequest|await /m.test(src), 'environments.js must be a pure module (no imports, DOM or network)');

  return { pass, fail, log };
}
