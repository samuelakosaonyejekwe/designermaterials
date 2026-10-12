// Builds data/spec.json from the engineering specification documents (converted to Markdown).
// Usage: node tools/build-spec.mjs <dir-with-md-files>
// Only structured engineering content (lists, tables, named models) is extracted.
import fs from 'node:fs';
import path from 'node:path';

const src = process.argv[2];
if (!src) { console.error('usage: node tools/build-spec.mjs <spec-md-dir>'); process.exit(1); }
const read = n => fs.readFileSync(path.join(src, n + '.md'), 'utf8').split(/\r?\n/);
const clean = s => s.replace(/\\([#.|*])/g, '$1').replace(/\*\*/g, '').trim();
const unbold = s => { const m = s.match(/^\*\*(.+)\*\*$/); return m ? clean(m[1]) : null; };

const SUITES = [
  ['Material Requirements and Functional Design', 'Target properties, functionality, specifications and engineering constraints'],
  ['Chemical Composition and Material Selection', 'Elements, compounds, dopants, constituent selection and property databases'],
  ['Atomic, Molecular and Crystal Structure', 'Atomic geometry, molecular structure, crystallography and crystal defects'],
  ['Nanomaterial Geometry and Morphology', 'Nanoparticles, nanorods, fibres, tubes, sheets, films, pores and coatings'],
  ['Metamaterial Architecture and Unit-Cell Design', 'Periodic structures, resonators, lattices, topology, metasurfaces and graded architectures'],
  ['Thermodynamics, Chemical Equilibrium and Phase Behaviour', 'Chemical potential, solubility, stability and phase transformation'],
  ['Chemical Reaction, Nucleation and Growth', 'Reaction kinetics, nucleation, crystal growth and particle population evolution'],
  ['Surface, Interface and Colloidal Engineering', 'Surface chemistry, adhesion, adsorption, wetting and particle interactions'],
  ['Fluid Dynamics and Porous-Media Transport', 'Flow, pressure, permeability, mixing and fluid–particle behaviour'],
  ['Heat Transfer and Thermal Engineering', 'Temperature fields, thermal conduction, heating, cooling and thermal stability'],
  ['Electrical, Electronic and Electrochemical Engineering', 'Charge transport, semiconductor behaviour, electrochemistry and dielectric properties'],
  ['Electromagnetic, Photonic and Optical Engineering', 'Electromagnetic response, optical properties, plasmonics and metasurfaces'],
  ['Magnetic and Spintronic Engineering', 'Magnetic properties, spin transport and magnetic-field-dependent behaviour'],
  ['Structural and Nanomechanical Engineering', 'Elasticity, plasticity, stress, fracture, fatigue, buckling and mechanical response'],
  ['Acoustic, Phononic and Elastic-Wave Engineering', 'Acoustic performance, wave propagation, vibration and dynamic properties'],
  ['Microstructure, Phase-Field and Defect Evolution', 'Grain growth, phase transformation, defect evolution and interfaces'],
  ['Quantum, Atomistic and Molecular Simulation', 'DFT, molecular dynamics, Monte Carlo and fundamental material properties'],
  ['Multiscale Modelling and Homogenisation', 'Linking atomic, nano, micro, component and system scales'],
  ['Large-Scale Component and System Analysis', 'Full structures, large arrays, device behaviour and industrial-scale models'],
  ['Active Tuning, Programmability and Control', 'Material switching, field-controlled response, actuators and feedback'],
  ['Fabrication, Synthesis and Process Engineering', 'Chemical reactors, deposition, lithography, printing and self-assembly'],
  ['Material Behaviour, Damage and Lifetime Prediction', 'Functional performance, degradation, failure and durability'],
  ['Experimental Characterisation and Metrology', 'Material property measurement, microscopy, spectroscopy and testing'],
  ['Verification, Calibration, Validation and Uncertainty', 'Numerical credibility, experimental benchmarking and prediction uncertainty'],
  ['Inverse Design, Optimisation and Artificial Intelligence', 'Generative design, topology optimisation, property targeting and ML'],
  ['Industrial Scale-Up, Automation and Quality', 'Production processes, monitoring, manufacturing yield and process control'],
  ['Safety, Standards, Environment and Device Integration', 'Exposure risk, regulatory compliance, product integration and qualification'],
  ['Economics, Sustainability, Data and Lifecycle', 'Financial feasibility, environmental impacts, digital records and end-of-life'],
  ['Nanomaterials and Metamaterials Reverse Engineering and Reconstruction', 'Reconstruct composition, structure, geometry, architecture, properties and manufacturing history from measured data'],
];
const GROUPS = [
  { name: 'Material definition and architecture', suites: [1, 2, 3, 4, 5] },
  { name: 'Chemical, thermal and transport physics', suites: [6, 7, 8, 9, 10] },
  { name: 'Electrical, electromagnetic, mechanical and wave physics', suites: [11, 12, 13, 14, 15, 16] },
  { name: 'Multiscale, large-scale, tunable and lifetime engineering', suites: [17, 18, 19, 20, 21, 22] },
  { name: 'Experiments, validation, optimisation and industrial engineering', suites: [23, 24, 25, 26, 27, 28] },
  { name: 'Reverse engineering', suites: [29] },
];
const suites = SUITES.map(([name, resp], i) => ({ id: i + 1, name, resp, eq: {}, ibc: {}, num: '', cvv: {}, io: { inputs: [], outputs: [] } }));

// split a comma list sentence into items, respecting parentheses
function splitList(s) {
  s = s.replace(/\.$/, '');
  const out = []; let depth = 0, cur = '';
  for (const ch of s) {
    if (ch === '(') depth++; if (ch === ')') depth--;
    if ((ch === ',' || ch === ';') && depth === 0) { out.push(cur); cur = ''; } else cur += ch;
  }
  out.push(cur);
  return out.map(x => x.trim().replace(/^and\s+/, '')).filter(Boolean).map(x => x[0].toUpperCase() + x.slice(1));
}
const personal = /\b(I recommend|I would|I strongly|your |Your |you |attached|previously supplied|previously discussed|supplied case)/;

// ---------- equations & models ----------
{
  const L = read('equations.models');
  let cur = null; const cross = []; let cc = null;
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      const m = b.match(/^Suite (\d+) —/);
      cur = m ? suites[+m[1] - 1] : null;
      const c = b.match(/^([A-E])\. (.+)$/);
      cc = c ? { key: c[1], title: c[2], items: [] } : null;
      if (cc) cross.push(cc);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    const kv = line.match(/^([A-Z][A-Za-z ,–-]+?):\s+(.+)$/);
    if (!kv) continue;
    const [, k, v] = kv;
    if (cur) {
      if (/^Typical equations/.test(k)) cur.eq.formulations = splitList(v);
      else if (/^Classical and contemporary/.test(k)) cur.eq.models = splitList(v);
      else if (/^Hybrid and computational/.test(k)) cur.eq.hybrid = splitList(v);
      else if (/^Typical output data/.test(k)) cur.eq.outputs = splitList(v);
      else if (/^Typical input data/.test(k)) cur.eq.inputs = splitList(v);
    } else if (cc) cc.items.push({ label: k, list: splitList(v) });
  }
  globalThis.crossEq = cross;
}

// ---------- initial & boundary conditions ----------
const bcClasses = []; const bcCapabilities = [];
{
  const L = read('ibcs');
  let cur = null, inTable = false, inCaps = false;
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      const m = b.match(/^(\d+)\. /);
      cur = m && +m[1] <= 29 ? suites[+m[1] - 1] : null;
      inTable = /^A\. General boundary/.test(b);
      inCaps = /^PART VII/.test(b);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    if (inTable && line.startsWith('|') && !/^\|[-|\s]+$/.test(line) && !/Boundary-condition class/.test(line)) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (c.length === 2) bcClasses.push({ cls: c[0], use: c[1] });
    }
    if (inCaps) { const m = line.match(/^\d+\.\s+([^:]+):\s+(.+)$/); if (m) bcCapabilities.push({ name: m[1], text: m[2] }); }
    if (!cur) continue;
    const kv = line.match(/^(Typical initial conditions|Typical boundary conditions[a-z ]*):\s+(.+)$/);
    if (kv) { if (/initial/.test(kv[1])) cur.ibc.ic = kv[2]; else { cur.ibc.bc = kv[2]; cur.ibc.bcLabel = kv[1]; } }
    else if (!personal.test(line)) cur.ibc.note = (cur.ibc.note ? cur.ibc.note + ' ' : '') + line;
  }
}

// ---------- numerical parameters ----------
const numGroups = []; const numModes = [];
{
  const L = read('num.simu.para');
  let g = null, part = '';
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      if (/^PART /.test(b)) { part = b.replace(/^PART [IVX]+ — /, ''); g = null; continue; }
      const m = b.match(/^(\d+)\. (.+)$/);
      g = m ? { id: +m[1], title: m[2], part, text: [] } : null;
      if (g) numGroups.push(g);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    if (line.startsWith('|')) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (c.length === 3 && /^\d+$/.test(c[0])) suites[+c[0] - 1].num = c[2];
      if (c.length === 2 && /^(Automatic|Engineering|Research)/.test(c[0])) numModes.push({ mode: c[0], text: c[1] });
      continue;
    }
    if (g && !personal.test(line)) g.text.push(line);
  }
}

// ---------- CVV ----------
const cvvGroups = []; const gates = []; const metrics = []; const credLevels = [];
{
  const L = read('cvv');
  let cur = null, g = null, part = '', pendingGate = null;
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      if (/^PART /.test(b)) { part = b.replace(/^PART [IVX]+ — /, ''); g = null; cur = null; continue; }
      const s = b.match(/^Suite (\d+) —/);
      cur = s ? suites[+s[1] - 1] : null;
      const m = b.match(/^(\d+)\. (.+)$/);
      g = m ? { id: +m[1], title: m[2], part, text: [] } : null;
      if (g) cvvGroups.push(g);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    const gm = line.match(/^Gate (\d) — (.+)$/);
    if (gm) { pendingGate = { id: +gm[1], name: gm[2], text: '' }; gates.push(pendingGate); continue; }
    if (pendingGate && !pendingGate.text) { pendingGate.text = line; continue; }
    if (line.startsWith('|')) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (c.length === 2 && !/^[-\s]+$/.test(c[0]) && c[0] !== 'Engineering output') metrics.push({ output: c[0], measures: c[1] });
      continue;
    }
    if (/^- /.test(raw) && g && g.id === 55) { credLevels.push(line.replace(/^- /, '').replace(/\.$/, '')); continue; }
    if (cur) {
      const kv = line.match(/^(Calibration|Verification|Validation|Expected results):\s+(.+)$/);
      if (kv) cur.cvv[{ Calibration: 'cal', Verification: 'ver', Validation: 'val', 'Expected results': 'exp' }[kv[1]]] = kv[2];
      else if (!personal.test(line)) cur.cvv.note = (cur.cvv.note ? cur.cvv.note + ' ' : '') + line;
    } else if (g && !personal.test(line)) g.text.push(line);
  }
}

// ---------- input / output ----------
const IOMAP = { 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 8, 7: 6, 8: 7, 9: 7, 10: 9, 11: 10, 12: 11, 13: 14, 14: 12, 15: 13, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 22, 23: 21, 24: 21, 25: 23, 26: 24, 27: 25, 28: 26, 29: 27, 30: 28 };
const units = []; const exportFormats = []; const inputMethods = []; const reportSections = [];
{
  const L = read('input.output');
  let cur = null, mode = null, sect = '', title = '';
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      const m = b.match(/^(\d+)\. (.+)$/);
      cur = null; mode = null; sect = b;
      if (m && IOMAP[+m[1]] && !/Representative input dataset|Relevant outputs/.test(b)) { cur = suites[IOMAP[+m[1]] - 1]; title = m[2]; }
      continue;
    }
    const line = clean(raw); if (!line) continue;
    if (line === 'Typical input data') { mode = 'inputs'; if (cur) cur.io.inputs.push({ group: title, items: [] }); continue; }
    if (line === 'Typical output data') { mode = 'outputs'; if (cur) cur.io.outputs.push({ group: title, items: [] }); continue; }
    if (cur && mode && /^- /.test(line)) { cur.io[mode].at(-1).items.push(line.slice(2).replace(/\.$/, '')); continue; }
    if (line.startsWith('|')) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (c.length !== 2 || /^[-\s]+$/.test(c[0])) continue;
      if (/^C\. Physical units/.test(sect) && c[0] !== 'Engineering variable') units.push({ q: c[0], u: c[1] });
      if (/^C\. CAD, manufacturing/.test(sect) && c[0] !== 'Output type') exportFormats.push({ type: c[0], formats: c[1] });
      if (/^PART X/.test(sect) && c[0] !== 'Input method') inputMethods.push({ method: c[0], use: c[1] });
      continue;
    }
    if (/^D\. Automatic engineering reports/.test(sect)) { const m = line.match(/^\d+\.\s+([^:]+):\s+(.+)$/); if (m) reportSections.push({ name: m[1], text: m[2] }); }
  }
}
// Suite 29 inputs/outputs
suites[28].io.inputs.push({ group: 'Reverse engineering and reconstruction', items: [
  'SEM, TEM, STEM and AFM images', 'XRD and electron diffraction measurements', 'XPS, EDS, Raman, FTIR and other chemical characterisation data',
  'X-ray and electron tomography', 'Mechanical load–displacement and stress–strain measurements', 'Thermal measurement data',
  'Electrical measurements and impedance spectra', 'Optical spectra and electromagnetic S-parameters', 'Acoustic transmission and resonance measurements',
  'Available manufacturing records and process histories', 'Known material compositions and property databases', 'Partial CAD, mesh or unit-cell geometry',
  'Manufacturing constraints and tolerances', ...(suites[28].eq.inputs || [])] });
suites[28].io.outputs.push({ group: 'Reverse engineering and reconstruction', items: [
  'Reconstructed nanostructure geometry', 'Crystal structures, lattice parameters and phase estimates', 'Composition, bonding and surface chemistry estimates',
  'Reconstructed 3D microstructures', 'Identified mechanical properties and constitutive parameters', 'Estimated thermal properties',
  'Identified electrical and dielectric parameters', 'Reconstructed optical and electromagnetic response models', 'Identified acoustic properties and resonator parameters',
  'Estimated or recovered fabrication conditions', 'Ranked candidate constituent materials', 'Reconstructed and completed geometry with uncertainty',
  'Reproduction feasibility and fabrication specifications'] });
suites[28].capabilities = [
  ['Material identification and composition reconstruction', 'Determine constituent elements, compounds, dopants, phases, impurities and chemical composition from experimental data.'],
  ['Atomic and crystal structure reconstruction', 'Recover crystal lattices, atomic arrangements, crystallographic orientations, grain structures and defects from diffraction, imaging and spectroscopy measurements.'],
  ['Nanogeometry and morphology reconstruction', 'Reconstruct nanoparticles, nanorods, nanotubes, nanosheets, nanopores, coatings and nanocomposite geometries from SEM, TEM, AFM and tomography data.'],
  ['Metamaterial architecture reconstruction', 'Identify and reconstruct unit-cell geometries, periodicity, lattice topology, resonator arrangements, chirality, graded architectures and spatial material distributions.'],
  ['Three-dimensional reconstruction', 'Generate physically meaningful 3D geometry, voxel structures, meshes and CAD representations from available 2D and 3D measurements, with uncertainty where the data are incomplete.'],
  ['Material property identification', 'Estimate unknown mechanical, thermal, electrical, electromagnetic, optical, acoustic, magnetic and chemical properties using experimental measurements and identification models.'],
  ['Constitutive model identification', 'Identify and calibrate mathematical material models that reproduce observed physical responses.'],
  ['Manufacturing process reconstruction', 'Estimate plausible synthesis conditions, growth mechanisms, processing histories and fabrication parameters from measured structures and known process constraints.'],
  ['Defect and imperfection reconstruction', 'Identify cracks, vacancies, grain boundaries, porosity, geometric defects, contamination, incomplete deposition and manufacturing irregularities.'],
  ['Multiscale reverse engineering', 'Reconstruct relationships between atomic structure, microstructure, metamaterial architecture and measured macroscopic properties.'],
  ['Performance reconstruction', 'Simulate the reconstructed material under measured test conditions to reproduce mechanical, thermal, electromagnetic, optical, acoustic or chemical responses.'],
  ['Inverse parameter estimation', 'Estimate unknown material parameters, geometric variables and processing parameters by fitting numerical model predictions to laboratory observations.'],
  ['AI-assisted reconstruction', 'Apply image processing, machine learning, statistical inference and optimisation to support geometry reconstruction and identification of material properties.'],
  ['Uncertainty and identifiability assessment', 'Quantify which material features or parameters can be determined reliably from the available observations and which remain ambiguous.'],
  ['Reproduction and redesign', 'Produce reconstructed geometry files, candidate material specifications, proposed fabrication recipes and recommendations for reproducing or modifying the original material.'],
].map(([name, text]) => ({ name, text }));

// ---------- environments ----------
const environments = []; const metaEnv = [];
{
  const L = read('environments'); let cat = '';
  for (const raw of L) {
    const b = unbold(raw);
    if (b) { const m = b.match(/^([A-J])\. (.+)$/); if (m) cat = m[2]; continue; }
    const line = clean(raw);
    if (!line.startsWith('|')) continue;
    const c = line.split('|').map(x => x.trim()).filter(Boolean);
    if (c.length === 3 && /^\d+$/.test(c[0])) environments.push({ id: +c[0], cat, name: c[1], makes: c[2] });
    if (c.length === 2 && /metamaterials|Metasurfaces/.test(c[0]) && c[0] !== 'Metamaterial class') metaEnv.push({ cls: c[0], env: c[1] });
  }
}

// ---------- engineering fields ----------
const fields = [];
{
  for (const raw of read('fields')) {
    const m = clean(raw).match(/^(\d+)\.\s+(.+?) — (.+)$/);
    if (m) fields.push({ id: +m[1], name: m[2], scope: m[3].replace(/\.$/, '') });
  }
}

// ---------- geometry formats ----------
const formatFamilies = []; const importerCaps = []; const priority = [];
{
  const L = read('geometry.formats'); let fam = null, sec = 0;
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      const m = b.match(/^(\d+)\. (.+)$/);
      sec = m ? +m[1] : 0;
      fam = m && sec <= 20 ? { id: sec, title: m[2], header: null, rows: [], notes: [] } : null;
      if (fam) formatFamilies.push(fam);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    if (line.startsWith('|')) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (!c.length || /^[-\s]+$/.test(c[0])) continue;
      if (fam) { if (!fam.header) fam.header = c; else fam.rows.push(c); }
      else if (sec === 21 && c.length === 3 && /^P[123]$/.test(c[2])) priority.push({ family: c[0], formats: c[1], p: c[2] });
      continue;
    }
    if (sec === 22) { const m = line.match(/^\d+\.\s+([^:]+):\s+(.+)$/); if (m) importerCaps.push({ name: m[1], text: m[2] }); continue; }
    if (fam && !personal.test(line) && line.length > 40) fam.notes.push(line);
  }
}

// ---------- framework (lifecycle areas) ----------
const framework = []; const modelFamilies = []; const mfgClasses = []; const principles = [];
{
  const L = read('info'); let part = '', area = null;
  for (const raw of L) {
    const b = unbold(raw);
    if (b) {
      if (/^PART /.test(b)) { part = b.replace(/^PART [IVX]+ — /, ''); area = null; continue; }
      const m = b.match(/^(\d+)\. (.+)$/);
      area = m ? { id: +m[1], title: m[2], part, items: [], notes: [] } : null;
      if (area) framework.push(area);
      continue;
    }
    const line = clean(raw); if (!line) continue;
    if (line.startsWith('|')) {
      const c = line.split('|').map(x => x.trim()).filter(Boolean);
      if (c.length === 2 && !/^[-\s]+$/.test(c[0])) {
        if (area && area.id === 35 && c[0] !== 'Modelling family') modelFamilies.push({ family: c[0], use: c[1] });
        else if (area && area.id === 36 && c[0] !== 'Manufacturing class') mfgClasses.push({ cls: c[0], methods: c[1] });
        else if (area && area.id === 5 && c[0] !== 'Property family') area.items.push(c[0] + ': ' + c[1]);
      }
      continue;
    }
    if (!area) { const p = line.match(/^(First|Second|Third|Fourth|Fifth), (.+)$/); if (p) principles.push(p[2][0].toUpperCase() + p[2].slice(1)); continue; }
    if (/^(-|\d+\.)\s/.test(line)) area.items.push(line.replace(/^(-|\d+\.)\s+/, '').replace(/\.$/, ''));
    else if (!personal.test(line) && line.length > 60) area.notes.push(line);
  }
}

const spec = {
  generated: new Date().toISOString().slice(0, 10),
  suites, groups: GROUPS,
  cross: { equations: globalThis.crossEq, bcClasses, bcCapabilities, numGroups, numModes, cvvGroups, gates, metrics, credLevels, units, exportFormats, inputMethods, reportSections },
  environments, metaEnv, fields, formats: { families: formatFamilies, priority, importerCaps },
  framework, modelFamilies, mfgClasses, principles,
};
fs.mkdirSync('data', { recursive: true });
fs.writeFileSync('data/spec.json', JSON.stringify(spec));
const rep = suites.map(s => `${s.id}: eq ${[s.eq.formulations, s.eq.models, s.eq.hybrid, s.eq.outputs].map(a => a ? a.length : 0).join('/')} ibc ${s.ibc.ic ? 1 : 0}${s.ibc.bc ? 1 : 0} num ${s.num ? 1 : 0} cvv ${Object.keys(s.cvv).length} io ${s.io.inputs.reduce((a, g) => a + g.items.length, 0)}/${s.io.outputs.reduce((a, g) => a + g.items.length, 0)}`);
console.log(rep.join('\n'));
console.log({ env: environments.length, metaEnv: metaEnv.length, fields: fields.length, fam: formatFamilies.length, prio: priority.length, caps: importerCaps.length, framework: framework.length, fwItems: framework.reduce((a, f) => a + f.items.length, 0), numGroups: numGroups.length, cvvGroups: cvvGroups.length, gates: gates.length, metrics: metrics.length, cred: credLevels.length, units: units.length, exp: exportFormats.length, bc: bcClasses.length, bcCaps: bcCapabilities.length, cross: globalThis.crossEq.length, mf: modelFamilies.length, mfg: mfgClasses.length, principles: principles.length, bytes: fs.statSync('data/spec.json').size });
