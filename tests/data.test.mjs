// Consistency tests for the materials library (data/materials.js).
// Usage: node -e "import('./tests/data.test.mjs').then(m=>m.default()).then(r=>{console.log(r.log.join('\n'));process.exit(r.fail?1:0)})"
import { CLASSES, MATERIALS, NEGATIVE_CLASS_IDS, findMaterial, findClass } from '../data/materials.js';
import * as chem from '../js/core/chem.js';

const NEGATIVE = ['negative-index', 'negative-permittivity', 'negative-permeability', 'auxetic', 'negative-thermal-expansion', 'negative-compressibility', 'negative-stiffness', 'negative-differential-resistance', 'negative-effective-mass', 'negative-differential-thermal-resistance'];
const REQUIRED_CLASSES = [...NEGATIVE, 'epsilon-near-zero', 'hyperbolic', 'chiral', 'photonic-crystal', 'phononic-crystal', 'acoustic-metasurface', 'optical-metasurface', 'thermal-metamaterial', 'seismic', 'mechanical-lattice', 'pentamode', 'origami-kirigami', 'tunable', 'nonlinear', 'topological', 'bianisotropic-willis',
  'metal-nanoparticle', 'oxide-nanoparticle', 'quantum-dot', 'magnetic-nanoparticle', 'fullerene', 'nanorod', 'nanowire', 'nanotube', 'nanofibre', 'graphene', 'tmd', 'mxene', 'hbn', 'nanosheet', 'nanofilm', 'nanoporous', 'nanocomposite', 'nanocrystalline', 'core-shell', 'nanolattice'];
const REQUIRED_MATERIALS = ['zno-nanorods', 'zno-nanowires', 'zno-nanoparticles', 'tio2-anatase', 'tio2-rutile', 'tio2-nanotubes', 'au-nanoparticles', 'ag-nanoparticles', 'cu-nanoparticles', 'pt-nanoparticles', 'pd-nanoparticles', 'fe3o4-nanoparticles', 'ceo2-nanoparticles', 'sio2-nanoparticles', 'al2o3-nanoparticles', 'cuo-nanoparticles', 'nio-nanoparticles', 'co3o4-nanoparticles', 'sno2-nanoparticles', 'wo3-nanoparticles', 'cdse-qd', 'cds-qd', 'zns-qd', 'pbs-qd', 'inp-qd', 'si-qd', 'cspbbr3-nc', 'carbon-dots', 'c60', 'swcnt', 'mwcnt', 'si-nanowires', 'gan-nanowires', 'ag-nanowires', 'graphene', 'graphene-oxide', 'mos2', 'ws2', 'hbn', 'ti3c2-mxene', 'silica-aerogel', 'zif-8', 'hydroxyapatite-nano', 'nanocellulose', 'plga-nanoparticles', 'liposomes', 'nanodiamond', 'batio3-nano', 'lifepo4-nano',
  'silicon', 'fused-silica', 'alumina', 'copper', 'aluminium-6061', 'ti6al4v', 'ss316l', 'pla', 'tpu', 'pdms', 'ip-dip', 'pmma', 'fr4', 'ptfe-laminate', 'gaas', 'vo2', 'gst', 'pzt', 'niti', 'epoxy'];
const KINDS = ['nanomaterial', 'metamaterial', 'nanostructured-metamaterial', 'constituent'];
const DIMS = ['0D', '1D', '2D', '3D', null];
const CRYSTALS = ['ZnO_wurtzite', 'ZnO_zincblende', 'TiO2_rutile', 'TiO2_anatase', 'Si', 'Ge', 'C_diamond', 'graphite', 'graphene', 'hBN', 'MoS2_2H', 'GaAs', 'GaN_wurtzite', 'CdSe_wurtzite', 'CdS_zincblende', 'ZnS_zincblende', 'Au', 'Ag', 'Cu', 'Pt', 'Pd', 'Ni', 'Al', 'Fe_bcc', 'W', 'NaCl', 'MgO', 'CeO2', 'SrTiO3', 'BaTiO3_cubic', 'CsCl', 'Fe3O4', 'SiO2_cristobalite', 'Al2O3_corundum'];
const ARCH = ['reentrant2D', 'chiral2D', 'rotatingSquares2D', 'honeycomb2D', 'octet', 'bcc', 'kelvin', 'gyroid', 'schwarzP', 'srr', 'doubleSrr', 'fishnet', 'wire', 'cross', 'nanopillar', 'nanobar', 'multilayer', 'massInMass', 'helmholtz', 'membrane', 'bistableBeam', 'wineRack', 'bimaterialLattice', 'kirigami', 'other'];
const PROP_KEYS = ['density', 'E', 'nu', 'yield', 'KIc', 'k_thermal', 'cp', 'cte', 'sigma_elec', 'eps_r', 'n_refr', 'bandgap_eV', 'Ms', 'surface_energy', 'c_sound', 'melting_K', 'cost_usd_kg', 'co2_kg_kg'];
// [min, max] inclusive sanity ranges (SI) for solids, aerogels and polymers
const RANGE = {
  density: [1, 23000], E: [1e3, 1.3e12], nu: [-1, 0.5], yield: [1e5, 2e11], KIc: [1e4, 3e8], k_thermal: [0.005, 6000], cp: [50, 5000], cte: [-1e-3, 1e-3],
  sigma_elec: [1e-20, 1e8], eps_r: [1, 1e5], n_refr: [1, 8], bandgap_eV: [0, 12], Ms: [1, 2e6], surface_energy: [0.005, 5], c_sound: [10, 2e4], melting_K: [100, 4500], cost_usd_kg: [0.01, 1e8], co2_kg_kg: [0.01, 1e5],
};
const FORBIDDEN = /\b(doi|et al\.?|arxiv|https?:)\b/i;

export default async function () {
  let pass = 0, fail = 0; const log = [];
  const ok = (cond, msg) => { if (cond) pass++; else { fail++; log.push('  ✗ ' + msg); } return cond; };
  const isNum = v => typeof v === 'number' && Number.isFinite(v);
  const str = (v, n = 1) => typeof v === 'string' && v.trim().length >= n;
  const range2 = (a, lo) => Array.isArray(a) && a.length === 2 && isNum(a[0]) && isNum(a[1]) && a[0] <= a[1] && a[0] >= lo;

  // ---- classes
  const classIds = new Set();
  for (const c of CLASSES) {
    ok(str(c.id) && !classIds.has(c.id), `class id missing or duplicated: ${c.id}`); classIds.add(c.id);
    ok(str(c.name, 4) && str(c.signature, 8) && str(c.mechanism, 80) && str(c.caveats, 60), `class ${c.id}: name/signature/mechanism/caveats incomplete`);
    ok(['nanomaterial', 'metamaterial', 'constituent'].includes(c.kind), `class ${c.id}: bad kind ${c.kind}`);
    ok(Array.isArray(c.suites) && c.suites.length > 0 && c.suites.every(s => Number.isInteger(s) && s >= 1 && s <= 29), `class ${c.id}: suites must be ids 1–29`);
    ok(findClass(c.id) === c, `findClass(${c.id}) failed`);
  }
  for (const id of REQUIRED_CLASSES) ok(classIds.has(id), `required class missing: ${id}`);
  ok(NEGATIVE.every(id => NEGATIVE_CLASS_IDS.includes(id)) && NEGATIVE_CLASS_IDS.length === NEGATIVE.length, 'NEGATIVE_CLASS_IDS does not match the ten negative classes');

  // ---- materials
  const ids = new Set(), perClass = {}, perKind = {}; let nRx = 0, nIonic = 0;
  for (const m of MATERIALS) {
    const t = `[${m.id}]`;
    ok(str(m.id) && /^[a-z0-9-]+$/.test(m.id) && !ids.has(m.id), `${t} id missing, malformed or duplicated`); ids.add(m.id);
    ok(findMaterial(m.id) === m, `${t} findMaterial failed`);
    ok(str(m.name, 3) && str(m.description, 60), `${t} name/description incomplete`);
    ok(KINDS.includes(m.kind), `${t} bad kind ${m.kind}`);
    ok(classIds.has(m.classId), `${t} unknown classId ${m.classId}`);
    ok(DIMS.includes(m.dimensionality), `${t} bad dimensionality ${m.dimensionality}`);
    perClass[m.classId] = (perClass[m.classId] || 0) + 1; perKind[m.kind] = (perKind[m.kind] || 0) + 1;
    const cls = findClass(m.classId), meta = m.kind === 'metamaterial' || m.kind === 'nanostructured-metamaterial';
    if (cls) {
      if (m.kind === 'constituent') ok(cls.kind === 'constituent', `${t} constituent must sit in the constituent class`);
      if (m.kind === 'nanomaterial') ok(cls.kind === 'nanomaterial', `${t} nanomaterial in non-nanomaterial class ${cls.id}`);
      if (m.kind === 'metamaterial') ok(cls.kind === 'metamaterial', `${t} metamaterial in non-metamaterial class ${cls.id}`);
    }
    if (m.kind === 'nanomaterial') ok(m.dimensionality !== null, `${t} nanomaterial needs a dimensionality`);
    // formula
    if (m.formula !== null) {
      try { const f = chem.parseFormula(m.formula); ok(Object.keys(f.atoms).length > 0 && chem.molarMass(m.formula) > 0, `${t} formula ${m.formula} parsed empty`); } catch (e) { ok(false, `${t} formula ${m.formula} does not parse: ${e.message}`); }
    } else pass++;
    // props
    ok(m.props && PROP_KEYS.every(k => k in m.props) && Object.keys(m.props).length === PROP_KEYS.length, `${t} props must hold exactly the 18 standard keys`);
    for (const k of PROP_KEYS) {
      const v = m.props?.[k]; if (v === null || v === undefined) continue;
      ok(isNum(v) && v >= RANGE[k][0] && v <= RANGE[k][1], `${t} props.${k} = ${v} outside sane range [${RANGE[k]}]`);
    }
    const p = m.props || {};
    if (p.E !== null && p.yield !== null) ok(p.yield < p.E, `${t} yield must be below E`);
    if (p.E !== null && p.density !== null) { const c = Math.sqrt(p.E / p.density); ok(c > 20 && c < 3e4, `${t} bar wave speed sqrt(E/ρ) = ${c.toFixed(0)} m/s implausible`); }
    if (p.c_sound !== null && p.E !== null && p.density !== null) { const r = p.c_sound / Math.sqrt(p.E / p.density); ok(r > 0.6 && r < 2.2, `${t} c_sound inconsistent with sqrt(E/ρ) (ratio ${r.toFixed(2)})`); }
    if (p.n_refr !== null && p.eps_r !== null) ok(p.eps_r >= 0.8 * p.n_refr ** 2 || p.sigma_elec > 1, `${t} static eps_r below optical n² `);
    // effective / architecture / crystal
    ok(m.crystal === null || CRYSTALS.includes(m.crystal), `${t} unknown crystal key ${m.crystal}`);
    if (meta) {
      ok(m.effective && typeof m.effective === 'object' && str(m.effective.note, 20), `${t} metamaterial needs effective values with a note`);
      ok(m.architecture && ARCH.includes(m.architecture.type) && m.architecture.dims && Object.keys(m.architecture.dims).length > 0, `${t} metamaterial needs a valid architecture {type, dims}`);
      ok(['architected', 'intrinsic', 'device', 'model'].includes(m.origin), `${t} bad origin ${m.origin}`);
    }
    if (m.architecture) for (const [k, v] of Object.entries(m.architecture.dims)) ok(isNum(v), `${t} architecture.dims.${k} must be a finite number`);
    if (m.effective) for (const [k, v] of Object.entries(m.effective)) ok(v === null || isNum(v) || typeof v === 'string' || (Array.isArray(v) && v.every(isNum)), `${t} effective.${k} must be number, string, numeric array or null`);
    if (m.effective?.nu != null) ok(m.effective.nu >= -1 && m.effective.nu <= 0.5, `${t} effective ν outside the stable isotropic range`);
    if (m.classId === 'auxetic' && m.effective?.nu != null) ok(m.effective.nu < 0, `${t} auxetic entry needs ν < 0`);
    for (const k of ['band_Hz', 'band_m', 'range_K', 'range_Pa']) if (m.effective?.[k]) ok(range2(m.effective[k], 0), `${t} effective.${k} must be [lo, hi] ≥ 0`);
    // synthesis
    ok(Array.isArray(m.synthesis) && m.synthesis.length >= 1 && m.synthesis.length <= 3, `${t} needs 1–3 synthesis routes`);
    for (const s of m.synthesis || []) {
      ok(str(s.route, 10) && typeof s.notes === 'string', `${t} synthesis route text missing`);
      ok(Array.isArray(s.envIds) && s.envIds.length > 0 && s.envIds.every(e => Number.isInteger(e) && e >= 1 && e <= 146) && new Set(s.envIds).size === s.envIds.length, `${t} envIds must be unique integers within 1–146 (${s.envIds})`);
      ok(range2(s.T_C, -273.15) && s.T_C[1] <= 4500, `${t} T_C window invalid (${s.T_C})`);
      ok(range2(s.P_bar, 0) && s.P_bar[0] > 0 && s.P_bar[1] <= 1e6, `${t} P_bar window invalid (${s.P_bar})`);
      ok(range2(s.time_h, 0) && s.time_h[0] > 0 && s.time_h[1] <= 1e4, `${t} time_h window invalid (${s.time_h})`);
      ok(Array.isArray(s.precursors) && s.precursors.length > 0 && s.precursors.every(x => str(x)), `${t} precursors missing`);
      if (s.reaction) {
        nRx++;
        const R = s.reaction.reactants, Pd = s.reaction.products;
        ok(Array.isArray(R) && Array.isArray(Pd) && R.length > 0 && Pd.length > 0 && [...R, ...Pd].every(x => str(x.f) && [null, 's', 'l', 'g', 'aq'].includes(x.phase)), `${t} reaction species malformed`);
        try {
          const rf = R.map(x => x.f), pf = Pd.map(x => x.f), b = chem.balance(rf, pf);
          ok(b.coefficients.every(c => Number.isInteger(c) && c > 0), `${t} non-integer coefficients`);
          const a = chem.auditEquation(rf, pf, b.coefficients);
          ok(a.massBalanced && a.chargeBalanced, `${t} reaction fails mass/charge audit: ${chem.equationString(rf, pf)}`);
          ok(Math.max(...b.coefficients) <= 100, `${t} suspiciously large coefficients: ${chem.equationString(rf, pf)}`);
        } catch (e) { ok(false, `${t} reaction does not balance (${R.map(x => x.f).join(' + ')} -> ${Pd.map(x => x.f).join(' + ')}): ${e.message}`); }
        if ([...R, ...Pd].some(x => x.phase === 'aq')) {
          try { const io = chem.ionicEquations(R, Pd); nIonic++; ok(str(io.molecular, 3) && str(io.completeIonic, 3) && str(io.netIonic, 3) && Array.isArray(io.spectators), `${t} ionic equations incomplete`); } catch (e) { ok(false, `${t} ionicEquations threw: ${e.message}`); }
        }
      }
    }
    // text fields
    ok(Array.isArray(m.applications) && m.applications.length > 0 && m.applications.every(x => str(x, 3)), `${t} applications missing`);
    ok(Array.isArray(m.hazards) && m.hazards.length > 0 && m.hazards.every(x => str(x, 8)), `${t} hazards missing`);
    ok(Array.isArray(m.suites) && m.suites.length > 0 && m.suites.every(s => Number.isInteger(s) && s >= 1 && s <= 29), `${t} suites must be ids 1–29`);
    ok(m.search && str(m.search.terms, 6) && (m.search.pubchem === null || str(m.search.pubchem, 2)), `${t} search {pubchem, terms} malformed`);
    ok(!FORBIDDEN.test(JSON.stringify(m)), `${t} contains a citation-like string`);
    if (m.constituent != null) ok(ids.has(m.constituent), `${t} constituent ${m.constituent} must be defined earlier`);
  }

  // ---- coverage and counts
  for (const id of REQUIRED_MATERIALS) ok(ids.has(id), `required material missing: ${id}`);
  for (const id of NEGATIVE) ok((perClass[id] || 0) >= 3, `negative class ${id} has only ${perClass[id] || 0} materials (need ≥ 3)`);
  for (const c of CLASSES) ok((perClass[c.id] || 0) >= 1, `class ${c.id} has no material`);
  const nNano = perKind.nanomaterial || 0, nMeta = (perKind.metamaterial || 0) + (perKind['nanostructured-metamaterial'] || 0), nCons = perKind.constituent || 0;
  ok(MATERIALS.length >= 110, `only ${MATERIALS.length} materials (need ≥ 110)`);
  ok(nNano >= 50, `only ${nNano} nanomaterials (need ≥ 50)`);
  ok(nMeta >= 45, `only ${nMeta} metamaterials (need ≥ 45)`);
  ok(nCons >= 15, `only ${nCons} constituents (need ≥ 15)`);
  ok(findMaterial('no-such-material') === null, 'findMaterial must return null for unknown ids');

  // spot checks against well-known reference values
  const zno = findMaterial('zno-nanorods');
  ok(zno && zno.props.density === 5606 && Math.abs(zno.props.bandgap_eV - 3.37) < 1e-9 && zno.crystal === 'ZnO_wurtzite', 'ZnO reference values changed');
  ok(findMaterial('zrw2o8')?.effective.cte < 0 && findMaterial('ag3cocn6')?.effective.beta_L < 0 && findMaterial('srr-wire-nim')?.effective.n < 0, 'negative reference values must be negative');
  const mim = findMaterial('mass-in-mass-lattice');
  if (mim) { const d = mim.architecture.dims, f0 = Math.sqrt(d.inner_spring_N_per_m / (d.inner_mass_g / 1000)) / (2 * Math.PI), f1 = f0 * Math.sqrt(1 + d.inner_mass_g / d.outer_mass_g); ok(Math.abs(f0 - mim.effective.band_Hz[0]) < 1 && Math.abs(f1 - mim.effective.band_Hz[1]) < 1, `mass-in-mass band inconsistent with dims (${f0.toFixed(1)}–${f1.toFixed(1)} Hz)`); }

  log.unshift(`Materials library: ${MATERIALS.length} materials (${nNano} nanomaterials, ${nMeta} metamaterials, ${nCons} constituents), ${CLASSES.length} classes, ${nRx} balanced reactions, ${nIonic} ionic equation sets`,
    'Negative classes: ' + NEGATIVE.map(id => `${id} ${perClass[id] || 0}`).join(', '));
  log.push(`${pass} checks passed, ${fail} failed`);
  return { pass, fail, log };
}
