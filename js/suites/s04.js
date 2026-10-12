// Suite 4 — Nanomaterial geometry and morphology.
// Shape descriptors of single nanostructures, particle-size populations, nanorod arrays on substrates, stochastic
// porous microstructures, aggregates and packing, and thin-film stacks with self-affine roughness.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import * as atomic from '../geometry/atomic.js';
import * as implicit from '../geometry/implicit.js';
import * as msh from '../geometry/mesh.js';

const F = N.fmt, PI = Math.PI, { sqrt, abs, min, max, log, exp, sin, cos, round, floor, cbrt } = Math;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const nm = v => v * 1e9;
const fin = (v, d = 0) => (Number.isFinite(v) ? v : d);
const toOBJ = m => `# units: ${m.units || 'nm'}\n` + m.vertices.map(p => `v ${+p[0].toPrecision(7)} ${+p[1].toPrecision(7)} ${+(p[2] || 0).toPrecision(7)}`).join('\n') + '\n' + m.faces.map(f => `f ${f[0] + 1} ${f[1] + 1} ${f[2] + 1}`).join('\n') + '\n';

// Handbook-typical bulk densities (kg/m³) and the crystal used for atomistic models.
const MAT = {
  ZnO: { f: 'ZnO', rho: 5606, x: 'ZnO_wurtzite', name: 'zinc oxide (wurtzite)' },
  TiO2: { f: 'TiO2', rho: 3890, x: 'TiO2_anatase', name: 'titanium dioxide (anatase)' },
  Au: { f: 'Au', rho: 19300, x: 'Au', name: 'gold' },
  Ag: { f: 'Ag', rho: 10490, x: 'Ag', name: 'silver' },
  Pt: { f: 'Pt', rho: 21450, x: 'Pt', name: 'platinum' },
  Cu: { f: 'Cu', rho: 8960, x: 'Cu', name: 'copper' },
  Si: { f: 'Si', rho: 2329, x: 'Si', name: 'silicon' },
  SiO2: { f: 'SiO2', rho: 2200, x: null, name: 'amorphous silica' },
  Fe3O4: { f: 'Fe3O4', rho: 5170, x: 'Fe3O4', name: 'magnetite' },
  CeO2: { f: 'CeO2', rho: 7215, x: 'CeO2', name: 'ceria' },
  Al2O3: { f: 'Al2O3', rho: 3987, x: 'Al2O3_corundum', name: 'α-alumina' },
  C: { f: 'C', rho: 2260, x: 'graphite', name: 'graphitic carbon' },
  CdSe: { f: 'CdSe', rho: 5816, x: 'CdSe_wurtzite', name: 'cadmium selenide (wurtzite)' },
};
const MAT_OPTS = Object.keys(MAT).map(k => ({ value: k, label: `${chem.prettyFormula(MAT[k].f)} — ${MAT[k].name}` }));
function matInfo(f, rho, name, x) {
  const p = chem.parseFormula(f), M = chem.molarMass(p) * 1e-3; let nAt = 0; for (const e in p.atoms) nAt += p.atoms[e];
  if (!(M > 0) || !(nAt > 0)) throw new Error('empty formula');
  const vAt = M / (rho * N.N_A * nAt);
  return { f, rho, name, x, M, nAt, vAt, delta: cbrt(vAt), pretty: chem.prettyFormula(f) };
}
function material(v, warnings) {
  if (v.material !== 'custom') { const m = MAT[v.material]; return matInfo(m.f, m.rho, m.name, m.x); }
  try { return matInfo(v.formula.trim(), v.rho, 'user-defined material', null); } catch (e) {
    warnings.push(`The custom formula "${v.formula}" could not be parsed (${e.message}); ZnO was used instead.`);
    return matInfo('ZnO', MAT.ZnO.rho, MAT.ZnO.name, MAT.ZnO.x);
  }
}
const matEq = m => eqn('formula', 'Material', `${m.pretty}: M = ${F(m.M * 1e3, 5)} g/mol, ρ = ${F(m.rho, 5)} kg/m³, ${m.nAt} atom${m.nAt > 1 ? 's' : ''} per formula unit`, `Mean atomic volume v_a = M/(ρ N_A n) = ${F(m.vAt * 1e30, 4)} Å³, atomic spacing δ = v_a^(1/3) = ${F(m.delta * 1e10, 4)} Å.`);

// ---------------------------------------------------------------- demonstration tables
const DEMO_SIZES = (() => { const r = N.rng(42), a = ['diameter_nm']; for (let i = 0; i < 160; i++) a.push((r.lognormal(log(15), log(1.25)) + r.normal(0, 0.3)).toFixed(2)); return a.join('\n'); })();
const DEMO_STACK = 'layer,formula,thickness_nm,density_kg_m3\nSeed layer,ZnO,30,5606\nBarrier,Al2O3,10,3100\nFunctional film,TiO2,60,3890';
function parseRows(text) {
  return String(text || '').split(/\r?\n/).map(l => l.trim()).filter(l => l && !/^[#%]/.test(l)).map(l => l.split(/[,;\t]/).map(c => c.trim()));
}

// ---------------------------------------------------------------- inputs
const SH = m => ({ mode: Array.isArray(m) ? m : [m] });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'single', group: 'Model', options: [
    { value: 'single', label: 'Single nanostructure descriptors and tolerances' }, { value: 'population', label: 'Particle population / size distribution' },
    { value: 'rods', label: 'Nanorod / nanowire array on a substrate' }, { value: 'porous', label: 'Porous microstructure (RVE) generation and analysis' },
    { value: 'aggregate', label: 'Aggregates and packing (DLA, RSA)' }, { value: 'film', label: 'Thin film / coating / multilayer' }], help: 'Selects the morphology model. Every analysis builds a geometry for the viewer and publishes its descriptors to the shared bus.' },
  { key: 'material', label: 'Material', type: 'select', value: 'ZnO', group: 'Material', options: [...MAT_OPTS, { value: 'custom', label: 'Custom (formula and density below)' }], showIf: SH(['single', 'population', 'rods', 'porous', 'aggregate']), help: 'Sets the bulk density and molar mass used for specific surface area, mass and atom counts, and the crystal used for atomistic models of small particles.' },
  { key: 'formula', label: 'Custom chemical formula', type: 'text', value: 'ZnO', group: 'Material', bus: 'material.formula', showIf: { mode: ['single', 'population', 'rods', 'porous', 'aggregate'], material: ['custom'] }, help: 'Formula of the solid, e.g. ZnO, TiO2, Fe3O4. Used for the molar mass and the number of atoms per formula unit.' },
  { key: 'rho', label: 'Custom bulk density', unit: 'kg/m³', value: 5606, min: 100, max: 23000, group: 'Material', bus: 'prop.density', showIf: { mode: ['single', 'population', 'rods', 'porous', 'aggregate'], material: ['custom'] }, help: 'Bulk density of the custom material. Specific surface area scales as 1/ρ.' },
  { key: 'shape', label: 'Shape', type: 'select', value: 'sphere', group: 'Nanostructure', showIf: SH('single'), options: [
    { value: 'sphere', label: 'Sphere' }, { value: 'spheroid', label: 'Spheroid (equatorial d, polar axis L)' }, { value: 'rod_hex', label: 'Hexagonal rod / wire (across flats d, length L)' },
    { value: 'rod_circ', label: 'Circular rod / wire / fibre' }, { value: 'tube', label: 'Tube (outer d, wall t, length L)' }, { value: 'sheet', label: 'Sheet / platelet (side d, thickness t)' },
    { value: 'coreshell', label: 'Core–shell sphere (core d, shell t)' }, { value: 'cube', label: 'Cube (edge d)' }, { value: 'octahedron', label: 'Octahedron (edge d)' }], help: 'Geometric model of the nanostructure. The meaning of d, L and t for each shape is given in the option label.' },
  { key: 'd', label: 'Diameter / edge d', unit: 'm', value: 1e-8, min: 1e-9, max: 1e-5, log: true, group: 'Nanostructure', bus: 'growth.diameter', showIf: SH(['single', 'population', 'rods', 'aggregate']), help: 'Characteristic transverse size: sphere, rod or tube diameter, cube edge, sheet side, core diameter, median particle diameter of a population, or primary-particle diameter of an aggregate. Surface-to-volume ratio scales as 1/d.' },
  { key: 'L', label: 'Length L', unit: 'm', value: 1e-6, min: 1e-9, max: 1e-3, log: true, group: 'Nanostructure', bus: 'growth.length', showIf: SH(['single', 'rods']), help: 'Rod, wire or tube length, or the polar axis of a spheroid. Sets the aspect ratio L/d.' },
  { key: 'tw', label: 'Wall / shell / sheet thickness t', unit: 'm', value: 2e-9, min: 2e-10, max: 1e-6, log: true, group: 'Nanostructure', showIf: SH('single'), help: 'Tube wall thickness (limited to just under d/2), shell thickness of a core–shell particle, or sheet thickness.' },
  { key: 'shellMat', label: 'Shell material', type: 'select', value: 'SiO2', group: 'Nanostructure', options: MAT_OPTS, showIf: { mode: ['single'], shape: ['coreshell'] }, help: 'Material of the shell of a core–shell particle; it contributes its own density to the particle mass.' },
  { key: 'tolPct', label: 'Dimensional tolerance (±)', unit: '%', value: 5, min: 0, max: 30, group: 'Tolerances', showIf: SH('single'), help: 'Relative tolerance applied to every dimension. First-order propagation gives the resulting spread of volume, surface area, specific surface area and aspect ratio.' },
  { key: 'popSource', label: 'Population source', type: 'select', value: 'generate', group: 'Population', showIf: SH('population'), options: [{ value: 'generate', label: 'Generate from a distribution (seeded sampling)' }, { value: 'fit', label: 'Fit distributions to measured diameters' }], help: 'Either sample a model distribution, or fit log-normal, normal and Weibull distributions to a table of measured diameters.' },
  { key: 'dist', label: 'Size distribution', type: 'select', value: 'lognormal', group: 'Population', showIf: SH('population'), options: [{ value: 'lognormal', label: 'Log-normal (median d, geometric std σ_g)' }, { value: 'normal', label: 'Normal (mean d, coefficient of variation)' }, { value: 'weibull', label: 'Weibull (median d, shape k)' }], help: 'Distribution sampled (generate) or reported as the selected fit (fit). Log-normal is typical of nucleation-and-growth and aerosol products.' },
  { key: 'gsd', label: 'Geometric standard deviation σ_g', unit: '–', value: 1.3, min: 1.01, max: 3, group: 'Population', showIf: SH('population'), help: 'Width of the log-normal distribution; 1.0 is monodisperse, 1.1–1.2 narrow, above 1.5 broad.' },
  { key: 'cvPct', label: 'Coefficient of variation (normal)', unit: '%', value: 15, min: 1, max: 60, group: 'Population', showIf: SH('population'), help: 'Standard deviation divided by the mean for the normal distribution. Negative samples are rejected (truncated normal).' },
  { key: 'wk', label: 'Weibull shape k', unit: '–', value: 4, min: 0.8, max: 20, group: 'Population', showIf: SH('population'), help: 'Shape parameter of the Weibull distribution; larger values give a narrower distribution.' },
  { key: 'nSample', label: 'Number of particles sampled', unit: '–', value: 2000, min: 50, max: 20000, step: 50, group: 'Population', sens: false, showIf: SH('population'), help: 'Sample size for the generated population. Volume-weighted statistics converge more slowly than number-weighted ones.' },
  { key: 'sizes', label: 'Measured diameters', type: 'table', columns: ['diameter_nm'], value: DEMO_SIZES, group: 'Population', showIf: SH('population'), help: 'One measured diameter per row in nm (e.g. from TEM/SEM image analysis). Used when the population source is "fit". The default is synthetic demonstration data.' },
  { key: 'nDens', label: 'Rod number density', unit: '1/m²', value: 5e13, min: 1e10, max: 1e16, log: true, group: 'Array', showIf: SH('rods'), help: 'Rods per unit substrate area (1e12 1/m² = 1 per µm²). Together with the diameter it sets the areal coverage and the inter-rod gap.' },
  { key: 'cvL', label: 'Length spread (CV)', unit: '%', value: 15, min: 0, max: 60, group: 'Array', showIf: SH('rods'), help: 'Coefficient of variation of rod length (log-normal).' },
  { key: 'cvD', label: 'Diameter spread (CV)', unit: '%', value: 12, min: 0, max: 60, group: 'Array', showIf: SH('rods'), help: 'Coefficient of variation of rod diameter (log-normal).' },
  { key: 'tilt', label: 'Tilt spread', unit: '°', value: 8, min: 0, max: 45, group: 'Array', showIf: SH('rods'), help: 'Standard deviation of the tilt from the substrate normal (half-normal distribution of the tilt angle, uniform azimuth). Zero gives perfectly vertical rods.' },
  { key: 'arrangement', label: 'Arrangement', type: 'select', value: 'random', group: 'Array', showIf: SH('rods'), options: [{ value: 'random', label: 'Random (hard-core, seeded)' }, { value: 'hexagonal', label: 'Hexagonal pattern (templated)' }, { value: 'square', label: 'Square pattern (templated)' }], help: 'Random nucleation on a seed layer, or patterned growth through a lithographic template.' },
  { key: 'section', label: 'Rod cross-section', type: 'select', value: 'hex', group: 'Array', showIf: SH('rods'), options: [{ value: 'hex', label: 'Hexagonal (d across flats)' }, { value: 'round', label: 'Circular' }], help: 'Hexagonal prisms are typical of wurtzite nanorods; circular for fibres and VLS wires.' },
  { key: 'pType', label: 'Microstructure model', type: 'select', value: 'randomSpheres', group: 'Porous medium', showIf: SH('porous'), options: [{ value: 'randomSpheres', label: 'Overlapping spheres (sintered powder)' }, { value: 'gaussianField', label: 'Level-cut random field (spinodal-like)' }, { value: 'fibres', label: 'Random fibres (nonwoven / nanofibre mat)' }, { value: 'voronoiFoam', label: 'Voronoi open-cell foam' }], help: 'Stochastic generator of the representative volume element. Solid is phase 1 and pore is phase 0.' },
  { key: 'phi', label: 'Target porosity', unit: '–', value: 0.25, min: 0.05, max: 0.9, group: 'Porous medium', showIf: SH('porous'), help: 'Pore volume fraction requested from the generator. The achieved value is reported and published.' },
  { key: 'feat', label: 'Feature size', unit: 'm', value: 2e-8, min: 2e-9, max: 1e-5, log: true, group: 'Porous medium', showIf: SH('porous'), help: 'Grain or fibre diameter, correlation length or foam cell size. Specific surface scales as 1/feature size and permeability as its square.' },
  { key: 'featVox', label: 'Voxels per feature', unit: '–', value: 7, min: 3, max: 12, step: 1, group: 'Porous medium', sens: false, showIf: SH('porous'), help: 'Spatial resolution of the voxel grid. Fewer than about five voxels per feature degrades the surface-area and tortuosity estimates.' },
  { key: 'pN', label: 'Grid size per edge', unit: 'voxels', value: 28, min: 16, max: 36, step: 2, group: 'Porous medium', sens: false, showIf: SH('porous'), help: 'Edge length of the cubic representative volume element in voxels. The RVE should span at least four to five features.' },
  { key: 'aggType', label: 'Aggregation / packing model', type: 'select', value: 'dla', group: 'Aggregate', showIf: SH('aggregate'), options: [{ value: 'dla', label: 'Diffusion-limited aggregation (3-D lattice, particle–cluster)' }, { value: 'rsa', label: 'Random sequential adsorption of discs on a surface' }], help: 'DLA grows a fractal cluster from random walkers; RSA deposits non-overlapping particles irreversibly on a surface up to the jamming limit.' },
  { key: 'nPart', label: 'Primary particles in the aggregate', unit: '–', value: 600, min: 20, max: 1500, step: 10, group: 'Aggregate', showIf: SH('aggregate'), help: 'Number of primary particles grown into the DLA cluster. More particles give a better estimate of the fractal dimension.' },
  { key: 'stick', label: 'Sticking probability', unit: '–', value: 1, min: 0.05, max: 1, group: 'Aggregate', showIf: SH('aggregate'), help: 'Probability that a walker touching the cluster sticks. Values below one let particles penetrate and give denser clusters (towards reaction-limited behaviour).' },
  { key: 'rsaL', label: 'RSA domain side', unit: 'diameters', value: 30, min: 10, max: 60, step: 1, group: 'Aggregate', sens: false, showIf: SH('aggregate'), help: 'Side of the periodic square substrate in particle diameters. Larger domains reduce the statistical scatter of the coverage.' },
  { key: 'rsaTau', label: 'RSA deposition attempts', unit: 'per disc area', value: 150, min: 10, max: 400, group: 'Aggregate', sens: false, showIf: SH('aggregate'), help: 'Dimensionless deposition time: attempts × disc area / substrate area. The jamming limit is approached as τ^(−1/2).' },
  { key: 'stack', label: 'Layer stack (bottom to top)', type: 'table', columns: ['layer', 'formula', 'thickness_nm', 'density_kg_m3'], value: DEMO_STACK, group: 'Film', showIf: SH('film'), help: 'One layer per row from the substrate upwards: name, chemical formula, thickness in nm and film density in kg/m³. The default is an example stack.' },
  { key: 'rmsT', label: 'Surface roughness (rms)', unit: 'm', value: 2e-9, min: 1e-11, max: 1e-6, log: true, group: 'Film', showIf: SH('film'), help: 'Root-mean-square roughness of the top surface; the generated self-affine surface is scaled to this value.' },
  { key: 'hurst', label: 'Hurst exponent H', unit: '–', value: 0.7, min: 0.1, max: 1, group: 'Film', showIf: SH('film'), help: 'Roughness exponent of the self-affine surface: H near 1 is smooth at short range, small H is jagged. The local fractal dimension is 3 − H.' },
  { key: 'xi', label: 'Lateral correlation length', unit: 'm', value: 5e-8, min: 1e-9, max: 1e-5, log: true, group: 'Film', showIf: SH('film'), help: 'Distance beyond which heights become uncorrelated (roll-off of the power spectrum).' },
  { key: 'scan', label: 'Scan size', unit: 'm', value: 1e-6, min: 1e-8, max: 1e-4, log: true, group: 'Film', sens: false, showIf: SH('film'), help: 'Side of the square height map (as an AFM scan). It should span several correlation lengths.' },
  { key: 'islandN', label: 'Island number density', unit: '1/m²', value: 1e15, min: 1e13, max: 1e17, log: true, group: 'Film', showIf: SH('film'), help: 'Nucleation density for the island-growth (Volmer–Weber) coverage model; higher density closes the film at a smaller thickness.' },
  { key: 'featAR', label: 'Feature aspect ratio (depth/width)', unit: '–', value: 20, min: 1, max: 200, log: true, group: 'Film', showIf: SH('film'), help: 'Aspect ratio of the trench or via to be coated, used for the step-coverage estimates.' },
  { key: 'sticking', label: 'Reactive sticking coefficient', unit: '–', value: 0.01, min: 1e-5, max: 1, log: true, group: 'Film', showIf: SH('film'), help: 'Probability that a precursor molecule reacts on striking the wall (CVD-type first-order model). Lower values give more conformal coatings.' },
  { key: 'seed', label: 'Random seed', unit: '–', value: 1, min: 1, max: 9999, step: 1, group: 'Model', sens: false, showIf: SH(['population', 'rods', 'porous', 'aggregate', 'film']), help: 'Seed of the pseudo-random generator. Identical inputs and seed reproduce the same geometry exactly.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
function clean(inp) {
  const v = {};
  for (const i of inputs) {
    const t = i.type || 'number', x = inp ? inp[i.key] : undefined;
    if (t === 'number') { const q = Number(x); v[i.key] = Number.isFinite(q) ? min(i.max, max(i.min, q)) : i.value; }
    else if (t === 'select') v[i.key] = i.options.some(o => o.value === x) ? x : i.value;
    else v[i.key] = x === undefined || x === null ? i.value : String(x);
  }
  v.seed = round(v.seed);
  return v;
}

// ---------------------------------------------------------------- single nanostructure
const SHAPE_NAME = { sphere: 'sphere', spheroid: 'spheroid', rod_hex: 'hexagonal rod', rod_circ: 'circular rod', tube: 'tube', sheet: 'sheet', coreshell: 'core–shell sphere', cube: 'cube', octahedron: 'octahedron' };
const SHAPE_VARS = { sphere: ['d'], spheroid: ['d', 'L'], rod_hex: ['d', 'L'], rod_circ: ['d', 'L'], tube: ['d', 'L', 'tw'], sheet: ['d', 'tw'], coreshell: ['d', 'tw'], cube: ['d'], octahedron: ['d'] };
const SHAPE_EQ = {
  sphere: 'V = πd³/6,  A = πd²',
  spheroid: 'V = (π/6)d²L;  prolate: A = 2πa²[1 + (c/(a e)) asin e], e² = 1 − a²/c²;  oblate: A = 2πa²[1 + ((1 − e²)/e) atanh e], e² = 1 − c²/a²  (a = d/2, c = L/2)',
  rod_hex: 'V = (√3/2)d²L,  A = 2√3·d·L + √3·d²  (d across flats)',
  rod_circ: 'V = πd²L/4,  A = πdL + πd²/2',
  tube: 'V = π(d² − d_i²)L/4,  A = π(d + d_i)L + π(d² − d_i²)/2,  d_i = d − 2t',
  sheet: 'V = d²t,  A = 2d² + 4dt',
  coreshell: 'D = d + 2t,  V = πD³/6,  A = πD²,  m = ρ_c πd³/6 + ρ_s π(D³ − d³)/6',
  cube: 'V = d³,  A = 6d²',
  octahedron: 'V = (√2/3)d³,  A = 2√3·d²',
};
function shapeGeom(shape, d, L, tw) {
  let V, A, dims, Vcore = null;
  if (shape === 'sphere') { V = PI * d ** 3 / 6; A = PI * d * d; dims = [d, d, d]; }
  else if (shape === 'spheroid') {
    const a = d / 2, c = L / 2; V = 4 / 3 * PI * a * a * c; dims = [d, d, L];
    if (abs(c - a) < 1e-6 * a) A = 4 * PI * a * a;
    else if (c > a) { const e = sqrt(1 - a * a / (c * c)); A = 2 * PI * a * a * (1 + c / (a * e) * Math.asin(e)); }
    else { const e = sqrt(1 - c * c / (a * a)); A = 2 * PI * a * a * (1 + (1 - e * e) / e * Math.atanh(e)); }
  }
  else if (shape === 'rod_hex') { V = sqrt(3) / 2 * d * d * L; A = 2 * sqrt(3) * d * L + sqrt(3) * d * d; dims = [d, 2 * d / sqrt(3), L]; }
  else if (shape === 'rod_circ') { V = PI * d * d * L / 4; A = PI * d * L + PI * d * d / 2; dims = [d, d, L]; }
  else if (shape === 'tube') { const t = min(tw, 0.49 * d), di = d - 2 * t; V = PI * (d * d - di * di) * L / 4; A = PI * (d + di) * L + PI * (d * d - di * di) / 2; dims = [d, d, L]; }
  else if (shape === 'sheet') { V = d * d * tw; A = 2 * d * d + 4 * d * tw; dims = [d, d, tw]; }
  else if (shape === 'coreshell') { const D = d + 2 * tw; V = PI * D ** 3 / 6; A = PI * D * D; dims = [D, D, D]; Vcore = PI * d ** 3 / 6; }
  else if (shape === 'cube') { V = d ** 3; A = 6 * d * d; dims = [d, d, d]; }
  else { V = Math.SQRT2 / 3 * d ** 3; A = 2 * sqrt(3) * d * d; dims = [d * Math.SQRT2, d * Math.SQRT2, d * Math.SQRT2]; }
  return { V, A, dims, Vcore, aspect: max(...dims) / min(...dims) };
}
function shapeMesh(shape, d, L, tw) {   // nm
  const D = nm(d), Ln = nm(L), T = nm(tw);
  let m;
  if (shape === 'sphere') m = msh.sphere({ radius: D / 2, segments: 24 });
  else if (shape === 'spheroid') m = msh.transformMesh(msh.sphere({ radius: 1, segments: 24 }), { scale: [D / 2, D / 2, Ln / 2] });
  else if (shape === 'rod_hex') m = msh.hexPrism({ acrossFlats: D, height: Ln });
  else if (shape === 'rod_circ') m = msh.cylinder({ radius: D / 2, height: Ln, segments: 24 });
  else if (shape === 'tube') m = msh.tube({ outer: D / 2, inner: D / 2 - min(T, 0.49 * D), height: Ln, segments: 24 });
  else if (shape === 'sheet') m = msh.box({ size: [D, D, T] });
  else if (shape === 'coreshell') m = msh.mergeMeshes([msh.sphere({ radius: D / 2 + T, segments: 24 }), msh.sphere({ radius: D / 2, segments: 16 })], 'nm');
  else if (shape === 'cube') m = msh.box({ size: [D, D, D] });
  else { const r = D / Math.SQRT2; m = { kind: 'mesh', vertices: [[r, 0, 0], [-r, 0, 0], [0, r, 0], [0, -r, 0], [0, 0, r], [0, 0, -r]], faces: [[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4], [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]] }; }
  return { kind: 'mesh', vertices: m.vertices, faces: m.faces, units: 'nm' };
}
function runSingle(v) {
  const warnings = [], mat = material(v, warnings), sh = v.shape, vars = SHAPE_VARS[sh];
  const shell = sh === 'coreshell' ? matInfo(MAT[v.shellMat].f, MAT[v.shellMat].rho, MAT[v.shellMat].name, MAT[v.shellMat].x) : null;
  if (sh === 'tube' && v.tw > 0.49 * v.d) warnings.push(`Wall thickness limited to 0.49 d = ${F(nm(0.49 * v.d), 3)} nm so that the tube keeps an open bore.`);
  const massOf = g => shell ? mat.rho * g.Vcore + shell.rho * (g.V - g.Vcore) : mat.rho * g.V;
  const desc = (d, L, tw) => { const g = shapeGeom(sh, d, L, tw), m = massOf(g); return { ...g, mass: m, ssa: g.A / m }; };
  const g = desc(v.d, v.L, v.tw);
  const dV = cbrt(6 * g.V / PI), dS = sqrt(g.A / PI), d32 = 6 * g.V / g.A, psi = PI ** (1 / 3) * (6 * g.V) ** (2 / 3) / g.A;
  const nAtoms = shell ? g.Vcore / mat.vAt + (g.V - g.Vcore) / shell.vAt : g.V / mat.vAt;
  const surfMat = shell || mat, delta = surfMat.delta;
  const disp = min(1, g.A * delta / g.V);                                  // first-order: one atomic layer over the whole surface
  const fShell = d32 > 2 * delta ? 1 - (1 - 2 * delta / d32) ** 3 : 1;     // exact shell fraction of the Sauter-equivalent sphere
  // first-order tolerance propagation
  const tol = v.tolPct / 100, base = { d: v.d, L: v.L, tw: v.tw }, outs = [['Volume V', 'V'], ['Surface area A', 'A'], ['Specific surface area', 'ssa'], ['Aspect ratio', 'aspect']];
  const sens = outs.map(([, k]) => vars.map(x => { const h = 0.01, p = { ...base }, q = { ...base }; p[x] *= 1 + h; q[x] *= 1 - h; return (log(desc(p.d, p.L, p.tw)[k]) - log(desc(q.d, q.L, q.tw)[k])) / (log(1 + h) - log(1 - h)); }));
  const worst = sens.map(r => r.reduce((s, c) => s + abs(c), 0) * tol), rss = sens.map(r => sqrt(r.reduce((s, c) => s + c * c, 0)) * tol);
  const second = vars[1] || null, nmap = 21, data = new Float64Array(nmap * nmap);
  for (let j = 0; j < nmap; j++) for (let i = 0; i < nmap; i++) {
    const ex = tol * (2 * i / (nmap - 1) - 1), ey = tol * (2 * j / (nmap - 1) - 1), p = { ...base }; p.d *= 1 + ex; let s;
    if (second) { p[second] *= 1 + ey; s = desc(p.d, p.L, p.tw).ssa; } else s = desc(p.d, p.L, p.tw).ssa / (1 + ey);
    data[j * nmap + i] = fin((s / g.ssa - 1) * 100);
  }
  const VN = { d: 'd', L: 'L', tw: 't' };
  // size scaling
  const ks = N.logspace(1e-9 / v.d, 1e-6 / v.d, 60), sx = ks.map(k => nm(k * v.d)), sS = ks.map(k => g.ssa / k / 1e3), sD = ks.map(k => min(1, disp / k) * 100);
  // geometry: atoms for small particles with a known crystal, otherwise a surface mesh
  let geometry = null, atomsNote = null, exportGeo;
  const aShape = { sphere: 'sphere', cube: 'cube', rod_circ: 'rod', rod_hex: 'hexrod', octahedron: 'octahedron' }[sh];
  if (nAtoms <= 5000 && mat.x && (aShape || (sh === 'coreshell' && shell.x))) {
    try {
      const a = sh === 'coreshell' ? atomic.coreShell(mat.x, shell.x, { core_nm: nm(v.d), shell_nm: nm(v.tw), maxAtoms: 6000 })
        : atomic.nanoparticle(mat.x, { shape: aShape, diameter_nm: nm(v.d) * (sh === 'octahedron' ? Math.SQRT2 : 1), length_nm: nm(v.L), maxAtoms: 6000 });
      if (a.atoms.length >= 4) {
        geometry = { kind: 'atoms', atoms: a.atoms.map(q => ({ el: q.el, x: q.x, y: q.y, z: q.z })), cell: null, title: a.title };
        atomsNote = { n: a.atoms.length, truncated: !!a.truncated };
        exportGeo = { name: `${mat.f}_${sh}.xyz`, mime: 'chemical/x-xyz', content: atomic.toXYZ(geometry), description: 'Atomistic model carved from the bulk crystal (unrelaxed surfaces), XYZ format.' };
      }
    } catch { geometry = null; }
  }
  if (!geometry) { geometry = shapeMesh(sh, v.d, v.L, v.tw); exportGeo = { name: `${mat.f}_${sh}.obj`, mime: 'model/obj', content: toOBJ(geometry), description: 'Surface mesh of the nanostructure in nm (Wavefront OBJ).' }; }
  let meshCheck = null;
  if (geometry.kind === 'mesh' && sh !== 'coreshell') { const st = msh.meshStats(geometry); meshCheck = { vol: st.volume * 1e-27, area: st.surfaceArea * 1e-27 * 1e9, water: st.watertight }; }
  const ssa_g = g.ssa / 1e3, what = `${F(nm(v.d), 3)} nm ${mat.pretty}${shell ? '@' + shell.pretty : ''} ${SHAPE_NAME[sh]}`;
  const checks = [
    chk('Isoperimetric inequality (sphericity ≤ 1)', psi <= 1 + 1e-9, `sphericity ψ = π^(1/3)(6V)^(2/3)/A = ${F(psi, 4)}`),
    chk('Equivalent-diameter ordering d₃₂ ≤ d_V ≤ d_S', d32 <= dV * (1 + 1e-9) && dV <= dS * (1 + 1e-9), `d₃₂ = ${F(nm(d32), 4)} nm, d_V = ${F(nm(dV), 4)} nm, d_S = ${F(nm(dS), 4)} nm`),
    chk('Continuum description meaningful (≥ 5 atomic spacings across the smallest dimension)', min(...g.dims) >= 5 * delta, `smallest dimension ${F(nm(min(...g.dims)), 3)} nm against atomic spacing ${F(delta * 1e10, 3)} Å`),
  ];
  if (meshCheck) checks.push(chk('Mesh volume against the analytical volume', abs(meshCheck.vol / g.V - 1) < 0.05, `faceted mesh ${F(meshCheck.vol, 4)} m³ against ${F(g.V, 4)} m³ (${F((meshCheck.vol / g.V - 1) * 100, 2)} %; faceting of curved surfaces); watertight: ${meshCheck.water}`));
  if (atomsNote) checks.push(chk('Atomistic count against the continuum estimate', abs(atomsNote.n / nAtoms - 1) < 0.25 || atomsNote.truncated, `${atomsNote.n} atoms in the carved model against ${F(nAtoms, 4)} from V·ρ·N_A·n/M${atomsNote.truncated ? ' (model truncated to the atom cap)' : ''}`));
  if (min(...g.dims) < 5 * delta) warnings.push('The smallest dimension is below five atomic spacings: continuum volume and surface area are only indicative; use the atomistic model.');
  if (disp >= 1) warnings.push('Every atom lies within one atomic layer of the surface: the surface-atom fraction is capped at 1.');
  return {
    summary: `A ${what} has volume ${F(g.V, 3)} m³ and surface area ${F(g.A, 3)} m², giving a specific surface area of ${F(ssa_g, 4)} m²/g and a Sauter diameter of ${F(nm(d32), 4)} nm. It holds about ${F(nAtoms, 3)} atoms, of which roughly ${F(disp * 100, 3)} % sit in the outermost atomic layer. A ±${F(v.tolPct, 3)} % dimensional tolerance changes the specific surface area by up to ±${F(worst[2] * 100, 3)} % (worst case, first order).`,
    kpis: [
      kpi('ssa', 'Specific surface area', ssa_g, 'm²/g', `Geometric surface area divided by mass, A/(ρV) = ${F(ssa_g, 4)} m²/g for a smooth, fully accessible, non-porous particle. Measured BET values are higher when the surface is rough or porous and lower when particles are necked or agglomerated.`),
      kpi('V', 'Volume', g.V, 'm³', `Solid volume of one ${SHAPE_NAME[sh]} from the closed-form shape formula; it sets the mass (${F(g.mass, 3)} kg per particle) and the atom count.`),
      kpi('A', 'Surface area', g.A, 'm²', `External surface area of one particle including end faces (and the bore of a tube) from the closed-form shape formula.`),
      kpi('SV', 'Surface-to-volume ratio', g.A / g.V, '1/m', `A/V = ${F(g.A / g.V, 4)} 1/m; it scales as 1/size, so halving every dimension doubles it. This ratio governs dissolution rate, catalytic activity per volume and surface-energy contributions.`),
      kpi('d32', 'Sauter diameter d₃₂', d32, 'm', `6V/A = ${F(nm(d32), 4)} nm: the diameter of the sphere with the same surface-to-volume ratio. It is the size that BET-derived "particle size" reports.`),
      kpi('dV', 'Volume-equivalent sphere diameter', dV, 'm', `(6V/π)^(1/3) = ${F(nm(dV), 4)} nm: the sphere of equal volume (and mass); relevant to sedimentation and counting by mass.`),
      kpi('dS', 'Surface-equivalent sphere diameter', dS, 'm', `(A/π)^(1/2) = ${F(nm(dS), 4)} nm: the sphere of equal surface area; it is never smaller than the volume-equivalent diameter.`),
      kpi('aspect', 'Aspect ratio', g.aspect, '–', `Largest over smallest principal dimension = ${F(g.aspect, 4)}. Values far from 1 indicate 1-D (rod, wire, tube) or 2-D (sheet) character and anisotropic properties.`),
      kpi('psi', 'Sphericity', psi, '–', `Wadell sphericity, the surface area of the volume-equivalent sphere divided by the actual area = ${F(psi, 4)}; 1 for a sphere and lower for elongated or flat shapes.`),
      kpi('nAtoms', 'Number of atoms', nAtoms, '–', `V·ρ·N_A·n/M with the bulk density and molar mass = ${F(nAtoms, 4)} atoms. Bulk density is assumed down to the nanoscale; lattice relaxation is neglected.`),
      kpi('disp', 'Surface-atom fraction (dispersion)', disp, '–', `Fraction of atoms in the outermost atomic layer, A·δ/V with δ = v_a^(1/3) = ${F(delta * 1e10, 3)} Å, = ${F(disp * 100, 3)} %. In catalysis this is the dispersion; the exact shell fraction of the Sauter-equivalent sphere is ${F(fShell * 100, 3)} %.`, null),
      kpi('ssaTol', 'SSA deviation at tolerance (worst case)', worst[2] * 100, '%', `First-order worst-case change of the specific surface area when every dimension deviates by ±${F(v.tolPct, 3)} %: Σ|∂ln SSA/∂ln xᵢ|·tol. The root-sum-square value for independent deviations is ±${F(rss[2] * 100, 3)} %.`, worst[2] < 0.1),
    ],
    series: [
      { key: 'ssa_size', title: 'Specific surface area vs size (all dimensions scaled together)', xlabel: 'characteristic size d (nm)', ylabel: 'specific surface area (m²/g)', logx: true, logy: true, lines: [{ name: 'SSA', x: sx, y: sS }], explain: 'Geometric similarity gives SSA ∝ 1/size: a straight line of slope −1 on logarithmic axes. Read off how much surface a given size class contributes.' },
      { key: 'disp_size', title: 'Surface-atom fraction vs size', xlabel: 'characteristic size d (nm)', ylabel: 'surface atoms (%)', logx: true, logy: false, lines: [{ name: 'first-order shell model', x: sx, y: sD }], explain: 'Share of atoms in the outermost atomic layer. Below a few nanometres most atoms are surface atoms and size-dependent properties (melting point, lattice parameter, reactivity) become pronounced.' },
    ],
    maps: [{ key: 'tolmap', title: 'Tolerance map: change of specific surface area', nx: nmap, ny: nmap, data, x0: -v.tolPct, x1: v.tolPct, y0: -v.tolPct, y1: v.tolPct, xlabel: 'deviation of d (%)', ylabel: second ? `deviation of ${VN[second]} (%)` : 'deviation of density (%)', unit: '%', diverging: true, explain: `Relative change of the specific surface area over the tolerance window, evaluated exactly from the shape formulae${second ? '' : ' (the second axis is the bulk-density deviation because this shape has a single dimension)'}. The gradient direction shows which dimension must be controlled most tightly.` }],
    tables: [{ key: 'tol', title: 'First-order tolerance propagation', columns: ['Output', ...vars.map(x => `∂ln/∂ln ${VN[x]}`), 'Worst case (± %)', 'Root-sum-square (± %)'], rows: outs.map(([n], i) => [n, ...sens[i].map(s => +s.toFixed(3)), +(worst[i] * 100).toFixed(3), +(rss[i] * 100).toFixed(3)]), explain: 'Logarithmic sensitivities from central finite differences of the closed-form shape functions, combined with the dimensional tolerance as a worst-case sum and as a root-sum-square for independent deviations.' }],
    equations: [
      matEq(mat), ...(shell ? [matEq(shell)] : []),
      eqn('formula', `Volume and surface area — ${SHAPE_NAME[sh]}`, SHAPE_EQ[sh]),
      eqn('formula', 'Specific surface area and equivalent diameters', 'SSA = A/(ρV);  d_V = (6V/π)^(1/3);  d_S = (A/π)^(1/2);  d₃₂ = 6V/A;  ψ = π^(1/3)(6V)^(2/3)/A'),
      eqn('formula', 'Atom count and dispersion', 'N = V ρ N_A n / M;  v_a = M/(ρ N_A n);  δ = v_a^(1/3);  D ≈ A δ / V'),
      eqn('governing', 'First-order tolerance propagation', 'Δy/y ≈ Σᵢ (∂ln y/∂ln xᵢ)(Δxᵢ/xᵢ);  worst case Σ|sᵢ|·tol;  independent: (Σ sᵢ²)^(1/2)·tol'),
    ],
    geometry, checks, warnings,
    recommendations: [
      rec('high', `To raise the specific surface area above ${F(ssa_g, 3)} m²/g, reduce the smallest dimension: SSA scales as 1/size, so ${F(nm(v.d) / 2, 3)} nm would give about ${F(ssa_g * 2, 3)} m²/g for the same shape.`, 'Geometric scaling SSA ∝ 1/d.'),
      rec('medium', `Control ${VN[vars[N.argmax(sens[2].map(abs))]]} most tightly: it has the largest logarithmic sensitivity (${F(max(...sens[2].map(abs)), 3)}) on the specific surface area.`, 'Tolerance sensitivity table.'),
      rec('medium', `Compare the geometric value with a BET measurement: a measured SSA well below ${F(ssa_g, 3)} m²/g indicates agglomeration or necking, a higher one indicates surface roughness or internal porosity.`, 'Smooth, fully accessible surface assumed in the geometric model.'),
      rec('low', disp > 0.2 ? `With ${F(disp * 100, 3)} % of atoms at the surface, less active material is needed per unit of surface function; dose by surface area rather than by mass to avoid over-use of ${mat.pretty}.` : `Only ${F(disp * 100, 3)} % of the atoms are at the surface; for surface-driven functions a smaller or higher-aspect-ratio morphology would use ${mat.pretty} more efficiently and reduce raw-material demand.`, 'Surface-atom fraction; material efficiency.', true),
      rec('low', 'Nanoscale powders with high specific surface area are respirable and can be pyrophoric or dust-explosive: handle as suspensions or in enclosed transfer wherever possible.', 'General nanomaterial handling practice for high surface-to-volume solids.', true),
    ],
    publish: { 'geom.diameter': v.d, 'geom.length': g.dims[2], 'geom.aspect': g.aspect, 'geom.ssa': g.ssa, 'geom.thickness': sh === 'sheet' || sh === 'tube' || sh === 'coreshell' ? v.tw : min(...g.dims), 'geom.volume': g.V, 'geom.surface_area': g.A, 'geom.d32': d32, 'prop.density': mat.rho, 'prop.molar_mass': mat.M, 'material.formula': mat.f },
    exports: [exportGeo, { name: 'descriptors.json', mime: 'application/json', content: JSON.stringify({ material: mat.f, shape: sh, d_m: v.d, L_m: v.L, t_m: v.tw, volume_m3: g.V, area_m2: g.A, ssa_m2_per_g: ssa_g, d_volume_m: dV, d_surface_m: dS, d_sauter_m: d32, aspect: g.aspect, sphericity: psi, atoms: nAtoms, dispersion: disp }, null, 1), description: 'Shape descriptors of this nanostructure.' }],
    numerics: { method: 'Closed-form shape geometry; central finite differences for tolerance sensitivities', dof: vars.length, steps: nmap * nmap, iterations: 1, residual: 0, notes: `Geometry: ${geometry.kind === 'atoms' ? 'atomistic model carved from the bulk crystal' : 'faceted surface mesh'}. Bulk density assumed at all sizes.` },
  };
}

// ---------------------------------------------------------------- particle population
const lnPdf = (x, mu, s) => x > 0 ? exp(-((log(x) - mu) ** 2) / (2 * s * s)) / (x * s * sqrt(2 * PI)) : 0;
const DISTS = {
  lognormal: { name: 'log-normal', cdf: (x, p) => x > 0 ? N.normCdf(log(x), p[0], p[1]) : 0, pdf: (x, p) => lnPdf(x, p[0], p[1]), text: p => `median ${F(exp(p[0]), 4)} nm, σ_g ${F(exp(p[1]), 4)}` },
  normal: { name: 'normal', cdf: (x, p) => N.normCdf(x, p[0], p[1]), pdf: (x, p) => exp(-((x - p[0]) ** 2) / (2 * p[1] * p[1])) / (p[1] * sqrt(2 * PI)), text: p => `mean ${F(p[0], 4)} nm, σ ${F(p[1], 4)} nm` },
  weibull: { name: 'Weibull', cdf: (x, p) => x > 0 ? 1 - exp(-((x / p[0]) ** p[1])) : 0, pdf: (x, p) => x > 0 ? p[1] / p[0] * (x / p[0]) ** (p[1] - 1) * exp(-((x / p[0]) ** p[1])) : 0, text: p => `scale ${F(p[0], 4)} nm, shape ${F(p[1], 4)}` },
};
function fitDists(s) {   // s sorted ascending, all > 0
  const n = s.length, ln = s.map(log), mu = N.mean(ln), sg = max(1e-6, sqrt(N.mean(ln.map(q => (q - mu) ** 2)))), m = N.mean(s), sd = max(1e-9 * m, sqrt(N.mean(s.map(q => (q - m) ** 2))));
  let k; const g = kk => { let a = 0, b = 0; const xm = s[n - 1]; for (const x of s) { const w = (x / xm) ** kk; a += w * log(x); b += w; } return a / b - 1 / kk - mu; };
  try { k = g(0.05) * g(200) < 0 ? N.brent(g, 0.05, 200, 1e-10) : N.weibullFit(s).shape; } catch { k = N.weibullFit(s).shape; }
  if (!(k > 0) || !Number.isFinite(k)) k = 1.2 * m / sd;
  const lam = N.mean(s.map(x => (x / s[n - 1]) ** k)) ** (1 / k) * s[n - 1];
  const P = { lognormal: [mu, sg], normal: [m, sd], weibull: [lam, k] }, out = {};
  for (const key in P) {
    const D = DISTS[key], p = P[key]; let ks = 0, ll = 0;
    for (let i = 0; i < n; i++) { const c = D.cdf(s[i], p); ks = max(ks, abs(c - i / n), abs(c - (i + 1) / n)); ll += log(max(1e-300, D.pdf(s[i], p))); }
    out[key] = { p, ks, aic: 4 - 2 * ll };
  }
  return out;
}
function popStats(s) {   // sorted ascending
  const n = s.length; let m1 = 0, m2 = 0, m3 = 0, m4 = 0; for (const x of s) { m1 += x; m2 += x * x; m3 += x ** 3; m4 += x ** 4; }
  const q = p => { const t = (n - 1) * p, i = floor(t); return s[i] + (s[min(i + 1, n - 1)] - s[i]) * (t - i); };
  let acc = 0, dv50 = s[n - 1], dv10 = null, dv90 = s[n - 1]; for (const x of s) { acc += x ** 3; if (dv10 === null && acc >= 0.1 * m3) dv10 = x; if (acc >= 0.5 * m3 && dv50 === s[n - 1]) dv50 = x; if (acc >= 0.9 * m3) { dv90 = x; break; } }
  const mean = m1 / n, sd = sqrt(max(0, m2 / n - mean * mean));
  return { n, mean, sd, cv: sd / mean, d32: m3 / m2, d43: m4 / m3, d10: q(0.1), d50: q(0.5), d90: q(0.9), span: (q(0.9) - q(0.1)) / q(0.5), pdi: (sd / mean) ** 2, dv10: dv10 ?? s[0], dv50, dv90 };
}
function runPopulation(v) {
  const warnings = [], mat = material(v, warnings), rand = N.rng(v.seed), fit = v.popSource === 'fit';
  let s = [], src;
  if (fit) {
    const rows = parseRows(v.sizes);
    for (const r of rows) { const x = Number(r[r.length - 1]); if (Number.isFinite(x) && x > 0) s.push(x); }
    if (String(v.sizes).trim() === DEMO_SIZES) warnings.push('The diameter table holds synthetic demonstration data — replace with measurements.');
    if (s.length < 8) { warnings.push(`Only ${s.length} valid positive diameters were found in the table (at least 8 are needed); the synthetic demonstration data were analysed instead.`); s = parseRows(DEMO_SIZES).map(r => Number(r[0])).filter(x => x > 0); }
    if (s.length > 20000) { s = s.slice(0, 20000); warnings.push('Only the first 20 000 diameters were used.'); }
    src = `${s.length} measured diameters`;
  } else {
    const dn = nm(v.d), n = round(v.nSample);
    if (v.dist === 'lognormal') for (let i = 0; i < n; i++) s.push(rand.lognormal(log(dn), log(v.gsd)));
    else if (v.dist === 'normal') { const sd = dn * v.cvPct / 100; while (s.length < n) { const x = rand.normal(dn, sd); if (x > 0.02 * dn) s.push(x); } }
    else { const lam = dn / log(2) ** (1 / v.wk); for (let i = 0; i < n; i++) s.push(lam * (-log(1 - rand() * 0.999999999)) ** (1 / v.wk)); }
    src = `${n} particles sampled from a ${DISTS[v.dist].name} distribution`;
  }
  s.sort((a, b) => a - b);
  const st = popStats(s), fits = fitDists(s), best = Object.keys(fits).sort((a, b) => fits[a].ks - fits[b].ks)[0], sel = fits[v.dist], D = DISTS[v.dist];
  const ksCrit = 1.36 / sqrt(st.n), ssa = 6 / (mat.rho * st.d32 * 1e-9), ssaMean = 6 / (mat.rho * st.mean * 1e-9);
  // histograms (number- and volume-weighted fractions) and cumulative curves
  const lo = s[0], hi = s[s.length - 1], bins = min(40, max(8, round(sqrt(st.n)))), w = (hi - lo) / bins || 1, hn = new Array(bins).fill(0), hv = new Array(bins).fill(0); let v3 = 0;
  for (const x of s) { const b = min(bins - 1, floor((x - lo) / w)); hn[b]++; hv[b] += x ** 3; v3 += x ** 3; }
  const hx = hn.map((_, i) => lo + (i + 0.5) * w), xs = N.linspace(lo, hi, 120), step = max(1, floor(st.n / 200)), cx = [], cn = [], cv = []; let acc = 0;
  s.forEach((x, i) => { acc += x ** 3; if (i % step === 0 || i === st.n - 1) { cx.push(x); cn.push((i + 1) / st.n * 100); cv.push(acc / v3 * 100); } });
  // Hatch–Choate check for the generated log-normal
  const checks = [
    chk('Moment-mean ordering D[1,0] ≤ D[3,2] ≤ D[4,3]', st.mean <= st.d32 * (1 + 1e-12) && st.d32 <= st.d43 * (1 + 1e-12), `D[1,0] = ${F(st.mean, 4)} nm, D[3,2] = ${F(st.d32, 4)} nm, D[4,3] = ${F(st.d43, 4)} nm`),
    chk(`Kolmogorov–Smirnov test of the ${D.name} fit (5 % level)`, sel.ks < ksCrit, `D = ${F(sel.ks, 3)} against the critical value 1.36/√n = ${F(ksCrit, 3)} (parameters estimated from the same data, so the test is conservative)`),
  ];
  if (!fit && v.dist === 'lognormal') { const l2 = log(v.gsd) ** 2, hc = nm(v.d) * exp(2.5 * l2); checks.push(chk('Sampled Sauter diameter against the Hatch–Choate relation', abs(st.d32 / hc - 1) < 0.1, `sampled D[3,2] = ${F(st.d32, 4)} nm against d_med·exp(2.5 ln²σ_g) = ${F(hc, 4)} nm (difference ${F((st.d32 / hc - 1) * 100, 2)} %; sampling error grows with σ_g)`)); }
  if (st.n < 200) warnings.push(`Only ${st.n} particles: volume-weighted statistics (D[4,3], Dv90) are dominated by the few largest particles and carry a large sampling error.`);
  if (st.d43 / st.mean > 2) warnings.push('Broad distribution: the volume-weighted mean exceeds twice the number mean, so number-based (TEM) and volume- or intensity-based (laser diffraction, DLS) sizes will disagree strongly.');
  // representative subset drawn to scale
  const pick = Array.from({ length: 36 }, (_, i) => s[min(st.n - 1, floor((i + 0.5) / 36 * st.n))]), pitch = pick[35] * 1.15;
  const mesh = msh.mergeMeshes(pick.map((dd, i) => msh.sphere({ radius: dd / 2, segments: 10, origin: [(i % 6) * pitch, floor(i / 6) * pitch, dd / 2] })), 'nm');
  const fitTable = Object.keys(fits).map(k => [DISTS[k].name, DISTS[k].text(fits[k].p), +fits[k].ks.toFixed(4), fits[k].ks < ksCrit ? 'not rejected' : 'rejected', +fits[k].aic.toFixed(1)]);
  return {
    summary: `${src} of ${mat.pretty}: number mean ${F(st.mean, 4)} nm, median ${F(st.d50, 4)} nm, Sauter mean D[3,2] ${F(st.d32, 4)} nm and volume-weighted mean D[4,3] ${F(st.d43, 4)} nm (span ${F(st.span, 3)}, polydispersity index ${F(st.pdi, 3)}). The population specific surface area is ${F(ssa / 1e3, 4)} m²/g. The best-fitting model by the Kolmogorov–Smirnov statistic is ${DISTS[best].name} (${DISTS[best].text(fits[best].p)}, D = ${F(fits[best].ks, 3)}).`,
    kpis: [
      kpi('d10m', 'Number mean D[1,0]', st.mean * 1e-9, 'm', `Arithmetic mean diameter Σd/N = ${F(st.mean, 4)} nm; this is what counting methods (TEM, SEM) report.`),
      kpi('d32', 'Sauter mean D[3,2]', st.d32 * 1e-9, 'm', `Σd³/Σd² = ${F(st.d32, 4)} nm: the surface-weighted mean that fixes the specific surface area of the whole population.`),
      kpi('d43', 'Volume-weighted mean D[4,3]', st.d43 * 1e-9, 'm', `Σd⁴/Σd³ = ${F(st.d43, 4)} nm (De Brouckere mean); it is what laser diffraction and mass-based methods emphasise and is sensitive to a few large particles.`),
      kpi('p10', 'Number percentile d10', st.d10 * 1e-9, 'm', `Ten per cent of the particles by number are smaller than ${F(st.d10, 4)} nm (volume-based Dv10 = ${F(st.dv10, 4)} nm).`),
      kpi('p50', 'Number median d50', st.d50 * 1e-9, 'm', `Half the particles by number are smaller than ${F(st.d50, 4)} nm; the volume-based median Dv50 is ${F(st.dv50, 4)} nm.`),
      kpi('p90', 'Number percentile d90', st.d90 * 1e-9, 'm', `Ninety per cent of the particles by number are smaller than ${F(st.d90, 4)} nm (volume-based Dv90 = ${F(st.dv90, 4)} nm).`),
      kpi('span', 'Span (d90 − d10)/d50', st.span, '–', `Width of the number distribution relative to its median = ${F(st.span, 3)}. Values below about 0.5 are narrow; above 1.5 broad.`, st.span < 1),
      kpi('pdi', 'Polydispersity index (σ/mean)²', st.pdi, '–', `Squared coefficient of variation of the number distribution = ${F(st.pdi, 3)} (CV = ${F(st.cv * 100, 3)} %). Below 0.05 is near-monodisperse; the DLS cumulant index is defined on the intensity distribution and is not identical.`, st.pdi < 0.1),
      kpi('ssa', 'Population specific surface area', ssa / 1e3, 'm²/g', `6/(ρ·D[3,2]) = ${F(ssa / 1e3, 4)} m²/g for dense spheres. Using the number mean instead would give ${F(ssaMean / 1e3, 4)} m²/g, an overestimate for any polydisperse powder.`),
      kpi('ks', `KS statistic of the ${D.name} fit`, sel.ks, '–', `Largest distance between the empirical and fitted cumulative distributions = ${F(sel.ks, 3)}; the 5 % critical value for ${st.n} particles is ${F(ksCrit, 3)}.`, sel.ks < ksCrit),
    ],
    series: [
      { key: 'hist', title: 'Size distribution: number- and volume-weighted', xlabel: 'diameter (nm)', ylabel: 'fraction per bin (%)', lines: [{ name: 'number', x: hx, y: hn.map(c => c / st.n * 100) }, { name: 'volume', x: hx, y: hv.map(c => c / v3 * 100) }, { name: `${D.name} fit (number)`, x: xs, y: xs.map(x => fin(D.pdf(x, sel.p)) * w * 100) }], explain: 'The volume-weighted histogram is shifted to larger sizes because volume scales as d³: a small number of coarse particles carries much of the mass. The fitted density is overlaid on the number histogram.' },
      { key: 'cum', title: 'Cumulative undersize curves', xlabel: 'diameter (nm)', ylabel: 'cumulative undersize (%)', lines: [{ name: 'number', x: cx, y: cn }, { name: 'volume', x: cx, y: cv }, { name: `${D.name} fit (number)`, x: xs, y: xs.map(x => fin(D.cdf(x, sel.p)) * 100) }], explain: 'Read d10, d50 and d90 where the curves cross 10, 50 and 90 %. The vertical gap between the empirical number curve and the fit is the Kolmogorov–Smirnov statistic.' },
    ],
    tables: [
      { key: 'fits', title: 'Distribution fits (maximum likelihood)', columns: ['Model', 'Parameters', 'KS statistic D', `KS at 5 % (critical ${F(ksCrit, 3)})`, 'AIC'], rows: fitTable, explain: 'Log-normal and normal parameters are the maximum-likelihood moments of ln d and d; the Weibull shape solves the likelihood equation by root finding. A lower KS statistic and a lower AIC indicate the better description.' },
      { key: 'means', title: 'Moment means and percentiles', columns: ['Quantity', 'Value (nm)'], rows: [['D[1,0] number mean', st.mean], ['D[3,2] Sauter mean', st.d32], ['D[4,3] volume-weighted mean', st.d43], ['d10 (number)', st.d10], ['d50 (number)', st.d50], ['d90 (number)', st.d90], ['Dv10 (volume)', st.dv10], ['Dv50 (volume)', st.dv50], ['Dv90 (volume)', st.dv90], ['standard deviation', st.sd]].map(r => [r[0], +r[1].toPrecision(5)]), explain: 'Number-based percentiles come from counting; volume-based percentiles weight each particle by d³, as in laser diffraction or sieving by mass.' },
    ],
    equations: [
      matEq(mat),
      eqn('formula', 'Moment means', 'D[p,q] = (Σ dᵢᵖ / Σ dᵢ^q)^(1/(p−q));  D[1,0] number mean,  D[3,2] Sauter mean,  D[4,3] volume-weighted mean'),
      eqn('formula', 'Population specific surface area (dense spheres)', 'SSA = 6 / (ρ · D[3,2])'),
      eqn('constitutive', 'Distributions', 'log-normal: f(d) = exp(−(ln d − μ)²/(2σ²)) / (d σ √(2π)), σ_g = e^σ;  Weibull: F(d) = 1 − exp(−(d/λ)^k);  normal: F(d) = Φ((d − m)/s)'),
      eqn('formula', 'Hatch–Choate relations (log-normal)', 'D[3,2] = d_med·exp(2.5 ln²σ_g);  D[4,3] = d_med·exp(3.5 ln²σ_g);  D[1,0] = d_med·exp(0.5 ln²σ_g)'),
      eqn('formula', 'Kolmogorov–Smirnov statistic', 'D = supₓ |F_n(x) − F(x)|;  critical value at 5 %: 1.36/√n'),
    ],
    geometry: { kind: 'mesh', vertices: mesh.vertices, faces: mesh.faces, units: 'nm' },
    checks, warnings,
    recommendations: [
      rec('high', `Quote D[3,2] = ${F(st.d32, 3)} nm when relating size to surface area or reactivity, and D[4,3] = ${F(st.d43, 3)} nm when comparing with laser-diffraction or mass-based data; the number mean ${F(st.mean, 3)} nm alone under-represents the coarse tail.`, 'Moment means of the distribution.'),
      rec(st.span > 1 ? 'high' : 'medium', st.span > 1 ? `The span of ${F(st.span, 3)} is broad: separate nucleation from growth (rapid mixing, seeded growth) or classify the product (centrifugation, size-selective precipitation) to narrow it.` : `The span of ${F(st.span, 3)} is already narrow; keep mixing time and temperature history constant between batches to preserve it.`, 'Span and polydispersity index of the population.'),
      rec('medium', `Count at least ${max(300, round(4 * (st.cv / 0.05) ** 2))} particles if the number mean must be known to ±5 % at 95 % confidence (CV = ${F(st.cv * 100, 3)} %); more are needed for volume-weighted statistics.`, 'Standard error of the mean, 1.96·CV/√n.'),
      rec('low', `A narrower distribution removes the oversize fraction that contributes mass but little surface: here the coarsest 10 % by number hold ${F((1 - N.interp1(cx, cv, st.d90) / 100) * 100, 3)} % of the volume. Tightening the distribution therefore saves ${mat.pretty} for the same surface function.`, 'Volume-weighted cumulative curve; material efficiency.', true),
    ],
    publish: { 'geom.diameter': st.d50 * 1e-9, 'geom.ssa': ssa, 'geom.d32': st.d32 * 1e-9, 'geom.d43': st.d43 * 1e-9, 'geom.span': st.span, 'geom.pdi': st.pdi, 'geom.aspect': 1, 'prop.density': mat.rho, 'material.formula': mat.f },
    exports: [
      { name: 'size_distribution.csv', mime: 'text/csv', content: 'diameter_nm,number_fraction,volume_fraction\n' + hx.map((x, i) => `${x.toPrecision(6)},${(hn[i] / st.n).toPrecision(5)},${(hv[i] / v3).toPrecision(5)}`).join('\n') + '\n', description: 'Number- and volume-weighted histogram of the population.' },
      { name: 'population_statistics.json', mime: 'application/json', content: JSON.stringify({ material: mat.f, n: st.n, statistics_nm: st, fits: Object.fromEntries(Object.keys(fits).map(k => [k, { parameters: fits[k].p, ks: fits[k].ks, aic: fits[k].aic }])), ssa_m2_per_g: ssa / 1e3 }, null, 1), description: 'Moment means, percentiles and fitted distribution parameters.' },
    ],
    numerics: { method: fit ? 'Maximum-likelihood fits (closed form for log-normal and normal, Brent root finding for the Weibull shape); empirical moments' : 'Seeded inverse-transform / Box–Muller sampling; empirical moments; maximum-likelihood re-fit', dof: st.n, steps: bins, iterations: 3, residual: sel.ks, notes: 'Percentiles by linear interpolation of the sorted sample. Spherical, dense particles assumed for the specific surface area.' },
  };
}

// ---------------------------------------------------------------- nanorod array
function runRods(v) {
  const warnings = [], mat = material(v, warnings), hex = v.section === 'hex';
  const side = min(1e6, max(60, sqrt(400 / v.nDens) * 1e9));   // nm; about 400 rods in view
  const arr = msh.nanorodArray({ density_per_um2: v.nDens * 1e-12, length_nm: nm(v.L), diameter_nm: nm(v.d), tiltSpread_deg: v.tilt, area_nm: [side, side], seed: v.seed, hex, lengthCV: v.cvL / 100, diameterCV: v.cvD / 100, arrangement: v.arrangement });
  let rods = arr.rods.map(r => ({ x: r.x, y: r.y, length: r.length, diameter: r.diameter, tilt: min(89, r.tilt), azimuth: r.azimuth }));
  if (!rods.length) { rods = [{ x: side / 2, y: side / 2, length: nm(v.L), diameter: nm(v.d), tilt: 0, azimuth: 0 }]; warnings.push('The requested density gives less than one rod in the modelled area; a single rod is shown.'); }
  if (arr.note) warnings.push(arr.note + '.'); if (arr.capped) warnings.push(arr.capped + '.');
  const n = rods.length, area = side * side, kA = hex ? sqrt(3) / 2 : PI / 4, kP = hex ? 2 * sqrt(3) : PI;   // cross-section = kA d², perimeter = kP d (hexagon: d across flats)
  const Ls = rods.map(r => r.length), Ds = rods.map(r => r.diameter), Ts = rods.map(r => r.tilt);
  let cover = 0, lateral = 0, vol = 0, p2 = 0, proj = 0;
  for (const r of rods) { const a = kA * r.diameter ** 2, c = cos(r.tilt * PI / 180); cover += a / max(0.02, c); lateral += kP * r.diameter * r.length; vol += a * r.length; p2 += (3 * c * c - 1) / 2; proj += a + r.diameter * r.length * sin(r.tilt * PI / 180); }
  const dens = n / (area * 1e-18), coverage = cover / area, sef = 1 + lateral / area, P2 = p2 / n, mL = N.mean(Ls), mD = N.mean(Ds);
  const massA = mat.rho * vol * 1e-27 / (area * 1e-18), ssa = (lateral + n * 0) * 1e-18 / (mat.rho * vol * 1e-27) + 0, ssaTop = (lateral + cover) * 1e-18 / (mat.rho * vol * 1e-27);
  // nearest-neighbour spacing (periodic)
  const nn = rods.map((a, i) => { let m = Infinity; for (let j = 0; j < n; j++) { if (j === i) continue; let dx = abs(a.x - rods[j].x), dy = abs(a.y - rods[j].y); if (dx > side / 2) dx = side - dx; if (dy > side / 2) dy = side - dy; const q = dx * dx + dy * dy; if (q < m) m = q; } return sqrt(m); }).filter(Number.isFinite);
  const nnMean = nn.length ? N.mean(nn) : side, pitch = sqrt(area / n), gap = nnMean - mD, nnPoisson = 0.5 * pitch, nnHex = pitch * sqrt(2 / sqrt(3));
  // volume fraction against height
  const hmax = N.maxOf(rods.map(r => r.length * cos(r.tilt * PI / 180))), zs = N.linspace(0, hmax, 50), phiZ = zs.map(z => { let s = 0; for (const r of rods) { const c = max(0.02, cos(r.tilt * PI / 180)); if (z <= r.length * c) s += kA * r.diameter ** 2 / c; } return min(1, s / area); });
  const phiMean = N.trapz(zs, phiZ) / (hmax || 1), porosity = 1 - phiMean;
  const open = exp(-min(50, proj / area)), alphaC = Math.atan2(max(0, gap), mL) * 180 / PI, hermans = P2;
  const hT = N.histogram(Ts, 15, 0, max(1, N.maxOf(Ts))), hS = N.histogram(nn.length ? nn : [side], 15);
  const mesh = msh.rodsToMesh({ kind: 'rods', rods, area: [side, side], units: 'nm', section: hex ? 'hex' : 'round' }, { maxRods: 300 });
  const checks = [
    chk('Areal coverage below unity (rods do not fill the plane)', coverage < 1, `coverage = ${F(coverage * 100, 3)} % of the substrate at the base`),
    chk('Positive mean gap between neighbouring rods', gap > 0, `mean nearest-neighbour centre spacing ${F(nnMean, 4)} nm against mean diameter ${F(mD, 4)} nm`),
    chk('Generated number density against the request', abs(dens / v.nDens - 1) < 0.15, `${F(dens, 3)} 1/m² generated (${n} rods on ${F(side, 4)} nm × ${F(side, 4)} nm) against ${F(v.nDens, 3)} 1/m² requested`),
    chk('Orientation order parameter within its bounds', P2 >= -0.5 - 1e-9 && P2 <= 1 + 1e-9, `⟨P₂⟩ = ${F(P2, 4)} (1 = perfectly vertical, 0 = isotropic)`),
  ];
  if (coverage > 0.6) warnings.push(`Areal coverage of ${F(coverage * 100, 3)} % exceeds what separate rods can reach in practice; neighbouring rods would coalesce into a film.`);
  if (mL / mD > 60) warnings.push(`Aspect ratio ${F(mL / mD, 3)}: slender rods of this kind bundle by capillary forces on drying and bend under their own adhesion; the straight-rod geometry is idealised.`);
  return {
    summary: `An array of ${mat.pretty} ${hex ? 'hexagonal' : 'circular'} nanorods at ${F(dens, 3)} rods/m² (${F(dens * 1e-12, 3)} per µm²), mean length ${F(mL, 4)} nm and diameter ${F(mD, 4)} nm, covers ${F(coverage * 100, 3)} % of the substrate and multiplies the surface area by ${F(sef, 4)}. Rods stand ${F(nnMean, 4)} nm apart centre to centre (gap ${F(gap, 3)} nm) with an alignment order parameter ⟨P₂⟩ of ${F(P2, 3)}. The forest is ${F(porosity * 100, 3)} % porous and carries ${F(massA * 100, 3)} mg/cm² of material.`,
    kpis: [
      kpi('dens', 'Number density', dens, '1/m²', `Rods generated per unit substrate area = ${F(dens * 1e-12, 4)} per µm²; the mean lattice pitch is ${F(pitch, 4)} nm.`),
      kpi('coverage', 'Areal coverage', coverage, '–', `Sum of rod base cross-sections (corrected for tilt by 1/cos θ) over the substrate area = ${F(coverage * 100, 3)} %. It is also the solid volume fraction at the base of the forest.`, coverage < 0.6),
      kpi('sef', 'Surface-area enhancement factor', sef, '–', `Total exposed area (side walls + tops + bare substrate) per unit substrate area = 1 + Σ(perimeter × length)/A = ${F(sef, 4)}. It multiplies sensing, catalytic or dye-loading capacity relative to a flat film.`, sef > 5),
      kpi('aspect', 'Mean aspect ratio L/d', mL / mD, '–', `Mean length ${F(mL, 4)} nm (std ${F(N.std(Ls), 3)} nm) over mean diameter ${F(mD, 4)} nm (std ${F(N.std(Ds), 3)} nm).`),
      kpi('spacing', 'Nearest-neighbour spacing', nnMean * 1e-9, 'm', `Mean centre-to-centre distance to the nearest rod = ${F(nnMean, 4)} nm, between the random (Poisson) value ${F(nnPoisson, 4)} nm and the hexagonal-lattice value ${F(nnHex, 4)} nm for this density.`),
      kpi('gap', 'Mean inter-rod gap', gap * 1e-9, 'm', `Nearest-neighbour spacing minus the mean diameter = ${F(gap, 4)} nm: the channel width available for electrolyte, gas or polymer infiltration.`, gap > 0),
      kpi('P2', 'Alignment order parameter ⟨P₂⟩', P2, '–', `⟨(3cos²θ − 1)/2⟩ over all rods (Hermans orientation factor) = ${F(hermans, 4)}; mean tilt ${F(N.mean(Ts), 3)}° from the substrate normal.`, P2 > 0.9),
      kpi('porosity', 'Forest porosity', porosity, '–', `One minus the height-averaged solid volume fraction over the forest height ${F(hmax, 4)} nm = ${F(porosity * 100, 3)} %.`),
      kpi('open', 'Line-of-sight open fraction (normal incidence)', open, '–', `Probability that a ray normal to the substrate reaches it unobstructed, exp(−n⟨A_proj⟩) with the projected area of tilted rods = ${F(open * 100, 3)} %. Oblique flux is shadowed beyond about ${F(alphaC, 3)}° from the normal (tan α = gap/L).`),
      kpi('mass', 'Mass loading', massA * 100, 'mg/cm²', `ρ·ΣV/A = ${F(massA, 3)} kg/m² of ${mat.pretty} in the rods (seed layer excluded).`),
      kpi('ssa', 'Specific surface area of the rods', ssaTop / 1e3, 'm²/g', `Side-wall plus top area divided by rod mass = ${F(ssaTop / 1e3, 4)} m²/g; for long rods this tends to ${hex ? '4/(ρd)' : '4/(ρd)'} and is independent of length.`),
    ],
    series: [
      { key: 'phi_z', title: 'Solid volume fraction vs height', xlabel: 'height above substrate (nm)', ylabel: 'solid volume fraction (–)', lines: [{ name: 'solid fraction', x: zs, y: phiZ }], explain: 'Fraction of each horizontal slice occupied by rods. It equals the areal coverage at the base and falls to zero over the range of rod heights; the fall-off width reflects the length spread.' },
      { key: 'tilt', title: 'Tilt-angle distribution', xlabel: 'tilt from the substrate normal (°)', ylabel: 'number of rods', lines: [{ name: 'rods', x: hT.x, y: hT.y }], explain: 'Histogram of the rod tilt angles. A narrow distribution near zero gives a high order parameter and a strong (002)-type texture in diffraction.' },
      { key: 'nn', title: 'Nearest-neighbour spacing distribution', xlabel: 'centre-to-centre distance (nm)', ylabel: 'number of rods', lines: [{ name: 'rods', x: hS.x, y: hS.y }], explain: 'Distribution of the distance from each rod to its nearest neighbour. Distances below one diameter indicate touching or coalesced rods.' },
    ],
    tables: [{ key: 'stats', title: 'Array statistics', columns: ['Quantity', 'Mean', 'Std', 'Min', 'Max'], rows: [['Length (nm)', Ls], ['Diameter (nm)', Ds], ['Tilt (°)', Ts], ['Nearest-neighbour spacing (nm)', nn.length ? nn : [side]]].map(([q, a]) => [q, +N.mean(a).toPrecision(5), +N.std(a).toPrecision(4), +N.minOf(a).toPrecision(5), +N.maxOf(a).toPrecision(5)]), explain: 'Statistics of the generated rods. Lengths and diameters are log-normally distributed with the requested coefficients of variation; tilt is half-normal.' }],
    equations: [
      matEq(mat),
      eqn('formula', 'Coverage and area enhancement', hex ? 'coverage = n⟨(√3/2)d²/cos θ⟩;  SEF = 1 + n⟨2√3·d·L⟩  (d across flats)' : 'coverage = n⟨πd²/(4 cos θ)⟩;  SEF = 1 + n⟨πdL⟩'),
      eqn('formula', 'Alignment (Hermans) order parameter', '⟨P₂⟩ = ⟨(3cos²θ − 1)/2⟩'),
      eqn('formula', 'Spacing references', 'random: ⟨r_nn⟩ = 1/(2√n);  hexagonal lattice: a = (2/(√3 n))^(1/2);  gap = ⟨r_nn⟩ − ⟨d⟩'),
      eqn('formula', 'Volume fraction, mass loading and shadowing', 'φ(z) = (1/A) Σᵢ aᵢ/cos θᵢ · [z ≤ Lᵢcos θᵢ];  m″ = ρ ΣVᵢ/A;  P_open = exp(−n⟨A_proj⟩)'),
    ],
    geometry: { kind: 'rods', rods, area: [side, side], units: 'nm' },
    checks, warnings,
    recommendations: [
      rec('high', gap < 0.5 * mD ? `The mean gap of ${F(gap, 3)} nm is less than half the rod diameter: reduce the seed density (thinner or sparser seed layer, patterned template) to keep rods separate and the side walls accessible.` : `The gap of ${F(gap, 3)} nm leaves the side walls accessible; raising the density towards a gap of about one diameter would increase the area enhancement (now ${F(sef, 3)}) without coalescence.`, 'Nearest-neighbour spacing against rod diameter.'),
      rec('medium', P2 < 0.9 ? `Improve vertical alignment (⟨P₂⟩ = ${F(P2, 3)}): use a textured, c-axis-oriented seed layer and avoid over-dense nucleation, which forces rods to tilt.` : `Alignment is good (⟨P₂⟩ = ${F(P2, 3)}); maintain seed-layer texture and annealing conditions to keep it.`, 'Orientation order parameter.'),
      rec('medium', `For conformal over-coating by line-of-sight deposition keep the incidence within ${F(alphaC, 3)}° of the normal or rotate the substrate; beyond that the lower rod walls are shadowed. Use ALD or solution routes for full side-wall coverage.`, 'Geometric shadowing angle tan α = gap/L.'),
      rec('low', `The array delivers ${F(sef, 3)} times the area of a flat film with only ${F(massA * 100, 3)} mg/cm² of ${mat.pretty}; low-temperature solution growth of such arrays uses far less material and energy than depositing an equivalent-area thick porous film.`, 'Surface-area enhancement per unit mass loading; material and energy efficiency.', true),
    ],
    publish: { 'geom.diameter': mD * 1e-9, 'geom.length': mL * 1e-9, 'geom.aspect': mL / mD, 'geom.number_density': dens, 'geom.ssa': ssaTop, 'geom.porosity': porosity, 'geom.coverage': coverage, 'geom.area_enhancement': sef, 'geom.spacing': nnMean * 1e-9, 'geom.order_parameter': P2, 'geom.mass_loading': massA, 'prop.density': mat.rho, 'material.formula': mat.f },
    exports: [
      { name: 'nanorod_array.csv', mime: 'text/csv', content: 'x_nm,y_nm,length_nm,diameter_nm,tilt_deg,azimuth_deg\n' + rods.map(r => [r.x, r.y, r.length, r.diameter, r.tilt, r.azimuth].map(q => +q.toPrecision(6)).join(',')).join('\n') + '\n', description: 'Position, size and orientation of every rod in the modelled area.' },
      { name: 'nanorod_array.obj', mime: 'model/obj', content: toOBJ(mesh), description: `Prism mesh of the first ${min(300, n)} rods in nm (Wavefront OBJ).` },
    ],
    numerics: { method: 'Seeded stochastic array generation (hard-core placement, log-normal sizes, half-normal tilt); direct summation of rod geometry; O(n²) periodic nearest-neighbour search', dof: n, steps: zs.length, iterations: 1, residual: abs(dens / v.nDens - 1), notes: `Modelled area ${F(side, 4)} nm square with ${n} rods. Rods are straight prisms; ${ssa >= 0 ? 'bundling, tapering and base coalescence are not modelled.' : ''}` },
  };
}

// ---------------------------------------------------------------- porous microstructure
function runPorous(v) {
  const warnings = [], mat = material(v, warnings), n = round(v.pN), vox = v.feat / round(v.featVox), fs = round(v.featVox);
  const gen = seed => implicit.porousMedium({ type: v.pType, size: [n, n, n], porosity: v.phi, featureSize: fs, seed, voxel_um: vox * 1e6 });
  const reps = [gen(v.seed), gen(v.seed + 101), gen(v.seed + 202)], g = reps[0], st = implicit.voxelStats(g);
  if (g.capped) warnings.push(g.capped + '.');
  const phi = g.porosity, Sv = max(1e-30, g.ssa * 1e6), solid = max(1e-9, 1 - phi);       // interfacial area per bulk volume, 1/m
  const ssa = Sv / (mat.rho * solid), perc = g.percolates, nPerc = ['x', 'y', 'z'].filter(a => perc[a]).length;
  const tauW = g.tortuosity, tauB = phi ** -0.5, tauG = ['x', 'y', 'z'].map(a => g.tortuosityGeometric[a]).filter(q => q > 0), tauGm = tauG.length ? N.mean(tauG) : null;
  const tau = tauW ?? tauB, kKC = phi ** 3 / (5 * Sv * Sv), dH = 4 * phi / Sv, dBrug = phi ** 1.5, dWalk = tauW ? phi / tauW : 0, chord = g.meanPoreChord * 1e-6;
  // chord-length (pore-size) distribution, all three axes pooled
  const pool = new Map(); for (const a of ['x', 'y', 'z']) for (const [len, c] of st.chords[a].histogram) pool.set(len, (pool.get(len) || 0) + c);
  const cl = [...pool.keys()].sort((a, b) => a - b), tot = cl.reduce((s, k) => s + pool.get(k) * k, 0) || 1;
  const cX = cl.map(k => k * 1e3), cY = cl.map(k => pool.get(k) * k / tot * 100);
  // representative-volume adequacy across three seeds
  const cvOf = a => { const m = N.mean(a); return m ? N.std(a) / abs(m) : 0; };
  const rP = reps.map(r => r.porosity), rS = reps.map(r => r.ssa * 1e6), rT = reps.map(r => r.tortuosity).filter(q => q > 0), cvP = cvOf(rP), cvS = cvOf(rS), cvT = rT.length > 1 ? cvOf(rT) : 0;
  const spanF = n / fs, adequate = spanF >= 4 && cvP < 0.03 && cvS < 0.05 && cvT < 0.2;
  const checks = [
    chk('Achieved porosity matches the target', abs(phi - v.phi) < 0.03, `achieved ${F(phi, 4)} against target ${F(v.phi, 4)}`),
    chk('Pore space percolates in all three directions', nPerc === 3, `percolation x/y/z: ${perc.x}/${perc.y}/${perc.z} (6-connected flood fill)`),
    chk('Solid skeleton is connected (self-supporting)', g.solidPercolates.x && g.solidPercolates.y && g.solidPercolates.z, `solid percolation x/y/z: ${g.solidPercolates.x}/${g.solidPercolates.y}/${g.solidPercolates.z}`),
    chk('Representative volume adequate (three seeds)', adequate, `RVE spans ${F(spanF, 3)} features; coefficient of variation: porosity ${F(cvP * 100, 3)} %, surface area ${F(cvS * 100, 3)} %, tortuosity ${F(cvT * 100, 3)} %`),
    chk('Tortuosity factor admissible (τ ≥ 1)', tauW === null || tauW >= 1 - 1e-9, tauW === null ? 'pore space does not percolate: no transport tortuosity' : `random-walk τ = ${F(tauW, 4)}`),
    chk('Resolution: at least five voxels per feature', fs >= 5, `${fs} voxels per feature (voxel ${F(nm(vox), 3)} nm)`),
  ];
  if (tauW === null) warnings.push('The pore space does not percolate: permeability and effective diffusivity of the connected network are zero in at least one direction; the Kozeny–Carman and Bruggeman values quoted assume an open network.');
  if (!adequate) warnings.push('The representative volume is small for this microstructure: properties vary noticeably between seeds. Increase the grid size or reduce the voxels per feature.');
  if (kKC < 1e-21) warnings.push('Pores of a few nanometres: gas flow is in the Knudsen/slip regime and liquid flow is affected by surface forces, so the continuum (Darcy) permeability is only a lower-bound indication.');
  return {
    summary: `A ${n}³-voxel ${DEF.pType.options.find(o => o.value === v.pType).label.toLowerCase()} microstructure of ${mat.pretty} with ${F(nm(v.feat), 3)} nm features has porosity ${F(phi, 3)} and ${F(ssa / 1e3, 4)} m²/g of internal surface. The pore space ${nPerc === 3 ? 'percolates in all directions' : `percolates in ${nPerc} of 3 directions`}, with a tortuosity factor of ${F(tau, 3)}${tauW ? '' : ' (Bruggeman estimate)'}, an effective diffusivity ratio D_eff/D₀ of ${F(tauW ? dWalk : dBrug, 3)} and a Kozeny–Carman permeability of ${F(kKC, 3)} m². Across three seeds the porosity varies by ${F(cvP * 100, 2)} % and the surface area by ${F(cvS * 100, 2)} %.`,
    kpis: [
      kpi('phi', 'Porosity', phi, '–', `Pore-voxel fraction of the generated volume = ${F(phi, 4)} (target ${F(v.phi, 3)}). It is the void volume available to fluid.`, abs(phi - v.phi) < 0.03),
      kpi('ssa', 'Specific surface area', ssa / 1e3, 'm²/g', `Interfacial area per bulk volume S_v = ${F(Sv, 3)} 1/m (voxel face count × 2/3) divided by ρ(1 − φ) = ${F(ssa / 1e3, 4)} m²/g. The 2/3 factor corrects the staircase overestimate for an isotropic surface.`),
      kpi('chord', 'Mean pore chord length', chord, 'm', `Mean length of line intercepts through the pore phase along x, y and z = ${F(nm(chord), 4)} nm, a direct measure of pore size; the hydraulic diameter 4φ/S_v is ${F(nm(dH), 4)} nm.`),
      kpi('tau', 'Tortuosity factor', tau, '–', tauW ? `τ = D₀/D_pore from a lattice random walk in the pore space = ${F(tauW, 4)} (statistical scatter about ±15 %). The Bruggeman estimate φ^(−1/2) is ${F(tauB, 4)}${tauGm ? ` and the geometric shortest-path tortuosity ${F(tauGm, 4)}` : ''}.` : `The pore space does not percolate, so the Bruggeman estimate φ^(−1/2) = ${F(tauB, 4)} is quoted for an open network of the same porosity.`, tau < 3),
      kpi('Deff', 'Effective diffusivity ratio D_eff/D₀', tauW ? dWalk : dBrug, '–', `φ/τ = ${F(tauW ? dWalk : dBrug, 4)}${tauW ? ` from the random-walk tortuosity; the Bruggeman relation φ^1.5 gives ${F(dBrug, 4)}` : ' from the Bruggeman relation φ^1.5'}. Multiply by the bulk (or Knudsen-corrected) diffusivity to obtain the effective value.`),
      kpi('k', 'Permeability (Kozeny–Carman)', kKC, 'm²', `k = φ³/(5 S_v²) = ${F(kKC, 3)} m² with the Kozeny constant 5. It is an order-of-magnitude estimate for a connected, roughly isotropic pore space.`),
      kpi('perc', 'Percolating directions', nPerc, 'of 3', `Number of axes along which the pore phase connects opposite faces (x: ${perc.x}, y: ${perc.y}, z: ${perc.z}). Transport across the sample needs percolation along the transport direction.`, nPerc === 3),
      kpi('rve', 'RVE porosity variation across seeds', cvP * 100, '%', `Coefficient of variation of porosity over three independently seeded volumes of the same size = ${F(cvP * 100, 3)} %; surface area varies by ${F(cvS * 100, 3)} %. Values below a few per cent indicate a representative volume.`, adequate),
    ],
    series: [
      { key: 'chords', title: 'Pore-size (chord-length) distribution', xlabel: 'chord length (nm)', ylabel: 'length-weighted fraction (%)', lines: [{ name: 'pore chords (x, y, z pooled)', x: cX.length ? cX : [0], y: cY.length ? cY : [0] }], explain: 'Distribution of the lengths of straight line segments lying wholly in the pore phase. The peak marks the typical pore width; the tail shows large connected channels. The resolution is one voxel.' },
      { key: 'transport', title: 'Transport estimates vs porosity (this microstructure marked)', xlabel: 'porosity (–)', ylabel: 'D_eff/D₀ (–)', lines: [{ name: 'Bruggeman φ^1.5', x: N.linspace(0.05, 0.95, 40), y: N.linspace(0.05, 0.95, 40).map(p => p ** 1.5) }, { name: 'this RVE (three seeds)', x: reps.map(r => r.porosity), y: reps.map(r => r.tortuosity ? r.porosity / r.tortuosity : 0) }], explain: 'Effective diffusivity relative to the free fluid. Points for the three seeded volumes lie below the Bruggeman curve when the pore space is more tortuous than a packed bed of spheres, and at zero when it does not percolate.' },
    ],
    tables: [{ key: 'rve', title: 'Representative-volume study (three seeds)', columns: ['Seed', 'Porosity', 'S_v (1/µm)', 'Tortuosity factor', 'Percolates x/y/z'], rows: reps.map(r => [r.seed, +r.porosity.toFixed(4), +r.ssa.toPrecision(4), r.tortuosity ? +r.tortuosity.toPrecision(4) : 'none', `${r.percolates.x ? 'y' : 'n'}/${r.percolates.y ? 'y' : 'n'}/${r.percolates.z ? 'y' : 'n'}`]), explain: 'The same generator run with three seeds on the same grid. Scatter between rows is the statistical uncertainty of a single volume of this size; it shrinks as the volume grows relative to the feature size.' }],
    equations: [
      matEq(mat),
      eqn('formula', 'Porosity and specific surface', 'φ = V_pore/V;  S_v = (2/3)·N_faces·h²/V;  SSA = S_v/(ρ(1 − φ));  d_h = 4φ/S_v'),
      eqn('constitutive', 'Kozeny–Carman permeability', 'k = φ³ / (c_K S_v²),  c_K = 5'),
      eqn('constitutive', 'Effective diffusivity', 'D_eff = D₀ φ/τ;  Bruggeman: τ = φ^(−1/2), D_eff/D₀ = φ^1.5'),
      eqn('governing', 'Tortuosity from a pore-space random walk', 'τ = D₀/D_pore = lim t/⟨r²(t)⟩  (lattice walkers, periodic wrapping)'),
      eqn('governing', 'Darcy flow through the medium', 'u = −(k/μ)∇p'),
    ],
    geometry: { kind: 'voxels', nx: g.nx, ny: g.ny, nz: g.nz, data: g.data, voxel: nm(vox), units: 'nm' },
    checks, warnings,
    recommendations: [
      rec('high', nPerc < 3 ? `The pores are not fully connected at φ = ${F(phi, 3)}: raise the porosity (pore former, lower sintering temperature) above the percolation threshold of this microstructure, or closed pores will not contribute to transport.` : `The pore network is connected; D_eff/D₀ = ${F(tauW ? dWalk : dBrug, 3)} means diffusion through the layer is ${F(1 / max(1e-9, tauW ? dWalk : dBrug), 3)} times slower than in the free fluid — size the layer thickness accordingly.`, 'Percolation and tortuosity of the generated volume.'),
      rec('medium', `Permeability scales with the square of the feature size: doubling the ${F(nm(v.feat), 3)} nm feature raises k from ${F(kKC, 2)} to about ${F(4 * kKC, 2)} m² but halves the specific surface area (${F(ssa / 1e3, 3)} m²/g). Choose the feature size from the transport-versus-surface trade-off.`, 'Kozeny–Carman scaling k ∝ φ³/S_v².'),
      rec(adequate ? 'low' : 'high', adequate ? 'The volume is representative for porosity and surface area; keep at least this size relative to the feature when importing tomography data.' : `Increase the RVE: it spans only ${F(spanF, 3)} features and properties vary between seeds. Average several realisations before passing values to transport suites.`, 'Variation across three seeded realisations.'),
      rec('low', `A hierarchical pore structure (large transport pores feeding ${F(nm(chord), 3)} nm pores) reaches the same accessible surface with lower pressure drop, cutting pumping energy and allowing thinner, lighter layers of ${mat.pretty}.`, 'Trade-off between Kozeny–Carman permeability and specific surface; energy and material efficiency.', true),
    ],
    publish: { 'geom.porosity': phi, 'geom.ssa': ssa, 'geom.tortuosity': tau, 'geom.pore_size': chord, 'eff.permeability': kKC, 'eff.D_ratio': tauW ? dWalk : dBrug, 'eff.density': mat.rho * solid, 'geom.diameter': v.feat, 'prop.density': mat.rho, 'material.formula': mat.f },
    exports: [
      { name: `porous_${n}x${n}x${n}_uint8.raw`, mime: 'application/octet-stream', bytes: new Uint8Array(g.data), description: `Voxel phase map (1 = solid, 0 = pore), x fastest, voxel ${F(nm(vox), 4)} nm.` },
      { name: 'porous_metadata.json', mime: 'application/json', content: JSON.stringify({ dims: [n, n, n], voxel_m: vox, model: g.model, seed: v.seed, porosity: phi, Sv_per_m: Sv, ssa_m2_per_g: ssa / 1e3, tortuosity_factor: tauW, permeability_m2: kKC, percolates: perc }, null, 1), description: 'Grid description and computed descriptors for the raw voxel file.' },
      { name: 'pore_chords.csv', mime: 'text/csv', content: 'chord_nm,length_weighted_fraction\n' + cX.map((x, i) => `${x.toPrecision(5)},${(cY[i] / 100).toPrecision(5)}`).join('\n') + '\n', description: 'Pore chord-length distribution.' },
    ],
    numerics: { method: `Stochastic voxel generation (${g.model}); face-count surface area; flood-fill percolation; lattice random walk for tortuosity`, dof: n ** 3, steps: 3, iterations: 3, residual: abs(phi - v.phi), notes: `Voxel ${F(nm(vox), 4)} nm, ${fs} voxels per feature. Three seeds for the representative-volume indicator.` },
  };
}

// ---------------------------------------------------------------- aggregates and packing
const PACK_REF = [['Random sequential adsorption of discs on a plane (jamming)', 0.547, '2-D'], ['Random close packing of discs', 0.84, '2-D'], ['Hexagonal close packing of discs', 0.9069, '2-D'], ['Random sequential adsorption of spheres in a volume', 0.38, '3-D'], ['Random loose packing of spheres', 0.55, '3-D'], ['Random close packing of spheres', 0.64, '3-D'], ['Face-centred cubic / hexagonal close packing', 0.7405, '3-D']];
function runDLA(v, mat, warnings) {
  const rand = N.rng(v.seed), G = 161, c = 80, occ = new Uint8Array(G * G * G), id = (x, y, z) => x + G * (y + G * z), target = round(v.nPart);
  const P = [[c, c, c]]; occ[id(c, c, c)] = 1;
  let rmax = 0, steps = 0, sx = c, sy = c, sz = c, s2 = 3 * c * c; const budget = 4e6, rgN = [], rgR = []; let nextRec = 8;
  const NB = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
  outer: while (P.length < target && rmax < 34) {
    const Rl = rmax + 4, Rk2 = min(78, 2 * Rl + 6) ** 2, near2 = (rmax + 2) ** 2; let x, y, z;
    const launch = () => { const u = 2 * rand() - 1, a = 2 * PI * rand(), s = sqrt(1 - u * u); x = round(c + Rl * s * cos(a)); y = round(c + Rl * s * sin(a)); z = round(c + Rl * u); };
    launch();
    for (;;) {
      if (++steps > budget) break outer;
      const m = NB[(rand() * 6) | 0], far = sqrt((x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2) - rmax - 2, st = far > 2 ? floor(far) : 1;   // long jumps while far from the cluster
      const X = x + st * m[0], Y = y + st * m[1], Z = z + st * m[2], r2 = (X - c) ** 2 + (Y - c) ** 2 + (Z - c) ** 2;
      if (r2 > Rk2) { launch(); continue; }
      if (occ[id(X, Y, Z)]) continue;
      x = X; y = Y; z = Z;
      if (r2 <= near2 && (occ[id(x + 1, y, z)] || occ[id(x - 1, y, z)] || occ[id(x, y + 1, z)] || occ[id(x, y - 1, z)] || occ[id(x, y, z + 1)] || occ[id(x, y, z - 1)]) && rand() < v.stick) {
        occ[id(x, y, z)] = 1; P.push([x, y, z]); rmax = max(rmax, sqrt(r2)); sx += x; sy += y; sz += z; s2 += x * x + y * y + z * z;
        if (P.length >= nextRec || P.length === target) { const n = P.length; rgN.push(n); rgR.push(sqrt(max(1e-12, s2 / n - (sx / n) ** 2 - (sy / n) ** 2 - (sz / n) ** 2))); nextRec = Math.ceil(nextRec * 1.25); }
        break;
      }
    }
  }
  const n = P.length; if (!rgN.length || rgN[rgN.length - 1] !== n) { rgN.push(n); rgR.push(sqrt(max(1e-12, s2 / n - (sx / n) ** 2 - (sy / n) ** 2 - (sz / n) ** 2))); }
  if (n < target) warnings.push(`The cluster stopped at ${n} of ${target} particles (${steps > budget ? 'walker step budget reached' : 'lattice radius limit reached'}).`);
  if (n < 150) warnings.push('Fewer than about 150 particles: the cluster is too small for a reliable fractal dimension; the fitted exponent is indicative only.');
  const Rg = rgR[rgR.length - 1], i0 = rgN.findIndex(q => q >= 12), fitN = rgN.slice(max(0, i0)), fitR = rgR.slice(max(0, i0));
  const lr = fitN.length >= 3 ? N.linreg(fitR.map(log), fitN.map(log)) : { slope: NaN, intercept: 0, r2: 0 }, Df = fin(lr.slope, 0), k0 = exp(fin(lr.intercept));
  if (Df > 3) warnings.push(`The fitted exponent ${F(Df, 3)} exceeds 3: at low sticking probability walkers fill the interior of a still-compact cluster, so it has not reached the fractal scaling regime. Grow more particles before quoting a fractal dimension.`);
  // box counting on the final cluster
  const bx = [], bc = []; for (let s = 1; s <= 32; s *= 2) { const set = new Set(); for (const p of P) set.add(floor(p[0] / s) + 512 * (floor(p[1] / s) + 512 * floor(p[2] / s))); bx.push(s); bc.push(set.size); }
  const use = bx.map((s, i) => i).filter(i => bx[i] <= max(2, rmax / 1.5) && bc[i] > 1), lb = use.length >= 2 ? N.linreg(use.map(i => log(bx[i])), use.map(i => log(bc[i]))) : { slope: NaN, r2: 0 }, Dbox = fin(-lb.slope, 0);
  let bonds = 0; for (const p of P) for (const m of NB) if (occ[id(p[0] + m[0], p[1] + m[1], p[2] + m[2])]) bonds++;
  const coord = n > 1 ? bonds / n : 0, Rs = sqrt(5 / 3) * Rg, phiIn = Rs > 0 ? min(1, n * 0.125 / Rs ** 3) : 1;
  // cropped voxel geometry
  const lo = [0, 1, 2].map(a => N.minOf(P.map(p => p[a]))), hi = [0, 1, 2].map(a => N.maxOf(P.map(p => p[a]))), dims = [0, 1, 2].map(a => hi[a] - lo[a] + 1), data = new Uint8Array(dims[0] * dims[1] * dims[2]);
  for (const p of P) data[(p[0] - lo[0]) + dims[0] * ((p[1] - lo[1]) + dims[1] * (p[2] - lo[2]))] = 1;
  const d = v.d, ssa = 6 / (mat.rho * d), fitLine = fitR.map(r => k0 * r ** Df);
  return {
    summary: `Diffusion-limited aggregation of ${n} primary particles of ${mat.pretty} (${F(nm(d), 3)} nm, sticking probability ${F(v.stick, 3)}) gives a cluster with radius of gyration ${F(nm(Rg * d), 4)} nm and a mass-fractal dimension of ${F(Df, 3)} from radius-of-gyration scaling (${F(Dbox, 3)} by box counting). Each particle touches ${F(coord, 3)} neighbours on average and the solid fills ${F(phiIn * 100, 3)} % of the equivalent sphere.`,
    kpis: [
      kpi('Df', 'Fractal dimension (R_g scaling)', Df, '–', `Slope of ln N against ln R_g during growth = ${F(Df, 3)} (R² = ${F(lr.r2, 4)}). Particle–cluster DLA in three dimensions gives about 2.5; cluster–cluster (soot-like) aggregates about 1.8; compact packing 3.`, Df > 1 && Df <= 3.05),
      kpi('Dbox', 'Fractal dimension (box counting)', Dbox, '–', `Slope of ln N_box against ln(1/box size) on the final cluster = ${F(Dbox, 3)}. For clusters of this size box counting is biased low by the finite scaling range.`),
      kpi('Rg', 'Radius of gyration', Rg * d, 'm', `Root-mean-square distance of the primary particles from the centre of mass = ${F(Rg, 4)} particle diameters = ${F(nm(Rg * d), 4)} nm. It is the size measured by static light or X-ray scattering (Guinier regime).`),
      kpi('coord', 'Mean coordination number', coord, '–', `Average number of touching nearest neighbours per primary particle on the cubic lattice = ${F(coord, 3)} (maximum 6). Values near 2 indicate chain-like, open structures.`),
      kpi('k0', 'Fractal prefactor k₀', k0, '–', `Prefactor in N = k₀ (R_g/d)^D_f fitted with R_g in units of the particle diameter = ${F(k0, 3)}.`),
      kpi('phiIn', 'Packing fraction inside the equivalent sphere', phiIn, '–', `Solid volume of ${n} spheres divided by the sphere of radius (5/3)^(1/2) R_g = ${F(phiIn * 100, 3)} %. It falls with size for a fractal (D_f < 3), which is why aggregates are so voluminous.`),
      kpi('n', 'Primary particles in the cluster', n, '–', `Number of particles aggregated out of ${target} requested, using ${F(steps, 3)} random-walk steps.`, n === target),
      kpi('ssa', 'Specific surface area (point contacts)', ssa / 1e3, 'm²/g', `6/(ρd) = ${F(ssa / 1e3, 4)} m²/g: an open aggregate with point contacts keeps the surface of its primary particles; sintered necks reduce it.`),
    ],
    series: [
      { key: 'rg', title: 'Radius-of-gyration scaling N(R_g)', xlabel: 'R_g / d (–)', ylabel: 'number of particles N', logx: true, logy: true, lines: [{ name: 'growing cluster', x: rgR, y: rgN }, { name: `fit N = ${F(k0, 3)}·(R_g/d)^${F(Df, 3)}`, x: fitR, y: fitLine }], explain: 'A straight line on logarithmic axes shows mass-fractal scaling; its slope is the fractal dimension. Early points (few particles) are excluded from the fit.' },
      { key: 'box', title: 'Box counting of the final cluster', xlabel: 'box size (particle diameters)', ylabel: 'occupied boxes', logx: true, logy: true, lines: [{ name: 'occupied boxes', x: bx, y: bc }], explain: 'Number of boxes of a given size needed to cover the cluster. The negative slope over the scaling range estimates the fractal dimension; it flattens once the box approaches the cluster size.' },
    ],
    tables: [{ key: 'pack', title: 'Reference packing fractions', columns: ['Structure', 'Packing / coverage fraction', 'Dimension'], rows: PACK_REF, explain: 'Standard literature values for comparison with the open fractal aggregate: random and ordered packings of equal discs and spheres.' }],
    equations: [
      matEq(mat),
      eqn('formula', 'Mass-fractal scaling law', 'N = k₀ (R_g/d)^D_f'),
      eqn('formula', 'Radius of gyration', 'R_g² = (1/N) Σᵢ |rᵢ − r_cm|²'),
      eqn('formula', 'Box-counting dimension', 'D_box = −d ln N_box / d ln s'),
      eqn('governing', 'Diffusion-limited aggregation', 'Random walkers (∂c/∂t = D∇²c in the continuum limit) stick irreversibly on first contact with probability p_s'),
    ],
    geometry: { kind: 'voxels', nx: dims[0], ny: dims[1], nz: dims[2], data, voxel: nm(d), units: 'nm' },
    checks: [
      chk('Fractal dimension within physical bounds (1 ≤ D_f ≤ 3)', Df >= 1 && Df <= 3.05, `D_f = ${F(Df, 3)} from ${fitN.length} points, R² = ${F(lr.r2, 4)}`),
      chk('Power-law scaling holds over the growth history', lr.r2 > 0.98, `R² = ${F(lr.r2, 4)} of the log–log regression`),
      chk('Coordination number admissible on the cubic lattice', coord >= (n > 1 ? 1 : 0) && coord <= 6, `mean coordination ${F(coord, 3)} (every particle has at least one contact)`),
      chk('Cluster completed within the lattice and step budget', n === target, `${n} of ${target} particles; maximum radius ${F(rmax, 3)} lattice units`),
    ],
    extra: { ssa, Rg: Rg * d, Df },
    recommendations: [
      rec('high', Df < 2.2 ? `The open structure (D_f = ${F(Df, 3)}) fills space at very low solid fraction: expect high viscosity and gelation at low loading; disperse with shear or ultrasound and stabilise electrostatically or sterically.` : `The compact structure (D_f = ${F(Df, 3)}) packs efficiently; if an open, high-porosity network is wanted (aerogel, filter cake), move towards fast, diffusion-limited aggregation by screening the surface charge.`, 'Fractal dimension of the aggregate.'),
      rec('medium', `Use the fitted law N = ${F(k0, 3)}·(R_g/d)^${F(Df, 3)} to convert a scattering-derived R_g into a particle count and effective density for sedimentation and filtration calculations.`, 'Mass-fractal scaling fitted to the simulated growth.'),
      rec('medium', `Lowering the sticking probability (reaction-limited conditions: higher surface charge, lower ionic strength) densifies clusters; here p_s = ${F(v.stick, 3)} gave D_f = ${F(Df, 3)}.`, 'Sticking-probability dependence of lattice DLA.'),
      rec('low', 'Fractal aggregates of nanoparticles settle and filter far more readily than the dispersed primaries; controlled flocculation before filtration reduces membrane fouling, wash-water use and nanoparticle release to effluent.', 'Aggregate size against primary size; waste-water and exposure control.', true),
    ],
    publish: { 'geom.diameter': d, 'geom.ssa': ssa, 'geom.fractal_dimension': Df, 'geom.radius_gyration': Rg * d, 'geom.coordination': coord, 'geom.porosity': 1 - phiIn },
    exports: [{ name: 'dla_aggregate.xyz', mime: 'chemical/x-xyz', content: `${n}\nDLA aggregate; coordinates in nm (primary diameter ${F(nm(d), 4)} nm)\n` + P.map(p => `X ${((p[0] - c) * nm(d)).toPrecision(6)} ${((p[1] - c) * nm(d)).toPrecision(6)} ${((p[2] - c) * nm(d)).toPrecision(6)}`).join('\n') + '\n', description: 'Centres of the primary particles of the aggregate.' }],
    numerics: { method: 'On-lattice particle–cluster diffusion-limited aggregation (simple cubic, 6-neighbour contact, launch and kill spheres)', dof: n, steps, iterations: n, residual: 1 - lr.r2, notes: 'Walkers take long lattice jumps while far from the cluster and are relaunched beyond the kill radius, which slightly biases the angular flux. Lattice anisotropy becomes visible only for much larger clusters.' },
  };
}
function runRSA(v, mat, warnings) {
  const rand = N.rng(v.seed), Lr = round(v.rsaL), cells = Array.from({ length: Lr * Lr }, () => []), P = [], d = v.d;
  const attempts = round(v.rsaTau * Lr * Lr / (PI / 4)), tx = [], ty = []; let next = 0.25;
  for (let a = 1; a <= attempts; a++) {
    const x = rand() * Lr, y = rand() * Lr, ci = floor(x), cj = floor(y); let ok = true;
    for (let di = -1; di <= 1 && ok; di++) for (let dj = -1; dj <= 1 && ok; dj++) for (const q of cells[(ci + di + Lr) % Lr + Lr * ((cj + dj + Lr) % Lr)]) { let dx = abs(q[0] - x), dy = abs(q[1] - y); if (dx > Lr / 2) dx = Lr - dx; if (dy > Lr / 2) dy = Lr - dy; if (dx * dx + dy * dy < 1) { ok = false; break; } }
    if (ok) { P.push([x, y]); cells[ci + Lr * cj].push([x, y]); }
    const tau = a * (PI / 4) / (Lr * Lr);
    if (tau >= next || a === attempts) { tx.push(tau); ty.push(P.length * (PI / 4) / (Lr * Lr)); next = tau * 1.15; }
  }
  const n = P.length, theta = n * (PI / 4) / (Lr * Lr), late = tx.map((t, i) => i).filter(i => tx[i] >= min(15, 0.4 * v.rsaTau));
  const fl = late.length >= 3 ? N.linreg(late.map(i => 1 / sqrt(tx[i])), late.map(i => ty[i])) : { intercept: theta, slope: 0, r2: 0 }, thInf = min(0.6, max(theta, fin(fl.intercept, theta)));
  const nn = P.map((a, i) => { let m = Infinity; for (let j = 0; j < n; j++) { if (j === i) continue; let dx = abs(a[0] - P[j][0]), dy = abs(a[1] - P[j][1]); if (dx > Lr / 2) dx = Lr - dx; if (dy > Lr / 2) dy = Lr - dy; m = min(m, dx * dx + dy * dy); } return sqrt(m); }).filter(Number.isFinite);
  const nnMean = nn.length ? N.mean(nn) : Lr, dens = n / (Lr * d) ** 2, massA = dens * mat.rho * PI * d ** 3 / 6, dn = nm(d), seg = 16;
  const polygons = P.map(p => ({ layer: 1, pts: Array.from({ length: seg }, (_, k) => [(p[0] + 0.5 * cos(2 * PI * k / seg)) * dn, (p[1] + 0.5 * sin(2 * PI * k / seg)) * dn]) }));
  const tol = Lr >= 25 ? 0.02 : 0.04;
  if (Lr < 20) warnings.push('Small RSA domain: the coverage of a single realisation scatters by a few per cent; increase the domain side for a tighter estimate.');
  return {
    summary: `Random sequential adsorption of ${F(dn, 3)} nm ${mat.pretty} particles on a ${Lr} × ${Lr} diameter periodic surface deposited ${n} particles after ${F(v.rsaTau, 3)} attempts per disc area, reaching a coverage of ${F(theta, 4)}. Extrapolating with Feder's law gives a jamming coverage of ${F(thInf, 4)}, to be compared with the accepted value 0.547. The mean centre-to-centre spacing is ${F(nnMean * dn, 4)} nm.`,
    kpis: [
      kpi('theta', 'Coverage reached', theta, '–', `Fraction of the surface covered by the projected discs after ${attempts} deposition attempts = ${F(theta, 4)}. Irreversible, non-overlapping adsorption cannot exceed the jamming limit.`, theta < 0.56),
      kpi('thInf', 'Extrapolated jamming coverage', thInf, '–', `Intercept of coverage against τ^(−1/2) (Feder's law) = ${F(thInf, 4)}; the accepted jamming limit for discs is 0.547.`, abs(thInf - 0.547) < tol),
      kpi('dens', 'Surface number density', dens, '1/m²', `${n} particles on (${Lr}·d)² = ${F(dens, 3)} 1/m² (${F(dens * 1e-12, 4)} per µm²) for d = ${F(dn, 3)} nm.`),
      kpi('nn', 'Mean nearest-neighbour spacing', nnMean * d, 'm', `Mean centre-to-centre distance to the closest particle = ${F(nnMean, 4)} diameters = ${F(nnMean * dn, 4)} nm; a hexagonal monolayer of the same particles would have 1.000 diameter.`),
      kpi('gapfrac', 'Uncovered area fraction', 1 - theta, '–', `Share of the substrate left bare = ${F((1 - theta) * 100, 3)} %; at jamming about 45 % of the surface is still exposed but no gap is large enough for another particle.`),
      kpi('mass', 'Areal mass of the monolayer', massA * 1e5, 'µg/cm²', `Number density × particle mass (ρπd³/6) = ${F(massA, 3)} kg/m².`),
      kpi('rel', 'Coverage relative to hexagonal close packing', theta / 0.9069, '–', `Coverage divided by π/(2√3) = 0.9069: random irreversible adsorption reaches only ${F(theta / 0.9069 * 100, 3)} % of an ordered monolayer.`),
    ],
    series: [
      { key: 'kin', title: 'Coverage vs deposition attempts', xlabel: 'attempts per disc area τ (–)', ylabel: 'coverage θ (–)', logx: true, lines: [{ name: 'simulation', x: tx, y: ty }, { name: 'jamming limit 0.547', x: [tx[0], tx[tx.length - 1]], y: [0.547, 0.547] }], explain: 'Adsorption is fast at first (Langmuir-like) and then slows sharply as the remaining gaps become rare; the approach to jamming follows a power law rather than an exponential.' },
      { key: 'feder', title: "Feder's law extrapolation", xlabel: 'τ^(−1/2) (–)', ylabel: 'coverage θ (–)', lines: [{ name: 'simulation', x: tx.map(t => 1 / sqrt(t)), y: ty }, { name: 'late-time fit', x: [0, late.length ? 1 / sqrt(tx[late[0]]) : 1], y: [fin(fl.intercept, theta), fin(fl.intercept, theta) + fin(fl.slope) * (late.length ? 1 / sqrt(tx[late[0]]) : 1)] }], explain: 'At late times θ(τ) = θ_J − k τ^(−1/2); the intercept of the straight line at τ^(−1/2) → 0 is the jamming coverage.' },
    ],
    tables: [{ key: 'pack', title: 'Reference packing fractions', columns: ['Structure', 'Packing / coverage fraction', 'Dimension'], rows: PACK_REF, explain: 'Standard literature values for random and ordered packings of equal discs and spheres, for comparison with the simulated coverage.' }],
    equations: [
      matEq(mat),
      eqn('governing', 'Random sequential adsorption', 'Each trial position is uniform on the surface; a particle is added irreversibly if it overlaps no earlier particle (|rᵢ − rⱼ| ≥ d)'),
      eqn('formula', 'Coverage and number density', 'θ = N πd²/(4A);  n = N/A'),
      eqn('formula', "Feder's law (approach to jamming)", 'θ(τ) = θ_J − k τ^(−1/2),  θ_J ≈ 0.547 for discs'),
    ],
    geometry: { kind: 'layout', polygons, units: 'nm' },
    checks: [
      chk('Coverage does not exceed the jamming limit', theta <= 0.56, `θ = ${F(theta, 4)} against θ_J ≈ 0.547`),
      chk('Extrapolated jamming coverage agrees with 0.547', abs(thInf - 0.547) < tol, `θ_J = ${F(thInf, 4)} (tolerance ±${tol} for a ${Lr} × ${Lr} domain, single realisation)`),
      chk('No overlapping particles', !nn.length || N.minOf(nn) >= 1 - 1e-9, `smallest centre spacing ${F(nn.length ? N.minOf(nn) : 1, 5)} diameters`),
    ],
    extra: { ssa: 6 / (mat.rho * d), dens, theta },
    recommendations: [
      rec('high', `Do not design for a monolayer coverage above about 0.55 from irreversible adsorption: this run saturates at ${F(thInf, 3)}. Higher coverage needs lateral mobility (annealing, capillary assembly, Langmuir–Blodgett transfer).`, 'RSA jamming limit.'),
      rec('medium', `Stop the deposition near τ ≈ ${F(min(v.rsaTau, 20), 3)} attempts per disc area: beyond that the coverage gain is marginal (θ grows as −τ^(−1/2)) while the exposure time and precursor use keep rising.`, "Feder's law kinetics."),
      rec('medium', `Electrostatic repulsion enlarges the effective exclusion diameter: with a Debye length λ the coverage falls roughly as θ_J (d/(d + 2λ))²; raise the ionic strength to approach ${F(thInf, 3)}.`, 'Effective hard-disc scaling of RSA coverage.'),
      rec('low', `A monolayer needs only ${F(massA * 1e5, 3)} µg/cm² of ${mat.pretty}; reuse the depleted suspension for subsequent substrates rather than discarding it, since less than a monolayer is consumed per run.`, 'Areal mass of the adsorbed layer; precursor efficiency.', true),
    ],
    publish: { 'geom.diameter': d, 'geom.number_density': dens, 'geom.coverage': theta, 'geom.spacing': nnMean * d, 'geom.ssa': 6 / (mat.rho * d) },
    exports: [{ name: 'rsa_positions.csv', mime: 'text/csv', content: 'x_nm,y_nm,diameter_nm\n' + P.map(p => `${(p[0] * dn).toPrecision(6)},${(p[1] * dn).toPrecision(6)},${dn.toPrecision(5)}`).join('\n') + '\n', description: 'Centres of the adsorbed particles on the periodic substrate.' }],
    numerics: { method: 'Random sequential adsorption with a cell list on a periodic square domain; linear regression of θ against τ^(−1/2)', dof: n, steps: attempts, iterations: attempts, residual: abs(thInf - 0.547), notes: 'Single seeded realisation; hard discs, no surface diffusion or desorption.' },
  };
}
function runAggregate(v) {
  const warnings = [], mat = material(v, warnings), r = v.aggType === 'dla' ? runDLA(v, mat, warnings) : runRSA(v, mat, warnings);
  delete r.extra;
  return { ...r, warnings, publish: { ...r.publish, 'prop.density': mat.rho, 'material.formula': mat.f } };
}

// ---------------------------------------------------------------- thin film / multilayer
function selfAffine(n, scan, xi, H, rms, rand) {   // spectral synthesis, PSD ∝ (1 + (qξ)²)^−(1+H)
  const re = Array.from({ length: n }, () => new Array(n).fill(0)), im = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    const kx = i <= n / 2 ? i : i - n, ky = j <= n / 2 ? j : j - n, q = 2 * PI * sqrt(kx * kx + ky * ky) / scan;
    if (!kx && !ky) continue;
    const a = (1 + (q * xi) ** 2) ** (-(1 + H) / 2); re[j][i] = a * rand.normal(); im[j][i] = a * rand.normal();
  }
  for (let j = 0; j < n; j++) N.fft(re[j], im[j], true);
  for (let i = 0; i < n; i++) { const cr = re.map(r => r[i]), ci = im.map(r => r[i]); N.fft(cr, ci, true); for (let j = 0; j < n; j++) re[j][i] = cr[j]; }
  const h = new Float64Array(n * n); let m = 0, s = 0;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { h[j * n + i] = re[j][i]; m += re[j][i]; }
  m /= n * n; for (let k = 0; k < h.length; k++) { h[k] -= m; s += h[k] * h[k]; }
  const f = rms / (sqrt(s / h.length) || 1); for (let k = 0; k < h.length; k++) h[k] *= f;
  return h;
}
function runFilm(v) {
  const warnings = [], rand = N.rng(v.seed);
  let rows = parseRows(v.stack); const layers = [];
  const read = rs => { for (const r of rs) { if (r.length < 3) continue; const t = Number(r[2]); if (!(t > 0)) continue; let info; try { const rho = Number(r[3]); info = matInfo(r[1], rho > 0 ? min(25000, max(50, rho)) : 3000, r[0], null); if (!(rho > 0)) warnings.push(`Layer "${r[0]}": no valid density given; 3000 kg/m³ assumed.`); } catch { warnings.push(`Layer "${r[0]}": formula "${r[1]}" could not be parsed and the layer was skipped.`); continue; } layers.push({ name: r[0] || info.pretty, t: min(1e-3, t * 1e-9), ...info }); if (layers.length >= 12) break; } };
  read(rows);
  if (String(v.stack).trim() === DEMO_STACK) warnings.push('The layer stack is an example (synthetic demonstration data — replace with measurements of your own stack).');
  if (!layers.length) { warnings.push('No valid layer rows were found (need: layer, formula, thickness_nm, density_kg_m3); the example stack was used.'); read(parseRows(DEMO_STACK)); }
  const total = N.sum(layers.map(l => l.t)), massA = N.sum(layers.map(l => l.rho * l.t)), top = layers[layers.length - 1], thin = N.minOf(layers.map(l => l.t));
  // island-growth (Volmer–Weber, hemispherical caps with Avrami overlap) coverage of the top layer
  const Ni = v.islandN, cov = h => { const r = cbrt(3 * h / (2 * PI * Ni)); return 1 - exp(-Ni * PI * r * r); }, r99 = sqrt(log(100) / (PI * Ni)), h99 = 2 * PI * Ni * r99 ** 3 / 3;
  const hs = N.logspace(1e-11, max(3 * h99, top.t, 1e-9), 80), covTop = cov(top.t), covLbL = min(1, top.t / top.delta);
  // self-affine surface
  const n = 64, xi = min(v.xi, v.scan / 2), h = selfAffine(n, v.scan, xi, v.hurst, v.rmsT, rand), dx = v.scan / n;
  let ra = 0, s3 = 0, s4 = 0, hmin = Infinity, hmax = -Infinity, slope2 = 0, areaR = 0;
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const q = h[j * n + i]; ra += abs(q); s3 += q ** 3; s4 += q ** 4; hmin = min(hmin, q); hmax = max(hmax, q); const gx = (h[j * n + (i + 1) % n] - q) / dx, gy = (h[((j + 1) % n) * n + i] - q) / dx; slope2 += gx * gx + gy * gy; areaR += sqrt(1 + gx * gx + gy * gy); }
  const M = n * n, rq = v.rmsT; ra /= M; const skew = s3 / M / rq ** 3, kurt = s4 / M / rq ** 4, rmsSlope = sqrt(slope2 / M), areaRatio = areaR / M;
  // height–height correlation g(r) = ⟨(h(x+r) − h(x))²⟩ along x and y → Hurst estimate from the short-range slope
  const lags = [], gr = []; for (let r = 1; r <= n / 2; r = r < 8 ? r + 1 : Math.ceil(r * 1.3)) { let s = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const q = h[j * n + i]; s += (h[j * n + (i + r) % n] - q) ** 2 + (h[((j + r) % n) * n + i] - q) ** 2; } lags.push(r * dx); gr.push(s / (2 * M)); }
  const short = lags.map((l, i) => i).filter(i => lags[i] <= max(xi / 2, 3.5 * dx)).slice(0, 8), hf = short.length >= 2 ? N.linreg(short.map(i => log(lags[i])), short.map(i => log(max(1e-300, gr[i])))) : { slope: 2 * v.hurst }, Hest = min(1, max(0, fin(hf.slope / 2, v.hurst)));
  const xiEst = (() => { for (let i = 0; i < lags.length; i++) if (gr[i] >= 2 * rq * rq * (1 - 1 / Math.E)) return lags[i]; return lags[lags.length - 1]; })();
  // conformality on a high-aspect-ratio feature
  const AR = v.featAR, thiele = AR * sqrt(3 * v.sticking), scCVD = 1 / Math.cosh(min(700, thiele)), scPVD = 1 / sqrt(1 + 4 * AR * AR), aldMult = 1 + 19 * AR / 4 + 1.5 * AR * AR;
  const ars = N.logspace(0.5, 300, 60);
  // voxel picture of the stack with the rough top surface (vertical scale exaggerated when the roughness is small)
  const nv = 32, nzB = 20, nz = 28, data = new Uint8Array(nv * nv * nz), bounds = []; let acc = 0; for (const l of layers) { acc += l.t; bounds.push(acc / total * nzB); }
  const amp = max(rq / total * nzB, 1.2) / (rq || 1) * 1.0;
  for (let j = 0; j < nv; j++) for (let i = 0; i < nv; i++) { const topZ = nzB + h[(j * 2) * n + i * 2] * amp; for (let k = 0; k < nz; k++) { if (k + 0.5 > topZ) break; let p = layers.length; for (let b = 0; b < bounds.length; b++) if (k + 0.5 <= bounds[b]) { p = b + 1; break; } data[i + nv * (j + nv * k)] = p; } }
  const stackEq = layers.map(l => eqn('formula', `Layer: ${l.name}`, `${l.pretty}: M = ${F(l.M * 1e3, 5)} g/mol, ρ = ${F(l.rho, 4)} kg/m³, t = ${F(nm(l.t), 4)} nm → ${F(l.rho * l.t * 1e5, 4)} µg/cm², ${F(l.rho * l.t / l.M * 1e5, 4)} nmol/cm² of formula units`));
  const checks = [
    chk('Roughness small against the thinnest layer (continuous layers)', rq < thin / 3, `rms roughness ${F(nm(rq), 3)} nm against thinnest layer ${F(nm(thin), 3)} nm`),
    chk('Top layer closed in the island-growth model', covTop > 0.99, `coverage ${F(covTop * 100, 4)} % at ${F(nm(top.t), 3)} nm; 99 % closure needs ${F(nm(h99), 3)} nm at this island density`),
    chk('Generated surface reproduces the Hurst exponent', abs(Hest - v.hurst) < 0.2, `H from the short-range structure function = ${F(Hest, 3)} against the input ${F(v.hurst, 3)} (64 × 64 grid, finite-size bias expected)`),
    chk('Scan spans several correlation lengths', v.scan >= 5 * v.xi, `scan ${F(nm(v.scan), 3)} nm against correlation length ${F(nm(v.xi), 3)} nm`),
    chk('Step coverage within [0, 1]', scCVD >= 0 && scCVD <= 1 && scPVD >= 0 && scPVD <= 1, `CVD-type ${F(scCVD, 3)}, line-of-sight ${F(scPVD, 3)}`),
  ];
  if (v.xi > v.scan / 2) warnings.push('The correlation length exceeds half the scan size and was limited to it: the height map does not sample the long-wavelength roughness.');
  if (rq >= thin / 3) warnings.push(`The rms roughness (${F(nm(rq), 3)} nm) is comparable with the thinnest layer (${F(nm(thin), 3)} nm): that layer is likely to be discontinuous or pin-holed.`);
  return {
    summary: `The ${layers.length}-layer stack (${layers.map(l => `${l.pretty} ${F(nm(l.t), 3)} nm`).join(' / ')}) is ${F(nm(total), 4)} nm thick with an areal mass of ${F(massA * 1e5, 4)} µg/cm². In the island-growth model the top layer covers ${F(covTop * 100, 4)} % of the surface and closes (99 %) at ${F(nm(h99), 3)} nm. The self-affine surface has rms roughness ${F(nm(rq), 3)} nm, average roughness ${F(nm(ra), 3)} nm and an area ratio of ${F(areaRatio, 5)}. On a feature of aspect ratio ${F(AR, 3)}, the step coverage is ${F(scCVD * 100, 3)} % for a CVD-type process and ${F(scPVD * 100, 3)} % for line-of-sight deposition.`,
    kpis: [
      kpi('total', 'Total thickness', total, 'm', `Sum of the layer thicknesses = ${F(nm(total), 4)} nm over ${layers.length} layers.`),
      kpi('massA', 'Areal mass', massA * 1e5, 'µg/cm²', `Σρᵢtᵢ = ${F(massA, 4)} kg/m²; this is what a quartz-crystal microbalance or weighing measures and what sets precursor consumption per unit area.`),
      kpi('cov', 'Top-layer coverage (island growth)', covTop, '–', `Hemispherical islands at ${F(Ni, 3)} 1/m² with random overlap (Avrami): θ = 1 − exp(−Nπr²) = ${F(covTop * 100, 4)} % at the nominal thickness ${F(nm(top.t), 3)} nm. Layer-by-layer growth would need only ${F(top.delta * 1e10, 3)} Å for a full monolayer.`, covTop > 0.99),
      kpi('h99', 'Thickness for 99 % closure', h99, 'm', `Nominal (mass-equivalent) thickness at which the island model reaches 99 % coverage = ${F(nm(h99), 4)} nm; it scales as N^(−1/2), so higher nucleation density gives thinner continuous films.`),
      kpi('rq', 'RMS roughness R_q', rq, 'm', `Root-mean-square height of the generated surface, scaled to the requested ${F(nm(rq), 4)} nm; it is ${F(rq / total * 100, 3)} % of the total thickness.`, rq < thin / 3),
      kpi('ra', 'Average roughness R_a', ra, 'm', `Mean absolute height = ${F(nm(ra), 4)} nm; for a Gaussian surface R_a/R_q = 0.80 (here ${F(ra / rq, 3)}). Peak-to-valley over the scan: ${F(nm(hmax - hmin), 4)} nm; skewness ${F(skew, 3)}, kurtosis ${F(kurt, 3)}.`),
      kpi('slope', 'RMS surface slope', rmsSlope, '–', `Root-mean-square gradient from forward differences at ${F(nm(dx), 3)} nm pixel size = ${F(rmsSlope, 4)} (${F(Math.atan(rmsSlope) * 180 / PI, 3)}°). It controls light scattering and contact mechanics and depends on the pixel size for a self-affine surface.`),
      kpi('areaR', 'True-to-projected area ratio', areaRatio, '–', `⟨(1 + |∇h|²)^(1/2)⟩ = ${F(areaRatio, 5)}: the roughness adds ${F((areaRatio - 1) * 100, 3)} % surface area at this resolution.`),
      kpi('Hest', 'Hurst exponent recovered', Hest, '–', `Half the log–log slope of the height–height correlation function at short range = ${F(Hest, 3)} (input ${F(v.hurst, 3)}); the lateral correlation length recovered from the 1 − 1/e point is ${F(nm(xiEst), 3)} nm.`, abs(Hest - v.hurst) < 0.2),
      kpi('scCVD', 'Step coverage, CVD-type', scCVD, '–', `Bottom-to-top thickness ratio in a cylindrical via from first-order Knudsen diffusion–reaction, 1/cosh(AR·√(3s)) with Thiele modulus ${F(thiele, 3)} = ${F(scCVD * 100, 3)} %.`, scCVD > 0.9),
      kpi('scPVD', 'Step coverage, line of sight', scPVD, '–', `Flux reaching the centre of a trench bottom from a diffuse (cosine) source relative to the open surface, 1/(1 + 4AR²)^(1/2) = ${F(scPVD * 100, 3)} %.`, scPVD > 0.5),
      kpi('ald', 'ALD exposure multiplier', aldMult, '–', `Precursor exposure needed to saturate the feature relative to a flat surface, 1 + (19/4)AR + (3/2)AR² = ${F(aldMult, 4)} (kinetic model for a hole in molecular flow). ALD stays conformal if the dose and purge times are scaled by this factor.`),
    ],
    series: [
      { key: 'cov', title: 'Coverage vs nominal thickness', xlabel: 'nominal thickness (nm)', ylabel: 'coverage (–)', logx: true, lines: [{ name: 'island growth (Volmer–Weber)', x: hs.map(nm), y: hs.map(cov) }, { name: 'layer by layer (Frank–van der Merwe)', x: hs.map(nm), y: hs.map(q => min(1, q / top.delta)) }], explain: 'Layer-by-layer growth closes within one monolayer, whereas island growth leaves bare substrate until islands coalesce. The offset between the curves is the thickness penalty of poor wetting.' },
      { key: 'hhcf', title: 'Height–height correlation function', xlabel: 'lateral distance r (nm)', ylabel: 'g(r) = ⟨(h(x+r) − h(x))²⟩ (nm²)', logx: true, logy: true, lines: [{ name: 'generated surface', x: lags.map(nm), y: gr.map(q => max(1e-30, q * 1e18)) }, { name: 'saturation 2R_q²', x: [nm(lags[0]), nm(lags[lags.length - 1])], y: [2 * rq * rq * 1e18, 2 * rq * rq * 1e18] }], explain: 'For a self-affine surface g(r) rises as r^(2H) and saturates at 2R_q² beyond the correlation length. The initial slope gives the Hurst exponent and the knee the correlation length.' },
      { key: 'sc', title: 'Step coverage vs aspect ratio', xlabel: 'aspect ratio (depth/width)', ylabel: 'bottom/top thickness ratio (–)', logx: true, lines: [{ name: `CVD-type, s = ${F(v.sticking, 2)}`, x: ars, y: ars.map(a => 1 / Math.cosh(min(700, a * sqrt(3 * v.sticking)))) }, { name: 'line of sight (PVD)', x: ars, y: ars.map(a => 1 / sqrt(1 + 4 * a * a)) }], explain: 'Conformality falls with aspect ratio. A low sticking coefficient keeps CVD-type coatings conformal to high aspect ratio; line-of-sight processes lose bottom coverage almost immediately.' },
    ],
    maps: [{ key: 'height', title: 'Surface height map (self-affine)', nx: n, ny: n, data: Float64Array.from(h, q => q * 1e9), x0: 0, x1: nm(v.scan), y0: 0, y1: nm(v.scan), xlabel: 'x (nm)', ylabel: 'y (nm)', unit: 'nm', diverging: true, explain: 'Synthetic topography with the requested rms roughness, Hurst exponent and correlation length, generated by Fourier filtering of seeded Gaussian noise (periodic). Compare its statistics with an AFM scan of the real film.' }],
    tables: [{ key: 'layers', title: 'Layer stack', columns: ['Layer', 'Formula', 'Thickness (nm)', 'Density (kg/m³)', 'Areal mass (µg/cm²)', 'Monolayers (t/δ)'], rows: layers.map(l => [l.name, l.pretty, +nm(l.t).toPrecision(5), l.rho, +(l.rho * l.t * 1e5).toPrecision(4), +(l.t / l.delta).toPrecision(4)]), explain: 'Each layer from the substrate upwards with its areal mass ρt and the equivalent number of atomic layers t/δ, where δ is the mean atomic spacing from the molar volume.' }],
    equations: [
      ...stackEq,
      eqn('formula', 'Stack totals', 't_tot = Σtᵢ;  m″ = Σρᵢtᵢ'),
      eqn('constitutive', 'Island-growth coverage (hemispherical caps, random overlap)', 'h = N(2/3)πr³;  θ = 1 − exp(−Nπr²)'),
      eqn('constitutive', 'Self-affine surface (power spectral density)', 'C(q) ∝ (1 + (qξ)²)^−(1+H);  g(r) ∝ r^(2H) for r ≪ ξ;  g(r) → 2R_q² for r ≫ ξ'),
      eqn('governing', 'Diffusion–reaction in a high-aspect-ratio via (Knudsen regime)', 'D_K d²n/dz² = (s v̄/d) n  →  n(H)/n(0) = 1/cosh(Φ),  Φ = AR·√(3s)  (D_K = v̄d/3)'),
      eqn('formula', 'Line-of-sight flux at a trench bottom', 'J_bottom/J_top = sin(arctan(1/(2AR))) = 1/(1 + 4AR²)^(1/2)'),
    ],
    geometry: { kind: 'voxels', nx: nv, ny: nv, nz, data, voxel: nm(v.scan) / nv, units: 'nm' },
    checks, warnings,
    recommendations: [
      rec('high', covTop < 0.99 ? `The ${F(nm(top.t), 3)} nm top layer is not closed in the island-growth model (${F(covTop * 100, 3)} %): raise the nucleation density (surface activation, seed or adhesion layer, lower temperature) or deposit at least ${F(nm(h99), 3)} nm.` : `The top layer is continuous with margin (closure at ${F(nm(h99), 3)} nm against ${F(nm(top.t), 3)} nm deposited); the thickness could be reduced towards ${F(nm(h99) * 1.5, 3)} nm if only continuity is required.`, 'Island-growth coverage model.'),
      rec('medium', scCVD < 0.9 ? `Step coverage on the aspect-ratio-${F(AR, 3)} feature is ${F(scCVD * 100, 3)} %: lower the sticking coefficient (lower temperature, less reactive precursor) or switch to ALD with doses scaled by about ${F(aldMult, 3)} × the flat-surface exposure.` : `The CVD-type process is conformal on this feature (${F(scCVD * 100, 3)} %); the sticking coefficient could rise to about ${F(min(1, (Math.acosh(1 / 0.9) / AR) ** 2 / 3), 2)} before coverage drops below 90 %.`, 'Thiele-modulus model of diffusion–reaction in a via.'),
      rec('medium', `Verify the roughness by AFM over a scan of at least ${F(nm(10 * v.xi), 3)} nm (ten correlation lengths) and compare R_q, the Hurst exponent and the correlation length with the values used here.`, 'Finite-size requirements of self-affine roughness statistics.'),
      rec('low', `The stack contains ${F(massA * 1e5, 3)} µg/cm² of material; because island-free (well-nucleated) growth closes films at a fraction of the thickness, improving nucleation is the most effective way to cut precursor consumption, deposition time and energy per coated area.`, 'Areal mass and closure thickness; resource efficiency.', true),
    ],
    publish: { 'geom.thickness': total, 'film.areal_mass': massA, 'geom.roughness_rms': rq, 'geom.coverage': covTop, 'film.step_coverage': scCVD, 'geom.hurst': v.hurst, 'material.formula': top.f, 'prop.density': top.rho },
    exports: [
      { name: 'height_map.csv', mime: 'text/csv', content: `# height in nm, ${n} x ${n} pixels, pixel ${nm(dx).toPrecision(5)} nm\n` + Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (h[j * n + i] * 1e9).toPrecision(5)).join(',')).join('\n') + '\n', description: 'Synthetic self-affine surface topography.' },
      { name: 'layer_stack.json', mime: 'application/json', content: JSON.stringify({ layers: layers.map(l => ({ name: l.name, formula: l.f, thickness_m: l.t, density_kg_m3: l.rho, areal_mass_kg_m2: l.rho * l.t })), total_thickness_m: total, areal_mass_kg_m2: massA, roughness: { rq_m: rq, ra_m: ra, hurst: v.hurst, correlation_length_m: xi } }, null, 1), description: 'Layer stack with areal masses and roughness parameters.' },
    ],
    numerics: { method: 'Spectral synthesis of a self-affine surface (two-dimensional inverse FFT of filtered seeded Gaussian noise, 64 × 64); closed-form coverage and step-coverage models', dof: n * n, steps: 2 * n, iterations: 1, residual: abs(Hest - v.hurst), notes: 'Periodic surface; the rms height is imposed by scaling. The voxel picture exaggerates the roughness when it is small relative to the thickness.' },
  };
}

const RUN = { single: runSingle, population: runPopulation, rods: runRods, porous: runPorous, aggregate: runAggregate, film: runFilm };

export default {
  id: 4,
  tagline: 'Shape descriptors, size distributions, nanorod arrays, porous microstructures, aggregates and thin-film morphology with generated geometry and tolerance maps.',
  inputs,
  presets: [
    { name: 'ZnO nanorod array on a seeded substrate', note: 'Hydrothermally grown wurtzite nanorods: 80 nm × 1.5 µm at 40 rods per µm² with a few degrees of tilt.', values: { mode: 'rods', material: 'ZnO', d: 8e-8, L: 1.5e-6, nDens: 4e13, cvL: 15, cvD: 12, tilt: 6, arrangement: 'random', section: 'hex' } },
    { name: 'Au nanoparticle population, 15 nm citrate sol', note: 'Log-normal gold colloid typical of citrate reduction; shows how number, surface and volume means differ.', values: { mode: 'population', popSource: 'generate', material: 'Au', dist: 'lognormal', d: 1.5e-8, gsd: 1.15, nSample: 3000 } },
    { name: 'Nanoporous medium, porosity 0.25', note: 'Sintered-powder type representative volume with 20 nm grains for transport suites.', values: { mode: 'porous', material: 'TiO2', pType: 'randomSpheres', phi: 0.25, feat: 2e-8, featVox: 7, pN: 28 } },
    { name: 'ZnO sphere, 10 nm (reference: 107 m²/g)', note: 'Hand-check case: SSA = 6/(ρd) = 6/(5606 × 10⁻⁸) m²/kg.', values: { mode: 'single', material: 'ZnO', shape: 'sphere', d: 1e-8, tolPct: 5 } },
    { name: 'Au@SiO₂ core–shell, 3 nm core', note: 'Small enough for an atomistic core with a continuum silica shell description.', values: { mode: 'single', material: 'Au', shape: 'coreshell', d: 3e-9, tw: 2e-9, shellMat: 'SiO2' } },
    { name: 'Fumed-silica-like DLA aggregate', note: 'Fractal aggregate of 12 nm primary particles.', values: { mode: 'aggregate', aggType: 'dla', material: 'SiO2', d: 1.2e-8, nPart: 1000, stick: 1 } },
    { name: 'Fit measured TEM diameters', note: 'Fits log-normal, normal and Weibull models to the diameter table and ranks them by KS statistic.', values: { mode: 'population', popSource: 'fit', material: 'Au', dist: 'lognormal' } },
    { name: 'Oxide multilayer on a high-aspect-ratio trench', note: 'Seed / barrier / functional film stack with self-affine roughness and conformality estimate.', values: { mode: 'film', rmsT: 2e-9, hurst: 0.7, xi: 5e-8, scan: 1e-6, featAR: 30, sticking: 0.005 } },
  ],
  run(inp) { const v = clean(inp); return RUN[v.mode](v); },
};
