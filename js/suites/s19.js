// Suite 19 — Large-scale component and system analysis.
// Finite lattice panels by 2-D frame finite elements (banded direct solver), homogenised Kirchhoff plates, N-DOF modal and
// harmonic analysis of resonator chains, finite phased-array / metasurface apertures, lumped reactor thermal networks with a
// 3-D simulation cost estimate, and Weibull reliability block diagrams.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => note ? { kind, label, text, note } : { kind, label, text };
const rec = (priority, text, basis, sustainability) => sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis };
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const line = (name, x, y) => { const X = [], Y = [], st = Math.max(1, Math.ceil(x.length / 400)); for (let i = 0; i < x.length; i += st) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) { X.push(x[i]); Y.push(y[i]); } if (!X.length) { X.push(0); Y.push(0); } return { name, x: X, y: Y }; };
const dB = v => 10 * Math.log10(Math.max(v, 1e-30));

const MODES = [
  { value: 'lattice_fe', label: 'Finite lattice panel (frame finite elements)' },
  { value: 'plate', label: 'Homogenised plate / panel component' },
  { value: 'modal', label: 'Resonator array: modes, harmonic response, finite band gap' },
  { value: 'array', label: 'Finite metasurface / antenna aperture' },
  { value: 'reactor', label: 'Scaled reactor: thermal network and 3-D simulation cost' },
  { value: 'reliability', label: 'System reliability block diagram' },
];
const WALLS = { steel: { name: 'stainless steel 316L', k: 16, rho: 8000, cp: 500 }, titanium: { name: 'titanium', k: 21, rho: 4500, cp: 520 }, nickel: { name: 'nickel alloy C-276', k: 10, rho: 8900, cp: 430 }, glass: { name: 'borosilicate glass', k: 1.1, rho: 2230, cp: 830 } };
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'lattice_fe', options: MODES, group: 'Model', help: 'Selects the component- or system-level model. Each analysis uses the inputs shown for it.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Model', showIf: { mode: ['array', 'reliability'] }, help: 'Seed for manufacturing-error realisations and Monte-Carlo sampling.' },
  // lattice FE
  { key: 'ctype', label: 'Cell type', type: 'select', value: 'reentrant', options: [{ value: 'reentrant', label: 'Re-entrant honeycomb (auxetic)' }, { value: 'honeycomb', label: 'Regular hexagonal honeycomb' }, { value: 'triangular', label: 'Triangular (stretch-dominated)' }, { value: 'square', label: 'Square grid' }], group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Unit-cell topology of the planar lattice. Re-entrant cells give a negative Poisson\'s ratio; triangular cells are stretch-dominated and stiff.' },
  { key: 'ncx', label: 'Cells along x', value: 8, min: 2, max: 20, step: 1, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Number of unit cells across the panel width. More cells reduce edge effects relative to the homogenised prediction.' },
  { key: 'ncy', label: 'Cells along y', value: 5, min: 2, max: 12, step: 1, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Number of unit cells along the panel height (loading direction for tension, compression and shear).' },
  { key: 'cell', label: 'Cell width', unit: 'm', value: 0.01, min: 1e-7, max: 1, log: true, group: 'Lattice panel', bus: 'geom.cell_size', showIf: { mode: ['lattice_fe'] }, help: 'Width of one unit cell along x; member lengths follow from the cell type and angle.' },
  { key: 't_strut', label: 'Member in-plane thickness', unit: 'm', value: 5e-4, min: 1e-9, max: 0.1, log: true, group: 'Lattice panel', bus: 'geom.strut', showIf: { mode: ['lattice_fe'] }, help: 'In-plane thickness of the cell walls; bending stiffness scales with its cube.' },
  { key: 'depth', label: 'Out-of-plane depth', unit: 'm', value: 0.01, min: 1e-8, max: 1, log: true, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Extrusion depth of the planar lattice (plane-stress members).' },
  { key: 'th_re', label: 'Re-entrant angle magnitude', unit: '°', value: 20, min: 5, max: 40, step: 1, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Inward inclination of the oblique walls of the re-entrant cell; larger angles increase the magnitude of the negative Poisson\'s ratio along y.' },
  { key: 'h_l', label: 'Wall ratio h/l (re-entrant)', value: 2, min: 1.2, max: 4, step: 0.1, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Ratio of vertical to oblique wall length of the re-entrant cell; it must exceed 2·sin θ for the cell not to self-intersect.' },
  { key: 'E_s', label: 'Member Young\'s modulus', unit: 'Pa', value: 7e10, min: 1e6, max: 1.2e12, log: true, group: 'Lattice panel', bus: 'prop.E', showIf: { mode: ['lattice_fe'] }, help: 'Young\'s modulus of the solid from which the lattice is made.' },
  { key: 'sy', label: 'Member yield strength', unit: 'Pa', value: 2.5e8, min: 1e5, max: 2e10, log: true, group: 'Lattice panel', bus: 'prop.yield', showIf: { mode: ['lattice_fe', 'plate'] }, help: 'Yield (or allowable) stress used for the safety factor.' },
  { key: 'lcase', label: 'Load case', type: 'select', value: 'tension', options: [{ value: 'tension', label: 'Uniaxial tension along y' }, { value: 'compression', label: 'Uniaxial compression along y' }, { value: 'shear', label: 'Simple shear (top edge displaced along x)' }, { value: 'bend3', label: 'Three-point bending' }], group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Supports and loading. Tension, compression and shear are displacement-controlled; bending applies a central force.' },
  { key: 'strain', label: 'Applied nominal strain', unit: '%', value: 0.2, min: 0.001, max: 5, log: true, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Edge displacement divided by panel height for tension, compression and shear (linear analysis, so results scale with it).' },
  { key: 'F_load', label: 'Central force (bending)', unit: 'N', value: 50, min: 1e-9, max: 1e6, log: true, group: 'Lattice panel', showIf: { mode: ['lattice_fe'] }, help: 'Downward force at mid-span for the three-point bending case.' },
  // plate
  { key: 'pa', label: 'Plate length a', unit: 'm', value: 0.5, min: 1e-6, max: 20, log: true, group: 'Plate', showIf: { mode: ['plate'] }, help: 'Plate dimension along x (simply supported on all four edges).' },
  { key: 'pb', label: 'Plate width b', unit: 'm', value: 0.4, min: 1e-6, max: 20, log: true, group: 'Plate', showIf: { mode: ['plate'] }, help: 'Plate dimension along y.' },
  { key: 'ph', label: 'Plate thickness', unit: 'm', value: 0.01, min: 1e-9, max: 0.5, log: true, group: 'Plate', bus: 'geom.thickness', showIf: { mode: ['plate'] }, help: 'Thickness of the (homogenised) panel; bending stiffness scales with its cube.' },
  { key: 'E_p', label: 'Effective Young\'s modulus', unit: 'Pa', value: 2e9, min: 1e4, max: 1.2e12, log: true, group: 'Plate', bus: 'eff.E', showIf: { mode: ['plate'] }, help: 'Homogenised in-plane modulus of the panel material (for example from the lattice or homogenisation suites).' },
  { key: 'nu_p', label: 'Effective Poisson\'s ratio', value: -0.4, min: -0.95, max: 0.49, step: 0.01, group: 'Plate', bus: 'eff.nu', showIf: { mode: ['plate'] }, help: 'Homogenised Poisson\'s ratio. Negative (auxetic) values are admissible down to −1 and produce dome-shaped (synclastic) bending.' },
  { key: 'rho_p', label: 'Effective density', unit: 'kg/m³', value: 400, min: 1, max: 23000, log: true, group: 'Plate', bus: 'eff.density', showIf: { mode: ['plate'] }, help: 'Homogenised mass density of the panel.' },
  { key: 'q_p', label: 'Uniform pressure', unit: 'Pa', value: 2000, min: 1e-3, max: 1e8, log: true, group: 'Plate', showIf: { mode: ['plate'] }, help: 'Uniform transverse pressure on the plate.' },
  { key: 'Nx_p', label: 'In-plane compressive load', unit: 'N/m', value: 1e4, min: 0, max: 1e9, group: 'Plate', showIf: { mode: ['plate'] }, help: 'Uniform edge compression along x per unit width, compared with the buckling load.' },
  // modal
  { key: 'n_cells', label: 'Number of unit cells', value: 10, min: 2, max: 40, step: 1, group: 'Resonator chain', showIf: { mode: ['modal'] }, help: 'Number of two-mass unit cells in the finite chain; attenuation in the gap grows exponentially with it.' },
  { key: 'm1', label: 'Mass 1', unit: 'kg', value: 0.01, min: 1e-15, max: 1e3, log: true, group: 'Resonator chain', showIf: { mode: ['modal'] }, help: 'First mass of the diatomic unit cell.' },
  { key: 'm2', label: 'Mass 2', unit: 'kg', value: 0.03, min: 1e-15, max: 1e3, log: true, group: 'Resonator chain', showIf: { mode: ['modal'] }, help: 'Second mass of the unit cell; the mass contrast opens the Bragg band gap.' },
  { key: 'k_s', label: 'Coupling stiffness', unit: 'N/m', value: 1e5, min: 1e-6, max: 1e12, log: true, group: 'Resonator chain', showIf: { mode: ['modal'] }, help: 'Stiffness of the springs connecting neighbouring masses.' },
  { key: 'eta', label: 'Structural loss factor', value: 0.01, min: 1e-5, max: 0.5, log: true, group: 'Resonator chain', showIf: { mode: ['modal'] }, help: 'Hysteretic damping of the springs, k(1 + iη); the modal quality factor is about 1/η.' },
  { key: 'gap_lo', label: 'Bloch gap lower edge (from unit-cell analysis)', unit: 'Hz', value: 0, min: 0, max: 1e12, group: 'Resonator chain', bus: 'res.bandgap_lo', showIf: { mode: ['modal'] }, help: 'Lower band-gap edge of the infinite structure from another suite; 0 uses the analytical diatomic-chain value.' },
  { key: 'gap_hi', label: 'Bloch gap upper edge (from unit-cell analysis)', unit: 'Hz', value: 0, min: 0, max: 1e12, group: 'Resonator chain', bus: 'res.bandgap_hi', showIf: { mode: ['modal'] }, help: 'Upper band-gap edge of the infinite structure; 0 uses the analytical diatomic-chain value.' },
  // array
  { key: 'nx', label: 'Elements along x', value: 16, min: 2, max: 40, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Number of radiating or scattering elements along the steering axis.' },
  { key: 'ny', label: 'Elements along y', value: 16, min: 1, max: 40, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Number of elements along the orthogonal axis (affects gain and error averaging).' },
  { key: 'pitch', label: 'Element pitch', unit: 'm', value: 0.015, min: 1e-8, max: 1, log: true, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Centre-to-centre element spacing. It must stay below λ/(1 + |sin θ₀|) to avoid grating lobes.' },
  { key: 'freq', label: 'Operating frequency', unit: 'Hz', value: 1e10, min: 1e6, max: 1e15, log: true, group: 'Aperture', bus: 'res.f0', showIf: { mode: ['array'] }, help: 'Operating frequency; sets the free-space wavelength.' },
  { key: 'steer', label: 'Steering angle', unit: '°', value: 20, min: -75, max: 75, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Target beam direction from broadside in the x–z plane, produced by a progressive phase gradient.' },
  { key: 'taper', label: 'Amplitude taper', type: 'select', value: 'uniform', options: [{ value: 'uniform', label: 'Uniform' }, { value: 'cosine', label: 'Cosine on a 0.2 pedestal' }, { value: 'triangular', label: 'Triangular' }], group: 'Aperture', showIf: { mode: ['array'] }, help: 'Aperture illumination: tapering lowers sidelobes at the cost of beamwidth and gain.' },
  { key: 'phase_err', label: 'Random phase error (rms)', unit: '°', value: 15, min: 0, max: 90, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Element-to-element phase error from manufacturing variation (geometry, material, bias).' },
  { key: 'amp_err', label: 'Random amplitude error (rms)', unit: '%', value: 10, min: 0, max: 50, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Relative amplitude error of each element.' },
  { key: 'bits', label: 'Phase quantisation bits', value: 3, min: 0, max: 8, step: 1, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Number of phase-control bits per element (0 = continuous phase). One- and two-bit metasurfaces lose several dB.' },
  { key: 'fail', label: 'Failed elements', unit: '%', value: 2, min: 0, max: 50, step: 0.5, group: 'Aperture', showIf: { mode: ['array'] }, help: 'Fraction of dead elements (fabrication yield or in-service failures).' },
  // reactor
  { key: 'V_r', label: 'Reactor liquid volume', unit: 'm³', value: 0.05, min: 1e-5, max: 50, log: true, group: 'Reactor', bus: 'proc.volume', showIf: { mode: ['reactor'] }, help: 'Liquid volume of the scaled hydrothermal or solution reactor.' },
  { key: 'T_set', label: 'Process temperature', unit: '°C', value: 90, min: 30, max: 350, group: 'Reactor', bus: 'proc.T_C', showIf: { mode: ['reactor'] }, help: 'Target liquid temperature.' },
  { key: 'P_heat', label: 'Heater power', unit: 'W', value: 5000, min: 1, max: 5e6, log: true, group: 'Reactor', bus: 'proc.power', showIf: { mode: ['reactor'] }, help: 'Installed heating power applied to the vessel wall.' },
  { key: 'wall_mat', label: 'Vessel wall material', type: 'select', value: 'steel', options: Object.entries(WALLS).map(([value, w]) => ({ value, label: w.name })), group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Wall material: sets wall conductivity and heat capacity.' },
  { key: 'wall_mm', label: 'Wall thickness', unit: 'mm', value: 8, min: 0.5, max: 200, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Vessel wall thickness; it adds thermal mass and carries the pressure load.' },
  { key: 'ins_mm', label: 'Insulation thickness', unit: 'mm', value: 50, min: 0, max: 500, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Thickness of external insulation (conductivity 0.04 W/m/K assumed).' },
  { key: 'h_int', label: 'Internal heat-transfer coefficient', unit: 'W/m²/K', value: 500, min: 20, max: 10000, log: true, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Wall-to-liquid coefficient: about 200–1000 for natural convection in water, several thousand when stirred.' },
  { key: 'h_ext', label: 'External heat-transfer coefficient', unit: 'W/m²/K', value: 10, min: 1, max: 200, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Combined convection and radiation coefficient on the outer surface.' },
  { key: 'T_amb', label: 'Ambient temperature', unit: '°C', value: 25, min: -20, max: 60, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Surrounding air temperature and initial temperature of the charge.' },
  { key: 'HD', label: 'Height / diameter ratio', value: 1.5, min: 0.3, max: 10, step: 0.1, group: 'Reactor', showIf: { mode: ['reactor'] }, help: 'Aspect ratio of the cylindrical vessel; it sets the surface-to-volume ratio.' },
  { key: 'mesh_mm', label: '3-D mesh size', unit: 'mm', value: 2, min: 0.05, max: 100, log: true, group: 'Simulation cost', showIf: { mode: ['reactor'] }, help: 'Cell size of the equivalent full 3-D simulation (boundary layers usually need about 1 mm or finer).' },
  { key: 't_proc', label: 'Simulated process time', unit: 's', value: 3600, min: 1, max: 1e6, log: true, group: 'Simulation cost', bus: 'proc.time', showIf: { mode: ['reactor'] }, help: 'Physical time the transient 3-D simulation would have to cover.' },
  { key: 'cores', label: 'Compute cores', value: 32, min: 1, max: 100000, log: true, sens: false, group: 'Simulation cost', showIf: { mode: ['reactor'] }, help: 'Number of processor cores available (70 % parallel efficiency assumed).' },
  { key: 'cpu_us', label: 'Solver cost per DOF per step', unit: 'µs', value: 2, min: 0.05, max: 200, log: true, sens: false, group: 'Simulation cost', showIf: { mode: ['reactor'] }, help: 'Core time per degree of freedom per time step of the 3-D solver (about 1–10 µs for implicit finite-volume codes); calibrate on a short test run.' },
  // reliability
  { key: 'n_ser', label: 'Stages in series', value: 100, min: 1, max: 500, step: 1, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Number of stages that must all work (unit cells in a load path, elements in a chain).' },
  { key: 'm_par', label: 'Redundant components per stage', value: 1, min: 1, max: 5, step: 1, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Identical components in parallel within each stage; the stage fails only when all of them have failed.' },
  { key: 'beta_w', label: 'Component Weibull shape β', value: 2.5, min: 0.5, max: 8, step: 0.1, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Shape parameter of the component life distribution: below 1 infant mortality, about 1 random failures, above 1 wear-out.' },
  { key: 'eta_h', label: 'Component characteristic life η', unit: 'h', value: 2e5, min: 10, max: 1e8, log: true, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Weibull scale: time at which 63.2 % of components have failed.' },
  { key: 'beta_c', label: 'Common element shape β', value: 1, min: 0.5, max: 8, step: 0.1, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Weibull shape of a single shared element in series with all stages (controller, joint, power supply).' },
  { key: 'eta_c', label: 'Common element characteristic life', unit: 'h', value: 5e5, min: 10, max: 1e9, log: true, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Weibull scale of the shared element.' },
  { key: 't_miss', label: 'Mission time', unit: 'h', value: 8760, min: 1, max: 1e7, log: true, group: 'Reliability', showIf: { mode: ['reliability'] }, help: 'Service time at which reliability is evaluated (8760 h = one year of continuous operation).' },
];
function read(inp) {
  const p = {};
  for (const d of inputs) {
    let v = inp ? inp[d.key] : undefined;
    if (d.type === 'select') p[d.key] = d.options.some(o => o.value === v) ? v : d.value;
    else { v = Number(v); p[d.key] = N.clamp(Number.isFinite(v) ? v : d.value, d.min, d.max); }
  }
  return p;
}
function finish(r) {
  r.warnings = r.warnings || [];
  for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { k.value = 0; r.warnings.push(`KPI "${k.label}" could not be evaluated for these inputs and is reported as 0.`); }
  for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
  if (!Object.keys(r.publish).length) r.publish['mech.status'] = 'not evaluated';
  return r;
}

// ---------- lattice panel geometry ----------
function buildLattice(type, Nx, Ny, a, th, hl) {
  let nodes = [], edges = []; const id = new Map(), nd = (x, y) => { const k = Math.round(x / a * 1e5) + '_' + Math.round(y / a * 1e5); let i = id.get(k); if (i === undefined) { i = nodes.length; id.set(k, i); nodes.push([x, y]); } return i; };
  let info;
  if (type === 'square') { for (let j = 0; j <= Ny; j++) for (let i = 0; i <= Nx; i++) { const n0 = nd(i * a, j * a); if (i < Nx) edges.push([n0, nd((i + 1) * a, j * a)]); if (j < Ny) edges.push([n0, nd(i * a, (j + 1) * a)]); } info = { l: a, cw: a, ch: a }; }
  else if (type === 'triangular') { const rh = a * Math.sqrt(3) / 2; for (let j = 0; j <= Ny; j++) for (let i = 0; i <= Nx - (j % 2); i++) nd(i * a + (j % 2) * a / 2, j * rh); for (let i = 0; i < nodes.length; i++) for (let j = i + 1; j < nodes.length; j++) if (Math.abs(Math.hypot(nodes[i][0] - nodes[j][0], nodes[i][1] - nodes[j][1]) - a) < 1e-6 * a) edges.push([i, j]); info = { l: a, cw: a, ch: rh }; }
  else {
    const c = Math.cos(th), s = Math.sin(th), l = a / (2 * c), h = hl * l, ch = 2 * (h + l * s);
    for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) { const x = i * a, y = j * ch, A = nd(x, y), B = nd(x, y + h), C = nd(x + l * c, y + h + l * s), D = nd(x + l * c, y + 2 * h + l * s); edges.push([A, B], [B, C], [C, D], [B, nd(x - l * c, y + h + l * s)], [D, nd(x, y + ch)], [D, nd(x + a, y + ch)]); }
    const key = e => Math.min(e[0], e[1]) + '_' + Math.max(e[0], e[1]), seen = new Set(); edges = edges.filter(e => { const k = key(e); if (seen.has(k)) return false; seen.add(k); return true; });
    const yTop = Ny * ch; let changed = true, alive = nodes.map(() => true);
    while (changed) { changed = false; const deg = new Int32Array(nodes.length); for (const e of edges) { deg[e[0]]++; deg[e[1]]++; } for (let i = 0; i < nodes.length; i++) if (alive[i] && deg[i] <= 1 && nodes[i][1] > 1e-9 * a && nodes[i][1] < yTop - 1e-9 * a) { alive[i] = false; changed = true; } if (changed) edges = edges.filter(e => alive[e[0]] && alive[e[1]]); }
    const deg = new Int32Array(nodes.length); for (const e of edges) { deg[e[0]]++; deg[e[1]]++; }
    const map = new Int32Array(nodes.length).fill(-1), nn = []; nodes.forEach((q, i) => { if (alive[i] && deg[i] > 0) { map[i] = nn.length; nn.push(q); } }); edges = edges.map(e => [map[e[0]], map[e[1]]]); nodes = nn; info = { l, h, cw: a, ch, c, s };
  }
  return { nodes, edges, info };
}
// ---------- 2-D frame FE with banded Cholesky ----------
function frameSolve(nodes, edges, E, A, I, fix, load) { // non-dimensionalised by a reference member length so that translational and rotational stiffness are commensurate
  const e0 = edges[0], Lr = Math.hypot(nodes[e0[0]][0] - nodes[e0[1]][0], nodes[e0[0]][1] - nodes[e0[1]][1]) || 1;
  const S = frameCore(nodes.map(q => [q[0] / Lr, q[1] / Lr]), edges, E, A / Lr ** 2, I / Lr ** 4, fix.map(([n, c, v]) => [n, c, c < 2 ? v / Lr : v]), load.map(([n, c, v]) => [n, c, c < 2 ? v / Lr ** 2 : v / Lr ** 3]));
  S.disp = S.disp.map(d => [d[0] * Lr, d[1] * Lr, d[2]]); S.reac = S.reac.map(r => [r[0] * Lr ** 2, r[1] * Lr ** 2, r[2] * Lr ** 3]); S.res = S.res.map(r => ({ N: r.N * Lr ** 2, M1: r.M1 * Lr ** 3, M2: r.M2 * Lr ** 3, L: r.L * Lr }));
  return S;
}
function frameCore(nodes, edges, E, A, I, fix, load) {
  const nN = nodes.length, xs = nodes.map(q => q[0]), ys = nodes.map(q => q[1]), sc = 1e-7 * (N.maxOf(xs) - N.minOf(xs) + N.maxOf(ys) - N.minOf(ys));
  const ranking = byX => { const o = nodes.map((_, i) => i).sort((p, q) => byX ? (Math.abs(xs[p] - xs[q]) > sc ? xs[p] - xs[q] : ys[p] - ys[q]) : (Math.abs(ys[p] - ys[q]) > sc ? ys[p] - ys[q] : xs[p] - xs[q])), r = new Int32Array(nN); o.forEach((nd, k) => { r[nd] = k; }); let w = 0; for (const e of edges) w = Math.max(w, 3 * Math.abs(r[e[0]] - r[e[1]]) + 2); return { r, w }; };
  const rkX = ranking(true), rkY = ranking(false), rank = rkX.w <= rkY.w ? rkX.r : rkY.r, bw = Math.min(rkX.w, rkY.w), n = 3 * nN;
  const w = bw + 1, Kb = new Float64Array(n * w), f = new Float64Array(n), at = (i, j) => i * w + (j - i + bw);
  const ke = new Float64Array(36), dofs = new Int32Array(6), elem = [];
  for (const e of edges) {
    const dx = xs[e[1]] - xs[e[0]], dy = ys[e[1]] - ys[e[0]], L = Math.hypot(dx, dy), c = dx / L, s = dy / L, a = E * A / L, b = 12 * E * I / L ** 3, g = 6 * E * I / L ** 2, d4 = 4 * E * I / L, d2 = 2 * E * I / L;
    const k11 = a * c * c + b * s * s, k12 = (a - b) * c * s, k22 = a * s * s + b * c * c, rows = [[k11, k12, -g * s, -k11, -k12, -g * s], [k12, k22, g * c, -k12, -k22, g * c], [-g * s, g * c, d4, g * s, -g * c, d2], [-k11, -k12, g * s, k11, k12, g * s], [-k12, -k22, -g * c, k12, k22, -g * c], [-g * s, g * c, d2, g * s, -g * c, d4]];
    for (let q = 0; q < 3; q++) { dofs[q] = 3 * rank[e[0]] + q; dofs[3 + q] = 3 * rank[e[1]] + q; }
    for (let r = 0; r < 6; r++) for (let q = 0; q < 6; q++) { ke[r * 6 + q] = rows[r][q]; if (dofs[q] <= dofs[r]) Kb[at(dofs[r], dofs[q])] += rows[r][q]; }
    elem.push({ L, c, s, a, g, d4, d2, d: Array.from(dofs), ke: Float64Array.from(ke) });
  }
  for (const [nd, comp, val] of load) f[3 * rank[nd] + comp] += val;
  const fext = Float64Array.from(f), isFix = new Uint8Array(n), uc = new Float64Array(n);
  for (const [nd, comp, val] of fix) { const d = 3 * rank[nd] + comp; isFix[d] = 1; uc[d] = val; }
  let dmax = 0; const d0 = new Float64Array(n); for (let i = 0; i < n; i++) { d0[i] = Kb[at(i, i)]; dmax = Math.max(dmax, d0[i]); }
  for (let d = 0; d < n; d++) if (isFix[d]) { const v = uc[d]; for (let j = Math.max(0, d - bw); j < d; j++) { f[j] -= Kb[at(d, j)] * v; Kb[at(d, j)] = 0; } for (let i = d + 1; i <= Math.min(n - 1, d + bw); i++) { f[i] -= Kb[at(i, d)] * v; Kb[at(i, d)] = 0; } Kb[at(d, d)] = dmax; }
  for (let d = 0; d < n; d++) if (isFix[d]) f[d] = dmax * uc[d];
  let singular = false;
  for (let j = 0; j < n; j++) { // banded Cholesky, in place
    const rj = j * w + bw - j; let sdiag = Kb[rj + j]; const k0 = Math.max(0, j - bw); for (let k = k0; k < j; k++) { const v = Kb[rj + k]; sdiag -= v * v; }
    if (!(sdiag > 1e-11 * (isFix[j] ? dmax : d0[j]))) { singular = true; sdiag = 1e-9 * dmax; }
    const piv = Math.sqrt(sdiag); Kb[rj + j] = piv; const iMax = Math.min(n - 1, j + bw);
    for (let i = j + 1; i <= iMax; i++) { const ri = i * w + bw - i; let sv = Kb[ri + j]; if (sv === 0 && i - j > 2) { let any = false; for (let k = Math.max(0, i - bw); k < j; k++) if (Kb[ri + k] !== 0) { any = true; break; } if (!any) continue; } for (let k = Math.max(0, i - bw); k < j; k++) sv -= Kb[ri + k] * Kb[rj + k]; Kb[ri + j] = sv / piv; }
  }
  const u = Float64Array.from(f);
  for (let i = 0; i < n; i++) { let sv = u[i]; for (let k = Math.max(0, i - bw); k < i; k++) sv -= Kb[at(i, k)] * u[k]; u[i] = sv / Kb[at(i, i)]; }
  for (let i = n - 1; i >= 0; i--) { let sv = u[i]; for (let k = i + 1; k <= Math.min(n - 1, i + bw); k++) sv -= Kb[at(k, i)] * u[k]; u[i] = sv / Kb[at(i, i)]; }
  const fint = new Float64Array(n), res = [];
  for (const el of elem) { const ue = el.d.map(d => u[d]); for (let r = 0; r < 6; r++) { let sv = 0; for (let q = 0; q < 6; q++) sv += el.ke[r * 6 + q] * ue[q]; fint[el.d[r]] += sv; } const u1 = el.c * ue[0] + el.s * ue[1], v1 = -el.s * ue[0] + el.c * ue[1], u2 = el.c * ue[3] + el.s * ue[4], v2 = -el.s * ue[3] + el.c * ue[4]; res.push({ N: el.a * (u2 - u1), M1: el.g * (v1 - v2) + el.d4 * ue[2] + el.d2 * ue[5], M2: el.g * (v1 - v2) + el.d2 * ue[2] + el.d4 * ue[5], L: el.L }); }
  let rnorm = 0, fnorm = 0; const R = new Float64Array(n); for (let d = 0; d < n; d++) { if (isFix[d]) { R[d] = fint[d] - fext[d]; fnorm += R[d] * R[d]; } else { rnorm += (fint[d] - fext[d]) ** 2; fnorm += fext[d] * fext[d]; } }
  const disp = nodes.map((_, i) => [u[3 * rank[i]], u[3 * rank[i] + 1], u[3 * rank[i] + 2]]), reac = nodes.map((_, i) => [R[3 * rank[i]], R[3 * rank[i] + 1], R[3 * rank[i] + 2]]);
  return { disp, reac, res, dof: n, bw, residual: Math.sqrt(rnorm) / (Math.sqrt(fnorm) || 1), singular };
}

function latticeFE(p) {
  const warnings = [], Nx = Math.round(p.ncx), Ny = Math.round(p.ncy), a = p.cell, Es = p.E_s, b = p.depth;
  let th = 0, hl = 1; if (p.ctype === 'honeycomb') { th = N.PI / 6; hl = 1; } else if (p.ctype === 'reentrant') { th = -p.th_re * N.PI / 180; hl = p.h_l; if (hl < 2 * Math.sin(-th) + 0.3) { hl = 2 * Math.sin(-th) + 0.3; warnings.push(`h/l raised to ${F(hl, 3)} so that the re-entrant vertices do not overlap (h/l > 2 sin θ).`); } }
  const G = buildLattice(p.ctype, Nx, Ny, a, th, hl), nodes = G.nodes, edges = G.edges, L0 = G.info.l; let t = p.t_strut;
  if (t > 0.3 * L0) { t = 0.3 * L0; warnings.push(`Member thickness limited to ${F(t, 3)} m (30 % of the member length): thicker members are not beams and need a continuum model.`); }
  const A = t * b, I = b * t ** 3 / 12, xs = nodes.map(q => q[0]), ys = nodes.map(q => q[1]), xmin = N.minOf(xs), xmax = N.maxOf(xs), ymin = N.minOf(ys), hex = p.ctype === 'honeycomb' || p.ctype === 'reentrant', ymax = hex ? ymin + Ny * G.info.ch : N.maxOf(ys) /* loaded line = whole number of cell periods (re-entrant vertices protrude above it) */, H = ymax - ymin, Wd = xmax - xmin, tol = 1e-6 * a;
  const bot = [], top = []; nodes.forEach((q, i) => { if (q[1] < ymin + tol) bot.push(i); if (Math.abs(q[1] - ymax) < tol) top.push(i); });
  const span = arr => N.maxOf(arr.map(i => xs[i])) - N.minOf(arr.map(i => xs[i])), trib = arr => arr.length > 1 ? span(arr) * arr.length / (arr.length - 1) : a, Wt = p.ctype === 'square' ? trib(top) : Nx * a; // load-bearing width: Nx + 1 full columns for the square grid, Nx cell widths otherwise (end nodes of hexagonal and triangular rows are half-connected)
  const mid = arr => arr.reduce((bst, i) => Math.abs(xs[i] - (xmin + xmax) / 2) < Math.abs(xs[bst] - (xmin + xmax) / 2) ? i : bst, arr[0]);
  const eps = p.strain / 100, fix = [], load = []; let delta = 0, Fapp = 0;
  if (p.lcase === 'bend3') { const lft = bot.reduce((m, i) => xs[i] < xs[m] ? i : m, bot[0]), rgt = bot.reduce((m, i) => xs[i] > xs[m] ? i : m, bot[0]); fix.push([lft, 0, 0], [lft, 1, 0], [rgt, 1, 0]); Fapp = p.F_load; load.push([mid(top), 1, -Fapp]); }
  else if (p.lcase === 'shear') { delta = eps * H; for (const i of bot) fix.push([i, 0, 0], [i, 1, 0], [i, 2, 0]); for (const i of top) fix.push([i, 0, delta], [i, 1, 0], [i, 2, 0]); }
  else { delta = (p.lcase === 'compression' ? -1 : 1) * eps * H; for (const i of bot) fix.push([i, 1, 0]); fix.push([mid(bot), 0, 0]); for (const i of top) fix.push([i, 1, delta]); }
  const S = frameSolve(nodes, edges, Es, A, I, fix, load);
  if (S.singular) warnings.push('The stiffness matrix was numerically singular (a mechanism or extremely slender members); a small stabilisation was applied and results are unreliable.');
  // member results
  const sig = S.res.map(r => Math.abs(r.N) / A + Math.max(Math.abs(r.M1), Math.abs(r.M2)) * (t / 2) / I), iMax = N.argmax(sig), sMax = sig[iMax], eM = edges[iMax], locX = (xs[eM[0]] + xs[eM[1]]) / 2, locY = (ys[eM[0]] + ys[eM[1]]) / 2;
  const nAbs = N.maxOf(S.res.map(r => Math.abs(r.N))); let buck = Infinity, iB = -1; S.res.forEach((r, i) => { if (r.N < -1e-9 * nAbs) { const m = N.PI ** 2 * Es * I / (r.L * r.L) / -r.N; if (m < buck) { buck = m; iB = i; } } });
  const bendShare = N.mean(S.res.map((r, i) => sig[i] > 0 ? 1 - Math.abs(r.N) / A / sig[i] : 0));
  let Rx = 0, Ry = 0, Rtx = 0, Rty = 0; S.reac.forEach(r => { Rx += r[0]; Ry += r[1]; }); for (const i of top) { Rtx += S.reac[i][0]; Rty += S.reac[i][1]; }
  const eqErr = Math.hypot(Rx, Ry - (p.lcase === 'bend3' ? Fapp : 0)) / (Math.abs(p.lcase === 'bend3' ? Fapp : p.lcase === 'shear' ? Rtx : Rty) || 1);
  // homogenised predictions
  const c = Math.cos(th), s = Math.sin(th), hh = hl * L0, tl = t / L0; let E2h, nu21h, E1h, G12h, relD;
  if (p.ctype === 'square') { E2h = E1h = Es * t / a; nu21h = 0; G12h = Es * tl ** 3 / 2; relD = 2 * t / a; }
  else if (p.ctype === 'triangular') { E2h = E1h = 2 / Math.sqrt(3) * Es * tl; nu21h = 1 / 3; G12h = Math.sqrt(3) / 4 * Es * tl; relD = 2 * Math.sqrt(3) * tl; }
  else { const db = c * L0 ** 3 / t ** 3, da = s * L0 / t, dh = 2 * hh / t; E2h = Es * (hh + L0 * s) / (L0 * c * (db * c + da * s + dh)); nu21h = ((db * s - da * c) / (L0 * c)) / ((db * c + da * s + dh) / (hh + L0 * s)); E1h = Es * L0 * c / ((hh + L0 * s) * (s * s * L0 ** 3 / t ** 3 + c * c * L0 / t)); G12h = Es * tl ** 3 * (hl + s) / (hl * hl * (1 + 2 * hl) * c); relD = tl * (hl + 2) / (2 * c * (hl + s)); }
  // FE effective response
  const midNodes = nodes.map((_, i) => i).filter(i => ys[i] > ymin + 0.25 * H && ys[i] < ymax - 0.25 * H), xmR = midNodes.length ? N.maxOf(midNodes.map(i => xs[i])) : xmax, xmL = midNodes.length ? N.minOf(midNodes.map(i => xs[i])) : xmin, Rg = midNodes.filter(i => xs[i] > xmR - 0.3 * a), Lg = midNodes.filter(i => xs[i] < xmL + 0.3 * a);
  const mU = (arr, cmp) => arr.length ? N.mean(arr.map(i => S.disp[i][cmp])) : 0, mX = arr => arr.length ? N.mean(arr.map(i => xs[i])) : 0;
  let effFE, effHom, effName, effUnit = 'Pa', nuFE = 0, kGlobal, indent = 0, spanDepth = 0, dMax = N.maxOf(S.disp.map(d => Math.hypot(d[0], d[1])));
  if (p.lcase === 'bend3') { const lft = N.minOf(bot.map(i => xs[i])), rgt = N.maxOf(bot.map(i => xs[i])), Ls = rgt - lft, mb = mid(bot), dTop = -S.disp[mid(top)][1], dBot = -S.disp[mb][1], dMid = xs[mb] > lft + 1e-9 * a && xs[mb] < rgt - 1e-9 * a && dBot > 1e-12 * dTop ? dBot : dTop, Ht = p.ctype === 'square' || p.ctype === 'triangular' ? H * (Ny + 1) / Ny : H, Ib = b * Ht ** 3 / 12; indent = dTop - dMid; spanDepth = Ls / H; kGlobal = Fapp / dMid; effFE = kGlobal; effHom = 1 / (Ls ** 3 / (48 * E1h * Ib) + Ls / (4 * (5 / 6) * G12h * b * Ht)); effName = 'Bending stiffness F/δ'; effUnit = 'N/m'; }
  else if (p.lcase === 'shear') { effFE = (Rtx / (Wt * b)) / (delta / H); effHom = G12h / (1 + G12h * H * H / (E2h * Wt * Wt)); effName = 'Apparent shear modulus'; kGlobal = Rtx / delta; }
  else { effFE = (Rty / (Wt * b)) / (delta / H); effHom = E2h; effName = 'Effective modulus along y'; kGlobal = Rty / delta; nuFE = Rg.length && Lg.length ? -((mU(Rg, 0) - mU(Lg, 0)) / (mX(Rg) - mX(Lg))) / (delta / H) : 0; }
  const fsEff = 100 * (effFE / effHom - 1), SF = p.sy / sMax, mass = N.sum(S.res.map(r => r.L)) * A, Ftot = p.lcase === 'bend3' ? Fapp : Math.abs(p.lcase === 'shear' ? Rtx : Rty);
  if (SF < 1) warnings.push(`Peak member stress ${F(sMax / 1e6, 4)} MPa exceeds the yield strength: the linear-elastic result over-estimates load capacity (plastic hinges form at the joints).`);
  if (buck < 1) warnings.push(`Member buckling margin ${F(buck, 3)} is below 1: compressed members buckle before this load is reached.`);
  if (p.lcase === 'bend3' && spanDepth < 4) warnings.push(`Span/depth = ${F(spanDepth, 3)}: the panel is a deep beam, so the comparison with beam theory is indicative only; use at least four times more cells along x than along y for a bending test.`);
  if (Math.min(Nx, Ny) < 5) warnings.push('Fewer than five cells across: the response is dominated by edge effects and differs from the homogenised material.');
  if (p.ctype === 'square' && p.lcase !== 'shear' && p.lcase !== 'bend3') warnings.push('The square grid loaded along its members has no lateral coupling (ν ≈ 0) and is extremely shear-compliant.');
  const scale = dMax > 0 ? 0.06 * Math.max(Wd, H) / dMax : 0, mm = 1e3, gNodes = nodes.map(q => [q[0] * mm, q[1] * mm, 0]), frames = [0, 0.2, 0.4, 0.6, 0.8, 1].map(f => ({ t: f, nodes: nodes.map((q, i) => [(q[0] + f * scale * S.disp[i][0]) * mm, (q[1] + f * scale * S.disp[i][1]) * mm, 0]) }));
  const sorted = [...sig].sort((x, y) => y - x), ranks = sorted.map((_, i) => (i + 1) / sorted.length), order = sig.map((_, i) => i).sort((x, y) => sig[y] - sig[x]).slice(0, 6);
  const prof = p.lcase === 'bend3' ? bot.slice().sort((x, y) => xs[x] - xs[y]) : nodes.map((_, i) => i).filter(i => xs[i] > xmR - 0.3 * a).sort((x, y) => ys[x] - ys[y]);
  const aux = nuFE < -0.02 && (p.lcase === 'tension' || p.lcase === 'compression');
  return {
    summary: `${Nx}×${Ny}-cell ${p.ctype} lattice panel (${nodes.length} nodes, ${edges.length} members, ${S.dof} DOF) under ${p.lcase === 'bend3' ? `three-point bending with ${F(Fapp)} N` : `${p.lcase} at ${F(p.strain)} % nominal strain`}: ${effName.toLowerCase()} ${F(effFE, 4)} ${effUnit} against ${F(effHom, 4)} ${effUnit} from homogenisation (finite-size effect ${F(fsEff, 3)} %)${p.lcase === 'tension' || p.lcase === 'compression' ? `, Poisson's ratio ${F(nuFE, 3)} (homogenised ${F(nu21h, 3)})${aux ? ' — auxetic' : ''}` : ''}. Peak member stress ${F(sMax / 1e6, 4)} MPa at (${F(locX * 1e3, 3)}, ${F(locY * 1e3, 3)}) mm gives a yield safety factor of ${F(SF, 3)}; buckling margin ${Number.isFinite(buck) ? F(buck, 3) : 'not applicable (no member in compression)'}.`,
    kpis: [
      kpi('eff_fe', effName + ' (finite panel)', effFE, effUnit, `From the summed support reactions of the frame model divided by the tributary edge area and the applied strain${p.lcase === 'bend3' ? ' (here: applied force over mid-span deflection)' : ''}; includes clamped-edge and free-edge effects.`, null),
      kpi('eff_hom', effName + ' (homogenised)', effHom, effUnit, p.lcase === 'bend3' ? `Timoshenko beam with the homogenised modulus E₁ and shear modulus G₁₂ of the infinite lattice (shear factor 5/6, tributary section depth), span/depth = ${F(spanDepth, 3)}. Deflection is taken at the lower chord; local indentation under the load adds ${F(indent, 3)} m.` : p.lcase === 'shear' ? `Homogenised G₁₂ = ${F(G12h, 4)} Pa of the infinite lattice, reduced by the bending compliance of the clamped panel: G/(1 + G·H²/(E·W²)).` : 'Infinite-lattice prediction from unit-cell beam theory (cell-wall bending plus axial deformation).', null),
      kpi('fs_effect', 'Finite-size / edge effect', fsEff, '%', 'Relative difference between the finite panel and the homogenised prediction. It decays as more cells are added; clamped edges stiffen and free edges soften the panel.', Math.abs(fsEff) < 15),
      kpi('nu_fe', 'Poisson\'s ratio of the panel', nuFE, '–', p.lcase === 'tension' || p.lcase === 'compression' ? `Negative of lateral over axial strain measured on the free edges in the middle half of the panel; the unit-cell value is ${F(nu21h, 3)}. Negative values mean the panel widens when stretched.` : 'Not evaluated for this load case (reported as 0); run tension or compression to measure it.', null),
      kpi('s_max', 'Maximum member stress', sMax, 'Pa', `Largest combined axial-plus-bending fibre stress |N|/A + |M|(t/2)/I over all member ends, located at x = ${F(locX * 1e3, 3)} mm, y = ${F(locY * 1e3, 3)} mm.`, SF >= 1.5),
      kpi('SF', 'Yield safety factor', SF, '–', 'Yield strength divided by the peak member stress. Because the analysis is linear, the load at first yield is the applied load times this factor.', SF >= 1.5),
      kpi('buck', 'Member buckling margin', Number.isFinite(buck) ? buck : 1e6, '–', Number.isFinite(buck) ? 'Smallest ratio of the pinned-end Euler load π²EI/L² to the compressive force in any member (conservative: rigid joints give some end restraint).' : 'No member is in compression, so member buckling does not arise; reported as a large number.', buck >= 1.5),
      kpi('k_glob', 'Global stiffness', kGlobal, 'N/m', 'Total reaction (or applied) force divided by the imposed (or resulting) displacement of the loaded edge or load point.', null),
      kpi('F_tot', 'Total load carried', Ftot, 'N', 'Resultant force on the loaded edge at the applied deformation (or the applied force in bending).', null),
      kpi('bend_share', 'Bending share of member stress', bendShare, '–', 'Mean fraction of the peak fibre stress that comes from bending rather than axial force: near 1 for bending-dominated cells, near 0 for stretch-dominated ones.', null),
      kpi('mass', 'Panel mass per unit density', mass, 'm³', 'Solid volume of all members; multiply by the material density to obtain the mass. Relative density of the cell is ' + F(relD, 3) + '.', null),
      kpi('eq_err', 'Global equilibrium error', eqErr, '–', 'Resultant of all support reactions plus applied loads, relative to the load carried; it should be at round-off level.', eqErr < 1e-6),
    ],
    series: [
      { key: 'stress_rank', title: 'Member stress distribution', xlabel: 'fraction of members (ranked)', ylabel: 'peak fibre stress (MPa)', lines: [line('members', ranks, sorted.map(v => v / 1e6)), line('yield strength', [0, 1], [p.sy / 1e6, p.sy / 1e6])], explain: 'Members ranked from most to least stressed. A steep initial drop indicates stress concentration at edges or load points; a flat curve indicates uniform utilisation.' },
      p.lcase === 'bend3' ? { key: 'defl', title: 'Deflection of the lower chord', xlabel: 'x (mm)', ylabel: 'vertical displacement (µm)', lines: [line('bottom nodes', prof.map(i => xs[i] * 1e3), prof.map(i => S.disp[i][1] * 1e6))], explain: 'Vertical displacement of the bottom boundary nodes between the supports; it should resemble a simply supported beam with additional shear deflection.' }
        : { key: 'edge', title: 'Lateral displacement of the right free edge', xlabel: 'y (mm)', ylabel: 'x-displacement (µm)', lines: [line('edge nodes', prof.map(i => ys[i] * 1e3), prof.map(i => S.disp[i][0] * 1e6))], explain: 'Horizontal displacement of the right-edge nodes along the panel height. A uniform plateau away from the supports shows the bulk Poisson response; the ends show the grip constraint.' },
    ],
    tables: [{ key: 'critical', title: 'Most highly stressed members', columns: ['x (mm)', 'y (mm)', 'Axial force (N)', 'Max end moment (N·m)', 'Stress (MPa)', 'Utilisation'], rows: order.map(i => [F((xs[edges[i][0]] + xs[edges[i][1]]) / 2 * 1e3, 4), F((ys[edges[i][0]] + ys[edges[i][1]]) / 2 * 1e3, 4), F(S.res[i].N, 4), F(Math.max(Math.abs(S.res[i].M1), Math.abs(S.res[i].M2)), 4), F(sig[i] / 1e6, 4), F(sig[i] / p.sy, 3)]), explain: 'Location and internal forces of the six most highly stressed members: the probable failure initiation sites of the panel.' }],
    equations: [
      eqn('governing', 'Direct stiffness equilibrium', 'K u = f,  K = Σₑ TₑᵀkₑTₑ  (Euler–Bernoulli frame elements, 3 DOF per node)'),
      eqn('constitutive', 'Member forces and stress', 'N = EA(u₂′ − u₁′)/L;  M = EI κ;  σ = |N|/A + |M|(t/2)/I'),
      eqn('formula', 'Homogenised cell (bending + axial)', 'E₂ = E_s(h + l sin θ) / {l cos θ [l³cos²θ/t³ + (l sin²θ + 2h)/t]};  ν₂₁ from the same wall deflections'),
      eqn('formula', 'Member buckling', 'P_cr = π²EI/L²  (pinned ends, conservative)'),
    ],
    geometry: { kind: 'lattice', nodes: gNodes, edges, radius: t * mm / 2, units: 'mm', frames, note: `deformation magnified ${F(scale, 3)}×` },
    checks: [
      chk('Global force equilibrium', eqErr < 1e-6, `|ΣR + ΣF| / load = ${eqErr.toExponential(2)}`),
      chk('Solver residual at free DOF', S.residual < 1e-6, `relative residual ${S.residual.toExponential(2)}, bandwidth ${S.bw}`),
      chk('Stiffness matrix positive definite', !S.singular, S.singular ? 'pivot below tolerance' : 'all Cholesky pivots positive'),
      chk('Members are slender beams', t <= 0.2 * L0, `t/l = ${F(t / L0, 3)}`),
      chk('Stress below yield', SF >= 1, `safety factor ${F(SF, 3)}`),
      chk('No member buckling', buck >= 1, Number.isFinite(buck) ? `margin ${F(buck, 3)}` : 'no compression members'),
      chk('Small-displacement assumption', dMax < 0.05 * Math.max(Wd, H), `max displacement ${F(dMax, 3)} m vs panel ${F(Math.max(Wd, H), 3)} m`),
    ],
    warnings,
    recommendations: [
      rec('high', SF < 1.5 ? `Peak stress is ${F(100 / SF, 3)} % of yield. Round the joints (fillet radius about one wall thickness) or thicken the walls near (${F(locX * 1e3, 3)}, ${F(locY * 1e3, 3)}) mm, where failure will start.` : `Safety factor ${F(SF, 3)} leaves margin; the wall thickness could be reduced by up to ${F(100 * (1 - Math.pow(1.5 / SF, bendShare > 0.5 ? 0.5 : 1)), 3)} % (keeping a factor of 1.5) to save mass.`, 'Peak fibre stress from frame analysis.'),
      rec('high', Math.abs(fsEff) > 15 ? `The finite panel differs from the homogenised model by ${F(fsEff, 3)} %: do not use effective properties for this size; either add cells or analyse the discrete lattice.` : 'Edge effects are below 15 %: homogenised properties are adequate for component-level design of this panel.', 'Comparison of discrete and homogenised stiffness.'),
      rec('medium', Number.isFinite(buck) && buck < 3 ? `Compressed members are within a factor ${F(buck, 3)} of Euler buckling; shorten or thicken them, since elastic buckling triggers sudden collapse of the cell row.` : 'Member buckling is not critical for this load case; check it again under compression or reversed loading.', 'Euler load of individual members.'),
      rec('medium', aux ? 'Exploit the auxetic response where conformability, indentation resistance or shear stiffness matters; verify the measured ν with digital image correlation on a panel of the same cell count.' : 'Measure the full-field strain on a printed coupon with the same number of cells to validate joint stiffness assumptions (rigid, fillet-free joints here).', 'Model assumes ideal rigid joints and perfect geometry.'),
      rec('low', `The lattice uses ${F(100 * relD, 3)} % of the solid volume; additive manufacture of only the load-bearing members avoids machining waste and reduces embodied energy roughly in proportion to mass.`, 'Relative density of the cell.', true),
    ],
    publish: { 'mech.max_stress': sMax, 'mech.safety_factor': SF, 'mech.buckling_load': Number.isFinite(buck) ? buck * Ftot : Ftot * 1e6, 'mech.stiffness': kGlobal, 'geom.rel_density': relD, ...(p.lcase === 'tension' || p.lcase === 'compression' ? { 'eff.E': effFE, 'eff.nu': nuFE } : p.lcase === 'shear' ? { 'eff.G': effFE } : {}) },
    numerics: { method: '2-D Euler–Bernoulli frame finite elements, direct stiffness method, banded Cholesky factorisation', dof: S.dof, steps: 1, iterations: 1, residual: S.residual, notes: `Half bandwidth ${S.bw}; linear elastic, small displacements, rigid joints.` },
  };
}

// ---------- homogenised Kirchhoff plate ----------
function plate(p) {
  const a = p.pa, b = p.pb, h = p.ph, E = p.E_p, nu = p.nu_p, D = E * h ** 3 / (12 * (1 - nu * nu)), q = p.q_p, warnings = [];
  const term = (m, n) => 16 * q / (N.PI ** 6 * D * m * n * (m * m / (a * a) + n * n / (b * b)) ** 2);
  const centre = M => { let w = 0, mx = 0, my = 0; for (let m = 1; m <= M; m += 2) for (let n = 1; n <= M; n += 2) { const W = term(m, n), sg = Math.sin(m * N.PI / 2) * Math.sin(n * N.PI / 2); w += W * sg; mx += D * W * N.PI ** 2 * (m * m / (a * a) + nu * n * n / (b * b)) * sg; my += D * W * N.PI ** 2 * (n * n / (b * b) + nu * m * m / (a * a)) * sg; } return { w, mx, my }; };
  const c1 = centre(39), c0 = centre(19), wMax = c1.w, Mmax = Math.max(Math.abs(c1.mx), Math.abs(c1.my)), sMax = 6 * Mmax / (h * h), conv = Math.max(Math.abs(c1.w - c0.w) / Math.abs(c1.w), Math.abs(c1.mx - c0.mx) / Math.abs(c1.mx));
  const fr = []; for (let m = 1; m <= 6; m++) for (let n = 1; n <= 6; n++) fr.push({ m, n, f: N.PI / 2 * (m * m / (a * a) + n * n / (b * b)) * Math.sqrt(D / (p.rho_p * h)) }); fr.sort((x, y) => x.f - y.f);
  let Ncr = Infinity, mcr = 1; for (let m = 1; m <= 12; m++) { const v = N.PI ** 2 * D / (b * b) * (m * b / a + a / (m * b)) ** 2; if (v < Ncr) { Ncr = v; mcr = m; } }
  const mass = p.rho_p * a * b * h, SF = p.sy / sMax, bSF = p.Nx_p > 0 ? Ncr / p.Nx_p : 1e6, coef = wMax * D / (q * Math.min(a, b) ** 4), square = Math.abs(a - b) < 1e-9 * a, thin = Math.min(a, b) / h;
  // curvature character under pure bending: κ_y = −ν κ_x
  const synclastic = nu < 0, nx = 41, ny = 41, map = new Float64Array(nx * ny), sx = [], sy = [];
  for (let m = 1; m <= 19; m += 2) { sx.push(N.linspace(0, 1, nx).map(x => Math.sin(m * N.PI * x))); sy.push(N.linspace(0, 1, ny).map(y => Math.sin(m * N.PI * y))); }
  for (let mi = 0; mi < 10; mi++) for (let ni = 0; ni < 10; ni++) { const W = term(2 * mi + 1, 2 * ni + 1); for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) map[i + nx * j] += W * sx[mi][i] * sy[ni][j] * 1e3; }
  const Dref = E * h ** 3 / (12 * (1 - 0.09));
  if (wMax > 0.5 * h) warnings.push(`Centre deflection is ${F(wMax / h, 3)} plate thicknesses: beyond about 0.5 h membrane stretching stiffens the plate and linear Kirchhoff theory over-predicts deflection.`);
  if (thin < 20) warnings.push(`Span/thickness = ${F(thin, 3)} is below 20: transverse shear (Mindlin) effects are significant, particularly for low-shear-modulus lattice cores.`);
  if (nu < -0.8) warnings.push('Poisson\'s ratio is close to the −1 stability limit, where the bending stiffness diverges and small errors in ν change the result strongly.');
  const xs = N.linspace(0, a, 81), wline = xs.map(x => { let w = 0; for (let m = 1; m <= 29; m += 2) for (let n = 1; n <= 29; n += 2) w += term(m, n) * Math.sin(m * N.PI * x / a) * Math.sin(n * N.PI / 2); return w * 1e3; });
  const nus = N.linspace(-0.9, 0.45, 55);
  return {
    summary: `Simply supported ${F(a)} × ${F(b)} m plate, ${F(h * 1e3)} mm thick, with effective E = ${F(E / 1e9, 4)} GPa and ν = ${F(nu, 3)}${synclastic ? ' (auxetic: synclastic, dome-like bending)' : ''}: flexural rigidity D = ${F(D, 4)} N·m, centre deflection ${F(wMax * 1e3, 4)} mm under ${F(q)} Pa, peak bending stress ${F(sMax / 1e6, 4)} MPa (safety factor ${F(SF, 3)}). First natural frequency ${F(fr[0].f, 4)} Hz, buckling load ${F(Ncr, 4)} N/m (m = ${mcr}), mass ${F(mass, 4)} kg.`,
    kpis: [
      kpi('D', 'Flexural rigidity', D, 'N·m', `D = Eh³/(12(1 − ν²)). With ν = ${F(nu, 3)} it is ${F(D / Dref, 4)} times that of a plate with ν = 0.3 and the same modulus: |ν| near 1 stiffens the plate in bending.`, null),
      kpi('w_max', 'Centre deflection', wMax, 'm', `Navier double sine series (odd terms to 39) for uniform pressure; equals ${F(coef, 4)}·q·b_min⁴/D${square ? ' (classical value for a square plate: 0.00406)' : ''}.`, wMax <= 0.5 * h),
      kpi('s_max', 'Maximum bending stress', sMax, 'Pa', 'Surface stress 6M/h² from the larger of the two centre bending moments M_x and M_y of the homogenised plate. Local member stresses in a lattice core are higher.', SF >= 1.5),
      kpi('SF', 'Stress safety factor', SF, '–', 'Allowable (yield) stress divided by the maximum homogenised bending stress.', SF >= 1.5),
      kpi('f1', 'Fundamental frequency', fr[0].f, 'Hz', 'f₁₁ = (π/2)(1/a² + 1/b²)√(D/ρh) for simply supported edges.', null),
      kpi('f2', 'Second natural frequency', fr[1].f, 'Hz', `Mode (${fr[1].m}, ${fr[1].n}); the spacing to the fundamental indicates how isolated the first resonance is.`, null),
      kpi('N_cr', 'Buckling load (uniaxial)', Ncr, 'N/m', `N_cr = (π²D/b²)·min_m(mb/a + a/mb)² with ${mcr} half-wave(s) along the loading direction.`, bSF >= 1.5),
      kpi('buck_SF', 'Buckling safety factor', bSF, '–', p.Nx_p > 0 ? 'Critical edge load divided by the applied in-plane compression.' : 'No in-plane compression applied; reported as a large number.', bSF >= 1.5),
      kpi('mass', 'Panel mass', mass, 'kg', 'Effective (homogenised) density multiplied by the plate volume a·b·h.', null),
      kpi('w_h', 'Deflection / thickness', wMax / h, '–', 'Linear plate theory is accurate while this ratio stays below about 0.5.', wMax / h <= 0.5),
      kpi('anticl', 'Curvature coupling κ_y/κ_x in pure bending', -nu, '–', 'Under a pure bending moment the transverse curvature is −ν times the principal one: negative for conventional plates (saddle, anticlastic), positive for auxetic plates (dome, synclastic).', null),
    ],
    series: [
      { key: 'w_x', title: 'Deflection along the centre line', xlabel: 'x (m)', ylabel: 'deflection (mm)', lines: [line('w(x, b/2)', xs, wline)], explain: 'Deflected shape along the mid-width line under uniform pressure, from the Navier series.' },
      { key: 'D_nu', title: 'Flexural rigidity versus Poisson\'s ratio', xlabel: 'Poisson\'s ratio', ylabel: 'D / D(ν = 0)', lines: [line('1/(1 − ν²)', nus, nus.map(v => 1 / (1 - v * v)))], explain: 'For a fixed Young\'s modulus the plate rigidity grows as |ν| approaches its limits; strongly auxetic panels are stiffer in bending than conventional ones of the same modulus.' },
    ],
    maps: [{ key: 'w', title: 'Deflection field', nx, ny, data: map, x0: 0, x1: a, y0: 0, y1: b, xlabel: 'x (m)', ylabel: 'y (m)', unit: 'mm', explain: 'Transverse deflection of the plate under uniform pressure; the maximum is at the centre and the edges are simply supported.' }],
    tables: [{ key: 'modes', title: 'Natural frequencies', columns: ['Mode (m, n)', 'Frequency (Hz)', 'Ratio to fundamental'], rows: fr.slice(0, 8).map(r => [`(${r.m}, ${r.n})`, F(r.f, 5), F(r.f / fr[0].f, 4)]), explain: 'First eight flexural modes of the simply supported plate with half-wave numbers m (along a) and n (along b).' }],
    equations: [
      eqn('governing', 'Kirchhoff plate equation', 'D∇⁴w = q,  D = Eh³/[12(1 − ν²)]'),
      eqn('formula', 'Navier solution', 'w = Σ_{m,n odd} 16q / [π⁶ D m n (m²/a² + n²/b²)²] · sin(mπx/a) sin(nπy/b)'),
      eqn('formula', 'Natural frequencies', 'f_mn = (π/2)(m²/a² + n²/b²)√(D/ρh)'),
      eqn('formula', 'Uniaxial buckling', 'N_cr = (π²D/b²)(mb/a + a/(mb))²'),
      eqn('constitutive', 'Moment–curvature', 'M_x = −D(w,ₓₓ + ν w,ᵧᵧ);  pure bending: κ_y = −ν κ_x'),
    ],
    checks: [
      chk('Thermodynamic admissibility of ν', nu > -1 && nu < 0.5, `ν = ${F(nu, 3)} within (−1, 0.5)`),
      chk('Navier series converged', conv < 1e-3, `change from 19 to 39 terms: ${conv.toExponential(2)}`),
      ...(square ? [chk('Square-plate benchmark coefficient', Math.abs(coef - 0.00406) < 5e-5, `w·D/(q·a⁴) = ${F(coef, 4)} vs 0.00406`)] : []),
      chk('Small-deflection theory valid', wMax <= 0.5 * h, `w/h = ${F(wMax / h, 3)}`),
      chk('Thin-plate assumption', thin >= 20, `span/thickness = ${F(thin, 3)}`),
      chk('Stress below allowable', SF >= 1, `safety factor ${F(SF, 3)}`),
      chk('Stable against buckling', bSF >= 1, `buckling factor ${F(bSF, 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', SF < 1.5 || wMax > 0.5 * h ? `Increase thickness: deflection scales with h⁻³ and stress with h⁻²; a thickness of ${F(h * Math.max(Math.cbrt(wMax / (0.5 * h)), Math.sqrt(1.5 / SF), 1) * 1e3, 3)} mm would meet both the small-deflection limit and a stress factor of 1.5.` : `The panel meets stress (factor ${F(SF, 3)}) and deflection limits; thickness could be reduced to about ${F(h * Math.max(Math.sqrt(1.5 / SF), Math.cbrt(wMax / (0.5 * h))) * 1e3, 3)} mm before either governs.`, 'Kirchhoff scaling of deflection and stress with thickness.'),
      rec('medium', synclastic ? 'The auxetic panel bends into a dome: use it for doubly curved shells and conformal panels without the forming stresses a conventional sheet would need.' : 'For doubly curved (dome) shapes consider an auxetic core, which bends synclastically instead of forming a saddle.', 'Curvature coupling κ_y = −ν κ_x.'),
      rec('medium', `Keep excitation away from ${F(fr[0].f, 3)} Hz or add damping; verify boundary conditions, because clamped edges raise the fundamental frequency by roughly a factor 1.8 for a square plate.`, 'First-mode frequency of the simply supported plate.'),
      rec('low', `At ${F(mass, 3)} kg the homogenised panel is ${F(p.rho_p, 3)} kg/m³ dense; choosing a core architecture with higher specific bending stiffness reduces material use and transport energy for the same deflection.`, 'Mass of the panel.', true),
    ],
    publish: { 'mech.max_stress': sMax, 'mech.safety_factor': Math.min(SF, bSF), 'mech.buckling_load': Ncr * b, 'res.f0': fr[0].f, 'mech.deflection': wMax },
    numerics: { method: 'Navier double Fourier series for the Kirchhoff plate (closed-form modal and buckling expressions)', dof: 400, steps: 1, iterations: 400, residual: conv, notes: 'Simply supported edges, uniform pressure, isotropic effective properties.' },
  };
}

// ---------- resonator chain ----------
function modal(p) {
  const nc = Math.round(p.n_cells), n = 2 * nc, k = p.k_s, eta = p.eta, warnings = [], ms = Array.from({ length: n }, (_, i) => i % 2 === 0 ? p.m1 : p.m2);
  const K = N.zeros(n, n), M = N.zeros(n, n); for (let i = 0; i < n; i++) { M[i][i] = ms[i]; K[i][i] = i < n - 1 ? 2 * k : k; if (i > 0) { K[i][i - 1] = -k; K[i - 1][i] = -k; } }
  const eg = N.eigGen(K, M), fn = eg.values.map(v => Math.sqrt(Math.max(v, 0)) / (2 * N.PI));
  const mHi = Math.max(p.m1, p.m2), mLo = Math.min(p.m1, p.m2), fA = Math.sqrt(2 * k / mHi) / (2 * N.PI), fO = Math.sqrt(2 * k / mLo) / (2 * N.PI), fC = Math.sqrt(2 * k * (1 / p.m1 + 1 / p.m2)) / (2 * N.PI), hasGap = mHi / mLo > 1.0001;
  const ext = p.gap_lo > 0 && p.gap_hi > p.gap_lo, gLo = ext ? p.gap_lo : fA, gHi = ext ? p.gap_hi : fO;
  const trans = f => { // complex Thomas for (K(1+iη) − ω²M) x = f, base displacement 1 through the first spring
    const w2 = (2 * N.PI * f) ** 2, br = new Float64Array(n), bi = new Float64Array(n), dr = new Float64Array(n), di = new Float64Array(n), or = -k, oi = -k * eta;
    for (let i = 0; i < n; i++) { const kd = i < n - 1 ? 2 * k : k; br[i] = kd - w2 * ms[i]; bi[i] = kd * eta; } dr[0] = k; di[0] = k * eta;
    for (let i = 1; i < n; i++) { const den = br[i - 1] ** 2 + bi[i - 1] ** 2, mr = (or * br[i - 1] + oi * bi[i - 1]) / den, mi = (oi * br[i - 1] - or * bi[i - 1]) / den; br[i] -= mr * or - mi * oi; bi[i] -= mr * oi + mi * or; dr[i] -= mr * dr[i - 1] - mi * di[i - 1]; di[i] -= mr * di[i - 1] + mi * dr[i - 1]; }
    const den = br[n - 1] ** 2 + bi[n - 1] ** 2; return Math.hypot((dr[n - 1] * br[n - 1] + di[n - 1] * bi[n - 1]) / den, (di[n - 1] * br[n - 1] - dr[n - 1] * bi[n - 1]) / den);
  };
  const fs = N.linspace(0.02 * fC, 1.4 * fC, 360), T = fs.map(trans);
  const kappa = f => { const w2 = (2 * N.PI * f) ** 2, cq = ((2 * k - p.m1 * w2) * (2 * k - p.m2 * w2) - 2 * k * k) / (2 * k * k); return Math.abs(cq) > 1 ? Math.acosh(Math.abs(cq)) : 0; };
  const fMid = hasGap ? Math.sqrt(fA * fO) : 1.2 * fC, Tmid = trans(fMid), kap = kappa(fMid), Tpred = Math.exp(-nc * kap), inGapIdx = fs.map((f, i) => f > gLo && f < gHi ? i : -1).filter(i => i >= 0), Tmin = inGapIdx.length ? N.minOf(inGapIdx.map(i => T[i])) : Tmid;
  const modesInGap = hasGap ? fn.filter(f => f > fA * 1.001 && f < fO * 0.999).length : 0, Tpeak = N.maxOf(T), overlap = ext && hasGap ? Math.max(0, Math.min(gHi, fO) - Math.max(gLo, fA)) / (gHi - gLo) : 1;
  const needCells = kap > 0 ? Math.ceil(Math.log(100) / kap) : 0; // cells for 40 dB
  if (!hasGap) warnings.push('Equal masses: the chain is monatomic and has no Bragg band gap, only a cut-off frequency. Introduce a mass (or stiffness) contrast to open a gap.');
  if (ext && hasGap && overlap < 0.5) warnings.push(`The supplied Bloch gap (${F(gLo, 4)}–${F(gHi, 4)} Hz) overlaps the analytical gap of this chain (${F(fA, 4)}–${F(fO, 4)} Hz) by only ${F(100 * overlap, 3)} %: the chain parameters do not represent that unit cell.`);
  if (modesInGap > 0) warnings.push(`${modesInGap} natural frequency(ies) of the finite chain lie inside the Bloch gap: these are boundary (edge) modes localised at the free end, absent from the infinite-lattice analysis.`);
  const v1 = eg.vectors[0], vmax = N.maxOf(v1.map(Math.abs)) || 1, nodes = [[-1, 0, 0], ...ms.map((_, i) => [i, 0, 0])], edges = ms.map((_, i) => [i, i + 1]), frames = N.linspace(0, 1, 9).map(t => ({ t, nodes: [[-1, 0, 0], ...ms.map((_, i) => [i + 0.35 * v1[i] / vmax * Math.sin(2 * N.PI * t), 0, 0])] }));
  return {
    summary: `Finite chain of ${nc} two-mass cells (${n} DOF): natural frequencies from ${F(fn[0], 4)} to ${F(fn[n - 1], 4)} Hz. ${hasGap ? `The infinite-lattice Bloch gap is ${F(fA, 4)}–${F(fO, 4)} Hz (relative width ${F(100 * (fO - fA) / fMid, 3)} %); the finite chain transmits ${F(2 * dB(Tmid), 4)} dB at mid-gap against ${F(2 * dB(Tpred), 4)} dB from the Bloch decay rate (${F(8.686 * kap, 3)} dB per cell).` : 'With equal masses there is no band gap; transmission falls only above the cut-off at ' + F(fC, 4) + ' Hz.'} Peak resonant amplification is ${F(Tpeak, 3)} for a loss factor of ${F(eta)}.`,
    kpis: [
      kpi('f1', 'Fundamental frequency', fn[0], 'Hz', 'Lowest eigenfrequency of the fixed–free chain from the generalised eigenproblem K v = ω² M v.', null),
      kpi('f_max', 'Highest natural frequency', fn[n - 1], 'Hz', `Top of the optical branch; the infinite-chain cut-off is ${F(fC, 4)} Hz.`, null),
      kpi('gap_lo', 'Band-gap lower edge', hasGap ? fA : 0, 'Hz', 'Top of the acoustic branch of the infinite diatomic chain, √(2k/m_heavy)/2π. Zero when no gap exists.', hasGap),
      kpi('gap_hi', 'Band-gap upper edge', hasGap ? fO : 0, 'Hz', 'Bottom of the optical branch, √(2k/m_light)/2π. Zero when no gap exists.', hasGap),
      kpi('T_mid', 'Mid-gap transmission (finite chain)', 2 * dB(Tmid), 'dB', `Displacement transmissibility |x_N/x₀| of the damped finite chain at ${F(fMid, 4)} Hz, from a complex tridiagonal solve.`, 2 * dB(Tmid) < -20),
      kpi('T_bloch', 'Mid-gap transmission (Bloch estimate)', 2 * dB(Tpred), 'dB', 'exp(−N·κ) with the evanescent Bloch wavenumber κ = acosh|cos qa| per cell: the infinite-lattice prediction for N cells, without boundary reflections or damping.', null),
      kpi('att_cell', 'Attenuation per cell at mid-gap', 8.686 * kap, 'dB', 'Spatial decay of the evanescent wave per unit cell; multiply by the number of cells for the total isolation.', kap > 0),
      kpi('cells_40', 'Cells for 40 dB isolation', needCells, '–', 'Number of unit cells needed for 40 dB amplitude attenuation at mid-gap according to the Bloch decay rate (0 if there is no gap).', null),
      kpi('T_min', 'Minimum transmission inside the gap', 2 * dB(Tmin), 'dB', `Lowest transmissibility found between ${F(gLo, 4)} and ${F(gHi, 4)} Hz in the frequency sweep.`, null),
      kpi('edge_modes', 'Finite-chain modes inside the gap', modesInGap, '–', 'Number of natural frequencies of the finite chain that fall inside the infinite-lattice gap: localised boundary modes that can leak energy.', modesInGap === 0),
      kpi('Q', 'Modal quality factor', 1 / eta, '–', 'For hysteretic damping every mode has Q ≈ 1/η, which sets the resonant amplification in the pass bands.', null),
    ],
    series: [
      { key: 'T_f', title: 'Transmission of the finite chain', xlabel: 'frequency (Hz)', ylabel: 'transmissibility |x_N / x₀|', logy: true, lines: [line(`${nc} cells`, fs, T.map(v => Math.max(v, 1e-30))), ...(hasGap || ext ? [line('gap lower edge', [gLo, gLo], [1e-6, 10]), line('gap upper edge', [gHi, gHi], [1e-6, 10])] : [])], explain: 'Harmonic response of the last mass to base excitation. Pass bands show resonance peaks; inside the Bloch gap the transmission drops exponentially with the number of cells.' },
      { key: 'modes', title: 'Natural frequencies', xlabel: 'mode number', ylabel: 'frequency (Hz)', lines: [line('finite chain', fn.map((_, i) => i + 1), fn)], explain: 'Sorted eigenfrequencies. The jump between mode N and N + 1 is the finite-structure image of the band gap; isolated points inside the jump are edge modes.' },
      { key: 'att', title: 'Bloch attenuation per cell', xlabel: 'frequency (Hz)', ylabel: 'attenuation (dB per cell)', lines: [line('8.686·κ', fs, fs.map(f => 8.686 * kappa(f)))], explain: 'Evanescent decay rate of the infinite lattice: zero in pass bands, positive in the gap and above the cut-off.' },
    ],
    tables: [{ key: 'mlist', title: 'Lowest modes', columns: ['Mode', 'Frequency (Hz)', 'Branch'], rows: fn.slice(0, Math.min(10, n)).map((f, i) => [i + 1, F(f, 5), !hasGap ? 'pass band' : f <= fA * 1.001 ? 'acoustic' : f < fO * 0.999 ? 'gap (edge mode)' : 'optical']), explain: 'First ten natural frequencies with the branch of the infinite-lattice dispersion relation to which each belongs.' }],
    equations: [
      eqn('governing', 'Modal eigenproblem', 'K v = ω² M v'),
      eqn('governing', 'Harmonic response', '[K(1 + iη) − ω²M] x = f(ω)'),
      eqn('formula', 'Diatomic Bloch dispersion', 'cos(qa) = [(2k − m₁ω²)(2k − m₂ω²) − 2k²] / (2k²)'),
      eqn('formula', 'Band-gap edges', 'ω_lo = √(2k/m_heavy),  ω_hi = √(2k/m_light),  ω_cut = √(2k(1/m₁ + 1/m₂))'),
    ],
    geometry: { kind: 'lattice', nodes, edges, radius: 0.08, units: 'cell', frames },
    checks: [
      chk('All eigenvalues positive', eg.values[0] > 0, `lowest ω² = ${F(eg.values[0], 4)} rad²/s²`),
      chk('Spectrum bounded by the Bloch cut-off', fn[n - 1] <= fC * 1.0001, `${F(fn[n - 1], 5)} ≤ ${F(fC, 5)} Hz`),
      chk('Finite-chain attenuation consistent with Bloch decay', !hasGap || Math.abs(2 * dB(Tmid) - 2 * dB(Tpred)) < 12 + 0.25 * Math.abs(2 * dB(Tpred)), `${F(2 * dB(Tmid), 4)} dB vs ${F(2 * dB(Tpred), 4)} dB`),
      chk('Quasi-static transmissibility tends to 1', Math.abs(trans(0.05 * fn[0]) - 1) < 0.02, `|x_N/x₀| = ${F(trans(0.05 * fn[0]), 5)} at ${F(0.05 * fn[0], 3)} Hz (5 % of the fundamental)`),
      chk('Supplied Bloch gap consistent with this chain', overlap >= 0.5, ext ? `overlap ${F(100 * overlap, 3)} %` : 'no external gap supplied; analytical edges used'),
    ],
    warnings,
    recommendations: [
      rec('high', hasGap ? `Use at least ${needCells} cells for 40 dB isolation at mid-gap; the present ${nc} cells give ${F(-2 * dB(Tmid), 3)} dB.` : 'Introduce a mass ratio of at least 2 between alternating resonators to open a usable band gap.', 'Exponential decay of evanescent Bloch waves.'),
      rec('medium', modesInGap > 0 ? 'Edge modes lie inside the gap: terminate the chain on the heavy mass or add a tuned absorber at the free end to remove them.' : 'No edge modes pollute the gap for this termination; keep the same termination in the built structure.', `${modesInGap} mode(s) inside the gap.`),
      rec('medium', `Measure the transmission of a short prototype (3–4 cells) first: the gap edges ${F(fA, 3)}–${F(fO, 3)} Hz and the per-cell attenuation of ${F(8.686 * kap, 3)} dB are already visible and validate k and the masses before the full array is built.`, 'Finite-array convergence to Bloch behaviour.'),
      rec('low', 'Passive periodic isolation needs no power or control hardware; prefer it to active cancellation where the target band is fixed, and choose the lightest mass contrast that gives the required gap.', 'Passive band-gap isolation.', true),
    ],
    publish: { 'res.f0': fn[0], 'res.bandgap_lo': hasGap ? fA : 0, 'res.bandgap_hi': hasGap ? fO : 0, 'res.transmission': Tmid * Tmid, 'res.Q': 1 / eta },
    numerics: { method: 'Generalised symmetric eigenproblem (Cholesky reduction + Jacobi) and complex tridiagonal (Thomas) harmonic solves', dof: n, steps: fs.length, iterations: fs.length, residual: 0, notes: 'Fixed–free chain, base excitation, hysteretic damping.' },
  };
}

// ---------- finite aperture ----------
function aperture(p) {
  const nx = Math.round(p.nx), ny = Math.round(p.ny), d = p.pitch, lam = N.C0 / p.freq, kd = 2 * N.PI * d / lam, th0 = p.steer * N.PI / 180, s0 = Math.sin(th0), bits = Math.round(p.bits), rnd = N.rng(p.seed), warnings = [];
  const tap = (i, n) => { const u = n > 1 ? (i - (n - 1) / 2) / ((n - 1) / 2) : 0; return p.taper === 'cosine' ? 0.2 + 0.8 * Math.cos(N.PI * u / 2) : p.taper === 'triangular' ? 1 - Math.abs(u) * (n - 1) / (n + 1) : 1; };
  const sp = p.phase_err * N.PI / 180, sa = p.amp_err / 100, pf = p.fail / 100, q = bits > 0 ? 2 * N.PI / 2 ** bits : 0, sq = q / Math.sqrt(12);
  const wr = new Float64Array(nx * ny), wi = new Float64Array(nx * ny), ir = new Float64Array(nx), cr = new Float64Array(nx), ci = new Float64Array(nx); let pw = 0, pw0 = 0, nFail = 0;
  for (let m = 0; m < nx; m++) for (let n = 0; n < ny; n++) {
    const a0 = tap(m, nx) * tap(n, ny), ph0 = -kd * m * s0; let ph = q > 0 ? Math.round(ph0 / q) * q : ph0; const dead = rnd() < pf, a = dead ? 0 : Math.max(0, a0 * (1 + sa * rnd.normal())); ph += sp * rnd.normal(); if (dead) nFail++;
    wr[m + nx * n] = a * Math.cos(ph); wi[m + nx * n] = a * Math.sin(ph); cr[m] += wr[m + nx * n]; ci[m] += wi[m + nx * n]; ir[m] += a0; pw += a * a; pw0 += a0 * a0;
  }
  const nT = 721, ths = N.linspace(-90, 90, nT), pat = new Float64Array(nT), pat0 = new Float64Array(nT);
  for (let t = 0; t < nT; t++) { const s = Math.sin(ths[t] * N.PI / 180); let re = 0, im = 0, r0 = 0, i0 = 0; for (let m = 0; m < nx; m++) { const ang = kd * m * s, c = Math.cos(ang), sn = Math.sin(ang); re += cr[m] * c - ci[m] * sn; im += cr[m] * sn + ci[m] * c; const a2 = kd * m * (s - s0); r0 += ir[m] * Math.cos(a2); i0 += ir[m] * Math.sin(a2); } pat[t] = re * re + im * im; pat0[t] = r0 * r0 + i0 * i0; }
  const analyse = P => { const ip = N.argmax(P), pk = P[ip]; let l = ip, r = ip; while (l > 0 && P[l - 1] < P[l]) l--; while (r < nT - 1 && P[r + 1] < P[r]) r++; let sl = 0; for (let t = 0; t < nT; t++) if (t < l || t > r) sl = Math.max(sl, P[t]); const cross = dir => { let t = ip; while (t > 0 && t < nT - 1 && P[t] > pk / 2) t += dir; if (P[t] > pk / 2 || t === ip) return ths[t]; const t2 = t - dir; return P[t2] === P[t] ? ths[t] : ths[t2] + (ths[t] - ths[t2]) * (P[t2] - pk / 2) / (P[t2] - P[t]); }; return { peak: ths[ip], pk, hpbw: Math.abs(cross(1) - cross(-1)), sll: sl / pk }; };
  const A = analyse(pat), A0 = analyse(pat0), gl = d / lam > 1 / (1 + Math.abs(s0)), dMaxScan = Math.asin(N.clamp(lam / d - 1, 0, 1)) * 180 / N.PI;
  const glAngles = []; for (const m of [-2, -1, 1, 2]) { const sg = s0 + m * lam / d; if (Math.abs(sg) <= 1) glAngles.push(Math.asin(sg) * 180 / N.PI); }
  const sumA = N.sum(ir), etaTaper = sumA * sumA / (nx * ny * pw0), area = nx * ny * d * d, dirIdeal = 4 * N.PI * area / (lam * lam) * etaTaper * Math.cos(th0), dirCap = Math.min(dirIdeal, N.PI * nx * ny * 2);
  const ruze = (1 - pf) * Math.exp(-(sp * sp + sq * sq)) / (1 + sa * sa), numEff = pw > 0 ? (A.pk / pw) / (A0.pk / pw0) : 0, lossNum = -dB(numEff), lossRuze = -dB(ruze), gain = dB(dirCap) - lossNum;
  const hpbwTh = 0.886 * lam / (nx * d * Math.cos(th0)) * 180 / N.PI, slFloor = dB((sp * sp + sq * sq + sa * sa + pf) / Math.max(ruze * nx * ny, 1e-9));
  const phaseGrad = -kd * s0 * 180 / N.PI;
  if (gl) warnings.push(`Element pitch ${F(d / lam, 3)} λ exceeds the grating-lobe limit ${F(1 / (1 + Math.abs(s0)), 3)} λ for ${F(p.steer)}° steering: a grating lobe appears at ${glAngles.map(v => F(v, 3) + '°').join(', ') || 'the horizon'} and removes power from the main beam.`);
  if (A.hpbw > 120) warnings.push('The aperture is electrically small: the main beam fills most of the visible half-space, so beamwidth, sidelobe level and steering are not meaningful.');
  if (d / lam < 0.2) warnings.push('Elements are closer than 0.2 λ: strong mutual coupling makes the independent-element array factor optimistic; the directivity is limited by the aperture area, not the element count.');
  if (bits === 1) warnings.push('One-bit phase control produces a symmetric mirror beam for a plane-wave feed (quantisation lobe) and about 3–4 dB loss.');
  warnings.push('Array-factor model: isotropic, uncoupled elements; element pattern, mutual coupling, feed loss and edge diffraction are not included.');
  const nu = 41, map = new Float64Array(nu * nu), us = N.linspace(-1, 1, nu); let mapMax = 0;
  const ex = us.map(u => { const c = new Float64Array(nx), s = new Float64Array(nx); for (let m = 0; m < nx; m++) { c[m] = Math.cos(kd * m * u); s[m] = Math.sin(kd * m * u); } return [c, s]; }), ey = us.map(v => { const c = new Float64Array(ny), s = new Float64Array(ny); for (let n = 0; n < ny; n++) { c[n] = Math.cos(kd * n * v); s[n] = Math.sin(kd * n * v); } return [c, s]; });
  for (let j = 0; j < nu; j++) for (let i = 0; i < nu; i++) { if (us[i] ** 2 + us[j] ** 2 > 1) { map[i + nu * j] = -1; continue; } let re = 0, im = 0; const [cx, sx] = ex[i], [cy, sy] = ey[j]; for (let n = 0; n < ny; n++) { let rr = 0, ri = 0; const o = nx * n; for (let m = 0; m < nx; m++) { rr += wr[o + m] * cx[m] - wi[o + m] * sx[m]; ri += wr[o + m] * sx[m] + wi[o + m] * cx[m]; } re += rr * cy[n] - ri * sy[n]; im += rr * sy[n] + ri * cy[n]; } const v = re * re + im * im; map[i + nu * j] = v; mapMax = Math.max(mapMax, v); }
  for (let i = 0; i < map.length; i++) map[i] = map[i] < 0 ? -50 : Math.max(-50, dB(map[i] / (mapMax || 1)));
  const mm = 1e3, sz = 0.4 * d * mm, polys = []; for (let n = 0; n < ny; n++) for (let m = 0; m < nx; m++) { const x = m * d * mm, y = n * d * mm; polys.push({ layer: wr[m + nx * n] === 0 && wi[m + nx * n] === 0 ? 2 : 1, pts: [[x - sz, y - sz], [x + sz, y - sz], [x + sz, y + sz], [x - sz, y + sz]] }); }
  return {
    summary: `${nx}×${ny}-element aperture at ${F(p.freq / 1e9, 4)} GHz (pitch ${F(d / lam, 3)} λ, ${F(Math.sqrt(area) / lam, 3)} λ across): beam at ${F(A.peak, 4)}° (commanded ${F(p.steer)}°), half-power beamwidth ${F(A.hpbw, 3)}°, peak sidelobe ${F(dB(A.sll), 3)} dB (error-free ${F(dB(A0.sll), 3)} dB). Phase, amplitude, quantisation and failure errors cost ${F(lossNum, 3)} dB of gain (Ruze-type estimate ${F(lossRuze, 3)} dB), leaving about ${F(gain, 4)} dBi. Grating lobes are ${gl ? 'PRESENT' : 'absent'}.`,
    kpis: [
      kpi('beam', 'Realised beam direction', A.peak, '°', `Peak of the computed array factor in the steering plane with the seeded error realisation; pointing error ${F(A.peak - p.steer, 3)}°.`, Math.abs(A.peak - p.steer) < 0.25 * A.hpbw),
      kpi('hpbw', 'Half-power beamwidth', A.hpbw, '°', `Width of the main lobe at −3 dB from the pattern; the uniform-aperture estimate 0.886λ/(Nd·cos θ₀) is ${F(hpbwTh, 3)}°.`, null),
      kpi('sll', 'Peak sidelobe level', dB(A.sll), 'dB', `Highest lobe outside the main beam relative to the peak, including random-error sidelobes (design value without errors: ${F(dB(A0.sll), 3)} dB).`, dB(A.sll) < -10),
      kpi('gain', 'Estimated gain', gain, 'dBi', `Aperture directivity 4πA/λ²·η_taper·cos θ₀ = ${F(dB(dirCap), 4)} dBi minus the computed error loss; element pattern and ohmic loss are not included.`, null),
      kpi('loss_num', 'Gain loss from errors (computed)', lossNum, 'dB', 'Reduction of peak directive gain of this realisation relative to the error-free array, normalised to radiated power.', lossNum < 1),
      kpi('loss_ruze', 'Gain loss from errors (Ruze-type)', lossRuze, 'dB', `Expected value (1 − P_fail)·exp(−σ_φ²)/(1 + σ_a²) with rms phase error ${F(Math.sqrt(sp * sp + sq * sq) * 180 / N.PI, 3)}° including quantisation.`, lossRuze < 1),
      kpi('eta_taper', 'Taper efficiency', etaTaper, '–', 'Aperture efficiency of the amplitude taper, |Σa|²/(NΣa²): 1 for uniform illumination.', null),
      kpi('d_lam', 'Pitch in wavelengths', d / lam, '–', `Grating lobes are avoided for d/λ < 1/(1 + |sin θ₀|) = ${F(1 / (1 + Math.abs(s0)), 3)}.`, !gl),
      kpi('scan_max', 'Maximum grating-lobe-free scan angle', dMaxScan, '°', 'Largest steering angle for which no grating lobe enters visible space at this pitch, asin(λ/d − 1).', dMaxScan >= Math.abs(p.steer)),
      kpi('grad', 'Phase gradient', phaseGrad, '° per element', 'Progressive inter-element phase −k·d·sin θ₀ that steers the beam (the generalised Snell law for a reflective or transmissive metasurface).', null),
      kpi('floor', 'Mean error sidelobe floor', slFloor, 'dB', 'Average sidelobe level produced by random errors, about (σ_φ² + σ_a² + P_fail)/(η·N) relative to the peak; it limits how low designed sidelobes can be.', null),
      kpi('n_fail', 'Failed elements', nFail, '–', `Number of dead elements in this seeded realisation (${F(100 * nFail / (nx * ny), 3)} % of ${nx * ny}).`, null),
    ],
    series: [
      { key: 'pat', title: 'Radiation pattern in the steering plane', xlabel: 'angle from broadside (°)', ylabel: 'relative power (dB)', lines: [line('with errors', ths, Array.from(pat, v => Math.max(-60, dB(v / A0.pk * pw0 / (pw || 1))))), line('error-free', ths, Array.from(pat0, v => Math.max(-60, dB(v / A0.pk))))], explain: 'Array factor normalised to the error-free peak at equal radiated power. Errors lower the peak, fill the nulls and raise the far sidelobes; a grating lobe shows as a second full-height beam.' },
      { key: 'ruze', title: 'Gain loss versus rms phase error', xlabel: 'rms phase error (°)', ylabel: 'gain loss (dB)', lines: [line('exp(−σ²)', N.linspace(0, 90, 46), N.linspace(0, 90, 46).map(v => 4.343 * (v * N.PI / 180) ** 2))], explain: 'Ruze relation: loss grows with the square of the rms phase error (about 1 dB at 27°). Use it to set fabrication tolerances.' },
    ],
    maps: [{ key: 'uv', title: 'Pattern in direction-cosine (u, v) space', nx: nu, ny: nu, data: map, x0: -1, x1: 1, y0: -1, y1: 1, xlabel: 'u = sin θ cos φ', ylabel: 'v = sin θ sin φ', unit: 'dB', explain: 'Normalised far-field power of the full two-dimensional aperture with errors; directions outside the unit circle are not visible and shown at the floor value.' }],
    equations: [
      eqn('formula', 'Array factor', 'AF(θ) = Σₘ Σₙ aₘₙ exp[j(k d m sin θ + φₘₙ)],  φₘₙ = −k d m sin θ₀'),
      eqn('formula', 'Grating-lobe condition', 'd/λ < 1/(1 + |sin θ₀|)'),
      eqn('formula', 'Beamwidth and directivity', 'HPBW ≈ 0.886 λ/(N d cos θ₀);  D ≈ (4πA/λ²)·η_taper·cos θ₀'),
      eqn('formula', 'Random-error gain loss (Ruze)', 'G/G₀ ≈ (1 − P_fail)·exp(−σ_φ²)/(1 + σ_a²),  σ_q = (2π/2ᵇ)/√12'),
    ],
    geometry: { kind: 'layout', polygons: polys, units: 'mm' },
    checks: [
      chk('No grating lobes in visible space', !gl, `d/λ = ${F(d / lam, 3)}, limit ${F(1 / (1 + Math.abs(s0)), 3)}`),
      chk('Beam points where commanded', Math.abs(A.peak - p.steer) < Math.max(0.5, 0.25 * A.hpbw), `error ${F(A.peak - p.steer, 3)}° (beamwidth ${F(A.hpbw, 3)}°)`),
      chk('Computed error loss consistent with Ruze estimate', Math.abs(lossNum - lossRuze) < 1 + 0.5 * lossRuze, `${F(lossNum, 3)} dB vs ${F(lossRuze, 3)} dB (single realisation)`),
      chk('Error-free beamwidth matches aperture theory', p.taper !== 'uniform' || gl || Math.abs(A0.hpbw - hpbwTh) < 0.15 * hpbwTh + 0.5, `${F(A0.hpbw, 3)}° vs ${F(hpbwTh, 3)}°`),
      chk('Error-free uniform sidelobe near −13.3 dB', p.taper !== 'uniform' || gl || nx < 6 || Math.abs(dB(A0.sll) + 13.3) < 1.2, `${F(dB(A0.sll), 3)} dB`),
    ],
    warnings,
    recommendations: [
      rec('high', gl ? `Reduce the pitch below ${F(lam / (1 + Math.abs(s0)) * 1e3, 3)} mm or restrict steering to ±${F(dMaxScan, 3)}° to remove the grating lobe.` : `The pitch supports grating-lobe-free steering to ±${F(dMaxScan, 3)}°.`, 'Grating-lobe condition.'),
      rec('high', lossRuze > 1 ? `Errors cost ${F(lossRuze, 3)} dB: hold the rms phase error below about 20° (0.5 dB) by tightening the critical dimension tolerance or adding phase bits (currently ${bits || 'continuous'}).` : 'Manufacturing tolerances are adequate: error loss stays below 1 dB.', 'Ruze-type tolerance budget.'),
      rec('medium', dB(A.sll) > -13 && p.taper === 'uniform' ? 'Apply an amplitude taper if sidelobes matter; a cosine-on-pedestal taper lowers them by about 8 dB for roughly 1 dB of gain.' : 'Sidelobes are governed by the taper and the random-error floor; further taper gives diminishing returns once the floor is reached.', `Peak sidelobe ${F(dB(A.sll), 3)} dB, error floor ${F(slFloor, 3)} dB.`),
      rec('low', `A ${F(100 * pf, 3)} % element failure rate costs only ${F(-dB(1 - pf), 3)} dB: graceful degradation means panels can stay in service rather than being scrapped for isolated defects, which improves fabrication yield and reduces waste.`, 'Failure tolerance of large arrays.', true),
    ],
    publish: { 'ant.gain_dBi': gain, 'ant.hpbw_deg': A.hpbw, 'ant.sll_dB': dB(A.sll), 'ant.efficiency': ruze * etaTaper, 'fab.yield': 1 - nFail / (nx * ny), 'res.f0': p.freq },
    numerics: { method: 'Direct summation of the array factor with seeded random excitation errors', dof: nx * ny, steps: nT + nu * nu, iterations: 1, residual: 0, notes: 'Single error realisation; isotropic uncoupled elements.' },
  };
}

// ---------- reactor thermal network + simulation cost ----------
function reactor(p) {
  const warnings = [], V = p.V_r, Dm = Math.cbrt(4 * V / (N.PI * p.HD)), Hh = p.HD * Dm, A = N.PI * Dm * Hh + N.PI * Dm * Dm / 2, wl = WALLS[p.wall_mat], tw = p.wall_mm * 1e-3, ti = p.ins_mm * 1e-3;
  const Cf = 1000 * 4186 * V, Cw = wl.rho * wl.cp * A * tw, Ri = 1 / (p.h_int * A) + tw / (2 * wl.k * A), Ro = tw / (2 * wl.k * A) + ti / (0.04 * A) + 1 / (p.h_ext * A);
  const Ta = p.T_amb, dTset = p.T_set - Ta, Tinf = Ta + p.P_heat * Ro, reach = Tinf > p.T_set + 0.5 && dTset > 0, tau1 = (Cf + Cw) * Ro + Cf * Ri, tau2 = Ri * Cf * Cw / (Cf + Cw);
  const rhs = (_, y) => [(y[1] - y[0]) / (Ri * Cf), (p.P_heat - (y[1] - y[0]) / Ri - (y[1] - Ta) / Ro) / Cw];
  const tAd = (Cf + Cw) * Math.max(dTset, 1) / p.P_heat, tEnd = Math.min(5 * tau1, Math.max(4 * tAd, 20 * tau2)), nSt = N.clamp(Math.ceil(tEnd / (0.4 * tau2)), 400, 6000), sol = N.rk4(rhs, [Ta, Ta], 0, tEnd, nSt);
  const Tf = sol.y.map(v => v[0]), Tw = sol.y.map(v => v[1]); let tHeat = NaN; for (let i = 1; i < Tf.length; i++) if (Tf[i] >= p.T_set) { tHeat = sol.t[i - 1] + (sol.t[i] - sol.t[i - 1]) * (p.T_set - Tf[i - 1]) / (Tf[i] - Tf[i - 1]); break; }
  if (dTset <= 0) tHeat = 0;
  if (!Number.isFinite(tHeat)) { tHeat = reach ? -tau1 * Math.log(1 - dTset / (p.P_heat * Ro)) : 3 * tau1; }
  let dTwf = 0; for (let i = 0; i < Tf.length && sol.t[i] <= tHeat + 1e-9; i++) dTwf = Math.max(dTwf, Tw[i] - Tf[i]);
  const TwMax = dTset > 0 ? N.interp1(sol.t, Tw, tHeat) : Ta, Phold = Math.max(0, dTset) / Ro, Eheat = p.P_heat * tHeat / 3.6e6;
  const kw = 0.65, alw = 1.6e-7, nuw = 4e-7, betaw = 5e-4, Lc = V / A, Bi = p.h_int * Lc / kw, Fo = alw * tHeat / (Dm / 2) ** 2, tCond = (Dm / 2) ** 2 / (N.PI ** 2 * alw), Ra = 9.81 * betaw * Math.max(dTwf, 0.1) * Hh ** 3 / (nuw * alw), uConv = Math.sqrt(9.81 * betaw * Math.max(dTwf, 0.1) * Hh), tMix = 5 * Hh / uConv;
  const Tc = p.T_set, ant = Tc < 100 ? [8.07131, 1730.63, 233.426] : [8.14019, 1810.94, 244.485], Psat = 133.322 * 10 ** (ant[0] - ant[1] / (ant[2] + Tc)), Pg = Math.max(0, Psat - 101325), hoop = Pg * Dm / (2 * tw), allow = p.wall_mat === 'glass' ? 7e6 : 1.15e8;
  // 3-D simulation cost
  const dxm = p.mesh_mm * 1e-3, cost = r => { const dx = dxm / r, cells = V / dx ** 3, dof = 6 * cells, dt = 0.5 * dx / uConv, steps = p.t_proc / dt, coreS = dof * steps * p.cpu_us * 1e-6, wall = coreS / (0.7 * p.cores); return { dx, cells, dof, mem: dof * 1000, dt, steps, coreH: coreS / 3600, wallH: wall / 3600 }; }, c1 = cost(1), c2 = cost(2), c4 = cost(4);
  const blayer = Hh / Math.pow(Math.max(Ra, 1e4), 0.25), resolved = dxm < blayer / 3;
  if (!reach && dTset > 0) warnings.push(`Heater power is insufficient: steady-state temperature is ${F(Tinf, 4)} °C, below the ${F(p.T_set)} °C set-point. At least ${F(Phold, 4)} W is needed just to hold it.`);
  if (Tc > 100) warnings.push(`Above 100 °C the vessel is pressurised to about ${F(Psat / 1e5, 3)} bar (water saturation); a pressure-rated vessel and relief device are required.`);
  if (hoop > allow) warnings.push(`Hoop stress ${F(hoop / 1e6, 3)} MPa exceeds the indicative allowable of ${F(allow / 1e6, 3)} MPa for ${wl.name}: increase wall thickness. Verify against the applicable pressure-vessel code.`);
  if (!resolved) warnings.push(`The 3-D mesh size (${F(p.mesh_mm)} mm) does not resolve the natural-convection boundary layer (≈ ${F(blayer * 1e3, 3)} mm): wall heat transfer would have to be modelled, not resolved.`);
  let zn = null; try { zn = chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'H2O', phase: 'l' }, { f: 'NaNO3', phase: 'aq' }]); } catch { zn = null; }
  const vols = N.logspace(1e-4, 10, 40), heatV = v => { const d = Math.cbrt(4 * v / (N.PI * p.HD)), a = N.PI * d * d * p.HD + N.PI * d * d / 2, ro = tw / (2 * wl.k * a) + ti / (0.04 * a) + 1 / (p.h_ext * a), P = p.P_heat * v / V, c = 1000 * 4186 * v + wl.rho * wl.cp * a * tw, x = Math.max(dTset, 0) / (P * ro); return x < 1 ? -c * ro * Math.log(1 - x) / 3600 : NaN; };
  const sub = Math.max(1, Math.floor(sol.t.length / 300)), tt = sol.t.filter((_, i) => i % sub === 0), speedup = c1.coreH * 3600 / 1e-2;
  return {
    summary: `${F(V * 1e3, 4)} L reactor (D = ${F(Dm, 3)} m, H = ${F(Hh, 3)} m, ${wl.name} wall) with ${F(p.P_heat)} W: heat-up from ${F(Ta)} to ${F(p.T_set)} °C takes ${F(tHeat / 60, 4)} min (${F(Eheat, 3)} kWh), holding needs ${F(Phold, 4)} W, and the wall runs up to ${F(dTwf, 3)} K hotter than the liquid during the ramp. The equivalent transient 3-D simulation at ${F(p.mesh_mm)} mm resolution has ${F(c1.dof, 3)} DOF, needs ${F(c1.mem / 1e9, 3)} GB and about ${F(c1.wallH, 3)} h on ${F(p.cores)} cores.`,
    kpis: [
      kpi('t_heat', 'Heat-up time', tHeat, 's', `Time for the liquid node of the two-node (liquid + wall) network to reach the set-point at full heater power${reach || dTset <= 0 ? '' : ' — not reachable; three time constants reported'}.`, reach || dTset <= 0),
      kpi('P_hold', 'Holding power', Phold, 'W', 'Steady heat loss through wall, insulation and external film at the set-point: (T_set − T_amb)/R_loss.', Phold < p.P_heat),
      kpi('E_heat', 'Heat-up energy', Eheat, 'kWh', 'Heater energy consumed during the ramp at full power; the adiabatic minimum is ' + F((Cf + Cw) * Math.max(dTset, 0) / 3.6e6, 3) + ' kWh.', null),
      kpi('dT_wall', 'Wall–liquid temperature difference', dTwf, 'K', 'Largest superheat of the wall over the bulk liquid during the ramp. It drives wall nucleation and scaling and measures spatial temperature non-uniformity.', dTwf < 10),
      kpi('T_wall', 'Wall temperature at end of ramp', TwMax + 273.15, 'K', 'Wall-node temperature when the liquid reaches the set-point; the hottest surface the precursor solution sees.', null),
      kpi('tau', 'Dominant thermal time constant', tau1, 's', 'Slow time constant (C_liquid + C_wall)·R_loss + C_liquid·R_wall of the network; the open-loop response time used for controller tuning.', null),
      kpi('Bi', 'Biot number of the charge', Bi, '–', 'h_int·(V/A)/k_liquid: values far above 1 mean that, without convection or stirring, conduction inside the liquid (not the wall film) limits uniformity.', null),
      kpi('t_mix', 'Natural-convection turnover time', tMix, 's', `About five circulation times H/u with u = √(gβΔT·H) = ${F(uConv, 3)} m/s (Rayleigh number ${F(Ra, 3)}); pure conduction would need ${F(tCond, 3)} s (Fourier number at heat-up ${F(Fo, 3)}).`, tMix < 0.2 * tHeat || dTset <= 0),
      kpi('P_sat', 'Vessel pressure at set-point', Psat, 'Pa', 'Saturation pressure of water from the Antoine equation; the autogenous pressure of a closed hydrothermal vessel filled with dilute aqueous solution.', null),
      kpi('hoop', 'Wall hoop stress', hoop, 'Pa', `Thin-wall estimate p·D/(2t) for the gauge pressure, compared with an indicative allowable of ${F(allow / 1e6, 3)} MPa.`, hoop <= allow),
      kpi('dof', '3-D simulation size', c1.dof, 'DOF', `${F(c1.cells, 3)} cells of ${F(p.mesh_mm)} mm with six unknowns each (three velocities, pressure, temperature, one species).`, null),
      kpi('mem', '3-D simulation memory', c1.mem, 'B', 'About 1 kB per degree of freedom for an implicit unstructured finite-volume solver (assumption; state-of-practice range 0.5–2 kB).', c1.mem < 2.56e11),
      kpi('wall_h', '3-D simulation wall-clock time', c1.wallH, 'h', `${F(c1.steps, 3)} time steps of ${F(c1.dt, 3)} s (CFL 0.5 on the convection velocity) at ${F(p.cpu_us)} µs per DOF-step on ${F(p.cores)} cores with 70 % parallel efficiency.`, c1.wallH < 168),
    ],
    series: [
      { key: 'T_t', title: 'Heat-up transient', xlabel: 'time (min)', ylabel: 'temperature (°C)', lines: [line('liquid', tt.map(v => v / 60), Tf.filter((_, i) => i % sub === 0)), line('wall', tt.map(v => v / 60), Tw.filter((_, i) => i % sub === 0)), line('set-point', [0, tEnd / 60], [p.T_set, p.T_set])], explain: 'Open-loop response at full heater power. The wall leads the liquid; the gap between the curves is the driving temperature difference across the internal film.' },
      { key: 'scale', title: 'Heat-up time versus scale at constant specific power', xlabel: 'reactor volume (m³)', ylabel: 'heat-up time (h)', logx: true, lines: [...(vols.some(v => Number.isFinite(heatV(v))) ? [line(`${F(p.P_heat / V / 1e3, 3)} kW/m³`, vols.filter(v => Number.isFinite(heatV(v))), vols.filter(v => Number.isFinite(heatV(v))).map(heatV))] : []), line('adiabatic limit (liquid only)', vols, vols.map(v => 4.186e6 * v * Math.max(dTset, 0) / (p.P_heat * v / V) / 3600))], explain: 'Scale-up at the same power per unit volume: larger vessels lose relatively less heat (smaller surface-to-volume ratio) so the heat-up time approaches the adiabatic value; small vessels may never reach the set-point.' },
      { key: 'cost', title: '3-D simulation cost versus mesh refinement', xlabel: 'cell size (mm)', ylabel: 'core-hours', logx: true, logy: true, lines: [line('core-hours', [c1, c2, c4].map(c => c.dx * 1e3), [c1, c2, c4].map(c => c.coreH))], explain: 'Each halving of the cell size multiplies cells by 8 and time steps by 2, so the cost rises 16-fold: the basis for choosing a reduced-order model for design sweeps.' },
    ],
    tables: [{ key: 'cost', title: '3-D simulation cost estimate', columns: ['Cell size (mm)', 'Cells', 'DOF', 'Memory (GB)', 'Time steps', 'Core-hours', 'Wall-clock (h)'], rows: [c1, c2, c4].map(c => [F(c.dx * 1e3, 3), F(c.cells, 3), F(c.dof, 3), F(c.mem / 1e9, 3), F(c.steps, 3), F(c.coreH, 3), F(c.wallH, 3)]), explain: 'Order-of-magnitude resource estimate for a transient conjugate heat-transfer and flow simulation of the full vessel at three mesh levels, to be calibrated with a short test run of the actual solver.' }],
    equations: [
      eqn('governing', 'Lumped thermal network', 'C_f dT_f/dt = (T_w − T_f)/R_i;   C_w dT_w/dt = P − (T_w − T_f)/R_i − (T_w − T_amb)/R_o'),
      eqn('formula', 'Thermal resistances', 'R_i = 1/(h_i A) + t_w/(2k_w A);   R_o = t_w/(2k_w A) + t_ins/(k_ins A) + 1/(h_o A)'),
      eqn('formula', 'Dimensionless groups', 'Bi = h L_c/k;  Fo = αt/R²;  Ra = gβΔT H³/(να)'),
      eqn('formula', 'Simulation cost scaling', 'N_cells = V/Δ³;  Δt = CFL·Δ/u;  cost ∝ Δ⁻⁴'),
      ...(zn ? [eqn('chemical', 'Reference reaction (ZnO hydrothermal)', zn.molecular, 'process carried out in the reactor'), eqn('ionic', 'Net ionic equation', zn.netIonic)] : []),
    ],
    checks: [
      chk('Set-point reachable with installed power', reach || dTset <= 0, `steady-state limit ${F(Tinf, 4)} °C`),
      chk('Energy balance of the transient', Math.abs((Cf * (Tf[Tf.length - 1] - Ta) + Cw * (Tw[Tw.length - 1] - Ta) + N.trapz(sol.t, Tw.map(v => (v - Ta) / Ro))) / (p.P_heat * tEnd) - 1) < 1e-3, 'stored + lost heat equals heater input within 0.1 %'),
      chk('Time-scale separation for the two-node model', tau2 < 0.3 * tau1, `τ₂/τ₁ = ${F(tau2 / tau1, 3)}`),
      chk('Liquid mixes faster than it heats', tMix < 0.2 * tHeat || dTset <= 0, `turnover ${F(tMix, 3)} s vs heat-up ${F(tHeat, 3)} s`),
      chk('Wall stress within indicative allowable', hoop <= allow, `${F(hoop / 1e6, 3)} MPa vs ${F(allow / 1e6, 3)} MPa`),
      chk('3-D mesh resolves the boundary layer', resolved, `Δ = ${F(p.mesh_mm)} mm, boundary layer ≈ ${F(blayer * 1e3, 3)} mm`),
    ],
    warnings,
    recommendations: [
      rec('high', !reach && dTset > 0 ? `Install at least ${F(1.5 * Phold, 3)} W or add insulation: the present heater cannot reach the set-point.` : `Heat-up takes ${F(tHeat / 60, 3)} min; if batch time matters, more power shortens it almost proportionally, but the wall superheat (now ${F(dTwf, 3)} K) rises too.`, 'Two-node thermal network.'),
      rec('high', dTwf > 10 ? `Wall superheat of ${F(dTwf, 3)} K promotes heterogeneous nucleation and deposits on the wall: stir, ramp at lower power, or use jacket heating with a larger area.` : 'Wall superheat is small; particle nucleation should be governed by the bulk temperature.', 'Wall-to-liquid temperature difference.'),
      rec('medium', `Use the lumped model for design sweeps and control tuning (plant gain ${F(Ro, 3)} K/W, time constant ${F(tau1, 3)} s, apparent delay ${F(tau2 / 2, 3)} s) and reserve the ${F(c1.coreH, 3)} core-hour 3-D simulation for final verification of mixing and hot spots.`, 'Cost scaling of resolved simulations.'),
      rec('medium', `Insulation: holding losses are ${F(Phold, 3)} W (${F(Phold * 24 / 1000, 3)} kWh per day). Doubling the insulation thickness would nearly halve them; recover heat from the hot effluent to pre-heat the next batch.`, 'Steady-state loss through the insulation.', true),
    ],
    publish: { 'proc.time': tHeat, 'proc.hold_power': Phold, 'proc.UA': 1 / Ro, 'proc.P': Psat, 'therm.T_max': TwMax + 273.15, 'proc.volume': V, 'ctrl.plant_gain': Ro, 'ctrl.plant_tau': tau1 + tau2 / 2, 'ctrl.plant_delay': tau2 / 2, 'fab.energy_kWh_batch': Eheat },
    numerics: { method: 'Two-node lumped thermal network integrated by classical RK4; algebraic cost model', dof: 2, steps: nSt, iterations: nSt, residual: 0, notes: `Constant water properties (ρ 1000 kg/m³, c_p 4186 J/kg/K); speed-up of the lumped model over the 3-D estimate ≈ ${F(speedup, 2)}.` },
  };
}

// ---------- reliability block diagram ----------
function reliability(p) {
  const ns = Math.round(p.n_ser), mp = Math.round(p.m_par), b = p.beta_w, eta = p.eta_h, t = p.t_miss, warnings = [];
  const Rc = x => Math.exp(-Math.pow(x / eta, b)), Rst = x => 1 - Math.pow(1 - Rc(x), mp), Rcom = x => Math.exp(-Math.pow(x / p.eta_c, p.beta_c)), Rsys = x => Math.exp(ns * Math.log(Math.max(Rst(x), 1e-300))) * Rcom(x);
  const R = Rsys(t); let tHi = Math.min(eta, p.eta_c); for (let i = 0; i < 60 && Rsys(tHi) > 1e-6; i++) tHi *= 1.5; for (let i = 0; i < 80 && Rsys(tHi) < 1e-9; i++) tHi /= 1.3;
  const ts = N.linspace(0, tHi, 400), Rs = ts.map(Rsys), mttf = N.trapz(ts, Rs);
  const life = q => { try { return N.brent(x => Rsys(x) - q, 1e-12 * tHi, tHi); } catch { return 0; } }, B10 = life(0.9), B1 = life(0.99), B50 = life(0.5);
  const haz = x => { const h = Math.max(x * 1e-4, 1e-9); return -(Math.log(Math.max(Rsys(x + h), 1e-300)) - Math.log(Math.max(Rsys(Math.max(x - h, 0)), 1e-300))) / (x + h - Math.max(x - h, 0)); };
  const etaSer = eta * Math.pow(ns, -1 / b), weakest = mp === 1 ? Math.exp(-Math.pow(t / etaSer, b)) * Rcom(t) : NaN;
  const rnd = N.rng(p.seed), nMC = N.clamp(Math.floor(3e5 / (ns * mp)), 300, 4000); let surv = 0; const lives = [];
  for (let s = 0; s < nMC; s++) { let sysL = p.eta_c * Math.pow(-Math.log(1 - rnd()), 1 / p.beta_c); for (let i = 0; i < ns && sysL > 0; i++) { let st = 0; for (let j = 0; j < mp; j++) { const l = eta * Math.pow(-Math.log(1 - rnd()), 1 / b); if (l > st) st = l; } if (st < sysL) sysL = st; } lives.push(sysL); if (sysL > t) surv++; }
  const Rmc = surv / nMC, seMC = Math.sqrt(Math.max(R * (1 - R), 1e-12) / nMC), mttfMC = N.mean(lives), hzSt = -ns * Math.log(Math.max(Rst(t), 1e-300)), hzCom = Math.pow(t / p.eta_c, p.beta_c), share = hzSt + hzCom > 0 ? hzSt / (hzSt + hzCom) : 0.5;
  const gainPar = mp > 1 ? Rsys(t) - Math.exp(ns * Math.log(Math.max(Rc(t), 1e-300))) * Rcom(t) : 0, RnoRed = Math.exp(ns * Math.log(Math.max(Rc(t), 1e-300))) * Rcom(t), Rred = (1 - Math.pow(1 - Rc(t), mp + 1)) ** ns * Rcom(t);
  if (R < 0.5) warnings.push(`System reliability at the mission time is only ${F(R, 3)}: with ${ns} stages in series the weakest of many components governs. Add redundancy or shorten the service interval.`);
  if (b < 1) warnings.push('Shape parameter below 1 indicates infant mortality: burn-in screening is more effective than redundancy or preventive replacement.');
  const tl = N.linspace(0, tHi, 200);
  return {
    summary: `${ns} stage(s) in series, each with ${mp} parallel Weibull component(s) (β = ${F(b)}, η = ${F(eta)} h), plus a common element: system reliability at ${F(t)} h is ${F(R, 4)} (Monte-Carlo ${F(Rmc, 4)} from ${nMC} seeded histories). MTTF = ${F(mttf, 4)} h, B10 life = ${F(B10, 4)} h (${F(B10 / 8760, 3)} years). The stages account for ${F(100 * share, 3)} % of the failure risk and the common element for the rest.`,
    kpis: [
      kpi('R', 'System reliability at mission time', R, '–', 'Probability that every stage has at least one working component and the common element survives: R = [1 − (1 − R_c)^m]^n · R_common.', R >= 0.9),
      kpi('R_mc', 'Reliability by Monte-Carlo', Rmc, '–', `Fraction of ${nMC} seeded simulated system lives exceeding the mission time; the sampling standard error is ${F(seMC, 2)}.`, Math.abs(Rmc - R) < 4 * seMC + 0.01),
      kpi('mttf', 'Mean time to failure', mttf, 'h', `Integral of the system reliability function (trapezoidal, 400 points); Monte-Carlo mean ${F(mttfMC, 4)} h.`, null),
      kpi('B10', 'B10 life', B10, 'h', 'Time by which 10 % of systems are expected to have failed (R = 0.9): the usual design-life measure.', B10 >= t),
      kpi('B1', 'B1 life', B1, 'h', 'Time by which 1 % of systems have failed; relevant for safety-critical service.', null),
      kpi('B50', 'Median life', B50, 'h', 'Time at which half of the systems have failed.', null),
      kpi('haz', 'Hazard rate at mission time', haz(t), '1/h', 'Instantaneous failure rate −d ln R/dt of the system. A rising hazard indicates wear-out, for which scheduled replacement is effective.', null),
      kpi('eta_ser', 'Series-chain characteristic life', etaSer, 'h', `Weakest-link scaling η·n^(−1/β): ${ns} identical components in series behave as one Weibull unit with this scale (size effect).`, null),
      kpi('R_c', 'Single-component reliability', Rc(t), '–', 'Weibull survival probability exp(−(t/η)^β) of one component at the mission time.', null),
      kpi('share', 'Risk share of the stages', share, '–', 'Fraction of the system log-unreliability contributed by the repeated stages (the remainder is the common element): shows where improvement pays.', null),
      kpi('R_plus', 'Reliability with one more redundant unit', Rred, '–', `System reliability if every stage had ${mp + 1} parallel components instead of ${mp}; without any redundancy it would be ${F(RnoRed, 4)}.`, null),
    ],
    series: [
      { key: 'R_t', title: 'Reliability functions', xlabel: 'time (h)', ylabel: 'reliability', lines: [line('system', tl, tl.map(Rsys)), line('single component', tl, tl.map(Rc)), line('one stage', tl, tl.map(Rst)), line('common element', tl, tl.map(Rcom))], explain: 'Survival probability of the system and of its building blocks. The series system always lies below its weakest stage; redundancy flattens the early part of the stage curve.' },
      { key: 'haz_t', title: 'System hazard rate', xlabel: 'time (h)', ylabel: 'hazard rate (1/h)', logy: true, lines: [line('h(t)', tl.slice(1), tl.slice(1).map(x => Math.max(haz(x), 1e-300)))], explain: 'Failure rate of the surviving population. Decreasing: infant mortality; flat: random failures; increasing: wear-out.' },
    ],
    equations: [
      eqn('formula', 'Weibull component reliability', 'R_c(t) = exp[−(t/η)^β]'),
      eqn('formula', 'Parallel stage and series system', 'R_stage = 1 − (1 − R_c)^m;   R_sys = R_stage^n · R_common'),
      eqn('formula', 'Weakest-link size effect', 'η_n = η · n^(−1/β)'),
      eqn('formula', 'Mean life and hazard', 'MTTF = ∫₀^∞ R(t) dt;   h(t) = −d ln R/dt'),
    ],
    checks: [
      chk('Analytical and Monte-Carlo reliability agree', Math.abs(Rmc - R) < 4 * seMC + 0.01, `difference ${F(Rmc - R, 3)}, 4σ = ${F(4 * seMC, 3)}`),
      chk('Weakest-link formula reproduced (no redundancy)', mp > 1 || Math.abs(weakest - R) < 1e-9, mp > 1 ? 'not applicable with redundancy' : `${F(weakest, 6)} vs ${F(R, 6)}`),
      chk('Reliability within [0, 1] and decreasing', R >= 0 && R <= 1 && Rs.every((v, i) => i === 0 || v <= Rs[i - 1] + 1e-12), `R(0) = ${F(Rs[0], 4)}, R(t_end) = ${F(Rs[Rs.length - 1], 3)}`),
      chk('MTTF integral converged', Rs[Rs.length - 1] < 1e-4, `tail reliability ${Rs[Rs.length - 1].toExponential(1)}`),
      chk('Mission within B10 life', B10 >= t, `B10 = ${F(B10, 4)} h, mission ${F(t)} h`),
    ],
    warnings,
    recommendations: [
      rec('high', share > 0.7 ? (mp === 1 ? `The ${ns} series stages dominate the risk: one redundant component per stage would raise reliability from ${F(R, 3)} to ${F(Rred, 3)}.` : `Stages still dominate; a further redundant unit would raise reliability to ${F(Rred, 3)}, but check common-cause failures, which redundancy does not remove.`) : `The shared element dominates the risk (${F(100 * (1 - share), 3)} %): improve or duplicate it before adding redundancy to the stages.`, 'Risk apportionment of the block diagram.'),
      rec('medium', b > 1.5 ? `Wear-out behaviour (β = ${F(b)}): scheduled replacement at about the B10 life (${F(B10, 3)} h) prevents most failures.` : 'With β near or below 1 scheduled replacement does not help; use condition monitoring or burn-in instead.', 'Shape of the hazard function.'),
      rec('medium', 'Obtain β and η from at least 10–20 component failures (or suspended tests) fitted by median-rank regression; the series size effect η·n^(−1/β) is very sensitive to β.', 'Weibull parameter uncertainty.'),
      rec('low', 'Design for replaceable stages or modules so that a single failed unit is exchanged instead of scrapping the whole array; this extends service life and reduces material waste.', 'Series systems fail at their weakest element.', true),
    ],
    publish: { 'life.years': B10 / 8760, 'life.reliability': R, 'life.mttf_h': mttf, 'risk.index': 1 - R },
    numerics: { method: 'Closed-form reliability block diagram with Brent root finding for life quantiles, trapezoidal MTTF integration and seeded Monte-Carlo verification', dof: ns * mp + 1, steps: 400, iterations: nMC, residual: Math.abs(Rmc - R), notes: 'Independent, identically distributed component lives; no repair, no common-cause failure.' },
  };
}

export default {
  id: 19,
  tagline: 'Analyses full components and systems: finite lattice panels, homogenised plates, resonator arrays, finite apertures, scaled reactors and system reliability.',
  inputs,
  presets: [
    { name: 'Re-entrant auxetic panel in tension', note: 'Aluminium re-entrant honeycomb, 8 × 5 cells: negative Poisson\'s ratio and edge effects versus homogenisation.', values: { mode: 'lattice_fe', ctype: 'reentrant', ncx: 8, ncy: 5, lcase: 'tension', th_re: 20, h_l: 2 } },
    { name: 'Triangular lattice beam in three-point bending', note: 'Stretch-dominated lattice loaded as a beam; compare with Timoshenko beam on homogenised properties.', values: { mode: 'lattice_fe', ctype: 'triangular', ncx: 16, ncy: 4, lcase: 'bend3', F_load: 200, t_strut: 4e-4 } },
    { name: 'Polymer nanolattice honeycomb in compression', note: 'Micro-printed hexagonal honeycomb with 2 µm cells: member buckling governs.', values: { mode: 'lattice_fe', ctype: 'honeycomb', ncx: 10, ncy: 6, lcase: 'compression', cell: 2e-6, t_strut: 1e-7, depth: 2e-6, E_s: 3e9, sy: 7e7, strain: 0.5 } },
    { name: 'Auxetic sandwich panel under pressure', note: 'Homogenised auxetic plate (ν = −0.4): deflection, frequencies and buckling.', values: { mode: 'plate', nu_p: -0.4, E_p: 2e9, rho_p: 400, pa: 0.5, pb: 0.4, ph: 0.01, q_p: 2000 } },
    { name: 'Phononic resonator chain, mass ratio 3', note: 'Ten-cell diatomic chain: finite-array transmission against the infinite Bloch gap.', values: { mode: 'modal', n_cells: 10, m1: 0.01, m2: 0.03, k_s: 1e5, eta: 0.01 } },
    { name: '10 GHz beam-steering metasurface, 3-bit', note: 'Programmable 16 × 16 metasurface with quantised phase and manufacturing errors.', values: { mode: 'array', nx: 16, ny: 16, pitch: 0.015, freq: 1e10, steer: 20, bits: 3, phase_err: 15 } },
    { name: 'ZnO hydrothermal reactor, 50 L at 90 °C', note: 'Scale-up of the reference ZnO nanorod synthesis: heat-up, losses and 3-D simulation cost.', values: { mode: 'reactor', V_r: 0.05, T_set: 90, P_heat: 5000 } },
    { name: 'Metamaterial array reliability, 100 cells', note: 'One hundred unit cells in series with wear-out life and a shared controller.', values: { mode: 'reliability', n_ser: 100, m_par: 1, beta_w: 2.5, eta_h: 2e5, t_miss: 8760 } },
  ],
  run(inp) {
    const p = read(inp);
    const r = p.mode === 'plate' ? plate(p) : p.mode === 'modal' ? modal(p) : p.mode === 'array' ? aperture(p) : p.mode === 'reactor' ? reactor(p) : p.mode === 'reliability' ? reliability(p) : latticeFE(p);
    return finish(r);
  },
};
