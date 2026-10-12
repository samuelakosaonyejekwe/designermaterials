// Suite 10 — Heat transfer and thermal engineering.
// Transient conduction with Joule and reaction heating, lumped / composite-wall / fin models with
// ramp–soak schedules, radiative furnace balance, nanoscale conduction and thermal metamaterials
// (anisotropic laminates, bilayer cloak, negative thermal expansion, thermal rectification / NDTR).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, PI = Math.PI, R = N.R_GAS, SB = N.SIGMA_SB;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
const deg = x => x * 180 / PI, rad = x => x * PI / 180;

const BC = [{ value: 'robin', label: 'Convective (Robin): h, ambient temperature' }, { value: 'fixed', label: 'Fixed temperature (= ambient)' }, { value: 'insulated', label: 'Insulated (zero flux)' }];
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'slab', group: 'Model', options: [
    { value: 'slab', label: 'Transient conduction with Joule and reaction heating' }, { value: 'lumped', label: 'Lumped body, composite wall, fin and ramp–soak schedule' }, { value: 'furnace', label: 'Radiation + convection furnace / wafer balance' },
    { value: 'nano', label: 'Nanoscale conduction and interface resistance' }, { value: 'meta', label: 'Thermal metamaterials: laminate tensor, cloak, concentrator' }, { value: 'nte', label: 'Negative thermal expansion composites and lattices' }, { value: 'rect', label: 'Thermal rectification and negative differential thermal resistance' }], help: 'Selects the thermal model. Results (temperatures, time constants, effective properties) are published for the kinetics, mechanics and control suites.' },
  // ---- slab + lumped shared
  { key: 'k', label: 'Thermal conductivity', unit: 'W/m/K', value: 0.6, min: 0.01, max: 2000, log: true, group: 'Material', bus: 'eff.k', showIf: { mode: ['slab', 'lumped'] }, help: 'Conductivity of the slab or body (0.6 W/m/K for water or a water-saturated porous bed).' },
  { key: 'rho', label: 'Density', unit: 'kg/m³', value: 1000, min: 1, max: 25000, log: true, group: 'Material', showIf: { mode: ['slab', 'lumped'] }, help: 'Bulk density of the medium being heated.' },
  { key: 'cp', label: 'Specific heat capacity', unit: 'J/kg/K', value: 4186, min: 50, max: 20000, log: true, group: 'Material', showIf: { mode: ['slab', 'lumped'] }, help: 'Heat capacity per unit mass (4186 J/kg/K for water).' },
  { key: 'T0', label: 'Initial temperature', unit: 'K', value: 298, min: 50, max: 2500, group: 'Conditions', showIf: { mode: ['slab', 'lumped'] }, help: 'Uniform temperature at the start.' },
  { key: 'h', label: 'Convection coefficient', unit: 'W/m²/K', value: 50, min: 0.1, max: 1e5, log: true, group: 'Conditions', showIf: { mode: ['slab', 'lumped'] }, help: 'Heat-transfer coefficient to the surroundings: about 5–25 for natural convection in air, 50–1000 for forced liquids.' },
  { key: 'Tamb', label: 'Ambient temperature', unit: 'K', value: 298, min: 50, max: 2500, group: 'Conditions', showIf: { mode: ['slab', 'lumped'] }, help: 'Temperature of the surroundings or coolant.' },
  // ---- slab
  { key: 'dim', label: 'Dimensionality', type: 'select', value: '1d', group: 'Geometry', showIf: { mode: ['slab'] }, options: [{ value: '1d', label: '1-D through the thickness' }, { value: '2d', label: '2-D cross-section (width × thickness)' }], help: 'The 1-D model resolves the through-thickness profile. The 2-D model adds the width, edge losses and partial electrode coverage.' },
  { key: 'Lz', label: 'Slab thickness', unit: 'm', value: 0.01, min: 1e-6, max: 5, log: true, group: 'Geometry', bus: 'geom.thickness', showIf: { mode: ['slab'] }, help: 'Thickness between the two cooled faces.' },
  { key: 'Lxw', label: 'Slab width', unit: 'm', value: 1, min: 1e-5, max: 50, log: true, group: 'Geometry', showIf: { mode: ['slab'] }, help: 'In-plane width (resolved in the 2-D model).' },
  { key: 'Lyd', label: 'Slab depth', unit: 'm', value: 1, min: 1e-5, max: 50, log: true, group: 'Geometry', showIf: { mode: ['slab'] }, help: 'Second in-plane dimension, used for total power and energy.' },
  { key: 'Ef', label: 'Electric field', unit: 'V/m', value: 1000, min: 0, max: 1e6, group: 'Heat sources', bus: 'proc.E_field', showIf: { mode: ['slab'] }, help: 'Field in the conducting medium. Joule heating is σE² per unit volume.' },
  { key: 'sig', label: 'Electrical conductivity', unit: 'S/m', value: 0.25, min: 0, max: 1e7, group: 'Heat sources', bus: 'eff.sigma', showIf: { mode: ['slab'] }, help: 'Conductivity of the medium at the initial temperature (about 0.1–1 S/m for dilute electrolytes).' },
  { key: 'alS', label: 'Temperature coefficient of conductivity', unit: '1/K', value: 0.02, min: -0.01, max: 0.1, group: 'Heat sources', showIf: { mode: ['slab'] }, help: 'σ(T) = σ₀[1 + α(T − T₀)]. Electrolytes gain about 2 % per kelvin, which feeds back positively on Joule heating; metals have a negative coefficient.' },
  { key: 'qR', label: 'Reaction heat release at the initial temperature', unit: 'W/m³', value: 100, min: -1e7, max: 1e7, group: 'Heat sources', showIf: { mode: ['slab'] }, help: 'Volumetric heat of reaction at T₀ (positive = exothermic), scaled with an Arrhenius factor as temperature changes.' },
  { key: 'EaR', label: 'Activation energy of the reaction heat source', unit: 'kJ/mol', value: 65, min: 0, max: 300, group: 'Heat sources', showIf: { mode: ['slab'] }, help: 'Arrhenius activation energy of the heat-releasing reaction (65 kJ/mol in the ZnO reference case).' },
  { key: 'bcA', label: 'Boundary at face z = 0', type: 'select', value: 'robin', group: 'Conditions', showIf: { mode: ['slab'] }, options: BC, help: 'Thermal boundary condition on the lower face.' },
  { key: 'bcB', label: 'Boundary at face z = L', type: 'select', value: 'robin', group: 'Conditions', showIf: { mode: ['slab'] }, options: BC, help: 'Thermal boundary condition on the upper face.' },
  { key: 'tEnd', label: 'Simulated time', unit: 's', value: 3600, min: 1e-3, max: 1e7, log: true, group: 'Conditions', showIf: { mode: ['slab'] }, help: 'Duration of the transient.' },
  { key: 'Tset', label: 'Setpoint temperature', unit: 'K', value: 320, min: 50, max: 2500, group: 'Conditions', showIf: { mode: ['slab'] }, help: 'Target mean temperature; the time to reach it is reported.' },
  { key: 'hfrac', label: 'Electrode coverage of the width (2-D)', unit: '', value: 0.5, min: 0.05, max: 1, group: 'Heat sources', showIf: { mode: ['slab'] }, help: 'Central fraction of the width in which the field (and Joule heating) is present in the 2-D model.' },
  // ---- lumped
  { key: 'Vl', label: 'Body volume', unit: 'm³', value: 1e-3, min: 1e-9, max: 100, log: true, group: 'Geometry', bus: 'proc.volume', showIf: { mode: ['lumped'] }, help: 'Volume of the heated contents (10⁻³ m³ = 1 L).' },
  { key: 'Al', label: 'Heat-loss surface area', unit: 'm²', value: 0.06, min: 1e-6, max: 1e3, log: true, group: 'Geometry', showIf: { mode: ['lumped'] }, help: 'External area through which heat is lost.' },
  { key: 'Ph', label: 'Heater power', unit: 'W', value: 300, min: 0, max: 1e7, group: 'Schedule', bus: 'proc.power', showIf: { mode: ['lumped'] }, help: 'Constant power delivered to the contents during the ramp.' },
  { key: 'Tsoak', label: 'Soak temperature', unit: 'K', value: 363, min: 50, max: 2500, group: 'Schedule', showIf: { mode: ['lumped'] }, help: 'Hold temperature of the ramp–soak schedule (363 K = 90 °C).' },
  { key: 'tSoak', label: 'Soak time', unit: 's', value: 10800, min: 0, max: 1e7, group: 'Schedule', bus: 'proc.time', showIf: { mode: ['lumped'] }, help: 'Duration of the hold at the soak temperature.' },
  { key: 'Tcool', label: 'Cool-down target temperature', unit: 'K', value: 313, min: 50, max: 2500, group: 'Schedule', showIf: { mode: ['lumped'] }, help: 'Temperature at which the vessel can be opened or unloaded.' },
  { key: 'Lw', label: 'Wall thickness', unit: 'm', value: 0.005, min: 0, max: 1, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: 'Thickness of the vessel or furnace wall.' },
  { key: 'kw', label: 'Wall conductivity', unit: 'W/m/K', value: 16, min: 0.01, max: 500, log: true, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: '16 W/m/K for stainless steel, about 1 for glass or ceramics.' },
  { key: 'Lins', label: 'Insulation thickness', unit: 'm', value: 0.02, min: 0, max: 1, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: 'Thickness of the insulation layer outside the wall (zero for none).' },
  { key: 'kins', label: 'Insulation conductivity', unit: 'W/m/K', value: 0.04, min: 0.005, max: 5, log: true, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: '0.04 W/m/K for mineral wool, 0.02 for aerogel blankets.' },
  { key: 'dfin', label: 'Pin-fin diameter', unit: 'm', value: 0.005, min: 1e-7, max: 0.1, log: true, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: 'Diameter of a cylindrical fin (or of a nanowire acting as a fin).' },
  { key: 'Lfin', label: 'Pin-fin length', unit: 'm', value: 0.05, min: 1e-7, max: 1, log: true, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: 'Length of the fin from its base.' },
  { key: 'kfin', label: 'Fin conductivity', unit: 'W/m/K', value: 200, min: 0.1, max: 2500, log: true, group: 'Wall and fin', showIf: { mode: ['lumped'] }, help: 'Conductivity of the fin material (about 200 W/m/K for aluminium).' },
  // ---- furnace
  { key: 'Th', label: 'Heater (lamp / susceptor) temperature', unit: 'K', value: 1273, min: 300, max: 3500, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Temperature of the radiating heater facing the wafer.' },
  { key: 'epsW', label: 'Wafer emissivity', unit: '', value: 0.65, min: 0.02, max: 1, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Total hemispherical emissivity of the wafer or sample (about 0.65–0.7 for hot silicon; films and metasurfaces change it strongly).' },
  { key: 'epsH', label: 'Heater emissivity', unit: '', value: 0.85, min: 0.02, max: 1, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Emissivity of the heater surface (graphite or SiC susceptors: 0.8–0.9).' },
  { key: 'Twall', label: 'Chamber wall temperature', unit: 'K', value: 300, min: 50, max: 2500, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Temperature of the walls seen by the back of the wafer (cold-wall reactor: near ambient).' },
  { key: 'hg', label: 'Gas convection coefficient', unit: 'W/m²/K', value: 20, min: 0, max: 2000, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Convective coefficient between the wafer faces and the process gas.' },
  { key: 'Tg', label: 'Gas temperature', unit: 'K', value: 600, min: 50, max: 2500, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Bulk temperature of the process gas.' },
  { key: 'tw', label: 'Wafer thickness', unit: 'm', value: 7.25e-4, min: 1e-6, max: 0.05, log: true, group: 'Furnace', bus: 'geom.thickness', showIf: { mode: ['furnace'] }, help: 'Thickness of the wafer or sample (725 µm for a 200 mm silicon wafer).' },
  { key: 'rcW', label: 'Wafer volumetric heat capacity', unit: 'J/m³/K', value: 2.0e6, min: 1e4, max: 1e7, log: true, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'ρ·c_p of the wafer (about 2×10⁶ J/m³/K for silicon at high temperature).' },
  { key: 'kWf', label: 'Wafer thermal conductivity', unit: 'W/m/K', value: 30, min: 0.1, max: 2000, log: true, group: 'Furnace', showIf: { mode: ['furnace'] }, help: 'Used only for the Biot-number check (silicon: 150 W/m/K at 300 K falling to about 30 at 1000 K).' },
  // ---- nano
  { key: 'kb', label: 'Bulk thermal conductivity', unit: 'W/m/K', value: 148, min: 0.05, max: 3000, log: true, group: 'Carrier transport', bus: 'prop.k_thermal', showIf: { mode: ['nano'] }, help: 'Bulk conductivity of the material (148 W/m/K for silicon, about 50 for ZnO single crystal).' },
  { key: 'Cv', label: 'Volumetric heat capacity of the carriers', unit: 'J/m³/K', value: 1.66e6, min: 1e3, max: 1e7, log: true, group: 'Carrier transport', showIf: { mode: ['nano'] }, help: 'Phonon (or electron) heat capacity per unit volume (1.66×10⁶ J/m³/K for silicon at 300 K).' },
  { key: 'vs', label: 'Carrier group velocity', unit: 'm/s', value: 6400, min: 100, max: 2e6, log: true, group: 'Carrier transport', bus: 'prop.c_sound', showIf: { mode: ['nano'] }, help: 'Average phonon group (sound) velocity; for metals use the Fermi velocity with the electronic heat capacity.' },
  { key: 'dch', label: 'Film thickness / wire diameter', unit: 'm', value: 5e-8, min: 1e-9, max: 1e-3, log: true, group: 'Carrier transport', bus: 'geom.diameter', showIf: { mode: ['nano'] }, help: 'Characteristic size of the nanostructure (5×10⁻⁸ m = 50 nm).' },
  { key: 'spec', label: 'Boundary specularity p', unit: '', value: 0, min: 0, max: 0.99, group: 'Carrier transport', showIf: { mode: ['nano'] }, help: 'Fraction of carriers reflected specularly at the surface: 0 for rough (fully diffuse) surfaces, towards 1 for atomically smooth ones.' },
  { key: 'RK', label: 'Kapitza (interface) resistance', unit: 'm²K/W', value: 1e-8, min: 1e-10, max: 1e-5, log: true, group: 'Nanocomposite', showIf: { mode: ['nano'] }, help: 'Thermal boundary resistance between particle and matrix: about 10⁻⁹–10⁻⁸ for solid–solid, 10⁻⁸–10⁻⁷ for solid–liquid or weakly bonded interfaces.' },
  { key: 'kpar', label: 'Particle conductivity', unit: 'W/m/K', value: 50, min: 0.01, max: 3000, log: true, group: 'Nanocomposite', showIf: { mode: ['nano'] }, help: 'Conductivity of the dispersed particles (use the size-reduced value for very small particles).' },
  { key: 'kmat', label: 'Matrix conductivity', unit: 'W/m/K', value: 0.6, min: 0.01, max: 500, log: true, group: 'Nanocomposite', showIf: { mode: ['nano'] }, help: 'Conductivity of the continuous phase (0.6 W/m/K for water, about 0.2 for polymers).' },
  { key: 'phiP', label: 'Particle volume fraction', unit: '', value: 0.05, min: 0, max: 0.5, group: 'Nanocomposite', showIf: { mode: ['nano'] }, help: 'Volume fraction of spherical particles; effective-medium formulae are reliable up to about 0.2–0.3.' },
  { key: 'aP', label: 'Particle radius', unit: 'm', value: 2e-8, min: 5e-10, max: 1e-3, log: true, group: 'Nanocomposite', showIf: { mode: ['nano'] }, help: 'Radius of the dispersed particles; compared with the Kapitza radius R_K·k_matrix.' },
  // ---- meta
  { key: 'kA', label: 'Laminate layer A conductivity', unit: 'W/m/K', value: 400, min: 0.01, max: 3000, log: true, group: 'Laminate', showIf: { mode: ['meta'] }, help: 'High-conductivity layer of the laminate (copper: 400 W/m/K).' },
  { key: 'kBl', label: 'Laminate layer B conductivity', unit: 'W/m/K', value: 0.15, min: 0.005, max: 3000, log: true, group: 'Laminate', showIf: { mode: ['meta'] }, help: 'Low-conductivity layer (polymers: 0.15–0.3 W/m/K).' },
  { key: 'fA', label: 'Volume fraction of layer A', unit: '', value: 0.5, min: 0.01, max: 0.99, group: 'Laminate', showIf: { mode: ['meta'] }, help: 'Thickness fraction of material A in the periodic stack.' },
  { key: 'ang', label: 'Layer orientation angle', unit: '°', value: 30, min: 0, max: 90, group: 'Laminate', showIf: { mode: ['meta'] }, help: 'Angle between the layers and the x-axis (or the radial direction in a shell). It tilts the conductivity tensor and deflects or rotates the heat flux.' },
  { key: 'kbg', label: 'Background conductivity', unit: 'W/m/K', value: 2.3, min: 0.01, max: 500, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Conductivity of the host material in which the cloak is embedded.' },
  { key: 'kin', label: 'Inner (insulating) layer conductivity', unit: 'W/m/K', value: 0.03, min: 0.001, max: 100, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Conductivity of the inner shell of the bilayer cloak; ideally zero (expanded polystyrene: 0.03 W/m/K).' },
  { key: 'kout', label: 'Outer layer conductivity', unit: 'W/m/K', value: 9.8, min: 0.01, max: 3000, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Conductivity of the outer compensating shell; the design value is reported for comparison.' },
  { key: 'kcore', label: 'Cloaked object conductivity', unit: 'W/m/K', value: 200, min: 0.001, max: 3000, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Conductivity of the object inside the cloak.' },
  { key: 'ra', label: 'Object (core) radius a', unit: 'm', value: 0.006, min: 1e-8, max: 10, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Radius of the cloaked region and inner radius of the shell.' },
  { key: 'rb', label: 'Inner-layer outer radius b', unit: 'm', value: 0.0095, min: 1e-8, max: 10, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Interface between the insulating and conducting layers.' },
  { key: 'rc', label: 'Outer radius c', unit: 'm', value: 0.012, min: 1e-8, max: 10, log: true, group: 'Cloak', showIf: { mode: ['meta'] }, help: 'Outer radius of the cloak or shell.' },
  // ---- nte
  { key: 'aM', label: 'Matrix thermal expansion coefficient', unit: '1/K', value: 17e-6, min: -1e-4, max: 3e-4, group: 'Composite', bus: 'prop.cte', showIf: { mode: ['nte'] }, help: 'Linear CTE of the matrix (copper 17×10⁻⁶, aluminium 23×10⁻⁶, epoxy about 60×10⁻⁶ 1/K).' },
  { key: 'KM', label: 'Matrix bulk modulus', unit: 'GPa', value: 140, min: 0.1, max: 600, log: true, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Bulk modulus of the matrix (copper 140 GPa, epoxy about 5 GPa).' },
  { key: 'GM', label: 'Matrix shear modulus', unit: 'GPa', value: 48, min: 0.03, max: 300, log: true, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Shear modulus of the matrix; it enters the Hashin–Shtrikman bounds and the misfit stress.' },
  { key: 'aN', label: 'Filler thermal expansion coefficient', unit: '1/K', value: -9e-6, min: -1e-4, max: 3e-4, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Linear CTE of the filler; negative for NTE compounds such as ZrW₂O₈ (about −9×10⁻⁶ 1/K) or ScF₃.' },
  { key: 'KN', label: 'Filler bulk modulus', unit: 'GPa', value: 72, min: 0.1, max: 600, log: true, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Bulk modulus of the filler (about 70 GPa for ZrW₂O₈). A stiff filler imposes its expansion more effectively.' },
  { key: 'GN', label: 'Filler shear modulus', unit: 'GPa', value: 34, min: 0.03, max: 300, log: true, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Shear modulus of the filler.' },
  { key: 'Vf', label: 'Filler volume fraction', unit: '', value: 0.4, min: 0, max: 0.95, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Volume fraction of the NTE filler.' },
  { key: 'dTn', label: 'Temperature excursion', unit: 'K', value: 100, min: -1500, max: 1500, group: 'Composite', showIf: { mode: ['nte'] }, help: 'Temperature change for thermal strain and thermal-mismatch stress.' },
  { key: 'aHi', label: 'Lattice: high-CTE member (base)', unit: '1/K', value: 23e-6, min: 1e-7, max: 3e-4, group: 'Lattice', showIf: { mode: ['nte'] }, help: 'CTE of the base member of the bi-material triangle (aluminium: 23×10⁻⁶ 1/K).' },
  { key: 'aLo', label: 'Lattice: low-CTE member (legs)', unit: '1/K', value: 1.2e-6, min: 0, max: 3e-4, group: 'Lattice', showIf: { mode: ['nte'] }, help: 'CTE of the two inclined members (Invar: 1.2×10⁻⁶, titanium: 8.6×10⁻⁶ 1/K).' },
  { key: 'phiL', label: 'Lattice: half apex angle', unit: '°', value: 30, min: 1, max: 85, group: 'Lattice', showIf: { mode: ['nte'] }, help: 'Half the angle between the two legs. Wider triangles amplify the pull-in of the apex by the expanding base.' },
  // ---- rect
  { key: 'kA0', label: 'Segment A conductivity at T_ref', unit: 'W/m/K', value: 1, min: 1e-3, max: 3000, log: true, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Conductivity of segment A at the mean temperature.' },
  { key: 'nA', label: 'Segment A temperature exponent', unit: '', value: 1.5, min: -3, max: 3, group: 'Junction', showIf: { mode: ['rect'] }, help: 'k_A = k_A0 (T/T_ref)^n. Positive below the conductivity peak (boundary scattering), negative above it (Umklapp scattering).' },
  { key: 'LA', label: 'Segment A length', unit: 'm', value: 1e-6, min: 1e-9, max: 1, log: true, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Length of segment A along the heat-flow direction.' },
  { key: 'kB0', label: 'Segment B conductivity at T_ref', unit: 'W/m/K', value: 1, min: 1e-3, max: 3000, log: true, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Conductivity of segment B at the mean temperature.' },
  { key: 'nB', label: 'Segment B temperature exponent', unit: '', value: -1.5, min: -3, max: 3, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Opposite signs of the two exponents are the classic recipe for a bulk thermal diode.' },
  { key: 'LB', label: 'Segment B length', unit: 'm', value: 1e-6, min: 1e-9, max: 1, log: true, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Length of segment B.' },
  { key: 'Gi0', label: 'Interface conductance at perfect spectral overlap', unit: 'W/m²/K', value: 2e5, min: 1, max: 1e10, log: true, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Maximum conductance of the weak link between the segments. It must be small compared with the segment conductances k/L for interface effects to dominate.' },
  { key: 'del0', label: 'Spectral mismatch at zero bias', unit: 'bandwidths', value: 2, min: -5, max: 5, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Offset between the dominant phonon bands of the two sides, in units of the band width, when both are at the mean temperature.' },
  { key: 'sA', label: 'Band shift coefficient of side A', unit: '', value: 6, min: -30, max: 30, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Change of the mismatch per unit fractional temperature change of the A-side interface (anharmonic softening or hardening of its spectrum).' },
  { key: 'sB', label: 'Band shift coefficient of side B', unit: '', value: 6, min: -30, max: 30, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Same for the B side; it enters with the opposite sign.' },
  { key: 'Tm', label: 'Mean temperature', unit: 'K', value: 300, min: 20, max: 1500, group: 'Junction', showIf: { mode: ['rect'] }, help: 'The two terminals are held at T_mean ± ΔT/2.' },
  { key: 'dTmax', label: 'Bias sweep range ±ΔT', unit: 'K', value: 200, min: 1, max: 2000, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Largest temperature difference of the sweep (limited to keep both terminals above absolute zero).' },
  { key: 'dTop', label: 'Operating bias ΔT', unit: 'K', value: 100, min: 1, max: 2000, group: 'Junction', showIf: { mode: ['rect'] }, help: 'Temperature difference at which heat current and rectification ratio are reported.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
const num = (inp, k) => { const d = DEF[k]; let v = +inp[k]; if (!Number.isFinite(v)) v = d.value; return N.clamp(v, d.min, d.max); };
const sel = (inp, k) => (DEF[k].options.some(o => o.value === inp[k]) ? inp[k] : DEF[k].value);

// ---------------------------------------------------------------- transient slab conduction
// Conjugate gradient for the 5-point conduction matrix with a line (tridiagonal block) preconditioner taken along the
// more strongly coupled grid direction; the kernel cg offers only diagonal preconditioning, which stalls on thin slabs.
function pcgLine(nx, nz, diag, cx, cz, b, x0, tol, maxIter) {
  const n = nx * nz, x = Float64Array.from(x0), r = new Float64Array(n), z = new Float64Array(n), p = new Float64Array(n), Ap = new Float64Array(n);
  const alongZ = cz >= cx, m = alongZ ? nz : nx, lines = alongZ ? nx : nz, off = alongZ ? cz : cx, st = alongZ ? nx : 1, cp = new Float64Array(m), dp = new Float64Array(m);
  const mv = (v, o) => { for (let j = 0; j < nz; j++) for (let i = 0; i < nx; i++) { const q = j * nx + i; let s = diag[q] * v[q]; if (j > 0) s -= cz * v[q - nx]; if (j < nz - 1) s -= cz * v[q + nx]; if (i > 0) s -= cx * v[q - 1]; if (i < nx - 1) s -= cx * v[q + 1]; o[q] = s; } };
  const prec = () => { for (let l = 0; l < lines; l++) { const q0 = alongZ ? l : l * nx; cp[0] = -off / diag[q0]; dp[0] = r[q0] / diag[q0]; for (let k = 1; k < m; k++) { const q = q0 + k * st, mm = diag[q] + off * cp[k - 1]; cp[k] = -off / mm; dp[k] = (r[q] + off * dp[k - 1]) / mm; } z[q0 + (m - 1) * st] = dp[m - 1]; for (let k = m - 2; k >= 0; k--) z[q0 + k * st] = dp[k] - cp[k] * z[q0 + (k + 1) * st]; } };
  let bn = 0, rn = 0, rz = 0, it = 0; mv(x, Ap); for (let q = 0; q < n; q++) { r[q] = b[q] - Ap[q]; bn += b[q] * b[q]; rn += r[q] * r[q]; } bn = Math.sqrt(bn) || 1; prec(); for (let q = 0; q < n; q++) { p[q] = z[q]; rz += r[q] * z[q]; }
  let res = Math.sqrt(rn) / bn;
  while (res > tol && it < maxIter) { mv(p, Ap); let pAp = 0; for (let q = 0; q < n; q++) pAp += p[q] * Ap[q]; const al = rz / pAp; rn = 0; for (let q = 0; q < n; q++) { x[q] += al * p[q]; r[q] -= al * Ap[q]; rn += r[q] * r[q]; } prec(); let rz2 = 0; for (let q = 0; q < n; q++) rz2 += r[q] * z[q]; const be = rz2 / rz; rz = rz2; for (let q = 0; q < n; q++) p[q] = z[q] + be * p[q]; res = Math.sqrt(rn) / bn; it++; }
  return { x, iterations: it, residual: res };
}
function runSlab(inp) {
  const k = num(inp, 'k'), rho = num(inp, 'rho'), cp = num(inp, 'cp'), T0 = num(inp, 'T0'), h = num(inp, 'h'), Ta = num(inp, 'Tamb'), dim = sel(inp, 'dim'), Lz = num(inp, 'Lz'), Lx = num(inp, 'Lxw'), Ly = num(inp, 'Lyd');
  const E = num(inp, 'Ef'), s0 = num(inp, 'sig'), al = num(inp, 'alS'), qR0 = num(inp, 'qR'), EaR = num(inp, 'EaR') * 1e3, bcA = sel(inp, 'bcA'), bcB = sel(inp, 'bcB'), tEnd = num(inp, 'tEnd'), Tset = num(inp, 'Tset'), hf = num(inp, 'hfrac'), warnings = [], checks = [];
  const qJ = T => s0 * N.clamp(1 + al * (T - T0), 0.02, 50) * E * E, qRx = T => qR0 * Math.exp(N.clamp((EaR / R) * (1 / T0 - 1 / Math.max(T, 1)), -50, 50));
  const two = dim === '2d', nz = two ? 16 : 60, nx = two ? 40 : 1, n = nx * nz, dz = Lz / nz, dx = Lx / nx, steps = two ? 180 : 240, dt = tEnd / steps, TCAP = 4000;
  const gOf = (bc, half) => (bc === 'robin' ? 1 / (1 / h + half / k) : bc === 'fixed' ? k / half : 0);   // W/m²/K to ambient
  const GA = gOf(bcA, dz / 2), GB = gOf(bcB, dz / 2), GS = two ? 1 / (1 / h + dx / (2 * k)) : 0;
  const cz = k * dx / dz, cx = two ? k * dz / dx : 0, Cc = rho * cp * dx * dz, id = (i, j) => j * nx + i;
  const heated = i => !two || Math.abs((i + 0.5) * dx - Lx / 2) <= hf * Lx / 2 + 1e-12;
  let T = new Float64Array(n).fill(T0), Ein = 0, Eout = 0, runaway = false, cgIt = 0, cgRes = 0;
  const diag = new Float64Array(n), gb = new Float64Array(n);   // boundary conductance per cell (per unit depth)
  for (let j = 0; j < nz; j++) for (let i = 0; i < nx; i++) { const q = id(i, j); gb[q] = (j === 0 ? GA * dx : 0) + (j === nz - 1 ? GB * dx : 0) + (two && i === 0 ? GS * dz : 0) + (two && i === nx - 1 ? GS * dz : 0); diag[q] = Cc / dt + gb[q] + (j > 0 ? cz : 0) + (j < nz - 1 ? cz : 0) + (i > 0 ? cx : 0) + (i < nx - 1 ? cx : 0); }
  const ta = new Float64Array(nz), tb = new Float64Array(nz), tc = new Float64Array(nz);
  if (!two) for (let j = 0; j < nz; j++) { ta[j] = j > 0 ? -cz : 0; tc[j] = j < nz - 1 ? -cz : 0; tb[j] = diag[j]; }
  const hist = { t: [0], mean: [T0], max: [T0], min: [T0], ein: [0], eout: [0] }, snaps = [{ t: 0, T: Float64Array.from(T) }], snapAt = new Set([0.05, 0.15, 0.4, 1].map(f => Math.max(1, Math.round(f * steps))));
  const tmap = two ? null : new Float64Array(60 * nz), src = new Float64Array(n);
  for (let st = 1; st <= steps; st++) {
    const rhs = new Float64Array(n); let qs = 0;
    for (let j = 0; j < nz; j++) for (let i = 0; i < nx; i++) { const q = id(i, j); src[q] = ((heated(i) ? qJ(T[q]) : 0) + qRx(T[q])) * dx * dz; qs += src[q]; rhs[q] = Cc / dt * T[q] + src[q] + gb[q] * Ta; }
    let Tn;
    if (two) { const r = pcgLine(nx, nz, diag, cx, cz, rhs, T, 1e-11, 2000); Tn = r.x; cgIt += r.iterations; cgRes = Math.max(cgRes, r.residual); } else Tn = N.solveTridiag(ta, tb, tc, rhs);
    let lo = 0; for (let q = 0; q < n; q++) { if (!(Tn[q] < TCAP)) { Tn[q] = TCAP; runaway = true; } if (Tn[q] < 1) Tn[q] = 1; lo += gb[q] * (Tn[q] - Ta); }
    Ein += qs * dt * Ly; Eout += lo * dt * Ly; T = Tn;
    let m = 0, mx = -Infinity, mn = Infinity; for (let q = 0; q < n; q++) { m += T[q]; if (T[q] > mx) mx = T[q]; if (T[q] < mn) mn = T[q]; }
    hist.t.push(st * dt); hist.mean.push(m / n); hist.max.push(mx); hist.min.push(mn); hist.ein.push(Ein); hist.eout.push(Eout);
    if (snapAt.has(st)) snaps.push({ t: st * dt, T: Float64Array.from(T) });
    if (tmap) { const col = Math.round(st * 60 / steps) - 1; if (col >= 0 && st === Math.round((col + 1) * steps / 60)) for (let j = 0; j < nz; j++) tmap[j * 60 + col] = T[j]; }
  }
  let stored = 0; for (let q = 0; q < n; q++) stored += Cc * (T[q] - T0); stored *= Ly;
  const ebRes = Math.abs(stored - (Ein - Eout)) / Math.max(Math.abs(Ein), Math.abs(Eout), Math.abs(stored), 1e-12);
  const iHot = N.argmax(T), Tmax = T[iHot], Tmin = N.minOf(T), Tmean = hist.mean[steps], zHot = (Math.floor(iHot / nx) + 0.5) * dz, xHot = ((iHot % nx) + 0.5) * dx;
  let tSet = null; const up = Tset >= T0; for (let i = 1; i <= steps; i++) { const a = hist.mean[i - 1], b = hist.mean[i]; if ((up && a < Tset && b >= Tset) || (!up && a > Tset && b <= Tset)) { tSet = hist.t[i - 1] + (Tset - a) / (b - a) * dt; break; } } if (tSet === null && Math.abs(T0 - Tset) < 1e-9) tSet = 0;
  const Gsum = GA + GB, tau = Gsum > 0 ? rho * cp * Lz / Gsum : Infinity, Bi = h * (Lz / 2) / k, vol = Lx * Ly * Lz, qJ0 = qJ(T0), Pj0 = qJ0 * vol * (two ? hf : 1), lossEnd = (hist.eout[steps] - hist.eout[steps - 1]) / dt, powEnd = (hist.ein[steps] - hist.ein[steps - 1]) / dt;
  const symRobin = bcA === 'robin' && bcB === 'robin', TssLin = symRobin ? Ta + (qJ0 + qRx(T0)) * Lz / (2 * h) + (qJ0 + qRx(T0)) * Lz * Lz / (8 * k) : null;
  const gainLoop = Gsum > 0 ? (s0 * al * E * E + qR0 * EaR / (R * T0 * T0)) * Lz / Gsum : Infinity, steady = Math.abs(hist.mean[steps] - hist.mean[steps - 1]) / dt < 1e-3 * Math.max(Math.abs(Tmean - T0), 1e-6) / Math.max(tau, dt);
  const fo = k / (rho * cp) * dt / (dz * dz);
  checks.push(chk('Energy balance (stored = generated − lost)', ebRes < 1e-6, `relative residual ${F(ebRes, 2)}; stored ${F(stored, 4)} J, in ${F(Ein, 4)} J, out ${F(Eout, 4)} J`));
  checks.push(chk('No thermal runaway (temperature below the model cap)', !runaway, runaway ? `temperature reached the ${TCAP} K cap` : `maximum ${F(Tmax, 5)} K`));
  checks.push(chk('Positive-feedback loop gain below one (stable against electrothermal/reaction runaway)', gainLoop < 1, `(σ₀αE² + q_R E_a/RT₀²)·L/ΣG = ${F(Math.min(gainLoop, 1e300), 3)}`));
  if (!two && symRobin && al === 0 && qR0 === 0 && tEnd > 6 * tau) { const err = Math.abs(Tmax - TssLin) / Math.max(Math.abs(TssLin - Ta), 1e-9); checks.push(chk('Numerical steady state matches the analytical solution', err < 0.01, `T_max ${F(Tmax, 6)} K vs analytical ${F(TssLin, 6)} K`)); }
  checks.push(chk('Temporal resolution (at least 20 steps per time constant)', !Number.isFinite(tau) || dt < tau / 20 || steady, `Δt = ${F(dt, 3)} s, τ = ${Number.isFinite(tau) ? F(tau, 3) + ' s' : '∞'}; grid Fourier number ${F(fo, 3)} (implicit scheme, unconditionally stable)`));
  if (runaway) warnings.push('Thermal runaway: the temperature-dependent sources outgrow the heat removal and the solution was capped. Reduce the field, improve cooling, or use a current-limited supply.');
  else if (gainLoop >= 1) warnings.push('The linearised feedback gain of the temperature-dependent sources exceeds one: the steady state is unstable and the temperature will keep rising with time.');
  if (Tmax > 373 && rho > 900 && rho < 1100 && cp > 3800) warnings.push('Temperature exceeds 373 K in what appears to be a water-based medium: boiling (or pressurisation in a sealed vessel) is not modelled.');
  if (two) warnings.push('2-D cross-section model: the field is taken as uniform under the electrodes and zero elsewhere (no current spreading), and the third dimension is uniform.');
  warnings.push('Constant thermal properties and a uniform electric field; natural convection inside the liquid, which flattens gradients, is represented only through the conductivity value.');
  const zc = Array.from({ length: nz }, (_, j) => (j + 0.5) * dz * 1e3), xc = Array.from({ length: nx }, (_, i) => (i + 0.5) * dx);
  const ti = N.linspace(0, steps, Math.min(steps + 1, 241)).map(Math.round), th = ti.map(i => hist.t[i]);
  const series = [
    { key: 'Tt', title: 'Temperature vs time', xlabel: 'time (s)', ylabel: 'temperature (K)', lines: [{ name: 'mean', x: th, y: ti.map(i => hist.mean[i]) }, { name: 'maximum', x: th, y: ti.map(i => hist.max[i]) }, { name: 'minimum', x: th, y: ti.map(i => hist.min[i]) }, { name: 'setpoint', x: [0, tEnd], y: [Tset, Tset] }], explain: 'Heating curve of the slab. The gap between maximum and minimum is the internal non-uniformity; the approach to a plateau has the time constant ρc_pL/Σh.' },
    two ? { key: 'prof', title: 'Mid-plane temperature across the width', xlabel: 'x (m)', ylabel: 'temperature (K)', lines: snaps.map(s => ({ name: `t = ${F(s.t, 3)} s`, x: xc, y: xc.map((_, i) => s.T[id(i, Math.floor(nz / 2))]) })), explain: 'Temperature along the width at mid-thickness: hot under the electrodes, cooler in the unheated margins and at the edges.' }
      : { key: 'prof', title: 'Through-thickness temperature profiles', xlabel: 'z (mm)', ylabel: 'temperature (K)', lines: snaps.map(s => ({ name: `t = ${F(s.t, 3)} s`, x: zc, y: Array.from(s.T) })), explain: 'Profiles at several times. With volumetric heating and cooled faces the profile becomes parabolic with the hotspot in the interior.' },
    { key: 'E', title: 'Energy budget', xlabel: 'time (s)', ylabel: 'energy (J)', lines: [{ name: 'generated', x: th, y: ti.map(i => hist.ein[i]) }, { name: 'lost to surroundings', x: th, y: ti.map(i => hist.eout[i]) }, { name: 'stored', x: th, y: ti.map(i => hist.ein[i] - hist.eout[i]) }], explain: 'Cumulative Joule plus reaction heat, heat lost through the boundaries and their difference (stored sensible heat). At steady state generation and loss grow at the same rate.' },
  ];
  const maps = two ? [{ key: 'T', title: 'Temperature field at the end time', nx, ny: nz, data: Float64Array.from(T), x0: 0, x1: Lx, y0: 0, y1: Lz, xlabel: 'x (m)', ylabel: 'z (m)', unit: 'K', explain: 'Cross-section of the slab (thickness exaggerated). The hotspot sits under the electrodes at mid-thickness; edges and unheated margins are cooler.', diverging: false }]
    : [{ key: 'Tzt', title: 'Temperature through the thickness vs time', nx: 60, ny: nz, data: tmap, x0: 0, x1: tEnd, y0: 0, y1: Lz, xlabel: 'time (s)', ylabel: 'z (m)', unit: 'K', explain: 'Space–time map of the through-thickness temperature: heat builds up in the interior while the cooled faces lag behind.', diverging: false }];
  const ie = safe(() => chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }]), null);
  return {
    summary: `After ${F(tEnd, 4)} s the slab reaches a mean temperature of ${F(Tmean, 5)} K and a maximum of ${F(Tmax, 5)} K at z = ${F(zHot * 1e3, 3)} mm${two ? `, x = ${F(xHot, 3)} m` : ''}, with an initial Joule power density of ${F(qJ0, 3)} W/m³ (${F(Pj0, 3)} W in total). ${tSet === null ? `The setpoint of ${F(Tset, 4)} K is not reached` : `The setpoint of ${F(Tset, 4)} K is reached after ${F(tSet, 3)} s`}. The energy balance closes to ${F(ebRes, 2)} (relative) and the thermal time constant is ${Number.isFinite(tau) ? F(tau, 3) + ' s' : 'infinite (adiabatic)'}.`,
    kpis: [
      kpi('Tmax', 'Maximum temperature', Tmax, 'K', `Highest cell temperature at the end time (${F(Tmax - 273.15, 4)} °C); the value to compare with boiling, decomposition or material limits.`, !runaway),
      kpi('Tmean', 'Mean temperature', Tmean, 'K', 'Volume-averaged temperature at the end time, which governs the mean reaction rate; published as the process temperature.'),
      kpi('zhot', 'Hotspot position (thickness)', zHot * 1e3, 'mm', 'Distance of the hottest cell from the z = 0 face. Mid-thickness for symmetric cooling; it moves towards an insulated face.'),
      kpi('dT', 'Internal temperature non-uniformity', Tmax - Tmin, 'K', 'Maximum minus minimum temperature in the slab. With an activation energy of tens of kJ/mol even a few kelvin changes local reaction rates by tens of per cent.', Tmax - Tmin < 5),
      kpi('tset', 'Time to reach the setpoint', tSet === null ? 'not reached' : tSet, 's', 'First time at which the mean temperature crosses the setpoint (linear interpolation between time steps).'),
      kpi('qJ', 'Joule power density (initial)', qJ0, 'W/m³', 'σ₀E² at the initial temperature. It rises with temperature when the conductivity coefficient is positive.'),
      kpi('Pj', 'Total Joule power (initial)', Pj0, 'W', 'Power density times the heated volume; this is the electrical power the supply must deliver at the start.'),
      kpi('Pend', 'Heat generation at the end time', powEnd, 'W', `Joule plus reaction heat at the final temperature field; heat loss at the same instant is ${F(lossEnd, 4)} W (equal at steady state).`),
      kpi('tau', 'Thermal time constant', Number.isFinite(tau) ? tau : 'infinite (adiabatic)', 's', 'Lumped estimate ρc_pL/ΣG of the time to approach steady state; published as the heat-up time constant for the growth suite.'),
      kpi('Bi', 'Biot number', Bi, '', 'h(L/2)/k: below 0.1 the slab is nearly isothermal through its thickness; above 1 internal conduction limits heat removal.', Bi < 0.1 ? true : null),
      kpi('Tss', 'Analytical steady maximum (constant sources)', TssLin === null ? 'not applicable' : TssLin, 'K', 'T_amb + qL/2h + qL²/8k for uniform constant heating with both faces convectively cooled; the numerical value is higher when sources grow with temperature.'),
      kpi('eb', 'Energy balance residual', ebRes, '', 'Relative difference between stored energy and generated minus lost energy: a conservation check of the finite-volume solution.', ebRes < 1e-6),
    ],
    series, maps,
    equations: [
      eqn('governing', 'Transient heat conduction with sources', 'ρ c_p ∂T/∂t = ∇·(k ∇T) + σ(T) E² + q_R(T)'),
      eqn('constitutive', 'Joule heating', 'q_J = σ E² = J·E,   σ(T) = σ₀ [1 + α (T − T₀)]'),
      eqn('constitutive', 'Reaction heat (Arrhenius)', 'q_R(T) = q_R(T₀) exp[ (E_a/R)(1/T₀ − 1/T) ]'),
      eqn('governing', 'Boundary conditions', '−k ∂T/∂n = h (T − T_∞)  (Robin);   T = T_∞  (fixed);   ∂T/∂n = 0  (insulated)'),
      ...(ie ? [eqn('chemical', 'Example heat-releasing reaction (ZnO formation)', ie.molecular, 'ΔH ≈ −22 kJ per mole ZnO from standard formation enthalpies'), eqn('ionic', 'Net ionic equation', ie.netIonic)] : []),
    ],
    checks, warnings,
    recommendations: [
      rec('high', runaway || gainLoop >= 1 ? 'The operating point is electrothermally unstable. Use a constant-current (not constant-voltage) supply, lower the field, or raise the heat-transfer coefficient so that the feedback gain falls below one.' : tSet === null ? `The setpoint is not reached within ${F(tEnd, 3)} s (mean ${F(Tmean, 4)} K). ${Tset > Tmean ? `Raise the field: Joule power scales with E², so about ${F(E * Math.sqrt(Math.max((Tset - Ta) / Math.max(Tmean - Ta, 1e-9), 1)), 3)} V/m is needed, or insulate one face.` : 'Increase cooling or reduce the source.'}` : `The setpoint is reached in ${F(tSet, 3)} s; hold it with closed-loop control of the field since the steady temperature would otherwise drift to ${F(Tmax, 4)} K.`, 'Heating curve, setpoint crossing and feedback gain.'),
      rec('medium', Tmax - Tmin > 5 ? `Internal non-uniformity is ${F(Tmax - Tmin, 3)} K. Reduce the thickness (the conductive part scales with L²), pulse the field, or stir/recirculate to even out reaction rates.` : `Internal non-uniformity is only ${F(Tmax - Tmin, 3)} K: the slab is thermally thin (Bi = ${F(Bi, 3)}), so a lumped model is adequate for control design.`, 'Biot number and temperature spread.'),
      rec('medium', `Place the control thermocouple at the hotspot (z ≈ ${F(zHot * 1e3, 3)} mm) or correct a wall reading by the predicted offset of ${F(Tmax - T[0], 3)} K.`, 'Location of the maximum temperature.'),
      rec('low', `At the end of the run ${F(lossEnd, 3)} W is being lost to the surroundings. Insulating the faces to a tenth of the present heat-transfer coefficient would hold the same temperature with roughly a tenth of the electrical power (${F(Ein / 3.6e6, 3)} kWh were used in this run).`, 'Steady-state balance between generation and loss.', true),
    ],
    publish: { 'therm.T_max': Tmax, 'proc.T_C': Tmean - 273.15, 'proc.power': powEnd, 'therm.tau_heat': Number.isFinite(tau) ? tau : 1e30, 'therm.dT_internal': Tmax - Tmin, 'elec.joule_W_m3': qJ0 },
    numerics: { method: two ? 'Cell-centred finite volumes, backward Euler in time, conjugate gradient with a line (tridiagonal) preconditioner each step; temperature-dependent sources lagged one step' : 'Cell-centred finite volumes, backward Euler in time, tridiagonal (Thomas) solve each step; temperature-dependent sources lagged one step', dof: n, steps, iterations: two ? cgIt : steps, residual: two ? cgRes : ebRes, notes: `Grid ${nx} × ${nz}; time step ${F(dt, 3)} s.` },
  };
}

// ---------------------------------------------------------------- lumped / wall / fin / schedule
function runLumped(inp) {
  const k = num(inp, 'k'), rho = num(inp, 'rho'), cp = num(inp, 'cp'), T0 = num(inp, 'T0'), h = num(inp, 'h'), Ta = num(inp, 'Tamb'), V = num(inp, 'Vl'), A = num(inp, 'Al'), P = num(inp, 'Ph'), Ts = num(inp, 'Tsoak'), tS = num(inp, 'tSoak'), Tc = num(inp, 'Tcool');
  const Lw = num(inp, 'Lw'), kw = num(inp, 'kw'), Li = num(inp, 'Lins'), ki = num(inp, 'kins'), df = num(inp, 'dfin'), Lf = num(inp, 'Lfin'), kf = num(inp, 'kfin'), warnings = [], checks = [];
  const Rw = Lw / kw, Ri = Li / ki, Rh = 1 / h, U = 1 / (Rw + Ri + Rh), UA = U * A, C = rho * V * cp, tau = C / UA, Lc = V / A, Bi = U * Lc / k, Tss = Ta + P / UA;
  const heatOK = Ts > T0 ? Tss > Ts : true, tRamp = Ts <= T0 ? 0 : heatOK ? tau * Math.log((Tss - T0) / (Tss - Ts)) : null, Phold = Math.max(UA * (Ts - Ta), 0);
  const coolOK = Tc > Ta && Ts > Tc, tCool = Ts <= Tc ? 0 : coolOK ? tau * Math.log((Ts - Ta) / (Tc - Ta)) : null;
  const Eramp = tRamp === null ? 0 : P * tRamp, Esoak = Phold * tS, Etot = Eramp + Esoak, Euse = Math.max(C * (Ts - T0), 0), eff = Etot > 0 ? Euse / Etot : 0, Pmin = UA * (Ts - Ta);
  const qw = U * (Ts - Ta), Tsurf = Ta + qw * Rh, Tint = Ts - qw * Rw;
  const m = Math.sqrt(4 * h / (kf * df)), mL = m * Lf, etaF = mL > 1e-9 ? Math.tanh(mL) / mL : 1, Ac = PI * df * df / 4, qFin = Math.sqrt(h * PI * df * kf * Ac) * (Ts - Ta) * Math.tanh(mL), epsF = qFin / Math.max(h * Ac * (Ts - Ta), 1e-300), BiF = h * df / (2 * kf);
  checks.push(chk('Lumped-capacitance validity (Bi < 0.1)', Bi < 0.1, `Bi = U(V/A)/k = ${F(Bi, 3)}`));
  checks.push(chk('Heater can reach the soak temperature', heatOK, `steady temperature at full power ${F(Tss, 5)} K vs soak ${F(Ts, 5)} K`));
  checks.push(chk('Cool-down target above ambient', Ts <= Tc || coolOK, `target ${F(Tc, 4)} K, ambient ${F(Ta, 4)} K`));
  checks.push(chk('One-dimensional fin assumption (h d / 2k_fin < 0.1)', BiF < 0.1, `fin Biot number ${F(BiF, 3)}`));
  checks.push(chk('Outer surface safe to touch (< 333 K)', Tsurf < 333, `outer surface ${F(Tsurf, 4)} K at soak`));
  if (Bi >= 0.1) warnings.push(`Biot number ${F(Bi, 3)} ≥ 0.1: the body is not isothermal and the lumped model underestimates the time for its centre to follow; use the transient conduction mode. (For a stirred liquid the effective conductivity is much higher and the lumped model is acceptable.)`);
  if (!heatOK) warnings.push(`The heater cannot reach the soak temperature: at ${F(P, 3)} W the vessel settles at ${F(Tss, 4)} K. At least ${F(Pmin, 3)} W is needed just to hold ${F(Ts, 4)} K.`);
  if (Ts > Tc && !coolOK) warnings.push('The cool-down target is at or below ambient and cannot be reached by natural cooling.');
  warnings.push('Uniform heat-loss coefficient over the stated area, constant properties, wall and insulation heat capacity neglected, and radiation included only through the value of h.');
  const tr = tRamp === null ? 3 * tau : tRamp, tc2 = tCool === null ? 3 * tau : tCool, x1 = N.linspace(0, tr, 60), x3 = N.linspace(0, tc2, 60);
  const sched = { name: 'ramp–soak–cool', x: [...x1, tr + tS, ...x3.slice(1).map(t => tr + tS + t)], y: [...x1.map(t => Tss + (T0 - Tss) * Math.exp(-t / tau)), tRamp === null ? Tss : Ts, ...x3.slice(1).map(t => Ta + ((tRamp === null ? Tss : Ts) - Ta) * Math.exp(-t / tau))] };
  const xf = N.linspace(0, Lf, 50);
  return {
    summary: `The ${F(V * 1e3, 3)} L body has an overall loss coefficient of ${F(U, 3)} W/m²/K (UA = ${F(UA, 3)} W/K), a time constant of ${F(tau, 3)} s and a Biot number of ${F(Bi, 3)}. ${tRamp === null ? `The heater cannot reach ${F(Ts, 4)} K` : `Heating to ${F(Ts, 4)} K takes ${F(tRamp, 3)} s`}, holding needs ${F(Phold, 3)} W, and natural cooling to ${F(Tc, 4)} K takes ${tCool === null ? 'indefinitely long' : F(tCool, 3) + ' s'}. The cycle uses ${F(Etot / 3.6e6, 3)} kWh, of which ${F(100 * eff, 3)} % is sensible heat of the contents.`,
    kpis: [
      kpi('Bi', 'Biot number', Bi, '', 'U·(V/A)/k: the lumped-capacitance model is valid below 0.1, where internal gradients are under about 5 % of the overall temperature difference.', Bi < 0.1),
      kpi('tau', 'Thermal time constant', tau, 's', 'ρVc_p/(UA): time to cover 63 % of a step change; heating and cooling take three to five time constants.'),
      kpi('U', 'Overall heat-loss coefficient', U, 'W/m²/K', `Series resistances of wall (${F(Rw, 2)}), insulation (${F(Ri, 2)}) and external film (${F(Rh, 2)} m²K/W).`),
      kpi('tRamp', 'Ramp time to soak temperature', tRamp === null ? 'not reachable' : tRamp, 's', 'τ ln[(T_ss − T₀)/(T_ss − T_soak)] at constant heater power, where T_ss is the steady temperature at full power.'),
      kpi('Phold', 'Holding power at soak', Phold, 'W', 'UA(T_soak − T_amb): the power that just balances the losses during the soak.'),
      kpi('tCool', 'Natural cool-down time', tCool === null ? 'not reachable' : tCool, 's', 'τ ln[(T_soak − T_amb)/(T_cool − T_amb)] with the heater off.'),
      kpi('E', 'Cycle energy (ramp + soak)', Etot / 3.6e6, 'kWh', `Ramp ${F(Eramp / 3.6e6, 3)} kWh plus soak ${F(Esoak / 3.6e6, 3)} kWh of heater energy.`),
      kpi('eff', 'Thermal efficiency of the cycle', eff, '', 'Sensible heat stored in the contents divided by the heater energy; long soaks in poorly insulated vessels give low values.', eff > 0.5 ? true : eff < 0.2 ? false : null),
      kpi('Tsurf', 'Outer surface temperature at soak', Tsurf, 'K', `From the composite-wall resistance network; the wall–insulation interface is at ${F(Tint, 4)} K. Above about 333 K (60 °C) surfaces cause burns.`, Tsurf < 333),
      kpi('qw', 'Wall heat flux at soak', qw, 'W/m²', 'U(T_soak − T_amb): heat loss per unit wall area during the hold.'),
      kpi('etaF', 'Fin efficiency', etaF, '', `tanh(mL)/(mL) with mL = ${F(mL, 3)} for a pin fin with an adiabatic tip: the fraction of the ideal (isothermal-fin) heat flow actually delivered.`),
      kpi('qFin', 'Heat flow through one fin', qFin, 'W', `√(hPkA_c)·θ_b·tanh(mL) at the soak-to-ambient temperature difference; fin effectiveness ${F(epsF, 3)} (heat flow relative to the bare base area).`),
    ],
    series: [
      { key: 'sched', title: 'Ramp–soak–cool temperature schedule', xlabel: 'time (s)', ylabel: 'temperature (K)', lines: [sched, { name: 'ambient', x: [0, sched.x[sched.x.length - 1]], y: [Ta, Ta] }], explain: 'Exponential approach to the soak temperature at constant power, the controlled hold and the natural cool-down. If the heater is undersized the ramp levels off below the soak temperature.' },
      { key: 'wall', title: 'Temperature through the composite wall', xlabel: 'distance from the inner surface (mm)', ylabel: 'temperature (K)', lines: [{ name: 'wall + insulation', x: [0, Lw * 1e3, (Lw + Li) * 1e3], y: [Ts, Tint, Tsurf] }, { name: 'ambient', x: [0, (Lw + Li) * 1e3 + 1e-9], y: [Ta, Ta] }], explain: 'Piecewise-linear steady profile: most of the temperature drop occurs across the layer with the largest thermal resistance.' },
      { key: 'fin', title: 'Temperature along the fin', xlabel: 'distance from the base (mm)', ylabel: 'temperature (K)', lines: [{ name: 'fin', x: xf.map(x => x * 1e3), y: xf.map(x => Ta + (Ts - Ta) * Math.cosh(Math.min(m * (Lf - x), 700)) / Math.cosh(Math.min(mL, 700))) }], explain: 'Adiabatic-tip fin profile cosh[m(L−x)]/cosh(mL). A fin much longer than 1/m adds material without adding heat flow.' },
    ],
    equations: [
      eqn('governing', 'Lumped energy balance', 'ρ V c_p dT/dt = P − U A (T − T_∞)   ⇒   T(t) = T_ss + (T₀ − T_ss) e^(−t/τ),  T_ss = T_∞ + P/UA,  τ = ρVc_p/UA'),
      eqn('constitutive', 'Composite wall (series resistances)', '1/U = L_w/k_w + L_ins/k_ins + 1/h'),
      eqn('formula', 'Biot number', 'Bi = U L_c / k,  L_c = V/A'),
      eqn('governing', 'Pin fin with adiabatic tip', 'd²θ/dx² = m² θ,  m = (4h / k d)^½;   q = (h P k A_c)^½ θ_b tanh(mL);   η = tanh(mL)/(mL)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', !heatOK ? `Raise the heater power above ${F(Pmin, 3)} W or add insulation: the present heater saturates at ${F(Tss, 4)} K.` : `Heater sizing is adequate (steady temperature at full power ${F(Tss, 4)} K). Use full power for the ramp and switch to about ${F(Phold, 3)} W (${F(100 * Phold / Math.max(P, 1e-9), 3)} % duty) for the hold to avoid overshoot.`, 'Lumped energy balance at full power.'),
      rec('medium', Tsurf >= 333 ? `Outer surface reaches ${F(Tsurf, 4)} K. Increase the insulation to at least ${F(ki * ((Ts - Ta) / Math.max(h * (333 - Ta), 1e-9) - Rw - Rh), 3)} m, or guard the surface.` : `Outer surface stays at ${F(Tsurf, 4)} K, safe to touch.`, 'Composite-wall temperature profile.'),
      rec('medium', mL > 3 ? `The fin is longer than useful (mL = ${F(mL, 3)}): beyond about ${F(2.65 / m * 1e3, 3)} mm it adds under 1 % more heat flow.` : `Fin efficiency is ${F(100 * etaF, 3)} %; lengthening it up to about ${F(2.65 / m * 1e3, 3)} mm still adds heat flow.`, 'Fin parameter mL.'),
      rec('low', `Losses during the soak account for ${F(100 * Esoak / Math.max(Etot, 1e-9), 3)} % of the cycle energy. Doubling the insulation thickness would cut the holding power to about ${F(A * (Ts - Ta) / (Rw + 2 * Ri + Rh), 3)} W; shortening the soak to the kinetically necessary time saves in proportion.`, 'Split of cycle energy between ramp and soak.', true),
    ],
    publish: { 'therm.tau_heat': tau, 'proc.power': Phold, 'fab.energy_kWh_batch': Etot / 3.6e6, 'ctrl.settling_s': tRamp === null ? 4 * tau : tRamp, 'proc.T_C': Ts - 273.15, 'therm.UA': UA },
    numerics: { method: 'Closed-form lumped-capacitance, series-resistance and fin solutions', dof: 1, steps: sched.x.length, iterations: 1, residual: 0, notes: 'Wall and insulation thermal mass neglected.' },
  };
}

// ---------------------------------------------------------------- furnace / wafer radiation balance
function runFurnace(inp) {
  const Th = num(inp, 'Th'), eW = num(inp, 'epsW'), eH = num(inp, 'epsH'), Tw = num(inp, 'Twall'), hg = num(inp, 'hg'), Tg = num(inp, 'Tg'), tw = num(inp, 'tw'), rc = num(inp, 'rcW'), kW = num(inp, 'kWf'), warnings = [], checks = [];
  const eEff = 1 / (1 / eW + 1 / eH - 1), qIn = T => eEff * SB * (Th ** 4 - T ** 4), qBack = T => eW * SB * (T ** 4 - Tw ** 4), qConv = T => 2 * hg * (T - Tg), net = T => qIn(T) - qBack(T) - qConv(T);
  const lo = Math.min(Th, Tw, Tg), hi = Math.max(Th, Tw, Tg), Tss = hi - lo < 1e-9 ? lo : safe(() => N.brent(net, lo, hi, 1e-9), 0.5 * (lo + hi));
  const hr = 4 * SB * Tss ** 3 * (eEff + eW), htot = hr + 2 * hg, tau = rc * tw / Math.max(htot, 1e-12), Ti = Math.min(300, lo), rate0 = net(Ti) / (rc * tw), Bi = htot / 2 * (tw / 2) / kW;
  const tEnd = 8 * tau, sol = N.rk45((_, y) => [net(y[0]) / (rc * tw)], [Ti], 0, tEnd, { rtol: 1e-7, atol: 1e-6, hmax: tEnd / 100, maxSteps: 4000 });
  let t99 = null; const band = 0.01 * Math.abs(Tss - Ti) + 1e-9; for (let i = 0; i < sol.t.length; i++) if (Math.abs(sol.y[i][0] - Tss) <= band) { t99 = sol.t[i]; break; }
  const dTde = (() => { const e2 = Math.min(eW + 0.01, 1), ee = 1 / (1 / e2 + 1 / eH - 1), f = T => ee * SB * (Th ** 4 - T ** 4) - e2 * SB * (T ** 4 - Tw ** 4) - 2 * hg * (T - Tg); return safe(() => N.brent(f, lo, hi, 1e-9), Tss) - Tss; })();
  const TssRad = Math.pow((eEff * Th ** 4 + eW * Tw ** 4) / (eEff + eW), 0.25), resid = Math.abs(net(Tss)) / Math.max(Math.abs(qIn(Tss)), 1e-9);
  checks.push(chk('Steady heat balance closed', resid < 1e-6, `net flux ${F(net(Tss), 2)} W/m² of ${F(qIn(Tss), 4)} W/m² absorbed`));
  checks.push(chk('Wafer isothermal through its thickness (Bi < 0.1)', Bi < 0.1, `Bi = ${F(Bi, 3)}`));
  checks.push(chk('Steady temperature between the surrounding temperatures (second law)', Tss >= lo - 1e-6 && Tss <= hi + 1e-6, `${F(lo, 4)} K ≤ ${F(Tss, 5)} K ≤ ${F(hi, 4)} K`));
  checks.push(chk('Transient integration reached steady state', t99 !== null, t99 === null ? 'not within 8 time constants' : `within 1 % after ${F(t99, 3)} s`));
  if (Bi >= 0.1) warnings.push('The sample is not thermally thin: front and back surfaces differ in temperature and a through-thickness conduction model is needed.');
  warnings.push('Grey, diffuse surfaces with temperature-independent emissivities; infinite parallel plates for the heater side and a large cold enclosure for the back side. Real wafers are spectrally selective and partly transparent below about 900 K.');
  const Ts = N.linspace(lo, hi, 80), ix = N.linspace(0, sol.t.length - 1, Math.min(sol.t.length, 200)).map(Math.round);
  return {
    summary: `With the heater at ${F(Th, 4)} K (effective exchange emissivity ${F(eEff, 3)}) the wafer settles at ${F(Tss, 5)} K (${F(Tss - 273.15, 4)} °C), absorbing ${F(qIn(Tss) / 1e3, 4)} kW/m² and losing ${F(qBack(Tss) / 1e3, 4)} kW/m² by radiation and ${F(qConv(Tss) / 1e3, 4)} kW/m² by convection. The thermal time constant is ${F(tau, 3)} s, with an initial heating rate of ${F(rate0, 3)} K/s. A 0.01 change of wafer emissivity shifts the temperature by ${F(dTde, 3)} K.`,
    kpis: [
      kpi('Tss', 'Steady wafer temperature', Tss, 'K', 'Root of the heat balance between radiation from the heater, radiation to the cold wall and convection to the gas (Brent\'s method).'),
      kpi('qin', 'Net radiative flux from the heater', qIn(Tss) / 1e3, 'kW/m²', 'ε_eff σ(T_h⁴ − T⁴) between two grey parallel plates with ε_eff = 1/(1/ε_w + 1/ε_h − 1).'),
      kpi('qback', 'Radiative loss to the chamber wall', qBack(Tss) / 1e3, 'kW/m²', 'ε_w σ(T⁴ − T_wall⁴) from the back face to a large cold enclosure.'),
      kpi('qconv', 'Convective loss to the gas', qConv(Tss) / 1e3, 'kW/m²', '2h(T − T_gas) from both faces. Negative values mean the gas heats the wafer.'),
      kpi('hr', 'Linearised radiative coefficient', hr, 'W/m²/K', '4σT³(ε_eff + ε_w): the equivalent heat-transfer coefficient of radiation at the steady temperature, to compare with convection.'),
      kpi('tau', 'Thermal time constant', tau, 's', 'ρc_p·t/(h_rad + 2h): response time of the wafer to heater changes near the steady state.'),
      kpi('t99', 'Time to within 1 % of steady state', t99 === null ? 'not reached' : t99, 's', 'From integrating the nonlinear lumped energy balance from the initial temperature.'),
      kpi('rate0', 'Initial heating rate', rate0, 'K/s', 'Net absorbed flux divided by the areal heat capacity at the start; rapid thermal processing reaches 50–200 K/s.'),
      kpi('dTde', 'Temperature shift per +0.01 emissivity', dTde, 'K', 'Sensitivity of the steady temperature to wafer emissivity: patterned films or metasurfaces with different emissivity run at different temperatures.'),
      kpi('TssRad', 'Radiation-only steady temperature', TssRad, 'K', 'Balance without convection, [(ε_eff T_h⁴ + ε_w T_wall⁴)/(ε_eff + ε_w)]^¼; the difference from the full result shows the role of the gas.'),
      kpi('Bi', 'Wafer Biot number', Bi, '', 'Mean surface coefficient times half-thickness over conductivity; below 0.1 the wafer is isothermal through its thickness.', Bi < 0.1),
    ],
    series: [
      { key: 'Tt', title: 'Wafer heat-up', xlabel: 'time (s)', ylabel: 'temperature (K)', lines: [{ name: 'wafer', x: ix.map(i => sol.t[i]), y: ix.map(i => sol.y[i][0]) }, { name: 'steady state', x: [0, tEnd], y: [Tss, Tss] }], explain: 'Nonlinear heat-up from the initial temperature: fast at first when the radiative driving force T_h⁴ − T⁴ is large, slowing as the losses grow.' },
      { key: 'bal', title: 'Heat-balance terms vs wafer temperature', xlabel: 'wafer temperature (K)', ylabel: 'flux (kW/m²)', lines: [{ name: 'from heater', x: Ts, y: Ts.map(T => qIn(T) / 1e3) }, { name: 'radiative loss', x: Ts, y: Ts.map(T => qBack(T) / 1e3) }, { name: 'convective loss', x: Ts, y: Ts.map(T => qConv(T) / 1e3) }, { name: 'net', x: Ts, y: Ts.map(T => net(T) / 1e3) }], explain: 'The steady temperature is where the net curve crosses zero. Its slope there is the restoring coefficient that sets the time constant.' },
    ],
    equations: [
      eqn('governing', 'Wafer energy balance', 'ρ c_p t dT/dt = ε_eff σ (T_h⁴ − T⁴) − ε_w σ (T⁴ − T_wall⁴) − 2 h (T − T_gas)'),
      eqn('constitutive', 'Stefan–Boltzmann exchange between grey parallel plates', 'q = σ (T_h⁴ − T⁴) / (1/ε_w + 1/ε_h − 1)'),
      eqn('formula', 'Linearised radiative coefficient', 'h_rad = 4 σ T³ (ε_eff + ε_w)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', `Temperature uniformity is set by emissivity: ${F(dTde, 3)} K per 0.01 of emissivity. Calibrate pyrometers for the actual film stack and expect patterned or metasurface-coated regions to run ${F(Math.abs(dTde) * 10, 3)} K apart for an emissivity contrast of 0.1.`, 'Emissivity sensitivity of the radiative balance.'),
      rec('medium', hr > 5 * 2 * hg ? `Radiation dominates (h_rad = ${F(hr, 3)} against ${F(2 * hg, 3)} W/m²/K for convection): gas flow changes have little thermal effect, so control the heater and the view factors.` : `Convection is significant (${F(2 * hg, 3)} against h_rad = ${F(hr, 3)} W/m²/K): gas flow rate and composition shift the wafer temperature and must be held constant.`, 'Comparison of radiative and convective coefficients.'),
      rec('medium', `The wafer follows heater changes with a time constant of ${F(tau, 3)} s; ramp-rate control faster than this needs feed-forward lamp power rather than feedback alone.`, 'Lumped time constant at the operating point.'),
      rec('low', `A cold-wall design loses ${F(qBack(Tss) / 1e3, 3)} kW/m² from the back of the wafer to the chamber. A reflective (low-emissivity) shield behind the wafer returns most of this and reduces heater power for the same temperature.`, 'Radiative loss term of the balance.', true),
    ],
    publish: { 'therm.T_max': Tss, 'proc.T_C': Tss - 273.15, 'therm.tau_heat': tau, 'therm.heat_flux': qIn(Tss) },
    numerics: { method: 'Brent root finding on the steady balance; adaptive Dormand–Prince RK45 for the heat-up transient', dof: 1, steps: sol.steps, iterations: sol.steps + sol.rejected, residual: resid, notes: 'Lumped (thermally thin) wafer.' },
  };
}

// ---------------------------------------------------------------- nanoscale conduction
function fuchs(delta, p) { const m = 200; let s = 0; for (let i = 1; i <= m; i++) { const u = i / m, w = i === m ? 1 : i % 2 ? 4 : 2, e = Math.exp(-Math.min(delta / u, 700)); s += w * (u - u ** 3) * (1 - e) / (1 - p * e); } return N.clamp(1 - 3 * (1 - p) / (2 * delta) * s / (3 * m), 0, 1); }
function runNano(inp) {
  const kb = num(inp, 'kb'), Cv = num(inp, 'Cv'), v = num(inp, 'vs'), d = num(inp, 'dch'), p = num(inp, 'spec'), RK = num(inp, 'RK'), kp = num(inp, 'kpar'), km = num(inp, 'kmat'), phi = num(inp, 'phiP'), a = num(inp, 'aP'), warnings = [], checks = [];
  const Lam = 3 * kb / (Cv * v), Kn = Lam / d, wire = x => kb / (1 + Lam / (x * (1 + p) / (1 - p))), film = x => kb * fuchs(x / Lam, p), cross = x => kb / (1 + 4 * Lam / (3 * x));
  const kWire = wire(d), kFilm = film(d), kCross = cross(d), regime = Kn < 0.1 ? 'diffusive (Fourier)' : Kn < 10 ? 'quasi-ballistic (transition)' : 'ballistic';
  const maxw = (f, kpp) => km * (kpp + 2 * km + 2 * f * (kpp - km)) / (kpp + 2 * km - f * (kpp - km)), hj = (f, rad_) => { const al = RK * km / rad_, n1 = kp * (1 + 2 * al) + 2 * km, n2 = kp * (1 - al) - km; return km * (n1 + 2 * f * n2) / (n1 - f * n2); };
  const kMax = maxw(phi, kp), kHJ = hj(phi, a), aK = RK * km, alpha = aK / a, aCrit = kp > km ? aK / (1 - km / kp) : Infinity, G = 1 / RK;
  const wLo = 1 / ((1 - phi) / km + phi / kp), wHi = (1 - phi) * km + phi * kp;
  checks.push(chk('Fuchs–Sondheimer factor recovers the bulk limit for thick films', Math.abs(fuchs(1000, p) - (1 - 3 * (1 - p) / 8000)) < 1e-4, `k/k_bulk at d = 1000Λ: ${F(fuchs(1000, p), 6)}`));
  checks.push(chk('Size-reduced conductivities do not exceed the bulk value', kWire <= kb && kFilm <= kb && kCross <= kb, `wire ${F(kWire, 4)}, film ${F(kFilm, 4)}, cross-plane ${F(kCross, 4)} W/m/K`));
  checks.push(chk('Composite conductivity within the Wiener (series/parallel) bounds', kMax >= wLo * (1 - 1e-9) && kMax <= wHi * (1 + 1e-9), `${F(wLo, 4)} ≤ ${F(kMax, 4)} ≤ ${F(wHi, 4)} W/m/K (Maxwell)`));
  checks.push(chk('Interface resistance never increases the composite conductivity', kHJ <= kMax * (1 + 1e-12), `Hasselman–Johnson ${F(kHJ, 4)} vs Maxwell ${F(kMax, 4)} W/m/K`));
  checks.push(chk('Dilute-dispersion validity (φ ≤ 0.25)', phi <= 0.25, `φ = ${F(phi, 3)}`));
  if (Kn > 10) warnings.push('Ballistic regime (Kn > 10): temperature is not well defined inside the structure and Fourier\'s law with a reduced conductivity is only an effective description; the phonon Boltzmann equation is needed for detail.');
  if (phi > 0.25) warnings.push('Above about 25 vol% particle interactions and percolation are significant; Maxwell-type effective-medium formulae underestimate the conductivity of chained or aggregated particles.');
  if (d < 5 * 5e-10 * 2) warnings.push('Below about 5 nm phonon confinement changes the dispersion relation itself; the bulk heat capacity and velocity are no longer appropriate.');
  warnings.push('Grey (single mean-free-path) kinetic model: real phonon spectra span mean free paths over orders of magnitude, so the grey model underestimates the size effect at large sizes and overestimates it at small sizes.');
  const ds = N.logspace(1e-9, 1e-4, 90), fs = N.linspace(0, 0.3, 60);
  return {
    summary: `The grey mean free path is ${F(Lam * 1e9, 3)} nm, so a ${F(d * 1e9, 3)} nm structure has a Knudsen number of ${F(Kn, 3)}: ${regime}. The conductivity falls from ${F(kb, 4)} W/m/K in bulk to ${F(kWire, 4)} (nanowire), ${F(kFilm, 4)} (film, in-plane) and ${F(kCross, 4)} W/m/K (film, cross-plane). For ${F(100 * phi, 3)} vol% of ${F(a * 1e9, 3)} nm-radius particles the effective conductivity is ${F(kHJ, 4)} W/m/K with interface resistance against ${F(kMax, 4)} W/m/K without (Kapitza radius ${F(aK * 1e9, 3)} nm).`,
    kpis: [
      kpi('Lam', 'Carrier mean free path (grey)', Lam * 1e9, 'nm', 'Λ = 3k/(C·v) from kinetic theory; the length scale below which boundaries limit conduction.'),
      kpi('Kn', 'Knudsen number', Kn, '', `Λ/d. ${regime[0].toUpperCase() + regime.slice(1)} transport; Fourier's law with bulk conductivity holds only for Kn below about 0.1.`, Kn < 0.1 ? true : null),
      kpi('kWire', 'Nanowire conductivity', kWire, 'W/m/K', `Matthiessen's rule with a Casimir boundary mean free path d(1+p)/(1−p); ${F(100 * kWire / kb, 3)} % of bulk.`),
      kpi('kFilm', 'Thin-film in-plane conductivity', kFilm, 'W/m/K', `Fuchs–Sondheimer integral for a film of thickness d with specularity p; ${F(100 * kFilm / kb, 3)} % of bulk.`),
      kpi('kCross', 'Thin-film cross-plane conductivity', kCross, 'W/m/K', `k/(1 + 4Kn/3), the ballistic–diffusive interpolation between Fourier conduction and the radiative (Casimir) limit; ${F(100 * kCross / kb, 3)} % of bulk.`),
      kpi('kHJ', 'Nanocomposite conductivity with interface resistance', kHJ, 'W/m/K', `Hasselman–Johnson model with α = R_K k_m/a = ${F(alpha, 3)}; ${F(100 * (kHJ / km - 1), 3)} % change relative to the matrix.`, kHJ > km),
      kpi('kMax', 'Nanocomposite conductivity, perfect interfaces', kMax, 'W/m/K', 'Maxwell (Maxwell-Garnett) effective-medium result for dilute spheres without interface resistance: the upper limit for well-dispersed particles.'),
      kpi('aK', 'Kapitza radius', aK * 1e9, 'nm', 'R_K·k_matrix: particles of about this radius or smaller are thermally masked by their interface.'),
      kpi('aCrit', 'Critical particle radius for enhancement', Number.isFinite(aCrit) ? aCrit * 1e9 : 'none (particle less conductive than matrix)', 'nm', 'Radius a_K/(1 − k_m/k_p) below which adding particles lowers the composite conductivity.'),
      kpi('G', 'Interface conductance', G / 1e6, 'MW/m²/K', 'Reciprocal of the Kapitza resistance; solid–solid interfaces span roughly 10–1000 MW/m²/K.'),
    ],
    series: [
      { key: 'kd', title: 'Conductivity vs characteristic size', xlabel: 'size (nm)', ylabel: 'k / k_bulk', logx: true, lines: [{ name: 'nanowire (diameter)', x: ds.map(x => x * 1e9), y: ds.map(x => wire(x) / kb) }, { name: 'film, in-plane (thickness)', x: ds.map(x => x * 1e9), y: ds.map(x => film(x) / kb) }, { name: 'film, cross-plane', x: ds.map(x => x * 1e9), y: ds.map(x => cross(x) / kb) }], explain: 'Classical size effect: conductivity falls once the size approaches the mean free path. Wires are affected more than films because carriers are confined in two directions.' },
      { key: 'kphi', title: 'Nanocomposite conductivity vs particle loading', xlabel: 'particle volume fraction', ylabel: 'k_eff / k_matrix', lines: [{ name: 'Maxwell (no interface resistance)', x: fs, y: fs.map(f => maxw(f, kp) / km) }, { name: `Hasselman–Johnson, a = ${F(a * 1e9, 3)} nm`, x: fs, y: fs.map(f => hj(f, a) / km) }, { name: `Hasselman–Johnson, a = ${F(a * 1e9 / 5, 3)} nm`, x: fs, y: fs.map(f => hj(f, a / 5) / km) }], explain: 'Interface resistance reduces the benefit of conductive fillers, and for particles below the critical radius the conductivity falls with loading.' },
    ],
    equations: [
      eqn('governing', 'Kinetic theory of conduction', 'k = ⅓ C v Λ'),
      eqn('constitutive', 'Matthiessen\'s rule with boundary scattering', '1/Λ_eff = 1/Λ_bulk + 1/Λ_b,   Λ_b = d (1 + p)/(1 − p)'),
      eqn('constitutive', 'Fuchs–Sondheimer (film, in-plane)', 'k/k_bulk = 1 − [3(1 − p)/(2δ)] ∫₁^∞ (t⁻³ − t⁻⁵)(1 − e^(−δt))/(1 − p e^(−δt)) dt,   δ = d/Λ'),
      eqn('constitutive', 'Maxwell and Hasselman–Johnson effective conductivity', 'k_eff/k_m = [k_p(1+2α) + 2k_m + 2φ(k_p(1−α) − k_m)] / [k_p(1+2α) + 2k_m − φ(k_p(1−α) − k_m)],   α = R_K k_m / a   (α = 0: Maxwell)'),
      eqn('formula', 'Knudsen number', 'Kn = Λ / d'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', kHJ < km ? `The particles lower the conductivity (radius ${F(a * 1e9, 3)} nm is below the critical ${Number.isFinite(aCrit) ? F(aCrit * 1e9, 3) + ' nm' : 'radius'}). For a heat-transfer fluid or thermal interface material use larger or high-aspect-ratio fillers and improve interfacial bonding (functionalisation) to lower R_K.` : `The filler raises conductivity by ${F(100 * (kHJ / km - 1), 3)} % at ${F(100 * phi, 3)} vol%; interface resistance costs ${F(100 * (1 - (kHJ - km) / Math.max(kMax - km, 1e-300)), 3)} % of the ideal gain, so interface engineering is ${alpha > 0.3 ? 'worthwhile' : 'of minor benefit'}.`, 'Hasselman–Johnson against Maxwell conductivity.'),
      rec('medium', Kn > 0.1 ? `Use the size-reduced conductivity (${F(kWire, 3)} W/m/K for the wire, ${F(kFilm, 3)} in-plane for the film) in device thermal models; using the bulk value would underestimate temperature rises by a factor of up to ${F(kb / Math.min(kWire, kCross), 3)}.` : 'The structure is large compared with the mean free path: bulk conductivity can be used.', 'Knudsen number.'),
      rec('medium', 'For thermoelectric or thermal-barrier use, the same size effect is an advantage: nanostructure at the scale of the phonon mean free path while keeping electron paths (much shorter mean free path) unaffected.', 'Different mean free paths of heat and charge carriers.'),
      rec('low', `High filler loadings raise viscosity, pumping power and material cost; at ${F(100 * phi, 3)} vol% the gain is ${F(100 * (kHJ / km - 1), 3)} %. Check that the heat-transfer benefit outweighs the added pumping energy and the embodied impact of the nanoparticles.`, 'Modest effective-medium enhancement at dilute loading.', true),
    ],
    publish: { 'eff.k': kHJ, 'therm.k_wire': kWire, 'therm.k_film': kFilm, 'therm.Kn': Kn, 'therm.mfp': Lam },
    numerics: { method: 'Closed-form kinetic-theory and effective-medium relations; Simpson integration of the Fuchs–Sondheimer integral', dof: 1, steps: 200, iterations: ds.length, residual: Math.abs(fuchs(1000, p) - (1 - 3 * (1 - p) / 8000)), notes: 'Grey mean-free-path model.' },
  };
}

// ---------------------------------------------------------------- thermal metamaterials: laminate, cloak, concentrator, rotator
function runMeta(inp) {
  const kA = num(inp, 'kA'), kB = num(inp, 'kBl'), fA = num(inp, 'fA'), th = rad(num(inp, 'ang')), kbg = num(inp, 'kbg'), k1 = num(inp, 'kin'), k2 = num(inp, 'kout'), k0 = num(inp, 'kcore'), warnings = [], checks = [];
  const a = num(inp, 'ra'), b = Math.max(num(inp, 'rb'), a * 1.02), c = Math.max(num(inp, 'rc'), b * 1.02);
  if (b !== num(inp, 'rb') || c !== num(inp, 'rc')) warnings.push('Radii were re-ordered so that a < b < c (each layer at least 2 % thick).');
  const kpar = fA * kA + (1 - fA) * kB, kperp = 1 / (fA / kA + (1 - fA) / kB), cs = Math.cos(th), sn = Math.sin(th);
  const kxx = kpar * cs * cs + kperp * sn * sn, kyy = kpar * sn * sn + kperp * cs * cs, kxy = (kpar - kperp) * sn * cs, defl = deg(Math.atan2(kxy, kxx)), aniso = kpar / kperp, detK = kxx * kyy - kxy * kxy;
  // bilayer cloak: exact 2-D solution for a core + two shells in a uniform gradient (radii scaled by c)
  const solveShell = (kc, ki, ko, A_, B_) => { const M = [[A_, -A_, -1 / A_, 0, 0, 0], [kc, -ki, ki / (A_ * A_), 0, 0, 0], [0, B_, 1 / B_, -B_, -1 / B_, 0], [0, ki, -ki / (B_ * B_), -ko, ko / (B_ * B_), 0], [0, 0, 0, 1, 1, -1], [0, 0, 0, ko, -ko, kbg]]; return safe(() => N.solveLinear(M, [0, 0, 0, 0, 1, kbg]), null); };
  const A_ = a / c, B_ = b / c, s = solveShell(k0, k1, k2, A_, B_) || [1, 1, 0, 1, 0, 0], dist = Math.abs(s[5]), shield = Math.abs(s[0]);
  const bare = Math.abs((kbg - k0) / (kbg + k0)) * A_ * A_, k2req = kbg * (c * c + b * b) / (c * c - b * b), k2req3 = kbg * (2 * c ** 3 + b ** 3) / (2 * (c ** 3 - b ** 3)), cOpt = k2 > kbg ? b * Math.sqrt((k2 + kbg) / (k2 - kbg)) : null;
  const res = (() => { const M = [[A_, -A_, -1 / A_], [k0, -k1, k1 / (A_ * A_)]]; return Math.abs(M[0][0] * s[0] + M[0][1] * s[1] + M[0][2] * s[2]) + Math.abs(M[1][0] * s[0] + M[1][1] * s[1] + M[1][2] * s[2]) / Math.max(k0, k1); })();
  // concentrator / rotator from the laminate used as an anisotropic shell between a and c
  const l = Math.sqrt(kperp / kpar), conc = Math.pow(c / a, 1 - l), shieldRing = Math.pow(c / a, 1 - 1 / l), match = Math.sqrt(kpar * kperp) / kbg, rot = deg((kxy / kxx) * Math.log(c / a));
  // temperature field map (unit gradient along x)
  const nm = 61, half = 2 * c, data = new Float64Array(nm * nm);
  for (let j = 0; j < nm; j++) for (let i = 0; i < nm; i++) { const x = (-1 + 2 * i / (nm - 1)) * 2, y = (-1 + 2 * j / (nm - 1)) * 2, r2 = x * x + y * y, r = Math.sqrt(r2); let T; if (r < A_) T = s[0] * x; else if (r < B_) T = s[1] * x + s[2] * x / r2; else if (r < 1) T = s[3] * x + s[4] * x / r2; else T = x + s[5] * x / r2; data[j * nm + i] = T; }
  checks.push(chk('Conductivity tensor positive definite (second law)', kxx > 0 && kyy > 0 && detK > 0, `det = ${F(detK, 4)} (W/m/K)², eigenvalues ${F(kpar, 4)} and ${F(kperp, 4)} W/m/K`));
  checks.push(chk('Laminate conductivities within the Wiener bounds', kperp <= kpar * (1 + 1e-12) && kperp >= Math.min(kA, kB) * (1 - 1e-12) && kpar <= Math.max(kA, kB) * (1 + 1e-12), `series ${F(kperp, 4)} ≤ parallel ${F(kpar, 4)} W/m/K`));
  checks.push(chk('Interface conditions of the cloak solution satisfied', res < 1e-8, `residual ${F(res, 2)}`));
  checks.push(chk('Cloak leaves the external field undisturbed (dipole < 2 % of applied)', dist < 0.02, `relative dipole strength ${F(dist, 3)}; bare object ${F(bare, 3)}`));
  checks.push(chk('Outer-layer conductivity above the background (required for a bilayer cloak)', k2 > kbg, `k₂ = ${F(k2, 4)}, background ${F(kbg, 4)}, design value ${F(k2req, 4)} W/m/K`));
  if (dist >= 0.02) warnings.push(`The cloak is detuned: the outer layer should have ${F(k2req, 4)} W/m/K for these radii (or outer radius ${cOpt ? F(cOpt * 1e3, 4) + ' mm' : 'not achievable'} for this material). The remaining distortion is ${F(100 * dist, 3)} % of the applied gradient at the cloak surface.`);
  warnings.push('Steady-state, two-dimensional, perfect thermal contact between layers. In transients the cloak is imperfect unless ρc_p is also matched, and interface (contact) resistances shift the design value of the outer layer.');
  warnings.push('The laminate is treated as a homogeneous anisotropic medium, valid only when the layer thickness is much smaller than the shell thickness and the length scale of the temperature field.');
  const ths = N.linspace(0, 90, 91), cr = N.linspace(b * 1.02, Math.max(3 * b, c * 1.2), 80);
  return {
    summary: `The A/B laminate has conductivities of ${F(kpar, 4)} W/m/K along and ${F(kperp, 4)} W/m/K across the layers (anisotropy ${F(aniso, 3)}); tilted by ${F(deg(th), 3)}° it deflects the heat flux by ${F(defl, 3)}° from the applied gradient. A bilayer cloak with these radii needs an outer layer of ${F(k2req, 4)} W/m/K (2-D); with ${F(k2, 4)} W/m/K the external disturbance is ${F(100 * dist, 3)} % (bare object: ${F(100 * bare, 3)} %) and the gradient inside the core is ${F(100 * shield, 3)} % of the applied one. As a radial-spoke shell the laminate concentrates the gradient by a factor ${F(conc, 3)}; as a spiral it rotates the flux by ${F(rot, 3)}°.`,
    kpis: [
      kpi('kpar', 'In-plane (parallel) conductivity', kpar, 'W/m/K', 'Arithmetic (parallel) average along the layers: the upper Wiener bound.'),
      kpi('kperp', 'Cross-plane (series) conductivity', kperp, 'W/m/K', 'Harmonic (series) average across the layers: the lower Wiener bound.'),
      kpi('aniso', 'Anisotropy ratio', aniso, '', 'k_∥/k_⊥. Large ratios are what make flux guiding, concentration and rotation possible.'),
      kpi('kxy', 'Off-diagonal conductivity k_xy', kxy, 'W/m/K', `Tensor component in the rotated frame (k_xx = ${F(kxx, 4)}, k_yy = ${F(kyy, 4)} W/m/K). It can take either sign with the tilt direction; the tensor stays positive definite.`),
      kpi('defl', 'Heat-flux deflection angle', defl, '°', 'Angle between the heat flux and −∇T for a gradient along x, atan(k_xy/k_xx): the flux bends towards the layer direction.'),
      kpi('k2req', 'Required outer-layer conductivity (2-D cloak)', k2req, 'W/m/K', `κ_b(c² + b²)/(c² − b²) for a perfectly insulating inner layer; the 3-D (spherical) value is ${F(k2req3, 4)} W/m/K.`),
      kpi('dist', 'External field disturbance with the cloak', dist, '', 'Dipole coefficient of the exterior solution relative to the applied gradient at r = c. Zero is a perfect cloak.', dist < 0.02),
      kpi('bare', 'External disturbance of the bare object', bare, '', 'Same measure for the object alone in the background, |κ_b − κ₀|/(κ_b + κ₀)·(a/c)²: what the cloak has to cancel.'),
      kpi('shield', 'Core gradient / applied gradient', shield, '', 'Temperature gradient inside the cloaked object relative to the applied one. Near zero means the object is thermally protected (isothermal).', shield < 0.05),
      kpi('conc', 'Concentration factor (radial-spoke shell)', conc, '', `(c/a)^(1−l) with l = √(k_⊥/k_∥) = ${F(l, 3)}: gradient amplification in the core of an anisotropic shell with layers along the radius. Impedance match √(k_∥k_⊥)/κ_b = ${F(match, 3)} (one for no external disturbance).`),
      kpi('ring', 'Gradient reduction (concentric-ring shell)', shieldRing, '', '(c/a)^(1−1/l): the same laminate wound as concentric rings shields the core instead.'),
      kpi('rot', 'Flux rotation in a spiral shell', rot, '°', '(k_rθ/k_rr)·ln(c/a): rotation of the temperature-gradient direction inside the core when the laminate spirals at the stated angle to the radius.'),
    ],
    series: [
      { key: 'tens', title: 'Conductivity tensor vs layer angle', xlabel: 'layer angle (°)', ylabel: 'conductivity (W/m/K)', lines: [{ name: 'k_xx', x: ths, y: ths.map(t => kpar * Math.cos(rad(t)) ** 2 + kperp * Math.sin(rad(t)) ** 2) }, { name: 'k_yy', x: ths, y: ths.map(t => kpar * Math.sin(rad(t)) ** 2 + kperp * Math.cos(rad(t)) ** 2) }, { name: 'k_xy', x: ths, y: ths.map(t => (kpar - kperp) * Math.sin(rad(t)) * Math.cos(rad(t))) }], explain: 'Components of the effective tensor as the laminate is rotated. The off-diagonal term peaks at 45° and is responsible for flux deflection and rotation.' },
      { key: 'cloak', title: 'Cloak detuning vs outer radius', xlabel: 'outer radius c (mm)', ylabel: 'external disturbance', logy: true, lines: [{ name: 'relative dipole strength', x: cr.map(x => x * 1e3), y: cr.map(x => { const q = solveShell(k0, k1, k2, a / x, b / x); return Math.max(q ? Math.abs(q[5]) : 1, 1e-8); }) }, { name: '2 % criterion', x: [cr[0] * 1e3, cr[cr.length - 1] * 1e3], y: [0.02, 0.02] }], explain: 'External distortion for the chosen layer materials as the outer radius varies. The sharp minimum is the cloaking condition; its width shows the manufacturing tolerance.' },
    ],
    maps: [{ key: 'T', title: 'Temperature field around the cloak', nx: nm, ny: nm, data, x0: -half, x1: half, y0: -half, y1: half, xlabel: 'x (m)', ylabel: 'y (m)', unit: 'K per (K/m · c)', explain: 'Analytical 2-D solution for a unit gradient along x (temperature in units of gradient × outer radius). Straight, evenly spaced isotherms outside r = c indicate cloaking; a flat region in the centre indicates shielding.', diverging: true }],
    equations: [
      eqn('constitutive', 'Layered-medium effective conductivities', 'k_∥ = f k_A + (1 − f) k_B;   1/k_⊥ = f/k_A + (1 − f)/k_B'),
      eqn('constitutive', 'Rotated conductivity tensor', 'k_xx = k_∥ cos²θ + k_⊥ sin²θ;   k_yy = k_∥ sin²θ + k_⊥ cos²θ;   k_xy = (k_∥ − k_⊥) sin θ cos θ'),
      eqn('governing', 'Steady anisotropic conduction', '∇·(κ ∇T) = 0,   q = −κ ∇T'),
      eqn('formula', 'Bilayer cloak condition', 'κ₂ = κ_b (c² + b²)/(c² − b²)  (2-D);   κ₂ = κ_b (2c³ + b³)/[2(c³ − b³)]  (3-D);   κ₁ → 0'),
      eqn('formula', 'Anisotropic shell: concentration and rotation', 'T ∝ r^(±l) cos θ,  l = (k_θ/k_r)^½;   gain = (c/a)^(1−l);   rotation = (k_rθ/k_rr) ln(c/a)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', dist >= 0.02 ? `Retune the cloak: choose an outer-layer material near ${F(k2req, 4)} W/m/K, or keep the material and set the outer radius to ${cOpt ? F(cOpt * 1e3, 4) + ' mm' : 'a value that does not exist because the outer layer must be more conductive than the background'}.` : `The bilayer cloak is tuned (disturbance ${F(100 * dist, 3)} %). Hold the outer radius within the width of the minimum in the detuning curve and ensure good thermal contact (solder or thermal grease) at the outer interface.`, 'Exact solution of the layered cylinder.'),
      rec('medium', shield > 0.05 ? `The core still sees ${F(100 * shield, 3)} % of the applied gradient; use a lower-conductivity or thicker inner layer (aerogel, evacuated gap) to protect the object.` : `The core is shielded to ${F(100 * shield, 3)} % of the applied gradient, which protects temperature-sensitive components.`, 'Core field coefficient.'),
      rec('medium', aniso < 10 ? `Anisotropy is only ${F(aniso, 3)}; increase the conductivity contrast of the layers or move the fraction of A towards 0.5 for stronger concentration or rotation.` : `Anisotropy of ${F(aniso, 3)} gives a concentration factor of ${F(conc, 3)}; for no external disturbance choose the background so that √(k_∥k_⊥) = κ_b (currently ${F(match, 3)}).`, 'Laminate anisotropy and shell impedance.'),
      rec('low', 'Thermal metamaterials built from two common bulk materials (e.g. copper and a polymer) are passive, need no rare elements and can redirect waste heat to where it can be recovered or away from sensitive parts without active cooling power.', 'Passive nature of effective-medium thermal devices.', true),
    ],
    publish: { 'eff.k': kpar, 'therm.k_perp': kperp, 'therm.anisotropy': aniso, 'therm.cloak_disturbance': dist, 'therm.concentration': conc },
    numerics: { method: 'Closed-form laminate homogenisation; direct solution of the 6 × 6 interface system for the layered cylinder', dof: 6, steps: cr.length, iterations: cr.length + 1, residual: res, notes: 'Radii scaled by the outer radius to keep the system well conditioned.' },
  };
}

// ---------------------------------------------------------------- negative thermal expansion
function runNTE(inp) {
  const aM = num(inp, 'aM'), KM = num(inp, 'KM') * 1e9, GM = num(inp, 'GM') * 1e9, aN = num(inp, 'aN'), KN = num(inp, 'KN') * 1e9, GN = num(inp, 'GN') * 1e9, V = num(inp, 'Vf'), dT = num(inp, 'dTn'), aHi = num(inp, 'aHi'), aLo = Math.min(num(inp, 'aLo'), aHi), ph = rad(num(inp, 'phiL')), warnings = [], checks = [];
  const turner = v => (aM * KM * (1 - v) + aN * KN * v) / (KM * (1 - v) + KN * v), rom = v => aM * (1 - v) + aN * v;
  const hsA = v => (Math.abs(KN - KM) < 1e-3 * KM ? KM : KM + v / (1 / (KN - KM) + 3 * (1 - v) / (3 * KM + 4 * GM))), hsB = v => (Math.abs(KN - KM) < 1e-3 * KM ? KM : KN + (1 - v) / (1 / (KM - KN) + 3 * v / (3 * KN + 4 * GN)));
  const levin = (v, Ks) => (Math.abs(KN - KM) < 1e-3 * KM ? rom(v) : aM + (aN - aM) * (1 / Ks - 1 / KM) / (1 / KN - 1 / KM));
  const Ka = hsA(V), Kb = hsB(V), Klo = Math.min(Ka, Kb), Khi = Math.max(Ka, Kb), a1 = levin(V, Ka), a2 = levin(V, Kb), aT = turner(V), alo = Math.min(a1, a2), ahi = Math.max(a1, a2);
  const zero = f => { if (f(0) === 0) return 0; if (f(0) * f(0.999) > 0) return null; return safe(() => N.brent(f, 0, 0.999, 1e-12), null); };
  const V0T = zero(turner), V0a = zero(v => levin(v, hsA(v))), V0b = zero(v => levin(v, hsB(v))), V0R = zero(rom);
  const eps = (aN - aM) * dT, pMis = 12 * GM * KN * eps / (3 * KN + 4 * GM), strain = a1 * dT;
  const lat = p => (aLo - aHi * Math.sin(p) ** 2) / Math.cos(p) ** 2, aLat = lat(ph), ph0 = deg(Math.asin(Math.sqrt(N.clamp(aLo / aHi, 0, 1))));
  checks.push(chk('Hashin–Shtrikman bulk-modulus bounds ordered and within the phase values', Klo <= Khi * (1 + 1e-12) && Klo >= Math.min(KM, KN) * (1 - 1e-9) && Khi <= Math.max(KM, KN) * (1 + 1e-9), `${F(Klo / 1e9, 4)} ≤ K ≤ ${F(Khi / 1e9, 4)} GPa`));
  checks.push(chk('Composite CTE between the constituent values', aT >= Math.min(aM, aN) - 1e-15 && aT <= Math.max(aM, aN) + 1e-15, `Turner ${F(aT * 1e6, 4)} ppm/K; constituents ${F(aM * 1e6, 4)} and ${F(aN * 1e6, 4)} ppm/K`));
  checks.push(chk('Positive bulk moduli (mechanical stability; the sign of the CTE is not restricted by thermodynamics)', KM > 0 && KN > 0 && Klo > 0, `K_eff ≥ ${F(Klo / 1e9, 4)} GPa`));
  checks.push(chk('Thermal-mismatch pressure in the filler below 200 MPa', Math.abs(pMis) < 2e8, `${F(pMis / 1e6, 4)} MPa for ΔT = ${F(dT, 4)} K (dilute-sphere estimate)`));
  checks.push(chk('Lattice members remain a valid triangle (half angle < 85°)', deg(ph) < 85, `half apex angle ${F(deg(ph), 3)}°`));
  if (Math.abs(pMis) >= 2e8) warnings.push(`The estimated mismatch pressure of ${F(Math.abs(pMis) / 1e6, 3)} MPa is large. Many framework NTE compounds transform or amorphise under pressure (for ZrW₂O₈ a transition to a denser, positive-expansion phase is reported near 0.2 GPa), which destroys the NTE effect; matrix yielding and interface debonding are also likely.`);
  if (aN >= 0 && aM >= 0) warnings.push('Both constituents have non-negative expansion: the composite CTE cannot be negative, only reduced. The bi-material lattice result is independent of this.');
  warnings.push('Isotropic two-phase composite with perfect bonding and temperature-independent properties; Turner assumes equal strain in both phases, Levin\'s relation with Hashin–Shtrikman moduli gives rigorous bounds. Porosity and microcracking, common in NTE ceramics, lower the effective stiffness of the filler.');
  warnings.push('The lattice value is for an ideal pin-jointed bi-material triangle; bonded (rigid) joints reduce the effect and add thermal stresses.');
  const vs = N.linspace(0, 0.95, 96), ps = N.linspace(1, 80, 80), fmtV = v => (v === null ? 'none in 0–1' : v);
  return {
    summary: `At ${F(100 * V, 3)} vol% filler the effective linear CTE is ${F(aT * 1e6, 4)} ppm/K by Turner's model and between ${F(alo * 1e6, 4)} and ${F(ahi * 1e6, 4)} ppm/K by the Schapery (Levin/Hashin–Shtrikman) bounds (rule of mixtures: ${F(rom(V) * 1e6, 4)} ppm/K). Zero expansion needs at least ${V0T === null ? 'a filler fraction outside 0–1' : F(100 * V0T, 3) + ' vol%'} (Turner) and ${V0a === null ? 'is not reached below 100 vol%' : F(100 * V0a, 3) + ' vol%'} for a matrix-continuous microstructure. A ${F(dT, 3)} K excursion produces a mismatch pressure of about ${F(pMis / 1e6, 3)} MPa in the filler. The bi-material triangular lattice has a height-direction CTE of ${F(aLat * 1e6, 4)} ppm/K and zero expansion at a half angle of ${F(ph0, 3)}°.`,
    kpis: [
      kpi('aT', 'Effective CTE (Turner)', aT * 1e6, 'ppm/K', 'Bulk-modulus-weighted average Σα_iK_iV_i/ΣK_iV_i (uniform-strain assumption). It is the most optimistic estimate of the filler effect and can lie outside the rigorous bounds when the stiffness contrast is large. Negative values mean the composite contracts on heating.'),
      kpi('aU', 'Effective CTE, upper Schapery bound (Levin + Hashin–Shtrikman)', ahi * 1e6, 'ppm/K', 'Levin\'s exact relation between effective CTE and effective bulk modulus, evaluated at one Hashin–Shtrikman bound. For particles in a continuous softer matrix the measured value usually lies close to this bound.'),
      kpi('aL', 'Effective CTE, lower Schapery bound (Levin + Hashin–Shtrikman)', alo * 1e6, 'ppm/K', 'The other bound (stiffer phase continuous); the true CTE of a well-bonded isotropic composite lies between the two bounds.'),
      kpi('aR', 'Rule-of-mixtures CTE', rom(V) * 1e6, 'ppm/K', 'Volume average, valid only when both phases have equal stiffness; shown to illustrate the effect of stiffness mismatch.'),
      kpi('V0', 'Filler fraction for zero expansion (Turner)', V0T === null ? 'none in 0–1' : V0T, '', 'V = α_mK_m/(α_mK_m − α_fK_f): the smallest filler fraction that could give zero expansion (uniform-strain limit). It exists only when the constituent CTEs have opposite signs.', V0T !== null && V0T < 0.6 ? true : null),
      kpi('V0HS', 'Filler fraction for zero expansion (bounds)', V0a === null && V0b === null ? 'none in 0–1' : `${F(Math.min(fmtV(V0a) === 'none in 0–1' ? 1 : V0a, fmtV(V0b) === 'none in 0–1' ? 1 : V0b), 3)}–${F(Math.max(V0a === null ? 0 : V0a, V0b === null ? 0 : V0b), 3)}`, '', 'Range of zero-expansion fractions from the two Levin/Hashin–Shtrikman bounds; the experimental value should fall inside it.'),
      kpi('K', 'Effective bulk modulus (Hashin–Shtrikman range)', 0.5 * (Klo + Khi) / 1e9, 'GPa', `Mid-point of the bounds ${F(Klo / 1e9, 4)}–${F(Khi / 1e9, 4)} GPa at the stated filler fraction.`),
      kpi('p', 'Thermal-mismatch pressure in the filler', pMis / 1e6, 'MPa', '12G_mK_fε*/(3K_f + 4G_m) for an isolated sphere with misfit strain ε* = (α_f − α_m)ΔT. Positive is compression of the filler.', Math.abs(pMis) < 2e8),
      kpi('strain', 'Thermal strain over the excursion', strain * 1e6, 'µε', `α_eff·ΔT with the matrix-continuous Levin/Hashin–Shtrikman value (${F(a1 * 1e6, 4)} ppm/K, the value published as the effective CTE); compare with the dimensional tolerance of the application.`),
      kpi('aLat', 'Lattice CTE (height of bi-material triangle)', aLat * 1e6, 'ppm/K', '(α_legs − α_base sin²φ)/cos²φ: the expanding base widens the triangle and pulls the apex down, giving negative expansion from two positive-expansion materials.'),
      kpi('ph0', 'Lattice half angle for zero expansion', ph0, '°', 'sin²φ₀ = α_legs/α_base. Wider triangles give negative, narrower ones positive expansion.'),
    ],
    series: [
      { key: 'av', title: 'Effective CTE vs filler fraction', xlabel: 'filler volume fraction', ylabel: 'CTE (ppm/K)', lines: [{ name: 'Turner', x: vs, y: vs.map(v => turner(v) * 1e6) }, { name: 'Levin + HS (matrix-continuous)', x: vs, y: vs.map(v => levin(v, hsA(v)) * 1e6) }, { name: 'Levin + HS (filler-continuous)', x: vs, y: vs.map(v => levin(v, hsB(v)) * 1e6) }, { name: 'rule of mixtures', x: vs, y: vs.map(v => rom(v) * 1e6) }, { name: 'zero', x: [0, 0.95], y: [0, 0] }], explain: 'Models of the composite expansion. Crossing the zero line gives a zero-expansion composite; the spread between the curves shows the sensitivity to microstructure.' },
      { key: 'lat', title: 'Lattice CTE vs half apex angle', xlabel: 'half apex angle (°)', ylabel: 'CTE (ppm/K)', lines: [{ name: 'height direction', x: ps, y: ps.map(x => N.clamp(lat(rad(x)) * 1e6, -2000, 2000)) }, { name: 'zero', x: [1, 80], y: [0, 0] }], explain: 'Geometric tuning of the bi-material triangle (clipped at ±2000 ppm/K): the expansion passes through zero and becomes strongly negative for flat triangles, at the cost of stiffness.' },
    ],
    equations: [
      eqn('constitutive', 'Turner model', 'α_eff = Σ α_i K_i V_i / Σ K_i V_i'),
      eqn('constitutive', 'Levin relation (two-phase isotropic composite)', 'α_eff = α_m + (α_f − α_m)(1/K_eff − 1/K_m)/(1/K_f − 1/K_m)'),
      eqn('constitutive', 'Hashin–Shtrikman bulk-modulus bound', 'K = K₁ + V₂ / [1/(K₂ − K₁) + 3V₁/(3K₁ + 4G₁)]'),
      eqn('formula', 'Misfit pressure of a spherical inclusion', 'p = 12 G_m K_f ε* / (3K_f + 4G_m),   ε* = (α_f − α_m) ΔT'),
      eqn('formula', 'Bi-material triangle', 'H² = l² − b²/4   ⇒   α_H = (α_l − α_b sin²φ) / cos²φ,   sin φ = b/2l'),
      eqn('formula', 'Example NTE filler', chem.prettyFormula('ZrW2O8'), `zirconium tungstate, M = ${F(chem.molarMass('ZrW2O8'), 5)} g/mol; isotropic negative expansion over a wide temperature range`),
    ],
    checks, warnings,
    recommendations: [
      rec('high', V0T === null ? 'Zero expansion cannot be reached with these constituents; choose a filler with negative CTE or use the bi-material lattice route.' : (V0a === null || V0a > 0.6) ? `Zero expansion needs between ${F(100 * V0T, 3)} vol% (Turner, optimistic) and ${V0a === null ? 'more than 100' : F(100 * V0a, 3)} vol% (matrix-continuous bound) of filler, beyond practical packing for a particulate composite. Use a stiffer or more strongly NTE filler, a lower-expansion matrix, or a filler-continuous architecture (sintered or infiltrated NTE preform), which moves the composite towards the Turner limit.` : `Zero expansion is expected between ${F(100 * V0T, 3)} and ${F(100 * V0a, 3)} vol% filler depending on microstructure; make a small composition series across this range and measure by dilatometry, since porosity and particle contiguity shift the optimum.`, 'Turner and Schapery (Levin/Hashin–Shtrikman) zero crossings.'),
      rec('medium', Math.abs(pMis) >= 1e8 ? `Mismatch pressure of ${F(Math.abs(pMis) / 1e6, 3)} MPa over ${F(Math.abs(dT), 3)} K risks filler phase change, matrix yielding or debonding. Limit the processing temperature excursion, use a compliant interlayer, or reduce the CTE contrast.` : `Mismatch pressure of ${F(Math.abs(pMis) / 1e6, 3)} MPa is moderate; verify dimensional stability by thermal cycling, as stress relaxation in the matrix makes the CTE hysteretic.`, 'Misfit-inclusion stress estimate.'),
      rec('medium', aLat < 0 ? `The lattice gives ${F(aLat * 1e6, 3)} ppm/K from ordinary alloys; use pinned or flexure joints, and note that in-plane stiffness falls as the triangle is flattened.` : `Open the half angle beyond ${F(ph0, 3)}° to make the lattice expansion negative.`, 'Bi-material triangle kinematics.'),
      rec('low', 'Architected NTE lattices achieve tunable or zero expansion with common, recyclable alloys instead of tungsten- or zirconium-based compounds; where a particulate filler is used, zero expansion avoids active thermal control and its energy use in precision equipment.', 'Material choice for NTE function.', true),
    ],
    publish: { 'eff.cte': a1, 'eff.K': Ka, 'therm.cte_turner': aT, 'therm.cte_lattice': aLat, 'mech.mismatch_pressure': pMis },
    numerics: { method: 'Closed-form micromechanics (Turner, Levin, Hashin–Shtrikman, Eshelby misfit sphere); Brent root finding for zero-expansion fractions', dof: 1, steps: vs.length, iterations: 4, residual: 0, notes: 'Linear thermoelasticity.' },
  };
}

// ---------------------------------------------------------------- thermal rectification / NDTR
function runRect(inp) {
  const kA0 = num(inp, 'kA0'), nA = num(inp, 'nA'), LA = num(inp, 'LA'), kB0 = num(inp, 'kB0'), nB = num(inp, 'nB'), LB = num(inp, 'LB'), Gi0 = num(inp, 'Gi0'), d0 = num(inp, 'del0'), sA = num(inp, 'sA'), sB = num(inp, 'sB'), Tm = num(inp, 'Tm'), warnings = [], checks = [];
  const dTmax = Math.min(num(inp, 'dTmax'), 1.8 * Tm), dTop = Math.min(num(inp, 'dTop'), dTmax);
  if (dTmax < num(inp, 'dTmax')) warnings.push(`The sweep range was limited to ±${F(dTmax, 4)} K so that both terminals stay above 10 % of the mean temperature.`);
  // Kirchhoff transform of k = k0 (T/Tm)^n and its inverse
  const phi = (T, k0, n) => (Math.abs(n + 1) < 1e-9 ? k0 * Tm * Math.log(T / Tm) : k0 * Tm / (n + 1) * Math.pow(T / Tm, n + 1));
  const inv = (P, k0, n) => { if (Math.abs(n + 1) < 1e-9) return Tm * Math.exp(N.clamp(P / (k0 * Tm), -50, 50)); const x = P * (n + 1) / (k0 * Tm); return x <= 0 ? (n + 1 > 0 ? 1e-3 : 1e9) : N.clamp(Tm * Math.pow(x, 1 / (n + 1)), 1e-3, 1e9); };
  const Gi = (T1, T2) => { const del = d0 + sA * (T1 - Tm) / Tm - sB * (T2 - Tm) / Tm; return Gi0 * Math.exp(-Math.min(del * del / 4, 700)); };
  let multi = false, evals = 0;
  const solve = dT => {
    if (Math.abs(dT) < 1e-12) return { J: 0, T1: Tm, T2: Tm, G: Gi(Tm, Tm) };
    const s = dT > 0 ? 1 : -1, TA = Tm + dT / 2, TB = Tm - dT / 2, t1 = j => inv(phi(TA, kA0, nA) - s * j * LA, kA0, nA), t2 = j => inv(phi(TB, kB0, nB) + s * j * LB, kB0, nB), gap = j => s * (t1(j) - t2(j));
    let jHi = Math.abs(dT) * Math.min(kA0 / LA, kB0 / LB, Gi0) + 1e-30; for (let i = 0; i < 80 && gap(jHi) > 0; i++) jHi *= 2;
    const jMax = safe(() => N.brent(gap, 0, jHi, 1e-14 * jHi), jHi), f = j => { evals++; const a = t1(j), b = t2(j); return s * Gi(a, b) * (a - b) - j; };
    const m = 60; let prev = f(0), lo = 0, root = null, changes = 0;
    for (let i = 1; i <= m; i++) { const j = jMax * i / m, cur = f(j); if (prev > 0 && cur <= 0) { changes++; if (root === null) root = safe(() => N.brent(f, lo, j, 1e-13 * jMax), 0.5 * (lo + j)); } else if (prev <= 0 && cur > 0) changes++; prev = cur; lo = j; }
    if (changes > 1) multi = true; if (root === null) root = jMax;
    return { J: s * root, T1: t1(root), T2: t2(root), G: Gi(t1(root), t2(root)) };
  };
  const np = 161, dTs = N.linspace(-dTmax, dTmax, np), sols = dTs.map(solve), Js = sols.map(q => q.J), dJ = Js.map((_, i) => (i === 0 ? (Js[1] - Js[0]) / (dTs[1] - dTs[0]) : i === np - 1 ? (Js[i] - Js[i - 1]) / (dTs[i] - dTs[i - 1]) : (Js[i + 1] - Js[i - 1]) / (dTs[i + 1] - dTs[i - 1])));
  const mid = (np - 1) / 2, region = (from, to, step) => { let on = null, off = null; for (let i = from; step > 0 ? i <= to : i >= to; i += step) { if (dJ[i] < 0) { if (on === null) on = dTs[i]; off = dTs[i]; } else if (on !== null) break; } return on === null ? null : [on, off]; };
  const fwd = region(mid + 1, np - 1, 1), rev = region(mid - 1, 0, -1), jf = solve(dTop), jr = solve(-dTop), Jf = Math.abs(jf.J), Jr = Math.abs(jr.J), Rr = Math.max(Jf, Jr) / Math.max(Math.min(Jf, Jr), 1e-300), easy = Jf >= Jr ? 'forward (A hot)' : 'reverse (B hot)';
  const iMin = N.argmin(dJ), gMin = dJ[iMin], G0 = dJ[mid], second = Js.every((J, i) => J * dTs[i] >= -1e-12 * Math.abs(J * dTs[i]) - 1e-30);
  let pvr = 1; { const side = fwd ? [mid, np - 1, 1] : rev ? [mid, 0, -1] : null; if (side) { let pk = 0, vl = Infinity, seen = false; for (let i = side[0]; side[2] > 0 ? i <= side[1] : i >= side[1]; i += side[2]) { const v = Math.abs(Js[i]); if (!seen) { if (v >= pk) pk = v; else { seen = true; vl = v; } } else vl = Math.min(vl, v); } if (seen && vl > 0) pvr = pk / vl; } }
  const fracI = dT => { const q = solve(dT); return Math.abs(q.T1 - q.T2) / Math.abs(dT); }, fi = fracI(dTop), ndtr = fwd || rev, GA = kA0 / LA, GB = kB0 / LB, Gser = 1 / (1 / GA + 1 / GB + 1 / Gi(Tm, Tm));
  checks.push(chk('Second law: heat flows from hot to cold at every bias (J·ΔT ≥ 0)', second, 'entropy production J·ΔT/(T_h T_c) is non-negative over the whole sweep'));
  checks.push(chk('Zero-bias conductance equals the series combination of segments and interface', Math.abs(G0 - Gser) / Gser < 0.05, `numerical ${F(G0, 4)} vs series ${F(Gser, 4)} W/m²/K`));
  checks.push(chk('Single-valued solution (no thermal bistability in the sweep)', !multi, multi ? 'more than one steady state exists at some bias: the smallest-current branch is shown' : 'one steady state at every bias'));
  checks.push(chk('Interface-dominated junction (interface carries most of the temperature drop)', fi > 0.5, `${F(100 * fi, 3)} % of ΔT falls across the interface at the operating bias`));
  if (multi) warnings.push('The nonlinear interface gives more than one steady solution at some biases (thermal bistability, a known companion of negative differential thermal resistance). The branch with the smallest heat current is reported; a real device would show hysteresis.');
  if (!ndtr) warnings.push('No negative differential thermal resistance occurs in the swept range with these parameters; it requires an interface-dominated junction and a spectral mismatch that grows quickly with bias.');
  warnings.push('Reduced phenomenological model: segment conductivities follow power laws in temperature and the interface conductance is a Gaussian function of a temperature-dependent spectral mismatch. It reproduces the mechanism (bias-dependent overlap of vibrational spectra) but the parameters are not material constants; quantitative design needs lattice-dynamics or molecular-dynamics input.');
  return {
    summary: `At a bias of ${F(dTop, 4)} K about a mean of ${F(Tm, 4)} K the heat current is ${F(jf.J, 3)} W/m² forward and ${F(jr.J, 3)} W/m² in reverse: a rectification ratio of ${F(Rr, 3)} with the easy direction ${easy}. ${ndtr ? `Negative differential thermal resistance (dJ/dΔT < 0) occurs ${fwd ? `for forward bias between ${F(fwd[0], 3)} and ${F(fwd[1], 3)} K` : ''}${fwd && rev ? ' and ' : ''}${rev ? `for reverse bias between ${F(Math.abs(rev[0]), 3)} and ${F(Math.abs(rev[1]), 3)} K` : ''}, with a peak-to-valley current ratio of ${F(pvr, 3)}.` : 'No negative differential thermal resistance occurs in the swept range.'} This is a reduced phenomenological model.`,
    kpis: [
      kpi('Jf', 'Forward heat current', jf.J, 'W/m²', 'Heat flux with the A terminal hot, from solving flux continuity through segment A, the interface and segment B.'),
      kpi('Jr', 'Reverse heat current', jr.J, 'W/m²', 'Heat flux with the B terminal hot at the same magnitude of bias (negative sign: flow from B to A).'),
      kpi('R', 'Rectification ratio', Rr, '', `Larger over smaller of the two current magnitudes; one means no rectification. Easy direction: ${easy}.`, Rr > 2 ? true : null),
      kpi('on', 'NDTR onset bias', fwd ? fwd[0] : rev ? rev[0] : 'none in range', 'K', 'Smallest bias (signed: positive = A hot) at which the differential thermal conductance dJ/dΔT turns negative, i.e. more temperature difference gives less heat flow.'),
      kpi('off', 'NDTR end bias', fwd ? fwd[1] : rev ? rev[1] : 'none in range', 'K', 'Bias at which the differential conductance becomes positive again, or the end of the sweep if it stays negative.'),
      kpi('gmin', 'Minimum differential thermal conductance', gMin, 'W/m²/K', `Most negative slope dJ/dΔT in the sweep (at ΔT = ${F(dTs[iMin], 3)} K). A negative value is the thermal analogue of the tunnel-diode region and the basis of thermal transistors and logic.`),
      kpi('pvr', 'Peak-to-valley current ratio', pvr, '', 'Heat current at the local maximum divided by the minimum beyond it; the usable swing of a thermal switching element.'),
      kpi('G0', 'Zero-bias thermal conductance', G0, 'W/m²/K', `Linear-response conductance; segments contribute ${F(GA, 3)} and ${F(GB, 3)} W/m²/K and the interface ${F(Gi(Tm, Tm), 3)} W/m²/K in series.`),
      kpi('fi', 'Interface share of the temperature drop', fi, '', 'Fraction of the applied bias that falls across the interface at the operating point; interface-controlled behaviour needs a large share.'),
      kpi('Ti', 'Interface temperatures at forward bias', `${F(jf.T1, 4)} / ${F(jf.T2, 4)} K`, '', 'Temperatures on the A and B sides of the interface, which set the spectral mismatch.'),
    ],
    series: [
      { key: 'J', title: 'Heat current vs temperature bias', xlabel: 'ΔT = T_A − T_B (K)', ylabel: 'heat current (W/m²)', lines: [{ name: 'J(ΔT)', x: dTs, y: Js }, { name: 'linear response', x: [-dTmax, dTmax], y: [-G0 * dTmax, G0 * dTmax] }], explain: 'Thermal current–bias characteristic. Asymmetry between positive and negative bias is rectification; a region where the curve falls with increasing bias is negative differential thermal resistance.' },
      { key: 'dJ', title: 'Differential thermal conductance', xlabel: 'ΔT (K)', ylabel: 'dJ/dΔT (W/m²/K)', lines: [{ name: 'dJ/dΔT', x: dTs, y: dJ }, { name: 'zero', x: [-dTmax, dTmax], y: [0, 0] }], explain: 'Slope of the characteristic. Values below zero mark the NDTR window.' },
      { key: 'Gi', title: 'Interface conductance and temperature drop', xlabel: 'ΔT (K)', ylabel: 'G_i / G_i0  and interface share of ΔT', lines: [{ name: 'interface conductance / maximum', x: dTs, y: sols.map(q => q.G / Gi0) }, { name: 'interface share of |ΔT|', x: dTs, y: sols.map((q, i) => (Math.abs(dTs[i]) > 1e-9 ? Math.abs(q.T1 - q.T2) / Math.abs(dTs[i]) : fracI(1e-3))) }], explain: 'Bias moves the two spectra into or out of overlap: conductance peaks where they match and collapses where they separate, which produces both rectification and NDTR.' },
    ],
    equations: [
      eqn('governing', 'Flux continuity through the junction', 'J = (1/L_A) ∫_{T₁}^{T_A} k_A dT = G_i(T₁, T₂)(T₁ − T₂) = (1/L_B) ∫_{T_B}^{T₂} k_B dT'),
      eqn('constitutive', 'Segment conductivities', 'k_A = k_A0 (T/T_m)^(n_A),   k_B = k_B0 (T/T_m)^(n_B)'),
      eqn('constitutive', 'Interface conductance from spectral overlap (phenomenological)', 'G_i = G_i0 exp(−δ²/4),   δ = δ₀ + s_A (T₁ − T_m)/T_m − s_B (T₂ − T_m)/T_m'),
      eqn('formula', 'Rectification ratio and NDTR criterion', 'R = max(|J₊|, |J₋|) / min(|J₊|, |J₋|);   NDTR:  dJ/dΔT < 0'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', ndtr ? `Operate inside the NDTR window (${fwd ? `${F(fwd[0], 3)}–${F(fwd[1], 3)} K forward` : `${F(Math.abs(rev[0]), 3)}–${F(Math.abs(rev[1]), 3)} K reverse`}) for switching or amplification, and outside it for a stable diode; the peak-to-valley ratio of ${F(pvr, 3)} sets the achievable gain.` : `To obtain NDTR, weaken the interface relative to the segments (G_i0 ≪ k/L, currently ${F(Gi0 / Math.min(GA, GB), 3)} of the weaker segment) and increase the band-shift coefficients so that the mismatch grows by more than one bandwidth within the bias range.`, 'Differential conductance over the sweep.'),
      rec('medium', Rr < 2 ? `Rectification is weak (${F(Rr, 3)}). Increase the zero-bias mismatch magnitude or use segments with opposite temperature exponents and comparable resistance to the interface.` : `Rectification of ${F(Rr, 3)} at ${F(dTop, 3)} K is useful for passive thermal management; verify experimentally with reversed-bias measurements on the same sample, correcting for contact resistances and radiation losses.`, 'Forward and reverse currents at the operating bias.'),
      rec('medium', multi ? 'Bistability is predicted: sweep the bias in both directions when measuring and expect hysteresis; use it deliberately for thermal memory.' : 'The characteristic is single-valued in this range; no hysteresis is expected.', 'Number of steady solutions found.'),
      rec('low', 'Passive thermal diodes need no electrical power: they can retain heat in a reactor or building element during cold periods and shed it during hot ones, reducing active heating and cooling energy.', 'Direction-dependent heat flow without external power.', true),
    ],
    publish: { 'therm.rectification': Rr, 'therm.ndtr_onset_K': fwd ? fwd[0] : rev ? rev[0] : 0, 'therm.conductance': G0, 'therm.heat_flux': jf.J },
    numerics: { method: 'Kirchhoff transform of the segment conductivities (closed form) and Brent root finding on the interface flux balance at each bias; central differences for dJ/dΔT', dof: 2, steps: np, iterations: evals, residual: 0, notes: 'Smallest-current root selected when several exist.' },
  };
}

const RUN = { slab: runSlab, lumped: runLumped, furnace: runFurnace, nano: runNano, meta: runMeta, nte: runNTE, rect: runRect };

export default {
  id: 10,
  tagline: 'Transient Joule and reaction heating, lumped and wall models with ramp–soak schedules, radiative furnace balance, nanoscale conduction, and thermal metamaterials including NTE composites and thermal rectifiers.',
  inputs,
  presets: [
    { name: 'ZnO hydrothermal reference slab: 1 m × 1 m × 0.01 m, 1000 V/m', note: 'Reference case: k = 0.6 W/m/K, ρ = 1000 kg/m³, c_p = 4186 J/kg/K, T₀ = 298 K, h = 50 W/m²/K, ambient 298 K, electric field 1000 V/m, reaction activation energy 65 kJ/mol.', values: { mode: 'slab', dim: '1d', Lz: 0.01, Lxw: 1, Lyd: 1, k: 0.6, rho: 1000, cp: 4186, T0: 298, h: 50, Tamb: 298, Ef: 1000, sig: 0.25, alS: 0.02, qR: 100, EaR: 65, tEnd: 3600, Tset: 320 } },
    { name: 'Reference slab in 2-D with half electrode coverage', note: 'Same slab with the field applied to the central half of the width: lateral gradients and edge losses.', values: { mode: 'slab', dim: '2d', hfrac: 0.5, Ef: 1000, sig: 0.25, tEnd: 3600 } },
    { name: 'Laboratory autoclave: ramp to 90 °C, 3 h soak', note: 'Lumped energy and time budget for a 1 L hydrothermal batch.', values: { mode: 'lumped', Vl: 1e-3, Al: 0.06, Ph: 300, Tsoak: 363, tSoak: 10800, Tcool: 313, h: 10, k: 0.6, rho: 1000, cp: 4186, T0: 298, Tamb: 298 } },
    { name: 'Cold-wall CVD / RTP wafer at 1000 °C heater', note: 'Radiation-dominated heat balance of a silicon wafer.', values: { mode: 'furnace', Th: 1273, epsW: 0.65, epsH: 0.85, Twall: 300, hg: 20, Tg: 600 } },
    { name: 'Silicon nanowire, 50 nm, and ZnO nanofluid', note: 'Classical size effect and Kapitza-limited nanocomposite conductivity.', values: { mode: 'nano', kb: 148, Cv: 1.66e6, vs: 6400, dch: 5e-8, spec: 0, RK: 1e-8, kpar: 50, kmat: 0.6, phiP: 0.05, aP: 2e-8 } },
    { name: 'Bilayer thermal cloak (polystyrene / alloy shell)', note: 'Thermal metamaterial: insulating inner layer and conducting outer layer tuned to the background.', values: { mode: 'meta', kbg: 2.3, kin: 0.03, kout: 9.8, kcore: 200, ra: 0.006, rb: 0.0095, rc: 0.012, kA: 400, kBl: 0.15, fA: 0.5, ang: 30 } },
    { name: 'Low-expansion epoxy / ZrW₂O₈ composite and NTE lattice', note: 'Negative-thermal-expansion filler offsetting the expansion of a polymer matrix; also evaluates the bi-material NTE lattice.', values: { mode: 'nte', aM: 60e-6, KM: 5, GM: 1.5, aN: -9e-6, KN: 72, GN: 34, Vf: 0.5, dTn: 50 } },
    { name: 'Low-expansion Cu / ZrW₂O₈ heat-sink composite', note: 'Metal-matrix composite for CTE matching to semiconductors.', values: { mode: 'nte', aM: 17e-6, KM: 140, GM: 48, aN: -9e-6, KN: 72, GN: 34, Vf: 0.5, dTn: 50 } },
    { name: 'Thermal diode with negative differential thermal resistance', note: 'Two-segment junction with a weak, spectrally mismatched interface.', values: { mode: 'rect', kA0: 1, nA: 1.5, LA: 1e-6, kB0: 1, nB: -1.5, LB: 1e-6, Gi0: 2e5, del0: 2, sA: 6, sB: 6, Tm: 300, dTmax: 200, dTop: 100 } },
  ],
  run(inp) { return RUN[sel(inp, 'mode')](inp); },
};
