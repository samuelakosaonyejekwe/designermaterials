// Suite 9 — Fluid dynamics and porous-media transport.
// Packed-bed (Darcy–Forchheimer/Ergun) hydraulics, 2-D heterogeneous Darcy flow with reactive
// advection–diffusion, channel and coil microreactor flow with residence-time distribution,
// stirred-vessel correlations and particle–fluid interaction.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, PI = Math.PI, kB = N.K_B, G = N.G0;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
// closed-vessel axial-dispersion conversion for a first-order reaction (Danckwerts)
const danckwerts = (Pe, Da) => { if (Da <= 0) return 0; if (Pe > 1e6) return 1 - Math.exp(-Da); if (Pe < 1e-6) return Da / (1 + Da); const a = Math.sqrt(1 + 4 * Da / Pe), ex = Math.exp(Math.max(Pe * (1 - a) / 2, -700)), den = (1 + a) ** 2 - (1 - a) ** 2 * Math.exp(Math.max(-a * Pe, -700)); return N.clamp(1 - 4 * a * ex / den, 0, 1); };
const ZNO_EQ = () => { const ie = safe(() => chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'NaNO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }]), null); return ie ? [eqn('chemical', 'Example precursor-consuming reaction (ZnO deposition)', ie.molecular, 'represented in the transport model by a first-order sink −k c'), eqn('ionic', 'Net ionic equation', ie.netIonic)] : []; };

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'bed', group: 'Model', options: [
    { value: 'bed', label: 'Packed bed / porous reactor (Ergun, Kozeny–Carman)' }, { value: 'darcy2d', label: '2-D heterogeneous Darcy flow with reactive transport' },
    { value: 'channel', label: 'Channel, pipe and coil microreactor flow' }, { value: 'stirred', label: 'Stirred vessel' }, { value: 'particle', label: 'Particle–fluid interaction' }], help: 'Selects the flow model. Each analysis publishes hydraulic results (flow, permeability, pressure drop, conversion) for the other suites.' },
  { key: 'rho', label: 'Fluid density', unit: 'kg/m³', value: 997, min: 0.05, max: 3000, log: true, group: 'Fluid', help: 'Density of the flowing fluid (997 kg/m³ for water at 25 °C, about 1.2 for air).' },
  { key: 'mu', label: 'Fluid viscosity', unit: 'Pa·s', value: 8.9e-4, min: 1e-6, max: 100, log: true, group: 'Fluid', help: 'Dynamic viscosity (8.9×10⁻⁴ Pa·s for water at 25 °C, 1.85×10⁻⁵ for air).' },
  { key: 'Dmol', label: 'Molecular diffusivity of the species', unit: 'm²/s', value: 7.5e-10, min: 1e-12, max: 1e-3, log: true, group: 'Transport', bus: 'prop.D_diff', showIf: { mode: ['bed', 'darcy2d', 'channel'] }, help: 'Diffusion coefficient of the transported precursor in the fluid (7.5×10⁻¹⁰ m²/s for Zn²⁺ in water).' },
  { key: 'k1', label: 'First-order reaction rate constant', unit: '1/s', value: 1e-3, min: 0, max: 1e3, group: 'Transport', bus: 'kin.k', showIf: { mode: ['bed', 'darcy2d', 'channel'] }, help: 'Rate constant of precursor consumption (per unit fluid volume). Zero gives a passive tracer.' },
  { key: 'Q', label: 'Volumetric flow rate', unit: 'm³/s', value: 1e-7, min: 1e-12, max: 10, log: true, group: 'Operation', bus: 'proc.flow', showIf: { mode: ['bed', 'channel'] }, help: 'Throughput (10⁻⁷ m³/s = 0.1 mL/s = 6 mL/min).' },
  // packed bed
  { key: 'por', label: 'Porosity', unit: '', value: 0.25, min: 0.05, max: 0.95, group: 'Porous medium', bus: 'geom.porosity', showIf: { mode: ['bed', 'darcy2d'] }, help: 'Void fraction of the bed. Permeability depends on it very steeply (ε³/(1−ε)²).' },
  { key: 'dg', label: 'Grain (particle) diameter', unit: 'm', value: 1e-4, min: 1e-8, max: 0.05, log: true, group: 'Porous medium', showIf: { mode: ['bed'] }, help: 'Sauter-mean diameter of the packing grains (10⁻⁴ m = 100 µm).' },
  { key: 'Lbed', label: 'Bed length', unit: 'm', value: 0.1, min: 1e-4, max: 20, log: true, group: 'Porous medium', showIf: { mode: ['bed'] }, help: 'Length of the porous bed in the flow direction.' },
  { key: 'Abed', label: 'Bed cross-sectional area', unit: 'm²', value: 1e-3, min: 1e-8, max: 100, log: true, group: 'Porous medium', showIf: { mode: ['bed'] }, help: 'Area normal to the flow; the superficial velocity is Q/A.' },
  // 2-D Darcy
  { key: 'Lx', label: 'Domain length (flow direction)', unit: 'm', value: 1, min: 1e-4, max: 100, log: true, group: 'Domain', showIf: { mode: ['darcy2d'] }, help: 'Distance between the inlet (left) and outlet (right) faces.' },
  { key: 'Ly', label: 'Domain width', unit: 'm', value: 1, min: 1e-4, max: 100, log: true, group: 'Domain', showIf: { mode: ['darcy2d'] }, help: 'Extent normal to the flow; the top and bottom sides are no-flow boundaries.' },
  { key: 'thick', label: 'Slab thickness', unit: 'm', value: 0.01, min: 1e-5, max: 10, log: true, group: 'Domain', showIf: { mode: ['darcy2d'] }, help: 'Out-of-plane thickness used to convert the 2-D solution into a flow rate.' },
  { key: 'Kg', label: 'Geometric-mean permeability', unit: 'm²', value: 1e-12, min: 1e-18, max: 1e-6, log: true, group: 'Porous medium', bus: 'eff.permeability', showIf: { mode: ['darcy2d'] }, help: 'Median of the log-normal permeability field (10⁻¹² m² ≈ 1 darcy).' },
  { key: 'sigK', label: 'Standard deviation of ln K', unit: '', value: 1, min: 0, max: 3, group: 'Porous medium', showIf: { mode: ['darcy2d'] }, help: 'Heterogeneity of the permeability field. Zero is homogeneous; above about 2 flow concentrates in a few channels.' },
  { key: 'corr', label: 'Correlation length (grid cells)', unit: 'cells', value: 3, min: 0, max: 10, step: 1, group: 'Porous medium', showIf: { mode: ['darcy2d'] }, help: 'Number of smoothing passes applied to the random field: larger values give larger patches of high and low permeability. The filter is scaled with the cell size so that the correlation length is the same physical distance along x and y, about √(2n/3) cells of the finer spacing for n passes.' },
  { key: 'dP', label: 'Applied pressure difference', unit: 'Pa', value: 1e5, min: 1, max: 1e8, log: true, group: 'Operation', showIf: { mode: ['darcy2d'] }, help: 'Inlet pressure minus outlet pressure.' },
  { key: 'pout', label: 'Outlet pressure', unit: 'Pa', value: 101325, min: 0, max: 1e8, group: 'Operation', bus: 'proc.P', showIf: { mode: ['darcy2d'] }, help: 'Absolute pressure at the outlet face; it only shifts the pressure map.' },
  { key: 'cin', label: 'Inlet precursor concentration', unit: 'mol/m³', value: 10, min: 1e-6, max: 1e5, log: true, group: 'Transport', bus: 'proc.c_precursor', showIf: { mode: ['darcy2d'] }, help: 'Concentration entering at the inlet face (10 mol/m³ = 0.01 mol/L).' },
  { key: 'grid', label: 'Grid resolution', type: 'select', value: 'medium', group: 'Numerics', showIf: { mode: ['darcy2d'] }, options: [{ value: 'coarse', label: '32 × 24 cells' }, { value: 'medium', label: '48 × 32 cells' }, { value: 'fine', label: '72 × 48 cells' }], help: 'Finite-volume grid. Finer grids resolve channelling better at higher cost; compare two grids to judge discretisation error.' },
  { key: 'nens', label: 'Realisations in the permeability ensemble', unit: '', value: 8, min: 1, max: 16, step: 1, sens: false, group: 'Numerics', showIf: { mode: ['darcy2d'] }, help: 'Number of independent seeded fields solved for the ensemble-mean effective permeability and its standard error. The maps and the transport solution always show the first realisation.' },
  { key: 'seed', label: 'Random seed of the permeability field', unit: '', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Numerics', showIf: { mode: ['darcy2d'] }, help: 'Seed of the pseudo-random generator; the same seed always gives the same field.' },
  // channel
  { key: 'geo', label: 'Channel geometry', type: 'select', value: 'pipe', group: 'Channel', showIf: { mode: ['channel'] }, options: [{ value: 'pipe', label: 'Straight tube' }, { value: 'coil', label: 'Helical coil' }], help: 'A coil adds secondary (Dean) flow that raises pressure drop and delays transition to turbulence.' },
  { key: 'dh', label: 'Tube inner diameter', unit: 'm', value: 1e-3, min: 1e-5, max: 2, log: true, group: 'Channel', showIf: { mode: ['channel'] }, help: 'Inner (hydraulic) diameter of the tube or channel.' },
  { key: 'Lch', label: 'Tube length', unit: 'm', value: 5, min: 1e-3, max: 1e4, log: true, group: 'Channel', showIf: { mode: ['channel'] }, help: 'Developed length of the tube; with the flow rate it gives the mean residence time.' },
  { key: 'rough', label: 'Wall roughness', unit: 'm', value: 1.5e-6, min: 0, max: 1e-2, group: 'Channel', showIf: { mode: ['channel'] }, help: 'Absolute roughness (about 1.5 µm for drawn tubing, 45 µm for commercial steel). Only matters in turbulent flow.' },
  { key: 'Dc', label: 'Coil diameter', unit: 'm', value: 0.05, min: 1e-3, max: 10, log: true, group: 'Channel', showIf: { mode: ['channel'] }, help: 'Diameter of the helix (centre line) for coiled tubes.' },
  // stirred vessel
  { key: 'imp', label: 'Impeller type', type: 'select', value: 'pbt', group: 'Vessel', showIf: { mode: ['stirred'] }, options: [{ value: 'rushton', label: 'Rushton turbine (radial)' }, { value: 'pbt', label: 'Pitched-blade turbine, 45° (mixed flow)' }, { value: 'hydrofoil', label: 'Hydrofoil (axial)' }], help: 'Sets the turbulent power number, flow number and Metzner–Otto constant for a baffled tank.' },
  { key: 'Dimp', label: 'Impeller diameter', unit: 'm', value: 0.05, min: 0.005, max: 3, log: true, group: 'Vessel', showIf: { mode: ['stirred'] }, help: 'Usually one third to one half of the tank diameter.' },
  { key: 'Tt', label: 'Tank diameter', unit: 'm', value: 0.15, min: 0.01, max: 10, log: true, group: 'Vessel', showIf: { mode: ['stirred'] }, help: 'Inner diameter of the vessel; the liquid height is taken equal to it.' },
  { key: 'Nrpm', label: 'Stirrer speed', unit: 'rpm', value: 800, min: 1, max: 5000, log: true, group: 'Vessel', showIf: { mode: ['stirred'] }, help: 'Rotational speed of the impeller. Power scales with its cube in turbulent flow.' },
  { key: 'dpS', label: 'Particle / agglomerate diameter', unit: 'm', value: 2e-5, min: 1e-9, max: 1e-2, log: true, group: 'Particles', bus: 'geom.diameter', showIf: { mode: ['stirred', 'particle'] }, help: 'Size of the suspended particles or agglomerates (2×10⁻⁵ m = 20 µm). Nanoparticles are usually suspended as micron-size agglomerates.' },
  { key: 'rhoS', label: 'Particle density', unit: 'kg/m³', value: 5606, min: 100, max: 25000, group: 'Particles', bus: 'prop.density', showIf: { mode: ['stirred', 'particle'] }, help: 'Density of the particle (use the effective density for porous agglomerates).' },
  { key: 'Xs', label: 'Solids loading', unit: 'mass %', value: 1, min: 0.001, max: 50, log: true, group: 'Particles', showIf: { mode: ['stirred'] }, help: 'Mass of solids per 100 mass units of liquid, as used in the Zwietering correlation.' },
  { key: 'Szw', label: 'Zwietering geometry constant S', unit: '', value: 6, min: 1, max: 20, group: 'Particles', sens: false, showIf: { mode: ['stirred'] }, help: 'Dimensionless constant of the just-suspended-speed correlation. Typical for D/T = 1/3 and clearance T/4: about 8 for a Rushton turbine, 5–6 for a down-pumping pitched-blade turbine, 6–7 for hydrofoils.' },
  // particle–fluid
  { key: 'Tp', label: 'Fluid temperature', unit: '°C', value: 25, min: -50, max: 1500, group: 'Fluid', bus: 'proc.T_C', showIf: { mode: ['particle'] }, help: 'Sets the thermal energy for Brownian motion.' },
  { key: 'phiS', label: 'Particle volume fraction', unit: '', value: 0.05, min: 0, max: 0.55, group: 'Particles', showIf: { mode: ['particle'] }, help: 'Solids volume fraction for hindered settling (Richardson–Zaki).' },
  { key: 'Uc', label: 'Characteristic flow velocity', unit: 'm/s', value: 0.1, min: 1e-6, max: 300, log: true, group: 'Operation', showIf: { mode: ['particle'] }, help: 'Velocity of the carrying flow used in the Stokes number.' },
  { key: 'Lc', label: 'Characteristic obstacle / channel size', unit: 'm', value: 1e-3, min: 1e-7, max: 10, log: true, group: 'Operation', showIf: { mode: ['particle'] }, help: 'Length over which the flow changes direction (bend radius, fibre diameter, nozzle size).' },
  { key: 'gradT', label: 'Temperature gradient', unit: 'K/m', value: 1e4, min: 0, max: 1e8, group: 'Operation', showIf: { mode: ['particle'] }, help: 'Magnitude of the temperature gradient driving thermophoresis (10⁴ K/m = 10 K per mm).' },
  { key: 'DT', label: 'Thermophoretic mobility D_T', unit: 'm²/s/K', value: 1e-11, min: 0, max: 1e-6, group: 'Particles', showIf: { mode: ['particle'] }, help: 'Thermophoretic velocity per unit temperature gradient. Colloids in liquids: 10⁻¹²–10⁻¹¹ m²/s/K. Nanoparticles in a gas (free-molecular, Waldmann): about 0.55·ν/T, e.g. 3×10⁻⁸ in air at 300 K.' },
  { key: 'lam', label: 'Gas mean free path (0 for liquids)', unit: 'm', value: 0, min: 0, max: 1e-5, group: 'Fluid', showIf: { mode: ['particle'] }, help: 'Used for the Cunningham slip correction of drag on small particles in gases (6.8×10⁻⁸ m for air at ambient conditions). Leave at zero for liquids.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
const num = (inp, k) => { const d = DEF[k]; let v = +inp[k]; if (!Number.isFinite(v)) v = d.value; return N.clamp(v, d.min, d.max); };
const sel = (inp, k) => (DEF[k].options.some(o => o.value === inp[k]) ? inp[k] : DEF[k].value);

// ---------------------------------------------------------------- packed bed
function runBed(inp) {
  const rho = num(inp, 'rho'), mu = num(inp, 'mu'), D = num(inp, 'Dmol'), k = num(inp, 'k1'), Q = num(inp, 'Q'), e = num(inp, 'por'), d = num(inp, 'dg'), L = num(inp, 'Lbed'), A = num(inp, 'Abed'), warnings = [], checks = [];
  const U = Q / A, ui = U / e, Kkc = e ** 3 * d * d / (180 * (1 - e) ** 2), Kerg = e ** 3 * d * d / (150 * (1 - e) ** 2), beta = 1.75 * (1 - e) / (e ** 3 * d);
  const gradOf = u => mu * u / Kerg + beta * rho * u * u, visc = mu * U / Kerg, inert = beta * rho * U * U, dP = (visc + inert) * L, dPdarcy = mu * U / Kkc * L;
  const Rep = rho * U * d / mu, Rem = Rep / (1 - e), regime = Rep < 1 ? 'Darcy (creeping) flow' : Rep < 150 ? 'Forchheimer (steady inertial) flow' : Rep < 300 ? 'unsteady laminar flow' : 'turbulent flow';
  const tort = Math.pow(e, -0.5), Deff = e * D / tort, Dp = Deff / e, PeD = ui * d / D, DL = Dp + 0.5 * ui * d, PeL = ui * L / DL, tau = L / ui, Da = k * tau, X = danckwerts(PeL, Da), Xpfr = 1 - Math.exp(-Da);
  const brink = Math.sqrt(Kkc), DaII = k * d * d / Math.max(Deff, 1e-300), power = dP * Q, fo = inert / (visc + inert);
  checks.push(chk('Continuum pore flow (grain diameter ≥ 100 nm in liquid)', d >= 1e-7, `grain diameter ${F(d, 3)} m, Brinkman length ${F(brink, 3)} m`));
  checks.push(chk('Ergun correlation range (modified Reynolds number < 1000)', Rem < 1000, `Re_p/(1−ε) = ${F(Rem, 3)}`));
  checks.push(chk('Porosity in the range of random packings (0.3–0.6) for which Ergun constants were fitted', e >= 0.3 && e <= 0.6, `ε = ${F(e, 3)}`));
  checks.push(chk('Bed long compared with the grain size (L/d > 50) for a continuum description', L / d > 50, `L/d = ${F(L / d, 3)}`));
  checks.push(chk('Conversion bounded by plug flow and complete mixing', X <= Xpfr + 1e-9 && X >= Da / (1 + Da) - 1e-9, `X = ${F(X, 4)}, plug flow ${F(Xpfr, 4)}, mixed ${F(Da / (1 + Da), 4)}`));
  if (e < 0.3 || e > 0.6) warnings.push('Porosity is outside 0.3–0.6: Kozeny–Carman and Ergun constants (180, 150, 1.75) were fitted to random packings of near-spherical grains and are only indicative for consolidated or very open media.');
  if (d < 1e-7) warnings.push('Sub-100 nm grains: slip, electroviscous and surface-force effects make the continuum permeability estimate unreliable.');
  if (dP > 1e7) warnings.push(`The predicted pressure drop of ${F(dP / 1e5, 3)} bar is beyond normal laboratory equipment; the bed would compact or the pump stall.`);
  warnings.push('One-dimensional uniform bed: wall channelling, compaction, and clogging by deposited product (permeability decline) are not modelled.');
  const us = N.logspace(Math.max(U / 100, 1e-9), U * 100, 80), es = N.linspace(0.1, 0.9, 80);
  return {
    summary: `At a superficial velocity of ${F(U * 1e3, 3)} mm/s through ${F(d * 1e6, 3)} µm grains at porosity ${F(e, 3)}, the Ergun pressure drop over ${F(L, 3)} m is ${F(dP / 1e3, 4)} kPa (${F(100 * fo, 3)} % inertial) and the Kozeny–Carman permeability is ${F(Kkc, 3)} m². The pore Reynolds number is ${F(Rep, 3)}: ${regime}. With a residence time of ${F(tau, 3)} s (Da = ${F(Da, 3)}, Pe = ${F(PeL, 3)}) the first-order conversion is ${F(100 * X, 4)} %.`,
    kpis: [
      kpi('dP', 'Pressure drop (Ergun)', dP / 1e3, 'kPa', 'Darcy–Forchheimer form of the Ergun equation, ΔP/L = 150μ(1−ε)²U/(ε³d²) + 1.75ρ(1−ε)U²/(ε³d). It sets pump head and the mechanical load on the bed.', dP < 1e6),
      kpi('K', 'Permeability (Kozeny–Carman)', Kkc, 'm²', `K = ε³d²/(180(1−ε)²) = ${F(Kkc / 9.869e-13, 3)} darcy; published as the effective permeability for other suites.`),
      kpi('beta', 'Forchheimer coefficient β', beta, '1/m', 'Inertial resistance coefficient 1.75(1−ε)/(ε³d) in ΔP/L = μU/K + βρU².'),
      kpi('Rep', 'Particle Reynolds number', Rep, '', `ρUd/μ. ${regime}; Darcy's law alone holds below about 1.`),
      kpi('ui', 'Interstitial velocity', ui, 'm/s', 'Mean velocity in the pores, U/ε; it sets the residence time.'),
      kpi('Deff', 'Effective diffusivity (Bruggeman)', Deff, 'm²/s', `D_eff = εD/τ with tortuosity τ = ε^(−½) = ${F(tort, 3)}; the diffusive flux per unit bed area.`),
      kpi('brink', 'Brinkman screening length', brink, 'm', '√K: distance over which the velocity adjusts from a free-fluid or wall boundary to the Darcy value.'),
      kpi('PeL', 'Axial Péclet number', PeL, '', `u_i L/D_L with D_L = D_eff/ε + ½u_i d. Above about 100 the bed behaves as plug flow (grain Péclet number ${F(PeD, 3)}).`, PeL > 100),
      kpi('Da', 'Damköhler number', Da, '', 'kτ: reaction rate over throughflow rate. Above about 3 the precursor is largely consumed before the outlet.'),
      kpi('X', 'Precursor conversion', X, '', 'First-order conversion from the axial-dispersion model with closed (Danckwerts) boundaries.'),
      kpi('DaII', 'Second Damköhler number (grain scale)', DaII, '', 'k d²/D_eff: reaction against diffusion on the grain scale. Above one, deposition is limited by transport into the pores and is non-uniform.', DaII < 1),
      kpi('power', 'Hydraulic power', power, 'W', 'ΔP·Q, the ideal pumping power dissipated in the bed.'),
    ],
    series: [
      { key: 'dpu', title: 'Pressure gradient vs superficial velocity', xlabel: 'superficial velocity (m/s)', ylabel: 'pressure gradient (Pa/m)', logx: true, logy: true, lines: [{ name: 'Ergun (Darcy–Forchheimer)', x: us, y: us.map(gradOf) }, { name: 'Darcy with Kozeny–Carman', x: us, y: us.map(u => mu * u / Kkc) }, { name: 'operating point', x: [U, U], y: [gradOf(U) / 2, gradOf(U) * 2] }], explain: 'Linear Darcy behaviour at low velocity with the quadratic Forchheimer term taking over at high velocity. The small constant offset reflects the 150 versus 180 coefficients of the two correlations.' },
      { key: 'ke', title: 'Permeability vs porosity', xlabel: 'porosity', ylabel: 'permeability (m²)', logy: true, lines: [{ name: 'Kozeny–Carman', x: es, y: es.map(x => x ** 3 * d * d / (180 * (1 - x) ** 2)) }, { name: 'operating point', x: [e, e], y: [Kkc / 2, Kkc * 2] }], explain: 'Steep ε³/(1−ε)² dependence: partial pore filling by deposited product lowers permeability far faster than it lowers porosity.' },
    ],
    equations: [...ZNO_EQ(),
      eqn('governing', 'Darcy–Forchheimer (Ergun) momentum balance', '−dp/dx = (μ/K) U + β ρ U²,   K = ε³d²/[150(1−ε)²],   β = 1.75(1−ε)/(ε³d)'),
      eqn('constitutive', 'Kozeny–Carman permeability', 'K = ε³ d² / [180 (1 − ε)²]'),
      eqn('constitutive', 'Bruggeman effective diffusivity', 'D_eff = ε D / τ,   τ = ε^(−1/2)'),
      eqn('governing', 'Axial dispersion with first-order reaction', 'u_i dc/dx = D_L d²c/dx² − k c,   X = 1 − 4a e^(Pe/2) / [(1+a)² e^(aPe/2) − (1−a)² e^(−aPe/2)],  a = (1 + 4Da/Pe)^½')],
    checks, warnings,
    recommendations: [
      rec('high', dP > 5e5 ? `Pressure drop is ${F(dP / 1e5, 3)} bar. Use larger grains (ΔP ∝ d⁻²: doubling d cuts it fourfold), a shorter and wider bed at the same volume, or raise porosity; ${F(100 * fo, 2)} % of the loss is inertial.` : `Pressure drop of ${F(dP / 1e3, 3)} kPa is moderate; there is hydraulic margin to raise throughput by a factor of about ${F(Math.max(1, 5e5 / Math.max(dP, 1e-9)), 2)} before reaching 5 bar.`, 'Ergun pressure drop.'),
      rec('medium', DaII > 1 ? `Grain-scale Damköhler number is ${F(DaII, 3)}: deposition is diffusion limited and concentrates near grain surfaces and the inlet. Use smaller grains or a lower rate (temperature) for uniform coating.` : `Grain-scale Damköhler number is ${F(DaII, 3)}: reaction is slower than pore diffusion, so deposition is uniform inside the pore space.`, 'Second Damköhler number.'),
      rec('medium', Da > 3 ? `The precursor is ${F(100 * X, 3)} % consumed along the bed, giving a strong inlet-to-outlet gradient of deposit. Raise the flow rate, reverse the flow periodically or lower the temperature for uniform product.` : `Conversion per pass is ${F(100 * X, 3)} %; recycle the outlet stream to use the remaining precursor rather than discarding it.`, 'Damköhler number and conversion.'),
      rec('low', `Hydraulic power is ${F(power, 3)} W (${F(power / Math.max(Q, 1e-300) / 3.6e6, 3)} kWh per m³ pumped). Operating in the Darcy regime (Re_p < 1) avoids the quadratic inertial loss and minimises pumping energy per unit throughput.`, 'Split of viscous and inertial pressure loss.', true),
    ],
    publish: { 'eff.permeability': Kkc, 'flow.dP': dP, 'flow.velocity': ui, 'rxn.conversion': X, 'geom.porosity': e, 'proc.flow': Q, 'flow.D_eff': Deff, 'proc.time': tau },
    numerics: { method: 'Closed-form Ergun, Kozeny–Carman and Bruggeman relations; analytical Danckwerts solution of the axial-dispersion model', dof: 1, steps: us.length, iterations: 1, residual: 0, notes: 'Spherical grains assumed (sphericity 1).' },
  };
}

// ---------------------------------------------------------------- 2-D heterogeneous Darcy flow + reactive transport
function runDarcy2D(inp) {
  const mu = num(inp, 'mu'), D = num(inp, 'Dmol'), k = num(inp, 'k1'), e = num(inp, 'por'), Lx = num(inp, 'Lx'), Ly = num(inp, 'Ly'), th = num(inp, 'thick'), Kg = num(inp, 'Kg'), sig = num(inp, 'sigK'), corr = Math.round(num(inp, 'corr')), dP = num(inp, 'dP'), pout = num(inp, 'pout'), cin = num(inp, 'cin'), seed = Math.round(num(inp, 'seed'));
  const [nx, ny] = { coarse: [32, 24], medium: [48, 32], fine: [72, 48] }[sel(inp, 'grid')], n = nx * ny, dx = Lx / nx, dy = Ly / ny, warnings = [], checks = [], id = (i, j) => j * nx + i;
  // correlated log-normal permeability field (seeded). Each smoothing pass is a separable three-point filter whose weights are scaled
  // with the cell size, so that the correlation length is the same physical distance along x and y (statistically isotropic field)
  // even though the cells are not square; 'corr' passes give a Gaussian-like correlation of about √(2·corr/3) cells of the finer spacing.
  const hmin = Math.min(dx, dy), wx = (hmin / dx) ** 2 / 3, wy = (hmin / dy) ** 2 / 3;
  const field = sd0 => {
    const rg = N.rng(sd0); let z = new Float64Array(n); for (let q = 0; q < n; q++) z[q] = rg.normal();
    for (let p = 0; p < corr; p++) {
      const w = new Float64Array(n), w2 = new Float64Array(n);
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { let s = (1 - 2 * wx) * z[id(i, j)], c = 1 - 2 * wx; if (i > 0) { s += wx * z[id(i - 1, j)]; c += wx; } if (i < nx - 1) { s += wx * z[id(i + 1, j)]; c += wx; } w[id(i, j)] = s / c; }
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { let s = (1 - 2 * wy) * w[id(i, j)], c = 1 - 2 * wy; if (j > 0) { s += wy * w[id(i, j - 1)]; c += wy; } if (j < ny - 1) { s += wy * w[id(i, j + 1)]; c += wy; } w2[id(i, j)] = s / c; }
      z = w2;
    }
    let m = 0, v = 0; for (let q = 0; q < n; q++) m += z[q]; m /= n; for (let q = 0; q < n; q++) v += (z[q] - m) ** 2; const sd = Math.sqrt(v / n) || 1;
    const K = new Float64Array(n); for (let q = 0; q < n; q++) K[q] = Kg * Math.exp(sig * (z[q] - m) / sd);
    return K;
  };
  const ax = dy * th / (dx * mu), ay = dx * th / (dy * mu), hm = (a, b) => 2 * a * b / (a + b), Tref = Kg * ax;
  // pressure solution of one realisation: returns the transmissibilities, the scaled pressure and the solver record
  const flow = (K, tol = 1e-10) => {
  const Tx = new Float64Array((nx + 1) * ny), Ty = new Float64Array(nx * (ny + 1));   // transmissibilities (harmonic means)
  for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) Tx[j * (nx + 1) + i] = i === 0 ? 2 * K[id(0, j)] * ax : i === nx ? 2 * K[id(nx - 1, j)] * ax : hm(K[id(i - 1, j)], K[id(i, j)]) * ax;
  for (let j = 1; j < ny; j++) for (let i = 0; i < nx; i++) Ty[j * nx + i] = hm(K[id(i, j - 1)], K[id(i, j)]) * ay;
  const diag = new Float64Array(n), b = new Float64Array(n);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const q = id(i, j); diag[q] = (Tx[j * (nx + 1) + i] + Tx[j * (nx + 1) + i + 1] + Ty[j * nx + i] + Ty[(j + 1) * nx + i]) / Tref; if (i === 0) b[q] = Tx[j * (nx + 1)] / Tref; }
  // Conjugate gradient with a line (tridiagonal block) preconditioner along the more strongly coupled grid direction.
  // The kernel cg has diagonal preconditioning only, which stalls on strongly heterogeneous or stretched grids.
  const mv = (p, o) => { for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const q = id(i, j); let s = diag[q] * p[q]; if (i > 0) s -= Tx[j * (nx + 1) + i] / Tref * p[q - 1]; if (i < nx - 1) s -= Tx[j * (nx + 1) + i + 1] / Tref * p[q + 1]; if (j > 0) s -= Ty[j * nx + i] / Tref * p[q - nx]; if (j < ny - 1) s -= Ty[(j + 1) * nx + i] / Tref * p[q + nx]; o[q] = s; } };
  const alongX = ax >= ay, mL = alongX ? nx : ny, cpL = new Float64Array(mL), dpL = new Float64Array(mL);
  const prec = (r, z) => { const lines = alongX ? ny : nx; for (let l = 0; l < lines; l++) { for (let k2 = 0; k2 < mL; k2++) { const q = alongX ? id(k2, l) : id(l, k2), lo = k2 === 0 ? 0 : (alongX ? Tx[l * (nx + 1) + k2] : Ty[k2 * nx + l]) / Tref, up = k2 === mL - 1 ? 0 : (alongX ? Tx[l * (nx + 1) + k2 + 1] : Ty[(k2 + 1) * nx + l]) / Tref, mm = diag[q] + (k2 ? lo * cpL[k2 - 1] : 0); cpL[k2] = -up / mm; dpL[k2] = (r[q] + (k2 ? lo * dpL[k2 - 1] : 0)) / mm; } let zn = 0; for (let k2 = mL - 1; k2 >= 0; k2--) { const q = alongX ? id(k2, l) : id(l, k2); zn = dpL[k2] - cpL[k2] * zn; z[q] = zn; } } };
  const sol = (() => { const x = new Float64Array(n), r = new Float64Array(n), z = new Float64Array(n), pv = new Float64Array(n), Ap = new Float64Array(n), history = []; for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) x[id(i, j)] = 1 - (i + 0.5) / nx;
    let bn = 0, rn = 0, rz = 0, it = 0; mv(x, Ap); for (let q = 0; q < n; q++) { r[q] = b[q] - Ap[q]; bn += b[q] * b[q]; rn += r[q] * r[q]; } bn = Math.sqrt(bn) || 1; prec(r, z); for (let q = 0; q < n; q++) { pv[q] = z[q]; rz += r[q] * z[q]; }
    let res = Math.sqrt(rn) / bn; history.push(res); const maxIter = 3000;
    while (res > tol && it < maxIter) { mv(pv, Ap); let pAp = 0; for (let q = 0; q < n; q++) pAp += pv[q] * Ap[q]; const al = rz / pAp; rn = 0; for (let q = 0; q < n; q++) { x[q] += al * pv[q]; r[q] -= al * Ap[q]; rn += r[q] * r[q]; } prec(r, z); let rz2 = 0; for (let q = 0; q < n; q++) rz2 += r[q] * z[q]; const be = rz2 / rz; rz = rz2; for (let q = 0; q < n; q++) pv[q] = z[q] + be * pv[q]; res = Math.sqrt(rn) / bn; history.push(res); it++; }
    return { x, iterations: it, residual: res, history, converged: res <= tol }; })();
    let Qi = 0, Qo = 0; for (let j = 0; j < ny; j++) { Qi += Tx[j * (nx + 1)] * (1 - sol.x[id(0, j)]) * dP; Qo += Tx[j * (nx + 1) + nx] * sol.x[id(nx - 1, j)] * dP; }
    return { Tx, Ty, sol, Keff: 0.5 * (Qi + Qo) * mu * Lx / (Ly * th * dP) };
  };
  const K = field(seed), { Tx, Ty, sol, Keff: Keff0 } = flow(K), ph = sol.x;   // ph = (p − p_out)/ΔP
  // ensemble of independent seeded realisations (pressure solution only, tolerance 10⁻⁶): mean ratio K_eff/K_g and its standard error.
  // The ensemble stops early once a fixed work budget (conjugate-gradient iterations × cells) is spent, so the run time stays bounded.
  const nreq = Math.round(N.clamp(Number(inp.nens) || 1, 1, 16)), k0 = Keff0 / Kg, ensR = [k0]; let work = 0;
  for (let q = 1; q < nreq && work < 3e6; q++) { if (sig > 0) { const f = flow(field(seed + 10007 * q), 1e-6); work += (f.sol.iterations + 1) * n; ensR.push(f.Keff / Kg); } else ensR.push(k0); }
  const nens = ensR.length, ensMean = N.sum(ensR) / nens, ensSE = nens > 1 ? Math.sqrt(N.sum(ensR.map(r => (r - ensMean) ** 2)) / (nens - 1) / nens) : 0;
  if (nens < nreq) warnings.push(`Only ${nens} of the ${nreq} requested realisations were solved: the pressure solver needs many iterations for this field and the ensemble was stopped at its work budget. Use a coarser grid or a smoother field for a larger ensemble.`);
  // face fluxes (m³/s), positive in +x / +y
  const Fx = new Float64Array((nx + 1) * ny), Fy = new Float64Array(nx * (ny + 1));
  for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) { const pl = i === 0 ? 1 : ph[id(i - 1, j)], pr = i === nx ? 0 : ph[id(i, j)]; Fx[j * (nx + 1) + i] = Tx[j * (nx + 1) + i] * (pl - pr) * dP; }
  for (let j = 1; j < ny; j++) for (let i = 0; i < nx; i++) Fy[j * nx + i] = Ty[j * nx + i] * (ph[id(i, j - 1)] - ph[id(i, j)]) * dP;
  let Qin = 0, Qout = 0; const qo = []; for (let j = 0; j < ny; j++) { Qin += Fx[j * (nx + 1)]; Qout += Fx[j * (nx + 1) + nx]; qo.push(Fx[j * (nx + 1) + nx]); }
  const Qm = 0.5 * (Qin + Qout), mbFlow = Math.abs(Qin - Qout) / Math.max(Math.abs(Qm), 1e-300), Keff = Qm * mu * Lx / (Ly * th * dP), qbar = Qout / ny, unif = 1 - N.sum(qo.map(q => Math.abs(q - qbar))) / (2 * ny * Math.max(Math.abs(qbar), 1e-300));
  let Ka = 0, Kh = 0; for (let q = 0; q < n; q++) { Ka += K[q]; Kh += 1 / K[q]; } Ka /= n; Kh = n / Kh;
  const speed = new Float64Array(n), pmap = new Float64Array(n), lk = new Float64Array(n); let umax = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const q = id(i, j), ux = 0.5 * (Fx[j * (nx + 1) + i] + Fx[j * (nx + 1) + i + 1]) / (dy * th), uy = 0.5 * (Fy[j * nx + i] + Fy[(j + 1) * nx + i]) / (dx * th); speed[q] = Math.hypot(ux, uy); umax = Math.max(umax, speed[q]); pmap[q] = pout + ph[q] * dP; lk[q] = Math.log10(K[q]); }
  const Um = Qm / (Ly * th);
  // steady advection–diffusion–reaction, first-order upwind, column line Gauss–Seidel
  const De = Math.pow(e, 1.5) * D, gx = De * dy * th / dx, gy = De * dx * th / dy, Vc = dx * dy * th, c = new Float64Array(n), aP = new Float64Array(n), sweepsMax = 500, omega = Um * dx / Math.max(De, 1e-300) > 2 ? 1 : 1.6;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const q = id(i, j), fw = Fx[j * (nx + 1) + i], fe = Fx[j * (nx + 1) + i + 1], fs = Fy[j * nx + i], fn = Fy[(j + 1) * nx + i]; aP[q] = Math.max(fe, 0) + Math.max(-fw, 0) + Math.max(fn, 0) + Math.max(-fs, 0) + (i > 0 ? gx : 2 * gx) + (i < nx - 1 ? gx : 0) + (j > 0 ? gy : 0) + (j < ny - 1 ? gy : 0) + e * k * Vc; }
  let sweeps = 0, change = 1; const la = new Float64Array(ny), lb = new Float64Array(ny), lc = new Float64Array(ny), ld = new Float64Array(ny);
  for (; sweeps < sweepsMax && change > 1e-10; sweeps++) {
    change = 0;
    for (let ii = 0; ii < nx; ii++) {
      const i = sweeps % 2 ? nx - 1 - ii : ii;
      for (let j = 0; j < ny; j++) { const q = id(i, j), fw = Fx[j * (nx + 1) + i], fe = Fx[j * (nx + 1) + i + 1], fs = Fy[j * nx + i], fn = Fy[(j + 1) * nx + i]; const cw = i > 0 ? c[q - 1] : 1, aw = Math.max(fw, 0) + (i > 0 ? gx : 2 * gx), ae = Math.max(-fe, 0) + (i < nx - 1 ? gx : 0); la[j] = j > 0 ? -(Math.max(fs, 0) + gy) : 0; lc[j] = j < ny - 1 ? -(Math.max(-fn, 0) + gy) : 0; lb[j] = aP[q]; ld[j] = aw * cw + (i < nx - 1 ? ae * c[q + 1] : 0); }
      const cn = N.solveTridiag(la, lb, lc, ld); for (let j = 0; j < ny; j++) { const q = id(i, j), val = Math.max(c[q] + omega * (cn[j] - c[q]), 0); change = Math.max(change, Math.abs(val - c[q])); c[q] = val; }
    }
  }
  let outC = 0, react = 0, cmin = Infinity, diffIn = 0; for (let j = 0; j < ny; j++) { outC += Fx[j * (nx + 1) + nx] * c[id(nx - 1, j)]; diffIn += 2 * gx * (1 - c[id(0, j)]); } for (let q = 0; q < n; q++) { react += e * k * Vc * c[q]; cmin = Math.min(cmin, c[q]); }
  const inC = Qin + diffIn, mbC = Math.abs(inC - outC - react) / Math.max(inC, 1e-300), cOut = outC / Math.max(Qout, 1e-300), X = N.clamp(1 - cOut, 0, 1), cmap = Float64Array.from(c, vv => vv * cin);
  const tau = e * Lx / Math.max(Um, 1e-300), Da = k * tau, Pe = Um * Lx / Math.max(De, 1e-300), cellPe = Um * dx / Math.max(De, 1e-300), Xpfr = 1 - Math.exp(-Da);
  const prof = Array.from({ length: nx }, (_, i) => { let s = 0; for (let j = 0; j < ny; j++) s += c[id(i, j)]; return s / ny * cin; }), xc = Array.from({ length: nx }, (_, i) => (i + 0.5) * dx);
  checks.push(chk('Pressure solver converged (conjugate gradient)', sol.converged, `${sol.iterations} iterations, relative residual ${F(sol.residual, 2)}`));
  checks.push(chk('Global fluid mass balance (inflow = outflow)', mbFlow < 1e-6, `relative imbalance ${F(mbFlow, 2)}`));
  checks.push(chk('Effective permeability within the Wiener bounds (harmonic ≤ K_eff ≤ arithmetic mean)', Keff >= Kh * (1 - 1e-6) && Keff <= Ka * (1 + 1e-6), `harmonic ${F(Kh, 3)}, effective ${F(Keff, 3)}, arithmetic ${F(Ka, 3)} m²`));
  checks.push(chk('Species balance (inflow = outflow + reacted)', mbC < 1e-5, `relative residual ${F(mbC, 2)} after ${sweeps} sweeps`));
  checks.push(chk('Concentration bounded (0 ≤ c ≤ c_in, positivity-preserving upwind scheme)', cmin >= 0 && N.maxOf(c) <= 1 + 1e-9, `min ${F(cmin * cin, 3)}, max ${F(N.maxOf(c) * cin, 3)} mol/m³`));
  if (cellPe > 2) warnings.push(`Cell Péclet number is ${F(cellPe, 3)}: the first-order upwind scheme adds numerical diffusion of order u·Δx/2 (${F(Um * dx / 2, 2)} m²/s), which smears the concentration front more than physical diffusion does.`);
  if (change > 1e-8) warnings.push(`The transport iteration stopped after ${sweeps} sweeps with a remaining change of ${F(change, 2)} of the inlet concentration (diffusion-dominated case converges slowly); the species balance residual shows the resulting error.`);
  if (sig > 0.5 && corr < 2) warnings.push('The permeability varies from cell to cell (correlation length below about one cell): the harmonic face averaging of the finite-volume scheme then underestimates the effective permeability systematically. Use at least two or three smoothing passes or a finer grid so that each permeability patch spans several cells.');
  if (sig > 2) warnings.push('Strong heterogeneity (σ_lnK > 2): the result depends on the particular realisation; average several seeds and refine the grid before drawing conclusions.');
  if (!sol.converged) warnings.push('The pressure solver did not reach its tolerance; flow results are approximate.');
  warnings.push('Single-phase steady Darcy flow on a statistically generated permeability field (not a measured microstructure); mechanical dispersion and permeability change by deposition are not included.');
  const mk = (key, title, data, unit, explain, diverging = false) => ({ key, title, nx, ny, data, x0: 0, x1: Lx, y0: 0, y1: Ly, xlabel: 'x (m), flow →', ylabel: 'y (m)', unit, explain, diverging });
  return {
    summary: `Darcy flow through the ${F(Lx, 3)} m × ${F(Ly, 3)} m × ${F(th, 3)} m slab under ${F(dP / 1e3, 4)} kPa gives ${F(Qm * 1e6, 3)} mL/s and an effective permeability of ${F(Keff, 3)} m² (${F(Keff / Kg, 3)} times the geometric mean). The outlet flow-uniformity index is ${F(unif, 3)} and the peak velocity is ${F(umax / Math.max(Um, 1e-300), 3)} times the mean. The precursor is ${F(100 * X, 4)} % consumed (outlet ${F(cOut * cin, 3)} mol/m³, minimum ${F(cmin * cin, 3)} mol/m³) at Da = ${F(Da, 3)}.`,
    kpis: [
      kpi('Q', 'Volumetric flow rate', Qm, 'm³/s', 'Sum of the face fluxes through the inlet (equal to the outlet within the solver tolerance); published for reactor and scale-up suites.'),
      kpi('Keff', 'Effective permeability', Keff, 'm²', 'QμL/(AΔP) of the heterogeneous slab. For a 2-D log-normal field theory predicts the geometric mean; departures of this single realisation reflect finite size and sampling (see the ensemble mean).'),
      kpi('Krat', 'K_eff / geometric mean', Keff / Kg, '', 'Ratio to the geometric-mean permeability; values well above one indicate a connected high-permeability channel across the domain.'),
      kpi('Krat_ens', 'Ensemble-mean K_eff / geometric mean', ensMean, '', `Mean of K_eff/K_g over ${nens} independent seeded realisations of the same statistics (seeds ${seed}, ${seed} + 10007, …). For a statistically isotropic two-dimensional log-normal field the exact infinite-medium result is one; the harmonic face averaging of the finite-volume scheme lowers it by a few per cent when the correlation length is below about two cells.`),
      kpi('Krat_se', 'Standard error of the ensemble mean', ensSE, '', nens > 1 ? `Sample standard deviation of the ${nens} ratios divided by √${nens}: the sampling uncertainty of the ensemble mean. A single realisation scatters by about ${F(ensSE * Math.sqrt(nens), 2)}.` : 'Only one realisation was requested, so no sampling uncertainty can be estimated; raise the number of realisations.'),
      kpi('unif', 'Outlet flow-uniformity index', unif, '', '1 − Σ|q_j − q̄|/(2Σq_j): one means perfectly uniform outflow, lower values mean channelling.', unif > 0.9 ? true : unif < 0.7 ? false : null),
      kpi('chan', 'Peak-to-mean velocity ratio', umax / Math.max(Um, 1e-300), '', 'Largest cell velocity divided by the mean Darcy velocity: a measure of preferential flow paths.'),
      kpi('U', 'Mean Darcy velocity', Um, 'm/s', 'Flow rate per unit cross-section; the mean pore velocity is this value divided by the porosity.'),
      kpi('X', 'Precursor conversion', X, '', `1 − flux-weighted outlet concentration over inlet concentration. Plug flow at the same Damköhler number would give ${F(Xpfr, 4)}; channelling lowers it.`),
      kpi('cmin', 'Minimum precursor concentration', cmin * cin, 'mol/m³', 'Lowest concentration in the domain: stagnant low-permeability pockets are starved of precursor and receive little deposit.'),
      kpi('Da', 'Damköhler number', Da, '', 'k·εL/U: reaction rate relative to throughflow.'),
      kpi('Pe', 'Péclet number', Pe, '', 'UL/D_eff: advection relative to diffusion over the domain length.'),
      kpi('mb', 'Species balance residual', mbC, '', 'Relative mismatch between inflow, outflow and reacted amount: a conservation check of the finite-volume solution.', mbC < 1e-5),
    ],
    series: [
      { key: 'prof', title: 'Width-averaged concentration along the flow', xlabel: 'x (m)', ylabel: 'concentration (mol/m³)', lines: [{ name: '2-D heterogeneous', x: xc, y: prof }, { name: 'plug flow', x: xc, y: xc.map(x => cin * Math.exp(-Da * x / Lx)) }], explain: 'Depletion of the precursor from inlet to outlet compared with ideal plug flow at the same mean velocity.' },
      { key: 'out', title: 'Outlet flux distribution', xlabel: 'y (m)', ylabel: 'local outflow / mean', lines: [{ name: 'outlet flux', x: qo.map((_, j) => (j + 0.5) * dy), y: qo.map(q => q / Math.max(Math.abs(qbar), 1e-300)) }, { name: 'uniform', x: [0, Ly], y: [1, 1] }], explain: 'Normalised outflow across the outlet face; peaks mark preferential channels.' },
      { key: 'cg', title: 'Pressure-solver convergence', xlabel: 'iteration', ylabel: 'relative residual', logy: true, lines: [{ name: 'conjugate gradient', x: N.linspace(0, sol.iterations, Math.min(sol.history.length, 300)).map(Math.round), y: N.linspace(0, sol.history.length - 1, Math.min(sol.history.length, 300)).map(q => Math.max(sol.history[Math.round(q)], 1e-300)) }], explain: 'Residual history of the line-preconditioned conjugate-gradient solve; strong heterogeneity slows convergence.' },
    ],
    maps: [
      mk('lk', 'Permeability field', lk, 'log₁₀(m²)', 'Seeded log-normal permeability with the chosen variance and correlation; bright patches conduct flow preferentially.'),
      mk('p', 'Pressure', pmap, 'Pa', 'Pressure falls from the inlet (left) to the outlet (right); contours bend around low-permeability patches.'),
      mk('u', 'Darcy velocity magnitude', speed, 'm/s', 'Cell-centred magnitude of the Darcy flux; streaks of high velocity are the preferential flow channels.'),
      mk('c', 'Precursor concentration', cmap, 'mol/m³', 'Steady concentration with first-order consumption; fingers of high concentration follow fast channels while slow zones are depleted.'),
    ],
    equations: [...ZNO_EQ(),
      eqn('governing', 'Darcy flow (mass conservation)', '∇·( (K/μ) ∇p ) = 0,   u = −(K/μ) ∇p'),
      eqn('governing', 'Steady advection–diffusion–reaction', '∇·(u c) − ∇·(D_eff ∇c) + ε k c = 0,   D_eff = ε^(3/2) D'),
      eqn('constitutive', 'Permeability field', 'ln K ~ N(ln K_g, σ²), spatially correlated by box smoothing'),
      eqn('formula', 'Effective permeability and uniformity', 'K_eff = Q μ L_x / (L_y t Δp);   γ = 1 − Σ|q_j − q̄| / (2 Σ q_j)')],
    checks, warnings,
    recommendations: [
      rec('high', unif < 0.8 ? `Flow is channelled (uniformity ${F(unif, 3)}, peak velocity ${F(umax / Math.max(Um, 1e-300), 3)}× mean). Improve packing homogeneity, add a distributor plate or a high-resistance layer at the inlet/outlet so the bed resistance no longer controls the distribution.` : `Flow distribution is acceptable (uniformity ${F(unif, 3)}); keep the permeability variance at or below the present σ_lnK = ${F(sig, 2)} when scaling up.`, 'Outlet flow-uniformity index.'),
      rec('medium', Da > 1 ? `Da = ${F(Da, 3)}: the precursor is depleted along the flow path (minimum ${F(cmin * cin, 3)} of ${F(cin, 3)} mol/m³). Raise the pressure difference, shorten the flow path or lower the temperature to keep the deposition uniform.` : `Da = ${F(Da, 3)}: reaction is slow compared with throughflow, so deposition is nearly uniform but only ${F(100 * X, 3)} % of the precursor is used per pass; recirculate the outlet stream.`, 'Damköhler number and concentration map.'),
      rec('medium', `Check discretisation by repeating on the finer grid and with several seeds; the first-order upwind scheme is bounded but diffusive (cell Péclet number ${F(cellPe, 3)}).`, 'Numerical diffusion of the upwind scheme.'),
      rec('low', `Pumping power is ${F(dP * Qm, 3)} W. Since flow scales linearly with pressure in the Darcy regime, a wider, shorter bed delivers the same residence time at a fraction of the pressure and pumping energy, and recirculation avoids discharging unreacted precursor.`, 'Darcy scaling of flow with pressure gradient.', true),
    ],
    publish: { 'eff.permeability': Keff, 'proc.flow': Qm, 'rxn.conversion': X, 'flow.uniformity': unif, 'flow.velocity': Um / e, 'flow.dP': dP, 'geom.porosity': e, 'proc.time': tau },
    numerics: { method: 'Cell-centred finite volumes with harmonic transmissibilities; conjugate gradient with a line (tridiagonal) preconditioner for pressure; first-order upwind transport solved by alternating-direction column line Gauss–Seidel (tridiagonal solves, over-relaxed when diffusion dominates)', dof: n, steps: sweeps, iterations: sol.iterations, residual: sol.residual, notes: `Grid ${nx} × ${ny}; transport stopped at a maximum change of ${F(change, 2)} of the inlet concentration.` },
  };
}

// ---------------------------------------------------------------- channel / pipe / coil
function runChannel(inp) {
  const rho = num(inp, 'rho'), mu = num(inp, 'mu'), D = num(inp, 'Dmol'), k = num(inp, 'k1'), Q = num(inp, 'Q'), geo = sel(inp, 'geo'), d = num(inp, 'dh'), L = num(inp, 'Lch'), rough = num(inp, 'rough'), Dc = Math.max(num(inp, 'Dc'), 2 * d), warnings = [], checks = [];
  const A = PI * d * d / 4, U = Q / A, Re = rho * U * d / mu, tau = L / U, coil = geo === 'coil', De = coil ? Re * Math.sqrt(d / Dc) : 0, ReCrit = coil ? 2100 * (1 + 12 * Math.sqrt(d / Dc)) : 2300;
  const fTurb = re => { const rr = rough / d; let f = Math.pow(-1.8 * Math.log10(Math.pow(rr / 3.7, 1.11) + 6.9 / re), -2); for (let i = 0; i < 30; i++) f = Math.pow(-2 * Math.log10(rr / 3.7 + 2.51 / (re * Math.sqrt(f))), -2); return f; };
  const coilFac = coil && De > 11.6 ? 1 + 0.033 * Math.pow(Math.log10(De), 4) : 1, fLam = re => 64 / re * coilFac;
  const hiT = ReCrit * 4000 / 2300, regime = Re < ReCrit ? 'laminar' : Re < hiT ? 'transitional' : 'turbulent';
  const f = regime === 'laminar' ? fLam(Re) : regime === 'turbulent' ? fTurb(Re) : fLam(ReCrit) + (fTurb(hiT) - fLam(ReCrit)) * (Re - ReCrit) / (hiT - ReCrit), fHaal = Math.pow(-1.8 * Math.log10(Math.pow(rough / d / 3.7, 1.11) + 6.9 / Math.max(Re, 1)), -2);
  const dP = f * (L / d) * rho * U * U / 2, tw = f * rho * U * U / 8, Le = regime === 'laminar' ? 0.06 * Re * d : 4.4 * Math.pow(Re, 1 / 6) * d, power = dP * Q, gw = 8 * U / d;
  const Fo = D * tau / (d * d / 4), taylorOK = regime !== 'laminar' || Fo >= 0.8;
  const Dax = regime === 'laminar' ? D + U * U * d * d / (192 * D) : 10.1 * (d / 2) * U * Math.sqrt(f / 8) + D, Pe = U * L / Dax, s2 = 2 / Pe - 2 / (Pe * Pe) * (1 - Math.exp(-Math.min(Pe, 700))), Nt = 1 / Math.max(s2, 1e-12), Da = k * tau;
  const Xpfr = 1 - Math.exp(-Da), Xcstr = Da / (1 + Da), Xdisp = danckwerts(Pe, Da), Xtis = 1 - Math.pow(1 + Da / Nt, -Nt);
  const lamX = da => { let s = 0; const m = 400; for (let i = 0; i <= m; i++) { const u = i / m, w = i === 0 || i === m ? 1 : i % 2 ? 4 : 2; s += w * (u > 0 ? 2 * u * Math.exp(-Math.min(da / (2 * u), 700)) : 0); } return 1 - s / (3 * m); };
  const Xlam = lamX(Da), model = regime === 'laminar' ? (taylorOK ? 'Taylor–Aris dispersion' : 'laminar convection (segregated flow)') : 'turbulent axial dispersion', X = model.startsWith('laminar convection') ? Xlam : Xdisp;
  const th = N.linspace(0.02, 3, 200), Etis = th.map(t => Math.exp(Math.min(Math.log(Nt) + (Nt - 1) * Math.log(Nt * t) - Nt * t - N.gammaLn(Nt), 50)));
  checks.push(chk('Flow regime unambiguous (not transitional)', regime !== 'transitional', `Re = ${F(Re, 4)}, critical ${F(ReCrit, 4)}`));
  checks.push(chk('Tube long compared with the hydrodynamic entrance length', L > 10 * Le || regime !== 'laminar' ? L > Le : false, `entrance length ${F(Le, 3)} m vs tube length ${F(L, 3)} m`));
  checks.push(chk('Taylor–Aris dispersion regime reached (Dτ/a² ≥ 0.8)', taylorOK, regime === 'laminar' ? `Dτ/a² = ${F(Fo, 3)}` : 'turbulent or transitional flow: turbulent dispersion correlation used'));
  checks.push(chk('Conversion between the mixed-tank and plug-flow limits', X <= Xpfr + 1e-9 && X >= Math.min(Xcstr, Xlam) - 1e-9, `X = ${F(X, 4)} (plug ${F(Xpfr, 4)}, mixed ${F(Xcstr, 4)})`));
  checks.push(chk('Incompressible-flow assumption (pressure drop small against 1 bar for gases, any for liquids)', rho > 50 || dP < 2e4, `ΔP = ${F(dP, 3)} Pa, density ${F(rho, 3)} kg/m³`));
  if (regime === 'transitional') warnings.push('Transitional Reynolds number: friction factor is interpolated between laminar and turbulent values and is uncertain by tens of per cent; avoid operating here.');
  if (!taylorOK) warnings.push(`Radial diffusion is too slow to average the parabolic profile within the residence time (Dτ/a² = ${F(Fo, 3)}): the residence-time distribution is the broad laminar-convection form E(θ) = 1/(2θ³), and conversion is computed for segregated laminar flow. Taylor–Aris values are shown for reference only.`);
  if (coil) warnings.push('In coils the secondary Dean flow narrows the residence-time distribution relative to a straight tube; this beneficial effect is not included in the dispersion estimate, which is therefore conservative.');
  if (rho <= 50 && dP >= 2e4) warnings.push('Gas flow with a pressure drop comparable to the absolute pressure: compressibility changes density and velocity along the tube.');
  warnings.push('Fully developed Newtonian single-phase flow; nanoparticle suspensions above a few volume per cent are shear-thinning and wall deposition (fouling) narrows the channel over time.');
  return {
    summary: `Flow of ${F(Q * 1e6, 3)} mL/s in a ${F(d * 1e3, 3)} mm ${coil ? 'coiled' : 'straight'} tube is ${regime} (Re = ${F(Re, 4)}${coil ? `, Dean number ${F(De, 3)}` : ''}) with a pressure drop of ${F(dP / 1e3, 4)} kPa over ${F(L, 3)} m and a mean residence time of ${F(tau, 3)} s. Axial mixing follows ${model} (Pe = ${F(Pe, 3)}, equivalent to ${F(Nt, 3)} tanks in series); first-order conversion is ${F(100 * X, 4)} % against ${F(100 * Xpfr, 4)} % for ideal plug flow.`,
    kpis: [
      kpi('Re', 'Reynolds number', Re, '', `ρUd/μ. ${regime[0].toUpperCase() + regime.slice(1)} flow; the transition starts near ${F(ReCrit, 4)} for this geometry.`),
      kpi('dP', 'Pressure drop', dP / 1e3, 'kPa', `Darcy–Weisbach ΔP = f(L/d)ρU²/2 with f = ${F(f, 4)} (${regime === 'laminar' ? 'Hagen–Poiseuille 64/Re' + (coilFac > 1 ? ' × Mishra–Gupta coil factor ' + F(coilFac, 3) : '') : 'Colebrook; Haaland gives ' + F(fHaal, 4)}).`),
      kpi('tau', 'Mean residence time', tau, 's', 'Volume over flow rate, L/U; the first moment of the residence-time distribution.'),
      kpi('U', 'Mean velocity', U, 'm/s', 'Flow rate divided by the cross-section; the centre-line velocity is twice this in laminar flow.'),
      kpi('Dax', 'Axial dispersion coefficient', Dax, 'm²/s', regime === 'laminar' ? 'Taylor–Aris D + U²d²/(192D): shear dispersion from the parabolic profile and radial diffusion.' : 'Taylor turbulent dispersion 10.1·a·u* plus molecular diffusion.'),
      kpi('Pe', 'Axial Péclet (Bodenstein) number', Pe, '', 'UL/D_ax. Above about 100 the reactor is close to plug flow; below about 10 back-mixing is strong.', Pe > 100 ? true : Pe < 10 ? false : null),
      kpi('Nt', 'Equivalent tanks in series', Nt, '', 'Reciprocal of the dimensionless variance of the closed-vessel dispersion model; one is a single stirred tank, infinity is plug flow.'),
      kpi('sig', 'Standard deviation of residence time', tau * Math.sqrt(s2), 's', 'Spread of residence times from the dispersion model, τ√(σ_θ²); particles formed in a wide distribution of residence times have a wide size distribution.'),
      kpi('X', 'First-order conversion', X, '', `From ${model}. Tanks-in-series gives ${F(Xtis, 4)}, laminar segregated flow ${F(Xlam, 4)}, plug flow ${F(Xpfr, 4)}.`),
      kpi('tw', 'Wall shear stress', tw, 'Pa', `fρU²/8, with a nominal wall shear rate 8U/d = ${F(gw, 3)} 1/s; relevant to deposit removal and to shear-sensitive agglomerates.`),
      kpi('Le', 'Hydrodynamic entrance length', Le, 'm', regime === 'laminar' ? '0.06·Re·d for laminar flow; beyond it the velocity profile is fully developed.' : '4.4·Re^(1/6)·d for turbulent flow.'),
      kpi('De', 'Dean number', De, '', coil ? 'Re√(d/D_c). Above about 10 secondary vortices enhance radial mixing and heat transfer.' : 'Zero for a straight tube (no secondary flow).'),
      kpi('power', 'Pumping power', power, 'W', 'ΔP·Q, the ideal hydraulic power delivered to the fluid; divide by the pump efficiency for the electrical demand.'),
    ],
    series: [
      { key: 'rtd', title: 'Residence-time distribution', xlabel: 'dimensionless time θ = t/τ', ylabel: 'E(θ)', lines: [{ name: `tanks in series, N = ${F(Nt, 3)}`, x: th, y: Etis }, { name: 'laminar convection', x: th, y: th.map(t => (t >= 0.5 ? 1 / (2 * t ** 3) : 0)) }, { name: 'single stirred tank', x: th, y: th.map(t => Math.exp(-t)) }], explain: 'Distribution of residence times for the dispersion-equivalent tanks-in-series model, for pure laminar convection (first fluid arrives at θ = 0.5) and for a single stirred tank. A narrow peak at θ = 1 is the plug-flow ideal.' },
      { key: 'xda', title: 'Conversion vs Damköhler number', xlabel: 'Da = kτ', ylabel: 'conversion', logx: true, lines: (() => { const das = N.logspace(0.01, 100, 80); return [{ name: 'plug flow', x: das, y: das.map(v => 1 - Math.exp(-v)) }, { name: 'axial dispersion', x: das, y: das.map(v => danckwerts(Pe, v)) }, { name: 'laminar segregated', x: das, y: das.map(lamX) }, { name: 'stirred tank', x: das, y: das.map(v => v / (1 + v)) }]; })(), explain: 'First-order conversion of the reactor models. The gap between the curves at high Damköhler number is the penalty of back-mixing or of the laminar velocity profile.' },
    ],
    equations: [...ZNO_EQ(),
      eqn('governing', 'Darcy–Weisbach / Hagen–Poiseuille', 'Δp = f (L/d) ρ U²/2;   laminar f = 64/Re  ⇔  Q = π d⁴ Δp / (128 μ L)'),
      eqn('constitutive', 'Turbulent friction (Colebrook, Haaland start)', '1/√f = −2 log₁₀[ ε/(3.7d) + 2.51/(Re √f) ];   1/√f ≈ −1.8 log₁₀[ (ε/3.7d)^1.11 + 6.9/Re ]'),
      eqn('constitutive', 'Taylor–Aris dispersion', 'D_ax = D + U² d² / (192 D)   (laminar);   D_ax ≈ 10.1 a u*   (turbulent)'),
      eqn('governing', 'Residence-time distribution', 'E(θ) = N (Nθ)^(N−1) e^(−Nθ) / Γ(N),   σ_θ² = 2/Pe − 2(1 − e^(−Pe))/Pe² = 1/N'),
      eqn('formula', 'Dean number and coil friction', 'De = Re (d/D_c)^½;   f_coil/f_straight = 1 + 0.033 (log₁₀ De)⁴   (laminar, Mishra–Gupta)')],
    checks, warnings,
    recommendations: [
      rec('high', Pe < 50 ? `Back-mixing is significant (Pe = ${F(Pe, 3)}, ${F(Nt, 3)} tanks): expect broad particle-size distributions. ${regime === 'laminar' ? 'Use a smaller diameter (D_ax ∝ d²), a coiled or segmented (gas–liquid slug) flow, or a longer tube at the same residence time.' : 'Lengthen the tube at constant residence time.'}` : `Axial dispersion is small (Pe = ${F(Pe, 3)}); the tube behaves close to plug flow, which favours narrow particle-size distributions.`, 'Bodenstein number.'),
      rec('medium', regime === 'transitional' ? 'Move out of the transitional regime by changing diameter or flow rate; pressure drop and mixing are not reproducible there.' : regime === 'laminar' && !taylorOK ? `Radial mixing is incomplete (Dτ/a² = ${F(Fo, 3)}). Reduce the diameter to below ${F(2 * Math.sqrt(D * tau / 0.8) * 1e3, 3)} mm or add static mixing/coiling so that all fluid elements see the same residence time.` : `Pressure drop of ${F(dP / 1e3, 3)} kPa is ${dP > 5e5 ? 'high: consider a larger diameter (ΔP ∝ d⁻⁴ in laminar flow)' : 'within the range of syringe and HPLC pumps'}.`, 'Flow regime and radial-diffusion criterion.'),
      rec('medium', `Wall shear stress is ${F(tw, 3)} Pa; ${tw < 0.5 ? 'this is low, so particles depositing on the wall are unlikely to be swept away — expect fouling and plan periodic flushing or a segmented flow.' : 'this helps to limit wall deposits.'}`, 'Wall shear from the friction factor.'),
      rec('low', `Continuous flow uses ${F(power, 3)} W of pumping for ${F(Q * 3.6e6, 3)} L/h; number up identical tubes rather than enlarging the diameter to keep the same residence-time distribution, heat transfer and product quality with minimal re-optimisation and waste.`, 'Scaling rules of laminar flow reactors.', true),
    ],
    publish: { 'flow.Re': Re, 'flow.dP': dP, 'rxn.conversion': X, 'proc.time': tau, 'flow.velocity': U, 'proc.flow': Q, 'flow.peclet': Pe, 'proc.volume': A * L },
    numerics: { method: 'Closed-form friction and dispersion correlations; fixed-point Colebrook iteration; Simpson integration of the laminar segregated-flow conversion', dof: 1, steps: 400, iterations: 30, residual: 0, notes: `Critical Reynolds number ${F(ReCrit, 4)}${coil ? ' (Srinivasan correlation for coils)' : ''}.` },
  };
}

// ---------------------------------------------------------------- stirred vessel
const IMP = { rushton: { name: 'Rushton turbine', Np: 5.0, Kp: 70, Nq: 0.72, ks: 11.5 }, pbt: { name: 'pitched-blade turbine', Np: 1.27, Kp: 50, Nq: 0.79, ks: 13 }, hydrofoil: { name: 'hydrofoil', Np: 0.3, Kp: 40, Nq: 0.56, ks: 10 } };
function runStirred(inp) {
  const rho = num(inp, 'rho'), mu = num(inp, 'mu'), im = IMP[sel(inp, 'imp')], Tt = num(inp, 'Tt'), Dm = Math.min(num(inp, 'Dimp'), 0.9 * Tt), Ns = num(inp, 'Nrpm') / 60, dp = num(inp, 'dpS'), rs = num(inp, 'rhoS'), Xs = num(inp, 'Xs'), S = num(inp, 'Szw'), warnings = [], checks = [];
  const nu = mu / rho, Re = rho * Ns * Dm * Dm / mu, Np = im.Kp / Re + im.Np, V = PI * Tt ** 3 / 4, P = Np * rho * Ns ** 3 * Dm ** 5, eps = P / (rho * V), eta = Math.pow(nu ** 3 / eps, 0.25), tip = PI * Ns * Dm;
  const regime = Re < 10 ? 'laminar' : Re < 1e4 ? 'transitional' : 'turbulent', crit = Math.cbrt(Np) * Re;
  const tMix = crit > 6370 ? 5.2 * Tt * Tt / (Math.cbrt(Np) * Ns * Dm * Dm) : rho * Tt * Tt / mu * (183 / crit) ** 2;
  const gMO = im.ks * Ns, gK = Math.sqrt(eps / nu), tauK = Math.sqrt(nu / eps), drho = Math.abs(rs - rho), Njs = S * Math.pow(nu, 0.1) * Math.pow(dp, 0.2) * Math.pow(G * drho / rho, 0.45) * Math.pow(Xs, 0.13) / Math.pow(Dm, 0.85);
  const Fr = Ns * Ns * Dm / G, tCirc = V / (im.Nq * Ns * Dm ** 3), stress = mu * gK, PV = P / V, PofN = nn => (im.Kp / (rho * nn * Dm * Dm / mu) + im.Np) * rho * nn ** 3 * Dm ** 5;
  checks.push(chk('Turbulent regime for the power-number and mixing-time correlations (Re > 10⁴)', Re > 1e4, `Re = ${F(Re, 4)}`));
  checks.push(chk('Solids fully suspended (N ≥ N_js)', Ns >= Njs, `N = ${F(Ns * 60, 4)} rpm, N_js = ${F(Njs * 60, 4)} rpm`));
  checks.push(chk('Impeller-to-tank diameter ratio in the usual range (0.25–0.6)', Dm / Tt >= 0.25 && Dm / Tt <= 0.6, `D/T = ${F(Dm / Tt, 3)}`));
  checks.push(chk('Particle size within the range of the Zwietering correlation (≳ 100 µm)', dp >= 1e-4, `d_p = ${F(dp * 1e6, 3)} µm`));
  checks.push(chk('Tip speed below 10 m/s (seal and erosion limit)', tip < 10, `tip speed ${F(tip, 3)} m/s`));
  if (Re <= 1e4) warnings.push(`Reynolds number ${F(Re, 3)} is below the fully turbulent range: the power number is estimated as K_p/Re + N_p,turb, which ignores the dip of the real power curve in the transitional region, and the transitional mixing-time correlation is used where applicable.`);
  if (dp < 1e-4) warnings.push('The Zwietering correlation was established for particles of roughly 0.1–2 mm; for fine agglomerates it is an extrapolation and tends to overestimate the just-suspended speed (fine particles follow the liquid).');
  if (Fr > 1) warnings.push('High Froude number: surface vortexing and air entrainment are likely unless the tank is fully baffled.');
  warnings.push('Standard baffled tank with liquid height equal to tank diameter and one impeller; local energy dissipation near the impeller is typically 10–50 times the mean used here.');
  const ns = N.logspace(Math.max(Ns / 10, 0.02), Ns * 5, 80);
  return {
    summary: `A ${F(Dm * 1e3, 3)} mm ${im.name} at ${F(Ns * 60, 4)} rpm in a ${F(Tt * 1e3, 3)} mm tank (${F(V * 1e3, 3)} L) runs at Re = ${F(Re, 3)} (${regime}) and draws ${F(P, 3)} W (${F(PV, 3)} W/m³). The 95 % mixing time is ${F(tMix, 3)} s, the Kolmogorov scale ${F(eta * 1e6, 3)} µm and the mean shear rate ${F(gMO, 3)} 1/s. The just-suspended speed for ${F(dp * 1e6, 3)} µm particles is ${F(Njs * 60, 4)} rpm, so solids are ${Ns >= Njs ? 'fully suspended' : 'not fully suspended'}.`,
    kpis: [
      kpi('Re', 'Impeller Reynolds number', Re, '', `ρND²/μ; ${regime}. Correlations with constant power number require more than about 10⁴.`, Re > 1e4),
      kpi('P', 'Shaft power', P, 'W', `P = N_p ρ N³ D⁵ with N_p = ${F(Np, 3)}. Power rises with the cube of speed and the fifth power of impeller diameter in turbulent flow.`),
      kpi('PV', 'Specific power input', PV, 'W/m³', 'Power per unit liquid volume; 100–1000 W/m³ is typical for solid suspension and precipitation, and it is a common scale-up criterion.'),
      kpi('tmix', 'Mixing time (95 %)', tMix, 's', crit > 6370 ? 'Grenville–Nienow turbulent correlation θ = 5.2 T²/(N_p^⅓ N D²) for H = T.' : 'Grenville transitional correlation N_p^⅓ Re = 183 Fo^(−½).'),
      kpi('eta', 'Kolmogorov length scale', eta * 1e6, 'µm', '(ν³/ε)^¼ from the mean dissipation. Agglomerates smaller than this experience viscous shear only; micromixing occurs on this scale.'),
      kpi('shear', 'Mean shear rate (Metzner–Otto)', gMO, '1/s', `γ̇ = k_s N with k_s = ${im.ks}; the turbulent (Kolmogorov) shear rate (ε/ν)^½ is ${F(gK, 3)} 1/s.`),
      kpi('tip', 'Impeller tip speed', tip, 'm/s', 'πND; sets maximum shear and is limited to a few m/s for fragile agglomerates.', tip < 10),
      kpi('Njs', 'Just-suspended speed (Zwietering)', Njs * 60, 'rpm', 'N_js = S ν^0.1 d_p^0.2 (gΔρ/ρ)^0.45 X^0.13 / D^0.85: lowest speed at which no particle rests on the bottom for more than about a second.'),
      kpi('ratio', 'N / N_js', Ns / Njs, '', 'Above one all solids are in motion; operating far above it wastes power.', Ns >= Njs),
      kpi('tcirc', 'Circulation time', tCirc, 's', `Tank volume over impeller pumping rate N_q N D³ (N_q = ${im.Nq}); mixing takes about four to five circulations.`),
      kpi('stress', 'Turbulent viscous stress on small agglomerates', stress, 'Pa', 'μ(ε/ν)^½ from the mean dissipation: agglomerates with a cohesive strength below this are broken up (higher near the impeller).'),
      kpi('Fr', 'Froude number', Fr, '', 'N²D/g; values above about one indicate surface vortexing in unbaffled tanks.'),
    ],
    series: [
      { key: 'pn', title: 'Power and just-suspended speed', xlabel: 'speed (rpm)', ylabel: 'power (W)', logx: true, logy: true, lines: [{ name: 'shaft power', x: ns.map(x => x * 60), y: ns.map(PofN) }, { name: 'just-suspended speed', x: [Njs * 60, Njs * 60], y: [PofN(ns[0]), PofN(ns[ns.length - 1])] }, { name: 'operating point', x: [Ns * 60], y: [P] }], explain: 'Power rises as N³ in turbulent flow; the vertical line marks the minimum speed for complete suspension.' },
      { key: 'tm', title: 'Mixing time and Kolmogorov scale vs speed', xlabel: 'speed (rpm)', ylabel: 'mixing time (s)  /  η (µm)', logx: true, logy: true, lines: [{ name: 'mixing time (s)', x: ns.map(x => x * 60), y: ns.map(nn => { const re = rho * nn * Dm * Dm / mu, np = im.Kp / re + im.Np, cr = Math.cbrt(np) * re; return cr > 6370 ? 5.2 * Tt * Tt / (Math.cbrt(np) * nn * Dm * Dm) : rho * Tt * Tt / mu * (183 / cr) ** 2; }) }, { name: 'Kolmogorov scale (µm)', x: ns.map(x => x * 60), y: ns.map(nn => Math.pow(nu ** 3 / (PofN(nn) / (rho * V)), 0.25) * 1e6) }], explain: 'Faster stirring shortens the mixing time roughly as 1/N and reduces the smallest eddy size as N^(−¾).' },
    ],
    equations: [...ZNO_EQ(),
      eqn('governing', 'Power draw', 'P = N_p ρ N³ D⁵,   Re = ρ N D² / μ'),
      eqn('constitutive', 'Mixing time (Grenville–Nienow)', 'θ₉₅ = 5.2 T^1.5 H^0.5 / (N_p^⅓ N D²)   (turbulent);   N_p^⅓ Re = 183 Fo^(−½)   (transitional)'),
      eqn('formula', 'Kolmogorov scales', 'ε = P / (ρ V);   η = (ν³/ε)^¼;   γ̇_K = (ε/ν)^½'),
      eqn('constitutive', 'Just-suspended speed (Zwietering)', 'N_js = S ν^0.1 d_p^0.2 (g Δρ/ρ)^0.45 X^0.13 / D^0.85'),
      eqn('constitutive', 'Metzner–Otto shear rate', 'γ̇ = k_s N')],
    checks, warnings,
    recommendations: [
      rec('high', Ns < Njs ? `Solids are not fully suspended: raise the speed to at least ${F(Njs * 60, 4)} rpm, lower the impeller clearance or use a larger down-pumping impeller.` : Ns > 1.5 * Njs ? `Speed is ${F(Ns / Njs, 3)} times the just-suspended speed. Reducing to about ${F(1.1 * Njs * 60, 4)} rpm would cut power to roughly ${F(100 * Math.pow(1.1 * Njs / Ns, 3), 3)} % of the present value if mixing time (then about ${F(tMix * Ns / (1.1 * Njs), 3)} s) remains acceptable.` : 'Speed is just above the suspension limit: an efficient operating point.', 'Zwietering just-suspended speed.'),
      rec('medium', `Mixing time is ${F(tMix, 3)} s. For precipitation, the reagent feed time should be at least ten times longer (${F(10 * tMix, 3)} s), or feed close to the impeller, so that local supersaturation peaks do not broaden the size distribution.`, 'Mixing time against reaction and nucleation time scales.'),
      rec('medium', dp > 5 * eta ? `Agglomerates (${F(dp * 1e6, 3)} µm) are larger than the Kolmogorov scale (${F(eta * 1e6, 3)} µm): they are exposed to inertial turbulent stresses and will be broken down towards that scale.` : `Agglomerates (${F(dp * 1e6, 3)} µm) are below or near the Kolmogorov scale (${F(eta * 1e6, 3)} µm): stirring will not disperse them further; use a rotor–stator or ultrasound if smaller sizes are required.`, 'Agglomerate size relative to the Kolmogorov scale.'),
      rec('low', `Specific power is ${F(PV, 3)} W/m³. An axial hydrofoil suspends solids at a fraction of the power of a Rushton turbine (turbulent power numbers 0.3 against 5); choosing it and running near N_js minimises energy per batch.`, 'Power numbers of the impeller types.', true),
    ],
    publish: { 'proc.power': P, 'proc.volume': V, 'mix.time_s': tMix, 'mix.shear_rate': gMO, 'mix.kolmogorov_m': eta, 'mix.Njs_rpm': Njs * 60, 'flow.Re': Re },
    numerics: { method: 'Closed-form stirred-tank correlations (power number, Grenville–Nienow mixing time, Zwietering, Metzner–Otto)', dof: 1, steps: ns.length, iterations: 1, residual: 0, notes: `Kolmogorov time ${F(tauK, 3)} s; laminar power constant K_p ≈ ${im.Kp} (approximate).` },
  };
}

// ---------------------------------------------------------------- particle–fluid
function runParticle(inp) {
  const rho = num(inp, 'rho'), mu = num(inp, 'mu'), d = num(inp, 'dpS'), rs = num(inp, 'rhoS'), T = num(inp, 'Tp') + 273.15, phi = num(inp, 'phiS'), Uc = num(inp, 'Uc'), Lc = num(inp, 'Lc'), gT = num(inp, 'gradT'), DT = num(inp, 'DT'), lam = num(inp, 'lam'), warnings = [], checks = [];
  const ccOf = x => { if (lam <= 0) return 1; const Kn = 2 * lam / x; return 1 + Kn * (1.257 + 0.4 * Math.exp(-1.1 / Kn)); };
  const vtOf = x => { const dr = Math.abs(rs - rho), W = (PI / 6) * x ** 3 * dr * G, cc = ccOf(x); if (W <= 0) return 0; const drag = v => { const re = rho * v * x / mu, cd = re < 1000 ? 24 / Math.max(re, 1e-300) * (1 + 0.15 * Math.pow(re, 0.687)) / cc : 0.44; return cd * (PI / 4) * x * x * rho * v * v / 2 - W; }; return safe(() => Math.exp(N.brent(lv => drag(Math.exp(lv)), Math.log(1e-30), Math.log(1e4), 1e-12)), cc * dr * G * x * x / (18 * mu)); };
  const Cc = ccOf(d), sgn = rs >= rho ? 1 : -1, vt = vtOf(d), vSt = Cc * Math.abs(rs - rho) * G * d * d / (18 * mu), Ret = rho * vt * d / mu;
  const nRZ = Ret < 0.2 ? 4.65 : Ret < 1 ? 4.4 * Math.pow(Ret, -0.03) : Ret < 500 ? 4.4 * Math.pow(Ret, -0.1) : 2.39, vh = vt * Math.pow(1 - phi, nRZ);
  const tauP = rs * d * d * Cc / (18 * mu), Stk = tauP * Uc / Lc, Dof = x => kB * T * ccOf(x) / (3 * PI * mu * x), D = Dof(d), vT = DT * gT, ST = DT / D, rms = Math.sqrt(2 * D), PeS = vt * (d / 2) / D;
  const dom = [['settling', vt], ['thermophoresis', vT], ['Brownian motion', rms]].sort((a, b) => b[1] - a[1])[0][0];
  checks.push(chk('Drag law converged to the force balance', Math.abs(vt - vSt) / Math.max(vSt, 1e-300) < 0.05 || Ret > 0.05, `terminal velocity ${F(vt, 3)} m/s, Stokes value ${F(vSt, 3)} m/s, Re_t = ${F(Ret, 3)}`));
  checks.push(chk('Hindered-settling correlation range (volume fraction ≤ 0.4)', phi <= 0.4, `φ = ${F(phi, 3)}`));
  checks.push(chk('Particle follows the flow (Stokes number < 0.1)', Stk < 0.1, `Stk = ${F(Stk, 3)}`));
  checks.push(chk('Hindered velocity does not exceed the single-particle value', vh <= vt * (1 + 1e-12), `factor (1−φ)^n = ${F(Math.pow(1 - phi, nRZ), 3)}`));
  if (phi > 0.4) warnings.push('Above about 40 vol% the Richardson–Zaki correlation is unreliable; the suspension approaches a packed, network-forming state.');
  if (Stk >= 1) warnings.push(`Stokes number ${F(Stk, 3)}: particles leave the streamlines at bends and obstacles and deposit by inertial impaction.`);
  if (lam > 0 && 2 * lam / d > 10) warnings.push('Free-molecular regime (Kn > 10): the Cunningham-corrected Stokes drag is an interpolation formula; thermophoresis follows the Waldmann limit.');
  warnings.push('Isolated rigid spheres; agglomerates have a lower effective density and a larger drag than the equivalent solid sphere, and electrostatic or DLVO forces near walls are not included.');
  const ds = N.logspace(1e-9, 1e-3, 90), ph = N.linspace(0, 0.5, 60);
  return {
    summary: `A ${F(d * 1e6, 3)} µm particle of density ${F(rs, 4)} kg/m³ ${sgn > 0 ? 'settles' : 'rises'} at ${F(vt, 3)} m/s (Re = ${F(Ret, 3)}); at ${F(100 * phi, 3)} vol% the hindered velocity is ${F(vh, 3)} m/s (Richardson–Zaki exponent ${F(nRZ, 3)}). Its Stokes number in the stated flow is ${F(Stk, 3)}. Over one second it moves ${F(vt * 1e6, 3)} µm by settling, ${F(vT * 1e6, 3)} µm by thermophoresis and ${F(rms * 1e6, 3)} µm (rms) by Brownian motion: ${dom} dominates.`,
    kpis: [
      kpi('vt', 'Terminal velocity (single particle)', sgn * vt, 'm/s', `Force balance with the Schiller–Naumann drag law${lam > 0 ? ' and Cunningham slip correction ' + F(Cc, 4) : ''}; negative means the particle rises. The Stokes value is ${F(vSt, 3)} m/s.`),
      kpi('Ret', 'Particle Reynolds number', Ret, '', 'ρ v d/μ at terminal velocity; Stokes drag is accurate below about 0.1.'),
      kpi('vh', 'Hindered settling velocity', sgn * vh, 'm/s', `Richardson–Zaki v = v_t(1 − φ)^n with n = ${F(nRZ, 3)}: return flow of displaced liquid and crowding slow the settling of concentrated suspensions.`),
      kpi('tauP', 'Particle relaxation time', tauP, 's', 'ρ_p d² C_c/(18μ): time for a particle to adjust to a change of fluid velocity.'),
      kpi('Stk', 'Stokes number', Stk, '', 'τ_p U/L. Well below one the particle is a tracer; of order one or more it impacts on obstacles and walls.', Stk < 0.1),
      kpi('D', 'Brownian diffusivity', D, 'm²/s', 'Stokes–Einstein (with slip correction in gases).'),
      kpi('vT', 'Thermophoretic velocity', vT, 'm/s', 'D_T·|∇T|, directed from hot to cold for most systems; it keeps particles off hot walls and drives them to cold ones.'),
      kpi('ST', 'Soret coefficient', ST, '1/K', 'D_T/D; multiplied by a temperature difference it gives the logarithmic concentration ratio at steady state.'),
      kpi('PeS', 'Sedimentation Péclet number', PeS, '', 'v_t·a/D; below one Brownian motion keeps the particle suspended.'),
      kpi('dom', 'Dominant transport mechanism over one second', dom, '', 'Largest of the settling distance, thermophoretic drift and Brownian root-mean-square displacement in one second.'),
    ],
    series: [
      { key: 'mech', title: 'Displacement in one second vs particle diameter', xlabel: 'diameter (µm)', ylabel: 'displacement (m)', logx: true, logy: true, lines: [{ name: 'settling', x: ds.map(x => x * 1e6), y: ds.map(x => Math.max(vtOf(x), 1e-30)) }, { name: 'Brownian (rms)', x: ds.map(x => x * 1e6), y: ds.map(x => Math.sqrt(2 * Dof(x))) }, { name: 'thermophoresis', x: [ds[0] * 1e6, ds[ds.length - 1] * 1e6], y: [Math.max(vT, 1e-30), Math.max(vT, 1e-30)] }], explain: 'Brownian motion dominates for the smallest particles, gravity for the largest; thermophoresis (taken as size independent) can dominate in between when gradients are steep.' },
      { key: 'rz', title: 'Hindered settling factor', xlabel: 'solids volume fraction', ylabel: 'v / v_t', lines: [{ name: `Richardson–Zaki, n = ${F(nRZ, 3)}`, x: ph, y: ph.map(x => Math.pow(1 - x, nRZ)) }, { name: 'solids flux φ·v/v_t', x: ph, y: ph.map(x => x * Math.pow(1 - x, nRZ)) }], explain: 'Settling slows steeply with concentration. The solids-flux curve has a maximum that sets the capacity of thickeners and the limit of batch settling.' },
    ],
    equations: [
      eqn('governing', 'Force balance at terminal velocity', '(π/6) d³ (ρ_p − ρ) g = C_D (π/4) d² ρ v²/2,   C_D = (24/Re)(1 + 0.15 Re^0.687)   (Re < 1000)'),
      eqn('constitutive', 'Richardson–Zaki hindered settling', 'v_h = v_t (1 − φ)ⁿ,   n = 4.65 (Re < 0.2) … 2.39 (Re > 500)'),
      eqn('formula', 'Stokes number and relaxation time', 'Stk = τ_p U / L,   τ_p = ρ_p d² C_c / (18 μ)'),
      eqn('governing', 'Brownian and thermophoretic transport', 'D = k T C_c / (3 π μ d);   v_T = −D_T ∇T;   S_T = D_T / D'),
      eqn('constitutive', 'Cunningham slip correction', 'C_c = 1 + Kn [1.257 + 0.4 exp(−1.1/Kn)],   Kn = 2λ/d')],
    checks, warnings,
    recommendations: [
      rec('high', dom === 'settling' ? `Settling dominates: in ${F(Lc / Math.max(vh, 1e-300), 3)} s a particle crosses the characteristic length of ${F(Lc * 1e3, 3)} mm. Keep the suspension agitated or the flow velocity above the settling velocity by a wide margin to avoid deposits in horizontal sections.` : dom === 'thermophoresis' ? `Thermophoresis dominates (${F(vT * 1e6, 3)} µm/s): use it deliberately — a cold collection surface gathers the particles, and hot reactor walls stay clean.` : 'Brownian motion dominates: particles reach walls by diffusion, so deposition scales with residence time and surface-to-volume ratio; coat or cool the walls, or shorten the residence time, to limit losses.', 'Comparison of displacements per second.'),
      rec('medium', Stk > 0.1 ? `Stokes number ${F(Stk, 3)}: expect inertial deposition in bends and at sampling probes; use gentle bends (larger radius) or isokinetic sampling.` : `Stokes number ${F(Stk, 3)}: particles follow the streamlines, so sampling is representative and inertial separators (cyclones, impactors) will not collect them — use filtration, electrostatic or thermophoretic collection.`, 'Stokes number.'),
      rec('medium', phi > 0.1 ? `At ${F(100 * phi, 3)} vol% the settling rate is only ${F(100 * Math.pow(1 - phi, nRZ), 3)} % of the dilute value: dilute before gravity separation, or accept longer settling times.` : 'The suspension is dilute enough for single-particle settling rates to apply.', 'Richardson–Zaki factor.'),
      rec('low', 'Airborne or suspended nanoparticles with Stokes numbers far below one are not removed by inertia: enclose the process and use HEPA filtration or wet scrubbing for exhaust, and capture product by thermophoresis or filtration rather than venting.', 'Low inertia of nanoparticles.', true),
    ],
    publish: { 'flow.settling_velocity': sgn * vh, 'prop.D_diff': D, 'flow.stokes': Stk, 'flow.thermophoretic_velocity': vT },
    numerics: { method: 'Brent root finding on the drag–weight balance (Schiller–Naumann drag); closed-form Richardson–Zaki, Stokes–Einstein and slip relations', dof: 1, steps: ds.length, iterations: ds.length + 1, residual: 0, notes: 'Thermophoretic mobility is a user input and is taken as size independent.' },
  };
}

const RUN = { bed: runBed, darcy2d: runDarcy2D, channel: runChannel, stirred: runStirred, particle: runParticle };

export default {
  id: 9,
  tagline: 'Flow and transport in porous reactors, heterogeneous Darcy fields, channels and coils, stirred vessels and particle-laden fluids.',
  inputs,
  presets: [
    { name: 'ZnO hydrothermal reference: porous slab, porosity 0.25', note: 'Reference case: 1 m × 1 m × 0.01 m porous reactor, D(Zn²⁺) = 7.5×10⁻¹⁰ m²/s, inlet 10 mol/m³ (0.01 mol/L), heterogeneous permeability.', values: { mode: 'darcy2d', por: 0.25, Lx: 1, Ly: 1, thick: 0.01, Kg: 1e-12, sigK: 1, corr: 3, dP: 1e5, pout: 101325, cin: 10, Dmol: 7.5e-10, k1: 1e-3, rho: 997, mu: 8.9e-4 } },
    { name: 'Packed bed of 100 µm grains, porosity 0.25', note: 'One-dimensional hydraulics of the reference porous reactor.', values: { mode: 'bed', por: 0.25, dg: 1e-4, Lbed: 0.1, Abed: 1e-3, Q: 1e-7, Dmol: 7.5e-10, k1: 1e-3 } },
    { name: 'Coiled-tube microreactor for nanoparticle synthesis', note: '1 mm PTFE coil, 60 mL/h: laminar flow with Taylor–Aris dispersion.', values: { mode: 'channel', geo: 'coil', dh: 1e-3, Lch: 10, Q: 1.67e-8, Dc: 0.05, Dmol: 1e-9, k1: 0.01 } },
    { name: 'Stirred hydrothermal precipitation vessel, 5 L', note: 'Pitched-blade turbine suspending ZnO agglomerates.', values: { mode: 'stirred', imp: 'pbt', Dimp: 0.06, Tt: 0.18, Nrpm: 700, dpS: 2e-5, rhoS: 5606, Xs: 1, Szw: 6 } },
    { name: 'Nanoparticles in a CVD gas stream (thermophoresis)', note: '50 nm oxide particles in hot gas near a cooled substrate; metamaterial film deposition.', values: { mode: 'particle', rho: 0.6, mu: 3e-5, Tp: 300, dpS: 5e-8, rhoS: 5606, phiS: 0, Uc: 1, Lc: 0.01, gradT: 1e5, DT: 5e-8, lam: 1.3e-7 } },
  ],
  run(inp) { return RUN[sel(inp, 'mode')](inp); },
};
