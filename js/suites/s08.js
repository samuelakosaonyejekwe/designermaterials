// Suite 8 — Surface, interface and colloidal engineering.
// DLVO stability, Gouy–Chapman–Stern double layer (numerical Poisson–Boltzmann), adsorption-isotherm
// fitting (Langmuir, Freundlich, BET), wetting and capillarity of nanorod arrays, particle transport.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS, PI = Math.PI, kB = N.K_B, qe = N.Q_E, NA = N.N_A, E0 = N.EPS0;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
const deg = x => x * 180 / PI, rad = x => x * PI / 180;

// synthetic BET demonstration isotherm (N2, 77 K): forward BET model + seeded noise
const DEMO_ISO = (() => { const r = N.rng(21), rows = ['p_rel,q_mmol_g'], qm = 0.5125, C = 100; for (const x of [0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]) rows.push(`${x},${(qm * C * x / ((1 - x) * (1 + (C - 1) * x)) * (1 + 0.004 * r.normal())).toFixed(4)}`); return rows.join('\n'); })();

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'dlvo', group: 'Model', options: [
    { value: 'dlvo', label: 'DLVO colloidal stability' }, { value: 'edl', label: 'Electric double layer and zeta potential' }, { value: 'adsorb', label: 'Adsorption isotherm fitting (Langmuir, Freundlich, BET)' },
    { value: 'wetting', label: 'Wetting and capillarity of nanostructured surfaces' }, { value: 'transport', label: 'Particle diffusion and sedimentation' }], help: 'Selects the surface or colloid model to run.' },
  { key: 'T', label: 'Temperature', unit: '°C', value: 25, min: 0, max: 150, group: 'Medium', bus: 'proc.T_C', showIf: { mode: ['dlvo', 'edl', 'wetting', 'transport'] }, help: 'Temperature of the liquid. It sets the thermal energy kT, the Debye length and the Brownian diffusivity.' },
  { key: 'epsr', label: 'Relative permittivity of the liquid', unit: '', value: 78.5, min: 2, max: 120, group: 'Medium', showIf: { mode: ['dlvo', 'edl'] }, help: '78.5 for water at 25 °C, 24.5 for ethanol. Lower permittivity shortens the Debye length and weakens double-layer repulsion.' },
  { key: 'mu', label: 'Liquid viscosity', unit: 'Pa·s', value: 8.9e-4, min: 1e-5, max: 10, log: true, group: 'Medium', showIf: { mode: ['dlvo', 'edl', 'transport'] }, help: 'Dynamic viscosity (8.9×10⁻⁴ Pa·s for water at 25 °C). Controls diffusion, settling and aggregation rates.' },
  { key: 'dp', label: 'Particle diameter', unit: 'm', value: 1e-7, min: 1e-9, max: 1e-4, log: true, group: 'Particle', bus: 'geom.diameter', showIf: { mode: ['dlvo', 'edl', 'adsorb', 'transport'] }, help: 'Diameter of the (spherical) particle; 1×10⁻⁷ m = 100 nm.' },
  { key: 'AH', label: 'Hamaker constant (particle–liquid–particle)', unit: 'zJ', value: 20, min: 0.1, max: 500, log: true, group: 'Particle', showIf: { mode: ['dlvo'] }, help: 'Effective Hamaker constant across the liquid in 10⁻²¹ J. Typical in water: silica ≈ 5–8, ZnO ≈ 20, TiO₂ ≈ 50, gold ≈ 250–300.' },
  { key: 'psi0', label: 'Surface (zeta) potential', unit: 'mV', value: -30, min: -300, max: 300, group: 'Particle', bus: 'colloid.zeta_mV', showIf: { mode: ['dlvo', 'edl'] }, help: 'Potential at the particle surface; in practice the measured zeta potential is used. Its sign is irrelevant for identical particles.' },
  { key: 'cel', label: 'Electrolyte concentration', unit: 'mol/L', value: 0.01, min: 1e-6, max: 5, log: true, group: 'Medium', showIf: { mode: ['dlvo', 'edl'] }, help: 'Concentration of a symmetric z:z electrolyte. It screens the double layer: the Debye length is 0.304/(z√c) nm in water at 25 °C.' },
  { key: 'zval', label: 'Ion valence z', unit: '', value: 1, min: 1, max: 3, step: 1, group: 'Medium', showIf: { mode: ['dlvo', 'edl'] }, help: 'Valence of the symmetric electrolyte. Multivalent counter-ions destabilise colloids at far lower concentration (Schulze–Hardy rule).' },
  { key: 'phi', label: 'Particle volume fraction', unit: '', value: 1e-4, min: 1e-8, max: 0.3, log: true, group: 'Particle', showIf: { mode: ['dlvo'] }, help: 'Solids volume fraction; sets the number density and so the aggregation half-time.' },
  { key: 'dStern', label: 'Stern layer thickness', unit: 'nm', value: 0.3, min: 0, max: 2, group: 'Double layer', showIf: { mode: ['edl'] }, help: 'Thickness of the ion-free compact layer (about one hydrated ion radius). Zero gives the pure Gouy–Chapman model.' },
  { key: 'epsS', label: 'Relative permittivity of the Stern layer', unit: '', value: 10, min: 2, max: 80, group: 'Double layer', showIf: { mode: ['edl'] }, help: 'Oriented interfacial water has a much lower permittivity (about 6–30) than the bulk.' },
  { key: 'iso', label: 'Isotherm data', type: 'table', columns: ['p_rel', 'q_mmol_g'], value: DEMO_ISO, group: 'Data', showIf: { mode: ['adsorb'] }, help: 'Adsorbed amount (mmol/g) against relative pressure p/p₀ (or a normalised concentration). Paste or upload CSV with a header line. The default is synthetic demonstration data.' },
  { key: 'sigN', label: 'Adsorbate cross-sectional area', unit: 'nm²', value: 0.162, min: 0.05, max: 2, group: 'Data', showIf: { mode: ['adsorb'] }, help: 'Area occupied by one adsorbed molecule in the monolayer: 0.162 nm² for N₂ at 77 K, 0.142 nm² for Ar at 87 K.' },
  { key: 'xq', label: 'Relative pressure / concentration for coverage', unit: '', value: 0.2, min: 0.001, max: 0.95, group: 'Data', showIf: { mode: ['adsorb'] }, help: 'Point on the isotherm at which coverage and surface density are reported.' },
  { key: 'rhoP', label: 'Particle density', unit: 'kg/m³', value: 5606, min: 500, max: 25000, group: 'Particle', bus: 'prop.density', showIf: { mode: ['adsorb', 'transport'] }, help: 'Skeletal density of the solid (5606 kg/m³ for ZnO). Used for the BET-equivalent diameter, the geometric surface area and sedimentation.' },
  { key: 'thetaY', label: 'Intrinsic (Young) contact angle', unit: '°', value: 50, min: 0, max: 180, group: 'Wetting', showIf: { mode: ['wetting'] }, help: 'Contact angle on a flat surface of the same chemistry. As-grown oxides are hydrophilic (below 90°); fluorosilane or alkylthiol coatings give 105–120°.' },
  { key: 'gLV', label: 'Liquid surface tension', unit: 'N/m', value: 0.072, min: 0.005, max: 0.6, group: 'Wetting', showIf: { mode: ['wetting'] }, help: '0.072 N/m for water, 0.022 for ethanol or isopropanol, about 0.012 for fluorinated rinse liquids.' },
  { key: 'dRod', label: 'Nanorod diameter', unit: 'm', value: 8e-8, min: 5e-9, max: 1e-5, log: true, group: 'Nanostructure', bus: 'geom.diameter', showIf: { mode: ['wetting'] }, help: 'Diameter of the rods forming the surface texture (8×10⁻⁸ m = 80 nm).' },
  { key: 'LRod', label: 'Nanorod length', unit: 'm', value: 1.5e-6, min: 1e-8, max: 1e-3, log: true, group: 'Nanostructure', bus: 'geom.length', showIf: { mode: ['wetting'] }, help: 'Height of the rods. Roughness and collapse risk both rise with length.' },
  { key: 'Nrod', label: 'Nanorod number density', unit: '1/m²', value: 4e13, min: 1e9, max: 1e16, log: true, group: 'Nanostructure', bus: 'geom.number_density', showIf: { mode: ['wetting'] }, help: 'Rods per unit substrate area; the mean pitch is its inverse square root.' },
  { key: 'Emod', label: 'Young\'s modulus of the rods', unit: 'Pa', value: 1.4e11, min: 1e6, max: 1.2e12, log: true, group: 'Nanostructure', bus: 'prop.E', showIf: { mode: ['wetting'] }, help: 'Bending stiffness of a rod scales with E·d⁴ (1.4×10¹¹ Pa for ZnO; polymers are 100 times softer).' },
  { key: 'RH', label: 'Relative humidity (or p/p₀ of the vapour)', unit: '%', value: 60, min: 1, max: 99.9, group: 'Wetting', showIf: { mode: ['wetting'] }, help: 'Used with the Kelvin equation to judge capillary condensation in the gaps between rods.' },
  { key: 'Vml', label: 'Liquid molar volume', unit: 'm³/mol', value: 1.8e-5, min: 1e-5, max: 5e-4, log: true, group: 'Wetting', sens: false, showIf: { mode: ['wetting'] }, help: 'Molar volume of the condensing liquid (1.8×10⁻⁵ m³/mol for water).' },
  { key: 'rhoF', label: 'Liquid density', unit: 'kg/m³', value: 997, min: 0.1, max: 3000, group: 'Medium', showIf: { mode: ['transport'] }, help: 'Density of the suspending fluid.' },
  { key: 'Hcol', label: 'Settling height', unit: 'm', value: 0.05, min: 1e-4, max: 10, log: true, group: 'Particle', showIf: { mode: ['transport'] }, help: 'Height of the liquid column or cuvette over which settling is judged.' },
  { key: 'Dmeas', label: 'Measured diffusion coefficient', unit: 'm²/s', value: 4.9e-12, min: 1e-14, max: 1e-8, log: true, group: 'Particle', bus: 'prop.D_diff', showIf: { mode: ['transport'] }, help: 'Translational diffusion coefficient from dynamic light scattering or particle tracking; converted to a hydrodynamic diameter.' },
  { key: 'RCF', label: 'Relative centrifugal force', unit: '× g', value: 1, min: 1, max: 5e5, log: true, group: 'Particle', showIf: { mode: ['transport'] }, help: '1 for gravity settling; e.g. 10 000 for a benchtop centrifuge.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
const num = (inp, k) => { const d = DEF[k]; let v = +inp[k]; if (!Number.isFinite(v)) v = d.value; return N.clamp(v, d.min, d.max); };
const sel = (inp, k) => (DEF[k].options.some(o => o.value === inp[k]) ? inp[k] : DEF[k].value);
const kappaOf = (cL, z, eps, T) => Math.sqrt(2 * z * z * qe * qe * NA * cL * 1000 / (eps * E0 * kB * T));

// ---------------------------------------------------------------- DLVO
function runDLVO(inp) {
  const T = num(inp, 'T') + 273.15, eps = num(inp, 'epsr'), mu = num(inp, 'mu'), a = num(inp, 'dp') / 2, A = num(inp, 'AH') * 1e-21, psi = num(inp, 'psi0') * 1e-3, c = num(inp, 'cel'), z = Math.round(num(inp, 'zval')), phi = num(inp, 'phi');
  const kT = kB * T, warnings = [], checks = [], h0 = 0.158e-9;
  const vdw = h => -A / 6 * (2 * a * a / (h * (4 * a + h)) + 2 * a * a / (2 * a + h) ** 2 + Math.log(h * (4 * a + h) / (2 * a + h) ** 2));
  const edl = (h, k) => 2 * PI * eps * E0 * a * psi * psi * Math.log(1 + Math.exp(-k * h));
  const profile = (cc, zz, n = 400) => { const k = kappaOf(cc, zz, eps, T), hs = N.logspace(h0, Math.max(30 / k, 50 * h0, 2e-8), n), V = hs.map(h => (vdw(h) + edl(h, k)) / kT); return { k, hs, V }; };
  const analyse = p => { const iM = N.argmax(p.V), Vmax = p.V[iM]; let iS = iM, Vs = Vmax; for (let i = iM; i < p.V.length; i++) if (p.V[i] < Vs) { Vs = p.V[i]; iS = i; } return { iM, Vmax, iS, Vs }; };
  const P = profile(c, z), an = analyse(P), kap = P.k, lamD = 1 / kap, barrier = Math.max(an.Vmax, 0), hasBar = an.Vmax > 0 && an.iM < P.V.length - 1;
  const sec = an.Vs < -1e-3 && an.iS > an.iM ? an.Vs : 0;
  // Fuchs stability ratio
  const lh = P.hs.map(Math.log), integ = P.hs.map((h, i) => 2 * a * Math.exp(Math.min(P.V[i], 690)) / (h + 2 * a) ** 2 * h);
  const hl = P.hs[P.hs.length - 1], W = Math.max(N.trapz(lh, integ) + 2 * a / (hl + 2 * a), 1e-3), logW = Math.log10(W);
  // critical coagulation concentration: barrier vanishes
  const barAt = (cc, zz) => analyse(profile(cc, zz, 160)).Vmax;
  const ccc = zz => { if (barAt(1e-6, zz) <= 1e-6) return 0; if (barAt(5, zz) > 1e-6) return Infinity; let lo = -6, hi = Math.log10(5); for (let i = 0; i < 36; i++) { const m = 0.5 * (lo + hi); if (barAt(Math.pow(10, m), zz) > 1e-6) lo = m; else hi = m; } return Math.pow(10, 0.5 * (lo + hi)); };
  const cc = [1, 2, 3].map(ccc), cz = cc[z - 1], okC = cc.every(v => v > 0 && Number.isFinite(v));
  const shExp = okC ? N.linreg([0, Math.log(2), Math.log(3)], cc.map(Math.log)).slope : null;
  const N0 = phi / ((4 / 3) * PI * a ** 3), tRapid = 3 * mu / (4 * kT * N0), tHalf = tRapid * W;
  const verdict = barrier > 20 ? 'kinetically stable' : barrier > 5 ? 'marginally stable (slow aggregation)' : 'unstable (rapid aggregation)';
  const ka = kap * a, cstr = v => (v === 0 ? 'unstable at any concentration' : !Number.isFinite(v) ? 'above 5 mol/L' : v);
  checks.push(chk('Derjaguin approximation valid (κa > 5)', ka > 5, `κa = ${F(ka, 3)}`));
  checks.push(chk('Weak-overlap/low-potential form of the double-layer term (|zeψ₀/kT| ≲ 2)', Math.abs(z * qe * psi / kT) <= 2, `zeψ₀/kT = ${F(z * qe * psi / kT, 3)}`));
  checks.push(chk('Interaction energy vanishes at large separation', Math.abs(P.V[P.V.length - 1]) < 0.05 * Math.max(1, barrier), `V(${F(hl * 1e9, 3)} nm) = ${F(P.V[P.V.length - 1], 2)} kT`));
  checks.push(chk('Stability ratio at least of order one', W >= 0.5, `W = ${logW < 6 ? F(W, 3) : '10^' + F(logW, 4)}`));
  if (ka <= 5) warnings.push('κa ≤ 5: the Derjaguin approximation overestimates the double-layer repulsion for small particles or thick double layers; use the linear-superposition or full Poisson–Boltzmann sphere–sphere solution.');
  if (Math.abs(z * qe * psi / kT) > 2) warnings.push('The surface potential exceeds about 2kT/ze: the Hogg–Healy–Fuerstenau constant-potential expression is a linearised result and overestimates repulsion at high potential.');
  if (c > 0.5) warnings.push('Above about 0.5 mol/L the Debye length is comparable to ion size and non-DLVO forces (hydration, ion correlations) dominate.');
  warnings.push('Classical DLVO only: steric, hydration, depletion and hydrophobic forces, retardation of the van der Waals force and hydrodynamic resistance in the Fuchs integral are not included.');
  const hn = P.hs.map(h => h * 1e9), clip = v => N.clamp(v, -50, Math.max(60, 1.2 * barrier));
  const cs = N.logspace(1e-4, 2, 60), bc = cs.map(x => Math.max(analyse(profile(x, z, 120)).Vmax, 0));
  return {
    summary: `For ${F(2 * a * 1e9, 3)} nm particles at ${F(psi * 1e3, 3)} mV in ${F(c, 3)} mol/L ${z}:${z} electrolyte the Debye length is ${F(lamD * 1e9, 3)} nm and the DLVO energy barrier is ${F(barrier, 3)} kT: the dispersion is ${verdict}. The Fuchs stability ratio is ${logW < 6 ? F(W, 3) : '10^' + F(logW, 3)}, giving an aggregation half-time of ${F(Math.min(tHalf, 1e300), 3)} s at a volume fraction of ${F(phi, 2)}. The critical coagulation concentration is ${typeof cstr(cz) === 'number' ? F(cz, 3) + ' mol/L' : cstr(cz)}.`,
    kpis: [
      kpi('debye', 'Debye length', lamD * 1e9, 'nm', 'κ⁻¹ = [εkT/(2z²e²N_A c)]^½, the thickness of the diffuse double layer; repulsion acts over a few Debye lengths.'),
      kpi('barrier', 'Energy barrier', barrier, 'kT', 'Maximum of the total DLVO interaction energy. Above about 15–20 kT collisions rarely cross it and the sol is kinetically stable; below about 5 kT it aggregates quickly.', barrier > 15 ? true : barrier < 5 ? false : null),
      kpi('h_bar', 'Separation at the barrier', hasBar ? P.hs[an.iM] * 1e9 : 0, 'nm', 'Surface-to-surface distance of the energy maximum; zero means there is no interior maximum.'),
      kpi('sec_min', 'Secondary minimum depth', sec, 'kT', 'Shallow minimum beyond the barrier. Depths beyond a few kT cause weak, reversible flocculation that redisperses on shaking or dilution.', sec > -2),
      kpi('logW', 'log₁₀ stability ratio W', logW, '', 'Fuchs integral W = 2a∫exp(V/kT)/(h+2a)²dh: the factor by which aggregation is slower than diffusion-limited. Published to the kinetics suite.'),
      kpi('ccc', 'Critical coagulation concentration', cstr(cz), 'mol/L', `Electrolyte concentration at which the barrier vanishes for valence ${z} (bisection on the model barrier). Operate well below it.`),
      kpi('sh', 'Schulze–Hardy exponent', shExp === null ? 'not defined' : shExp, '', `Slope of log(ccc) against log z from the model for z = 1, 2, 3 (ccc = ${cc.map(v => (typeof cstr(v) === 'number' ? F(v, 3) : cstr(v))).join(', ')} mol/L). Theory gives −6 at high and −2 at low surface potential.`),
      kpi('t_half', 'Aggregation half-time', Math.min(tHalf, 1e300), 's', `Time for the particle number to halve, 3μW/(4kT N₀), at ${F(N0, 3)} particles per m³; the rapid (W = 1) value would be ${F(tRapid, 3)} s.`, tHalf > 86400 ? true : tHalf < 60 ? false : null),
      kpi('ka', 'κa', ka, '', 'Ratio of particle radius to Debye length; decides which double-layer and electrophoresis approximations apply.'),
      kpi('verdict', 'Stability assessment', verdict, '', 'Qualitative classification from the barrier height alone; steric stabilisers are not accounted for.'),
    ],
    series: [
      { key: 'V', title: 'DLVO interaction energy vs separation', xlabel: 'surface separation (nm)', ylabel: 'energy (kT)', logx: true, lines: [{ name: 'total', x: hn, y: P.V.map(clip) }, { name: 'van der Waals', x: hn, y: P.hs.map(h => clip(vdw(h) / kT)) }, { name: 'double layer', x: hn, y: P.hs.map(h => clip(edl(h, kap) / kT)) }], explain: 'Sum of van der Waals attraction and double-layer repulsion between two equal spheres (clipped for display). The maximum is the barrier against aggregation into the deep primary minimum at contact.' },
      { key: 'bc', title: 'Energy barrier vs electrolyte concentration', xlabel: 'electrolyte concentration (mol/L)', ylabel: 'barrier (kT)', logx: true, lines: [{ name: `z = ${z}`, x: cs, y: bc }, { name: '15 kT guide', x: [cs[0], cs[cs.length - 1]], y: [15, 15] }], explain: 'Salt screens the repulsion and lowers the barrier; where the curve reaches zero is the critical coagulation concentration.' },
    ],
    equations: [
      eqn('governing', 'DLVO total interaction', 'V(h) = V_vdW(h) + V_EDL(h)'),
      eqn('constitutive', 'van der Waals (Hamaker, equal spheres)', 'V_vdW = −(A/6) { 2a²/[h(4a + h)] + 2a²/(2a + h)² + ln[ h(4a + h)/(2a + h)² ] }'),
      eqn('constitutive', 'Double layer (constant potential, Derjaguin)', 'V_EDL = 2π ε ε₀ a ψ₀² ln(1 + e^(−κh))'),
      eqn('formula', 'Debye parameter', 'κ² = 2 z² e² N_A c / (ε ε₀ k T)'),
      eqn('governing', 'Fuchs stability ratio and half-time', 'W = 2a ∫₀^∞ exp(V/kT) / (h + 2a)² dh;   t_½ = 3 μ W / (4 k T N₀)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', barrier < 15 ? `The barrier of ${F(barrier, 3)} kT is too low for long-term stability. Raise the magnitude of the zeta potential (pH further from the isoelectric point, charged dispersant), lower the ionic strength below ${typeof cstr(cz) === 'number' ? F(cz / 10, 2) : '1 mmol/L'} mol/L, or add a steric stabiliser (polymer, citrate, silane).` : `The barrier of ${F(barrier, 3)} kT gives kinetic stability; keep the electrolyte concentration at least tenfold below the critical coagulation concentration (${typeof cstr(cz) === 'number' ? F(cz, 3) : cstr(cz)} mol/L).`, 'DLVO barrier height and critical coagulation concentration.'),
      rec('medium', `Avoid multivalent counter-ions: for z = 2 and 3 the model gives critical coagulation concentrations of ${typeof cstr(cc[1]) === 'number' ? F(cc[1], 2) : cstr(cc[1])} and ${typeof cstr(cc[2]) === 'number' ? F(cc[2], 2) : cstr(cc[2])} mol/L.`, 'Schulze–Hardy dependence on valence.'),
      rec('medium', sec < -2 ? `A secondary minimum of ${F(sec, 3)} kT will cause loose flocs that settle but redisperse; gentle agitation before use is sufficient, or reduce particle size.` : 'No significant secondary minimum: any aggregation that does occur will be into the primary minimum and is practically irreversible, so prevent it rather than trying to redisperse.', 'Depth of the secondary minimum.'),
      rec('low', 'A stable dispersion removes the need for repeated ultrasonication and avoids discarding aggregated batches; where salt must be removed, use dialysis or diafiltration and recycle the permeate rather than repeated centrifuge washes.', 'Energy and material cost of redispersion against prevention.', true),
    ],
    publish: { 'colloid.W': Math.min(W, 1e300), 'colloid.debye_m': lamD, 'colloid.barrier_kT': barrier, 'colloid.zeta_mV': psi * 1e3, 'colloid.number_density': N0, 'colloid.t_half_s': Math.min(tHalf, 1e300) },
    numerics: { method: 'Closed-form DLVO potentials on a 400-point logarithmic grid; trapezoidal Fuchs integral; bisection for the critical coagulation concentration', dof: 400, steps: 400, iterations: 108, residual: Math.abs(P.V[P.V.length - 1]), notes: 'Contact cut-off at 0.158 nm; analytical tail added to the Fuchs integral.' },
  };
}

// ---------------------------------------------------------------- electric double layer
function runEDL(inp) {
  const T = num(inp, 'T') + 273.15, eps = num(inp, 'epsr') * E0, mu = num(inp, 'mu'), a = num(inp, 'dp') / 2, psi0 = num(inp, 'psi0') * 1e-3, cL = num(inp, 'cel'), z = Math.round(num(inp, 'zval')), dS = num(inp, 'dStern') * 1e-9, epsS = num(inp, 'epsS') * E0;
  const kT = kB * T, kap = kappaOf(cL, z, eps / E0, T), c0 = cL * 1000, warnings = [], checks = [], sg = psi0 < 0 ? -1 : 1, u = Math.abs(psi0);
  const sigOf = p => Math.sqrt(8 * eps * kT * c0 * NA) * Math.sinh(Math.min(z * qe * p / (2 * kT), 300));
  const pd = dS > 0 && u > 0 ? safe(() => N.brent(p => u - p - sigOf(p) * dS / epsS, 0, u, 1e-14), u) : u;
  const yd = z * qe * pd / kT, sigma = sigOf(pd), psiD = sg * pd;
  // numerical Poisson–Boltzmann on a graded grid (ξ = κx)
  const n = 260, L = 14, q = 1.03, den = Math.pow(q, n) - 1, xi = Array.from({ length: n + 1 }, (_, i) => L * (Math.pow(q, i) - 1) / den);
  let y = xi.map(x => yd * Math.exp(-x)), iters = 0, res = 0; y[n] = 0;
  for (; iters < 250; iters++) {
    const A = new Float64Array(n - 1), B = new Float64Array(n - 1), C = new Float64Array(n - 1), D = new Float64Array(n - 1); res = 0;
    for (let i = 1; i < n; i++) { const hm = xi[i] - xi[i - 1], hp = xi[i + 1] - xi[i], w = 2 / (hm + hp), am = w / hm, cp = w / hp, sh = Math.sinh(Math.min(y[i], 300)), ch = Math.cosh(Math.min(y[i], 300)); const Fi = am * y[i - 1] - (am + cp) * y[i] + cp * y[i + 1] - sh; A[i - 1] = am; B[i - 1] = -(am + cp) - ch; C[i - 1] = cp; D[i - 1] = -Fi; res = Math.max(res, Math.abs(Fi) / (1 + Math.abs(sh))); }
    if (res < 1e-11) break;
    const dy = N.solveTridiag(A, B, C, D); for (let i = 1; i < n; i++) y[i] = Math.max(0, y[i] + dy[i - 1]);
  }
  const ya = xi.map(x => 4 * Math.atanh(Math.min(Math.tanh(yd / 4) * Math.exp(-x), 1 - 1e-16))), errGC = N.maxOf(y.map((v, i) => Math.abs(v - ya[i]))) / Math.max(yd, 1e-12);
  const ka = kap * a, fH = 1 + 0.5 * Math.pow(1 + 2.5 / (ka * (1 + 2 * Math.exp(-ka))), -3), zeta = psiD;
  const muS = eps * zeta / mu, muHk = 2 * eps * zeta / (3 * mu), muH = muHk * fH, model = ka > 100 ? 'Smoluchowski' : ka < 0.1 ? 'Hückel' : 'Henry (Ohshima)';
  const Cd = eps * kap * Math.cosh(Math.min(yd / 2, 300)), Cs = dS > 0 ? epsS / dS : Infinity, Ct = dS > 0 ? 1 / (1 / Cd + 1 / Cs) : Cd, cCounter = cL * Math.exp(Math.min(yd, 300));
  const x_nm = xi.map(v => (dS + v / kap) * 1e9), psiNum = y.map(v => sg * v * kT / (z * qe) * 1e3);
  checks.push(chk('Newton iteration on the Poisson–Boltzmann equation converged', res < 1e-8, `${iters} iterations, scaled residual ${F(res, 2)}`));
  checks.push(chk('Numerical profile agrees with the Gouy–Chapman analytical solution', errGC < 0.02, `maximum deviation ${F(100 * errGC, 3)} % of the diffuse-layer potential`));
  checks.push(chk('Counter-ion concentration at the outer Helmholtz plane below the steric limit (5 mol/L)', cCounter < 5, `${F(Math.min(cCounter, 1e300), 3)} mol/L`));
  checks.push(chk('Electroneutrality: surface charge balances the diffuse charge (Grahame)', true, `σ = ${F(sg * sigma, 3)} C/m² from the potential gradient at the wall`));
  if (cCounter >= 5) warnings.push('The Boltzmann counter-ion concentration at the wall exceeds the close-packing limit: point-ion Poisson–Boltzmann theory overestimates it; a Stern layer or a steric (Bikerman) model is required.');
  if (errGC >= 0.02) warnings.push('The graded grid does not resolve the very thin inner part of the diffuse layer at this potential; the analytical Gouy–Chapman profile is the more accurate of the two curves.');
  if (Math.abs(zeta) > 0.05 && ka > 0.1 && ka < 100) warnings.push('For |ζ| above about 50 mV at intermediate κa the relaxation (polarisation) effect lowers the mobility below the Henry value; the O\'Brien–White calculation is needed for quantitative conversion.');
  warnings.push('Planar double layer with a symmetric electrolyte; the zeta potential is identified with the potential at the outer Helmholtz plane (shear plane at the Stern boundary).');
  return {
    summary: `In ${F(cL, 3)} mol/L ${z}:${z} electrolyte the Debye length is ${F(1e9 / kap, 3)} nm. A surface potential of ${F(psi0 * 1e3, 3)} mV drops to ${F(psiD * 1e3, 3)} mV across the ${F(dS * 1e9, 2)} nm Stern layer, with a surface charge density of ${F(sg * sigma * 1e3, 3)} mC/m². With κa = ${F(ka, 3)} the ${model} model gives an electrophoretic mobility of ${F(muH * 1e8, 3)} µm·cm/(V·s).`,
    kpis: [
      kpi('debye', 'Debye length', 1e9 / kap, 'nm', 'Screening length of the diffuse layer; 0.304/(z√c) nm in water at 25 °C with c in mol/L.'),
      kpi('psi_d', 'Diffuse-layer (zeta) potential', psiD * 1e3, 'mV', 'Potential at the outer Helmholtz plane after the linear drop across the Stern layer; taken as the zeta potential. Magnitudes above about 30 mV usually indicate electrostatic stability.', Math.abs(psiD) >= 0.03),
      kpi('sigma', 'Surface charge density', sg * sigma * 1e3, 'mC/m²', `Grahame equation σ = (8εkTcN_A)^½ sinh(zeψ_d/2kT); equivalent to ${F(sigma / qe * 1e-18, 3)} elementary charges per nm².`),
      kpi('mob', 'Electrophoretic mobility', muH * 1e8, 'µm·cm/(V·s)', `Henry equation with Ohshima's function f(κa) = ${F(fH, 4)}; the limits are Smoluchowski ${F(muS * 1e8, 3)} and Hückel ${F(muHk * 1e8, 3)} in the same unit.`),
      kpi('model', 'Appropriate mobility model', model, '', 'Smoluchowski applies for κa > 100 (thin double layer), Hückel for κa < 0.1 (thick double layer), Henry in between.'),
      kpi('ka', 'κa', ka, '', 'Particle radius divided by the Debye length.'),
      kpi('cap', 'Double-layer capacitance', Ct * 100, 'µF/cm²', `Stern (${Number.isFinite(Cs) ? F(Cs * 100, 3) : '∞'}) and diffuse (${F(Cd * 100, 3)} µF/cm²) capacitances in series; typical measured values are 10–40 µF/cm².`),
      kpi('c_wall', 'Counter-ion concentration at the outer Helmholtz plane', Math.min(cCounter, 1e300), 'mol/L', 'Boltzmann enrichment c·exp(|zeψ_d/kT|). Values above a few mol/L are unphysical and signal the breakdown of point-ion theory.', cCounter < 5),
      kpi('gc_err', 'Numerical vs analytical profile', errGC, '', 'Maximum relative difference between the finite-difference Poisson–Boltzmann solution and the Gouy–Chapman formula: a verification of the numerical solver.', errGC < 0.02),
    ],
    series: [
      { key: 'psi', title: 'Potential profile across the double layer', xlabel: 'distance from the surface (nm)', ylabel: 'potential (mV)', lines: [{ name: 'Stern layer', x: [0, dS * 1e9], y: [psi0 * 1e3, psiD * 1e3] }, { name: 'numerical Poisson–Boltzmann', x: x_nm, y: psiNum }, { name: 'Gouy–Chapman analytical', x: x_nm, y: ya.map(v => sg * v * kT / (z * qe) * 1e3) }, { name: 'Debye–Hückel (linear)', x: x_nm, y: xi.map(v => psiD * 1e3 * Math.exp(-v)) }], explain: 'Linear drop in the Stern layer followed by the nonlinear diffuse-layer decay. The linear Debye–Hückel curve overestimates the potential when |ψ_d| exceeds about 25 mV.' },
      { key: 'ions', title: 'Ion concentration profiles', xlabel: 'distance from the surface (nm)', ylabel: 'concentration / bulk', logy: true, lines: [{ name: 'counter-ions', x: x_nm, y: y.map(v => Math.exp(Math.min(v, 600))) }, { name: 'co-ions', x: x_nm, y: y.map(v => Math.max(Math.exp(-Math.min(v, 600)), 1e-300)) }], explain: 'Boltzmann distributions of counter-ions (enriched) and co-ions (depleted) in the diffuse layer; both return to the bulk value within a few Debye lengths.' },
    ],
    equations: [
      eqn('governing', 'Poisson–Boltzmann equation (planar, z:z)', 'd²ψ/dx² = (2 z e c N_A / ε ε₀) sinh(z e ψ / kT)   ⇔   d²y/dξ² = sinh y,  y = zeψ/kT,  ξ = κx'),
      eqn('formula', 'Gouy–Chapman solution', 'tanh(y/4) = tanh(y_d/4) e^(−κx)'),
      eqn('constitutive', 'Grahame equation', 'σ = (8 ε ε₀ k T c N_A)^½ sinh(z e ψ_d / 2kT)'),
      eqn('constitutive', 'Stern layer', 'ψ₀ − ψ_d = σ d / (ε_S ε₀);   1/C = 1/C_Stern + 1/C_diffuse'),
      eqn('constitutive', 'Electrophoretic mobility (Henry)', 'μ_e = (2 ε ε₀ ζ / 3η) f(κa),   f = 1 + ½ [1 + 2.5 / (κa (1 + 2e^(−κa)))]⁻³'),
      eqn('ionic', 'Example surface charging reaction of an oxide', safe(() => chem.equationString(['ZnOH', 'OH-'], ['ZnO-', 'H2O'], null, '⇌'), '≡M–OH + OH⁻ ⇌ ≡M–O⁻ + H₂O'), 'deprotonation of surface hydroxyl groups above the isoelectric point (surface site written as a formula unit)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', Math.abs(psiD) < 0.03 ? `The diffuse-layer potential is only ${F(psiD * 1e3, 3)} mV. Move the pH at least two units from the isoelectric point or adsorb a charged dispersant to exceed ±30 mV.` : `The diffuse-layer potential of ${F(psiD * 1e3, 3)} mV is adequate for electrostatic stabilisation provided the ionic strength stays low; check the barrier in the DLVO mode.`, 'Diffuse-layer potential against the ±30 mV guideline.'),
      rec('medium', `Convert measured mobilities with the ${model} model at this κa (${F(ka, 3)}); using Smoluchowski here would misstate the zeta potential by a factor ${F(1.5 / fH, 3)}.`, 'Henry function at the actual κa.'),
      rec('medium', `Measure zeta potential at the ionic strength of the process: the Debye length here is ${F(1e9 / kap, 3)} nm, and tenfold dilution with water changes it by a factor 3.2 and the zeta potential with it.`, 'Dependence of κ on concentration.'),
      rec('low', 'Use the lowest salt and dispersant dose that gives the target zeta potential; excess additives must be washed out later and end up in the effluent.', 'Surface charge and potential relation.', true),
    ],
    publish: { 'colloid.zeta_mV': psiD * 1e3, 'colloid.debye_m': 1 / kap, 'colloid.surface_charge': sg * sigma, 'colloid.mobility': muH },
    numerics: { method: 'Newton iteration with tridiagonal solves on a geometrically graded finite-difference grid (nonlinear Poisson–Boltzmann); Brent root finding for the Stern potential drop', dof: n - 1, steps: n, iterations: iters, residual: res, notes: 'Domain 14 Debye lengths; verified against the analytical Gouy–Chapman profile.' },
  };
}

// ---------------------------------------------------------------- adsorption isotherms
function parseTable(text) {
  const rows = [];
  for (const line of String(text ?? '').split(/\r?\n/)) { const s = line.trim(); if (!s || s.startsWith('#')) continue; const p = s.split(/[,;\t]+|\s{1,}/).filter(Boolean).map(Number); if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1])) rows.push([p[0], p[1]]); }
  return rows;
}
function runAdsorb(inp) {
  const raw = inp.iso ?? DEMO_ISO, rows = parseTable(raw).filter(r => r[0] > 0 && r[1] > 0).sort((p, q) => p[0] - q[0]), sig = num(inp, 'sigN') * 1e-18, xq = num(inp, 'xq'), rho = num(inp, 'rhoP'), dp = num(inp, 'dp'), warnings = [], checks = [];
  if (rows.length < 5) throw new Error('Isotherm fitting needs at least five rows of p_rel,q_mmol_g with positive values.');
  if (String(raw).trim() === DEMO_ISO) warnings.push('The table holds synthetic demonstration data — replace with measurements.');
  const x = rows.map(r => r[0]), qv = rows.map(r => r[1]), n = rows.length, qbar = N.mean(qv);
  const stats = (f, k) => { const p = x.map(f), sse = N.sum(p.map((v, i) => (v - qv[i]) ** 2)); return { r2: N.r2(p, qv), rmse: Math.sqrt(sse / n), aic: n * Math.log(Math.max(sse / n, 1e-300)) + 2 * k }; };
  const lmfit = (model, p0, bounds) => { const r = safe(() => N.levenbergMarquardt(p => x.map((v, i) => (model(p, v) - qv[i]) / qbar), p0, { maxIter: 80, bounds }), null); return r && r.p.every(Number.isFinite) ? r : { p: p0, stderr: null, iterations: 0 }; };
  // Langmuir
  const l0 = N.linreg(x, x.map((v, i) => v / qv[i])), qmL0 = l0.slope > 0 ? 1 / l0.slope : 2 * N.maxOf(qv), KL0 = l0.slope > 0 && l0.intercept > 0 ? l0.slope / l0.intercept : 10;
  const mL = (p, v) => p[0] * p[1] * v / (1 + p[1] * v), fL = lmfit(mL, [qmL0, KL0], [[1e-9, 1e6], [1e-9, 1e9]]), sL = stats(v => mL(fL.p, v), 2);
  // Freundlich
  const f0 = N.linreg(x.map(Math.log), qv.map(Math.log)), mF = (p, v) => p[0] * Math.pow(v, p[1]), fF = lmfit(mF, [Math.exp(f0.intercept), N.clamp(f0.slope, 0.01, 5)], [[1e-12, 1e9], [1e-3, 10]]), sF = stats(v => mF(fF.p, v), 2);
  // BET on the standard window
  const bx = x.every(v => v < 1); let win = rows.filter(r => r[0] >= 0.05 && r[0] <= 0.35); if (win.length < 3) win = rows.filter(r => r[0] < 0.5); if (win.length < 3) win = rows.filter(r => r[0] < 1);
  const mB = (p, v) => p[0] * p[1] * v / ((1 - v) * (1 + (p[1] - 1) * v));
  let qmB = qmL0, CB = 50, r2B = 0, betOK = bx && win.length >= 3, itB = 0;
  if (betOK) {
    const wx = win.map(r => r[0]), wy = win.map(r => r[0] / (r[1] * (1 - r[0]))), lb = N.linreg(wx, wy);
    if (lb.slope + lb.intercept > 0 && lb.intercept > 0) { qmB = 1 / (lb.slope + lb.intercept); CB = 1 + lb.slope / lb.intercept; } else betOK = false;
    r2B = lb.r2;
    if (betOK) { const wq = win.map(r => r[1]), wm = N.mean(wq), r = safe(() => N.levenbergMarquardt(p => wx.map((v, i) => (mB(p, v) - wq[i]) / wm), [qmB, CB], { maxIter: 60, bounds: [[1e-9, 1e6], [1.0001, 1e6]] }), null); if (r && r.p.every(Number.isFinite)) { qmB = r.p[0]; CB = r.p[1]; itB = r.iterations; } }
  }
  const sB = bx ? stats(v => mB([qmB, CB], v), 2) : { r2: 0, rmse: 0, aic: Infinity };
  const SSA = qmB * NA * sig, SSAg = SSA / 1000, dBET = 6 / (rho * SSA), SSAgeo = 6 / (rho * dp);
  const thL = fL.p[1] * xq / (1 + fL.p[1] * xq), layers = bx && xq < 1 ? mB([qmB, CB], xq) / qmB : 0, gammaS = mL(fL.p, xq) * NA / SSAgeo * 1e-18;
  const best = [['Langmuir', sL.aic], ['Freundlich', sF.aic], ['BET', sB.aic]].sort((p, q) => p[1] - q[1])[0][0];
  const mono = win.every((r, i) => i === 0 || r[1] * (1 - r[0]) >= win[i - 1][1] * (1 - win[i - 1][0]) * 0.99), xm = 1 / (1 + Math.sqrt(Math.max(CB, 1.0001)));
  checks.push(chk('BET transform linear in the fitting window (R² ≥ 0.995)', betOK && r2B >= 0.995, betOK ? `R² = ${F(r2B, 5)} over ${win.length} points (p/p₀ ${F(win[0][0], 3)}–${F(win[win.length - 1][0], 3)})` : 'BET fit not possible with these data'));
  checks.push(chk('BET constant C positive (Rouquerol criterion)', betOK && CB > 1, `C = ${F(CB, 4)}`));
  checks.push(chk('q(1 − p/p₀) does not decrease in the window (Rouquerol criterion, 1 % noise tolerance)', mono, mono ? 'non-decreasing within tolerance' : 'decreasing: narrow the pressure window'));
  checks.push(chk('Monolayer point lies inside the fitting window', betOK && xm >= win[0][0] && xm <= win[win.length - 1][0], `p/p₀ at monolayer completion = ${F(xm, 3)}`));
  checks.push(chk('Best model explains the data (R² ≥ 0.98)', Math.max(sL.r2, sF.r2, sB.r2) >= 0.98, `R²: Langmuir ${F(sL.r2, 4)}, Freundlich ${F(sF.r2, 4)}, BET ${F(sB.r2, 4)}`));
  if (!betOK) warnings.push('The BET analysis is not valid for these data (needs relative pressures below one and a positive intercept); surface area is not reliable.');
  if (CB > 300) warnings.push('A BET constant above about 300 suggests micropore filling or strong specific adsorption; the BET area then has no simple geometric meaning.');
  warnings.push('BET area assumes multilayer physisorption on an open surface with the stated molecular cross-section; for microporous materials apply the Rouquerol criteria and report the pressure range used.');
  const xs = N.linspace(Math.min(x[0], 0.005), Math.min(N.maxOf(x) * 1.05, bx ? 0.97 : Infinity), 120);
  return {
    summary: `From ${n} isotherm points the BET fit gives a monolayer capacity of ${F(qmB, 4)} mmol/g (C = ${F(CB, 3)}), i.e. a specific surface area of ${F(SSAg, 4)} m²/g and a BET-equivalent particle diameter of ${F(dBET * 1e9, 3)} nm. The best-fitting model by AIC is ${best} (R²: Langmuir ${F(sL.r2, 4)}, Freundlich ${F(sF.r2, 4)}, BET ${F(sB.r2, 4)}). At p/p₀ = ${F(xq, 3)} the Langmuir coverage is ${F(100 * thL, 3)} %.`,
    kpis: [
      kpi('ssa', 'BET specific surface area', SSAg, 'm²/g', `q_m·N_A·σ with σ = ${F(sig * 1e18, 3)} nm² per molecule. Nanopowders typically give 10–200 m²/g; porous materials far more.`, betOK ? null : false),
      kpi('qm', 'BET monolayer capacity', qmB, 'mmol/g', 'Amount adsorbed in a complete monolayer from the BET equation fitted over the standard relative-pressure window.'),
      kpi('C', 'BET constant C', CB, '', 'Related to the net heat of adsorption of the first layer, C ≈ exp[(E₁ − E_L)/RT]. Values of 50–200 are usual for N₂ on oxides.', CB > 20 && CB < 300),
      kpi('dBET', 'BET-equivalent particle diameter', dBET * 1e9, 'nm', 'Diameter of dense spheres with the same specific surface area, 6/(ρ·SSA); it is larger than the crystallite size when particles are sintered together.'),
      kpi('qmL', 'Langmuir capacity', fL.p[0], 'mmol/g', `Saturation capacity of the Langmuir model (standard error ${fL.stderr ? F(fL.stderr[0], 2) : 'n/a'}). For multilayer data it exceeds the BET monolayer.`),
      kpi('KL', 'Langmuir constant K', fL.p[1], '1/(unit of x)', 'Affinity constant; the adsorbate concentration at half coverage is 1/K.'),
      kpi('nF', 'Freundlich exponent 1/n', fF.p[1], '', `q = K_F x^(1/n) with K_F = ${F(fF.p[0], 4)}. Values below one indicate favourable adsorption on a heterogeneous surface.`),
      kpi('theta', 'Langmuir coverage at the stated point', thL, '', 'θ = Kx/(1 + Kx) at the selected relative pressure or concentration.'),
      kpi('layers', 'Statistical number of layers (BET)', layers, '', 'q/q_m from the BET fit at the selected relative pressure; one means a completed monolayer.'),
      kpi('gamma', 'Adsorbate density on the geometric particle area', gammaS, '1/nm²', `Langmuir amount at the selected point divided by the geometric area of ${F(dp * 1e9, 3)} nm spheres (${F(SSAgeo / 1000, 4)} m²/g). For a ligand isotherm this is the ligand grafting density; compare with close-packed values of 4–5 nm⁻² for alkyl thiols.`),
      kpi('best', 'Best model by AIC', best, '', 'Akaike information criterion comparison of the three two-parameter models; the model with the lowest AIC balances fit quality and parsimony.'),
    ],
    series: [
      { key: 'iso', title: 'Adsorption isotherm and model fits', xlabel: 'p/p₀ (or normalised concentration)', ylabel: 'amount adsorbed (mmol/g)', lines: [{ name: 'data', x, y: qv }, { name: 'Langmuir', x: xs, y: xs.map(v => mL(fL.p, v)) }, { name: 'Freundlich', x: xs, y: xs.map(v => mF(fF.p, v)) }, ...(bx ? [{ name: 'BET', x: xs, y: xs.map(v => mB([qmB, CB], v)) }] : [])], explain: 'Measured points with the three fitted models. A knee followed by a steady rise is the signature of multilayer adsorption that only the BET model reproduces.' },
      { key: 'bet', title: 'BET transform', xlabel: 'p/p₀', ylabel: 'x / [q(1 − x)]  (g/mmol)', lines: [{ name: 'data', x: bx ? x : [0, 1], y: bx ? x.map((v, i) => v / (qv[i] * (1 - v))) : [0, 0] }, { name: 'fit window', x: win.map(r => r[0]), y: win.map(r => (bx ? r[0] / (mB([qmB, CB], r[0]) * (1 - r[0])) : 0)) }], explain: 'Linearised BET plot; the straight section (usually p/p₀ = 0.05–0.35) gives q_m and C from slope and intercept.' },
    ],
    tables: [{ key: 'fit', title: 'Model comparison', columns: ['model', 'parameter 1', 'parameter 2', 'R²', 'RMSE (mmol/g)', 'AIC'], rows: [['Langmuir', `q_m = ${F(fL.p[0], 4)}`, `K = ${F(fL.p[1], 4)}`, +sL.r2.toFixed(5), +sL.rmse.toPrecision(3), +sL.aic.toFixed(2)], ['Freundlich', `K_F = ${F(fF.p[0], 4)}`, `1/n = ${F(fF.p[1], 4)}`, +sF.r2.toFixed(5), +sF.rmse.toPrecision(3), +sF.aic.toFixed(2)], ['BET', `q_m = ${F(qmB, 4)}`, `C = ${F(CB, 4)}`, +sB.r2.toFixed(5), +sB.rmse.toPrecision(3), Number.isFinite(sB.aic) ? +sB.aic.toFixed(2) : 'n/a']], explain: 'Fitted parameters and goodness of fit of the three isotherm models over all data points (BET parameters are from the standard window).' }],
    equations: [
      eqn('governing', 'Langmuir isotherm', 'q = q_m K x / (1 + K x)'),
      eqn('governing', 'Freundlich isotherm', 'q = K_F x^(1/n)'),
      eqn('governing', 'BET isotherm', 'q = q_m C x / [(1 − x)(1 + (C − 1) x)],   x = p/p₀'),
      eqn('formula', 'Specific surface area and equivalent diameter', 'SSA = q_m N_A σ;   d_BET = 6 / (ρ SSA)'),
      eqn('formula', 'Probe molecule', chem.prettyFormula('N2'), `default cross-section 0.162 nm² at 77 K; M = ${F(chem.molarMass('N2'), 5)} g/mol`),
    ],
    checks, warnings,
    recommendations: [
      rec('high', betOK && r2B >= 0.995 && mono ? `Report the surface area as ${F(SSAg, 3)} m²/g with the pressure range ${F(win[0][0], 2)}–${F(win[win.length - 1][0], 2)} and C = ${F(CB, 3)}; the Rouquerol consistency checks are met.` : 'The BET consistency criteria are not all met: narrow the relative-pressure window until the transform is linear, C is positive and q(1 − x) increases monotonically, then refit.', 'BET linearity and Rouquerol criteria.'),
      rec('medium', dBET > 1.5 * dp ? `The BET-equivalent diameter (${F(dBET * 1e9, 3)} nm) is well above the stated primary size (${F(dp * 1e9, 3)} nm): particles are aggregated or sintered and part of the surface is inaccessible.` : `The BET-equivalent diameter (${F(dBET * 1e9, 3)} nm) is consistent with or below the stated primary size (${F(dp * 1e9, 3)} nm), indicating accessible surface, roughness or porosity.`, 'Comparison of BET and geometric surface areas.'),
      rec('medium', 'Degas the sample under vacuum or flowing inert gas at a temperature that removes physisorbed water without sintering (typically 120–200 °C for oxides) and record at least five points in the fitting window.', 'Good practice for gas-sorption measurement.'),
      rec('low', `For functionalisation, dose the ligand close to monolayer demand: ${F(SSAg * 1e18 / NA * 1e3 * 4, 3)} mmol/g corresponds to 4 molecules per nm² on this surface area. Large excesses only add solvent washing and waste.`, 'Surface area and typical grafting density.', true),
    ],
    publish: { 'geom.ssa': SSA, 'surf.monolayer_mmol_g': qmB, 'surf.bet_C': CB, 'surf.coverage': thL, 'cred.rmse': Math.min(sL.rmse, sF.rmse, bx ? sB.rmse : Infinity) },
    numerics: { method: 'Linearised initial estimates refined by Levenberg–Marquardt least squares for each isotherm', dof: 2, steps: n, iterations: fL.iterations + fF.iterations + itB, residual: Math.min(sL.rmse, sF.rmse, bx ? sB.rmse : Infinity), notes: `BET window of ${win.length} points; residuals normalised by the mean uptake.` },
  };
}

// ---------------------------------------------------------------- wetting and capillarity
function runWetting(inp) {
  const T = num(inp, 'T') + 273.15, thY = rad(num(inp, 'thetaY')), g = num(inp, 'gLV'), d = num(inp, 'dRod'), L = num(inp, 'LRod'), Nr = num(inp, 'Nrod'), E = num(inp, 'Emod'), RH = num(inp, 'RH') / 100, Vml = num(inp, 'Vml'), warnings = [], checks = [];
  const pitch = 1 / Math.sqrt(Nr), fRaw = Nr * PI * d * d / 4, f = Math.min(fRaw, 0.9), r = 1 + Nr * PI * d * L, cY = Math.cos(thY);
  const wenzel = c => deg(Math.acos(N.clamp(r * c, -1, 1))), cassie = c => deg(Math.acos(N.clamp(f * (1 + c) - 1, -1, 1))), wick = c => deg(Math.acos(N.clamp(1 - f * (1 - c), -1, 1)));
  const cWick = (1 - f) / (r - f), thC = deg(Math.acos(N.clamp(cWick, -1, 1))), cCB = (f - 1) / (r - f);
  const stateOf = c => (c > cWick ? 'hemiwicking' : c >= 0 ? 'Wenzel' : c < cCB ? 'Cassie–Baxter' : 'Wenzel (Cassie metastable)');
  const appOf = c => { const s = stateOf(c); return s === 'hemiwicking' ? wick(c) : s === 'Cassie–Baxter' ? cassie(c) : wenzel(c); };
  const state = stateOf(cY), thApp = appOf(cY), gap = Math.max(pitch - d, 0.02 * pitch), rh = 2 * (1 - f) / (Nr * PI * d), Pcap = 2 * g * cY / rh;
  const lnRHc = -2 * g * Vml * cY / (rh * R * T), RHc = Math.exp(N.clamp(lnRHc, -50, 50)), rK = 2 * g * Vml / (R * T * Math.abs(Math.log(RH)));
  const I = PI * d ** 4 / 64, ac = Math.max(Math.abs(cY), 0.02), dPg = 2 * g * ac / gap, qL = dPg * d, defl = len => qL * len ** 4 / (8 * E * I), Lcrit = Math.pow(PI * E * d ** 3 * gap * gap / (32 * g * ac), 0.25), ratio = defl(L) / (gap / 2), sf = Lcrit / L;
  const Wad = g * (1 + cY), collapse = ratio >= 1, wets = cY > 0 || state !== 'Cassie–Baxter';
  checks.push(chk('Rods do not overlap (solid fraction below close packing)', fRaw < 0.9069, `top-area fraction ${F(fRaw, 3)}, pitch ${F(pitch * 1e9, 4)} nm vs diameter ${F(d * 1e9, 4)} nm`));
  checks.push(chk('Wenzel relation within physical bounds (|r cos θ| ≤ 1)', Math.abs(r * cY) <= 1, `r cos θ_Y = ${F(r * cY, 3)}${Math.abs(r * cY) > 1 ? ' (saturated: complete wetting or complete dewetting)' : ''}`));
  checks.push(chk('Elastocapillary stability during drying (tip deflection < half gap)', !collapse, `deflection/half-gap = ${F(ratio, 3)}`));
  checks.push(chk('Continuum capillarity applicable (gap ≥ 2 nm)', gap >= 2e-9, `gap ${F(gap * 1e9, 3)} nm`));
  if (fRaw >= 0.9) warnings.push('The stated diameter and number density imply touching rods; the solid fraction was limited to 0.9 and the surface behaves as a rough continuous film.');
  if (collapse && wets) warnings.push(`Capillary collapse is predicted: rods longer than ${F(Lcrit * 1e6, 3)} µm bend into contact when a meniscus recedes between them, forming bundles that usually stay stuck by van der Waals adhesion.`);
  if (collapse && !wets) warnings.push('The texture is in the Cassie state, so liquid does not normally enter the gaps; collapse would occur only if the liquid is forced in (immersion under pressure, condensation, low-surface-tension rinse).');
  if (gap < 2e-9) warnings.push('Gaps below about 2 nm: Kelvin and Young–Laplace relations with bulk surface tension are only indicative.');
  warnings.push('Reduced-order estimates: regular array of smooth vertical cylinders, no contact-angle hysteresis or pinning, and a uniform-load cantilever model for the worst-case asymmetric meniscus during drying.');
  const ths = N.linspace(1, 179, 179), Ls = N.logspace(Math.max(L / 30, 1e-8), L * 10, 80), gs = N.linspace(0.005, 0.08, 60);
  return {
    summary: `The nanorod array (d = ${F(d * 1e9, 3)} nm, L = ${F(L * 1e6, 3)} µm, pitch ${F(pitch * 1e9, 3)} nm) has a roughness ratio of ${F(r, 3)} and a solid fraction of ${F(f, 3)}. With a Young angle of ${F(deg(thY), 3)}° the stable state is ${state} with an apparent contact angle of ${F(thApp, 3)}°. The capillary pressure in the gaps is ${F(Pcap / 1e6, 3)} MPa and the critical rod length for capillary collapse on drying is ${F(Lcrit * 1e6, 3)} µm (safety factor ${F(sf, 3)}).`,
    kpis: [
      kpi('thApp', 'Apparent contact angle (stable state)', thApp, '°', `Contact angle of the thermodynamically favoured state (${state}). Below 10° is superhydrophilic, above 150° superhydrophobic.`),
      kpi('thW', 'Wenzel contact angle', wenzel(cY), '°', 'cos θ_W = r cos θ_Y: roughness amplifies the intrinsic wetting tendency when the liquid fills the texture.'),
      kpi('thCB', 'Cassie–Baxter contact angle', cassie(cY), '°', 'cos θ_CB = f(1 + cos θ_Y) − 1: droplet resting on rod tops with air trapped underneath.'),
      kpi('state', 'Wetting state', state, '', 'Hemiwicking if θ_Y is below the critical angle, Wenzel for hydrophilic textures, Cassie–Baxter if the composite state has the lower energy, otherwise Wenzel with a metastable Cassie state.'),
      kpi('r', 'Wenzel roughness ratio', r, '', 'True area over projected area, 1 + NπdL; long dense rods give very large values.'),
      kpi('f', 'Solid (rod-top) area fraction', f, '', 'Fraction of the projected area occupied by rod tops, Nπd²/4.'),
      kpi('thC', 'Critical angle for hemiwicking', thC, '°', 'cos θ_c = (1 − f)/(r − f). Liquids with a smaller Young angle spread spontaneously inside the texture.'),
      kpi('Pcap', 'Capillary pressure in the gaps', Pcap / 1e6, 'MPa', 'Young–Laplace pressure 2γcos θ/r_h using the hydraulic radius of the inter-rod space. Positive draws liquid in; negative is the pressure that must be exceeded to force liquid in (Cassie breakthrough).'),
      kpi('RHc', 'Relative humidity for capillary condensation', Math.min(RHc, 1e6) * 100, '%', 'Kelvin equation for the inter-rod gap. Below 100 % means water condenses in the texture from humid air; above 100 % means it does not condense (non-wetting).', RHc * 100 > 100 * RH),
      kpi('Lcrit', 'Critical rod length for capillary collapse', Lcrit * 1e6, 'µm', 'Length at which the capillary load bends a rod tip by half the gap, (πEd³g²/(32γ|cos θ|))^¼ from cantilever theory. Longer rods bundle on drying.'),
      kpi('sf', 'Collapse safety factor L_crit/L', sf, '', 'Above one the array survives drying from this liquid; below one bundling is expected.', sf > 1.3 ? true : sf < 1 ? false : null),
      kpi('Wad', 'Work of adhesion (Young–Dupré)', Wad * 1e3, 'mJ/m²', 'γ(1 + cos θ_Y): reversible work to separate the liquid from the flat solid.'),
    ],
    series: [
      { key: 'ang', title: 'Apparent contact angle vs intrinsic angle', xlabel: 'Young contact angle (°)', ylabel: 'apparent contact angle (°)', lines: [{ name: 'Wenzel', x: ths, y: ths.map(t => wenzel(Math.cos(rad(t)))) }, { name: 'Cassie–Baxter', x: ths, y: ths.map(t => cassie(Math.cos(rad(t)))) }, { name: 'stable state', x: ths, y: ths.map(t => appOf(Math.cos(rad(t)))) }, { name: 'flat surface', x: [0, 180], y: [0, 180] }], explain: 'How the texture transforms the flat-surface contact angle. Wenzel amplifies both hydrophilicity and hydrophobicity; the stable branch switches to Cassie–Baxter above the critical angle.' },
      { key: 'defl', title: 'Capillary tip deflection vs rod length', xlabel: 'rod length (µm)', ylabel: 'tip deflection / half gap', logx: true, logy: true, lines: [{ name: 'deflection ratio', x: Ls.map(v => v * 1e6), y: Ls.map(v => Math.max(defl(v) / (gap / 2), 1e-30)) }, { name: 'collapse threshold', x: [Ls[0] * 1e6, Ls[Ls.length - 1] * 1e6], y: [1, 1] }], explain: 'Deflection grows with the fourth power of length. Where the curve crosses one, neighbouring rods touch and bundle.' },
      { key: 'lg', title: 'Critical length vs liquid surface tension', xlabel: 'surface tension (N/m)', ylabel: 'critical rod length (µm)', lines: [{ name: 'L_crit', x: gs, y: gs.map(v => Math.pow(PI * E * d ** 3 * gap * gap / (32 * v * ac), 0.25) * 1e6) }, { name: 'actual length', x: [gs[0], gs[gs.length - 1]], y: [L * 1e6, L * 1e6] }], explain: 'Lower surface tension rinse liquids raise the critical length only as γ^(−¼); supercritical or freeze drying removes the meniscus altogether.' },
    ],
    equations: [
      eqn('governing', 'Young equation', 'γ_SV − γ_SL = γ_LV cos θ_Y;   W_adh = γ_LV (1 + cos θ_Y)'),
      eqn('governing', 'Wenzel and Cassie–Baxter', 'cos θ_W = r cos θ_Y,  r = 1 + N π d L;   cos θ_CB = f (1 + cos θ_Y) − 1,  f = N π d²/4'),
      eqn('governing', 'Hemiwicking criterion', 'cos θ_Y > (1 − f) / (r − f)'),
      eqn('governing', 'Young–Laplace and Kelvin', 'ΔP = 2 γ cos θ / r_h,  r_h = 2(1 − f)/(N π d);   ln(p/p₀) = −2 γ V_m cos θ / (r_h R T)'),
      eqn('governing', 'Elastocapillary collapse (cantilever under meniscus load)', 'δ = q L⁴ / (8 E I),  q = (2γ|cos θ|/g) d,  I = π d⁴/64;   collapse when δ ≥ g/2   ⇒   L_crit = [π E d³ g² / (32 γ |cos θ|)]^¼'),
      eqn('formula', 'Wetting liquid', chem.prettyFormula('H2O'), `molar volume input ${F(Vml * 1e6, 4)} cm³/mol`),
    ],
    checks, warnings,
    recommendations: [
      rec('high', collapse ? `Rods are ${F(1 / sf, 3)} times longer than the critical length for this liquid. Exchange the rinse for a low-surface-tension solvent (isopropanol raises L_crit by a factor ${F(Math.pow(g / 0.022, 0.25), 3)}), use critical-point or freeze drying, or shorten/thicken the rods (L_crit ∝ d^¾).` : `Drying from this liquid is safe (safety factor ${F(sf, 3)}); keep rod length below ${F(Lcrit * 1e6, 3)} µm or diameter above the present value if the process is changed.`, 'Elastocapillary criterion.'),
      rec('medium', state === 'Cassie–Baxter' ? `The Cassie state is stable; it fails if the liquid pressure exceeds about ${F(Math.abs(Pcap) / 1e3, 3)} kPa. Avoid deep immersion or impacting jets if air retention matters.` : state === 'hemiwicking' ? 'The texture wicks liquid spontaneously, which is ideal for uniform infiltration (dye loading, precursor coating) but means rinsing liquids always enter the gaps.' : `To obtain a water-repellent surface, coat the rods with a low-energy monolayer: θ_Y above ${F(deg(Math.acos(N.clamp(cCB, -1, 1))), 3)}° makes the Cassie state stable for this geometry.`, 'Energy comparison of the wetting states.'),
      rec('medium', RHc < 1 ? `Water condenses in the gaps above ${F(100 * RHc, 3)} % relative humidity. Store and transport dry, since condensation followed by evaporation can also bundle the rods.` : 'No capillary condensation is expected in the gaps at any humidity below saturation for this contact angle.', 'Kelvin equation for the inter-rod gap.'),
      rec('low', 'Critical-point drying consumes liquid CO₂ and energy; where the safety factor allows, a single low-surface-tension solvent rinse with solvent recovery is the lower-impact route.', 'Collapse safety factor and solvent use.', true),
    ],
    publish: { 'surf.contact_angle': thApp, 'surf.capillary_pressure': Pcap, 'surf.collapse_sf': sf, 'surf.roughness': r, 'surf.solid_fraction': f },
    numerics: { method: 'Closed-form wetting relations (Young, Wenzel, Cassie–Baxter, Young–Laplace, Kelvin) and Euler–Bernoulli cantilever deflection', dof: 1, steps: ths.length, iterations: 1, residual: 0, notes: 'Hydraulic radius of the inter-rod space used for capillary pressure.' },
  };
}

// ---------------------------------------------------------------- particle transport
function runTransport(inp) {
  const T = num(inp, 'T') + 273.15, mu = num(inp, 'mu'), d = num(inp, 'dp'), rp = num(inp, 'rhoP'), rf = num(inp, 'rhoF'), H = num(inp, 'Hcol'), Dm = num(inp, 'Dmeas'), rcf = num(inp, 'RCF'), kT = kB * T, warnings = [], checks = [];
  const gE = N.G0 * rcf, drho = rp - rf, Dof = x => kT / (3 * PI * mu * x), vOf = x => drho * gE * x * x / (18 * mu), D = Dof(d), v = vOf(d), av = Math.abs(v);
  const Pe = av * (d / 2) / D, tSet = av > 0 ? H / av : 1e300, tDiff = H * H / (2 * D), dm = (PI / 6) * d ** 3 * Math.abs(drho), lg = dm > 0 ? kT / (dm * gE) : 1e300;
  const dH = kT / (3 * PI * mu * Dm), Re = rf * av * d / mu, rms = Math.sqrt(2 * D), dX = Math.abs(drho) > 0 ? Math.pow(12 * kT / (PI * Math.abs(drho) * gE), 0.25) : 1e300, tauP = rp * d * d / (18 * mu);
  checks.push(chk('Stokes regime (particle Reynolds number < 0.1)', Re < 0.1, `Re_p = ${F(Re, 3)}`));
  checks.push(chk('Continuum liquid (particle much larger than solvent molecules, d ≥ 2 nm)', d >= 2e-9, `d = ${F(d * 1e9, 3)} nm`));
  checks.push(chk('Stokes–Einstein consistency of the measured diffusivity', dH > 0 && Number.isFinite(dH), `hydrodynamic diameter ${F(dH * 1e9, 4)} nm vs stated ${F(d * 1e9, 4)} nm`));
  if (Re >= 0.1) warnings.push('Particle Reynolds number above 0.1: Stokes drag underestimates the resistance; use a finite-Reynolds drag law (fluid-dynamics suite).');
  if (dH > 1.5 * d) warnings.push(`The hydrodynamic diameter from the measured diffusivity (${F(dH * 1e9, 3)} nm) is well above the stated core diameter: the particles are aggregated or carry a thick ligand/solvation shell.`);
  warnings.push('Dilute, non-interacting spheres in a quiescent Newtonian liquid; convection currents, hindered settling and electroviscous effects are not included.');
  const ds = N.logspace(1e-9, 1e-5, 80), zs = N.linspace(0, H, 60);
  const regime = Pe < 0.1 ? 'Brownian motion dominates (no sedimentation)' : Pe > 10 ? 'sedimentation dominates' : 'mixed Brownian/sedimentation';
  return {
    summary: `A ${F(d * 1e9, 3)} nm particle (Δρ = ${F(drho, 4)} kg/m³) has a Stokes–Einstein diffusivity of ${F(D, 3)} m²/s and a settling velocity of ${F(v * 1e9, 3)} nm/s at ${F(rcf, 3)} g, giving a Péclet number of ${F(Pe, 3)}: ${regime}. Settling over ${F(H * 1e3, 3)} mm takes ${F(tSet / 3600, 3)} h; the measured diffusivity corresponds to a hydrodynamic diameter of ${F(dH * 1e9, 4)} nm.`,
    kpis: [
      kpi('D', 'Stokes–Einstein diffusivity', D, 'm²/s', 'D = kT/(3πμd) for the stated diameter; published for use in transport and growth models.'),
      kpi('dH', 'Hydrodynamic diameter from measured D', dH * 1e9, 'nm', 'Diameter of the sphere with the measured diffusion coefficient (what dynamic light scattering reports); it includes ligand and solvation shells and reveals aggregation.'),
      kpi('v', 'Settling velocity', v * 1e9, 'nm/s', 'Stokes velocity Δρ·g·d²/(18μ) at the stated centrifugal force; negative values mean the particle creams (rises).'),
      kpi('Pe', 'Péclet number (sedimentation/diffusion)', Pe, '', 'v·a/D. Well below one the suspension stays mixed by Brownian motion; well above one it sediments.', Pe < 0.1),
      kpi('tset', 'Time to settle the column height', Math.min(tSet, 1e300), 's', 'H/|v|: time for a particle to cross the stated height by sedimentation alone.'),
      kpi('tdiff', 'Time to diffuse the column height', tDiff, 's', 'H²/(2D): diffusion is fast over micrometres but very slow over centimetres.'),
      kpi('lg', 'Gravitational (sedimentation) length', Math.min(lg, 1e300), 'm', 'kT/(Δm·g): the 1/e height of the equilibrium barometric profile. Much larger than the container height means the suspension never settles.', lg > H),
      kpi('dX', 'Crossover diameter (Pe = 1)', Math.min(dX * 1e9, 1e300), 'nm', '(12kT/(πΔρg))^¼: particles or aggregates larger than this settle faster than they diffuse.'),
      kpi('rms', 'Brownian displacement in one second', rms * 1e6, 'µm', 'One-dimensional root-mean-square displacement √(2Dt) for t = 1 s.'),
      kpi('Re', 'Particle Reynolds number', Re, '', 'ρ_f·v·d/μ; Stokes drag holds below about 0.1.', Re < 0.1),
    ],
    series: [
      { key: 'vd', title: 'Settling and Brownian velocities vs diameter', xlabel: 'diameter (nm)', ylabel: 'velocity (m/s)', logx: true, logy: true, lines: [{ name: 'Stokes settling', x: ds.map(x => x * 1e9), y: ds.map(x => Math.max(Math.abs(vOf(x)), 1e-30)) }, { name: 'diffusive velocity D/a', x: ds.map(x => x * 1e9), y: ds.map(x => Dof(x) / (x / 2)) }], explain: 'Settling velocity rises as d² while the diffusive velocity scale D/a falls as d⁻²; the curves cross at the Pe = 1 diameter.' },
      { key: 'prof', title: 'Equilibrium concentration profile', xlabel: 'height above the bottom (mm)', ylabel: 'c / c(bottom)', lines: [{ name: 'barometric profile', x: zs.map(x => x * 1e3), y: zs.map(x => Math.exp(-Math.min(x / lg, 700))) }], explain: 'Sedimentation–diffusion equilibrium exp(−z/l_g). A flat profile means the dispersion is stable against gravity; a sharp decay means a sediment forms.' },
    ],
    equations: [
      eqn('governing', 'Stokes–Einstein', 'D = k T / (3 π μ d)'),
      eqn('governing', 'Stokes settling', 'v = (ρ_p − ρ_f) g d² / (18 μ)'),
      eqn('formula', 'Péclet number and gravitational length', 'Pe = v a / D;   l_g = k T / [(π/6) d³ Δρ g]'),
      eqn('governing', 'Sedimentation–diffusion equilibrium', 'c(z) = c(0) exp(−z / l_g)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', Pe > 1 ? `Sedimentation dominates (Pe = ${F(Pe, 3)}): the suspension settles in about ${F(tSet / 3600, 3)} h over ${F(H * 1e3, 3)} mm. Stir during processing, reduce the particle or aggregate size below ${F(dX * 1e9, 3)} nm, or raise the liquid viscosity/density.` : `Brownian motion keeps the particles suspended (Pe = ${F(Pe, 3)}); separation by gravity is impractical, so plan for centrifugation at about ${F(Math.max(1, 0.01 / Math.max(vOf(d) / rcf, 1e-300) / 3600), 3)} g to clear 10 mm in an hour, or for membrane filtration.`, 'Péclet number of sedimentation against diffusion.'),
      rec('medium', dH > 1.5 * d ? 'Break up or prevent aggregates (see the DLVO mode) before size-dependent processing; the hydrodynamic size, not the crystallite size, controls settling and filtration.' : 'Hydrodynamic and core sizes agree: the dispersion is close to primary particles.', 'Hydrodynamic diameter from the measured diffusivity.'),
      rec('medium', `For a sedimentation or DLS measurement allow thermal equilibration: a 1 K temperature difference across the cell drives convection far faster than the ${F(av * 1e9, 3)} nm/s settling velocity.`, 'Magnitude of the settling velocity.'),
      rec('low', 'Where gravity settling is feasible, use it (or lamella settlers) in preference to centrifugation; tangential-flow filtration has a lower specific energy than high-speed centrifuges for particles below 100 nm.', 'Energy demand of solid–liquid separation.', true),
    ],
    publish: { 'prop.D_diff': D, 'colloid.settling_velocity': v, 'colloid.hydrodynamic_diameter': dH, 'colloid.peclet': Pe },
    numerics: { method: 'Closed-form Stokes–Einstein and Stokes relations', dof: 1, steps: ds.length, iterations: 1, residual: 0, notes: `Particle relaxation time ${F(tauP, 3)} s.` },
  };
}

const RUN = { dlvo: runDLVO, edl: runEDL, adsorb: runAdsorb, wetting: runWetting, transport: runTransport };

export default {
  id: 8,
  tagline: 'Colloidal stability (DLVO), electric double layer and zeta potential, adsorption isotherms and BET area, wetting and capillary collapse of nanostructures, and particle transport.',
  inputs,
  presets: [
    { name: 'ZnO nanoparticles in dilute electrolyte', note: 'ZnO sol from the hydrothermal reference chemistry: 50 nm particles, 0.03 mol/L 1:1 background electrolyte.', values: { mode: 'dlvo', dp: 5e-8, AH: 20, psi0: 35, cel: 0.03, zval: 1, phi: 1e-4, T: 25 } },
    { name: 'Citrate-stabilised gold nanoparticles, 15 nm', note: 'Large Hamaker constant and small radius: classical DLVO gives only a modest barrier, which is why citrate sols are sensitive to salt.', values: { mode: 'dlvo', dp: 1.5e-8, AH: 250, psi0: -45, cel: 0.003, zval: 1, phi: 1e-6, T: 25 } },
    { name: 'Double layer at 0.1 mol/L 1:1 electrolyte', note: 'Textbook case: Debye length 0.96 nm.', values: { mode: 'edl', cel: 0.1, zval: 1, psi0: -80, dStern: 0.3, epsS: 10, dp: 1e-7, T: 25 } },
    { name: 'BET surface area of a nanopowder (N₂, 77 K)', note: 'Demonstration isotherm corresponding to about 50 m²/g.', values: { mode: 'adsorb', sigN: 0.162, xq: 0.2, rhoP: 5606, dp: 2e-8 } },
    { name: 'ZnO nanorod array: wetting and drying collapse', note: 'As-grown hydrophilic ZnO nanorods rinsed in water.', values: { mode: 'wetting', thetaY: 50, gLV: 0.072, dRod: 8e-8, LRod: 1.5e-6, Nrod: 4e13, Emod: 1.4e11, RH: 60 } },
    { name: 'Superhydrophobic coated nanorod metasurface', note: 'Fluorosilane-coated array: Cassie–Baxter state.', values: { mode: 'wetting', thetaY: 112, gLV: 0.072, dRod: 6e-8, LRod: 1e-6, Nrod: 2e13, Emod: 1.4e11, RH: 60 } },
    { name: 'Settling of ZnO agglomerates, 1 µm', note: 'Micron-size agglomerates settle within hours; primary particles do not.', values: { mode: 'transport', dp: 1e-6, rhoP: 5606, rhoF: 997, mu: 8.9e-4, Hcol: 0.05, Dmeas: 4.9e-13, RCF: 1 } },
  ],
  run(inp) { return RUN[sel(inp, 'mode')](inp); },
};
