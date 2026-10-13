// Suite 7 — Chemical reaction, nucleation and growth.
// Seeded ZnO nanorod growth (HMTA / hydroxide bath), LaMer burst nucleation with a cohort population
// balance, Ostwald ripening and Smoluchowski aggregation, JMAK transformation kinetics, Arrhenius fitting.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt, R = N.R_GAS, LN10 = Math.LN10, T25 = 298.15, PI = Math.PI;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
const muWater = T => 2.414e-5 * Math.pow(10, 247.8 / (Math.max(T, 250) - 140));   // Pa·s (Vogel form)
const thin = (arr, n = 300) => { if (arr.length <= n) return arr.slice(); const o = []; for (let i = 0; i < n; i++) o.push(arr[Math.round(i * (arr.length - 1) / (n - 1))]); return o; };

// ---- chemistry sets: reactions generated and balanced with chem.js
const CHEMS = {
  zn_nitrate: { label: 'ZnO from zinc nitrate', solid: 'ZnO', rho: 5606, salt: 'Zn(NO3)2', by: { NaOH: 'NaNO3', NH4OH: 'NH4NO3' } },
  zn_chloride: { label: 'ZnO from zinc chloride', solid: 'ZnO', rho: 5606, salt: 'ZnCl2', by: { NaOH: 'NaCl', NH4OH: 'NH4Cl' } },
  au: { label: 'Au by citrate reduction', solid: 'Au', rho: 19300, rx: [[{ f: 'HAuCl4', phase: 'aq' }, { f: 'C6H8O7', phase: 'aq' }], [{ f: 'Au', phase: 's' }, { f: 'C5H6O5', phase: 'aq' }, { f: 'CO2', phase: 'g' }, { f: 'HCl', phase: 'aq' }]], note: 'citric acid/citrate oxidised to acetonedicarboxylic acid' },
  ag: { label: 'Ag by borohydride reduction', solid: 'Ag', rho: 10490, rx: [[{ f: 'AgNO3', phase: 'aq' }, { f: 'NaBH4', phase: 'aq' }, { f: 'H2O', phase: 'l' }], [{ f: 'Ag', phase: 's' }, { f: 'H2', phase: 'g' }, { f: 'B(OH)3', phase: 'aq' }, { f: 'NaNO3', phase: 'aq' }]], note: 'borohydride hydrolysed to boric acid' },
  tio2: { label: 'TiO₂ by hydrolysis of TiCl₄', solid: 'TiO2', rho: 3900, rx: [[{ f: 'TiCl4', phase: 'aq' }, { f: 'H2O', phase: 'l' }], [{ f: 'TiO2', phase: 's' }, { f: 'HCl', phase: 'aq' }]], note: 'forced hydrolysis' },
};
function chemEquations(key, base = 'NaOH') {
  const c = CHEMS[key], out = [];
  const rx = c.rx || [[{ f: c.salt, phase: 'aq' }, { f: base, phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: c.by[base], phase: 'aq' }, { f: 'H2O', phase: 'l' }]];
  const ie = safe(() => chem.ionicEquations(rx[0], rx[1]), null);
  if (ie) out.push(eqn('chemical', 'Balanced chemical equation', ie.molecular, c.note), eqn('ionic', 'Complete ionic equation', ie.completeIonic), eqn('ionic', 'Net ionic equation', ie.netIonic, ie.spectators.length ? `spectators: ${ie.spectators.join(', ')}` : undefined));
  else out.push(eqn('chemical', 'Reaction', `${rx[0].map(s => s.f).join(' + ')} → ${rx[1].map(s => s.f).join(' + ')}`, 'could not be balanced automatically'));
  out.push(eqn('formula', 'Solid product', chem.prettyFormula(c.solid), `M = ${F(chem.molarMass(c.solid), 5)} g/mol, ρ = ${c.rho} kg/m³, V_m = ${F(chem.molarMass(c.solid) * 1e-3 / c.rho * 1e6, 4)} cm³/mol`));
  return out;
}
const molarVolume = key => chem.molarMass(CHEMS[key].solid) * 1e-3 / CHEMS[key].rho;

// ---- synthetic Arrhenius demonstration table (forward model + seeded noise)
const DEMO_TABLE = (() => { const r = N.rng(7), rows = ['T_C,k']; for (let t = 50; t <= 90; t += 5) rows.push(`${t},${(1.5e6 * Math.exp(-65000 / (R * (t + 273.15))) * (1 + 0.04 * r.normal())).toExponential(4)}`); return rows.join('\n'); })();

const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'zno', group: 'Model', options: [
    { value: 'zno', label: 'Seeded ZnO nanorod growth (hydrothermal / CBD)' }, { value: 'lamer', label: 'LaMer burst nucleation and growth' },
    { value: 'ripen', label: 'Ostwald ripening and aggregation' }, { value: 'jmak', label: 'JMAK transformation kinetics / TTT' }, { value: 'arrfit', label: 'Arrhenius parameter estimation' }], help: 'Selects the kinetic model. All modes list the balanced reaction of the selected chemistry and the rate laws solved.' },
  { key: 'chem', label: 'Chemistry', type: 'select', value: 'zn_nitrate', group: 'Chemistry', options: Object.entries(CHEMS).map(([value, c]) => ({ value, label: c.label })), help: 'Reaction system. It sets the balanced chemical and ionic equations and the molar volume of the solid product. The nanorod mode always grows ZnO and uses the zinc salt selected here.' },
  { key: 'T', label: 'Process temperature', unit: '°C', value: 90, min: 5, max: 300, group: 'Conditions', bus: 'proc.T_C', showIf: { mode: ['zno', 'lamer', 'ripen', 'arrfit'] }, help: 'Bath or reaction temperature (also the prediction temperature of the Arrhenius fit). Rate constants follow the Arrhenius equation; diffusivities scale with T/μ.' },
  { key: 'tmin', label: 'Process time', unit: 'min', value: 180, min: 0.1, max: 6000, log: true, group: 'Conditions', showIf: { mode: ['zno', 'lamer', 'ripen', 'jmak'] }, help: 'Duration of growth, ageing or annealing.' },
  // ZnO growth
  { key: 'cZn', label: 'Zinc salt concentration', unit: 'mol/m³', value: 25, min: 0.1, max: 500, log: true, group: 'Chemistry', bus: 'proc.c_precursor', showIf: { mode: ['zno'] }, help: 'Initial Zn²⁺ concentration (25 mol/m³ = 25 mmol/L). It sets the amount of material available; the rod length depends on it only weakly because most of the zinc precipitates in the bulk and the buffered pH, not the concentration, fixes the supersaturation.' },
  { key: 'cHMTA', label: 'HMTA concentration', unit: 'mol/m³', value: 25, min: 0, max: 500, group: 'Chemistry', showIf: { mode: ['zno'] }, help: 'Hexamethylenetetramine: a weak base (pKa 4.89) that buffers the bath near pH 5–7, and a slow base source, each mole releasing four moles of ammonia on hydrolysis. Set to zero for a directly added base.' },
  { key: 'cOH', label: 'Hydroxide added directly', unit: 'mol/m³', value: 0, min: 0, max: 1000, group: 'Chemistry', bus: 'proc.c_base', showIf: { mode: ['zno'] }, help: 'Strong base (e.g. NaOH) added at the start. It raises the pH immediately, giving high supersaturation and more bulk precipitation than HMTA.' },
  { key: 'Ea', label: 'Activation energy of surface growth', unit: 'kJ/mol', value: 65, min: 10, max: 200, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'Arrhenius activation energy of the ZnO surface-incorporation step, applied to both facets and to bulk precipitation.' },
  // Defaults of k_ax, the bulk particle yield and k_h are calibrated against Mosalagae, Murape & Lepodise, Heliyon 6 (2020) e04458, Tables 1–3
  // (see the comment above runZnO and the notes in data/benchmarks/s07.js for values and residuals).
  { key: 'kax', label: 'Axial (0001) surface rate constant at 90 °C', unit: 'm/s', value: 1.3e-6, min: 1e-10, max: 1e-3, log: true, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'First-order surface rate constant of the polar top facet at the 90 °C reference. It sets how quickly rods and bulk particles follow the supersaturation (and with the facet ratio the diameter); the final length depends on it only weakly. Calibrated on one published series; seed layer and additives change it.' },
  { key: 'aniso', label: 'Facet rate ratio (0001)/(10-10)', unit: '', value: 30, min: 1, max: 500, log: true, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'Ratio of axial to lateral surface rate constants. Capping agents (citrate, PEI, chloride) change it and with it the aspect ratio. 30 gives diameters of 60–70 nm for rods 0.7 µm long, as commonly observed.' },
  { key: 'kh', label: 'HMTA hydrolysis rate constant at 90 °C', unit: '1/s', value: 2e-5, min: 1e-8, max: 1e-2, log: true, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'First-order rate constant of HMTA decomposition at 90 °C; it paces the slow growth after the initial burst. Calibrated from the growth between 2 and 3 h of one published series; activation energy 90 kJ/mol is assumed. It depends on pH in reality.' },
  { key: 'npm', label: 'Bulk ZnO particles formed per mole of zinc', unit: '1/mol', value: 4.4e12, min: 1e8, max: 1e18, log: true, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'Number of ZnO particles that nucleate in the liquid per mole of zinc charged (4.4×10¹² per mole = 1.1×10⁸ per mL at 25 mmol/L). They compete with the nanorods for the supersaturation: the rod length scales as this number to the power −⅓. Calibration parameter; dust, filtration, heating rate and additives change it.' },
  { key: 'rhoB', label: 'Bulk-particle to axial rate-constant ratio', unit: '', value: 1, min: 0.01, max: 100, log: true, group: 'Kinetics', showIf: { mode: ['zno'] }, help: 'Surface rate constant of the bulk particles (equivalent spheres) relative to the axial facet. The rod length is inversely proportional to it when both surfaces are reaction limited; it cannot be separated from the particle yield by length data alone, so it is kept at 1 unless the precipitate size is measured.' },
  { key: 'Nseed', label: 'Nanorod (seed) number density', unit: '1/m²', value: 4e13, min: 1e11, max: 1e16, log: true, group: 'Substrate', bus: 'geom.number_density', showIf: { mode: ['zno'] }, help: 'Rods per unit substrate area (4×10¹³ m⁻² = 40 per µm²), fixed by the seed layer. It bounds the diameter through the pitch.' },
  { key: 'd0', label: 'Seed diameter', unit: 'nm', value: 20, min: 2, max: 500, group: 'Substrate', showIf: { mode: ['zno'] }, help: 'Initial diameter (and height) of the seed crystallites.' },
  { key: 'Vb', label: 'Bath volume', unit: 'm³', value: 5e-5, min: 1e-7, max: 10, log: true, group: 'Reactor', bus: 'proc.volume', showIf: { mode: ['zno', 'lamer'] }, help: 'Liquid volume (5×10⁻⁵ m³ = 50 mL). A small volume per unit substrate area depletes faster.' },
  { key: 'Asub', label: 'Seeded substrate area', unit: 'cm²', value: 4, min: 0.01, max: 1e5, log: true, group: 'Reactor', showIf: { mode: ['zno'] }, help: 'Area carrying the seed layer; together with the bath volume it sets precursor depletion.' },
  { key: 'delta', label: 'Diffusion boundary-layer thickness', unit: 'µm', value: 200, min: 1, max: 5000, log: true, group: 'Reactor', showIf: { mode: ['zno'] }, help: 'Stagnant film through which species reach the substrate. About 100–300 µm in an unstirred bath; stirring makes it thinner.' },
  { key: 'DZn', label: 'Diffusivity of Zn²⁺ at 25 °C', unit: 'm²/s', value: 7.5e-10, min: 1e-11, max: 1e-8, log: true, group: 'Transport', bus: 'prop.D_diff', showIf: { mode: ['zno'] }, help: 'Aqueous diffusion coefficient at 25 °C; scaled to the bath temperature with the Stokes–Einstein T/μ rule.' },
  { key: 'DOH', label: 'Diffusivity of OH⁻ at 25 °C', unit: 'm²/s', value: 5.3e-9, min: 1e-11, max: 2e-8, log: true, group: 'Transport', showIf: { mode: ['zno'] }, help: 'Used to test whether base supply through the boundary layer could limit growth.' },
  { key: 'E', label: 'Applied electric field (towards substrate)', unit: 'V/m', value: 0, min: -5000, max: 5000, group: 'Transport', bus: 'proc.E_field', showIf: { mode: ['zno'] }, help: 'Field in the electrolyte normal to the substrate. Positive values drive cations to the substrate (Nernst–Planck drift); the sign of the effect follows the mean charge of the dissolved zinc species.' },
  { key: 'T0', label: 'Initial bath temperature', unit: '°C', value: 25, min: 5, max: 300, group: 'Conditions', showIf: { mode: ['zno'] }, help: 'Temperature at the start; the bath approaches the process temperature exponentially with the heat-up time constant.' },
  { key: 'tauH', label: 'Heat-up time constant', unit: 's', value: 600, min: 1, max: 20000, log: true, group: 'Conditions', bus: 'therm.tau_heat', showIf: { mode: ['zno'] }, help: 'First-order time constant of bath heating (from the heat-transfer suite). Slow heating delays HMTA hydrolysis and growth.' },
  // LaMer / ripening shared
  { key: 'cP', label: 'Precursor concentration', unit: 'mol/m³', value: 10, min: 1e-3, max: 2000, log: true, group: 'Chemistry', showIf: { mode: ['lamer'] }, help: 'Initial concentration of the molecular precursor that is converted to growth units (monomer).' },
  { key: 'kf', label: 'Precursor conversion rate constant', unit: '1/s', value: 0.02, min: 1e-6, max: 1e3, log: true, group: 'Kinetics', showIf: { mode: ['lamer'] }, help: 'First-order rate constant of precursor → monomer (reduction or hydrolysis). Fast conversion relative to growth gives a short burst and narrow sizes.' },
  { key: 'gam', label: 'Solid–liquid interfacial energy', unit: 'J/m²', value: 0.25, min: 0.01, max: 2, group: 'Particle', bus: 'prop.surface_energy', showIf: { mode: ['lamer', 'ripen'] }, help: 'Interfacial energy of the particle in its solution. It controls the nucleation barrier (cubic dependence) and the ripening rate.' },
  { key: 'cInf', label: 'Bulk solubility of the solid (monomer)', unit: 'mol/m³', value: 2e-3, min: 1e-9, max: 100, log: true, group: 'Particle', bus: 'thermo.solubility', showIf: { mode: ['lamer', 'ripen'] }, help: 'Equilibrium monomer concentration over a flat surface. Supersaturation is the monomer concentration divided by this value.' },
  { key: 'Dm', label: 'Monomer diffusivity', unit: 'm²/s', value: 1e-9, min: 1e-12, max: 1e-7, log: true, group: 'Particle', showIf: { mode: ['lamer', 'ripen'] }, help: 'Diffusion coefficient of the growth unit in the solvent.' },
  { key: 'ks', label: 'Surface-reaction rate constant', unit: 'm/s', value: 1e-5, min: 1e-9, max: 1, log: true, group: 'Particle', showIf: { mode: ['lamer'] }, help: 'Interface attachment coefficient in series with diffusion. Small values give reaction-limited growth, large values diffusion-limited growth (size focusing).' },
  { key: 'logA', label: 'log₁₀ nucleation prefactor', unit: 'log(1/m³/s)', value: 34, min: 22, max: 42, group: 'Particle', sens: false, showIf: { mode: ['lamer'] }, help: 'Kinetic prefactor of the classical nucleation rate; uncertain by orders of magnitude and best treated as a fitting parameter.' },
  { key: 'r0', label: 'Initial mean particle radius', unit: 'nm', value: 5, min: 0.5, max: 500, log: true, group: 'Particle', showIf: { mode: ['ripen'] }, help: 'Mean radius at the start of ageing.' },
  { key: 'N0', label: 'Initial particle number density', unit: '1/m³', value: 1e18, min: 1e10, max: 1e24, log: true, group: 'Particle', bus: 'colloid.number_density', showIf: { mode: ['ripen'] }, help: 'Particles per cubic metre of suspension; aggregation is second order in this number.' },
  { key: 'W', label: 'Stability ratio W', unit: '', value: 1e4, min: 1, max: 1e15, log: true, group: 'Particle', bus: 'colloid.W', showIf: { mode: ['ripen'] }, help: 'Fuchs stability ratio from the colloid suite: 1 is rapid (diffusion-limited) aggregation, large values mean a repulsive barrier slows aggregation by that factor.' },
  { key: 'mu', label: 'Solvent viscosity', unit: 'Pa·s', value: 3.2e-4, min: 1e-5, max: 10, log: true, group: 'Particle', showIf: { mode: ['ripen'] }, help: 'Dynamic viscosity of the suspending liquid (3.2×10⁻⁴ Pa·s for water at 90 °C, 8.9×10⁻⁴ at 25 °C).' },
  { key: 'df', label: 'Aggregate fractal dimension', unit: '', value: 1.8, min: 1.5, max: 3, group: 'Particle', showIf: { mode: ['ripen'] }, help: 'About 1.8 for diffusion-limited and 2.1 for reaction-limited cluster aggregation; 3 corresponds to coalescence into compact spheres.' },
  // JMAK
  { key: 'Tj', label: 'Annealing / transformation temperature', unit: '°C', value: 400, min: 0, max: 1500, group: 'Conditions', showIf: { mode: ['jmak'] }, help: 'Isothermal hold temperature of the transformation (e.g. crystallisation of an amorphous oxide).' },
  { key: 'nAv', label: 'Avrami exponent n', unit: '', value: 2.5, min: 0.5, max: 4, group: 'Kinetics', showIf: { mode: ['jmak'] }, help: 'Depends on nucleation mode and growth dimensionality: 4 = continuous nucleation with 3-D interface-controlled growth, 3 = site-saturated 3-D, 2.5 = continuous nucleation with 3-D diffusion-controlled growth, 1.5 = site-saturated diffusion-controlled.' },
  { key: 'k0', label: 'Rate pre-exponential factor', unit: '1/s', value: 5e8, min: 1, max: 1e20, log: true, group: 'Kinetics', showIf: { mode: ['jmak'] }, help: 'Pre-exponential factor of the JMAK rate constant k = k₀ exp(−Ea/RT).' },
  { key: 'EaJ', label: 'Activation energy of transformation', unit: 'kJ/mol', value: 150, min: 10, max: 600, group: 'Kinetics', showIf: { mode: ['jmak'] }, help: 'Effective activation energy of the overall transformation (nucleation and growth combined).' },
  { key: 'drive', label: 'Driving-force model', type: 'select', value: 'arrhenius', group: 'Kinetics', showIf: { mode: ['jmak'] }, options: [{ value: 'arrhenius', label: 'Arrhenius only (crystallisation on heating)' }, { value: 'undercool', label: 'Undercooling-limited nucleation (C-curve)' }], help: 'Arrhenius gives a monotonic TTT line. The undercooling option multiplies k by a nucleation factor that vanishes at the equilibrium temperature, giving the classic C-shaped TTT curve.' },
  { key: 'Teq', label: 'Equilibrium transformation temperature', unit: 'K', value: 1200, min: 300, max: 3500, group: 'Kinetics', showIf: { mode: ['jmak'] }, help: 'Temperature above which the product phase is not stable (used only with the undercooling option).' },
  { key: 'Bn', label: 'Nucleation barrier parameter', unit: 'K', value: 300, min: 1, max: 1e5, log: true, group: 'Kinetics', showIf: { mode: ['jmak'] }, help: 'Strength of the undercooling term exp[−B T_eq²/(T ΔT²)]; larger values push the nose of the C-curve to lower temperature.' },
  // Arrhenius fit
  { key: 'data', label: 'Rate-constant data', type: 'table', columns: ['T_C', 'k'], value: DEMO_TABLE, group: 'Data', showIf: { mode: ['arrfit'] }, help: 'Measured rate constants versus temperature (°C). Paste or upload CSV with a header line. The default is synthetic demonstration data.' },
];
const DEF = Object.fromEntries(inputs.map(i => [i.key, i]));
const num = (inp, k) => { const d = DEF[k]; let v = +inp[k]; if (!Number.isFinite(v)) v = d.value; return N.clamp(v, d.min, d.max); };
const sel = (inp, k) => (DEF[k].options.some(o => o.value === inp[k]) ? inp[k] : DEF[k].value);

// ---------------------------------------------------------------- ZnO nanorod growth
// Seeded growth from a zinc salt + HMTA (or hydroxide) bath. Physical basis of the growth law:
//  * Driving force. Facets advance at v = V_m k (c − c_eq) = V_m k c_eq (S − 1): first order in the excess of dissolved Zn(II) over the ZnO
//    solubility c_eq(pH, T) (hydrolysis constants of Baes & Mesmer, as in the speciation suite). ZnO solubility is retrograde, so heating a bath
//    whose pH is held near 5–7 by the HMTA/HMTAH⁺ (pKa 4.89) and NH₃/NH₄⁺ pairs creates the supersaturation; no base release is needed to start.
//  * Competition with homogeneous precipitation. ZnO particles in the bulk liquid grow under the same driving force. Their number is taken
//    proportional to the zinc in the bath (n_p particles per mole of zinc, i.e. a concentration-independent final particle size) and their
//    area soon exceeds that of the substrate by orders of magnitude, so they, not the rods, consume the supersaturation. Since rods and
//    particles share (c − c_eq), the rod length follows the particle radius: L − L₀ ≈ (k_ax/k_bulk)(r_p − r_p0), r_p = [3 V_m P/(4π N_p)]^⅓ with
//    P the precipitated zinc. Length therefore scales with the cube root of the precipitated FRACTION, not with the nominal concentration.
//  * Self-limitation. Each ZnO formed releases 2 H⁺; the buffer lets only about a quarter of the zinc precipitate in the initial burst (pH
//    falls to about 5, where S → 1). Further growth is paced by HMTA hydrolysis (4 NH₃ per HMTA) and ends when the zinc is exhausted.
//  * Supply to the array through a stagnant film in series with the facet reaction (rate per rod falls with rod density when diffusion limits).
// Calibration (declared): Mosalagae, Murape & Lepodise, Heliyon 6 (2020) e04458, Tables 1–3 (seeded glass, 100 mL, equimolar zinc nitrate +
// HMTA), ten lengths. Fitted by least absolute deviation of ln(length): k_ax(90 °C) = 1.3×10⁻⁶ m/s and n_p = 4.4×10¹² particles per mole of
// zinc; k_h(90 °C) = 2×10⁻⁵ 1/s was scanned over 0.6–4.5×10⁻⁵ with n_p refitted and is only weakly determined (the misfit changes by under 6 %).
// Residuals (model − measured): Table 1 (25 mM, 90 °C; 1.5, 2, 2.5, 3 h) +62, −1, 0, −8 %; Table 2 (25 mM, 2 h; 60, 70, 80 °C) 0, −22, −12 %;
// Table 3 (90 °C, 2 h; 10, 15, 30 mM) +25, +28, +62 %; mean |ln ratio| 0.19. The 1.5 h and 30 mM points are not reproduced.
// Everything else (E_a, facet ratio, seed density, film thickness, heat-up constant, equilibrium constants, E_a of HMTA hydrolysis) is not
// fitted. Independent validation (not used in the fit; data/benchmarks/s07.js): Mason et al. 25 mM, 70 °C/12 h −26 % and 85 °C/1 h −17 %;
// Ridhuan et al. 0.1 M, 80 °C, 4/8/16/24 h −9, −10, −7, −16 %; Chae et al. 0.1 M, 90 °C, 10 h −23 %. Not reproduced: Ridhuan 0.05 M, 4 h
// (321 nm measured, 680 nm modelled) and Chae 0.1 M, 60 °C, 24 h (about 400 nm measured, 655 nm modelled).
const R_P0 = 5e-9;   // radius of a bulk nucleus when the bath first becomes supersaturated, m
let BC_T = NaN, BC_OBJ = null;
function bathChem(T) {
  if (T === BC_T) return BC_OBJ;
  const vh = (lk, dH) => lk - (dH / (R * LN10)) * (1 / T - 1 / T25), p10 = x => Math.pow(10, x);
  const Kw = p10(vh(-14, 55.8e3)), KS = p10(vh(11.14, -89.3e3)), B1 = p10(vh(-8.96, 55.8e3)), B2 = p10(vh(-16.9, 111.6e3)), B3 = p10(vh(-28.4, 167.4e3)), B4 = p10(vh(-41.2, 223.2e3)), Ka = p10(vh(-9.25, 52.2e3)), KaH = p10(-4.89);
  const at = pH => { const h = Math.exp(-LN10 * pH), i = 1 / h, r1 = B1 * i, r2 = B2 * i * i, r3 = B3 * i * i * i, r4 = B4 * i * i * i * i, D = 1 + r1 + r2 + r3 + r4; return { zbar: (2 + r1 - r3 - 2 * r4) / D, ceq: Math.min(KS * h * h, 1e5) * D * 1000 }; };
  // charge balance (concentrations in mol/m³): Zn species + NH₄⁺ + HMTAH⁺ + H⁺ − OH⁻ + Na⁺ − anion = 0, monotonically decreasing in pH
  const F0 = (p, Zn, Nt, Na, An, Ht) => { const h = Math.exp(-LN10 * p), i = 1 / h, r1 = B1 * i, r2 = B2 * i * i, r3 = B3 * i * i * i, r4 = B4 * i * i * i * i; return (2 + r1 - r3 - 2 * r4) / (1 + r1 + r2 + r3 + r4) * Zn / 1000 + (Nt / 1000) * h / (h + Ka) + (Ht / 1000) * h / (h + KaH) + h - Kw / h + (Na - An) / 1000; };
  const pH = (Zn, Nt, Na, An, Ht = 0) => { let a = 0, b = 14.5, fa = F0(a, Zn, Nt, Na, An, Ht), fb = F0(b, Zn, Nt, Na, An, Ht); if (fa <= 0) return a; if (fb >= 0) return b;
    for (let i = 0; i < 80 && b - a > 1e-10; i++) { let x = i % 3 === 2 ? 0.5 * (a + b) : (a * fb - b * fa) / (fb - fa); if (!(x > a && x < b)) x = 0.5 * (a + b); const fx = F0(x, Zn, Nt, Na, An, Ht); if (fx > 0) { a = x; fa = fx; } else { b = x; fb = fx; } } return 0.5 * (a + b); };
  BC_T = T; BC_OBJ = { Kw, Ka, KaH, at, pH }; return BC_OBJ;
}
const bern = x => (Math.abs(x) < 1e-6 ? 1 + x / 2 : x > 700 ? x : x < -700 ? 0 : x / (1 - Math.exp(-x)));

function runZnO(inp) {
  const ck0 = sel(inp, 'chem'), ck = ck0.startsWith('zn_') ? ck0 : 'zn_nitrate', warnings = [], checks = [];
  if (ck !== ck0) warnings.push('The nanorod mode grows ZnO; the selected chemistry is not a zinc salt, so zinc nitrate was used.');
  const Tp = num(inp, 'T') + 273.15, T0 = num(inp, 'T0') + 273.15, tau = num(inp, 'tauH'), tEnd = num(inp, 'tmin') * 60;
  const Zn0 = num(inp, 'cZn'), H0 = num(inp, 'cHMTA'), Na = num(inp, 'cOH'), Ea = num(inp, 'Ea') * 1e3, kax90 = num(inp, 'kax'), an = num(inp, 'aniso'), kh90 = num(inp, 'kh'), npm = num(inp, 'npm'), rhoB = num(inp, 'rhoB');
  const Ns = num(inp, 'Nseed'), d0 = num(inp, 'd0'), V = num(inp, 'Vb'), As = num(inp, 'Asub') * 1e-4, del = num(inp, 'delta') * 1e-6, DZ25 = num(inp, 'DZn'), DO25 = num(inp, 'DOH'), E = num(inp, 'E');
  const Vm = molarVolume(ck), Tref = 363.15, EaH = 90e3, pitch = 1 / Math.sqrt(Ns) * 1e9, Np = npm * Zn0, rp0n = R_P0 * 1e9;   // pitch nm, Np 1/m³
  const d00 = Math.min(d0, 0.9 * pitch); if (d00 < d0) warnings.push(`Seed diameter exceeds the pitch (${F(pitch, 3)} nm) of the chosen number density; it was reduced to 90 % of the pitch.`);
  const Tt = t => Tp + (T0 - Tp) * Math.exp(-t / tau), arr = (k, EaX, T) => k * Math.exp(-(EaX / R) * (1 / T - 1 / Tref));
  const Dscale = T => (T / T25) * muWater(T25) / muWater(T);
  const volRod = (Lx, dx) => PI * dx * dx / 4 * Lx * 1e-27, vol0 = volRod(d00, d00);
  // state: [HMTA (mol/m³), bulk particle radius r_p (nm), L (nm), d (nm)]; ammonia, bulk precipitate and dissolved zinc follow from balances
  const bulkOf = y => Np * (4 * PI / 3) * ((y[1] * 1e-9) ** 3 - R_P0 ** 3) / Vm, depOf = y => Ns * As * (volRod(y[2], y[3]) - vol0) / Vm / V;   // mol per m³ of bath
  const znOf = y => Zn0 - bulkOf(y) - depOf(y);
  const local = (t, y) => {
    const T = Tt(t), bc = bathChem(T), Zn = Math.max(znOf(y), 0), NT = 4 * Math.max(H0 - y[0], 0), pH = bc.pH(Zn, NT, Na, 2 * Zn0, y[0]), s = bc.at(pH);
    const kax = arr(kax90, Ea, T), klat = kax / an * N.clamp((pitch - y[3]) / (0.05 * pitch), 0, 1), D = DZ25 * Dscale(T);
    const Pe = s.zbar * N.FARADAY * E * del / (R * T), km = (D / del) * bern(Pe);
    const dm = y[3] * 1e-9, Lm = y[2] * 1e-9, Ksurf = Ns * (PI * dm * dm / 4 * kax + PI * dm * Lm * klat);
    const drv = Zn - s.ceq, drive = Math.max(0, drv), ds = drive * km / (km + Ksurf + 1e-300);
    const rp = y[1] * 1e-9, kbulk = 1 / (1 / (rhoB * kax) + rp / D); let vp = Vm * kbulk * drv; if (drv < 0 && y[1] <= rp0n) vp = 0;   // bulk particles may redissolve, not below the nucleus size
    return { T, pH, s, kax, klat, km, Ksurf, drive, ds, Pe, kh: arr(kh90, EaH, T), vp, rp, D, bc, Zn, NT, ab: Np * 4 * PI * rp * rp };
  };
  const rhs = (t, y) => { const q = local(t, y); return [-q.kh * y[0], q.vp * 1e9, Vm * q.kax * q.ds * 1e9, 2 * Vm * q.klat * q.ds * 1e9]; };
  const segs = [0, 1e-3, 1e-2, 3e-2, 1e-1, 0.2, 0.35, 0.5, 0.75, 1].map(f => f * tEnd), per = 24, nSteps = per * (segs.length - 1), sol = { t: [0], y: [[H0, rp0n, d00, d00]], newtonIterations: 0 };
  for (let k = 1; k < segs.length; k++) { const y0 = sol.y[sol.y.length - 1].slice(); if (y0[1] < rp0n) y0[1] = rp0n; const part = N.bdf1(rhs, y0, segs[k - 1], segs[k], per, { tol: 1e-7, maxNewton: 30, nonNegative: true }); for (let i = 1; i < part.t.length; i++) { sol.t.push(part.t[i]); sol.y.push(part.y[i]); } sol.newtonIterations += part.newtonIterations; }
  const idx = thin(sol.t.map((_, i) => i), 400), ts = idx.map(i => sol.t[i]), ys = idx.map(i => sol.y[i]), qs = ts.map((t, i) => local(t, ys[i]));
  const yE = sol.y[sol.y.length - 1], qE = local(tEnd, yE), L = yE[2], d = Math.min(yE[3], pitch), AR = L / d;
  const dep = depOf(yE) * V, bulk = bulkOf(yE) * V, left = qE.Zn * V, tot = Zn0 * V, znMin = N.minOf(sol.y.map(znOf)), mbRes = Math.max(0, -znMin) / Zn0;
  const conv = 1 - qE.Zn / Zn0, Ssat = qs.map(q => N.clamp(q.Zn / Math.max(q.s.ceq, 1e-12), 1e-6, 1e12)), iS = N.argmax(Ssat), Send = Ssat[Ssat.length - 1];
  const vAx = (L - d00) / tEnd, vLat = (yE[3] - d00) / 2 / tEnd, Mz = chem.molarMass('ZnO') * 1e-3, massSub = dep * Mz, selSub = dep / Math.max(dep + bulk, 1e-300);
  const DaS = qE.Ksurf / qE.km, enh = bern(qE.Pe), heat = 22.4e3 * (dep + bulk) / tEnd;
  const cover = Math.min(1, Ns * PI * (d * 1e-9) ** 2 / 4), fPrec = bulk / tot, areaRatio = qE.ab * V / Math.max(Ns * As * (PI * (d * 1e-9) ** 2 / 4 + PI * d * L * 1e-18), 1e-300);
  const hf = Math.pow(10, -qE.pH), OHfree = qE.bc.Kw / hf * 1000, NH3 = qE.NT * qE.bc.Ka / (hf + qE.bc.Ka);
  const JZn = qE.Ksurf * qE.ds, Jbase = (DO25 * Dscale(qE.T) / del) * (OHfree + NH3 + yE[0] * qE.bc.KaH / (hf + qE.bc.KaH));
  // geometry with frames
  const rg = N.rng(11), nR = 36, side = Math.sqrt(nR / Ns) * 1e9, g = Math.round(Math.sqrt(nR)), proto = [];
  for (let i = 0; i < g; i++) for (let j = 0; j < g; j++) proto.push({ x: (i + 0.5 + 0.5 * (rg() - 0.5)) * side / g, y: (j + 0.5 + 0.5 * (rg() - 0.5)) * side / g, fl: rg.lognormal(0, 0.12), fd: rg.lognormal(0, 0.08), tilt: Math.abs(rg.normal(0, 4)), az: 360 * rg() });
  const rodsAt = (Lx, dx) => proto.map(p => ({ x: +p.x.toFixed(2), y: +p.y.toFixed(2), length: +(Lx * p.fl).toFixed(2), diameter: +Math.min(dx * p.fd, pitch).toFixed(2), tilt: +p.tilt.toFixed(2), azimuth: +p.az.toFixed(1) }));
  const fi = Array.from({ length: 24 }, (_, k) => { const tt = k / 23 * tEnd; let j = 0; while (j < sol.t.length - 1 && sol.t[j] < tt) j++; return j; }), frames = fi.map(i => ({ t: sol.t[i], rods: rodsAt(sol.y[i][2], sol.y[i][3]) }));
  const geometry = { kind: 'rods', rods: rodsAt(L, yE[3]), area: [side, side], units: 'nm', title: 'ZnO nanorod array (representative patch)', frames };
  // checks
  checks.push(chk('Zinc balance: dissolved zinc from the balance stays non-negative', mbRes < 1e-3, `dissolved Zn = initial − bulk precipitate − nanorods (closed by construction); largest overshoot ${F(mbRes, 2)} of the initial zinc`));
  checks.push(chk('Non-negative concentrations and sizes', sol.y.every(y => y.every(v => v >= 0)), 'all state variables ≥ 0 throughout'));
  checks.push(chk('Base supply through the boundary layer not limiting', Jbase >= 2 * JZn || JZn === 0, `proton-acceptor flux capacity ${F(Jbase, 3)} vs demand ${F(2 * JZn, 3)} mol/m²/s at the end`));
  checks.push(chk('Rods not coalesced into a film', yE[3] < 0.98 * pitch, `diameter ${F(yE[3], 4)} nm vs pitch ${F(pitch, 4)} nm`));
  checks.push(chk('Bath relaxed towards saturation at the end (S < 2 or zinc exhausted)', Send < 2 || qE.Zn < 0.05 * Zn0, `final S = ${F(Send, 3)} with ${F(100 * qE.Zn / Zn0, 3)} % of the zinc still dissolved; peak S = ${F(Ssat[iS], 3)} after ${F(ts[iS] / 60, 3)} min`));
  if (!(Send < 2 || qE.Zn < 0.05 * Zn0) && Ssat[iS] > 1) warnings.push(`The bath is still supersaturated at the end (S = ${F(Send, 3)}): at this temperature, concentration or duration the precipitation burst is not complete, and the length is then sensitive to the surface rate constant and its activation energy rather than to the buffer chemistry.`);
  checks.push(chk('Rod length follows the bulk-particle radius (shared driving force)', DaS > 0.3 || Ssat[iS] <= 1 || Math.abs((L - d00) / Math.max(rhoB > 0 ? (yE[1] - rp0n) / rhoB : 1, 1e-9) - 1) < 0.35, `L − L₀ = ${F(L - d00, 4)} nm vs (k_ax/k_bulk)(r_p − r_p0) = ${F((yE[1] - rp0n) / rhoB, 4)} nm (they differ when diffusion limits either surface)`));
  const inCal = H0 > 0 && Na === 0 && Zn0 >= 10 && Zn0 <= 30 && Math.abs(H0 / Zn0 - 1) < 0.05 && Tp >= 333.15 && Tp <= 363.15 && tEnd >= 5400 && tEnd <= 10800;
  if (Ssat[iS] <= 1) warnings.push('The bath never becomes supersaturated with respect to ZnO: no growth occurs. Add a base (HMTA or hydroxide), raise the temperature or the zinc concentration.');
  if (yE[3] >= 0.98 * pitch) warnings.push('Rod diameter has reached the pitch: neighbouring rods merge into a dense columnar film and lateral growth stops.');
  if (E !== 0) warnings.push('The electromigration term assumes the stated field exists in the electrolyte film next to the substrate. In a conducting bath the field is set by current density over conductivity and is largely screened in the double layer; electrode reactions (local OH⁻ generation) are not modelled.');
  if (Na > 0) warnings.push('Directly added hydroxide: the bath is strongly supersaturated from the start and nearly all zinc precipitates within seconds to minutes, mostly in the bulk. The rate constants were calibrated for HMTA baths; zincate-rich (pH > 11) growth solutions and Zn(OH)₂ intermediates are outside the calibrated range.');
  if (H0 > 0 && Na === 0 && !inCal) warnings.push(`Outside the calibration range of the rate parameters (equimolar zinc nitrate + HMTA, 10–30 mmol/L, 60–90 °C, 1.5–3 h). Independent published lengths at 25–100 mmol/L and 1–24 h are reproduced within about 25 %, with the model on the low side; seed layer, vessel and heating rate change the result between laboratories by a similar amount.`);
  if (H0 > 0 && Na === 0 && tEnd < 5400) warnings.push('Growth shorter than 1.5 h: the length depends on when the bath reaches temperature and bulk precipitation sets in, which the first-order heat-up model captures only roughly; the calibration series itself is over-predicted by 60 % at 1.5 h.');
  warnings.push('Reduced-order model: uniform bath, one size class of bulk particles whose number is proportional to the zinc charged, 1-D film transport to the array, first-order surface kinetics per facet, ideal-dilute hydrolysis equilibria (no ammine complexes, no ammonia loss to the gas phase), HMTA protonation constant taken as independent of temperature. Rate constants and the particle yield are calibration parameters.');
  const baseForm = Na > 0 && H0 === 0 ? 'NaOH' : 'NH4OH', eqs = [];
  if (H0 > 0) eqs.push(eqn('chemical', 'HMTA hydrolysis (slow base release)', safe(() => chem.equationString(['C6H12N4', 'H2O'], ['HCHO', 'NH3']), 'C₆H₁₂N₄ + 6 H₂O → 6 HCHO + 4 NH₃')), eqn('ionic', 'HMTA as a proton acceptor (pH buffer)', '(CH₂)₆N₄ + H⁺ ⇌ (CH₂)₆N₄H⁺', 'pKa = 4.89 (25 °C value, used at all temperatures)'), eqn('ionic', 'Ammonia protonation', safe(() => chem.equationString(['NH3', 'H2O'], ['NH4+', 'OH-'], null, '⇌'), 'NH₃ + H₂O ⇌ NH₄⁺ + OH⁻'), `pKa(NH₄⁺) = ${F(-Math.log10(bathChem(Tp).Ka), 4)} at ${F(Tp - 273.15, 4)} °C`));
  eqs.push(...chemEquations(ck, baseForm));
  eqs.push(eqn('ionic', 'Zinc hydrolysis (speciation)', safe(() => chem.equationString(['Zn2+', 'H2O'], ['Zn(OH)2', 'H+'], null, '⇌'), 'Zn²⁺ + 2 H₂O ⇌ Zn(OH)₂ + 2 H⁺'), 'with ZnOH⁺, Zn(OH)₃⁻ and Zn(OH)₄²⁻ included in the charge balance'),
    eqn('constitutive', 'Supersaturation with respect to ZnO', 'S = [Zn]_T / c_eq(pH, T),   c_eq = *K_s0 [H⁺]² (1 + Σ *βₙ/[H⁺]ⁿ),   ZnO(s) + 2 H⁺ ⇌ Zn²⁺ + H₂O', `log *K_s0 = ${F(Math.log10(bathChem(Tp).at(0).ceq / 1000), 4)} at ${F(Tp - 273.15, 4)} °C (11.14 at 25 °C; retrograde solubility)`),
    eqn('governing', 'Facet growth velocities (first order in relative supersaturation)', 'dL/dt = V_m k_ax (c_s − c_eq) = V_m k_ax c_eq (S_s − 1);   dd/dt = 2 V_m k_lat (c_s − c_eq)'),
    eqn('governing', 'Surface reaction and film diffusion in series (array)', 'k_m (c − c_s) = K_s (c_s − c_eq),   K_s = N [ (π d²/4) k_ax + π d L k_lat ],   k_m = (D/δ) · Pe / (1 − e^(−Pe)),   Pe = z̄ F E δ / (R T)'),
    eqn('governing', 'Competing homogeneous precipitation (bulk particles)', 'dr_p/dt = V_m (c − c_eq) / (1/k_bulk + r_p/D),   N_p = n_p [Zn]₀,   P = N_p (4π/3)(r_p³ − r_p0³)/V_m', 'equal driving force on rods and particles gives L − L₀ ≈ (k_ax/k_bulk)(r_p − r_p0)'),
    eqn('governing', 'Balances', 'd[HMTA]/dt = −k_h [HMTA];   [NH₃]_T = 4([HMTA]₀ − [HMTA]);   [Zn] = [Zn]₀ − P − (A/V) N (π d² L/4 − π d₀³/4)/V_m;   charge balance → pH'),
    eqn('constitutive', 'Arrhenius kinetics', 'k(T) = k(363 K) exp[−(E_a/R)(1/T − 1/363.15)]'));
  const good = AR > 5, ctrl = DaS < 0.1 ? 'surface-reaction limited' : DaS > 10 ? 'diffusion limited' : 'under mixed control';
  return {
    summary: `After ${F(tEnd / 60, 4)} min at ${F(Tp - 273.15, 4)} °C the model gives ZnO nanorods ${F(L, 4)} nm long and ${F(d, 4)} nm in diameter (aspect ratio ${F(AR, 3)}), a mean axial rate of ${F(vAx * 60, 3)} nm/min. Supersaturation peaks at S = ${F(Ssat[iS], 3)} after ${F(ts[iS] / 60, 3)} min and has relaxed to ${F(Send, 3)} at pH ${F(qE.pH, 3)}: ${F(100 * conv, 3)} % of the zinc has left solution, ${F(100 * fPrec, 3)} % as bulk precipitate (particles of ${F(yE[1], 3)} nm radius) and ${F(100 * dep / tot, 2)} % on the substrate. The array itself is ${ctrl}; the length is set by the competition with bulk precipitation, not by the nominal concentration.`,
    kpis: [
      kpi('L_final', 'Final nanorod length', L, 'nm', 'Length along the c-axis from integrating dL/dt = V_m k_ax (c_s − c_eq) while bulk particles compete for the same supersaturation; compare with cross-section SEM. Published seeded-growth lengths for 10–100 mmol/L baths lie between about 0.4 and 1.2 µm.'),
      kpi('d_final', 'Final nanorod diameter', d, 'nm', `Diameter from lateral (10-10) facet growth, bounded by the pitch of ${F(pitch, 4)} nm set by the seed density.`),
      kpi('aspect', 'Aspect ratio L/d', AR, '', 'Length over diameter. It follows from the facet rate ratio and the seed size; values above about 10 are typical targets for piezoelectric and sensing arrays.', good),
      kpi('v_ax', 'Mean axial growth rate', vAx * 60, 'nm/min', 'Average c-axis growth velocity over the run. The instantaneous rate is far higher during the precipitation burst and falls to a few nm/min or less once the bath is near saturation.'),
      kpi('v_lat', 'Mean lateral growth rate', vLat * 60, 'nm/min', 'Average radial growth velocity of the side facets; small compared with the axial rate when the anisotropy is high.'),
      kpi('S_peak', 'Peak supersaturation', Ssat[iS], '', `Largest ratio of dissolved zinc to the ZnO solubility at the bath pH, reached after ${F(ts[iS] / 60, 3)} min while the bath heats up and before bulk particles have grown enough to consume it.`),
      kpi('S_end', 'Final supersaturation', Send, '', 'Ratio of dissolved zinc to the ZnO solubility at the end. Values close to 1 mean that further growth is paced by base release (HMTA hydrolysis), not by the surface kinetics.', Send < 2 && Send >= 1),
      kpi('pH_end', 'Final bath pH', qE.pH, '', 'pH from the charge balance including zinc hydroxo species and the HMTA/HMTAH⁺ and ammonia/ammonium buffers at the final temperature.'),
      kpi('conv', 'Zinc conversion', conv, '', 'Fraction of the initial dissolved zinc removed from solution by growth and bulk precipitation.'),
      kpi('f_bulk', 'Zinc precipitated in the bulk', fPrec, '', 'Fraction of the zinc charged that has formed ZnO particles in the liquid. The rod length scales with its cube root; in an equimolar HMTA bath the buffer limits it to roughly a quarter until HMTA hydrolysis supplies more base.'),
      kpi('r_bulk', 'Bulk precipitate particle radius', yE[1], 'nm', `Equivalent-sphere radius of the ${F(Np, 3)} particles per m³ assumed to form in the liquid (${F(npm, 3)} per mole of zinc). Their surface is ${F(areaRatio, 3)} times that of the nanorod array, which is why they control the supersaturation.`),
      kpi('sel', 'Fraction of converted zinc on the substrate', selSub, '', 'Share of the consumed zinc that ends up as nanorods rather than bulk precipitate; a measure of precursor utilisation.', selSub > 0.3 ? true : selSub < 0.05 ? false : null),
      kpi('mass', 'ZnO deposited on the substrate', massSub * 1e6, 'mg', 'Mass of the nanorod film from the geometric volume of the rods on the seeded area.'),
      kpi('hmta', 'HMTA converted', H0 > 0 ? 1 - yE[0] / H0 : 0, '', 'Fraction of HMTA hydrolysed; the remainder is unused base reservoir (and proton buffer).'),
      kpi('DaS', 'Surface-to-diffusion ratio K_s/k_m', DaS, '', 'Second Damköhler number of the array. Below 0.1 growth is surface-reaction limited (uniform, facet controlled); above 10 it is diffusion limited (edge effects, tapered rods, rate per rod falling with rod density).'),
      kpi('enh', 'Electromigration enhancement of mass transfer', enh, '', 'Factor Pe/(1 − e^(−Pe)) multiplying the film mass-transfer coefficient due to Nernst–Planck drift; 1 means no field effect.'),
      kpi('cover', 'Substrate coverage by rod tops', cover, '', 'Projected area fraction N·πd²/4; it approaches π/4 when rods touch.'),
    ],
    series: [
      { key: 'len', title: 'Nanorod length and diameter vs time', xlabel: 'time (min)', ylabel: 'size (nm)', lines: [{ name: 'length', x: ts.map(t => t / 60), y: ys.map(y => y[2]) }, { name: 'diameter', x: ts.map(t => t / 60), y: ys.map(y => Math.min(y[3], pitch)) }, { name: 'bulk particle radius', x: ts.map(t => t / 60), y: ys.map(y => y[1]) }], explain: 'Growth curves of both rod dimensions and of the bulk precipitate particles. After a lag while the bath heats, a burst of growth ends when the bath has acidified to saturation; the slow rise afterwards is paced by HMTA hydrolysis.' },
      { key: 'conc', title: 'Bath concentrations vs time', xlabel: 'time (min)', ylabel: 'concentration (mol/m³)', lines: [{ name: 'dissolved Zn', x: ts.map(t => t / 60), y: qs.map(q => q.Zn) }, { name: 'HMTA', x: ts.map(t => t / 60), y: ys.map(y => y[0]) }, { name: 'total ammonia', x: ts.map(t => t / 60), y: qs.map(q => q.NT) }, { name: 'bulk precipitate', x: ts.map(t => t / 60), y: ys.map(y => Math.max(bulkOf(y), 0)) }, { name: 'ZnO solubility c_eq', x: ts.map(t => t / 60), y: qs.map(q => Math.min(q.s.ceq, 2 * Zn0)) }], explain: 'Depletion of zinc in the finite bath, HMTA decomposition, the zinc lost to bulk precipitate and the solubility at the bath pH (clipped at twice the initial zinc). Dissolved zinc follows the solubility once the burst is over.' },
      { key: 'S', title: 'Supersaturation and pH vs time', xlabel: 'time (min)', ylabel: 'log₁₀ S  /  pH', lines: [{ name: 'log₁₀ S', x: ts.map(t => t / 60), y: Ssat.map(v => Math.log10(v)) }, { name: 'pH', x: ts.map(t => t / 60), y: qs.map(q => q.pH) }, { name: 'temperature / 100 (°C)', x: ts.map(t => t / 60), y: qs.map(q => (q.T - 273.15) / 100) }], explain: 'Driving force for ZnO formation (log scale, positive when supersaturated) with the bath pH and the heating curve. Heating raises S because ZnO is less soluble when hot; precipitation releases acid and brings S back to 1.' },
      { key: 'rate', title: 'Facet growth rates vs time', xlabel: 'time (min)', ylabel: 'growth rate (nm/min)', lines: [{ name: 'axial (0001)', x: ts.map(t => t / 60), y: qs.map(q => Vm * q.kax * q.ds * 6e10) }, { name: 'lateral (10-10)', x: ts.map(t => t / 60), y: qs.map(q => Vm * q.klat * q.ds * 6e10) }], explain: 'Instantaneous facet velocities from the series surface-reaction/diffusion model; their ratio sets the evolving aspect ratio. The peak marks the precipitation burst.' },
    ],
    tables: [{ key: 'bal', title: 'Zinc balance at the end of the run', columns: ['pool', 'amount (µmol)', 'share (%)'], rows: [['dissolved', +(left * 1e6).toPrecision(4), +(100 * left / tot).toPrecision(4)], ['nanorods on substrate', +(dep * 1e6).toPrecision(4), +(100 * dep / tot).toPrecision(4)], ['bulk precipitate', +(bulk * 1e6).toPrecision(4), +(100 * bulk / tot).toPrecision(4)]], explain: 'Where the zinc initially in the bath ends up. A large bulk-precipitate share means poor precursor utilisation.' }],
    equations: eqs, geometry, checks, warnings,
    recommendations: [
      rec('high', Ssat[iS] <= 1 ? 'No supersaturation is reached: raise the temperature, add HMTA or another base, or increase the zinc concentration.' : `Raising the precursor concentration will not lengthen the rods in proportion: the length follows the cube root of the precipitated fraction (now ${F(100 * fPrec, 3)} % of the zinc). For longer rods extend the time (HMTA has released ${F(100 * (H0 > 0 ? 1 - yE[0] / H0 : 0), 3)} % of its base), refresh the bath for a second growth step, or suppress bulk nucleation (filtered solutions, lower heating rate, a complexing agent such as PEI or ammonia) so that fewer, larger bulk particles compete.`, 'Shared driving force between nanorods and bulk precipitate.'),
      rec('medium', AR < 10 ? `Aspect ratio is ${F(AR, 3)}. Increase the facet rate ratio with a side-facet capping agent (e.g. polyethyleneimine) or use smaller seeds; reducing diameter growth also delays coalescence at ${F(pitch, 3)} nm pitch.` : `Aspect ratio ${F(AR, 3)} is high; check the collapse risk during drying in the surface suite before rinsing in water.`, 'Facet rate ratio and seed geometry.'),
      rec('medium', DaS > 1 ? `Transport is significant (K_s/k_m = ${F(DaS, 3)}): stir or use a thinner boundary layer to avoid tapered rods and edge-to-centre non-uniformity.` : `The array is surface-reaction limited (K_s/k_m = ${F(DaS, 3)}), so bath uniformity matters more than stirring. The final length is only weakly temperature dependent above about 70 °C, but the time at which the precipitation burst occurs shifts with the heating rate: record the bath temperature history.`, 'Series resistance of surface reaction and film diffusion.'),
      rec('medium', `Only ${F(100 * dep / tot, 2)} % of the zinc in the bath ends up on the substrate. Raise the substrate area per bath volume (now ${F(As / V, 3)} m²/m³), stack substrates, or reuse the spent bath for a second growth step; recover the bulk ZnO precipitate by filtration rather than discharging zinc to drain.`, 'Zinc balance of the finite bath.', true),
      rec('low', 'Formaldehyde and ammonia are released by HMTA hydrolysis: run in a closed vessel or fume hood and neutralise the spent bath before disposal.', 'Stoichiometry of HMTA hydrolysis.', true),
    ],
    publish: { 'growth.rate': Math.max(vAx, 0) * 1e-9, 'growth.length': L * 1e-9, 'growth.diameter': d * 1e-9, 'geom.length': L * 1e-9, 'geom.diameter': d * 1e-9, 'geom.aspect': AR, 'geom.number_density': Ns, 'rxn.conversion': conv, 'rxn.yield_kg': massSub, 'rxn.heat': heat, 'proc.supersaturation': Ssat[iS], 'proc.pH': qE.pH, 'proc.time': tEnd, 'material.formula': 'ZnO' },
    numerics: { method: 'Backward Euler + Newton (finite-difference Jacobian) on four coupled ODEs (HMTA, bulk particle radius, rod length, rod diameter), nine graded time segments; dissolved zinc and ammonia from the balances; bracketed regula falsi on the charge balance for pH at every evaluation', dof: 4, steps: nSteps, iterations: sol.newtonIterations, residual: mbRes, notes: 'First-order implicit scheme; halving the step changes the final length by well under 1 %.' },
  };
}

// ---------------------------------------------------------------- LaMer burst nucleation (cohort population balance)
function runLaMer(inp) {
  const ck = sel(inp, 'chem'), T = num(inp, 'T') + 273.15, tEnd = num(inp, 'tmin') * 60, P0 = num(inp, 'cP'), kf = num(inp, 'kf'), g = num(inp, 'gam'), Cinf = num(inp, 'cInf'), D = num(inp, 'Dm'), ks = num(inp, 'ks'), A = Math.pow(10, num(inp, 'logA')), Vb = num(inp, 'Vb');
  const Vm = molarVolume(ck), v = Vm / N.N_A, kT = N.K_B * T, lcap = 2 * g * Vm / (R * T), B = 16 * PI * g ** 3 * v * v / (3 * kT ** 3), rmin = Math.cbrt(3 * 2 * v / (4 * PI)), aS = D / ks, warnings = [], checks = [];
  let t = 0, P = P0, C = 0, Nn = [], rr = [], dt = Math.min(1 / kf, tEnd) * 1e-4, steps = 0, rejected = 0, clampedNuclei = false, evals = 0, work = 0;
  const hist = { t: [], S: [], J: [], N: [], r: [], P: [], sol: [] }, maxSteps = 9000;
  const solidOf = (nArr, rArr) => { let s = 0; for (let i = 0; i < nArr.length; i++) s += nArr[i] * rArr[i] ** 3; return s * 4 * PI / (3 * Vm); };
  const Jof = c => { const S = c / Cinf; if (S <= 1) return 0; const l = Math.log(S); return A * Math.exp(-Math.min(B / (l * l), 700)); };
  const record = (J) => { let n = 0, m = 0; for (let i = 0; i < Nn.length; i++) { n += Nn[i]; m += Nn[i] * rr[i]; } hist.t.push(t); hist.S.push(C / Cinf); hist.J.push(J); hist.N.push(n); hist.r.push(n > 0 ? m / n : 0); hist.P.push(P); hist.sol.push(solidOf(Nn, rr)); };
  record(0);
  let Cr = [], rn = [];
  while (t < tEnd * (1 - 1e-12) && steps < maxSteps && work < 4.5e6) {
    dt = Math.min(dt, tEnd - t);
    const S = C / Cinf, lnS = S > 1 ? Math.log(S) : 0, m = rr.length; let J = Jof(C), ntot = 0; for (let i = 0; i < m; i++) ntot += Nn[i];
    if (J * dt < Math.max(1e-7 * ntot, 1)) J = 0;   // negligible nucleation is not tracked as a cohort
    let rs = lnS > 0 ? 1.1 * lcap / lnS : 0, clampNow = false; if (J > 0 && rs < rmin) { rs = rmin; clampNow = true; }
    const nMol = (4 / 3) * PI * rs ** 3 / Vm, dCn = J * dt * nMol, Pn = P * Math.exp(-kf * dt), VDt = Vm * D * dt;
    if (Cr.length !== m) { Cr = new Array(m); rn = new Array(m); }
    let rbar = 0, nt = 0; for (let i = 0; i < m; i++) { Cr[i] = Cinf * Math.exp(Math.min(lcap / rr[i], 60)); rbar += Nn[i] * rr[i]; nt += Nn[i]; } rbar = nt > 0 ? rbar / nt : 0;
    // implicit monomer balance: g(c) = total − precursor − solid(c) − c = 0, growth of each cohort integrated exactly at fixed c
    const g = c => { let sld = 0; for (let i = 0; i < m; i++) { const q = rr[i] * rr[i] / 2 + aS * rr[i] + VDt * (c - Cr[i]), x0 = q <= 0 ? 0 : -aS + Math.sqrt(aS * aS + 2 * q), x = x0 < 0.8 * rmin ? 0 : x0; rn[i] = x; sld += Nn[i] * x * x * x; } evals++; work += m + 4; return P0 - Pn - sld * 4 * PI / (3 * Vm) - dCn - c; };
    const hi = P0 - Pn; let Cn;
    let lo = Math.max(0, 0.7 * C), up = Math.min(hi, 1.4 * C + 1e-3 * Cinf); if (!(up > lo) || g(lo) <= 0 || g(up) >= 0) { lo = 0; up = hi; }
    if (lo === 0 && g(0) <= 0) Cn = 0; else if (up === hi && g(hi) >= 0) Cn = hi; else Cn = safe(() => N.brent(g, lo, up, 1e-9 * Math.max(C, Cinf), 60), N.bisect(c => -g(c), lo, up, 1e-9 * Math.max(C, Cinf)));
    Cn = Math.max(0, Cn + g(Cn));   // exact closure of the mass balance with the radii just computed
    let bad = dCn > 0.05 * Math.max(C, 1e-300) && dCn > 1e-9 * P0;
    if (!bad && m) { let stot = 0; for (let i = 0; i < m; i++) stot += Nn[i] * rr[i] ** 3; for (let i = 0; i < m; i++) if (Nn[i] * rr[i] ** 3 > 1e-5 * stot && rr[i] >= 0.75 * rbar && Math.abs(rn[i] - rr[i]) > 0.2 * rr[i] + 0.03 * rbar) { bad = true; break; } }
    if (!bad) { const Jn = Jof(Cn), jm = Math.max(J, Jn); if (jm * dt * nMol > 1e-7 * P0 && Math.abs(Math.log((Jn + 1e-300) / (J + 1e-300))) > 0.6) bad = true; if (Math.abs(Cn - C) > 0.3 * Math.max(C, Cinf) && t > 0) bad = true; }
    if (bad && dt > 1e-13 * tEnd) { dt *= 0.5; rejected++; steps++; continue; }
    if (clampNow) clampedNuclei = true;
    P = Pn; C = Cn; t += dt;
    const n2 = [], r2 = []; for (let i = 0; i < m; i++) if (rn[i] > 0 && Nn[i] > 0) { n2.push(Nn[i]); r2.push(rn[i]); }
    if (J * dt > 0) { n2.push(J * dt); r2.push(rs); }
    Nn = n2; rr = r2;
    if (Nn.length > 120) { const ord = rr.map((_, i) => i).sort((i, j) => rr[i] - rr[j]); Nn = ord.map(i => Nn[i]); rr = ord.map(i => rr[i]); const a = [], b = []; for (let i = 0; i + 1 < Nn.length; i += 2) { const nn = Nn[i] + Nn[i + 1]; a.push(nn); b.push(Math.cbrt((Nn[i] * rr[i] ** 3 + Nn[i + 1] * rr[i + 1] ** 3) / nn)); } if (Nn.length % 2) { a.push(Nn[Nn.length - 1]); b.push(rr[rr.length - 1]); } Nn = a; rr = b; }
    record(J); steps++;
    dt = Math.min(dt * 1.3, tEnd / 150);
  }
  const done = t >= tEnd * (1 - 1e-9), Ntot = N.sum(Nn), solid = solidOf(Nn, rr), mb = Math.abs(P0 - P - C - solid) / P0;
  let rMean = 0, r2m = 0, r3m = 0; for (let i = 0; i < Nn.length; i++) { rMean += Nn[i] * rr[i]; r2m += Nn[i] * rr[i] ** 2; r3m += Nn[i] * rr[i] ** 3; }
  rMean = Ntot > 0 ? rMean / Ntot : 0; const sd = Ntot > 0 ? Math.sqrt(Math.max(0, r2m / Ntot - rMean * rMean)) : 0, cv = rMean > 0 ? sd / rMean : 0, rVol = Ntot > 0 ? Math.cbrt(r3m / Ntot) : 0;
  const iJ = N.argmax(hist.J), Jmax = hist.J[iJ], tInd = Jmax > 0 ? hist.t[iJ] : 0, Smax = N.maxOf(hist.S);
  let tA = tInd, tB = tInd; if (Jmax > 0) { for (let i = iJ; i >= 0 && hist.J[i] > 0.01 * Jmax; i--) tA = hist.t[i]; for (let i = iJ; i < hist.J.length && hist.J[i] > 0.01 * Jmax; i++) tB = hist.t[i]; }
  const conv = solid / P0, Ms = chem.molarMass(CHEMS[ck].solid) * 1e-3;
  // size distribution
  const dia = rr.map(x => 2e9 * x); let psd = { x: [0, 1], y: [0, 0] };
  if (Ntot > 0) { const lo = N.minOf(dia), hi = N.maxOf(dia), span = Math.max(hi - lo, 0.02 * hi + 1e-6), nb = 30, y = new Array(nb).fill(0), w = (span * 1.2) / nb, x0 = lo - 0.1 * span; for (let i = 0; i < dia.length; i++) y[N.clamp(Math.floor((dia[i] - x0) / w), 0, nb - 1)] += Nn[i] / Ntot; psd = { x: y.map((_, i) => x0 + (i + 0.5) * w), y }; }
  const H = k => thin(hist[k], 300), ht = H('t'), pos = a => a.map(vv => Math.max(vv, 1e-30));
  checks.push(chk('Monomer mass balance (precursor + monomer + solid; conservative closure each step)', mb < 1e-9, `relative residual ${F(mb, 2)}`));
  checks.push(chk('Integration reached the end time', done, `${steps} steps (${rejected} rejected), final time ${F(t, 4)} s of ${F(tEnd, 4)} s`));
  checks.push(chk('Critical nucleus at least two formula units', !clampedNuclei, clampedNuclei ? 'nuclei were clamped to a dimer: classical nucleation theory is at its limit' : 'r* above the dimer size throughout'));
  checks.push(chk('Nucleation and growth separated in time (burst shorter than 20 % of the run)', (tB - tA) < 0.2 * tEnd, `nucleation window ${F(tB - tA, 3)} s`));
  if (!done) warnings.push(`The adaptive integrator used its computational budget at t = ${F(t, 3)} s, before the end time (stiff ripening stage with many small particles); results refer to that time. Shorten the process time or use the ripening mode for the ageing stage.`);
  if (clampedNuclei) warnings.push('At the peak supersaturation the classical critical nucleus is smaller than a dimer; nuclei were introduced as dimers. Treat the number density as an order-of-magnitude estimate.');
  if (Ntot === 0) warnings.push('No nucleation occurred: the monomer concentration never produced a measurable nucleation rate. Raise the precursor concentration or conversion rate, or lower the interfacial energy/solubility.');
  warnings.push('Aggregative growth and coalescence, which dominate many metal-nanoparticle syntheses, are not included; particles grow by monomer addition only. The classical nucleation prefactor and interfacial energy are effective fitting parameters.');
  const regime = aS > 10 * rMean ? 'surface-reaction limited' : aS < 0.1 * rMean ? 'diffusion limited' : 'mixed diffusion/reaction control';
  return {
    summary: `Burst nucleation peaks at ${F(tInd, 3)} s with a maximum supersaturation of ${F(Smax, 3)} and lasts about ${F(tB - tA, 3)} s. After ${F(tEnd / 60, 4)} min the model predicts ${F(Ntot, 3)} particles per m³ with a mean diameter of ${F(2e9 * rMean, 4)} nm and a relative standard deviation of ${F(100 * cv, 3)} %; ${F(100 * conv, 4)} % of the precursor is in the solid. Growth is ${regime}.`,
    kpis: [
      kpi('d_mean', 'Number-mean particle diameter', 2e9 * rMean, 'nm', 'Mean of the number-weighted size distribution from the cohort population balance; compare with TEM counting statistics.'),
      kpi('cv', 'Relative standard deviation of size', cv, '', `Standard deviation over mean radius (σ = ${F(2e9 * sd, 3)} nm in diameter). Below about 0.1 is usually called monodisperse.`, cv < 0.1 ? true : cv > 0.25 ? false : null),
      kpi('pdi', 'Polydispersity index (σ/mean)²', cv * cv, '', 'Square of the relative standard deviation, comparable in definition with the cumulant polydispersity index reported by dynamic light scattering.'),
      kpi('N', 'Particle number density', Ntot, '1/m³', 'Total number of particles per cubic metre surviving at the end; fixed essentially during the nucleation burst.'),
      kpi('t_ind', 'Induction time (peak nucleation)', tInd, 's', 'Time at which the nucleation rate is highest, i.e. when the monomer concentration peaks in the LaMer diagram.'),
      kpi('burst', 'Nucleation window', tB - tA, 's', 'Duration over which the nucleation rate exceeds 1 % of its maximum. A short window relative to the growth time gives narrow distributions.'),
      kpi('S_max', 'Maximum supersaturation', Smax, '', 'Peak monomer concentration divided by the bulk solubility; the critical supersaturation of the LaMer picture.'),
      kpi('J_max', 'Peak nucleation rate', Jmax, '1/m³/s', 'Maximum classical nucleation rate J = A exp(−B/ln²S) reached during the burst.'),
      kpi('conv', 'Precursor converted to solid', conv, '', 'Fraction of the precursor present as particles at the end; the rest is unreacted precursor and dissolved monomer.', conv > 0.9),
      kpi('d_vol', 'Volume-equivalent mean diameter', 2e9 * rVol, 'nm', 'Diameter of the sphere of mean particle volume; follows directly from yield and number density.'),
      kpi('yield', 'Solid yield in the batch', solid * Ms * Vb * 1e6, 'mg', 'Mass of particles in the stated bath volume from the converted amount and the molar mass of the product.'),
    ],
    series: [
      { key: 'lamer', title: 'LaMer diagram: supersaturation vs time', xlabel: 'time (s)', ylabel: 'supersaturation S', logx: true, logy: true, lines: [{ name: 'S = C/C∞', x: ht.map(x => Math.max(x, 1e-9)), y: pos(H('S')) }, { name: 'saturation', x: [Math.max(ht[1] || 1e-9, 1e-9), tEnd], y: [1, 1] }], explain: 'Monomer concentration rises as the precursor converts, overshoots the critical level where nucleation bursts, then relaxes towards saturation as particles grow.' },
      { key: 'J', title: 'Nucleation rate and number density', xlabel: 'time (s)', ylabel: 'J (1/m³/s)  /  N (1/m³)', logx: true, logy: true, lines: [{ name: 'nucleation rate J', x: ht.map(x => Math.max(x, 1e-9)), y: pos(H('J')).map(x => Math.max(x, 1)) }, { name: 'number density N', x: ht.map(x => Math.max(x, 1e-9)), y: pos(H('N')).map(x => Math.max(x, 1)) }], explain: 'The nucleation rate is a sharp pulse; the number density steps up during it and then stays constant or falls slightly as the smallest particles redissolve (values floored at 1 for the log axis).' },
      { key: 'r', title: 'Mean particle diameter vs time', xlabel: 'time (s)', ylabel: 'diameter (nm)', logx: true, lines: [{ name: 'number-mean diameter', x: ht.map(x => Math.max(x, 1e-9)), y: H('r').map(x => 2e9 * x) }], explain: 'Growth by monomer addition after the burst: diffusion-limited growth follows r ∝ t^½ while monomer lasts, then levels off.' },
      { key: 'psd', title: 'Final particle-size distribution', xlabel: 'diameter (nm)', ylabel: 'number fraction per bin', lines: [{ name: 'number distribution', x: psd.x, y: psd.y }], explain: 'Histogram of particle diameters from the nucleation cohorts. Particles born early are largest; diffusion-limited growth narrows the relative width.' },
      { key: 'mass', title: 'Distribution of material vs time', xlabel: 'time (s)', ylabel: 'concentration (mol/m³)', logx: true, lines: [{ name: 'precursor', x: ht.map(x => Math.max(x, 1e-9)), y: H('P') }, { name: 'solid', x: ht.map(x => Math.max(x, 1e-9)), y: H('sol') }], explain: 'Precursor consumed by first-order conversion and solid formed; the difference is dissolved monomer.' },
    ],
    equations: [...chemEquations(ck),
      eqn('governing', 'Precursor and monomer balance', 'dP/dt = −k_f P;   dC/dt = k_f P − J n* − Σᵢ Nᵢ 4π rᵢ² (drᵢ/dt) / V_m'),
      eqn('governing', 'Classical nucleation rate', 'J = A exp[−16π γ³ v² / (3 (kT)³ ln²S)],   S = C / C_∞,   r* = 2 γ V_m / (R T ln S)'),
      eqn('governing', 'Growth law (diffusion and surface reaction in series)', 'dr/dt = V_m D (C − C_r) / (r + D/k_s),   C_r = C_∞ exp[2 γ V_m / (r R T)]'),
      eqn('governing', 'Population balance (method of characteristics)', '∂n/∂t + ∂(G n)/∂r = J δ(r − r*);   each nucleation cohort (Nᵢ, rᵢ) is advected with G(rᵢ)')],
    checks, warnings,
    recommendations: [
      rec('high', cv > 0.15 ? `Size spread is ${F(100 * cv, 3)} %. Shorten the nucleation window (now ${F(tB - tA, 3)} s): inject the reductant/base rapidly at temperature, raise the conversion rate constant, or use seeded growth so that nucleation and growth are fully separated.` : `Size spread is ${F(100 * cv, 3)} %: nucleation and growth are well separated. Keep mixing times below the nucleation window of ${F(tB - tA, 3)} s when scaling up, otherwise the burst is smeared out.`, 'Width of the nucleation pulse relative to growth.'),
      rec('medium', `Particle size scales as (converted amount / number density)^⅓. To change the diameter from ${F(2e9 * rMean, 3)} nm, change the number of nuclei (stronger or weaker burst) rather than the reaction time.`, 'Mass balance between yield and number density.'),
      rec('medium', regime === 'diffusion limited' ? 'Growth is diffusion limited, which focuses the size distribution; keep monomer supply up (slow continued precursor feed) to stay in the focusing regime.' : `Growth is ${regime}; reaction-limited growth does not focus sizes. A higher temperature or a less strongly bound ligand moves growth towards diffusion control.`, 'Ratio D/k_s to particle radius.'),
      rec('low', conv < 0.95 ? `Only ${F(100 * conv, 3)} % of the precursor is converted; extend the time or temperature before work-up so that unreacted metal precursor does not end up in the waste stream.` : `Conversion is ${F(100 * conv, 4)} %: stop the reaction now to avoid ripening and save heating energy; recover the metal from washings.`, 'Precursor conversion at the end of the run.', true),
    ],
    publish: { 'nucl.rate': Jmax, 'geom.diameter': 2 * rMean, 'colloid.number_density': Ntot, 'rxn.conversion': conv, 'growth.rate': tEnd > 0 ? rMean / tEnd : 0, 'rxn.yield_kg': solid * Ms * Vb, 'proc.supersaturation': Math.max(Smax, 1e-30), 'colloid.cv': cv, 'material.formula': CHEMS[ck].solid },
    numerics: { method: 'Cohort (method-of-characteristics) population balance; implicit monomer balance solved by Brent iteration with exact per-cohort growth integration; adaptive step limited by relative change of radius, monomer and nucleation rate', dof: Nn.length + 2, steps, iterations: evals, residual: mb, notes: `${rejected} rejected steps; size-adjacent cohorts merged pairwise (number and volume conserving) above 120; cohorts below dimer size dissolved; nucleation below 10⁻⁷ of the population per step neglected; nuclei introduced at 1.1 r*; Gibbs–Thomson term evaluated at the start of each step; cohorts smaller than 0.75 of the mean radius are not used for step control.` },
  };
}

// ---------------------------------------------------------------- Ostwald ripening + Smoluchowski aggregation
function runRipen(inp) {
  const ck = sel(inp, 'chem'), T = num(inp, 'T') + 273.15, tEnd = num(inp, 'tmin') * 60, g = num(inp, 'gam'), Cinf = num(inp, 'cInf'), D = num(inp, 'Dm'), r0 = num(inp, 'r0') * 1e-9, N0 = num(inp, 'N0'), W = num(inp, 'W'), mu = num(inp, 'mu'), df = num(inp, 'df');
  const Vm = molarVolume(ck), warnings = [], checks = [];
  const K = 8 * g * D * Cinf * Vm * Vm / (9 * R * T), rOf = t => Math.cbrt(r0 ** 3 + K * t), rT = rOf(tEnd), tDoubleR = 7 * r0 ** 3 / K;
  const kS = 4 * N.K_B * T / (3 * mu * W), tau = 1 / (kS * N0), NofT = t => N0 / (1 + t / tau), Nend = NofT(tEnd), nAgg = N0 / Nend, rAgg = r0 * Math.pow(nAgg, 1 / df), tDoubleA = tau * (Math.pow(2, df) - 1);
  const phi = N0 * (4 / 3) * PI * r0 ** 3, lcap = 2 * g * Vm / (R * T);
  // discrete Smoluchowski with the Brownian fractal kernel, dimensionless time θ = t/τ
  const nc = 60, thEnd = Math.min(tEnd / tau, 6), inv = 1 / df, pw = Array.from({ length: nc + 1 }, (_, i) => [Math.pow(i || 1, inv), Math.pow(i || 1, -inv)]);
  const ker = (i, j) => 0.25 * (pw[i][0] + pw[j][0]) * (pw[i][1] + pw[j][1]);
  const rhs = (_, y) => { const dy = new Array(nc).fill(0); for (let i = 1; i <= nc; i++) { const yi = y[i - 1]; if (yi <= 0) continue; for (let j = i; j <= nc; j++) { const yj = y[j - 1]; if (yj <= 0) continue; const rate = 2 * ker(i, j) * yi * yj * (i === j ? 0.5 : 1); dy[i - 1] -= rate; dy[j - 1] -= rate; if (i + j <= nc) dy[i + j - 1] += rate; } } return dy; };
  const y0 = new Array(nc).fill(0); y0[0] = 1;
  const sm = thEnd > 1e-9 ? N.rk45(rhs, y0, 0, thEnd, { rtol: 1e-5, atol: 1e-9, maxSteps: 3000 }) : { t: [0], y: [y0], steps: 0 };
  const yE = sm.y[sm.y.length - 1], massIn = yE.reduce((s, v, i) => s + v * (i + 1), 0), numIn = N.sum(yE);
  const ts = N.logspace(Math.max(tEnd * 1e-4, 1e-3), tEnd, 120);
  // LSW stationary distribution
  const us = N.linspace(0.01, 1.49, 150), gl = us.map(u => u * u * Math.pow(3 / (3 + u), 7 / 3) * Math.pow(1.5 / (1.5 - u), 11 / 3) * Math.exp(-u / (1.5 - u))), nrm = N.trapz(us, gl);
  const dominant = tDoubleA < tDoubleR ? 'aggregation' : 'Ostwald ripening';
  checks.push(chk('Cluster mass retained in the 60-class Smoluchowski window', massIn > 0.97, `${F(100 * massIn, 4)} % of mass in classes 1–60 at θ = ${F(thEnd, 3)}`));
  checks.push(chk('Numerical cluster count agrees with the constant-kernel solution within 25 %', Math.abs(numIn - 1 / (1 + thEnd)) < 0.25 / (1 + thEnd) + (1 - massIn), `numerical N/N₀ = ${F(numIn, 4)}, analytical ${F(1 / (1 + thEnd), 4)}`));
  checks.push(chk('Dilute suspension (volume fraction < 5 %) for LSW and Smoluchowski theory', phi < 0.05, `solid volume fraction ${F(phi, 3)}`));
  checks.push(chk('Particles larger than the capillary length (linearised Gibbs–Thomson valid)', r0 > lcap, `r₀ = ${F(r0 * 1e9, 3)} nm, capillary length ${F(lcap * 1e9, 3)} nm`));
  if (phi >= 0.05) warnings.push('Volume fraction above 5 %: LSW underestimates ripening (diffusion fields overlap) and binary-collision aggregation theory is inaccurate.');
  if (r0 <= lcap) warnings.push('Initial radius is below the capillary length: the linearised Gibbs–Thomson relation used by LSW theory is not valid and the cube law underestimates early coarsening.');
  if (tEnd / tau > 6) warnings.push('The cluster-size distribution is shown at t = 6τ, where the 60-class window still holds the mass; the analytical number decay is used for the end time.');
  warnings.push('LSW theory gives the asymptotic diffusion-limited coarsening rate of a dilute system; aggregation uses the Brownian kernel with a size-independent stability ratio and no break-up or sedimentation.');
  return {
    summary: `Over ${F(tEnd / 60, 4)} min at ${F(T - 273.15, 4)} °C, LSW ripening takes the mean radius from ${F(r0 * 1e9, 3)} to ${F(rT * 1e9, 4)} nm (rate constant ${F(K * 1e27, 3)} nm³/s). Brownian aggregation with W = ${F(W, 3)} has a half-time of ${F(tau, 3)} s, leaving ${F(100 * Nend / N0, 3)} % of the initial particle number and a mean aggregate radius of ${F(rAgg * 1e9, 4)} nm. The faster coarsening route is ${dominant}.`,
    kpis: [
      kpi('K_lsw', 'LSW ripening rate constant', K * 1e27, 'nm³/s', 'K = 8γDC∞V_m²/(9RT): slope of mean radius cubed versus time for diffusion-limited Ostwald ripening.'),
      kpi('r_end', 'Mean radius after ripening', rT * 1e9, 'nm', `Cube law r³ = r₀³ + K t at the end time; a growth factor of ${F(rT / r0, 3)} relative to the starting radius.`),
      kpi('t2_rip', 'Time to double the radius by ripening', tDoubleR, 's', 'Time 7r₀³/K for the mean radius to double by dissolution of small particles and growth of large ones.'),
      kpi('k_smol', 'Aggregation rate constant', kS, 'm³/s', 'k = 4kT/(3μW), the Smoluchowski rate constant for the decay of total particle number, reduced by the stability ratio.'),
      kpi('t_half', 'Aggregation half-time', tau, 's', 'Time 1/(kN₀) for the particle number to halve. Seconds or less means an unstable sol; days or more means practically stable.', tau > 10 * tEnd ? true : tau < tEnd ? false : null),
      kpi('N_rel', 'Remaining particle number N/N₀', Nend / N0, '', 'Smoluchowski solution N = N₀/(1 + t/τ) at the end time.'),
      kpi('r_agg', 'Mean aggregate radius', rAgg * 1e9, 'nm', `Radius of gyration scale r₀ (N₀/N)^(1/d_f) for fractal aggregates of ${F(nAgg, 3)} primary particles on average.`),
      kpi('t2_agg', 'Time to double the size by aggregation', tDoubleA, 's', 'Time for the mean aggregate radius to reach twice the primary radius, τ(2^d_f − 1).'),
      kpi('dominant', 'Dominant coarsening mechanism', dominant, '', 'Whichever of ripening and aggregation doubles the size sooner under the stated conditions.'),
      kpi('phi', 'Solid volume fraction', phi, '', 'N₀·(4/3)πr₀³; both theories assume a dilute suspension.'),
    ],
    series: [
      { key: 'rt', title: 'Size evolution: ripening and aggregation', xlabel: 'time (s)', ylabel: 'radius (nm)', logx: true, logy: true, lines: [{ name: 'LSW mean radius', x: ts, y: ts.map(x => rOf(x) * 1e9) }, { name: 'mean aggregate radius', x: ts, y: ts.map(x => r0 * Math.pow(N0 / NofT(x), 1 / df) * 1e9) }], explain: 'LSW ripening follows t^⅓ at long times; aggregation follows t^(1/d_f). The upper curve at a given time is the mechanism that controls the measured hydrodynamic size.' },
      { key: 'nt', title: 'Particle number decay', xlabel: 'time (s)', ylabel: 'N/N₀', logx: true, logy: true, lines: [{ name: 'aggregation, N₀/(1 + t/τ)', x: ts, y: ts.map(x => NofT(x) / N0) }, { name: 'ripening, (r₀/r)³', x: ts, y: ts.map(x => (r0 / rOf(x)) ** 3) }], explain: 'Both mechanisms reduce the particle count: aggregation by collisions, ripening by dissolution of the smallest particles at constant solid volume.' },
      { key: 'lsw', title: 'LSW stationary size distribution', xlabel: 'r / mean radius', ylabel: 'probability density', lines: [{ name: 'LSW', x: us, y: gl.map(x => x / nrm) }], explain: 'Self-similar distribution reached in diffusion-limited ripening: skewed towards small sizes with a sharp cut-off at 1.5 times the mean radius.' },
      { key: 'clus', title: 'Cluster-size distribution from the Smoluchowski equations', xlabel: 'primary particles per cluster', ylabel: 'n_k / N₀', logx: true, logy: true, lines: [{ name: `numerical, t = ${F(thEnd, 3)} τ`, x: yE.map((_, i) => i + 1), y: yE.map(x => Math.max(x, 1e-30)) }, { name: 'constant-kernel analytical', x: yE.map((_, i) => i + 1), y: yE.map((_, i) => Math.max(Math.pow(thEnd, i) / Math.pow(1 + thEnd, i + 2), 1e-30)) }], explain: 'Number of k-mers from integrating the discrete coagulation equations with the Brownian fractal kernel, compared with the classical constant-kernel solution.' },
    ],
    equations: [...chemEquations(ck),
      eqn('governing', 'LSW ripening (cube law)', 'r̄³(t) − r̄₀³ = K t,   K = 8 γ D C_∞ V_m² / (9 R T)'),
      eqn('governing', 'Smoluchowski coagulation', 'dn_k/dt = ½ Σ_{i+j=k} β_ij n_i n_j − n_k Σ_i β_ik n_i'),
      eqn('constitutive', 'Brownian kernel with stability ratio', 'β_ij = (2kT / 3μW) (i^(1/d_f) + j^(1/d_f)) (i^(−1/d_f) + j^(−1/d_f))'),
      eqn('formula', 'Number decay and half-time', 'N(t) = N₀ / (1 + t/τ),   τ = 3 μ W / (4 kT N₀)')],
    checks, warnings,
    recommendations: [
      rec('high', dominant === 'aggregation' ? `Aggregation is the faster route (half-time ${F(tau, 3)} s). Raise the stability ratio: move the pH away from the isoelectric point, lower the ionic strength or add a steric stabiliser; a tenfold increase in W gives a tenfold longer shelf life.` : `Ripening is the faster route (radius doubles in ${F(tDoubleR, 3)} s). Lower the solubility (pH of minimum solubility, poorer solvent), cool the suspension or cap the surface with ligands to lower the interfacial energy.`, 'Comparison of size-doubling times of the two mechanisms.'),
      rec('medium', `Diluting the sol by a factor of 10 extends the aggregation half-time to ${F(10 * tau, 3)} s but does not affect the ripening rate per particle.`, 'Second-order number dependence of Smoluchowski kinetics.'),
      rec('medium', 'Follow the mean size by dynamic light scattering against time: a slope of ⅓ on a log–log plot indicates ripening, a slope of 1/d_f (about 0.5) indicates aggregation.', 'Different power laws of the two mechanisms.'),
      rec('low', 'Storing concentrates cold and stabilised avoids re-dispersion by sonication and re-synthesis of off-specification batches, which saves energy and precursor.', 'Arrhenius and solubility dependence of the coarsening rates.', true),
    ],
    publish: { 'colloid.number_density': Nend, 'colloid.t_half_s': tau, 'geom.diameter': 2 * Math.max(rT, rAgg), 'growth.rate': (rT - r0) / tEnd, 'kin.K_lsw': K },
    numerics: { method: 'Closed-form LSW and Smoluchowski solutions; Dormand–Prince RK45 on 60 discrete cluster classes', dof: nc, steps: sm.steps || 0, iterations: sm.steps || 0, residual: Math.abs(1 - massIn), notes: 'Cluster classes above 60 primary particles are truncated; the retained mass is reported as a check.' },
  };
}

// ---------------------------------------------------------------- JMAK
function runJMAK(inp) {
  const ck = sel(inp, 'chem'), T = num(inp, 'Tj') + 273.15, tEnd = num(inp, 'tmin') * 60, n = num(inp, 'nAv'), k0 = num(inp, 'k0'), Ea = num(inp, 'EaJ') * 1e3, drive = sel(inp, 'drive'), Teq = num(inp, 'Teq'), Bn = num(inp, 'Bn'), warnings = [], checks = [];
  const kOf = Tk => { let k = k0 * Math.exp(-Ea / (R * Tk)); if (drive === 'undercool') k = Tk >= Teq ? 0 : k * Math.exp(-Math.min(Bn * Teq * Teq / (Tk * (Teq - Tk) ** 2), 700)); return k; };
  const k = kOf(T), X = tt => 1 - Math.exp(-Math.pow(k * tt, n)), tX = (x, kk) => (kk > 0 ? Math.pow(-Math.log(1 - x), 1 / n) / kk : Infinity), Xend = X(tEnd);
  const t50 = tX(0.5, k), t99 = tX(0.99, k), t01 = tX(0.01, k), tPeak = n > 1 && k > 0 ? Math.pow((n - 1) / n, 1 / n) / k : 0;
  const cl = v => N.clamp(Math.log10(Math.max(v, 1e-30)), -6, 20);
  const Ts = drive === 'undercool' ? N.linspace(0.3 * Teq, 0.995 * Teq, 120) : N.linspace(Math.max(250, T - 200), T + 200, 120);
  const ttt = [0.01, 0.5, 0.99].map(x => ({ name: `${100 * x} % transformed`, x: Ts.map(Tk => cl(tX(x, kOf(Tk)))), y: Ts.map(Tk => Tk - 273.15) }));
  const iNose = N.argmin(ttt[1].x), Tnose = Ts[iNose];
  const tHi = Number.isFinite(t99) ? Math.max(1.5 * t99, 1e-9) : tEnd, tt = N.linspace(0, Math.min(Math.max(tHi, tEnd), 1e12), 200);
  const lt = tt.slice(1).filter(x => { const xx = X(x); return xx > 1e-6 && xx < 1 - 1e-9; });
  const slope = lt.length > 2 ? N.linreg(lt.map(Math.log), lt.map(x => Math.log(-Math.log(1 - X(x))))).slope : n;
  checks.push(chk('Transformed fraction within [0, 1]', Xend >= 0 && Xend <= 1, `X = ${F(Xend, 5)}`));
  checks.push(chk('Avrami slope recovered from ln[−ln(1−X)] vs ln t', Math.abs(slope - n) < 1e-6 * n + 1e-6, `slope ${F(slope, 6)} vs n = ${F(n, 6)}`));
  checks.push(chk('Hold temperature below the equilibrium temperature', drive !== 'undercool' || T < Teq, drive === 'undercool' ? `T = ${F(T, 4)} K, T_eq = ${F(Teq, 4)} K` : 'not applicable (Arrhenius-only option)'));
  if (k <= 0) warnings.push('The hold temperature is at or above the equilibrium temperature: there is no driving force and nothing transforms.');
  if (!Number.isFinite(t50) || t50 > 1e15) warnings.push('The transformation is effectively frozen at this temperature (half-time beyond 10¹⁵ s).');
  warnings.push('JMAK assumes random nucleation in an infinite sample, isotropic growth and a constant rate law; in thin films and nanoparticles the finite size lowers the effective exponent. Parameters must be fitted to isothermal data (e.g. in-situ XRD or DSC).');
  const fin = v => (Number.isFinite(v) && v < 1e30 ? v : 'not reached');
  const mech = n >= 3.5 ? 'continuous nucleation, 3-D interface-controlled growth' : n >= 2.75 ? 'site-saturated 3-D growth or continuous nucleation with 2-D growth' : n >= 2.25 ? 'continuous nucleation with 3-D diffusion-controlled growth' : n >= 1.75 ? 'site-saturated 2-D growth' : n >= 1.25 ? 'site-saturated diffusion-controlled 3-D growth' : 'one-dimensional or strongly impinged growth';
  return {
    summary: `At ${F(T - 273.15, 4)} °C the JMAK rate constant is ${F(k, 3)} 1/s, so ${F(100 * Xend, 4)} % is transformed after ${F(tEnd / 60, 4)} min; the half-transformation time is ${typeof fin(t50) === 'number' ? F(t50, 3) + ' s' : 'not reached'}. An Avrami exponent of ${F(n, 3)} is consistent with ${mech}.${drive === 'undercool' ? ` The nose of the C-curve lies near ${F(Tnose - 273.15, 4)} °C.` : ''}`,
    kpis: [
      kpi('X', 'Transformed fraction at the end of the hold', Xend, '', 'X = 1 − exp[−(kt)ⁿ] at the stated time and temperature; the predicted crystalline (product) fraction.', Xend > 0.95),
      kpi('k', 'JMAK rate constant', k, '1/s', 'k = k₀ exp(−Ea/RT), multiplied by the undercooling factor when that option is selected. Its reciprocal is the characteristic transformation time.'),
      kpi('t50', 'Time to 50 % transformation', fin(t50), 's', 'Half-transformation time (ln 2)^(1/n)/k; the usual measure of transformation speed on a TTT diagram.'),
      kpi('t99', 'Time to 99 % transformation', fin(t99), 's', 'Hold time needed for practically complete transformation at this temperature.'),
      kpi('t01', 'Incubation time (1 % transformed)', fin(t01), 's', 'Time at which 1 % has transformed; the apparent start of the transformation.'),
      kpi('tpeak', 'Time of maximum transformation rate', tPeak, 's', 'Maximum of dX/dt, at ((n−1)/n)^(1/n)/k for n > 1; zero means the rate is highest at the start.'),
      kpi('mech', 'Mechanism suggested by the exponent', mech, '', 'Textbook interpretation of the Avrami exponent in terms of nucleation mode and growth dimensionality; it is not unique and should be confirmed by microscopy.'),
      kpi('Tnose', drive === 'undercool' ? 'Nose temperature of the C-curve' : 'Temperature for 50 % in the stated time', drive === 'undercool' ? Tnose - 273.15 : (() => { const kreq = Math.pow(Math.LN2, 1 / n) / tEnd; return Ea / (R * Math.log(k0 / kreq)) - 273.15; })(), '°C', drive === 'undercool' ? 'Temperature of fastest transformation, where the gain in mobility on heating balances the loss of driving force.' : 'Temperature at which half the material transforms within the stated hold time, from inverting the Arrhenius rate constant.'),
    ],
    series: [
      { key: 'xt', title: 'Transformed fraction vs time', xlabel: 'time (s)', ylabel: 'fraction transformed', lines: [{ name: `X(t) at ${F(T - 273.15, 4)} °C`, x: tt, y: tt.map(X) }], explain: 'Sigmoidal JMAK curve: incubation, rapid transformation, then slowing by impingement of growing regions.' },
      { key: 'ttt', title: 'Time–temperature–transformation diagram', xlabel: 'log₁₀ time (s)', ylabel: 'temperature (°C)', lines: ttt, explain: 'Iso-conversion lines (1, 50 and 99 %); log time is clipped to −6…20. With the undercooling option the lines form a C-curve whose nose marks the fastest transformation.' },
      { key: 'av', title: 'Avrami plot', xlabel: 'ln t', ylabel: 'ln[−ln(1 − X)]', lines: [{ name: 'model', x: lt.length ? lt.map(Math.log) : [0, 1], y: lt.length ? lt.map(x => Math.log(-Math.log(1 - X(x)))) : [0, 0] }], explain: 'Linearised form used to extract n (slope) and k (intercept n ln k) from measured isothermal data.' },
    ],
    equations: [...chemEquations(ck),
      eqn('governing', 'Kolmogorov–Johnson–Mehl–Avrami equation', 'X(t) = 1 − exp[−(k t)ⁿ]'),
      eqn('constitutive', 'Rate constant', drive === 'undercool' ? 'k(T) = k₀ exp(−E_a/RT) · exp[−B T_eq² / (T (T_eq − T)²)]' : 'k(T) = k₀ exp(−E_a/RT)'),
      eqn('formula', 'Iso-conversion time', 't_X = [−ln(1 − X)]^(1/n) / k(T)')],
    checks, warnings,
    recommendations: [
      rec('high', Xend < 0.95 ? `Only ${F(100 * Xend, 3)} % is transformed. ${typeof fin(t99) === 'number' ? `Hold for ${F(t99, 3)} s at this temperature` : 'Raise the temperature'}${drive === 'undercool' ? `, or move towards the nose at ${F(Tnose - 273.15, 4)} °C` : ` or raise the temperature: every 10 K multiplies k by ${F(Math.exp(10 * Ea / (R * T * T)), 3)}`}.` : `Transformation is essentially complete (${F(100 * Xend, 4)} %); the hold could be shortened to about ${F(t99, 3)} s.`, 'Transformed fraction and iso-conversion times.'),
      rec('medium', 'For nanostructures choose the lowest temperature that completes the transformation in an acceptable time: grain growth and sintering have their own (usually higher) activation energies and are suppressed more strongly by lower temperature.', 'Arrhenius temperature dependence of competing processes.'),
      rec('medium', 'Calibrate n, k₀ and E_a from isothermal measurements at three or more temperatures (Avrami plots, then an Arrhenius plot of k) with the fitting mode of this suite.', 'Model parameters are empirical.'),
      rec('low', Xend >= 0.99 && Number.isFinite(t99) && t99 < 0.5 * tEnd ? `The hold is ${F(tEnd / t99, 3)} times longer than needed for 99 % transformation; shortening it saves furnace energy in proportion.` : 'Use rapid thermal annealing or a tighter ramp–soak schedule once kinetics are calibrated, to minimise furnace energy per batch.', 'Ratio of hold time to t₉₉.', true),
    ],
    publish: { 'rxn.conversion': Xend, 'kin.k': k, 'kin.t50_s': Number.isFinite(t50) ? Math.min(t50, 1e30) : 1e30, 'proc.T_C': T - 273.15 },
    numerics: { method: 'Closed-form JMAK kinetics with Arrhenius rate constant', dof: 1, steps: tt.length, iterations: 1, residual: Math.abs(slope - n), notes: 'Isothermal; non-isothermal paths require the additivity rule.' },
  };
}

// ---------------------------------------------------------------- Arrhenius fit
function parseTable(text) {
  const rows = [];
  for (const line of String(text ?? '').split(/\r?\n/)) { const s = line.trim(); if (!s || s.startsWith('#')) continue; const p = s.split(/[,;\t]+|\s{1,}/).filter(Boolean).map(Number); if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1])) rows.push([p[0], p[1]]); }
  return rows;
}
function runArrFit(inp) {
  const ck = sel(inp, 'chem'), Tp = num(inp, 'T') + 273.15, raw = inp.data ?? DEMO_TABLE, rows = parseTable(raw).filter(r => r[0] > -273.15 && r[1] > 0), warnings = [], checks = [];
  if (rows.length < 3 || new Set(rows.map(r => r[0])).size < 3) throw new Error('Arrhenius fit needs at least three rows of T_C,k with positive rate constants at three different temperatures.');
  if (String(raw).trim() === DEMO_TABLE) warnings.push('The table holds synthetic demonstration data — replace with measurements.');
  const Tk = rows.map(r => r[0] + 273.15), kk = rows.map(r => r[1]), x = Tk.map(t => 1 / t), y = kk.map(Math.log), n = rows.length;
  const lr = N.linreg(x, y), mx = N.mean(x), Sxx = N.sum(x.map(v => (v - mx) ** 2)), res = y.map((v, i) => v - (lr.intercept + lr.slope * x[i])), s2 = N.sum(res.map(v => v * v)) / Math.max(1, n - 2), s = Math.sqrt(s2);
  const seSl = Math.sqrt(s2 / Sxx), seInt = s * Math.sqrt(1 / n + mx * mx / Sxx), Ea = -lr.slope * R, seEa = seSl * R, lnA = lr.intercept, A = Math.exp(N.clamp(lnA, -700, 700));
  const z = 1.959964, nu = Math.max(1, n - 2), tq = z + (z ** 3 + z) / (4 * nu) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * nu * nu);
  const pred = xx => lr.intercept + lr.slope * xx, sePred = xx => s * Math.sqrt(1 / n + (xx - mx) ** 2 / Sxx), lnkP = pred(1 / Tp), kP = Math.exp(N.clamp(lnkP, -700, 700)), ciF = Math.exp(Math.min(tq * sePred(1 / Tp), 700));
  // nonlinear least squares in k-space
  const km = N.mean(kk), lm = safe(() => N.levenbergMarquardt(p => Tk.map((t, i) => (Math.exp(N.clamp(p[0] - p[1] * 1000 / (R * t), -700, 700)) - kk[i]) / km), [lnA, Ea / 1000], { maxIter: 60 }), null);
  const EaNL = lm ? lm.p[1] * 1000 : Ea, seNL = lm && lm.stderr ? lm.stderr[1] * 1000 : 0;
  // Eyring
  const ey = N.linreg(x, kk.map((v, i) => Math.log(v / Tk[i]))), dHe = -ey.slope * R, dSe = (ey.intercept - Math.log(N.K_B / N.H_PLANCK)) * R;
  const Tm = N.mean(Tk), span = N.maxOf(Tk) - N.minOf(Tk), extrap = Tp < N.minOf(Tk) || Tp > N.maxOf(Tk);
  checks.push(chk('Linearity of ln k vs 1/T (R² ≥ 0.98)', lr.r2 >= 0.98, `R² = ${F(lr.r2, 5)}`));
  checks.push(chk('Activation energy positive and significant (E_a > 2 standard errors)', Ea > 2 * seEa, `E_a = ${F(Ea / 1000, 4)} ± ${F(seEa / 1000, 3)} kJ/mol`));
  checks.push(chk('Temperature span at least 20 K', span >= 20, `span ${F(span, 4)} K over ${n} points`));
  checks.push(chk('Linearised and nonlinear fits agree within two standard errors', Math.abs(EaNL - Ea) <= 2 * seEa + 1e-9 * Math.abs(Ea), `nonlinear E_a = ${F(EaNL / 1000, 4)} kJ/mol`));
  checks.push(chk('Prediction temperature inside the measured range', !extrap, `${F(Tp - 273.15, 4)} °C vs data ${F(N.minOf(Tk) - 273.15, 4)}–${F(N.maxOf(Tk) - 273.15, 4)} °C`));
  if (extrap) warnings.push('The prediction temperature lies outside the measured range: the confidence interval widens and a change of mechanism cannot be excluded.');
  if (n < 5) warnings.push('Fewer than five data points: standard errors are poorly determined.');
  if (Ea <= 0) warnings.push('The fitted activation energy is not positive: the data do not follow simple Arrhenius behaviour (transport limitation, equilibrium shift or measurement scatter).');
  warnings.push('Eyring activation entropy assumes a first-order rate constant in 1/s; for other units only ΔH‡ is meaningful.');
  const xs = N.linspace(N.minOf(x) - 0.05 * (N.maxOf(x) - N.minOf(x) + 1e-9), N.maxOf(x) + 0.05 * (N.maxOf(x) - N.minOf(x) + 1e-9), 60), X3 = xs.map(v => 1000 * v);
  return {
    summary: `A linear Arrhenius fit of ${n} rate constants between ${F(N.minOf(Tk) - 273.15, 4)} and ${F(N.maxOf(Tk) - 273.15, 4)} °C gives E_a = ${F(Ea / 1000, 4)} ± ${F(seEa / 1000, 3)} kJ/mol and ln A = ${F(lnA, 4)} ± ${F(seInt, 3)} (R² = ${F(lr.r2, 4)}). The predicted rate constant at ${F(Tp - 273.15, 4)} °C is ${F(kP, 3)} with a 95 % confidence factor of ×/÷ ${F(ciF, 3)}.`,
    kpis: [
      kpi('Ea', 'Activation energy', Ea / 1000, 'kJ/mol', 'From the slope −E_a/R of ln k against 1/T (ordinary least squares). Typical values: 10–25 kJ/mol for diffusion control, 40–100 kJ/mol for surface-reaction control.', Ea > 0),
      kpi('Ea_se', 'Standard error of E_a', seEa / 1000, 'kJ/mol', `One standard error from the regression; the 95 % confidence half-width is ${F(tq * seEa / 1000, 3)} kJ/mol with ${nu} degrees of freedom.`),
      kpi('lnA', 'ln A (pre-exponential factor)', lnA, 'ln(units of k)', `Intercept of the Arrhenius plot; A = ${F(A, 3)} in the units of k. Strongly correlated with E_a, so always quote both together.`),
      kpi('lnA_se', 'Standard error of ln A', seInt, '', 'One standard error of the intercept. It is large when the temperature span is narrow because the intercept is extrapolated to 1/T = 0.'),
      kpi('r2', 'Coefficient of determination R²', lr.r2, '', 'Fraction of the variance of ln k explained by the straight line; curvature (low R²) signals a change of mechanism or transport limitation.', lr.r2 >= 0.98),
      kpi('k_pred', 'Predicted rate constant at the process temperature', kP, 'units of k', `Value of the fitted line at ${F(Tp - 273.15, 4)} °C; 95 % confidence range ${F(kP / ciF, 3)} to ${F(kP * ciF, 3)}.`),
      kpi('Ea_nl', 'E_a from nonlinear fit in k', EaNL / 1000, 'kJ/mol', `Levenberg–Marquardt fit of k = A exp(−E_a/RT) without the log transform (standard error ${F(seNL / 1000, 3)} kJ/mol). It weights high-temperature points more heavily; a large difference from the linear fit indicates heteroscedastic data.`),
      kpi('dH', 'Eyring activation enthalpy ΔH‡', dHe / 1000, 'kJ/mol', 'From the slope of ln(k/T) against 1/T; equals E_a − RT at the mean temperature of the data.'),
      kpi('dS', 'Eyring activation entropy ΔS‡', dSe, 'J/mol/K', 'From the intercept of the Eyring plot assuming a first-order rate constant in 1/s. Negative values indicate an ordered (associative) transition state.'),
      kpi('q10', 'Rate increase per 10 K', Math.exp(N.clamp(10 * Ea / (R * Tp * (Tp + 10)), -700, 700)), '', 'Factor by which the rate constant grows when the temperature is raised by 10 K from the process temperature.'),
    ],
    series: [
      { key: 'arr', title: 'Arrhenius plot', xlabel: '1000/T (1/K)', ylabel: 'ln k', lines: [{ name: 'data', x: x.map(v => 1000 * v), y }, { name: 'linear fit', x: X3, y: xs.map(pred) }, { name: '95 % upper', x: X3, y: xs.map(v => pred(v) + tq * sePred(v)) }, { name: '95 % lower', x: X3, y: xs.map(v => pred(v) - tq * sePred(v)) }], explain: 'Measured ln k against reciprocal temperature with the least-squares line and its 95 % confidence band. Systematic curvature of the points means a single activation energy is not adequate.' },
      { key: 'res', title: 'Residuals of the fit', xlabel: 'temperature (°C)', ylabel: 'ln k (data) − ln k (fit)', lines: [{ name: 'residual', x: Tk.map(t => t - 273.15), y: res }, { name: 'zero', x: [N.minOf(Tk) - 273.15, N.maxOf(Tk) - 273.15], y: [0, 0] }], explain: 'Residuals should scatter randomly about zero; a trend (e.g. a bow shape) indicates model inadequacy rather than noise.' },
    ],
    tables: [{ key: 'fit', title: 'Data and fitted values', columns: ['T (°C)', 'k (data)', 'k (fit)', 'relative deviation (%)'], rows: rows.map((r, i) => [r[0], +r[1].toExponential(4), +Math.exp(pred(x[i])).toExponential(4), +(100 * (r[1] / Math.exp(pred(x[i])) - 1)).toFixed(2)]), explain: 'Point-by-point comparison of the measurements with the fitted Arrhenius law.' }],
    equations: [...chemEquations(ck),
      eqn('governing', 'Arrhenius equation', 'k = A exp(−E_a / R T);   ln k = ln A − (E_a/R)(1/T)'),
      eqn('governing', 'Eyring transition-state equation', 'k = (k_B T / h) exp(ΔS‡/R) exp(−ΔH‡/RT)'),
      eqn('formula', 'Least-squares standard errors', 'se(slope) = s / √S_xx;   se(intercept) = s √(1/n + x̄²/S_xx);   s² = Σ residual² / (n − 2)')],
    checks, warnings,
    recommendations: [
      rec('high', seEa / Math.abs(Ea || 1) > 0.1 ? `E_a is known only to ±${F(100 * seEa / Math.abs(Ea || 1), 3)} % (one standard error). Add measurements at the extremes of the temperature range: the standard error of the slope falls with the spread of 1/T, so widening the span from ${F(span, 3)} K helps more than repeating mid-range points.` : `E_a is determined to ±${F(100 * seEa / Math.abs(Ea || 1), 3)} %; the fit is adequate for design within the measured range.`, 'Standard error of the regression slope.'),
      rec('medium', Ea / 1000 < 25 && Ea > 0 ? 'An activation energy below about 25 kJ/mol suggests mass-transfer control: check for a stirring-speed dependence before treating it as intrinsic kinetics.' : 'An activation energy in the chemical range suggests intrinsic kinetics; verify by confirming that the rate is independent of stirring speed and particle size.', 'Magnitude of the fitted activation energy.'),
      rec('medium', `Temperature control matters: at ${F(Tp - 273.15, 4)} °C a 1 K error changes the rate by ${F(100 * Ea / (R * Tp * Tp), 3)} %.`, 'Derivative d ln k/dT = E_a/RT².'),
      rec('low', `Use the fitted law to trade temperature against time: lowering the temperature by 10 K lengthens the batch by a factor ${F(Math.exp(N.clamp(10 * Ea / (R * Tp * (Tp - 10)), -700, 700)), 3)} but reduces heat losses; choose the minimum-energy combination for the reactor's heat-loss coefficient.`, 'Arrhenius time–temperature equivalence.', true),
    ],
    publish: { 'kin.Ea': Ea, 'kin.lnA': lnA, 'kin.k': kP, 'cred.rmse': s },
    numerics: { method: 'Ordinary least squares on ln k vs 1/T with analytical standard errors; Levenberg–Marquardt nonlinear fit in k for comparison', dof: 2, steps: n, iterations: lm ? lm.iterations : 0, residual: s, notes: `Student-t quantile ${F(tq, 4)} (Cornish–Fisher approximation) for ${nu} degrees of freedom; mean data temperature ${F(Tm, 4)} K.` },
  };
}

const RUN = { zno: runZnO, lamer: runLaMer, ripen: runRipen, jmak: runJMAK, arrfit: runArrFit };

export default {
  id: 7,
  tagline: 'Reaction kinetics, nucleation and growth: seeded ZnO nanorod growth, LaMer burst nucleation with a population balance, ripening and aggregation, JMAK kinetics and Arrhenius fitting.',
  inputs,
  presets: [
    { name: 'ZnO nanorods, zinc nitrate + HMTA, 90 °C, 3 h', note: 'Standard equimolar 25 mmol/L chemical-bath growth on a seeded substrate.', values: { mode: 'zno', chem: 'zn_nitrate', cZn: 25, cHMTA: 25, cOH: 0, T: 90, tmin: 180, E: 0 } },
    { name: 'ZnO hydrothermal reference: 0.01 M Zn²⁺ + 0.02 M OH⁻, 1000 V/m', note: 'Reference case: direct hydroxide, E_a = 65 kJ/mol, D(Zn²⁺) = 7.5×10⁻¹⁰ and D(OH⁻) = 5.3×10⁻⁹ m²/s, heated from 298 K, with an applied field.', values: { mode: 'zno', chem: 'zn_nitrate', cZn: 10, cHMTA: 0, cOH: 20, T: 90, T0: 25, tmin: 180, Ea: 65, DZn: 7.5e-10, DOH: 5.3e-9, E: 1000 } },
    { name: 'Gold nanoparticles by citrate reduction (Turkevich type)', note: 'Burst nucleation of Au at the boil; 0.25 mmol/L HAuCl₄.', values: { mode: 'lamer', chem: 'au', T: 100, tmin: 15, cP: 0.25, kf: 0.02, gam: 0.45, cInf: 1e-6, Dm: 1.5e-9, ks: 2e-2, logA: 34, Vb: 1e-4 } },
    { name: 'Silver nanoparticles by borohydride reduction', note: 'Very fast reduction at room temperature gives many small nuclei.', values: { mode: 'lamer', chem: 'ag', T: 25, tmin: 5, cP: 1, kf: 5, gam: 0.4, cInf: 1e-6, Dm: 1e-9, ks: 1e-2, logA: 34, Vb: 1e-4 } },
    { name: 'Ageing of a ZnO nanoparticle sol', note: 'Ripening against aggregation of 8 nm particles during storage at 60 °C.', values: { mode: 'ripen', chem: 'zn_nitrate', T: 60, tmin: 600, r0: 4, gam: 0.25, cInf: 2e-3, Dm: 1e-9, N0: 1e20, W: 1e6, mu: 4.7e-4, df: 1.8 } },
    { name: 'Crystallisation of amorphous TiO₂ to anatase at 400 °C', note: 'Isothermal annealing of a sol–gel film described by JMAK kinetics.', values: { mode: 'jmak', chem: 'tio2', Tj: 400, tmin: 60, nAv: 2.5, k0: 5e8, EaJ: 150, drive: 'arrhenius' } },
  ],
  run(inp) { return RUN[sel(inp, 'mode')](inp); },
};
