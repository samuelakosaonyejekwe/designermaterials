// Suite 20 — Active tuning, programmability and control.
// Closed-loop synthesis control (PID with tuning rules, anti-windup, margins; constrained MPC of nanorod growth), LQR/Kalman
// design for a resonator, and reduced-order models of physical tuning mechanisms (varactor, phase-change, piezoelectric,
// liquid-crystal, electrostatic MEMS, ferrite, bistable arrays, quasi-zero-stiffness isolators).
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';

const F = N.fmt;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => note ? { kind, label, text, note } : { kind, label, text };
const rec = (priority, text, basis, sustainability) => sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis };
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const line = (name, x, y) => { const X = [], Y = [], st = Math.max(1, Math.ceil(x.length / 400)); for (let i = 0; i < x.length; i += st) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) { X.push(x[i]); Y.push(y[i]); } if (!X.length) { X.push(0); Y.push(0); } return { name, x: X, y: Y }; };
const mechTable = rows => ({ key: 'mechanism', title: 'Tuning-mechanism assessment', columns: ['Aspect', 'Assessment for this design'], rows, explain: 'Engineering assessment of the physical tuning mechanism: what drives it, what it costs in energy and time, and what limits its reversibility and life.' });
const formulaEq = (f, note) => { try { return [eqn('formula', 'Active material', `${chem.prettyFormula(f)}  (M = ${F(chem.molarMass(f), 5)} g/mol)`, note)]; } catch { return []; } };
const GAMMA = 2 * N.PI * 28.025e9; // gyromagnetic ratio, rad/s/T

const MODES = [
  { value: 'pid', label: 'PID temperature control of a synthesis reactor' },
  { value: 'mpc', label: 'Model-predictive control of nanorod growth' },
  { value: 'lqg', label: 'LQR + Kalman filter for a resonator' },
  { value: 'varactor', label: 'Varactor-loaded split-ring resonator' },
  { value: 'pcm', label: 'Phase-change metasurface (VO₂ / GST)' },
  { value: 'piezo', label: 'Piezoelectric stack and shunt stiffness tuning' },
  { value: 'lc', label: 'Liquid-crystal index tuning' },
  { value: 'mems', label: 'Electrostatic MEMS tuning and pull-in' },
  { value: 'ferrite', label: 'Magnetically biased ferrite' },
  { value: 'bistable', label: 'Programmable bistable mechanical array' },
  { value: 'qzs', label: 'Quasi-zero-stiffness isolator tuning' },
];
const S = m => ({ mode: Array.isArray(m) ? m : [m] });
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'pid', options: MODES, group: 'Model', help: 'Selects the control problem or tuning mechanism. Each analysis uses the inputs shown for it.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 9999, step: 1, sens: false, group: 'Model', showIf: S(['pid', 'mpc', 'lqg', 'bistable']), help: 'Seed for sensor noise, process noise and manufacturing variation.' },
  // PID
  { key: 'V_reac', label: 'Reactor liquid volume', unit: 'm³', value: 1e-3, min: 1e-5, max: 50, log: true, group: 'Plant', bus: 'proc.volume', showIf: S('pid'), help: 'Liquid volume; with the vessel heat capacity it sets the thermal mass and time constant.' },
  { key: 'P_heater', label: 'Heater power', unit: 'W', value: 400, min: 1, max: 5e6, log: true, group: 'Plant', showIf: S('pid'), help: 'Maximum heater power (controller output 100 %).' },
  { key: 'UA', label: 'Heat-loss conductance UA', unit: 'W/K', value: 2, min: 0.01, max: 1e5, log: true, group: 'Plant', bus: 'proc.UA', showIf: S('pid'), help: 'Overall loss conductance from liquid to ambient. Steady-state gain is P/UA and the time constant C/UA.' },
  { key: 'C_vessel', label: 'Vessel heat capacity', unit: 'J/K', value: 1500, min: 0, max: 1e8, group: 'Plant', showIf: S('pid'), help: 'Heat capacity of the vessel wall and internals that heat up with the liquid.' },
  { key: 'dead', label: 'Dead time', unit: 's', value: 45, min: 1, max: 5000, log: true, group: 'Plant', bus: 'ctrl.plant_delay', showIf: S('pid'), help: 'Apparent delay between a heater change and the sensor response (wall lag, mixing, sensor).' },
  { key: 'T_set', label: 'Temperature set-point', unit: '°C', value: 90, min: 30, max: 350, group: 'Plant', bus: 'proc.T_C', showIf: S('pid'), help: 'Target reaction temperature.' },
  { key: 'T_amb', label: 'Ambient / start temperature', unit: '°C', value: 25, min: -20, max: 60, group: 'Plant', showIf: S('pid'), help: 'Ambient temperature, also the initial temperature of the charge.' },
  { key: 'tuning', label: 'Tuning rule', type: 'select', value: 'simc', options: [{ value: 'simc', label: 'SIMC (Skogestad) PI' }, { value: 'imc', label: 'IMC PID' }, { value: 'zn', label: 'Ziegler–Nichols open-loop PID' }, { value: 'cc', label: 'Cohen–Coon PID' }, { value: 'manual', label: 'Manual gains' }], group: 'Controller', showIf: S('pid'), help: 'Rule used to compute the gains from the first-order-plus-dead-time plant model.' },
  { key: 'lam', label: 'Closed-loop time constant / dead time', value: 1, min: 0.25, max: 20, step: 0.25, group: 'Controller', showIf: S('pid'), help: 'Desired closed-loop time constant as a multiple of the dead time (SIMC and IMC only). Larger values are slower and more robust.' },
  { key: 'Kc_man', label: 'Manual proportional gain', unit: '1/K', value: 0.05, min: 1e-5, max: 100, log: true, group: 'Controller', showIf: S('pid'), help: 'Heater fraction per kelvin of error (manual tuning only).' },
  { key: 'Ti_man', label: 'Manual integral time', unit: 's', value: 600, min: 1, max: 1e6, log: true, group: 'Controller', showIf: S('pid'), help: 'Integral (reset) time (manual tuning only).' },
  { key: 'Td_man', label: 'Manual derivative time', unit: 's', value: 0, min: 0, max: 1e4, group: 'Controller', showIf: S('pid'), help: 'Derivative time (manual tuning only); 0 gives a PI controller.' },
  { key: 'aw', label: 'Anti-windup', type: 'select', value: 'on', options: [{ value: 'on', label: 'Back-calculation anti-windup' }, { value: 'off', label: 'None (integrator winds up)' }], group: 'Controller', showIf: S('pid'), help: 'Back-calculation stops the integrator accumulating while the heater is saturated during heat-up.' },
  { key: 'Ts', label: 'Sampling interval', unit: 's', value: 2, min: 0.05, max: 600, log: true, sens: false, group: 'Controller', showIf: S('pid'), help: 'Controller sampling period; it should be below about one fifth of the dead time.' },
  { key: 'noise', label: 'Sensor noise (rms)', unit: 'K', value: 0.1, min: 0, max: 5, step: 0.05, group: 'Controller', showIf: S('pid'), help: 'Standard deviation of the temperature measurement noise.' },
  { key: 'dist', label: 'Load disturbance', unit: '% of heater', value: 20, min: 0, max: 80, step: 1, group: 'Controller', showIf: S('pid'), help: 'Step increase in heat loss (for example cold reagent dosing) applied at 60 % of the simulated time.' },
  { key: 't_sim', label: 'Simulated time', unit: 's', value: 7200, min: 60, max: 1e6, log: true, sens: false, group: 'Controller', showIf: S('pid'), help: 'Duration of the closed-loop simulation.' },
  // MPC
  { key: 'L_target', label: 'Target nanorod length', unit: 'm', value: 2e-6, min: 2e-8, max: 1e-4, log: true, group: 'Growth control', bus: 'growth.length', showIf: S('mpc'), help: 'Final length the growth run must reach.' },
  { key: 't_growth', label: 'Growth time', unit: 's', value: 14400, min: 600, max: 4e5, log: true, group: 'Growth control', showIf: S('mpc'), help: 'Duration of the growth run; the reference length rises linearly over it.' },
  { key: 'G_ref', label: 'Model growth rate at 90 °C', unit: 'm/s', value: 1.5e-10, min: 1e-13, max: 1e-7, log: true, group: 'Growth control', bus: 'growth.rate', showIf: S('mpc'), help: 'Axial growth rate predicted by the controller model at the reference temperature and nominal precursor concentration.' },
  { key: 'Ea', label: 'Growth activation energy', unit: 'kJ/mol', value: 60, min: 10, max: 200, group: 'Growth control', showIf: S('mpc'), help: 'Arrhenius activation energy of the growth rate (temperature manipulation only).' },
  { key: 'manip', label: 'Manipulated variable', type: 'select', value: 'temp', options: [{ value: 'temp', label: 'Reactor temperature' }, { value: 'dose', label: 'Precursor dosing (relative concentration)' }], group: 'Growth control', showIf: S('mpc'), help: 'Input the controller adjusts: temperature (Arrhenius rate) or precursor concentration (first-order rate).' },
  { key: 'T_lo', label: 'Minimum temperature', unit: '°C', value: 65, min: 25, max: 150, group: 'Growth control', showIf: S('mpc'), help: 'Lower temperature constraint (growth morphology and nucleation limits).' },
  { key: 'T_hi', label: 'Maximum temperature', unit: '°C', value: 97, min: 40, max: 250, group: 'Growth control', showIf: S('mpc'), help: 'Upper temperature constraint (boiling, pressure rating, secondary nucleation).' },
  { key: 'bias', label: 'Plant–model rate mismatch', unit: '%', value: -15, min: -60, max: 100, step: 1, group: 'Growth control', showIf: S('mpc'), help: 'True growth rate relative to the model (for example seed-layer quality or pH offset). The controller does not know it.' },
  { key: 'deplete', label: 'Rate loss by precursor depletion', value: 0.2, min: 0, max: 0.8, step: 0.05, group: 'Growth control', showIf: S('mpc'), help: 'Fractional fall of the true growth rate when the target length is reached, not included in the controller model.' },
  { key: 'len_noise', label: 'Length measurement noise (rms)', unit: 'nm', value: 15, min: 0, max: 500, group: 'Growth control', showIf: S('mpc'), help: 'Noise of the in-situ length estimate (optical or QCM-based).' },
  { key: 'n_ctrl', label: 'Control intervals', value: 48, min: 8, max: 200, step: 1, sens: false, group: 'Growth control', showIf: S('mpc'), help: 'Number of control moves over the growth run.' },
  { key: 'N_h', label: 'Prediction horizon', value: 10, min: 2, max: 30, step: 1, group: 'Growth control', showIf: S('mpc'), help: 'Number of future intervals optimised at each step.' },
  { key: 'rho_u', label: 'Move-suppression weight', value: 1, min: 0, max: 100, step: 0.5, group: 'Growth control', showIf: S('mpc'), help: 'Penalty on changes of the manipulated variable; larger values give smoother, slower action.' },
  { key: 'tau_act', label: 'Actuator time constant', unit: 's', value: 600, min: 1, max: 2e4, log: true, group: 'Growth control', showIf: S('mpc'), help: 'First-order lag with which the reactor temperature (or concentration) follows its command; not modelled by the controller.' },
  // LQG
  { key: 'f0_r', label: 'Resonance frequency', unit: 'Hz', value: 1e4, min: 0.1, max: 1e10, log: true, group: 'Resonator', bus: 'res.f0', showIf: S('lqg'), help: 'Natural frequency of the resonator (unit cell, cantilever or isolator) to be controlled.' },
  { key: 'Q_r', label: 'Quality factor', value: 200, min: 0.6, max: 1e5, log: true, group: 'Resonator', bus: 'res.Q', showIf: S('lqg'), help: 'Open-loop quality factor; the free ring-down time is Q/(πf₀).' },
  { key: 'm_r', label: 'Effective mass', unit: 'kg', value: 1e-6, min: 1e-18, max: 1e3, log: true, group: 'Resonator', showIf: S('lqg'), help: 'Modal mass, used to convert the normalised control force to newtons.' },
  { key: 'x0_r', label: 'Initial displacement', unit: 'm', value: 1e-7, min: 1e-12, max: 0.1, log: true, group: 'Resonator', showIf: S('lqg'), help: 'Initial displacement disturbance to be damped out.' },
  { key: 'q_pos', label: 'Position weight', value: 10, min: 0, max: 1e4, group: 'Resonator', showIf: S('lqg'), help: 'LQR state weight on normalised displacement; raises the closed-loop stiffness.' },
  { key: 'q_vel', label: 'Velocity weight', value: 1, min: 0.001, max: 1e4, log: true, group: 'Resonator', showIf: S('lqg'), help: 'LQR state weight on normalised velocity; raises the closed-loop damping.' },
  { key: 'r_u', label: 'Control weight', value: 1, min: 1e-3, max: 1e3, log: true, group: 'Resonator', showIf: S('lqg'), help: 'LQR weight on the normalised actuator force; larger values economise actuation.' },
  { key: 'w_proc', label: 'Process noise (normalised rms)', value: 0.3, min: 1e-4, max: 1, log: true, group: 'Resonator', showIf: S('lqg'), help: 'Square root of the force-noise intensity in units of k·x₀ (per √(radian of phase)).' },
  { key: 'v_meas', label: 'Measurement noise (normalised rms)', value: 0.02, min: 1e-4, max: 1, log: true, group: 'Resonator', showIf: S('lqg'), help: 'Square root of the displacement-sensor noise intensity in units of x₀.' },
  // varactor
  { key: 'L_srr', label: 'Ring inductance', unit: 'nH', value: 5, min: 0.01, max: 1000, log: true, group: 'Varactor SRR', showIf: S('varactor'), help: 'Self-inductance of the split ring.' },
  { key: 'C_gap', label: 'Fixed gap capacitance', unit: 'pF', value: 0.2, min: 0.001, max: 100, log: true, group: 'Varactor SRR', showIf: S('varactor'), help: 'Geometric gap and parasitic capacitance in parallel with the varactor; it dilutes the tuning range.' },
  { key: 'Cj0', label: 'Zero-bias junction capacitance', unit: 'pF', value: 2, min: 0.01, max: 100, log: true, group: 'Varactor SRR', showIf: S('varactor'), help: 'Varactor capacitance at zero bias.' },
  { key: 'Vbi', label: 'Built-in potential', unit: 'V', value: 0.7, min: 0.2, max: 3, step: 0.05, group: 'Varactor SRR', showIf: S('varactor'), help: 'Junction built-in potential in C(V) = Cj0/(1 + V/Vbi)^M.' },
  { key: 'Mj', label: 'Grading coefficient M', value: 0.5, min: 0.3, max: 2, step: 0.05, group: 'Varactor SRR', showIf: S('varactor'), help: '0.5 for an abrupt junction, about 0.33 for linearly graded, above 1 for hyper-abrupt varactors.' },
  { key: 'R_s', label: 'Series resistance', unit: 'Ω', value: 1.5, min: 0.01, max: 100, log: true, group: 'Varactor SRR', showIf: S('varactor'), help: 'Ring ohmic plus varactor series resistance; it limits the quality factor.' },
  { key: 'V_bias', label: 'Reverse bias', unit: 'V', value: 5, min: 0, max: 30, step: 0.5, group: 'Varactor SRR', showIf: S('varactor'), help: 'Operating reverse bias voltage.' },
  { key: 'V_max', label: 'Maximum reverse bias', unit: 'V', value: 20, min: 1, max: 60, group: 'Varactor SRR', showIf: S('varactor'), help: 'Largest usable bias (below breakdown).' },
  { key: 'R_bias', label: 'Bias-line resistance', unit: 'kΩ', value: 100, min: 0.1, max: 1e4, log: true, group: 'Varactor SRR', showIf: S('varactor'), help: 'Isolation resistor of the bias network; with the junction capacitance it sets the tuning speed.' },
  // PCM
  { key: 'pcm_mat', label: 'Phase-change material', type: 'select', value: 'vo2', options: [{ value: 'vo2', label: 'VO₂ (volatile insulator–metal transition)' }, { value: 'gst', label: 'Ge₂Sb₂Te₅ (non-volatile amorphous–crystalline)' }], group: 'Phase-change metasurface', showIf: S('pcm'), help: 'VO₂ switches reversibly near 68 °C and needs holding power; GST retains its state without power and is reset by melt-quenching.' },
  { key: 'area_um2', label: 'Patch area', unit: 'µm²', value: 100, min: 0.01, max: 1e8, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'In-plane area of the phase-change patch in one unit cell.' },
  { key: 'thick_nm', label: 'Film thickness', unit: 'nm', value: 100, min: 5, max: 5000, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Thickness of the phase-change film; with conductivity it gives the sheet resistance that shorts the resonator gap.' },
  { key: 'G_th', label: 'Thermal conductance to heat sink', unit: 'W/K', value: 2e-5, min: 1e-9, max: 1, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Conductance from the patch to the substrate heat sink; sets holding power and cooling time.' },
  { key: 'Cth_mult', label: 'Heated-volume factor', value: 10, min: 1, max: 1000, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Ratio of the total heated thermal mass (film plus adjacent substrate and heater) to that of the film alone.' },
  { key: 'P_mw', label: 'Heater power per cell', unit: 'mW', value: 3, min: 1e-4, max: 1e4, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Electrical or optical heating power delivered to one cell during switching.' },
  { key: 'T_base', label: 'Base temperature', unit: '°C', value: 25, min: -50, max: 60, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Heat-sink temperature.' },
  { key: 'f0_ms', label: 'Metasurface resonance', unit: 'Hz', value: 1e12, min: 1e8, max: 1e15, log: true, group: 'Phase-change metasurface', bus: 'res.f0', showIf: S('pcm'), help: 'Resonance of the unit cell with the patch insulating.' },
  { key: 'X0', label: 'Resonator reactance √(L/C)', unit: 'Ω', value: 300, min: 10, max: 5000, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Characteristic reactance of the equivalent series RLC sheet; larger values give a narrower resonance.' },
  { key: 'R_ms', label: 'Resonator loss resistance', unit: 'Ω', value: 8, min: 0.01, max: 500, log: true, group: 'Phase-change metasurface', showIf: S('pcm'), help: 'Ohmic resistance of the metal resonator in the sheet-impedance model.' },
  // piezo
  { key: 'd33', label: 'Piezoelectric coefficient d₃₃', unit: 'pm/V', value: 400, min: 1, max: 3000, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Strain per unit field along the poling direction (soft PZT ≈ 400–650 pm/V, AlN ≈ 5 pm/V).' },
  { key: 'n_lay', label: 'Number of layers', value: 100, min: 1, max: 2000, step: 1, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Number of electrically parallel, mechanically series layers.' },
  { key: 't_lay', label: 'Layer thickness', unit: 'µm', value: 100, min: 0.1, max: 2000, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Thickness of one ceramic layer; the field is V/t.' },
  { key: 'A_pz', label: 'Cross-section', unit: 'mm²', value: 25, min: 1e-6, max: 1e4, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Active cross-sectional area of the stack.' },
  { key: 'Y_pz', label: 'Short-circuit modulus', unit: 'GPa', value: 50, min: 1, max: 400, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Elastic modulus along the stack with electrodes shorted.' },
  { key: 'k33', label: 'Coupling coefficient k₃₃', value: 0.7, min: 0.02, max: 0.95, step: 0.01, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Electromechanical coupling; the open- to short-circuit stiffness ratio is 1/(1 − k²).' },
  { key: 'epsr', label: 'Relative permittivity (free)', value: 1800, min: 5, max: 20000, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Stress-free relative permittivity of the ceramic.' },
  { key: 'V_pz', label: 'Drive voltage', unit: 'V', value: 100, min: 0.1, max: 1000, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Applied voltage for the stroke and energy estimates.' },
  { key: 'm_pz', label: 'Supported mass', unit: 'kg', value: 0.1, min: 1e-9, max: 1e3, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Mass carried by the stack; with the stack stiffness it forms the tunable resonator.' },
  { key: 'alpha_sh', label: 'Shunt capacitance / blocked capacitance', value: 1, min: -0.9, max: 100, step: 0.1, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Capacitive shunt ratio: 0 is open circuit, large values approach short circuit; negative values (active negative-capacitance circuit) stiffen beyond open circuit.' },
  { key: 'R_drv', label: 'Driver output resistance', unit: 'Ω', value: 50, min: 0.1, max: 1e5, log: true, group: 'Piezoelectric stack', showIf: S('piezo'), help: 'Source resistance charging the stack capacitance.' },
  // LC
  { key: 'n_o', label: 'Ordinary index', value: 1.52, min: 1.3, max: 2, step: 0.01, group: 'Liquid crystal', showIf: S('lc'), help: 'Ordinary refractive index of the nematic.' },
  { key: 'n_e', label: 'Extraordinary index', value: 1.74, min: 1.31, max: 2.4, step: 0.01, group: 'Liquid crystal', showIf: S('lc'), help: 'Extraordinary index; birefringence n_e − n_o is the tuning resource.' },
  { key: 'K11', label: 'Splay elastic constant', unit: 'pN', value: 11, min: 1, max: 50, group: 'Liquid crystal', showIf: S('lc'), help: 'Frank splay constant; sets threshold voltage and relaxation time.' },
  { key: 'd_eps', label: 'Dielectric anisotropy Δε', value: 11, min: 0.5, max: 60, group: 'Liquid crystal', showIf: S('lc'), help: 'Low-frequency dielectric anisotropy driving the reorientation.' },
  { key: 'gamma1', label: 'Rotational viscosity', unit: 'Pa·s', value: 0.1, min: 0.005, max: 2, log: true, group: 'Liquid crystal', showIf: S('lc'), help: 'Rotational viscosity γ₁ of the nematic.' },
  { key: 'd_lc', label: 'Cell gap', unit: 'µm', value: 5, min: 0.2, max: 200, log: true, group: 'Liquid crystal', showIf: S('lc'), help: 'Liquid-crystal layer thickness; response time scales with its square.' },
  { key: 'V_lc', label: 'Drive voltage (rms)', unit: 'V', value: 5, min: 0, max: 100, step: 0.5, group: 'Liquid crystal', showIf: S('lc'), help: 'Applied AC voltage across the cell.' },
  { key: 'lam_nm', label: 'Operating wavelength', unit: 'nm', value: 1550, min: 300, max: 3e6, log: true, group: 'Liquid crystal', showIf: S('lc'), help: 'Free-space wavelength of the resonance being tuned.' },
  { key: 'overlap', label: 'Mode overlap with the liquid crystal', value: 0.5, min: 0.01, max: 1, step: 0.01, group: 'Liquid crystal', showIf: S('lc'), help: 'Fraction of the resonant mode energy located in the liquid crystal.' },
  // MEMS
  { key: 'k_m', label: 'Suspension stiffness', unit: 'N/m', value: 10, min: 0.01, max: 1e4, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Mechanical spring constant of the suspended plate.' },
  { key: 'A_m', label: 'Electrode area', unit: 'µm²', value: 1e4, min: 1, max: 1e8, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Overlap area of the parallel-plate electrodes.' },
  { key: 'g0', label: 'Initial air gap', unit: 'µm', value: 2, min: 0.05, max: 100, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Rest gap between the plates.' },
  { key: 'f_mech', label: 'Mechanical resonance', unit: 'Hz', value: 5e4, min: 10, max: 1e9, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Unbiased mechanical resonance of the plate.' },
  { key: 't_d', label: 'Dielectric thickness', unit: 'nm', value: 100, min: 1, max: 5000, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Insulating layer on the fixed electrode.' },
  { key: 'eps_d', label: 'Dielectric permittivity', value: 7.5, min: 1, max: 100, group: 'MEMS actuator', showIf: S('mems'), help: 'Relative permittivity of the insulating layer (silicon nitride ≈ 7.5).' },
  { key: 'V_m', label: 'Applied voltage', unit: 'V', value: 5, min: 0, max: 500, group: 'MEMS actuator', showIf: S('mems'), help: 'Bias voltage across the gap.' },
  { key: 'Q_m', label: 'Mechanical quality factor', value: 5, min: 0.2, max: 1e5, log: true, group: 'MEMS actuator', showIf: S('mems'), help: 'Quality factor including squeeze-film damping (low in air, high in vacuum).' },
  // ferrite
  { key: 'Ms', label: 'Saturation magnetisation', unit: 'A/m', value: 1.4e5, min: 1e3, max: 2e6, log: true, group: 'Ferrite', bus: 'prop.Ms', showIf: S('ferrite'), help: 'Saturation magnetisation (YIG ≈ 1.4×10⁵ A/m, i.e. 4πM_s ≈ 1750 G).' },
  { key: 'H_b', label: 'Internal bias field', unit: 'A/m', value: 1e5, min: 1e3, max: 2e6, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Static internal magnetic field (1 Oe = 79.6 A/m).' },
  { key: 'alpha_g', label: 'Gilbert damping', value: 0.003, min: 1e-5, max: 0.2, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Dimensionless damping (single-crystal YIG ≈ 10⁻⁴, polycrystalline ferrites 10⁻³–10⁻²).' },
  { key: 'f_op', label: 'Operating frequency', unit: 'Hz', value: 6e9, min: 1e7, max: 1e11, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Frequency at which the effective permeability is evaluated.' },
  { key: 'N_coil', label: 'Coil turns', value: 200, min: 1, max: 20000, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Turns of the bias electromagnet.' },
  { key: 'gap_mm', label: 'Magnet gap', unit: 'mm', value: 5, min: 0.1, max: 200, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Air gap of the bias magnet (core reluctance neglected).' },
  { key: 'A_pole', label: 'Pole area', unit: 'cm²', value: 4, min: 0.01, max: 1000, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'Pole-face area, used for the coil inductance.' },
  { key: 'R_coil', label: 'Coil resistance', unit: 'Ω', value: 2, min: 0.01, max: 1000, log: true, group: 'Ferrite', showIf: S('ferrite'), help: 'DC resistance of the bias coil.' },
  // bistable
  { key: 'n_units', label: 'Number of bistable units', value: 8, min: 1, max: 40, step: 1, group: 'Bistable array', showIf: S('bistable'), help: 'Units stacked in series; each can be in one of two stable states.' },
  { key: 'k_bar', label: 'Bar axial stiffness', unit: 'N/m', value: 2e4, min: 1, max: 1e8, log: true, group: 'Bistable array', showIf: S('bistable'), help: 'Axial stiffness of each inclined bar of the von Mises truss unit.' },
  { key: 'w_b', label: 'Half-span', unit: 'mm', value: 10, min: 0.001, max: 1000, log: true, group: 'Bistable array', showIf: S('bistable'), help: 'Horizontal half-width of the truss unit.' },
  { key: 'h_b', label: 'Initial apex height', unit: 'mm', value: 1.2, min: 0.0001, max: 500, log: true, group: 'Bistable array', showIf: S('bistable'), help: 'Rise of the apex above the supports; snap-through travel is twice this height.' },
  { key: 'k_v', label: 'Parallel linear spring', unit: 'N/m', value: 0, min: 0, max: 1e7, group: 'Bistable array', showIf: S('bistable'), help: 'Linear spring acting on the apex; it biases the unit and removes bistability when large.' },
  { key: 'var_pct', label: 'Unit-to-unit variation', unit: '%', value: 5, min: 0, max: 30, step: 0.5, group: 'Bistable array', showIf: S('bistable'), help: 'Manufacturing scatter of unit stiffness; it fixes the order in which units snap.' },
  { key: 'm_unit', label: 'Moving mass per unit', unit: 'kg', value: 1e-3, min: 1e-12, max: 100, log: true, group: 'Bistable array', showIf: S('bistable'), help: 'Mass moved during snap-through; sets the switching time.' },
  // QZS
  { key: 'kq_v', label: 'Vertical spring stiffness', unit: 'N/m', value: 2000, min: 1, max: 1e7, log: true, group: 'QZS isolator', showIf: S('qzs'), help: 'Load-bearing positive-stiffness spring.' },
  { key: 'kq_h', label: 'Oblique spring stiffness', unit: 'N/m', value: 1500, min: 1, max: 1e7, log: true, group: 'QZS isolator', showIf: S('qzs'), help: 'Stiffness of each of the two pre-compressed lateral springs supplying negative stiffness.' },
  { key: 'L0_q', label: 'Free length of oblique springs', unit: 'mm', value: 100, min: 1, max: 2000, group: 'QZS isolator', showIf: S('qzs'), help: 'Unloaded length of the lateral springs.' },
  { key: 'a_q', label: 'Compressed length at equilibrium', unit: 'mm', value: 62, min: 0.5, max: 2000, group: 'QZS isolator', showIf: S('qzs'), help: 'Length of the lateral springs when horizontal (the tuning parameter): shorter means more pre-compression and more negative stiffness.' },
  { key: 'm_q', label: 'Payload mass', unit: 'kg', value: 5, min: 0.001, max: 1e4, log: true, group: 'QZS isolator', showIf: S('qzs'), help: 'Isolated mass.' },
  { key: 'zeta_q', label: 'Damping ratio (of the linear mount)', value: 0.02, min: 0.001, max: 0.5, log: true, group: 'QZS isolator', showIf: S('qzs'), help: 'Damping ratio referred to the vertical spring alone; the damper itself is unchanged by tuning.' },
  { key: 'f_exc', label: 'Excitation frequency', unit: 'Hz', value: 2, min: 0.01, max: 1000, log: true, group: 'QZS isolator', showIf: S('qzs'), help: 'Disturbance frequency to be isolated.' },
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
  if (!Object.keys(r.publish).length) r.publish['ctrl.status'] = 'not evaluated';
  return r;
}
const znoEqs = () => { try { const z = chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NaOH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'H2O', phase: 'l' }, { f: 'NaNO3', phase: 'aq' }]); return [eqn('chemical', 'Controlled reaction (ZnO hydrothermal growth)', z.molecular), eqn('ionic', 'Net ionic equation', z.netIonic, 'spectators: ' + z.spectators.join(', '))]; } catch { return []; } };

// ---------- PID ----------
function pid(p) {
  const warnings = [], C = 4.186e6 * p.V_reac + p.C_vessel, tau = C / p.UA, K = p.P_heater / p.UA, th = p.dead, step = p.T_set - p.T_amb;
  let Kc, Ti, Td, rule;
  if (p.tuning === 'zn') { Kc = 1.2 * tau / (K * th); Ti = 2 * th; Td = 0.5 * th; rule = 'Ziegler–Nichols open-loop'; }
  else if (p.tuning === 'cc') { Kc = tau / (K * th) * (4 / 3 + th / (4 * tau)); Ti = th * (32 + 6 * th / tau) / (13 + 8 * th / tau); Td = 4 * th / (11 + 2 * th / tau); rule = 'Cohen–Coon'; }
  else if (p.tuning === 'imc') { const tc = p.lam * th; Kc = (tau + th / 2) / (K * (tc + th / 2)); Ti = tau + th / 2; Td = tau * th / (2 * tau + th); rule = 'IMC PID'; }
  else if (p.tuning === 'manual') { Kc = p.Kc_man; Ti = p.Ti_man; Td = p.Td_man; rule = 'manual'; }
  else { const tc = p.lam * th; Kc = tau / (K * (tc + th)); Ti = Math.min(tau, 4 * (tc + th)); Td = 0; rule = 'SIMC PI'; }
  let Ts = p.Ts, n = Math.round(p.t_sim / Ts); if (n > 20000) { n = 20000; Ts = p.t_sim / n; warnings.push(`Sampling interval increased to ${F(Ts, 3)} s to limit the simulation to 20 000 steps.`); } if (n < 20) { n = 20; Ts = p.t_sim / n; }
  if (Ts > th / 2) warnings.push(`Sampling interval (${F(Ts, 3)} s) is long compared with the dead time (${F(th, 3)} s): the sampling delay erodes the phase margin.`);
  const d = Math.max(0, Math.round(th / Ts)), a = Math.exp(-Ts / tau), Nf = 10, rnd = N.rng(p.seed), kD = Math.round(0.6 * n), dis = p.dist / 100;
  const T = new Float64Array(n + 1), U = new Float64Array(n + 1), t = new Float64Array(n + 1), ubuf = new Float64Array(d + 1); let x = 0, I = 0, Dt = 0, yPrev = p.T_amb, satN = 0, iae = 0, ise = 0, en = 0, bi = 0;
  T[0] = p.T_amb;
  for (let k = 0; k < n; k++) {
    const y = p.T_amb + x + p.noise * rnd.normal(), e = p.T_set - y;
    if (Td > 0) Dt = Td / (Td + Nf * Ts) * Dt - Kc * Td * Nf / (Td + Nf * Ts) * (y - yPrev); yPrev = y;
    const v = Kc * e + I + Dt, u = N.clamp(v, 0, 1); if (u !== v) satN++;
    I += Kc * Ts / Ti * e + (p.aw === 'on' ? Ts / Math.max(Ti, Ts) * (u - v) : 0);
    ubuf[bi] = u; bi = (bi + 1) % (d + 1); const ud = ubuf[bi]; // value written d steps ago
    x = a * x + (1 - a) * K * (ud - (k >= kD ? dis : 0));
    U[k] = u; T[k + 1] = p.T_amb + x; t[k + 1] = (k + 1) * Ts; const et = p.T_set - T[k + 1]; iae += Math.abs(et) * Ts; ise += et * et * Ts; en += u * p.P_heater * Ts;
  }
  U[n] = U[n - 1];
  const band = 0.02 * Math.abs(step) || 0.5; let tr10 = NaN, tr90 = NaN, ovs = 0, tSet = 0, devD = 0, tRec = 0;
  for (let k = 0; k <= kD; k++) { const fr = step !== 0 ? (T[k] - p.T_amb) / step : 1; if (!Number.isFinite(tr10) && fr >= 0.1) tr10 = t[k]; if (!Number.isFinite(tr90) && fr >= 0.9) tr90 = t[k]; ovs = Math.max(ovs, step !== 0 ? (T[k] - p.T_set) / step * Math.sign(step) : 0); if (Math.abs(T[k] - p.T_set) > band) tSet = t[Math.min(k + 1, kD)]; }
  for (let k = kD; k <= n; k++) { devD = Math.max(devD, Math.abs(T[k] - p.T_set)); if (Math.abs(T[k] - p.T_set) > band) tRec = t[Math.min(k + 1, n)] - t[kD]; }
  const reached = Number.isFinite(tr90), settled = tSet < t[kD] - 5 * Ts, rise = reached ? tr90 - tr10 : t[kD], endErr = Math.abs(T[n] - p.T_set), feasible = K * (1 - dis) >= step;
  // frequency response of the loop (continuous controller + plant + half-sample delay)
  const thE = th + Ts / 2, Lw = w => { const wd = w * Td, den = 1 + (wd / Nf) ** 2, cr = Kc * (1 + wd * wd / Nf / den), ci = Kc * (-1 / (w * Ti) + wd / den), g = K / Math.sqrt(1 + (w * tau) ** 2); return { mag: Math.hypot(cr, ci) * g, ph: Math.atan2(ci, cr) - Math.atan(w * tau) - w * thE }; };
  const ws = N.logspace(1e-3 * Math.min(1 / tau, 1 / Ti), N.PI / Ts, 500), LL = ws.map(Lw); let wgc = 0, wpc = 0, Ms = 0;
  for (let i = 1; i < ws.length; i++) { if (!wgc && LL[i - 1].mag >= 1 && LL[i].mag < 1) wgc = Math.exp(N.brent(lw => Math.log(Lw(Math.exp(lw)).mag), Math.log(ws[i - 1]), Math.log(ws[i]))); if (!wpc && LL[i - 1].ph >= -N.PI && LL[i].ph < -N.PI) wpc = Math.exp(N.brent(lw => Lw(Math.exp(lw)).ph + N.PI, Math.log(ws[i - 1]), Math.log(ws[i]))); }
  for (const l of LL) Ms = Math.max(Ms, 1 / Math.sqrt(1 + 2 * l.mag * Math.cos(l.ph) + l.mag * l.mag));
  const PM = wgc ? (N.PI + Lw(wgc).ph) * 180 / N.PI : 180, GM = wpc ? 1 / Lw(wpc).mag : 1000, GMdB = 20 * Math.log10(GM), delayM = wgc ? PM * N.PI / 180 / wgc : 1e9, stable = GM > 1 && PM > 0;
  const w180 = N.brent(w => -Math.atan(w * tau) - w * th + N.PI, 1e-9 / th, N.PI / th), Ku = Math.sqrt(1 + (w180 * tau) ** 2) / K, Pu = 2 * N.PI / w180;
  const Phold = step / K * p.P_heater, eKWh = en / 3.6e6;
  if (!feasible) warnings.push(`Heater too small: full power raises the temperature by only ${F(K * (1 - dis), 3)} K against the ${F(step, 3)} K required (including the disturbance).`);
  if (!reached) warnings.push('The set-point was not reached (90 % of the step) before the disturbance; extend the simulated time or increase heater power.');
  if (!stable) warnings.push('The loop is unstable according to the frequency-domain margins: reduce the gain or increase the integral time.');
  if (p.aw === 'off' && satN > 0.05 * n) warnings.push('The heater saturated without anti-windup: the integrator accumulates during heat-up, which produces the overshoot seen in the response.');
  if (Td > 0 && p.noise > 0) warnings.push('Derivative action amplifies sensor noise; the derivative is filtered (N = 10) and acts on the measurement only.');
  return {
    summary: `${rule} control of a ${F(p.V_reac * 1e3, 3)} L reactor (FOPDT plant: gain ${F(K, 4)} K per unit duty, time constant ${F(tau, 4)} s, dead time ${F(th, 3)} s) with K_c = ${F(Kc, 3)} 1/K, T_i = ${F(Ti, 4)} s, T_d = ${F(Td, 3)} s: the set-point ${F(p.T_set)} °C is ${reached ? `reached with ${F(100 * Math.max(ovs, 0), 3)} % overshoot, rise time ${F(rise, 4)} s and settling after ${F(tSet, 4)} s` : 'NOT reached in the simulated time'}. Gain margin ${F(GMdB, 3)} dB, phase margin ${F(PM, 3)}°, peak sensitivity ${F(Ms, 3)}; the ${F(p.dist)} % load step causes a ${F(devD, 3)} K deviation. Heating energy ${F(eKWh, 3)} kWh.`,
    kpis: [
      kpi('overshoot', 'Overshoot', 100 * Math.max(ovs, 0), '%', 'Largest excursion beyond the set-point before the disturbance, as a percentage of the set-point step. Overshoot matters in synthesis because it can trigger secondary nucleation.', ovs < 0.05),
      kpi('t_rise', 'Rise time (10–90 %)', rise, 's', `Time for the temperature to go from 10 % to 90 % of the step. During heat-up the heater is saturated for ${F(100 * satN / n, 3)} % of the samples, so the rise is limited by installed power, not tuning.`, reached),
      kpi('t_settle', 'Settling time (2 %)', tSet, 's', 'Time after which the temperature stays within ±2 % of the step around the set-point (evaluated before the disturbance).', settled),
      kpi('iae', 'Integral absolute error', iae, 'K·s', 'Time integral of |error| over the whole run, including heat-up and disturbance rejection; lower is better for a fixed scenario.', null),
      kpi('ise', 'Integral squared error', ise, 'K²·s', 'Time integral of the squared error; it weights large deviations (heat-up) more heavily than IAE.', null),
      kpi('GM', 'Gain margin', GMdB, 'dB', `Factor by which the loop gain could rise before instability, from the loop frequency response including a half-sample delay. Ultimate gain of the plant alone: ${F(Ku, 3)} 1/K at a period of ${F(Pu, 4)} s.`, GMdB >= 6),
      kpi('PM', 'Phase margin', PM, '°', `Additional phase lag tolerated at the gain crossover (${F(wgc, 3)} rad/s); equivalent to a delay margin of ${F(delayM, 3)} s.`, PM >= 30),
      kpi('Ms', 'Peak sensitivity M_s', Ms, '–', 'Maximum of |1/(1 + L)|: the worst-case amplification of disturbances. Values of 1.2–2 indicate a robust loop.', Ms <= 2),
      kpi('dev_dist', 'Peak deviation after load step', devD, 'K', `Largest temperature error after the ${F(p.dist)} % load disturbance; recovery to the 2 % band takes ${F(tRec, 4)} s.`, devD < 0.1 * Math.abs(step) + 0.5),
      kpi('energy', 'Heating energy', eKWh, 'kWh', `Integral of heater power over the run; steady holding needs ${F(Phold, 4)} W (${F(100 * step / K, 3)} % duty).`, null),
      kpi('sat', 'Actuator saturation', 100 * satN / n, '%', 'Fraction of samples in which the controller output was clipped at 0 or 100 % heater power.', null),
      kpi('Kc', 'Proportional gain', Kc, '1/K', `Gain from the ${rule} rule, in heater fraction per kelvin; integral time ${F(Ti, 4)} s and derivative time ${F(Td, 3)} s.`, null),
    ],
    series: [
      { key: 'T_t', title: 'Closed-loop temperature response', xlabel: 'time (s)', ylabel: 'temperature (°C)', lines: [line('temperature', Array.from(t), Array.from(T)), line('set-point', [0, t[n]], [p.T_set, p.T_set])], explain: 'Heat-up from ambient followed by the load disturbance at 60 % of the run. Look for overshoot after the saturated heat-up phase and the size and duration of the dip after the disturbance.' },
      { key: 'u_t', title: 'Heater duty', xlabel: 'time (s)', ylabel: 'heater output (%)', lines: [line('duty', Array.from(t), Array.from(U, v => 100 * v))], explain: 'Controller output. Full saturation during heat-up is normal; rapid chatter indicates noise amplification by excessive gain or derivative action.' },
      { key: 'bode', title: 'Loop gain (Bode magnitude)', xlabel: 'frequency (rad/s)', ylabel: '|L(jω)|', logx: true, logy: true, lines: [line('|L|', ws, LL.map(l => Math.max(l.mag, 1e-12))), line('unity', [ws[0], ws[ws.length - 1]], [1, 1])], explain: 'Open-loop gain of controller times plant. The crossing of unity defines the bandwidth; the gain at the −180° frequency defines the gain margin.' },
      { key: 'phase', title: 'Loop phase', xlabel: 'frequency (rad/s)', ylabel: 'phase (°)', logx: true, lines: [line('∠L', ws, LL.map(l => Math.max(l.ph * 180 / N.PI, -720))), line('−180°', [ws[0], ws[ws.length - 1]], [-180, -180])], explain: 'Phase of the loop transfer function; the distance above −180° at the gain crossover is the phase margin. Dead time makes the phase fall without limit.' },
    ],
    equations: [
      eqn('governing', 'Plant (first order plus dead time)', 'C dT/dt = P·u(t − θ) − UA·(T − T_amb) − d(t);   G(s) = K e^(−θs)/(τs + 1),  K = P/UA, τ = C/UA'),
      eqn('constitutive', 'PID law', 'u = K_c[e + (1/T_i)∫e dt − T_d dy_f/dt],  0 ≤ u ≤ 1,  derivative filter N = 10'),
      eqn('formula', 'Anti-windup (back-calculation)', 'dI/dt = (K_c/T_i)e + (u_sat − u)/T_t'),
      eqn('formula', 'Tuning rule used', p.tuning === 'zn' ? 'K_c = 1.2τ/(Kθ), T_i = 2θ, T_d = 0.5θ' : p.tuning === 'cc' ? 'K_c = (τ/Kθ)(4/3 + θ/4τ), T_i = θ(32 + 6θ/τ)/(13 + 8θ/τ), T_d = 4θ/(11 + 2θ/τ)' : p.tuning === 'imc' ? 'K_c = (τ + θ/2)/[K(τ_c + θ/2)], T_i = τ + θ/2, T_d = τθ/(2τ + θ)' : p.tuning === 'manual' ? 'gains entered manually' : 'K_c = τ/[K(τ_c + θ)], T_i = min(τ, 4(τ_c + θ))'),
      eqn('formula', 'Stability margins', 'GM = 1/|L(jω₁₈₀)|;  PM = 180° + ∠L(jω_c);  M_s = max|1/(1 + L)|'),
      ...znoEqs(),
    ],
    tables: [mechTable([['Physical mechanism', 'Resistive heating against conductive and convective loss; first-order thermal lag with transport delay'], ['Energy demand', `${F(eKWh, 3)} kWh for this run; ${F(Phold, 4)} W to hold the set-point`], ['Response time', `rise ${F(rise, 4)} s (power-limited), disturbance recovery ${F(tRec, 4)} s`], ['Loss', `${F(100 * step / K, 3)} % of heater rating is continuously lost to ambient at the set-point`], ['Reversibility', 'Heating only: cooling relies on passive loss (time constant ' + F(tau, 4) + ' s), so overshoot is slow to correct'], ['Hysteresis', 'None in the model; real heaters with relay drive add a dead band'], ['Cycle durability', 'Solid-state relay switching at the sampling rate; heater life limited by thermal cycling of the element']])],
    checks: [
      chk('Closed loop stable (margins)', stable, `GM = ${F(GMdB, 3)} dB, PM = ${F(PM, 3)}°`),
      chk('Robustness guideline GM ≥ 6 dB, PM ≥ 30°', GMdB >= 6 && PM >= 30, `M_s = ${F(Ms, 3)}`),
      chk('Set-point reached and held', reached && endErr < Math.max(band, 3 * p.noise + 0.2), `final error ${F(endErr, 3)} K`),
      chk('Heater able to hold set-point under load', feasible, `available rise ${F(K * (1 - dis), 4)} K vs required ${F(step, 4)} K`),
      chk('Sampling fast enough', Ts <= th / 2, `T_s/θ = ${F(Ts / th, 3)}`),
      chk('Temperature bounded', Array.from(T).every(v => v > p.T_amb - 5 && v < p.T_amb + K + 5), `range ${F(N.minOf(T), 4)}–${F(N.maxOf(T), 4)} °C`),
    ],
    warnings,
    recommendations: [
      rec('high', !stable || PM < 30 || GMdB < 6 ? `Detune the loop: margins (${F(GMdB, 3)} dB, ${F(PM, 3)}°) are below the 6 dB / 30° guideline. Use SIMC with τ_c ≥ θ, which gives about 60° phase margin for this plant.` : `Margins are adequate (${F(GMdB, 3)} dB, ${F(PM, 3)}°); the loop tolerates a dead-time increase of ${F(delayM, 3)} s or a gain change by a factor ${F(GM, 3)}.`, 'Frequency-domain stability margins.'),
      rec('high', ovs > 0.05 ? `Overshoot of ${F(100 * ovs, 3)} % risks secondary nucleation and morphology changes: enable anti-windup, ramp the set-point, or lengthen the integral time.` : 'Overshoot is below 5 %; keep anti-windup enabled because heat-up is always power-limited.', 'Step-response overshoot.'),
      rec('medium', `Identify K, τ and θ from an open-loop step test on the real reactor (here ${F(K, 3)} K, ${F(tau, 3)} s, ${F(th, 3)} s were derived from thermal mass and loss) and re-run the tuning; dead time is the parameter that most affects stability.`, 'Model-based tuning rules.'),
      rec('medium', `Holding ${F(p.T_set)} °C costs ${F(Phold, 3)} W continuously; improving insulation (halving UA) halves that loss and doubles the plant gain, so re-tune afterwards. Avoiding overshoot also avoids wasted heat and off-specification batches.`, 'Steady-state energy balance.', true),
    ],
    publish: { 'ctrl.settling_s': tSet, 'ctrl.overshoot': Math.max(ovs, 0), 'ctrl.gain_margin_dB': GMdB, 'ctrl.phase_margin_deg': PM, 'proc.T_C': p.T_set, 'proc.hold_power': Phold, 'ctrl.plant_tau': tau, 'ctrl.plant_delay': th, 'ctrl.plant_gain': K / p.P_heater },
    numerics: { method: 'Exact discretisation of the FOPDT plant with sampled PID (derivative filter, back-calculation anti-windup); margins by Brent root finding on the analytical loop response', dof: 1, steps: n, iterations: n, residual: endErr, notes: `Sampling ${F(Ts, 3)} s, delay ${d} samples, seeded sensor noise.` },
  };
}

// ---------- MPC of growth ----------
function mpc(p) {
  const warnings = [], n = Math.round(p.n_ctrl), Nh = Math.round(p.N_h), dt = p.t_growth / n, Lt = p.L_target, Tref = 363.15, R = N.R_GAS, Ea = p.Ea * 1e3, temp = p.manip === 'temp', rnd = N.rng(p.seed);
  let Tlo = p.T_lo + 273.15, Thi = p.T_hi + 273.15; if (Thi < Tlo + 2) { Thi = Tlo + 2; warnings.push('Maximum temperature raised to 2 K above the minimum.'); }
  const vOfZ = z => temp ? p.G_ref * Math.exp(-Ea / R * (1 / z - 1 / Tref)) : p.G_ref * z, zOfV = v => temp ? 1 / (1 / Tref - R * Math.log(Math.max(v, 1e-300) / p.G_ref) / Ea) : v / p.G_ref;
  const zLo = temp ? Tlo : 0, zHi = temp ? Thi : 2.5, vLo = vOfZ(zLo), vHi = vOfZ(zHi), vNom = Lt / p.t_growth, feasible = vNom >= vLo && vNom <= vHi, nuLo = vLo / vNom, nuHi = vHi / vNom, g = 1 / n;
  const zNom = N.clamp(zOfV(N.clamp(vNom, Math.max(vLo, 1e-300), vHi)), zLo, zHi), mult = 1 + p.bias / 100, sub = 10, aAct = 1 - Math.exp(-dt / sub / p.tau_act);
  const plant = (L, z, zc) => { for (let s = 0; s < sub; s++) { z += (zc - z) * aAct; L += dt / sub * vOfZ(z) * mult * (1 - p.deplete * Math.min(L / Lt, 1.2)); } return [L, z]; };
  const ref = k => Lt * Math.min(k / n, 1);
  // open loop: constant nominal input
  const Lol = [0]; { let L = 0, z = zNom; for (let k = 0; k < n; k++) { [L, z] = plant(L, z, zNom); Lol.push(L); } }
  // closed loop
  const Lcl = [0], Zc = [], Za = [], Mh = []; let L = 0, z = zNom, mhat = 1, yPrev = 0, nuPrev = N.clamp(1, nuLo, nuHi), plan = new Float64Array(Nh).fill(nuPrev), itTot = 0, satN = 0, pgMax = 0, zPrevAct = z;
  for (let k = 0; k < n; k++) {
    const y = L + p.len_noise * 1e-9 * rnd.normal();
    if (k > 0) { const pred = vOfZ(0.5 * (zPrevAct + z)) * dt; if (pred > 1e-30) mhat = N.clamp(mhat + 0.3 * ((y - yPrev) / pred - mhat), 0.2, 5); }
    yPrev = y; zPrevAct = z; const nLive = Math.min(Nh, n - k), l0 = y / Lt, gm = g * mhat, rr = Array.from({ length: Nh }, (_, j) => ref(k + j + 1) / Lt), Lip = 2 * (gm * gm * Nh * (Nh + 1) / 2 + 4 * p.rho_u * g * g) + 1e-30;
    for (let j = 0; j < Nh - 1; j++) plan[j] = plan[j + 1]; // warm start by shifting
    let xk = Float64Array.from(plan), yk = Float64Array.from(plan), tk = 1; const grad = new Float64Array(Nh), cost = v => { let l = l0, J = 0, pr = nuPrev; for (let j = 0; j < nLive; j++) { l += gm * v[j]; J += (l - rr[j]) ** 2 + p.rho_u * g * g * (v[j] - pr) ** 2; pr = v[j]; } return J; };
    const gradient = v => { let l = l0, acc = 0; const e = new Float64Array(Nh); for (let j = 0; j < nLive; j++) { l += gm * v[j]; e[j] = l - rr[j]; } for (let j = nLive - 1; j >= 0; j--) { acc += e[j]; grad[j] = 2 * gm * acc + 2 * p.rho_u * g * g * ((v[j] - (j > 0 ? v[j - 1] : nuPrev)) - (j < nLive - 1 ? v[j + 1] - v[j] : 0)); } };
    let it = 0, pg = 1;
    for (; it < 80 && pg > 1e-9; it++) { gradient(yk); pg = 0; const xn = Float64Array.from(xk); for (let j = 0; j < nLive; j++) { xn[j] = N.clamp(yk[j] - grad[j] / Lip, nuLo, nuHi); pg = Math.max(pg, Math.abs(xn[j] - yk[j])); } const tn = 0.5 * (1 + Math.sqrt(1 + 4 * tk * tk)); for (let j = 0; j < nLive; j++) yk[j] = xn[j] + (tk - 1) / tn * (xn[j] - xk[j]); if (cost(xn) > cost(xk)) { yk = Float64Array.from(xn); tk = 1; } else tk = tn; xk = xn; }
    itTot += it; pgMax = Math.max(pgMax, pg); plan = xk; const nu = plan[0]; if (nu <= nuLo + 1e-9 || nu >= nuHi - 1e-9) satN++;
    const zc = N.clamp(zOfV(nu * vNom), zLo, zHi); nuPrev = nu; Zc.push(zc); Mh.push(mhat);
    [L, z] = plant(L, z, zc); Za.push(z); Lcl.push(L);
  }
  const tt = N.linspace(0, p.t_growth, n + 1), refs = tt.map((_, k) => ref(k)), rm = arr => Math.sqrt(N.mean(arr.map((v, k) => (v - refs[k]) ** 2)));
  const rmCl = rm(Lcl), rmOl = rm(Lol), eCl = Lcl[n] - Lt, eOl = Lol[n] - Lt, zDisp = v => temp ? v - 273.15 : v, dz = Zc.map((v, k) => k ? Math.abs(v - Zc[k - 1]) : 0), maxMove = N.maxOf(dz), zMean = N.mean(Za), improve = Math.abs(eOl) > 1e-30 ? Math.abs(eCl) / Math.abs(eOl) : 0;
  const Eheat = temp ? N.sum(Za.map(v => v - 298.15)) * dt : 0, EheatOl = temp ? (zNom - 298.15) * p.t_growth : 0;
  if (!feasible) warnings.push(`The nominal growth rate ${F(vNom, 3)} m/s needed for the target lies outside the range achievable within the input constraints (${F(vLo, 3)}–${F(vHi, 3)} m/s): the target cannot be met.`);
  if (satN > 0.3 * n) warnings.push(`The manipulated variable sat on a constraint for ${F(100 * satN / n, 3)} % of the run: the mismatch exceeds the available control authority.`);
  if (p.tau_act > 3 * dt) warnings.push('The actuator lag is longer than three control intervals and is not in the controller model; expect sluggish correction. Use longer intervals or add the lag to the model.');
  const unit = temp ? '°C' : 'relative concentration';
  return {
    summary: `Model-predictive control of nanorod growth by ${temp ? 'temperature' : 'precursor dosing'} (${n} intervals of ${F(dt, 4)} s, horizon ${Nh}): final length ${F(Lcl[n] * 1e9, 5)} nm against a target of ${F(Lt * 1e9, 5)} nm (error ${F(eCl * 1e9, 3)} nm, ${F(100 * eCl / Lt, 3)} %), compared with ${F(100 * eOl / Lt, 3)} % in open loop under a ${F(p.bias)} % rate mismatch and ${F(100 * p.deplete, 3)} % depletion. Tracking RMSE ${F(rmCl * 1e9, 4)} nm (open loop ${F(rmOl * 1e9, 4)} nm); the controller identified a rate multiplier of ${F(Mh[n - 1], 3)}.`,
    kpis: [
      kpi('L_final', 'Final length (closed loop)', Lcl[n], 'm', 'Length at the end of the run with the MPC acting on noisy length measurements.', Math.abs(eCl) < 0.05 * Lt),
      kpi('err_cl', 'Final length error (closed loop)', 100 * eCl / Lt, '%', 'Relative deviation of the final length from the target with feedback.', Math.abs(eCl) < 0.05 * Lt),
      kpi('err_ol', 'Final length error (open loop)', 100 * eOl / Lt, '%', 'Deviation when the nominal recipe (constant input computed from the model) is applied without feedback: the cost of plant–model mismatch and depletion.', Math.abs(eOl) < 0.05 * Lt),
      kpi('rmse_cl', 'Tracking RMSE (closed loop)', rmCl, 'm', 'Root-mean-square deviation from the linear reference trajectory over the run.', rmCl < 0.05 * Lt),
      kpi('rmse_ol', 'Tracking RMSE (open loop)', rmOl, 'm', 'Root-mean-square deviation from the reference when the constant nominal recipe is applied without feedback.', null),
      kpi('improve', 'Residual error ratio', improve, '–', 'Closed-loop final error divided by the open-loop final error; values well below 1 show the benefit of feedback.', improve < 0.5),
      kpi('m_hat', 'Identified rate multiplier', Mh[n - 1], '–', `Controller estimate of true-to-model growth-rate ratio from measured length increments (true initial value ${F(mult, 3)}, falling with depletion). It provides offset-free prediction.`, null),
      kpi('z_mean', temp ? 'Mean reactor temperature' : 'Mean relative concentration', zDisp(zMean), unit, `Average of the actual manipulated variable; the open-loop recipe uses ${F(zDisp(zNom), 4)} ${unit}.`, null),
      kpi('move', 'Largest control move', maxMove, temp ? 'K' : '–', 'Largest change of the commanded input between consecutive intervals; limited by the move-suppression weight.', null),
      kpi('sat', 'Time on constraints', 100 * satN / n, '%', 'Share of control intervals in which the optimal input lay on a bound.', satN < 0.3 * n),
      kpi('qp_it', 'Mean optimiser iterations', itTot / n, '–', `Average accelerated projected-gradient iterations per control step; largest final step ${pgMax.toExponential(1)} (normalised).`, pgMax < 1e-4),
    ],
    series: [
      { key: 'L_t', title: 'Nanorod length trajectory', xlabel: 'time (h)', ylabel: 'length (nm)', lines: [line('reference', tt.map(v => v / 3600), refs.map(v => v * 1e9)), line('MPC (closed loop)', tt.map(v => v / 3600), Lcl.map(v => v * 1e9)), line('open-loop recipe', tt.map(v => v / 3600), Lol.map(v => v * 1e9))], explain: 'Length against time. The open-loop run drifts away from the reference because the true rate differs from the model and falls as precursor is consumed; the MPC corrects the drift.' },
      { key: 'u_t', title: temp ? 'Reactor temperature' : 'Relative precursor concentration', xlabel: 'time (h)', ylabel: unit, lines: [line('command', tt.slice(0, n).map(v => v / 3600), Zc.map(zDisp)), line('actual', tt.slice(1).map(v => v / 3600), Za.map(zDisp)), line('lower limit', [0, p.t_growth / 3600], [zDisp(zLo), zDisp(zLo)]), line('upper limit', [0, p.t_growth / 3600], [zDisp(zHi), zDisp(zHi)])], explain: 'Manipulated variable chosen by the optimiser and the lagged value actually reached, with the input constraints. A command resting on a limit shows lack of authority.' },
      { key: 'm_t', title: 'Identified growth-rate multiplier', xlabel: 'time (h)', ylabel: 'estimated true/model rate', lines: [line('estimate', tt.slice(0, n).map(v => v / 3600), Mh), line('true', tt.slice(0, n).map(v => v / 3600), Lcl.slice(0, n).map(v => mult * (1 - p.deplete * Math.min(v / Lt, 1.2))))], explain: 'Online estimate of the ratio between real and modelled growth rate, used as a disturbance model. It should follow the true multiplier, which falls as precursor depletes.' },
    ],
    equations: [
      eqn('governing', 'Growth model used by the controller', temp ? 'dL/dt = m̂·G_ref·exp[−E_a/R·(1/T − 1/T_ref)]' : 'dL/dt = m̂·G_ref·(c/c_ref)'),
      eqn('governing', 'Finite-horizon optimal control problem', 'min Σⱼ (L_{k+j} − r_{k+j})² + ρ Σⱼ (Δv_{k+j})²   s.t.  v_min ≤ v ≤ v_max,  L_{j+1} = L_j + Δt·m̂·v_j;  the horizon shrinks to the intervals remaining in the batch'),
      eqn('formula', 'Disturbance estimator', 'm̂ ← m̂ + 0.3·[ΔL_measured/(v_model·Δt) − m̂]'),
      eqn('formula', 'Solver', 'accelerated projected gradient on the box-constrained QP, step 1/L_Lipschitz, warm start'),
      ...znoEqs(),
    ],
    tables: [mechTable([['Physical mechanism', temp ? 'Arrhenius temperature dependence of the surface-reaction-limited growth rate' : 'First-order dependence of growth rate on precursor concentration (dosing pump)'], ['Energy demand', temp ? `integrated excess temperature ${F(Eheat / 3600, 4)} K·h (open loop ${F(EheatOl / 3600, 4)} K·h); heater energy is this times the loss conductance` : 'pump power only; reagent consumption follows the concentration profile'], ['Response time', `actuator lag ${F(p.tau_act, 3)} s; control interval ${F(dt, 3)} s`], ['Loss', 'Reagent consumed by homogeneous nucleation at high supersaturation is not available for rod growth'], ['Reversibility', 'Growth is irreversible: excess length cannot be removed, so overshoot must be avoided by constraints'], ['Hysteresis', 'None in the rate law; dissolution–regrowth asymmetry is not modelled'], ['Cycle durability', 'Not applicable (single batch); sensor fouling limits repeated in-situ measurement']])],
    checks: [
      chk('Target reachable within input constraints', feasible, `required rate ${F(vNom, 3)} m/s, achievable ${F(vLo, 3)}–${F(vHi, 3)} m/s`),
      chk('Constraints respected', Zc.every(v => v >= zLo - 1e-9 && v <= zHi + 1e-9), `command range ${F(zDisp(N.minOf(Zc)), 4)}–${F(zDisp(N.maxOf(Zc)), 4)} ${unit}`),
      chk('Optimiser converged', pgMax < 1e-4, `largest final projected step ${pgMax.toExponential(2)}`),
      chk('Feedback improves on open loop', Math.abs(eCl) <= Math.abs(eOl) + 3 * p.len_noise * 1e-9, `|error| ${F(Math.abs(eCl) * 1e9, 3)} nm vs ${F(Math.abs(eOl) * 1e9, 3)} nm`),
      chk('Length monotonic (growth irreversible)', Lcl.every((v, k) => k === 0 || v >= Lcl[k - 1]), 'no decrease of length in the simulated plant'),
    ],
    warnings,
    recommendations: [
      rec('high', Math.abs(eOl) > 0.05 * Lt ? `Open-loop recipes miss the target by ${F(100 * eOl / Lt, 3)} % under this mismatch; feedback from an in-situ length or mass signal reduces the error to ${F(100 * eCl / Lt, 3)} %.` : 'The open-loop recipe is already within 5 %; feedback mainly protects against batch-to-batch variation.', 'Closed- versus open-loop final error.'),
      rec('medium', satN > 0.3 * n ? 'Widen the admissible temperature or concentration window, or lengthen the growth time: the controller is constraint-limited.' : `Control authority is sufficient (${F(100 * satN / n, 3)} % of time on limits); the constraints could be tightened to protect morphology.`, 'Constraint activity.'),
      rec('medium', `Re-identify G_ref and E_a from the logged data after each batch: the estimated multiplier (${F(Mh[n - 1], 3)}) shows the present model is ${Mh[n - 1] < 1 ? 'optimistic' : 'pessimistic'}.`, 'Online disturbance estimate.'),
      rec('low', 'Hitting the length target in one batch avoids re-growth runs and rejected substrates, saving precursor, water and heating energy; prefer the lowest temperature profile that still meets the schedule.', 'First-pass yield.', true),
    ],
    publish: { 'growth.length': Lcl[n], 'ctrl.tracking_rmse': rmCl, 'growth.rate': Lcl[n] / p.t_growth, 'proc.time': p.t_growth, ...(temp ? { 'proc.T_C': zMean - 273.15 } : {}) },
    numerics: { method: 'Receding-horizon MPC: box-constrained QP by accelerated projected gradient (FISTA with restart); plant integrated with 10 sub-steps per interval', dof: Nh, steps: n, iterations: itTot, residual: pgMax, notes: 'Controller model omits actuator lag and depletion (deliberate mismatch); offset removed by a multiplicative disturbance estimate.' },
  };
}

// ---------- LQR + Kalman ----------
function care(A, B, Q, R) { // continuous algebraic Riccati equation for a 2-state, single-input system by Kleinman–Newton iteration
  let K = [0, 0], P = [[0, 0], [0, 0]], it = 0, dK = 1;
  for (; it < 60 && dK > 1e-13; it++) {
    const a = [[A[0][0] - B[0] * K[0], A[0][1] - B[0] * K[1]], [A[1][0] - B[1] * K[0], A[1][1] - B[1] * K[1]]], q = [[Q[0][0] + R * K[0] * K[0], Q[0][1] + R * K[0] * K[1]], [Q[0][1] + R * K[0] * K[1], Q[1][1] + R * K[1] * K[1]]];
    const s = N.solveLinear([[2 * a[0][0], 2 * a[1][0], 0], [a[0][1], a[0][0] + a[1][1], a[1][0]], [0, 2 * a[0][1], 2 * a[1][1]]], [-q[0][0], -q[0][1], -q[1][1]]);
    P = [[s[0], s[1]], [s[1], s[2]]]; const Kn = [(B[0] * s[0] + B[1] * s[1]) / R, (B[0] * s[1] + B[1] * s[2]) / R]; dK = Math.abs(Kn[0] - K[0]) + Math.abs(Kn[1] - K[1]); K = Kn;
  }
  const PA = (i, j) => A[0][i] * P[0][j] + A[1][i] * P[1][j] + P[i][0] * A[0][j] + P[i][1] * A[1][j], PB = [P[0][0] * B[0] + P[0][1] * B[1], P[1][0] * B[0] + P[1][1] * B[1]];
  let res = 0; for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) res = Math.max(res, Math.abs(PA(i, j) - PB[i] * PB[j] / R + Q[i][j]));
  return { K, P, it, res };
}
function lqg(p) {
  const warnings = [], w0 = 2 * N.PI * p.f0_r, Qf = p.Q_r, kS = p.m_r * w0 * w0, A = [[0, 1], [-1, -1 / Qf]], B = [0, 1], Cm = [1, 0], W = p.w_proc ** 2, V = p.v_meas ** 2;
  const lq = care(A, B, [[p.q_pos, 0], [0, p.q_vel]], p.r_u), K = lq.K, kf = care([[A[0][0], A[1][0]], [A[0][1], A[1][1]]], Cm, [[0, 0], [0, W]], V), Lk = kf.K;
  const wcl = Math.sqrt(1 + K[0]), zcl = (1 / Qf + K[1]) / (2 * wcl), b1 = Lk[0] + 1 / Qf, c1 = Lk[0] / Qf + 1 + Lk[1], west = Math.sqrt(Math.max(c1, 1e-12)), zest = b1 / (2 * west);
  const rankOf = M => Math.abs(M[0][0] * M[1][1] - M[0][1] * M[1][0]) > 1e-12 ? 2 : (Math.abs(M[0][0]) + Math.abs(M[0][1]) + Math.abs(M[1][0]) + Math.abs(M[1][1]) > 1e-12 ? 1 : 0);
  const rc = rankOf([[B[0], A[0][0] * B[0] + A[0][1] * B[1]], [B[1], A[1][0] * B[0] + A[1][1] * B[1]]]), ro = rankOf([[Cm[0], Cm[1]], [Cm[0] * A[0][0] + Cm[1] * A[1][0], Cm[0] * A[0][1] + Cm[1] * A[1][1]]]);
  const poles = (b, c) => { const d = b * b - 4 * c; return d >= 0 ? [[(-b + Math.sqrt(d)) / 2, 0], [(-b - Math.sqrt(d)) / 2, 0]] : [[-b / 2, Math.sqrt(-d) / 2], [-b / 2, -Math.sqrt(-d) / 2]]; }, pc = poles(1 / Qf + K[1], 1 + K[0]), pe = poles(b1, c1);
  const tauOpen = 2 * Qf, tauCl = 1 / Math.max(zcl * wcl, 1e-9), tsOpen = tauOpen * Math.log(50), tsClTh = tauCl * Math.log(50);
  const lamMax = Math.max(...pc.map(q => Math.hypot(q[0], q[1])), ...pe.map(q => Math.hypot(q[0], q[1])), 1);
  let Tend = Math.max(3 * tsClTh, 30), dt = Math.min(0.05, 0.5 / lamMax), ns = Math.ceil(Tend / dt); if (ns > 6000) { ns = 6000; Tend = ns * dt; }
  const rnd = N.rng(p.seed), sw = Math.sqrt(W / dt), sv = Math.sqrt(V / dt), tt = [0], X = [1], Xo = [1], Xh = [0], Uu = [0]; let s = [1, 0, 0, 0], so = [1, 0], uRms = 0, uPk = 0, estErr = 0, tsCl = 0;
  const f = (y, wn, yn) => { const u = -(K[0] * y[2] + K[1] * y[3]), inn = yn - y[2]; return [y[1], -y[0] - y[1] / Qf + u + wn, y[3] + Lk[0] * inn, -y[2] - y[3] / Qf + u + Lk[1] * inn]; };
  const rk = (y, fn) => { const k1 = fn(y), k2 = fn(y.map((v, i) => v + dt / 2 * k1[i])), k3 = fn(y.map((v, i) => v + dt / 2 * k2[i])), k4 = fn(y.map((v, i) => v + dt * k3[i])); return y.map((v, i) => v + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])); };
  for (let k = 0; k < ns; k++) {
    const wn = sw * rnd.normal(), yn = s[0] + sv * rnd.normal(), u = -(K[0] * s[2] + K[1] * s[3]); uRms += u * u; uPk = Math.max(uPk, Math.abs(u));
    s = rk(s, y => f(y, wn, yn)); so = rk(so, y => [y[1], -y[0] - y[1] / Qf + wn]); estErr += (s[0] - s[2]) ** 2;
    if (Math.abs(s[0]) > 0.02 + 3 * Math.sqrt(W / Math.max(4 * zcl * wcl ** 3, 1e-9))) tsCl = (k + 1) * dt;
    tt.push((k + 1) * dt); X.push(s[0]); Xo.push(so[0]); Xh.push(s[2]); Uu.push(u);
  }
  uRms = Math.sqrt(uRms / ns); estErr = Math.sqrt(estErr / ns); const Fpk = uPk * kS * p.x0_r, Frms = uRms * kS * p.x0_r, Pmech = Frms * (p.x0_r * w0) * 0.5, tReal = v => v / w0, Qcl = 1 / (2 * zcl), stableC = pc.every(q => q[0] < 0), stableE = pe.every(q => q[0] < 0);
  const noiseFloor = Math.sqrt(W / Math.max(4 * zcl * wcl ** 3, 1e-12)); // stationary rms displacement under LQR with full state
  if (Qcl < 0.5) warnings.push('The closed loop is over-damped (ζ > 1): position recovery becomes slow; lower the velocity weight or raise the position weight.');
  if (west < 2 * wcl) warnings.push('The estimator is not much faster than the regulator: the separation principle still guarantees stability, but robustness margins of LQG are not guaranteed.');
  if (ns === 6000) warnings.push('Simulation length capped at 6000 steps; the response plot covers only the first part of the transient.');
  return {
    summary: `LQR for a ${F(p.f0_r, 4)} Hz resonator (Q = ${F(Qf, 4)}): feedback gains K = [${F(K[0], 4)}, ${F(K[1], 4)}] (normalised) move the poles to ω = ${F(wcl, 4)}·ω₀ with damping ratio ${F(zcl, 3)} (effective Q ${F(Qcl, 3)}), shortening the 2 % settling time from ${F(tReal(tsOpen), 3)} s to ${F(tReal(tsClTh), 3)} s. The Kalman filter (gains [${F(Lk[0], 4)}, ${F(Lk[1], 4)}]) has poles at ω = ${F(west, 4)}·ω₀, ζ = ${F(zest, 3)}. Peak actuator force ${F(Fpk, 3)} N for a ${F(p.x0_r, 3)} m disturbance.`,
    kpis: [
      kpi('zeta_cl', 'Closed-loop damping ratio', zcl, '–', `From the characteristic polynomial s² + (1/Q + k₂)s + (1 + k₁) of A − BK; the open-loop value is ${F(1 / (2 * Qf), 3)}.`, zcl > 0.2 && zcl <= 1),
      kpi('Q_cl', 'Closed-loop quality factor', Qcl, '–', 'Effective Q = 1/(2ζ) of the actively damped resonator; the reduction relative to the passive Q is the active damping achieved.', null),
      kpi('f_cl', 'Closed-loop resonance', p.f0_r * wcl, 'Hz', `Position feedback adds active stiffness k₁ = ${F(K[0], 3)} times the passive stiffness, shifting the resonance by ${F(100 * (wcl - 1), 3)} %.`, null),
      kpi('ts_cl', 'Settling time (closed loop)', tReal(tsClTh), 's', `ln(50)/(ζω) for the closed-loop poles; in the noisy simulation the response stayed above the noise-adjusted 2 % band until ${F(tReal(tsCl), 3)} s.`, tsClTh < tsOpen),
      kpi('ts_open', 'Settling time (open loop)', tReal(tsOpen), 's', 'Free ring-down time to 2 % of the initial amplitude, 2Q·ln(50)/ω₀.', null),
      kpi('F_peak', 'Peak control force', Fpk, 'N', `Largest actuator force in the simulation, |u|·k·x₀ with k = mω₀² = ${F(kS, 3)} N/m. Compare with the actuator saturation limit.`, null),
      kpi('F_rms', 'RMS control force', Frms, 'N', 'Root-mean-square actuator force over the simulated transient, including the response to process and sensor noise.', null),
      kpi('P_ctrl', 'Indicative control power', Pmech, 'W', 'Order-of-magnitude mechanical power ½·F_rms·x₀ω₀ exchanged by the actuator; electrical drive power is higher by the inverse actuator efficiency.', null),
      kpi('est_err', 'RMS estimation error', estErr, 'x₀', 'Root-mean-square difference between true and estimated displacement in the simulation, including the initial transient.', estErr < 0.3),
      kpi('rank_c', 'Controllability rank', rc, '–', 'Rank of [B, AB]; 2 means both states can be driven by the single force actuator.', rc === 2),
      kpi('rank_o', 'Observability rank', ro, '–', 'Rank of [C; CA]; 2 means displacement measurement alone suffices to reconstruct velocity.', ro === 2),
      kpi('are_res', 'Riccati residual', Math.max(lq.res, kf.res), '–', `Largest element of AᵀP + PA − PBR⁻¹BᵀP + Q after ${lq.it} (control) and ${kf.it} (filter) Kleinman–Newton iterations.`, Math.max(lq.res, kf.res) < 1e-8),
    ],
    series: [
      { key: 'x_t', title: 'Displacement response', xlabel: 'time (s)', ylabel: 'displacement / x₀', lines: [line('closed loop (LQG)', tt.map(tReal), X), line('open loop', tt.map(tReal), Xo), line('estimate', tt.map(tReal), Xh)], explain: 'Response to an initial displacement with seeded process and sensor noise. The controlled trace decays within a few periods; the estimate converges to the true state at the estimator rate.' },
      { key: 'u_t', title: 'Control force', xlabel: 'time (s)', ylabel: 'force (N)', lines: [line('actuator force', tt.map(tReal), Uu.map(v => v * kS * p.x0_r))], explain: 'Actuator force commanded by the state feedback acting on the estimated state. The noise-driven part at late times is the price of using a noisy sensor.' },
    ],
    tables: [
      { key: 'poles', title: 'Poles (in units of ω₀)', columns: ['System', 'Real part', 'Imaginary part'], rows: [['open loop', F(-1 / (2 * Qf), 4), F(Math.sqrt(Math.max(1 - 1 / (4 * Qf * Qf), 0)), 4)], ['regulator 1', F(pc[0][0], 4), F(pc[0][1], 4)], ['regulator 2', F(pc[1][0], 4), F(pc[1][1], 4)], ['estimator 1', F(pe[0][0], 4), F(pe[0][1], 4)], ['estimator 2', F(pe[1][0], 4), F(pe[1][1], 4)]], explain: 'Eigenvalues of A (open loop), A − BK (regulator) and A − LC (estimator), normalised by the natural frequency. All real parts must be negative.' },
      mechTable([['Physical mechanism', 'Force actuator (piezoelectric, electrostatic or electromagnetic) driven by state feedback: synthetic stiffness and damping'], ['Energy demand', `≈ ${F(Pmech, 3)} W mechanical while rejecting this disturbance; zero when at rest apart from noise-driven action`], ['Response time', `${F(tReal(tsClTh), 3)} s to settle (open loop ${F(tReal(tsOpen), 3)} s)`], ['Loss', `effective Q reduced from ${F(Qf, 3)} to ${F(Qcl, 3)}; sensor noise is injected as force noise`], ['Reversibility', 'Fully reversible: gains can be reprogrammed electronically at any time'], ['Hysteresis', 'None in the linear design; actuator hysteresis would appear as a nonlinearity inside the loop'], ['Cycle durability', 'Limited by the actuator (for example piezoelectric fatigue), not by the control law']]),
    ],
    equations: [
      eqn('governing', 'State-space model (normalised)', 'ẋ = A x + B u + G w,  y = C x + v;  A = [[0, 1], [−1, −1/Q]], B = G = [0, 1]ᵀ, C = [1, 0]'),
      eqn('formula', 'LQR', 'min ∫(xᵀQx + r u²)dt  ⇒  u = −Kx̂,  K = r⁻¹BᵀP,  AᵀP + PA − PBr⁻¹BᵀP + Q = 0'),
      eqn('formula', 'Kalman–Bucy filter', 'dx̂/dt = A x̂ + B u + L(y − C x̂),  L = P_f Cᵀ V⁻¹,  AP_f + P_f Aᵀ − P_f CᵀV⁻¹C P_f + GWGᵀ = 0'),
      eqn('formula', 'Closed-loop poles', 's² + (1/Q + k₂)s + (1 + k₁) = 0'),
    ],
    checks: [
      chk('Regulator poles in the left half-plane', stableC, `Re = ${F(pc[0][0], 4)}, ${F(pc[1][0], 4)}`),
      chk('Estimator poles in the left half-plane', stableE, `Re = ${F(pe[0][0], 4)}, ${F(pe[1][0], 4)}`),
      chk('System controllable and observable', rc === 2 && ro === 2, `ranks ${rc} and ${ro}`),
      chk('Riccati equations solved', Math.max(lq.res, kf.res) < 1e-8, `residuals ${lq.res.toExponential(1)} and ${kf.res.toExponential(1)}`),
      chk('Riccati solutions positive definite', lq.P[0][0] > 0 && lq.P[0][0] * lq.P[1][1] - lq.P[0][1] ** 2 > 0 && kf.P[0][0] > 0, `P₁₁ = ${F(lq.P[0][0], 4)}, det P = ${F(lq.P[0][0] * lq.P[1][1] - lq.P[0][1] ** 2, 4)}`),
      chk('Simulated response bounded', X.every(v => Math.abs(v) < 5), `max |x|/x₀ = ${F(N.maxOf(X.map(Math.abs)), 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', `Check that the actuator can deliver ${F(Fpk, 3)} N: saturation inside an LQG loop degrades damping. If it cannot, raise the control weight r (now ${F(p.r_u)}).`, 'Peak control force in the simulation.'),
      rec('medium', west < 2 * wcl ? 'Speed up the estimator (lower the assumed measurement noise or raise the process noise) so that its poles are two to five times faster than the regulator poles.' : 'Estimator and regulator are well separated; verify robustness with a 20 % error in resonance frequency before deployment.', `Estimator/regulator bandwidth ratio ${F(west / wcl, 3)}.`),
      rec('medium', 'Full state feedback has guaranteed margins (at least 60° phase), but the observer-based controller does not: measure the loop transfer function on the hardware and add a notch or roll-off for unmodelled higher modes.', 'Known limitation of LQG robustness.'),
      rec('low', `Active damping replaces lossy passive treatments: the same vibration attenuation is obtained with about ${F(Pmech, 2)} W only while a disturbance is present, and without added mass or viscoelastic materials that are hard to recycle.`, 'Control power estimate.', true),
    ],
    publish: { 'res.Q': Qcl, 'res.f0': p.f0_r * wcl, 'ctrl.settling_s': tReal(tsClTh), 'ctrl.force_peak': Fpk },
    numerics: { method: 'Kleinman–Newton iteration on the algebraic Riccati equations (3×3 Lyapunov solves); RK4 simulation with seeded band-limited white noise', dof: 4, steps: ns, iterations: lq.it + kf.it, residual: Math.max(lq.res, kf.res), notes: `Time normalised by ω₀; step ${F(dt, 3)} rad.` },
  };
}

// ---------- varactor SRR ----------
function varactor(p) {
  const warnings = [], L = p.L_srr * 1e-9, Cg = p.C_gap * 1e-12, Cj0 = p.Cj0 * 1e-12, Vmax = Math.max(p.V_max, p.V_bias), Cj = V => Cj0 / Math.pow(1 + V / p.Vbi, p.Mj), Ct = V => Cg + Cj(V), f = V => 1 / (2 * N.PI * Math.sqrt(L * Ct(V))), Qv = V => Math.sqrt(L / Ct(V)) / p.R_s;
  const fB = f(p.V_bias), f0 = f(0), fM = f(Vmax), range = 2 * (fM - f0) / (fM + f0), sens = (f(p.V_bias + 0.01) - f(Math.max(p.V_bias - 0.01, 0))) / (p.V_bias >= 0.01 ? 0.02 : 0.01), Q = Qv(p.V_bias), bw = fB / Q, fom = (fM - f0) / bw;
  const tau = p.R_bias * 1e3 * Cj(p.V_bias), vs = N.linspace(0, Vmax, 121), Esw = N.trapz(vs, vs.map(v => v * Cj(v))), Estore = 0.5 * Ct(p.V_bias) * p.V_bias ** 2, Phold = p.V_bias * 1e-8, capRatio = Cj(0) / Cj(Vmax), maxRange = Math.sqrt(Ct(0) / Ct(Vmax)) - 1;
  if (Q < 5) warnings.push(`Quality factor ${F(Q, 3)} is low: the resonance is broad and the tunable shift is less than a few linewidths. Reduce series resistance.`);
  if (Cg > Cj(Vmax)) warnings.push('The fixed gap capacitance exceeds the minimum varactor capacitance and dilutes the tuning range; reduce parasitics or use a larger varactor.');
  warnings.push('Small-signal model: at high incident power the RF voltage across the junction modulates the capacitance (self-tuning, harmonic generation).');
  return {
    summary: `Varactor-loaded split-ring resonator: resonance ${F(fB / 1e9, 4)} GHz at ${F(p.V_bias)} V reverse bias, tunable from ${F(f0 / 1e9, 4)} GHz (0 V) to ${F(fM / 1e9, 4)} GHz (${F(Vmax)} V), a ${F(100 * range, 3)} % range equal to ${F(fom, 3)} linewidths (Q = ${F(Q, 3)}). Sensitivity ${F(sens / 1e6, 3)} MHz/V; tuning speed limited by the bias network to ${F(tau * 1e6, 3)} µs; switching energy ${F(Esw * 1e12, 3)} pJ per cell.`,
    kpis: [
      kpi('f_bias', 'Resonance at operating bias', fB, 'Hz', 'f = 1/(2π√(L(C_gap + C_j(V)))) with the reverse-biased junction capacitance C_j = C_j0/(1 + V/V_bi)^M.', null),
      kpi('range', 'Tuning range', 100 * range, '%', `Fractional frequency span between 0 V and ${F(Vmax)} V relative to the centre; the junction capacitance ratio is ${F(capRatio, 3)}.`, range > 0.1),
      kpi('sens', 'Tuning sensitivity', sens, 'Hz/V', 'Slope df/dV at the operating bias; it converts bias noise and drift into frequency error.', null),
      kpi('Q', 'Quality factor', Q, '–', '√(L/C)/R_s of the loaded resonator; the series resistance of the varactor usually dominates the loss.', Q >= 10),
      kpi('fom', 'Tuning figure of merit', fom, 'linewidths', 'Tuning span divided by the resonance linewidth f/Q: the number of resolvable states (above about 1 is needed for useful switching).', fom > 1),
      kpi('tau', 'Bias response time', tau, 's', 'R_bias·C_j time constant of the bias network charging the junction.', null),
      kpi('E_sw', 'Switching energy per cell', Esw, 'J', 'Energy drawn from the bias supply to charge the junction from 0 V to the maximum bias, ∫V dQ.', null),
      kpi('P_hold', 'Holding power per cell', Phold, 'W', 'Bias voltage times an assumed 10 nA reverse leakage; essentially zero static power.', null),
      kpi('C_bias', 'Total capacitance at bias', Ct(p.V_bias), 'F', 'Gap plus junction capacitance at the operating point.', null),
    ],
    series: [
      { key: 'f_V', title: 'Resonance versus bias', xlabel: 'reverse bias (V)', ylabel: 'resonance (GHz)', lines: [line('f(V)', vs, vs.map(v => f(v) / 1e9))], explain: 'Tuning curve: steep at low bias, flattening as the junction capacitance approaches the fixed gap capacitance.' },
      { key: 'Q_V', title: 'Quality factor and capacitance versus bias', xlabel: 'reverse bias (V)', ylabel: 'Q and C (pF)', lines: [line('Q', vs, vs.map(Qv)), line('C total (pF)', vs, vs.map(v => Ct(v) * 1e12))], explain: 'Q rises as the capacitance falls because the characteristic impedance √(L/C) increases at fixed series resistance.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Depletion-width modulation of a reverse-biased p–n (or Schottky) junction changes the capacitance loading the resonator gap'], ['Energy demand', `${F(Esw * 1e12, 3)} pJ per switching event; holding power ≈ ${F(Phold * 1e9, 3)} nW (leakage)`], ['Response time', `${F(tau * 1e6, 3)} µs (bias network); intrinsic junction response is sub-nanosecond`], ['Loss', `series resistance ${F(p.R_s)} Ω gives Q = ${F(Q, 3)} (insertion loss rises towards low bias)`], ['Reversibility', 'Fully reversible and continuous (analogue tuning)'], ['Hysteresis', 'None for a junction varactor (ferroelectric varactors differ)'], ['Cycle durability', 'Practically unlimited electrical cycling; limits are reverse breakdown and electrostatic discharge']])],
    equations: [
      eqn('constitutive', 'Junction capacitance', 'C_j(V) = C_j0 / (1 + V/V_bi)^M'),
      eqn('formula', 'Resonance and quality factor', 'f = 1/[2π√(L(C_gap + C_j))];  Q = √(L/C)/R_s'),
      eqn('formula', 'Maximum tuning ratio', 'f_max/f_min = √[(C_gap + C_j0)/(C_gap + C_j(V_max))]', `= ${F(1 + maxRange, 4)}`),
    ],
    checks: [
      chk('Tuning span exceeds one linewidth', fom > 1, `${F(fom, 3)} linewidths`),
      chk('Frequency increases monotonically with bias', vs.every((v, i) => i === 0 || f(v) >= f(vs[i - 1])), 'capacitance falls with reverse bias'),
      chk('Passive (positive resistance and capacitance)', p.R_s > 0 && Ct(Vmax) > 0, `R_s = ${F(p.R_s)} Ω`),
      chk('Bias within stated maximum', p.V_bias <= p.V_max, `${F(p.V_bias)} V ≤ ${F(p.V_max)} V`),
    ],
    warnings,
    recommendations: [
      rec('high', fom < 2 ? 'Lower the series resistance (higher-Q varactor, thicker metal) or the parasitic gap capacitance: the tuning span is under two linewidths.' : `The design resolves ${F(fom, 3)} linewidths; a ${F(Math.abs(bw / sens) * 1e3, 3)} mV bias resolution is enough to place the resonance within one linewidth.`, 'Tuning figure of merit.'),
      rec('medium', 'Use hyper-abrupt varactors (M > 1) if a wider or more linear tuning curve is required; they trade Q for capacitance ratio.', 'Grading coefficient dependence.'),
      rec('medium', `Bias-line RC of ${F(tau * 1e6, 3)} µs limits reconfiguration to about ${F(0.35 / tau / 1e3, 3)} kHz; lower the bias resistance if faster modulation is needed, at the cost of RF leakage.`, 'Bias network time constant.'),
      rec('low', 'Varactor tuning draws only leakage current in the static state, so large programmable surfaces can run from milliwatt-level supplies — far less than thermal or mechanical tuning.', 'Holding power estimate.', true),
    ],
    publish: { 'res.f0': fB, 'res.Q': Q, 'tune.range': range, 'tune.energy_J': Esw, 'tune.response_s': tau },
    numerics: { method: 'Closed-form lumped RLC with depletion-capacitance law; trapezoidal charge–voltage integral', dof: 1, steps: 121, iterations: 1, residual: 0, notes: 'Quasi-static, small-signal, lumped-element model of one unit cell.' },
  };
}

// ---------- phase-change metasurface ----------
const PCMS = { vo2: { f: 'VO2', name: 'VO₂', Tup: 68, Tdn: 62, w: 1.5, si: 50, sm: 3e5, rho: 4340, cp: 690, lat: 5.15e4, volatile: true, Treset: 0, latReset: 0 }, gst: { f: 'Ge2Sb2Te5', name: 'Ge₂Sb₂Te₅', Tup: 150, Tdn: -1e9, w: 5, si: 3, sm: 5e4, rho: 6150, cp: 218, lat: 0, volatile: false, Treset: 615, latReset: 1e5 } };
function pcm(p) {
  const warnings = [], M = PCMS[p.pcm_mat], A = p.area_um2 * 1e-12, th = p.thick_nm * 1e-9, Vol = A * th, P = p.P_mw * 1e-3, Tb = p.T_base, Cth = p.Cth_mult * M.rho * M.cp * Vol, tauTh = Cth / p.G_th, Tpk = Tb + P / p.G_th;
  const frUp = T => 1 / (1 + Math.exp(-(T - M.Tup) / M.w)), frDn = T => M.volatile ? 1 / (1 + Math.exp(-(T - M.Tdn) / M.w)) : 1, sig = fr => Math.exp((1 - fr) * Math.log(M.si) + fr * Math.log(M.sm));
  const w0 = 2 * N.PI * p.f0_ms, Lr = p.X0 / w0, Cr = 1 / (p.X0 * w0), Z0 = 376.73;
  const Tr = (f, fr) => { const w = 2 * N.PI * f, Rv = 1 / (sig(fr) * th), g = 1 / Rv, b = w * Cr, d = g * g + b * b, zr = p.R_ms + g / d, zi = w * Lr - b / d; return (4 * (zr * zr + zi * zi)) / ((2 * zr + Z0) ** 2 + 4 * zi * zi); }; // shunt sheet between free-space half-spaces
  const fs = N.linspace(0.4 * p.f0_ms, 1.6 * p.f0_ms, 241), Tins = fs.map(f => Tr(f, 0)), Tmet = fs.map(f => Tr(f, 1)), iMin = N.argmin(Tins), fRes = fs[iMin], T0 = Tr(fRes, 0), T1 = Tr(fRes, 1), depth = Math.abs(T1 - T0), contrast = 10 * Math.log10(Math.max(T1, 1e-12) / Math.max(T0, 1e-12));
  const Tsw = M.Tup + 3 * M.w, reach = Tpk > Tsw, dTs = Tsw - Tb, Eset = Cth * dTs + M.rho * Vol * M.lat, tOn = reach ? -tauTh * Math.log(1 - dTs / (P / p.G_th)) + M.rho * Vol * M.lat / Math.max(P - p.G_th * dTs, 1e-30) : tauTh * 5;
  const Phold = M.volatile ? p.G_th * dTs : 0, tOff = M.volatile ? (reach ? tauTh * Math.log(Math.max((Tpk - Tb) / Math.max(M.Tdn - 3 * M.w - Tb, 1e-6), 1.0001)) : tauTh) : tauTh * Math.log(Math.max((M.Treset - Tb) / (M.Tup - Tb), 1.0001));
  const Ereset = M.volatile ? 0 : Cth * (M.Treset - Tb) + M.rho * Vol * M.latReset, Preset = M.volatile ? 0 : p.G_th * (M.Treset - Tb), canReset = M.volatile || Tpk > M.Treset, quench = M.volatile ? 0 : (M.Treset - M.Tup) / Math.max(tOff, 1e-30);
  const Tmax = Math.max(M.Tup + 30, Math.min(Tpk, M.Tup + 60)), Tu = N.linspace(Tb, Tmax, 80), Td = N.linspace(Tmax, Tb, 80), loopArea = N.trapz(Tu, Tu.map(T => frDn(T) - frUp(T)));
  if (!reach) warnings.push(`Heater power is insufficient: the steady-state temperature ${F(Tpk, 4)} °C stays below the switching temperature ${F(Tsw, 4)} °C. At least ${F(p.G_th * dTs * 1e3, 3)} mW per cell is needed.`);
  if (!M.volatile && !canReset) warnings.push(`This heater cannot reach the melt-quench temperature (${F(M.Treset)} °C): the cell can be set (crystallised) but not reset to the amorphous state; ${F(Preset * 1e3, 3)} mW would be needed.`);
  if (!M.volatile && quench < 1e9) warnings.push(`Cooling rate after reset ≈ ${F(quench, 2)} K/s; amorphisation of GST typically needs of the order of 10⁹ K/s or more, so reduce the heated thermal mass or raise the thermal conductance.`);
  if (M.volatile && Tb > M.Tdn - 5) warnings.push('The base temperature is within 5 K of the cooling transition: the cell may not return to the insulating state.');
  if (p.f0_ms > 3e13) warnings.push('Above about 30 THz the phase-change film acts mainly through its refractive-index contrast rather than its DC-like sheet conductance; the conductive-short model used here is then only qualitative.');
  warnings.push('Material constants (transition temperatures, conductivities, latent heat) are typical thin-film values and vary strongly with stoichiometry, strain and deposition; calibrate on measured films.');
  return {
    summary: `${M.name} phase-change metasurface cell (${F(p.area_um2)} µm² × ${F(p.thick_nm)} nm): transmission at the ${F(fRes / 1e12, 4)} THz resonance switches from ${F(T0, 3)} (insulating/amorphous) to ${F(T1, 3)} (metallic/crystalline), a modulation depth of ${F(depth, 3)} (${F(contrast, 3)} dB). Switching on needs ${F(Eset * 1e9, 3)} nJ and ${F(tOn * 1e6, 3)} µs at ${F(p.P_mw)} mW; ${M.volatile ? `holding the metallic state costs ${F(Phold * 1e3, 3)} mW per cell, and the hysteresis is ${F(M.Tup - M.Tdn, 3)} K wide` : `the state is retained without power, and reset by melt-quench needs ${F(Ereset * 1e9, 3)} nJ`}.`,
    kpis: [
      kpi('T_off', 'Transmission, insulating state', T0, '–', 'Power transmission at resonance through the equivalent shunt RLC sheet with the patch in its low-conductivity state (deep resonance dip).', null),
      kpi('T_on', 'Transmission, conducting state', T1, '–', `Transmission at the same frequency when the film (sheet resistance ${F(1 / (M.sm * th), 3)} Ω per square) shorts the capacitive gap and quenches the resonance.`, null),
      kpi('depth', 'Modulation depth', depth, '–', 'Absolute change of power transmission between the two states at the resonance frequency.', depth > 0.3),
      kpi('E_set', 'Switching-on energy', Eset, 'J', `Sensible heat of the heated volume (${F(Cth, 3)} J/K) from ${F(Tb)} °C to ${F(Tsw, 4)} °C${M.lat ? ' plus the latent heat of the transition' : ''}.`, null),
      kpi('t_on', 'Switch-on time', tOn, 's', `Thermal RC heating time (τ = ${F(tauTh, 3)} s) to the switching temperature at the stated heater power${reach ? '' : ' — not reachable; five time constants reported'}.`, reach),
      kpi('t_off', M.volatile ? 'Switch-off (cooling) time' : 'Quench time after reset', tOff, 's', M.volatile ? 'Time to cool passively from the peak temperature to below the cooling branch of the hysteresis loop.' : 'Time to cool from the melt temperature to the crystallisation temperature; it must be short enough to freeze the amorphous state.', null),
      kpi('P_hold', 'Holding power per cell', Phold, 'W', M.volatile ? 'Steady heater power needed to keep the film above the transition, G_th·ΔT: the continuous cost of a volatile mechanism.' : 'Zero: the amorphous and crystalline states are both stable at the base temperature (non-volatile).', Phold < 1e-3),
      kpi('E_reset', 'Reset energy', Ereset, 'J', M.volatile ? 'Not applicable to the volatile VO₂ transition (cooling restores the insulating state); reported as 0.' : `Heat needed to melt the film (${F(M.Treset)} °C plus latent heat) before quenching to the amorphous state.`, null),
      kpi('hyst', 'Thermal hysteresis width', M.volatile ? M.Tup - M.Tdn : M.Treset - M.Tup, 'K', M.volatile ? 'Separation of the heating and cooling transition temperatures (two-branch sigmoid model); it provides memory within the loop but demands temperature overdrive.' : 'Separation between the crystallisation and melt temperatures: the window that makes the state non-volatile.', null),
      kpi('T_peak', 'Steady-state temperature at heater power', Tpk + 273.15, 'K', 'Base temperature plus P/G_th; it must exceed the switching temperature with margin but stay below damage limits.', reach),
    ],
    series: [
      { key: 'T_f', title: 'Transmission spectra of the two states', xlabel: 'frequency (THz)', ylabel: 'power transmission', lines: [line('insulating / amorphous', fs.map(v => v / 1e12), Tins), line('metallic / crystalline', fs.map(v => v / 1e12), Tmet)], explain: 'Resonant dip of the metasurface with the patch insulating, and its suppression when the patch conducts. The vertical distance at the dip is the modulation depth.' },
      { key: 'loop', title: 'Hysteresis loop: transmission at resonance versus temperature', xlabel: 'temperature (°C)', ylabel: 'power transmission', lines: [line('heating', Tu, Tu.map(T => Tr(fRes, frUp(T)))), line('cooling', Td, Td.map(T => Tr(fRes, frDn(T))))], explain: M.volatile ? 'Heating and cooling branches of the insulator–metal transition. The loop width is the thermal hysteresis; inside it the state depends on history.' : 'Heating crystallises the film; on cooling it stays crystalline (non-volatile memory). Only melt-quenching, outside this plot, restores the initial state.' },
    ],
    tables: [mechTable([['Physical mechanism', M.volatile ? 'Thermally driven insulator–metal (monoclinic–rutile) transition of VO₂ near 68 °C' : 'Thermally driven amorphous–crystalline transition of Ge₂Sb₂Te₅ (set by annealing, reset by melt-quench)'], ['Energy demand', `set ${F(Eset * 1e9, 3)} nJ per cell; ${M.volatile ? 'hold ' + F(Phold * 1e3, 3) + ' mW per cell' : 'reset ' + F(Ereset * 1e9, 3) + ' nJ, no holding power'}`], ['Response time', `on ${F(tOn * 1e6, 3)} µs, ${M.volatile ? 'off' : 'quench'} ${F(tOff * 1e6, 3)} µs (thermal RC ${F(tauTh * 1e6, 3)} µs); the intrinsic transition is far faster`], ['Loss', `conducting-state sheet resistance ${F(1 / (M.sm * th), 3)} Ω/sq adds absorption; insulating-state conductivity ${F(M.si)} S/m`], ['Reversibility', M.volatile ? 'Reversible and volatile: returns to the insulating state on cooling' : 'Reversible but requires two different thermal pulses; non-volatile between pulses'], ['Hysteresis', M.volatile ? `${F(M.Tup - M.Tdn, 3)} K thermal hysteresis (two-branch model, area ${F(loopArea, 3)} K)` : 'Memory-type: state persists until the opposite pulse is applied'], ['Cycle durability', M.volatile ? 'Thin films commonly survive very large cycle numbers; cracking from the ~1 % transformation strain limits thick films' : 'Typically 10⁵–10⁹ cycles in memory cells; large-area films degrade sooner through segregation and void formation']])],
    equations: [
      ...formulaEq(M.f, M.volatile ? 'monoclinic (insulating) ⇌ rutile (metallic)' : 'amorphous ⇌ crystalline'),
      eqn('constitutive', 'Two-branch hysteresis', 'f_m↑(T) = 1/{1 + exp[−(T − T_c↑)/w]},  f_m↓(T) = 1/{1 + exp[−(T − T_c↓)/w]};  ln σ = (1 − f_m) ln σ_i + f_m ln σ_m'),
      eqn('governing', 'Thermal RC', 'C_th dT/dt = P − G_th(T − T_base);  τ = C_th/G_th'),
      eqn('formula', 'Sheet transmission', 'T = |2Z_s/(2Z_s + Z₀)|²,  Z_s = R + jωL + (jωC + 1/R_film)⁻¹'),
    ],
    checks: [
      chk('Switching temperature reachable', reach, `T_peak = ${F(Tpk, 4)} °C vs ${F(Tsw, 4)} °C`),
      chk('Reset achievable', canReset, M.volatile ? 'volatile: reset by cooling' : `needs ${F(M.Treset)} °C, heater reaches ${F(Tpk, 4)} °C`),
      chk('Passivity of the sheet model (0 ≤ T ≤ 1)', Tins.every(v => v >= 0 && v <= 1 + 1e-9) && Tmet.every(v => v >= 0 && v <= 1 + 1e-9), `T range ${F(Math.min(N.minOf(Tins), N.minOf(Tmet)), 3)}–${F(Math.max(N.maxOf(Tins), N.maxOf(Tmet)), 3)}`),
      chk('Useful modulation depth', depth > 0.3, `ΔT = ${F(depth, 3)}`),
      chk('Peak temperature below damage limit', Tpk < (M.volatile ? 400 : 900), `${F(Tpk, 4)} °C`),
    ],
    warnings,
    recommendations: [
      rec('high', !reach ? `Increase heater power above ${F(p.G_th * dTs * 1.5e3, 3)} mW per cell or reduce the thermal conductance to the substrate.` : M.volatile ? `Holding power of ${F(Phold * 1e3, 3)} mW per cell scales with array size: for a static configuration prefer a non-volatile material; for fast modulation VO₂ is appropriate.` : 'The state is retained without power: ideal for reconfigure-rarely surfaces; verify that reset pulses reach the melt temperature with the required quench rate.', 'Thermal power balance.'),
      rec('medium', depth < 0.3 ? 'Increase film thickness or resonator reactance so that the conducting film shorts the gap more effectively; the present modulation depth is small.' : `Modulation depth ${F(depth, 3)} is useful; thin the film until the depth starts to fall, which lowers the switching energy proportionally.`, 'Sheet-impedance model.'),
      rec('medium', `Thermal isolation trades speed for power: τ = ${F(tauTh * 1e6, 3)} µs at G_th = ${F(p.G_th, 3)} W/K. Suspended or low-conductivity supports cut holding power but slow cooling.`, 'Thermal RC time constant.'),
      rec('low', M.volatile ? 'Pulse the heater only while the function is needed and operate the array close to (but safely below) the transition to minimise overdrive energy.' : 'Non-volatile switching eliminates standby power entirely, the dominant energy term for rarely reconfigured surfaces; tellurium and germanium are critical raw materials, so minimise film volume.', 'Energy and materials demand of the mechanism.', true),
    ],
    publish: { 'res.transmission': T1, 'res.f0': fRes, 'tune.energy_J': Eset, 'tune.response_s': tOn, 'tune.hold_power_W': Phold, 'therm.T_max': Tpk + 273.15, 'material.formula': M.f },
    numerics: { method: 'Two-branch sigmoid hysteresis, lumped thermal RC and equivalent-circuit sheet transmission (closed form)', dof: 1, steps: 241, iterations: 1, residual: 0, notes: 'Logarithmic mixing of conductivity; uniform patch temperature.' },
  };
}

// ---------- piezoelectric ----------
function piezo(p) {
  const warnings = [], n = Math.round(p.n_lay), t = p.t_lay * 1e-6, A = p.A_pz * 1e-6, Ls = n * t, Y = p.Y_pz * 1e9, k2 = p.k33 ** 2, ksc = Y * A / Ls, koc = ksc / (1 - k2), stroke = n * p.d33 * 1e-12 * p.V_pz, Fb = ksc * stroke, CT = n * N.EPS0 * p.epsr * A / t, CS = CT * (1 - k2), al = p.alpha_sh;
  const kOf = a => ksc * (1 + k2 / (1 - k2) / (1 + a)), kSh = kOf(al), fOf = k => Math.sqrt(k / p.m_pz) / (2 * N.PI), fsc = fOf(ksc), foc = fOf(koc), fsh = fOf(kSh), field = p.V_pz / t, Est = 0.5 * CT * p.V_pz ** 2, tanD = 0.02, Eloss = N.PI * tanD * CT * p.V_pz ** 2 / 4, tauE = p.R_drv * CT, fAct = Math.sqrt(Y / 7600) / (2 * Ls), tResp = Math.max(2.2 * tauE, 1 / (3 * fAct));
  const Ropt = Math.sqrt(1 - k2) / (2 * N.PI * fsc * CS), etaMax = k2 / (2 * Math.sqrt(1 - k2)), stableNC = 1 + al > 0, range = foc / fsc - 1, Wmech = 0.5 * Fb * stroke / 2;
  if (field > 2e6) warnings.push(`Electric field ${F(field / 1e6, 3)} MV/m exceeds the usual 2 MV/m limit of soft PZT: depolarisation or breakdown risk.`);
  if (al < 0) warnings.push(`Negative-capacitance shunt (ratio ${F(al)}): an active circuit is required and the system is stable only while C_shunt > −C_blocked; margin here ${F(1 + al, 3)}.`);
  if (fsh > 0.3 * fAct) warnings.push('The mass–stack resonance approaches the stack\'s own longitudinal resonance; the lumped-spring model loses accuracy.');
  const als = N.logspace(0.01, 100, 60);
  return {
    summary: `Piezoelectric stack (${n} layers × ${F(p.t_lay)} µm, ${F(p.A_pz)} mm²): free stroke ${F(stroke * 1e6, 4)} µm and blocked force ${F(Fb, 4)} N at ${F(p.V_pz)} V. Electrical boundary conditions tune its stiffness from ${F(ksc / 1e6, 4)} (short circuit) to ${F(koc / 1e6, 4)} N/µm (open circuit), shifting the resonance with a ${F(p.m_pz)} kg mass from ${F(fsc, 4)} to ${F(foc, 4)} Hz (${F(100 * range, 3)} %); with the selected shunt the resonance is ${F(fsh, 4)} Hz. Energy per actuation ${F(Est * 1e3, 3)} mJ, response time ${F(tResp * 1e6, 3)} µs.`,
    kpis: [
      kpi('stroke', 'Free stroke', stroke, 'm', `n·d₃₃·V for the unloaded stack at ${F(field / 1e6, 3)} MV/m; real stacks reach about 0.1 % strain.`, field <= 2e6),
      kpi('F_block', 'Blocked force', Fb, 'N', 'Short-circuit stiffness times free stroke: the force when displacement is fully prevented.', null),
      kpi('k_sc', 'Short-circuit stiffness', ksc, 'N/m', 'Y·A/L with shorted electrodes (lowest stiffness state).', null),
      kpi('k_oc', 'Open-circuit stiffness', koc, 'N/m', `k_sc/(1 − k₃₃²): the charge generated by strain cannot flow and stiffens the stack by a factor ${F(1 / (1 - k2), 3)}.`, null),
      kpi('k_sh', 'Stiffness with selected shunt', kSh, 'N/m', `k_sc[1 + k²/(1 − k²)·C_S/(C_S + C_sh)] for a capacitive shunt ratio ${F(al)}; negative ratios exceed the open-circuit value.`, stableNC && kSh > 0),
      kpi('f_range', 'Resonance tuning range', 100 * range, '%', `Passive short–open tuning 1/√(1 − k²) − 1 of the mass–stack resonance (${F(fsc, 4)}–${F(foc, 4)} Hz).`, range > 0.05),
      kpi('f_sh', 'Resonance with selected shunt', fsh, 'Hz', 'Resonance of the supported mass on the shunted stack.', null),
      kpi('C', 'Stack capacitance', CT, 'F', `Free capacitance n·ε·A/t; the blocked value is ${F(CS, 3)} F.`, null),
      kpi('E_act', 'Energy per actuation', Est, 'J', `½CV² drawn per full-stroke cycle unless recovered; about ${F(Eloss, 3)} J of it is dissipated as heat per cycle (tan δ = 0.02 assumed).`, null),
      kpi('t_resp', 'Response time', tResp, 's', `Larger of the electrical charging time 2.2·R·C (${F(2.2 * tauE, 3)} s) and one third of the stack's mechanical period (resonance ${F(fAct, 3)} Hz).`, null),
      kpi('eta_max', 'Maximum shunt-damping loss factor', etaMax, '–', `Peak loss factor k²/(2√(1 − k²)) of a resistive shunt, reached at R = ${F(Ropt, 3)} Ω near the short-circuit resonance.`, null),
    ],
    series: [
      { key: 'k_alpha', title: 'Stiffness versus shunt capacitance', xlabel: 'C_shunt / C_blocked', ylabel: 'stiffness (N/µm)', logx: true, lines: [line('shunted', als, als.map(a => kOf(a) / 1e6)), line('open circuit', [0.01, 100], [koc / 1e6, koc / 1e6]), line('short circuit', [0.01, 100], [ksc / 1e6, ksc / 1e6])], explain: 'A switchable capacitor bank moves the stiffness continuously between the open- and short-circuit limits; the range is set by the coupling coefficient.' },
      { key: 'F_x', title: 'Actuator load line', xlabel: 'displacement (µm)', ylabel: 'force (N)', lines: [line(`${F(p.V_pz)} V`, [0, stroke * 1e6], [Fb, 0]), line(`${F(p.V_pz / 2)} V`, [0, stroke * 1e6 / 2], [Fb / 2, 0])], explain: 'Force–displacement capability at constant voltage: blocked force at zero displacement, free stroke at zero force. A spring load intersects the line at the operating point.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Converse piezoelectric effect (actuation) and electromechanical coupling to a shunt circuit (stiffness and damping tuning)'], ['Energy demand', `${F(Est * 1e3, 3)} mJ per stroke (reactive, largely recoverable); ${F(Eloss * 1e3, 3)} mJ dissipated per cycle; none to hold a position or a passive shunt state`], ['Response time', `${F(tResp * 1e6, 3)} µs`], ['Loss', `dielectric loss tan δ ≈ 0.02; resistive shunt can add loss factor up to ${F(etaMax, 3)}`], ['Reversibility', 'Fully reversible within the field limit'], ['Hysteresis', '10–15 % of stroke in open-loop voltage drive (ferroelectric domain switching) plus creep; charge drive or feedback removes most of it'], ['Cycle durability', 'Typically above 10⁹ cycles at rated field in dry conditions; humidity-driven electrode migration and tensile stress are the usual failure causes']])],
    equations: [
      ...formulaEq('PbZr0.52Ti0.48O3', 'lead zirconate titanate near the morphotropic phase boundary (typical actuator ceramic)'),
      eqn('constitutive', 'Piezoelectric relations (33 mode)', 'S₃ = s₃₃ᴱ T₃ + d₃₃ E₃;   D₃ = d₃₃ T₃ + ε₃₃ᵀ E₃'),
      eqn('formula', 'Stack actuator', 'Δ = n d₃₃ V;  F_b = k_sc Δ;  k_sc = Y A/L'),
      eqn('formula', 'Shunt stiffness tuning', 'k = k_sc[1 + k₃₃²/(1 − k₃₃²)·C_S/(C_S + C_sh)];  k_oc = k_sc/(1 − k₃₃²)'),
      eqn('formula', 'Resistive shunt damping', 'η_max = k²/(2√(1 − k²)) at R = √(1 − k²)/(ωC_S)'),
    ],
    checks: [
      chk('Electric field within material limit', field <= 2e6, `${F(field / 1e6, 3)} MV/m`),
      chk('Shunted system stable', stableNC && kSh > 0, `1 + C_sh/C_S = ${F(1 + al, 3)}`),
      chk('Coupling coefficient admissible', p.k33 > 0 && p.k33 < 1, `k₃₃ = ${F(p.k33)}`),
      chk('Stiffness ordering k_sc ≤ k ≤ k_oc for passive shunt', al < 0 || (kSh >= ksc * (1 - 1e-12) && kSh <= koc * (1 + 1e-12)), `${F(ksc, 4)} ≤ ${F(kSh, 4)} ≤ ${F(koc, 4)} N/m`),
      chk('Lumped model valid (below stack resonance)', fsh <= 0.3 * fAct, `f/f_stack = ${F(fsh / fAct, 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', range < 0.05 ? 'Passive shunt tuning gives under 5 % frequency range with this coupling; use a higher-k material (single crystal), or a negative-capacitance shunt with a stability margin.' : `Passive short/open switching already gives ${F(100 * range, 3)} % tuning without power; a capacitor bank provides intermediate states.`, 'Coupling-coefficient limit on stiffness change.'),
      rec('medium', 'For positioning, use charge control or closed-loop strain feedback to remove the 10–15 % hysteresis and creep of voltage-driven stacks.', 'Ferroelectric hysteresis.'),
      rec('medium', `Keep the stack under compressive preload of about ${F(0.1 * Fb, 3)}–${F(0.2 * Fb, 3)} N (10–20 % of blocked force) to avoid tensile cracking in dynamic operation.`, 'Brittle multilayer ceramics tolerate little tension.'),
      rec('low', 'Use energy-recovery drivers: most of the ½CV² per stroke is reactive and can be returned to the supply. Lead-based PZT is restricted under substance regulations — consider lead-free ceramics where the lower d₃₃ is acceptable.', 'Reactive drive energy and material compliance.', true),
    ],
    publish: { 'eff.E': kSh * Ls / A, 'res.f0': fsh, 'tune.range': range, 'tune.energy_J': Est, 'tune.response_s': tResp, 'mech.stiffness': kSh },
    numerics: { method: 'Closed-form linear piezoelectric stack and shunt relations', dof: 1, steps: 60, iterations: 1, residual: 0, notes: 'Quasi-static stiffness; ideal lossless shunt capacitor; stack density 7600 kg/m³ assumed for its own resonance.' },
  };
}

// ---------- liquid crystal ----------
function lcTune(p) {
  const warnings = [], ne = Math.max(p.n_e, p.n_o + 0.001), no = p.n_o, K = p.K11 * 1e-12, d = p.d_lc * 1e-6, lam = p.lam_nm * 1e-9, Vth = N.PI * Math.sqrt(K / (N.EPS0 * p.d_eps));
  const tilt = V => V <= Vth ? 0 : N.PI / 2 - 2 * Math.atan(Math.exp(-(V - Vth) / Vth)), nEff = th => no * ne / Math.sqrt(ne * ne * Math.sin(th) ** 2 + no * no * Math.cos(th) ** 2);
  const th = tilt(p.V_lc), nV = nEff(th), dn = ne - nV, dnMax = ne - no, phase = 2 * N.PI * dn * d / lam, nAvg = 0.5 * (ne + no), shift = p.overlap * dn / nAvg, f0 = N.C0 / lam, dLam = -shift * lam, tOff = p.gamma1 * d * d / (K * N.PI ** 2), tOn = p.V_lc > Vth * 1.001 ? tOff / ((p.V_lc / Vth) ** 2 - 1) : tOff * 1000;
  const epsPar = 5 + p.d_eps, Carea = N.EPS0 * epsPar / d, Earea = 0.5 * Carea * p.V_lc ** 2, Parea = 2 * N.PI * 1000 * Carea * p.V_lc ** 2, usedFrac = dnMax > 0 ? dn / dnMax : 0;
  if (p.n_e <= p.n_o) warnings.push('Extraordinary index must exceed the ordinary index for a positive-birefringence nematic; a minimal birefringence of 0.001 is used.');
  if (p.V_lc <= Vth) warnings.push(`Drive voltage ${F(p.V_lc)} V is below the Fréedericksz threshold ${F(Vth, 3)} V: the director does not reorient and no tuning occurs.`);
  if (lam > 1e-4) warnings.push('At millimetre and longer wavelengths liquid-crystal birefringence is retained but absorption (tan δ ≈ 0.01–0.03) limits the resonator Q.');
  const vs = N.linspace(0, Math.max(4 * Vth, p.V_lc, 1), 120);
  return {
    summary: `Nematic liquid-crystal tuning (Δn = ${F(dnMax, 3)}, ${F(p.d_lc)} µm cell): at ${F(p.V_lc)} V rms (threshold ${F(Vth, 3)} V) the effective index falls from ${F(ne, 4)} to ${F(nV, 4)}, giving a phase change of ${F(phase / N.PI, 3)}π and a resonance blue shift of ${F(100 * shift, 3)} % (${F(dLam * 1e9, 4)} nm in wavelength) for a mode overlap of ${F(p.overlap)}. Response times: ${F(tOn * 1e3, 3)} ms on, ${F(tOff * 1e3, 3)} ms off.`,
    kpis: [
      kpi('V_th', 'Fréedericksz threshold', Vth, 'V', 'π√(K₁₁/(ε₀Δε)): minimum rms voltage for director reorientation in a planar cell, independent of thickness.', p.V_lc > Vth),
      kpi('n_eff', 'Effective index at drive voltage', nV, '–', `Index seen by light polarised along the rubbing direction for a mid-layer tilt of ${F(th * 180 / N.PI, 3)}° (single-angle approximation of the director profile).`, null),
      kpi('dn', 'Index change', -dn, '–', `Change from the zero-field value n_e; ${F(100 * usedFrac, 3)} % of the full birefringence ${F(dnMax, 3)} is used.`, usedFrac > 0.5),
      kpi('phase', 'Phase modulation', phase, 'rad', '2πΔn·d/λ for a single pass through the cell; a full 2π is needed for arbitrary beam steering.', phase >= 2 * N.PI),
      kpi('shift', 'Relative resonance shift', 100 * shift, '%', 'First-order perturbation Δf/f ≈ −Γ·δn/n̄ with δn = n(V) − n_e ≤ 0: the resonance moves to higher frequency (blue shift) as the index falls with voltage.', Math.abs(shift) > 0.005),
      kpi('d_lambda', 'Resonance wavelength shift', dLam, 'm', 'Equivalent wavelength shift of the tuned resonance.', null),
      kpi('t_on', 'Turn-on time', tOn, 's', 'τ_off/((V/V_th)² − 1): overdrive speeds the field-driven reorientation.', null),
      kpi('t_off', 'Relaxation time', tOff, 's', 'γ₁d²/(K₁₁π²): elastic relaxation after the field is removed, the speed-limiting step, proportional to the square of the cell gap.', tOff < 0.05),
      kpi('E_area', 'Stored energy per unit area', Earea, 'J/m²', `½CV² with C/A = ε₀ε∥/d (ε∥ ≈ ${F(epsPar, 3)} assumed). The 1 kHz AC drive circulates about ${F(Parea, 3)} W/m² of reactive power, of which only a small fraction is dissipated.`, null),
    ],
    series: [
      { key: 'n_V', title: 'Effective index versus voltage', xlabel: 'rms voltage (V)', ylabel: 'effective index', lines: [line('n(V)', vs, vs.map(v => nEff(tilt(v))))], explain: 'Electro-optic tuning curve: no change below threshold, then a fall from n_e towards n_o that saturates at a few times the threshold voltage.' },
      { key: 'shift_V', title: 'Resonance shift versus voltage', xlabel: 'rms voltage (V)', ylabel: 'relative shift (%)', lines: [line('Δf/f', vs, vs.map(v => 100 * p.overlap * (ne - nEff(tilt(v))) / nAvg))], explain: 'Predicted fractional shift of the metasurface resonance. Its saturation value is set by birefringence times mode overlap.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Field-induced reorientation of the nematic director (Fréedericksz transition) changes the index seen by the resonant mode'], ['Energy demand', `${F(Earea * 1e6, 3)} µJ/m² stored; AC holding drive is almost purely reactive (≈ µW/cm² dissipated)`], ['Response time', `on ${F(tOn * 1e3, 3)} ms, off ${F(tOff * 1e3, 3)} ms (scales with d²)`], ['Loss', 'Low in the visible and near infrared; scattering and absorption grow towards THz and microwave bands'], ['Reversibility', 'Fully reversible and continuous'], ['Hysteresis', 'Negligible for nematics with stable alignment layers; anchoring memory can appear in very thin cells'], ['Cycle durability', 'Display-grade materials sustain > 10⁹ cycles with AC drive; DC bias causes ion migration and image sticking']])],
    equations: [
      eqn('formula', 'Fréedericksz threshold', 'V_th = π√(K₁₁/(ε₀Δε))'),
      eqn('constitutive', 'Index ellipsoid', 'n(θ) = n_o n_e / √(n_e² sin²θ + n_o² cos²θ)'),
      eqn('formula', 'Director tilt (single-angle approximation)', 'θ(V) = π/2 − 2·atan{exp[−(V − V_th)/V_th]},  V > V_th'),
      eqn('formula', 'Dynamics', 'τ_off = γ₁d²/(K₁₁π²);  τ_on = τ_off/[(V/V_th)² − 1]'),
      eqn('formula', 'Resonance perturbation', 'Δf/f ≈ −Γ·[n(V) − n_e]/n̄'),
    ],
    checks: [
      chk('Drive above threshold', p.V_lc > Vth, `V/V_th = ${F(p.V_lc / Vth, 3)}`),
      chk('Index between n_o and n_e', nV >= no - 1e-12 && nV <= ne + 1e-12, `${F(no, 4)} ≤ ${F(nV, 4)} ≤ ${F(ne, 4)}`),
      chk('Positive birefringence', p.n_e > p.n_o, `Δn = ${F(p.n_e - p.n_o, 3)}`),
      chk('Relaxation faster than 50 ms', tOff < 0.05, `τ_off = ${F(tOff * 1e3, 3)} ms`),
    ],
    warnings,
    recommendations: [
      rec('high', tOff > 0.05 ? `Relaxation of ${F(tOff * 1e3, 3)} ms is slow: halve the cell gap for a fourfold speed-up, or use a dual-frequency or polymer-network liquid crystal.` : 'Response is display-grade; thinner cells would be faster still at the cost of phase range.', 'τ_off ∝ d².'),
      rec('medium', usedFrac < 0.9 ? `Only ${F(100 * usedFrac, 3)} % of the birefringence is used at ${F(p.V_lc)} V; about ${F(4 * Vth, 3)} V reaches saturation.` : 'The drive saturates the index change; higher voltage only shortens the turn-on time.', 'Electro-optic saturation.'),
      rec('medium', `Increase the mode overlap (now ${F(p.overlap)}) by infiltrating the resonator gaps: the shift scales linearly with the fraction of field energy in the liquid crystal.`, 'First-order perturbation theory.'),
      rec('low', 'Liquid-crystal tuning is field-driven and consumes microwatts per square centimetre, making it the lowest-power continuous mechanism for large optical apertures; seal cells properly, as liquid-crystal mixtures require controlled end-of-life treatment.', 'Reactive AC drive.', true),
    ],
    publish: { 'eff.n': nV, 'res.f0': f0 * (1 + shift), 'tune.range': Math.abs(shift), 'tune.response_s': Math.max(tOn, tOff), 'tune.energy_J_m2': Earea },
    numerics: { method: 'Closed-form Fréedericksz threshold, index-ellipsoid and relaxation-time relations', dof: 1, steps: 120, iterations: 1, residual: 0, notes: 'Single effective tilt angle instead of the full director profile; strong anchoring.' },
  };
}

// ---------- electrostatic MEMS ----------
function mems(p) {
  const warnings = [], k = p.k_m, A = p.A_m * 1e-12, g0 = p.g0 * 1e-6, gd = p.t_d * 1e-9 / p.eps_d, ge0 = g0 + gd, Vpi = Math.sqrt(8 * k * ge0 ** 3 / (27 * N.EPS0 * A)), Vh = Math.sqrt(2 * k * g0 * gd * gd / (N.EPS0 * A)), w0 = 2 * N.PI * p.f_mech;
  const xMax = Math.min(ge0 / 3, g0), gapAt = V => { if (V >= Vpi) return 0; const fn = x => k * x - N.EPS0 * A * V * V / (2 * (ge0 - x) ** 2); try { return g0 - N.brent(fn, 0, xMax); } catch { return g0; } };
  const pulled = p.V_m >= Vpi, g = gapAt(p.V_m), x = g0 - g, ge = g + gd, C0 = N.EPS0 * A / ge0, Cv = N.EPS0 * A / ge, Cdown = N.EPS0 * A / gd, ke = pulled ? k : N.EPS0 * A * p.V_m ** 2 / ge ** 3, keff = k - ke, fV = pulled ? 0 : p.f_mech * Math.sqrt(Math.max(keff / k, 0));
  const Cratio = Cv / C0, Esw = 0.5 * (pulled ? Cdown : Cv) * p.V_m ** 2, tSw = pulled ? 3.67 * Vpi / (p.V_m * w0) : Math.max(1 / (2 * p.f_mech), p.Q_m / (N.PI * p.f_mech)), Edown = p.V_m / (p.t_d * 1e-9), onOff = Cdown / C0, mass = k / (w0 * w0);
  if (pulled) warnings.push(`Applied voltage exceeds the pull-in voltage ${F(Vpi, 4)} V: the plate collapses onto the dielectric (switch closed). It releases only below the hold-down voltage ${F(Vh, 3)} V.`);
  if (pulled && Edown > 3e8) warnings.push(`Field in the dielectric when closed is ${F(Edown / 1e8, 3)}×10⁸ V/m: dielectric charging and breakdown limit lifetime; reduce the hold voltage after actuation.`);
  if (gd > g0 / 3) warnings.push('The dielectric contributes a large part of the electrical gap; pull-in occurs after more than one third of the air gap has closed.');
  const vs = N.linspace(0, 0.999 * Vpi, 100);
  return {
    summary: `Electrostatic parallel-plate actuator (${F(p.A_m)} µm², ${F(p.g0)} µm gap, ${F(k)} N/m): pull-in at ${F(Vpi, 4)} V, hold-down at ${F(Vh, 3)} V. At ${F(p.V_m)} V the plate ${pulled ? 'is pulled in (switch closed), capacitance ratio ' + F(onOff, 3) : `deflects ${F(x * 1e9, 4)} nm, raising the capacitance by ${F(100 * (Cratio - 1), 3)} % and softening the resonance from ${F(p.f_mech, 4)} to ${F(fV, 4)} Hz`}. Switching energy ${F(Esw * 1e12, 3)} pJ, response time ${F(tSw * 1e6, 3)} µs.`,
    kpis: [
      kpi('V_pi', 'Pull-in voltage', Vpi, 'V', '√(8k·g³/(27ε₀A)) with the electrical gap g = g₀ + t_d/ε_d: beyond it the electrostatic force grows faster than the spring force and the plate snaps down.', p.V_m < Vpi),
      kpi('V_hold', 'Hold-down (release) voltage', Vh, 'V', 'Voltage below which the spring pulls the collapsed plate off the dielectric; the difference from the pull-in voltage is the electromechanical hysteresis window.', null),
      kpi('x', 'Static deflection', pulled ? g0 : x, 'm', pulled ? 'Plate in contact with the dielectric: full gap closed.' : `Stable equilibrium of k·x = ε₀AV²/(2(g − x)²); the stable travel limit is ${F(xMax * 1e9, 4)} nm (one third of the electrical gap).`, !pulled),
      kpi('C_ratio', 'Capacitance tuning ratio', pulled ? onOff : Cratio, '–', pulled ? 'Down-state over up-state capacitance of the switch.' : 'Capacitance at bias over zero-bias capacitance; an analogue parallel-plate varactor cannot exceed 1.5 before pull-in.', null),
      kpi('f_V', 'Resonance at bias', fV, 'Hz', pulled ? 'No free resonance in the collapsed state; reported as 0.' : `f₀√(1 − k_e/k) with electrostatic spring k_e = ε₀AV²/g³ = ${F(ke, 3)} N/m (a negative stiffness).`, !pulled),
      kpi('k_eff', 'Effective stiffness', pulled ? 0 : keff, 'N/m', 'Mechanical stiffness minus the electrostatic negative stiffness; it reaches zero exactly at pull-in.', keff > 0 && !pulled),
      kpi('t_sw', 'Response time', tSw, 's', pulled ? 'Switching-time estimate 3.67·V_pi/(V·ω₀) for an inertia-limited plate.' : 'Mechanical settling time of the analogue-tuned plate, about Q/(πf₀) (or half a period if over-damped).', null),
      kpi('E_sw', 'Energy per actuation', Esw, 'J', '½CV² stored on the actuator capacitance; no static current flows, so holding power is limited to leakage.', null),
      kpi('m_eff', 'Effective mass', mass, 'kg', 'k/ω₀² implied by the stiffness and resonance frequency.', null),
    ],
    series: [
      { key: 'x_V', title: 'Deflection versus voltage', xlabel: 'voltage (V)', ylabel: 'deflection (nm)', lines: [line('stable branch', vs, vs.map(v => (g0 - gapAt(v)) * 1e9)), line('travel limit', [0, Vpi], [xMax * 1e9, xMax * 1e9])], explain: 'Stable equilibrium deflection up to pull-in. The curve becomes vertical at the pull-in voltage, where control of position is lost.' },
      { key: 'f_V', title: 'Resonance and capacitance versus voltage', xlabel: 'voltage (V)', ylabel: 'f/f₀ and C/C₀', lines: [line('f/f₀', vs, vs.map(v => { const gg = gapAt(v) + gd; return Math.sqrt(Math.max(1 - N.EPS0 * A * v * v / gg ** 3 / k, 0)); })), line('C/C₀', vs, vs.map(v => ge0 / (gapAt(v) + gd)))], explain: 'Electrostatic spring softening lowers the resonance to zero at pull-in while the capacitance rises by at most 50 %: the two usable analogue tuning effects.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Electrostatic attraction between a suspended and a fixed electrode: gap tuning and negative electrostatic stiffness'], ['Energy demand', `${F(Esw * 1e12, 3)} pJ per actuation; near-zero holding power (leakage only)`], ['Response time', `${F(tSw * 1e6, 3)} µs (mechanical, set by resonance ${F(p.f_mech, 3)} Hz and damping)`], ['Loss', 'Very low electrical loss (air gap, metal plates); mechanical Q limited by squeeze-film damping'], ['Reversibility', 'Reversible; analogue below pull-in, bistable (digital) above it'], ['Hysteresis', `Pull-in / release window ${F(Vpi, 3)} V / ${F(Vh, 3)} V`], ['Cycle durability', 'Commonly 10⁹–10¹¹ cycles when switched without RF power; dielectric charging and stiction are the dominant failure modes']])],
    equations: [
      eqn('governing', 'Static equilibrium', 'k x = ε₀ A V² / [2 (g₀ + t_d/ε_d − x)²]'),
      eqn('formula', 'Pull-in and hold-down', 'V_pi = √[8k g³/(27 ε₀ A)],  x_pi = g/3;   V_h = √[2k g₀ (t_d/ε_d)²/(ε₀A)]'),
      eqn('formula', 'Electrostatic spring softening', 'k_eff = k − ε₀AV²/g³;  f = f₀√(k_eff/k)'),
      eqn('formula', 'Switching time', 't_s ≈ 3.67 V_pi/(V_s ω₀)'),
    ],
    checks: [
      chk('Operating below pull-in (analogue regime)', !pulled, `V/V_pi = ${F(p.V_m / Vpi, 3)}`),
      chk('Effective stiffness positive (stable)', pulled || keff > 0, `k_eff = ${F(pulled ? 0 : keff, 4)} N/m`),
      chk('Deflection within stable travel', pulled || x <= xMax * (1 + 1e-9), `${F(x * 1e9, 4)} nm of ${F(xMax * 1e9, 4)} nm`),
      chk('Hold-down below pull-in (hysteresis window exists)', Vh < Vpi, `${F(Vh, 3)} V < ${F(Vpi, 3)} V`),
      chk('Dielectric field acceptable when closed', !pulled || Edown <= 3e8, `${F(Edown, 3)} V/m`),
    ],
    warnings,
    recommendations: [
      rec('high', pulled ? `After pull-in reduce the voltage towards ${F(1.3 * Vh, 3)} V (30 % above hold-down): this lowers dielectric charging, the main lifetime limiter.` : `Analogue tuning is limited to ${F(100 * (Cratio - 1), 3)} % capacitance change here and 50 % at most; use separate actuation and signal electrodes or a comb drive for a wider range.`, 'Pull-in instability at one third of the gap.'),
      rec('medium', `Bias noise converts to frequency noise through spring softening; at this operating point the resonance shifts by ${F(pulled ? 0 : 100 * (1 - fV / p.f_mech), 3)} %. Stabilise the bias supply accordingly.`, 'Electrostatic negative stiffness.'),
      rec('medium', 'Package hermetically: humidity promotes stiction and dielectric charging, and pressure sets the damping (and thus the settling time).', 'Typical MEMS failure mechanisms.'),
      rec('low', 'Electrostatic actuation holds its state with essentially no power, so large programmable arrays consume energy only when reconfigured.', 'Capacitive actuation energy.', true),
    ],
    publish: { 'res.f0': fV, 'tune.range': pulled ? 1 : 1 - fV / p.f_mech, 'tune.energy_J': Esw, 'tune.response_s': tSw, 'mech.stiffness': pulled ? 0 : keff },
    numerics: { method: 'Parallel-plate electrostatic equilibrium by Brent root finding; closed-form pull-in relations', dof: 1, steps: 100, iterations: 100, residual: 0, notes: 'Rigid plate, no fringing fields, linear suspension.' },
  };
}

// ---------- ferrite ----------
function ferrite(p) {
  const warnings = [], mu0 = N.MU0, wm = GAMMA * mu0 * p.Ms, al = p.alpha_g, w0of = H => GAMMA * mu0 * H, fresOf = H => GAMMA * mu0 * Math.sqrt(H * (H + p.Ms)) / (2 * N.PI), farOf = H => GAMMA * mu0 * (H + p.Ms) / (2 * N.PI);
  const muEff = (f, H) => { const w = 2 * N.PI * f, C = N.C, w0c = { re: w0of(H), im: al * w }, den = C.sub(C.mul(w0c, w0c), { re: w * w, im: 0 }), mu = C.add({ re: 1, im: 0 }, C.div(C.scale(w0c, wm), den)), ka = C.div({ re: w * wm, im: 0 }, den); return { mu, ka, eff: C.div(C.sub(C.mul(mu, mu), C.mul(ka, ka)), mu) }; };
  const fres = fresOf(p.H_b), far = farOf(p.H_b), f0 = w0of(p.H_b) / (2 * N.PI), m = muEff(p.f_op, p.H_b), neg = m.eff.re < 0, tanM = Math.abs(m.eff.im) / Math.max(Math.abs(m.eff.re), 1e-12), sens = (fresOf(p.H_b * 1.001) - fresOf(p.H_b * 0.999)) / (0.002 * p.H_b), dH = 2 * al * 2 * N.PI * fres / (GAMMA * mu0), Qf = 1 / (2 * al);
  const gap = p.gap_mm * 1e-3, I = p.H_b * gap / p.N_coil, Pc = I * I * p.R_coil, Lc = mu0 * p.N_coil ** 2 * p.A_pole * 1e-4 / gap, tauC = Lc / p.R_coil, Emag = 0.5 * Lc * I * I, Hneed = (() => { try { return N.brent(H => fresOf(H) - p.f_op, 1, 5e7); } catch { return 0; } })();
  const inBand = p.f_op > fres && p.f_op < far, saturated = p.H_b > 0.05 * p.Ms;
  if (!saturated) warnings.push('The bias field is small compared with the magnetisation: the ferrite may not be saturated, and domain (low-field) losses and hysteresis, not included here, will dominate.');
  if (Math.abs(p.f_op - fres) < 3 * al * fres + 1e-3 * fres) warnings.push('Operating within a few linewidths of ferromagnetic resonance: absorption is very high and the effective permeability varies rapidly with field and temperature.');
  if (Pc > 100) warnings.push(`Bias coil dissipates ${F(Pc, 3)} W: use a permanent-magnet bias with a small trim coil instead.`);
  const fsw = N.linspace(Math.max(0.3 * fres, 1e6), 1.5 * far, 300), mm = fsw.map(f => muEff(f, p.H_b).eff), Hs = N.linspace(0.1 * p.H_b, 3 * p.H_b, 80);
  return {
    summary: `Transversely magnetised ferrite (M_s = ${F(p.Ms, 3)} A/m) biased at ${F(p.H_b, 3)} A/m: ferromagnetic resonance at ${F(fres / 1e9, 4)} GHz and anti-resonance at ${F(far / 1e9, 4)} GHz, between which the effective permeability is negative. At ${F(p.f_op / 1e9, 4)} GHz μ_eff = ${F(m.eff.re, 4)} ${m.eff.im < 0 ? '−' : '+'} j${F(Math.abs(m.eff.im), 3)}${neg ? ' (negative-permeability band)' : ''}. Tuning sensitivity ${F(sens * 79.577 / 1e6, 4)} MHz/Oe; the bias coil needs ${F(I, 3)} A and ${F(Pc, 3)} W with a ${F(tauC * 1e3, 3)} ms time constant.`,
    kpis: [
      kpi('f_res', 'Ferromagnetic resonance', fres, 'Hz', '(γμ₀/2π)√(H(H + M_s)): pole of the effective permeability for a wave propagating perpendicular to the bias (Kittel frequency of an in-plane-magnetised slab).', null),
      kpi('f_ar', 'Anti-resonance', far, 'Hz', '(γμ₀/2π)(H + M_s): zero of the effective permeability. Between resonance and anti-resonance μ_eff is negative and waves are evanescent.', null),
      kpi('mu_re', 'Re μ_eff at operating frequency', m.eff.re, '–', 'Real part of (μ² − κ²)/μ from the Polder tensor with Gilbert damping; negative values provide the magnetic half of a negative-index medium.', null),
      kpi('mu_im', 'Im μ_eff at operating frequency', m.eff.im, '–', `Loss part (negative for the e^{jωt} convention); magnetic loss tangent ${F(tanM, 3)}.`, m.eff.im <= 1e-9),
      kpi('sens', 'Field tuning sensitivity', sens, 'Hz·m/A', `Slope of the resonance with bias field, equal to ${F(sens * 79.577 / 1e6, 4)} MHz/Oe.`, null),
      kpi('H_need', 'Bias field for resonance at operating frequency', Hneed, 'A/m', 'Internal field that would place the ferromagnetic resonance exactly at the operating frequency.', null),
      kpi('dH', 'Resonance linewidth', dH, 'A/m', `Full linewidth 2αω/(γμ₀) = ${F(dH / 79.577, 3)} Oe; the resonance quality factor is about 1/(2α) = ${F(Qf, 3)}.`, null),
      kpi('I_coil', 'Bias coil current', I, 'A', 'H·gap/N for an electromagnet with negligible core reluctance.', null),
      kpi('P_coil', 'Bias holding power', Pc, 'W', 'Ohmic dissipation I²R in the coil: the continuous cost of electromagnet tuning.', Pc < 10),
      kpi('tau', 'Tuning response time', tauC, 's', `Coil time constant L/R with L = μ₀N²A/gap = ${F(Lc, 3)} H; the spin system itself responds in nanoseconds.`, null),
      kpi('E_mag', 'Stored magnetic energy', Emag, 'J', '½LI² in the bias magnet, supplied at each field change.', null),
    ],
    series: [
      { key: 'mu_f', title: 'Effective permeability versus frequency', xlabel: 'frequency (GHz)', ylabel: 'μ_eff', lines: [line('Re μ_eff', fsw.map(v => v / 1e9), mm.map(v => N.clamp(v.re, -30, 30))), line('Im μ_eff', fsw.map(v => v / 1e9), mm.map(v => N.clamp(v.im, -30, 30)))], explain: 'Dispersion of the effective permeability (clipped at ±30). The real part is negative between resonance and anti-resonance; the imaginary part peaks at resonance.' },
      { key: 'f_H', title: 'Resonance and anti-resonance versus bias field', xlabel: 'bias field (kA/m)', ylabel: 'frequency (GHz)', lines: [line('resonance', Hs.map(v => v / 1e3), Hs.map(h => fresOf(h) / 1e9)), line('anti-resonance', Hs.map(v => v / 1e3), Hs.map(h => farOf(h) / 1e9)), line('operating frequency', [Hs[0] / 1e3, Hs[Hs.length - 1] / 1e3], [p.f_op / 1e9, p.f_op / 1e9])], explain: 'Magnetic tuning map. The negative-permeability band lies between the two curves and moves up with bias field.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Precession of the saturated magnetisation about the bias field (gyromagnetic resonance); the bias field sets the resonance frequency'], ['Energy demand', `${F(Pc, 3)} W continuous in the coil (zero with a permanent magnet); ${F(Emag, 3)} J per field change`], ['Response time', `${F(tauC * 1e3, 3)} ms (coil L/R); intrinsic spin response in nanoseconds`], ['Loss', `magnetic loss tangent ${F(tanM, 3)} at the operating point; linewidth ${F(dH / 79.577, 3)} Oe`], ['Reversibility', 'Fully reversible and continuous above saturation'], ['Hysteresis', 'Negligible when saturated; coercivity and remanence matter below saturation (latching devices exploit this)'], ['Cycle durability', 'No wear-out mechanism in the ferrite; limits are coil heating and temperature drift of M_s']])],
    equations: [
      ...formulaEq('Y3Fe5O12', 'yttrium iron garnet, the reference low-loss microwave ferrite'),
      eqn('constitutive', 'Polder permeability tensor', 'μ = 1 + ω₀ω_m/(ω₀² − ω²),  κ = ωω_m/(ω₀² − ω²),  ω₀ → ω₀ + jαω;  ω₀ = γμ₀H, ω_m = γμ₀M_s'),
      eqn('formula', 'Effective permeability (wave ⟂ bias)', 'μ_eff = (μ² − κ²)/μ'),
      eqn('formula', 'Resonance and anti-resonance', 'ω_res = γμ₀√(H(H + M_s));   ω_ar = γμ₀(H + M_s)'),
      eqn('formula', 'Bias electromagnet', 'H ≈ N I / gap;  P = I²R;  τ = L/R'),
    ],
    checks: [
      chk('Passivity (Im μ_eff ≤ 0)', mm.every(v => v.im <= 1e-9), `Im μ_eff at operating point ${F(m.eff.im, 3)}`),
      chk('Negative-μ band consistent with sign of Re μ_eff', inBand === neg || Math.abs(p.f_op - fres) < 5 * al * fres || Math.abs(p.f_op - far) < 5 * al * far, `operating point ${inBand ? 'inside' : 'outside'} [f_res, f_ar], Re μ_eff = ${F(m.eff.re, 3)}`),
      chk('Ferrite magnetically saturated', saturated, `H/M_s = ${F(p.H_b / p.Ms, 3)}`),
      chk('Static limit μ_eff(ω→0) = 1 + M_s/H', Math.abs(muEff(1, p.H_b).eff.re - (1 + p.Ms / p.H_b)) < 1e-3 * (1 + p.Ms / p.H_b), `${F(muEff(1, p.H_b).eff.re, 5)} vs ${F(1 + p.Ms / p.H_b, 5)}`),
      chk('Coil power moderate', Pc < 100, `${F(Pc, 3)} W`),
    ],
    warnings,
    recommendations: [
      rec('high', neg ? `The operating point lies in the negative-permeability band (${F(fres / 1e9, 3)}–${F(far / 1e9, 3)} GHz): combine with a negative-permittivity wire medium for a tunable negative-index band; loss tangent here is ${F(tanM, 3)}.` : `To obtain negative permeability at ${F(p.f_op / 1e9, 3)} GHz set the bias between ${F(Math.max(p.f_op * 2 * N.PI / (GAMMA * mu0) - p.Ms, 0), 3)} and ${F(Hneed, 3)} A/m.`, 'Position relative to resonance and anti-resonance.'),
      rec('medium', Pc > 10 ? 'Provide the main bias with permanent magnets and use the coil only for trimming: the present holding power is high.' : 'Coil power is modest; a soft-iron yoke reduces it further.', `Coil dissipation ${F(Pc, 3)} W.`),
      rec('medium', 'Stabilise temperature or compensate the bias: saturation magnetisation falls with temperature, shifting the band edges.', 'Temperature dependence of M_s.'),
      rec('low', 'A latched (remanent-state) ferrite or permanent-magnet bias holds the tuned state without power, removing the continuous coil loss.', 'Holding power of electromagnet bias.', true),
    ],
    publish: { 'eff.mu': m.eff.re, 'res.f0': fres, 'res.bandgap_lo': fres, 'res.bandgap_hi': far, 'res.Q': Qf, 'tune.response_s': tauC, 'tune.hold_power_W': Pc },
    numerics: { method: 'Closed-form Polder tensor with Gilbert damping (complex arithmetic); Brent inversion for the bias field', dof: 1, steps: 300, iterations: 1, residual: 0, notes: 'Uniformly saturated ferrite, internal field given, demagnetising and anisotropy fields neglected.' },
  };
}

// ---------- bistable array ----------
function bistable(p) {
  const warnings = [], kb = p.k_bar, w = p.w_b * 1e-3, h = p.h_b * 1e-3, L0 = Math.hypot(w, h), nU = Math.round(p.n_units), rnd = N.rng(p.seed);
  const Fu = x => { const L = Math.hypot(w, h - x); return 2 * kb * (L0 - L) * (h - x) / L + p.k_v * x; }, dF = x => (Fu(x + 1e-6 * h) - Fu(x - 1e-6 * h)) / (2e-6 * h);
  const xPk = N.goldenSection(x => -Fu(x), 0, h, 1e-10 * h), xVal = N.goldenSection(Fu, h, 2.2 * h, 1e-10 * h), Fpk = Fu(xPk), Fval = Fu(xVal), bist = Fval < 0;
  let xU = h, x2 = 2 * h; try { xU = N.brent(Fu, xPk, xVal); } catch { xU = h; } if (bist) { try { x2 = N.brent(Fu, xVal, 3 * h); } catch { x2 = 2 * h; } }
  const integ = (a, b) => { const xs = N.linspace(a, b, 200); return N.trapz(xs, xs.map(Fu)); }, Ebar = integ(0, xU), Eback = bist ? -integ(xU, x2) : 0, k1 = dF(0), k2 = bist ? dF(x2) : k1, kNeg = dF(h), strain = (L0 - w) / L0;
  const sc = Array.from({ length: nU }, () => Math.max(0.3, 1 + p.var_pct / 100 * rnd.normal())), Fhi = 1.15 * N.maxOf(sc) * Fpk, Flo = bist ? 1.15 * N.maxOf(sc) * Fval : 0, nStep = 60;
  const inv = (F, s, st) => { const fn = x => s * Fu(x) - F; try { if (st === 0) return N.brent(fn, -3 * h, xPk); let hi = 3 * h; for (let i = 0; i < 20 && fn(hi) < 0; i++) hi *= 1.5; return N.brent(fn, xVal, hi); } catch { return st === 0 ? (F > 0 ? xPk : 0) : xVal; } };
  const state = new Uint8Array(nU), Fpath = [...N.linspace(0, Fhi, nStep), ...N.linspace(Fhi, 0, 30), ...(bist ? [...N.linspace(0, Flo, 30).slice(1), ...N.linspace(Flo, 0, 20).slice(1)] : [])], Xp = [], map = new Float64Array(Fpath.length * nU); let nOnMax = 0, retained = 0;
  Fpath.forEach((F, q) => { let X = 0; for (let i = 0; i < nU; i++) { if (state[i] === 0 && F > sc[i] * Fpk) state[i] = 1; else if (state[i] === 1 && F < sc[i] * Fval) state[i] = 0; X += inv(F, sc[i], state[i]); map[q + Fpath.length * i] = state[i]; } Xp.push(X); nOnMax = Math.max(nOnMax, N.sum(state)); if (q === nStep + 29) retained = N.sum(state); });
  let loop = 0; for (let q = 1; q < Fpath.length; q++) loop += 0.5 * (Fpath[q] + Fpath[q - 1]) * (Xp[q] - Xp[q - 1]);
  const fn1 = Math.sqrt(Math.max(k1, 1e-30) / p.m_unit) / (2 * N.PI), tSw = 1 / (4 * fn1), kChain = nOn => 1 / ((nU - nOn) / k1 + nOn / Math.max(k2, 1e-30)), kRange = bist ? Math.max(kChain(0), kChain(nU)) / Math.min(kChain(0), kChain(nU)) : 1;
  const order = sc.map((s, i) => [s, i]).sort((a, b) => a[0] - b[0]).map(v => v[1] + 1);
  if (!bist) warnings.push(`The unit is monostable (valley force ${F(Fval, 3)} N > 0): it snaps through under load but returns when unloaded, so it stores no state. Reduce the parallel spring or increase the apex height.`);
  if (strain > 0.02) warnings.push(`The bars are compressed by ${F(100 * strain, 3)} % when the unit passes through the flat configuration: beyond the elastic range of most stiff materials; use compliant flexures or elastomers.`);
  if (p.var_pct === 0 && nU > 1) warnings.push('With zero variation all units reach their threshold simultaneously: the snapping order is not defined in practice and is set by imperfections.');
  const xs = N.linspace(-0.2 * h, 2.4 * h, 160);
  return {
    summary: `Von Mises truss unit (half-span ${F(p.w_b)} mm, rise ${F(p.h_b)} mm): snap-through at ${F(Fpk, 4)} N, ${bist ? `snap-back at ${F(Fval, 4)} N — bistable with an energy barrier of ${F(Ebar * 1e3, 4)} mJ forward and ${F(Eback * 1e3, 4)} mJ back` : 'monostable (no stored state)'}; most negative tangent stiffness ${F(kNeg, 4)} N/m. A chain of ${nU} units has ${bist ? `2^${nU} = ${F(2 ** nU, 4)} stable configurations` : 'a single stable configuration'}; loading above all thresholds and unloading to zero leaves ${retained} of ${nU} units switched (non-volatile), and a full set–reset force cycle dissipates ${F(loop * 1e3, 4)} mJ.`,
    kpis: [
      kpi('F_snap', 'Snap-through force', Fpk, 'N', `Local maximum of the unit force–displacement curve at ${F(xPk * 1e3, 4)} mm: the force that programmes a unit from state 0 to state 1.`, null),
      kpi('F_back', 'Snap-back force', Fval, 'N', 'Local minimum of the force curve. A negative value means a reverse force is needed to reset the unit, i.e. the unit is bistable.', bist),
      kpi('E_bar', 'Forward energy barrier', Ebar, 'J', 'Work ∫F dx from the first stable state to the unstable equilibrium: the minimum energy to switch one unit, and its protection against accidental switching.', null),
      kpi('E_back', 'Reverse energy barrier', Eback, 'J', bist ? 'Work needed to return from the second stable state over the barrier; the difference between the two barriers is the energy stored in the switched state.' : 'No second stable state exists; reported as 0.', bist),
      kpi('k_neg', 'Most negative tangent stiffness', kNeg, 'N/m', 'Slope dF/dx at mid-travel. Negative stiffness is stable only under displacement control or when embedded in a sufficiently stiff surrounding structure.', null),
      kpi('n_states', 'Stable configurations (log₁₀)', bist ? nU * Math.log10(2) : 0, '–', bist ? `Each bistable unit stores one bit: 2^${nU} configurations and ${nU + 1} distinct stack heights for identical units.` : 'Monostable units give a single configuration.', bist),
      kpi('E_loop', 'Energy dissipated per load cycle', loop, 'J', 'Area enclosed by the force–displacement loop of the chain under quasi-static force control: energy released in the snaps (reusable energy absorption).', null),
      kpi('k_range', 'Programmable stiffness ratio', kRange, '–', `Ratio of chain tangent stiffness between all-off (${F(kChain(0), 4)} N/m) and all-on (${F(kChain(nU), 4)} N/m) configurations at zero load.`, null),
      kpi('t_sw', 'Switching time', tSw, 's', `About a quarter period of the unit resonance (${F(fn1, 4)} Hz) for the stated moving mass; inertia-limited snap.`, null),
      kpi('strain', 'Peak bar strain', strain, '–', 'Compressive strain of the inclined bars in the flat configuration, (L₀ − w)/L₀; it governs fatigue life of the flexures.', strain < 0.01),
      kpi('retained', 'Units left switched after the cycle', retained, '–', 'Number of units remaining in state 1 after loading above all thresholds and unloading to zero force: the programmed state held without power. A reverse force then resets them.', null),
    ],
    series: [
      { key: 'F_x', title: 'Unit force–displacement curve', xlabel: 'apex displacement (mm)', ylabel: 'force (N)', lines: [line('F(x)', xs.map(v => v * 1e3), xs.map(Fu)), line('zero', [xs[0] * 1e3, xs[xs.length - 1] * 1e3], [0, 0])], explain: 'Non-monotonic curve of one unit: positive stiffness, a negative-stiffness branch between the force maximum and minimum, then positive stiffness again. Zero crossings with positive slope are stable states.' },
      { key: 'chain', title: 'Chain response under a force cycle', xlabel: 'total displacement (mm)', ylabel: 'force (N)', lines: [line('load–unload path', Xp.map(v => v * 1e3), Fpath)], explain: 'Sequential snapping produces displacement jumps on loading; after unloading to zero the bistable units stay switched, and a reverse force snaps them back. The enclosed area is the energy dissipated per set–reset cycle.' },
    ],
    maps: [{ key: 'states', title: 'Programmable state map', nx: Fpath.length, ny: nU, data: map, x0: 0, x1: Fpath.length, y0: 0, y1: nU, xlabel: 'load step', ylabel: 'unit index', unit: 'state (0 or 1)', explain: 'State of every unit along the load cycle. Units switch in order of their thresholds (weakest first) and reset in reverse order of their snap-back forces.' }],
    tables: [mechTable([['Physical mechanism', 'Geometric snap-through of a shallow two-bar truss (elastic instability with two energy minima)'], ['Energy demand', `${F(Ebar * 1e3, 4)} mJ per unit to set, ${F(Eback * 1e3, 4)} mJ to reset; none to hold a state`], ['Response time', `${F(tSw * 1e3, 4)} ms per snap (inertia-limited)`], ['Loss', `${F(loop * 1e3, 4)} mJ dissipated per full cycle of the chain (useful for energy absorption, a cost for actuation)`], ['Reversibility', bist ? 'Reversible by applying the opposite force; states are non-volatile' : 'Self-resetting (monostable)'], ['Hysteresis', `Intrinsic: force window ${F(Fval, 3)} to ${F(Fpk, 3)} N`], ['Cycle durability', `Governed by flexure fatigue at ${F(100 * strain, 3)} % strain; snapping order ${order.slice(0, 8).join(', ')}${nU > 8 ? ', …' : ''} is set by manufacturing scatter`]])],
    equations: [
      eqn('constitutive', 'Von Mises truss unit', 'F(x) = 2k_b(L₀ − L)(h − x)/L + k_v x,  L = √(w² + (h − x)²),  L₀ = √(w² + h²)'),
      eqn('formula', 'Bistability condition', 'min F(x) < 0 on h < x < 2h'),
      eqn('formula', 'Energy barrier', 'E_b = ∫₀^{x_u} F dx,  F(x_u) = 0, dF/dx < 0'),
      eqn('formula', 'Series chain', 'X = Σᵢ xᵢ(F, sᵢ);   1/K = Σᵢ 1/kᵢ(sᵢ)'),
    ],
    checks: [
      chk('Units bistable', bist, `snap-back force ${F(Fval, 4)} N`),
      chk('Force curve has a negative-stiffness branch', kNeg < 0, `dF/dx at mid-travel = ${F(kNeg, 4)} N/m`),
      chk('Stable states have positive stiffness', k1 > 0 && k2 > 0, `k₁ = ${F(k1, 4)} N/m, k₂ = ${F(k2, 4)} N/m`),
      chk('Dissipated energy non-negative', loop >= -1e-9 * Math.abs(Fpk * h), `${F(loop, 3)} J per cycle`),
      chk('Bar strain within elastic range', strain < 0.02, `${F(100 * strain, 3)} %`),
      chk('All units switched at peak load', nOnMax === nU, `${nOnMax} of ${nU}`),
    ],
    warnings,
    recommendations: [
      rec('high', bist ? `Barrier of ${F(Ebar * 1e3, 3)} mJ per unit: check it against the expected vibration or impact energy so that states are not lost in service; raise h/w for a higher barrier.` : 'Make the unit bistable (remove or soften the parallel spring, or increase the rise) if non-volatile programmable states are required.', 'Energy landscape of the unit.'),
      rec('medium', p.var_pct < 3 && nU > 1 ? 'Grade the unit stiffness deliberately (a few per cent per unit) to obtain a deterministic snapping sequence instead of relying on random imperfections.' : `Unit scatter of ${F(p.var_pct)} % fixes the snapping order; tighten tolerances if individual addressing by force level is intended.`, 'Sequential snapping governed by threshold spread.'),
      rec('medium', strain > 0.01 ? `Peak bar strain ${F(100 * strain, 3)} % is high for metals and stiff polymers: use flexure hinges so that deformation is carried in bending.` : 'Bar strain is low enough for high-cycle operation in most engineering polymers and metals.', 'Strain in the flat configuration.'),
      rec('low', 'Bistable units hold their programmed state with zero power and release absorbed impact energy elastically, so the structure is reusable rather than sacrificial like crushable foams.', 'Non-volatile mechanical states and recoverable snap-through.', true),
    ],
    publish: { 'mech.snap_force': Fpk, 'mech.stiffness': kChain(retained), 'tune.energy_J': Ebar, 'tune.response_s': tSw, 'tune.range': kRange - 1, 'mech.energy_absorbed': loop },
    numerics: { method: 'Closed-form truss force law; golden-section extrema, Brent roots, quasi-static force-controlled path following of the series chain', dof: nU, steps: Fpath.length, iterations: Fpath.length * nU, residual: 0, notes: 'Rigid supports, pin joints, no inertia during the path; seeded stiffness scatter.' },
  };
}

// ---------- quasi-zero-stiffness isolator ----------
function qzs(p) {
  const warnings = [], kv = p.kq_v, kh = p.kq_h, L0 = p.L0_q * 1e-3; let a = p.a_q * 1e-3; if (a >= L0) { a = 0.999 * L0; warnings.push('Compressed length limited to the free length: the oblique springs must be pre-compressed to provide negative stiffness.'); }
  const K0 = kv - 2 * kh * (L0 / a - 1), aStar = L0 * 2 * kh / (kv + 2 * kh), k3 = kh * L0 / a ** 3, Fx = x => kv * x + 2 * kh * x * (1 - L0 / Math.hypot(x, a)), Kx = x => kv + 2 * kh * (1 - L0 * a * a / Math.pow(x * x + a * a, 1.5));
  const stable = K0 > 0, m = p.m_q, fLin = Math.sqrt(kv / m) / (2 * N.PI), fn = stable ? Math.sqrt(K0 / m) / (2 * N.PI) : 0, c = 2 * p.zeta_q * Math.sqrt(kv * m), zEff = stable ? c / (2 * Math.sqrt(K0 * m)) : 1;
  const Tr = (f, f0, z) => { if (f0 <= 0) return 1; const r = f / f0; return Math.sqrt((1 + (2 * z * r) ** 2) / ((1 - r * r) ** 2 + (2 * z * r) ** 2)); }, Tq = stable ? Tr(p.f_exc, fn, zEff) : 1, Tl = Tr(p.f_exc, fLin, p.zeta_q);
  const xRange = K0 < 0.1 * kv && stable ? Math.sqrt((0.1 * kv - K0) / (3 * k3)) : 0, dKda = 2 * kh * L0 / (a * a), tolA = 0.05 * kv / dKda, Epre = kh * (L0 - a) ** 2, Fpre = kh * (L0 - a), sag = m * 9.81 / kv, tSettle = stable ? 4 / (zEff * 2 * N.PI * fn) : 0, xBi = !stable ? a * Math.sqrt(Math.max(Math.pow(2 * kh * L0 / ((kv + 2 * kh) * a), 2) - 1, 0)) : 0;
  if (!stable) warnings.push(`Net stiffness is negative (${F(K0, 4)} N/m): the central position is unstable and the mount snaps to ±${F(xBi * 1e3, 3)} mm (bistable). Increase the compressed length above ${F(aStar * 1e3, 4)} mm.`);
  if (stable && K0 < 0.02 * kv) warnings.push(`Stiffness is within 2 % of zero: a length error of ${F(tolA * 1e6 * 0.4, 3)} µm or a payload change would make the mount unstable; keep a margin.`);
  if (sag > 0.5 * a) warnings.push('Static sag under the payload is large compared with the oblique-spring geometry; the working point leaves the low-stiffness region unless the vertical spring is preloaded to carry the weight at x = 0.');
  const xs = N.linspace(-0.8 * a, 0.8 * a, 161), fsw = N.logspace(Math.max(0.02 * fLin, 1e-3), 20 * fLin, 200), as = N.linspace(0.3 * L0, 0.999 * L0, 100);
  return {
    summary: `Quasi-zero-stiffness mount (k_v = ${F(kv)} N/m, 2 × ${F(kh)} N/m oblique springs compressed from ${F(p.L0_q)} to ${F(a * 1e3, 4)} mm): net stiffness ${F(K0, 4)} N/m (${F(100 * K0 / kv, 3)} % of the vertical spring), ${stable ? `natural frequency ${F(fn, 4)} Hz instead of ${F(fLin, 4)} Hz, and transmissibility ${F(Tq, 3)} at ${F(p.f_exc)} Hz against ${F(Tl, 3)} for the linear mount` : 'i.e. NEGATIVE — unstable'}. Exact zero stiffness occurs at ${F(aStar * 1e3, 4)} mm; the tuning tolerance for 5 % residual stiffness is ${F(tolA * 1e6, 3)} µm.`,
    kpis: [
      kpi('K0', 'Net stiffness at equilibrium', K0, 'N/m', 'k_v − 2k_h(L₀/a − 1): the pre-compressed oblique springs contribute a negative stiffness that cancels part of the load-bearing spring.', stable),
      kpi('ratio', 'Stiffness ratio K/k_v', K0 / kv, '–', 'Residual fraction of the linear-mount stiffness; values of a few per cent give a tenfold lower natural frequency at unchanged static load capacity.', K0 > 0 && K0 < 0.2 * kv),
      kpi('f_n', 'Natural frequency', fn, 'Hz', stable ? `√(K/m)/2π for the ${F(m)} kg payload; isolation begins above ${F(Math.SQRT2 * fn, 4)} Hz.` : 'Not defined for negative stiffness; reported as 0.', stable),
      kpi('f_lin', 'Natural frequency of the linear mount', fLin, 'Hz', 'Reference value with the vertical spring alone.', null),
      kpi('T_q', 'Transmissibility at excitation frequency', Tq, '–', `Linear single-degree-of-freedom transmissibility with the physical damper unchanged (effective damping ratio ${F(zEff, 3)}); the linear mount gives ${F(Tl, 3)}.`, Tq < 1),
      kpi('a_star', 'Compressed length for zero stiffness', aStar, 'm', 'a* = L₀·2k_h/(k_v + 2k_h): shorter than this the mount is bistable.', null),
      kpi('tol', 'Tuning tolerance', tolA, 'm', 'Change of compressed length that alters the stiffness by 5 % of k_v, from dK/da = 2k_hL₀/a²; it indicates the adjustment resolution required.', null),
      kpi('x_range', 'Low-stiffness travel', xRange, 'm', 'Displacement amplitude within which the tangent stiffness stays below 10 % of k_v (cubic hardening K = K₀ + 3k₃x²); 0 if the net stiffness already exceeds that.', xRange > 0),
      kpi('E_pre', 'Pre-compression energy', Epre, 'J', `Elastic energy stored in the two oblique springs, k_h(L₀ − a)²; the tuning actuator must hold ${F(Fpre, 4)} N per spring.`, null),
      kpi('t_settle', 'Settling time', tSettle, 's', '4/(ζω_n) of the tuned mount: low stiffness makes the mount slow to recover from shocks, the price of low-frequency isolation.', null),
      kpi('sag', 'Static deflection under payload', sag, 'm', 'm·g/k_v carried by the vertical spring; the static load capacity is not reduced by the negative-stiffness elements.', null),
    ],
    series: [
      { key: 'F_x', title: 'Restoring force and tangent stiffness', xlabel: 'displacement (mm)', ylabel: 'force (N) and stiffness (N/mm)', lines: [line('force', xs.map(v => v * 1e3), xs.map(Fx)), line('stiffness (N/mm)', xs.map(v => v * 1e3), xs.map(v => Kx(v) / 1e3)), line('linear spring force', xs.map(v => v * 1e3), xs.map(v => kv * v))], explain: 'The force curve is flat around the equilibrium (low stiffness) and hardens at larger displacement. A negative slope at the centre would indicate over-tuning.' },
      { key: 'T_f', title: 'Transmissibility', xlabel: 'frequency (Hz)', ylabel: 'transmissibility', logx: true, logy: true, lines: [line('tuned mount', fsw, fsw.map(f => stable ? Tr(f, fn, zEff) : 1)), line('linear mount', fsw, fsw.map(f => Tr(f, fLin, p.zeta_q)))], explain: 'Small-amplitude transmissibility. The tuned mount shifts the resonance down and begins to isolate at a much lower frequency than the linear mount of equal load capacity.' },
      { key: 'K_a', title: 'Net stiffness versus tuning length', xlabel: 'compressed length (mm)', ylabel: 'K / k_v', lines: [line('K(a)/k_v', as.map(v => v * 1e3), as.map(v => (kv - 2 * kh * (L0 / v - 1)) / kv)), line('zero', [as[0] * 1e3, as[as.length - 1] * 1e3], [0, 0])], explain: 'Tuning curve: the stiffness passes through zero at a*. Operating slightly to the right of the crossing gives quasi-zero, positive stiffness.' },
    ],
    tables: [mechTable([['Physical mechanism', 'Pre-compressed oblique springs act as a negative-stiffness element in parallel with a positive load-bearing spring'], ['Energy demand', `${F(Epre, 4)} J stored once in pre-compression; none to hold if the adjustment is self-locking`], ['Response time', `re-tuning is quasi-static (actuator speed); the mount itself settles in ${F(tSettle, 4)} s`], ['Loss', `unchanged physical damper; effective damping ratio rises from ${F(p.zeta_q, 3)} to ${F(zEff, 3)} as stiffness falls`], ['Reversibility', 'Fully reversible and continuous by changing the compressed length'], ['Hysteresis', 'None ideally; friction in pivots and guides adds a dead band that limits small-amplitude isolation'], ['Cycle durability', `Springs remain at high static compression (${F(100 * (1 - a / L0), 3)} % of free length): check relaxation and buckling of the oblique springs`]])],
    equations: [
      eqn('constitutive', 'Restoring force', 'F(x) = k_v x + 2k_h x[1 − L₀/√(x² + a²)]'),
      eqn('formula', 'Equilibrium stiffness and QZS condition', 'K₀ = k_v − 2k_h(L₀/a − 1);   K₀ = 0 at a* = 2k_hL₀/(k_v + 2k_h)'),
      eqn('formula', 'Cubic hardening', 'K(x) ≈ K₀ + 3k₃x²,  k₃ = k_hL₀/a³'),
      eqn('formula', 'Transmissibility', 'T = √{[1 + (2ζr)²]/[(1 − r²)² + (2ζr)²]},  r = f/f_n'),
    ],
    checks: [
      chk('Net stiffness positive (stable equilibrium)', stable, `K₀ = ${F(K0, 4)} N/m`),
      chk('Stability margin above 2 % of k_v', K0 >= 0.02 * kv, `K₀/k_v = ${F(K0 / kv, 4)}`),
      chk('Isolation at the excitation frequency', stable && Tq < 1, `T = ${F(Tq, 3)} (linear mount ${F(Tl, 3)})`),
      chk('Improves on the linear mount', stable && Tq <= Tl * (1 + 1e-9), `ratio ${F(Tq / Tl, 3)}`),
      chk('Oblique springs pre-compressed', p.a_q < p.L0_q, `a/L₀ = ${F(a / L0, 3)}`),
    ],
    warnings,
    recommendations: [
      rec('high', !stable ? `Lengthen the compressed length to at least ${F(aStar * 1.02e3, 4)} mm to restore positive stiffness with a 2 % margin.` : K0 > 0.2 * kv ? `Shorten the compressed length towards ${F(aStar * 1.03e3, 4)} mm: the present ${F(100 * K0 / kv, 3)} % residual stiffness leaves most of the benefit unused.` : `Tuning is effective (${F(100 * K0 / kv, 3)} % residual stiffness); provide an adjustment with resolution better than ${F(tolA * 1e6 / 5, 3)} µm to hold it.`, 'Stiffness-cancellation condition.'),
      rec('medium', `Keep vibration amplitude below ${F(xRange * 1e3, 3)} mm: outside this range cubic hardening raises the stiffness and resonance frequency (jump phenomena at large amplitude).`, 'Cubic stiffness of the geometry.'),
      rec('medium', 'Re-tune when the payload mass changes, because the low-stiffness point must coincide with the static equilibrium; an automatic levelling actuator makes this practical.', 'Static equilibrium must sit at x = 0.'),
      rec('low', 'A passive quasi-zero-stiffness mount isolates low frequencies without continuous power, where an active mount of equal performance would draw power permanently.', 'Passive stiffness cancellation.', true),
    ],
    publish: { 'mech.stiffness': K0, 'res.f0': fn, 'res.transmission': Tq, 'tune.range': 1 - Math.max(K0, 0) / kv, 'tune.energy_J': Epre, 'ctrl.settling_s': tSettle },
    numerics: { method: 'Closed-form geometric-nonlinear spring model and linearised transmissibility', dof: 1, steps: 200, iterations: 1, residual: 0, notes: 'Linear springs, frictionless pivots, small-amplitude (linearised) dynamics about the equilibrium.' },
  };
}

export default {
  id: 20,
  tagline: 'Designs feedback control for synthesis and resonators and evaluates physical tuning mechanisms for programmable materials: energy, speed, loss, hysteresis and durability.',
  inputs,
  presets: [
    { name: 'ZnO hydrothermal reactor, 1 L at 90 °C (SIMC PI)', note: 'Reference synthesis: closed-loop temperature control with anti-windup and a 20 % load disturbance.', values: { mode: 'pid', tuning: 'simc', V_reac: 1e-3, P_heater: 400, UA: 2, dead: 45, T_set: 90 } },
    { name: 'Same reactor, Ziegler–Nichols without anti-windup', note: 'Shows the overshoot caused by aggressive tuning and integrator wind-up during power-limited heat-up.', values: { mode: 'pid', tuning: 'zn', aw: 'off', V_reac: 1e-3, P_heater: 400, UA: 2, dead: 45, T_set: 90 } },
    { name: 'ZnO nanorod length control (MPC, temperature)', note: 'Tracks a 2 µm length target despite a 15 % slower real growth rate and precursor depletion.', values: { mode: 'mpc', manip: 'temp', L_target: 2e-6, t_growth: 14400, bias: -15, deplete: 0.2 } },
    { name: 'Active damping of a MEMS resonator (LQG)', note: 'State-feedback damping of a 10 kHz, Q = 200 resonator with a Kalman filter on displacement.', values: { mode: 'lqg', f0_r: 1e4, Q_r: 200, q_pos: 10, q_vel: 1, r_u: 1 } },
    { name: 'Varactor-tuned microwave metasurface cell', note: 'Split-ring resonator loaded with a 2 pF abrupt-junction varactor.', values: { mode: 'varactor', L_srr: 5, C_gap: 0.2, Cj0: 2, V_bias: 5, V_max: 20 } },
    { name: 'VO₂ terahertz switch metasurface', note: 'Volatile insulator–metal transition shorting a 1 THz resonator gap.', values: { mode: 'pcm', pcm_mat: 'vo2', f0_ms: 1e12, P_mw: 3, G_th: 2e-5 } },
    { name: 'GST non-volatile metasurface pixel', note: 'Terahertz Ge₂Sb₂Te₅ pixel on a micro-heater: set by annealing, reset by melt-quench, no holding power.', values: { mode: 'pcm', pcm_mat: 'gst', area_um2: 25, thick_nm: 100, G_th: 1e-5, Cth_mult: 3, P_mw: 8, f0_ms: 1e12, X0: 300, R_ms: 8 } },
    { name: 'Programmable bistable lattice, 8 units', note: 'Mechanical memory and reusable energy absorption from snap-through units.', values: { mode: 'bistable', n_units: 8, k_bar: 2e4, w_b: 10, h_b: 1.2, var_pct: 5 } },
    { name: 'Quasi-zero-stiffness isolator', note: 'Negative-stiffness springs tuned to cancel 92 % of the mount stiffness.', values: { mode: 'qzs', kq_v: 2000, kq_h: 1500, L0_q: 100, a_q: 62, m_q: 5, f_exc: 2 } },
    { name: 'YIG ferrite, negative-permeability band', note: 'Field-tunable negative effective permeability between resonance and anti-resonance.', values: { mode: 'ferrite', Ms: 1.4e5, H_b: 1e5, f_op: 6e9, alpha_g: 0.003 } },
  ],
  run(inp) {
    const p = read(inp), fnMap = { pid, mpc, lqg, varactor, pcm, piezo, lc: lcTune, mems, ferrite, bistable, qzs };
    return finish((fnMap[p.mode] || pid)(p));
  },
};
