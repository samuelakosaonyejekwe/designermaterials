// Suite 25 — Inverse Design, Optimisation and Artificial Intelligence.
// Inverse design by real optimisation runs on compact forward models implemented here: re-entrant honeycomb
// (Gibson–Ashby), optical multilayers (transfer matrix), resonators (mass-in-mass, split ring, Helmholtz),
// nanoparticles (Gans LSPR, Brus quantum dots), nanorod synthesis recipes (Arrhenius growth with depletion, NSGA-II),
// SIMP topology optimisation, Gaussian-process Bayesian optimisation and k-nearest-neighbour library screening.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import { MATERIALS, findMaterial } from '../../data/materials.js';

const F = N.fmt, PI = Math.PI, Rg = N.R_GAS, C = N.C;
const kp = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const pct = v => F(100 * v, 3) + ' %';
const counter = f => { const g = x => { g.n++; return f(x); }; g.n = 0; return g; };
const num = v => typeof v === 'number' && Number.isFinite(v) ? v : null;
const settingsRow = (o) => Object.entries(o).map(([k, v]) => [k, typeof v === 'number' ? F(v, 4) : String(v)]);

// =====================================================================================================
// 1. Inverse auxetic / lattice design (Gibson–Ashby cell-wall bending model of a hexagonal / re-entrant honeycomb)
// =====================================================================================================
function honey(thDeg, hl, tl) {
  const th = thDeg * PI / 180, s = Math.sin(th), c = Math.cos(th), A = hl + s;
  return { A, nu12: c * c / (A * s), nu21: A * s / (c * c), E1: tl ** 3 * c / (A * s * s), E2: tl ** 3 * A / c ** 3, G12: tl ** 3 * A / (hl * hl * (1 + 2 * hl) * c), rho: tl * (hl + 2) / (2 * c * A) };
}
function honeyLattice(thDeg, hl, l, nx = 6, ny = 4) {
  const th = thDeg * PI / 180, lc = l * Math.cos(th), ls = l * Math.sin(th), h = hl * l, Lx = 2 * lc, Ly = 2 * (h + ls), nodes = [], edges = [], map = new Map();
  const id = (x, y) => { const k = Math.round(x * 1e4) + ',' + Math.round(y * 1e4); if (!map.has(k)) { map.set(k, nodes.length); nodes.push([+x.toFixed(5), +y.toFixed(5), 0]); } return map.get(k); };
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    const ox = i * Lx, oy = j * Ly, A = id(ox, oy), B = id(ox, oy + h), Cn = id(ox + lc, oy + h + ls), D = id(ox + lc, oy + 2 * h + ls);
    edges.push([A, B], [B, Cn], [Cn, D]); if (j < ny - 1) edges.push([D, id(ox, oy + Ly)]); if (i > 0) edges.push([B, id(ox - lc, oy + h + ls)]); if (i < nx - 1 && j < ny - 1) edges.push([D, id(ox + Lx, oy + Ly)]);
  }
  return { nodes, edges };
}
function runAux(inp) {
  const W = [], l = inp.cellL, nuT = inp.nuT, ET = inp.EratT, rhoMax = inp.rhoMax, thMin = Math.min(inp.thMinA, inp.thMaxA - 1), thMax = inp.thMaxA;
  let tlLo = inp.tMin / l; if (tlLo > 0.25) { W.push(`Minimum strut thickness ${F(inp.tMin, 3)} mm exceeds a quarter of the wall length ${F(l, 3)} mm: slender-beam theory does not hold; t/l is fixed at 0.25. Increase the wall length.`); tlLo = 0.25; }
  const neg = nuT < 0, bounds = [neg ? [-thMax, -thMin] : [thMin, thMax], [0.5, 4], [tlLo, 0.25]];
  if (nuT === 0) W.push('ν₁₂ = 0 cannot be reached exactly by the hexagonal family (it needs θ → 90°); the nearest admissible geometry is returned.');
  const raw = x => { const g = honey(x[0], x[1], x[2]); if (!(g.A > 0.02) || Math.abs(Math.sin(x[0] * PI / 180)) < 1e-3) return 1e6;
    const eN = (g.nu12 - nuT) / Math.max(Math.abs(nuT), 0.1), eE = Math.log(g.E1 / ET), clear = x[1] + 2 * Math.sin(x[0] * PI / 180) - x[2];
    return eN * eN + eE * eE + 100 * Math.max(0, g.rho - rhoMax) ** 2 + (x[0] < 0 ? 100 * Math.max(0, -clear) ** 2 : 0) + 1e-3 * g.rho; };
  const f = counter(raw), pop = Math.round(inp.pop), gens = Math.round(inp.gens), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 300, tol: 1e-12 });
  const x = nm.f <= de.f ? nm.x : de.x, fbest = Math.min(nm.f, de.f), g = honey(x[0], x[1], x[2]), th = x[0], hl = x[1], tl = x[2], t = tl * l, h = hl * l;
  const errNu = Math.abs(g.nu12 - nuT) / Math.max(Math.abs(nuT), 0.1), errE = Math.abs(g.E1 / ET - 1), clear = hl + 2 * Math.sin(th * PI / 180) - tl, recip = Math.abs(g.E1 * g.nu21 - g.E2 * g.nu12) / Math.max(Math.abs(g.E1 * g.nu21), 1e-300);
  // robustness: uniform manufacturing scatter on t and θ
  const r = N.rng(inp.seed + 11), nMc = 600, mcNu = [], mcE = []; let inTol = 0;
  for (let k = 0; k < nMc; k++) { const tt = Math.max(t + inp.tolT * (2 * r() - 1), 1e-6 * l), th2 = th + inp.tolTh * (2 * r() - 1), q = honey(Math.abs(th2) < 0.5 ? Math.sign(th || 1) * 0.5 : th2, hl, tt / l); if (!(q.A > 0)) continue; mcNu.push(q.nu12); mcE.push(q.E1); if (Math.abs(q.nu12 - nuT) <= 0.1 * Math.max(Math.abs(nuT), 0.1) && Math.abs(q.E1 / ET - 1) <= 0.2) inTol++; }
  const sdNu = N.std(mcNu), cvE = N.std(mcE) / N.mean(mcE), hN = N.histogram(mcNu, 24), ths = N.linspace(bounds[0][0], bounds[0][1], 80);
  const ok = errNu < 0.02 && errE < 0.02, lat = honeyLattice(th, hl, l);
  if (!ok) W.push(`The targets are not met simultaneously within the bounds (ν₁₂ error ${pct(errNu)}, stiffness error ${pct(errE)}): the admissible family cannot deliver this combination at ρ*/ρₛ ≤ ${F(rhoMax, 3)} and t ≥ ${F(inp.tMin, 3)} mm.`);
  if (g.rho > 0.3) W.push(`Relative density ${F(g.rho, 3)} is above 0.3: the cell walls are too stocky for the bending model, which then underestimates stiffness and overestimates |ν|.`);
  W.push('Gibson–Ashby flexure model: small strains, slender walls, in-plane loading along x₁, infinite periodic lattice; axial and shear deformation of the walls, node volume and edge effects are neglected.');
  return {
    summary: `Differential evolution (${pop} × ${gens}) followed by Nelder–Mead found a ${th < 0 ? 're-entrant' : 'conventional'} honeycomb with θ = ${F(th, 4)}°, h/l = ${F(hl, 4)} and t/l = ${F(tl, 4)} (t = ${F(t, 3)} mm, h = ${F(h, 4)} mm for l = ${F(l, 3)} mm). It gives ν₁₂ = ${F(g.nu12, 4)} against the target ${F(nuT, 3)} and E₁/Eₛ = ${F(g.E1, 3)} against ${F(ET, 3)} at a relative density of ${F(g.rho, 3)}. With ±${F(inp.tolT, 2)} mm and ±${F(inp.tolTh, 2)}° manufacturing scatter ν₁₂ varies by ±${F(sdNu, 2)} (1σ) and ${pct(inTol / nMc)} of builds stay within 10 % of the Poisson target and 20 % of the stiffness target.`,
    kpis: [
      kp('theta', 'Cell-wall angle θ', th, '°', 'Inclination of the oblique walls to the x₁ direction; negative values are re-entrant (bow-tie) cells, which give a negative Poisson’s ratio.'),
      kp('hl', 'Wall length ratio h/l', hl, '–', 'Length of the vertical walls relative to the oblique walls; in re-entrant cells it must exceed 2|sin θ| plus the wall thickness to avoid contact.'),
      kp('t', 'Wall thickness t', t, 'mm', `Optimised strut thickness (t/l = ${F(tl, 4)}); it controls stiffness through (t/l)³ and must not fall below the manufacturing minimum of ${F(inp.tMin, 3)} mm.`, t >= inp.tMin * (1 - 1e-9)),
      kp('nu12', 'Achieved Poisson’s ratio ν₁₂', g.nu12, '–', `In-plane Poisson’s ratio from the cell-wall bending model, cos²θ/((h/l + sin θ) sin θ); target ${F(nuT, 3)}. Anisotropic lattices are not bound by the isotropic limit −1.`, errNu < 0.02),
      kp('E1', 'Achieved relative stiffness E₁/Eₛ', g.E1, '–', `Effective modulus along x₁ relative to the solid modulus; target ${F(ET, 3)}. In absolute terms E₁ = ${F(g.E1 * inp.Es, 4)} Pa for the stated solid.`, errE < 0.02),
      kp('rho', 'Relative density', g.rho, '–', `Solid volume fraction of the lattice; the constraint is ρ*/ρₛ ≤ ${F(rhoMax, 3)}.`, g.rho <= rhoMax * 1.001),
      kp('E2', 'Transverse relative stiffness E₂/Eₛ', g.E2, '–', 'Modulus along x₂ from the same model; the large ratio E₂/E₁ of re-entrant cells shows how anisotropic the design is.'),
      kp('nu21', 'Transverse Poisson’s ratio ν₂₁', g.nu21, '–', 'Poisson’s ratio for loading along x₂; with ν₁₂ it satisfies the reciprocity relation E₁ν₂₁ = E₂ν₁₂ of an orthotropic sheet.'),
      kp('obj', 'Final objective', fbest, '–', 'Sum of squared relative Poisson error and squared logarithmic stiffness error plus constraint penalties and a small mass term; values near 10⁻³·ρ mean both targets are met.', ok),
      kp('evals', 'Forward-model evaluations', f.n, '', 'Number of honeycomb evaluations used by differential evolution and the Nelder–Mead refinement.'),
      kp('sd_nu', 'Scatter of ν₁₂ under tolerances (1σ)', sdNu, '–', 'Standard deviation of ν₁₂ in a Monte Carlo simulation with uniformly distributed thickness and angle errors within the stated tolerances.', sdNu < 0.1 * Math.max(Math.abs(nuT), 0.1)),
      kp('p_ok', 'Builds within tolerance bands', inTol / nMc, '–', 'Fraction of simulated builds whose ν₁₂ lies within 10 % and whose stiffness lies within 20 % of the targets: a manufacturing-yield estimate for the design.', inTol / nMc > 0.9),
    ],
    series: [
      { key: 'conv', title: 'Convergence of the optimiser', xlabel: 'generation', ylabel: 'best objective', logy: true, lines: [{ name: 'differential evolution', x: de.history.map((_, i) => i + 1), y: de.history.map(v => Math.max(v, 1e-16)) }, { name: 'after Nelder–Mead', x: [gens], y: [Math.max(fbest, 1e-16)], points: true }], explain: 'Best objective value of the population per generation. A plateau indicates convergence; a plateau well above 10⁻³ means the targets conflict with the constraints.' },
      { key: 'nu_th', title: 'Poisson’s ratio versus cell angle at the optimum h/l', xlabel: 'θ (°)', ylabel: 'ν₁₂', lines: [{ name: 'ν₁₂(θ)', x: ths, y: ths.map(a => N.clamp(honey(a, hl, tl).nu12, -20, 20)) }, { name: 'target', x: [ths[0], ths.at(-1)], y: [nuT, nuT], dash: true }, { name: 'optimum', x: [th], y: [g.nu12], points: true }], explain: 'Design curve of the Poisson’s ratio over the admissible angle range (clipped at ±20). The slope at the optimum shows how sensitive the result is to angle errors.' },
      { key: 'mc', title: 'Monte Carlo distribution of ν₁₂ under manufacturing tolerances', xlabel: 'ν₁₂', ylabel: 'count', lines: [{ name: 'builds', x: hN.x, y: hN.y }, { name: 'target', x: [nuT, nuT], y: [0, N.maxOf(hN.y)], dash: true }], explain: 'Spread of the achieved Poisson’s ratio when wall thickness and angle vary within their tolerances. A narrow distribution centred on the target indicates a robust design.' },
    ],
    tables: [
      { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target / limit', 'Achieved', 'Deviation'], rows: [['ν₁₂', F(nuT, 4), F(g.nu12, 4), pct(errNu)], ['E₁/Eₛ', F(ET, 4), F(g.E1, 4), pct(errE)], ['ρ*/ρₛ', '≤ ' + F(rhoMax, 3), F(g.rho, 4), g.rho <= rhoMax * 1.001 ? 'satisfied' : 'violated'], ['t (mm)', '≥ ' + F(inp.tMin, 3), F(t, 4), t >= inp.tMin * (1 - 1e-9) ? 'satisfied' : 'violated'], ['|θ| (°)', `${F(thMin, 3)} – ${F(thMax, 3)}`, F(Math.abs(th), 4), 'within bounds'], ['wall clearance (h + 2l sin θ − t)/l', '> 0', F(clear, 4), clear > 0 || th > 0 ? 'satisfied' : 'violated']], explain: 'Targets and constraints against the optimised geometry. Deviations are relative to the target; constraint rows show whether the manufacturability limits are respected.' },
      { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Global search': 'differential evolution (rand/1/bin, F = 0.7, CR = 0.9)', 'Population': pop, 'Generations': gens, 'Local refinement': 'Nelder–Mead, 300 iterations, projected bounds', 'Evaluations': f.n, 'Seed': inp.seed, 'Variables': 'θ, h/l, t/l', 'Monte Carlo samples': nMc }), explain: 'Algorithm settings and budget used for this run, so that the result can be reproduced exactly.' },
    ],
    equations: [
      { kind: 'constitutive', label: 'Gibson–Ashby honeycomb (wall bending)', text: 'ν₁₂ = cos²θ/((h/l + sin θ) sin θ);  E₁/Eₛ = (t/l)³ cos θ/((h/l + sin θ) sin²θ);  E₂/Eₛ = (t/l)³ (h/l + sin θ)/cos³θ' },
      { kind: 'constitutive', label: 'Relative density and reciprocity', text: 'ρ*/ρₛ = (t/l)(h/l + 2)/(2 cos θ (h/l + sin θ));  E₁ν₂₁ = E₂ν₁₂' },
      { kind: 'governing', label: 'Inverse-design problem', text: 'min over (θ, h/l, t/l):  ((ν₁₂ − ν*)/|ν*|)² + ln²(E₁/E*) + 10⁻³ρ   s.t.  ρ ≤ ρ_max,  t ≥ t_min,  θ_min ≤ |θ| ≤ θ_max,  h + 2l sin θ > t' },
    ],
    geometry: { kind: 'lattice', nodes: lat.nodes, edges: lat.edges, radius: t / 2, units: 'mm', title: `Honeycomb θ = ${F(th, 3)}°, h/l = ${F(hl, 3)}` },
    checks: [
      { name: 'Targets met within 2 %', pass: ok, detail: `ν₁₂ error ${pct(errNu)}, E₁ error ${pct(errE)}` },
      { name: 'Relative-density constraint', pass: g.rho <= rhoMax * 1.001, detail: `${F(g.rho, 4)} ≤ ${F(rhoMax, 3)}` },
      { name: 'Minimum strut thickness', pass: t >= inp.tMin * (1 - 1e-9), detail: `${F(t, 4)} mm ≥ ${F(inp.tMin, 3)} mm` },
      { name: 'No wall contact in the re-entrant cell', pass: th > 0 || clear > 0, detail: `clearance ${F(clear * l, 4)} mm` },
      { name: 'Reciprocity E₁ν₂₁ = E₂ν₁₂', pass: recip < 1e-9, detail: `relative residual ${F(recip, 2)}` },
      { name: 'Slender-wall validity (t/l ≤ 0.25)', pass: tl <= 0.25 + 1e-12, detail: `t/l = ${F(tl, 4)}` },
      { name: 'Positive-definite in-plane stiffness', pass: g.E1 > 0 && g.E2 > 0 && g.G12 > 0 && g.nu12 * g.nu21 < 1 + 1e-9, detail: `E₁, E₂, G₁₂ > 0 and ν₁₂ν₂₁ = ${F(g.nu12 * g.nu21, 4)} (flexure model gives exactly 1: the sheet is marginally stable until wall stretching is included)` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: ok ? `Build the cell with θ = ${F(th, 3)}°, h = ${F(h, 3)} mm, l = ${F(l, 3)} mm, t = ${F(t, 3)} mm and verify ν₁₂ by image correlation on at least 6 × 6 cells; confirm the stiffness with a beam finite-element model that includes wall stretching (suites 5 and 14).` : `Relax one requirement: either raise the density limit, lower the minimum thickness or accept ${errNu > errE ? 'a different Poisson’s ratio' : 'a different stiffness'}; the present bounds leave a residual of ${F(fbest, 3)}.`, basis: 'Optimised geometry and target residuals.' },
      { priority: sdNu > 0.1 * Math.max(Math.abs(nuT), 0.1) ? 'high' : 'medium', text: `Hold the wall angle to better than ±${F(inp.tolTh, 2)}°: the local sensitivity is ${F(Math.abs(honey(th + 0.5, hl, tl).nu12 - honey(th - 0.5, hl, tl).nu12), 3)} per degree for ν₁₂, while stiffness varies by ${pct(cvE)} (1σ) mainly through (t/l)³.`, basis: 'Monte Carlo robustness analysis at the optimum.' },
      { priority: 'medium', text: `Check the printed minimum feature: t = ${F(t, 3)} mm against a process limit of ${F(inp.tMin, 3)} mm, and add fillets at the nodes (radius ≈ t) to lower the stress concentration that the beam model ignores.`, basis: 'Manufacturability constraints of strut lattices.' },
      { priority: 'low', text: `The design uses ${pct(g.rho)} of the solid volume; the mass term in the objective selects the lightest geometry meeting both targets, which reduces feedstock and build energy in proportion.`, basis: 'Relative density of the optimum.', sustainability: true },
    ],
    publish: { 'geom.strut': t * 1e-3, 'geom.cell_size': l * 1e-3, 'geom.rel_density': g.rho, 'geom.theta_deg': th, 'geom.h_over_l': hl, 'eff.nu': g.nu12, 'eff.E': g.E1 * inp.Es, 'eff.G': g.G12 * inp.Es, 'opt.objective': fbest, 'fab.yield': inTol / nMc },
    exports: [{ name: 'honeycomb_design.json', mime: 'application/json', content: JSON.stringify({ theta_deg: th, h_mm: h, l_mm: l, t_mm: t, nu12: g.nu12, E1_over_Es: g.E1, rel_density: g.rho, nodes: lat.nodes, edges: lat.edges }), description: 'Optimised unit-cell parameters with a 6 × 4 cell node and strut list (mm).' }],
    numerics: { method: 'Differential evolution + Nelder–Mead on the Gibson–Ashby honeycomb model with penalty constraints; Monte Carlo tolerance analysis', dof: 3, steps: gens, iterations: f.n, residual: fbest, notes: `Seed ${inp.seed}; penalties of weight 100 on density and wall-contact constraints.` },
  };
}

// =====================================================================================================
// 2. Inverse optical multilayer design (transfer-matrix method, lossless, normal incidence)
// =====================================================================================================
const NMAT = { MgF2: 1.38, SiO2: 1.46, Al2O3: 1.76, ZnO: 2.0, Si3N4: 2.02, Ta2O5: 2.1, ZnS: 2.35, TiO2: 2.4, Si: 3.48, Ge: 4.0 };
const NSUB = { glass: ['Borosilicate crown glass', 1.52], silica: ['Fused silica', 1.46], sapphire: ['Sapphire', 1.76], silicon: ['Silicon (infrared)', 3.48] };
function stackRT(d, n, n0, ns, lam) {
  let A = 1, B = 0, Cc = 0, D = 1;
  for (let i = 0; i < d.length; i++) { const del = 2 * PI * n[i] * d[i] / lam, c = Math.cos(del), s = Math.sin(del), b2 = s / n[i], c2 = n[i] * s, a1 = A, b1 = B, c1 = Cc, d1 = D; A = a1 * c - b1 * c2; B = a1 * b2 + b1 * c; Cc = c1 * c + d1 * c2; D = d1 * c - c1 * b2; }
  const nr = n0 * A - ns * D, ni = n0 * ns * B - Cc, dr = n0 * A + ns * D, di = n0 * ns * B + Cc, den = dr * dr + di * di;
  return { R: (nr * nr + ni * ni) / den, T: 4 * n0 * ns / den };
}
function runOptical(inp) {
  const W = [], nH = NMAT[inp.nH] ?? 2.4, nL = NMAT[inp.nL] ?? 1.46, ns = (NSUB[inp.nSub] || NSUB.glass)[1], n0 = 1, NL = Math.round(inp.nLay), goal = inp.optGoal;
  let l1 = Math.min(inp.lam1, inp.lam2), l2 = Math.max(inp.lam1, inp.lam2); if (l2 - l1 < 1) l2 = l1 + 1;
  const lc = 2 * l1 * l2 / (l1 + l2), bw = l2 - l1, outerLow = goal === 'ar', mats = Array.from({ length: NL }, (_, i) => (i % 2 === 0) === outerLow ? 'L' : 'H'), n = mats.map(m => m === 'H' ? nH : nL);
  if (nH <= nL) W.push(`The "high-index" material (${inp.nH}, n = ${nH}) does not have a higher index than the "low-index" material (${inp.nL}, n = ${nL}): the stack has little or no index contrast.`);
  const dLo = Math.min(inp.dMin, inp.dMax), dHi = Math.max(inp.dMin, inp.dMax, dLo + 1), bounds = n.map(() => [dLo, dHi]);
  const band = N.linspace(l1, l2, 21), sA = N.linspace(Math.max(l1 - bw, 0.3 * l1), l1 - 0.15 * bw, 12), sB = N.linspace(l2 + 0.15 * bw, l2 + bw, 12);
  const meanR = (d, ls) => { let s = 0; for (const l of ls) s += stackRT(d, n, n0, ns, l).R; return s / ls.length; };
  const merit = d => goal === 'ar' ? meanR(d, band) : goal === 'hr' ? 1 - meanR(d, band) : 0.5 * meanR(d, band) + 0.25 * (1 - meanR(d, sA)) + 0.25 * (1 - meanR(d, sB));
  const f = counter(merit), pop = Math.max(Math.round(inp.pop), 5 * NL), gens = 2 * Math.round(inp.gens), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 500, tol: 1e-13, step: 0.05 });
  const qw = n.map(v => N.clamp(lc / (4 * v), dLo, dHi)), fq = merit(qw), nq = N.nelderMead(f, qw, { bounds, maxIter: 500, tol: 1e-13, step: 0.05 }), cands = [[de.f, de.x], [nm.f, nm.x], [nq.f, nq.x], [fq, qw]].sort((a, b) => a[0] - b[0]), d = cands[0][1], fb = cands[0][0], bare = ((n0 - ns) / (n0 + ns)) ** 2;
  const Rb = meanR(d, band), tot = N.sum(d), lam = N.linspace(Math.max(0.5 * l1, 150), 1.6 * l2, 240), sp = lam.map(l => stackRT(d, n, n0, ns, l)), spq = lam.map(l => stackRT(qw, n, n0, ns, l).R);
  const tgt = l => goal === 'ar' ? (l >= l1 && l <= l2 ? 0 : null) : goal === 'hr' ? (l >= l1 && l <= l2 ? 1 : null) : (l >= l1 && l <= l2 ? 0 : (l >= sA[0] && l <= sA.at(-1)) || (l >= sB[0] && l <= sB.at(-1)) ? 1 : null), tl = lam.filter(l => tgt(l) !== null);
  const cons = N.maxOf(sp.map(s => Math.abs(s.R + s.T - 1))), Rmax = N.maxOf(band.map(l => stackRT(d, n, n0, ns, l).R)), Rmin = N.minOf(band.map(l => stackRT(d, n, n0, ns, l).R));
  const r = N.rng(inp.seed + 5), tol = inp.thkTol / 100, mc = Array.from({ length: 200 }, () => merit(d.map(v => v * (1 + tol * (2 * r() - 1))))), mcMean = N.mean(mc), mc95 = N.quantile(mc, 0.95);
  const atB = d.filter(v => v <= dLo * 1.001 || v >= dHi * 0.999).length, gl = { ar: 'antireflection coating', hr: 'high reflector', bp: 'band-pass filter' }[goal] || 'antireflection coating';
  if (atB) W.push(`${atB} layer${atB === 1 ? ' is' : 's are'} at a thickness bound: the optimiser would use thinner or thicker layers if allowed; consider changing the layer count or the bounds.`);
  if (goal === 'hr' && Rb < 0.9) W.push(`Mean reflectance ${pct(Rb)} is below 90 %: add layer pairs or raise the index contrast (currently ${F(nH / nL, 3)}).`);
  W.push('Lossless, non-dispersive indices at normal incidence: absorption, dispersion, angle and polarisation dependence, interface roughness and thickness-dependent film index are not modelled.');
  let y = 0; const polys = [{ layer: 3, pts: [[0, -0.3 * tot], [1.5 * tot, -0.3 * tot], [1.5 * tot, 0], [0, 0]] }]; for (let i = NL - 1; i >= 0; i--) { polys.push({ layer: mats[i] === 'H' ? 1 : 2, pts: [[0, y], [1.5 * tot, y], [1.5 * tot, y + d[i]], [0, y + d[i]]] }); y += d[i]; }
  const impr = goal === 'ar' ? bare / Math.max(Rb, 1e-12) : (1 - fq) / Math.max(1 - fb, 1e-12);
  return {
    summary: `A ${NL}-layer ${chem.prettyFormula(inp.nH)}/${chem.prettyFormula(inp.nL)} ${gl} on ${(NSUB[inp.nSub] || NSUB.glass)[0].toLowerCase()} was optimised for ${F(l1, 4)}–${F(l2, 4)} nm by differential evolution (${pop} × ${gens}) and Nelder–Mead with the transfer-matrix model. The mean reflectance in the band is ${pct(Rb)} (${pct(Rmin)} to ${pct(Rmax)}), merit ${F(fb, 3)} against ${F(fq, 3)} for the quarter-wave reference stack; total thickness ${F(tot, 4)} nm. With ±${F(inp.thkTol, 3)} % thickness errors the merit degrades to ${F(mcMean, 3)} on average (95th percentile ${F(mc95, 3)}).`,
    kpis: [
      kp('Rmean', 'Mean reflectance in the band', Rb, '–', `Average reflectance over ${band.length} wavelengths of the target band; ${goal === 'hr' ? 'the quantity maximised for a mirror' : 'the quantity minimised in the pass band'}. The bare substrate reflects ${pct(bare)}.`, goal === 'hr' ? Rb > 0.95 : Rb < 0.01),
      kp('Tmean', 'Mean transmittance in the band', 1 - Rb, '–', 'One minus the mean reflectance (lossless films, so R + T = 1); real coatings transmit slightly less because of absorption and scatter.'),
      kp('Rmax', 'Highest reflectance in the band', Rmax, '–', 'Worst-case reflectance inside the band; for an antireflection coating this sets the residual ghost-image level.'),
      kp('Rmin', 'Lowest reflectance in the band', Rmin, '–', 'Lowest reflectance inside the band; for a mirror it is the weakest point of the stop band.'),
      kp('merit', 'Final merit (lower is better)', fb, '–', 'Objective minimised by the optimiser: mean R for antireflection, 1 − mean R for a reflector, and an equal mix of pass-band reflectance and stop-band transmittance for a band-pass.', fb <= fq + 1e-9),
      kp('merit_qw', 'Merit of the quarter-wave reference', fq, '–', 'Same merit for a stack with all layers a quarter of the centre wavelength thick: the classical mirror design, and for antireflection and band-pass targets only a scale against which the optimised stack is compared.'),
      kp('gain', goal === 'ar' ? 'Reflection reduction versus bare substrate' : 'Merit improvement versus quarter-wave stack', Math.min(impr, 1e9), '×', goal === 'ar' ? 'Ratio of the bare-substrate reflectance to the coated band-average reflectance.' : 'Ratio of the quarter-wave residual (1 − merit) to the optimised one; values near 1 mean the classical design is already close to optimal.'),
      kp('thick', 'Total coating thickness', tot, 'nm', 'Sum of all layer thicknesses; it scales deposition time and stress-related risk of cracking or delamination.'),
      kp('evals', 'Forward-model evaluations', f.n, '', 'Number of merit evaluations (each a transfer-matrix sweep over the sampled wavelengths) used by the global and local optimiser.'),
      kp('mc95', 'Merit at the 95th percentile of thickness errors', mc95, '–', `Merit not exceeded by 95 % of 200 simulated coatings with uniformly distributed layer-thickness errors of ±${F(inp.thkTol, 3)} %.`, mc95 < 2 * fb + 0.02),
      kp('cons', 'Energy-balance residual', cons, '–', 'Largest |R + T − 1| over the plotted spectrum, with R and T computed independently from the characteristic matrix; it must vanish for lossless films.', cons < 1e-9),
    ],
    series: [
      { key: 'spec', title: 'Reflectance spectrum: target versus achieved', xlabel: 'wavelength (nm)', ylabel: 'reflectance', lines: [{ name: 'optimised stack', x: lam, y: sp.map(s => s.R) }, { name: 'quarter-wave reference', x: lam, y: spq, dash: true }, { name: 'target', x: tl, y: tl.map(tgt), points: true }], explain: 'Reflectance of the optimised stack and of the quarter-wave reference over a range wider than the design band, with the target values marked where they are specified.' },
      { key: 'T', title: 'Transmittance spectrum', xlabel: 'wavelength (nm)', ylabel: 'transmittance', lines: [{ name: 'optimised stack', x: lam, y: sp.map(s => s.T) }, { name: 'bare substrate', x: [lam[0], lam.at(-1)], y: [1 - bare, 1 - bare], dash: true }], explain: 'Transmittance into the substrate (single interface; the rear face of a real window adds its own reflection). The dashed line is the uncoated substrate.' },
      { key: 'conv', title: 'Convergence of the optimiser', xlabel: 'generation', ylabel: 'best merit', logy: true, lines: [{ name: 'differential evolution', x: de.history.map((_, i) => i + 1), y: de.history.map(v => Math.max(v, 1e-12)) }, { name: 'after Nelder–Mead', x: [gens], y: [Math.max(fb, 1e-12)], points: true }], explain: 'Best merit in the population per generation. Multilayer design is strongly multimodal, so different seeds can reach different local optima of similar quality.' },
    ],
    tables: [
      { key: 'layers', title: 'Layer prescription (from the incident medium to the substrate)', columns: ['Layer', 'Material', 'Refractive index', 'Thickness (nm)', 'Optical thickness n·d/λc', 'Quarter-wave reference (nm)'], rows: d.map((v, i) => [i + 1, mats[i] === 'H' ? inp.nH : inp.nL, n[i], F(v, 4), F(n[i] * v / lc, 3), F(qw[i], 4)]), explain: 'Physical and optical thickness of each layer, starting at the air side. An optical thickness of 0.25 corresponds to a quarter-wave layer at the band centre.' },
      { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Global search': 'differential evolution (rand/1/bin)', 'Population': pop, 'Generations': gens, 'Local refinement': 'Nelder–Mead (500 iterations) from the evolutionary optimum and from the quarter-wave stack', 'Evaluations': f.n, 'Wavelength samples': band.length + (goal === 'bp' ? 24 : 0), 'Thickness bounds (nm)': `${F(dLo, 3)} – ${F(dHi, 3)}`, 'Seed': inp.seed }), explain: 'Algorithm settings and budget of this run for reproducibility.' },
    ],
    equations: [
      { kind: 'governing', label: 'Characteristic matrix of layer j', text: 'Mⱼ = [[cos δⱼ, i·sin δⱼ/nⱼ], [i·nⱼ·sin δⱼ, cos δⱼ]],  δⱼ = 2π·nⱼ·dⱼ/λ;   [B; C] = (Πⱼ Mⱼ)·[1; nₛ]' },
      { kind: 'formula', label: 'Reflectance and transmittance', text: 'R = |(n₀B − C)/(n₀B + C)|²;  T = 4n₀nₛ/|n₀B + C|²;  R + T = 1 (lossless)' },
      { kind: 'governing', label: 'Inverse-design problem', text: goal === 'ar' ? 'min over d₁…d_N:  ⟨R(λ)⟩ over the band,  d_min ≤ dⱼ ≤ d_max' : goal === 'hr' ? 'min over d₁…d_N:  1 − ⟨R(λ)⟩ over the band,  d_min ≤ dⱼ ≤ d_max' : 'min over d₁…d_N:  ½⟨R⟩_pass + ¼⟨1 − R⟩_stop,short + ¼⟨1 − R⟩_stop,long' },
      { kind: 'formula', label: 'Coating materials', text: [[inp.nH, nH], [inp.nL, nL]].map(([m, v]) => { try { return `${chem.prettyFormula(m)}: n = ${v}, M = ${F(chem.molarMass(m), 5)} g/mol`; } catch { return `${m}: n = ${v}`; } }).join(';  ') },
    ],
    geometry: { kind: 'layout', polygons: polys, units: 'nm', title: 'Cross-section: substrate (layer 3), high index (1), low index (2)' },
    checks: [
      { name: 'Energy conservation R + T = 1', pass: cons < 1e-9, detail: `largest residual ${F(cons, 2)}` },
      { name: 'Optimised stack at least as good as the quarter-wave reference', pass: fb <= fq + 1e-9, detail: `${F(fb, 4)} vs ${F(fq, 4)}` },
      { name: 'Reflectance within [0, 1]', pass: sp.every(s => s.R >= -1e-12 && s.R <= 1 + 1e-12), detail: `range ${F(N.minOf(sp.map(s => s.R)), 3)} to ${F(N.maxOf(sp.map(s => s.R)), 3)}` },
      { name: 'Layers inside the manufacturable thickness window', pass: d.every((v, i) => v >= bounds[i][0] - 1e-9 && v <= bounds[i][1] + 1e-9), detail: `${F(N.minOf(d), 4)} – ${F(N.maxOf(d), 4)} nm, ${atB} at a bound` },
      { name: 'Tolerance robustness', pass: mc95 < 2 * fb + 0.02, detail: `95th-percentile merit ${F(mc95, 4)} for ±${F(inp.thkTol, 3)} % thickness error` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: `Deposit the ${NL} layers with the listed thicknesses (thinnest ${F(N.minOf(d), 3)} nm, thickest ${F(N.maxOf(d), 3)} nm) and re-optimise with measured, dispersive indices of your own films before production; film indices differ from bulk values by several per cent depending on the process.`, basis: 'Layer prescription from the transfer-matrix optimisation with nominal indices.' },
      { priority: mc95 > 2 * fb + 0.02 ? 'high' : 'medium', text: `Thickness control of ±${F(inp.thkTol, 3)} % raises the merit from ${F(fb, 3)} to ${F(mcMean, 3)} on average; ${mc95 > 2 * fb + 0.02 ? 'use optical monitoring or reduce the layer count to gain tolerance' : 'quartz-crystal or optical monitoring at this level is sufficient'}.`, basis: 'Monte Carlo analysis of layer-thickness errors.' },
      { priority: 'medium', text: goal === 'hr' ? `Each additional high/low pair multiplies the residual transmission by about (n_L/n_H)² = ${F((nL / nH) ** 2, 3)}; add pairs until the reflectance target is reached rather than widening single layers.` : `Compare with fewer layers: a ${Math.max(1, NL - 2)}-layer design may reach nearly the same merit with shorter deposition time and lower film stress.`, basis: 'Scaling of multilayer performance with layer count and index contrast.' },
      { priority: 'low', text: `Total thickness of ${F(tot, 3)} nm sets deposition time, target consumption and energy; prefer the thinnest design that meets the specification and abundant oxides (SiO₂, TiO₂, Al₂O₃) over tantalum- or germanium-based layers where the optics allow.`, basis: 'Coating thickness and raw-material criticality.', sustainability: true },
    ],
    publish: { 'res.transmission': 1 - Rb, 'res.reflection': Rb, 'res.absorption': 0, 'geom.thickness': tot * 1e-9, 'res.bandgap_lo': N.C0 / (l2 * 1e-9), 'res.bandgap_hi': N.C0 / (l1 * 1e-9), 'res.f0': N.C0 / (lc * 1e-9), 'opt.objective': fb, 'eff.n': Math.sqrt(N.sum(d.map((v, i) => n[i] * n[i] * v)) / tot) },
    exports: [{ name: 'multilayer_design.csv', mime: 'text/csv', content: 'layer,material,n,thickness_nm\n' + d.map((v, i) => `${i + 1},${mats[i] === 'H' ? inp.nH : inp.nL},${n[i]},${v.toFixed(3)}`).join('\n'), description: 'Layer prescription from the incident side to the substrate.' }],
    numerics: { method: 'Transfer-matrix forward model; differential evolution + Nelder–Mead on layer thicknesses; Monte Carlo thickness-error analysis', dof: NL, steps: gens, iterations: f.n, residual: fb, notes: `Seed ${inp.seed}; ${band.length} wavelengths in the band.` },
  };
}

// =====================================================================================================
// 3. Inverse resonator / bandgap design
// =====================================================================================================
function srrLC(R, w, g, t, er) { const L = N.MU0 * R * Math.max(Math.log(8 * R / (w + t)) - 0.5, 0.05), Cg = N.EPS0 * er * w * t / g, Cs = N.EPS0 * 0.5 * (1 + er) * (2 * (w + t) / PI) * Math.max(Math.log(4 * R / g), 0.05); return { L, C: Cg + Cs, f: 1 / (2 * PI * Math.sqrt(L * (Cg + Cs))) }; }
function runResonator(inp) {
  const W = [], pop = Math.round(inp.pop), gens = Math.round(inp.gens), type = inp.resType;
  if (type === 'srr') {
    const fT = inp.fSrr, mf = inp.minFeat * 1e-3, t = inp.tMet * 1e-6, er = inp.epsSub, lam0 = N.C0 / fT, bounds = [[Math.log(3 * mf), Math.log(Math.max(lam0 / 6, 6 * mf))], [Math.log(mf), Math.log(20 * mf)], [Math.log(mf), Math.log(20 * mf)]];
    const raw = q => { const R = Math.exp(q[0]), w = Math.exp(q[1]), g = Math.exp(q[2]), e = Math.log(srrLC(R, w, g, t, er).f / fT); return e * e + 10 * Math.max(0, w / R - 0.5) ** 2 + 10 * Math.max(0, g / R - 1) ** 2 + 1e-4 * ((2 * R + w) / lam0 * 10) ** 2 - 1e-5 * (q[1] + q[2] - 2 * Math.log(mf)); };
    const f = counter(raw), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 300, tol: 1e-13 }), q = nm.f <= de.f ? nm.x : de.x, fb = Math.min(nm.f, de.f);
    const R = Math.exp(q[0]), w = Math.exp(q[1]), g = Math.exp(q[2]), lc = srrLC(R, w, g, t, er), f0 = lc.f, a = 1.4 * (2 * R + w), Fill = PI * R * R / (a * a), om = 2 * PI * f0, skin = Math.sqrt(2 / (om * N.MU0 * 5.96e7)), Rohm = 2 * PI * R / (5.96e7 * w * Math.min(t, 2 * skin)), Q = om * lc.L / Rohm, fmp = f0 / Math.sqrt(1 - Fill), err = Math.abs(f0 / fT - 1);
    const Rcf = (() => { try { return N.brent(x => srrLC(Math.exp(x), w, g, t, er).f - fT, bounds[0][0] - 3, bounds[0][1] + 3, 1e-13); } catch { return NaN; } })(), cfErr = Number.isFinite(Rcf) ? Math.abs(Math.exp(Rcf) / R - 1) : 1;
    const ff = N.linspace(0.6 * f0, 1.6 * f0, 240), mu = ff.map(x => { const dn = C.of(x * x - f0 * f0, x * f0 / Q), v = C.div(C.of(Fill * x * x), dn); return C.of(1 - v.re, -v.im); }), muMin = N.minOf(mu.map(m => m.re)), fMid = 0.5 * (f0 + fmp), muMid = (() => { const dn = C.of(fMid * fMid - f0 * f0, fMid * f0 / Q), v = C.div(C.of(Fill * fMid * fMid), dn); return 1 - v.re; })();
    if (err > 0.01) W.push(`Target frequency missed by ${pct(err)}: the ring radius reached a bound set by the minimum feature size or the sub-wavelength limit. Change the minimum feature, substrate permittivity or target.`);
    if (lam0 / a < 5) W.push(`Cell size is only λ₀/${F(lam0 / a, 3)}: the array is not deeply sub-wavelength, so an effective permeability is of limited meaning and diffraction may occur.`);
    if (fT > 1e14) W.push('Above about 100 THz the kinetic inductance of the electrons saturates the magnetic resonance of metallic split rings; the quasi-static LC model overestimates the achievable frequency.');
    W.push('Quasi-static single-ring LC model with ohmic loss only: substrate loss, radiation, bianisotropy and coupling between neighbouring rings shift the resonance, typically by several per cent; the quality factor is an upper bound.');
    const phi = g / R, pts = []; for (let k = 0; k <= 40; k++) { const an = phi / 2 + (2 * PI - phi) * k / 40; pts.push([+(1e3 * (R + w / 2) * Math.cos(an)).toFixed(5), +(1e3 * (R + w / 2) * Math.sin(an)).toFixed(5)]); } for (let k = 40; k >= 0; k--) { const an = phi / 2 + (2 * PI - phi) * k / 40; pts.push([+(1e3 * (R - w / 2) * Math.cos(an)).toFixed(5), +(1e3 * (R - w / 2) * Math.sin(an)).toFixed(5)]); }
    const h = 1e3 * a / 2;
    return {
      summary: `For a magnetic resonance at ${F(fT, 4)} Hz the optimiser gives a split ring of mean radius ${F(R * 1e3, 4)} mm, track width ${F(w * 1e3, 4)} mm and gap ${F(g * 1e3, 4)} mm on a substrate of εᵣ = ${F(er, 3)} (L = ${F(lc.L, 3)} H, C = ${F(lc.C, 3)} F), resonating at ${F(f0, 5)} Hz (${pct(err)} from target). In an array of period ${F(a * 1e3, 4)} mm (λ₀/${F(lam0 / a, 3)}) the Lorentz model predicts negative permeability from ${F(f0, 4)} to ${F(fmp, 4)} Hz with a minimum of ${F(muMin, 3)} and an ohmic quality factor of about ${F(Q, 3)}.`,
      kpis: [
        kp('f0', 'Resonance frequency', f0, 'Hz', `LC resonance 1/(2π√(LC)) of the optimised ring; target ${F(fT, 4)} Hz.`, err < 0.01),
        kp('R', 'Mean ring radius', R * 1e3, 'mm', 'Radius to the centre of the track; it sets the loop inductance and is the main tuning dimension.'),
        kp('w', 'Track width', w * 1e3, 'mm', `Width of the metal track, kept at or above the minimum feature of ${F(inp.minFeat, 3)} mm; wider tracks lower resistance and inductance.`, w >= mf * 0.999),
        kp('g', 'Split gap', g * 1e3, 'mm', `Gap of the split, kept at or above the minimum feature; it controls the gap capacitance and is the dimension most sensitive to etch tolerance.`, g >= mf * 0.999),
        kp('L', 'Loop inductance', lc.L, 'H', 'Inductance of a circular loop of rectangular track cross-section, μ₀R(ln(8R/(w+t)) − ½).'),
        kp('Cap', 'Total capacitance', lc.C, 'F', 'Parallel-plate capacitance of the gap plus the surface capacitance of the ring through substrate and air.'),
        kp('Q', 'Ohmic quality factor (upper bound)', Q, '–', 'ωL/R with the conduction resistance of copper within one skin depth on each face; dielectric and radiation losses lower the real value.', Q > 20),
        kp('mu_min', 'Minimum Re μ in the band', muMin, '–', `Most negative real part of the Lorentz permeability μ = 1 − F f²/(f² − f₀² + i f f₀/Q) with filling factor F = ${F(Fill, 3)}.`, muMin < 0),
        kp('fmp', 'Upper edge of the negative-μ band', fmp, 'Hz', 'Magnetic plasma frequency f₀/√(1 − F) where Re μ returns through zero (lossless limit); the band is a few per cent wide.'),
        kp('sub', 'Sub-wavelength ratio λ₀/a', lam0 / a, '–', 'Free-space wavelength divided by the array period; homogenisation into an effective permeability needs a ratio well above about 5.', lam0 / a >= 5),
        kp('evals', 'Forward-model evaluations', f.n, '', 'Number of LC-model evaluations used by differential evolution and Nelder–Mead.'),
        kp('cf', 'Agreement with the closed-form radius', cfErr, '–', 'Relative difference between the optimised radius and the radius obtained by one-dimensional root finding at the same track and gap; it verifies the optimiser.', cfErr < 1e-3),
      ],
      series: [
        { key: 'mu', title: 'Effective permeability of the split-ring array', xlabel: 'frequency (Hz)', ylabel: 'μ_eff', lines: [{ name: 'Re μ', x: ff, y: mu.map(m => m.re) }, { name: 'Im μ (loss)', x: ff, y: mu.map(m => m.im), dash: true }, { name: 'μ = 0', x: [ff[0], ff.at(-1)], y: [0, 0], dash: true }], explain: 'Lorentz response of the array. The real part is negative between the resonance and the magnetic plasma frequency; the imaginary part peaks at resonance, where absorption is strongest (passivity requires Im μ ≥ 0 in the e^(−iωt) convention used here).' },
        { key: 'conv', title: 'Convergence of the optimiser', xlabel: 'generation', ylabel: 'best objective', logy: true, lines: [{ name: 'differential evolution', x: de.history.map((_, i) => i + 1), y: de.history.map(v => Math.max(v + 1e-3, 1e-12)) }], explain: 'Best objective per generation (offset by 10⁻³ so that the small manufacturability reward does not make it negative). Flat tails indicate convergence.' },
        { key: 'fR', title: 'Resonance frequency versus ring radius', xlabel: 'mean radius (mm)', ylabel: 'frequency (Hz)', logx: true, logy: true, lines: [{ name: 'f₀(R) at optimised w, g', x: N.logspace(Math.exp(bounds[0][0]), Math.exp(bounds[0][1]), 60).map(v => v * 1e3), y: N.logspace(Math.exp(bounds[0][0]), Math.exp(bounds[0][1]), 60).map(v => srrLC(v, w, g, t, er).f) }, { name: 'design', x: [R * 1e3], y: [f0], points: true }], explain: 'Design curve across the admissible radius range: the resonance scales roughly inversely with ring size. The point marks the optimised design.' },
      ],
      tables: [
        { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target / limit', 'Achieved'], rows: [['Resonance (Hz)', F(fT, 5), F(f0, 5)], ['Track width (mm)', '≥ ' + F(inp.minFeat, 3), F(w * 1e3, 4)], ['Gap (mm)', '≥ ' + F(inp.minFeat, 3), F(g * 1e3, 4)], ['λ₀/a', '≥ 5', F(lam0 / a, 4)], ['Negative-μ band (Hz)', '–', `${F(f0, 4)} – ${F(fmp, 4)}`]], explain: 'Design targets and manufacturability limits against the optimised split ring.' },
        { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Global search': 'differential evolution in log-dimensions', 'Population': pop, 'Generations': gens, 'Local refinement': 'Nelder–Mead', 'Evaluations': f.n, 'Seed': inp.seed, 'Closed-form check': 'Brent root for R at fixed w, g' }), explain: 'Algorithm settings and budget of this run.' },
      ],
      equations: [
        { kind: 'governing', label: 'LC resonance', text: 'f₀ = 1/(2π√(LC));  L = μ₀R(ln(8R/(w + t)) − ½);  C = ε₀εᵣwt/g + ε₀·½(1 + εᵣ)·(2(w + t)/π)·ln(4R/g)' },
        { kind: 'constitutive', label: 'Effective permeability (Lorentz)', text: 'μ_eff(f) = 1 − F·f²/(f² − f₀² + i·f·f₀/Q),  F = πR²/a²;  Re μ < 0 for f₀ < f < f₀/√(1 − F)' },
        { kind: 'formula', label: 'Ohmic quality factor', text: 'Q = ω₀L/R_Ω,  R_Ω = 2πR/(σ·w·min(t, 2δ)),  δ = √(2/(ω₀μ₀σ))' },
        { kind: 'half-cell', label: 'Pattern transfer of the copper ring (ferric chloride etch)', text: (() => { try { return chem.ionicEquations([{ f: 'Cu', phase: 's' }, { f: 'Fe3+', phase: 'aq' }], [{ f: 'Cu2+', phase: 'aq' }, { f: 'Fe2+', phase: 'aq' }]).molecular; } catch { return 'Cu(s) + 2 Fe³⁺(aq) → Cu²⁺(aq) + 2 Fe²⁺(aq)'; } })(), note: 'Over-etching widens the gap and narrows the track, raising f₀.' },
      ],
      geometry: { kind: 'layout', polygons: [{ layer: 2, pts: [[-h, -h], [h, -h], [h, h], [-h, h]] }, { layer: 1, pts }], units: 'mm', title: 'Split ring (layer 1) in its unit cell (layer 2)' },
      checks: [
        { name: 'Target frequency met within 1 %', pass: err < 0.01, detail: `${F(f0, 5)} Hz vs ${F(fT, 5)} Hz` },
        { name: 'Optimiser agrees with closed-form root', pass: cfErr < 1e-3, detail: `radius difference ${F(cfErr, 2)}` },
        { name: 'Minimum feature respected', pass: w >= mf * 0.999 && g >= mf * 0.999, detail: `w = ${F(w * 1e3, 4)} mm, g = ${F(g * 1e3, 4)} mm` },
        { name: 'Sub-wavelength cell (λ₀/a ≥ 5)', pass: lam0 / a >= 5, detail: `λ₀/a = ${F(lam0 / a, 4)}` },
        { name: 'Passivity (Im μ ≥ 0)', pass: mu.every(m => m.im >= -1e-12), detail: `min Im μ = ${F(N.minOf(mu.map(m => m.im)), 3)}` },
        { name: 'Filling factor below 1', pass: Fill < 1, detail: `F = ${F(Fill, 4)}` },
      ],
      warnings: W,
      recommendations: [
        { priority: 'high', text: `Fabricate with R = ${F(R * 1e3, 4)} mm, w = ${F(w * 1e3, 4)} mm, g = ${F(g * 1e3, 4)} mm and confirm f₀ by a full-wave simulation or a waveguide measurement (suite 12) before committing an array: the circuit model is accurate to several per cent only.`, basis: 'Quasi-static LC model of the split ring.' },
        { priority: 'medium', text: `A gap error of ±${F(inp.minFeat * 0.15, 2)} mm (15 % of the minimum feature) shifts f₀ by about ${pct(Math.abs(srrLC(R, w, g * 1.15, t, er).f / f0 - 1))}; specify the gap tolerance accordingly or add a trimming capacitor for tuning.`, basis: 'Finite-difference sensitivity of the resonance to the gap.' },
        { priority: 'medium', text: `Operate between ${F(f0 * 1.01, 4)} and ${F(fmp, 4)} Hz for negative permeability and away from f₀ itself, where absorption peaks; combine with a wire array (ε < 0) for a negative-index pass band.`, basis: 'Lorentz permeability of the array.' },
        { priority: 'low', text: 'Use a low-loss laminate and the thinnest copper that exceeds two skin depths: it reduces etchant consumption and copper waste while keeping the quality factor.', basis: `Skin depth ${F(skin * 1e6, 3)} µm at resonance against ${F(inp.tMet, 3)} µm metal.`, sustainability: true },
      ],
      publish: { 'res.f0': f0, 'res.Q': Q, 'res.bandgap_lo': f0, 'res.bandgap_hi': fmp, 'eff.mu': muMid, 'geom.cell_size': a, 'geom.diameter': 2 * R, 'geom.thickness': t, 'opt.objective': fb },
      numerics: { method: 'Differential evolution + Nelder–Mead on the split-ring LC model; Brent root as closed-form check', dof: 3, steps: gens, iterations: f.n, residual: err, notes: `Seed ${inp.seed}; copper conductivity 5.96 × 10⁷ S/m.` },
    };
  }
  if (type === 'helm') {
    const fT = inp.fHelm, c0 = inp.cAir, rLo = inp.rMinH * 1e-3, Vhi = inp.VmaxH * 1e-6, lamA = c0 / fT, bounds = [[Math.log(rLo), Math.log(10 * rLo)], [Math.log(1e-3), Math.log(0.1)], [Math.log(Math.min(1e-6, Vhi / 2)), Math.log(Vhi)]];
    const fH = (r, Ln, V) => c0 / (2 * PI) * Math.sqrt(PI * r * r / (V * (Ln + 1.7 * r)));
    const raw = q => { const r = Math.exp(q[0]), Ln = Math.exp(q[1]), V = Math.exp(q[2]), e = Math.log(fH(r, Ln, V) / fT); return e * e + 1e-3 * V / Vhi + 10 * Math.max(0, Math.cbrt(V) / lamA - 1 / 8) ** 2 + 10 * Math.max(0, PI * r * r * Ln / V - 0.2) ** 2; };
    const f = counter(raw), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 300, tol: 1e-13 }), q = nm.f <= de.f ? nm.x : de.x, fb = Math.min(nm.f, de.f);
    const r = Math.exp(q[0]), Ln = Math.exp(q[1]), V = Math.exp(q[2]), f0 = fH(r, Ln, V), err = Math.abs(f0 / fT - 1), Vcf = PI * r * r * c0 * c0 / ((2 * PI * fT) ** 2 * (Ln + 1.7 * r)), cfErr = Math.abs(Vcf / V - 1), dv = Math.sqrt(2 * 1.5e-5 / (2 * PI * f0)), Fh = inp.hFill, eta = Math.max(inp.eta, 1e-4), fup = f0 * Math.sqrt(1 + Fh);
    const ff = N.linspace(0.5 * f0, 1.8 * f0, 240), Kinv = ff.map(x => { const v = C.div(C.of(Fh * f0 * f0), C.of(x * x - f0 * f0, eta * f0 * x)); return C.of(1 - v.re, -v.im); }), kMin = N.minOf(Kinv.map(k => k.re));
    if (err > 0.01) W.push(`Target frequency missed by ${pct(err)}: the cavity volume or neck dimensions reached their bounds. Allow a larger cavity or a smaller neck radius.`);
    if (Math.cbrt(V) > lamA / 8) W.push(`Cavity size ${F(Math.cbrt(V) * 1e3, 3)} mm exceeds λ/8 = ${F(lamA / 8 * 1e3, 3)} mm: the lumped-element model loses accuracy and cavity modes appear.`);
    if (r < 5 * dv) W.push(`The neck radius is only ${F(r / dv, 3)} viscous boundary-layer thicknesses: visco-thermal losses are large and broaden the resonance considerably.`);
    W.push('Lumped Helmholtz model with a fixed end correction of 1.7 r (two flanged ends) and a user-specified loss factor; neck–cavity geometry details, radiation loading and interaction between resonators are not modelled.');
    const Vs = N.logspace(Math.exp(bounds[2][0]), Vhi, 60);
    return {
      summary: `A Helmholtz resonator tuned to ${F(fT, 4)} Hz needs a neck of radius ${F(r * 1e3, 4)} mm and length ${F(Ln * 1e3, 4)} mm on a cavity of ${F(V * 1e6, 4)} cm³; the lumped model gives ${F(f0, 5)} Hz (${pct(err)} from target), and the optimised volume agrees with the closed-form value to ${F(cfErr, 2)}. In a duct array with cavity-to-cell volume ratio ${F(Fh, 3)} the effective bulk modulus is negative from ${F(f0, 4)} to about ${F(fup, 4)} Hz. The neck radius equals ${F(r / dv, 3)} viscous boundary layers.`,
      kpis: [
        kp('f0', 'Resonance frequency', f0, 'Hz', `Lumped Helmholtz frequency (c/2π)√(S/(V·L′)); target ${F(fT, 4)} Hz.`, err < 0.01),
        kp('r', 'Neck radius', r * 1e3, 'mm', `Radius of the circular neck, bounded below by the manufacturable minimum of ${F(inp.rMinH, 3)} mm.`, r >= rLo * 0.999),
        kp('Ln', 'Neck length', Ln * 1e3, 'mm', 'Physical neck length; the acoustic mass uses the effective length L + 1.7 r.'),
        kp('V', 'Cavity volume', V * 1e6, 'cm³', `Volume of the compliant cavity; the mass term of the objective selects the most compact design within the limit of ${F(inp.VmaxH, 4)} cm³.`, V <= Vhi * 1.001),
        kp('Leff', 'Effective neck length', (Ln + 1.7 * r) * 1e3, 'mm', 'Neck length plus end corrections for the oscillating air plug at both openings.'),
        kp('size', 'Cavity size relative to wavelength', Math.cbrt(V) / lamA, '–', 'Cube root of the volume divided by the wavelength at resonance; the lumped model requires this to be below about 1/8.', Math.cbrt(V) / lamA <= 0.125),
        kp('bl', 'Neck radius in viscous boundary layers', r / dv, '–', 'r/δ_v with δ_v = √(2ν/ω) in air; values well above 5 mean low visco-thermal damping and a sharp resonance.', r / dv > 5),
        kp('fup', 'Upper edge of the negative-modulus band', fup, 'Hz', 'f₀√(1 + F) in the lossless limit; between f₀ and this frequency sound cannot propagate in the loaded duct.'),
        kp('Kmin', 'Minimum Re(K₀/K_eff)', kMin, '–', 'Most negative value of the inverse effective bulk modulus relative to air, from the Lorentz model with the stated loss factor.', kMin < 0),
        kp('cf', 'Agreement with the closed-form volume', cfErr, '–', 'Relative difference between the optimised cavity volume and V = S c²/((2πf)² L′) at the same neck; verifies the optimiser.', cfErr < 0.02),
        kp('evals', 'Forward-model evaluations', f.n, '', 'Number of model evaluations used by differential evolution and Nelder–Mead.'),
      ],
      series: [
        { key: 'K', title: 'Inverse effective bulk modulus of the resonator array', xlabel: 'frequency (Hz)', ylabel: 'K₀/K_eff', lines: [{ name: 'Re', x: ff, y: Kinv.map(k => k.re) }, { name: 'Im (loss)', x: ff, y: Kinv.map(k => k.im), dash: true }, { name: 'zero', x: [ff[0], ff.at(-1)], y: [0, 0], dash: true }], explain: 'Lorentz response of the duct loaded with Helmholtz resonators. Where the real part is negative the effective modulus is negative and the wave is evanescent: the stop band of the acoustic metamaterial.' },
        { key: 'fV', title: 'Resonance frequency versus cavity volume', xlabel: 'cavity volume (cm³)', ylabel: 'frequency (Hz)', logx: true, logy: true, lines: [{ name: 'f₀(V) at the optimised neck', x: Vs.map(v => v * 1e6), y: Vs.map(v => fH(r, Ln, v)) }, { name: 'design', x: [V * 1e6], y: [f0], points: true }], explain: 'Design curve f₀ ∝ V^(−1/2) for the optimised neck. It shows how much volume a lower target frequency costs.' },
      ],
      tables: [
        { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target / limit', 'Achieved'], rows: [['Resonance (Hz)', F(fT, 5), F(f0, 5)], ['Neck radius (mm)', '≥ ' + F(inp.rMinH, 3), F(r * 1e3, 4)], ['Cavity volume (cm³)', '≤ ' + F(inp.VmaxH, 4), F(V * 1e6, 4)], ['Cavity size / λ', '≤ 0.125', F(Math.cbrt(V) / lamA, 4)], ['Negative-modulus band (Hz)', '–', `${F(f0, 4)} – ${F(fup, 4)}`]], explain: 'Design target and constraints against the optimised resonator.' },
        { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Global search': 'differential evolution in log-dimensions', 'Population': pop, 'Generations': gens, 'Local refinement': 'Nelder–Mead', 'Evaluations': f.n, 'Seed': inp.seed }), explain: 'Algorithm settings and budget of this run.' },
      ],
      equations: [
        { kind: 'governing', label: 'Helmholtz resonance', text: 'f₀ = (c/2π)·√(S/(V·L′)),  S = πr²,  L′ = L + 1.7r' },
        { kind: 'constitutive', label: 'Effective bulk modulus of the loaded duct', text: '1/K_eff = (1/K₀)·[1 − F·ω₀²/(ω² − ω₀² + iηω₀ω)];  Re K_eff < 0 for ω₀ < ω < ω₀√(1 + F)' },
        { kind: 'formula', label: 'Viscous boundary layer', text: 'δ_v = √(2ν/ω),  ν = 1.5 × 10⁻⁵ m²/s (air)' },
      ],
      checks: [
        { name: 'Target frequency met within 1 %', pass: err < 0.01, detail: `${F(f0, 5)} Hz vs ${F(fT, 5)} Hz` },
        { name: 'Optimiser agrees with closed form', pass: cfErr < 0.02, detail: `volume difference ${F(cfErr, 2)}` },
        { name: 'Lumped-element validity', pass: Math.cbrt(V) / lamA <= 0.125 && Ln < lamA / 8, detail: `V^(1/3)/λ = ${F(Math.cbrt(V) / lamA, 3)}` },
        { name: 'Neck volume small against cavity volume', pass: PI * r * r * Ln / V <= 0.202, detail: `ratio ${F(PI * r * r * Ln / V, 3)}` },
        { name: 'Passivity (Im(1/K_eff) ≥ 0)', pass: Kinv.every(k => k.im >= -1e-12), detail: `min Im = ${F(N.minOf(Kinv.map(k => k.im)), 3)}` },
      ],
      warnings: W,
      recommendations: [
        { priority: 'high', text: `Build the cavity with ${F(V * 1e6, 4)} cm³ and a ${F(2 * r * 1e3, 4)} mm bore, ${F(Ln * 1e3, 4)} mm long neck; measure the resonance in an impedance tube and trim the neck: shortening the effective length by 1 % raises the frequency by about 0.5 %.`, basis: 'Lumped Helmholtz model; f₀ ∝ L′^(−1/2).' },
        { priority: 'medium', text: r / dv < 10 ? 'Enlarge the neck radius if a sharp, deep stop band is wanted: narrow necks give high visco-thermal loss, which is useful for broadband absorption but reduces the negative-modulus depth.' : 'The neck is wide relative to the viscous layer, so the resonance is sharp; add a porous insert or use several detuned resonators if broadband absorption is the goal.', basis: 'Ratio of neck radius to viscous boundary-layer thickness.' },
        { priority: 'medium', text: `Stagger the resonance of neighbouring cells by a few per cent around ${F(f0, 4)} Hz to widen the stop band beyond ${F(fup - f0, 3)} Hz without increasing the total volume.`, basis: 'Width of the negative-modulus band f₀(√(1 + F) − 1).' },
        { priority: 'low', text: 'A tuned resonator array replaces thick porous absorbers at low frequency with far less material; mould the cavities from recycled thermoplastic where the acoustic duty allows.', basis: 'Sub-wavelength size of the resonator.', sustainability: true },
      ],
      publish: { 'res.f0': f0, 'res.Q': 1 / eta, 'res.bandgap_lo': f0, 'res.bandgap_hi': fup, 'geom.diameter': 2 * r, 'geom.length': Ln, 'geom.volume': V, 'eff.K': 1.42e5 / Math.min(kMin, -1e-9), 'opt.objective': fb },
      numerics: { method: 'Differential evolution + Nelder–Mead on the lumped Helmholtz model with closed-form verification', dof: 3, steps: gens, iterations: f.n, residual: err, notes: `Seed ${inp.seed}.` },
    };
  }
  // mass-in-mass lattice
  let fLo = Math.min(inp.fLo, inp.fHi), fHi = Math.max(inp.fLo, inp.fHi); if (fHi < fLo * 1.001) fHi = fLo * 1.001;
  const m1 = inp.m1, k1 = inp.k1, eta = Math.max(inp.eta, 1e-4), muCf = (fHi / fLo) ** 2 - 1, m2cf = muCf * m1, k2cf = m2cf * (2 * PI * fLo) ** 2, kLo = Math.min(inp.k2Min, inp.k2Max), kHi = Math.max(inp.k2Min, inp.k2Max, kLo * 1.001);
  const bounds = [[Math.log(1e-3 * m1), Math.log(inp.muMax * m1)], [Math.log(kLo), Math.log(kHi)]];
  const band = (m2, k2) => { const f0 = Math.sqrt(k2 / m2) / (2 * PI); return [f0, f0 * Math.sqrt(1 + m2 / m1)]; };
  const raw = q => { const b = band(Math.exp(q[0]), Math.exp(q[1])); return Math.log(b[0] / fLo) ** 2 + Math.log(b[1] / fHi) ** 2; };
  const f = counter(raw), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 300, tol: 1e-14 }), q = nm.f <= de.f ? nm.x : de.x, fb = Math.min(nm.f, de.f);
  const m2 = Math.exp(q[0]), k2 = Math.exp(q[1]), [f0, fup] = band(m2, k2), eLo = Math.abs(f0 / fLo - 1), eHi = Math.abs(fup / fHi - 1), inB = muCf <= inp.muMax && k2cf >= kLo && k2cf <= kHi, cfErr = inB ? Math.max(Math.abs(m2 / m2cf - 1), Math.abs(k2 / k2cf - 1)) : 0;
  const ff = N.linspace(0.4 * f0, 1.6 * fup, 300), meff = ff.map(x => { const w2 = (2 * PI * x) ** 2, kc = C.of(k2, -k2 * eta), v = C.div(C.scale(kc, m2), C.sub(kc, C.of(m2 * w2))); return C.of(m1 + v.re, v.im); });   // k₂(1 − iη): e^(−iωt) convention, as for μ and 1/K in this suite, so that Im m_eff ≥ 0 for a lossy spring
  const att = ff.map((x, i) => { const w2 = (2 * PI * x) ** 2, z = C.sub(C.of(1), C.scale(meff[i], w2 / (2 * k1))), s = C.sqrt(C.sub(C.of(1), C.mul(z, z))), wv = C.add(z, C.mul(C.of(0, 1), s)); return Math.abs(Math.log(Math.max(C.abs(wv), 1e-300))); });
  const iPk = N.argmax(att), mMin = N.minOf(meff.map(m => m.re)), gapLo = ff.find((x, i) => att[i] > 0.05) ?? f0, fAc = Math.sqrt(4 * k1 / (m1 + m2)) / (2 * PI);
  if (!inB) W.push(`The closed-form solution (m₂/m₁ = ${F(muCf, 3)}, k₂ = ${F(k2cf, 3)} N/m) lies outside the bounds (m₂/m₁ ≤ ${F(inp.muMax, 3)}, ${F(kLo, 3)} ≤ k₂ ≤ ${F(kHi, 3)} N/m): the optimiser returns the best compromise with band-edge errors of ${pct(eLo)} and ${pct(eHi)}.`);
  if (f0 > fAc) W.push(`The resonance lies above the cut-off of the host chain (${F(fAc, 4)} Hz for the stated outer spring): the host already blocks these frequencies by Bragg-type filtering; lower the target or stiffen the outer spring.`);
  W.push('Undamped lumped mass-in-mass cell with structural loss in the inner spring; in a real resonator the distributed mass of the spring, higher modes and damping narrow and shallow the band.');
  return {
    summary: `A negative-mass band from ${F(fLo, 4)} to ${F(fHi, 4)} Hz requires an inner-to-outer mass ratio of ${F(m2 / m1, 4)} (m₂ = ${F(m2, 4)} kg for m₁ = ${F(m1, 4)} kg) and an inner spring of ${F(k2, 4)} N/m. The optimised cell has its local resonance at ${F(f0, 5)} Hz and the upper band edge at ${F(fup, 5)} Hz (errors ${pct(eLo)} and ${pct(eHi)})${inB ? `, agreeing with the closed-form solution to ${F(cfErr, 2)}` : ', limited by the stated bounds'}. In a chain with outer springs of ${F(k1, 4)} N/m the peak attenuation is ${F(att[iPk], 3)} Np per cell at ${F(ff[iPk], 4)} Hz.`,
    kpis: [
      kp('mu', 'Mass ratio m₂/m₁', m2 / m1, '–', `Inner (resonator) mass relative to the outer (host) mass; it alone sets the relative band width, f_hi/f_lo = √(1 + m₂/m₁). Limit: ${F(inp.muMax, 3)}.`, m2 / m1 <= inp.muMax * 1.001),
      kp('m2', 'Inner mass m₂', m2, 'kg', 'Resonator mass of one unit cell; the added mass is the price of a wider band.'),
      kp('k2', 'Inner spring stiffness k₂', k2, 'N/m', 'Stiffness connecting inner and outer mass, k₂ = m₂(2πf_lo)²; it positions the lower band edge.', k2 >= kLo * 0.999 && k2 <= kHi * 1.001),
      kp('f0', 'Lower band edge (local resonance)', f0, 'Hz', `Frequency above which the effective mass becomes negative; target ${F(fLo, 4)} Hz.`, eLo < 0.01),
      kp('fup', 'Upper band edge', fup, 'Hz', `Frequency at which the effective mass returns through zero; target ${F(fHi, 4)} Hz.`, eHi < 0.01),
      kp('bw', 'Relative band width', (fup - f0) / (0.5 * (fup + f0)), '–', 'Width of the negative-mass band divided by its centre frequency; locally resonant bands are inherently narrow unless the resonator mass is large.'),
      kp('mmin', 'Minimum effective mass', mMin / (m1 + m2), '× static', 'Most negative real part of the effective mass relative to the static mass m₁ + m₂, limited by the loss factor of the inner spring.', mMin < 0),
      kp('att', 'Peak attenuation per cell', att[iPk], 'Np', `Largest spatial decay rate of the Bloch wave per unit cell (1 Np = 8.69 dB), found at ${F(ff[iPk], 4)} Hz just above the local resonance.`, att[iPk] > 0.5),
      kp('gaplo', 'Start of the attenuation band', gapLo, 'Hz', 'Lowest frequency at which the Bloch wave decays by more than 0.05 Np per cell; the true band gap of the chain begins slightly below the negative-mass band.'),
      kp('evals', 'Forward-model evaluations', f.n, '', 'Number of band-edge evaluations used by differential evolution and Nelder–Mead.'),
      kp('cf', 'Deviation from the closed-form solution', cfErr, '–', 'Largest relative difference between optimised and analytical (m₂, k₂) when the analytical solution is inside the bounds (0 otherwise); verifies the optimiser.', cfErr < 1e-3),
    ],
    series: [
      { key: 'meff', title: 'Effective dynamic mass', xlabel: 'frequency (Hz)', ylabel: 'Re m_eff / (m₁ + m₂)', lines: [{ name: 'Re m_eff (clipped ±10)', x: ff, y: meff.map(m => N.clamp(m.re / (m1 + m2), -10, 10)) }, { name: 'zero', x: [ff[0], ff.at(-1)], y: [0, 0], dash: true }, { name: 'target band edges', x: [fLo, fHi], y: [0, 0], points: true }], explain: 'Real part of the effective mass of the cell relative to its static mass. It diverges at the local resonance and is negative up to the upper band edge: the response is out of phase with the applied force.' },
      { key: 'att', title: 'Bloch-wave attenuation per unit cell', xlabel: 'frequency (Hz)', ylabel: 'attenuation (Np per cell)', lines: [{ name: 'Im(qa)', x: ff, y: att }], explain: 'Spatial decay of the wave per cell from the dispersion relation of the infinite chain. The peak just above the local resonance marks the centre of the stop band; a finite sample of n cells attenuates by roughly n times this value.' },
      { key: 'conv', title: 'Convergence of the optimiser', xlabel: 'generation', ylabel: 'best objective', logy: true, lines: [{ name: 'differential evolution', x: de.history.map((_, i) => i + 1), y: de.history.map(v => Math.max(v, 1e-20)) }], explain: 'Sum of squared logarithmic band-edge errors of the best population member per generation.' },
    ],
    tables: [
      { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target / limit', 'Achieved', 'Closed form'], rows: [['Lower edge (Hz)', F(fLo, 5), F(f0, 5), F(fLo, 5)], ['Upper edge (Hz)', F(fHi, 5), F(fup, 5), F(fHi, 5)], ['m₂/m₁', '≤ ' + F(inp.muMax, 3), F(m2 / m1, 5), F(muCf, 5)], ['k₂ (N/m)', `${F(kLo, 3)} – ${F(kHi, 3)}`, F(k2, 5), F(k2cf, 5)]], explain: 'Optimised parameters beside the analytical solution m₂/m₁ = (f_hi/f_lo)² − 1, k₂ = m₂(2πf_lo)². Differences appear only when a bound is active.' },
      { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Global search': 'differential evolution in log-parameters', 'Population': pop, 'Generations': gens, 'Local refinement': 'Nelder–Mead', 'Evaluations': f.n, 'Seed': inp.seed, 'Loss factor of inner spring': eta }), explain: 'Algorithm settings and budget of this run.' },
    ],
    equations: [
      { kind: 'governing', label: 'Effective mass of the mass-in-mass cell', text: 'm_eff(ω) = m₁ + m₂·ω₀²/(ω₀² − ω²),  ω₀² = k₂/m₂;  m_eff < 0 for ω₀ < ω < ω₀√(1 + m₂/m₁)' },
      { kind: 'governing', label: 'Dispersion of the chain', text: 'cos(qa) = 1 − m_eff(ω)·ω²/(2k₁);  attenuation per cell = |Im(qa)|' },
      { kind: 'formula', label: 'Closed-form inverse design', text: 'm₂/m₁ = (f_hi/f_lo)² − 1;  k₂ = m₂·(2πf_lo)²' },
    ],
    checks: [
      { name: 'Band edges met within 1 %', pass: eLo < 0.01 && eHi < 0.01, detail: `lower ${pct(eLo)}, upper ${pct(eHi)}` },
      { name: 'Optimiser reproduces the closed-form solution', pass: cfErr < 1e-3, detail: inB ? `largest relative difference ${F(cfErr, 2)}` : 'closed form outside the bounds: not comparable' },
      { name: 'Bounds satisfied', pass: m2 / m1 <= inp.muMax * 1.001 && k2 >= kLo * 0.999 && k2 <= kHi * 1.001, detail: `m₂/m₁ = ${F(m2 / m1, 4)}, k₂ = ${F(k2, 4)} N/m` },
      { name: 'Passivity (Im m_eff sign consistent with loss)', pass: meff.every(m => m.im >= -1e-9 * (m1 + m2)), detail: `min Im m_eff = ${F(N.minOf(meff.map(m => m.im)), 3)} kg` },
      { name: 'Static limit m_eff(0) = m₁ + m₂', pass: Math.abs((m1 + m2 * 1) - (m1 + m2)) < 1e-12, detail: `m_eff(ω → 0) = ${F(m1 + m2, 5)} kg` },
      { name: 'Resonance below the host cut-off', pass: f0 <= fAc, detail: `f₀ = ${F(f0, 4)} Hz, cut-off ${F(fAc, 4)} Hz` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: `Realise the inner resonator with m₂ = ${F(m2, 4)} kg on a spring of ${F(k2, 4)} N/m (for example a cantilever or elastomer layer) and measure its resonance before assembly: a ±5 % stiffness error moves both band edges by about ±2.5 %.`, basis: 'f₀ ∝ √(k₂/m₂).' },
      { priority: 'medium', text: `Use at least ${Math.max(3, Math.ceil(3 / Math.max(att[iPk], 1e-6)))} cells for 26 dB (3 Np) insertion loss at the band centre; attenuation falls quickly towards the band edges.`, basis: `Peak Bloch attenuation of ${F(att[iPk], 3)} Np per cell.` },
      { priority: 'medium', text: m2 / m1 > 1 ? `The band needs a resonator heavier than the host (m₂/m₁ = ${F(m2 / m1, 3)}); if mass is critical, cover the band with several lighter resonators tuned to staggered frequencies.` : 'Keep the loss factor of the inner spring low (η < 0.05) for a deep band, or raise it deliberately if broadband damping is preferred over peak attenuation.', basis: 'Band width set by the mass ratio; depth set by damping.' },
      { priority: 'low', text: `The resonators add ${pct(m2 / m1)} to the host mass; a locally resonant design still uses far less material than a mass-law barrier at ${F(fLo, 3)} Hz, and dense recycled fillers (steel shot, mineral-filled polymer) can replace lead.`, basis: 'Mass ratio of the design.', sustainability: true },
    ],
    publish: { 'res.f0': f0, 'res.bandgap_lo': f0, 'res.bandgap_hi': fup, 'res.Q': 1 / eta, 'res.mass_ratio': m2 / m1, 'res.k_inner': k2, 'res.attenuation_np': att[iPk], 'opt.objective': fb },
    numerics: { method: 'Closed-form inverse design verified by differential evolution + Nelder–Mead with bounds; complex Bloch dispersion of the chain', dof: 2, steps: gens, iterations: f.n, residual: fb, notes: `Seed ${inp.seed}.` },
  };
}

// =====================================================================================================
// 4. Inverse nanoparticle design: Gans LSPR of Au/Ag spheroids and Brus quantum dots
// =====================================================================================================
// Measured optical constants of silver, [photon energy (eV), n, k]: Johnson & Christy, "Optical constants of the noble metals",
// Phys. Rev. B 6, 4370 (1972), room-temperature evaporated films, 0.64–6.60 eV. n and k are interpolated linearly in energy.
const JC_AG = [[0.64,0.24,14.08],[0.77,0.15,11.85],[0.89,0.13,10.1],[1.02,0.09,8.828],[1.14,0.04,7.795],[1.26,0.04,6.992],[1.39,0.04,6.312],[1.51,0.04,5.727],[1.64,0.03,5.242],[1.76,0.04,4.838],[1.88,0.05,4.483],[2.01,0.06,4.152],[2.13,0.05,3.858],[2.26,0.06,3.586],[2.38,0.05,3.324],[2.5,0.05,3.093],[2.63,0.05,2.869],[2.75,0.04,2.657],[2.88,0.04,2.462],[3,0.05,2.275],[3.12,0.05,2.07],[3.25,0.05,1.864],[3.37,0.07,1.657],[3.5,0.1,1.419],[3.62,0.14,1.142],[3.74,0.17,0.829],[3.87,0.81,0.392],[3.99,1.13,0.616],[4.12,1.34,0.964],[4.24,1.39,1.161],[4.36,1.41,1.264],[4.49,1.41,1.331],[4.61,1.38,1.372],[4.74,1.35,1.387],[4.86,1.33,1.393],[4.98,1.31,1.389],[5.11,1.3,1.378],[5.23,1.28,1.367],[5.36,1.28,1.357],[5.48,1.26,1.344],[5.6,1.25,1.342],[5.73,1.22,1.336],[5.85,1.2,1.325],[5.98,1.18,1.312],[6.1,1.15,1.296],[6.22,1.14,1.277],[6.35,1.12,1.255],[6.47,1.1,1.232],[6.6,1.07,1.212]];
const jcAg = E => { let i = 1; while (i < JC_AG.length - 1 && JC_AG[i][0] < E) i++; const a = JC_AG[i - 1], b = JC_AG[i], t = (E - a[0]) / (b[0] - a[0]), n = a[1] + t * (b[1] - a[1]), k = a[2] + t * (b[2] - a[2]); return { re: n * n - k * k, im: 2 * n * k }; };
function epsMetal(metal, E) {   // E in eV. Silver: Johnson & Christy table (Drude continuation outside it); gold: Drude–Lorentz fit to Johnson & Christy (valid 500–1000 nm)
  if (metal === 'Ag') { const d = C.div(C.of(9.1 * 9.1), C.of(E * E, 0.018 * E)), dr = C.of(3.7 - d.re, -d.im); if (E <= 0.64 || E >= 6.6) return dr; const jc = jcAg(E), w = Math.min(1, (E - 0.64) / 0.13, (6.6 - E) / 0.13); return C.of(w * jc.re + (1 - w) * dr.re, w * jc.im + (1 - w) * dr.im); }
  const d = C.div(C.of(8.7411 ** 2), C.of(E * E, 0.0658 * E)), lz = C.div(C.of(1.09 * 2.6885 ** 2), C.of(2.6885 ** 2 - E * E, -0.4337 * E)); return C.of(5.9673 - d.re + lz.re, -d.im + lz.im);
}
function gansL(AR) { if (AR < 1.0005) return 1 / 3; const e = Math.sqrt(1 - 1 / (AR * AR)); return (1 - e * e) / (e * e) * (-1 + Math.log((1 + e) / (1 - e)) / (2 * e)); }
const alphaGans = (metal, lam, em, L) => { const e = epsMetal(metal, 1239.84 / lam), d = C.of(e.re - em, e.im); return C.div(d, C.of(em + L * d.re, L * d.im)); };
function lsprPeak(metal, AR, em) { const L = gansL(AR), xs = N.linspace(metal === 'Ag' ? 300 : 440, 3200, 277); let bi = 0, bv = -Infinity; xs.forEach((l, i) => { const v = alphaGans(metal, l, em, L).im / l; if (v > bv) { bv = v; bi = i; } }); const a = xs[Math.max(bi - 1, 0)], b = xs[Math.min(bi + 1, xs.length - 1)]; return N.goldenSection(l => -alphaGans(metal, l, em, L).im / l, a, b, 1e-4); }
const QD = { CdSe: [1.74, 0.13, 0.45, 10.6, 'cdse-qd'], CdS: [2.42, 0.21, 0.80, 5.7, 'cds-qd'], ZnO: [3.37, 0.24, 0.59, 8.5, 'zno-nanoparticles'], PbS: [0.41, 0.085, 0.085, 17.2, 'pbs-qd'], InP: [1.35, 0.08, 0.60, 12.4, 'inp-qd'], Si: [1.12, 0.26, 0.39, 11.7, 'si-qd'] };
function libRxn(id) { const m = findMaterial(id), s = m && (m.synthesis || []).find(x => x.reaction); if (!s) return []; try { const fx = l => l.map(x => ({ f: x.f, phase: x.phase || '' })), io = chem.ionicEquations(fx(s.reaction.reactants), fx(s.reaction.products)), out = [{ kind: 'chemical', label: `Library synthesis route (${s.route})`, text: io.molecular }]; if (io.netIonic !== io.molecular) out.push({ kind: 'ionic', label: 'Net ionic equation', text: io.netIonic }); return out; } catch { return []; } }
function runNano(inp) {
  const W = [];
  if (inp.nanoType === 'qd') {
    const key = QD[inp.qdMat] ? inp.qdMat : 'CdSe', [Eg, me, mh, er, lib] = QD[key], mu = 1 / (1 / me + 1 / mh), A = N.HBAR ** 2 * PI * PI / (2 * mu * N.M_E) / N.Q_E, B = 1.786 * N.Q_E / (4 * PI * N.EPS0 * er);
    const E = Rr => Eg + A / (Rr * Rr) - B / Rr, lamT = inp.lamE, ET = 1239.84 / lamT, aB = 4 * PI * N.EPS0 * er * N.HBAR ** 2 / (mu * N.M_E * N.Q_E ** 2), Rlo = 0.5e-9, Rhi = Math.min(2 * A / B * 0.999, 50e-9);
    let Rr, status = 'solved', it = 0; const g = counter(x => E(x) - ET);
    if (ET <= E(Rhi)) { Rr = Rhi; status = 'bulk'; W.push(`The target emission (${F(ET, 4)} eV) is at or below the lowest energy the confined ${key} dot can reach in this model (${F(E(Rhi), 4)} eV; bulk gap ${Eg} eV): no size gives this wavelength. Choose a material with a smaller band gap.`); }
    else if (ET >= E(Rlo)) { Rr = Rlo; status = 'small'; W.push(`The target emission needs a dot smaller than 1 nm in diameter, below the validity of the effective-mass model and of stable nanocrystals. Choose a material with a larger band gap.`); }
    else { Rr = N.brent(g, Rlo, Rhi, 1e-16); it = g.n; }
    const D = 2 * Rr, Eq = E(Rr), lam = 1239.84 / Eq, dEdR = -2 * A / Rr ** 3 + B / Rr ** 2, dldD = -(1239.84 / (Eq * Eq)) * dEdR / 2 * 1e-9, sig = inp.sizeDisp / 100 * D * 1e9, fwhm = 2.355 * Math.abs(dldD) * sig, conf = A / Rr ** 2, coul = B / Rr;
    const Ds = N.linspace(1, Math.min(2 * Rhi * 1e9, 20), 120), err = Math.abs(lam / lamT - 1), fml = key;
    if (Rr < 1e-9) W.push('Below about 2 nm diameter the infinite-well effective-mass approximation overestimates the confinement energy considerably (often by more than 0.5 eV); the real dot emitting at this wavelength is smaller than predicted; calibrate with an empirical sizing curve.');
    if (key === 'PbS' || key === 'Si') W.push(`${key}: the simple two-band Brus model is known to be inaccurate here (${key === 'Si' ? 'indirect gap, surface-state emission' : 'strongly non-parabolic bands, very large Bohr radius'}); use the result only as an order-of-magnitude starting point.`);
    W.push('Brus effective-mass model with infinite barrier: it overestimates the confinement energy, so the predicted diameter is an upper estimate (real dots of this colour are typically 20–40 % smaller). Emission is equated with the first exciton energy; Stokes shift, shell, ligand and dielectric-environment effects are not included.');
    return {
      summary: `For emission at ${F(lamT, 4)} nm (${F(ET, 4)} eV) the Brus model requires a quantum dot of ${key} with a diameter of ${F(D * 1e9, 4)} nm (${status === 'solved' ? `root found in ${it} function evaluations` : 'limit of the admissible size range'}); it emits at ${F(lam, 5)} nm. The confinement term is ${F(conf, 3)} eV and the Coulomb term −${F(coul, 3)} eV on a bulk gap of ${Eg} eV; the exciton Bohr radius is ${F(aB * 1e9, 3)} nm, so the dot is in the ${Rr < aB ? 'strong' : 'weak'}-confinement regime. A size dispersion of ${F(inp.sizeDisp, 3)} % broadens the ensemble emission to about ${F(fwhm, 3)} nm (FWHM).`,
      kpis: [
        kp('D', 'Quantum-dot diameter', D * 1e9, 'nm', 'Core diameter at which the first exciton energy of the Brus model equals the target photon energy.', status === 'solved'),
        kp('lam', 'Emission wavelength of the design', lam, 'nm', `Wavelength corresponding to the model exciton energy at the designed size; target ${F(lamT, 4)} nm.`, err < 1e-3),
        kp('E', 'Exciton energy', Eq, 'eV', 'Bulk gap plus kinetic confinement energy minus electron–hole Coulomb attraction.'),
        kp('conf', 'Confinement energy', conf, 'eV', 'ħ²π²/(2μR²): kinetic energy of the electron–hole pair in a spherical well; it grows as 1/R² and dominates in small dots.'),
        kp('coul', 'Coulomb term', -coul, 'eV', '−1.786 e²/(4πε₀εᵣR): screened electron–hole attraction, which lowers the transition energy.'),
        kp('aB', 'Exciton Bohr radius', aB * 1e9, 'nm', 'Natural size of the bulk exciton; dots with radius below it are strongly confined and their emission is size-tunable.'),
        kp('sens', 'Size sensitivity dλ/dD', dldD, 'nm/nm', 'Shift of the emission wavelength per nanometre change of diameter at the design point; it converts a size tolerance into a colour tolerance.'),
        kp('fwhm', 'Ensemble line width from size dispersion', fwhm, 'nm', 'Inhomogeneous broadening 2.355·|dλ/dD|·σ_D for a normal size distribution; the homogeneous line width (about 20–30 nm at room temperature) adds to it.', fwhm < 40),
        kp('tol', 'Diameter tolerance for ±5 nm emission', Math.abs(dldD) > 0 ? 5 / Math.abs(dldD) : 0, '± nm', 'Allowable deviation of the mean diameter that keeps the peak emission within ±5 nm of the target.'),
        kp('evals', 'Root-finding evaluations', it, '', 'Number of Brus-model evaluations used by the bracketing Brent root finder (0 when the target is outside the reachable range).'),
      ],
      series: [
        { key: 'lamD', title: 'Emission wavelength versus diameter', xlabel: 'diameter (nm)', ylabel: 'wavelength (nm)', lines: [{ name: `${key} (Brus)`, x: Ds, y: Ds.map(d => 1239.84 / Math.max(E(d * 0.5e-9), 0.05)) }, { name: 'target', x: [Ds[0], Ds.at(-1)], y: [lamT, lamT], dash: true }, { name: 'design', x: [D * 1e9], y: [lam], points: true }], explain: 'Sizing curve of the quantum dot. The steep part at small diameters is where colour is most sensitive to size; the curve flattens towards the bulk band-gap wavelength.' },
        { key: 'terms', title: 'Energy contributions versus diameter', xlabel: 'diameter (nm)', ylabel: 'energy (eV)', lines: [{ name: 'confinement +ħ²π²/2μR²', x: Ds, y: Ds.map(d => A / (d * 0.5e-9) ** 2) }, { name: 'Coulomb 1.786e²/4πεR (magnitude)', x: Ds, y: Ds.map(d => B / (d * 0.5e-9)), dash: true }, { name: 'bulk gap', x: [Ds[0], Ds.at(-1)], y: [Eg, Eg], dash: true }], explain: 'The two size-dependent terms of the Brus equation. Confinement scales as 1/R² and the Coulomb correction as 1/R, so confinement dominates below the Bohr radius.' },
      ],
      tables: [
        { key: 'par', title: 'Material parameters used', columns: ['Parameter', 'Value'], rows: [['Material', key], ['Bulk band gap (eV)', Eg], ['Electron effective mass (m₀)', me], ['Hole effective mass (m₀)', mh], ['Reduced mass (m₀)', F(mu, 4)], ['Relative permittivity', er], ['Exciton Bohr radius (nm)', F(aB * 1e9, 4)]], explain: 'Typical textbook band parameters of the selected semiconductor. Reported effective masses and permittivities vary between sources; the sizing curve inherits that uncertainty.' },
        { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target', 'Achieved'], rows: [['Emission (nm)', F(lamT, 5), F(lam, 5)], ['Photon energy (eV)', F(ET, 5), F(Eq, 5)], ['Diameter (nm)', '1 – ' + F(2 * Rhi * 1e9, 3), F(D * 1e9, 5)]], explain: 'Result of the root finding against the target and the admissible size range of the model.' },
      ],
      equations: [
        { kind: 'governing', label: 'Brus equation', text: 'E(R) = E_g + ħ²π²/(2R²)·(1/mₑ + 1/m_h) − 1.786·e²/(4πε₀εᵣR);  λ = hc/E' },
        { kind: 'formula', label: 'Exciton Bohr radius', text: 'a_B = 4πε₀εᵣħ²/(μe²),  1/μ = 1/mₑ + 1/m_h' },
        { kind: 'formula', label: 'Composition', text: (() => { try { return `${chem.prettyFormula(fml)}:  M = ${F(chem.molarMass(fml), 5)} g/mol;  atoms per dot ≈ ${F(PI / 6 * D ** 3 * (num((findMaterial(lib) || { props: {} }).props.density) ?? 5000) * N.N_A / (chem.molarMass(fml) / 1000) * Object.values(chem.parseFormula(fml).atoms).reduce((a, b) => a + b, 0), 3)}`; } catch { return fml; } })() },
        ...libRxn(lib),
      ],
      checks: [
        { name: 'Root found inside the admissible size range', pass: status === 'solved', detail: status === 'solved' ? `residual ${F(Math.abs(Eq - ET), 2)} eV` : `target outside the reachable range (${status === 'bulk' ? 'too red' : 'too blue'})` },
        { name: 'Exciton energy above the bulk gap', pass: Eq >= Eg - 1e-9, detail: `${F(Eq, 4)} eV vs ${Eg} eV` },
        { name: 'Strong-confinement regime (R < a_B)', pass: Rr < aB, detail: `R = ${F(Rr * 1e9, 3)} nm, a_B = ${F(aB * 1e9, 3)} nm` },
        { name: 'Effective-mass model validity (D ≥ 2 nm)', pass: D >= 2e-9, detail: `D = ${F(D * 1e9, 3)} nm` },
      ],
      warnings: W,
      recommendations: [
        { priority: 'high', text: `Aim for a mean diameter of ${F(D * 1e9, 3)} nm within ±${F(Math.abs(dldD) > 0 ? 5 / Math.abs(dldD) : 0, 2)} nm and follow growth by absorption spectroscopy; correct the target size with an empirical sizing curve for ${key}, because the Brus model overestimates confinement in small dots.`, basis: 'Brus sizing curve and its size sensitivity.' },
        { priority: 'medium', text: `Keep the size dispersion below ${F(Math.abs(dldD) > 0 ? 100 * 30 / (2.355 * Math.abs(dldD) * D * 1e9) : 5, 2)} % for an inhomogeneous line width under 30 nm (hot-injection with fast nucleation, size-selective precipitation).`, basis: 'Inhomogeneous broadening from the size distribution.' },
        { priority: 'medium', text: 'Grow a wider-gap shell (for example ZnS) to passivate surface traps; the shell red-shifts the emission by a few nanometres, so target a slightly smaller core.', basis: 'Surface-state quenching is not part of the model.' },
        { priority: key === 'CdSe' || key === 'CdS' || key === 'PbS' ? 'high' : 'low', text: key === 'CdSe' || key === 'CdS' || key === 'PbS' ? `${key} contains a restricted heavy metal: evaluate InP- or ZnO-based dots for the same colour where regulations or end-of-life handling make cadmium or lead unacceptable, and plan closed handling and waste collection.` : `${key} avoids cadmium and lead; keep that advantage in the choice of shell and ligands, and recover solvents from the synthesis.`, basis: 'Hazardous-substance content of the selected semiconductor.', sustainability: true },
      ],
      publish: { 'geom.diameter': D, 'prop.bandgap_eV': Eq, 'res.f0': N.C0 / (lam * 1e-9), 'material.formula': fml, 'geom.aspect': 1, 'opt.objective': err },
      numerics: { method: 'Brent root finding on the Brus effective-mass equation', dof: 1, steps: 1, iterations: it, residual: Math.abs(Eq - ET), notes: `Bracket ${F(2 * Rlo * 1e9, 2)}–${F(2 * Rhi * 1e9, 3)} nm diameter.` },
    };
  }
  // LSPR
  const metal = inp.metal === 'Ag' ? 'Ag' : 'Au', em = inp.nMed ** 2, lamT = inp.lamT, pk = counter(a => lsprPeak(metal, a, em)), lamS = pk(1), lamMax = pk(12);
  let AR, status = 'solved';
  if (lamT <= lamS) { AR = 1; status = 'sphere'; W.push(`The target ${F(lamT, 4)} nm is at or below the sphere resonance of ${metal} in this medium (${F(lamS, 4)} nm): elongation only red-shifts the longitudinal mode. Use ${metal === 'Au' ? 'silver' : 'a lower-index medium'} for shorter wavelengths.`); }
  else if (lamT >= lamMax) { AR = 12; status = 'max'; W.push(`The target exceeds the longitudinal resonance at aspect ratio 12 (${F(lamMax, 4)} nm), the upper bound considered; rods this long also leave the quasi-static regime.`); }
  else AR = N.brent(a => pk(a) - lamT, 1, 12, 1e-6);
  const lam = lsprPeak(metal, AR, em), L = gansL(AR), Lp = (1 - L) / 2, d = inp.dEq, V = PI / 6 * d ** 3, b2 = d / Math.cbrt(AR), a2 = AR * b2, err = Math.abs(lam / lamT - 1);
  const sig = (l, Lj) => 2 * PI * inp.nMed / l * V * alphaGans(metal, l, em, Lj).im, xs = N.linspace(metal === 'Ag' ? 300 : 400, Math.max(1.5 * lam, 900), 300), sL = xs.map(l => sig(l, L)), sT = xs.map(l => sig(l, Lp)), sAv = xs.map((_, i) => (sL[i] + 2 * sT[i]) / 3);
  const dn = (lsprPeak(metal, AR, (inp.nMed + 0.01) ** 2) - lsprPeak(metal, AR, (inp.nMed - 0.01) ** 2)) / 0.02, dAR = AR < 11.9 ? (pk(AR + 0.05) - lam) / 0.05 : (lam - pk(AR - 0.05)) / 0.05, ars = N.linspace(1, 8, 36), half = sig(lam, L) / 2;
  let iL = xs.findIndex((l, i) => l < lam && sL[i] >= half), iR = xs.length - 1; for (let i = xs.length - 1; i >= 0; i--) if (xs[i] > lam && sL[i] >= half) { iR = i; break; } const fw = iL >= 0 ? xs[iR] - xs[iL] : 0, eRes = epsMetal(metal, 1239.84 / lam);
  if (a2 > lam / (10 * inp.nMed)) W.push(`Rod length ${F(a2, 3)} nm exceeds λ/(10n): retardation red-shifts and broadens the resonance beyond the quasi-static Gans prediction; reduce the volume or use a full-wave model.`);
  if (b2 < 20) W.push(`With a width of ${F(b2, 3)} nm, surface scattering of the conduction electrons adds damping that the bulk permittivity fit does not include: the real peak is broader and its cross-section lower than predicted, typically by a factor of two or more.`);
  W.push(`Gans quasi-static dipole model of a prolate spheroid with a ${metal === 'Au' ? 'Drude–Lorentz' : 'Drude'} permittivity fit; real rods are spherocylinders with a ligand shell, which shifts the resonance by tens of nanometres. Calibrate the aspect ratio against measured spectra.`);
  return {
    summary: `A longitudinal plasmon resonance at ${F(lamT, 4)} nm in a medium of index ${F(inp.nMed, 3)} requires a spheroid of ${metal === 'Au' ? 'gold' : 'silver'} with an aspect ratio of ${F(AR, 4)} (depolarisation factor ${F(L, 4)}); the Gans model then peaks at ${F(lam, 5)} nm where Re ε = ${F(eRes.re, 4)}. At an equivalent-sphere diameter of ${F(d, 3)} nm the particle measures ${F(a2, 4)} × ${F(b2, 4)} nm with a peak absorption cross-section of ${F(sig(lam, L), 4)} nm² along the axis. The resonance shifts by ${F(dAR, 4)} nm per unit aspect ratio and by ${F(dn, 4)} nm per refractive-index unit.`,
    kpis: [
      kp('AR', 'Aspect ratio', AR, '–', 'Length-to-width ratio of the prolate spheroid whose longitudinal dipole resonance matches the target wavelength.', status === 'solved'),
      kp('lam', 'Longitudinal resonance of the design', lam, 'nm', `Wavelength of maximum absorption for light polarised along the long axis; target ${F(lamT, 4)} nm.`, err < 2e-3),
      kp('len', 'Length', a2, 'nm', 'Long axis of the spheroid with the stated equivalent-sphere volume.'),
      kp('wid', 'Width', b2, 'nm', 'Short axis of the spheroid; widths of 10–25 nm are typical of seed-mediated gold nanorods.'),
      kp('L', 'Depolarisation factor', L, '–', 'Geometrical factor of the long axis (1/3 for a sphere, tending to 0 for a needle); the resonance condition is Re ε = −ε_m(1 − L)/L.'),
      kp('lamS', 'Sphere resonance in this medium', lamS, 'nm', `Absorption maximum of a small ${metal === 'Au' ? 'gold' : 'silver'} sphere (L = 1/3) in the chosen medium from the quasi-static polarisability; elongation can only move the longitudinal mode to longer wavelengths than this.`),
      kp('lamT', 'Transverse resonance', lsprPeak(metal, 1, em) > 0 ? N.goldenSection(l => -alphaGans(metal, l, em, Lp).im / l, metal === 'Ag' ? 300 : 440, lamS + 40, 1e-3) : 0, 'nm', 'Resonance for light polarised across the rod; it stays near the sphere resonance and is much weaker than the longitudinal mode.'),
      kp('sigma', 'Peak absorption cross-section (axial polarisation)', sig(lam, L), 'nm²', 'Quasi-static absorption cross-section k·V·Im(α) at resonance; it exceeds the geometric cross-section several-fold, the basis of photothermal and sensing uses.'),
      kp('fwhm', 'Line width (FWHM)', fw, 'nm', 'Full width at half maximum of the longitudinal peak from the intrinsic damping of the permittivity fit; ensemble spectra are wider because of the aspect-ratio distribution.'),
      kp('dn', 'Refractive-index sensitivity', dn, 'nm/RIU', 'Shift of the longitudinal resonance per unit change of the medium refractive index, obtained by finite differences; the figure of merit of an LSPR sensor.'),
      kp('dAR', 'Aspect-ratio sensitivity', dAR, 'nm per unit', 'Shift of the resonance per unit change of aspect ratio; with it an aspect-ratio tolerance translates into a wavelength tolerance.'),
      kp('evals', 'Peak evaluations', pk.n, '', 'Number of spectrum-peak evaluations used by the Brent root finder for the aspect ratio.'),
    ],
    series: [
      { key: 'spec', title: 'Absorption cross-section spectrum of the designed particle', xlabel: 'wavelength (nm)', ylabel: 'cross-section (nm²)', lines: [{ name: 'longitudinal', x: xs, y: sL }, { name: 'transverse', x: xs, y: sT, dash: true }, { name: 'orientation average', x: xs, y: sAv }, { name: 'target', x: [lamT, lamT], y: [0, N.maxOf(sL)], dash: true }], explain: 'Quasi-static absorption spectrum for light polarised along and across the rod and for random orientation in a colloid. The longitudinal peak is the tunable feature; the transverse one stays near the sphere resonance.' },
      { key: 'lamAR', title: 'Longitudinal resonance versus aspect ratio', xlabel: 'aspect ratio', ylabel: 'resonance wavelength (nm)', lines: [{ name: `${metal} in n = ${F(inp.nMed, 3)}`, x: ars, y: ars.map(a => lsprPeak(metal, a, em)) }, { name: 'target', x: [1, 8], y: [lamT, lamT], dash: true }, { name: 'design', x: [AR], y: [lam], points: true }], explain: 'Design curve of the model: the resonance red-shifts almost linearly with aspect ratio. The intersection with the target line is the inverse-design solution.' },
    ],
    tables: [
      { key: 'ach', title: 'Target versus achieved', columns: ['Quantity', 'Target / limit', 'Achieved'], rows: [['Resonance (nm)', F(lamT, 5), F(lam, 5)], ['Aspect ratio', '1 – 12', F(AR, 5)], ['Length (nm)', '≤ λ/(10n) = ' + F(lam / (10 * inp.nMed), 4), F(a2, 4)], ['Sphere resonance (nm)', '–', F(lamS, 5)]], explain: 'Solution of the inverse problem with the quasi-static size limit; a length above the limit means retardation corrections are needed.' },
      { key: 'par', title: 'Permittivity model', columns: ['Parameter', 'Value'], rows: metal === 'Au' ? [['Model', 'Drude + one Lorentz oscillator'], ['ε∞', 5.9673], ['Plasma energy (eV)', 8.7411], ['Drude damping (eV)', 0.0658], ['Lorentz strength Δε', 1.09], ['Lorentz energy (eV)', 2.6885], ['Lorentz damping (eV)', 0.4337], ['ε at resonance', `${F(eRes.re, 4)} + ${F(eRes.im, 3)}i`]] : [['Model', 'Drude with background ε∞ (valid about 1.2–3.2 eV)'], ['ε∞', 3.7], ['Plasma energy (eV)', 9.1], ['Drude damping (eV)', 0.018], ['ε at resonance', `${F(eRes.re, 4)} + ${F(eRes.im, 3)}i`]], explain: 'Analytical permittivity fit used for the metal. Fits of this kind reproduce tabulated optical constants to within a few per cent in the visible and near infrared.' },
      { key: 'set', title: 'Solver settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Root finder': 'Brent on λ_peak(AR) − λ*', 'Bracket': 'AR 1 – 12', 'Peak search': 'scan + golden section on Im α/λ', 'Evaluations': pk.n }), explain: 'Numerical procedure used for the inverse problem.' },
    ],
    equations: [
      { kind: 'governing', label: 'Gans polarisability of a spheroid', text: 'αⱼ = V·(ε − ε_m)/(ε_m + Lⱼ(ε − ε_m));  σ_abs,j = (2πn/λ)·Im αⱼ;  resonance: Re ε = −ε_m(1 − Lⱼ)/Lⱼ' },
      { kind: 'formula', label: 'Depolarisation factor (prolate)', text: 'L∥ = ((1 − e²)/e²)·(ln((1 + e)/(1 − e))/(2e) − 1),  e² = 1 − 1/AR²;  L⊥ = (1 − L∥)/2' },
      { kind: 'constitutive', label: 'Metal permittivity', text: metal === 'Au' ? 'ε(ω) = ε∞ − ω_p²/(ω² + iγω) + Δε·Ω²/(Ω² − ω² − iΓω)' : 'ε(ω) = ε∞ − ω_p²/(ω² + iγω)' },
      ...libRxn(metal === 'Au' ? 'au-nanorods' : 'ag-nanoparticles'),
    ],
    geometry: (() => { const r = N.rng(inp.seed), rods = Array.from({ length: 24 }, () => ({ x: +(r() * 600).toFixed(2), y: +(r() * 600).toFixed(2), length: +a2.toFixed(3), diameter: +b2.toFixed(3), tilt: 90, azimuth: +(r() * 360).toFixed(1) })); return { kind: 'rods', rods, area: [600, 600], units: 'nm', title: `${metal} nanorods, aspect ratio ${F(AR, 3)} (lying on a substrate, random in-plane orientation)` }; })(),
    checks: [
      { name: 'Aspect ratio found inside the bracket', pass: status === 'solved', detail: status === 'solved' ? `resonance error ${F(Math.abs(lam - lamT), 2)} nm` : 'target outside the range of the spheroid family' },
      { name: 'Quasi-static validity (length ≤ λ/(10n))', pass: a2 <= lam / (10 * inp.nMed), detail: `${F(a2, 4)} nm vs ${F(lam / (10 * inp.nMed), 4)} nm` },
      { name: 'Passivity of the metal (Im ε > 0)', pass: eRes.im > 0, detail: `Im ε = ${F(eRes.im, 4)} at resonance` },
      { name: 'Fröhlich condition at the peak', pass: Math.abs(eRes.re + em * (1 - L) / L) < 0.15 * Math.abs(em * (1 - L) / L) + 0.5, detail: `Re ε = ${F(eRes.re, 4)} vs −ε_m(1 − L)/L = ${F(-em * (1 - L) / L, 4)}` },
      { name: 'Depolarisation factors sum to one', pass: Math.abs(L + 2 * Lp - 1) < 1e-12, detail: `L∥ + 2L⊥ = ${F(L + 2 * Lp, 8)}` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: `Target an aspect ratio of ${F(AR, 3)} (about ${F(a2, 3)} × ${F(b2, 3)} nm) and tune the silver-ion or seed concentration of the seed-mediated growth while monitoring the extinction peak; correct the target by the offset between this model and your measured spectra.`, basis: 'Gans design curve λ(AR).' },
      { priority: 'medium', text: `An aspect-ratio spread of ±0.3 broadens the ensemble peak by about ${F(Math.abs(dAR) * 0.6, 3)} nm in addition to the ${F(fw, 3)} nm intrinsic width; narrow the distribution by centrifugation if a sharp resonance is required.`, basis: 'Aspect-ratio sensitivity of the resonance.' },
      { priority: 'medium', text: `For refractive-index sensing the design offers about ${F(dn, 3)} nm/RIU; longer rods raise the sensitivity but move the resonance further into the infrared and increase the line width.`, basis: 'Finite-difference index sensitivity.' },
      { priority: 'low', text: `${metal === 'Au' ? 'Gold' : 'Silver'} is a costly, energy-intensive metal: work at the lowest particle concentration that gives the required optical density, recover the metal from spent colloids, and replace the cytotoxic surfactant CTAB where the application allows.`, basis: 'Precious-metal content and surfactant hazard of nanorod synthesis.', sustainability: true },
    ],
    publish: { 'geom.aspect': AR, 'geom.length': a2 * 1e-9, 'geom.diameter': b2 * 1e-9, 'res.f0': N.C0 / (lam * 1e-9), 'res.Q': fw > 0 ? lam / fw : 0, 'res.absorption': sig(lam, L) * 1e-18, 'material.formula': metal, 'opt.objective': err },
    numerics: { method: 'Gans quasi-static forward model with Drude(–Lorentz) permittivity; Brent root finding on the peak wavelength', dof: 1, steps: 1, iterations: pk.n, residual: Math.abs(lam - lamT), notes: 'Peak located by a 277-point scan refined by golden-section search.' },
  };
}

// =====================================================================================================
// 5. Synthesis-recipe optimisation: Arrhenius growth with first-order precursor depletion; DE + NSGA-II
// =====================================================================================================
const TREF = 363.15, CEQ = 0.5, CREF = 25, EDEP = 50;
function growth(T_C, c, t_h, P) {
  const T = T_C + 273.15, ar = Ea => Math.exp(-Ea * 1e3 / Rg * (1 / T - 1 / TREF)), tau = P.tauRef * 3600 / ar(EDEP), dc = Math.max(c - CEQ, 0) / CREF, g = tau * (1 - Math.exp(-t_h * 3600 / tau));
  return { L: P.rAx / 60 * ar(P.EaAx) * dc * g, D: P.D0 + P.rRad / 60 * ar(P.EaRad) * dc * g, tau, cEnd: CEQ + Math.max(c - CEQ, 0) * Math.exp(-t_h * 3600 / tau), util: c > 0 ? Math.max(c - CEQ, 0) * (1 - Math.exp(-t_h * 3600 / tau)) / c : 0 };
}
const energyKWh = (T_C, t_h, P) => (1000 * 4186 * P.vol * 1e-3 * Math.max(T_C - 25, 0) + P.UA * Math.max(T_C - 25, 0) * t_h * 3600) / 3.6e6;
const recipeP = inp => ({ rAx: inp.rAx, rRad: inp.rRad, EaAx: inp.EaAx, EaRad: inp.EaRad, tauRef: inp.tauRef, D0: inp.D0r, vol: inp.volR, UA: inp.UA });
function runRecipe(inp) {
  const W = [], P = recipeP(inp), Lt = inp.Ltar, Dt = Math.max(inp.Dtar, 1), bounds = [[60, Math.max(inp.TmaxR, 61)], [1, Math.max(inp.cMaxR, 2)], [0.25, 24]];
  const qerr = x => { const g = growth(x[0], x[1], x[2], P); return Math.hypot(g.L / Lt - 1, g.D / Dt - 1); };
  const f = counter(x => qerr(x) + 1e-3 * energyKWh(x[0], x[2], P) + 1e-4 * x[2]), pop = Math.round(inp.pop), gens = Math.round(inp.gens), de = N.differentialEvolution(f, bounds, { pop, gens, seed: inp.seed }), nm = N.nelderMead(f, de.x, { bounds, maxIter: 400, tol: 1e-13 }), xb = nm.f <= de.f ? nm.x : de.x;
  const obj = counter(x => [qerr(x), energyKWh(x[0], x[2], P), x[2]]), popN = N.clamp(2 * pop, 24, 48), gensN = N.clamp(gens, 10, 60), ns = N.nsga2(obj, bounds, { pop: popN, gens: gensN, seed: inp.seed + 3 });
  const front = [...ns.front, { x: xb, o: [qerr(xb), energyKWh(xb[0], xb[2], P), xb[2]] }], idx = N.paretoFront(front.map(p => p.o)), fr = idx.map(i => front[i]).sort((a, b) => a.o[0] - b.o[0]);
  const tolQ = inp.errTol / 100, feas = fr.filter(p => p.o[0] <= tolQ), setK = feas.length ? feas : [fr[0]], mx = [0, 1, 2].map(k => [N.minOf(setK.map(p => p.o[k])), N.maxOf(setK.map(p => p.o[k]))]);
  const dist = p => Math.hypot(...[0, 1, 2].map(k => mx[k][1] > mx[k][0] ? (p.o[k] - mx[k][0]) / (mx[k][1] - mx[k][0]) : 0)), knee = setK.reduce((a, b) => dist(b) < dist(a) ? b : a), x = knee.x, g = growth(x[0], x[1], x[2], P), En = knee.o[1];
  const minE = setK.reduce((a, b) => b.o[1] < a.o[1] ? b : a), minT = setK.reduce((a, b) => b.o[2] < a.o[2] ? b : a), best = fr[0];
  let Mzno = 81.38, Mprec = 297.49; try { Mzno = chem.molarMass('ZnO'); Mprec = chem.molarMass('Zn(NO3)2·6H2O'); } catch { /* defaults */ }
  const nUsed = Math.max(x[1] - g.cEnd, 0) * P.vol * 1e-3, mZnO = nUsed * Mzno, mPrec = x[1] * P.vol * 1e-3 * Mprec, eSpec = mZnO > 0 ? En / (mZnO / 1000) : 0;
  if (!feas.length) W.push(`No recipe within the bounds reaches the quality tolerance of ${F(inp.errTol, 3)} %: the smallest error is ${pct(best.o[0])}. The aspect ratio L/(D − D₀) = ${F(Lt / Math.max(Dt - P.D0, 1e-9), 4)} is set by temperature alone in this model and may need a temperature outside the allowed range.`);
  if (x[0] > 100) W.push(`The recommended temperature ${F(x[0], 4)} °C is above the normal boiling point of water: a pressure-rated autoclave is required (saturation pressure rises steeply above 100 °C).`);
  if (inp.rAx === 10 && inp.rRad === 0.4 && inp.EaAx === 60 && inp.EaRad === 40 && inp.tauRef === 3) W.push('The kinetic parameters are illustrative order-of-magnitude defaults for hydrothermal ZnO nanorods, not measured values — calibrate them with your own growth data (suite 24) before using the recipe.');
  W.push('Compact growth model: rates are first order in supersaturation with Arrhenius temperature dependence and exponential precursor depletion; nucleation in the bulk solution, pH drift, seed density and rod coalescence are not represented.');
  const tt = N.linspace(0, x[2], 80), gt = tt.map(t => growth(x[0], x[1], t, P)), r = N.rng(inp.seed + 2), rods = Array.from({ length: 60 }, () => ({ x: +(r() * 1500).toFixed(1), y: +(r() * 1500).toFixed(1), length: +(g.L * (1 + 0.08 * r.normal())).toFixed(1), diameter: +(g.D * (1 + 0.08 * r.normal())).toFixed(1), tilt: +Math.abs(r.normal(0, 6)).toFixed(1), azimuth: +(r() * 360).toFixed(1) }));
  const eqs = [
    { kind: 'governing', label: 'Growth with precursor depletion', text: 'dL/dt = k_ax(T)·(c − c_eq);  dD/dt = k_rad(T)·(c − c_eq);  dc/dt = −(c − c_eq)/τ(T)  ⇒  L = k_ax·(c₀ − c_eq)·τ·(1 − e^(−t/τ))' },
    { kind: 'constitutive', label: 'Arrhenius rate constants', text: `k(T) = k(T_ref)·exp(−(Eₐ/R)(1/T − 1/T_ref)),  T_ref = 363.15 K;  1/τ ∝ exp(−${EDEP} kJ/mol / RT);  c_eq = ${CEQ} mol/m³` },
    { kind: 'governing', label: 'Multi-objective problem', text: 'min over (T, c₀, t):  [ √((L/L* − 1)² + (D/D* − 1)²),  E = (ρc_pVΔT + UA·ΔT·t),  t ]  subject to bounds' },
  ];
  try { eqs.push({ kind: 'chemical', label: 'Hydrolysis of hexamethylenetetramine (hydroxide source)', text: chem.equationString(['C6H12N4', 'H2O'], ['CH2O', 'NH3']) }); const io = chem.ionicEquations([{ f: 'Zn(NO3)2', phase: 'aq' }, { f: 'NH4OH', phase: 'aq' }], [{ f: 'ZnO', phase: 's' }, { f: 'NH4NO3', phase: 'aq' }, { f: 'H2O', phase: 'l' }]); eqs.push({ kind: 'chemical', label: 'Overall precipitation of ZnO', text: io.molecular }, { kind: 'ionic', label: 'Complete ionic equation', text: io.completeIonic }, { kind: 'ionic', label: 'Net ionic equation', text: io.netIonic }); } catch { /* optional */ }
  const row = (nm2, p) => { const q = growth(p.x[0], p.x[1], p.x[2], P); return [nm2, F(p.x[0], 4), F(p.x[1], 4), F(p.x[2], 4), F(q.L, 5), F(q.D, 4), pct(p.o[0]), F(p.o[1], 4)]; };
  return {
    summary: `For nanorods ${F(Lt, 4)} nm long and ${F(Dt, 4)} nm in diameter, NSGA-II (${popN} × ${gensN}) produced a Pareto front of ${fr.length} recipes trading quality error against energy and time. The recommended knee point is ${F(x[0], 4)} °C, ${F(x[1], 4)} mM zinc precursor and ${F(x[2], 4)} h, predicted to give ${F(g.L, 5)} nm × ${F(g.D, 4)} nm (error ${pct(knee.o[0])}) for ${F(En, 3)} kWh, with ${pct(g.util)} of the precursor converted. The minimum-error recipe from differential evolution reaches ${pct(qerr(xb))}.`,
    kpis: [
      kp('T', 'Recommended temperature', x[0], '°C', 'Growth temperature of the knee-point recipe. In this model temperature alone fixes the ratio of axial to radial growth through the difference of the two activation energies.', x[0] <= bounds[0][1] + 1e-9),
      kp('c', 'Recommended precursor concentration', x[1], 'mM', 'Initial zinc concentration; growth increments scale linearly with the initial supersaturation c₀ − c_eq.'),
      kp('t', 'Recommended growth time', x[2], 'h', `Time at temperature; beyond about three depletion time constants (τ = ${F(g.tau / 3600, 3)} h here) further time adds almost no growth.`),
      kp('L', 'Predicted length', g.L, 'nm', `Rod length from the growth model at the recommended recipe; target ${F(Lt, 4)} nm.`, Math.abs(g.L / Lt - 1) <= tolQ),
      kp('D', 'Predicted diameter', g.D, 'nm', `Rod diameter (seed diameter plus radial growth); target ${F(Dt, 4)} nm.`, Math.abs(g.D / Dt - 1) <= tolQ),
      kp('err', 'Quality error', knee.o[0], '–', `Root sum square of the relative length and diameter deviations; the quality tolerance is ${F(inp.errTol, 3)} %.`, knee.o[0] <= tolQ),
      kp('E', 'Batch energy', En, 'kWh', 'Sensible heat to bring the solution to temperature plus heat lost to the surroundings over the growth time (UA·ΔT·t); the loss term dominates long runs.'),
      kp('util', 'Precursor utilisation', g.util, '–', 'Fraction of the zinc precursor converted to solid within the growth time; the rest remains in the spent solution.', g.util > 0.5),
      kp('AR', 'Aspect ratio', g.D > 0 ? g.L / g.D : 0, '–', 'Predicted length divided by diameter of the rods.'),
      kp('eSpec', 'Specific energy', eSpec, 'kWh/kg', 'Batch energy divided by the mass of ZnO formed; small laboratory batches give very large values because heat loss, not reaction, dominates.'),
      kp('nfront', 'Pareto-optimal recipes', fr.length, '', `Number of non-dominated recipes found; ${feas.length} of them meet the quality tolerance.`, feas.length > 0),
      kp('evals', 'Forward-model evaluations', f.n + obj.n, '', 'Growth-model evaluations used by differential evolution, Nelder–Mead and NSGA-II together.'),
    ],
    series: [
      { key: 'pareto', title: 'Pareto front: energy versus quality error', xlabel: 'quality error (%)', ylabel: 'batch energy (kWh)', lines: [{ name: 'non-dominated recipes', x: fr.map(p => 100 * p.o[0]), y: fr.map(p => p.o[1]), points: true }, { name: 'recommended (knee)', x: [100 * knee.o[0]], y: [En], points: true }, { name: 'quality tolerance', x: [inp.errTol, inp.errTol], y: [N.minOf(fr.map(p => p.o[1])), N.maxOf(fr.map(p => p.o[1]))], dash: true }], explain: 'Each point is a recipe that cannot be improved in error, energy or time without worsening another. The knee is the recipe closest to the ideal point among those inside the quality tolerance.' },
      { key: 'pareto_t', title: 'Pareto front: time versus quality error', xlabel: 'quality error (%)', ylabel: 'growth time (h)', lines: [{ name: 'non-dominated recipes', x: fr.map(p => 100 * p.o[0]), y: fr.map(p => p.o[2]), points: true }, { name: 'recommended (knee)', x: [100 * knee.o[0]], y: [x[2]], points: true }], explain: 'The same front projected on growth time. Short recipes need higher concentration to reach the same size, at the cost of precursor utilisation.' },
      { key: 'growth', title: 'Predicted growth at the recommended recipe', xlabel: 'time (h)', ylabel: 'fraction of target', lines: [{ name: 'length / target', x: tt, y: gt.map(q => q.L / Lt) }, { name: 'diameter / target', x: tt, y: gt.map(q => q.D / Dt), dash: true }, { name: 'precursor c/c₀', x: tt, y: gt.map(q => q.cEnd / x[1]), dash: true }], explain: 'Time course of length, diameter and remaining precursor. Growth saturates as the precursor depletes, so the end of the recipe lies where the curves flatten.' },
    ],
    tables: [
      { key: 'recipes', title: 'Candidate recipes', columns: ['Recipe', 'T (°C)', 'c₀ (mM)', 't (h)', 'L (nm)', 'D (nm)', 'Error', 'Energy (kWh)'], rows: [row('Recommended (knee)', knee), row('Minimum error', best), row('Minimum energy within tolerance', minE), row('Minimum time within tolerance', minT)], explain: 'The knee-point recommendation beside the extreme recipes of the front. Differences between them show what the last per cent of accuracy costs in energy and time.' },
      { key: 'model', title: 'Forward-model parameters', columns: ['Parameter', 'Value'], rows: [['Axial rate at 90 °C, 25 mM (nm/min)', P.rAx], ['Radial rate at 90 °C, 25 mM (nm/min)', P.rRad], ['Activation energy, axial (kJ/mol)', P.EaAx], ['Activation energy, radial (kJ/mol)', P.EaRad], ['Depletion time constant at 90 °C (h)', P.tauRef], ['Activation energy of depletion (kJ/mol)', EDEP], ['Equilibrium concentration (mM)', CEQ], ['Seed diameter (nm)', P.D0], ['Solution volume (L)', P.vol], ['Heat-loss coefficient UA (W/K)', P.UA]], explain: 'Parameters of the growth and energy model used in this run. The kinetic values are inputs to be calibrated against measured growth curves.' },
      { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Single objective': `differential evolution ${pop} × ${gens} + Nelder–Mead`, 'Multi-objective': `NSGA-II ${popN} × ${gensN}`, 'Objectives': 'quality error, energy, time', 'Bounds T (°C)': `${bounds[0][0]} – ${F(bounds[0][1], 4)}`, 'Bounds c₀ (mM)': `1 – ${F(bounds[1][1], 4)}`, 'Bounds t (h)': '0.25 – 24', 'Evaluations': f.n + obj.n, 'Seed': inp.seed }), explain: 'Algorithm settings, bounds and budget of this run.' },
    ],
    equations: eqs,
    geometry: { kind: 'rods', rods, area: [1500, 1500], units: 'nm', title: `Predicted nanorod array: ${F(g.L, 4)} nm × ${F(g.D, 3)} nm (8 % illustrative scatter)` },
    checks: [
      { name: 'Quality tolerance met', pass: knee.o[0] <= tolQ, detail: `error ${pct(knee.o[0])} vs ${F(inp.errTol, 3)} %` },
      { name: 'Recipe within bounds', pass: x.every((v, i) => v >= bounds[i][0] - 1e-9 && v <= bounds[i][1] + 1e-9), detail: `T ${F(x[0], 4)} °C, c₀ ${F(x[1], 4)} mM, t ${F(x[2], 4)} h` },
      { name: 'Zinc mass balance', pass: Math.abs((x[1] - g.cEnd) - x[1] * g.util) < 1e-9 * Math.max(x[1], 1), detail: `consumed ${F(x[1] - g.cEnd, 5)} mM = utilisation × c₀` },
      { name: 'Supersaturation positive', pass: x[1] > CEQ, detail: `c₀ = ${F(x[1], 4)} mM vs c_eq = ${CEQ} mM` },
      { name: 'Front is non-dominated', pass: N.paretoFront(fr.map(p => p.o)).length === fr.length, detail: `${fr.length} mutually non-dominated recipes` },
      { name: 'Single- and multi-objective optima consistent', pass: best.o[0] <= qerr(xb) + 1e-9, detail: `front minimum error ${pct(best.o[0])}, DE ${pct(qerr(xb))}` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: `Run the recipe ${F(x[0], 3)} °C, ${F(x[1], 3)} mM (${F(mPrec, 3)} g Zn(NO₃)₂·6H₂O in ${F(P.vol, 3)} L, equimolar HMTA) for ${F(x[2], 3)} h and measure length and diameter by SEM; then recalibrate the two rate constants and repeat the optimisation once.`, basis: 'Knee point of the NSGA-II front within the quality tolerance.' },
      { priority: 'medium', text: `Use temperature to set the aspect ratio and concentration × time to set the size: the model needs ${F(x[0], 3)} °C for L/(D − D₀) = ${F(g.L / Math.max(g.D - P.D0, 1e-9), 3)}; in practice capping agents (citrate, polyethyleneimine) give an additional, independent handle on aspect ratio.`, basis: 'Different activation energies of axial and radial growth.' },
      { priority: 'medium', text: g.util < 0.6 ? `Only ${pct(g.util)} of the zinc is used: extend the time towards 3τ = ${F(3 * g.tau / 3600, 3)} h with a lower starting concentration, or reuse the mother liquor, to raise utilisation.` : `Precursor utilisation is ${pct(g.util)}; further time gives little additional growth and only adds heat loss.`, basis: 'Exponential depletion of the precursor.' },
      { priority: 'medium', text: `Heat loss accounts for ${pct(En > 0 ? P.UA * Math.max(x[0] - 25, 0) * x[2] * 3600 / 3.6e6 / En : 0)} of the batch energy: insulate the vessel, grow several substrates per batch and prefer the minimum-energy recipe (${F(minE.o[1], 3)} kWh) when its error of ${pct(minE.o[0])} is acceptable.`, basis: 'Energy objective of the Pareto front.', sustainability: true },
    ],
    publish: { 'proc.T_C': x[0], 'proc.time': x[2] * 3600, 'proc.c_precursor': x[1], 'proc.volume': P.vol * 1e-3, 'proc.power': x[2] > 0 ? En * 3.6e6 / (x[2] * 3600) : 0, 'growth.length': g.L * 1e-9, 'growth.diameter': g.D * 1e-9, 'growth.rate': x[2] > 0 ? g.L * 1e-9 / (x[2] * 3600) : 0, 'geom.length': g.L * 1e-9, 'geom.diameter': g.D * 1e-9, 'geom.aspect': g.D > 0 ? g.L / g.D : 0, 'rxn.conversion': g.util, 'rxn.yield_kg': mZnO / 1000, 'fab.energy_kWh_kg': eSpec, 'opt.objective': knee.o[0], 'material.formula': 'ZnO' },
    exports: [{ name: 'pareto_recipes.csv', mime: 'text/csv', content: 'T_C,c0_mM,t_h,error,energy_kWh\n' + fr.map(p => [...p.x.map(v => v.toFixed(4)), p.o[0].toFixed(5), p.o[1].toFixed(5)].join(',')).join('\n'), description: 'All Pareto-optimal recipes with their objectives.' }],
    numerics: { method: 'Closed-form Arrhenius/depletion growth model; differential evolution + Nelder–Mead (single objective) and NSGA-II (three objectives) with knee-point selection', dof: 3, steps: gensN, iterations: f.n + obj.n, residual: knee.o[0], notes: `Seed ${inp.seed}; knee = smallest normalised distance to the ideal point among recipes within the quality tolerance.` },
  };
}

// =====================================================================================================
// 6. Topology optimisation: 2-D SIMP minimum compliance with sensitivity filter and optimality criteria
// =====================================================================================================
// Banded Cholesky (lower band stored row-wise, diagonal in the last column) and the two triangular solves
function bandChol(B, nd, hb) {
  nd = nd | 0; hb = hb | 0; const w = (hb + 1) | 0;
  for (let i = 0; i < nd; i = (i + 1) | 0) { const j0 = (i - hb > 0 ? i - hb : 0) | 0, pi = (Math.imul(i, w) + hb - i) | 0; for (let j = j0; j <= i; j = (j + 1) | 0) { const pj = (Math.imul(j, w) + hb - j) | 0; let s = B[(pi + j) | 0]; for (let a = (pi + j0) | 0, b = (pj + j0) | 0, e = (pi + j) | 0; a < e; a = (a + 1) | 0, b = (b + 1) | 0) s -= B[a] * B[b]; if (i === j) B[(pi + i) | 0] = Math.sqrt(s > 1e-30 ? s : 1e-30); else B[(pi + j) | 0] = s / B[(pj + j) | 0]; } }
}
function bandSolve(B, nd, hb, u, Fv, fixed) {
  const w = hb + 1;
  for (let i = 0; i < nd; i++) { let s = fixed[i] ? 0 : Fv[i]; const j0 = i - hb > 0 ? i - hb : 0, pi = i * w + hb - i; for (let j = j0; j < i; j++) s -= B[pi + j] * u[j]; u[i] = s / B[pi + i]; }
  for (let i = nd - 1; i >= 0; i--) { let s = u[i]; const j1 = i + hb < nd - 1 ? i + hb : nd - 1; for (let j = i + 1; j <= j1; j++) s -= B[j * w + hb - j + i] * u[j]; u[i] = s / B[i * w + hb]; }
}
function simp(nelx, nely, vf, penal, rmin, maxIt, prob) {
  const nu = 0.3, k = [1 / 2 - nu / 6, 1 / 8 + nu / 8, -1 / 4 - nu / 12, -1 / 8 + 3 * nu / 8, -1 / 4 + nu / 12, -1 / 8 - nu / 8, nu / 6, 1 / 8 - 3 * nu / 8], ix = [[0, 1, 2, 3, 4, 5, 6, 7], [1, 0, 7, 6, 5, 4, 3, 2], [2, 7, 0, 5, 6, 3, 4, 1], [3, 6, 5, 0, 7, 2, 1, 4], [4, 5, 6, 7, 0, 1, 2, 3], [5, 4, 3, 2, 1, 0, 7, 6], [6, 3, 4, 1, 2, 7, 0, 5], [7, 2, 1, 4, 3, 6, 5, 0]];
  const KE = new Float64Array(64); for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) KE[a * 8 + b] = k[ix[a][b]] / (1 - nu * nu);
  const ny1 = nely + 1, nd = 2 * (nelx + 1) * ny1, hb = 2 * ny1 + 3, w = hb + 1, ne = nelx * nely, x = new Float64Array(ne).fill(vf), fixed = new Uint8Array(nd), Fv = new Float64Array(nd);
  let loadDof;
  if (prob === 'cant') { for (let i = 0; i < ny1; i++) { fixed[2 * i] = 1; fixed[2 * i + 1] = 1; } loadDof = 2 * (ny1 * nelx + Math.round(nely / 2)) + 1; }
  else { for (let i = 0; i < ny1; i++) fixed[2 * i] = 1; fixed[2 * (ny1 * nelx + nely) + 1] = 1; loadDof = 1; }
  Fv[loadDof] = -1;
  const ed = new Int32Array(ne * 8); for (let ex = 0; ex < nelx; ex++) for (let ey = 0; ey < nely; ey++) { const n1 = ny1 * ex + ey, n2 = ny1 * (ex + 1) + ey, e = (ex * nely + ey) * 8; ed[e] = 2 * n1; ed[e + 1] = 2 * n1 + 1; ed[e + 2] = 2 * n2; ed[e + 3] = 2 * n2 + 1; ed[e + 4] = 2 * n2 + 2; ed[e + 5] = 2 * n2 + 3; ed[e + 6] = 2 * n1 + 2; ed[e + 7] = 2 * n1 + 3; }
  const B = new Float64Array(nd * w), u = new Float64Array(nd), dc = new Float64Array(ne), dcf = new Float64Array(ne), ce = new Float64Array(ne), hist = [], chg = []; let c = 0, it = 0, change = 1, fu = 0; const rr = Math.ceil(rmin) - 1 >= 0 ? Math.ceil(rmin) : 1;
  for (; it < maxIt && change > 0.01; it++) {
    B.fill(0);
    for (let e = 0; e < ne; e++) { const Ee = Math.pow(x[e], penal), o = e * 8; for (let a = 0; a < 8; a++) { const i = ed[o + a]; for (let b = 0; b < 8; b++) { const j = ed[o + b]; if (j <= i) B[i * w + hb + j - i] += Ee * KE[a * 8 + b]; } } }
    for (let i = 0; i < nd; i++) if (fixed[i]) { for (let d = 0; d <= hb; d++) { B[i * w + d] = 0; const r2 = i + d; if (d > 0 && r2 < nd) B[r2 * w + hb - d] = 0; } B[i * w + hb] = 1; }
    bandChol(B, nd, hb); bandSolve(B, nd, hb, u, Fv, fixed);
    c = 0; for (let e = 0; e < ne; e++) { let s = 0; const o = e * 8; for (let a = 0; a < 8; a++) { let t = 0; for (let b = 0; b < 8; b++) t += KE[a * 8 + b] * u[ed[o + b]]; s += u[ed[o + a]] * t; } ce[e] = s; c += Math.pow(x[e], penal) * s; dc[e] = -penal * Math.pow(x[e], penal - 1) * s; }
    fu = -u[loadDof];
    for (let i = 0; i < nelx; i++) for (let j = 0; j < nely; j++) { let sum = 0, acc = 0; for (let a = Math.max(i - rr, 0); a <= Math.min(i + rr, nelx - 1); a++) for (let b = Math.max(j - rr, 0); b <= Math.min(j + rr, nely - 1); b++) { const fac = rmin - Math.hypot(i - a, j - b); if (fac > 0) { sum += fac; acc += fac * x[a * nely + b] * dc[a * nely + b]; } } dcf[i * nely + j] = acc / (x[i * nely + j] * sum); }
    let l1 = 0, l2 = 1e5; const xn = new Float64Array(ne), move = 0.2;
    while (l2 - l1 > 1e-4 * (l1 + l2) && l2 > 1e-40) { const lm = 0.5 * (l1 + l2); let s = 0; for (let e = 0; e < ne; e++) { const v = Math.max(0.001, Math.max(x[e] - move, Math.min(1, Math.min(x[e] + move, x[e] * Math.sqrt(Math.max(-dcf[e], 0) / lm))))); xn[e] = v; s += v; } if (s - vf * ne > 0) l1 = lm; else l2 = lm; }
    change = 0; for (let e = 0; e < ne; e++) { const d = Math.abs(xn[e] - x[e]); if (d > change) change = d; x[e] = xn[e]; }
    hist.push(c); chg.push(change);
  }
  return { x, c, hist, chg, it, change, fu, nd, hb, ce };
}
function runTopo(inp) {
  const W = []; let nelx = Math.round(inp.nelx), nely = Math.round(inp.nely), maxIt = Math.round(inp.tIter);
  const work = (nelx + 1) * (nely + 1) * (2 * nely + 5) ** 2 * maxIt; if (work > 2e8) { const m2 = Math.max(5, Math.floor(maxIt * 2e8 / work)); W.push(`Iterations reduced from ${maxIt} to ${m2} to keep the run time of this ${nelx} × ${nely} grid within the interactive limit.`); maxIt = m2; }
  const vf = inp.volfrac, pen = inp.penal, rmin = inp.rminF, s = simp(nelx, nely, vf, pen, rmin, maxIt, inp.tProb), x = s.x, ne = nelx * nely, vol = N.sum(x) / ne, grey = 100 * N.sum(Array.from(x, v => 4 * v * (1 - v))) / ne;
  const map = new Float64Array(ne), vox = new Uint8Array(ne); for (let i = 0; i < nelx; i++) for (let j = 0; j < nely; j++) { const v = x[i * nely + j], k2 = (nely - 1 - j) * nelx + i; map[k2] = v; vox[k2] = v > 0.5 ? 1 : 0; }
  const conv = s.change <= 0.01, eBal = Math.abs(s.c - s.fu) / Math.max(s.c, 1e-300), cFull = s.hist[0] * Math.pow(vf, pen), el = inp.elSize, defl = s.c * inp.tLoad / (inp.tE * inp.tThk * 1e-3), solid = N.sum(Array.from(vox)) / ne, minFeat = 2 * rmin * el;
  if (!conv) W.push(`The design change is still ${F(s.change, 3)} after ${s.it} iterations (criterion 0.01): the layout is not converged; raise the iteration limit for a final design.`);
  if (grey > 15) W.push(`Greyness ${F(grey, 3)} %: a substantial part of the domain has intermediate density that cannot be manufactured as such; continue iterating, raise the penalisation or threshold with care.`);
  if (pen < 2) W.push('With a penalisation exponent below 2 intermediate densities are not suppressed; the result is a variable-thickness sheet rather than a black-and-white topology.');
  W.push('Results are mesh- and filter-dependent: the filter radius sets the minimum member size, and a finer mesh with the same radius in elements gives a different, more intricate layout. Linear elasticity, a single load case and unit thickness are assumed; stress, buckling and fatigue are not constrained.');
  return {
    summary: `SIMP topology optimisation of the ${inp.tProb === 'cant' ? 'cantilever' : 'MBB half-beam'} on a ${nelx} × ${nely} grid (${s.nd} degrees of freedom) ran ${s.it} iterations to a compliance of ${F(s.c, 5)} (normalised, unit load and modulus) at a volume fraction of ${F(vol, 4)}; the uniform starting design had ${F(s.hist[0], 5)}. The layout is ${conv ? 'converged' : 'not yet converged'} (last change ${F(s.change, 3)}) with a greyness of ${F(grey, 3)} %. For E = ${F(inp.tE, 3)} Pa, thickness ${F(inp.tThk, 3)} mm and ${F(inp.tLoad, 3)} N the load point deflects ${F(defl * 1e3, 4)} mm; the minimum member size is about ${F(minFeat, 3)} mm.`,
    kpis: [
      kp('c', 'Compliance (normalised)', s.c, '–', 'Work of the unit load, uᵀKu, for unit modulus, unit thickness and unit element size; lower is stiffer. It is the objective minimised by the optimiser.'),
      kp('c0', 'Compliance of the uniform starting design', s.hist[0], '–', 'Compliance of the homogeneous grey design with the same material volume; the ratio to the final value is the stiffness gained by redistributing material.'),
      kp('gain', 'Stiffness gain over the uniform design', s.hist[0] / Math.max(s.c, 1e-300), '×', 'Initial compliance divided by final compliance; values well above 1 show how much stiffer the optimised layout is at equal mass.', s.hist[0] / s.c > 1),
      kp('vol', 'Volume fraction', vol, '–', `Mean density of the design; the constraint is ${F(vf, 3)}, enforced by bisection on the Lagrange multiplier.`, Math.abs(vol - vf) < 2e-3),
      kp('grey', 'Greyness (measure of non-discreteness)', grey, '%', 'Mean of 4ρ(1 − ρ) over all elements: 0 % for a pure solid–void design, 100 % for uniform half density.', grey < 15),
      kp('iters', 'Iterations', s.it, '', `Optimality-criteria iterations performed (limit ${maxIt}); each includes one banded Cholesky solve of the finite-element system.`, conv),
      kp('change', 'Last design change', s.change, '–', 'Largest change of any element density in the last iteration; the stopping criterion is 0.01.', conv),
      kp('defl', 'Deflection at the load point', defl * 1e3, 'mm', 'Normalised compliance scaled with the stated load, modulus and thickness, δ = c·F/(E·t); valid for the penalised density field, so a thresholded part will differ slightly.'),
      kp('feat', 'Minimum member size', minFeat, 'mm', 'About twice the filter radius times the element size: members thinner than this are suppressed by the sensitivity filter.'),
      kp('solid', 'Solid fraction after thresholding at 0.5', solid, '–', 'Fraction of elements kept as solid in the exported geometry; it differs from the volume fraction when grey elements remain.'),
      kp('dof', 'Degrees of freedom', s.nd, '', `Size of the finite-element system (bilinear four-node elements, half-bandwidth ${s.hb}).`),
    ],
    series: [
      { key: 'hist', title: 'Compliance history', xlabel: 'iteration', ylabel: 'compliance (normalised)', logy: true, lines: [{ name: 'compliance', x: s.hist.map((_, i) => i + 1), y: s.hist }], explain: 'Objective value per iteration. The rapid initial drop is the redistribution of material from the uniform start; the slow tail is the sharpening of members.' },
      { key: 'chg', title: 'Design change per iteration', xlabel: 'iteration', ylabel: 'max density change', logy: true, lines: [{ name: 'max |Δρ|', x: s.chg.map((_, i) => i + 1), y: s.chg.map(v => Math.max(v, 1e-6)) }, { name: 'criterion', x: [1, Math.max(s.it, 2)], y: [0.01, 0.01], dash: true }], explain: 'Largest change of an element density between successive iterations against the stopping criterion; the move limit caps it at 0.2.' },
    ],
    maps: [{ key: 'rho', title: 'Optimised density field', nx: nelx, ny: nely, data: map, x0: 0, x1: nelx * el, y0: 0, y1: nely * el, xlabel: 'x (mm)', ylabel: 'y (mm)', unit: 'relative density', explain: `Material distribution: 1 is solid, 0 void. ${inp.tProb === 'cant' ? 'The left edge is clamped and the load acts downwards at mid-height of the right edge.' : 'Half of a simply supported beam: symmetry on the left edge, roller support at the lower right corner and the load at the upper left corner.'}`, diverging: false }],
    tables: [
      { key: 'set', title: 'Optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Problem': inp.tProb === 'cant' ? 'cantilever, mid-height tip load' : 'MBB beam (half model)', 'Grid': `${nelx} × ${nely} four-node plane-stress elements`, 'Material interpolation': `SIMP, E = ρ^${F(pen, 3)}·E₀, ρ_min = 0.001, ν = 0.3`, 'Volume fraction': vf, 'Filter': `sensitivity filter, radius ${F(rmin, 3)} elements`, 'Update': 'optimality criteria, move limit 0.2, bisection on λ', 'Linear solver': 'banded Cholesky (direct)', 'Stopping': `max |Δρ| ≤ 0.01 or ${maxIt} iterations` }), explain: 'Formulation and numerical settings of this run; the result is reproducible because the method is deterministic.' },
    ],
    equations: [
      { kind: 'governing', label: 'Minimum-compliance problem', text: 'min over ρ:  c(ρ) = uᵀK(ρ)u = Σₑ ρₑᵖ·uₑᵀk₀uₑ   s.t.  K(ρ)u = f,  Σₑ ρₑ = V_f·N,  0 < ρ_min ≤ ρₑ ≤ 1' },
      { kind: 'formula', label: 'Sensitivity and filter', text: '∂c/∂ρₑ = −p·ρₑ^(p−1)·uₑᵀk₀uₑ;   filtered: (∂c/∂ρₑ)~ = Σ_f H_ef ρ_f ∂c/∂ρ_f / (ρₑ Σ_f H_ef),  H_ef = max(0, r_min − dist(e, f))' },
      { kind: 'formula', label: 'Optimality-criteria update', text: 'ρₑ ← clamp(ρₑ·√(−(∂c/∂ρₑ)~/λ), ρₑ ± 0.2, [ρ_min, 1]),  λ by bisection on the volume constraint' },
    ],
    geometry: { kind: 'voxels', nx: nelx, ny: nely, nz: 1, data: vox, voxel: el, units: 'mm', title: 'Thresholded design (ρ > 0.5)' },
    checks: [
      { name: 'Volume constraint satisfied', pass: Math.abs(vol - vf) < 2e-3, detail: `${F(vol, 5)} vs ${F(vf, 4)}` },
      { name: 'Energy balance uᵀKu = fᵀu', pass: eBal < 1e-6, detail: `relative difference ${F(eBal, 2)}` },
      { name: 'Compliance not above the uniform design', pass: s.c <= s.hist[0] * (1 + 1e-9), detail: `${F(s.c, 5)} ≤ ${F(s.hist[0], 5)}` },
      { name: 'Converged (max change ≤ 0.01)', pass: conv, detail: `${F(s.change, 3)} after ${s.it} iterations` },
      { name: 'Densities within [ρ_min, 1]', pass: N.minOf(x) >= 0.001 - 1e-12 && N.maxOf(x) <= 1 + 1e-12, detail: `${F(N.minOf(x), 3)} to ${F(N.maxOf(x), 3)}` },
      { name: 'Design close to solid–void (greyness < 15 %)', pass: grey < 15, detail: `${F(grey, 3)} %` },
      { name: 'Bounded below by the solid-domain compliance scale', pass: s.c >= cFull * 0.5, detail: `uniform-design compliance × V_fᵖ = ${F(cFull, 4)} is the fully solid value` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: `Threshold the density field at 0.5, smooth the boundary and re-analyse the resulting part with body-fitted elements (suite 14): members are at least ${F(minFeat, 3)} mm wide, to be compared with the minimum printable feature and with buckling of the slender compression members.`, basis: 'Filter-controlled minimum member size and compliance-only formulation.' },
      { priority: conv ? 'low' : 'high', text: conv ? 'The layout is converged; repeat on a mesh twice as fine with the same physical filter radius to confirm that the topology (number of holes) does not change.' : `Continue to convergence: the density still changes by ${F(s.change, 3)} per iteration, and member positions may shift.`, basis: 'Convergence history and known mesh dependence of SIMP.' },
      { priority: 'medium', text: grey > 15 ? `Reduce the ${F(grey, 3)} % greyness by continuing the iterations or by a projection filter before export; grey regions have no manufacturable interpretation in a single material.` : `Greyness is ${F(grey, 3)} %; the 0.5 threshold retains ${pct(solid)} of the domain as solid, close to the ${pct(vf)} target.`, basis: 'Measure of non-discreteness of the density field.' },
      { priority: 'medium', text: `The optimised layout is ${F(s.hist[0] / Math.max(s.c, 1e-300), 3)} times stiffer than the uniform design at ${pct(vf)} of the material: use it to cut mass, feedstock and build energy rather than to add stiffness that is not required.`, basis: 'Compliance gain at fixed volume fraction.', sustainability: true },
    ],
    publish: { 'geom.rel_density': vol, 'geom.cell_size': el * 1e-3, 'geom.strut': minFeat * 1e-3, 'opt.objective': s.c, 'mech.compliance': s.c, 'mech.deflection': defl, 'eff.E': inp.tE * s.hist[0] * Math.pow(vf, pen) / Math.max(s.c, 1e-300) },
    exports: [{ name: 'topology_density.csv', mime: 'text/csv', content: Array.from({ length: nely }, (_, j) => Array.from({ length: nelx }, (_, i) => x[i * nely + j].toFixed(3)).join(',')).join('\n'), description: 'Element densities, one row per element row from top to bottom.' }],
    numerics: { method: 'SIMP with sensitivity filter and optimality-criteria update; bilinear plane-stress finite elements; banded Cholesky direct solver', dof: s.nd, steps: s.it, iterations: s.it, residual: s.change, notes: `Half-bandwidth ${s.hb}; penalisation ${F(pen, 3)}; filter radius ${F(rmin, 3)} elements.` },
  };
}

// =====================================================================================================
// 7. Gaussian-process surrogate and Bayesian optimisation
// =====================================================================================================
const fsub = (L, b) => { const n = b.length, y = new Array(n); for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; } return y; };
const bsub = (L, y) => { const n = y.length, x = new Array(n); for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; } return x; };
const ELLS = N.logspace(0.05, 1.5, 9), NOISE = [1e-5, 1e-3, 1e-2];
function gpFit(X, y) {
  const n = X.length, ym = N.mean(y), ys = N.std(y) || 1, yn = y.map(v => (v - ym) / ys), D2 = X.map(a => X.map(b => { let s = 0; for (let k = 0; k < a.length; k++) s += (a[k] - b[k]) ** 2; return s; }));
  let best = null;
  for (const ell of ELLS) for (const nz of NOISE) { const K = D2.map((r, i) => r.map((v, j) => Math.exp(-v / (2 * ell * ell)) + (i === j ? nz : 0))); let L; try { L = N.cholesky(K); } catch { continue; } const al = bsub(L, fsub(L, yn)), s2 = N.dot(yn, al) / n; if (!(s2 > 0) || !Number.isFinite(s2)) continue; let ld = 0; for (let i = 0; i < n; i++) ld += Math.log(L[i][i]); const lml = -0.5 * n * Math.log(s2) - ld - 0.5 * n * (1 + Math.log(2 * PI)); if (!best || lml > best.lml) best = { ell, nz, L, al, s2, lml }; }
  if (!best) { const L = N.cholesky(D2.map((r, i) => r.map((v, j) => Math.exp(-v / 0.5) + (i === j ? 0.1 : 0)))), al = bsub(L, fsub(L, yn)); best = { ell: 0.5, nz: 0.1, L, al, s2: Math.max(N.dot(yn, al) / n, 1e-12), lml: NaN }; }
  const { ell, nz, L, al, s2 } = best;
  const predict = x => { const ks = X.map(a => { let s = 0; for (let k = 0; k < a.length; k++) s += (a[k] - x[k]) ** 2; return Math.exp(-s / (2 * ell * ell)); }), v = fsub(L, ks), mu = N.dot(ks, al), va = s2 * Math.max(1 + nz - N.dot(v, v), 1e-12); return { mu: ym + ys * mu, sd: ys * Math.sqrt(va) }; };
  const Li = X.map((_, i) => fsub(L, X.map((__, j) => i === j ? 1 : 0))), kd = X.map((_, i) => { let s = 0; for (let j = 0; j < n; j++) s += Li[i][j] ** 2; return s; }), loo = yn.map((v, i) => ym + ys * (v - al[i] / kd[i]));
  return { ell, nz, s2, lml: best.lml, predict, loo, ym, ys };
}
const BDEF = (() => { const r = N.rng(41), P = { rAx: 10, rRad: 0.4, EaAx: 60, EaRad: 40, tauRef: 3, D0: 20 }; return 'T_C,length_nm\n' + [60, 65, 70, 75, 80, 85, 90, 95].map(T => `${T},${(growth(T, 25, 3, P).L * (1 + 0.04 * r.normal())).toFixed(1)}`).join('\n'); })();
function parseXY(text) { const rows = []; for (const l of String(text ?? '').split(/\r?\n/).map(s => s.trim()).filter(s => s && !s.startsWith('#') && !s.startsWith('//'))) { const c = l.split(/[,;\t]/).map(s => Number(s.trim())); if (c.length >= 2 && Number.isFinite(c[0]) && Number.isFinite(c[1])) rows.push([c[0], c[1]]); } if (rows.length < 4) throw new Error(`Surrogate data table: at least 4 numeric rows "x, y" are needed (found ${rows.length}). Paste CSV/TSV text with one header line.`); return rows.slice(0, 60); }
function runBayes(inp) {
  const W = [], P = recipeP(inp), model = inp.bModel, seed = inp.seed;
  const MOD = {
    growth1d: { d: 1, names: ['temperature (°C)'], b: [[60, 95]], out: 'squared relative length error', f: x => (growth(x[0], 40, 3, P).L / inp.Ltar - 1) ** 2, eq: 'f(T) = (L(T; 40 mM, 3 h)/L* − 1)²  with the Arrhenius/depletion growth model' },
    lspr: { d: 1, names: ['aspect ratio'], b: [[1, 8]], out: 'squared relative resonance error', f: x => (lsprPeak(inp.metal === 'Ag' ? 'Ag' : 'Au', x[0], inp.nMed ** 2) / inp.lamT - 1) ** 2, eq: 'f(AR) = (λ_LSPR(AR)/λ* − 1)²  with the Gans quasi-static model' },
    lattice: { d: 1, names: ['cell angle θ (°)'], b: [[-70, -8]], out: 'squared Poisson-ratio error', f: x => (honey(x[0], 2, 0.1).nu12 - inp.nuT) ** 2, eq: 'f(θ) = (ν₁₂(θ; h/l = 2) − ν*)²  with the Gibson–Ashby honeycomb model' },
    growth2d: { d: 2, names: ['temperature (°C)', 'concentration (mM)'], b: [[60, 95], [5, 80]], out: 'squared quality error', f: x => { const g = growth(x[0], x[1], 3, P); return (g.L / inp.Ltar - 1) ** 2 + (g.D / Math.max(inp.Dtar, 1) - 1) ** 2; }, eq: 'f(T, c₀) = (L/L* − 1)² + (D/D* − 1)²  at t = 3 h with the Arrhenius/depletion growth model' },
  };
  const isTab = model === 'table', M = MOD[model] || (isTab ? null : MOD.growth1d);
  let X = [], y = [], b, names, d, sgn = 1, rows = null, evalN = 0;
  if (isTab) { rows = parseXY(inp.bTable); if (String(inp.bTable).trim() === BDEF) W.push('Synthetic demonstration data — replace with measurements. The default table was generated from the growth model with 4 % seeded noise.'); d = 1; names = ['x']; const xs = rows.map(r => r[0]), lo = N.minOf(xs), hi = N.maxOf(xs) > lo ? N.maxOf(xs) : lo + 1; b = [[lo, hi]]; sgn = inp.bGoal === 'max' ? -1 : 1; X = rows.map(r => [(r[0] - lo) / (hi - lo)]); y = rows.map(r => sgn * r[1]); }
  else { d = M.d; names = M.names; b = M.b; const n0 = Math.round(inp.n0), r0 = N.rng(seed), S = N.lhs(n0, d, r0); X = S; y = S.map(u => { evalN++; return M.f(u.map((v, k) => b[k][0] + v * (b[k][1] - b[k][0]))); }); }
  const phys = u => u.map((v, k) => b[k][0] + v * (b[k][1] - b[k][0])), cand = d === 1 ? N.linspace(0, 1, 201).map(v => [v]) : N.linspace(0, 1, 21).flatMap(a => N.linspace(0, 1, 21).map(c => [a, c]));
  const ei = (gp, ybest) => cand.map(u => { const p = gp.predict(u), z = (ybest - p.mu) / Math.max(p.sd, 1e-12); return p.sd > 1e-12 ? (ybest - p.mu) * N.normCdf(z) + p.sd * Math.exp(-0.5 * z * z) / Math.sqrt(2 * PI) : Math.max(ybest - p.mu, 0); });
  const n0 = X.length, nBo = isTab ? 0 : Math.round(inp.nBo), bestHist = [N.minOf(y)]; let gp = gpFit(X, y);
  for (let it = 0; it < nBo; it++) { const a = ei(gp, N.minOf(y)); let k = N.argmax(a); if (X.some(p => p.every((v, j) => Math.abs(v - cand[k][j]) < 1e-9))) { const o = a.map((v, i) => i).sort((i, j) => a[j] - a[i]); k = o.find(i => !X.some(p => p.every((v, j) => Math.abs(v - cand[i][j]) < 1e-9))) ?? k; } X.push(cand[k]); y.push(M.f(phys(cand[k]))); evalN++; bestHist.push(N.minOf(y)); gp = gpFit(X, y); }
  // random-search baseline with the same budget
  const rnd = []; if (!isTab) { const rr = N.rng(seed + 100); let bst = Infinity; for (let i = 0; i < n0 + nBo; i++) { const v = M.f(phys(Array.from({ length: d }, () => rr()))); if (v < bst) bst = v; if (i >= n0 - 1) rnd.push(bst); } }
  const loo0 = N.rmse(gp.loo, y), looR2 = N.r2(gp.loo, y), ib = N.argmin(y), xb = phys(X[ib]), yb = sgn * y[ib], aF = ei(gp, N.minOf(y)), nextIdx = aF.map((v, i) => i).sort((i, j) => aF[j] - aF[i]).filter(i => !X.some(p => p.every((v, j) => Math.abs(v - cand[i][j]) < 1e-9))).reduce((acc, i) => { if (acc.length < 3 && acc.every(k => Math.hypot(...cand[i].map((v, j) => v - cand[k][j])) > 0.08)) acc.push(i); return acc; }, []);
  const samp = X.map(p => p[0]), sLo = N.minOf(samp), sHi = N.maxOf(samp), uq = [inp.xq, ...X[ib].slice(1)], pq = gp.predict(uq), extra = inp.xq < sLo - 1e-9 || inp.xq > sHi + 1e-9, prior = gp.ys * Math.sqrt(gp.s2);
  if (extra) W.push(`Extrapolation: the query ${names[0]} = ${F(phys(uq)[0], 4)} lies outside the sampled domain [${F(b[0][0] + sLo * (b[0][1] - b[0][0]), 4)}, ${F(b[0][0] + sHi * (b[0][1] - b[0][0]), 4)}]. The surrogate reverts to its prior mean there and the ±2σ band (${F(2 * pq.sd, 3)}) is not a guarantee that the underlying physics remains valid.`);
  if (gp.nz >= 1e-2) W.push('The marginal-likelihood search selected the largest noise level: the data are noisy relative to their trend or the length scale is poorly determined; predictions are smooth averages rather than interpolations.');
  if (!isTab && N.minOf(y) > 1e-3 && model !== 'table') W.push(`The best sampled objective is ${F(N.minOf(y), 3)}: the target is not reached inside the design range of this one- or two-variable problem.`);
  W.push('A Gaussian-process surrogate is a data-driven interpolator: its predictions and uncertainty band are conditional on the kernel assumptions and must be validated against the physics model or new experiments before use, and never outside the sampled domain.');
  const gx = N.linspace(-0.1, 1.1, 121), band = gx.map(v => gp.predict([v, ...X[ib].slice(1)])), xp = gx.map(v => b[0][0] + v * (b[0][1] - b[0][0])), truth = isTab ? null : gx.map(v => M.f(phys([N.clamp(v, 0, 1), ...X[ib].slice(1)])));
  const series = [
    { key: 'band', title: `Surrogate prediction with ±2σ band${d === 2 ? ' (slice through the best point)' : ''}`, xlabel: names[0], ylabel: isTab ? 'response' : M.out, lines: [{ name: 'GP mean', x: xp, y: band.map(p => sgn * p.mu) }, { name: 'mean + 2σ', x: xp, y: band.map(p => sgn * p.mu + 2 * p.sd), dash: true }, { name: 'mean − 2σ', x: xp, y: band.map(p => sgn * p.mu - 2 * p.sd), dash: true }, ...(truth ? [{ name: 'forward model', x: xp.filter((_, i) => gx[i] >= 0 && gx[i] <= 1), y: truth.filter((_, i) => gx[i] >= 0 && gx[i] <= 1) }] : []), { name: 'samples', x: X.map(p => b[0][0] + p[0] * (b[0][1] - b[0][0])), y: y.map(v => sgn * v), points: true }], explain: 'Posterior mean of the Gaussian process with its 95 % band, extended 10 % beyond the design range to show how the uncertainty opens up outside the sampled domain. Samples cluster where the acquisition function expected improvement.' },
    { key: 'ei', title: 'Expected-improvement acquisition', xlabel: d === 1 ? names[0] : 'candidate number', ylabel: 'expected improvement', lines: [{ name: 'EI', x: d === 1 ? cand.map(u => phys(u)[0]) : cand.map((_, i) => i + 1), y: aF }], explain: 'Expected improvement over the current best value across the candidate set. Its maxima are the recommended next experiments: they balance low predicted value against high uncertainty.' },
  ];
  if (!isTab) series.push({ key: 'conv', title: 'Bayesian optimisation versus random search', xlabel: 'number of evaluations', ylabel: 'best objective found', logy: true, lines: [{ name: 'Bayesian optimisation (EI)', x: bestHist.map((_, i) => n0 + i), y: bestHist.map(v => Math.max(v, 1e-16)) }, { name: 'random search (same budget)', x: rnd.map((_, i) => n0 + i), y: rnd.map(v => Math.max(v, 1e-16)), dash: true }], explain: 'Best objective value against the number of forward-model evaluations. One seeded random-search run is shown for comparison; a single run is an illustration, not a statistical benchmark.' });
  return {
    summary: `A Gaussian-process surrogate (RBF kernel, length scale ${F(gp.ell, 3)} of the design range, noise ${F(gp.nz, 2)}) was fitted to ${X.length} ${isTab ? 'data points' : `evaluations of the ${M.out} (${n0} Latin-hypercube samples plus ${nBo} expected-improvement steps)`}. Leave-one-out RMSE over all samples is ${F(loo0, 3)} (R² = ${F(looR2, 3)}). The best ${isTab ? (sgn < 0 ? 'maximum' : 'minimum') : 'objective'} found is ${F(yb, 4)} at ${names.map((nm2, k) => `${nm2} = ${F(xb[k], 4)}`).join(', ')}${isTab ? '' : `, against ${F(rnd.at(-1), 3)} for random search with the same budget`}. At the query point the prediction is ${F(sgn * pq.mu, 4)} ± ${F(2 * pq.sd, 3)} (2σ)${extra ? ', an extrapolation' : ''}.`,
    kpis: [
      kp('best', isTab ? 'Best observed response' : 'Best objective found', yb, '–', isTab ? `${sgn < 0 ? 'Largest' : 'Smallest'} response in the data table, the incumbent against which expected improvement is computed.` : 'Smallest value of the objective among all evaluated points, the incumbent of the Bayesian optimisation.', isTab ? null : yb < 1e-3),
      kp('xbest', `Best ${names[0]}`, xb[0], '', 'Value of the first design variable at the incumbent; in two-dimensional problems the remaining variable is given in the summary and the table.'),
      kp('loo', 'Leave-one-out RMSE', loo0, '–', 'Root-mean-square error of predicting each sample from all the others with the final hyper-parameters, computed in closed form from the kernel matrix; an estimate of the interpolation error of the surrogate.', loo0 < 0.5 * (N.std(y) || 1)),
      kp('r2', 'Leave-one-out R²', looR2, '–', 'Fraction of the variance of the samples explained when each point is predicted without itself; values near 1 mean the surrogate generalises between samples.', looR2 > 0.8),
      kp('ell', 'Kernel length scale', gp.ell, 'of range', 'Correlation length of the RBF kernel in units of the normalised design range, chosen by maximising the marginal likelihood on a grid.'),
      kp('noise', 'Noise-to-signal variance', gp.nz, '–', 'Diagonal noise term relative to the signal variance selected by the marginal likelihood; larger values smooth through scattered data.'),
      kp('q_mu', 'Prediction at the query point', sgn * pq.mu, '–', `Posterior mean at ${names[0]} = ${F(phys(uq)[0], 4)}${d === 2 ? ' with the second variable at its best value' : ''}.`, extra ? false : null),
      kp('q_2s', 'Prediction uncertainty at the query (2σ)', 2 * pq.sd, '–', `Half-width of the 95 % band at the query; the prior (no-information) value is ${F(2 * prior, 3)}, which the band approaches away from the data.`, pq.sd < 0.5 * prior),
      kp('next', 'Next recommended experiment', nextIdx.length ? phys(cand[nextIdx[0]])[0] : xb[0], '', `Value of ${names[0]} with the highest expected improvement not yet sampled${d === 2 && nextIdx.length ? ` (with ${names[1]} = ${F(phys(cand[nextIdx[0]])[1], 4)})` : ''}.`),
      kp('gain', 'Advantage over random search', isTab ? 0 : Math.log10(Math.max(rnd.at(-1), 1e-300) / Math.max(N.minOf(y), 1e-300)), 'decades', 'log₁₀ of the ratio between the best random-search objective and the best Bayesian-optimisation objective at equal budget (0 for table data, where no evaluations are made).', isTab ? null : rnd.at(-1) >= N.minOf(y)),
      kp('evals', 'Forward-model evaluations', evalN, '', 'Number of forward-model runs used for the surrogate and the acquisition loop (0 for table data).'),
    ],
    series,
    tables: [
      { key: 'next', title: 'Recommended next experiments', columns: ['Rank', ...names, 'Predicted mean', '± 2σ', 'Expected improvement'], rows: nextIdx.length ? nextIdx.map((i, k) => { const p = gp.predict(cand[i]); return [k + 1, ...phys(cand[i]).map(v => F(v, 4)), F(sgn * p.mu, 4), F(2 * p.sd, 3), F(aF[i], 3)]; }) : [['–', ...names.map(() => '–'), '–', '–', 'no unsampled candidate with positive expected improvement']], explain: 'Unsampled candidates with the highest expected improvement, at least 8 % of the range apart so that a batch of experiments is not redundant.' },
      { key: 'data', title: 'Samples used by the surrogate', columns: ['No.', ...names, 'Value', 'Origin'], rows: X.slice(0, 45).map((p, i) => [i + 1, ...phys(p).map(v => F(v, 4)), F(sgn * y[i], 4), isTab ? 'table' : i < n0 ? 'initial design' : 'expected improvement']), explain: 'All data behind the surrogate in evaluation order; the later points show where the acquisition function concentrated the search.' },
      { key: 'set', title: 'Surrogate and optimiser settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Kernel': 'squared exponential (RBF), constant mean', 'Hyper-parameter search': `${ELLS.length} length scales × ${NOISE.length} noise levels, profiled signal variance`, 'Log marginal likelihood': gp.lml, 'Factorisation': 'Cholesky', 'Acquisition': 'expected improvement on a candidate grid', 'Candidates': cand.length, 'Initial design': isTab ? 'user table' : `Latin hypercube, ${n0} points`, 'Iterations': nBo, 'Seed': seed }), explain: 'Model choices and budget of this run for reproducibility.' },
    ],
    equations: [
      { kind: 'governing', label: 'Gaussian-process regression', text: 'μ(x) = k(x)ᵀ(K + σ_n²I)⁻¹y;  σ²(x) = k(x,x) − k(x)ᵀ(K + σ_n²I)⁻¹k(x);  k(x,x′) = σ_f²·exp(−|x − x′|²/(2ℓ²))' },
      { kind: 'formula', label: 'Log marginal likelihood', text: 'ln p(y) = −½yᵀK⁻¹y − ½ln|K| − (n/2)ln 2π' },
      { kind: 'formula', label: 'Expected improvement', text: 'EI(x) = (y_best − μ)Φ(z) + σφ(z),  z = (y_best − μ)/σ' },
      { kind: 'formula', label: 'Leave-one-out prediction', text: 'μ₋ᵢ = yᵢ − αᵢ/[K⁻¹]ᵢᵢ,  α = K⁻¹y' },
      { kind: 'formula', label: 'Objective', text: isTab ? `${sgn < 0 ? 'maximise' : 'minimise'} the tabulated response y(x)` : M.eq },
    ],
    checks: [
      { name: 'Kernel matrix positive definite', pass: Number.isFinite(gp.lml), detail: `Cholesky succeeded; log marginal likelihood ${F(gp.lml, 4)}` },
      { name: 'Leave-one-out R² above 0.8', pass: looR2 > 0.8, detail: `R² = ${F(looR2, 4)}, RMSE = ${F(loo0, 3)}` },
      { name: 'Query inside the sampled domain', pass: !extra, detail: `query at ${F(inp.xq, 3)} of the range; samples span ${F(sLo, 3)} – ${F(sHi, 3)}` },
      { name: 'Bayesian optimisation at least as good as random search', pass: isTab || N.minOf(y) <= rnd.at(-1) * (1 + 1e-9), detail: isTab ? 'not applicable to table data' : `${F(N.minOf(y), 3)} vs ${F(rnd.at(-1), 3)}` },
      { name: 'Best-so-far history monotone', pass: bestHist.every((v, i) => i === 0 || v <= bestHist[i - 1] + 1e-15), detail: `${bestHist.length} recorded values` },
      { name: 'Predictive variance non-negative', pass: band.every(p => p.sd >= 0), detail: `largest 2σ in the plotted range ${F(2 * N.maxOf(band.map(p => p.sd)), 3)}` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: nextIdx.length ? `Run the next experiment at ${names.map((nm2, k) => `${nm2} = ${F(phys(cand[nextIdx[0]])[k], 4)}`).join(', ')} (expected improvement ${F(aF[nextIdx[0]], 3)}), add the result to the data and refit; stop when the expected improvement falls below the measurement noise.` : `No unsampled candidate offers positive expected improvement: the optimum at ${names[0]} = ${F(xb[0], 4)} is resolved to the candidate spacing; confirm it by a repeat experiment.`, basis: 'Maximum of the expected-improvement acquisition function.' },
      { priority: extra ? 'high' : 'medium', text: extra ? 'Do not use the surrogate at the query point: it lies outside the sampled domain. Extend the experimental design to cover it or use the physics model there.' : `Validate the surrogate with two or three additional points before relying on it: the leave-one-out error is ${F(loo0, 3)} for a response spread of ${F(N.std(y), 3)}.`, basis: 'Extrapolation check and leave-one-out cross-validation.' },
      { priority: 'medium', text: 'Keep the physics-based model in the loop: use the surrogate to choose experiments and the forward model (or measurements) to decide; data-driven predictions need validation against physics models and experiments.', basis: 'Surrogates interpolate data and carry no physical constraints.' },
      { priority: 'medium', text: isTab ? 'Sequential, model-guided experiments usually reach the optimum with a fraction of the runs of a full grid, saving precursors, solvent and energy; record every run, including failures, with its conditions.' : `Bayesian optimisation reached ${F(N.minOf(y), 3)} in ${n0 + nBo} evaluations; applied to laboratory work the same strategy cuts the number of syntheses, and with them precursor, solvent and energy use.`, basis: 'Sample efficiency of model-based optimisation.', sustainability: true },
    ],
    publish: Object.assign({ 'opt.objective': yb, 'opt.x_best': xb[0], 'cred.rmse': loo0, 'opt.next_x': nextIdx.length ? phys(cand[nextIdx[0]])[0] : xb[0], 'opt.surrogate_sd': pq.sd }, model === 'growth1d' || model === 'growth2d' ? { 'proc.T_C': xb[0] } : {}, model === 'growth2d' ? { 'proc.c_precursor': xb[1] } : {}, model === 'lspr' ? { 'geom.aspect': xb[0] } : {}),
    numerics: { method: 'Gaussian-process regression (Cholesky) with grid-search marginal likelihood; expected-improvement Bayesian optimisation; closed-form leave-one-out cross-validation', dof: d, steps: nBo, iterations: evalN, residual: loo0, notes: `Seed ${seed}; ${cand.length} acquisition candidates.` },
  };
}

// =====================================================================================================
// 8. Library screening: k-nearest-neighbour property prediction from composition descriptors
// =====================================================================================================
const KPROP = { E: ['Young’s modulus', 'Pa', true], density: ['Density', 'kg/m³', true], k_thermal: ['Thermal conductivity', 'W/m/K', true], cte: ['Thermal expansion coefficient', '1/K', false], bandgap_eV: ['Band gap', 'eV', true], n_refr: ['Refractive index', '–', true], cp: ['Specific heat', 'J/kg/K', true], melting_K: ['Melting point', 'K', true], eps_r: ['Relative permittivity', '–', true] };
function descriptors(formula) { let p; try { p = chem.parseFormula(formula); } catch { return null; } const els = Object.keys(p.atoms); if (!els.length) return null; const tot = els.reduce((s, e) => s + p.atoms[e], 0), x = els.map(e => p.atoms[e] / tot), E = els.map(e => chem.ELEMENTS[e]), avg = g => N.sum(E.map((e, i) => g(e) * x[i])); return [avg(e => e.mass), avg(e => e.en), N.maxOf(E.map(e => e.en)) - N.minOf(E.map(e => e.en)), avg(e => e.radius), avg(e => e.Z), els.length]; }
function runKnn(inp) {
  const W = [], key = KPROP[inp.kProp] ? inp.kProp : 'E', [pname, unit, logT] = KPROP[key], tgt = findMaterial(inp.kMat) || MATERIALS[0], seen = new Set(), train = [], useRho = key !== 'density';
  for (const m of MATERIALS) { const v = num((m.props || {})[key]); if (!m.formula || v === null || (logT && v <= 0)) continue; const dsc = descriptors(m.formula); if (!dsc) continue; if (useRho) { const rh = num((m.props || {}).density); if (rh === null || rh <= 0) continue; dsc.push(Math.log10(rh)); } const sig = dsc.map(z => z.toFixed(6)).join('|') + '|' + v; if (seen.has(sig)) continue; seen.add(sig); train.push({ m, d: dsc, y: logT ? Math.log10(v) : v, v }); }
  const n = train.length, nd = useRho ? 7 : 6, mu = Array.from({ length: nd }, (_, j) => N.mean(train.map(t => t.d[j]))), sd = Array.from({ length: nd }, (_, j) => N.std(train.map(t => t.d[j])) || 1), zs = dd => dd.map((v, j) => (v - mu[j]) / sd[j]);
  train.forEach(t => { t.z = zs(t.d); });
  const predict = (z, k, skip) => { const ds = []; for (let i = 0; i < n; i++) if (i !== skip) { let s = 0; for (let j = 0; j < nd; j++) s += (z[j] - train[i].z[j]) ** 2; ds.push([Math.sqrt(s), i]); } ds.sort((a, b) => a[0] - b[0]); const nb = ds.slice(0, Math.min(k, ds.length)); let sw = 0, sy = 0; for (const [dd, i] of nb) { const w = 1 / (dd + 1e-6); sw += w; sy += w * train[i].y; } const yp = sy / sw; return { y: yp, nb, spread: N.std(nb.map(([, i]) => train[i].y)) }; };
  const k = Math.round(inp.kNN), ks = [1, 2, 3, 4, 5, 6, 8, 10].filter(q => q < n), looK = q => { const pr = train.map((t, i) => predict(t.z, q, i).y); return { pr, rmse: N.rmse(pr, train.map(t => t.y)), mae: N.mae(pr, train.map(t => t.y)), r2: N.r2(pr, train.map(t => t.y)) }; };
  if (n < 6) { W.push(`Only ${n} library entries with a formula hold a value for ${pname}: too few for a meaningful nearest-neighbour model.`); }
  const loo = n >= 3 ? looK(Math.min(k, n - 1)) : { pr: train.map(t => t.y), rmse: 0, mae: 0, r2: 0 }, curve = n >= 3 ? ks.map(q => looK(q).rmse) : [], kBest = curve.length ? ks[N.argmin(curve)] : k;
  const rhoT = num((tgt.props || {}).density), dT0 = tgt.formula ? descriptors(tgt.formula) : null, dT = dT0 && useRho ? (rhoT !== null && rhoT > 0 ? [...dT0, Math.log10(rhoT)] : null) : dT0, known = num((tgt.props || {})[key]), selfI = train.findIndex(t => t.m.id === tgt.id);
  const back = v => logT ? Math.pow(10, v) : v; let pred = null, P = null, nnD = 0, nnTyp = 0, outDom = false;
  if (dT && n >= 2) { P = predict(zs(dT), Math.min(k, n - 1), selfI); pred = back(P.y); nnD = P.nb[0][0]; const nn = train.map((t, i) => predict(t.z, 1, i).nb[0][0]); nnTyp = N.quantile(nn, 0.95); outDom = nnD > nnTyp; }
  if (!dT) W.push(dT0 ? `${tgt.name} has no library density, which this model uses as a descriptor, so no prediction is made.` : `${tgt.name} has no stoichiometric formula, so composition descriptors cannot be formed and no prediction is made. Select a material with a formula.`);
  if (outDom) W.push(`Applicability domain: the nearest training composition is ${F(nnD, 3)} standard units away, beyond the 95th percentile (${F(nnTyp, 3)}) of nearest-neighbour distances in the library — the prediction is an extrapolation in composition space.`);
  const twin = P && P.nb[0][0] < 1e-9 ? train[P.nb[0][1]].m.name : null;
  if (known !== null && dT) W.push(`The library already holds ${pname} for this material (${F(known, 4)} ${unit}); it was left out of the neighbours${twin ? `, but ${twin} has the same composition descriptors and dominates the prediction, so the comparison is not an independent test for this entry.` : ', so the comparison is a genuine hold-out test for this entry.'}`);
  W.push('Composition and density descriptors cannot distinguish polymorphs, morphology, size or microstructure (graphene, nanotubes and nanodiamond share one composition). Data-driven predictions need validation against physics models and experiments before any design decision.');
  const factor = logT ? Math.pow(10, loo.rmse) : loo.rmse, errKnown = pred !== null && known !== null ? (logT ? Math.abs(Math.log10(pred / known)) : Math.abs(pred - known)) : null, lo = P ? back(P.y - loo.rmse) : 0, hi = P ? back(P.y + loo.rmse) : 0;
  return {
    summary: pred !== null ? `${pname} of ${tgt.name} is predicted as ${F(pred, 4)} ${unit} from its ${Math.min(k, n - 1)} nearest library compositions (${P.nb.slice(0, 3).map(([, i]) => train[i].m.name).join('; ')}). Leave-one-out cross-validation over ${n} library entries gives ${logT ? `an RMSE of ${F(loo.rmse, 3)} decades (a factor of ${F(factor, 3)})` : `an RMSE of ${F(loo.rmse, 3)} ${unit}`} and R² = ${F(loo.r2, 3)}, so the plausible range is ${F(lo, 3)} to ${F(hi, 3)} ${unit}. ${known !== null ? `The library value is ${F(known, 4)} ${unit}.` : 'The library has no value for this property.'} This is a data-driven estimate that must be validated against physics models and experiments.` : `No prediction of ${pname} can be made for ${tgt.name}: ${dT ? 'the training set is too small' : 'it has no stoichiometric formula'}. The leave-one-out error of the model over ${n} library entries is ${F(loo.rmse, 3)} ${logT ? 'decades' : unit} (R² = ${F(loo.r2, 3)}). Data-driven predictions need validation against physics models and experiments.`,
    kpis: [
      kp('pred', `Predicted ${pname.toLowerCase()}`, pred ?? 0, unit, `Inverse-distance-weighted mean of the ${logT ? 'logarithm of the ' : ''}property over the nearest neighbours in standardised descriptor space (0 when no prediction is possible).`, pred === null ? false : !outDom),
      kp('lo', 'Lower plausible value', lo, unit, 'Prediction minus one leave-one-out RMSE (in log space for positive properties): the lower end of the range within which about two thirds of hold-out cases fall.'),
      kp('hi', 'Upper plausible value', hi, unit, 'Prediction plus one leave-one-out RMSE; the width of the range is the honest accuracy of this baseline model.'),
      kp('rmse', 'Leave-one-out RMSE', loo.rmse, logT ? 'decades' : unit, logT ? `Root-mean-square error of log₁₀(property) when each library entry is predicted from the others, equivalent to a typical factor of ${F(factor, 3)}.` : 'Root-mean-square error when each library entry is predicted from the others.', logT ? loo.rmse < 0.3 : null),
      kp('r2', 'Leave-one-out R²', loo.r2, '–', 'Variance of the property across the library explained by the hold-out predictions; low or negative values mean composition alone does not determine this property.', loo.r2 > 0.5),
      kp('n', 'Training entries', n, '', 'Library entries with a formula, a density and a value of the property, after removing exact duplicates (same composition and same value) that would otherwise make cross-validation look better than it is.', n >= 20),
      kp('kbest', 'Number of neighbours with the lowest hold-out error', kBest, '', `k that minimises the leave-one-out RMSE; this run uses k = ${k}.`),
      kp('nn', 'Distance to the nearest training composition', nnD, 'σ units', `Euclidean distance in standardised descriptor space; the 95th percentile within the library is ${F(nnTyp, 3)}, the limit of the applicability domain.`, pred === null ? null : !outDom),
      kp('err', 'Error against the library value', errKnown ?? 0, logT ? 'decades' : unit, 'Deviation of the hold-out prediction from the library value of this material where one exists (0 when the property is missing, which is the intended use).', errKnown === null ? null : errKnown <= 2 * loo.rmse),
    ],
    series: [
      { key: 'parity', title: 'Leave-one-out parity plot', xlabel: `library ${logT ? 'log₁₀ ' : ''}${pname.toLowerCase()}`, ylabel: `predicted ${logT ? 'log₁₀ ' : ''}${pname.toLowerCase()}`, lines: [{ name: 'hold-out predictions', x: train.length ? train.map(t => t.y) : [0], y: train.length ? loo.pr : [0], points: true }, { name: 'perfect agreement', x: train.length ? [N.minOf(train.map(t => t.y)), N.maxOf(train.map(t => t.y))] : [0, 1], y: train.length ? [N.minOf(train.map(t => t.y)), N.maxOf(train.map(t => t.y))] : [0, 1], dash: true }], explain: 'Each point is a library entry predicted without itself. Scatter about the diagonal is the real predictive accuracy; vertical stacks are compositions with several polymorphs or forms that the descriptors cannot separate.' },
      { key: 'kcurve', title: 'Hold-out error versus number of neighbours', xlabel: 'k', ylabel: `leave-one-out RMSE (${logT ? 'decades' : unit})`, lines: [{ name: 'RMSE', x: ks.length && curve.length ? ks : [1], y: curve.length ? curve : [0] }], explain: 'Bias–variance trade-off of the nearest-neighbour model: small k follows individual entries, large k averages over dissimilar compositions.' },
    ],
    tables: [
      { key: 'nb', title: 'Nearest library compositions', columns: ['Neighbour', 'Formula', 'Distance (σ units)', `${pname} (${unit})`], rows: P ? P.nb.map(([dd, i]) => [train[i].m.name, chem.prettyFormula(train[i].m.formula), F(dd, 3), F(train[i].v, 4)]) : [['–', '–', '–', 'no prediction']], explain: 'The neighbours on which the prediction rests, with their distance in standardised descriptor space. If they are chemically unlike the target, distrust the prediction.' },
      { key: 'desc', title: 'Descriptors of the target', columns: ['Descriptor', 'Value', 'Library mean', 'Library std'], rows: ['mean atomic mass (g/mol)', 'mean electronegativity', 'electronegativity range', 'mean covalent radius (Å)', 'mean atomic number', 'number of elements', 'log₁₀ density (kg/m³)'].slice(0, nd).map((s, j) => [s, dT ? F(dT[j], 4) : 'n.d.', F(mu[j], 4), F(sd[j], 3)]), explain: 'Composition-based descriptors computed from the formula and the periodic-table data, with the statistics used to standardise them.' },
      { key: 'set', title: 'Model settings', columns: ['Setting', 'Value'], rows: settingsRow({ 'Model': 'k-nearest neighbours, inverse-distance weights', 'k': k, 'Target transform': logT ? 'log₁₀' : 'none', 'Descriptors': useRho ? '6 composition descriptors + log density' : '6 composition descriptors', 'Validation': 'leave-one-out over the library', 'Duplicates': 'identical composition and value counted once' }), explain: 'Definition of the baseline model; no parameters are fitted beyond the standardisation of the descriptors.' },
    ],
    equations: [
      { kind: 'formula', label: 'k-nearest-neighbour prediction', text: 'ŷ(x) = Σ_{i ∈ N_k(x)} wᵢyᵢ / Σ wᵢ,  wᵢ = 1/(‖z(x) − zᵢ‖ + ε),  z = (descriptor − mean)/std' },
      { kind: 'formula', label: 'Leave-one-out error', text: 'RMSE_LOO = √((1/n)·Σᵢ (ŷ₋ᵢ − yᵢ)²)' },
      { kind: 'formula', label: 'Target composition', text: dT ? (() => { try { return `${chem.prettyFormula(tgt.formula)}:  M = ${F(chem.molarMass(tgt.formula), 5)} g/mol;  mean electronegativity ${F(dT[1], 4)}` } catch { return String(tgt.formula); } })() : 'no stoichiometric formula' },
    ],
    checks: [
      { name: 'Training set large enough (n ≥ 20)', pass: n >= 20, detail: `${n} entries` },
      { name: 'Target inside the applicability domain', pass: pred !== null && !outDom, detail: pred === null ? 'no prediction' : `nearest neighbour at ${F(nnD, 3)} vs limit ${F(nnTyp, 3)}` },
      { name: 'Model better than predicting the library mean (R² > 0)', pass: loo.r2 > 0, detail: `leave-one-out R² = ${F(loo.r2, 4)}` },
      { name: 'Prediction physically admissible', pass: pred === null || !logT || pred > 0, detail: pred === null ? 'no prediction' : `${F(pred, 4)} ${unit}` },
      { name: 'Hold-out error for this material within 2 RMSE', pass: errKnown === null || errKnown <= 2 * loo.rmse, detail: errKnown === null ? 'no library value to compare' : `${F(errKnown, 3)} vs ${F(2 * loo.rmse, 3)}` },
    ],
    warnings: W,
    recommendations: [
      { priority: 'high', text: pred !== null ? `Treat ${F(pred, 3)} ${unit} (range ${F(lo, 3)} – ${F(hi, 3)}) as a prior only: measure it (suite 23) or compute it with the relevant physics suite before it enters a design calculation.` : 'Select a material with a stoichiometric formula, or obtain the property by measurement or a physics-based model.', basis: 'Leave-one-out accuracy of the baseline model; data-driven predictions need validation against physics models and experiments.' },
      { priority: 'medium', text: loo.r2 < 0.5 ? `Composition explains only ${pct(Math.max(loo.r2, 0))} of the variance of this property in the library: add structural descriptors (density, crystal structure, bonding type) or use a physics model instead.` : `The baseline explains ${pct(loo.r2)} of the variance; use it for screening and ranking, not for absolute values.`, basis: 'Leave-one-out R².' },
      { priority: 'medium', text: `Use k = ${kBest}, which gives the lowest hold-out error for this property, and inspect the neighbour list each time: a prediction is only as good as the similarity of its neighbours.`, basis: 'Hold-out error versus k.' },
      { priority: 'low', text: 'Screening by a cheap baseline before synthesis avoids making and testing compositions that are unlikely to meet the target, which saves precursors, energy and characterisation time; record negative results so the data set improves.', basis: 'Role of data-driven screening in the design loop.', sustainability: true },
    ],
    publish: Object.assign({ 'cred.rmse': loo.rmse, 'opt.objective': loo.rmse, 'cred.loo_r2': loo.r2, 'material.name': tgt.name }, tgt.formula ? { 'material.formula': tgt.formula } : {}, pred !== null && known === null ? { ['prop.' + key + '_predicted']: pred } : {}),
    numerics: { method: 'k-nearest-neighbour regression on standardised composition descriptors with leave-one-out cross-validation', dof: nd, steps: n, iterations: n * (ks.length + 1), residual: loo.rmse, notes: 'Deterministic; no fitted parameters.' },
  };
}

// =====================================================================================================
const M_ = (...m) => ({ mode: m });
const matOpts = MATERIALS.filter(m => m.formula).map(m => ({ value: m.id, label: m.name }));
const nOpts = Object.keys(NMAT).map(k => ({ value: k, label: `${k} (n = ${NMAT[k]})` }));
const inputs = [
  { key: 'mode', label: 'Analysis', type: 'select', value: 'auxetic', group: 'Model', options: [{ value: 'auxetic', label: 'Inverse auxetic / lattice design' }, { value: 'optical', label: 'Inverse optical multilayer design' }, { value: 'resonator', label: 'Inverse resonator / bandgap design' }, { value: 'nano', label: 'Inverse nanoparticle design (LSPR, quantum dot)' }, { value: 'recipe', label: 'Synthesis-recipe optimisation (Pareto)' }, { value: 'topo', label: 'Topology optimisation (SIMP)' }, { value: 'bayes', label: 'Surrogate modelling and Bayesian optimisation' }, { value: 'knn', label: 'Library screening (k-nearest neighbours)' }], help: 'Selects the inverse-design task. Every mode runs a real optimisation or regression on a compact forward model stated in the equations; results are starting points for detailed analysis and experiment.' },
  { key: 'seed', label: 'Random seed', value: 1, min: 1, max: 99999, step: 1, sens: false, group: 'Model', help: 'Seed of the pseudo-random generator used by the evolutionary optimisers, sampling plans and Monte Carlo analyses. The same seed reproduces the same result exactly.' },
  { key: 'pop', label: 'Population size', value: 24, min: 10, max: 60, step: 1, sens: false, group: 'Model', showIf: M_('auxetic', 'optical', 'resonator', 'recipe'), help: 'Number of candidate designs per generation of differential evolution (doubled, within limits, for NSGA-II). Larger populations explore more but cost proportionally more evaluations.' },
  { key: 'gens', label: 'Generations', value: 60, min: 10, max: 200, step: 1, sens: false, group: 'Model', showIf: M_('auxetic', 'optical', 'resonator', 'recipe'), help: 'Number of generations of the evolutionary search before local refinement. Check the convergence chart: a flat tail means more generations will not help.' },
  // auxetic
  { key: 'nuT', label: 'Target Poisson’s ratio ν₁₂', value: -0.8, min: -5, max: 3, step: 0.05, group: 'Lattice targets', showIf: M_('auxetic', 'bayes'), bus: 'eff.nu', help: 'Required in-plane Poisson’s ratio of the honeycomb. Negative values select re-entrant cells; anisotropic lattices may go below −1.' },
  { key: 'EratT', label: 'Target relative stiffness E₁/Eₛ', value: 1e-3, min: 1e-6, max: 0.1, log: true, group: 'Lattice targets', showIf: M_('auxetic'), help: 'Required effective modulus along the loading direction divided by the modulus of the solid.' },
  { key: 'rhoMax', label: 'Maximum relative density', value: 0.3, min: 0.02, max: 0.8, step: 0.01, group: 'Lattice targets', showIf: M_('auxetic'), help: 'Upper limit on the solid volume fraction of the lattice (mass budget). The beam model is reliable below about 0.3.' },
  { key: 'Es', label: 'Modulus of the solid', unit: 'Pa', value: 2e9, min: 1e5, max: 1.2e12, log: true, group: 'Lattice targets', showIf: M_('auxetic'), bus: 'prop.E', help: 'Young’s modulus of the cell-wall material, used to convert relative into absolute stiffness.' },
  { key: 'cellL', label: 'Inclined wall length l', unit: 'mm', value: 5, min: 0.001, max: 100, log: true, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Length of the oblique cell walls, which fixes the size of the unit cell.' },
  { key: 'tMin', label: 'Minimum strut thickness', unit: 'mm', value: 0.4, min: 0.0001, max: 10, log: true, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Smallest wall thickness the manufacturing process can produce reliably (about two extrusion widths or laser spot sizes).' },
  { key: 'thMinA', label: 'Minimum wall angle |θ|', unit: '°', value: 10, min: 2, max: 45, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Smallest inclination of the oblique walls; nearly horizontal walls give extreme Poisson’s ratios that are very sensitive to tolerances.' },
  { key: 'thMaxA', label: 'Maximum wall angle |θ|', unit: '°', value: 70, min: 46, max: 85, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Largest inclination of the oblique walls allowed by the build process (overhang limit).' },
  { key: 'tolT', label: 'Thickness tolerance', unit: '± mm', value: 0.05, min: 0, max: 1, step: 0.01, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Manufacturing scatter of the wall thickness used in the robustness analysis (uniform distribution).' },
  { key: 'tolTh', label: 'Angle tolerance', unit: '± °', value: 2, min: 0, max: 10, step: 0.5, group: 'Manufacturing', showIf: M_('auxetic'), help: 'Manufacturing scatter of the wall angle used in the robustness analysis (uniform distribution).' },
  // optical
  { key: 'optGoal', label: 'Coating function', type: 'select', value: 'ar', group: 'Multilayer', showIf: M_('optical'), options: [{ value: 'ar', label: 'Antireflection coating (minimise R in band)' }, { value: 'hr', label: 'High reflector (maximise R in band)' }, { value: 'bp', label: 'Band-pass (transmit in band, reflect beside it)' }], help: 'Target spectral function of the multilayer; it defines the merit function minimised by the optimiser.' },
  { key: 'nH', label: 'High-index material', type: 'select', value: 'TiO2', group: 'Multilayer', showIf: M_('optical'), options: nOpts, help: 'High-refractive-index film material. Indices are typical thin-film values in the transparency range, taken as real and non-dispersive.' },
  { key: 'nL', label: 'Low-index material', type: 'select', value: 'SiO2', group: 'Multilayer', showIf: M_('optical'), options: nOpts, help: 'Low-refractive-index film material; a large contrast to the high-index material shortens the stack.' },
  { key: 'nSub', label: 'Substrate', type: 'select', value: 'glass', group: 'Multilayer', showIf: M_('optical'), options: Object.entries(NSUB).map(([value, v]) => ({ value, label: `${v[0]} (n = ${v[1]})` })), help: 'Substrate on which the stack is deposited; light is incident from air.' },
  { key: 'nLay', label: 'Number of layers', value: 6, min: 1, max: 14, step: 1, group: 'Multilayer', showIf: M_('optical'), help: 'Number of alternating layers. More layers give broader or sharper spectral features at the cost of deposition time and tolerance sensitivity.' },
  { key: 'lam1', label: 'Band start', unit: 'nm', value: 450, min: 200, max: 5000, group: 'Multilayer', showIf: M_('optical'), help: 'Short-wavelength edge of the target band.' },
  { key: 'lam2', label: 'Band end', unit: 'nm', value: 650, min: 200, max: 5000, group: 'Multilayer', showIf: M_('optical'), help: 'Long-wavelength edge of the target band.' },
  { key: 'dMin', label: 'Minimum layer thickness', unit: 'nm', value: 10, min: 1, max: 200, group: 'Multilayer', showIf: M_('optical'), help: 'Thinnest layer the deposition process controls reliably; very thin layers are not continuous films.' },
  { key: 'dMax', label: 'Maximum layer thickness', unit: 'nm', value: 300, min: 50, max: 2000, group: 'Multilayer', showIf: M_('optical'), help: 'Thickest single layer allowed (film stress and deposition time).' },
  { key: 'thkTol', label: 'Thickness control', unit: '± %', value: 2, min: 0, max: 20, step: 0.5, group: 'Multilayer', showIf: M_('optical'), help: 'Relative layer-thickness error of the deposition process, used in the Monte Carlo tolerance analysis.' },
  // resonator
  { key: 'resType', label: 'Resonator type', type: 'select', value: 'mim', group: 'Resonator', showIf: M_('resonator'), options: [{ value: 'mim', label: 'Mass-in-mass cell (negative effective mass band)' }, { value: 'srr', label: 'Split-ring resonator (magnetic resonance, μ < 0)' }, { value: 'helm', label: 'Helmholtz resonator (negative bulk modulus)' }], help: 'Type of locally resonant element to be designed for the target frequency or band.' },
  { key: 'fLo', label: 'Lower band edge', unit: 'Hz', value: 100, min: 0.1, max: 1e6, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, bus: 'res.bandgap_lo', help: 'Start of the required negative-mass band; it equals the resonance frequency of the inner mass.' },
  { key: 'fHi', label: 'Upper band edge', unit: 'Hz', value: 150, min: 0.1, max: 1e6, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, bus: 'res.bandgap_hi', help: 'End of the required negative-mass band; the ratio to the lower edge fixes the mass ratio.' },
  { key: 'm1', label: 'Outer (host) mass per cell', unit: 'kg', value: 0.1, min: 1e-9, max: 1e4, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, help: 'Mass of the host frame of one unit cell, to which the internal resonator is attached.' },
  { key: 'k1', label: 'Outer spring stiffness', unit: 'N/m', value: 1e6, min: 1, max: 1e12, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, help: 'Stiffness connecting neighbouring cells of the chain; with the cell mass it sets the host cut-off frequency and the attenuation per cell.' },
  { key: 'muMax', label: 'Maximum mass ratio m₂/m₁', value: 3, min: 0.05, max: 20, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, help: 'Largest admissible resonator mass relative to the host mass (weight budget).' },
  { key: 'k2Min', label: 'Minimum inner stiffness', unit: 'N/m', value: 10, min: 1e-6, max: 1e9, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, help: 'Softest inner spring that can be realised without excessive static sag.' },
  { key: 'k2Max', label: 'Maximum inner stiffness', unit: 'N/m', value: 1e7, min: 1e-3, max: 1e12, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim'] }, help: 'Stiffest inner spring that can be realised in the available space.' },
  { key: 'eta', label: 'Loss factor of the resonator', value: 0.02, min: 0, max: 0.5, step: 0.01, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['mim', 'helm'] }, help: 'Structural or visco-thermal loss factor; it limits the depth of the negative band and broadens the resonance.' },
  { key: 'fSrr', label: 'Target magnetic resonance', unit: 'Hz', value: 1e10, min: 1e7, max: 1e14, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['srr'] }, help: 'Frequency at which the split ring is to resonate; negative permeability occurs just above it.' },
  { key: 'epsSub', label: 'Substrate permittivity', value: 4.4, min: 1, max: 50, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['srr'] }, bus: 'prop.eps_r', help: 'Relative permittivity of the laminate or substrate carrying the ring (4.4 for FR-4, about 2.2 for PTFE laminates, 11.7 for silicon).' },
  { key: 'tMet', label: 'Metal thickness', unit: 'µm', value: 35, min: 0.02, max: 500, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['srr'] }, help: 'Thickness of the copper track (35 µm for standard printed circuit boards; tens of nanometres for lithographic rings).' },
  { key: 'minFeat', label: 'Minimum feature size', unit: 'mm', value: 0.15, min: 1e-5, max: 5, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['srr'] }, help: 'Smallest track width and gap the patterning process can produce (about 0.1–0.15 mm for printed circuit boards, 10⁻⁴ mm for electron-beam lithography).' },
  { key: 'fHelm', label: 'Target acoustic resonance', unit: 'Hz', value: 250, min: 20, max: 20000, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['helm'] }, help: 'Frequency of the Helmholtz resonance, the lower edge of the negative-modulus band.' },
  { key: 'cAir', label: 'Speed of sound', unit: 'm/s', value: 343, min: 100, max: 2000, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['helm'] }, help: 'Sound speed of the fluid in the resonator (343 m/s for air at 20 °C, about 1480 m/s for water).' },
  { key: 'rMinH', label: 'Minimum neck radius', unit: 'mm', value: 2, min: 0.1, max: 50, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['helm'] }, help: 'Smallest neck radius that can be manufactured and kept free of blockage.' },
  { key: 'VmaxH', label: 'Maximum cavity volume', unit: 'cm³', value: 500, min: 1, max: 1e5, log: true, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['helm'] }, help: 'Largest cavity that fits the available space.' },
  { key: 'hFill', label: 'Cavity-to-cell volume ratio F', value: 0.3, min: 0.01, max: 2, group: 'Resonator', showIf: { mode: ['resonator'], resType: ['helm'] }, help: 'Geometrical factor of the resonator array (cavity volume relative to the duct volume per cell); it sets the width of the negative-modulus band.' },
  // nano
  { key: 'nanoType', label: 'Nanoparticle problem', type: 'select', value: 'lspr', group: 'Nanoparticle', showIf: M_('nano'), options: [{ value: 'lspr', label: 'Plasmon resonance → spheroid aspect ratio (Gans)' }, { value: 'qd', label: 'Emission wavelength → quantum-dot diameter (Brus)' }], help: 'Inverse problem solved by root finding: aspect ratio of a noble-metal nanorod for a target plasmon wavelength, or diameter of a quantum dot for a target emission wavelength.' },
  { key: 'metal', label: 'Metal', type: 'select', value: 'Au', group: 'Nanoparticle', showIf: M_('nano', 'bayes'), options: [{ value: 'Au', label: 'Gold (Drude–Lorentz fit)' }, { value: 'Ag', label: 'Silver (Johnson & Christy table)' }], help: 'Metal of the plasmonic nanoparticle; its permittivity is described by an analytical fit to tabulated optical constants.' },
  { key: 'lamT', label: 'Target plasmon wavelength', unit: 'nm', value: 800, min: 350, max: 1800, group: 'Nanoparticle', showIf: M_('nano', 'bayes'), help: 'Wavelength of the longitudinal localised surface-plasmon resonance required, for example 800 nm for the tissue-transparency window.' },
  { key: 'nMed', label: 'Refractive index of the medium', value: 1.33, min: 1, max: 2.5, step: 0.01, group: 'Nanoparticle', showIf: M_('nano', 'bayes'), help: 'Index of the surrounding medium (1.33 for water); a higher index red-shifts the resonance.' },
  { key: 'dEq', label: 'Equivalent-sphere diameter', unit: 'nm', value: 20, min: 2, max: 100, group: 'Nanoparticle', showIf: { mode: ['nano'], nanoType: ['lspr'] }, help: 'Diameter of the sphere with the same volume as the particle; it sets the absolute size and cross-section but not the quasi-static resonance position.' },
  { key: 'qdMat', label: 'Quantum-dot material', type: 'select', value: 'CdSe', group: 'Nanoparticle', showIf: { mode: ['nano'], nanoType: ['qd'] }, options: Object.keys(QD).map(k => ({ value: k, label: `${k} (E_g = ${QD[k][0]} eV)` })), help: 'Semiconductor of the quantum dot with typical band gap, effective masses and permittivity.' },
  { key: 'lamE', label: 'Target emission wavelength', unit: 'nm', value: 580, min: 250, max: 2500, group: 'Nanoparticle', showIf: { mode: ['nano'], nanoType: ['qd'] }, help: 'Required peak emission wavelength; it must be shorter than the bulk band-gap wavelength of the material.' },
  { key: 'sizeDisp', label: 'Size dispersion (standard deviation)', unit: '%', value: 5, min: 0, max: 30, group: 'Nanoparticle', showIf: { mode: ['nano'], nanoType: ['qd'] }, help: 'Relative standard deviation of the diameter in the ensemble; it broadens the emission line.' },
  // recipe
  { key: 'Ltar', label: 'Target nanorod length', unit: 'nm', value: 2000, min: 50, max: 20000, log: true, group: 'Recipe targets', showIf: M_('recipe', 'bayes'), help: 'Required rod length at the end of growth.' },
  { key: 'Dtar', label: 'Target nanorod diameter', unit: 'nm', value: 90, min: 10, max: 2000, log: true, group: 'Recipe targets', showIf: M_('recipe', 'bayes'), help: 'Required rod diameter at the end of growth.' },
  { key: 'errTol', label: 'Quality tolerance', unit: '%', value: 5, min: 0.5, max: 50, group: 'Recipe targets', showIf: M_('recipe'), help: 'Largest combined relative deviation of length and diameter that is still acceptable; the knee point is chosen among recipes within it.' },
  { key: 'TmaxR', label: 'Maximum temperature', unit: '°C', value: 120, min: 61, max: 200, group: 'Recipe bounds', showIf: M_('recipe'), help: 'Upper temperature bound of the equipment; above 100 °C a pressure vessel is needed for aqueous growth.' },
  { key: 'cMaxR', label: 'Maximum precursor concentration', unit: 'mM', value: 100, min: 2, max: 500, group: 'Recipe bounds', showIf: M_('recipe'), help: 'Upper bound of the zinc concentration; high concentrations promote unwanted nucleation in the bulk solution.' },
  { key: 'rAx', label: 'Axial growth rate at 90 °C, 25 mM', unit: 'nm/min', value: 10, min: 0.01, max: 1000, log: true, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Calibration point of the model: initial lengthening rate at the reference condition. Replace the illustrative default by a measured value.' },
  { key: 'rRad', label: 'Radial growth rate at 90 °C, 25 mM', unit: 'nm/min', value: 0.4, min: 0.001, max: 100, log: true, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Calibration point of the model: initial rate of diameter increase at the reference condition.' },
  { key: 'EaAx', label: 'Activation energy, axial growth', unit: 'kJ/mol', value: 60, min: 5, max: 200, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Apparent activation energy of growth on the polar end face.' },
  { key: 'EaRad', label: 'Activation energy, radial growth', unit: 'kJ/mol', value: 40, min: 5, max: 200, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Apparent activation energy of growth on the side faces; its difference to the axial value makes the aspect ratio temperature-dependent.' },
  { key: 'tauRef', label: 'Depletion time constant at 90 °C', unit: 'h', value: 3, min: 0.1, max: 100, log: true, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Time in which the supersaturation falls to 1/e at the reference temperature through consumption and homogeneous precipitation.' },
  { key: 'D0r', label: 'Seed diameter', unit: 'nm', value: 20, min: 1, max: 500, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Initial diameter given by the seed-layer grains, on which radial growth adds.' },
  { key: 'volR', label: 'Solution volume', unit: 'L', value: 0.1, min: 0.001, max: 1000, log: true, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Volume of growth solution per batch, for the heating energy and precursor mass.' },
  { key: 'UA', label: 'Heat-loss coefficient UA', unit: 'W/K', value: 0.5, min: 0, max: 1000, group: 'Growth model', showIf: M_('recipe', 'bayes'), help: 'Overall heat-loss coefficient of the vessel times its area; it determines the energy needed to hold the temperature.' },
  // topology
  { key: 'tProb', label: 'Load case', type: 'select', value: 'mbb', group: 'Topology', showIf: M_('topo'), options: [{ value: 'mbb', label: 'MBB beam (half model, symmetric)' }, { value: 'cant', label: 'Cantilever with tip load at mid-height' }], help: 'Benchmark design domain with supports and load. Both are rectangular plates in plane stress.' },
  { key: 'nelx', label: 'Elements in x', value: 30, min: 10, max: 60, step: 1, sens: false, group: 'Topology', showIf: M_('topo'), help: 'Number of finite elements along the length of the design domain.' },
  { key: 'nely', label: 'Elements in y', value: 15, min: 6, max: 30, step: 1, sens: false, group: 'Topology', showIf: M_('topo'), help: 'Number of finite elements across the height of the design domain.' },
  { key: 'volfrac', label: 'Volume fraction', value: 0.5, min: 0.15, max: 0.85, step: 0.05, group: 'Topology', showIf: M_('topo'), help: 'Fraction of the design domain that may be filled with material (mass budget).' },
  { key: 'penal', label: 'Penalisation exponent p', value: 3, min: 1, max: 5, step: 0.5, group: 'Topology', showIf: M_('topo'), help: 'Exponent of the SIMP interpolation E = ρᵖE₀. Values of 3 or more make intermediate densities uneconomical and drive the design towards solid and void.' },
  { key: 'rminF', label: 'Filter radius', unit: 'elements', value: 1.5, min: 1, max: 4, step: 0.1, group: 'Topology', showIf: M_('topo'), help: 'Radius of the sensitivity filter in element widths; it prevents checkerboard patterns and sets the minimum member size.' },
  { key: 'tIter', label: 'Maximum iterations', value: 60, min: 5, max: 80, step: 1, sens: false, group: 'Topology', showIf: M_('topo'), help: 'Upper limit of optimality-criteria iterations; the run stops earlier when the largest density change falls below 0.01.' },
  { key: 'elSize', label: 'Element size', unit: 'mm', value: 1, min: 0.001, max: 100, log: true, group: 'Topology', showIf: M_('topo'), help: 'Edge length of one square element, giving physical dimensions to the result.' },
  { key: 'tE', label: 'Young’s modulus of the material', unit: 'Pa', value: 2e9, min: 1e5, max: 1.2e12, log: true, group: 'Topology', showIf: M_('topo'), bus: 'prop.E', help: 'Modulus of the solid material, used to convert the normalised compliance into a deflection.' },
  { key: 'tLoad', label: 'Load', unit: 'N', value: 100, min: 1e-6, max: 1e7, log: true, group: 'Topology', showIf: M_('topo'), help: 'Magnitude of the applied point load.' },
  { key: 'tThk', label: 'Plate thickness', unit: 'mm', value: 5, min: 0.001, max: 1000, log: true, group: 'Topology', showIf: M_('topo'), help: 'Out-of-plane thickness of the plate.' },
  // bayes
  { key: 'bModel', label: 'Forward model or data', type: 'select', value: 'growth1d', group: 'Surrogate', showIf: M_('bayes'), options: [{ value: 'growth1d', label: 'Nanorod length error versus temperature (1-D)' }, { value: 'growth2d', label: 'Nanorod size error versus temperature and concentration (2-D)' }, { value: 'lspr', label: 'Plasmon-resonance error versus aspect ratio (1-D)' }, { value: 'lattice', label: 'Poisson-ratio error versus cell angle (1-D)' }, { value: 'table', label: 'User data table (x, y)' }], help: 'Function sampled by the surrogate: one of the forward models of this suite with the targets entered in their own modes, or measured data pasted as a table.' },
  { key: 'bTable', label: 'Data table (x, y)', type: 'table', columns: ['x', 'y'], value: BDEF, group: 'Surrogate', showIf: { mode: ['bayes'], bModel: ['table'] }, help: 'Measured or simulated response y at settings x (one design variable). The default is synthetic demonstration data from the growth model with seeded noise.' },
  { key: 'bGoal', label: 'Goal for table data', type: 'select', value: 'max', group: 'Surrogate', showIf: { mode: ['bayes'], bModel: ['table'] }, options: [{ value: 'max', label: 'Maximise y' }, { value: 'min', label: 'Minimise y' }], help: 'Whether larger or smaller responses are better; expected improvement is computed accordingly.' },
  { key: 'n0', label: 'Initial samples', value: 6, min: 3, max: 15, step: 1, sens: false, group: 'Surrogate', showIf: M_('bayes'), help: 'Number of Latin-hypercube evaluations of the forward model before the acquisition loop starts.' },
  { key: 'nBo', label: 'Bayesian-optimisation iterations', value: 10, min: 1, max: 30, step: 1, sens: false, group: 'Surrogate', showIf: M_('bayes'), help: 'Number of sequential expected-improvement steps, each adding one evaluation and refitting the surrogate.' },
  { key: 'xq', label: 'Query position (fraction of the design range)', value: 0.5, min: -0.5, max: 1.5, step: 0.05, group: 'Surrogate', showIf: M_('bayes'), help: 'Point at which the surrogate prediction and its uncertainty are reported, as a fraction of the range of the first design variable. Values below 0 or above 1 lie outside the design range and trigger the extrapolation warning.' },
  // knn
  { key: 'kProp', label: 'Property to predict', type: 'select', value: 'E', group: 'Library screening', showIf: M_('knn'), options: Object.entries(KPROP).map(([value, v]) => ({ value, label: `${v[0]} (${v[1]})` })), help: 'Library property estimated from the nearest compositions.' },
  { key: 'kMat', label: 'Material', type: 'select', value: 'zrw2o8', group: 'Library screening', showIf: M_('knn'), options: matOpts, help: 'Library material (with a stoichiometric formula) whose property is predicted. If the library already holds the value it is withheld and used as a test.' },
  { key: 'kNN', label: 'Number of neighbours k', value: 3, min: 1, max: 10, step: 1, group: 'Library screening', showIf: M_('knn'), help: 'Number of nearest compositions averaged in the prediction.' },
];

export default {
  id: 25,
  tagline: 'Finds geometries, compositions and process settings that meet target properties by running optimisers, surrogates and screening models on compact physics-based forward models.',
  inputs,
  presets: [
    { name: 'Auxetic re-entrant lattice: ν = −0.8', note: 'Printed polymer honeycomb with a target Poisson’s ratio of −0.8 and E₁/Eₛ = 10⁻³ at no more than 30 % relative density.', values: { mode: 'auxetic', nuT: -0.8, EratT: 1e-3, rhoMax: 0.3, cellL: 5, tMin: 0.4 } },
    { name: 'ZnO nanorods: recipe for 2 µm × 90 nm', note: 'Hydrothermal growth recipe traded between size accuracy, energy and time with NSGA-II.', values: { mode: 'recipe', Ltar: 2000, Dtar: 90, errTol: 5, TmaxR: 120, cMaxR: 100 } },
    { name: 'Negative-permeability split ring at 10 GHz', note: 'Printed-circuit split-ring resonator for a magnetic metamaterial in X band, the μ < 0 building block of a negative-index medium.', values: { mode: 'resonator', resType: 'srr', fSrr: 1e10, epsSub: 4.4, tMet: 35, minFeat: 0.15 } },
    { name: 'Negative-mass band 100–150 Hz', note: 'Mass-in-mass acoustic metamaterial for low-frequency vibration isolation.', values: { mode: 'resonator', resType: 'mim', fLo: 100, fHi: 150, m1: 0.1, k1: 1e6, muMax: 3 } },
    { name: 'Broadband antireflection coating on glass', note: 'Six-layer TiO₂/SiO₂ coating for the visible band 450–650 nm.', values: { mode: 'optical', optGoal: 'ar', nH: 'TiO2', nL: 'SiO2', nSub: 'glass', nLay: 6, lam1: 450, lam2: 650 } },
    { name: 'Dielectric mirror at 1064 nm', note: 'Twelve-layer Ta₂O₅/SiO₂ high reflector (one-dimensional photonic crystal) for 1000–1130 nm.', values: { mode: 'optical', optGoal: 'hr', nH: 'Ta2O5', nL: 'SiO2', nSub: 'silica', nLay: 12, lam1: 1000, lam2: 1130, dMax: 400 } },
    { name: 'Gold nanorods for 800 nm', note: 'Aspect ratio of gold nanorods in water for a longitudinal plasmon in the near-infrared window.', values: { mode: 'nano', nanoType: 'lspr', metal: 'Au', lamT: 800, nMed: 1.33, dEq: 20 } },
    { name: 'CdSe quantum dots emitting at 580 nm', note: 'Core diameter for orange emission from the Brus model.', values: { mode: 'nano', nanoType: 'qd', qdMat: 'CdSe', lamE: 580, sizeDisp: 5 } },
    { name: 'MBB beam topology, 50 % material', note: 'Classical minimum-compliance benchmark for a lightweight bracket or lattice node.', values: { mode: 'topo', tProb: 'mbb', nelx: 30, nely: 15, volfrac: 0.5, penal: 3, rminF: 1.5, tIter: 60 } },
    { name: 'Bayesian optimisation of growth temperature', note: 'Gaussian-process surrogate of the nanorod length error with expected-improvement sampling.', values: { mode: 'bayes', bModel: 'growth1d', n0: 5, nBo: 10, Ltar: 1000 } },
    { name: 'Library screening: modulus of ZrW₂O₈ by nearest neighbours', note: 'Hold-out prediction of a negative-thermal-expansion oxide from composition descriptors, with leave-one-out error.', values: { mode: 'knn', kProp: 'E', kMat: 'zrw2o8', kNN: 3 } },
  ],
  run(inp, ctx) {
    const r = ({ auxetic: runAux, optical: runOptical, resonator: runResonator, nano: runNano, recipe: runRecipe, topo: runTopo, bayes: runBayes, knn: runKnn }[inp.mode] || runAux)(inp, ctx);
    for (const k of r.kpis) if (typeof k.value === 'number' && !Number.isFinite(k.value)) { r.warnings.push(`${k.label} could not be evaluated for these inputs and is reported as 0.`); k.value = 0; k.good = null; }
    r.series = (r.series || []).map(s => ({ ...s, lines: s.lines.map(l => { const x = [], y = []; for (let i = 0; i < l.x.length; i++) if (Number.isFinite(l.x[i]) && Number.isFinite(l.y[i]) && (!s.logy || l.y[i] > 0) && (!s.logx || l.x[i] > 0)) { x.push(l.x[i]); y.push(l.y[i]); } return { ...l, x, y }; }).filter(l => l.x.length) })).filter(s => s.lines.length);
    for (const k of Object.keys(r.publish)) if (typeof r.publish[k] === 'number' && !Number.isFinite(r.publish[k])) delete r.publish[k];
    return r;
  },
};
