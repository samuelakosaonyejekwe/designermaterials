// Suite 17 — Quantum, atomistic and molecular simulation.
// Small in-browser atomistic calculations: pair-potential molecular dynamics, Metropolis and kinetic Monte Carlo,
// model-Hamiltonian electronic structure (tight binding, quantum confinement), lattice dynamics, nudged-elastic-band
// pathways and transition-state rates — plus input decks for external first-principles and molecular-dynamics codes.
// This suite does not run density-functional theory itself.
import * as N from '../core/num.js';
import * as chem from '../core/chem.js';
import * as A from '../geometry/atomic.js';

const F = N.fmt, PI = Math.PI, KB = N.K_B / N.Q_E /* eV/K */, TU = 10.180505 /* fs per internal time unit √(u Å²/eV) */, EVA3_GPA = 160.21766, HB2M = 0.0380998 /* ħ²/2m₀, eV nm² */, AMU = 1.66053907e-27;
const kpi = (key, label, value, unit, explain, good = null) => ({ key, label, value, unit, explain, good });
const eqn = (kind, label, text, note) => (note ? { kind, label, text, note } : { kind, label, text });
const rec = (priority, text, basis, sustainability) => (sustainability ? { priority, text, basis, sustainability: true } : { priority, text, basis });
const chk = (name, pass, detail) => ({ name, pass: !!pass, detail });
const safe = (fn, fallback) => { try { return fn(); } catch { return fallback; } };
const fin = (v, d = 0) => (Number.isFinite(v) ? v : d);
const num = (inp, k, d, lo, hi) => N.clamp(fin(+inp[k], d), lo, hi);
const int = (inp, k, d, lo, hi) => Math.round(num(inp, k, d, lo, hi));
const pick = (v, list, d) => (list.includes(v) ? v : d);
const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a || 1; };
const slug = s => String(s).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 40) || 'structure';
const thin = (x, y, max = 400) => { const n = x.length, k = Math.ceil(n / max); if (k <= 1) return { x: Array.from(x), y: Array.from(y) }; const X = [], Y = []; for (let i = 0; i < n; i += k) { X.push(x[i]); Y.push(y[i]); } return { x: X, y: Y }; };
const KEYS = Object.keys(A.CRYSTALS);

// Pair-potential presets. LJ: ε (eV), σ (Å). Morse: D (eV), α (1/Å), r0 (Å) — classic lattice-property-fitted parameters for
// cubic metals. a: experimental fcc lattice parameter (Å); nuMax: highest measured phonon frequency (THz); thD: Debye temperature (K).
const POT = {
  Ar: { name: 'Argon — Lennard-Jones', type: 'lj', e: 0.0104, r: 3.405, a: 5.31, mass: 39.948, nuMax: 2.0, thD: 92 },
  Cu: { name: 'Copper — Morse', type: 'morse', e: 0.3429, al: 1.3588, r: 2.866, a: 3.615, mass: 63.546, nuMax: 7.4, thD: 343 },
  Ag: { name: 'Silver — Morse', type: 'morse', e: 0.3323, al: 1.3690, r: 3.115, a: 4.086, mass: 107.868, nuMax: 5.1, thD: 225 },
  Au: { name: 'Gold — Morse', type: 'morse', e: 0.4753, al: 1.5830, r: 3.024, a: 4.078, mass: 196.967, nuMax: 4.7, thD: 165 },
  Ni: { name: 'Nickel — Morse', type: 'morse', e: 0.4205, al: 1.4199, r: 2.780, a: 3.524, mass: 58.693, nuMax: 8.9, thD: 450 },
  Al: { name: 'Aluminium — Morse', type: 'morse', e: 0.2703, al: 1.1646, r: 3.253, a: 4.050, mass: 26.982, nuMax: 9.7, thD: 428 },
};
const phiOf = pt => pt.type === 'lj' ? r => { const q = (pt.r / r) ** 6; return 4 * pt.e * (q * q - q); } : r => { const u = Math.exp(-pt.al * (r - pt.r)); return pt.e * (u * u - 2 * u); };
const dphiOf = pt => pt.type === 'lj' ? r => { const q = (pt.r / r) ** 6; return -24 * pt.e * (2 * q * q - q) / r; } : r => { const u = Math.exp(-pt.al * (r - pt.r)); return -2 * pt.al * pt.e * (u * u - u); };
const d2phiOf = pt => pt.type === 'lj' ? r => { const q = (pt.r / r) ** 6; return 24 * pt.e * (26 * q * q - 7 * q) / (r * r); } : r => { const u = Math.exp(-pt.al * (r - pt.r)); return 2 * pt.al * pt.al * pt.e * (2 * u * u - u); };
const fccDef = (el, a) => ({ name: `${el} fcc`, formula: el, system: 'cubic', a, b: a, c: a, alpha: 90, beta: 90, gamma: 90, basis: [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]].map(f => ({ el, f })) });
const potLine = pt => pt.type === 'lj' ? `φ(r) = 4ε[(σ/r)¹² − (σ/r)⁶],  ε = ${F(pt.e, 4)} eV, σ = ${F(pt.r, 4)} Å` : `φ(r) = D[e^(−2α(r − r₀)) − 2e^(−α(r − r₀))],  D = ${F(pt.e, 4)} eV, α = ${F(pt.al, 5)} Å⁻¹, r₀ = ${F(pt.r, 4)} Å`;

// Pair forces over a neighbour list (shifted-force truncation). L > 0 selects the cubic minimum-image convention.
// Returns the potential energy; the virial Σ r·F is accumulated in acc[0] inside the loop (no property access after the
// hot loop, which keeps the JIT-compiled loop entry valid on every call).
function mdForces(lst, cnt, x, y, z, fx, fy, fz, L, rc, lj, E, p1, r0, phiC, dC, acc) {
  const rc2 = rc * rc, iL = 1 / L, per = L < 1e29; let ep = 0; fx.fill(0); fy.fill(0); fz.fill(0); acc[0] = 0;
  for (let q = 0; q < cnt; q += 2) {
    const i = lst[q], j = lst[q + 1]; let dx = x[i] - x[j], dy = y[i] - y[j], dz = z[i] - z[j];
    if (per) { dx -= L * Math.floor(dx * iL + 0.5); dy -= L * Math.floor(dy * iL + 0.5); dz -= L * Math.floor(dz * iL + 0.5); }
    const r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < rc2) {
      const r = Math.sqrt(r2); let phi = 0.5, dphi = 0.5;
      if (lj) { const s2 = p1 / r2, s6 = s2 * s2 * s2; phi = 4 * E * (s6 * s6 - s6); dphi = -24 * E * (2 * s6 * s6 - s6) / r; } else { const u = Math.exp(-p1 * (r - r0)); phi = E * (u * u - 2 * u); dphi = -2 * p1 * E * (u * u - u); }
      ep += phi - phiC - (r - rc) * dC; const fr = (dC - dphi) / r; acc[0] += fr * r2;
      fx[i] += fr * dx; fy[i] += fr * dy; fz[i] += fr * dz; fx[j] -= fr * dx; fy[j] -= fr * dy; fz[j] -= fr * dz;
    }
  }
  return ep;
}
// Verlet list into a pre-allocated buffer; returns the number of entries used, or −1 if the buffer is too small.
function mdList(buf, x, y, z, n, L, rl2) {
  const iL = 1 / L, per = L < 1e29, cap = buf.length - 2; let c = 0;
  for (let i = 0; i < n - 1; i++) {
    const xi = x[i], yi = y[i], zi = z[i];
    for (let j = i + 1; j < n; j++) {
      let dx = xi - x[j], dy = yi - y[j], dz = zi - z[j];
      if (per) { dx -= L * Math.floor(dx * iL + 0.5); dy -= L * Math.floor(dy * iL + 0.5); dz -= L * Math.floor(dz * iL + 0.5); }
      if (dx * dx + dy * dy + dz * dz < rl2) { if (c > cap) return -1; buf[c] = i; buf[c + 1] = j; c += 2; }
    }
  }
  return c;
}

// ───────────────────────────── mode 1: molecular dynamics ─────────────────────────────
function runMD(inp) {
  const warnings = [], checks = [], pk = pick(inp.mdPot, Object.keys(POT), 'Ar'), pt = POT[pk], lj = pt.type === 'lj', sys = pick(inp.mdSys, ['bulk', 'cluster'], 'bulk'), ens = pick(inp.mdEns, ['nve', 'berendsen', 'langevin'], 'nve'), pbc = sys === 'bulk';
  const T = num(inp, 'T', 94.4, 1, 3000), dtfs = num(inp, 'dt', 10, 0.2, 20), dt = dtfs / TU, tau = num(inp, 'tau', 300, 10, 5000) / TU, seed = int(inp, 'seed', 1, 1, 99999), rand = N.rng(seed), m = pt.mass;
  let nEq = int(inp, 'nEq', 300, 50, 5000), nProd = int(inp, 'nProd', 1000, 100, 10000);
  const a0 = lj ? 1.0902 * Math.SQRT2 * pt.r : pt.a, dnn0 = a0 / Math.SQRT2, rcPot = lj ? 2.5 * pt.r : pt.r + 4.5 / pt.al, skin = 0.25 * dnn0 / 1.5;
  let L = 0, pos, init, label;
  if (pbc) {
    const nc = int(inp, 'ncell', 3, 2, 5), rel = num(inp, 'dens', 0.75, 0.4, 1.15), a = a0 / Math.cbrt(rel); L = nc * a; pos = [];
    for (let i = 0; i < nc; i++) for (let j = 0; j < nc; j++) for (let k = 0; k < nc; k++) for (const b of [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) pos.push([(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a]);
    label = `${pos.length} ${pk} atoms in a periodic cubic box of ${F(L, 4)} Å (initial fcc lattice, a = ${F(a, 4)} Å, ${F(100 * rel, 3)} % of the static-crystal density)`;
  } else {
    const st = A.nanoparticle(fccDef(pk, a0), { shape: 'sphere', diameter_nm: num(inp, 'mdD', 1.8e-9, 8e-10, 4e-9) * 1e9, maxAtoms: 300 }); pos = st.atoms.map(q => [q.x, q.y, q.z]);
    if (st.truncated) warnings.push(`The cluster was scaled down to ${pos.length} atoms (${F(st.diameter_nm, 3)} nm) to respect the 300-atom limit of the in-browser integrator.`);
    if (pos.length < 4) throw new Error('The cluster contains fewer than four atoms: increase its diameter.');
    label = `free ${pk} cluster of ${pos.length} atoms (${F(st.diameter_nm, 3)} nm sphere carved from fcc, a = ${F(a0, 4)} Å)`;
  }
  const n = pos.length, x = new Float64Array(n), y = new Float64Array(n), z = new Float64Array(n), vx = new Float64Array(n), vy = new Float64Array(n), vz = new Float64Array(n), fx = new Float64Array(n), fy = new Float64Array(n), fz = new Float64Array(n);
  pos.forEach((p, i) => { x[i] = p[0]; y[i] = p[1]; z[i] = p[2]; });
  init = { kind: 'atoms', atoms: pos.map(p => ({ el: pk, x: p[0], y: p[1], z: p[2] })), cell: pbc ? [[L, 0, 0], [0, L, 0], [0, 0, L]] : null, title: `${pk} MD start` };
  let rc = rcPot; if (pbc && rc > L / 2 - skin) { rc = L / 2 - skin; warnings.push(`The potential cut-off was reduced from ${F(rcPot, 3)} to ${F(rc, 3)} Å so that it fits the minimum-image convention in this small box; energies and pressure shift accordingly. Use more cells for a converged cut-off.`); }
  const rl = rc + skin, rl2 = rl * rl, rc2 = rc * rc, phiC = phiOf(pt)(rc), dC = dphiOf(pt)(rc), V = pbc ? L ** 3 : 0, dof = Math.max(1, ens === 'langevin' ? 3 * n : 3 * n - 3);
  const pairsEst = pbc ? Math.min(n * (n - 1) / 2, n * (4 / 3 * PI * rl ** 3 * n / V) / 2) : n * (n - 1) / 2, budget = 1.6e7;
  if (pairsEst * (nEq + nProd) > budget) { const k = budget / (pairsEst * (nEq + nProd)); nEq = Math.max(50, Math.floor(nEq * k)); nProd = Math.max(100, Math.floor(nProd * k)); warnings.push(`Step counts were reduced to ${nEq} equilibration + ${nProd} production steps to keep the run within the in-browser time budget.`); }
  // velocities
  const sv = Math.sqrt(KB * T / m); let px = 0, py = 0, pz = 0; for (let i = 0; i < n; i++) { vx[i] = rand.normal(0, sv); vy[i] = rand.normal(0, sv); vz[i] = rand.normal(0, sv); px += vx[i]; py += vy[i]; pz += vz[i]; }
  for (let i = 0; i < n; i++) { vx[i] -= px / n; vy[i] -= py / n; vz[i] -= pz / n; }
  const kin = () => { let k = 0; for (let i = 0; i < n; i++) k += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]; return 0.5 * m * k; };
  { const s = Math.sqrt(T / Math.max(2 * kin() / (dof * KB), 1e-12)); for (let i = 0; i < n; i++) { vx[i] *= s; vy[i] *= s; vz[i] *= s; } }
  // neighbour list
  let list = new Int32Array(Math.min(n * (n - 1), Math.max(256, Math.ceil(3 * pairsEst) * 2)) + 4), nl = 0, rebuilds = 0; const x0 = new Float64Array(n), y0 = new Float64Array(n), z0 = new Float64Array(n);
  const build = () => { nl = mdList(list, x, y, z, n, Lbox, rl2); if (nl < 0) { list = new Int32Array(n * (n - 1) + 4); nl = mdList(list, x, y, z, n, Lbox, rl2); } x0.set(x); y0.set(y); z0.set(z); rebuilds++; };
  const stale = () => { const lim = (skin / 2) ** 2; for (let i = 0; i < n; i++) { const dx = x[i] - x0[i], dy = y[i] - y0[i], dz = z[i] - z0[i]; if (dx * dx + dy * dy + dz * dz > lim) return true; } return false; };
  const acc = new Float64Array(1), Lbox = pbc ? L : 1e30, p1 = lj ? pt.r * pt.r : pt.al; let epot = 0; const forces = () => { epot = mdForces(list, nl, x, y, z, fx, fy, fz, Lbox, rc, lj, pt.e, p1, pt.r, phiC, dC, acc); };
  build(); forces();
  const total = nEq + nProd, nS = Math.min(200, nProd), sEvery = Math.max(1, Math.floor(nProd / nS)), tS = [], TS = [], EtS = [], EpS = [], PS = [], msdS = [];
  const xs = new Float64Array(n), ys = new Float64Array(n), zs = new Float64Array(n);
  const nb = 120, rmaxR = pbc ? L / 2 : 3.2 * dnn0, drR = rmaxR / nb, hist = new Float64Array(nb); let nR = 0; const rEvery = Math.max(1, Math.floor(nProd / 40));
  const np = n * (n - 1) / 2, lind = !pbc, sr = lind ? new Float64Array(np) : null, sr2 = lind ? new Float64Array(np) : null; let nLs = 0;
  const frames = [], fEvery = Math.max(1, Math.floor(nProd / 30)), wrap = v => (pbc ? v - L * Math.floor(v / L) : v);
  const snap = () => { const at = new Array(n); for (let i = 0; i < n; i++) at[i] = { el: pk, x: +wrap(x[i]).toFixed(3), y: +wrap(y[i]).toFixed(3), z: +wrap(z[i]).toFixed(3) }; return at; };
  const cL = Math.exp(-dt / tau), sL = Math.sqrt((1 - cL * cL) * KB * T / m), hdt = 0.5 * dt / m; let Tsum = 0, Psum = 0, nAvg = 0, maxT = 0, blown = false;
  for (let s = 0; s < total && !blown; s++) {
    for (let i = 0; i < n; i++) { vx[i] += hdt * fx[i]; vy[i] += hdt * fy[i]; vz[i] += hdt * fz[i]; x[i] += dt * vx[i]; y[i] += dt * vy[i]; z[i] += dt * vz[i]; }
    if (stale()) build(); forces();
    for (let i = 0; i < n; i++) { vx[i] += hdt * fx[i]; vy[i] += hdt * fy[i]; vz[i] += hdt * fz[i]; }
    let ke = kin(), Ti = 2 * ke / (dof * KB); if (!Number.isFinite(Ti) || Ti > 1e6) { blown = true; break; }
    const prod = s >= nEq;
    if (!prod || ens === 'berendsen') { const lam = Math.sqrt(Math.max(0.25, 1 + dt / (prod ? tau : Math.min(tau, 20 * dt)) * (T / Math.max(Ti, 1e-9) - 1))); for (let i = 0; i < n; i++) { vx[i] *= lam; vy[i] *= lam; vz[i] *= lam; } ke *= lam * lam; Ti *= lam * lam; }
    else if (ens === 'langevin') { for (let i = 0; i < n; i++) { vx[i] = cL * vx[i] + sL * rand.normal(); vy[i] = cL * vy[i] + sL * rand.normal(); vz[i] = cL * vz[i] + sL * rand.normal(); } ke = kin(); Ti = 2 * ke / (dof * KB); }
    if (s === nEq - 1) { xs.set(x); ys.set(y); zs.set(z); }
    if (!prod) continue;
    const k = s - nEq, ep = epot, P = pbc ? (n * KB * Ti + acc[0] / 3) / V * EVA3_GPA : 0; Tsum += Ti; Psum += P; nAvg++; if (Ti > maxT) maxT = Ti;
    if (k % sEvery === 0) {
      let cx = 0, cy = 0, cz = 0; for (let i = 0; i < n; i++) { cx += x[i] - xs[i]; cy += y[i] - ys[i]; cz += z[i] - zs[i]; } cx /= n; cy /= n; cz /= n;
      let ms = 0; for (let i = 0; i < n; i++) ms += (x[i] - xs[i] - cx) ** 2 + (y[i] - ys[i] - cy) ** 2 + (z[i] - zs[i] - cz) ** 2;
      tS.push((k + 1) * dtfs * 1e-3); TS.push(Ti); EtS.push((ep + ke) / n); EpS.push(ep / n); PS.push(P); msdS.push(ms / n);
    }
    if (k % rEvery === 0) {
      nR++; for (let i = 0; i < n - 1; i++) for (let j = i + 1; j < n; j++) { let dx = x[i] - x[j], dy = y[i] - y[j], dz = z[i] - z[j]; if (pbc) { dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L); } const r = Math.sqrt(dx * dx + dy * dy + dz * dz); if (r < rmaxR) hist[Math.floor(r / drR)] += 2; if (lind) { const q = i * n - i * (i + 1) / 2 + j - i - 1; sr[q] += r; sr2[q] += r * r; } }
      if (lind) nLs++;
    }
    if (k % fEvery === 0 && frames.length < 40) frames.push({ t: +((k + 1) * dtfs * 1e-3).toFixed(4), atoms: snap() });
  }
  if (blown || tS.length < 4) throw new Error('The integration became unstable (temperature diverged). Reduce the time step, lower the temperature or lower the density.');
  const Tm = Tsum / nAvg, Pm = Psum / nAvg, Et0 = N.mean(EtS.slice(0, Math.ceil(EtS.length / 10))), Et1 = N.mean(EtS.slice(-Math.ceil(EtS.length / 10))), Em = N.mean(EtS), drift = (Et1 - Et0) / Math.max(Math.abs(Em), 1e-12), fluct = N.std(EtS) / Math.max(Math.abs(Em), 1e-12);
  const half = Math.floor(tS.length / 2), fit = N.linreg(tS.slice(half), msdS.slice(half)), slope = fit.slope, D = Math.max(slope, 0) / 6 * 1e-8, tProd = nProd * dtfs * 1e-3;
  // radial distribution
  let rho; if (pbc) rho = n / V; else { let cx = 0, cy = 0, cz = 0; for (const p of pos) { cx += p[0]; cy += p[1]; cz += p[2]; } cx /= n; cy /= n; cz /= n; let s2 = 0; for (const p of pos) s2 += (p[0] - cx) ** 2 + (p[1] - cy) ** 2 + (p[2] - cz) ** 2; rho = n / (4 / 3 * PI * Math.max(Math.sqrt(5 / 3 * s2 / n), 1) ** 3); }
  const rR = [], gR = []; for (let b = 0; b < nb; b++) { const r1 = b * drR, r2 = r1 + drR; rR.push(r1 + drR / 2); gR.push(hist[b] / (nR * n * rho * 4 / 3 * PI * (r2 ** 3 - r1 ** 3))); }
  const ipk = N.argmax(gR), rPeak = rR[ipk], gPeak = gR[ipk];
  let lindV; if (lind) { let s = 0, c = 0; for (let q = 0; q < np; q++) { const mr = sr[q] / nLs, v = sr2[q] / nLs - mr * mr; if (mr > 0) { s += Math.sqrt(Math.max(v, 0)) / mr; c++; } } lindV = c ? s / c : 0; } else lindV = Math.sqrt(Math.max(N.mean(msdS.slice(half)), 0) / 2) / dnn0;
  const diffusive = D > 1e-10 && slope * tProd > 0.5 * dnn0 * dnn0 * 0.1, state = lind ? (lindV > 0.1 ? 'liquid-like (melted or evaporating)' : 'solid-like') : (diffusive && lindV > 0.15 ? 'liquid' : 'solid');
  checks.push(chk('Total-energy conservation (NVE)', ens !== 'nve' || Math.abs(drift) < 2e-3, ens === 'nve' ? `relative drift ${F(drift, 2)} over ${F(tProd, 3)} ps, rms fluctuation ${F(fluct, 2)}` : `not applicable: ${ens} thermostat exchanges energy with a heat bath (drift ${F(drift, 2)})`),
    chk('Temperature near the target', Math.abs(Tm - T) / T < (ens === 'nve' ? 0.25 : 0.1), `⟨T⟩ = ${F(Tm, 4)} K, target ${F(T, 4)} K${ens === 'nve' ? ' (NVE: the temperature drifts from the equilibration value as potential and kinetic energy equilibrate)' : ''}`),
    chk('Momentum conservation', ens === 'langevin' || (() => { let a = 0, b = 0, c = 0; for (let i = 0; i < n; i++) { a += vx[i]; b += vy[i]; c += vz[i]; } return Math.hypot(a, b, c) / n < 1e-6 * sv + 1e-12; })(), ens === 'langevin' ? 'not conserved by construction with a Langevin bath' : 'centre-of-mass velocity stays zero to round-off'),
    chk('Time step resolves the fastest vibration', dtfs < 0.04 * 2 * PI * Math.sqrt(m / Math.max(d2phiOf(pt)(lj ? 1.1225 * pt.r : pt.r), 1e-9)) * TU, `Δt = ${F(dtfs, 3)} fs; pair oscillation period ≈ ${F(2 * PI * Math.sqrt(m / (2 * d2phiOf(pt)(lj ? 1.1225 * pt.r : pt.r))) * TU, 4)} fs`));
  if (ens === 'langevin') warnings.push('Langevin friction damps momentum, so the self-diffusion coefficient is lowered relative to Newtonian dynamics; use NVE production runs for transport coefficients.');
  if (ens === 'berendsen') warnings.push('The Berendsen thermostat does not sample the canonical ensemble exactly (it suppresses kinetic-energy fluctuations); it is intended for equilibration.');
  if (!pbc) warnings.push('For a free cluster the mean-squared displacement includes surface-atom mobility and slow rigid rotation of the whole cluster, so a non-zero "diffusion coefficient" does not by itself indicate melting; use the Lindemann index.');
  if (!lj) warnings.push('Morse pair potentials are qualitative for metals (no many-body embedding): melting points, vacancy energies and surface relaxations are not reliable. Use an embedded-atom potential in LAMMPS for quantitative work.');
  warnings.push(`Small system (${n} atoms, ${F(tProd, 3)} ps): statistical uncertainty of D and P is of the order of 10–30 %, and finite-size effects lower D in periodic cells.`);
  const final = snap(), geom = { kind: 'atoms', atoms: final, cell: init.cell, title: `${pk} MD, ${F(Tm, 4)} K, t = ${F(tProd, 3)} ps`, frames, timeUnit: 'ps' };
  const traj = frames.map(f => A.toXYZ({ atoms: f.atoms }, { title: `${pk} MD frame t = ${f.t} ps` })).join('');
  const rcL = F(rc, 5), lmp = ['# LAMMPS input matching the in-browser run (shifted-force pair potential)', 'units           metal', 'atom_style      atomic', `boundary        ${pbc ? 'p p p' : 'f f f'}`, `read_data       ${slug(pk)}_md.lammps-data`, lj ? `pair_style      lj/smooth/linear ${rcL}` : `pair_style      morse/smooth/linear ${rcL}`, lj ? `pair_coeff      1 1 ${pt.e} ${pt.r}` : `pair_coeff      1 1 ${pt.e} ${pt.al} ${pt.r}`, `mass            1 ${m}`, `neighbor        ${F(skin, 3)} bin`, `velocity        all create ${F(T, 5)} ${seed} mom yes rot yes dist gaussian`, `timestep        ${F(dtfs * 1e-3, 4)}`, '# equilibration', 'fix             eq all nve', `fix             th all temp/berendsen ${F(T, 5)} ${F(T, 5)} ${F(20 * dtfs * 1e-3, 3)}`, 'thermo          100', `run             ${nEq}`, 'unfix           th', '# production', ens === 'berendsen' ? `fix             th all temp/berendsen ${F(T, 5)} ${F(T, 5)} ${F(tau * TU * 1e-3, 3)}` : ens === 'langevin' ? `fix             th all langevin ${F(T, 5)} ${F(T, 5)} ${F(tau * TU * 1e-3, 3)} ${seed}` : '# NVE: no thermostat', 'compute         msd all msd com yes', 'compute         rdf all rdf 120', 'fix             rdfout all ave/time 10 100 1000 c_rdf[*] file rdf.dat mode vector', 'thermo_style    custom step temp pe ke etotal press c_msd[4]', `dump            trj all xyz ${fEvery} traj.xyz`, `run             ${nProd}`].join('\n') + '\n';
  const msdT = thin(tS, msdS, 200);
  return {
    summary: `Velocity-Verlet molecular dynamics of a ${label} with the ${pt.name} potential (${ens.toUpperCase()} production, Δt = ${F(dtfs, 3)} fs, ${nEq} + ${nProd} steps). Mean temperature ${F(Tm, 4)} K${pbc ? `, virial pressure ${F(Pm * 1e3, 3)} MPa` : ''}, potential energy ${F(N.mean(EpS), 4)} eV/atom; the system is ${state} with a self-diffusion coefficient of ${F(D, 3)} m²/s and a Lindemann ${lind ? 'index' : 'ratio'} of ${F(lindV, 3)}. Relative total-energy drift ${F(drift, 2)}. Level of theory: classical empirical pair potential, truncated at ${F(rc, 3)} Å (shifted-force).`,
    kpis: [
      kpi('T_mean', 'Mean temperature', Tm, 'K', `Time average of the kinetic temperature 2E_kin/(N_dof k_B) over the production run (${dof} degrees of freedom); the maximum instantaneous value was ${F(maxT, 4)} K.`),
      kpi('E_pot', 'Potential energy', N.mean(EpS), 'eV/atom', `Mean pair-interaction energy per atom with the truncated, shifted-force potential; for the static ${pk} crystal the full-range value is about ${lj ? F(-8.61 * pt.e, 3) : '−E_coh'} eV.`),
      kpi('drift', 'Relative total-energy drift', drift, '–', `Change of the total energy between the first and last tenth of the production run divided by its mean. In NVE it measures integration accuracy; with a thermostat energy is exchanged by design.`, ens === 'nve' ? Math.abs(drift) < 2e-3 : null),
      kpi('P', 'Virial pressure', pbc ? Pm * 1e9 : 0, 'Pa', pbc ? `P = (N k_B T + ⅓Σ r·F)/V averaged over the run; no long-range tail correction is applied for the truncated potential${lj ? ` (for a uniform fluid it would add about ${F(16 / 3 * PI * (n / V) ** 2 * pt.e * pt.r ** 3 * (2 / 3 * (pt.r / rc) ** 9 - (pt.r / rc) ** 3) * EVA3_GPA * 1e3, 3)} MPa)` : ''}. Standard deviation ${F(N.std(PS) * 1e3, 3)} MPa.` : 'A free cluster in vacuum has no defined bulk pressure; the value is reported as zero.'),
      kpi('D', 'Self-diffusion coefficient', D, 'm²/s', `Einstein relation D = (1/6) d⟨|Δr|²⟩/dt from a linear fit to the second half of the mean-squared displacement (slope ${F(slope, 3)} Å²/ps, R² = ${F(fit.r2 ?? 0, 3)}). ${diffusive ? 'The displacement grows linearly: liquid-like diffusion.' : 'No net diffusion is resolved on this time scale: the atoms vibrate about fixed sites.'}`),
      kpi('lindemann', lind ? 'Lindemann index' : 'Lindemann ratio', lindV, '–', lind ? 'Mean relative root-mean-square fluctuation of all pair distances; values above about 0.1 indicate a melted (or evaporating) cluster.' : `Root-mean-square displacement √(MSD/2) over the nearest-neighbour distance ${F(dnn0, 4)} Å; crystals melt when it reaches roughly 0.1–0.15, and it grows without bound in a liquid.`, null),
      kpi('g_peak', 'First RDF peak height', gPeak, '–', `Maximum of g(r) at r = ${F(rPeak, 4)} Å, the mean nearest-neighbour distance. Liquids show a peak of about 2.5–3 followed by damped oscillations; crystals show sharp, higher peaks.`),
      kpi('r_peak', 'Nearest-neighbour distance', rPeak, 'Å', 'Position of the first maximum of the radial distribution function averaged over the production run.'),
      kpi('state', 'Aggregate state', state, '', 'Classification from the mean-squared displacement (diffusive or bounded) and the Lindemann criterion; small periodic cells can stay superheated or supercooled, so treat it as indicative.'),
      kpi('natoms', 'Atoms', n, '–', `Number of atoms integrated; neighbour list rebuilt ${rebuilds} times with a skin of ${F(skin, 3)} Å.`),
    ],
    series: [
      { key: 'T_t', title: 'Temperature vs time', xlabel: 'time (ps)', ylabel: 'temperature (K)', lines: [{ name: 'instantaneous T', x: tS, y: TS }, { name: 'target', x: [tS[0], tS[tS.length - 1]], y: [T, T] }], explain: 'Kinetic temperature during production. Fluctuations of relative size √(2/3N) are intrinsic to a small system; a systematic drift means the system was not equilibrated.' },
      { key: 'E_t', title: 'Energy vs time', xlabel: 'time (ps)', ylabel: 'energy (eV/atom)', lines: [{ name: 'total', x: tS, y: EtS }, { name: 'potential', x: tS, y: EpS }], explain: 'Total and potential energy per atom. In an NVE run the total must be flat: a visible slope signals too large a time step or a discontinuous potential.' },
      { key: 'rdf', title: 'Radial distribution function', xlabel: 'r (Å)', ylabel: 'g(r)', lines: [{ name: 'g(r)', x: rR, y: gR }], explain: pbc ? 'Pair distribution averaged over the production run. Sharp separated peaks indicate crystalline order; a broad first peak and g → 1 at long range indicate a liquid.' : 'Pair distribution of the cluster normalised to the density of its equivalent sphere; it decays with r because surface atoms have fewer neighbours.' },
      { key: 'msd', title: 'Mean-squared displacement', xlabel: 'time (ps)', ylabel: 'MSD (Å²)', lines: [{ name: 'MSD', x: msdT.x, y: msdT.y }, { name: 'linear fit (2nd half)', x: [tS[half], tS[tS.length - 1]], y: [fit.intercept + fit.slope * tS[half], fit.intercept + fit.slope * tS[tS.length - 1]] }], explain: 'Displacement from the positions at the start of production with centre-of-mass motion removed. A plateau means a solid; a straight line means diffusion with slope 6D.' },
    ],
    equations: [
      eqn('formula', 'System', `${chem.prettyFormula(pk)}, M = ${F(chem.molarMass(pk), 6)} g/mol, N = ${n}`),
      eqn('governing', 'Newtonian equations of motion', 'm d²r_i/dt² = F_i = −∇_i Σ_{j≠i} φ(r_ij)'), eqn('constitutive', 'Pair potential', potLine(pt), `Shifted-force truncation at r_c = ${F(rc, 4)} Å: φ_sf(r) = φ(r) − φ(r_c) − (r − r_c) φ′(r_c).`),
      eqn('governing', 'Velocity-Verlet integrator', 'v(t+½Δt) = v(t) + ½Δt F(t)/m;  r(t+Δt) = r(t) + Δt v(t+½Δt);  v(t+Δt) = v(t+½Δt) + ½Δt F(t+Δt)/m'),
      eqn('governing', ens === 'langevin' ? 'Langevin thermostat' : ens === 'berendsen' ? 'Berendsen thermostat' : 'Microcanonical ensemble', ens === 'langevin' ? 'm dv/dt = F − (m/τ) v + √(2 m k_BT/τ) ξ(t)' : ens === 'berendsen' ? 'v ← v √[1 + (Δt/τ)(T₀/T − 1)]' : 'E = Σ ½mv² + Σ φ = constant (no thermostat during production)'),
      eqn('governing', 'Einstein relation and virial pressure', 'D = lim (1/6t) ⟨|r(t) − r(0)|²⟩;   P = [N k_BT + ⅓ Σ_{i<j} r_ij·F_ij]/V'),
    ],
    geometry: geom, checks, warnings,
    recommendations: [
      rec('high', ens === 'nve' ? (Math.abs(drift) < 2e-3 ? `Energy is conserved to ${F(Math.abs(drift), 2)} relative: the ${F(dtfs, 3)} fs time step is adequate. Keep it when moving this set-up to LAMMPS with the exported script.` : `The energy drift of ${F(drift, 2)} is too large: halve the time step (currently ${F(dtfs, 3)} fs) before trusting any property.`) : 'Switch production to NVE (after thermostatted equilibration) to verify energy conservation and to measure diffusion without thermostat artefacts.', 'NVE energy-conservation test.'),
      rec('medium', `Statistics from ${n} atoms and ${F(tProd, 3)} ps are indicative only. For publishable D or melting data run the exported LAMMPS deck with at least a few thousand atoms and nanoseconds, and average several seeds.`, 'Finite-size and sampling uncertainty.'),
      rec('medium', lj ? 'Lennard-Jones argon is a good model of rare-gas liquids and a generic model of simple fluids; map results to other noble gases by rescaling ε and σ (corresponding states).' : 'For metals replace the Morse pair potential by an embedded-atom (EAM) potential from a curated repository in the LAMMPS script; pair potentials give wrong vacancy and surface energies.', 'Transferability limits of the potential.'),
      rec('low', 'Prototype time step, cut-off and thermostat settings here on a few hundred atoms before submitting large production runs: failed or mis-configured cluster jobs are a major avoidable energy cost in atomistic simulation.', 'Compute-energy stewardship.', true),
    ],
    publish: { 'material.formula': pk, 'prop.D_diff': D, 'proc.T_C': Tm - 273.15, 'sim.pressure': pbc ? Pm * 1e9 : 0, 'sim.E_pot_eV': N.mean(EpS), 'sim.lindemann': lindV, 'prop.density': pbc ? n * m * AMU / (V * 1e-30) : n * m * AMU / (n / rho * 1e-30) },
    exports: [
      { name: `${slug(pk)}_md_traj.xyz`, mime: 'chemical/x-xyz', content: traj, description: `XYZ trajectory, ${frames.length} frames over ${F(tProd, 3)} ps (coordinates in Å${pbc ? ', wrapped into the box' : ''}).` },
      { name: `${slug(pk)}_md.lammps-data`, mime: 'text/plain', content: A.toLAMMPSData(init), description: 'LAMMPS data file of the starting configuration.' },
      { name: `in.${slug(pk)}_md.lammps`, mime: 'text/plain', content: lmp, description: 'LAMMPS input script reproducing this run (same potential, cut-off, time step, ensemble and step counts).' },
    ],
    numerics: { method: `Velocity Verlet, Verlet neighbour list (skin ${F(skin, 3)} Å), ${ens === 'nve' ? 'NVE' : ens === 'berendsen' ? 'Berendsen NVT' : 'Langevin NVT'} production after Berendsen equilibration`, dof: 3 * n, steps: total, iterations: rebuilds, residual: Math.abs(drift), notes: `Seed ${seed}; internal units Å, eV, u (time unit ${TU} fs).` },
  };
}

// ───────────────────────────── mode 2: Metropolis Monte Carlo ─────────────────────────────
const agm = (a, b) => { for (let i = 0; i < 40 && Math.abs(a - b) > 1e-15; i++) [a, b] = [(a + b) / 2, Math.sqrt(a * b)]; return a; };
const TC2D = 2 / Math.log(1 + Math.SQRT2);
function isingExactU(T) { const K = 1 / T, k = 2 * Math.sinh(2 * K) / Math.cosh(2 * K) ** 2, K1 = PI / (2 * agm(1, Math.sqrt(Math.max(1 - k * k, 1e-300)))); return -1 / Math.tanh(2 * K) * (1 + 2 / PI * (2 * Math.tanh(2 * K) ** 2 - 1) * K1); }
function runIsing(inp) {
  const L = int(inp, 'mcL', 16, 8, 48), n = L * L, Jm = num(inp, 'Jmev', 10, 0.1, 500), Tu = num(inp, 'mcT', 2.0, 0.5, 5), H = num(inp, 'mcH', 0, -1, 1), seed = int(inp, 'seed', 1, 1, 99999), rand = N.rng(seed), warnings = [], checks = [];
  let sweeps = int(inp, 'mcSweeps', 600, 100, 4000); const nT = 17; if (n * sweeps * 1.4 * (nT + 2) > 2.5e7) { sweeps = Math.max(100, Math.floor(2.5e7 / (n * 1.4 * (nT + 2)))); warnings.push(`Sweeps per temperature were reduced to ${sweeps} to keep the run within the in-browser time budget.`); }
  const eq = Math.max(50, Math.floor(0.4 * sweeps)), s = new Int8Array(n).fill(1), up = new Int32Array(n), dn = new Int32Array(n), lf = new Int32Array(n), rt = new Int32Array(n);
  for (let i = 0; i < n; i++) { const xx = i % L, yy = (i - xx) / L; rt[i] = yy * L + (xx + 1) % L; lf[i] = yy * L + (xx + L - 1) % L; up[i] = ((yy + 1) % L) * L + xx; dn[i] = ((yy + L - 1) % L) * L + xx; }
  const sweep = T => { let acc = 0; for (let q = 0; q < n; q++) { const i = rand.int(n), dE = 2 * s[i] * (s[rt[i]] + s[lf[i]] + s[up[i]] + s[dn[i]] + H); if (dE <= 0 || rand() < Math.exp(-dE / T)) { s[i] = -s[i]; acc++; } } return acc / n; };
  const measure = () => { let e = 0, mg = 0; for (let i = 0; i < n; i++) { e -= s[i] * (s[rt[i]] + s[up[i]]) + H * s[i]; mg += s[i]; } return [e, mg]; };
  const Ts = N.linspace(1.2, 3.6, nT), mT = [], cT = [], xT = [], eT = [];
  for (const T of Ts) { for (let k = 0; k < eq; k++) sweep(T); let e1 = 0, e2 = 0, m1 = 0, m2 = 0; for (let k = 0; k < sweeps; k++) { sweep(T); const [e, mg] = measure(), am = Math.abs(mg); e1 += e; e2 += e * e; m1 += am; m2 += am * am; } e1 /= sweeps; e2 /= sweeps; m1 /= sweeps; m2 /= sweeps; mT.push(m1 / n); eT.push(e1 / n); cT.push((e2 - e1 * e1) / (n * T * T)); xT.push((m2 - m1 * m1) / (n * T)); }
  // diagnostics at the chosen temperature from a random (hot) start
  for (let i = 0; i < n; i++) s[i] = rand() < 0.5 ? 1 : -1;
  const eS = [], mS = [], tS = []; let acc = 0; const nd = Math.min(sweeps + eq, 2000);
  for (let k = 0; k < nd; k++) { acc += sweep(Tu); const [e, mg] = measure(); eS.push(e / n); mS.push(Math.abs(mg) / n); tS.push(k + 1); }
  acc /= nd; const tail = mS.slice(Math.floor(nd / 2)), mm = N.mean(tail), vv = N.variance(tail); let tauI = 0.5;
  if (vv > 1e-12) for (let lag = 1; lag < Math.min(200, tail.length / 4); lag++) { let c = 0; for (let i = 0; i + lag < tail.length; i++) c += (tail[i] - mm) * (tail[i + lag] - mm); c /= (tail.length - lag) * vv; if (c <= 0) break; tauI += c; }
  const eqT = (() => { const ef = N.mean(eS.slice(Math.floor(nd / 2))), sd = N.std(eS.slice(Math.floor(nd / 2))) + 1e-9; for (let k = 0; k < nd; k++) if (Math.abs(eS[k] - ef) < 2 * sd) return k + 1; return nd; })();
  const iC = N.argmax(cT), iX = N.argmax(xT), TcC = Ts[iC], TcX = Ts[iX], TcK = TC2D * Jm * 1e-3 / KB, onsM = T => (T < TC2D ? (1 - Math.sinh(2 / T) ** -4) ** 0.125 : 0), map = new Float64Array(n); for (let i = 0; i < n; i++) map[i] = s[i];
  const uEx = Ts.map(isingExactU), iLow = 0, iHigh = nT - 1;
  checks.push(chk('Susceptibility peak near the Onsager critical temperature', Math.abs(TcX - TC2D) < 0.35, `χ peaks at ${F(TcX, 4)} J/k_B on a ${L}×${L} lattice; exact infinite-lattice value ${F(TC2D, 5)} (finite-size shift ∝ 1/L)`),
    chk('Energy against the exact Onsager solution (away from T_c)', H !== 0 || (Math.abs(eT[iLow] - uEx[iLow]) < 0.03 && Math.abs(eT[iHigh] - uEx[iHigh]) < 0.05), H !== 0 ? 'not applicable with an applied field' : `u(${Ts[iLow]}) = ${F(eT[iLow], 4)} vs ${F(uEx[iLow], 4)};  u(${Ts[iHigh]}) = ${F(eT[iHigh], 4)} vs ${F(uEx[iHigh], 4)} J per spin`),
    chk('Acceptance ratio is physical', acc > 0 && acc < 1, `${F(100 * acc, 3)} % of proposed flips accepted at T = ${F(Tu, 3)} J/k_B`), chk('Run longer than the autocorrelation time', nd / 2 > 20 * tauI, `τ_int ≈ ${F(tauI, 3)} sweeps, ${Math.floor(nd / 2)} measurement sweeps`));
  if (Math.abs(Tu - TC2D) < 0.2) warnings.push('The chosen temperature is close to the critical point: critical slowing down makes single-spin-flip sampling inefficient (long autocorrelation time). Cluster algorithms are needed for accurate critical properties.');
  if (H !== 0) warnings.push('With an applied field the transition is smeared; the Onsager comparison applies only at zero field.');
  return {
    summary: `Metropolis Monte Carlo of the 2-D Ising model on a ${L}×${L} periodic lattice (${sweeps} sweeps per temperature). The susceptibility peaks at ${F(TcX, 4)} J/k_B and the specific heat at ${F(TcC, 4)} J/k_B, against the exact Onsager value ${F(TC2D, 5)} J/k_B; for J = ${F(Jm, 3)} meV this corresponds to T_c = ${F(TcK, 4)} K. At T = ${F(Tu, 3)} J/k_B the magnetisation per spin is ${F(mm, 3)} with ${F(100 * acc, 3)} % acceptance.`,
    kpis: [
      kpi('Tc_chi', 'Critical temperature from χ peak', TcX, 'J/k_B', `Temperature of the maximum magnetic susceptibility on the simulated grid (spacing ${F(Ts[1] - Ts[0], 2)}). Finite lattices shift the peak above the exact value by an amount ∝ 1/L.`, Math.abs(TcX - TC2D) < 0.35),
      kpi('Tc_c', 'Critical temperature from specific-heat peak', TcC, 'J/k_B', 'Temperature of the maximum of c = (⟨E²⟩ − ⟨E⟩²)/(N k_B T²); it converges to T_c faster than the susceptibility peak.'),
      kpi('Tc_exact', 'Onsager critical temperature', TC2D, 'J/k_B', 'Exact result k_B T_c/J = 2/ln(1 + √2) for the infinite square lattice — the benchmark for the simulation.'),
      kpi('Tc_K', 'Critical temperature for the given J', TcK, 'K', `Exact T_c converted with the exchange constant J = ${F(Jm, 3)} meV. A real 2-D magnet needs anisotropy for this Ising behaviour to apply.`),
      kpi('m', 'Magnetisation per spin at the chosen T', mm, '–', `Mean |M|/N over the second half of the diagnostic run at T = ${F(Tu, 3)} J/k_B; the exact spontaneous magnetisation is ${F(onsM(Tu), 4)}.`),
      kpi('accept', 'Acceptance ratio', acc, '–', 'Fraction of proposed single-spin flips accepted at the chosen temperature. It falls towards zero deep in the ordered phase, where the dynamics freeze.'),
      kpi('tau', 'Integrated autocorrelation time', tauI, 'sweeps', 'Number of sweeps between statistically independent magnetisation samples; the effective sample size is the number of sweeps divided by 2τ.'),
      kpi('t_eq', 'Equilibration time', eqT, 'sweeps', 'First sweep at which the energy of the run started from random spins enters the ±2σ band of its final value; discard at least this many sweeps before measuring.'),
    ],
    series: [
      { key: 'm_T', title: 'Magnetisation vs temperature', xlabel: 'k_B T / J', ylabel: '|M| / N', lines: [{ name: `Monte Carlo ${L}×${L}`, x: Ts, y: mT }, { name: 'Onsager (exact, infinite lattice)', x: N.linspace(1.2, 3.6, 121), y: N.linspace(1.2, 3.6, 121).map(onsM) }], explain: 'Order parameter across the ferromagnetic transition. The finite lattice rounds the sharp exact curve and leaves a tail above T_c of order 1/√N.' },
      { key: 'c_T', title: 'Specific heat and susceptibility', xlabel: 'k_B T / J', ylabel: 'per spin', lines: [{ name: 'specific heat c/k_B', x: Ts, y: cT }, { name: 'susceptibility χ J', x: Ts, y: xT }], explain: 'Fluctuation estimators of the response functions; both peak near the critical temperature and the peaks sharpen and grow with lattice size.' },
      { key: 'u_T', title: 'Internal energy vs temperature', xlabel: 'k_B T / J', ylabel: 'u / J per spin', lines: [{ name: 'Monte Carlo', x: Ts, y: eT }, { name: 'Onsager (exact)', x: Ts, y: uEx }], explain: 'Energy per spin against the exact infinite-lattice result; agreement away from T_c verifies the sampling.' },
      { key: 'eq', title: 'Equilibration at the chosen temperature', xlabel: 'sweep', ylabel: 'per spin', lines: [{ name: 'energy u/J', ...thin(tS, eS, 300) }, { name: '|M|/N', ...thin(tS, mS, 300) }], explain: 'Relaxation from a random spin configuration. Measurements are valid only after both curves have stopped drifting.' },
    ],
    maps: [{ key: 'spins', title: `Spin configuration at T = ${F(Tu, 3)} J/k_B`, nx: L, ny: L, data: map, x0: 0, x1: L, y0: 0, y1: L, xlabel: 'x (sites)', ylabel: 'y (sites)', unit: 'spin', diverging: true, explain: 'Final spin configuration (+1 / −1). Below T_c one orientation dominates; near T_c domains of all sizes coexist; above T_c the pattern is random.' }],
    equations: [
      eqn('governing', 'Ising Hamiltonian', 'H = −J Σ_⟨ij⟩ s_i s_j − h Σ_i s_i,  s_i = ±1'), eqn('governing', 'Metropolis acceptance', 'P(accept) = min[1, exp(−ΔE/k_BT)]', 'Satisfies detailed balance with the Boltzmann distribution.'),
      eqn('constitutive', 'Fluctuation relations', 'c = (⟨E²⟩ − ⟨E⟩²)/(N k_BT²),   χ = (⟨M²⟩ − ⟨|M|⟩²)/(N k_BT)'), eqn('formula', 'Onsager solution', 'k_BT_c/J = 2/ln(1 + √2) ≈ 2.269;   m = [1 − sinh⁻⁴(2J/k_BT)]^(1/8) for T < T_c'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', `Discard at least ${Math.max(eqT, Math.ceil(10 * tauI))} sweeps and space measurements by about ${Math.ceil(2 * tauI)} sweeps at this temperature to obtain independent samples.`, 'Equilibration time and integrated autocorrelation time of this run.'),
      rec('medium', `Locate T_c precisely by repeating for several lattice sizes and extrapolating the peak position against 1/L (finite-size scaling); the present ${L}×${L} estimate is ${F(100 * (TcX - TC2D) / TC2D, 3)} % from the exact value.`, 'Finite-size shift of the susceptibility peak.'),
      rec('medium', 'Map exchange constants from first-principles total-energy differences onto J to predict ordering temperatures of magnetic nanostructures; the same lattice model describes order–disorder alloys and lattice-gas adsorption.', 'Ising model as an effective Hamiltonian.'),
      rec('low', 'Use cluster (Wolff) updates or parallel tempering near T_c instead of brute-force longer single-flip runs: they cut the number of sweeps, and therefore the compute energy, by orders of magnitude.', 'Critical slowing down of local updates.', true),
    ],
    publish: { 'mag.Tc': TcK, 'sim.Tc_reduced': TcX, 'sim.acceptance': acc, 'sim.magnetisation': mm },
    numerics: { method: 'Metropolis single-spin-flip Monte Carlo, random-site updates, periodic boundaries', dof: n, steps: (sweeps + eq) * nT + nd, iterations: n, residual: Math.abs(TcX - TC2D), notes: `Seed ${seed}; ${eq} equilibration + ${sweeps} measurement sweeps at each of ${nT} temperatures (annealed upward from the ordered state).` },
  };
}
function runLatticeGas(inp) {
  const L = int(inp, 'mcL', 16, 8, 48), n = L * L, Ea = num(inp, 'Eads', 0.5, 0.01, 3), w = num(inp, 'wmev', 20, -200, 200) * 1e-3, T = num(inp, 'Tads', 300, 20, 2000), mass = num(inp, 'mAds', 28, 1, 400), P = num(inp, 'Pgas', 1e3, 1e-10, 1e8), seed = int(inp, 'seed', 1, 1, 99999), rand = N.rng(seed), warnings = [], checks = [];
  let sweeps = int(inp, 'mcSweeps', 600, 100, 4000); const nP = 15; if (n * sweeps * 1.5 * (nP + 1) > 2.5e7) { sweeps = Math.max(100, Math.floor(2.5e7 / (n * 1.5 * (nP + 1)))); warnings.push(`Sweeps per pressure were reduced to ${sweeps} to keep the run within the in-browser time budget.`); }
  const kT = KB * T, wt = w / kT, lam = N.H_PLANCK / Math.sqrt(2 * PI * mass * AMU * N.K_B * T), lnP0 = Math.log(N.K_B * T / lam ** 3), lnK = Ea / kT - lnP0, z = 4;
  const occ = new Uint8Array(n), nbr = new Int32Array(4 * n); for (let i = 0; i < n; i++) { const xx = i % L, yy = (i - xx) / L; nbr[4 * i] = yy * L + (xx + 1) % L; nbr[4 * i + 1] = yy * L + (xx + L - 1) % L; nbr[4 * i + 2] = ((yy + 1) % L) * L + xx; nbr[4 * i + 3] = ((yy + L - 1) % L) * L + xx; }
  let cov = 0;
  const sweep = mu => { let acc = 0; for (let q = 0; q < n; q++) { const i = rand.int(n), nn = occ[nbr[4 * i]] + occ[nbr[4 * i + 1]] + occ[nbr[4 * i + 2]] + occ[nbr[4 * i + 3]], dE = (occ[i] ? 1 : -1) * (wt * nn + mu); if (dE <= 0 || rand() < Math.exp(-dE)) { cov += occ[i] ? -1 : 1; occ[i] ^= 1; acc++; } } return acc / n; };   // dE in k_BT: adsorb −(w·nn + μ), desorb +(w·nn + μ)
  const run = (mu, eq, ms) => { for (let k = 0; k < eq; k++) sweep(mu); let th = 0, a = 0, th2 = 0; for (let k = 0; k < ms; k++) { a += sweep(mu); const t = cov / n; th += t; th2 += t * t; } return { theta: th / ms, acc: a / ms, sd: Math.sqrt(Math.max(th2 / ms - (th / ms) ** 2, 0)) }; };
  const mid = -z * wt / 2, mus = N.linspace(mid - 7, mid + 7, nP), eq = Math.max(50, Math.floor(0.4 * sweeps)), thMC = mus.map(mu => run(mu, eq, sweeps).theta);
  const log10P = mu => (mu - lnK) / Math.LN10, muU = Math.log(P) + lnK, ur = run(muU, 2 * eq, sweeps), thU = ur.theta, map = new Float64Array(n); for (let i = 0; i < n; i++) map[i] = occ[i];
  const thG = N.linspace(0.005, 0.995, 100), lang = mu => 1 / (1 + Math.exp(-mu)), bwMu = th => Math.log(th / (1 - th)) - z * wt * th, kTc = w > 0 ? w / 4 * TC2D : 0, Tc = kTc / KB, iHalf = Math.floor(nP / 2);
  let pHalf = log10P(mid); for (let i = 1; i < nP; i++) if ((thMC[i - 1] - 0.5) * (thMC[i] - 0.5) <= 0 && thMC[i] !== thMC[i - 1]) { pHalf = log10P(mus[i - 1] + (0.5 - thMC[i - 1]) / (thMC[i] - thMC[i - 1]) * (mus[i] - mus[i - 1])); break; }
  const qst = Ea + z * w * thU, sub = w > 0 && T < Tc;
  checks.push(chk('Coverage within [0, 1] and monotonic in pressure', thMC.every(v => v >= 0 && v <= 1) && thMC.every((v, i) => i === 0 || v >= thMC[i - 1] - 0.05), `θ from ${F(thMC[0], 3)} to ${F(thMC[nP - 1], 3)}`),
    chk('Particle–hole symmetry point θ = ½', sub || Math.abs(thMC[iHalf] - 0.5) < 0.06, sub ? 'not tested below the lattice-gas critical temperature (two-phase coexistence)' : `θ = ${F(thMC[iHalf], 4)} at ln(KP) = −z w/2k_BT`),
    chk('Langmuir limits recovered', Math.abs(thMC[0] - lang(mus[0] + z * wt * thMC[0])) < 0.02 && thMC[nP - 1] > 0.97, `low-pressure θ = ${F(thMC[0], 3)}, high-pressure θ = ${F(thMC[nP - 1], 4)}`), chk('Acceptance ratio is physical', ur.acc >= 0 && ur.acc <= 1, `${F(100 * ur.acc, 3)} % at the chosen pressure`));
  if (sub) warnings.push(`T = ${F(T, 4)} K is below the lattice-gas critical temperature ${F(Tc, 4)} K: the adsorbate condenses into islands (first-order transition), the isotherm has a near-vertical step and grand-canonical sampling shows hysteresis.`);
  if (w < 0 && -wt > 1.5) warnings.push('Strong lateral repulsion: an ordered c(2×2) half-coverage phase appears as a plateau at θ = ½, which mean-field theory does not describe.');
  warnings.push('The gas-phase reference is a monatomic ideal gas and the adsorbate is immobile and structureless: rotational, vibrational and configurational-entropy terms of real molecules shift the pressure scale by orders of magnitude.');
  return {
    summary: `Grand-canonical Metropolis Monte Carlo of a square lattice gas (${L}×${L} sites, adsorption energy ${F(Ea, 3)} eV, lateral interaction ${F(w * 1e3, 3)} meV, ${F(T, 4)} K). At ${F(P, 3)} Pa the coverage is ${F(thU, 3)} against ${F(lang(muU), 3)} for the non-interacting Langmuir isotherm; half coverage is reached at about 10^${F(pHalf, 3)} Pa. ${w > 0 ? `Attractive interactions steepen the isotherm (lattice-gas T_c = ${F(Tc, 4)} K).` : w < 0 ? 'Repulsive interactions broaden the isotherm.' : 'Without lateral interactions the simulation reproduces the Langmuir isotherm.'}`,
    kpis: [
      kpi('theta', 'Coverage at the chosen pressure', thU, 'ML', `Mean fractional site occupation at ${F(P, 3)} Pa and ${F(T, 4)} K from ${sweeps} Monte Carlo sweeps (standard deviation of the instantaneous coverage ${F(ur.sd, 2)}).`),
      kpi('theta_L', 'Langmuir coverage', lang(muU), 'ML', 'θ = KP/(1 + KP) for independent sites at the same pressure; the difference from the simulated value is the effect of lateral interactions.'),
      kpi('logP_half', 'Pressure of half coverage (log₁₀)', pHalf, 'log₁₀ Pa', `Pressure at which θ = ½ on the simulated isotherm; for independent sites it equals 1/K = 10^${F(-lnK / Math.LN10, 3)} Pa.`),
      kpi('K', 'Langmuir equilibrium constant (log₁₀)', lnK / Math.LN10, 'log₁₀ 1/Pa', `K = exp(E_ads/k_BT)·Λ³/k_BT with the thermal wavelength Λ = ${F(lam * 1e12, 4)} pm of a ${F(mass, 4)} u adsorbate.`),
      kpi('Tc', 'Lattice-gas critical temperature', Tc, 'K', w > 0 ? 'k_BT_c = (w/4)·2/ln(1 + √2) from the exact mapping to the 2-D Ising model; below it the adlayer separates into dense islands and a dilute gas.' : 'No condensation transition exists for zero or repulsive nearest-neighbour interactions; the value is reported as zero.'),
      kpi('qst', 'Isosteric heat estimate', qst, 'eV', 'Mean-field differential heat of adsorption E_ads + z·w·θ at the simulated coverage; it rises with coverage for attractive and falls for repulsive interactions.'),
      kpi('accept', 'Acceptance ratio', ur.acc, '–', 'Fraction of proposed adsorption/desorption moves accepted at the chosen pressure.'),
      kpi('w_kT', 'Lateral interaction in thermal units', wt, 'k_BT', 'w/k_BT; above about 1.76 (attractive) the system is below the critical temperature of the square lattice gas.'),
    ],
    series: [
      { key: 'iso', title: 'Adsorption isotherm', xlabel: 'log₁₀ pressure (Pa)', ylabel: 'coverage θ', lines: [{ name: 'Monte Carlo', x: mus.map(log10P), y: thMC }, { name: 'Langmuir (w = 0)', x: mus.map(log10P), y: mus.map(lang) }, { name: 'Bragg–Williams mean field', x: thG.map(t => log10P(bwMu(t))), y: thG }], explain: 'Coverage against gas pressure. Lateral attraction shifts the isotherm to lower pressure and steepens it; the mean-field curve develops an unphysical S-shaped loop below its own critical temperature, where the true isotherm has a vertical step.' },
      { key: 'dev', title: 'Deviation from the Langmuir isotherm', xlabel: 'log₁₀ pressure (Pa)', ylabel: 'θ_MC − θ_Langmuir', lines: [{ name: 'difference', x: mus.map(log10P), y: mus.map((mu, i) => thMC[i] - lang(mu)) }], explain: 'Excess coverage caused by adsorbate–adsorbate interactions; zero everywhere for an ideal Langmuir adlayer.' },
    ],
    maps: [{ key: 'occ', title: `Site occupation at ${F(P, 3)} Pa`, nx: L, ny: L, data: map, x0: 0, x1: L, y0: 0, y1: L, xlabel: 'x (sites)', ylabel: 'y (sites)', unit: 'occupied', explain: 'Final adsorbate configuration (1 = occupied). Attractive interactions produce compact islands, repulsive interactions a checkerboard-like arrangement, and no interaction a random pattern.' }],
    equations: [
      eqn('governing', 'Lattice-gas Hamiltonian', 'H = −E_ads Σ n_i − w Σ_⟨ij⟩ n_i n_j − μ Σ n_i,  n_i ∈ {0, 1}'), eqn('governing', 'Grand-canonical Metropolis move', 'P(accept) = min{1, exp[∓(w Σ_nn n_j + μ + E_ads − μ°)/k_BT]}  (− adsorption, + desorption)'),
      eqn('constitutive', 'Langmuir isotherm', 'θ = KP/(1 + KP),  K = exp(E_ads/k_BT) Λ³/(k_BT),  Λ = h/√(2π m k_BT)'), eqn('constitutive', 'Bragg–Williams (Fowler–Guggenheim) isotherm', 'KP = θ/(1 − θ) · exp(−z w θ/k_BT),  z = 4'),
      eqn('chemical', 'Adsorption equilibrium', 'A(g) + ∗ ⇌ A∗', 'One adsorbate A per surface site ∗; non-dissociative adsorption.'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', `Operate near 10^${F(pHalf, 3)} Pa (at ${F(T, 4)} K) for half coverage; one decade above or below moves the surface to nearly full or nearly empty for these parameters.`, 'Simulated isotherm.'),
      rec('medium', Math.abs(wt) > 0.3 ? 'Lateral interactions are significant here: a Langmuir fit to measured isotherms would give a coverage-dependent, apparent adsorption energy. Fit a Fowler–Guggenheim form or the lattice-gas model instead.' : 'Lateral interactions are weak at this temperature: the Langmuir model is an adequate description for process calculations.', 'Size of w/k_BT.'),
      rec('medium', 'Take the adsorption energy and nearest-neighbour interaction from DFT (adsorbate pair at neighbouring sites versus isolated) to make this model predictive for a specific surface.', 'DFT-informed lattice-gas modelling.'),
      rec('low', 'Use the isotherm to choose the lowest regeneration temperature or mildest vacuum that empties the surface: avoiding over-heating in sorbent or sensor regeneration saves energy and limits sintering of nanostructured surfaces.', 'Temperature dependence of K.', true),
    ],
    publish: { 'sim.coverage': thU, 'sim.acceptance': ur.acc, 'proc.P': P, 'proc.T_C': T - 273.15, 'sim.q_st_eV': qst },
    numerics: { method: 'Grand-canonical Metropolis Monte Carlo (single-site insertion/removal), periodic square lattice', dof: n, steps: (sweeps + eq) * nP + sweeps + 2 * eq, iterations: n, residual: Math.abs(thMC[iHalf] - 0.5), notes: `Seed ${seed}; pressure sweep upward from the empty surface.` },
  };
}

// ───────────────────────────── mode 3: kinetic Monte Carlo (solid-on-solid growth) ─────────────────────────────
function sosKMC(L, T, Es, En, Ees, nu0, Fdep, covTarget, maxEvents, rand, sampleN = 0) {
  const n = L * L, h = new Int16Array(n), cls = new Int8Array(n).fill(-1), posIn = new Int32Array(n), lists = Array.from({ length: 5 }, () => new Int32Array(n)), cnt = [0, 0, 0, 0, 0], kT = KB * T;
  const rate = [0, 1, 2, 3, 4].map(c => nu0 * Math.exp(-(Es + c * En) / kT)), pES = Math.exp(-Ees / kT), nbr = new Int32Array(4 * n);
  for (let i = 0; i < n; i++) { const xx = i % L, yy = (i - xx) / L; nbr[4 * i] = yy * L + (xx + 1) % L; nbr[4 * i + 1] = yy * L + (xx + L - 1) % L; nbr[4 * i + 2] = ((yy + 1) % L) * L + xx; nbr[4 * i + 3] = ((yy + L - 1) % L) * L + xx; }
  const upd = i => {
    let c = -1; if (h[i] > 0) { c = 0; const hi = h[i]; for (let q = 0; q < 4; q++) if (h[nbr[4 * i + q]] >= hi) c++; }
    const o = cls[i]; if (o === c) return;
    if (o >= 0) { const last = lists[o][--cnt[o]], p = posIn[i]; lists[o][p] = last; posIn[last] = p; }
    if (c >= 0) { posIn[i] = cnt[c]; lists[c][cnt[c]++] = i; } cls[i] = c;
  };
  const touch = i => { upd(i); upd(nbr[4 * i]); upd(nbr[4 * i + 1]); upd(nbr[4 * i + 2]); upd(nbr[4 * i + 3]); };
  let t = 0, dep = 0, ev = 0, hops = 0; const target = Math.round(covTarget * n), Rdep = Fdep * n, samples = [], every = sampleN ? Math.max(1, Math.floor(target / sampleN)) : 0;
  const stats = () => { let s1 = 0, s2 = 0, ab = 0, st = 0; for (let i = 0; i < n; i++) { const v = h[i]; s1 += v; s2 += v * v; ab += v & 1 ? -1 : 1; st += Math.abs(v - h[nbr[4 * i]]) + Math.abs(v - h[nbr[4 * i + 2]]); } const mean = s1 / n; return { cov: dep / n, t, W: Math.sqrt(Math.max(s2 / n - mean * mean, 0)), I: (ab / n) ** 2, step: st / (2 * n), h: sampleN ? Int16Array.from(h) : null }; };
  if (sampleN) samples.push(stats());
  while (dep < target && ev < maxEvents) {
    const R = Rdep + cnt[0] * rate[0] + cnt[1] * rate[1] + cnt[2] * rate[2] + cnt[3] * rate[3] + cnt[4] * rate[4]; let u = rand() * R; t += -Math.log(1 - rand()) / R; ev++;
    if (u < Rdep) { const i = rand.int(n); h[i]++; dep++; touch(i); if (every && dep % every === 0) samples.push(stats()); continue; }
    u -= Rdep; let c = 0; while (c < 4 && u >= cnt[c] * rate[c]) { u -= cnt[c] * rate[c]; c++; } if (!cnt[c]) continue;
    const i = lists[c][Math.min(cnt[c] - 1, Math.floor(u / rate[c]))], j = nbr[4 * i + rand.int(4)];
    if (h[j] < h[i] - 1 && pES < 1 && rand() >= pES) continue;   // step-edge (Ehrlich–Schwoebel) barrier for descending hops
    h[i]--; h[j]++; hops++; touch(i); touch(j);
  }
  if (sampleN && samples[samples.length - 1].cov < dep / n) samples.push(stats());
  return { h, t, dep, ev, hops, samples, capped: dep < target, rate, stats };
}
function countIslands(h, L) {   // 4-connected clusters of columns with height ≥ 1
  const n = L * L, seen = new Uint8Array(n), stack = []; let c = 0;
  for (let s = 0; s < n; s++) { if (seen[s] || h[s] < 1) continue; c++; stack.push(s); seen[s] = 1; while (stack.length) { const i = stack.pop(), xx = i % L, yy = (i - xx) / L; for (const j of [yy * L + (xx + 1) % L, yy * L + (xx + L - 1) % L, ((yy + 1) % L) * L + xx, ((yy + L - 1) % L) * L + xx]) if (!seen[j] && h[j] >= 1) { seen[j] = 1; stack.push(j); } } }
  return c;
}
function runKMC(inp) {
  const L = int(inp, 'kL', 32, 16, 48), T = num(inp, 'T_C', 227, -150, 1500) + 273.15, Es = num(inp, 'Es', 0.75, 0.05, 2.5), En = num(inp, 'En', 0.25, 0, 1), Ees = num(inp, 'Ees', 0, 0, 0.5), Fd = num(inp, 'flux', 1, 1e-4, 1e3), cov = num(inp, 'kCov', 2, 0.2, 5), nu0 = num(inp, 'nu0', 1e13, 1e10, 1e15), seed = int(inp, 'seed', 1, 1, 99999), warnings = [], checks = [];
  const a = num(inp, 'aLat', 2.9e-10, 1e-10, 1e-9), rand = N.rng(seed), res = sosKMC(L, T, Es, En, Ees, nu0, Fd, cov, 6e5, rand, 60), S = res.samples, fin = S[S.length - 1], n = L * L, kT = KB * T, k0 = res.rate[0], DF = k0 / (4 * Fd);
  if (res.capped) warnings.push(`The event budget (600 000 events) was exhausted at ${F(res.dep / n, 3)} ML of the requested ${F(cov, 3)} ML: surface diffusion is very fast relative to deposition (D/F = ${F(DF, 3)}). Lower the temperature, raise the flux or reduce the lattice size to simulate thicker films.`);
  // island density vs D/F at 0.15 ML
  const fMul = [0.1, 0.3, 1, 3, 10, 30], isl = fMul.map((f, q) => { const r = sosKMC(L, T, Es, En, Ees, nu0, Fd * f, 0.15, 8e4, N.rng(seed + 17 * (q + 1))); return { DF: k0 / (4 * Fd * f), N: countIslands(r.h, L) / n, cov: r.dep / n, capped: r.capped }; });
  const okI = isl.filter(q => q.N > 0 && q.cov > 0.1), fitI = okI.length >= 3 ? N.linreg(okI.map(q => Math.log10(q.DF)), okI.map(q => Math.log10(q.N))) : { slope: 0, intercept: 0, r2: 0 }, chi = -fitI.slope;
  if (isl.some(q => q.capped)) warnings.push('Some island-density runs hit their event budget before reaching 0.15 ML and were excluded from the scaling fit.');
  const hmax = N.maxOf(Array.from(res.h)), hmean = res.dep / n, Wf = fin.W, lastML = S.filter(q => q.cov >= Math.max(0, fin.cov - 1)), Wmax = N.maxOf(lastML.map(q => q.W)), Wpoisson = Math.sqrt(fin.cov);
  const modeG = Wmax < 0.62 ? 'layer-by-layer (2-D)' : Wf > 0.75 * Wpoisson ? 'rough / multilayer (3-D, kinetically limited)' : 'intermediate (multilayer with partial interlayer transport)';
  const exposed = (() => { const c = {}; for (let i = 0; i < n; i++) c[res.h[i]] = (c[res.h[i]] || 0) + 1; return Object.keys(c).length; })(), nIsl = countIslands(res.h, L);
  const map = new Float64Array(n); for (let i = 0; i < n; i++) map[i] = res.h[i];
  const nz = Math.min(hmax, 14) + 2, vox = hh => { const d = new Uint8Array(L * L * nz); for (let i = 0; i < n; i++) { d[i] = 1; const top = Math.min(hh[i], nz - 1); for (let zz = 1; zz <= top; zz++) d[i + n * zz] = 2; } return d; };
  const fr = []; const stepF = Math.max(1, Math.floor(S.length / 20)); for (let q = 0; q < S.length; q += stepF) fr.push({ t: +S[q].t.toPrecision(4), data: vox(S[q].h) });
  const D = k0 * a * a / 4, covS = S.map(q => q.cov);
  checks.push(chk('Mass conservation', Math.abs(N.sum(Array.from(res.h)) - res.dep) < 0.5, `Σ heights = ${N.sum(Array.from(res.h))}, atoms deposited = ${res.dep}`), chk('Deposition time consistent with the flux', Math.abs(res.t * Fd / (res.dep / n) - 1) < 0.15, `${F(res.dep / n, 4)} ML in ${F(res.t, 4)} s at ${F(Fd, 3)} ML/s (stochastic arrival times)`),
    chk('Roughness between ideal layer-by-layer and random deposition', Wf <= 1.25 * Wpoisson + 0.1, `W = ${F(Wf, 3)} ML; Poisson limit √θ = ${F(Wpoisson, 3)} ML`), chk('Requested coverage reached', !res.capped, `${F(res.dep / n, 4)} of ${F(cov, 3)} ML in ${res.ev} events`));
  if (En / kT > 14) warnings.push('Lateral bonds are effectively irreversible at this temperature (E_n ≫ k_BT): islands grow fractal-like and cannot compact.');
  warnings.push('Solid-on-solid model on a square lattice: no overhangs, vacancies, strain or reconstruction; rates follow bond counting of the initial site only. Barriers should come from DFT or experiment for a specific surface.');
  return {
    summary: `Solid-on-solid kinetic Monte Carlo growth on a ${L}×${L} lattice at ${F(T, 4)} K: ${F(res.dep / n, 3)} ML deposited at ${F(Fd, 3)} ML/s in ${F(res.t, 3)} s (${res.ev} events, ${res.hops} hops). Adatom hop rate ${F(k0, 3)} s⁻¹ gives D/F = ${F(DF, 3)}. The final rms roughness is ${F(Wf, 3)} ML (random deposition would give ${F(Wpoisson, 3)} ML): growth mode ${modeG}. Island density at 0.15 ML scales as (D/F)^−${F(chi, 2)}.`,
    kpis: [
      kpi('W', 'Final rms roughness', Wf, 'ML', `Standard deviation of the column heights at ${F(fin.cov, 3)} ML (${F(Wf * a * 1e9, 3)} nm for a ${F(a * 1e10, 3)} Å layer). Ideal layer-by-layer growth oscillates between 0 and 0.5 ML; random deposition grows as √θ.`, Wf < 0.75 * Wpoisson),
      kpi('mode', 'Growth mode', modeG, '', `Classified from the largest roughness over the last monolayer (${F(Wmax, 3)} ML) against the layer-by-layer bound 0.5 ML and the random-deposition limit ${F(Wpoisson, 3)} ML.`),
      kpi('DF', 'D/F ratio', DF, '–', `Adatom diffusion coefficient in lattice units, k₀/4 with k₀ = ν₀ exp(−E_s/k_BT) = ${F(k0, 3)} s⁻¹, divided by the deposition flux. Large values (≳ 10⁵) favour few large islands and smooth layers.`),
      kpi('D', 'Adatom surface diffusion coefficient', D, 'm²/s', `D = k₀a²/4 for uncorrelated hops of length a = ${F(a * 1e10, 3)} Å on a square lattice (terrace diffusion of an isolated adatom).`),
      kpi('chi', 'Island-density scaling exponent', chi, '–', `Slope of log N_island against log(D/F) at 0.15 ML from ${okI.length} flux values (R² = ${F(fitI.r2 ?? 0, 3)}). Rate-equation theory gives 1/3 for irreversible attachment (critical nucleus of one atom) and larger values when dimers can dissociate.`, null),
      kpi('N_isl', 'Island density at the reference flux', (isl[2] || isl[0]).N, '1/site', `Number of first-layer islands per lattice site at 0.15 ML and the chosen flux (${F((isl[2] || isl[0]).N / (a * a), 3)} m⁻²).`),
      kpi('layers', 'Exposed layers', exposed, '–', `Number of different height levels present on the final surface (maximum height ${hmax} ML, mean ${F(hmean, 3)} ML); 1–2 for layer-by-layer growth.`),
      kpi('I_ab', 'Anti-Bragg intensity (final)', fin.I, '–', 'Kinematic diffraction intensity at the anti-phase condition, [Σ_h (−1)^h θ_h]², normalised to a flat surface; it oscillates with monolayer period during layer-by-layer growth (RHEED / XRD growth oscillations).'),
      kpi('t', 'Deposition time', res.t, 's', 'Simulated physical time from the residence-time (BKL) algorithm, Δt = −ln(u)/R_total per event.'),
    ],
    series: [
      { key: 'W_cov', title: 'Surface roughness vs coverage', xlabel: 'coverage (ML)', ylabel: 'rms roughness (ML)', lines: [{ name: 'KMC', x: covS, y: S.map(q => q.W) }, { name: 'random deposition √θ', x: covS, y: covS.map(Math.sqrt) }, { name: 'ideal layer-by-layer', x: covS, y: covS.map(c => { const f = c - Math.floor(c); return Math.sqrt(f * (1 - f)); }) }], explain: 'Roughness between the two limiting cases. Curves that follow the oscillating lower bound indicate layer-by-layer growth; curves that approach √θ indicate kinetic roughening because adatoms cannot reach step edges or cross them.' },
      { key: 'I_cov', title: 'Anti-Bragg intensity and step density', xlabel: 'coverage (ML)', ylabel: 'normalised', lines: [{ name: 'anti-Bragg intensity', x: covS, y: S.map(q => q.I) }, { name: 'step density', x: covS, y: S.map(q => q.step) }], explain: 'Growth oscillations: the anti-phase diffraction intensity has maxima at completed layers and the step density has minima there. Damping of the oscillations measures the transition to multilayer growth.' },
      { key: 'N_DF', title: 'Island density vs D/F', xlabel: 'log₁₀ D/F', ylabel: 'log₁₀ islands per site', lines: [{ name: 'KMC at 0.15 ML', x: isl.map(q => Math.log10(q.DF)), y: isl.map(q => Math.log10(Math.max(q.N, 0.5 / n))) }, { name: `fit, slope −${F(chi, 2)}`, x: isl.map(q => Math.log10(q.DF)), y: isl.map(q => fitI.intercept + fitI.slope * Math.log10(q.DF)) }], explain: 'Nucleation scaling N ∝ (D/F)^−χ in the sub-monolayer regime. The exponent identifies the critical nucleus size; the lattice is small, so each point carries counting noise of 1/√(number of islands).' },
    ],
    maps: [{ key: 'height', title: 'Final surface height', nx: L, ny: L, data: map, x0: 0, x1: L * a * 1e9, y0: 0, y1: L * a * 1e9, xlabel: 'x (nm)', ylabel: 'y (nm)', unit: 'ML', explain: 'Height of every lattice column in monolayers. Compact single-level islands indicate two-dimensional growth; stacked mounds indicate hindered interlayer transport.' }],
    equations: [
      eqn('governing', 'Arrhenius hop rate with bond counting', 'k_n = ν₀ exp[−(E_s + n E_n)/k_BT],  n = number of lateral neighbours', `Descending a step is additionally accepted with probability exp(−E_ES/k_BT), E_ES = ${F(Ees, 3)} eV.`),
      eqn('governing', 'Residence-time (BKL) algorithm', 'R = F L² + Σ_n N_n k_n;  event i chosen with probability r_i/R;  Δt = −ln(u)/R'), eqn('governing', 'Master equation sampled', 'dP(σ,t)/dt = Σ_σ′ [W(σ′→σ) P(σ′,t) − W(σ→σ′) P(σ,t)]'),
      eqn('constitutive', 'Nucleation scaling', 'N_island ∝ (D/F)^−χ,  χ = i*/(i* + 2);  D = k₀ a²/4'), eqn('formula', 'Roughness and anti-Bragg intensity', 'W² = ⟨h²⟩ − ⟨h⟩²;   I = [Σ_h (−1)^h θ_h]²'),
    ],
    geometry: { kind: 'voxels', nx: L, ny: L, nz, data: vox(res.h), voxel: a * 1e9, units: 'nm', frames: fr, timeUnit: 's' },
    checks, warnings,
    recommendations: [
      rec('high', modeG.startsWith('layer') ? `Conditions give layer-by-layer growth (D/F = ${F(DF, 3)}): stop deposition at a maximum of the growth oscillation to obtain atomically complete layers.` : `Growth is not layer-by-layer at D/F = ${F(DF, 3)}: raise the substrate temperature (each +50 K multiplies the hop rate by about ${F(Math.exp(Es / KB * (1 / T - 1 / (T + 50))), 3)}) or lower the flux until the roughness curve follows the oscillating lower bound.`, 'Roughness evolution against the limiting cases.'),
      rec('medium', Ees > 0 ? `The step-edge barrier of ${F(Ees, 3)} eV suppresses downward transport (probability ${F(Math.exp(-Ees / kT), 3)} per attempt) and promotes mounds; surfactants or pulsed deposition can restore smooth growth.` : 'No step-edge barrier is included; add an Ehrlich–Schwoebel barrier to study mound formation typical of metal homoepitaxy.', 'Ehrlich–Schwoebel interlayer transport.'),
      rec('medium', `For nanostructure (island or quantum-dot) arrays the island density ${F((isl[2] || isl[0]).N / (a * a), 3)} m⁻² sets the spacing; it follows (D/F)^−${F(chi, 2)}, so a tenfold higher flux raises the density about ${F(10 ** chi, 3)}-fold.`, 'Sub-monolayer nucleation scaling.'),
      rec('low', 'Use the lowest substrate temperature that still gives the required smoothness: thermal budget dominates the energy use of vacuum deposition, and lower temperatures also limit interdiffusion.', 'Temperature dependence of D/F.', true),
    ],
    publish: { 'prop.D_diff': D, 'geom.thickness': hmean * a, 'growth.rate': Fd * a, 'proc.T_C': T - 273.15, 'geom.number_density': (isl[2] || isl[0]).N / (a * a), 'sim.roughness_ML': Wf },
    numerics: { method: 'Rejection-free (Bortz–Kalos–Lebowitz) kinetic Monte Carlo with rate classes by lateral coordination', dof: n, steps: res.ev, iterations: res.hops, residual: 0, notes: `Seed ${seed}; event budget 6×10⁵ (main run) + 8×10⁴ per island-density run.` },
  };
}

// ───────────────────────────── mode 4: electronic structure (model Hamiltonians) ─────────────────────────────
const ACC = 1.42, AG = ACC * Math.sqrt(3);   // carbon–carbon bond and graphene lattice constant, Å
const gA1 = [AG, 0], gA2 = [AG / 2, AG * Math.sqrt(3) / 2], gB1 = [2 * PI / AG, -2 * PI / (AG * Math.sqrt(3))], gB2 = [0, 4 * PI / (AG * Math.sqrt(3))];
const fAbs = (kx, ky) => { const p1 = kx * gA1[0] + ky * gA1[1], p2 = kx * gA2[0] + ky * gA2[1]; return Math.hypot(1 + Math.cos(p1) + Math.cos(p2), Math.sin(p1) + Math.sin(p2)); };
const piBands = (w, g0, s) => [-g0 * w / (1 + s * w), g0 * w / Math.max(1 - s * w, 1e-6)];
function dosHist(vals, lo, hi, nb, norm) { const h = new Float64Array(nb), de = (hi - lo) / nb; for (const v of vals) { const b = Math.floor((v - lo) / de); if (b >= 0 && b < nb) h[b]++; } return { x: Array.from({ length: nb }, (_, i) => lo + (i + 0.5) * de), y: Array.from(h, v => v * norm / de) }; }
function runGraphene(inp, base) {
  const g0 = num(inp, 'gamma0', 2.7, 1, 4), s = num(inp, 'tbS', 0, 0, 0.2), checks = [], warnings = [];
  const Kp = [4 * PI / (3 * AG), 0], Mp = [PI / AG, PI / (Math.sqrt(3) * AG)], path = [[0, 0], Kp, Mp, [0, 0]], names = ['Γ', 'K', 'M', 'Γ'], kx = [], ev = [], ec = [], ticks = [0]; let dist = 0;
  for (let sgm = 0; sgm < 3; sgm++) { const p = path[sgm], q = path[sgm + 1], len = Math.hypot(q[0] - p[0], q[1] - p[1]), np = Math.max(20, Math.round(60 * len / (4 * PI / (3 * AG)))); for (let i = sgm ? 1 : 0; i <= np; i++) { const t = i / np, b = piBands(fAbs(p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])), g0, s); kx.push(dist + t * len); ev.push(b[0]); ec.push(b[1]); } dist += len; ticks.push(dist); }
  const nk = 180, vals = []; for (let i = 0; i < nk; i++) for (let j = 0; j < nk; j++) { const u = (i + 0.5) / nk, v = (j + 0.5) / nk, b = piBands(fAbs(u * gB1[0] + v * gB2[0], u * gB1[1] + v * gB2[1]), g0, s); vals.push(b[0], b[1]); }
  const lo = -3.3 * g0, hi = 3 * g0 / Math.max(1 - 3 * s, 0.3) * 1.05, dos = dosHist(vals, lo, hi, 200, 1 / (nk * nk)), eK = piBands(fAbs(Kp[0], Kp[1]), g0, s), eM = piBands(fAbs(Mp[0], Mp[1]), g0, s), eG = piBands(3, g0, s);
  const vF = 1.5 * g0 * N.Q_E * ACC * 1e-10 / N.HBAR, q = 0.01, slope = piBands(fAbs(Kp[0] + q, Kp[1]), g0, s)[1] / q, vFn = slope * N.Q_E * 1e-10 / N.HBAR, iz = N.argmin(dos.x.map(Math.abs));
  checks.push(chk('Dirac point at K', Math.abs(eK[0]) < 1e-9 && Math.abs(eK[1]) < 1e-9, `E(K) = ${F(eK[1], 2)} eV: valence and conduction bands touch`), chk('Fermi velocity: numerical slope vs 3γ₀a_cc/2ħ', Math.abs(vFn - vF) / vF < 0.02, `${F(vFn, 4)} vs ${F(vF, 4)} m/s`), chk('Density of states integrates to two π states per cell', Math.abs(N.trapz(dos.x, dos.y) - 2) < 0.03, `∫ g(E) dE = ${F(N.trapz(dos.x, dos.y), 4)}`), chk('Van Hove energy at M equals ±γ₀', s > 0 || Math.abs(eM[1] - g0) < 1e-9, `E(M) = ${F(eM[0], 4)} / ${F(eM[1], 4)} eV`));
  warnings.push('Nearest-neighbour π-orbital tight binding: γ₀ ≈ 2.7–3.0 eV is an empirical fit; σ bands, second-neighbour hopping (electron–hole asymmetry beyond the overlap s), spin–orbit coupling and many-body velocity renormalisation are not included.');
  return {
    ...base,
    summary: `Nearest-neighbour tight-binding π bands of graphene (γ₀ = ${F(g0, 3)} eV, overlap s = ${F(s, 3)}): the valence and conduction bands touch at K (zero gap, E(K) = ${F(eK[1], 2)} eV) with Fermi velocity ${F(vF, 4)} m/s; the van Hove singularities lie at ${F(eM[0], 4)} and ${F(eM[1], 4)} eV (M point) and the π band width is ${F(eG[1] - eG[0], 4)} eV.`,
    kpis: [
      kpi('gap', 'Band gap', 2 * Math.abs(eK[1]), 'eV', 'Energy difference between conduction and valence band at the K point; zero for ideal graphene, which is a semimetal with a point-like Fermi surface.', null),
      kpi('vF', 'Fermi velocity', vF, 'm/s', `v_F = 3γ₀a_cc/2ħ with a_cc = ${ACC} Å; the bands are linear (massless Dirac fermions) within about 1 eV of the Dirac point. Measured values are close to 1×10⁶ m/s.`),
      kpi('E_M', 'Van Hove singularity (conduction)', eM[1], 'eV', 'Saddle point of the π* band at M, where the density of states diverges logarithmically; it controls the optical absorption peak near 2γ₀.'),
      kpi('E_Mv', 'Van Hove singularity (valence)', eM[0], 'eV', 'Saddle point of the π band at M; equal and opposite to the conduction value when the overlap integral is zero.'),
      kpi('width', 'π band width', eG[1] - eG[0], 'eV', 'Energy span from the bottom of the π band to the top of the π* band at Γ (6γ₀ for s = 0).'),
      kpi('dos0', 'Density of states at the Dirac point', dos.y[iz], '1/eV per cell', 'Vanishes linearly, g(E) = 2A_c|E|/(πħ²v_F²) per spin; the histogram value is limited by the energy resolution of the k-mesh.'),
      kpi('a', 'Lattice constant', AG, 'Å', 'Graphene lattice constant √3·a_cc used for the Brillouin zone; K = 4π/3a.'),
    ],
    series: [
      { key: 'bands', title: 'π bands along Γ–K–M–Γ', xlabel: `k path (Å⁻¹): ${names.map((nm, i) => `${nm} = ${F(ticks[i], 3)}`).join(', ')}`, ylabel: 'energy (eV)', lines: [{ name: 'π (valence)', x: kx, y: ev }, { name: 'π* (conduction)', x: kx, y: ec }], explain: 'Valence and conduction π bands. They cross linearly at K (the Dirac point, placed at E = 0) and have saddle points at M.' },
      { key: 'dos', title: 'Density of states', xlabel: 'energy (eV)', ylabel: 'states / eV per unit cell (per spin)', lines: [{ name: 'g(E)', x: dos.x, y: dos.y }], explain: 'Histogram over a 180×180 k-mesh: linear V-shape around the Dirac point and van Hove peaks at the M-point energies.' },
    ],
    equations: [eqn('governing', 'Tight-binding secular equation', 'det[H(k) − E S(k)] = 0,   H_AB = −γ₀ f(k),  S_AB = s f(k),  f(k) = 1 + e^(ik·a₁) + e^(ik·a₂)'), eqn('constitutive', 'π-band dispersion', 'E±(k) = ± γ₀ |f(k)| / (1 ∓ s |f(k)|)'), eqn('formula', 'Fermi velocity', 'v_F = 3 γ₀ a_cc / 2ħ'), eqn('formula', 'Material', `${chem.prettyFormula('C')} (graphene), 2 atoms per cell, M = ${F(chem.molarMass('C'), 5)} g/mol`)],
    geometry: (() => { const sh = A.nanosheet('graphene', { nx: 8, ny: 8 }); return { kind: 'atoms', atoms: sh.atoms.map(q2 => ({ el: q2.el, x: q2.x, y: q2.y, z: q2.z })), cell: sh.cell, title: 'Graphene 8×8 sheet' }; })(),
    checks, warnings,
    recommendations: [rec('high', 'Graphene has no band gap: for switching devices open one by quantum confinement (nanoribbons, see the nanotube mode for the analogous zone-folding rule), bilayer gating or chemical functionalisation.', 'Zero gap at the K point.'), rec('medium', `Fit γ₀ to a measured Fermi velocity (v_F = ${F(vF, 3)} m/s for γ₀ = ${F(g0, 3)} eV) or to a DFT band structure before using the model quantitatively.`, 'Empirical nature of the hopping parameter.'), rec('medium', 'Within ±1 eV of the Dirac point the linear Dirac cone is accurate; optical transitions near the M-point saddle need excitonic corrections that this single-particle model lacks.', 'Validity range of the nearest-neighbour model.'), rec('low', 'Model Hamiltonians answer band-topology questions in milliseconds; reserve energy-intensive first-principles runs for parameters the model cannot supply.', 'Method hierarchy.', true)],
    publish: { 'material.formula': 'C', 'prop.bandgap_eV': 0, 'elec.fermi_velocity': vF, 'prop.lattice_a': AG * 1e-10 },
    numerics: { method: 'Analytic 2×2 tight-binding Hamiltonian; Brillouin-zone histogram for the density of states', dof: 2, steps: nk * nk, iterations: 1, residual: Math.abs(eK[1]), notes: `${kx.length} k-points on the path, 180×180 mesh for the DOS.` },
  };
}
function runCNT(inp, base) {
  const g0 = num(inp, 'gamma0', 2.7, 1, 4), s = num(inp, 'tbS', 0, 0, 0.2), checks = [], warnings = []; let n = int(inp, 'tn', 10, 1, 20), m = int(inp, 'tm', 10, 0, 20); if (m > n) [n, m] = [m, n];
  const N2 = n * n + n * m + m * m, dR = gcd(2 * n + m, 2 * m + n), Nc = 2 * N2 / dR, d = AG * Math.sqrt(N2) / PI, Tl = Math.sqrt(3) * AG * Math.sqrt(N2) / dR, t1 = (2 * m + n) / dR, t2 = -(2 * n + m) / dR;
  const K1 = [(-t2 * gB1[0] + t1 * gB2[0]) / Nc, (-t2 * gB1[1] + t1 * gB2[1]) / Nc], K2 = [(m * gB1[0] - n * gB2[0]) / Nc, (m * gB1[1] - n * gB2[1]) / Nc], k2 = Math.hypot(K2[0], K2[1]), ux = K2[0] / k2, uy = K2[1] / k2;
  const nkF = 600, kMax = PI / Tl, minE = new Float64Array(Nc), vals = [], stride = Math.max(1, Math.ceil(Nc * (nkF + 1) / 2.4e5)); let gapC = Infinity;
  for (let q = 0; q < Nc; q++) { let mn = Infinity; for (let i = 0; i <= nkF; i++) { const k = -kMax + 2 * kMax * i / nkF, w = fAbs(k * ux + q * K1[0], k * uy + q * K1[1]), b = piBands(w, g0, s); if (b[1] - b[0] < mn) mn = b[1] - b[0]; if ((i + q) % stride === 0 && i < nkF) vals.push(b[0], b[1]); } minE[q] = mn; if (mn < gapC) gapC = mn; }
  const gap = gapC < 2e-3 ? 0 : gapC, metallic = gap === 0, rule = (n - m) % 3 === 0, Erule = 2 * ACC * g0 / d, order = Array.from({ length: Nc }, (_, q) => q).sort((p, q) => minE[p] - minE[q]).slice(0, Math.min(Nc, 8)), nkP = 120, ks = N.linspace(-kMax, kMax, nkP + 1), lines = [];
  for (const q of order) { const v = [], c = []; for (const k of ks) { const b = piBands(fAbs(k * ux + q * K1[0], k * uy + q * K1[1]), g0, s); v.push(b[0]); c.push(b[1]); } lines.push({ name: `q = ${q} (π)`, x: ks, y: v }, { name: `q = ${q} (π*)`, x: ks, y: c }); }
  const lo = -3.2 * g0, hi = 3.2 * g0 / Math.max(1 - 3 * s, 0.3), dos = dosHist(vals, lo, hi, 300, 2 / vals.length * Nc), sub = [...new Set(Array.from(minE).map(v => +(v).toFixed(3)))].sort((p, q) => p - q), E11 = sub.find(v => v > 2e-3) ?? 0, E22 = sub.find(v => v > E11 + 1e-3) ?? 0;
  const tube = safe(() => A.nanotube({ n, m, length_nm: 3, maxAtoms: 3000 }), null);
  checks.push(chk('Metallic / semiconducting class agrees with the (n − m) mod 3 rule', metallic === rule, `computed gap ${F(gap, 3)} eV; (n − m) mod 3 = ${(n - m) % 3}`), chk('Semiconducting gap against 2a_ccγ₀/d', metallic || Math.abs(gap - Erule) / Erule < 0.25, metallic ? 'not applicable to a metallic tube' : `zone folding ${F(gap, 4)} eV, rule ${F(Erule, 4)} eV (trigonal warping accounts for the difference)`), chk('Sub-band count', order.length > 0 && Nc === 2 * N2 / dR, `${Nc} cutting lines → ${2 * Nc} π sub-bands, ${4 * N2 / dR} atoms per translational cell`));
  if (rule && n !== m) warnings.push('Zone folding predicts this non-armchair tube to be metallic; in reality wall curvature opens a small gap of order 10–100 meV (∝ 1/d²). Only armchair (n,n) tubes are truly metallic.');
  if (d < 7) warnings.push(`For this narrow tube (d = ${F(d / 10, 3)} nm) curvature-induced σ–π rehybridisation changes the bands appreciably; zone folding is reliable above about 1 nm diameter.`);
  warnings.push('Single-particle zone-folded tight binding: optical transition energies of real tubes are shifted by electron–electron and excitonic effects (several tenths of an eV).');
  return {
    ...base,
    summary: `Zone-folded tight-binding bands of the (${n},${m}) carbon nanotube (diameter ${F(d / 10, 4)} nm, ${Nc} cutting lines, translational period ${F(Tl, 4)} Å, γ₀ = ${F(g0, 3)} eV): the tube is ${metallic ? 'metallic (bands cross at the Fermi level)' : `semiconducting with a band gap of ${F(gap, 4)} eV`}; the simple rule 2a_ccγ₀/d gives ${F(Erule, 4)} eV for a semiconducting tube of this diameter. First van Hove transition E₁₁ = ${F(E11, 4)} eV.`,
    kpis: [
      kpi('gap', 'Band gap', gap, 'eV', metallic ? 'A cutting line passes through the K point, so valence and conduction sub-bands touch: metallic within the zone-folding model.' : `Smallest direct separation between π and π* sub-bands; it scales as 1/d (rule value ${F(Erule, 4)} eV).`, null),
      kpi('class', 'Electronic type', metallic ? 'metallic' : 'semiconducting', '', `(n − m) mod 3 = ${(n - m) % 3}: tubes with a remainder of 0 are metallic, 1 or 2 semiconducting. Statistically one third of all chiralities are metallic.`),
      kpi('d', 'Diameter', d * 1e-10, 'm', `d = a√(n² + nm + m²)/π with a = ${F(AG, 4)} Å; chiral angle ${F(Math.atan2(Math.sqrt(3) * m, 2 * n + m) * 180 / PI, 4)}°.`),
      kpi('rule', 'Gap from the 2a_ccγ₀/d rule', Erule, 'eV', 'Linear-band (Dirac-cone) estimate for a semiconducting tube of this diameter; the full zone-folding gap deviates by trigonal warping.'),
      kpi('E11', 'First van Hove transition E₁₁', E11, 'eV', `Separation of the first pair of van Hove singularities${metallic ? ' above the crossing bands (M₁₁)' : ' (S₁₁, equal to the gap)'}; it dominates optical absorption, with an emission wavelength near ${F(E11 > 0 ? 1239.84 / E11 : 0, 4)} nm before excitonic corrections.`),
      kpi('E22', 'Second van Hove transition E₂₂', E22, 'eV', 'Separation of the second pair of van Hove singularities; the ratio E₂₂/E₁₁ is close to 2 for semiconducting tubes in this model.'),
      kpi('Nc', 'Hexagons per translational cell', Nc, '–', `N = 2(n² + nm + m²)/d_R with d_R = ${dR}: the number of one-dimensional cutting lines and of sub-band pairs.`),
      kpi('T', 'Translational period', Tl * 1e-10, 'm', 'Length of the one-dimensional unit cell along the tube axis; armchair tubes have the shortest (one graphene lattice constant).'),
    ],
    series: [
      { key: 'bands', title: `(${n},${m}) sub-bands nearest the Fermi level`, xlabel: 'axial wave vector k (Å⁻¹)', ylabel: 'energy (eV)', lines, explain: 'The eight cutting lines with the smallest direct gap, each giving a π and a π* sub-band over the one-dimensional Brillouin zone −π/T…π/T. Bands meeting at E = 0 make the tube metallic.' },
      { key: 'dos', title: 'Density of states', xlabel: 'energy (eV)', ylabel: 'states / eV per translational cell (per spin)', lines: [{ name: 'g(E)', x: dos.x, y: dos.y }], explain: 'One-dimensional van Hove singularities (1/√E spikes) at every sub-band edge. A finite plateau at E = 0 marks a metallic tube, a gap a semiconducting one.' },
      { key: 'gap_d', title: 'Band gap vs diameter for (n,0) zigzag tubes', xlabel: 'diameter (nm)', ylabel: 'gap (eV)', lines: (() => { const xs = [], ys = [], yr = []; for (let nn = 5; nn <= 25; nn++) { if (nn % 3 === 0) continue; const dd = AG * nn / PI; let mn = Infinity; for (let q = 0; q < 2 * nn; q++) { const w = fAbs(q * 2 * PI / (nn * AG), 0), b = piBands(w, g0, s); mn = Math.min(mn, b[1] - b[0]); } xs.push(dd / 10); ys.push(mn); yr.push(2 * ACC * g0 / dd); } return [{ name: 'zone folding', x: xs, y: ys }, { name: '2a_ccγ₀/d', x: xs, y: yr }]; })(), explain: 'Semiconducting zigzag tubes: the gap falls as 1/d. The two families (n mod 3 = 1 or 2) scatter on either side of the linear-band rule because of trigonal warping.' },
    ],
    equations: [eqn('governing', 'Zone folding (periodic boundary condition around the circumference)', 'k·C_h = 2πq,  q = 0 … N − 1;   E_q(k) = ± γ₀ |f(k K₂/|K₂| + q K₁)|'), eqn('constitutive', 'Metallicity rule and gap', '(n − m) mod 3 = 0 ⇒ metallic;  otherwise E_g ≈ 2 a_cc γ₀ / d'), eqn('formula', 'Geometry', `C_h = n a₁ + m a₂,  d = a√(n² + nm + m²)/π = ${F(d, 4)} Å,  N = 2(n² + nm + m²)/d_R = ${Nc}`), eqn('formula', 'Material', `${chem.prettyFormula('C')}, ${4 * N2 / dR} atoms per translational cell`)],
    geometry: tube ? { kind: 'atoms', atoms: tube.atoms.map(q => ({ el: q.el, x: q.x, y: q.y, z: q.z })), cell: tube.cell, title: tube.title } : undefined,
    checks, warnings,
    recommendations: [rec('high', metallic ? `(${n},${m}) is metallic: suitable for interconnects and transparent electrodes; for transistor channels choose a chirality with (n − m) mod 3 ≠ 0.` : `(${n},${m}) is semiconducting with E_g ≈ ${F(gap, 3)} eV: suitable for transistor channels and near-infrared emitters (E₁₁ ≈ ${F(E11 > 0 ? 1239.84 / E11 : 0, 4)} nm before excitonic shifts).`, 'Zone-folding classification.'), rec('medium', `Select diameter to set the gap: E_g ≈ ${F(2 * ACC * g0 / 10, 3)} eV·nm / d. Tubes near 1 nm give gaps near 0.8 eV; wider tubes approach the semimetallic graphene limit.`, 'E_g ≈ 2a_ccγ₀/d.'), rec('medium', 'Assign chirality experimentally from the pair (E₁₁, E₂₂) by photoluminescence-excitation or resonant Raman spectroscopy, using the ratios computed here as the starting pattern.', 'Van Hove transition energies.'), rec('low', 'As-grown material is a mixture of metallic and semiconducting tubes; chirality-selective growth or low-energy sorting (gel chromatography, aqueous two-phase extraction) avoids wasting two thirds or one third of the product.', 'One third of chiralities are metallic.', true)],
    publish: { 'material.formula': 'C', 'prop.bandgap_eV': gap, 'geom.diameter': d * 1e-10, 'elec.E11_eV': E11 },
    exports: tube ? [{ name: `cnt_${n}_${m}.xyz`, mime: 'chemical/x-xyz', content: A.toXYZ(tube), description: 'Nanotube coordinates (Å).' }, { name: `cnt_${n}_${m}.cif`, mime: 'chemical/x-cif', content: A.toCIF(tube), description: 'Nanotube in a periodic box (CIF).' }] : [],
    numerics: { method: 'Zone-folded nearest-neighbour tight binding; minimum-gap search on a uniform axial k-grid', dof: 2 * Nc, steps: Nc * (nkF + 1), iterations: 1, residual: metallic ? gapC : 0, notes: `${nkF + 1} k-points per cutting line (grid contains the K-point projection for metallic tubes); gap threshold 2 meV.` },
  };
}
function runSSH(inp, base) {
  const t1 = num(inp, 't1', 2.4, 0.05, 5), t2 = num(inp, 't2', 3.0, 0.05, 5), a = num(inp, 'aCh', 2.8, 1, 10), checks = [], warnings = [], nk = 161, ks = N.linspace(-PI / a, PI / a, nk), Ek = ks.map(k => Math.sqrt(t1 * t1 + t2 * t2 + 2 * t1 * t2 * Math.cos(k * a)));
  const gap = 2 * Math.abs(t1 - t2), width = 2 * (t1 + t2), Nsite = 60, H = N.zeros(Nsite, Nsite); for (let i = 0; i < Nsite - 1; i++) { const t = i % 2 === 0 ? t1 : t2; H[i][i + 1] = H[i + 1][i] = -t; }
  const eg = N.eigSym(H), mid = eg.values.filter(v => Math.abs(v) < 0.25 * Math.max(gap / 2, 1e-9)).length, topo = t2 > t1, iEdge = N.argmin(eg.values.map(Math.abs)), psi = eg.vectors[iEdge].map(v => v * v);
  const vals = []; for (let i = 0; i < 4000; i++) { const k = -PI / a + 2 * PI / a * (i + 0.5) / 4000, e = Math.sqrt(t1 * t1 + t2 * t2 + 2 * t1 * t2 * Math.cos(k * a)); vals.push(e, -e); } const dos = dosHist(vals, -1.05 * (t1 + t2), 1.05 * (t1 + t2), 200, 1 / 4000);
  const mEff = gap > 1e-6 ? 2 * HB2M * 100 * (gap / 2) / (t1 * t2 * a * a) : 0;   // m*/m₀ = ħ²Δ/(t₁t₂a²) with Δ = |t₁ − t₂|, in units ħ²/2m₀ = 3.81 eV Å²
  const bulkMin = N.minOf(eg.values.filter(v => v > 0.25 * Math.max(gap / 2, 1e-9)).concat([t1 + t2]));
  checks.push(chk('Finite-chain bulk gap matches 2|t₁ − t₂|', gap < 1e-6 || Math.abs(2 * bulkMin - gap) / gap < 0.15, `finite chain ${F(2 * bulkMin, 4)} eV, analytic ${F(gap, 4)} eV`), chk('Edge states present only in the topological phase', (mid === 2) === (topo && gap > 0.05), `${mid} mid-gap states; t₂ ${topo ? '>' : '≤'} t₁`), chk('Spectrum is particle–hole symmetric', Math.abs(eg.values[0] + eg.values[Nsite - 1]) < 1e-8, `E_min + E_max = ${F(eg.values[0] + eg.values[Nsite - 1], 2)} eV`));
  warnings.push('One-orbital-per-site chain with nearest-neighbour hopping; the parameters are illustrative of trans-polyacetylene (t ≈ 2.5 eV, dimerisation gap ≈ 1.4 eV). Lattice relaxation (solitons, polarons) and electron correlation are not included.');
  return {
    ...base,
    summary: `${Math.abs(t1 - t2) < 1e-9 ? 'Uniform' : 'Dimerised (Su–Schrieffer–Heeger)'} chain with hoppings t₁ = ${F(t1, 3)} eV and t₂ = ${F(t2, 3)} eV: band gap ${F(gap, 4)} eV, total band width ${F(width, 4)} eV. A ${Nsite}-site open chain starting with a t₁ bond has ${mid} mid-gap edge state${mid === 1 ? '' : 's'} (${topo ? 'topological phase, winding number 1' : 'trivial phase, winding number 0'}).`,
    kpis: [
      kpi('gap', 'Band gap', gap, 'eV', 'E_g = 2|t₁ − t₂| at the zone boundary: bond-length alternation (Peierls distortion) opens a gap in the half-filled chain; it vanishes for equal hoppings (metallic monatomic chain).'),
      kpi('width', 'Total band width', width, 'eV', 'Energy span 2(t₁ + t₂) of the two bands together; 4t for the uniform chain.'),
      kpi('edge', 'Mid-gap edge states', mid, '–', `Number of states of the finite ${Nsite}-site chain inside the central quarter of the gap. Two (one per end) appear when the weaker bond terminates the chain.`),
      kpi('winding', 'Winding number', topo && gap > 0 ? 1 : 0, '–', 'Topological invariant of the bulk Hamiltonian h(k) = t₁ + t₂e^(ika): it counts how often h(k) encircles the origin and equals the number of protected zero modes per end.'),
      kpi('m_eff', 'Band-edge effective mass', mEff, 'm₀', 'm* = ħ²|t₁ − t₂|/(t₁t₂a²) from the curvature of the band at the gap edge; lighter carriers for weaker dimerisation.'),
      kpi('loc', 'Edge-state weight on the first ten sites', N.sum(psi.slice(0, 10)) + 0, '–', 'Probability of the state closest to E = 0 on the first ten sites; about 0.5 for a topological edge-state pair (half on each end), small for an extended bulk state.'),
    ],
    series: [
      { key: 'bands', title: 'Band structure', xlabel: 'k (Å⁻¹)', ylabel: 'energy (eV)', lines: [{ name: 'valence', x: ks, y: Ek.map(v => -v) }, { name: 'conduction', x: ks, y: Ek }], explain: 'Two bands of the dimerised chain. The gap opens at the zone boundary k = ±π/a and equals twice the difference of the hopping integrals.' },
      { key: 'spec', title: `Spectrum of the open ${Nsite}-site chain`, xlabel: 'state index', ylabel: 'energy (eV)', lines: [{ name: 'eigenvalues', x: eg.values.map((_, i) => i + 1), y: eg.values }], explain: 'Eigenvalues from direct diagonalisation. Two states pinned at zero energy in the middle of the gap are the topological edge states.' },
      { key: 'edge', title: 'Probability density of the state nearest E = 0', xlabel: 'site', ylabel: '|ψ|²', lines: [{ name: '|ψ|²', x: psi.map((_, i) => i + 1), y: psi }], explain: 'An edge state decays exponentially into the chain with localisation length a/ln(t₂/t₁) and lives on one sublattice only; a bulk state is spread over the whole chain.' },
      { key: 'dos', title: 'Density of states', xlabel: 'energy (eV)', ylabel: 'states / eV per cell (per spin)', lines: [{ name: 'g(E)', x: dos.x, y: dos.y }], explain: 'One-dimensional van Hove singularities at all four band edges and the dimerisation gap between them.' },
    ],
    equations: [eqn('governing', 'SSH Hamiltonian', 'H = −Σ_n [t₁ c†_{A,n} c_{B,n} + t₂ c†_{B,n} c_{A,n+1} + h.c.]'), eqn('constitutive', 'Dispersion and gap', 'E±(k) = ± √(t₁² + t₂² + 2 t₁ t₂ cos ka),   E_g = 2|t₁ − t₂|'), eqn('formula', 'Winding number', 'ν = (1/2πi) ∮ dk ∂_k ln h(k),  h(k) = t₁ + t₂ e^(ika);  ν = 1 for t₂ > t₁')],
    checks, warnings,
    recommendations: [rec('high', gap > 0 ? `The dimerisation gap of ${F(gap, 3)} eV makes the chain a semiconductor; doping creates mid-gap soliton or polaron states, which is the conduction mechanism of conjugated polymers.` : 'With equal hoppings the chain is metallic, but a real half-filled chain is unstable to dimerisation (Peierls theorem).', 'Dispersion relation.'), rec('medium', 'Terminate the chain on the weak bond to obtain protected end states; the same design rule gives edge modes in photonic, acoustic and mechanical SSH metamaterial lattices.', 'Bulk–boundary correspondence.'), rec('medium', 'Extract t₁ and t₂ for a real polymer or atomic chain from DFT bands (band width and gap give both parameters) before quantitative use.', 'Two-parameter model.'), rec('low', 'A two-parameter chain model reproduces the gap physics of conjugated polymers at negligible compute cost; use it to screen substitution patterns before running large quantum-chemistry jobs.', 'Method hierarchy.', true)],
    publish: { 'prop.bandgap_eV': gap, 'elec.bandwidth_eV': width },
    numerics: { method: 'Analytic Bloch dispersion + Jacobi diagonalisation of the open-chain Hamiltonian', dof: Nsite, steps: nk, iterations: 1, residual: Math.abs(eg.values[0] + eg.values[Nsite - 1]), notes: `${Nsite} sites, open boundaries, first bond t₁.` },
  };
}
// Quantum-dot materials: bulk gap (eV, 300 K), electron and hole effective masses (m₀), relative permittivity used in the Brus Coulomb term.
const QD = {
  CdSe: { Eg: 1.74, me: 0.13, mh: 0.45, eps: 10.6, f: 'CdSe' }, CdS: { Eg: 2.42, me: 0.21, mh: 0.80, eps: 5.7, f: 'CdS' }, ZnO: { Eg: 3.37, me: 0.24, mh: 0.59, eps: 8.5, f: 'ZnO' },
  Si: { Eg: 1.12, me: 0.26, mh: 0.39, eps: 11.7, f: 'Si', indirect: true }, PbS: { Eg: 0.41, me: 0.085, mh: 0.085, eps: 17.2, f: 'PbS' }, InP: { Eg: 1.34, me: 0.08, mh: 0.60, eps: 12.4, f: 'InP' },
};
// lowest s-state of a finite spherical well (radius R nm, depth V eV, masses in m₀; BenDaniel–Duke matching)
function sphereWell(R, V, mIn, mOut) {
  const Einf = HB2M * PI * PI / (mIn * R * R), g = E => { const k = Math.sqrt(E * mIn / HB2M), kap = Math.sqrt((V - E) * mOut / HB2M); return (k / Math.tan(k * R) - 1 / R) / mIn + (kap + 1 / R) / mOut; };
  const hi = Math.min(V, Einf) * (1 - 1e-9), lo = hi * 1e-6; let a = lo, fa = g(a);
  if (!(hi > lo)) return null; for (let i = 1; i <= 200; i++) { const b = lo + (hi - lo) * i / 200, fb = g(b); if (Number.isFinite(fa) && Number.isFinite(fb) && fa * fb <= 0) return safe(() => N.brent(g, a, b, 1e-12), null); a = b; fa = fb; }
  return null;
}
function runBrus(inp, base) {
  const mk = pick(inp.qdMat, Object.keys(QD), 'CdSe'), q = QD[mk], d = num(inp, 'qdD', 4e-9, 1e-9, 2e-8) * 1e9, V0 = num(inp, 'qdV0', 4, 0.5, 10), checks = [], warnings = [];
  const brus = D => { const R = D / 2, conf = HB2M * PI * PI / (R * R) * (1 / q.me + 1 / q.mh), coul = 1.786 * 1.439964 / (q.eps * R); return { conf, coul, E: q.Eg + conf - coul }; };
  const fw = D => { const R = D / 2, e = sphereWell(R, V0, q.me, 1), h = sphereWell(R, V0, q.mh, 1); if (e === null || h === null) return null; return q.Eg + e + h - 1.786 * 1.439964 / (q.eps * R); };
  const b = brus(d), Efw = fw(d), lam = E => (E > 0 ? 1239.841984 / E : 0), mu = 1 / (1 / q.me + 1 / q.mh), aB = 0.0529177 * q.eps / mu, Ds = N.linspace(1.2, 12, 80), EB = Ds.map(D => brus(D).E), EF = Ds.map(D => fw(D) ?? brus(D).E);
  // CdSe: empirical sizing curve D(λ) of first-exciton absorption (polynomial calibration for colloidal dots, valid for about 1.5–8 nm)
  const yu = l => 1.6122e-9 * l ** 4 - 2.6575e-6 * l ** 3 + 1.6242e-3 * l * l - 0.4277 * l + 41.57; let lamEmp = null;
  if (mk === 'CdSe' && d >= 1.5 && d <= 8) lamEmp = safe(() => N.brent(l => yu(l) - d, 400, 690, 1e-9), null);
  const regime = d / 2 < aB ? 'strong confinement (R < a_B)' : d / 2 < 3 * aB ? 'intermediate confinement' : 'weak confinement (R ≫ a_B)';
  checks.push(chk('Confinement raises the gap above the bulk value', b.E > q.Eg - b.coul, `E = ${F(b.E, 4)} eV, bulk ${q.Eg} eV`), chk('Finite-barrier gap lies between bulk and the infinite-barrier (Brus) value', Efw === null || (Efw <= b.E + 1e-9 && Efw >= q.Eg - b.coul - 1e-9), Efw === null ? 'no bound state for this barrier' : `${F(Efw, 4)} eV vs ${F(b.E, 4)} eV`), chk('Effective-mass approximation applicable (dot larger than about two lattice constants)', d > 1.5, `d = ${F(d, 3)} nm`), chk('Coulomb term treated as a perturbation', b.coul < b.conf || d / 2 > aB, `confinement ${F(b.conf, 3)} eV, Coulomb ${F(b.coul, 3)} eV`));
  if (lamEmp) checks.push(chk('Finite-barrier estimate within 0.3 eV of the empirical CdSe sizing curve', Efw !== null && Math.abs(Efw - 1239.841984 / lamEmp) < 0.3, `model ${F(Efw ?? 0, 4)} eV, empirical first-exciton energy ${F(1239.841984 / lamEmp, 4)} eV`));
  warnings.push(`The Brus equation assumes parabolic bands and infinite barriers: it overestimates the gap of small dots (below about ${F(2 * aB, 3)} nm here) by several tenths of an eV. The finite-barrier value (barrier ${F(V0, 3)} eV, free-electron mass outside) is the better estimate; neither includes band non-parabolicity or surface states.`);
  if (q.indirect) warnings.push('Silicon has an indirect gap with anisotropic valleys: the single-band effective-mass estimate is qualitative only, and emission is phonon-assisted and weak.');
  if (mk === 'PbS') warnings.push('PbS has strongly non-parabolic bands; the parabolic Brus formula overestimates confinement markedly for dots below about 6 nm (hyperbolic-band models are used in practice).');
  const Euse = Efw ?? b.E;
  return {
    ...base,
    summary: `${chem.prettyFormula(q.f)} quantum dot of ${F(d, 3)} nm diameter (${regime}, exciton Bohr radius ${F(aB, 3)} nm): the Brus equation gives a gap of ${F(b.E, 4)} eV (bulk ${q.Eg} eV + confinement ${F(b.conf, 3)} eV − Coulomb ${F(b.coul, 3)} eV), i.e. emission near ${F(lam(b.E), 4)} nm. With a finite ${F(V0, 3)} eV barrier the estimate is ${Efw === null ? 'not bound' : `${F(Efw, 4)} eV (${F(lam(Efw), 4)} nm)`}${lamEmp ? `; the empirical sizing curve for colloidal CdSe gives ${F(1239.841984 / lamEmp, 4)} eV (${F(lamEmp, 4)} nm)` : ''}.`,
    kpis: [
      kpi('Eg_brus', 'Band gap (Brus, infinite barrier)', b.E, 'eV', `E = E_g + ħ²π²/(2R²)(1/m_e + 1/m_h) − 1.786e²/(4πε₀ε_rR) with R = ${F(d / 2, 3)} nm, m_e = ${q.me}, m_h = ${q.mh}, ε_r = ${q.eps}. An upper estimate: hard walls over-confine the carriers.`),
      kpi('Eg_fw', 'Band gap (finite barrier)', Efw ?? b.E, 'eV', Efw === null ? 'No bound s-state exists for this barrier height; the Brus value is repeated.' : `Lowest electron and hole s-levels of a spherical well of depth ${F(V0, 3)} eV with mass-discontinuity (BenDaniel–Duke) matching, plus the same Coulomb term. Wave-function leakage lowers the gap by ${F(b.E - Efw, 3)} eV relative to Brus.`),
      kpi('lambda', 'Emission wavelength (finite barrier)', lam(Euse), 'nm', `λ = hc/E for the finite-barrier gap; the Brus value corresponds to ${F(lam(b.E), 4)} nm. Real emission is further red-shifted by the Stokes shift (tens of meV).`),
      kpi('conf', 'Confinement energy', b.conf, 'eV', 'Kinetic-energy cost of squeezing electron and hole into the dot, ∝ 1/R²; it dominates below the exciton Bohr radius.'),
      kpi('coul', 'Electron–hole Coulomb energy', b.coul, 'eV', 'First-order electron–hole attraction 1.786e²/(4πε₀ε_rR), ∝ 1/R; it partly cancels the confinement energy.'),
      kpi('aB', 'Exciton Bohr radius', aB * 1e-9, 'm', `a_B = ε_r a₀ m₀/μ with reduced mass μ = ${F(mu, 3)} m₀; confinement effects are strong when the dot radius (${F(d / 2, 3)} nm) is below it.`),
      kpi('shift', 'Blue shift from bulk', Euse - q.Eg, 'eV', `Increase of the optical gap relative to bulk ${chem.prettyFormula(q.f)} (${q.Eg} eV, ${F(lam(q.Eg), 4)} nm).`),
      kpi('emp', 'Empirical first-exciton energy (CdSe only)', lamEmp ? 1239.841984 / lamEmp : 'not available', lamEmp ? 'eV' : '', lamEmp ? 'From the widely used polynomial sizing curve relating the first absorption peak of colloidal CdSe dots to their diameter; an experimental calibration to judge the two model estimates.' : 'An experimental sizing curve is included only for CdSe between 1.5 and 8 nm.'),
    ],
    series: [
      { key: 'Eg_d', title: 'Band gap vs dot diameter', xlabel: 'diameter (nm)', ylabel: 'energy gap (eV)', lines: [{ name: 'Brus (infinite barrier)', x: Ds, y: EB.map(v => Math.min(v, q.Eg + 6)) }, { name: `finite barrier ${F(V0, 3)} eV`, x: Ds, y: EF.map(v => Math.min(v, q.Eg + 6)) }, { name: 'bulk', x: [Ds[0], Ds[Ds.length - 1]], y: [q.Eg, q.Eg] }, ...(mk === 'CdSe' ? [{ name: 'empirical sizing curve (CdSe)', x: N.linspace(450, 660, 40).map(yu), y: N.linspace(450, 660, 40).map(l => 1239.841984 / l) }] : [])], explain: 'Size dependence of the optical gap (Brus curve clipped at 6 eV above bulk). The gap rises steeply below the exciton Bohr radius; the finite-barrier curve flattens for small dots because the wave function leaks out.' },
      { key: 'lam_d', title: 'Emission wavelength vs dot diameter', xlabel: 'diameter (nm)', ylabel: 'wavelength (nm)', lines: [{ name: 'finite barrier', x: Ds, y: EF.map(lam) }, { name: 'Brus', x: Ds, y: EB.map(lam) }], explain: 'Colour tuning by size: read off the diameter needed for a target emission wavelength, then refine against an experimental sizing curve for the material.' },
    ],
    equations: [eqn('formula', 'Material', `${chem.prettyFormula(q.f)}, M = ${F(chem.molarMass(q.f), 5)} g/mol, E_g(bulk) = ${q.Eg} eV`), eqn('constitutive', 'Brus equation', 'E(R) = E_g + (ħ²π²/2R²)(1/m_e* + 1/m_h*) − 1.786 e²/(4πε₀ε_r R)'), eqn('governing', 'Finite spherical well (s-state)', '(1/m_in)[k cot(kR) − 1/R] = −(1/m_out)[κ + 1/R],  k = √(2m_inE)/ħ,  κ = √(2m_out(V₀ − E))/ħ'), eqn('formula', 'Exciton Bohr radius', 'a_B = 4πε₀ε_rħ²/(μe²),  1/μ = 1/m_e* + 1/m_h*')],
    checks, warnings,
    recommendations: [rec('high', `For emission at a target wavelength choose the diameter from the finite-barrier curve and then calibrate against measured absorption: a ±5 % size dispersion at ${F(d, 3)} nm broadens the emission by about ${F(Math.abs(((fw(d * 1.05) ?? brus(d * 1.05).E) - (fw(d * 0.95) ?? brus(d * 0.95).E))) * 1000, 3)} meV.`, 'Slope dE/dD of the size-dependence curve.'), rec('medium', 'Grow a wider-gap shell (e.g. ZnS on CdSe) to passivate surface traps and raise the confining barrier; the barrier-height input shows how the gap responds.', 'Finite-barrier model.'), rec('medium', `Dots smaller than about ${F(2 * aB, 3)} nm are in strong confinement, where the effective-mass approximation degrades; use atomistic tight-binding or DFT for quantitative gaps there.`, 'Exciton Bohr radius.'), rec('low', mk === 'CdSe' || mk === 'CdS' || mk === 'PbS' ? 'Cadmium and lead chalcogenides are toxic and restricted in many products; InP- or ZnO-based dots cover much of the same spectral range with lower hazard and should be assessed first.' : 'This material avoids cadmium and lead; keep the synthesis solvent and ligand choice equally benign to retain the environmental advantage.', 'Hazard and regulatory profile of quantum-dot materials.', true)],
    publish: { 'material.formula': q.f, 'prop.bandgap_eV': Euse, 'geom.diameter': d * 1e-9, 'elec.emission_nm': lam(Euse), 'prop.eps_r': q.eps },
    numerics: { method: 'Closed-form Brus equation; Brent root of the finite spherical-well matching condition', dof: 2, steps: Ds.length, iterations: 1, residual: 0, notes: 'Single parabolic band per carrier, barrier mass m₀.' },
  };
}
function runWell(inp, base) {
  const W = num(inp, 'wW', 6e-9, 5e-10, 3e-8) * 1e9, V0 = num(inp, 'wV0', 0.3, 0.01, 5), mw = num(inp, 'wMw', 0.067, 0.01, 2), mb = num(inp, 'wMb', 0.092, 0.01, 2), checks = [], warnings = [];
  const kap0 = Math.sqrt(V0 * mb / HB2M), pad0 = N.clamp(Math.max(6 / kap0, 0.5 * W), 0.5, 60), nTot = 101, nW = N.clamp(Math.round(nTot * W / (W + 2 * pad0)), 10, 81), dx = W / nW, nP = Math.floor((nTot - nW) / 2), ng = nW + 2 * nP, pad = nP * dx, Ltot = W + 2 * pad;
  // cell-centred grid: the interfaces lie exactly midway between two grid points; 1/m at a midpoint is the mean of its two neighbours
  const xs = Array.from({ length: ng }, (_, i) => -Ltot / 2 + (i + 0.5) * dx), inW = x => Math.abs(x) < W / 2, H = N.zeros(ng, ng), im = xs.map(x => 1 / (inW(x) ? mw : mb)), imAt = (i, j) => 0.5 * (im[i] + im[N.clamp(j, 0, ng - 1)]);
  for (let i = 0; i < ng; i++) { const mp = imAt(i, i + 1), mm = imAt(i, i - 1); H[i][i] = HB2M / (dx * dx) * (mp + mm) + (inW(xs[i]) ? 0 : V0); if (i + 1 < ng) H[i][i + 1] = H[i + 1][i] = -HB2M / (dx * dx) * mp; }
  const eg = N.eigSym(H), bound = eg.values.filter(v => v < V0), nB = bound.length, Einf = k => HB2M * PI * PI * k * k / (mw * W * W);
  // analytic ground state (even parity): (k/m_w) tan(kW/2) = κ/m_b
  const gE = E => { const k = Math.sqrt(E * mw / HB2M), kp = Math.sqrt((V0 - E) * mb / HB2M); return k / mw * Math.tan(k * W / 2) - kp / mb; }, top = Math.min(V0, Einf(1)) * (1 - 1e-9), E1a = safe(() => N.brent(gE, 1e-12, top, 1e-13), null);
  const E1 = eg.values[0], nExp = 1 + Math.floor(Math.sqrt(V0 * mw / HB2M) * W / PI + 1e-9), lines = [{ name: 'potential V(x)', x: xs, y: xs.map(x => (inW(x) ? 0 : V0)) }], sc = 0.6 * (nB > 1 ? (bound[1] - bound[0]) : V0 / 2);
  for (let k = 0; k < Math.min(nB, 4); k++) { const v = eg.vectors[k], mx = N.maxOf(v.map(Math.abs)) || 1, sgn = v[N.argmax(v.map(Math.abs))] > 0 ? 1 : -1; lines.push({ name: `ψ${k + 1} (E = ${F(bound[k], 3)} eV)`, x: xs, y: v.map(u => bound[k] + sgn * sc * u / mx) }); }
  const pOut = (() => { const v = eg.vectors[0]; let o = 0, t = 0; v.forEach((u, i) => { t += u * u; if (!inW(xs[i])) o += u * u; }); return o / t; })(), E12 = nB > 1 ? bound[1] - bound[0] : 0;
  checks.push(chk('Finite-difference ground state against the analytic matching condition', E1a === null || Math.abs(E1 - E1a) / E1a < 0.03, E1a === null ? 'analytic root not bracketed' : `numerical ${F(E1, 5)} eV, analytic ${F(E1a, 5)} eV (grid spacing ${F(dx, 3)} nm)`), chk('Bound-state count agrees with 1 + ⌊√(2m_wV₀)W/πħ⌋', mw !== mb ? Math.abs(nB - nExp) <= 1 : nB === nExp, `${nB} found, ${nExp} expected (exact for equal masses)`), chk('Levels below the infinite-well values', nB === 0 || E1 < Einf(1), `E₁ = ${F(E1, 4)} eV, infinite well ${F(Einf(1), 4)} eV`), chk('Wave function vanishes at the box edge', Math.abs(eg.vectors[0][0]) < 1e-3 * N.maxOf(eg.vectors[0].map(Math.abs)), `edge amplitude ${F(Math.abs(eg.vectors[0][0]) / N.maxOf(eg.vectors[0].map(Math.abs)), 2)} of the maximum`));
  if (nB === 0) warnings.push('No level lies below the barrier on this grid; a symmetric one-dimensional well always has one bound state, but it is too weakly bound to resolve here.');
  warnings.push('Single-band envelope-function model with position-dependent effective mass; band non-parabolicity, strain, electric fields and many-body effects are not included. Defaults correspond to a GaAs well in Al₀.₃Ga₀.₇As barriers (conduction band).');
  return {
    ...base,
    summary: `Finite quantum well of width ${F(W, 3)} nm and depth ${F(V0, 3)} eV (masses ${F(mw, 3)} / ${F(mb, 3)} m₀): ${nB} bound state${nB === 1 ? '' : 's'}, ground level E₁ = ${F(E1, 4)} eV (infinite-well value ${F(Einf(1), 4)} eV${E1a !== null ? `, analytic finite-well value ${F(E1a, 4)} eV` : ''}); ${F(100 * pOut, 3)} % of the ground-state probability lies in the barriers.${nB > 1 ? ` The E₁→E₂ intersubband transition is ${F(E12, 4)} eV (${F(1.239841984 / E12, 4)} µm).` : ''}`,
    kpis: [
      kpi('E1', 'Ground-state energy', E1, 'eV', 'Lowest eigenvalue of the finite-difference Hamiltonian, measured from the well bottom; it is the confinement shift of the band edge.'),
      kpi('E1_inf', 'Infinite-well ground state', Einf(1), 'eV', 'ħ²π²/(2m_wW²) for hard walls; the finite barrier lowers the level because the wave function penetrates the barriers.'),
      kpi('E1_an', 'Analytic finite-well ground state', E1a ?? E1, 'eV', 'Root of the even-parity matching condition (k/m_w) tan(kW/2) = κ/m_b; the benchmark for the numerical eigenvalue.'),
      kpi('n_bound', 'Bound states', nB, '–', 'Number of eigenvalues below the barrier height; states above are box-quantised continuum states of the finite grid.'),
      kpi('E12', 'E₁ → E₂ transition', E12, 'eV', nB > 1 ? `Intersubband spacing, corresponding to a wavelength of ${F(1.239841984 / E12, 4)} µm — the design quantity of quantum-well infrared photodetectors and cascade lasers.` : 'Only one bound state exists, so there is no bound-to-bound transition.'),
      kpi('leak', 'Probability in the barriers', pOut, '–', 'Fraction of |ψ₁|² outside the well; it controls tunnel coupling to neighbouring wells and the sensitivity to interface roughness.'),
      kpi('decay', 'Barrier decay length of ψ₁', nB ? 1e-9 / Math.sqrt(Math.max(V0 - E1, 1e-9) * mb / HB2M) : 0, 'm', '1/κ with κ = √(2m_b(V₀ − E₁))/ħ: barriers thinner than a few decay lengths couple adjacent wells into a superlattice miniband.'),
    ],
    series: [
      { key: 'psi', title: 'Potential, levels and eigenfunctions', xlabel: 'x (nm)', ylabel: 'energy (eV)', lines, explain: 'Bound-state wave functions (arbitrary amplitude) drawn at their eigenenergies on top of the well potential. The number of nodes increases by one per level and the tails decay exponentially in the barriers.' },
      { key: 'levels', title: 'Energy levels: finite vs infinite well', xlabel: 'level index n', ylabel: 'energy (eV)', lines: [{ name: 'finite well (numerical)', x: (nB ? bound : [E1]).slice(0, 8).map((_, i) => i + 1), y: (nB ? bound : [E1]).slice(0, 8) }, { name: 'infinite well n²E₁', x: [1, 2, 3, 4, 5, 6, 7, 8].slice(0, Math.max(nB, 1)), y: [1, 2, 3, 4, 5, 6, 7, 8].slice(0, Math.max(nB, 1)).map(Einf) }], explain: 'Finite-well levels lie below the n² ladder of the infinite well, increasingly so near the top of the barrier.' },
    ],
    equations: [eqn('governing', 'Envelope-function Schrödinger equation', '−(ħ²/2) d/dx [ (1/m*(x)) dψ/dx ] + V(x) ψ = E ψ'), eqn('governing', 'Finite-difference Hamiltonian', 'H_ii = (ħ²/2Δ²)(1/m_{i+½} + 1/m_{i−½}) + V_i,   H_{i,i±1} = −ħ²/(2Δ² m_{i±½})'), eqn('constitutive', 'Analytic matching (even states)', '(k/m_w) tan(kW/2) = κ/m_b'), eqn('formula', 'Infinite well', 'E_n = ħ²π²n²/(2 m_w W²)')],
    checks, warnings,
    recommendations: [rec('high', nB > 1 ? `Tune the well width to place the E₁→E₂ transition (${F(E12, 3)} eV) at the target infrared wavelength; it scales roughly as 1/W².` : 'Widen or deepen the well to obtain a second bound state if an intersubband transition is required.', 'Computed level spacing.'), rec('medium', `Keep barriers thicker than about ${F(nB ? 5 / Math.sqrt(Math.max(V0 - E1, 1e-9) * mb / HB2M) : 5, 3)} nm (five decay lengths) to isolate the well, or thinner to couple wells deliberately.`, 'Barrier decay length of the ground state.'), rec('medium', 'A one-monolayer width fluctuation shifts E₁ by 2E₁·ΔW/W: specify growth-rate control accordingly (see the kinetic Monte Carlo mode for layer-by-layer conditions).', 'E ∝ 1/W² sensitivity.'), rec('low', 'Design layer sequences numerically before epitaxy: each avoided calibration growth run saves hours of reactor time, precursor gases and energy.', 'Design-before-build practice.', true)],
    publish: { 'elec.E1_eV': E1, 'geom.thickness': W * 1e-9, 'elec.E12_eV': E12 },
    numerics: { method: 'Finite differences (BenDaniel–Duke operator) + Jacobi eigenvalue solver', dof: ng, steps: 1, iterations: 1, residual: E1a !== null ? Math.abs(E1 - E1a) : 0, notes: `${ng} cell-centred grid points over ${F(Ltot, 3)} nm (${nW} in the well), Dirichlet walls ${F(pad, 3)} nm beyond the well.` },
  };
}
function runTB(inp) {
  const which = pick(inp.tbModel, ['graphene', 'cnt', 'ssh', 'brus', 'well'], 'graphene');
  return { graphene: runGraphene, cnt: runCNT, ssh: runSSH, brus: runBrus, well: runWell }[which](inp, {});
}

// ───────────────────────────── mode 5: lattice dynamics ─────────────────────────────
const THZ = 15.633302;   // √(eV Å⁻² u⁻¹) / 2π in THz
const cvMode = x => (x < 1e-6 ? 1 : x > 500 ? 0 : x * x * Math.exp(x) / (Math.exp(x) - 1) ** 2);   // Einstein function, x = hν/k_BT
const debyeCv = (T, th) => { const xm = th / Math.max(T, 1e-9), n = 200; let s = 0; for (let i = 1; i <= n; i++) { const x = xm * (i - 0.5) / n; s += x * x * cvMode(x); } return 3 * 3 * (s * xm / n) / xm ** 3; };   // per atom, units of k_B
const debye1D = (T, th) => { const xm = th / Math.max(T, 1e-9), n = 200; let s = 0; for (let i = 1; i <= n; i++) s += cvMode(xm * (i - 0.5) / n); return s / n; };   // per mode
const H_K = N.H_PLANCK * 1e12 / N.K_B;   // K per THz
function runChain(inp, diatomic) {
  const m1 = num(inp, 'phM1', 65.38, 1, 250), m2 = diatomic ? num(inp, 'phM2', 15.999, 1, 250) : m1, K = num(inp, 'phK', 50, 0.5, 500), a = num(inp, 'phA', 3.0, 1, 10), T = num(inp, 'T', 94.4, 1, 3000), checks = [], warnings = [];
  const u = AMU, aM = a * 1e-10, nq = 121, qs = N.linspace(0, PI / aM, nq), S = 1 / (m1 * u) + 1 / (m2 * u), toTHz = w => w / (2 * PI) * 1e-12;
  const wAc = q => (diatomic ? Math.sqrt(Math.max(K * S - K * Math.sqrt(Math.max(S * S - 4 * Math.sin(q * aM / 2) ** 2 / (m1 * u * m2 * u), 0)), 0)) : 2 * Math.sqrt(K / (m1 * u)) * Math.abs(Math.sin(q * aM / 2)));
  const wOp = q => Math.sqrt(K * S + K * Math.sqrt(Math.max(S * S - 4 * Math.sin(q * aM / 2) ** 2 / (m1 * u * m2 * u), 0)));
  const vs = diatomic ? aM * Math.sqrt(K / (2 * (m1 + m2) * u)) : aM * Math.sqrt(K / (m1 * u)), wMaxAc = wAc(PI / aM), wOpG = diatomic ? wOp(0) : 0, wOpB = diatomic ? wOp(PI / aM) : 0, gap = diatomic ? toTHz(wOpB - wMaxAc) : 0, wMax = diatomic ? wOpG : wMaxAc;
  const nAt = diatomic ? 2 : 1, thD = N.HBAR * vs * PI * nAt / aM / N.K_B, thE = N.HBAR * wMax / N.K_B / Math.SQRT2, vg = qs.map((q, i) => { const h = PI / aM * 1e-4, f = wAc; return (f(Math.min(q + h, PI / aM)) - f(Math.max(q - h, 0))) / (Math.min(q + h, PI / aM) - Math.max(q - h, 0)); });
  // heat capacity per atom from the exact chain spectrum (1-D), in units of k_B
  const nS = 600, cvChain = Tt => { let s = 0; for (let i = 0; i < nS; i++) { const q = PI / aM * (i + 0.5) / nS; s += cvMode(N.HBAR * wAc(q) / (N.K_B * Tt)); if (diatomic) s += cvMode(N.HBAR * wOp(q) / (N.K_B * Tt)); } return s / (nS * nAt); };
  const Ts = N.linspace(Math.max(thD / 40, 1), Math.max(1.6 * thD, 1.2 * T), 60), vsNum = wAc(PI / aM * 1e-4) / (PI / aM * 1e-4);
  checks.push(chk('Acoustic sum rule', wAc(0) < 1e-3 * wMax, `ω_acoustic(q = 0) = ${F(toTHz(wAc(0)), 2)} THz`), chk('Long-wavelength sound velocity', Math.abs(vsNum - vs) / vs < 1e-3, `numerical slope ${F(vsNum, 5)} m/s, analytic ${F(vs, 5)} m/s`), chk('Dulong–Petit limit', Math.abs(cvChain(20 * thD) - 1) < 0.01, `C_v(20 θ_D) = ${F(cvChain(20 * thD), 4)} k_B per atom (one vibrational degree of freedom per atom in 1-D)`));
  if (diatomic) checks.push(chk('Band-edge frequencies', Math.abs(wMaxAc - Math.sqrt(2 * K / (Math.max(m1, m2) * u))) / wMaxAc < 1e-6 && Math.abs(wOpB - Math.sqrt(2 * K / (Math.min(m1, m2) * u))) / wOpB < 1e-6, `zone-boundary acoustic ${F(toTHz(wMaxAc), 4)} THz = √(2K/M)/2π, optical ${F(toTHz(wOpB), 4)} THz = √(2K/m)/2π`));
  warnings.push('One-dimensional harmonic chain with nearest-neighbour springs: a pedagogical model of phonon branches and gaps, also the reference for mass-spring acoustic metamaterials. Real crystals need the three-dimensional dynamical matrix (see the fcc option).');
  return {
    summary: `${diatomic ? 'Diatomic' : 'Monatomic'} harmonic chain (${diatomic ? `masses ${F(m1, 4)} and ${F(m2, 4)} u` : `mass ${F(m1, 4)} u`}, spring constant ${F(K, 3)} N/m, period ${F(a, 3)} Å): sound velocity ${F(vs, 4)} m/s, maximum frequency ${F(toTHz(wMax), 4)} THz${diatomic ? `, acoustic branch up to ${F(toTHz(wMaxAc), 4)} THz and an acoustic–optical gap of ${F(gap, 4)} THz` : ''}. Debye temperature estimate ${F(thD, 4)} K; heat capacity at ${F(T, 4)} K is ${F(cvChain(T), 3)} k_B per atom.`,
    kpis: [
      kpi('vs', 'Sound velocity', vs, 'm/s', diatomic ? 'Long-wavelength slope of the acoustic branch, a√(K/2(m₁ + m₂)): both atoms move in phase as an elastic continuum.' : 'Long-wavelength group velocity a√(K/m) of the acoustic branch.'),
      kpi('nu_max', 'Maximum frequency', toTHz(wMax), 'THz', diatomic ? 'Optical-branch frequency at q = 0, √(2K(1/m₁ + 1/m₂))/2π: the two sublattices vibrate against each other (infrared-active in ionic crystals).' : 'Zone-boundary frequency 2√(K/m)/2π where neighbouring atoms move in antiphase and the group velocity vanishes.'),
      kpi('gap', 'Acoustic–optical gap', gap, 'THz', diatomic ? `Forbidden frequency band between √(2K/M) and √(2K/m) at the zone boundary; it closes when the masses are equal (ratio here ${F(Math.max(m1, m2) / Math.min(m1, m2), 3)}).` : 'A monatomic chain has a single acoustic branch and no gap.'),
      kpi('thD', 'Debye temperature', thD, 'K', 'θ_D = ħ v_s k_D/k_B with the one-dimensional Debye wave vector k_D = π × (atoms per unit length); above it all modes are excited.'),
      kpi('cv', `Heat capacity at ${F(T, 4)} K`, cvChain(T), 'k_B per atom', `Harmonic heat capacity summed over the exact chain spectrum; the classical limit is 1 k_B per atom in one dimension. The 1-D Debye model gives ${F(debye1D(T, thD), 3)} k_B.`),
      kpi('zpe', 'Zero-point energy', (() => { let s = 0; for (let i = 0; i < nS; i++) { const q = PI / aM * (i + 0.5) / nS; s += wAc(q) + (diatomic ? wOp(q) : 0); } return 0.5 * N.HBAR * s / (nS * nAt) / N.Q_E * 1e3; })(), 'meV per atom', 'Mean ½ħω over all modes: the quantum vibrational energy present even at 0 K.'),
      kpi('thE', 'Einstein temperature (comparison model)', thE, 'K', 'Single-frequency model with ω_E = ω_max/√2; it fails at low temperature, where long-wavelength acoustic modes dominate.'),
    ],
    series: [
      { key: 'disp', title: 'Phonon dispersion', xlabel: 'q a / π', ylabel: 'frequency (THz)', lines: [{ name: 'acoustic', x: qs.map(q => q * aM / PI), y: qs.map(q => toTHz(wAc(q))) }, ...(diatomic ? [{ name: 'optical', x: qs.map(q => q * aM / PI), y: qs.map(q => toTHz(wOp(q))) }] : [])], explain: 'Frequency against wave vector over half the Brillouin zone. The acoustic branch starts linearly (sound) and flattens at the zone boundary; the optical branch of a two-atom cell is separated by a gap that grows with the mass contrast.' },
      { key: 'vg', title: 'Acoustic group velocity', xlabel: 'q a / π', ylabel: 'group velocity (m/s)', lines: [{ name: 'dω/dq', x: qs.map(q => q * aM / PI), y: vg }], explain: 'Energy-transport velocity of acoustic phonons: equal to the sound velocity at long wavelength and zero at the zone boundary (standing wave).' },
      { key: 'cv', title: 'Heat capacity vs temperature', xlabel: 'temperature (K)', ylabel: 'C_v (k_B per atom)', lines: [{ name: 'exact chain spectrum', x: Ts, y: Ts.map(cvChain) }, { name: 'Debye (1-D)', x: Ts, y: Ts.map(t => debye1D(t, thD)) }, { name: 'Einstein', x: Ts, y: Ts.map(t => cvMode(thE / t)) }], explain: 'Harmonic heat capacity rising to the classical limit. The Debye model captures the low-temperature power law (linear in T in one dimension); the Einstein model falls off exponentially and underestimates it.' },
    ],
    equations: [eqn('governing', 'Equations of motion', diatomic ? 'm₁ ü_n = K(v_n + v_{n−1} − 2u_n),   m₂ v̈_n = K(u_{n+1} + u_n − 2v_n)' : 'm ü_n = K(u_{n+1} + u_{n−1} − 2u_n)'), eqn('constitutive', 'Dispersion relation', diatomic ? 'ω² = K(1/m₁ + 1/m₂) ± K√[(1/m₁ + 1/m₂)² − 4 sin²(qa/2)/(m₁m₂)]' : 'ω(q) = 2√(K/m) |sin(qa/2)|'), eqn('constitutive', 'Harmonic heat capacity', 'C_v = k_B Σ_modes x² eˣ/(eˣ − 1)²,  x = ħω/k_BT'), eqn('formula', 'Debye temperature (1-D)', 'θ_D = ħ v_s k_D / k_B,  k_D = π N_atoms/L')],
    checks, warnings,
    recommendations: [rec('high', diatomic ? `The gap between ${F(toTHz(wMaxAc), 3)} and ${F(toTHz(wOpB), 3)} THz blocks vibrational transport at those frequencies; widen it by increasing the mass contrast, the same principle used in phononic-crystal and locally resonant metamaterial band gaps.` : 'Introduce a second mass or spring per cell to open a band gap; a uniform chain transmits all frequencies up to the cut-off.', 'Dispersion relation.'), rec('medium', `Fit the spring constant to a measured sound velocity or zone-boundary frequency (currently ${F(vs, 4)} m/s, ${F(toTHz(wMax), 3)} THz) before using the model for thermal estimates.`, 'Two-parameter model calibration.'), rec('medium', `Below about ${F(thD / 10, 3)} K (θ_D/10) only long-wavelength acoustic modes are populated; nanostructure boundaries then dominate phonon scattering and thermal conductivity.`, 'Debye temperature.'), rec('low', 'Phonon band-gap and low-group-velocity designs lower lattice thermal conductivity, improving thermoelectric efficiency and thermal insulation with no added material.', 'Group-velocity curve.', true)],
    publish: { 'prop.c_sound': vs, 'therm.debye_T': thD, 'res.bandgap_lo': diatomic ? toTHz(wMaxAc) * 1e12 : 0, 'res.bandgap_hi': diatomic ? toTHz(wOpB) * 1e12 : 0, 'res.f0': toTHz(wMax) * 1e12 },
    numerics: { method: 'Closed-form harmonic dispersion; mode sums on a uniform q-grid', dof: nAt, steps: nS, iterations: 1, residual: 0, notes: `${nS} q-points for thermodynamic sums.` },
  };
}
function runFccPhonon(inp) {
  const pk = pick(inp.mdPot, Object.keys(POT), 'Ar'), pt = POT[pk], lj = pt.type === 'lj', T = num(inp, 'T', 94.4, 1, 3000), checks = [], warnings = [], Rc = 5.0, S2 = Math.SQRT2, vec = [];
  for (let i = -4; i <= 4; i++) for (let j = -4; j <= 4; j++) for (let k = -4; k <= 4; k++) for (const b of [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) { const x = (i + b[0]) * S2, y = (j + b[1]) * S2, z = (k + b[2]) * S2, r = Math.hypot(x, y, z); if (r > 1e-9 && r <= Rc) vec.push([x, y, z, r]); }
  const phi = phiOf(pt), dphi = dphiOf(pt), d2phi = d2phiOf(pt), Elat = d => { let e = 0; for (const v of vec) e += phi(d * v[3]); return e / 2; }, dG = lj ? 1.09 * pt.r : pt.r;
  const scan = N.linspace(0.5 * dG, 1.5 * dG, 41), im = N.argmin(scan.map(Elat)), d0 = N.goldenSection(Elat, scan[Math.max(im - 1, 0)], scan[Math.min(im + 1, 40)], 1e-10 * dG), a = S2 * d0, m = pt.mass;
  const Aq = vec.map(v => { const r = d0 * v[3]; return d2phi(r) - dphi(r) / r; }), Bq = vec.map(v => dphi(d0 * v[3]) / (d0 * v[3]));
  const freqs = (qx, qy, qz) => { let xx = 0, yy = 0, zz = 0, xy = 0, xz = 0, yz = 0; for (let i = 0; i < vec.length; i++) { const v = vec[i], c = 1 - Math.cos(d0 * (qx * v[0] + qy * v[1] + qz * v[2])), A1 = Aq[i] * c / (v[3] * v[3]), B1 = Bq[i] * c; xx += A1 * v[0] * v[0] + B1; yy += A1 * v[1] * v[1] + B1; zz += A1 * v[2] * v[2] + B1; xy += A1 * v[0] * v[1]; xz += A1 * v[0] * v[2]; yz += A1 * v[1] * v[2]; } return N.eigSym([[xx / m, xy / m, xz / m], [xy / m, yy / m, yz / m], [xz / m, yz / m, zz / m]]).values.map(l => Math.sign(l) * Math.sqrt(Math.abs(l)) * THZ); };
  const g = 2 * PI / a, segs = [[[0, 0, 0], [1, 0, 0], 'Γ', 'X'], [[1, 1, 0], [0, 0, 0], 'X', 'Γ (via K)'], [[0, 0, 0], [.5, .5, .5], 'Γ', 'L']], xs = [], br = [[], [], []], ticks = []; let dist = 0;
  for (const [p, q, n1, n2] of segs) { const len = Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]), np = Math.round(40 * len); ticks.push(`${n1} = ${F(dist, 3)}`); for (let i = 0; i <= np; i++) { const t = i / np, f = freqs(g * (p[0] + t * (q[0] - p[0])), g * (p[1] + t * (q[1] - p[1])), g * (p[2] + t * (q[2] - p[2]))); xs.push(dist + t * len); for (let b = 0; b < 3; b++) br[b].push(f[b]); } dist += len; ticks.push(`${n2} = ${F(dist, 3)}`); }
  // Brillouin-zone sampling (uniform grid over the cubic reciprocal cell of side 4π/a, which holds two zones)
  const ng = 12, all = []; let minF = Infinity; for (let i = 0; i < ng; i++) for (let j = 0; j < ng; j++) for (let k = 0; k < ng; k++) { const f = freqs(2 * g * (i + 0.5) / ng, 2 * g * (j + 0.5) / ng, 2 * g * (k + 0.5) / ng); for (const v of f) { all.push(v); if (v < minF) minF = v; } }
  const nuMax = N.maxOf(all), dos = dosHist(all.filter(v => v > 0), 0, nuMax * 1.02, 60, 3 / all.length), fX = freqs(g, 0, 0), fL = freqs(g / 2, g / 2, g / 2);
  // elastic constants from sound velocities
  const rho = 4 * m * AMU / (a * 1e-10) ** 3, eps = 0.02, vel = (dx, dy, dz) => freqs(g * eps * dx, g * eps * dy, g * eps * dz).map(f => 2 * PI * f * 1e12 / (g * 1e10 * eps * Math.hypot(dx, dy, dz)));
  const v100 = vel(1, 0, 0), v110 = vel(1, 1, 0), C11 = rho * v100[2] ** 2, C44 = rho * v100[0] ** 2, C12 = 2 * rho * v110[2] ** 2 - C11 - 2 * C44, Bm = (C11 + 2 * C12) / 3, rd = N.rng(7); let s3 = 0; const nd = 80;
  for (let i = 0; i < nd; i++) { const zc = 2 * rd() - 1, ph = 2 * PI * rd(), sq = Math.sqrt(1 - zc * zc), v = vel(sq * Math.cos(ph), sq * Math.sin(ph), zc); s3 += (1 / v[0] ** 3 + 1 / v[1] ** 3 + 1 / v[2] ** 3) / 3; }
  const vD = (s3 / nd) ** (-1 / 3), nAt = 4 / (a * 1e-10) ** 3, thD = N.HBAR * vD / N.K_B * Math.cbrt(6 * PI * PI * nAt), cv = Tt => { let s = 0; for (const f of all) if (f > 0) s += cvMode(H_K * f / Tt); return 3 * s / all.length; }, Ts = N.linspace(Math.max(thD / 30, 1), Math.max(1.5 * thD, 1.2 * T), 50);
  const zpe = 1.5 * N.mean(all.filter(v => v > 0)) * 1e12 * N.H_PLANCK / N.Q_E * 1e3, lindT = (() => { let s = 0; for (const f of all) if (f > 0) { const x = H_K * f / T; s += (1 / Math.tanh(x / 2)) / f; } const u2 = 3 * N.HBAR / (2 * m * AMU * 2 * PI * 1e12) * s / all.length; return Math.sqrt(u2) / (d0 * 1e-10); })();
  checks.push(chk('Dynamical stability (no imaginary modes)', minF > -1e-6 * nuMax, `lowest sampled frequency ${F(minF, 3)} THz over ${ng ** 3} q-points`), chk('Acoustic sum rule at Γ', Math.abs(freqs(0, 0, 0)[2]) < 1e-6, 'all three branches vanish at q = 0'), chk('Cauchy relation C₁₂ = C₄₄', Math.abs(C12 - C44) / C44 < 0.03, `C₁₂ = ${F(C12 / 1e9, 4)} GPa, C₄₄ = ${F(C44 / 1e9, 4)} GPa (exact for central pair forces at zero pressure)`), chk('Dulong–Petit limit', Math.abs(cv(30 * thD) - 3) < 0.02, `C_v(30 θ_D) = ${F(cv(30 * thD), 4)} k_B per atom`), chk('Maximum frequency within 30 % of experiment', Math.abs(nuMax - pt.nuMax) / pt.nuMax < 0.3, `${F(nuMax, 4)} THz vs about ${pt.nuMax} THz measured`));
  warnings.push(`Harmonic lattice dynamics with a central pair potential summed to ${Rc} nearest-neighbour distances at its own static equilibrium (a = ${F(a, 5)} Å): the Cauchy relation C₁₂ = C₄₄ is forced, which real metals violate, so transverse branches and elastic anisotropy are approximate.`);
  if (lj) warnings.push('Anharmonicity and zero-point expansion, which are large in rare-gas solids, are not included.');
  return {
    summary: `Harmonic phonons of fcc ${pk} from the ${pt.name} potential at a = ${F(a, 5)} Å: maximum frequency ${F(nuMax, 4)} THz (measured ≈ ${pt.nuMax} THz), X-point modes ${F(fX[0], 3)} (T) and ${F(fX[2], 3)} (L) THz. Elastic constants C₁₁ = ${F(C11 / 1e9, 4)}, C₁₂ = ${F(C12 / 1e9, 4)}, C₄₄ = ${F(C44 / 1e9, 4)} GPa; Debye temperature ${F(thD, 4)} K (measured ≈ ${pt.thD} K); C_v(${F(T, 4)} K) = ${F(cv(T), 3)} k_B per atom.`,
    kpis: [
      kpi('nu_max', 'Maximum phonon frequency', nuMax, 'THz', `Highest frequency in the sampled Brillouin zone (longitudinal mode near the zone boundary); the measured value for ${pk} is about ${pt.nuMax} THz.`, Math.abs(nuMax - pt.nuMax) / pt.nuMax < 0.3),
      kpi('thD', 'Debye temperature', thD, 'K', `θ_D = (ħ v_D/k_B)(6π²n)^⅓ with the Debye-averaged sound velocity v_D = ${F(vD, 4)} m/s from ${nd} propagation directions; the low-temperature calorimetric value is about ${pt.thD} K.`, Math.abs(thD - pt.thD) / pt.thD < 0.3),
      kpi('cv', `Heat capacity at ${F(T, 4)} K`, cv(T), 'k_B per atom', `Harmonic C_v from the sampled phonon spectrum (classical limit 3 k_B); the Debye model with the same θ_D gives ${F(debyeCv(T, thD), 3)} k_B. Molar value ${F(cv(T) * N.R_GAS, 4)} J/mol/K.`),
      kpi('C11', 'Elastic constant C₁₁', C11, 'Pa', 'From the longitudinal sound velocity along [100], C₁₁ = ρ v_L².'),
      kpi('C44', 'Elastic constant C₄₄', C44, 'Pa', 'From the transverse sound velocity along [100], C₄₄ = ρ v_T²; equal to C₁₂ for any central pair potential.'),
      kpi('B', 'Bulk modulus', Bm, 'Pa', 'B = (C₁₁ + 2C₁₂)/3 from the long-wavelength limit of the dynamical matrix; consistent with the curvature of the lattice energy.'),
      kpi('vD', 'Debye sound velocity', vD, 'm/s', `Direction- and branch-averaged velocity, (⟨Σ v⁻³⟩/3)^(−⅓); longitudinal [100] velocity ${F(v100[2], 4)} m/s, transverse ${F(v100[0], 4)} m/s.`),
      kpi('zpe', 'Zero-point energy', zpe, 'meV per atom', 'Mean (3/2)hν over the spectrum; for argon it is about a tenth of the cohesive energy, which explains the expanded lattice of the real crystal.'),
      kpi('lind', `Lindemann ratio at ${F(T, 4)} K`, lindT, '–', 'Harmonic root-mean-square displacement √⟨u²⟩ divided by the nearest-neighbour distance; melting is expected when it reaches roughly 0.1–0.15.', lindT < 0.12),
    ],
    series: [
      { key: 'disp', title: `Phonon dispersion of fcc ${pk}`, xlabel: `reduced wave vector (2π/a): ${ticks.filter((t, i) => i % 2 === 0 || i === ticks.length - 1).join(', ')}`, ylabel: 'frequency (THz)', lines: [{ name: 'transverse 1', x: xs, y: br[0] }, { name: 'transverse 2', x: xs, y: br[1] }, { name: 'longitudinal', x: xs, y: br[2] }], explain: 'Three acoustic branches along Γ–X [ξ00], X–K–Γ [ξξ0] and Γ–L [ξξξ]. Transverse branches are degenerate along the four-fold and three-fold axes and split along [ξξ0].' },
      { key: 'dos', title: 'Phonon density of states', xlabel: 'frequency (THz)', ylabel: 'modes per THz per atom', lines: [{ name: 'g(ν)', x: dos.x, y: dos.y }], explain: 'Histogram over the Brillouin zone: quadratic (Debye) rise at low frequency, then peaks from the flat transverse and longitudinal zone-boundary modes.' },
      { key: 'cv', title: 'Heat capacity vs temperature', xlabel: 'temperature (K)', ylabel: 'C_v (k_B per atom)', lines: [{ name: 'lattice dynamics', x: Ts, y: Ts.map(cv) }, { name: `Debye, θ_D = ${F(thD, 3)} K`, x: Ts, y: Ts.map(t => debyeCv(t, thD)) }, { name: 'Einstein, θ_E = 0.806 θ_D', x: Ts, y: Ts.map(t => 3 * cvMode(0.806 * thD / t)) }], explain: 'Harmonic heat capacity from the full spectrum against the Debye and Einstein models; all converge to 3 k_B (Dulong–Petit) at high temperature.' },
    ],
    equations: [eqn('formula', 'Crystal', `${chem.prettyFormula(pk)}, fcc, a = ${F(a, 5)} Å, M = ${F(chem.molarMass(pk), 6)} g/mol`), eqn('governing', 'Dynamical matrix', 'D_αβ(q) = (1/m) Σ_R Φ_αβ(R) [1 − cos(q·R)],   Φ_αβ = (φ″ − φ′/R) R_αR_β/R² + δ_αβ φ′/R'), eqn('governing', 'Eigenvalue problem', 'D(q) e = ω² e'), eqn('constitutive', 'Pair potential', potLine(pt)), eqn('constitutive', 'Debye temperature and heat capacity', 'θ_D = (ħ v_D/k_B)(6π²n)^⅓;   C_v = k_B Σ_{q,s} x² eˣ/(eˣ − 1)²,  x = ħω_s(q)/k_BT')],
    geometry: (() => { const sc = A.supercell(A.buildCell(fccDef(pk, a)), [3, 3, 3]); return { kind: 'atoms', atoms: sc.atoms.map(q => ({ el: q.el, x: q.x, y: q.y, z: q.z })), cell: sc.cell, title: `${pk} fcc 3×3×3 at the potential's equilibrium` }; })(),
    checks, warnings,
    recommendations: [rec('high', `Use θ_D = ${F(thD, 3)} K as the scale for thermal properties: below about ${F(thD / 10, 3)} K the heat capacity follows T³ and boundary scattering in nanostructures controls heat transport.`, 'Computed Debye temperature.'), rec('medium', lj ? 'For rare-gas solids add quasi-harmonic zero-point and thermal expansion (repeat at the expanded lattice parameter) before comparing with measured frequencies.' : 'For quantitative metal phonons use an embedded-atom or DFT force-constant model: the pair potential enforces C₁₂ = C₄₄ and cannot reproduce the measured elastic anisotropy.', 'Cauchy-relation check.'), rec('medium', 'Phonon confinement in particles or films thinner than a few nanometres cuts off wavelengths longer than the size and softens surface modes; lower the effective Debye temperature accordingly in nanoscale thermal models.', 'Long-wavelength acoustic branches.'), rec('low', 'Harmonic lattice dynamics from a potential costs milliseconds and already gives θ_D, C_v and elastic constants; reserve DFT phonon calculations for final validation of selected materials.', 'Method hierarchy.', true)],
    publish: { 'material.formula': pk, 'therm.debye_T': thD, 'prop.c_sound': vD, 'prop.K': Bm, 'prop.cp': cv(T) * N.R_GAS / (chem.molarMass(pk) * 1e-3), 'prop.lattice_a': a * 1e-10, 'prop.density': rho },
    numerics: { method: 'Born–von Kármán lattice dynamics: 3×3 dynamical matrix from pair-potential force constants, Jacobi diagonalisation', dof: 3, steps: ng ** 3 + xs.length, iterations: vec.length, residual: Math.abs(freqs(0, 0, 0)[2]), notes: `${vec.length} neighbours within ${Rc} nearest-neighbour distances; ${ng}³ q-point grid for the spectrum.` },
  };
}
function runPhonon(inp) { const w = pick(inp.phModel, ['mono', 'diatomic', 'fcc'], 'diatomic'); return w === 'fcc' ? runFccPhonon(inp) : runChain(inp, w === 'diatomic'); }

// ───────────────────────────── mode 6: reaction and diffusion pathways ─────────────────────────────
const MB = { A: [-200, -100, -170, 15], a: [-1, -1, -6.5, 0.7], b: [0, 0, 11, 0.6], c: [-10, -10, -6.5, 0.7], x0: [1, 0, -0.5, -1], y0: [0, 0.5, 1.5, 1] };
function pesMB(x, y) { let v = 0, gx = 0, gy = 0; for (let k = 0; k < 4; k++) { const dx = x - MB.x0[k], dy = y - MB.y0[k], e = MB.A[k] * Math.exp(Math.min(MB.a[k] * dx * dx + MB.b[k] * dx * dy + MB.c[k] * dy * dy, 50)); v += e; gx += e * (2 * MB.a[k] * dx + MB.b[k] * dy); gy += e * (MB.b[k] * dx + 2 * MB.c[k] * dy); } return [v, gx, gy]; }
function pesHop(x, y) { const t = 2 * PI; return [50 * (2 - Math.cos(t * x) - Math.cos(t * y)), 50 * t * Math.sin(t * x), 50 * t * Math.sin(t * y)]; }
function hess(V, x, y) { const h = 1e-4, f = (p, q) => V(p, q), xx = (f(x + h, y)[1] - f(x - h, y)[1]) / (2 * h), yy = (f(x, y + h)[2] - f(x, y - h)[2]) / (2 * h), xy = (f(x + h, y)[2] - f(x - h, y)[2]) / (2 * h), tr = xx + yy, dt = Math.sqrt(Math.max((xx - yy) ** 2 / 4 + xy * xy, 0)); return [tr / 2 - dt, tr / 2 + dt]; }
function runNEB(inp) {
  const which = pick(inp.pes, ['mb', 'hop'], 'mb'), V = which === 'mb' ? pesMB : pesHop, pathSel = pick(inp.nebPath, ['AC', 'AB', 'BC'], 'AC'), nImg = int(inp, 'nImg', 15, 7, 25), sE = num(inp, 'pesE', 0.005, 1e-4, 0.05), sL = num(inp, 'pesL', 1, 0.2, 5), mass = num(inp, 'nebMass', 63.546, 1, 250), nu0 = num(inp, 'nebNu', 1e13, 1e10, 1e15), dim = int(inp, 'dim', 2, 1, 3), T = num(inp, 'T_C', 300, -150, 1500) + 273.15, checks = [], warnings = [];
  const relax = p0 => { const r = N.nelderMead(p => V(p[0], p[1])[0], p0, { maxIter: 800, tol: 1e-15, step: 0.03 }), p = r.x ?? r.p ?? r; let x = p[0], y = p[1]; for (let i = 0; i < 200; i++) { const [, gx, gy] = V(x, y), h = 1e-4, hxx = (V(x + h, y)[1] - V(x - h, y)[1]) / (2 * h), hyy = (V(x, y + h)[2] - V(x, y - h)[2]) / (2 * h), hxy = (V(x + h, y)[2] - V(x - h, y)[2]) / (2 * h), det = hxx * hyy - hxy * hxy; if (!(det > 0)) break; const dx = (hyy * gx - hxy * gy) / det, dy = (hxx * gy - hxy * gx) / det; x -= dx; y -= dy; if (Math.hypot(dx, dy) < 1e-12) break; } return [x, y]; };
  let P0, P1, names;
  if (which === 'mb') { const mins = { A: relax([-0.558, 1.442]), B: relax([-0.05, 0.467]), C: relax([0.623, 0.028]) }; P0 = mins[pathSel[0]]; P1 = mins[pathSel[1]]; names = [`minimum ${pathSel[0]}`, `minimum ${pathSel[1]}`]; } else { P0 = [0, 0]; P1 = [1, 0]; names = ['site (0, 0)', 'site (1, 0)']; }
  const n = nImg, X = new Float64Array(n), Y = new Float64Array(n), vx = new Float64Array(n), vy = new Float64Array(n), E = new Float64Array(n), GX = new Float64Array(n), GY = new Float64Array(n), FX = new Float64Array(n), FY = new Float64Array(n);
  for (let i = 0; i < n; i++) { const t = i / (n - 1); X[i] = P0[0] + t * (P1[0] - P0[0]); Y[i] = P0[1] + t * (P1[1] - P0[1]) + (which === 'hop' ? 0.18 * Math.sin(PI * t) : 0); }
  const kS = 300, dt = 0.0015, maxIt = 5000, tol = 0.05, hist = []; let it = 0, fmax = Infinity, climb = -1;
  const evalF = () => {
    for (let i = 0; i < n; i++) { const r = V(X[i], Y[i]); E[i] = r[0]; GX[i] = r[1]; GY[i] = r[2]; }
    fmax = 0; for (let i = 1; i < n - 1; i++) {
      const px = X[i + 1] - X[i], py = Y[i + 1] - Y[i], mx = X[i] - X[i - 1], my = Y[i] - Y[i - 1]; let tx, ty;
      if (E[i + 1] > E[i] && E[i] > E[i - 1]) { tx = px; ty = py; } else if (E[i + 1] < E[i] && E[i] < E[i - 1]) { tx = mx; ty = my; } else { const d1 = Math.abs(E[i + 1] - E[i]), d2 = Math.abs(E[i - 1] - E[i]), hi = Math.max(d1, d2), lo = Math.min(d1, d2); if (E[i + 1] > E[i - 1]) { tx = px * hi + mx * lo; ty = py * hi + my * lo; } else { tx = px * lo + mx * hi; ty = py * lo + my * hi; } }
      const tn = Math.hypot(tx, ty) || 1; tx /= tn; ty /= tn; const gp = GX[i] * tx + GY[i] * ty;
      if (i === climb) { FX[i] = -GX[i] + 2 * gp * tx; FY[i] = -GY[i] + 2 * gp * ty; } else { const sp = kS * (Math.hypot(px, py) - Math.hypot(mx, my)); FX[i] = -GX[i] + gp * tx + sp * tx; FY[i] = -GY[i] + gp * ty + sp * ty; }
      const perp = i === climb ? Math.hypot(FX[i], FY[i]) : Math.hypot(-GX[i] + gp * tx, -GY[i] + gp * ty); if (perp > fmax) fmax = perp;
    }
  };
  evalF();
  for (it = 0; it < maxIt; it++) {
    if (it === 400) { climb = 1; for (let i = 2; i < n - 1; i++) if (E[i] > E[climb]) climb = i; evalF(); }
    if (it % 10 === 0) hist.push([it, fmax]); if (fmax < tol && it > 420) break;
    for (let i = 1; i < n - 1; i++) { vx[i] += dt * FX[i]; vy[i] += dt * FY[i]; const vf = vx[i] * FX[i] + vy[i] * FY[i], ff = FX[i] * FX[i] + FY[i] * FY[i]; if (vf <= 0 || ff === 0) { vx[i] = 0; vy[i] = 0; } else { vx[i] = vf / ff * FX[i]; vy[i] = vf / ff * FY[i]; } let sx = dt * vx[i], sy = dt * vy[i]; const sn = Math.hypot(sx, sy); if (sn > 0.02) { sx *= 0.02 / sn; sy *= 0.02 / sn; } X[i] += sx; Y[i] += sy; }
    evalF();
  }
  hist.push([it, fmax]);
  const iS = N.argmax(Array.from(E)), Es = E[iS], bF = Es - E[0], bR = Es - E[n - 1], dE = E[n - 1] - E[0], hMin = hess(V, P0[0], P0[1]), hSad = hess(V, X[iS], Y[iS]);
  const conv = sE * N.Q_E / (mass * AMU * (sL * 1e-10) ** 2), nuV = hMin[0] > 0 && hSad[1] > 0 ? Math.sqrt(hMin[0] * hMin[1] / hSad[1] * conv) / (2 * PI) : nu0, wIm = hSad[0] < 0 ? Math.sqrt(-hSad[0] * conv) : 0, Ea = bF * sE, kT = KB * T;
  const kV = nuV * Math.exp(-Ea / kT), kU = nu0 * Math.exp(-Ea / kT), jump = Math.hypot(P1[0] - P0[0], P1[1] - P0[1]) * sL * 1e-10, D = kV * jump * jump / (2 * dim), D0 = nuV * jump * jump / (2 * dim), wig = 1 + (N.HBAR * wIm / (N.K_B * T)) ** 2 / 24;
  const arc = [0]; for (let i = 1; i < n; i++) arc.push(arc[i - 1] + Math.hypot(X[i] - X[i - 1], Y[i] - Y[i - 1]) * sL);
  // potential-energy map
  const dom = which === 'mb' ? [-1.5, 1.2, -0.3, 2.0] : [-0.3, 1.3, -0.5, 0.5], nx = 90, ny = 80, data = new Float64Array(nx * ny), clip = (Es + 1.2 * Math.max(bF, bR, 1) - E[0]) * sE;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) data[j * nx + i] = Math.min((V(dom[0] + (dom[1] - dom[0]) * i / (nx - 1), dom[2] + (dom[3] - dom[2]) * j / (ny - 1))[0] - E[0]) * sE, clip);
  const invT = N.linspace(1000 / Math.max(T * 2, 400), 1000 / Math.max(T * 0.5, 150), 40), sadOK = hSad[0] < 0 && hSad[1] > 0;
  checks.push(chk('NEB converged (perpendicular force below tolerance)', fmax < tol, `max force ${F(fmax, 3)} energy-units per length after ${it} iterations (tolerance ${tol})`), chk('Climbing image is a first-order saddle', sadOK, `Hessian eigenvalues ${F(hSad[0], 4)} and ${F(hSad[1], 4)} (one negative required)`), chk('Gradient at the saddle', Math.hypot(GX[iS], GY[iS]) < 1, `|∇V| = ${F(Math.hypot(GX[iS], GY[iS]), 3)} at (${F(X[iS], 4)}, ${F(Y[iS], 4)})`));
  if (which === 'mb' && pathSel !== 'BC') checks.push(chk('Müller–Brown saddle energy against the reference', Math.abs(Es + 40.665) < 0.3, `E‡ = ${F(Es, 6)} (reference −40.665 at (−0.822, 0.624))`));
  if (which === 'mb' && pathSel === 'BC') checks.push(chk('Müller–Brown saddle energy against the reference', Math.abs(Es + 72.249) < 0.3, `E‡ = ${F(Es, 6)} (reference −72.249 at (0.212, 0.293))`));
  if (which === 'hop') checks.push(chk('Barrier against the analytic value', Math.abs(bF - 100) < 0.5, `NEB ${F(bF, 6)} vs 100 energy units at the bridge site (0.5, 0)`));
  if (Ea / kT < 3) warnings.push(`The barrier is only ${F(Ea / kT, 3)} k_BT: harmonic transition-state theory assumes rare events (E_a ≳ 5 k_BT) and overestimates the rate here.`);
  if (wig > 1.5) warnings.push(`Tunnelling correction ${F(wig, 3)} is large: the Wigner formula is not valid; use instanton or path-integral methods for light atoms at low temperature.`);
  warnings.push('Model two-dimensional potential-energy surface scaled to physical units by the chosen energy and length units; it demonstrates the method. For a real process take the images from DFT or a force field and use the same NEB procedure in an external code.');
  return {
    summary: `Climbing-image nudged elastic band with ${n} images on the ${which === 'mb' ? 'Müller–Brown' : 'periodic surface-hopping'} surface from ${names[0]} to ${names[1]}: saddle at (${F(X[iS], 4)}, ${F(Y[iS], 4)}) with a forward barrier of ${F(bF, 5)} energy units = ${F(Ea, 4)} eV (reverse ${F(bR * sE, 4)} eV, reaction energy ${F(dE * sE, 4)} eV). Harmonic TST with the Vineyard prefactor ${F(nuV, 3)} Hz gives k = ${F(kV, 3)} s⁻¹ at ${F(T, 4)} K and a diffusion coefficient of ${F(D, 3)} m²/s for a ${F(jump * 1e10, 3)} Å jump.`,
    kpis: [
      kpi('Ea', 'Forward activation energy', Ea, 'eV', `Energy of the climbing image minus the initial minimum, ${F(bF, 5)} surface units × ${F(sE, 3)} eV; the rate depends exponentially on it (E_a/k_BT = ${F(Ea / kT, 3)}).`),
      kpi('Ea_rev', 'Reverse activation energy', bR * sE, 'eV', 'Barrier seen from the final state; the difference of forward and reverse barriers is the reaction energy.'),
      kpi('dE', 'Reaction energy', dE * sE, 'eV', 'Energy of the final minimum relative to the initial one (negative = exothermic step).'),
      kpi('nu', 'Vineyard attempt frequency', nuV, 'Hz', `ν* = Π ν_i(minimum)/Π′ ν_j(saddle) from finite-difference Hessians, for a mass of ${F(mass, 4)} u; the commonly assumed value is 10¹²–10¹³ Hz.`),
      kpi('k', 'Rate constant (harmonic TST)', kV, '1/s', `k = ν* exp(−E_a/k_BT) at ${F(T, 4)} K; with the user attempt frequency ${F(nu0, 3)} Hz it would be ${F(kU, 3)} s⁻¹. Mean residence time ${F(1 / Math.max(kV, 1e-300), 3)} s.`),
      kpi('D', 'Diffusion coefficient', D, 'm²/s', `D = k ℓ²/(2d) for uncorrelated hops of length ℓ = ${F(jump * 1e10, 3)} Å in d = ${dim} dimensions; pre-exponential factor D₀ = ${F(D0, 3)} m²/s.`),
      kpi('wigner', 'Wigner tunnelling factor', wig, '–', `κ = 1 + (ħω‡/k_BT)²/24 with the imaginary saddle frequency ${F(wIm / (2 * PI) * 1e-12, 3)} THz; values near 1 mean classical barrier crossing.`, wig < 1.2),
      kpi('saddle_x', 'Saddle position along the path', arc[iS] / arc[n - 1], '–', `Fractional arc length of the transition state; coordinates (${F(X[iS], 4)}, ${F(Y[iS], 4)}) in surface units. An early saddle indicates an exothermic, reactant-like transition state (Hammond postulate).`),
    ],
    series: [
      { key: 'profile', title: 'Minimum-energy-path profile', xlabel: 'reaction coordinate (Å)', ylabel: 'energy (eV)', lines: [{ name: 'NEB images', x: arc, y: Array.from(E, v => (v - E[0]) * sE) }], explain: 'Energy of each image relative to the initial state along the converged path. The highest point is the climbing image at the saddle; intermediate minima appear as dips.' },
      { key: 'path', title: 'Converged path in configuration space', xlabel: 'x (surface units)', ylabel: 'y (surface units)', lines: [{ name: 'minimum-energy path', x: Array.from(X), y: Array.from(Y) }, { name: 'initial straight-line guess', x: [P0[0], P1[0]], y: [P0[1], P1[1]] }], explain: 'Image positions after relaxation compared with the linear interpolation used as the starting guess; the path bends to follow the valley floor through the saddle.' },
      { key: 'arrh', title: 'Arrhenius plot', xlabel: '1000 / T (1/K)', ylabel: 'log₁₀ value', lines: [{ name: 'log₁₀ k (1/s)', x: invT, y: invT.map(v => Math.log10(nuV) - Ea / (KB * 1000 / v) / Math.LN10) }, { name: 'log₁₀ D (m²/s)', x: invT, y: invT.map(v => Math.log10(D0) - Ea / (KB * 1000 / v) / Math.LN10) }], explain: 'Straight lines of slope −E_a/(k_B ln 10): compare with measured rate or diffusion data to validate the barrier.' },
      { key: 'conv', title: 'NEB convergence', xlabel: 'iteration', ylabel: 'log₁₀ max force', lines: [{ name: 'max perpendicular force', ...(() => { const t = thin(hist.map(h => h[0]), hist.map(h => Math.log10(Math.max(h[1], 1e-12))), 300); return t; })() }], explain: 'Largest force component perpendicular to the path (full force on the climbing image). The jump at iteration 400 marks the switch to the climbing image.' },
    ],
    maps: [{ key: 'pes', title: 'Potential-energy surface with the NEB path', nx, ny, data, x0: dom[0], x1: dom[1], y0: dom[2], y1: dom[3], xlabel: 'x (surface units)', ylabel: 'y (surface units)', unit: 'eV', overlay: [{ name: 'NEB path', x: Array.from(X), y: Array.from(Y) }, { name: 'saddle', x: [X[iS]], y: [Y[iS]] }], explain: 'Energy relative to the initial minimum (clipped above the barrier region). The overlaid path follows the valley from the initial to the final minimum and crosses the ridge at its lowest point, the saddle.' }],
    equations: [
      eqn('governing', 'Nudged elastic band force', 'F_i = −∇V(R_i)|_⊥ + k(|R_{i+1} − R_i| − |R_i − R_{i−1}|) τ̂_i;   climbing image: F = −∇V + 2(∇V·τ̂)τ̂'),
      eqn('governing', 'Harmonic transition-state theory', 'k = ν* exp(−E_a/k_BT),   ν* = Π_i^N ν_i^min / Π_j^(N−1) ν_j^‡  (Vineyard)'),
      eqn('constitutive', 'Arrhenius diffusion', 'D = D₀ exp(−E_a/k_BT),  D₀ = ν* ℓ²/(2d)'), eqn('constitutive', 'Wigner correction', 'κ = 1 + (1/24)(ħω‡/k_BT)²'),
      eqn('formula', 'Model surface', which === 'mb' ? 'V(x,y) = Σ_{k=1}^{4} A_k exp[a_k(x − x_k)² + b_k(x − x_k)(y − y_k) + c_k(y − y_k)²]  (Müller–Brown)' : 'V(x,y) = 50[2 − cos 2πx − cos 2πy]  (square-lattice adatom corrugation, bridge-site saddle)'),
    ],
    checks, warnings,
    recommendations: [
      rec('high', `With E_a = ${F(Ea, 3)} eV the process has a rate of ${F(kV, 3)} s⁻¹ at ${F(T, 4)} K; a rate of 1 s⁻¹ is reached near ${F(Ea / (KB * Math.log(Math.max(nuV, 10))), 4)} K — use this as the onset temperature for annealing or growth.`, 'Harmonic TST rate.'),
      rec('medium', `Always converge the climbing image (force ${F(fmax, 2)} here) and confirm exactly one negative Hessian eigenvalue; an unconverged band underestimates the barrier by tens of per cent.`, 'Saddle-point verification checks.'),
      rec('medium', 'Feed the barrier and attempt frequency into the kinetic Monte Carlo mode (surface diffusion barrier) or continuum diffusion models: this is the standard DFT → NEB → KMC multiscale chain.', 'Published activation energy and D.'),
      rec('low', 'A catalyst or surfactant that lowers the barrier by 0.1 eV reduces the temperature needed for the same rate by tens of kelvin, a direct saving in process heat.', 'Exponential dependence of k on E_a.', true),
    ],
    publish: { 'prop.D_diff': D, 'rxn.Ea_eV': Ea, 'rxn.rate_const': kV, 'rxn.attempt_freq': nuV, 'proc.T_C': T - 273.15 },
    numerics: { method: 'Climbing-image nudged elastic band (improved tangent), projected velocity-Verlet (quick-min) optimiser; finite-difference Hessians', dof: 2 * (n - 2), steps: it, iterations: it, residual: fmax, notes: `Spring constant ${kS}, time step ${dt}, step cap 0.02, climbing image from iteration 400.` },
  };
}

// ───────────────────────────── mode 7: input decks for external solvers ─────────────────────────────
const METALS = ['Au', 'Ag', 'Cu', 'Pt', 'Pd', 'Ni', 'Al', 'Fe_bcc', 'W', 'graphite', 'graphene'], MAGNETIC = ['Fe_bcc', 'Ni', 'Fe3O4'];
function runDeck(inp) {
  const key = pick(inp.dkCrystal, KEYS, 'ZnO_wurtzite'), def = A.CRYSTALS[key], rep = int(inp, 'dkRep', 1, 1, 3), ks = num(inp, 'kspacing', 0.25, 0.05, 1), ecut = num(inp, 'ecut', 60, 20, 200), calc = pick(inp.dkCalc, ['scf', 'relax', 'vc-relax'], 'scf'), sm = pick(inp.dkSmear, ['auto', 'metal', 'insulator'], 'auto'), checks = [], warnings = [];
  const cell = rep > 1 ? A.supercell(A.buildCell(key), [rep, rep, rep]) : A.buildCell(key), M = cell.cell, nat = cell.atoms.length, Mi = N.inverse(M), vol = Math.abs(N.det(M)), els = [...new Set(cell.atoms.map(q => q.el))], metal = sm === 'metal' || (sm === 'auto' && METALS.includes(key)), mag = MAGNETIC.includes(key);
  const bLen = [0, 1, 2].map(k => 2 * PI * Math.hypot(Mi[0][k], Mi[1][k], Mi[2][k])), vac = key === 'graphene', mesh = sp => bLen.map((b, k) => (vac && k === 2 ? 1 : Math.max(1, Math.ceil(b / sp)))), km = mesh(ks), nk = km[0] * km[1] * km[2];
  const BOHR = 0.529177, npw = E => vol / BOHR ** 3 * E ** 1.5 / (6 * PI * PI), fr = cell.atoms.map(q => [0, 1, 2].map(k => q.x * Mi[0][k] + q.y * Mi[1][k] + q.z * Mi[2][k])), f9 = v => v.toFixed(9).padStart(14), base = slug(key) + (rep > 1 ? `_${rep}x${rep}x${rep}` : '');
  const qe = [`! Quantum ESPRESSO pw.x input for ${def.name} (${def.formula}) — generated starting point, NOT a converged set-up`, '! Replace every pseudopotential file name and re-converge ecutwfc / ecutrho / k-points for YOUR pseudopotentials.', '&CONTROL', `  calculation = '${calc}'`, `  prefix = '${base}'`, "  pseudo_dir = './pseudo'", "  outdir = './tmp'", '  tprnfor = .true.', '  tstress = .true.', '/', '&SYSTEM', '  ibrav = 0', `  nat = ${nat}`, `  ntyp = ${els.length}`, `  ecutwfc = ${F(ecut, 4)}   ! Ry — PLACEHOLDER: use the value recommended for your pseudopotentials, then converge`, `  ecutrho = ${F(8 * ecut, 4)}   ! Ry — 4×ecutwfc for norm-conserving, 8–12× for ultrasoft/PAW`, ...(metal ? ["  occupations = 'smearing'", "  smearing = 'marzari-vanderbilt'", '  degauss = 0.02'] : ["  occupations = 'fixed'"]), ...(mag ? ['  nspin = 2', ...els.map((e, i) => `  starting_magnetization(${i + 1}) = ${['Fe', 'Ni', 'Co', 'Mn'].includes(e) ? '0.5' : '0.0'}`)] : []), '/', '&ELECTRONS', '  conv_thr = 1.0d-8', '  mixing_beta = 0.5', '/', ...(calc !== 'scf' ? ['&IONS', "  ion_dynamics = 'bfgs'", '/'] : []), ...(calc === 'vc-relax' ? ['&CELL', "  cell_dynamics = 'bfgs'", '  press_conv_thr = 0.5', '/'] : []),
    'ATOMIC_SPECIES', ...els.map(e => `  ${e.padEnd(3)} ${String(chem.ELEMENTS[e].mass).padEnd(9)} ${e}.REPLACE_WITH_YOUR_PSEUDOPOTENTIAL.UPF`), 'CELL_PARAMETERS angstrom', ...M.map(v => '  ' + v.map(f9).join(' ')), 'ATOMIC_POSITIONS crystal', ...cell.atoms.map((q, i) => `  ${q.el.padEnd(3)} ${fr[i].map(f9).join(' ')}`), 'K_POINTS automatic', `  ${km.join(' ')}  0 0 0`].join('\n') + '\n';
  const encut = Math.round(ecut * 13.605693 / 10) * 10, incar = [`SYSTEM = ${def.formula} (${def.name})`, '# Generated starting point — re-converge ENCUT and KPOINTS for your POTCAR set', `ENCUT  = ${encut}        # eV, PLACEHOLDER: at least 1.3 × the largest ENMAX in your POTCAR`, 'PREC   = Accurate', 'EDIFF  = 1E-6', 'LREAL  = .FALSE.', ...(metal ? ['ISMEAR = 1', 'SIGMA  = 0.2'] : ['ISMEAR = 0', 'SIGMA  = 0.05']), ...(mag ? ['ISPIN  = 2', `MAGMOM = ${els.map(e => `${cell.atoms.filter(q => q.el === e).length}*${['Fe', 'Ni', 'Co', 'Mn'].includes(e) ? '3.0' : '0.0'}`).join(' ')}`] : ['ISPIN  = 1']), ...(calc === 'scf' ? ['IBRION = -1', 'NSW    = 0'] : ['IBRION = 2', 'NSW    = 100', 'EDIFFG = -0.01', `ISIF   = ${calc === 'vc-relax' ? 3 : 2}`]), 'LWAVE  = .FALSE.', 'LCHARG = .FALSE.'].join('\n') + '\n';
  const kpts = [`Automatic mesh, target spacing ${F(ks, 3)} 1/Angstrom (2*pi included)`, '0', hexLike(M) ? 'Gamma' : 'Monkhorst-Pack', km.join(' '), '0 0 0'].join('\n') + '\n';
  const pp = POT[key] && def.system === 'cubic' ? POT[key] : null;
  const lmp = ['# LAMMPS input for ' + def.name, 'units           metal', 'atom_style      atomic', 'boundary        p p p', `read_data       ${base}.lammps-data`, ...(pp ? [`pair_style      morse ${F(pp.r + 4.5 / pp.al, 4)}`, `pair_coeff      1 1 ${pp.e} ${pp.al} ${pp.r}   # qualitative pair potential; prefer an EAM file`] : ['# PLACEHOLDER: supply an interatomic potential validated for this material', '# (e.g. EAM/MEAM for metals, Tersoff/Stillinger-Weber for covalent solids, Buckingham+Coulomb or ReaxFF for oxides)', 'pair_style      REPLACE_WITH_PAIR_STYLE', 'pair_coeff      * * REPLACE_WITH_POTENTIAL_FILE ' + els.join(' ')]), 'neighbor        2.0 bin', 'thermo          10', 'fix             rel all box/relax aniso 0.0 vmax 0.001', 'min_style       cg', 'minimize        1.0e-10 1.0e-10 5000 50000', `write_data      ${base}_relaxed.lammps-data`].join('\n') + '\n';
  const cuts = [0.6, 0.8, 1, 1.2, 1.4, 1.7].map(f => Math.round(ecut * f)), sps = [0.5, 0.4, 0.3, 0.25, 0.2, 0.15], plan = cuts.map((c, i) => [i + 1, c, 8 * c, Math.round(c * 13.605693), Math.round(npw(c)), mesh(sps[i]).join('×'), sps[i], mesh(sps[i]).reduce((p, q) => p * q, 1)]);
  const kLen = km.map((k, i) => k * 2 * PI / bLen[i]), cost = nk * npw(ecut) * nat * nat;
  checks.push(chk('Fractional coordinates inside the cell', fr.every(f => f.every(v => v > -1e-6 && v < 1 + 1e-6)), `${nat} atoms`), chk('k-mesh meets the target spacing', km.every((k, i) => (vac && i === 2) || bLen[i] / k <= ks * 1.0001), `actual spacings ${km.map((k, i) => F(bLen[i] / k, 3)).join(', ')} Å⁻¹`), chk('Real-space sampling length above 15 Å for metals', !metal || Math.min(...kLen.filter((_, i) => !(vac && i === 2))) > 15, `k × a = ${kLen.map(v => F(v, 3)).join(', ')} Å${metal ? '' : ' (insulator: 10 Å is usually sufficient)'}`), chk('Charge-neutral stoichiometric cell', Number.isInteger(nat / N.sum(Object.values(chem.parseFormula(def.formula).atoms))), `${nat} atoms = ${nat / N.sum(Object.values(chem.parseFormula(def.formula).atoms))} formula units`));
  warnings.push('This application does not run density-functional theory. The decks are syntactically complete starting points: pseudopotential files, cut-offs, smearing, Hubbard U, van der Waals corrections and magnetic order must be chosen and converged by the user.');
  if (mag) warnings.push(`${def.name} is magnetic: a spin-polarised ferromagnetic start is written; antiferro- or ferrimagnetic order (e.g. magnetite) needs explicit per-atom moments.`);
  if (/TiO2|Fe3O4|CeO2|ZnO/.test(key)) warnings.push('Semi-local functionals underestimate the band gap of this oxide and misplace localised d/f states; plan for DFT+U or a hybrid functional if electronic properties matter.');
  if (/graph|hBN|MoS2/.test(key)) warnings.push('Layered material: interlayer binding requires a van der Waals correction (e.g. DFT-D3 or a vdW functional); none is set in the decks.');
  const rows2 = N.linspace(20, 200, 37);
  return {
    summary: `Input decks for ${def.name} (${chem.prettyFormula(def.formula)}, ${nat} atoms, ${els.length} species, cell volume ${F(vol, 5)} Å³): Quantum ESPRESSO pw.x (${calc}), VASP POSCAR/KPOINTS/INCAR and LAMMPS data + input. A k-point spacing of ${F(ks, 3)} Å⁻¹ gives a ${km.join('×')} mesh (${nk} points before symmetry reduction); a ${F(ecut, 4)} Ry cut-off corresponds to about ${F(npw(ecut), 3)} plane waves per band. This application does not run DFT itself — the files are for external solvers.`,
    kpis: [
      kpi('nat', 'Atoms in the cell', nat, '–', `Atoms in the ${rep > 1 ? `${rep}×${rep}×${rep} supercell` : 'conventional cell'} written to every deck (${els.map(e => `${e} ${cell.atoms.filter(q => q.el === e).length}`).join(', ')}).`),
      kpi('kmesh', 'k-point mesh', km.join('×'), '', `Smallest mesh whose spacing |b_i|/n_i does not exceed ${F(ks, 3)} Å⁻¹ along each reciprocal axis${vac ? ' (one point along the vacuum direction)' : ''}.`),
      kpi('nk', 'k-points (full grid)', nk, '–', 'Number of mesh points before symmetry reduction; the solver reduces it by the point-group order (up to 48 for cubic crystals).'),
      kpi('npw', 'Plane waves per band', npw(ecut), '–', `N_pw ≈ Ω k_cut³/(6π²) for a ${F(ecut, 4)} Ry wave-function cut-off; memory and time per band scale linearly with it.`),
      kpi('encut', 'Equivalent VASP cut-off', encut, 'eV', 'The Rydberg cut-off converted to eV and rounded to 10 eV; a placeholder to be replaced by at least 1.3 × the largest ENMAX of the chosen POTCAR files.'),
      kpi('vol', 'Cell volume', vol, 'Å³', 'Volume of the simulation cell; plane-wave count is proportional to it, so vacuum in slab or sheet cells is not free.'),
      kpi('cost', 'Relative cost index', cost, '–', 'N_k × N_pw × N_atoms²: a scaling guide (not a time prediction) for comparing set-ups; doubling the supercell edge multiplies it by about 64 at fixed k-density reduction.'),
      kpi('smear', 'Occupation scheme', metal ? 'smearing (metal)' : 'fixed (insulator)', '', metal ? 'Marzari–Vanderbilt cold smearing 0.02 Ry (QE) / Methfessel–Paxton 0.2 eV (VASP): needed for Fermi-surface sampling; check that the entropy term stays below 1 meV/atom.' : 'Fixed occupations (QE) / Gaussian smearing 0.05 eV (VASP), appropriate for a gapped material.'),
    ],
    series: [
      { key: 'npw_E', title: 'Basis size vs plane-wave cut-off', xlabel: 'ecutwfc (Ry)', ylabel: 'plane waves per band', lines: [{ name: 'N_pw', x: rows2, y: rows2.map(npw) }], explain: 'The number of plane waves grows as E_cut^(3/2); together with the cut-off ladder in the table it shows the price of each convergence step.' },
      { key: 'nk_s', title: 'k-points vs target spacing', xlabel: 'k-point spacing (1/Å)', ylabel: 'k-points in the full grid', logy: true, lines: [{ name: 'N_k', x: N.linspace(0.08, 0.8, 37), y: N.linspace(0.08, 0.8, 37).map(s => mesh(s).reduce((p, q) => p * q, 1)) }], explain: 'Mesh size required for a given reciprocal-space resolution. Metals typically need 0.1–0.2 Å⁻¹, insulators 0.2–0.4 Å⁻¹.' },
    ],
    tables: [{ key: 'plan', title: 'Convergence-test plan', columns: ['step', 'ecutwfc (Ry)', 'ecutrho (Ry)', 'ENCUT (eV)', 'plane waves', 'k-mesh', 'spacing (1/Å)', 'k-points'], rows: plan, explain: 'Run the cut-off ladder at a fixed moderate k-mesh, then the k-mesh ladder at the converged cut-off. Accept the first setting beyond which total energy changes by less than 1 meV/atom, forces by less than 0.01 eV/Å and stress by less than 0.5 kbar.' }],
    equations: [
      eqn('formula', 'Composition', `${chem.prettyFormula(def.formula)},  M = ${F(chem.molarMass(def.formula), 6)} g/mol,  ${def.spacegroup}`), eqn('governing', 'Kohn–Sham equations (solved by the external code)', '[−(ħ²/2m)∇² + v_ext(r) + v_H[n](r) + v_xc[n](r)] ψ_i(r) = ε_i ψ_i(r),   n(r) = Σ_i f_i |ψ_i(r)|²'),
      eqn('formula', 'Plane-wave basis size', 'N_pw ≈ Ω (2m E_cut/ħ²)^(3/2) / (6π²)'), eqn('formula', 'k-mesh from spacing', 'n_i = ⌈ |b_i| / Δk ⌉,  b_i = 2π (a_j × a_k)/Ω'),
    ],
    geometry: { kind: 'atoms', atoms: cell.atoms.map(q => ({ el: q.el, x: q.x, y: q.y, z: q.z })), cell: M, title: cell.title },
    checks, warnings,
    recommendations: [
      rec('high', 'Before any production run complete the convergence plan in the table and record the results: cut-off and k-mesh convergence must be demonstrated for your pseudopotentials, not assumed.', 'Verification requirement for first-principles calculations.'),
      rec('medium', `Start from the experimental structure exported here, relax ions first (${calc === 'scf' ? 'switch the calculation type to relax' : 'as configured'}), then the cell; compare the relaxed lattice parameters with experiment (PBE typically overestimates by 1–2 %).`, 'Standard structural-relaxation workflow.'),
      rec('medium', metal ? 'For this metal converge the k-mesh together with the smearing width; use the tetrahedron method for the final density of states.' : 'For this gapped material a coarser mesh suffices for total energies, but band gaps need a hybrid functional or GW correction.', 'Occupation scheme selected for the deck.'),
      rec('low', `Cost scales as about N_atoms³: use the smallest cell and the converged — not the maximum — cut-off and mesh. Halving an over-converged k-mesh in each direction saves roughly ${nk > 1 ? '8×' : 'little'} compute energy with no loss of accuracy.`, 'Cost index and convergence ladder.', true),
    ],
    publish: { 'material.formula': def.formula, 'sim.n_atoms': nat, 'sim.n_kpoints': nk, 'sim.ecut_Ry': ecut, 'prop.lattice_a': A.cellParams(A.buildCell(key).cell).a * 1e-10 },
    exports: [
      { name: `${base}.pw.in`, mime: 'text/plain', content: qe, description: 'Quantum ESPRESSO pw.x input (ibrav = 0, explicit cell and crystal coordinates); pseudopotential names and cut-offs are marked placeholders.' },
      { name: `${base}.POSCAR`, mime: 'text/plain', content: A.toPOSCAR(cell), description: 'VASP POSCAR (direct coordinates, species order as in the file).' },
      { name: `${base}.KPOINTS`, mime: 'text/plain', content: kpts, description: `VASP KPOINTS, ${km.join('×')} automatic mesh.` },
      { name: `${base}.INCAR`, mime: 'text/plain', content: incar, description: 'VASP INCAR starting point (ENCUT is a placeholder; POTCAR must be assembled by the user under their licence).' },
      { name: `${base}.lammps-data`, mime: 'text/plain', content: A.toLAMMPSData(cell), description: 'LAMMPS data file (atom_style atomic).' },
      { name: `in.${base}.lammps`, mime: 'text/plain', content: lmp, description: pp ? 'LAMMPS relaxation input with a qualitative Morse pair potential.' : 'LAMMPS relaxation input; the pair style and potential file are marked placeholders.' },
      { name: `${base}.cif`, mime: 'chemical/x-cif', content: A.toCIF(cell), description: 'CIF of the same cell for visual checking.' },
    ],
    numerics: { method: 'Deterministic input-deck generation from the crystal library (no electronic-structure calculation performed)', dof: 3 * nat, steps: 1, iterations: 1, residual: 0, notes: 'k-mesh from reciprocal-lattice lengths; plane-wave count from the free-electron sphere estimate.' },
  };
}
function hexLike(M) { const p = A.cellParams(M); return Math.abs(p.gamma - 120) < 0.5 || Math.abs(p.gamma - 60) < 0.5; }

const Mo = (...m) => ({ mode: m });
const potOpts = Object.keys(POT).map(k => ({ value: k, label: POT[k].name }));
const inputs = [
  { key: 'mode', label: 'Simulation', type: 'select', value: 'tb', group: 'Model', help: 'Selects the atomistic method; each mode shows only its own inputs.', options: [{ value: 'md', label: 'Molecular dynamics (pair potential)' }, { value: 'mc', label: 'Metropolis Monte Carlo (Ising / lattice gas)' }, { value: 'kmc', label: 'Kinetic Monte Carlo film growth' }, { value: 'tb', label: 'Electronic structure (model Hamiltonians)' }, { value: 'phonon', label: 'Lattice dynamics (phonons)' }, { value: 'neb', label: 'Reaction / diffusion pathway (NEB, TST)' }, { value: 'deck', label: 'Input decks for DFT and MD codes' }] },
  { key: 'mdSys', label: 'System', type: 'select', value: 'bulk', group: 'Molecular dynamics', showIf: Mo('md'), options: [{ value: 'bulk', label: 'Periodic bulk cell' }, { value: 'cluster', label: 'Free nanocluster' }], help: 'Periodic box of fcc cells (bulk solid or liquid) or a free spherical cluster in vacuum.' },
  { key: 'mdPot', label: 'Potential', type: 'select', value: 'Ar', group: 'Potential', showIf: Mo('md', 'phonon'), options: potOpts, help: 'Lennard-Jones argon or a Morse potential for an fcc metal (published lattice-property-fitted parameters).' },
  { key: 'mdEns', label: 'Production ensemble', type: 'select', value: 'nve', group: 'Molecular dynamics', showIf: Mo('md'), options: [{ value: 'nve', label: 'NVE (microcanonical)' }, { value: 'berendsen', label: 'NVT, Berendsen thermostat' }, { value: 'langevin', label: 'NVT, Langevin thermostat' }], help: 'Equilibration always uses a Berendsen thermostat; this sets the ensemble of the production phase.' },
  { key: 'ncell', label: 'fcc cells per box edge', unit: '–', value: 3, min: 2, max: 5, step: 1, group: 'Molecular dynamics', showIf: Mo('md'), sens: false, help: 'The periodic box holds 4n³ atoms (32 to 500).' },
  { key: 'dens', label: 'Density relative to the static crystal', unit: '–', value: 0.75, min: 0.4, max: 1.15, group: 'Molecular dynamics', showIf: Mo('md'), help: 'Number density divided by that of the equilibrium fcc crystal: about 0.75 for liquid argon near its triple point, 1.0 for the solid.' },
  { key: 'mdD', label: 'Cluster diameter', unit: 'm', value: 1.8e-9, min: 8e-10, max: 4e-9, log: true, group: 'Molecular dynamics', bus: 'geom.diameter', showIf: Mo('md'), help: 'Diameter of the free cluster carved from the fcc lattice (limited to 300 atoms).' },
  { key: 'T', label: 'Temperature', unit: 'K', value: 94.4, min: 1, max: 3000, group: 'Conditions', showIf: Mo('md', 'phonon'), help: 'Target temperature of the molecular-dynamics run, or the temperature at which the harmonic heat capacity is evaluated.' },
  { key: 'dt', label: 'Time step', unit: 'fs', value: 10, min: 0.2, max: 20, group: 'Molecular dynamics', showIf: Mo('md'), help: 'Integration step: about 10 fs for argon, 1–3 fs for metals. Too large a step shows up as energy drift.' },
  { key: 'tau', label: 'Thermostat time constant', unit: 'fs', value: 300, min: 10, max: 5000, group: 'Molecular dynamics', showIf: Mo('md'), help: 'Relaxation time of the Berendsen or Langevin thermostat in production.' },
  { key: 'nEq', label: 'Equilibration steps', unit: '–', value: 300, min: 50, max: 5000, step: 50, group: 'Molecular dynamics', showIf: Mo('md'), sens: false, help: 'Thermostatted steps discarded before measurement.' },
  { key: 'nProd', label: 'Production steps', unit: '–', value: 1000, min: 100, max: 10000, step: 100, group: 'Molecular dynamics', showIf: Mo('md'), sens: false, help: 'Steps over which averages, RDF and MSD are collected.' },
  { key: 'seed', label: 'Random seed', unit: '–', value: 1, min: 1, max: 99999, step: 1, group: 'Model', showIf: Mo('md', 'mc', 'kmc'), sens: false, help: 'Seed of the random-number generator (velocities, Monte Carlo moves); identical seeds reproduce identical runs.' },
  { key: 'mcModel', label: 'Monte Carlo model', type: 'select', value: 'ising', group: 'Monte Carlo', showIf: Mo('mc'), options: [{ value: 'ising', label: '2-D Ising model (magnetic ordering)' }, { value: 'latticegas', label: 'Lattice-gas adsorption isotherm' }], help: 'Spin model with the exact Onsager benchmark, or grand-canonical adsorption with lateral interactions.' },
  { key: 'mcL', label: 'Lattice size L', unit: 'sites', value: 16, min: 8, max: 48, step: 1, group: 'Monte Carlo', showIf: Mo('mc'), sens: false, help: 'The periodic square lattice has L × L sites.' },
  { key: 'mcSweeps', label: 'Measurement sweeps', unit: '–', value: 600, min: 100, max: 4000, step: 100, group: 'Monte Carlo', showIf: Mo('mc'), sens: false, help: 'Monte Carlo sweeps (one attempted move per site) per temperature or pressure point; 40 % as many are used for equilibration.' },
  { key: 'Jmev', label: 'Exchange constant J', unit: 'meV', value: 10, min: 0.1, max: 500, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['ising'] }, help: 'Nearest-neighbour coupling; it converts the reduced critical temperature into kelvin.' },
  { key: 'mcT', label: 'Reduced temperature k_BT/J', unit: '–', value: 2.0, min: 0.5, max: 5, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['ising'] }, help: 'Temperature of the diagnostic run and the displayed spin configuration (the critical value is 2.269).' },
  { key: 'mcH', label: 'Reduced field h/J', unit: '–', value: 0, min: -1, max: 1, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['ising'] }, help: 'External magnetic field in units of J.' },
  { key: 'Tads', label: 'Adsorption temperature', unit: 'K', value: 300, min: 20, max: 2000, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['latticegas'] }, help: 'Temperature of surface and gas in the lattice-gas model.' },
  { key: 'Eads', label: 'Adsorption energy', unit: 'eV', value: 0.5, min: 0.01, max: 3, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['latticegas'] }, help: 'Binding energy of an isolated adsorbate (positive = bound).' },
  { key: 'wmev', label: 'Lateral interaction w', unit: 'meV', value: 20, min: -200, max: 200, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['latticegas'] }, help: 'Nearest-neighbour adsorbate interaction: positive attractive, negative repulsive.' },
  { key: 'mAds', label: 'Adsorbate mass', unit: 'u', value: 28, min: 1, max: 400, group: 'Monte Carlo', showIf: { mode: ['mc'], mcModel: ['latticegas'] }, help: 'Molecular mass for the ideal-gas translational partition function (28 u for CO or N₂).' },
  { key: 'Pgas', label: 'Gas pressure', unit: 'Pa', value: 1e3, min: 1e-10, max: 1e8, log: true, group: 'Monte Carlo', bus: 'proc.P', showIf: { mode: ['mc'], mcModel: ['latticegas'] }, help: 'Pressure at which coverage, acceptance and the displayed configuration are evaluated.' },
  { key: 'T_C', label: 'Process temperature', unit: '°C', value: 227, min: -150, max: 1500, group: 'Conditions', bus: 'proc.T_C', showIf: Mo('kmc', 'neb'), help: 'Substrate temperature during growth, or the temperature at which the transition-state rate is evaluated.' },
  { key: 'kL', label: 'Lattice size L', unit: 'sites', value: 32, min: 16, max: 48, step: 1, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), sens: false, help: 'The periodic square substrate has L × L columns.' },
  { key: 'Es', label: 'Surface diffusion barrier E_s', unit: 'eV', value: 0.75, min: 0.05, max: 2.5, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Activation energy for a free adatom hop on a terrace.' },
  { key: 'En', label: 'Lateral bond energy E_n', unit: 'eV', value: 0.25, min: 0, max: 1, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Additional barrier per in-plane nearest neighbour; it controls detachment from island edges.' },
  { key: 'Ees', label: 'Step-edge (Ehrlich–Schwoebel) barrier', unit: 'eV', value: 0, min: 0, max: 0.5, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Extra barrier for hopping down a step; positive values hinder interlayer transport and cause mounds.' },
  { key: 'flux', label: 'Deposition flux', unit: 'ML/s', value: 1, min: 1e-4, max: 1e3, log: true, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Arrival rate of atoms in monolayers per second.' },
  { key: 'kCov', label: 'Deposited thickness', unit: 'ML', value: 2, min: 0.2, max: 5, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Total coverage to deposit (subject to the event budget).' },
  { key: 'nu0', label: 'Attempt frequency', unit: '1/s', value: 1e13, min: 1e10, max: 1e15, log: true, group: 'Kinetic Monte Carlo', showIf: Mo('kmc'), help: 'Pre-exponential factor of the hop rate.' },
  { key: 'aLat', label: 'Surface lattice spacing', unit: 'm', value: 2.9e-10, min: 1e-10, max: 1e-9, group: 'Kinetic Monte Carlo', bus: 'prop.lattice_a', showIf: Mo('kmc'), help: 'Hop distance and layer height used to convert lattice units into metres.' },
  { key: 'tbModel', label: 'Electronic model', type: 'select', value: 'graphene', group: 'Electronic structure', showIf: Mo('tb'), options: [{ value: 'graphene', label: 'Graphene π bands (tight binding)' }, { value: 'cnt', label: 'Carbon nanotube (n,m), zone folding' }, { value: 'ssh', label: '1-D chain / SSH dimerised chain' }, { value: 'brus', label: 'Quantum dot (Brus equation)' }, { value: 'well', label: 'Finite quantum well (numerical Schrödinger)' }], help: 'Model Hamiltonian to solve.' },
  { key: 'gamma0', label: 'Hopping integral γ₀', unit: 'eV', value: 2.7, min: 1, max: 4, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['graphene', 'cnt'] }, help: 'Nearest-neighbour π hopping energy; empirical values lie between 2.7 and 3.0 eV.' },
  { key: 'tbS', label: 'Overlap integral s', unit: '–', value: 0, min: 0, max: 0.2, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['graphene', 'cnt'] }, help: 'Nearest-neighbour orbital overlap (0 for the symmetric model, about 0.13 for an asymmetric π/π* fit).' },
  { key: 'tn', label: 'Chiral index n', unit: '–', value: 10, min: 1, max: 20, step: 1, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['cnt'] }, sens: false, help: 'First chiral index of the nanotube.' },
  { key: 'tm', label: 'Chiral index m', unit: '–', value: 10, min: 0, max: 20, step: 1, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['cnt'] }, sens: false, help: 'Second chiral index (m ≤ n; swapped automatically otherwise).' },
  { key: 't1', label: 'Hopping t₁ (intra-cell)', unit: 'eV', value: 2.4, min: 0.05, max: 5, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['ssh'] }, help: 'Hopping on the first bond of each cell; set t₁ = t₂ for a uniform chain.' },
  { key: 't2', label: 'Hopping t₂ (inter-cell)', unit: 'eV', value: 3.0, min: 0.05, max: 5, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['ssh'] }, help: 'Hopping on the second bond; t₂ > t₁ gives the topological phase for a chain ending on t₁ bonds.' },
  { key: 'aCh', label: 'Chain period', unit: 'Å', value: 2.8, min: 1, max: 10, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['ssh'] }, help: 'Length of the two-site unit cell.' },
  { key: 'qdMat', label: 'Quantum-dot material', type: 'select', value: 'CdSe', group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['brus'] }, options: Object.keys(QD).map(k => ({ value: k, label: `${chem.prettyFormula(QD[k].f)} (E_g ${QD[k].Eg} eV)` })), help: 'Semiconductor with its bulk gap, effective masses and permittivity.' },
  { key: 'qdD', label: 'Dot diameter', unit: 'm', value: 4e-9, min: 1e-9, max: 2e-8, log: true, group: 'Electronic structure', bus: 'geom.diameter', showIf: { mode: ['tb'], tbModel: ['brus'] }, help: 'Diameter of the spherical nanocrystal.' },
  { key: 'qdV0', label: 'Confining barrier', unit: 'eV', value: 4, min: 0.5, max: 10, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['brus'] }, help: 'Band offset to the surrounding ligand shell or matrix for the finite-barrier estimate (several eV for organic ligands, about 1 eV for a semiconductor shell).' },
  { key: 'wW', label: 'Well width', unit: 'm', value: 6e-9, min: 5e-10, max: 3e-8, log: true, group: 'Electronic structure', bus: 'geom.thickness', showIf: { mode: ['tb'], tbModel: ['well'] }, help: 'Thickness of the quantum-well layer.' },
  { key: 'wV0', label: 'Well depth (band offset)', unit: 'eV', value: 0.3, min: 0.01, max: 5, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['well'] }, help: 'Barrier height; about 0.3 eV for the conduction band of GaAs in Al₀.₃Ga₀.₇As.' },
  { key: 'wMw', label: 'Effective mass in the well', unit: 'm₀', value: 0.067, min: 0.01, max: 2, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['well'] }, help: 'Carrier effective mass in the well material (0.067 for electrons in GaAs).' },
  { key: 'wMb', label: 'Effective mass in the barrier', unit: 'm₀', value: 0.092, min: 0.01, max: 2, group: 'Electronic structure', showIf: { mode: ['tb'], tbModel: ['well'] }, help: 'Carrier effective mass in the barrier material.' },
  { key: 'phModel', label: 'Lattice-dynamics model', type: 'select', value: 'diatomic', group: 'Lattice dynamics', showIf: Mo('phonon'), options: [{ value: 'mono', label: '1-D monatomic chain' }, { value: 'diatomic', label: '1-D diatomic chain' }, { value: 'fcc', label: '3-D fcc crystal from the pair potential' }], help: 'Analytic chain models or the full dynamical matrix of an fcc crystal with the selected potential.' },
  { key: 'phM1', label: 'Mass 1', unit: 'u', value: 65.38, min: 1, max: 250, group: 'Lattice dynamics', showIf: { mode: ['phonon'], phModel: ['mono', 'diatomic'] }, help: 'Atomic mass of the (first) atom of the chain.' },
  { key: 'phM2', label: 'Mass 2', unit: 'u', value: 15.999, min: 1, max: 250, group: 'Lattice dynamics', showIf: { mode: ['phonon'], phModel: ['diatomic'] }, help: 'Atomic mass of the second atom of the diatomic chain.' },
  { key: 'phK', label: 'Spring constant', unit: 'N/m', value: 50, min: 0.5, max: 500, group: 'Lattice dynamics', showIf: { mode: ['phonon'], phModel: ['mono', 'diatomic'] }, help: 'Nearest-neighbour force constant (tens of N/m for covalent and ionic bonds).' },
  { key: 'phA', label: 'Chain period', unit: 'Å', value: 3.0, min: 1, max: 10, group: 'Lattice dynamics', showIf: { mode: ['phonon'], phModel: ['mono', 'diatomic'] }, help: 'Repeat distance of the chain (one atom per period for the monatomic, two for the diatomic chain).' },
  { key: 'pes', label: 'Potential-energy surface', type: 'select', value: 'mb', group: 'Pathway', showIf: Mo('neb'), options: [{ value: 'mb', label: 'Müller–Brown model surface' }, { value: 'hop', label: 'Adatom hop on a square lattice' }], help: 'Model two-dimensional surface on which the minimum-energy path is searched.' },
  { key: 'nebPath', label: 'End states (Müller–Brown)', type: 'select', value: 'AC', group: 'Pathway', showIf: { mode: ['neb'], pes: ['mb'] }, options: [{ value: 'AC', label: 'Deep minimum A → minimum C (via intermediate B)' }, { value: 'AB', label: 'Minimum A → intermediate B' }, { value: 'BC', label: 'Intermediate B → minimum C' }], help: 'Initial and final minima of the path.' },
  { key: 'nImg', label: 'Number of images', unit: '–', value: 15, min: 7, max: 25, step: 1, group: 'Pathway', showIf: Mo('neb'), sens: false, help: 'Images along the band including both end points.' },
  { key: 'pesE', label: 'Energy unit of the surface', unit: 'eV', value: 0.005, min: 1e-4, max: 0.05, log: true, group: 'Pathway', showIf: Mo('neb'), help: 'Converts the dimensionless model surface into electron-volts (0.005 eV per unit makes the main Müller–Brown barrier about 0.53 eV).' },
  { key: 'pesL', label: 'Length unit of the surface', unit: 'Å', value: 1, min: 0.2, max: 5, group: 'Pathway', showIf: Mo('neb'), help: 'Converts the dimensionless coordinates into ångström; it sets the jump length and the vibrational frequencies.' },
  { key: 'nebMass', label: 'Mass of the moving atom', unit: 'u', value: 63.546, min: 1, max: 250, group: 'Pathway', showIf: Mo('neb'), help: 'Mass used for the Vineyard attempt frequency and the tunnelling correction.' },
  { key: 'nebNu', label: 'Reference attempt frequency', unit: 'Hz', value: 1e13, min: 1e10, max: 1e15, log: true, group: 'Pathway', showIf: Mo('neb'), help: 'Conventional prefactor for comparison with the computed Vineyard value.' },
  { key: 'dim', label: 'Diffusion dimensionality', unit: '–', value: 2, min: 1, max: 3, step: 1, group: 'Pathway', showIf: Mo('neb'), sens: false, help: 'Dimension d in D = kℓ²/(2d): 1 channel, 2 surface, 3 bulk.' },
  { key: 'dkCrystal', label: 'Crystal', type: 'select', value: 'ZnO_wurtzite', group: 'Input decks', showIf: Mo('deck'), options: KEYS.map(k => ({ value: k, label: `${A.CRYSTALS[k].name} — ${chem.prettyFormula(A.CRYSTALS[k].formula)}` })), help: 'Library structure written to the Quantum ESPRESSO, VASP and LAMMPS files.' },
  { key: 'dkRep', label: 'Supercell repeats', unit: '–', value: 1, min: 1, max: 3, step: 1, group: 'Input decks', showIf: Mo('deck'), sens: false, help: 'Repeats of the conventional cell along each axis (1 = unit cell).' },
  { key: 'kspacing', label: 'k-point spacing', unit: '1/Å', value: 0.25, min: 0.05, max: 1, group: 'Input decks', showIf: Mo('deck'), help: 'Target reciprocal-space resolution (2π included): 0.1–0.2 for metals, 0.2–0.4 for insulators.' },
  { key: 'ecut', label: 'Wave-function cut-off', unit: 'Ry', value: 60, min: 20, max: 200, group: 'Input decks', showIf: Mo('deck'), help: 'Starting plane-wave cut-off — a placeholder to be converged for the chosen pseudopotentials.' },
  { key: 'dkCalc', label: 'Calculation type', type: 'select', value: 'scf', group: 'Input decks', showIf: Mo('deck'), options: [{ value: 'scf', label: 'Single-point (scf)' }, { value: 'relax', label: 'Ionic relaxation' }, { value: 'vc-relax', label: 'Variable-cell relaxation' }], help: 'Task written to the DFT decks.' },
  { key: 'dkSmear', label: 'Occupations', type: 'select', value: 'auto', group: 'Input decks', showIf: Mo('deck'), options: [{ value: 'auto', label: 'Automatic (metal or insulator from the library)' }, { value: 'metal', label: 'Metal: smearing' }, { value: 'insulator', label: 'Insulator: fixed occupations' }], help: 'Occupation / smearing scheme for the electronic states.' },
];

const RUN = { md: runMD, mc: inp => (pick(inp.mcModel, ['ising', 'latticegas'], 'ising') === 'ising' ? runIsing(inp) : runLatticeGas(inp)), kmc: runKMC, tb: runTB, phonon: runPhonon, neb: runNEB, deck: runDeck };

export default {
  id: 17,
  tagline: 'Small real atomistic simulations in the browser — molecular dynamics, Monte Carlo, kinetic Monte Carlo growth, tight-binding and quantum-confinement models, phonons and reaction pathways — plus ready-to-run input decks for external DFT and MD codes.',
  inputs,
  presets: [
    { name: 'Liquid argon near the triple point (Lennard-Jones MD)', note: 'Classic benchmark: 94.4 K, 1.37 g/cm³, D ≈ 2.4 × 10⁻⁹ m²/s.', values: { mode: 'md', mdSys: 'bulk', mdPot: 'Ar', mdEns: 'nve', ncell: 3, dens: 0.75, T: 94.4, dt: 10, nEq: 300, nProd: 1000 } },
    { name: 'Gold nanocluster heating (Morse MD, 1.6 nm)', note: 'Free Au cluster with a Langevin bath; the Lindemann index indicates surface pre-melting.', values: { mode: 'md', mdSys: 'cluster', mdPot: 'Au', mdEns: 'langevin', mdD: 1.6e-9, T: 600, dt: 3, tau: 500, nEq: 300, nProd: 1000 } },
    { name: 'Graphene π bands and Dirac point', note: 'Nearest-neighbour tight binding with γ₀ = 2.7 eV.', values: { mode: 'tb', tbModel: 'graphene', gamma0: 2.7, tbS: 0 } },
    { name: '(10,10) armchair carbon nanotube', note: 'Metallic tube from zone folding; try (10,0) for a semiconducting gap near 1 eV.', values: { mode: 'tb', tbModel: 'cnt', tn: 10, tm: 10, gamma0: 2.7 } },
    { name: 'CdSe quantum dot, 4 nm', note: 'Size-tunable emission: Brus and finite-barrier estimates against the empirical sizing curve.', values: { mode: 'tb', tbModel: 'brus', qdMat: 'CdSe', qdD: 4e-9, qdV0: 4 } },
    { name: 'ZnO-like diatomic chain phonons', note: 'Zn and O masses: acoustic and optical branches with a wide gap.', values: { mode: 'phonon', phModel: 'diatomic', phM1: 65.38, phM2: 15.999, phK: 50, phA: 3.0, T: 300 } },
    { name: 'Epitaxial film growth, 500 K, 1 ML/s (KMC)', note: 'Solid-on-solid growth with bond counting: roughness and growth oscillations.', values: { mode: 'kmc', T_C: 227, Es: 0.75, En: 0.25, flux: 1, kCov: 2, kL: 32 } },
    { name: 'TiO₂ anatase DFT input decks', note: 'Quantum ESPRESSO, VASP and LAMMPS starting files with a convergence plan.', values: { mode: 'deck', dkCrystal: 'TiO2_anatase', dkRep: 1, kspacing: 0.25, ecut: 60, dkCalc: 'relax' } },
    { name: 'Müller–Brown reaction path (NEB + TST)', note: 'Benchmark surface: saddle at (−0.822, 0.624), E‡ = −40.66.', values: { mode: 'neb', pes: 'mb', nebPath: 'AC', nImg: 15, T_C: 300 } },
    { name: '2-D Ising model and the Onsager transition', note: 'Monte Carlo against the exact critical temperature 2.269 J/k_B.', values: { mode: 'mc', mcModel: 'ising', mcL: 16, mcT: 2.0, mcSweeps: 600 } },
  ],
  run(inp) { return RUN[pick(inp.mode, Object.keys(RUN), 'tb')](inp); },
};
