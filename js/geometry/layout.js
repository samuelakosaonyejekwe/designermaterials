// Lithography / metasurface layout engine: resonator primitives, periodic arrays, phase-gradient layouts and
// layout statistics. Coordinates in nm; polygons are simple (no holes), counter-clockwise, on integer layers.
const { sin, cos, sqrt, abs, min, max, PI, round, hypot } = Math;
const num = (v, d) => Number.isFinite(Number(v)) && v !== null && v !== '' ? Number(v) : d;

export const polyArea = pts => { let a = 0; for (let i = 0, n = pts.length; i < n; i++) { const p = pts[i], q = pts[(i + 1) % n]; a += p[0] * q[1] - q[0] * p[1]; } return a / 2; };
const clean = pts => { const o = []; for (const p of pts) { const l = o[o.length - 1]; if (!l || hypot(p[0] - l[0], p[1] - l[1]) > 1e-9) o.push([p[0], p[1]]); } if (o.length > 1 && hypot(o[0][0] - o[o.length - 1][0], o[0][1] - o[o.length - 1][1]) < 1e-9) o.pop(); return polyArea(o) < 0 ? o.reverse() : o; };
const rect = (x0, y0, x1, y1) => [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
const circle = (r, n = 48, cx = 0, cy = 0) => Array.from({ length: n }, (_, i) => [cx + r * cos(2 * PI * i / n), cy + r * sin(2 * PI * i / n)]);
const rot = (pts, a, cx = 0, cy = 0) => { const c = cos(a), s = sin(a); return pts.map(([x, y]) => [cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c]); };
const crossPoly = (L, w) => { const a = L / 2, b = w / 2; return [[b, b], [a, b], [a, -b], [b, -b], [b, -a], [-b, -a], [-b, -b], [-a, -b], [-a, b], [-b, b], [-b, a], [b, a]]; };
const cRing = (a, w, g) => [[g / 2, a / 2], [a / 2, a / 2], [a / 2, -a / 2], [-a / 2, -a / 2], [-a / 2, a / 2], [-g / 2, a / 2], [-g / 2, a / 2 - w], [-a / 2 + w, a / 2 - w], [-a / 2 + w, -a / 2 + w], [a / 2 - w, -a / 2 + w], [a / 2 - w, a / 2 - w], [g / 2, a / 2 - w]];
const cArc = (ro, w, g, n = 40) => { const ri = ro - w, h = Math.asin(min(0.95, g / 2 / max(1e-9, ro - w / 2))), a0 = PI / 2 + h, span = 2 * PI - 2 * h, o = [], i = []; for (let k = 0; k <= n; k++) { const t = a0 + span * k / n; o.push([ro * cos(t), ro * sin(t)]); i.push([ri * cos(t), ri * sin(t)]); } return [...o, ...i.reverse()]; };

// Stroke an open polyline into a simple outline polygon (mitre joins limited to 3× half-width; flat or extended ends).
export function pathToPolygon(pts, width, { extend = 0 } = {}) {
  const P = pts.filter((p, i) => i === 0 || hypot(p[0] - pts[i - 1][0], p[1] - pts[i - 1][1]) > 1e-12), n = P.length, hw = width / 2;
  if (n < 2) return [];
  const dirs = []; for (let i = 0; i < n - 1; i++) { const d = [P[i + 1][0] - P[i][0], P[i + 1][1] - P[i][1]], l = hypot(d[0], d[1]); dirs.push([d[0] / l, d[1] / l]); }
  const L = [], R = [];
  for (let i = 0; i < n; i++) {
    const d0 = dirs[max(0, i - 1)], d1 = dirs[min(n - 2, i)]; let p = P[i];
    if (i === 0) p = [p[0] - d1[0] * extend, p[1] - d1[1] * extend]; if (i === n - 1) p = [p[0] + d0[0] * extend, p[1] + d0[1] * extend];
    let mx = -(d0[1] + d1[1]), my = d0[0] + d1[0]; const ml = hypot(mx, my) || 1; mx /= ml; my /= ml;
    const k = min(3, 1 / max(1e-6, mx * -d1[1] + my * d1[0]));
    L.push([p[0] + mx * hw * k, p[1] + my * hw * k]); R.push([p[0] - mx * hw * k, p[1] - my * hw * k]);
  }
  const out = [...L, ...R.reverse()];
  return polyArea(out) < 0 ? out.reverse() : out;
}

export const RESONATOR_DEFAULTS = {
  srr: { size: 400, width: 60, gap: 60, shape: 'square' }, doubleSrr: { size: 400, width: 50, gap: 50, spacing: 30, shape: 'square' },
  elc: { size: 400, width: 40, gap: 30, plate: 160 }, cross: { length: 300, width: 60 }, patch: { width: 300, height: 300 },
  fishnet: { px: 600, py: 600, holeX: 300, holeY: 300 }, nanodisc: { diameter: 200, segments: 48 }, nanopillar: { diameter: 200, segments: 48, shape: 'circle', height: 600 },
  nanobar: { length: 250, width: 90, rotation_deg: 0 }, vAntenna: { armLength: 200, width: 40, angle_deg: 90, rotation_deg: 0 }, cWire: { length: 300, width: 60, pair: false, gap: 40 },
  jerusalemCross: { length: 300, width: 40, capLength: 140, capWidth: 40 }, helmholtzNeck2D: { cavityWidth: 300, cavityHeight: 300, wall: 30, neckWidth: 60, neckLength: 80 },
};

// One resonator centred on the origin. Returns an array of polygons [{layer, pts}] (nm).
export function resonator({ type = 'srr', params = {}, layer = 1 } = {}) {
  const D = RESONATOR_DEFAULTS[type];
  if (!D) throw new Error(`Unknown resonator type "${type}". Available: ${Object.keys(RESONATOR_DEFAULTS).join(', ')}`);
  const p = {}; for (const k in D) p[k] = typeof D[k] === 'number' ? max(k.endsWith('_deg') ? -1e9 : 1e-6, num(params[k], D[k])) : (params[k] ?? D[k]);
  let out;
  if (type === 'srr') { const w = min(p.width, p.size * 0.45), g = min(p.gap, p.size - 2 * w); out = [p.shape === 'circular' ? cArc(p.size / 2, w, g) : cRing(p.size, w, g)]; }
  else if (type === 'doubleSrr') {
    const w = min(p.width, p.size * 0.2), a2 = p.size - 2 * w - 2 * p.spacing, g = min(p.gap, max(1, a2 - 2 * w));
    if (a2 <= 2 * w) throw new Error('doubleSrr: size too small for two rings with this width and spacing');
    out = p.shape === 'circular' ? [cArc(p.size / 2, w, g), rot(cArc(a2 / 2, w, g), PI)] : [cRing(p.size, w, g), rot(cRing(a2, w, g), PI)];
  } else if (type === 'elc') {
    const a = p.size / 2, w = min(p.width, a * 0.4), g = min(p.gap, p.size - 4 * w), pl = min(p.plate, p.size - 2 * w - 2) / 2, tw = min(w, (p.size - 2 * w - g) / 4);
    out = [rect(-a, -a, -a + w, a), rect(a - w, -a, a, a), rect(-a + w, a - w, a - w, a), rect(-a + w, -a, a - w, -a + w),
      rect(-w / 2, g / 2 + tw, w / 2, a - w), rect(-w / 2, -a + w, w / 2, -g / 2 - tw), rect(-pl, g / 2, pl, g / 2 + tw), rect(-pl, -g / 2 - tw, pl, -g / 2)];
  } else if (type === 'cross') out = [crossPoly(p.length, min(p.width, p.length))];
  else if (type === 'patch') out = [rect(-p.width / 2, -p.height / 2, p.width / 2, p.height / 2)];
  else if (type === 'fishnet') { const X = p.px / 2, Y = p.py / 2, hx = min(p.holeX, p.px * 0.98) / 2, hy = min(p.holeY, p.py * 0.98) / 2; out = [rect(-X, -Y, -hx, Y), rect(hx, -Y, X, Y), rect(-hx, hy, hx, Y), rect(-hx, -Y, hx, -hy)]; }
  else if (type === 'nanodisc') out = [circle(p.diameter / 2, max(8, round(p.segments)))];
  else if (type === 'nanopillar') out = [p.shape === 'square' ? rect(-p.diameter / 2, -p.diameter / 2, p.diameter / 2, p.diameter / 2) : circle(p.diameter / 2, max(8, round(p.segments)))];
  else if (type === 'nanobar') out = [rot(rect(-p.length / 2, -p.width / 2, p.length / 2, p.width / 2), num(params.rotation_deg, 0) * PI / 180)];
  else if (type === 'vAntenna') { const h = min(179, max(10, num(params.angle_deg, 90))) * PI / 360, L = p.armLength; out = [rot(pathToPolygon([[L * sin(h), L * cos(h)], [0, 0], [-L * sin(h), L * cos(h)]], p.width), num(params.rotation_deg, 0) * PI / 180)]; }
  else if (type === 'cWire') { const r = rect(-p.length / 2, -p.width / 2, p.length / 2, p.width / 2), s = (p.width + p.gap) / 2; out = p.pair ? [r.map(q => [q[0], q[1] + s]), r.map(q => [q[0], q[1] - s])] : [r]; }
  else if (type === 'jerusalemCross') { const a = p.length / 2, cw = min(p.capWidth, a * 0.5), w = min(p.width, p.capLength), c = p.capLength / 2; out = [crossPoly(p.length - 2 * cw, w), rect(a - cw, -c, a, c), rect(-a, -c, -a + cw, c), rect(-c, a - cw, c, a), rect(-c, -a, c, -a + cw)]; }
  else { const a = p.cavityWidth / 2, b = p.cavityHeight / 2, A = a + p.wall, B = b + p.wall, g = min(p.neckWidth / 2, a * 0.95), t = min(p.wall, A - g), n = num(params.neckLength, D.neckLength) > 0 ? p.neckLength : 0; out = [[[g, B + n], [g + t, B + n], [g + t, B], [A, B], [A, -B], [-A, -B], [-A, B], [-g - t, B], [-g - t, B + n], [-g, B + n], [-g, b], [-a, b], [-a, -b], [a, -b], [a, b], [g, b]]]; }
  return out.map(pts => ({ layer, pts: clean(pts) })).filter(q => q.pts.length >= 3);
}

// Replicate a unit cell. rotateFn(i, j) → rotation in radians about the cell centre. Hexagonal lattice: odd rows
// shifted by px/2 with row pitch py (default px·√3/2).
export function array(polys, { nx = 5, ny = 5, px, py, rotateFn, lattice = 'square' } = {}) {
  const hex = lattice === 'hex';
  let Nx = max(1, round(num(nx, 5))), Ny = max(1, round(num(ny, 5))), capped = null;
  let W = 0, H = 0; for (const q of polys) for (const p of q.pts) { W = max(W, 2 * abs(p[0])); H = max(H, 2 * abs(p[1])); }
  const Px = max(1e-9, num(px, W * 1.5 || 1)), Py = max(1e-9, num(py, hex ? Px * sqrt(3) / 2 : Px));
  while (Nx * Ny * polys.length > 40000) { if (Nx >= Ny) Nx--; else Ny--; capped = `array reduced to ${Nx} × ${Ny} cells to keep the layout below 40 000 polygons`; }
  const polygons = [];
  for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) {
    const cx = (i + 0.5 + (hex && j % 2 ? 0.5 : 0)) * Px, cy = (j + 0.5) * Py, a = rotateFn ? Number(rotateFn(i, j)) || 0 : 0;
    for (const q of polys) polygons.push({ layer: q.layer ?? 1, pts: (a ? rot(q.pts, a) : q.pts).map(p => [p[0] + cx, p[1] + cy]) });
  }
  const out = { kind: 'layout', polygons, units: 'nm', period: [Px, Py], cells: [Nx, Ny], lattice: hex ? 'hex' : 'square', extent: [(Nx + (hex && Ny > 1 ? 0.5 : 0)) * Px, Ny * Py] };
  if (!rotateFn && !hex) out.hierarchy = { unit: polys.map(q => ({ layer: q.layer ?? 1, pts: q.pts.map(p => [p[0], p[1]]) })), nx: Nx, ny: Ny, px: Px, py: Py, origin: [Px / 2, Py / 2] };
  if (capped) out.capped = capped;
  return out;
}

// Phase map → meta-atom layout.
//  'nanopillar': propagation phase of a truncated-waveguide pillar, φ(d) = (2π/λ)·n_eff(d)·H with the zeroth-order
//    effective-medium index n_eff = √(f·n_p² + (1 − f)·n_m²), f = πd²/(4p²). This is a first-pass estimate: a real
//    library must be calibrated by RCWA/FDTD because n_eff of a guided HE11 mode differs from the volume average.
//  'nanobar-PB': geometric (Pancharatnam–Berry) phase, exact for an ideal half-wave-plate element under circularly
//    polarised light: rotation θ = φ/2.
export function metasurfacePhaseLayout({ phases, period_nm = 400, library = 'nanopillar', params = {}, layer = 1 } = {}) {
  if (!Array.isArray(phases) || !phases.length || !Array.isArray(phases[0])) throw new Error('metasurfacePhaseLayout needs phases as a 2-D array [row][column] in radians');
  const p = max(1, num(period_nm, 400)), ny = phases.length, nx = phases[0].length, polygons = [], notes = [], TWO = 2 * PI, wrap = v => ((v % TWO) + TWO) % TWO;
  let capped = null; const rows = min(ny, 200), cols = min(nx, 200);
  if (rows < ny || cols < nx) capped = `phase map cropped to ${cols} × ${rows} cells`;
  const design = [];
  if (library === 'nanopillar') {
    const lam = max(1, num(params.wavelength_nm, 633)), np = max(1.01, num(params.n_pillar, 2.4)), nm = max(1, num(params.n_medium, 1)), Hh = max(1, num(params.height_nm, 600));
    const dmin = max(0.05, min(0.9, num(params.dmin_frac, 0.2))) * p, dmax = max(dmin + 1e-6, min(0.95, num(params.dmax_frac, 0.8)) * p);
    const neff = d => sqrt((PI * d * d / (4 * p * p)) * np * np + (1 - PI * d * d / (4 * p * p)) * nm * nm), ph = d => TWO / lam * neff(d) * Hh;
    const cover = ph(dmax) - ph(dmin);
    if (cover < TWO) notes.push(`phase coverage is ${(cover / PI).toFixed(2)}π < 2π for height ${Hh} nm: phases beyond the range are clamped (increase the pillar height or index)`);
    if (p > lam / nm) notes.push('period exceeds the wavelength in the surrounding medium: diffraction orders will appear');
    for (let j = 0; j < rows; j++) { const row = []; for (let i = 0; i < cols; i++) {
      const target = min(wrap(num(phases[j][i], 0)), cover); let lo = dmin, hi = dmax;
      for (let it = 0; it < 40; it++) { const m = (lo + hi) / 2; if (ph(m) - ph(dmin) < target) lo = m; else hi = m; }
      const d = (lo + hi) / 2; row.push(d);
      polygons.push({ layer, pts: circle(d / 2, 24, (i + .5) * p, (j + .5) * p) });
    } design.push(row); }
    return { kind: 'layout', polygons, units: 'nm', period: [p, p], cells: [cols, rows], extent: [cols * p, rows * p], library, diameters_nm: design, coverage_rad: cover, height_nm: Hh, model: 'φ = (2π/λ)·n_eff·H, n_eff = √(f n_p² + (1 − f) n_m²) (zeroth-order effective medium; calibrate against full-wave simulation)', notes, ...(capped ? { capped } : {}) };
  }
  if (library !== 'nanobar-PB') throw new Error(`Unknown meta-atom library "${library}". Use nanopillar or nanobar-PB.`);
  const L = min(num(params.length_nm, 0.65 * p), p * 0.95), w = min(num(params.width_nm, 0.25 * p), L), base = rect(-L / 2, -w / 2, L / 2, w / 2);
  if (hypot(L, w) > p) notes.push('bar diagonal exceeds the period: rotated bars can touch their neighbours');
  for (let j = 0; j < rows; j++) { const row = []; for (let i = 0; i < cols; i++) { const th = wrap(num(phases[j][i], 0)) / 2; row.push(th); polygons.push({ layer, pts: rot(base, th).map(q => [q[0] + (i + .5) * p, q[1] + (j + .5) * p]) }); } design.push(row); }
  notes.push('assumes each bar acts as a half-wave plate at the design wavelength; conversion efficiency must be checked by full-wave simulation');
  return { kind: 'layout', polygons, units: 'nm', period: [p, p], cells: [cols, rows], extent: [cols * p, rows * p], library, rotations_rad: design, model: 'Pancharatnam–Berry phase φ = 2σθ (σ = ±1 for circular handedness)', notes, ...(capped ? { capped } : {}) };
}

function segDist(a, b, c, d) { // distance between segments ab and cd (2-D)
  const pd = (p, q, r) => { const vx = r[0] - q[0], vy = r[1] - q[1], l = vx * vx + vy * vy, t = l ? max(0, min(1, ((p[0] - q[0]) * vx + (p[1] - q[1]) * vy) / l)) : 0; return hypot(p[0] - q[0] - t * vx, p[1] - q[1] - t * vy); };
  const o = (p, q, r) => (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
  if (o(a, b, c) * o(a, b, d) < 0 && o(c, d, a) * o(c, d, b) < 0) return 0;
  return min(pd(a, c, d), pd(b, c, d), pd(c, a, b), pd(d, a, b));
}
// Polygon counts and area per layer, bounding box, fill factor, and estimates of the smallest feature and gap.
// minFeature: smallest of the corner-to-corner edge lengths and the distances between facing (anti-parallel) edges of
// one polygon (line width, or an internal slot such as an SRR gap); facets of discretised curves are ignored. minGap: smallest clearance between different polygons on the same layer
// (touching polygons excluded). Overlapping polygons are counted twice in the fill factor.
export function layoutStats(layout) {
  const P = layout.polygons || [], perLayer = {}, lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
  let total = 0, minEdge = Infinity, minWidth = Infinity, verts = 0;
  const boxes = P.map(q => {
    const L = perLayer[q.layer] || (perLayer[q.layer] = { polygons: 0, area: 0 }), a = abs(polyArea(q.pts)); L.polygons++; L.area += a; total += a; verts += q.pts.length;
    const b = [Infinity, Infinity, -Infinity, -Infinity];
    for (const p of q.pts) { b[0] = min(b[0], p[0]); b[1] = min(b[1], p[1]); b[2] = max(b[2], p[0]); b[3] = max(b[3], p[1]); }
    lo[0] = min(lo[0], b[0]); lo[1] = min(lo[1], b[1]); hi[0] = max(hi[0], b[2]); hi[1] = max(hi[1], b[3]);
    return b;
  });
  const sample = P.length > 400 ? Array.from({ length: 400 }, (_, i) => Math.floor(i * P.length / 400)) : P.map((_, i) => i);
  for (const i of sample) {
    const pts = P[i].pts, n = pts.length; if (n > 400) continue;
    const dir = pts.map((p, a) => { const q = pts[(a + 1) % n], l = hypot(q[0] - p[0], q[1] - p[1]) || 1; return [(q[0] - p[0]) / l, (q[1] - p[1]) / l, l]; });
    const corner = a => { const u = dir[(a + n - 1) % n], v = dir[a]; return u[0] * v[0] + u[1] * v[1] < 0.866; }; // turn of 30° or more
    for (let a = 0; a < n; a++) {
      if (dir[a][2] > 1e-9 && dir[a][2] < minEdge && corner(a) && corner((a + 1) % n)) minEdge = dir[a][2]; // facets of a discretised curve are not features
      for (let b = a + 2; b < n; b++) {
        if ((a === 0 && b === n - 1) || dir[a][0] * dir[b][0] + dir[a][1] * dir[b][1] > -0.7) continue; // only facing (anti-parallel) edges measure a width or a slot
        const d = segDist(pts[a], pts[(a + 1) % n], pts[b], pts[(b + 1) % n]); if (d > 1e-9 && d < minWidth) minWidth = d;
      }
    }
  }
  let minGap = Infinity, touching = 0;
  const ext = max(hi[0] - lo[0], hi[1] - lo[1]) || 1, cs = ext / max(1, round(sqrt(P.length))), grid = new Map();
  boxes.forEach((b, i) => { const k = Math.floor((b[0] + b[2]) / 2 / cs) + ',' + Math.floor((b[1] + b[3]) / 2 / cs); (grid.get(k) || grid.set(k, []).get(k)).push(i); });
  for (const i of sample) {
    const b = boxes[i], gx = Math.floor((b[0] + b[2]) / 2 / cs), gy = Math.floor((b[1] + b[3]) / 2 / cs), A = P[i].pts; if (A.length > 200) continue;
    for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (const j of grid.get((gx + dx) + ',' + (gy + dy)) || []) {
      if (j === i || P[j].layer !== P[i].layer) continue; const c = boxes[j], B = P[j].pts; if (B.length > 200) continue;
      const bd = max(0, max(b[0] - c[2], c[0] - b[2]), max(b[1] - c[3], c[1] - b[3])); if (bd >= minGap) continue;
      let d = Infinity; for (let p = 0; p < A.length && d > 0; p++) for (let q = 0; q < B.length; q++) { const v = segDist(A[p], A[(p + 1) % A.length], B[q], B[(q + 1) % B.length]); if (v < d) d = v; }
      if (d < 1e-9) touching++; else if (d < minGap) minGap = d;
    }
  }
  const has = P.length > 0, bboxArea = has ? (hi[0] - lo[0]) * (hi[1] - lo[1]) : 0, ref = layout.extent ? layout.extent[0] * layout.extent[1] : bboxArea;
  return {
    nPolygons: P.length, nVertices: verts, layers: perLayer, bbox: has ? [lo, hi] : [[0, 0], [0, 0]], size: has ? [hi[0] - lo[0], hi[1] - lo[1]] : [0, 0], area: total,
    fillFactor: ref > 0 ? min(1, total / ref) : 0, fillReference: layout.extent ? 'array extent (cells × period)' : 'bounding box',
    minEdge: Number.isFinite(minEdge) ? minEdge : null, minFeature: Number.isFinite(min(minEdge, minWidth)) ? min(minEdge, minWidth) : null, minGap: Number.isFinite(minGap) ? minGap : null,
    touchingPairs: touching / 2, sampled: sample.length < P.length ? `feature and gap estimates from ${sample.length} of ${P.length} polygons` : null, units: layout.units || 'nm',
  };
}
