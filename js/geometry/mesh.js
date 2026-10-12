// Mesh / CAD-surface engine: triangle-mesh diagnostics and repair, parametric solids, reactor builder,
// nanorod arrays and image-based microstructure reconstruction. Pure ES module.
import { rng } from '../core/num.js';

const { sin, cos, sqrt, abs, min, max, PI, floor, round, atan2 } = Math;
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len3 = a => sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
const num = (v, d) => Number.isFinite(Number(v)) && v !== null && v !== '' ? Number(v) : d;

export function bboxOf(points) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const p of points) for (let a = 0; a < 3; a++) { const v = p[a] || 0; if (v < lo[a]) lo[a] = v; if (v > hi[a]) hi[a] = v; }
  return points.length ? [lo, hi] : [[0, 0, 0], [0, 0, 0]];
}

// ---------------------------------------------------------------- diagnostics
export function meshStats(mesh) {
  const V = mesh.vertices, F = mesh.faces, nV = V.length, nF = F.length, bbox = bboxOf(V), diag = len3(sub(bbox[1], bbox[0])) || 1;
  let area = 0, vol = 0, degenerate = 0;
  const und = new Map(), dir = new Map(), parent = new Int32Array(nV).map((_, i) => i), usedV = new Uint8Array(nV);
  const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
  for (const f of F) {
    const [a, b, c] = f, A = V[a], B = V[b], C = V[c];
    if (!A || !B || !C) { degenerate++; continue; }
    const n = cross(sub(B, A), sub(C, A)), ar = len3(n) / 2;
    if (a === b || b === c || a === c || ar < 1e-14 * diag * diag) degenerate++;
    area += ar; vol += (A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0])) / 6;
    usedV[a] = usedV[b] = usedV[c] = 1;
    for (const [p, q] of [[a, b], [b, c], [c, a]]) {
      if (p === q) continue;
      const k = p < q ? p * nV + q : q * nV + p; und.set(k, (und.get(k) || 0) + 1);
      const d = p * nV + q; dir.set(d, (dir.get(d) || 0) + 1);
      const rp = find(p), rq = find(q); if (rp !== rq) parent[rp] = rq;
    }
  }
  let boundary = 0, nonManifold = 0, inconsistent = 0;
  for (const c of und.values()) { if (c === 1) boundary++; else if (c > 2) nonManifold++; }
  for (const c of dir.values()) if (c > 1) inconsistent++;
  const seen = new Map(); let dup = 0, used = 0; const q = 1e-9 * diag, comps = new Set();
  for (let i = 0; i < nV; i++) {
    if (usedV[i]) { used++; comps.add(find(i)); }
    const k = round(V[i][0] / q) + ',' + round(V[i][1] / q) + ',' + round((V[i][2] || 0) / q);
    if (seen.has(k)) dup++; else seen.set(k, i);
  }
  const watertight = nF > 0 && boundary === 0 && nonManifold === 0, euler = used - und.size + nF, components = comps.size;
  return {
    nVertices: nV, nFaces: nF, bbox, size: sub(bbox[1], bbox[0]), surfaceArea: area, volume: vol, watertight, nonManifoldEdges: nonManifold, boundaryEdges: boundary,
    duplicateVertices: dup, degenerateFaces: degenerate, unreferencedVertices: nV - used, consistentNormals: inconsistent === 0 && nonManifold === 0, inconsistentEdges: inconsistent,
    outwardNormals: watertight ? vol > 0 : null, eulerCharacteristic: euler, components, genus: watertight && inconsistent === 0 ? components - euler / 2 : null, nEdges: und.size, units: mesh.units || '',
  };
}

// Conservative repair: weld coincident vertices, drop degenerate and duplicate faces, make winding consistent by
// propagation across manifold edges, orient closed components outwards. Holes are reported, never filled.
export function repairMesh(mesh, { tolerance } = {}) {
  const V0 = mesh.vertices, bbox = bboxOf(V0), diag = len3(sub(bbox[1], bbox[0])) || 1, tol = tolerance ?? 1e-7 * diag, report = [];
  const map = new Map(), remap = new Int32Array(V0.length), V = [];
  for (let i = 0; i < V0.length; i++) {
    const p = V0[i], k = round(p[0] / tol) + ',' + round(p[1] / tol) + ',' + round((p[2] || 0) / tol); let j = map.get(k);
    if (j === undefined) { j = V.length; V.push([p[0], p[1], p[2] || 0]); map.set(k, j); }
    remap[i] = j;
  }
  if (V.length < V0.length) report.push(`welded ${V0.length - V.length} coincident vertices (tolerance ${tol.toExponential(2)} ${mesh.units || 'units'})`);
  let F = [], degen = 0, dupF = 0; const fseen = new Set(), nV = V.length;
  for (const f of mesh.faces) {
    const a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
    if (a === undefined || b === undefined || c === undefined || a === b || b === c || a === c || len3(cross(sub(V[b], V[a]), sub(V[c], V[a]))) < 1e-14 * diag * diag) { degen++; continue; }
    const s = [a, b, c].sort((x, y) => x - y), k = (s[0] * nV + s[1]) * nV + s[2];
    if (fseen.has(k)) { dupF++; continue; }
    fseen.add(k); F.push([a, b, c]);
  }
  if (degen) report.push(`removed ${degen} degenerate (zero-area or repeated-index) faces`);
  if (dupF) report.push(`removed ${dupF} duplicate faces`);
  // winding propagation
  const adj = new Map();
  F.forEach((f, i) => { for (let e = 0; e < 3; e++) { const p = f[e], q = f[(e + 1) % 3], k = p < q ? p * nV + q : q * nV + p; let l = adj.get(k); if (!l) adj.set(k, l = []); l.push(i); } });
  const comp = new Int32Array(F.length).fill(-1); let nc = 0, flips = 0;
  const hasDir = (f, p, q) => (f[0] === p && f[1] === q) || (f[1] === p && f[2] === q) || (f[2] === p && f[0] === q);
  for (let s = 0; s < F.length; s++) {
    if (comp[s] >= 0) continue;
    const stack = [s]; comp[s] = nc;
    while (stack.length) {
      const i = stack.pop(), f = F[i];
      for (let e = 0; e < 3; e++) {
        const p = f[e], q = f[(e + 1) % 3], l = adj.get(p < q ? p * nV + q : q * nV + p);
        if (l.length !== 2) continue;
        const j = l[0] === i ? l[1] : l[0];
        if (comp[j] >= 0) continue;
        if (hasDir(F[j], p, q)) { F[j] = [F[j][0], F[j][2], F[j][1]]; flips++; }
        comp[j] = nc; stack.push(j);
      }
    }
    nc++;
  }
  if (flips) report.push(`reversed ${flips} faces to make winding consistent`);
  const cv = new Float64Array(nc), open = new Uint8Array(nc);
  F.forEach((f, i) => { const A = V[f[0]], B = V[f[1]], C = V[f[2]]; cv[comp[i]] += A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0]); });
  for (const [, l] of adj) if (l.length !== 2) for (const i of l) open[comp[i]] = 1;
  let inv = 0;
  // a closed component nested inside another (a cavity) is legitimately inward-facing; only flip when the whole mesh volume is negative
  const total = cv.reduce((a, b) => a + b, 0);
  if (total < 0) { F = F.map((f, i) => open[comp[i]] ? f : (inv++, [f[0], f[2], f[1]])); }
  if (inv) report.push(`net enclosed volume was negative: reversed ${inv} faces of closed components so that normals point outwards`);
  const usedIdx = new Int32Array(nV).fill(-1), V2 = [];
  F = F.map(f => f.map(i => { if (usedIdx[i] < 0) { usedIdx[i] = V2.length; V2.push(V[i]); } return usedIdx[i]; }));
  if (V2.length < V.length) report.push(`dropped ${V.length - V2.length} unreferenced vertices`);
  const out = { ...mesh, kind: 'mesh', vertices: V2, faces: F, units: mesh.units }, st = meshStats(out);
  if (st.boundaryEdges) report.push(`${st.boundaryEdges} open boundary edges remain (holes are not filled automatically)`);
  if (st.nonManifoldEdges) report.push(`${st.nonManifoldEdges} non-manifold edges remain (not altered)`);
  if (!report.length) report.push('no defects found; mesh unchanged');
  return { mesh: out, report, stats: st };
}

// scale (number or [sx,sy,sz]) → rotate ([rx,ry,rz] degrees, applied X then Y then Z) → translate.
export function transformMesh(mesh, { scale = 1, translate = [0, 0, 0], rotate = [0, 0, 0] } = {}) {
  const s = typeof scale === 'number' ? [scale, scale, scale] : scale, [a, b, c] = rotate.map(d => d * PI / 180);
  const ca = cos(a), sa = sin(a), cb = cos(b), sb = sin(b), cc = cos(c), sc = sin(c);
  const R = [[cc * cb, cc * sb * sa - sc * ca, cc * sb * ca + sc * sa], [sc * cb, sc * sb * sa + cc * ca, sc * sb * ca - cc * sa], [-sb, cb * sa, cb * ca]];
  const vertices = mesh.vertices.map(p => { const q = [p[0] * s[0], p[1] * s[1], (p[2] || 0) * s[2]]; return [R[0][0] * q[0] + R[0][1] * q[1] + R[0][2] * q[2] + translate[0], R[1][0] * q[0] + R[1][1] * q[1] + R[1][2] * q[2] + translate[1], R[2][0] * q[0] + R[2][1] * q[1] + R[2][2] * q[2] + translate[2]]; });
  const flip = s[0] * s[1] * s[2] < 0;
  return { ...mesh, vertices, faces: flip ? mesh.faces.map(f => [f[0], f[2], f[1]]) : mesh.faces.map(f => [...f]) };
}

export function mergeMeshes(list, units) {
  const vertices = [], faces = [], ranges = [];
  for (const m of list) { const b = vertices.length, f0 = faces.length; for (const p of m.vertices) vertices.push(p); for (const f of m.faces) faces.push([f[0] + b, f[1] + b, f[2] + b]); ranges.push([f0, m.faces.length]); }
  return { kind: 'mesh', vertices, faces, units: units || list[0]?.units || 'mm', ranges };
}
const signedVolume = m => { let v = 0; for (const f of m.faces) { const A = m.vertices[f[0]], B = m.vertices[f[1]], C = m.vertices[f[2]]; v += A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0]); } return v / 6; };
const areaOf = m => { let a = 0; for (const f of m.faces) a += len3(cross(sub(m.vertices[f[1]], m.vertices[f[0]]), sub(m.vertices[f[2]], m.vertices[f[0]]))) / 2; return a; };
const flipped = m => ({ ...m, faces: m.faces.map(f => [f[0], f[2], f[1]]) });
const outward = m => signedVolume(m) < 0 ? flipped(m) : m;
const moved = (m, o) => o && (o[0] || o[1] || o[2]) ? { ...m, vertices: m.vertices.map(p => [p[0] + o[0], p[1] + o[1], p[2] + o[2]]) } : m;

// ---------------------------------------------------------------- polygon helpers
export const polygonArea = pts => { let a = 0; for (let i = 0, n = pts.length; i < n; i++) { const p = pts[i], q = pts[(i + 1) % n]; a += p[0] * q[1] - q[0] * p[1]; } return a / 2; };
// Ear clipping for a simple polygon (no holes). Returns counter-clockwise index triples into pts.
export function triangulate(pts) {
  const n = pts.length; if (n < 3) return [];
  const idx = pts.map((_, i) => i); if (polygonArea(pts) < 0) idx.reverse();
  const tris = [], cr = (a, b, c) => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  let guard = 0;
  while (idx.length > 3 && guard++ < 4 * n * n) {
    let cut = false;
    for (let i = 0; i < idx.length; i++) {
      const ia = idx[(i + idx.length - 1) % idx.length], ib = idx[i], ic = idx[(i + 1) % idx.length], A = pts[ia], B = pts[ib], C = pts[ic];
      if (cr(A, B, C) <= 1e-14) { if (abs(cr(A, B, C)) <= 1e-14 && idx.length > 3) { idx.splice(i, 1); cut = true; break; } continue; }
      let ear = true;
      for (const k of idx) { if (k === ia || k === ib || k === ic) continue; const P = pts[k]; if (cr(A, B, P) >= 0 && cr(B, C, P) >= 0 && cr(C, A, P) >= 0) { ear = false; break; } }
      if (ear) { tris.push([ia, ib, ic]); idx.splice(i, 1); cut = true; break; }
    }
    if (!cut) break; // self-intersecting or degenerate input: stop with what was found
  }
  if (idx.length === 3) tris.push([idx[0], idx[1], idx[2]]);
  return tris;
}
// Prism from a 2-D triangulation (CCW triangles): top, bottom and a wall on every boundary edge. Watertight by construction.
export function extrude(pts, tris, z0, z1, units = 'mm') {
  const n = pts.length, vertices = [...pts.map(p => [p[0], p[1], z0]), ...pts.map(p => [p[0], p[1], z1])], faces = [], d = new Set();
  for (const [a, b, c] of tris) { faces.push([a + n, b + n, c + n], [a, c, b]); d.add(a * n + b); d.add(b * n + c); d.add(c * n + a); }
  for (const [a, b, c] of tris) for (const [p, q] of [[a, b], [b, c], [c, a]]) if (!d.has(q * n + p)) faces.push([p, q, q + n], [p, q + n, p + n]);
  return { kind: 'mesh', vertices, faces, units };
}
// Solid of revolution about z from an (r, z) profile: an open polyline whose ends lie on the axis, or a closed loop.
export function revolve(profile, segments = 32, closedLoop = false, units = 'mm') {
  const seg = max(3, round(segments)), vertices = [], faces = [], ring = [];
  for (const [r, z] of profile) {
    if (r < 1e-12) { ring.push([vertices.length]); vertices.push([0, 0, z]); }
    else { const b = vertices.length; ring.push(Array.from({ length: seg }, (_, s) => b + s)); for (let s = 0; s < seg; s++) vertices.push([r * cos(2 * PI * s / seg), r * sin(2 * PI * s / seg), z]); }
  }
  const m = profile.length, last = closedLoop ? m : m - 1;
  for (let i = 0; i < last; i++) {
    const A = ring[i], B = ring[(i + 1) % m];
    if (A.length === 1 && B.length === 1) continue;
    for (let s = 0; s < seg; s++) {
      const t = (s + 1) % seg;
      if (A.length === 1) faces.push([A[0], B[s], B[t]]);
      else if (B.length === 1) faces.push([A[s], B[0], A[t]]);
      else faces.push([A[s], B[s], B[t]], [A[s], B[t], A[t]]);
    }
  }
  return outward({ kind: 'mesh', vertices, faces, units });
}

// ---------------------------------------------------------------- parametric solids (mm; outward normals; watertight)
export function box({ size = [10, 10, 10], origin = [0, 0, 0] } = {}) {
  const [a, b, c] = size.map(v => max(1e-9, num(v, 1)));
  return moved(extrude([[0, 0], [a, 0], [a, b], [0, b]], [[0, 1, 2], [0, 2, 3]], 0, c), origin);
}
export function cylinder({ radius = 5, height = 10, segments = 32, origin } = {}) { return moved(revolve([[0, 0], [radius, 0], [radius, height], [0, height]], segments), origin); }
export function cone({ r1 = 5, r2 = 0, height = 10, segments = 32, origin } = {}) { return moved(revolve([[0, 0], [r1, 0], ...(r2 > 1e-12 ? [[r2, height]] : []), [0, height]], segments), origin); }
export function sphere({ radius = 5, segments = 24, origin } = {}) {
  const n = max(4, round(segments / 2)), prof = Array.from({ length: n + 1 }, (_, i) => { const t = -PI / 2 + PI * i / n; return [i === 0 || i === n ? 0 : radius * cos(t), radius * sin(t)]; });
  return moved(revolve(prof, segments), origin);
}
export function tube({ outer = 5, inner = 3, height = 10, segments = 32, origin } = {}) { const ri = min(inner, outer * 0.999); return moved(revolve([[ri, 0], [outer, 0], [outer, height], [ri, height]], segments, true), origin); }
export function torus({ R = 10, r = 2, segments = 32, tubeSegments = 16, origin } = {}) {
  const rr = min(r, R * 0.999), prof = Array.from({ length: tubeSegments }, (_, i) => [R + rr * cos(2 * PI * i / tubeSegments), rr * sin(2 * PI * i / tubeSegments)]);
  return moved(revolve(prof, segments, true), origin);
}
export function hexPrism({ acrossFlats = 10, height = 10, origin } = {}) {
  const R = acrossFlats / sqrt(3), pts = Array.from({ length: 6 }, (_, i) => [R * cos(PI / 3 * i), R * sin(PI / 3 * i)]);
  return moved(extrude(pts, [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]], 0, height), origin);
}
// Rectangular plate with an nx × ny array of through holes (one hole per pitch cell).
export function plateWithHoles({ size = [40, 40], thickness = 2, holes = {}, origin } = {}) {
  const nx = max(1, round(num(holes.nx, 4))), ny = max(1, round(num(holes.ny, 4))), px = size[0] / nx, py = size[1] / ny;
  const d = min(num(holes.diameter, 0.5 * min(px, py)), 0.95 * min(px, py)), seg = max(8, 4 * round(num(holes.segments, 16) / 4)), q = seg / 4;
  const pts = [], tris = [], key = new Map();
  const add = (x, y) => { const k = round(x / px * 1e6) + ',' + round(y / py * 1e6); let i = key.get(k); if (i === undefined) { i = pts.length; pts.push([x, y]); key.set(k, i); } return i; };
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const cx = (i + .5) * px, cy = (j + .5) * py, cor = [[px / 2, -py / 2], [px / 2, py / 2], [-px / 2, py / 2], [-px / 2, -py / 2]], P = [], C = [];
    for (let s = 0; s < 4; s++) for (let t = 0; t < q; t++) {
      const a = cor[s], b = cor[(s + 1) % 4], x = a[0] + (b[0] - a[0]) * t / q, y = a[1] + (b[1] - a[1]) * t / q, l = sqrt(x * x + y * y);
      P.push(add(cx + x, cy + y)); C.push(add(cx + d / 2 * x / l, cy + d / 2 * y / l));
    }
    for (let s = 0; s < seg; s++) { const t = (s + 1) % seg; tris.push([P[s], P[t], C[t]], [P[s], C[t], C[s]]); }
  }
  return moved(extrude(pts, tris, 0, thickness), origin);
}


export const plate = plateWithHoles, plateWithHoleArray = plateWithHoles;
// ---------------------------------------------------------------- reactor builder
// Dimensions in mm. Each solid part is a closed shell; `regions[i].faces = [firstFace, faceCount]` indexes the merged mesh.
// Fluid, inlet and outlet regions are analytical (no faces). Volumes of solids come from the meshed shells (faceted
// cylinders: slightly below the ideal value); fluid volume and wetted area use the exact analytical dimensions.
export function reactor({ type = 'autoclave', dims = {}, features = {} } = {}) {
  const D = (k, d) => max(1e-6, num(dims[k], d)), F = (k, d = true) => features[k] === undefined ? d : features[k];
  const parts = [], regions = [], notes = [], seg = 40, mm3 = 1e-9, mm2 = 1e-6;
  const solid = (name, role, mesh) => parts.push({ name, role, mesh });
  const virt = (name, role, volume_mm3, area_mm2) => regions.push({ name, role, volume_m3: volume_mm3 * mm3, area_m2: area_mm2 * mm2, faces: null });
  const hollow = (ro, ri, zb, zt, floor_, lid) => mergeMeshes([revolve([[0, zb - floor_], [ro, zb - floor_], [ro, zt + lid], ...(lid > 0 ? [[0, zt + lid]] : [[ri, zt], [ri, zb], [0, zb]])], seg), ...(lid > 0 ? [flipped(revolve([[0, zb], [ri, zb], [ri, zt], [0, zt]], seg))] : [])]);
  const alongX = (m, x0 = 0) => transformMesh(m, { rotate: [0, 90, 0], translate: [x0, 0, 0] });
  let fluidV = 0, wetted = 0;
  if (type === 'autoclave') {
    const Di = D('innerDiameter', 30), H = D('innerHeight', 50), tw = D('wall', 6), fill = min(1, max(0.05, num(dims.fill, 0.7))), tl = F('liner') ? D('linerThickness', 3) : 0, Ri = Di / 2, hf = fill * H;
    solid('pressure vessel (steel shell with lid)', 'wall', hollow(Ri + tl + tw, Ri + tl, -tl, H + tl, tw, tw));
    if (tl) solid('liner (PTFE cup with cap)', 'wall', hollow(Ri + tl, Ri, 0, H, tl, tl));
    fluidV = PI * Ri * Ri * hf; wetted = PI * Ri * Ri + PI * Di * hf;
    if (F('substrate')) { const s = min(D('substrateSize', 10), Di * 0.6), ts = D('substrateThickness', 0.5), z = min(hf * 0.4, H - ts); solid('substrate (seeded wafer piece)', 'substrate', box({ size: [s, s, ts], origin: [-s / 2, -s / 2, z] })); fluidV -= s * s * ts; wetted += 2 * s * s + 4 * s * ts; }
    if (F('heater')) solid('heater band', 'wall', tube({ outer: Ri + tl + tw + 2, inner: Ri + tl + tw, height: H * 0.6, segments: seg, origin: [0, 0, H * 0.2] }));
    virt('reaction solution', 'fluid', fluidV, wetted);
    virt('head space (vapour)', 'fluid', PI * Ri * Ri * (H - hf), PI * Di * (H - hf) + PI * Ri * Ri);
    notes.push(`fill fraction ${(fill * 100).toFixed(0)} % of the liner volume; autogenous pressure depends on fill and temperature`);
  } else if (type === 'tubular' || type === 'cvdTube') {
    const cvd = type === 'cvdTube', Di = D(cvd ? 'tubeInnerDiameter' : 'innerDiameter', cvd ? 50 : 10), L = D(cvd ? 'tubeLength' : 'length', cvd ? 600 : 500), tw = D('wall', cvd ? 2.5 : 1.5), Ri = Di / 2;
    solid(cvd ? 'quartz process tube' : 'tube wall', 'wall', alongX(tube({ outer: Ri + tw, inner: Ri, height: L, segments: seg })));
    fluidV = PI * Ri * Ri * L; wetted = PI * Di * L;
    if (F('heater', cvd)) { const Lh = min(L, D('hotZoneLength', L * (cvd ? 0.35 : 0.6))); solid(cvd ? 'furnace hot zone (heater sleeve)' : 'heater sleeve', 'wall', alongX(tube({ outer: Ri + tw + (cvd ? 12 : 3), inner: Ri + tw + (cvd ? 4 : 0.5), height: Lh, segments: seg }), (L - Lh) / 2)); }
    if (F('substrate', cvd)) { const s = min(D('substrateSize', cvd ? 20 : Di * 0.5), Di * 0.65), ts = D('substrateThickness', 0.5); solid('substrate', 'substrate', box({ size: [s, s, ts], origin: [L / 2 - s / 2, -s / 2, -Ri * 0.5] })); fluidV -= s * s * ts; wetted += 2 * s * s + 4 * s * ts; }
    virt(cvd ? 'process gas' : 'process fluid', 'fluid', fluidV, wetted);
    if (F('inlet')) virt('inlet', 'inlet', 0, PI * Ri * Ri);
    if (F('outlet')) virt('outlet', 'outlet', 0, PI * Ri * Ri);
  } else if (type === 'cstr') {
    const Di = D('diameter', 100), H = D('height', 120), tw = D('wall', 3), fill = min(1, max(0.05, num(dims.fill, 0.8))), Ri = Di / 2, hf = fill * H, Dimp = min(D('impellerDiameter', Di * 0.4), Di * 0.9), nb = max(0, round(num(dims.baffles, 4)));
    solid('tank (open-top vessel)', 'wall', hollow(Ri + tw, Ri, 0, H, tw, 0));
    fluidV = PI * Ri * Ri * hf; wetted = PI * Ri * Ri + PI * Di * hf;
    const rs = max(1, Dimp * 0.06), zi = H / 3, bw = Dimp / 5, bt = max(0.5, Dimp / 40), bl = Dimp / 2 - rs;
    solid('impeller shaft', 'wall', cylinder({ radius: rs, height: H - zi + 10, segments: 16, origin: [0, 0, zi] }));
    fluidV -= PI * rs * rs * max(0, hf - zi); wetted += 2 * PI * rs * max(0, hf - zi);
    for (let i = 0; i < 4; i++) solid(`impeller blade ${i + 1}`, 'wall', transformMesh(box({ size: [bl, bt, bw], origin: [rs, -bt / 2, zi] }), { rotate: [0, 0, 90 * i] }));
    fluidV -= 4 * bl * bt * bw; wetted += 4 * 2 * (bl * bw + bl * bt);
    for (let i = 0; i < nb; i++) { const w = Di / 12, t = max(1, Di / 60); solid(`baffle ${i + 1}`, 'wall', transformMesh(box({ size: [w, t, H * 0.9], origin: [Ri - w, -t / 2, H * 0.05] }), { rotate: [0, 0, 360 / nb * i + 45] })); fluidV -= w * t * min(hf, H * 0.9); wetted += 2 * w * min(hf, H * 0.9); }
    const dp = D('portDiameter', Di / 10);
    if (F('inlet')) { solid('inlet nozzle', 'inlet', transformMesh(tube({ outer: dp / 2 + 1, inner: dp / 2, height: 15, segments: 16 }), { rotate: [0, 90, 0], translate: [Ri + tw, 0, H * 0.85] })); }
    if (F('outlet')) { solid('outlet nozzle', 'outlet', transformMesh(tube({ outer: dp / 2 + 1, inner: dp / 2, height: 15, segments: 16 }), { rotate: [0, -90, 0], translate: [-Ri - tw, 0, H * 0.1] })); }
    if (F('heater', false)) solid('heating jacket', 'wall', tube({ outer: Ri + tw + 6, inner: Ri + tw + 1, height: hf, segments: seg }));
    virt('reaction mixture', 'fluid', fluidV, wetted);
    notes.push(`liquid level ${hf.toFixed(1)} mm; impeller clearance H/3 and D_imp/D_tank = ${(Dimp / Di).toFixed(2)} (standard stirred-tank proportions are 0.3–0.5)`);
  } else if (type === 'porousSlab') {
    const L = D('length', 50), W = D('width', 50), T = D('thickness', 10), phi = min(0.99, max(0.01, num(dims.porosity, 0.4))), dp = D('particleDiameter', 0.5), tw = D('wall', 2);
    solid('porous slab (bulk envelope)', 'porous', box({ size: [L, W, T] }));
    if (F('liner', true)) { solid('housing plate (bottom side)', 'wall', box({ size: [L, tw, T], origin: [0, -tw, 0] })); solid('housing plate (top side)', 'wall', box({ size: [L, tw, T], origin: [0, W, 0] })); }
    fluidV = phi * L * W * T; wetted = 6 * (1 - phi) / dp * L * W * T;
    virt('pore fluid', 'fluid', fluidV, wetted);
    if (F('inlet')) virt('inlet face', 'inlet', 0, W * T);
    if (F('outlet')) virt('outlet face', 'outlet', 0, W * T);
    notes.push(`wetted area from the packed-sphere relation S_v = 6(1 − φ)/d_p with d_p = ${dp} mm; replace with a measured specific surface when available`);
  } else if (type === 'electrochemCell') {
    const L = D('length', 60), W = D('width', 40), H = D('height', 50), tw = D('wall', 3), fill = min(1, max(0.05, num(dims.fill, 0.8))), hf = fill * H;
    const ew = min(D('electrodeWidth', 20), W * 0.9), eh = min(D('electrodeHeight', 40), H), et = D('electrodeThickness', 1), gap = min(D('gap', 20), L - 4 * et), clr = min(5, H * 0.1);
    solid('cell base', 'wall', box({ size: [L + 2 * tw, W + 2 * tw, tw], origin: [-tw, -tw, -tw] }));
    solid('cell wall −x', 'wall', box({ size: [tw, W + 2 * tw, H], origin: [-tw, -tw, 0] })); solid('cell wall +x', 'wall', box({ size: [tw, W + 2 * tw, H], origin: [L, -tw, 0] }));
    solid('cell wall −y', 'wall', box({ size: [L, tw, H], origin: [0, -tw, 0] })); solid('cell wall +y', 'wall', box({ size: [L, tw, H], origin: [0, W, 0] }));
    fluidV = L * W * hf; wetted = L * W + 2 * (L + W) * hf;
    const ne = F('electrodes', 2) === 3 ? 3 : F('electrodes', 2) ? 2 : 0, imm = max(0, min(eh, hf - clr));
    if (ne) {
      solid(F('substrate') ? 'working electrode (substrate)' : 'working electrode', F('substrate') ? 'substrate' : 'electrode', box({ size: [et, ew, eh], origin: [L / 2 - gap / 2 - et, (W - ew) / 2, clr] }));
      solid('counter electrode', 'electrode', box({ size: [et, ew, eh], origin: [L / 2 + gap / 2, (W - ew) / 2, clr] }));
      fluidV -= 2 * et * ew * imm; wetted += 2 * (2 * ew * imm + 2 * et * imm + et * ew);
      if (ne === 3) { const rr = 1.5; solid('reference electrode', 'electrode', cylinder({ radius: rr, height: eh, segments: 12, origin: [L / 2, W * 0.2, clr] })); fluidV -= PI * rr * rr * imm; wetted += 2 * PI * rr * imm + PI * rr * rr; }
      notes.push(`immersed electrode height ${imm.toFixed(1)} mm; facing electrode area ${(ew * imm).toFixed(1)} mm² per electrode; inter-electrode gap ${gap} mm`);
    }
    virt('electrolyte', 'fluid', fluidV, wetted);
  } else throw new Error(`Unknown reactor type "${type}". Use autoclave, tubular, cstr, porousSlab, cvdTube or electrochemCell.`);
  const mesh = mergeMeshes(parts.map(p => p.mesh), 'mm'), solids = parts.map((p, i) => ({ name: p.name, role: p.role, volume_m3: abs(signedVolume(p.mesh)) * mm3, area_m2: areaOf(p.mesh) * mm2, faces: mesh.ranges[i] }));
  delete mesh.ranges;
  return { type, mesh, regions: [...solids, ...regions], volume_m3: max(0, fluidV) * mm3, wettedArea_m2: wetted * mm2, notes };
}

// ---------------------------------------------------------------- nanorod arrays
// Rods standing on a substrate. tilt = polar angle from the surface normal (deg), azimuth in deg. `hex` marks a
// hexagonal cross-section (wurtzite habit; diameter = across corners); placement is random sequential with a
// minimum centre spacing of one diameter where the density allows, or a jittered hexagonal grid (arrangement:'hexagonal').
export function nanorodArray({ density_per_um2 = 50, length_nm = 1000, diameter_nm = 80, tiltSpread_deg = 5, area_nm = [1000, 1000], seed = 1, hex = true, lengthCV = 0.1, diameterCV = 0.1, arrangement = 'random' } = {}) {
  const Lx = max(1, num(area_nm[0], 1000)), Ly = max(1, num(area_nm[1], 1000)), rand = rng(seed ?? 1), dens = max(0, num(density_per_um2, 50));
  let n = round(dens * Lx * Ly * 1e-6), capped = null;
  if (n > 4000) { capped = `requested ${n} rods; limited to 4000 for interactive use (reduce the area to see the full density)`; n = 4000; }
  const d0 = max(1, num(diameter_nm, 80)), l0 = max(1, num(length_nm, 1000)), ts = max(0, num(tiltSpread_deg, 5)), rods = [];
  const ln = (m, cv) => cv > 0 ? m * Math.exp(rand.normal(-0.5 * Math.log(1 + cv * cv), sqrt(Math.log(1 + cv * cv)))) : m;
  const cellN = max(1, floor(Lx / d0)), cellM = max(1, floor(Ly / d0)), grid = new Map(); let overlaps = 0;
  const free = (x, y) => { const i = floor(x / Lx * cellN), j = floor(y / Ly * cellM); for (let a = i - 1; a <= i + 1; a++) for (let b = j - 1; b <= j + 1; b++) for (const q of grid.get(a * 100003 + b) || []) if ((q[0] - x) ** 2 + (q[1] - y) ** 2 < d0 * d0) return false; return true; };
  const put = (x, y) => { const k = floor(x / Lx * cellN) * 100003 + floor(y / Ly * cellM); if (!grid.has(k)) grid.set(k, []); grid.get(k).push([x, y]); rods.push({ x, y, length: ln(l0, lengthCV), diameter: ln(d0, diameterCV), tilt: abs(rand.normal(0, ts)), azimuth: rand() * 360 }); };
  if (arrangement === 'hexagonal' || arrangement === 'square') {
    const hx = arrangement === 'hexagonal', a = sqrt(Lx * Ly / max(1, n) * (hx ? 2 / sqrt(3) : 1)), dy = hx ? a * sqrt(3) / 2 : a;
    for (let j = 0, y = dy / 2; y < Ly && rods.length < n + 64; j++, y += dy) for (let x = hx && j % 2 ? a : a / 2; x < Lx; x += a) put(min(Lx, max(0, x + rand.normal(0, 0.03 * a))), min(Ly, max(0, y + rand.normal(0, 0.03 * a))));
  } else for (let i = 0; i < n; i++) {
    let x = rand() * Lx, y = rand() * Ly, tries = 0;
    while (!free(x, y) && tries++ < 30) { x = rand() * Lx; y = rand() * Ly; }
    if (tries > 30) overlaps++;
    put(x, y);
  }
  const cover = rods.reduce((s, r) => s + (hex ? 3 * sqrt(3) / 8 : PI / 4) * r.diameter ** 2, 0) / (Lx * Ly), N = rods.length || 1;
  const out = {
    kind: 'rods', rods, area: [Lx, Ly], units: 'nm', section: hex ? 'hex' : 'round', arrangement, count: rods.length,
    density_per_um2: rods.length / (Lx * Ly * 1e-6), coverage: min(1, cover), meanLength_nm: rods.reduce((s, r) => s + r.length, 0) / N, meanDiameter_nm: rods.reduce((s, r) => s + r.diameter, 0) / N,
    meanSpacing_nm: sqrt(Lx * Ly / N), meanTilt_deg: rods.reduce((s, r) => s + r.tilt, 0) / N,
    surfaceAreaEnhancement: 1 + rods.reduce((s, r) => s + (hex ? 3 : PI) * r.diameter * r.length, 0) / (Lx * Ly),
  };
  if (overlaps) out.note = `${overlaps} rods overlap a neighbour: the requested density exceeds what one-diameter spacing allows (coalescence would occur in growth)`;
  if (capped) out.capped = capped;
  return out;
}
export function rodsToMesh(geom, { maxRods = 2500 } = {}) {
  const list = geom.rods.slice(0, maxRods), sides = geom.section === 'round' ? 12 : 6, V = [], Fc = [];
  for (const r of list) {
    const t = (r.tilt || 0) * PI / 180, az = (r.azimuth || 0) * PI / 180, w = [sin(t) * cos(az), sin(t) * sin(az), cos(t)], R = r.diameter / 2, b = V.length;
    const u = abs(w[2]) < 0.999 ? (() => { const q = cross([0, 0, 1], w), l = len3(q); return q.map(v => v / l); })() : [1, 0, 0], v = cross(w, u);
    for (let s = 0; s < sides; s++) { const a = 2 * PI * s / sides, q = [R * (cos(a) * u[0] + sin(a) * v[0]), R * (cos(a) * u[1] + sin(a) * v[1]), R * (cos(a) * u[2] + sin(a) * v[2])]; V.push([r.x + q[0], r.y + q[1], q[2]], [r.x + q[0] + r.length * w[0], r.y + q[1] + r.length * w[1], q[2] + r.length * w[2]]); }
    for (let s = 0; s < sides; s++) { const a0 = b + 2 * s, a1 = b + 2 * ((s + 1) % sides); Fc.push([a0, a1, a1 + 1], [a0, a1 + 1, a0 + 1]); }
    for (let s = 1; s < sides - 1; s++) Fc.push([b, b + 2 * (s + 1), b + 2 * s], [b + 1, b + 2 * s + 1, b + 2 * (s + 1) + 1]);
  }
  const out = { kind: 'mesh', vertices: V, faces: Fc, units: geom.units || 'nm' };
  if (list.length < geom.rods.length) out.truncated = `only the first ${maxRods} of ${geom.rods.length} rods were meshed`;
  return out;
}

// ---------------------------------------------------------------- image-based reconstruction
export function toGrey(data, width, height) {
  const n = width * height, g = new Uint8Array(n), ch = round(data.length / n);
  if (ch !== 1 && ch !== 3 && ch !== 4) throw new Error(`Image data length ${data.length} does not match ${width}×${height} with 1, 3 or 4 channels`);
  if (ch === 1) g.set(data.subarray ? data.subarray(0, n) : data.slice(0, n)); else for (let i = 0; i < n; i++) g[i] = round(0.299 * data[i * ch] + 0.587 * data[i * ch + 1] + 0.114 * data[i * ch + 2]);
  return g;
}
// Otsu's method: threshold maximising the between-class variance of the grey-level histogram.
export function otsu(grey) {
  const h = new Float64Array(256); for (let i = 0; i < grey.length; i++) h[grey[i]]++;
  let sum = 0; for (let i = 0; i < 256; i++) sum += i * h[i];
  let wB = 0, sB = 0, best = 0, thr = 127; const N = grey.length;
  for (let t = 0; t < 256; t++) { wB += h[t]; if (!wB) continue; const wF = N - wB; if (!wF) break; sB += t * h[t]; const d = sB / wB - (sum - sB) / wF, v = wB * wF * d * d; if (v > best) { best = v; thr = t; } }
  return thr;
}
// Segment a micrograph and extrude it `depth` voxels (a 2-D section carries no depth information: depth > 1 is a
// columnar assumption). Phase 1 = brighter than the threshold (or darker with invert). Image row 0 (top) maps to the
// highest y so that the geometry is displayed upright.
export function voxelsFromImage({ width, height, data }, { threshold = 'otsu', invert = false, depth = 1, pixelSize = 1, units = 'px' } = {}) {
  if (!(width > 0 && height > 0) || !data) throw new Error('voxelsFromImage needs {width, height, data}');
  const g = toGrey(data, width, height), thr = typeof threshold === 'number' && Number.isFinite(threshold) ? threshold : otsu(g), nz = max(1, min(512, round(depth))), out = new Uint8Array(width * height * nz);
  let c = 0;
  for (let y = 0; y < height; y++) for (let x = 0; x < width; x++) { const v = (g[x + width * (height - 1 - y)] > thr) !== !!invert ? 1 : 0; c += v; if (v) for (let z = 0; z < nz; z++) out[x + width * (y + height * z)] = 1; }
  return { kind: 'voxels', nx: width, ny: height, nz, data: out, voxel: pixelSize, units, threshold: thr, thresholdMethod: typeof threshold === 'number' ? 'user value' : 'Otsu (maximum between-class variance)', phaseFraction: c / (width * height), assumption: nz > 1 ? 'section extruded along z (columnar microstructure assumed)' : 'single 2-D section' };
}
// Connected-component labelling (8-connectivity) of pixels brighter than the threshold; shape descriptors from
// second central moments (equivalent ellipse).
export function segmentParticles(grey, width, height, { threshold = 'otsu', invert = false, minArea = 1 } = {}) {
  const thr = typeof threshold === 'number' ? threshold : otsu(grey), n = width * height, lab = new Int32Array(n), stack = new Int32Array(n), out = [];
  const on = i => (grey[i] > thr) !== !!invert;
  for (let s = 0; s < n; s++) {
    if (lab[s] || !on(s)) continue;
    let top = 0, area = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, x0 = width, x1 = 0, y0 = height, y1 = 0, edge = false;
    stack[top++] = s; lab[s] = out.length + 1;
    while (top) {
      const p = stack[--top], x = p % width, y = (p / width) | 0;
      area++; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
      if (x === 0 || y === 0 || x === width - 1 || y === height - 1) edge = true;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) { const X = x + dx, Y = y + dy; if (X < 0 || Y < 0 || X >= width || Y >= height) continue; const q = X + width * Y; if (!lab[q] && on(q)) { lab[q] = lab[s]; stack[top++] = q; } }
    }
    if (area < minArea) continue;
    const cx = sx / area, cy = sy / area, mxx = sxx / area - cx * cx + 1 / 12, myy = syy / area - cy * cy + 1 / 12, mxy = sxy / area - cx * cy;
    const tr = mxx + myy, df = sqrt(max(0, (mxx - myy) ** 2 + 4 * mxy * mxy)), l1 = (tr + df) / 2, l2 = max(1e-12, (tr - df) / 2);
    out.push({ area_px: area, centroid: [cx, cy], bbox: [x0, y0, x1, y1], equivDiameter_px: 2 * sqrt(area / PI), aspect: sqrt(l1 / l2), orientation_deg: 0.5 * atan2(2 * mxy, mxx - myy) * 180 / PI, majorAxis_px: 4 * sqrt(l1), minorAxis_px: 4 * sqrt(l2), touchesEdge: edge });
  }
  return out;
}
