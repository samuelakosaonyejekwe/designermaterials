// On-demand geometry kernels (WebAssembly, vendored under /vendor): OpenCASCADE (occt-import-js) for STEP / IGES / BREP,
// openNURBS (rhino3dm) for Rhino .3dm and web-ifc for IFC. Nothing here is fetched or instantiated until loadKernel() or
// tessellate() is called. The loaders run under a strict Content-Security-Policy (script-src 'self' 'wasm-unsafe-eval'):
// the kernels are imported as same-origin ES modules, the .wasm files are fetched from the application's own vendor/
// directory, and no code is ever built from strings. In Node (tests) the same files are read from disk.
const { min, max, abs, ceil, hypot, PI } = Math;
const IS_NODE = typeof globalThis.process === 'object' && !!globalThis.process?.versions?.node && typeof globalThis.document === 'undefined' && typeof globalThis.WorkerGlobalScope === 'undefined';
const now = () => (globalThis.performance ? globalThis.performance.now() : Date.now());

// Pack registry. `bytes` is what a first use downloads (loader + WebAssembly binary); licence texts sit next to the files.
export const KERNELS = {
  occt: { id: 'occt', name: 'occt-import-js', engine: 'Open CASCADE Technology', version: '0.0.23', licence: 'LGPL-2.1 (Open CASCADE Technology: LGPL-2.1 with the Open CASCADE exception)', formats: ['step', 'iges', 'brep'], bytes: 7699271,
    loader: new URL('../../vendor/occt-import-js/occt-import-js.js', import.meta.url), wasm: new URL('../../vendor/occt-import-js/occt-import-js.wasm', import.meta.url) },
  rhino: { id: 'rhino', name: 'rhino3dm', engine: 'openNURBS', version: '8.35.0', licence: 'MIT', formats: ['3dm'], bytes: 2795907,
    loader: new URL('../../vendor/rhino3dm/rhino3dm.module.js', import.meta.url), wasm: new URL('../../vendor/rhino3dm/rhino3dm.wasm', import.meta.url) },
  ifc: { id: 'ifc', name: 'web-ifc', engine: 'web-ifc', version: '0.0.78', licence: 'MPL-2.0', formats: ['ifc'], bytes: 7504318,
    loader: new URL('../../vendor/web-ifc/web-ifc-api.js', import.meta.url), wasm: new URL('../../vendor/web-ifc/web-ifc.wasm', import.meta.url) },
};
const ALIAS = { occt: 'occt', 'occt-import-js': 'occt', opencascade: 'occt', step: 'occt', stp: 'occt', iges: 'occt', igs: 'occt', brep: 'occt', brp: 'occt', rhino: 'rhino', rhino3dm: 'rhino', opennurbs: 'rhino', '3dm': 'rhino', ifc: 'ifc', 'web-ifc': 'ifc', webifc: 'ifc' };
// Kernel id for a kernel name or a format id ('step', 'iges', 'brep', '3dm', 'ifc'); null when no kernel covers it.
export const kernelFor = name => ALIAS[String(name || '').toLowerCase()] || null;
export const DEFAULT_TESSELLATION = { linearDeflection: 0.001, deflectionType: 'relative', angularDeflection: 0.5 };

const loaded = {}, pending = {}, failed = {};
async function readBytes(url) {
  if (IS_NODE) { const fs = await import('node:fs'), { fileURLToPath } = await import('node:url'); return new Uint8Array(fs.readFileSync(fileURLToPath(url))); }
  const r = await fetch(url.href, { credentials: 'same-origin' }); if (!r.ok) throw new Error(`${url.pathname.split('/').pop()}: HTTP ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
}
// True when the kernel is already instantiated or its files can be reached (on disk in Node; from the network or the
// offline cache in a browser). Never instantiates anything.
export async function kernelAvailable(name) {
  const id = kernelFor(name), k = KERNELS[id]; if (!k) return false; if (loaded[id]) return true;
  try {
    if (IS_NODE) { const fs = await import('node:fs'), { fileURLToPath } = await import('node:url'); return fs.existsSync(fileURLToPath(k.loader)) && fs.existsSync(fileURLToPath(k.wasm)); }
    for (const u of [k.loader, k.wasm]) {
      let hit = false; if (globalThis.caches) { try { hit = !!(await globalThis.caches.match(u.href, { ignoreSearch: true })); } catch (e) { hit = false; } }
      if (!hit) { const r = await fetch(u.href, { method: 'HEAD', credentials: 'same-origin' }); if (!r.ok) return false; }
    }
    return true;
  } catch (e) { return false; }
}
export const kernelStatus = name => { const id = kernelFor(name); return !id ? 'unknown' : loaded[id] ? 'loaded' : pending[id] ? 'loading' : failed[id] ? 'failed: ' + failed[id] : 'not loaded'; };

async function instantiate(id) {
  const k = KERNELS[id], log = [], sink = s => { if (log.length < 200) log.push(String(s)); };
  if (id === 'ifc') {
    // web-ifc's browser build fetches its own binary; it is pointed at the vendored file. In Node the same build is used:
    // the file is served to it from disk and the browser-only environment check is satisfied for the duration of start-up.
    const W = await import(k.loader.href), api = new W.IfcAPI();
    if (!IS_NODE) { await api.Init(() => k.wasm.href, true); }
    else {
      const bin = await readBytes(k.wasm), url = 'https://vendor.invalid/web-ifc.wasm', realFetch = globalThis.fetch, proc = globalThis.process, hadWindow = 'window' in globalThis; let init;
      globalThis.fetch = async (u, o) => String(u) === url ? new Response(bin, { headers: { 'content-type': 'application/wasm' } }) : realFetch(u, o);
      try { globalThis.process = undefined; if (!hadWindow) globalThis.window = globalThis; init = api.Init(() => url, true); } finally { globalThis.process = proc; if (!hadWindow) delete globalThis.window; }
      try { await init; } finally { globalThis.fetch = realFetch; }
    }
    if (W.LogLevel && api.SetLogLevel) { try { api.SetLogLevel(W.LogLevel.LOG_LEVEL_OFF); } catch (e) { /* optional */ } }
    return { id, api, module: W, log, version: k.version };
  }
  const [factory, wasmBinary] = await Promise.all([import(k.loader.href).then(m => m.default), readBytes(k.wasm)]);
  const arg = { wasmBinary, locateFile: p => p, print: sink, printErr: sink };
  if (id === 'occt' && IS_NODE) { arg.nodeFs = await import('node:fs'); arg.nodePath = await import('node:path'); }
  return { id, api: await factory(arg), log, version: k.version };
}
// loadKernel('occt' | 'rhino' | 'ifc' or a format id) → the instantiated kernel (cached). Rejects when the files cannot be
// fetched or instantiated; a later call tries again.
export function loadKernel(name) {
  const id = kernelFor(name); if (!id) return Promise.reject(new Error(`no geometry kernel is registered for "${name}"`));
  if (loaded[id]) return Promise.resolve(loaded[id]);
  if (!pending[id]) pending[id] = instantiate(id).then(k => { loaded[id] = k; delete pending[id]; delete failed[id]; return k; }, e => { delete pending[id]; failed[id] = (e && e.message ? e.message : String(e)).slice(0, 300); throw new Error(`${KERNELS[id].name} kernel could not be loaded (${failed[id]})${/CompileError|Content Security Policy|unsafe-eval/i.test(failed[id]) ? " — the page's Content-Security-Policy must allow 'wasm-unsafe-eval' in script-src for WebAssembly kernels" : ''}`); });
  return pending[id];
}

// ---------------------------------------------------------------- body assembly
// Weld coincident vertices of one body (tolerance-safe: neighbouring hash cells are checked), drop collapsed triangles,
// measure area / signed volume / closedness and orient closed bodies outwards.
function weldBody(P, I, tolRel = 1e-7) {
  const n = P.length / 3; let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) { const v = P[3 * i + a]; if (v < lo[a]) lo[a] = v; if (v > hi[a]) hi[a] = v; }
  const diag = n ? hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1 : 1, tol = tolRel * diag, q = 4 * tol, cells = new Map(), remap = new Int32Array(n), V = [];
  for (let i = 0; i < n; i++) {
    const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2], ax = Math.floor((x - tol) / q), bx = Math.floor((x + tol) / q), ay = Math.floor((y - tol) / q), by = Math.floor((y + tol) / q), az = Math.floor((z - tol) / q), bz = Math.floor((z + tol) / q);
    let found = -1;
    for (let cx = ax; cx <= bx && found < 0; cx++) for (let cy = ay; cy <= by && found < 0; cy++) for (let cz = az; cz <= bz && found < 0; cz++) { const l = cells.get(cx + ',' + cy + ',' + cz); if (l) for (const j of l) { const p = V[j]; if (abs(p[0] - x) <= tol && abs(p[1] - y) <= tol && abs(p[2] - z) <= tol) { found = j; break; } } }
    if (found < 0) { found = V.length; V.push([x, y, z]); const key = Math.floor(x / q) + ',' + Math.floor(y / q) + ',' + Math.floor(z / q), l = cells.get(key); if (l) l.push(found); else cells.set(key, [found]); }
    remap[i] = found;
  }
  const F = [], nV = V.length, und = new Map(), dir = new Set(); let area = 0, vol = 0, dropped = 0, clash = 0; const eps = 1e-14 * diag * diag;
  for (let t = 0; t + 2 < I.length; t += 3) {
    const a = remap[I[t]], b = remap[I[t + 1]], c = remap[I[t + 2]]; if (a === undefined || b === undefined || c === undefined || a === b || b === c || a === c) { dropped++; continue; }
    const A = V[a], B = V[b], C = V[c], ux = B[0] - A[0], uy = B[1] - A[1], uz = B[2] - A[2], vx = C[0] - A[0], vy = C[1] - A[1], vz = C[2] - A[2], ar = hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx) / 2;
    if (!(ar >= eps)) { dropped++; continue; }
    area += ar; vol += (A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0])) / 6; F.push([a, b, c]);
    for (const [p, r] of [[a, b], [b, c], [c, a]]) { const k = p < r ? p * nV + r : r * nV + p; und.set(k, (und.get(k) || 0) + 1); const d = p * nV + r; if (dir.has(d)) clash++; else dir.add(d); }
  }
  let open = 0, multi = 0; for (const c of und.values()) { if (c === 1) open++; else if (c > 2) multi++; }
  const closed = F.length > 0 && open === 0 && multi === 0, oriented = clash === 0; let flipped = false;
  if (closed && oriented && vol < 0) { for (const f of F) { const t = f[1]; f[1] = f[2]; f[2] = t; } vol = -vol; flipped = true; }
  // connected sets of triangles (shells)
  const par = new Int32Array(nV).map((_, i) => i), find = i => { while (par[i] !== i) { par[i] = par[par[i]]; i = par[i]; } return i; };
  for (const f of F) { const r = find(f[0]); for (const j of [f[1], f[2]]) { const s = find(j); if (s !== r) par[s] = r; } }
  const roots = new Set(); for (const f of F) roots.add(find(f[0]));
  return { vertices: V, faces: F, area, volume: vol, closed, oriented, flipped, openEdges: open, nonManifoldEdges: multi, shells: roots.size, dropped, bbox: [lo, hi] };
}
function assemble(raw, t0) {
  const vertices = [], faces = [], parts = [], warnings = [...(raw.warnings || [])]; let solids = 0, shells = 0, brepFaces = 0, emptyFaces = 0, flips = 0, dropped = 0, skipped = 0;
  for (const b of raw.bodies) {
    if (!b.positions || !b.indices || b.indices.length < 3) { skipped++; emptyFaces += b.emptyFaces || 0; brepFaces += b.brepFaces || 0; continue; }
    const w = weldBody(b.positions, b.indices), base = vertices.length, start = faces.length;
    if (!w.faces.length) { skipped++; continue; }
    for (const p of w.vertices) vertices.push(p); for (const f of w.faces) faces.push([f[0] + base, f[1] + base, f[2] + base]);
    const part = { name: b.name || `body ${parts.length + 1}`, colour: b.colour || null, faceRange: [start, faces.length], triangles: w.faces.length, area: w.area, closed: w.closed && w.oriented, shells: w.shells };
    if (part.closed) { part.volume = w.volume; solids++; } else part.openEdges = w.openEdges;
    if (b.brepFaces !== undefined) part.brepFaces = b.brepFaces; if (b.type) part.type = b.type; if (b.source) part.source = b.source;
    parts.push(part); shells += w.shells; brepFaces += b.brepFaces || 0; emptyFaces += b.emptyFaces || 0; if (w.flipped) flips++; dropped += w.dropped;
  }
  if (flips) warnings.push(`${flips} closed ${flips === 1 ? 'body was' : 'bodies were'} stored inside-out and ${flips === 1 ? 'was' : 'were'} re-oriented so that normals point outwards`);
  if (dropped) warnings.push(`${dropped} collapsed (zero-area) triangles produced by the tessellation were dropped`);
  if (emptyFaces) warnings.push(`${emptyFaces} B-rep ${emptyFaces === 1 ? 'face' : 'faces'} could not be tessellated by the kernel and ${emptyFaces === 1 ? 'is' : 'are'} missing from the mesh`);
  if (skipped) warnings.push(`${skipped} ${skipped === 1 ? 'body' : 'bodies'} produced no triangles (wireframe, construction or empty shapes)`);
  return { mesh: { kind: 'mesh', vertices, faces, units: raw.units || '' }, parts, units: raw.units || '', warnings,
    stats: { bodies: parts.length, solids, shells, faces: brepFaces, untessellatedFaces: emptyFaces, triangles: faces.length, vertices: vertices.length, seconds: (now() - t0) / 1000, kernel: raw.kernel, kernelVersion: raw.version, ...(raw.stats || {}) }, ...(raw.meta ? { meta: raw.meta } : {}) };
}

// ---------------------------------------------------------------- OpenCASCADE: STEP, IGES, BREP
function occtBodies(k, format, bytes, o) {
  const fn = { step: 'ReadStepFile', iges: 'ReadIgesFile', brep: 'ReadBrepFile' }[format], lin = o.linearDeflection ?? DEFAULT_TESSELLATION.linearDeflection, ang = o.angularDeflection ?? DEFAULT_TESSELLATION.angularDeflection, absolute = o.deflectionType === 'absolute';
  if (!(lin > 0) || !(ang > 0)) throw new Error('tessellate: linearDeflection and angularDeflection must be positive');
  k.log.length = 0;
  const r = k.api[fn](bytes, { linearUnit: 'millimeter', linearDeflectionType: absolute ? 'absolute_value' : 'bounding_box_ratio', linearDeflection: lin, angularDeflection: ang });
  if (!r || !r.success) throw new Error(`OpenCASCADE could not read the ${format.toUpperCase()} data${k.log.length ? ` (${k.log.slice(-2).join('; ').slice(0, 200)})` : ''}`);
  const names = new Map(), walk = (node, path) => { const p = node.name ? [...path, node.name] : path; for (const i of node.meshes || []) if (!names.has(i)) names.set(i, p); for (const c of node.children || []) walk(c, p); }; if (r.root) walk(r.root, []);
  const bodies = (r.meshes || []).map((m, i) => {
    const path = names.get(i) || [], bf = m.brep_faces || [], empty = bf.filter(f => !(f.last >= f.first)).length;
    return { name: m.name || path[path.length - 1] || '', source: path.join(' / '), colour: m.color ? [m.color[0], m.color[1], m.color[2]] : null, positions: m.attributes && m.attributes.position ? m.attributes.position.array : null, indices: m.index ? m.index.array : null, brepFaces: bf.length, emptyFaces: empty };
  });
  return { bodies, units: format === 'brep' ? '' : 'mm', kernel: 'occt-import-js', version: k.version, warnings: [], stats: { deflection: { type: absolute ? 'absolute' : 'relative', linear: lin, angular: ang } } };
}

// ---------------------------------------------------------------- openNURBS: Rhino 3DM
const RH_UNITS = { Angstroms: 'Å', Nanometers: 'nm', Microns: 'um', Millimeters: 'mm', Centimeters: 'cm', Meters: 'm', Inches: 'in', Feet: 'ft' };
const enumName = (E, v) => { for (const key of Object.keys(E)) if (key !== 'values' && key !== 'argCount' && E[key] === v) return key; return v && v.value !== undefined ? String(v.value) : ''; };
const free = x => { try { if (x && typeof x.delete === 'function') x.delete(); } catch (e) { /* already released */ } };
// An untrimmed surface sampled on a parameter grid; the grid is refined until the surface at cell centres lies within
// `dev` of the bilinear facets (or the cap is reached).
function sampleSurface(s, dev, cap = 192) {
  const du = s.domain(0), dv = s.domain(1); if (!du || !dv || !(du[1] > du[0]) || !(dv[1] > dv[0])) return null;
  const grid = n => { const P = new Float64Array((n + 1) * (n + 1) * 3); for (let j = 0; j <= n; j++) for (let i = 0; i <= n; i++) { const p = s.pointAt(du[0] + (du[1] - du[0]) * i / n, dv[0] + (dv[1] - dv[0]) * j / n), o = 3 * (j * (n + 1) + i); P[o] = p[0]; P[o + 1] = p[1]; P[o + 2] = p[2]; } return P; };
  let n = 8, P = grid(n), worst = Infinity;
  for (;;) {
    worst = 0; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const c = s.pointAt(du[0] + (du[1] - du[0]) * (i + 0.5) / n, dv[0] + (dv[1] - dv[0]) * (j + 0.5) / n), a = 3 * (j * (n + 1) + i), b = a + 3, d = 3 * ((j + 1) * (n + 1) + i), e = d + 3; let m = 0; for (let x = 0; x < 3; x++) { const mid = (P[a + x] + P[b + x] + P[d + x] + P[e + x]) / 4; m += (c[x] - mid) ** 2; } if (m > worst) worst = m; }
    worst = Math.sqrt(worst); if (worst <= dev || n * 2 > cap) break; n *= 2; P = grid(n);
  }
  const I = []; for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const a = j * (n + 1) + i, b = a + 1, d = a + n + 1, e = d + 1; I.push(a, b, e, a, e, d); }
  return { positions: P, indices: I, n, deviation: worst };
}
function rhinoBodies(k, bytes, o) {
  const rh = k.api, doc = rh.File3dm.fromByteArray(bytes); if (!doc) throw new Error('openNURBS could not read the 3DM archive');
  const warnings = [], bodies = [], inv = {}, T = rh.ObjectType, lin = o.linearDeflection ?? DEFAULT_TESSELLATION.linearDeflection;
  let units = ''; try { const u = enumName(rh.UnitSystem, doc.settings().modelUnitSystem); units = RH_UNITS[u] || ''; if (!units && u && u !== 'None' && u !== 'Unset') warnings.push(`model unit system "${u}" has no equivalent here: coordinates are left in file units`); } catch (e) { /* no settings table */ }
  const objs = doc.objects(), n = objs.count, byId = new Map(), refs = [], counts = { meshes: 0, brepCached: 0, brepPartial: 0, brepNone: 0, surfaces: 0, extrusionCached: 0, extrusionNone: 0, subd: 0, instances: 0, hidden: 0 };
  const meshArrays = m => { const j = m.toThreejsJSON(), d = j && j.data; return d && d.attributes && d.attributes.position && d.index ? { positions: d.attributes.position.array, indices: d.index.array } : null; };
  const merge = list => { let nv = 0, ni = 0; for (const a of list) { nv += a.positions.length; ni += a.indices.length; } const P = new Float64Array(nv), I = new Array(ni); let pv = 0, pi = 0; for (const a of list) { const base = pv / 3; for (let i = 0; i < a.positions.length; i++) P[pv++] = a.positions[i]; for (let i = 0; i < a.indices.length; i++) I[pi++] = a.indices[i] + base; } return { positions: P, indices: I }; };
  const sizeOf = g => { try { const b = g.getBoundingBox(), lo = b.min || b.Min, hi = b.max || b.Max, c = p => Array.isArray(p) ? p : [p.X ?? p.x, p.Y ?? p.y, p.Z ?? p.z], a = c(lo), z = c(hi); return hypot(z[0] - a[0], z[1] - a[1], z[2] - a[2]) || 1; } catch (e) { return 1; } };
  for (let i = 0; i < n; i++) {
    const obj = objs.get(i); if (!obj) continue; const g = obj.geometry(), at = obj.attributes(); if (!g) continue;
    const kind = enumName(T, g.objectType); inv[kind] = (inv[kind] || 0) + 1;
    let colour = null; try { const c = at.drawColor(doc); if (c) colour = [c.r / 255, c.g / 255, c.b / 255]; } catch (e) { colour = null; }
    let layer = ''; try { const L = doc.layers().get(at.layerIndex); if (L) layer = L.fullPath || L.name || ''; } catch (e) { layer = ''; }
    const body = { name: at.name || (layer ? `${kind} on ${layer}` : `${kind} ${i + 1}`), source: layer, colour, type: kind, id: at.id, definition: !!at.isInstanceDefinitionObject };
    let geo = null;
    try {
      if (kind === 'Mesh') { geo = meshArrays(g); if (geo) counts.meshes++; }
      else if (kind === 'Brep') {
        const fs = g.faces(), nf = fs.count, got = []; let missing = 0;
        for (let f = 0; f < nf; f++) { const face = fs.get(f), m = face.getMesh(rh.MeshType.Any); const a = m ? meshArrays(m) : null; if (a && a.indices.length) got.push(a); else missing++; free(m); free(face); }
        body.brepFaces = nf;
        if (got.length && !missing) { geo = merge(got); counts.brepCached++; body.type = 'Brep (render mesh stored in the file)'; }
        else if (!got.length && nf === 1 && g.isSurface) { const face = fs.get(0), s = face.underlyingSurface(), sm = s ? sampleSurface(s, lin * sizeOf(g)) : null; if (sm) { geo = sm; counts.surfaces++; body.type = `Brep, untrimmed surface sampled on a ${sm.n} × ${sm.n} parameter grid`; } free(s); free(face); }
        else if (got.length) { geo = merge(got); counts.brepPartial++; body.emptyFaces = missing; body.type = 'Brep (render mesh stored for some faces only)'; }
        else { counts.brepNone++; body.emptyFaces = nf; }
      } else if (kind === 'Surface') { const sm = sampleSurface(g, lin * sizeOf(g)); if (sm) { geo = sm; counts.surfaces++; body.type = `untrimmed surface sampled on a ${sm.n} × ${sm.n} parameter grid`; } }
      else if (kind === 'Extrusion') { const m = g.getMesh(rh.MeshType.Any), a = m ? meshArrays(m) : null; if (a && a.indices.length) { geo = a; counts.extrusionCached++; body.type = 'Extrusion (render mesh stored in the file)'; } else counts.extrusionNone++; free(m); }
      else if (kind === 'SubD') { const m = rh.Mesh.createFromSubDControlNet(g, false), a = m ? meshArrays(m) : null; if (a && a.indices.length) { geo = a; counts.subd++; body.type = 'SubD control net (not the limit surface)'; } free(m); }
      else if (kind === 'InstanceReference') { const X = g.xform; refs.push({ body, def: g.parentIdefId, xform: X && X.toFloatArray ? Array.from(X.toFloatArray(true)) : X }); free(X); counts.instances++; }
    } catch (e) { warnings.push(`object ${i + 1} (${kind}) could not be converted: ${e && e.message ? e.message : e}`); }
    if (geo) { Object.assign(body, geo); byId.set(String(body.id).toLowerCase(), body); if (!body.definition) { if (at.visible === false) counts.hidden++; bodies.push(body); } }
    free(g); free(at); free(obj);
  }
  // block instances: definition geometry placed by the reference transform (one level; nested references are reported)
  if (refs.length) {
    let placed = 0, nested = 0; const defs = doc.instanceDefinitions();
    for (const r of refs) {
      let d = null; try { const ix = defs.findId ? defs.findId(r.def) : null; d = ix && typeof ix === 'object' ? ix : null; if (!d) for (let i = 0; i < defs.count && !d; i++) { const c = defs.get(i); if (c && String(c.id).toLowerCase() === String(r.def).toLowerCase()) d = c; } } catch (e) { d = null; }
      const ids = d ? d.getObjectIds() : null, M = r.xform; if (!ids || !M || M.length !== 16) { nested++; continue; }
      for (const id of ids) { const src = byId.get(String(id).toLowerCase()); if (!src || !src.positions) { nested++; continue; } const P = new Float64Array(src.positions.length); for (let v = 0; v < P.length; v += 3) { const x = src.positions[v], y = src.positions[v + 1], z = src.positions[v + 2]; P[v] = M[0] * x + M[1] * y + M[2] * z + M[3]; P[v + 1] = M[4] * x + M[5] * y + M[6] * z + M[7]; P[v + 2] = M[8] * x + M[9] * y + M[10] * z + M[11]; }
        const det = M[0] * (M[5] * M[10] - M[6] * M[9]) - M[1] * (M[4] * M[10] - M[6] * M[8]) + M[2] * (M[4] * M[9] - M[5] * M[8]); let I = src.indices; if (det < 0) { I = Array.from(src.indices); for (let t = 0; t + 2 < I.length; t += 3) { const s = I[t + 1]; I[t + 1] = I[t + 2]; I[t + 2] = s; } }
        bodies.push({ ...src, name: `${src.name} (instance of ${d.name || 'block'})`, positions: P, indices: I, definition: false }); placed++; }
    }
    if (placed) warnings.push(`${placed} block-instance ${placed === 1 ? 'body was' : 'bodies were'} placed from ${refs.length} instance reference(s)`);
    if (nested) warnings.push(`${nested} instance reference(s) or members could not be expanded (nested blocks, linked files or definitions without a mesh)`);
  }
  if (counts.brepNone) warnings.push(`${counts.brepNone} Brep object(s) carry no stored render mesh: their trimmed NURBS faces were not tessellated (openNURBS reads and evaluates NURBS but contains no surface mesher)`);
  if (counts.brepPartial) warnings.push(`${counts.brepPartial} Brep object(s) had a stored mesh for only some of their faces`);
  if (counts.extrusionNone) warnings.push(`${counts.extrusionNone} extrusion object(s) carry no stored render mesh and were not tessellated`);
  if (counts.surfaces) warnings.push(`${counts.surfaces} untrimmed surface(s) were sampled on a regular parameter grid refined to a chord deviation of about ${lin} × the object size`);
  if (counts.subd) warnings.push(`${counts.subd} SubD object(s) are represented by their control net, which lies outside the smooth limit surface`);
  if (counts.hidden) warnings.push(`${counts.hidden} hidden object(s) are included`);
  const other = Object.entries(inv).filter(([t]) => !['Mesh', 'Brep', 'Surface', 'Extrusion', 'SubD', 'InstanceReference'].includes(t)); if (other.length) warnings.push(`not converted to mesh: ${other.map(([t, c]) => `${c} ${t}`).join(', ')}`);
  const meta = { objects: n, objectTypes: inv, layers: 0, applicationName: '', archiveVersion: null, converted: counts };
  try { meta.layers = doc.layers().count; meta.applicationName = doc.applicationName || ''; meta.archiveVersion = doc.archiveVersion; meta.createdBy = doc.createdBy || ''; meta.lastEditedBy = doc.lastEditedBy || ''; meta.revision = doc.revision; } catch (e) { /* optional */ }
  free(doc);
  return { bodies, units, kernel: 'rhino3dm (openNURBS)', version: k.version, warnings, meta, stats: { deflection: { type: 'relative', linear: lin, angular: null } } };
}

// ---------------------------------------------------------------- web-ifc: IFC
function ifcBodies(k, bytes, o) {
  const api = k.api, ang = o.angularDeflection ?? 0.2, seg = max(12, min(180, o.circleSegments ?? ceil(2 * PI / ang))), warnings = [], bodies = [];
  let id = -1; try { id = api.OpenModel(bytes, { COORDINATE_TO_ORIGIN: false, CIRCLE_SEGMENTS: seg }); } catch (e) { throw new Error(`web-ifc could not open the IFC model (${e && e.message ? e.message : e})`); }
  if (!(id >= 0)) throw new Error('web-ifc could not open the IFC model');
  const meta = { schema: '', elementTypes: {} }; let failedEl = 0;
  try {
    try { meta.schema = api.GetModelSchema(id) || ''; } catch (e) { meta.schema = ''; }
    api.StreamAllMeshes(id, fm => {
      const pgs = fm.geometries, list = []; let colour = null;
      for (let i = 0; i < pgs.size(); i++) {
        const pg = pgs.get(i), T = pg.flatTransformation; let g = null;
        try { g = api.GetGeometry(id, pg.geometryExpressID); const v = api.GetVertexArray(g.GetVertexData(), g.GetVertexDataSize()), ix = api.GetIndexArray(g.GetIndexData(), g.GetIndexDataSize()), P = new Float64Array(v.length / 2);
          // 4 × 4 placement (column-major), then web-ifc's y-up axes back to the IFC z-up convention: (x, y, z) ← (x, −z, y)
          for (let s = 0, d = 0; s + 5 < v.length; s += 6, d += 3) { const x = v[s], y = v[s + 1], z = v[s + 2], X = T[0] * x + T[4] * y + T[8] * z + T[12], Y = T[1] * x + T[5] * y + T[9] * z + T[13], Z = T[2] * x + T[6] * y + T[10] * z + T[14]; P[d] = X; P[d + 1] = -Z; P[d + 2] = Y; }
          const det = T[0] * (T[5] * T[10] - T[6] * T[9]) - T[4] * (T[1] * T[10] - T[2] * T[9]) + T[8] * (T[1] * T[6] - T[2] * T[5]), I = Array.from(ix); if (det < 0) for (let t = 0; t + 2 < I.length; t += 3) { const s = I[t + 1]; I[t + 1] = I[t + 2]; I[t + 2] = s; }
          if (I.length) list.push({ positions: P, indices: I }); if (!colour && pg.color) colour = [pg.color.x, pg.color.y, pg.color.z];
        } catch (e) { failedEl++; } finally { free(g); }
      }
      if (!list.length) return;
      let type = '', name = ''; try { type = api.GetNameFromTypeCode(api.GetLineType(id, fm.expressID)) || ''; } catch (e) { type = ''; }
      try { const L = api.GetLine(id, fm.expressID); name = L && L.Name && L.Name.value ? String(L.Name.value) : ''; } catch (e) { name = ''; }
      meta.elementTypes[type || 'unknown'] = (meta.elementTypes[type || 'unknown'] || 0) + 1;
      let nv = 0, ni = 0; for (const a of list) { nv += a.positions.length; ni += a.indices.length; } const P = new Float64Array(nv), I = new Array(ni); let pv = 0, pi = 0; for (const a of list) { const base = pv / 3; P.set(a.positions, pv); pv += a.positions.length; for (let i = 0; i < a.indices.length; i++) I[pi++] = a.indices[i] + base; }
      bodies.push({ name: name || `${type || 'element'} #${fm.expressID}`, source: `#${fm.expressID}`, type, colour, positions: P, indices: I });
    });
  } finally { try { api.CloseModel(id); } catch (e) { /* already closed */ } }
  if (failedEl) warnings.push(`${failedEl} placed geometry item(s) raised an error in the kernel and are missing`);
  return { bodies, units: 'm', kernel: 'web-ifc', version: k.version, warnings, meta, stats: { deflection: { type: 'segments', linear: null, angular: 2 * PI / seg, circleSegments: seg } } };
}

// tessellate({format: 'step' | 'iges' | 'brep' | '3dm' | 'ifc', bytes, options}) →
//   {mesh: {kind: 'mesh', vertices, faces, units}, parts: [{name, colour: [r, g, b] | null, faceRange: [start, end), triangles,
//    area, closed, volume (closed bodies only), shells, brepFaces?, type?, source?}], units, warnings: [], stats: {bodies, solids,
//    shells, faces, untessellatedFaces, triangles, vertices, seconds, kernel, kernelVersion, deflection}}.
// options: linearDeflection (default 0.001; a fraction of the mean bounding-box edge, or a length in millimetres when
// deflectionType is 'absolute'), angularDeflection in radians (default 0.5; for IFC it sets the segments per full circle,
// default 0.2 rad → 32). `solids` counts bodies whose welded triangulation is closed and consistently oriented, `shells`
// connected sets of triangles, `faces` the B-rep faces the kernel reported. STEP and IGES come back in millimetres whatever
// unit the file declares, IFC in metres, 3DM in the model unit, BREP in the file's own (undeclared) unit.
export async function tessellate({ format, bytes, options = {} } = {}) {
  const fmt = String(format || '').toLowerCase().replace(/^\./, ''), f = ({ stp: 'step', p21: 'step', igs: 'iges', brp: 'brep' })[fmt] || fmt, id = ({ step: 'occt', iges: 'occt', brep: 'occt', '3dm': 'rhino', ifc: 'ifc' })[f];
  if (!id) throw new Error(`tessellate: no kernel handles the format "${format}"`);
  if (!bytes || !bytes.length) throw new Error('tessellate: no data');
  const k = await loadKernel(id), t0 = now(), b = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes.buffer || bytes, bytes.byteOffset || 0, bytes.byteLength);
  let raw; try { raw = id === 'occt' ? occtBodies(k, f, b, options) : id === 'rhino' ? rhinoBodies(k, b, options) : ifcBodies(k, b, options); }
  catch (e) {
    // a trapped or aborted WebAssembly instance cannot be reused: drop it so that the next import starts a fresh one
    if ((globalThis.WebAssembly && e instanceof globalThis.WebAssembly.RuntimeError) || /Aborted|out of memory|unreachable|memory access/i.test(String(e && e.message))) delete loaded[id];
    if (e instanceof Error && !/^\d+$/.test(e.message)) throw e; throw new Error(`${KERNELS[id].name} stopped while reading the ${f.toUpperCase()} data (${e && e.message ? e.message : e})`); }
  return assemble(raw, t0);
}
