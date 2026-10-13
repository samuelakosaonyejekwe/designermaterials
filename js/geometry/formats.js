// Writers and parsers for non-atomic geometry formats (surface meshes, CAD exchange, FE/CFD meshes, lithography
// layouts, voxel volumes, graphs, point clouds, toolpaths and procedural descriptions). Pure ES module: no DOM, no
// file system, no code execution — procedural formats are parsed, never run.
// Parsers return a geometry object (CONTRACT schemas) that may carry `meta` (header data), `notes` (assumptions made)
// and `lost` (information present in the file that was not preserved). They throw Error with a specific message on
// malformed input.
import { sdf, sdfToVoxels, marchingCubes } from './implicit.js';
import { triangulate, otsu } from './mesh.js';
import { pathToPolygon, polyArea } from './layout.js';

const { abs, min, max, floor, round, sqrt, PI, cos, sin, hypot } = Math;
const ENC = new TextEncoder();
export const toBytes = x => typeof x === 'string' ? ENC.encode(x) : x instanceof Uint8Array ? x : new Uint8Array(x.buffer || x, x.byteOffset || 0, x.byteLength);
export function toText(x, limit) {
  if (typeof x === 'string') return limit ? x.slice(0, limit) : x;
  const b = limit && x.length > limit ? x.subarray(0, limit) : x;
  if (b.length >= 2 && ((b[0] === 0xFF && b[1] === 0xFE) || (b[0] === 0xFE && b[1] === 0xFF))) return new TextDecoder(b[0] === 0xFF ? 'utf-16le' : 'utf-16be').decode(b);
  return new TextDecoder('utf-8').decode(b);
}
const latin = (b, o, n) => { let s = ''; for (let i = o, e = min(b.length, o + n); i < e; i++) s += String.fromCharCode(b[i]); return s; };
const N9 = v => { const r = +Number(v).toPrecision(9); return r === 0 ? '0' : String(r); };
const nums = s => { const o = []; for (const t of s.split(/[\s,;]+/)) if (t !== '') o.push(Number(t)); return o; };
const lines = t => t.split(/\r\n|\r|\n/);
const view = b => new DataView(b.buffer, b.byteOffset, b.byteLength);
const fail = (fmt, msg) => { throw new Error(`${fmt}: ${msg}`); };
const UNIT_M = { m: 1, mm: 1e-3, cm: 1e-2, um: 1e-6, 'µm': 1e-6, 'μm': 1e-6, micron: 1e-6, nm: 1e-9, 'Å': 1e-10, A: 1e-10, angstrom: 1e-10, in: 0.0254, inch: 0.0254, ft: 0.3048 };
export const unitToMetres = u => UNIT_M[u] ?? UNIT_M[String(u || '').toLowerCase()] ?? null;

class Writer {
  constructor(n = 1024) { this.b = new Uint8Array(n); this.n = 0; this.v = new DataView(this.b.buffer); }
  need(k) { if (this.n + k > this.b.length) { const nb = new Uint8Array(max(this.b.length * 2, this.n + k)); nb.set(this.b); this.b = nb; this.v = new DataView(nb.buffer); } }
  u8(v) { this.need(1); this.b[this.n++] = v; return this; }
  u16(v, le = true) { this.need(2); this.v.setUint16(this.n, v, le); this.n += 2; return this; }
  i16(v, le = true) { this.need(2); this.v.setInt16(this.n, v, le); this.n += 2; return this; }
  u32(v, le = true) { this.need(4); this.v.setUint32(this.n, v, le); this.n += 4; return this; }
  i32(v, le = true) { this.need(4); this.v.setInt32(this.n, v, le); this.n += 4; return this; }
  u64(v, le = true) { const lo = v % 4294967296, hi = floor(v / 4294967296); if (le) { this.u32(lo); this.u32(hi); } else { this.u32(hi, false); this.u32(lo, false); } return this; }
  f32(v, le = true) { this.need(4); this.v.setFloat32(this.n, v, le); this.n += 4; return this; }
  f64(v, le = true) { this.need(8); this.v.setFloat64(this.n, v, le); this.n += 8; return this; }
  bytes(a) { this.need(a.length); this.b.set(a, this.n); this.n += a.length; return this; }
  str(s) { return this.bytes(ENC.encode(s)); }
  zero(k) { this.need(k); this.n += k; return this; }
  align(k) { const r = this.n % k; return r ? this.zero(k - r) : this; }
  out() { return this.b.slice(0, this.n); }
}

// ---------------------------------------------------------------- CRC-32, inflate, zip
const CRC_T = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
export function crc32(b) { let c = 0xFFFFFFFF; for (let i = 0; i < b.length; i++) c = CRC_T[(c ^ b[i]) & 255] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }

const LBASE = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258], LEXT = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DBASE = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577], DEXT = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
const CLORD = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
function huff(lens) { const count = new Uint16Array(16), offs = new Uint16Array(16), symbol = new Uint16Array(lens.length); for (const l of lens) count[l]++; count[0] = 0; for (let i = 1; i < 16; i++) offs[i] = offs[i - 1] + count[i - 1]; for (let s = 0; s < lens.length; s++) if (lens[s]) symbol[offs[lens[s]]++] = s; return { count, symbol }; }
let FIXED = null;
// DEFLATE (RFC 1951) decoder. Returns {data, end} where end is the offset just after the compressed stream.
export function inflateRaw(src, start = 0, sizeHint = 0) {
  let pos = start, bb = 0, bc = 0, out = new Uint8Array(min(sizeHint, (src.length - start) * 1032 + 1024) || max(4096, (src.length - start) * 4)), op = 0; // DEFLATE cannot expand by more than ~1032:1, so a larger declared size is not trusted
  const grow = n => { if (op + n > out.length) { const b = new Uint8Array(max(out.length * 2, op + n)); b.set(out); out = b; } };
  const bits = n => { while (bc < n) { if (pos >= src.length) fail('inflate', 'unexpected end of compressed data'); bb |= src[pos++] << bc; bc += 8; } const v = bb & ((1 << n) - 1); bb >>>= n; bc -= n; return v; };
  const sym = h => { let code = 0, first = 0, index = 0; for (let len = 1; len < 16; len++) { code |= bits(1); const c = h.count[len]; if (code - c < first) return h.symbol[index + (code - first)]; index += c; first += c; first <<= 1; code <<= 1; } fail('inflate', 'invalid Huffman code'); };
  for (let last = 0; !last;) {
    last = bits(1); const type = bits(2);
    if (type === 0) {
      bb = 0; bc = 0; const len = src[pos] | src[pos + 1] << 8; pos += 4;
      if (pos + len > src.length) fail('inflate', 'stored block exceeds input');
      grow(len); out.set(src.subarray(pos, pos + len), op); op += len; pos += len;
    } else if (type === 1 || type === 2) {
      let lh, dh;
      if (type === 1) { if (!FIXED) { const l = new Uint8Array(288); l.fill(8, 0, 144); l.fill(9, 144, 256); l.fill(7, 256, 280); l.fill(8, 280, 288); FIXED = [huff(l), huff(new Uint8Array(30).fill(5))]; } [lh, dh] = FIXED; }
      else {
        const nl = bits(5) + 257, nd = bits(5) + 1, nc = bits(4) + 4, cl = new Uint8Array(19);
        for (let i = 0; i < nc; i++) cl[CLORD[i]] = bits(3);
        const ch = huff(cl), lens = new Uint8Array(nl + nd);
        for (let i = 0; i < nl + nd;) { const s = sym(ch); if (s < 16) lens[i++] = s; else { let rep, v = 0; if (s === 16) { if (!i) fail('inflate', 'repeat with no previous length'); v = lens[i - 1]; rep = 3 + bits(2); } else rep = s === 17 ? 3 + bits(3) : 11 + bits(7); if (i + rep > nl + nd) fail('inflate', 'too many code lengths'); while (rep--) lens[i++] = v; } }
        lh = huff(lens.subarray(0, nl)); dh = huff(lens.subarray(nl));
      }
      for (; ;) {
        const s = sym(lh);
        if (s < 256) { grow(1); out[op++] = s; } else if (s === 256) break; else {
          if (s > 285) fail('inflate', 'invalid length symbol');
          const len = LBASE[s - 257] + bits(LEXT[s - 257]), ds = sym(dh); if (ds > 29) fail('inflate', 'invalid distance symbol');
          const d = DBASE[ds] + bits(DEXT[ds]); if (d > op) fail('inflate', 'distance beyond start of output');
          grow(len); for (let i = 0; i < len; i++, op++) out[op] = out[op - d];
        }
      }
    } else fail('inflate', 'invalid block type');
  }
  return { data: out.subarray(0, op), end: pos };
}
export function gunzip(b) {
  if (b[0] !== 0x1F || b[1] !== 0x8B) fail('gzip', 'missing gzip signature'); if (b[2] !== 8) fail('gzip', 'unsupported compression method');
  const flg = b[3]; let p = 10;
  if (flg & 4) p += 2 + (b[p] | b[p + 1] << 8); if (flg & 8) { while (b[p++]); } if (flg & 16) { while (b[p++]); } if (flg & 2) p += 2;
  return inflateRaw(b, p).data;
}
export const unzlib = (b, hint) => inflateRaw(b, 2, hint).data;

// ZIP writer (stored, no compression) with CRC-32 — sufficient for OPC packages such as 3MF.
export function zipStore(files) {
  const w = new Writer(4096), central = [];
  for (const f of files) {
    const name = ENC.encode(f.name), data = toBytes(f.data), crc = crc32(data), off = w.n;
    w.u32(0x04034b50).u16(20).u16(0x0800).u16(0).u16(0).u16(0x21).u32(crc).u32(data.length).u32(data.length).u16(name.length).u16(0).bytes(name).bytes(data);
    central.push({ name, crc, size: data.length, off });
  }
  const cd = w.n;
  for (const c of central) w.u32(0x02014b50).u16(20).u16(20).u16(0x0800).u16(0).u16(0).u16(0x21).u32(c.crc).u32(c.size).u32(c.size).u16(c.name.length).u16(0).u16(0).u16(0).u16(0).u32(0).u32(c.off).bytes(c.name);
  const cdSize = w.n - cd;
  w.u32(0x06054b50).u16(0).u16(0).u16(central.length).u16(central.length).u32(cdSize).u32(cd).u16(0);
  return w.out();
}
// ZIP reader: central directory, stored (0) and DEFLATE (8) entries. `only` (RegExp) limits which entries are inflated.
// Optional decoders registered on demand (e.g. Zstandard for recent Fusion archives).
export const codecs = {};
export function unzip(bytes, only) {
  const b = toBytes(bytes), v = view(b); let e = -1;
  for (let i = b.length - 22; i >= max(0, b.length - 66000); i--) if (v.getUint32(i, true) === 0x06054b50) { e = i; break; }
  if (e < 0) fail('zip', 'end-of-central-directory record not found (truncated or not a zip archive)');
  const n = v.getUint16(e + 10, true); let p = v.getUint32(e + 16, true); const out = [];
  if (p === 0xFFFFFFFF) fail('zip', 'ZIP64 archives are not supported');
  for (let i = 0; i < n && p + 46 <= b.length; i++) {
    if (v.getUint32(p, true) !== 0x02014b50) break;
    const method = v.getUint16(p + 10, true), csize = v.getUint32(p + 20, true), usize = v.getUint32(p + 24, true), nl = v.getUint16(p + 28, true), xl = v.getUint16(p + 30, true), cl = v.getUint16(p + 32, true), lo = v.getUint32(p + 42, true);
    const name = toText(b.subarray(p + 46, p + 46 + nl)), entry = { name, size: usize, method, data: null };
    if (!only || only.test(name)) {
      const ds = lo + 30 + v.getUint16(lo + 26, true) + v.getUint16(lo + 28, true);
      if (method === 0) entry.data = b.subarray(ds, ds + csize); else if (method === 8) entry.data = inflateRaw(b, ds, usize || 0).data; else if (method === 93 && codecs.zstd) { try { entry.data = codecs.zstd(b.subarray(ds, ds + csize)); } catch (err) { entry.error = 'Zstandard stream could not be decoded: ' + err.message; } } else entry.error = `compression method ${method} not supported`;
    }
    out.push(entry); p += 46 + nl + xl + cl;
  }
  return out;
}

// ---------------------------------------------------------------- shared conversion helpers
const RX = new Map();
const attr = (s, name) => { let r = RX.get(name); if (!r) RX.set(name, r = new RegExp('(?:^|[\\s<])' + name + '\\s*=\\s*(?:"([^"]*)"|\'([^\']*)\')')); const m = r.exec(s); return m ? (m[1] ?? m[2]) : undefined; };
const xmlText = s => s.replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&amp;/g, '&');
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const fan = (idx, faces) => { for (let i = 1; i < idx.length - 1; i++) faces.push([idx[0], idx[i], idx[i + 1]]); };
const b64 = s => Uint8Array.from(atob(s.replace(/\s+/g, '')), c => c.charCodeAt(0));
const toB64 = b => { let s = ''; for (let i = 0; i < b.length; i += 0x8000) s += String.fromCharCode.apply(null, b.subarray(i, i + 0x8000)); return btoa(s); };

// Typed reads from a byte buffer.
const DT = { uint8: [1, 'getUint8'], int8: [1, 'getInt8'], uint16: [2, 'getUint16'], int16: [2, 'getInt16'], uint32: [4, 'getUint32'], int32: [4, 'getInt32'], float32: [4, 'getFloat32'], float64: [8, 'getFloat64'], uint64: [8, 'getBigUint64'], int64: [8, 'getBigInt64'] };
export function readTyped(b, offset, count, dtype, le = true) {
  const d = DT[dtype]; if (!d) fail('binary data', `unsupported data type "${dtype}"`);
  if (offset + count * d[0] > b.length) fail('binary data', `expected ${count} × ${dtype} (${count * d[0]} bytes) at offset ${offset}, but only ${b.length - offset} bytes are present`);
  if (dtype === 'uint8') return b.subarray(offset, offset + count);
  const v = view(b), out = dtype === 'float64' || d[0] === 8 ? new Float64Array(count) : dtype === 'float32' ? new Float32Array(count) : new Float64Array(count), big = d[0] === 8 && dtype !== 'float64';
  for (let i = 0, p = offset; i < count; i++, p += d[0]) out[i] = big ? Number(v[d[1]](p, le)) : v[d[1]](p, le);
  return out;
}
// Scalar array → voxels. Small non-negative integer label sets are kept as phase ids; anything else is treated as a
// grey-level field: it is rescaled to 8 bits (kept as `grey`) and binarised with Otsu's threshold (phase 1 = brighter).
export function scalarToVoxels(a, nx, ny, nz, { voxel = 1, units = '', origin, spacing, threshold } = {}) {
  const N = nx * ny * nz; if (!(N > 0)) fail('volume', `invalid dimensions ${nx} × ${ny} × ${nz}`);
  if (a.length < N) fail('volume', `data holds ${a.length} values but ${nx} × ${ny} × ${nz} = ${N} are required`);
  let lo = Infinity, hi = -Infinity, integer = true; const seen = new Uint8Array(256); let distinct = 0;
  for (let i = 0; i < N; i++) { const v = a[i]; if (v < lo) lo = v; if (v > hi) hi = v; if (integer) { if (v !== (v | 0) || v < 0 || v > 255) integer = false; else if (!seen[v]) { seen[v] = 1; distinct++; } } }
  const out = { kind: 'voxels', nx, ny, nz, data: null, voxel, units, notes: [] };
  if (origin) out.origin = origin; if (spacing) { out.spacing = spacing; if (abs(spacing[0] - spacing[1]) > 1e-9 * abs(spacing[0]) || abs(spacing[0] - spacing[2]) > 1e-9 * abs(spacing[0])) out.notes.push(`anisotropic voxel spacing ${spacing.join(' × ')}: \`voxel\` holds the x spacing, \`spacing\` the full triple`); }
  if (integer && distinct <= 32 && threshold === undefined) { out.data = a instanceof Uint8Array && a.length === N ? a.slice() : Uint8Array.from(a.subarray ? a.subarray(0, N) : a.slice(0, N)); out.valueMode = 'labels'; out.phases = distinct; }
  else {
    // 8-bit data are kept as stored; wider or floating-point data are rescaled linearly from [min, max] to 0–255
    const g = new Uint8Array(N), keep = integer, sc = keep ? 1 : hi > lo ? 255 / (hi - lo) : 0, l0 = keep ? 0 : lo; for (let i = 0; i < N; i++) g[i] = (a[i] - l0) * sc;
    const thr = threshold !== undefined ? max(0, min(255, (threshold - l0) * sc)) : otsu(g), d = new Uint8Array(N); for (let i = 0; i < N; i++) d[i] = g[i] > thr ? 1 : 0;
    Object.assign(out, { data: d, grey: g, greyScale: keep ? 'as stored (8-bit)' : `linear from [${N9(lo)}, ${N9(hi)}]`, valueMode: 'grey', valueRange: [lo, hi], threshold: sc ? l0 + thr / sc : lo });
    out.notes.push(`grey-level data (${lo} … ${hi}) binarised at ${N9(out.threshold)} (${threshold !== undefined ? 'user threshold' : 'Otsu'}); phase 1 = brighter. The 8-bit field is kept as \`grey\` for re-segmentation`);
  }
  return out;
}
const TETF = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]], HEXF = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]], WEDF = [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]], PYRF = [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
const CELLF = { tet: TETF, hex: HEXF, wedge: WEDF, pyr: PYRF }, CORNERS = { line: 2, tri: 3, quad: 4, tet: 4, hex: 8, wedge: 6, pyr: 5 };
// FE/CFD connectivity → geometry. Shell elements become triangles; solid elements contribute their boundary faces
// (faces used by exactly one cell); beam/truss-only models become a lattice. Mid-side nodes are dropped.
export function femToGeometry(pts, elems, { units = '', sets, meta } = {}) {
  const counts = {}, faces = [], edges = [], fmap = new Map(), notes = [], lost = []; let higher = 0;
  const key = f => [...f].sort((a, b) => a - b).join(',');
  for (const e of elems) {
    counts[e.type] = (counts[e.type] || 0) + 1; const n = e.n; if (n.length > CORNERS[e.type]) higher++;
    if (e.type === 'line') edges.push([n[0], n[1]]); else if (e.type === 'tri') faces.push([n[0], n[1], n[2]]); else if (e.type === 'quad') faces.push([n[0], n[1], n[2]], [n[0], n[2], n[3]]);
    else for (const f of CELLF[e.type]) { const g = f.map(i => n[i]), k = key(g), o = fmap.get(k); if (o) o.c++; else fmap.set(k, { c: 1, g }); }
  }
  let nb = 0; for (const { c, g } of fmap.values()) if (c === 1) { nb++; fan(g, faces); }
  if (higher) lost.push(`${higher} higher-order elements reduced to their corner nodes (mid-side nodes dropped)`);
  const solid = (counts.tet || 0) + (counts.hex || 0) + (counts.wedge || 0) + (counts.pyr || 0);
  if (solid) notes.push(`${solid} solid elements: the displayed surface is their boundary (${nb} exterior faces); interior connectivity is summarised in elementCounts`);
  const used = new Int32Array(pts.length).fill(-1), V = [], ix = i => { if (i === undefined || i < 0 || i >= pts.length) fail('mesh', 'element references a node that is not defined'); if (used[i] < 0) { used[i] = V.length; V.push(pts[i]); } return used[i]; };
  if (!faces.length && edges.length) return { kind: 'lattice', nodes: V, edges: edges.map(e => [ix(e[0]), ix(e[1])]), radius: 0, units, elementCounts: counts, nNodesFile: pts.length, sets, meta, notes: [...notes, 'beam/truss model: strut radius is not stored with the connectivity (radius set to 0; assign a section)'], lost };
  const F = faces.map(f => [ix(f[0]), ix(f[1]), ix(f[2])]);
  let vol = 0; for (const f of F) { const A = V[f[0]], B = V[f[1]], C = V[f[2]]; vol += A[0] * (B[1] * C[2] - B[2] * C[1]) + A[1] * (B[2] * C[0] - B[0] * C[2]) + A[2] * (B[0] * C[1] - B[1] * C[0]); }
  if (solid && vol < 0) for (const f of F) { const t = f[1]; f[1] = f[2]; f[2] = t; }
  const out = { kind: 'mesh', vertices: V, faces: F, units, elementCounts: counts, nNodesFile: pts.length, solidElements: solid, sets, meta, notes, lost };
  if (edges.length) out.beams = edges.map(e => [ix(e[0]), ix(e[1])]);
  if (!F.length && !edges.length) notes.push('no supported elements found (nodes only)');
  return out;
}
const needMesh = (g, fmt) => { if (!g || g.kind !== 'mesh' || !g.vertices || !g.faces) fail(fmt, `needs a mesh geometry (got ${g && g.kind})`); return g; };
const faceNormal = (A, B, C) => { const u = [B[0] - A[0], B[1] - A[1], B[2] - A[2]], v = [C[0] - A[0], C[1] - A[1], C[2] - A[2]], n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], l = hypot(n[0], n[1], n[2]) || 1; return [n[0] / l, n[1] / l, n[2] / l]; };
function weld(tri) { // flat xyz triples → indexed mesh (exact coordinate match)
  const map = new Map(), V = [], F = []; let f = [];
  for (let i = 0; i < tri.length; i += 3) { const k = tri[i] + ',' + tri[i + 1] + ',' + tri[i + 2]; let j = map.get(k); if (j === undefined) { j = V.length; V.push([tri[i], tri[i + 1], tri[i + 2]]); map.set(k, j); } f.push(j); if (f.length === 3) { F.push(f); f = []; } }
  return { V, F };
}

// ---------------------------------------------------------------- STL
export function toSTL(mesh, { binary = true, name = 'designer_materials' } = {}) {
  const { vertices: V, faces: F } = needMesh(mesh, 'STL');
  if (!binary) { const o = [`solid ${name}`]; for (const f of F) { const n = faceNormal(V[f[0]], V[f[1]], V[f[2]]); o.push(` facet normal ${N9(n[0])} ${N9(n[1])} ${N9(n[2])}`, '  outer loop'); for (const i of f) o.push(`   vertex ${N9(V[i][0])} ${N9(V[i][1])} ${N9(V[i][2] || 0)}`); o.push('  endloop', ' endfacet'); } o.push(`endsolid ${name}`); return o.join('\n') + '\n'; }
  const w = new Writer(84 + 50 * F.length); w.str(`Binary STL, units ${mesh.units || 'unspecified'}`.slice(0, 79)); w.n = 80; w.u32(F.length);
  for (const f of F) { const n = faceNormal(V[f[0]], V[f[1]], V[f[2]]); w.f32(n[0]).f32(n[1]).f32(n[2]); for (const i of f) w.f32(V[i][0]).f32(V[i][1]).f32(V[i][2] || 0); w.u16(0); }
  return w.out();
}
// Binary vs ASCII is decided by the size rule 84 + 50·n = file length, never by the leading word "solid".
export function parseSTL(input) {
  const b = toBytes(input), notes = ['STL carries no units: lengths are taken as stored (commonly mm)'], lost = [];
  let tri, enc, name = '';
  const nb = b.length >= 84 ? view(b).getUint32(80, true) : -1, exact = nb >= 0 && 84 + 50 * nb === b.length;
  const head = latin(b, 0, 512), looksAscii = /^\s*solid\b/.test(head) && /facet|endsolid/.test(latin(b, 0, 4096));
  if (exact || (!looksAscii && nb > 0 && 84 + 50 * nb <= b.length)) {
    const v = view(b); tri = new Float32Array(nb * 9); enc = 'binary'; name = latin(b, 0, 80).replace(/\0.*$/, '').trim();
    for (let i = 0, p = 84; i < nb; i++, p += 50) for (let k = 0; k < 9; k++) tri[i * 9 + k] = v.getFloat32(p + 12 + 4 * k, true);
    if (!exact) notes.push(`${b.length - 84 - 50 * nb} trailing bytes after the declared ${nb} facets were ignored`);
    lost.push('per-facet attribute bytes (colour extensions) ignored');
  } else {
    const t = toText(b), re = /vertex\s+(\S+)\s+(\S+)\s+(\S+)/g, a = []; let m; enc = 'ASCII';
    while ((m = re.exec(t))) a.push(+m[1], +m[2], +m[3]);
    if (!a.length) fail('STL', 'neither a binary STL (size ≠ 84 + 50·n) nor an ASCII STL with "vertex" records');
    if (a.length % 9) fail('STL', `ASCII file has ${a.length / 3} vertices, not a multiple of 3 (truncated facet)`);
    if (a.some(x => !Number.isFinite(x))) fail('STL', 'non-numeric vertex coordinate');
    tri = Float32Array.from(a); name = (/^\s*solid[ \t]+([^\r\n]*)/.exec(t) || [])[1] || '';
    const ns = (t.match(/^\s*solid\b/gm) || []).length; if (ns > 1) notes.push(`${ns} solids merged into one mesh`);
  }
  const { V, F } = weld(tri);
  notes.push(`facet soup of ${F.length} triangles welded to ${V.length} shared vertices (exact coordinate match)`);
  return { kind: 'mesh', vertices: V, faces: F, units: 'mm', meta: { encoding: enc, name, facets: F.length }, notes, lost: [...lost, 'stored facet normals discarded (recomputed from winding)'] };
}

// ---------------------------------------------------------------- OBJ / PLY / OFF / GTS
export function toOBJ(mesh) { const { vertices: V, faces: F } = needMesh(mesh, 'OBJ'); return `# units: ${mesh.units || 'unspecified'}\no mesh\n` + V.map(p => `v ${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n') + '\n' + F.map(f => `f ${f[0] + 1} ${f[1] + 1} ${f[2] + 1}`).join('\n') + '\n'; }
export function parseOBJ(input) {
  const V = [], F = [], groups = [], L = []; let cur = null, polys = 0, curves = 0, mtl = null;
  for (const raw of lines(toText(input))) {
    const l = raw.trim(); if (!l || l[0] === '#') continue; const t = l.split(/\s+/), k = t[0];
    if (k === 'v') { const p = [+t[1], +t[2], +t[3] || 0]; if (!Number.isFinite(p[0] + p[1] + p[2])) fail('OBJ', `bad vertex line "${l}"`); V.push(p); }
    else if (k === 'f') { const idx = t.slice(1).map(s => { const i = parseInt(s, 10); return i < 0 ? V.length + i : i - 1; }); if (idx.length < 3) continue; if (idx.some(i => !(i >= 0 && i < V.length))) fail('OBJ', `face references an undefined vertex: "${l}"`); if (idx.length > 3) polys++; const f0 = F.length; fan(idx, F); if (cur) cur.faces += F.length - f0; }
    else if (k === 'l') { const idx = t.slice(1).map(s => { const i = parseInt(s, 10); return i < 0 ? V.length + i : i - 1; }); for (let i = 0; i + 1 < idx.length; i++) L.push([idx[i], idx[i + 1]]); }
    else if (k === 'g' || k === 'o') groups.push(cur = { name: t.slice(1).join(' '), faces: 0, material: mtl });
    else if (k === 'usemtl') { mtl = t[1]; if (cur) cur.material = mtl; }
    else if (k === 'curv' || k === 'surf' || k === 'cstype') curves++;
  }
  if (!V.length) fail('OBJ', 'no vertex ("v") records found');
  const notes = ['OBJ carries no units: lengths taken as stored'], lost = ['texture coordinates, vertex normals and material library (.mtl) not imported'];
  if (polys) notes.push(`${polys} polygons with more than three vertices fan-triangulated`); if (curves) lost.push('free-form curve/surface statements ignored');
  if (!F.length && L.length) return { kind: 'lattice', nodes: V, edges: L, radius: 0, units: '', notes: [...notes, 'line elements only: imported as a node–edge network'], lost };
  if (!F.length) return { kind: 'lattice', nodes: V, edges: [], radius: 0, units: '', pointCloud: true, notes: [...notes, 'no faces: imported as a point cloud'], lost };
  return { kind: 'mesh', vertices: V, faces: F, units: '', groups: groups.filter(g => g.faces), notes, lost };
}
export function toPLY(mesh, { binary = false } = {}) {
  const { vertices: V, faces: F } = needMesh(mesh, 'PLY');
  const head = `ply\nformat ${binary ? 'binary_little_endian' : 'ascii'} 1.0\ncomment units ${mesh.units || 'unspecified'}\nelement vertex ${V.length}\nproperty float64 x\nproperty float64 y\nproperty float64 z\nelement face ${F.length}\nproperty list uchar int vertex_indices\nend_header\n`;
  if (!binary) return head + V.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n') + '\n' + F.map(f => `3 ${f[0]} ${f[1]} ${f[2]}`).join('\n') + '\n';
  const w = new Writer(head.length + V.length * 24 + F.length * 13); w.str(head); for (const p of V) w.f64(p[0]).f64(p[1]).f64(p[2] || 0); for (const f of F) w.u8(3).i32(f[0]).i32(f[1]).i32(f[2]); return w.out();
}
const PLYT = { char: 'int8', int8: 'int8', uchar: 'uint8', uint8: 'uint8', short: 'int16', int16: 'int16', ushort: 'uint16', uint16: 'uint16', int: 'int32', int32: 'int32', uint: 'uint32', uint32: 'uint32', float: 'float32', float32: 'float32', double: 'float64', float64: 'float64' };
export function parsePLY(input) {
  const b = toBytes(input), head = latin(b, 0, min(b.length, 65536)), m = /end_header\r?\n/.exec(head);
  if (!/^ply\r?\n/.test(head) || !m) fail('PLY', 'missing "ply" magic or "end_header"');
  const els = []; let fmt = '';
  for (const l of lines(head.slice(0, m.index))) { const t = l.trim().split(/\s+/); if (t[0] === 'format') fmt = t[1]; else if (t[0] === 'element') els.push({ name: t[1], count: +t[2], props: [] }); else if (t[0] === 'property' && els.length) els[els.length - 1].props.push(t[1] === 'list' ? { list: PLYT[t[2]], type: PLYT[t[3]], name: t[4] } : { type: PLYT[t[1]], name: t[2] }); }
  if (!['ascii', 'binary_little_endian', 'binary_big_endian'].includes(fmt)) fail('PLY', `unknown format "${fmt}"`);
  const V = [], F = []; let polys = 0, colours = false;
  const take = (el, rec) => {
    if (el.name === 'vertex') V.push([rec.x, rec.y, rec.z || 0]);
    else if (el.name === 'face' || el.name === 'tristrips') { const idx = rec.vertex_indices || rec.vertex_index; if (idx && idx.length >= 3) { if (idx.length > 3) polys++; fan(idx, F); } }
  };
  if (fmt === 'ascii') {
    const tk = toText(b.subarray(m.index + m[0].length)).split(/\s+/); let p = tk[0] === '' ? 1 : 0;
    for (const el of els) for (let i = 0; i < el.count; i++) { if (p >= tk.length) fail('PLY', `the file ends before the declared ${el.count} ${el.name} records`); const rec = {}; for (const pr of el.props) { if (pr.list) { const n = +tk[p++]; if (!(n >= 0 && n <= tk.length - p)) fail('PLY', 'bad list length in an ASCII record'); const a = new Array(n); for (let k = 0; k < n; k++) a[k] = +tk[p++]; rec[pr.name] = a; } else rec[pr.name] = +tk[p++]; } take(el, rec); }
  } else {
    const le = fmt === 'binary_little_endian', v = view(b); let p = m.index + m[0].length;
    const rd = t => { const d = DT[t]; if (!d) fail('PLY', 'unsupported property type'); if (p + d[0] > b.length) fail('PLY', 'binary payload shorter than the header declares'); const x = v[d[1]](p, le); p += d[0]; return x; };
    for (const el of els) for (let i = 0; i < el.count; i++) { const rec = {}; for (const pr of el.props) { if (pr.list) { const n = rd(pr.list), a = new Array(n); for (let k = 0; k < n; k++) a[k] = rd(pr.type); rec[pr.name] = a; } else rec[pr.name] = rd(pr.type); } take(el, rec); }
  }
  const ve = els.find(e => e.name === 'vertex'); if (!ve || !V.length) fail('PLY', 'no vertex element');
  if (V.some(p => !Number.isFinite(p[0] + p[1] + p[2]))) fail('PLY', 'vertex element lacks numeric x, y, z properties');
  colours = ve.props.some(p => /red|diffuse/.test(p.name));
  const lost = []; if (colours) lost.push('vertex colours not imported'); if (ve.props.some(p => p.name === 'nx')) lost.push('vertex normals not imported');
  const notes = ['PLY carries no units: lengths taken as stored']; if (polys) notes.push(`${polys} polygons fan-triangulated`);
  const meta = { format: fmt, elements: els.map(e => `${e.name} × ${e.count}`) };
  if (!F.length) return { kind: 'lattice', nodes: V, edges: [], radius: 0, units: '', pointCloud: true, meta, notes: [...notes, 'no faces: imported as a point cloud'], lost };
  if (F.some(f => f.some(i => !(i >= 0 && i < V.length)))) fail('PLY', 'face references an undefined vertex');
  return { kind: 'mesh', vertices: V, faces: F, units: '', meta, notes, lost };
}
export function toOFF(mesh) { const { vertices: V, faces: F } = needMesh(mesh, 'OFF'); return `OFF\n${V.length} ${F.length} 0\n` + V.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n') + '\n' + F.map(f => `3 ${f[0]} ${f[1]} ${f[2]}`).join('\n') + '\n'; }
export function parseOFF(input) {
  const L = lines(toText(input)).map(l => l.replace(/#.*/, '').trim()).filter(Boolean);
  const h = /^(ST)?(C)?(N)?(4)?(n)?OFF\b(.*)$/.exec(L[0] || ''); if (!h) fail('OFF', 'missing OFF header keyword');
  let li = 1, c = nums(h[6]); if (c.length < 2) c = nums(L[li++] || ''); const [nv, nf] = c;
  if (!(nv >= 0 && nf >= 0) || L.length < li + nv + nf) fail('OFF', `header declares ${nv} vertices and ${nf} faces but the file has ${L.length - li} data lines`);
  const V = [], F = []; let polys = 0;
  for (let i = 0; i < nv; i++) { const a = nums(L[li++]); if (a.length < 3 || !Number.isFinite(a[0] + a[1] + a[2])) fail('OFF', `bad vertex line ${i + 1}`); V.push([a[0], a[1], a[2]]); }
  for (let i = 0; i < nf; i++) { const a = nums(L[li++]), n = a[0], idx = a.slice(1, 1 + n); if (idx.length < n || idx.some(k => !(k >= 0 && k < nv))) fail('OFF', `bad face line ${i + 1}`); if (n > 3) polys++; if (n >= 3) fan(idx, F); }
  const lost = []; if (h[2]) lost.push('per-vertex colours (COFF) not imported'); if (h[3]) lost.push('per-vertex normals (NOFF) not imported');
  return { kind: 'mesh', vertices: V, faces: F, units: '', meta: { variant: L[0].split(/\s/)[0] }, notes: ['OFF carries no units: lengths taken as stored', ...(polys ? [`${polys} polygons fan-triangulated`] : [])], lost };
}
export function parseGTS(input) {
  const L = lines(toText(input)).map(l => l.replace(/[#!].*/, '').trim()).filter(Boolean), [nv, ne, nf] = nums(L[0] || '');
  if (!(nv > 0 && ne >= 0 && nf >= 0) || L.length < 1 + nv + ne + nf) fail('GTS', 'header "nVertices nEdges nFaces" does not match the file length');
  const V = [], E = [], F = [];
  for (let i = 0; i < nv; i++) { const a = nums(L[1 + i]); V.push([a[0], a[1], a[2] || 0]); }
  for (let i = 0; i < ne; i++) { const a = nums(L[1 + nv + i]); E.push([a[0] - 1, a[1] - 1]); }
  for (let i = 0; i < nf; i++) { const a = nums(L[1 + nv + ne + i]), e1 = E[a[0] - 1], e2 = E[a[1] - 1]; if (!e1 || !e2) fail('GTS', `face ${i + 1} references an undefined edge`); const o = e2[0] === e1[0] || e2[0] === e1[1] ? e2[1] : e2[0]; F.push(e2.includes(e1[1]) ? [e1[0], e1[1], o] : [e1[1], e1[0], o]); }
  return { kind: 'mesh', vertices: V, faces: F, units: '', notes: ['GTS carries no units: lengths taken as stored'], lost: [] };
}

// ---------------------------------------------------------------- 3MF / AMF
const mat12 = s => { const a = nums(s || ''); return a.length === 12 ? a : null; };
const ap12 = (T, p) => T ? [p[0] * T[0] + p[1] * T[3] + p[2] * T[6] + T[9], p[0] * T[1] + p[1] * T[4] + p[2] * T[7] + T[10], p[0] * T[2] + p[1] * T[5] + p[2] * T[8] + T[11]] : p;
const mul12 = (A, B) => !A ? B : !B ? A : [0, 1, 2, 3].flatMap(r => [0, 1, 2].map(c => A[r * 3] * B[c] + A[r * 3 + 1] * B[3 + c] + A[r * 3 + 2] * B[6 + c] + (r === 3 ? B[9 + c] : 0)));
export function to3MF(mesh, { name = 'part' } = {}) {
  const { vertices: V, faces: F } = needMesh(mesh, '3MF'), u = ({ mm: 'millimeter', um: 'micron', 'µm': 'micron', cm: 'centimeter', m: 'meter', in: 'inch', ft: 'foot' })[mesh.units] || 'millimeter';
  const model = `<?xml version="1.0" encoding="UTF-8"?>\n<model unit="${u}" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">\n <metadata name="Application">Designer Materials</metadata>\n <resources>\n  <object id="1" name="${esc(name)}" type="model">\n   <mesh>\n    <vertices>\n` +
    V.map(p => `     <vertex x="${N9(p[0])}" y="${N9(p[1])}" z="${N9(p[2] || 0)}"/>`).join('\n') + '\n    </vertices>\n    <triangles>\n' + F.map(f => `     <triangle v1="${f[0]}" v2="${f[1]}" v3="${f[2]}"/>`).join('\n') + '\n    </triangles>\n   </mesh>\n  </object>\n </resources>\n <build>\n  <item objectid="1"/>\n </build>\n</model>\n';
  return zipStore([
    { name: '[Content_Types].xml', data: '<?xml version="1.0" encoding="UTF-8"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>' },
    { name: '_rels/.rels', data: '<?xml version="1.0" encoding="UTF-8"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>' },
    { name: '3D/3dmodel.model', data: model }]);
}
export function parse3MF(input) {
  const files = unzip(input, /\.model$|\.rels$/i), mf = files.find(f => /^\/?3D\/3dmodel\.model$/i.test(f.name)) || files.find(f => /\.model$/i.test(f.name));
  if (!mf || !mf.data) fail('3MF', 'archive holds no 3D model part (3D/3dmodel.model)');
  const x = toText(mf.data), unit = attr((/<model\b[^>]*>/.exec(x) || [''])[0], 'unit') || 'millimeter', objs = new Map(), lost = [], notes = [];
  const reO = /<object\b([^>]*)>([\s\S]*?)<\/object>/g; let m;
  while ((m = reO.exec(x))) {
    const o = { V: [], F: [], comp: [], name: attr(m[1], 'name') }, body = m[2]; let q;
    const rv = /<vertex\b([^>]*)\/?>/g; while ((q = rv.exec(body))) o.V.push([+attr(q[1], 'x'), +attr(q[1], 'y'), +attr(q[1], 'z')]);
    const rt = /<triangle\b([^>]*)\/?>/g; while ((q = rt.exec(body))) o.F.push([+attr(q[1], 'v1'), +attr(q[1], 'v2'), +attr(q[1], 'v3')]);
    const rc = /<component\b([^>]*)\/?>/g; while ((q = rc.exec(body))) o.comp.push({ id: attr(q[1], 'objectid'), T: mat12(attr(q[1], 'transform')) });
    objs.set(attr(m[1], 'id'), o);
  }
  if (!objs.size) fail('3MF', 'model part contains no <object> resources');
  const V = [], F = [], parts = [];
  const emit = (id, T, depth = 0) => { const o = objs.get(id); if (!o || depth > 8) return; if (o.F.length) { const b = V.length; for (const p of o.V) V.push(ap12(T, p)); for (const f of o.F) F.push([f[0] + b, f[1] + b, f[2] + b]); parts.push({ name: o.name || `object ${id}`, faces: o.F.length }); } for (const c of o.comp) emit(c.id, mul12(c.T, T), depth + 1); };
  const items = [...x.matchAll(/<item\b([^>]*)\/?>/g)];
  if (items.length) for (const it of items) emit(attr(it[1], 'objectid'), mat12(attr(it[1], 'transform'))); else { notes.push('no <build> items: all mesh objects imported untransformed'); for (const id of objs.keys()) emit(id, null); }
  if (!F.length) fail('3MF', 'no triangles found in the build items');
  if (F.some(f => f.some(i => !(i >= 0 && i < V.length)))) fail('3MF', 'triangle references an undefined vertex');
  if (/<(m:)?basematerials|<m:colorgroup|pid=/.test(x)) lost.push('materials, colours and property groups not imported'); if (/beamlattice/i.test(x)) lost.push('beam-lattice extension data not imported'); if (/slicestack/i.test(x)) lost.push('slice-extension data not imported');
  const u = ({ millimeter: 'mm', micron: 'µm', centimeter: 'cm', meter: 'm', inch: 'in', foot: 'ft' })[unit] || 'mm';
  return { kind: 'mesh', vertices: V, faces: F, units: u, meta: { unit, objects: objs.size, parts, title: (/<metadata name="Title">([^<]*)/.exec(x) || [])[1] }, notes, lost };
}
export function toAMF(mesh) {
  const { vertices: V, faces: F } = needMesh(mesh, 'AMF'), u = ({ mm: 'millimeter', um: 'micron', 'µm': 'micron', m: 'meter', in: 'inch', ft: 'feet' })[mesh.units] || 'millimeter';
  return `<?xml version="1.0" encoding="UTF-8"?>\n<amf unit="${u}" version="1.1">\n <object id="0">\n  <mesh>\n   <vertices>\n` + V.map(p => `    <vertex><coordinates><x>${N9(p[0])}</x><y>${N9(p[1])}</y><z>${N9(p[2] || 0)}</z></coordinates></vertex>`).join('\n') +
    '\n   </vertices>\n   <volume>\n' + F.map(f => `    <triangle><v1>${f[0]}</v1><v2>${f[1]}</v2><v3>${f[2]}</v3></triangle>`).join('\n') + '\n   </volume>\n  </mesh>\n </object>\n</amf>\n';
}
export function parseAMF(input) {
  let b = toBytes(input); const notes = [];
  if (b[0] === 0x50 && b[1] === 0x4B) { const f = unzip(b).find(e => e.data && /<amf\b/.test(latin(e.data, 0, 2048))); if (!f) fail('AMF', 'zip archive holds no AMF document'); b = f.data; notes.push('compressed (zip) AMF container unpacked'); }
  const x = toText(b); if (!/<amf\b/.test(x)) fail('AMF', 'no <amf> root element');
  const unit = attr((/<amf\b[^>]*>/.exec(x) || [''])[0], 'unit') || 'millimeter', V = [], F = []; let nObj = 0, nVol = 0;
  for (const o of x.split(/<object\b/).slice(1)) {
    nObj++; const base = V.length; let m;
    const rv = /<coordinates>\s*<x>([^<]+)<\/x>\s*<y>([^<]+)<\/y>\s*<z>([^<]+)<\/z>/g; while ((m = rv.exec(o))) V.push([+m[1], +m[2], +m[3]]);
    const rt = /<triangle>\s*<v1>\s*(\d+)\s*<\/v1>\s*<v2>\s*(\d+)\s*<\/v2>\s*<v3>\s*(\d+)\s*<\/v3>/g; while ((m = rt.exec(o))) F.push([+m[1] + base, +m[2] + base, +m[3] + base]);
    nVol += (o.match(/<volume\b/g) || []).length;
  }
  if (!F.length) fail('AMF', 'no triangles found'); if (F.some(f => f.some(i => !(i >= 0 && i < V.length)))) fail('AMF', 'triangle references an undefined vertex');
  const lost = []; if (/<material\b/.test(x)) lost.push('material definitions and per-volume material assignment not imported'); if (/<edge\b|<normal\b/.test(x)) lost.push('curved-triangle data (edge tangents / normals) ignored: facets treated as planar'); if (/<constellation\b/.test(x)) lost.push('constellation instances not expanded');
  return { kind: 'mesh', vertices: V, faces: F, units: ({ millimeter: 'mm', micron: 'µm', meter: 'm', inch: 'in', feet: 'ft' })[unit] || 'mm', meta: { unit, objects: nObj, volumes: nVol }, notes, lost };
}

// ---------------------------------------------------------------- glTF / GLB, VRML / X3D
const GLC = { 5120: 'int8', 5121: 'uint8', 5122: 'int16', 5123: 'uint16', 5125: 'uint32', 5126: 'float32' }, GLN = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4, MAT4: 16 };
export function parseGLTF(input) {
  const b = toBytes(input); let json, bin = null;
  if (latin(b, 0, 4) === 'glTF') { const v = view(b); let p = 12; while (p + 8 <= b.length) { const len = v.getUint32(p, true), type = v.getUint32(p + 4, true); if (type === 0x4E4F534A) json = JSON.parse(toText(b.subarray(p + 8, p + 8 + len))); else if (type === 0x004E4942) bin = b.subarray(p + 8, p + 8 + len); p += 8 + len + (len % 4 ? 4 - len % 4 : 0); } }
  else { try { json = JSON.parse(toText(b)); } catch (e) { fail('glTF', 'not valid JSON: ' + e.message); } }
  if (!json || !json.asset) fail('glTF', 'missing "asset" object');
  const lost = [], notes = ['glTF convention: metres, +Y up (kept as stored)'], bufs = (json.buffers || []).map((bf, i) => bf.uri ? (bf.uri.startsWith('data:') ? b64(bf.uri.slice(bf.uri.indexOf(',') + 1)) : null) : (i === 0 ? bin : null));
  const acc = i => { const a = json.accessors[i], bv = json.bufferViews[a.bufferView], buf = bufs[bv.buffer]; if (!buf) return null; const nc = GLN[a.type], t = GLC[a.componentType]; if (!nc || !t) fail('glTF', 'unsupported accessor type'); const sz = DT[t][0], stride = bv.byteStride || nc * sz, o = (bv.byteOffset || 0) + (a.byteOffset || 0); if (!(a.count >= 0) || o + max(0, a.count - 1) * stride + nc * sz > buf.length) fail('glTF', 'an accessor points outside its buffer'); const out = new Float64Array(a.count * nc), v = view(buf); for (let k = 0; k < a.count; k++) for (let c = 0; c < nc; c++) out[k * nc + c] = v[DT[t][1]](o + k * stride + c * sz, true); return out; };
  const V = [], F = []; let skipped = 0, missing = 0;
  const mesh = (mi, M) => { for (const pr of json.meshes[mi].primitives || []) { if ((pr.mode ?? 4) !== 4) { skipped++; continue; } const P = pr.attributes && pr.attributes.POSITION !== undefined ? acc(pr.attributes.POSITION) : null; if (!P) { missing++; continue; } const base = V.length; for (let i = 0; i < P.length; i += 3) { const x = P[i], y = P[i + 1], z = P[i + 2]; V.push(M ? [M[0] * x + M[4] * y + M[8] * z + M[12], M[1] * x + M[5] * y + M[9] * z + M[13], M[2] * x + M[6] * y + M[10] * z + M[14]] : [x, y, z]); } const I = pr.indices !== undefined ? acc(pr.indices) : Float64Array.from({ length: P.length / 3 }, (_, i) => i); if (!I) { missing++; continue; } for (let i = 0; i + 2 < I.length; i += 3) F.push([base + I[i], base + I[i + 1], base + I[i + 2]]); } };
  const mm = (A, B) => { const o = new Array(16); for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) o[c * 4 + r] = A[r] * B[c * 4] + A[4 + r] * B[c * 4 + 1] + A[8 + r] * B[c * 4 + 2] + A[12 + r] * B[c * 4 + 3]; return o; };
  const local = n => { if (n.matrix) return n.matrix; const t = n.translation || [0, 0, 0], [x, y, z, w] = n.rotation || [0, 0, 0, 1], s = n.scale || [1, 1, 1]; return [(1 - 2 * (y * y + z * z)) * s[0], 2 * (x * y + z * w) * s[0], 2 * (x * z - y * w) * s[0], 0, 2 * (x * y - z * w) * s[1], (1 - 2 * (x * x + z * z)) * s[1], 2 * (y * z + x * w) * s[1], 0, 2 * (x * z + y * w) * s[2], 2 * (y * z - x * w) * s[2], (1 - 2 * (x * x + y * y)) * s[2], 0, t[0], t[1], t[2], 1]; };
  const walk = (ni, M, d) => { const n = json.nodes[ni]; if (!n || d > 32) return; const W = M ? mm(M, local(n)) : local(n); if (n.mesh !== undefined) mesh(n.mesh, W); for (const c of n.children || []) walk(c, W, d + 1); };
  const scene = json.scenes && json.scenes[json.scene ?? 0];
  if (scene && json.nodes) for (const r of scene.nodes || []) walk(r, null, 0); else (json.meshes || []).forEach((_, i) => mesh(i, null));
  if (missing) lost.push(`${missing} primitives reference external or missing buffers (supply a .glb or embedded data URIs)`); if (skipped) lost.push(`${skipped} non-triangle primitives (points/lines/strips) skipped`);
  if (json.materials) lost.push('materials and textures not imported'); if (json.animations) lost.push('animations and skins not imported');
  if (!F.length) fail('glTF', 'no triangle geometry could be read' + (missing ? ' (buffers are external to this file)' : ''));
  return { kind: 'mesh', vertices: V, faces: F, units: 'm', meta: { generator: json.asset.generator, version: json.asset.version, meshes: (json.meshes || []).length, nodes: (json.nodes || []).length }, notes, lost };
}
export function parseVRML(input) {
  const t = toText(input), V = [], F = [], x3d = /<X3D\b/.test(t); let sets = 0, polys = 0;
  const add = (ptS, idxS) => { const p = nums(ptS), base = V.length; for (let i = 0; i + 2 < p.length; i += 3) V.push([p[i], p[i + 1], p[i + 2]]); let cur = []; for (const k of [...nums(idxS), -1]) { if (k < 0) { if (cur.length >= 3) { if (cur.length > 3) polys++; fan(cur.map(i => i + base), F); } cur = []; } else cur.push(k); } sets++; };
  if (x3d) { const re = /<IndexedFaceSet\b([^>]*)>([\s\S]*?)<\/IndexedFaceSet>/g; let m; while ((m = re.exec(t))) { const c = /<Coordinate\b([^>]*)\/?>/.exec(m[2]); if (c) add(attr(c[1], 'point') || '', attr(m[1], 'coordIndex') || ''); } }
  else { const re = /IndexedFaceSet\s*\{/g; let m; while ((m = re.exec(t))) { const seg = t.slice(m.index, t.indexOf('IndexedFaceSet', m.index + 10) > 0 ? t.indexOf('IndexedFaceSet', m.index + 10) : undefined), p = /point\s*\[([^\]]*)\]/.exec(seg), c = /coordIndex\s*\[([^\]]*)\]/.exec(seg); if (p && c) add(p[1], c[1]); } }
  if (!F.length) fail(x3d ? 'X3D' : 'VRML', 'no IndexedFaceSet with inline coordinates found (primitive nodes such as Box/Sphere and binary .x3db are not tessellated)');
  if (F.some(f => f.some(i => !(i >= 0 && i < V.length)))) fail(x3d ? 'X3D' : 'VRML', 'coordIndex references an undefined point');
  return { kind: 'mesh', vertices: V, faces: F, units: 'm', meta: { faceSets: sets, dialect: x3d ? 'X3D XML' : (/#VRML V(\S+)/.exec(t) || [])[1] || 'VRML' }, notes: ['scene units assumed to be metres (VRML/X3D default)', ...(polys ? [`${polys} polygons fan-triangulated`] : [])], lost: ['Transform nodes are not applied: each face set is placed in its local coordinates', 'appearance, primitive nodes (Box, Sphere, Cylinder), extrusions and scripts not imported'] };
}

// ---------------------------------------------------------------- STEP (ISO 10303-21)
// Faceted B-rep, AP214: one planar FACE_SURFACE with a POLY_LOOP per triangle inside a CLOSED_SHELL / FACETED_BREP.
export function toSTEP(mesh, { name = 'part', timestamp } = {}) {
  const { vertices: V, faces: F } = needMesh(mesh, 'STEP'), o = [], q = s => String(s).replace(/'/g, "''").replace(/[^\x20-\x7E]/g, '_'), R = v => { const r = +Number(v).toPrecision(12); let t = String(r === 0 ? 0 : r); if (!/[.eE]/.test(t)) t += '.'; return t.replace(/e\+?/, 'E'); };
  const pre = ({ mm: '.MILLI.', um: '.MICRO.', 'µm': '.MICRO.', nm: '.NANO.', cm: '.CENTI.', m: '$' })[mesh.units] || '.MILLI.';
  let id = 0; const E = t => { o.push(`#${++id}=${t};`); return id; };
  E("APPLICATION_CONTEXT('core data for automotive mechanical design processes')"); E("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#1)"); E("PRODUCT_CONTEXT('',#1,'mechanical')");
  E(`PRODUCT('${q(name)}','${q(name)}','',(#3))`); E("PRODUCT_DEFINITION_FORMATION('','',#4)"); E("PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design')"); E("PRODUCT_DEFINITION('design','',#5,#6)"); E("PRODUCT_DEFINITION_SHAPE('','',#7)");
  E(`(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(${pre},.METRE.))`); E('(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))'); E('(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())');
  E("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#9,'distance_accuracy_value','')");
  E("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#12))GLOBAL_UNIT_ASSIGNED_CONTEXT((#9,#10,#11))REPRESENTATION_CONTEXT('',''))");
  E("CARTESIAN_POINT('',(0.,0.,0.))"); E("DIRECTION('',(0.,0.,1.))"); E("DIRECTION('',(1.,0.,0.))"); E("AXIS2_PLACEMENT_3D('',#14,#15,#16)");
  const pid = V.map(p => E(`CARTESIAN_POINT('',(${R(p[0])},${R(p[1])},${R(p[2] || 0)}))`)), fids = [];
  for (const f of F) {
    const A = V[f[0]], B = V[f[1]], n = faceNormal(A, B, V[f[2]]), ul = hypot(B[0] - A[0], B[1] - A[1], (B[2] || 0) - (A[2] || 0)); if (!(ul > 0) || !(hypot(...n) > 0.5)) continue;
    const dn = E(`DIRECTION('',(${R(n[0])},${R(n[1])},${R(n[2])}))`), dr = E(`DIRECTION('',(${R((B[0] - A[0]) / ul)},${R((B[1] - A[1]) / ul)},${R(((B[2] || 0) - (A[2] || 0)) / ul)}))`), ax = E(`AXIS2_PLACEMENT_3D('',#${pid[f[0]]},#${dn},#${dr})`), pl = E(`PLANE('',#${ax})`);
    const lp = E(`POLY_LOOP('',(#${pid[f[0]]},#${pid[f[1]]},#${pid[f[2]]}))`), bd = E(`FACE_OUTER_BOUND('',#${lp},.T.)`); fids.push(E(`FACE_SURFACE('',(#${bd}),#${pl},.T.)`));
  }
  const sh = E(`CLOSED_SHELL('',(${fids.map(i => '#' + i).join(',')}))`), br = E(`FACETED_BREP('${q(name)}',#${sh})`), rep = E(`FACETED_BREP_SHAPE_REPRESENTATION('',(#17,#${br}),#13)`); E(`SHAPE_DEFINITION_REPRESENTATION(#8,#${rep})`);
  const ts = timestamp || new Date().toISOString().slice(0, 19);
  return `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('faceted boundary representation'),'2;1');\nFILE_NAME('${q(name)}.step','${ts}',(''),(''),'Designer Materials','Designer Materials','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n${o.join('\n')}\nENDSEC;\nEND-ISO-10303-21;\n`;
}
// Reads what a Part 21 file states explicitly: header, products, units, entity inventory, and geometry that needs no
// surface kernel — faceted B-rep (POLY_LOOP), AP242 tessellation (TRIANGULATED_FACE_SET) and planar ADVANCED_FACEs
// bounded only by straight edges without inner loops. Curved or trimmed faces are counted, not tessellated; the
// topological vertices and edges are then returned as a straight-chord wireframe.
export function parseSTEP(input) {
  const t = toText(input); if (!/ISO-10303-21/.test(t.slice(0, 2000))) fail('STEP', 'missing "ISO-10303-21;" header');
  const di = t.search(/\bDATA\s*;/), head = di > 0 ? t.slice(0, di) : t.slice(0, 5000), strs = s => [...s.matchAll(/'((?:[^']|'')*)'/g)].map(m => m[1]);
  const meta = { schema: strs((/FILE_SCHEMA\s*\(([\s\S]*?)\)\s*;/.exec(head) || ['', ''])[1])[0] || '', description: strs((/FILE_DESCRIPTION\s*\(([\s\S]*?)\)\s*;/.exec(head) || ['', ''])[1])[0] || '' };
  const fn = strs((/FILE_NAME\s*\(([\s\S]*?)\)\s*;/.exec(head) || ['', ''])[1]); Object.assign(meta, { fileName: fn[0], timestamp: fn[1], originatingSystem: fn[fn.length - 2], preprocessor: fn[fn.length - 3] });
  meta.protocol = /IFC/i.test(meta.schema) ? 'IFC' : /AP242|MANAGED_MODEL_BASED/i.test(meta.schema) ? 'AP242' : /AUTOMOTIVE_DESIGN/i.test(meta.schema) ? 'AP214' : /CONFIG_CONTROL_DESIGN/i.test(meta.schema) ? 'AP203' : meta.schema.split(/[\s{]/)[0];
  const counts = {}, keep = new Map(), KEEP = /^(CARTESIAN_POINT|VERTEX_POINT|EDGE_CURVE|ORIENTED_EDGE|EDGE_LOOP|POLY_LOOP|FACE_OUTER_BOUND|FACE_BOUND|ADVANCED_FACE|FACE_SURFACE|PLANE|LINE|COORDINATES_LIST|TRIANGULATED_FACE_SET|PRODUCT)$/;
  // single linear pass over "#id = body ;" instances (semicolons inside quoted strings do not terminate an instance)
  let lenUnit = null, nEnt = 0, pos = max(0, di); const tl = t.length, m = [0, 0, ''];
  for (; ;) {
    const hsh = t.indexOf('#', pos); if (hsh < 0) break; let q = hsh + 1, idn = 0, nd = 0; for (; q < tl; q++) { const c = t.charCodeAt(q); if (c < 48 || c > 57) break; idn = idn * 10 + c - 48; nd++; }
    while (q < tl && t.charCodeAt(q) <= 32) q++; if (!nd || t[q] !== '=') { pos = hsh + 1; continue; }
    let e = q + 1, inq = false; for (; e < tl; e++) { const c = t.charCodeAt(e); if (c === 39) inq = !inq; else if (c === 59 && !inq) break; } if (e >= tl) break;
    m[1] = idn; m[2] = t.slice(q + 1, e); pos = e + 1;
    nEnt++; const body = m[2].trim();
    if (body[0] === '(') { for (const k of body.matchAll(/([A-Z_0-9]+)\s*\(/g)) counts[k[1]] = (counts[k[1]] || 0) + 1; if (/LENGTH_UNIT/.test(body)) { const si = /SI_UNIT\s*\(\s*(\$|\.[A-Z]+\.)\s*,\s*\.METRE\./.exec(body); if (si && !lenUnit) lenUnit = ({ $: 'm', '.MILLI.': 'mm', '.MICRO.': 'µm', '.NANO.': 'nm', '.CENTI.': 'cm' })[si[1]] || null; const cb = /CONVERSION_BASED_UNIT\s*\(\s*'([^']*)'/.exec(body); if (cb && !lenUnit) lenUnit = /INCH/i.test(cb[1]) ? 'in' : /FOOT|FEET/i.test(cb[1]) ? 'ft' : cb[1].toLowerCase(); } continue; }
    const p = body.indexOf('('), nm = (p < 0 ? body : body.slice(0, p)).trim(); counts[nm] = (counts[nm] || 0) + 1;
    if (KEEP.test(nm)) keep.set(+m[1], [nm, body.slice(p + 1, body.lastIndexOf(')'))]);
  }
  if (!nEnt) fail('STEP', 'DATA section holds no entity instances');
  const refs = s => [...s.matchAll(/#(\d+)/g)].map(x => +x[1]), pt = new Map(), V = [], F = [], vix = new Map(), products = [], notes = [], lost = [];
  for (const [id, [nm, a]] of keep) { if (nm === 'CARTESIAN_POINT') { const c = nums((/\(([^()]*)\)\s*$/.exec(a) || ['', ''])[1]); pt.set(id, [c[0] || 0, c[1] || 0, c[2] || 0]); } else if (nm === 'PRODUCT') { const s = strs(a); products.push(s[1] || s[0] || ''); } }
  const vert = pidx => { let i = vix.get(pidx); if (i === undefined) { const p = pt.get(pidx); if (!p) return -1; i = V.length; V.push(p); vix.set(pidx, i); } return i; };
  const vpoint = id => { const k = keep.get(id); return k && k[0] === 'VERTEX_POINT' ? refs(k[1])[0] : undefined; };
  let facet = 0, tess = 0, planar = 0, curved = 0, holed = 0;
  const poly3 = ids => { // triangulate a planar polygon given as vertex indices
    if (ids.length < 3 || ids.includes(-1)) return false; if (ids.length === 3) { F.push(ids); return true; }
    const n = [0, 0, 0]; for (let i = 0; i < ids.length; i++) { const a = V[ids[i]], b = V[ids[(i + 1) % ids.length]]; n[0] += (a[1] - b[1]) * (a[2] + b[2]); n[1] += (a[2] - b[2]) * (a[0] + b[0]); n[2] += (a[0] - b[0]) * (a[1] + b[1]); }
    const ax = abs(n[0]) > abs(n[1]) && abs(n[0]) > abs(n[2]) ? 0 : abs(n[1]) > abs(n[2]) ? 1 : 2, u = (ax + 1) % 3, w = (ax + 2) % 3, sgn = n[ax] < 0;
    const tr = triangulate(ids.map(i => sgn ? [V[i][w], V[i][u]] : [V[i][u], V[i][w]])); for (const q of tr) F.push([ids[q[0]], ids[q[1]], ids[q[2]]]); return tr.length > 0;
  };
  for (const [, [nm, a]] of keep) {
    if (nm === 'POLY_LOOP') { if (poly3(refs(a).map(vert))) facet++; }
    else if (nm === 'TRIANGULATED_FACE_SET') {
      const r = refs(a), cl = keep.get(r[0]); if (!cl || cl[0] !== 'COORDINATES_LIST') continue;
      const co = nums(cl[1].slice(cl[1].indexOf('((')).replace(/[()]/g, ' ')), base = V.length; for (let i = 0; i + 2 < co.length; i += 3) V.push([co[i], co[i + 1], co[i + 2]]);
      const grp = [...a.matchAll(/\(\s*\(([\d\s,()]*?)\)\s*\)|\(([\d\s,]*)\)/g)], tri = grp.length ? grp[grp.length - 1] : null, pn = /\$\s*,\s*\(([\d\s,]*)\)\s*,\s*\(\s*\(/.exec(a), map = pn && pn[1].trim() ? nums(pn[1]) : null;
      const ti = tri ? nums((tri[1] ?? tri[2]).replace(/[()]/g, ' ')) : []; for (let i = 0; i + 2 < ti.length; i += 3) { const f = [ti[i], ti[i + 1], ti[i + 2]].map(k => base + (map ? map[k - 1] : k) - 1); if (f.every(k => k >= base && k < V.length)) { F.push(f); tess++; } }
    }
  }
  const lineIds = new Set([...keep].filter(([, v]) => v[0] === 'LINE').map(([k]) => k)), wire = new Map(), wn = [], we = [];
  const wv = pid => { let i = wire.get(pid); if (i === undefined) { const p = pt.get(pid); if (!p) return -1; i = wn.length; wn.push(p); wire.set(pid, i); } return i; };
  for (const [, [nm, a]] of keep) if (nm === 'EDGE_CURVE') { const r = refs(a), p1 = vpoint(r[0]), p2 = vpoint(r[1]); if (p1 !== undefined && p2 !== undefined && p1 !== p2) { const i = wv(p1), j = wv(p2); if (i >= 0 && j >= 0) we.push([i, j]); } }
  for (const [, [nm, a]] of keep) {
    if (nm !== 'ADVANCED_FACE' && nm !== 'FACE_SURFACE') continue;
    const r = refs(a), surf = keep.get(r[r.length - 1]), bounds = r.slice(0, -1).map(b => keep.get(b)).filter(Boolean);
    if (bounds.some(b => keep.get(refs(b[1])[0])?.[0] === 'POLY_LOOP')) continue;
    if (!surf || surf[0] !== 'PLANE') { curved++; continue; }
    if (bounds.length !== 1) { holed++; continue; }
    const lp = keep.get(refs(bounds[0][1])[0]); if (!lp || lp[0] !== 'EDGE_LOOP') { curved++; continue; }
    const ids = []; let ok = true;
    for (const oe of refs(lp[1])) { const k = keep.get(oe); if (!k || k[0] !== 'ORIENTED_EDGE') { ok = false; break; } const er = refs(k[1]), ec = keep.get(er[er.length - 1]); if (!ec || ec[0] !== 'EDGE_CURVE') { ok = false; break; } const cr = refs(ec[1]); if (!lineIds.has(cr[2])) { ok = false; break; } const fwd = !/\.F\.\s*$/.test(k[1]); ids.push(vert(vpoint(fwd ? cr[0] : cr[1]))); }
    if (!ok) { curved++; continue; }
    if (/\.F\.\s*$/.test(bounds[0][1])) ids.reverse();
    if (poly3(ids)) planar++; else curved++;
  }
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (const p of pt.values()) for (let k = 0; k < 3; k++) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
  const top = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  Object.assign(meta, { products: [...new Set(products)].slice(0, 50), entities: nEnt, entityCounts: Object.fromEntries(top.slice(0, 40)), solids: (counts.MANIFOLD_SOLID_BREP || 0) + (counts.FACETED_BREP || 0) + (counts.BREP_WITH_VOIDS || 0), brepFaces: (counts.ADVANCED_FACE || 0) + (counts.FACE_SURFACE || 0), pointBounds: pt.size ? [lo, hi] : null,
    surfaces: Object.fromEntries(top.filter(([k]) => /SURFACE$|^PLANE$|SURFACE_WITH_KNOTS/.test(k) && !/FACE_SURFACE|CURVE/.test(k))), tessellated: { polyLoops: facet, triangulatedFaceSetTriangles: tess, planarFaces: planar, untessellatedCurvedFaces: curved, planarFacesWithHolesSkipped: holed } });
  const units = lenUnit || 'mm'; if (!lenUnit) notes.push('no length unit found in the file: millimetres assumed');
  if (counts.NEXT_ASSEMBLY_USAGE_OCCURRENCE) lost.push(`assembly structure (${counts.NEXT_ASSEMBLY_USAGE_OCCURRENCE} component placements) not applied: every body is shown in its own definition coordinates`);
  if (counts.COLOUR_RGB || counts.STYLED_ITEM) lost.push('colours, layers and presentation styles not imported');
  if (counts.COMPLEX_TRIANGULATED_FACE_SET) lost.push(`${counts.COMPLEX_TRIANGULATED_FACE_SET} complex (strip/fan) tessellated face sets not decoded`);
  if (curved || holed) lost.push(`${curved} curved or curve-bounded faces and ${holed} planar faces with inner loops were not tessellated: exact B-rep surfaces need a CAD kernel. Re-export as AP242 with tessellation, or as STL/3MF, for the full shape`);
  if (F.length) { notes.push(`surface built from ${facet} faceted loops, ${tess} tessellated triangles and ${planar} planar straight-edged B-rep faces`); notes.push('face orientation flags are not evaluated: winding is made consistent by the mesh repair step'); return { kind: 'mesh', vertices: V, faces: F, units, complete: curved + holed === 0, meta, notes, lost, ...(we.length && (curved || holed) ? { wireframe: { nodes: wn, edges: we } } : {}) }; }
  if (we.length) { notes.push('no tessellatable faces: the B-rep vertices and edges are shown as a wireframe with curved edges drawn as straight chords'); return { kind: 'lattice', nodes: wn, edges: we, radius: 0, units, wireframe: true, meta, notes, lost }; }
  notes.push('no boundary topology readable: only the header, product and entity inventory were extracted');
  return { kind: null, units, meta, notes, lost: [...lost, 'geometry not reconstructed'] };
}

// ---------------------------------------------------------------- IGES
const IGES_NAMES = { 100: 'circular arc', 102: 'composite curve', 104: 'conic arc', 106: 'copious data / polyline', 108: 'plane', 110: 'line', 112: 'parametric spline curve', 114: 'parametric spline surface', 116: 'point', 118: 'ruled surface', 120: 'surface of revolution', 122: 'tabulated cylinder', 124: 'transformation matrix', 126: 'rational B-spline curve', 128: 'rational B-spline surface', 130: 'offset curve', 140: 'offset surface', 141: 'boundary', 142: 'curve on parametric surface', 143: 'bounded surface', 144: 'trimmed parametric surface', 186: 'manifold solid B-rep', 190: 'plane surface', 192: 'right circular cylindrical surface', 194: 'right circular conical surface', 196: 'spherical surface', 198: 'toroidal surface', 308: 'subfigure definition', 314: 'colour definition', 402: 'associativity instance', 406: 'property', 408: 'singular subfigure instance', 502: 'vertex list', 504: 'edge list', 508: 'loop', 510: 'face', 514: 'shell' };
export function parseIGES(input) {
  const L = lines(toText(input)).filter(l => l.length >= 73), sec = { S: [], G: [], D: [], P: [], T: [] };
  for (const l of L) { const s = sec[l[72]]; if (s) s.push(l); }
  if (!sec.D.length || !sec.G.length) fail('IGES', 'no fixed-format Global/Directory sections found (columns 73 must hold S, G, D, P, T); compressed-ASCII and binary IGES are not supported');
  const g = sec.G.map(l => l.slice(0, 72)).join(''); let pd = ',', rd = ';', p = 0; const gf = [];
  const tok = () => { while (g[p] === ' ') p++; const h = /^(\d+)H/.exec(g.slice(p, p + 12)); let v; if (h) { v = g.substr(p + h[0].length, +h[1]); p += h[0].length + +h[1]; } else { let e = p; while (e < g.length && g[e] !== pd && g[e] !== rd) e++; v = g.slice(p, e).trim(); p = e; } const end = g[p] === rd; p++; return [v, end]; };
  if (g.startsWith('1H')) { pd = g[2]; p = 4; gf.push(pd); if (g.startsWith('1H', p)) { rd = g[p + 2]; p += 4; gf.push(rd); } else { p++; gf.push(';'); } } else { p = g.indexOf(',') + 1; gf.push(','); if (g[p] === ',') { p++; gf.push(';'); } }
  while (p < g.length && gf.length < 30) { const [v, end] = tok(); gf.push(v); if (end) break; }
  const UF = { 1: 'in', 2: 'mm', 4: 'ft', 5: 'mi', 6: 'm', 7: 'km', 8: 'mil', 9: 'µm', 10: 'cm', 11: 'µin' }, flag = +gf[13] || 1, units = UF[flag] || (gf[14] || '').toLowerCase() || 'in';
  const meta = { productId: gf[2], fileName: gf[3], system: gf[4], preprocessor: gf[5], modelScale: +gf[12] || 1, unitsFlag: flag, unitsName: gf[14], created: gf[17], author: gf[20], organisation: gf[21], igesVersionFlag: +gf[22] || null, startSection: sec.S.map(l => l.slice(0, 72).trim()).join(' ').slice(0, 300) };
  const counts = {}, ents = []; let xform = 0;
  for (let i = 0; i + 1 < sec.D.length; i += 2) { const a = sec.D[i], b = sec.D[i + 1], type = parseInt(a.slice(0, 8), 10); if (!type) continue; counts[type] = (counts[type] || 0) + 1; const tm = parseInt(a.slice(48, 56), 10) || 0; if ([100, 106, 110, 116].includes(type)) { ents.push({ type, ptr: parseInt(a.slice(8, 16), 10), n: parseInt(b.slice(24, 32), 10) || 1, form: parseInt(b.slice(32, 40), 10) || 0 }); if (tm) xform++; } }
  const nodes = [], edges = [], np = [], num1 = s => Number(String(s).trim().replace(/[dD]/, 'e'));
  const rec = e => nums(sec.P.slice(e.ptr - 1, e.ptr - 1 + e.n).map(l => l.slice(0, 64)).join('').split(rd)[0].replace(/[dD]([+-]?\d)/g, 'e$1'));
  const add = q => { nodes.push(q); return nodes.length - 1; };
  for (const e of ents) {
    if (e.ptr < 1 || e.ptr > sec.P.length) continue; const a = rec(e);
    if (e.type === 116) np.push(add([a[1], a[2], a[3] || 0]));
    else if (e.type === 110) edges.push([add([a[1], a[2], a[3]]), add([a[4], a[5], a[6]])]);
    else if (e.type === 100) { const [, zt, cx, cy, x1, y1, x2, y2] = a, r = hypot(x1 - cx, y1 - cy), a0 = Math.atan2(y1 - cy, x1 - cx); let a1 = Math.atan2(y2 - cy, x2 - cx); if (a1 <= a0 + 1e-12) a1 += 2 * PI; const n = max(2, Math.ceil((a1 - a0) / (PI / 12))); let prev = add([x1, y1, zt]); for (let k = 1; k <= n; k++) { const t = a0 + (a1 - a0) * k / n, c = add([cx + r * cos(t), cy + r * sin(t), zt]); edges.push([prev, c]); prev = c; } }
    else if (e.type === 106 && [1, 2, 11, 12, 63].includes(e.form)) { const ip = a[1], n = a[2]; let prev = -1; for (let k = 0; k < n; k++) { const q = ip === 1 ? [a[4 + 2 * k], a[5 + 2 * k], a[3]] : [a[3 + 3 * k], a[4 + 3 * k], a[5 + 3 * k]]; if (!Number.isFinite(q[0] + q[1] + q[2])) break; const c = add(q); if (prev >= 0) edges.push([prev, c]); prev = c; } }
  }
  void num1;
  meta.entityCounts = Object.fromEntries(Object.entries(counts).map(([k, v]) => [`${k} ${IGES_NAMES[k] || 'entity'}`, v])); meta.entities = Object.values(counts).reduce((a, b) => a + b, 0);
  const surf = [108, 114, 118, 120, 122, 128, 140, 143, 144, 186, 190, 192, 194, 196, 198, 510, 514].reduce((s, k) => s + (counts[k] || 0), 0), crv = [102, 104, 112, 126, 130].reduce((s, k) => s + (counts[k] || 0), 0);
  const notes = [`units flag ${flag} → ${units}`], lost = [];
  if (surf) lost.push(`${surf} surface / solid entities (NURBS, trimmed, analytic) not tessellated: needs a CAD kernel — re-export as STEP AP242 with tessellation, or STL`);
  if (crv) lost.push(`${crv} spline / conic / composite curve entities not sampled`); if (xform) lost.push(`${xform} wireframe entities reference a transformation matrix (type 124) that was not applied`);
  if (!nodes.length) return { kind: null, units, meta, notes: [...notes, 'entity inventory only: no points, lines, arcs or polylines present'], lost: [...lost, 'geometry not reconstructed'] };
  notes.push(`wireframe from ${counts[110] || 0} lines, ${counts[100] || 0} arcs, ${counts[106] || 0} polylines and ${counts[116] || 0} points`);
  return { kind: 'lattice', nodes, edges, radius: 0, units, wireframe: true, pointCloud: !edges.length, meta, notes, lost };
}

// ---------------------------------------------------------------- VTK legacy
const VTKC = { 3: 'line', 5: 'tri', 9: 'quad', 10: 'tet', 12: 'hex', 13: 'wedge', 14: 'pyr', 21: 'line', 22: 'tri', 23: 'quad', 24: 'tet', 25: 'hex', 26: 'wedge', 27: 'pyr', 28: 'quad', 29: 'hex' };
function vtkCell(type, ids, elems, skip) {
  if (type === 8) elems.push({ type: 'quad', n: [ids[0], ids[1], ids[3], ids[2]] }); else if (type === 11) elems.push({ type: 'hex', n: [ids[0], ids[1], ids[3], ids[2], ids[4], ids[5], ids[7], ids[6]] });
  else if (type === 4) for (let i = 0; i + 1 < ids.length; i++) elems.push({ type: 'line', n: [ids[i], ids[i + 1]] });
  else if (type === 6) for (let i = 0; i + 2 < ids.length; i++) elems.push({ type: 'tri', n: i % 2 ? [ids[i + 1], ids[i], ids[i + 2]] : [ids[i], ids[i + 1], ids[i + 2]] });
  else if (type === 7) for (let i = 1; i + 1 < ids.length; i++) elems.push({ type: 'tri', n: [ids[0], ids[i], ids[i + 1]] });
  else if (type === 13 || type === 26) elems.push({ type: 'wedge', n: [ids[0], ids[2], ids[1], ids[3], ids[5], ids[4]] }); // VTK wedge winding is opposite to the Abaqus/Gmsh convention used here
  else if (VTKC[type]) elems.push({ type: VTKC[type], n: ids }); else if (type !== 1 && type !== 2) skip[type] = (skip[type] || 0) + 1;
}
export function toVTK(g) {
  const h = `# vtk DataFile Version 3.0\nDesigner Materials ${g.kind} (units: ${g.units || 'unspecified'})\nASCII\n`;
  if (g.kind === 'voxels') { const N = g.nx * g.ny * g.nz, o = g.origin || [0, 0, 0], sp = g.spacing || [g.voxel || 1, g.voxel || 1, g.voxel || 1], rows = []; for (let i = 0; i < N; i += 40) rows.push(Array.prototype.join.call(g.data.subarray(i, min(N, i + 40)), ' ')); return h + `DATASET STRUCTURED_POINTS\nDIMENSIONS ${g.nx} ${g.ny} ${g.nz}\nORIGIN ${o.map(N9).join(' ')}\nSPACING ${sp.map(N9).join(' ')}\nPOINT_DATA ${N}\nSCALARS phase unsigned_char 1\nLOOKUP_TABLE default\n${rows.join('\n')}\n`; }
  if (g.kind === 'lattice') return h + `DATASET POLYDATA\nPOINTS ${g.nodes.length} double\n${g.nodes.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n')}\nLINES ${g.edges.length} ${3 * g.edges.length}\n${g.edges.map(e => `2 ${e[0]} ${e[1]}`).join('\n')}\n`;
  needMesh(g, 'VTK'); return h + `DATASET POLYDATA\nPOINTS ${g.vertices.length} double\n${g.vertices.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n')}\nPOLYGONS ${g.faces.length} ${4 * g.faces.length}\n${g.faces.map(f => `3 ${f[0]} ${f[1]} ${f[2]}`).join('\n')}\n`;
}
export function parseVTK(input) {
  const b = toBytes(input), headTxt = latin(b, 0, 600), hl = lines(headTxt);
  if (!/^# vtk DataFile Version/.test(hl[0] || '')) fail('VTK', 'missing "# vtk DataFile Version" header');
  if (/^BINARY/i.test((hl[2] || '').trim())) fail('VTK', 'legacy BINARY encoding is not supported by this reader: re-save as legacy ASCII or as XML (.vtu/.vtp/.vti)');
  const text = toText(b), nl2 = text.indexOf('\n', text.indexOf('\n') + 1), tk = text.slice(nl2 + 1).split(/\s+/).filter(Boolean); let p = 0;
  const meta = { version: (hl[0].match(/Version\s+(\S+)/) || [])[1], title: (hl[1] || '').trim() }, pts = [], elems = [], skip = {}, fields = []; let dataset = '', dims = null, spacing = [1, 1, 1], origin = [0, 0, 0], scal = null, scalAt = '', where = '', nData = 0;
  const take = n => { if (!(n >= 0) || n > tk.length - p) fail('VTK', `a section declares ${n} values but only ${tk.length - p} remain in the file (truncated or malformed)`); const a = new Float64Array(n); for (let i = 0; i < n; i++) a[i] = +tk[p++]; return a; };
  const conn = (handler) => { const n = +tk[p++], size = +tk[p++]; if (tk[p] === 'OFFSETS') { p += 2; const off = take(n); if (tk[p] !== 'CONNECTIVITY') fail('VTK', 'OFFSETS without CONNECTIVITY'); p += 2; const c = take(size); const cells = []; for (let i = 0; i + 1 < n; i++) cells.push(Array.from(c.subarray(off[i], off[i + 1]))); handler(cells); return cells.length; } const cells = []; for (let i = 0; i < n; i++) { const k = +tk[p++]; cells.push(Array.from(take(k))); } handler(cells); return n; };
  let cellList = null, types = null;
  while (p < tk.length) {
    const k = tk[p++].toUpperCase();
    if (k === 'ASCII') continue;
    if (k === 'DATASET') dataset = tk[p++].toUpperCase();
    else if (k === 'POINTS') { const n = +tk[p++]; p++; const a = take(3 * n); for (let i = 0; i < n; i++) pts.push([a[3 * i], a[3 * i + 1], a[3 * i + 2]]); }
    else if (k === 'POLYGONS') conn(c => c.forEach(ids => vtkCell(ids.length === 3 ? 5 : ids.length === 4 ? 9 : 7, ids, elems, skip)));
    else if (k === 'TRIANGLE_STRIPS') conn(c => c.forEach(ids => vtkCell(6, ids, elems, skip)));
    else if (k === 'LINES') conn(c => c.forEach(ids => vtkCell(4, ids, elems, skip)));
    else if (k === 'VERTICES') conn(() => { });
    else if (k === 'CELLS') conn(c => { cellList = c; });
    else if (k === 'CELL_TYPES') { const n = +tk[p++]; types = take(n); }
    else if (k === 'DIMENSIONS') dims = Array.from(take(3));
    else if (k === 'SPACING' || k === 'ASPECT_RATIO') spacing = Array.from(take(3));
    else if (k === 'ORIGIN') origin = Array.from(take(3));
    else if (k === 'POINT_DATA' || k === 'CELL_DATA') { where = k; nData = +tk[p++]; }
    else if (k === 'SCALARS') { const name = tk[p++]; p++; let nc = 1; if (/^\d+$/.test(tk[p])) nc = +tk[p++]; if ((tk[p] || '').toUpperCase() === 'LOOKUP_TABLE') p += 2; const a = take(nData * nc); fields.push(`${where} ${name}`); if (!scal && nc === 1) { scal = a; scalAt = where; meta.scalar = name; } }
    else if (k === 'VECTORS' || k === 'NORMALS') { fields.push(`${where} ${tk[p]}`); p += 2; p += 3 * nData; }
    else if (k === 'TENSORS') { fields.push(`${where} ${tk[p]}`); p += 2; p += 9 * nData; }
    else if (k === 'TEXTURE_COORDINATES') { const d = +tk[p + 1]; p += 3; p += d * nData; }
    else if (k === 'COLOR_SCALARS') { const d = +tk[p + 1]; p += 2; p += d * nData; }
    else if (k === 'LOOKUP_TABLE') { const n = +tk[p + 1]; p += 2 + 4 * n; }
    else if (k === 'FIELD') { p++; const n = +tk[p++]; for (let i = 0; i < n; i++) { const nm = tk[p++], nc = +tk[p++], nt = +tk[p++]; p++; p += nc * nt; fields.push(`field ${nm}`); } }
    else if (k === 'METADATA' || k === 'INFORMATION') { /* v4+ metadata blocks: tokens are skipped one by one */ }
  }
  meta.dataset = dataset; meta.fields = fields;
  const lost = fields.length ? [`data arrays not attached to the geometry: ${fields.slice(0, 12).join(', ')}`] : [];
  if (dataset === 'STRUCTURED_POINTS') {
    if (!dims) fail('VTK', 'STRUCTURED_POINTS without DIMENSIONS'); if (!scal) fail('VTK', 'STRUCTURED_POINTS without a scalar array to interpret as phases');
    const d = scalAt === 'CELL_DATA' ? dims.map(v => max(1, v - 1)) : dims, vx = scalarToVoxels(scal, d[0], d[1], d[2], { voxel: spacing[0], spacing, origin, units: '' });
    return Object.assign(vx, { meta, lost: fields.length > 1 ? [`only the first scalar array ("${meta.scalar}") was imported; others: ${fields.slice(1).join(', ')}`] : [], notes: [...vx.notes, 'VTK stores no length unit: spacing taken as stored'] });
  }
  if (dataset !== 'POLYDATA' && dataset !== 'UNSTRUCTURED_GRID') fail('VTK', `dataset type ${dataset || '(none)'} is not supported (supported: POLYDATA, UNSTRUCTURED_GRID, STRUCTURED_POINTS)`);
  if (cellList) { if (!types || types.length < cellList.length) fail('VTK', 'CELLS without matching CELL_TYPES'); cellList.forEach((ids, i) => vtkCell(types[i], ids, elems, skip)); }
  if (!pts.length) fail('VTK', 'no POINTS section');
  const g = femToGeometry(pts, elems, { meta }); g.notes.push('VTK stores no length unit: coordinates taken as stored'); g.lost.push(...lost);
  if (Object.keys(skip).length) g.lost.push(`unsupported VTK cell types skipped: ${Object.entries(skip).map(([k, v]) => `type ${k} × ${v}`).join(', ')}`);
  if (scal && scalAt === 'CELL_DATA' && scal.length === elems.length) g.cellData = { name: meta.scalar, values: scal };
  return g;
}

// ---------------------------------------------------------------- VTK XML (.vtu / .vtp / .vti)
function vtkXml(input, type) {
  const b = toBytes(input), t = toText(b), root = /<VTKFile\b([^>]*)>/.exec(t);
  if (!root) fail('VTK XML', 'no <VTKFile> root element'); const ft = attr(root[1], 'type'); if (type && ft !== type) fail('VTK XML', `file type is ${ft}, expected ${type}`);
  const le = (attr(root[1], 'byte_order') || 'LittleEndian') === 'LittleEndian', ht = (attr(root[1], 'header_type') || 'UInt32').toLowerCase(), hs = DT[ht][0], comp = attr(root[1], 'compressor') || '';
  if (comp && !/ZLib/i.test(comp)) fail('VTK XML', `compressor ${comp} is not supported (only zlib); re-save uncompressed or ASCII`);
  const ap = /<AppendedData\b([^>]*)>\s*_/.exec(t); let rawStart = -1;
  if (ap && attr(ap[1], 'encoding') === 'raw') { const key = ENC.encode('<AppendedData'); outer: for (let i = 0; i < b.length - key.length; i++) { for (let k = 0; k < key.length; k++) if (b[i + k] !== key[k]) continue outer; let j = i; while (b[j] !== 0x3E) j++; while (b[j] !== 0x5F) j++; rawStart = j + 1; break; } }
  const apText = ap && rawStart < 0 ? t.slice(ap.index + ap[0].length) : null;
  const hdr = (bytes, n) => readTyped(bytes, 0, n, ht, le), c4 = n => Math.ceil(n / 3) * 4;
  const inflateBlocks = (h0, cs, at) => { const parts = []; let off = 0, total = 0; for (let i = 0; i < cs.length; i++) { const u = unzlib(at(off, cs[i]), i === cs.length - 1 && h0[2] ? h0[2] : h0[1]); parts.push(u); total += u.length; off += cs[i]; } const out = new Uint8Array(total); let o = 0; for (const q of parts) { out.set(q, o); o += q.length; } return out; };
  const fromRaw = base => { // [header][data] as raw bytes
    if (!comp) { const n = hdr(b.subarray(base, base + hs), 1)[0]; return b.subarray(base + hs, base + hs + n); }
    const h0 = hdr(b.subarray(base, base + 3 * hs), 3), nb = h0[0], d0 = base + (3 + nb) * hs, cs = hdr(b.subarray(base + 3 * hs, d0), nb);
    return inflateBlocks(h0, cs, (o, n) => b.subarray(d0 + o, d0 + o + n));
  };
  const fromB64 = s => { // uncompressed: one base64 stream holding [size][data] (or two padded streams); compressed: header stream, then data stream
    if (!comp) { const sep = s.slice(0, c4(hs)).includes('='), n = hdr(b64(s.slice(0, c4(hs))), 1)[0]; return sep ? b64(s.slice(c4(hs), c4(hs) + c4(n))) : b64(s.slice(0, c4(hs + n))).subarray(hs, hs + n); }
    const h0 = hdr(b64(s.slice(0, c4(3 * hs))), 3), nb = h0[0], hb = (3 + nb) * hs, cs = hdr(b64(s.slice(0, c4(hb))), 3 + nb).subarray(3); let tot = 0; for (const v of cs) tot += v;
    const data = b64(s.slice(c4(hb), c4(hb) + c4(tot))); return inflateBlocks(h0, cs, (o, n) => data.subarray(o, o + n));
  };
  const arr = (attrs, inner) => {
    const fmt = attr(attrs, 'format') || 'ascii', ty = (attr(attrs, 'type') || 'Float32').toLowerCase();
    if (fmt === 'ascii') return Float64Array.from(nums(inner));
    if (!DT[ty]) fail('VTK XML', `unsupported DataArray type ${ty}`);
    let bytes;
    if (fmt === 'binary') bytes = fromB64(inner.replace(/\s+/g, ''));
    else if (fmt === 'appended') { const off = +attr(attrs, 'offset') || 0; if (rawStart >= 0) bytes = fromRaw(rawStart + off); else if (apText) bytes = fromB64(apText.slice(off)); else fail('VTK XML', 'appended DataArray without <AppendedData>'); }
    else fail('VTK XML', `unknown DataArray format "${fmt}"`);
    return readTyped(bytes, 0, floor(bytes.length / DT[ty][0]), ty, le);
  };
  const section = (name, scope = t) => { const m = new RegExp('<' + name + '\\b([^>]*?)(?:/>|>([\\s\\S]*?)</' + name + '>)').exec(scope); return m ? { attrs: m[1], body: m[2] || '' } : null; };
  const arrays = body => [...body.matchAll(/<DataArray\b([^>]*?)(?:\/>|>([\s\S]*?)<\/DataArray>)/g)].map(m => ({ name: attr(m[1], 'Name'), nc: +attr(m[1], 'NumberOfComponents') || 1, get: () => arr(m[1], m[2] || '') }));
  return { t, ft, section, arrays, version: attr(root[1], 'version'), encoding: comp ? 'zlib-compressed binary' : ap ? `appended ${rawStart >= 0 ? 'raw' : 'base64'}` : /format="binary"/.test(t) ? 'inline base64' : 'ASCII' };
}
function vtkXmlMesh(input, type) {
  const x = vtkXml(input, type), piece = x.section('Piece'); if (!piece) fail('VTK XML', 'no <Piece>');
  const P = x.section('Points', piece.body); if (!P) fail('VTK XML', 'no <Points>'); const pa = x.arrays(P.body)[0].get(), pts = []; for (let i = 0; i + 2 < pa.length; i += 3) pts.push([pa[i], pa[i + 1], pa[i + 2]]);
  const elems = [], skip = {}, named = (sec, n) => { const a = sec && x.arrays(sec.body).find(q => q.name === n); return a ? a.get() : null; };
  if (x.ft === 'UnstructuredGrid') { const C = x.section('Cells', piece.body), cn = named(C, 'connectivity'), of = named(C, 'offsets'), ty = named(C, 'types'); if (!cn || !of || !ty) fail('VTU', 'Cells need connectivity, offsets and types arrays'); for (let i = 0, s = 0; i < of.length; i++) { vtkCell(ty[i], Array.from(cn.subarray(s, of[i])), elems, skip); s = of[i]; } }
  else for (const [sec, vt] of [['Polys', 7], ['Strips', 6], ['Lines', 4]]) { const S = x.section(sec, piece.body), cn = named(S, 'connectivity'), of = named(S, 'offsets'); if (cn && of) for (let i = 0, s = 0; i < of.length; i++) { const ids = Array.from(cn.subarray(s, of[i])); vtkCell(vt === 7 && ids.length === 4 ? 9 : vt, ids, elems, skip); s = of[i]; } }
  const pieces = (x.t.match(/<Piece\b/g) || []).length, fields = [...piece.body.matchAll(/<(PointData|CellData)\b[^>]*>([\s\S]*?)<\/\1>/g)].flatMap(m => x.arrays(m[2]).map(a => `${m[1]} ${a.name}`));
  const g = femToGeometry(pts, elems, { meta: { type: x.ft, version: x.version, encoding: x.encoding, fields } }); g.notes.push('VTK stores no length unit: coordinates taken as stored');
  if (fields.length) g.lost.push(`data arrays not attached to the geometry: ${fields.slice(0, 12).join(', ')}`); if (pieces > 1) g.lost.push(`only the first of ${pieces} pieces was imported`);
  if (Object.keys(skip).length) g.lost.push(`unsupported cell types skipped: ${Object.entries(skip).map(([k, v]) => `type ${k} × ${v}`).join(', ')}`);
  return g;
}
export const parseVTU = input => vtkXmlMesh(input, 'UnstructuredGrid');
export const parseVTP = input => vtkXmlMesh(input, 'PolyData');
export function parseVTI(input) {
  const x = vtkXml(input, 'ImageData'), im = x.section('ImageData'); if (!im) fail('VTI', 'no <ImageData> element');
  const ext = nums(attr(im.attrs, 'WholeExtent') || ''), sp = nums(attr(im.attrs, 'Spacing') || '1 1 1'), og = nums(attr(im.attrs, 'Origin') || '0 0 0'); if (ext.length < 6) fail('VTI', 'WholeExtent missing');
  const cd = x.section('CellData', im.body), pd = x.section('PointData', im.body), ca = cd ? x.arrays(cd.body) : [], pa = pd ? x.arrays(pd.body) : [], pick = ca.find(a => a.nc === 1) || pa.find(a => a.nc === 1);
  if (!pick) fail('VTI', 'no single-component data array to interpret as phases or grey levels');
  const cell = ca.includes(pick), d = [0, 1, 2].map(i => max(1, ext[2 * i + 1] - ext[2 * i] + (cell ? 0 : 1)));
  const vx = scalarToVoxels(pick.get(), d[0], d[1], d[2], { voxel: sp[0], spacing: sp, origin: og, units: '' }), all = [...ca, ...pa].map(a => a.name);
  return Object.assign(vx, { meta: { type: 'ImageData', array: pick.name, centring: cell ? 'cell' : 'point', encoding: x.encoding, arrays: all }, lost: all.length > 1 ? [`only array "${pick.name}" imported; others present: ${all.filter(n => n !== pick.name).join(', ')}`] : [], notes: [...vx.notes, 'VTK stores no length unit: spacing taken as stored'] });
}
export function toVTU(g) {
  const lat = g.kind === 'lattice'; if (!lat) needMesh(g, 'VTU');
  const P = lat ? g.nodes : g.vertices, C = lat ? g.edges : g.faces, k = lat ? 2 : 3;
  return `<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n<UnstructuredGrid>\n<Piece NumberOfPoints="${P.length}" NumberOfCells="${C.length}">\n<Points>\n<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n${P.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n')}\n</DataArray>\n</Points>\n<Cells>\n<DataArray type="Int32" Name="connectivity" format="ascii">\n${C.map(c => c.join(' ')).join('\n')}\n</DataArray>\n<DataArray type="Int32" Name="offsets" format="ascii">\n${C.map((_, i) => (i + 1) * k).join(' ')}\n</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii">\n${C.map(() => lat ? 3 : 5).join(' ')}\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n`;
}
export function toVTI(g) {
  if (g.kind !== 'voxels') fail('VTI', `needs a voxels geometry (got ${g.kind})`);
  const N = g.nx * g.ny * g.nz, sp = g.spacing || [g.voxel || 1, g.voxel || 1, g.voxel || 1], o = g.origin || [0, 0, 0], rows = []; for (let i = 0; i < N; i += 40) rows.push(Array.prototype.join.call(g.data.subarray(i, min(N, i + 40)), ' '));
  const e = `0 ${g.nx} 0 ${g.ny} 0 ${g.nz}`;
  return `<?xml version="1.0"?>\n<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n<ImageData WholeExtent="${e}" Origin="${o.map(N9).join(' ')}" Spacing="${sp.map(N9).join(' ')}">\n<Piece Extent="${e}">\n<CellData Scalars="phase">\n<DataArray type="UInt8" Name="phase" format="ascii">\n${rows.join('\n')}\n</DataArray>\n</CellData>\n</Piece>\n</ImageData>\n</VTKFile>\n`;
}

// ---------------------------------------------------------------- Gmsh MSH
const GMSH = { 1: ['line', 2], 2: ['tri', 3], 3: ['quad', 4], 4: ['tet', 4], 5: ['hex', 8], 6: ['wedge', 6], 7: ['pyr', 5], 8: ['line', 3], 9: ['tri', 6], 10: ['quad', 9], 11: ['tet', 10], 12: ['hex', 27], 13: ['wedge', 18], 14: ['pyr', 14], 15: ['point', 1], 16: ['quad', 8], 17: ['hex', 20], 18: ['wedge', 15], 19: ['pyr', 13], 20: ['tri', 9], 21: ['tri', 10], 26: ['line', 4], 29: ['tet', 20] };
export function toMSH(g) {
  const lat = g.kind === 'lattice'; if (!lat) needMesh(g, 'MSH'); const P = lat ? g.nodes : g.vertices, C = lat ? g.edges : g.faces;
  return `$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n${lat ? 1 : 2} 1 "${lat ? 'struts' : 'surface'}"\n$EndPhysicalNames\n$Nodes\n${P.length}\n${P.map((p, i) => `${i + 1} ${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n')}\n$EndNodes\n$Elements\n${C.length}\n${C.map((c, i) => `${i + 1} ${lat ? 1 : 2} 2 1 1 ${c.map(k => k + 1).join(' ')}`).join('\n')}\n$EndElements\n`;
}
export function parseMSH(input) {
  const t = toText(input), sec = name => { const a = t.indexOf('$' + name), e = t.indexOf('$End' + name); return a >= 0 && e > a ? lines(t.slice(t.indexOf('\n', a) + 1, e)).map(l => l.trim()).filter(Boolean) : null; };
  const mf = sec('MeshFormat'); if (!mf) fail('MSH', 'no $MeshFormat section (not a Gmsh mesh; Fluent and ANSYS also use .msh)');
  const [ver, bin] = nums(mf[0]); if (bin) fail('MSH', 'binary Gmsh files are not supported by this reader: save with "Mesh.Binary = 0"');
  const N = sec('Nodes'), E = sec('Elements'); if (!N || !E) fail('MSH', 'missing $Nodes or $Elements section');
  const idx = new Map(), pts = [], elems = [], skipped = {}, phys = new Map(), names = {};
  for (const l of (sec('PhysicalNames') || []).slice(1)) { const m = /^(\d+)\s+(\d+)\s+"(.*)"/.exec(l); if (m) names[m[2]] = { dim: +m[1], name: m[3] }; }
  const el = (type, ids, tag) => { const g = GMSH[type]; if (!g || g[0] === 'point') { if (!g) skipped[type] = (skipped[type] || 0) + 1; return; } const n = ids.map(i => idx.get(i)); if (n.some(v => v === undefined)) fail('MSH', 'element references an undefined node'); elems.push({ type: g[0], n, tag }); phys.set(tag, (phys.get(tag) || 0) + 1); };
  if (ver < 4) {
    const n = +N[0]; if (N.length < n + 1) fail('MSH', `$Nodes declares ${n} nodes but holds ${N.length - 1}`);
    for (let i = 1; i <= n; i++) { const a = nums(N[i]); idx.set(a[0], pts.length); pts.push([a[1], a[2], a[3] || 0]); }
    const ne = +E[0]; for (let i = 1; i <= ne && i < E.length; i++) { const a = nums(E[i]), nt = a[2]; el(a[1], a.slice(3 + nt), nt ? a[3] : 0); }
  } else {
    let li = 1; const nbk = nums(N[0])[0], v41 = ver >= 4.1;
    for (let k = 0; k < nbk; k++) { const h = nums(N[li++]), cnt = h[3]; if (v41) { const tags = []; for (let i = 0; i < cnt; i++) tags.push(+N[li++]); for (let i = 0; i < cnt; i++) { const a = nums(N[li++]); idx.set(tags[i], pts.length); pts.push([a[0], a[1], a[2] || 0]); } } else for (let i = 0; i < cnt; i++) { const a = nums(N[li++]); idx.set(a[0], pts.length); pts.push([a[1], a[2], a[3] || 0]); } }
    li = 1; const ebk = nums(E[0])[0];
    for (let k = 0; k < ebk; k++) { const h = nums(E[li++]), tag = v41 ? h[1] : h[0], type = h[2], cnt = h[3]; for (let i = 0; i < cnt; i++) { const a = nums(E[li++]); el(type, a.slice(1), tag); } }
  }
  const dim = Math.max(...elems.map(e => e.type === 'line' ? 1 : ['tri', 'quad'].includes(e.type) ? 2 : 3), 0);
  const use = dim === 3 ? elems.filter(e => !['tri', 'quad', 'line'].includes(e.type)) : dim === 2 ? elems.filter(e => e.type !== 'line') : elems;
  const g = femToGeometry(pts, use, { meta: { version: ver, dimension: dim } });
  g.regions = [...phys].map(([tag, n]) => ({ tag, name: v41Name(names, tag, ver), elements: n }));
  if (use.length < elems.length) g.notes.push(`${elems.length - use.length} lower-dimensional boundary elements (physical surfaces/curves) kept as regions, not as geometry`);
  g.notes.push('Gmsh stores no length unit: coordinates taken as stored'); if (ver >= 4) g.notes.push('MSH 4: the entity tag is used as region id (physical groups are defined in $Entities)');
  if (Object.keys(skipped).length) g.lost.push(`unsupported element types skipped: ${Object.entries(skipped).map(([k, v]) => `type ${k} × ${v}`).join(', ')}`);
  return g;
}
const v41Name = (names, tag, ver) => (ver < 4 && names[tag] ? names[tag].name : `tag ${tag}`);

// ---------------------------------------------------------------- Abaqus INP
function abqType(t) {
  let m; t = t.toUpperCase();
  if ((m = /^(?:D|DC|AC|Q)?C3D(\d+)/.exec(t))) return ({ 4: 'tet', 10: 'tet', 8: 'hex', 20: 'hex', 27: 'hex', 6: 'wedge', 15: 'wedge', 5: 'pyr', 13: 'pyr' })[m[1]] || null;
  if (/^(B\d\d|T\dD\d|PIPE\d\d|FRAME\dD|CONN\dD2|SPRINGA|DASHPOTA)/.test(t)) return 'line';
  if ((m = /^(?:S|SC|STRI|M3D|CPS|CPE|CPEG|CAX|DS|R3D|DC2D|DCAX|SFM3D|AC2D|ACAX|CGAX)(\d+)/.exec(t))) { const n = m[1] === '65' ? 6 : +m[1][0]; return n === 3 || n === 6 ? 'tri' : n === 4 || n === 8 || n === 9 ? 'quad' : null; }
  return null;
}
export function toINP(g) {
  const lat = g.kind === 'lattice'; if (!lat) needMesh(g, 'INP'); const P = lat ? g.nodes : g.vertices, C = lat ? g.edges : g.faces;
  return `*HEADING\nDesigner Materials ${lat ? 'strut lattice (B31 beams)' : 'surface mesh (S3 shells)'}\n** length unit: ${g.units || 'unspecified'}${lat ? `; strut radius ${N9(g.radius || 0)}` : ''}\n*NODE\n${P.map((p, i) => `${i + 1}, ${N9(p[0])}, ${N9(p[1])}, ${N9(p[2] || 0)}`).join('\n')}\n*ELEMENT, TYPE=${lat ? 'B31' : 'S3'}, ELSET=${lat ? 'STRUTS' : 'SURFACE'}\n${C.map((c, i) => `${i + 1}, ${c.map(k => k + 1).join(', ')}`).join('\n')}\n`;
}
export function parseINP(input) {
  const L = lines(toText(input)), pts = [], idx = new Map(), elems = [], sets = {}, unknown = {}; let mode = '', etype = null, part = '', cur = null, gen = false, pend = null, inst = 0, heading = '', wantHead = false, ename = '';
  for (const raw of L) {
    const l = raw.trim(); if (!l || l.startsWith('**')) continue;
    if (l[0] === '*') {
      const f = l.split(',').map(s => s.trim()), kw = f[0].toUpperCase().replace(/\s+/g, ''), par = {}; for (const q of f.slice(1)) { const [k, v] = q.split('='); par[k.trim().toUpperCase()] = (v ?? '').trim(); }
      mode = ''; pend = null; wantHead = false;
      if (kw === '*NODE' && !('FILE' in par)) mode = 'node';
      else if (kw === '*ELEMENT' && par.TYPE) { mode = 'elem'; ename = par.TYPE.toUpperCase(); etype = abqType(ename); if (!etype) unknown[ename] = 0; if (par.ELSET) sets[`elset ${par.ELSET}`] = sets[`elset ${par.ELSET}`] || 0; cur = par.ELSET ? `elset ${par.ELSET}` : null; }
      else if (kw === '*NSET' || kw === '*ELSET') { mode = 'set'; cur = `${kw.slice(1).toLowerCase()} ${par.NSET || par.ELSET || '?'}`; sets[cur] = sets[cur] || 0; gen = 'GENERATE' in par; }
      else if (kw === '*PART') part = par.NAME || `part${Object.keys(sets).length}`; else if (kw === '*ENDPART') part = ''; else if (kw === '*INSTANCE') inst++; else if (kw === '*HEADING') wantHead = true;
      continue;
    }
    if (wantHead) { heading = heading || l; continue; }
    if (mode === 'node') { const a = l.split(',').map(Number); if (!Number.isFinite(a[0]) || !Number.isFinite(a[1])) fail('INP', `bad *NODE line "${l}"`); idx.set(part + '|' + a[0], pts.length); pts.push([a[1], a[2] || 0, a[3] || 0]); }
    else if (mode === 'elem') {
      const a = l.split(',').map(s => s.trim()).filter(s => s !== '').map(Number); pend = pend ? pend.concat(a) : a; if (l.endsWith(',')) continue;
      const ids = pend.slice(1); pend = null; if (cur) sets[cur]++;
      if (!etype) { unknown[ename]++; continue; }
      const n = ids.map(i => idx.get(part + '|' + i)); if (n.some(v => v === undefined)) fail('INP', `element ${ids.length ? '' : '(empty) '}on line "${l}" references a node that is not defined before it`);
      elems.push({ type: etype, n });
    } else if (mode === 'set') { const a = l.split(',').map(s => s.trim()).filter(Boolean); if (gen) { const [s, e, st = 1] = a.map(Number); sets[cur] += max(0, floor((e - s) / (st || 1)) + 1); } else sets[cur] += a.length; }
  }
  if (!pts.length) fail('INP', 'no *NODE data found (if this is a CP2K input, it is an atomistic file, not an Abaqus mesh)');
  const g = femToGeometry(pts, elems, { sets, meta: { heading, instances: inst } });
  g.notes.push('Abaqus input carries no unit system: coordinates taken as stored');
  if (inst) g.lost.push(`${inst} *INSTANCE placements not applied: parts are shown in their own coordinates`);
  const uk = Object.entries(unknown).filter(([, v]) => v); if (uk.length) g.lost.push(`element types not mapped to geometry: ${uk.map(([k, v]) => `${k} × ${v}`).join(', ')}`);
  g.lost.push('materials, sections, boundary conditions, loads and steps are not imported');
  return g;
}

// ---------------------------------------------------------------- NASTRAN bulk data
const nasNum = s => { s = String(s).trim(); if (s === '') return NaN; return Number(s.replace(/[dD]/, 'e').replace(/([0-9.])([+-]\d+)$/, '$1e$2')); };
const NASEL = { CTRIA3: ['tri', 3, 2], CTRIA6: ['tri', 6, 2], CTRIAR: ['tri', 3, 2], CQUAD4: ['quad', 4, 2], CQUAD8: ['quad', 8, 2], CQUADR: ['quad', 4, 2], CTETRA: ['tet', 10, 2], CHEXA: ['hex', 20, 2], CPENTA: ['wedge', 15, 2], CPYRAM: ['pyr', 13, 2], CBAR: ['line', 2, 2], CBEAM: ['line', 2, 2], CROD: ['line', 2, 2], CTUBE: ['line', 2, 2], CONROD: ['line', 2, 1] };
export function toBDF(g) {
  const lat = g.kind === 'lattice'; if (!lat) needMesh(g, 'BDF'); const P = lat ? g.nodes : g.vertices, C = lat ? g.edges : g.faces;
  const f8 = v => String(v).padEnd(8).slice(0, 8), i16 = v => String(v).padEnd(16), e16 = v => { const s = Number(v).toExponential(9).toUpperCase(); return s.padStart(16); };
  const o = [`$ Designer Materials ${lat ? 'strut lattice (CROD)' : 'surface mesh (CTRIA3)'}; length unit: ${g.units || 'unspecified'}`, 'BEGIN BULK'];
  P.forEach((p, i) => o.push(`GRID*   ${i16(i + 1)}${i16('')}${e16(p[0])}${e16(p[1])}*`, `*       ${e16(p[2] || 0)}`));
  C.forEach((c, i) => o.push(f8(lat ? 'CROD' : 'CTRIA3') + f8(i + 1) + f8(1) + c.map(k => f8(k + 1)).join('')));
  o.push('ENDDATA'); return o.join('\n') + '\n';
}
export function parseBDF(input) {
  const raw = lines(toText(input).replace(/\t/g, '        ')), bb = raw.findIndex(l => /^\s*BEGIN\s+BULK/i.test(l)), cards = []; let cur = null;
  for (let i = bb + 1; i < raw.length; i++) {
    let l = raw[i]; if (!l.trim() || l[0] === '$') continue; if (/^ENDDATA/i.test(l)) break; l = l.replace(/\$.*$/, '');
    let f, cont;
    if (l.includes(',')) { f = l.split(',').map(s => s.trim()); cont = f[0] === '' || /^[+*]/.test(f[0]); if (!cont) { cur = [f[0].toUpperCase().replace('*', ''), ...f.slice(1)]; cards.push(cur); } else if (cur) cur.push(...f.slice(1)); while (cur && /^[+*]/.test(cur[cur.length - 1] || '')) cur.pop(); continue; }
    const name = l.slice(0, 8).trim(), large = name.endsWith('*') || (name.startsWith('*')); cont = name === '' || /^[+*]/.test(name);
    f = []; if (large) for (let c = 8; c < 72; c += 16) f.push(l.slice(c, c + 16).trim()); else for (let c = 8; c < 72; c += 8) f.push(l.slice(c, c + 8).trim());
    if (!cont) { cur = [name.toUpperCase().replace('*', ''), ...f]; cards.push(cur); } else if (cur) cur.push(...f);
  }
  const pts = [], idx = new Map(), elems = [], other = {}; let cp = 0;
  for (const c of cards) if (c[0] === 'GRID') { const id = parseInt(c[1], 10), x = [nasNum(c[3]), nasNum(c[4]), nasNum(c[5])].map(v => Number.isFinite(v) ? v : 0); if (!Number.isFinite(id)) fail('BDF', `bad GRID card "${c.join(',')}"`); if (c[2] && +c[2]) cp++; idx.set(id, pts.length); pts.push(x); }
  for (const c of cards) { const d = NASEL[c[0]]; if (!d) { if (c[0] !== 'GRID') other[c[0]] = (other[c[0]] || 0) + 1; continue; } const ids = c.slice(1 + d[2], 1 + d[2] + d[1]).map(s => parseInt(s, 10)).filter(Number.isFinite), need = CORNERS[d[0]]; if (ids.length < need) fail('BDF', `${c[0]} ${c[1]} has ${ids.length} grid points, needs ${need}`); const n = ids.slice(0, d[0] === 'line' ? 2 : ids.length).map(i => idx.get(i)); if (n.slice(0, need).some(v => v === undefined)) fail('BDF', `${c[0]} ${c[1]} references an undefined GRID`); elems.push({ type: d[0], n: n.filter(v => v !== undefined) }); }
  if (!pts.length) fail('BDF', 'no GRID cards found');
  const g = femToGeometry(pts, elems, { meta: { cards: cards.length, otherCards: other } });
  g.notes.push('NASTRAN bulk data carries no unit system: coordinates taken as stored'); if (cp) g.lost.push(`${cp} GRID points are defined in a local coordinate system (CP ≠ 0) that was not transformed`);
  const props = Object.keys(other).filter(k => /^(P|MAT|SPC|MPC|FORCE|LOAD|RBE|CORD)/.test(k)); if (props.length) g.lost.push(`property, material, constraint and load cards not imported (${props.slice(0, 10).join(', ')})`);
  return g;
}

// ---------------------------------------------------------------- I-DEAS universal file, SU2
export function parseUNV(input) {
  const L = lines(toText(input)), pts = [], idx = new Map(), elems = [], sets = []; let i = 0, units = '', unknown = 0;
  const isSep = l => l !== undefined && l.trim() === '-1';
  if (!L.some(isSep)) fail('UNV', 'no "    -1" dataset delimiters found');
  const dnum = s => Number(s.replace(/[dD]/g, 'e'));
  while (i < L.length) {
    if (!isSep(L[i])) { i++; continue; } const ds = parseInt(L[i + 1], 10); i += 2; const s = i; while (i < L.length && !isSep(L[i])) i++; const blk = L.slice(s, i); i++;
    if (ds === 2411) for (let k = 0; k + 1 < blk.length; k += 2) { const id = parseInt(blk[k], 10), c = blk[k + 1].trim().split(/\s+/).map(dnum); idx.set(id, pts.length); pts.push([c[0], c[1], c[2] || 0]); }
    else if (ds === 2412) for (let k = 0; k < blk.length;) {
      const h = nums(blk[k++]), fe = h[1], nn = h[5], beam = fe >= 11 && fe <= 24; if (beam) k++; const ids = []; while (ids.length < nn && k < blk.length) ids.push(...nums(blk[k++]));
      let type = beam ? 'line' : ({ 3: 'tri', 6: fe === 112 ? 'wedge' : 'tri', 4: fe === 111 ? 'tet' : 'quad', 8: fe === 115 ? 'hex' : 'quad', 9: 'quad', 10: 'tet', 15: 'wedge', 20: 'hex', 5: 'pyr' })[nn]; if (!type) { unknown++; continue; }
      let n = ids; if (!beam && nn === 6 && type === 'tri') n = [ids[0], ids[2], ids[4]]; else if (nn === 8 && type === 'quad') n = [ids[0], ids[2], ids[4], ids[6]]; else if (nn === 10 && type === 'tet') n = [ids[0], ids[2], ids[4], ids[9]]; else if (beam) n = [ids[0], ids[ids.length - 1]];
      const m = n.map(q => idx.get(q)); if (m.some(v => v === undefined)) fail('UNV', 'element references an undefined node (dataset 2412 before 2411?)'); elems.push({ type, n: m });
    } else if (ds === 164 && blk.length) units = blk[0].slice(10, 30).trim() || units;
    else if ((ds === 2467 || ds === 2477 || ds === 2452) && blk.length > 1) for (let k = 0; k + 1 < blk.length;) { const h = nums(blk[k]), n = h[h.length - 1], name = (blk[k + 1] || '').trim(); sets.push({ name, entities: n }); k += 2 + Math.ceil(n / 2); }
  }
  if (!pts.length) fail('UNV', 'no node dataset (2411) found');
  const g = femToGeometry(pts, elems, { meta: { unitsLine: units, groups: sets } });
  g.notes.push(units ? `unit system line from dataset 164: "${units}"` : 'no unit dataset (164): coordinates taken as stored'); if (unknown) g.lost.push(`${unknown} elements of unsupported shape skipped`);
  return g;
}
export function parseSU2(input) {
  const L = lines(toText(input)).map(l => l.replace(/%.*/, '').trim()).filter(Boolean), pts = [], elems = [], markers = [], NN = { 3: 2, 5: 3, 9: 4, 10: 4, 12: 8, 13: 6, 14: 5 }; let dim = 3, skip = {};
  for (let i = 0; i < L.length; i++) {
    const m = /^(\w+)\s*=\s*(.*)$/.exec(L[i]); if (!m) continue; const k = m[1].toUpperCase(), v = m[2].trim();
    if (k === 'NDIME') dim = +v;
    else if (k === 'NPOIN') { const n = parseInt(v, 10); for (let j = 0; j < n; j++) { const a = nums(L[++i] || ''); if (a.length < dim) fail('SU2', `point line ${j + 1} is incomplete`); pts.push(dim === 2 ? [a[0], a[1], 0] : [a[0], a[1], a[2]]); } }
    else if (k === 'NELEM') { const n = parseInt(v, 10); for (let j = 0; j < n; j++) { const a = nums(L[++i] || ''), nn = NN[a[0]]; if (!nn) { skip[a[0]] = (skip[a[0]] || 0) + 1; continue; } vtkCell(a[0], a.slice(1, 1 + nn), elems, skip); } }
    else if (k === 'MARKER_TAG') markers.push({ name: v, elements: 0 });
    else if (k === 'MARKER_ELEMS') { const n = parseInt(v, 10); if (markers.length) markers[markers.length - 1].elements = n; i += n; }
  }
  if (!pts.length) fail('SU2', 'no NPOIN block found'); if (elems.some(e => e.n.some(q => !(q >= 0 && q < pts.length)))) fail('SU2', 'element references an undefined point');
  const g = femToGeometry(pts, elems, { meta: { dimension: dim } }); g.regions = markers.map(m => ({ name: m.name, role: 'boundary marker', elements: m.elements }));
  g.notes.push('SU2 stores no length unit: coordinates taken as stored'); if (dim === 2) g.notes.push('2-D mesh placed in the z = 0 plane');
  return g;
}

// ---------------------------------------------------------------- XDMF / PVD
export function toXDMF(g) {
  const head = '<?xml version="1.0"?>\n<Xdmf Version="3.0">\n<Domain>\n';
  if (g.kind === 'voxels') { const sp = g.spacing || [g.voxel || 1, g.voxel || 1, g.voxel || 1], o = g.origin || [0, 0, 0], N = g.nx * g.ny * g.nz, rows = []; for (let i = 0; i < N; i += 40) rows.push(Array.prototype.join.call(g.data.subarray(i, min(N, i + 40)), ' '));
    return head + `<Grid Name="voxels" GridType="Uniform">\n<Topology TopologyType="3DCoRectMesh" Dimensions="${g.nz + 1} ${g.ny + 1} ${g.nx + 1}"/>\n<Geometry GeometryType="ORIGIN_DXDYDZ">\n<DataItem Dimensions="3" NumberType="Float" Precision="8" Format="XML">${N9(o[2])} ${N9(o[1])} ${N9(o[0])}</DataItem>\n<DataItem Dimensions="3" NumberType="Float" Precision="8" Format="XML">${N9(sp[2])} ${N9(sp[1])} ${N9(sp[0])}</DataItem>\n</Geometry>\n<Attribute Name="phase" AttributeType="Scalar" Center="Cell">\n<DataItem Dimensions="${g.nz} ${g.ny} ${g.nx}" NumberType="UChar" Precision="1" Format="XML">\n${rows.join('\n')}\n</DataItem>\n</Attribute>\n</Grid>\n</Domain>\n</Xdmf>\n`; }
  const lat = g.kind === 'lattice'; if (!lat) needMesh(g, 'XDMF'); const P = lat ? g.nodes : g.vertices, C = lat ? g.edges : g.faces;
  return head + `<Grid Name="${lat ? 'lattice' : 'mesh'}" GridType="Uniform">\n<Topology TopologyType="${lat ? 'Polyline" NodesPerElement="2' : 'Triangle'}" NumberOfElements="${C.length}">\n<DataItem Dimensions="${C.length} ${lat ? 2 : 3}" NumberType="Int" Format="XML">\n${C.map(c => c.join(' ')).join('\n')}\n</DataItem>\n</Topology>\n<Geometry GeometryType="XYZ">\n<DataItem Dimensions="${P.length} 3" NumberType="Float" Precision="8" Format="XML">\n${P.map(p => `${N9(p[0])} ${N9(p[1])} ${N9(p[2] || 0)}`).join('\n')}\n</DataItem>\n</Geometry>\n</Grid>\n</Domain>\n</Xdmf>\n`;
}
// Inventory of grids, time steps and referenced heavy-data files. Inline (Format="XML") data are reconstructed;
// HDF5-backed items are read when the companion file is supplied in `files` ({name: bytes}).
export function parseXDMF(input, { files = {} } = {}) {
  const t = toText(input); if (!/<Xdmf\b/.test(t)) fail('XDMF', 'no <Xdmf> root element');
  const grids = [...t.matchAll(/<Grid\b([^>]*)>/g)].map(m => ({ name: attr(m[1], 'Name'), type: attr(m[1], 'GridType') || 'Uniform', collection: attr(m[1], 'CollectionType') }));
  const times = [...t.matchAll(/<Time\b([^>]*)\/?>/g)].map(m => +attr(m[1], 'Value')).filter(Number.isFinite), refs = new Set(), items = [];
  const reI = /<DataItem\b([^>]*?)(?:\/>|>([\s\S]*?)<\/DataItem>)/g; let m;
  while ((m = reI.exec(t))) { const fmt = (attr(m[1], 'Format') || 'XML').toUpperCase(), body = (m[2] || '').trim(); if (fmt !== 'XML' && body) refs.add(body.split(':')[0].trim()); items.push({ at: m.index, fmt, dims: nums(attr(m[1], 'Dimensions') || ''), body }); }
  const meta = { version: attr((/<Xdmf\b[^>]*>/.exec(t) || [''])[0], 'Version'), grids, timeSteps: times, referencedFiles: [...refs], attributes: [...t.matchAll(/<Attribute\b([^>]*)>/g)].map(a => attr(a[1], 'Name')) };
  const notes = [], lost = [], h5 = {};
  const load = it => { if (!it) return null; if (it.fmt === 'XML') return it.body ? Float64Array.from(nums(it.body)) : null; if (it.fmt === 'HDF') { const [fn, path] = it.body.split(':').map(s => s.trim()), key = Object.keys(files).find(k => k === fn || k.endsWith('/' + fn)); if (!key) return null; h5[key] = h5[key] || h5Open(toBytes(files[key])); const d = h5[key].datasets.find(q => q.path === path || q.path === '/' + path.replace(/^\//, '')); return d && d.readable ? d.read() : null; } return null; };
  const first = (tag, from = 0) => { const mm = new RegExp('<' + tag + '\\b([^>]*?)(/?)>', 'g'); mm.lastIndex = from; const r = mm.exec(t); if (!r) return null; const end = r[2] ? r.index + r[0].length : t.indexOf('</' + tag + '>', r.index); return { attrs: r[1], at: r.index, end: end < 0 ? t.length : end }; };
  const topo = first('Topology'), geo = first('Geometry'), inside = (s) => s ? items.filter(it => it.at > s.at && it.at < s.end) : [];
  let g = null;
  if (topo && geo) {
    const tt = (attr(topo.attrs, 'TopologyType') || attr(topo.attrs, 'Type') || '').toLowerCase(), gt = (attr(geo.attrs, 'GeometryType') || attr(geo.attrs, 'Type') || 'XYZ').toUpperCase(), gi = inside(geo), ti = inside(topo);
    if (/corectmesh/.test(tt)) {
      const d = nums(attr(topo.attrs, 'Dimensions') || attr(topo.attrs, 'NumberOfElements') || ''), og = load(gi[0]), sp = load(gi[1]), at = first('Attribute'), ai = at ? inside(at)[0] : null, vals = load(ai);
      if (d.length === 3 && vals) { const cell = (attr(at.attrs, 'Center') || 'Node').toLowerCase() === 'cell', n = cell ? d.map(v => v - 1) : d; g = scalarToVoxels(vals, n[2], n[1], n[0], { voxel: sp ? sp[2] : 1, spacing: sp ? [sp[2], sp[1], sp[0]] : undefined, origin: og ? [og[2], og[1], og[0]] : undefined }); }
    } else {
      const NPE = { triangle: ['tri', 3], quadrilateral: ['quad', 4], tetrahedron: ['tet', 4], hexahedron: ['hex', 8], wedge: ['wedge', 6], pyramid: ['pyr', 5], polyline: ['line', 2], triangle_6: ['tri', 6], tetrahedron_10: ['tet', 10], edge: ['line', 2] }[tt], xyz = load(gi[0]), cn = load(ti[0]);
      if (NPE && xyz && cn) { const dimG = gt === 'XY' ? 2 : 3, pts = []; for (let i = 0; i + dimG - 1 < xyz.length; i += dimG) pts.push([xyz[i], xyz[i + 1], dimG === 3 ? xyz[i + 2] : 0]); const elems = []; for (let i = 0; i + NPE[1] - 1 < cn.length; i += NPE[1]) elems.push({ type: NPE[0], n: Array.from(cn.subarray(i, i + NPE[1])) }); if (elems.some(e => e.n.some(q => !(q >= 0 && q < pts.length)))) fail('XDMF', 'topology references an undefined point'); g = femToGeometry(pts, elems, {}); }
      else if (!NPE) notes.push(`topology type "${tt}" is not reconstructed`);
    }
  }
  if (grids.length > 1) lost.push(`${grids.length} grids / ${times.length} time steps listed: only the first grid is reconstructed`); if (meta.attributes.length) lost.push(`attribute arrays not attached: ${meta.attributes.slice(0, 10).join(', ')}`);
  if (!g) { if (refs.size) notes.push(`heavy data live in ${[...refs].join(', ')}: load that file together with the .xdmf to reconstruct the grid`); return { kind: null, meta, notes: [...notes, 'inventory only (grids, time steps, referenced files)'], lost: [...lost, 'geometry not reconstructed'] }; }
  return Object.assign(g, { meta, notes: [...(g.notes || []), ...notes, 'XDMF stores no length unit: coordinates taken as stored'], lost: [...(g.lost || []), ...lost] });
}
// ParaView collection: list of (time, file). Referenced .vtu/.vtp/.vti files are parsed when supplied in `files`.
export function parsePVD(input, { files = {}, maxFrames = 40 } = {}) {
  const t = toText(input); if (!/<VTKFile\b[^>]*type\s*=\s*"Collection"/.test(t)) fail('PVD', 'not a VTK Collection file');
  const sets = [...t.matchAll(/<DataSet\b([^>]*)\/?>/g)].map(m => ({ t: +attr(m[1], 'timestep') || 0, file: attr(m[1], 'file') || '', part: attr(m[1], 'part') }));
  const frames = [], missing = []; let first = null;
  for (const s of sets.slice(0, maxFrames)) { const key = Object.keys(files).find(k => k === s.file || k.endsWith('/' + s.file) || s.file.endsWith('/' + k)); if (!key) { missing.push(s.file); continue; } const g = /\.vti$/i.test(s.file) ? parseVTI(files[key]) : /\.vtp$/i.test(s.file) ? parseVTP(files[key]) : parseVTU(files[key]); first = first || g; frames.push(g.kind === 'voxels' ? { t: s.t, data: g.data } : g.kind === 'lattice' ? { t: s.t, nodes: g.nodes, edges: g.edges } : { t: s.t, vertices: g.vertices, faces: g.faces }); }
  const meta = { dataSets: sets, timeSteps: [...new Set(sets.map(s => s.t))], referencedFiles: sets.map(s => s.file) }, lost = [];
  if (sets.length > maxFrames) lost.push(`${sets.length - maxFrames} time steps beyond the first ${maxFrames} not loaded`);
  if (!first) return { kind: null, meta, frames: [], notes: [`collection of ${sets.length} data sets; the referenced files (${missing.slice(0, 5).join(', ')}${missing.length > 5 ? ', …' : ''}) must be loaded together with the .pvd`], lost: [...lost, 'geometry not reconstructed'] };
  return Object.assign(first, { frames, meta, notes: [...(first.notes || []), `${frames.length} of ${sets.length} time steps loaded as animation frames`], lost: [...(first.lost || []), ...lost, ...(missing.length ? [`${missing.length} referenced files were not supplied`] : [])] });
}

// ---------------------------------------------------------------- GDSII stream
function gdsReal(v) { const b = new Uint8Array(8); if (!v) return b; const neg = v < 0; v = abs(v); let e = 64; while (v >= 1) { v /= 16; e++; } while (v < 1 / 16) { v *= 16; e--; } b[0] = (neg ? 0x80 : 0) | e; for (let i = 1; i < 8; i++) { v *= 256; const d = floor(v); b[i] = d; v -= d; } return b; }
function gdsRealRead(b, o) { let m = 0; for (let i = 1; i < 8; i++) m = m * 256 + b[o + i]; return (b[o] & 0x80 ? -1 : 1) * m / 2 ** 56 * 16 ** ((b[o] & 0x7F) - 64); }
const layoutScale = (layout, dbu_m) => (unitToMetres(layout.units) ?? 1e-9) / dbu_m;
function gridPolys(layout, k) { // integer-grid polygons, consecutive duplicates removed
  const out = [];
  for (const q of layout.polygons || []) { let pts = q.pts.map(p => [round(p[0] * k), round(p[1] * k)]).filter((p, i, a) => i === 0 || p[0] !== a[i - 1][0] || p[1] !== a[i - 1][1]); while (pts.length > 1 && pts[0][0] === pts[pts.length - 1][0] && pts[0][1] === pts[pts.length - 1][1]) pts.pop(); if (pts.length > 8000) { const st = Math.ceil(pts.length / 8000); pts = pts.filter((_, i) => i % st === 0); } if (pts.length >= 3) out.push({ layer: max(0, min(32767, round(q.layer ?? 1))), datatype: max(0, min(32767, round(q.datatype ?? 0))), pts }); }
  return out;
}
// user unit `unit_nm` (default 1 µm) and database unit `dbu_nm` (default 1 nm). A layout produced by array() without
// per-cell rotation is written hierarchically: structure UNIT placed by one AREF in structure TOP.
export function toGDS(layout, { unit_nm = 1000, dbu_nm = 1, name = 'TOP', libname = 'DESIGNER_MATERIALS', hierarchical = true, timestamp } = {}) {
  if (!layout || layout.kind !== 'layout') fail('GDSII', `needs a layout geometry (got ${layout && layout.kind})`);
  const k = layoutScale(layout, dbu_nm * 1e-9), w = new Writer(1 << 16), d = timestamp ? new Date(timestamp) : new Date(), date = [d.getFullYear(), d.getMonth() + 1, d.getDate(), d.getHours(), d.getMinutes(), d.getSeconds()];
  const rec = (type, dt, pl = new Uint8Array(0)) => { const pad = pl.length % 2; w.u16(4 + pl.length + pad, false).u8(type).u8(dt).bytes(pl); if (pad) w.u8(0); };
  const i2 = v => { const b = new Uint8Array(v.length * 2), dv = new DataView(b.buffer); v.forEach((x, i) => dv.setInt16(2 * i, x, false)); return b; }, i4 = v => { const b = new Uint8Array(v.length * 4), dv = new DataView(b.buffer); v.forEach((x, i) => dv.setInt32(4 * i, x, false)); return b; };
  const sname = s => ENC.encode(String(s).replace(/[^A-Za-z0-9_?$]/g, '_').slice(0, 32) || 'CELL');
  const boundary = q => { rec(0x08, 0); rec(0x0D, 2, i2([q.layer])); rec(0x0E, 2, i2([q.datatype])); rec(0x10, 3, i4([...q.pts.flat(), q.pts[0][0], q.pts[0][1]])); rec(0x11, 0); };
  const struct = (nm, body) => { rec(0x05, 2, i2([...date, ...date])); rec(0x06, 6, sname(nm)); body(); rec(0x07, 0); };
  rec(0x00, 2, i2([600])); rec(0x01, 2, i2([...date, ...date])); rec(0x02, 6, sname(libname));
  const un = new Uint8Array(16); un.set(gdsReal(dbu_nm / unit_nm), 0); un.set(gdsReal(dbu_nm * 1e-9), 8); rec(0x03, 5, un);
  const h = hierarchical && layout.hierarchy;
  if (h && h.nx <= 32767 && h.ny <= 32767) {
    struct('UNIT', () => gridPolys({ polygons: h.unit, units: layout.units }, k).forEach(boundary));
    const ox = round(h.origin[0] * k), oy = round(h.origin[1] * k);
    struct(name, () => { rec(0x0B, 0); rec(0x12, 6, sname('UNIT')); rec(0x13, 2, i2([h.nx, h.ny])); rec(0x10, 3, i4([ox, oy, ox + round(h.nx * h.px * k), oy, ox, oy + round(h.ny * h.py * k)])); rec(0x11, 0); });
  } else struct(name, () => gridPolys(layout, k).forEach(boundary));
  rec(0x04, 0); return w.out();
}
const affMul = (A, B) => [A[0] * B[0] + A[2] * B[1], A[1] * B[0] + A[3] * B[1], A[0] * B[2] + A[2] * B[3], A[1] * B[2] + A[3] * B[3], A[0] * B[4] + A[2] * B[5] + A[4], A[1] * B[4] + A[3] * B[5] + A[5]];
const affPt = (T, p) => [T[0] * p[0] + T[2] * p[1] + T[4], T[1] * p[0] + T[3] * p[1] + T[5]];
const placeT = (x, y, angDeg, mag, flip) => { const a = angDeg * PI / 180, c = cos(a) * mag, s = sin(a) * mag; return flip ? [c, s, s, -c, x, y] : [c, s, -s, c, x, y]; };
// Flatten a cell hierarchy {name → {polys, refs}} from its top cells; caps the output polygon count.
function flattenCells(cells, scale, cap = 200000) {
  const referenced = new Set(); for (const c of cells.values()) for (const r of c.refs) referenced.add(r.cell);
  const tops = [...cells.keys()].filter(k => !referenced.has(k)), polygons = [], layers = {}; let truncated = false, missing = new Set();
  const emit = (key, T, depth) => {
    const c = cells.get(key); if (!c) { missing.add(String(key)); return; } if (depth > 24) return;
    for (const q of c.polys) { if (polygons.length >= cap) { truncated = true; return; } const lk = `${q.layer}/${q.datatype}`; layers[lk] = (layers[lk] || 0) + 1; polygons.push({ layer: q.layer, datatype: q.datatype, pts: q.pts.map(p => { const r = affPt(T, p); return [r[0] * scale, r[1] * scale]; }) }); }
    for (const r of c.refs) for (const o of r.offs) { if (polygons.length >= cap) { truncated = true; return; } emit(r.cell, affMul(T, placeT(r.x + o[0], r.y + o[1], r.angle, r.mag, r.flip)), depth + 1); }
  };
  for (const tcell of (tops.length ? tops : [...cells.keys()].slice(0, 1))) emit(tcell, [1, 0, 0, 1, 0, 0], 0);
  return { polygons, layers, tops, truncated, missing: [...missing] };
}
export function parseGDS(input) {
  const b = toBytes(input), v = view(b); if (b.length < 6 || b[2] !== 0x00 || b[3] !== 0x02) fail('GDSII', 'file does not start with a HEADER record (00 06 00 02 …)');
  const cells = new Map(), meta = { version: v.getInt16(4, false) }, labels = []; let p = 0, cell = null, el = null, uu = 1e-3, dbu = 1e-9, npath = 0, roundEnds = 0, nref = 0, bigRef = false;
  const str = (o, n) => latin(b, o, n).replace(/\0+$/, '');
  while (p + 4 <= b.length) {
    const len = v.getUint16(p, false), type = b[p + 2], o = p + 4, n = len - 4; if (len < 4) break; p += len;
    if (type === 0x02) meta.libname = str(o, n);
    else if (type === 0x03) { uu = gdsRealRead(b, o); dbu = gdsRealRead(b, o + 8); }
    else if (type === 0x06) { cell = { polys: [], refs: [] }; cells.set(str(o, n), cell); }
    else if (type === 0x07) cell = null;
    else if ([0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x2D, 0x15].includes(type)) el = { type, layer: 0, datatype: 0, xy: [], width: 0, ptype: 0, mag: 1, angle: 0, flip: false, cols: 1, rows: 1, name: '', text: '' };
    else if (!el) continue;
    else if (type === 0x0D) el.layer = v.getInt16(o, false);
    else if (type === 0x0E || type === 0x2E || type === 0x16) el.datatype = v.getInt16(o, false);
    else if (type === 0x0F) el.width = v.getInt32(o, false);
    else if (type === 0x21) el.ptype = v.getInt16(o, false);
    else if (type === 0x10) for (let i = 0; i + 8 <= n; i += 8) el.xy.push([v.getInt32(o + i, false), v.getInt32(o + i + 4, false)]);
    else if (type === 0x12) el.name = str(o, n);
    else if (type === 0x13) { el.cols = v.getInt16(o, false); el.rows = v.getInt16(o + 2, false); }
    else if (type === 0x1A) el.flip = !!(v.getUint16(o, false) & 0x8000);
    else if (type === 0x1B) el.mag = gdsRealRead(b, o);
    else if (type === 0x1C) el.angle = gdsRealRead(b, o);
    else if (type === 0x19) el.text = str(o, n);
    else if (type === 0x11) {
      if (cell) {
        if ((el.type === 0x08 || el.type === 0x2D) && el.xy.length >= 4) cell.polys.push({ layer: el.layer, datatype: el.datatype, pts: el.xy.slice(0, -1) });
        else if (el.type === 0x09 && el.xy.length >= 2) { npath++; if (el.ptype === 1) roundEnds++; const wd = abs(el.width); if (wd > 0) { const pg = pathToPolygon(el.xy, wd, { extend: el.ptype === 2 || el.ptype === 1 ? wd / 2 : 0 }); if (pg.length >= 3) cell.polys.push({ layer: el.layer, datatype: el.datatype, pts: pg }); } }
        else if (el.type === 0x0A && el.xy.length) { nref++; cell.refs.push({ cell: el.name, x: el.xy[0][0], y: el.xy[0][1], angle: el.angle, mag: el.mag, flip: el.flip, offs: [[0, 0]] }); }
        else if (el.type === 0x0B && el.xy.length >= 3) { nref++; const [p0, p1, p2] = el.xy, c = max(1, el.cols), r = max(1, el.rows), offs = []; if (c * r > 200000) bigRef = true; for (let j = 0; j < r && offs.length < 200000; j++) for (let i = 0; i < c && offs.length < 200000; i++) offs.push([i * (p1[0] - p0[0]) / c + j * (p2[0] - p0[0]) / r, i * (p1[1] - p0[1]) / c + j * (p2[1] - p0[1]) / r]); cell.refs.push({ cell: el.name, x: p0[0], y: p0[1], angle: el.angle, mag: el.mag, flip: el.flip, offs }); }
        else if (el.type === 0x0C && el.xy.length) labels.push({ text: el.text, layer: el.layer, x: el.xy[0][0] * dbu * 1e9, y: el.xy[0][1] * dbu * 1e9 });
      }
      el = null;
    } else if (type === 0x04) break;
  }
  if (!cells.size) fail('GDSII', 'no structures (BGNSTR/STRNAME) found');
  const flat = flattenCells(cells, dbu * 1e9); // AREF lattice vectors are already expressed in parent coordinates
  Object.assign(meta, { userUnit_m: dbu / uu, dbUnit_m: dbu, cells: [...cells.keys()], topCells: flat.tops, layers: flat.layers, references: nref, labels: labels.slice(0, 200) });
  const notes = [`database unit ${N9(dbu * 1e9)} nm; coordinates converted to nm`, `hierarchy of ${cells.size} structures flattened from top cell(s) ${flat.tops.join(', ') || '(cyclic)'}`], lost = ['cell hierarchy is flattened (structure names kept in meta); properties and node/box types not imported'];
  if (npath) notes.push(`${npath} PATH elements converted to outline polygons`); if (roundEnds) lost.push(`${roundEnds} round-ended paths drawn with square extended ends`);
  if (flat.truncated || bigRef) lost.push('output limited to 200 000 polygons / array instances: the layout is larger than shown'); if (flat.missing.length) lost.push(`references to undefined structures ignored: ${flat.missing.slice(0, 8).join(', ')}`);
  if (labels.length) notes.push(`${labels.length} text labels listed in meta.labels`);
  return { kind: 'layout', polygons: flat.polygons, units: 'nm', meta, notes, lost };
}

// ---------------------------------------------------------------- OASIS (SEMI P39)
const OCT = [[1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]];
// CTRAPEZOID shapes (SEMI P39 compact trapezoids): vertex = (a·w + b·h, c·w + d·h) for each [a, b, c, d].
const CTRAP = ['0000 0001 1-101 1000', '0000 0001 1001 1-100', '0000 0101 1001 1000', '0100 0001 1001 1000', '0000 0101 1-101 1000', '0100 0001 1001 1-100', '0000 0101 1001 1-100', '0100 0001 1-101 1000',
  '0000 0001 10-11 1000', '0000 00-11 1001 1000', '0000 0001 1001 1010', '0010 0001 1001 1000', '0000 0001 10-11 1010', '0010 00-11 1001 1000', '0000 00-11 1001 1010', '0010 0001 10-11 1000',
  '0000 0010 1000', '0000 0010 1010', '0000 1010 1000', '0010 1010 1000', '0000 0101 0200', '0001 0201 0100', '0000 0020 1010', '1000 0010 1020', '0000 0001 1001 1000', '0000 0010 1010 1000'].map(r => r.split(' ').map(q => (q.match(/-?\d/g) || []).map(Number)));
// Minimal, strict-mode-free OASIS: START (table offsets in START, all zero), one or two CELL records by name,
// RECTANGLE / POLYGON (type-4 point lists) / PLACEMENT with a regular repetition, END padded to 256 bytes, no validation.
export function toOASIS(layout, { dbu_nm = 1, name = 'TOP', hierarchical = true } = {}) {
  if (!layout || layout.kind !== 'layout') fail('OASIS', `needs a layout geometry (got ${layout && layout.kind})`);
  const k = layoutScale(layout, dbu_nm * 1e-9), w = new Writer(1 << 16);
  const uint = x => { x = max(0, round(x)); while (x >= 128) { w.u8((x % 128) | 128); x = floor(x / 128); } w.u8(x); }, sint = x => uint(abs(x) * 2 + (x < 0 ? 1 : 0)), str = s => { const e = ENC.encode(s); uint(e.length); w.bytes(e); };
  const gdelta = (dx, dy) => { uint(abs(dx) * 4 + (dx < 0 ? 2 : 0) + 1); sint(dy); };
  const nstr = s => String(s).replace(/[^\x21-\x7E]/g, '_') || 'CELL';
  const shape = q => {
    const P = q.pts, xs = P.map(p => p[0]), ys = P.map(p => p[1]), x0 = Math.min(...xs), y0 = Math.min(...ys), x1 = Math.max(...xs), y1 = Math.max(...ys);
    if (P.length === 4 && P.every(p => (p[0] === x0 || p[0] === x1) && (p[1] === y0 || p[1] === y1)) && x1 > x0 && y1 > y0) { w.u8(20).u8(0x7B); uint(q.layer); uint(q.datatype); uint(x1 - x0); uint(y1 - y0); sint(x0); sint(y0); return; }
    w.u8(21).u8(0x3B); uint(q.layer); uint(q.datatype); uint(4); uint(P.length - 1); for (let i = 1; i < P.length; i++) gdelta(P[i][0] - P[i - 1][0], P[i][1] - P[i - 1][1]); sint(P[0][0]); sint(P[0][1]);
  };
  w.str('%SEMI-OASIS\r\n'); w.u8(1); str('1.0'); const gpm = 1000 / dbu_nm; if (Number.isInteger(gpm)) { uint(0); uint(gpm); } else { uint(7); w.f64(gpm); } uint(0); for (let i = 0; i < 12; i++) uint(0);
  const h = hierarchical && layout.hierarchy;
  if (h && h.nx >= 2 && h.ny >= 2) {
    w.u8(14); str('UNIT'); gridPolys({ polygons: h.unit, units: layout.units }, k).forEach(shape);
    w.u8(14); str(nstr(name)); w.u8(17).u8(0xB8); str('UNIT'); sint(round(h.origin[0] * k)); sint(round(h.origin[1] * k)); uint(1); uint(h.nx - 2); uint(h.ny - 2); uint(round(h.px * k)); uint(round(h.py * k));
  } else { w.u8(14); str(nstr(name)); gridPolys(layout, k).forEach(shape); }
  w.u8(2); uint(252); w.zero(252); uint(0);
  return w.out();
}
export function parseOASIS(input) {
  let b = toBytes(input); if (latin(b, 0, 13) !== '%SEMI-OASIS\r\n') fail('OASIS', 'missing "%SEMI-OASIS" magic');
  let p = 13, unit = 1000, version = '', cell = null, stop = '', ctrap = 0, cblocks = 0, texts = 0, big = false;
  const eof = () => fail('OASIS', 'unexpected end of file inside a record');
  const uint = () => { let x = 0, m = 1; for (; ;) { if (p >= b.length) eof(); const c = b[p++]; x += (c & 127) * m; if (!(c & 128)) return x; m *= 128; } };
  const sint = () => { const u = uint(); return u % 2 ? -floor(u / 2) : floor(u / 2); };
  const real = (t = uint()) => { if (t === 0) return uint(); if (t === 1) return -uint(); if (t === 2) return 1 / uint(); if (t === 3) return -1 / uint(); if (t === 4 || t === 5) { const a = uint(), c = uint(); return (t === 5 ? -1 : 1) * a / c; } if (t === 6) { const x = view(b).getFloat32(p, true); p += 4; return x; } if (t === 7) { const x = view(b).getFloat64(p, true); p += 8; return x; } fail('OASIS', `invalid real type ${t}`); };
  const str = () => { const n = uint(), s = latin(b, p, n); p += n; return s; };
  const d2 = () => { const u = uint(), m = floor(u / 4), d = OCT[u % 4]; return [d[0] * m, d[1] * m]; }, d3 = () => { const u = uint(), m = floor(u / 8), d = OCT[u % 8]; return [d[0] * m, d[1] * m]; };
  const dg = () => { const u = uint(); if (u % 2 === 0) { const d = OCT[floor(u / 2) % 8], m = floor(u / 16); return [d[0] * m, d[1] * m]; } const dx = floor(u / 4) * (floor(u / 2) % 2 ? -1 : 1); return [dx, sint()]; };
  const plist = poly => { const t = uint(), n = uint(), pts = [[0, 0]]; let x = 0, y = 0;
    if (t < 2) { let hz = t === 0; for (let i = 0; i < n; i++) { const d = sint(); if (hz) x += d; else y += d; pts.push([x, y]); hz = !hz; } if (poly) pts.push(hz ? [0, y] : [x, 0]); }
    else if (t < 5) for (let i = 0; i < n; i++) { const d = t === 2 ? d2() : t === 3 ? d3() : dg(); x += d[0]; y += d[1]; pts.push([x, y]); }
    else if (t === 5) { let ax = 0, ay = 0; for (let i = 0; i < n; i++) { const d = dg(); ax += d[0]; ay += d[1]; x += ax; y += ay; pts.push([x, y]); } } else fail('OASIS', `invalid point-list type ${t}`);
    return pts; };
  const M = { abs: true, px: 0, py: 0, gx: 0, gy: 0, tx: 0, ty: 0, layer: 0, dt: 0, w: 0, h: 0, pl: null, hw: 0, se: 0, ee: 0, ppl: null, rep: [[0, 0]], pcell: null, rad: 0, ct: 0 };
  const repetition = () => { const t = uint(); if (t === 0) return M.rep; let o = [[0, 0]]; const lim = 2e5;
    const line = (n, f) => { let x = 0, y = 0; for (let i = 1; i < n; i++) { const d = f(); x += d[0]; y += d[1]; if (o.length < lim) o.push([x, y]); else big = true; } };
    if (t === 1 || t === 2 || t === 3) { const nx = t === 3 ? 1 : uint() + 2, ny = t === 2 ? 1 : uint() + 2, sx = t === 3 ? 0 : uint(), sy = t === 2 ? 0 : uint(); o = []; if (nx * ny > lim) big = true; for (let j = 0; j < ny && o.length < lim; j++) for (let i = 0; i < nx && o.length < lim; i++) o.push([i * sx, j * sy]); }
    else if (t === 4 || t === 5) { const n = uint() + 2, g = t === 5 ? uint() : 1; line(n, () => [uint() * g, 0]); } else if (t === 6 || t === 7) { const n = uint() + 2, g = t === 7 ? uint() : 1; line(n, () => [0, uint() * g]); }
    else if (t === 8) { const n = uint() + 2, m = uint() + 2, a = dg(), c = dg(); o = []; if (n * m > lim) big = true; for (let j = 0; j < m && o.length < lim; j++) for (let i = 0; i < n && o.length < lim; i++) o.push([i * a[0] + j * c[0], i * a[1] + j * c[1]]); }
    else if (t === 9) { const n = uint() + 2, a = dg(); o = []; if (n > lim) big = true; for (let i = 0; i < n && i < lim; i++) o.push([i * a[0], i * a[1]]); }
    else if (t === 10 || t === 11) { const n = uint() + 2, g = t === 11 ? uint() : 1; line(n, () => { const d = dg(); return [d[0] * g, d[1] * g]; }); } else fail('OASIS', `invalid repetition type ${t}`);
    return (M.rep = o); };
  const xy = (info, kx, ky) => { if (info & 16) { const x = sint(); M[kx] = M.abs ? x : M[kx] + x; } if (info & 8) { const y = sint(); M[ky] = M.abs ? y : M[ky] + y; } };
  const ld = info => { if (info & 1) M.layer = uint(); if (info & 2) M.dt = uint(); };
  const put = (pts, offs) => { for (const o of offs) cell.polys.push({ layer: M.layer, datatype: M.dt, pts: pts.map(q => [q[0] + M.gx + o[0], q[1] + M.gy + o[1]]) }); };
  const cells = new Map(), names = new Map(); let autoName = 0;
  const newCell = key => { cell = { polys: [], refs: [] }; cells.set(key, cell); Object.assign(M, { abs: true, px: 0, py: 0, gx: 0, gy: 0, tx: 0, ty: 0 }); };
  const need = () => { if (!cell) fail('OASIS', 'element record outside a CELL'); };
  try {
    for (; p < b.length;) {
      const at = p, id = uint();
      if (id === 0) continue;
      if (id === 1) { version = str(); unit = real(); if (uint() === 0) for (let i = 0; i < 12; i++) uint(); }
      else if (id === 2) break;
      else if (id === 3) names.set(autoName++, str()); else if (id === 4) { const s = str(); names.set(uint(), s); }
      else if (id === 5 || id === 7 || id === 9) str(); else if (id === 6 || id === 8 || id === 10) { str(); uint(); }
      else if (id === 11 || id === 12) { str(); for (let k = 0; k < 2; k++) { const t = uint(); if (t === 4) { uint(); uint(); } else if (t > 0) uint(); } }
      else if (id === 13) newCell(uint()); else if (id === 14) newCell(str());
      else if (id === 15) M.abs = true; else if (id === 16) M.abs = false;
      else if (id === 17 || id === 18) { need(); const info = b[p++]; if (info & 128) M.pcell = info & 64 ? uint() : str(); let mag = 1, ang = (info >> 1 & 3) * 90; if (id === 18) { if (info & 4) mag = real(); ang = info & 2 ? real() : 0; } if (info & 32) { const x = sint(); M.px = M.abs ? x : M.px + x; } if (info & 16) { const y = sint(); M.py = M.abs ? y : M.py + y; } const offs = info & 8 ? repetition() : [[0, 0]]; cell.refs.push({ cell: M.pcell, x: M.px, y: M.py, angle: ang, mag, flip: !!(info & 1), offs }); }
      else if (id === 19) { need(); const info = b[p++]; if (info & 64) { if (info & 32) uint(); else str(); } if (info & 1) uint(); if (info & 2) uint(); xy(info, 'tx', 'ty'); if (info & 4) repetition(); texts++; }
      else if (id === 20) { need(); const info = b[p++]; ld(info); if (info & 64) M.w = uint(); if (info & 32) M.h = uint(); if (info & 128) M.h = M.w; xy(info, 'gx', 'gy'); put([[0, 0], [M.w, 0], [M.w, M.h], [0, M.h]], info & 4 ? repetition() : [[0, 0]]); }
      else if (id === 21) { need(); const info = b[p++]; ld(info); if (info & 32) M.pl = plist(true); xy(info, 'gx', 'gy'); const offs = info & 4 ? repetition() : [[0, 0]]; if (!M.pl) fail('OASIS', 'POLYGON reuses a point list that was never defined'); put(M.pl, offs); }
      else if (id === 22) { need(); const info = b[p++]; ld(info); if (info & 64) M.hw = uint(); if (info & 128) { const e = uint(), ss = e >> 2 & 3, ee = e & 3; if (ss === 1) M.se = 0; else if (ss === 2) M.se = M.hw; else if (ss === 3) M.se = sint(); if (ee === 1) M.ee = 0; else if (ee === 2) M.ee = M.hw; else if (ee === 3) M.ee = sint(); } if (info & 32) M.ppl = plist(false); xy(info, 'gx', 'gy'); const offs = info & 4 ? repetition() : [[0, 0]]; if (M.ppl && M.hw > 0) { const pg = pathToPolygon(M.ppl, 2 * M.hw, { extend: (M.se + M.ee) / 2 }); if (pg.length >= 3) put(pg, offs); } }
      else if (id >= 23 && id <= 25) { need(); const info = b[p++]; ld(info); if (info & 64) M.w = uint(); if (info & 32) M.h = uint(); const da = id !== 25 ? sint() : 0, db = id !== 24 ? sint() : 0; xy(info, 'gx', 'gy'); const offs = info & 4 ? repetition() : [[0, 0]], w = M.w, h = M.h;
        put(info & 128 ? [[0, max(da, 0)], [0, h + min(db, 0)], [w, h - max(db, 0)], [w, -min(da, 0)]] : [[max(da, 0), h], [w + min(db, 0), h], [w - max(db, 0), 0], [-min(da, 0), 0]], offs); }
      else if (id === 26) { need(); const info = b[p++]; ld(info); if (info & 128) M.ct = uint(); if (info & 64) M.w = uint(); if (info & 32) M.h = uint(); xy(info, 'gx', 'gy'); const offs = info & 4 ? repetition() : [[0, 0]], tb = CTRAP[M.ct]; if (!tb) { ctrap++; continue; }
        const w = M.ct === 20 || M.ct === 21 ? 2 * M.h : M.w, h = (M.ct >= 16 && M.ct <= 19) || M.ct === 25 ? w : M.ct === 22 || M.ct === 23 ? 2 * w : M.h; put(tb.map(q => [q[0] * w + q[1] * h, q[2] * w + q[3] * h]), offs); }
      else if (id === 27) { need(); const info = b[p++]; ld(info); if (info & 32) M.rad = uint(); xy(info, 'gx', 'gy'); const offs = info & 4 ? repetition() : [[0, 0]]; put(Array.from({ length: 32 }, (_, i) => [M.rad * cos(PI * i / 16), M.rad * sin(PI * i / 16)]), offs); }
      else if (id === 28) { const info = b[p++]; if (info & 4) { if (info & 2) uint(); else str(); } if (!(info & 8)) { let n = info >> 4; if (n === 15) n = uint(); for (let i = 0; i < n; i++) { const t = uint(); if (t <= 7) real(t); else if (t === 8) uint(); else if (t === 9) sint(); else if (t <= 12) str(); else uint(); } } }
      else if (id === 29) { /* PROPERTY repeat */ }
      else if (id === 30) { uint(); str(); } else if (id === 31) { uint(); str(); uint(); } else if (id === 32) { uint(); str(); }
      else if (id === 33) { need(); const info = b[p++]; uint(); ld(info); str(); xy(info, 'gx', 'gy'); if (info & 4) repetition(); }
      else if (id === 34) { const ct = uint(), un = uint(), cn = uint(); if (ct !== 0) fail('OASIS', `CBLOCK compression type ${ct} is not supported`); const inf = inflateRaw(b, p, un).data, nb = new Uint8Array(at + inf.length + b.length - p - cn); nb.set(b.subarray(0, at)); nb.set(inf, at); nb.set(b.subarray(p + cn), at + inf.length); b = nb; p = at; cblocks++; }
      else { stop = `unknown record id ${id} at byte ${at}`; break; }
    }
  } catch (e) { if (!cells.size) throw e; stop = e.message; }
  if (!cells.size) fail('OASIS', 'no CELL records found');
  const res = k => typeof k === 'number' ? (names.get(k) ?? `cell#${k}`) : k, named = new Map([...cells].map(([k, c]) => [res(k), { polys: c.polys, refs: c.refs.map(r => ({ ...r, cell: res(r.cell) })) }]));
  const nmPerGrid = 1000 / unit, flat = flattenCells(named, nmPerGrid), notes = [`grid ${N9(unit)} steps per µm (${N9(nmPerGrid)} nm); coordinates converted to nm`, `hierarchy of ${cells.size} cells flattened from top cell(s) ${flat.tops.join(', ') || '(cyclic)'}`], lost = ['cell hierarchy is flattened; properties, text and layer-name tables not imported'];
  if (cblocks) notes.push(`${cblocks} compressed CBLOCK records inflated`); if (ctrap) lost.push(`${ctrap} CTRAPEZOID records with an invalid type code were skipped`); if (texts) lost.push(`${texts} TEXT records not imported`);
  if (stop) lost.push(`parsing stopped early (${stop}): shapes after that point are missing`); if (big || flat.truncated) lost.push('very large repetitions or more than 200 000 polygons: output truncated');
  return { kind: 'layout', polygons: flat.polygons, units: 'nm', complete: !stop && !ctrap && !flat.truncated, meta: { version, gridPerMicron: unit, cells: [...named.keys()], topCells: flat.tops, layers: flat.layers }, notes, lost };
}

// ---------------------------------------------------------------- DXF
const INSUNITS = { 1: 'in', 2: 'ft', 4: 'mm', 5: 'cm', 6: 'm', 8: 'µin', 9: 'mil', 11: 'Å', 12: 'nm', 13: 'µm' };
// Minimal DXF (HEADER + ENTITIES). Layouts → closed LWPOLYLINEs (or R12 POLYLINEs with r12:true), lattices → LINEs.
// The file has no TABLES/BLOCKS/OBJECTS sections: layout and CAM tools read it; AutoCAD itself is strict about
// post-R12 files and may ask for recovery — use r12:true for it.
export function toDXF(g, { r12 = false } = {}) {
  const code = Object.entries(INSUNITS).find(([, u]) => u === g.units || (g.units === 'um' && u === 'µm'))?.[0] || 0, o = ['0', 'SECTION', '2', 'HEADER', '9', '$ACADVER', '1', r12 ? 'AC1009' : 'AC1015', '9', '$INSUNITS', '70', String(code), '0', 'ENDSEC', '0', 'SECTION', '2', 'ENTITIES'];
  let hnd = 0x100; const ent = (type, layer, sub) => { o.push('0', type); if (!r12) o.push('5', (hnd++).toString(16).toUpperCase(), '100', 'AcDbEntity'); o.push('8', layer); if (!r12 && sub) o.push('100', sub); };
  if (g.kind === 'layout') for (const q of g.polygons) {
    if (r12) { ent('POLYLINE', 'L' + (q.layer ?? 1)); o.push('66', '1', '70', '1'); for (const p of q.pts) { ent('VERTEX', 'L' + (q.layer ?? 1)); o.push('10', N9(p[0]), '20', N9(p[1]), '30', '0'); } ent('SEQEND', 'L' + (q.layer ?? 1)); }
    else { ent('LWPOLYLINE', 'L' + (q.layer ?? 1), 'AcDbPolyline'); o.push('90', String(q.pts.length), '70', '1'); for (const p of q.pts) o.push('10', N9(p[0]), '20', N9(p[1])); }
  } else if (g.kind === 'lattice') for (const e of g.edges) { const a = g.nodes[e[0]], b = g.nodes[e[1]]; ent('LINE', 'STRUTS', 'AcDbLine'); o.push('10', N9(a[0]), '20', N9(a[1]), '30', N9(a[2] || 0), '11', N9(b[0]), '21', N9(b[1]), '31', N9(b[2] || 0)); }
  else fail('DXF', `needs a layout or lattice geometry (got ${g.kind})`);
  o.push('0', 'ENDSEC', '0', 'EOF'); return o.join('\n') + '\n';
}
const arcPts = (cx, cy, r, a0, a1, close) => { const n = max(2, Math.ceil(abs(a1 - a0) / (PI / 24))), o = []; for (let k = 0; k <= n - (close ? 1 : 0); k++) { const t = a0 + (a1 - a0) * k / n; o.push([cx + r * cos(t), cy + r * sin(t)]); } return o; };
function bulge(p, q, bg) { const th = 4 * Math.atan(bg), d = hypot(q[0] - p[0], q[1] - p[1]); if (!(d > 0) || abs(th) < 1e-9) return []; const r = d / (2 * sin(th / 2)), hgt = r * cos(th / 2), cx = (p[0] + q[0]) / 2 - (q[1] - p[1]) / d * hgt, cy = (p[1] + q[1]) / 2 + (q[0] - p[0]) / d * hgt, a0 = Math.atan2(p[1] - cy, p[0] - cx); return arcPts(cx, cy, abs(r), a0, a0 + th, false).slice(1, -1); }
export function parseDXF(input) {
  const b = toBytes(input); if (latin(b, 0, 18) === 'AutoCAD Binary DXF') fail('DXF', 'binary DXF is not supported: save as ASCII DXF');
  const L = lines(toText(b)), pairs = []; for (let i = 0; i + 1 < L.length; i += 2) { const c = parseInt(L[i], 10); if (!Number.isFinite(c)) { if (L[i].trim() === '' && i + 2 >= L.length) break; fail('DXF', `line ${i + 1}: expected a numeric group code, found "${L[i].trim().slice(0, 30)}"`); } pairs.push([c, L[i + 1].trim()]); }
  if (!pairs.some(q => q[0] === 0 && q[1] === 'SECTION') && !pairs.some(q => q[0] === 0 && q[1] === 'EOF')) fail('DXF', 'no SECTION / EOF markers found');
  const meta = { version: '', units: '' }, layerIds = new Map(), polygons = [], nodes = [], edges = [], V = [], F = [], other = {}, points = []; let sec = '', i = 0, unitCode = 0, inserts = 0, bulges = 0, openPaths = 0;
  const lid = n => { if (!layerIds.has(n)) layerIds.set(n, layerIds.size + 1); return layerIds.get(n); };
  const nkey = new Map(), nadd = q => { const k = N9(q[0]) + ',' + N9(q[1]) + ',' + N9(q[2] || 0); let i = nkey.get(k); if (i === undefined) { i = nodes.length; nodes.push([q[0], q[1], q[2] || 0]); nkey.set(k, i); } return i; };
  const seg = pts => { let prev = -1; for (const q of pts) { const c = nadd(q); if (prev >= 0 && prev !== c) edges.push([prev, c]); prev = c; } }; // coincident end points are merged so that LINE soups become a connected network
  const poly = (verts, closed, layer) => { const pts = []; verts.forEach((v, k) => { pts.push([v.x, v.y]); const nx = verts[(k + 1) % verts.length]; if (v.b && (closed || k < verts.length - 1)) { bulges++; pts.push(...bulge([v.x, v.y], [nx.x, nx.y], v.b)); } }); if (pts.length > 2 && hypot(pts[0][0] - pts[pts.length - 1][0], pts[0][1] - pts[pts.length - 1][1]) < 1e-12) { pts.pop(); closed = true; } if (closed && pts.length >= 3) polygons.push({ layer: lid(layer), pts: polyArea(pts) < 0 ? pts.reverse() : pts }); else if (pts.length >= 2) { openPaths++; seg(pts); } };
  while (i < pairs.length) {
    const [c, v] = pairs[i];
    if (c === 0 && v === 'SECTION') { sec = pairs[i + 1] ? pairs[i + 1][1] : ''; i += 2; continue; }
    if (sec === 'HEADER') { if (c === 9 && v === '$ACADVER') meta.version = pairs[i + 1][1]; if (c === 9 && v === '$INSUNITS') unitCode = +pairs[i + 1][1]; if (c === 9 && v === '$MEASUREMENT' && !unitCode) meta.measurement = +pairs[i + 1][1] ? 'metric' : 'imperial'; i++; continue; }
    if (sec === 'BLOCKS') { if (c === 0 && v === 'BLOCK') other.BLOCK = (other.BLOCK || 0) + 1; i++; continue; }
    if (sec !== 'ENTITIES' || c !== 0) { i++; continue; }
    const type = v, e = {}; let j = i + 1; const multi = { 10: [], 20: [], 30: [], 42: [] };
    for (; j < pairs.length && pairs[j][0] !== 0; j++) { const [k, val] = pairs[j]; if (k in multi) { multi[k].push(+val); if (k === 10) multi[42].length = multi[10].length - 1; } if (!(k in e)) e[k] = val; }
    const layer = e[8] || '0', f = k => +e[k] || 0;
    if (type === 'LINE') seg([[f(10), f(20), f(30)], [f(11), f(21), f(31)]]);
    else if (type === 'LWPOLYLINE') { const bl = pairs.slice(i + 1, j); let vs = [], curv = null; for (const [k, val] of bl) { if (k === 10) { curv = { x: +val, y: 0, b: 0 }; vs.push(curv); } else if (k === 20 && curv) curv.y = +val; else if (k === 42 && curv) curv.b = +val; } poly(vs, !!(f(70) & 1), layer); }
    else if (type === 'POLYLINE') { const closed = !!(f(70) & 1), mesh3 = !!(f(70) & 64), vs = []; let k = j; while (k < pairs.length && !(pairs[k][0] === 0 && pairs[k][1] === 'SEQEND')) { if (pairs[k][0] === 0 && pairs[k][1] === 'VERTEX') { const vv = { x: 0, y: 0, z: 0, b: 0 }; k++; for (; k < pairs.length && pairs[k][0] !== 0; k++) { const [cc, val] = pairs[k]; if (cc === 10) vv.x = +val; else if (cc === 20) vv.y = +val; else if (cc === 30) vv.z = +val; else if (cc === 42) vv.b = +val; } vs.push(vv); } else k++; } j = k + 1; if (mesh3) other['POLYLINE (polyface mesh)'] = (other['POLYLINE (polyface mesh)'] || 0) + 1; else poly(vs, closed, layer); }
    else if (type === 'CIRCLE') polygons.push({ layer: lid(layer), pts: arcPts(f(10), f(20), f(40), 0, 2 * PI, true) });
    else if (type === 'ARC') { let a0 = f(50) * PI / 180, a1 = f(51) * PI / 180; if (a1 <= a0) a1 += 2 * PI; seg(arcPts(f(10), f(20), f(40), a0, a1, false).map(q => [q[0], q[1], f(30)])); }
    else if (type === 'ELLIPSE') { const mx = f(11), my = f(21), a = hypot(mx, my), bb = a * (f(40) || 1), rot = Math.atan2(my, mx), t0 = f(41), t1 = e[42] !== undefined ? +e[42] : 2 * PI, full = abs(t1 - t0 - 2 * PI) < 1e-6, n = 64, pts = []; for (let k = 0; k <= n - (full ? 1 : 0); k++) { const t = t0 + (t1 - t0) * k / n, x = a * cos(t), y = bb * sin(t); pts.push([f(10) + x * cos(rot) - y * sin(rot), f(20) + x * sin(rot) + y * cos(rot)]); } if (full) polygons.push({ layer: lid(layer), pts }); else seg(pts); }
    else if (type === '3DFACE') { const q = [[f(10), f(20), f(30)], [f(11), f(21), f(31)], [f(12), f(22), f(32)], [f(13), f(23), f(33)]], base = V.length; V.push(...q); F.push([base, base + 1, base + 2]); if (hypot(q[3][0] - q[2][0], q[3][1] - q[2][1], q[3][2] - q[2][2]) > 1e-12) F.push([base, base + 2, base + 3]); }
    else if (type === 'SOLID' || type === 'TRACE') polygons.push({ layer: lid(layer), pts: [[f(10), f(20)], [f(11), f(21)], [f(13), f(23)], [f(12), f(22)]] });
    else if (type === 'POINT') points.push([f(10), f(20), f(30)]);
    else if (type === 'INSERT') inserts++;
    else if (!['ENDSEC', 'EOF', 'SEQEND', 'VERTEX'].includes(type)) other[type] = (other[type] || 0) + 1;
    i = j;
  }
  meta.units = INSUNITS[unitCode] || ''; meta.insunits = unitCode; meta.layers = Object.fromEntries(layerIds); meta.otherEntities = other;
  const notes = [meta.units ? `$INSUNITS = ${unitCode} (${meta.units})` : 'no $INSUNITS in the header: the drawing is unitless and the length unit must be confirmed'], lost = [];
  if (bulges) notes.push(`${bulges} polyline arc segments (bulges) discretised at 7.5° steps`); if (openPaths) notes.push(`${openPaths} open polylines kept as wire segments`);
  if (inserts) lost.push(`${inserts} block references (INSERT) not expanded`); const oth = Object.entries(other).filter(([k]) => k !== 'BLOCK'); if (oth.length) lost.push(`entities not converted: ${oth.map(([k, n]) => `${k} × ${n}`).join(', ')}`);
  if (F.length) { const w = weld(F.flatMap(f => f.flatMap(k => V[k]))); if (polygons.length || edges.length) lost.push(`${polygons.length} closed profiles and ${edges.length} wire segments ignored in favour of the 3DFACE mesh`); return { kind: 'mesh', vertices: w.V, faces: w.F, units: meta.units, meta, notes, lost }; }
  if (polygons.length) { const g = { kind: 'layout', polygons, units: meta.units || 'unitless', layerNames: Object.fromEntries([...layerIds].map(([n, k]) => [k, n])), meta, notes, lost }; if (edges.length) g.wire = { nodes, edges }; return g; }
  if (edges.length) return { kind: 'lattice', nodes, edges, radius: 0, units: meta.units || 'unitless', wireframe: true, meta, notes, lost };
  if (points.length) return { kind: 'lattice', nodes: points, edges: [], radius: 0, units: meta.units || 'unitless', pointCloud: true, meta, notes, lost };
  fail('DXF', `no convertible entities in the ENTITIES section${oth.length ? ` (found only ${oth.map(q => q[0]).join(', ')})` : ''}`);
}

// ---------------------------------------------------------------- SVG
export function toSVG(layout, { fill } = {}) {
  if (!layout || layout.kind !== 'layout') fail('SVG', `needs a layout geometry (got ${layout && layout.kind})`);
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity; for (const q of layout.polygons) for (const p of q.pts) { x0 = min(x0, p[0]); y0 = min(y0, p[1]); x1 = max(x1, p[0]); y1 = max(y1, p[1]); }
  if (!layout.polygons.length) { x0 = y0 = 0; x1 = y1 = 1; }
  const cols = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#17becf', '#8c564b', '#7f7f7f'], byL = new Map(); for (const q of layout.polygons) { const k = q.layer ?? 1; if (!byL.has(k)) byL.set(k, []); byL.get(k).push(q); }
  // y is flipped (SVG y points down) so the pattern is displayed as designed
  const body = [...byL].map(([k, ps], i) => `<g id="layer${k}" data-layer="${k}" fill="${fill || cols[i % cols.length]}">\n${ps.map(q => `<polygon points="${q.pts.map(p => `${N9(p[0] - x0)},${N9(y1 - p[1])}`).join(' ')}"/>`).join('\n')}\n</g>`).join('\n');
  return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${N9(x1 - x0)} ${N9(y1 - y0)}" data-units="${esc(layout.units || 'nm')}" data-origin="${N9(x0)} ${N9(y0)}">\n<title>layout (${esc(layout.units || 'nm')} per user unit)</title>\n${body}\n</svg>\n`;
}
const svgT = s => { let M = [1, 0, 0, 1, 0, 0]; for (const m of (s || '').matchAll(/(matrix|translate|scale|rotate|skewX|skewY)\s*\(([^)]*)\)/g)) { const a = nums(m[2]); let T; if (m[1] === 'matrix') T = a.slice(0, 6); else if (m[1] === 'translate') T = [1, 0, 0, 1, a[0] || 0, a[1] || 0]; else if (m[1] === 'scale') T = [a[0], 0, 0, a[1] ?? a[0], 0, 0]; else if (m[1] === 'rotate') { const c = cos(a[0] * PI / 180), sn = sin(a[0] * PI / 180), cx = a[1] || 0, cy = a[2] || 0; T = [c, sn, -sn, c, cx - c * cx + sn * cy, cy - sn * cx - c * cy]; } else if (m[1] === 'skewX') T = [1, 0, Math.tan(a[0] * PI / 180), 1, 0, 0]; else T = [1, Math.tan(a[0] * PI / 180), 0, 1, 0, 0]; M = affMul(M, T); } return M; };
function svgPath(d) { // → [{pts, closed}], curves flattened
  const tk = d.match(/[MmLlHhVvCcSsQqTtAaZz]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g) || [], subs = []; let i = 0, cmd = '', x = 0, y = 0, sx = 0, sy = 0, cur = null, px = 0, py = 0, prevC = '';
  const n = () => +tk[i++], more = () => i < tk.length && !/^[A-Za-z]$/.test(tk[i]), bez = (P, steps = 12) => { for (let k = 1; k <= steps; k++) { const t = k / steps, u = 1 - t; cur.pts.push(P.length === 4 ? [u * u * u * P[0][0] + 3 * u * u * t * P[1][0] + 3 * u * t * t * P[2][0] + t * t * t * P[3][0], u * u * u * P[0][1] + 3 * u * u * t * P[1][1] + 3 * u * t * t * P[2][1] + t * t * t * P[3][1]] : [u * u * P[0][0] + 2 * u * t * P[1][0] + t * t * P[2][0], u * u * P[0][1] + 2 * u * t * P[1][1] + t * t * P[2][1]]); } };
  while (i < tk.length) {
    if (/^[A-Za-z]$/.test(tk[i])) cmd = tk[i++]; else if (!cmd) break; const rel = cmd === cmd.toLowerCase(), C = cmd.toUpperCase();
    if (C === 'Z') { if (cur) cur.closed = true; x = sx; y = sy; cur = null; prevC = C; continue; }
    if (!more()) continue;
    if (C === 'M') { x = n() + (rel ? x : 0); y = n() + (rel ? y : 0); sx = x; sy = y; cur = { pts: [[x, y]], closed: false }; subs.push(cur); cmd = rel ? 'l' : 'L'; prevC = 'M'; continue; }
    if (!cur) { cur = { pts: [[x, y]], closed: false }; subs.push(cur); }
    if (C === 'L') { x = n() + (rel ? x : 0); y = n() + (rel ? y : 0); cur.pts.push([x, y]); }
    else if (C === 'H') { x = n() + (rel ? x : 0); cur.pts.push([x, y]); } else if (C === 'V') { y = n() + (rel ? y : 0); cur.pts.push([x, y]); }
    else if (C === 'C' || C === 'S') { const ox = rel ? x : 0, oy = rel ? y : 0; let c1 = C === 'C' ? [n() + ox, n() + oy] : (/[CS]/.test(prevC) ? [2 * x - px, 2 * y - py] : [x, y]); const c2 = [n() + ox, n() + oy], e = [n() + ox, n() + oy]; bez([[x, y], c1, c2, e]); px = c2[0]; py = c2[1]; x = e[0]; y = e[1]; }
    else if (C === 'Q' || C === 'T') { const ox = rel ? x : 0, oy = rel ? y : 0, c1 = C === 'Q' ? [n() + ox, n() + oy] : (/[QT]/.test(prevC) ? [2 * x - px, 2 * y - py] : [x, y]), e = [n() + ox, n() + oy]; bez([[x, y], c1, e]); px = c1[0]; py = c1[1]; x = e[0]; y = e[1]; }
    else if (C === 'A') { let rx = abs(n()), ry = abs(n()); const phi = n() * PI / 180, fa = n(), fs = n(), ex = n() + (rel ? x : 0), ey = n() + (rel ? y : 0);
      if (rx > 0 && ry > 0 && (ex !== x || ey !== y)) { const cp = cos(phi), sp = sin(phi), dx = (x - ex) / 2, dy = (y - ey) / 2, x1 = cp * dx + sp * dy, y1 = -sp * dx + cp * dy, lam = x1 * x1 / (rx * rx) + y1 * y1 / (ry * ry); if (lam > 1) { rx *= sqrt(lam); ry *= sqrt(lam); }
        const num0 = max(0, rx * rx * ry * ry - rx * rx * y1 * y1 - ry * ry * x1 * x1), co = (fa !== fs ? 1 : -1) * sqrt(num0 / (rx * rx * y1 * y1 + ry * ry * x1 * x1 || 1)), cxp = co * rx * y1 / ry, cyp = -co * ry * x1 / rx, cx = cp * cxp - sp * cyp + (x + ex) / 2, cy = sp * cxp + cp * cyp + (y + ey) / 2;
        const t0 = Math.atan2((y1 - cyp) / ry, (x1 - cxp) / rx); let dt = Math.atan2((-y1 - cyp) / ry, (-x1 - cxp) / rx) - t0; if (fs && dt < 0) dt += 2 * PI; if (!fs && dt > 0) dt -= 2 * PI; const st = max(2, Math.ceil(abs(dt) / (PI / 18)));
        for (let k = 1; k <= st; k++) { const t = t0 + dt * k / st; cur.pts.push([cx + rx * cos(t) * cp - ry * sin(t) * sp, cy + rx * cos(t) * sp + ry * sin(t) * cp]); } } else cur.pts.push([ex, ey]); x = ex; y = ey; }
    prevC = C;
  }
  return subs;
}
// String-based SVG reader (no DOM): path, polygon, polyline, rect, circle, ellipse and line elements with nested
// group transforms. Closed shapes become polygons; open strokes are kept as wire segments.
export function parseSVG(input) {
  let b = toBytes(input); if (b[0] === 0x1F && b[1] === 0x8B) b = gunzip(b);
  const t = toText(b), root = /<svg\b([^>]*)>/.exec(t); if (!root) fail('SVG', 'no <svg> root element');
  const vb = nums(attr(root[1], 'viewBox') || ''), wAttr = attr(root[1], 'width') || '', hAttr = attr(root[1], 'height') || '', wm = /^([\d.eE+-]+)\s*(mm|cm|in|pt|pc|px|um|µm|nm)?$/.exec(wAttr.trim());
  const polygons = [], nodes = [], edges = [], stack = [[1, 0, 0, 1, 0, 0]], lstack = [1], skipTags = /^(defs|clipPath|mask|pattern|symbol|marker|metadata|style|script)$/; let skip = 0, layerN = 0, uses = 0, textN = 0, images = 0, rounded = 0, openN = 0;
  const emit = (pts, closed, T, layer) => { const P = pts.map(p => affPt(T, p)); if (closed && P.length >= 3) polygons.push({ layer, pts: P }); else if (P.length >= 2) { openN++; let prev = -1; for (const p of P) { nodes.push([p[0], p[1], 0]); if (prev >= 0) edges.push([prev, nodes.length - 1]); prev = nodes.length - 1; } } };
  const re = /<(\/?)([A-Za-z][\w:.-]*)\b((?:"[^"]*"|'[^']*'|[^>"'])*?)(\/?)>/g; let m;
  while ((m = re.exec(t))) {
    const [, close, tag, at, self] = m, nm = tag.replace(/^svg:/, '');
    if (close) { if (skipTags.test(nm)) skip = max(0, skip - 1); else if (stack.length > 1) { stack.pop(); lstack.pop(); } continue; }
    if (skipTags.test(nm)) { if (!self) skip++; continue; } if (skip) { continue; }
    const T = affMul(stack[stack.length - 1], svgT(attr(at, 'transform'))); let layer = lstack[lstack.length - 1];
    if (nm === 'g' && /groupmode\s*=\s*"layer"/.test(at)) layer = ++layerN; else if (nm === 'g' && attr(at, 'data-layer') !== undefined) layer = +attr(at, 'data-layer') || layer;
    const f = k => parseFloat(attr(at, k)) || 0;
    if (nm === 'path') for (const s of svgPath(attr(at, 'd') || '')) emit(s.pts, s.closed, T, layer);
    else if (nm === 'polygon' || nm === 'polyline') { const a = nums(attr(at, 'points') || ''), pts = []; for (let k = 0; k + 1 < a.length; k += 2) pts.push([a[k], a[k + 1]]); emit(pts, nm === 'polygon', T, layer); }
    else if (nm === 'rect') { const x = f('x'), y = f('y'), w = f('width'), h = f('height'); if (attr(at, 'rx') || attr(at, 'ry')) rounded++; if (w > 0 && h > 0) emit([[x, y], [x + w, y], [x + w, y + h], [x, y + h]], true, T, layer); }
    else if (nm === 'circle' || nm === 'ellipse') { const rx = nm === 'circle' ? f('r') : f('rx'), ry = nm === 'circle' ? f('r') : f('ry'); if (rx > 0 && ry > 0) emit(Array.from({ length: 64 }, (_, k) => [f('cx') + rx * cos(PI * k / 32), f('cy') + ry * sin(PI * k / 32)]), true, T, layer); }
    else if (nm === 'line') emit([[f('x1'), f('y1')], [f('x2'), f('y2')]], false, T, layer);
    else if (nm === 'use') uses++; else if (nm === 'text') textN++; else if (nm === 'image') images++;
    if (!self) { stack.push(T); lstack.push(layer); }
  }
  if (!polygons.length && !edges.length) fail('SVG', 'no path, polygon, rect, circle, ellipse or line elements found' + (images ? ' (the file embeds raster images only)' : ''));
  // user units → physical units when width carries a unit and a viewBox exists; then flip y to point upwards
  let scale = 1, units = attr(root[1], 'data-units') || 'px', H = vb.length === 4 ? vb[3] : parseFloat(hAttr) || 0; const notes = [], lost = [];
  if (attr(root[1], 'data-units')) notes.push(`units "${units}" taken from the data-units attribute written by this application`);
  else if (wm && wm[2] && wm[2] !== 'px' && vb.length === 4 && vb[2] > 0) { scale = +wm[1] / vb[2]; units = wm[2] === 'um' ? 'µm' : wm[2]; if (units === 'pt') { scale *= 25.4 / 72; units = 'mm'; } if (units === 'pc') { scale *= 25.4 / 6; units = 'mm'; } notes.push(`physical size from width="${wAttr}" and the viewBox: 1 user unit = ${N9(scale)} ${units}`); }
  else notes.push('no physical size declared: coordinates are SVG user units (px); the scale must be set by the user');
  if (!H) for (const q of polygons) for (const p of q.pts) H = max(H, p[1]);
  const og = nums(attr(root[1], 'data-origin') || '0 0'), y0 = vb.length === 4 ? vb[1] : 0, x0 = vb.length === 4 ? vb[0] : 0, fx = p => [(p[0] - x0) * scale + (og[0] || 0), (H + y0 - p[1]) * scale + (og[1] || 0)];
  for (const q of polygons) { q.pts = q.pts.map(fx); if (polyArea(q.pts) < 0) q.pts.reverse(); }
  notes.push('y axis flipped so that +y points up; Bézier and elliptical-arc segments flattened to polylines');
  if (uses) lost.push(`${uses} <use> references not expanded`); if (textN) lost.push(`${textN} text elements ignored`); if (images) lost.push(`${images} embedded raster images ignored`); if (rounded) lost.push(`${rounded} rounded rectangles imported with sharp corners`);
  lost.push('fill rules, strokes widths, clipping and compound-path holes are not evaluated: each closed sub-path is a separate filled polygon');
  const meta = { viewBox: vb.length === 4 ? vb : null, width: wAttr, height: hAttr, layers: layerN || 1 };
  if (!polygons.length) return { kind: 'lattice', nodes: nodes.map(p => { const q = fx(p); return [q[0], q[1], 0]; }), edges, radius: 0, units, wireframe: true, meta, notes: [...notes, 'only open strokes found: imported as a node–edge network'], lost };
  const g = { kind: 'layout', polygons, units, meta, notes, lost }; if (edges.length) { g.wire = { nodes: nodes.map(p => { const q = fx(p); return [q[0], q[1], 0]; }), edges }; notes.push(`${openN} open strokes kept as wire segments`); }
  return g;
}

// ---------------------------------------------------------------- Gerber RS-274X, Excellon, CLI, G-code
export function parseGerber(input) {
  const t = toText(input); if (!/%FS|%MO|%ADD|G04|D0?[123]\*/.test(t)) fail('Gerber', 'no RS-274X commands found (%FS, %MO, %ADD, D01–D03)');
  const macros = (t.match(/%AM[^%]*%/g) || []).length, blocks = t.replace(/%AM[^%]*%/g, '').replace(/[%\r\n]/g, '').split('*').map(s => s.trim()).filter(Boolean);
  let ni = 2, nd = 4, trailing = false, units = 'mm', x = 0, y = 0, mode = 1, region = false, dark = true, ap = null, absC = true, contour = null, arcs = 0, unk = 0, sr = false;
  const aps = new Map(), polygons = [];
  const coord = s => { const neg = s[0] === '-'; let d = s.replace(/^[+-]/, ''); if (d.includes('.')) return +s; if (trailing) d = d.padEnd(ni + nd, '0'); else d = d.padStart(ni + nd, '0'); return (neg ? -1 : 1) * +d / 10 ** nd; };
  const add = pts => { if (pts.length >= 3) polygons.push({ layer: 1, pts: polyArea(pts) < 0 ? pts.slice().reverse() : pts, ...(dark ? {} : { clear: true }) }); };
  const flash = () => { if (!ap) return; const [a, bq, c] = ap.p; if (ap.s === 'C') add(arcPts(x, y, a / 2, 0, 2 * PI, true)); else if (ap.s === 'R') add([[x - a / 2, y - bq / 2], [x + a / 2, y - bq / 2], [x + a / 2, y + bq / 2], [x - a / 2, y + bq / 2]]);
    else if (ap.s === 'O') { const r = min(a, bq) / 2, hz = a >= bq, l = abs(a - bq) / 2; add(hz ? [...arcPts(x + l, y, r, -PI / 2, PI / 2, false), ...arcPts(x - l, y, r, PI / 2, 3 * PI / 2, false)] : [...arcPts(x, y + l, r, 0, PI, false), ...arcPts(x, y - l, r, PI, 2 * PI, false)]); }
    else if (ap.s === 'P') { const n = max(3, round(bq)), r0 = (c || 0) * PI / 180; add(Array.from({ length: n }, (_, k) => [x + a / 2 * cos(r0 + 2 * PI * k / n), y + a / 2 * sin(r0 + 2 * PI * k / n)])); } else unk++; };
  for (const bk of blocks) {
    let m;
    if ((m = /^FS([LT])([AI])X(\d)(\d)Y\d\d/.exec(bk))) { trailing = m[1] === 'T'; absC = m[2] === 'A'; ni = +m[3]; nd = +m[4]; continue; }
    if (/^MO(MM|IN)/.test(bk)) { units = bk[2] === 'M' ? 'mm' : 'in'; continue; }
    if ((m = /^ADD(\d+)([A-Za-z_.$][\w.$]*),?(.*)$/.exec(bk))) { aps.set(+m[1], { s: ['C', 'R', 'O', 'P'].includes(m[2]) ? m[2] : '?', p: m[3].split('X').map(Number), name: m[2] }); continue; }
    if (/^LP[DC]/.test(bk)) { dark = bk[2] === 'D'; continue; } if (/^SR/.test(bk)) { if (/X\d+Y\d+/.test(bk) && !/^SRX1Y1/.test(bk)) sr = true; continue; }
    if (/^(G04|LN|IN|IP|TF|TA|TO|TD|OF|AS|MI|SF|LM|LR|LS|AB)/.test(bk)) continue;
    let s = bk, g; while ((g = /^G(\d+)/.exec(s))) { const c = +g[1]; if (c === 1 || c === 2 || c === 3) mode = c; else if (c === 36) { region = true; contour = null; } else if (c === 37) { if (contour) add(contour); contour = null; region = false; } else if (c === 70) units = 'in'; else if (c === 71) units = 'mm'; else if (c === 90) absC = true; else if (c === 91) absC = false; s = s.slice(g[0].length); }
    if (!s || /^M0?[02]/.test(s)) continue;
    const X = /X([+-]?[\d.]+)/.exec(s), Y = /Y([+-]?[\d.]+)/.exec(s), I = /I([+-]?[\d.]+)/.exec(s), J = /J([+-]?[\d.]+)/.exec(s), Dm = /D(\d+)$/.exec(s), dc = Dm ? +Dm[1] : null;
    if (dc !== null && dc >= 10 && !X && !Y) { ap = aps.get(dc) || null; continue; }
    const nx = X ? (absC ? coord(X[1]) : x + coord(X[1])) : x, ny = Y ? (absC ? coord(Y[1]) : y + coord(Y[1])) : y, op = dc ?? 1;
    if (op === 1) {
      let path = [[x, y], [nx, ny]];
      if (mode !== 1 && (I || J)) { arcs++; const cx = x + (I ? coord(I[1]) : 0), cy = y + (J ? coord(J[1]) : 0), r = hypot(x - cx, y - cy), a0 = Math.atan2(y - cy, x - cx); let a1 = Math.atan2(ny - cy, nx - cx); if (mode === 2) { if (a1 >= a0 - 1e-12) a1 -= 2 * PI; } else if (a1 <= a0 + 1e-12) a1 += 2 * PI; path = arcPts(cx, cy, r, a0, a1, false); path[path.length - 1] = [nx, ny]; }
      if (region) { if (!contour) contour = [[x, y]]; contour.push(...path.slice(1)); }
      else if (ap && ap.p[0] > 0) { const pg = pathToPolygon(path, ap.p[0], { extend: ap.p[0] / 2 }); add(pg); } else unk++;
    } else if (op === 2) { if (region && contour) { add(contour); contour = null; } }
    x = nx; y = ny; if (op === 3) flash();
  }
  if (contour) add(contour);
  if (!polygons.length) fail('Gerber', 'no draws, flashes or regions could be converted');
  const clearN = polygons.filter(q => q.clear).length, lost = ['stroke ends are squared (extended by half the aperture) instead of rounded'];
  if (clearN) lost.push(`${clearN} clear-polarity (LPC) shapes are flagged "clear" but not subtracted from the dark image`); if (macros) lost.push(`${macros} aperture macros (AM) not interpreted: shapes using them are missing`); if (unk) lost.push(`${unk} operations used an undefined or unsupported aperture and were skipped`); if (sr) lost.push('step-and-repeat (SR) blocks not replicated');
  return { kind: 'layout', polygons, units, meta: { format: `${ni}.${nd}`, zeroOmission: trailing ? 'trailing' : 'leading', apertures: aps.size, arcs }, notes: [`coordinates in ${units} (format ${ni}.${nd})`, ...(arcs ? [`${arcs} circular interpolations discretised at 7.5° steps`] : [])], lost };
}
export function parseExcellon(input) {
  const L = lines(toText(input)).map(l => l.trim()).filter(Boolean); if (!L.some(l => /^M48/.test(l)) && !L.some(l => /^T\d+C/.test(l))) fail('Excellon', 'no M48 header or tool definitions found');
  let units = 'in', tz = false, tool = 0, x = 0, y = 0, implicit = 0, slots = 0; const tools = new Map(), polygons = [];
  const co = s => { if (s.includes('.')) return +s; implicit++; const neg = s[0] === '-', d = s.replace(/^[+-]/, ''), dec = units === 'mm' ? 3 : 4, tot = units === 'mm' ? 6 : 6; return (neg ? -1 : 1) * +(tz ? d.padStart(tot, '0') : d.padEnd(tot, '0')) / 10 ** dec; };
  for (const l of L) { let m;
    if (/^(METRIC|M71)/.test(l)) { units = 'mm'; tz = /TZ/.test(l) ? true : /LZ/.test(l) ? false : tz; } else if (/^(INCH|M72)/.test(l)) { units = 'in'; tz = /TZ/.test(l) ? true : /LZ/.test(l) ? false : tz; }
    else if ((m = /^T(\d+).*?C([\d.]+)/.exec(l))) tools.set(+m[1], +m[2]); else if ((m = /^T(\d+)$/.exec(l))) tool = +m[1];
    else if (/^[XY]/.test(l)) { if (/G85/.test(l)) slots++; const X = /X([+-]?[\d.]+)/.exec(l), Y = /Y([+-]?[\d.]+)/.exec(l); if (X) x = co(X[1]); if (Y) y = co(Y[1]); const d = tools.get(tool); if (d > 0) polygons.push({ layer: tool, pts: arcPts(x, y, d / 2, 0, 2 * PI, true) }); } }
  if (!polygons.length) fail('Excellon', 'no drill hits with a defined tool diameter found');
  return { kind: 'layout', polygons, units, meta: { tools: Object.fromEntries(tools), hits: polygons.length }, notes: [`drill hits as circles, layer = tool number, units ${units}`, ...(implicit ? [`${implicit} coordinates without a decimal point decoded with the assumed ${units === 'mm' ? '3.3' : '2.4'} format (${tz ? 'trailing' : 'leading'} zeros kept) — verify against the fabrication drawing`] : [])], lost: slots ? [`${slots} routed slots (G85) imported as single hits`] : [] };
}
export function parseCLI(input) {
  const b = toBytes(input), head = latin(b, 0, min(b.length, 4096)); if (!/\$\$HEADERSTART/.test(head)) fail('CLI', 'no $$HEADERSTART keyword');
  if (/\$\$BINARY/.test(head)) fail('CLI', 'binary CLI is not supported: export the slice file as ASCII CLI');
  const t = toText(b), u = +((/\$\$UNITS\s*\/\s*([\d.eE+-]+)/.exec(t) || [])[1]) || 1, polygons = [], layers = []; let cur = null, hatch = 0, openN = 0;
  for (const m of t.matchAll(/\$\$(LAYER|POLYLINE|HATCHES)\s*\/\s*([^$]*)/g)) { const a = nums(m[2]);
    if (m[1] === 'LAYER') { cur = { z: a[0] * u, polygons: 0, hatchLength: 0 }; layers.push(cur); }
    else if (!cur) continue;
    else if (m[1] === 'POLYLINE') { const n = a[2], pts = []; for (let k = 0; k < n; k++) pts.push([a[3 + 2 * k] * u, a[4 + 2 * k] * u]); if (pts.some(p => !Number.isFinite(p[0] + p[1]))) fail('CLI', `POLYLINE in layer ${layers.length} is shorter than its declared ${n} points`); if (a[1] === 2) { openN++; continue; } if (pts.length > 1 && hypot(pts[0][0] - pts[n - 1][0], pts[0][1] - pts[n - 1][1]) < 1e-9) pts.pop(); if (pts.length >= 3) { polygons.push({ layer: layers.length, pts, z: cur.z, ...(a[1] === 0 ? { hole: true } : {}) }); cur.polygons++; } }
    else { const n = a[1]; for (let k = 0; k < n; k++) { const q = a.slice(2 + 4 * k, 6 + 4 * k); cur.hatchLength += hypot(q[2] - q[0], q[3] - q[1]) * u; } hatch += n; } }
  if (!layers.length) fail('CLI', 'no $$LAYER records found');
  const dz = layers.length > 1 ? (layers[layers.length - 1].z - layers[0].z) / (layers.length - 1) : 0;
  return { kind: 'layout', polygons, units: 'mm', layers, meta: { unitsFactor_mm: u, nLayers: layers.length, layerThickness_mm: dz, hatches: hatch }, notes: [`${layers.length} slices; layout layer number = slice index, slice height in polygons[].z (mm)`, 'clockwise (internal) contours are flagged hole:true and are not subtracted'], lost: [...(hatch ? [`${hatch} hatch vectors summarised by length per layer, not stored`] : []), ...(openN ? [`${openN} open polylines ignored`] : [])] };
}
export function parseGcode(input, { maxSegments = 60000 } = {}) {
  const L = lines(toText(input)); let x = 0, y = 0, z = 0, e = 0, f = 0, absXY = true, absE = true, inch = 1, hasE = false, mode = 0, nMove = 0;
  const nodes = [], edges = [], zs = new Map(); let dep = 0, travel = 0, time = 0, eTot = 0, prevNode = -1, kept = 0, total = 0, lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const raw of L) if (/\bE[-+]?[\d.]/.test(raw.split(';')[0]) && /^\s*(N\d+\s+)?G0?1\b/i.test(raw)) { hasE = true; break; }
  const seg = (nx, ny, nz, work, len) => { if (work) { dep += len; total++; const zk = N9(nz); zs.set(zk, (zs.get(zk) || 0) + len); for (const [k, v] of [[0, nx], [1, ny], [2, nz]]) { if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; } if (kept < maxSegments) { if (prevNode < 0) { nodes.push([x, y, z]); prevNode = nodes.length - 1; } nodes.push([nx, ny, nz]); edges.push([prevNode, nodes.length - 1]); prevNode = nodes.length - 1; kept++; } } else { travel += len; prevNode = -1; } if (f > 0) time += len / f * 60; };
  for (const raw of L) {
    const l = raw.replace(/\(.*?\)/g, '').split(';')[0].trim().toUpperCase(); if (!l) continue; const w = {}; for (const m of l.matchAll(/([A-Z])\s*([-+]?(?:\d+\.?\d*|\.\d+))/g)) if (!(m[1] in w) || m[1] === 'G') { if (m[1] === 'G') (w.G = w.G || []).push(+m[2]); else w[m[1]] = +m[2]; }
    for (const g of w.G || []) { if (g === 90) { absXY = true; absE = true; } else if (g === 91) { absXY = false; absE = false; } else if (g === 20) inch = 25.4; else if (g === 21) inch = 1; else if (g === 92) { if ('E' in w) e = w.E; if ('X' in w) x = w.X * inch; if ('Y' in w) y = w.Y * inch; if ('Z' in w) z = w.Z * inch; } else if (g <= 3) mode = g; }
    if (w.M === 82) absE = true; else if (w.M === 83) absE = false; if ((w.G || []).includes(92) || (w.G || []).includes(28)) continue;
    if ('F' in w) f = w.F * inch; if (!('X' in w || 'Y' in w || 'Z' in w || 'E' in w)) continue; if (!(w.G || []).length && nMove === 0 && !('X' in w)) continue;
    const nx = 'X' in w ? (absXY ? w.X * inch : x + w.X * inch) : x, ny = 'Y' in w ? (absXY ? w.Y * inch : y + w.Y * inch) : y, nz = 'Z' in w ? (absXY ? w.Z * inch : z + w.Z * inch) : z, ne = 'E' in w ? (absE ? w.E : e + w.E) : e, de = ne - e; nMove++;
    const work = mode !== 0 && (hasE ? de > 1e-9 && (nx !== x || ny !== y || nz !== z) : true); if (de > 0) eTot += de;
    if ((mode === 2 || mode === 3) && ('I' in w || 'J' in w)) { const cx = x + (w.I || 0) * inch, cy = y + (w.J || 0) * inch, r = hypot(x - cx, y - cy), a0 = Math.atan2(y - cy, x - cx); let a1 = Math.atan2(ny - cy, nx - cx); if (mode === 2) { if (a1 >= a0 - 1e-12) a1 -= 2 * PI; } else if (a1 <= a0 + 1e-12) a1 += 2 * PI; const pts = arcPts(cx, cy, r, a0, a1, false); let px = x, py = y; pts.slice(1).forEach((q, k) => { const zz = z + (nz - z) * (k + 1) / (pts.length - 1), ln = hypot(q[0] - px, q[1] - py); const sx = x, sy = y; x = px; y = py; seg(q[0], q[1], zz, work, ln); x = sx; y = sy; px = q[0]; py = q[1]; }); }
    else { const len = hypot(nx - x, ny - y, nz - z); if (len > 0) seg(nx, ny, nz, work, len); }
    x = nx; y = ny; z = nz; e = ne;
  }
  if (!total) fail('G-code', hasE ? 'no extruding moves found' : 'no G1/G2/G3 feed moves found');
  const zl = [...zs.keys()].map(Number).sort((a, c) => a - c), dzs = zl.slice(1).map((v, i) => v - zl[i]).sort((a, c) => a - c), layerHeight = dzs.length ? dzs[floor(dzs.length / 2)] : 0;
  const stats = { mode: hasE ? 'additive (extrusion)' : 'subtractive / motion (feed moves)', workMoves: total, nLayers: zl.length, layerHeight_mm: layerHeight, workLength_mm: dep, travelLength_mm: travel, extrusion_E: eTot, estimatedTime_s: time, bbox: [lo, hi] };
  const g = { kind: 'lattice', nodes, edges, radius: 0, units: 'mm', toolpath: true, stats, layers: zl.slice(0, 2000).map(v => ({ z: v, length_mm: zs.get(N9(v)) })), meta: stats, notes: [`toolpath of ${total} ${hasE ? 'extruding' : 'feed'} moves imported as a node–edge polyline network (mm)`, 'time estimate = Σ length / programmed feed rate (accelerations ignored: a lower bound)'], lost: ['temperatures, fan, tool changes and machine-specific M-codes are not interpreted', 'a toolpath is a set of manufacturing instructions, not a solid: bead width and overlap are not reconstructed'] };
  if (total > kept) g.truncated = `only the first ${kept} of ${total} work moves are kept as geometry (statistics cover the whole file)`;
  return g;
}

// ---------------------------------------------------------------- TIFF
const CM_PER = { nm: 1e-7, 'µm': 1e-4, um: 1e-4, micron: 1e-4, mm: 0.1, cm: 1, m: 100, 'Å': 1e-8, in: 2.54, inch: 2.54 };
// Multi-page, uncompressed, 8-bit greyscale, little-endian classic TIFF: one page per z slice, phase ids stored as
// grey values (contrast:true stretches them to 0–255 for viewing, which is lossy for the ids). ImageJ-style
// ImageDescription carries the z spacing and unit; X/YResolution hold pixels per that unit (ImageJ/Fiji convention,
// ResolutionUnit = none), so non-ImageJ readers see an uncalibrated image.
export function toTIFFStack(vox, { contrast = false } = {}) {
  if (!vox || vox.kind !== 'voxels') fail('TIFF', `needs a voxels geometry (got ${vox && vox.kind})`);
  const { nx, ny, nz, data } = vox, N = nx * ny, sp = vox.spacing || [vox.voxel || 1, vox.voxel || 1, vox.voxel || 1], cm = CM_PER[vox.units], un = vox.units === 'um' ? 'micron' : vox.units === 'µm' ? 'micron' : (vox.units || 'pixel');
  const desc = ENC.encode(`ImageJ=1.53t\nimages=${nz}\nslices=${nz}\nunit=${un}\nspacing=${N9(sp[2])}\nloop=false\n\0`);
  const rat = v => { if (!cm || !(v > 0)) return [1, 1]; const ppu = 1 / v; let den = 1; while (ppu * den < 1e6 && den < 1e6) den *= 10; return [min(4294967295, round(ppu * den)), den]; }, rx = rat(sp[0]), ry = rat(sp[1]); // ImageJ convention: pixels per `unit`, ResolutionUnit = none
  let mx = 1; if (contrast) for (let i = 0; i < data.length; i++) if (data[i] > mx) mx = data[i];
  const w = new Writer(8 + nz * (N + 200) + desc.length); w.str('II').u16(42).u32(8);
  for (let k = 0; k < nz; k++) {
    const first = k === 0, ne = first ? 13 : 12, ifd = w.n, ext = ifd + 2 + 12 * ne + 4, dOff = ext + 16, pix = dOff + (first ? desc.length + (desc.length % 2) : 0), next = k === nz - 1 ? 0 : pix + N + (N % 2);
    const e = (tag, type, cnt, val) => { w.u16(tag).u16(type).u32(cnt); if (type === 3 && cnt === 1) w.u16(val).u16(0); else w.u32(val); };
    w.u16(ne); e(256, 4, 1, nx); e(257, 4, 1, ny); e(258, 3, 1, 8); e(259, 3, 1, 1); e(262, 3, 1, 1); if (first) e(270, 2, desc.length, dOff); e(273, 4, 1, pix); e(277, 3, 1, 1); e(278, 4, 1, ny); e(279, 4, 1, N); e(282, 5, 1, ext); e(283, 5, 1, ext + 8); e(296, 3, 1, 1);
    w.u32(next).u32(rx[0]).u32(rx[1]).u32(ry[0]).u32(ry[1]); if (first) { w.bytes(desc); if (desc.length % 2) w.u8(0); }
    w.need(N + 1); for (let r = 0; r < ny; r++) { const src = nx * ((ny - 1 - r) + ny * k); if (contrast) for (let i = 0; i < nx; i++) w.b[w.n + r * nx + i] = round(data[src + i] * 255 / mx); else w.b.set(data.subarray(src, src + nx), w.n + r * nx); } w.n += N; if (N % 2) w.u8(0);
  }
  return w.out();
}
function lzwTiff(src, expected) {
  const out = new Uint8Array(expected), prefix = new Uint16Array(4098), suffix = new Uint8Array(4098), stack = new Uint8Array(4100); let op = 0, bitBuf = 0, bitCnt = 0, sp = 0, width = 9, next = 258, old = -1, first = 0;
  for (let i = 0; i < 256; i++) suffix[i] = i;
  while (op < expected) {
    while (bitCnt < width) { if (sp >= src.length) return out; bitBuf = (bitBuf << 8 | src[sp++]) >>> 0; bitCnt += 8; }
    const code = bitBuf >>> (bitCnt - width) & ((1 << width) - 1); bitCnt -= width; bitBuf &= (1 << bitCnt) - 1;
    if (code === 257) break; if (code === 256) { width = 9; next = 258; old = -1; continue; }
    if (old < 0) { out[op++] = code; old = first = code; continue; }
    let cur = code, n = 0; if (code >= next) { stack[n++] = first; cur = old; }
    while (cur >= 258) { stack[n++] = suffix[cur]; cur = prefix[cur]; }
    first = cur; stack[n++] = cur; while (n && op < expected) out[op++] = stack[--n];
    if (next < 4096) { prefix[next] = old; suffix[next] = first; next++; if (next >= (1 << width) - 1 && width < 12) width++; }
    old = code;
  }
  return out;
}
function packBits(src, expected) { const out = new Uint8Array(expected); let op = 0, p = 0; while (p < src.length && op < expected) { const n = src[p++] << 24 >> 24; if (n >= 0) { for (let i = 0; i <= n && op < expected; i++) out[op++] = src[p++]; } else if (n !== -128) { const v = src[p++]; for (let i = 0; i < 1 - n && op < expected; i++) out[op++] = v; } } return out; }
// Classic and BigTIFF, strips or tiles; uncompressed, PackBits, LZW (with horizontal predictor) and Deflate;
// 1/8/16/32-bit integer or 32-bit float samples; grey or RGB(A) (converted to luminance). Pages → z slices.
export function parseTIFF(input, { threshold } = {}) {
  const b = toBytes(input), le = b[0] === 0x49 && b[1] === 0x49; if (!le && !(b[0] === 0x4D && b[1] === 0x4D)) fail('TIFF', 'missing II/MM byte-order mark');
  const v = view(b), u16 = o => v.getUint16(o, le), u32 = o => v.getUint32(o, le), u64 = o => Number(v.getBigUint64(o, le)), magic = u16(2), big = magic === 43; if (magic !== 42 && !big) fail('TIFF', `unexpected magic number ${magic}`);
  const TS = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8, 13: 4, 16: 8, 17: 8, 18: 8 }, pages = []; let ifd = big ? u64(8) : u32(4);
  while (ifd && ifd + 2 < b.length && pages.length < 8192) { const n = big ? u64(ifd) : u16(ifd), es = big ? 20 : 12, base = ifd + (big ? 8 : 2), tags = {}; if (base + n * es > b.length) break; for (let i = 0; i < n; i++) { const o = base + i * es, type = u16(o + 2), cnt = big ? u64(o + 4) : u32(o + 4), vo = o + (big ? 12 : 8), inl = (TS[type] || 1) * cnt <= (big ? 8 : 4); tags[u16(o)] = { type, cnt, at: inl ? vo : (big ? u64(vo) : u32(vo)) }; } pages.push(tags); ifd = big ? u64(base + n * es) : u32(base + n * es); }
  if (!pages.length) fail('TIFF', 'no image file directory found');
  const val = (t, i = 0) => { if (!t) return undefined; const o = t.at + i * (TS[t.type] || 1); switch (t.type) { case 3: return u16(o); case 4: case 13: return u32(o); case 16: case 18: return u64(o); case 5: return u32(o) / (u32(o + 4) || 1); case 8: return v.getInt16(o, le); case 9: return v.getInt32(o, le); case 11: return v.getFloat32(o, le); case 12: return v.getFloat64(o, le); default: return b[o]; } };
  const list = t => { if (!t) return []; if (t.cnt > 4e6 || t.at + t.cnt * (TS[t.type] || 1) > b.length) fail('TIFF', 'a strip / tile table points outside the file (truncated or corrupt)'); return Array.from({ length: t.cnt }, (_, i) => val(t, i)); };
  const T0 = pages[0], W = val(T0[256]), H = val(T0[257]); if (!(W > 0 && H > 0)) fail('TIFF', 'ImageWidth/ImageLength missing');
  const same = pages.filter(t => val(t[256]) === W && val(t[257]) === H && !(val(t[254]) & 1)), nz = same.length; if (W * H * nz > 6e7) fail('TIFF', `stack of ${W} × ${H} × ${nz} voxels is too large for in-browser import: crop or bin the volume first`);
  if (W * H * nz * (val(T0[258]) || 1) * (val(T0[277]) || 1) / 8 > b.length * ((val(T0[259]) || 1) === 1 ? 1 : 1040)) fail('TIFF', `the directory declares ${W} × ${H} × ${nz} pixels, more than a ${b.length}-byte file can hold (corrupt header)`);
  const vol = new Float32Array(W * H * nz), notes = [], lost = []; let info = null;
  same.forEach((t, k) => {
    const bps = val(t[258]) || 1, comp = val(t[259]) || 1, photo = val(t[262]) ?? 1, spp = val(t[277]) || 1, sf = val(t[339]) || 1, pred = val(t[317]) || 1, planar = val(t[284]) || 1;
    if (![1, 8, 16, 32].includes(bps)) fail('TIFF', `${bps} bits per sample is not supported`); if (planar === 2 && spp > 1) fail('TIFF', 'planar (separate-plane) RGB storage is not supported');
    const tiled = !!t[322], tw = tiled ? val(t[322]) : W, th = tiled ? val(t[323]) : (val(t[278]) || H), offs = list(t[tiled ? 324 : 273]), cnts = list(t[tiled ? 325 : 279]), rowB = Math.ceil(tw * spp * bps / 8), B = bps / 8, img = new Uint8Array(Math.ceil(W * spp * bps / 8) * H), across = tiled ? Math.ceil(W / tw) : 1, iRow = Math.ceil(W * spp * bps / 8);
    offs.forEach((off, ci) => {
      const rows = tiled ? th : min(th, H - ci * th), exp = rowB * rows, raw = b.subarray(off, off + (cnts[ci] ?? exp)); let d;
      if (comp === 1) d = raw; else if (comp === 5) d = lzwTiff(raw, exp); else if (comp === 32773) d = packBits(raw, exp); else if (comp === 8 || comp === 32946) d = unzlib(raw, exp); else fail('TIFF', `compression scheme ${comp} is not supported (supported: none, LZW, PackBits, Deflate)`);
      if (pred === 2 && bps >= 8) { if (d === raw) d = raw.slice(); const dv = view(d); for (let r = 0; r < rows; r++) { const o = r * rowB; if (B === 1) for (let i = spp; i < rowB && o + i < d.length; i++) d[o + i] = d[o + i] + d[o + i - spp]; else if (B === 2) for (let i = spp; i < tw * spp && o + 2 * i + 1 < d.length; i++) dv.setUint16(o + 2 * i, dv.getUint16(o + 2 * i, le) + dv.getUint16(o + 2 * (i - spp), le), le); } }
      const x0 = tiled ? (ci % across) * tw : 0, y0 = tiled ? floor(ci / across) * th : ci * th, cb = x0 * spp * bps / 8;
      for (let r = 0; r < rows && y0 + r < H; r++) { const n = min(rowB, iRow - cb, d.length - r * rowB); if (n > 0) img.set(d.subarray(r * rowB, r * rowB + n), (y0 + r) * iRow + cb); }
    });
    const iv = view(img), samp = (o, s) => bps === 8 ? img[o + s] : bps === 16 ? (sf === 2 ? iv.getInt16(o + 2 * s, le) : iv.getUint16(o + 2 * s, le)) : sf === 3 ? iv.getFloat32(o + 4 * s, le) : sf === 2 ? iv.getInt32(o + 4 * s, le) : iv.getUint32(o + 4 * s, le), top = bps === 1 ? 1 : bps === 8 ? 255 : bps === 16 ? 65535 : 1;
    for (let y = 0; y < H; y++) { const dst = W * ((H - 1 - y) + H * k); for (let x = 0; x < W; x++) { let g; if (bps === 1) g = img[y * iRow + (x >> 3)] >> (7 - (x & 7)) & 1; else { const o = y * iRow + x * spp * B; g = spp >= 3 && photo === 2 ? 0.299 * samp(o, 0) + 0.587 * samp(o, 1) + 0.114 * samp(o, 2) : samp(o, 0); } vol[dst + x] = photo === 0 ? top - g : g; } }
    if (!info) info = { bitsPerSample: bps, samplesPerPixel: spp, compression: ({ 1: 'none', 5: 'LZW', 32773: 'PackBits', 8: 'Deflate', 32946: 'Deflate' })[comp], photometric: photo, sampleFormat: sf === 3 ? 'float' : sf === 2 ? 'int' : 'uint', tiled };
    if (k === 0) { if (photo === 3) lost.push('palette-colour image: colour-map indices are used as values'); if (spp >= 3 && photo === 2) notes.push('RGB samples converted to luminance (0.299 R + 0.587 G + 0.114 B)'); if (photo === 0) notes.push('WhiteIsZero photometric interpretation inverted'); }
  });
  let px = 1, pz = 1, units = 'px'; const descr = T0[270] ? latin(b, T0[270].at, min(T0[270].cnt, 200000)) : '', ru = val(T0[296]) ?? 2, xr = val(T0[282]);
  const ome = /PhysicalSizeX="([\d.eE+-]+)"/.exec(descr), ij = /unit=(\S+)/.exec(descr), ijs = /spacing=([\d.eE+-]+)/.exec(descr);
  if (ome) { px = +ome[1]; units = ((/PhysicalSizeXUnit="([^"]+)"/.exec(descr) || [])[1] || 'µm'); pz = +((/PhysicalSizeZ="([\d.eE+-]+)"/.exec(descr) || [])[1] || px); notes.push('voxel size from the OME-XML header (PhysicalSizeX/Z)'); lost.push('OME channel, time and plane metadata not interpreted: all pages are stacked along z'); }
  else if (ij && xr > 0 && !/^pixels?$/i.test(ij[1])) { units = /^(micron|um|µm|\\u00B5m)$/i.test(ij[1]) ? 'µm' : ij[1]; px = 1 / xr; pz = ijs ? +ijs[1] : px; notes.push('voxel size from the ImageJ description (unit, spacing) and X resolution'); }
  else if (xr > 0 && (ru === 2 || ru === 3) && xr !== 72 && xr !== 96) { px = (ru === 2 ? 25.4e3 : 1e4) / xr; pz = px; units = 'µm'; notes.push(`pixel size ${N9(px)} µm from the resolution tags; slice spacing assumed equal (not stored)`); }
  else notes.push('no physical calibration in the file: voxel size 1 px (set the scale manually)');
  if (same.length < pages.length) lost.push(`${pages.length - same.length} directories of a different size (thumbnails / pyramid levels) skipped`);
  const g = scalarToVoxels(vol, W, H, nz, { voxel: px, units, spacing: [px, px, pz], threshold });
  return Object.assign(g, { meta: { ...info, pages: pages.length, width: W, height: H, bigTiff: big, byteOrder: le ? 'little-endian' : 'big-endian', description: descr.slice(0, 300).replace(/\0/g, '') }, notes: [...notes, ...g.notes, 'image row 0 (top) is mapped to the highest y'], lost });
}

// ---------------------------------------------------------------- HDF5
// Minimal file: version-0 superblock, root group (symbol table: v1 B-tree + local heap + one symbol node) and one
// contiguous dataset with version-1 object header. Shape is (nz, ny, nx) in C order, i.e. flat index x + nx·(y + ny·z).
export function toHDF5(g, { name = 'data' } = {}) {
  let dims, raw, f64 = false; const attrs = [];
  if (g && g.kind === 'voxels') { dims = [g.nz, g.ny, g.nx]; raw = g.data; const sp = g.spacing || [g.voxel || 1, g.voxel || 1, g.voxel || 1]; attrs.push(['voxel_size', g.voxel || 1], ['spacing_xyz', sp], ['origin_xyz', g.origin || [0, 0, 0]], ['units', String(g.units || '')], ['axis_order', 'z,y,x (C order): flat index = x + nx*(y + ny*z)']); }
  else if (g && g.values && g.dims) { dims = [...g.dims]; f64 = !(g.values instanceof Uint8Array); if (f64) { const w = new Writer(g.values.length * 8); for (let i = 0; i < g.values.length; i++) w.f64(g.values[i]); raw = w.out(); } else raw = g.values; for (const [k, v] of Object.entries(g.attrs || {})) attrs.push([k, v]); }
  else fail('HDF5', 'needs a voxels geometry or a field {dims:[…], values:TypedArray}');
  const pad8 = n => (n + 7) & ~7, nm = ENC.encode(String(name).replace(/[\/\0]/g, '_') || 'data');
  const dtU8 = () => new Writer(16).u8(0x10).u8(0).u8(0).u8(0).u32(1).u16(0).u16(8).out(), dtF64 = () => new Writer(24).u8(0x11).u8(0x20).u8(0x3F).u8(0).u32(8).u16(0).u16(64).u8(52).u8(11).u8(0).u8(52).u32(1023).out(), dtStr = n => new Writer(8).u8(0x13).u8(0x10).u8(0).u8(0).u32(n).out(); // null-terminated, UTF-8
  const space = d => { const w = new Writer(8 + 8 * d.length).u8(1).u8(d.length).u8(0).u8(0).u32(0); d.forEach(v => w.u64(v)); return w.out(); };
  const msg = (type, data, flags = 0) => { const w = new Writer(8 + pad8(data.length)).u16(type).u16(pad8(data.length)).u8(flags).zero(3).bytes(data); return w.align(8).out(); };
  const attrMsg = ([k, val]) => { const kn = ENC.encode(k + '\0'); let dt, ds, data; if (typeof val === 'string') { const s = ENC.encode(val + '\0'); dt = dtStr(s.length); ds = space([]); data = s; } else { const a = Array.isArray(val) ? val : [val]; dt = dtF64(); ds = space(Array.isArray(val) ? [a.length] : []); const w = new Writer(8 * a.length); a.forEach(x => w.f64(x)); data = w.out(); }
    const w = new Writer(64).u8(1).u8(0).u16(kn.length).u16(dt.length).u16(ds.length).bytes(kn).align(8).bytes(dt).align(8).bytes(ds).align(8).bytes(data); return msg(0x0C, w.out()); };
  const ohdr = msgs => { let size = 0; msgs.forEach(m => size += m.length); const w = new Writer(16 + size).u8(1).u8(0).u16(msgs.length).u32(1).u32(size).zero(4); msgs.forEach(m => w.bytes(m)); return w.out(); };
  const U = 0xFFFFFFFF, rootOH = 96, bt = rootOH + 40, heap = bt + 544, heapData = heap + 32, segSize = 8 + pad8(nm.length + 1), snod = heapData + segSize, dsOH = snod + 328;
  const dsMsgsNoLayout = [msg(0x01, space(dims)), msg(0x03, f64 ? dtF64() : dtU8(), 1), msg(0x05, new Writer(8).u8(2).u8(2).u8(2).u8(1).u32(0).out(), 1)], attrM = attrs.map(attrMsg);
  let dsSize = 16 + 8 + 24; [...dsMsgsNoLayout, ...attrM].forEach(m => dsSize += m.length); const dataAt = pad8(dsOH + dsSize), total = dataAt + raw.length;
  const layout = msg(0x08, new Writer(18).u8(3).u8(1).u64(dataAt).u64(raw.length).out()), w = new Writer(total);
  w.bytes([0x89, 0x48, 0x44, 0x46, 0x0D, 0x0A, 0x1A, 0x0A]).u8(0).u8(0).u8(0).u8(0).u8(0).u8(8).u8(8).u8(0).u16(4).u16(16).u32(0).u64(0).u32(U).u32(U).u64(total).u32(U).u32(U).u64(0).u64(rootOH).u32(1).u32(0).u64(bt).u64(heap);
  w.bytes(ohdr([msg(0x11, new Writer(16).u64(bt).u64(heap).out())]));
  w.str('TREE').u8(0).u8(0).u16(1).u32(U).u32(U).u32(U).u32(U).u64(0).u64(snod).u64(8); w.n = heap;
  w.str('HEAP').u8(0).zero(3).u64(segSize).u64(1).u64(heapData); w.zero(8).bytes(nm); w.n = snod;
  w.str('SNOD').u8(1).u8(0).u16(1).u64(8).u64(dsOH).u32(0).u32(0).zero(16); w.n = dsOH;
  w.bytes(ohdr([...dsMsgsNoLayout, layout, ...attrM])); w.n = dataAt; w.bytes(raw);
  return w.out();
}
// Reader for the classic HDF5 file layout: superblock v0–v3, old-style groups (symbol tables) and compact link
// messages, version 1/2 object headers, contiguous / compact / chunked (v1 B-tree; deflate, shuffle, fletcher32) datasets
// of fixed-point and IEEE float types. Everything else is listed in the inventory with the reason it was not read.
function h5Open(b) {
  let sig = -1; for (let o = 0; o < min(b.length, 1 << 20); o = o ? o * 2 : 512) if (b[o] === 0x89 && latin(b, o + 1, 3) === 'HDF') { sig = o; break; }
  if (sig < 0) fail('HDF5', 'signature \\x89HDF not found'); const v = view(b), sbv = b[sig + 8]; let so, sl, rootOH = null, rootBT = null, rootHeap = null, base = 0;
  const rd = (o, n) => n === 8 ? Number(v.getBigUint64(o, true)) : n === 4 ? v.getUint32(o, true) : n === 2 ? v.getUint16(o, true) : b[o];
  if (sbv <= 1) { so = b[sig + 13]; sl = b[sig + 14]; const p0 = sig + (sbv === 1 ? 28 : 24); base = rd(p0, so); const re = p0 + 4 * so; rootOH = rd(re + so, so); if (v.getUint32(re + 2 * so, true) === 1) { rootBT = rd(re + 2 * so + 8, so); rootHeap = rd(re + 2 * so + 8 + so, so); } }
  else if (sbv <= 3) { so = b[sig + 9]; sl = b[sig + 10]; base = rd(sig + 12, so); rootOH = rd(sig + 12 + 3 * so, so); } else fail('HDF5', `superblock version ${sbv} is not supported`);
  const A = a => a + base, undef = a => !Number.isFinite(a) || a >= b.length || a === 0xFFFFFFFF || a > 2 ** 52;
  const messages = addr => { // → [{type, o, size}]
    const out = [], q = []; let a = A(addr);
    if (latin(b, a, 4) === 'OHDR') { const fl = b[a + 5]; let p = a + 6; if (fl & 0x20) p += 16; if (fl & 0x10) p += 4; const cs = rd(p, 1 << (fl & 3)); p += 1 << (fl & 3); q.push([p, p + cs]); const co = fl & 4 ? 2 : 0;
      while (q.length) { let [s, e] = q.shift(); while (s + 4 + co <= e - 0) { const type = b[s], size = v.getUint16(s + 1, true), o = s + 4 + co; if (o + size > e) break; if (type === 0x10) { const ca = A(rd(o, so)), cl = rd(o + so, sl); if (latin(b, ca, 4) === 'OCHK') q.push([ca + 4, ca + cl - 4]); } else if (type) out.push({ type, o, size, v2: true }); s = o + size; } } return out; }
    if (b[a] !== 1) return out; const nm = v.getUint16(a + 2, true), hs = v.getUint32(a + 8, true); q.push([a + 16, a + 16 + hs]);
    while (q.length && out.length < nm + 64) { let [s, e] = q.shift(); while (s + 8 <= e) { const type = v.getUint16(s, true), size = v.getUint16(s + 2, true), o = s + 8; if (type === 0x10) q.push([A(rd(o, so)), A(rd(o, so)) + rd(o + so, sl)]); else if (type) out.push({ type, o, size }); s = o + size; } }
    return out;
  };
  const dtype = m => { const cls = b[m.o] & 15, bits = b[m.o + 1], size = v.getUint32(m.o + 4, true), be = !!(bits & 1); if (cls === 0 && [1, 2, 4, 8].includes(size)) return { name: (bits & 8 ? 'int' : 'uint') + size * 8, size, be }; if (cls === 1 && (size === 4 || size === 8)) return { name: 'float' + size * 8, size, be }; return { name: ['fixed', 'float', 'time', 'string', 'bitfield', 'opaque', 'compound', 'reference', 'enum', 'vlen', 'array'][cls] || `class${cls}`, size, be, unsupported: cls !== 0 && cls !== 1, string: cls === 3 }; };
  const dspace = m => { const ver = b[m.o], rank = b[m.o + 1], off = ver === 1 ? 8 : 4, d = []; for (let i = 0; i < rank; i++) d.push(rd(m.o + off + i * sl, sl)); return d; };
  const attrOf = m => { try { const ver = b[m.o]; let p = m.o + 2; const ns = v.getUint16(p, true), ts = v.getUint16(p + 2, true), ss = v.getUint16(p + 4, true); p += 6; if (ver === 3) p++; const al = n => ver === 1 ? (n + 7) & ~7 : n, name = latin(b, p, ns).replace(/\0.*$/, ''); p += al(ns); const dt = dtype({ o: p }); const dsp = dspace({ o: p + al(ts) }); p += al(ts) + al(ss); const n = dsp.reduce((x, y) => x * y, 1); if (dt.string) return [name, toText(b.subarray(p, p + dt.size)).replace(/\0[\s\S]*$/, '')]; if (dt.unsupported || n > 64) return [name, `<${dt.name}>`]; const a = Array.from(readTyped(b, p, n, dt.name, !dt.be)); return [name, dsp.length ? a : a[0]]; } catch (e) { return null; } };
  const datasets = [], groups = [], notes = [], seen = new Set();
  const visit = (addr, path, depth) => {
    if (undef(A(addr)) || seen.has(addr) || depth > 12 || datasets.length > 2000) return; seen.add(addr); const ms = messages(addr), find = t => ms.find(m => m.type === t);
    const st = find(0x11), links = ms.filter(m => m.type === 0x06), lay = find(0x08), dsm = find(0x01), dtm = find(0x03);
    if (lay && dsm && dtm) {
      const shape = dspace(dsm), dt = dtype(dtm), count = shape.reduce((x, y) => x * y, 1), d = { path: path || '/', shape, dtype: dt.name, layout: '', filters: [], readable: false, why: '', attrs: {} };
      for (const m of ms) if (m.type === 0x0C) { const kv = attrOf(m); if (kv) d.attrs[kv[0]] = kv[1]; }
      const fm = find(0x0B); if (fm) { const ver = b[fm.o], nf = b[fm.o + 1]; let p = fm.o + (ver === 1 ? 8 : 2); for (let i = 0; i < nf; i++) { const id = v.getUint16(p, true); let nl = 0; if (ver === 1 || id >= 256) { nl = v.getUint16(p + 2, true); p += 2; } const ncd = v.getUint16(p + 4, true); p += 6 + (ver === 1 ? (nl + 7) & ~7 : nl) + 4 * ncd + (ver === 1 && ncd % 2 ? 4 : 0); d.filters.push(id); } }
      const lv = b[lay.o], cls = lv >= 3 ? b[lay.o + 1] : b[lay.o + 2]; d.layout = ['compact', 'contiguous', 'chunked', 'virtual'][cls] || `class ${cls}`;
      const bytes = count * dt.size, FN = { 1: 'deflate', 2: 'shuffle', 3: 'fletcher32', 4: 'szip', 32000: 'lzf', 32001: 'blosc', 32004: 'lz4', 32015: 'zstd' }; d.filterNames = d.filters.map(f => FN[f] || `filter ${f}`);
      if (dt.unsupported) d.why = `datatype ${dt.name} is not numeric`; else if (bytes > 2.5e8) d.why = 'dataset larger than 250 MB (too large for in-browser import)'; else if (!count) d.why = 'empty dataset';
      else if (cls === 1) { const ad = lv >= 3 ? rd(lay.o + 2, so) : rd(lay.o + 8, so); if (undef(A(ad))) d.why = 'no storage allocated'; else { d.readable = true; d.read = () => readTyped(b, A(ad), count, dt.name, !dt.be); } }
      else if (cls === 0 && lv >= 3) { d.readable = true; d.read = () => readTyped(b, lay.o + 4, count, dt.name, !dt.be); }
      else if (cls === 2 && lv === 3) {
        const nd = b[lay.o + 2], bta = rd(lay.o + 3, so), cd = []; for (let i = 0; i < nd; i++) cd.push(v.getUint32(lay.o + 3 + so + 4 * i, true));
        const bad = d.filters.filter(f => ![1, 2, 3].includes(f)); if (bad.length) d.why = `chunk filter(s) ${d.filterNames.join(', ')} not supported (only deflate, shuffle, fletcher32)`; else if (undef(A(bta))) d.why = 'no chunks allocated';
        else { d.readable = true; d.read = () => { const rank = shape.length, es = dt.size, out = new Uint8Array(bytes), csz = cd.slice(0, rank), cn = csz.reduce((x, y) => x * y, 1) * es;
          const node = a => { a = A(a); if (latin(b, a, 4) !== 'TREE' || b[a + 4] !== 1) fail('HDF5', 'chunk B-tree node is corrupt'); const lvl = b[a + 5], n = v.getUint16(a + 6, true), ks = 8 + 8 * nd; let p = a + 8 + 2 * so;
            for (let i = 0; i < n; i++, p += ks + so) { const child = rd(p + ks, so); if (lvl > 0) { node(child); continue; } const size = v.getUint32(p, true), mask = v.getUint32(p + 4, true), off = []; for (let k = 0; k < rank; k++) off.push(rd(p + 8 + 8 * k, 8)); let c = b.subarray(A(child), A(child) + size);
              for (let fi = d.filters.length - 1; fi >= 0; fi--) { if (mask >> fi & 1) continue; const f = d.filters[fi]; if (f === 1) c = unzlib(c, cn); else if (f === 3) c = c.subarray(0, c.length - 4); else if (f === 2 && es > 1) { const u = new Uint8Array(c.length), ne = floor(c.length / es); for (let k = 0; k < ne; k++) for (let q = 0; q < es; q++) u[k * es + q] = c[q * ne + k]; c = u; } }
              const idx = new Array(rank).fill(0), rowLen = min(csz[rank - 1], shape[rank - 1] - off[rank - 1]); if (rowLen <= 0) continue;
              for (; ;) { let inside = true, dst = 0, src = 0; for (let k = 0; k < rank - 1; k++) { if (off[k] + idx[k] >= shape[k]) inside = false; dst = dst * shape[k] + off[k] + idx[k]; src = src * csz[k] + idx[k]; } if (rank === 1) { dst = 0; src = 0; } if (inside) out.set(c.subarray(src * csz[rank - 1] * es, src * csz[rank - 1] * es + rowLen * es), (dst * shape[rank - 1] + off[rank - 1]) * es); let k = rank - 2; while (k >= 0 && ++idx[k] >= csz[k]) idx[k--] = 0; if (k < 0) break; } } };
          node(bta); return readTyped(out, 0, count, dt.name, !dt.be); }; }
      } else d.why = cls === 2 ? `chunked layout version ${lv} (new chunk indexing) is not supported` : `${d.layout} storage is not supported`;
      datasets.push(d); return;
    }
    groups.push(path || '/');
    if (st) { const bta = rd(st.o, so), ha = A(rd(st.o + so, so)); if (latin(b, ha, 4) !== 'HEAP') return; const hd = A(rd(ha + 8 + 2 * sl, so));
      const nodeG = a => { a = A(a); if (undef(a)) return; const s4 = latin(b, a, 4); if (s4 === 'TREE') { const n = v.getUint16(a + 6, true); for (let i = 0; i < n; i++) nodeG(rd(a + 8 + 2 * so + sl + i * (sl + so), so)); } else if (s4 === 'SNOD') { const n = v.getUint16(a + 6, true); for (let i = 0; i < n; i++) { const e = a + 8 + i * (2 * so + 24), no = rd(e, so); let nm = ''; for (let k = hd + no; b[k] && k < b.length; k++) nm += String.fromCharCode(b[k]); visit(rd(e + so, so), path + '/' + nm, depth + 1); } } };
      nodeG(bta); }
    for (const m of links) { const fl = b[m.o + 1]; let p = m.o + 2, lt = 0; if (fl & 8) lt = b[p++]; if (fl & 4) p += 8; if (fl & 16) p++; const ls = 1 << (fl & 3), nl = rd(p, ls); p += ls; const nm = toText(b.subarray(p, p + nl)); p += nl; if (lt === 0) visit(rd(p, so), path + '/' + nm, depth + 1); }
    const li = find(0x02); if (li && !links.length && !st) { const fl = b[li.o + 1], p = li.o + 2 + (fl & 1 ? 8 : 0); if (!undef(A(rd(p, so)))) notes.push(`group ${path || '/'} uses dense link storage (fractal heap): its members are not listed`); }
  };
  if (rootBT !== null && undef(A(rootOH))) fail('HDF5', 'root group object header address is invalid'); visit(rootOH, '', 0);
  return { superblock: sbv, offsetSize: so, datasets, groups, notes };
}
export function parseHDF5(input, { dataset, threshold } = {}) {
  const h = h5Open(toBytes(input)), inv = h.datasets.map(d => ({ path: d.path, shape: d.shape, dtype: d.dtype, layout: d.layout, filters: d.filterNames, readable: d.readable, ...(d.why ? { notRead: d.why } : {}), ...(Object.keys(d.attrs).length ? { attrs: d.attrs } : {}) }));
  const meta = { superblockVersion: h.superblock, groups: h.groups.slice(0, 200), datasets: inv.slice(0, 400), nDatasets: h.datasets.length }, notes = [...h.notes], lost = [];
  const cand = h.datasets.filter(d => d.readable && d.shape.filter(n => n > 1).length >= 2 && d.shape.filter(n => n > 1).length <= 3), size = d => d.shape.reduce((x, y) => x * y, 1), best = l => l.reduce((a, d) => !a || size(d) > size(a) ? d : a, null), pick = dataset ? h.datasets.find(d => d.path === dataset || d.path === '/' + dataset) : (best(cand.filter(d => d.shape.filter(n => n > 1).length === 3)) || best(cand));
  if (dataset && !pick) fail('HDF5', `dataset "${dataset}" not found; available: ${h.datasets.map(d => d.path).slice(0, 20).join(', ')}`);
  if (!pick || !pick.readable) { const why = h.datasets.filter(d => !d.readable).slice(0, 6).map(d => `${d.path}: ${d.why}`); return { kind: null, meta, notes: [...notes, h.datasets.length ? `no 2-D/3-D numeric dataset could be read${why.length ? ' (' + why.join('; ') + ')' : ''}` : 'no datasets found in the groups that could be traversed', 'inventory only (groups, datasets, shapes, types, filters)'], lost: ['array data not imported'] }; }
  const sh = pick.shape.filter(n => n > 1), [nz, ny, nx] = sh.length === 3 ? sh : [1, sh[0], sh[1]], a = pick.attrs, spx = Array.isArray(a.spacing_xyz) ? a.spacing_xyz : null, vs = typeof a.voxel_size === 'number' ? a.voxel_size : Array.isArray(a.element_size_um) ? a.element_size_um[2] : spx ? spx[0] : 1;
  const g = scalarToVoxels(pick.read(), nx, ny, nz, { voxel: vs, units: typeof a.units === 'string' ? a.units : Array.isArray(a.element_size_um) ? 'µm' : '', spacing: spx || undefined, origin: Array.isArray(a.origin_xyz) ? a.origin_xyz : undefined, threshold });
  meta.dataset = pick.path; notes.push(`dataset ${pick.path} (${pick.dtype}, shape ${pick.shape.join(' × ')}, ${pick.layout}${pick.filterNames.length ? ', ' + pick.filterNames.join('+') : ''}) read as (z, y, x) in C order — confirm the axis order for your data`);
  if (!('voxel_size' in a) && !spx && !a.element_size_um) notes.push('no voxel-size attribute found: voxel size 1 (unitless)'); if (h.datasets.length > 1) lost.push(`${h.datasets.length - 1} other datasets listed in meta.datasets were not imported (the largest 3-D dataset is chosen by default; select another with the dataset option)`);
  return Object.assign(g, { meta, notes: [...notes, ...g.notes], lost });
}

// ---------------------------------------------------------------- NRRD, MetaImage, NIfTI, MRC, RAW
const NRRDT = { 'signed char': 'int8', int8: 'int8', int8_t: 'int8', uchar: 'uint8', 'unsigned char': 'uint8', uint8: 'uint8', uint8_t: 'uint8', short: 'int16', 'short int': 'int16', 'signed short': 'int16', int16: 'int16', int16_t: 'int16', ushort: 'uint16', 'unsigned short': 'uint16', uint16: 'uint16', uint16_t: 'uint16', int: 'int32', 'signed int': 'int32', int32: 'int32', int32_t: 'int32', uint: 'uint32', 'unsigned int': 'uint32', uint32: 'uint32', uint32_t: 'uint32', longlong: 'int64', 'long long': 'int64', int64: 'int64', int64_t: 'int64', ulonglong: 'uint64', 'unsigned long long': 'uint64', uint64: 'uint64', uint64_t: 'uint64', float: 'float32', double: 'float64' };
const splitHead = (b, re) => { const s = latin(b, 0, min(b.length, 1 << 16)), m = re.exec(s); return m ? { text: s.slice(0, m.index), end: m.index + m[0].length } : null; };
export function toNRRD(vox) {
  if (!vox || vox.kind !== 'voxels') fail('NRRD', 'needs a voxels geometry'); const sp = vox.spacing || [vox.voxel || 1, vox.voxel || 1, vox.voxel || 1];
  const h = ENC.encode(`NRRD0004\n# units: ${vox.units || 'unspecified'}\ntype: uint8\ndimension: 3\nsizes: ${vox.nx} ${vox.ny} ${vox.nz}\nspacings: ${sp.map(N9).join(' ')}\nencoding: raw\n\n`), o = new Uint8Array(h.length + vox.data.length); o.set(h); o.set(vox.data, h.length); return o;
}
export function parseNRRD(input, { raw, threshold } = {}) {
  const b = toBytes(input); if (latin(b, 0, 4) !== 'NRRD') fail('NRRD', 'missing NRRD magic'); const hd = splitHead(b, /\r?\n\r?\n/) || { text: latin(b, 0, b.length), end: b.length }, f = {};
  for (const l of lines(hd.text).slice(1)) { if (l[0] === '#') { const u = /units?:\s*(\S+)/i.exec(l); if (u) f._units = u[1]; continue; } const m = /^([^:]+):=?\s*(.*)$/.exec(l); if (m) f[m[1].trim().toLowerCase()] = m[2].trim(); }
  const dim = +f.dimension, sizes = nums(f.sizes || ''), type = NRRDT[(f.type || '').toLowerCase()], encn = (f.encoding || 'raw').toLowerCase(), le = (f.endian || 'little').toLowerCase() !== 'big';
  if (!type) fail('NRRD', `unsupported sample type "${f.type}"`); if (!(dim >= 2 && dim <= 3) || sizes.length !== dim) fail('NRRD', `dimension ${f.dimension} with sizes "${f.sizes}" is not a 2-D or 3-D scalar volume`);
  let d = b.subarray(hd.end); const df = f['data file'] || f.datafile;
  if (df) { if (!raw) fail('NRRD', `detached header: the samples are in "${df}", which must be supplied together with the header`); d = toBytes(raw); }
  const [nx, ny, nz = 1] = sizes, N = nx * ny * nz; let a;
  if (encn === 'raw') a = readTyped(d, +f['byte skip'] > 0 ? +f['byte skip'] : 0, N, type, le); else if (encn === 'gzip' || encn === 'gz') a = readTyped(gunzip(d), 0, N, type, le); else if (/^(ascii|text|txt)$/.test(encn)) a = Float64Array.from(nums(toText(d))); else fail('NRRD', `encoding "${encn}" is not supported (raw, gzip, ascii)`);
  let sp = f.spacings ? nums(f.spacings.replace(/nan/gi, '1')) : null; const sdir = f['space directions']; if (!sp && sdir) sp = [...sdir.matchAll(/\(([^)]*)\)/g)].map(m => hypot(...nums(m[1])));
  const units = (f['space units'] ? (/"([^"]*)"/.exec(f['space units']) || [])[1] : null) || f._units || '', g = scalarToVoxels(a, nx, ny, nz, { voxel: sp ? sp[0] : 1, spacing: sp && sp.length >= 3 ? sp.slice(0, 3) : undefined, units: units === 'unspecified' ? '' : units, threshold });
  return Object.assign(g, { meta: { type, encoding: encn, endian: le ? 'little' : 'big', space: f.space, fields: Object.keys(f).filter(k => k[0] !== '_') }, notes: [...g.notes, ...(sp ? [] : ['no spacing in the header: voxel size 1']), ...(units ? [] : ['no unit in the header'])], lost: sdir && /[1-9]/.test(sdir.replace(/\(([^,]*),0,0\)|\(0,([^,]*),0\)|\(0,0,([^)]*)\)/g, '')) ? ['oblique space directions reduced to their lengths (orientation not applied)'] : [] });
}
const METT = { MET_UCHAR: 'uint8', MET_CHAR: 'int8', MET_USHORT: 'uint16', MET_SHORT: 'int16', MET_UINT: 'uint32', MET_INT: 'int32', MET_ULONG: 'uint32', MET_LONG: 'int32', MET_FLOAT: 'float32', MET_DOUBLE: 'float64', MET_ULONG_LONG: 'uint64', MET_LONG_LONG: 'int64' };
export function toMHA(vox) {
  if (!vox || vox.kind !== 'voxels') fail('MetaImage', 'needs a voxels geometry'); const sp = vox.spacing || [vox.voxel || 1, vox.voxel || 1, vox.voxel || 1], o = vox.origin || [0, 0, 0];
  const h = ENC.encode(`ObjectType = Image\nNDims = 3\nBinaryData = True\nBinaryDataByteOrderMSB = False\nCompressedData = False\nOffset = ${o.map(N9).join(' ')}\nElementSpacing = ${sp.map(N9).join(' ')}\nDimSize = ${vox.nx} ${vox.ny} ${vox.nz}\nElementType = MET_UCHAR\nElementDataFile = LOCAL\n`), out = new Uint8Array(h.length + vox.data.length); out.set(h); out.set(vox.data, h.length); return out;
}
export function parseMHD(input, { raw, threshold } = {}) {
  const b = toBytes(input), hd = splitHead(b, /ElementDataFile\s*=\s*([^\r\n]*)\r?\n?/); if (!hd || !/ObjectType|NDims/.test(hd.text)) fail('MetaImage', 'no "NDims … ElementDataFile" header found');
  const f = {}; for (const l of lines(hd.text)) { const m = /^\s*(\w+)\s*=\s*(.*)$/.exec(l); if (m) f[m[1]] = m[2].trim(); } const file = (/ElementDataFile\s*=\s*([^\r\n]*)/.exec(latin(b, 0, hd.end)) || [])[1].trim();
  const dims = nums(f.DimSize || ''), type = METT[f.ElementType], le = !/true/i.test(f.BinaryDataByteOrderMSB || f.ElementByteOrderMSB || 'False'); if (!type) fail('MetaImage', `unsupported ElementType "${f.ElementType}"`); if (dims.length < 2 || dims.length > 3) fail('MetaImage', `NDims = ${f.NDims} is not a 2-D or 3-D image`); if (+f.ElementNumberOfChannels > 1) fail('MetaImage', 'multi-channel images are not supported');
  let d; if (file === 'LOCAL') d = b.subarray(hd.end); else { if (!raw) fail('MetaImage', `the samples are in the companion file "${file}", which must be supplied together with the .mhd header`); d = toBytes(raw); }
  const [nx, ny, nz = 1] = dims, N = nx * ny * nz; if (/true/i.test(f.CompressedData || '')) d = unzlib(d, N * DT[type][0]);
  const sp = nums(f.ElementSpacing || f.ElementSize || ''), g = scalarToVoxels(readTyped(d, +f.HeaderSize > 0 ? +f.HeaderSize : 0, N, type, le), nx, ny, nz, { voxel: sp[0] || 1, spacing: sp.length >= 3 ? sp.slice(0, 3) : undefined, origin: nums(f.Offset || f.Position || '').slice(0, 3), units: 'mm', threshold });
  return Object.assign(g, { meta: { elementType: f.ElementType, dataFile: file, compressed: /true/i.test(f.CompressedData || '') }, notes: [...g.notes, 'MetaImage spacing is conventionally in millimetres (assumed)'], lost: f.TransformMatrix && !/^1 0 0 0 1 0 0 0 1$/.test(f.TransformMatrix) ? ['non-identity TransformMatrix (orientation) not applied'] : [] });
}
export function parseNIfTI(input, { raw, threshold } = {}) {
  let b = toBytes(input); const notes = []; if (b[0] === 0x1F && b[1] === 0x8B) { b = gunzip(b); notes.push('gzip container (.nii.gz) inflated'); }
  if (b.length < 348) fail('NIfTI', 'file shorter than the 348-byte header'); const v = view(b); let le = true; if (v.getInt32(0, true) !== 348) { if (v.getInt32(0, false) === 348) le = false; else if (v.getInt32(0, true) === 540 || v.getInt32(0, false) === 540) fail('NIfTI', 'NIfTI-2 (540-byte header) is not supported: convert to NIfTI-1'); else fail('NIfTI', 'sizeof_hdr is not 348'); }
  const magic = latin(b, 344, 3), nd = v.getInt16(40, le), dim = [1, 2, 3, 4].map(i => max(1, v.getInt16(40 + 2 * i, le))), code = v.getInt16(70, le), type = ({ 2: 'uint8', 4: 'int16', 8: 'int32', 16: 'float32', 64: 'float64', 256: 'int8', 512: 'uint16', 768: 'uint32' })[code];
  if (magic !== 'n+1' && magic !== 'ni1') fail('NIfTI', 'magic string is neither "n+1" nor "ni1"'); if (!type) fail('NIfTI', `datatype code ${code} is not supported`);
  const pix = [1, 2, 3].map(i => abs(v.getFloat32(76 + 4 * i, le)) || 1), vo = v.getFloat32(108, le), slope = v.getFloat32(112, le), inter = v.getFloat32(116, le), xu = b[123] & 7, units = ({ 1: 'm', 2: 'mm', 3: 'µm' })[xu] || '', N = dim[0] * dim[1] * dim[2];
  let d = b, off = floor(vo) || 352; if (magic === 'ni1') { if (!raw) fail('NIfTI', 'header/image pair (.hdr + .img): the .img file must be supplied together with the header'); d = toBytes(raw); off = 0; }
  let a = readTyped(d, off, N, type, le); if (slope && (slope !== 1 || inter)) { a = Float64Array.from(a, x => x * slope + inter); notes.push(`intensity scaling applied: value × ${N9(slope)} + ${N9(inter)}`); }
  const g = scalarToVoxels(a, dim[0], dim[1], dim[2], { voxel: pix[0], spacing: pix, units, threshold });
  return Object.assign(g, { meta: { datatype: type, dims: nd, shape: dim, description: latin(b, 148, 80).replace(/\0.*$/, ''), littleEndian: le }, notes: [...notes, ...g.notes, ...(units ? [] : ['spatial unit not declared in xyzt_units'])], lost: [...(dim[3] > 1 ? [`${dim[3]} time points / volumes present: only the first was imported`] : []), 'qform/sform orientation matrices not applied (voxel axes used as stored)'] });
}
export function parseMRC(input, { threshold } = {}) {
  const b = toBytes(input); if (b.length < 1024) fail('MRC', 'file shorter than the 1024-byte header'); const v = view(b), le = !(b[212] === 0x11 && b[213] === 0x11), i32 = o => v.getInt32(o, le), f32 = o => v.getFloat32(o, le);
  const nx = i32(0), ny = i32(4), nz = i32(8), mode = i32(12), type = ({ 0: 'int8', 1: 'int16', 2: 'float32', 6: 'uint16' })[mode], ext = i32(92), map = latin(b, 208, 4);
  if (!(nx > 0 && ny > 0 && nz > 0 && nx < 1e5 && ny < 1e5 && nz < 1e5)) fail('MRC', 'implausible dimensions in the header (not an MRC/CCP4 map, or wrong byte order)'); if (!type) fail('MRC', `mode ${mode} is not supported (supported: 0 int8, 1 int16, 2 float32, 6 uint16)`);
  const mx = [i32(28), i32(32), i32(36)], cell = [f32(40), f32(44), f32(48)], sp = cell.map((c, i) => c > 0 && mx[i] > 0 ? c / mx[i] : 1), cal = cell[0] > 0 && mx[0] > 0, order = [i32(64), i32(68), i32(72)];
  const g = scalarToVoxels(readTyped(b, 1024 + max(0, ext), nx * ny * nz, type, le), nx, ny, nz, { voxel: sp[0], spacing: sp, units: cal ? 'Å' : '', threshold });
  return Object.assign(g, { meta: { mode, type, mapStamp: map, extendedHeaderBytes: ext, cell_A: cell, density: { min: f32(76), max: f32(80), mean: f32(84) }, axisOrder: order, spaceGroup: i32(88) }, notes: [...g.notes, cal ? 'voxel size = cell length / grid sampling (Å)' : 'cell dimensions are zero: voxel size 1 (uncalibrated)', ...(map === 'MAP ' ? [] : ['"MAP " stamp absent (pre-2000 or non-standard MRC variant)'])], lost: order[0] === 1 && order[1] === 2 && order[2] === 3 || !order[0] ? [] : [`non-standard axis order (mapc, mapr, maps) = ${order.join(', ')} not permuted`] });
}
export const toRAW = vox => { if (!vox || vox.kind !== 'voxels') fail('RAW', 'needs a voxels geometry'); return vox.data.slice(); };
// Headerless volume. Dimensions and sample type must come from the user; a perfect-cube byte count is the only guess made.
export function parseRAW(input, { nx, ny, nz, dtype = 'uint8', littleEndian = true, headerBytes = 0, voxel = 1, units = '', threshold } = {}) {
  const b = toBytes(input), es = (DT[dtype] || fail('RAW', `unknown dtype "${dtype}"`))[0], n = floor((b.length - headerBytes) / es), notes = [];
  if (!(nx > 0 && ny > 0 && nz > 0)) { const c = round(Math.cbrt(n)); if (c > 1 && c * c * c === n) { nx = ny = nz = c; notes.push(`dimensions not given: ${n} samples form a perfect cube, so ${c} × ${c} × ${c} was assumed`); } else fail('RAW', `a headerless volume needs nx, ny, nz and dtype: the file holds ${n} ${dtype} samples, which is not a perfect cube`); }
  if (nx * ny * nz > n) fail('RAW', `${nx} × ${ny} × ${nz} ${dtype} needs ${nx * ny * nz * es} bytes but the file has ${b.length - headerBytes}`);
  if (nx * ny * nz < n) notes.push(`${(n - nx * ny * nz) * es} trailing bytes ignored`);
  const g = scalarToVoxels(readTyped(b, headerBytes, nx * ny * nz, dtype, littleEndian), nx, ny, nz, { voxel, units, threshold });
  return Object.assign(g, { meta: { dtype, littleEndian, headerBytes }, notes: [...notes, ...g.notes, 'raw volumes carry no metadata: dimensions, sample type, byte order, voxel size and x-fastest ordering are user-supplied or assumed'], lost: [] });
}

// ---------------------------------------------------------------- EBSD orientation maps (Oxford CTF, EDAX/TSL ANG)
export function parseCTF(input) {
  const L = lines(toText(input)); if (!/^Channel Text File/i.test(L[0] || '')) fail('CTF', 'first line is not "Channel Text File"');
  const hi = L.findIndex(l => /^Phase\s+X\s+Y\b/.test(l)); if (hi < 0) fail('CTF', 'column header line "Phase X Y …" not found'); const f = {}, phases = [];
  for (let i = 1; i < hi; i++) { const t = L[i].split(/\t/); if (/^Phases$/i.test(t[0])) { const n = +t[1]; for (let k = 1; k <= n && i + k < hi; k++) { const q = L[i + k].split(/\t/); phases.push({ id: k, name: q[2] || `phase ${k}`, cell: q[0], angles: q[1], laueGroup: q[3], spaceGroup: q[4] }); } i += n; } else if (t.length >= 2) f[t[0].trim()] = t[1].trim(); }
  const cols = L[hi].trim().split(/\s+/), ci = n => cols.indexOf(n), nx = +f.XCells, ny = +f.YCells, xs = +f.XStep || 1, ys = +f.YStep || xs; if (!(nx > 0 && ny > 0)) fail('CTF', 'XCells / YCells missing from the header'); if (nx * ny > 4e7) fail('CTF', `map of ${nx} × ${ny} points is too large for in-browser import`);
  const N = nx * ny, data = new Uint8Array(N), e1 = new Float32Array(N), e2 = new Float32Array(N), e3 = new Float32Array(N), q = new Float32Array(N), cP = ci('Phase'), cX = ci('X'), cY = ci('Y'), cE = ci('Euler1'), cQ = ci('BC'); let n = 0, zero = 0;
  for (let i = hi + 1; i < L.length; i++) { if (!L[i].trim()) continue; const t = L[i].trim().split(/\s+/), ix = round(+t[cX] / xs), iy = round(+t[cY] / ys); if (!(ix >= 0 && ix < nx && iy >= 0 && iy < ny)) continue; const p = ix + nx * iy; data[p] = +t[cP]; if (!data[p]) zero++; e1[p] = +t[cE]; e2[p] = +t[cE + 1]; e3[p] = +t[cE + 2]; if (cQ >= 0) q[p] = +t[cQ]; n++; }
  if (!n) fail('CTF', 'no data rows could be placed on the XCells × YCells grid');
  return { kind: 'voxels', nx, ny, nz: 1, data, voxel: xs, spacing: [xs, ys, xs], units: 'µm', orientation: { euler1: e1, euler2: e2, euler3: e3, units: 'deg', convention: 'Bunge ZXZ (φ1, Φ, φ2), acquisition frame' }, quality: { name: 'band contrast', values: q }, phases,
    meta: { project: f.Prj, xCells: nx, yCells: ny, xStep: xs, yStep: ys, kV: +f.KV || null, magnification: +f.Mag || null, indexedFraction: 1 - zero / n, points: n }, notes: ['voxel value = phase id (0 = not indexed); Euler angles in degrees kept in `orientation`', 'map row order kept as acquired (y index increases down the scan)', `step size ${xs} × ${ys} µm`], lost: ['sample/acquisition reference-frame rotations (AcqE1–3) are not applied', 'band count, error and MAD columns not imported'] };
}
export function parseANG(input) {
  const L = lines(toText(input)), f = {}, phases = []; let i = 0;
  for (; i < L.length && (L[i].trim() === '' || L[i].trim()[0] === '#'); i++) { const m = /^#\s*([\w-]+):?\s+(.*)$/.exec(L[i].trim()); if (m) { if (m[1] === 'MaterialName') phases.push({ id: phases.length + 1, name: m[2].trim() }); else if (!(m[1] in f)) f[m[1].toUpperCase()] = m[2].trim(); } }
  if (i === 0 || i >= L.length) fail('ANG', 'no "#" header block followed by data rows'); const xs = +f.XSTEP || 1, ys = +f.YSTEP || xs, hex = /hex/i.test(f.GRID || ''), rows = [];
  for (; i < L.length; i++) { const t = L[i].trim(); if (!t) continue; const a = t.split(/\s+/).map(Number); if (a.length < 5 || !Number.isFinite(a[4])) fail('ANG', `data row ${rows.length + 1} has fewer than five numeric columns`); rows.push(a); }
  let nx = +f.NCOLS_ODD || 0, ny = +f.NROWS || 0; if (!(nx > 0 && ny > 0)) { let mxx = 0, mxy = 0; for (const a of rows) { mxx = max(mxx, a[3]); mxy = max(mxy, a[4]); } nx = round(mxx / xs) + 1; ny = round(mxy / ys) + 1; }
  if (!(nx * ny <= 4e7)) fail('ANG', `map of ${nx} × ${ny} points is too large or the grid header is corrupt`);
  const N = nx * ny, data = new Uint8Array(N), e1 = new Float32Array(N), e2 = new Float32Array(N), e3 = new Float32Array(N), q = new Float32Array(N), D = 180 / PI; let bad = 0, placed = 0;
  for (const a of rows) { const iy = round(a[4] / ys), ix = hex ? floor(a[3] / xs + 1e-6) : round(a[3] / xs); if (!(ix >= 0 && ix < nx && iy >= 0 && iy < ny)) continue; const p = ix + nx * iy, ok = !(a[6] < 0); data[p] = ok ? (Number.isFinite(a[7]) ? a[7] : 0) + 1 : 0; if (!ok) bad++; e1[p] = a[0] * D; e2[p] = a[1] * D; e3[p] = a[2] * D; q[p] = a[5] || 0; placed++; }
  if (!placed) fail('ANG', 'no data rows could be placed on the scan grid');
  return { kind: 'voxels', nx, ny, nz: 1, data, voxel: xs, spacing: [xs, ys, xs], units: 'µm', orientation: { euler1: e1, euler2: e2, euler3: e3, units: 'deg', convention: 'Bunge ZXZ (φ1, Φ, φ2), converted from radians' }, quality: { name: 'image quality (IQ)', values: q }, phases,
    meta: { grid: f.GRID || 'SqrGrid', xStep: xs, yStep: ys, nCols: nx, nRows: ny, points: placed, indexedFraction: 1 - bad / placed, operator: f.OPERATOR, sampleId: f.SAMPLEID }, notes: ['voxel value = phase index + 1 (0 = confidence index < 0, i.e. not indexed); Euler angles converted to degrees in `orientation`', `step size ${xs} × ${ys} µm`, ...(hex ? ['hexagonal scan grid resampled by nearest column onto a rectangular array (odd rows shifted by half a step)'] : [])], lost: ['confidence index, fit and detector-signal columns not imported'] };
}

// ---------------------------------------------------------------- graphs: GraphML, GEXF, GML, node/edge CSV
function graphOut(ids, pos, edges, radii, meta, fmt) {
  const idx = new Map(ids.map((id, i) => [String(id), i])), notes = [], lost = []; let have = pos.filter(p => p && Number.isFinite(p[0]) && Number.isFinite(p[1])).length, bad = 0;
  if (!ids.length) fail(fmt, 'no nodes found');
  const nodes = pos.map((p, i) => p && Number.isFinite(p[0]) && Number.isFinite(p[1]) ? [p[0], p[1], Number.isFinite(p[2]) ? p[2] : 0] : null), synthetic = have < ids.length;
  if (synthetic) { const R = max(1, ids.length / (2 * PI)); nodes.forEach((p, i) => { if (!p) nodes[i] = [R * cos(2 * PI * i / ids.length), R * sin(2 * PI * i / ids.length), 0]; }); notes.push(`${ids.length - have} of ${ids.length} nodes carry no coordinates: they were placed on a circle. Only the connectivity is meaningful for those nodes`); }
  const E = [], rr = []; edges.forEach((e, k) => { const a = idx.get(String(e[0])), b = idx.get(String(e[1])); if (a === undefined || b === undefined) { bad++; return; } if (a !== b) { E.push([a, b]); rr.push(radii ? radii[k] : undefined); } });
  if (bad) lost.push(`${bad} edges reference undefined nodes and were dropped`);
  const g = { kind: 'lattice', nodes, edges: E, radius: 0, units: '', meta: { ...meta, nodes: ids.length, edges: E.length }, coordinatesSynthetic: synthetic, notes: [...notes, 'graph files carry no length unit and no strut section: radius 0, units unspecified'], lost };
  if (rr.some(v => Number.isFinite(v))) { g.radii = rr.map(v => Number.isFinite(v) ? v : 0); g.radius = g.radii.reduce((a, b) => a + b, 0) / g.radii.length; g.notes.push('per-edge radius attribute imported as strut radii'); }
  return g;
}
export function parseGraphML(input) {
  const t = toText(input); if (!/<graphml\b/i.test(t)) fail('GraphML', 'no <graphml> root element');
  const keys = {}; for (const m of t.matchAll(/<key\b([^>]*)\/?>/g)) keys[attr(m[1], 'id')] = (attr(m[1], 'attr.name') || attr(m[1], 'id') || '').toLowerCase();
  const dataOf = body => { const o = {}; for (const m of body.matchAll(/<data\b([^>]*)>([^<]*)<\/data>/g)) o[keys[attr(m[1], 'key')] || attr(m[1], 'key')] = xmlText(m[2].trim()); return o; };
  const ids = [], pos = [], edges = [], radii = [];
  for (const m of t.matchAll(/<node\b([^>]*?)(?:\/>|>([\s\S]*?)<\/node>)/g)) { const d = dataOf(m[2] || ''), yx = /<y:Geometry\b([^>]*)\/?>/.exec(m[2] || ''); ids.push(attr(m[1], 'id')); pos.push(yx ? [+attr(yx[1], 'x'), +attr(yx[1], 'y'), 0] : [parseFloat(d.x), parseFloat(d.y), parseFloat(d.z)]); }
  for (const m of t.matchAll(/<edge\b([^>]*?)(?:\/>|>([\s\S]*?)<\/edge>)/g)) { const d = dataOf(m[2] || ''); edges.push([attr(m[1], 'source'), attr(m[1], 'target')]); radii.push(parseFloat(d.radius ?? (d.diameter !== undefined ? d.diameter / 2 : NaN))); }
  return graphOut(ids, pos, edges, radii, { attributes: Object.values(keys), directed: /edgedefault\s*=\s*"directed"/.test(t) }, 'GraphML');
}
export function parseGEXF(input) {
  const t = toText(input); if (!/<gexf\b/i.test(t)) fail('GEXF', 'no <gexf> root element'); const ids = [], pos = [], edges = [];
  for (const m of t.matchAll(/<node\b([^>]*?)(?:\/>|>([\s\S]*?)<\/node>)/g)) { const p = /<(?:viz:)?position\b([^>]*)\/?>/.exec(m[2] || ''); ids.push(attr(m[1], 'id')); pos.push(p ? [+attr(p[1], 'x'), +attr(p[1], 'y'), +attr(p[1], 'z') || 0] : null); }
  for (const m of t.matchAll(/<edge\b([^>]*?)\/?>/g)) edges.push([attr(m[1], 'source'), attr(m[1], 'target')]);
  return graphOut(ids, pos, edges, null, { version: attr((/<gexf\b[^>]*>/.exec(t) || [''])[0], 'version') }, 'GEXF');
}
export function parseGML(input) {
  const t = toText(input); if (!/\bgraph\s*\[/.test(t)) fail('GML', 'no "graph [" block'); const ids = [], pos = [], edges = [], num1 = (s, k) => { const m = new RegExp('\\b' + k + '\\s+(-?[\\d.eE+-]+)').exec(s); return m ? +m[1] : NaN; };
  for (const m of t.matchAll(/\bnode\s*\[((?:[^\[\]]|\[[^\[\]]*\])*)\]/g)) { ids.push(num1(m[1], 'id')); pos.push([num1(m[1], 'x'), num1(m[1], 'y'), num1(m[1], 'z')]); }
  for (const m of t.matchAll(/\bedge\s*\[((?:[^\[\]]|\[[^\[\]]*\])*)\]/g)) edges.push([num1(m[1], 'source'), num1(m[1], 'target')]);
  return graphOut(ids, pos, edges, null, { directed: /\bdirected\s+1/.test(t) }, 'GML');
}
// Delimited table → {header, rows}; comment lines (#, %, //) and blank lines skipped; delimiter inferred.
export function parseTable(text) {
  const L = lines(text).map(l => l.trim()).filter(l => l && !/^(#|%|\/\/)/.test(l)); if (!L.length) return { header: null, rows: [] };
  const delim = [',', ';', '\t'].find(d => L[0].includes(d)), split = l => (delim ? l.split(delim) : l.split(/\s+/)).map(s => s.trim().replace(/^"(.*)"$/, '$1')), first = split(L[0]), hasH = first.some(s => s !== '' && !Number.isFinite(Number(s)));
  return { header: hasH ? first.map(s => s.toLowerCase()) : null, rows: L.slice(hasH ? 1 : 0).map(split) };
}
export const EDGE_HEADER = /^\s*["']?(source|from|src|node1|n1|start)\b["']?\s*[,;\t ]\s*["']?(target|to|dst|node2|n2|end)\b/i; // header line of an edge table
// Node/edge tables. Accepted: (a) one text with a nodes table (id,x,y[,z]) and an edges table (source,target[,radius])
// separated by a blank line or a second header; (b) the two tables as separate arguments; (c) a segment table
// x1,y1,z1,x2,y2,z2[,radius].
export function parseNodeEdgeCSV(input, edgesInput) {
  const t = toText(input); let nt = t, et = edgesInput !== undefined ? toText(edgesInput) : null;
  if (et === null) { const L = lines(t), k = L.findIndex((l, i) => i > 0 && EDGE_HEADER.test(l)); if (k > 0) { nt = L.slice(0, k).join('\n'); et = L.slice(k).join('\n'); } }
  const A = parseTable(nt), col = (h, names, d) => { if (!h) return d; const i = h.findIndex(c => names.includes(c)); return i; };
  if (et === null) {
    const h = A.header, x1 = col(h, ['x1', 'xa', 'x_start'], 0); if (h && x1 < 0) fail('node–edge CSV', 'expected a nodes table followed by an edges table (header source,target), or a segment table with columns x1,y1,z1,x2,y2,z2');
    if (!A.rows.length || A.rows[0].length < 4) fail('node–edge CSV', 'a segment table needs at least x1,y1,x2,y2 per row');
    const six = A.rows[0].length >= 6, key = new Map(), nodes = [], edges = [], radii = [], add = p => { const k = p.map(v => N9(v)).join(','); let i = key.get(k); if (i === undefined) { i = nodes.length; nodes.push(p); key.set(k, i); } return i; };
    for (const r of A.rows) { const a = r.map(Number), p = six ? [a[0], a[1], a[2]] : [a[0], a[1], 0], q = six ? [a[3], a[4], a[5]] : [a[2], a[3], 0]; if (![...p, ...q].every(Number.isFinite)) fail('node–edge CSV', `non-numeric coordinate in row "${r.join(',')}"`); edges.push([add(p), add(q)]); radii.push(a[six ? 6 : 4]); }
    const g = { kind: 'lattice', nodes, edges, radius: 0, units: '', meta: { layout: 'segment table', segments: edges.length }, notes: ['segment end points merged into shared nodes (9 significant digits)', 'CSV carries no length unit'], lost: [] }; if (radii.every(Number.isFinite)) { g.radii = radii; g.radius = radii.reduce((s, v) => s + v, 0) / radii.length; } return g;
  }
  const B = parseTable(et), h = A.header, ci = col(h, ['id', 'node', 'label', 'name', 'n'], 0), cx = col(h, ['x'], 1), cy = col(h, ['y'], 2), cz = col(h, ['z'], 3), hb = B.header, cs = col(hb, ['source', 'from', 'n1', 'node1', 'start', 'i'], 0), ct = col(hb, ['target', 'to', 'n2', 'node2', 'end', 'j'], 1), cr = hb ? hb.findIndex(c => ['radius', 'r'].includes(c)) : 2, cd = hb ? hb.indexOf('diameter') : -1;
  if (cx < 0 || cy < 0) fail('node–edge CSV', 'nodes table needs x and y columns'); if (cs < 0 || ct < 0) fail('node–edge CSV', 'edges table needs source and target columns');
  return graphOut(A.rows.map((r, i) => ci >= 0 && h ? r[ci] : (h ? i : r[0])), A.rows.map(r => [+r[cx], +r[cy], cz >= 0 && r[cz] !== undefined ? +r[cz] : 0]), B.rows.map(r => [r[cs], r[ct]]), B.rows.map(r => cr >= 0 && r[cr] !== undefined && r[cr] !== '' ? +r[cr] : cd >= 0 ? +r[cd] / 2 : NaN), { layout: 'nodes + edges tables' }, 'node–edge CSV');
}

// ---------------------------------------------------------------- point clouds (XYZ, PTS, PTX, PCD, CSV)
export function parsePointCloud(input, { maxPoints = 300000 } = {}) {
  const b = toBytes(input), head = latin(b, 0, min(b.length, 4096)), meta = {}, notes = [], lost = []; let pts = [], total = 0, fmt = 'xyz';
  if (/^#\s*\.PCD|^VERSION\b/m.test(head) && /\bDATA\s+(ascii|binary|binary_compressed)/.test(head)) {
    fmt = 'pcd'; const dm = /\bDATA\s+(\w+)\r?\n/.exec(head), hl = lines(head.slice(0, dm.index)), get = k => (hl.find(l => l.startsWith(k + ' ')) || '').split(/\s+/).slice(1), fields = get('FIELDS'), size = get('SIZE').map(Number), type = get('TYPE'), cnt = get('COUNT').map(Number), n = +get('POINTS')[0] || (+get('WIDTH')[0] * +get('HEIGHT')[0]);
    const fi = ['x', 'y', 'z'].map(k => fields.indexOf(k)); if (fi[0] < 0 || fi[1] < 0) fail('PCD', 'FIELDS must contain x and y'); meta.fields = fields; total = n;
    if (dm[1] === 'ascii') { const L = lines(toText(b.subarray(dm.index + dm[0].length))); const st = max(1, Math.ceil(n / maxPoints)); for (let i = 0; i < L.length && i < n; i += st) { const a = L[i].trim().split(/\s+/); if (a.length > fi[1]) pts.push([+a[fi[0]], +a[fi[1]], fi[2] >= 0 ? +a[fi[2]] : 0]); } }
    else if (dm[1] === 'binary') { const offs = []; let o = 0; fields.forEach((_, i) => { offs.push(o); o += size[i] * (cnt[i] || 1); }); const v = view(b), base = dm.index + dm[0].length, rd = (p, i) => { if (i < 0) return 0; const q = p + offs[i]; return type[i] === 'F' ? (size[i] === 8 ? v.getFloat64(q, true) : v.getFloat32(q, true)) : type[i] === 'U' ? (size[i] === 4 ? v.getUint32(q, true) : size[i] === 2 ? v.getUint16(q, true) : b[q]) : (size[i] === 4 ? v.getInt32(q, true) : size[i] === 2 ? v.getInt16(q, true) : v.getInt8(q)); }; if (base + n * o > b.length) fail('PCD', 'binary payload shorter than POINTS × point size'); const st = max(1, Math.ceil(n / maxPoints)); for (let i = 0; i < n; i += st) pts.push([rd(base + i * o, fi[0]), rd(base + i * o, fi[1]), rd(base + i * o, fi[2])]); }
    else fail('PCD', 'binary_compressed (LZF) point data are not supported: save as ascii or binary');
    if (fields.some(f => /rgb|intensity|normal/.test(f))) lost.push('colour, intensity and normal fields not imported');
  } else {
    const L = lines(toText(b)); let i = 0, hdr = null; while (i < L.length && (!L[i].trim() || /^(#|%|\/\/)/.test(L[i].trim()))) i++;
    const f0 = (L[i] || '').trim().split(/[\s,;]+/);
    if (f0.length === 1 && /^\d+$/.test(f0[0]) && L[i + 1] !== undefined) { const f1 = L[i + 1].trim().split(/[\s,;]+/); if (f1.length === 1 && /^\d+$/.test(f1[0])) { fmt = 'ptx'; meta.columns = +f0[0]; meta.rows = +f1[0]; i += 10; lost.push('scanner pose and registration matrix in the PTX header are not applied (points kept in scanner coordinates)'); } else { fmt = 'pts'; meta.declaredPoints = +f0[0]; i++; } }
    else if (f0.some(s => s !== '' && !Number.isFinite(Number(s)))) { hdr = f0.map(s => s.toLowerCase().replace(/"/g, '')); i++; fmt = 'csv'; }
    const cx = hdr ? hdr.indexOf('x') : 0, cy = hdr ? hdr.indexOf('y') : 1, cz = hdr ? hdr.indexOf('z') : 2; if (cx < 0 || cy < 0) fail('point cloud', `header "${hdr.join(',')}" has no x and y columns`);
    const rest = L.length - i, st = max(1, Math.ceil(rest / maxPoints)); let extra = 0, skipped = 0;
    for (let k = i; k < L.length; k++) { const l = L[k].trim(); if (!l) continue; total++; if ((k - i) % st) continue; const a = l.split(/[\s,;]+/), p = [+a[cx], +a[cy], cz >= 0 && a[cz] !== undefined ? +a[cz] : 0]; if (!Number.isFinite(p[0] + p[1] + p[2])) { skipped++; continue; } if (fmt === 'ptx' && p[0] === 0 && p[1] === 0 && p[2] === 0) { skipped++; continue; } if (a.length > 3) extra = max(extra, a.length - 3); pts.push(p); }
    if (extra) lost.push(`${extra} extra column(s) per point (intensity / colour / normals) not imported`); if (skipped) notes.push(`${skipped} empty or non-numeric rows skipped`);
  }
  if (!pts.length) fail('point cloud', 'no numeric x y z rows found');
  meta.format = fmt; meta.points = total; const g = { kind: 'lattice', nodes: pts, edges: [], radius: 0, units: '', pointCloud: true, meta, notes: [...notes, 'point clouds carry no surface: reconstruction (meshing) is needed before use as a solid', 'no length unit stored: coordinates taken as written'], lost };
  if (pts.length < total) g.truncated = `uniformly subsampled to ${pts.length} of ${total} points for display`;
  return g;
}

// ---------------------------------------------------------------- safe expression evaluation (no eval)
const TOK = /\s*(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|[A-Za-z_$][\w$]*|"(?:[^"\\]|\\.)*"|==|!=|<=|>=|&&|\|\||[-+*\/%^()\[\]{};,=:.<>!?#])/y;
function tokenize(src, fmt) { const s = src.replace(/\/\*[\s\S]*?\*\//g, ' ').replace(/\/\/[^\n]*/g, ' '), out = []; TOK.lastIndex = 0; let m, p = 0; while (p < s.length) { TOK.lastIndex = p; m = TOK.exec(s); if (!m) { if (/^\s*$/.test(s.slice(p))) break; fail(fmt, `unexpected character "${s.slice(p).trim()[0]}" near "${s.slice(p, p + 20).trim()}"`); } out.push(m[1]); p = TOK.lastIndex; } return out; }
const FN = { sqrt: Math.sqrt, abs: Math.abs, exp: Math.exp, ln: Math.log, log: Math.log10, floor: Math.floor, ceil: Math.ceil, round: Math.round, min: Math.min, max: Math.max, pow: Math.pow, sign: Math.sign, hypot: Math.hypot, fabs: Math.abs, fmod: (a, b) => a % b };
// Arithmetic over numbers and vectors with variables and a fixed function table. `deg` selects degree-based trigonometry.
function makeEval(tk, env, { deg = false, fmt = 'expression' } = {}) {
  const st = { i: 0 }, k = deg ? PI / 180 : 1, trig = { sin: x => sin(x * k), cos: x => cos(x * k), tan: x => Math.tan(x * k), asin: x => Math.asin(x) / k, acos: x => Math.acos(x) / k, atan: x => Math.atan(x) / k, atan2: (a, b) => Math.atan2(a, b) / k };
  const peek = () => tk[st.i], next = () => tk[st.i++], expect = c => { if (tk[st.i] !== c) fail(fmt, `expected "${c}" but found "${tk[st.i] ?? 'end of file'}"`); st.i++; };
  const vop = (a, b, f) => Array.isArray(a) ? (Array.isArray(b) ? a.map((x, i) => vop(x, b[i] ?? 0, f)) : a.map(x => vop(x, b, f))) : Array.isArray(b) ? b.map(x => vop(a, x, f)) : f(a, b);
  const prim = () => {
    const t = next(); if (t === undefined) fail(fmt, 'unexpected end of expression');
    if (t === '(') { const v = expr(); expect(')'); return v; }
    if (t === '[' || t === '{') { const end = t === '[' ? ']' : '}', a = []; while (peek() !== end) { a.push(expr()); if (peek() === ':') fail(fmt, 'range expressions are not supported'); if (peek() === ',') next(); else if (peek() !== end) fail(fmt, `expected "," or "${end}"`); } next(); return a; }
    if (t === '-') return vop(unary(), 0, a => -a); if (t === '+') return unary(); if (t === '!') return unary() ? 0 : 1;
    if (/^[\d.]/.test(t)) return Number(t); if (t[0] === '"') return t.slice(1, -1);
    if (/^[A-Za-z_$]/.test(t)) {
      const lc = t.toLowerCase();
      if (peek() === '(' || peek() === '[' && (FN[lc] || trig[lc]) && !(t in env)) { const close = next() === '(' ? ')' : ']', a = []; while (peek() !== close) { a.push(expr()); if (peek() === ',') next(); } next(); const f = trig[lc] || FN[lc]; if (!f) fail(fmt, `function "${t}" is not supported`); return f(...a.flat()); }
      if (t in env) { let v = env[t]; while (peek() === '[' && Array.isArray(v)) { next(); const ix = expr(); expect(']'); v = v[ix]; } return v; }
      if (lc === 'pi') return PI; if (lc === 'true') return 1; if (lc === 'false') return 0; if (lc === 'undef') return undefined;
      fail(fmt, `unknown variable "${t}"`);
    }
    fail(fmt, `unexpected token "${t}"`);
  };
  const unary = () => { const v = prim(); if (peek() === '^') { next(); return Math.pow(v, unary()); } return v; };
  const term = () => { let v = unary(); for (; ;) { const o = peek(); if (o === '*') { next(); const r = unary(); v = Array.isArray(v) && Array.isArray(r) ? v.reduce((s, x, i) => s + x * r[i], 0) : vop(v, r, (a, b) => a * b); } else if (o === '/') { next(); v = vop(v, unary(), (a, b) => a / b); } else if (o === '%') { next(); v = vop(v, unary(), (a, b) => a % b); } else return v; } };
  const sum = () => { let v = term(); for (; ;) { const o = peek(); if (o === '+') { next(); v = vop(v, term(), (a, b) => a + b); } else if (o === '-') { next(); v = vop(v, term(), (a, b) => a - b); } else return v; } };
  const cmp = () => { let v = sum(); for (; ;) { const o = peek(); if (['<', '>', '<=', '>=', '==', '!='].includes(o)) { next(); const r = sum(); v = (o === '<' ? v < r : o === '>' ? v > r : o === '<=' ? v <= r : o === '>=' ? v >= r : o === '==' ? v === r : v !== r) ? 1 : 0; } else return v; } };
  const expr = () => { let v = cmp(); while (peek() === '&&' || peek() === '||') { const o = next(), r = cmp(); v = (o === '&&' ? v && r : v || r) ? 1 : 0; } if (peek() === '?') { next(); const a = expr(); expect(':'); const c = expr(); return v ? a : c; } return v; };
  return { st, peek, next, expect, expr };
}

// ---------------------------------------------------------------- OpenSCAD (parsed, never executed)
// Subset: numeric/vector variables, cube, sphere, cylinder, translate, rotate, scale, mirror, union, difference,
// intersection (plus group/render/color as pass-through) and explicit polyhedron(). The CSG tree is evaluated as a signed-distance field,
// voxelised and iso-surfaced. Unsupported modules are reported and contribute no geometry.
export function parseSCAD(input, { resolution = 64 } = {}) {
  const src = toText(input), tk = tokenize(src, 'OpenSCAD'), env = { $fn: 0, $fa: 12, $fs: 2 }, ev = makeEval(tk, env, { deg: true, fmt: 'OpenSCAD' }), { st, peek, next, expect } = ev, unsupported = {}, prims = {}, polys = [], PASS = ['union', 'group', 'render', 'color', 'colour'];
  const skipBalanced = () => { let d = 0; for (; st.i < tk.length;) { const t = next(); if (t === '{' || t === '(' || t === '[') d++; else if (t === '}' || t === ')' || t === ']') { d--; if (d === 0 && t === '}') return; } else if (t === ';' && d === 0) return; } };
  const vec3 = (v, d) => v === undefined ? [d, d, d] : Array.isArray(v) ? [v[0] ?? d, v[1] ?? d, v[2] ?? d].map(Number) : [Number(v), Number(v), Number(v)];
  const box8 = (lo, hi, f) => { const o = [[Infinity, Infinity, Infinity], [-Infinity, -Infinity, -Infinity]]; for (let m = 0; m < 8; m++) { const p = f([m & 1 ? hi[0] : lo[0], m & 2 ? hi[1] : lo[1], m & 4 ? hi[2] : lo[2]]); for (let a = 0; a < 3; a++) { o[0][a] = min(o[0][a], p[a]); o[1][a] = max(o[1][a], p[a]); } } return o; };
  const hull2 = (a, c) => !a ? c : !c ? a : { lo: a.lo.map((v, i) => min(v, c.lo[i])), hi: a.hi.map((v, i) => max(v, c.hi[i])) };
  const U = list => { const L = list.filter(Boolean); if (!L.length) return null; return { f: L.length === 1 ? L[0].f : sdf.union(...L.map(n => n.f)), ...L.reduce((a, c) => hull2(a, c), null) }; };
  function statement() {
    let t = peek(); while (['*', '!', '#', '%'].includes(t)) { const mod = next(); if (mod === '*') { skipBalanced(); return null; } t = peek(); }
    if (t === ';') { next(); return null; } if (t === '{') { next(); const c = block('}'); return U(c); }
    if (t === 'module' || t === 'function') { unsupported[`${t} definitions`] = (unsupported[`${t} definitions`] || 0) + 1; skipBalanced(); return null; }
    if (['for', 'if', 'else', 'let', 'include', 'use', 'intersection_for', 'assert', 'echo'].includes(t)) { unsupported[t] = (unsupported[t] || 0) + 1; next(); if (t === 'include' || t === 'use') { while (st.i < tk.length && peek() !== '>') next(); next(); return null; } skipBalanced(); return null; }
    if (!/^[A-Za-z_$]/.test(t || '')) fail('OpenSCAD', `unexpected token "${t}"`);
    if (tk[st.i + 1] === '=') { const name = next(); next(); env[name] = ev.expr(); expect(';'); return null; }
    const op = next(), pos = [], named = {}; expect('(');
    while (peek() !== ')') { if (/^[A-Za-z_$]/.test(peek()) && tk[st.i + 1] === '=') { const n = next(); next(); named[n] = ev.expr(); } else pos.push(ev.expr()); if (peek() === ',') next(); else if (peek() !== ')') fail('OpenSCAD', `expected "," or ")" in the arguments of ${op}`); }
    next(); const A = (i, n) => named[n] !== undefined ? named[n] : pos[i];
    let kids = []; if (peek() === ';') next(); else if (peek() === '{') { next(); kids = block('}'); } else { const c = statement(); if (c) kids = [c]; }
    if (op === 'cube') { prims.cube = (prims.cube || 0) + 1; const s = vec3(A(0, 'size'), 1).map(v => abs(v) || 1e-9), c = A(1, 'center') ? [0, 0, 0] : s.map(v => v / 2); return { f: sdf.box(s.map(v => v / 2), c), lo: c.map((v, i) => v - s[i] / 2), hi: c.map((v, i) => v + s[i] / 2) }; }
    if (op === 'sphere') { prims.sphere = (prims.sphere || 0) + 1; const r = abs(named.d !== undefined ? named.d / 2 : A(0, 'r') ?? 1); return { f: sdf.sphere(r), lo: [-r, -r, -r], hi: [r, r, r] }; }
    if (op === 'cylinder') { prims.cylinder = (prims.cylinder || 0) + 1; const h = abs(A(0, 'h') ?? 1), r = named.d !== undefined ? named.d / 2 : named.r, r1 = abs(named.d1 !== undefined ? named.d1 / 2 : named.r1 ?? r ?? pos[1] ?? 1), r2 = abs(named.d2 !== undefined ? named.d2 / 2 : named.r2 ?? r ?? pos[2] ?? pos[1] ?? 1), cz = (named.center ?? pos[3]) ? 0 : h / 2, R = max(r1, r2); return { f: abs(r1 - r2) < 1e-12 ? sdf.cylinder(r1, h, [0, 0, cz]) : sdf.cone(r1, r2, h, [0, 0, cz]), lo: [-R, -R, cz - h / 2], hi: [R, R, cz + h / 2] }; }
    if (op === 'polyhedron') { const P = A(0, 'points'), Fc = A(1, 'faces') ?? named.triangles; if (!Array.isArray(P) || !Array.isArray(Fc)) fail('OpenSCAD', 'polyhedron needs points and faces lists'); polys.push({ P, Fc }); return null; }
    const child = U(kids);
    if (PASS.includes(op)) return child;
    if (op === 'difference') { const L = kids.filter(Boolean); if (!L.length) return null; return { f: L.length > 1 ? sdf.subtract(L[0].f, ...L.slice(1).map(n => n.f)) : L[0].f, lo: L[0].lo, hi: L[0].hi }; }
    if (op === 'intersection') { const L = kids.filter(Boolean); if (!L.length) return null; return { f: sdf.intersect(...L.map(n => n.f)), lo: L[0].lo.map((v, i) => max(v, ...L.map(n => n.lo[i]))), hi: L[0].hi.map((v, i) => min(v, ...L.map(n => n.hi[i]))) }; }
    if (['translate', 'rotate', 'scale', 'mirror'].includes(op)) {
      if (!child) return null; const v = A(0, op === 'rotate' ? 'a' : 'v');
      if (op === 'translate') { const d = vec3(v, 0); return { f: sdf.translate(child.f, d), lo: child.lo.map((x, i) => x + d[i]), hi: child.hi.map((x, i) => x + d[i]) }; }
      if (op === 'scale') { const s = vec3(v, 1).map(x => x || 1e-9), bb = box8(child.lo, child.hi, p => p.map((x, i) => x * s[i])); return { f: sdf.scale(child.f, s), lo: bb[0], hi: bb[1] }; }
      if (op === 'mirror') { const n = vec3(v, 0), l = hypot(...n) || 1, u = n.map(x => x / l), refl = p => { const d = 2 * (p[0] * u[0] + p[1] * u[1] + p[2] * u[2]); return [p[0] - d * u[0], p[1] - d * u[1], p[2] - d * u[2]]; }, bb = box8(child.lo, child.hi, refl), f = child.f; return { f: (x, y, z) => { const q = refl([x, y, z]); return f(q[0], q[1], q[2]); }, lo: bb[0], hi: bb[1] }; }
      let R; const axis = A(1, 'v');
      if (Array.isArray(axis) && !Array.isArray(v)) { const l = hypot(...axis) || 1, [ux, uy, uz] = axis.map(x => x / l), a = v * PI / 180, c = cos(a), s = sin(a), C = 1 - c; R = [[c + ux * ux * C, ux * uy * C - uz * s, ux * uz * C + uy * s], [uy * ux * C + uz * s, c + uy * uy * C, uy * uz * C - ux * s], [uz * ux * C - uy * s, uz * uy * C + ux * s, c + uz * uz * C]]; }
      else R = sdf.rotationMatrix(Array.isArray(v) ? vec3(v, 0) : [0, 0, Number(v) || 0]);
      const bb = box8(child.lo, child.hi, p => [R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2], R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2], R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2]]); return { f: sdf.rotate(child.f, R), lo: bb[0], hi: bb[1] };
    }
    unsupported[op] = (unsupported[op] || 0) + 1; return null;
  }
  function block(end) { const out = []; while (st.i < tk.length && peek() !== end) { const n = statement(); if (n) out.push(n); } if (end) expect(end); return out; }
  const root = U(block(undefined)), vars = Object.fromEntries(Object.entries(env).filter(([k, v]) => k[0] !== '$' && (typeof v === 'number' || Array.isArray(v))).slice(0, 60));
  const meta = { primitives: prims, unsupported, variables: vars, resolution }, lost = Object.keys(unsupported).length ? [`modules / statements not evaluated (they contribute no geometry): ${Object.entries(unsupported).map(([k, n]) => `${k} × ${n}`).join(', ')}`] : [];
  if (polys.length) { // explicit polyhedra are imported directly (faces are clockwise seen from outside in OpenSCAD, so they are reversed)
    const V = [], Fo = []; for (const q of polys) { const base = V.length; for (const p of q.P) V.push([+p[0], +p[1], +p[2] || 0]); for (const f of q.Fc) { const idx = f.map(i => i + base).reverse(); if (idx.some(i => !(i >= base && i < V.length))) fail('OpenSCAD', 'polyhedron face references an undefined point'); fan(idx, Fo); } }
    if (root) lost.push('the script mixes polyhedron() with CSG primitives: only the explicit polyhedra were imported');
    return { kind: 'mesh', vertices: V, faces: Fo, units: 'mm', meta: { ...meta, polyhedra: polys.length }, notes: ['script parsed (not executed): explicit polyhedron() point and face lists imported as written (polyhedra inside transforms or Booleans are not transformed)', 'OpenSCAD is unitless: millimetres assumed'], lost };
  }
  if (!root || !root.lo.every(Number.isFinite) || root.hi.some((v, i) => !(v > root.lo[i]))) return { kind: null, units: 'mm', meta, notes: ['the script was parsed but produced no solid from the supported subset (cube, sphere, cylinder with transforms and Booleans)'], lost: [...lost, 'geometry not reconstructed'] };
  const padv = 0.04 * max(...root.hi.map((v, i) => v - root.lo[i])), vx = sdfToVoxels(root.f, [root.lo.map(v => v - padv), root.hi.map(v => v + padv)], resolution), mesh = marchingCubes(vx);
  return Object.assign(mesh, { units: 'mm', meta: { ...meta, grid: [vx.nx, vx.ny, vx.nz], voxel: vx.voxel }, notes: ['script parsed (not executed); CSG evaluated as a signed-distance field and iso-surfaced', `surface resolution ${N9(vx.voxel)} mm: sharp edges and features below this size are rounded or lost`, 'OpenSCAD is unitless: millimetres assumed'], lost });
}
export function toSCAD(g) {
  const v = p => `[${N9(p[0])},${N9(p[1])},${N9(p[2] || 0)}]`, h = `// Designer Materials ${g.kind}; units: ${g.units || 'unspecified'}\n`;
  if (g.kind === 'mesh') return h + `polyhedron(points=[${g.vertices.map(v).join(',')}],\n  faces=[${g.faces.map(f => `[${f[0]},${f[2]},${f[1]}]`).join(',')}], convexity=10);\n`; // OpenSCAD faces are clockwise seen from outside
  if (g.kind === 'lattice') return h + `r = ${N9(g.radius || 0.05)};\n$fn = 12;\nunion() {\n${g.edges.map((e, i) => ` hull() { translate(${v(g.nodes[e[0]])}) sphere(${g.radii ? N9(g.radii[i]) : 'r'}); translate(${v(g.nodes[e[1]])}) sphere(${g.radii ? N9(g.radii[i]) : 'r'}); }`).join('\n')}\n}\n`;
  if (g.kind === 'rods') return h + `union() {\n${g.rods.map(r => ` translate([${N9(r.x)},${N9(r.y)},0]) rotate([0,${N9(r.tilt || 0)},${N9(r.azimuth || 0)}]) cylinder(h=${N9(r.length)}, d=${N9(r.diameter)}, $fn=${g.section === 'round' ? 24 : 6});`).join('\n')}\n}\n`;
  if (g.kind === 'layout') return h + `linear_extrude(height=${N9(g.thickness || 1)}) union() {\n${g.polygons.map(q => ` polygon([${q.pts.map(p => `[${N9(p[0])},${N9(p[1])}]`).join(',')}]);`).join('\n')}\n}\n`;
  fail('OpenSCAD', `cannot write a ${g.kind} geometry`);
}

// ---------------------------------------------------------------- Gmsh .geo (parsed, never executed)
export function parseGEO(input) {
  const tk = tokenize(toText(input), 'Gmsh GEO'), env = {}, ev = makeEval(tk, env, { fmt: 'Gmsh GEO' }), { st, peek, next } = ev, P = new Map(), inv = {}, unsupported = {}, nodes = [], edges = [], phys = [], prim = [], idx = new Map();
  const node = id => { if (!idx.has(id)) { const p = P.get(id); if (!p) return -1; idx.set(id, nodes.length); nodes.push(p); } return idx.get(id); };
  const skip = () => { let d = 0; while (st.i < tk.length) { const t = next(); if (t === '{' || t === '(') d++; else if (t === '}' || t === ')') d--; else if (t === ';' && d <= 0) return; } };
  const pending = [];
  while (st.i < tk.length) {
    const s0 = st.i, t = next();
    if (t === ';') continue;
    if (/^[A-Za-z_]\w*$/.test(t) && peek() === '=') { next(); try { env[t] = ev.expr(); } catch (e) { unsupported['unevaluated assignment'] = (unsupported['unevaluated assignment'] || 0) + 1; } skip(); continue; }
    let kw = t; while (/^[A-Za-z]/.test(peek() || '') && !['(', '=', '{'].includes(peek())) kw += ' ' + next();
    if (peek() !== '(') { unsupported[kw.split(' ')[0]] = (unsupported[kw.split(' ')[0]] || 0) + 1; st.i = s0 + 1; skip(); continue; }
    next(); let id; try { id = peek() === '"' || /^"/.test(peek()) ? next() : ev.expr(); } catch (e) { id = null; } while (st.i < tk.length && peek() !== ')') next(); next();
    if (peek() !== '=') { unsupported[kw] = (unsupported[kw] || 0) + 1; skip(); continue; } next();
    let list; try { list = ev.expr(); } catch (e) { unsupported[kw] = (unsupported[kw] || 0) + 1; skip(); continue; } skip(); if (!Array.isArray(list)) list = [list]; list = list.flat();
    const K = kw.toLowerCase(); inv[kw] = (inv[kw] || 0) + 1;
    if (K === 'point') P.set(id, [list[0] || 0, list[1] || 0, list[2] || 0]);
    else if (K === 'line' || K === 'spline' || K === 'bspline' || K === 'bezier') pending.push(['poly', list]);
    else if (K === 'circle' && list.length === 3) pending.push(['arc', list]); else if (K === 'ellipse' && list.length >= 4) pending.push(['poly', [list[0], list[list.length - 1]]]);
    else if (/^physical /.test(K)) phys.push({ kind: kw.replace(/^Physical /i, ''), name: typeof id === 'string' ? id.replace(/"/g, '') : `tag ${id}`, entities: list.length });
    else if (['box', 'sphere', 'cylinder', 'cone', 'torus', 'rectangle', 'disk', 'wedge'].includes(K)) prim.push({ type: kw, params: list });
  }
  for (const [kind, l] of pending) {
    if (kind === 'poly') for (let i = 0; i + 1 < l.length; i++) { const a = node(l[i]), b = node(l[i + 1]); if (a >= 0 && b >= 0) edges.push([a, b]); }
    else { const s = P.get(l[0]), c = P.get(l[1]), e = P.get(l[2]); if (!s || !c || !e) continue; const u = [s[0] - c[0], s[1] - c[1], s[2] - c[2]], w = [e[0] - c[0], e[1] - c[1], e[2] - c[2]], r = hypot(...u), ang = Math.acos(max(-1, min(1, (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / (r * hypot(...w) || 1)))), n = max(2, Math.ceil(ang / (PI / 16))); let prev = node(l[0]);
      for (let k = 1; k <= n; k++) { let cur; if (k === n) cur = node(l[2]); else { const t = k / n, sa = Math.sin(ang) || 1, a = Math.sin((1 - t) * ang) / sa, bq = Math.sin(t * ang) / sa; nodes.push([c[0] + a * u[0] + bq * w[0], c[1] + a * u[1] + bq * w[1], c[2] + a * u[2] + bq * w[2]]); cur = nodes.length - 1; } edges.push([prev, cur]); prev = cur; } }
  }
  for (const id of P.keys()) node(id);
  const meta = { inventory: inv, physicalGroups: phys, occPrimitives: prim.slice(0, 50), variables: Object.fromEntries(Object.entries(env).filter(([, v]) => typeof v === 'number').slice(0, 60)), unsupported }, lost = [];
  const surf = Object.entries(inv).filter(([k]) => /surface|volume|loop/i.test(k)).reduce((s, [, n]) => s + n, 0); if (surf) lost.push(`${surf} curve-loop / surface / volume definitions are inventoried but not meshed: run Gmsh and import the .msh for a mesh`);
  if (prim.length) lost.push(`${prim.length} OpenCASCADE primitives listed in meta.occPrimitives are not built`); if (Object.keys(unsupported).length) lost.push(`statements not interpreted: ${Object.entries(unsupported).map(([k, n]) => `${k} × ${n}`).join(', ')}`);
  if (!nodes.length) return { kind: null, units: '', meta, notes: ['script parsed (not executed): no Point entities found'], lost: [...lost, 'geometry not reconstructed'] };
  return { kind: 'lattice', nodes, edges, radius: 0, units: '', wireframe: true, pointCloud: !edges.length, meta, notes: ['script parsed (not executed): points, lines and circle arcs shown as a wireframe; splines drawn through their control points', 'Gmsh geometry is unitless: coordinates taken as written'], lost };
}
export function toGEO(g, { lc } = {}) {
  const o = [`// Designer Materials ${g.kind}; units: ${g.units || 'unspecified'}`, 'SetFactory("Built-in");'], pt = (i, p, l) => o.push(`Point(${i}) = {${N9(p[0])}, ${N9(p[1])}, ${N9(p[2] || 0)}, ${N9(l)}};`);
  if (g.kind === 'lattice') { const L = lc ?? max(1e-9, (g.radius || 0.1) * 2); g.nodes.forEach((p, i) => pt(i + 1, p, L)); g.edges.forEach((e, i) => o.push(`Line(${i + 1}) = {${e[0] + 1}, ${e[1] + 1}};`)); o.push(`Physical Curve("struts") = {1:${g.edges.length}};`); }
  else if (g.kind === 'layout') { let np = 0, nl = 0, ns = 0; const surf = []; for (const q of g.polygons) { const ext = Math.max(...q.pts.map(p => p[0])) - Math.min(...q.pts.map(p => p[0])), L = lc ?? max(1e-9, ext / 8), p0 = np + 1, l0 = nl + 1; q.pts.forEach(p => pt(++np, p, L)); q.pts.forEach((_, i) => o.push(`Line(${++nl}) = {${p0 + i}, ${p0 + (i + 1) % q.pts.length}};`)); ns++; o.push(`Curve Loop(${ns}) = {${q.pts.map((_, i) => l0 + i).join(', ')}};`, `Plane Surface(${ns}) = {${ns}};`); surf.push(ns); } if (surf.length) o.push(`Physical Surface("pattern") = {${surf.join(', ')}};`); }
  else if (g.kind === 'mesh') { const em = new Map(); let nl = 0; const L = lc ?? 1; g.vertices.forEach((p, i) => pt(i + 1, p, L)); const line = (a, b) => { const k = a < b ? a + ',' + b : b + ',' + a; let id = em.get(k); if (!id) { id = ++nl; em.set(k, id); o.push(`Line(${id}) = {${min(a, b) + 1}, ${max(a, b) + 1}};`); } return a < b ? id : -id; }; g.faces.forEach((f, i) => { const l = [line(f[0], f[1]), line(f[1], f[2]), line(f[2], f[0])]; o.push(`Curve Loop(${i + 1}) = {${l.join(', ')}};`, `Plane Surface(${i + 1}) = {${i + 1}};`); }); o.push(`Surface Loop(1) = {1:${g.faces.length}};`, 'Volume(1) = {1};', 'Physical Volume("solid") = {1};'); }
  else fail('Gmsh GEO', `cannot write a ${g.kind} geometry`);
  return o.join('\n') + '\n';
}

// ---------------------------------------------------------------- JSON / YAML geometry and descriptor files
const TYPED = { Uint8Array, Float32Array, Float64Array, Int32Array, Uint16Array, Uint32Array, Int16Array };
export function toJSON(g, { pretty = false } = {}) {
  return JSON.stringify({ schema: 'designer-materials.geometry', version: 1, geometry: g }, (k, v) => { if (ArrayBuffer.isView(v)) { const nm = v.constructor.name; return { __typed: nm, base64: toB64(new Uint8Array(v.buffer, v.byteOffset, v.byteLength)) }; } return typeof v === 'function' ? undefined : v; }, pretty ? 1 : undefined);
}
// Minimal YAML subset: block mappings and sequences by indentation, flow sequences [a, b], scalars, # comments.
export function parseYAML(text) {
  const L = lines(toText(text)).map(l => l.replace(/\s+#.*$/, '').replace(/^#.*$/, '')).filter(l => l.trim() && l.trim() !== '---').map(l => ({ ind: l.match(/^ */)[0].length, s: l.trim() }));
  const scalar = s => { s = s.trim(); if (s === '' || s === '~' || s === 'null') return null; if (/^\[.*\]$/.test(s)) { const inner = s.slice(1, -1), out = []; let d = 0, cur = ''; for (const ch of inner) { if (ch === '[') d++; if (ch === ']') d--; if (ch === ',' && d === 0) { out.push(scalar(cur)); cur = ''; } else cur += ch; } if (cur.trim()) out.push(scalar(cur)); return out; } if (/^\{.*\}$/.test(s)) { const o = {}; for (const kv of s.slice(1, -1).split(',')) { const i = kv.indexOf(':'); if (i > 0) o[kv.slice(0, i).trim()] = scalar(kv.slice(i + 1)); } return o; } if (/^(true|false)$/i.test(s)) return /^true$/i.test(s); if (/^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(s)) return Number(s); return s.replace(/^(["'])(.*)\1$/, '$2'); };
  let i = 0;
  const node = ind => {
    if (i >= L.length) return null;
    if (L[i].s.startsWith('- ') || L[i].s === '-') { const a = []; while (i < L.length && L[i].ind === ind && (L[i].s.startsWith('- ') || L[i].s === '-')) { const rest = L[i].s.slice(1).trim(); if (!rest) { i++; a.push(i < L.length && L[i].ind > ind ? node(L[i].ind) : null); } else if (/^[^\[{"'][^:]*:(\s|$)/.test(rest)) { L[i] = { ind: ind + 2, s: rest }; a.push(node(ind + 2)); } else { i++; a.push(scalar(rest)); } } return a; }
    const o = {}; while (i < L.length && L[i].ind === ind && !L[i].s.startsWith('- ')) { const m = /^([^:]+):\s*(.*)$/.exec(L[i].s); if (!m) fail('YAML', `cannot parse line "${L[i].s}" (only simple block YAML is supported)`); i++; if (m[2] === '') o[m[1].trim()] = i < L.length && L[i].ind > ind ? node(L[i].ind) : (i < L.length && L[i].ind === ind && L[i].s.startsWith('- ') ? node(ind) : null); else o[m[1].trim()] = scalar(m[2]); } return o;
  };
  if (!L.length) fail('YAML', 'empty document'); return node(L[0].ind);
}
// Structured (JSON/YAML) geometry: this application's own export; a plain geometry object; node/edge lattices;
// vertex/face meshes; point lists; generator definitions (`{generator:'tpms'|'lattice'|'porous', …}`); and layer-stack tables.
export function fromStructured(obj, generators = {}) {
  const revive = v => { if (v && typeof v === 'object') { if (v.__typed && TYPED[v.__typed] && typeof v.base64 === 'string') { const by = b64(v.base64); return new TYPED[v.__typed](by.buffer.slice(by.byteOffset, by.byteOffset + by.byteLength)); } for (const k of Object.keys(v)) v[k] = revive(v[k]); } return v; };
  if (!obj || typeof obj !== 'object') fail('structured geometry', 'the document is not an object');
  if (obj.schema === 'designer-materials.geometry' && obj.geometry) { const g = revive(obj.geometry); return { ...g, notes: [...(g.notes || []), 'native geometry export of this application: read back without conversion'], lost: g.lost || [] }; }
  const xyz = p => Array.isArray(p) ? [+p[0], +p[1], +p[2] || 0] : [+p.x, +p.y, +p.z || 0], kinds = ['atoms', 'lattice', 'mesh', 'voxels', 'layout', 'rods'];
  if (kinds.includes(obj.kind)) { const g = revive(obj); if (g.kind === 'voxels' && Array.isArray(g.data)) g.data = Uint8Array.from(g.data); return { ...g, notes: [...(g.notes || []), `geometry object of kind "${g.kind}" read as stored`], lost: g.lost || [] }; }
  const gen = obj.generator || obj.generate;
  if (gen) { const f = generators[String(gen).toLowerCase()]; if (!f) fail('structured geometry', `unknown generator "${gen}" (available: ${Object.keys(generators).join(', ')})`); const g = f(obj.params || obj); return { ...g, meta: { generator: gen, params: obj.params || obj }, notes: [...(g.notes || []), `geometry generated from the parametric definition (generator "${gen}")`], lost: [] }; }
  if (Array.isArray(obj.nodes) && (Array.isArray(obj.edges) || Array.isArray(obj.struts) || Array.isArray(obj.beams) || Array.isArray(obj.links))) {
    const E = obj.edges || obj.struts || obj.beams || obj.links, ids = obj.nodes.map((n, i) => (n && !Array.isArray(n) && n.id !== undefined) ? n.id : i), g = graphOut(ids, obj.nodes.map(xyz), E.map(e => Array.isArray(e) ? [e[0], e[1]] : [e.source ?? e.from ?? e.i ?? e.n1, e.target ?? e.to ?? e.j ?? e.n2]), E.map(e => Array.isArray(e) ? NaN : +(e.radius ?? e.r ?? (e.diameter !== undefined ? e.diameter / 2 : NaN))), { name: obj.name, cell: obj.cell || obj.period || obj.periodicity }, 'structured lattice');
    if (obj.units) g.units = obj.units; if (Number.isFinite(+obj.radius)) g.radius = +obj.radius; else if (Number.isFinite(+obj.strut)) g.radius = +obj.strut / 2; if (obj.period || obj.cell) g.period = obj.period || obj.cell; if (obj.units) g.notes = g.notes.filter(n => !/no length unit/.test(n)); return g; }
  if (Array.isArray(obj.vertices) && (Array.isArray(obj.faces) || Array.isArray(obj.triangles))) { const F = []; for (const f of obj.faces || obj.triangles) fan(Array.from(f), F); const V = obj.vertices.map(xyz); if (F.some(f => f.some(i => !(i >= 0 && i < V.length)))) fail('structured geometry', 'a face references an undefined vertex'); return { kind: 'mesh', vertices: V, faces: F, units: obj.units || '', notes: obj.units ? [] : ['no "units" field: lengths taken as stored'], lost: [] }; }
  if (Array.isArray(obj.polygons)) return { kind: 'layout', polygons: obj.polygons.map(q => ({ layer: q.layer ?? 1, pts: (q.pts || q.points || q).map(p => [+p[0], +p[1]]) })), units: obj.units || 'nm', notes: obj.units ? [] : ['no "units" field: nanometres assumed for a layout'], lost: [] };
  const stack = obj.layers || obj.stack || obj.layerStack;
  if (Array.isArray(stack) && stack.length && typeof stack[0] === 'object') { let z = 0; const layers = stack.map((l, i) => { const t = +(l.thickness ?? l.thickness_nm ?? l.t ?? l.d ?? 0), o = { index: i + 1, material: l.material ?? l.name ?? l.mat ?? `layer ${i + 1}`, thickness: t, z0: z, z1: z + t, ...(l.n !== undefined ? { n: l.n } : {}), ...(l.layer !== undefined ? { gdsLayer: l.layer } : {}) }; z += t; return o; }); return { kind: null, units: obj.units || 'nm', layerStack: layers, meta: { layers: layers.length, totalThickness: z }, notes: [`layer-stack table with ${layers.length} layers, total thickness ${N9(z)} ${obj.units || 'nm (assumed)'}: this is vertical process information to combine with a 2-D layout, not a geometry by itself`], lost: [] }; }
  if (Array.isArray(obj.points) || Array.isArray(obj.coordinates)) { const P = (obj.points || obj.coordinates).map(xyz); return { kind: 'lattice', nodes: P, edges: [], radius: 0, units: obj.units || '', pointCloud: true, notes: ['point list imported as a point cloud'], lost: [] }; }
  fail('structured geometry', `no recognised geometry schema (top-level keys: ${Object.keys(obj).slice(0, 12).join(', ')}). Expected one of: kind, generator, nodes+edges, vertices+faces, polygons, layers, points`);
}
export { h5Open };
