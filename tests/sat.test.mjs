// Self-test of the ACIS / ShapeManager path: the SAT text and SAB binary readers (old and new record layouts, bare
// and "$" references, 4- and 8-byte binary variants), the ACIS → STEP translator and the importer wiring (.sat / .sab
// files, an Inventor compound file, a Fusion archive, a JT file carrying an XT B-rep segment, SOLIDWORKS
// multi-configuration parts and assemblies). Default export: async () => {pass, fail, log}. Every fixture is
// synthesised here, so the test runs offline; the tessellation checks are skipped (with a logged note) when the
// OpenCASCADE kernel pack is absent.
import * as F from '../js/geometry/formats.js';
import * as K from '../js/geometry/kernels.js';
import * as N from '../js/geometry/native.js';
import * as S from '../js/geometry/sat.js';
import * as IM from '../js/geometry/importer.js';
import { fixtures as X } from './xt.test.mjs';

const { enc, cat, u32, zlibStored, swFile, xtBinary, solid: xtSolid, block: xtBlock } = X, PI = Math.PI;
// abstract record tokens: reference, number, integer, logical (with its SAT word pair), vector, unbounded interval end, quoted string, bare word
const p = n => ['p', n], f = x => ['f', x], d = n => ['d', n], L = (v, a, b) => ['l', v, a, b], v = a => ['v', a], I = ['i'], str = s => ['s', s], w = s => ['w', s];
const fwd = r => L(!!r, 'forward', 'reversed');
// a model: {faces: [{surface: [class, ...tokens], rev, loops: [[{edge, rev}]]}], edges: [{curve: [class, ...tokens] | null, a, b, t0, t1, rev}], points}
// → records [[class name, ...tokens]] in the layout of ACIS 7.0 and later (style 'new') or of ACIS 4.0 (style 'old')
function records({ faces, edges, points, transform }, style = 'new') {
  const R = [], H = style === 'new' ? [p(-1), d(-1), p(-1)] : [p(-1)]; let n = 0; const id = () => n++, BODY = id(), LUMP = id(), SHELL = id(), TF = transform ? id() : -1;
  const fa = faces.map(() => id()), su = faces.map(() => id()), lp = faces.map(fc => fc.loops.map(() => id())), co = faces.map(fc => fc.loops.map(l => l.map(() => id()))), ed = edges.map(() => id()), cu = edges.map(e => (e.curve ? id() : -1)), vx = points.map(() => id()), pt = points.map(() => id());
  const uses = edges.map(() => []); faces.forEach((fc, i) => fc.loops.forEach((l, j) => l.forEach((h, k) => uses[h.edge].push(co[i][j][k]))));
  const vEdge = points.map((_, i) => ed[edges.findIndex(e => e.a === i || e.b === i)]);
  R[BODY] = ['body', ...H, p(LUMP), p(-1), p(TF)]; R[LUMP] = ['lump', ...H, p(-1), p(SHELL), p(BODY)]; R[SHELL] = ['shell', ...H, p(-1), p(-1), p(fa[0]), p(-1), p(LUMP)];
  if (transform) R[TF] = ['transform', ...(style === 'new' ? [p(-1), d(-1)] : [p(-1)]), ...transform.R.map(f), ...transform.t.map(f), f(transform.s || 1), L(true, 'no_rotate', 'rotate'), L(false, 'no_reflect', 'reflect'), L(false, 'no_shear', 'shear')];
  faces.forEach((fc, i) => {
    R[fa[i]] = ['face', ...H, p(fa[i + 1] ?? -1), p(lp[i][0] ?? -1), p(SHELL), p(-1), p(su[i]), fwd(fc.rev), L(false, 'single', 'double')]; R[su[i]] = [fc.surface[0], ...H, ...fc.surface.slice(1)];
    fc.loops.forEach((l, j) => { R[lp[i][j]] = ['loop', ...H, p(lp[i][j + 1] ?? -1), p(co[i][j][0]), p(fa[i])]; l.forEach((h, k) => { const me = co[i][j][k], us = uses[h.edge], other = us.length > 1 ? us[(us.indexOf(me) + 1) % us.length] : -1; R[me] = ['coedge', ...H, p(co[i][j][(k + 1) % l.length]), p(co[i][j][(k + l.length - 1) % l.length]), p(other), p(ed[h.edge]), fwd(h.rev), p(lp[i][j]), ...(style === 'new' ? [d(0)] : []), p(-1)]; }); });
  });
  edges.forEach((e, i) => { R[ed[i]] = ['edge', ...H, p(vx[e.a]), ...(style === 'new' ? [f(e.t0)] : []), p(vx[e.b]), ...(style === 'new' ? [f(e.t1)] : []), p(uses[i][0]), p(cu[i]), fwd(e.rev), ...(style === 'new' ? [str('unknown')] : [])]; if (e.curve) R[cu[i]] = [e.curve[0], ...H, ...e.curve.slice(1)]; });
  points.forEach((q, i) => { R[vx[i]] = ['vertex', ...H, p(vEdge[i]), ...(style === 'new' ? [d(0)] : []), p(pt[i])]; R[pt[i]] = ['point', ...H, v(q)]; });
  return R;
}
// SAT text: version line, (2.0 and later) product / kernel / date strings and the units line, one record per line
function satText(R, { version = 700, scale = 1, bare = false, explicit = false } = {}) {
  const tok = t => t[0] === 'p' ? (bare ? String(t[1]) : '$' + t[1]) : t[0] === 'f' || t[0] === 'd' ? String(t[1]) : t[0] === 'l' ? (version < 200 && t[2] === 'forward' ? (t[1] ? '1' : '0') : t[1] ? t[3] : t[2]) : t[0] === 'v' ? t[1].join(' ') : t[0] === 'i' ? 'I' : t[0] === 's' ? `@${t[1].length} ${t[1]}` : t[1];
  let s = `${version} 0 1 0\n`; if (version >= 200) s += `7 fixture 8 ACIS 7.0 24 Mon Jan 01 00:00:00 2024\n${scale} 9.9999999999999995e-007 1e-010\n`;
  R.forEach((r, i) => { s += `${explicit ? '-' + i + ' ' : ''}${r[0]} ${r.slice(1).map(tok).join(' ')} #\n`; });
  return enc(s + (version >= 200 ? 'End-of-ACIS-data\n' : ''));
}
// SAB binary: signature, four counts, tagged header data, tagged records. wide: the 8-byte "ASM BinaryFile8" variant.
function sabBinary(R, { sig = 'ACIS BinaryFile', version = 700, scale = 1, wide = false, history = false, product = 'fixture' } = {}) {
  const out = [], W = wide ? 8 : 4, int = (n, k = 4) => { const b = new Uint8Array(k), q = new DataView(b.buffer); if (k === 8) q.setBigInt64(0, BigInt(n), true); else q.setInt32(0, n, true); out.push(...b); }, dbl = x => { const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, x, true); out.push(...b); }, s8 = (tag, s) => { out.push(tag, s.length, ...enc(s)); };
  const ident = name => { const parts = name.split('-'); parts.forEach((q, i) => s8(i === parts.length - 1 ? 13 : 14, q)); };
  out.push(...enc(sig)); for (const n of [version, 0, 1, 0]) int(n, W); s8(7, product); s8(7, 'ASM 223.0'); s8(7, 'Mon Jan 01 00:00:00 2024'); out.push(6); dbl(scale); out.push(6); dbl(1e-6); out.push(6); dbl(1e-10);
  for (const r of R) { ident(r[0]); for (const t of r.slice(1)) { if (t[0] === 'p') { out.push(12); int(t[1], W); } else if (t[0] === 'd') { out.push(4); int(t[1], W); } else if (t[0] === 'f') { out.push(6); dbl(t[1]); } else if (t[0] === 'l') out.push(t[1] ? 10 : 11); else if (t[0] === 'v') { out.push(19); t[1].forEach(dbl); } else if (t[0] === 'i') out.push(11); else if (t[0] === 's') s8(7, t[1]); else if (t[1] === '{') out.push(15); else if (t[1] === '}') out.push(16); else if (t[1] in { open: 1, closed: 1, periodic: 1, none: 1, full: 1 }) { out.push(21); int({ open: 0, closed: 1, periodic: 2, none: 0, full: 0 }[t[1]], W); } else s8(13, t[1]); } out.push(17); }
  const kernel = sig.startsWith('ASM') ? 'ASM' : 'ACIS'; if (history) { ident(`Begin-of-${kernel}-History-Data`); ident('history_stream'); out.push(4); int(1, W); out.push(17); } else ident(`End-of-${kernel}-data`);
  return Uint8Array.from(out);
}
const plane = (o, n, u) => ['plane-surface', v(o), v(n), v(u), L(false, 'forward_v', 'reverse_v'), I, I, I, I];
const line = (a, b) => { const q = b.map((x, i) => x - a[i]), l = Math.hypot(...q); return { curve: ['straight-curve', v(a), v(q.map(x => x / l)), I, I], t0: 0, t1: l }; };
// block [0, a] × [0, b] × [0, c]: six planes with outward normals, twelve straight edges, eight vertices
function block(a, b, c, opts = {}) {
  const size = [a, b, c], points = [], edges = [], key = new Map(), faces = [], unitv = k => [0, 1, 2].map(i => (i === k ? 1 : 0));
  for (let i = 0; i < 8; i++) points.push([i & 1 ? a : 0, i & 2 ? b : 0, i & 4 ? c : 0]);
  const edge = (m, n) => { const lo = Math.min(m, n), hi = Math.max(m, n), k = lo * 8 + hi; if (!key.has(k)) { key.set(k, edges.length); edges.push({ ...line(points[lo], points[hi]), a: lo, b: hi, rev: false }); } return { edge: key.get(k), rev: m > n }; };
  for (let k = 0; k < 3; k++) for (const s of [0, 1]) { const uu = (k + 1) % 3, vv = (k + 2) % 3, at = (i, j) => (s << k) | (i << uu) | (j << vv), ring = s ? [at(0, 0), at(1, 0), at(1, 1), at(0, 1)] : [at(0, 0), at(0, 1), at(1, 1), at(1, 0)], o = [0, 0, 0]; o[k] = s ? size[k] : 0;
    faces.push({ surface: plane(o, unitv(k).map(x => (s ? x : -x)), unitv(uu)), rev: false, loops: [ring.map((m, i) => edge(m, ring[(i + 1) % 4]))] }); }
  if (opts.spline) {
    // the top face on an exact bilinear B-spline surface (control points with u fastest), a second face through "{ ref 0 }" is not needed here;
    // one vertical edge on an exact B-spline curve, saved reversed
    const top = faces[5]; top.surface = ['spline-surface', fwd(false), w('{'), w('exactsur'), w('full'), w('nubs'), d(1), d(1), w('open'), w('open'), w('none'), w('none'), d(2), d(2), f(0), d(1), f(a), d(1), f(0), d(1), f(b), d(1), v([0, 0, c]), v([a, 0, c]), v([0, b, c]), v([a, b, c]), f(0), w('}'), I, I, I, I];
    const e = edges[key.get(0 * 8 + 4)]; e.curve = ['intcurve-curve', fwd(true), w('{'), w('exactcur'), w('full'), w('nubs'), d(1), w('open'), d(2), f(-c), d(1), f(0), d(1), v([0, 0, c]), v([0, 0, 0]), f(0), w('}'), I, I];
  }
  return { faces, edges, points };
}
// cylinder of radius r and height h about z: one cylindrical and two planar faces, two closed circular edges with one vertex each
function cylinder(r, h, cone = 0) {
  const r1 = r + h * cone, circle = (z, rad) => ['ellipse-curve', v([0, 0, z]), v([0, 0, 1]), v([rad, 0, 0]), f(1), I, I], ang = Math.atan(Math.abs(cone));
  // a cone's sine and cosine have opposite signs when its radius shrinks along the axis
  return { points: [[r, 0, 0], [r1, 0, h]], edges: [{ curve: circle(0, r), a: 0, b: 0, t0: 0, t1: 2 * PI, rev: false }, { curve: circle(h, r1), a: 1, b: 1, t0: 0, t1: 2 * PI, rev: false }],
    faces: [{ surface: ['cone-surface', v([0, 0, 0]), v([0, 0, 1]), v([r, 0, 0]), f(1), I, I, f(cone ? (cone > 0 ? 1 : -1) * Math.sin(ang) : 0), f(cone ? Math.cos(ang) : 1), f(r), L(false, 'forward', 'reversed'), I, I, I, I], rev: false, loops: [[{ edge: 0, rev: false }], [{ edge: 1, rev: true }]] },
      { surface: plane([0, 0, 0], [0, 0, -1], [1, 0, 0]), rev: false, loops: [[{ edge: 0, rev: true }]] }, { surface: plane([0, 0, h], [0, 0, 1], [1, 0, 0]), rev: false, loops: [[{ edge: 1, rev: false }]] }] };
}
// closed surfaces without any loop: a sphere of radius r (negative: hollow sense) and a torus
const sphere = r => ({ points: [], edges: [], faces: [{ surface: ['sphere-surface', v([0, 0, 0]), f(r), v([1, 0, 0]), v([0, 0, 1]), L(false, 'forward_v', 'reverse_v'), I, I, I, I], rev: r < 0, loops: [] }] });
const torus = (R, r) => ({ points: [], edges: [], faces: [{ surface: ['torus-surface', v([0, 0, 0]), v([0, 0, 1]), f(R), f(r), v([1, 0, 0]), L(false, 'forward_v', 'reverse_v'), I, I, I, I], rev: false, loops: [] }] });
// OLE2 compound file (512-byte sectors, one FAT sector) with one storage holding the given streams
function cfbStorage(storage, streams) {
  const SS = 512, small = streams.filter(s => s.data.length < 4096), big = streams.filter(s => s.data.length >= 4096), sec = n => Math.ceil(n / SS);
  const mini = new Uint8Array(Math.max(SS, sec(small.reduce((n, s) => n + Math.ceil(s.data.length / 64) * 64, 0)) * SS)), miniFat = new Uint32Array(SS / 4 * 4).fill(0xFFFFFFFF), fat = new Uint32Array(SS / 4).fill(0xFFFFFFFF); let mp = 0;
  for (const s of small) { s.start = mp; const n = Math.ceil(s.data.length / 64); mini.set(s.data, mp * 64); for (let i = 0; i < n; i++) miniFat[mp + i] = i === n - 1 ? 0xFFFFFFFE : mp + i + 1; mp += n; }
  const nMF = Math.ceil(mp / (SS / 4)) || 1, nDir = Math.ceil((streams.length + 2) / 4); let next = 1; const chain = n => { const s0 = next; for (let i = 0; i < n; i++) fat[next + i] = i === n - 1 ? 0xFFFFFFFE : next + i + 1; next += n; return s0; };
  fat[0] = 0xFFFFFFFD; const dirStart = chain(nDir), mfStart = chain(nMF), miniStart = chain(mini.length / SS), blobs = []; for (const s of big) { s.start = chain(sec(s.data.length)); const q = new Uint8Array(sec(s.data.length) * SS); q.set(s.data); blobs.push(q); }
  if (next > SS / 4) throw new Error('fixture too large for one FAT sector');
  const dir = new Uint8Array(nDir * SS), dv = new DataView(dir.buffer), ent = (i, name, type, start, size, child, right) => { const o = 128 * i; for (let k = 0; k < name.length; k++) dv.setUint16(o + 2 * k, name.charCodeAt(k), true); dv.setUint16(o + 64, 2 * name.length + 2, true); dir[o + 66] = type; dir[o + 67] = 1; dv.setUint32(o + 68, 0xFFFFFFFF, true); dv.setUint32(o + 72, right, true); dv.setUint32(o + 76, child, true); dv.setUint32(o + 116, start, true); dv.setUint32(o + 120, size, true); };
  ent(0, 'Root Entry', 5, miniStart, mini.length, 1, 0xFFFFFFFF); ent(1, storage, 1, 0, 0, 2, 0xFFFFFFFF); streams.forEach((s, i) => ent(i + 2, s.name, 2, s.start, s.data.length, 0xFFFFFFFF, i + 1 < streams.length ? i + 3 : 0xFFFFFFFF));
  const head = new Uint8Array(SS), hv = new DataView(head.buffer); head.set([0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1]); hv.setUint16(24, 0x3E, true); hv.setUint16(26, 3, true); hv.setUint16(28, 0xFFFE, true); hv.setUint16(30, 9, true); hv.setUint16(32, 6, true); hv.setUint32(44, 1, true); hv.setUint32(48, dirStart, true); hv.setUint32(56, 4096, true); hv.setUint32(60, mfStart, true); hv.setUint32(64, nMF, true); hv.setUint32(68, 0xFFFFFFFE, true); for (let i = 0; i < 109; i++) hv.setUint32(76 + 4 * i, i ? 0xFFFFFFFF : 0, true);
  return cat(head, new Uint8Array(fat.buffer), dir, new Uint8Array(miniFat.buffer, 0, nMF * SS), mini, ...blobs);
}
const utf16le = s => { const o = new Uint8Array(2 * s.length); for (let i = 0; i < s.length; i++) { o[2 * i] = s.charCodeAt(i) & 255; o[2 * i + 1] = s.charCodeAt(i) >> 8; } return o; };
// JT 9.5 file with one XT B-rep segment: 80-byte version line, byte order, table of contents, zlib-wrapped segment data
function jtFile(xt) {
  const le32 = n => u32(n), guid = k => new Uint8Array(16).fill(k), data = cat(u32(xt.length + 40), guid(7), Uint8Array.from([1]), u32(1, 0, 0, 0, xt.length), xt), z = zlibStored(data), seg = cat(guid(2), le32(17), le32(24 + 9 + z.length), le32(2), le32(z.length + 1), Uint8Array.from([2]), z);
  const head = cat(enc('Version 9.5 JT  DM 8.0.0.0'.padEnd(75, ' ') + ' \n\r\n '), Uint8Array.from([0]), le32(0), le32(0), guid(1)), tocOff = head.length; new DataView(head.buffer).setInt32(85, tocOff, true);
  const segOff = tocOff + 4 + 28; return cat(head, le32(1), guid(2), le32(segOff), le32(seg.length), le32(17 << 24), seg);
}

export default async function run() {
  const log = []; let pass = 0, fail = 0; const t00 = Date.now();
  const ok = (cond, msg) => { if (cond) pass++; else { fail++; log.push('✗ ' + msg); } };
  const near = (a, b, rel) => Math.abs(a - b) <= rel * Math.abs(b), count = (s, re) => (s.match(re) || []).length;
  const have = await K.kernelAvailable('occt'), A = 40, B = 30, C = 20, Rr = 10, H = 50;
  const blockNew = satText(records(block(A, B, C))), blockOld = satText(records(block(A, B, C), 'old'), { version: 400 }), blockBin = sabBinary(records(block(A, B, C))), cylNew = satText(records(cylinder(Rr, H)));
  const bbox = g => [0, 1, 2].map(a => Math.max(...g.vertices.map(q => q[a])) - Math.min(...g.vertices.map(q => q[a])));

  // ---- reader: text (7.0 layout, 4.0 layout, 1.x bare references, explicit indices) and binary give the same records
  {
    const mn = S.readSAT(blockNew), mo = S.readSAT(blockOld), mb = S.readSAT(blockBin);
    ok(mn.format === 'sat' && mn.version === 700 && mn.product === 'fixture' && mn.acis === 'ACIS 7.0' && mn.scale === 1 && mn.resabs === 1e-6 && mn.terminated && mb.format === 'sab' && mb.signature === 'ACIS BinaryFile' && mb.version === 700 && mb.scale === 1 && mb.terminated, 'readSAT: header of a text and a binary save file (version, product, kernel string, units, resabs, end marker)');
    ok(['body', 'lump', 'shell', 'face', 'loop', 'coedge', 'edge', 'vertex', 'point', 'straight-curve', 'plane-surface'].every(k => mn.counts[k] === mb.counts[k] && mn.counts[k] === mo.counts[k]) && mn.counts.face === 6 && mn.counts.coedge === 24 && mn.counts.edge === 12 && mn.counts.vertex === 8 && mn.entities.size === 85, 'readSAT: text (7.0 and 4.0 layouts) and binary decode to the same record inventory (6 faces, 24 coedges, 12 edges, 8 vertices)');
    const e = [...mn.entities.values()].find(r => r.name === 'edge'), eb = mb.entities.get(e.i); ok(e.t[0] instanceof S.Ptr && e.t[0].i === -1 && e.t.some(x => x instanceof S.Str && x.s === 'unknown') && e.t.includes('forward') && eb.t.includes(false) && eb.t.filter(x => x instanceof S.Ptr).length === e.t.filter(x => x instanceof S.Ptr).length, 'readSAT: references, quoted strings and sense words in text; references and logical tags in binary');
    const pl = [...mb.entities.values()].find(r => r.name === 'plane-surface'); ok(pl.base === 'surface' && pl.t.filter(x => typeof x === 'number').length === 10, 'readSAT: binary position / vector tags are flattened to numbers and dashed class names are rebuilt');
    const mx = S.readSAT(satText(records(block(A, B, C)), { explicit: true })); ok(mx.entities.size === 85 && mx.entities.get(0).name === 'body', 'readSAT: explicit "-n" record indices');
    const m8 = S.readSAT(sabBinary(records(block(A, B, C)), { sig: 'ASM BinaryFile8', version: 22300, wide: true, scale: 10 })), m4 = S.readSAT(sabBinary(records(block(A, B, C)), { sig: 'ASM BinaryFile4', version: 22300, scale: 10, history: true }));
    ok(m8.signature === 'ASM BinaryFile8' && m8.version === 22300 && m8.entities.size === 85 && m8.scale === 10 && m4.entities.size === 85 && m4.history && m4.terminated, 'readSAT: ShapeManager "ASM BinaryFile4" and 8-byte "ASM BinaryFile8" variants; a roll-back history section ends the model');
    const junk = new Uint8Array(300).map((_, i) => (i * 29) & 255), hits = S.findACIS(cat(junk, blockBin, junk, sabBinary(records(cylinder(Rr, H)), { sig: 'ASM BinaryFile4', version: 22300 }))); ok(hits.length === 2 && hits[0].offset === 300 && hits[1].signature === 'ASM BinaryFile4' && S.readSAT(cat(junk, blockBin), { offset: 300 }).entities.size === 85, 'findACIS: binary streams are located inside another buffer by signature and read from their offset');
    let msg = ''; try { S.readSAT(enc('hello world')); } catch (err) { msg = err.message; } ok(/SAT/.test(msg), 'readSAT: data without a save-file header is rejected');
    msg = ''; try { S.readSAT(blockBin.subarray(0, blockBin.length - 200)); } catch (err) { msg = err.message; } let cut = ''; try { const tr = S.readSAT(blockNew.subarray(0, blockNew.length - 400)); cut = tr.warnings.join(' '); } catch (err) { cut = err.message; } ok(/end of the binary data/.test(msg) && /end-of-data marker|SAT: malformed/.test(cut), 'readSAT: truncation is reported (binary: explicit error; text: error or a missing-end-marker warning)');
  }

  // ---- translator
  {
    const tb = S.satToStep(S.readSAT(blockNew), { name: 'block.sat' }), to = S.satToStep(S.readSAT(blockOld)), tx = S.satToStep(S.readSAT(blockBin)), tc = S.satToStep(S.readSAT(cylNew), { name: 'cyl.sat' });
    ok(tb.stats.faces === 6 && tb.stats.translated === 6 && tb.stats.exactEdges === 12 && tb.stats.approximatedEdges === 0 && tb.bodies.length === 1 && tb.bodies[0].closed && tb.bodies[0].type === 'solid' && !tb.stats.offSurface && !tb.stats.paramMismatch, 'satToStep(block): 6 of 6 faces and 12 exact edges in one closed solid; boundary and parameter self-checks clean');
    ok(/^ISO-10303-21;/.test(tb.step) && /AUTOMOTIVE_DESIGN/.test(tb.step) && count(tb.step, /=MANIFOLD_SOLID_BREP\(/g) === 1 && count(tb.step, /=ADVANCED_FACE\(/g) === 6 && count(tb.step, /=PLANE\(/g) === 6 && count(tb.step, /=EDGE_CURVE\(/g) === 12 && count(tb.step, /=ORIENTED_EDGE\(/g) === 24 && count(tb.step, /=VERTEX_POINT\(/g) === 8 && /SI_UNIT\(\.MILLI\.,\.METRE\.\)/.test(tb.step) && /PRODUCT\('block'/.test(tb.step) && /END-ISO-10303-21;\s*$/.test(tb.step), 'satToStep(block): well-formed AP214 text in millimetres with the expected entity counts');
    ok(to.stats.translated === 6 && to.bodies[0].closed && tx.stats.translated === 6 && tx.bodies[0].closed && count(to.step, /=EDGE_CURVE\(/g) === 12, 'satToStep: the 4.0 text layout and the binary encoding translate to the same solid');
    ok(tc.stats.translated === 3 && count(tc.step, /=CYLINDRICAL_SURFACE\(/g) === 1 && count(tc.step, /=CIRCLE\(/g) === 2 && count(tc.step, /=FACE_BOUND\(/g) === 2 && tc.bodies[0].closed, 'satToStep(cylinder): cone-surface with zero sine becomes a cylindrical surface with two circular bounds');
    const sc = S.satToStep(S.readSAT(satText(records(block(1, 2, 3)), { scale: 25.4 }))); ok(/CARTESIAN_POINT\('',\(25\.4,50\.8,76\.(2|19999\d*)\)\)/.test(sc.step), 'satToStep: lengths are multiplied by the millimetres-per-unit of the header (inch file)');
    const un = S.satToStep(S.readSAT(satText(records(block(A, B, C), 'old'), { version: 105 }))); ok(un.stats.translated === 6 && un.notes.some(s => /predates the units record/.test(s)), 'satToStep: a 1.x file (numeric senses, no units record) translates, with the unit assumption stated');
    const bare = S.satToStep(S.readSAT(satText(records(block(A, B, C), 'old'), { version: 102, bare: true }))); ok(bare.stats.translated === 6 && bare.bodies[0].closed, 'satToStep: bare integer references of the earliest save files');
    const sp = S.satToStep(S.readSAT(satText(records(block(A, B, C, { spline: true })))));
    ok(sp.stats.translated === 6 && sp.stats.surfaces['spline (exactsur)'] === 1 && sp.stats.edgeCurves['intcurve (exactcur)'] === 1 && sp.stats.approximatedFaces === 0 && sp.stats.approximatedEdges === 0 && count(sp.step, /=B_SPLINE_SURFACE_WITH_KNOTS\(/g) === 1 && count(sp.step, /=B_SPLINE_CURVE_WITH_KNOTS\(/g) === 1 && /B_SPLINE_SURFACE_WITH_KNOTS\('',1,1,[^;]*\(2,2\),\(2,2\)/.test(sp.step), 'satToStep: exact spline surface and intcurve data are written as B-splines with their clamping knots restored');
    // a procedural curve with only a summary of its approximation costs the faces it bounds
    const bad = block(A, B, C); bad.edges[0].curve = ['intcurve-curve', fwd(false), w('{'), w('bldcur'), w('summary'), d(2), f(0), f(1), w('}'), I, I]; const tl = S.satToStep(S.readSAT(satText(records(bad))));
    ok(tl.stats.translated === 4 && Object.keys(tl.stats.skipped).some(k => /bldcur/.test(k)) && !tl.bodies[0].closed && /=OPEN_SHELL\(/.test(tl.step) && !/MANIFOLD_SOLID_BREP/.test(tl.step), 'satToStep: an edge curve without stored control points costs exactly the two faces it bounds, reported by subtype, and the rest is an open shell');
    const mv = block(A, B, C); mv.transform = { R: [0, 1, 0, -1, 0, 0, 0, 0, 1], t: [100, 0, 0] }; const tm = S.satToStep(S.readSAT(satText(records(mv)))); ok(/CARTESIAN_POINT\('',\(70\.,40\.,20\.\)\)/.test(tm.step) && tm.stats.transformed === 1, 'satToStep: a body transform acts on row vectors (p · R + t)');
    ok(S.satToStep(S.readSAT(enc('700 0 0 0\n7 fixture 8 ACIS 7.0 24 Mon Jan 01 00:00:00 2024\n1 1e-06 1e-10\nEnd-of-ACIS-data\n'))).step === null, 'satToStep: a save file without bodies yields no STEP text');
  }

  // ---- containers
  const asm = sabBinary(records(block(A / 10, B / 10, C / 10)), { sig: 'ASM BinaryFile4', version: 22600, scale: 10, history: true, product: 'Inventor' });
  const ipt = cfbStorage('RSeStorage', [{ name: 'Mabc', data: cat(u32(25), enc('RSe Meta Stream Version 8'), new Uint8Array(15), u32(13), utf16le('PmBRepSegment'), new Uint8Array(40)) }, { name: 'Babc', data: cat(new Uint8Array(18).fill(1), zlibStored(cat(new Uint8Array(77).fill(3), asm))) }, { name: 'Mdef', data: cat(u32(25), enc('RSe Meta Stream Version 8'), new Uint8Array(15), u32(11), utf16le('PmDCSegment'), new Uint8Array(40)) }, { name: 'Bdef', data: cat(new Uint8Array(18).fill(1), zlibStored(sabBinary(records(cylinder(1, 5)), { sig: 'ASM BinaryFile4', version: 22600, scale: 10 }))) }]);
  const smbh = sabBinary(records(block(A / 10, B / 10, C / 10)), { sig: 'ASM BinaryFile8', version: 22300, wide: true, scale: 30, history: true, product: 'Autodesk Neutron' }), smb = sabBinary(records(cylinder(1, 5)), { sig: 'ASM BinaryFile8', version: 22300, wide: true, scale: 10, product: 'Autodesk Neutron' });
  const f3d = F.zipStore([{ name: 'FusionAssetName[Active]/Breps.BlobParts/BREP.1.smb', data: smb }, { name: 'FusionAssetName[Active]/Breps.BlobParts/BREP.2.smbh', data: smbh }, { name: 'Manifest.dat', data: new Uint8Array(16) }]);
  const xtPart = (a, b, c) => xtBinary(xtSolid(xtBlock(a, b, c), true), true, false), partition = xt => cat(u32(1), new Uint8Array(16), u32(xt.length, 0), zlibStored(xt));
  {
    const r = N.parseOleCad(ipt, 'Inventor'), e = r.extra.embedded || []; ok(e.length === 1 && e[0].format === 'acis' && e[0].segment === 'PmBRepSegment' && S.readSAT(e[0].bytes).counts.face === 6 && r.recovered.some(s => /1 further stream/.test(s)) && r.meta.segmentNames.length === 2, 'parseOleCad(Inventor): segments are paired by key and inflated; the B-rep segment\'s ShapeManager stream is offered and the sketch-data stream left aside');
    const rf = N.parseFusion(f3d), ef = rf.extra.embedded || []; ok(ef.length === 1 && /smbh$/.test(ef[0].name) && ef[0].unitScale === 10 && rf.recovered.some(s => /1 further \.smb blob/.test(s)), 'parseFusion: the history blob (.smbh) holds the current bodies; the .smb beside it is left aside');
    const rj = N.parseJT(jtFile(xtPart(0.04, 0.03, 0.02))), ej = rj.extra.embedded || []; ok(ej.length === 1 && ej[0].format === 'parasolid' && ej[0].schema === 'SCH_1200000_12006' && rj.meta.segmentCounts['XT B-rep (Parasolid)'] === 1, 'parseJT: the Parasolid transmit stream inside a zlib-wrapped XT B-rep segment is unpacked');
    const tree = N.swInstanceTree('<swSolidWorks><swHeader><swFile id="3" swDocType="ASSEMBLY" swPath="C:\\a\\top.SLDASM"/><swFile id="6" swDocType="PART" swPath="C:\\a\\blk.SLDPRT"/><swFile id="9" swDocType="ASSEMBLY" swPath="C:\\a\\sub.SLDASM"/></swHeader><swModelList><swModel id="5" swName="blk" swFileRef="6"/><swModel id="8" swName="sub" swFileRef="9"><swReference swName="blk" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="5" swTransform="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0.05 1"/></swModel><swModel id="2" swName="top" swFileRef="3" swBoundingBox="0 0 0 0.1 0.04 0.07"><swReference swName="blk" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="5" swTransform="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/><swReference swName="sub" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="8" swTransform="0 1 0 0 -1 0 0 0 0 0 1 0 0.1 0 0 1"/></swModel></swModelList><swConfigurationList><swConfiguration swName="Default" swMostRecentConfiguration="YES" swModelRef="2"/></swConfigurationList></swSolidWorks>');
    ok(tree && tree.name === 'top' && tree.components.length === 2 && tree.components[1].path.join('/') === 'sub<1>/blk<1>' && tree.components[1].depth === 1 && tree.components[1].matrix.map(x => +x.toFixed(6)).join() === '0,-1,0,0.1,1,0,0,0,0,0,1,0.05' && tree.bbox[3] === 0.1 && /blk\.SLDPRT$/.test(tree.components[1].file), 'swInstanceTree: components, files and row-vector placements, with a sub-assembly expanded and its placement composed');
  }

  // ---- importer: registry, detection, fallback without the kernel, never throwing
  {
    const by = Object.fromEntries(IM.FORMATS.map(x => [x.id, x]));
    ok(by.sat.support === 'partial' && by.sab.support === 'partial' && by.inventor.support === 'partial' && by.fusion.support === 'partial' && /translated/.test(by.sat.note) && by.sab.ext.includes('.smb') && by.catia.support === 'inventory' && by.creo.support === 'inventory', 'registry: SAT, SAB, Inventor and Fusion are "partial" with what is and is not translated in the note; CATIA and Creo stay at inventory');
    ok(IM.detectFormat('a.sat', blockNew).id === 'sat' && IM.detectFormat('a.sab', blockBin).id === 'sab' && IM.detectFormat('x.smb', smb).id === 'sab' && IM.detectFormat('p.ipt', ipt).id === 'inventor' && IM.detectFormat('d.f3d', f3d).id === 'fusion', 'detectFormat: SAT text, ACIS and ShapeManager binary signatures, Inventor and Fusion containers');
    const off = await IM.importFile({ name: 'block.sab', bytes: blockBin }, { kernel: false });
    ok(off.format.id === 'sab' && off.geometry === null && off.extra.meta.entityCounts.face === 6 && off.report.lost.some(s => /kernel was switched off/.test(s)) && off.report.validation.some(x => x.name === 'Geometry reconstructed' && !x.pass), 'importFile(sab, kernel off): record inventory only, and the report says why there is no geometry');
    const offT = await IM.importFile({ name: 'block.sat', bytes: blockNew }, { kernel: false }); ok(offT.geometry && offT.geometry.kind === 'lattice' && offT.geometry.wireframe && offT.report.assumptions.some(s => /kernel was switched off/.test(s)), 'importFile(sat, kernel off): falls back to the vertex-and-chord wireframe, stated as such');
    let bad = null, threw = ''; try { bad = await IM.importFile({ name: 'broken.sab', bytes: blockBin.subarray(0, 300) }); } catch (err) { threw = err.message; } ok(!threw && bad && bad.geometry === null, 'importFile(sab, truncated): falls back to the inventory without throwing');
    threw = ''; try { bad = await IM.importFile({ name: 'junk.ipt', bytes: cfbStorage('RSeStorage', [{ name: 'Mabc', data: new Uint8Array(80) }, { name: 'Babc', data: new Uint8Array(200).fill(0x78) }]) }); } catch (err) { threw = err.message; } ok(!threw && bad && bad.geometry === null && bad.report.lost.some(s => /no ShapeManager B-rep stream/.test(s)), 'importFile(Inventor without a readable B-rep): inventory, with the reason');
    const asmNo = await IM.importFile({ name: 'top.SLDASM', bytes: swFile([['swXmlContents/COMPINSTANCETREE', enc('<swSolidWorks><swHeader><swFile id="3" swDocType="ASSEMBLY" swPath="C:\\a\\top.SLDASM"/><swFile id="6" swDocType="PART" swPath="C:\\a\\blk.SLDPRT"/></swHeader><swModelList><swModel id="5" swName="blk" swFileRef="6"/><swModel id="2" swName="top" swFileRef="3"><swReference swName="blk" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="5" swTransform="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/></swModel></swModelList></swSolidWorks>')]]) });
    ok(asmNo.geometry === null && asmNo.extra.assembly.components.length === 1 && asmNo.report.lost.some(s => /import them together with the assembly/.test(s)), 'importFile(SLDASM alone): the referenced part files are named and asked for');
  }

  // ---- end to end through OpenCASCADE
  if (!have) log.push('note: OpenCASCADE kernel pack absent — ACIS → STEP → mesh checks skipped');
  else {
    const vol = A * B * C, area = 2 * (A * B + B * C + A * C);
    const rb = await IM.importFile({ name: 'block.sat', bytes: blockNew }), g = rb.geometry, pb = (rb.extra.parts || [])[0];
    ok(g && g.kind === 'mesh' && g.units === 'mm' && pb && pb.closed && near(pb.volume, vol, 1e-9) && near(pb.area, area, 1e-9) && rb.extra.acis.translated === 6 && bbox(g).every((x, i) => near(x, [A, B, C][i], 1e-9)), 'importFile(block.sat): closed solid of 40 × 30 × 20 mm with exact volume, area and bounding box');
    ok(rb.report.assumptions.some(s => /ACIS \/ ShapeManager B-rep translated to STEP and tessellated by OpenCASCADE/.test(s) && /6 of 6 faces translated/.test(s)) && rb.report.assumptions.some(s => /1 mm per unit/.test(s)) && rb.report.lost.some(s => /feature history, sketches, parameters and constraints are not recovered/.test(s)) && rb.report.provenance.geometryPath.startsWith('kernel') && (rb.extra.embedded || []).some(e => e.format === 'step' && /ADVANCED_FACE/.test(F.toText(e.bytes))), 'importFile(block.sat): the report names the path, the coverage, the unit factor and what is not recovered; the STEP text is kept');
    for (const [nm, bytes] of [['4.0 text', blockOld], ['binary', blockBin], ['explicit indices', satText(records(block(A, B, C)), { explicit: true })], ['spline top face and B-spline edge', satText(records(block(A, B, C, { spline: true })))]]) { const r = await IM.importFile({ name: nm.startsWith('bin') ? 'b.sab' : 'b.sat', bytes }); ok(r.geometry && r.extra.parts[0].closed && near(r.extra.parts[0].volume, vol, 1e-9), `importFile(block, ${nm}): the same closed solid`); }
    const fine = { tessellation: { linearDeflection: 0.0002, angularDeflection: 0.1 } }, mid = { tessellation: { linearDeflection: 0.001, angularDeflection: 0.2 } };
    const rc = await IM.importFile({ name: 'cyl.sat', bytes: cylNew }, fine), pc = (rc.extra.parts || [])[0], cv = PI * Rr * Rr * H, ca = 2 * PI * Rr * (Rr + H);
    ok(pc && pc.closed && near(pc.volume, cv, 4e-3) && pc.volume < cv && near(pc.area, ca, 3e-3) && rc.extra.acis.translated === 3, 'importFile(cyl.sat): closed cylinder, volume and area within the tessellation tolerance of π r² h and 2 π r (r + h)');
    // frusta: radius growing (sine and cosine of equal sign) and shrinking (opposite signs) along the axis
    const fv = (r0, r1) => PI * H * (r0 * r0 + r0 * r1 + r1 * r1) / 3, ru = await IM.importFile({ name: 'up.sat', bytes: satText(records(cylinder(Rr, H, 0.2))) }, mid), rd = await IM.importFile({ name: 'down.sat', bytes: satText(records(cylinder(Rr, H, -0.1))) }, mid), pu = (ru.extra.parts || [])[0], pd = (rd.extra.parts || [])[0];
    ok(pu && pu.closed && near(pu.volume, fv(Rr, Rr + 0.2 * H), 1e-2) && pd && pd.closed && near(pd.volume, fv(Rr, Rr - 0.1 * H), 1e-2) && count(F.toText(ru.extra.embedded.find(e => e.format === 'step').bytes), /=CONICAL_SURFACE\(/g) === 1, 'importFile(cone frusta): widening and narrowing cones give the frustum volume π h (r₀² + r₀ r₁ + r₁²) / 3');
    const rs = await IM.importFile({ name: 'ball.sat', bytes: satText(records(sphere(7))) }), ps = (rs.extra.parts || [])[0], rt = await IM.importFile({ name: 'ring.sat', bytes: satText(records(torus(12, 3))) }), ptor = (rt.extra.parts || [])[0];
    ok(ps && ps.closed && near(ps.volume, 4 / 3 * PI * 343, 1e-2) && ptor && ptor.closed && near(ptor.volume, 2 * PI * PI * 12 * 9, 1e-2), 'importFile(sphere, torus): faces without loops close into solids of volume 4/3 π r³ and 2 π² R r²');
    log.push(`ACIS cylinder r = 10 mm, h = 50 mm: volume ${pc ? pc.volume.toFixed(2) : '—'} vs ${cv.toFixed(2)} mm³; frusta ${pu ? pu.volume.toFixed(1) : '—'} vs ${fv(Rr, Rr + 0.2 * H).toFixed(1)} and ${pd ? pd.volume.toFixed(1) : '—'} vs ${fv(Rr, Rr - 0.1 * H).toFixed(1)} mm³; sphere ${ps ? ps.volume.toFixed(1) : '—'} vs ${(4 / 3 * PI * 343).toFixed(1)}; torus ${ptor ? ptor.volume.toFixed(1) : '—'} vs ${(2 * PI * PI * 108).toFixed(1)} mm³`);
    const rin = await IM.importFile({ name: 'inch.sat', bytes: satText(records(block(1, 2, 3)), { scale: 25.4 }) }); ok(rin.geometry && near(rin.extra.parts[0].volume, 6 * 25.4 ** 3, 1e-9) && rin.report.assumptions.some(s => /25\.4 mm per unit/.test(s)), 'importFile(inch file): millimetres through the header factor 25.4');
    // Inventor compound file and Fusion archive (lengths in centimetres)
    const ri = await IM.importFile({ name: 'block.ipt', bytes: ipt }); ok(ri.format.id === 'inventor' && ri.geometry && near(ri.extra.parts[0].volume, vol, 1e-9) && ri.report.assumptions.some(s => /PmBRepSegment/.test(s)) && ri.report.lost.some(s => /roll-back history/.test(s)) && ri.report.lost.some(s => /feature tree/.test(s)), 'importFile(Inventor): the B-rep segment\'s body is translated and tessellated (cm → mm); history and feature tree reported as not read');
    const rf = await IM.importFile({ name: 'block.f3d', bytes: f3d }); ok(rf.format.id === 'fusion' && rf.geometry && rf.extra.parts.length === 1 && near(rf.extra.parts[0].volume, vol, 1e-9) && rf.report.assumptions.some(s => /centimetres/.test(s)), 'importFile(Fusion): the history blob\'s current body is used, in centimetres whatever its header field says');
    // JT with an XT B-rep segment → Parasolid path
    const rj = await IM.importFile({ name: 'block.jt', bytes: jtFile(xtPart(0.04, 0.03, 0.02)) }); ok(rj.format.id === 'jt' && rj.geometry && near(rj.extra.parts[0].volume, vol, 1e-9) && rj.report.lost.some(s => /scene graph/.test(s)), 'importFile(JT with an XT B-rep segment): routed through the Parasolid path; the scene graph is reported as not decoded');
    // SOLIDWORKS part with two configurations
    const sw2 = swFile([['Contents/Config-0-Partition', partition(xtPart(0.04, 0.03, 0.02))], ['Contents/Config-1-Partition', partition(xtPart(0.04, 0.03, 0.01))], ['Contents/Config-2-Partition', partition(xtPart(0.04, 0.03, 0.02))]]), rc2 = await IM.importFile({ name: 'two.SLDPRT', bytes: sw2 }), cf = rc2.extra.configurations || [];
    ok(rc2.geometry && near(rc2.extra.parts[0].volume, vol, 1e-9) && cf.length === 3 && cf[0].active && cf[1].geometry && cf[1].geometry.kind === 'mesh' && near(cf[1].parts[0].volume, vol / 2, 1e-9) && cf[1].stats.translated === 6 && /configuration 0/.test(cf[2].sameAs || '') && rc2.report.imported.some(s => /^configurations —/.test(s)), 'importFile(SOLIDWORKS, three configurations): the first is the geometry, the second is tessellated under extra.configurations, an identical third is marked as such');
    // SOLIDWORKS assembly: two instances of the block, the second rotated a quarter turn about z and shifted, inside a sub-assembly
    const tree = '<swSolidWorks><swHeader><swFile id="3" swDocType="ASSEMBLY" swPath="C:\\a\\top.SLDASM"/><swFile id="6" swDocType="PART" swPath="C:\\a\\blk.SLDPRT"/><swFile id="9" swDocType="ASSEMBLY" swPath="C:\\a\\sub.SLDASM"/></swHeader><swModelList><swModel id="5" swName="blk" swFileRef="6"/><swModel id="8" swName="sub" swFileRef="9"><swReference swName="blk" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="5" swTransform="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0.05 1"/></swModel><swModel id="2" swName="top" swFileRef="3" swBoundingBox="0 0 0 0.1 0.04 0.07"><swReference swName="blk" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="5" swTransform="1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"/><swReference swName="sub" swReferenceNumber="1" swConfigurationId="0" swSuppressed="NO" swHidden="NO" swModelRef="8" swTransform="0 1 0 0 -1 0 0 0 0 0 1 0 0.1 0 0 1"/></swModel></swModelList><swConfigurationList><swConfiguration swName="Default" swMostRecentConfiguration="YES" swModelRef="2"/></swConfigurationList></swSolidWorks>';
    const part = swFile([['Contents/Config-0-Partition', partition(xtPart(0.04, 0.03, 0.02))]]), ra = await IM.importFile({ name: 'top.SLDASM', bytes: swFile([['swXmlContents/COMPINSTANCETREE', enc(tree)]]) }, { siblings: [{ name: 'BLK.sldprt', bytes: part }, { name: 'other.SLDPRT', bytes: part }] }), ga = ra.geometry;
    const ext = ga ? [0, 1, 2].map(a => [Math.min(...ga.vertices.map(q => q[a])), Math.max(...ga.vertices.map(q => q[a]))]) : [];
    ok(ga && ra.extra.parts.length === 2 && ra.extra.assemblyResult.placed === 2 && ra.extra.assemblyResult.bboxCheck.pass && near(ext[0][1], 100, 1e-9) && Math.abs(ext[0][0]) < 1e-9 && near(ext[1][1], 40, 1e-9) && near(ext[2][1], 70, 1e-9) && near(ra.extra.parts.reduce((s, q) => s + q.volume, 0), 2 * vol, 1e-9) && ra.report.assumptions.some(s => /placement check: the assembled instances lie inside/.test(s)), 'importFile(SLDASM + parts): instances are placed by the stored transforms (sub-assembly composed) and agree with the extents the assembly stores');
    const rm = await IM.importFile({ name: 'top.SLDASM', bytes: swFile([['swXmlContents/COMPINSTANCETREE', enc(tree.replace('blk.SLDPRT', 'gone.SLDPRT'))]]) }, { siblings: [{ name: 'blk.SLDPRT', bytes: part }] }); ok(rm.geometry === null && rm.report.lost.some(s => /none of which was supplied/.test(s)), 'importFile(SLDASM, parts absent): no geometry is invented; the missing files are named');
  }
  log.push(`sat self-test: ${pass} passed, ${fail} failed in ${((Date.now() - t00) / 1000).toFixed(1)} s (OpenCASCADE kernel ${have ? 'present' : 'absent'})`);
  return { pass, fail, log };
}
