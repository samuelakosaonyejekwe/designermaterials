// Self-test of the Parasolid XT path: the node reader (text and neutral binary, base schema and embedded schema
// differences), the XT → STEP translator and the importer wiring (transmit files, a SOLIDWORKS partition, an NX
// container). Default export: async () => {pass, fail, log}. Every fixture is synthesised here from the published
// transmit layout, so the test runs offline; the tessellation checks are skipped (with a logged note) when the
// OpenCASCADE kernel pack is absent.
import * as F from '../js/geometry/formats.js';
import * as K from '../js/geometry/kernels.js';
import * as N from '../js/geometry/native.js';
import * as X from '../js/geometry/xt.js';
import * as IM from '../js/geometry/importer.js';

const enc = s => Uint8Array.from(s, c => c.charCodeAt(0)), PI = Math.PI;
const cat = (...a) => { const o = new Uint8Array(a.reduce((n, x) => n + x.length, 0)); let p = 0; for (const x of a) { o.set(x, p); p += x.length; } return o; };
const u32 = (...v) => { const o = new Uint8Array(4 * v.length), d = new DataView(o.buffer); v.forEach((x, i) => d.setUint32(4 * i, x >>> 0, true)); return o; };
const stored = d => { const parts = []; for (let i = 0; i < d.length || !parts.length; i += 65535) { const c = d.subarray(i, i + 65535), last = i + 65535 >= d.length; parts.push(Uint8Array.from([last ? 1 : 0, c.length & 255, c.length >> 8, ~c.length & 255, (~c.length >> 8) & 255]), c); } return cat(...parts); };
const zlibStored = d => cat(Uint8Array.from([0x78, 0x01]), stored(d), u32(0));
const swFile = streams => cat(u32(0x12345678), Uint8Array.from([0, 0, 0, 4]), new Uint8Array(12), ...streams.flatMap(([name, data]) => { const d = stored(data), nm = Uint8Array.from(enc(name), c => ((c << 4) & 0xF0) | (c >> 4)); return [Uint8Array.from([0x14, 0, 6, 0, 8, 0]), u32(0xAABBCCDD, F.crc32(data), d.length, data.length, nm.length), nm, d]; }));

// typed field values of a node record: pointer, int, byte, double (null = unset), char, vector
const p = n => ['p', n], d = n => ['d', n], u = n => ['u', n], f = x => ['f', x], c = s => ['c', s], v = a => ['v', a];
const HEAD = (fmt, sch) => `**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**************************\n**PARASOLID !"#$%&'()*+,-./:;<=>?@[\\]^_\`{|}~0123456789**************************\n**PART1;\nMC=test;\nAPPL=fixture;\nFORMAT=${fmt};\nGUISE=transmit;\nDATE=1-jan-2024;\n**PART2;\nSCH=${sch};\nUSFLD_SIZE=0;\n**PART3;\n**END_OF_HEADER*****************************************************************\n`;
const MOD = partition => `: TRANSMIT FILE ${partition ? '(partition) ' : ''}created by modeller version 1200000`, SCH = 'SCH_1200000_12006';
// text transmit: numbers are followed by one space, chars are not, an unset value is "?"; records are wrapped at 80 columns
function xtText(nodes, partition = false) {
  const num = x => (x === null ? '?' : String(x) + ' '); let s = `T${MOD(partition).length} ${MOD(partition)}${SCH.length} ${SCH}0 `;
  for (const [type, index, ...vals] of nodes) { s += `${type} ${index} `; for (const [k, x] of vals) s += k === 'c' ? x : k === 'v' ? (x ? x.map(num).join('') : '?') : num(x); }
  s += '1 0 '; let out = ''; for (let i = 0; i < s.length; i += 79) out += s.slice(i, i + 79) + '\n';
  return enc(HEAD('text', SCH) + out);
}
// neutral binary transmit: big-endian, pointers and node indices as index + 1 in two bytes
function xtBinary(nodes, partition = false, textHeader = true) {
  const out = [], i16 = n => out.push((n >> 8) & 255, n & 255), i32 = n => out.push((n >>> 24) & 255, (n >>> 16) & 255, (n >>> 8) & 255, n & 255), f64 = x => { const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, x === null ? -3.14158e13 : x, false); out.push(...b); };
  out.push(0x50, 0x53, 0, 0); i16(MOD(partition).length); out.push(...enc(MOD(partition))); i32(SCH.length); out.push(...enc(SCH)); i32(0);
  for (const [type, index, ...vals] of nodes) { i16(type); i16(index + 1); for (const [k, x] of vals) { if (k === 'p') i16(x + 1); else if (k === 'd') i32(x); else if (k === 'u') out.push(x); else if (k === 'f') f64(x); else if (k === 'c') out.push(x.charCodeAt(0)); else for (const q of x || [null, null, null]) f64(q); } }
  i16(1); i16(1); return textHeader ? cat(enc(HEAD('binary', SCH)), Uint8Array.from(out)) : Uint8Array.from(out);
}
const curveHead = (id, sense = '+') => [d(id), p(0), p(0), p(0), p(0), p(0), c(sense)];
// a solid from faces [{surface: [type, ...fields], sense, loops: [[{edge, sense, vertex}]]}], edges [{curve: [type, ...fields]}], points
function solid({ faces, edges, points }, partition) {
  const nodes = []; let next = partition ? 2 : 1; const id = () => next++, BODY = id(), VOID = id(), SOLID = id(), SHELL = id();
  const vx = points.map(() => id()), pt = points.map(() => id()), ed = edges.map(() => id()), cu = edges.map(() => id()), fa = faces.map(() => id()), su = faces.map(() => id());
  const lp = faces.map(fc => fc.loops.map(() => id())), fin = faces.map(fc => fc.loops.map(l => l.map(() => id()))), uses = edges.map(() => []);
  faces.forEach((fc, i) => fc.loops.forEach((l, j) => l.forEach((h, k) => uses[h.edge].push(fin[i][j][k]))));
  if (partition) nodes.push([101, 1, p(0), p(0), p(BODY), p(0), p(0), p(0), p(0), ['c', 'T'], p(0), d(0), d(0)]);
  nodes.push([12, BODY, d(next), p(0), p(0), p(0), p(0), p(0), p(0), f(1000), f(1e-8), p(0), p(0), p(0), u(1), p(0), u(1), u(1), p(SHELL), p(su[0]), p(cu[0] || 0), p(pt[0] || 0), p(VOID), p(ed[0] || 0), p(vx[0] || 0)]);
  nodes.push([19, VOID, d(1), p(0), p(BODY), p(SOLID), p(0), p(0), c('V')], [19, SOLID, d(2), p(0), p(BODY), p(0), p(VOID), p(SHELL), c('S')], [13, SHELL, d(3), p(0), p(BODY), p(0), p(fa[0]), p(0), p(0), p(SOLID), p(0)]);
  points.forEach((q, i) => nodes.push([18, vx[i], d(100 + i), p(0), p(0), p(0), p(0), p(pt[i]), f(null), p(BODY)], [29, pt[i], d(200 + i), p(0), p(vx[i]), p(0), p(0), v(q)]));
  edges.forEach((e, i) => nodes.push([16, ed[i], d(300 + i), p(0), f(null), p(uses[i][0]), p(0), p(0), p(cu[i]), p(0), p(0), p(BODY)], [e.curve[0], cu[i], ...curveHead(400 + i), ...e.curve.slice(1)]));
  faces.forEach((fc, i) => {
    nodes.push([14, fa[i], d(500 + i), p(0), f(null), p(fa[i + 1] || 0), p(fa[i - 1] || 0), p(lp[i][0] || 0), p(SHELL), p(su[i]), c(fc.sense || '+'), p(0), p(0), p(0), p(0), p(0)], [fc.surface[0], su[i], ...curveHead(600 + i), ...fc.surface.slice(1)]);
    fc.loops.forEach((l, j) => { nodes.push([15, lp[i][j], d(700 + 10 * i + j), p(0), p(fin[i][j][0]), p(fa[i]), p(lp[i][j + 1] || 0)]); l.forEach((h, k) => { const us = uses[h.edge], me = fin[i][j][k]; nodes.push([17, me, p(0), p(lp[i][j]), p(fin[i][j][(k + 1) % l.length]), p(fin[i][j][(k + l.length - 1) % l.length]), p(h.vertex === undefined ? 0 : vx[h.vertex]), p(us[(us.indexOf(me) + 1) % us.length]), p(ed[h.edge]), p(0), p(0), c(h.sense)]); }); });
  });
  return nodes;
}
// block [0, a] × [0, b] × [0, c]: six planes with outward normals, twelve straight edges, eight vertices
function block(a, b, cz) {
  const size = [a, b, cz], points = [], edges = [], key = new Map(), faces = [], unitv = k => [0, 1, 2].map(i => (i === k ? 1 : 0));
  for (let i = 0; i < 8; i++) points.push([i & 1 ? a : 0, i & 2 ? b : 0, i & 4 ? cz : 0]);
  const edge = (m, n) => { const lo = Math.min(m, n), hi = Math.max(m, n), k = lo * 8 + hi; if (!key.has(k)) { const q = points[hi].map((x, i) => x - points[lo][i]), l = Math.hypot(...q); key.set(k, edges.length); edges.push({ curve: [30, v(points[lo]), v(q.map(x => x / l))] }); } return { edge: key.get(k), sense: m < n ? '+' : '-', vertex: n }; };
  for (let k = 0; k < 3; k++) for (const s of [0, 1]) { const uu = (k + 1) % 3, vv = (k + 2) % 3, at = (i, j) => (s << k) | (i << uu) | (j << vv), ring = s ? [at(0, 0), at(1, 0), at(1, 1), at(0, 1)] : [at(0, 0), at(0, 1), at(1, 1), at(1, 0)], o = [0, 0, 0]; o[k] = s ? size[k] : 0;
    faces.push({ surface: [50, v(o), v(unitv(k).map(x => (s ? x : -x))), v(unitv(uu))], loops: [ring.map((m, i) => edge(m, ring[(i + 1) % 4]))] }); }
  return { faces, edges, points };
}
// cylinder of radius r and height h about z: one cylindrical and two planar faces, two vertex-free circular edges
const cylinder = (r, h) => ({ points: [], edges: [{ curve: [31, v([0, 0, 0]), v([0, 0, 1]), v([1, 0, 0]), f(r)] }, { curve: [31, v([0, 0, h]), v([0, 0, 1]), v([1, 0, 0]), f(r)] }],
  faces: [{ surface: [51, v([0, 0, 0]), v([0, 0, 1]), f(r), v([1, 0, 0])], loops: [[{ edge: 0, sense: '+' }], [{ edge: 1, sense: '-' }]] }, { surface: [50, v([0, 0, 0]), v([0, 0, -1]), v([1, 0, 0])], loops: [[{ edge: 0, sense: '-' }]] }, { surface: [50, v([0, 0, h]), v([0, 0, 1]), v([1, 0, 0])], loops: [[{ edge: 1, sense: '+' }]] }] });

// fixture builders shared with the other container tests (sat.test.mjs)
export const fixtures = { enc, cat, u32, stored, zlibStored, swFile, xtBinary, xtText, solid, block, cylinder };

export default async function run() {
  const log = []; let pass = 0, fail = 0; const t00 = Date.now();
  const ok = (cond, msg) => { if (cond) pass++; else { fail++; log.push('✗ ' + msg); } };
  const near = (a, b, rel) => Math.abs(a - b) <= rel * Math.abs(b);
  const have = await K.kernelAvailable('occt'), A = 0.04, B = 0.03, C = 0.02, R = 0.01, H = 0.05;
  const blockT = xtText(solid(block(A, B, C))), blockB = xtBinary(solid(block(A, B, C))), cylT = xtText(solid(cylinder(R, H)));

  // ---- reader: text and neutral binary give the same nodes
  {
    const mt = X.readXT(blockT), mb = X.readXT(blockB);
    ok(mt.encoding === 'text' && mb.encoding === 'binary' && mt.schema === 'SCH_1200000_12006' && mt.modellerVersion === '1200000' && !mt.partition && mt.header.appl === 'fixture' && mt.root.name === 'BODY', 'readXT: header, modeller version, schema and root node of a text and a binary transmit');
    ok(['FACE', 'EDGE', 'VERTEX', 'FIN', 'LOOP', 'PLANE', 'LINE', 'POINT'].every(k => mt.counts[k] === mb.counts[k]) && mt.counts.FACE === 6 && mt.counts.EDGE === 12 && mt.counts.VERTEX === 8 && mt.counts.FIN === 24 && mt.nodes.size === mb.nodes.size && !mt.trailing && !mb.trailing, 'readXT: text and binary encodings decode to the same node inventory (6 faces, 12 edges, 8 vertices, 24 fins)');
    const same = [...mt.nodes.values()].every(n => JSON.stringify(n) === JSON.stringify(mb.nodes.get(n.index))); ok(same, 'readXT: every field of every node agrees between the two encodings');
    const pt = [...mt.nodes.values()].filter(n => n.name === 'POINT').map(n => n.pvec), e = [...mt.nodes.values()].find(n => n.name === 'EDGE');
    ok(pt.some(q => q[0] === A && q[1] === B && q[2] === C) && e.tolerance === null && mt.nodes.get(e.curve).name === 'LINE' && mt.nodes.get(e.curve).sense === '+', 'readXT: vectors, unset ("?") values, pointers and sense characters');
    // embedded schema: POINT edited by Copy / Insert / Append instructions, CYLINDER declared identical (255), and an unknown node type described in full
    const emb = HEAD('text', 'SCH_2700142_26105') + `T51 : TRANSMIT FILE created by modeller version 270014223 SCH_2700142_26105_13006196 0 12 255 1 9 0 0 0 0 0 0 1e3 1e-8 0 0 0 1 0 1 1 0 0 0 0 0 0 0 29 8 CCI5 extra0 0 1 dCCCCA4 tail29 0 Z2 7 0 42 1 0 0 1 2 3 0 51 255 3 5 0 0 0 0 0 +0 0 0 0 0 1 .25 1 0 0 200 2 4 demo9 test node3 val0 0 1 f4 text0 1 1 cT4 4 1.5 abcd1 0 `;
    const me = X.readXT(enc(emb)), q = me.nodes.get(2), cy = me.nodes.get(3), ex = me.nodes.get(4);
    ok(me.baseSchema === '13006' && q.name === 'POINT' && q.extra === 42 && q.owner === 1 && q.pvec.join() === '1,2,3' && q.tail === 0 && cy.name === 'CYLINDER' && cy.radius === 0.25 && cy.sense === '+', 'readXT: embedded schema differences (copy, insert, append; unchanged flag 255) drive the field layout');
    ok(ex && ex.name === 'DEMO' && ex.val === 1.5 && ex.text === 'abcd' && me.extensions.includes('DEMO'), 'readXT: a node type the base schema lacks is read from the description the file carries');
    let msg = ''; try { X.readXT(blockT.subarray(0, blockT.length - 300)); } catch (err) { msg = err.message; } ok(/end of the transmit data|expected/.test(msg), 'readXT: a truncated file is rejected with an explicit message');
    msg = ''; try { X.readXT(enc('hello')); } catch (err) { msg = err.message; } ok(/Parasolid/.test(msg), 'readXT: data without a transmit flag sequence is rejected');
  }

  // ---- translator
  {
    const tb = X.xtToStep(X.readXT(blockT), { name: 'block.x_t' }), tc = X.xtToStep(X.readXT(cylT), { name: 'cyl.x_t' }), n = (s, re) => (s.match(re) || []).length;
    ok(tb.stats.faces === 6 && tb.stats.translated === 6 && tb.stats.exactEdges === 12 && tb.stats.approximatedEdges === 0 && tb.bodies.length === 1 && tb.bodies[0].closed && tb.bodies[0].type === 'solid', 'xtToStep(block): 6 of 6 faces and 12 exact edges translated into one closed solid');
    ok(/^ISO-10303-21;/.test(tb.step) && /AUTOMOTIVE_DESIGN/.test(tb.step) && n(tb.step, /=MANIFOLD_SOLID_BREP\(/g) === 1 && n(tb.step, /=CLOSED_SHELL\(/g) === 1 && n(tb.step, /=ADVANCED_FACE\(/g) === 6 && n(tb.step, /=PLANE\(/g) === 6 && n(tb.step, /=EDGE_CURVE\(/g) === 12 && n(tb.step, /=ORIENTED_EDGE\(/g) === 24 && n(tb.step, /=VERTEX_POINT\(/g) === 8 && n(tb.step, /=LINE\(/g) === 12 && /END-ISO-10303-21;\s*$/.test(tb.step), 'xtToStep(block): AP214 text with the expected entity counts');
    ok(/SI_UNIT\(\$,\.METRE\.\)/.test(tb.step) && /UNCERTAINTY_MEASURE_WITH_UNIT\(LENGTH_MEASURE\(1\.E-8\)/.test(tb.step) && /ADVANCED_BREP_SHAPE_REPRESENTATION/.test(tb.step) && /PRODUCT\('block'/.test(tb.step), 'xtToStep(block): metres declared, uncertainty taken from the body resolution, product named after the file');
    ok(tc.stats.translated === 3 && n(tc.step, /=CYLINDRICAL_SURFACE\(/g) === 1 && n(tc.step, /=CIRCLE\(/g) === 2 && n(tc.step, /=FACE_BOUND\(/g) === 2 && n(tc.step, /=FACE_OUTER_BOUND\(/g) === 2 && n(tc.step, /=VERTEX_POINT\(/g) === 2, 'xtToStep(cylinder): vertex-free circular edges get a vertex on the curve; the cylindrical face keeps its two bounds');
    // a face whose surface has no STEP form is left out and the body becomes an open shell
    let swapped = false; const nodes = solid(block(A, B, C)).map(x => (x[0] === 50 && !swapped ? (swapped = true, [59, x[1], ...x.slice(2, 9), d(0), p(0)]) : x)); // one plane → a blend-boundary surface (boundary: short, blend: pointer)
    const to = X.xtToStep(X.readXT(xtText(nodes)));
    ok(to.stats.faces === 6 && to.stats.translated === 5 && to.stats.skipped['BLEND_BOUND surface'] === 1 && !to.bodies[0].closed && /=OPEN_SHELL\(/.test(to.step) && /SHELL_BASED_SURFACE_MODEL/.test(to.step) && !/MANIFOLD_SOLID_BREP/.test(to.step), 'xtToStep: an untranslatable surface costs exactly its face (5 of 6), reported by type, and the rest is written as an open shell');
    ok(X.xtToStep(X.readXT(enc(HEAD('text', SCH) + 'T51 : TRANSMIT FILE created by modeller version 120000017 SCH_1200000_120060 1 0 '))).step === null, 'xtToStep: a transmit without bodies yields no STEP text');
  }

  // ---- containers: Parasolid streams located by signature
  {
    const part = xtBinary(solid(block(A, B, C), true), true, false), junk = new Uint8Array(700).map((_, i) => (i * 37) & 255), nx = cat(enc('SPLMSSTR'), junk, enc('/Root/UG_PART/\0'), zlibStored(part), junk, part);
    const found = N.findParasolid(nx); ok(found.length === 2 && found[0].compressed && found[0].partition && found[0].bytes.length === part.length && !found[1].compressed && found[1].schema === 'SCH_1200000_12006', 'findParasolid: deflated and raw neutral-binary transmit streams are located inside another file');
    const r = N.parseNX(nx); ok(r.extra.embedded.length === 2 && r.extra.embedded[0].format === 'parasolid' && /SPLM/.test(r.meta.layout) && r.geometry === null && r.notRecovered.some(s => /not interpreted/.test(s)), 'parseNX: container recognised and its Parasolid partitions offered');
    ok(IM.detectFormat('a.prt', nx).id === 'nx' && IM.detectFormat('b.x_b', blockB).id === 'parasolid' && IM.detectFormat('b.x_b', blockB).confidence > 0.9 && IM.detectFormat('c.xmt_txt', blockT).id === 'parasolid', 'detectFormat: NX container and Parasolid transmit files (text header before binary data included)');
  }

  // ---- importer: registry, fallback without the kernel, never throwing
  {
    const by = Object.fromEntries(IM.FORMATS.map(x => [x.id, x]));
    ok(by.parasolid.support === 'partial' && by.solidworks.support === 'partial' && by.nx.support === 'partial' && /translated/.test(by.parasolid.note) && /faces/.test(by.parasolid.note) && /not/.test(by.parasolid.note), 'registry: Parasolid, SOLIDWORKS and NX are "partial" with the measured coverage in the note');
    const off = await IM.importFile({ name: 'block.x_t', bytes: blockT }, { kernel: false });
    ok(off.format.id === 'parasolid' && off.geometry === null && off.extra.meta.schema === 'SCH_1200000_12006' && off.report.lost.some(s => /kernel was switched off/.test(s)) && off.report.validation.some(x => x.name === 'Geometry reconstructed' && !x.pass), 'importFile(x_t, kernel off): header inventory only, and the report says why there is no geometry');
    let bad = null, threw = ''; try { bad = await IM.importFile({ name: 'broken.x_t', bytes: blockT.subarray(0, blockT.length - 400) }); } catch (err) { threw = err.message; }
    ok(!threw && bad && bad.geometry === null && bad.report.lost.some(s => /node records could not be translated/.test(s)), 'importFile(x_t, truncated): falls back to the inventory without throwing');
  }

  // ---- end to end through OpenCASCADE
  if (!have) log.push('note: OpenCASCADE kernel pack absent — XT → STEP → mesh checks skipped');
  else {
    const rb = await IM.importFile({ name: 'block.x_t', bytes: blockT }), g = rb.geometry, pb = (rb.extra.parts || [])[0], vol = A * B * C * 1e9, area = 2 * (A * B + B * C + A * C) * 1e6;
    ok(g && g.kind === 'mesh' && g.units === 'mm' && pb && pb.closed && near(pb.volume, vol, 1e-9) && near(pb.area, area, 1e-9) && rb.extra.parasolid.translated === 6 && rb.extra.parasolid.faces === 6, 'importFile(block.x_t): closed solid of 40 × 30 × 20 mm with exact volume and area');
    const [lo, hi] = [0, 1, 2].map(a => [Math.min(...g.vertices.map(q => q[a])), Math.max(...g.vertices.map(q => q[a]))]).reduce((s, x) => [[...s[0], x[0]], [...s[1], x[1]]], [[], []]);
    ok(near(hi[0] - lo[0], 40, 1e-9) && near(hi[1] - lo[1], 30, 1e-9) && near(hi[2] - lo[2], 20, 1e-9), 'importFile(block.x_t): bounding box in millimetres (metres in the file)');
    ok(rb.report.assumptions.some(s => /Parasolid B-rep translated to STEP and tessellated by OpenCASCADE/.test(s) && /6 of 6 faces translated/.test(s)) && rb.report.assumptions.some(s => /metres/.test(s) && /millimetres/.test(s)) && rb.report.lost.some(s => /feature history, sketches, PMI and assembly constraints are not recovered/.test(s)) && rb.report.provenance.geometryPath.startsWith('kernel') && (rb.extra.embedded || []).some(e => e.format === 'step' && /ISO-10303-21/.test(F.toText(e.bytes.subarray(0, 40)))), 'importFile(block.x_t): report states the path, coverage, units and what is not recovered; the STEP text is offered');
    const rbin = await IM.importFile({ name: 'block.x_b', bytes: blockB }); ok(rbin.geometry && near(rbin.extra.parts[0].volume, vol, 1e-9), 'importFile(block.x_b): neutral binary gives the same solid');
    const rc = await IM.importFile({ name: 'cyl.x_t', bytes: cylT }, { tessellation: { linearDeflection: 0.0002, angularDeflection: 0.1 } }), pc = (rc.extra.parts || [])[0], cv = PI * R * R * H * 1e9, ca = (2 * PI * R * R + 2 * PI * R * H) * 1e6;
    ok(pc && pc.closed && near(pc.volume, cv, 4e-3) && pc.volume < cv && near(pc.area, ca, 3e-3) && rc.extra.parasolid.translated === 3, 'importFile(cyl.x_t): closed cylinder, volume and area within the tessellation tolerance of π r² h and 2 π r (r + h)');
    log.push(`XT cylinder r = 10 mm, h = 50 mm: volume ${pc ? pc.volume.toFixed(2) : '—'} vs ${cv.toFixed(2)} mm³, area ${pc ? pc.area.toFixed(2) : '—'} vs ${ca.toFixed(2)} mm²; block 40 × 30 × 20 mm: ${pb ? pb.volume.toFixed(3) : '—'} mm³`);
    // SOLIDWORKS 2015+ container holding the block as a Parasolid partition; NX container holding the same partition
    const part = xtBinary(solid(block(A, B, C), true), true, false), sw = swFile([['Contents/Config-0-Partition', cat(u32(1), new Uint8Array(16), u32(part.length, 0), zlibStored(part))]]);
    const rs = await IM.importFile({ name: 'block.SLDPRT', bytes: sw }); ok(rs.format.id === 'solidworks' && rs.geometry && near(rs.extra.parts[0].volume, vol, 1e-9) && rs.report.assumptions.some(s => /partition of configuration 0/.test(s)) && rs.report.lost.some(s => /feature tree/.test(s)), 'importFile(SOLIDWORKS): the embedded partition is translated and tessellated; the feature tree is reported as not interpreted');
    const rn = await IM.importFile({ name: 'block.prt', bytes: cat(enc('SPLMSSTR'), new Uint8Array(64), zlibStored(part)) }); ok(rn.format.id === 'nx' && rn.geometry && near(rn.extra.parts[0].volume, vol, 1e-9) && rn.report.lost.some(s => /NX object model/.test(s)), 'importFile(NX): the partition located in the container is translated and tessellated');
  }
  log.push(`xt self-test: ${pass} passed, ${fail} failed in ${((Date.now() - t00) / 1000).toFixed(1)} s (OpenCASCADE kernel ${have ? 'present' : 'absent'})`);
  return { pass, fail, log };
}
