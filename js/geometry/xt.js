// Parasolid XT transmit files (.x_t / .x_b, and the partition streams embedded in SOLIDWORKS documents) read node by
// node from the schema the file declares, and translated entity by entity into an ISO 10303-21 (STEP AP214) advanced
// B-rep text. Nothing is evaluated or tessellated here: trimming, healing and meshing are left to the OpenCASCADE
// kernel that reads the STEP text. Field layouts follow the Parasolid XT Format Reference (base schema 13006 plus the
// schema differences each file embeds).
const { abs, min, max, hypot, cos, sin, atan2, sqrt, PI, floor } = Math;

// ---------------------------------------------------------------- base schema (SCH_13006)
// type|NAME|field:code[:count] …   count: omitted = scalar, 1 = variable length, n > 1 = array of n.
// codes: d int, n short, u byte, l logical, c char, w unicode char, f double, p pointer, v vector, i interval, b box, h hvec, t tag
const CU = 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c ';
const BASE = `
10|ASSEMBLY|highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p
11|INSTANCE|node_id:d attributes_groups:p type:u part:p transform:p assembly:p next_in_part:p prev_in_part:p next_of_part:p prev_of_part:p
12|BODY|highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p
13|SHELL|node_id:d attributes_groups:p body:p next:p face:p edge:p vertex:p region:p front_face:p
14|FACE|node_id:d attributes_groups:p tolerance:f next:p previous:p loop:p shell:p surface:p sense:c next_on_surface:p previous_on_surface:p next_front:p previous_front:p front_shell:p
15|LOOP|node_id:d attributes_groups:p fin:p face:p next:p
16|EDGE|node_id:d attributes_groups:p tolerance:f fin:p previous:p next:p curve:p next_on_curve:p previous_on_curve:p owner:p
17|FIN|attributes_groups:p loop:p forward:p backward:p vertex:p other:p edge:p curve:p next_at_vx:p sense:c
18|VERTEX|node_id:d attributes_groups:p fin:p previous:p next:p point:p tolerance:f owner:p
19|REGION|node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c
29|POINT|node_id:d attributes_groups:p owner:p next:p previous:p pvec:v
30|LINE|${CU}pvec:v direction:v
31|CIRCLE|${CU}centre:v normal:v x_axis:v radius:f
32|ELLIPSE|${CU}centre:v normal:v x_axis:v major_radius:f minor_radius:f
38|INTERSECTION|${CU}surface:p:2 chart:p start:p end:p
40|CHART|base_parameter:f base_scale:f chart_count:d chordal_error:f angular_error:f parameter_error:f:2 hvec:h:1
41|LIMIT|type:c hvec:h:1
45|BSPLINE_VERTICES|vertices:f:1
50|PLANE|${CU}pvec:v normal:v x_axis:v
51|CYLINDER|${CU}pvec:v axis:v radius:f x_axis:v
52|CONE|${CU}pvec:v axis:v radius:f sin_half_angle:f cos_half_angle:f x_axis:v
53|SPHERE|${CU}centre:v radius:f axis:v x_axis:v
54|TORUS|${CU}centre:v axis:v major_radius:f minor_radius:f x_axis:v
56|BLENDED_EDGE|${CU}blend_type:c surface:p:2 spine:p range:f:2 thumb_weight:f:2 boundary:p:2 start:p end:p
59|BLEND_BOUND|${CU}boundary:n blend:p
60|OFFSET_SURF|${CU}check:c true_offset:l surface:p offset:f scale:f
67|SWEPT_SURF|${CU}section:p sweep:v scale:f
68|SPUN_SURF|${CU}profile:p base:v axis:v start:v end:v start_param:f end_param:f x_axis:v scale:f
70|LIST|node_id:d owner:p next:p previous:p list_type:d list_length:d block_length:d size_of_entry:d list_block:p finger_block:p finger_index:d notransmit:l
74|POINTER_LIS_BLOCK|n_entries:d next_block:p entries:p:1
79|ATT_DEF_ID|string:c:1
80|ATTRIB_DEF|next:p identifier:p type_id:d actions:u:8 field_names:p legal_owners:l:14 fields:u:1
81|ATTRIBUTE|node_id:d definition:p owner:p next:p previous:p next_of_type:p previous_of_type:p fields:p:1
82|INT_VALUES|values:d:1
83|REAL_VALUES|values:f:1
84|CHAR_VALUES|values:c:1
85|POINT_VALUES|values:v:1
86|VECTOR_VALUES|values:v:1
87|AXIS_VALUES|values:v:1
88|TAG_VALUES|values:t:1
89|DIRECTION_VALUES|values:v:1
90|GROUP|node_id:d attributes_groups:p owner:p next:p previous:p type:u first_member:p
91|MEMBER_OF_GROUP|dummy_node_id:d owning_group:p owner:p next:p previous:p next_member:p previous_member:p
98|UNICODE_VALUES|values:w:1
99|FIELD_NAMES|names:p:1
100|TRANSFORM|node_id:d owner:p next:p previous:p rotation_matrix:f:9 translation_vector:v scale:f flag:d perspective_vector:v
101|WORLD|assembly:p attribute:p body:p transform:p surface:p curve:p point:p alive:l attrib_def:p highest_id:d current_id:d
102|KEY|string:c:1
120|PE_SURF|${CU}type:c data:p tf:p internal_geom:p:1
121|INT_PE_DATA|geom_type:d real_array:p int_array:p
122|EXT_PE_DATA|key:p real_array:p int_array:p
124|B_SURFACE|${CU}nurbs:p data:p
125|SURFACE_DATA|original_uint:i original_vint:i extended_uint:i extended_vint:i self_int:u original_u_start:c original_u_end:c original_v_start:c original_v_end:c extended_u_start:c extended_u_end:c extended_v_start:c extended_v_end:c analytic_form_type:c swept_form_type:c spun_form_type:c blend_form_type:c analytic_form:p swept_form:p spun_form:p blend_form:p
126|NURBS_SURF|u_periodic:l v_periodic:l u_degree:n v_degree:n n_u_vertices:d n_v_vertices:d u_knot_type:u v_knot_type:u n_u_knots:d n_v_knots:d rational:l u_closed:l v_closed:l surface_form:u vertex_dim:n bspline_vertices:p u_knot_mult:p v_knot_mult:p u_knots:p v_knots:p
127|KNOT_MULT|mult:n:1
128|KNOT_SET|knots:f:1
130|PE_CURVE|${CU}type:c data:p tf:p internal_geom:p:1
133|TRIMMED_CURVE|${CU}basis_curve:p point_1:v point_2:v parm_1:f parm_2:f
134|B_CURVE|${CU}nurbs:p data:p
135|CURVE_DATA|self_int:u analytic_form:p
136|NURBS_CURVE|degree:n n_vertices:d vertex_dim:n n_knots:d knot_type:u periodic:l closed:l rational:l curve_form:u bspline_vertices:p knot_mult:p knots:p
137|SP_CURVE|${CU}surface:p b_curve:p original:p tolerance_to_original:f
141|GEOMETRIC_OWNER|owner:p next:p previous:p shared_geometry:p
163|HELIX_SU_FORM|axis_pt:v axis_dir:v hand:c turns:i pitch:f gap:f tol:f
184|HELIX_CU_FORM|axis_pt:v axis_dir:v point:v hand:c turns:i pitch:f tol:f
`;
const SCHEMA = new Map(BASE.trim().split('\n').map(l => { const [t, name, f] = l.split('|'); return [+t, { type: +t, name, fields: f.trim().split(/\s+/).map(s => { const [n, c, k] = s.split(':'); return { name: n, code: c, n: k ? +k : 0 }; }) }]; }));
export const NODE_NAMES = Object.fromEntries([...SCHEMA].map(([t, s]) => [t, s.name]));
const NULL_F = -3.14158e13, NULL_D = -32764;

// ---------------------------------------------------------------- token readers
class Bin {
  constructor(b, o) { this.b = b; this.v = new DataView(b.buffer, b.byteOffset, b.byteLength); this.o = o; }
  need(n) { if (this.o + n > this.b.length) throw new Error('unexpected end of the transmit data'); }
  byte() { this.need(1); return this.b[this.o++]; }
  short() { this.need(2); const x = this.v.getInt16(this.o, false); this.o += 2; return x; }
  int() { this.need(4); const x = this.v.getInt32(this.o, false); this.o += 4; return x; }
  dbl() { this.need(8); const x = this.v.getFloat64(this.o, false); this.o += 8; return x === NULL_F ? null : x; }
  ptr() { let r = this.short(), q = 0; if (r < 0) { q = this.short(); r = -r; } return q * 32767 + r - 1; }
  type() { return this.short(); }
  count() { return this.int(); }
  nfields() { return this.byte(); }
  ch() { return String.fromCharCode(this.byte()); }
  chars(n) { this.need(n); let s = ''; for (let i = 0; i < n; i++) s += String.fromCharCode(this.b[this.o + i]); this.o += n; return s; }
  sstr() { return this.chars(this.byte()); }
  logical() { return this.byte() !== 0; }
  pos() { return this.ptr(); }
  d() { const x = this.int(); return x === NULL_D ? null : x; }
  n() { const x = this.short(); return x === NULL_D ? null : x; }
  u() { return this.byte(); }
  vec() { const x = this.dbl(), y = this.dbl(), z = this.dbl(); return x === null || y === null || z === null ? null : [x, y, z]; }
  uf() { return this.int(); }
  done() { return this.o >= this.b.length; }
}
class Txt {
  constructor(s, o) { this.s = s; this.o = o; }
  num() {
    const s = this.s; let o = this.o; while (s.charCodeAt(o) === 32) o++;
    if (o >= s.length) throw new Error('unexpected end of the transmit data');
    if (s[o] === '?') { this.o = o + 1; return null; }
    let e = o; while (e < s.length && s.charCodeAt(e) !== 32) e++;
    const t = s.slice(o, e), x = Number(t); if (t === '' || Number.isNaN(x)) throw new Error(`"${t.slice(0, 24)}" found where a number was expected (offset ${o})`);
    this.o = e + 1; return x;
  }
  int() { const x = this.num(); if (x !== null && !Number.isInteger(x)) throw new Error(`${x} found where an integer was expected (offset ${this.o})`); return x; }
  byte() { return this.int(); }
  short() { return this.int(); }
  dbl() { return this.num(); }
  ptr() { return this.int() || 0; }
  type() { return this.int(); }
  count() { return this.int(); }
  nfields() { return this.int(); }
  raw() { if (this.o >= this.s.length) throw new Error('unexpected end of the transmit data'); return this.s[this.o++]; }
  ch() { const c = this.raw(); if (c !== '\\') return c; const d = this.raw(); return d === '0' ? '\0' : d === 'n' ? '\n' : d === 'r' ? '\r' : d === '9' ? '         ' : d; }
  chars(n) { let s = ''; while (s.length < n) s += this.ch(); return s; }
  sstr() { const n = this.int(); return this.s.slice(this.o, this.o += n); }
  logical() { const c = this.raw(); if (c !== 'T' && c !== 'F') throw new Error(`"${c}" found where a logical was expected (offset ${this.o - 1})`); return c === 'T'; }
  pos() { return this.int(); }
  d() { return this.int(); }
  n() { return this.int(); }
  u() { return this.int(); }
  vec() { const x = this.num(); if (x === null) return null; return [x, this.num(), this.num()]; }
  uf() { return this.int(); }
  done() { let o = this.o; while (this.s.charCodeAt(o) === 32) o++; return o >= this.s.length; }
}

function readField(r, code, n) {
  switch (code) {
    case 'p': if (n === 1) return r.ptr(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.ptr(); return a; }
    case 'f': if (n === 1) return r.dbl(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.dbl(); return a; }
    case 'c': return n === 1 ? r.ch() : r.chars(n);
    case 'd': case 't': if (n === 1) return r.d(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.d(); return a; }
    case 'n': case 'w': if (n === 1) return r.n(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.n(); return a; }
    case 'u': if (n === 1) return r.u(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.u(); return a; }
    case 'l': if (n === 1) return r.logical(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.logical(); return a; }
    case 'v': case 'h': if (n === 1) return r.vec(); { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = r.vec(); return a; }
    case 'i': { const a = []; for (let i = 0; i < 2 * n; i++) a.push(r.dbl()); return a; }
    case 'b': { const a = []; for (let i = 0; i < 6 * n; i++) a.push(r.dbl()); return a; }
    default: throw new Error(`field type "${code}" is not part of the XT schema language`);
  }
}
// One field description of an embedded schema: name, pointer class, element count, type (absent for pointers), transmit flag (variable-length fields only)
function fieldSpec(r, text) {
  const name = r.sstr(), cls = text ? r.int() : r.short(), n = text ? r.int() : r.pos(), code = cls ? 'p' : r.sstr()[0], xmt = n === 1 ? r.logical() : true;
  if (!/^[\w ]{1,64}$/.test(name) || !(n >= 0 && n < 1e6)) throw new Error(`implausible embedded field description "${name.slice(0, 20)}"`);
  return { name, code, n, cls, xmt };
}
// The field list of a node type in the file's schema: the base list edited by the Copy / Delete / Insert / Append sequence
// that precedes the first node of that type, or a complete description for a type the base schema does not have.
function embeddedSchema(r, text, type) {
  const base = SCHEMA.get(type), nf = r.nfields();
  if (nf === 255) { if (!base) throw new Error(`node type ${type} is declared identical to the base schema, which has no such type in this reader`); return base; }
  if (!base) { const name = r.sstr(), desc = r.sstr(), fields = []; if (!/^[\w ]{1,64}$/.test(name)) throw new Error(`node type ${type} is not in the base schema and is not described by the file`); for (let i = 0; i < nf; i++) fields.push(fieldSpec(r, text)); return { type, name: name.toUpperCase(), desc, fields, extension: true }; }
  const fields = []; let bi = 0;
  for (let g = 0; g < 4096; g++) {
    const op = r.raw ? r.raw() : r.ch();
    if (op === 'Z') break;
    if (op === 'C') { if (bi >= base.fields.length) throw new Error(`embedded schema of ${base.name} copies more fields than the base schema holds`); fields.push(base.fields[bi++]); }
    else if (op === 'D') bi++;
    else if (op === 'I' || op === 'A') fields.push(fieldSpec(r, text));
    else throw new Error(`"${op}" is not an embedded-schema instruction (node type ${type})`);
  }
  if (fields.length !== nf) throw new Error(`embedded schema of ${base.name} lists ${fields.length} fields where ${nf} were announced`);
  return { type, name: base.name, fields };
}

// ---------------------------------------------------------------- reader
const latin = (b, o, n) => { let s = ''; for (let i = o, e = min(b.length, o + n); i < e; i++) s += String.fromCharCode(b[i]); return s; };
function findSeq(b, seq, from, to) { outer: for (let i = from; i + seq.length <= to; i++) { for (let k = 0; k < seq.length; k++) if (b[i + k] !== seq[k]) continue outer; return i; } return -1; }
// The key = value pairs of the textual header that a frustrum writes in front of the transmit data
function headerKeys(head) {
  const out = {}, h = head.replace(/\r?\n/g, '');
  for (const m of h.matchAll(/(?:^|;)(?:\*\*PART\d;)?([A-Z_0-9]{1,80})=((?:\^.|[^;^])*)/g)) out[m[1].toLowerCase()] = m[2].replace(/\^(.)/g, (q, c) => c === 'n' ? '\n' : c === '_' ? ' ' : c).trim();
  return out;
}

// readXT(bytes) → {encoding, partition, modellerVersion, schema, baseSchema, userFieldSize, header, nodes: Map(index → node),
// root, counts: {NAME: n}, extensions: [names of node types the file itself describes], warnings}. A node is
// {index, type, name, …fields by schema name}; pointer fields hold node indices (0 = null).
export function readXT(bytes) {
  const b = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes), head = latin(b, 0, 16384), warnings = [];
  let start = 0, header = {};
  if (head.startsWith('**')) { const e = head.indexOf('**END_OF_HEADER'); if (e < 0) throw new Error('Parasolid: the text header has no END_OF_HEADER record'); let nl = head.indexOf('\n', e); if (nl < 0) throw new Error('Parasolid: nothing follows the header'); start = nl + 1; header = headerKeys(head.slice(0, e)); }
  while (start < b.length && (b[start] === 10 || b[start] === 13)) start++;
  let r, text, modeller, schema;
  if (b[start] === 0x50 && b[start + 1] === 0x53 && b[start + 2] === 0 && b[start + 3] === 0) {
    text = false; r = new Bin(b, start + 4); modeller = r.chars(r.short()); schema = r.chars(r.int());
  } else if (b[start] === 0x54) {
    text = true; let s = ''; const CH = 8192; for (let i = start + 1; i < b.length; i += CH) s += String.fromCharCode.apply(null, b.subarray(i, min(b.length, i + CH))); s = s.replace(/[\r\n]+/g, '');
    r = new Txt(s, 0); modeller = r.sstr(); schema = r.sstr();
  } else if (b[start] === 0x42) throw new Error('Parasolid: bare (machine-dependent) binary transmit files are not read; save as neutral binary or text');
  else throw new Error('Parasolid: no transmit flag sequence ("T" or "PS") follows the header');
  const mv = /TRANSMIT FILE\s*(\(([a-z ]+)\))?\s*created by modeller version (\d+)/.exec(modeller); if (!mv) throw new Error('Parasolid: the modeller identification string is missing');
  const sm = /^SCH_(\d+)_(\d+)(?:_(\d+))?$/.exec(schema); if (!sm) throw new Error(`Parasolid: unreadable schema name "${schema.slice(0, 40)}"`);
  const embedded = sm[3] !== undefined; if (embedded && sm[3] !== '13006') warnings.push(`the file's base schema is ${sm[3]}, not 13006: fields copied from the base schema may be misread`);
  if (!embedded && +sm[2] < 9000) warnings.push(`schema ${sm[2]} predates the layout this reader knows (12006 / 13006): node fields may be misread`);
  if (embedded) r.short();
  const usfld = text ? r.int() : r.int(); if (!(usfld >= 0 && usfld <= 16)) throw new Error(`Parasolid: implausible user-field size ${usfld}`);
  const nodes = new Map(), seen = new Map(), counts = {}, extensions = [];
  for (let guard = 0; guard < 5e7; guard++) {
    const type = r.type();
    if (type === 1) { r.ptr(); break; }
    let sc = seen.get(type);
    if (!sc) { if (embedded) sc = embeddedSchema(r, text, type); else sc = SCHEMA.get(type); if (!sc) throw new Error(`Parasolid: node type ${type} is not in schema ${sm[2]} as known to this reader`); seen.set(type, sc); if (sc.extension) extensions.push(sc.name); sc.variable = sc.fields.some(f => f.n === 1); }
    const len = sc.variable ? r.count() : 0, index = r.ptr();
    if (!(index > 0) || len < 0 || len > 1e8) throw new Error(`Parasolid: corrupt node record (type ${type}, index ${index}, length ${len})`);
    const node = { index, type, name: sc.name };
    for (const f of sc.fields) { if (f.xmt === false) continue; node[f.name] = f.n === 1 ? (len === 0 ? (f.code === 'c' ? '' : []) : f.code === 'c' ? r.chars(len) : len === 1 ? [readField(r, f.code, 1)] : readField(r, f.code, len)) : readField(r, f.code, f.n || 1); }
    for (let i = 0; i < usfld; i++) r.uf();
    nodes.set(index, node); counts[sc.name] = (counts[sc.name] || 0) + 1;
  }
  return { encoding: text ? 'text' : 'binary', partition: /partition/.test(mv[2] || ''), guise: mv[2] || 'part', modellerVersion: mv[3], schema, schemaNumber: +sm[2], baseSchema: embedded ? sm[3] : null, userFieldSize: usfld, header, nodes, root: nodes.get(1) || null, counts, extensions, warnings, trailing: !r.done() };
}

// ---------------------------------------------------------------- small vector and NURBS kit (sampling of curves that have no STEP form)
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = a => hypot(a[0], a[1], a[2]), dist = (a, b) => hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]), unit = a => { const l = len(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const isVec = v => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);
// de Boor evaluation of a B-spline with flat knot vector U, control points P (flat, `dim` numbers each) at parameter u
function deBoor(deg, U, P, dim, n, u) {
  let k = deg; const hi = n - 1; if (u >= U[hi + 1]) k = hi; else if (u > U[deg]) { let lo = deg, up = hi + 1; while (up - lo > 1) { const mid = (lo + up) >> 1; if (U[mid] <= u) lo = mid; else up = mid; } k = lo; }
  const d = new Array((deg + 1) * dim); for (let j = 0; j <= deg; j++) for (let c = 0; c < dim; c++) d[j * dim + c] = P[(k - deg + j) * dim + c];
  for (let r = 1; r <= deg; r++) for (let j = deg; j >= r; j--) { const i = k - deg + j, den = U[i + deg - r + 1] - U[i], a = den > 0 ? (u - U[i]) / den : 0; for (let c = 0; c < dim; c++) d[j * dim + c] = (1 - a) * d[(j - 1) * dim + c] + a * d[j * dim + c]; }
  return d.slice(deg * dim, deg * dim + dim);
}
const flatKnots = (knots, mult) => { const U = []; for (let i = 0; i < knots.length; i++) for (let m = 0; m < mult[i]; m++) U.push(knots[i]); return U; };
// A knot vector / control-point row in the non-periodic form STEP wants. Periodic data (the knots span one period and the
// wrapped control points are not stored) is unwrapped: knots are continued by the period and control points repeated.
function unwrap(deg, n, knots, mult, periodic) {
  const total = mult.reduce((a, c) => a + c, 0);
  if (total === n + deg + 1) return { knots, mult, n, map: null };
  if (!periodic) throw new Error(`B-spline data inconsistent: ${n} control points of degree ${deg} with ${total} knots`);
  const m = knots.length - 1, mu = mult[m], T = knots[m] - knots[0], per = total - mu; if (per !== n || mult[0] !== mu || !(T > 0) || mu > deg) throw new Error(`periodic B-spline data not understood: ${n} control points, degree ${deg}, ${total} knots`);
  const U = flatKnots(knots, mult), extra = deg + 1 - mu, F = [];
  for (let i = extra; i >= 1; i--) F.push(U[per - i] - T); for (const k of U) F.push(k); for (let i = 0; i < extra; i++) F.push(U[mu + i] + T);
  const nk = [], nm = []; for (const k of F) { if (nk.length && abs(k - nk[nk.length - 1]) <= 1e-14 * max(1, abs(k))) nm[nm.length - 1]++; else { nk.push(k); nm.push(1); } }
  const nn = F.length - deg - 1, map = []; for (let i = 0; i < nn; i++) map.push((((i - extra) % n) + n) % n);
  return { knots: nk, mult: nm, n: nn, map };
}

// Boehm knot insertion: parameter t inserted once into the flat knot vector U / control points P (in place)
function insertKnot(deg, U, P, dim, t) {
  let k = deg; while (k + 1 < U.length - deg - 1 && U[k + 1] <= t) k++;
  const Q = new Array(P.length + dim);
  for (let i = 0; i <= P.length / dim; i++) for (let c = 0; c < dim; c++) {
    if (i <= k - deg) Q[i * dim + c] = P[i * dim + c]; else if (i > k) Q[i * dim + c] = P[(i - 1) * dim + c];
    else { const den = U[i + deg] - U[i], a = den > 0 ? (t - U[i]) / den : 0; Q[i * dim + c] = (1 - a) * P[(i - 1) * dim + c] + a * P[i * dim + c]; }
  }
  U.splice(k + 1, 0, t); P.length = 0; for (const x of Q) P.push(x);
}
// The exact piece [t0, t1] of a B-spline as a clamped B-spline; t1 may run past the end of a periodic curve, whose knots
// and control points are then continued by one period first. Returns null when the data does not allow it.
function subCurve(q, t0, t1) {
  const p = q.degree, dim = q.dim, a = q.U[p], b = q.U[q.n], T = b - a, tol = 1e-9 * max(1, abs(a), abs(b), T); let U = q.U.slice(), P = q.P.slice(0, q.n * dim), n = q.n;
  if (!(t1 > t0 + tol) || t1 - t0 > T + tol) return null;
  if (q.closed) { const sh = floor((t0 - a) / T + 1e-9) * T; t0 -= sh; t1 -= sh; }
  if (t0 < a - tol || t0 > b + tol) return null; t0 = max(a, t0);
  if (t1 > b + tol) {
    const m = n - p; if (!q.closed || m < 1) return null;
    for (let i = 0; i + m <= n + p; i++) if (abs(U[i + m] - U[i] - T) > 1e-7 * T) return null;
    for (let i = 0; i < p * dim; i++) if (abs(P[m * dim + i] - P[i]) > 1e-7 * max(1, abs(P[i]))) return null;
    for (let j = n + p + 1; j <= n + p + m; j++) U.push(U[j - m] + T); for (let j = n; j < n + m; j++) for (let c = 0; c < dim; c++) P.push(P[(j - m) * dim + c]); n += m;
  } else t1 = min(b, t1);
  for (const t of [t0, t1]) { let have = 0; for (const u of U) if (abs(u - t) <= tol) have++; const at = U.find(u => abs(u - t) <= tol); for (let i = have; i < p; i++) insertKnot(p, U, P, dim, at !== undefined ? at : t); }
  const i0 = U.findIndex(u => abs(u - t0) <= tol), i1 = U.length - 1 - U.slice().reverse().findIndex(u => abs(u - t1) <= tol); if (i0 < 0 || i1 >= U.length || i1 <= i0) return null;
  // knots from the first t0 to the last t1, each end raised to multiplicity p + 1
  let m0 = 0, m1 = 0; for (const u of U) { if (abs(u - t0) <= tol) m0++; if (abs(u - t1) <= tol) m1++; }
  const inner = U.slice(i0 + m0, i1 - m1 + 1), flat = [...new Array(p + 1).fill(t0), ...inner, ...new Array(p + 1).fill(t1)], cnt = flat.length - p - 1, first = i0 + m0 - 1 - p + (p + 1 - m0) - (p + 1 - m0);
  const start = i0 + m0 - p - 1; if (start < 0 || (start + cnt) * dim > P.length || cnt < p + 1) return null;
  const knots = [], mult = []; for (const k of flat) { if (knots.length && abs(k - knots[knots.length - 1]) <= tol) mult[mult.length - 1]++; else { knots.push(k); mult.push(1); } }
  return { degree: p, dim, rational: q.rational, n: cnt, P: P.slice(start * dim, (start + cnt) * dim), knots, mult, U: flat, closed: false };
}
// Parameter of the point of an evaluable curve nearest to p (coarse scan, then interval refinement)
function nearestParam(e, p, lo, hi) {
  const cuts = [lo]; for (const k of e.breaks || []) if (k > lo && k < hi) cuts.push(k); cuts.push(hi); let best = lo, bd = Infinity, step = 0;
  for (let i = 0; i + 1 < cuts.length; i++) for (let j = 0; j <= 8; j++) { const t = cuts[i] + (cuts[i + 1] - cuts[i]) * j / 8, d = dist(e.f(t), p); if (d < bd) { bd = d; best = t; step = (cuts[i + 1] - cuts[i]) / 8; } }
  let a = max(lo, best - step), b = min(hi, best + step);
  for (let g = 0; g < 60 && b - a > 1e-14 * max(1, abs(b)); g++) { const m1 = a + (b - a) * 0.381966, m2 = a + (b - a) * 0.618034; if (dist(e.f(m1), p) < dist(e.f(m2), p)) b = m2; else a = m1; }
  return (a + b) / 2;
}

// ---------------------------------------------------------------- geometry access
const SURFACES = new Set(['PLANE', 'CYLINDER', 'CONE', 'SPHERE', 'TORUS', 'BLENDED_EDGE', 'BLEND_BOUND', 'OFFSET_SURF', 'SWEPT_SURF', 'SPUN_SURF', 'PE_SURF', 'B_SURFACE']);
function nurbsCurve(M, bc) {
  const nu = M.get(bc.nurbs), V = nu && M.get(nu.bspline_vertices), K = nu && M.get(nu.knots), Q = nu && M.get(nu.knot_mult); if (!nu || !V || !K || !Q) throw new Error('B-curve without NURBS data');
  const nk = nu.n_knots > 0 ? nu.n_knots : K.knots.length, w = unwrap(nu.degree, nu.n_vertices, K.knots.slice(0, nk), Q.mult.slice(0, nk), nu.periodic), dim = nu.vertex_dim, src = V.vertices; if (src.length < nu.n_vertices * dim) throw new Error('B-curve control points are incomplete');
  let P = src; if (w.map) { P = []; for (const i of w.map) for (let c = 0; c < dim; c++) P.push(src[i * dim + c]); }
  return { degree: nu.degree, dim, rational: !!nu.rational, n: w.n, P, knots: w.knots, mult: w.mult, U: flatKnots(w.knots, w.mult), closed: !!nu.closed || !!nu.periodic };
}
function nurbsSurface(M, bs) {
  const nu = M.get(bs.nurbs), V = nu && M.get(nu.bspline_vertices); if (!nu || !V) throw new Error('B-surface without NURBS data');
  const ku = M.get(nu.u_knots), kv = M.get(nu.v_knots), mu = M.get(nu.u_knot_mult), mv = M.get(nu.v_knot_mult); if (!ku || !kv || !mu || !mv) throw new Error('B-surface without knot data');
  const nku = nu.n_u_knots > 0 ? nu.n_u_knots : ku.knots.length, nkv = nu.n_v_knots > 0 ? nu.n_v_knots : kv.knots.length;
  const wu = unwrap(nu.u_degree, nu.n_u_vertices, ku.knots.slice(0, nku), mu.mult.slice(0, nku), nu.u_periodic), wv = unwrap(nu.v_degree, nu.n_v_vertices, kv.knots.slice(0, nkv), mv.mult.slice(0, nkv), nu.v_periodic), dim = nu.vertex_dim, src = V.vertices, nv0 = nu.n_v_vertices;
  if (src.length < nu.n_u_vertices * nv0 * dim) throw new Error('B-surface control points are incomplete');
  // control points are stored with the v index varying fastest: vertex (i, j) is at i * n_v + j
  const P = new Array(wu.n * wv.n * dim);
  for (let i = 0; i < wu.n; i++) for (let j = 0; j < wv.n; j++) { const si = wu.map ? wu.map[i] : i, sj = wv.map ? wv.map[j] : j; for (let c = 0; c < dim; c++) P[(i * wv.n + j) * dim + c] = src[(si * nv0 + sj) * dim + c]; }
  return { du: nu.u_degree, dv: nu.v_degree, nu: wu.n, nv: wv.n, dim, rational: !!nu.rational, P, ku: wu.knots, mu: wu.mult, kv: wv.knots, mv: wv.mult, Uu: flatKnots(wu.knots, wu.mult), Uv: flatKnots(wv.knots, wv.mult), uClosed: !!nu.u_closed || !!nu.u_periodic, vClosed: !!nu.v_closed || !!nu.v_periodic };
}
const dehom = (p, dim, rational) => rational ? [p[0] / p[dim - 1], p[1] / p[dim - 1], dim > 3 ? p[2] / p[dim - 1] : 0] : [p[0], p[1], dim > 2 ? p[2] : 0];
// Evaluator of a curve node: {f(t) → point, range: [t0, t1] | null, closed}, or null when the curve has no closed form here
function curveEval(M, c, memo) {
  if (!c) return null; if (memo.has(c.index)) return memo.get(c.index); let e = null;
  try {
    if (c.name === 'LINE' && isVec(c.pvec) && isVec(c.direction)) e = { f: t => add(c.pvec, mul(c.direction, t)), range: null };
    else if ((c.name === 'CIRCLE' || c.name === 'ELLIPSE') && isVec(c.centre) && isVec(c.normal) && isVec(c.x_axis)) { const Y = cross(c.normal, c.x_axis), a = c.name === 'CIRCLE' ? c.radius : c.major_radius, b = c.name === 'CIRCLE' ? c.radius : c.minor_radius; e = { f: t => add(c.centre, add(mul(c.x_axis, a * cos(t)), mul(Y, b * sin(t)))), range: [0, 2 * PI], closed: true }; }
    else if (c.name === 'B_CURVE') { const q = nurbsCurve(M, c); e = { f: t => dehom(deBoor(q.degree, q.U, q.P, q.dim, q.n, t), q.dim, q.rational), range: [q.U[q.degree], q.U[q.n]], closed: q.closed, breaks: q.knots, degree: q.degree }; }
    else if (c.name === 'TRIMMED_CURVE') { const b = curveEval(M, M.get(c.basis_curve), memo); if (b) e = { ...b, range: [min(c.parm_1, c.parm_2), max(c.parm_1, c.parm_2)], closed: false }; }
    else if (c.name === 'SP_CURVE') { const s = surfEval(M, M.get(c.surface), memo), bc = M.get(c.b_curve); if (s && bc) { const q = nurbsCurve(M, bc); e = { f: t => { const h = deBoor(q.degree, q.U, q.P, q.dim, q.n, t); return q.rational ? s.f(h[0] / h[q.dim - 1], h[1] / h[q.dim - 1]) : s.f(h[0], h[1]); }, range: [q.U[q.degree], q.U[q.n]], closed: q.closed, breaks: q.knots, degree: max(3, q.degree) }; } }
  } catch (err) { e = null; }
  memo.set(c.index, e); return e;
}
// Evaluator of a surface node: {f(u, v) → point} or null
function surfEval(M, s, memo) {
  if (!s) return null; const key = 's' + s.index; if (memo.has(key)) return memo.get(key); let e = null;
  try {
    const ax = s.axis, X = s.x_axis, Y = isVec(ax) && isVec(X) ? cross(ax, X) : null;
    if (s.name === 'PLANE') { const Yp = cross(s.normal, s.x_axis); e = { f: (u, v) => add(s.pvec, add(mul(s.x_axis, u), mul(Yp, v))) }; }
    else if (s.name === 'CYLINDER') e = { f: (u, v) => add(s.pvec, add(add(mul(X, s.radius * cos(u)), mul(Y, s.radius * sin(u))), mul(ax, v))) };
    else if (s.name === 'CONE') { const t = s.sin_half_angle / s.cos_half_angle; e = { f: (u, v) => { const r = s.radius + v * t; return add(add(s.pvec, mul(ax, v)), add(mul(X, r * cos(u)), mul(Y, r * sin(u)))); } }; }
    else if (s.name === 'SPHERE') e = { f: (u, v) => add(s.centre, add(mul(add(mul(X, cos(u)), mul(Y, sin(u))), s.radius * cos(v)), mul(ax, s.radius * sin(v)))) };
    else if (s.name === 'TORUS') e = { f: (u, v) => add(s.centre, add(mul(add(mul(X, cos(u)), mul(Y, sin(u))), s.major_radius + s.minor_radius * cos(v)), mul(ax, s.minor_radius * sin(v)))) };
    else if (s.name === 'B_SURFACE') { const q = nurbsSurface(M, s), row = new Array(q.nu * q.dim); e = { f: (u, v) => { for (let i = 0; i < q.nu; i++) { const p = deBoor(q.dv, q.Uv, q.P.slice(i * q.nv * q.dim, (i + 1) * q.nv * q.dim), q.dim, q.nv, v); for (let c = 0; c < q.dim; c++) row[i * q.dim + c] = p[c]; } return dehom(deBoor(q.du, q.Uu, row, q.dim, q.nu, u), q.dim, q.rational); } }; }
    else if (s.name === 'SWEPT_SURF') { const c = curveEval(M, M.get(s.section), memo); if (c) e = { f: (u, v) => add(c.f(u), mul(s.sweep, v)) }; }
    else if (s.name === 'SPUN_SURF') { const c = curveEval(M, M.get(s.profile), memo), A = s.axis, B = s.base; if (c) e = { f: (u, v) => { const p = c.f(u), z = add(B, mul(A, dot(sub(p, B), A))), d = sub(p, z); return add(z, add(mul(d, cos(v)), mul(cross(A, d), sin(v)))); } }; }
    else if (s.name === 'OFFSET_SURF') { const b = M.get(s.surface), be = surfEval(M, b, memo); if (be && b) { const sg = b.sense === '-' ? -1 : 1; e = { f: (u, v) => { const h = 1e-6 * max(1, abs(u), abs(v)), p = be.f(u, v), n = unit(cross(sub(be.f(u + h, v), be.f(u - h, v)), sub(be.f(u, v + h), be.f(u, v - h)))); return add(p, mul(n, s.offset * sg)); } }; } }
  } catch (err) { e = null; }
  memo.set(key, e); return e;
}

// ---------------------------------------------------------------- polylines → interpolating B-spline (approximated edges)
// Nearest position of point p on polyline Q: {i, t, d} with the point at Q[i] + t (Q[i+1] − Q[i])
function locate(Q, p) {
  let best = { i: 0, t: 0, d: Infinity };
  for (let i = 0; i + 1 < Q.length; i++) { const a = Q[i], s = sub(Q[i + 1], a), l2 = dot(s, s), t = l2 > 0 ? max(0, min(1, dot(sub(p, a), s) / l2)) : 0, d = dist(p, add(a, mul(s, t))); if (d < best.d) best = { i, t, d }; }
  return best;
}
// The part of polyline Q (closed when Q's ends coincide) that runs from point a to point b in the direction of increasing index
function cutPolyline(Q, closed, a, b) {
  if (!a && !b) return Q.slice();
  if (!a || !b) { const p = a || b, l = locate(Q, p); return a ? [a, ...Q.slice(l.i + 1)] : [...Q.slice(0, l.i + 1), b]; }
  const la = locate(Q, a), lb = locate(Q, b), pa = la.i + la.t, pb = lb.i + lb.t, same = dist(a, b) <= 1e-12 * max(1, len(a));
  if (!same && pb > pa) return [a, ...Q.slice(la.i + 1, lb.i + 1), b];
  if (closed) return [a, ...Q.slice(la.i + 1, Q.length - 1), ...Q.slice(0, lb.i + 1), b];
  if (same) return [a, ...Q.slice(1, Q.length - 1), b];
  return null;
}
function tidy(Q) {
  let L = 0; for (let i = 0; i + 1 < Q.length; i++) L += dist(Q[i], Q[i + 1]); const eps = 1e-7 * L, out = [Q[0]];
  for (let i = 1; i < Q.length - 1; i++) if (dist(Q[i], out[out.length - 1]) > eps && dist(Q[i], Q[Q.length - 1]) > eps) out.push(Q[i]);
  out.push(Q[Q.length - 1]); return L > 0 && out.length >= 2 ? out : null;
}
// C1 cubic through the points (chord-length parameters unless `t` is given, Catmull-Rom tangents) as B-spline poles with
// triple interior knots
function interpolate(Q, t, closedFlag) {
  const m = Q.length - 1; if (!t) { t = [0]; for (let i = 0; i < m; i++) t.push(t[i] + dist(Q[i], Q[i + 1])); }
  if (m === 1) return { degree: 1, P: [Q[0], Q[1]], knots: [0, t[1]], mult: [2, 2] };
  const closed = closedFlag !== undefined ? closedFlag : dist(Q[0], Q[m]) <= 1e-9 * t[m], T = [], slope = (a, b) => mul(sub(Q[b], Q[a]), 1 / ((t[b] - t[a]) || 1));
  for (let i = 0; i <= m; i++) {
    if (i === 0) { if (closed) T.push(mul(sub(Q[1], Q[m - 1]), 1 / (t[1] + t[m] - t[m - 1]))); else { const w = (t[1] - t[0]) / (t[2] - t[0]); T.push(add(mul(slope(0, 1), 1 + w), mul(slope(1, 2), -w))); } }
    else if (i === m) { if (closed) T.push(T[0]); else { const w = (t[m] - t[m - 1]) / (t[m] - t[m - 2]); T.push(add(mul(slope(m - 1, m), 1 + w), mul(slope(m - 2, m - 1), -w))); } }
    else T.push(slope(i - 1, i + 1));
  }
  const P = [Q[0]], knots = [t[0]], mult = [4];
  for (let i = 0; i < m; i++) { const h = (t[i + 1] - t[i]) / 3; P.push(add(Q[i], mul(T[i], h)), sub(Q[i + 1], mul(T[i + 1], h)), Q[i + 1]); knots.push(t[i + 1]); mult.push(i === m - 1 ? 4 : 3); }
  return { degree: 3, P, knots, mult };
}
// Nearest-point map onto an analytic surface node (used for the contact curves of rolling-ball blends), or null
function projector(s) {
  if (!s) return null; const A = s.axis, C = s.name === 'SPHERE' || s.name === 'TORUS' ? s.centre : s.pvec;
  const radial = p => { const d = sub(p, C), h = dot(d, A), e = sub(d, mul(A, h)); return { h, rho: len(e), e: unit(e) }; };
  if (s.name === 'PLANE') return p => sub(p, mul(s.normal, dot(sub(p, C), s.normal)));
  if (s.name === 'CYLINDER') return p => { const q = radial(p); return add(C, add(mul(A, q.h), mul(q.e, abs(s.radius)))); };
  if (s.name === 'SPHERE') return p => add(C, mul(unit(sub(p, C)), abs(s.radius)));
  if (s.name === 'CONE') { const t = s.sin_half_angle / s.cos_half_angle; return p => { const q = radial(p), v = ((q.rho - s.radius) * t + q.h) / (t * t + 1); return add(C, add(mul(A, v), mul(q.e, s.radius + v * t))); }; }
  if (s.name === 'TORUS') return p => { const q = radial(p), c = add(C, mul(q.e, s.major_radius)); return add(c, mul(unit(sub(p, c)), abs(s.minor_radius))); };
  return null;
}
// Points along an evaluable curve between two parameters, dense enough for the interpolant to follow it
function sampleCurve(e, t0, t1) {
  const lo = min(t0, t1), hi = max(t0, t1), cuts = [lo]; for (const k of e.breaks || []) if (k > lo + 1e-12 && k < hi - 1e-12) cuts.push(k); cuts.push(hi);
  const per = e.breaks ? max(2, min(8, 2 * (e.degree || 3))) : 0, Q = [];
  if (per) { for (let i = 0; i + 1 < cuts.length; i++) for (let j = 0; j < per; j++) Q.push(e.f(cuts[i] + (cuts[i + 1] - cuts[i]) * j / per)); Q.push(e.f(hi)); }
  else { const n = max(8, Math.ceil(48 * (hi - lo) / (2 * PI))); for (let j = 0; j <= n; j++) Q.push(e.f(lo + (hi - lo) * j / n)); }
  if (Q.length > 4000) return null; if (t1 < t0) Q.reverse(); return Q.every(isVec) ? Q : null;
}

// ---------------------------------------------------------------- attributes
// Attributes attached to a node: [{id: 'SDL/TYSA_COLOUR', values: [[…], …]}]
function attributesOf(M, node) {
  const out = [];
  for (let a = M.get(node.attributes_groups), g = 0; a && g < 1000; a = M.get(a.next), g++) {
    if (a.name !== 'ATTRIBUTE') continue; const def = M.get(a.definition), id = def && M.get(def.identifier); if (!id) continue;
    out.push({ id: id.string, values: (a.fields || []).map(i => { const v = M.get(i); return !v ? null : v.name === 'UNICODE_VALUES' ? String.fromCharCode(...v.values.filter(c => c > 0)) : v.values; }) });
  }
  return out;
}
const attr = (M, node, id) => { const a = attributesOf(M, node).find(x => x.id === id); return a ? a.values : null; };

// ---------------------------------------------------------------- XT → STEP
const real = x => { if (!Number.isFinite(x)) throw new Error('non-finite coordinate'); if (x === 0) return '0.'; const s = String(x), i = s.indexOf('e'); if (i < 0) return s.includes('.') ? s : s + '.'; const m = s.slice(0, i); return (m.includes('.') ? m : m + '.') + 'E' + (+s.slice(i + 1)); };
const str = s => "'" + String(s).replace(/[^\x20-\x7e]/g, c => c.charCodeAt(0) > 0xff ? '\\X2\\' + c.charCodeAt(0).toString(16).toUpperCase().padStart(4, '0') + '\\X0\\' : '_').replace(/\\(?!X2\\|X0\\)/g, '\\\\').replace(/'/g, "''") + "'";
const B = v => v ? '.T.' : '.F.';
const BODY_TYPES = { 1: 'solid', 2: 'wire', 3: 'sheet', 6: 'general' };

// xtToStep(model, opts) → {step: AP214 text | null, bodies: [{name, type, faces, translated, closed, colour}], stats, notes, lost}.
// `model` is the result of readXT. Every face whose surface and edge curves have a STEP form is written as an ADVANCED_FACE;
// a body whose faces all translate becomes a MANIFOLD_SOLID_BREP (BREP_WITH_VOIDS for cavities), otherwise its faces go into
// a SHELL_BASED_SURFACE_MODEL. Lengths are written in metres, the unit of every XT file. opts: name (product name stem).
export function xtToStep(model, opts = {}) {
  const M = model.nodes, L = [], E = s => '#' + L.push(s), memo = new Map();
  const stats = { bodies: 0, solids: 0, sheets: 0, wireBodies: 0, faces: 0, translated: 0, skipped: {}, edges: 0, exactEdges: 0, approximatedEdges: 0, approximated: {}, edgeCurves: {}, surfaces: {}, instances: 0 }, bodies = [], notes = [], lost = [];
  const skip = (why, n = 1) => { stats.skipped[why] = (stats.skipped[why] || 0) + n; };
  let tf = null, tfKey = '';
  const P = p => { if (!isVec(p)) throw new Error('unset point'); if (!tf) return p; const R = tf.R; return [(R[0] * p[0] + R[1] * p[1] + R[2] * p[2] + tf.t[0]) * tf.s, (R[3] * p[0] + R[4] * p[1] + R[5] * p[2] + tf.t[1]) * tf.s, (R[6] * p[0] + R[7] * p[1] + R[8] * p[2] + tf.t[2]) * tf.s]; };
  const D = d => { if (!isVec(d)) throw new Error('unset direction'); if (!tf) return d; const R = tf.R; return [R[0] * d[0] + R[1] * d[1] + R[2] * d[2], R[3] * d[0] + R[4] * d[1] + R[5] * d[2], R[6] * d[0] + R[7] * d[1] + R[8] * d[2]]; };
  const S = r => tf ? r * tf.s : r;
  const pt = p => E(`CARTESIAN_POINT('',(${real(p[0])},${real(p[1])},${real(p[2])}))`);
  const dir = d => { const l = len(d); if (!(l > 1e-12)) throw new Error('zero direction'); return E(`DIRECTION('',(${real(d[0] / l)},${real(d[1] / l)},${real(d[2] / l)}))`); };
  const ax2 = (o, z, x) => { const zz = unit(D(z)); let xx = D(x); xx = sub(xx, mul(zz, dot(xx, zz))); return E(`AXIS2_PLACEMENT_3D('',${pt(P(o))},${dir(zz)},${dir(xx)})`); };
  const list = a => '(' + a.join(',') + ')';
  const bspCurve = (deg, pts, knots, mult, weights, closed) => {
    const ids = list(pts.map(pt)), k = `${list(mult)},${list(knots.map(real))},.UNSPECIFIED.`;
    return weights ? E(`(BOUNDED_CURVE()B_SPLINE_CURVE(${deg},${ids},.UNSPECIFIED.,${B(closed)},.F.)B_SPLINE_CURVE_WITH_KNOTS(${k})CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE(${list(weights.map(real))})REPRESENTATION_ITEM(''))`) : E(`B_SPLINE_CURVE_WITH_KNOTS('',${deg},${ids},.UNSPECIFIED.,${B(closed)},.F.,${k})`);
  };
  // exact curves: {id} or null
  const curve3 = c => {
    const key = 'c' + c.index + tfKey; if (memo.has(key)) return memo.get(key); let id = null;
    if (c.name === 'LINE') id = E(`LINE('',${pt(P(c.pvec))},${E(`VECTOR('',${dir(D(c.direction))},1.)`)})`);
    else if (c.name === 'CIRCLE') id = E(`CIRCLE('',${ax2(c.centre, c.normal, c.x_axis)},${real(S(c.radius))})`);
    else if (c.name === 'ELLIPSE') id = E(`ELLIPSE('',${ax2(c.centre, c.normal, c.x_axis)},${real(S(c.major_radius))},${real(S(c.minor_radius))})`);
    else if (c.name === 'B_CURVE') { const q = nurbsCurve(M, c); if (q.dim !== (q.rational ? 4 : 3)) throw new Error('B-curve is not three-dimensional'); const pts = [], w = []; for (let i = 0; i < q.n; i++) { const h = q.P.slice(i * q.dim, (i + 1) * q.dim); pts.push(P(dehom(h, q.dim, q.rational))); if (q.rational) w.push(h[3]); } id = bspCurve(q.degree, pts, q.knots, q.mult, q.rational ? w : null, q.closed); }
    memo.set(key, id); return id;
  };
  // surfaces: {id, flip} (flip: the STEP surface normal is opposite to the XT natural normal) or {why}
  const surface = s => {
    const key = 'S' + s.index + tfKey; if (memo.has(key)) return memo.get(key); let r;
    try {
      if (s.name === 'PLANE') r = { id: E(`PLANE('',${ax2(s.pvec, s.normal, s.x_axis)})`), flip: false };
      else if (s.name === 'CYLINDER') r = { id: E(`CYLINDRICAL_SURFACE('',${ax2(s.pvec, s.axis, s.x_axis)},${real(S(abs(s.radius)))})`), flip: false };
      else if (s.name === 'SPHERE') r = { id: E(`SPHERICAL_SURFACE('',${ax2(s.centre, s.axis, s.x_axis)},${real(S(abs(s.radius)))})`), flip: false };
      else if (s.name === 'TORUS') { const a = s.major_radius, b = abs(s.minor_radius), pl = ax2(s.centre, s.axis, s.x_axis); r = a > b ? { id: E(`TOROIDAL_SURFACE('',${pl},${real(S(a))},${real(S(b))})`), flip: false } : { id: E(`DEGENERATE_TOROIDAL_SURFACE('',${pl},${real(S(abs(a)))},${real(S(b))},${B(a > 0)})`), flip: false }; }
      else if (s.name === 'CONE') {
        // as transmitted the radius grows along +A: R(u, v) = P + vA + (X cos u + Y sin u)(r + v tan a), which is STEP's own form;
        // a negative half angle (radius shrinking along +A) is written with the axis reversed, which also reverses the normal
        const t = s.sin_half_angle / s.cos_half_angle, ang = atan2(abs(s.sin_half_angle), abs(s.cos_half_angle)); if (!(ang > 1e-9 && ang < PI / 2 - 1e-9)) throw new Error('degenerate cone angle');
        let o = s.pvec, rad = s.radius; if (rad < 0) { o = add(o, mul(s.axis, -rad / t)); rad = 0; }
        r = { id: E(`CONICAL_SURFACE('',${ax2(o, mul(s.axis, t > 0 ? 1 : -1), s.x_axis)},${real(S(rad))},${real(ang)})`), flip: t < 0 };
      } else if (s.name === 'B_SURFACE') {
        const q = nurbsSurface(M, s); if (q.dim !== (q.rational ? 4 : 3)) throw new Error('B-surface is not three-dimensional'); const rows = [], wr = [];
        for (let i = 0; i < q.nu; i++) { const row = [], w = []; for (let j = 0; j < q.nv; j++) { const h = q.P.slice((i * q.nv + j) * q.dim, (i * q.nv + j + 1) * q.dim); row.push(pt(P(dehom(h, q.dim, q.rational)))); if (q.rational) w.push(real(h[3])); } rows.push(list(row)); wr.push(list(w)); }
        const k = `${list(q.mu)},${list(q.mv)},${list(q.ku.map(real))},${list(q.kv.map(real))},.UNSPECIFIED.`, head = `${q.du},${q.dv},${list(rows)},.UNSPECIFIED.,${B(q.uClosed)},${B(q.vClosed)},.F.`;
        r = { id: q.rational ? E(`(BOUNDED_SURFACE()B_SPLINE_SURFACE(${head})B_SPLINE_SURFACE_WITH_KNOTS(${k})GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(${list(wr)})REPRESENTATION_ITEM('')SURFACE())`) : E(`B_SPLINE_SURFACE_WITH_KNOTS('',${head},${k})`), flip: false };
      } else if (s.name === 'SWEPT_SURF') { const c = M.get(s.section), id = c && curve3(c); r = id ? { id: E(`SURFACE_OF_LINEAR_EXTRUSION('',${id},${E(`VECTOR('',${dir(D(s.sweep))},1.)`)})`), flip: false } : { why: `swept surface on ${c ? c.name : 'a missing'} section curve` }; }
      else if (s.name === 'SPUN_SURF') { const c = M.get(s.profile), id = c && curve3(c); r = id ? { id: E(`SURFACE_OF_REVOLUTION('',${id},${E(`AXIS1_PLACEMENT('',${pt(P(s.base))},${dir(D(s.axis))})`)})`), flip: true } : { why: `spun surface on ${c ? c.name : 'a missing'} profile curve` }; }
      else if (s.name === 'OFFSET_SURF') { const b = M.get(s.surface), br = b ? surface(b) : { why: 'offset of a missing surface' }; r = br.id ? { id: E(`OFFSET_SURFACE('',${br.id},${real(S(s.offset) * ((b.sense === '-') !== br.flip ? -1 : 1))},.F.)`), flip: br.flip } : { why: `offset of: ${br.why}` }; }
      else r = { why: `${s.name} surface` };
    } catch (e) { r = { why: `${s.name} surface (${e.message})` }; }
    if (r.id) stats.surfaces[s.name] = (stats.surfaces[s.name] || 0) + 1;
    memo.set(key, r); return r;
  };
  // A rolling-ball blend has no STEP form: the part of it under one face is written as a bicubic B-spline surface through
  // circular sections taken along the spine (exact contact points on the two analytic supporting surfaces; each section a
  // chain of cubic arcs). An approximation, counted in stats.approximatedFaces.
  const blend = (fa, s) => {
    if (s.blend_type !== 'R') return { why: 'cliff-edge blend surface' };
    const sp = M.get(s.spine), s1 = M.get(s.surface[0]), s2 = M.get(s.surface[1]), p1 = projector(s1), p2 = projector(s2), r = abs(s.range[0]);
    if (!p1 || !p2) return { why: `rolling-ball blend between ${s1 ? s1.name : 'a missing surface'} and ${s2 ? s2.name : 'a missing surface'}` };
    if (!sp || !(r > 0) || abs(abs(s.range[1]) - r) > 1e-9 * r) return { why: 'rolling-ball blend with unequal or unset radii' };
    let Q = null; if (sp.name === 'INTERSECTION') { const ch = M.get(sp.chart); Q = ch && ch.hvec ? ch.hvec.filter(isVec) : null; } else { const ev = curveEval(M, sp, memo); if (ev && ev.range) Q = ev.breaks ? sampleCurve(ev, ev.range[0], ev.range[1]) : Array.from({ length: 145 }, (_, i) => ev.f(ev.range[0] + (ev.range[1] - ev.range[0]) * i / 144)); }
    Q = Q && tidy(Q); if (!Q || Q.length < 3) return { why: `rolling-ball blend on a ${sp.name} spine that cannot be sampled` };
    let closed = dist(Q[0], Q[Q.length - 1]) <= 1e-9 * max(1, len(Q[0])); const n = Q.length - 1, pos = [];
    for (let lp = M.get(fa.loop), g = 0; lp && g < 1e4; lp = M.get(lp.next), g++) { const f0 = M.get(lp.fin); for (let fin = f0, h = 0; fin && h < 1e4; h++) { const v = vpos(M.get(fin.vertex)); if (v) { const l = locate(Q, v); pos.push(l.i + l.t); } fin = M.get(fin.forward); if (fin === f0) break; } }
    if (pos.length) {
      pos.sort((a, b) => a - b); let lo = pos[0], hi = pos[pos.length - 1];
      if (closed) { let gap = lo + n - hi, at = -1; for (let i = 0; i + 1 < pos.length; i++) if (pos[i + 1] - pos[i] > gap) { gap = pos[i + 1] - pos[i]; at = i; } if (gap > 4) { if (at >= 0) { lo = pos[at + 1]; hi = pos[at] + n; } const a = floor(lo) - 1, b = Math.ceil(hi) + 1, R2 = []; for (let i = a; i <= b; i++) R2.push(Q[((i % n) + n) % n]); Q = R2; closed = false; } }
      else Q = Q.slice(max(0, floor(lo) - 1), min(n, Math.ceil(hi) + 1) + 1);
    }
    const [pa, pb] = sp.sense === '-' ? [p2, p1] : [p1, p2], secs = []; let amax = 0, prevY = null;
    for (const c of Q) { const a = pa(c), b = pb(c); if (!isVec(a) || !isVec(b) || abs(dist(a, c) - r) > 1e-4 * r || abs(dist(b, c) - r) > 1e-4 * r) return { why: 'rolling-ball blend whose spine is not at the blend radius from its supporting surfaces' }; const X = unit(sub(a, c)), w = sub(b, c), wy = sub(w, mul(X, dot(w, X))), ang = atan2(len(wy), dot(w, X)), Y = len(wy) > 1e-9 * r ? unit(wy) : prevY; if (!Y) return { why: 'rolling-ball blend with a degenerate first section' }; prevY = Y; amax = max(amax, ang); secs.push({ c, X, Y, ang }); }
    const m = max(1, Math.ceil(1.06 * amax / (PI / 4))), rows = secs.map(q => { const a0 = -0.03 * q.ang, th = 1.06 * q.ang / m, k = 4 / 3 * Math.tan(th / 4), at = a => add(q.c, add(mul(q.X, r * cos(a)), mul(q.Y, r * sin(a)))), tan = a => add(mul(q.X, -r * sin(a)), mul(q.Y, r * cos(a))), row = [at(a0)]; for (let j = 0; j < m; j++) { const b0 = a0 + j * th, b1 = b0 + th; row.push(add(at(b0), mul(tan(b0), k)), sub(at(b1), mul(tan(b1), k)), at(b1)); } return row; });
    const tU = [0]; for (let i = 0; i + 1 < Q.length; i++) tU.push(tU[i] + dist(Q[i], Q[i + 1])); const cols = []; let ku = null, mu = null;
    for (let j = 0; j <= 3 * m; j++) { const c = interpolate(rows.map(row => row[j]), tU, closed); cols.push(c.P); ku = c.knots; mu = c.mult; if (c.degree !== 3) return { why: 'rolling-ball blend with too few sections' }; }
    const net = cols[0].map((_, i) => list(cols.map(col => pt(P(col[i]))))), kv = [], mv = []; for (let j = 0; j <= m; j++) { kv.push(real(j / m)); mv.push(j === 0 || j === m ? 4 : 3); }
    stats.approximatedFaces = (stats.approximatedFaces || 0) + 1;
    return { id: E(`B_SPLINE_SURFACE_WITH_KNOTS('',3,3,${list(net)},.UNSPECIFIED.,${B(closed)},.F.,.F.,${list(mu)},${list(mv)},${list(ku.map(real))},${list(kv)},.UNSPECIFIED.)`), flip: false };
  };
  const vertex = v => { const key = 'v' + v.index + tfKey; if (memo.has(key)) return memo.get(key); const p = M.get(v.point); if (!p || !isVec(p.pvec)) throw new Error('vertex without a point'); const id = E(`VERTEX_POINT('',${pt(P(p.pvec))})`); memo.set(key, id); return id; };
  const vpos = v => { const p = v && M.get(v.point); return p && isVec(p.pvec) ? p.pvec : null; };
  const approx = (kind, e, Q, vs, ve) => {
    Q = Q && tidy(Q); if (!Q) return { why: `${kind} with no usable sample points` };
    const c = interpolate(Q), id = bspCurve(c.degree, c.P.map(P), c.knots, c.mult, null, false), a = vs ? vertex(vs) : E(`VERTEX_POINT('',${pt(P(Q[0]))})`), b = ve ? (ve === vs ? a : vertex(ve)) : a;
    stats.approximatedEdges++; stats.approximated[kind] = (stats.approximated[kind] || 0) + 1; return { id: E(`EDGE_CURVE('',${a},${b},${id},.T.)`) };
  };
  const edge = e => {
    const key = 'e' + e.index + tfKey; if (memo.has(key)) return memo.get(key); let r; stats.edges++;
    try {
      let vs = null, ve = null, fc = null; const f0 = M.get(e.fin);
      for (let f = f0, g = 0; f && g < 4096; g++) { const v = M.get(f.vertex); if (f.sense === '-') { if (v) vs = v; } else if (v) ve = v; if (!fc && f.curve) fc = f; f = M.get(f.other); if (f === f0) break; }
      if (!vs) vs = ve; if (!ve) ve = vs;
      let c = M.get(e.curve), lo = null, hi = null, viaFin = false; if (!c && fc) { c = M.get(fc.curve); viaFin = true; }
      if (!c) r = { why: 'edge without a curve' };
      else {
        for (let g = 0; c.name === 'TRIMMED_CURVE' && g < 8; g++) { if (lo === null) { lo = c.parm_1; hi = c.parm_2; } const b = M.get(c.basis_curve); if (!b) break; c = b; }
        let fwd = c.sense !== '-'; const a = vpos(vs), b = vpos(ve); stats.edgeCurves[c.name] = (stats.edgeCurves[c.name] || 0) + 1;
        if (viaFin) {
          // tolerant edge: the fin's trimmed SP-curve, sampled between its trimming parameters, runs in the direction of the fin
          const ev = curveEval(M, c, memo); let Q = ev && lo !== null ? sampleCurve(ev, lo, hi) : null; if (Q && fc.sense === '-') Q.reverse();
          if (Q && a && b) { Q[0] = a; Q[Q.length - 1] = b; }
          r = Q ? approx('tolerant edge (surface curve)', e, Q, vs, ve) : { why: `tolerant edge on a ${c.name} that cannot be evaluated` };
        } else if (c.name === 'LINE' || c.name === 'CIRCLE' || c.name === 'ELLIPSE' || c.name === 'B_CURVE') {
          if (c.name === 'LINE' && (!vs || vs === ve)) r = { why: 'unbounded straight edge' };
          else {
            let id = null;
            if (c.name === 'B_CURVE' && vs) {
              // an edge on a closed B-curve is written as the exact piece of the curve it spans, so that it never straddles the curve's seam
              const q = nurbsCurve(M, c);
              if (q.closed) {
                const ev = curveEval(M, c, memo), A0 = q.U[q.degree], B0 = q.U[q.n], T = B0 - A0; let t0 = lo !== null ? min(lo, hi) : null, t1 = lo !== null ? max(lo, hi) : null;
                if (t0 === null) { const ts = nearestParam(ev, a, A0, B0), te = vs === ve ? ts : nearestParam(ev, b, A0, B0); if (fwd) { t0 = ts; t1 = te > ts + 1e-9 * T ? te : te + T; } else { t0 = te; t1 = ts > te + 1e-9 * T ? ts : ts + T; } }
                if (!(vs === ve && abs(t0 - A0) < 1e-9 * T)) { const sc = subCurve(q, t0, t1); if (sc) { const pts = [], w = []; for (let i = 0; i < sc.n; i++) { const h = sc.P.slice(i * sc.dim, (i + 1) * sc.dim); pts.push(P(dehom(h, sc.dim, sc.rational))); if (sc.rational) w.push(h[3]); } id = bspCurve(sc.degree, pts, sc.knots, sc.mult, sc.rational ? w : null, false); stats.trimmedBCurves = (stats.trimmedBCurves || 0) + 1; } }
              }
            }
            if (!id) id = curve3(c); let va, vb; if (vs) { va = vertex(vs); vb = ve === vs ? va : vertex(ve); } else { const ev = curveEval(M, c, memo); va = vb = E(`VERTEX_POINT('',${pt(P(ev.f(ev.range[0])))})`); } r = { id: E(`EDGE_CURVE('',${va},${vb},${id},${B(fwd)})`) }; stats.exactEdges++;
          }
        } else if (c.name === 'INTERSECTION') {
          const ch = M.get(c.chart), s = M.get(c.start), t = M.get(c.end); let Q = ch && ch.hvec ? ch.hvec.filter(isVec) : [];
          if (s && s.type === 'T' && isVec(s.hvec[0])) Q = [s.hvec[0], ...Q]; if (t && t.type === 'T' && isVec(t.hvec[0])) Q = [...Q, t.hvec[0]];
          if (Q.length < 2) r = { why: 'intersection curve without chart points' };
          else { const closed = dist(Q[0], Q[Q.length - 1]) <= 1e-9; if (!fwd) Q = Q.slice().reverse(); const cut = cutPolyline(Q, closed, a, b); r = cut ? approx('intersection curve (chart points)', e, cut, vs, ve) : { why: 'intersection curve whose chart does not span the edge' }; }
        } else if (c.name === 'SP_CURVE') {
          const ev = curveEval(M, c, memo); let Q = ev ? sampleCurve(ev, lo !== null ? lo : ev.range[0], lo !== null ? hi : ev.range[1]) : null;
          if (Q && lo === null && !fwd) Q.reverse(); const cut = Q && (lo !== null ? (a && b ? [a, ...Q.slice(1, -1), b] : Q) : cutPolyline(Q, !!ev.closed, a, b));
          r = cut ? approx('surface curve (SP-curve samples)', e, cut, vs, ve) : { why: `SP-curve on a ${(M.get(c.surface) || {}).name || 'missing surface'} that cannot be evaluated` };
        } else r = { why: `${c.name} edge curve` };
      }
    } catch (err) { r = { why: `edge (${err.message})` }; }
    memo.set(key, r); return r;
  };
  const loop = lp => {
    const f0 = M.get(lp.fin); if (!f0) return { why: 'loop without fins' };
    if (!f0.edge) { const v = M.get(f0.vertex); return v ? { id: E(`VERTEX_LOOP('',${vertex(v)})`) } : { why: 'vertex loop without a vertex' }; }
    const oe = [];
    for (let f = f0, g = 0; f; g++) { if (g > 1e5) return { why: 'loop that does not close' }; const en = M.get(f.edge); if (!en) return { why: 'fin without an edge' }; const e = edge(en); if (!e.id) return e; oe.push(E(`ORIENTED_EDGE('',*,*,${e.id},${B(f.sense !== '-')})`)); f = M.get(f.forward); if (f === f0) break; }
    return { id: E(`EDGE_LOOP('',${list(oe)})`) };
  };
  const face = fa => {
    if (opts.filter && !opts.filter(fa)) return null;
    stats.faces++; const s = M.get(fa.surface); if (!s) { skip('face without a surface'); return null; }
    try {
      const sr = s.name === 'BLENDED_EDGE' ? blend(fa, s) : surface(s); if (!sr.id) { skip(sr.why); return null; }
      const bounds = []; for (let lp = M.get(fa.loop), g = 0; lp && g < 1e5; lp = M.get(lp.next), g++) { const r = loop(lp); if (!r.id) { skip(`face bounded by: ${r.why}`); return null; } bounds.push(r.id); }
      if (!bounds.length) {
        if (s.name === 'SPHERE') bounds.push(E(`VERTEX_LOOP('',${E(`VERTEX_POINT('',${pt(P(sub(s.centre, mul(s.axis, abs(s.radius)))))})`)})`));
        else if (s.name === 'TORUS') bounds.push(E(`VERTEX_LOOP('',${E(`VERTEX_POINT('',${pt(P(add(s.centre, mul(s.x_axis, s.major_radius + abs(s.minor_radius)))))})`)})`));
        else { skip(`unbounded face on a ${s.name}`); return null; }
      }
      const same = ((fa.sense === '-') === (s.sense === '-')) !== sr.flip, id = E(`ADVANCED_FACE('',${list(bounds.map((b, i) => E(`${bounds.length === 1 ? 'FACE_OUTER_BOUND' : 'FACE_BOUND'}('',${b},.T.)`)))},${sr.id},${B(same)})`);
      stats.translated++; return id;
    } catch (e) { skip(`face (${e.message})`); return null; }
  };
  const shellFaces = sh => { const out = [], seen = new Set(); for (const head of [sh.face, sh.front_face]) for (let f = M.get(head), g = 0; f && g < 1e6; f = M.get(head === sh.face ? f.next : f.next_front), g++) { if (seen.has(f.index)) break; seen.add(f.index); out.push(f); } return out; };
  const extent = fs => { let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; const grow = (p, r = 0) => { for (let a = 0; a < 3; a++) { lo[a] = min(lo[a], p[a] - r); hi[a] = max(hi[a], p[a] + r); } };
    for (const f of fs) { const s = M.get(f.surface); if (s && s.name === 'SPHERE' && isVec(s.centre)) grow(s.centre, abs(s.radius)); if (s && s.name === 'TORUS' && isVec(s.centre)) grow(s.centre, abs(s.major_radius) + abs(s.minor_radius)); for (let lp = M.get(f.loop), g = 0; lp && g < 1e4; lp = M.get(lp.next), g++) { const f0 = M.get(lp.fin); for (let fin = f0, h = 0; fin && h < 1e4; h++) { const p = vpos(M.get(fin.vertex)); if (p) grow(p); fin = M.get(fin.forward); if (fin === f0) break; } } }
    return hi[0] >= lo[0] ? hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) : 0; };
  // shared context
  const ctxApp = E("APPLICATION_CONTEXT('core data for automotive mechanical design processes')"); E(`APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,${ctxApp})`);
  const ctxProd = E(`PRODUCT_CONTEXT('',${ctxApp},'mechanical')`), ctxDef = E(`PRODUCT_DEFINITION_CONTEXT('part definition',${ctxApp},'design')`);
  const uLen = E('(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.))'), uAng = E('(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))'), uSol = E('(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())');
  let res = 1e-8; for (const n of M.values()) if ((n.name === 'BODY' || n.name === 'ASSEMBLY') && n.res_linear > 0 && n.res_linear < 1e-3) { res = n.res_linear; break; }
  const unc = E(`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(${real(res)}),${uLen},'distance_accuracy_value','confusion accuracy')`);
  const ctx = E(`(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((${unc}))GLOBAL_UNIT_ASSIGNED_CONTEXT((${uLen},${uAng},${uSol}))REPRESENTATION_CONTEXT('Context #1','3D'))`);
  const origin = () => E(`AXIS2_PLACEMENT_3D('',${pt([0, 0, 0])},${dir([0, 0, 1])},${dir([1, 0, 0])})`);
  const stem = opts.name ? String(opts.name).replace(/\.[^.]+$/, '') : 'body';
  const body = (b, label) => {
    const type = BODY_TYPES[b.body_type] || `type ${b.body_type}`, nm = attr(M, b, 'SDL/TYSA_NAME'), name = (nm && typeof nm[0] === 'string' && nm[0].trim()) || label, rec = { name, type, faces: 0, translated: 0, closed: false, colour: null }; stats.bodies++;
    const regions = []; for (let r = M.get(b.region), g = 0; r && g < 1e5; r = M.get(r.next), g++) regions.push(r);
    const shellsOf = r => { const out = []; for (let s = M.get(r.shell), g = 0; s && g < 1e5; s = M.get(s.next), g++) out.push(s); return out; };
    const solid = regions.filter(r => r.type === 'S'), items = [], open = [], f0 = stats.faces, t0 = stats.translated, colours = new Map(); let complete = true, any = false;
    const tally = fs => { for (const f of fs) { const c = attr(M, f, 'SDL/TYSA_COLOUR'); if (c && c[0] && c[0].length >= 3) { const k = c[0].slice(0, 3).join(','); colours.set(k, (colours.get(k) || 0) + 1); } } };
    if (solid.length) {
      for (const r of solid) {
        const sh = shellsOf(r).map(s => { const fs = []; for (let f = M.get(s.face), g = 0; f && g < 1e6; f = M.get(f.next), g++) fs.push(f); return fs; }).filter(fs => fs.length); if (!sh.length) continue;
        if (sh.length > 1) { const ex = sh.map(extent); let o = 0; ex.forEach((x, i) => { if (x > ex[o]) o = i; }); if (o) sh.unshift(sh.splice(o, 1)[0]); }
        const done = sh.map(fs => { tally(fs); const ids = fs.map(face).filter(Boolean); return { ids, all: ids.length === fs.length }; });
        if (done.every(d => d.all)) { const cs = done.map(d => E(`CLOSED_SHELL('',${list(d.ids)})`)); items.push(cs.length === 1 ? E(`MANIFOLD_SOLID_BREP(${str(name)},${cs[0]})`) : E(`BREP_WITH_VOIDS(${str(name)},${cs[0]},${list(cs.slice(1).map(c => E(`ORIENTED_CLOSED_SHELL('',*,${c},.F.)`)))})`)); any = true; }
        else { complete = false; for (const d of done) if (d.ids.length) open.push(E(`OPEN_SHELL('',${list(d.ids)})`)); }
      }
    } else { const fs = []; for (const r of regions) for (const s of shellsOf(r)) for (const f of shellFaces(s)) if (!fs.includes(f)) fs.push(f); tally(fs); const ids = fs.map(face).filter(Boolean); if (ids.length) open.push(E(`OPEN_SHELL('',${list(ids)})`)); complete = false; if (!fs.length) { stats.wireBodies++; } }
    rec.faces = stats.faces - f0; rec.translated = stats.translated - t0; rec.closed = !!solid.length && complete && any;
    const bc = attr(M, b, 'SDL/TYSA_COLOUR'); let best = bc && bc[0] && bc[0].length >= 3 ? bc[0].slice(0, 3).join(',') : null; if (!best) { let n = 0; for (const [k, c] of colours) if (c > n) { n = c; best = k; } } if (best) rec.colour = best.split(',').map(Number);
    bodies.push(rec); if (!items.length && !open.length) return;
    if (solid.length && rec.closed) stats.solids++; else stats.sheets++;
    const prod = E(`PRODUCT(${str(name)},${str(name)},'',(${ctxProd}))`), pdf = E(`PRODUCT_DEFINITION_FORMATION('','',${prod})`), pd = E(`PRODUCT_DEFINITION('design','',${pdf},${ctxDef})`), pds = E(`PRODUCT_DEFINITION_SHAPE('','',${pd})`);
    E(`PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(${prod}))`);
    const all = [...items]; if (open.length) all.push(E(`SHELL_BASED_SURFACE_MODEL(${str(name)},${list(open)})`));
    const rep = E(`${open.length ? 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION' : 'ADVANCED_BREP_SHAPE_REPRESENTATION'}('',${list([origin(), ...all])},${ctx})`); E(`SHAPE_DEFINITION_REPRESENTATION(${pds},${rep})`);
    rec.written = true;
  };
  // parts reachable from the root: a body, a partition's world, an assembly with instances, or a list of parts
  const compose = (outer, t) => { if (!t) return outer; const R = t.rotation_matrix, v = isVec(t.translation_vector) ? t.translation_vector : [0, 0, 0], s = t.scale > 0 ? t.scale : 1, inner = { R, t: v, s }; if (!outer) return inner; const A = outer.R, C = new Array(9); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[3 * i + j] = A[3 * i] * R[j] + A[3 * i + 1] * R[3 + j] + A[3 * i + 2] * R[6 + j]; const sv = mul(v, s); return { R: C, t: [A[0] * sv[0] + A[1] * sv[1] + A[2] * sv[2] + outer.t[0], A[3] * sv[0] + A[4] * sv[1] + A[5] * sv[2] + outer.t[1], A[6] * sv[0] + A[7] * sv[1] + A[8] * sv[2] + outer.t[2]], s: outer.s * s }; };
  const det = R => R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]);
  let serial = 0, mirrored = 0; const visiting = new Set();
  const part = (n, t, depth) => {
    if (!n || depth > 32) return;
    if (n.name === 'BODY') { tf = t; tfKey = t ? '@' + (++serial) : ''; body(n, `${stem}${stats.bodies ? ' ' + (stats.bodies + 1) : ''}`); tf = null; tfKey = ''; }
    else if (n.name === 'ASSEMBLY') { if (visiting.has(n.index)) return; visiting.add(n.index); for (let i = M.get(n.sub_instance), g = 0; i && g < 1e5; i = M.get(i.next_in_part), g++) { const x = compose(t, M.get(i.transform)); stats.instances++; if (x && det(x.R) < 0) { mirrored++; continue; } part(M.get(i.part), x, depth + 1); } visiting.delete(n.index); }
  };
  const root = model.root;
  if (root && root.name === 'WORLD') { for (let b = M.get(root.body), g = 0; b && g < 1e5; b = M.get(b.next), g++) part(b, null, 0); if (root.assembly) for (let a = M.get(root.assembly), g = 0; a && g < 1e5; a = M.get(a.next), g++) if (!a.ref_instance) part(a, null, 0); }
  else if (root && root.name === 'POINTER_LIS_BLOCK') { for (let blk = root, g = 0; blk && g < 1e5; blk = M.get(blk.next_block), g++) for (const i of blk.entries || []) part(M.get(i), null, 0); }
  else if (root) part(root, null, 0);
  if (!stats.bodies) for (const n of M.values()) if (n.name === 'BODY') part(n, null, 0);
  if (mirrored) lost.push(`${mirrored} mirrored assembly instance(s) were left out (a reflection cannot be applied to the analytic surface placements)`);
  const written = bodies.some(b => b.written);
  const head = `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('Parasolid XT B-rep translated entity by entity'),'2;1');\nFILE_NAME(${str(stem)},'',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n`;
  const step = written ? head + L.map((s, i) => `#${i + 1}=${s};`).join('\n') + '\nENDSEC;\nEND-ISO-10303-21;\n' : null;
  return { step, bodies, stats, notes, lost };
}
// The NURBS and polyline kit and the STEP literal writers, shared with the ACIS translator (sat.js)
export const kit = { deBoor, flatKnots, insertKnot, subCurve, nearestParam, interpolate, tidy, locate, cutPolyline, sampleCurve, real, str };
