// ACIS / Autodesk ShapeManager save files: SAT text (.sat), SAB binary (.sab) and the "ASM BinaryFile" streams that
// Inventor and Fusion embed (.smb / .smbh and the B-rep segments of .ipt), read record by record and translated entity
// by entity into an ISO 10303-21 (STEP AP214) advanced B-rep text. Nothing is evaluated or tessellated here: trimming,
// healing and meshing are left to the OpenCASCADE kernel that reads the STEP text. Record layouts follow the published
// SAT save-file description (topology, analytic geometry, B-spline data blocks) and what the files themselves show.
import { kit } from './xt.js';

const { abs, min, max, hypot, atan2, sqrt, PI, cos, sin } = Math;
const { deBoor, flatKnots, subCurve, nearestParam, real, str } = kit;

// ---------------------------------------------------------------- tokens
// A record is {i, name, base, t}: `name` is the full dash-joined class name ("plane-surface"), `base` its last part
// and `t` the tokens after it: numbers, Ptr (entity reference), Str (quoted string), Enum (binary enumeration value),
// true / false (binary logicals), and bare words ("forward", "nubs", "{", "}" …) as plain strings.
export class Ptr { constructor(i) { this.i = i; } }
export class Str { constructor(s) { this.s = s; } }
export class Enum { constructor(v) { this.v = v; } }
const NUM = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const latin = (b, o, n) => { let s = ''; for (let i = o, e = min(b.length, o + n); i < e; i++) s += String.fromCharCode(b[i]); return s; };
const BIN_SIGS = ['ACIS BinaryFile', 'ASM BinaryFile4', 'ASM BinaryFile8'];
const MARK_END = /^End-of-[A-Za-z]+-data$/i, MARK_HIST = /^Begin-of-[A-Za-z]+-History/i;

// Offsets of binary ACIS / ShapeManager streams inside a larger buffer, found by signature → [{offset, signature}]
export function findACIS(bytes, limit = 256) {
  const b = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes), out = [];
  for (let i = 0; i + 40 < b.length && out.length < limit; i++) {
    if (b[i] !== 0x41) continue; // 'A'
    if (b[i + 1] === 0x43 && b[i + 2] === 0x49 && b[i + 3] === 0x53 && latin(b, i, 15) === BIN_SIGS[0]) out.push({ offset: i, signature: BIN_SIGS[0] });
    else if (b[i + 1] === 0x53 && b[i + 2] === 0x4D && b[i + 3] === 0x20) { const s = latin(b, i, 15); if (s === BIN_SIGS[1] || s === BIN_SIGS[2]) out.push({ offset: i, signature: s }); }
  }
  return out;
}

function readText(b, model) {
  let s = ''; const CH = 8192; for (let i = 0; i < b.length; i += CH) s += String.fromCharCode.apply(null, b.subarray(i, min(b.length, i + CH)));
  if (s.charCodeAt(0) === 0xEF && s.charCodeAt(1) === 0xBB) s = s.slice(3);
  let o = 0; const n = s.length, ws = c => c === 32 || c === 9 || c === 10 || c === 13;
  const raw = () => { while (o < n && ws(s.charCodeAt(o))) o++; if (o >= n) return null; const a = o; while (o < n && !ws(s.charCodeAt(o))) o++; return s.slice(a, o); };
  const hnum = what => { const t = raw(); if (t === null || !NUM.test(t)) throw new Error(`SAT: the header ${what} is missing`); return +t; };
  model.version = hnum('version'); model.records = hnum('record count'); model.bodyCount = hnum('entity count'); model.flags = hnum('history flag');
  if (!(model.version >= 100 && model.version < 1e6 && Number.isInteger(model.version))) throw new Error(`SAT: implausible save-file version ${model.version}`);
  if (model.version >= 200) {
    const lstr = () => { const t = raw(); if (t === null || !/^@?\d+$/.test(t)) throw new Error('SAT: the header strings are missing'); const len = +t.replace('@', ''); o++; const r = s.slice(o, o + len); o += len; return r; };
    model.product = lstr(); model.acis = lstr(); model.date = lstr(); model.scale = hnum('unit scale'); model.resabs = hnum('resabs'); model.resnor = hnum('resnor');
  }
  // records
  const ents = model.entities, bare = model.version < 104; let t = [], seq = 0, first = true, done = false;
  const flush = () => {
    let k = 0, idx = null;
    if (first) { const ref = x => x instanceof Ptr || (bare && typeof x === 'number'); while (k < t.length && !(typeof t[k] === 'string' && /^[A-Za-z_]/.test(t[k]) && ref(t[k + 1])) && !(typeof t[k] === 'number' && typeof t[k + 1] === 'string' && ref(t[k + 2]))) k++; first = false; }
    if (typeof t[k] === 'number' && t[k] <= 0 && Number.isInteger(t[k])) { idx = -t[k]; k++; }
    const name = t[k]; if (typeof name !== 'string' || !/^[A-Za-z_]/.test(name)) { if (t.length) model.warnings.push(`a record without a class name was skipped (near record ${seq})`); t = []; return; }
    if (idx === null) idx = seq; seq = idx + 1; ents.set(idx, { i: idx, name, base: name.slice(name.lastIndexOf('-') + 1), t: t.slice(k + 1) }); t = [];
  };
  while (o < n && !done) {
    const c = s.charCodeAt(o);
    if (ws(c)) { o++; continue; }
    if (c === 35) { o++; flush(); continue; } // #
    if (c === 123 || c === 125) { t.push(s[o++]); continue; } // { }
    if (c === 36) { const a = ++o; while (o < n && !ws(s.charCodeAt(o)) && s.charCodeAt(o) !== 35) o++; const v = +s.slice(a, o); if (!Number.isInteger(v)) throw new Error(`SAT: malformed entity reference near offset ${a}`); t.push(new Ptr(v)); continue; }
    if (c === 64) { const a = ++o; while (o < n && !ws(s.charCodeAt(o))) o++; const len = +s.slice(a, o); if (!Number.isInteger(len) || len < 0) throw new Error(`SAT: malformed string length near offset ${a}`); o++; t.push(new Str(s.slice(o, o + len))); o += len; continue; }
    const a = o; while (o < n && !ws(s.charCodeAt(o))) o++; let w = s.slice(a, o); let hash = false;
    if (w.length > 1 && w.endsWith('#')) { w = w.slice(0, -1); hash = true; }
    if (NUM.test(w)) t.push(+w); else if (!t.length && MARK_END.test(w)) done = true; else if (!t.length && MARK_HIST.test(w)) { done = true; model.history = true; } else t.push(w);
    if (hash) flush();
  }
  model.terminated = done;
}

function readBinary(b, start, model) {
  const v = new DataView(b.buffer, b.byteOffset, b.byteLength), sig = latin(b, start, 15); let o = start + 15, W = sig === BIN_SIGS[2] ? 8 : 4;
  const need = k => { if (o + k > b.length) throw new Error('SAB: unexpected end of the binary data'); };
  const long = () => { need(W); const x = W === 8 ? Number(v.getBigInt64(o, true)) : v.getInt32(o, true); o += W; return x; };
  // the four header counts are as wide as the file's integers (eight bytes in the "…File8" variant)
  if (W === 8) { need(32); model.version = Number(v.getBigInt64(o, true)); model.records = Number(v.getBigInt64(o + 8, true)); model.bodyCount = Number(v.getBigInt64(o + 16, true)); model.flags = Number(v.getBigInt64(o + 24, true)); o += 32; }
  else { need(16); model.version = v.getInt32(o, true); model.records = v.getInt32(o + 4, true); model.bodyCount = v.getInt32(o + 8, true); model.flags = v.getInt32(o + 12, true); o += 16; }
  if (!(model.version >= 100 && model.version < 1e6)) throw new Error(`SAB: implausible save-file version ${model.version}`);
  const ents = model.entities; let t = [], name = '', pre = '', seq = 0, head = model.version >= 200 ? 6 : 0; const hd = [];
  for (let guard = 0; o < b.length; guard++) {
    const tag = b[o++]; let x;
    switch (tag) {
      case 2: need(1); x = v.getInt8(o); o += 1; break;
      case 3: need(2); x = v.getInt16(o, true); o += 2; break;
      case 4: x = long(); break;
      case 5: need(4); x = v.getFloat32(o, true); o += 4; break;
      case 6: need(8); x = v.getFloat64(o, true); o += 8; break;
      case 7: { need(1); const l = b[o]; need(1 + l); x = new Str(latin(b, o + 1, l)); o += 1 + l; break; }
      case 8: { need(2); const l = v.getUint16(o, true); need(2 + l); x = new Str(latin(b, o + 2, l)); o += 2 + l; break; }
      case 9: case 18: { need(4); const l = v.getUint32(o, true); need(4 + l); x = new Str(latin(b, o + 4, l)); o += 4 + l; break; }
      case 10: x = true; break;
      case 11: x = false; break;
      case 12: x = new Ptr(long()); break;
      case 13: case 14: { need(1); const l = b[o]; need(1 + l); const w = latin(b, o + 1, l); o += 1 + l; if (tag === 14) { pre += w + '-'; continue; } x = pre + w; pre = ''; break; }
      case 15: x = '{'; break;
      case 16: x = '}'; break;
      case 17: if (name) { ents.set(seq, { i: seq, name, base: name.slice(name.lastIndexOf('-') + 1), t }); seq++; } else if (t.length) model.warnings.push(`a record without a class name was skipped (near record ${seq})`); t = []; name = ''; continue;
      case 19: case 20: need(24); t.push(v.getFloat64(o, true), v.getFloat64(o + 8, true), v.getFloat64(o + 16, true)); o += 24; continue;
      case 21: x = new Enum(long()); break;
      case 22: need(16); t.push(v.getFloat64(o, true), v.getFloat64(o + 8, true)); o += 16; continue;
      case 23: need(8); x = Number(v.getBigInt64(o, true)); o += 8; break;
      default: throw new Error(`SAB: unknown data tag ${tag} at offset ${o - 1} (record ${seq})`);
    }
    if (head > 0) { hd.push(x); if (--head === 0) { model.product = hd[0] instanceof Str ? hd[0].s : ''; model.acis = hd[1] instanceof Str ? hd[1].s : ''; model.date = hd[2] instanceof Str ? hd[2].s : ''; model.scale = +hd[3]; model.resabs = +hd[4]; model.resnor = +hd[5]; } continue; }
    if (!name && !t.length && typeof x === 'string' && tag === 13) {
      if (MARK_END.test(x)) { model.terminated = true; break; }
      if (MARK_HIST.test(x)) { model.terminated = true; model.history = true; break; }
      name = x; continue;
    }
    t.push(x);
  }
  model.end = o;
}

// readSAT(bytes, {offset}) → {format: 'sat' | 'sab', signature, version, product, acis, date, scale (millimetres per
// model unit, null when the file predates the units record), resabs, resnor, entities: Map(index → record), counts,
// history (a roll-back history section follows the model and was left unread), warnings}.
export function readSAT(bytes, opts = {}) {
  const b = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes), start = opts.offset || 0, sig = latin(b, start, 15);
  const model = { format: 'sat', signature: null, version: 0, records: 0, bodyCount: 0, flags: 0, product: '', acis: '', date: '', scale: null, resabs: null, resnor: null, entities: new Map(), counts: {}, history: false, terminated: false, warnings: [] };
  if (BIN_SIGS.includes(sig)) { model.format = 'sab'; model.signature = sig; readBinary(b, start, model); }
  else readText(start ? b.subarray(start) : b, model);
  if (!(model.scale > 0)) model.scale = null;
  for (const e of model.entities.values()) model.counts[e.name] = (model.counts[e.name] || 0) + 1;
  if (!model.terminated && model.version >= 200) model.warnings.push('the end-of-data marker was not found: the file may be truncated');
  return model;
}

// ---------------------------------------------------------------- record access
class Cur {
  constructor(t, i = 0, bare = false) { this.t = t; this.i = i; this.bare = bare; }
  peek() { return this.t[this.i]; }
  // an entity reference; the earliest save files (1.x) write it as a bare integer
  ptr() { const x = this.t[this.i]; if (x instanceof Ptr) { this.i++; return x.i; } if (this.bare && typeof x === 'number' && Number.isInteger(x)) { this.i++; return x; } throw new Error('an entity reference is missing'); }
  isPtr() { const x = this.t[this.i]; return x instanceof Ptr; }
  num() { const x = this.t[this.i]; if (typeof x !== 'number') throw new Error('a number is missing'); this.i++; return x; }
  vec() { return [this.num(), this.num(), this.num()]; }
  // a logical saved as a word pair in text ("forward" / "reversed") and as a true / false tag in binary
  flag(f, t) { const x = this.t[this.i]; let r; if (x === true || x === false) r = x; else if (x === t) r = true; else if (x === f) r = false; else if (x instanceof Enum && (x.v === 0 || x.v === 1)) r = x.v === 1; else if (x === 0 || x === 1) r = x === 1; else throw new Error(`"${f}" or "${t}" is missing`); this.i++; return r; }
  isFlag(f, t) { const x = this.t[this.i]; return x === true || x === false || x === f || x === t; }
  // one end of a parameter interval: unbounded ("I") or a finite value ("F" value)
  ival() { const x = this.t[this.i]; if (x === 'I' || x === false) { this.i++; return null; } if (x === 'F' || x === true) { this.i++; return this.num(); } throw new Error('a parameter bound is missing'); }
  isIval() { const x = this.t[this.i]; return x === 'I' || x === 'F' || x === true || x === false; }
  skipNums() { while (typeof this.t[this.i] === 'number') this.i++; }
}
const CLOSURE = { open: 0, closed: 1, periodic: 2 };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2], cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = a => hypot(a[0], a[1], a[2]), unit = a => { const l = len(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const isVec = v => Array.isArray(v) && v.length === 3 && v.every(Number.isFinite);

// The save-file model with decoding of the records the translator needs. Records are decoded on demand and memoised.
export function access(model) {
  const E = model.entities, newStyle = model.version >= 700, memo = new Map();
  // attrib pointer, then (7.0 and later) the entity id number(s) and one further reference, absent from attribute records
  const bare = model.version < 104;
  const head = (e, attrib) => { const c = new Cur(e.t, 0, bare); const attr = c.ptr(); if (newStyle) { if (typeof c.peek() === 'number') c.i++; if (typeof c.peek() === 'number' && c.t[c.i + 1] instanceof Ptr) c.i++; if (!attrib && c.isPtr()) c.i++; } c.attr = attr; return c; };
  // the run of entity references that follows; `n` is how many a record of this class holds when they are bare integers
  const ptrs = (c, n) => { const p = []; if (bare) { for (let i = 0; i < n; i++) p.push(c.ptr()); return p; } while (c.isPtr()) p.push(c.ptr()); return p; };
  // subtype objects ("{ name … }") in file order: a later "{ ref n }" names the n-th of them
  let table = null, byId = null, useId = false;
  const subtypes = () => {
    if (table) return table; table = []; byId = new Map(); let withId = 0, same = 0;
    for (const e of E.values()) { const t = e.t; for (let i = 0; i < t.length; i++) if (t[i] === '{' && typeof t[i + 1] === 'string' && t[i + 1] !== 'ref') { const id = typeof t[i + 2] === 'number' && Number.isInteger(t[i + 2]) ? t[i + 2] : null; if (id !== null) { withId++; if (id === table.length) same++; if (!byId.has(id)) byId.set(id, table.length); } table.push({ t, i }); } }
    useId = table.length > 0 && withId === table.length && same < table.length && byId.size === table.length;
    return table;
  };
  // the token span [a, b) of the subtype object that opens at t[i], following "ref" indirections
  const span = (t, i) => {
    for (let g = 0; g < 8; g++) {
      if (t[i] !== '{') throw new Error('a subtype block is missing');
      let d = 0, j = i; for (; j < t.length; j++) { if (t[j] === '{') d++; else if (t[j] === '}' && --d === 0) break; } if (j >= t.length) throw new Error('a subtype block does not close');
      if (t[i + 1] !== 'ref') return { t, a: i + 1, b: j, end: j + 1 };
      const tb = subtypes(), n = t[i + 2], k = useId ? byId.get(n) : n, r = tb[k]; if (!r) throw new Error(`subtype reference ${n} has no target`);
      const s = span(r.t, r.i); return { t: s.t, a: s.a, b: s.b, end: j + 1 };
    }
    throw new Error('subtype references nest too deeply');
  };
  // first B-spline block at the top level of a subtype span: the index of its "nubs" / "nurbs" word, or -1
  const findBs = (s, wantSurface) => {
    let d = 0, k0 = s.a + 1; if (typeof s.t[k0] === 'number' && Number.isInteger(s.t[k0]) && typeof s.t[k0 + 1] !== 'number') k0++;
    // the approximation level directly after the subtype name: "summary" and "none" carry no control points
    const lv = s.t[k0]; if (lv === 'summary' || lv === 'none' || (lv instanceof Enum && lv.v !== 0 && (s.t[k0 + 1] !== 'nubs' && s.t[k0 + 1] !== 'nurbs'))) return -1;
    for (let k = s.a; k < s.b; k++) {
      const x = s.t[k]; if (x === '{') d++; else if (x === '}') d--;
      else if (d === 0 && (x === 'nubs' || x === 'nurbs') && typeof s.t[k + 1] === 'number') {
        const surf = typeof s.t[k + 2] === 'number'; if (surf !== wantSurface) { if (!wantSurface) return -1; continue; }
        const p = s.t[k - 1]; if (p === 'none' || p === 'summary' || (p instanceof Enum && p.v !== 0)) return -1;
        return k;
      }
    }
    return -1;
  };
  const knots = (c, n) => { const k = [], m = []; for (let i = 0; i < n; i++) { k.push(c.num()); const q = c.num(); if (!(Number.isInteger(q) && q >= 1 && q <= 64)) throw new Error('implausible knot multiplicity'); m.push(q); } for (let i = 1; i < n; i++) if (!(k[i] > k[i - 1])) throw new Error('knots do not increase'); return [k, m]; };
  const closure = c => { const x = c.t[c.i++]; if (x instanceof Enum) return x.v; if (typeof x === 'string' && x in CLOSURE) return CLOSURE[x]; throw new Error('a closure word is missing'); };
  const total = m => m.reduce((a, x) => a + x, 0);
  // B-spline curve block → {degree, rational, n, P: [[x, y, z]], W | null, knots, mult (clamping knots restored), closure, fitol}
  const bsCurve = (t, k) => {
    const c = new Cur(t, k + 1), rational = t[k] === 'nurbs', degree = c.num(), cl = closure(c), nk = c.num(); if (!(degree >= 1 && degree <= 25 && nk >= 2 && nk < 1e6)) throw new Error('implausible B-spline curve header');
    const [kn, mu] = knots(c, nk), n = total(mu) - degree + 1; if (n < degree + 1) throw new Error('B-spline curve with too few control points');
    const P = [], W = rational ? [] : null; for (let i = 0; i < n; i++) { P.push(c.vec()); if (rational) W.push(c.num()); }
    const fitol = typeof c.peek() === 'number' ? c.num() : null, mult = mu.slice(); mult[0]++; mult[nk - 1]++;
    return { degree, rational, n, P, W, knots: kn, mult, closure: cl, fitol };
  };
  // B-spline surface block → {du, dv, nu, nv, P[u][v], W[u][v] | null, ku, mu, kv, mv, cu, cv, fitol}
  const bsSurface = (t, k) => {
    const c = new Cur(t, k + 1), rational = t[k] === 'nurbs', du = c.num(), dv = c.num(), fl = []; while (c.i < t.length && typeof c.peek() !== 'number' && fl.length < 6) fl.push(c.t[c.i++]);
    if (fl.length !== 4 && fl.length !== 5) throw new Error('unrecognised B-spline surface header'); const o = fl.length - 4, clo = x => x instanceof Enum ? x.v : (x in CLOSURE ? CLOSURE[x] : 0);
    const nku = c.num(), nkv = c.num(); if (!(du >= 1 && du <= 25 && dv >= 1 && dv <= 25 && nku >= 2 && nkv >= 2 && nku < 1e6 && nkv < 1e6)) throw new Error('implausible B-spline surface header');
    const [ku, mu0] = knots(c, nku), [kv, mv0] = knots(c, nkv), nu = total(mu0) - du + 1, nv = total(mv0) - dv + 1; if (nu < du + 1 || nv < dv + 1 || nu * nv > 4e6) throw new Error('B-spline surface with an implausible control net');
    // control points are saved with the u index varying fastest
    const P = Array.from({ length: nu }, () => new Array(nv)), W = rational ? Array.from({ length: nu }, () => new Array(nv)) : null;
    for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) { P[i][j] = c.vec(); if (rational) W[i][j] = c.num(); }
    const fitol = typeof c.peek() === 'number' ? c.num() : null, mu = mu0.slice(), mv = mv0.slice(); mu[0]++; mu[nku - 1]++; mv[0]++; mv[nkv - 1]++;
    return { du, dv, rational, nu, nv, P, W, ku, mu, kv, mv, cu: clo(fl[o]), cv: clo(fl[o + 1]), fitol };
  };
  // a curve written inline inside a subtype block ("straight …", "ellipse …", "intcurve forward { … }") → [curve, next index]
  const inlineCurve = (t, k) => {
    const c = new Cur(t, k + 1), w = t[k]; let r;
    if (w === 'straight') r = { kind: 'line', p: c.vec(), d: c.vec() };
    else if (w === 'ellipse') r = { kind: 'ellipse', c: c.vec(), n: c.vec(), m: c.vec(), ratio: c.num() };
    else if (w === 'intcurve') { const rev = c.flag('forward', 'reversed'), s = span(c.t, c.i), name = String(s.t[s.a]), kk = findBs(s, false); c.i = s.end; r = kk < 0 ? { kind: 'none', why: `intcurve (${name}) without a B-spline` } : { kind: 'bspline', rev, sub: name, exact: EXACT_CUR.has(name), q: bsCurve(s.t, kk) }; }
    else throw new Error(`inline ${typeof w === 'string' ? w : 'curve'} definition not understood`);
    while (c.isIval()) c.ival(); return [r, c.i];
  };
  const EXACT_CUR = new Set(['exactcur', 'exact_int_cur']), EXACT_SUR = new Set(['exactsur', 'exact_spl_sur']);
  const decode = e => {
    switch (e.base) {
      case 'body': { const c = head(e), p = ptrs(c, 3); return { attr: c.attr, lump: p[0], wire: p[1], transform: p[2] }; }
      case 'lump': { const c = head(e), p = ptrs(c, 3); return { attr: c.attr, next: p[0], shell: p[1], body: p[2] }; }
      case 'shell': { const c = head(e), p = ptrs(c, 4); return p.length >= 5 ? { attr: c.attr, next: p[0], subshell: p[1], face: p[2], wire: p[3], lump: p[4] } : { attr: c.attr, next: p[0], subshell: p[1], face: p[2], wire: -1, lump: p[3] }; }
      case 'subshell': { const c = head(e), p = ptrs(c, 4); return { attr: c.attr, parent: p[0], next: p[1], child: p[2], face: p[3], wire: p.length > 4 ? p[4] : -1 }; }
      case 'face': { const c = head(e), next = c.ptr(), loop = c.ptr(), shell = c.ptr(), subshell = c.ptr(), surface = c.ptr(), rev = c.flag('forward', 'reversed'), dbl = c.isFlag('single', 'double') ? c.flag('single', 'double') : false; return { attr: c.attr, next, loop, shell, subshell, surface, rev, double: dbl }; }
      case 'loop': { const c = head(e), p = ptrs(c, 3); return { attr: c.attr, next: p[0], coedge: p[1], face: p[2] }; }
      case 'coedge': { const c = head(e), next = c.ptr(), prev = c.ptr(), partner = c.ptr(), edge = c.ptr(), rev = c.flag('forward', 'reversed'), owner = c.ptr(); return { attr: c.attr, next, prev, partner, edge, rev, loop: owner }; }
      case 'edge': { const c = head(e), start = c.ptr(), t0 = !bare && typeof c.peek() === 'number' ? c.num() : null, end = c.ptr(), t1 = !bare && typeof c.peek() === 'number' ? c.num() : null, coedge = c.ptr(), curve = c.ptr(), rev = c.flag('forward', 'reversed'); return { attr: c.attr, start, end, t0, t1, coedge, curve, rev }; }
      case 'vertex': { const c = head(e), edge = c.ptr(); if (!bare) c.skipNums(); return { attr: c.attr, edge, point: c.ptr() }; }
      case 'point': { const c = head(e); return { p: c.vec() }; }
      case 'transform': { const c = head(e), R = []; for (let i = 0; i < 9; i++) R.push(c.num()); const t = c.vec(), s = c.num(); let reflect = false; if (c.isFlag('no_rotate', 'rotate')) { c.flag('no_rotate', 'rotate'); if (c.isFlag('no_reflect', 'reflect')) reflect = c.flag('no_reflect', 'reflect'); } return { R, t, s, reflect }; }
      case 'curve': {
        const c = head(e), kind = e.name.slice(0, -6);
        if (kind === 'straight') return { kind: 'line', p: c.vec(), d: c.vec() };
        if (kind === 'ellipse') return { kind: 'ellipse', c: c.vec(), n: c.vec(), m: c.vec(), ratio: c.num() };
        if (kind === 'intcurve') {
          const rev = c.flag('forward', 'reversed'), s = span(c.t, c.i), name = String(s.t[s.a]), k = findBs(s, false);
          if (k < 0) return { kind: 'none', why: `intcurve (${name}) saved without its approximating B-spline` };
          return { kind: 'bspline', rev, sub: name, exact: EXACT_CUR.has(name), q: bsCurve(s.t, k) };
        }
        return { kind: 'none', why: `${e.name}` };
      }
      case 'surface': {
        const c = head(e), kind = e.name.slice(0, -8);
        if (kind === 'plane') return { kind: 'plane', p: c.vec(), n: c.vec(), u: typeof c.peek() === 'number' && typeof c.t[c.i + 2] === 'number' ? c.vec() : null };
        if (kind === 'cone') { const o = c.vec(), n = c.vec(), m = c.vec(), ratio = c.num(); while (c.isIval()) c.ival(); const sn = c.num(), cs = c.num(); return { kind: 'cone', c: o, n, m, ratio, sin: sn, cos: cs }; }
        if (kind === 'sphere') { const o = c.vec(), r = c.num(), u = typeof c.peek() === 'number' ? c.vec() : null, pole = typeof c.peek() === 'number' ? c.vec() : null; return { kind: 'sphere', c: o, r, u, pole }; }
        if (kind === 'torus') { const o = c.vec(), n = c.vec(), R = c.num(), r = c.num(), u = typeof c.peek() === 'number' ? c.vec() : null; return { kind: 'torus', c: o, n, R, r, u }; }
        if (kind === 'spline') {
          const rev = c.flag('forward', 'reversed'), s = span(c.t, c.i), name = String(s.t[s.a]), k = findBs(s, true);
          if (k < 0 && (name === 'sumsur' || name === 'sum_spl_sur')) {
            // a sum surface c1(u) + c2(v) − origin of a straight line and another curve is the linear extrusion of that curve
            let i = s.a + 1; if (typeof s.t[i] === 'number') i++; const [c1, n1] = inlineCurve(s.t, i), [c2, n2] = inlineCurve(s.t, n1), line = c1.kind === 'line' ? c1 : c2.kind === 'line' ? c2 : null, prof = line === c1 ? c2 : c1;
            if (line && (prof.kind === 'ellipse' || prof.kind === 'bspline')) {
              // exact when the saved origin (or, failing that, an end of the profile) lies on the straight line
              const d = unit(line.d), onLine = q => { const x = sub(line.p, q); return len(sub(x, mul(d, dot(x, d)))) <= 1e-6 * max(1, len(q), len(line.p)); }, o = typeof s.t[n2] === 'number' && typeof s.t[n2 + 1] === 'number' && typeof s.t[n2 + 2] === 'number' ? [s.t[n2], s.t[n2 + 1], s.t[n2 + 2]] : null;
              if (o ? onLine(o) : prof.kind === 'bspline' ? onLine(prof.q.P[0]) || onLine(prof.q.P[prof.q.n - 1]) : onLine(add(prof.c, prof.m))) return { kind: 'extrusion', rev, sub: name, profile: prof, d: line.d, lineFirst: line === c1 };
            }
            return { kind: 'none', why: `spline surface (${name}) that is not a straight extrusion and has no approximating B-spline` };
          }
          if (k < 0) return { kind: 'none', why: `spline surface (${name}) saved without its approximating B-spline` };
          return { kind: 'bspline', rev, sub: name, exact: EXACT_SUR.has(name), q: bsSurface(s.t, k) };
        }
        return { kind: 'none', why: `${e.name}` };
      }
      default: return null;
    }
  };
  const get = i => { if (!(i >= 0)) return null; if (memo.has(i)) return memo.get(i); const e = E.get(i); let r = null; if (e) { try { r = decode(e); if (r) { r.i = i; r.name = e.name; r.base = e.base; } } catch (err) { r = { i, name: e.name, base: e.base, bad: err.message, kind: 'none', why: `${e.name} record that could not be decoded (${err.message})` }; } } memo.set(i, r); return r; };
  // attribute records chained from an entity: [{name, t (tokens after the owner reference)}]
  const attributes = attr => { const out = []; for (let i = attr, g = 0; i >= 0 && g < 4096; g++) { const e = E.get(i); if (!e || e.base !== 'attrib') break; let c; try { c = head(e, true); const next = c.ptr(); c.ptr(); c.ptr(); out.push({ name: e.name, t: e.t.slice(c.i) }); i = next; } catch (err) { break; } } return out; };
  return { get, attributes, E, subtypes, span, findBs, bsSurface, bsCurve };
}
const utf8 = s => { if (!/[\x80-\xff]/.test(s)) return s; try { return new TextDecoder('utf-8', { fatal: true }).decode(Uint8Array.from(s, c => c.charCodeAt(0))); } catch (e) { return s; } };

// ---------------------------------------------------------------- ACIS → STEP
const B = v => v ? '.T.' : '.F.';
// satToStep(model, opts) → {step: AP214 text | null, bodies: [{name, type, faces, translated, closed, colour}], stats, notes, lost}.
// `model` is the result of readSAT. Every face whose surface and edge curves have a STEP form is written as an
// ADVANCED_FACE; a lump whose shells are topologically closed and whose faces all translate becomes a
// MANIFOLD_SOLID_BREP (BREP_WITH_VOIDS for cavities), otherwise its faces go into a SHELL_BASED_SURFACE_MODEL. Lengths
// are written in millimetres: model units × the header's millimetres-per-unit (1 when the file has no units record).
// opts: name (product name stem), bodies (indices of the body records to translate; default all), scale (millimetres
// per model unit, overriding the header), filter (face record → boolean).
export function satToStep(model, opts = {}) {
  const A = access(model), M = A.get, L = [], E = s => '#' + L.push(s), memo = new Map(), scale = opts.scale > 0 ? opts.scale : model.scale || 1;
  const stats = { bodies: 0, solids: 0, sheets: 0, wireBodies: 0, faces: 0, translated: 0, skipped: {}, edges: 0, exactEdges: 0, approximatedEdges: 0, approximated: {}, approximatedFaces: 0, approximatedSurfaces: {}, edgeCurves: {}, surfaces: {}, coneVotes: [0, 0] }, bodies = [], notes = [], lost = [];
  const skip = (why, n = 1) => { stats.skipped[why] = (stats.skipped[why] || 0) + n; };
  let tf = null, tfKey = '';
  // ACIS transforms act on row vectors: p' = (p · R) s + t
  const P = p => { if (!isVec(p)) throw new Error('unset point'); if (!tf) return [p[0] * scale, p[1] * scale, p[2] * scale]; const R = tf.R, s = tf.s * scale; return [(p[0] * R[0] + p[1] * R[3] + p[2] * R[6]) * s + tf.t[0] * scale, (p[0] * R[1] + p[1] * R[4] + p[2] * R[7]) * s + tf.t[1] * scale, (p[0] * R[2] + p[1] * R[5] + p[2] * R[8]) * s + tf.t[2] * scale]; };
  const D = d => { if (!isVec(d)) throw new Error('unset direction'); if (!tf) return d; const R = tf.R; return [d[0] * R[0] + d[1] * R[3] + d[2] * R[6], d[0] * R[1] + d[1] * R[4] + d[2] * R[7], d[0] * R[2] + d[1] * R[5] + d[2] * R[8]]; };
  const S = r => (tf ? r * tf.s : r) * scale;
  const pt = p => E(`CARTESIAN_POINT('',(${real(p[0])},${real(p[1])},${real(p[2])}))`);
  const dir = d => { const l = len(d); if (!(l > 1e-12)) throw new Error('zero direction'); return E(`DIRECTION('',(${real(d[0] / l)},${real(d[1] / l)},${real(d[2] / l)}))`); };
  const perp = z => { const a = abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return unit(cross(cross(z, a), z)); };
  const ax2 = (o, z, x) => { const zz = unit(D(z)); let xx = x ? D(x) : perp(zz); xx = sub(xx, mul(zz, dot(xx, zz))); if (!(len(xx) > 1e-9 * max(1, len(x || xx)))) xx = perp(zz); return E(`AXIS2_PLACEMENT_3D('',${pt(P(o))},${dir(zz)},${dir(xx)})`); };
  const list = a => '(' + a.join(',') + ')';
  const bspCurve = (deg, pts, knots, mult, weights, closed) => {
    const ids = list(pts.map(pt)), k = `${list(mult)},${list(knots.map(real))},.UNSPECIFIED.`;
    return weights ? E(`(BOUNDED_CURVE()B_SPLINE_CURVE(${deg},${ids},.UNSPECIFIED.,${B(closed)},.F.)B_SPLINE_CURVE_WITH_KNOTS(${k})CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE(${list(weights.map(real))})REPRESENTATION_ITEM(''))`) : E(`B_SPLINE_CURVE_WITH_KNOTS('',${deg},${ids},.UNSPECIFIED.,${B(closed)},.F.,${k})`);
  };
  // homogeneous form of a decoded B-spline curve, as the knot-insertion kit wants it
  const homog = q => { const Pn = []; for (let i = 0; i < q.n; i++) { const w = q.W ? q.W[i] : 1; Pn.push(q.P[i][0] * w, q.P[i][1] * w, q.P[i][2] * w, w); } const U = flatKnots(q.knots, q.mult); return { degree: q.degree, dim: 4, rational: true, n: q.n, P: Pn, U, knots: q.knots, mult: q.mult, closed: q.closure !== 0 }; };
  const writeBs = (h, wasRational) => { const pts = [], w = []; for (let i = 0; i < h.n; i++) { const ww = h.P[4 * i + 3]; pts.push(P([h.P[4 * i] / ww, h.P[4 * i + 1] / ww, h.P[4 * i + 2] / ww])); w.push(ww); } return bspCurve(h.degree, pts, h.knots, h.mult, wasRational ? w : null, false); };
  const clamped = q => q.mult[0] === q.degree + 1 && q.mult[q.mult.length - 1] === q.degree + 1;
  // exact curves → STEP id (memoised per transform)
  const curve3 = (c) => {
    const key = 'c' + c.i + tfKey; if (memo.has(key)) return memo.get(key); let id = null;
    if (c.kind === 'line') id = E(`LINE('',${pt(P(c.p))},${E(`VECTOR('',${dir(D(c.d))},1.)`)})`);
    else if (c.kind === 'ellipse') {
      const a = len(c.m), r = abs(c.ratio); if (!(a > 0 && r > 0)) throw new Error('degenerate ellipse');
      if (abs(r - 1) <= 1e-12) id = E(`CIRCLE('',${ax2(c.c, c.n, c.m)},${real(S(a))})`);
      else if (r < 1) id = E(`ELLIPSE('',${ax2(c.c, c.n, c.m)},${real(S(a))},${real(S(a * r))})`);
      else id = E(`ELLIPSE('',${ax2(c.c, c.n, cross(c.n, c.m))},${real(S(a * r))},${real(S(a))})`);
    } else if (c.kind === 'bspline') { const q = c.q; if (clamped(q)) id = bspCurve(q.degree, q.P.map(P), q.knots, q.mult, q.W, false); else { const h = homog(q), sc = subCurve({ ...h, closed: false }, h.U[h.degree], h.U[h.n]); if (!sc) throw new Error('unclamped B-spline curve could not be clamped'); id = writeBs(sc, q.rational); } }
    memo.set(key, id); return id;
  };
  // the radius of a cone changes along its axis by ± |sin / cos| per unit length; the sign is taken from the face's own
  // vertices when they decide it and from the save-file rule otherwise (radius shrinks along the axis when sin · cos < 0)
  const coneSlope = (s, fa) => {
    const k = abs(s.sin / s.cos), R = len(s.m), ax = unit(s.n); let up = 0, down = 0;
    if (fa) for (const p of facePoints(fa)) { const d = sub(p, s.c), h = dot(d, ax), rho = len(sub(d, mul(ax, h))); if (abs(h) < 1e-7 * max(1, R)) continue; const eu = abs(rho - (R + h * k)), ed = abs(rho - (R - h * k)); if (eu < ed * 0.5) up++; else if (ed < eu * 0.5) down++; }
    const rule = s.sin * s.cos < 0 ? -1 : 1; if (up || down) { stats.coneVotes[(up > down ? 1 : -1) === rule ? 0 : 1]++; return up > down ? k : -k; }
    return rule * k;
  };
  // surfaces → {id, flip} (flip: the STEP surface normal is opposite to the ACIS natural normal) or {why}
  const surface = (s, fa) => {
    const perFace = s.kind === 'cone' && abs(s.sin) > 1e-12, key = 'S' + s.i + tfKey + (perFace && fa ? '/' + fa.i : ''); if (memo.has(key)) return memo.get(key); let r;
    try {
      if (s.kind === 'plane') r = { id: E(`PLANE('',${ax2(s.p, s.n, s.u)})`), flip: false };
      else if (s.kind === 'cone') {
        if (abs(s.ratio - 1) > 1e-9) r = { why: 'elliptical cone or cylinder' };
        else if (abs(s.sin) <= 1e-12) r = { id: E(`CYLINDRICAL_SURFACE('',${ax2(s.c, s.n, s.m)},${real(S(len(s.m)))})`), flip: s.cos < 0 };
        else {
          const k = coneSlope(s, fa), ang = atan2(abs(s.sin), abs(s.cos)); if (!(ang > 1e-9 && ang < PI / 2 - 1e-9)) throw new Error('degenerate cone angle');
          // STEP's cone widens along its axis: a cone that narrows along the saved axis is written with the axis reversed
          r = { id: E(`CONICAL_SURFACE('',${ax2(s.c, k > 0 ? s.n : mul(s.n, -1), s.m)},${real(S(len(s.m)))},${real(ang)})`), flip: s.cos < 0, slope: k };
        }
      } else if (s.kind === 'sphere') r = { id: E(`SPHERICAL_SURFACE('',${ax2(s.c, s.pole || [0, 0, 1], s.u)},${real(S(abs(s.r)))})`), flip: s.r < 0 };
      else if (s.kind === 'torus') { const a = s.R, b = abs(s.r), pl = ax2(s.c, s.n, s.u); r = a > b ? { id: E(`TOROIDAL_SURFACE('',${pl},${real(S(a))},${real(S(b))})`), flip: s.r < 0 } : { id: E(`DEGENERATE_TOROIDAL_SURFACE('',${pl},${real(S(abs(a)))},${real(S(b))},${B(a > 0)})`), flip: s.r < 0 }; }
      else if (s.kind === 'bspline') {
        const q = s.q, rows = [], wr = []; if (q.mu[0] !== q.du + 1 || q.mu[q.mu.length - 1] !== q.du + 1 || q.mv[0] !== q.dv + 1 || q.mv[q.mv.length - 1] !== q.dv + 1) throw new Error('unclamped (periodic) knot vector');
        for (let i = 0; i < q.nu; i++) { rows.push(list(q.P[i].map(p => pt(P(p))))); if (q.W) wr.push(list(q.W[i].map(real))); }
        const k = `${list(q.mu)},${list(q.mv)},${list(q.ku.map(real))},${list(q.kv.map(real))},.UNSPECIFIED.`, hd = `${q.du},${q.dv},${list(rows)},.UNSPECIFIED.,.F.,.F.,.F.`;
        r = { id: q.W ? E(`(BOUNDED_SURFACE()B_SPLINE_SURFACE(${hd})B_SPLINE_SURFACE_WITH_KNOTS(${k})GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(${list(wr)})REPRESENTATION_ITEM('')SURFACE())`) : E(`B_SPLINE_SURFACE_WITH_KNOTS('',${hd},${k})`), flip: s.rev, approx: !s.exact, sub: s.sub };
      } else if (s.kind === 'extrusion') {
        const pc = s.profile; pc.i = 'x' + s.i; const id = curve3(pc);
        // natural normal c1' × c2': the STEP extrusion's is profile' × direction
        r = { id: E(`SURFACE_OF_LINEAR_EXTRUSION('',${id},${E(`VECTOR('',${dir(D(s.d))},1.)`)})`), flip: s.rev !== (s.lineFirst !== (pc.kind === 'bspline' && pc.rev)), sub: s.sub };
      } else r = { why: s.why || `${s.name} surface` };
    } catch (e) { r = { why: `${s.name} (${e.message})` }; }
    if (r.id) { const kd = s.kind === 'bspline' ? `spline (${s.sub})` : s.kind === 'extrusion' ? `extrusion (${s.sub})` : s.kind === 'cone' && abs(s.sin) <= 1e-12 ? 'cylinder' : s.kind; stats.surfaces[kd] = (stats.surfaces[kd] || 0) + 1; }
    memo.set(key, r); return r;
  };
  const vpos = v => { const p = v && M(v.point); return p && isVec(p.p) ? p.p : null; };
  const vertex = v => { const key = 'v' + v.i + tfKey; if (memo.has(key)) return memo.get(key); const p = vpos(v); if (!p) throw new Error('vertex without a point'); const id = E(`VERTEX_POINT('',${pt(P(p))})`); memo.set(key, id); return id; };
  const ring = lp => { const out = [], c0 = M(lp.coedge); for (let c = c0, g = 0; c && c.base === 'coedge' && g < 1e5; g++) { out.push(c); c = M(c.next); if (c === c0) return out; } return out.length === 1 || !c0 ? out : null; };
  const loopsOf = fa => { const out = []; for (let lp = M(fa.loop), g = 0; lp && lp.base === 'loop' && g < 1e5; lp = M(lp.next), g++) out.push(lp); return out; };
  function facePoints(fa) { const out = []; for (const lp of loopsOf(fa)) for (const c of ring(lp) || []) { const e = M(c.edge); if (!e || e.base !== 'edge') continue; for (const v of [M(e.start), M(e.end)]) { const p = vpos(v); if (p) out.push(p); } const cu = M(e.curve); if (cu && cu.kind === 'ellipse' && isVec(cu.m) && isVec(cu.c) && e.t0 !== null && e.t1 !== null) { const t = (e.rev ? -1 : 1) * (e.t0 + e.t1) / 2; out.push(add(cu.c, add(mul(cu.m, cos(t)), mul(cross(unit(cu.n), cu.m), cu.ratio * sin(t))))); } } return out; }
  // edge → {id} | {degenerate: vertex} | {why}
  const edge = e => {
    const key = 'e' + e.i + tfKey; if (memo.has(key)) return memo.get(key); let r; stats.edges++;
    try {
      const vs = M(e.start), ve = M(e.end), c = M(e.curve);
      if (!vs || !ve || vs.base !== 'vertex' || ve.base !== 'vertex') r = { why: 'edge without vertices' };
      else if (!c) r = { degenerate: vs };
      else if (c.base !== 'curve' || c.kind === 'none') r = { why: c.why || `edge on a ${c.name}` };
      else {
        const va = vertex(vs), vb = ve === vs ? va : vertex(ve); let id = null, same = !e.rev; stats.edgeCurves[c.kind === 'bspline' ? `intcurve (${c.sub})` : c.kind] = (stats.edgeCurves[c.kind === 'bspline' ? `intcurve (${c.sub})` : c.kind] || 0) + 1;
        if (c.kind === 'line' && vs === ve) r = { why: 'unbounded straight edge' };
        else {
          if (c.kind === 'bspline') {
            same = e.rev === c.rev; const q = c.q, h = homog(q), a = h.U[h.degree], b = h.U[h.n], T = b - a;
            // the stretch of the B-spline that the edge covers, in the spline's own parameter (edge and curve senses undone)
            let t0 = e.t0, t1 = e.t1; if (t0 !== null && t1 !== null) { if (e.rev) [t0, t1] = [-t1, -t0]; if (c.rev) [t0, t1] = [-t1, -t0]; }
            const whole = t0 === null || t1 === null || (abs(t0 - a) <= 1e-7 * T && abs(t1 - b) <= 1e-7 * T);
            if (!whole && t1 > t0 + 1e-9 * T && (q.closure !== 0 || (t0 >= a - 1e-7 * T && t1 <= b + 1e-7 * T))) { const sc = subCurve(h, max(t0, q.closure !== 0 ? -Infinity : a), q.closure !== 0 ? t1 : min(t1, b)); if (sc) { id = writeBs(sc, q.rational); stats.trimmedBCurves = (stats.trimmedBCurves || 0) + 1; } }
          }
          if (!id) id = curve3(c);
          // the saved edge parameters must reproduce the end vertices on the decoded curve (a check of the decoding, not a repair)
          if (e.t0 !== null && e.t1 !== null && c.kind !== 'bspline') { const sg = e.rev ? -1 : 1, at = t => c.kind === 'line' ? add(c.p, mul(c.d, t)) : add(c.c, add(mul(c.m, cos(t)), mul(cross(unit(c.n), c.m), c.ratio * sin(t)))), a = vpos(vs), b = vpos(ve), ext = max(len(sub(a, b)), c.kind === 'ellipse' ? len(c.m) : 0, 1e-9); if (len(sub(at(sg * e.t0), a)) > 1e-3 * ext + 1e3 * (model.resabs || 1e-6) || len(sub(at(sg * e.t1), b)) > 1e-3 * ext + 1e3 * (model.resabs || 1e-6)) stats.paramMismatch = (stats.paramMismatch || 0) + 1; }
          r = { id: E(`EDGE_CURVE('',${same ? va : vb},${same ? vb : va},${id},.T.)`), flipped: !same };
          if (c.kind === 'bspline' && !c.exact) { stats.approximatedEdges++; const kd = `intcurve (${c.sub}) through its stored approximating B-spline`; stats.approximated[kd] = (stats.approximated[kd] || 0) + 1; } else stats.exactEdges++;
        }
      }
    } catch (err) { r = { why: `edge (${err.message})` }; }
    memo.set(key, r); return r;
  };
  const loop = lp => {
    const cs = ring(lp); if (!cs || !cs.length) return { why: 'loop that does not close' };
    const oe = []; let apex = null;
    for (const c of cs) { const en = M(c.edge); if (!en || en.base !== 'edge') { return { why: 'coedge without an edge' }; } const e = edge(en); if (e.degenerate) { apex = e.degenerate; continue; } if (!e.id) return e; oe.push(E(`ORIENTED_EDGE('',*,*,${e.id},${B(c.rev === !!e.flipped)})`)); }
    if (!oe.length) return apex ? { id: E(`VERTEX_LOOP('',${vertex(apex)})`) } : { why: 'empty loop' };
    return { id: E(`EDGE_LOOP('',${list(oe)})`) };
  };
  // distance of a point from a decoded analytic surface (model units), used to check each face's boundary against its surface
  const offSurface = (s, p, slope) => {
    if (s.kind === 'plane') return abs(dot(sub(p, s.p), unit(s.n)));
    if (s.kind === 'sphere') return abs(len(sub(p, s.c)) - abs(s.r));
    const ax = unit(s.n), d = sub(p, s.c), h = dot(d, ax), rho = len(sub(d, mul(ax, h)));
    if (s.kind === 'torus') return abs(hypot(rho - s.R, h) - abs(s.r));
    if (s.kind === 'cone') { if (abs(s.ratio - 1) > 1e-9) return 0; const k = abs(s.sin) > 1e-12 ? slope : 0; return abs(rho - abs(len(s.m) + h * k)) / sqrt(1 + k * k); }
    return 0;
  };
  const face = fa => {
    if (opts.filter && !opts.filter(fa)) return null;
    stats.faces++; const s = M(fa.surface); if (!s || s.base !== 'surface') { skip('face without a surface'); return null; }
    try {
      const sr = surface(s, fa); if (!sr.id) { skip(sr.why); return null; }
      const bounds = []; for (const lp of loopsOf(fa)) { const r = loop(lp); if (!r.id) { skip(`face bounded by: ${r.why}`); return null; } bounds.push(r.id); }
      if (!bounds.length) {
        if (s.kind === 'sphere') bounds.push(E(`VERTEX_LOOP('',${E(`VERTEX_POINT('',${pt(P(sub(s.c, mul(unit(s.pole || [0, 0, 1]), abs(s.r)))))})`)})`));
        else if (s.kind === 'torus') bounds.push(E(`VERTEX_LOOP('',${E(`VERTEX_POINT('',${pt(P(add(s.c, mul(unit(s.u || perp(unit(s.n))), s.R + abs(s.r)))))})`)})`));
        else { skip(`unbounded face on a ${s.kind} surface`); return null; }
      }
      if (s.kind !== 'bspline' && s.kind !== 'extrusion') { const pts = facePoints(fa); let dmax = 0, ext = 0; const k = s.kind === 'cone' && abs(s.sin) > 1e-12 ? sr.slope : 0; for (const p of pts) { dmax = max(dmax, offSurface(s, p, k)); ext = max(ext, len(sub(p, pts[0]))); } if (dmax > max(1e-4 * ext, 1e3 * (model.resabs || 1e-6))) { stats.offSurface = (stats.offSurface || 0) + 1; stats.offSurfaceKinds = stats.offSurfaceKinds || {}; stats.offSurfaceKinds[s.kind] = max(stats.offSurfaceKinds[s.kind] || 0, dmax * scale); } }
      // edge orientation was written explicitly (EDGE_CURVE always runs with its curve), so ORIENTED_EDGE carries coedge × edge sense
      const same = fa.rev === sr.flip, id = E(`ADVANCED_FACE('',${list(bounds.map(b => E(`${bounds.length === 1 ? 'FACE_OUTER_BOUND' : 'FACE_BOUND'}('',${b},.T.)`)))},${sr.id},${B(same)})`);
      if (sr.approx) { stats.approximatedFaces++; stats.approximatedSurfaces[sr.sub] = (stats.approximatedSurfaces[sr.sub] || 0) + 1; }
      stats.translated++; return id;
    } catch (e) { skip(`face (${e.message})`); return null; }
  };
  const facesOf = sh => {
    const out = [], seen = new Set(), chain = f0 => { for (let f = M(f0), g = 0; f && f.base === 'face' && g < 1e6; f = M(f.next), g++) { if (seen.has(f.i)) break; seen.add(f.i); out.push(f); } };
    chain(sh.face); const walk = (i, depth) => { for (let s = M(i), g = 0; s && s.base === 'subshell' && g < 1e5 && depth < 64; s = M(s.next), g++) { chain(s.face); walk(s.child, depth + 1); } }; walk(sh.subshell, 0);
    return out;
  };
  const closedTopo = fs => { for (const f of fs) { if (f.double) return false; for (const lp of loopsOf(f)) for (const c of ring(lp) || []) { const p = M(c.partner); if (!p || p.base !== 'coedge') { const e = M(c.edge); if (e && e.curve >= 0) return false; } } } return true; };
  const extent = fs => { let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (const f of fs) { const s = M(f.surface); const grow = (p, r = 0) => { for (let a = 0; a < 3; a++) { lo[a] = min(lo[a], p[a] - r); hi[a] = max(hi[a], p[a] + r); } }; if (s && s.kind === 'sphere') grow(s.c, abs(s.r)); if (s && s.kind === 'torus') grow(s.c, abs(s.R) + abs(s.r)); for (const p of facePoints(f)) grow(p); } return hi[0] >= lo[0] ? hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) : 0; };
  // shared context
  const ctxApp = E("APPLICATION_CONTEXT('core data for automotive mechanical design processes')"); E(`APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,${ctxApp})`);
  const ctxProd = E(`PRODUCT_CONTEXT('',${ctxApp},'mechanical')`), ctxDef = E(`PRODUCT_DEFINITION_CONTEXT('part definition',${ctxApp},'design')`);
  const uLen = E('(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))'), uAng = E('(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))'), uSol = E('(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())');
  const res = model.resabs > 0 && model.resabs * scale < 0.1 ? model.resabs * scale : 1e-5;
  const unc = E(`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(${real(res)}),${uLen},'distance_accuracy_value','confusion accuracy')`);
  const ctx = E(`(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((${unc}))GLOBAL_UNIT_ASSIGNED_CONTEXT((${uLen},${uAng},${uSol}))REPRESENTATION_CONTEXT('Context #1','3D'))`);
  const origin = () => E(`AXIS2_PLACEMENT_3D('',${pt([0, 0, 0])},${dir([0, 0, 1])},${dir([1, 0, 0])})`);
  const stem = opts.name ? String(opts.name).replace(/\.[^.]+$/, '') : 'body';
  const nameOf = attr => { let nm = ''; for (const a of A.attributes(attr)) { if (!/name_attrib|gen-attrib/.test(a.name)) continue; const ss = a.t.filter(x => x instanceof Str).map(x => utf8(x.s)); if (ss.length >= 2 && /NAME$/i.test(ss[0]) && ss[1].trim()) return ss[1].trim(); if (ss.length >= 2 && /^(name|ATTRIB_XACIS_ID)$/i.test(ss[0]) && ss[1].trim() && !nm) nm = ss[1].trim(); } return nm; };
  const colourOf = attr => { for (const a of A.attributes(attr)) { if (/^rgb_color-/.test(a.name) || /^(colour|color)-/.test(a.name)) { const n = a.t.filter(x => typeof x === 'number'); if (/^rgb/.test(a.name) && n.length >= 3) { const c = n.slice(-3); if (c.every(x => x >= 0 && x <= 1)) return c; } } if (/^adesk-(truecolor|color)-/.test(a.name) || /truecolor/.test(a.name)) { const n = a.t.filter(x => typeof x === 'number'); if (n.length) { const u = n[n.length - 1] >>> 0; return [((u >>> 16) & 255) / 255, ((u >>> 8) & 255) / 255, (u & 255) / 255]; } } } return null; };
  const body = (b, label) => {
    const name = nameOf(b.attr) || label, rec = { index: b.i, name, type: 'sheet', faces: 0, translated: 0, closed: false, colour: null }; stats.bodies++;
    const t = M(b.transform); tf = null; tfKey = '';
    if (t && t.base === 'transform' && !t.bad) { const idn = t.R.every((x, i) => abs(x - (i % 4 === 0 ? 1 : 0)) < 1e-14) && t.t.every(x => x === 0) && abs(t.s - 1) < 1e-14; if (!idn) { const R = t.R, det = R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]); if (t.reflect || det < 0) { lost.push(`body "${name}" carries a mirroring transform, which cannot be applied to analytic surface placements: it was left out`); bodies.push(rec); return; } tf = { R, t: t.t, s: t.s > 0 ? t.s : 1 }; tfKey = '@' + b.i; stats.transformed = (stats.transformed || 0) + 1; } }
    const items = [], open = [], f0 = stats.faces, t0 = stats.translated, colours = new Map(); let complete = true, any = false, lumps = 0, wires = 0;
    if (M(b.wire)) wires++;
    for (let lu = M(b.lump), g = 0; lu && lu.base === 'lump' && g < 1e5; lu = M(lu.next), g++) {
      lumps++; const sh = []; for (let s = M(lu.shell), h = 0; s && s.base === 'shell' && h < 1e5; s = M(s.next), h++) { const fs = facesOf(s); if (fs.length) sh.push(fs); else if (M(s.wire)) wires++; }
      if (!sh.length) continue;
      if (sh.length > 1) { const ex = sh.map(extent); let o = 0; ex.forEach((x, i) => { if (x > ex[o]) o = i; }); if (o) sh.unshift(sh.splice(o, 1)[0]); }
      const done = sh.map(fs => { for (const f of fs) { const c = colourOf(f.attr); if (c) { const k = c.join(','); colours.set(k, (colours.get(k) || 0) + 1); } } const ids = fs.map(face).filter(Boolean); return { ids, all: ids.length === fs.length && closedTopo(fs) }; });
      if (done.every(d => d.all)) { const cs = done.map(d => E(`CLOSED_SHELL('',${list(d.ids)})`)); items.push(cs.length === 1 ? E(`MANIFOLD_SOLID_BREP(${str(name)},${cs[0]})`) : E(`BREP_WITH_VOIDS(${str(name)},${cs[0]},${list(cs.slice(1).map(c => E(`ORIENTED_CLOSED_SHELL('',*,${c},.F.)`)))})`)); any = true; }
      else { complete = false; for (const d of done) if (d.ids.length) open.push(E(`OPEN_SHELL('',${list(d.ids)})`)); }
    }
    tf = null; tfKey = '';
    rec.faces = stats.faces - f0; rec.translated = stats.translated - t0; rec.closed = complete && any; rec.type = rec.closed ? 'solid' : rec.faces ? 'sheet' : 'wire';
    let best = colourOf(b.attr); if (!best) { let n = 0; for (const [k, c] of colours) if (c > n) { n = c; best = k.split(',').map(Number); } } rec.colour = best;
    bodies.push(rec); if (!rec.faces) { stats.wireBodies++; return; }
    if (!items.length && !open.length) return;
    if (rec.closed) stats.solids++; else stats.sheets++;
    const prod = E(`PRODUCT(${str(name)},${str(name)},'',(${ctxProd}))`), pdf = E(`PRODUCT_DEFINITION_FORMATION('','',${prod})`), pd = E(`PRODUCT_DEFINITION('design','',${pdf},${ctxDef})`), pds = E(`PRODUCT_DEFINITION_SHAPE('','',${pd})`);
    E(`PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(${prod}))`);
    const all = [...items]; if (open.length) all.push(E(`SHELL_BASED_SURFACE_MODEL(${str(name)},${list(open)})`));
    const rep = E(`${open.length ? 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION' : 'ADVANCED_BREP_SHAPE_REPRESENTATION'}('',${list([origin(), ...all])},${ctx})`); E(`SHAPE_DEFINITION_REPRESENTATION(${pds},${rep})`);
    rec.written = true;
  };
  const want = opts.bodies ? new Set(opts.bodies) : null; let nb = 0;
  for (const e of model.entities.values()) { if (e.base !== 'body') continue; if (want && !want.has(e.i)) continue; const b = M(e.i); if (!b || b.bad) { lost.push(`a body record could not be decoded (${b ? b.bad : 'missing'})`); continue; } nb++; body(b, `${stem}${nb > 1 ? ' ' + nb : ''}`); }
  if (model.scale === null) notes.push('the file predates the units record: one model unit was taken as one millimetre');
  const written = bodies.some(b => b.written);
  const headTxt = `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('ACIS B-rep translated entity by entity'),'2;1');\nFILE_NAME(${str(stem)},'',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n`;
  const step = written ? headTxt + L.map((s, i) => `#${i + 1}=${s};`).join('\n') + '\nENDSEC;\nEND-ISO-10303-21;\n' : null;
  return { step, bodies, stats, notes, lost };
}
