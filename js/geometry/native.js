// Native / proprietary CAD containers: everything that can be read from the container structure without the vendor's
// modelling kernel — OLE2 compound files (property sets, thumbnails, stream inventory), the SOLIDWORKS 2015+ chunk
// container (preview, document properties, configuration names, embedded Parasolid partitions), Siemens NX parts (their
// Parasolid partitions), Inventor, Fusion and FreeCAD archives, the JT table of contents, the Parasolid text header and
// ACIS SAT topology points. Nothing here evaluates a B-rep: each reader states exactly what it recovered, and the
// Parasolid and ShapeManager streams it unpacks are translated elsewhere (xt.js, sat.js).
import * as F from './formats.js';
import { findACIS, readSAT } from './sat.js';

const { min, max, abs, hypot } = Math;
const latin = (b, o, n) => { let s = ''; for (let i = o, e = min(b.length, o + n); i < e; i++) s += String.fromCharCode(b[i]); return s; };
const dv = b => new DataView(b.buffer, b.byteOffset, b.byteLength);
const startsWith = (b, sig, o = 0) => { if (b.length < o + sig.length) return false; for (let i = 0; i < sig.length; i++) if (b[o + i] !== sig[i]) return false; return true; };
const PNG = [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A], JPG = [0xFF, 0xD8, 0xFF];
const utf16 = (b, o, n) => new TextDecoder('utf-16le').decode(b.subarray(o, min(b.length, o + n)));
const xmlText = s => s.replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&#(\d+);/g, (m, d) => String.fromCharCode(+d)).replace(/&amp;/g, '&');

// ---------------------------------------------------------------- images
// A packed DIB (BITMAPINFOHEADER + palette + pixels), as stored in OLE thumbnails, wrapped into a .bmp file.
export function dibToBmp(dib) {
  if (dib.length < 40) return null; const v = dv(dib), hs = v.getUint32(0, true); if (hs < 12 || hs > 124 || hs > dib.length) return null;
  const bpp = hs === 12 ? v.getUint16(10, true) : v.getUint16(14, true), comp = hs >= 20 ? v.getUint32(16, true) : 0; let colours = hs >= 36 ? v.getUint32(32, true) : 0; if (!colours && bpp <= 8) colours = 1 << bpp;
  const off = 14 + hs + colours * (hs === 12 ? 3 : 4) + (comp === 3 && hs === 40 ? 12 : 0), out = new Uint8Array(14 + dib.length), o = dv(out);
  out[0] = 0x42; out[1] = 0x4D; o.setUint32(2, out.length, true); o.setUint32(10, min(off, out.length), true); out.set(dib, 14); return out;
}
// Recognise an image at (or shortly after) the start of a stream: PNG, JPEG, BMP or a bare DIB.
export function sniffImage(b, scan = 64) {
  for (let o = 0; o <= min(scan, max(0, b.length - 8)); o++) {
    if (startsWith(b, PNG, o)) return { mime: 'image/png', bytes: b.subarray(o) }; if (startsWith(b, JPG, o)) return { mime: 'image/jpeg', bytes: b.subarray(o) };
    if (b[o] === 0x42 && b[o + 1] === 0x4D && b.length - o > 54 && dv(b).getUint32(o + 14, true) === 40) return { mime: 'image/bmp', bytes: b.subarray(o) };
    if (b.length - o > 40 && dv(b).getUint32(o, true) === 40 && dv(b).getUint16(o + 12, true) === 1 && [1, 4, 8, 16, 24, 32].includes(dv(b).getUint16(o + 14, true))) { const bmp = dibToBmp(b.subarray(o)); if (bmp) return { mime: 'image/bmp', bytes: bmp }; }
  }
  return null;
}
function imageSize(img) { try { const b = img.bytes, v = dv(b); if (img.mime === 'image/png' && b.length > 24) return [v.getUint32(16, false), v.getUint32(20, false)]; if (img.mime === 'image/bmp' && b.length > 26) return [abs(v.getInt32(18, true)), abs(v.getInt32(22, true))]; } catch (e) { /* unknown */ } return null; }

// ---------------------------------------------------------------- OLE2 compound file (MS-CFB)
// cfbOpen(bytes) → {version, sectorSize, entries: [{name, path, type: 'storage' | 'stream' | 'root', size, clsid}], read(entry | path) → Uint8Array, find(regexp)}
export function cfbOpen(bytes) {
  const b = F.toBytes(bytes), v = dv(b); if (b.length < 512 || !startsWith(b, [0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1])) throw new Error('CFB: missing OLE2 compound-file signature');
  const major = v.getUint16(26, true), ss = 1 << v.getUint16(30, true), ms = 1 << v.getUint16(32, true), nFat = v.getUint32(44, true), dirStart = v.getUint32(48, true), cutoff = v.getUint32(56, true), miniFatStart = v.getUint32(60, true), difStart = v.getUint32(68, true), nDif = v.getUint32(72, true);
  if (ss !== 512 && ss !== 4096) throw new Error(`CFB: unsupported sector size ${ss}`); if (ms !== 64) throw new Error(`CFB: unsupported mini-sector size ${ms}`);
  const END = 0xFFFFFFFE, FREE = 0xFFFFFFFF, off = s => (s + 1) * ss, maxSec = Math.floor(b.length / ss);
  const fatSecs = []; for (let i = 0; i < 109 && fatSecs.length < nFat; i++) { const s = v.getUint32(76 + 4 * i, true); if (s < END) fatSecs.push(s); }
  for (let s = difStart, g = 0; s < END && g < nDif + 2 && fatSecs.length < nFat; g++) { const o = off(s); if (o + ss > b.length) break; for (let i = 0; i < ss / 4 - 1 && fatSecs.length < nFat; i++) { const t = v.getUint32(o + 4 * i, true); if (t < END) fatSecs.push(t); } s = v.getUint32(o + ss - 4, true); }
  const per = ss / 4, fat = new Uint32Array(fatSecs.length * per).fill(FREE); fatSecs.forEach((s, k) => { const o = off(s); if (o + ss <= b.length) for (let i = 0; i < per; i++) fat[k * per + i] = v.getUint32(o + 4 * i, true); });
  const chain = (start, limitBytes) => { const out = []; for (let s = start, g = 0; s < END && s < fat.length && g <= maxSec; g++) { out.push(s); if (limitBytes !== undefined && out.length * ss >= limitBytes) break; s = fat[s]; } return out; };
  const readChain = (start, size) => { const secs = chain(start, size), out = new Uint8Array(size === undefined ? secs.length * ss : size); let p = 0; for (const s of secs) { const o = off(s), n = min(ss, out.length - p, max(0, b.length - o)); if (n <= 0) break; out.set(b.subarray(o, o + n), p); p += n; } return out; };
  const dir = readChain(dirStart), dvd = dv(dir), raw = [];
  for (let o = 0; o + 128 <= dir.length; o += 128) { const nl = dvd.getUint16(o + 64, true), type = dir[o + 66]; raw.push({ name: nl >= 2 && nl <= 64 ? utf16(dir, o, nl - 2) : '', type, left: dvd.getUint32(o + 68, true), right: dvd.getUint32(o + 72, true), child: dvd.getUint32(o + 76, true), clsid: Array.from(dir.subarray(o + 80, o + 96)).map(x => x.toString(16).padStart(2, '0')).join(''), start: dvd.getUint32(o + 116, true), size: major >= 4 ? dvd.getUint32(o + 120, true) + dvd.getUint32(o + 124, true) * 4294967296 : dvd.getUint32(o + 120, true) }); }
  if (!raw.length || raw[0].type !== 5) throw new Error('CFB: root directory entry not found');
  const mini = readChain(raw[0].start, min(raw[0].size, b.length)), miniFatBytes = miniFatStart < END ? readChain(miniFatStart) : new Uint8Array(0), mfv = dv(miniFatBytes), miniFat = i => 4 * i + 4 <= miniFatBytes.length ? mfv.getUint32(4 * i, true) : END;
  const entries = [], seen = new Set();
  const walk = (id, path) => { const stack = [id]; while (stack.length) { const i = stack.pop(); if (i >= raw.length || seen.has(i)) continue; seen.add(i); const e = raw[i]; if (e.type !== 1 && e.type !== 2) continue; const p = path ? path + '/' + e.name : e.name; entries.push({ name: e.name, path: p, type: e.type === 1 ? 'storage' : 'stream', size: e.type === 2 ? e.size : 0, clsid: e.clsid, _i: i }); if (e.left < raw.length) stack.push(e.left); if (e.right < raw.length) stack.push(e.right); if (e.type === 1 && e.child < raw.length) walk(e.child, p); } };
  seen.add(0); if (raw[0].child < raw.length) walk(raw[0].child, '');
  const read = x => { const e = typeof x === 'string' ? entries.find(q => q.path === x) : x; if (!e || e.type !== 'stream') return null; const r = raw[e._i], size = min(r.size, b.length);
    if (r.size < cutoff) { const out = new Uint8Array(size); let p = 0; for (let s = r.start, g = 0; s < END && p < size && g < 1e6; g++) { const o = s * 64, n = min(64, size - p); if (o + n > mini.length) break; out.set(mini.subarray(o, o + n), p); p += n; s = miniFat(s); } return out; }
    return readChain(r.start, size); };
  return { version: major, sectorSize: ss, rootClsid: raw[0].clsid, entries, read, find: re => entries.filter(e => e.type === 'stream' && re.test(e.path)) };
}

// ---------------------------------------------------------------- OLE property sets (MS-OLEPS)
const FMTID = { 'e0859ff2f94f6810ab9108002b27b3d9': 'SummaryInformation', '02d5cdd59c2e1b10939708002b2cf9ae': 'DocumentSummaryInformation', '05d5cdd59c2e1b10939708002b2cf9ae': 'UserDefinedProperties', '0f3f8532443411d19e930060b03c1ca6': 'Inventor Design Tracking Properties' };
const PID = { SummaryInformation: { 2: 'title', 3: 'subject', 4: 'author', 5: 'keywords', 6: 'comments', 7: 'template', 8: 'lastSavedBy', 9: 'revisionNumber', 12: 'created', 13: 'lastSaved', 18: 'application' }, DocumentSummaryInformation: { 2: 'category', 14: 'manager', 15: 'company' }, 'Inventor Design Tracking Properties': { 4: 'creationDate', 5: 'partNumber', 7: 'project', 20: 'material', 29: 'description', 41: 'designer', 42: 'engineer' } };
const fileTime = (lo, hi) => { const t = (hi * 4294967296 + lo) / 1e4 - 11644473600000; return t > 0 && t < 4e12 ? new Date(t).toISOString() : null; };
// parsePropertySets(stream) → [{name, fmtid, properties: {key: value}, thumbnail?: {mime, bytes}}]
export function parsePropertySets(bytes) {
  const b = F.toBytes(bytes), out = []; if (b.length < 48) return out; const v = dv(b); if (v.getUint16(0, true) !== 0xFFFE) return out;
  const n = min(2, v.getUint32(24, true));
  for (let s = 0; s < n; s++) {
    const fid = Array.from(b.subarray(28 + 20 * s, 44 + 20 * s)).map(x => x.toString(16).padStart(2, '0')).join(''), base = v.getUint32(44 + 20 * s, true); if (base + 8 > b.length) continue;
    const count = min(4096, v.getUint32(base + 4, true)), name = FMTID[fid] || 'property set', names = PID[name] || {}, set = { name, fmtid: fid, properties: {} }, items = [];
    for (let i = 0; i < count && base + 16 + 8 * i <= b.length; i++) items.push([v.getUint32(base + 8 + 8 * i, true), base + v.getUint32(base + 12 + 8 * i, true)]);
    let cp = 1252; const cpi = items.find(x => x[0] === 1); if (cpi && cpi[1] + 6 <= b.length && v.getUint16(cpi[1], true) === 2) cp = v.getUint16(cpi[1] + 4, true);
    const str = (o, len) => { if (cp === 1200) return utf16(b, o, len * 2).replace(/\0+$/, ''); const raw = b.subarray(o, min(b.length, o + len)); let t; try { t = new TextDecoder(cp === 65001 ? 'utf-8' : 'windows-' + cp).decode(raw); } catch (e) { t = latin(raw, 0, raw.length); } return t.replace(/\0+$/, ''); };
    const dict = {}; const di = items.find(x => x[0] === 0);
    if (di) { let o = di[1]; const m = o + 4 <= b.length ? min(4096, v.getUint32(o, true)) : 0; o += 4; for (let i = 0; i < m && o + 8 <= b.length; i++) { const id = v.getUint32(o, true), len = v.getUint32(o + 4, true); o += 8; if (len > 4096) break; dict[id] = str(o, len); o += cp === 1200 ? ((len * 2 + 3) & ~3) : len; } }
    for (const [id, o] of items) {
      if (id === 0 || id === 1 || o + 8 > b.length) continue; const type = v.getUint16(o, true), key = dict[id] || names[id] || `pid${id}`; let val;
      if (type === 2) val = v.getInt16(o + 4, true); else if (type === 3 || type === 22) val = v.getInt32(o + 4, true); else if (type === 19 || type === 23) val = v.getUint32(o + 4, true); else if (type === 4) val = v.getFloat32(o + 4, true); else if (type === 5 && o + 12 <= b.length) val = v.getFloat64(o + 4, true);
      else if (type === 11) val = v.getInt16(o + 4, true) !== 0; else if (type === 30) val = str(o + 8, min(65536, v.getUint32(o + 4, true))); else if (type === 31) val = utf16(b, o + 8, min(131072, v.getUint32(o + 4, true) * 2)).replace(/\0+$/, '');
      else if (type === 64 && o + 12 <= b.length) val = fileTime(v.getUint32(o + 4, true), v.getUint32(o + 8, true));
      else if (type === 71) { const size = v.getUint32(o + 4, true), tag = v.getInt32(o + 8, true), fmt = tag === -1 || tag === -2 ? v.getUint32(o + 12, true) : null, data = b.subarray(o + (fmt === null ? 12 : 16), min(b.length, o + 8 + size)); let img = sniffImage(data, 32); if (!img && fmt === 8) { const bmp = dibToBmp(data); if (bmp) img = { mime: 'image/bmp', bytes: bmp }; } if (!img && fmt === 3 && data.length > 8) img = { mime: 'image/wmf', bytes: data.subarray(8) }; if (!img && fmt === 14) img = { mime: 'image/emf', bytes: data }; if (img) set.thumbnail = img; continue; }
      if (val !== undefined && val !== null && val !== '') set.properties[key] = val;
    }
    out.push(set);
  }
  return out;
}
const INV_SETS = { 'Inventor Summary Information': 'SummaryInformation', 'Inventor Document Summary Information': 'DocumentSummaryInformation', 'Design Tracking Properties': 'Inventor Design Tracking Properties' };
const DISPLAYABLE = ['image/png', 'image/jpeg', 'image/bmp', 'image/gif'];
const preview = (img, source) => { const size = imageSize(img); return { mime: img.mime, bytes: img.bytes, source, displayable: DISPLAYABLE.includes(img.mime), ...(size ? { width: size[0], height: size[1] } : {}) }; };
const describePreview = p => `preview image (${p.mime}${p.width ? `, ${p.width} × ${p.height} px` : ''}, ${p.bytes.length} bytes) from ${p.source}${p.displayable ? '' : ' — a vector metafile that browsers cannot display'}`;
const parasolidInfo = b => { const head = latin(b, 0, 2048), m = /TRANSMIT FILE\s*(\(partition\))?\s*created by modeller version (\d+)/.exec(head), s = /SCH_(\d+)_(\d+)/.exec(head); return { encoding: head.startsWith('PS') ? 'binary' : 'text', partition: !!(m && m[1]), modellerVersion: m ? m[2] : null, schema: s ? s[0] : null }; };

// Property sets, thumbnail and stream inventory common to every OLE2-based CAD document.
function cfbDocument(b, what) {
  const cfb = cfbOpen(b), meta = { container: `OLE2 compound file v${cfb.version} (${cfb.sectorSize}-byte sectors)`, streams: cfb.entries.filter(e => e.type === 'stream').length, storages: cfb.entries.filter(e => e.type === 'storage').length, properties: {} }, recovered = [], extra = {}; let pv = null;
  for (const e of cfb.entries) {
    if (e.type !== 'stream' || e.name.charCodeAt(0) !== 5 || e.size > 4e6) continue;
    let sets = []; try { sets = parsePropertySets(cfb.read(e)); } catch (err) { sets = []; }
    for (const s of sets) {
      // Inventor names its sets in a "Property Set Name" property and reuses the standard property ids
      const own = s.properties['Property Set Name']; if (own) { delete s.properties['Property Set Name']; s.name = String(own); const ids = PID[INV_SETS[s.name]] || {}; for (const k of Object.keys(s.properties)) { const m = /^pid(\d+)$/.exec(k); if (m && ids[m[1]]) { s.properties[ids[m[1]]] = s.properties[k]; delete s.properties[k]; } } }
      delete s.properties.pid2147483648; if (/^_/.test(s.name)) continue;
      if (Object.keys(s.properties).length) { const key = s.name === 'property set' ? e.name.slice(1) : s.name; meta.properties[key] = { ...(meta.properties[key] || {}), ...s.properties }; } if (s.thumbnail && (!pv || (!pv.displayable && DISPLAYABLE.includes(s.thumbnail.mime)))) pv = preview(s.thumbnail, `the ${s.name} thumbnail property`); }
  }
  for (const e of cfb.find(/(^|\/)(PreviewPNG|Preview|Thumbnail|RevitPreview[\d.]*|JPreview)$/i)) { if (e.size < 64 || e.size > 2e7) continue; const img = sniffImage(cfb.read(e), 4096); if (img && (!pv || !pv.displayable || img.bytes.length > pv.bytes.length)) pv = preview(img, `stream "${e.path}"`); }
  const nprops = Object.values(meta.properties).reduce((a, s) => a + Object.keys(s).length, 0);
  if (nprops) recovered.push(`${nprops} document properties from ${Object.keys(meta.properties).length} OLE property set(s)`);
  if (pv) { extra.preview = pv; recovered.push(describePreview(pv)); }
  meta.streamList = cfb.entries.filter(e => e.type === 'stream').sort((a, c) => c.size - a.size).slice(0, 60).map(e => `${e.path.replace(/[\x00-\x1f]/g, '·')} (${e.size} bytes)`);
  recovered.push(`${what} stream inventory: ${meta.streams} streams in ${meta.storages} storages`);
  return { cfb, meta, recovered, extra };
}
const flatProps = meta => { const o = {}; for (const s of Object.values(meta.properties || {})) for (const [k, v] of Object.entries(s)) if (typeof v !== 'object' && !/^pid\d+$/.test(k) && String(v).length < 200 && !(k in o)) o[k] = v; return o; };

// ---------------------------------------------------------------- SOLIDWORKS
// 2015+ files are a sequence of chunks: marker 14 00 06 00 08 00, u32 tag, u32 CRC-32 of the inflated data, u32 compressed
// size, u32 inflated size, u32 name length, the stream name with the two nibbles of every byte swapped, raw DEFLATE data.
const swap = c => ((c << 4) & 0xF0) | (c >> 4);
export function swStreams(bytes) {
  const b = F.toBytes(bytes), v = dv(b), out = new Map();
  for (let i = 0; i + 26 < b.length; i++) {
    if (b[i] !== 0x14 || b[i + 1] !== 0 || b[i + 2] !== 6 || b[i + 3] !== 0 || b[i + 4] !== 8 || b[i + 5] !== 0) continue;
    const crc = v.getUint32(i + 10, true), cs = v.getUint32(i + 14, true), us = v.getUint32(i + 18, true), nl = v.getUint32(i + 22, true); if (!nl || nl > 255 || i + 26 + nl + cs > b.length || us > 1.5e9) continue;
    let name = ''; for (let k = 0; k < nl; k++) name += String.fromCharCode(swap(b[i + 26 + k])); if (!/^[\x20-\x7e]+$/.test(name)) continue;
    const prev = out.get(name); if (prev && prev.size >= us) continue; // the same stream can be present in several saved revisions; the largest is kept
    out.set(name, { name, offset: i + 26 + nl, compressed: cs, size: us, crc });
  }
  const read = name => { const e = out.get(name); if (!e) return null; if (!e.data) { let d; try { d = e.size === 0 ? new Uint8Array(0) : F.inflateRaw(b, e.offset, e.size).data; } catch (err) { return null; } if (d.length !== e.size || (e.size && F.crc32(d) !== e.crc)) return null; e.data = d; } return e.data; };
  return { entries: [...out.values()], read };
}
// The Parasolid body of a configuration: 4-byte length, 16-byte class id, two sizes, then a zlib stream holding a binary
// Parasolid transmit file. Returned only when the inflated data carries the Parasolid transmit signature.
function swPartition(d) { if (!d || d.length < 40) return null; for (const o of [28, 24, 32, 20, 36]) { if ((d[o] & 0x0F) !== 8 || ((d[o] << 8) | d[o + 1]) % 31) continue; try { const r = F.inflateRaw(d, o + 2).data; if (/^PS\0\0\0.: TRANSMIT FILE|^\*\*ABCDEFGHIJKLM/.test(latin(r, 0, 40))) return r; } catch (e) { /* not at this offset */ } } return null; }
const wideStrings = (d, minLen = 3, limit = 400) => { const out = []; let cur = ''; for (let i = 0; i + 1 < d.length && out.length < limit; i += 2) { const c = d[i] | (d[i + 1] << 8); if (c >= 32 && c < 0xD800 && c !== 0x7F && (c < 0x80 || c > 0xA0)) cur += String.fromCharCode(c); else { if (cur.length >= minLen) out.push(cur); cur = ''; if (d[i + 1] !== 0 && d[i] === 0) i--; } } if (cur.length >= minLen) out.push(cur); return out; };
// The component instance tree of an assembly → {name, bbox (metres) | null, files: [path], components: [{name, instance,
// path: [names from the top], depth, file, type, configuration, configurationId, suppressed, hidden, virtual,
// matrix: 3 × 4 row-major acting on column vectors (metres) | null}]}. Every model (the assembly and each sub-assembly)
// lists its own swReference children; a reference carries a 4 × 4 "swTransform" whose first three rows are the rotation
// acting on row vectors and whose fourth row is the translation (p' = p · R + t) relative to the model that owns it.
// Sub-assemblies are expanded here, so each part instance comes out with its placement in the top assembly.
export function swInstanceTree(xml) {
  if (!/<swReference\b/.test(xml)) return null;
  const at = (tag, k) => { const m = new RegExp('\\b' + k + '="([^"]*)"').exec(tag); return m ? xmlText(m[1]) : ''; }, files = new Map(), models = new Map(); let cur = null;
  for (const m of xml.matchAll(/<swFile\b[^>]*>|<swModel\b[^>]*>|<\/swModel>|<swReference\b[^>]*>/g)) {
    const tag = m[0];
    if (tag.startsWith('<swFile')) files.set(at(tag, 'id'), { path: at(tag, 'swPath'), type: at(tag, 'swDocType') });
    else if (tag === '</swModel>') cur = null;
    else if (tag.startsWith('<swModel')) { const mo = { name: at(tag, 'swName'), file: at(tag, 'swFileRef'), bbox: at(tag, 'swBoundingBox').trim().split(/\s+/).map(Number), refs: [] }; models.set(at(tag, 'id'), mo); cur = /\/>$/.test(tag) ? null : mo; }
    else if (cur) {
      const a = at(tag, 'swTransform').trim().split(/\s+/).map(Number); let matrix = null;
      if (a.length === 16 && a.every(Number.isFinite) && abs(a[3]) + abs(a[7]) + abs(a[11]) < 1e-9) { const sc = a[15] > 0 ? a[15] : 1; matrix = [a[0] * sc, a[4] * sc, a[8] * sc, a[12], a[1] * sc, a[5] * sc, a[9] * sc, a[13], a[2] * sc, a[6] * sc, a[10] * sc, a[14]]; }
      cur.refs.push({ name: at(tag, 'swName'), instance: +at(tag, 'swReferenceNumber') || 1, model: at(tag, 'swModelRef'), configuration: at(tag, 'swConfigurationName'), configurationId: at(tag, 'swConfigurationId'), suppressed: at(tag, 'swSuppressed') === 'YES', hidden: at(tag, 'swHidden') === 'YES', virtual: at(tag, 'swIsVirtualComponent') === 'YES', matrix });
    }
  }
  let top = null; for (const m of xml.matchAll(/<swConfiguration\b[^>]*>/g)) if (!top || /swMostRecentConfiguration="YES"/.test(m[0])) top = models.get(at(m[0], 'swModelRef')) || top;
  if (!top) { const used = new Set(); for (const mo of models.values()) for (const r of mo.refs) used.add(r.model); for (const [id, mo] of models) if (mo.refs.length && !used.has(id)) top = mo; }
  if (!top) return null;
  const components = [], walk = (mo, M, path, off, depth) => {
    if (depth > 24) return;
    for (const r of mo.refs) {
      const child = models.get(r.model), file = child ? files.get(child.file) : null, here = r.matrix ? (M ? matMul(M, r.matrix) : r.matrix) : null, gone = off || r.suppressed, p = [...path, `${r.name}<${r.instance}>`];
      if (child && child.refs.length && child !== mo) { walk(child, M === null && path.length ? null : here, p, gone, depth + 1); continue; }
      components.push({ name: r.name, instance: r.instance, path: p, depth, file: file ? file.path : '', type: file ? file.type : '', configuration: r.configuration, configurationId: r.configurationId, suppressed: gone, hidden: r.hidden, virtual: r.virtual, matrix: path.length && M === null ? null : here });
    }
  };
  walk(top, null, [], false, 0);
  return { name: top.name, bbox: top.bbox.length === 6 && top.bbox.every(Number.isFinite) ? top.bbox : null, files: [...new Set(components.map(c => c.file).filter(Boolean))], components };
}
export function parseSolidWorks(bytes) {
  const b = F.toBytes(bytes), recovered = [], notRecovered = [], extra = {}, meta = {}; let get, names;
  if (startsWith(b, [0xD0, 0xCF, 0x11, 0xE0])) { const d = cfbDocument(b, 'SOLIDWORKS (OLE2 layout, 2014 and earlier)'); Object.assign(meta, d.meta); recovered.push(...d.recovered); Object.assign(extra, d.extra); names = d.cfb.entries.filter(e => e.type === 'stream').map(e => e.path); get = n => d.cfb.read(n); meta.layout = 'OLE2 compound file (SOLIDWORKS 2014 and earlier)'; }
  else { const s = swStreams(b); if (!s.entries.length) throw new Error('SOLIDWORKS: neither an OLE2 compound file nor the 2015+ chunk container (no stream markers found)'); names = s.entries.map(e => e.name); get = s.read; meta.layout = 'SOLIDWORKS 2015+ chunk container (DEFLATE-compressed named streams)'; meta.streams = names.length; meta.streamList = s.entries.slice().sort((a, c) => c.size - a.size).slice(0, 60).map(e => `${e.name} (${e.size} bytes)`); recovered.push(`stream inventory: ${names.length} named streams, each inflated and CRC-checked on demand`); }
  const text = n => { const d = get(n); return d ? F.toText(d[0] === 0x3C || d[0] === 0xFF || d[0] === 0xEF ? d : d.subarray(1)) : ''; };
  // previews
  if (!extra.preview || !extra.preview.displayable) for (const n of ['PreviewPNG', 'Preview', ...names.filter(x => /^ModelBasedImages\/ModelBasedImage/.test(x))]) { const d = names.includes(n) ? get(n) : null, img = d ? sniffImage(d, 64) : null; if (img) { extra.preview = preview(img, `stream "${n}"`); recovered.push(describePreview(extra.preview)); break; } }
  // document properties (OPC-style XML parts in the 2015+ layout)
  const props = {}; for (const [n, tags] of [['docProps/core.xml', ['dc:title', 'dc:subject', 'dc:creator', 'dc:description', 'cp:keywords', 'dc:lastModifiedBy', 'dcterms:created', 'dcterms:modified']], ['docProps/app.xml', ['Application', 'AppVersion', 'Company']]]) { const t = names.includes(n) ? text(n) : ''; for (const tag of tags) { const m = new RegExp('<' + tag + '\\b[^>]*>([^<]*)</' + tag + '>').exec(t); if (m && m[1].trim()) props[tag.replace(/^\w+:/, '')] = xmlText(m[1].trim()); } }
  const custom = {}; for (const n of names.filter(x => /^docProps\/(custom|Config-\d+-Properties)\.xml$/.test(x))) { const t = text(n), scope = /Config-(\d+)/.exec(n); for (const m of t.matchAll(/<property\b[^>]*\bname="([^"]+)"[^>]*>\s*<vt:\w+>([^<]*)<\/vt:\w+>/g)) { const k = xmlText(m[1]).trim(), val = xmlText(m[2]).trim(); if (k && val) custom[scope ? `${k} [configuration ${scope[1]}]` : k] = val; } }
  if (Object.keys(props).length) { meta.documentProperties = props; recovered.push(`document properties — ${Object.entries(props).slice(0, 8).map(([k, x]) => `${k}: ${String(x).slice(0, 80)}`).join('; ')}`); }
  if (Object.keys(custom).length) { meta.customProperties = custom; recovered.push(`${Object.keys(custom).length} custom properties (${Object.entries(custom).slice(0, 6).map(([k, x]) => `${k} = ${String(x).replace(/\s+/g, ' ').slice(0, 40)}`).join('; ')})`); }
  if (props.AppVersion) { const major = parseInt(props.AppVersion, 10); if (major >= 8 && major < 60) meta.release = `SOLIDWORKS ${1992 + major} (application version ${props.AppVersion})`; }
  const mo = names.map(x => /^_MO_VERSION_(\d+)/.exec(x)).find(Boolean); if (mo) meta.fileSchema = `_MO_VERSION_${mo[1]}`;
  // configurations
  const cfgIdx = [...new Set(names.map(x => /(?:^|\/)Config-(\d+)(?:$|-)/.exec(x)).filter(Boolean).map(m => +m[1]))].sort((a, c) => a - c);
  const hdr = names.find(x => /(^|\/)CMgrHdr2?$/.test(x)), hd = hdr ? get(hdr) : null; let cfgNames = [];
  if (hd) { const seen = new Set(); for (let i = 0; i + 6 < hd.length; i++) if (hd[i] === 0xFF && hd[i + 1] === 0xFE && hd[i + 2] === 0xFF) { const n = hd[i + 3]; if (n && i + 4 + 2 * n <= hd.length) { const s = utf16(hd, i + 4, 2 * n); if (/^[^\x00-\x1f]+$/.test(s) && !seen.has(s)) { seen.add(s); cfgNames.push(s); } i += 3 + 2 * n; } } }
  if (cfgNames.length || cfgIdx.length) { meta.configurations = cfgNames.length ? cfgNames : cfgIdx.map(i => `configuration ${i} (name not readable)`); recovered.push(`${meta.configurations.length} configuration(s): ${meta.configurations.slice(0, 12).join(', ')}`); }
  const kw = names.includes('swXmlContents/KeyWords') ? text('swXmlContents/KeyWords') : ''; if (kw) { const feats = [...kw.matchAll(/<(?:swFeature|Feature)\b[^>]*\b(?:Name|name)="([^"]+)"/g)].map(m => xmlText(m[1])); if (feats.length) { meta.featureNames = feats.slice(0, 200); recovered.push(`${feats.length} feature and item names from the keyword index (names only — no parameters or sketches)`); } }
  // assembly structure: the component instance tree SOLIDWORKS writes as XML (files, configurations, placements, extents)
  if (names.includes('swXmlContents/COMPINSTANCETREE')) { const asm = swInstanceTree(text('swXmlContents/COMPINSTANCETREE')); if (asm && asm.components.length) { extra.assembly = asm; meta.components = asm.components.slice(0, 200).map(c => `${c.path.join(' / ')} → ${c.file || 'virtual component'}${c.configuration ? ` [${c.configuration}]` : ''}${c.suppressed ? ' (suppressed)' : ''}`); recovered.push(`assembly structure: ${asm.components.length} part instance(s) (sub-assemblies expanded) of ${asm.files.length} referenced file(s) with their placement transforms (${asm.files.slice(0, 8).map(f => f.split(/[\\/]/).pop()).join(', ')}${asm.files.length > 8 ? ', …' : ''})`); } }
  // embedded Parasolid bodies
  const embedded = []; for (const n of names.filter(x => /(^|\/)Config-\d+-Partition$/.test(x))) { const ps = swPartition(get(n)); if (ps) { const info = parasolidInfo(ps), c = /Config-(\d+)/.exec(n)[1]; embedded.push({ name: `configuration-${c}.x_b`, mime: 'application/octet-stream', format: 'parasolid', bytes: ps, description: `Parasolid ${info.encoding} transmit${info.partition ? ' (partition)' : ''} of configuration ${c}${cfgNames[+c] ? ` "${cfgNames[+c]}"` : ''}, modeller version ${info.modellerVersion || 'unknown'}, schema ${info.schema || 'unknown'}`, ...info }); } }
  if (embedded.length) { extra.embedded = embedded; meta.parasolid = embedded.map(e => e.description); recovered.push(`${embedded.length} embedded Parasolid body stream(s) unpacked intact (${embedded.map(e => `${e.name}: ${e.bytes.length} bytes, ${e.schema || 'schema not read'}`).join('; ')}) — the exact B-rep of the part in Parasolid partition form, kept for saving`); }
  else if (names.some(x => /Partition$/.test(x))) notRecovered.push('the body partition streams are present but did not unpack to a Parasolid transmit file');
  if (names.some(x => /DisplayLists/.test(x))) notRecovered.push('the stored display lists (the graphics tessellation SOLIDWORKS draws on screen) use an undocumented encoding and were not decoded');
  notRecovered.push('the feature tree, sketches, dimensions and mates are stored as serialised SOLIDWORKS objects and were not interpreted');
  if (!extra.preview) notRecovered.push('no preview image stream was found');
  return { meta, recovered, notRecovered, extra, geometry: null };
}

// ---------------------------------------------------------------- Inventor, Solid Edge, Revit and other OLE2 documents
export function parseOleCad(bytes, product) {
  const d = cfbDocument(F.toBytes(bytes), product), notRecovered = [], flat = flatProps(d.meta);
  if (Object.keys(flat).length) d.recovered.unshift(`properties — ${Object.entries(flat).slice(0, 10).map(([k, v]) => `${k}: ${String(v).slice(0, 60)}`).join('; ')}`);
  if (product === 'Inventor') {
    const seg = d.cfb.entries.filter(e => /RSeStorage\//.test(e.path) && e.type === 'stream'); if (seg.length) d.meta.segments = seg.length;
    const info = d.cfb.find(/RSeDbRevisionInfo$|RSeDb$/)[0], s = info ? wideStrings(d.cfb.read(info), 6, 60).filter(x => /Schema|Saved|Build|\d{4} \(/.test(x)) : []; if (s.length) { d.meta.revisionInfo = s.slice(0, 8); d.recovered.push(`save history — ${s.slice(0, 4).join('; ')}`); }
    // files this document references (the parts and sub-assemblies of an assembly, a derived part's base): full paths in UFRxDoc
    const ufr = d.cfb.entries.find(e => e.type === 'stream' && e.name === 'UFRxDoc'), ud = ufr ? d.cfb.read(ufr) : null;
    if (ud) { const refs = []; for (const off of [0, 1]) for (const m of utf16(ud, off, ud.length - off).matchAll(/(?:[A-Za-z]:\\|\\\\)[^\x00-\x1f"*<>?|:]{1,240}?\.(?:ipt|iam)/gi)) { const pth = m[0]; if (!/\\Templates\\/i.test(pth) && !refs.includes(pth)) refs.push(pth); }
      const own = refs.length ? refs[0] : ''; const others = refs.slice(1); if (own) d.meta.savedAs = own; if (others.length) { d.extra.references = others; d.meta.references = others.slice(0, 200); d.recovered.push(`${others.length} referenced file(s): ${others.slice(0, 8).map(x => x.split(/[\\/]/).pop()).join(', ')}${others.length > 8 ? ', …' : ''}`); } }
    // model segments: "M…" meta streams name each segment, the "B…" stream of the same key holds its data as one zlib stream
    const segs = inventorSegments(d.cfb), found = [], zstd = segs.filter(g => g.zstd).length;
    if (segs.length) { d.meta.segmentNames = segs.map(g => `${g.name} (${g.data.length} bytes${g.inflated ? ' inflated' : g.zstd ? ', Zstandard-compressed, not inflated' : ''})`); d.recovered.push(zstd ? `${segs.length} model segment(s) named: ${segs.map(g => g.name).join(', ')}` : `${segs.length} model segment(s) unpacked: ${segs.map(g => g.name).join(', ')}`); }
    // assembly: occurrences, their files and their placements
    if (segs.some(g => g.name === 'AmRxSegment')) { let asm = null; try { asm = inventorOccurrences(d.cfb, segs); } catch (e) { asm = null; }
      if (asm && asm.components.length) { d.extra.assembly = asm; d.meta.components = asm.components.slice(0, 200).map(c => `${c.file.split(/[\\/]/).pop()}:${c.instance}${c.subAssembly ? ' (sub-assembly)' : ''}`); d.recovered.push(`assembly structure: ${asm.components.length} occurrence(s) of ${new Set(asm.components.map(c => c.file)).size} file(s), each with its placement transform (rotation and translation in centimetres) from the assembly's result segment`); }
      else if (!zstd) notRecovered.push('the occurrence list of the assembly result segment is not in the layout this reader knows: component placements were not decoded'); }
    // an assembly segment's save without any body record (no assembly-level bodies) is not offered
    const bodiless = f => { if (!/^Am/.test(f.segment)) return false; try { return !readSAT(f.bytes).counts.body; } catch (e) { return false; } };
    for (const g of segs) for (const h of findACIS(g.data)) { const f = { ...h, segment: g.name, brep: /b-?rep/i.test(g.name), bytes: g.data.subarray(h.offset) }; if (!bodiless(f)) found.push(f); }
    const use = found.some(f => f.brep) ? found.filter(f => f.brep) : found, aside = found.length - use.length;
    if (use.length) {
      d.extra.embedded = use.map((f, i) => ({ name: `brep-${i + 1}.sab`, mime: 'application/octet-stream', format: 'acis', bytes: f.bytes, segment: f.segment, description: `Autodesk ShapeManager binary B-rep ("${f.signature}") ${i + 1} of ${use.length} from segment ${f.segment}` }));
      d.meta.shapeManager = d.extra.embedded.map(e => e.description);
      d.recovered.push(`${use.length} ShapeManager (ACIS-derived) B-rep stream(s) located in ${[...new Set(use.map(f => f.segment))].join(', ')}${aside ? `; ${aside} further stream(s) in other segments (sketch and construction geometry) left aside` : ''}`);
    } else if (d.extra.assembly) { /* an assembly: its geometry is in the part files it references */ } else notRecovered.push(zstd ? `the ${zstd} model segment(s) are stored with Zstandard compression (written by Inventor 2026; the saving release is in the properties above) and the Zstandard decoder was not available for this import: neither the ShapeManager B-rep nor the assembly occurrences were unpacked` : segs.length ? 'no ShapeManager B-rep stream was found in the model segments' : 'the model segments could not be unpacked (no meta / data stream pairs in the layout this reader knows): no B-rep located');
    notRecovered.push('the feature tree, sketches, parameters, browser, graphics and design-view segments are stored in Inventor\'s proprietary object serialisation and were not interpreted');
  } else {
    if (product === 'Solid Edge') { const found = []; for (const e of d.cfb.entries) if (e.type === 'stream' && e.size > 512) { const data = d.cfb.read(e); if (data) for (const f of findParasolid(data)) if (f.guise !== 'deltas') found.push({ ...f, stream: e.path.replace(/[\x00-\x1f]/g, '·') }); }
      if (found.length) { d.extra.embedded = found.map((f, i) => ({ name: `parasolid-${i + 1}.x_b`, mime: 'application/octet-stream', format: 'parasolid', bytes: f.bytes, description: `Parasolid binary transmit (${f.guise}) from stream "${f.stream}", modeller version ${f.modellerVersion}, schema ${f.schema || 'unknown'}`, encoding: f.encoding, partition: f.partition, modellerVersion: f.modellerVersion, schema: f.schema })); d.recovered.push(`${found.length} embedded Parasolid stream(s) located by signature`); } }
    notRecovered.push(d.extra.embedded ? `the ${product} feature and document streams are in the vendor's proprietary encoding and were not decoded` : `the ${product} model streams are in the vendor's proprietary encoding and were not decoded`);
  }
  if (!d.extra.preview) notRecovered.push('no thumbnail was found in the property sets or preview streams');
  return { meta: d.meta, recovered: d.recovered, notRecovered, extra: d.extra, geometry: null };
}

// Inventor model segments → [{name, key, data, inflated, zstd}]. A meta stream "M<key>" carries the segment name; the data
// stream "B<key>" is a short header followed by one zlib stream (older releases store the data raw; Inventor 2026 writes a
// Zstandard frame instead, recognised by its magic number and inflated when the optional decoder F.codecs.zstd is loaded;
// `zstd` is true when such a frame was left packed).
export function inventorSegments(cfb) {
  const out = [];
  for (const e of cfb.entries) {
    if (e.type !== 'stream' || !/(^|\/)RSeStorage\/B[^/]*$/.test(e.path) || e.size < 8) continue;
    const key = e.name.slice(1), me = cfb.entries.find(x => x.type === 'stream' && x.path === e.path.slice(0, e.path.length - e.name.length) + 'M' + key); let name = '';
    if (me) { const m = cfb.read(me); if (m) { const t = latin(m, 0, min(m.length, 512)).replace(/\0/g, ''), q = /[A-Za-z][A-Za-z0-9]{1,40}Segment/.exec(t); if (q) name = q[0]; } }
    const d = cfb.read(e); if (!d) continue; let data = null, zstd = false;
    for (let i = 0; i < min(64, d.length - 2) && !data; i++) if (d[i] === 0x78 && (d[i + 1] === 0x01 || d[i + 1] === 0x5E || d[i + 1] === 0x9C || d[i + 1] === 0xDA)) { try { const r = F.inflateRaw(d, i + 2); if (r.data.length > 0) data = r.data; } catch (err) { data = null; } }
    if (!data) for (let i = 0; i < min(64, d.length - 4) && !zstd; i++) if (d[i] === 0x28 && d[i + 1] === 0xB5 && d[i + 2] === 0x2F && d[i + 3] === 0xFD) { zstd = true; if (F.codecs.zstd) { try { const r = F.codecs.zstd(d.subarray(i)); if (r && r.length) data = r; } catch (err) { data = null; } } } // Zstandard frame magic
    out.push({ name: name || `segment ${key.slice(0, 8)}`, key, data: data || d, inflated: !!data, zstd: zstd && !data, codec: data ? (zstd ? 'Zstandard' : 'zlib') : '' });
  }
  return out;
}

// ---- Inventor assemblies: occurrences and their placements
// The object records of an inflated segment, in file order → [{type, at, len}]. Each record is a zero byte, a u32 class
// index, the payload and a u32 payload length; the stream ends with a zero byte and FF FF FF FF, so it is walked backwards.
export function inventorRecords(data) {
  const d = F.toBytes(data), v = dv(d), out = []; let end = -1; for (let o = d.length - 5, g = 0; o > 0 && g < 4096; o--, g++) if (d[o] === 0 && v.getUint32(o + 1, true) === 0xFFFFFFFF) { end = o; break; } if (end < 0) return out;
  for (let p = end; p > 9;) { const len = v.getUint32(p - 4, true), s = p - 9 - len; if (s < 0 || d[s] !== 0) break; out.push({ type: v.getUint32(s + 1, true), at: s + 5, len }); p = s; }
  return out.reverse();
}
// A packed 4 × 4 transform: a u32 mask, then one float64 for every cell (row-major) that the mask does not fix. Cell k
// is fixed to 1 by bit k alone, to -1 by bits k and k + 16, to 0 by bit k + 16 alone. The rows are [R | t]: the matrix
// acts on column vectors, lengths in centimetres. → {matrix: 12 numbers (3 × 4 row-major), end} or null when the cells
// do not form a rigid motion.
export function inventorTransform(data, offset, limit) {
  const d = F.toBytes(data), end = limit === undefined ? d.length : min(limit, d.length); if (offset + 4 > end) return null; const v = dv(d), m = v.getUint32(offset, true), a = []; let p = offset + 4;
  for (let k = 0; k < 16; k++) { const lo = (m >>> k) & 1, hi = (m >>> (k + 16)) & 1; if (lo) a.push(hi ? -1 : 1); else if (hi) a.push(0); else { if (p + 8 > end) return null; const x = v.getFloat64(p, true); if (!Number.isFinite(x) || abs(x) > 1e7) return null; a.push(x); p += 8; } }
  if (abs(a[12]) + abs(a[13]) + abs(a[14]) > 1e-9 || abs(a[15] - 1) > 1e-9) return null;
  for (let i = 0; i < 3; i++) { const j = (i + 1) % 3; if (abs(hypot(a[4 * i], a[4 * i + 1], a[4 * i + 2]) - 1) > 1e-6 || abs(a[4 * i] * a[4 * j] + a[4 * i + 1] * a[4 * j + 1] + a[4 * i + 2] * a[4 * j + 2]) > 1e-6) return null; }
  return { matrix: a.slice(0, 12), end: p };
}
// The reference table of UFRxDoc → {files: [{path, key, count}], occurrences: [file | null]}. Every referenced file is
// a length-prefixed UTF-16 path followed by a 70-byte record (16-byte document id at +22, number of occurrences at +58);
// the first path is the document itself. After the table comes one property block per occurrence, in occurrence order,
// whose property 3 repeats the document id of its file (copies of a file share that id: property 5, the file's base
// name, then tells them apart).
export function inventorReferences(bytes) {
  const ud = F.toBytes(bytes), v = dv(ud), all = [];
  for (let o = 0; o + 12 < ud.length; o++) { const n = v.getUint32(o, true); if (n < 8 || n > 400 || o + 4 + 2 * n + 70 > ud.length || ud[o + 5] !== 0 || ud[o + 4 + 2 * n - 1] !== 0) continue; const s = utf16(ud, o + 4, 2 * n); if (!/^(?:[A-Za-z]:\\|\\\\)[^\x00-\x1f"*<>?|]+\.(?:ipt|iam)$/i.test(s)) continue; const e = o + 4 + 2 * n; all.push({ path: s, key: Array.from(ud.subarray(e + 22, e + 38)).join(','), count: v.getUint32(e + 58, true), end: e + 70 }); o += 3 + 2 * n; }
  const files = all.slice(1), occurrences = [], hits = []; if (!files.length) return { files, occurrences };
  const prop = (o, id) => ud[o] <= 1 && ud[o + 1] === id && ud[o + 2] === 1 && !ud[o + 3] && !ud[o + 4] && !ud[o + 5] && ud[o + 6] === id, stem = pth => pth.split(/[\\/]/).pop().replace(/\.[^.]+$/, '').toLowerCase();
  for (let o = all[all.length - 1].end; o + 23 <= ud.length; o++) if (prop(o, 3)) { const key = Array.from(ud.subarray(o + 7, o + 23)).join(','); if (files.some(f => f.key === key)) { hits.push({ o, key }); o += 22; } }
  hits.forEach((h, k) => { let cand = files.filter(f => f.key === h.key); if (cand.length > 1) { const stop = k + 1 < hits.length ? hits[k + 1].o : ud.length; let name = null; for (let o = h.o + 23; o + 11 < stop && name === null; o++) if (prop(o, 5)) { const n = v.getUint32(o + 7, true); if (n > 0 && n < 260 && o + 11 + 2 * n <= ud.length) name = utf16(ud, o + 11, 2 * n).toLowerCase(); } cand = cand.filter(f => stem(f.path) === name); } occurrences.push(cand.length === 1 ? cand[0] : null); });
  return { files: files.map(f => ({ path: f.path, key: f.key, count: f.count })), occurrences: occurrences.map(f => (f ? files.indexOf(f) : -1)) };
}
// Occurrences of an assembly → {components: [{id, name, file, instance, subAssembly, matrix: 3 × 4 row-major acting on
// column vectors, centimetres, relative to this assembly}], files: [path]} or null when the layout is not the one this
// reader knows (nothing is guessed). The result segment (AmRxSegment) holds, as three record classes that are recognised
// by their layout: one placement record per occurrence (object reference, packed transform, occurrence id 12 bytes
// before the end), one 24-byte handle per occurrence (its id) and one record per component definition listing the
// object numbers of its handles. The component definitions follow the order of the occurrence blocks of UFRxDoc, which
// name the files; the run lengths of both must agree, and every handle, id, placement and file must pair off exactly.
export function inventorOccurrences(cfb, segments) {
  const rx = (segments || inventorSegments(cfb)).find(g => g.name === 'AmRxSegment'), ue = cfb.entries.find(e => e.type === 'stream' && e.name === 'UFRxDoc'); if (!rx || !rx.inflated || !ue) return null;
  const d = rx.data, v = dv(d), recs = inventorRecords(d), ub = cfb.read(ue), refs = ub ? inventorReferences(ub) : null; if (!recs.length || !refs || !refs.files.length) return null;
  const byType = new Map(); recs.forEach((r, i) => { r.i = i; if (!byType.has(r.type)) byType.set(r.type, []); byType.get(r.type).push(r); });
  const placement = r => { if (r.len < 60) return null; for (const off of [10, 14]) { const q = inventorTransform(d, r.at + off, r.at + r.len - 12); if (q) return { occ: v.getUint32(r.at + r.len - 12, true), matrix: q.matrix }; } return null; };
  const handle = r => (r.len === 24 && v.getUint32(r.at, true) === 0x30000002 && v.getUint32(r.at + 4, true) === 1 ? { occ: v.getUint32(r.at + 16, true) } : null);
  const definition = r => { if (r.len < 40 || v.getUint32(r.at, true) !== 0x30000002) return null; const n = v.getUint32(r.at + 4, true); if (n < 1 || n > 100000 || v.getUint32(r.at + 8, true) < n || 20 + 4 * n > r.len) return null; let p = r.at + 16; const ids = []; for (let k = 0; k < n; k++, p += 4) ids.push(v.getUint32(p, true)); const nl = v.getUint32(p, true); p += 4; if (nl > 260 || p + 2 * nl + 8 > r.at + r.len || d[p + 2 * nl] !== 1) return null; return { ids, name: utf16(d, p, 2 * nl), kind: d[p + 2 * nl + 1] }; };
  const all = (f, skip) => { const res = []; for (const [t, list] of byType) { if (skip.includes(t)) continue; const got = list.map(f); if (got.every(Boolean)) res.push({ type: t, list: got.map((x, k) => ({ ...x, i: list[k].i })) }); } return res; };
  const H = all(handle, [])[0], X = H && all(placement, [H.type])[0]; if (!X) return null;
  // object numbers count the records from 2 (number 1 is the segment header); Inventor 2026 writes a second set of
  // definitions and handles (same ids), so the first definition class that pairs off completely is taken
  const handles = new Map(H.list.map(h => [h.i + 2, h])), where = new Map(X.list.map(x => [x.occ, x])); if (where.size !== X.list.length) return null;
  for (const C of all(definition, [H.type, X.type])) {
    const ids = C.list.flatMap(c => c.ids); if (ids.length !== X.list.length || new Set(ids).size !== ids.length || !ids.every(id => handles.has(id) && where.has(handles.get(id).occ)) || new Set(ids.map(id => handles.get(id).occ)).size !== ids.length || refs.occurrences.length !== ids.length || refs.occurrences.includes(-1)) continue;
    const out = []; let q = 0, good = true;
    for (const c of C.list) { const fi = refs.occurrences[q]; for (let k = 1; k < c.ids.length; k++) if (refs.occurrences[q + k] !== fi) good = false; q += c.ids.length; const file = refs.files[fi].path;
      for (let k = 0; k < c.ids.length; k++) { const occ = handles.get(c.ids[k]).occ; out.push({ id: occ, name: c.name, file, instance: k + 1, subAssembly: /\.iam$/i.test(file), matrix: where.get(occ).matrix }); } }
    if (good) return { components: out, files: refs.files.map(f => f.path) };
  }
  return null;
}

// ---------------------------------------------------------------- Autodesk Fusion (.f3d / .f3z)
export function parseFusion(bytes) {
  const b = F.toBytes(bytes), WANT = /\.(png|jpe?g|smb|smbh|smt|sab|sat|f3d)$|(^|\/)(Manifest|Properties)\.dat$|(^|\/)Manifest\.json$|(Meta|Bulk|Redirections)Stream\.dat$/i, all = F.unzip(b, WANT), recovered = [], notRecovered = [], extra = {}, meta = { entries: all.length };
  const breps = all.filter(e => /\.(smb|smbh|smt|sab|sat)$/i.test(e.name)), imgs = all.filter(e => e.data && /\.(png|jpe?g)$/i.test(e.name) && sniffImage(e.data, 0)).sort((a, c) => c.size - a.size), nested = all.filter(e => /\.f3d$/i.test(e.name));
  meta.archiveEntries = all.slice(0, 60).map(e => `${e.name} (${e.size} bytes)`); recovered.push(`archive listing: ${all.length} entries`);
  if (imgs.length) { extra.preview = preview(sniffImage(imgs[0].data, 0), `archive entry "${imgs[0].name}"`); recovered.push(describePreview(extra.preview)); }
  const man = all.find(e => e.data && /Manifest\.dat$/i.test(e.name)); if (man) { const s = wideStrings(man.data, 4, 40).concat((latin(man.data, 0, min(man.data.length, 65536)).match(/[\x20-\x7e]{6,}/g) || []).slice(0, 40)); if (s.length) meta.manifestStrings = [...new Set(s)].slice(0, 30); }
  if (nested.length) { meta.components = nested.map(e => e.name); recovered.push(`${nested.length} nested design archive(s) (.f3z assembly): ${nested.slice(0, 8).map(e => e.name).join(', ')}`); }
  // B-rep blobs: a ".smbh" blob holds the current bodies of one component followed by their roll-back history; the ".smb"
  // blobs beside it hold sketch profiles and base-feature bodies, and are the design itself only when no .smbh exists
  // (direct-modelling designs). A .f3z holds one archive per component document.
  const blobs = []; let aside = 0, zstd = 0;
  const collect = (entries, owner) => { zstd += entries.filter(e => !e.data && e.method === 93 && /\.(smb|smbh|smt)$/i.test(e.name)).length; const hist = entries.filter(e => e.data && /\.smbh$/i.test(e.name)), cur = entries.filter(e => e.data && /\.(smb|smt|sab|sat)$/i.test(e.name)), use = hist.length ? hist : cur; if (hist.length) aside += cur.length; for (const e of use) blobs.push({ name: e.name, owner, data: e.data }); };
  collect(all, '');
  for (const n of nested.slice(0, 64)) { if (!n.data) continue; let inner = []; try { inner = F.unzip(n.data, /\.(smb|smbh|smt|sab|sat)$/i); } catch (e) { inner = []; } collect(inner, n.name); }
  const good = blobs.filter(q => findACIS(q.data.subarray(0, min(q.data.length, 64))).some(h => h.offset === 0) || /^\s*\d{3,5}\s+\d+\s+\d+\s+\d+/.test(latin(q.data, 0, 40)));
  if (breps.length) { meta.brepSegments = breps.map(e => `${e.name} (${e.size} bytes)`); recovered.push(`${breps.length} B-rep entr${breps.length === 1 ? 'y' : 'ies'} located: ${breps.slice(0, 6).map(e => e.name.split('/').pop()).join(', ')}`); }
  if (good.length) { extra.embedded = good.map((q, i) => ({ name: `${q.name.split('/').pop()}`, mime: 'application/octet-stream', format: 'acis', bytes: q.data, component: q.owner, unitScale: 10, description: `Autodesk ShapeManager B-rep ${i + 1} of ${good.length}${q.owner ? ` of component archive ${q.owner}` : ''} (${q.name.split('/').pop()}, ${q.data.length} bytes)` })); meta.shapeManager = extra.embedded.map(e => e.description); recovered.push(`${good.length} ShapeManager (ACIS-derived) B-rep stream(s) unpacked intact${aside ? `; ${aside} further .smb blob(s) holding sketch profiles and base-feature bodies left aside` : ''}`); }
  else if (zstd) notRecovered.push(`the ${zstd} B-rep entr${zstd === 1 ? 'y is' : 'ies are'} stored with Zstandard compression (zip method 93, used by recent Fusion releases), which this reader does not inflate: the ShapeManager data was not unpacked`);
  else if (breps.length || nested.length) notRecovered.push('the B-rep entries did not unpack to ShapeManager binary data with a known signature');
  if (zstd && good.length) notRecovered.push(`${zstd} further B-rep entr${zstd === 1 ? 'y is' : 'ies are'} stored with Zstandard compression (zip method 93), which this reader does not inflate: those bodies are missing`);
  else notRecovered.push('no B-rep segment entries were recognised in the archive');
  // component occurrences and which blob holds each component's bodies, from the design stream
  let design = null; try { design = fusionAssembly(all, { unzip: F.unzip }); } catch (e) { design = null; }
  if (design) { const occ = design.components.length - 1, moved = design.components.filter(c => c.matrix && c.matrix.some((x, i) => abs(x - IDENT[i]) > 1e-12)).length; extra.design = design; meta.designStream = { version: design.streamVersion, root: design.root, occurrences: occ, documents: (design.documents || []).length || 1 }; meta.occurrences = design.components.slice(1, 201).map(c => `${c.path.slice(1).join(' / ')}${c.bodies.length ? ` (${c.bodies.length} blob${c.bodies.length === 1 ? '' : 's'})` : ''}${c.external ? ' (external design)' : ''}`);
    recovered.push(`design stream: root component "${design.root}"${occ ? `, ${occ} component occurrence(s), ${moved} of them with a placement other than the identity (4 × 4 matrices in centimetres, in the root frame)` : ', no component occurrences'}; each component's bodies are tied to their B-rep blob`); }
  else if (all.some(e => /(Meta|Bulk)Stream\.dat$/.test(e.name))) notRecovered.push('the design stream is not in the layout this reader knows: component occurrences and their placements were not decoded, so bodies stay in their own stored coordinates');
  notRecovered.push('the design history (timeline, sketches, parameters, joints) is stored in Fusion\'s proprietary serialisation and was not interpreted');
  return { meta, recovered, notRecovered, extra, geometry: null };
}

// ---- Fusion design stream: component occurrences, their placement transforms and the B-rep blob(s) that hold each
// component's bodies.
//
// Layout (all little-endian), as established on Fusion archives written between 2017 and 2025 (stream versions
// 0x08000EDB to 0x08005404):
//   FusionAssetName[Active]/Design1/MetaStream.dat
//     astr "Design", u32, wstr GUID, u32 stream version, u32, astr "Design", u32 × 3, u32 nTypes,
//     nTypes × { astr class GUID, astr base-class GUID, u32 class version, astr module, u32 n, n × u64 object id },
//     u32 nRoots, nRoots × u64, u32 n1, n1 × { u64 object id, u64 offset of the object in BulkStream.dat },
//     u32 n2, n2 × { u64 object id, u64 offset of the object's second (base-class) chunk }, …
//     (astr = u32 length + 8-bit characters, wstr = u32 length + UTF-16LE characters)
//   FusionAssetName[Active]/Design1/BulkStream.dat
//     objects in id order; each chunk starts with astr class number (decimal text, = 256 + index in the type table)
//     and the u64 object id. A reference to another object is u8 flag (0 = null, nothing follows), u64 id,
//     [astr class GUID of the target — streams written before about 2019 only], u8 k, k links into another design
//     document (external reference; newer streams: u64 key, u8 0, astr class GUID; older: u32 0, wstr GUID, u8 0,
//     wstr role key), u8 0. So a plain reference is 11 bytes in newer streams.
//   Classes used (GUIDs are the same in every file; the class numbers are not):
//     CE2913AA-…  occurrence ("instance"): … ref component … 4 × 4 float64 matrix, row-major, acting on column
//                 vectors, translation in the fourth column, centimetres (omitted when the preceding flag byte is 1 =
//                 identity) … wstr name, [wstr role key of an external design] … ref occurrence list
//     904F885E-…  occurrence list of one component: … refs to its occurrences …, ref owning component
//     E03784ED-… + 42054630-…  component (two chunks); its name is the last wstr of the second chunk
//     549FBB80-…  occurrence path (absent from 2017 streams): zero bytes, u32 n, n refs to the occurrences from the root
//                 down, …, 4 × 4 float64 matrix (same conventions) as the last 128 bytes unless it is the identity
//     CD57BC48-…  body container: …, ref owning component
//     74CA4562-… (history blob, ".smbh") / 40DCC3F4-… (plain blob, ".smb"): u32 n, n × { u64 key, u64 body object id },
//                 ref body container, wstr "BREP.<guid>.smbh" = the name of the archive entry under Breps.BlobParts/
//
// Result: { units: 'cm', streamVersion, root, components: [{ name, component, path, depth, instance, bodies, bodyKeys,
//           bodyObjects, matrix, placement, occurrenceMatrix?, external? }], notes }. components[0] is the root component
//           (identity); the others are its occurrences in tree order, sub-assemblies expanded (path = names from the root).
//   external { key, role } on an occurrence of another design document (.f3z); it has no bodies of its own here.
//   matrix   12 numbers, 3 × 4 row-major, acting on column vectors, translation (elements 3, 7, 11) in CENTIMETRES — the
//            unit of the ShapeManager data in the blobs. It maps the blob's coordinates to the root design's frame; it is
//            a stored matrix, NOT composed with the parent occurrence (the file stores root-frame matrices).
//   placement 'path-record' (matrix taken from the occurrence-path record, class 549FBB80-…) | 'occurrence-record' (no
//            path record for this occurrence; the matrix of the occurrence record is used) | 'unknown' (matrix: null —
//            no path record and the occurrence lies in a sub-assembly that is placed more than once).
//            occurrenceMatrix is present when the two records disagree or the placement is unknown.
//   bodies   names of the archive entries (final path segment, e.g. "BREP.<guid>.smbh") that hold this component's
//            bodies: the history blob if the component has one, otherwise its plain blobs.
//   bodyKeys { blob name: [key, …] } — for a ".smbh" blob the key is the entity id number written in each body record
//            (the integer that follows the attribute pointer). In a ".smb" blob the body records carry no id (-1) and
//            the keys run 0 … n-1 for n body records: they are taken to be the 0-based order of the body records. A body record whose key is not listed is not a body of the design (the history
//            blob of a PCB board carries two unlisted sheets, for instance) and is to be left out.
//   bodyObjects { blob name: [object id, …] } — the design-stream object of each listed body, parallel to bodyKeys.
// Returns null when the streams are missing or their layout is not the one described above.

const FUSION_CLASS = {
  instance: 'CE2913AA-CFE0-4F04-9102-24424ED3BCFA', list: '904F885E-D859-4755-A7EF-6EBED8732FAD',
  component: 'E03784ED-5E19-4E14-B9F2-3B07017018CD', componentBase: '42054630-20A0-40E1-B969-CFE9E742F5C9',
  bodies: 'CD57BC48-50EC-47DC-975A-FB6DEA72F4DA', path: '549FBB80-B890-473E-A5C0-415D3D9BF4E6', blobHistory: '74CA4562-59F8-4C97-8AD1-F8297C21F9AA', blob: '40DCC3F41D-E94F-B3B1-3E95-9B6ABE572C',
};
class FusionLayout extends Error { }
const need = (c, why) => { if (!c) throw new FusionLayout(why); };

function fusionReader(b) {
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), n = b.length;
  const r = {
    b, n,
    u8: p => { need(p >= 0 && p < n, 'read past the end'); return b[p]; },
    u16: p => { need(p >= 0 && p + 2 <= n, 'read past the end'); return dv.getUint16(p, true); },
    u32: p => { need(p >= 0 && p + 4 <= n, 'read past the end'); return dv.getUint32(p, true); },
    u64: p => { need(p >= 0 && p + 8 <= n, 'read past the end'); const hi = dv.getUint32(p + 4, true); need(hi < 0x200000, 'a 64-bit value is out of range'); return hi * 4294967296 + dv.getUint32(p, true); },
    f64: p => { need(p >= 0 && p + 8 <= n, 'read past the end'); return dv.getFloat64(p, true); },
    astr: (p, max = 4096) => { const l = r.u32(p); need(l <= max && p + 4 + l <= n, 'a text field is too long'); let s = ''; for (let i = 0; i < l; i++) s += String.fromCharCode(b[p + 4 + i]); return { s, next: p + 4 + l }; },
    wstr: (p, max = 4096) => { const l = r.u32(p); need(l <= max && p + 4 + 2 * l <= n, 'a text field is too long'); let s = ''; for (let i = 0; i < l; i++) s += String.fromCharCode(dv.getUint16(p + 4 + 2 * i, true)); return { s, next: p + 4 + 2 * l }; },
  };
  return r;
}
const isGuid = s => /^[0-9A-Fa-f-]{36,37}$/.test(s);

// MetaStream.dat → { version, types, typeOf: Map(GUID → type), first: Map(id → offset), second: Map(id → offset), order }
function fusionMeta(bytes) {
  const r = fusionReader(bytes); let p = 0;
  const a0 = r.astr(p, 64); need(a0.s === 'Design', 'the stream is not a design stream'); p = a0.next + 4;
  const g0 = r.wstr(p, 64); p = g0.next; const version = r.u32(p); p += 8;
  const a1 = r.astr(p, 64); p = a1.next + 12; const nTypes = r.u32(p); p += 4; need(nTypes > 0 && nTypes < 65536, 'the type table is implausible');
  const types = [], typeOf = new Map();
  for (let i = 0; i < nTypes; i++) {
    const g = r.astr(p, 64); p = g.next; const base = r.astr(p, 64); p = base.next; const ver = r.u32(p); p += 4; const mod = r.astr(p, 256); p = mod.next; const k = r.u32(p); p += 4;
    need(isGuid(g.s) && (base.s === '' || isGuid(base.s)) && k < 1e7 && p + 8 * k <= r.n, 'a type table entry is malformed');
    const ids = new Array(k); for (let j = 0; j < k; j++) { ids[j] = r.u64(p); p += 8; }
    const t = { index: i, guid: g.s.toUpperCase(), base: base.s.toUpperCase(), version: ver, module: mod.s, ids }; types.push(t); typeOf.set(t.guid, t);
  }
  const nRoots = r.u32(p); p += 4; need(nRoots < 1e6, 'the root list is implausible'); p += 8 * nRoots;
  const table = () => { const k = r.u32(p); p += 4; need(k < 5e7 && p + 16 * k <= r.n, 'an offset table is malformed'); const m = new Map(); let last = -1; for (let j = 0; j < k; j++) { const id = r.u64(p), off = r.u64(p + 8); p += 16; need(off > last, 'an offset table is not ascending'); last = off; m.set(id, off); } return m; };
  const first = table(), second = p + 4 <= r.n ? table() : new Map();
  return { version, types, typeOf, first, second, order: [...first].sort((x, y) => x[1] - y[1]) };
}

// fusionDesign(metaBytes, bulkBytes, opts) → result | null (see the head of the file). opts.onReject(reason) is called
// with the reason when the layout is not recognised.
export function fusionDesign(metaBytes, bulkBytes, opts = {}) {
  try { return fusionDecode(F.toBytes(metaBytes), F.toBytes(bulkBytes)); } catch (e) { if (e instanceof FusionLayout || e instanceof RangeError) { if (typeof opts.onReject === 'function') opts.onReject(e.message); return null; } throw e; }
}

function fusionDecode(mb, bb) {
  const meta = fusionMeta(mb), r = fusionReader(bb), notes = [];
  const T = k => meta.typeOf.get(FUSION_CLASS[k]); for (const k of ['instance', 'list', 'component', 'componentBase']) need(T(k), `class ${FUSION_CLASS[k]} is not in the type table`);
  // extent of an object: [start of its first chunk, start of the next object); second chunk start if it has one
  const nextOff = new Map(); meta.order.forEach(([id], i) => nextOff.set(id, i + 1 < meta.order.length ? meta.order[i + 1][1] : r.n));
  // chunk header at p for object id of class t: astr class number, u64 id → offset of the chunk's fields
  const head = (p, id, t) => { const c = r.astr(p, 8); need(c.s === String(256 + t.index), `object ${id}: the chunk is not of the expected class`); need(r.u64(c.next) === id, `object ${id}: the chunk carries another id`); return c.next + 8; };
  const chunk = (id, t) => { const s = meta.first.get(id); need(s !== undefined, `object ${id} is not in the offset table`); const mid = meta.second.get(id), end = nextOff.get(id); return { p: head(s, id, t), end: mid !== undefined ? mid : end }; };
  // a reference read forwards → {id | null, next}
  const ref = p => {
    const f = r.u8(p); if (f === 0) return { id: null, next: p + 1 }; need(f === 1, 'a reference flag is neither 0 nor 1'); const id = r.u64(p + 1); let q = p + 9;
    const l = q + 4 <= r.n ? r.u32(q) : 0;
    if ((l === 36 || l === 37) && isGuid(r.astr(q, 64).s)) {                         // before about 2019: the class GUID of the target follows
      q += 4 + l; const k = r.u8(q++); need(k <= 1, 'a reference has an implausible number of external links'); let role = '';
      if (k) { need(r.u32(q) === 0, 'an external link is malformed'); const g = r.wstr(q + 4, 64); need(isGuid(g.s) && r.u8(g.next) === 0, 'an external link is malformed'); const w = r.wstr(g.next + 1, 1024); role = w.s; q = w.next; }
      need(r.u8(q) === 0, 'a reference is not closed by a zero byte'); return k ? { id, ext: role, role, next: q + 1 } : { id, next: q + 1 };
    }
    // k further links into another design document (an external reference): each a 64-bit key, a zero byte and the class GUID
    const k = r.u8(q++); need(k <= 4, 'a reference has an implausible number of external links'); const ext = [];
    for (let i = 0; i < k; i++) { const lo = r.u32(q), hi = r.u32(q + 4); need(r.u8(q + 8) === 0, 'an external link is malformed'); const g = r.astr(q + 9, 64); need(isGuid(g.s), 'an external link is malformed'); ext.push(hi.toString(16).padStart(8, '0') + lo.toString(16).padStart(8, '0')); q = g.next; }
    need(r.u8(q) === 0, 'a reference is not closed by a zero byte'); return ext.length ? { id, ext: ext.join('/'), next: q + 1 } : { id, next: q + 1 };
  };
  // a non-null reference that ends at `end`, read backwards → {id, start}
  const refBack = end => {
    need(end >= 11 && r.u16(end - 2) === 0, 'a trailing reference is missing');
    for (const l of [36, 37]) { const s = end - 2 - l - 4 - 9; if (s >= 0 && r.u32(s + 9) === l && r.u8(s) === 1 && isGuid(r.astr(s + 9, 64).s)) return { id: r.u64(s + 1), start: s }; }
    const s = end - 11; need(r.u8(s) === 1, 'a trailing reference is missing'); return { id: r.u64(s + 1), start: s };
  };

  // ---- occurrences
  const inst = new Map();
  for (const id of T('instance').ids) {
    const c = chunk(id, T('instance')); let p = c.p;
    const tag = r.astr(p, 64); p = tag.next; need(r.u8(p) <= 1 && r.u32(p + 1) <= 255, `occurrence ${id}: unexpected header`); p += 5;
    const comp = ref(p); need(comp.id !== null, `occurrence ${id}: no component`); p = comp.next;
    // counters (and, from about 2023, a list of property-group names, a self reference and a short label); the last u32
    // is the 1-based instance number; then two GUIDs as wstr
    let g1 = null, g2 = null; for (let k = 4; k <= 400 && p + k + 160 <= c.end && !g2; k++) { if (r.u32(p + k) !== 36) continue; const x = r.wstr(p + k, 64); if (!isGuid(x.s) || r.u32(x.next) !== 36) continue; const y = r.wstr(x.next, 64); if (isGuid(y.s)) { g1 = x; g2 = y; p += k; } }
    need(g2, `occurrence ${id}: the two GUIDs are missing`); const index = Math.max(0, r.u32(p - 4) - 1); p = g2.next;
    if (r.u8(p) === 0) p++;                                                   // one zero byte in streams written from about 2019
    const comp2 = ref(p); need(comp2.id === comp.id && comp2.ext === comp.ext, `occurrence ${id}: the two component references differ`); p = comp2.next;
    const f = r.u8(p++); need(f <= 1, `occurrence ${id}: unexpected transform flag`);
    let matrix = IDENT.slice();
    if (f === 0) {
      const m = []; for (let i = 0; i < 16; i++) m.push(r.f64(p + 8 * i)); p += 128;
      need(m.every(Number.isFinite) && Math.abs(m[12]) + Math.abs(m[13]) + Math.abs(m[14]) < 1e-9 && Math.abs(m[15] - 1) < 1e-9, `occurrence ${id}: the matrix is not affine`);
      // rows of the 3 × 3 block must be orthonormal (a rigid motion, possibly mirrored)
      for (let i = 0; i < 3; i++) for (let j = i; j < 3; j++) { const d = m[4 * i] * m[4 * j] + m[4 * i + 1] * m[4 * j + 1] + m[4 * i + 2] * m[4 * j + 2]; need(Math.abs(d - (i === j ? 1 : 0)) < 1e-6, `occurrence ${id}: the matrix is not a rigid motion`); }
      matrix = m.slice(0, 12);
    }
    const name = r.wstr(p, 1024); p = name.next;
    // for an occurrence of an external design: the role key under which RedirectionsStream.dat names the referenced archive
    let role = ''; if (c.end - p >= 15) { try { const w = r.wstr(p, 512); if (w.s.length >= 36 && isGuid(w.s.slice(0, 36))) { role = w.s; p = w.next; } } catch (e) { if (!(e instanceof FusionLayout)) throw e; } }
    let list = null; if (c.end - p >= 11) { const l = refBack(c.end); need(l.start >= p && l.start - p <= 8, `occurrence ${id}: unexpected bytes after the name`); list = l.id; }
    inst.set(id, { id, tag: tag.s, component: comp.id, ext: comp.ext || '', role: role || comp.role || '', index, matrix, name: name.s, list });
  }
  // ---- occurrence lists → owning component
  const listOwner = new Map(), listOrder = new Map();
  for (const id of T('list').ids) {
    const c = chunk(id, T('list')); listOwner.set(id, refBack(c.end).id);
    try { let p = c.p + 5; const k = r.u32(p); p += 6; const ids = []; if (k === 1) { const cnt = r.u32(p); p += 4; need(cnt < 1e6, ''); for (let i = 0; i < cnt; i++) { const x = ref(p); ids.push(x.id); p = x.next; } } if (k <= 1 && ref(p).next === c.end) listOrder.set(id, ids); } catch (e) { if (!(e instanceof FusionLayout)) throw e; }
  }
  // ---- components: name
  const comps = new Map();
  for (const id of T('component').ids) {
    const s = meta.first.get(id), mid = meta.second.get(id); need(s !== undefined, `component ${id} is not in the offset table`); head(s, id, T('component'));
    let name = ''; const end = nextOff.get(id);
    if (mid !== undefined) {
      head(mid, id, T('componentBase'));
      // tail: wstr GUID, wstr name, u8 0, reference to the root occurrence (null for every component but the root)
      const tails = [2]; try { tails.push(end - refBack(end).start + 1); } catch (e) { if (!(e instanceof FusionLayout)) throw e; }
      find: for (const t of tails) for (let l = 0; l <= 512; l++) { const q = end - t - 2 * l - 4; if (q - 76 < mid) break; if (r.u32(q) !== l || r.u32(q - 76) !== 36) continue; if (!isGuid(r.wstr(q - 76, 64).s)) continue; name = r.wstr(q, 1024).s; break find; }
    }
    comps.set(id, { id, name, blobs: [], children: [] });
  }
  // ---- B-rep blobs → body container → component
  const bodyOwner = new Map(); if (T('bodies')) for (const id of T('bodies').ids) { const c = chunk(id, T('bodies')); try { bodyOwner.set(id, refBack(c.end).id); } catch (e) { if (!(e instanceof FusionLayout)) throw e; } }
  for (const k of ['blobHistory', 'blob']) {
    if (!T(k)) continue;
    for (const id of T(k).ids) {
      const c = chunk(id, T(k)); let nameAt = -1, name = '';
      for (let l = 10; l <= 200 && nameAt < 0; l++) { const q = c.end - 2 * l - 4; if (q < c.p) break; if (r.u32(q) === l) { const w = r.wstr(q, 256).s; if (/^BREP\.[0-9A-Fa-f-]{36}\.smbh?$/.test(w)) { nameAt = q; name = w; } } }
      if (nameAt < 0) continue;                                                     // a blob record without a stored blob
      // the reference to the body container ends one byte before the name
      let cont = null; for (const gap of [1, 0, 2]) { try { const c2 = refBack(nameAt - gap); if (bodyOwner.has(c2.id)) { cont = c2; break; } } catch (e) { if (!(e instanceof FusionLayout)) throw e; } }
      if (!cont) { notes.push(`blob ${name}: its body container reference was not understood`); continue; }
      const owner = bodyOwner.get(cont.id); if (!comps.has(owner)) { notes.push(`blob ${name}: its owning component was not found`); continue; }
      const keys = [], objs = []; for (let n = Math.floor((cont.start - c.p - 4) / 16); n >= 0; n--) { const s = cont.start - 16 * n - 4; if (r.u32(s) !== n) continue; const kv = [], ov = []; let ok = true; for (let i = 0; i < n && ok; i++) { const hiK = r.u32(s + 4 + 16 * i + 4), hiV = r.u32(s + 4 + 16 * i + 12); if (hiK || hiV) { ok = false; break; } const key = r.u32(s + 4 + 16 * i), v = r.u32(s + 4 + 16 * i + 8); if (!meta.first.has(v)) ok = false; kv.push(key); ov.push(v); } if (ok && new Set(kv).size === n) { keys.push(...kv); objs.push(...ov); break; } }
      comps.get(owner).blobs.push({ name, history: k === 'blobHistory', keys, objs });
    }
  }
  // ---- occurrence-path records (streams from about 2018): zero bytes, u32 n, n references to the occurrences along the
  // path from the root, …, and — unless it is the identity — a 4 × 4 matrix as the last 128 bytes
  const rigid = m => m.every(Number.isFinite) && Math.abs(m[12]) + Math.abs(m[13]) + Math.abs(m[14]) < 1e-9 && Math.abs(m[15] - 1) < 1e-9 && [0, 1, 2].every(i => [0, 1, 2].every(j => Math.abs(m[4 * i] * m[4 * j] + m[4 * i + 1] * m[4 * j + 1] + m[4 * i + 2] * m[4 * j + 2] - (i === j ? 1 : 0)) < 1e-6));
  const pathRec = new Map();
  if (T('path')) for (const id of T('path').ids) {
    try {
      const c = chunk(id, T('path')); let p = c.p; while (p < c.end && p - c.p < 24 && r.u8(p) === 0) p++; const n = r.u32(p); if (n < 1 || n > 32) continue; p += 4;
      const ids = []; let local = true; for (let i = 0; i < n; i++) { const x = ref(p); p = x.next; if (x.id === null || x.ext || !inst.has(x.id)) local = false; ids.push(x.id); } if (!local) continue;
      let M = IDENT.slice(); if (c.end - p >= 128) { const m = []; for (let i = 0; i < 16; i++) m.push(r.f64(c.end - 128 + 8 * i)); if (rigid(m)) M = m.slice(0, 12); else continue; }
      const key = ids.join('/'); if (!pathRec.has(key)) pathRec.set(key, M);
    } catch (e) { if (!(e instanceof FusionLayout)) throw e; }
  }
  const sameM = (A, B) => A.every((v, i) => Math.abs(v - B[i]) < 1e-7);
  // The path record gives the placement of an occurrence in the root frame. The occurrence record stores the same matrix
  // unless a later feature (a joint) has moved things: then only the path records are mutually consistent, and they are
  // used. Streams without path records (2017) have the root-frame matrix in the occurrence record.
  const check = (ids, M, o, deep) => { const pm = pathRec.get(ids.join('/')); if (!pm && deep) { o.placement = 'unknown'; o.occurrenceMatrix = M.slice(); o.matrix = null; shared++; } else if (!pm) { o.placement = 'occurrence-record'; unconfirmed++; } else { o.placement = 'path-record'; if (!sameM(pm, M)) { o.occurrenceMatrix = M.slice(); conflicts++; } o.matrix = pm.slice(); } return o; };
  let conflicts = 0, unconfirmed = 0, shared = 0;
  // ---- tree
  let root = null; const loose = [];
  for (const i of inst.values()) {
    need(i.ext || comps.has(i.component), `occurrence ${i.id}: its component ${i.component} is not a component object`);
    const parent = i.list === null ? undefined : listOwner.get(i.list);
    if (i.tag === 'rootInstance' || i.list === null) { need(root === null, 'the design has two root occurrences'); root = i; continue; }
    if (parent === undefined || !comps.has(parent)) { loose.push(i); continue; }
    i.parent = parent; comps.get(parent).children.push(i);
  }
  need(root !== null, 'the design has no root occurrence');
  if (loose.length) notes.push(`${loose.length} occurrence(s) are not attached to a component and were left out`);
  for (const c of comps.values()) { const ord = [...listOwner].filter(([, o]) => o === c.id).flatMap(([l]) => listOrder.get(l) || []); c.children.sort((x, y) => { const a = ord.indexOf(x.id), b = ord.indexOf(y.id); return (a < 0 ? 1e9 : a) - (b < 0 ? 1e9 : b) || x.id - y.id; }); }
  const useCount = new Map(); for (const i of inst.values()) if (i !== root && !i.ext) useCount.set(i.component, (useCount.get(i.component) || 0) + 1);
  const out = []; let external = 0;
  const bodiesOf = c => { const h = c.blobs.filter(q => q.history && q.keys.length), use = h.length ? h : c.blobs.filter(q => q.keys.length); return use; };
  const walk = (c, path, depth, matrix, instance, viaShared, guard, ids) => {
    need(depth < 64 && !guard.has(c.id), 'the occurrence tree is cyclic');
    const use = bodiesOf(c), bodyKeys = {}, bodyObjects = {}; for (const q of use) { bodyKeys[q.name] = q.keys.slice(); bodyObjects[q.name] = q.objs.slice(); }
    const o = { name: path[path.length - 1], component: c.name, path, depth, instance, bodies: use.map(q => q.name), bodyKeys, bodyObjects, matrix }; out.push(ids.length ? check(ids, matrix, o, viaShared) : o);
    guard.add(c.id);
    for (const i of c.children) {
      const deep = viaShared || (depth > 0 && (useCount.get(c.id) || 0) > 1);
      // an occurrence of another design document: its bodies are in that document's archive (see fusionAssembly)
      if (i.ext) { external++; out.push(check([...ids, i.id], i.matrix, { name: i.name || `external:${i.index + 1}`, component: '', path: [...path, i.name || `external:${i.index + 1}`], depth: depth + 1, instance: i.index + 1, bodies: [], bodyKeys: {}, bodyObjects: {}, matrix: i.matrix.slice(), external: { key: i.ext, role: i.role } }, deep)); continue; }
      const k = comps.get(i.component), nm = i.name || `${k.name}:${i.index + 1}`;
      walk(k, [...path, nm], depth + 1, i.matrix.slice(), i.index + 1, deep, guard, [...ids, i.id]);
    }
    guard.delete(c.id);
  };
  walk(comps.get(root.component), [comps.get(root.component).name || 'root'], 0, IDENT.slice(), 1, false, new Set(), []);
  if (conflicts) notes.push(`${conflicts} occurrence(s) have a path record whose matrix differs from the one in the occurrence record (kept as occurrenceMatrix): the path-record matrix is used`);
  if (unconfirmed && pathRec.size) notes.push(`${unconfirmed} occurrence(s) have no path record: the matrix of the occurrence record is used for them`);
  const risky = conflicts ? out.filter(o => o.placement === 'occurrence-record' && o.bodies.length).length : 0;
  if (risky) notes.push(`${risky} of them hold bodies while other occurrences of this design were moved after insertion: their placement relative to those is not established`);
  if (shared) notes.push(`${shared} occurrence(s) lie inside a sub-assembly that is placed more than once and have no path record: the occurrence record holds one root-frame matrix only, so their placement is not established (placement: 'unknown', matrix: null)`);
  if (external) notes.push(`${external} occurrence(s) refer to another design document (external reference)`);
  notes.push('matrices are the stored ones: row-major, acting on column vectors, translation in centimetres, expressed in the root design frame (not relative to the parent occurrence)');
  return { units: 'cm', streamVersion: meta.version, root: comps.get(root.component).name, components: out, notes };
}

const textOf = d => { let t = ''; for (let i = 0; i < d.length; i++) t += String.fromCharCode(d[i]); return t; };
const jsonOf = e => { if (!e || !e.data) return null; try { return JSON.parse(textOf(e.data)); } catch (x) { return null; } };
const FUSION_STREAM = n => new RegExp('(^|/)FusionAssetName\\[Active\\]/Design1/' + n + '\\.dat$');

// One design document: entries of one .f3d archive → result of fusionDesign with every component given
// document: docName, and the blob names checked against the archive (result.missing).
function fusionDocument(entries, docName, opts) {
  const m = entries.find(e => e && e.data && FUSION_STREAM('MetaStream').test(e.name)), b = entries.find(e => e && e.data && FUSION_STREAM('BulkStream').test(e.name)); if (!m || !b) return null;
  const res = fusionDesign(m.data, b.data, opts); if (!res) return null;
  const have = new Set(entries.filter(e => e && /(^|\/)Breps\.BlobParts\/[^/]+$/.test(e.name)).map(e => e.name.split('/').pop())), missing = new Set();
  for (const c of res.components) { c.document = docName; for (const n of c.bodies) if (!have.has(n)) missing.add(n); }
  res.missing = [...missing]; if (missing.size) res.notes.push(`${missing.size} blob(s) named by the design stream are not in the archive: ${[...missing].slice(0, 6).join(', ')}`);
  // external references: role key → file name of the referenced archive
  res.redirect = new Map(); const red = jsonOf(entries.find(e => e && /(^|\/)RedirectionsStream\.dat$/.test(e.name) && !/Active\]/.test(e.name)));
  if (red && Array.isArray(red.references)) for (const x of red.references) { if (!x || typeof x.relativePath !== 'string') continue; for (const p of Array.isArray(x.properties) ? x.properties : []) { const v = p && p.neutronRole && p.neutronRole.value; if (typeof v === 'string' && v) res.redirect.set(v, x.relativePath); } }
  return res;
}

// fusionAssembly(entries, opts) — entries: the unzipped archive, [{name, data: Uint8Array | null}].
//   .f3d: the entries must include FusionAssetName[Active]/Design1/MetaStream.dat and BulkStream.dat, inflated.
//   .f3z: the entries are Manifest.json ({"root": "<guid>.f3d"}) and one .f3d per design document; opts.unzip(bytes,
//         regexp) → entries is then needed to open them (formats.js unzip). Occurrences of an external design are
//         resolved through RedirectionsStream.dat of the referring document (role key → archive name) and the
//         referenced document's own occurrences are appended with matrix = (matrix of the referring occurrence) ·
//         (their matrix), document = the name of the nested .f3d entry.
// Returns the result of fusionDesign (each component with document: '' or the nested .f3d name), or null.
export function fusionAssembly(entries, opts = {}) {
  if (!Array.isArray(entries)) return null;
  try {
    if (entries.some(e => e && e.data && FUSION_STREAM('MetaStream').test(e.name))) return fusionDocument(entries, '', opts);
    const nested = entries.filter(e => e && e.data && /\.f3d$/i.test(e.name)); if (!nested.length || typeof opts.unzip !== 'function') return null;
    const man = jsonOf(entries.find(e => e && /(^|\/)Manifest\.json$/.test(e.name))); let rootName = man && typeof man.root === 'string' ? man.root : nested.length === 1 ? nested[0].name : null; if (!rootName) return null;
    const WANT = /(Meta|Bulk)Stream\.dat$|RedirectionsStream\.dat$|Breps\.BlobParts\//, docs = new Map();
    const open = name => { if (docs.has(name)) return docs.get(name); const e = nested.find(x => x.name === name || x.name.split('/').pop() === name); let d = null; if (e) { let inner = null; try { inner = opts.unzip(e.data, WANT); } catch (x) { inner = null; } if (Array.isArray(inner)) d = fusionDocument(inner, e.name, opts); } docs.set(name, d); return d; };
    const root = open(rootName); if (!root) return null;
    const out = [], notes = root.notes.slice(), missing = new Set(root.missing); let unresolved = 0;
    const add = (doc, M, prefix, depth0, stack) => {
      for (const c of doc.components) {
        if (prefix && c.depth === 0 && !c.bodies.length) continue;                    // the referenced document's root is the occurrence already listed
        const k = { ...c, path: prefix ? [...prefix, ...c.path.slice(1)] : c.path, depth: c.depth + depth0, matrix: !c.matrix || M === undefined ? null : M ? matMul(M, c.matrix) : c.matrix.slice() };
        if (!c.external) { out.push(k); continue; }
        const file = doc.redirect.get(c.external.role) || (doc.redirect.size === 1 && !c.external.role ? [...doc.redirect.values()][0] : null), child = file && !stack.includes(file) ? open(file) : null;
        k.external = { ...c.external, document: file || '' };
        if (child && /^external:\d+$/.test(c.name)) { k.name = `${child.root}:${c.instance}`; k.component = child.root; k.path = [...k.path.slice(0, -1), k.name]; }
        out.push(k); if (!child) { unresolved++; continue; }
        for (const n of child.missing) missing.add(n);
        add(child, k.matrix === null ? undefined : k.matrix, k.path, k.depth, [...stack, file]);
      }
    };
    add(root, null, null, 0, [rootName]);
    if (unresolved) notes.push(`${unresolved} occurrence(s) of an external design could not be resolved to a design document in the archive: their bodies are not placed`);
    return { units: 'cm', streamVersion: root.streamVersion, root: root.root, rootDocument: rootName, documents: [...docs].filter(([, d]) => d).map(([n]) => n), components: out, notes, missing: [...missing] };
  } catch (e) { if (e instanceof FusionLayout || e instanceof RangeError) return null; throw e; }
}

// the 3 × 4 matrix in millimetres (rotation unchanged, translation × 10), for bodies already scaled from cm to mm
export const matrixToMm = m => m.map((x, i) => (i % 4 === 3 ? x * 10 : x));

// ---------------------------------------------------------------- FreeCAD (.FCStd)
const quatMat = (q, p) => { const [x, y, z, w] = q, n = hypot(x, y, z, w) || 1, a = x / n, b = y / n, c = z / n, d = w / n; return [1 - 2 * (b * b + c * c), 2 * (a * b - c * d), 2 * (a * c + b * d), p[0], 2 * (a * b + c * d), 1 - 2 * (a * a + c * c), 2 * (b * c - a * d), p[1], 2 * (a * c - b * d), 2 * (b * c + a * d), 1 - 2 * (a * a + b * b), p[2]]; };
const matMul = (A, B) => { const o = new Array(12); for (let r = 0; r < 3; r++) { for (let c = 0; c < 3; c++) o[4 * r + c] = A[4 * r] * B[c] + A[4 * r + 1] * B[4 + c] + A[4 * r + 2] * B[8 + c]; o[4 * r + 3] = A[4 * r] * B[3] + A[4 * r + 1] * B[7] + A[4 * r + 2] * B[11] + A[4 * r + 3]; } return o; };
const IDENT = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0];
// readFCStd(bytes) → {meta, shapes: [{name, label, type, file, brep: Uint8Array, visible, colour, transform: 3 × 4 row-major | null}], notes, extra}.
// Each shape is the OpenCASCADE BREP that FreeCAD saved for an object; `transform` is the placement of the enclosing
// App::Part containers (the object's own placement is already inside its BREP).
export function readFCStd(bytes) {
  const entries = F.unzip(F.toBytes(bytes), /(^|\/)(Document|GuiDocument)\.xml$|\.brp$|\.brep$|thumbnails\/.*\.png$/i), doc = entries.find(e => e.name === 'Document.xml' && e.data); if (!doc) throw new Error('FreeCAD: the archive holds no Document.xml');
  const t = F.toText(doc.data), gui = entries.find(e => e.name === 'GuiDocument.xml' && e.data), gt = gui ? F.toText(gui.data) : '', notes = [], extra = {};
  const meta = { programVersion: (/ProgramVersion="([^"]+)"/.exec(t) || [])[1], archiveEntries: entries.slice(0, 60).map(e => `${e.name} (${e.size} bytes)`) };
  const docProp = k => { const m = new RegExp('<Property name="' + k + '"[^>]*>\\s*<String value="([^"]*)"').exec(t.slice(0, t.indexOf('<Objects') > 0 ? t.indexOf('<Objects') : 20000)); return m ? xmlText(m[1]) : ''; };
  for (const k of ['Label', 'CreatedBy', 'CreationDate', 'LastModifiedBy', 'LastModifiedDate', 'Company', 'Comment', 'License']) { const x = docProp(k); if (x) meta[k[0].toLowerCase() + k.slice(1)] = x; }
  const types = new Map(); for (const m of t.matchAll(/<Object type="([^"]+)" name="([^"]+)"/g)) types.set(m[2], m[1]); meta.objects = [...types].slice(0, 200).map(([n, ty]) => `${n} : ${ty}`);
  const od = t.indexOf('<ObjectData'), body = od >= 0 ? t.slice(od) : '', objs = new Map();
  for (const m of body.matchAll(/<Object name="([^"]+)"[^>]*>([\s\S]*?)<\/Object>/g)) {
    const p = m[2], prop = k => { const i = p.indexOf(`<Property name="${k}"`); if (i < 0) return ''; const j = p.indexOf('</Property>', i); return p.slice(i, j < 0 ? i + 2000 : j); };
    const o = { name: m[1], type: types.get(m[1]) || '', label: xmlText((/<String value="([^"]*)"/.exec(prop('Label')) || [])[1] || m[1]), file: (/<Part\b[^>]*\bfile="([^"]+)"/.exec(prop('Shape')) || [])[1] || '', group: [...prop('Group').matchAll(/<Link value="([^"]+)"/g)].map(x => x[1]), placement: null };
    const pl = /<PropertyPlacement\b([^>]*)\/>/.exec(prop('Placement')); if (pl) { const a = k => +((new RegExp('\\b' + k + '="([^"]+)"').exec(pl[1]) || [])[1] || 0); o.placement = quatMat([a('Q0'), a('Q1'), a('Q2'), a('Q3')], [a('Px'), a('Py'), a('Pz')]); }
    objs.set(o.name, o);
  }
  const vis = new Map(), col = new Map();
  for (const m of gt.matchAll(/<ViewProvider name="([^"]+)"[^>]*>([\s\S]*?)<\/ViewProvider>/g)) { const v = /<Property name="Visibility"[^>]*>\s*<Bool value="(\w+)"/.exec(m[2]); if (v) vis.set(m[1], v[1] === 'true'); const c = /<Property name="ShapeColor"[^>]*>\s*<PropertyColor value="(\d+)"/.exec(m[2]); if (c) { const u = +c[1] >>> 0; col.set(m[1], [(u >>> 24) / 255, ((u >>> 16) & 255) / 255, ((u >>> 8) & 255) / 255]); } }
  const parent = new Map(), swallowed = new Set();
  for (const o of objs.values()) { if (!o.group.length) continue; if (o.file) for (const g of o.group) swallowed.add(g); else if (/^App::(Part|DocumentObjectGroup)$/.test(o.type) || /Part$|Group$/.test(o.type)) for (const g of o.group) if (!parent.has(g)) parent.set(g, o.name); }
  const chain = n => { let M = null, hidden = false; for (let p = parent.get(n), g = 0; p && g < 64; p = parent.get(p), g++) { const o = objs.get(p); if (!o) break; if (vis.get(p) === false) hidden = true; if (o.placement && /^App::Part$/.test(o.type)) M = M ? matMul(o.placement, M) : o.placement; } return { M, hidden }; };
  const data = new Map(entries.filter(e => e.data && /\.(brp|brep)$/i.test(e.name)).map(e => [e.name, e.data])), shapes = []; let hidden = 0, swall = 0, missing = 0;
  for (const o of objs.values()) {
    if (!o.file) continue; const brep = data.get(o.file); if (!brep || brep.length < 40) { missing++; continue; }
    if (/Sketch|Datum|^App::(Origin|Plane|Line)|^PartDesign::(Plane|Line|Point|CoordinateSystem|ShapeBinder)/.test(o.type)) continue;
    const c = chain(o.name); if (swallowed.has(o.name)) { swall++; continue; } if (gui && (vis.get(o.name) === false || c.hidden)) { hidden++; continue; }
    shapes.push({ name: o.name, label: o.label, type: o.type, file: o.file, brep, visible: gui ? vis.get(o.name) !== false : null, colour: col.get(o.name) || null, transform: c.M && c.M.some((x, i) => abs(x - IDENT[i]) > 1e-12) ? c.M : null });
  }
  if (!gui) notes.push('the archive has no GuiDocument.xml, so visibility is unknown: every stored shape is imported, including intermediate ones');
  if (hidden) notes.push(`${hidden} hidden object(s) (for example the inputs of Boolean operations) were left out, as in the FreeCAD 3-D view`);
  if (swall) notes.push(`${swall} feature shape(s) inside a Body or compound were left out because the enclosing object's shape already contains them`);
  if (missing) notes.push(`${missing} object(s) reference a shape file that is absent from the archive`);
  if ([...types.values()].some(ty => /^App::Link|Assembly/.test(ty))) notes.push('App::Link / assembly instances are not expanded: only objects that store their own shape are imported');
  const th = entries.find(e => e.data && /thumbnails\/.*\.png$/i.test(e.name)), img = th ? sniffImage(th.data, 0) : null; if (img) extra.preview = preview(img, `archive entry "${th.name}"`);
  return { meta, shapes, notes, extra };
}

// ---------------------------------------------------------------- Parasolid transmit data inside other files
// Neutral-binary transmit streams found by signature: stored raw ("PS\0\0\0" + modeller string) or as zlib streams
// (each candidate zlib header is inflated and kept when the result starts with the transmit signature).
// → [{offset, bytes, compressed, encoding, partition, guise, modellerVersion, schema}]
export function findParasolid(bytes, limit = 64) {
  const b = F.toBytes(bytes), out = [];
  const info = d => { const head = latin(d, 0, 200), m = /^PS\0\0\0.: TRANSMIT FILE\s*(?:\(([a-z ]+)\))?\s*created by modeller version (\d+)/.exec(head), sc = /SCH_\d+_\d+/.exec(head); return m ? { encoding: 'binary', guise: m[1] || 'part', partition: m[1] === 'partition', modellerVersion: m[2], schema: sc ? sc[0] : null } : null; };
  for (let i = 0; i + 24 < b.length && out.length < limit; i++) {
    if (b[i] === 0x50 && b[i + 1] === 0x53 && b[i + 2] === 0 && b[i + 3] === 0 && b[i + 4] === 0 && b[i + 6] === 0x3A && b[i + 7] === 0x20 && b[i + 8] === 0x54) { const d = b.subarray(i), m = info(d); if (m) out.push({ offset: i, bytes: d, compressed: false, ...m }); }
    else if (b[i] === 0x78 && (b[i + 1] === 0x01 || b[i + 1] === 0x5E || b[i + 1] === 0x9C || b[i + 1] === 0xDA)) { let r = null, plausible = true; try { F.inflateRaw(b.subarray(i + 2, i + 130)); } catch (e) { plausible = /unexpected end|exceeds input/.test(String(e && e.message)); } /* a short window rejects most false headers cheaply */ if (plausible) { try { r = F.inflateRaw(b, i + 2); } catch (e) { r = null; } } const m = r && r.data.length > 64 ? info(r.data) : null; if (m) { out.push({ offset: i, bytes: r.data, compressed: true, ...m }); i = r.end - 1; } }
  }
  return out;
}
// Siemens NX part files: an OLE2 compound file (through NX 12) or the "SPLMSSTR" container (NX 1847 and later). The solid
// model is a Parasolid partition, stored raw in the former and zlib-deflated in the latter; it is located by signature.
export function parseNX(bytes) {
  const b = F.toBytes(bytes), recovered = [], notRecovered = [], extra = {}; let meta = {}, areas = [b];
  if (startsWith(b, [0xD0, 0xCF, 0x11, 0xE0])) { const d = cfbDocument(b, 'NX (OLE2 layout)'); meta = d.meta; recovered.push(...d.recovered); Object.assign(extra, d.extra); areas = d.cfb.entries.filter(e => e.type === 'stream' && e.size > 512).map(e => d.cfb.read(e)).filter(Boolean); meta.layout = 'OLE2 compound file (NX 12 and earlier)'; }
  else if (latin(b, 0, 8) === 'SPLMSSTR') { meta.layout = 'SPLM structured storage (NX 1847 and later)'; const names = [...new Set((latin(b, 0, min(b.length, 16384)).match(/\/Root\/[\x20-\x7e]{1,60}/g) || []))]; if (names.length) { meta.entries = names.slice(0, 60); recovered.push(`container directory: ${names.length} entries (${names.slice(0, 8).join(', ')}${names.length > 8 ? ', …' : ''})`); } }
  else throw new Error('NX: neither an OLE2 compound file nor the SPLMSSTR container');
  const found = []; for (const a of areas) found.push(...findParasolid(a));
  const whole = found.filter(f => f.guise !== 'deltas'), parts = whole.some(f => f.partition) ? whole.filter(f => f.partition) : whole, deltas = found.length - whole.length, aside = whole.length - parts.length;
  if (parts.length) { extra.embedded = parts.map((f, i) => ({ name: `partition-${i + 1}.x_b`, mime: 'application/octet-stream', format: 'parasolid', bytes: f.bytes, description: `Parasolid binary transmit (${f.guise}) ${i + 1} of ${parts.length}, modeller version ${f.modellerVersion}, schema ${f.schema || 'unknown'}${f.compressed ? ', stored deflated' : ''}`, encoding: f.encoding, partition: f.partition, modellerVersion: f.modellerVersion, schema: f.schema })); meta.parasolid = extra.embedded.map(e => e.description); recovered.push(`${parts.length} embedded Parasolid ${parts[0].partition ? 'partition' : 'transmit'} stream(s) located (${extra.embedded.map(e => `${e.bytes.length} bytes, ${e.schema || 'schema not read'}`).join('; ')})${aside ? `; ${aside} further single-body transmit stream(s) (sketch, section and reference bodies) left aside` : ''}${deltas ? `; ${deltas} roll-back delta stream(s) left aside` : ''}`); }
  else notRecovered.push('no Parasolid transmit stream was found in the file');
  notRecovered.push('the NX object model (features, expressions, sketches, PMI, drawings) is stored in a proprietary serialisation and was not interpreted');
  return { meta, recovered, notRecovered, extra, geometry: null };
}

// ---------------------------------------------------------------- XZ / LZMA2 (JT 10 segment compression)
// unxz(bytes, start) → Uint8Array: the blocks of one .xz stream, each a chain with the single LZMA2 filter.
export function unxz(bytes, start = 0, limit = 1 << 30) {
  const b = F.toBytes(bytes); let o = start; const fail = m => { throw new Error('xz: ' + m); };
  if (!startsWith(b, [0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00], o)) fail('missing stream signature'); const check = b[o + 7] & 15, checkSize = check === 0 ? 0 : check < 4 ? 4 : check < 7 ? 8 : check < 10 ? 16 : check < 13 ? 32 : 64; o += 12;
  const vli = () => { let x = 0, sh = 0; for (let i = 0; i < 9; i++) { const c = b[o++]; x += (c & 127) * 2 ** sh; sh += 7; if (!(c & 128)) return x; } fail('bad integer'); };
  let out = new Uint8Array(65536), n = 0; const need = k => { if (n + k > limit) fail('output larger than allowed'); if (n + k > out.length) { const g = new Uint8Array(max(out.length * 2, n + k)); g.set(out.subarray(0, n)); out = g; } };
  // LZMA state (probabilities are reset on request of each chunk)
  let lc = 0, lp = 0, pb = 0, probs = null, state = 0, r0 = 0, r1 = 0, r2 = 0, r3 = 0, range = 0, code = 0, dictStart = 0;
  const IS_MATCH = 0, IS_REP = 192, IS_G0 = 204, IS_G1 = 216, IS_G2 = 228, IS_R0L = 240, SLOT = 432, SPEC = 688, ALIGN = 802, LEN = 818, RLEN = 1332, LIT = 1846;
  const reset = () => { probs = new Uint16Array(LIT + (0x300 << (lc + lp))).fill(1024); state = 0; r0 = r1 = r2 = r3 = 0; };
  const bit = i => { const pr = probs[i], bound = (range >>> 11) * pr; let r; if ((code >>> 0) < bound) { range = bound; probs[i] = pr + ((2048 - pr) >> 5); r = 0; } else { range = (range - bound) >>> 0; code = (code - bound) >>> 0; probs[i] = pr - (pr >> 5); r = 1; } if (range < 16777216) { range = (range << 8) >>> 0; code = ((code << 8) | b[o++]) >>> 0; } return r; };
  const tree = (base, bits) => { let m = 1; for (let i = 0; i < bits; i++) m = (m << 1) | bit(base + m); return m - (1 << bits); };
  const rtree = (base, bits) => { let m = 1, sym = 0; for (let i = 0; i < bits; i++) { const q = bit(base + m); m = (m << 1) | q; sym |= q << i; } return sym; };
  const direct = bits => { let r = 0; for (let i = 0; i < bits; i++) { range >>>= 1; let q = 0; if ((code >>> 0) >= range) { code = (code - range) >>> 0; q = 1; } r = (r * 2) + q; if (range < 16777216) { range = (range << 8) >>> 0; code = ((code << 8) | b[o++]) >>> 0; } } return r; };
  const length = (base, ps) => (bit(base) === 0 ? tree(base + 2 + ps * 8, 3) : bit(base + 1) === 0 ? 8 + tree(base + 130 + ps * 8, 3) : 16 + tree(base + 258, 8)) + 2;
  const lzma = (unpacked, end) => {
    range = 0xFFFFFFFF; code = 0; o++; for (let i = 0; i < 4; i++) code = ((code << 8) | b[o++]) >>> 0; const stop = n + unpacked; need(unpacked);
    while (n < stop) {
      if (o > end + 8) fail('compressed chunk overrun'); const pos = n - dictStart, ps = pos & ((1 << pb) - 1);
      if (bit(IS_MATCH + (state << 4) + ps) === 0) {
        const prev = pos > 0 ? out[n - 1] : 0, base = LIT + 0x300 * (((pos & ((1 << lp) - 1)) << lc) + (prev >> (8 - lc))); let sym = 1;
        if (state >= 7) { let mb = out[n - r0 - 1]; do { const q = (mb >> 7) & 1; mb <<= 1; const x = bit(base + ((1 + q) << 8) + sym); sym = (sym << 1) | x; if (q !== x) break; } while (sym < 256); }
        while (sym < 256) sym = (sym << 1) | bit(base + sym);
        out[n++] = sym & 255; state = state < 4 ? 0 : state < 10 ? state - 3 : state - 6; continue;
      }
      let len;
      if (bit(IS_REP + state) === 0) {
        r3 = r2; r2 = r1; r1 = r0; len = length(LEN, ps); state = state < 7 ? 7 : 10; const slot = tree(SLOT + (min(len - 2, 3) << 6), 6);
        if (slot < 4) r0 = slot; else { const nb = (slot >> 1) - 1; r0 = (2 | (slot & 1)) * 2 ** nb; if (slot < 14) r0 += rtree(SPEC + r0 - slot - 1, nb); else { r0 += direct(nb - 4) * 16; r0 += rtree(ALIGN, 4); } if (r0 >= 0xFFFFFFFF) fail('unexpected end marker'); }
      } else {
        if (bit(IS_G0 + state) === 0) { if (bit(IS_R0L + (state << 4) + ps) === 0) { if (n - dictStart < 1) fail('corrupt data'); state = state < 7 ? 9 : 11; out[n] = out[n - r0 - 1]; n++; continue; } }
        else { let dist; if (bit(IS_G1 + state) === 0) dist = r1; else { if (bit(IS_G2 + state) === 0) dist = r2; else { dist = r3; r3 = r2; } r2 = r1; } r1 = r0; r0 = dist; }
        len = length(RLEN, ps); state = state < 7 ? 8 : 11;
      }
      if (r0 >= n - dictStart) fail('distance beyond the data decoded so far'); if (n + len > stop) fail('match runs past the chunk');
      for (let i = 0, src = n - r0 - 1; i < len; i++) out[n++] = out[src++];
    }
  };
  for (let blk = 0; blk < 4096 && o < b.length && b[o] !== 0; blk++) {
    const h0 = o, hs = (b[o] + 1) * 4, flags = b[o + 1]; o += 2; if (flags & 3) fail('filter chains other than a single LZMA2 filter are not supported'); if (flags & 0x40) vli(); if (flags & 0x80) vli(); if (vli() !== 0x21) fail('filter is not LZMA2'); o = h0 + hs;
    let propsSet = false; dictStart = n;
    for (let g = 0; g < 1e7; g++) {
      const c = b[o++]; if (c === undefined) fail('unexpected end of data'); if (c === 0) break;
      if (c === 1 || c === 2) { const len = ((b[o] << 8) | b[o + 1]) + 1; o += 2; if (c === 1) dictStart = n; need(len); out.set(b.subarray(o, o + len), n); n += len; o += len; continue; }
      if (c < 0x80) fail('bad LZMA2 control byte');
      const unpacked = (((c & 31) << 16) | (b[o] << 8) | b[o + 1]) + 1, packed = ((b[o + 2] << 8) | b[o + 3]) + 1, mode = (c >> 5) & 3; o += 4;
      if (mode >= 2) { let d = b[o++]; if (d > 224) fail('bad LZMA properties'); lc = d % 9; d = (d - lc) / 9; lp = d % 5; pb = (d - lp) / 5; propsSet = true; } if (mode === 3) dictStart = n; if (mode >= 1) { if (!propsSet) fail('LZMA chunk before its properties'); reset(); } else if (!probs) fail('LZMA chunk without a state');
      const end = o + packed; lzma(unpacked, end); o = end;
    }
    o = h0 + Math.ceil((o - h0) / 4) * 4 + checkSize;
  }
  return out.subarray(0, n);
}

// ---------------------------------------------------------------- JT (ISO 14306)
const JT_SEG = { 1: 'logical scene graph', 2: 'JT B-rep', 3: 'PMI data', 4: 'meta data', 6: 'shape', 7: 'shape LOD 0', 8: 'shape LOD 1', 9: 'shape LOD 2', 10: 'shape LOD 3', 11: 'shape LOD 4', 12: 'shape LOD 5', 13: 'shape LOD 6', 14: 'shape LOD 7', 15: 'shape LOD 8', 16: 'shape LOD 9', 17: 'XT B-rep (Parasolid)', 18: 'wireframe representation', 20: 'ULP', 24: 'LWPA', 30: 'MultiXT B-rep', 31: 'InfoSegment', 33: 'STEP B-rep' };
export function parseJT(bytes) {
  const b = F.toBytes(bytes), v = dv(b), head = latin(b, 0, 80), m = /^Version (\d+)\.(\d+)(?:\s*JT|\s{8})/.exec(head); if (!m) throw new Error('JT: missing "Version n.n JT" header');
  const major = +m[1], meta = { version: `${m[1]}.${m[2]}`, header: head.trim().replace(/\s+/g, ' ') }, recovered = [`JT version ${meta.version}`], notRecovered = [], extra = {}, le = b[80] === 0; meta.byteOrder = le ? 'little-endian' : 'big-endian';
  try {
    const wide = major >= 10, toc = wide ? Number(v.getBigUint64(85, le)) : v.getInt32(85, le); if (!(toc > 0 && toc + 4 <= b.length)) throw new Error('table of contents offset outside the file');
    const n = v.getInt32(toc, le), es = wide ? 16 + 8 + 4 + 4 : 16 + 4 + 4 + 4; if (!(n >= 0 && n < 1e6 && toc + 4 + n * es <= b.length)) throw new Error('table of contents does not fit in the file');
    const counts = {}, sizes = {}, xt = []; let lods = 0, tri = 0, xtSegs = 0, lzma = 0;
    for (let i = 0; i < n; i++) {
      const o = toc + 4 + i * es, so = wide ? Number(v.getBigUint64(o + 16, le)) : v.getInt32(o + 16, le), len = v.getInt32(o + (wide ? 24 : 20), le), attr = v.getUint32(o + (wide ? 28 : 24), le), type = attr >>> 24, name = JT_SEG[type] || `segment type ${type}`; counts[name] = (counts[name] || 0) + 1; sizes[name] = (sizes[name] || 0) + max(0, len); if (type >= 6 && type <= 16) { lods++; tri += max(0, len); }
      // XT B-rep segments carry a Parasolid transmit stream, stored raw or zlib-deflated behind the segment header: found by signature
      if ((type === 17 || type === 30) && so > 0 && len > 24 && so + len <= b.length) {
        xtSegs++; let data = b.subarray(so, so + len), packed = false;
        // segment header (id, type, length), then compression flag 2, compressed length and algorithm 2 (zlib) or 3 (LZMA)
        if (len > 40 && (v.getInt32(so + 24, le) === 2 || v.getInt32(so + 24, le) === 3)) { const algo = b[so + 32]; if (algo === 2 && b[so + 33] === 0x78) { try { data = F.inflateRaw(b, so + 35).data; packed = true; } catch (e) { data = b.subarray(so, so + len); } } else if (algo === 3 && b[so + 33] === 0xFD && b[so + 34] === 0x37) { try { data = unxz(b, so + 33, 4e8); packed = true; } catch (e) { lzma++; continue; } } }
        for (const f of findParasolid(data)) if (f.guise !== 'deltas') xt.push({ ...f, compressed: packed || f.compressed, segment: xtSegs, kind: name });
      }
    }
    if (xt.length) { extra.embedded = xt.map((f, i) => ({ name: `xt-brep-${i + 1}.x_b`, mime: 'application/octet-stream', format: 'parasolid', bytes: f.bytes, description: `Parasolid binary transmit (${f.guise}) ${i + 1} of ${xt.length} from JT ${f.kind} segment ${f.segment}, modeller version ${f.modellerVersion}, schema ${f.schema || 'unknown'}${f.compressed ? ', stored deflated' : ''}`, encoding: f.encoding, partition: f.partition, modellerVersion: f.modellerVersion, schema: f.schema })); meta.parasolid = extra.embedded.map(e => e.description); recovered.push(`${xt.length} Parasolid transmit stream(s) unpacked from ${xtSegs} XT B-rep segment(s) — the exact B-rep of the parts, each in its own part coordinates`); }
    else if (xtSegs) notRecovered.push(`${xtSegs} XT B-rep segment(s) are present but no Parasolid transmit stream could be unpacked from them`);
    if (lzma) notRecovered.push(`${lzma} LZMA-compressed XT B-rep segment(s) could not be inflated: their Parasolid data was not unpacked`);
    meta.segments = n; meta.segmentCounts = counts; meta.segmentBytes = sizes;
    recovered.push(`table of contents: ${n} data segments — ${Object.entries(counts).map(([k, c]) => `${c} × ${k} (${sizes[k]} bytes)`).join(', ')}`);
    if (lods) notRecovered.push(`${lods} tessellated shape / LOD segment(s) (${tri} bytes) are present but were not decoded: JT shape data use quantised, entropy-coded vertex and topology streams whose decoder is not implemented here`);
    if (counts[JT_SEG[2]] || counts[JT_SEG[33]]) notRecovered.push('the JT B-rep / STEP B-rep segments were inventoried but not decoded');
    if (xt.length) notRecovered.push('the logical scene graph (assembly structure, instance placements, names, materials) was not decoded: bodies are not placed by their instance transforms, and a part instanced several times appears once');
  } catch (e) { notRecovered.push(`the table of contents could not be read (${e.message}): only the version string is reported`); }
  return { meta, recovered, notRecovered, extra, geometry: null };
}

// ---------------------------------------------------------------- Parasolid transmit files (.x_t / .x_b)
export function parseParasolid(bytes) {
  const b = F.toBytes(bytes), head = latin(b, 0, 4096), meta = { ...parasolidInfo(b) }, recovered = [], notRecovered = [];
  for (const k of ['MC', 'MC_MODEL', 'MC_ID', 'OS', 'OS_RELEASE', 'FRU', 'APPL', 'SITE', 'USER', 'FORMAT', 'GUISE', 'KEY', 'FILE', 'DATE', 'SCH', 'USFLD_SIZE']) { const r = new RegExp('(?:^|[\\r\\n;])' + k + '=([^;\\r\\n]*)').exec(head); if (r && r[1].trim()) meta[k.toLowerCase()] = r[1].trim(); }
  if (meta.sch && !meta.schema) meta.schema = meta.sch;
  recovered.push(`Parasolid ${meta.encoding} transmit file${meta.partition ? ' (partition)' : ''}${meta.modellerVersion ? `, written by modeller version ${meta.modellerVersion}` : ''}${meta.schema ? `, schema ${meta.schema}` : ''}`);
  const hk = ['appl', 'date', 'key', 'file', 'format', 'guise'].filter(k => meta[k]).map(k => `${k}: ${meta[k]}`); if (hk.length) recovered.push(`header — ${hk.join('; ')}`);
  return { meta, recovered, notRecovered, extra: {}, geometry: null };
}

// ---------------------------------------------------------------- ACIS SAT (text)
// Header, entity inventory and the topological vertex points; edges are drawn as straight chords between their end
// vertices. No surface or curve is evaluated.
export function parseSAT(bytes) {
  const t = F.toText(F.toBytes(bytes)), L = t.split(/\r?\n/), h = (L[0] || '').trim().split(/\s+/).map(Number); if (!(h.length >= 3 && h[0] >= 100 && h.slice(0, 3).every(Number.isFinite))) throw new Error('SAT: the first line is not an ACIS header (version, record count, body count)');
  const meta = { version: h[0] / 100, records: h[1], bodies: h[2] }, recovered = [], notRecovered = []; let line = 1, scale = 1;
  if (h[0] >= 200) { const p = /^\s*(\d+) (.*)$/.exec(L[1] || ''); if (p) { let rest = L[1].trim(); const take = () => { const q = /^(\d+) /.exec(rest); if (!q) return ''; const n = +q[1], s = rest.slice(q[0].length, q[0].length + n); rest = rest.slice(q[0].length + n).trimStart(); return s; }; meta.product = take(); meta.acis = take(); meta.date = take(); } const u = (L[2] || '').trim().split(/\s+/).map(Number); if (u.length >= 2 && u[0] > 0) { meta.unitScale_mm = u[0]; meta.resabs = u[1]; meta.resnor = u[2]; scale = u[0]; } line = 3; }
  const body = L.slice(line).join('\n'), recs = body.split('#'), ent = [], counts = {};
  for (const r of recs) { const s = r.trim(); if (!s) continue; if (/^End-of-ACIS-data|^End of ACIS/i.test(s)) break; const m = /^(?:-(\d+)\s+)?([A-Za-z_][\w-]*)\s*([\s\S]*)$/.exec(s); if (!m) continue; const idx = m[1] !== undefined ? +m[1] : ent.length; ent[idx] = { type: m[2], tok: m[3].split(/\s+/).filter(Boolean) }; const base = m[2].split('-')[0] === m[2] ? m[2] : m[2]; counts[base] = (counts[base] || 0) + 1; }
  meta.entityCounts = Object.fromEntries(Object.entries(counts).sort((a, c) => c[1] - a[1]).slice(0, 40));
  recovered.push(`ACIS ${meta.version}${meta.product ? ` written by ${meta.product}` : ''}${meta.acis ? ` (${meta.acis})` : ''}${meta.date ? `, ${meta.date}` : ''}`);
  recovered.push(`entity inventory — ${Object.entries(meta.entityCounts).slice(0, 14).map(([k, c]) => `${k}: ${c}`).join(', ')}`);
  const isNum = s => /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(s), ptr = s => /^\$-?\d+$/.test(s) ? +s.slice(1) : null, nodes = [], pid = new Map();
  ent.forEach((e, i) => { if (!e || e.type !== 'point') return; const nums = e.tok.filter(isNum), n = nums.length; if (n < 3) return; pid.set(i, nodes.length); nodes.push([+nums[n - 3] * scale, +nums[n - 2] * scale, +nums[n - 1] * scale]); });
  const vpt = new Map(); ent.forEach((e, i) => { if (!e || e.type !== 'vertex') return; for (const s of e.tok) { const p = ptr(s); if (p !== null && pid.has(p)) { vpt.set(i, pid.get(p)); break; } } });
  const edges = [], seen = new Set(); let curved = 0;
  ent.forEach(e => { if (!e || e.type !== 'edge') return; const vs = []; let crv = null; for (const s of e.tok) { const p = ptr(s); if (p === null || p < 0) continue; if (vpt.has(p)) vs.push(vpt.get(p)); else if (ent[p] && /-curve$/.test(ent[p].type)) crv = ent[p].type; } if (crv && crv !== 'straight-curve') curved++; if (vs.length >= 2 && vs[0] !== vs[1]) { const k = min(vs[0], vs[1]) + ',' + max(vs[0], vs[1]); if (!seen.has(k)) { seen.add(k); edges.push([vs[0], vs[1]]); } } });
  const moved = ent.filter(e => e && e.type === 'transform').some(e => { const n = e.tok.filter(isNum).map(Number); return n.length >= 13 && (abs(n[0] - 1) > 1e-12 || abs(n[4] - 1) > 1e-12 || abs(n[8] - 1) > 1e-12 || abs(n[9]) + abs(n[10]) + abs(n[11]) > 1e-12 || abs(n[12] - 1) > 1e-12); });
  notRecovered.push('no face was tessellated: the ACIS surfaces and curves (planes, cones, spheres, tori, splines) are counted but not evaluated');
  if (!nodes.length) return { meta, recovered, notRecovered, extra: {}, geometry: null };
  const notes = [`${nodes.length} topological vertex points${edges.length ? ` and ${edges.length} edges drawn as straight chords between their end vertices` : ''}${meta.unitScale_mm ? `, scaled to millimetres by the header factor ${meta.unitScale_mm}` : ''}`];
  if (curved) notes.push(`${curved} edges lie on curved geometry (ellipse / spline curves): their chords cut across the true curve, and closed circular edges with a single vertex are not drawn`);
  if (moved) notes.push('bodies carry placement transforms that were not applied: points are in each body\'s local coordinates, so the extents are not those of the placed model');
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (const p of nodes) for (let a = 0; a < 3; a++) { lo[a] = min(lo[a], p[a]); hi[a] = max(hi[a], p[a]); }
  meta.vertexBoundingBox = { min: lo, max: hi, note: 'bounding box of the topological vertices only: curved faces can extend beyond it' };
  recovered.push(`vertex bounding box ${lo.map(x => +x.toPrecision(6)).join(', ')} → ${hi.map(x => +x.toPrecision(6)).join(', ')} ${meta.unitScale_mm ? 'mm' : 'file units'} (vertices only: curved faces can extend beyond it)`);
  return { meta, recovered, notRecovered, extra: {}, geometry: { kind: 'lattice', nodes, edges, radius: 0, units: meta.unitScale_mm ? 'mm' : '', wireframe: true, pointCloud: !edges.length, notes, lost: [] } };
}

// ---------------------------------------------------------------- CATIA V5 and PTC Creo: container facts and preview only
// Complete JPEG images embedded anywhere in a file → [{offset, bytes, width, height}] (marker segments walked to the frame header and to the end marker)
export function findJPEG(bytes, limit = 16) {
  const b = F.toBytes(bytes), out = [];
  for (let i = 0; i + 12 < b.length && out.length < limit; i++) {
    if (b[i] !== 0xFF || b[i + 1] !== 0xD8 || b[i + 2] !== 0xFF || b[i + 3] < 0xC0) continue;
    let p = i + 2, w = 0, h = 0, okay = false;
    for (let g = 0; g < 256 && p + 4 <= b.length && b[p] === 0xFF; g++) { const mk = b[p + 1], len = (b[p + 2] << 8) | b[p + 3]; if (mk >= 0xC0 && mk <= 0xC2 && p + 9 <= b.length) { h = (b[p + 5] << 8) | b[p + 6]; w = (b[p + 7] << 8) | b[p + 8]; } if (mk === 0xDA) { okay = w > 0 && h > 0; break; } if (len < 2) break; p += 2 + len; }
    if (!okay) continue; let e = -1; for (let q = p; q + 1 < b.length && q < p + 2e7; q++) if (b[q] === 0xFF && b[q + 1] === 0xD9) { e = q + 2; break; } if (e < 0) continue;
    out.push({ offset: i, bytes: b.subarray(i, e), width: w, height: h }); i = e - 1;
  }
  return out;
}
// ---- CATIA V5: container directory and the display tessellation of the CGR cache
// Directory of a "V5_CFV2" container → [{name, size, ext: [{off, len}]}] or null. The header gives the big-endian offset
// of the directory at the end of the file; every stream has one record there (length 0x6A + 20 per further extent, its
// UTF-16BE name in front of it, each extent as offset / stored / logical length / position). Streams are not compressed
// but can be split over several extents.
export function catiaDirectory(bytes) {
  const b = F.toBytes(bytes); if (b.length < 64 || latin(b, 0, 7) !== 'V5_CFV2') return null; const v = dv(b), dirOff = v.getUint32(8), out = []; if (!(dirOff > 0 && dirOff + 16 <= b.length)) return null;
  for (let o = dirOff; o + 0x6A <= b.length; o++) {
    const L = v.getUint32(o); if (L < 0x6A || L > 0x6A + 20 * 65536 || (L - 0x6A) % 20 || o + L > b.length) continue; const t = v.getUint16(o + 0x30); if (t !== 0x12 && t !== 0x11) continue;
    const n = v.getUint32(o + 0x52); if (n !== (L - 0x6A) / 20 + 1) continue; const size = v.getUint32(o + 14), ext = []; let okay = true, sum = 0;
    for (let i = 0; i < n; i++) { const q = o + 0x56 + 20 * i, off = v.getUint32(q), len = v.getUint32(q + 4); if (off + len > b.length || v.getUint32(q + 8) !== len || v.getUint32(q + 12) !== sum) { okay = false; break; } sum += len; ext.push({ off, len }); }
    if (!okay || sum !== size) continue; const e = o - 2; let st = e; while (st - 2 >= dirOff && b[st - 2] === 0 && b[st - 1] >= 0x20 && b[st - 1] < 0x7F) st -= 2;
    let name = ''; for (let i = st; i < e; i += 2) name += String.fromCharCode(b[i + 1]); out.push({ name, size, ext }); o += L - 1;
  }
  return out;
}
export function catiaStream(bytes, entry) { const b = F.toBytes(bytes); if (entry.ext.length === 1) return b.subarray(entry.ext[0].off, entry.ext[0].off + entry.ext[0].len); const d = new Uint8Array(entry.size); let q = 0; for (const x of entry.ext) { d.set(b.subarray(x.off, x.off + x.len), q); q += x.len; } return d; }
// One tessellated body block of the cache's "MainDataStream", starting at its marker FF FF 02 00 01. Counts are one byte,
// or FF and a little-endian u32. Layout: A (position entries), B (entries with a normal), u32 0; a 2-bit code per entry
// (0: x y z follow, 3: y z follow, 2: z follows, 1: the same point again); a count and that many float32 coordinates;
// a count and that many int16 normal components; B normal nibbles; float32 chordal tolerance, 00 04 02 0A, face count;
// per face 01, flags, the numbers of single triangles / strips / fans (flag bits 1, 2, 4), FF and the u32 index total,
// a plane normal of three float32 (flag bit 8), the strip and fan lengths, then big-endian u16 indices into the entries.
// Every count has to agree with the others, or the block is refused (null).
function catiaBlock(m, o) {
  const v = dv(m), end = m.length; let p = o + 5, bad = false;
  const cnt = () => { if (p >= end) { bad = true; return 0; } if (m[p] !== 0xFF) return m[p++]; if (p + 5 > end) { bad = true; return 0; } const x = v.getUint32(p + 1, true); p += 5; return x; };
  const A = cnt(), B = cnt(); if (bad || A < 3 || B > A || A > 65535 || p + 4 > end || v.getUint32(p, true) !== 0) return null; p += 4;
  const sym = p; p += Math.ceil(A / 4); if (p >= end) return null; const code = i => (m[sym + (i >> 2)] >> ((i & 3) * 2)) & 3, n = [0, 0, 0, 0]; for (let i = 0; i < A; i++) n[code(i)]++;
  const nf = cnt(); if (bad || nf !== 3 * n[0] + 2 * n[3] + n[2] || p + 4 * nf > end) return null;
  const P = new Float32Array(3 * A); let x = 0, y = 0, z = 0;
  for (let i = 0; i < A; i++) { const c = code(i); if (c === 0) { x = v.getFloat32(p, true); y = v.getFloat32(p + 4, true); z = v.getFloat32(p + 8, true); p += 12; } else if (c === 3) { y = v.getFloat32(p, true); z = v.getFloat32(p + 4, true); p += 8; } else if (c === 2) { z = v.getFloat32(p, true); p += 4; } if (!(Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z))) return null; P[3 * i] = x; P[3 * i + 1] = y; P[3 * i + 2] = z; }
  const M = cnt(); if (bad || M & 1 || p + 2 * M + Math.ceil(B / 2) + 9 > end) return null; p += 2 * M;
  let stored = 0; for (let i = 0; i < B; i++) { const c = m[p + (i >> 1)], k = i & 1 ? c >> 4 : c & 15; if (k < 6) stored++; else if (k > 11) return null; } if (stored * 2 !== M) return null; p += Math.ceil(B / 2);
  const sag = v.getFloat32(p, true); if (m[p + 4] !== 0 || m[p + 5] !== 4 || m[p + 6] !== 2 || m[p + 7] !== 0x0A) return null; p += 8;
  const faces = cnt(); if (bad) return null; const tris = [];
  for (let f = 0; f < faces; f++) {
    if (p + 7 > end || m[p] !== 1) return null; const flags = m[p + 1]; p += 2; const nt = flags & 1 ? cnt() : 0, ns = flags & 2 ? cnt() : 0, nfan = flags & 4 ? cnt() : 0; if (bad || m[p] !== 0xFF || p + 5 > end) return null;
    const total = v.getUint32(p + 1, true); p += 5; if (flags & 8) p += 12; const sl = [], fl = []; let sum = 3 * nt; for (let i = 0; i < ns; i++) { const c = cnt(); sl.push(c); sum += c; } for (let i = 0; i < nfan; i++) { const c = cnt(); fl.push(c); sum += c; }
    if (bad || sum !== total || p + 2 * total > end) return null; const ix = k => v.getUint16(p + 2 * k, false); for (let k = 0; k < total; k++) if (ix(k) >= A) return null;
    let q = 0; for (let i = 0; i < nt; i++, q += 3) tris.push(ix(q), ix(q + 1), ix(q + 2));
    for (const len of sl) { for (let i = 0; i + 2 < len; i++) { if (i & 1) tris.push(ix(q + i + 1), ix(q + i), ix(q + i + 2)); else tris.push(ix(q + i), ix(q + i + 1), ix(q + i + 2)); } q += len; }
    for (const len of fl) { for (let i = 1; i + 1 < len; i++) tris.push(ix(q), ix(q + i), ix(q + i + 1)); q += len; }
    p += 2 * total;
  }
  return { A, P, tris, faces, sag, end: p };
}
// The display tessellation CATIA saved with a part → {vertices, faces, units: 'mm', blocks: [{triangles, faces, sag,
// faceRange, area, volume, closed}], skipped, dropped} or null (no cache, or no block that passes its count checks).
// Entries with identical coordinates are welded; the stored winding is kept. This is the visualisation mesh (float32,
// millimetres, at the chordal tolerance it was saved with), not the exact B-rep.
export function catiaDisplayMesh(bytes) {
  const b = F.toBytes(bytes), dir = catiaDirectory(b); if (!dir) return null; const mains = [];
  for (const e of dir) { if (e.size < 128) continue; if (e.name === 'MainDataStream') { mains.push(catiaStream(b, e)); continue; } const st = catiaStream(b, e), d2 = catiaDirectory(st); if (d2) for (const x of d2) if (x.name === 'MainDataStream' && x.size > 128) mains.push(catiaStream(st, x)); }
  const vertices = [], faces = [], blocks = []; let skipped = 0, dropped = 0;
  for (const m of mains) for (let o = 0; o + 32 < m.length; o++) {
    if (m[o] !== 0xFF || m[o + 1] !== 0xFF || m[o + 2] !== 2 || m[o + 3] !== 0 || m[o + 4] !== 1) continue; const r = catiaBlock(m, o); if (!r || !r.tris.length) { if (!r && m[o + 5] === 0xFF) skipped++; continue; } o = r.end - 1;
    const seen = new Map(), map = new Int32Array(r.A), bits = new Int32Array(r.P.buffer), f0 = faces.length;
    for (let i = 0; i < r.A; i++) { const key = bits[3 * i] + ',' + bits[3 * i + 1] + ',' + bits[3 * i + 2]; let j = seen.get(key); if (j === undefined) { j = vertices.length; vertices.push([r.P[3 * i], r.P[3 * i + 1], r.P[3 * i + 2]]); seen.set(key, j); } map[i] = j; }
    const edges = new Map(); let area = 0, vol = 0;
    for (let t = 0; t < r.tris.length; t += 3) { const f = [map[r.tris[t]], map[r.tris[t + 1]], map[r.tris[t + 2]]]; if (f[0] === f[1] || f[1] === f[2] || f[0] === f[2]) { dropped++; continue; } faces.push(f);
      const p0 = vertices[f[0]], p1 = vertices[f[1]], p2 = vertices[f[2]], ux = p1[0] - p0[0], uy = p1[1] - p0[1], uz = p1[2] - p0[2], wx = p2[0] - p0[0], wy = p2[1] - p0[1], wz = p2[2] - p0[2]; area += hypot(uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx) / 2; vol += (p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])) / 6;
      for (let k = 0; k < 3; k++) { const i = f[k], j = f[(k + 1) % 3], key = i < j ? i * 4294967296 + j : j * 4294967296 + i, e = edges.get(key) || [0, 0]; e[0]++; e[1] += i < j ? 1 : -1; edges.set(key, e); } }
    let closed = edges.size > 0; for (const e of edges.values()) if (e[0] !== 2 || e[1] !== 0) { closed = false; break; }
    blocks.push({ triangles: faces.length - f0, faces: r.faces, sag: r.sag, faceRange: [f0, faces.length], area, volume: closed ? vol : undefined, closed });
  }
  return blocks.length ? { vertices, faces, units: 'mm', blocks, skipped, dropped } : null;
}
const WHY_CATIA = 'no geometry: this file holds no display tessellation that passes the reader\'s checks (a CATProduct, or a part saved without its CGR display cache), and the exact part body in the CGMGeom container is in the binary encoding of Dassault\'s CGM kernel, which is not decoded';
export function parseCatia(bytes) {
  const b = F.toBytes(bytes), recovered = [], notRecovered = [], extra = {}, meta = { container: 'CATIA V5 "V5_CFV2" container' }; if (latin(b, 0, 7) !== 'V5_CFV2') throw new Error('CATIA: missing V5_CFV2 signature');
  const names = new Map(); { const CH = 1 << 20; for (let o = 0; o < b.length; o += CH) { const t = latin(b, o, CH + 64); for (const m of t.matchAll(/CAT[A-Za-z0-9]{3,30}Cont(?:ainer)?(?![A-Za-z])|CGMGeom|CATMFBRP|CATGEO3D/g)) names.set(m[0], (names.get(m[0]) || 0) + 1); } }
  let dir = null; try { dir = catiaDirectory(b); } catch (e) { dir = null; }
  if (dir && dir.length) { meta.streams = dir.length; meta.streamList = dir.slice().sort((x, y) => y.size - x.size).slice(0, 60).map(e => `${e.name || '(unnamed)'} (${e.size} bytes${e.ext.length > 1 ? `, ${e.ext.length} extents` : ''})`); recovered.push(`container directory: ${dir.length} streams, stored uncompressed (${dir.filter(e => e.ext.length > 1).length} of them split over several extents)`); }
  if (names.size) { meta.containers = [...names.keys()].sort(); recovered.push(`${names.size} internal container names: ${meta.containers.slice(0, 16).join(', ')}${names.size > 16 ? ', …' : ''}`); meta.hasGeometryContainer = names.has('CGMGeom'); meta.hasDisplayCache = names.has('CATCGRCont') || names.has('CATSeeBodyCont'); }
  const imgs = findJPEG(b).sort((x, y) => y.bytes.length - x.bytes.length); if (imgs.length) { extra.preview = { mime: 'image/jpeg', bytes: imgs[0].bytes, source: `an embedded JPEG image at byte ${imgs[0].offset}`, displayable: true, width: imgs[0].width, height: imgs[0].height }; recovered.push(describePreview(extra.preview)); } else notRecovered.push('no preview image was found');
  let mesh = null; try { mesh = catiaDisplayMesh(b); } catch (e) { mesh = null; }
  let geometry = null;
  if (mesh) {
    const nf = mesh.blocks.reduce((a, q) => a + q.faces, 0), sags = [...new Set(mesh.blocks.map(q => +q.sag.toPrecision(3)))].join(' / '), closed = mesh.blocks.filter(q => q.closed).length;
    recovered.push(`display tessellation saved with the part (CGR cache): ${mesh.faces.length} triangles over ${nf} B-rep faces in ${mesh.blocks.length} block(s), ${closed} of them closed, in millimetres`);
    meta.displayMesh = { blocks: mesh.blocks.length, brepFaces: nf, triangles: mesh.faces.length, points: mesh.vertices.length, chordalTolerance: mesh.blocks.map(q => q.sag), closedBlocks: closed };
    extra.parts = mesh.blocks.map((q, i) => ({ name: mesh.blocks.length > 1 ? `display block ${i + 1}` : 'display tessellation', faceRange: q.faceRange, triangles: q.triangles, area: q.area, closed: q.closed, volume: q.volume, type: 'CATIA display tessellation' }));
    const notes = [`geometry path: the display tessellation CATIA saved with the part (the CGR cache in the CATCGRCont container) — triangles, strips and fans of ${nf} B-rep faces read directly, ${mesh.faces.length} triangles on ${mesh.vertices.length} distinct points; this is the visualisation mesh stored in the file, not the exact B-rep`,
      `coordinates are single-precision millimetres and curved faces are facetted at the chordal tolerance the cache was saved with (${sags} mm): areas and volumes carry that facetting error, and a cache saved at a coarse setting gives a coarse mesh`,
      'the cache reflects the part as it was last displayed and saved: whether it is up to date with the feature tree is not checked'];
    if (mesh.dropped) notes.push(`${mesh.dropped} degenerate strip or fan triangles (repeated points) were dropped`);
    const lost = ['the exact part body (CGMGeom container, the encoding of Dassault\'s CGM kernel) is not decoded: there are no analytic surfaces, edges or face identities in the result, only triangles']; if (mesh.skipped) lost.push(`${mesh.skipped} further tessellation block(s) did not pass the count checks and were left out: the mesh may be incomplete`);
    geometry = { kind: 'mesh', vertices: mesh.vertices, faces: mesh.faces, units: 'mm', notes, lost, meta: {} }; if (mesh.blocks.length > 1) geometry.shells = `${mesh.blocks.length} tessellation blocks kept as separate shells`; if (mesh.skipped) geometry.complete = false;
  } else notRecovered.push(WHY_CATIA);
  notRecovered.push('the specification tree, sketches, parameters, PMI and material are stored as serialised CATIA feature objects and were not interpreted');
  return { meta, recovered, notRecovered, extra, geometry };
}
// ---- PTC Creo: header, table of contents and what the sections are
// creoSections(bytes) → {objectType, release, headerLength, sections: [{name, offset, stored, size, kind, flag, body,
// bodyLength, compressed}]} or null. The last field of the first header line is the header length, to which every
// table-of-contents offset is relative; the table is a run of fixed-width lines (name, offset, stored size, expanded
// size, release, kind, flag, two checksums; hexadecimal) that a NEXT_TOC_ENTRY line chains to a second table. A section
// begins with a "#name" line; flag -1 means the body is stored as it is, a positive flag that it is a Unix-compress
// (LZW) stream. Pro/ENGINEER files of the 1990s have no table: null.
export function creoSections(bytes) {
  const b = F.toBytes(bytes), head = latin(b, 0, 4096), h = /^#UGC:(\d+) (\w+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) ([0-9a-f]{8})/.exec(head); if (!h) return null;
  const base = parseInt(h[11], 16), sections = []; let pos = head.indexOf('#UGC_TOC');
  for (let g = 0; pos >= 0 && pos < b.length && g < 64; g++) {
    const m = /^#UGC_TOC (\d+) (\d+) (\d+) (\d+)/.exec(latin(b, pos, 81)); if (!m) break; const n = min(+m[2], 100000), w = +m[3]; let next = -1; if (w < 20 || w > 512) break;
    for (let i = 1; i <= n && pos + (i + 1) * w <= b.length; i++) {
      const f = latin(b, pos + i * w, w).replace(/#+\s*$/, '').split(' '); if (f.length < 8) continue; const offset = parseInt(f[1], 16) + base, stored = parseInt(f[2], 16), size = parseInt(f[3], 16); if (f[0] === 'NEXT_TOC_ENTRY') { next = offset; continue; }
      if (!(offset >= 0 && stored >= 0 && offset + stored <= b.length)) continue; let nl = offset; while (nl < offset + stored && b[nl] !== 10) nl++; const body = nl + 1;
      sections.push({ name: f[0].replace(/^ND:\d+:/, ''), offset, stored, size, kind: parseInt(f[5], 16), flag: +f[6], body, bodyLength: max(0, offset + stored - body), compressed: +f[6] > 0 && b[body] === 0x1F && b[body + 1] === 0x9D });
    }
    pos = next;
  }
  return sections.length ? { objectType: h[2].toLowerCase(), release: +h[3], headerLength: base, sections } : null;
}
const WHY_CREO = 'no geometry: the file holds no tessellated or neutral data (its only display data are the JPEG view thumbnails and one display table in an undecoded packing), and the surface table of the solid (section VisibGeom: planes, cylinders, cones, tori with their local systems and radii), although present, is in PTC\'s tagged object encoding — field names spelled out once, later records positional, real numbers byte-coded — whose record structure is not decoded, so no face can be rebuilt and checked';
export function parseCreo(bytes) {
  const b = F.toBytes(bytes), head = latin(b, 0, min(b.length, 65536)), recovered = [], notRecovered = [], extra = {}, meta = {}; if (!head.startsWith('#UGC:')) throw new Error('Creo: missing #UGC header');
  const h1 = /^#UGC:(\d+)\s+(\w+)\s+(\d+)\s+(\d+)/.exec(head); if (h1) { meta.headerVersion = +h1[1]; meta.objectType = h1[2].toLowerCase(); meta.releaseCode = +h1[3]; meta.minimumReleaseCode = +h1[4]; }
  for (const [k, key] of [['MACH', 'platform'], ['CMNM', 'originalName'], ['HOST', 'host']]) { const m = new RegExp('#- ' + k + ' +([^\\\\\\r\\n]*?)\\s*\\\\').exec(head); if (m && m[1].trim()) meta[key] = m[1].trim().replace(/^[0-9a-f]{3}/, ''); }
  const rel = /#((?:Creo|Pro\/ENGINEER|Wildfire)[^\r\n]{0,90})/.exec(head); if (rel) meta.release = rel[1].replace(/\s+/g, ' ').trim();
  recovered.push(`${meta.release || 'PTC model file'}${meta.objectType ? `, ${meta.objectType}` : ''}${meta.originalName ? `, saved as ${meta.originalName}` : ''}${meta.platform ? ` on ${meta.platform.replace(/^_/, '')}` : ''}`);
  let toc = null; try { toc = creoSections(b); } catch (e) { toc = null; }
  if (toc) { const S = toc.sections, packed = S.filter(x => x.compressed), plain = S.filter(x => x.flag < 0), bytesPlain = plain.reduce((a, x) => a + x.stored, 0), total = S.reduce((a, x) => a + x.stored, 0), vg = S.find(x => x.name === 'VisibGeom');
    meta.sections = S.slice(0, 120).map(t => `${t.name} (${t.stored} bytes stored, ${t.size} expanded${t.compressed ? ', LZW' : ''})`); meta.sectionCount = S.length; meta.hasSurfaceTable = !!vg;
    recovered.push(`table of contents: ${S.length} model sections (${S.slice(0, 8).map(t => t.name).join(', ')}${S.length > 8 ? ', …' : ''}); ${plain.length} are stored as they are (${total ? (100 * bytesPlain / total).toFixed(0) : 0} % of the section bytes) and ${packed.length} as Unix-compress (LZW) streams${vg ? `; the surface table VisibGeom is ${vg.stored} bytes${vg.compressed ? '' : ', uncompressed'}` : ''}`); }
  else { const n = (head.match(/(?:^|\n)#[A-Z][A-Za-z_]{3,40}\r?\n/g) || []).length; if (n) recovered.push(`${n} section heading(s) in the first 64 kB (a file without a table of contents, as written by Pro/ENGINEER releases of the 1990s)`); }
  const imgs = findJPEG(b).sort((x, y) => y.width * y.height - x.width * x.height); if (imgs.length) { extra.preview = { mime: 'image/jpeg', bytes: imgs[0].bytes, source: `an embedded JPEG image at byte ${imgs[0].offset} (${imgs.length} view thumbnail(s) in the file)`, displayable: true, width: imgs[0].width, height: imgs[0].height }; recovered.push(describePreview(extra.preview)); } else notRecovered.push('no preview image was found');
  notRecovered.push(WHY_CREO, 'features, sketches, parameters, relations and annotations were not interpreted');
  return { meta, recovered, notRecovered, extra, geometry: null };
}

// ---------------------------------------------------------------- ACIS SAB / ShapeManager binary
// Header and record inventory of a binary save file (the B-rep itself is translated by sat.js).
export function parseSAB(bytes) {
  const m = readSAT(F.toBytes(bytes)), meta = { signature: m.signature, version: m.version / 100, product: m.product, acis: m.acis, date: m.date, unitScale_mm: m.scale, resabs: m.resabs, records: m.entities.size }, recovered = [], notRecovered = [];
  meta.entityCounts = Object.fromEntries(Object.entries(m.counts).sort((a, c) => c[1] - a[1]).slice(0, 40));
  recovered.push(`${m.signature} save file, version ${meta.version}${m.product ? ` written by ${m.product}` : ''}${m.acis ? ` (${m.acis})` : ''}${m.date ? `, ${m.date}` : ''}`);
  recovered.push(`record inventory — ${Object.entries(meta.entityCounts).slice(0, 14).map(([k, c]) => `${k}: ${c}`).join(', ')}`);
  if (m.history) notRecovered.push('the roll-back history section that follows the model was not read');
  return { meta, recovered, notRecovered, extra: {}, geometry: null };
}

// ---------------------------------------------------------------- Unix compress (.Z, LZW) — used for .stp.Z
export function unLZW(bytes) {
  const b = F.toBytes(bytes); if (b.length < 3 || b[0] !== 0x1F || b[1] !== 0x9D) throw new Error('compress: missing 1F 9D signature');
  const maxBits = b[2] & 0x1F, block = !!(b[2] & 0x80); if (maxBits < 9 || maxBits > 16) throw new Error(`compress: unsupported code width ${maxBits}`);
  const prefix = new Uint16Array(1 << maxBits), suffix = new Uint8Array(1 << maxBits), stack = new Uint8Array(1 << maxBits); let out = new Uint8Array(max(4096, b.length * 3)), op = 0;
  const put = c => { if (op >= out.length) { const g = new Uint8Array(out.length * 2); g.set(out); out = g; } out[op++] = c; };
  let bits = 9, free = block ? 257 : 256, bitPos = 24, groupStart = 24, old = -1, fin = 0;
  const maxCode = () => (bits === maxBits ? 1 << maxBits : (1 << bits) - 1);
  for (;;) {
    if (bitPos + bits > b.length * 8) break; let code = 0; for (let i = 0; i < bits; i++) { const p = bitPos + i; code |= ((b[p >> 3] >> (p & 7)) & 1) << i; } bitPos += bits;
    if (old < 0) { old = fin = code; if (code > 255) throw new Error('compress: corrupt stream'); put(code); continue; }
    if (code === 256 && block) { free = 256; const span = bits * 8, used = bitPos - groupStart; bitPos = groupStart + Math.ceil(used / span) * span; groupStart = bitPos; bits = 9; continue; }
    const incode = code; let sp = 0; if (code >= free) { if (code > free) throw new Error('compress: corrupt stream'); stack[sp++] = fin; code = old; }
    while (code >= 256) { stack[sp++] = suffix[code]; code = prefix[code]; } fin = code; stack[sp++] = fin; while (sp) put(stack[--sp]);
    if (free < (1 << maxBits)) { prefix[free] = old; suffix[free] = fin; free++; }
    old = incode;
    if (free > maxCode() && bits < maxBits) { const span = bits * 8, used = bitPos - groupStart; bitPos = groupStart + Math.ceil(used / span) * span; groupStart = bitPos; bits++; }
  }
  return out.subarray(0, op);
}
