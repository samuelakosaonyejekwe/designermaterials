// Native / proprietary CAD containers: everything that can be read from the container structure without the vendor's
// modelling kernel — OLE2 compound files (property sets, thumbnails, stream inventory), the SOLIDWORKS 2015+ chunk
// container (preview, document properties, configuration names, embedded Parasolid partitions), Siemens NX parts (their
// Parasolid partitions), Inventor, Fusion and FreeCAD archives, the JT table of contents, the Parasolid text header and
// ACIS SAT topology points. Nothing here evaluates a B-rep: each reader states exactly what it recovered, and the
// Parasolid streams it unpacks are translated elsewhere (xt.js).
import * as F from './formats.js';

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
    notRecovered.push('the model segments (B-rep, graphics, browser, design-view data) are stored in Inventor\'s proprietary segment encoding, whose B-rep part is an Autodesk ShapeManager (ACIS-derived) body: not decoded');
  } else {
    if (product === 'Solid Edge') { const found = []; for (const e of d.cfb.entries) if (e.type === 'stream' && e.size > 512) { const data = d.cfb.read(e); if (data) for (const f of findParasolid(data)) if (f.guise !== 'deltas') found.push({ ...f, stream: e.path.replace(/[\x00-\x1f]/g, '·') }); }
      if (found.length) { d.extra.embedded = found.map((f, i) => ({ name: `parasolid-${i + 1}.x_b`, mime: 'application/octet-stream', format: 'parasolid', bytes: f.bytes, description: `Parasolid binary transmit (${f.guise}) from stream "${f.stream}", modeller version ${f.modellerVersion}, schema ${f.schema || 'unknown'}`, encoding: f.encoding, partition: f.partition, modellerVersion: f.modellerVersion, schema: f.schema })); d.recovered.push(`${found.length} embedded Parasolid stream(s) located by signature`); } }
    notRecovered.push(d.extra.embedded ? `the ${product} feature and document streams are in the vendor's proprietary encoding and were not decoded` : `the ${product} model streams are in the vendor's proprietary encoding and were not decoded`);
  }
  if (!d.extra.preview) notRecovered.push('no thumbnail was found in the property sets or preview streams');
  return { meta: d.meta, recovered: d.recovered, notRecovered, extra: d.extra, geometry: null };
}

// ---------------------------------------------------------------- Autodesk Fusion (.f3d / .f3z)
export function parseFusion(bytes) {
  const b = F.toBytes(bytes), all = F.unzip(b, /\.(png|jpe?g)$|(^|\/)(Manifest|Properties)\.dat$/i), recovered = [], notRecovered = [], extra = {}, meta = { entries: all.length };
  const breps = all.filter(e => /\.(smb|smbh|smt|sab|sat)$/i.test(e.name)), imgs = all.filter(e => e.data && /\.(png|jpe?g)$/i.test(e.name) && sniffImage(e.data, 0)).sort((a, c) => c.size - a.size), nested = all.filter(e => /\.f3d$/i.test(e.name));
  meta.archiveEntries = all.slice(0, 60).map(e => `${e.name} (${e.size} bytes)`); recovered.push(`archive listing: ${all.length} entries`);
  if (imgs.length) { extra.preview = preview(sniffImage(imgs[0].data, 0), `archive entry "${imgs[0].name}"`); recovered.push(describePreview(extra.preview)); }
  const man = all.find(e => e.data && /Manifest\.dat$/i.test(e.name)); if (man) { const s = wideStrings(man.data, 4, 40).concat((latin(man.data, 0, min(man.data.length, 65536)).match(/[\x20-\x7e]{6,}/g) || []).slice(0, 40)); if (s.length) meta.manifestStrings = [...new Set(s)].slice(0, 30); }
  if (nested.length) { meta.components = nested.map(e => e.name); recovered.push(`${nested.length} nested design archive(s) (.f3z assembly): ${nested.slice(0, 8).map(e => e.name).join(', ')}`); }
  if (breps.length) { meta.brepSegments = breps.map(e => `${e.name} (${e.size} bytes)`); recovered.push(`${breps.length} B-rep segment(s) located: ${breps.slice(0, 6).map(e => e.name).join(', ')}`); notRecovered.push('the B-rep segments are Autodesk ShapeManager binary (.smb / .smbh, an ACIS derivative): they are located and sized, but no open kernel reads them, so no geometry was reconstructed'); }
  else notRecovered.push('no B-rep segment entries were recognised in the archive');
  notRecovered.push('the design history (timeline, sketches, parameters, joints) is stored in Fusion\'s proprietary serialisation and was not interpreted');
  return { meta, recovered, notRecovered, extra, geometry: null };
}

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

// ---------------------------------------------------------------- JT (ISO 14306)
const JT_SEG = { 1: 'logical scene graph', 2: 'JT B-rep', 3: 'PMI data', 4: 'meta data', 6: 'shape', 7: 'shape LOD 0', 8: 'shape LOD 1', 9: 'shape LOD 2', 10: 'shape LOD 3', 11: 'shape LOD 4', 12: 'shape LOD 5', 13: 'shape LOD 6', 14: 'shape LOD 7', 15: 'shape LOD 8', 16: 'shape LOD 9', 17: 'XT B-rep (Parasolid)', 18: 'wireframe representation', 20: 'ULP', 24: 'LWPA', 30: 'MultiXT B-rep', 31: 'InfoSegment', 33: 'STEP B-rep' };
export function parseJT(bytes) {
  const b = F.toBytes(bytes), v = dv(b), head = latin(b, 0, 80), m = /^Version (\d+)\.(\d+)\s*JT/.exec(head); if (!m) throw new Error('JT: missing "Version n.n JT" header');
  const major = +m[1], meta = { version: `${m[1]}.${m[2]}`, header: head.trim().replace(/\s+/g, ' ') }, recovered = [`JT version ${meta.version}`], notRecovered = [], le = b[80] === 0; meta.byteOrder = le ? 'little-endian' : 'big-endian';
  try {
    const wide = major >= 10, toc = wide ? Number(v.getBigUint64(85, le)) : v.getInt32(85, le); if (!(toc > 0 && toc + 4 <= b.length)) throw new Error('table of contents offset outside the file');
    const n = v.getInt32(toc, le), es = wide ? 16 + 8 + 4 + 4 : 16 + 4 + 4 + 4; if (!(n >= 0 && n < 1e6 && toc + 4 + n * es <= b.length)) throw new Error('table of contents does not fit in the file');
    const counts = {}, sizes = {}; let lods = 0, tri = 0;
    for (let i = 0; i < n; i++) { const o = toc + 4 + i * es, len = v.getInt32(o + (wide ? 24 : 20), le), attr = v.getUint32(o + (wide ? 28 : 24), le), type = attr >>> 24, name = JT_SEG[type] || `segment type ${type}`; counts[name] = (counts[name] || 0) + 1; sizes[name] = (sizes[name] || 0) + max(0, len); if (type >= 6 && type <= 16) { lods++; tri += max(0, len); } }
    meta.segments = n; meta.segmentCounts = counts; meta.segmentBytes = sizes;
    recovered.push(`table of contents: ${n} data segments — ${Object.entries(counts).map(([k, c]) => `${c} × ${k} (${sizes[k]} bytes)`).join(', ')}`);
    if (lods) notRecovered.push(`${lods} tessellated shape / LOD segment(s) (${tri} bytes) are present but were not decoded: JT shape data use quantised, entropy-coded vertex and topology streams whose decoder is not implemented here`);
    if (counts[JT_SEG[2]] || counts[JT_SEG[17]] || counts[JT_SEG[30]] || counts[JT_SEG[33]]) notRecovered.push('the exact B-rep segments (JT B-rep / Parasolid XT / STEP) were inventoried but not extracted or evaluated');
  } catch (e) { notRecovered.push(`the table of contents could not be read (${e.message}): only the version string is reported`); }
  return { meta, recovered, notRecovered, extra: {}, geometry: null };
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
