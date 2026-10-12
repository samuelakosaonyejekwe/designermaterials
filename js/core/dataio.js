// Universal tabular / structured data readers for the data portal and table inputs:
// CSV, TSV, whitespace tables, JSON, NDJSON, XML, YAML (flat), XLSX (first sheets), instrument text files.
import { parseCSV } from './ui.js';

const td = new TextDecoder();
async function inflateRaw(bytes) { const ds = new DecompressionStream('deflate-raw'), w = ds.writable.getWriter(); w.write(bytes); w.close(); return new Uint8Array(await new Response(ds.readable).arrayBuffer()); }
export async function unzip(bytes) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength), out = {}; let p = bytes.length - 22; while (p >= 0 && dv.getUint32(p, true) !== 0x06054b50) p--; if (p < 0) throw new Error('Not a ZIP container');
  let off = dv.getUint32(p + 16, true); const n = dv.getUint16(p + 10, true);
  for (let i = 0; i < n; i++) { const method = dv.getUint16(off + 10, true), csz = dv.getUint32(off + 20, true), nl = dv.getUint16(off + 28, true), el = dv.getUint16(off + 30, true), cl = dv.getUint16(off + 32, true), lo = dv.getUint32(off + 42, true), name = td.decode(bytes.subarray(off + 46, off + 46 + nl)), lnl = dv.getUint16(lo + 26, true), lel = dv.getUint16(lo + 28, true), data = bytes.subarray(lo + 30 + lnl + lel, lo + 30 + lnl + lel + csz); out[name] = method === 0 ? data : method === 8 ? await inflateRaw(data) : null; off += 46 + nl + el + cl; }
  return out;
}
const xmlText = s => s.replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&amp;/g, '&');
async function readXLSX(bytes) {
  const z = await unzip(bytes), get = n => z[n] ? td.decode(z[n]) : '', shared = [...get('xl/sharedStrings.xml').matchAll(/<si>([\s\S]*?)<\/si>/g)].map(m => xmlText([...m[1].matchAll(/<t[^>]*>([\s\S]*?)<\/t>/g)].map(t => t[1]).join('')));
  const sheets = Object.keys(z).filter(k => /^xl\/worksheets\/sheet\d+\.xml$/.test(k)).sort((a, b) => parseInt(a.match(/\d+/)[0]) - parseInt(b.match(/\d+/)[0])); if (!sheets.length) throw new Error('No worksheets found in the workbook');
  const names = [...get('xl/workbook.xml').matchAll(/<sheet [^>]*name="([^"]*)"/g)].map(m => xmlText(m[1])), tables = [];
  for (const [k, sh] of sheets.entries()) {
    const grid = []; for (const rm of td.decode(z[sh]).matchAll(/<row[^>]*>([\s\S]*?)<\/row>/g)) { const row = []; for (const c of rm[1].matchAll(/<c ([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)) { const ref = (c[1].match(/r="([A-Z]+)\d+"/) || [])[1] || 'A', t = (c[1].match(/t="(\w+)"/) || [])[1], v = (c[2] || '').match(/<v>([\s\S]*?)<\/v>/), is = (c[2] || '').match(/<t[^>]*>([\s\S]*?)<\/t>/); let col = 0; for (const ch of ref) col = col * 26 + ch.charCodeAt(0) - 64; let val = v ? v[1] : is ? is[1] : ''; val = t === 's' ? shared[+val] : t === 'inlineStr' || t === 'str' ? xmlText(val) : val === '' ? '' : Number(val); row[col - 1] = val; } grid.push(Array.from(row, x => x ?? '')); }
    const rows = grid.filter(r => r.some(c => c !== '')); if (!rows.length) continue; const w = Math.max(...rows.map(r => r.length)), hasHeader = rows[0].some(c => typeof c === 'string' && c !== '');
    tables.push({ name: names[k] || `Sheet ${k + 1}`, columns: hasHeader ? Array.from({ length: w }, (_, i) => String(rows[0][i] ?? `col${i + 1}`)) : Array.from({ length: w }, (_, i) => `col${i + 1}`), rows: rows.slice(hasHeader ? 1 : 0).map(r => Array.from({ length: w }, (_, i) => r[i] ?? '')) });
  }
  return tables;
}
function flatten(o, pre = '', out = {}) { if (o !== null && typeof o === 'object' && !Array.isArray(o)) for (const k in o) flatten(o[k], pre ? pre + '.' + k : k, out); else out[pre] = Array.isArray(o) ? JSON.stringify(o) : o; return out; }
function fromRecords(recs) { const cols = [...new Set(recs.flatMap(r => Object.keys(r)))]; return { columns: cols, rows: recs.map(r => cols.map(c => r[c] ?? '')) }; }
function fromJSON(o) {
  if (Array.isArray(o)) { if (o.every(r => Array.isArray(r))) return parseCSV(o.map(r => r.join(',')).join('\n')); if (o.every(r => r && typeof r === 'object')) return fromRecords(o.map(r => flatten(r))); return { columns: ['value'], rows: o.map(v => [v]) }; }
  const arrs = Object.entries(o).filter(([, v]) => Array.isArray(v) && v.length && v.every(x => typeof x !== 'object'));
  if (arrs.length >= 2) { const n = Math.max(...arrs.map(a => a[1].length)); return { columns: arrs.map(a => a[0]), rows: Array.from({ length: n }, (_, i) => arrs.map(a => a[1][i] ?? '')) }; }
  const rec = Object.values(o).find(v => Array.isArray(v) && v.length && typeof v[0] === 'object'); if (rec) return fromRecords(rec.map(r => flatten(r)));
  const f = flatten(o); return { columns: ['parameter', 'value'], rows: Object.entries(f) };
}
function fromYAML(text) { const rows = []; const stack = []; for (const line of text.split(/\r?\n/)) { if (!line.trim() || /^\s*#/.test(line)) continue; const m = line.match(/^(\s*)([^:#]+):\s*(.*)$/); if (!m) continue; const depth = m[1].length; while (stack.length && stack.at(-1).d >= depth) stack.pop(); const key = m[2].trim(); if (m[3] === '') stack.push({ d: depth, k: key }); else { let v = m[3].replace(/\s+#.*$/, '').replace(/^["']|["']$/g, ''); rows.push([[...stack.map(s => s.k), key].join('.'), v !== '' && !Number.isNaN(Number(v)) ? Number(v) : v]); } } return { columns: ['parameter', 'value'], rows }; }
function fromXML(text) { const rows = []; const walk = (s, path) => { let any = false; for (const m of s.matchAll(/<([\w:.-]+)([^>]*)>([\s\S]*?)<\/\1>/g)) { any = true; const p = path ? path + '.' + m[1] : m[1]; for (const a of m[2].matchAll(/([\w:.-]+)="([^"]*)"/g)) rows.push([p + '@' + a[1], Number.isNaN(Number(a[2])) || a[2] === '' ? xmlText(a[2]) : Number(a[2])]); if (!walk(m[3], p)) { const v = xmlText(m[3].trim()); if (v !== '') rows.push([p, Number.isNaN(Number(v)) ? v : Number(v)]); } } return any; }; walk(text.replace(/<\?[\s\S]*?\?>|<!--[\s\S]*?-->/g, ''), ''); return { columns: ['parameter', 'value'], rows }; }
// Instrument exports often carry a metadata preamble before the numeric block: keep the longest numeric run.
function numericBlock(text) { const lines = text.replace(/\r/g, '').split('\n'); let best = [0, 0], start = -1; const isNum = l => { const p = l.trim().split(/[\s,;\t]+/); return p.length >= 2 && p.every(x => x !== '' && !Number.isNaN(Number(x))); }; lines.forEach((l, i) => { if (isNum(l)) { if (start < 0) start = i; if (i - start + 1 > best[1]) best = [start, i - start + 1]; } else start = -1; }); if (best[1] < 3) return null; const head = best[0] > 0 ? lines[best[0] - 1].trim().split(/[\t,;]+|\s{2,}/).filter(Boolean) : []; const rows = lines.slice(best[0], best[0] + best[1]).map(l => l.trim().split(/[\s,;\t]+/).map(Number)), w = rows[0].length; return { columns: head.length === w && head.some(c => Number.isNaN(Number(c))) ? head : Array.from({ length: w }, (_, i) => w === 2 ? ['x', 'y'][i] : `col${i + 1}`), rows, preamble: lines.slice(0, best[0]).filter(l => l.trim()).slice(0, 40) }; }

export const DATA_FORMATS = [
  ['CSV / TSV / delimited text', '.csv .tsv .tab .txt .dat .asc', 'Tables of measurements, process logs, spectra; delimiter and header detected automatically'],
  ['Instrument XY exports', '.xy .xye .xrdml(text) .dpt .spc(text) .jdx', 'Diffraction, spectroscopy and analyser exports with a metadata preamble; the numeric block is extracted and the preamble kept'],
  ['Excel workbook', '.xlsx', 'Every worksheet is read as a table (values and shared strings; formulas use their cached values)'],
  ['JSON / NDJSON', '.json .ndjson .jsonl', 'Arrays of records, column arrays, or nested parameter sets (flattened to parameter/value)'],
  ['XML', '.xml', 'Structured material, layer-stack and instrument records flattened to parameter/value paths'],
  ['YAML / INI-style', '.yaml .yml .ini .cfg .toml', 'Process recipes and configuration parameter sets'],
  ['Touchstone S-parameters', '.s1p .s2p', 'Frequency-dependent S-parameters for electromagnetic retrieval'],
  ['Images', '.png .jpg .jpeg .bmp .gif .webp .tif .tiff', 'SEM / TEM / AFM / optical micrographs for image-based reconstruction (opened in Geometry Studio)'],
  ['Geometry and structure files', 'see Geometry import', 'CAD, meshes, crystal/atomic structures, lithography layouts, voxel volumes — routed to the geometry importer'],
];
export async function readTabular(name, bytes) {
  const ext = (name.toLowerCase().match(/\.[a-z0-9]+$/) || [''])[0];
  if (ext === '.xlsx' || (bytes[0] === 0x50 && bytes[1] === 0x4b && /\.xls[xm]$/.test(name.toLowerCase()))) { const t = await readXLSX(bytes); return { ...t[0], sheets: t, kind: 'Excel workbook' }; }
  if (ext === '.xls') throw new Error('Legacy binary .xls is not readable in the browser: save the workbook as .xlsx or .csv and import again.');
  const text = td.decode(bytes).replace(/^﻿/, '');
  if (ext === '.json' || /^\s*[\[{]/.test(text)) { try { return { ...fromJSON(JSON.parse(text)), kind: 'JSON' }; } catch (e) { if (ext === '.json') throw new Error('Invalid JSON: ' + e.message); } }
  if (ext === '.ndjson' || ext === '.jsonl') return { ...fromRecords(text.split('\n').filter(l => l.trim()).map(l => flatten(JSON.parse(l)))), kind: 'NDJSON' };
  if (/^\.s\d+p$/.test(ext)) { const opt = (text.match(/^#.*$/m) || [''])[0], rows = text.split(/\r?\n/).filter(l => l.trim() && !/^[!#]/.test(l)).map(l => l.trim().split(/\s+/).map(Number)), n = rows[0].length; return { columns: ['frequency', ...Array.from({ length: n - 1 }, (_, i) => `S${Math.floor(i / 2) < 2 ? ['11', '21', '12', '22'][Math.floor(i / 2)] || i : i}_${i % 2 ? 'b' : 'a'}`)], rows, kind: 'Touchstone', preamble: [opt] }; }
  if (ext === '.xml' || /^\s*<\?xml|^\s*<[\w:]+[\s>]/.test(text)) return { ...fromXML(text), kind: 'XML' };
  if (['.yaml', '.yml', '.ini', '.cfg', '.toml'].includes(ext)) return { ...fromYAML(text.replace(/^\s*\[.*\]\s*$/gm, '').replace(/^(\s*[\w.-]+)\s*=\s*/gm, '$1: ')), kind: 'Parameter file' };
  const nb = numericBlock(text), csv = parseCSV(text);
  if (nb && (nb.preamble.length > 1 || csv.columns.length < 2 || nb.rows.length > csv.rows.length)) return { ...nb, kind: 'Numeric table' };
  if (!csv.rows.length) throw new Error('No tabular data found in this file.'); return { ...csv, kind: 'Delimited text' };
}
// Suggest which bus quantities a parameter/value table can fill (name + unit heuristics)
const ALIASES = [['prop.density', /dens/i], ['prop.E', /young|elastic mod|^e$|modulus of elast/i], ['prop.nu', /poisson/i], ['prop.k_thermal', /therm.*cond|^k$/i], ['prop.cp', /specific heat|heat cap|^cp$/i], ['prop.sigma_elec', /elec.*cond/i], ['prop.bandgap_eV', /band ?gap/i], ['prop.eps_r', /permittiv|dielectric const/i], ['prop.n_refr', /refractive/i], ['geom.porosity', /poros/i], ['geom.diameter', /diamet/i], ['geom.length', /length/i], ['geom.thickness', /thick/i], ['proc.T_C', /temp/i], ['proc.P', /press/i], ['proc.pH', /^ph$/i], ['proc.time', /time|duration/i], ['proc.E_field', /electric field|e-?field/i], ['proc.c_precursor', /zn.*conc|precursor conc/i], ['proc.c_base', /oh.*conc|base conc|hydroxide/i], ['material.formula', /formula|material$/i], ['prop.D_diff', /diffus/i]];
export function suggestBus(name) { const hit = ALIASES.find(([, re]) => re.test(String(name).trim())); return hit ? hit[0] : ''; }
