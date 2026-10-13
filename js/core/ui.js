// Tiny DOM toolkit shared by every page.
// DOM insertion methods accept nested arrays and skip null / undefined / false, so conditional children can be written inline.
if (typeof Element !== 'undefined') for (const m of ['replaceChildren', 'append', 'prepend', 'after', 'before']) { const orig = Element.prototype[m]; Element.prototype[m] = function (...k) { return orig.apply(this, k.flat(Infinity).filter(x => x !== null && x !== undefined && x !== false)); }; }
export function h(tag, attrs, ...kids) {
  const m = tag.match(/^([a-z0-9]+)?((?:[.#][\w-]+)*)$/i) || [0, tag, ''];
  const el = document.createElement(m[1] || 'div');
  (m[2].match(/[.#][\w-]+/g) || []).forEach(t => t[0] === '.' ? el.classList.add(t.slice(1)) : (el.id = t.slice(1)));
  if (attrs && (typeof attrs !== 'object' || attrs instanceof Node || Array.isArray(attrs))) { kids.unshift(attrs); attrs = null; }
  for (const k in attrs || {}) {
    const v = attrs[k]; if (v === null || v === undefined || v === false) continue;
    if (k.startsWith('on')) el.addEventListener(k.slice(2).toLowerCase(), v);
    else if (k === 'class') el.className += ' ' + v;
    else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
    else if (k === 'html') el.innerHTML = v;
    else if (k in el && k !== 'list' && k !== 'form') { try { el[k] = v; } catch { el.setAttribute(k, v); } }
    else el.setAttribute(k, v === true ? '' : v);
  }
  const add = c => { if (c === null || c === undefined || c === false) return; if (Array.isArray(c)) c.forEach(add); else el.append(c instanceof Node ? c : document.createTextNode(String(c))); };
  kids.forEach(add);
  return el;
}
export const svg = (tag, attrs = {}, ...kids) => { const el = document.createElementNS('http://www.w3.org/2000/svg', tag); for (const k in attrs) if (attrs[k] !== null && attrs[k] !== undefined) el.setAttribute(k, attrs[k]); kids.flat().forEach(c => c !== null && c !== undefined && el.append(c instanceof Node ? c : document.createTextNode(String(c)))); return el; };
export const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
export function toast(msg, kind = 'info', ms = 3600) {
  let box = document.getElementById('toasts'); if (!box) { box = h('div#toasts', { role: 'status', 'aria-live': 'polite' }); document.body.append(box); }
  const t = h('div.toast.' + kind, msg); box.append(t); setTimeout(() => { t.classList.add('out'); setTimeout(() => t.remove(), 300); }, ms);
}
export function download(name, content, mime = 'application/octet-stream') {
  const blob = content instanceof Blob ? content : new Blob([content], { type: mime });
  const a = h('a', { href: URL.createObjectURL(blob), download: name }); document.body.append(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 800);
}
export const icon = (name, size = 20) => { const p = ICONS[name] || ICONS.dot; const s = svg('svg', { viewBox: '0 0 24 24', width: size, height: size, fill: 'none', stroke: 'currentColor', 'stroke-width': 1.8, 'stroke-linecap': 'round', 'stroke-linejoin': 'round', 'aria-hidden': 'true', class: 'ico' }); s.innerHTML = p; return s; };
const ICONS = {
  dot: '<circle cx="12" cy="12" r="3"/>',
  back: '<path d="M15 5l-7 7 7 7"/>', fwd: '<path d="M9 5l7 7-7 7"/>', home: '<path d="M3 11l9-8 9 8M5 10v10h5v-6h4v6h5V10"/>',
  grid: '<rect x="3" y="3" width="7" height="7" rx="1.5"/><rect x="14" y="3" width="7" height="7" rx="1.5"/><rect x="3" y="14" width="7" height="7" rx="1.5"/><rect x="14" y="14" width="7" height="7" rx="1.5"/>',
  book: '<path d="M4 5a2 2 0 012-2h13v16H6a2 2 0 00-2 2zM4 21V5M9 7h6"/>', flask: '<path d="M9 3h6M10 3v6L4.5 19a1.5 1.5 0 001.3 2.2h12.4a1.5 1.5 0 001.3-2.2L14 9V3M7.5 15h9"/>',
  cube: '<path d="M12 2l9 5v10l-9 5-9-5V7zM3 7l9 5 9-5M12 12v10"/>', upload: '<path d="M12 16V4M7 9l5-5 5 5M4 20h16"/>', download: '<path d="M12 4v12M7 11l5 5 5-5M4 20h16"/>',
  flow: '<circle cx="5" cy="6" r="2.5"/><circle cx="19" cy="6" r="2.5"/><circle cx="12" cy="18" r="2.5"/><path d="M7.5 6h9M6.2 8.2l4.6 7.6M17.8 8.2l-4.6 7.6"/>',
  globe: '<circle cx="12" cy="12" r="9"/><path d="M3 12h18M12 3a14 14 0 010 18M12 3a14 14 0 000 18"/>', doc: '<path d="M6 2h8l5 5v15H6zM14 2v5h5M9 13h7M9 17h7"/>',
  gear: '<circle cx="12" cy="12" r="3"/><path d="M19.4 15a1.7 1.7 0 00.3 1.8l.1.1a2 2 0 11-2.8 2.8l-.1-.1a1.7 1.7 0 00-2.8 1.2V21a2 2 0 11-4 0v-.2a1.7 1.7 0 00-2.9-1.1l-.1.1a2 2 0 11-2.8-2.8l.1-.1a1.7 1.7 0 00-1.2-2.9H3a2 2 0 110-4h.2a1.7 1.7 0 001.1-2.9l-.1-.1a2 2 0 112.8-2.8l.1.1a1.7 1.7 0 002.9-1.2V3a2 2 0 114 0v.2a1.7 1.7 0 002.9 1.1l.1-.1a2 2 0 112.8 2.8l-.1.1a1.7 1.7 0 001.2 2.9H21a2 2 0 110 4h-.2a1.7 1.7 0 00-1.4 1z"/>',
  play: '<path d="M7 4l13 8-13 8z"/>', search: '<circle cx="11" cy="11" r="7"/><path d="M21 21l-4.5-4.5"/>', sun: '<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M2 12h2M20 12h2M5 5l1.5 1.5M17.5 17.5L19 19M19 5l-1.5 1.5M6.5 17.5L5 19"/>',
  menu: '<path d="M4 6h16M4 12h16M4 18h16"/>', x: '<path d="M6 6l12 12M18 6L6 18"/>', check: '<path d="M5 12l5 5 9-10"/>', warn: '<path d="M12 3l10 18H2zM12 10v5M12 18v.5"/>', info: '<circle cx="12" cy="12" r="9"/><path d="M12 11v6M12 7.5v.5"/>',
  link: '<path d="M10 14a5 5 0 007 0l3-3a5 5 0 00-7-7l-1.5 1.5M14 10a5 5 0 00-7 0l-3 3a5 5 0 007 7L12.5 18.5"/>', leaf: '<path d="M5 19c0-9 5-14 15-14 0 10-5 15-14 15M5 19c2-5 5-8 9-10"/>', atom: '<circle cx="12" cy="12" r="1.6"/><ellipse cx="12" cy="12" rx="10" ry="4"/><ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(60 12 12)"/><ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(120 12 12)"/>',
  chart: '<path d="M4 20V4M4 20h16M8 16v-5M12 16V8M16 16v-8"/>', layers: '<path d="M12 3l9 5-9 5-9-5zM3 13l9 5 9-5M3 17l9 5 9-5"/>', refresh: '<path d="M20 11a8 8 0 10-2.3 5.7M20 4v7h-7"/>', install: '<path d="M12 3v11M7.5 10l4.5 4.5 4.5-4.5M5 21h14a2 2 0 002-2v-3M3 16v3a2 2 0 002 2"/>',
  wifi: '<path d="M2 9a15 15 0 0120 0M5.5 12.5a10 10 0 0113 0M9 16a5 5 0 016 0"/><circle cx="12" cy="19.5" r="1"/>', star: '<path d="M12 3l2.8 5.7 6.2.9-4.5 4.4 1.1 6.2L12 17.3 6.4 20.2l1.1-6.2L3 9.6l6.2-.9z"/>', target: '<circle cx="12" cy="12" r="9"/><circle cx="12" cy="12" r="5"/><circle cx="12" cy="12" r="1"/>', swap: '<path d="M7 4L3 8l4 4M3 8h14M17 20l4-4-4-4M21 16H7"/>',
};
export function debounce(fn, ms = 200) { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }
export const fileBytes = f => new Promise((res, rej) => { const r = new FileReader(); r.onload = () => res(new Uint8Array(r.result)); r.onerror = () => rej(r.error); r.readAsArrayBuffer(f); });
export function tabs(defs, { initial = 0, onChange } = {}) {
  const bar = h('div.tabs', { role: 'tablist' }), body = h('div.tabbody'); let cur = -1; const cache = [];
  const show = i => { if (i === cur) return; cur = i;[...bar.children].forEach((b, j) => { b.classList.toggle('on', j === i); b.setAttribute('aria-selected', j === i); }); if (!cache[i] || defs[i].fresh) cache[i] = h('div.tabpane', defs[i].render()); body.replaceChildren(cache[i]); onChange && onChange(i); };
  defs.forEach((d, i) => bar.append(h('button.tab', { role: 'tab', onclick: () => show(i) }, d.label, d.count != null ? h('span.count', d.count) : null)));
  show(initial); return h('div.tabset', bar, body);
}
export function table(columns, rows, { cls = '', max = 400 } = {}) {
  const fmtc = v => typeof v === 'number' ? (Number.isInteger(v) && Math.abs(v) < 1e6 ? String(v) : (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(3) : String(+v.toPrecision(5)))) : (v ?? '—');
  return h('div.tablewrap', h('table.' + (cls || 'data'), h('thead', h('tr', columns.map(c => h('th', c)))), h('tbody', rows.slice(0, max).map(r => h('tr', r.map(c => h('td', { class: typeof c === 'number' ? 'num' : '' }, c instanceof Node ? c : fmtc(c))))))), rows.length > max ? h('p.muted', `Showing ${max} of ${rows.length} rows — export for the full table.`) : null);
}
export function toCSV(columns, rows) { const q = v => { const s = v === null || v === undefined ? '' : String(v); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }; return [columns.map(q).join(','), ...rows.map(r => r.map(q).join(','))].join('\n'); }
export function parseCSV(text) {
  const lines = text.replace(/\r/g, '').split('\n').filter(l => l.trim() && !/^\s*[#%]/.test(l)); if (!lines.length) return { columns: [], rows: [] };
  const delim = [',', '\t', ';', '|'].map(d => [d, lines[0].split(d).length]).sort((a, b) => b[1] - a[1])[0]; const d = delim[1] > 1 ? delim[0] : /\s+/;
  const split = l => { if (typeof d !== 'string') return l.trim().split(d); const out = []; let cur = '', q = false; for (let i = 0; i < l.length; i++) { const c = l[i]; if (q) { if (c === '"' && l[i + 1] === '"') { cur += '"'; i++; } else if (c === '"') q = false; else cur += c; } else if (c === '"') q = true; else if (c === d) { out.push(cur); cur = ''; } else cur += c; } out.push(cur); return out.map(s => s.trim()); };
  const first = split(lines[0]), hasHeader = first.some(c => c !== '' && Number.isNaN(Number(c)));
  const columns = hasHeader ? first : first.map((_, i) => 'col' + (i + 1));
  const rows = lines.slice(hasHeader ? 1 : 0).map(l => split(l).map(c => c !== '' && !Number.isNaN(Number(c)) ? Number(c) : c));
  return { columns, rows };
}

// Background worker for long computations. In the single-file edition the worker is started from the page's own bundle.
export function makeWorker() {
  if (typeof Worker === 'undefined') throw new Error('Background computation is not available in this browser.');
  if (globalThis.__DM_STANDALONE__) return new Worker(URL.createObjectURL(new Blob([document.getElementById('dm-bundle').textContent], { type: 'text/javascript' })), { type: 'module' });
  const u = new URL('./worker.js', import.meta.url);
  return u.origin === location.origin ? new Worker(u, { type: 'module' }) : new Worker(URL.createObjectURL(new Blob([`import ${JSON.stringify(u.href)};`], { type: 'text/javascript' })), { type: 'module' });
}
