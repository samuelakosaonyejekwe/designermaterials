// SVG / canvas charts: line, scatter, bar (tornado), heat map. Colours come from CSS custom properties.
import { h, svg } from './ui.js';
import { fmt } from './num.js';

const NS = 8;
const col = i => `var(--s${(i % NS) + 1})`;
function ticks(lo, hi, n = 5) {
  if (!(hi > lo)) { hi = lo + (Math.abs(lo) || 1); }
  const span = hi - lo, step0 = 10 ** Math.floor(Math.log10(span / n)), err = span / n / step0, step = step0 * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1), out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}
function logTicks(lo, hi) { const out = []; for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) out.push(10 ** e); return out.filter(v => v >= lo * 0.999 && v <= hi * 1.001); }
const tickFmt = v => { const a = Math.abs(v); if (v === 0) return '0'; if (a >= 1e4 || a < 1e-2) { const e = Math.floor(Math.log10(a)), m = v / 10 ** e; return (Math.abs(m - 1) < 1e-9 ? '' : +m.toPrecision(2) + '×') + '10' + String(e).replace(/[-\d]/g, c => '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'['-0123456789'.indexOf(c)]); } return String(+v.toPrecision(4)); };

export function legend(names) { return names.length < 2 ? null : h('div.legend', names.map((n, i) => h('span.lg', h('i', { style: { background: col(i) } }), n))); }

// spec: {title,xlabel,ylabel,logx,logy,lines:[{name,x,y,points,dash}]}
export function lineChart(spec, { height = 260 } = {}) {
  const W = 640, H = height, m = { l: 58, r: 14, t: 10, b: 40 }, pw = W - m.l - m.r, ph = H - m.t - m.b;
  const lines = (spec.lines || []).filter(l => l.x && l.x.length);
  const pos = (a, log) => log ? Array.from(a).filter(v => v > 0) : Array.from(a);
  let xs = lines.flatMap(l => pos(l.x, spec.logx)), ys = lines.flatMap(l => pos(l.y, spec.logy));
  if (!xs.length || !ys.length) return h('div.chart.empty', 'No plottable data');
  let x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  if (x0 === x1) { x0 -= 1; x1 += 1; }
  if (!spec.logy) { const pad = (y1 - y0) * 0.06 || Math.abs(y1) * 0.1 || 1; y0 -= pad; y1 += pad; if (Math.min(...ys) >= 0 && y0 < 0) y0 = 0; } else if (y0 === y1) { y0 /= 2; y1 *= 2; }
  const tx = v => spec.logx ? Math.log10(v) : v, ty = v => spec.logy ? Math.log10(v) : v;
  const X = v => m.l + (tx(v) - tx(x0)) / (tx(x1) - tx(x0)) * pw, Y = v => m.t + ph - (ty(v) - ty(y0)) / (ty(y1) - ty(y0)) * ph;
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot', role: 'img', 'aria-label': spec.title || 'chart' });
  const g = svg('g', { class: 'axes' });
  for (const t of spec.logy ? logTicks(y0, y1) : ticks(y0, y1, 5)) g.append(svg('line', { x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: 'grid' }), svg('text', { x: m.l - 7, y: Y(t) + 4, 'text-anchor': 'end' }, tickFmt(t)));
  for (const t of spec.logx ? logTicks(x0, x1) : ticks(x0, x1, 6)) g.append(svg('line', { x1: X(t), x2: X(t), y1: m.t + ph, y2: m.t + ph + 4, class: 'tick' }), svg('text', { x: X(t), y: m.t + ph + 17, 'text-anchor': 'middle' }, tickFmt(t)));
  if (!spec.logy && y0 < 0 && y1 > 0) g.append(svg('line', { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0), class: 'zero' }));
  g.append(svg('line', { x1: m.l, x2: W - m.r, y1: m.t + ph, y2: m.t + ph, class: 'axis' }));
  g.append(svg('text', { x: m.l + pw / 2, y: H - 5, 'text-anchor': 'middle', class: 'lab' }, spec.xlabel || ''));
  g.append(svg('text', { x: 13, y: m.t + ph / 2, 'text-anchor': 'middle', class: 'lab', transform: `rotate(-90 13 ${m.t + ph / 2})` }, spec.ylabel || ''));
  root.append(g);
  lines.forEach((l, i) => {
    const ok = k => Number.isFinite(l.x[k]) && Number.isFinite(l.y[k]) && (!spec.logx || l.x[k] > 0) && (!spec.logy || l.y[k] > 0);
    if (l.points) { for (let k = 0; k < l.x.length; k++) if (ok(k)) root.append(svg('circle', { cx: X(l.x[k]), cy: Y(l.y[k]), r: 3.6, fill: col(i), class: 'pt' })); return; }
    let d = '', pen = false; for (let k = 0; k < l.x.length; k++) { if (!ok(k)) { pen = false; continue; } d += (pen ? 'L' : 'M') + X(l.x[k]).toFixed(1) + ' ' + Y(l.y[k]).toFixed(1); pen = true; }
    root.append(svg('path', { d, fill: 'none', stroke: col(i), 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-dasharray': l.dash ? '6 4' : null }));
  });
  // hover layer
  const cross = svg('line', { y1: m.t, y2: m.t + ph, class: 'cross', visibility: 'hidden' }), dots = lines.map((_, i) => svg('circle', { r: 4.5, fill: col(i), class: 'hov', visibility: 'hidden' }));
  root.append(cross, ...dots);
  const tip = h('div.tip', { hidden: true }), wrap = h('div.chart', root, tip, legend(lines.map(l => l.name || '')));
  const move = ev => {
    const r = root.getBoundingClientRect(), px = ((ev.touches ? ev.touches[0].clientX : ev.clientX) - r.left) / r.width * W; if (px < m.l || px > W - m.r) return leave();
    const xv = spec.logx ? 10 ** (tx(x0) + (px - m.l) / pw * (tx(x1) - tx(x0))) : x0 + (px - m.l) / pw * (x1 - x0); let xn = null; const rowsT = [];
    lines.forEach((l, i) => { let best = 0, bd = Infinity; for (let k = 0; k < l.x.length; k++) { const dd = Math.abs(tx(Math.max(l.x[k], spec.logx ? 1e-300 : -Infinity)) - tx(xv)); if (dd < bd) { bd = dd; best = k; } } const yv = l.y[best]; if (!Number.isFinite(yv) || (spec.logy && yv <= 0)) { dots[i].setAttribute('visibility', 'hidden'); return; } if (xn === null) xn = l.x[best]; dots[i].setAttribute('cx', X(l.x[best])); dots[i].setAttribute('cy', Y(yv)); dots[i].setAttribute('visibility', 'visible'); rowsT.push(h('div', h('i', { style: { background: col(i) } }), (l.name ? l.name + ': ' : ''), h('b', fmt(yv)))); });
    if (xn === null) return leave(); cross.setAttribute('x1', X(xn)); cross.setAttribute('x2', X(xn)); cross.setAttribute('visibility', 'visible');
    tip.replaceChildren(h('div.th', `${spec.xlabel || 'x'} = ${fmt(xn)}`), ...rowsT); tip.hidden = false; const left = (X(xn) / W) * r.width; tip.style.left = Math.min(Math.max(left + 12, 4), r.width - 170) + 'px'; tip.style.top = '8px';
  };
  const leave = () => { cross.setAttribute('visibility', 'hidden'); dots.forEach(d => d.setAttribute('visibility', 'hidden')); tip.hidden = true; };
  root.addEventListener('pointermove', move); root.addEventListener('pointerleave', leave); root.addEventListener('touchmove', move, { passive: true });
  return wrap;
}

// Horizontal bar chart (signed), used for tornado / sensitivity rankings. items:[{label,value,lo?,hi?}]
export function barChart(items, { unit = '', signed = true, max = 14 } = {}) {
  items = items.slice(0, max); const W = 640, rowH = 26, m = { l: 190, r: 60, t: 6, b: 6 }, H = m.t + m.b + rowH * items.length, pw = W - m.l - m.r;
  const ext = Math.max(1e-300, ...items.map(i => Math.max(Math.abs(i.lo ?? i.value), Math.abs(i.hi ?? i.value)))), hasNeg = signed && items.some(i => (i.lo ?? i.value) < 0), zero = hasNeg ? m.l + pw / 2 : m.l, sc = (hasNeg ? pw / 2 : pw) / ext;
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot bars', role: 'img' });
  items.forEach((it, k) => {
    const y = m.t + k * rowH + 4, lab = it.label.length > 30 ? it.label.slice(0, 29) + '…' : it.label;
    root.append(svg('text', { x: m.l - 8, y: y + 13, 'text-anchor': 'end' }, lab));
    const seg = (v, c) => { const w = Math.abs(v) * sc; const r = svg('rect', { x: v < 0 ? zero - w : zero, y, width: Math.max(w, 1), height: rowH - 9, rx: 3, fill: c }); r.append(svg('title', {}, `${it.label}: ${fmt(v)} ${unit}`)); root.append(r); };
    if (it.lo !== undefined) { seg(it.lo, 'var(--s1)'); seg(it.hi, 'var(--s2)'); } else seg(it.value, it.value < 0 ? 'var(--s2)' : 'var(--s1)');
    root.append(svg('text', { x: W - m.r + 6, y: y + 13, class: 'val' }, fmt(it.lo !== undefined ? Math.max(Math.abs(it.lo), Math.abs(it.hi)) : it.value, 3)));
  });
  root.append(svg('line', { x1: zero, x2: zero, y1: m.t, y2: H - m.b, class: 'axis' }));
  return h('div.chart', root);
}

const SEQ = [[205, 226, 251], [134, 182, 239], [57, 135, 229], [28, 92, 171], [13, 54, 107]], DIV = [[42, 120, 214], [240, 239, 236], [227, 73, 72]];
function ramp(stops, t) { t = Math.min(1, Math.max(0, t)) * (stops.length - 1); const i = Math.min(stops.length - 2, Math.floor(t)), f = t - i; return stops[i].map((c, k) => Math.round(c + (stops[i + 1][k] - c) * f)); }
// spec: {nx,ny,data,x0,x1,y0,y1,xlabel,ylabel,unit,diverging}
export function heatmap(spec, { height = 260 } = {}) {
  const { nx, ny, data } = spec; let lo = Infinity, hi = -Infinity; for (const v of data) { if (v < lo) lo = v; if (v > hi) hi = v; }
  const div = spec.diverging ?? (lo < 0 && hi > 0), ext = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
  const cv = h('canvas', { width: nx, height: ny, class: 'heat' }), ctx = cv.getContext('2d'), img = ctx.createImageData(nx, ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const v = data[i + nx * j], c = div ? ramp(DIV, 0.5 + 0.5 * v / ext) : ramp(SEQ, hi > lo ? (v - lo) / (hi - lo) : 0.5), p = 4 * (i + nx * (ny - 1 - j)); img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255; }
  ctx.putImageData(img, 0, 0);
  const tip = h('div.tip', { hidden: true }), box = h('div.heatbox', { style: { aspectRatio: `${Math.max(0.4, Math.min(3, nx / ny))}`, maxHeight: height + 'px' } }, cv, tip);
  const xr = [spec.x0 ?? 0, spec.x1 ?? nx], yr = [spec.y0 ?? 0, spec.y1 ?? ny];
  box.addEventListener('pointermove', ev => { const r = cv.getBoundingClientRect(), i = Math.min(nx - 1, Math.max(0, Math.floor((ev.clientX - r.left) / r.width * nx))), j = Math.min(ny - 1, Math.max(0, ny - 1 - Math.floor((ev.clientY - r.top) / r.height * ny))); tip.hidden = false; tip.replaceChildren(h('div.th', `${spec.xlabel || 'x'} ${fmt(xr[0] + (i + 0.5) / nx * (xr[1] - xr[0]), 3)}, ${spec.ylabel || 'y'} ${fmt(yr[0] + (j + 0.5) / ny * (yr[1] - yr[0]), 3)}`), h('div', h('b', fmt(data[i + nx * j])), ' ', spec.unit || '')); tip.style.left = Math.min(ev.clientX - r.left + 12, r.width - 150) + 'px'; tip.style.top = Math.max(ev.clientY - r.top - 44, 0) + 'px'; });
  box.addEventListener('pointerleave', () => tip.hidden = true);
  const grad = div ? `linear-gradient(90deg, rgb(${DIV[0]}), rgb(${DIV[1]}), rgb(${DIV[2]}))` : `linear-gradient(90deg, ${SEQ.map(c => `rgb(${c})`).join(',')})`;
  return h('div.chart', box, h('div.cbar', h('span', fmt(div ? -ext : lo, 3)), h('i', { style: { background: grad } }), h('span', fmt(div ? ext : hi, 3) + (spec.unit ? ' ' + spec.unit : ''))), h('div.axnote', `${spec.xlabel || 'x'}: ${fmt(xr[0], 3)} → ${fmt(xr[1], 3)}   ·   ${spec.ylabel || 'y'}: ${fmt(yr[0], 3)} → ${fmt(yr[1], 3)}`));
}
export function histChart(values, { bins = 24, xlabel = '', height = 220 } = {}) {
  let lo = Math.min(...values), hi = Math.max(...values); if (lo === hi) hi = lo + 1; const w = (hi - lo) / bins, c = new Array(bins).fill(0); for (const v of values) c[Math.min(bins - 1, Math.floor((v - lo) / w))]++;
  const W = 640, H = height, m = { l: 48, r: 12, t: 8, b: 38 }, pw = W - m.l - m.r, ph = H - m.t - m.b, mx = Math.max(...c), root = svg('svg', { viewBox: `0 0 ${W} ${H}`, class: 'plot', role: 'img' });
  c.forEach((n, i) => { const bh = n / mx * ph, r = svg('rect', { x: m.l + i * pw / bins + 1, y: m.t + ph - bh, width: pw / bins - 2, height: bh, rx: 2, fill: 'var(--s1)' }); r.append(svg('title', {}, `${fmt(lo + i * w, 3)} – ${fmt(lo + (i + 1) * w, 3)}: ${n}`)); root.append(r); });
  for (const t of ticks(lo, hi, 6)) root.append(svg('text', { x: m.l + (t - lo) / (hi - lo) * pw, y: m.t + ph + 16, 'text-anchor': 'middle' }, tickFmt(t)));
  root.append(svg('line', { x1: m.l, x2: W - m.r, y1: m.t + ph, y2: m.t + ph, class: 'axis' }), svg('text', { x: m.l + pw / 2, y: H - 4, 'text-anchor': 'middle', class: 'lab' }, xlabel), svg('text', { x: m.l - 6, y: m.t + 10, 'text-anchor': 'end' }, String(mx)));
  return h('div.chart', root);
}
