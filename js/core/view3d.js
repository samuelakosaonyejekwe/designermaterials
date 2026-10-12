// Lightweight canvas viewer for every geometry representation (atoms, lattice, mesh, voxels, layout, rods).
// Drag to rotate, wheel / pinch to zoom, plays `frames` so geometry can evolve with simulation time.
import { h, icon } from './ui.js';
import { elementColor, ELEMENTS } from './chem.js';

const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim() || '#888';

function voxelFaces(g, cap = 26000) {
  const { nx, ny, nz, data } = g, step = Math.max(1, Math.ceil(Math.cbrt(nx * ny * nz / 260000))), at = (x, y, z) => x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz ? 0 : data[x + nx * (y + ny * z)];
  const quads = [], D = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
  const C = { 0: [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], 1: [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0]], 2: [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]], 3: [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]], 4: [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], 5: [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]] };
  for (let z = 0; z < nz; z += step) for (let y = 0; y < ny; y += step) for (let x = 0; x < nx; x += step) { const p = at(x, y, z); if (!p) continue; for (let f = 0; f < 6; f++) { if (at(x + D[f][0] * step, y + D[f][1] * step, z + D[f][2] * step)) continue; quads.push({ p: C[f].map(c => [x + c[0] * step, y + c[1] * step, z + c[2] * step]), n: D[f], ph: p }); if (quads.length >= cap) return quads; } }
  return quads;
}

export function viewer(geom, { height = 340 } = {}) {
  const cv = h('canvas.v3d', { style: { height: height + 'px' } }), info = h('span.muted'), wrap = h('div.viewer');
  let g = geom, rx = -1.05, rz = 0.6, zoom = 1, frame = -1, timer = null, prims = null, center = [0, 0, 0], radius = 1, flat = false;
  const frames = geom.frames || [];
  const merged = () => frame >= 0 && frames[frame] ? { ...geom, ...frames[frame] } : geom;

  function build() {
    g = merged(); const P = []; let pts = [];
    if (g.kind === 'atoms') { pts = g.atoms.map(a => [a.x, a.y, a.z]); const n = g.atoms.length, stride = Math.max(1, Math.ceil(n / 9000)); for (let i = 0; i < n; i += stride) { const a = g.atoms[i]; P.push({ t: 's', c: [a.x, a.y, a.z], r: 0.42 * (ELEMENTS[a.el]?.radius || 1.2), col: elementColor(a.el) }); } if (g.cell) { const [a, b, c] = g.cell, o = [0, 0, 0], ad = (u, v) => u.map((x, i) => x + v[i]), V = [o, a, b, c, ad(a, b), ad(a, c), ad(b, c), ad(ad(a, b), c)]; for (const [i, j] of [[0, 1], [0, 2], [0, 3], [1, 4], [1, 5], [2, 4], [2, 6], [3, 5], [3, 6], [4, 7], [5, 7], [6, 7]]) P.push({ t: 'l', a: V[i], b: V[j], w: 1, col: css('--muted') }); pts = pts.concat(V); } info.textContent = `${n} atoms${stride > 1 ? ` (showing 1 in ${stride})` : ''}`; }
    else if (g.kind === 'lattice') { pts = g.nodes; const w = Math.max(1.2, Math.min(5, 220 / Math.sqrt(g.edges.length + 1))); for (const [i, j] of g.edges) P.push({ t: 'l', a: g.nodes[i], b: g.nodes[j], w, col: css('--s1') }); if (g.ghost) for (const [i, j] of g.edges) P.push({ t: 'l', a: g.ghost[i], b: g.ghost[j], w: 1, col: css('--muted') }); info.textContent = `${g.nodes.length} nodes · ${g.edges.length} struts`; flat = g.nodes.every(p => !p[2]); }
    else if (g.kind === 'mesh') { pts = g.vertices; const nf = g.faces.length, stride = Math.max(1, Math.ceil(nf / 30000)); for (let i = 0; i < nf; i += stride) { const f = g.faces[i]; P.push({ t: 'p', p: [g.vertices[f[0]], g.vertices[f[1]], g.vertices[f[2]]] }); } info.textContent = `${g.vertices.length} vertices · ${nf} faces`; }
    else if (g.kind === 'voxels') { const q = voxelFaces(g); for (const f of q) P.push({ t: 'p', p: f.p, n: f.n, ph: f.ph }); pts = [[0, 0, 0], [g.nx, g.ny, g.nz]]; info.textContent = `${g.nx}×${g.ny}×${g.nz} voxels`; }
    else if (g.kind === 'layout') { flat = true; for (const pg of g.polygons) { P.push({ t: 'p', p: pg.pts.map(p => [p[0], p[1], 0]), layer: pg.layer, flat: true }); for (const p of pg.pts) pts.push([p[0], p[1], 0]); } info.textContent = `${g.polygons.length} polygons`; }
    else if (g.kind === 'rods') { const [Lx, Ly] = g.area || [1, 1]; P.push({ t: 'p', p: [[0, 0, 0], [Lx, 0, 0], [Lx, Ly, 0], [0, Ly, 0]], sub: true }); pts = [[0, 0, 0], [Lx, Ly, 0]]; const n = g.rods.length, stride = Math.max(1, Math.ceil(n / 2500)); let zmax = 0; for (let i = 0; i < n; i += stride) { const r = g.rods[i], t = (r.tilt || 0) * Math.PI / 180, az = (r.azimuth || 0) * Math.PI / 180, tip = [r.x + r.length * Math.sin(t) * Math.cos(az), r.y + r.length * Math.sin(t) * Math.sin(az), r.length * Math.cos(t)]; P.push({ t: 'l', a: [r.x, r.y, 0], b: tip, wd: r.diameter, col: css('--s3'), cap: true }); zmax = Math.max(zmax, tip[2]); } pts.push([Lx / 2, Ly / 2, Math.max(zmax, (geom.frames || []).reduce((m, f) => Math.max(m, ...(f.rods || []).slice(0, 200).map(r => r.length)), 0))]); info.textContent = `${n} rods`; }
    if (frame < 0 || !prims) { let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (const p of pts) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k] || 0); hi[k] = Math.max(hi[k], p[k] || 0); } if (!pts.length) { lo = [0, 0, 0]; hi = [1, 1, 1]; } center = lo.map((v, k) => (v + hi[k]) / 2); radius = Math.max(1e-12, Math.hypot(...hi.map((v, k) => v - lo[k])) / 2); if (flat && frame < 0) { rx = 0; rz = 0; } }
    prims = P;
  }
  function draw() {
    const dpr = Math.min(2, window.devicePixelRatio || 1), W = cv.clientWidth, H = cv.clientHeight; if (!W) return; if (cv.width !== W * dpr) { cv.width = W * dpr; cv.height = H * dpr; }
    const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    const s = Math.min(W, H) * 0.46 / radius * zoom, cz = Math.cos(rz), sz = Math.sin(rz), cx = Math.cos(rx), sx = Math.sin(rx);
    const pr = p => { const x = p[0] - center[0], y = p[1] - center[1], z = (p[2] || 0) - center[2], x1 = x * cz - y * sz, y1 = x * sz + y * cz, y2 = y1 * cx - z * sx, z2 = y1 * sx + z * cx; return [W / 2 + x1 * s, H / 2 - y2 * s, z2]; };
    const list = prims.map(q => { if (q.t === 's') { const c = pr(q.c); return { q, d: c[2], c }; } if (q.t === 'l') { const a = pr(q.a), b = pr(q.b); return { q, d: (a[2] + b[2]) / 2, a, b }; } const pp = q.p.map(pr); return { q, d: pp.reduce((m, p) => m + p[2], 0) / pp.length - (q.sub ? 1e9 : 0), pp }; });
    list.sort((a, b) => a.d - b.d);
    const base = css('--s1'), sub = css('--surface-2'), line = css('--line'), L = [0.35, 0.5, 0.79];
    for (const it of list) { const q = it.q;
      if (q.t === 's') { const r = Math.max(1.2, q.r * s); ctx.beginPath(); ctx.arc(it.c[0], it.c[1], r, 0, 7); ctx.fillStyle = q.col; ctx.fill(); if (r > 3) { ctx.lineWidth = 0.8; ctx.strokeStyle = 'rgba(0,0,0,.45)'; ctx.stroke(); } }
      else if (q.t === 'l') { ctx.beginPath(); ctx.moveTo(it.a[0], it.a[1]); ctx.lineTo(it.b[0], it.b[1]); ctx.lineWidth = q.wd ? Math.max(1, q.wd * s) : q.w; ctx.lineCap = q.cap ? 'round' : 'butt'; ctx.strokeStyle = q.col; ctx.stroke(); }
      else { ctx.beginPath(); it.pp.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.closePath();
        if (q.sub) { ctx.fillStyle = sub; ctx.fill(); ctx.strokeStyle = line; ctx.lineWidth = 1; ctx.stroke(); continue; }
        if (q.flat) { ctx.fillStyle = css('--s' + (((q.layer || 1) - 1) % 8 + 1)); ctx.globalAlpha = 0.78; ctx.fill(); ctx.globalAlpha = 1; continue; }
        let n = q.n; if (!n) { const [a, b, c] = q.p, u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]]; n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]; const m = Math.hypot(...n) || 1; n = n.map(x => x / m); }
        const nx1 = n[0] * cz - n[1] * sz, ny1 = n[0] * sz + n[1] * cz, ny2 = ny1 * cx - n[2] * sx, nz2 = ny1 * sx + n[2] * cx, sh = 0.34 + 0.66 * Math.abs(nx1 * L[0] + ny2 * L[1] + nz2 * L[2]);
        ctx.fillStyle = q.ph > 1 ? css('--s' + ((q.ph - 1) % 8 + 1)) : base; ctx.globalAlpha = 1; ctx.fill(); ctx.fillStyle = `rgba(0,0,0,${(1 - sh) * 0.62})`; ctx.fill(); }
    }
  }
  const redraw = () => requestAnimationFrame(draw);
  let drag = null, pinch = 0; const ptrs = new Map();
  cv.addEventListener('pointerdown', e => { cv.setPointerCapture(e.pointerId); ptrs.set(e.pointerId, [e.clientX, e.clientY]); drag = [e.clientX, e.clientY]; });
  cv.addEventListener('pointermove', e => { if (!ptrs.has(e.pointerId)) return; ptrs.set(e.pointerId, [e.clientX, e.clientY]); if (ptrs.size === 2) { const [a, b] = [...ptrs.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (pinch) zoom = Math.min(40, Math.max(0.2, zoom * d / pinch)); pinch = d; redraw(); return; } if (!drag) return; rz += (e.clientX - drag[0]) * 0.01; rx += (e.clientY - drag[1]) * 0.01; drag = [e.clientX, e.clientY]; redraw(); });
  const up = e => { ptrs.delete(e.pointerId); drag = null; pinch = 0; }; cv.addEventListener('pointerup', up); cv.addEventListener('pointercancel', up);
  cv.addEventListener('wheel', e => { e.preventDefault(); zoom = Math.min(40, Math.max(0.2, zoom * (e.deltaY < 0 ? 1.12 : 0.89))); redraw(); }, { passive: false });
  const bar = h('div.vbar', info);
  if (frames.length) {
    const tl = h('span.muted'), sl = h('input', { type: 'range', min: -1, max: frames.length - 1, value: frames.length - 1, 'aria-label': 'Simulation time' });
    const set = i => { frame = i; sl.value = i; tl.textContent = i < 0 ? 'final' : `t = ${Number(frames[i].t).toPrecision(4)}${geom.timeUnit ? ' ' + geom.timeUnit : ''}`; build(); redraw(); };
    sl.addEventListener('input', () => set(+sl.value));
    const btn = h('button.btn.sm', { onclick: () => { if (timer) { clearInterval(timer); timer = null; btn.replaceChildren(icon('play', 14), ' Play growth'); return; } btn.replaceChildren('Pause'); let i = 0; timer = setInterval(() => { set(i); if (++i >= frames.length) { clearInterval(timer); timer = null; btn.replaceChildren(icon('play', 14), ' Play growth'); } }, 160); } }, icon('play', 14), ' Play growth');
    bar.append(btn, sl, tl); frame = frames.length - 1;
  }
  bar.append(h('button.btn.sm.ghost', { title: 'Reset view', onclick: () => { rx = flat ? 0 : -1.05; rz = flat ? 0 : 0.6; zoom = 1; redraw(); } }, 'Reset view'));
  wrap.append(cv, bar, h('p.hint', 'Drag to rotate · scroll or pinch to zoom' + (geom.units ? ` · units: ${geom.units}` : geom.kind === 'atoms' ? ' · units: Å' : '')));
  build(); new ResizeObserver(redraw).observe(cv); setTimeout(redraw, 0);
  return wrap;
}
