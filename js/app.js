// Application shell: navigation (with back / forward controls), routing, theme, install, offline and updates.
import { h, icon, toast, download } from './core/ui.js';
import { loadProject, prefs, project, onProject } from './core/store.js';
import { loadSpec } from './core/engine.js';
import { startAutoRefresh } from './core/live.js';

const main = document.getElementById('main'), top = document.getElementById('top'), side = document.getElementById('side'), bottom = document.getElementById('bottom');
const NAV = [['#/', 'home', 'Home'], ['#/case', 'flow', 'Case pipeline'], ['#/suites', 'grid', 'Suites'], ['#/library', 'book', 'Materials library'], ['#/environments', 'flask', 'Environments'], ['#/geometry', 'cube', 'Geometry Studio'], ['#/data', 'upload', 'Data portal'], ['#/live', 'globe', 'Live data'], ['#/reference', 'doc', 'Reference'], ['#/settings', 'install', 'Install & settings']];
const routes = {
  '': async r => (await import('./pages/main.js')).homePage(r),
  suites: async r => (await import('./pages/main.js')).suitesPage(r),
  suite: async (r, a) => (await import('./pages/suite.js')).suitePage(r, Math.min(29, Math.max(1, parseInt(a) || 1))),
  case: async r => (await import('./pages/main.js')).casePage(r),
  library: async (r, a) => (await import('./pages/library.js')).libraryPage(r, a),
  environments: async (r, a) => (await import('./pages/library.js')).environmentsPage(r, a),
  geometry: async r => (await import('./pages/geometry.js')).geometryPage(r),
  data: async r => (await import('./pages/main.js')).dataPage(r),
  live: async r => (await import('./pages/main.js')).livePage(r),
  reference: async r => (await import('./pages/main.js')).referencePage(r),
  settings: async r => (await import('./pages/main.js')).settingsPage(r, ctx),
};

// ---- theme ----
function applyTheme() { const t = prefs.get('theme', 'auto'); if (t === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.dataset.theme = t; }
applyTheme();

// ---- history-aware back / forward ----
let idx = history.state?.i ?? 0, maxIdx = +(sessionStorage.getItem('dm.max') || idx);
if (history.state?.i === undefined) history.replaceState({ i: idx }, '');
const backBtn = h('button', { title: 'Back', 'aria-label': 'Back', onclick: () => history.back() }, icon('back', 18), h('span', 'Back')), fwdBtn = h('button', { title: 'Forward', 'aria-label': 'Forward', onclick: () => history.forward() }, h('span', 'Forward'), icon('fwd', 18));
const syncNav = () => { backBtn.disabled = idx <= 0; fwdBtn.disabled = idx >= maxIdx; };

// ---- install ----
let deferred = null;
addEventListener('beforeinstallprompt', e => { e.preventDefault(); deferred = e; installBtn.hidden = false; });
addEventListener('appinstalled', () => { deferred = null; installBtn.hidden = true; toast('Installed — open it from your home screen or app list', 'ok', 6000); });
const install = async () => { if (deferred) { deferred.prompt(); await deferred.userChoice; deferred = null; installBtn.hidden = true; } else location.hash = '#/settings'; };
const installBtn = h('button.btn.sm.primary', { hidden: true, onclick: install, title: 'Install on this device' }, icon('install', 15), ' Install');

// ---- service worker ----
let reg = null, swState = 'Offline support is starting…';
async function initSW() {
  if (!('serviceWorker' in navigator)) { swState = 'This browser does not support offline installation; the app still works online.'; return; }
  const hadController = !!navigator.serviceWorker.controller;
  try {
    reg = await navigator.serviceWorker.register('sw.js'); swState = 'Offline support is active on this device.';
    const offer = w => { const b = h('div#update', h('button.btn.primary.sm', { onclick: () => { w.postMessage('skipWaiting'); } }, icon('refresh', 14), ' Update available — reload')); document.body.append(b); };
    if (reg.waiting && navigator.serviceWorker.controller) offer(reg.waiting);
    reg.addEventListener('updatefound', () => { const w = reg.installing; w && w.addEventListener('statechange', () => { if (w.state === 'installed' && navigator.serviceWorker.controller) offer(w); }); });
    let reloaded = false; navigator.serviceWorker.addEventListener('controllerchange', () => { if (hadController && !reloaded) { reloaded = true; location.reload(); } });
    setInterval(() => reg.update().catch(() => { }), 30 * 60e3); addEventListener('online', () => reg.update().catch(() => { }));
    try { if ('periodicSync' in reg) { const st = await navigator.permissions.query({ name: 'periodic-background-sync' }); if (st.state === 'granted') await reg.periodicSync.register('refresh-live', { minInterval: 6 * 3600e3 }); } } catch { /* optional capability */ }
  } catch (e) { swState = 'Offline support could not start: ' + e.message; }
}
const swMsg = m => new Promise(res => { if (!navigator.serviceWorker?.controller) return res(null); const ch = new MessageChannel(); ch.port1.onmessage = e => res(e.data); navigator.serviceWorker.controller.postMessage(m, [ch.port2]); setTimeout(() => res(null), 60000); });
// Store-only ZIP of every application file, for copying to another machine or an internal server.
async function downloadBundle() {
  toast('Collecting application files…'); const list = (await (await fetch('assets.json')).json()).files, enc = new TextEncoder(), parts = [], central = []; let off = 0;
  const crcT = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })(), crc = b => { let c = ~0; for (let i = 0; i < b.length; i++) c = crcT[(c ^ b[i]) & 255] ^ (c >>> 8); return ~c >>> 0; };
  for (const f of list) { const data = new Uint8Array(await (await fetch(f)).arrayBuffer()), name = enc.encode('designer-materials/' + f), c = crc(data), hd = new DataView(new ArrayBuffer(30)); hd.setUint32(0, 0x04034b50, true); hd.setUint16(4, 20, true); hd.setUint32(14, c, true); hd.setUint32(18, data.length, true); hd.setUint32(22, data.length, true); hd.setUint16(26, name.length, true); parts.push(new Uint8Array(hd.buffer), name, data);
    const cd = new DataView(new ArrayBuffer(46)); cd.setUint32(0, 0x02014b50, true); cd.setUint16(4, 20, true); cd.setUint16(6, 20, true); cd.setUint32(16, c, true); cd.setUint32(20, data.length, true); cd.setUint32(24, data.length, true); cd.setUint16(28, name.length, true); cd.setUint32(42, off, true); central.push(new Uint8Array(cd.buffer), name); off += 30 + name.length + data.length; }
  const cs = central.reduce((a, b) => a + b.length, 0), end = new DataView(new ArrayBuffer(22)); end.setUint32(0, 0x06054b50, true); end.setUint16(8, list.length, true); end.setUint16(10, list.length, true); end.setUint32(12, cs, true); end.setUint32(16, off, true);
  download('designer-materials-offline.zip', new Blob([...parts, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' })); toast('Bundle saved. Unzip and serve the folder with any static web server (for example: python3 -m http.server).', 'ok', 9000);
}
const ctx = { applyTheme, install, canInstall: () => !!deferred, swStatus: () => swState, downloadBundle,
  checkUpdate: async () => { if (!reg) return toast(swState); await reg.update(); toast(reg.waiting || reg.installing ? 'A new version is being installed' : 'You have the latest version', 'ok'); },
  precache: async () => { toast('Downloading all suites and data for offline use…'); const r = await swMsg('precacheAll'); toast(r ? `Ready for offline use: ${r.cached} files stored` : 'Offline storage will complete after the next reload', 'ok', 6000); } };

// ---- shell ----
const net = h('span#net');
const syncNet = () => { const on = navigator.onLine !== false; net.textContent = on ? 'Online · live data' : 'Offline · cached data'; net.classList.toggle('off', !on); };
function shell(spec) {
  top.replaceChildren(h('button.iconbtn', { 'aria-label': 'Menu', id: 'menuBtn', onclick: () => side.classList.toggle('open') }, icon('menu', 18)), h('a.brand', { href: '#/' }, h('img', { src: 'icons/icon.svg', alt: '', width: 28, height: 28 }), h('span', 'Designer Materials')), h('div.navpair', { role: 'group', 'aria-label': 'History' }, backBtn, fwdBtn), h('div.grow'), net, installBtn,
    h('button.iconbtn', { 'aria-label': 'Switch light / dark', title: 'Switch light / dark', onclick: () => { const cur = document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'); prefs.set('theme', cur === 'dark' ? 'light' : 'dark'); applyTheme(); route(true); } }, icon('sun', 18)));
  side.replaceChildren(...NAV.map(([href, ic, label]) => h('a', { href, 'data-nav': href }, icon(ic, 18), label)), h('div.sep', '29 engineering suites'), ...spec.groups.map((g, k) => h('div.sgroup', { class: 'g' + (k + 1) }, h('div.sghead', h('i'), g.name), g.suites.map(i => spec.suites[i - 1]).map(s => h('a.mini', { href: `#/suite/${s.id}`, 'data-nav': `#/suite/${s.id}`, title: s.resp }, h('span.sn', String(s.id).padStart(2, '0')), h('span.st', s.name))))));
  bottom.replaceChildren(...[NAV[0], NAV[1], NAV[2], NAV[3], NAV[5]].map(([href, ic, label]) => h('a', { href, 'data-nav': href }, icon(ic, 20), label.split(' ')[0])));
  if (matchMedia('(min-width: 901px)').matches) document.getElementById('menuBtn').hidden = true;
  syncNet(); addEventListener('online', syncNet); addEventListener('offline', syncNet);
}
let token = 0;
async function route(keepScroll) {
  const [, name = '', arg] = (location.hash || '#/').split('/'), fn = routes[name] || routes[''], my = ++token, cur = name === 'suite' ? `#/suite/${parseInt(arg) || 1}` : '#/' + name;
  document.querySelectorAll('[data-nav]').forEach(a => a.classList.toggle('on', a.dataset.nav === cur)); side.classList.remove('open'); syncNav();
  const root = h('div.page'); try { await fn(root, arg && decodeURIComponent(arg)); } catch (e) { console.error(e); root.replaceChildren(h('div.card.error', h('h2', 'This page could not be opened'), h('p', navigator.onLine === false ? 'It has not been stored for offline use yet. Reconnect once, open Install & settings, and choose “Download everything for offline now”.' : 'An unexpected error occurred.'), h('pre', String(e.message || e)))); }
  if (my !== token) return; main.replaceChildren(root); if (!keepScroll) { scrollTo(0, 0); main.focus({ preventScroll: true }); }
}
addEventListener('hashchange', () => { if (history.state?.i === undefined) { idx++; maxIdx = idx; history.replaceState({ i: idx }, ''); sessionStorage.setItem('dm.max', maxIdx); } else idx = history.state.i; route(); });
addEventListener('popstate', () => { if (history.state?.i !== undefined) { idx = history.state.i; syncNav(); } });
addEventListener('keydown', e => { if (e.altKey && e.key === 'ArrowLeft') history.back(); if (e.altKey && e.key === 'ArrowRight') history.forward(); });

(async () => { try { const spec = await loadSpec(); await loadProject(); shell(spec); await route(); initSW(); startAutoRefresh(); } catch (e) { main.replaceChildren(h('div.card.error', h('h2', 'Designer Materials could not start'), h('pre', String(e.stack || e)))); } })();
