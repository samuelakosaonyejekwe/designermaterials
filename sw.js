// Offline support. The application shell and every suite are cached on install and refreshed in the background;
// live reference data is cached per request so the latest answers remain available without a connection.
const VERSION = 'cb4d18c1ba3e';
const SHELL = 'dm-shell-' + VERSION, LIVE = 'dm-live-v1';
const LIVE_HOSTS = ['api.openalex.org', 'api.crossref.org', 'www.ebi.ac.uk', 'zenodo.org', 'api.datacite.org', 'alexandria.icams.rub.de', 'optimade.odbx.science', 'nomad-lab.eu', 'pubchem.ncbi.nlm.nih.gov', 'en.wikipedia.org', 'api.github.com', 'api.worldbank.org', 'api.frankfurter.dev', 'api.carbonintensity.org.uk', 'api.db.nomics.world'];
const FEEDS = ['https://api.frankfurter.dev/v1/latest?base=EUR&symbols=USD', 'https://api.carbonintensity.org.uk/intensity', ...['PZINC', 'PCOPP', 'PALUM', 'PNICK', 'PTIN', 'PGOLD', 'PSILVER', 'PPLAT', 'PCOBA', 'PLITH'].map(c => `https://api.db.nomics.world/v22/series/IMF/PCPS/M.W00.${c}.USD?observations=1`)];

async function precacheAll() {
  const cache = await caches.open(SHELL), list = (await (await fetch('assets.json', { cache: 'no-store' })).json()).files; let cached = 0;
  await Promise.all(list.map(async f => { try { const r = await fetch(f, { cache: 'no-cache' }); if (r.ok) { await cache.put(new Request(f), r); cached++; } } catch { /* retried on next activation */ } }));
  return cached;
}
self.addEventListener('install', e => e.waitUntil(precacheAll()));
self.addEventListener('activate', e => e.waitUntil((async () => { for (const k of await caches.keys()) if (k.startsWith('dm-shell-') && k !== SHELL) await caches.delete(k); await self.clients.claim(); })()));
self.addEventListener('message', e => { if (e.data === 'skipWaiting') self.skipWaiting(); if (e.data === 'precacheAll') e.waitUntil(precacheAll().then(n => e.ports[0] && e.ports[0].postMessage({ cached: n }))); });
const refreshFeeds = async () => { const c = await caches.open(LIVE); try { const r = await fetch('data/live/snapshot.json', { cache: 'no-store' }); if (r.ok) await c.put(new URL('data/live/snapshot.json', location.href).pathname, r); } catch { /* offline */ } await Promise.all(FEEDS.map(async u => { try { const r = await fetch(u); if (r.ok) await c.put(u, r); } catch { /* offline */ } })); };
self.addEventListener('periodicsync', e => { if (e.tag === 'refresh-live') e.waitUntil(refreshFeeds()); });
self.addEventListener('sync', e => { if (e.tag === 'refresh-live') e.waitUntil(refreshFeeds()); });

self.addEventListener('fetch', e => {
  const req = e.request; if (req.method !== 'GET') return; const url = new URL(req.url);
  if (url.origin === location.origin && url.pathname.includes('/data/live/')) {
    // Server-refreshed data snapshot: always try the network, keep the last copy for offline use.
    e.respondWith((async () => { const cache = await caches.open(LIVE); try { const r = await fetch(req, { cache: 'no-store' }); if (r.ok) cache.put(url.pathname, r.clone()); return r; } catch (err) { return (await cache.match(url.pathname)) || Response.error(); } })());
  } else if (url.origin === location.origin) {
    // Stale-while-revalidate for the application: instant load, silent refresh, full offline fallback.
    e.respondWith((async () => { const cache = await caches.open(SHELL), key = url.pathname.endsWith('/') ? new Request(url.pathname + 'index.html') : new Request(url.pathname), hit = await cache.match(key) || await cache.match(req, { ignoreSearch: true });
      const net = fetch(req).then(r => { if (r.ok && r.type === 'basic') cache.put(key, r.clone()); return r; }).catch(() => null);
      if (hit) { e.waitUntil(net); return hit; } const r = await net; if (r) return r; if (req.mode === 'navigate') return (await cache.match('index.html')) || Response.error(); return Response.error(); })());
  } else if (LIVE_HOSTS.includes(url.hostname)) {
    // Network first for live data; fall back to the last stored answer when offline.
    e.respondWith((async () => { const cache = await caches.open(LIVE); try { const r = await fetch(req); if (r.ok) cache.put(req, r.clone()); return r; } catch (err) { const hit = await cache.match(req); if (hit) return hit; throw err; } })());
  }
});
