// Collects reference data on a server schedule so every device (including phones that suspend closed web apps)
// receives a fresh same-origin snapshot in one request, and so sources that refuse direct browser requests are covered.
// Usage: node tools/fetch-live.mjs   → writes data/live/snapshot.json (never throws: failed sources keep their last value)
import fs from 'node:fs';
import { pathToFileURL } from 'node:url';
import path from 'node:path';

const out = 'data/live/snapshot.json', prev = fs.existsSync(out) ? JSON.parse(fs.readFileSync(out, 'utf8')) : {};
const sleep = ms => new Promise(r => setTimeout(r, ms));
const get = async (url, { text = false, timeout = 25000 } = {}) => { const c = new AbortController(), t = setTimeout(() => c.abort(), timeout); try { const r = await fetch(url, { signal: c.signal, headers: { 'User-Agent': 'designer-materials-data-refresh/1.0', Accept: text ? '*/*' : 'application/json' } }); if (!r.ok) throw new Error('HTTP ' + r.status); return text ? r.text() : r.json(); } finally { clearTimeout(t); } };
const snap = { at: new Date().toISOString(), feeds: { ...(prev.feeds || {}) }, arxiv: { ...(prev.arxiv || {}) }, cod: { ...(prev.cod || {}) }, sources: {} };
const mark = (name, ok, n) => { snap.sources[name] = { ok, records: n, at: new Date().toISOString() }; console.log(`${ok ? 'ok ' : 'FAIL'} ${name}${n !== undefined ? ' (' + n + ')' : ''}`); };

// 1. scalar feeds
const pcps = { 'price.zinc_usd_t': 'PZINC', 'price.copper_usd_t': 'PCOPP', 'price.aluminium_usd_t': 'PALUM', 'price.nickel_usd_t': 'PNICK', 'price.tin_usd_t': 'PTIN', 'price.gold_usd_oz': 'PGOLD', 'price.silver_usd_oz': 'PSILVER', 'price.platinum_usd_oz': 'PPLAT', 'price.cobalt_usd_t': 'PCOBA', 'price.lithium_usd_t': 'PLITH' };
let nf = 0;
for (const [k, code] of Object.entries(pcps)) { try { const s = (await get(`https://api.db.nomics.world/v22/series/IMF/PCPS/M.W00.${code}.USD?observations=1`)).series.docs[0]; let i = s.value.length - 1; while (i > 0 && typeof s.value[i] !== 'number') i--; snap.feeds[k] = { value: s.value[i], period: s.period[i], history: { x: s.period.slice(-60), y: s.value.slice(-60) } }; nf++; } catch { /* keep previous */ } }
try { const d = await get('https://api.frankfurter.dev/v1/latest?base=EUR&symbols=USD'); snap.feeds['fx.USD_per_EUR'] = { value: d.rates.USD, period: d.date }; nf++; } catch { /* keep */ }
try { const x = (await get('https://api.carbonintensity.org.uk/intensity')).data[0]; snap.feeds['grid.co2_g_kWh'] = { value: x.intensity.actual ?? x.intensity.forecast, period: x.from, note: x.intensity.index }; nf++; } catch { /* keep */ }
mark('feeds', nf > 0, nf);

// 2. arXiv preprints per suite topic (arXiv does not allow direct browser requests)
const topics = { 1: 'materials selection multi-objective design', 2: 'dopant phase stability composition', 3: 'crystal structure defects first-principles', 4: 'nanoparticle morphology control', 5: 'metamaterial unit cell lattice design', 6: 'solubility speciation precipitation thermodynamics', 7: 'nucleation growth kinetics nanocrystal', 8: 'colloidal stability surface functionalization', 9: 'porous media reactive transport', 10: 'nanoscale heat transfer thermal metamaterial', 11: 'electrodeposition charge transport nanostructure', 12: 'metasurface plasmonic negative index', 13: 'magnetic nanoparticles spintronics', 14: 'mechanical metamaterial auxetic', 15: 'acoustic metamaterial phononic bandgap', 16: 'phase-field microstructure grain growth', 17: 'machine-learning interatomic potential molecular dynamics', 18: 'homogenization multiscale effective properties', 19: 'architected material structural analysis', 20: 'tunable reconfigurable metamaterial', 21: 'hydrothermal synthesis nanorod process', 22: 'fatigue degradation lifetime reliability materials', 23: 'materials characterization metrology uncertainty', 24: 'verification validation uncertainty quantification', 25: 'inverse design topology optimization machine learning materials', 26: 'scale-up continuous manufacturing nanomaterials', 27: 'nanomaterial safety exposure risk', 28: 'techno-economic life cycle assessment materials', 29: 'inverse problem reconstruction tomography materials' };
let na = 0;
for (const [id, q] of Object.entries(topics)) {
  try { const xml = await get(`https://export.arxiv.org/api/query?search_query=${encodeURIComponent('all:' + q.split(' ').join(' AND all:'))}&sortBy=submittedDate&sortOrder=descending&max_results=6`, { text: true });
    const items = [...xml.matchAll(/<entry>([\s\S]*?)<\/entry>/g)].map(m => { const g = t => (m[1].match(new RegExp(`<${t}[^>]*>([\\s\\S]*?)</${t}>`)) || [])[1]?.replace(/\s+/g, ' ').trim(); return { title: g('title'), url: g('id'), year: +(g('published') || '').slice(0, 4), date: (g('published') || '').slice(0, 10), venue: 'arXiv preprint' }; }).filter(x => x.title);
    if (items.length) { snap.arxiv[id] = items; na++; } } catch { /* keep */ }
  await sleep(3200);
}
mark('arxiv', na > 0, na);

// 3. Crystallography Open Database: experimental structures for the library's compositions
try {
  const { MATERIALS } = await import(pathToFileURL(path.resolve('data/materials.js'))), chem = await import(pathToFileURL(path.resolve('js/core/chem.js')));
  const hill = f => { const a = chem.parseFormula(f).atoms, els = Object.keys(a), hasC = els.includes('C'); els.sort((x, y) => hasC ? ((x === 'C' ? -2 : x === 'H' ? -1 : 0) - (y === 'C' ? -2 : y === 'H' ? -1 : 0) || x.localeCompare(y)) : x.localeCompare(y)); return els.map(e => e + (a[e] === 1 ? '' : a[e])).join(' '); };
  const formulas = [...new Set(MATERIALS.filter(m => m.formula).map(m => { try { const a = chem.parseFormula(m.formula).atoms; return Object.keys(a).length <= 4 && Object.values(a).every(Number.isInteger) ? m.formula : null; } catch { return null; } }).filter(Boolean))];
  let nc = 0;
  for (const f of formulas) {
    try { const rows = await get(`https://www.crystallography.net/cod/result?formula=${encodeURIComponent(hill(f))}&format=json`);
      snap.cod[f] = { total: rows.length, entries: rows.slice(0, 6).map(r => ({ id: r.file, sg: r.sg, a: +r.a, b: +r.b, c: +r.c, alpha: +r.alpha, beta: +r.beta, gamma: +r.gamma, vol: +r.vol, year: r.year, mineral: r.mineral || null, title: (r.title || '').slice(0, 140) })) }; nc++; } catch { /* keep */ }
    await sleep(600);
  }
  mark('cod', nc > 0, nc);
} catch (e) { mark('cod', false, 0); console.log(e.message); }

fs.mkdirSync('data/live', { recursive: true }); fs.writeFileSync(out, JSON.stringify(snap));
console.log(`snapshot written: ${(fs.statSync(out).size / 1024).toFixed(0)} kB`);
