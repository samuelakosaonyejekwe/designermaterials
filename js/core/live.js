// Live data layer. Every request is made directly from the user's browser to public scientific,
// engineering and economic services; responses are cached locally (stale-while-revalidate) so the
// last known data stays available offline and is refreshed whenever a connection exists.
import { kv, prefs } from './store.js';

const J = async (url, opt = {}) => { const c = new AbortController(), t = setTimeout(() => c.abort(), opt.timeout || 15000); try { const r = await fetch(url, { signal: c.signal, keepalive: true, headers: opt.headers || { Accept: 'application/json' } }); if (!r.ok) throw new Error(`HTTP ${r.status}`); return opt.text ? r.text() : r.json(); } finally { clearTimeout(t); } };
const enc = encodeURIComponent, yearsAgo = n => `${new Date().getFullYear() - n}-01-01`;
const dbn = code => async () => { const d = await J(`https://api.db.nomics.world/v22/series/IMF/PCPS/M.W00.${code}.USD?observations=1`), s = d.series.docs[0]; let i = s.value.length - 1; while (i > 0 && typeof s.value[i] !== 'number') i--; return { value: s.value[i], period: s.period[i], history: { x: s.period.slice(-60), y: s.value.slice(-60) } }; };

// Scalar feeds that can fill suite inputs (keys referenced by `live:` on inputs)
export const FEEDS = {
  'fx.USD_per_EUR': { label: 'USD per EUR', unit: 'USD/EUR', source: 'European Central Bank reference rates (Frankfurter)', url: 'https://www.frankfurter.dev', get: async () => { const d = await J('https://api.frankfurter.dev/v1/latest?base=EUR&symbols=USD'); return { value: d.rates.USD, period: d.date }; } },
  'grid.co2_g_kWh': { label: 'Grid carbon intensity (Great Britain)', unit: 'g CO₂/kWh', source: 'National Energy System Operator Carbon Intensity API', url: 'https://carbonintensity.org.uk', get: async () => { const d = await J('https://api.carbonintensity.org.uk/intensity'), x = d.data[0]; return { value: x.intensity.actual ?? x.intensity.forecast, period: x.from, note: x.intensity.index }; } },
  'price.zinc_usd_t': { label: 'Zinc price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PZINC') },
  'price.copper_usd_t': { label: 'Copper price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PCOPP') },
  'price.aluminium_usd_t': { label: 'Aluminium price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PALUM') },
  'price.nickel_usd_t': { label: 'Nickel price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PNICK') },
  'price.tin_usd_t': { label: 'Tin price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PTIN') },
  'price.gold_usd_oz': { label: 'Gold price', unit: 'USD/troy oz', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PGOLD') },
  'price.silver_usd_oz': { label: 'Silver price', unit: 'USD/troy oz', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PSILVER') },
  'price.platinum_usd_oz': { label: 'Platinum price', unit: 'USD/troy oz', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PPLAT') },
  'price.cobalt_usd_t': { label: 'Cobalt price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PCOBA') },
  'price.lithium_usd_t': { label: 'Lithium price', unit: 'USD/t', source: 'IMF Primary Commodity Prices (DBnomics)', url: 'https://db.nomics.world/IMF/PCPS', get: dbn('PLITH') },
};

const work = w => ({ title: w.title || w.display_name, year: w.publication_year, date: w.publication_date, venue: w.primary_location?.source?.display_name || '', doi: w.doi ? w.doi.replace('https://doi.org/', '') : null, url: w.doi || w.id, cited: w.cited_by_count, oa: !!w.open_access?.is_oa });
// Query connectors (search-type sources)
export const CONNECTORS = {
  openalex: { name: 'OpenAlex', org: 'OurResearch', kind: 'Literature', home: 'https://openalex.org', about: 'Open index of scholarly works; used for the latest peer-reviewed practice on each topic.',
    search: async q => (await J(`https://api.openalex.org/works?search=${enc(q)}&filter=from_publication_date:${yearsAgo(3)},type:article&sort=relevance_score:desc&per-page=8&select=id,title,publication_year,publication_date,doi,cited_by_count,primary_location,open_access`)).results.map(work) },
  crossref: { name: 'Crossref', org: 'Crossref', kind: 'Literature', home: 'https://www.crossref.org', about: 'DOI registration agency metadata for journals, standards and reports.',
    search: async q => (await J(`https://api.crossref.org/works?query=${enc(q)}&rows=8&filter=from-pub-date:${yearsAgo(3)}&select=title,DOI,issued,container-title,is-referenced-by-count,type`)).message.items.map(i => ({ title: (i.title || ['(untitled)'])[0], year: i.issued?.['date-parts']?.[0]?.[0], venue: (i['container-title'] || [''])[0], doi: i.DOI, url: 'https://doi.org/' + i.DOI, cited: i['is-referenced-by-count'], type: i.type })) },
  europepmc: { name: 'Europe PMC', org: 'EMBL-EBI', kind: 'Literature (life sciences, toxicology)', home: 'https://europepmc.org', about: 'Biomedical and nanotoxicology literature; used by safety and biomedical assessments.',
    search: async q => ((await J(`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=${enc(q)}%20sort_date:y&format=json&pageSize=8`)).resultList?.result || []).map(i => ({ title: i.title, year: +i.pubYear, venue: i.journalTitle || '', doi: i.doi || null, url: i.doi ? 'https://doi.org/' + i.doi : `https://europepmc.org/article/${i.source}/${i.id}`, cited: i.citedByCount })) },
  zenodo: { name: 'Zenodo', org: 'CERN', kind: 'Datasets and software', home: 'https://zenodo.org', about: 'Open research datasets, simulation inputs and software records.',
    search: async q => ((await J(`https://zenodo.org/api/records?q=${enc(q)}&size=8&sort=mostrecent&type=dataset`)).hits?.hits || []).map(i => ({ title: i.metadata?.title || i.title, year: +(i.metadata?.publication_date || '').slice(0, 4), venue: 'Zenodo dataset', doi: i.doi || null, url: i.links?.self_html || i.links?.html || ('https://doi.org/' + i.doi) })) },
  datacite: { name: 'DataCite', org: 'DataCite', kind: 'Datasets', home: 'https://datacite.org', about: 'DOI metadata for research datasets across repositories.',
    search: async q => ((await J(`https://api.datacite.org/dois?query=${enc(q)}&resource-type-id=dataset&page[size]=8&sort=-created`)).data || []).map(i => ({ title: i.attributes.titles?.[0]?.title, year: i.attributes.publicationYear, venue: i.attributes.publisher, doi: i.attributes.doi, url: 'https://doi.org/' + i.attributes.doi })) },
  alexandria: { name: 'Alexandria (OPTIMADE)', org: 'Ruhr-Universität Bochum', kind: 'Computed crystal structures (DFT, PBE)', home: 'https://alexandria.icams.rub.de', about: 'Millions of DFT-relaxed crystal structures with band gap and formation energy, served through the OPTIMADE standard.',
    structures: async els => ((await J(`https://alexandria.icams.rub.de/pbe/v1/structures?filter=${enc(`elements HAS ALL ${els.map(e => `"${e}"`).join(',')} AND nelements=${els.length}`)}&page_limit=10`)).data || []).map(optimade('alexandria')) },
  odbx: { name: 'odbx (OPTIMADE)', org: 'Open Database of Xtals', kind: 'Computed crystal structures', home: 'https://odbx.science', about: 'Open crystal-structure database exposed through the OPTIMADE standard.',
    structures: async els => ((await J(`https://optimade.odbx.science/v1/structures?filter=${enc(`elements HAS ALL ${els.map(e => `"${e}"`).join(',')} AND nelements=${els.length}`)}&page_limit=10`)).data || []).map(optimade('odbx')) },
  nomad: { name: 'NOMAD', org: 'FAIRmat / Max Planck', kind: 'Materials simulation archive', home: 'https://nomad-lab.eu', about: 'Archive of electronic-structure and experiment records; used to gauge how much reference data exists for a composition.',
    count: async reduced => { const d = await J(`https://nomad-lab.eu/prod/v1/api/v1/entries?json_query=${enc(JSON.stringify({ 'results.material.chemical_formula_reduced': reduced }))}&page_size=1`); return { total: d.pagination?.total ?? 0, url: `https://nomad-lab.eu/prod/v1/gui/search/entries?formula=${enc(reduced)}` }; } },
  pubchem: { name: 'PubChem', org: 'US National Library of Medicine', kind: 'Chemical identity, properties and hazards', home: 'https://pubchem.ncbi.nlm.nih.gov', about: 'Compound identity, computed properties and GHS hazard classification.',
    compound: async name => { const p = (await J(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${enc(name)}/property/MolecularFormula,MolecularWeight,IUPACName/JSON`)).PropertyTable.Properties[0]; let ghs = []; try { const v = await J(`https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/${p.CID}/JSON?heading=GHS+Classification`), out = new Set(); const walk = n => { if (!n) return; if (Array.isArray(n)) return n.forEach(walk); if (typeof n === 'object') { if (n.Name === 'GHS Hazard Statements' || n.Name === 'Signal') (n.Value?.StringWithMarkup || []).forEach(s => out.add(s.String)); Object.values(n).forEach(walk); } }; walk(v.Record); ghs = [...out].slice(0, 14); } catch { /* no GHS record */ } return { cid: p.CID, formula: p.MolecularFormula, mw: +p.MolecularWeight, iupac: p.IUPACName, ghs, url: `https://pubchem.ncbi.nlm.nih.gov/compound/${p.CID}` }; } },
  wikipedia: { name: 'Wikipedia', org: 'Wikimedia', kind: 'Background summaries', home: 'https://en.wikipedia.org', about: 'Plain-language background for concepts; not used as a data source for calculations.',
    summary: async title => { const d = await J(`https://en.wikipedia.org/api/rest_v1/page/summary/${enc(title)}`); return { title: d.title, extract: d.extract, url: d.content_urls?.desktop?.page }; } },
  github: { name: 'GitHub releases', org: 'GitHub', kind: 'Scientific software versions', home: 'https://github.com', about: 'Current releases of open scientific codes this platform exchanges files with.',
    release: async repo => { const d = await J(`https://api.github.com/repos/${repo}/releases/latest`); return { repo, tag: d.tag_name, date: (d.published_at || '').slice(0, 10), url: d.html_url }; } },
  worldbank: { name: 'World Bank Open Data', org: 'World Bank', kind: 'Economic indicators', home: 'https://data.worldbank.org', about: 'Macro-economic indicators used for discount-rate and cost-escalation context.',
    indicator: async (code, country = 'WLD') => { const d = await J(`https://api.worldbank.org/v2/country/${country}/indicator/${code}?format=json&per_page=8`), rows = (d[1] || []).filter(r => r.value !== null); return rows.length ? { value: rows[0].value, period: rows[0].date, name: rows[0].indicator.value } : null; } },
};
function optimade(src) { return d => { const a = d.attributes; return { id: d.id, source: src, formula: a.chemical_formula_reduced, nsites: a.nsites, spacegroup: a._alexandria_space_group ?? null, bandgap: a._alexandria_band_gap ?? null, eform: a._alexandria_formation_energy_per_atom ?? null, hull: a._alexandria_hull_distance ?? null, geometry: { kind: 'atoms', title: `${a.chemical_formula_reduced} (${src} ${d.id})`, cell: a.lattice_vectors, atoms: (a.cartesian_site_positions || []).map((p, i) => ({ el: a.species_at_sites[i], x: p[0], y: p[1], z: p[2] })) } }; }; }

// Reference services without a browser-accessible API: opened as deep links instead of being queried in-app.
export const PORTALS = [
  { name: 'Materials Project', kind: 'Computed materials properties', url: q => `https://next-gen.materialsproject.org/materials?formula=${enc(q)}` },
  { name: 'Crystallography Open Database', kind: 'Experimental crystal structures (CIF)', url: q => `https://www.crystallography.net/cod/result.php?formula=${enc(q)}` },
  { name: 'AFLOW', kind: 'High-throughput computed materials', url: q => `https://aflow.org/search/?search=${enc(q)}` },
  { name: 'OQMD', kind: 'Formation energies and phase stability', url: q => `https://oqmd.org/materials/composition/${enc(q)}` },
  { name: 'NIST Chemistry WebBook', kind: 'Thermochemical and spectral data', url: q => `https://webbook.nist.gov/cgi/cbook.cgi?Formula=${enc(q)}&NoIon=on&Units=SI` },
  { name: 'JARVIS (NIST)', kind: 'DFT, force-field and ML materials data', url: q => `https://jarvis.nist.gov/jarvisdft/?q=${enc(q)}` },
  { name: 'Materials Cloud', kind: 'Curated computational datasets', url: q => `https://www.materialscloud.org/discover/mc3d/dashboard/ptable?q=${enc(q)}` },
  { name: 'RRUFF', kind: 'Raman, XRD and chemistry of minerals', url: q => `https://rruff.info/${enc(q)}` },
  { name: 'ECHA substance information', kind: 'EU chemical registration and classification', url: q => `https://chem.echa.europa.eu/substance-search?searchText=${enc(q)}` },
  { name: 'NIOSH nanotechnology', kind: 'Occupational exposure guidance', url: () => 'https://www.cdc.gov/niosh/nanotech/' },
  { name: 'ISO/TC 229 Nanotechnologies', kind: 'Standards catalogue', url: () => 'https://www.iso.org/committee/381983/x/catalogue/' },
  { name: 'refractiveindex.info', kind: 'Optical constants', url: q => `https://refractiveindex.info/?search=${enc(q)}` },
  { name: 'arXiv', kind: 'Preprints', url: q => `https://arxiv.org/search/?query=${enc(q)}&searchtype=all&order=-announced_date_first` },
  { name: 'USGS Mineral Commodity Summaries', kind: 'Raw-material supply and criticality', url: () => 'https://www.usgs.gov/centers/national-minerals-information-center/mineral-commodity-summaries' },
];
export const TOOL_REPOS = ['materialsproject/pymatgen', 'lammps/lammps', 'QEF/q-e', 'cp2k/cp2k', 'OpenFOAM/OpenFOAM-dev', 'FEniCS/dolfinx', 'NanoComp/meep', 'KLayout/klayout', 'openscad/openscad', 'ovito-org/ovito', 'Materials-Consortia/OPTIMADE'];

// ---- cache ----
const TTL = { feed: 6 * 3600e3, search: 12 * 3600e3, struct: 7 * 86400e3 };
export const online = () => typeof navigator === 'undefined' || navigator.onLine !== false;
async function cached(key, ttl, fn, { force = false } = {}) {
  const hit = await kv.get('live:' + key), fresh = hit && Date.now() - hit.at < ttl;
  if (hit && (fresh || !online()) && !force) return { ...hit, cached: true, stale: !fresh };
  try { const data = await fn(); const rec = { data, at: Date.now() }; await kv.set('live:' + key, rec); return { ...rec, cached: false, stale: false }; }
  catch (e) { if (hit) return { ...hit, cached: true, stale: true, error: e.message }; throw e; }
}
export const feed = (key, opt) => cached('feed:' + key, TTL.feed, FEEDS[key].get, opt);
export const literature = (q, src = 'openalex', opt) => cached(`lit:${src}:${q}`, TTL.search, () => CONNECTORS[src].search(q), opt);
export const structures = (els, src = 'alexandria', opt) => cached(`str:${src}:${[...els].sort().join('-')}`, TTL.struct, () => CONNECTORS[src].structures(els), opt);
export const compound = (name, opt) => cached('cmp:' + name.toLowerCase(), TTL.struct, () => CONNECTORS.pubchem.compound(name), opt);
export const nomadCount = (f, opt) => cached('nomad:' + f, TTL.struct, () => CONNECTORS.nomad.count(f), opt);
export const summary = (t, opt) => cached('wiki:' + t, TTL.struct, () => CONNECTORS.wikipedia.summary(t), opt);
export const release = (r, opt) => cached('rel:' + r, TTL.search, () => CONNECTORS.github.release(r), opt);
export const indicator = (c, opt) => cached('wb:' + c, TTL.struct, () => CONNECTORS.worldbank.indicator(c), opt);
export async function refreshFeeds(force = false) {
  const out = {}; await Promise.all(Object.keys(FEEDS).map(async k => { try { out[k] = await feed(k, { force }); } catch (e) { out[k] = { error: e.message }; } }));
  prefs.set('lastFeedRefresh', Date.now()); return out;
}
export const age = at => { const s = (Date.now() - at) / 1000; return s < 90 ? 'just now' : s < 3600 ? `${Math.round(s / 60)} min ago` : s < 86400 ? `${Math.round(s / 3600)} h ago` : `${Math.round(s / 86400)} d ago`; };
// Keep data fresh while the app is open and whenever connectivity returns.
export function startAutoRefresh(onUpdate) {
  const go = async () => { if (!online()) return; if (Date.now() - prefs.get('lastFeedRefresh', 0) > 30 * 60e3) { await refreshFeeds(); onUpdate && onUpdate(); } };
  go(); setInterval(go, 10 * 60e3); addEventListener('online', go);
  // Refresh when the app comes to the foreground, and once more as it is sent to the background (requests are kept alive),
  // so the stored copy is as recent as the last moment the app was in use on platforms that suspend closed web apps.
  const leave = () => { if (online() && Date.now() - prefs.get('lastFeedRefresh', 0) > 5 * 60e3) refreshFeeds(true).then(() => onUpdate && onUpdate()).catch(() => { }); };
  document.addEventListener('visibilitychange', () => document.visibilityState === 'visible' ? go() : leave()); addEventListener('pagehide', leave);
}
// Search phrases per suite for "current practice" literature look-ups
export const SUITE_TOPICS = { 1: 'materials selection multi-objective design requirements', 2: 'composition design dopant phase stability', 3: 'crystal structure defects first-principles', 4: 'nanoparticle morphology control shape', 5: 'metamaterial unit cell design lattice', 6: 'solubility speciation thermodynamics precipitation', 7: 'nucleation growth kinetics nanocrystal', 8: 'colloidal stability surface functionalization', 9: 'porous media reactor flow mixing nanoparticle synthesis', 10: 'thermal transport nanoscale heat transfer', 11: 'electrodeposition charge transport nanostructure', 12: 'metasurface plasmonic optical response', 13: 'magnetic nanoparticles spintronic', 14: 'mechanical metamaterial auxetic nanomechanics', 15: 'acoustic metamaterial phononic bandgap', 16: 'phase-field microstructure grain growth', 17: 'molecular dynamics machine-learning interatomic potential', 18: 'homogenization multiscale effective properties', 19: 'architected material component large-scale analysis', 20: 'tunable reconfigurable metamaterial control', 21: 'synthesis process optimization scalable fabrication', 22: 'degradation fatigue lifetime reliability', 23: 'characterization metrology measurement uncertainty', 24: 'verification validation uncertainty quantification simulation', 25: 'inverse design machine learning topology optimization', 26: 'scale-up continuous manufacturing quality control', 27: 'nanomaterial safety exposure risk assessment', 28: 'techno-economic life cycle assessment sustainability', 29: 'reverse engineering reconstruction tomography inverse problem' };
