// Contract and robustness tests for every engineering suite, the geometry engines and the kernels.
// Usage: node tests/run.mjs [suiteId ...]
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const only = process.argv.slice(2).map(Number).filter(Boolean);
let fail = 0, pass = 0;
const bad = (id, msg) => { fail++; console.log(`  ✗ [${id}] ${msg}`); };
const ok = () => pass++;
const isNum = v => typeof v === 'number' && Number.isFinite(v);
const KINDS = ['chemical', 'ionic', 'half-cell', 'governing', 'constitutive', 'formula'];
const GEOM = ['atoms', 'lattice', 'mesh', 'voxels', 'layout', 'rods'];

function defaults(suite, over = {}) { const o = {}; for (const i of suite.inputs) o[i.key] = i.value; return { ...o, ...over }; }
function checkResult(id, r, label) {
  const t = `${label}:`;
  if (!r || typeof r !== 'object') return bad(id, `${t} run() returned nothing`);
  if (typeof r.summary !== 'string' || r.summary.length < 40) bad(id, `${t} summary missing/too short`); else ok();
  if (!Array.isArray(r.kpis) || r.kpis.length < 4) bad(id, `${t} needs >= 4 kpis`);
  for (const k of r.kpis || []) {
    if (!k.key || !k.label) bad(id, `${t} kpi without key/label`);
    if (!(isNum(k.value) || typeof k.value === 'string')) bad(id, `${t} kpi ${k.key} value not finite: ${k.value}`); else ok();
    if (typeof k.explain !== 'string' || k.explain.length < 40) bad(id, `${t} kpi ${k.key} lacks a real explanation`);
  }
  for (const s of r.series || []) {
    if (!s.title || !s.explain || s.explain.length < 40) bad(id, `${t} series "${s.title}" lacks explanation`);
    for (const l of s.lines || []) { if (!l.x || !l.y || l.x.length !== l.y.length || !l.x.length) bad(id, `${t} series "${s.title}" line "${l.name}" malformed`); else if (!Array.from(l.y).every(v => isNum(v)) || !Array.from(l.x).every(v => isNum(v))) bad(id, `${t} series "${s.title}" line "${l.name}" has non-finite values`); else ok(); }
  }
  for (const m of r.maps || []) { if (!m.explain || m.explain.length < 40) bad(id, `${t} map "${m.title}" lacks explanation`); if (!m.data || m.data.length !== m.nx * m.ny) bad(id, `${t} map "${m.title}" size mismatch`); else if (!Array.from(m.data).every(isNum)) bad(id, `${t} map "${m.title}" non-finite`); else ok(); }
  for (const tb of r.tables || []) { if (!tb.explain || tb.explain.length < 40) bad(id, `${t} table "${tb.title}" lacks explanation`); if (!tb.columns || !tb.rows || tb.rows.some(row => row.length !== tb.columns.length)) bad(id, `${t} table "${tb.title}" malformed`); else ok(); }
  if (!Array.isArray(r.equations) || r.equations.length < 2) bad(id, `${t} needs >= 2 equations`);
  for (const e of r.equations || []) if (!KINDS.includes(e.kind) || !e.label || !e.text) bad(id, `${t} equation malformed: ${JSON.stringify(e).slice(0, 80)}`);
  if (!Array.isArray(r.recommendations) || r.recommendations.length < 3) bad(id, `${t} needs >= 3 recommendations`);
  for (const c of r.recommendations || []) if (!c.text || !c.basis || !['high', 'medium', 'low'].includes(c.priority)) bad(id, `${t} recommendation malformed`);
  if (!(r.recommendations || []).some(c => c.sustainability)) bad(id, `${t} needs at least one sustainability recommendation`);
  if (!Array.isArray(r.checks) || !r.checks.length) bad(id, `${t} needs physical/numerical checks`);
  for (const c of r.checks || []) if (!c.name || typeof c.pass !== 'boolean' || !c.detail) bad(id, `${t} check malformed`);
  if (!r.publish || !Object.keys(r.publish).length) bad(id, `${t} must publish values to the shared bus`);
  for (const [k, v] of Object.entries(r.publish || {})) if (!(isNum(v) || typeof v === 'string')) bad(id, `${t} publish ${k} invalid`);
  if (r.geometry && !GEOM.includes(r.geometry.kind)) bad(id, `${t} geometry kind invalid`);
  for (const x of r.exports || []) if (!x.name || !(typeof x.content === 'string' || x.bytes instanceof Uint8Array)) bad(id, `${t} export malformed`);
  if (!r.numerics || !r.numerics.method) bad(id, `${t} numerics.method missing`);
}

const spec = JSON.parse(fs.readFileSync(path.join(root, 'data/spec.json'), 'utf8'));
console.log('Suites');
for (let id = 1; id <= 29; id++) {
  if (only.length && !only.includes(id)) continue;
  const file = path.join(root, 'js/suites', `s${String(id).padStart(2, '0')}.js`);
  if (!fs.existsSync(file)) { bad(id, 'module missing'); continue; }
  let suite; try { suite = (await import(pathToFileURL(file))).default; } catch (e) { bad(id, 'import failed: ' + e.message); continue; }
  const before = fail;
  if (suite.id !== id) bad(id, 'id mismatch');
  if (!suite.tagline) bad(id, 'tagline missing');
  if (!Array.isArray(suite.inputs) || suite.inputs.length < 6) bad(id, 'needs >= 6 inputs');
  const keys = new Set();
  for (const i of suite.inputs || []) {
    if (keys.has(i.key)) bad(id, `duplicate input ${i.key}`); keys.add(i.key);
    if (!i.label || !i.help) bad(id, `input ${i.key} needs label and help`);
    const type = i.type || 'number';
    if (type === 'number' && !(isNum(i.value) && isNum(i.min) && isNum(i.max) && i.min <= i.value && i.value <= i.max)) bad(id, `input ${i.key} needs finite value within [min,max]`);
    if (type === 'select' && !(i.options && i.options.length)) bad(id, `input ${i.key} needs options`);
  }
  if (!Array.isArray(suite.presets) || suite.presets.length < 2) bad(id, 'needs >= 2 presets');
  const ctx = { bus: {}, spec: spec.suites[id - 1] };
  const runs = [['defaults', defaults(suite)], ...(suite.presets || []).map(p => ['preset ' + p.name, defaults(suite, p.values)])];
  // every select option
  for (const i of suite.inputs || []) if (i.type === 'select') for (const o of i.options) { const v = typeof o === 'object' ? o.value : o; if (v !== i.value) runs.push([`${i.key}=${v}`, defaults(suite, { [i.key]: v })]); }
  let maxT = 0;
  for (const [label, inp] of runs) {
    try { const t0 = performance.now(); const r = await suite.run(inp, ctx); maxT = Math.max(maxT, performance.now() - t0); checkResult(id, r, label); }
    catch (e) { bad(id, `${label}: threw ${e.message}\n${(e.stack || '').split('\n').slice(1, 3).join('\n')}`); }
  }
  if (maxT > 4000) bad(id, `too slow: ${maxT.toFixed(0)} ms`);
  // robustness: perturb every numeric input to its bounds and ±10 %
  let threw = 0, nonFinite = 0, n = 0;
  for (const i of suite.inputs || []) {
    if ((i.type || 'number') !== 'number') continue;
    for (const v of [i.min, i.max, Math.min(i.max, i.value * 1.1 + 1e-12), Math.max(i.min, i.value * 0.9)]) {
      n++; try { const r = await suite.run(defaults(suite, { [i.key]: v }), ctx); if (!r.kpis.every(k => isNum(k.value) || typeof k.value === 'string')) { nonFinite++; if (nonFinite <= 3) bad(id, `non-finite KPI at ${i.key}=${v}`); } } catch (e) { threw++; if (threw <= 3) bad(id, `threw at ${i.key}=${v}: ${e.message}`); }
    }
  }
  // determinism
  try { const a = await suite.run(defaults(suite), ctx), b = await suite.run(defaults(suite), ctx); if (JSON.stringify(a.kpis.map(k => k.value)) !== JSON.stringify(b.kpis.map(k => k.value))) bad(id, 'not deterministic (seed your RNG)'); } catch { /* reported above */ }
  console.log(`${fail === before ? '✓' : '✗'} Suite ${String(id).padStart(2)} ${spec.suites[id - 1].name}  (${runs.length} runs, ${n} perturbations, max ${maxT.toFixed(0)} ms)`);
}

// optional module self-tests
for (const f of ['tests/atomic.test.mjs', 'tests/geometry.test.mjs', 'tests/kernels.test.mjs', 'tests/data.test.mjs', 'tests/environments.test.mjs']) {
  if (only.length || !fs.existsSync(path.join(root, f))) continue;
  console.log('\n' + f);
  try { const m = await import(pathToFileURL(path.join(root, f))); const r = await m.default(); pass += r.pass; fail += r.fail; for (const l of r.log || []) console.log('  ' + l); } catch (e) { bad(f, 'threw ' + e.stack); }
}
console.log(`\n${pass} checks passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
