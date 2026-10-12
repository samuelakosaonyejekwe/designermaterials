// Runs the reference benchmarks. Usage: node tests/benchmarks.mjs [suiteId ...]
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..'), only = process.argv.slice(2).map(Number).filter(Boolean);
const imp = f => import(pathToFileURL(path.join(root, f)));
const f4 = v => typeof v === 'number' ? (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(4) : String(+v.toPrecision(6))) : String(v);
let pass = 0, fail = 0, missing = []; const ids = new Set();
for (let id = 1; id <= 29; id++) {
  if (only.length && !only.includes(id)) continue;
  const file = `data/benchmarks/s${String(id).padStart(2, '0')}.js`; if (!fs.existsSync(path.join(root, file))) { missing.push(id); continue; }
  const list = (await imp(file)).default, suite = (await imp(`js/suites/s${String(id).padStart(2, '0')}.js`)).default, def = Object.fromEntries(suite.inputs.map(i => [i.key, i.value]));
  console.log(`Suite ${id} — ${list.length} benchmarks`);
  for (const b of list) {
    let got = NaN, err = '';
    if (ids.has(b.id)) err = 'duplicate id'; ids.add(b.id);
    for (const k of ['id', 'title', 'inputs', 'expected', 'kind', 'reference']) if (b[k] === undefined) err = `missing field ${k}`;
    if (b.suite !== id) err = 'suite mismatch'; if (b.tol === undefined && b.abs === undefined) err = 'needs tol or abs';
    for (const k in b.inputs || {}) if (!(k in def)) err = `unknown input "${k}"`;
    try { const r = await suite.run({ ...def, ...b.inputs }, { bus: {}, benchmark: true }); got = b.pub ? r.publish?.[b.pub] : r.kpis.find(k => k.key === b.kpi)?.value; } catch (e) { err = 'threw: ' + e.message; }
    const abs = Math.abs(got - b.expected), dev = b.abs !== undefined ? abs : abs / Math.abs(b.expected || 1), ok = !err && Number.isFinite(got) && (b.abs !== undefined ? abs <= b.abs : dev <= b.tol);
    ok ? pass++ : fail++;
    console.log(`  ${ok ? 'PASS' : 'FAIL'}  ${b.id}: expected ${f4(b.expected)}, got ${f4(got)}, ${b.abs !== undefined ? 'abs dev ' + f4(dev) + ' (≤ ' + b.abs + ')' : 'dev ' + (dev * 100).toFixed(3) + ' % (≤ ' + (b.tol * 100) + ' %)'} ${err}`);
  }
}
console.log(`\n${pass} benchmarks passed, ${fail} failed${missing.length ? `; no benchmark file for suites ${missing.join(', ')}` : ''}`);
process.exit(fail ? 1 : 0);
