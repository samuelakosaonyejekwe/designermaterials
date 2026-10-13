// Brings the Internet Archive copy up to date incrementally. For every application file it compares the SHA-1 digest
// of the newest archived capture with the published file and re-archives only those that differ, within a time budget.
// Run on every scheduled publish, so the archive copy converges to the current release and is verified by content.
// Usage: node tools/archive-catchup.mjs [minutes]
import fs from 'node:fs';
import crypto from 'node:crypto';
import { execSync } from 'node:child_process';

const budget = (Number(process.argv[2]) || 20) * 60e3, t0 = Date.now();
const remote = execSync('git remote get-url origin').toString().trim().replace(/\.git$/, ''), [, user, repo] = remote.match(/github\.com[:/]([^/]+)\/(.+)$/);
const pages = `https://${user}.github.io/${repo}/`, A = JSON.parse(fs.readFileSync('assets.json', 'utf8')), files = ['standalone.html', '', ...A.files.filter(f => f !== 'index.html'), 'assets.json', 'sw.js'];
// standalone.html comes first: it is the whole application in one document, so the archive copy is complete and
// consistent as soon as that single file has been saved.
const sleep = ms => new Promise(r => setTimeout(r, ms));
const b32 = buf => { const abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ234567'; let bits = 0, val = 0, out = ''; for (const b of buf) { val = (val << 8) | b; bits += 8; while (bits >= 5) { out += abc[(val >>> (bits - 5)) & 31]; bits -= 5; } } if (bits) out += abc[(val << (5 - bits)) & 31]; return out; };
const digest = f => fs.existsSync(f || 'index.html') ? b32(crypto.createHash('sha1').update(fs.readFileSync(f || 'index.html')).digest()) : 'missing';
const archived = async url => { try { const r = await fetch(`https://web.archive.org/cdx/search/cdx?url=${encodeURIComponent(url)}&output=json&limit=-1&fl=digest,statuscode&filter=statuscode:200`); if (!r.ok) return null; const j = await r.json(); return j.length > 1 ? j[j.length - 1][0] : ''; } catch { return null; } };

let current = 0, saved = 0, pending = [];
for (const f of files) {
  if (Date.now() - t0 > budget) { pending.push(f); continue; }
  const have = await archived(pages + f); await sleep(1200);
  if (have !== null && have === digest(f)) { current++; continue; }
  let s; try { s = (await fetch('https://web.archive.org/save/' + pages + f)).status; } catch { s = 'ERR'; }
  if (s === 200) saved++; else { pending.push(f); if (s === 429) await sleep(60000); }
  await sleep(8000);
}
console.log(`Internet Archive copy: ${current} files already current, ${saved} re-archived, ${pending.length} still pending (picked up by the next run).`);
