// Publishes independent copies of a tagged release so the app stays reachable if the primary host is down.
// Usage: node tools/publish-mirrors.mjs <tag> [--skip-archive]
//   1. cdn.jsdelivr.net/gh/<user>/<repo>@<tag>/ — CDN copy stored permanently once each file has been requested
//      (this script requests every file, including the on-demand kernels). launcher.html runs the app from it.
//   2. web.archive.org — Internet Archive snapshot of every application file.
//   3. archive.softwareheritage.org — permanent archive of the source repository.
// No account is needed for any of them.
import fs from 'node:fs';
import { execSync } from 'node:child_process';

const tag = process.argv[2]; if (!tag) { console.error('usage: node tools/publish-mirrors.mjs <tag>'); process.exit(1); }
const remote = execSync('git remote get-url origin').toString().trim().replace(/\.git$/, ''), [, user, repo] = remote.match(/github\.com[:/]([^/]+)\/(.+)$/);
const A = JSON.parse(fs.readFileSync('assets.json', 'utf8')), files = A.files.concat(A.packs || [], ['assets.json', 'sw.js']);
const pages = `https://${user}.github.io/${repo}/`, cdn = `https://cdn.jsdelivr.net/gh/${user}/${repo}@${tag}/`;
const sleep = ms => new Promise(r => setTimeout(r, ms));
const get = async (url, opt) => { try { const r = await fetch(url, opt); return r.status; } catch (e) { return 'ERR ' + e.message; } };

let ok = 0; for (const f of files) { const s = await get(cdn + f); if (s === 200) ok++; else console.log('  cdn', f, s); }
console.log(`CDN copy stored: ${ok}/${files.length} files at ${cdn}`);

// Tell the CDN that the release pointer changed, so the portable launcher resolves the new version at once.
console.log('Launcher pointer refreshed:', await get(`https://purge.jsdelivr.net/gh/${user}/${repo}@main/assets.json`));

console.log('Software Heritage:', await get(`https://archive.softwareheritage.org/api/1/origin/save/git/url/${remote}/`, { method: 'POST' }));

if (!process.argv.includes('--skip-archive')) {
  let saved = 0;
  for (const f of ['', ...files.filter(f => f !== 'index.html')]) {
    let s; for (let attempt = 0; attempt < 4; attempt++) { s = await get('https://web.archive.org/save/' + pages + f, { redirect: 'follow' }); if (s === 200) break; await sleep(40000); }
    if (s === 200) saved++; else console.log('  archive', f || '(root)', s); await sleep(7000);
  }
  console.log(`Internet Archive: ${saved}/${files.length} files saved — https://web.archive.org/web/${pages}`);
}
