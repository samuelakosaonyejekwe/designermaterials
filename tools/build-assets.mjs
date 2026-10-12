// Lists every file the app needs offline (assets.json) and stamps the service worker with a content hash,
// so any change produces a new cache version. Usage: node tools/build-assets.mjs
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';

const roots = ['index.html', 'manifest.webmanifest', 'mirrors.json', 'css', 'js', 'data', 'icons'], files = [];
const walk = p => { const st = fs.statSync(p); if (st.isDirectory()) fs.readdirSync(p).sort().forEach(f => walk(path.posix.join(p, f))); else files.push(p); };
roots.filter(r => fs.existsSync(r)).forEach(walk);
const hash = crypto.createHash('sha256'); for (const f of files) { hash.update(f); hash.update(fs.readFileSync(f)); }
const build = hash.digest('hex').slice(0, 12);
const packs = []; if (fs.existsSync('vendor')) { const all = files.length; walk('vendor'); packs.push(...files.splice(all).filter(f => /\.(js|mjs|wasm)$/.test(f))); }
fs.writeFileSync('assets.json', JSON.stringify({ build, files, packs }));
const sw = fs.readFileSync('sw.js', 'utf8').replace(/const VERSION = '[^']*';/, `const VERSION = '${build}';`);
fs.writeFileSync('sw.js', sw);
console.log(`build ${build}: ${files.length} files, ${(files.reduce((a, f) => a + fs.statSync(f).size, 0) / 1e6).toFixed(2)} MB`);
