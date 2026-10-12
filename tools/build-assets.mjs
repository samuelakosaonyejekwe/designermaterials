// Lists every file the app needs offline (assets.json) and stamps the service worker with a content hash,
// so any change produces a new cache version. Usage: node tools/build-assets.mjs
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';

const roots = ['index.html', 'launcher.html', 'manifest.webmanifest', 'mirrors.json', 'css', 'js', 'data', 'icons'], files = [];
const walk = p => { const st = fs.statSync(p); if (st.isDirectory()) fs.readdirSync(p).sort().forEach(f => walk(path.posix.join(p, f))); else files.push(p); };
roots.filter(r => fs.existsSync(r)).forEach(walk);
const hash = crypto.createHash('sha256'); for (const f of files) { hash.update(f); hash.update(fs.readFileSync(f)); }
const build = hash.digest('hex').slice(0, 12);
const packs = []; if (fs.existsSync('vendor')) { const all = files.length; walk('vendor'); packs.push(...files.splice(all).filter(f => /\.(js|mjs|wasm)$/.test(f))); }
fs.writeFileSync('assets.json', JSON.stringify({ build, files, packs }));
const sw = fs.readFileSync('sw.js', 'utf8').replace(/const VERSION = '[^']*';/, `const VERSION = '${build}';`);
fs.writeFileSync('sw.js', sw);
console.log(`build ${build}: ${files.length} files, ${(files.reduce((a, f) => a + fs.statSync(f).size, 0) / 1e6).toFixed(2)} MB`);

// Portable launcher: a single HTML file that runs the app from the permanently cached CDN copy of the tagged release.
// It works when opened from disk or from any other host, so the app stays reachable if the primary host is down.
const pkg = JSON.parse(fs.readFileSync('package.json', 'utf8')), remote = (() => { try { return fs.readFileSync('.git/config', 'utf8').match(/github\.com[:/]([^/\s]+)\/([^\s.]+)/); } catch { return null; } })();
if (remote) {
  const cdn = `https://cdn.jsdelivr.net/gh/${remote[1]}/${remote[2]}@v${pkg.version}/`;
  const html = fs.readFileSync('index.html', 'utf8')
    .replace('<meta charset="utf-8">', `<meta charset="utf-8">\n<base href="${cdn}">`)
    .replace(/default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; style-src 'self'; img-src 'self' data: blob:; font-src 'self'; worker-src 'self'; manifest-src 'self'; object-src 'none'; base-uri 'self'; form-action 'none'; connect-src 'self'/, `default-src https://cdn.jsdelivr.net; script-src https://cdn.jsdelivr.net 'wasm-unsafe-eval'; style-src https://cdn.jsdelivr.net; img-src https://cdn.jsdelivr.net data: blob:; font-src https://cdn.jsdelivr.net; worker-src https://cdn.jsdelivr.net blob:; object-src 'none'; form-action 'none'; connect-src https://cdn.jsdelivr.net`)
    .replace('<link rel="manifest" href="manifest.webmanifest">\n', '');
  if (!html.includes('<base href=') || !html.includes("default-src https://cdn.jsdelivr.net")) throw new Error('launcher template out of date');
  fs.writeFileSync('launcher.html', html); console.log('launcher.html → ' + cdn);
}
