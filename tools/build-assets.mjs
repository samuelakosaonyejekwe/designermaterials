// Lists every file the app needs offline (assets.json) and stamps the service worker with a content hash,
// so any change produces a new cache version. Usage: node tools/build-assets.mjs
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';

const roots = ['index.html', 'launcher.html', 'boot.js', 'CHANGELOG.md', 'manifest.webmanifest', 'mirrors.json', 'css', 'js', 'data', 'icons'], files = [];
const walk = p => { const st = fs.statSync(p); if (st.isDirectory()) fs.readdirSync(p).sort().forEach(f => walk(path.posix.join(p, f))); else if (!p.startsWith('data/live/')) files.push(p); };
roots.filter(r => fs.existsSync(r)).forEach(walk);
const hash = crypto.createHash('sha256'); for (const f of files) { hash.update(f); hash.update(fs.readFileSync(f)); }
const build = hash.digest('hex').slice(0, 12);
const packs = []; if (fs.existsSync('vendor')) { const all = files.length; walk('vendor'); packs.push(...files.splice(all).filter(f => /\.(js|mjs|wasm)$/.test(f))); }
const version = JSON.parse(fs.readFileSync('package.json', 'utf8')).version;
fs.writeFileSync('assets.json', JSON.stringify({ build, version, files, packs }));
const sw = fs.readFileSync('sw.js', 'utf8').replace(/const VERSION = '[^']*';/, `const VERSION = '${build}';`);
fs.writeFileSync('sw.js', sw);
console.log(`build ${build}: ${files.length} files, ${(files.reduce((a, f) => a + fs.statSync(f).size, 0) / 1e6).toFixed(2)} MB`);

// Portable launcher: one small HTML file that runs the app from the content-delivery copy of the newest tagged release
// (resolved at start-up by boot.js, with the release it was built from as fallback). It works from disk or any host.
const pkg = JSON.parse(fs.readFileSync('package.json', 'utf8')), remote = (() => { try { return fs.readFileSync('.git/config', 'utf8').match(/github\.com[:/]([^/\s]+)\/([^\s.]+)/); } catch { return null; } })();
if (remote) {
  const root = `https://cdn.jsdelivr.net/gh/${remote[1]}/${remote[2]}`, body = fs.readFileSync('index.html', 'utf8').match(/<body>([\s\S]*?)<noscript>/)[1].trim();
  fs.writeFileSync('boot.js', `// Resolves the newest published release and starts the app from the content-delivery network.
(function () {
  var root = '${root}', fallback = 'v${pkg.version}', done = false;
  function start(tag) {
    if (done) return; done = true; var base = root + '@' + tag + '/';
    var b = document.createElement('base'); b.href = base; document.head.appendChild(b);
    var l = document.createElement('link'); l.rel = 'stylesheet'; l.href = base + 'css/app.css'; document.head.appendChild(l);
    var s = document.createElement('script'); s.type = 'module'; s.src = base + 'js/app.js'; document.head.appendChild(s);
  }
  var timer = setTimeout(function () { start(fallback); }, 4000);
  fetch(root + '@main/assets.json', { cache: 'no-store' }).then(function (r) { return r.json(); }).then(function (a) { clearTimeout(timer); start(/^\\d+\\.\\d+\\.\\d+$/.test(a.version) ? 'v' + a.version : fallback); }).catch(function () { clearTimeout(timer); start(fallback); });
})();
`);
  fs.writeFileSync('launcher.html', `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<title>Designer Materials</title>
<meta http-equiv="Content-Security-Policy" content="default-src https://cdn.jsdelivr.net; script-src https://cdn.jsdelivr.net 'wasm-unsafe-eval'; style-src https://cdn.jsdelivr.net; img-src https://cdn.jsdelivr.net data: blob:; font-src https://cdn.jsdelivr.net; worker-src https://cdn.jsdelivr.net blob:; object-src 'none'; form-action 'none'; connect-src ${fs.readFileSync('index.html', 'utf8').match(/connect-src 'self' ([^"]+)"/)[1]} https://cdn.jsdelivr.net https://${remote[1]}.github.io">
<meta name="referrer" content="no-referrer">
<meta name="color-scheme" content="light dark">
</head>
<body>
${body}
<script src="${root}@v${pkg.version}/boot.js"></script>
</body>
</html>
`);
  console.log('launcher.html + boot.js → ' + root + '@v' + pkg.version);
}
