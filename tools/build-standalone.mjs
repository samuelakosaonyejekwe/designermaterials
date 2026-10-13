// Builds standalone.html: the whole application in ONE self-contained file (code, styles, specification, libraries of
// materials / environments / benchmarks). It runs from any web address or straight from disk with no other files, so it
// can be mirrored, archived or passed on as a single document. Uses esbuild only as a build-time tool.
// Usage: node tools/build-standalone.mjs
import fs from 'node:fs';
import { execFileSync } from 'node:child_process';

const tmp = fs.mkdtempSync('/tmp/dm-standalone-');
execFileSync('npx', ['--yes', 'esbuild@0.21.5', 'js/standalone.js', '--bundle', '--format=esm', '--target=es2022', '--minify-whitespace', '--minify-syntax', '--external:node:*', '--log-level=warning', `--outfile=${tmp}/bundle.js`], { stdio: 'inherit' });
const bundle = fs.readFileSync(`${tmp}/bundle.js`, 'utf8').replace(/<\/script/gi, '<\\/script');
const index = fs.readFileSync('index.html', 'utf8'), body = index.match(/<body>([\s\S]*?)<noscript>/)[1].trim(), connect = index.match(/connect-src 'self' ([^"]+)"/)[1];
const remote = (() => { try { return fs.readFileSync('.git/config', 'utf8').match(/github\.com[:/]([^/\s]+)\/([^\s.]+)/); } catch { return null; } })();
const primary = remote ? `https://${remote[1]}.github.io/${remote[2]}/` : './';
const embed = {}; for (const f of ['data/spec.json', 'assets.json', 'mirrors.json', 'CHANGELOG.md']) embed[f] = fs.readFileSync(f, 'utf8');
const icon = 'data:image/svg+xml;base64,' + fs.readFileSync('icons/icon.svg').toString('base64');
const json = JSON.stringify(embed).replace(/</g, '\\u003c');
const shim = `
globalThis.__DM_STANDALONE__ = true; globalThis.__DM_ICON__ = ${JSON.stringify(icon)}; globalThis.__DM_PRIMARY__ = ${JSON.stringify(primary)};
(() => { const files = JSON.parse(document.getElementById('dm-files').textContent), real = globalThis.fetch.bind(globalThis);
  globalThis.fetch = (input, init) => { const u = typeof input === 'string' ? input : input && input.url; if (typeof u === 'string' && !/^[a-z][a-z0-9+.-]*:/i.test(u)) { const key = u.replace(/^\\.\\//, '').split('?')[0]; if (key in files) return Promise.resolve(new Response(files[key], { status: 200, headers: { 'Content-Type': key.endsWith('.json') ? 'application/json' : 'text/plain' } })); return real(new URL(key, globalThis.__DM_PRIMARY__).href, init); } return real(input, init); }; })();`;
const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<title>Designer Materials</title>
<meta name="description" content="Designer Materials — single-file edition: the complete nanomaterials and metamaterials engineering platform in one document.">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline' 'wasm-unsafe-eval' blob: ${primary.replace(/\/[^/]*\/$/, '')}; style-src 'unsafe-inline'; img-src data: blob: ${primary.replace(/\/[^/]*\/$/, '')}; worker-src blob:; object-src 'none'; form-action 'none'; connect-src ${connect} ${primary.replace(/\/[^/]*\/$/, '')}">
<meta name="referrer" content="no-referrer">
<meta name="color-scheme" content="light dark">
<link rel="icon" href="${icon}">
<style>${fs.readFileSync('css/app.css', 'utf8').replace(/<\/style/gi, '<\\/style')}</style>
</head>
<body>
${body}
<script type="application/json" id="dm-files">${json}</script>
<script>${shim}</script>
<script type="module" id="dm-bundle">${bundle}</script>
</body>
</html>
`;
fs.writeFileSync('standalone.html', html); fs.rmSync(tmp, { recursive: true, force: true });
console.log(`standalone.html: ${(html.length / 1e6).toFixed(2)} MB`);
