import { readFileSync, writeFileSync, readdirSync, rmSync, cpSync } from 'node:fs';
import { resolve, join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { patchNavigation } from './patch-donnieplus-navigation.mjs';
import { patchAccountOnly } from './patch-donnieplus-account.mjs';
import { patchLogin } from './patch-donnieplus-login.mjs';

const upstream = resolve(process.argv[2] || '');
if (!process.argv[2]) throw new Error('Usage: node scripts/update-donnieplus.mjs /path/to/leotv-web');
const destination = resolve('public/donnieplus');
const html = readFileSync(join(upstream, 'index.html'), 'utf8');
const main = html.match(/<script type="module" crossorigin src="\.\/assets\/([^"]+)"/)[1];
const css = html.match(/<link rel="stylesheet" crossorigin href="\.\/assets\/([^"]+)"/)[1];
const preloads = [...html.matchAll(/<link rel="modulepreload" crossorigin href="\.\/assets\/([^"]+)"/g)].map(m => m[1]);
const polyfill = html.match(/id="vite-legacy-polyfill" src="\.\/assets\/([^"]+)"/)[1];
const legacy = html.match(/id="vite-legacy-entry" data-src="\.\/assets\/([^"]+)"/)[1];
const commit = execFileSync('git', ['-C', upstream, 'rev-parse', 'HEAD'], { encoding: 'utf8' }).trim();
const version = JSON.parse(readFileSync(join(upstream, 'version.json'), 'utf8')).build;
const revision = version + '-nav2';
const suffix = '-dp' + revision;
const renames = Object.fromEntries(readdirSync(join(upstream, 'assets')).map(name => [name, name.replace(/(\.[^.]+)$/, suffix + '$1')]));
rmSync(join(destination, 'assets'), { recursive: true, force: true });
cpSync(join(upstream, 'assets'), join(destination, 'assets'), { recursive: true });
for (const name of readdirSync(join(destination, 'assets'))) {
  const path = join(destination, 'assets', name);
  let source = readFileSync(path, 'utf8');
  if (name.endsWith('.js')) {
    source = source.replaceAll('LeoTv', 'DonniePlus').replaceAll('LeoTV', 'DonniePlus')
      .replaceAll('DonniePlusNative', 'LeoTvNative').replaceAll('X-DonniePlus-Key', 'X-LeoTv-Key');
    source = source.replace(/\{"X-LeoTv-Key":[^}]+\}/g, '{}');
    source = source.replace(/access:\(\)=>[\w$]+\(`\/api\/access`\)/g, 'access:()=>Promise.resolve({key:null})');
    source = source.replace(/access:function\(\)\{return [\w$]+\(`\/api\/access`\)\}/g, 'access:function(){return Promise.resolve({key:null})}');
    source = source.replaceAll('`https://leoz1000.github.io/leotv-web/`', '(location.origin+`/donnieplus/`)')
      .replaceAll('"https://leoz1000.github.io/leotv-web/"', '(location.origin+"/donnieplus/")');
  }
  if (/^Account-/.test(name)) source = patchLogin(source, /^Account-legacy-/.test(name));
  source = patchAccountOnly(source, name, main, legacy);
  if (name === main || name === legacy) source = patchNavigation(source, name === legacy);
  for (const [oldName, newName] of Object.entries(renames)) source = source.replaceAll(oldName, newName);
  writeFileSync(join(destination, 'assets', renames[name]), source);
  if (renames[name] !== name) rmSync(path);
}
const bootstrap = `try{var prior=localStorage.getItem('leotv.server');var target='https://macleo.tail1ed5a3.ts.net';if(prior&&prior!==target){var account=localStorage.getItem('leotv.account@'+prior);if(account&&!localStorage.getItem('leotv.account@'+target))localStorage.setItem('leotv.account@'+target,account)}localStorage.setItem('leotv.server',target);localStorage.removeItem('leotv.accessKey');localStorage.removeItem('leotv.moveToMain');sessionStorage.setItem('leotv.intro','1')}catch(e){}window.__leotvIntroDone=1;`;
writeFileSync(join(destination, 'index.html'), `<!doctype html>
<html lang="it"><head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover,interactive-widget=resizes-content">
<meta name="color-scheme" content="dark"><meta name="theme-color" content="#050507">
<meta name="description" content="DonniePlus. Film, serie e un'interfaccia tutta tua."><title>DonniePlus</title>
<link rel="icon" type="image/svg+xml" href="./donnieplus-icon.svg"><link rel="manifest" href="./manifest.webmanifest">
<link rel="apple-touch-icon" href="./apple-touch-icon.png?v=4"><meta name="mobile-web-app-capable" content="yes"><meta name="apple-mobile-web-app-capable" content="yes"><meta name="apple-mobile-web-app-title" content="DonniePlus"><meta name="apple-mobile-web-app-status-bar-style" content="black-translucent">
<style>html,body{background:#050507;margin:0}</style>
<script data-cfasync="false">${bootstrap}</script>
<link rel="stylesheet" href="./assets/${renames[css]}"><link rel="stylesheet" href="./donnieplus-theme.css?v=${version}"><link rel="stylesheet" href="./donnieplus-settings.css?v=${version}"><link rel="stylesheet" href="./donnieplus-navigation.css?v=${revision}">
${preloads.map(file => `<link rel="modulepreload" crossorigin href="./assets/${renames[file]}">`).join('\n')}
<script data-cfasync="false" type="module" src="./donnieplus-settings.js?v=${revision}"></script>
<script data-cfasync="false" type="module" crossorigin src="./assets/${renames[main]}"></script>
</head><body><div id="root"></div>
<script data-cfasync="false" nomodule crossorigin id="vite-legacy-polyfill" src="./assets/${renames[polyfill]}"></script>
<script data-cfasync="false" nomodule crossorigin id="vite-legacy-entry" data-src="./assets/${renames[legacy]}">System.import(document.getElementById('vite-legacy-entry').getAttribute('data-src'))</script>
</body></html>`);
for (const file of ['version.json', 'apk-version.json']) cpSync(join(upstream, file), join(destination, file));
const manifest = JSON.parse(readFileSync(join(upstream, 'manifest.webmanifest'), 'utf8'));
Object.assign(manifest, { name: 'DonniePlus', short_name: 'DonniePlus', id: '/donnieplus/', start_url: '/donnieplus/', scope: '/donnieplus/', description: 'Film, serie e un’interfaccia tutta tua.' });
writeFileSync(join(destination, 'manifest.webmanifest'), JSON.stringify(manifest, null, 2) + '\n');
const settingsPath = join(destination, 'donnieplus-settings.js');
writeFileSync(settingsPath, readFileSync(settingsPath, 'utf8').replace(/(donnieplus-(?:preferences\.js|wordmark\.svg)\?v=)[^'"]+/g, '$1' + version));
const custom = [`donnieplus-navigation.css?v=${revision}`, `donnieplus-theme.css?v=${version}`, `donnieplus-settings.css?v=${version}`, `donnieplus-settings.js?v=${revision}`, `donnieplus-preferences.js?v=${version}`, `donnieplus-wordmark.svg?v=${version}`];
writeFileSync(join(destination, 'precache.json'), JSON.stringify({ files: [...Object.values(renames).map(name => 'assets/' + name), ...custom] }));
let sw = readFileSync(join(upstream, 'sw.js'), 'utf8').replace(/const CACHE = '[^']+';/, `const CACHE = 'donnieplus-pages-${revision}';`).replaceAll('LeoTv', 'DonniePlus');
sw = sw.replace('keys.filter((k) => k !== CACHE)', "keys.filter((k) => k !== CACHE && k.startsWith('donnieplus-pages-'))");
sw = sw.replace("path.includes('/assets/')", "(path.includes('/assets/') || /\\/donnieplus-(settings|preferences|theme|wordmark|navigation)[^/]*\\.(css|js|svg)$/.test(path))");
// Cloudflare caches the bare entry URL for four hours. Fetch a build-specific
// shell URL both at install and navigation so an older HTML cannot replace it.
sw = sw.replace("fetch(f, { cache: f === './' ? 'no-store' : 'default' })", "fetch(f === './' ? './?dp-build=" + revision + "' : f, { cache: f === './' ? 'no-store' : 'default' })");
sw = sw.replace('self.registration.navigationPreload.enable()', 'self.registration.navigationPreload.disable()');
sw = sw.replace('Promise.resolve(preloaded).then((early) => early || fetch(request))', "Promise.resolve().then(() => { const fresh = new URL(request.url); fresh.searchParams.set('dp-build', '" + revision + "'); return fetch(fresh, { cache: 'no-store', credentials: 'same-origin' }); })");
writeFileSync(join(destination, 'sw.js'), sw);
writeFileSync(join(destination, 'upstream.json'), JSON.stringify({ repository: 'Leoz1000/leotv-web', commit, build: version }, null, 2) + '\n');
console.log(`Updated DonniePlus from LeoTV ${version} (${commit}), ${Object.keys(renames).length} assets.`);
