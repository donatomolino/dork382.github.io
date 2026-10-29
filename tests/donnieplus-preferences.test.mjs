import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { DEFAULTS, THEMES, normalize, loadPreferences } from '../public/donnieplus/donnieplus-preferences.js';

test('corrupt or inaccessible storage safely restores defaults', () => {
  for (const getItem of [() => '{broken', () => null, () => { throw new Error('blocked'); }]) assert.deepEqual(loadPreferences({ getItem }), DEFAULTS);
});
test('imports allow only supported preferences and clamp numeric ranges', () => {
  const value = normalize({ theme: 'ocean', blur: 999, volume: -5, opacity: NaN, bold: 'true', speed: 1.5, accent: 'url(unsafe)', token: 'secret' });
  assert.equal(value.theme, 'ocean'); assert.equal(value.blur, 32); assert.equal(value.volume, 0); assert.equal(value.opacity, DEFAULTS.opacity); assert.equal(value.bold, false); assert.equal(value.speed, '1.5'); assert.equal(value.accent, 'auto'); assert.equal('token' in value, false);
  for (const input of [null, [], 'bad']) assert.deepEqual(normalize(input), DEFAULTS);
  for (const [theme] of THEMES) assert.equal(normalize({ theme }).theme, theme);
});
test('all precached assets exist and both API bundles use account-only authentication', () => {
  const base = new URL('../public/donnieplus/', import.meta.url);
  const manifest = JSON.parse(readFileSync(new URL('precache.json', base)));
  for (const path of manifest.files) assert.ok(existsSync(new URL(path.split('?')[0], base)), path);
  const assets = readdirSync(new URL('assets/', base));
  let accountBundles = 0;
  for (const file of assets.filter(f => f.endsWith('.js'))) {
    const source = readFileSync(new URL('assets/' + file, base), 'utf8');
    assert.ok(!/\b[\w$]+\.access\(\)/.test(source), file + ' contains server key consumer'); assert.ok(!source.includes('Codice di accesso'), file + ' contains server key UI'); assert.ok(!source.includes('set-srv-key'), file); assert.ok(!source.includes('X-LeoTv-Key'), file); assert.ok(!source.includes('/api/access'), file);
    if (source.includes('/api/account/login')) { accountBundles++; assert.ok(source.includes('X-Account-Token')); assert.ok(source.includes('X-Profile-Token')); }
    for (const ref of source.matchAll(/(?:\.\/|\/assets\/)([\w-]+\.(?:js|css))(?![a-z])/g)) assert.ok(assets.includes(ref[1]) || existsSync(new URL(ref[1], base)), file + ' references ' + ref[1]);
  }
  assert.equal(accountBundles, 2);
});

test('shell refresh bypasses stale CDN HTML and isolates cache deletion', () => {
  const sw = readFileSync(new URL('../public/donnieplus/sw.js', import.meta.url), 'utf8');
  assert.ok(sw.includes("fresh.searchParams.set('dp-build'"));
  assert.ok(sw.includes("./?dp-build="));
  assert.ok(sw.includes("k.startsWith('donnieplus-pages-')"));
  assert.ok(!sw.includes('early || fetch(request)'));
});
