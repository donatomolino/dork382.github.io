import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Window } from 'happy-dom';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import { dpCreateNavigationStage, dpScreenReady } from '../public/donnieplus/donnieplus-navigation.js';
const window = new Window();
for (const key of ['window', 'document', 'HTMLElement', 'MutationObserver']) globalThis[key] = key === 'window' ? window : window[key];
globalThis.IS_REACT_ACT_ENVIRONMENT = true;
const h = React.createElement;
const Stage = dpCreateNavigationStage(React);
const settle = () => new Promise(resolve => setTimeout(resolve, 0));

test('keeps outgoing catalog, commits only the latest request, and restores cached screen/scroll without remounting', async () => {
  const container = document.createElement('div'); document.body.append(container);
  const root = createRoot(container), pending = new Map(), mounts = new Map();
  function Entry({ entry, visible }) {
    const [loaded, setLoaded] = React.useState(entry.name === 'home');
    React.useEffect(() => { mounts.set(entry.name, (mounts.get(entry.name) || 0) + 1); pending.set(entry.name, () => setLoaded(true)); }, []);
    return h('div', { className: 'screen-slot' + (visible ? '' : ' screen-slot--hidden') }, loaded ? h('div', { className: 'screen', 'data-catalog': entry.name }, entry.name + ' catalog') : h('div', { className: 'home--skeleton', 'aria-busy': true }));
  }
  async function route(name, id) { await act(async () => { const current = { name, id }; root.render(h(Stage, { entries: [current], current, Entry })); await settle(); }); }
  const visible = () => container.querySelector('.screen-slot:not(.screen-slot--hidden)');
  await route('home', 1);
  const home = visible().querySelector('.screen'); home.scrollTop = 230;
  await route('movies', 2);
  assert.equal(visible().textContent, 'home catalog');
  assert.match(container.querySelector('[role=status]').textContent, /Caricamento Film/);
  assert.equal(container.querySelector('.dp-navigation').getAttribute('aria-busy'), 'true');
  await route('tv', 3);
  await act(async () => { pending.get('movies')(); await settle(); });
  assert.equal(visible().textContent, 'home catalog', 'late Film response must not replace the requested TV section');
  assert.match(container.querySelector('[role=status]').textContent, /Serie TV/);
  await act(async () => { pending.get('tv')(); await settle(); });
  assert.equal(visible().textContent, 'tv catalog'); assert.equal(container.querySelector('[role=status]'), null);
  await route('home', 4);
  assert.equal(visible().querySelector('.screen'), home); assert.equal(home.scrollTop, 230); assert.equal(mounts.get('home'), 1);
  await route('movies', 5);
  assert.equal(visible().textContent, 'movies catalog'); assert.equal(mounts.get('movies'), 1);
  await act(async () => root.unmount()); container.remove();
});

test('waiting and failed lazy screens keep a visible status or error instead of an empty screen', () => {
  const node = document.createElement('div');
  node.innerHTML = '<div class="screen-slot"><div class="screen screen--center"><div class="spinner"></div></div></div>';
  assert.equal(dpScreenReady(node), false);
  node.innerHTML = '<div class="screen-slot"><div class="screen"><div role="alert">Riprova</div></div></div>';
  assert.equal(dpScreenReady(node), true);
});
