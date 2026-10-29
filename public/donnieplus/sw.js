// DonniePlus service worker. It makes the web app installable and lets it open
// quickly, also with no network at all (on a plane, to watch downloads). It
// never stores the catalog, the account data or the video:
// - pages: from the network, the saved copy when offline or after 3 s;
// - hashed files in assets/: the whole build is saved when the worker
//   installs (precache.json lists it), so every screen opens offline;
// - api, video and version checks: never touched.
// It also shows the "è pronto" notices of downloads (push, below).
const CACHE = 'donnieplus-pages-20261029160502-shell2';
const MAX_FILES = 60;
// the page, its icon and every file of the current build
const SHELL = ['./', 'icon.png?v=3', 'manifest.webmanifest'];

async function saveBuild() {
  const cache = await caches.open(CACHE);
  let files = [];
  try {
    const list = await fetch('precache.json', { cache: 'no-store' });
    if (list.ok) files = ((await list.json()).files || []).filter((f) => typeof f === 'string');
  } catch (e) { /* offline or an old server: files are saved as they are used */ }
  // one missing file must not stop the others (or the worker) from installing
  const have = new Set((await cache.keys()).map((r) => r.url));
  // the page itself is always fetched again: it names the current build
  const wanted = SHELL.concat(files).filter((f) => f === './' || !have.has(new URL(f, self.registration.scope).href));
  await Promise.all(wanted.map((f) =>
    fetch(f === './' ? './?dp-build=20261029160502-ui2' : f, { cache: f === './' ? 'no-store' : 'default' })
      .then((r) => (r.ok ? cache.put(f, r) : null))
      .catch(() => null)));
  if (files.length) await cache.put('precache.json', new Response(JSON.stringify({ files })));
}

self.addEventListener('install', (event) => {
  event.waitUntil(saveBuild().then(() => self.skipWaiting()));
});

// The worker itself rarely changes, but every new build of the app does: the
// page asks for its files to be saved each time it opens (see pwa.ts), and
// only the files not saved yet are downloaded.
self.addEventListener('message', (event) => {
  if (event.data && event.data.type === 'save-build') {
    event.waitUntil(saveBuild().then(() => caches.open(CACHE)).then(trim));
  }
});

self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.filter((k) => k !== CACHE && k.startsWith('donnieplus-pages-')).map((k) => caches.delete(k))))
      // the page request starts while the worker itself is waking up
      .then(() => (self.registration.navigationPreload ? self.registration.navigationPreload.disable() : null))
      .then(() => self.clients.claim()),
  );
});

async function trim(cache) {
  // files of older builds go; the current build always stays complete
  let keep = null;
  try {
    const saved = await cache.match('precache.json');
    if (saved) keep = new Set(((await saved.json()).files || []).map((f) => new URL(f, self.registration.scope).href));
  } catch (e) { keep = null; }
  const assets = (await cache.keys()).filter((r) => r.url.includes('/assets/'));
  if (keep && keep.size) {
    const old = assets.filter((r) => !keep.has(r.url));
    // keep a few recent extras (a page still running the previous build)
    for (let i = 0; i < old.length - 20; i++) await cache.delete(old[i]);
    return;
  }
  for (let i = 0; i < assets.length - MAX_FILES; i++) await cache.delete(assets[i]);
}

// A slow network (a weak Wi-Fi, a sleeping server) must not leave a black
// screen: after this long the saved copy of the page opens instead, and the
// network answer, when it comes, still refreshes that copy for next time.
const NAVIGATION_TIMEOUT = 3000;

async function fromNetworkFirst(request, patient, preloaded) {
  const cache = await caches.open(CACHE);
  const network = Promise.resolve().then(() => { const fresh = new URL(request.url); fresh.searchParams.set('dp-build', '20261029160502-ui2'); return fetch(fresh, { cache: 'no-store', credentials: 'same-origin' }); }).then(async (response) => {
    // only the app's own page is kept for offline use, never another HTML
    // file of the site (an .html page would replace the app when offline)
    const appPage = !/\.[a-z0-9]+$/i.test(new URL(request.url).pathname);
    if (appPage && response.ok && (response.headers.get('content-type') || '').includes('text/html')) await cache.put('./', response.clone());
    return response;
  });
  network.catch(() => {});   // a late failure after the saved page answered is not an error
  // right after "Aggiorna" the page waits for the network: the saved copy is the old build
  const slow = patient ? new Promise(() => {}) : new Promise((resolve) => setTimeout(resolve, NAVIGATION_TIMEOUT)).then(() => cache.match('./'));
  try {
    // the first useful answer wins: the network, or the saved page after 3 s
    const first = await Promise.race([network, slow.then((saved) => saved || network)]);
    return first;
  } catch (err) {
    const saved = await cache.match('./');
    if (saved) return saved;
    throw err;
  }
}

async function fromCacheFirst(request) {
  const cache = await caches.open(CACHE);
  const saved = await cache.match(request);
  if (saved) return saved;
  const response = await fetch(request);
  if (response.ok) {
    await cache.put(request, response.clone());
    trim(cache);
  }
  return response;
}

self.addEventListener('fetch', (event) => {
  const request = event.request;
  if (request.method !== 'GET' || request.headers.has('range')) return;
  const url = new URL(request.url);
  if (url.origin !== self.location.origin) return;
  const path = url.pathname;
  if (path.includes('/api/') || path.endsWith('/version.json') || path.endsWith('/sw.js')) return;
  // the app download (and its short page) go straight to the network
  if (path.endsWith('.apk') || /\/apk\/?$/.test(path)) return;
  if (request.mode === 'navigate') {
    event.respondWith(fromNetworkFirst(request, url.searchParams.has('upd'), event.preloadResponse));
  } else if ((path.includes('/assets/') || /\/donnieplus-(settings|preferences|theme|wordmark)[^/]*\.(css|js|svg)$/.test(path))) {
    event.respondWith(fromCacheFirst(request));
  }
});

// "È pronto" notices for downloads, sent by the DonniePlus at home (Web Push).
// iOS requires every push to show a notification, so one is always shown.
self.addEventListener('push', (event) => {
  let data = {};
  try { data = event.data ? event.data.json() : {}; } catch (e) { data = { body: event.data ? event.data.text() : '' }; }
  event.waitUntil(self.registration.showNotification(data.title || 'DonniePlus', {
    body: data.body || 'Apri DonniePlus per finire il download.',
    tag: data.tag || 'leotv',
    icon: 'icon.png',
    badge: 'icon.png',
    data: { open: data.open || 'downloads' },
  }));
});

// A tap on the notice opens DonniePlus on the downloads, where the copy resumes by itself.
self.addEventListener('notificationclick', (event) => {
  event.notification.close();
  const screen = (event.notification.data && event.notification.data.open) || 'downloads';
  const url = new URL('./?open=' + encodeURIComponent(screen), self.registration.scope).href;
  event.waitUntil(
    self.clients.matchAll({ type: 'window', includeUncontrolled: true }).then((windows) => {
      for (const client of windows) {
        if ('focus' in client) {
          client.postMessage({ type: 'open', screen });
          return client.focus();
        }
      }
      return self.clients.openWindow(url);
    }),
  );
});
