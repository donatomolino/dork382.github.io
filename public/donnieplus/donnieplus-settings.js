import { DEFAULTS, THEMES, STORAGE_KEY, normalize, loadPreferences } from './donnieplus-preferences.js?v=20261029160502';

const root = document.documentElement;
let prefs = loadPreferences({ getItem: key => localStorage.getItem(key) });
let dialog, opener, activeTab = 'appearance', saveTimer;
const media = new WeakSet();
const accents = { violet: '#7b5cff', blue: '#5caaff', mint: '#56d9b2', amber: '#ffc46b', pink: '#f78cbe' };
const icon = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" aria-hidden="true"><path d="m12 2 2.5 6.5L21 11l-6.5 2.5L12 20l-2.5-6.5L3 11l6.5-2.5L12 2Z"/><path d="m20 16 1 2.5 2.5 1-2.5 1L20 23l-1-2.5-2.5-1 2.5-1L20 16Z"/></svg>';
const escape = value => String(value).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);

function apply() {
  for (const key of ['theme','quality','density','text','motion','fit']) root.dataset['dp' + key[0].toUpperCase() + key.slice(1)] = prefs[key];
  for (const key of ['bold','contrast','hover','ratings','synopsis','captions','tabLabels']) root.dataset['dp' + key[0].toUpperCase() + key.slice(1)] = String(prefs[key]);
  const theme = THEMES.find(t => t[0] === prefs.theme);
  const accent = prefs.accent === 'auto' ? theme[3] : accents[prefs.accent];
  root.style.setProperty('--dp-accent', accent);
  const rgb = accent.slice(1).match(/../g).map(c => parseInt(c, 16) / 255).map(c => c <= .04045 ? c / 12.92 : ((c + .055) / 1.055) ** 2.4);
  root.style.setProperty('--dp-on-accent', rgb[0] * .2126 + rgb[1] * .7152 + rgb[2] * .0722 > .179 ? '#10121c' : '#fff');
  root.style.setProperty('--dp-blur', prefs.blur + 'px');
  root.style.setProperty('--dp-opacity', prefs.opacity / 100);
  root.style.setProperty('--dp-radius', prefs.radius + 'px');
  root.style.setProperty('--dp-brightness', prefs.brightness / 100);
  root.style.setProperty('--dp-saturation', prefs.saturation / 100);
  root.style.setProperty('--dp-text-scale', { normal: 1, large: 1.12, xlarge: 1.25 }[prefs.text]);
  root.style.colorScheme = prefs.theme === 'light' ? 'light' : 'dark';
  document.querySelector('meta[name="theme-color"]')?.setAttribute('content', theme[4]);
  document.querySelectorAll('video').forEach(video => { video.style.objectFit = prefs.fit; });
}
function persist() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => { try { localStorage.setItem(STORAGE_KEY, JSON.stringify(prefs)); } catch {} }, 100);
}
function announce(message) { const status = dialog?.querySelector('#dp-status'); if (status) status.textContent = message; }
function update(key, value) {
  prefs = normalize({ ...prefs, [key]: value });
  apply(); persist();
  if (['volume','speed','muted','fit'].includes(key)) document.querySelectorAll('video').forEach(video => {
    try {
      if (key === 'volume') video.volume = prefs.volume / 100;
      if (key === 'speed') video.playbackRate = Number(prefs.speed);
      if (key === 'muted') video.muted = prefs.muted;
      if (key === 'fit') video.style.objectFit = prefs.fit;
    } catch {}
  });
  announce('Preferenze salvate su questo dispositivo.');
}
function toggle(key, title, description) { return `<label class="dp-control dp-toggle"><span><strong>${title}</strong><small>${description}</small></span><input data-pref="${key}" type="checkbox" ${prefs[key] ? 'checked' : ''}><span class="dp-switch" aria-hidden="true"></span></label>`; }
function select(key, title, description, options) { return `<label class="dp-control"><span><strong>${title}</strong><small>${description}</small></span><select data-pref="${key}" aria-label="${title}">${options.map(([value,label]) => `<option value="${value}" ${prefs[key] === value ? 'selected' : ''}>${label}</option>`).join('')}</select></label>`; }
function range(key, title, description, min, max, unit = '%') { return `<label class="dp-control dp-range"><span><strong>${title}</strong><small>${description}</small></span><span class="dp-range-input"><input data-pref="${key}" type="range" min="${min}" max="${max}" value="${prefs[key]}" aria-label="${title}"><output data-unit="${unit}">${prefs[key]}${unit}</output></span></label>`; }
function group(title, body) { return `<section class="dp-group"><h3>${title}</h3><div class="dp-controls">${body}</div></section>`; }

function render() {
  dialog.innerHTML = `<div class="dp-panel"><header class="dp-panel-header"><div class="dp-panel-brand">${icon}<span>DonniePlus <small>Il tuo spazio, a modo tuo.</small></span></div><button type="button" class="dp-close" data-action="close" aria-label="Chiudi personalizzazione">×</button></header>
  <div class="dp-panel-intro"><span class="dp-eyebrow">PERSONALIZZAZIONE</span><h2 id="dp-title">Rendilo tuo.</h2><p>Scegli l’atmosfera e i dettagli della tua esperienza.</p></div>
  <nav class="dp-tabs" role="tablist" aria-label="Categorie delle preferenze">${[['appearance','Aspetto'],['interface','Interfaccia'],['video','Video'],['data','Dati e preferenze']].map(([id,label]) => `<button type="button" role="tab" id="dp-tab-${id}" data-tab="${id}" aria-controls="dp-${id}" aria-selected="${activeTab === id}" tabindex="${activeTab === id ? 0 : -1}">${label}</button>`).join('')}</nav>
  <div class="dp-panel-scroll">
  <section role="tabpanel" id="dp-appearance" aria-labelledby="dp-tab-appearance" ${activeTab !== 'appearance' ? 'hidden' : ''}>
  <section class="dp-group"><div class="dp-group-heading"><h3>Un’atmosfera per ogni momento</h3><span>6 temi</span></div><div class="dp-theme-grid">${THEMES.map(([id,name,description,accent,bg]) => `<button type="button" class="dp-theme-card" data-theme="${id}" aria-pressed="${prefs.theme === id}" style="--sample-accent:${accent};--sample-bg:${bg}"><span class="dp-theme-preview" aria-hidden="true"><i></i><span class="dp-sample-bar"></span><span class="dp-sample-poster"></span><span class="dp-sample-poster"></span><span class="dp-sample-poster"></span></span><strong>${name}<span class="dp-theme-check" aria-hidden="true">✓</span></strong><small>${description}</small></button>`).join('')}</div></section>
  ${group('Vetro e dettagli', select('accent','Colore d’accento','Pulsanti, selezioni e dettagli.',[['auto','Dal tema'],['violet','Viola'],['blue','Blu'],['mint','Menta'],['amber','Ambra'],['pink','Rosa']]) + range('blur','Sfocatura del vetro','Da trasparente a morbido. Disponibile nei temi in vetro.',0,32,' px') + range('opacity','Opacità dei pannelli','Aumentala per leggere meglio sopra le immagini.',40,100) + range('radius','Angoli','Dai bordi netti alle forme arrotondate.',0,28,' px') + select('quality','Effetti grafici','Riduci gli effetti sui dispositivi meno potenti.',[['full','Completi'],['balanced','Bilanciati'],['light','Leggeri']]))}
  </section>
  <section role="tabpanel" id="dp-interface" aria-labelledby="dp-tab-interface" ${activeTab !== 'interface' ? 'hidden' : ''}>
  ${group('Layout e leggibilità', select('density','Dimensione delle copertine','Più titoli a schermo o immagini più grandi.',[['compact','Compatte'],['normal','Standard'],['wide','Grandi']]) + select('text','Dimensione dei testi','Ingrandisce etichette, pulsanti e descrizioni.',[['normal','Standard'],['large','Grande'],['xlarge','Molto grande']]) + toggle('bold','Testi più marcati','Rende più evidenti titoli ed etichette.') + toggle('contrast','Contrasto elevato','Testi secondari più chiari e bordi più visibili.') + select('motion','Animazioni','Controlla transizioni e movimenti dell’interfaccia.',[['system','Come il dispositivo'],['reduced','Ridotte'],['full','Complete']]) + toggle('hover','Effetto sulle copertine','Piccolo ingrandimento al passaggio del puntatore.'))}
  ${group('Contenuti visibili', toggle('ratings','Mostra valutazioni','Visualizza i voti accanto ai contenuti.') + toggle('synopsis','Mostra descrizioni','Mantiene visibili le trame nella home e nei dettagli.') + toggle('captions','Mostra titoli delle copertine','Visualizza le etichette sotto le immagini.') + toggle('tabLabels','Etichette della navigazione','Mostra il nome sotto le icone su mobile.') + range('brightness','Luminosità dello sfondo','Regola le immagini di copertina della home, non i video.',50,120) + range('saturation','Colore delle copertine','Da toni neutri a immagini più vivaci.',0,150))}
  </section>
  <section role="tabpanel" id="dp-video" aria-labelledby="dp-tab-video" ${activeTab !== 'video' ? 'hidden' : ''}>
  <p class="dp-note">Queste preferenze si applicano al player su questo dispositivo. I controlli del player restano disponibili durante la visione.</p>
  ${group('Audio e riproduzione', range('volume','Volume iniziale','Il volume del player; su iPhone può essere gestito dal sistema.',0,100) + toggle('rememberVolume','Ricorda il volume','Memorizza le modifiche fatte dal player.') + toggle('muted','Avvia senza audio','Inizia la riproduzione con l’audio disattivato.') + select('speed','Velocità iniziale','Scegli la velocità predefinita dei video.',[['0.75','0,75×'],['1','Normale · 1×'],['1.25','1,25×'],['1.5','1,5×'],['2','2×']]) + select('fit','Formato video','Adatta conserva tutto il fotogramma; Riempi può ritagliarlo.',[['contain','Adatta allo schermo'],['cover','Riempi lo schermo']]))}
  </section>
  <section role="tabpanel" id="dp-data" aria-labelledby="dp-tab-data" ${activeTab !== 'data' ? 'hidden' : ''}>
  ${group('Le tue preferenze', `<div class="dp-data-row"><div><strong>Porta con te il tuo stile</strong><p>Esporta o importa solo queste preferenze grafiche e video.</p></div><div class="dp-data-actions"><button type="button" data-action="export">Esporta</button><button type="button" data-action="import">Importa</button><input type="file" id="dp-import" accept="application/json,.json" hidden></div></div><div class="dp-data-row"><div><strong>Cache dell’interfaccia</strong><p>Rimuove le copie locali dei file di DonniePlus. Il primo caricamento sarà più lento.</p></div><button type="button" data-action="cache">Svuota cache</button></div><div class="dp-data-row"><div><strong>Ripristina le preferenze</strong><p>Torna al tema Liquid Glass e ai valori iniziali. Account, lista e progressi vengono conservati.</p></div><button type="button" class="dp-danger" data-action="reset">Ripristina</button></div>`)}
  <div class="dp-about"><strong>DonniePlus</strong><span id="dp-version">Preferenze di questo dispositivo</span></div>
  </section></div>
  <footer class="dp-panel-footer"><span role="status" aria-live="polite" id="dp-status">Le modifiche vengono salvate automaticamente.</span><button type="button" data-action="close">Fatto</button></footer></div>`;
  dialog.querySelector('#dp-import').addEventListener('change', async event => {
    const file = event.target.files?.[0]; if (!file) return;
    if (file.size > 50000) { announce('Il file è troppo grande. Seleziona un’esportazione di DonniePlus.'); return; }
    try {
      const data = JSON.parse(await file.text());
      if (data.app !== 'DonniePlus' || data.version !== 1 || !data.preferences || typeof data.preferences !== 'object') throw new Error('invalid');
      prefs = normalize(data.preferences); apply(); persist(); render(); announce('Preferenze importate.');
    } catch { announce('Questo file non contiene preferenze DonniePlus valide.'); }
  });
  fetch(new URL('./upstream.json', import.meta.url)).then(r => r.json()).then(data => { const el = dialog.querySelector('#dp-version'); if (el) el.textContent = 'LeoTV · build ' + escape(data.build) + ' · personalizzazione DonniePlus'; }).catch(() => {});
}
function chooseTab(tab) {
  activeTab = tab;
  dialog.querySelectorAll('[data-tab]').forEach(button => { const on = button.dataset.tab === tab; button.setAttribute('aria-selected', String(on)); button.tabIndex = on ? 0 : -1; });
  dialog.querySelectorAll('[role="tabpanel"]').forEach(panel => { panel.hidden = panel.id !== 'dp-' + tab; });
  dialog.querySelector('.dp-panel-scroll').scrollTop = 0;
}
function openPanel(event) {
  opener = event?.currentTarget || document.activeElement;
  if (!dialog) {
    dialog = document.createElement('dialog'); dialog.id = 'dp-preferences'; dialog.setAttribute('aria-labelledby', 'dp-title');
    document.body.append(dialog);
    dialog.addEventListener('close', () => { if (opener?.isConnected) opener.focus({ preventScroll: true }); });
    dialog.addEventListener('click', event => {
      if (event.target === dialog) { dialog.close(); return; }
      const theme = event.target.closest('[data-theme]');
      if (theme) { update('theme', theme.dataset.theme); dialog.querySelectorAll('[data-theme]').forEach(button => button.setAttribute('aria-pressed', String(button === theme))); }
      const tab = event.target.closest('[data-tab]'); if (tab) chooseTab(tab.dataset.tab);
      const action = event.target.closest('[data-action]')?.dataset.action;
      if (action === 'close') dialog.close();
      if (action === 'reset') { prefs = { ...DEFAULTS }; apply(); persist(); render(); announce('Preferenze ripristinate.'); }
      if (action === 'export') {
        const blob = new Blob([JSON.stringify({ app: 'DonniePlus', version: 1, preferences: prefs }, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob), link = document.createElement('a'); link.href = url; link.download = 'DonniePlus-preferenze.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000); announce('Preferenze esportate.');
      }
      if (action === 'import') dialog.querySelector('#dp-import').click();
      if (action === 'cache') clearInterfaceCache(event.target.closest('button'));
    });
    dialog.addEventListener('input', event => {
      const input = event.target, key = input.dataset.pref; if (!key) return;
      const value = input.type === 'checkbox' ? input.checked : input.type === 'range' ? Number(input.value) : input.value;
      if (input.type === 'range') { const output = input.parentElement.querySelector('output'); output.textContent = input.value + output.dataset.unit; }
      update(key, value);
    });
    dialog.addEventListener('keydown', event => {
      event.stopPropagation();
      if (event.key === 'Escape') { event.preventDefault(); dialog.close(); }
      if (event.target.matches('[role="tab"]') && ['ArrowLeft','ArrowRight','Home','End'].includes(event.key)) {
        event.preventDefault(); const buttons = [...dialog.querySelectorAll('[data-tab]')], index = buttons.indexOf(event.target);
        const next = event.key === 'Home' ? 0 : event.key === 'End' ? buttons.length - 1 : (index + (event.key === 'ArrowRight' ? 1 : -1) + buttons.length) % buttons.length;
        chooseTab(buttons[next].dataset.tab); buttons[next].focus();
      }
    });
  }
  render(); if (!dialog.open) dialog.showModal();
}
async function clearInterfaceCache(button) {
  button.disabled = true;
  try {
    if (!('caches' in window)) throw new Error('unavailable');
    const names = (await caches.keys()).filter(name => name.startsWith('donnieplus-pages-'));
    await Promise.all(names.map(name => caches.delete(name)));
    announce('Cache dell’interfaccia svuotata.');
  } catch { announce('La cache non è disponibile in questo browser.'); }
  finally { button.disabled = false; }
}
function initializeVideo(video) {
  if (media.has(video)) return; media.add(video);
  let initializing = false, initialized = false;
  function initialize() {
    initializing = true;
    try { video.volume = prefs.volume / 100; video.muted = prefs.muted; video.playbackRate = Number(prefs.speed); video.style.objectFit = prefs.fit; } catch {}
    initialized = true; setTimeout(() => { initializing = false; }, 0);
  }
  video.addEventListener('loadedmetadata', initialize);
  video.addEventListener('play', () => { if (!initialized) initialize(); });
  video.addEventListener('volumechange', () => {
    if (prefs.rememberVolume && initialized && !initializing && !video.muted) { prefs.volume = Math.round(video.volume * 100); persist(); }
  });
  video.style.objectFit = prefs.fit;
}
function button(className, label) {
  const button = document.createElement('button'); button.type = 'button'; button.className = className; button.innerHTML = icon + `<span>${label}</span>`;
  button.addEventListener('click', openPanel); return button;
}
function enhance() {
  document.querySelectorAll('.auth__brand,.auth__icon--app').forEach(brand => {
    if (!brand.querySelector('.dp-wordmark')) {
      const image = document.createElement('img'); image.className = 'dp-wordmark'; image.src = new URL('./donnieplus-wordmark.svg?v=20261029160502', import.meta.url).href; image.alt = 'DonniePlus'; image.width = 260; image.height = 50; brand.replaceChildren(image); brand.removeAttribute('aria-hidden'); brand.setAttribute('aria-label', 'DonniePlus');
    }
  });
  document.querySelectorAll('.auth__panel').forEach(panel => { if (!panel.querySelector('.dp-login-preferences')) panel.append(button('dp-login-preferences', 'Aspetto e preferenze')); });
  document.querySelectorAll('.settings__content').forEach(content => {
    if (content.querySelector('.dp-settings-entry')) return;
    const section = document.createElement('section'); section.className = 'dp-settings-entry';
    section.innerHTML = '<h2>PERSONALIZZAZIONE</h2>';
    const entry = button('dp-entry', '<strong>Aspetto e preferenze</strong><small>Temi, vetro, interfaccia e video</small>');
    section.append(entry); const title = content.querySelector('.settings__title'); content.insertBefore(section, title ? title.nextSibling : content.firstChild);
  });
  document.querySelectorAll('video').forEach(initializeVideo);
}
apply();
let scheduled = false;
new MutationObserver(() => { if (!scheduled) { scheduled = true; requestAnimationFrame(() => { scheduled = false; enhance(); }); } }).observe(document.getElementById('root'), { childList: true, subtree: true });
enhance();
window.addEventListener('storage', event => { if (event.key === STORAGE_KEY) { prefs = loadPreferences({ getItem: key => localStorage.getItem(key) }); apply(); if (dialog?.open) render(); } });
window.addEventListener('pagehide', () => { try { localStorage.setItem(STORAGE_KEY, JSON.stringify(prefs)); } catch {} });
// The upstream load handler registers the bare worker URL. Register a build-
// specific URL afterwards so the CDN cannot keep serving an older worker.
if ('serviceWorker' in navigator && window.isSecureContext) {
  const registerWorker = () => setTimeout(() => {
    const build = new URL(import.meta.url).searchParams.get('v') || '20261029160502';
    const worker = new URL('./sw.js', import.meta.url); worker.searchParams.set('dp-build', build + '-shell2');
    navigator.serviceWorker.register(worker.href, { scope: './', updateViaCache: 'none' }).catch(() => {});
  }, 0);
  if (document.readyState === 'complete') registerWorker();
  else window.addEventListener('load', registerWorker, { once: true });
}
