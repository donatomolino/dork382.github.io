// Shared by the modern/legacy bundles and the lifecycle tests. No account data
// is copied: mounted screens belong to the current profile's app instance.
export const DP_ROOT_SECTIONS = new Set(['home', 'movies', 'tv', 'categories', 'search', 'mylist', 'history', 'channels', 'settings', 'downloads', 'party']);
export function dpRouteKey(entry) {
  return DP_ROOT_SECTIONS.has(entry.name) ? entry.name : 'screen:' + entry.id;
}
export function dpScreenReady(node) {
  if (!node || !node.querySelector('.screen-slot')) return false;
  return !node.querySelector('.home--skeleton, .screen--center .spinner, .page > .spinner, [aria-busy="true"]');
}
export function dpCreateNavigationStage(React) {
  const h = React.createElement;
  return function DonniePlusNavigation({ entries, current, Entry, homeButton }) {
    const retained = React.useRef(new Map());
    const ready = React.useRef(new Set());
    const nodes = React.useRef(new Map());
    const [displayed, setDisplayed] = React.useState(null);
    const requested = dpRouteKey(current);
    // Root screens keep their component identity, scroll and data when the
    // upstream router replaces its stack with a new entry/id.
    entries.forEach(entry => {
      const key = dpRouteKey(entry);
      if (!retained.current.has(key) || !DP_ROOT_SECTIONS.has(entry.name)) retained.current.set(key, entry);
    });
    const present = new Set(entries.map(dpRouteKey));
    for (const [key, entry] of retained.current) {
      if (!DP_ROOT_SECTIONS.has(entry.name) && !present.has(key)) retained.current.delete(key);
    }
    const previous = retained.current.get(displayed);
    const keepPrevious = previous && DP_ROOT_SECTIONS.has(previous.name) && DP_ROOT_SECTIONS.has(current.name);
    const pending = DP_ROOT_SECTIONS.has(current.name) && !ready.current.has(requested);
    const shown = pending && keepPrevious ? displayed : requested;
    React.useLayoutEffect(() => {
      const node = nodes.current.get(requested);
      if (!node) return;
      let stopped = false;
      const check = () => {
        if (!stopped && dpScreenReady(node)) {
          ready.current.add(requested);
          setDisplayed(requested);
        }
      };
      const observer = new MutationObserver(check);
      observer.observe(node, { childList: true, subtree: true, attributes: true, attributeFilter: ['class', 'aria-busy'] });
      check();
      return () => { stopped = true; observer.disconnect(); };
    }, [requested, current.id]);
    const labels = { home: 'Home', movies: 'Film', tv: 'Serie TV', categories: 'Categorie', search: 'Cerca', mylist: 'La mia lista', history: 'Cronologia', channels: 'Canali', settings: 'Impostazioni', downloads: 'Download', party: 'Guarda insieme' };
    return h('div', { className: 'dp-navigation', 'data-dp-section': current.name, 'data-dp-pending': String(pending), 'aria-busy': pending },
      Array.from(retained.current, ([key, entry]) => h('div', {
        key, className: 'dp-navigation__entry',
        ref: node => {
          if (!node) { nodes.current.delete(key); return; }
          node.inert = shown !== key;
          node.querySelectorAll('.home__spotlight,.home__rows,.page > :not(.home__head):not(.home__bar)').forEach(part => { part.inert = pending && shown === key; });
          node.querySelectorAll('.home__kinds > .kind-chip').forEach((chip, index) => chip.setAttribute('aria-pressed', String(['home', 'movies', 'tv'][index] === current.name)));
          nodes.current.set(key, node);
        },
        'data-dp-entry': key, 'aria-hidden': shown !== key || undefined,
        // The outgoing page stays painted but cannot trigger an old title
        // while the newly selected section is loading. Navigation stays usable.
      }, h(Entry, { entry, active: key === requested, visible: key === shown, homeButton }))),
      pending ? h('div', { className: 'dp-navigation__loading', role: 'status', 'aria-live': 'polite', 'aria-atomic': 'true' },
        h('span', { className: 'dp-navigation__progress', 'aria-hidden': 'true' }),
        h('span', null, 'Caricamento ' + (labels[current.name] || 'della sezione') + '…'),
        h('small', null, 'La pagina attuale resta visibile.')) : null);
  };
}
