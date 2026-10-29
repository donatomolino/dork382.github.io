export const STORAGE_KEY = 'donnieplus.preferences.v1';
export const THEMES = [
  ['glass', 'Liquid Glass', 'Il vetro originale, con accenti viola.', '#7b5cff', '#141419'],
  ['ocean', 'Oceano', 'Vetro blu e riflessi acquamarina.', '#60d8fa', '#091b29'],
  ['minimal', 'Minimal', 'Superfici opache e dettagli essenziali.', '#d6dce8', '#171a20'],
  ['cinema', 'Cinema', 'Nero profondo e accenti rossi.', '#ff5c68', '#090909'],
  ['light', 'Alba', 'Superfici chiare e vetro luminoso.', '#6154df', '#edf0f8'],
  ['neon', 'Neon', 'Viola elettrico su blu notte.', '#dc8aff', '#100c25'],
];
export const DEFAULTS = Object.freeze({ theme: 'glass', accent: 'auto', blur: 20, opacity: 62, radius: 16, quality: 'full', density: 'normal', text: 'normal', bold: false, contrast: false, motion: 'system', hover: true, ratings: true, synopsis: true, captions: true, tabLabels: true, brightness: 100, saturation: 100, volume: 80, rememberVolume: true, muted: false, speed: '1', fit: 'contain' });
const ENUMS = { theme: THEMES.map(t => t[0]), accent: ['auto', 'violet', 'blue', 'mint', 'amber', 'pink'], quality: ['full', 'balanced', 'light'], density: ['normal', 'compact', 'wide'], text: ['normal', 'large', 'xlarge'], motion: ['system', 'reduced', 'full'], speed: ['0.75', '1', '1.25', '1.5', '2'], fit: ['contain', 'cover'] };
const RANGES = { blur: [0, 32], opacity: [40, 100], radius: [0, 28], brightness: [50, 120], saturation: [0, 150], volume: [0, 100] };
export function normalize(input) {
  const result = { ...DEFAULTS };
  if (!input || typeof input !== 'object' || Array.isArray(input)) return result;
  for (const key of Object.keys(DEFAULTS)) {
    if (ENUMS[key] && ENUMS[key].includes(String(input[key]))) result[key] = String(input[key]);
    else if (RANGES[key] && typeof input[key] === 'number' && Number.isFinite(input[key])) result[key] = Math.round(Math.max(RANGES[key][0], Math.min(RANGES[key][1], input[key])));
    else if (typeof DEFAULTS[key] === 'boolean' && typeof input[key] === 'boolean') result[key] = input[key];
  }
  return result;
}
export function loadPreferences(storage) {
  try { return normalize(JSON.parse(storage.getItem(STORAGE_KEY))); } catch { return { ...DEFAULTS }; }
}
