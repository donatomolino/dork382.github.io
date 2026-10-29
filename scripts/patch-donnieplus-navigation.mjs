import { readFileSync } from 'node:fs';
import { transformSync } from 'esbuild';

// Fail closed when a future upstream build changes these component contracts.
function once(source, before, after) {
  if (source.split(before).length !== 2) throw new Error('Navigation patch no longer matches upstream: ' + before.slice(0, 90));
  return source.replace(before, after);
}
export function patchNavigation(source, legacy = false) {
  const factory = readFileSync(new URL('../public/donnieplus/donnieplus-navigation.js', import.meta.url), 'utf8').replace(/^export /gm, '');
  if (legacy) {
    source = once(source, 'function Lc(){let e=Pe()', 'var dpNavigationStage=dpCreateNavigationStage(V);function Lc(){let e=Pe()');
    source = once(source, 'children:u.map(e=>(0,H.jsx)(Ui,{entry:e,active:e.id===r.id,homeButton:e.name!==`home`||t.length!==1},e.id))', 'children:(0,H.jsx)(dpNavigationStage,{entries:u,current:r,Entry:Ui,homeButton:r.name!==`home`||t.length!==1})');
    source = once(source, 'function({entry:e,active:t,homeButton:n})', 'function({entry:e,active:t,visible:dpVisible=t,homeButton:n})');
    source = once(source, '(Gi,{active:t,nodeRef:i,', '(Gi,{active:dpVisible,nodeRef:i,');
    return transformSync(factory, { target: 'es2015' }).code + source;
  }
  source = once(source, 'function Oc(){let e=Ne()', 'var dpNavigationStage=dpCreateNavigationStage(B);function Oc(){let e=Ne()');
  source = once(source, 'children:u.map(e=>(0,V.jsx)(Ac,{entry:e,active:e.id===r.id,homeButton:e.name!==`home`||t.length!==1},e.id))', 'children:(0,V.jsx)(dpNavigationStage,{entries:u,current:r,Entry:Ac,homeButton:r.name!==`home`||t.length!==1})');
  source = once(source, 'function({entry:e,active:t,homeButton:n})', 'function({entry:e,active:t,visible:dpVisible=t,homeButton:n})');
  source = once(source, '(Mc,{active:t,nodeRef:i,', '(Mc,{active:dpVisible,nodeRef:i,');
  return factory + source;
}
