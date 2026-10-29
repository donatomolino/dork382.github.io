import { readFileSync } from 'node:fs';
function once(source, before, after) {
  if (source.split(before).length !== 2) throw new Error('Login patch no longer matches upstream: ' + before.slice(0, 70));
  return source.replace(before, after);
}
export function patchLogin(source, legacy) {
  const helper = readFileSync(new URL('../public/donnieplus/donnieplus-auth.js', import.meta.url), 'utf8').replace(/^export /gm, '');
  const submitStart = legacy ? 'G=async()=>{if(!k){' : 'K=async()=>{if(!j){';
  const submitEnd = legacy ? '},[q,J]=' : '},[te,q]=';
  const begin = source.indexOf(submitStart), end = source.indexOf(submitEnd, begin);
  if (begin < 0 || end < 0) throw new Error('Login submit handler changed');
  let body = source.slice(begin, end);
  body = body.replaceAll('b.trim()', 'dpFields.email').replace(legacy ? /\bS\b/g : /\bC\b/g, 'dpFields.password');
  const read = legacy ? 'const dpFields=dpReadCredentials(L.current,b,S);x(dpFields.email);C(dpFields.password);' : 'const dpFields=dpReadCredentials(I.current,b,C);S(dpFields.email);E(dpFields.password);';
  body = body.replace(submitStart, submitStart + read);
  source = source.slice(0, begin) + body + source.slice(end);
  if (legacy) {
    const old = '(0,K.jsx)(j,{focusKey:`auth-submit`,primary:!0,className:`auth__submit ${k?`is-busy`:``} ${U?``:`is-waiting`}`,label:k?g===`login`?`Accesso in corso…`:`Creazione in corso…`:g===`login`?`Accedi`:`Crea account`,onEnter:G})';
    source = once(source, old, 'l?' + old + ':(0,K.jsx)(`button`,{type:`submit`,className:`btn btn--primary auth__submit ${k?`is-busy`:``}`,disabled:k,"aria-busy":k,children:(0,K.jsx)(`span`,{className:`btn__label`,children:k?g===`login`?`Accesso in corso…`:`Creazione in corso…`:g===`login`?`Accedi`:`Crea account`})})');
    source = once(source, 'id:e,className:`field__input', 'id:e,name:e,className:`field__input');
    source = once(source, '(0,K.jsxs)(`form`,{onSubmit:', '(0,K.jsxs)(`form`,{onFocusCapture:()=>{const dpFields=dpReadCredentials(L.current,b,S);x(dpFields.email);C(dpFields.password)},onBlurCapture:()=>{const dpFields=dpReadCredentials(L.current,b,S);x(dpFields.email);C(dpFields.password)},onSubmit:');
  } else {
    const old = '(0,U.jsx)(O,{focusKey:`auth-submit`,primary:!0,className:`auth__submit ${j?`is-busy`:``} ${G?``:`is-waiting`}`,label:j?v===`login`?`Accesso in corso…`:`Creazione in corso…`:v===`login`?`Accedi`:`Crea account`,onEnter:K})';
    source = once(source, old, 'o?' + old + ':(0,U.jsx)(`button`,{type:`submit`,className:`btn btn--primary auth__submit ${j?`is-busy`:``}`,disabled:j,"aria-busy":j,children:(0,U.jsx)(`span`,{className:`btn__label`,children:j?v===`login`?`Accesso in corso…`:`Creazione in corso…`:v===`login`?`Accedi`:`Crea account`})})');
    source = once(source, 'id:e,className:`field__input', 'id:e,name:e,className:`field__input');
    source = once(source, '(0,U.jsxs)(`form`,{onSubmit:', '(0,U.jsxs)(`form`,{onFocusCapture:()=>{const dpFields=dpReadCredentials(I.current,b,C);S(dpFields.email);E(dpFields.password)},onBlurCapture:()=>{const dpFields=dpReadCredentials(I.current,b,C);S(dpFields.email);E(dpFields.password)},onSubmit:');
  }
  return helper + source;
}
