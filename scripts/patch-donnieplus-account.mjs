// Remove the original server-key features rather than retaining dead key UI.
function replaceFunction(source, start, end, replacement) {
  const a = source.indexOf(start), b = source.indexOf(end, a + start.length);
  if (a < 0 || b < 0) throw new Error('Account-only patch no longer matches upstream: ' + start);
  return source.slice(0, a) + replacement + source.slice(b);
}
export function patchAccountOnly(source, name, main, legacyMain) {
  if (name === main) {
    source = replaceFunction(source, 'async function Ic(', 'function Lc()', 'async function Ic(e,t=`leotv`){if(t===`leomusic`)return{link:`${Pc}#invito=${encodeURIComponent(e)}`,direct:!0};return{link:location.origin+`/donnieplus/?invito=${encodeURIComponent(e)}`,direct:!1}}');
    source = replaceFunction(source, 'function Fc(){try{let e=location.hash', 'async function Ic(', 'function Fc(){}');
  } else if (name === legacyMain) {
    source = replaceFunction(source, 'async function Bc(', 'function $()', 'async function Bc(e,t=`leotv`){if(t===`leomusic`)return{link:`${qi}#invito=${encodeURIComponent(e)}`,direct:!0};return{link:location.origin+`/donnieplus/?invito=${encodeURIComponent(e)}`,direct:!1}}');
    source = replaceFunction(source, 'function zc(){try{let e=location.hash', 'async function Bc(', 'function zc(){}');
  } else if (/^Library-legacy-/.test(name)) {
    source = replaceFunction(source, 'function Ge({onFocusNode:e,brief:t})', 'function Ke({onFocusNode:e})', 'function Ge({onFocusNode:e}){return(0,Z.jsx)(J,{focusKey:`set-srv`,title:`Collegamento`,footer:`Accesso tramite account DonniePlus.`,children:(0,Z.jsx)(q,{k:`set-srv-address`,icon:Y.link,tint:K.teal,title:`Server`,value:(F()||location.origin).replace(/\\/+$/,``),onFocusNode:e})})}');
  } else if (/^Library-/.test(name)) {
    source = replaceFunction(source, 'function ze({onFocusNode:t,brief:n})', 'function Be({onFocusNode:e})', 'function ze({onFocusNode:t}){return(0,Z.jsx)(J,{focusKey:`set-srv`,title:`Collegamento`,footer:`Accesso tramite account DonniePlus.`,children:(0,Z.jsx)(q,{k:`set-srv-address`,icon:Y.link,tint:G.teal,title:`Server`,value:(P()||location.origin).replace(/\\/+$/,``),onFocusNode:t})})}');
  }
  return source;
}
