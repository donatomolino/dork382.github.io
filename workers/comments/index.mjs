const origins = new Set(["https://donatomolino.it", "https://www.donatomolino.it"]);
const hostnames = new Set(["donatomolino.it", "www.donatomolino.it"]);
const MAX_BODY = 16384;
const threadPattern = /^[a-z0-9][a-z0-9-]{0,95}$/;
const idPattern = /^[0-9a-f-]{36}$/i;
const columns = "id, parent_id, name, body, lang, created_at";
function response(request, data, status = 200) {
  const headers = {"Content-Type":"application/json; charset=utf-8", "Cache-Control":"no-store", "X-Content-Type-Options":"nosniff", "Vary":"Origin"};
  const origin = request.headers.get("Origin");
  if (origins.has(origin)) {
    headers["Access-Control-Allow-Origin"] = origin;
    headers["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS";
    headers["Access-Control-Allow-Headers"] = "Content-Type";
  }
  return new Response(status === 204 ? null : JSON.stringify(data), {status, headers});
}
async function readBody(request) {
  if (!request.body || !request.headers.get("Content-Type")?.startsWith("application/json")) throw new Error("invalid");
  const reader = request.body.getReader();
  const decoder = new TextDecoder();
  let size = 0, body = "";
  try {
    for (;;) {
      const {done,value} = await reader.read();
      if (done) break;
      size += value.byteLength;
      if (size > MAX_BODY) { await reader.cancel(); throw new Error("large"); }
      body += decoder.decode(value,{stream:true});
    }
    return JSON.parse(body + decoder.decode());
  } finally { reader.releaseLock(); }
}
async function dailyKey(secret, day, ip) {
  const enc = new TextEncoder();
  const key = await crypto.subtle.importKey("raw",enc.encode(secret),{name:"HMAC",hash:"SHA-256"},false,["sign"]);
  const bytes = await crypto.subtle.sign("HMAC",key,enc.encode("comment-rate:"+day+":"+ip));
  return Array.from(new Uint8Array(bytes),b=>b.toString(16).padStart(2,"0")).join("");
}
export default {
  async fetch(request, env) {
    const reply = (data,status) => response(request,data,status);
    const url = new URL(request.url);
    if (url.pathname !== "/v1/comments") return reply({error:"not_found"},404);
    const origin = request.headers.get("Origin");
    if (origin && !origins.has(origin)) return reply({error:"origin"},403);
    if (request.method === "OPTIONS") return origins.has(origin) ? reply(null,204) : reply({error:"origin"},403);
    if (!["GET","POST"].includes(request.method)) return reply({error:"method"},405);
    if (!env.DB || !env.TURNSTILE_SITE_KEY || !env.TURNSTILE_SECRET || !env.RATE_SECRET) return reply({error:"unavailable"},503);
    try {
      if (request.method === "GET") {
        const thread = url.searchParams.get("thread") || "";
        if (!threadPattern.test(thread)) return reply({error:"invalid_thread"},400);
        const roots = await env.DB.prepare("SELECT "+columns+" FROM comments WHERE thread=? AND status='approved' AND parent_id IS NULL ORDER BY created_at DESC,id DESC LIMIT 30").bind(thread).all();
        const replies = await env.DB.prepare("SELECT "+columns+" FROM comments WHERE thread=? AND status='approved' AND parent_id IN (SELECT id FROM comments WHERE thread=? AND status='approved' AND parent_id IS NULL ORDER BY created_at DESC,id DESC LIMIT 30) ORDER BY created_at ASC,id ASC LIMIT 100").bind(thread,thread).all();
        return reply({comments:[...roots.results,...replies.results],siteKey:env.TURNSTILE_SITE_KEY});
      }
      if (!origins.has(origin)) return reply({error:"origin"},403);
      let data;
      try { data = await readBody(request); } catch { return reply({error:"invalid_body"},400); }
      if (!data || typeof data !== "object" || Array.isArray(data)) return reply({error:"invalid_body"},400);
      const {thread,lang,token} = data;
      const name = typeof data.name === "string" ? data.name.trim() : "";
      const body = typeof data.body === "string" ? data.body.trim() : "";
      const parent = data.parentId ?? null;
      if (typeof thread !== "string" || !threadPattern.test(thread) || !["it","en"].includes(lang) || !name || name.length>60 || !body || body.length>2000 || typeof token !== "string" || token.length<1 || token.length>2048 || (parent!==null && (typeof parent!=="string" || !idPattern.test(parent))) || data.website) return reply({error:"invalid_fields"},400);
      const verify = await fetch("https://challenges.cloudflare.com/turnstile/v0/siteverify",{
        method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({secret:env.TURNSTILE_SECRET,response:token})
      });
      const result = await verify.json();
      if (!verify.ok || result.success!==true || !hostnames.has(result.hostname) || result.action!=="blog_comment") return reply({error:"verification"},403);
      const now = Date.now(), day = new Date(now).toISOString().slice(0,10);
      const ip = request.headers.get("CF-Connecting-IP");
      if (!ip) return reply({error:"unavailable"},503);
      const rateKey = await dailyKey(env.RATE_SECRET,day,ip);
      if (parent) {
        const root = await env.DB.prepare("SELECT id FROM comments WHERE id=? AND thread=? AND status='approved' AND parent_id IS NULL").bind(parent,thread).first();
        if (!root) return reply({error:"invalid_parent"},400);
      }
      // Atomic conditional insert: no race between checking the daily cap and writing.
      const write = await env.DB.prepare(`INSERT INTO comments (id,thread,parent_id,name,body,lang,status,created_at,day,rate_key)
        SELECT ?,?,?,?,?,?,'pending',?,?,?
        WHERE (SELECT COUNT(*) FROM comments WHERE day=?) < 200
        AND (SELECT COUNT(*) FROM comments WHERE day=? AND rate_key=?) < 5`)
        .bind(crypto.randomUUID(),thread,parent,name,body,lang,now,day,rateKey,day,day,rateKey).run();
      if (!write.meta.changes) return reply({error:"rate_limit"},429);
      return reply({status:"pending"},201);
    } catch {
      // Never leak credentials, submitted content or SQL diagnostics.
      return reply({error:"unavailable"},503);
    }
  },
  async scheduled(_event,env,ctx) {
    const day = new Date().toISOString().slice(0,10);
    ctx.waitUntil(env.DB.prepare("UPDATE comments SET rate_key=NULL WHERE day<? AND rate_key IS NOT NULL").bind(day).run());
  }
};
