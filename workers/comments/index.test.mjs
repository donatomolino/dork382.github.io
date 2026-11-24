import test from "node:test";
import assert from "node:assert/strict";
import {DatabaseSync} from "node:sqlite";
import {readFileSync} from "node:fs";
import worker from "./index.mjs";

test("comments are moderated, validated, rate limited and scoped", async () => {
  const db = new DatabaseSync(":memory:");
  db.exec(readFileSync(new URL("./schema.sql",import.meta.url),"utf8"));
  const DB = {prepare(sql) {return {bind(...values) {
    const q=db.prepare(sql);
    return {all:async()=>({results:q.all(...values)}),first:async()=>q.get(...values)??null,run:async()=>({meta:{changes:Number(q.run(...values).changes)}})};
  }};}};
  const env={DB,TURNSTILE_SITE_KEY:"public",TURNSTILE_SECRET:"secret",RATE_SECRET:"rate-secret"};
  const realFetch=globalThis.fetch;
  let valid=true,hostname="donatomolino.it",action="blog_comment";
  globalThis.fetch=async()=>Response.json({success:valid,hostname,action});
  const url="https://comments.donatomolino.it/v1/comments";
  const origin="https://donatomolino.it";
  const fields={thread:"future-tech",lang:"it",name:"Reader",body:"<script>plain text</script>",token:"test",website:""};
  const post=(patch={},headers={})=>worker.fetch(new Request(url,{method:"POST",headers:{"Content-Type":"application/json","Origin":origin,"CF-Connecting-IP":"192.0.2.1",...headers},body:JSON.stringify({...fields,...patch})}),env);
  const get=async(thread="future-tech")=>(await worker.fetch(new Request(url+"?thread="+thread),env)).json();
  try{
    assert.equal((await post({},{"Origin":"https://evil.test"})).status,403);
    assert.equal((await post({website:"bot"})).status,400);
    assert.equal((await post({body:"x".repeat(2001)})).status,400);
    valid=false;assert.equal((await post()).status,403);valid=true;
    hostname="evil.test";assert.equal((await post()).status,403);hostname="donatomolino.it";
    action="other";assert.equal((await post()).status,403);action="blog_comment";
    assert.equal((await post()).status,201);
    assert.equal((await get()).comments.length,0);
    const id=db.prepare("SELECT id FROM comments").get().id;
    assert.equal((await post({parentId:id})).status,400);
    db.prepare("UPDATE comments SET status='approved' WHERE id=?").run(id);
    const shown=(await get()).comments;
    assert.equal(shown.length,1);assert.equal(shown[0].body,fields.body);assert.equal(shown[0].rate_key,undefined);
    assert.equal((await get("another-thread")).comments.length,0);
    assert.equal((await post({parentId:id,thread:"another-thread"})).status,400);
    assert.equal((await post({parentId:id,lang:"en"})).status,201);
    db.exec("UPDATE comments SET status='approved'");
    assert.equal((await get()).comments.length,2);
    for(let i=0;i<3;i++)assert.equal((await post()).status,201);
    assert.equal((await post()).status,429);
    assert.equal((await post({},{"CF-Connecting-IP":"192.0.2.2"})).status,201);
    const preflight=await worker.fetch(new Request(url,{method:"OPTIONS",headers:{Origin:origin}}),env);
    assert.equal(preflight.status,204);assert.equal(preflight.headers.get("Access-Control-Allow-Origin"),origin);
    assert.equal((await worker.fetch(new Request(url+"?thread=future-tech"),{...env,RATE_SECRET:undefined})).status,503);
    db.exec("UPDATE comments SET day='2020-01-01'");
    let cleanup;
    await worker.scheduled(null,env,{waitUntil(p){cleanup=p;}});
    await cleanup;
    assert.equal(db.prepare("SELECT COUNT(*) AS count FROM comments WHERE rate_key IS NOT NULL").get().count,0);
  }finally{globalThis.fetch=realFetch;db.close();}
});
