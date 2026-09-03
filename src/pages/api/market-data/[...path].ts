import type { APIRoute } from "astro";

export const prerender=false;
export const GET:APIRoute=async({params,request})=>{
  const base=import.meta.env.MARKET_DATA_API_URL;
  const key=import.meta.env.MARKET_DATA_API_KEY;
  if(!base)return new Response(JSON.stringify({error:"Remote provider is not configured"}),{status:503,headers:{"content-type":"application/json"}});
  const url=new URL(params.path??"",base.endsWith("/")?base:`${base}/`);url.search=new URL(request.url).search;
  const controller=new AbortController();const timer=setTimeout(()=>controller.abort(),8000);
  try{const upstream=await fetch(url,{headers:{Accept:"application/json",...(key?{Authorization:`Bearer ${key}`}:{})},signal:controller.signal});return new Response(upstream.body,{status:upstream.status,headers:{"content-type":"application/json","cache-control":"s-maxage=55, stale-while-revalidate=240"}})}catch(error){return new Response(JSON.stringify({error:error instanceof Error?error.message:"Provider unavailable"}),{status:502,headers:{"content-type":"application/json"}})}finally{clearTimeout(timer)}
};
