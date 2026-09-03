import { companiesSchema, historySchema } from "./schema";
import { mockCompanies } from "./mock-data";
import type { HistoricalMarketPoint, MarketDataProvider, TimeRange } from "./types";

const delay = (ms:number) => new Promise(resolve=>setTimeout(resolve,ms));
export class MockMarketDataProvider implements MarketDataProvider {
  async getCompanies(){ await delay(180); return companiesSchema.parse(mockCompanies); }
  async getHistoricalData(companyId:string, range:TimeRange):Promise<HistoricalMarketPoint[]>{
    const company=mockCompanies.find(c=>c.id===companyId); if(!company?.marketCapUsd) return [];
    const points=range==="24H"?24:range==="7D"?28:range==="1M"?30:52; const swing=(company.priceChangePercent1y??4)/100;
    return historySchema.parse(Array.from({length:points},(_,i)=>({date:new Date(Date.now()-(points-i)*86400000).toISOString(),marketCapUsd:company.marketCapUsd!*(1-swing+swing*i/points)*(1+Math.sin(i*.8)*.012)})));
  }
}
export class RemoteMarketDataProvider implements MarketDataProvider {
  async request(path:string){ const controller=new AbortController(); const timer=setTimeout(()=>controller.abort(),8000); try { const response=await fetch(`/api/market-data${path}`,{signal:controller.signal}); if(response.status===429) throw new Error("Rate limited"); if(!response.ok) throw new Error(`Provider error ${response.status}`); return await response.json(); } finally { clearTimeout(timer); } }
  async getCompanies(){ return companiesSchema.parse(await this.request("/companies")); }
  async getHistoricalData(companyId:string,range:TimeRange){ return historySchema.parse(await this.request(`/history/${encodeURIComponent(companyId)}?range=${range}`)); }
}
export const getProvider = () => import.meta.env.PUBLIC_MARKET_DATA_PROVIDER === "remote" ? new RemoteMarketDataProvider() : new MockMarketDataProvider();
