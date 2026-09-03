import type { CompanyMarketData, Country, TimeRange } from "./types";

export const deduplicateCompanies = (items: CompanyMarketData[]) => {
  const seen = new Map<string, CompanyMarketData>();
  for (const company of items) {
    const key = company.canonicalId ?? company.id;
    if (!seen.has(key) || (seen.get(key)?.status === "private" && company.status === "public")) seen.set(key, company);
  }
  return [...seen.values()];
};
export const publicCompanies = (items: CompanyMarketData[]) => deduplicateCompanies(items).filter(c => c.status === "public" && c.marketCapUsd !== undefined);
export const selectUniverse = (items: CompanyMarketData[], count: number | "all") => {
  const ranked = publicCompanies(items).sort((a,b) => (b.marketCapUsd ?? 0) - (a.marketCapUsd ?? 0));
  return count === "all" ? ranked : ranked.slice(0, count);
};
export const countryTotal = (items: CompanyMarketData[], country: Country) => publicCompanies(items).filter(c=>c.country===country).reduce((sum,c)=>sum+(c.marketCapUsd ?? 0),0);
export const dominance = (items: CompanyMarketData[]) => { const usa=countryTotal(items,"USA"), china=countryTotal(items,"CHINA"), total=usa+china; return { usa, china, usaPercent: total ? usa/total*100 : 50, chinaPercent: total ? china/total*100 : 50, total }; };
export const performance = (c: CompanyMarketData, range: TimeRange) => ({"24H":c.priceChangePercent24h,"7D":c.priceChangePercent7d,"1M":c.priceChangePercent1m,"1Y":c.priceChangePercent1y}[range] ?? 0);
export const formatMarketCap = (value?: number) => value === undefined ? "Unranked" : value >= 1e12 ? `$${(value/1e12).toFixed(2)}T` : `$${(value/1e9).toFixed(value>=1e11?0:1)}B`;
export const unitScale = (cap=0, reference=3.5e12) => Math.max(.65, Math.min(1.65, .55 + Math.sqrt(cap/reference)*1.1));
export const territoryX = (usaPercent:number, width=24) => -width/2 + width*Math.max(0,Math.min(100,usaPercent))/100;
