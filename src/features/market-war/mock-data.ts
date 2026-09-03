import type { CompanyMarketData } from "./types";

const now = "2026-07-16T16:00:00.000Z";
const p = (id: string, ticker: string, name: string, country: "USA"|"CHINA", cap: number, changes: number[], sector: string, unitType: CompanyMarketData["unitType"], exchange: string): CompanyMarketData => ({
  id, ticker, name, country, exchange, status: "public", unitType, marketCapUsd: cap * 1e9,
  priceChangePercent24h: changes[0], priceChangePercent7d: changes[1], priceChangePercent1m: changes[2], priceChangePercent1y: changes[3],
  sector, updatedAt: now, sourceName: "Deterministic demo dataset", sourceUrl: "https://example.com/demo-market-data", canonicalId: id,
});

export const mockCompanies: CompanyMarketData[] = [
  p("apple","AAPL","Apple","USA",3410,[1.2,2.4,6.8,18.4],"Technology","aircraft","NASDAQ"),
  p("nvidia","NVDA","Nvidia","USA",3320,[2.8,5.1,11.6,47.2],"Semiconductors","drone","NASDAQ"),
  p("microsoft","MSFT","Microsoft","USA",3010,[-.4,1.1,4.2,14.8],"Technology","command_center","NASDAQ"),
  p("alphabet","GOOGL","Alphabet","USA",2180,[.8,-1.2,2.1,11.6],"Communication","tank","NASDAQ"),
  p("amazon","AMZN","Amazon","USA",2110,[-1.1,.6,5.5,22.3],"Consumer","ship","NASDAQ"),
  p("meta","META","Meta","USA",1740,[1.6,3.8,8.2,29.4],"Communication","aircraft","NASDAQ"),
  p("broadcom","AVGO","Broadcom","USA",1380,[.3,2.8,9.4,38.6],"Semiconductors","drone","NASDAQ"),
  p("berkshire","BRK.B","Berkshire Hathaway","USA",1050,[-.2,.4,1.8,9.1],"Financials","command_center","NYSE"),
  p("tesla","TSLA","Tesla","USA",1010,[-2.4,-4.1,7.2,-6.8],"Automotive","armored_vehicle","NASDAQ"),
  p("eli-lilly","LLY","Eli Lilly","USA",790,[.9,1.9,-2.3,15.2],"Healthcare","armored_vehicle","NYSE"),
  p("tencent","0700.HK","Tencent","CHINA",610,[1.7,3.1,8.4,24.2],"Communication","tank","HKEX"),
  p("alibaba","9988.HK","Alibaba","CHINA",365,[-.8,2.2,10.1,18.6],"Consumer","ship","HKEX"),
  p("xiaomi","1810.HK","Xiaomi","CHINA",190,[2.3,5.8,13.2,42.1],"Technology","aircraft","HKEX"),
  p("china-mobile","0941.HK","China Mobile","CHINA",185,[.2,1.4,3.9,12.2],"Telecom","command_center","HKEX"),
  p("pdd","PDD","PDD Holdings","CHINA",160,[-1.7,-3.2,4.6,-8.1],"Consumer","armored_vehicle","NASDAQ"),
  p("byd","1211.HK","BYD","CHINA",145,[1.4,2.7,6.1,28.3],"Automotive","armored_vehicle","HKEX"),
  p("netease","9999.HK","NetEase","CHINA",92,[-.4,.8,2.5,6.9],"Communication","drone","HKEX"),
  p("meituan","3690.HK","Meituan","CHINA",88,[.7,-2.1,-5.2,-12.4],"Consumer","armored_vehicle","HKEX"),
  p("jd","9618.HK","JD.com","CHINA",61,[-1.2,.3,3.1,4.8],"Consumer","ship","HKEX"),
  p("baidu","9888.HK","Baidu","CHINA",33,[1.1,2.6,-1.7,-4.2],"Technology","drone","HKEX"),
  { id:"deepseek", name:"DeepSeek", country:"CHINA", status:"private", unitType:"drone", sector:"Artificial Intelligence", estimatedPrivateValuationUsd:8e9, updatedAt:now, sourceName:"Illustrative estimate — demo only", sourceUrl:"https://example.com/private-valuation", canonicalId:"deepseek" },
  { id:"moonshot", name:"Moonshot AI", country:"CHINA", status:"private", unitType:"ship", sector:"Artificial Intelligence", updatedAt:now, sourceName:"Unranked — no reliable public valuation", canonicalId:"moonshot" },
];
