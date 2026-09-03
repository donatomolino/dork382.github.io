import { describe,expect,it } from "vitest";
import { companiesSchema } from "./schema";
import { countryTotal,deduplicateCompanies,dominance,formatMarketCap,selectUniverse,territoryX,unitScale } from "./calculations";
import { mockCompanies } from "./mock-data";

describe("market-war calculations",()=>{
 it("calculates dominance from public companies only",()=>{const d=dominance(mockCompanies);expect(d.usaPercent+d.chinaPercent).toBeCloseTo(100);expect(d.usa).toBeGreaterThan(d.china)});
 it("separates private valuation from totals",()=>{const total=countryTotal(mockCompanies,"CHINA");expect(total).not.toBeGreaterThanOrEqual(mockCompanies.filter(c=>c.country==="CHINA").reduce((s,c)=>s+(c.marketCapUsd??0)+(c.estimatedPrivateValuationUsd??0),0))});
 it("prevents duplicate canonical listings",()=>{const duplicate={...mockCompanies[10],id:"tencent-adr"};expect(deduplicateCompanies([...mockCompanies,duplicate])).toHaveLength(mockCompanies.length)});
 it("formats market caps",()=>{expect(formatMarketCap(2.15e12)).toBe("$2.15T");expect(formatMarketCap(81e9)).toBe("$81.0B")});
 it("sorts and selects a universe",()=>{const selected=selectUniverse(mockCompanies,10);expect(selected).toHaveLength(10);expect(selected[0].marketCapUsd).toBeGreaterThanOrEqual(selected[9].marketCapUsd!)});
 it("validates the dataset with Zod",()=>expect(companiesSchema.parse(mockCompanies)).toHaveLength(22));
 it("maps territory and clamps unit scale",()=>{expect(territoryX(50,24)).toBe(0);expect(territoryX(100,24)).toBe(12);expect(unitScale(1e15)).toBe(1.65)});
});
