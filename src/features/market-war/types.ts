export type Country = "USA" | "CHINA";
export type CompanyStatus = "public" | "private";
export type UnitType = "tank" | "aircraft" | "ship" | "drone" | "armored_vehicle" | "command_center";
export type TimeRange = "24H" | "7D" | "1M" | "1Y";

export interface CompanyMarketData {
  id: string; ticker?: string; name: string; country: Country; exchange?: string;
  status: CompanyStatus; unitType: UnitType; marketCapUsd?: number;
  estimatedPrivateValuationUsd?: number; priceChangePercent24h?: number;
  priceChangePercent7d?: number; priceChangePercent1m?: number; priceChangePercent1y?: number;
  rank?: number; sector?: string; updatedAt: string; sourceName?: string; sourceUrl?: string;
  canonicalId?: string;
}

export interface HistoricalMarketPoint { date: string; marketCapUsd: number; }
export interface MarketDataProvider {
  getCompanies(): Promise<CompanyMarketData[]>;
  getHistoricalData(companyId: string, range: TimeRange): Promise<HistoricalMarketPoint[]>;
}
