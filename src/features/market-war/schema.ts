import { z } from "zod";

export const companySchema = z.object({
  id: z.string().min(1), ticker: z.string().optional(), name: z.string().min(1),
  country: z.enum(["USA", "CHINA"]), exchange: z.string().optional(),
  status: z.enum(["public", "private"]),
  unitType: z.enum(["tank", "aircraft", "ship", "drone", "armored_vehicle", "command_center"]),
  marketCapUsd: z.number().nonnegative().optional(), estimatedPrivateValuationUsd: z.number().nonnegative().optional(),
  priceChangePercent24h: z.number().optional(), priceChangePercent7d: z.number().optional(),
  priceChangePercent1m: z.number().optional(), priceChangePercent1y: z.number().optional(),
  rank: z.number().int().positive().optional(), sector: z.string().optional(), updatedAt: z.string().datetime(),
  sourceName: z.string().optional(), sourceUrl: z.string().url().optional(), canonicalId: z.string().optional(),
}).superRefine((value, context) => {
  if (value.status === "private" && value.marketCapUsd !== undefined) context.addIssue({ code: "custom", message: "Private companies cannot have public market cap", path: ["marketCapUsd"] });
});

export const companiesSchema = z.array(companySchema);
export const historySchema = z.array(z.object({ date: z.string(), marketCapUsd: z.number().nonnegative() }));
