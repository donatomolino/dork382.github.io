# donatomolino.it

Personal Astro website, including **USA vs China — Market Cap Frontline** at `/usavschina`.

## Market Cap Frontline

An original low-poly financial diorama built with React, TypeScript, React Three Fiber, Drei, Zustand, TanStack Query, Zod, Recharts, and Framer Motion. The position of the gold frontline represents the market-cap share of the selected public-company universe. Unit scale uses square-root normalization and movement represents selected-period performance.

This is a symbolic, non-violent visualization. It compares the supported company universe—not every listed company in either country—and is not financial advice.

## Local development

```bash
npm install
npm run dev
```

Open `http://localhost:4321/usavschina`. Other scripts:

```bash
npm run check
npm test
npm run build
npm run preview
```

## Market data configuration

Mock mode is the zero-configuration default and is clearly labeled as demo data. It uses deterministic values and never claims to be live.

```env
PUBLIC_MARKET_DATA_PROVIDER=mock
MARKET_DATA_API_URL=
MARKET_DATA_API_KEY=
```

Set `PUBLIC_MARKET_DATA_PROVIDER=remote` to use the remote adapter. The browser calls `/api/market-data/*`; the Astro server proxy adds `MARKET_DATA_API_KEY`, so credentials are not shipped to clients. The remote response must conform to the normalized `CompanyMarketData` schema. Requests time out after eight seconds, retry twice through TanStack Query, are cached for 55 seconds, respect HTTP 429 errors, and fall back gracefully to mock data.

For static-only hosting, keep mock mode or deploy the API route on an Astro server adapter/serverless platform. No proprietary assets are required.

## Architecture

- `src/features/market-war`: normalized types, Zod schemas, deterministic dataset, providers, and pure financial calculations.
- `src/components/market-war/MarketWarScene.tsx`: procedural terrain, units, effects, camera, labels, and animated frontline.
- `src/components/market-war/MarketWarApp.tsx`: query lifecycle, dashboard UI, filters, settings, rankings, chart, timeline, fallback, and accessible data table.
- `src/stores/marketWarStore.ts`: persisted user preferences and interaction state, separate from raw data.
- `src/pages/api/market-data/[...path].ts`: server-side credential proxy.
- `src/features/market-war/calculations.test.ts`: financial, deduplication, formatting, schema, territory, and unit-scale tests.

Rendering never owns financial calculations, and raw provider data is validated before entering the application. Chinese ADR/HK dual listings use `canonicalId`; only one canonical economic entity is retained, preferring a public record over a private one. Private-company estimates are stored separately as `estimatedPrivateValuationUsd` and never enter public totals.

## Performance and accessibility

The 3D scene is route-local and lazy-loaded, uses shared procedural geometry, constrained orbit controls, adaptive DPR, limited shadows, memoized placement, and no React state updates in frame loops. Mobile keeps WebGL with a simplified overlay layout. If WebGL is unavailable, rankings and the keyboard-accessible data table remain usable. The interface includes semantic labels, visible focus states, non-color country labels, reduced-motion support, a screen-reader market summary, and persistent settings.

## Deployment and production data

Run `npm run build`, then start the standalone Node output with `node dist/server/entry.mjs`. Set the two server environment variables in the deployment platform for remote mode. Before declaring real-time status in production, connect an appropriately licensed market-data feed, provide exchange-calendar metadata, historical snapshot coverage, and source timestamps. Private valuations require dated, attributable sources and should be removed when reliability cannot be established.
