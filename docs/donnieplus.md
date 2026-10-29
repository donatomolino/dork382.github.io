# DonniePlus on GitHub Pages

Live URL: https://donatomolino.it/donnieplus/

Astro copies `public/donnieplus/` unchanged into the static build. The portfolio routes remain independent. The existing GitHub Pages workflow builds and deploys pushes to `master`.

## Upstream updates

Current upstream: LeoTV build `20261029160502`, commit `5c251237a0d3e9032a276cd6df00f30379a75ef3` from `Leoz1000/leotv-web`. See `public/donnieplus/upstream.json`.

To import another checked-out upstream version:

```sh
node scripts/update-donnieplus.mjs /absolute/path/to/leotv-web
node --test tests/donnieplus-preferences.test.mjs
npm run build
```

The importer copies the complete modern and legacy bundles, preserves the original Liquid Glass CSS, rebrands visible text, removes access-key requests, renames bundled assets, and rebuilds the entry page and precache. Review its transformations whenever upstream changes authentication or entry-point structure. It deliberately excludes upstream APK downloads.

## Authentication and backend

The static client connects directly to `https://macleo.tail1ed5a3.ts.net`. Account login and profile selection supply `X-Account-Token`, `X-Profile` and `X-Profile-Token`. No server key or password is embedded. `/api/access` is never called and `X-LeoTv-Key` is never sent.

The backend must allow origin `https://donatomolino.it`, handle OPTIONS preflights, accept account/profile headers and Range, and expose Content-Range, Accept-Ranges and Content-Length. Preserve signed URLs, streaming and 206 responses. CORS for health and login preflight was verified on the live domain during the original deployment. Full login and video playback require an authenticated account.

## Personalization

`donnieplus-settings.js` adds the DonniePlus login wordmark and an “Aspetto e preferenze” entry in the login and settings screens. The separate, keyboard-accessible dialog offers six themes (original Liquid Glass, Oceano, Minimal, Cinema, Alba, Neon) and 23 device-local preferences covering glass effects, layout, readability, motion, content labels and video defaults.

Preferences are stored in `donnieplus.preferences.v1`, validated and clamped when imported, and synchronized across tabs. Exports contain only these preferences. Reset preserves account data. Cache clearing and service-worker activation touch only `donnieplus-pages-` caches, preserving portfolio caches and downloaded media. iOS can retain system control over video volume.

Custom styles layer over upstream CSS; the default remains Liquid Glass. Cloudflare Rocket Loader is disabled for app scripts through `data-cfasync="false"`.

## Validation

The Node tests cover malformed storage, import validation and ranges, precache completeness, account-only authentication in both bundles and local asset references. Also verify the login wordmark, theme persistence, tabs and settings controls in a browser. After a deployment, confirm the GitHub Pages workflow succeeds and check the live entry page.
