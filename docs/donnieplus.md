# DonniePlus on GitHub Pages

URL: https://donatomolino.it/donnieplus/

`public/donnieplus/` contains the DonniePlus static frontend with the original LeoTV Liquid Glass CSS. Astro copies it to `dist/donnieplus/` without changing portfolio routes or requiring a server adapter. The app has a separate PWA scope and cache.

The client calls `https://macleo.tail1ed5a3.ts.net` directly. GitHub Pages does not run the Sites Worker proxy. Existing LeoTV account login and profile selection supply `X-Account-Token`, `X-Profile`, and `X-Profile-Token`. No access key or password is included. `/api/access` is not called.

## Backend prerequisite — blocking publication

At preparation time, a GET to `/api/health` with `Origin: https://donatomolino.it` returned HTTP 200 **without** `Access-Control-Allow-Origin`. The static frontend therefore cannot use this backend from the intended domain yet.

Before merging, the backend owner must:

- Permit the exact origin `https://donatomolino.it` on API, signed video and artwork responses, including errors.
- Handle browser OPTIONS preflights for the methods used by the app (GET, HEAD, POST, PUT, PATCH, DELETE, OPTIONS).
- Allow request headers `Content-Type`, `X-Account-Token`, `X-Profile`, `X-Profile-Token`, `Range`, and `If-Range`.
- Expose `Content-Range`, `Accept-Ranges`, and `Content-Length`; preserve 206 responses and streaming.
- Keep account authentication and signed media URL validation enabled. No server access key is needed.

Test `/api/health`, account login, profile selection and video playback from the real domain before calling the migration complete. Existing sessions on the ChatGPT Site cannot transfer across origins, so visitors log in again.

## Validation

Run `npm run build`. Check that `dist/donnieplus/index.html` and all referenced assets exist. The GitHub Pages workflow publishes the Astro build on pushes to `master`. This integration is a snapshot; changes to the upstream LeoTV build must be imported deliberately.
