# Blog comments — Cloudflare Free

Astro remains on GitHub Pages. The Worker only serves the comments API at `https://comments.donatomolino.it/v1/comments`. Italian and English articles use the same translation key. No email, account, cookies set by our code, or raw IP storage. Public comments are plain text. All submissions are pending until approved. The reader loads the list explicitly; Turnstile only loads when interacting with the form.

## Production setup (once, from your authenticated terminal)

1. Verify **Workers Free** in Cloudflare. Create a managed Turnstile widget with hostnames `donatomolino.it` and `www.donatomolino.it`; copy its public site key and keep the secret private.
2. From this repository, in `workers/comments`:

```sh
npx wrangler login
npx wrangler d1 create donatomolino-comments
cp wrangler.example.toml wrangler.toml
```

Edit the local `wrangler.toml` with the returned database ID and public Turnstile site key. Keep the Free plan; do not add paid products. The local config is gitignored. The custom domain requires a Cloudflare-managed DNS zone and no conflicting DNS record.

```sh
npx wrangler d1 execute donatomolino-comments --remote --file=schema.sql
npx wrangler secret put TURNSTILE_SECRET
npx wrangler secret put RATE_SECRET
npx wrangler deploy
```

Paste the Turnstile secret at the first prompt, and a separate randomly generated secret (e.g. `openssl rand -hex 32`) at the second. Never commit these values or paste them into chat. If Wrangler offers to create a paid resource, stop. If prompted for custom domain DNS changes, confirm only the dedicated comments subdomain.

```sh
curl 'https://comments.donatomolino.it/v1/comments?thread=future-tech'
```

Expect JSON containing an empty comments array and the public site key. Once healthy, merge the frontend PR and test a real submission from the published blog; verify it is pending, then approve it and refresh. A Worker/API test is not proof of live Turnstile or DNS.

## Moderate

In Cloudflare D1 Studio, run:

```sql
SELECT id, thread, name, body, lang, created_at
FROM comments WHERE status='pending' ORDER BY created_at ASC;
UPDATE comments SET status='approved' WHERE id='PASTE_REVIEWED_ID' AND status='pending';
-- Or reject:
UPDATE comments SET status='rejected' WHERE id='PASTE_REVIEWED_ID' AND status='pending';
```

No public admin route. A rejected parent hides its replies. Do not approve impersonation, spam or sensitive personal information. For deletion requests delete replies first, then their parent, if needed. Review rejected comments periodically; keep daily counts intact until the next UTC day. If adding new articles, use lowercase slug or translationKey (letters, numbers, hyphens).

## Limits and privacy

Workers Free: 100,000 requests/day, shared with other Workers on the account. D1 Free: 5 million rows read and 100,000 rows written/day, 5 GB total account storage (individual database limits also apply). These quotas are account-wide; service may stop when limits are reached. No guarantee of availability during abuse. The application limits submissions to 200/day total and five/day per daily HMAC identifier derived from the connecting IP. The raw IP is not written to the database, and a scheduled cleanup removes old pseudonymous rate keys at 00:15 UTC (subject to successful execution). Cloudflare still processes network/challenge data. See the privacy disclosure in the form; add this processing to the site's privacy notice if necessary. The list returns the latest 30 approved top-level comments and up to 100 approved replies. No external paid services or email notifications are configured.

Tests use local SQLite and mocked Turnstile, not a live Cloudflare deployment:

```sh
node --test workers/comments/index.test.mjs
```
