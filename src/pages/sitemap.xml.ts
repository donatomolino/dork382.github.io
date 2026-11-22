import type { APIRoute } from "astro";
import { getCollection } from "astro:content";

export const prerender = true;

const SITE = "https://donatomolino.it";
const toUrl = (path: string) => new URL(path, SITE).toString();
const escapeXml = (value: string) =>
  value.replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

type Entry = {
  loc: string;
  lastmod?: string;
  alternates?: { lang: "it" | "en"; url: string }[];
};

function staticRoutes(): string[] {
  // Vite enumerates Astro pages at build time; no manual sitemap maintenance needed
  // when a new ordinary .astro route is added. Parametric blog paths are added below.
  const files = Object.keys(import.meta.glob("./**/*.astro"));
  return files
    .filter((file) =>
      !file.includes("[") &&
      !/(^|\/)(404|500|_.*)\.astro$/.test(file) &&
      !file.startsWith("./api/") &&
      !file.startsWith("./admin/") &&
      !file.startsWith("./dev/") &&
      !file.startsWith("./test/") &&
      file !== "./blog/index.astro" // language detector is noindex
    )
    .map((file) => {
      const name = file.slice(2, -".astro".length);
      if (name === "index") return "/";
      const route = name.endsWith("/index")
        ? name.slice(0, -"/index".length)
        : name;
      return "/" + route + "/";
    });
}

function generateEntry({ loc, lastmod, alternates }: Entry): string {
  const languages = (alternates ?? [])
    .map(({ lang, url }) =>
      '<xhtml:link rel="alternate" hreflang="' + lang + '" href="' + escapeXml(url) + '"/>'
    )
    .join("");
  return "<url><loc>" + escapeXml(loc) + "</loc>" +
    (lastmod ? "<lastmod>" + lastmod + "</lastmod>" : "") +
    languages + "</url>";
}

export const GET: APIRoute = async () => {
  const posts = await getCollection("blog", ({ data }) =>
    !data.draft && data.pubDate <= new Date()
  );

  const indexes = [
    { lang: "it" as const, url: toUrl("/blog/it/") },
    { lang: "en" as const, url: toUrl("/blog/en/") },
  ];

  const byTranslationKey = new Map<string, typeof posts>();
  for (const post of posts) {
    if (!post.data.translationKey) continue;
    const siblings = byTranslationKey.get(post.data.translationKey) ?? [];
    siblings.push(post);
    byTranslationKey.set(post.data.translationKey, siblings);
  }

  const entries = new Map<string, Entry>();
  for (const path of staticRoutes()) entries.set(toUrl(path), { loc: toUrl(path) });

  for (const { url } of indexes) {
    entries.set(url, { loc: url, alternates: indexes });
  }

  for (const post of posts) {
    const loc = toUrl("/blog/" + post.data.lang + "/" + post.data.slug + "/");
    const alternates = post.data.translationKey
      ? (byTranslationKey.get(post.data.translationKey) ?? [])
          .map((sibling) => ({
            lang: sibling.data.lang,
            url: toUrl("/blog/" + sibling.data.lang + "/" + sibling.data.slug + "/"),
          }))
      : [];
    entries.set(loc, {
      loc,
      lastmod: (post.data.updatedDate ?? post.data.pubDate).toISOString().slice(0, 10),
      alternates: alternates.length > 1 ? alternates : undefined,
    });
  }

  const xml = '<?xml version="1.0" encoding="UTF-8"?>\n' +
    '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9" ' +
    'xmlns:xhtml="http://www.w3.org/1999/xhtml">\n' +
    Array.from(entries.values())
      .sort((a, b) => a.loc.localeCompare(b.loc))
      .map(generateEntry)
      .join("\n") +
    '\n</urlset>';

  return new Response(xml, {
    headers: {
      "Content-Type": "application/xml; charset=utf-8",
      "Cache-Control": "public, max-age=3600",
    },
  });
};
