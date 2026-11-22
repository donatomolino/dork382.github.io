import type { APIRoute } from "astro";
import { getCollection } from "astro:content";

export const prerender = true;

const escapeXml = (value: string) =>
  value.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");

export const GET: APIRoute = async () => {
  const posts = (await getCollection("blog", ({ data }) =>
    !data.draft && data.pubDate <= new Date()
  )).sort((a, b) => b.data.pubDate.getTime() - a.data.pubDate.getTime());

  const items = posts.map((post) => {
    const url = "https://donatomolino.it/blog/" + post.data.lang + "/" + post.data.slug + "/";
    const { title, description, pubDate, lang } = post.data;
    return [
      "<item>",
      "<title>" + escapeXml(title) + "</title>",
      "<link>" + escapeXml(url) + "</link>",
      '<guid isPermaLink="true">' + escapeXml(url) + "</guid>",
      "<description>" + escapeXml(description) + "</description>",
      "<pubDate>" + pubDate.toUTCString() + "</pubDate>",
      '<dc:language xmlns:dc="http://purl.org/dc/elements/1.1/">' + lang + "</dc:language>",
      "</item>",
    ].join("");
  }).join("\n");

  const rss = '<?xml version="1.0" encoding="UTF-8"?>\n' +
    '<rss version="2.0"><channel>' +
    "<title>Donato Molino — Blog (EN &amp; IT)</title>" +
    "<link>https://donatomolino.it/blog/</link>" +
    "<description>Original writing on software, AI, product engineering and the future of technology.</description>" +
    "<language>en</language>" +
    items + "</channel></rss>";

  return new Response(rss, {
    headers: { "Content-Type": "application/rss+xml; charset=utf-8" },
  });
};
