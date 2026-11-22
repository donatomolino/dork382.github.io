import type { APIRoute, GetStaticPaths } from "astro";
import { getCollection } from "astro:content";
import sharp from "sharp";

export const prerender = true;

const escapeXml = (text: string) =>
  text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&apos;");

function splitTitle(title: string): string[] {
  const words = title.split(/\s+/);
  const lines: string[] = [];
  let current = "";
  for (const word of words) {
    const next = current ? current + " " + word : word;
    if (next.length > 25 && current) {
      lines.push(current);
      current = word;
    } else {
      current = next;
    }
  }
  if (current) lines.push(current);
  return lines.slice(0, 3);
}

export const getStaticPaths = (async () => {
  const posts = await getCollection("blog", ({ data }) =>
    !data.draft && data.pubDate <= new Date()
  );
  return posts.map((post) => ({
    params: { lang: post.data.lang, slug: post.data.slug },
    props: { title: post.data.title, lang: post.data.lang },
  }));
}) satisfies GetStaticPaths;

export const GET: APIRoute = async ({ props }) => {
  const title = String(props.title);
  const lang = String(props.lang);
  const lines = splitTitle(title);
  const fontSize = lines.length > 2 ? 69 : 83;
  const gap = fontSize * 1.12;
  const startY = lines.length > 2 ? 234 : 277;
  const text = lines.map((line, index) => (
    '<text x="72" y="' + (startY + index * gap) +
    '" font-family="Arial,Helvetica,sans-serif" font-size="' + fontSize +
    '" letter-spacing="-2" font-weight="800" fill="#11110f">' +
    escapeXml(line) + "</text>"
  )).join("");

  const svg = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="630" viewBox="0 0 1200 630">',
    '<rect width="1200" height="630" fill="#f2f1ed"/>',
    '<rect x="72" y="72" width="22" height="22" fill="#ed1c24"/>',
    '<text x="112" y="91" font-family="Arial,Helvetica,sans-serif" font-size="23" font-weight="800" fill="#11110f">DONATO MOLINO</text>',
    '<text x="1128" y="92" font-family="Arial,Helvetica,sans-serif" font-size="19" font-weight="700" fill="#11110f" text-anchor="end">' + escapeXml(lang.toUpperCase()) + ' / BLOG</text>',
    '<line x1="72" y1="122" x2="1128" y2="122" stroke="#11110f" stroke-width="3"/>',
    text,
    '<line x1="72" y1="539" x2="1128" y2="539" stroke="#11110f" stroke-width="3"/>',
    '<text x="72" y="583" font-family="Arial,Helvetica,sans-serif" font-size="23" font-weight="700" fill="#11110f">PRODUCT ENGINEERING / AI / SOFTWARE</text>',
    '<text x="1128" y="583" font-family="Arial,Helvetica,sans-serif" font-size="23" font-weight="800" fill="#ed1c24" text-anchor="end">DONATOMOLINO.IT</text>',
    '</svg>',
  ].join("");

  const png = await sharp(Buffer.from(svg)).png({ compressionLevel: 9 }).toBuffer();
  return new Response(new Uint8Array(png), {
    headers: { "Content-Type": "image/png", "Cache-Control": "public, max-age=31536000" },
  });
};
