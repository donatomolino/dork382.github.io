export interface TelegramPost {
  id: number
  date: string
  text?: string
  imageUrl?: string
  videoThumb?: string
  videoUrl?: string
}

function formatDate(isoDate: string): string {
  return new Date(isoDate).toLocaleDateString("it-IT", {
    year: "numeric",
    month: "long",
    day: "numeric",
    hour: "2-digit",
    minute: "2-digit",
    timeZone: "Europe/Rome",
  })
}

function decode(text: string): string {
  return text
    .replace(/<br\s*\/?>/gi, "\n")
    .replace(/<[^>]+>/g, "")
    .replace(/&amp;/g, "&")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'")
    .replace(/&nbsp;/g, " ")
    .trim()
}

export function parsePosts(html: string): TelegramPost[] {
  const posts: TelegramPost[] = []
  const messageRegex =
    /data-post="[^"]+\/(\d+)"[^>]*>([\s\S]*?)(?=<div class="tgme_widget_message_wrap|$)/g

  let match: RegExpExecArray | null
  while ((match = messageRegex.exec(html)) !== null) {
    const [, idStr, body] = match
    if (body.includes("service_message")) continue

    const id = parseInt(idStr, 10)

    const textMatch = body.match(
      /<div class="tgme_widget_message_text[^"]*"[^>]*dir="auto">([\s\S]*?)<\/div>/,
    )
    const text = textMatch ? decode(textMatch[1]) : undefined

    let imageUrl: string | undefined
    let videoThumb: string | undefined
    let videoUrl: string | undefined

    const bgMatch = body.match(
      /style="background-image:\s*url\(['"]?([^'"')]+)['"]?\)/,
    )
    if (bgMatch) {
      const url = bgMatch[1].replace(/&amp;/g, "&").replace(/&#39;/g, "'")
      if (/photo_wrap|image/.test(bgMatch.input?.substring(0, bgMatch.index) || "")) {
        imageUrl = url
      } else if (/video_wrap|video|animation|roundvideo/.test(bgMatch.input?.substring(0, bgMatch.index) || "")) {
        videoThumb = url
      }
    }

    const srcMatch = body.match(/<video\s[^>]*src="([^"]+)"/)
    if (srcMatch) {
      videoUrl = srcMatch[1]
    }

    const dateMatch = body.match(/<time[^>]*datetime="([^"]+)"[^>]*>/)
    const date = dateMatch ? formatDate(dateMatch[1]) : ""

    posts.push({ id, date, text, imageUrl, videoThumb, videoUrl })
  }
  return posts.reverse()
}
