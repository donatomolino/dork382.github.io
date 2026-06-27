import type { MoodId } from "./thinky-moods";

export type ThinkyReactionId =
  | "cookie"
  | "pizza"
  | "coffee"
  | "sleep"
  | "music"
  | "code"
  | "idea"
  | "love"
  | "rain"
  | "rocket"
  | "bug"
  | "paint"
  | "book"
  | "money"
  | "fire"
  | "hug"
  | "magic"
  | "default";

export type ThinkyReaction = {
  id: ThinkyReactionId;
  mood: MoodId;
  bubble: string;
  keywords: string[];
};

export const THINKY_REACTIONS: ThinkyReaction[] = [
  {
    id: "cookie",
    mood: "felice",
    bubble: "Gnam.",
    keywords: ["biscotto", "cookie", "cookies", "cracker"],
  },
  {
    id: "pizza",
    mood: "entusiasta",
    bubble: "Fetta.",
    keywords: ["pizza", "margherita", "pepperoni"],
  },
  {
    id: "coffee",
    mood: "focus",
    bubble: "Sorso.",
    keywords: ["coffee", "caffe", "caffè", "espresso", "latte"],
  },
  {
    id: "sleep",
    mood: "pensieroso",
    bubble: "Pausa.",
    keywords: ["sleep", "sleepy", "stanco", "tired", "zzz", "nap"],
  },
  {
    id: "music",
    mood: "entusiasta",
    bubble: "Ballo.",
    keywords: ["music", "song", "canzone", "musica", "dance", "beat"],
  },
  {
    id: "code",
    mood: "focus",
    bubble: "Ci sono.",
    keywords: ["code", "debug", "bugfix", "react", "astro", "css", "deploy"],
  },
  {
    id: "idea",
    mood: "wow",
    bubble: "Aha.",
    keywords: ["idea", "light", "lamp", "genius", "think"],
  },
  {
    id: "love",
    mood: "felice",
    bubble: "Morbido.",
    keywords: ["love", "heart", "cuore", "amore", "cute"],
  },
  {
    id: "rain",
    mood: "pensieroso",
    bubble: "Plin.",
    keywords: ["rain", "pioggia", "sad", "cry", "triste", "worried"],
  },
  {
    id: "rocket",
    mood: "entusiasta",
    bubble: "Via.",
    keywords: ["rocket", "launch", "lancio", "go", "fly"],
  },
  {
    id: "bug",
    mood: "determinato",
    bubble: "Lo sistemo.",
    keywords: ["bug", "errore", "error", "broken", "fix"],
  },
  {
    id: "paint",
    mood: "curioso",
    bubble: "Colore.",
    keywords: ["paint", "draw", "disegna", "color", "art"],
  },
  {
    id: "book",
    mood: "pensieroso",
    bubble: "Leggo.",
    keywords: ["book", "libro", "read", "study", "learn"],
  },
  {
    id: "money",
    mood: "wow",
    bubble: "Oh.",
    keywords: ["money", "cash", "soldi", "budget", "price"],
  },
  {
    id: "fire",
    mood: "determinato",
    bubble: "Caldo.",
    keywords: ["fire", "fuoco", "hot", "burn", "intense"],
  },
  {
    id: "hug",
    mood: "felice",
    bubble: "Abbraccio.",
    keywords: ["hug", "abbraccio", "comfort", "friend"],
  },
  {
    id: "magic",
    mood: "wow",
    bubble: "Magia.",
    keywords: ["magic", "magia", "spell", "sparkle", "wizard"],
  },
];

export const DEFAULT_REACTION: ThinkyReaction = {
  id: "default",
  mood: "curioso",
  bubble: "Mh.",
  keywords: [],
};

export function detectThinkyReaction(text: string) {
  const normalized = text.trim().toLowerCase();

  if (!normalized) {
    return DEFAULT_REACTION;
  }

  return (
    THINKY_REACTIONS.find((reaction) =>
      reaction.keywords.some((keyword) => normalized.includes(keyword)),
    ) ?? DEFAULT_REACTION
  );
}
