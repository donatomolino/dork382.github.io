export type MoodId =
  | "curioso"
  | "pensieroso"
  | "wow"
  | "determinato"
  | "felice"
  | "entusiasta"
  | "focus"
  | "neutrale";

export type Mood = {
  id: MoodId;
  label: string;
  colorName: string;
  color: string;
  personality: string;
  reply: string;
  keywords: string[];
};

export const MOODS: Mood[] = [
  {
    id: "curioso",
    label: "Curious",
    colorName: "Blue",
    color: "#4da3ff",
    personality: "tilted, observant, playful",
    reply: "That makes me curious. Tell me one more detail.",
    keywords: ["why", "how", "what", "maybe", "idea", "explore", "wonder", "?"],
  },
  {
    id: "pensieroso",
    label: "Thoughtful",
    colorName: "Violet",
    color: "#b56cff",
    personality: "slow, reflective, soft",
    reply: "I am thinking about it. There is something interesting in that.",
    keywords: ["think", "consider", "reflect", "unsure", "perhaps", "complex", "question"],
  },
  {
    id: "wow",
    label: "Wow",
    colorName: "Sky",
    color: "#65d9ff",
    personality: "surprised, bright, lifted",
    reply: "Whoa. That got my glow going.",
    keywords: ["wow", "amazing", "surprise", "incredible", "huge", "wild", "unbelievable"],
  },
  {
    id: "determinato",
    label: "Determined",
    colorName: "Red",
    color: "#ff4d5e",
    personality: "stable, sharp, focused",
    reply: "Locked in. We can handle this.",
    keywords: ["must", "ship", "fix", "urgent", "goal", "win", "hard", "challenge", "fight"],
  },
  {
    id: "felice",
    label: "Happy",
    colorName: "Green",
    color: "#6dff79",
    personality: "warm, buoyant, smiling",
    reply: "That feels good. My forehead is basically smiling too.",
    keywords: ["happy", "great", "nice", "love", "thanks", "good", "yay", "beautiful", "success"],
  },
  {
    id: "entusiasta",
    label: "Excited",
    colorName: "Orange",
    color: "#ffad42",
    personality: "bouncy, sparkly, expressive",
    reply: "Yes. I am absolutely here for that.",
    keywords: ["excited", "launch", "party", "awesome", "let's go", "cool", "energy", "!"],
  },
  {
    id: "focus",
    label: "Focus",
    colorName: "Deep blue",
    color: "#336dff",
    personality: "precise, quiet, technical",
    reply: "Focus mode. Give me the task and I will stay sharp.",
    keywords: ["code", "debug", "build", "task", "work", "focus", "terminal", "logic", "system"],
  },
  {
    id: "neutrale",
    label: "Neutral",
    colorName: "White",
    color: "#ffffff",
    personality: "calm, clean, balanced",
    reply: "I am listening.",
    keywords: ["neutral", "calm", "ok", "hello", "hi", "hey"],
  },
];

export const IDLE_MOOD: MoodId = "curioso";

const moodPriority: MoodId[] = [
  "wow",
  "entusiasta",
  "determinato",
  "focus",
  "felice",
  "pensieroso",
  "curioso",
  "neutrale",
];

export function getMoodById(id: MoodId) {
  return MOODS.find((mood) => mood.id === id) ?? MOODS[0];
}

export function getRandomMood(current: MoodId) {
  const nextMoods = MOODS.filter((mood) => mood.id !== current);
  return nextMoods[Math.floor(Math.random() * nextMoods.length)].id;
}

export function detectMoodFromText(text: string): MoodId {
  const normalized = text.trim().toLowerCase();

  if (!normalized) {
    return "neutrale";
  }

  const scores = new Map<MoodId, number>();

  for (const mood of MOODS) {
    let score = 0;

    for (const keyword of mood.keywords) {
      if (normalized.includes(keyword)) {
        score += keyword.length > 2 ? 2 : 1;
      }
    }

    scores.set(mood.id, score);
  }

  if (normalized.includes("sad") || normalized.includes("worried") || normalized.includes("scared")) {
    scores.set("pensieroso", (scores.get("pensieroso") ?? 0) + 3);
  }

  if (normalized.includes("angry") || normalized.includes("serious")) {
    scores.set("determinato", (scores.get("determinato") ?? 0) + 3);
  }

  if (normalized.length > 120) {
    scores.set("pensieroso", (scores.get("pensieroso") ?? 0) + 1);
  }

  const best = moodPriority.reduce(
    (winner, moodId) =>
      (scores.get(moodId) ?? 0) > (scores.get(winner) ?? 0) ? moodId : winner,
    "neutrale" as MoodId,
  );

  return (scores.get(best) ?? 0) > 0 ? best : "curioso";
}
