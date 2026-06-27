import type { ComponentProps, CSSProperties, PointerEvent } from "react";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import Thinky from "./Thinky";
import {
  DEFAULT_REACTION,
  detectThinkyReaction,
  type ThinkyReactionId,
} from "./thinky-reactions";
import {
  IDLE_MOOD,
  detectMoodFromText,
  getMoodById,
  getRandomMood,
  type MoodId,
} from "./thinky-moods";

type ChatMessage = {
  id: number;
  speaker: "you" | "thinky";
  text: string;
  mood?: MoodId;
};

type ThinkyApiLine = {
  text: string;
  mood: MoodId;
  action: ThinkyReactionId | "none";
};

type FormSubmitEvent = Parameters<
  NonNullable<ComponentProps<"form">["onSubmit"]>
>[0];

const THINKY_API_URL = "https://ai.dork3802.workers.dev/";
const MAX_MESSAGE_LENGTH = 255;
const IDLE_TIMEOUT = 7600;
const MAX_EYE_OFFSET = 8;
const NEAR_DISTANCE = 260;
const SUGGESTED_PROMPTS = [
  "Have a cookie.",
  "Eat the pizza!",
  "Take a tiny coffee break.",
  "Show me a bright idea.",
  "Debug this little bug.",
  "Launch the rocket.",
  "Do a tiny dance.",
  "Read this book.",
  "Make some magic.",
  "Paint something bright.",
];
const INITIAL_MESSAGES: ChatMessage[] = [
  {
    id: 1,
    speaker: "you",
    text: "Say something to Thinky.",
    mood: IDLE_MOOD,
  },
];

function clamp(value: number, min: number, max: number) {
  return Math.min(Math.max(value, min), max);
}

function getRandomPrompts(count: number) {
  return [...SUGGESTED_PROMPTS]
    .sort(() => Math.random() - 0.5)
    .slice(0, count);
}

function isMoodId(value: string): value is MoodId {
  return [
    "curioso",
    "pensieroso",
    "wow",
    "determinato",
    "felice",
    "entusiasta",
    "focus",
    "neutrale",
  ].includes(value);
}

function isReactionId(value: string): value is ThinkyReactionId {
  return [
    "cookie",
    "pizza",
    "coffee",
    "sleep",
    "music",
    "code",
    "idea",
    "love",
    "rain",
    "rocket",
    "bug",
    "paint",
    "book",
    "money",
    "fire",
    "hug",
    "magic",
    "default",
  ].includes(value);
}

function normalizeApiLines(value: unknown): ThinkyApiLine[] {
  const candidate = value as { lines?: unknown };
  const lines = Array.isArray(candidate.lines) ? candidate.lines : [];

  return lines
    .slice(0, 3)
    .map((line) => {
      const item = line as { text?: unknown; mood?: unknown; action?: unknown };
      const rawMood = String(item.mood ?? "");
      const rawAction = String(item.action ?? "none");
      const text = String(item.text ?? "")
        .replace(/\s+/g, " ")
        .trim()
        .slice(0, 64);

      return {
        text,
        mood: isMoodId(rawMood) ? rawMood : "curioso",
        action:
          rawAction === "none"
            ? "none"
            : isReactionId(rawAction)
              ? rawAction
              : "none",
      };
    })
    .filter((line) => line.text);
}

function combineThinkyLines(lines: ThinkyApiLine[]): ThinkyApiLine {
  const firstActionLine = lines.find((line) => line.action !== "none");
  const firstMoodLine = lines.find((line) => line.mood !== "neutrale");

  return {
    text: lines.map((line) => line.text).join(" ").trim().slice(0, 190),
    mood: firstMoodLine?.mood ?? lines[0]?.mood ?? "curioso",
    action: firstActionLine?.action ?? "none",
  };
}

function getLocalThinkyLine(text: string): ThinkyApiLine {
  const mood = detectMoodFromText(text);
  const reaction = detectThinkyReaction(text);

  return {
    text: reaction.bubble,
    mood: reaction.id === "default" ? mood : reaction.mood,
    action: reaction.id === "default" ? "none" : reaction.id,
  };
}

export default function ThinkyPlayground() {
  const [activeMood, setActiveMood] = useState<MoodId>(IDLE_MOOD);
  const [eyeOffset, setEyeOffset] = useState({ x: 0, y: 0 });
  const [look, setLook] = useState({ x: 0, y: 0 });
  const [isNear, setIsNear] = useState(false);
  const [isPressed, setIsPressed] = useState(false);
  const [isBlinking, setIsBlinking] = useState(false);
  const [activity, setActivity] = useState(0);
  const [burstKey, setBurstKey] = useState(0);
  const [reactionKey, setReactionKey] = useState(0);
  const [reactionId, setReactionId] = useState<ThinkyReactionId>(
    DEFAULT_REACTION.id,
  );
  const [thought, setThought] = useState("");
  const [isThinking, setIsThinking] = useState(false);
  const [chatValue, setChatValue] = useState("");
  const [suggestedPrompts, setSuggestedPrompts] = useState<string[]>([]);
  const [messages, setMessages] = useState<ChatMessage[]>(INITIAL_MESSAGES);
  const stageRef = useRef<HTMLDivElement | null>(null);
  const targetEyeRef = useRef({ x: 0, y: 0 });
  const eyeRef = useRef({ x: 0, y: 0 });
  const targetLookRef = useRef({ x: 0, y: 0 });
  const lookRef = useRef({ x: 0, y: 0 });
  const idleTimerRef = useRef<number | undefined>(undefined);
  const pressTimerRef = useRef<number | undefined>(undefined);
  const blinkTimerRef = useRef<number | undefined>(undefined);
  const activityTimerRef = useRef<number | undefined>(undefined);
  const thoughtTimerRef = useRef<number | undefined>(undefined);
  const responseTimersRef = useRef<number[]>([]);
  const messageIdRef = useRef(2);

  const mood = useMemo(
    () => getMoodById(activeMood),
    [activeMood],
  );

  const resetIdleTimer = useCallback(() => {
    window.clearTimeout(idleTimerRef.current);
    idleTimerRef.current = window.setTimeout(() => {
      setActiveMood(IDLE_MOOD);
      setIsNear(false);
    }, IDLE_TIMEOUT);
  }, []);

  const changeMood = useCallback(
    (nextMood: MoodId) => {
      setActiveMood(nextMood);
      setBurstKey((key) => key + 1);
      resetIdleTimer();
    },
    [resetIdleTimer],
  );

  const handleRandomMood = useCallback(() => {
    changeMood(getRandomMood(activeMood));
    setIsPressed(true);
    window.clearTimeout(pressTimerRef.current);
    pressTimerRef.current = window.setTimeout(() => setIsPressed(false), 420);
  }, [activeMood, changeMood]);

  const wakeThinky = useCallback(() => {
    setActivity(1);
    window.clearTimeout(activityTimerRef.current);
    activityTimerRef.current = window.setTimeout(() => setActivity(0), 950);
  }, []);

  const playThinkyLine = useCallback(
    (line: ThinkyApiLine) => {
      setReactionId(line.action === "none" ? DEFAULT_REACTION.id : line.action);
      setThought(line.text);
      setReactionKey((key) => key + 1);
      changeMood(line.mood);
      setIsPressed(true);
      wakeThinky();
      window.clearTimeout(pressTimerRef.current);
      window.clearTimeout(thoughtTimerRef.current);
      pressTimerRef.current = window.setTimeout(() => setIsPressed(false), 420);
      thoughtTimerRef.current = window.setTimeout(() => {
        setReactionId(DEFAULT_REACTION.id);
      }, 2600);
    },
    [changeMood, wakeThinky],
  );

  const submitThinkyMessage = useCallback(
    async (rawText: string) => {
      const text = rawText.trim().slice(0, MAX_MESSAGE_LENGTH);

      if (!text || isThinking) {
        return;
      }

      const userMessage: ChatMessage = {
        id: messageIdRef.current++,
        speaker: "you",
        text,
      };

      responseTimersRef.current.forEach((timer) => window.clearTimeout(timer));
      responseTimersRef.current = [];
      setMessages((current) => [...current.slice(-5), userMessage]);
      setChatValue("");
      setSuggestedPrompts([]);
      setIsThinking(true);
      setReactionId(DEFAULT_REACTION.id);
      setThought("...");
      setReactionKey((key) => key + 1);
      changeMood("pensieroso");

      try {
        const response = await fetch(THINKY_API_URL, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            message: text,
            language: "it",
            locale: "it-IT",
          }),
        });

        if (!response.ok) {
          throw new Error(`Thinky API returned ${response.status}`);
        }

        const lines = normalizeApiLines(await response.json());

        if (!lines.length) {
          throw new Error("Thinky API returned no lines");
        }

        playThinkyLine(combineThinkyLines(lines));
      } catch {
        const fallbackLine = getLocalThinkyLine(text);
        playThinkyLine(fallbackLine);
      } finally {
        setIsThinking(false);
      }
    },
    [changeMood, isThinking, playThinkyLine],
  );

  const handleChatSubmit = useCallback(
    (event: FormSubmitEvent) => {
      event.preventDefault();
      void submitThinkyMessage(chatValue);
    },
    [chatValue, submitThinkyMessage],
  );

  useEffect(() => {
    let frame = 0;

    const tick = () => {
      const current = eyeRef.current;
      const target = targetEyeRef.current;
      const lookCurrent = lookRef.current;
      const lookTarget = targetLookRef.current;
      const next = {
        x: current.x + (target.x - current.x) * 0.12,
        y: current.y + (target.y - current.y) * 0.12,
      };
      const nextLook = {
        x: lookCurrent.x + (lookTarget.x - lookCurrent.x) * 0.1,
        y: lookCurrent.y + (lookTarget.y - lookCurrent.y) * 0.1,
      };

      eyeRef.current = next;
      lookRef.current = nextLook;
      setEyeOffset(next);
      setLook(nextLook);
      frame = window.requestAnimationFrame(tick);
    };

    frame = window.requestAnimationFrame(tick);
    resetIdleTimer();

    return () => {
      window.cancelAnimationFrame(frame);
      window.clearTimeout(idleTimerRef.current);
      window.clearTimeout(pressTimerRef.current);
      window.clearTimeout(blinkTimerRef.current);
      window.clearTimeout(activityTimerRef.current);
      window.clearTimeout(thoughtTimerRef.current);
      responseTimersRef.current.forEach((timer) => window.clearTimeout(timer));
    };
  }, [resetIdleTimer]);

  useEffect(() => {
    setSuggestedPrompts(getRandomPrompts(2));
  }, []);

  useEffect(() => {
    const scheduleBlink = () => {
      const delay = activeMood === "focus" ? 5200 : activeMood === "wow" ? 3200 : 4100;

      blinkTimerRef.current = window.setTimeout(() => {
        setIsBlinking(true);
        window.setTimeout(() => {
          setIsBlinking(false);
          scheduleBlink();
        }, activeMood === "pensieroso" ? 190 : 130);
      }, delay + Math.random() * 1400);
    };

    scheduleBlink();

    return () => window.clearTimeout(blinkTimerRef.current);
  }, [activeMood]);

  const handlePointerMove = useCallback(
    (event: PointerEvent<HTMLDivElement>) => {
      const stage = stageRef.current;

      if (!stage) {
        return;
      }

      const rect = stage.getBoundingClientRect();
      const centerX = rect.left + rect.width / 2;
      const centerY = rect.top + rect.height / 2;
      const deltaX = event.clientX - centerX;
      const deltaY = event.clientY - centerY;
      const distance = Math.hypot(deltaX, deltaY);
      const pull = clamp(distance / NEAR_DISTANCE, 0, 1);
      const strength = distance < NEAR_DISTANCE ? 1 : 0.62;

      targetEyeRef.current = {
        x: clamp(
          (deltaX / Math.max(distance, 1)) * MAX_EYE_OFFSET * pull * strength,
          -MAX_EYE_OFFSET,
          MAX_EYE_OFFSET,
        ),
        y: clamp(
          (deltaY / Math.max(distance, 1)) * MAX_EYE_OFFSET * pull * strength,
          -MAX_EYE_OFFSET,
          MAX_EYE_OFFSET,
        ),
      };
      targetLookRef.current = {
        x: clamp(deltaX / NEAR_DISTANCE, -1, 1),
        y: clamp(deltaY / NEAR_DISTANCE, -1, 1),
      };

      setIsNear(distance < NEAR_DISTANCE);
      setActivity(distance < NEAR_DISTANCE ? 1 : 0.45);
      resetIdleTimer();
    },
    [resetIdleTimer],
  );

  const handlePointerLeave = useCallback(() => {
    targetEyeRef.current = { x: 0, y: 0 };
    targetLookRef.current = { x: 0, y: 0 };
    setIsNear(false);
    setActivity(0);
  }, []);

  return (
    <main
      className="thinky-page"
      style={{ "--mood-color": mood.color } as CSSProperties}
      onPointerMove={handlePointerMove}
      onPointerLeave={handlePointerLeave}
    >
      <div className="thinky-ambient" aria-hidden="true" />

      <header className="thinky-header">
        <p className="thinky-kicker">Interactive mascot / 01</p>
        <h1>Thinky</h1>
      </header>

      <section className="thinky-stage" ref={stageRef} aria-live="polite">
        <Thinky
          mood={mood.id}
          color={mood.color}
          eyeOffset={eyeOffset}
          look={look}
          isNear={isNear}
          isPressed={isPressed}
          isBlinking={isBlinking}
          activity={activity}
          burstKey={burstKey}
          reactionId={reactionId}
          reactionKey={reactionKey}
          thought={thought}
          isThinking={isThinking}
          onClick={handleRandomMood}
        />
      </section>

      <section className="thinky-console" aria-label="Talk to Thinky">
        <div className="thinky-chat-log" aria-live="polite">
          {messages
            .filter((message) => message.speaker === "you")
            .slice(-1)
            .map((message) => (
            <p
              className={`thinky-message thinky-message--${message.speaker}`}
              key={message.id}
            >
              <span>{message.speaker === "you" ? "You" : "Thinky"}</span>
              {message.text}
            </p>
          ))}
        </div>

        {suggestedPrompts.length > 0 && (
          <div className="thinky-suggestions" aria-label="Suggested prompts">
            {suggestedPrompts.map((prompt) => (
              <button
                key={prompt}
                type="button"
                disabled={isThinking}
                onClick={() => void submitThinkyMessage(prompt)}
              >
                {prompt}
              </button>
            ))}
          </div>
        )}

        <form className="thinky-chat-form" onSubmit={handleChatSubmit}>
          <input
            type="text"
            value={chatValue}
            maxLength={MAX_MESSAGE_LENGTH}
            placeholder="Tell Thinky something..."
            aria-label="Message Thinky"
            disabled={isThinking}
            onChange={(event) => {
              setSuggestedPrompts([]);
              setChatValue(event.target.value.slice(0, MAX_MESSAGE_LENGTH));
            }}
          />
          <span className="thinky-char-counter" aria-live="polite">
            {chatValue.length}/{MAX_MESSAGE_LENGTH}
          </span>
          <button type="submit" disabled={isThinking}>
            {isThinking ? "Thinking" : "Send"}
          </button>
        </form>
      </section>
    </main>
  );
}
