import type { ComponentProps, CSSProperties, PointerEvent } from "react";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import Thinky from "./Thinky";
import {
  IDLE_MOOD,
  MOODS,
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

type FormSubmitEvent = Parameters<
  NonNullable<ComponentProps<"form">["onSubmit"]>
>[0];

const IDLE_TIMEOUT = 7600;
const MAX_EYE_OFFSET = 8;
const NEAR_DISTANCE = 260;
const INITIAL_MESSAGES: ChatMessage[] = [
  {
    id: 1,
    speaker: "thinky",
    text: "Say something and I will glow with the mood I feel.",
    mood: IDLE_MOOD,
  },
];

function clamp(value: number, min: number, max: number) {
  return Math.min(Math.max(value, min), max);
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
  const [chatValue, setChatValue] = useState("");
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

  const handleChatSubmit = useCallback(
    (event: FormSubmitEvent) => {
      event.preventDefault();

      const text = chatValue.trim();

      if (!text) {
        return;
      }

      const nextMood = detectMoodFromText(text);
      const nextMoodConfig = getMoodById(nextMood);
      const userMessage: ChatMessage = {
        id: messageIdRef.current++,
        speaker: "you",
        text,
      };
      const thinkyMessage: ChatMessage = {
        id: messageIdRef.current++,
        speaker: "thinky",
        text: nextMoodConfig.reply,
        mood: nextMood,
      };

      setMessages((current) => [...current.slice(-4), userMessage, thinkyMessage]);
      setChatValue("");
      changeMood(nextMood);
      setIsPressed(true);
      wakeThinky();
      window.clearTimeout(pressTimerRef.current);
      pressTimerRef.current = window.setTimeout(() => setIsPressed(false), 420);
    },
    [changeMood, chatValue, wakeThinky],
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
    };
  }, [resetIdleTimer]);

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
        <p className="thinky-kicker">Mood-aware companion</p>
        <h1>Thinky</h1>
        <p>Talk to him. Move your mouse. Watch the glow react.</p>
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
          onClick={handleRandomMood}
        />
      </section>

      <section className="thinky-console" aria-label="Talk to Thinky">
        <div className="thinky-chat-log" aria-live="polite">
          {messages.map((message) => (
            <p
              className={`thinky-message thinky-message--${message.speaker}`}
              key={message.id}
            >
              <span>{message.speaker === "you" ? "You" : "Thinky"}</span>
              {message.text}
            </p>
          ))}
        </div>

        <form className="thinky-chat-form" onSubmit={handleChatSubmit}>
          <input
            type="text"
            value={chatValue}
            maxLength={180}
            placeholder="Tell Thinky something..."
            aria-label="Message Thinky"
            onChange={(event) => setChatValue(event.target.value)}
          />
          <button type="submit">Send</button>
        </form>

        <nav className="thinky-mood-bar" aria-label="Thinky moods">
          {MOODS.map((item) => (
            <button
              className={item.id === mood.id ? "is-active" : ""}
              key={item.id}
              type="button"
              aria-label={`Set Thinky mood to ${item.label}`}
              aria-pressed={item.id === mood.id}
              style={{ "--mood-button-color": item.color } as CSSProperties}
              onClick={() => changeMood(item.id)}
            >
              <span className="thinky-mood-dot" aria-hidden="true" />
              <span className="thinky-mood-copy">
                <strong>{item.label}</strong>
                <small>{item.colorName}</small>
              </span>
            </button>
          ))}
        </nav>
      </section>
    </main>
  );
}
