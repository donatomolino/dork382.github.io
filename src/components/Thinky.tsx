import type { CSSProperties } from "react";
import type { MoodId } from "./thinky-moods";
import type { ThinkyReactionId } from "./thinky-reactions";

type ThinkyProps = {
  mood: MoodId;
  color: string;
  eyeOffset: { x: number; y: number };
  look: { x: number; y: number };
  isNear: boolean;
  isPressed: boolean;
  isBlinking: boolean;
  activity: number;
  burstKey: number;
  reactionId: ThinkyReactionId;
  reactionKey: number;
  thought: string;
  isThinking: boolean;
  onClick: () => void;
};

const expressiveMoods = new Set<MoodId>([
  "curioso",
  "pensieroso",
  "wow",
  "determinato",
  "felice",
  "entusiasta",
  "focus",
]);

export default function Thinky({
  mood,
  color,
  eyeOffset,
  look,
  isNear,
  isPressed,
  isBlinking,
  activity,
  burstKey,
  reactionId,
  reactionKey,
  thought,
  isThinking,
  onClick,
}: ThinkyProps) {
  const eyeX = eyeOffset.x;
  const eyeY = eyeOffset.y;
  const showLeftEye = mood !== "entusiasta";
  const showCodePanel = mood === "focus";
  const showAccentLines = mood === "wow" || mood === "entusiasta";
  const showThinkingPaw = mood === "curioso" || mood === "pensieroso";

  return (
    <button
      className={[
        "thinky-shell",
        `thinky-shell--${mood}`,
        isNear ? "is-near" : "",
        isPressed ? "is-pressed" : "",
        isBlinking ? "is-blinking" : "",
        isThinking ? "is-thinking" : "",
        reactionId !== "default" ? `thinky-shell--react-${reactionId}` : "",
      ]
        .filter(Boolean)
        .join(" ")}
      style={
        {
          "--mood-color": color,
          "--look-x": look.x,
          "--look-y": look.y,
          "--tilt-x": `${look.y * -7}deg`,
          "--tilt-y": `${look.x * 9}deg`,
          "--glow-shift-x": `${look.x * 9}px`,
          "--glow-shift-y": `${look.y * 6}px`,
          "--orbit-activity": activity,
        } as CSSProperties
      }
      type="button"
      aria-label="Change Thinky's mood"
      onClick={onClick}
    >
      {thought && (
        <span className="thinky-thought" key={`thought-${reactionKey}`} aria-live="polite">
          {thought}
        </span>
      )}
      {reactionId !== "default" && (
        <span
          className={`thinky-prop thinky-prop--${reactionId}`}
          key={`prop-${reactionKey}`}
          aria-hidden="true"
        >
          <span />
          <span />
          <span />
        </span>
      )}
      {burstKey > 0 && (
        <span className="thinky-burst" key={burstKey} aria-hidden="true">
          <span />
          <span />
          <span />
          <span />
          <span />
          <span />
        </span>
      )}
      <svg
        className="thinky-svg"
        viewBox="0 0 420 360"
        role="img"
        aria-labelledby="thinky-title thinky-desc"
      >
        <title id="thinky-title">Thinky</title>
        <desc id="thinky-desc">
          A soft white mascot with a mood-colored forehead glow.
        </desc>

        <defs>
          <radialGradient id="body-shine" cx="38%" cy="24%" r="72%">
            <stop offset="0%" stopColor="#ffffff" />
            <stop offset="38%" stopColor="#fbfdff" />
            <stop offset="74%" stopColor="#eef3fb" />
            <stop offset="100%" stopColor="#d6dde9" />
          </radialGradient>
          <radialGradient id="ear-shine" cx="42%" cy="22%" r="80%">
            <stop offset="0%" stopColor="#ffffff" />
            <stop offset="100%" stopColor="#e9edf5" />
          </radialGradient>
          <radialGradient id="soft-paw" cx="40%" cy="28%" r="70%">
            <stop offset="0%" stopColor="#ffffff" />
            <stop offset="100%" stopColor="#e5e9f2" />
          </radialGradient>
          <filter id="body-shadow" x="-25%" y="-25%" width="150%" height="150%">
            <feDropShadow
              dx="0"
              dy="10"
              stdDeviation="12"
              floodColor="#000000"
              floodOpacity="0.32"
            />
          </filter>
          <filter id="inner-glow" x="-40%" y="-60%" width="180%" height="170%">
            <feGaussianBlur stdDeviation="16" />
          </filter>
          <filter id="line-glow" x="-80%" y="-80%" width="260%" height="260%">
            <feGaussianBlur stdDeviation="2.6" result="blur" />
            <feMerge>
              <feMergeNode in="blur" />
              <feMergeNode in="SourceGraphic" />
            </feMerge>
          </filter>
          <linearGradient id="body-rim" x1="76" y1="88" x2="342" y2="320">
            <stop offset="0%" stopColor="#ffffff" stopOpacity="0.88" />
            <stop offset="52%" stopColor="#ffffff" stopOpacity="0.08" />
            <stop offset="100%" stopColor="#8390a6" stopOpacity="0.42" />
          </linearGradient>
        </defs>

        <ellipse className="thinky-ground" cx="210" cy="325" rx="152" ry="20" />

        <g className="thinky-body-group" filter="url(#body-shadow)">
          <path
            className="thinky-ear thinky-ear--left"
            d="M103 104C84 82 88 42 116 34C144 26 164 56 160 92C157 122 126 132 103 104Z"
            fill="url(#ear-shine)"
          />
          <path
            className="thinky-ear thinky-ear--right"
            d="M317 104C336 82 332 42 304 34C276 26 256 56 260 92C263 122 294 132 317 104Z"
            fill="url(#ear-shine)"
          />
          <path
            className="thinky-body"
            d="M72 211C72 125 126 82 210 82C294 82 348 125 348 211C348 292 303 326 210 326C117 326 72 292 72 211Z"
            fill="url(#body-shine)"
          />
          <path
            className="thinky-body-rim"
            d="M72 211C72 125 126 82 210 82C294 82 348 125 348 211C348 292 303 326 210 326C117 326 72 292 72 211Z"
            fill="none"
            stroke="url(#body-rim)"
          />
          <path
            className="thinky-body-soft-shadow"
            d="M90 245C111 298 154 318 211 318C268 318 312 298 330 246C319 301 280 334 210 334C141 334 103 302 90 245Z"
          />
          <ellipse
            className="thinky-ear-occlusion thinky-ear-occlusion--left"
            cx="132"
            cy="109"
            rx="32"
            ry="14"
          />
          <ellipse
            className="thinky-ear-occlusion thinky-ear-occlusion--right"
            cx="288"
            cy="109"
            rx="32"
            ry="14"
          />
          <ellipse
            className="thinky-forehead-glow"
            cx={210 + look.x * 9}
            cy={132 + look.y * 6}
            rx="110"
            ry="58"
            fill="var(--mood-color)"
            filter="url(#inner-glow)"
          />
          <path
            className="thinky-highlight"
            d="M111 139C124 109 156 92 196 89C156 98 131 119 119 151C104 191 111 229 132 256C96 230 91 185 111 139Z"
          />
          <ellipse
            className="thinky-cheek-light"
            cx={276 + look.x * 4}
            cy={229 + look.y * 3}
            rx="44"
            ry="32"
          />

          <path
            className="thinky-foot thinky-foot--left"
            d="M86 303C92 281 116 271 134 286C150 300 142 323 116 326C93 329 81 321 86 303Z"
            fill="url(#soft-paw)"
          />
          <path
            className="thinky-foot thinky-foot--right"
            d="M334 303C328 281 304 271 286 286C270 300 278 323 304 326C327 329 339 321 334 303Z"
            fill="url(#soft-paw)"
          />

          {showThinkingPaw && (
            <g className="thinky-paw thinky-paw--thinking">
              <path
                d="M105 286C88 273 91 244 110 225C125 210 143 214 146 230C151 254 135 291 117 293C113 294 109 291 105 286Z"
                fill="url(#soft-paw)"
              />
              <path
                d="M115 234C104 246 99 261 102 275"
                className="thinky-paw-line"
              />
            </g>
          )}

          <g className="thinky-face">
            {mood === "felice" ? (
              <>
                <path className="thinky-smile-eye" d="M154 205C164 190 178 190 188 205" />
                <path className="thinky-smile-eye" d="M232 205C242 190 256 190 266 205" />
              </>
            ) : (
              <>
                {showLeftEye ? (
                  <ellipse
                    className="thinky-eye thinky-eye--left"
                    cx={165 + eyeX}
                    cy={202 + eyeY}
                    rx={mood === "focus" ? 13 : 16}
                    ry={mood === "focus" ? 19 : 28}
                  />
                ) : (
                  <path className="thinky-wink" d="M149 199C160 191 175 192 186 203" />
                )}
                <ellipse
                  className="thinky-eye thinky-eye--right"
                  cx={255 + eyeX}
                  cy={202 + eyeY}
                  rx={mood === "focus" ? 13 : 16}
                  ry={mood === "focus" ? 19 : 28}
                />
                <ellipse
                  className="thinky-eye-shine"
                  cx={158 + eyeX}
                  cy={190 + eyeY}
                  rx="4"
                  ry="6"
                />
                {showLeftEye && (
                  <ellipse
                    className="thinky-eye-shine"
                    cx={248 + eyeX}
                    cy={190 + eyeY}
                    rx="4"
                    ry="6"
                  />
                )}
              </>
            )}

            {expressiveMoods.has(mood) && (
              <g className="thinky-brows">
                {mood === "determinato" ? (
                  <>
                    <path d="M143 164L184 174" />
                    <path d="M277 164L236 174" />
                  </>
                ) : mood === "pensieroso" ? (
                  <>
                    <path d="M144 170C153 158 166 157 177 166" />
                    <path d="M241 166C252 157 265 158 274 170" />
                  </>
                ) : mood === "focus" ? (
                  <>
                    <path d="M145 170H178" />
                    <path d="M242 170H275" />
                  </>
                ) : mood !== "felice" ? (
                  <>
                    <path d="M148 170C158 160 172 160 182 170" />
                    <path d="M238 170C248 160 262 160 272 170" />
                  </>
                ) : null}
              </g>
            )}

            <g className="thinky-mouth">
              {mood === "wow" ? (
                <ellipse cx="210" cy="249" rx="19" ry="28" />
              ) : mood === "felice" || mood === "entusiasta" ? (
                <path d="M186 244C197 268 224 269 237 244C226 255 200 255 186 244Z" />
              ) : mood === "determinato" ? (
                <path d="M190 251C201 240 219 240 230 251" />
              ) : mood === "pensieroso" || mood === "curioso" ? (
                <path d="M191 247C201 239 215 239 225 247" />
              ) : mood === "focus" ? (
                <path d="M194 246H226" />
              ) : (
                <path d="M194 250H226" />
              )}
            </g>
          </g>

          {showCodePanel && (
            <g className="thinky-code-panel">
              <rect x="122" y="260" width="176" height="72" rx="10" />
              <path d="M178 296L160 281L178 266" />
              <path d="M242 266L260 281L242 296" />
              <path d="M219 263L202 299" />
              <path d="M277 303H287" />
            </g>
          )}
        </g>

        {showAccentLines && (
          <g
            className="thinky-accent-lines"
            stroke="var(--mood-color)"
            filter="url(#line-glow)"
          >
            {mood === "wow" ? (
              <>
                <path d="M267 35L263 67" />
                <path d="M298 48L280 74" />
                <path d="M322 76L292 87" />
              </>
            ) : (
              <>
                <path d="M313 32L303 65" />
                <path d="M345 50L323 74" />
              </>
            )}
          </g>
        )}
      </svg>
    </button>
  );
}
