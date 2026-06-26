import type { CSSProperties } from "react";
import type { MoodId } from "./thinky-moods";

type OrbitProps = {
  mood: MoodId;
  color: string;
  isNear: boolean;
  activity: number;
};

export default function Orbit({ mood, color, isNear, activity }: OrbitProps) {
  return (
    <span
      className={[
        "thinky-orbit",
        `thinky-orbit--${mood}`,
        isNear ? "is-near" : "",
      ]
        .filter(Boolean)
        .join(" ")}
      style={
        {
          "--mood-color": color,
          "--orbit-activity": activity,
        } as CSSProperties
      }
      aria-hidden="true"
    >
      <span className="thinky-orbit-core" />
      <span className="thinky-orbit-particle thinky-orbit-particle--one" />
      <span className="thinky-orbit-particle thinky-orbit-particle--two" />
    </span>
  );
}
