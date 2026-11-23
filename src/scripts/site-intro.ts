/**
 * One 3000ms timeline: code rises and scrambles, "Donato Molino" resolves,
 * compiles into the cloned red header square, two red waves expand, and the
 * logo travels to the header while a circular iris reveals the static page.
 * There are no animation dependencies; every track uses Web Animations API.
 */
type DonatoWindow = Window & {
  donatoIntro?: Promise<void>;
  __donatoIntroEligible?: boolean;
  __donatoIntroBootstrapTime?: number;
  __donatoIntroResolve?: () => void;
};

const win = window as DonatoWindow;
const root = document.documentElement;
const DURATION = 3000;
const glyphs = "{}[]/=*#$01";
const codeTokens = [
  "Astro", "getCollection", "getStaticPaths", "render", "BlogLayout",
  "projects", "latestPosts", "journalDate", "showcase", "pubDate",
  "canonical", "lang", "translationKey", "sitemap", "build",
  "interface", "deploy", "TypeScript", "React", "Cloudflare",
  "AI", "Unity", "Product", "const", "await", "return",
  "filter", "sort", "map", "true", "()", "=>", "{}", "[]",
];

function randomFrom<T>(values: readonly T[]): T {
  return values[Math.floor(Math.random() * values.length)];
}

function makeTokens(): string {
  return Array.from({ length: 5 }, () => randomFrom(codeTokens)).join("  ");
}

function startIntro(): void {
  if (!win.__donatoIntroEligible || !root.classList.contains("intro")) return;

  const overlay = document.querySelector<HTMLElement>("[data-dm-intro-overlay]");
  const code = overlay?.querySelector<HTMLElement>("[data-dm-intro-code]");
  const word = overlay?.querySelector<HTMLElement>("[data-dm-intro-word]");
  const lettersRoot = overlay?.querySelector<HTMLElement>("[data-dm-intro-letters]");
  const cursor = overlay?.querySelector<HTMLElement>("[data-dm-intro-cursor]");
  const logoHolder = overlay?.querySelector<HTMLElement>("[data-dm-intro-logo]");
  const wavesHolder = overlay?.querySelector<HTMLElement>("[data-dm-intro-waves]");
  const logo = document.querySelector<HTMLElement>(".home-header .home-brand .signal-square");
  const main = document.querySelector<HTMLElement>("main");

  const noIntro = () => {
    root.classList.remove("intro", "intro-running");
    win.__donatoIntroResolve?.();
  };

  if (!overlay || !code || !word || !lettersRoot || !cursor ||
      !logoHolder || !wavesHolder || !logo || !main ||
      !Element.prototype.animate ||
      matchMedia("(prefers-reduced-motion: reduce)").matches ||
      // No late intro after the flash guard has nearly completed.
      performance.now() - (win.__donatoIntroBootstrapTime ?? 0) > 1500) {
    noIntro();
    return;
  }

  const vw = document.documentElement.clientWidth;
  const vh = window.innerHeight;
  const centerX = vw / 2;
  const centerY = vh / 2;
  const logoBounds = logo.getBoundingClientRect();
  const logoStyle = getComputedStyle(logo);
  const visible =
    logoStyle.visibility !== "hidden" && logoStyle.display !== "none" &&
    logoBounds.width > 0 && logoBounds.height > 0 &&
    logoBounds.top >= 0 && logoBounds.bottom <= vh &&
    logoBounds.left >= 0 && logoBounds.right <= vw;
  const landingX = visible ? logoBounds.left + logoBounds.width / 2 : centerX;
  const landingY = visible ? logoBounds.top + logoBounds.height / 2 : centerY;
  const sourceVisibility = logo.style.visibility;
  const originalOrigin = main.style.transformOrigin;
  const mainBounds = main.getBoundingClientRect();

  let completed = false;
  let raf = 0;
  let fallbackTimer = 0;
  const animations: Animation[] = [];
  let sourceRevealed = !visible;
  const stopResize = () => finish();
  const stopPagehide = () => finish();

  function animate(el: Element, frames: Keyframe[]) {
    const animation = el.animate(frames, { duration: DURATION, fill: "forwards" });
    animations.push(animation);
    return animation;
  }

  function finish() {
    if (completed) return;
    completed = true;
    window.cancelAnimationFrame(raf);
    window.clearTimeout(fallbackTimer);
    window.removeEventListener("resize", stopResize);
    window.removeEventListener("pagehide", stopPagehide);
    for (const animation of animations) animation.cancel();
    logo.style.visibility = sourceVisibility;
    main.style.transformOrigin = originalOrigin;
    root.classList.remove("intro", "intro-running");
    overlay.remove();
    win.__donatoIntroResolve?.();
  }

  // Copy the real header mark. The size, background and roundness come
  // from the live DOM rather than from a second hand-drawn logo.
  const size = Math.max(1, logoBounds.width);
  const clone = logo.cloneNode(true) as HTMLElement;
  clone.removeAttribute("id");
  clone.style.width = "100%";
  clone.style.height = "100%";
  clone.style.background = logoStyle.backgroundColor;
  clone.style.borderRadius = "inherit";
  logoHolder.style.left = centerX - size / 2 + "px";
  logoHolder.style.top = centerY - size / 2 + "px";
  logoHolder.style.width = size + "px";
  logoHolder.style.height = size + "px";
  logoHolder.style.borderRadius = logoStyle.borderRadius;
  logoHolder.appendChild(clone);

  // The existing Astro content is present from the first HTML render.
  // Only its visual presentation is eased in while the iris opens.
  main.style.transformOrigin =
    landingX - mainBounds.left + "px " + (landingY - mainBounds.top) + "px";

  const lines = Array.from({ length: 15 }, () => {
    const line = document.createElement("div");
    line.className = "dm-intro-code-line";
    line.textContent = makeTokens();
    code.appendChild(line);
    return line;
  });

  // For a 13-character name, compress the per-character interval so all
  // letters settle before the requested 0.37 compile phase.
  const name = Array.from("Donato Molino");
  const glyphElements = name.map((letter) => {
    const span = document.createElement("span");
    span.className = "dm-intro-letter";
    span.textContent = letter === " " ? "\u00a0" : randomFrom(Array.from(glyphs));
    lettersRoot.appendChild(span);
    return span;
  });
  const count = name.filter((letter) => letter !== " ").length;
  const step = Math.min(.055, .265 / Math.max(1, count - 1));
  let characterIndex = 0;
  const revealAt = name.map((letter) =>
    letter === " " ? 0 : .09 + characterIndex++ * step
  );

  function addWave(x: number, y: number, from: number, until: number, diameter: number) {
    const ring = document.createElement("span");
    ring.className = "dm-intro-wave";
    ring.style.left = x + "px";
    ring.style.top = y + "px";
    wavesHolder.appendChild(ring);
    animate(ring, [
      { width: "0px", height: "0px", opacity: 0, offset: 0 },
      { width: "0px", height: "0px", opacity: .76, offset: from },
      { width: diameter + "px", height: diameter + "px", opacity: 0, offset: until },
      { width: diameter + "px", height: diameter + "px", opacity: 0, offset: 1 },
    ]);
  }

  try {
    // Finish building the overlay before replacing the first-paint guard.
    root.classList.add("intro-running");
    if (visible) logo.style.visibility = "hidden";

    const timingLines: Keyframe[] = [
      { transform: "translateY(24px)", opacity: 0, offset: 0 },
      { transform: "translateY(0px)", opacity: 1, offset: .10 },
      { transform: "translateY(0px)", opacity: 1, offset: .36 },
      { transform: "translateY(-30px)", opacity: 0, offset: .48 },
      { transform: "translateY(-30px)", opacity: 0, offset: 1 },
    ];
    lines.forEach((line) => animate(line, timingLines));

    animate(word, [
      { opacity: 0, filter: "blur(15px)", transform: "translateY(24px) scale(.96)", offset: 0 },
      { opacity: 1, filter: "blur(0px)", transform: "translateY(0) scale(1)", offset: .08 },
      { opacity: 1, filter: "blur(0px)", transform: "translateY(0) scale(1)", offset: .37, easing: "ease-out" },
      { opacity: 1, filter: "blur(0px)", transform: "translateY(0) scale(1.05)", color: "var(--dm-intro-accent)", offset: .40 },
      { opacity: 0, filter: "blur(20px)", transform: "translateY(0) scale(.25)", color: "var(--dm-intro-accent)", offset: .47 },
      { opacity: 0, transform: "translateY(0) scale(.25)", offset: 1 },
    ]);

    const large = Math.max(1, .34 * Math.min(vw, vh) / size);
    const dx = landingX - centerX;
    const dy = landingY - centerY;
    if (visible) {
      animate(logoHolder, [
        { transform: "translate(0px, 0px) scale(0)", opacity: 0, offset: 0 },
        { transform: "translate(0px, 0px) scale(0)", opacity: 0, offset: .40, easing: "cubic-bezier(.34,1.56,.64,1)" },
        { transform: "translate(0px, 0px) scale(" + large + ")", opacity: 1, offset: .58 },
        { transform: "translate(0px, 0px) scale(" + large + ")", opacity: 1, offset: .64, easing: "cubic-bezier(.55,0,.3,1)" },
        { transform: "translate(" + dx + "px, " + dy + "px) scale(1.2)", opacity: 1, offset: .90, easing: "cubic-bezier(.2,.8,.2,1)" },
        { transform: "translate(" + dx + "px, " + dy + "px) scale(1)", opacity: 1, offset: 1 },
      ]);
    } else {
      // No jump to an off-screen header on restored, scrolled pages.
      animate(logoHolder, [
        { transform: "scale(0)", opacity: 0, offset: 0 },
        { transform: "scale(0)", opacity: 0, offset: .40, easing: "cubic-bezier(.34,1.56,.64,1)" },
        { transform: "scale(" + large + ")", opacity: 1, offset: .58 },
        { transform: "scale(" + large + ")", opacity: 1, offset: .64 },
        { transform: "scale(" + (large * 1.3) + ")", opacity: 0, offset: .90 },
        { transform: "scale(" + (large * 1.3) + ")", opacity: 0, offset: 1 },
      ]);
    }

    const diagonal = Math.hypot(vw, vh) * 1.1;
    addWave(centerX, centerY, .44, .80, diagonal);
    addWave(centerX, centerY, .50, .86, diagonal);
    if (visible) addWave(landingX, landingY, .89, 1, size * 3);

    const radius = Math.hypot(vw, vh) * 1.1;
    const full = "circle(" + radius + "px at " + landingX + "px " + landingY + "px)";
    const closed = "circle(0px at " + landingX + "px " + landingY + "px)";
    const curtain = animate(overlay, [
      { clipPath: full, offset: 0 },
      { clipPath: full, offset: .64, easing: "cubic-bezier(.55,0,.3,1)" },
      { clipPath: closed, offset: .90 },
      { clipPath: closed, offset: 1 },
    ]);

    animate(main, [
      { transform: "scale(.96)", offset: 0 },
      { transform: "scale(.96)", offset: .64, easing: "cubic-bezier(.55,0,.3,1)" },
      { transform: "scale(1)", offset: .94 },
      { transform: "scale(1)", offset: 1 },
    ]);

    // Use one shared wall-clock progress for scrambles and source hand-off.
    const started = performance.now();
    const tick = (now: number) => {
      if (completed) return;
      const elapsed = Math.max(0, now - started);
      const phase = Math.min(1, elapsed / DURATION);

      if (phase <= .40) {
        name.forEach((letter, i) => {
          glyphElements[i].textContent = letter === " " ? "\u00a0"
            : phase >= revealAt[i] ? letter : randomFrom(Array.from(glyphs));
        });
        cursor.style.visibility = phase < .38 && Math.floor(elapsed / 160) % 2 === 0
          ? "visible" : "hidden";
      } else {
        cursor.style.visibility = "hidden";
        name.forEach((letter, i) => {
          glyphElements[i].textContent = letter === " " ? "\u00a0" : letter;
        });
      }

      if (phase < .48) {
        const first = Math.floor(Math.random() * lines.length);
        const second = (first + 1 + Math.floor(Math.random() * (lines.length - 1))) % lines.length;
        lines[first].textContent = makeTokens();
        lines[second].textContent = makeTokens();
      }

      if (phase >= .90 && !sourceRevealed) {
        sourceRevealed = true;
        logo.style.visibility = sourceVisibility;
      }

      if (phase < 1) raf = requestAnimationFrame(tick);
    };

    window.addEventListener("resize", stopResize, { once: true });
    window.addEventListener("pagehide", stopPagehide, { once: true });
    raf = window.requestAnimationFrame(tick);
    void curtain.finished.then(finish, finish);
    fallbackTimer = window.setTimeout(finish, DURATION + 300);
  } catch (_error) {
    finish();
  }
}

startIntro();
