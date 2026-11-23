/*
 * Progressive enhancements for the portfolio home page.
 * The content renders normally with JavaScript disabled; WebGL is optional.
 */
const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
const finePointer = window.matchMedia("(hover: hover) and (pointer: fine)");
const intro = document.querySelector<HTMLElement>("[data-code-intro]");

function initIntro() {
  if (!intro || !document.documentElement.classList.contains("intro-playing")) return;

  intro.setAttribute("aria-hidden", "false");
  const button = intro.querySelector<HTMLButtonElement>("[data-intro-skip]");
  const background = Array.from(document.querySelectorAll<HTMLElement>(
    "body > .home-header, body > main, body > .home-footer, body > .home-skip"
  ));
  background.forEach((element) => { element.inert = true; });
  button?.focus({ preventScroll: true });
  let finished = false;
  let startExit: number;
  let endExit: number;

  const complete = (skip = false) => {
    if (finished) return;
    finished = true;
    window.clearTimeout(startExit);
    window.clearTimeout(endExit);
    if (skip) {
      document.documentElement.classList.remove("intro-playing", "intro-leaving");
      intro.setAttribute("aria-hidden", "true");
    } else {
      document.documentElement.classList.add("intro-leaving");
      endExit = window.setTimeout(() => {
        document.documentElement.classList.remove("intro-playing", "intro-leaving");
        intro.setAttribute("aria-hidden", "true");
      }, 260);
    }
    background.forEach((element) => { element.inert = false; });
    if (document.activeElement === button) {
      document.querySelector<HTMLElement>("main#top")?.focus({ preventScroll: true });
    }
    document.removeEventListener("keydown", onKeyDown);
  };

  function onKeyDown(event: KeyboardEvent) {
    if (event.key === "Escape") complete(true);
  }

  document.addEventListener("keydown", onKeyDown);
  button?.addEventListener("click", () => complete(true), { once: true });
  startExit = window.setTimeout(() => complete(false), 1750);
  if (reduceMotion.matches) complete(true);
}

function initCardTilt() {
  if (reduceMotion.matches || !finePointer.matches) return;
  const cards = document.querySelectorAll<HTMLElement>(".project-card");

  cards.forEach((card) => {
    let frame = 0;
    let eventX = 0;
    let eventY = 0;

    const update = () => {
      frame = 0;
      const bounds = card.getBoundingClientRect();
      const nx = (eventX - bounds.left) / bounds.width;
      const ny = (eventY - bounds.top) / bounds.height;
      const x = Math.max(-1, Math.min(1, (ny - 0.5) * -2));
      const y = Math.max(-1, Math.min(1, (nx - 0.5) * 2));

      card.style.setProperty("--card-rx", (x * 4.2).toFixed(2) + "deg");
      card.style.setProperty("--card-ry", (y * 5).toFixed(2) + "deg");
      card.style.setProperty("--pointer-x", (nx * 100).toFixed(1) + "%");
      card.style.setProperty("--pointer-y", (ny * 100).toFixed(1) + "%");
    };

    card.addEventListener("pointermove", (event) => {
      eventX = event.clientX;
      eventY = event.clientY;
      if (!frame) frame = window.requestAnimationFrame(update);
    }, { passive: true });

    card.addEventListener("pointerleave", () => {
      window.cancelAnimationFrame(frame);
      frame = 0;
      card.style.removeProperty("--card-rx");
      card.style.removeProperty("--card-ry");
      card.style.removeProperty("--pointer-x");
      card.style.removeProperty("--pointer-y");
    });
  });
}

function initReveals() {
  if (reduceMotion.matches || !("IntersectionObserver" in window)) return;
  const sections = document.querySelectorAll<HTMLElement>(
    ".section-header, .project-card, .writing-card, .manifesto-callout, .capabilities-grid article, .about-section, .contact-section"
  );
  const observer = new IntersectionObserver((entries, current) => {
    entries.forEach((entry) => {
      if (entry.isIntersecting) {
        entry.target.classList.add("motion-in-view");
        current.unobserve(entry.target);
      }
    });
  }, { threshold: 0.08, rootMargin: "0px 0px -26px 0px" });

  sections.forEach((section, i) => {
    section.classList.add("motion-reveal");
    section.style.setProperty("--reveal-delay", (i % 4) * 35 + "ms");
    observer.observe(section);
  });
}

async function initWebGLArtifact() {
  const artifact = document.querySelector<HTMLElement>("[data-hero-artifact]");
  const hero = document.querySelector<HTMLElement>(".home-hero");
  const canvas = artifact?.querySelector<HTMLCanvasElement>("[data-hero-canvas]");
  const saveData = (navigator as Navigator & { connection?: { saveData?: boolean } }).connection?.saveData;

  if (!artifact || !hero || !canvas ||
    reduceMotion.matches || !finePointer.matches ||
    window.innerWidth < 960 || saveData) return;

  // Defer WebGL and the Three.js bundle until the art is actually visible.
  if (!("IntersectionObserver" in window)) return;

  let started = false;
  let inView = false;
  let running = false;
  let disposed = false;
  let frame = 0;
  let renderer: import("three").WebGLRenderer | undefined;
  let stageObserver: ResizeObserver | undefined;
  let disposeResources: (() => void) | undefined;
  let play: (() => void) | undefined;
  let pause: (() => void) | undefined;

  const onVisibility = () => {
    if (inView && !document.hidden) play?.();
    else pause?.();
  };

  const onEntry = new IntersectionObserver((entries) => {
    inView = entries[0]?.isIntersecting ?? false;
    if (inView && !started) {
      started = true;
      void boot();
    } else {
      onVisibility();
    }
  }, { threshold: 0.02 });
  onEntry.observe(artifact);

  async function boot() {
    if (!canvas || !artifact || !hero) return;
    try {
      const THREE = await import("three");
      if (disposed || reduceMotion.matches) return;

      renderer = new THREE.WebGLRenderer({
        canvas,
        antialias: true,
        alpha: true,
        powerPreference: "low-power",
      });
      renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.5));
      renderer.outputColorSpace = THREE.SRGBColorSpace;

      const scene = new THREE.Scene();
      const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 50);
      camera.position.set(0, 0, 7.5);

      const root = new THREE.Group();
      scene.add(root);

      const crystalGeometry = new THREE.IcosahedronGeometry(1.13, 1);
      const crystalMaterial = new THREE.MeshPhysicalMaterial({
        color: 0x3157dc,
        metalness: 0.32,
        roughness: 0.24,
        flatShading: true,
        clearcoat: 0.9,
        clearcoatRoughness: 0.12,
      });
      const crystal = new THREE.Mesh(crystalGeometry, crystalMaterial);
      root.add(crystal);

      const edgesGeometry = new THREE.EdgesGeometry(crystalGeometry);
      const edgesMaterial = new THREE.LineBasicMaterial({
        color: 0xe7edff,
        transparent: true,
        opacity: 0.68,
      });
      const edges = new THREE.LineSegments(edgesGeometry, edgesMaterial);
      crystal.add(edges);

      const rings = [
        { radius: 1.61, x: 1.1, y: 0.12, color: 0x3157dc },
        { radius: 1.78, x: 0.48, y: 1.2, color: 0xcc3045 },
        { radius: 1.9, x: -0.36, y: 1.05, color: 0x91aae7 },
      ].map((def) => {
        const geometry = new THREE.TorusGeometry(def.radius, 0.014, 7, 128);
        const material = new THREE.MeshBasicMaterial({ color: def.color, transparent: true, opacity: 0.86 });
        const mesh = new THREE.Mesh(geometry, material);
        mesh.rotation.set(def.x, def.y, 0);
        root.add(mesh);
        return { geometry, material, mesh };
      });

      const satelliteGeometry = new THREE.SphereGeometry(0.095, 14, 12);
      const satelliteMaterial = new THREE.MeshBasicMaterial({ color: 0xcc3045 });
      const satellite = new THREE.Mesh(satelliteGeometry, satelliteMaterial);
      satellite.position.set(1.6, 0.8, 0.18);
      root.add(satellite);

      scene.add(new THREE.AmbientLight(0xffffff, 1.9));
      const keyLight = new THREE.DirectionalLight(0xeaf1ff, 3.3);
      keyLight.position.set(-3, 4, 5);
      scene.add(keyLight);
      const rimLight = new THREE.PointLight(0xd4ef8d, 40);
      rimLight.position.set(3, -3, 4);
      scene.add(rimLight);

      const aim = { x: 0, y: 0 };
      let last = 0;

      const onPointer = (event: PointerEvent) => {
        const b = hero.getBoundingClientRect();
        aim.x = Math.max(-1, Math.min(1, ((event.clientX - b.left) / b.width - 0.5) * 2));
        aim.y = Math.max(-1, Math.min(1, ((event.clientY - b.top) / b.height - 0.5) * 2));
      };
      const onPointerLeave = () => { aim.x = 0; aim.y = 0; };
      hero.addEventListener("pointermove", onPointer, { passive: true });
      hero.addEventListener("pointerleave", onPointerLeave);

      const resize = () => {
        if (!renderer || !artifact) return;
        const width = Math.max(1, artifact.clientWidth);
        const height = Math.max(1, artifact.clientHeight);
        camera.aspect = width / height;
        camera.updateProjectionMatrix();
        renderer.setSize(width, height, false);
      };
      resize();
      stageObserver = new ResizeObserver(resize);
      stageObserver.observe(artifact);

      const tick = (timestamp: number) => {
        if (!running || !renderer) return;
        const dt = Math.min(40, timestamp - (last || timestamp));
        last = timestamp;

        root.rotation.x += ((-aim.y * 0.38) - root.rotation.x) * 0.065;
        root.rotation.y += ((aim.x * 0.52) - root.rotation.y) * 0.065;
        crystal.rotation.y += dt * 0.00022;
        crystal.rotation.x += dt * 0.0001;
        rings[0].mesh.rotation.z += dt * 0.00011;
        rings[1].mesh.rotation.y += dt * 0.00008;
        rings[2].mesh.rotation.x += dt * 0.0001;

        renderer.render(scene, camera);
        frame = window.requestAnimationFrame(tick);
      };

      play = () => {
        if (running || disposed) return;
        running = true;
        last = 0;
        frame = window.requestAnimationFrame(tick);
      };
      pause = () => {
        running = false;
        window.cancelAnimationFrame(frame);
      };

      disposeResources = () => {
        hero.removeEventListener("pointermove", onPointer);
        hero.removeEventListener("pointerleave", onPointerLeave);
        stageObserver?.disconnect();
        crystalGeometry.dispose();
        crystalMaterial.dispose();
        edgesGeometry.dispose();
        edgesMaterial.dispose();
        satelliteGeometry.dispose();
        satelliteMaterial.dispose();
        for (const ring of rings) { ring.geometry.dispose(); ring.material.dispose(); }
        renderer?.dispose();
        artifact.classList.remove("hero-artifact--webgl");
      };

      // The CSS model stays visible until the first successful WebGL frame.
      renderer.render(scene, camera);
      artifact.classList.add("hero-artifact--webgl");
      onVisibility();
    } catch (error) {
      // WebGL is an enhancement. An available CSS model remains visible on failure.
      console.warn("3D scene unavailable; using CSS fallback.", error);
      renderer?.dispose();
      artifact?.classList.remove("hero-artifact--webgl");
    }
  }

  window.addEventListener("pagehide", () => {
    disposed = true;
    onEntry.disconnect();
    pause?.();
    disposeResources?.();
    document.removeEventListener("visibilitychange", onVisibility);
  }, { once: true });
  document.addEventListener("visibilitychange", onVisibility);
}

initIntro();
initCardTilt();
initReveals();
void initWebGLArtifact();
