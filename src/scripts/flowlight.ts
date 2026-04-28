export {};

type FluidBufferKey =
  | "velocityX"
  | "velocityY"
  | "nextVelocityX"
  | "nextVelocityY"
  | "dyeR"
  | "dyeG"
  | "dyeB"
  | "nextDyeR"
  | "nextDyeG"
  | "nextDyeB"
  | "pressure"
  | "nextPressure"
  | "divergence"
  | "curl";

interface Viewport {
  width: number;
  height: number;
  dpr: number;
  cellSize: number;
}

interface FluidField {
  width: number;
  height: number;
  size: number;
  velocityX: Float32Array;
  velocityY: Float32Array;
  nextVelocityX: Float32Array;
  nextVelocityY: Float32Array;
  dyeR: Float32Array;
  dyeG: Float32Array;
  dyeB: Float32Array;
  nextDyeR: Float32Array;
  nextDyeG: Float32Array;
  nextDyeB: Float32Array;
  pressure: Float32Array;
  nextPressure: Float32Array;
  divergence: Float32Array;
  curl: Float32Array;
  imageData: ImageData;
}

interface PointerState {
  id: number;
  type: string;
  x: number;
  y: number;
  prevX: number;
  prevY: number;
  down: boolean;
  hue: number;
  seed: number;
}

const canvasElement = document.getElementById("flowlight-canvas");
const resetButtonElement = document.getElementById("reset-field");

if (
  !(canvasElement instanceof HTMLCanvasElement) ||
  !(resetButtonElement instanceof HTMLButtonElement)
) {
  throw new Error("Unable to find the Flowlight UI elements.");
}

const canvas = canvasElement;
const resetButton = resetButtonElement;
const contextCandidate = canvas.getContext("2d", { alpha: false });
const fluidCanvas = document.createElement("canvas");
const fluidContextCandidate = fluidCanvas.getContext("2d");

if (!contextCandidate || !fluidContextCandidate) {
  throw new Error("Unable to initialize the Flowlight canvas.");
}

const context: CanvasRenderingContext2D = contextCandidate;
const fluidContext: CanvasRenderingContext2D = fluidContextCandidate;

const viewport: Viewport = { width: 0, height: 0, dpr: 1, cellSize: 1 };
const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
const field: FluidField = {
  width: 0,
  height: 0,
  size: 0,
  velocityX: new Float32Array(0),
  velocityY: new Float32Array(0),
  nextVelocityX: new Float32Array(0),
  nextVelocityY: new Float32Array(0),
  dyeR: new Float32Array(0),
  dyeG: new Float32Array(0),
  dyeB: new Float32Array(0),
  nextDyeR: new Float32Array(0),
  nextDyeG: new Float32Array(0),
  nextDyeB: new Float32Array(0),
  pressure: new Float32Array(0),
  nextPressure: new Float32Array(0),
  divergence: new Float32Array(0),
  curl: new Float32Array(0),
  imageData: fluidContext.createImageData(1, 1),
};
const pointers = new Map<number, PointerState>();
let lastFrame = performance.now();
let ambientCooldown = 0;

const clamp = (value: number, min: number, max: number): number =>
  Math.min(Math.max(value, min), max);
const indexOf = (x: number, y: number): number => x + y * field.width;

function swap(keyA: FluidBufferKey, keyB: FluidBufferKey) {
  const next = field[keyA];
  field[keyA] = field[keyB];
  field[keyB] = next;
}

function allocateField() {
  field.size = field.width * field.height;
  field.velocityX = new Float32Array(field.size);
  field.velocityY = new Float32Array(field.size);
  field.nextVelocityX = new Float32Array(field.size);
  field.nextVelocityY = new Float32Array(field.size);
  field.dyeR = new Float32Array(field.size);
  field.dyeG = new Float32Array(field.size);
  field.dyeB = new Float32Array(field.size);
  field.nextDyeR = new Float32Array(field.size);
  field.nextDyeG = new Float32Array(field.size);
  field.nextDyeB = new Float32Array(field.size);
  field.pressure = new Float32Array(field.size);
  field.nextPressure = new Float32Array(field.size);
  field.divergence = new Float32Array(field.size);
  field.curl = new Float32Array(field.size);
  field.imageData = fluidContext.createImageData(field.width, field.height);
}

function resize() {
  viewport.width = window.innerWidth;
  viewport.height = window.innerHeight;
  viewport.dpr = Math.min(window.devicePixelRatio || 1, 2);

  canvas.width = Math.floor(viewport.width * viewport.dpr);
  canvas.height = Math.floor(viewport.height * viewport.dpr);
  canvas.style.width = `${viewport.width}px`;
  canvas.style.height = `${viewport.height}px`;

  context.setTransform(viewport.dpr, 0, 0, viewport.dpr, 0, 0);

  const densityBias = reducedMotion.matches ? 0.88 : 1;
  const cellSize = clamp(
    Math.min(viewport.width, viewport.height) / (170 * viewport.dpr * densityBias),
    2.2,
    5,
  );

  viewport.cellSize = cellSize;
  field.width = Math.max(180, Math.round(viewport.width / cellSize));
  field.height = Math.max(110, Math.round(viewport.height / cellSize));
  fluidCanvas.width = field.width;
  fluidCanvas.height = field.height;
  allocateField();
  resetField();
}

function sample(fieldBuffer: Float32Array, x: number, y: number): number {
  const sampleX = clamp(x, 0, field.width - 1);
  const sampleY = clamp(y, 0, field.height - 1);
  const x0 = Math.floor(sampleX);
  const y0 = Math.floor(sampleY);
  const x1 = Math.min(x0 + 1, field.width - 1);
  const y1 = Math.min(y0 + 1, field.height - 1);
  const tx = sampleX - x0;
  const ty = sampleY - y0;

  const top =
    fieldBuffer[indexOf(x0, y0)] * (1 - tx) +
    fieldBuffer[indexOf(x1, y0)] * tx;
  const bottom =
    fieldBuffer[indexOf(x0, y1)] * (1 - tx) +
    fieldBuffer[indexOf(x1, y1)] * tx;

  return top * (1 - ty) + bottom * ty;
}

function hslToRgb(
  hue: number,
  saturation: number,
  lightness: number,
): [number, number, number] {
  if (saturation === 0) {
    return [lightness, lightness, lightness];
  }

  const hueToChannel = (p: number, q: number, t: number): number => {
    let channel = t;

    if (channel < 0) {
      channel += 1;
    }

    if (channel > 1) {
      channel -= 1;
    }

    if (channel < 1 / 6) {
      return p + (q - p) * 6 * channel;
    }

    if (channel < 1 / 2) {
      return q;
    }

    if (channel < 2 / 3) {
      return p + (q - p) * (2 / 3 - channel) * 6;
    }

    return p;
  };

  const q =
    lightness < 0.5
      ? lightness * (1 + saturation)
      : lightness + saturation - lightness * saturation;
  const p = 2 * lightness - q;

  return [
    hueToChannel(p, q, hue + 1 / 3),
    hueToChannel(p, q, hue),
    hueToChannel(p, q, hue - 1 / 3),
  ];
}

function splat(
  screenX: number,
  screenY: number,
  velocityX: number,
  velocityY: number,
  energy: number,
  hue: number,
  strength = 1,
) {
  const centerX = (screenX / viewport.width) * (field.width - 1);
  const centerY = (screenY / viewport.height) * (field.height - 1);
  const radius = clamp((18 + energy * 0.38) / viewport.cellSize, 4, 24);
  const radiusSquared = radius * radius;
  const [red, green, blue] = hslToRgb(
    (((hue % 360) + 360) % 360) / 360,
    0.9,
    0.56,
  );
  const minX = Math.max(1, Math.floor(centerX - radius * 2.2));
  const maxX = Math.min(field.width - 2, Math.ceil(centerX + radius * 2.2));
  const minY = Math.max(1, Math.floor(centerY - radius * 2.2));
  const maxY = Math.min(field.height - 2, Math.ceil(centerY + radius * 2.2));

  for (let y = minY; y <= maxY; y += 1) {
    for (let x = minX; x <= maxX; x += 1) {
      const offsetX = x - centerX;
      const offsetY = y - centerY;
      const distanceSquared = offsetX * offsetX + offsetY * offsetY;
      const influence = Math.exp(-distanceSquared / radiusSquared);

      if (influence < 0.0035) {
        continue;
      }

      const index = indexOf(x, y);
      field.velocityX[index] += velocityX * influence * 0.78;
      field.velocityY[index] += velocityY * influence * 0.78;
      field.dyeR[index] = Math.min(
        field.dyeR[index] + red * influence * (0.38 + strength * 0.22),
        3.2,
      );
      field.dyeG[index] = Math.min(
        field.dyeG[index] + green * influence * (0.38 + strength * 0.22),
        3.2,
      );
      field.dyeB[index] = Math.min(
        field.dyeB[index] + blue * influence * (0.38 + strength * 0.22),
        3.2,
      );
    }
  }
}

function resetField() {
  const keys: FluidBufferKey[] = [
    "velocityX",
    "velocityY",
    "nextVelocityX",
    "nextVelocityY",
    "dyeR",
    "dyeG",
    "dyeB",
    "nextDyeR",
    "nextDyeG",
    "nextDyeB",
    "pressure",
    "nextPressure",
    "divergence",
    "curl",
  ];

  keys.forEach((key) => field[key].fill(0));

  for (let bloom = 0; bloom < 7; bloom += 1) {
    const seedX = viewport.width * (0.12 + Math.random() * 0.76);
    const seedY = viewport.height * (0.14 + Math.random() * 0.72);
    const impulseX = (Math.random() - 0.5) * 3.5;
    const impulseY = (Math.random() - 0.5) * 3.5;
    splat(
      seedX,
      seedY,
      impulseX,
      impulseY,
      14 + Math.random() * 8,
      90 + Math.random() * 250,
      0.72,
    );
  }
}

function advectScalar(
  target: Float32Array,
  source: Float32Array,
  dissipation: number,
  dt: number,
) {
  for (let y = 0; y < field.height; y += 1) {
    for (let x = 0; x < field.width; x += 1) {
      const index = indexOf(x, y);
      const backX = x - dt * field.velocityX[index];
      const backY = y - dt * field.velocityY[index];
      target[index] = sample(source, backX, backY) * dissipation;
    }
  }
}

function advectVelocity(dt: number) {
  for (let y = 0; y < field.height; y += 1) {
    for (let x = 0; x < field.width; x += 1) {
      const index = indexOf(x, y);
      const backX = x - dt * field.velocityX[index];
      const backY = y - dt * field.velocityY[index];
      field.nextVelocityX[index] = sample(field.velocityX, backX, backY) * 0.996;
      field.nextVelocityY[index] = sample(field.velocityY, backX, backY) * 0.996;
    }
  }

  swap("velocityX", "nextVelocityX");
  swap("velocityY", "nextVelocityY");
}

function computeCurl() {
  for (let y = 1; y < field.height - 1; y += 1) {
    for (let x = 1; x < field.width - 1; x += 1) {
      const left = field.velocityY[indexOf(x - 1, y)];
      const right = field.velocityY[indexOf(x + 1, y)];
      const top = field.velocityX[indexOf(x, y - 1)];
      const bottom = field.velocityX[indexOf(x, y + 1)];
      field.curl[indexOf(x, y)] = 0.5 * (right - left - bottom + top);
    }
  }
}

function applyVorticity(dt: number) {
  computeCurl();

  for (let y = 2; y < field.height - 2; y += 1) {
    for (let x = 2; x < field.width - 2; x += 1) {
      const index = indexOf(x, y);
      const forceX =
        Math.abs(field.curl[indexOf(x, y + 1)]) -
        Math.abs(field.curl[indexOf(x, y - 1)]);
      const forceY =
        Math.abs(field.curl[indexOf(x - 1, y)]) -
        Math.abs(field.curl[indexOf(x + 1, y)]);
      const length = Math.hypot(forceX, forceY) + 0.00001;
      const normalizedX = forceX / length;
      const normalizedY = forceY / length;
      const vorticity = field.curl[index];
      field.velocityX[index] += normalizedX * vorticity * dt * 18;
      field.velocityY[index] += normalizedY * vorticity * dt * 18;
    }
  }
}

function projectPressure() {
  for (let y = 1; y < field.height - 1; y += 1) {
    for (let x = 1; x < field.width - 1; x += 1) {
      const index = indexOf(x, y);
      const divergence =
        -0.5 *
        (field.velocityX[indexOf(x + 1, y)] -
          field.velocityX[indexOf(x - 1, y)] +
          field.velocityY[indexOf(x, y + 1)] -
          field.velocityY[indexOf(x, y - 1)]);
      field.divergence[index] = divergence;
      field.pressure[index] = 0;
    }
  }

  for (let iteration = 0; iteration < 14; iteration += 1) {
    for (let y = 1; y < field.height - 1; y += 1) {
      for (let x = 1; x < field.width - 1; x += 1) {
        const index = indexOf(x, y);
        field.nextPressure[index] =
          (field.divergence[index] +
            field.pressure[indexOf(x + 1, y)] +
            field.pressure[indexOf(x - 1, y)] +
            field.pressure[indexOf(x, y + 1)] +
            field.pressure[indexOf(x, y - 1)]) *
          0.25;
      }
    }

    swap("pressure", "nextPressure");
  }

  for (let y = 1; y < field.height - 1; y += 1) {
    for (let x = 1; x < field.width - 1; x += 1) {
      const index = indexOf(x, y);
      field.velocityX[index] -=
        0.5 *
        (field.pressure[indexOf(x + 1, y)] - field.pressure[indexOf(x - 1, y)]);
      field.velocityY[index] -=
        0.5 *
        (field.pressure[indexOf(x, y + 1)] - field.pressure[indexOf(x, y - 1)]);
    }
  }
}

function applyPointerForces(now: number, dt: number) {
  pointers.forEach((pointer) => {
    const deltaX = pointer.x - pointer.prevX;
    const deltaY = pointer.y - pointer.prevY;
    const distance = Math.hypot(deltaX, deltaY);
    const orbit = now * 0.0012 + pointer.seed;
    const driftX = Math.cos(orbit) * (pointer.down ? 0.7 : 0.25);
    const driftY = Math.sin(orbit * 1.2) * (pointer.down ? 0.7 : 0.25);
    const impulseX = (deltaX / viewport.width) * field.width * 26 + driftX;
    const impulseY = (deltaY / viewport.height) * field.height * 26 + driftY;
    const energy = Math.max(distance, pointer.down ? 7 : 0);

    if (energy > 0.1) {
      pointer.hue = (pointer.hue + dt * 38 + distance * 0.42) % 360;
      splat(
        pointer.x,
        pointer.y,
        impulseX,
        impulseY,
        energy,
        pointer.hue,
        pointer.down ? 1.12 : 0.82,
      );
    }

    pointer.prevX = pointer.x;
    pointer.prevY = pointer.y;

    if (!pointer.down && pointer.type !== "mouse") {
      pointers.delete(pointer.id);
    }
  });
}

function spawnAmbient(now: number, dt: number) {
  if (reducedMotion.matches) {
    return;
  }

  ambientCooldown -= dt;

  if (ambientCooldown > 0) {
    return;
  }

  ambientCooldown = 0.06 + Math.random() * 0.12;

  const orbit = now * 0.00014;
  const x =
    viewport.width *
    (0.5 + Math.sin(orbit * 2.2) * 0.26 + (Math.random() - 0.5) * 0.04);
  const y =
    viewport.height *
    (0.5 + Math.cos(orbit * 1.7) * 0.22 + (Math.random() - 0.5) * 0.05);
  const impulseX = Math.cos(orbit * 3.1) * 2.5;
  const impulseY = Math.sin(orbit * 2.4) * 2.5;
  const hue = (now * 0.02 + 90 + Math.random() * 110) % 360;
  splat(x, y, impulseX, impulseY, 18 + Math.random() * 12, hue, 0.8);
}

function step(dt: number, now: number) {
  applyPointerForces(now, dt);
  spawnAmbient(now, dt);
  applyVorticity(dt);
  projectPressure();
  advectVelocity(dt);
  projectPressure();
  advectScalar(field.nextDyeR, field.dyeR, 0.988, dt);
  advectScalar(field.nextDyeG, field.dyeG, 0.988, dt);
  advectScalar(field.nextDyeB, field.dyeB, 0.988, dt);
  swap("dyeR", "nextDyeR");
  swap("dyeG", "nextDyeG");
  swap("dyeB", "nextDyeB");
}

function renderField(now: number) {
  const pixels = field.imageData.data;

  for (let index = 0; index < field.size; index += 1) {
    const base = index * 4;
    const speed = Math.hypot(field.velocityX[index], field.velocityY[index]) * 0.06;
    const baseRed = 1 - Math.exp(-field.dyeR[index] * 0.84);
    const baseGreen = 1 - Math.exp(-field.dyeG[index] * 0.84);
    const baseBlue = 1 - Math.exp(-field.dyeB[index] * 0.84);
    const tone = 1 / (1 + Math.max(0, baseRed + baseGreen + baseBlue - 0.9) * 2.4);
    const red = baseRed * tone;
    const green = baseGreen * tone;
    const blue = baseBlue * tone;
    const shimmer = 0.07 + 0.05 * Math.sin(now * 0.0015 + index * 0.0014);
    const intensity = Math.max(red, green, blue);

    pixels[base] = clamp((red + speed * 0.18 + shimmer * 0.08) * 255, 0, 255);
    pixels[base + 1] = clamp(
      (green + speed * 0.14 + shimmer * 0.06) * 255,
      0,
      255,
    );
    pixels[base + 2] = clamp((blue + speed * 0.22 + shimmer * 0.1) * 255, 0, 255);
    pixels[base + 3] = clamp((intensity * 0.5 + speed * 0.18) * 255, 0, 180);
  }

  fluidContext.putImageData(field.imageData, 0, 0);
}

function drawBackdrop(now: number) {
  const width = viewport.width;
  const height = viewport.height;
  const gradient = context.createLinearGradient(0, 0, width, height);
  gradient.addColorStop(0, "#03060b");
  gradient.addColorStop(0.55, "#07111a");
  gradient.addColorStop(1, "#020308");
  context.fillStyle = gradient;
  context.fillRect(0, 0, width, height);

  const haloX = width * (0.48 + Math.sin(now * 0.00012) * 0.09);
  const haloY = height * (0.42 + Math.cos(now * 0.00015) * 0.08);
  const halo = context.createRadialGradient(haloX, haloY, 0, haloX, haloY, width * 0.52);
  halo.addColorStop(0, "rgba(42, 120, 88, 0.2)");
  halo.addColorStop(0.38, "rgba(12, 24, 33, 0.1)");
  halo.addColorStop(1, "rgba(0, 0, 0, 0)");
  context.fillStyle = halo;
  context.fillRect(0, 0, width, height);
}

function render(now: number) {
  drawBackdrop(now);
  renderField(now);

  const width = viewport.width;
  const height = viewport.height;
  const blurLarge = clamp(viewport.cellSize * 2, 5, 10);
  const blurSmall = clamp(viewport.cellSize * 0.9, 2, 4);

  context.save();
  context.globalCompositeOperation = "screen";
  context.globalAlpha = 0.52;
  context.filter = `blur(${blurLarge}px) saturate(138%) brightness(96%)`;
  context.drawImage(fluidCanvas, 0, 0, width, height);
  context.globalAlpha = 0.74;
  context.filter = `blur(${blurSmall}px) saturate(126%) brightness(94%)`;
  context.drawImage(fluidCanvas, 0, 0, width, height);
  context.restore();

  context.save();
  context.globalAlpha = 0.28;
  context.globalCompositeOperation = "source-over";
  context.filter = "none";
  context.drawImage(fluidCanvas, 0, 0, width, height);
  context.restore();
}

function frame(now: number) {
  const dt = clamp((now - lastFrame) / 1000, 0.008, 0.024);
  lastFrame = now;
  step(dt, now);
  render(now);
  window.requestAnimationFrame(frame);
}

function updatePointer(event: PointerEvent): PointerState {
  const existing = pointers.get(event.pointerId);

  if (existing) {
    existing.x = event.clientX;
    existing.y = event.clientY;

    if (event.pointerType !== "mouse") {
      existing.down = true;
    }

    return existing;
  }

  const pointer: PointerState = {
    id: event.pointerId,
    type: event.pointerType,
    x: event.clientX,
    y: event.clientY,
    prevX: event.clientX,
    prevY: event.clientY,
    down: event.pointerType === "mouse" || event.buttons > 0,
    hue: Math.random() * 360,
    seed: Math.random() * Math.PI * 2,
  };

  pointers.set(event.pointerId, pointer);
  return pointer;
}

function handlePointerDown(event: PointerEvent) {
  event.preventDefault();
  const pointer = updatePointer(event);
  pointer.down = true;
  pointer.prevX = pointer.x = event.clientX;
  pointer.prevY = pointer.y = event.clientY;
  splat(pointer.x, pointer.y, 0, 0, 24, pointer.hue, 1.1);
}

function handlePointerMove(event: PointerEvent) {
  if (event.pointerType !== "mouse" && !pointers.has(event.pointerId)) {
    return;
  }

  updatePointer(event);
}

function handlePointerUp(event: PointerEvent) {
  const pointer = pointers.get(event.pointerId);

  if (!pointer) {
    return;
  }

  if (pointer.type === "mouse") {
    pointer.down = false;
    return;
  }

  pointers.delete(event.pointerId);
}

function handlePointerLeave(event: PointerEvent) {
  if (event.pointerType === "mouse") {
    pointers.delete(event.pointerId);
  }
}

canvas.addEventListener("pointerdown", handlePointerDown, { passive: false });
canvas.addEventListener("pointermove", handlePointerMove, { passive: true });
canvas.addEventListener("pointerup", handlePointerUp, { passive: true });
canvas.addEventListener("pointercancel", handlePointerUp, { passive: true });
canvas.addEventListener("pointerleave", handlePointerLeave, { passive: true });
resetButton.addEventListener("click", resetField);
window.addEventListener("resize", resize, { passive: true });
resize();
window.requestAnimationFrame(frame);
