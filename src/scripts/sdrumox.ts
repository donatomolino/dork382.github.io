export {};

type MoveKey = "succhiapolli" | "finocchio" | "mystery" | "lavica";
type ItemKey = "vitaBonus" | "crtShield" | "turboPatch";
type FloaterKind = "damage" | "heal" | "charge";
type BattleOutcome = "idle" | "win" | "lose";

interface EntityState {
  hp: number;
  maxHp: number;
  displayHp: number;
  hurtUntil: number;
  healUntil: number;
}

interface BattleState {
  round: number;
  bossCharge: number;
  displayCharge: number;
  locked: boolean;
  ended: boolean;
  outcome: BattleOutcome;
  player: EntityState & {
    bonusLives: number;
    shieldTurns: number;
  };
  boss: EntityState & {
    burnTurns: number;
    stunned: boolean;
    damageDownTurns: number;
  };
  items: Record<ItemKey, number>;
  log: string[];
  soundEnabled: boolean;
}

type Palette = Record<string, string>;

const PLAYER_SPRITE = [
  ".......nn.......",
  "......nnnn......",
  ".....nssssn.....",
  ".....nssssn.....",
  "......swws......",
  ".....bbbrr......",
  "....bbbbrrr.....",
  "....bbyyyrr.....",
  ".....byyyr......",
  ".....bbyyb......",
  ".....d..dd......",
  "....dd..ddd.....",
  "....d....dd.....",
  "...dd....ddd....",
  "................",
  "................",
];

const BOSS_SPRITE = [
  "....................",
  "........rr..........",
  ".....rrrppprrr......",
  "....rpppppppppr.....",
  "...rpppccppccpppr....",
  "..rppppppppppppppr...",
  ".rpppppyyyyyypppppr..",
  ".rppppyyyyyyyyppppr..",
  ".rpppbyyyyyyyy bpppr.",
  ".rpppbyyyyyyyy bpppr.",
  ".rppppyyyyyyyyppppr..",
  ".rrpppooooooooppprr..",
  "..rppooorrrrooopppr..",
  "..rppoo......oooppr..",
  "...rppo......ooppr...",
  "....rrr......rrr.....",
  "......d......d.......",
  ".....dd......dd......",
  "....................",
  "....................",
].map((row) => row.replaceAll(" ", "o"));

const PLAYER_PALETTE: Palette = {
  n: "#2b1c18",
  s: "#f2d4b0",
  w: "#f8f5ef",
  b: "#5ec4ff",
  y: "#ffd166",
  r: "#ff6b6b",
  d: "#202835",
};

const BOSS_PALETTE: Palette = {
  r: "#4d1022",
  p: "#8b1e4f",
  c: "#ffe082",
  y: "#ff8e4d",
  b: "#c44368",
  o: "#361016",
  d: "#1b0f12",
};

const playerLines = [
  "Pixel Don entra in arena con la faccia di chi non ha alcuna intenzione di perdere.",
  "Il boss carica il suo meter. La tua finestra è corta.",
  "Un colpo sbagliato qui si paga in round persi.",
];

function expectInstance<T extends Element>(
  value: Element | null,
  ctor: { new (): T; prototype: T },
  message: string,
): T {
  if (!(value instanceof ctor)) {
    throw new Error(message);
  }

  return value;
}

const moveButtons = Array.from(document.querySelectorAll<HTMLButtonElement>("[data-move]"));
const itemButtons = Array.from(document.querySelectorAll<HTMLButtonElement>("[data-item]"));
const canvas = expectInstance(
  document.getElementById("sdrumox-canvas"),
  HTMLCanvasElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const battleStage = expectInstance(
  document.getElementById("battle-stage"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const floaters = expectInstance(
  document.getElementById("floaters"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const stageMessage = expectInstance(
  document.getElementById("stage-message"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const playerHpFill = expectInstance(
  document.getElementById("player-hp-fill"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const playerHpText = expectInstance(
  document.getElementById("player-hp-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bossHpFill = expectInstance(
  document.getElementById("boss-hp-fill"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bossHpText = expectInstance(
  document.getElementById("boss-hp-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bossChargeFill = expectInstance(
  document.getElementById("boss-charge-fill"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bossChargeText = expectInstance(
  document.getElementById("boss-charge-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bossPhaseText = expectInstance(
  document.getElementById("boss-phase-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const roundText = expectInstance(
  document.getElementById("round-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const statusBadge = expectInstance(
  document.getElementById("status-badge"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const bonusLifeText = expectInstance(
  document.getElementById("bonus-life-text"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const dialogueLog = expectInstance(
  document.getElementById("dialogue-log"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const restartButton = expectInstance(
  document.getElementById("restart-battle"),
  HTMLButtonElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const soundToggleButton = expectInstance(
  document.getElementById("sound-toggle"),
  HTMLButtonElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const jumpscare = expectInstance(
  document.getElementById("jumpscare"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const itemVitaBonusCount = expectInstance(
  document.getElementById("item-vitaBonus-count"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const itemCrtShieldCount = expectInstance(
  document.getElementById("item-crtShield-count"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);
const itemTurboPatchCount = expectInstance(
  document.getElementById("item-turboPatch-count"),
  HTMLElement,
  "Unable to initialize the Sdrumox battle UI.",
);

const contextCandidate = canvas.getContext("2d");

if (!contextCandidate) {
  throw new Error("Unable to create the Sdrumox canvas context.");
}

const context: CanvasRenderingContext2D = contextCandidate;
let devicePixelRatioValue = Math.min(window.devicePixelRatio || 1, 2);
let animationFrameId = 0;
let audioContext: AudioContext | null = null;

const state: BattleState = {
  round: 1,
  bossCharge: 0,
  displayCharge: 0,
  locked: false,
  ended: false,
  outcome: "idle",
  player: {
    hp: 128,
    maxHp: 128,
    displayHp: 128,
    hurtUntil: 0,
    healUntil: 0,
    bonusLives: 1,
    shieldTurns: 0,
  },
  boss: {
    hp: 220,
    maxHp: 220,
    displayHp: 220,
    hurtUntil: 0,
    healUntil: 0,
    burnTurns: 0,
    stunned: false,
    damageDownTurns: 0,
  },
  items: {
    vitaBonus: 1,
    crtShield: 2,
    turboPatch: 1,
  },
  log: [],
  soundEnabled: true,
};

const clamp = (value: number, min: number, max: number): number =>
  Math.min(Math.max(value, min), max);

const randomInt = (min: number, max: number): number =>
  Math.floor(Math.random() * (max - min + 1)) + min;

function delay(duration: number) {
  return new Promise<void>((resolve) => {
    window.setTimeout(resolve, duration);
  });
}

function getAudioContext(): AudioContext | null {
  if (!state.soundEnabled) {
    return null;
  }

  if (!audioContext) {
    audioContext = new AudioContext();
  }

  if (audioContext.state === "suspended") {
    void audioContext.resume();
  }

  return audioContext;
}

function playTone(
  frequency: number,
  duration: number,
  type: OscillatorType,
  volume: number,
  slideTo?: number,
  delayTime = 0,
) {
  const currentAudio = getAudioContext();

  if (!currentAudio) {
    return;
  }

  const oscillator = currentAudio.createOscillator();
  const gain = currentAudio.createGain();
  const startAt = currentAudio.currentTime + delayTime;

  oscillator.type = type;
  oscillator.frequency.setValueAtTime(Math.max(1, frequency), startAt);

  if (slideTo) {
    oscillator.frequency.exponentialRampToValueAtTime(Math.max(1, slideTo), startAt + duration);
  }

  gain.gain.setValueAtTime(0.0001, startAt);
  gain.gain.exponentialRampToValueAtTime(volume, startAt + 0.02);
  gain.gain.exponentialRampToValueAtTime(0.0001, startAt + duration);

  oscillator.connect(gain);
  gain.connect(currentAudio.destination);
  oscillator.start(startAt);
  oscillator.stop(startAt + duration + 0.04);
}

function playNoise(duration: number, volume: number) {
  const currentAudio = getAudioContext();

  if (!currentAudio) {
    return;
  }

  const buffer = currentAudio.createBuffer(1, currentAudio.sampleRate * duration, currentAudio.sampleRate);
  const channel = buffer.getChannelData(0);

  for (let index = 0; index < channel.length; index += 1) {
    channel[index] = (Math.random() * 2 - 1) * (1 - index / channel.length);
  }

  const source = currentAudio.createBufferSource();
  const gain = currentAudio.createGain();
  source.buffer = buffer;
  gain.gain.value = volume;
  source.connect(gain);
  gain.connect(currentAudio.destination);
  source.start();
}

function playClickSound() {
  playTone(420, 0.08, "square", 0.03, 560);
}

function playAttackSound() {
  playTone(180, 0.14, "sawtooth", 0.07, 90);
  playNoise(0.09, 0.018);
}

function playHealSound() {
  playTone(440, 0.12, "triangle", 0.04, 660);
  playTone(660, 0.13, "triangle", 0.03, 920, 0.08);
}

function playChargeSound() {
  playTone(110, 0.22, "square", 0.05, 180);
  playTone(220, 0.14, "triangle", 0.025, 330, 0.06);
}

function playBossHitSound() {
  playTone(130, 0.2, "sawtooth", 0.06, 70);
  playNoise(0.12, 0.022);
}

function playJumpscareSound() {
  playNoise(0.24, 0.04);
  playTone(920, 0.12, "square", 0.04, 220);
  playTone(70, 0.22, "sawtooth", 0.05, 40);
}

function playWinSound() {
  playTone(330, 0.14, "triangle", 0.04, 440);
  playTone(440, 0.14, "triangle", 0.04, 660, 0.1);
  playTone(660, 0.18, "triangle", 0.04, 880, 0.22);
}

function playLoseSound() {
  playTone(180, 0.18, "square", 0.05, 120);
  playTone(120, 0.22, "square", 0.04, 75, 0.12);
}

function pushDialogue(message: string) {
  state.log = [...state.log, message].slice(-6);
  dialogueLog.innerHTML = "";

  state.log.forEach((line) => {
    const row = document.createElement("div");
    row.className = "dialogue-line";
    row.textContent = line;
    dialogueLog.append(row);
  });
}

function setStageMessage(message: string) {
  stageMessage.textContent = message;
}

function setStageShake() {
  battleStage.classList.remove("stage-shake");
  void battleStage.offsetWidth;
  battleStage.classList.add("stage-shake");
}

function spawnFloater(text: string, kind: FloaterKind, side: "left" | "right") {
  const floater = document.createElement("div");
  floater.className = `floater ${kind}`;
  floater.textContent = text;
  floater.style.left = side === "left" ? `${randomInt(28, 42)}%` : `${randomInt(60, 78)}%`;
  floater.style.top = `${randomInt(28, 42)}%`;
  floaters.append(floater);
  window.setTimeout(() => floater.remove(), 900);
}

async function tweenNumber(
  from: number,
  to: number,
  duration: number,
  update: (value: number) => void,
) {
  return new Promise<void>((resolve) => {
    const startedAt = performance.now();

    const frame = (now: number) => {
      const progress = clamp((now - startedAt) / duration, 0, 1);
      const eased = 1 - (1 - progress) * (1 - progress);
      update(from + (to - from) * eased);

      if (progress < 1) {
        window.requestAnimationFrame(frame);
        return;
      }

      update(to);
      resolve();
    };

    window.requestAnimationFrame(frame);
  });
}

function phaseLabel(charge: number): string {
  if (charge >= 90) return "Catastrofe imminente";
  if (charge >= 70) return "Overclock feroce";
  if (charge >= 50) return "Boss in furia";
  if (charge >= 30) return "Sdrumox si scalda";
  return "Ingresso Boss";
}

function renderUi() {
  const playerHpRatio = (state.player.displayHp / state.player.maxHp) * 100;
  const bossHpRatio = (state.boss.displayHp / state.boss.maxHp) * 100;

  playerHpFill.style.width = `${clamp(playerHpRatio, 0, 100)}%`;
  bossHpFill.style.width = `${clamp(bossHpRatio, 0, 100)}%`;
  bossChargeFill.style.width = `${clamp(state.displayCharge, 0, 100)}%`;

  playerHpText.textContent = `${Math.round(state.player.displayHp)} / ${state.player.maxHp}`;
  bossHpText.textContent = `${Math.round(state.boss.displayHp)} / ${state.boss.maxHp}`;
  bossChargeText.textContent = `${Math.round(state.displayCharge)}%`;
  bossPhaseText.textContent = `Sdrumox ${state.bossCharge}%`;
  roundText.textContent = `Round ${state.round}`;
  bonusLifeText.textContent = `Bonus Vita x${state.player.bonusLives}`;
  statusBadge.textContent =
    state.outcome === "win"
      ? "Sdrumox KO"
      : state.outcome === "lose"
        ? "Run finita"
        : phaseLabel(state.bossCharge);
  soundToggleButton.textContent = state.soundEnabled ? "Audio ON" : "Audio OFF";

  itemVitaBonusCount.textContent = `x${state.items.vitaBonus}`;
  itemCrtShieldCount.textContent = `x${state.items.crtShield}`;
  itemTurboPatchCount.textContent = `x${state.items.turboPatch}`;

  moveButtons.forEach((button) => {
    button.disabled = state.locked || state.ended;
  });

  itemButtons.forEach((button) => {
    const key = button.dataset.item as ItemKey | undefined;
    button.disabled = state.locked || state.ended || !key || state.items[key] <= 0;
  });
}

function drawSprite(
  sprite: string[],
  palette: Palette,
  x: number,
  y: number,
  pixelSize: number,
  flip = false,
) {
  sprite.forEach((row, rowIndex) => {
    Array.from(row).forEach((token, columnIndex) => {
      if (token === ".") {
        return;
      }

      context.fillStyle = palette[token] ?? "#ffffff";
      const drawX = flip
        ? x + (row.length - 1 - columnIndex) * pixelSize
        : x + columnIndex * pixelSize;
      context.fillRect(drawX, y + rowIndex * pixelSize, pixelSize, pixelSize);
    });
  });
}

function drawBackground(width: number, height: number, now: number) {
  context.clearRect(0, 0, width, height);
  context.imageSmoothingEnabled = false;

  context.fillStyle = "#111827";
  context.fillRect(0, 0, width, height);

  context.fillStyle = "#22304a";
  context.fillRect(0, 0, width, height * 0.52);

  context.fillStyle = "#2b3c59";
  for (let index = 0; index < 14; index += 1) {
    context.fillRect(index * 64, 24 + (index % 3) * 18, 34, 12);
  }

  context.fillStyle = "#ffcf6c";
  context.fillRect(width - 168, 42, 72, 72);

  context.fillStyle = "#314864";
  context.beginPath();
  context.moveTo(0, height * 0.58);
  context.lineTo(width * 0.18, height * 0.34);
  context.lineTo(width * 0.34, height * 0.58);
  context.closePath();
  context.fill();

  context.beginPath();
  context.moveTo(width * 0.18, height * 0.58);
  context.lineTo(width * 0.34, height * 0.26);
  context.lineTo(width * 0.48, height * 0.58);
  context.closePath();
  context.fill();

  context.fillStyle = "#4d6043";
  context.fillRect(0, height * 0.58, width, height * 0.42);

  context.fillStyle = "#3b2a22";
  context.fillRect(0, height * 0.72, width, height * 0.28);

  context.fillStyle = "#6a4e3d";
  for (let index = 0; index < 24; index += 1) {
    context.fillRect(index * 48, height * 0.78 + ((index + 1) % 2) * 10, 20, 44);
  }

  const chargeGlow = Math.max(0, state.displayCharge / 100);
  context.fillStyle = `rgba(185, 140, 255, ${0.08 + chargeGlow * 0.24})`;
  context.fillRect(width * 0.58, height * 0.18, width * 0.28, height * 0.36);

  context.fillStyle = "#1d2430";
  context.fillRect(width * 0.07, height * 0.64, width * 0.24, 18);
  context.fillRect(width * 0.64, height * 0.54, width * 0.2, 18);

  context.fillStyle = "#6fd8ff";
  const sparkleX = width * 0.68 + Math.sin(now * 0.003) * 16;
  const sparkleY = height * 0.18 + Math.cos(now * 0.0035) * 10;
  context.fillRect(sparkleX, sparkleY, 8, 8);
  context.fillRect(sparkleX - 12, sparkleY + 6, 4, 4);
  context.fillRect(sparkleX + 12, sparkleY - 6, 4, 4);
}

function drawScene(now: number) {
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  const playerPixelSize = Math.max(4, Math.floor(width / 82));
  const bossPixelSize = Math.max(4, Math.floor(width / 98));
  const playerBob = Math.sin(now * 0.004) * 2;
  const bossBob = Math.sin(now * 0.003 + 2) * 3;

  drawBackground(width, height, now);

  const playerScale = now < state.player.hurtUntil ? -4 : now < state.player.healUntil ? 3 : 0;
  const bossScale = now < state.boss.hurtUntil ? -4 : now < state.boss.healUntil ? 3 : 0;

  drawSprite(PLAYER_SPRITE, PLAYER_PALETTE, width * 0.18, height * 0.46 + playerBob + playerScale, playerPixelSize, false);
  drawSprite(BOSS_SPRITE, BOSS_PALETTE, width * 0.62, height * 0.29 + bossBob + bossScale, bossPixelSize, true);

  if (!state.ended) {
    context.fillStyle = "rgba(255,255,255,0.12)";
    context.fillRect(width * 0.62, height * 0.24, width * 0.01 + state.displayCharge * 1.6, 6);
  }

  animationFrameId = window.requestAnimationFrame(drawScene);
}

function resizeCanvas() {
  devicePixelRatioValue = Math.min(window.devicePixelRatio || 1, 2);
  const width = Math.floor(canvas.clientWidth * devicePixelRatioValue);
  const height = Math.floor(canvas.clientHeight * devicePixelRatioValue);

  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }

  context.setTransform(devicePixelRatioValue, 0, 0, devicePixelRatioValue, 0, 0);
}

async function animateEntityHp(target: "player" | "boss", nextHp: number) {
  const entity = state[target];
  const start = entity.displayHp;
  entity.hp = clamp(nextHp, 0, entity.maxHp);
  await tweenNumber(start, entity.hp, 380, (value) => {
    entity.displayHp = value;
    renderUi();
  });
}

async function animateCharge(nextCharge: number) {
  const start = state.displayCharge;
  state.bossCharge = clamp(nextCharge, 0, 100);
  await tweenNumber(start, state.bossCharge, 420, (value) => {
    state.displayCharge = value;
    renderUi();
  });
}

async function applyDamage(target: "player" | "boss", amount: number) {
  const entity = state[target];
  const nextHp = entity.hp - amount;
  spawnFloater(`-${amount}`, "damage", target === "player" ? "left" : "right");
  entity.hurtUntil = performance.now() + 220;
  setStageShake();
  if (target === "player") {
    playBossHitSound();
  } else {
    playAttackSound();
  }
  await animateEntityHp(target, nextHp);
}

async function applyHeal(target: "player" | "boss", amount: number) {
  const entity = state[target];
  const nextHp = entity.hp + amount;
  spawnFloater(`+${amount}`, "heal", target === "player" ? "left" : "right");
  entity.healUntil = performance.now() + 220;
  playHealSound();
  await animateEntityHp(target, nextHp);
}

async function changeCharge(delta: number) {
  if (delta === 0) {
    return;
  }

  const nextCharge = clamp(state.bossCharge + delta, 0, 100);
  spawnFloater(`${delta > 0 ? "+" : ""}${delta}%`, "charge", "right");
  playChargeSound();
  await animateCharge(nextCharge);
}

async function triggerJumpscare() {
  jumpscare.classList.add("active");
  playJumpscareSound();
  await delay(260);
  jumpscare.classList.remove("active");
  await delay(90);
}

async function useMove(move: MoveKey) {
  switch (move) {
    case "succhiapolli": {
      const damage = randomInt(18, 24);
      const heal = randomInt(12, 18);
      pushDialogue(`Succhiapolli aggancia il boss e gli strappa ${damage} HP.`);
      setStageMessage("Succhiapolli in corso...");
      await applyDamage("boss", damage);
      await delay(120);
      await applyHeal("player", heal);
      pushDialogue(`Ti rimetti in piedi con ${heal} HP rubati a caldo.`);
      break;
    }

    case "finocchio": {
      const damage = randomInt(14, 20);
      pushDialogue(`Finocchio colpisce di taglio e sballa il meter di Sdrumox.`);
      setStageMessage("Finocchio: controllo del charge.");
      state.boss.damageDownTurns = 1;
      await applyDamage("boss", damage);
      await delay(120);
      await changeCharge(-10);
      pushDialogue("La prossima botta del boss arriva smorzata.");
      break;
    }

    case "mystery": {
      const damage = randomInt(16, 30);
      pushDialogue("??? lacera lo schermo. Qualcosa entra nel match.");
      setStageMessage("Jumpscare in arrivo.");
      await triggerJumpscare();
      await applyDamage("boss", damage);

      if (Math.random() < 0.55) {
        state.boss.stunned = true;
        pushDialogue("Sdrumox resta in tilt e salta il turno.");
      } else {
        const backlash = randomInt(6, 10);
        pushDialogue(`Il jumpscare si ritorce: ti becchi ${backlash} di rinculo.`);
        await delay(100);
        await applyDamage("player", backlash);
      }
      break;
    }

    case "lavica": {
      const damage = randomInt(26, 34);
      pushDialogue(`La Vica esplode in faccia al boss per ${damage} HP.`);
      setStageMessage("La Vica incendia l'arena.");
      state.boss.burnTurns = Math.max(state.boss.burnTurns, 2);
      await applyDamage("boss", damage);
      pushDialogue("Sdrumox prende fuoco e continuerà a perdere HP.");
      break;
    }
  }
}

async function useItem(item: ItemKey) {
  state.items[item] -= 1;
  playClickSound();

  switch (item) {
    case "vitaBonus": {
      state.player.bonusLives += 1;
      pushDialogue("Vita Bonus attivata: guadagni una 1-UP supplementare.");
      setStageMessage("Vita Bonus online.");
      await applyHeal("player", 18);
      break;
    }

    case "crtShield": {
      state.player.shieldTurns += 1;
      pushDialogue("Schermo CRT attivo: il prossimo colpo viene quasi annullato.");
      setStageMessage("Scudo CRT armato.");
      break;
    }

    case "turboPatch": {
      pushDialogue("Turbo Patch: il meter del boss perde terreno.");
      setStageMessage("Turbo Patch applicata.");
      await changeCharge(-10);
      await applyHeal("player", 8);
      break;
    }
  }
}

async function reviveIfPossible() {
  if (state.player.hp > 0 || state.player.bonusLives <= 0 || state.ended) {
    return false;
  }

  state.player.bonusLives -= 1;
  pushDialogue("La Bonus Vita scatta al millimetro: torni in piedi per un altro round.");
  setStageMessage("1-UP attiva.");
  await applyHeal("player", 54);
  return true;
}

async function bossTurn() {
  if (state.boss.stunned) {
    state.boss.stunned = false;
    pushDialogue("Sdrumox è ancora stordito e perde il turno.");
    setStageMessage("Boss stunnato.");
    return;
  }

  const phase = Math.floor(state.bossCharge / 10);
  const rawDamage = randomInt(10 + phase * 2, 16 + phase * 4);
  const weakenedDamage = state.boss.damageDownTurns > 0 ? Math.floor(rawDamage * 0.65) : rawDamage;
  const shieldedDamage = state.player.shieldTurns > 0 ? Math.max(2, Math.floor(weakenedDamage * 0.25)) : weakenedDamage;
  const damage = shieldedDamage;

  if (state.boss.damageDownTurns > 0) {
    state.boss.damageDownTurns -= 1;
  }

  if (state.player.shieldTurns > 0) {
    state.player.shieldTurns -= 1;
    pushDialogue(`Lo Schermo CRT assorbe quasi tutto: Sdrumox passa da ${weakenedDamage} a ${damage}.`);
  } else {
    pushDialogue(`Sdrumox contrattacca e ti scarica addosso ${damage} danni reali.`);
  }

  setStageMessage(`Contrattacco boss: ${damage} danni.`);
  await applyDamage("player", damage);

  if (state.player.hp <= 0) {
    const revived = await reviveIfPossible();

    if (!revived) {
      await loseSequence("Sdrumox ti schiaccia prima che tu riesca a respirare.");
    }
  }
}

async function burnTick() {
  if (state.boss.burnTurns <= 0 || state.boss.hp <= 0) {
    return;
  }

  const damage = randomInt(8, 12);
  state.boss.burnTurns -= 1;
  pushDialogue(`La Vica continua a friggere il boss: ${damage} danni da burn.`);
  setStageMessage("Burn attivo.");
  await applyDamage("boss", damage);
}

async function instantKillSequence() {
  state.locked = true;
  state.ended = true;
  state.outcome = "lose";
  renderUi();
  pushDialogue("Sdrumox raggiunge il 100%. Scatta la kill certa e la run viene cancellata.");
  setStageMessage("100% raggiunto. Fine run.");
  playLoseSound();
  await delay(180);
  state.player.bonusLives = 0;
  await animateEntityHp("player", 0);
}

async function winSequence() {
  state.locked = true;
  state.ended = true;
  state.outcome = "win";
  renderUi();
  pushDialogue("Sdrumox esplode in pixel e detriti. Hai chiuso la run prima del 100%.");
  setStageMessage("Vittoria. Boss cancellato.");
  playWinSound();
}

async function loseSequence(message: string) {
  state.locked = true;
  state.ended = true;
  state.outcome = "lose";
  renderUi();
  pushDialogue(message);
  setStageMessage("Run fallita. Premi Ricomincia.");
  playLoseSound();
}

async function advanceRound() {
  state.round += 1;
  renderUi();
  await changeCharge(10);
  pushDialogue(`Sdrumox sale al ${state.bossCharge}% e il suo livello cresce.`);

  if (state.bossCharge >= 100 && state.boss.hp > 0) {
    await instantKillSequence();
  }
}

async function takeTurn(kind: "move" | "item", key: MoveKey | ItemKey) {
  if (state.locked || state.ended) {
    return;
  }

  getAudioContext();
  state.locked = true;
  renderUi();

  if (kind === "move") {
    await useMove(key as MoveKey);
  } else {
    await useItem(key as ItemKey);
  }

  if (state.boss.hp <= 0) {
    await winSequence();
    return;
  }

  await delay(180);

  await bossTurn();

  if (state.ended) {
    return;
  }

  if (state.boss.hp > 0) {
    await delay(140);
    await burnTick();
  }

  if (state.boss.hp <= 0) {
    await winSequence();
    return;
  }

  await delay(180);
  await advanceRound();

  if (!state.ended) {
    state.locked = false;
    renderUi();
    setStageMessage("Scegli la prossima mossa.");
  }
}

function initializeBattle() {
  state.round = 1;
  state.bossCharge = 0;
  state.displayCharge = 0;
  state.locked = false;
  state.ended = false;
  state.outcome = "idle";
  state.player.hp = 128;
  state.player.maxHp = 128;
  state.player.displayHp = 128;
  state.player.hurtUntil = 0;
  state.player.healUntil = 0;
  state.player.bonusLives = 1;
  state.player.shieldTurns = 0;
  state.boss.hp = 220;
  state.boss.maxHp = 220;
  state.boss.displayHp = 220;
  state.boss.hurtUntil = 0;
  state.boss.healUntil = 0;
  state.boss.burnTurns = 0;
  state.boss.stunned = false;
  state.boss.damageDownTurns = 0;
  state.items.vitaBonus = 1;
  state.items.crtShield = 2;
  state.items.turboPatch = 1;
  state.log = [];

  playerLines.forEach((line) => pushDialogue(line));
  setStageMessage("Sdrumox emerge dal rumore. Fai la prima mossa.");
  renderUi();
}

moveButtons.forEach((button) => {
  button.addEventListener("click", () => {
    const key = button.dataset.move as MoveKey | undefined;

    if (!key) {
      return;
    }

    void takeTurn("move", key);
  });
});

itemButtons.forEach((button) => {
  button.addEventListener("click", () => {
    const key = button.dataset.item as ItemKey | undefined;

    if (!key || state.items[key] <= 0) {
      return;
    }

    void takeTurn("item", key);
  });
});

restartButton.addEventListener("click", () => {
  playClickSound();
  initializeBattle();
});

soundToggleButton.addEventListener("click", () => {
  state.soundEnabled = !state.soundEnabled;

  if (state.soundEnabled) {
    getAudioContext();
    playClickSound();
  }

  renderUi();
});

window.addEventListener("resize", resizeCanvas, { passive: true });

resizeCanvas();
initializeBattle();
animationFrameId = window.requestAnimationFrame(drawScene);

window.addEventListener("beforeunload", () => {
  window.cancelAnimationFrame(animationFrameId);
});