/**
 * Browser QA for the 3000ms Swiss intro. Captures requested timeline frames at
 * 375/800/1440 px in light and dark, checks landing/overflow, reduced motion,
 * blocked storage, daily gating, the blog header, and scrolled navigation.
 * Run in CI against an Astro preview server on localhost:4322.
 */
import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { chromium } from "playwright";

const BASE = "http://127.0.0.1:4322";
const OUT = "artifacts/intro";
await fs.mkdir(OUT, { recursive: true });

const browser = await chromium.launch({ headless: true });
let pass = 0;

async function openContext(width, dark = false) {
  const context = await browser.newContext({
    viewport: { width, height: width === 375 ? 760 : 900 },
    deviceScaleFactor: 1,
    colorScheme: dark ? "dark" : "light",
  });
  if (dark) {
    await context.addInitScript(() => {
      // Ensure dark mode is set before any page bootstrap script can run.
      const apply = () => {
        if (!document.documentElement) return false;
        document.documentElement.dataset.theme = "dark";
        return true;
      };
      if (!apply()) {
        const observer = new MutationObserver(() => {
          if (apply()) observer.disconnect();
        });
        observer.observe(document, { childList: true });
      }
    });
  }
  return context;
}

async function checkNoOverflow(page, label) {
  const metrics = await page.evaluate(() => ({
    scrollWidth: document.documentElement.scrollWidth,
    viewport: document.documentElement.clientWidth,
  }));
  assert.ok(
    metrics.scrollWidth <= metrics.viewport + 2,
    label + ": horizontal overflow " + JSON.stringify(metrics)
  );
}

async function snapshotTimeline(width, dark) {
  const label = (dark ? "dark" : "light") + "-" + width;
  const context = await openContext(width, dark);
  const page = await context.newPage();
  const pageErrors = [];
  page.on("pageerror", (error) => pageErrors.push(error.message));

  await page.goto(BASE + "/?intro", { waitUntil: "domcontentloaded" });
  await page.waitForFunction(
    () => document.documentElement.classList.contains("intro-running"),
    null,
    { timeout: 3500 }
  );

  const surface = await page.evaluate(() => {
    const overlay = document.querySelector("[data-dm-intro-overlay]");
    return getComputedStyle(overlay).backgroundColor;
  });
  assert.equal(
    surface,
    dark ? "rgb(17, 17, 15)" : "rgb(242, 241, 237)",
    label + ": wrong themed intro background"
  );

  const points = [0, .4, .64, .9, 1];
  const started = Date.now();
  let landing;
  for (const fraction of points) {
    const remaining = Math.max(0, Math.round(fraction * 3000 - (Date.now() - started)));
    if (remaining) await page.waitForTimeout(remaining);
    await page.screenshot({ path: path.join(OUT, label + "-" + fraction.toFixed(2) + ".png") });
    await checkNoOverflow(page, label + " @ " + fraction);

    if (fraction === .9) {
      landing = await page.evaluate(() => {
        const source = document.querySelector(".home-header .signal-square");
        const copy = document.querySelector("[data-dm-intro-logo]");
        if (!source || !copy) return null;
        const a = source.getBoundingClientRect();
        const b = copy.getBoundingClientRect();
        return {
          distance: Math.hypot(
            a.left + a.width / 2 - b.left - b.width / 2,
            a.top + a.height / 2 - b.top - b.height / 2
          ),
          sourceWidth: a.width,
          copyWidth: b.width,
        };
      });
    }
  }

  await page.evaluate(() => window.donatoIntro);
  await page.waitForTimeout(100);
  const result = await page.evaluate(() => ({
    overlayRemoved: !document.querySelector("[data-dm-intro-overlay]"),
    htmlClean: !document.documentElement.classList.contains("intro"),
    logoVisible: getComputedStyle(document.querySelector(".home-header .signal-square")).visibility,
    mainTransform: getComputedStyle(document.querySelector("main")).transform,
  }));

  assert.ok(result.overlayRemoved && result.htmlClean, label + ": intro not cleaned up");
  assert.notEqual(result.logoVisible, "hidden", label + ": source mark hidden");
  assert.equal(result.mainTransform, "none", label + ": main transform not restored");
  assert.deepEqual(pageErrors, [], label + ": console errors");
  // This is a sanity bound at an imperfect wall-clock sampling point.
  if (landing) assert.ok(landing.distance < 60, label + ": logo flight missed header by " + landing.distance);
  pass += 1;
  console.log("PASS", label, "landing", landing, "cleanup", result);
  await context.close();
}

async function checkReducedMotion() {
  const context = await browser.newContext({ reducedMotion: "reduce" });
  const page = await context.newPage();
  await page.goto(BASE + "/?intro", { waitUntil: "domcontentloaded" });
  const result = await page.evaluate(async () => {
    await window.donatoIntro;
    return document.documentElement.classList.contains("intro");
  });
  assert.equal(result, false, "reduced motion should skip intro");
  await checkNoOverflow(page, "reduced-motion");
  console.log("PASS reduced motion");
  pass += 1;
  await context.close();
}

async function checkBlockedStorage() {
  const context = await browser.newContext();
  await context.addInitScript(() => {
    Object.defineProperty(window, "localStorage", {
      configurable: true,
      get() { throw new Error("Storage denied by browser"); },
    });
  });
  const page = await context.newPage();
  await page.goto(BASE + "/?intro", { waitUntil: "domcontentloaded" });
  const result = await page.evaluate(async () => {
    await window.donatoIntro;
    return document.documentElement.classList.contains("intro");
  });
  assert.equal(result, false, "blocked storage should skip intro");
  console.log("PASS blocked storage");
  pass += 1;
  await context.close();
}

async function checkOnceDailyAndBlog() {
  const context = await browser.newContext();
  const page = await context.newPage();
  await page.goto(BASE + "/?intro", { waitUntil: "domcontentloaded" });
  await page.evaluate(() => window.donatoIntro);
  await page.goto(BASE + "/blog/en/", { waitUntil: "domcontentloaded" });
  const skipped = await page.evaluate(() => ({
    intro: document.documentElement.classList.contains("intro"),
    promise: Boolean(window.donatoIntro),
  }));
  assert.ok(!skipped.intro && skipped.promise, "intro replayed on the same day");
  await page.goto(BASE + "/blog/it/?intro", { waitUntil: "domcontentloaded" });
  await page.waitForFunction(
    () => document.documentElement.classList.contains("intro-running"),
    null,
    { timeout: 3500 }
  );
  const header = await page.locator(".home-header .signal-square").count();
  assert.equal(header, 1, "blog header target missing");
  await page.evaluate(() => window.donatoIntro);
  console.log("PASS daily gate, forced blog intro");
  pass += 1;
  await context.close();
}

async function checkScrolledPage() {
  const context = await browser.newContext({ viewport: { width: 800, height: 750 } });
  const page = await context.newPage();

  // Inject a synchronous scroll at the very end of HTML parsing, before the
  // deferred Astro module starts the intro. This reproduces a restored page
  // with the shared header genuinely above the viewport.
  await page.route("**/*", async (route) => {
    const url = new URL(route.request().url());
    if (route.request().isNavigationRequest() &&
        url.pathname === "/" && url.searchParams.has("intro")) {
      const response = await route.fetch();
      const html = (await response.text()).replace(
        "</body>",
        '<script>window.scrollTo(0, 1250);document.documentElement.dataset.scrolledTest="yes";</script></body>'
      );
      await route.fulfill({ response, body: html });
    } else {
      await route.continue();
    }
  });

  await page.goto(BASE + "/?intro", { waitUntil: "domcontentloaded" });
  await page.waitForFunction(
    () => document.documentElement.classList.contains("intro-running"),
    null,
    { timeout: 3500 }
  );

  const startedOffscreen = await page.evaluate(() => {
    const r = document.querySelector(".home-header .signal-square").getBoundingClientRect();
    return document.documentElement.dataset.scrolledTest === "yes" && r.bottom < 0;
  });
  assert.ok(startedOffscreen, "test never actually started with the logo offscreen");

  await page.waitForTimeout(2300);
  const offscreenFlight = await page.evaluate(() => {
    const copy = document.querySelector("[data-dm-intro-logo]");
    if (!copy) return null;
    const r = copy.getBoundingClientRect();
    return { x: r.x + r.width / 2, y: r.y + r.height / 2 };
  });
  assert.ok(
    offscreenFlight && Math.abs(offscreenFlight.x - 400) < 8 &&
    Math.abs(offscreenFlight.y - 375) < 8,
    "offscreen header should leave the copy centered: " + JSON.stringify(offscreenFlight)
  );

  await page.evaluate(() => window.donatoIntro);
  await checkNoOverflow(page, "scrolled");
  assert.equal(
    await page.locator("[data-dm-intro-overlay]").count(),
    0, "scrolled intro overlay not cleaned up"
  );
  console.log("PASS genuinely scrolled page and center-only fallback");
  pass += 1;
  await context.close();
}

try {
  for (const width of [375, 800, 1440]) {
    for (const dark of [false, true]) await snapshotTimeline(width, dark);
  }
  await checkReducedMotion();
  await checkBlockedStorage();
  await checkOnceDailyAndBlog();
  await checkScrolledPage();
  console.log("INTRO QA PASSED:", pass, "checks. Screenshots:", OUT);
} finally {
  await browser.close();
}
