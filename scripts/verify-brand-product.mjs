/**
 * QA for the approved SVG identity and the single SwooshMania product.
 * Run against production Astro preview at http://127.0.0.1:4322.
 */
import assert from "node:assert/strict";
import fs from "node:fs/promises";
import { chromium } from "playwright";

const url = "http://127.0.0.1:4322";
const output = "artifacts/brand-product";
await fs.mkdir(output, { recursive: true });
// GitHub's Ubuntu runner already ships Chrome; avoid reinstalling OS dependencies.
const browser = await chromium.launch({ headless: true, channel: "chrome" });

async function checkOverflow(page, name) {
  const geometry = await page.evaluate(() => ({
    contentWidth: document.documentElement.scrollWidth,
    viewportWidth: document.documentElement.clientWidth,
  }));
  assert.ok(geometry.contentWidth <= geometry.viewportWidth + 2,
    name + ": unexpected horizontal scrolling " + JSON.stringify(geometry));
}

try {
  for (const width of [375, 800, 1440]) {
    const context = await browser.newContext({
      viewport: { width, height: width < 600 ? 780 : 900 },
      deviceScaleFactor: 1,
      reducedMotion: "no-preference"
    });
    const page = await context.newPage();
    const pageErrors = [];
    page.on("pageerror", error => pageErrors.push(error.message));

    await page.goto(url + "/", { waitUntil: "networkidle" });
    assert.equal(await page.locator(".home-brand .site-brand-mark").count(), 1,
      width + ": brand mark missing from homepage header");
    assert.equal(await page.locator(".featured-product").count(), 1,
      width + ": exactly one featured product expected");
    assert.equal(await page.locator("#swooshmania-title").textContent().then(x => x?.trim()), "SwooshMania.");
    for (const label of ["Unique demo players", "Steam wishlists", "Planned release"]) {
      assert.equal(await page.getByText(label, { exact: true }).count(), 1,
        width + ": missing product metric " + label);
    }
    assert.ok((await page.locator(".featured-product-stats").innerText()).includes("2027"));
    assert.ok((await page.locator(".featured-product-note").innerText()).includes("September 2026"));
    assert.equal(await page.locator(".featured-product-primary").getAttribute("href"),
      "https://store.steampowered.com/app/3560590/SwooshMania/");
    assert.ok((await page.locator("head link[rel=icon]").getAttribute("href")).includes("favicon.svg"));
    const dimensions = await page.locator(".site-brand-mark").evaluate(el => ({
      width: el.getBoundingClientRect().width,
      height: el.getBoundingClientRect().height,
    }));
    assert.ok(dimensions.width >= 24 && dimensions.height >= 24, "logo must be legible in header");
    await checkOverflow(page, width + " homepage");

    await page.screenshot({ path: output + "/home-" + width + ".png", fullPage: true });

    await page.locator(".home-brand").hover();
    await page.waitForTimeout(280);
    const hovered = await page.locator(".site-brand-mark").evaluate(el =>
      getComputedStyle(el).transform);
    assert.notEqual(hovered, "none", width + ": hover motion not active");
    await page.screenshot({ path: output + "/hover-" + width + ".png" });

    assert.deepEqual(pageErrors, [], width + ": runtime errors");
    console.log("PASS homepage + brand hover at", width, "px");
    await context.close();
  }

  const context = await browser.newContext({ viewport: { width: 375, height: 780 } });
  const page = await context.newPage();
  for (const route of ["/blog/en/", "/blog/it/", "/manifesto/"]) {
    await page.goto(url + route, { waitUntil: "domcontentloaded" });
    assert.equal(await page.locator(".site-brand-mark").count(), 1,
      route + ": shared brand mark missing");
    await checkOverflow(page, route);
    await page.screenshot({
      path: output + "/" + route.replaceAll("/", "_") + ".png",
      fullPage: false,
    });
    console.log("PASS", route);
  }
  await context.close();

  const reducedContext = await browser.newContext({ reducedMotion: "reduce" });
  const reducedPage = await reducedContext.newPage();
  await reducedPage.goto(url + "/", { waitUntil: "domcontentloaded" });
  const motion = await reducedPage.locator(".site-brand-mark").evaluate(el => ({
    transition: getComputedStyle(el).transitionDuration,
    viewTransitionName: getComputedStyle(el).viewTransitionName,
  }));
  assert.equal(motion.transition, "0s", "reduced-motion should suppress hover transition");
  assert.equal(motion.viewTransitionName, "none", "reduced-motion should suppress shared transition");
  console.log("PASS reduced-motion", motion);
  await reducedContext.close();
  console.log("IDENTITY & PRODUCT QA PASSED");
} finally {
  await browser.close();
}
