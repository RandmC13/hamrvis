// Opens every tab of a served bundle and fails on a page error, a lane that is not drawn or a lane reporting an error;
// prints each tab's "raw files (verify)" row.
// node viewer_smoke.js <url> <run id> ...   (HAMRVIS_CHROMIUM: browser to use instead of Playwright's own)
const { chromium } = require("playwright");
(async () => {
  const [url, ...runs] = process.argv.slice(2), problems = [];
  const browser = await chromium.launch({ executablePath: process.env.HAMRVIS_CHROMIUM || undefined, args: ["--no-sandbox"] });
  for (const run of runs) {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    page.on("pageerror", e => problems.push(`${run}: ${e.message}`));
    await page.goto(`${url}#${run}`);
    await page.waitForFunction(() => document.querySelectorAll(".lane canvas").length > 0, null, { timeout: 30000 });
    await page.waitForTimeout(500);
    const lanes = await page.evaluate(() => [...document.querySelectorAll(".lane")].filter(l => l.offsetParent).map(l => {
      const c = l.querySelector("canvas"), px = c.width && c.height ? c.getContext("2d").getImageData(0, 0, c.width, c.height).data : [];
      let drawn = false;
      for (let i = 3; i < px.length && !drawn; i += 4) drawn = px[i] !== 0;
      return { id: l.id, drawn };
    }));
    if (!lanes.length) problems.push(`${run}: no lanes`);
    for (const l of lanes) if (!l.drawn) problems.push(`${run}: ${l.id} is blank`);
    // smoothing is off when the viewer opens; on a 3 x 3 grid with an empty centre-right pixel, 2 x 2 boxes average
    // the non-empty pixels and keep the empty one empty
    const sm = await page.evaluate(() => { const T = window.HV_TEST, g = [1, 2, 3, 4, NaN, 6, 7, 8, 9];
      return { off: T.S.smooth === 1, out: Array.from(T.smoothed("test" + Math.random(), i => g[i], 3, 3, 2)) }; });
    const want = [1, 1.5, 2.5, 2.5, NaN, 11 / 3, 5.5, 19 / 3, 23 / 3];   // box of pixel (x, y): x - 1 .. x, y - 1 .. y
    if (!sm.off) problems.push(`${run}: smoothing is not off at start`);
    if (sm.out.some((v, i) => (v !== v) !== (want[i] !== want[i]) || (v === v && Math.abs(v - want[i]) > 1e-6))) problems.push(`${run}: smoothing gave ${sm.out}`);
    const error = await page.evaluate(() => { const e = document.getElementById("data-error"); return e && !e.hidden ? e.textContent : ""; });
    if (error) problems.push(`${run}: ${error}`);
    const verify = await page.evaluate(() => { const t = [...document.querySelectorAll("#details dt")].find(t => t.textContent === "raw files (verify)");
      return t ? t.nextElementSibling.textContent : "(no row)"; });
    console.log(`${run}: ${lanes.length} lanes; checked against raw files: ${verify}`);
    await page.close();
  }
  await browser.close();
  if (problems.length) { console.error(problems.join("\n")); process.exit(1); }
})().catch(e => { console.error(e); process.exit(1); });
