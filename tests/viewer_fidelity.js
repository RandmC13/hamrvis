// Reads back what the viewer computes and prints for each write run, at the inferred lag and with no lag: bit cells,
// stage facts, the grain under a grid of film points and the inspector rows of a few grains at the last frame.
// node viewer_fidelity.js <url> <run id> ...   prints JSON {run: {inferred: {...}, none: {...}}}
const { chromium } = require("playwright");
(async () => {
  const [url, ...runs] = process.argv.slice(2), out = {}, errors = [];
  const browser = await chromium.launch({ executablePath: process.env.HAMRVIS_CHROMIUM || undefined, args: ["--no-sandbox"] });
  for (const run of runs) {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    page.on("pageerror", e => errors.push(`${run}: ${e.message}`));
    await page.goto(`${url}#${run}`);
    await page.waitForFunction(r => window.HV_TEST && window.HV_TEST.D[r] && window.HV_TEST.S.stage === r, run, { timeout: 30000 });
    await page.waitForTimeout(300);
    out[run] = await page.evaluate(run => {
      const T = window.HV_TEST, d = T.D[run], e = d.e, strip = h => String(h).replace(/<sup class="src"[^>]*>[^<]*<\/sup>/g, "").replace(/<[^>]+>/g, "").trim(), res = {};
      const L = e.size_nm, pts = [];
      for (let i = 1; i < 40; i++) for (let j = 1; j < 20; j++) pts.push([L[0] * i / 40, L[1] * j / 20]);
      const heated = d.an.map((a, g) => [a, g]).filter(([a]) => a.decision_step != null).slice(0, 6).map(([, g]) => g);
      for (const mode of ["inferred", "none"]) {
        T.S.lag.mode = mode; T.applyLag();
        T.S.f = d.nF - 1;
        const insp = {};
        for (const g of heated) {
          T.S.sel = g; T.S.hov = -1; T.drawInspector();
          insp[g] = [...document.querySelectorAll("#inspector dl.facts dt")].map(dt => [strip(dt.innerHTML), strip(dt.nextElementSibling.innerHTML)]);
        }
        res[mode] = {
          lag: T.lagOf(e).lag, frame: T.S.f,
          cells: T.cellParts(d).map(c => [c.t, c.b, c.want, c.parts]),
          facts: T.stageFacts(run).filter(Boolean).map(([k, v]) => [strip(k), strip(v)]),
          hover: pts.map(([x, y]) => [x, y, T.grainAt(x, y)]), inspector: insp,
        };
      }
      return res;
    }, run);
    await page.close();
  }
  await browser.close();
  if (errors.length) { console.error(errors.join("\n")); process.exit(1); }
  console.log(JSON.stringify(out));
})().catch(e => { console.error(e); process.exit(1); });
