// Reads back what the viewer prints and plots for every run of a bundle, at a few frames: the clock, film tooltips
// (grain and atom views, with bit cells and polarity), bits lane and stray field tooltips, and the chart series it
// drew (scrubber, life story, space-time, score against lag, read-back waveform).
// node viewer_readouts.js <url> <run id> ...   prints JSON {run: {...}}
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
    out[run] = await page.evaluate(async run => {
      const T = window.HV_TEST, S = T.S, d = T.D[run], e = d.e, res = { frames: [], story: [] };
      const sleep = ms => new Promise(r => setTimeout(r, ms));
      const until = async cond => { for (let i = 0; i < 400 && !cond(); i++) await sleep(25); };
      const strip = h => String(h).replace(/<sup class="src"[^>]*>[^<]*<\/sup>/g, "").replace(/<[^>]+>/g, "").trim();
      const rows = html => { const el = document.createElement("div"); el.innerHTML = html || "";
        return [...el.querySelectorAll("dt")].map(dt => [strip(dt.innerHTML), strip(dt.nextElementSibling.innerHTML)]); };
      const draw = async () => { T.drawAll(); await until(() => T.D[run] && (!d.ov || d.g.blocks.get(Math.floor(S.f / d.g.frames)))); T.drawAll(); };
      const L = e.size_nm, pts = [];
      for (let i = 0; i < 9; i++) for (let j = 0; j < 5; j++) pts.push([L[0] * (i + .37) / 9, L[1] * (j + .41) / 5]);
      S.lag.mode = "inferred"; if (d.an0) T.applyLag();
      const lane = id => T.lanes().find(l => l.id === id);
      const nF = d.nF, frames = [...new Set([0, Math.floor(nF / 3), Math.floor(2 * nF / 3), nF - 1])];
      for (const f of frames) {
        S.f = S.fa = f; S.hov = S.sel = -1; await draw();
        const fr = { f, clock: [...document.querySelectorAll("#clock span")].map(s => { const b = s.querySelector("b"), k = s.cloneNode(true); k.querySelector("b").remove(); return [strip(k.innerHTML), strip(b.innerHTML)]; }), tips: {} };
        for (const view of ["grains", "atoms"]) {
          if (view === "atoms" && !d.snaps) continue;
          if (view === "grains" && !d.nG) continue;
          S.view = view; S.cells = S.pol = false;
          if (view === "atoms") await until(() => rows(T.filmTip(...pts[0])).some(r => r[0].startsWith("column")));
          fr.tips[view] = pts.map(([x, y]) => [x, y, rows(T.filmTip(x, y))]);
        }
        if (e.hamr && d.g && d.an) {   // bit cells and polarity on: the tooltip adds the bit's mean m_z
          S.view = "grains"; S.cells = S.pol = true;
          fr.tips.cells = pts.map(([x, y]) => [x, y, rows(T.filmTip(x, y))]);
          S.cells = S.pol = false;
        }
        const bits = lane("bits");
        if (bits && e.hamr.n_tracks <= 8) {   // track 0 is the bottom row
          const head = e.hamr.n_tracks > 1 ? 13 : 0, py = (e.hamr.n_tracks - 1) * (40 + head) + head + 5;
          fr.bits = Array.from({ length: e.hamr.bits_per_track + 2 }, (_, b) => { const x = (b - .5) * e.hamr.BL; return [x, rows(bits.hover(x, py))]; });
        }
        const bmap = lane("bmap");
        if (bmap) { await draw(); fr.bmap = pts.map(([x, y]) => [x, y, rows(bmap.hover(x, y))]); fr.wave = T.REC.wave; }
        S.hov = -1; S.view = "grains";
        res.frames.push(fr);
      }
      S.f = S.fa = nF - 1; await draw();
      res.scrub = T.REC.scrub || null;
      res.st = document.getElementById("st-panel").hidden ? null : T.REC.st;
      res.lagscan = e.analysis ? T.REC.lagscan : null;
      res.lag = d.an0 ? T.lagOf(e).lag : null;
      if (d.g && nF > 1) {   // life stories of a few grains: around their hottest moment, and the whole run
        const an = d.an, pick = an ? an.map((a, g) => [a, g]).filter(([a]) => a.decision_step != null).slice(0, 3).map(([, g]) => g) : [0, 1].filter(g => g < d.nG);
        for (const g of pick) for (const all of [false, true]) {
          if (all && !d.model) continue;
          S.storyAll = all; S.sel = g; S.hov = -1; T.drawInspector();
          await until(() => T.REC.story && T.REC.story.g === g && T.REC.story.ready);
          if (d.ov) T.drawInspector();
          res.story.push({ all, ...T.REC.story });
        }
        S.storyAll = false; S.sel = -1;
      }
      return res;
    }, run);
    await page.close();
  }
  await browser.close();
  if (errors.length) { console.error(errors.join("\n")); process.exit(1); }
  console.log(JSON.stringify(out));
})().catch(e => { console.error(e); process.exit(1); });
