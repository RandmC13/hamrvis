// Reads the pixels the viewer drew for every run of a bundle, in each view and colour mode it offers, so the test can
// decode them back to values: the film (grains, atom columns, one layer), the x-z section and the stray field map, at
// the last frame and a middle one; and the pixels of a selected grain's outline.
// node viewer_pixels.js <url> <run id> ...   prints JSON {run: {luts, captures: [...], outline}}
const { chromium } = require("playwright");
const STRIDE = 2;
(async () => {
  const [url, ...runs] = process.argv.slice(2), out = {}, errors = [];
  const browser = await chromium.launch({ executablePath: process.env.HAMRVIS_CHROMIUM || undefined, args: ["--no-sandbox"] });
  for (const run of runs) {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, deviceScaleFactor: 1 });
    page.on("pageerror", e => errors.push(`${run}: ${e.message}`));
    await page.goto(`${url}#${run}`);
    await page.waitForFunction(r => window.HV_TEST && window.HV_TEST.D[r] && window.HV_TEST.S.stage === r, run, { timeout: 30000 });
    await page.mouse.move(0, 999);
    await page.waitForTimeout(300);
    out[run] = await page.evaluate(async ([run, STRIDE]) => {
      const T = window.HV_TEST, S = T.S, d = T.D[run], e = d.e;
      const sleep = ms => new Promise(r => setTimeout(r, ms));
      const lane = id => T.lanes().find(l => l.id === id);
      const rgb = n => { const c = getComputedStyle(document.documentElement).getPropertyValue(n).trim(); return c.startsWith("rgb") ? c.match(/\d+/g).slice(0, 3).map(Number) : [1, 3, 5].map(i => parseInt(c.slice(i, i + 2), 16)); };
      const res = { luts: { m: Array.from(T.LUT.m), agree: Array.from(T.LUT.agree), b: Array.from(T.LUT.b), s: Array.from(T.LUT.s), void: Array.from(T.LUT.void), mid: Array.from(T.LUT.mid) }, ink: rgb("--ink"), captures: [] };
      S.lag.mode = "inferred"; if (d.an0) T.applyLag();
      S.cells = S.pol = false; S.hov = S.sel = -1; S.follow = false;
      if (e.dipole) S.rw = 1;   // the reader band shades the film: keep it to 1 nm, which the test leaves out
      // draw, waiting for the snapshot or block the views need (they redraw when it arrives)
      const loading = () => d.snaps && [...d.top.values(), ...d.vol.values()].some(a => a === null);
      const settle = async () => { T.drawAll(); await sleep(30); for (let i = 0; i < 200 && loading(); i++) await sleep(25); T.drawAll(); };
      const grab = (id, st) => {
        const L = lane(id), c = L.canvas, x = c.getContext("2d"), W = c.width, H = c.height, im = x.getImageData(0, 0, W, H).data, px = [];
        for (let r = 0; r < H; r += STRIDE) for (let q = 0; q < W; q += STRIDE) { const o = (r * W + q) * 4; px.push(im[o], im[o + 1], im[o + 2]); }
        return { lane: id, ...st, W, H, stride: STRIDE, win: { ...S.win }, y0: S.y0, ysec: S.ysec, z: S.z, rw: S.rw, px };
      };
      const lastF = d.nF - 1, midF = Math.floor(d.nF / 2);
      const states = [];
      if (d.nG && (e.grain_map || e.grains.outline)) {
        for (const colour of ["mz", "agree", "peak", "minm"]) {
          if (colour === "agree" && !(d.an || (e.ref && T.D[e.ref] && T.D[e.ref].an))) continue;
          if (colour === "peak" && !(d.model && d.an)) continue;
          if (colour === "minm" && !d.g) continue;
          states.push({ view: "grains", colour, f: lastF });
        }
        if (midF !== lastF && d.g) states.push({ view: "grains", colour: "mz", f: midF });
      }
      if (d.snaps) {
        states.push({ view: "atoms", colour: "mz", f: lastF });
        if (midF !== lastF) states.push({ view: "atoms", colour: "mz", f: midF });
        if (d.model) states.push({ view: "atoms", colour: "peak", f: lastF });
        states.push({ view: "atoms", colour: "minm", f: lastF });
        if (e.snapshots.volume) states.push({ view: "layer", colour: "mz", f: lastF, layer: Math.floor(d.grid.nz / 2) });
      }
      for (const st of states) {
        S.view = st.view; S.colour = st.colour; S.layer = st.layer || 0; S.f = S.fa = st.f;
        await settle();
        res.captures.push(grab("map", st));
        if (st.view === "layer" && lane("xz")) res.captures.push(grab("xz", { ...st, view: "xz" }));
      }
      if (lane("bmap")) { S.view = "grains"; S.colour = "mz"; S.f = S.fa = lastF; await settle(); res.captures.push(grab("bmap", { view: "bmap", f: lastF })); }
      if (e.grain_map && d.gmap) {   // outline of the grain with most pixels: the pixels that change when it is selected
        const n = new Map(); for (const g of d.gmap) if (g >= 0) n.set(g, (n.get(g) || 0) + 1);
        const g = [...n.entries()].sort((a, b) => b[1] - a[1])[0][0];
        S.view = "grains"; S.colour = "mz"; S.f = S.fa = lastF;
        const c = lane("map").canvas, pix = () => c.getContext("2d").getImageData(0, 0, c.width, c.height).data;
        S.sel = -1; await settle(); const before = pix();
        S.sel = g; await settle(); const after = pix(), W = c.width, H = c.height, pts = [];
        for (let i = 0; i < W * H; i++) if (before[4 * i] !== after[4 * i] || before[4 * i + 1] !== after[4 * i + 1] || before[4 * i + 2] !== after[4 * i + 2]) pts.push([i % W, Math.floor(i / W)]);
        res.outline = { g, W, H, win: { ...S.win }, y0: S.y0, ysec: S.ysec, rw: S.rw, dipole: !!e.dipole, pts };
        S.sel = -1;
      }
      return res;
    }, [run, STRIDE]);
    await page.close();
  }
  await browser.close();
  if (errors.length) { console.error(errors.join("\n")); process.exit(1); }
  console.log(JSON.stringify(out));
})().catch(e => { console.error(e); process.exit(1); });
