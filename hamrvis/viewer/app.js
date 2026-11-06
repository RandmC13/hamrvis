// hamrvis viewer: loads a bundle made by `hamrvis reduce` (docs/bundle-format.md) and draws it.
// Data comes from bundle/<path> over HTTP (`hamrvis view`) or from window.HV_EMBED (`hamrvis export`).
"use strict";
const HAMRVIS_VERSION = "0.9.0"; // same as hamrvis/__init__.py (a test checks)

// ------------------------------------------------------------------------------------------------ bundle access
const EMB = window.HV_EMBED || null;

async function gunzip(buf) {
  const s = new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Response(s).arrayBuffer();
}
async function getBuf(path) {
  if (EMB) {
    if (!(path in EMB)) throw new Error(`${path} is not in this exported page`);
    const bin = atob(EMB[path]), u = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
    return gunzip(u);                                          // export gzips every file
  }
  const r = await fetch("bundle/" + path);
  if (!r.ok) throw new Error(`bundle/${path}: HTTP ${r.status}`);
  const b = await r.arrayBuffer();
  return path.endsWith(".gz") ? gunzip(b) : b;
}
async function getRange(path, a, n) {                             // n bytes from byte a (long grain series)
  const r = await fetch("bundle/" + path, { headers: { Range: `bytes=${a}-${a + n - 1}` } });
  if (r.status !== 206) throw new Error(`bundle/${path}: HTTP ${r.status} to a byte range request (use hamrvis view)`);
  return r.arrayBuffer();
}
const has = path => EMB ? path in EMB : true;
async function tryBuf(path) { try { return await getBuf(path); } catch (e) { return null; } }

// ------------------------------------------------------------------------------------------------ state
const S = { M: null, stage: null, f: 0, fa: 0, view: "grains", colour: "mz", layer: 0, sel: -1, hov: -1, ysec: 0,
            z: 0, rw: 10, playing: false, storyAll: false, stAll: false, cells: false, pol: false, hx: null,
            win: { x0: 0, x1: 1 }, dspan: 1, y0: 0, follow: true };
const D = {};                                                    // per-stage data, filled by loadStage
const $ = id => document.getElementById(id);
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt = (v, d = 1) => (v == null || !isFinite(v)) ? "–" : v.toFixed(d);
const sgn = v => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v);
const kv = rows => `<dl class="facts">${rows.filter(Boolean).map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("")}</dl>`;

// ------------------------------------------------------------------------------------------------ HAMR model
// Port of hamrvis/vampire.py Hamr (itself src/hamr/hamr_continuous.cpp + field_profiles.cpp). Steps count from
// the start of the write run; lengths nm.
function Model(p) {
  const w = p.written ?? p.bits_per_track;                         // bits VAMPIRE's loop writes per track
  const trackY = t => p.TW * (Math.floor(t * w / Math.max(p.bits_per_track, 1)) + 0.5) + p.pad;
  const head = step => {
    const k = step - p.eq, on = k >= 1 && k <= p.laser_steps, kk = Math.max(k, 1);
    const track = Math.floor((kk - 1) / p.nt), j = (kk - 1) % p.nt + 1;
    return { x: -p.BL / 2 - p.NPS + (j + 1) * p.step_nm, y: trackY(track), on, track, j };
  };
  const field = step => {                                         // signed B_applied (T)
    const h = head(step); if (!h.on) return 0;
    const j0 = h.j - 1, b = Math.floor((j0 - p.ne) / p.nb), i = (j0 - p.ne) - b * p.nb;
    if (j0 < p.ne || b >= p.bits_per_track) return 0;
    let mag = Math.min((i + 1) / p.rise, 1);
    if (i >= p.nb - p.fall) mag = (p.nb - 1 - i) / p.fall;
    const pol = p.sequence[Math.min(h.track * p.bits_per_track + b, p.sequence.length - 1)] || 0;
    return mag * p.Hpeak * pol;
  };
  const T = (x, y, step, h = head(step)) => h.on ? p.Tmin + (p.Tmax - p.Tmin) *
    Math.exp(-((x - h.x) ** 2) / (2 * p.sigma_x ** 2) - ((y - h.y) ** 2) / (2 * p.sigma_y ** 2)) : p.Tmin;
  const inBox = (x, y, h) => h.on && Math.abs(x - (h.x - p.NPS)) <= p.box_x / 2 && Math.abs(y - h.y) <= p.box_y / 2;
  const bit = step => {                                           // [track, bit] being written, or null
    const h = head(step), b = Math.floor((h.j - 1 - p.ne) / p.nb);
    return h.on && h.j - 1 >= p.ne && b < w ? [h.track, b] : null;
  };
  const peakT = (x, y) => {                                       // hottest laser temperature ever at x, y
    const x0 = -p.BL / 2 - p.NPS + 2 * p.step_nm, x1 = -p.BL / 2 - p.NPS + (p.nt + 1) * p.step_nm;
    const dx = x < x0 ? x - x0 : x > x1 ? x - x1 : 0;
    let best = 0;
    for (let t = 0; t < p.n_tracks; t++) {
      const yt = trackY(t);
      best = Math.max(best, Math.exp(-(dx * dx) / (2 * p.sigma_x ** 2) - ((y - yt) ** 2) / (2 * p.sigma_y ** 2)));
    }
    return p.Tmin + (p.Tmax - p.Tmin) * best;
  };
  const passStep = (x, track) => p.eq + track * p.nt + Math.round((x + p.BL / 2 + p.NPS) / p.step_nm) - 1;
  return { p, head, field, T, inBox, bit, peakT, passStep };
}

// ------------------------------------------------------------------------------------------------ colour maps
const hex = c => c.startsWith("rgb") ? c.match(/\d+/g).map(Number) : [1, 3, 5].map(i => parseInt(c.slice(i, i + 2), 16));
// diverging lo - mid - hi, saturating at |v| = 0.85 (the old explorer's map); 256 entries for v in -1..1
function divLut(lo, mid, hi) {
  const a = hex(css(lo)), m = hex(css(mid)), b = hex(css(hi)), L = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const v = Math.max(-1, Math.min(1, (i / 127.5 - 1) / .85)), e = v > 0 ? b : a, t = Math.abs(v);
    for (let c = 0; c < 3; c++) L[3 * i + c] = m[c] + (e[c] - m[c]) * t;
  }
  return L;
}
function seqLut(names) {
  const st = names.map(n => hex(css(n))), L = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const t = i / 255 * (st.length - 1), k = Math.min(Math.floor(t), st.length - 2), f = t - k;
    for (let c = 0; c < 3; c++) L[3 * i + c] = st[k][c] + (st[k + 1][c] - st[k][c]) * f;
  }
  return L;
}
let LUT = {};
function buildLuts() {
  LUT.m = divLut("--down", "--mid", "--up");
  LUT.agree = divLut("--bad", "--mid", "--ok");
  LUT.b = divLut("--b-dn", "--b-mid", "--b-up");
  LUT.s = seqLut(["--s0", "--s1", "--s2", "--s3", "--s4"]);
  LUT.void = hex(css("--void"));
  LUT.mid = hex(css("--mid"));
}
const idxDiv = v => Math.max(0, Math.min(255, Math.round((v + 1) * 127.5)));     // -1..1
const idxSeq = v => Math.max(0, Math.min(255, Math.round(v * 255)));             // 0..1
const rgbOf = (L, i) => [L[3 * i], L[3 * i + 1], L[3 * i + 2]];
const rgb = (L, i) => `rgb(${L[3 * i]},${L[3 * i + 1]},${L[3 * i + 2]})`;

// ------------------------------------------------------------------------------------------------ loading
async function loadStage(st) {
  if (D[st]) return D[st];
  const e = S.M.byId[st], d = { e, st };
  d.nG = e.grains.id.length;
  if (e.grain_series) {
    const nF = e.grain_series.frames, nL = e.grain_series.lead.length;
    const gt = e.grain_series.dtype || "i16";                  // grains.i16 (x 1/32767) or grains.f32 (exact)
    d.gs = gt === "f32" ? 1 : 32767;
    if (e.grain_series.overview) {                             // too long to hold: blocks of frames on demand, overview for whole-run views
      d.g = { blocks: new Map(), T: gt === "f32" ? Float32Array : Int16Array, path: `${st}/grains.${gt}` };
      d.g.frames = Math.max(1, Math.ceil(16e6 / (d.nG * 4 * d.g.T.BYTES_PER_ELEMENT)));
      d.ov = new Int16Array(await getBuf(`${st}/grains_overview.i16`));
      d.ovB = e.grain_series.overview.bins;
    } else d.g = gt === "f32" ? new Float32Array(await getBuf(`${st}/grains.f32`)) : new Int16Array(await getBuf(`${st}/grains.i16`));
    d.lead = new Float64Array(await getBuf(`${st}/grain_lead.f64`));
    Object.assign(d, { nF, nL });
    d.time = Float64Array.from({ length: nF }, (_, f) => d.lead[f * nL]);
    if (d.ov) await loadFrame(d, nF - 1);                      // the final state is always wanted
    d.steps = Float64Array.from(d.time, t => Math.round(t / e.dt));
  } else if (e.snapshots && e.snapshots.frames.length > 1) {   // no grain series (e.g. a continuous film): one frame per snapshot
    const fr = e.snapshots.frames;
    Object.assign(d, { nF: fr.length, time: Float64Array.from(fr, s => s.time) });
    d.steps = Float64Array.from(d.time, t => Math.round(t / e.dt));
  } else Object.assign(d, { nF: 1, time: [0], steps: [0] });
  const f64 = e.output && e.output.dtype === "f64", out = await tryBuf(`${st}/output.${f64 ? "f64" : "f32"}`);
  if (out) { d.out = f64 ? new Float64Array(out) : new Float32Array(out); d.cols = e.output.columns; }
  if (e.grain_map) d.gmap = new Int16Array(await getBuf(`${st}/grain_map.i16`));
  if (e.snapshots) {
    d.snaps = e.snapshots.frames.map(s => ({ ...s, step: Math.round(s.time / e.dt) }));
    d.grid = e.snapshots.grid;
    d.topCount = new Uint8Array(await getBuf(`${st}/top_count.u8`));
    d.minm = new Uint8Array(await getBuf(`${st}/min_m.u8`));
    d.top = new Map(); d.vol = new Map();
  }
  if (e.dipole) d.B = new Float32Array(await getBuf(`${st}/dipole.f32`));
  if (e.hamr) d.model = Model(e.hamr);
  d.an = e.analysis ? e.analysis.grains : null;
  return (D[st] = d);
}
// snapshot grid encodings (reduce --precision): array type, scale, empty-pixel marker
const SNAP = { i8: [Int8Array, 127, -128], i16: [Int16Array, 32767, -32768], f32: [Float32Array, 1, NaN] };
const snapPath = (d, kind, idx) => `${d.st}/${kind}/${String(idx).padStart(8, "0")}.${d.e.snapshots.dtype || "i8"}.gz`;
const snapVal = (a, i) => { const v = a[i]; return v === a.E || v !== v ? NaN : v / a.S; };   // NaN: no atoms
// one snapshot array, cached; returns null while loading and redraws when it arrives
function snapArray(d, kind, idx) {
  const cache = d[kind];
  if (cache.has(idx)) return cache.get(idx);
  const path = snapPath(d, kind, idx);
  if (!has(path)) { cache.set(idx, undefined); return undefined; }
  cache.set(idx, null);
  const [T, sc, empty] = SNAP[d.e.snapshots.dtype || "i8"];
  getBuf(path).then(b => { const a = new T(b); a.S = sc; a.E = empty; cache.set(idx, a); if (cache.size > 80) cache.delete(cache.keys().next().value); redraw(); })
    .catch(err => { cache.set(idx, undefined); d.err = `snapshot ${idx}: ${err.message}`; redraw(); });
  return null;
}
function nearestSnap(d, step) {
  if (!d.snaps) return null;
  let best = null;
  for (const s of d.snaps) if (has(snapPath(d, "top", s.index)) &&
    (!best || Math.abs(s.step - step) < Math.abs(best.step - step))) best = s;
  return best;
}

// ------------------------------------------------------------------------------------------------ grain helpers
// Long series: a block of frames is fetched when first wanted and values read as NaN until it arrives. drawAll waits for
// the block of the current frame, so nothing drawn for the current frame is ever a placeholder.
const MAX_BLOCKS = 32;
function frameBlock(d, f) {
  const G = d.g, b = Math.floor(f / G.frames), have = G.blocks.get(b);
  if (have !== undefined) return have;
  G.blocks.set(b, null);
  loadBlock(d, b).then(redraw).catch(err => { d.err = `grain series: ${err.message}`; G.blocks.delete(b); redraw(); });
  return null;
}
async function loadBlock(d, b) {
  const G = d.g, sz = d.nG * 4 * G.T.BYTES_PER_ELEMENT, n = Math.min(G.frames, d.nF - b * G.frames);
  const a = new G.T(await getRange(G.path, b * G.frames * sz, n * sz));
  G.blocks.set(b, a);
  const pin = Math.floor((d.nF - 1) / G.frames);
  for (const k of G.blocks.keys()) { if (G.blocks.size <= MAX_BLOCKS) break; if (k !== pin && k !== b && G.blocks.get(k)) G.blocks.delete(k); }
  return a;
}
const loadFrame = (d, f) => loadBlock(d, Math.floor(f / d.g.frames));
const frameReady = (d, f) => !d.ov || !!frameBlock(d, f);
function gv(d, f, g, c) {
  if (!d.ov) return d.g[(f * d.nG + g) * 4 + c];
  const a = frameBlock(d, f);
  return a ? a[((f % d.g.frames) * d.nG + g) * 4 + c] : NaN;
}
const gmz = (d, f, g) => gv(d, f, g, 2) / d.gs * gv(d, f, g, 3) / d.gs;
const gml = (d, f, g) => gv(d, f, g, 3) / d.gs;
function frameAt(d, step) {                                       // last frame at or before step
  let lo = 0, hi = d.nF - 1;
  while (lo < hi) { const m = (lo + hi + 1) >> 1; if (d.steps[m] <= step) lo = m; else hi = m - 1; }
  return lo;
}
const _minGrain = {};
function minGrainM(d, g) {
  if (!_minGrain[d.st]) {
    const a = new Float32Array(d.nG).fill(1);
    if (d.ov) { for (let b = 0; b < d.ovB; b++) for (let k = 0; k < d.nG; k++) a[k] = Math.min(a[k], d.ov[(b * d.nG + k) * 4 + 2] / 32767); }
    else for (let f = 0; f < d.nF; f++) for (let k = 0; k < d.nG; k++) a[k] = Math.min(a[k], gml(d, f, k));
    _minGrain[d.st] = a;
  }
  return _minGrain[d.st][g];
}
// written as intended at frame f: melted at some point in the run, and pointing along its bit's field now
const okAt = (d, f, a, g) => a.intended !== 0 && Math.sign(gmz(d, f, g)) === a.intended;
const melted = a => a.category !== "unheated" && a.intended !== 0;
const refD = d => d.e.ref ? D[d.e.ref] : null;                    // the write before this run, if any
const wantOf = (d, g) => {                                        // intended z direction of grain g (write analysis);
  const an = d.an || (refD(d) && refD(d).an);                     // 0 for grains the laser never melts
  return an && an[g] && an[g].category !== "unheated" ? an[g].intended : 0;
};
function grainAt(x, y) {                                          // grain index at film position, from the atom map
  const d = cur(), M = d.e.grain_map;
  if (!M) return d.e.grains.outline ? d.e.grains.outline.findIndex(o => o && inPoly(o, x, y)) : -1;
  const ix = Math.floor(x / M.pitch_nm), iy = Math.floor(y / M.pitch_nm);
  return ix >= 0 && ix < M.nx && iy >= 0 && iy < M.ny ? d.gmap[iy * M.nx + ix] : -1;
}
function inPoly(o, x, y) {
  let c = false;
  for (let i = 0, j = o.length - 1; i < o.length; j = i++)
    if ((o[i][1] > y) !== (o[j][1] > y) && x < (o[j][0] - o[i][0]) * (y - o[i][1]) / (o[j][1] - o[i][1]) + o[i][0]) c = !c;
  return c;
}
const _edges = new Map();                                         // grain index -> its boundary pixels
function grainEdge(d, g) {
  const key = `${d.st}|${g}`;
  if (_edges.has(key)) return _edges.get(key);
  const M = d.e.grain_map, gm = d.gmap, px = [];
  for (let iy = 0; iy < M.ny; iy++) for (let ix = 0; ix < M.nx; ix++) {
    const i = iy * M.nx + ix; if (gm[i] !== g) continue;
    const n = (dx, dy) => gm[((iy + dy + M.ny) % M.ny) * M.nx + (ix + dx + M.nx) % M.nx];
    if (n(1, 0) !== g || n(-1, 0) !== g || n(0, 1) !== g || n(0, -1) !== g) px.push(i);
  }
  _edges.set(key, px);
  return px;
}

// colour of grain g at frame f in the current colour mode (RGB triple)
function grainColour(d, f, g) {
  const e = d.e;
  if (S.colour === "peak" && d.an) return rgbOf(LUT.s, idxSeq((d.an[g].peakT - e.hamr.Tmin) / (e.hamr.Tmax - e.hamr.Tmin)));
  if (S.colour === "minm" && d.g) return rgbOf(LUT.s, idxSeq(1 - minGrainM(d, g)));
  if (!d.g) return LUT.void;
  if (S.colour === "agree") { const w = wantOf(d, g); return w ? rgbOf(LUT.agree, idxDiv(gmz(d, f, g) * w)) : LUT.mid; }
  return rgbOf(LUT.m, idxDiv(gmz(d, f, g)));
}

// ------------------------------------------------------------------------------------------------ canvas helpers
const dpr = () => window.devicePixelRatio || 1;
function fit(c, h) {
  const w = c.clientWidth || c.parentElement.clientWidth, r = dpr();
  if (c.width !== Math.round(w * r) || c.height !== Math.round(h * r)) { c.width = Math.round(w * r); c.height = Math.round(h * r); }
  c.style.height = h + "px";
  const ctx = c.getContext("2d"); ctx.setTransform(r, 0, 0, r, 0, 0);
  return { ctx, w, h };
}
const off = (w, h) => { const c = document.createElement("canvas"); c.width = w; c.height = h; return c; };
// image of nx x ny pixels from a colour function (y flipped: row 0 is the bottom of the film)
function pixImage(nx, ny, colour) {
  const c = off(nx, ny), x = c.getContext("2d"), im = x.createImageData(nx, ny);
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) {
    const col = colour(iy * nx + ix), o = ((ny - 1 - iy) * nx + ix) * 4;
    im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
  }
  x.putImageData(im, 0, 0); return c;
}

// The part of an nx x ny grid (spanning the film in x) inside the x window and the rows fy0..fy1 (fractions), at no more
// than twice the screen resolution: canvases stay within the browser's limits on any film.
function gridImage(nx, ny, colour, fy0, fy1, w, h) {
  const L = Lx(), r = dpr();
  const ix0 = Math.max(0, Math.floor(S.win.x0 / L * nx)), ix1 = Math.min(nx, Math.ceil(S.win.x1 / L * nx));
  const iy0 = Math.max(0, Math.floor(fy0 * ny)), iy1 = Math.min(ny, Math.ceil(fy1 * ny));
  const stepx = Math.max(1, Math.floor((ix1 - ix0) / (2 * w * r))), stepy = Math.max(1, Math.floor((iy1 - iy0) / (2 * h * r)));
  const cw = Math.ceil((ix1 - ix0) / stepx), ch = Math.ceil((iy1 - iy0) / stepy);
  const c = off(cw, ch), x = c.getContext("2d"), im = x.createImageData(cw, ch);
  for (let jy = 0; jy < ch; jy++) for (let jx = 0; jx < cw; jx++) {
    const col = colour((iy0 + jy * stepy) * nx + ix0 + jx * stepx), o = ((ch - 1 - jy) * cw + jx) * 4;
    im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
  }
  x.putImageData(im, 0, 0);
  return { c, x0: ix0 / nx * L, x1: Math.min(ix0 + cw * stepx, nx) / nx * L, fy0: iy0 / ny, fy1: Math.min(iy0 + ch * stepy, ny) / ny, win: S.win };
}
// draw a gridImage: x from the window, y from the lane's screen positions of its rows
function putImage(ctx, im, w, yTop, yBot, smooth = false) {
  ctx.imageSmoothingEnabled = smooth;
  ctx.drawImage(im.c, X(im.x0, w), yTop, X(im.x1, w) - X(im.x0, w), yBot - yTop);
}

// ------------------------------------------------------------------------------------------------ film window
// One x window (S.win, nm) for every lane on the film axis. Lanes on the film show it at the true aspect ratio; a film
// taller than its lane scrolls in y (S.y0: y at the bottom of the lane), one too flat for its lane is stretched in y.
const MIN_FILM_PX = 120;
const span = () => S.win.x1 - S.win.x0;
const isDefaultWindow = () => S.win.x0 === 0 && S.win.x1 === S.dspan;
function fitWindow() {                                            // the view a run opens on: all of it if readable
  const w = $("axis-top").clientWidth || 1000, L = Lx();
  S.dspan = Math.min(L, w * Ly() / MIN_FILM_PX);
  S.win = { x0: 0, x1: S.dspan };
  S.y0 = 0;
}
function setWindow(x0, x1) {
  const L = Lx(), n = Math.min(L, Math.max(x1 - x0, 2 * pixelNm()));
  x0 = Math.max(0, Math.min(L - n, x0));
  S.win = { x0, x1: x0 + n };
  redraw();
}
const pixelNm = () => { const d = cur(); return d.grid ? d.grid.pitch_nm : d.e.grain_map ? d.e.grain_map.pitch_nm : 0.25; };
function zoomAt(x, factor) { setWindow(x - (x - S.win.x0) * factor, x + (S.win.x1 - x) * factor); }
function filmH(w) {
  const h = Math.round(w * Ly() / S.dspan);
  return Ly() > Lx() ? Math.min(h, Math.round(.6 * innerHeight)) : h;
}
function filmScale(w, h) {
  const sx = w / span(), full = Ly() * sx <= h + 1;
  const sy = full ? h / Ly() : sx, y0 = full ? 0 : Math.max(0, Math.min(Ly() - h / sy, S.y0));
  return { sx, sy, y0, y1: full ? Ly() : y0 + h / sy, scrolled: !full, py: y => h - (y - y0) * sy };
}
function filmImage(V, nx, ny, colour, w, h, smooth) {
  const im = gridImage(nx, ny, colour, V.y0 / Ly(), V.y1 / Ly(), w, h);
  im.top = V.scrolled ? V.py(im.fy1 * Ly()) : 0; im.bot = V.scrolled ? V.py(im.fy0 * Ly()) : h; im.smooth = smooth;
  return im;
}
const putFilm = (ctx, im, w) => putImage(ctx, im, w, im.top, im.bot, im.smooth);
function drawScrollBar(ctx, V, w, h) {                            // where the lane is in a film taller than it
  if (!V.scrolled) return;
  const y = h * (1 - V.y1 / Ly()), hh = Math.max(12, h * (V.y1 - V.y0) / Ly());
  ctx.fillStyle = css("--ink"); ctx.globalAlpha = .55; ctx.fillRect(w - 5, y, 4, hh); ctx.globalAlpha = 1;
}
const scrollNote = V => V.scrolled ? `y ${fmt(V.y0, 0)}–${fmt(V.y1, 0)} of ${fmt(Ly(), 0)} nm` : "";
// keep the head in view while it moves
function followHead(d) {
  if (!S.follow || !d.model) return;
  const h = d.model.head(d.steps[S.f]), w = $("axis-top").clientWidth || 1000, lane = filmH(w);
  if (!h.on) return;
  const n = span();
  if (n < Lx() && (h.x < S.win.x0 + .15 * n || h.x > S.win.x0 + .85 * n)) { const x0 = Math.max(0, Math.min(Lx() - n, h.x - .35 * n)); S.win = { x0, x1: x0 + n }; }
  const V = filmScale(w, lane);
  if (V.scrolled && (h.y < V.y0 + .15 * (V.y1 - V.y0) || h.y > V.y1 - .15 * (V.y1 - V.y0))) S.y0 = h.y - .5 * (V.y1 - V.y0);
}
function drawMinimap() {
  const d = cur(), c = $("minimap"), show = span() < Lx() - 1e-9;
  $("mini-row").hidden = !show;
  if (!show) return;
  const { ctx, w, h } = fit(c, 16), L = Lx();
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = css("--lane"); ctx.fillRect(0, 0, w, h);
  if (d.model) {
    const hd = d.model.head(d.steps[S.f]);
    if (hd.on) { ctx.fillStyle = css("--ink"); ctx.fillRect(hd.x / L * w - 1, 0, 2, h); }
  }
  ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5; ctx.strokeRect(S.win.x0 / L * w + .5, .5, Math.max(3, span() / L * w) - 1, h - 1);
  $("follow").checked = S.follow;
  $("follow").parentElement.hidden = !d.model;
}

// ------------------------------------------------------------------------------------------------ lanes
let LANES = [];
const cur = () => D[S.stage];
const Lx = () => cur().e.size_nm[0], Ly = () => cur().e.size_nm[1];
const X = (x, w) => (x - S.win.x0) / span() * w;
const XI = (px, w) => px / w * span() + S.win.x0;

function filmLane(title) {
  return {
    id: "map", group: "top", title,
    height: filmH, scroll: true,
    draw: drawFilm,
    hover: filmTip,
    click(x, y) { const g = grainAt(x, y); S.sel = g === S.sel ? -1 : g; S.ysec = y; redraw(); },   // clicking the selected grain deselects it
    drag: dragSection,
  };
}

const dragSection = (x, y) => { S.ysec = Math.min(Ly(), Math.max(0, y)); };   // section line follows a held mouse button
let _lastPix = null;
function drawFilm(ctx, w, h) {
  const d = cur(), e = d.e, V = filmScale(w, h), sx = V.sx, sy = V.sy, M = e.grain_map;
  ctx.fillStyle = `rgb(${LUT.void})`; ctx.fillRect(0, 0, w, h);
  const view = d.snaps ? S.view : "grains", world = () => { ctx.save(); ctx.translate(-S.win.x0 * sx, V.y0 * sy); };
  let note = "";
  if (view === "grains") {
    if (M) {
      const pal = Array.from({ length: d.nG }, (_, g) => grainColour(d, S.f, g));
      putFilm(ctx, filmImage(V, M.nx, M.ny, i => d.gmap[i] < 0 ? LUT.void : pal[d.gmap[i]], w, h), w);
    } else if (e.grains.outline) {   // no atom coordinates in the bundle: VAMPIRE's Voronoi outlines (approximate)
      world();
      e.grains.outline.forEach((o, g) => { if (!o) return; const c = grainColour(d, S.f, g); ctx.fillStyle = `rgb(${c})`; ctx.beginPath();
        o.forEach(([x, y], i) => i ? ctx.lineTo(x * sx, h - y * sy) : ctx.moveTo(x * sx, h - y * sy)); ctx.fill(); });
      ctx.restore();
      note = "approximate Voronoi outlines (no atom coordinates)";
    }
  } else {
    let img = pixelImage(d, view, V, w, h);
    const key = `${d.st}|${view}|${S.colour}|${S.layer}`;
    if (img) _lastPix = { img, key };
    else if (_lastPix && _lastPix.key === key) img = _lastPix.img;   // next snapshot still loading: no blank frame
    if (img) putFilm(ctx, img, w);
    const s = nearestSnap(d, d.steps[S.f]);
    if (s) note = snapLabel(d, s);
  }
  // overlays below are in film coordinates: x * sx, h - y * sy
  world();
  const fw = Lx() * sx;
  // hovered and selected grain outlines, from the atom map
  if (M) for (const [g, col, lw] of [[S.hov, "--muted", 1], [S.sel, "--ink", 1]]) {
    if (g < 0) continue;
    const cw = fw / M.nx, ch = Ly() * sy / M.ny;
    ctx.fillStyle = css(col); const pw = Math.max(cw * (lw + .2), cw < 1 ? 1 : 0), ph = Math.max(ch * (lw + .2), ch < 1 ? 1 : 0);
    for (const i of grainEdge(d, g)) ctx.fillRect((i % M.nx) * cw, h - (Math.floor(i / M.nx) + 1) * ch, pw, ph);
  }
  else if (e.grains.outline) for (const [g, col, lw] of [[S.hov, "--muted", 1.2], [S.sel, "--ink", 2]]) {   // Voronoi fallback
    const o = g >= 0 && e.grains.outline[g]; if (!o) continue;
    ctx.beginPath(); o.forEach(([x, y], i) => i ? ctx.lineTo(x * sx, h - y * sy) : ctx.moveTo(x * sx, h - y * sy)); ctx.closePath();
    haloStroke(ctx, css(col), lw);
  }
  ctx.setLineDash([5, 4]); ctx.lineWidth = 1; ctx.strokeStyle = css("--ink"); ctx.globalAlpha = .5;
  const yl = h - S.ysec * sy; ctx.beginPath(); ctx.moveTo(0, yl); ctx.lineTo(fw, yl); ctx.stroke();
  ctx.setLineDash([]); ctx.globalAlpha = 1;
  if (S.cells && cellGeom(d)) note = [note, drawCells(ctx, sx, sy, h)].filter(Boolean).join("; ");
  if (d.model) drawHead(ctx, d, sx, sy, h);
  if (e.dipole) drawReaderBand(ctx, fw, h, sy);
  ctx.restore();
  drawScrollBar(ctx, V, w, h);
  if (view !== "grains") {   // atom views show the nearest snapshot: say so on the image when it is not this time
    const s = nearestSnap(d, d.steps[S.f]);
    if (s && Math.abs(s.step - d.steps[S.f]) * e.dt > 1e-13) {
      const txt = `snapshot ${timeText(d, s.time)}`;
      ctx.font = `600 12px ${css("--font")}`; const tw = ctx.measureText(txt).width;
      ctx.fillStyle = css("--ink"); ctx.globalAlpha = .8; ctx.fillRect(6, h - 26, tw + 12, 20); ctx.globalAlpha = 1;
      ctx.fillStyle = css("--lane"); ctx.textAlign = "left"; ctx.fillText(txt, 12, h - 12);
    }
  }
  setNote("map", [note, `section y ${fmt(S.ysec)} nm`, scrollNote(V)]);
}
// bit cells: each bit's region (shifted upstream by the freeze lag) on every track; the band hotter than Tc
function cellGeom(d = cur()) {
  const W = d.e.hamr ? d.e : d.e.ref && S.M.byId[d.e.ref]; if (!W || !W.hamr) return null;
  const p = W.hamr, lag = W.freeze_lag_nm || 0, Tc = W.Tc;
  const half = Tc > p.Tmin && Tc < p.Tmax ? p.sigma_y * Math.sqrt(2 * Math.log((p.Tmax - p.Tmin) / (Tc - p.Tmin))) : 0;
  return { p, lag, Tc, half };
}
// atom-weighted grain shares of every bit cell, inside the writable band (|y - track centre| <= half, or the
// whole track without a Tc): [{t, b, want, parts: [[grain, weight], ...]}], from the grain map when there is one
const _cells = {};
function cellParts(d) {
  if (_cells[d.st]) return _cells[d.st];
  const C = cellGeom(d), p = C.p, out = [], M = d.e.grain_map, e = d.e;
  for (let t = 0; t < p.n_tracks; t++) for (let b = 0; b < p.bits_per_track; b++) {
    const yc = p.pad + (t + .5) * p.TW, hw = C.half || p.TW / 2, x0 = b * p.BL - C.lag, x1 = x0 + p.BL;
    const w = new Map();
    if (M && d.gmap) {
      const P = M.pitch_nm;
      for (let iy = Math.max(0, Math.floor((yc - hw) / P)); iy < Math.min(M.ny, Math.ceil((yc + hw) / P)); iy++)
        for (let ix = Math.floor(x0 / P); ix < Math.ceil(x1 / P); ix++) {
          const g = d.gmap[iy * M.nx + ((ix % M.nx) + M.nx) % M.nx]; if (g >= 0) w.set(g, (w.get(g) || 0) + 1);
        }
    } else e.grains.x.forEach((x, g) => { const y = e.grains.y[g]; if (x >= x0 && x < x1 && Math.abs(y - yc) <= hw) w.set(g, e.grains.atoms[g]); });
    out.push({ t, b, want: Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz), parts: [...w] });
  }
  return (_cells[d.st] = out);
}
const cellMz = (d, c, f) => { let s = 0, n = 0; for (const [g, w] of c.parts) { s += w * gmz(d, f, g); n += w; } return n ? s / n : 0; };
const bitOk = (d, c, f) => c.want !== 0 && c.parts.length > 0 && Math.sign(cellMz(d, c, f)) === c.want;
function drawCells(ctx, sx, sy, h) {
  const C = cellGeom(); if (!C) return "";
  const { p, lag, Tc, half } = C;
  ctx.font = `600 10.5px ${css("--font")}`; ctx.textAlign = "center";
  for (let t = 0; t < p.n_tracks; t++) {
    const y0 = p.pad + t * p.TW, y1 = y0 + p.TW, yc = y0 + p.TW / 2;
    ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1; ctx.globalAlpha = .7; ctx.setLineDash([3, 3]);
    for (let b = 0; b <= p.bits_per_track; b++) { const x = (b * p.BL - lag) * sx; ctx.beginPath(); ctx.moveTo(x, h - y1 * sy); ctx.lineTo(x, h - y0 * sy); ctx.stroke(); }
    ctx.setLineDash([]); ctx.strokeRect(-lag * sx, h - y1 * sy, p.bits_per_track * p.BL * sx, p.TW * sy);
    if (half) { ctx.lineWidth = 1.6; ctx.globalAlpha = .9;
      for (const yy of [yc - half, yc + half]) { ctx.beginPath(); ctx.moveTo(-lag * sx, h - yy * sy); ctx.lineTo((p.bits_per_track * p.BL - lag) * sx, h - yy * sy); ctx.stroke(); } }
    ctx.globalAlpha = 1; ctx.fillStyle = css("--ink");
    for (let b = 0; b < p.bits_per_track; b++) if (p.BL * sx > 14) ctx.fillText(b + 1, ((b + .5) * p.BL - lag) * sx, h - y1 * sy + 12);
  }
  const d = cur();
  if (S.pol && d.g) {   // outline each cell's writable part in the colour of its mean m_z; mark it against the bit
    const hw = half || p.TW / 2;
    for (const c of cellParts(d)) {
      const v = cellMz(d, c, S.f), yc = p.pad + (c.t + .5) * p.TW, x0 = (c.b * p.BL - lag) * sx;
      const rx = x0 + 2, ry = h - (yc + hw) * sy + 2, rw = p.BL * sx - 4, rh = 2 * hw * sy - 4;
      ctx.strokeStyle = css("--ink"); ctx.lineWidth = 5; ctx.strokeRect(rx, ry, rw, rh);              // halo
      ctx.strokeStyle = Math.abs(v) < .05 ? css("--mid") : css(v > 0 ? "--up" : "--down"); ctx.lineWidth = 3; ctx.strokeRect(rx, ry, rw, rh);
      if (c.want && p.BL * sx > 10) {
        const ok = bitOk(d, c, S.f), tx = x0 + p.BL * sx / 2, ty = h - (yc - hw) * sy + 16;
        ctx.font = `800 14px ${css("--font")}`; ctx.lineWidth = 3; ctx.strokeStyle = css("--lane");
        ctx.strokeText(ok ? "✓" : "✗", tx, ty); ctx.fillStyle = css(ok ? "--ok" : "--bad"); ctx.fillText(ok ? "✓" : "✗", tx, ty);
      }
    }
  }
  return `bit cells shifted ${fmt(lag)} nm upstream; solid: ${fmt(Tc, 0)} K band, ±${fmt(half)} nm` +
    (S.pol ? "; outline: sign of the cell's mean m_z inside the band (grey: |m_z| < 0.05)" : "");
}
// laser isotherms drawn on the film and keyed in the legend: freezing T solid, then dashed, dash-dot, dotted with falling T
function isotherms(d) {
  const p = d.model.p, Tc = d.e.Tc, dashes = [[10, 4], [6, 3, 1.5, 3], [1.5, 3]];
  const others = [600, 500, 400].filter(T => T > p.Tmin && T < p.Tmax && Math.abs(T - Tc) > 15);
  return [...(Tc > p.Tmin && Tc < p.Tmax ? [{ T: Tc, w: 2.4, dash: [] }] : []), ...others.map((T, i) => ({ T, w: 1.4, dash: dashes[i] }))];
}
// a line over the grain imagery: a wider stroke in the lane colour underneath keeps it readable on red and blue
function haloStroke(ctx, colour, w) {
  ctx.strokeStyle = css("--lane"); ctx.lineWidth = w + 3; ctx.globalAlpha = .85; ctx.stroke();
  ctx.strokeStyle = colour; ctx.lineWidth = w; ctx.globalAlpha = 1; ctx.stroke();
}
function drawHead(ctx, d, sx, sy, h) {
  const m = d.model, p = m.p, step = d.steps[S.f], hd = m.head(step);
  if (!hd.on) return;
  const px = x => x * sx, py = y => h - y * sy, ink = css("--ink");
  for (const { T, w, dash } of isotherms(d)) {
    const r = Math.sqrt(2 * Math.log((p.Tmax - p.Tmin) / (T - p.Tmin)));
    ctx.setLineDash(dash);
    ctx.beginPath(); ctx.ellipse(px(hd.x), py(hd.y), r * p.sigma_x * sx, r * p.sigma_y * sy, 0, 0, 2 * Math.PI); haloStroke(ctx, ink, w);
  }
  ctx.setLineDash([]);
  ctx.beginPath(); ctx.arc(px(hd.x), py(hd.y), 3, 0, 2 * Math.PI); haloStroke(ctx, ink, 1); ctx.fillStyle = ink; ctx.fill();
  const B = m.field(step) * p.field_uz, bx = hd.x - p.NPS, colour = css(B > 0 ? "--up" : B < 0 ? "--down" : "--muted");
  ctx.setLineDash(B === 0 ? [3, 3] : []);
  ctx.beginPath(); ctx.rect(px(bx - p.box_x / 2), py(hd.y + p.box_y / 2), p.box_x * sx, p.box_y * sy); haloStroke(ctx, colour, 2.4);
  ctx.setLineDash([]);
  const label = B === 0 ? "0 T" : `${B > 0 ? "↑" : "↓"} ${Math.abs(B).toFixed(2)} T`, ly = py(hd.y + p.box_y / 2) - 6;
  ctx.font = `600 12px ${css("--font")}`; ctx.textAlign = "center"; ctx.lineJoin = "round";
  ctx.strokeStyle = css("--lane"); ctx.lineWidth = 4; ctx.strokeText(label, px(bx), ly);
  ctx.fillStyle = colour; ctx.fillText(label, px(bx), ly); ctx.lineJoin = "miter";
}
function drawReaderBand(ctx, w, h, sy) {
  ctx.fillStyle = css("--ink"); ctx.globalAlpha = .1;
  ctx.fillRect(0, h - (S.ysec + S.rw / 2) * sy, w, S.rw * sy); ctx.globalAlpha = 1;
}

// pixel images (atom columns or one layer) at data resolution
function pixelImage(d, view, V, w, h) {
  const G = d.grid, vc = LUT.void;
  if (S.colour === "peak" && d.model) {
    const p = d.model.p;
    return filmImage(V, G.nx, G.ny, i => d.topCount[i] ? rgbOf(LUT.s, idxSeq((d.model.peakT((i % G.nx + .5) * G.pitch_nm, (Math.floor(i / G.nx) + .5) * G.pitch_nm) - p.Tmin) / (p.Tmax - p.Tmin))) : vc, w, h);
  }
  if (S.colour === "minm") return filmImage(V, G.nx, G.ny, i => d.topCount[i] ? rgbOf(LUT.s, idxSeq(1 - d.minm[i] / 255)) : vc, w, h);
  const s = nearestSnap(d, d.steps[S.f]);
  if (!s) return null;
  const agree = S.colour === "agree" && d.gmap;
  const col = (v, i) => {
    if (v !== v) return vc;
    if (agree) { const g = d.gmap[i], w = g >= 0 ? wantOf(d, g) : 0; return w ? rgbOf(LUT.agree, idxDiv(v * w)) : LUT.mid; }
    return rgbOf(LUT.m, idxDiv(v));
  };
  if (view === "atoms") {
    const a = snapArray(d, "top", s.index); if (!a) return null;
    return filmImage(V, G.nx, G.ny, i => col(snapVal(a, 2 * G.nx * G.ny + i), i), w, h);
  }
  const a = snapArray(d, "vol", s.index); if (!a) return null;
  const n = G.vx * G.vy, o = 2 * G.nz * n + S.layer * n, r = lp(G) / G.pitch_nm;
  return filmImage(V, G.vx, G.vy, i => col(snapVal(a, o + i), Math.min(G.ny - 1, Math.floor((Math.floor(i / G.vx) + .5) * r)) * G.nx + Math.min(G.nx - 1, Math.floor((i % G.vx + .5) * r))), w, h);
}
const lp = G => G.layer_pitch_nm || 2 * G.pitch_nm;
function snapLabel(d, s) {
  const dt = Math.abs(s.step - d.steps[S.f]) * d.e.dt;
  return `snapshot ${timeText(d, s.time)}` + (dt > 1e-13 ? ` (Δ ${fmt(dt * 1e12, 0)} ps)` : "");
}

function filmTip(x, y) {
  const d = cur(), e = d.e, g = grainAt(x, y);
  S.hov = g;
  const rows = [[`x, y`, `${fmt(x)}, ${fmt(y)} nm`]];
  if (d.snaps && S.view !== "grains") {
    const G = d.grid, ix = Math.floor(x / G.pitch_nm), iy = Math.floor(y / G.pitch_nm);
    const s = nearestSnap(d, d.steps[S.f]), a = s && snapArray(d, "top", s.index);
    if (a && ix >= 0 && ix < G.nx && iy >= 0 && iy < G.ny) {
      const v = snapVal(a, 2 * G.nx * G.ny + iy * G.nx + ix);
      rows.push(["column m<sub>z</sub>", v !== v ? "–" : fmt(v, 2)], ["atoms in column", d.topCount[iy * G.nx + ix]],
                ["lowest column |m|", fmt(d.minm[iy * G.nx + ix] / 255, 2)]);
    }
  }
  if (d.model) rows.push(["laser T now", `${fmt(d.model.T(x, y, d.steps[S.f]), 0)} K`], ["peak laser T", `${fmt(d.model.peakT(x, y), 0)} K`]);
  if (S.cells && S.pol && d.g && cellGeom(d)) {
    const C = cellGeom(d), p = C.p, t = Math.min(p.n_tracks - 1, Math.max(0, Math.floor((y - p.pad) / p.TW))), b = Math.floor((x + C.lag) / p.BL);
    const c = b >= 0 && b < p.bits_per_track ? cellParts(d)[t * p.bits_per_track + b] : null;
    if (c) rows.push([`bit ${b + 1} mean m<sub>z</sub>`, `${fmt(cellMz(d, c, S.f), 3)} ${c.want ? (bitOk(d, c, S.f) ? "✓" : "✗") : ""}`]);
  }
  if (g >= 0) {
    rows.push(["grain", e.grains.id[g]]);
    if (d.g) rows.push(["m<sub>z</sub>", fmt(gmz(d, S.f, g), 3)], ["|m|", fmt(gml(d, S.f, g), 3)]);
    if (d.an) rows.push(["outcome", outcome(d.an[g])]);
  }
  return kv(rows);
}
const CAT = { "full-field": "full field", "reversal": "field reversal", "zero-field": "field off", "unheated": "not melted", "undecided": "undecided" };
const outcome = a => a.category === "unheated" || a.category === "undecided" || !a.intended ? CAT[a.category]
  : `${a.correct ? "✓" : "✗"} ${CAT[a.category]}`;

// ---- cross-section lane (x–z through the section line, from the per-layer volume)
function xzLane() {
  return {
    id: "xz", group: "bottom", title: "x–z section",
    height: w => { const G = cur().grid; return Math.max(40, Math.round(w * G.nz * G.pitch_z_nm / S.dspan * 2)); },
    draw(ctx, w, h) {
      const d = cur(), G = d.grid;
      ctx.fillStyle = `rgb(${LUT.void})`; ctx.fillRect(0, 0, w, h);
      const s = nearestSnap(d, d.steps[S.f]), a = s ? snapArray(d, "vol", s.index) : null;
      if (!a) {   // still loading: keep showing the last section drawn instead of a blank frame
        if (this._last && this._last.st === d.st) { putImage(ctx, this._last.img, w, 0, h); setNote("xz", [`y ${fmt(S.ysec)} nm`, snapLabel(d, this._last.s), "loading next snapshot"]); }
        else setNote("xz", [s ? "loading" : "no snapshot"]);
        return;
      }
      const n = G.vx * G.vy, iy = Math.max(0, Math.min(G.vy - 1, Math.floor(S.ysec / lp(G))));
      const img = gridImage(G.vx, G.nz, i => { const v = snapVal(a, 2 * G.nz * n + Math.floor(i / G.vx) * n + iy * G.vx + (i % G.vx));
        return v !== v ? LUT.void : rgbOf(LUT.m, idxDiv(v)); }, 0, 1, w, h);
      putImage(ctx, img, w, 0, h);
      this._last = { img, s, st: d.st };
      setNote("xz", [`y ${fmt(S.ysec)} nm`, `z 0–${fmt(G.nz * G.pitch_z_nm)} nm, ×${+(h / (G.nz * G.pitch_z_nm) / (w / span())).toFixed(1)}`, snapLabel(d, s)]);
    },
    hover(x, y, ev, h) { return kv([["x", `${fmt(x)} nm`], ["z", `${fmt((1 - y / h) * cur().grid.nz * cur().grid.pitch_z_nm)} nm`]]); },
    yMode: "px",
  };
}

// ---- bits lane: intended bits shifted back by the freeze lag, with the per-bit outcome
const trackOfSection = p => Math.min(p.n_tracks - 1, Math.max(0, Math.floor((S.ysec - p.pad) / p.TW)));
// one row per track, stacked as on the film (last track on top)
function bitsLane() {
  const head = p => p.n_tracks > 1 ? 13 : 0, row = p => 40 + head(p);   // multi-track rows get a label line
  const MAX_ROWS = 8, shown = p => p.n_tracks > MAX_ROWS ? 4 : p.n_tracks;   // more tracks: 4 around the head or the section line
  const firstTrack = (d, p) => {
    if (p.n_tracks <= MAX_ROWS) return 0;
    const hd = d.model.head(d.steps[S.f]), t = hd.on ? hd.track : trackOfSection(p);
    return Math.max(0, Math.min(p.n_tracks - shown(p), t - 1));
  };
  const rowTrack = (d, p, py) => firstTrack(d, p) + Math.min(shown(p) - 1, Math.max(0, shown(p) - 1 - Math.floor(py / row(p))));
  const bitIndex = d => {                                               // melted grains of each bit
    if (d.bitIdx) return d.bitIdx;
    d.bitIdx = new Map();
    d.an.forEach((a, g) => { if (melted(a)) { const k = a.track * d.model.p.bits_per_track + a.bit; if (!d.bitIdx.has(k)) d.bitIdx.set(k, []); d.bitIdx.get(k).push(g); } });
    return d.bitIdx;
  };
  const bitGrains = (d, t, b) => bitIndex(d).get(t * d.model.p.bits_per_track + b) || [];
  // written once one of its grains has frozen: before that a grain "agrees" only because of the starting state
  const written = (d, f, gs) => f === d.nF - 1 || gs.some(g => d.an[g].decision_step != null && d.an[g].decision_step <= d.steps[f]);
  return {
    id: "bits", group: "top", title: "Bits", yMode: "px",
    height: () => row(cur().model.p) * shown(cur().model.p),
    draw(ctx, w, h) {
      const d = cur(), p = d.model.p, lag = d.e.freeze_lag_nm, bt = d.model.bit(d.steps[S.f]), f = S.f;
      ctx.clearRect(0, 0, w, h);
      ctx.font = `600 11.5px ${css("--font")}`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      const t0 = firstTrack(d, p);
      for (let t = t0; t < t0 + shown(p); t++) {
        const y0 = (shown(p) - 1 - (t - t0)) * row(p) + head(p);
        for (let b = 0; b < p.bits_per_track; b++) {
          const want = Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz), x0 = X(b * p.BL - lag, w), x1 = X((b + 1) * p.BL - lag, w);
          if (x1 < 0 || x0 > w) continue;
          ctx.fillStyle = want ? rgb(LUT.m, idxDiv(want * .45)) : css("--rule"); ctx.fillRect(x0 + .5, y0 + 2, x1 - x0 - 1, 18);
          const gs = bitGrains(d, t, b), ok = gs.filter(g => okAt(d, f, d.an[g], g)).length;
          if (gs.length && written(d, f, gs)) {
            const fw = (x1 - x0 - 1) * ok / gs.length;
            ctx.fillStyle = css("--ok"); ctx.fillRect(x0 + .5, y0 + 24, fw, 4);
            ctx.fillStyle = css("--bad"); ctx.fillRect(x0 + .5 + fw, y0 + 24, x1 - x0 - 1 - fw, 4);
            if (x1 - x0 > 26) { ctx.fillStyle = css("--muted"); ctx.fillText(`${ok}/${gs.length}`, (x0 + x1) / 2, y0 + 35); }
          } else if (gs.length) {
            ctx.fillStyle = css("--rule"); ctx.fillRect(x0 + .5, y0 + 24, x1 - x0 - 1, 4);
            if (x1 - x0 > 26) { ctx.globalAlpha = .45; ctx.fillStyle = css("--muted"); ctx.fillText(`${ok}/${gs.length}`, (x0 + x1) / 2, y0 + 35); ctx.globalAlpha = 1; }
          }
          if (x1 - x0 > 9) { ctx.fillStyle = css("--ink"); ctx.fillText(want > 0 ? "↑" : want < 0 ? "↓" : "0", (x0 + x1) / 2, y0 + 11); }
        }
        if (bt && bt[0] === t) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2; ctx.strokeRect(X(bt[1] * p.BL - lag, w) + 1, y0 + 1, X(p.BL, w) - 2, 20); }
        if (p.n_tracks > 1) {
          const label = `track ${t + 1}`;
          ctx.save(); ctx.font = `600 10.5px ${css("--font")}`; ctx.textAlign = "left";
          ctx.fillStyle = css(bt && bt[0] === t ? "--ink" : "--muted"); ctx.fillText(bt && bt[0] === t ? `${label} · writing` : label, 4, y0 - 6); ctx.restore();
        }
      }
      setNote("bits", [p.n_tracks > MAX_ROWS ? `tracks ${t0 + 1}–${t0 + shown(p)} of ${p.n_tracks}` : "", p.n_tracks > 1 ? `${p.n_tracks} tracks, as on the film; box: bit under the head` : "", `✓ / melted grains at ${timeText(d, d.time[S.f])}`, "grey: not yet written"]);
    },
    hover(x, py) {
      const d = cur(), p = d.model.p, lag = d.e.freeze_lag_nm, t = rowTrack(d, p, py), b = Math.floor((x + lag) / p.BL);
      if (b < 0 || b >= p.bits_per_track) return "";
      const want = Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz), gs = bitGrains(d, t, b);
      return kv([["bit", `${b + 1} / ${p.bits_per_track}`], p.n_tracks > 1 && ["track", t + 1], ["intended", want > 0 ? "up (+z)" : want < 0 ? "down (−z)" : "no field"],
        ["x", `${fmt(b * p.BL - lag)}–${fmt((b + 1) * p.BL - lag)} nm`], !written(d, S.f, gs) && ["written", "not yet"], [`✓ at ${timeText(d, d.time[S.f])}`, `${gs.filter(g => okAt(d, S.f, d.an[g], g)).length} / ${gs.length}`],
        ["✓ at end", `${gs.filter(g => d.an[g].correct).length} / ${gs.length}`]]);
    },
  };
}

// ---- read stage lanes
function dipoleLayer(d, iz) {
  const P = d.e.dipole, n = P.nx * P.ny, o = 2 * P.nz * n + iz * n;
  return d.B.subarray(o, o + n);
}
function bScale(arr) {
  const v = Array.from(arr).filter(isFinite).map(Math.abs).sort((a, b) => a - b);
  return v[Math.floor(v.length * .995)] || 1e-3;
}
const fieldColour = (P, val, sc) => i => rgbOf(LUT.b, idxDiv(val(Math.floor(i / (P.nx - 1)), i % (P.nx - 1)) / sc));
function sensorLane() {
  return {
    id: "bmap", group: "bottom", title: "Stray B<sub>z</sub>",
    height: filmH, scroll: true,
    draw(ctx, w, h) {
      const d = cur(), P = d.e.dipole, Bz = dipoleLayer(d, S.z), sc = bScale(Bz), V = filmScale(w, h);
      putFilm(ctx, filmImage(V, P.nx - 1, P.ny - 1, fieldColour(P, (iy, ix) => Bz[iy * P.nx + ix], sc), w, h, true), w);
      if (d.gmap) {   // grain boundaries (pixels without magnetic atoms) for orientation
        const M = d.e.grain_map, c = filmImage(V, M.nx, M.ny, i => d.gmap[i] < 0 ? [0, 0, 0] : [255, 255, 255], w, h);
        ctx.globalCompositeOperation = "multiply"; ctx.globalAlpha = .18; putFilm(ctx, c, w); ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";
      }
      ctx.save(); ctx.translate(-S.win.x0 * V.sx, V.y0 * V.sy);
      drawReaderBand(ctx, Lx() * V.sx, h, V.sy);
      ctx.restore();
      drawScrollBar(ctx, V, w, h);
      d.bscale = sc;
      setNote("bmap", [`z ${zText(S.z)}`, `±${fmt(sc * 1e3, 0)} mT (99.5 %)`, scrollNote(V)]);
    },
    hover(x, y) {
      const d = cur(), P = d.e.dipole, ix = Math.floor(x / cellNm(P, 0)), iy = Math.floor(y / cellNm(P, 1)), n = P.nx * P.ny;
      const B = c => d.B[c * P.nz * n + S.z * n + iy * P.nx + ix] * 1e3;
      return kv([["x, y, z", `${fmt(x)}, ${fmt(y)}, ${zText(S.z)}`], ["B<sub>x</sub>", `${fmt(B(0))} mT`], ["B<sub>y</sub>", `${fmt(B(1))} mT`], ["B<sub>z</sub>", `${fmt(B(2))} mT`]]);
    },
    click(x, y) { S.ysec = y; redraw(); },
    drag: dragSection,
  };
}
const cellNm = (P, a) => P.cell_nm_xyz ? P.cell_nm_xyz[a] : P.cell_nm;   // dipole cell size along x, y, z (nm)
const zText = iz => `${fmt((iz + .5) * cellNm(cur().e.dipole, 2), 1)} nm`;
function bxzLane() {
  return {
    id: "bxz", group: "bottom", title: "B<sub>z</sub> x–z",
    height: w => Math.round(w * cur().e.size_nm[2] / S.dspan * 2),
    draw(ctx, w, h) {
      const d = cur(), P = d.e.dipole, n = P.nx * P.ny, iy = Math.min(P.ny - 1, Math.floor(S.ysec / cellNm(P, 1)));
      putImage(ctx, gridImage(P.nx - 1, P.nz - 1, fieldColour(P, (iz, ix) => d.B[2 * P.nz * n + iz * n + iy * P.nx + ix], d.bscale || .05), 0, 1, w, h), w, 0, h, true);
      const yz = z => h - z / cur().e.size_nm[2] * h;
      ctx.strokeStyle = css("--ink"); ctx.setLineDash([4, 3]); ctx.globalAlpha = .6;
      for (const z of [10, (S.z + .5) * cellNm(P, 2)]) { ctx.beginPath(); ctx.moveTo(0, yz(z)); ctx.lineTo(w, yz(z)); ctx.stroke(); }
      ctx.setLineDash([]); ctx.globalAlpha = 1;
      setNote("bxz", [`y ${fmt(S.ysec)} nm`, `z 0–${fmt(d.e.size_nm[2], 0)} nm, ×${+(h / d.e.size_nm[2] / (w / span())).toFixed(1)}`]);
    },
    hover(x, y, ev, h) { return kv([["x", `${fmt(x)} nm`], ["z", `${fmt((1 - y / h) * cur().e.size_nm[2])} nm`]]); },
    yMode: "px",
  };
}
function waveLane() {
  return {
    id: "wave", group: "bottom", title: "Read-back",
    height: () => 160,
    draw(ctx, w, h) {
      const d = cur(), wv = waveform(d), lim = Math.max(...wv.map(Math.abs)) * 1.15 || 1e-3, P = d.e.dipole;
      ctx.clearRect(0, 0, w, h);
      const y = v => h / 2 - v / lim * (h / 2 - 8);
      ctx.strokeStyle = css("--rule"); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(0, y(0)); ctx.lineTo(w, y(0)); ctx.stroke();
      const wr = cellGeom(d) && (d.e.hamr ? d.e : S.M.byId[d.e.ref]);
      if (wr) {
        const p = wr.hamr, lag = wr.freeze_lag_nm, t = trackOfSection(p);
        ctx.strokeStyle = css("--muted"); ctx.setLineDash([2, 3]);
        for (let b = 1; b < p.bits_per_track; b++) if (p.sequence[t * p.bits_per_track + b] !== p.sequence[t * p.bits_per_track + b - 1]) {
          const xx = X(b * p.BL - lag, w); ctx.beginPath(); ctx.moveTo(xx, 4); ctx.lineTo(xx, h - 4); ctx.stroke();
        }
        ctx.setLineDash([]);
      }
      ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.6; ctx.beginPath();
      wv.forEach((v, ix) => { const xx = X((ix + .5) * cellNm(P, 0), w); ix ? ctx.lineTo(xx, y(v)) : ctx.moveTo(xx, y(v)); }); ctx.stroke();
      ctx.fillStyle = css("--muted"); ctx.font = `500 11px ${css("--font")}`; ctx.textAlign = "left";
      ctx.fillText(`+${fmt(lim * 1e3, 0)} mT`, 4, 12); ctx.fillText(`−${fmt(lim * 1e3, 0)} mT`, 4, h - 4);
      setNote("wave", [`mean B_z, y ${fmt(S.ysec)} ± ${fmt(S.rw / 2)} nm`, `z ${zText(S.z)}`, "dotted: intended transitions"]);
    },
    hover(x) { const d = cur(), wv = waveform(d), ix = Math.floor(x / cellNm(d.e.dipole, 0)); return kv([["x", `${fmt(x)} nm`], ["B<sub>z</sub>", `${fmt((wv[ix] || 0) * 1e3)} mT`]]); },
  };
}
function waveform(d) {
  const P = d.e.dipole, Bz = dipoleLayer(d, S.z), out = new Array(P.nx - 1).fill(0);
  let n = 0;
  for (let iy = 0; iy < P.ny - 1; iy++) {
    if (Math.abs((iy + .5) * cellNm(P, 1) - S.ysec) > S.rw / 2) continue;
    n++; for (let ix = 0; ix < P.nx - 1; ix++) out[ix] += Bz[iy * P.nx + ix];
  }
  return out.map(v => v / Math.max(n, 1));
}

// ------------------------------------------------------------------------------------------------ lane DOM
function setNote(id, parts) { const n = $("note-" + id); if (n) n.innerHTML = parts.filter(Boolean).map(p => `<span>${p}</span>`).join(""); }
function buildLanes() {
  const d = cur();
  for (const id of ["lanes-top", "lanes-bottom"]) $(id).querySelectorAll(".lane").forEach(n => n.remove());
  LANES = []; S.tip = null;
  if (d.e.hamr && d.an) LANES.push(bitsLane());
  LANES.push(filmLane(d.e.type === "dipole-read" ? "Written" : "Film"));
  if (d.snaps && d.e.snapshots.volume) LANES.push(xzLane());
  if (d.e.dipole) LANES.push(sensorLane(), bxzLane(), waveLane());
  for (const L of LANES) {
    const el = document.createElement("div"); el.className = "lane"; el.id = "lane-" + L.id;
    el.innerHTML = `<div class="lane-head">${L.title}</div><canvas aria-label="${L.title.replace(/<[^>]+>/g, "")}"></canvas><div class="note" id="note-${L.id}"></div>`;
    $("lanes-" + L.group).appendChild(el);
    L.canvas = el.querySelector("canvas");
    const pos = ev => { const r = L.canvas.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top, r.width, r.height]; };
    const xy = ev => { const [px, py, w, h] = pos(ev), V = filmScale(w, h); return [XI(px, w), L.yMode ? py : V.scrolled ? (h - py) / V.sy + V.y0 : (1 - py / h) * Ly(), h]; };
    let pan = null;                                                // shift + drag pans the window
    L.canvas.addEventListener("pointerdown", ev => { if (ev.shiftKey) { pan = { x: ev.clientX, y: ev.clientY, x0: S.win.x0, y0: S.y0 }; S.follow = false; L.canvas.setPointerCapture(ev.pointerId); } });
    L.canvas.addEventListener("pointerup", () => { pan = null; });
    L.canvas.addEventListener("wheel", ev => {
      const [px, , w, h] = pos(ev), V = filmScale(w, h);
      if (ev.ctrlKey) { ev.preventDefault(); S.follow = false; zoomAt(XI(px, w), Math.exp(ev.deltaY * .004)); return; }   // pinch, or ctrl + wheel
      if (Math.abs(ev.deltaX) > Math.abs(ev.deltaY)) { if (span() < Lx()) { ev.preventDefault(); S.follow = false; setWindow(S.win.x0 + ev.deltaX * span() / w, S.win.x1 + ev.deltaX * span() / w); } return; }
      if (L.scroll && V.scrolled) {                                // scroll the film inside the lane; at its ends the page scrolls
        const y = Math.max(0, Math.min(Ly() - (V.y1 - V.y0), V.y0 - ev.deltaY / V.sy));
        if (y !== V.y0) { ev.preventDefault(); S.follow = false; S.y0 = y; redraw(); }
      }
    }, { passive: false });
    L.canvas.addEventListener("mousemove", ev => {
      const [x, y, h] = xy(ev);
      if (pan) {
        const r = L.canvas.getBoundingClientRect(), n = span(), k = n / r.width;
        setWindow(pan.x0 - (ev.clientX - pan.x) * k, pan.x0 - (ev.clientX - pan.x) * k + n);
        const V = filmScale(r.width, r.height);
        S.y0 = Math.max(0, Math.min(Ly() - (V.y1 - V.y0), pan.y0 + (ev.clientY - pan.y) / V.sy));
        L.dragged = true;
        return;
      }
      if (L.drag && ev.buttons & 1) { L.dragged = true; L.drag(x, y); }
      S.tip = () => showTip(ev, L.hover ? L.hover(x, y, ev, h) : "");   // re-run by drawAll so the tip follows the timeline
      S.hx = x; placeXhair(ev.clientX); redraw();
    });
    L.canvas.addEventListener("mouseleave", () => { S.tip = null; hideTip(); S.hov = -1; S.hx = null; document.querySelectorAll(".xhair").forEach(x => x.hidden = true); redraw(); });
    if (L.drag) L.canvas.addEventListener("mousedown", ev => { L.dragged = false; ev.preventDefault(); });   // no text selection while dragging
    if (L.click) L.canvas.addEventListener("click", ev => { if (L.dragged) { L.dragged = false; return; } const [x, y, h] = xy(ev); L.click(x, y, ev, h); });
  }
  $("lanes-bottom").hidden = !LANES.some(L => L.group === "bottom");
}
function placeXhair(clientX) {
  for (const box of [$("lanes-top"), $("lanes-bottom")]) {
    const r = box.getBoundingClientRect(), xh = box.querySelector(".xhair");
    xh.hidden = box.hidden; xh.style.left = (clientX - r.left) + "px";
    xh.querySelector("span").textContent = `${fmt(S.hx)} nm`;
  }
}

// ------------------------------------------------------------------------------------------------ drawing
let _raf = 0;
function redraw() { if (!_raf) _raf = requestAnimationFrame(() => { _raf = 0; drawAll(); }); }
// a lane that cannot draw says why instead of staying blank; the others still draw
function laneError(L, ctx, w, h, err) {
  console.error(err);
  ctx.setTransform(dpr(), 0, 0, dpr(), 0, 0);
  ctx.fillStyle = css("--lane"); ctx.fillRect(0, 0, w, h);
  ctx.fillStyle = css("--bad"); ctx.font = `600 13px ${css("--font")}`; ctx.textAlign = "left";
  ctx.fillText(`This lane could not be drawn: ${err.message}`, 10, Math.min(20, h - 4));
  setNote(L.id, ["See the browser console for details"]);
}
function drawAll() {
  const d = cur();
  if (d.err) { $("data-error").hidden = false; $("data-error").textContent = `Part of the data could not be read: ${d.err}`; }
  else $("data-error").hidden = true;
  if (!frameReady(d, S.f)) { drawScrub(); drawClock(); return; }   // its block arrives shortly and redraws: keep the last picture
  followHead(d);
  for (const L of LANES) {
    const w = L.canvas.clientWidth; let g = null;
    if (!w) continue;                                              // not laid out yet
    try { g = fit(L.canvas, L.height(w)); L.draw(g.ctx, w, g.h); }
    catch (err) { if (g) laneError(L, g.ctx, w, g.h, err); else console.error(err); }
  }
  drawMinimap();
  drawAxis("axis-top"); if (!$("lanes-bottom").hidden) drawAxis("axis-bottom");
  drawScrub(); drawClock(); drawSt(); drawInspector();
  if (S.tip) S.tip();
}
function drawAxis(id) {
  const c = $(id), { ctx, w, h } = fit(c, 20), n = span();
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = css("--muted"); ctx.strokeStyle = css("--rule"); ctx.font = `500 11px ${css("--font")}`;
  const steps = n < Lx() - 1e-9 ? [.5, 1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000] : [5, 10, 20, 25, 50, 100, 200];
  const step = steps.find(s => w / (n / s) > 44) || (n > 1e4 ? 2000 : 500);
  ctx.textAlign = "center";
  for (let i = Math.ceil(S.win.x0 / step); i * step <= S.win.x1 + 1e-9; i++) { const x = +(i * step).toFixed(3), px = X(x, w); ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, 4); ctx.stroke(); ctx.fillText(x, Math.min(Math.max(px, 8), w - 12), 14); }
}
function timeText(d, t) {
  if (d.e.monte_carlo) return `MC step ${Math.round(t / d.e.dt).toLocaleString()}`;
  return t < 1e-9 ? `${fmt(t * 1e12, 0)} ps` : `${fmt(t * 1e9, 3)} ns`;
}
function drawClock() {
  const d = cur(), el = $("clock");
  if (d.nF <= 1) { el.innerHTML = ""; return; }
  const parts = [[`t`, `${timeText(d, d.time[S.f])} / ${timeText(d, d.time[d.nF - 1])}`]];
  if (d.model) {
    const step = d.steps[S.f], hd = d.model.head(step), bt = d.model.bit(step), B = d.model.field(step) * d.model.p.field_uz;
    parts.push(["laser x", hd.on ? `${fmt(hd.x)} nm` : "off"], ["bit", bt ? `${bt[1] + 1} / ${d.model.p.bits_per_track}` : "–"],
               ["B<sub>z</sub> head", `${B === 0 ? "0" : sgn(+B.toFixed(2))} T`]);
  }
  el.innerHTML = parts.map(([k, v]) => `<span>${k} <b>${v}</b></span>`).join("");
}

const SCRUB_ALL_ROWS = 16384;
// scrubber: system magnetisation under a draggable cursor, bit windows shaded
function drawScrub() {
  const d = cur(), c = $("scrub");
  $("transport").hidden = d.nF <= 1;
  if (d.nF <= 1) return;
  const { ctx, w, h } = fit(c, 46);
  ctx.clearRect(0, 0, w, h);
  const tx = f => f / (d.nF - 1) * w;
  if (d.model) {
    const p = d.model.p;
    for (let t = 0; t < p.n_tracks; t++) for (let b = 0; b < p.bits_per_track; b++) {
      const s0 = p.eq + t * p.nt + p.ne + b * p.nb + 1, x0 = tx(frameAt(d, s0)), x1 = tx(frameAt(d, s0 + p.nb - 1));
      ctx.fillStyle = rgb(LUT.m, idxDiv(Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz) * .3)); ctx.fillRect(x0, 0, x1 - x0, h);
    }
  }
  if (d.out) {
    const ci = n => d.cols.findIndex(c => c.endsWith(n)), iz = ci("M_norm_z"), il = ci("M_norm_l"), nc = d.cols.length, rows = d.out.length / nc;
    const tcol = d.cols.indexOf("Real_time"), t1 = d.time[d.nF - 1] || 1;
    const series = [];
    if (il >= 0) series.push({ v: r => d.out[r * nc + il], col: css("--muted"), lo: 0, hi: 1 });
    if (iz >= 0) series.push({ v: r => d.out[r * nc + iz] * (il >= 0 ? d.out[r * nc + il] : 1), col: css("--ink"), lo: -1, hi: 1 });
    for (const s of series) {
      ctx.strokeStyle = s.col; ctx.lineWidth = 1.5; ctx.beginPath();
      const yOf = v => h - 3 - (v - s.lo) / (s.hi - s.lo) * (h - 6), tOf = r => tcol >= 0 ? d.out[r * nc + tcol] : r / rows * t1;
      if (rows <= SCRUB_ALL_ROWS) {
        const stepR = Math.max(1, Math.floor(rows / (w * 2)));
        for (let r = 0; r < rows; r += stepR) {
          const yy = yOf(s.v(r));
          r ? ctx.lineTo(tOf(r) / t1 * w, yy) : ctx.moveTo(tOf(r) / t1 * w, yy);
        }
      } else {     // lowest and highest value in each pixel column, so a short event in a long run stays visible
        const lo = new Float64Array(Math.ceil(w) + 1).fill(Infinity), hi = new Float64Array(Math.ceil(w) + 1).fill(-Infinity);
        for (let r = 0; r < rows; r++) { const k = Math.max(0, Math.min(lo.length - 1, Math.floor(tOf(r) / t1 * w))), v = s.v(r); if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; }
        let first = true;
        for (let k = 0; k < lo.length; k++) if (hi[k] >= lo[k]) { first ? ctx.moveTo(k + .5, yOf(lo[k])) : ctx.lineTo(k + .5, yOf(lo[k])); ctx.lineTo(k + .5, yOf(hi[k])); first = false; }
      }
      ctx.stroke();
    }
  }
  const x = tx(S.f);
  ctx.fillStyle = css("--ink"); ctx.fillRect(x - 1, 0, 2, h);
  ctx.beginPath(); ctx.moveTo(x - 6, 0); ctx.lineTo(x + 6, 0); ctx.lineTo(x, 8); ctx.fill();
  c.setAttribute("aria-valuemin", 0); c.setAttribute("aria-valuemax", d.nF - 1); c.setAttribute("aria-valuenow", S.f);
  c.setAttribute("aria-valuetext", timeText(d, d.time[S.f]));
}

// ------------------------------------------------------------------------------------------------ space-time
// One row per grain (sorted by x, grouped by track), time left to right, as in the first explorer.
let _st = null;
const ST_MAX_FRAMES = 16384, ST_MAX_PIXELS = 64e6;
// Too many frames or grains for one image: one image pixel per screen pixel (per overview bin, if fewer). Each pixel holds
// the extreme of its frames and grains furthest from the pixel before it in time, so a short event stays visible.
function stResampled(d, rows, pw, ph) {
  const r = dpr(), nt = Math.min(Math.round(pw * r), d.ov ? d.ovB : d.nF), nr = Math.min(Math.round(ph * r), rows.length), c = off(nt, nr), x = c.getContext("2d"), im = x.createImageData(nt, nr);
  const agree = S.colour === "agree", N = d.ov ? d.ovB : d.nF;
  for (let py = 0; py < nr; py++) {
    const r0 = Math.floor(py * rows.length / nr), r1 = Math.max(r0 + 1, Math.floor((py + 1) * rows.length / nr));
    let prev = null;
    for (let px = 0; px < nt; px++) {
      const f0 = Math.floor(px * N / nt), f1 = Math.max(f0 + 1, Math.floor((px + 1) * N / nt));
      let lo = Infinity, hi = -Infinity, anyWant = !agree;
      for (let q = r0; q < r1; q++) {
        const g = rows[q], w = agree ? wantOf(d, g) : 1;
        if (!w) continue;
        anyWant = true;
        for (let f = f0; f < f1; f++) {
          let a, b;
          if (d.ov) { const o = (f * d.nG + g) * 4; a = d.ov[o] / 32767; b = d.ov[o + 1] / 32767; } else a = b = gmz(d, f, g);
          if (w < 0) { const t = a; a = -b; b = -t; }
          if (a < lo) lo = a; if (b > hi) hi = b;
        }
      }
      const v = !anyWant ? 0 : prev == null ? (lo + hi) / 2 : Math.abs(hi - prev) > Math.abs(lo - prev) ? hi : lo;
      prev = v;
      const col = !anyWant ? LUT.mid : rgbOf(agree ? LUT.agree : LUT.m, idxDiv(v)), o = (py * nt + px) * 4;
      im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
    }
  }
  x.putImageData(im, 0, 0);
  return c;
}
function stRows(d) {
  const an = d.an, e = d.e;
  let g = [...Array(d.nG).keys()];
  if (an && !S.stAll) g = g.filter(i => an[i].category !== "unheated");
  const tr = i => an ? an[i].track : 0;
  return g.sort((a, b) => tr(a) - tr(b) || e.grains.x[a] - e.grains.x[b]);
}
function drawSt() {
  const d = cur(), sec = $("st-panel");
  sec.hidden = !(d.g && d.nF > 1 && d.e.hamr);
  if (sec.hidden) return;
  const c = $("st"), rows = stRows(d), left = 44, top = 16, bot = 20;
  const rowH = Math.min(8, 520 / rows.length), H = top + bot + rowH * rows.length;
  const { ctx, w, h } = fit(c, H), pw = w - left - 4;
  const key = `${S.stage}|${S.stAll}|${S.colour}|${rows.length}|${document.documentElement.dataset.theme || ""}`;
  if (!_st || _st.key !== key) {
    const img = d.ov || d.nF > ST_MAX_FRAMES || d.nF * rows.length > ST_MAX_PIXELS ? stResampled(d, rows, pw, rowH * rows.length) : pixImage(d.nF, rows.length, i => { const f = i % d.nF, r = rows.length - 1 - Math.floor(i / d.nF), g = rows[r];
      if (S.colour === "agree") { const wv = wantOf(d, g); return wv ? rgbOf(LUT.agree, idxDiv(gmz(d, f, g) * wv)) : LUT.mid; }
      return rgbOf(LUT.m, idxDiv(gmz(d, f, g))); });
    _st = { key, img };
  }
  ctx.clearRect(0, 0, w, h);
  ctx.drawImage(_st.img, left, top, pw, rowH * rows.length);
  const t1 = d.steps[d.nF - 1], t0 = d.steps[0], tx = s => left + (s - t0) / (t1 - t0) * pw, ry = r => top + (r + .5) * rowH;
  if (d.model) {
    const p = d.model.p;
    for (let t = 0; t < p.n_tracks; t++) for (let b = 0; b < p.bits_per_track; b++) {   // bit windows along the top
      const s0 = p.eq + t * p.nt + p.ne + b * p.nb, want = Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz);
      ctx.fillStyle = css(want > 0 ? "--up" : want < 0 ? "--down" : "--rule"); ctx.fillRect(tx(s0) + .5, 4, tx(s0 + p.nb) - tx(s0) - 1, 6);
    }
    ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.1; ctx.setLineDash([4, 3]); ctx.beginPath();   // laser over each grain
    rows.forEach((g, r) => { const s = d.model.passStep(d.e.grains.x[g], d.an ? d.an[g].track : 0); r ? ctx.lineTo(tx(s), ry(r)) : ctx.moveTo(tx(s), ry(r)); });
    ctx.stroke(); ctx.setLineDash([]);
    if (d.an) {   // decision moments
      ctx.fillStyle = css("--ink");
      rows.forEach((g, r) => { const a = d.an[g]; if (a.decision_step != null) ctx.fillRect(tx(a.decision_step) - 1, ry(r) - rowH / 2, 2, Math.max(rowH, 1.5)); });
    }
  }
  if (d.an && d.model && d.model.p.n_tracks > 1) {   // track separators and names
    ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5; ctx.fillStyle = css("--ink"); ctx.font = `600 10.5px ${css("--font")}`; ctx.textAlign = "left";
    rows.forEach((g, r) => { if (r && d.an[g].track !== d.an[rows[r - 1]].track) { ctx.beginPath(); ctx.moveTo(left - 40, top + r * rowH); ctx.lineTo(left + pw, top + r * rowH); ctx.stroke(); } });
    for (let t = 0; t < d.model.p.n_tracks; t++) { const r = rows.findIndex(g => d.an[g].track === t); if (r >= 0) ctx.fillText(`track ${t + 1}`, left + 4, top + r * rowH + 11); }
  }
  ctx.fillStyle = css("--muted"); ctx.font = `500 10.5px ${css("--font")}`; ctx.textAlign = "right";
  const every = Math.ceil(14 / rowH);
  rows.forEach((g, r) => { if (r % every === 0) ctx.fillText(fmt(d.e.grains.x[g], 0), left - 5, ry(r) + 3.5); });
  const rs = rows.indexOf(S.sel);
  if (rs >= 0) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5; ctx.strokeRect(left - .5, ry(rs) - rowH / 2 - 1, pw + 1, rowH + 2); }
  ctx.textAlign = "center";
  const tEnd = d.time[d.nF - 1], nice = [1, 2, 5].flatMap(m => [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7].map(e => m * e)).sort((a, b) => a - b);
  const dtick = nice.find(v => v / tEnd * pw >= 80) || tEnd;          // labels at least 80 px apart
  for (let t = 0; t <= tEnd + 1e-15; t += dtick) ctx.fillText(timeText(d, t), Math.min(w - 20, tx(Math.round(t / d.e.dt))), h - 5);
  const xc = tx(d.steps[S.f]); ctx.fillStyle = css("--ink"); ctx.fillRect(xc - .75, top - 3, 1.5, rowH * rows.length + 6);
  c._map = { rows, left, top, rowH, pw, t0, t1 };
  $("st-meta").innerHTML = `<span>${rows.length} grains${S.stAll ? "" : " (melted)"}, sorted by x (nm)</span><span>dashed: laser centre passes</span><span>tick: decision</span>`;
}
function stPick(ev) {
  const c = $("st"), m = c._map, r = c.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
  const row = Math.floor((y - m.top) / m.rowH), s = m.t0 + (x - m.left) / m.pw * (m.t1 - m.t0);
  return { g: row >= 0 && row < m.rows.length ? m.rows[row] : -1, f: frameAt(cur(), s) };
}

// ------------------------------------------------------------------------------------------------ inspector
function drawInspector() {
  const d = cur(), e = d.e, el = $("inspector"), g = S.hov >= 0 ? S.hov : S.sel;
  if (g < 0) { el.innerHTML = `<h2>${stageTitle(S.stage)}</h2>` + kv(stageFacts(S.stage)) + (d.model ? legendHtml() : ""); return; }
  const a = d.an && d.an[g];
  let chip = "";
  if (a) chip = a.category === "unheated" ? `<span class="chip none">not melted</span>` : a.category === "undecided" ? `<span class="chip none">undecided</span>`
    : !a.intended ? `<span class="chip none">outside bits</span>` : a.correct ? `<span class="chip ok">✓ intended</span>` : `<span class="chip bad">✗ wrong</span>`;
  const rows = [["x, y", `${fmt(e.grains.x[g])}, ${fmt(e.grains.y[g])} nm`],
    ["magnetic atoms", (e.grains.magnetic_atoms || e.grains.atoms)[g].toLocaleString()]];
  if (d.g) rows.push(["m<sub>z</sub>", fmt(gmz(d, S.f, g), 3)], ["|m|", fmt(gml(d, S.f, g), 3)]);
  if (a) {
    if (a.intended) rows.push(["intended", a.intended > 0 ? "↑ +z" : "↓ −z"]);
    rows.push(["final m<sub>z</sub>", fmt(a.final, 3)]);
    if (a.decision_step != null) rows.push(["decided", timeText(d, a.decision_step * e.dt)], ["T then", `${fmt(a.T_dec, 0)} K`],
      ["B<sub>z</sub> then", `${sgn(+(a.B_dec * e.hamr.field_uz).toFixed(2))} T (${CAT[a.category]})`],
      ["behind laser", `${fmt(a.lag_dec)} nm`], ["off track centre", `${fmt(a.dy)} nm`]);
    rows.push(["peak T", `${fmt(a.peakT, 0)} K`], ["lowest |m|", fmt(a.minM, 3)]);
  }
  el.innerHTML = `<h2>Grain ${e.grains.id[g]} ${chip}</h2>` + kv(rows) +
    (d.g && d.nF > 1 ? `<svg id="story" role="img" aria-label="Grain history"></svg>` + (d.model ? `<label class="toggle"><input type="checkbox" id="story-all" ${S.storyAll ? "checked" : ""}> whole run</label>` : "") : "");
  if (d.g && d.nF > 1) drawStory(d, g, a);
  const cb = $("story-all"); if (cb) cb.onchange = () => { S.storyAll = cb.checked; redraw(); };
}
const stageTitle = st => S.M.byId[st].title;
function stageFacts(st) {
  const d = D[st], e = d.e;
  if (e.hamr && e.analysis) {
    const an = e.analysis.grains, heated = an.filter(melted), ok = heated.filter(a => a.correct).length;
    let now = 0; an.forEach((a, g) => { if (melted(a) && okAt(d, S.f, a, g)) now++; });
    return [[`✓ at ${timeText(d, d.time[S.f])}`, `${now} / ${heated.length} (${fmt(100 * now / heated.length, 0)} %)`],
      ["✓ at end", `${ok} / ${heated.length} (${fmt(100 * ok / heated.length, 0)} %)`],
      ...(() => { const cs = cellParts(d).filter(c => c.want && c.parts.length); if (!cs.length) return [];
        const n = f => cs.filter(c => bitOk(d, c, f)).length;
        return [[`bits ✓ at ${timeText(d, d.time[S.f])}`, `${n(S.f)} / ${cs.length}`], ["bits ✓ at end", `${n(d.nF - 1)} / ${cs.length}`]]; })(), ["freezing T (measured)", `${fmt(e.Tc, 0)} K`],
      ["freeze lag", `${fmt(e.freeze_lag_nm)} nm`], ["peak head field", `${fmt(e.hamr.Hpeak, 2)} T`], ["head speed", `${fmt(e.hamr.v, 0)} m/s`],
      ["bits", e.hamr.bits_per_track * e.hamr.n_tracks], ["bit length", `${fmt(e.hamr.BL)} nm`], ["laser T", `${fmt(e.hamr.Tmin, 0)}–${fmt(e.hamr.Tmax, 0)} K`]];
  }
  if (e.type === "monte-carlo") return [["integrator", e.input["sim:integrator"]], ["temperature", `${e.input["sim:temperature"] || "–"} K`],
    ["steps", d.nF > 1 ? Math.round(d.time[d.nF - 1] / e.dt).toLocaleString() : "–"], ["grains", d.nG]];
  const bits = d.g && cellGeom(d) ? cellParts(d).filter(c => c.want && c.parts.length) : [];
  const generic = e.type !== "dipole-read" && e.type !== "hamr";   // a kind of run with no specific views yet
  return [generic && ["program", e.input["sim:program"] || "–"], generic && ["integrator", e.input["sim:integrator"] || "–"],
    generic && ["temperature", `${e.input["sim:temperature"] || "–"} K`], generic && d.nF > 1 && ["steps", Math.round(d.time[d.nF - 1] / e.dt).toLocaleString()],
    generic && ["grains", d.nG], e.ref && !e.v1 && ["compared with", `${S.M.byId[e.ref].title} (${e.ref})`],
    bits.length && ["bits ✓ (written state)", `${bits.filter(c => bitOk(d, c, d.nF - 1)).length} / ${bits.length}`],
    (e.dipole || !generic) && ["dipole solver", e.input["dipole:solver"] || "–"], e.dipole && ["cells", `${e.dipole.nx - 1} × ${e.dipole.ny - 1} × ${e.dipole.nz - 1}`],
    e.dipole && ["cell size", e.dipole.cell_nm_xyz ? `${e.dipole.cell_nm_xyz.join(" x ")} nm` : `${e.dipole.cell_nm} nm`]];
}
function legendHtml() {
  const d = cur(), iso = isotherms(d).sort((a, b) => b.T - a.T).map(({ T, w, dash }) =>
    `<span><svg width="30" height="10" aria-hidden="true"><line x1="1" y1="5" x2="29" y2="5" stroke="currentColor" stroke-width="${w}" stroke-dasharray="${dash.join(" ")}"/></svg>${fmt(T, 0)} K${T === d.e.Tc ? " (freezing T)" : ""}</span>`).join("");
  return `<div class="legend-mini"><b>laser isotherms</b>${iso}<span><i class="box"></i>head-field box: red ↑, blue ↓, dashed 0 T</span></div>`;
}

// values of one grain at frames f0, f0 + step, ... up to f1, one short byte range per frame
const _samples = new Map();
function grainSamples(d, g, f0, f1, step) {
  const key = `${d.st}|${g}|${f0}|${f1}|${step}`;
  if (_samples.has(key)) return _samples.get(key);
  if (_samples.size > 64) _samples.delete(_samples.keys().next().value);
  const G = d.g, sz = G.T.BYTES_PER_ELEMENT, res = { ready: false, mz: new Map(), ml: new Map() }, jobs = [];
  for (let f = f0; f <= f1; f += step) jobs.push(getRange(G.path, (f * d.nG + g) * 4 * sz, 4 * sz).then(b => {
    const v = new G.T(b); res.ml.set(f, v[3] / d.gs); res.mz.set(f, v[2] / d.gs * v[3] / d.gs); }));
  _samples.set(key, res);
  Promise.all(jobs).then(() => { res.ready = true; redraw(); }).catch(err => { d.err = `grain series: ${err.message}`; _samples.delete(key); redraw(); });
  return res;
}
function drawStory(d, g, a) {
  const svg = $("story"), w = svg.clientWidth || 288, H = 290, NS = "http://www.w3.org/2000/svg";
  svg.setAttribute("viewBox", `0 0 ${w} ${H}`);
  const el = (tag, at) => { const n = document.createElementNS(NS, tag); for (const k in at) n.setAttribute(k, at[k]); svg.appendChild(n); return n; };
  let f0 = 0, f1 = d.nF - 1;
  const m = d.model, gx = d.e.grains.x[g], gy = d.e.grains.y[g];
  if (m && !S.storyAll) {                     // window around the grain's hottest moment
    let best = 0, bf = 0;
    for (let f = 0; f < d.nF; f += 2) { const T = m.T(gx, gy, d.steps[f]); if (T > best) { best = T; bf = f; } }
    const half = Math.round(3 * m.p.sigma_x / m.p.v * 1e-9 / (d.time[1] - d.time[0] || 1e-12));
    f0 = Math.max(0, bf - half); f1 = Math.min(d.nF - 1, bf + 2 * half);
  }
  // long series: only the values of this grain at the plotted frames are read from the file
  const pst = Math.max(1, Math.floor((f1 - f0) / 300));
  const sampled = d.ov ? grainSamples(d, g, f0, f1, pst) : null;
  const vmz = sampled ? f => sampled.mz.get(f) : f => gmz(d, f, g), vml = sampled ? f => sampled.ml.get(f) : f => gml(d, f, g);
  const panels = [];
  if (m) panels.push({ name: "T (K)", lo: m.p.Tmin, hi: m.p.Tmax, v: f => m.T(gx, gy, d.steps[f]), col: css("--ink"), ref: d.e.Tc });
  panels.push({ name: "|m|", lo: 0, hi: 1, v: vml, col: css("--muted"), series: true });
  panels.push({ name: "m<tspan baseline-shift='sub' font-size='9'>z</tspan>", lo: -1, hi: 1, v: vmz, col: css("--up"), ref: 0, series: true });
  if (m) panels.push({ name: "B<tspan baseline-shift='sub' font-size='9'>z</tspan> (T)", lo: -m.p.Hpeak, hi: m.p.Hpeak, v: f => { const hd = m.head(d.steps[f]); return m.inBox(gx, gy, hd) ? m.field(d.steps[f]) * m.p.field_uz : 0; }, col: css("--down"), ref: 0 });
  const L = 6, R = 6, top = 4, gap = 8, bot = 16, ph = (H - top - bot - gap * (panels.length - 1)) / panels.length;
  const tx = f => L + (f - f0) / Math.max(1, f1 - f0) * (w - L - R);
  panels.forEach((P, i) => {
    const y0 = top + i * (ph + gap), yv = v => y0 + ph - (v - P.lo) / (P.hi - P.lo) * ph;
    el("rect", { x: L, y: y0, width: w - L - R, height: ph, fill: "none", stroke: css("--rule") });
    if (P.ref != null) el("line", { x1: L, x2: w - R, y1: yv(P.ref), y2: yv(P.ref), stroke: css("--muted"), "stroke-dasharray": "3 3" });
    let dd = "";
    if (!(P.series && d.ov && !sampled.ready)) for (let f = f0; f <= f1; f += pst) dd += (f === f0 ? "M" : "L") + tx(f).toFixed(1) + " " + yv(P.v(f)).toFixed(1);
    el("path", { d: dd, fill: "none", stroke: P.col, "stroke-width": 1.5 });
    const t = el("text", { x: L + 5, y: y0 + 12, fill: css("--muted"), "font-size": 11, "font-family": css("--font") }); t.innerHTML = P.name;
  });
  const df = a && a.decision_step != null ? frameAt(d, a.decision_step) : null;
  if (df != null && df >= f0 && df <= f1) el("line", { x1: tx(df), x2: tx(df), y1: top, y2: H - bot, stroke: css("--ok"), "stroke-width": 1.2 });
  if (S.f >= f0 && S.f <= f1) el("line", { x1: tx(S.f), x2: tx(S.f), y1: top, y2: H - bot, stroke: css("--ink"), "stroke-width": 1 });
  for (const f of [f0, f1]) {
    const t = el("text", { x: tx(f), y: H - 3, fill: css("--muted"), "font-size": 11, "font-family": css("--font"), "text-anchor": f === f0 ? "start" : "end" });
    t.textContent = timeText(d, d.time[f]);
  }
  svg.onclick = ev => { const r = svg.getBoundingClientRect(); seek(Math.round(f0 + (ev.clientX - r.left - L) / (r.width - L - R) * (f1 - f0))); };
}

// ------------------------------------------------------------------------------------------------ sections
function drawAnalysis() {
  const sec = $("analysis"), d = cur(), e = d.e;
  if (!e.analysis) { sec.hidden = true; return; }
  sec.hidden = false;
  const an = e.analysis.grains, heated = an.filter(a => a.intended && a.category !== "unheated"), ok = heated.filter(a => a.correct).length;
  const cats = {}; an.forEach(a => { const k = a.category; cats[k] = cats[k] || { n: 0, ok: 0 }; cats[k].n++; if (a.correct) cats[k].ok++; });
  const th = e.analysis.thresholds;
  sec.innerHTML = `<h2>Write analysis</h2>
    <div class="stats">
      <div class="stat"><b>${ok} / ${heated.length}</b><span>melted grains ✓</span></div>
      <div class="stat"><b>${fmt(100 * ok / heated.length, 0)} %</b><span>success</span></div>
      <div class="stat"><b>${fmt(e.Tc, 0)} K</b><span>freezing T (measured)</span></div>
      <div class="stat"><b>${fmt(e.freeze_lag_nm)} nm</b><span>freeze lag</span></div>
      <div class="stat"><b>${fmt(e.hamr.Hpeak, 2)} T</b><span>peak head field</span></div>
    </div>
    <div class="figs">
      <div class="fig"><h3>Decision positions</h3><svg id="scatter" role="img" aria-label="Decision positions"></svg>
        <div class="figkey"><span><i style="background:var(--ok)"></i>✓</span><span><i style="background:var(--bad)"></i>✗</span><span><i class="line"></i>${fmt(e.Tc, 0)} K isotherm</span></div></div>
      <div class="fig"><h3>Outcome by field at decision</h3><div class="tablewrap"><table>
        <thead><tr><th>Field at decision</th><th class="num">Grains</th><th class="num">✓</th></tr></thead><tbody>
        ${Object.entries(cats).map(([k, v]) => `<tr><td>${CAT[k]}</td><td class="num">${v.n}</td><td class="num">${k === "unheated" || k === "undecided" ? "–" : v.ok}</td></tr>`).join("")}
        </tbody></table></div>
        ${kv([["decided", `m<sub>z</sub>|m| last below ${th.decision} in final direction`], ["melted", `|m| < ${th.heated} and lost direction`],
              ["undecided", `|final m<sub>z</sub>| < ${th.settled}`], ["reversal", `|B| < ${th.full_field} B<sub>peak</sub>`]])}</div>
    </div>`;
  drawScatter(d, an);
}
function drawScatter(d, an) {
  const svg = $("scatter"), w = svg.clientWidth || 480, h = 280, NS = "http://www.w3.org/2000/svg", e = d.e, p = e.hamr;
  svg.setAttribute("viewBox", `0 0 ${w} ${h}`); svg.style.height = h + "px";
  const el = (tag, at, txt) => { const n = document.createElementNS(NS, tag); for (const k in at) n.setAttribute(k, at[k]); if (txt != null) n.textContent = txt; svg.appendChild(n); return n; };
  const pts = an.filter(a => a.lag_dec != null && a.intended);
  const xmax = Math.max(5, ...pts.map(a => a.lag_dec)) * 1.1, ymax = Math.max(p.TW / 2, ...pts.map(a => Math.abs(a.dy))) * 1.05;
  const L = 44, R = 12, T = 10, B = 32, sx = v => L + v / xmax * (w - L - R), sy = v => T + (1 - (v + ymax) / (2 * ymax)) * (h - T - B);
  for (let v = 0; v <= xmax; v += xmax > 30 ? 10 : 5) { el("line", { x1: sx(v), x2: sx(v), y1: T, y2: h - B, stroke: css("--rule") }); el("text", { x: sx(v), y: h - B + 14, "text-anchor": "middle", fill: css("--muted"), "font-size": 11 }, v); }
  for (let v = -Math.floor(ymax / 10) * 10; v <= ymax; v += 10) { el("line", { x1: L, x2: w - R, y1: sy(v), y2: sy(v), stroke: css("--rule") }); el("text", { x: L - 6, y: sy(v) + 4, "text-anchor": "end", fill: css("--muted"), "font-size": 11 }, v); }
  el("text", { x: (L + w - R) / 2, y: h - 4, "text-anchor": "middle", fill: css("--muted"), "font-size": 11.5 }, "behind laser centre (nm)");
  el("text", { x: 12, y: (T + h - B) / 2, "text-anchor": "middle", fill: css("--muted"), "font-size": 11.5, transform: `rotate(-90 12 ${(T + h - B) / 2})` }, "off track centre (nm)");
  const k = e.Tc > p.Tmin && e.Tc < p.Tmax ? 2 * Math.log((p.Tmax - p.Tmin) / (e.Tc - p.Tmin)) : -1;   // no isotherm without a Tc
  let dd = "";
  for (let i = 0; i <= 120; i++) { const yy = -p.sigma_y * Math.sqrt(k) + i / 120 * 2 * p.sigma_y * Math.sqrt(k), r = k - yy * yy / p.sigma_y ** 2; if (!(r >= 0) || Math.abs(yy) > ymax) continue; dd += (dd ? "L" : "M") + sx(p.sigma_x * Math.sqrt(r)) + " " + sy(yy); }
  if (dd) el("path", { d: dd, fill: "none", stroke: css("--ink"), "stroke-width": 1.3, "stroke-dasharray": "5 3" });
  for (const a of pts) {
    const c = el("circle", { cx: sx(a.lag_dec), cy: sy(a.dy), r: 3.6, fill: css(a.correct ? "--ok" : "--bad"), "fill-opacity": .85, stroke: css("--lane"), "stroke-width": .8 });
    const g = an.indexOf(a);
    c.style.cursor = "pointer";
    c.addEventListener("mouseenter", ev => showTip(ev, kv([["grain", e.grains.id[g]], ["outcome", outcome(a)], ["behind laser", `${fmt(a.lag_dec)} nm`], ["off centre", `${fmt(a.dy)} nm`], ["T", `${fmt(a.T_dec, 0)} K`]])));
    c.addEventListener("mouseleave", hideTip);
    c.addEventListener("click", () => { S.sel = g; S.ysec = e.grains.y[g]; if (a.decision_step != null) seek(frameAt(d, a.decision_step)); window.scrollTo({ top: 0, behavior: "smooth" }); });
  }
}
function drawDetails() {
  const e = cur().e, sec = $("details"), ck = e.checks;
  const facts = [["source", e.source], ["reduced", `${S.M.created} with hamrvis ${S.M.hamrvis_version || "older than 0.9"}`], ["viewer", `hamrvis ${HAMRVIS_VERSION}`]];
  const from = e.grain_map && (e.grain_map.from_run || e.grain_map.from_stage);
  if (e.segments && e.segments.length > 1) facts.splice(1, 0, ["continued in", e.segments.slice(1).join(", ")]);
  for (const n of e.notes || []) facts.push(["note", n]);
  if (e.grain_map) facts.push(["grain shapes", `from atom coordinates${from ? ` (${from} ${e.grain_map.from_run ? "run" : "stage"})` : ""}, ${e.grain_map.mixed_columns} mixed columns`]);
  if (ck) {
    facts.push(["stored columns vs exact", `≤ ${fmt(ck.column_rounding_max, 4)}`]);
    if (ck.grain_txt_same_step && ck.grain_txt_same_step.n) facts.push(["grain.txt, same step", `${ck.grain_txt_same_step.n} snapshots, max diff ${ck.grain_txt_same_step.max_diff.toExponential(1)}`]);
    if (ck.grain_txt_one_step_earlier && ck.grain_txt_one_step_earlier.n) facts.push(["grain.txt, 1 step earlier", `${ck.grain_txt_one_step_earlier.n} snapshots, median ${ck.grain_txt_one_step_earlier.median_diff.toExponential(1)}, max ${ck.grain_txt_one_step_earlier.max_diff.toExponential(1)}`]);
  }
  const rows = Object.entries(e.input).map(([k, v]) => `<div><span>${k}</span><span>${v === "" ? "on" : v}${e.input_units[k] ? " " + e.input_units[k] : ""}</span></div>`).join("");
  sec.innerHTML = `<h2>Run details</h2><div class="detail-facts">${kv(facts)}</div><div class="params">${rows}</div>`;
}

// ------------------------------------------------------------------------------------------------ tooltip
function showTip(ev, html) {
  const t = $("tip"); if (!html) { t.hidden = true; return; }
  t.innerHTML = html; t.hidden = false;
  const r = t.getBoundingClientRect(), x = ev.clientX + 14 + r.width > innerWidth ? ev.clientX - r.width - 14 : ev.clientX + 14;
  t.style.left = x + "px"; t.style.top = Math.min(ev.clientY + 14, innerHeight - r.height - 8) + "px";
}
function hideTip() { $("tip").hidden = true; }

// ------------------------------------------------------------------------------------------------ controls
function seek(f) { const d = cur(); S.f = Math.max(0, Math.min(d.nF - 1, f)); S.fa = S.f; redraw(); }
function syncControls() {
  const d = cur();
  for (const [id, v] of [["v-grains", "grains"], ["v-atoms", "atoms"], ["v-layer", "layer"]]) {
    $(id).setAttribute("aria-pressed", S.view === v);
    $(id).disabled = v === "grains" ? !d.nG && !!d.snaps : !(d.snaps && (v === "atoms" || d.e.snapshots.volume));
  }
  const hasAn = !!(d.an || (refD(d) && refD(d).an));
  for (const [id, v] of [["c-mz", "mz"], ["c-agree", "agree"], ["c-peak", "peak"], ["c-minm", "minm"]]) {
    $(id).setAttribute("aria-pressed", S.colour === v);
    $(id).hidden = v === "agree" && !cellGeom(d);   // nothing written yet
    $(id).disabled = (v === "peak" && !d.model) || (v === "minm" && !(d.g || d.snaps)) || (v === "agree" && !hasAn);
  }
  $("layer-ctl").hidden = S.view !== "layer" || !d.snaps;
  if (d.grid) { $("layer").max = d.grid.nz - 1; $("layer-out").textContent = `z ${fmt((S.layer + .5) * d.grid.pitch_z_nm, 2)} nm`; }
  $("z-ctl").hidden = $("w-ctl").hidden = !d.e.dipole;
  if (d.e.dipole) { $("zsel").max = d.e.dipole.nz - 2; $("zsel").value = S.z; $("z-out").textContent = zText(S.z); $("w-out").textContent = `${S.rw} nm`; }
  $("st-all").checked = S.stAll;
  $("cells-ctl").hidden = !cellGeom(d);   // nothing written yet
  $("pol-ctl").hidden = $("cells-ctl").hidden || !S.cells || !d.g;
}
function bindControls() {
  [["v-grains", "grains"], ["v-atoms", "atoms"], ["v-layer", "layer"]].forEach(([id, v]) => $(id).onclick = () => { S.view = v; syncControls(); redraw(); });
  [["c-mz", "mz"], ["c-agree", "agree"], ["c-peak", "peak"], ["c-minm", "minm"]].forEach(([id, v]) => $(id).onclick = () => { S.colour = v; syncControls(); redraw(); });
  $("layer").oninput = ev => { S.layer = +ev.target.value; syncControls(); redraw(); };
  $("zsel").oninput = ev => { S.z = +ev.target.value; syncControls(); redraw(); };
  $("rw").oninput = ev => { S.rw = +ev.target.value; syncControls(); redraw(); };
  $("st-all").onchange = ev => { S.stAll = ev.target.checked; redraw(); };
  $("play").onclick = togglePlay;
  $("to-start").onclick = () => seek(0);
  $("to-end").onclick = () => seek(cur().nF - 1);
  $("cells").onchange = ev => { S.cells = ev.target.checked; syncControls(); redraw(); };
  $("pol").onchange = ev => { S.pol = ev.target.checked; redraw(); };
  $("theme").onclick = () => {   // header button; the artifact host may also set data-theme
    const root = document.documentElement, dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
    try { localStorage.setItem("hamrvis-theme", root.dataset.theme); } catch (e) { /* storage unavailable */ }
    themeChanged();
  };
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", themeChanged);
  const sc = $("scrub");
  const scrubTo = ev => { const r = sc.getBoundingClientRect(); seek(Math.round((ev.clientX - r.left) / r.width * (cur().nF - 1))); };
  sc.addEventListener("pointerdown", ev => { sc.setPointerCapture(ev.pointerId); scrubTo(ev); });
  sc.addEventListener("pointermove", ev => { if (ev.buttons) scrubTo(ev); });
  const st = $("st");
  st.addEventListener("mousemove", ev => { const { g, f } = stPick(ev), d = cur(); S.hov = g;
    showTip(ev, g >= 0 ? kv([["grain", d.e.grains.id[g]], ["x", `${fmt(d.e.grains.x[g])} nm`], ["t", timeText(d, d.time[f])], ["m<sub>z</sub>", fmt(gmz(d, f, g), 3)]]) : ""); redraw(); });
  st.addEventListener("mouseleave", () => { S.hov = -1; hideTip(); redraw(); });
  st.addEventListener("click", ev => { const { g, f } = stPick(ev); if (g >= 0) { S.sel = g; S.ysec = cur().e.grains.y[g]; } seek(f); });
  document.addEventListener("keydown", ev => {
    if (ev.target.closest("input, select, textarea")) return;
    const n = ev.shiftKey ? 50 : 1;
    if (ev.key === " ") { ev.preventDefault(); togglePlay(); }
    else if (ev.key === "ArrowRight") { ev.preventDefault(); seek(S.f + n); }
    else if (ev.key === "ArrowLeft") { ev.preventDefault(); seek(S.f - n); }
    else if (ev.key === "Home") { ev.preventDefault(); seek(0); }
    else if (ev.key === "End") { ev.preventDefault(); seek(cur().nF - 1); }
    else if (ev.key === "Escape") { S.sel = -1; redraw(); }
    else if (ev.key === "+" || ev.key === "=") $("zoom-in").click();
    else if (ev.key === "-") $("zoom-out").click();
    else if (ev.key === "0") $("zoom-fit").click();
  });
  $("zoom-in").onclick = () => { S.follow = false; zoomAt((S.win.x0 + S.win.x1) / 2, 1 / 1.6); };
  $("zoom-out").onclick = () => { S.follow = false; zoomAt((S.win.x0 + S.win.x1) / 2, 1.6); };
  $("zoom-fit").onclick = () => { S.follow = false; S.win = { x0: 0, x1: Lx() }; S.y0 = 0; redraw(); };
  $("follow").onchange = ev => { S.follow = ev.target.checked; redraw(); };
  const mini = $("minimap");
  const miniTo = ev => { const r = mini.getBoundingClientRect(), n = span(), x = (ev.clientX - r.left) / r.width * Lx(); S.follow = false; setWindow(x - n / 2, x + n / 2); };
  mini.addEventListener("pointerdown", ev => { mini.setPointerCapture(ev.pointerId); miniTo(ev); });
  mini.addEventListener("pointermove", ev => { if (ev.buttons) miniTo(ev); });
  $("save-png").hidden = !!EMB;
  $("save-png").onclick = savePng;
  window.addEventListener("resize", () => { if (isDefaultWindow()) fitWindow(); redraw(); });
}
function themeChanged() { buildLuts(); _st = null; redraw(); drawAnalysis(); }
let _last = 0;
function togglePlay() {
  S.playing = !S.playing; $("play").textContent = S.playing ? "❚❚" : "▶"; $("play").setAttribute("aria-label", S.playing ? "Pause" : "Play");
  if (S.playing) { if (S.f >= cur().nF - 1) S.f = S.fa = 0; _last = 0; requestAnimationFrame(tick); }
}
function tick(ts) {
  if (!S.playing) return;
  if (ts - _last > 33) {    // ~30 updates per second; speed = frames per update (fractions play in slow motion)
    _last = ts; const d = cur();
    S.fa += +$("speed").value; const f = Math.min(d.nF - 1, Math.floor(S.fa));
    if (f !== S.f) { S.f = f; prefetch(d); redraw(); }
    if (d.ov) frameBlock(d, Math.min(d.nF - 1, f + d.g.frames));   // the next block of a long series before playback reaches it
    if (S.f >= d.nF - 1) togglePlay();
  }
  requestAnimationFrame(tick);
}
// load the next snapshots the visible views will need before playback reaches them
function prefetch(d) {
  if (!d.snaps) return;
  const kinds = [];
  if (S.view === "atoms") kinds.push("top");
  if (S.view === "layer" || (d.e.snapshots.volume && LANES.some(L => L.id === "xz"))) kinds.push("vol");
  const step = d.steps[S.f];
  for (const s of d.snaps.filter(s => s.step > step).slice(0, 2)) for (const k of kinds) snapArray(d, k, s.index);
}
function savePng() {
  const r = dpr(), cs = LANES.map(L => L.canvas).concat([$("axis-top")]), head = 40 * r;
  const W = Math.max(...cs.map(c => c.width)), H = cs.reduce((s, c) => s + c.height + 6 * r, head);
  const out = off(W, H), x = out.getContext("2d");
  x.fillStyle = css("--bg"); x.fillRect(0, 0, W, H);
  x.fillStyle = css("--ink"); x.font = `700 ${18 * r}px ${css("--font")}`; x.fillText(`${S.M.name}: ${stageTitle(S.stage)}, ${timeText(cur(), cur().time[S.f])}`, 0, 24 * r);
  let y = head; for (const c of cs) { x.drawImage(c, 0, y); y += c.height + 6 * r; }
  const a = document.createElement("a"); a.href = out.toDataURL("image/png");
  a.download = `${S.M.name}-${S.stage}-${S.f}.png`; a.click();
}

// ------------------------------------------------------------------------------------------------ stages
async function setStage(st) {
  S.stage = st; S.sel = -1; S.hov = -1; _st = null;
  let d;
  try {
    d = await loadStage(st);
    if (d.e.ref) await loadStage(d.e.ref);                         // intended bits for "vs intended"
  } catch (err) {
    console.error(err);
    $("data-error").hidden = false; $("data-error").textContent = `This run could not be loaded: ${err.message}`;
    return;
  }
  const e = d.e;
  fitWindow(); S.follow = true;
  S.f = S.fa = e.hamr ? 0 : d.nF - 1;
  const C = cellGeom(d);
  S.ysec = C ? C.p.TW * .5 + C.p.pad : e.size_nm[1] / 2;
  if (e.dipole) S.z = Math.min(e.dipole.nz - 2, Math.round(e.size_nm[2] * .55 / cellNm(e.dipole, 2)));
  if (S.view !== "grains" && !d.snaps) S.view = "grains";
  if (S.view === "grains" && !d.nG && d.snaps) S.view = "atoms";      // no grains (continuous film): atom views only
  if ((S.colour === "peak" && !d.model) || (S.colour === "agree" && !C)) S.colour = "mz";
  document.querySelectorAll("#chain button").forEach(b => b.setAttribute("aria-pressed", b.dataset.stage === st));
  try { history.replaceState(null, "", "#" + st); } catch (err) { /* sandboxed */ }
  buildLanes(); syncControls(); drawAnalysis(); drawDetails(); redraw();
}
// one tab per run, in run order; the line under the title says what the run is
function runSub(e) {
  const T = e.input["sim:temperature"];
  if (e.hamr) return `${e.hamr.bits_per_track * e.hamr.n_tracks} bits, ${fmt(e.hamr.v, 0)} m/s`;
  if (e.type === "monte-carlo") return `${e.input["sim:integrator"] || "relax"}, ${T || ""} K`;
  if (e.dipole) return "dipole field";
  return e.v1 ? "1 step" : [e.input["sim:integrator"], T && `${T} K`].filter(Boolean).join(", ");
}
function buildChain() {
  const nav = $("chain");
  nav.innerHTML = "";
  for (const e of S.M.runs) {
    const b = document.createElement("button");
    b.dataset.stage = e.id;
    b.innerHTML = `<b>${e.title}</b><span>${e.v1 ? "" : e.id + " · "}${runSub(e)}</span>`;
    b.title = e.segments ? e.segments.join("\n") : e.source;
    b.onclick = () => setStage(e.id);
    nav.appendChild(b);
  }
}
// version 1 bundles had three fixed stages; give them the run list of version 2
function normalise(M) {
  if (!M.runs) {
    const T = { relax: "Relaxation", write: "Write", read: "Read-back" };
    M.runs = ["relax", "write", "read"].filter(s => M.stages[s]).map(s => Object.assign(M.stages[s], {
      id: s, title: T[s], v1: true, monte_carlo: s === "relax" && /monte-carlo/.test(M.stages[s].input["sim:integrator"] || ""),
      type: s === "write" ? "hamr" : s === "read" ? "dipole-read" : /monte-carlo/.test(M.stages[s].input["sim:integrator"] || "") ? "monte-carlo" : "generic",
      ref: s === "read" && M.stages.write ? "write" : undefined }));
  }
  M.byId = Object.fromEntries(M.runs.map(e => [e.id, e]));
  return M;
}

async function main() {
  try { const t = localStorage.getItem("hamrvis-theme"); if (t && !EMB) document.documentElement.dataset.theme = t; } catch (e) { /* storage unavailable */ }
  buildLuts(); bindControls();
  try {
    S.M = normalise(JSON.parse(new TextDecoder().decode(await getBuf("manifest.json"))));
  } catch (err) {
    $("run-sub").textContent = "";
    document.querySelector(".work").innerHTML = `<div class="error"><h2>No bundle</h2>${kv([["error", err.message], ["start with", "python3 -m hamrvis view &lt;bundle&gt;"]])}</div>`;
    return;
  }
  const runs = S.M.runs, any = runs.find(e => e.hamr) || runs[0];
  $("run-name").textContent = S.M.name;
  const nat = (runs.find(e => e.hamr && e.n_atoms) || runs.find(e => e.n_atoms) || {}).n_atoms;
  $("run-sub").textContent = [nat ? `${nat.toLocaleString()} atoms` : null, any.grains.id.length ? `${any.grains.id.length} grains` : "no grains",
    `${any.size_nm.map(v => fmt(v, 0)).join(" × ")} nm`].filter(Boolean).join(", ");
  document.title = `${S.M.name} HAMR run`;
  buildChain();
  const want = location.hash.slice(1);
  await setStage(S.M.byId[want] ? want : any.id);
  window.addEventListener("hashchange", () => { const h = location.hash.slice(1); if (S.M.byId[h] && h !== S.stage) setStage(h); });
}
main();
