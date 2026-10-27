// hamrvis viewer: loads a bundle made by `hamrvis reduce` (docs/bundle-format.md) and draws it.
// Data comes from bundle/<path> over HTTP (`hamrvis view`) or from window.HV_EMBED (`hamrvis export`).
"use strict";

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
const has = path => EMB ? path in EMB : true;
async function tryBuf(path) { try { return await getBuf(path); } catch (e) { return null; } }

// ------------------------------------------------------------------------------------------------ state
const S = { M: null, stage: null, f: 0, fa: 0, view: "grains", colour: "mz", layer: 0, sel: -1, hov: -1, ysec: 0,
            z: 0, rw: 10, playing: false, storyAll: false, stAll: false, hx: null };
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
  const head = step => {
    const k = step - p.eq, on = k >= 1 && k <= p.laser_steps, kk = Math.max(k, 1);
    const track = Math.floor((kk - 1) / p.nt), j = (kk - 1) % p.nt + 1;
    return { x: -p.BL / 2 - p.NPS + (j + 1) * p.step_nm, y: p.TW * (track + 0.5) + p.pad, on, track, j };
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
    return h.on && h.j - 1 >= p.ne && b < p.bits_per_track ? [h.track, b] : null;
  };
  const peakT = (x, y) => {                                       // hottest laser temperature ever at x, y
    const x0 = -p.BL / 2 - p.NPS + 2 * p.step_nm, x1 = -p.BL / 2 - p.NPS + (p.nt + 1) * p.step_nm;
    const dx = x < x0 ? x - x0 : x > x1 ? x - x1 : 0;
    let best = 0;
    for (let t = 0; t < p.n_tracks; t++) {
      const yt = p.TW * (t + 0.5) + p.pad;
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
}
const idxDiv = v => Math.max(0, Math.min(255, Math.round((v + 1) * 127.5)));     // -1..1
const idxSeq = v => Math.max(0, Math.min(255, Math.round(v * 255)));             // 0..1
const rgbOf = (L, i) => [L[3 * i], L[3 * i + 1], L[3 * i + 2]];
const rgb = (L, i) => `rgb(${L[3 * i]},${L[3 * i + 1]},${L[3 * i + 2]})`;

// ------------------------------------------------------------------------------------------------ loading
async function loadStage(st) {
  if (D[st]) return D[st];
  const e = S.M.stages[st], d = { e, st };
  d.nG = e.grains.id.length;
  if (e.grain_series) {
    const nF = e.grain_series.frames, nL = e.grain_series.lead.length;
    d.g = new Int16Array(await getBuf(`${st}/grains.i16`));
    d.lead = new Float64Array(await getBuf(`${st}/grain_lead.f64`));
    Object.assign(d, { nF, nL });
    d.time = Float64Array.from({ length: nF }, (_, f) => d.lead[f * nL]);
    d.steps = Float64Array.from(d.time, t => Math.round(t / e.dt));
  } else Object.assign(d, { nF: 1, time: [0], steps: [0] });
  const out = await tryBuf(`${st}/output.f32`);
  if (out) { d.out = new Float32Array(out); d.cols = e.output.columns; }
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
// one snapshot array, cached; returns null while loading and redraws when it arrives
function snapArray(d, kind, idx) {
  const cache = d[kind];
  if (cache.has(idx)) return cache.get(idx);
  const path = `${d.st}/${kind}/${String(idx).padStart(8, "0")}.i8.gz`;
  if (!has(path)) { cache.set(idx, undefined); return undefined; }
  cache.set(idx, null);
  getBuf(path).then(b => { cache.set(idx, new Int8Array(b)); if (cache.size > 80) cache.delete(cache.keys().next().value); redraw(); })
    .catch(() => cache.set(idx, undefined));
  return null;
}
function nearestSnap(d, step) {
  if (!d.snaps) return null;
  let best = null;
  for (const s of d.snaps) if (has(`${d.st}/top/${String(s.index).padStart(8, "0")}.i8.gz`) &&
    (!best || Math.abs(s.step - step) < Math.abs(best.step - step))) best = s;
  return best;
}

// ------------------------------------------------------------------------------------------------ grain helpers
const gmz = (d, f, g) => { const o = (f * d.nG + g) * 4; return d.g[o + 2] / 32767 * d.g[o + 3] / 32767; };
const gml = (d, f, g) => d.g[(f * d.nG + g) * 4 + 3] / 32767;
function frameAt(d, step) {                                       // last frame at or before step
  let lo = 0, hi = d.nF - 1;
  while (lo < hi) { const m = (lo + hi + 1) >> 1; if (d.steps[m] <= step) lo = m; else hi = m - 1; }
  return lo;
}
const _minGrain = {};
function minGrainM(d, g) {
  if (!_minGrain[d.st]) {
    const a = new Float32Array(d.nG).fill(1);
    for (let f = 0; f < d.nF; f++) for (let k = 0; k < d.nG; k++) a[k] = Math.min(a[k], gml(d, f, k));
    _minGrain[d.st] = a;
  }
  return _minGrain[d.st][g];
}
const wantOf = (d, g) => {                                        // intended z direction of grain g (write analysis);
  const an = d.an || (D.write && D.write.an);                     // 0 for grains the laser never melts
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
  const key = d.st + g;
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
  if (S.colour === "agree") { const w = wantOf(d, g); return w ? rgbOf(LUT.agree, idxDiv(gmz(d, f, g) * w)) : hex(css("--mid")); }
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
function blit(ctx, img, w, h, smooth = false) { ctx.imageSmoothingEnabled = smooth; ctx.drawImage(img, 0, 0, w, h); }
// image of nx x ny pixels from a colour function (y flipped: row 0 is the bottom of the film)
function pixImage(nx, ny, colour) {
  const c = off(nx, ny), x = c.getContext("2d"), im = x.createImageData(nx, ny);
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) {
    const col = colour(iy * nx + ix), o = ((ny - 1 - iy) * nx + ix) * 4;
    im.data[o] = col[0]; im.data[o + 1] = col[1]; im.data[o + 2] = col[2]; im.data[o + 3] = 255;
  }
  x.putImageData(im, 0, 0); return c;
}

// ------------------------------------------------------------------------------------------------ lanes
let LANES = [];
const cur = () => D[S.stage];
const Lx = () => cur().e.size_nm[0], Ly = () => cur().e.size_nm[1];
const X = (x, w) => x / Lx() * w;
const XI = (px, w) => px / w * Lx();

function filmLane(title) {
  return {
    id: "map", group: "top", title,
    height: w => Math.round(w * Ly() / Lx()),
    draw: drawFilm,
    hover: filmTip,
    click(x, y) { S.sel = grainAt(x, y); S.ysec = y; redraw(); },
  };
}

function drawFilm(ctx, w, h) {
  const d = cur(), e = d.e, sx = w / e.size_nm[0], sy = h / e.size_nm[1], M = e.grain_map;
  ctx.fillStyle = `rgb(${LUT.void})`; ctx.fillRect(0, 0, w, h);
  const view = d.snaps ? S.view : "grains";
  let note = "";
  if (view === "grains") {
    if (M) {
      const pal = Array.from({ length: d.nG }, (_, g) => grainColour(d, S.f, g));
      blit(ctx, pixImage(M.nx, M.ny, i => d.gmap[i] < 0 ? LUT.void : pal[d.gmap[i]]), w, h);
    } else if (e.grains.outline) {   // no atom coordinates in the bundle: VAMPIRE's Voronoi outlines (approximate)
      e.grains.outline.forEach((o, g) => { if (!o) return; const c = grainColour(d, S.f, g); ctx.fillStyle = `rgb(${c})`; ctx.beginPath();
        o.forEach(([x, y], i) => i ? ctx.lineTo(x * sx, h - y * sy) : ctx.moveTo(x * sx, h - y * sy)); ctx.fill(); });
      note = "approximate Voronoi outlines (no atom coordinates)";
    }
  } else {
    const img = pixelImage(d, view);
    if (img) blit(ctx, img, w, h);
    const s = nearestSnap(d, d.steps[S.f]);
    if (s) note = snapLabel(d, s);
  }
  // hovered and selected grain outlines, from the atom map
  if (M) for (const [g, col, lw] of [[S.hov, "--muted", 1], [S.sel, "--ink", 1]]) {
    if (g < 0) continue;
    ctx.fillStyle = css(col); const pw = w / M.nx * (lw + .2), ph = h / M.ny * (lw + .2);
    for (const i of grainEdge(d, g)) ctx.fillRect((i % M.nx) * w / M.nx, h - (Math.floor(i / M.nx) + 1) * h / M.ny, pw, ph);
  }
  ctx.setLineDash([5, 4]); ctx.lineWidth = 1; ctx.strokeStyle = css("--ink"); ctx.globalAlpha = .5;
  const yl = h - S.ysec * sy; ctx.beginPath(); ctx.moveTo(0, yl); ctx.lineTo(w, yl); ctx.stroke();
  ctx.setLineDash([]); ctx.globalAlpha = 1;
  if (d.model) drawHead(ctx, d, sx, sy, h);
  if (S.stage === "read") drawReaderBand(ctx, w, h, sy);
  setNote("map", [note, `section y ${fmt(S.ysec)} nm`]);
}
function drawHead(ctx, d, sx, sy, h) {
  const m = d.model, p = m.p, step = d.steps[S.f], hd = m.head(step);
  if (!hd.on) return;
  const px = x => x * sx, py = y => h - y * sy, Tc = d.e.Tc;
  ctx.strokeStyle = css("--ink");
  for (const T of [400, 500, 600, Tc].filter(T => T > p.Tmin && T < p.Tmax)) {
    const r = Math.sqrt(2 * Math.log((p.Tmax - p.Tmin) / (T - p.Tmin)));
    ctx.lineWidth = T === Tc ? 2.2 : 1; ctx.globalAlpha = T === Tc ? .95 : .45;
    ctx.beginPath(); ctx.ellipse(px(hd.x), py(hd.y), r * p.sigma_x * sx, r * p.sigma_y * sy, 0, 0, 2 * Math.PI); ctx.stroke();
  }
  ctx.globalAlpha = 1;
  ctx.fillStyle = css("--ink"); ctx.beginPath(); ctx.arc(px(hd.x), py(hd.y), 3, 0, 2 * Math.PI); ctx.fill();
  const B = m.field(step) * p.field_uz, bx = hd.x - p.NPS;
  ctx.strokeStyle = css(B > 0 ? "--up" : B < 0 ? "--down" : "--muted"); ctx.lineWidth = 2.2; ctx.setLineDash(B === 0 ? [3, 3] : []);
  ctx.strokeRect(px(bx - p.box_x / 2), py(hd.y + p.box_y / 2), p.box_x * sx, p.box_y * sy);
  ctx.setLineDash([]);
  ctx.fillStyle = ctx.strokeStyle; ctx.font = `600 12px ${css("--font")}`; ctx.textAlign = "center";
  ctx.fillText(B === 0 ? "0 T" : `${B > 0 ? "↑" : "↓"} ${Math.abs(B).toFixed(2)} T`, px(bx), py(hd.y + p.box_y / 2) - 5);
}
function drawReaderBand(ctx, w, h, sy) {
  ctx.fillStyle = css("--ink"); ctx.globalAlpha = .1;
  ctx.fillRect(0, h - (S.ysec + S.rw / 2) * sy, w, S.rw * sy); ctx.globalAlpha = 1;
}

// pixel images (atom columns or one layer) at data resolution
function pixelImage(d, view) {
  const G = d.grid, vc = LUT.void;
  if (S.colour === "peak" && d.model) {
    const p = d.model.p;
    return pixImage(G.nx, G.ny, i => d.topCount[i] ? rgbOf(LUT.s, idxSeq((d.model.peakT((i % G.nx + .5) * G.pitch_nm, (Math.floor(i / G.nx) + .5) * G.pitch_nm) - p.Tmin) / (p.Tmax - p.Tmin))) : vc);
  }
  if (S.colour === "minm") return pixImage(G.nx, G.ny, i => d.topCount[i] ? rgbOf(LUT.s, idxSeq(1 - d.minm[i] / 255)) : vc);
  const s = nearestSnap(d, d.steps[S.f]);
  if (!s) return null;
  const agree = S.colour === "agree" && d.gmap;
  const col = (v, i) => {
    if (v === -128) return vc;
    if (agree) { const g = d.gmap[i], w = g >= 0 ? wantOf(d, g) : 0; return w ? rgbOf(LUT.agree, idxDiv(v / 127 * w)) : hex(css("--mid")); }
    return rgbOf(LUT.m, idxDiv(v / 127));
  };
  if (view === "atoms") {
    const a = snapArray(d, "top", s.index); if (!a) return null;
    return pixImage(G.nx, G.ny, i => col(a[2 * G.nx * G.ny + i], i));
  }
  const a = snapArray(d, "vol", s.index); if (!a) return null;
  const n = G.vx * G.vy, o = 2 * G.nz * n + S.layer * n;
  return pixImage(G.vx, G.vy, i => col(a[o + i], ((Math.floor(i / G.vx) * 2) * G.nx + (i % G.vx) * 2)));
}
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
      const v = a[2 * G.nx * G.ny + iy * G.nx + ix];
      rows.push(["column m<sub>z</sub>", v === -128 ? "–" : fmt(v / 127, 2)], ["atoms in column", d.topCount[iy * G.nx + ix]],
                ["lowest column |m|", fmt(d.minm[iy * G.nx + ix] / 255, 2)]);
    }
  }
  if (d.model) rows.push(["laser T now", `${fmt(d.model.T(x, y, d.steps[S.f]), 0)} K`], ["peak laser T", `${fmt(d.model.peakT(x, y), 0)} K`]);
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
    height: w => { const G = cur().grid; return Math.max(40, Math.round(w * G.nz * G.pitch_z_nm / Lx() * 2)); },
    draw(ctx, w, h) {
      const d = cur(), G = d.grid;
      ctx.fillStyle = `rgb(${LUT.void})`; ctx.fillRect(0, 0, w, h);
      const s = nearestSnap(d, d.steps[S.f]), a = s ? snapArray(d, "vol", s.index) : null;
      if (!a) { setNote("xz", [s ? "loading" : "no snapshot"]); return; }
      const n = G.vx * G.vy, iy = Math.max(0, Math.min(G.vy - 1, Math.floor(S.ysec / (2 * G.pitch_nm))));
      blit(ctx, pixImage(G.vx, G.nz, i => { const v = a[2 * G.nz * n + Math.floor(i / G.vx) * n + iy * G.vx + (i % G.vx)];
        return v === -128 ? LUT.void : rgbOf(LUT.m, idxDiv(v / 127)); }), w, h);
      setNote("xz", [`y ${fmt(S.ysec)} nm`, `z 0–${fmt(G.nz * G.pitch_z_nm)} nm, ×2`, snapLabel(d, s)]);
    },
    hover(x, y, ev, h) { return kv([["x", `${fmt(x)} nm`], ["z", `${fmt((1 - y / h) * cur().grid.nz * cur().grid.pitch_z_nm)} nm`]]); },
    yMode: "px",
  };
}

// ---- bits lane: intended bits shifted back by the freeze lag, with the per-bit outcome
const trackOfSection = p => Math.min(p.n_tracks - 1, Math.max(0, Math.floor((S.ysec - p.pad) / p.TW)));
function bitsLane() {
  return {
    id: "bits", group: "top", title: "Bits",
    height: () => 40,
    draw(ctx, w, h) {
      const d = cur(), p = d.model.p, lag = d.e.freeze_lag_nm, t = trackOfSection(p);
      ctx.clearRect(0, 0, w, h);
      ctx.font = `600 11.5px ${css("--font")}`; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      for (let b = 0; b < p.bits_per_track; b++) {
        const want = Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz), x0 = X(b * p.BL - lag, w), x1 = X((b + 1) * p.BL - lag, w);
        ctx.fillStyle = rgb(LUT.m, idxDiv(want * .45)); ctx.fillRect(x0 + .5, 2, x1 - x0 - 1, 18);
        const gs = d.an.filter(a => a.track === t && a.bit === b && a.category !== "unheated"), ok = gs.filter(a => a.correct).length;
        if (gs.length) {
          const fw = (x1 - x0 - 1) * ok / gs.length;
          ctx.fillStyle = css("--ok"); ctx.fillRect(x0 + .5, 24, fw, 4);
          ctx.fillStyle = css("--bad"); ctx.fillRect(x0 + .5 + fw, 24, x1 - x0 - 1 - fw, 4);
          if (x1 - x0 > 26) { ctx.fillStyle = css("--muted"); ctx.fillText(`${ok}/${gs.length}`, (x0 + x1) / 2, 35); }
        }
        ctx.fillStyle = css("--ink"); ctx.fillText(want > 0 ? "↑" : "↓", (x0 + x1) / 2, 11);
      }
      const bt = d.model.bit(d.steps[S.f]);
      if (bt && bt[0] === t) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 2; ctx.strokeRect(X(bt[1] * p.BL - lag, w) + 1, 1, X(p.BL, w) - 2, 20); }
      setNote("bits", []);
    },
    hover(x) {
      const d = cur(), p = d.model.p, lag = d.e.freeze_lag_nm, t = trackOfSection(p), b = Math.floor((x + lag) / p.BL);
      if (b < 0 || b >= p.bits_per_track) return "";
      const want = Math.sign(p.sequence[t * p.bits_per_track + b] * p.field_uz);
      const gs = d.an.filter(a => a.track === t && a.bit === b && a.category !== "unheated");
      return kv([["bit", `${b + 1} / ${p.bits_per_track}`], p.n_tracks > 1 && ["track", t + 1], ["intended", want > 0 ? "up (+z)" : "down (−z)"],
        ["x", `${fmt(b * p.BL - lag)}–${fmt((b + 1) * p.BL - lag)} nm`], ["melted grains ✓", `${gs.filter(a => a.correct).length} / ${gs.length}`]]);
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
function fieldImage(P, nrows, val, sc) {
  return pixImage(P.nx - 1, nrows, i => rgbOf(LUT.b, idxDiv(val(Math.floor(i / (P.nx - 1)), i % (P.nx - 1)) / sc)));
}
function sensorLane() {
  return {
    id: "bmap", group: "bottom", title: "Stray B<sub>z</sub>",
    height: w => Math.round(w * Ly() / Lx()),
    draw(ctx, w, h) {
      const d = cur(), P = d.e.dipole, Bz = dipoleLayer(d, S.z), sc = bScale(Bz);
      blit(ctx, fieldImage(P, P.ny - 1, (iy, ix) => Bz[iy * P.nx + ix], sc), w, h, true);
      if (d.gmap) {   // grain boundaries (pixels without magnetic atoms) for orientation
        const M = d.e.grain_map, c = pixImage(M.nx, M.ny, i => d.gmap[i] < 0 ? [0, 0, 0] : [255, 255, 255]);
        ctx.globalCompositeOperation = "multiply"; ctx.globalAlpha = .18; blit(ctx, c, w, h); ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";
      }
      drawReaderBand(ctx, w, h, h / Ly());
      d.bscale = sc;
      setNote("bmap", [`z ${zText(S.z)}`, `±${fmt(sc * 1e3, 0)} mT (99.5 %)`]);
    },
    hover(x, y) {
      const d = cur(), P = d.e.dipole, ix = Math.floor(x / P.cell_nm), iy = Math.floor(y / P.cell_nm), n = P.nx * P.ny;
      const B = c => d.B[c * P.nz * n + S.z * n + iy * P.nx + ix] * 1e3;
      return kv([["x, y, z", `${fmt(x)}, ${fmt(y)}, ${zText(S.z)}`], ["B<sub>x</sub>", `${fmt(B(0))} mT`], ["B<sub>y</sub>", `${fmt(B(1))} mT`], ["B<sub>z</sub>", `${fmt(B(2))} mT`]]);
    },
    click(x, y) { S.ysec = y; redraw(); },
  };
}
const zText = iz => `${fmt((iz + .5) * cur().e.dipole.cell_nm, 1)} nm`;
function bxzLane() {
  return {
    id: "bxz", group: "bottom", title: "B<sub>z</sub> x–z",
    height: w => Math.round(w * cur().e.size_nm[2] / Lx() * 2),
    draw(ctx, w, h) {
      const d = cur(), P = d.e.dipole, n = P.nx * P.ny, iy = Math.min(P.ny - 1, Math.floor(S.ysec / P.cell_nm));
      blit(ctx, fieldImage(P, P.nz - 1, (iz, ix) => d.B[2 * P.nz * n + iz * n + iy * P.nx + ix], d.bscale || .05), w, h, true);
      const yz = z => h - z / cur().e.size_nm[2] * h;
      ctx.strokeStyle = css("--ink"); ctx.setLineDash([4, 3]); ctx.globalAlpha = .6;
      for (const z of [10, (S.z + .5) * P.cell_nm]) { ctx.beginPath(); ctx.moveTo(0, yz(z)); ctx.lineTo(w, yz(z)); ctx.stroke(); }
      ctx.setLineDash([]); ctx.globalAlpha = 1;
      setNote("bxz", [`y ${fmt(S.ysec)} nm`, `z 0–${fmt(d.e.size_nm[2], 0)} nm, ×2`]);
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
      const wr = S.M.stages.write;
      if (wr && wr.hamr) {
        const p = wr.hamr, lag = wr.freeze_lag_nm, t = trackOfSection(p);
        ctx.strokeStyle = css("--muted"); ctx.setLineDash([2, 3]);
        for (let b = 1; b < p.bits_per_track; b++) if (p.sequence[t * p.bits_per_track + b] !== p.sequence[t * p.bits_per_track + b - 1]) {
          const xx = X(b * p.BL - lag, w); ctx.beginPath(); ctx.moveTo(xx, 4); ctx.lineTo(xx, h - 4); ctx.stroke();
        }
        ctx.setLineDash([]);
      }
      ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.6; ctx.beginPath();
      wv.forEach((v, ix) => { const xx = X((ix + .5) * P.cell_nm, w); ix ? ctx.lineTo(xx, y(v)) : ctx.moveTo(xx, y(v)); }); ctx.stroke();
      ctx.fillStyle = css("--muted"); ctx.font = `500 11px ${css("--font")}`; ctx.textAlign = "left";
      ctx.fillText(`+${fmt(lim * 1e3, 0)} mT`, 4, 12); ctx.fillText(`−${fmt(lim * 1e3, 0)} mT`, 4, h - 4);
      setNote("wave", [`mean B_z, y ${fmt(S.ysec)} ± ${fmt(S.rw / 2)} nm`, `z ${zText(S.z)}`, "dotted: intended transitions"]);
    },
    hover(x) { const d = cur(), wv = waveform(d), ix = Math.floor(x / d.e.dipole.cell_nm); return kv([["x", `${fmt(x)} nm`], ["B<sub>z</sub>", `${fmt((wv[ix] || 0) * 1e3)} mT`]]); },
  };
}
function waveform(d) {
  const P = d.e.dipole, Bz = dipoleLayer(d, S.z), out = new Array(P.nx - 1).fill(0);
  let n = 0;
  for (let iy = 0; iy < P.ny - 1; iy++) {
    if (Math.abs((iy + .5) * P.cell_nm - S.ysec) > S.rw / 2) continue;
    n++; for (let ix = 0; ix < P.nx - 1; ix++) out[ix] += Bz[iy * P.nx + ix];
  }
  return out.map(v => v / Math.max(n, 1));
}

// ------------------------------------------------------------------------------------------------ lane DOM
function setNote(id, parts) { const n = $("note-" + id); if (n) n.innerHTML = parts.filter(Boolean).map(p => `<span>${p}</span>`).join(""); }
function buildLanes() {
  const d = cur();
  for (const id of ["lanes-top", "lanes-bottom"]) $(id).querySelectorAll(".lane").forEach(n => n.remove());
  LANES = [];
  if (S.stage === "write" && d.an) LANES.push(bitsLane());
  LANES.push(filmLane(S.stage === "read" ? "Written" : "Film"));
  if (d.snaps && d.e.snapshots.volume) LANES.push(xzLane());
  if (S.stage === "read" && d.e.dipole) LANES.push(sensorLane(), bxzLane(), waveLane());
  for (const L of LANES) {
    const el = document.createElement("div"); el.className = "lane"; el.id = "lane-" + L.id;
    el.innerHTML = `<div class="lane-head">${L.title}</div><canvas aria-label="${L.title.replace(/<[^>]+>/g, "")}"></canvas><div class="note" id="note-${L.id}"></div>`;
    $("lanes-" + L.group).appendChild(el);
    L.canvas = el.querySelector("canvas");
    const pos = ev => { const r = L.canvas.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top, r.width, r.height]; };
    const xy = ev => { const [px, py, w, h] = pos(ev); return [XI(px, w), L.yMode ? py : (1 - py / h) * Ly(), h]; };
    L.canvas.addEventListener("mousemove", ev => {
      const [x, y, h] = xy(ev);
      showTip(ev, L.hover ? L.hover(x, y, ev, h) : "");
      S.hx = x; placeXhair(ev.clientX); redraw();
    });
    L.canvas.addEventListener("mouseleave", () => { hideTip(); S.hov = -1; S.hx = null; document.querySelectorAll(".xhair").forEach(x => x.hidden = true); redraw(); });
    if (L.click) L.canvas.addEventListener("click", ev => { const [x, y, h] = xy(ev); L.click(x, y, ev, h); });
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
function drawAll() {
  for (const L of LANES) { const w = L.canvas.clientWidth; const { ctx, h } = fit(L.canvas, L.height(w)); L.draw(ctx, w, h); }
  drawAxis("axis-top"); if (!$("lanes-bottom").hidden) drawAxis("axis-bottom");
  drawScrub(); drawClock(); drawSt(); drawInspector();
}
function drawAxis(id) {
  const c = $(id), { ctx, w, h } = fit(c, 20), L = Lx();
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = css("--muted"); ctx.strokeStyle = css("--rule"); ctx.font = `500 11px ${css("--font")}`;
  const step = [5, 10, 20, 25, 50, 100, 200].find(s => w / (L / s) > 44) || 500;
  ctx.textAlign = "center";
  for (let x = 0; x <= L + 1e-9; x += step) { const px = X(x, w); ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, 4); ctx.stroke(); ctx.fillText(x, Math.min(Math.max(px, 8), w - 12), 14); }
}
function timeText(d, t) {
  if (d.st === "relax") return `MC step ${Math.round(t / d.e.dt).toLocaleString()}`;
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
      const stepR = Math.max(1, Math.floor(rows / (w * 2)));
      for (let r = 0; r < rows; r += stepR) {
        const t = tcol >= 0 ? d.out[r * nc + tcol] : r / rows * t1, yy = h - 3 - (s.v(r) - s.lo) / (s.hi - s.lo) * (h - 6);
        r ? ctx.lineTo(t / t1 * w, yy) : ctx.moveTo(t / t1 * w, yy);
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
function stRows(d) {
  const an = d.an, e = d.e;
  let g = [...Array(d.nG).keys()];
  if (an && !S.stAll) g = g.filter(i => an[i].category !== "unheated");
  const tr = i => an ? an[i].track : 0;
  return g.sort((a, b) => tr(a) - tr(b) || e.grains.x[a] - e.grains.x[b]);
}
function drawSt() {
  const d = cur(), sec = $("st-panel");
  sec.hidden = !(d.g && d.nF > 1 && S.stage === "write");
  if (sec.hidden) return;
  const c = $("st"), rows = stRows(d), left = 44, top = 16, bot = 20;
  const rowH = Math.max(1.5, Math.min(8, 460 / rows.length)), H = top + bot + rowH * rows.length;
  const { ctx, w, h } = fit(c, H), pw = w - left - 4;
  const key = `${S.stage}|${S.stAll}|${S.colour}|${rows.length}|${document.documentElement.dataset.theme || ""}`;
  if (!_st || _st.key !== key) {
    const img = pixImage(d.nF, rows.length, i => { const f = i % d.nF, r = rows.length - 1 - Math.floor(i / d.nF), g = rows[r];
      if (S.colour === "agree") { const wv = wantOf(d, g); return wv ? rgbOf(LUT.agree, idxDiv(gmz(d, f, g) * wv)) : hex(css("--mid")); }
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
      ctx.fillStyle = css(want > 0 ? "--up" : "--down"); ctx.fillRect(tx(s0) + .5, 4, tx(s0 + p.nb) - tx(s0) - 1, 6);
    }
    ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.1; ctx.setLineDash([4, 3]); ctx.beginPath();   // laser over each grain
    rows.forEach((g, r) => { const s = d.model.passStep(d.e.grains.x[g], d.an ? d.an[g].track : 0); r ? ctx.lineTo(tx(s), ry(r)) : ctx.moveTo(tx(s), ry(r)); });
    ctx.stroke(); ctx.setLineDash([]);
    if (d.an) {   // decision moments
      ctx.fillStyle = css("--ink");
      rows.forEach((g, r) => { const a = d.an[g]; if (a.decision_step != null) ctx.fillRect(tx(a.decision_step) - 1, ry(r) - rowH / 2, 2, Math.max(rowH, 1.5)); });
    }
  }
  ctx.fillStyle = css("--muted"); ctx.font = `500 10.5px ${css("--font")}`; ctx.textAlign = "right";
  const every = Math.ceil(14 / rowH);
  rows.forEach((g, r) => { if (r % every === 0) ctx.fillText(fmt(d.e.grains.x[g], 0), left - 5, ry(r) + 3.5); });
  const rs = rows.indexOf(S.sel);
  if (rs >= 0) { ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5; ctx.strokeRect(left - .5, ry(rs) - rowH / 2 - 1, pw + 1, rowH + 2); }
  ctx.textAlign = "center";
  const tEnd = d.time[d.nF - 1], dtick = tEnd > 2e-9 ? 0.5e-9 : tEnd > 4e-10 ? 1e-10 : 2e-11;
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
const stageTitle = st => ({ relax: "Relaxation", write: "Write", read: "Read-back" })[st];
function stageFacts(st) {
  const d = D[st], e = d.e;
  if (st === "write" && e.analysis) {
    const an = e.analysis.grains, heated = an.filter(a => a.intended && a.category !== "unheated"), ok = heated.filter(a => a.correct).length;
    return [["melted grains ✓", `${ok} / ${heated.length} (${fmt(100 * ok / heated.length, 0)} %)`], ["freezing T (measured)", `${fmt(e.Tc, 0)} K`],
      ["freeze lag", `${fmt(e.freeze_lag_nm)} nm`], ["peak head field", `${fmt(e.hamr.Hpeak, 2)} T`], ["head speed", `${fmt(e.hamr.v, 0)} m/s`],
      ["bits", e.hamr.bits_per_track * e.hamr.n_tracks], ["bit length", `${fmt(e.hamr.BL)} nm`], ["laser T", `${fmt(e.hamr.Tmin, 0)}–${fmt(e.hamr.Tmax, 0)} K`]];
  }
  if (st === "relax") return [["integrator", e.input["sim:integrator"]], ["temperature", `${e.input["sim:temperature"] || "–"} K`],
    ["steps", d.nF > 1 ? Math.round(d.time[d.nF - 1] / e.dt).toLocaleString() : "–"], ["grains", d.nG]];
  return [["dipole solver", e.input["dipole:solver"] || "–"], e.dipole && ["cells", `${e.dipole.nx - 1} × ${e.dipole.ny - 1} × ${e.dipole.nz - 1}`],
    e.dipole && ["cell size", `${e.dipole.cell_nm} nm`]];
}
function legendHtml() {
  return `<div class="legend-mini"><span><i class="ring"></i>laser isotherms, bold ${fmt(cur().e.Tc, 0)} K</span><span><i class="box"></i>head-field box</span></div>`;
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
  const panels = [];
  if (m) panels.push({ name: "T (K)", lo: m.p.Tmin, hi: m.p.Tmax, v: f => m.T(gx, gy, d.steps[f]), col: css("--ink"), ref: d.e.Tc });
  panels.push({ name: "|m|", lo: 0, hi: 1, v: f => gml(d, f, g), col: css("--muted") });
  panels.push({ name: "m<tspan baseline-shift='sub' font-size='9'>z</tspan>", lo: -1, hi: 1, v: f => gmz(d, f, g), col: css("--up"), ref: 0 });
  if (m) panels.push({ name: "B<tspan baseline-shift='sub' font-size='9'>z</tspan> (T)", lo: -m.p.Hpeak, hi: m.p.Hpeak, v: f => { const hd = m.head(d.steps[f]); return m.inBox(gx, gy, hd) ? m.field(d.steps[f]) * m.p.field_uz : 0; }, col: css("--down"), ref: 0 });
  const L = 6, R = 6, top = 4, gap = 8, bot = 16, ph = (H - top - bot - gap * (panels.length - 1)) / panels.length;
  const tx = f => L + (f - f0) / Math.max(1, f1 - f0) * (w - L - R);
  panels.forEach((P, i) => {
    const y0 = top + i * (ph + gap), yv = v => y0 + ph - (v - P.lo) / (P.hi - P.lo) * ph;
    el("rect", { x: L, y: y0, width: w - L - R, height: ph, fill: "none", stroke: css("--rule") });
    if (P.ref != null) el("line", { x1: L, x2: w - R, y1: yv(P.ref), y2: yv(P.ref), stroke: css("--muted"), "stroke-dasharray": "3 3" });
    let dd = ""; const st = Math.max(1, Math.floor((f1 - f0) / 300));
    for (let f = f0; f <= f1; f += st) dd += (f === f0 ? "M" : "L") + tx(f).toFixed(1) + " " + yv(P.v(f)).toFixed(1);
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
  if (S.stage !== "write" || !e.analysis) { sec.hidden = true; return; }
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
  const k = 2 * Math.log((p.Tmax - p.Tmin) / (e.Tc - p.Tmin));
  let dd = "";
  for (let i = 0; i <= 120; i++) { const yy = -p.sigma_y * Math.sqrt(k) + i / 120 * 2 * p.sigma_y * Math.sqrt(k), r = k - yy * yy / p.sigma_y ** 2; if (r < 0 || Math.abs(yy) > ymax) continue; dd += (dd ? "L" : "M") + sx(p.sigma_x * Math.sqrt(r)) + " " + sy(yy); }
  el("path", { d: dd, fill: "none", stroke: css("--ink"), "stroke-width": 1.3, "stroke-dasharray": "5 3" });
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
  const facts = [["source", e.source], ["reduced", S.M.created]];
  if (e.grain_map) facts.push(["grain shapes", `from atom coordinates${e.grain_map.from_stage ? ` (${e.grain_map.from_stage} stage)` : ""}, ${e.grain_map.mixed_columns} mixed columns`]);
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
    $(id).disabled = v !== "grains" && !(d.snaps && (v === "atoms" || d.e.snapshots.volume));
  }
  const hasAn = !!(d.an || (D.write && D.write.an));
  for (const [id, v] of [["c-mz", "mz"], ["c-agree", "agree"], ["c-peak", "peak"], ["c-minm", "minm"]]) {
    $(id).setAttribute("aria-pressed", S.colour === v);
    $(id).disabled = (v === "peak" && !d.model) || (v === "minm" && !(d.g || d.snaps)) || (v === "agree" && !hasAn);
  }
  $("layer-ctl").hidden = S.view !== "layer" || !d.snaps;
  if (d.grid) { $("layer").max = d.grid.nz - 1; $("layer-out").textContent = `z ${fmt((S.layer + .5) * d.grid.pitch_z_nm, 2)} nm`; }
  $("z-ctl").hidden = $("w-ctl").hidden = S.stage !== "read" || !d.e.dipole;
  if (d.e.dipole) { $("zsel").max = d.e.dipole.nz - 2; $("zsel").value = S.z; $("z-out").textContent = zText(S.z); $("w-out").textContent = `${S.rw} nm`; }
  $("st-all").checked = S.stAll;
}
function bindControls() {
  [["v-grains", "grains"], ["v-atoms", "atoms"], ["v-layer", "layer"]].forEach(([id, v]) => $(id).onclick = () => { S.view = v; syncControls(); redraw(); });
  [["c-mz", "mz"], ["c-agree", "agree"], ["c-peak", "peak"], ["c-minm", "minm"]].forEach(([id, v]) => $(id).onclick = () => { S.colour = v; syncControls(); redraw(); });
  $("layer").oninput = ev => { S.layer = +ev.target.value; syncControls(); redraw(); };
  $("zsel").oninput = ev => { S.z = +ev.target.value; syncControls(); redraw(); };
  $("rw").oninput = ev => { S.rw = +ev.target.value; syncControls(); redraw(); };
  $("st-all").onchange = ev => { S.stAll = ev.target.checked; redraw(); };
  $("play").onclick = togglePlay;
  $("theme").onclick = () => {
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
    else if (ev.key === "Escape") { S.sel = -1; redraw(); }
  });
  $("save-png").hidden = !!EMB;
  $("save-png").onclick = savePng;
  window.addEventListener("resize", redraw);
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
    if (S.f >= d.nF - 1) togglePlay();
  }
  requestAnimationFrame(tick);
}
function prefetch(d) {
  if (!d.snaps || S.view === "grains") return;
  const s = nearestSnap(d, d.steps[Math.min(d.nF - 1, S.f + 60)]);
  if (s) snapArray(d, S.view === "atoms" ? "top" : "vol", s.index);
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
  const d = await loadStage(st);
  if (st !== "write" && S.M.stages.write && S.M.stages.write.analysis) await loadStage("write");   // intended bits for "vs intended"
  S.f = S.fa = st === "write" ? 0 : d.nF - 1;
  const e = d.e;
  S.ysec = e.hamr ? e.hamr.TW * .5 + e.hamr.pad : e.size_nm[1] / 2;
  if (st === "read" && S.M.stages.write && S.M.stages.write.hamr) { const p = S.M.stages.write.hamr; S.ysec = p.TW * .5 + p.pad; }
  if (e.dipole) S.z = Math.min(e.dipole.nz - 2, Math.round(e.size_nm[2] * .55 / e.dipole.cell_nm));
  if (S.view !== "grains" && !d.snaps) S.view = "grains";
  if (S.colour === "peak" && !d.model) S.colour = "mz";
  document.querySelectorAll("#chain button").forEach(b => b.setAttribute("aria-pressed", b.dataset.stage === st));
  try { history.replaceState(null, "", "#" + st); } catch (err) { /* sandboxed */ }
  buildLanes(); syncControls(); drawAnalysis(); drawDetails(); redraw();
}
function buildChain() {
  const nav = $("chain");
  nav.innerHTML = "";
  for (const st of ["relax", "write", "read"]) {
    const e = S.M.stages[st], b = document.createElement("button");
    b.dataset.stage = st; b.disabled = !e;
    let sub = "not in bundle";
    if (e) sub = st === "write" && e.hamr ? `${e.hamr.bits_per_track * e.hamr.n_tracks} bits, ${fmt(e.hamr.v, 0)} m/s`
      : st === "relax" ? `${e.input["sim:integrator"] || "relax"}, ${e.input["sim:temperature"] || ""} K` : e.dipole ? "dipole field" : "1 step";
    b.innerHTML = `<b>${stageTitle(st)}</b><span>${sub}</span>`;
    b.onclick = () => setStage(st);
    nav.appendChild(b);
  }
}

async function main() {
  try { const t = localStorage.getItem("hamrvis-theme"); if (t && !EMB) document.documentElement.dataset.theme = t; } catch (e) { /* storage unavailable */ }
  buildLuts(); bindControls();
  try {
    S.M = JSON.parse(new TextDecoder().decode(await getBuf("manifest.json")));
  } catch (err) {
    $("run-sub").textContent = "";
    document.querySelector(".work").innerHTML = `<div class="error"><h2>No bundle</h2>${kv([["error", err.message], ["start with", "python3 -m hamrvis view &lt;bundle&gt;"]])}</div>`;
    return;
  }
  const st = S.M.stages, any = st.write || st.relax || st.read;
  $("run-name").textContent = S.M.name;
  const nat = (st.write || st.relax || {}).n_atoms;
  $("run-sub").textContent = [nat ? `${nat.toLocaleString()} atoms` : null, `${any.grains.id.length} grains`,
    `${any.size_nm.map(v => fmt(v, 0)).join(" × ")} nm`].filter(Boolean).join(", ");
  document.title = `${S.M.name} HAMR run`;
  buildChain();
  const want = location.hash.slice(1);
  await setStage(st[want] ? want : st.write ? "write" : st.relax ? "relax" : "read");
  window.addEventListener("hashchange", () => { const h = location.hash.slice(1); if (st[h] && h !== S.stage) setStage(h); });
}
main();
