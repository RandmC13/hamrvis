# Design decisions

Why hamrvis is built the way it is.

Measurements below come from two real runs. The *large run* is a 200 x 60 x 20 nm granular film with 8,724,701 atoms
including a 10 nm non-magnetic sensor layer, 416 grains, 33 bits at 50 m/s and a 4.29 ns write, with 430 binary
file-per-node snapshots (about 90 GB). The *two-track run* has 2 tracks of 11 bits with field-free edge bits, text
snapshots and a film cut by a geometry file (1.6M atoms).

## Architecture

- **Reduce where the data is, view anywhere.** Moving 90 GB to a laptop is the bottleneck, and a browser cannot hold
  it either, so a Python reducer runs on the cluster next to the output and writes a bundle; only the bundle is
  copied. Alternatives rejected: copying raw data (90 GB per run), a per-run generated web page (every old page
  freezes at the version that built it, and a 16 MB page cannot hold atom data of a large run), a desktop GUI
  (dependencies, no sharing).
- **One general viewer + a local server** (`hamrvis view`, Python stdlib `http.server`) instead of a page per run:
  one codebase, every bundle opens in the same, latest viewer. Snapshot frames are separate files the viewer
  fetches on demand, so a multi-GB bundle costs the browser one frame at a time.
- **`hamrvis export`** inlines the viewer and a thinned bundle into one HTML file for sharing; grain frames are thinned
  with `--grain-step` to stay under a size limit (16 MB by default).
- **Dependencies**: Python 3.8+ and numpy only; the viewer is vanilla JS + canvas/SVG, no build step, one `app.js`, so
  export can inline it. Fonts from Google Fonts with system fallbacks (works offline with the fallback).

## Reduction

- Snapshot binary layout (`src/config`): coordinates `u64 n, i32 type[n], i32 category[n], f64 xyz[3n],
  i32 grain[n]`; spins `u64 n, f64 xyz[3n]`; text variants line-based. The `.meta` file lists the data files, so
  file-per-node, file-per-process and single-file outputs are handled alike.
- Only magnetic atoms (materials without `non-magnetic` in the .mat file) are binned: about 45 % of a snapshot of the
  large run is the sensor layer and fill.
- **Top-down columns at the atomic pitch** (2.5 A: 800 x 240 for the large run) and a **per-layer volume at 2x the
  pitch** (400 x 120 x 40 layers); int8 components (resolution 1/127). Measured on the large run's final snapshot:
  top 0.35 MB, volume 4.6 MB per snapshot gzipped, so about 150 MB + 2.0 GB for all 430 snapshots. Output fidelity is
  preferred to saving disk, so the volume is on by default (`--no-volume` drops it).
- The volume serves the one-layer view, the x-z cross-section at any y, and (later) any other slice.
- `grain.txt` per-grain block only (int16); the per-material block (large run: 5 materials x 416 grains x 4) is skipped.
- Parallel over snapshots: about 2 s per snapshot of the large run per core (see Any size for the worker model).
- Read stage: `dipole-field.txt` cells gridded by ID (see bundle-format.md). The per-atom
  `atomistic_dipole_field.txt` is not reduced yet (marked `shortcut:` in reduce.py).

## Accuracy

- **Grain shapes come from the atoms.** `grain_shapes.txt` polygons were tried first. Checked against the large run's
  atoms, only 85.6 % of magnetic columns lay inside their own grain's polygon (5.1 % in another grain's, 9.4 % in
  none): the polygons are the Voronoi cells before `voronoi-elliptical-rounding`, they miss grains that wrap across
  the periodic boundaries, and some are matched to the wrong grain. The reducer builds a top-down grain-ID map
  (one value per 2.5 A column, from each atom's grain ID in the coordinate files; large run: 0 columns hold two grains)
  and the viewer colours those pixels. Polygons are used only when a bundle has no coordinates, and are labelled
  approximate.
- **Grain centres and atom counts come from the atoms** (periodic circular means). `grain-coordinates.txt` counts
  non-magnetic atoms assigned to a grain: in the large run all `fill` atoms belong to grain 0 (48,483 atoms listed,
  13,127 magnetic), which moved its listed centre 15 nm. Its magnetisation in `grain.txt` is not affected.
- **Built-in checks**, stored in the manifest and shown under Run details: (1) every stored column value against
  the exact float64 mean (int8 rounding, <= 1/254); (2) per-grain mean spin of every snapshot against `grain.txt`.
  VAMPIRE writes snapshot k at step `k x rate + 1`, one step after the grain row of the same output, so for a
  finite-temperature run the comparison is one step apart (differences = one step of thermal motion; about 2e-2 for
  the 290-atom grains of the test fixture). At 0 K and at a shared step the means agree to 0 (tests/fixture-t0).
- **What reduction keeps and drops** (90 GB -> about 2.2 GB for the large run): all grain-level time series at full time
  resolution (int16, 3e-5); every snapshot (none dropped by default); every magnetic column at full atomic in-plane
  resolution (through-thickness mean); every atomic layer at 2 x 2 sites (0.5 nm); all as int8 per component
  (resolution 0.008). Dropped: the non-magnetic atoms, the per-atom detail inside a 2 x 2 x 1 cell, and precision
  beyond 1/127. Analysis numbers (outcomes, freezing temperature) are computed from the full-precision files.
- **Randomly ordered media.** With `initial-spin-direction = random`, a Monte Carlo relaxation orders each grain up or
  down at random (47 % up in the large run), and every grain starts the write with the sign it ended the relaxation
  with. A saturated start makes writing easy to see; for random media the "m_z vs intended" colour mode shows agreement
  with the bit each grain sits in, and greys out grains the laser never melts (large run: only 156 of 416 grains, a
  band of about 24 nm, melt in a 58 nm track).

## Generality

- The head model copies VAMPIRE's track logic: `max_tracks = floor((Ly - 1 A - 2 padding) / track_size)`,
  `bits_per_track = floor((Lx - 1 A) / bit_size)`, tracks reduced to `round(number_of_bits / bits_per_track)` only
  if that is fewer; sequences trimmed or padded with field-free 0 bits to fill the tracks (user-defined sequences
  read from `hamr:bit-sequence`). Checked against the two-track run's `B_applied` (1e-16 T) and head positions in its
  log.
- `grain.txt` leading columns are taken from the input's `grain:` keywords in order (the two-track run has
  `time-steps`, the large run `real-time`); unknown-width keywords before `grain:magnetisation` stop the reducer with a
  message naming the table to extend.
- Grids adapt (TOP_MAX, LAYER_MAX in reduce.py); `--layer-factor 1` gives per-atom layers.
- Coordinate files without grain IDs: Voronoi outlines, labelled approximate; snapshot checks against grain.txt
  need the IDs and are skipped.
- Truncated snapshots (copies caught mid-transfer) are skipped with a warning and listed in the manifest.
- Space-time rows shrink below 1 px for large runs (canvas height capped), track separators for several tracks;
  time ticks spaced by the run length.

## Viewer behaviour

- The theme toggle sits next to the run name, so that it does not push the stage buttons on narrow screens.
- Start / end buttons (Home / End keys) compare the first and last frames.
- Correct-grain counts follow the timeline (bits lane, inspector: "at t" and "at end"); the analysis section
  keeps the final numbers.
- *Bit cells* overlay: bit regions on every track, and the band hotter than the freezing temperature.
- Atom views print the snapshot time on the image when it is not the timeline's time. The grains view is grain.txt at
  every grain frame (1 ps in the large run); atom views are the nearest spin snapshot (10 ps in the large run), so with
  few snapshots the two differ by the time between them.

## Bit polarity

- A better measure of a write than single grains. A bit's polarity = mean m_z over the bit cell (x shifted by the
  freeze lag, one bit length) inside the band hotter than Tc (|y - track centre| <= sigma_y sqrt(2 ln((Tmax - Tmin)/(Tc
  - Tmin))); the whole track if there is no Tc), weighted by the grain-map pixels each grain has in that region
  (pixels ~ atom columns), using the grain magnetisation of the current frame. Without a grain map: grains whose centre
  lies in the region, weighted by atom count. A bit is correct when the sign matches the intended field direction.
- Large run: 29 of 31 bits correct at the end (the first bit cell lies off the film once shifted by the lag, and
  the 33rd bit was never written), against 116 of 142 melted grains inside the bit cells; a cell holds parts of about 9 grains
  (8-11) of either sign, so its mean |m_z| is well below 1 (median 0.45, 0.15-0.56 for 80 % of the cells).

## Runs instead of fixed stages

- hamrvis must not assume relax -> write -> read, since workflows differ. It reads a directory, works out which runs
  are in it and what each is, and builds a tab for each. There is no code path per workflow: a registry of run types
  (`runs.py` TYPES) that can grow, with views otherwise following the data a run has.
- Types for now: hamr-simulation, Monte Carlo, dipole read; others get a generic tab.
- Order: runs chained with `load-checkpoint = restart` all start at simulated time 0, so simulated time cannot order
  them. The order of the log's first time stamp is used, then the folder name; `--order` overrides.
- Jobs of one run continued from checkpoints are joined: same input apart from the checkpoint keys, resumed at a step
  the previous job reached; rows and snapshots of the earlier job after the resume step are dropped. Checked: a write
  run as two jobs (first killed mid-write) reduces byte-identical to the same write in one job. A resumed job logs its
  bits again, so bits per track count distinct "Written bit N" lines.
- Bundle format 2: `runs` list, `<run id>/` folders, `output` as float64. Version 1 bundles open unchanged except
  for the tab bar (only tabs for runs present); checked by pixel comparison of screenshots against the old viewer.

## Physics choices

- The head field is rebuilt from the input with VAMPIRE's own trapezoid (`docs/analysis.md`), not read from
  `output`, as a run need not write `B_applied`. Note for reading inputs: with only `sim:applied-field-strength` set,
  Hmax and Hmin keep their defaults (+1 / -1 T) and the ramp climbs Hmax - Hmin, so the peak field is 2 T, not 1 T. This
  is a property of such inputs, not of hamrvis.
- **Bits lost to the lead-in.** In earlier VAMPIRE versions `hamr_continuous.cpp` starts its bit loop after the lead-in
  (half a bit + NPS, `ne` steps) but ends it at `bits_per_track * nb` steps into the track, so it writes
  `ceil((bits_per_track * nb - ne) / nb)` bits per track and leaves the field off over the rest. The large run (6 nm
  bits, 5 nm NPS) reports 33 bits and writes 32; its log has 32 "Written bit" lines. The head's track is
  `bits written so far / bits_per_track` (integer division), so with bits lost later tracks are written over earlier
  ones. hamrvis models both (`Hamr.written`, `Hamr.track_y`); cells never written have intended polarity 0, shown
  as "no field". In fixed versions every bit is written, each track at its own y. hamrvis takes the bits per track from
  the run's log ("Written bit" lines up to the first "end of track" line), so runs from both versions are modelled, and
  falls back to the old formula when the log has no whole track.
- Correctness uses the signed z field (`bit polarity x unit-vector z`), since a run may write with (0, 0, -1).
- The freezing temperature is inferred from the data (median centre temperature at the decision moments), 588 K for the
  large run, instead of requiring a Curie temperature on the command line.

## Visual design

- **Aligned lanes**: every view of the write shares the down-track axis (film x), stacked like oscilloscope
  channels, with one crosshair through all lanes; inspector column on the right; transport bar under the lanes.
  Down-track position is the coordinate everything in recording is about (bits, transitions, laser, read-back).
- Colour is reserved for data. m_z: blue #2a78d6 (down), grey #f0efec, red #e34948 (up), saturating at |m_z| = 0.85; dark
  theme #3987e5 / #383835 / #e66767. Fields use a separate green-purple map. Correct = green, wrong = ochre (red and
  blue already mean up and down). Laser isotherms in ink; the field box takes the colour of the field direction.
- Lines over the grain imagery (isotherms, field box and its label) get a halo in the lane colour so they read over red
  and blue; the freezing-T isotherm is solid and bold, the others dashed, dash-dot and dotted with falling T, and the
  legend keys each one by temperature.
- The bits lane shows every track at once (one labelled row per track, stacked as on the film, the track under the head
  marked "writing"), so no track switch is needed. A bit's counts are greyed until one of its grains has frozen (its
  decision step): before that a grain matches its bit only because of the starting magnetisation.
- The transport bar sits directly under the film (the main view), with slow-motion speeds (0.1x-0.5x). Space-time has one
  row per melted grain sorted by x, time left to right, laser pass dashed, decision ticks. Labelled data instead of
  prose everywhere.
- Chrome: cool aluminium / slate neutrals (the platter), one typeface (Instrument Sans, condensed bold for titles,
  tabular figures for numbers), sentence-case labels.

## Any size (target about 100M atoms, 10^4 grains, 10^5 frames, films to about 2 um)

Measured on a generated run (`tests/make_large.py`, 2000 x 60 nm, 99.84M atoms in one file, 9,809 grains, 100,199 frames, 12 GB grain.txt, 5 snapshots) on a machine with 7 GB and 12 CPUs: reduce 365 s, peak private memory 4.35 GB, of which 5 workers about 1 GB each; a 7.9 GB grain series opens in the viewer in 2 s with a 139 MB JS heap. On a 20M-atom film the previous reducer used 4.08 GB with 2 workers where this one uses 1.77 GB, and with its default worker count it lost a worker to the out-of-memory killer and hung. Generated grains are uniform, so gzip does far better on them than on real runs: compare bundle sizes with `inspect`'s estimate, not with this run.

- **No format change for size.** `grains.i16` is contiguous by frame, so byte ranges are enough to read parts of it; every
  old bundle opens as before, and default bundles of existing runs stay byte-identical. Only a series above 512 MB (64M
  frame-grains) gains one new file, the overview, and is read in blocks.
- **Whole-run views need extremes, not samples.** Space-time, the time line and lowest |m| read the overview: min and max per
  time bin, so a short reversal in a 10^5-frame run still shows. Sampled frames would miss it. Bins are 1024: 82 MB at 10^4
  grains, loaded whole. In each pixel the extreme furthest from the previous pixel is drawn.
- **A grain's life story reads only that grain.** 300 frames of one grain are 300 8-byte ranges, instead of a second copy of
  the series in grain order (8 GB more for the same numbers).
- **Grain analysis in one pass.** A grain's decision depends on the sign of its final m_z, known only at the last frame, so
  the accumulator follows both signs. This replaced a [frames, grains] array and gives the same numbers (byte-identical
  manifests; a test checks it against the whole-array formulas).
- **Workers use spawn on every platform**, not fork: one code path everywhere (fork is deprecated on macOS, absent on
  Windows). Per-atom indices (13 bytes per magnetic atom) are memory-mapped files, so workers share one copy.
- **Automatic resume, automatic worker count**: no flag to remember. The same command after a time limit carries on; the
  number of workers follows the job's CPUs and memory.
- **Viewer: one x window, true aspect ratio, y scrolls.** Overlays keep their film-coordinate drawing under a canvas
  translation, so they needed no rewrite. A plain wheel on a lane is left to the page unless the lane itself can scroll: a
  zoom on wheel would stop the page from scrolling whenever the pointer is over a lane, so zoom is pinch or ctrl + wheel.
  Plain drag stays the section line, pan is shift + drag.
- Images are built only for the visible part, at no more than twice the screen resolution by taking every n-th pixel: as
  before, an image shrunk on screen shows chosen pixels, not averages. Zoom in to see each column.
- Limits left: the dipole field map is loaded whole (29 MB at 2000 x 60 x 20 cells, about 1 GB at 2 x 2 um x 20 nm); a window
  snapshot run keeps one full state of the spins in memory (24 bytes per atom) since each snapshot updates the last.

## Testing

- `tests/fixture-run`: a 2,624-atom granular HAMR write made with VAMPIRE (regenerate: `tests/make_fixture.sh`).
- `tests/fixture-leadin`: the same film with a lead-in longer than a bit and two tracks (output and log only):
  VAMPIRE reports 6 bits and writes 4, the second track on the first (made before the fix; kept as recorded).
- `tests/fixture-leadin-fixed`: the same input run by the fixed VAMPIRE: 6 bits, two tracks.
- `tests/scenarios/`: folders of runs chained by checkpoints on the same film, made with VAMPIRE by
  `tests/make_scenarios.py` (about 40 s): MC relax / write / read; write only; one write over two jobs (the first killed
  mid-write); write / read / rewrite / read; an LLG anneal (generic tab) before a write. Machine-specific paths and host
  names in their logs are replaced with neutral ones.
- `tests/make_large.py`: a generated run of any size with a known answer per grain; `Scale` in `tests/test_hamrvis.py` uses a small one.
- `tests/test_hamrvis.py`: field model vs VAMPIRE's B_applied, default 2 T peak, column averages vs a direct numpy
  computation, analysis output, export, and the viewer's JavaScript model vs the Python one (node).
- `tests/test_viewer_smoke.py`: every scenario in a real browser through Playwright.

## Ideas not done yet

- Relaxation views: grain |m| and m_z histograms over MC steps; time to equilibrate (|m| reaching its plateau) per grain
  and overall; fraction of grains up against MC step; grain size against |m| (small grains are superparamagnetic);
  comparison of the relaxed state with the write's first frame (a checkpoint hand-over check: 416 of 416 grains of the
  large run keep their sign).
- Run comparison, for example a sweep of grain rounding.
- Movie export (canvas MediaRecorder to WebM, or ffmpeg from the reducer) and SVG export.
- Per-atom read-back field; peak temperature measured from snapshots rather than the laser model.
- Lossless workflow: a bundle that keeps all the simulation data, so the raw VAMPIRE files could be deleted. Today every
  bundle is lossy, even with `--layer-factor 1`: spin components are int8 (1/127), pixels hold means rather than atoms
  (no atom IDs or order), non-magnetic atoms are dropped, the top map coarsens above 4M pixels, counts cap at 255. Sketch:
  per-atom float32 spins per snapshot plus the coordinates, materials and grain IDs once (two-track run: about 18.6 MB per
  snapshot, about 2.7 GB vs 6 GB of text); the viewer keeps using the reduced grids.
