# How hamrvis works

Every number and picture the viewer shows, traced back to the VAMPIRE file it comes from, with the formula hamrvis
applies, a worked example from a real run, the test that checks it, and its known limits. How to use hamrvis is in
[user-guide.md](user-guide.md); the files of a bundle are in [bundle-format.md](bundle-format.md).

Each value is marked by where it comes from, as the viewer marks it on hovering a label:

- **recorded**: written by VAMPIRE, and copied, averaged or rounded by hamrvis;
- **reconstructed**: rebuilt by hamrvis from the run's input, as VAMPIRE computes it (the head model);
- **interpreted**: a definition hamrvis chose (decision, outcome, freezing temperature, freeze lag, bit correctness).

Paths such as `hamrvis/reduce.py` are hamrvis's; paths starting with `src/` are VAMPIRE's.

**Worked examples.** *chain* is the small test scenario in this repository, `tests/scenarios/chain` (12 × 6 × 1.2 nm,
2,624 atoms, 9 grains: a Monte Carlo relax, a write and a read-back made with VAMPIRE by `tests/make_scenarios.py`);
every chain number below can be reproduced with `hamrvis reduce tests/scenarios/chain`. *The large run* is a
200 × 60 × 20 nm granular FePt film of 8,724,701 atoms and 416 grains (Monte Carlo relax, a write of one 33-bit
track at 50 m/s, a dipole read), reduced with hamrvis 1.0. *The scale writes* are two later writes of the same system
from a new relax, one with `hamr:active-window` snapshots and one with full snapshots, 430 snapshots each. Their
bundles are not in the repository; their numbers were computed from the bundles while this was written.

Contents: [1 Overview](#1-from-vampire-files-to-the-screen) · [2 Runs and the input](#2-runs-tabs-and-the-input-file) ·
[3 Grain shape](#3-grain-shape-and-structure) · [4 Snapshots](#4-spin-snapshots) ·
[5 Grain series](#5-the-grain-series-recorded) · [6 Head model](#6-the-head-model-reconstructed) · [7 Write analysis](#7-write-analysis-interpreted) ·
[8 Read-back](#8-read-back-runs-recorded) · [9 Every readout](#9-every-readout) ·
[10 Tests, tolerances, limits](#10-tests-tolerances-and-limits)

## 1. From VAMPIRE files to the screen

hamrvis has two halves. `hamrvis reduce` runs next to the data: it reads the run folders, averages the atoms onto
grids, parses the tables, rebuilds the head, analyses the write and writes a *bundle*, a folder of small binary
arrays plus `manifest.json`. `hamrvis view` serves the bundle to a page (`hamrvis/viewer/app.js`) that loads arrays on
demand and computes the rest of what it shows. Only the bundle is copied home.

```
VAMPIRE run folders            hamrvis reduce (Python)                bundle            viewer (browser)
input, *.mat, log         ->   runs.find: kinds, order, joined jobs   manifest.json     head model (Model),
output, grain.txt              Input: keys and units as VAMPIRE       output.f64        scores at the lag in use,
grain-coordinates.txt          reduce_tables: output, grain.txt       grains.i16        bit cells and polarity,
atoms-coords.meta/.data        GrainAccumulator (one pass)            grain_map.i16     images, read-back waveform
spins-NNNNNNNN.meta/.data      coordinate pass: grid, grain map,      top/*.i8.gz
dipole-field.txt                 centres, per-atom bin indices        vol/*.i8.gz
                               snapshot workers: column and layer     top_count, vol_count,
                                 means, lowest |m|, checks            min_m
                               Hamr: head model and its check         dipole.f32,
                               write_analysis: Tc, lag, outcomes      dipole_atoms.u16
                               reduce_dipole: cells by ID
```

The bundle holds what VAMPIRE recorded, averaged and rounded, plus the reducer's analysis; the viewer recomputes the
head and the scores itself, so that the lag can be changed on screen. `hamrvis verify` re-reads the raw files with a
separate reader and checks the bundle against them ([section 10](#10-tests-tolerances-and-limits)).

## 2. Runs, tabs and the input file

### Finding runs and giving each a tab

A path given on the command line counts as a run if it holds a VAMPIRE input file. In a folder that is scanned, a
subfolder counts only if it holds an input file and an output table, so copies of inputs or source trees next to the
runs are skipped with a message (`runs.candidates`). Each run gets the kind of the first entry in the registry
(`runs.TYPES`) that matches its input:

| Kind | Tab title | Matches when the input has |
|---|---|---|
| `hamr` | Write | `sim:program = hamr-simulation` |
| `dipole-read` | Read-back | `dipole:output-dipole-field` or `dipole:output-atomistic-dipole-field` |
| `monte-carlo` | Monte Carlo | `sim:integrator` containing `monte-carlo` |
| `generic` | the program name | anything else |

**Order.** Runs chained with `sim:load-checkpoint = restart` all start at simulated time 0, so simulated time cannot
order them. hamrvis orders runs by the first time stamp of each `log` (`dd-mm-yyyy [hh:mm:ss]`, `vampire.log_start`),
then by folder name; `--order` overrides. **Reference write**: every run that is not a write takes the nearest write
before it as its reference, the bits it is compared with (a read-back's intended transitions, "m_z vs intended", bit
polarity). Tab ID = folder name with unsafe characters replaced, made unique with `-2`, `-3`.

*Large run*: three tabs in log order, relax (Monte Carlo), write (hamr), read (dipole-read, reference = write). The
relax wrote no coordinates, so its grain map is shared from the write (section 3).

### Jobs continued from checkpoints

Consecutive folders become one run when the later one has `sim:load-checkpoint = continue`, the same input apart from
the checkpoint keys and the keys that only say where it runs or how its files are laid out (`config:output-mode`,
`config:output-nodes`, `gpu:*`), and resumed (log line `Checkpoint file loaded at sim::time N`) at a step the earlier
job's output reached (`runs._continues`). Rows, grain frames and snapshots of the earlier job after the resume step
are dropped; the later job wins where they overlap (`reduce._join_rows`, `run_snapshots`, `_stream_grains`). A
resumed job replays its loop and logs its bits again, so bits are counted as distinct "Written bit N" numbers.

A job continued on another number of processes or GPUs, or with another output mode, lists its atoms in another order.
Each job that wrote its own coordinate files therefore gets its own atom layout: its coordinates are binned again, and
the reduce stops (`its atoms differ from those of the job it continues`) unless they fill the columns exactly as the
first job's did. Each snapshot is read with the layout of the job that wrote it; a job without coordinate files uses
the layout before it. A window snapshot cannot start a job with a new layout, since it updates the previous state
atom by atom.

Checked by `Bundles.test_joined_jobs_equal_one_job` (a write run as two jobs, the first killed mid-write, reduces
byte-identical to the same write in one job), the `mpi-continue` scenario (a write continued from 2 MPI ranks on 1,
every snapshot of both jobs checked against the raw files), `RunDiscovery.test_resumed_log_counts_bits_once` and
`test_order_follows_the_logs_not_the_names`.

### Reading the input as VAMPIRE does

hamrvis copies VAMPIRE's input rules (`vampire.read_keys`, `Input.num`):

- Lines are `key = value !unit`; `#` starts a comment; all whitespace inside a key is removed; the last line setting
  a key wins.
- Units come from VAMPIRE's own table (`src/utility/units.cpp`, copied into `_UNITS`): a value is converted to
  VAMPIRE's internal unit and then to hamrvis's (nm, s, T, m/s). An unknown unit, or one of the wrong kind for the key,
  stops with a message, as VAMPIRE would. A value without a unit is in VAMPIRE's internal unit (a bare
  `hamr:head-speed = 3e11` is 3e11 Å/s = 30 m/s).
- Integer keys ignore any unit. Integer lists (`hamr:bit-sequence`) are read as `vin::integers_from_string` reads
  them: comma separated, each field read with `stream >> int`, so `{1,0,-1}` reads 0, 0, −1.
- `dimensions:system-size` sets all three axes, `-x`/`-y`/`-z` one; the last line wins; default 10 nm.
- Defaults are VAMPIRE's where the input is silent (time step 1 fs, Tmax 300 K, Hmax +1 T, Hmin −1 T, head speed
  30 m/s, FWHM and field box 20 nm, ramps 1 ps).
- File names: `input`, else the name in the log line `Opening main input file "..."`; `output`, else the one file whose
  second line is VAMPIRE's header; `--input-file` / `--output-file` override.
- Magnetic materials: every `material[n]` of the `.mat` file without `:non-magnetic`. Without the `.mat` file every
  material counts as magnetic and the run gets a note.

*Large run, write*: `sim:time-step = 1.0 !fs` gives 1e-15 s; `hamr:head-speed = 50 !m/s` gives 50 m/s;
`hamr:bit-size = 6 !nm`; `cells:macro-cell-size = 1.0 !nm`. Materials 1–3 (FePt) are magnetic, 4 (fill) and 5
(sensor) are `non-magnetic = keep`, so the magnetic types are 0, 1, 2 (0-based, as in the coordinate files). Of the
8,724,701 atoms, 4,856,385 are magnetic; the sensor layer alone is 3,840,000.

Checked by `Parsing.test_units_as_vampire_converts_them`, `test_bad_units_say_what_is_wrong`,
`test_keys_as_vampire_reads_them`, `test_custom_names_and_units_give_the_same_bundle`,
`Bundles.test_defaults_are_vampires`, `test_system_size_sets_every_axis`, and verify's "head parameters" check
(input values in any units against an independent reading).

## 3. Grain shape and structure

Grain shapes come from the atoms. VAMPIRE's `grain_shapes.txt` polygons were tried first and rejected: on a large
run only 85.6 % of magnetic columns lay inside their own grain's polygon. The polygons are the Voronoi cells before
`voronoi-elliptical-rounding`, they miss grains that wrap across the periodic boundary, and some are matched to the
wrong grain ([decisions.md](decisions.md), Accuracy).

### Grain IDs (recorded)

The binary coordinate file of each data file is `u64 n, i32 type[n], i32 category[n], f64 xyz[3n], i32 grain[n]`
(`src/config`; xyz in Å). The grain array is present in newer files: hamrvis takes it when the file is at least
`8 + 36 n` bytes long (`vampire.coord_chunks`). Text coordinate files have the grain in column 6. The grain list
itself (IDs, listed centres, atom counts) comes from `grain-coordinates.txt`, which VAMPIRE writes only for granular
systems; without it a run has atom views only.

### The column grid

The grid is set by the atoms themselves (`reduce.Grid`):

```
a   = smallest spacing between distinct magnetic-atom x values      (coordinates within 1e-3 Å count as equal)
a_z = smallest spacing between distinct z values,   z0 = lowest z
nx0 = round(Lx / a),  ny0 = round(Ly / a),  nz = round((z_max − z0) / a_z) + 1
tf  = max(1, ceil(√(nx0 ny0 / 4e6)))            top map coarsened only above 4M pixels (--top-max)
lf  = max(2, ceil(√(nx0 ny0 nz / 2.5e6)))       layer cells: lf × lf sites per layer (--layer-max, --layer-factor)
column of an atom:  ix = rint(x / a) mod nx0,  iy = rint(y / a) mod ny0,  pixel (iy // tf, ix // tf)
layer cell:         (iz, iy // lf, ix // lf),  iz = clip(rint((z − z0) / a_z), 0, nz − 1)
```

The modulo wraps periodic coordinates. Only magnetic atoms are binned. Layer cells hold 2 × 2 sites because one site
per layer can be empty in non-cubic lattices, while 2 × 2 never is inside a grain.

**Pixel geometry.** Pixel i spans `o + i p` to `o + (i + 1) p` with `o = −a/2` and `p = f a` (viewer: `aGeo`,
`pixOf(g, v) = floor((v − o) / p)`), so a pixel is centred on its atoms. Its centre, used wherever a pixel must be
placed against a bit boundary, is the mean of its atom columns, `a (i f + (f − 1)/2)` (`pixC`).

*Large run*: the FePt lattice has x values every 2.5 Å, so a = 0.25 nm, nx0 = 800, ny0 = 240, tf = 1 (192,000
pixels), 40 atomic layers 0.25 nm apart from z0 = 0, layer grid 400 × 120 × 40 cells of 0.5 nm. Each column holds
25 to 40 magnetic atoms; 61,770 of the 192,000 columns hold none (grain boundaries).
*chain*: a = 3 Å, grid 40 × 20, layer grid 20 × 10 × 4. Atom 137 sits at (84, 27, 6) Å: rint(84/3) = 28,
rint(27/3) = 9, so it lands in column (ix 28, iy 9) with 3 other atoms (4 layers).

### The grain map (recorded)

While reading the coordinates the reducer counts atoms per (column, grain) pair. Each column takes the grain with the
most magnetic atoms in it; a column with no magnetic atom gets −1. `mixed_columns` counts columns holding atoms of
more than one grain, so a reader can see whether the majority rule ever had to choose. Grain IDs not in
`grain-coordinates.txt` fall into a catch-all bin and are mapped to −1.

*Large run*: 0 mixed columns, so the map is exact. 415 of the 416 grains appear on it; the missing one, grain 415,
has no magnetic atoms (section 7). *chain*: 0 mixed columns, 144 of 800 columns empty.

Checked by `Reduce.test_grain_map_matches_atoms`, verify's "grain map" and "mixed columns" checks, and the viewer
hover test (section 10).

### Grain centres and magnetic atom counts (recorded)

The listed values in `grain-coordinates.txt` include non-magnetic atoms VAMPIRE assigns to a grain, so the reducer
replaces them with values from the magnetic atoms. Centres are periodic (circular) means, so a grain wrapped across
the boundary gets its true centre:

```
x_c = (atan2(Σ sin θ, Σ cos θ) mod 2π) · Lx / 2π,     θ = 2π x / Lx     (likewise y)
magnetic_atoms = number of magnetic atoms with that grain ID
```

A grain with no magnetic atom keeps its listed centre. Values are rounded to 1 pm (3 decimals in nm).

*Large run, grain 0*: listed (110.820, 28.656) nm, 48,483 atoms; from the atoms (126.203, 26.908) nm, 13,127
magnetic atoms. The other 35,356 are fill and sensor atoms VAMPIRE assigned to grain 0, which pulled the listed centre
15.5 nm away. Its magnetisation in `grain.txt` is VAMPIRE's own and is not affected. Grain 415: 1,235,400 atoms, none
magnetic, so it keeps its listed centre.

### Without coordinates: shared maps, then Voronoi outlines

A run that wrote no coordinate files (often a relax or a read) takes the grain map, centres and magnetic atom counts
of the nearest run with the same system size and the same grain IDs, preferring earlier runs on a tie
(`share_grain_map`); Run details then says the outlines are drawn from the atom positions of that run. Only if no
such run exists does the viewer fall back to VAMPIRE's Voronoi polygons: each grain takes the polygon that contains its
centre, or the polygon with the nearest centre (`vampire.grain_outlines`; `grain_shapes.txt` lists every cell made
before a geometry cut, so its order is not the grain order). The film lane then says "approximate Voronoi outlines".
*Large run*: relax and read share the write's map.

### How the viewer draws grains

Grains view: each pixel of the grain map is painted in its grain's colour at the current grain frame; −1 pixels are
the void colour. Image rows count from the bottom of the film, so the image is flipped for display. Large films:
only the visible window is built, at no more than twice the screen resolution, by taking every n-th pixel, so a
zoomed-out view shows chosen pixels, not averages. Outlines of the hovered and selected grain are its boundary pixels
(pixels of the grain with a 4-neighbour of another value, with periodic wrap), so they lie on grain-map boundaries by
construction. Hover converts the pointer to a pixel with the same geometry and reads the grain map; with Voronoi
outlines it uses a point-in-polygon test.

## 4. Spin snapshots

### Column means and layer means (recorded)

For each snapshot a worker reads the spin files in the coordinate files' order, keeps the magnetic atoms, and sums
each spin component into its column, layer cell and grain, using per-atom bin indices written once by the coordinate
pass and memory-mapped by every worker. The stored value is the mean direction of the atoms in a pixel:

```
top[c, iy, ix]      = Σ s_c / N_column           c = x, y, z;   N_column = top_count
vol[c, iz, iy, ix]  = Σ s_c / N_cell             N_cell = vol_count
min_m[iy, ix]       = min over the bundle's snapshots of |Σ s / N_column|      (stored × 255, uint8)
```

The column value is a through-thickness average, so a column whose layers point different ways shows a small mean.
Atom columns draws top[z]; One layer draws vol[z, layer]; the x–z section draws vol[z] at the layer-grid row
containing the section line. *Smooth* averages these values over 2 × 2 to 8 × 8 neighbouring pixels on screen,
ignoring empty pixels; the stored values and the tooltip are unchanged.

### Quantisation

| `--precision` | Stored as | Value | Empty pixel | Largest rounding error |
|---|---|---|---|---|
| 8 (default) | int8 | rint(127 · mean), clipped to ±127 | −128 | 1/254 = 0.003937 |
| 16 | int16 | rint(32767 · mean) | −32768 | 1/65534 |
| 32 | float32 | mean | NaN | float rounding |

Every snapshot is checked as it is written: the largest difference between a decoded stored column and the exact
float64 mean is kept (`checks.column_rounding_max`) and shown under Run details.

*chain, write, snapshot 5*: column (28, 9) holds 4 atoms with exact mean m_z = −0.444764. 127 × −0.444764 = −56.49,
so the file stores −56, which decodes as −56/127 = −0.440945, a difference of 0.0038, inside 1/254. The run's
largest rounding is 0.003937, the bound itself, as in all three large-run writes.

The viewer's m_z colour map runs blue–grey–red and saturates at |v| = 0.85 (256 entries); values beyond ±0.85 all
draw the end colour. Exact values are in the tooltip.

### Window snapshots (`hamr:active-window`)

With `hamr:active-window`, VAMPIRE writes snapshots holding only the atoms of grains integrated since the previous
snapshot: `u64 n, u64 r`, r × (first, count) ranges, then the spins, marked `Spin data: window` in the `.meta` file.
The reducer reads them in order and changes the column, layer and grain sums by the binned difference of the spins in
the window (new minus old), so only atoms near the head are touched (`_window_series`). The first snapshot must be
complete. Memory: one full spin state (24 bytes per atom). Reading in order makes them slower to reduce: 351 s for
the window scale write against 45 s for its full-snapshot twin.

Checked by `WindowSnapshots.test_series_rebuilds_every_snapshot`, `test_reduce_same_as_full_snapshots`, and verify,
which rebuilds every window snapshot from the complete one before it with its own reader (the window scale write:
430 snapshots, every column and layer within 0.00394).

### Checks against grain.txt, and why they are one step apart

For every snapshot the reducer also forms the per-grain mean spin (Σ s / magnetic atoms) and compares it with
`grain.txt`'s direction × |m| for the same grain. VAMPIRE writes the two on different counters: `grain.txt` when
`sim::time % grain:output-rate == 0` (`src/vio/grain.cpp`), snapshots when
`sim::output_rate_counter % config:atoms-output-rate == 0` (`src/config/config.cpp`), which leads by one step. So a
snapshot is usually written one step after the grain row of the same output. The reducer compares at the same step
when a grain row exists there, else one step earlier, and stores both results separately. At finite temperature the
one-step comparison differs by one step of thermal motion; at 0 K, or at the same step, the two agree to rounding.

| Run | Same step | One step earlier |
|---|---|---|
| chain write (11 snapshots, T > 0) | 0 | 9, median 0.019, max 0.030 (290-atom grains) |
| chain read | 2, max 6.0e-7 (float32 arithmetic) | 0 |
| scale write, window | 329, max 0.0049 | 99, median 0.0028, max 0.0048 |
| scale write, full | 329, max 0.0052 | 99, median 0.0030, max 0.0048 |
| large-run write | none: its one snapshot is at the end of the cool-down, after the last grain row | |

Checked by `Reduce.test_snapshots_match_grain_txt`, `test_snapshot_reading_exact_at_same_step` (`tests/fixture-t0`:
0 K, same step, difference 0) and `test_column_average_matches_spins`.

### Which snapshot the atom views show

Atom views show the snapshot nearest to the current grain frame in steps, before or after it. When it is not at the
current time, the image says `snapshot <t>` and the lane note gives the gap (Δ ps). Snapshots are fetched on demand,
the last 80 kept and the next two prefetched during playback; while one loads, the previous image and its tooltip
values stay up. The grains view, by contrast, is `grain.txt` at every grain frame.

## 5. The grain series (recorded)

`grain.txt`'s columns follow the input's `grain:` keywords in their order (`src/vio/match_grain_list.cpp`): single
columns first (`time-steps` → step, `real-time` → time, `temperature`, field, ...), then blocks per grain, of which
`grain:magnetisation` gives 4 per grain: m_x, m_y, m_z of the direction and |m|. A keyword of unknown width before
`grain:magnetisation` leaves the grain views out with a note naming the table to extend (`vampire.grain_columns`). A
per-material block that may follow is skipped.

```
frame time = the time column; else step × dt; else (row + 1) × grain:output-rate × dt   (with a note)
stored     = rint(clip(v, −1, 1) × 32767) as int16, from grain.txt's float64 value      (error ≤ 1/65534)
m_z shown  = (m_z / 32767) × (|m| / 32767)      the z-magnetisation of the grain, not its direction
|m| shown  = |m| / 32767
```

`grain.txt` is read one row at a time (it can be tens of GB). A short last row from a job stopped mid-write is left
out with a warning; a damaged row anywhere else leaves the grain views out with a note. `--grain-every N` keeps every
N-th frame plus the last; `--grain-precision 32` stores float32. The write analysis uses the kept frames in float32.

**Long series.** Above 512 MB the series also gets an overview: 1,024 time bins, each holding the lowest and highest
m_z and |m| of each grain. Whole-run views (space–time, time bar, lowest |m|) read the overview, so a short reversal
still shows; the current frame and a grain's life story are read by HTTP byte ranges, in blocks of frames on demand.

*Large run, write*: `grain:real-time`, `grain:temperature`, `grain:magnetisation`, `grain:material-magnetisation`,
every 1,000 steps: 4,289 frames 1 ps apart, 416 grains, 4,289 × 416 × 4 × 2 bytes = 14.3 MB, so no overview. The raw
row has 9,986 columns (2 + 4 × 416 + 4 × 5 × 416); the material block is skipped.

Checked by verify's "grain series" check, `Scale.test_grain_series_and_overview_hold_the_known_answer`,
`Parsing.test_half_written_last_row`, `TrackRules.test_grain_columns_follow_the_input`.

## 6. The head model (reconstructed)

VAMPIRE does not write the laser position, laser temperature or head field per step to any file a run need have, so
hamrvis rebuilds them from the input exactly as `src/hamr/hamr_continuous.cpp`, `field_profiles.cpp`,
`temp_profiles.cpp` and `initialize.cpp` compute them (Python `vampire.Hamr`; the viewer's port is `Model` in
`app.js`; `hamrvis/verify/reference.py` has a third, independent `Head`). It is then checked against what VAMPIRE did
record.

### Steps, bits and tracks

```
step_nm        = v · dt                                   head travel per step
nb             = round(BL / step_nm)                      steps per bit
ne             = round(BL/2 / step_nm) + round(NPS / step_nm)    lead-in = run-out (field off)
bits_per_track = floor((Lx − 0.1 nm) / BL)
max_tracks     = floor((Ly − 0.1 nm − 2 pad) / TW);  tracks = round(n_bits / bits_per_track) if fewer, else max_tracks
nt             = bits_per_track · nb + 2 ne              steps per track
laser on       for 1 ≤ k ≤ tracks · nt,  k = step − equilibration steps;  then 10 ps cool-down at Tmin
track t        = (k − 1) // nt,   step in track j = (k − 1) % nt + 1
laser x        = −BL/2 − NPS + (j + 1) · step_nm
laser y        = TW · (floor(t · w / bits_per_track) + 0.5) + pad        w = bits written per track
sequence       single-tone: −1, +1, −1, ...;  user-defined: hamr:bit-sequence;  trimmed or padded with 0 (no field)
```

BL = `hamr:bit-size`, TW = `hamr:track-size`, pad = `hamr:track-padding`, NPS = the laser-to-pole distance. The
lead-in brings the laser from half a bit plus NPS before the film to x = 0; during bit b the laser centre crosses
[b·BL, (b + 1)·BL).

*Large run, from its input*: v = 50 m/s, dt = 1 fs, so step_nm = 5e-5 nm. nb = 6 / 5e-5 = 120,000.
ne = round(3 / 5e-5) + round(5 / 5e-5) = 160,000. bits_per_track = floor(199.9 / 6) = 33.
max_tracks = floor(59.9 / 58) = 1 and n_bits = 35, so 1 track. nt = 33 × 120,000 + 320,000 = 4,280,000 steps.
Track y = 29 nm. The input asks for 35 single-tone bits; the sequence is trimmed to 33. The laser is on for 4.28 ns,
then 10 ps at 300 K: the run ends at 4.29 ns.

### Laser temperature

```
T(x, y, step) = Tmin + (Tmax − Tmin) · exp(−(x − x_L)² / 2σx² − (y − y_L)² / 2σy²)   while the laser is on, else Tmin
σ = FWHM / √(8 ln 2)
peak T at a point (viewer): the same Gaussian at the closest laser position over the whole run, on any track
```

*Large run*: FWHM 30 nm gives σx = σy = 12.740 nm; Tmin 300 K, Tmax 700 K. Grain 0 lies 2.092 nm off the track
centre: peak T = 300 + 400 · exp(−2.092² / (2 · 12.740²)) = 694.6 K. The reducer's per-grain peak (the maximum over
grain frames at the grain centre) is also 694.6 K: at 50 m/s the laser moves 0.05 nm between frames.

There are two definitions of peak T: the grains view colours a grain by the reducer's per-grain maximum over frames at
its centre; the atom views and the film tooltip use the analytic maximum at each point. They agree to rounding when
frames are dense compared with σ/v.

### Field trapezoid and field box

VAMPIRE's bit loop sets |H| = |H_previous| + increment, with increment (Hmax − Hmin)/rise for the first rise steps and
−(Hmax − Hmin)/fall for the last fall steps of each bit (`update_field_time_trapz_profile`), starting from 0 after the
lead-in. So the plateau is **Hpeak = Hmax − Hmin**, not Hmax:

```
|B| = min((i + 1) / rise, 1) · Hpeak             step i = 0 ... nb − 1 of the bit
|B| = (nb − 1 − i) / fall · Hpeak                over the last fall steps
B   = |B| · sign(bit) · u                        u = sim:applied-field-unit-vector; the z field is B · u_z
```

An input that sets only `sim:applied-field-strength` keeps VAMPIRE's defaults Hmax = +1 T and Hmin = −1 T and writes
at 2 T; the large run is such an input. Its field unit vector is (0, 0, −1), so a sequence entry −1 means a field
along +z.

The field acts only inside a box around the pole, which trails the laser by NPS (`field_profiles.cpp`):

```
|x − (x_L − NPS)| ≤ head-field-x / 2    and    |y − y_L| ≤ head-field-y / 2     (box limited to the film)
```

*Large run*: the box is 14 × 22 nm, centred 5 nm behind the laser, so it spans 12 nm behind to 2 nm ahead of the
laser centre. The track centre cools through the freezing temperature 10.33 nm behind the laser (section 7), inside
the box, so a grain freezing on the track centre does so under the full field of the bit being written. Off the track
centre the freezing isotherm is closer to the laser (it is an ellipse).

### Bits written: from the log, or from the old loop

Earlier VAMPIRE versions started the bit loop after the lead-in but ended it at bits_per_track · nb steps into the
track, so with a lead-in of a bit or more they wrote `ceil((bits_per_track · nb − ne) / nb)` bits per track and left
the field off over the rest; the current loop writes every bit. hamrvis takes the count from the run's own log: the
distinct "Written bit N" lines before the first "Head position at end of track" (`written_source = "log"`). Without a
whole track in the log it uses the old formula and adds a note (`"formula"`). Bit b < w of track t takes sequence
entry t · w + b; cells b ≥ w get 0. VAMPIRE also places each track at floor(bits written so far / bits_per_track), so
with lost bits later tracks fall on earlier ones; hamrvis models that.

*Large run*: the log has 32 "Written bit" lines, and the old formula agrees: ceil((33 × 120,000 − 160,000) / 120,000)
= ceil(31.67) = 32. VAMPIRE printed "Bits to be written: 33". The viewer shows 32 of 33 and the tab "32 bits, 50 m/s";
the 33rd cell has no field. *Scale writes*, made with the current loop: 33 of 33.

### The head model against VAMPIRE's record

At reduce time each write is compared with whatever VAMPIRE recorded (`Hamr.check`), shown under Run details:

- `B_applied` of the output table, at every row whose printed time resolves a single step. VAMPIRE prints 6
  significant digits, so beyond about 10^6 steps a printed time no longer names one step; such rows are left out and
  counted. Agreement within a relative 1e-5.
- Logged head positions: after bit n the log says "Head moving to position: x nm, y nm"; the model gives
  x = −BL/2 − NPS + (j + 1) step_nm with j = ne + (n mod w + 1) nb + 1. In a job continued from a checkpoint only lines
  after "continued from checkpoint" count, since the replayed loop logs a stale head.

With neither the check says there is nothing to compare; any difference also becomes a note on the run and a warning
in the reduce log.

*Large run*: no `B_applied` column; all 32 logged positions agree, largest difference 1.0e-4 nm in x (the log's
printing precision), 0 in y. Scale writes: all 33 positions, the same. *chain*: all 1,953 `B_applied` rows agree,
largest difference 2.2e-16 T; 3 positions, 8.9e-16 nm.

Checked by `FieldModel.test_matches_vampire_output` (to 1e-12 T), `test_bits_lost_to_a_long_lead_in`
(`tests/fixture-leadin`: 6 reported, 4 written), `test_fixed_bit_loop_writes_every_bit`,
`test_default_peak_is_hmax_minus_hmin`, `TrackRules.*`, `HeadCheck.test_a_wrong_model_is_caught`,
`NewScenarios.test_head_field_matches_vampire_in_every_write`, `Reduce.test_viewer_model_matches_python` (`Model`
against `Hamr`, in node), and the viewer readout tests (clock against the independent `Head`).

## 7. Write analysis (interpreted)

The reducer turns each grain's m_z = m_z·|m| trace into an outcome (`reduce.write_analysis`; thresholds
`THRESHOLDS`, changed with `--decision`, `--heated`, `--settled`, `--full-field`). It is computed in one pass over the
frames (`GrainAccumulator`): a grain's decision depends on the sign of its final m_z, unknown until the last frame, so
both signs are followed. `Scale.test_accumulator_matches_the_whole_series` checks it against whole-array formulas.

### No magnetic atoms, melted, decided, undecided

| Quantity | Definition | Default |
|---|---|---|
| no magnetic atoms | a grain with 0 magnetic atoms (from the coordinates; without them, \|m\| = 0 in every frame), such as a non-magnetic layer VAMPIRE gives a grain ID, for which `grain.txt` holds 0. Category `empty`, left out of everything below | |
| final direction | sign of m_z at the last frame | |
| melted | \|m\| dropped below *heated* at some frame **and** m_z was below *decision* in its final direction at some frame (it lost its direction). Grains at the edge of the hot band can dip in \|m\| without losing their direction | 0.1, 0.05 |
| decision moment | the frame after the last frame with (final sign) · m_z < *decision*: after it the grain never comes back | 0.05 |
| undecided | melted, and \|final m_z\| < *settled*: still disordered at the end | 0.3 |

*Large run, grain 0*: it stays disordered while the laser passes, and its m_z last sits below 0.05 at 2.871 ns. The
decision moment is the next frame, 2.872 ns (frame 2871, step 2,872,000), when the laser centre was 9.397 nm ahead
of it and the laser temperature at its centre was 600.6 K. Final m_z +0.753.

### Categories: the head field at the decision

```
B_head = head field at the decision step (whether or not the box covers the grain)
zero-field ("field off")   if B_head = 0
reversal ("field reversal") if |B_head| < full_field · Hpeak  (0.4: during a ramp)
full-field ("full field")  otherwise
also stored: Bz at the grain (0 outside the box), T at the centre, distance behind the laser x_L − x_g, offset y_g − y_L
```

| Category | Large run | ✓ | Scale, window | ✓ | Scale, full | ✓ |
|---|---|---|---|---|---|---|
| not melted | 259 | | 262 | | 262 | |
| full field | 144 | 116 | 146 | 119 | 146 | 120 |
| field off | 11 | 0 | 6 | 0 | 7 | 0 |
| undecided | 1 | | 1 | | 0 | |
| no magnetic atoms | 1 | | 1 | | 1 | |
| field reversal | 0 | | 0 | | 0 | |

"Full field" means the head field, not the field at the grain: in the large run 14 of the 144 full-field grains
decided outside the box (5 of them ✓). The inspector shows both values and says "(outside the field box)". *chain*:
grain 1 is full field (head −1.5 T) with 0 T at the grain.

### Freezing temperature

```
Tc = median laser temperature at the grain centre at the decision moment,
     over melted grains with |final m_z| ≥ settled whose temperature then was above Tmin + 5 % (Tmax − Tmin)
     (decisions under the laser, not thermal flips far from it);   Tc_n = number of those grains
--tc T replaces it ("given with --tc")
```

It is an effective value for this head speed and thermal profile, read from the run, not the equilibrium Curie
temperature. If it does not lie between Tmin and Tmax, the grains are scored with no lag and the run gets a note.

| | Large run | Scale, window | Scale, full | chain |
|---|---|---|---|---|
| Tc | 587.9 K | 590.3 K | 592.4 K | 488.9 K |
| Tc_n | 155 | 151 | 153 | 1 |
| freeze lag | 10.334 nm | 10.202 nm | 10.085 nm | 2.749 nm |

A Tc from one grain (chain) is a weak estimate, which is why Run details gives the number of grains.

### Freeze lag and the intended bit

A grain freezes where the laser has cooled it through Tc, a distance behind the laser centre. On the track centre:

```
lag = σx · √(2 ln((Tmax − Tmin) / (Tc − Tmin)))         0 if Tc is not between Tmin and Tmax (with a note)
bit of grain g      b = floor((x_g + lag) / BL)       on its track t = floor((y_g − pad) / TW), clipped
intended direction  = sign(sequence[t · bits_per_track + b] · u_z),   0 outside 0 ≤ b < bits_per_track
```

So bit b's cell covers [b·BL − lag, (b + 1)·BL − lag): the cells sit a lag upstream of where the laser was when the
bit's field was on.

*Large run, grain 0*: lag = 12.740 · √(2 ln(400 / 287.86)) = 10.334 nm. Grain 0 at x = 126.203 nm:
b = floor(136.537 / 6) = 22 (bit 23). sequence[22] = −1 and u_z = −1, so intended +z. At its decision the head was on
the plateau of bit 22 with B = −2 T, so Bz = +2 T inside the box: full field. Final m_z = +0.753: ✓. Its distance
behind the laser at the decision, 9.40 nm, is close to the lag; over all melted grains the median is 8.37 nm, since
grains off the track centre freeze nearer the laser.

### The score

```
correct  = category not unheated, empty or undecided;  intended ≠ 0;  sign(final m_z) = intended
score    = correct / melted grains with intended ≠ 0          "melted grains correct"
left out = melted grains with intended = 0 (outside every bit cell), counted beside the score
```

| | Large run | Scale, window | Scale, full | chain |
|---|---|---|---|---|
| melted | 156 | 153 | 153 | 5 |
| left out (outside the bits) | 14 | 10 | 9 | 2 |
| score | 116 / 142 (81.7 %) | 119 / 143 (83.2 %) | 120 / 144 (83.3 %) | 2 / 3 |
| score with no lag | 85 / 151 | 87 / 153 | 87 / 153 | 1 / 4 |

In the large run the 14 left out lie in bit 33 (never written) or beyond the 33 bit slots. The two scale writes start
from the same relaxed state and differ only by thermal noise; their scores and lags agree within that.

### The lag control, rescoring and the score-against-lag chart

The viewer scores the write again at the lag chosen in the Write analysis section (`scoreGrains`, a port of the
reducer's assignment): category, decision and final sign do not depend on the lag; only the bit, the intended
direction and "correct" change. Options: *inferred* (the bundle's lag, or *from --tc*), *none* (0 nm), *custom*,
*from bit centres*: bit b is centred at (b + 0.5)·BL − lag, so given centres x_b the lag is the mean of
(b + 0.5)·BL − x_b. This is viewer state only; the bundle never changes. Every lag-dependent view follows it (grain
outcomes, "m_z vs intended", bit cells, bit polarity, the bits lane, the read-back's intended transitions, Save PNG) and
states the lag and its source.

The chart evaluates 151 lags from 0 to max(3 BL, 2 × inferred lag). *Large run*: highest 117 / 142 at 9.64 nm; the
inferred 10.33 nm gives 116 / 142. The denominator moves with the lag because cells move over the film's ends. An
alternating sequence repeats every 2 BL (12 nm here), so the curve does too: a high score cannot tell lag L from
L + 2 BL. The inferred lag comes from Tc, not from the curve, so it is not affected; a custom or bit-centre lag can
land on the wrong peak.

Checked by `Bundles.test_viewer_scores_grains_as_the_reducer` (every grain's bit, intended direction and outcome at
the bundle's lag), `test_viewer_rescores_at_another_lag`, `ViewerFidelity` (the score at the inferred lag and at none,
recomputed in Python) and `ViewerReadouts` (the chart's series).

### Bit cells and bit polarity

Bit polarity measures a bit, not a grain: the mean m_z over the bit cell inside the band hot enough to be written.

```
band half-width  h = σy · √(2 ln((Tmax − Tmin) / (Tc − Tmin)))       (whole track width TW if there is no Tc)
cell (t, b)      x in [b·BL − lag, (b + 1)·BL − lag),   |y − (pad + (t + 0.5) TW)| ≤ h
parts            grain-map pixels whose centre lies in the cell, counted per grain: on the film only
cell m_z         Σ pixels_g · m_z,g / Σ pixels_g                     at the current grain frame
bit matches      intended ≠ 0, the cell has pixels, sign(cell m_z) = intended
```

Without a grain map, grains whose centre lies in the cell count, weighted by atom count. The overlay draws the cell
borders shifted by the lag, the band in solid lines and, with *Bit polarity*, each cell's outline in the colour of its
mean (grey below |m_z| 0.05) with ✓ or ✗.

*Large run*: h = 10.334 nm (equal to the lag, as σx = σy). 29 of 31 cells match at the end; 31 rather than 33 because
cell 1 lies off the film (x −10.33 to −4.33 nm) and cell 33 was never written. At no lag: 23 of 32. Cell 23 holds
1,405 pixels of several grains (grain 0: 297, grain 25: 231, grain 256: 224, grain 102: 206, ...); mean m_z +0.469,
intended +z: ✓. *Scale writes*: 31 of 32 each (all 33 bits written; cell 1 off the film); at no lag 23 and 26 of 33.

Checked by `ViewerFidelity`: the page's cell list (every cell's grains and pixel counts) and bit count at the end
equal a Python recomputation from the bundle at both lags.

### Bits lane, "m_z vs intended", space–time, decision positions

- **Bits lane**: one row per track. Each bit at [b·BL − lag, (b+1)·BL − lag) shows the intended arrow and, under it,
  the share of its melted grains (by grain centre) whose current sign matches. Counts stay grey until one of its grains
  has passed its decision step, since before that a grain matches only because of its starting state. The box marks the
  bit under the head. Hover gives "✓ at t" (current sign) and "✓ at end" (the outcome).
- **m_z vs intended**: m_z × intended on a green–ochre map; grains with intended 0 or never melted are grey. In a
  later run it uses the reference write's analysis.
- **Space–time**: one row per melted grain (all grains with the box ticked), sorted by track then x, time left to
  right; the dashed line is the step when the laser centre passed each grain's x, ticks are decision moments, strips
  along the top are the bit windows. Above 16,384 frames or 64M cells it is resampled: each image pixel holds the
  extreme of its frames and grains furthest from the pixel before it.
- **Decision positions**: each decided grain in a cell at (distance behind the laser, offset from the track centre),
  ✓ or ✗, with the Tc isotherm x = σx √(2 ln((Tmax − Tmin)/(Tc − Tmin)) − y²/σy²). Grains decide a little inside it on
  average, which is why the median measured distance (8.37 nm in the large run) is below the lag.

## 8. Read-back runs (recorded)

A dipole run writes `dipole-field.txt`: one row per macrocell with ID, atoms in the cell, the atoms' centroid
(x, y, z), moment, and Bx, By, Bz in T. The centroid is not the cell centre, so hamrvis grids by cell ID
(`reduce_dipole`):

```
cell size per axis   cells:macro-cell-size or -x/-y/-z (last line wins; default 10 Å)
cells per axis       n = round(L / cell) + 1             IDs x-major, z fastest:  id = (ix · ny + iy) · nz + iz
stored               dipole.f32 [3, nz, ny, nx],  dipole_atoms.u16 [nz, ny, nx]    (0 atoms: a field-probe cell)
cell i spans         i · cell to (i + 1) · cell
```

VAMPIRE numbers L/cell + 1 cells per axis, so the last cell per axis starts at the film edge. The viewer draws
(nx − 1) × (ny − 1) cells and offers nz − 1 sensor heights. The reducer refuses a file whose row count or ID order does
not match.

- **Stray B_z**: B_z of one cell layer (*Sensor z*, default round(0.55 Lz / cell)), on a green–purple map scaled to
  the 99.5th percentile of |B_z| in that layer, with the grain boundaries multiplied over it at 18 % for orientation.
  The hover gives Bx, By, Bz of the cell under the pointer.
- **B_z x–z**: B_z of the cell row at the section line, z vertical (stretched), with dashed lines at the sensor height
  and the film top = z0 + (nz − 0.5) · a_z of the snapshot grid (this run's, or the reference write's).
- **Read-back**: for each x cell, the mean B_z over cell rows whose centre lies within half the reader width of the
  section line, at the sensor height. Dotted lines mark the intended transitions of the reference write under the
  section line: x = b·BL − lag wherever sequence entry b differs from entry b − 1, at the lag in use.

*Large run, read*: 201 × 61 × 21 cells of 1 nm. Cell layers 0–9 hold atoms, 10–20 none (probe cells). The default
sensor height is cell round(0.55 × 20) = 11, shown as 11.5 nm (the cell centre); the film top is 39.5 × 0.25 =
9.875 nm. Colour scale ±289 mT. The default reader band, y = 29 ± 5 nm, takes cell rows 24–33; the waveform runs from
−181.6 to +211.7 mT (scale reads: −203.8 to +213.6 mT after the window write, −209.5 to +206.4 mT after the full one).

Checked by `Parsing.test_dipole_cells_default_to_vampire_size`, `Bundles.test_reads_compare_with_the_write_before_them`,
verify's "dipole cells" checks, `ViewerReadouts` (stray-field tooltips, the waveform series) and `ViewerPixels` (the
stray-field map's pixels).

## 9. Every readout

"Tested" names the test that compares the printed or drawn value with an independent computation from the bundle.

| Where | Readout | Kind | Comes from | Shown to | Tested |
|---|---|---|---|---|---|
| Header | atoms, grains, size | recorded | atom count from `atoms-coords.meta`; grain count; system size | exact | |
| Tab | bits and speed / integrator and temperature | recorded / reconstructed | bits written × tracks, v; input keys | exact | |
| Film tooltip | x, y | | pointer position | 0.1 nm | ViewerReadouts |
| | column m_z, atoms in column, lowest column \|m\| | recorded | top[z] at the nearest snapshot; top_count; min_m / 255 | 0.01 | ViewerReadouts |
| | laser T now, peak laser T | reconstructed | head model at the pointer | 1 K | ViewerReadouts |
| | bit n mean m_z ✓/✗ | interpreted | cell m_z of the cell under the pointer | 0.001 | ViewerFidelity |
| | grain, m_z, \|m\|, outcome | recorded / interpreted | grain map; grain series at the frame; analysis | 0.001 | ViewerFidelity, ViewerReadouts |
| Bits lane | intended, ✓ counts, written at | interpreted | scoring at the lag in use | exact | ViewerReadouts |
| Inspector (grain) | x, y; magnetic atoms | recorded | centres and counts from the atoms | 0.1 nm | |
| | m_z, \|m\|, final m_z, lowest \|m\| | recorded | grain series (int16); analysis | 0.001 | ViewerFidelity |
| | intended, bit; freeze lag | interpreted | scoring at the lag in use | exact | scoring tests |
| | decided, T then, behind laser, off centre | interpreted / reconstructed | decision step; T, distance and offset at it | 1 K, 0.1 nm | |
| | head B_z then; B_z at grain then | reconstructed | head field · u_z; 0 outside the box | 0.01 T | ViewerFidelity |
| | peak T | reconstructed | reducer's maximum over frames at the centre | 1 K | |
| | life story: T, \|m\|, m_z, B_z | recorded / reconstructed | around the hottest frame (−3σx/v to +6σx/v) or the whole run, ≤ 300 points | chart | ViewerReadouts |
| Inspector (run) | grains correct now / at end | interpreted | now: current sign; at end: outcome | exact | ViewerFidelity |
| | bits correct now / at end | interpreted | bit polarity over cells with pixels and intended ≠ 0 | exact | ViewerFidelity |
| | freeze lag, freezing T | interpreted | lag in use; Tc and its source | 0.01 nm, 1 K | ViewerFidelity |
| | head field, speed, bits written, bit length, laser T | reconstructed | head parameters | as printed | bundle tests |
| Clock | t / end; laser x; bit; B_z head | reconstructed | frame time; head position, bit and field · u_z | ps or ns, 0.1 nm, 0.01 T | ViewerReadouts |
| Time bar | M_z, M_l, bit windows | recorded / reconstructed | the output table's `M_norm_z × M_norm_l` and `M_norm_l` placed by `Real_time`; bit windows from the head model | chart | ViewerReadouts |
| Space–time | grain, t, m_z | recorded | grain series (or its overview) | chart | ViewerReadouts |
| Read-back lanes | B_x, B_y, B_z, waveform | recorded | dipole cells (section 8) | 0.1 mT | ViewerReadouts |
| Write analysis | score, Tc, lag, head field; outcome table; charts | interpreted | section 7 | exact counts | ViewerFidelity, ViewerReadouts |
| Run details | source, versions, bits, checks, notes, input | recorded | the manifest; a warning when the bundle's hamrvis is older than the current analysis model (0.9.1) | as stored | `test_viewer_reads_version_1_manifests` |
| Film, section, stray map | every pixel | as its value | the colour map applied to the bundle's value at that film position | one colour step, one screen pixel | ViewerPixels |

The time bar places frames by index while the magnetisation trace is placed by time. That is exact for evenly spaced
frames; after `--grain-every` keeps an extra last frame the two drift slightly. Exported pages (`hamrvis export`) use
a thinned bundle (every 4th grain frame by default), so their grain frames, and any "now" count, are coarser than the
bundle's; the analysis numbers are the bundle's.

## 10. Tests, tolerances and limits

### How fidelity is tested

| Level | What is compared | Where |
|---|---|---|
| 1 · bundle against raw files | every array and number of the bundle against the raw run files, read by `hamrvis/verify/reference.py`, which shares no code with the rest of hamrvis (`hamrvis verify`; table below) | `tests/fidelity/`, CI, and any user's bundle |
| 2 · page against bundle | what the running page computes and prints (lag, cell list, scores, hover grain over a grid of points, inspector values, clock, tooltips, bits lane, chart series), read back from the page in headless Chromium and compared with values worked out from the bundle and the independent head model | `tests/test_viewer_smoke.py` (ViewerFidelity), `tests/test_viewer_readouts.py` |
| 3 · pixels | sampled canvas pixels of the film (every view and colour mode), the x–z section and the stray-field map decoded through the colour map against the bundle's value there; a selected grain's outline on its grain-map boundary | `tests/test_viewer_pixels.py` |
| 4 · known answers | 0 K same-step snapshot = `grain.txt` exactly (`tests/fixture-t0`); a generated large run with a known answer per grain (`tests/make_large.py`); `B_applied` to 1e-12 T | `tests/test_hamrvis.py` |
| 5 · coverage | the scenarios in `tests/scenarios/` (chain, continuous film, custom names, file per process, LLG relax, rewrite, split write, a write continued on fewer MPI ranks, text snapshots, two tracks, write only), version 1 and 2 bundles, window runs, and real bundles where present | all of the above |

Levels 2 and 3 run on every scenario and, where the bundles are present, on real runs. At full size the large run's
scale bundles (relax, write and read each) passed every one of their 48 `hamrvis verify` checks.

### Tolerances (`hamrvis verify`)

| Quantity | Tolerance |
|---|---|
| output table | exact (float64); joined jobs: every stored row is a row of one job, in time order |
| grain list | IDs and listed atom counts exact |
| grain centres | periodic mean of the grain's magnetic atoms, stored to 1e-3 nm; a grain without magnetic atoms keeps its listed centre exactly |
| atoms per column and layer cell | exact (capped at 255) |
| grid | pitch, columns, layers and layer cells from the atoms, exactly |
| grain map | the grain with most magnetic atoms in each column, −1 without any, exact |
| grain series | int16: within half a step (0.5/32767); float32: float32 of `grain.txt`'s value |
| snapshot columns and layers | int8: 0.5/127; int16: 0.5/32767; float32: float32 rounding; empty pixels must be empty; window snapshots rebuilt from the complete snapshot before them |
| snapshots left out | none (a bundle that left snapshots out fails) |
| lowest \|m\| | within one step of 1/255 (hamrvis keeps the column \|m\| in float32) |
| dipole cells | float32 of `dipole-field.txt`'s fields; atoms per cell exact |
| head parameters | input values in VAMPIRE's units, relative 1e-12 |

Printed numbers in the viewer tests: the reference value rounded to the digits shown (half a unit of the last digit).
Derived quantities (head, Tc, lag, outcomes, bit polarity): recomputed from their definitions; the head also against
VAMPIRE's `B_applied` and logged positions.

### Known limits

- One lag across the whole track: grains off the track centre freeze closer to the laser (the isotherm is an ellipse)
  and can take the previous bit's field. A y-dependent lag was tried on the large run and scored worse (97 / 145
  against 116 / 143 under the analysis of that time), so the uniform lag stays.
- A grain's centre temperature is a rough proxy for when it froze: the freezing front needs about 0.4 ns to cross a
  grain at 20 m/s.
- Grains deciding during a field reversal can go either way. "Full field" refers to the head field at the decision,
  whether or not the box covers the grain (section 7).
- The score against lag repeats every two bits for an alternating sequence (section 7).
- Lowest |m| in the atom views is the lowest over the bundle's snapshots. A run with one snapshot (the large run's
  write) shows only the final state there; the grains view of the same colour uses the full grain series.
- Multi-track runs from the old bit loop: the head model puts later passes on top of earlier tracks, as VAMPIRE did,
  but the analysis assigns each grain to its nominal track and that track's first pass, although a later pass wrote
  over it. For an alternating sequence with an even number of bits per track both passes ask for the same directions.
  No run of this kind has been tested.
- VAMPIRE prints 6 significant digits, so beyond about 10^6 steps printed times no longer resolve single steps:
  `B_applied` rows there are left out of the head check, and snapshot times lose step precision.
- The time bar places frames by index (section 9); the m_z colour map saturates at |v| = 0.85.
- `atomistic_dipole_field.txt` (the per-atom read-back field) is not reduced.
- Bundles are lossy by design: non-magnetic atoms are dropped, pixels hold means, spins are int8. The analysis uses
  the full-precision `grain.txt`, not the stored arrays.
