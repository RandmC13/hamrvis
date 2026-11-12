# hamrvis user guide

This guide covers everything hamrvis 1.0 does, from the VAMPIRE input to the last control of the viewer. How each
number is worked out, and where it comes from, is in [how-it-works.md](how-it-works.md); the files of a bundle are in
[bundle-format.md](bundle-format.md); the design choices behind it are in [decisions.md](decisions.md).

`hamrvis` below means the installed command; `python3 -m hamrvis` from a clone is the same.

Contents: [Install](#install) · [Quick start](#quick-start) · [Preparing VAMPIRE runs](#preparing-vampire-runs) ·
[inspect](#inspect-what-is-there) · [reduce](#reduce-make-a-bundle) · [On a cluster](#on-a-cluster) ·
[verify](#verify-check-a-bundle-against-the-raw-files) · [view](#view-the-viewer) · [export](#export-one-page-to-share) ·
[Runs and tabs](#runs-and-tabs) · [Workflows](#workflows) · [Runs it adapts to](#runs-it-adapts-to) ·
[How much is kept](#how-much-is-kept) · [Messages](#messages-and-troubleshooting)

## Install

hamrvis needs Python 3.8 or newer and numpy 1.17 or newer; nothing else. The viewer runs in any current browser.

```bash
pip install "git+https://github.com/RandmC13/hamrvis.git@hamrvis-v1.0.0"   # a release
# or from a clone:
git clone https://github.com/RandmC13/hamrvis.git && cd hamrvis && pip install .
```

Without installing, `python3 -m hamrvis ...` works from the top folder of a clone, as long as numpy is installed. On
a cluster without pip access, a clone and the cluster's Python module are enough (see [On a cluster](#on-a-cluster)).
`hamrvis --version` prints the version.

## Quick start

```bash
hamrvis inspect runs/                 # what was found, what each view needs, the estimated size
hamrvis reduce runs/ -o bundles/run   # one tab per run
hamrvis view bundles/run              # opens http://127.0.0.1:8737/
```

The repository's small example chain (a Monte Carlo relaxation, a write and a read-back) tries all of it in a few
seconds: `hamrvis reduce tests/scenarios/chain -o my-first-bundle`.

## Preparing VAMPIRE runs

hamrvis reads only what VAMPIRE writes, so the views a run gets depend on its input. To get every view from a HAMR
write, add:

```
output:column-headers
output:real-time
output:magnetisation
grain:real-time
grain:magnetisation
grain:output-rate = 1000
create:grain-shape-output
config:atoms
config:atoms-output-rate = 5000
config:output-format = binary
```

and `dipole:output-dipole-field` to a read-back run.

| Input keywords | Gives | Why |
|---|---|---|
| `output:column-headers` | Every run (required) | hamrvis finds the columns of `output` by their headers |
| `output:real-time`, `output:magnetisation` | Time line under the film, with the run's M_z and M_l | The time axis of a run; Monte Carlo runs use steps instead |
| `grain:real-time` (or `grain:time-steps`), `grain:magnetisation`, `grain:output-rate` | Grain view, grain life stories, space-time diagram, write analysis | `grain.txt` is the only per-grain time series; a time column is also what lets jobs continued from a checkpoint be joined |
| `create:grain-shape-output` | Grain list and outlines (`grain-coordinates.txt`) | VAMPIRE writes it for granular systems only; a continuous film gets atom views only |
| `config:atoms`, `config:atoms-output-rate`, `config:output-format = binary` | Atom-column view, one-layer view, x–z section, lowest-\|m\| map, exact grain outlines | Spin snapshots; binary is smaller and faster to read than text, and either works |
| `dipole:output-dipole-field` (read-back) | Stray-field map, field section, read-back waveform | The field VAMPIRE computed above the film |
| `output:applied-field-strength` | Not needed | The head field is rebuilt from the input exactly as VAMPIRE applies it; if the column is there it is used to check that |

**Choosing the rates.** `grain:output-rate` sets the time resolution of every grain view and of the write analysis
(the decision moment of a grain is found to one grain frame); `config:atoms-output-rate` sets how often there is an
atom picture. Snapshots dominate the disk use: one binary snapshot of an 8.7M-atom film is about 210 MB. `hamrvis inspect` on a short test job shows what a
choice costs. Keep `hamr:active-window` snapshots as they are: hamrvis reads them (see
[Runs it adapts to](#runs-it-adapts-to)).

**Chained runs.** Run each stage (relax, write, read, anneal, ...) in its own folder and chain them with checkpoints as
usual; hamrvis orders them by the time each started. Keep the same output keywords in jobs that continue one run.

## inspect: what is there

```bash
hamrvis inspect <runs...> [reduce options]
```

Lists every run it finds, in tab order, and for each: its kind and tab title, folder, start time, VAMPIRE program and
integrator, the input and output file names, materials, output columns, whether grain views, atom views, the head
model and the field map are possible (and the input keyword that would add a missing one), the grid and the snapshots
kept, the memory a reduce would need, and the estimated bundle size and reduce time. It reads only small files, so it
is quick on any run. It takes the same options as `reduce`, so the estimate follows them. Run it first, and include its
output when reporting a problem. Example (the example chain's write):

```
write: hamr tab "Write"
  folder        tests/scenarios/chain/write
  started       2026-11-04 13:27:27
  program       hamr-simulation, integrator llg-heun
  files         input file 'input', output table 'output'
  materials     m.mat
  output table  6 columns (Real_time, B_applied, s0M_norm_x, s0M_norm_y, s0M_norm_z, s0M_norm_l)
  grain views   yes: 9 grains x 390 frames in grain.txt (0.0 MB in the bundle)
  atom views    yes: 2624 atoms in 1 coordinate file(s), 11 snapshots, 0.7 MB of spin files
                grid ~40 x 20 columns of 0.300 nm, layers ~20 x 10 x 4; 11 snapshots kept, ~0.0 MB in the bundle
  memory        11 snapshot worker(s) of ~0.10 GB (one per CPU); up to ~1.5 GB for the coordinates; ...
  head          1 track(s) x 3 bits, 3 written per track, bit 3 nm, track 5 nm, speed 1500 m/s
```

A line reading `NOT USABLE: <reason>` or `no: <reason>` says what is missing and how to add it.

## reduce: make a bundle

```bash
hamrvis reduce <runs...> -o <bundle> [options]
```

`<runs>` are run directories, or folders whose subfolders are run directories, in any mix; every run becomes a tab
(see [Runs and tabs](#runs-and-tabs)). `<bundle>` must be a new or empty directory, or an existing bundle, which is
replaced. The reduce prints a log (`[hamrvis HH:MM:SS] ...`): per run, its kind and folders, the snapshot progress with
the time left, the checks it made on itself, and for a write the inferred freezing temperature, freeze lag and grain
outcomes. It ends with the bundle size.

### Options

All are optional; the defaults suit most runs.

| Option | Default | Effect |
|---|---|---|
| `-o`, `--out DIR` | (required) | Bundle directory |
| `--order DIR ...` | log start time, then name | Tab order by folder name, as `inspect` lists them |
| `--name TEXT` | the runs' common folder name | Name shown at the top of the viewer and in saved file names |
| `--every N` | 1 | Keep every N-th spin snapshot (the last is always kept) |
| `--no-snapshots` | off | Leave the spin snapshots out: grain views and grain shapes only |
| `--no-volume` | off | Leave the per-layer data out (one-layer view and x–z section): about 10 times smaller |
| `--precision 8\|16\|32` | 8 | Snapshot grids as int8 (steps of 1/127), int16 (1/32767) or float32 |
| `--grain-precision 16\|32` | 16 | Grain series as int16 (steps of 1/32767) or float32 (grain.txt's values) |
| `--grain-every N` | 1 | Keep every N-th grain.txt frame (and the last); the write analysis then uses the kept frames |
| `--layer-factor N` | 2, more for large systems | Layer cell size in atomic sites; 1 stores every atom (one site per layer can be empty in non-cubic lattices, so 1 can leave holes) |
| `--top-max P` | 4000000 | Pixels per snapshot before the atom-column map is coarsened |
| `--layer-max C` | 2500000 | Cells per snapshot before the layer grid is coarsened |
| `--tc K` | inferred from the data | Freezing temperature for the write analysis; sets the freeze lag through the laser profile |
| `--decision X` | 0.05 | Write analysis: m_z/\|m\| below this in the final direction counts as not yet decided |
| `--heated X` | 0.1 | \|m\| below this at some point means the grain was melted |
| `--settled X` | 0.3 | \|final m_z/\|m\|\| below this means still disordered at the end (undecided) |
| `--full-field X` | 0.4 | \|B\| / peak field at the decision above this is "full field", below it "field reversal" |
| `-j`, `--workers N` | CPUs − 1, limited by memory | Parallel snapshot workers |
| `--mem GB` | free memory, or what is left of the Slurm job's limit | Memory budget for the snapshot workers |
| `--input-file NAME` | `input`, else the name in the log | The runs' input file name, if VAMPIRE ran with `--input-file` and it is not found |
| `--output-file NAME` | `output`, else the file with VAMPIRE's header | The runs' output table name, if VAMPIRE ran with `--output-file` and it is not found |

The four analysis thresholds are explained in [how-it-works.md](how-it-works.md); they change the analysis only, not
the stored data.

### Size, time and memory

In one 8.7M-atom write, 430 snapshots (90 GB) became a 2.2 GB bundle. `--no-volume`
cuts most of the rest, `--every` and `--grain-every` cut proportionally. `inspect` estimates the size and time.

Snapshots are reduced in parallel, one worker per CPU with one left free, no more than there are snapshots, and no
more than fit the memory budget (the log says which limit applied). The reducer holds at most one frame of `grain.txt`
and 16.8M atoms of a spin file at a time, so a run's size is limited by disk, not memory. If a worker is killed for
lack of memory, the reduce stops with a message; run it again with a lower `-j` or `--mem`.

### Interruptions

An interrupted reduce (time limit, out of memory, Ctrl+C) carries on when the same command is run again: finished runs
and finished snapshots are kept (`carrying on the interrupted reduce ... k of n runs were complete`). The folder
`.reduce-work` inside the bundle holds what it needs meanwhile and is removed at the end. Running it with other runs
or other options starts again.

## On a cluster

Reduce where the data is, then copy only the bundle home. The repository's `scripts/slurm-reduce.job` does this:

```bash
cd hamrvis                                         # a clone on the cluster
sbatch --account=<project> scripts/slurm-reduce.job <runs> <bundle> [reduce options]
PY_MODULE=<python-module> sbatch --account=<project> scripts/slurm-reduce.job ...   # if python3 with numpy needs a module
```

It asks for 16 CPUs, 16 GB and 2 hours (change the `#SBATCH` lines or pass `sbatch` options), prints the Python and
numpy versions, the reduce log, the time taken, the bundle size and an `rsync` line to copy it home. Paths are relative
to the folder you submit from. Workers follow the job's CPUs and memory limit. If the job runs out of time, submit it
again with the same arguments: it carries on. As a guide, the 8.7M-atom write above took 45 s on 32 CPUs with
full snapshots, and 351 s with `hamr:active-window` snapshots (read one after another).

`hamrvis verify` (below) also belongs on the cluster, as it needs the raw files. To look at a bundle on the cluster
without copying it, run `hamrvis view <bundle> --no-browser` there and forward the port:
`ssh -L 8737:127.0.0.1:8737 <cluster>`, then open <http://127.0.0.1:8737/> at home. The viewer reads the bundle in
pieces as it needs them, so this works over a slow connection, but a local copy is faster.

## verify: check a bundle against the raw files

```bash
hamrvis verify <bundle> [--run ID ...] [--every N] [--json FILE]
```

Re-reads the raw run files the bundle was made from (at the paths recorded when it was reduced) with an independent
reader in `hamrvis/verify/` that imports nothing from the rest of hamrvis, and compares every stored array and number:
the output table, head parameters, grain list and series, grid, column and layer atom counts, grain map and centres,
every snapshot's columns and layers, the lowest-|m| map and the dipole cells. Each check prints one line:

```
  ok   snapshot columns (11 snapshots)          21,648 values  max diff 0.00394  (tolerance 0.003937008)
```

and the run ends with `<n> checks, <k> failed, <s> s`. The tolerances are the stored roundings (half a step of the
encoding), listed in `hamrvis/verify/check.py` and in [how-it-works.md](how-it-works.md). The result is stored in the
bundle and shown under Run details.

| Option | Effect |
|---|---|
| `--run ID ...` | Check these runs only (the IDs are the folder names, as in the tabs) |
| `--every N` | Compare every N-th snapshot (the lowest-\|m\| map is then not checked) |
| `--json FILE` | Write every check as JSON |

Exit status: 0 if everything agrees, 1 if anything differs, 2 if it cannot check (not a bundle, an unknown `--run`, or
the run directories are not where they were: a bundle copied away from its runs cannot be verified). It reads every
snapshot again with a simpler reader, so it takes longer than the reduce: about 285 s against 45 s for the 8.7M-atom write.

## view: the viewer

```bash
hamrvis view <bundle> [-p PORT] [--no-browser]
```

Serves the viewer and the bundle on `127.0.0.1:8737` (only this computer can reach it) and opens it in the browser.
`-p` chooses another port, `--no-browser` only prints the address. Ctrl+C stops it. Nothing leaves the computer;
the page loads its font from Google Fonts and falls back to a system font offline.

### Layout

From the top: the run name and size (atoms, grains, film size), the tabs (one per run), the controls bar (stays at
the top while scrolling), the film axis with the zoom buttons, the lanes, the transport bar, then the inspector beside
the space-time diagram, the write analysis and the Run details. The ◐ button switches between light and dark
(remembered; the default follows the system).

### Tabs

One tab per run, in run order. Each shows the kind of run (Write, Read-back, Monte Carlo, or the VAMPIRE program for
other runs) and a line with the folder name and a summary: bits and head speed for a write, integrator and temperature
for a Monte Carlo or other run, "dipole field" for a read-back. Hovering a tab shows its folders (several for a run
continued over jobs). The address ends in `#<folder>`, so a run can be bookmarked; without one the first write opens.
Bundles made before format version 2 use `#relax`, `#write`, `#read`.

### Controls bar

Controls that do not apply to the run or its data are hidden or greyed out.

| Control | Options | Notes |
|---|---|---|
| Film view | Grains · Atom columns · One layer | Grains: one colour per grain from `grain.txt`. Atom columns: the mean spin of each column of atoms through the film (about 40 atoms). One layer: one atomic layer, cells of 2 × 2 sites by default. Atom views show the nearest snapshot |
| Layer | slider, `z <x> nm` | One-layer view only: which atomic layer |
| Smooth | off · 2 × 2 · 4 × 4 · 8 × 8 | Atom views coloured by m_z or m_z vs intended. Averages neighbouring pixels on screen (empty pixels stay empty; the data and the tooltip are unchanged). A column averages only ~40 atoms, so a melted region looks speckled rather than white without it |
| Colour by | m_z · m_z vs intended · Peak T · Lowest \|m\| | m_z: blue (down) to red (up). m_z vs intended: green where it matches the bit written there, ochre where not, grey for grains never melted (writes only). Peak T: the highest laser temperature each point reaches in the run (head model needed). Lowest \|m\|: the lowest magnetisation length each grain or column reaches in the run, which shows where the film melted |
| Sensor z, Reader width | sliders; width 1–40 nm (default 10) | Read-back runs: the height of the stray-field map and the width of the band the read-back waveform averages |
| Bit cells | checkbox | Writes: every bit's region on every track, shifted upstream by the freeze lag, and the band hotter than the freezing temperature (the part of the track that can be written) |
| Bit polarity | checkbox (after Bit cells) | Outlines each cell's part inside that band in the colour of its mean m_z (atom-weighted over the grain pieces it holds; grey if \|m_z\| < 0.05) and marks ✓ / ✗ against the bit where the cell is wide enough; field-free `0` bits get no mark |

### Film axis and zoom

Every lane shows the same x range of the film. Zoom with ctrl + wheel on a lane (a trackpad pinch sends the same),
the − + buttons or the − + keys; ↔ or 0 shows the whole film (stretched in y if need be). Pan with shift + drag or a
horizontal wheel. When zoomed in, a strip under the axis shows where the window is (drag it to move) and, for a
write, *follow head* keeps the head in view (any zoom or pan turns it off). A film too flat to read at full width
opens on its first part at its true aspect ratio; a film taller than its lane scrolls with the wheel, and its y range
is printed under the lane. Hovering any lane draws one crosshair through all of them with its x.

### Lanes

- **Bits** (writes): the intended bits of every track, one row per track as on the film (more than 8 tracks: 4 at a
  time), shifted upstream by the freeze lag, each marked ↑, ↓ or 0. Under each bit a bar shows the share of its melted
  grains written correctly (green) and wrongly (ochre), with the count; grey until written. A box marks the bit under
  the head. Hover for the bit, track, intended direction, x and when it was written.
- **Film** ("Written" in a read-back): the film in the chosen view and colour. On it: the laser isotherms (solid:
  freezing temperature; dashed, dash-dotted, dotted: 600, 500, 400 K, where they fall between the laser's minimum
  and maximum), the head-field box with its direction and strength (`↑ 0.80 T`; dashed when the field is 0), the
  selected and hovered grain outlines, the dashed section line, the Bit cells and polarity overlays, and the reader
  band in a read-back. When the snapshot shown is not at the current time, its time is printed on the image. Hover for
  the column or grain values (see *Readouts* below). Click a grain to select it (click again to deselect); clicking or
  dragging also moves the section line to that y.
- **x–z section**: a slice through the section line from the nearest snapshot, height stretched by the printed factor
  (needs the layer data).
- **Stray B_z** (read-back): the field at the sensor height, green (down) to purple (up), scaled to the 99.5th
  percentile shown; click to move the section line. **B_z x–z**: the field in the plane of the section line, dashed
  lines at the film top and the sensor. **Read-back**: the mean B_z across the reader band along x, dotted lines at
  the intended bit transitions.

A lane that cannot draw says `This lane could not be drawn: ...` instead of staying blank; data that cannot be read
gives a banner `Part of the data could not be read: ...`.

### Transport

Directly under the film: start, play/pause, end, the speed (0.1× to 30×, default 3×; frames per update, about 30
updates a second) and the time bar, which shows the bit windows and the run's M_l (grey) and M_z. Click or drag it to
seek. The clock reads `t <time> / <end>` and, for a write, the laser x, the bit being written and the head field. Monte
Carlo runs count MC steps. The bar is hidden when a run has one frame.

### Keys and mouse

| Input | Does |
|---|---|
| Space | Play / pause |
| ← → | One frame back / on (shift: 50 frames) |
| Home, End | First, last frame |
| + (or =), −, 0 | Zoom in, out, whole film |
| Esc | Deselect the grain |
| Ctrl + wheel, trackpad pinch | Zoom about the pointer |
| Shift + drag, horizontal wheel | Pan |
| Click on the film | Select a grain and move the section line |
| Click on a chart (life story, space-time, scatter) | Go to that time (and select that grain) |

Keys are ignored while typing in a field. There is no touchscreen pinch.

### Readouts and marks

Every readout label carries a mark, explained when hovered: **recorded** (VAMPIRE wrote it: m_z, |m|, atom counts,
column values, fields), **reconstructed** (hamrvis rebuilt it from the run's input, as VAMPIRE computes it: laser
temperature, head field, head position) or **interpreted** (hamrvis's write analysis: outcomes, freezing temperature,
freeze lag, bit correctness). [how-it-works.md](how-it-works.md) gives the source of each.

### Inspector

With no grain selected it shows the run's facts: for a write, grains and bits correct now and at the end, freeze lag,
freezing temperature, head field and speed, bits written, bit length and laser temperatures, with a legend for the
isotherms and the head box; for other runs, program, integrator, temperature, steps and grains. Hover a grain or click
to keep it: its outcome (✓ intended, ✗ wrong, outside bits, undecided, not melted), position, magnetic atoms, m_z,
|m|, intended bit, final m_z, decision moment and what it felt then (temperature, head field, field at the grain, how
far behind the laser and off the track centre), peak temperature and lowest |m|. Below, its life story: laser
temperature, |m|, m_z and field felt against time around its hottest moment (*whole run* shows all of it); the green
line marks the decision, the other line the current time. Click the chart to go there.

### Space–time

One row per melted grain (*all grains* for every grain), sorted by x, time left to right, coloured by m_z; dashed: the
laser centre passes the grain; tick: its decision. Track and bit boundaries are marked. Hover for the grain, time and
m_z; click to select the grain and go to that time. With more frames or grains than pixels, each pixel shows the frame
or grain that differs most from the pixel before it, so a short event stays visible (the time bar does the same).

### Write analysis

For writes: melted grains matching their intended bit (grains outside the bits are not counted), the freezing
temperature (inferred, or "--tc"), the freeze lag and the head field. The **freeze lag** control rescores everything
in the viewer only (the bundle is unchanged): *inferred* (or *from --tc*), *none*, *custom* (a number in nm) or *from
bit centres* (the x of bit 1, 2, ... as measured elsewhere, in nm); *in use* shows the value. Charts: *Decision
positions* (where each grain froze relative to the laser, ✓ / ✗, with the freezing isotherm; click a point to select
that grain), *Score against lag* (the score for every lag, the one in use and the inferred one marked), and *Outcome by
head field at decision* (full field, field reversal, field off). *How these are worked out* (collapsed) states the
rules. A bundle reduced before hamrvis 0.9.1 shows a warning that its outcomes come from the older analysis: reduce it
again.

### Run details

The source folders, when and with which hamrvis it was reduced, the viewer version, the jobs a run was continued in,
bits written and reported, snapshots (and any left out as unreadable), the freezing temperature and lag with their
source, and every note the reduce made (see [Messages](#messages-and-troubleshooting)). Under *Checks made when
reducing*: the head model against VAMPIRE's head positions and `B_applied` column, the stored column values against the
exact means, the snapshots' grain averages against `grain.txt`, how the grain outlines were drawn, and the result of
`hamrvis verify` if it was run. Then every keyword and value of the run's input.

### Save PNG

Stacks the axis and all lanes, with a header line (run, stage, time, lag), into one PNG named
`<name>-<run>-<frame>.png`. Local viewer only.

### Colours

Magnetisation: blue (down) to red (up). Fields: green (down) to purple (up), so a field never reads as magnetisation.
Written as intended: green; wrong: ochre; not melted: grey.

## export: one page to share

```bash
hamrvis export <bundle> -o page.html [--grain-step N] [--snap-step N] [--no-volume] [--limit-mb MB] [--artifact]
```

Writes one self-contained HTML page (viewer and data) that opens without hamrvis. To keep it small, grain frames are
thinned to every 4th (`--grain-step`, default 4) and snapshots to every N-th (`--snap-step`, default 1; the last is
always kept); `--no-volume` drops the layer data. It refuses to write a page over `--limit-mb` (default 15.5, for
sites with a 16 MB limit) and says what to raise. `--artifact` leaves out the document skeleton, for hosts that supply
their own. Save PNG is not in exported pages.

**Before sharing** a bundle or page: both record the absolute paths of the run folders, shown under Run details.

## Runs and tabs

`reduce` finds every run in the folders given and works out what each one is from its own files; nothing is assumed
about the workflow.

- **Kinds of run** (`hamrvis/runs.py`, `TYPES`): *Write* (`sim:program = hamr-simulation`: head model, bits lane,
  write analysis, space-time), *Read-back* (dipole field output: stray-field map, field section, read-back waveform),
  *Monte Carlo* (Monte Carlo integrators: time in MC steps). Any other run gets a tab titled with its VAMPIRE program,
  showing what it wrote (grain and atom views, time line, run details). Field maps are shown for any run that wrote
  `dipole-field.txt`.
- **Order**: the order the runs started (the first time stamp of each `log`), since chained runs run one after
  another; then folder name. `--order` sets it by hand.
- **Runs continued over several jobs** (`sim:load-checkpoint = continue`, same input otherwise, resuming at a step the
  previous job reached) are one run and one tab; where jobs overlap the later one is used. A job may run on another
  number of processes or GPUs, or with another `config:output-mode` / `output-nodes`: its snapshots are read with its
  own atom coordinates.
- **Comparisons**: a run after a write (a read, an anneal, ...) is compared with the nearest write before it: bit
  cells, bit polarity and "m_z vs intended" use that write's bits.

## Workflows

Relax → write → read is one workflow among many. These were tested end to end (tabs listed in the order they appear):

| Workflow | Folders | Tabs |
|---|---|---|
| Relax, write, read | relax, write, read | Monte Carlo · Write · Read-back |
| Write first, then anneal and read | write, anneal, read | Write · Monte Carlo · Read-back (the anneal and read are compared with the write) |
| Other programs before a write | field-cool, loop, write, read | field-cool · hysteresis-loop · Write · Read-back |
| A write continued over two jobs on 2 then 1 MPI ranks | write-job1, write-job2, read | Write (one tab for both jobs) · Read-back |
| Several reads with different dipole solvers, folder names out of order | z-write, c-read-macrocell, b-read-tensor, a-read-hierarchical | Write · Read-back · Read-back · Read-back (in run order, not name order) |
| Overwrite a written track | write-tracks, overwrite-slow, read | Write · Write · Read-back (the read is compared with the overwrite) |

## Runs it adapts to

Nothing is tied to one run. The reducer reads everything from the run's own files:
- tracks, bits per track and the bit sequence follow VAMPIRE's rules (`hamr:number-of-bits`, `track-size`,
  `track-padding`, single-tone or user-defined sequences, field-free `0` bits);
- the head field follows `sim:min/maximum-applied-field-strength`, rise/fall times and the field direction;
- `grain.txt` columns follow the input's `grain:` keywords in order (`time-steps` or `real-time` etc.); a `grain.txt`
  hamrvis cannot read (no `grain:magnetisation`, or an unknown keyword before it) leaves out the grain views with a
  note instead of stopping;
- input files are read as VAMPIRE reads them: every unit VAMPIRE knows (`!um`, `!kOe`, `!ps`, `!nm/ns`, ...), values
  without a unit in VAMPIRE's internal units (A, s, T, A/s: a bare `hamr:head-speed = 3e11` is 30 m/s), integer lists
  such as `hamr:bit-sequence` read the way VAMPIRE reads them; an unknown unit or one of the wrong kind is an error
  naming the key;
- file names: the input is `input` or the one the log names (`Opening main input file`); the output table is `output`
  or the file starting with VAMPIRE's `# Output file for vampire simulation` header. A job continued from a checkpoint
  appends to its table without that header, so a custom output name there needs `--output-file`;
- a missing material file: every material counts as magnetic (with a warning); continuous films (no
  `grain-coordinates.txt`) get atom views only, with the snapshots as frames; a last row a stopped job left
  half-written in `output` or `grain.txt` is left out;
- snapshots in binary or text, single file or one file per node or process; truncated copies are skipped with a
  warning;
- window snapshots of `hamr:active-window` runs (only the atoms near the head, `Spin data: window` in the `.meta`
  file) are read in order, each updating the snapshot before it; a truncated window snapshot stops the reduction,
  since the snapshots after it depend on it;
- the grids adapt to the lattice and size: atom columns at the atomic pitch up to 4M pixels, then coarser; layers at
  2 × 2 sites up to 2.5M cells per snapshot, then coarser;
- grain shapes come from the atoms' grain IDs; coordinate files without grain IDs (older VAMPIRE text format) fall
  back to the Voronoi outlines, labelled approximate.

Tested on real runs including 1 track of 33 bits with binary file-per-node snapshots (8.7M atoms), the same with
`hamr:active-window` snapshots, and 2 tracks of 11 bits with field-free edge bits, text snapshots and a film cut by a
geometry file (1.6M atoms).

## How much is kept

Everything VAMPIRE writes per grain (`grain.txt`, every frame) and per run (`output`) is kept at full time resolution.
Every spin snapshot is kept (unless `--every` is given), reduced to the mean spin of the magnetic atoms in each column
through the film and in each 2 × 2-site cell of every atomic layer, stored to 1/127 per component. The reducer checks
itself on every run and records the result (Run details): each stored value against the exact mean, and each
snapshot's per-grain average against `grain.txt`. `hamrvis verify` checks everything again with an independent reader.
Keep the VAMPIRE output if you may need the exact values: a bundle holds averages, rounded.

## Messages and troubleshooting

Run `hamrvis inspect` first, and include its output when reporting a problem. Errors stop the command (exit status 1);
warnings are printed in the log and the run goes on; notes are stored in the bundle and shown under Run details.

### Errors

| Message | What to do |
|---|---|
| `no run directories (folders with an 'input' file)` / `no VAMPIRE runs in ...` | Point at a run directory or at the folder holding them. If VAMPIRE ran with `--input-file`, pass `--input-file NAME` |
| `<path>: not a directory` | Check the path |
| `no VAMPIRE input file (looked for 'input' and the name in the log ...)` | Pass `--input-file NAME` |
| `no VAMPIRE output table (looked for 'output' ...)` / `several VAMPIRE output tables` | Pass `--output-file NAME` |
| `... has no header row: add output:column-headers to the input` | Rerun with `output:column-headers` |
| `... has no rows: the run stopped before its first output` | Leave that folder out |
| `--order names ..., not among the runs found` | Use folder names exactly as `inspect` lists them |
| `--every must be 1 or more` (also `--grain-every`) | Give a positive number |
| `a HAMR write needs hamr:bit-size and hamr:track-size in its input` | VAMPIRE has no usable default for them: the run could not have written bits |
| `'<unit>' is not one of VAMPIRE's units` / `VAMPIRE reads <key> as a <kind>, but '<unit>' is a <kind> unit` | The input has a unit VAMPIRE would also reject; fix the input |
| `... is not empty and holds no hamrvis bundle` | Choose a new or empty directory for `-o` |
| `A snapshot worker was killed, most likely because memory ran out ...` | Run the same command again with a lower `-j` or `--mem`: it carries on |
| `snapshot numbers repeat across the jobs of this run` | The jobs joined into this run do not continue one another: reduce them as separate folders, or check `--order` |
| `the jobs of this run have different output columns` | Make the joined jobs' output keywords the same, or reduce them separately |
| `grain.txt needs grain:real-time or grain:time-steps to join the jobs` | Add one of them: without a time column the jobs cannot be joined |
| `its atoms differ from those of the job it continues` | The folders are not one run continued: reduce them separately |
| `no magnetic atoms in the coordinate files` | The material file marks no atoms as magnetic, or the coordinate files are not from this run |
| `snapshot N is a window snapshot with no complete snapshot before it` | A window run needs its first, complete snapshot: keep every snapshot of the run |
| `dipole-field.txt has N cells, expected ...` | The file does not belong to this run's cells, or is truncated |
| `<bundle> has no manifest.json` (view) / `is not a hamrvis bundle` (verify) | That folder is not a bundle: make one with `hamrvis reduce` |
| `cannot serve on port ...` | Another program holds the port: stop it or use `-p` |
| `page.html is X MB, over Y MB: raise --grain-step or --snap-step, or use --no-volume` (export) | As it says, or raise `--limit-mb` if the host allows |
| `run directories not found` (verify, exit 2) | Run verify where the raw runs are, at the paths recorded in the bundle |
| `no run ... in this bundle` (verify, exit 2) | Use run IDs as shown in the tabs |

### Warnings and notes

| Message | Meaning |
|---|---|
| `skipping <dir>: it has an input file but no output table` | That folder is not a run (or needs `--output-file`) |
| `last row incomplete (job stopped while writing it): left out` | A stopped job's half-written row in `output` or `grain.txt` is ignored |
| `material file ... not in <dir>` / `defines no material[n] keys`, note `every material counted as magnetic` | Copy the material file into the run folder, or non-magnetic atoms count in the atom views |
| `snapshot N skipped: ...` | A truncated or damaged snapshot (a copy caught mid-transfer?) is left out; listed under Run details |
| `grain.txt not used, so no grain magnetisation views: ...` | The reason follows (no `grain:magnetisation`, an unknown keyword before it, columns that do not match the grains, a damaged row) |
| `the coordinate files hold no grain IDs (older text format)` | Grain shapes stay VAMPIRE's approximate Voronoi outlines |
| note `no grain-coordinates.txt` | Not a granular system: atom views only |
| note `grain.txt has no time column: frame times assumed from grain:output-rate` | Add `grain:real-time` for exact times |
| `head model differs from VAMPIRE's record: ...` | hamrvis's head model does not match VAMPIRE's head positions or `B_applied`: report it, with the run's input |
| note `bits written per track (w of n) taken from VAMPIRE's bit loop, as the log has no whole track` | The log stopped before a track ended; the count follows VAMPIRE's code |
| note `no freezing temperature between the laser's minimum and maximum` | No settled grain decided under the laser (or `--tc` lies outside it): grains are scored with no freeze lag |

### In the viewer

| What you see | What to do |
|---|---|
| A view or control is missing or greyed | The run lacks that data: `inspect` names the keyword that adds it |
| `No bundle` | Start the viewer with `hamrvis view <bundle>` |
| `... HTTP <n> to a byte range request (use hamrvis view)` | The bundle is served by something other than `hamrvis view`, which the viewer needs for large files |
| `... is not in this exported page` | The page was exported without that data (`--no-volume`, thinning) |
| `This run could not be loaded` / `This lane could not be drawn` | Report it with the browser console's message and `inspect` output |
| The analysis warns of an older analysis model | Reduce the run again with this version |
