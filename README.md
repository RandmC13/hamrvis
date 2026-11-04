# hamrvis

Viewer for VAMPIRE heat-assisted magnetic recording (HAMR) simulations. Give it a folder of runs (for example a
Monte Carlo relaxation, one or more HAMR writes, dipole-field reads, in any number and order) and it builds one page
with a tab for each run, showing what that kind of simulation calls for, on one down-track axis.

It has two parts:

1. **`hamrvis reduce`** turns finished run directories into a compact *bundle*. It runs where the data is (on the
   cluster), so only the bundle is copied home. The 90 GB of spin snapshots of a 8.7M-atom write become about
   about 2.2 GB, or about 180 MB without the per-layer data (estimated from the final snapshot).
2. **`hamrvis view`** serves one general viewer for any bundle on `localhost` and opens it in the browser. The
   same viewer works for every run, so improvements apply to old runs too.

Only Python 3.8+ with numpy is needed; the viewer is plain HTML and JavaScript with no build step.

## Typical workflow

On cluster, after the runs have finished (`<chain>` is a folder holding one subfolder per run):

```bash
cd ~/hamrvis                        # a clone of this repository
sbatch --account=<project> scripts/slurm-reduce.job <chain> ~/scratch/bundles/<name>
```

On the laptop:

```bash
rsync -av --info=progress2 <cluster>:scratch/bundles/<name> ~/hamr-bundles/
python3 -m hamrvis view ~/hamr-bundles/<name>         # opens http://127.0.0.1:8737/
```

For a run on your own machine, reduce it directly: `python3 -m hamrvis reduce runs/my-run -o bundles/my-run`.

To share a run as one file (for example with a colleague, or as a localhost.ai artifact):

```bash
python3 -m hamrvis export bundles/<name> -o <name>.html --grain-step 8      # must stay under 16 MB for localhost.ai
```

## Commands

| Command | What it does |
|---|---|
| `python3 -m hamrvis inspect <runs...>` | Lists the runs found in tab order and, per run, its type, file names, what each view needs and whether it is there (with the input keyword that would add it), and the estimated bundle size and reduce time. Reads only small files: run it first, and include its output when reporting a problem. Takes the same options as `reduce`. |
| `python3 -m hamrvis reduce <runs...> -o <bundle>` | `<runs>`: run directories, or a folder whose subfolders are run directories; each run becomes a tab (see *Runs and tabs*). Options: `--order DIR ...` (tab order by folder name), `--no-volume` (skip per-layer data, ~10x smaller), `--every N` (keep every N-th snapshot), `--tc K` (freezing temperature for the analysis; measured by default), `--name`, `-j N` (parallel workers); see *Reduce options* for the rest. |
| `python3 -m hamrvis view <bundle>` | Serves the viewer and the bundle on `127.0.0.1:8737` (`-p` to change; a port in use is reported). |
| `python3 -m hamrvis export <bundle> -o page.html` | One self-contained page. `--grain-step N` keeps every N-th grain frame, `--snap-step N` every N-th snapshot, `--no-volume`, `--artifact` (no document skeleton, for localhost.ai). |
| `python3 -m unittest discover tests` | Tests (needs node for the viewer-model test). |

### Reduce options

All optional; without them the bundle is the same as before they existed.

| Option | Effect |
|---|---|
| `--precision 8\|16\|32` | Snapshot grids stored as int8 (1/127, default), int16 (1/32767) or float32 |
| `--grain-precision 16\|32` | Grain series stored as int16 (1/32767, default) or float32 (grain.txt's values) |
| `--grain-every N` | Keep every N-th grain.txt frame (and the last); the write analysis then uses the kept frames |
| `--no-snapshots` | Leave the spin snapshots out (grain views and grain shapes stay) |
| `--top-max P`, `--layer-max C` | Pixels per snapshot before the column map / layer grid is coarsened (4e6, 2.5e6) |
| `--decision`, `--heated`, `--settled`, `--full-field` | Write analysis thresholds (docs/analysis.md) |
| `--input-file NAME`, `--output-file NAME` | The runs' input and output file names, if VAMPIRE was run with `--input-file` / `--output-file` and they are not found automatically (see below) |

## What the run needs to write

The reducer reads what VAMPIRE writes; the more of it there is, the more the viewer shows.

| Input keywords | Gives |
|---|---|
| `output:column-headers`, `output:real-time`, `output:magnetisation` | The time line under the transport bar (required: `output` must have column headers) |
| `grain:real-time`, `grain:magnetisation`, `grain:output-rate` | Grain view, grain life stories, space-time diagram, write analysis |
| `create:grain-shape-output` | Grain outlines |
| `config:atoms`, `config:atoms-output-rate`, `config:output-format = binary` | Atom-column view, one-layer view, cross-section, lowest-\|m\| map |
| `dipole:output-dipole-field` (read stage) | Stray-field map, field section, read-back waveform |
| `output:applied-field-strength` | Not needed: the head field is rebuilt from the input, exactly as VAMPIRE applies it |

## Runs and tabs

`reduce` finds every run in the folders given and works out what each one is from its own files; nothing is
assumed about the workflow.

- **Kinds of run** (`hamrvis/runs.py`, `TYPES`): *Write* (`sim:program = hamr-simulation`: head model, bits lane,
  write analysis, space-time), *Read-back* (dipole field output: stray-field map, field section, read-back
  waveform), *Monte Carlo* (Monte Carlo integrators: time in MC steps). Any other run gets a tab titled with its
  VAMPIRE program showing what it wrote (grain and atom views, time line, run details). Field maps are shown for any
  run that wrote `dipole-field.txt`.
- **Order**: the order the runs started (the first time stamp of each `log`), since runs chained by checkpoints run
  one after another; then folder name. `--order` sets it by hand.
- **Runs continued over several jobs** (`sim:load-checkpoint = continue`, same input otherwise, resuming at a step
  the previous job reached) are one run and one tab; where jobs overlap the later one is used.
- **Comparisons**: a run after a write (a read, an anneal, ...) is compared with the nearest write before it: bit
  cells, bit polarity and "m_z vs intended" use that write's bits.

## Runs it adapts to

Nothing is tied to one run. The reducer reads everything from the run's own files:
- tracks, bits per track and the bit sequence follow VAMPIRE's rules (`hamr:number-of-bits`, `track-size`,
  `track-padding`, single-tone or user-defined sequences, field-free `0` bits);
- the head field follows `sim:min/maximum-applied-field-strength`, rise/fall times and the field direction;
- `grain.txt` columns follow the input's `grain:` keywords in order (`time-steps` or `real-time` etc.); a
  `grain.txt` hamrvis cannot read (no `grain:magnetisation`, or an unknown keyword before it) leaves out the grain
  views with a note instead of stopping;
- input files are read as VAMPIRE reads them: every unit VAMPIRE knows (`!um`, `!kOe`, `!ps`, `!nm/ns`, ...), values
  without a unit in VAMPIRE's internal units (A, s, T, A/s: a bare `hamr:head-speed = 3e11` is 30 m/s), integer lists
  such as `hamr:bit-sequence` read the way VAMPIRE reads them; an unknown unit or one of the wrong kind is an error
  naming the key;
- file names: the input is `input` or the one the log names (`Opening main input file`); the output table is
  `output` or the file starting with VAMPIRE's `# Output file for vampire simulation` header. A job continued from
  a checkpoint appends to its table without that header, so a custom output name there needs `--output-file`;
- a missing material file: every material counts as magnetic (with a warning); continuous films (no
  `grain-coordinates.txt`, which VAMPIRE writes only for granular systems) get atom views only, with the snapshots
  as frames; a last row a stopped job left half-written in `output` or `grain.txt` is left out;
- snapshots in binary or text, single file or one file per node or process; truncated copies are skipped with a
  warning;
- window snapshots of `hamr:active-window` runs (only the atoms near the head, `Spin data: window` in the `.meta`
  file) are read in order, each updating the snapshot before it; a truncated window snapshot stops the reduction,
  since the snapshots after it depend on it;
- the grids adapt to the lattice and size: atom columns at the atomic pitch up to 4M pixels, then coarser; layers
  at 2 x 2 sites (`--layer-factor 1` stores every atom) up to 2.5M cells per snapshot, then coarser;
- grain shapes come from the atoms' grain IDs; coordinate files without grain IDs (older VAMPIRE text format)
  fall back to the Voronoi outlines, labelled approximate.

Tested on the large run cluster chain (1 track, 33 bits, binary file-per-node, 8.7M atoms) and the 2026-10-15 two-track
run (2 tracks x 11 bits with field-free edge bits, text snapshots, film cut by a geometry file, 1.6M atoms).

## How much is kept

Everything VAMPIRE writes per grain (`grain.txt`, every frame) and per run (`output`) is kept at full time resolution.
Every spin snapshot is kept (none are dropped unless `--every` is given), reduced to the mean spin of the magnetic
atoms in each 0.25 nm column through the film and in each 0.5 nm x 0.5 nm cell of every atomic layer, stored to
1/127 per component. The reducer checks itself on every run and records the result (shown under Run details): each
stored value against the exact mean, and each snapshot's per-grain average against `grain.txt`. Grain shapes are
taken from the atom positions, not from VAMPIRE's Voronoi outlines. See `docs/decisions.md` (Accuracy).

## Viewer

- **Tabs**: one per run, in run order, at the top right (title from the kind of run, then the folder name). The
  address ends in `#<folder>`, so a run can be bookmarked; bundles made before version 2 keep `#relax`, `#write`,
  `#read`.
- **Lanes** share the film's x axis; hovering any lane draws one crosshair through all of them.
  - *Bits*: intended bits of every track (one row per track, stacked as on the film), shifted upstream by the
    freeze lag, with the melted grains pointing along their bit's field at the current time (green) under each; a box
    marks the bit under the head.
  - *Film*: grains (shapes from the atoms), atom columns (through-thickness average) or one atomic layer; coloured
    by m_z, m_z against the intended bit (never-melted grains grey), peak laser temperature or the lowest |m|
    reached. Laser isotherms (bold: freezing temperature) and the head-field box with its direction and strength.
    The dashed line is the section line: click the film to move it. *Bit cells* overlays every bit's region on
    every track and the band hotter than the freezing temperature (the part of the track that can be written).
    *Bit polarity* outlines each cell's part inside that band in the colour of its mean m_z (atom-weighted over
    the grain pieces it contains, at the current time) and marks it ✓ / ✗ against the bit; the side panel counts
    correct bits now and at the end.
    Atom views show the nearest snapshot; its time is printed on the image when it differs from the timeline.
  - *Transport* directly under the film: start / play / end (Home, space, End), 0.1x to 30x, step (arrow keys,
    shift for 50 frames).
  - *x–z section*: slice through the section line from the nearest snapshot (height x2).
  - Read stage: *stray field* B_z at a chosen sensor height, *field section* (x–z), *read-back* waveform (mean B_z
    across a reader band of chosen width) with the intended transitions.
- **Inspector**: hover a grain, or click to keep it: outcome, decision moment and circumstances, and its life story
  (laser temperature, |m|, m_z, field felt). Click the chart to jump there.
- **Space–time**: one row per melted grain (or all), sorted by x, time left to right; dashed: laser passes the
  grain; tick: decision moment. Click to select the grain and jump to that time.
- **Analysis** (write): success, freezing temperature and lag, where each grain froze, outcomes by circumstance.
- **Save PNG** stacks the lanes into one image (local viewer only).

Colours: magnetisation is blue (down) to red (up); fields are green (down) to purple (up), so a field never reads
as magnetisation; written-as-intended is green, wrong is ochre.

## More

- `docs/decisions.md`: why it is built this way, the choices made with you, and what was left out.
- `docs/bundle-format.md`: the bundle files.
- `docs/analysis.md`: how the write analysis classifies grains, and the head model.
