# Bundle format (version 2)

A bundle is a directory written by `hamrvis reduce`. All binary files are raw little-endian arrays with no header;
their shapes come from `manifest.json`. Files ending in `.gz` are gzip-compressed (the viewer inflates them with the
browser's `DecompressionStream`).

```
manifest.json
<run>/output.f64              float64 [rows, cols]          the run's 'output' file; column names in the manifest
<run>/grains.i16            int16   [frames, grains, 4]   grain.txt per-grain block: mx, my, mz (direction), |m|; value x 32767
<run>/grain_lead.f64        float64 [frames, n_lead]      grain.txt leading columns (time in s, temperature in K)
<run>/grains_overview.i16   int16   [bins, grains, 4]     series over 512 MB only: per time bin lowest and highest m_z (mz x |m|), lowest and highest |m|, x 32767
<run>/top/NNNNNNNN.i8.gz    int8    [3, ny, nx]           snapshot NNNNNNNN: mean spin of the magnetic atoms in each column, x 127
<run>/vol/NNNNNNNN.i8.gz    int8    [3, nz, vy, vx]       snapshot NNNNNNNN: mean spin per layer on 2x2-site pixels, x 127
<run>/top_count.u8          uint8   [ny, nx]              magnetic atoms per column (capped at 255)
<run>/vol_count.u8          uint8   [nz, vy, vx]          magnetic atoms per volume pixel
<run>/min_m.u8              uint8   [ny, nx]              lowest column |mean spin| over all snapshots, x 255
<run>/grain_map.i16         int16   [ny, nx]              grain index of the magnetic atoms in each column (-1: none)
<run>/dipole.f32              float32 [3, nz, ny, nx]       dipole-field.txt cell fields Bx, By, Bz (T)
<run>/dipole_atoms.u16        uint16  [nz, ny, nx]          atoms per dipole cell (0: a field-probe cell in non-magnetic material)
```

With `--precision 16` the snapshot files are `NNNNNNNN.i16.gz` (int16, x 32767, empty -32768), with `--precision 32`
`NNNNNNNN.f32.gz` (float32, empty NaN); `snapshots.dtype` says which. With `--grain-precision 32`, `grains.f32`
(float32, unscaled) replaces `grains.i16` and `grain_series.dtype` is `"f32"`.

`<run>` is the run's id (its folder name, made unique), one per tab. A run without grains (continuous film) has an
empty grain list and no grain files; the viewer then steps through its snapshots. Pixels with no magnetic atom (grain boundaries) hold -128 in the int8
arrays. The first array index of `[3, ...]` arrays is the x, y, z component. An atom at x goes to column
`rint(x / a) // f` (`a`: `atom_pitch_nm`, `f` atoms per pixel, so `pitch_nm = f a`), and likewise in y: column `ix`
spans `ix pitch - a/2` to `(ix + 1) pitch - a/2`, centred on its atoms, and row `iy` counts from the bottom of the film
(the viewer flips it for display). Dipole cell `i` spans `i cell` to `(i + 1) cell`.

A series with an overview is read in blocks by HTTP byte ranges (`hamrvis view` answers them); no other part of the
format depends on its size. While `reduce` runs, a folder `.reduce-work` holds its per-atom indices and progress; a bundle
that still has it is unfinished.

## manifest.json

```
format: "hamrvis-bundle", version: 2, created, name
verify               written by hamrvis verify: {hamrvis_version, date, checks, failed, runs (null: all), every}
runs: [ {                 in run order (the order the runs started), one per tab
   id, title            folder name (unique); tab title from the run type
   type                 hamr | monte-carlo | dipole-read | generic (hamrvis/runs.py TYPES)
   ref                  id of the nearest hamr run before this one (whose bits it is compared with), if any
   monte_carlo          true for Monte Carlo integrators (time axis in MC steps)
   source, segments     run directory; all directories of a run continued over several jobs, in order
   output               {columns, rows, dtype: "f64"}
   input                {key: value} of the input file (flags: ""), input_units {key: unit}
   dt                   time step (s); size_nm [x, y, z]; magnetic_types (0-based material indices); n_atoms
   grains               {id, x, y (nm), atoms, magnetic_atoms}; x, y, magnetic_atoms from the atoms when the run
                        wrote coordinates, else from grain-coordinates.txt plus outline: [[x, y] ...] (Voronoi, approximate)
   grain_map            {nx, ny, pitch_nm, atom_pitch_nm, mixed_columns, from_run (if copied from another run)}
   checks               {column_rounding_max, grain_txt_same_step {n, max_diff},
                         grain_txt_one_step_earlier {n, max_diff, median_diff}}
   source_files         {name: [bytes, mtime]} of input, output, grain.txt, dipole-field.txt
   grain_series         {frames, lead: ["time", "temperature"], dtype: "f32" (if not int16), every: N (if thinned),
                         overview: {bins} (if grains_overview.i16 exists; frames f of bin b: b*frames//bins <= f < (b+1)*frames//bins)}
   snapshots            {grid: {pitch_nm, pitch_z_nm, z0_nm, nx, ny, nz, vx, vy}, volume, frames: [{index, time}],
                         skipped: [indices], dtype: "i16" | "f32" (if not int8)}
   notes                [text] what was left out and why (no grain-coordinates.txt, unreadable grain.txt)
   hamr (hamr runs)         head model parameters (see docs/analysis.md), including Hpeak, field_uz, rise, fall
   Tc, Tc_source, freeze_lag_nm, analysis {grains: [...], thresholds} (hamr runs)
   dipole               {nx, ny, nz, cell_nm, cell_nm_xyz [x, y, z] (if the cells are not cubes)}
                        (runs that wrote dipole-field.txt)
} ]
```

## Grids

- Top-down columns use the smallest spacing between atom x coordinates (2.5 A in the rocksalt-supercell FePt
  model, so 800 x 240 columns for 200 x 60 nm). Periodic coordinates wrap.
- The per-layer volume uses twice that pitch in-plane (one site per layer per pixel can be empty in non-cubic
  lattices; 2 x 2 sites never are inside a grain) and one layer per distinct atom z.
- Dipole cells are gridded by cell ID: VAMPIRE numbers macrocells x-major with z fastest on
  (L / cell size + 1) cells per axis, the cell size per axis from cells:macro-cell-size or
  cells:macro-cell-size-x/y/z (the last line wins, as in VAMPIRE; default 10 A). The x, y, z columns of dipole-field.txt are atom centroids, not
  cell centres, so they are not used for gridding.

## Changing the format

Bump `FORMAT_VERSION` in `hamrvis/reduce.py` for incompatible changes and keep the viewer able to read older
bundles (or say in the viewer that the bundle must be reduced again).

## Version 1

Bundles reduced before the run list (version 1) hold `stages: {relax, write, read}` instead of `runs`, with each
stage's files in `relax/`, `write/`, `read/` and `output.f32` (float32). The viewer opens them as three tabs with
the same titles and links (`#relax`, `#write`, `#read`); `read` is compared with `write`.
