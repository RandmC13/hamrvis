# Bundle format (version 1)

A bundle is a directory written by `hamrvis reduce`. All binary files are raw little-endian arrays with no header;
their shapes come from `manifest.json`. Files ending in `.gz` are gzip-compressed (the viewer inflates them with the
browser's `DecompressionStream`).

```
manifest.json
<stage>/output.f32            float32 [rows, cols]          the run's 'output' file; column names in the manifest
<stage>/grains.i16            int16   [frames, grains, 4]   grain.txt per-grain block: mx, my, mz (direction), |m|; value x 32767
<stage>/grain_lead.f64        float64 [frames, n_lead]      grain.txt leading columns (time in s, temperature in K)
<stage>/top/NNNNNNNN.i8.gz    int8    [3, ny, nx]           snapshot NNNNNNNN: mean spin of the magnetic atoms in each column, x 127
<stage>/vol/NNNNNNNN.i8.gz    int8    [3, nz, vy, vx]       snapshot NNNNNNNN: mean spin per layer on 2x2-site pixels, x 127
<stage>/top_count.u8          uint8   [ny, nx]              magnetic atoms per column (capped at 255)
<stage>/vol_count.u8          uint8   [nz, vy, vx]          magnetic atoms per volume pixel
<stage>/min_m.u8              uint8   [ny, nx]              lowest column |mean spin| over all snapshots, x 255
<stage>/grain_map.i16         int16   [ny, nx]              grain index of the magnetic atoms in each column (-1: none)
read/dipole.f32               float32 [3, nz, ny, nx]       dipole-field.txt cell fields Bx, By, Bz (T)
read/dipole_atoms.u16         uint16  [nz, ny, nx]          atoms per dipole cell (0: a field-probe cell in non-magnetic material)
```

`<stage>` is `relax`, `write` or `read`. Pixels with no magnetic atom (grain boundaries) hold -128 in the int8
arrays. The first array index of `[3, ...]` arrays is the x, y, z component. Row `iy` of a map is y = (iy + 0.5) x
pitch from the bottom of the film; the viewer flips it for display.

## manifest.json

```
format: "hamrvis-bundle", version: 1, created, name
stages: { <stage>: {
   source               absolute path of the run directory that was reduced
   input                {key: value} of the input file (flags: ""), input_units {key: unit}
   dt                   time step (s); size_nm [x, y, z]; magnetic_types (0-based material indices); n_atoms
   output               {columns, rows}
   grains               {id, x, y (nm), atoms, magnetic_atoms}; x, y, magnetic_atoms from the atoms when the run
                        wrote coordinates, else from grain-coordinates.txt plus outline: [[x, y] ...] (Voronoi, approximate)
   grain_map            {nx, ny, pitch_nm, mixed_columns, from_stage (if copied from another stage)}
   checks               {column_rounding_max, grain_txt_same_step {n, max_diff},
                         grain_txt_one_step_earlier {n, max_diff, median_diff}}
   source_files         {name: [bytes, mtime]} of input, output, grain.txt, dipole-field.txt
   grain_series         {frames, lead: ["time", "temperature"]}
   snapshots            {grid: {pitch_nm, pitch_z_nm, z0_nm, nx, ny, nz, vx, vy}, volume, frames: [{index, time}]}
   hamr (write)         head model parameters (see docs/analysis.md), including Hpeak, field_uz, rise, fall
   Tc, Tc_source, freeze_lag_nm, analysis {grains: [...], thresholds} (write)
   dipole (read)        {nx, ny, nz, cell_nm}
}}
```

## Grids

- Top-down columns use the smallest spacing between atom x coordinates (2.5 A in the rocksalt-supercell FePt
  model, so 800 x 240 columns for 200 x 60 nm). Periodic coordinates wrap.
- The per-layer volume uses twice that pitch in-plane (one site per layer per pixel can be empty in non-cubic
  lattices; 2 x 2 sites never are inside a grain) and one layer per distinct atom z.
- Dipole cells are gridded by cell ID: VAMPIRE numbers macrocells x-major with z fastest on
  (L / cells:macro-cell-size + 1) cells per axis. The x, y, z columns of dipole-field.txt are atom centroids, not
  cell centres, so they are not used for gridding.

## Changing the format

Bump `FORMAT_VERSION` in `hamrvis/reduce.py` for incompatible changes and keep the viewer able to read older
bundles (or say in the viewer that the bundle must be reduced again).
