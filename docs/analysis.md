# Write analysis and head model

## Head, laser and field (`hamrvis/vampire.py` `Hamr`, ported to the viewer as `Model` in `hamrvis/viewer/app.js`)

From `src/hamr/hamr_continuous.cpp` and `src/hamr/field_profiles.cpp`. Steps count from the start of the run's
`sim::time` (equilibration included); `k = step - eq`.

- Steps per bit `nb = round(BL / v / dt)`; run-in and run-out `ne = round(BL / 2 / v / dt) + round(NPS / v / dt)`;
  bits per track `floor((Lx - 0.1 nm) / BL)`; steps per track `nt = bits_per_track nb + 2 ne`.
- Bits actually written per track `w = ceil((bits_per_track nb - ne) / nb)` (the bit loop runs from `ne` to
  `bits_per_track nb`), so `w < bits_per_track` when `ne >= nb`; bit `b < w` of track `t` takes sequence entry
  `t w + b`, cells `b >= w` get no field.
- Laser on for `1 <= k <= n_tracks nt`. Track `t = (k - 1) // nt`, step in track `j = (k - 1) % nt + 1`.
- Laser centre `x = -BL/2 - NPS + (j + 1) v dt`, `y = TW (floor(t w / bits_per_track) + 0.5) + track_padding`.
- Laser temperature `T = Tmin + (Tmax - Tmin) exp(-dx^2 / 2 sx^2 - dy^2 / 2 sy^2)`, `s = FWHM / sqrt(8 ln 2)`.
- Head field: inside bit `b` of a track (0-based step `i` in the bit) the magnitude rises linearly over `rise`
  steps and falls over the last `fall` steps; it is applied in the box `|x - (x_laser - NPS)| <= field_x / 2`,
  `|y - y_laser| <= field_y / 2`, along `sim:applied-field-unit-vector`, with the sign of the bit.
- **Peak field `Hpeak = Hmax - Hmin`**: the trapezoid adds `(Hmax - Hmin) / rise` per step. VAMPIRE's defaults are
  `sim:maximum-applied-field-strength = +1 T` and `sim:minimum-applied-field-strength = -1 T`, so an input that
  sets neither (only `sim:applied-field-strength`) writes at **2 T**. This was confirmed with a small serial run's
  `B_applied` column; `tests/test_hamrvis.py` checks the model against VAMPIRE's output to 1e-12 T.

## Grain outcomes (`hamrvis/reduce.py` `write_analysis`)

Per grain, from `grain.txt` (z-magnetisation `m_z |m|`):

| Quantity | Definition |
|---|---|
| melted | its \|m\| dropped below 0.1 **and** it lost its final direction (`m_z |m|` below 0.05 in the final direction) at some point. Grains at the edge of the hot band can dip below 0.1 in \|m\| without losing their direction. |
| decision moment | the frame after the last one where `m_z |m|` was below 0.05 in its final direction: after it the grain never changes sign again. |
| undecided | melted, but `|final m_z |m||` < 0.3: still disordered at the end. |
| circumstance | head field at the decision: off (`zero-field`), below 0.4 Hpeak (`reversal`, during a ramp), or full (`full-field`). |
| freezing temperature Tc | median laser temperature at the grain centre at the decision moment of all decided melted grains (or `--tc`). It is a measured, effective value for the head speed, not the equilibrium Curie temperature. |
| freeze lag | distance behind the laser centre where the track centre cools through Tc: `sx sqrt(2 ln((Tmax - Tmin) / (Tc - Tmin)))`. |
| intended direction | the field direction of the bit the grain lies in, with bit `b` covering `[b BL - lag, (b + 1) BL - lag)` on the nearest track: bits are written a lag behind the laser. The sign is `bit polarity x z-component of the field unit vector`. |
| written as intended | melted, decided, and its final sign equals the intended direction. |

Known limits (from the earlier earlier runs, still true): a grain's centre temperature is a rough proxy for when it
froze (the Tc front needs ~0.4 ns to cross a grain at 20 m/s); grains near the track edge freeze closer to the
laser (the Tc isotherm is an ellipse) and can take the previous bit's field; grains deciding during a reversal can
go either way.
