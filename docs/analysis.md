# Write analysis and head model

## Head, laser and field (`hamrvis/vampire.py` `Hamr`, ported to the viewer as `Model` in `hamrvis/viewer/app.js`)

From `src/hamr/hamr_continuous.cpp` and `src/hamr/field_profiles.cpp`. Steps count from the start of the run's
`sim::time` (equilibration included); `k = step - eq`.

- Steps per bit `nb = round(BL / v / dt)`; run-in and run-out `ne = round(BL / 2 / v / dt) + round(NPS / v / dt)`;
  bits per track `floor((Lx - 0.1 nm) / BL)`; steps per track `nt = bits_per_track nb + 2 ne`.
- Bits actually written per track `w`: the distinct "Written bit N" lines of the log before the first end of track
  (`written_source = "log"`). Without a whole track in the log, `w = ceil((bits_per_track nb - ne) / nb)`, VAMPIRE's
  bit loop before it was fixed (it ran from `ne` to `bits_per_track nb`, so `w < bits_per_track` when `ne >= nb`;
  `written_source = "formula"`, and the run gets a note). Bit `b < w` of track `t` takes sequence entry `t w + b`,
  cells `b >= w` get no field. The viewer shows both counts, VAMPIRE's "Bits to be written" and the bits written.
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

## Head check at reduce time (`Hamr.check`, stored as `head_check`)

Every write is checked against what VAMPIRE recorded for it, and the result is shown under Run details:
`B_applied` of the output table at every row whose printed time resolves a single step, and the head position
VAMPIRE logs after each bit ("Head moving to position"; in a job continued from a checkpoint only the lines after
"continued from checkpoint", as the replayed loop logs a stale head). VAMPIRE prints 6 significant digits, so values
agree when within a relative 1e-5. A run with neither is "not verifiable"; a difference is also a note on the run.

## Grain outcomes (`hamrvis/reduce.py` `write_analysis`)

Per grain, from `grain.txt` (z-magnetisation `m_z |m|`):

| Quantity | Definition |
|---|---|
| melted | its \|m\| dropped below 0.1 **and** it lost its final direction (`m_z |m|` below 0.05 in the final direction) at some point. Grains at the edge of the hot band can dip below 0.1 in \|m\| without losing their direction. |
| decision moment | the frame after the last one where `m_z |m|` was below 0.05 in its final direction: after it the grain never changes sign again. |
| undecided | melted, but `|final m_z |m||` < 0.3: still disordered at the end. |
| no magnetic atoms | a grain whose atoms are all non-magnetic (a non-magnetic layer VAMPIRE gives a grain ID): `grain.txt` holds 0 for it, which would read as melted. Taken from the magnetic atom count, or, without coordinate files, from \|m\| = 0 in every frame. Left out of the score, the freezing temperature and the views of melted grains (category `empty`). |
| circumstance | head field at the decision: off (`zero-field`), below 0.4 Hpeak (`reversal`, during a ramp), or full (`full-field`). |
| freezing temperature Tc | median laser temperature at the grain centre at the decision moment of the settled melted grains decided under the laser (`Tc_n` of them), or `--tc` (`Tc_kind` "inferred" or "given"). It is inferred from the run, an effective value for the head speed, not the equilibrium Curie temperature. |
| freeze lag | distance behind the laser centre where the track centre cools through Tc: `sx sqrt(2 ln((Tmax - Tmin) / (Tc - Tmin)))`. |
| intended direction | the field direction of the bit the grain lies in, with bit `b` covering `[b BL - lag, (b + 1) BL - lag)` on the nearest track: bits are written a lag behind the laser. The sign is `bit polarity x z-component of the field unit vector`. |
| matching its intended bit | melted, decided, and its final sign equals the intended direction. The score is the number of these over the melted grains inside a bit cell; melted grains outside every cell are left out and counted. |

## The lag in the viewer

The reducer scores the write at the inferred lag (or the one from `--tc`). The viewer's Write analysis section can
score it again at another lag: none, a custom value, or one taken from bit centre positions (bit `b` is centred at
`(b + 0.5) BL - lag`, so the lag is the mean of `(b + 0.5) BL - x_b` over the centres given). This is viewer state
only and never changes the bundle. Every lag-dependent view follows it (grain outcomes, colours against the intended
bit, bit cells, bit polarity, the bits lane, the read-back's intended transitions, Save PNG) and states the lag in
use and its source. A chart shows the score against the lag with the lag in use marked. The viewer's scoring
(`scoreGrains` in `app.js`) is tested to give the reducer's bit, intended direction and outcome for every grain at
the bundle's own lag.

Readouts are marked `r` (recorded by VAMPIRE: magnetisation, |m|, column values), `c` (reconstructed from the
input with the head model: laser temperature, head field, positions relative to the laser) or `i` (interpreted
by hamrvis: decision, intended bit, outcomes, freezing temperature, lag).

A bundle reduced before the current analysis model (hamrvis 0.9.1; bundles without `hamrvis_version` are older)
opens with a warning that its write analysis may differ from what this version computes.

Known limits: one lag is used across the whole track (see the next point); a grain's centre temperature is a rough proxy for when it
froze (the Tc front needs ~0.4 ns to cross a grain at 20 m/s); grains near the track edge freeze closer to the
laser (the Tc isotherm is an ellipse) and can take the previous bit's field; grains deciding during a reversal can
go either way.
