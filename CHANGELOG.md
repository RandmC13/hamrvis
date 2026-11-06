# Changelog

## 0.9.0 (unreleased)

First packaged version.

- One tab per run found, whatever the workflow: writes, Monte Carlo runs, reads and any other VAMPIRE run, with runs
  continued over several jobs joined into one tab.
- Input files read as VAMPIRE reads them, including its units and integer lists; `hamrvis inspect`; options for
  precision, thinning and thresholds.
- Reduction streams its input, so memory does not grow with the size of the run, and an interrupted reduce carries on
  when run again. The viewer handles films and series too large to load at once.
- Reads the window snapshots of `hamr:active-window` runs.
- `pip install` with a `hamrvis` command; the version is stored in each bundle and shown under Run details.

Fixed:

- A write that does not set `sim:maximum-temperature` is now modelled at VAMPIRE's default of 300 K, not 0 K.
- `dimensions:system-size` and VAMPIRE's default size of 10 nm are now understood, instead of stopping the reduce.
- A write without `hamr:bit-size` or `hamr:track-size`, or an output table with no rows, now gets a clear message.
- The local server can no longer be led outside the bundle by a Windows drive path.
