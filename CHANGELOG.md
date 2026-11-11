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
- What is recorded, reconstructed and interpreted is marked on every readout. The freeze lag can be changed in the
  viewer (inferred, none, custom, from bit centres), every lag-dependent view follows it and states it, and a chart
  shows the score against the lag. Each write's head model is checked against VAMPIRE's `B_applied` and logged head
  positions when it is reduced. Bits reported and bits written are shown with the source of the count; fallbacks are
  notes on the run; bundles reduced before the current analysis model open with a warning.

Fixed:

- A write that does not set `sim:maximum-temperature` is now modelled at VAMPIRE's default of 300 K, not 0 K.
- `dimensions:system-size` and VAMPIRE's default size of 10 nm are now understood, instead of stopping the reduce.
- A write without `hamr:bit-size` or `hamr:track-size`, or an output table with no rows, now gets a clear message.
- Bit cells no longer wrap round the edge of the film: a cell before the start of the film held grains from its far
  end, so a bit that lies off the film was counted (the large run's "bits matching intended" read 30 / 33, not 29 / 32).
- Atom grids are drawn centred on their atoms (they were half an atom spacing off), and every grid and dipole cell is
  placed by its own pitch, so the film, the hover readouts and the bit cells agree.
- The inspector gives the head field at the decision (which sets the category) and the field at the grain, which is 0
  outside the field box; it gave only the first, labelled as the second.
- The space-time chart follows the freeze lag; the dipole x–z lane marks the top of the magnetic atoms rather than a fixed
  10 nm; times under 1 ns show 0.01 ps; the x–z section gives the heights of its atom layers.
- The local server can no longer be led outside the bundle by a Windows drive path.
- `export` no longer depends on the system's text encoding, which stopped it on Windows with Python before 3.15.
