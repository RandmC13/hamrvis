# Changelog

## 1.0.0 (2026-11-12)

First release, including everything listed below. Checked at full size: an 8.7-million-atom relax, write (with and
without the active window) and read-back reduce to bundles that match the raw files in every stored value, within the
stored precision.

## 0.9.1 (not released)

Added:

- Viewer: *Smooth* in the atom column and layer views (off by default) averages 2 × 2 to 8 × 8 neighbouring pixels
  on screen, so a melted region, whose few atoms per pixel keep a random mean, shows as pale; the data and the
  tooltip are unchanged.
- `hamrvis verify <bundle>` checks every array and number of a bundle against the raw run files with an independent
  reader (formerly `tests/fidelity/bundle_check.py`, which still works), stores the result in the manifest and the
  viewer shows it under Run details.

Changed:

- Viewer: shorter labels in the stage boxes, tooltips and write analysis, without the r/c/i superscripts (the source
  of a readout shows on hovering its label); the explanation of the analysis is one collapsed list with the lag formula
  on its own line; the decision positions plot spans the page and fits the data; the controls stay at the top while
  scrolling, so the sensor height and reader width can be changed next to the field lanes; tooltips no longer
  flicker during playback (in the atom and layer views their column rows disappeared while the next snapshot loaded;
  they now keep the values shown until it arrives, as the film does); Run details gives the number of snapshots, lists
  the checks made when reducing as plain statements of what was compared, and shows a verify result only when there
  is one.

Fixed:

- A run continued in a job on another number of processes (or with another output mode) was split into two tabs,
  and once joined, the later job's snapshots were left out: they were read with the first job's atom coordinates,
  which list the atoms in another order. Each job's snapshots are now read with its own coordinates (checked to
  hold the same atoms), and `hamrvis verify` fails a bundle that left snapshots out.
- A grain with no magnetic atoms (a non-magnetic layer VAMPIRE gives a grain ID, for which `grain.txt` holds 0) is no
  longer counted as a melted, undecided grain: it was part of the score's denominator (the large run read 116 / 143,
  not 116 / 142). Such grains now have the category "no magnetic atoms" and are left out of the write analysis.
  Bundles reduced before this open with the warning that their analysis is older.

## 0.9.0 (not released)

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
