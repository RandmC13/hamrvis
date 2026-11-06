# hamrvis

Viewer for VAMPIRE heat-assisted magnetic recording (HAMR) simulations. Give it a folder of runs (for example a
Monte Carlo relaxation, one or more HAMR writes, dipole-field reads, in any number and order) and it builds one page
with a tab for each run, showing what that kind of simulation calls for, on one down-track axis.

It has two parts:

1. **`hamrvis reduce`** turns finished run directories into a compact *bundle*. It runs where the data is (on the
   cluster), so only the bundle is copied home. The 90 GB of spin snapshots of a 8.7M-atom write become
   about 2.2 GB, or about 180 MB without the per-layer data (estimated from the final snapshot).
2. **`hamrvis view`** serves one general viewer for any bundle on `localhost` and opens it in the browser. The
   same viewer works for every run, so improvements apply to old runs too.

Only Python 3.8+ with numpy is needed; the viewer is plain HTML and JavaScript with no build step.

## Install

```bash
pip install git+https://github.com/RandmC13/hamrvis.git
```

This gives the `hamrvis` command; `python3 -m hamrvis` from a clone without installing is the same thing. Python 3.8+ and numpy only.

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

To share a run as one file:

```bash
python3 -m hamrvis export bundles/<name> -o <name>.html --grain-step 8
```

The [user guide](docs/user-guide.md) has every command and option, what to put in your VAMPIRE input, what each view shows, and troubleshooting.

## More

- `docs/user-guide.md`: how to use it.
- `docs/decisions.md`: why it is built this way and what was left out.
- `docs/bundle-format.md`: the bundle files.
- `docs/analysis.md`: how the write analysis classifies grains, and the head model.
- `tests/make_large.py`: a generated run of any size, in VAMPIRE's file formats with a known answer, for scale tests.


## Licence and citation

BSD 2-Clause (`LICENSE`). To cite hamrvis, see `CITATION.cff`. Changes are listed in `CHANGELOG.md`.
