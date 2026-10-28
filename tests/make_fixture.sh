#!/bin/bash
# Regenerate tests/fixture-run and tests/fixture-leadin-fixed with a serial VAMPIRE build:
# tests/make_fixture.sh /path/to/vampire-serial
# tests/fixture-leadin was made by VAMPIRE before its bit loop was fixed (it lost bits to a long lead-in) and is kept
# as recorded, so that runs from those versions stay covered; do not regenerate it.
set -eu
cd "$(dirname "$0")/fixture-run"
rm -f atoms-coords* spins-* grain* output log stdout
"$1" > stdout
# fixture-leadin-fixed keeps only the output and log: a lead-in (half a bit + NPS) longer than a bit, two tracks
cd ../fixture-leadin-fixed
"$1" > /dev/null
rm -f atoms-coords* spins-* grain*
