#!/bin/bash
# Regenerate tests/fixture-run and tests/fixture-leadin with a serial VAMPIRE build:
# tests/make_fixture.sh /path/to/vampire-serial
set -eu
cd "$(dirname "$0")/fixture-run"
rm -f atoms-coords* spins-* grain* output log stdout
"$1" > stdout
# fixture-leadin keeps only the output and log: its lead-in is longer than a bit, so VAMPIRE writes fewer bits than it reports
cd ../fixture-leadin
"$1" > /dev/null
rm -f atoms-coords* spins-* grain*
