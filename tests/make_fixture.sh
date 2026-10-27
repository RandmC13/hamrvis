#!/bin/bash
# Regenerate tests/fixture-run with a serial VAMPIRE build: tests/make_fixture.sh /path/to/vampire-serial
set -eu
cd "$(dirname "$0")/fixture-run"
rm -f atoms-coords* spins-* grain* output log stdout
"$1" > stdout
