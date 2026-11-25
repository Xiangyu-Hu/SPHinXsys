#!/usr/bin/env bash
# Prepare CPU programmed-temperature bin/ so IOEnvironment won't choke on a
# half-broken output_backup (the failure mode of the first F3 1800 attempt).
set -euo pipefail
BIN="${1:-/home/ubuntu/software/SPHinXsys-SPDB/build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin}"
cd "$BIN"
if [[ -d output ]]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  mv output "output_saved_${stamp}"
  echo "Moved ./output -> ./output_saved_${stamp}"
fi
rm -rf output_backup
mkdir -p output restart reload
echo "Ready in $BIN"
