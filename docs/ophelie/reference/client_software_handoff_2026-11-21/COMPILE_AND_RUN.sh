#!/usr/bin/env bash
# Client NC + stirring handoff: compile and run recipes (Linux SYCL GPU).
# Usage: source oneAPI setvars first if needed, then bash COMPILE_AND_RUN.sh <step>
# Steps: audit | nc_relax | nc_iso | nc_dT | stir_relax | e0 | e1 | e2 | f0 | f1 | all_build
set -euo pipefail
SRC="${SRC:-/home/ubuntu/SPHinXsys}"
BUILD="${BUILD:-$SRC/build}"
OUT_NC="$SRC/docs/ophelie/reference/client_software_handoff_2026-11-21/nc"
OUT_STIR="$SRC/docs/ophelie/reference/elastic_paddle_runs_2026-11-21"
STEP="${1:-}"

cd "$BUILD"

case "$STEP" in
  all_build)
    cmake .
    cmake --build . --target \
      test_3d_ophelie_client_case_audit \
      test_3d_ophelie_client_nc_glass_relax \
      test_3d_ophelie_french_natural_convection_frozen_q \
      test_3d_ophelie_client_stirring_glass_relax \
      test_3d_ophelie_doc228_e0_paddle_hold \
      test_3d_ophelie_doc228_e1_paddle_load \
      test_3d_ophelie_doc228_e2_paddle_spin \
      test_3d_ophelie_client_stirring_f0_elastic_fsi \
      test_3d_ophelie_client_stirring_f1_heated_elastic \
      -j"$(nproc)"
    ;;
  audit)
    cmake --build . --target test_3d_ophelie_client_case_audit -j"$(nproc)"
    ./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_case_audit/bin/test_3d_ophelie_client_case_audit \
      --repo-root "$SRC" \
      --output-dir "$SRC/docs/ophelie/reference/client_case_audit_2026-11-20"
    ;;
  nc_relax)
    cmake --build . --target test_3d_ophelie_client_nc_glass_relax -j"$(nproc)"
    cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin
    ./test_3d_ophelie_client_nc_glass_relax \
      --dp=0.010 --melt-height=0.180 --relax-steps=1000 --relax-vtp-every=100 \
      --output-meta-dir=./client_nc_relax_meta
    ;;
  nc_iso)
    cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j"$(nproc)"
    FROZEN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin"
    RELOAD="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload"
    mkdir -p "$OUT_NC"
    cd "$FROZEN"
    ./test_3d_ophelie_french_natural_convection_frozen_q \
      --client-nc-geometry --client-nc-material --zero-q --zero-buoyancy \
      --dp=0.010 --melt-height=0.180 --reload-dir="$RELOAD" \
      --end-time=0.50 --hydro-settle-time=0.15 --state-record-interval=0.05 \
      --summary-json="$OUT_NC/client_nc_isothermal_hydro_summary.json"
    ;;
  nc_dT)
    cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j"$(nproc)"
    FROZEN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin"
    RELOAD="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload"
    mkdir -p "$OUT_NC"
    cd "$FROZEN"
    ./test_3d_ophelie_french_natural_convection_frozen_q \
      --client-nc-geometry --client-nc-material --zero-q \
      --dp=0.010 --melt-height=0.180 --reload-dir="$RELOAD" \
      --prescribe-t-radial=20 --end-time=0.10 --hydro-settle-time=0.05 \
      --summary-json="$OUT_NC/client_nc_dT20_boussinesq_summary.json"
    ./test_3d_ophelie_french_natural_convection_frozen_q \
      --client-nc-geometry --client-nc-material --zero-q --zero-buoyancy \
      --dp=0.010 --melt-height=0.180 --reload-dir="$RELOAD" \
      --prescribe-t-radial=20 --end-time=0.10 --hydro-settle-time=0.05 \
      --summary-json="$OUT_NC/client_nc_dT20_zero_buoyancy_summary.json"
    ;;
  stir_relax)
    cmake --build . --target test_3d_ophelie_client_stirring_glass_relax -j"$(nproc)"
    cd tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin
    ./test_3d_ophelie_client_stirring_glass_relax --bodies=all
    ;;
  e0)
    cmake --build . --target test_3d_ophelie_doc228_e0_paddle_hold -j"$(nproc)"
    BIN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e0_paddle_hold/bin"
    mkdir -p "$OUT_STIR/E0_hold"
    cd "$BIN"
    ./test_3d_ophelie_doc228_e0_paddle_hold \
      --dp=0.004 --solid-steps=100 --clamp-top-band=0.05 \
      --out-dir="$OUT_STIR/E0_hold" --summary-json=summary.json
    ;;
  e1)
    cmake --build . --target test_3d_ophelie_doc228_e1_paddle_load -j"$(nproc)"
    BIN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e1_paddle_load/bin"
    mkdir -p "$OUT_STIR/E1_load_response"
    cd "$BIN"
    ./test_3d_ophelie_doc228_e1_paddle_load \
      --dp=0.004 --solid-steps=5000 --clamp-top-band=0.05 \
      --force1=1.5e4 --force2=3e4 \
      --out-dir="$OUT_STIR/E1_load_response" --summary-json=summary.json
    ;;
  e2)
    cmake --build . --target test_3d_ophelie_doc228_e2_paddle_spin -j"$(nproc)"
    BIN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_doc228_e2_paddle_spin/bin"
    mkdir -p "$OUT_STIR/E2_driven_rotation"
    cd "$BIN"
    ./test_3d_ophelie_doc228_e2_paddle_spin \
      --dp=0.004 --end-time=0.05 --omega-ramp=0.05 --clamp-top-band=0.05 \
      --max-wall-hours=2 \
      --out-dir="$OUT_STIR/E2_driven_rotation" --summary-json=summary.json
    ;;
  f0)
    cmake --build . --target test_3d_ophelie_client_stirring_f0_elastic_fsi -j"$(nproc)"
    RELAX_RELOAD="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload"
    BIN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f0_elastic_fsi/bin"
    mkdir -p "$OUT_STIR/F0_isothermal_fsi" "$BIN/output"
    cd "$BIN"
    ./test_3d_ophelie_client_stirring_f0_elastic_fsi \
      --reload-dir="$RELAX_RELOAD" \
      --rotor-stl=./data/client_stir_paddle_solver_m.stl \
      --end-time=0.05 --omega-ramp=0.05 --c0=8 --clamp-top-band=0.05 \
      --max-wall-hours=2 \
      --out-dir="$OUT_STIR/F0_isothermal_fsi" --summary-json=summary.json
    ;;
  f1)
    cmake --build . --target test_3d_ophelie_client_stirring_f1_heated_elastic -j"$(nproc)"
    RELAX_RELOAD="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload"
    BIN="$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin"
    mkdir -p "$OUT_STIR/F1_heated_integration" "$BIN/output"
    cd "$BIN"
    ./test_3d_ophelie_client_stirring_f1_heated_elastic \
      --reload-dir="$RELAX_RELOAD" \
      --rotor-stl=./data/client_stir_paddle_solver_m.stl \
      --end-time=0.2 --omega-ramp=0.05 --c0=8 \
      --heat-power=1e5 --clamp-top-band=0.05 --max-wall-hours=2 \
      --out-dir="$OUT_STIR/F1_heated_integration" --summary-json=summary.json
    ;;
  *)
    echo "Usage: $0 {all_build|audit|nc_relax|nc_iso|nc_dT|stir_relax|e0|e1|e2|f0|f1}"
    exit 1
    ;;
esac
