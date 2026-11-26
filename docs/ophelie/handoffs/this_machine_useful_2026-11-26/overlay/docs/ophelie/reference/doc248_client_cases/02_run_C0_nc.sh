#!/usr/bin/env bash
# Doc248 C0 — client NC (flat slip + EM fixed-power).
# Power ledger:
#   generator_input = 30 kW (client-confirmed supply)
#   EM absorbed scaling target default = 18 kW (engineering ASSUMED; 17–20 kW OK via TARGET_ABSORBED)
#   NOT French Jacoutot 50 kW; NOT a claim that η·30 kW is physically computed.
# Requires 01_generate_particles first.
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C0_nc}
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload}
END_TIME=${END_TIME:-60}
WALL_S=${WALL_S:-36000}   # 10 h wall cap via timeout
# Absorbed EM scaling target [W]. Engineering ASSUMED from generator 30 kW ledger.
# Default 18000 (was 19000); override e.g. TARGET_ABSORBED=17000.
TARGET_ABSORBED=${TARGET_ABSORBED:-18000}
# Robin cooling: French default h_side=300 → P_loss≈57 kW @1000°C (net≪0).
# 0D retune for P_loss≈P_absorbed≈19 kW: h_side≈50 (bottom/free keep French defaults).
# Override e.g. H_SIDE=80 or H_SIDE=30.
H_SIDE=${H_SIDE:-50}
H_BOTTOM=${H_BOTTOM:-35}
H_FREE=${H_FREE:-20}
# Free-surface radiation (default French 0.8). Top-layer over-cool → lower EMISSIVITY first
# (rad ≫ conv on free surface); optionally also lower H_FREE.
EMISSIVITY=${EMISSIVITY:-0.8}
# Optional IC CSV. Default empty = isothermal T0.
# Hybrid thin tongue + axis layers:
#   INIT_T_RZ=$PACK/configs/client_nc_hybrid_tongue_H230_ASSUMED_T_rz.csv
# Pure Fig.8 form: INIT_T_RZ=$PACK/configs/client_nc_fig8_scaled_H230_ASSUMED_T_rz.csv
INIT_T_RZ=${INIT_T_RZ:-}
MELT_HEIGHT=${MELT_HEIGHT:-0.230}
# Joule-heat shape levers (Fig.6 tongue is upper-outer; current CAD remap is 15–92%H):
#   Q_UPPER_KAPPA>0 → Doc246 controlled reweight Q∝exp[κ(z/H−0.5)] (works on existing binary)
#   COIL_Z_FRAC_LO/HI → compress filament band (needs rebuild with --client-nc-coil-z-frac-*)
Q_UPPER_KAPPA=${Q_UPPER_KAPPA:-0}
COIL_Z_FRAC_LO=${COIL_Z_FRAC_LO:-}
COIL_Z_FRAC_HI=${COIL_Z_FRAC_HI:-}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q" ]] || { echo "FATAL: compile frozen_q"; exit 1; }

# Resolve relative paths before cd into FROZEN bin.
mkdir -p "$OUT/out"
OUT=$(cd "$OUT" && pwd)
RELOAD=$(cd "$RELOAD" && pwd)

INIT_T_RZ_ARGS=()
if [[ -n "$INIT_T_RZ" ]]; then
  [[ -f "$INIT_T_RZ" ]] || { echo "FATAL: INIT_T_RZ not found: $INIT_T_RZ"; exit 1; }
  INIT_T_RZ=$(cd "$(dirname "$INIT_T_RZ")" && pwd)/$(basename "$INIT_T_RZ")
  INIT_T_RZ_ARGS=(--init-T-rz-csv="$INIT_T_RZ")
fi
Q_SHAPE_ARGS=()
if [[ -n "$Q_UPPER_KAPPA" && "$Q_UPPER_KAPPA" != "0" ]]; then
  Q_SHAPE_ARGS+=(--doc246-q-upper-kappa="$Q_UPPER_KAPPA")
fi
if [[ -n "$COIL_Z_FRAC_LO" ]]; then
  Q_SHAPE_ARGS+=(--client-nc-coil-z-frac-lo="$COIL_Z_FRAC_LO")
fi
if [[ -n "$COIL_Z_FRAC_HI" ]]; then
  Q_SHAPE_ARGS+=(--client-nc-coil-z-frac-hi="$COIL_Z_FRAC_HI")
fi

cd "$FROZEN"
rm -rf output output_backup && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C0_effective.json" "$OUT/"
if [[ -n "$INIT_T_RZ" ]]; then
  cp -f "$INIT_T_RZ" "$OUT/" 2>/dev/null || true
  [[ -f "${INIT_T_RZ%.csv}.meta.json" ]] && cp -f "${INIT_T_RZ%.csv}.meta.json" "$OUT/" 2>/dev/null || true
fi

echo "[doc248][C0] start $(date -Is) end=$END_TIME wall=${WALL_S}s"
echo "[doc248][C0] melt_height=$MELT_HEIGHT m (NC fill; coil remapped in binary)"
echo "[doc248][C0] generator_ledger=30000 W; EM_absorbed_target=$TARGET_ABSORBED W (17–20 kW ASSUMED band; not French 50 kW)"
echo "[doc248][C0] Robin cooling: h_side=$H_SIDE h_bottom=$H_BOTTOM h_free=$H_FREE emissivity=$EMISSIVITY (aim P_loss≈P_absorbed near operating T)"
if [[ -n "$INIT_T_RZ" ]]; then
  echo "[doc248][C0] init_T_rz=$INIT_T_RZ (stratified ASSUMED; thaw path, not freeze-T)"
else
  echo "[doc248][C0] init_T=isothermal client T0 (~1000C)"
fi
if [[ -n "$Q_UPPER_KAPPA" && "$Q_UPPER_KAPPA" != "0" ]]; then
  echo "[doc248][C0] Q_upper_kappa=$Q_UPPER_KAPPA (Doc246 controlled shape; NOT paper Maxwell)"
fi
if [[ -n "$COIL_Z_FRAC_LO" || -n "$COIL_Z_FRAC_HI" ]]; then
  echo "[doc248][C0] coil_z_frac=[${COIL_Z_FRAC_LO:-0.15},${COIL_Z_FRAC_HI:-0.92}] (needs rebuilt binary)"
fi
set +e
timeout --signal=TERM --kill-after=120 "${WALL_S}s" \
  ./test_3d_ophelie_french_natural_convection_frozen_q \
    --client-nc-geometry --client-nc-material \
    --dp=0.010 --melt-height="$MELT_HEIGHT" \
    --reload-dir="$RELOAD" \
    --doc223-paper-flat-slip-top \
    --em-control=fixed-power --target-absorbed-power="$TARGET_ABSORBED" \
    --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
    --enable-thermal-diffusion \
    --h-side="$H_SIDE" --h-bottom="$H_BOTTOM" --h-free="$H_FREE" --emissivity="$EMISSIVITY" \
    "${INIT_T_RZ_ARGS[@]}" \
    "${Q_SHAPE_ARGS[@]}" \
    --doc224-u-blowup-gate=50 \
    --end-time="$END_TIME" --hydro-settle-time=1.0 \
    --state-record-interval=5.0 \
    --meridional-average-start=1.0 \
    --boundary-face-t-mode=local_dt_dz \
    --probes-timeseries="$OUT/probes_timeseries.csv" \
    --summary-json="$OUT/summary.json" \
    2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C0] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"

cp -f output/french_energy_budget.csv output/french_thermo_em_coupling.csv \
      output/french_spatial_stats.csv output/french_T_rz.csv \
      output/french_meridional_timeseries.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
exit "$rc"
