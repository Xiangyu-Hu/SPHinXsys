#!/usr/bin/env bash
# Build a relocatable (best-effort) Linux review bundle. GUI still needs a display to launch.
set -euo pipefail
LINUX_SRC="$(cd "$(dirname "$0")/.." && pwd)"
# Walk up to SPHinXsys repo root (directory with CMakeLists.txt + tests + src)
ROOT="$LINUX_SRC"
while [[ "$ROOT" != "/" ]]; do
  if [[ -f "$ROOT/CMakeLists.txt" && -d "$ROOT/tests" && -d "$ROOT/src" ]]; then
    break
  fi
  ROOT="$(dirname "$ROOT")"
done
VER="${OPHELIE_BUNDLE_VERSION:-0.1.0-dev}"
SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
OUT="${1:-$LINUX_SRC/dist/ophelie-linux-${VER}-${SHA}}"
STAGE="$OUT"
rm -rf "$STAGE"
mkdir -p "$STAGE"/{bin,examples,docs,share}

cp -a "$LINUX_SRC/bin/coldcrucible-solver" "$STAGE/bin/" || true
cp -a "$LINUX_SRC/coldcrucible_solver" "$STAGE/share/"
# Make launcher resolve package-local python module
cat > "$STAGE/bin/coldcrucible-solver" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
export PYTHONPATH="$ROOT/share${PYTHONPATH:+:$PYTHONPATH}"
export OPHELIE_BUILD="${OPHELIE_BUILD:-}"
exec python3 -m coldcrucible_solver "$@"
EOF
chmod +x "$STAGE/bin/coldcrucible-solver"

if [[ -x "$LINUX_SRC/desktop/build/ophelie_desktop" ]]; then
  cp -a "$LINUX_SRC/desktop/build/ophelie_desktop" "$STAGE/bin/"
fi

# Examples (copy, resolve symlinks for STL so bundle is self-contained)
for ex in natural-convection mechanical-stirring; do
  src="$LINUX_SRC/examples/$ex"
  dst="$STAGE/examples/$ex"
  mkdir -p "$dst/assets"
  cp -a "$src/case.json" "$dst/"
  if [[ -d "$src/assets" ]]; then
    for f in "$src/assets"/*; do
      [[ -e "$f" ]] || continue
      cp -aL "$f" "$dst/assets/" 2>/dev/null || cp -a "$f" "$dst/assets/"
    done
  fi
done

# Solver binaries + note: full SYCL runtime bundling is separate; record paths for now
mkdir -p "$STAGE/solver"
for key in nc_relax nc_run stir_relax stir_f0; do
  :
done
# Copy built test binaries if present
B="$ROOT/build/tests/extra_source_and_tests/3d_examples"
copy_bin() {
  local name="$1" rel="$2"
  local p="$B/$rel/bin/$name"
  if [[ -x "$p" ]]; then
    mkdir -p "$STAGE/solver/$name"
    cp -a "$p" "$STAGE/solver/$name/"
    # data dir if any
    if [[ -d "$B/$rel/bin/data" ]]; then
      cp -a "$B/$rel/bin/data" "$STAGE/solver/$name/"
    fi
  fi
}
copy_bin test_3d_ophelie_client_nc_glass_relax test_3d_ophelie_client_nc_glass_relax
copy_bin test_3d_ophelie_french_natural_convection_frozen_q test_3d_ophelie_french_natural_convection_frozen_q
copy_bin test_3d_ophelie_client_stirring_glass_relax test_3d_ophelie_client_stirring_glass_relax
copy_bin test_3d_ophelie_client_stirring_f0_elastic_fsi test_3d_ophelie_client_stirring_f0_elastic_fsi

# Bundle redistributable oneAPI/SYCL libs (NOT libcuda — that comes from NVIDIA driver).
mkdir -p "$STAGE/lib"
ONEAPI_LIB="${ONEAPI_ROOT:-/opt/intel/oneapi}/compiler/2025.2/lib"
copy_lib() {
  local f="$1"
  [[ -e "$f" ]] || return 0
  cp -a "$f" "$STAGE/lib/"
  # also copy soname symlink targets if this is a symlink name only
  if [[ -L "$f" ]]; then
    local real
    real="$(readlink -f "$f")"
    [[ -e "$real" ]] && cp -a "$real" "$STAGE/lib/" || true
  fi
}
if [[ -d "$ONEAPI_LIB" ]]; then
  for f in \
    "$ONEAPI_LIB"/libsycl.so* \
    "$ONEAPI_LIB"/libur_loader.so* \
    "$ONEAPI_LIB"/libur_adapter_cuda.so* \
    "$ONEAPI_LIB"/libimf.so \
    "$ONEAPI_LIB"/libsvml.so \
    "$ONEAPI_LIB"/libirng.so \
    "$ONEAPI_LIB"/libintlc.so*
  do
    copy_lib "$f"
  done
fi
# TBB often required by SYCL runtime
TBB_LIB="$(ls -d /opt/intel/oneapi/tbb/*/lib/intel64/gcc4.8 2>/dev/null | tail -1 || true)"
if [[ -n "${TBB_LIB}" ]]; then
  copy_lib "$TBB_LIB/libtbb.so.12"
  copy_lib "$TBB_LIB/libtbbmalloc.so.2"
fi
# Qt5 + VTK (+ xcb plugin) from packaging list if present
QT_VTK_LIST="$LINUX_SRC/packaging/deps/qt_vtk_libs.txt"
if [[ -f "$QT_VTK_LIST" ]]; then
  mkdir -p "$STAGE/plugins/platforms"
  while IFS= read -r libpath; do
    [[ -z "$libpath" || ! -e "$libpath" ]] && continue
    base="$(basename "$libpath")"
    if [[ "$libpath" == *"/qt5/plugins/"* ]]; then
      cp -a "$libpath" "$STAGE/plugins/platforms/"
    else
      copy_lib "$libpath"
    fi
  done < "$QT_VTK_LIST"
fi

# Record what was shipped
(
  echo "bundled_libs:"
  ls -1 "$STAGE/lib" | sed 's/^/  - /'
  echo "qt_plugins:"
  ls -1 "$STAGE/plugins/platforms" 2>/dev/null | sed 's/^/  - /' || true
  echo "excluded_by_policy:"
  echo "  - libcuda.so.1 (NVIDIA driver on target)"
) > "$STAGE/docs/BUNDLED_LIBS.txt"
cp -a "$LINUX_SRC/packaging/deps/ldd_snapshot.txt" "$STAGE/docs/" 2>/dev/null || true
cp -a "$LINUX_SRC/packaging/deps/RUNTIME_POLICY.txt" "$STAGE/docs/" 2>/dev/null || true
cp -a "$LINUX_SRC/packaging/deps/qt_vtk_libs.txt" "$STAGE/docs/" 2>/dev/null || true

cat > "$STAGE/bin/env.sh" <<EOF
# Source before running. Prefers bundled lib/ over developer oneAPI install.
export OPHELIE_BUNDLE_ROOT="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="\$OPHELIE_BUNDLE_ROOT/bin:\$PATH"
export PYTHONPATH="\$OPHELIE_BUNDLE_ROOT/share\${PYTHONPATH:+:\$PYTHONPATH}"
export LD_LIBRARY_PATH="\$OPHELIE_BUNDLE_ROOT/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
export QT_PLUGIN_PATH="\$OPHELIE_BUNDLE_ROOT/plugins\${QT_PLUGIN_PATH:+:\$QT_PLUGIN_PATH}"
export QT_QPA_PLATFORM_PLUGIN_PATH="\$OPHELIE_BUNDLE_ROOT/plugins/platforms"
# Do NOT override system libcuda; driver must be installed on target.
export ONEAPI_DEVICE_SELECTOR="\${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}"
export OPHELIE_BIN_NC_RELAX="\$OPHELIE_BUNDLE_ROOT/solver/test_3d_ophelie_client_nc_glass_relax/test_3d_ophelie_client_nc_glass_relax"
export OPHELIE_BIN_NC_RUN="\$OPHELIE_BUNDLE_ROOT/solver/test_3d_ophelie_french_natural_convection_frozen_q/test_3d_ophelie_french_natural_convection_frozen_q"
export OPHELIE_BIN_STIR_RELAX="\$OPHELIE_BUNDLE_ROOT/solver/test_3d_ophelie_client_stirring_glass_relax/test_3d_ophelie_client_stirring_glass_relax"
export OPHELIE_BIN_STIR_F0="\$OPHELIE_BUNDLE_ROOT/solver/test_3d_ophelie_client_stirring_f0_elastic_fsi/test_3d_ophelie_client_stirring_f0_elastic_fsi"
export OPHELIE_SOLVER="\$OPHELIE_BUNDLE_ROOT/bin/coldcrucible-solver"
EOF

# GUI launcher wrapper that loads env
cat > "$STAGE/bin/ophelie-desktop.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$HERE/env.sh"
exec "$HERE/ophelie_desktop" "$@"
EOF
chmod +x "$STAGE/bin/ophelie-desktop.sh"

cat > "$STAGE/docs/QUICKSTART_zh.txt" <<EOF
OPHELIE Linux 开发/审查包（$VER / git $SHA）

目录约定
- 本包（安装/解压目录）建议只读：bin、lib、solver、examples 模板
- 项目/结果请放在可写目录（桌面可选结果输出根目录）
- 不要把计算结果写回只读安装位置

CLI（无图形也可）
  source bin/env.sh
  coldcrucible-solver info --json
  coldcrucible-solver selfcheck
  coldcrucible-solver validate --config examples/natural-convection/case.json

桌面 GUI（需要 DISPLAY）
  ./bin/ophelie-desktop.sh

目标机仍需
- 兼容的 NVIDIA 驱动（提供 libcuda.so.1）
- 系统 Qt5/VTK（当前未完整随包；见 docs/RUNTIME_POLICY.txt）

已知限制
- 搅拌默认 elastic_driven；禁止刚体动画冒充
- F1 规定热 != 真实 EM；E2 短窗 != 整圈
- 30 kW 为电源输入
EOF

cp -a "$ROOT/docs/ophelie/software/ENVIRONMENT_AND_BASELINE.md" "$STAGE/docs/" 2>/dev/null || true
cp -a "$ROOT/docs/ophelie/software/230_OPHELIE_LINUX_DESKTOP_THIRD_PARTY_REVIEW_CURSOR_PLAN_2026-11-21.md" "$STAGE/docs/" 2>/dev/null || true

# acceptance snapshot
if [[ -x "$STAGE/bin/coldcrucible-solver" ]]; then
  ( cd "$STAGE" && source bin/env.sh >/dev/null && coldcrucible-solver selfcheck > "$STAGE/docs/selfcheck.json" || true )
fi

echo "Bundle written to: $STAGE"
du -sh "$STAGE"
ls "$STAGE/lib" | wc -l | awk '{print "bundled_lib_files:", $1}'
