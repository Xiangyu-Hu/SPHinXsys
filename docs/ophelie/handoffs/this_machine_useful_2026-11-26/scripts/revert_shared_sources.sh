#!/usr/bin/env bash
set -euo pipefail
PACK=$(cd "$(dirname "$0")/.." && pwd)
ROOT=$(cd "$PACK/../../../.." && pwd)
cd "$ROOT"
git checkout HEAD -- \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_client_geometry.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_client_case_params.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_client_stirring_geometry.h \
  tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_french_stirring_geometry.h \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/test_3d_ophelie_french_natural_convection_frozen_q.cpp \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/test_3d_ophelie_client_nc_glass_relax.cpp \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/test_3d_ophelie_french_stirring_em.cpp \
  tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/test_3d_ophelie_client_stirring_f1_heated_elastic.cpp \
  docs/ophelie/reference/doc248_client_cases/01_generate_particles.sh \
  docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh \
  docs/ophelie/reference/doc248_client_cases/03_run_C2_stir.sh \
  docs/ophelie/reference/doc248_client_cases/COMPILE.sh \
  docs/ophelie/reference/doc248_client_cases/README.md \
  docs/ophelie/reference/doc248_client_cases/configs/C0_effective.json \
  docs/ophelie/reference/doc248_client_cases/configs/C2_effective.json \
  docs/ophelie/reference/doc248_client_cases/configs/geometry_and_probe_mapping.json
echo "[pack] shared sources restored to HEAD; archive remains in $PACK"
