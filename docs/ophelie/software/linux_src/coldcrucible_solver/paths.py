from __future__ import annotations

import os
from pathlib import Path

_MODULE_DIR = Path(__file__).resolve().parent


def _find_repo_root(start: Path) -> Path:
    """Walk up until CMakeLists.txt + tests/ look like SPHinXsys repo root."""
    for p in [start, *start.parents]:
        if (p / "CMakeLists.txt").is_file() and (p / "tests").is_dir() and (p / "src").is_dir():
            return p
    # Fallback: linux_src/coldcrucible_solver -> parents until docs/ophelie/software/linux_src then up to repo
    # coldcrucible_solver(0)->linux_src(1)->software(2)->ophelie(3)->docs(4)->repo(5)
    return start.parents[5]


REPO_ROOT = Path(os.environ["OPHELIE_REPO_ROOT"]).resolve() if os.environ.get("OPHELIE_REPO_ROOT") else _find_repo_root(_MODULE_DIR)

# Packaged under docs/ophelie/software/linux_src (dev handoff layout)
LINUX_SRC_ROOT = _MODULE_DIR.parent  # .../linux_src
IN_BUNDLE = (_MODULE_DIR.parent.name == "share")

DEFAULT_BUILD = Path(os.environ.get("OPHELIE_BUILD", REPO_ROOT / "build"))
EXTRA_3D = DEFAULT_BUILD / "tests/extra_source_and_tests/3d_examples"

TARGETS = {
    "nc_relax": EXTRA_3D / "test_3d_ophelie_client_nc_glass_relax/bin/test_3d_ophelie_client_nc_glass_relax",
    "nc_run": EXTRA_3D
    / "test_3d_ophelie_french_natural_convection_frozen_q/bin/test_3d_ophelie_french_natural_convection_frozen_q",
    "stir_relax": EXTRA_3D / "test_3d_ophelie_client_stirring_glass_relax/bin/test_3d_ophelie_client_stirring_glass_relax",
    "stir_f0": EXTRA_3D / "test_3d_ophelie_client_stirring_f0_elastic_fsi/bin/test_3d_ophelie_client_stirring_f0_elastic_fsi",
    "stir_f1": EXTRA_3D / "test_3d_ophelie_client_stirring_f1_heated_elastic/bin/test_3d_ophelie_client_stirring_f1_heated_elastic",
}


def examples_root() -> Path:
    env = os.environ.get("OPHELIE_EXAMPLES_ROOT")
    if env:
        return Path(env).resolve()
    if IN_BUNDLE:
        return REPO_ROOT / "examples"
    cand = LINUX_SRC_ROOT / "examples"
    if cand.is_dir():
        return cand
    return REPO_ROOT / "software/examples"


def resolve_binary(key: str) -> Path:
    env_key = f"OPHELIE_BIN_{key.upper()}"
    override = os.environ.get(env_key)
    if override:
        return Path(override)
    if IN_BUNDLE:
        name_map = {
            "nc_relax": "test_3d_ophelie_client_nc_glass_relax",
            "nc_run": "test_3d_ophelie_french_natural_convection_frozen_q",
            "stir_relax": "test_3d_ophelie_client_stirring_glass_relax",
            "stir_f0": "test_3d_ophelie_client_stirring_f0_elastic_fsi",
            "stir_f1": "test_3d_ophelie_client_stirring_f1_heated_elastic",
        }
        bundled = REPO_ROOT / "solver" / name_map[key] / name_map[key]
        if bundled.is_file():
            return bundled
    return TARGETS[key]
