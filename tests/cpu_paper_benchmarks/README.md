# CPU paper benchmarks

Shared recorder/config live in `tests/paper_benchmarks_common` (identical schema with `tests/gpu_paper_benchmarks`).

This paper-only tree was copied from `tests/tests_sycl` at git commit
`19136915145d1d0040e9878c332e23806af2fd7a` and adapted with benchmark-only
configuration, timing, and result recording.
Enable these long-running executables with
`SPHINXSYS_BUILD_CPU_PAPER_BENCHMARKS=ON`; they are not registered with CTest.
Build these targets in a host tree (`SPHINXSYS_USE_SYCL=OFF`, typically `build-host`).

Use separate build roots for reproducible CPU/GPU comparisons, for example
`build-host` and `build-sycl`. The launchers in [`scripts/`](scripts/README.md)
consume those existing builds; they never configure, compile, or submit
benchmarks themselves. Actual benchmark runs must be started by the user after
checking the selected build, device, end time, output policy, and result path.

The shared parser removes only these options before SPHinXsys sees `argv`:
`--benchmark`, `--dp`, `--resolution`, `--end-time`, `--output`,
`--output-interval`, `--result-dir`, and `--run-id`. Options accept either
`--name=value` or `--name value`; boolean options also work as flags. Example:

```text
cpu_paper_dambreak_3d --benchmark --resolution=s3 --end-time=20 \
  --run-id=host-s3-01
```

If `--result-dir` is omitted, each executable defaults to
`<build-dir>/benchmark_results` (for example `build-host/benchmark_results`),
injected at CMake configure time for that build tree.

`standard` and `s1` currently map to `dp=0.050`; `s2` through `s6` map to
`0.040`, `0.032`, `0.025`, `0.020`, and `0.016`. These are explicit initial
values for later hardware calibration. `--dp` overrides the label mapping.

Without `--benchmark`, the dambreak template retains verification output and
regression checking. Benchmark mode suppresses heavy periodic VTP/restart output unless
`--output=on` is supplied. Lightweight feature/observer sampling still
runs so kinematics remain recorded. Final regression testing/database
generation still requires verification mode with output enabled.

Each run gets a non-overwriting directory containing `summary.csv`,
`environment.json`, and all SPHinXsys-generated `output/`, `restart/`, and log
files. Before constructing `SPHSystem`, the executable changes into this unique
run directory, so concurrent or repeated runs never write to the shared binary
directory. After `SPHSystem` initializes its run-local I/O environment and
before any I/O recorder is constructed, existing binary-directory `input/` and
`reload/` trees are recursively staged into the corresponding run directory.
Readers and writers therefore use only run-local paths, so regression
generation, relaxation, and reload writes cannot modify shared assets.
Git, build, and compiler metadata are injected at CMake
configure time; device/host details may be supplied with `SPH_BENCH_DEVICE`,
`SPH_BENCH_HOST`, and `SPH_BENCH_OS`. Raw per-run data belongs under the active
build tree in `<build-dir>/benchmark_results/` (for example
`build-host/benchmark_results/` or `build-host/benchmark_results/`); reviewed
aggregate material belongs in [`curated/`](curated/README.md).

The CPU verification/scaling scripts default to `build-host`, deliberately
remove `ONEAPI_DEVICE_SELECTOR`, and fix backend metadata to `host_tbb`; the CPU
baseline is the native host build. The GPU scripts default to `build-host`, fix
backend metadata to `sycl`, and always set `ONEAPI_DEVICE_SELECTOR` (for NVIDIA
CUDA SYCL builds use `cuda:gpu`; for Intel GPU use `level_zero:gpu`).
`BUILD_DIR`, `BIN_DIR`, multi-config `CONFIG`, `RESULT_ROOT`, `REPETITIONS`,
short/full end-time controls, and output controls are documented in the scripts
README.

Each launcher preserves the exact command and combined stdout/stderr in
`run.log`, including failed attempts, and isolates pre-directory failures from
real run directories. The collector validates ordered required columns, unions
case-specific extra columns, reports failed `run.status` entries, keeps
non-completed summaries only in raw combined output, and can write completed-run
timing statistics using only the Python standard library. Every summary records
the requested end time, benchmark/output mode, output interval, resolution, and
verification time as required columns. Repeat statistics group by requested end
time and these configured run settings, while final physical time remains a
measured result rather than a grouping key.


## CPU build note

CMake wiring and case copies are ready. Prefer a separate host tree (`SPHINXSYS_USE_SYCL=OFF`, e.g. `build-host`) with `SPHINXSYS_BUILD_CPU_PAPER_BENCHMARKS=ON`. GPU compilation can proceed independently; CPU binaries do not need to be built until baseline campaigns start.
