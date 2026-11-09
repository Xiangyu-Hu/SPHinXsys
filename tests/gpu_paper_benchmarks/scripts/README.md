# Benchmark scripts

These launchers only run already-built executables. They do not configure or
compile SPHinXsys. Keep two independent build trees:

- `build-host`: `SPHINXSYS_USE_SYCL=OFF`, used by the CPU scripts.
- `build-sycl`: `SPHINXSYS_USE_SYCL=ON`, used by the GPU scripts.

The launcher checks these non-recursive executable locations in order:

```text
<build>/tests/gpu_paper_benchmarks/verification/<case>/bin/<target>
<build>/tests/gpu_paper_benchmarks/verification/<case>/bin/<CONFIG>/<target>
```

Set `BUILD_DIR` to the corresponding build-tree root. If all executables were
copied to one directory, set `SPH_BENCH_BIN_DIR` to that absolute directory
instead. `CONFIG` defaults to `Release`; set it to another multi-config
configuration or to an empty string to check only the direct path. The launcher
does not recursively guess among multiple executables.

## Run suites

These commands are examples only; benchmark execution is left to the user:

```bash
BUILD_DIR=/path/to/build-host REPETITIONS=3 \
  ./tests/gpu_paper_benchmarks/scripts/run_verification_cpu.sh

BUILD_DIR=/path/to/build-sycl REPETITIONS=3 \
ONEAPI_DEVICE_SELECTOR=cuda:gpu SPH_BENCH_DEVICE=rtx-4090 \
  ./tests/gpu_paper_benchmarks/scripts/run_verification_gpu.sh

BUILD_DIR=/path/to/build-host \
  ./tests/gpu_paper_benchmarks/scripts/run_dambreak_scaling_cpu.sh

BUILD_DIR=/path/to/build-sycl ONEAPI_DEVICE_SELECTOR=cuda:gpu \
  ./tests/gpu_paper_benchmarks/scripts/run_dambreak_scaling_gpu.sh
```

Verification runs Poiseuille and diffusion at `coarse`, `standard`, and
`fine`, and dam-break, fish FSI, oscillating beam, and twisting column at
`standard`. Scaling runs dam-break at `s1` through `s6`.

CPU launchers execute an `SPHINXSYS_USE_SYCL=OFF` host binary with
`ONEAPI_DEVICE_SELECTOR` removed from its environment and always record backend
`host_tbb`. GPU launchers always record backend `sycl` and set
`ONEAPI_DEVICE_SELECTOR`, defaulting to `cuda:gpu` for NVIDIA CUDA SYCL builds
(use `level_zero:gpu` for Intel GPU); override the selector to match the
installed oneAPI runtime. Environment
variables cannot switch a CPU launcher onto the SYCL branch or a GPU launcher
onto the host branch.

## Environment

- `BUILD_DIR`: build root (`build-host` or `build-sycl` by default).
- `SPH_BENCH_BIN_DIR`: optional absolute directory containing all benchmark
  executables, overriding `BUILD_DIR`. Do not use oneAPI's `BIN_DIR=bin64`.
- `BIN_DIR`: accepted only if it is an existing absolute path; relative values
  such as oneAPI `bin64` are ignored.
- `CONFIG`: optional multi-config subdirectory, default `Release`; set empty to
  disable the configured-path check.
- `RESULT_ROOT`: raw result root; defaults to `${BUILD_DIR}/benchmark_results`
  (for example `build-sycl/benchmark_results`).
- `REPETITIONS`: positive integer, default `3`.
- `END_TIME_MODE`: `short` (default) or `full`.
- `SHORT_END_TIME`: end time used in short mode, default `0.01`.
- `FULL_END_TIME`: optional common override in full mode. Without it, each
  verification case uses its established full end time.
- `END_TIME`: highest-priority end-time override for every selected run.
- `OUTPUT`: value passed to `--output`, default `off`.
- `OUTPUT_INTERVAL`: optional value passed to `--output-interval`.
- `SPH_BENCH_DEVICE`: optional device label only; defaults to `host_cpu` for CPU
  launchers and the oneAPI selector for GPU launchers.
- `ONEAPI_DEVICE_SELECTOR`: GPU runtime selection; default `cuda:gpu` in
  `run_verification_gpu.sh` / `run_dambreak_scaling_gpu.sh` (use
  `level_zero:gpu` for Intel GPU).

Every run ID includes a nanosecond UTC timestamp, launcher PID, short Git
commit, backend/device, case, resolution, and repeat number. The command and
combined stdout/stderr are first captured in a temporary file because the C++
recorder creates `RESULT_ROOT/case/run-id` itself. Logs and `run.status` are
written without suffixing or replacing existing files. If the executable fails
before owning its run directory, or reports that the directory was concurrently
created, the launcher stores `run.log` and a failed `run.status` in a unique
directory below `RESULT_ROOT/.launcher_failures/<case>/`; it never writes those
artifacts into the conflicting run directory or fabricates a summary row.
The executable recursively stages any executable-directory `input/` and
`reload/` assets into its unique run directory before constructing I/O
recorders. All subsequent reads and writes use those run-local copies.

## Collect results

```bash
python3 tests/gpu_paper_benchmarks/scripts/collect_benchmark_results.py \
  build-sycl/benchmark_results \
  tests/gpu_paper_benchmarks/curated/all_runs.csv \
  --stats-output tests/gpu_paper_benchmarks/curated/repeat_stats.csv
```

The collector recursively finds `summary.csv`, requires every file to begin
with all required columns in their defined order, rejects duplicate/reserved
columns, and permits case-specific extra columns. It emits a deterministic,
alphabetically ordered union of extra columns, reports that union, leaves
missing case-specific values empty, and adds `source_summary`. Duplicate
`(case, run)` identities are errors.

It also scans `run.status`. By default, any failed run directory without a
summary is reported and makes collection fail. Pass `--allow-failures` to keep
collecting while still listing each such directory on stderr. A summary whose
CSV status or `run.status` is not `completed` remains in the combined raw CSV,
is reported, and is excluded from repeated-run statistics.

Optional statistics use the standard library and group by
`case/git/build/precision/backend/device/benchmark_mode/output_enabled/`
`output_interval/resolution/dp/requested_end_time`. These configured fields are
required in every summary. Final `physical_time` is retained in the combined
data but is not a grouping key, so small CPU/GPU end-time overshoots do not
split repeats. Statistics report count, mean, sample standard deviation,
minimum, and median for the standard timing columns. Repeat `--stat-column
NAME` to choose other numeric columns.

## Summary schema notes

Every summary includes shared throughput/memory columns from
`tests/paper_benchmarks_common`: `solid_steps`, `particle_updates`,
`particle_update_definition`, `particle_updates_per_second`, `mpps`,
`neighbor_interactions`, `mnips`, `gpips`, `peak_rss_kb`,
`peak_gpu_memory_kb`, and reserved component-timing columns.
Lightweight observer/feature sampling remains on even when
`--benchmark` disables heavy VTP/restart output.
