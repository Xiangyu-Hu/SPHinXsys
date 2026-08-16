# Benchmark scripts (CPU paper)

These launchers only run already-built executables. They do not configure or
compile SPHinXsys. Prefer a host build tree:

- `build-host`: `SPHINXSYS_USE_SYCL=OFF`, `SPHINXSYS_BUILD_CPU_PAPER_BENCHMARKS=ON`

GPU campaigns belong under `tests/gpu_paper_benchmarks/scripts/` with a separate
`build-sycl` tree. Both trees share the recorder schema in
`tests/paper_benchmarks_common`.

The launcher checks these non-recursive executable locations in order:

```text
<build>/tests/cpu_paper_benchmarks/verification/<case>/bin/<target>
<build>/tests/cpu_paper_benchmarks/verification/<case>/bin/<CONFIG>/<target>
```

Set `BUILD_DIR` to the host build-tree root. If all executables were copied to
one directory, set `SPH_BENCH_BIN_DIR` to that absolute directory instead.
`CONFIG` defaults to `Release`; set it to another multi-config configuration or
to an empty string to check only the direct path.

## Run suites

Examples only; compilation and execution remain user-driven:

```bash
BUILD_DIR=/path/to/build-host REPETITIONS=3 \
  ./tests/cpu_paper_benchmarks/scripts/run_verification_cpu.sh

BUILD_DIR=/path/to/build-host \
  ./tests/cpu_paper_benchmarks/scripts/run_dambreak_scaling_cpu.sh
```

Verification runs Poiseuille and diffusion at `coarse`, `standard`, and
`fine`, and dam-break, fish FSI, oscillating beam, and twisting column at
`standard`. Scaling runs dam-break at `s1` through `s6`.

CPU launchers execute an `SPHINXSYS_USE_SYCL=OFF` host binary with
`ONEAPI_DEVICE_SELECTOR` removed from its environment and always record backend
`host_tbb`.

## Summary schema notes

Every summary includes shared throughput/memory columns from
`tests/paper_benchmarks_common`: `solid_steps`, `particle_updates`,
`particle_update_definition`, `particle_updates_per_second`, `mpps`,
`neighbor_interactions`, `mnips`, `gpips`, `peak_rss_kb`,
`peak_gpu_memory_kb`, and reserved component-timing columns.
Lightweight observer/feature sampling remains on even when
`--benchmark` disables heavy VTP/restart output.

Use the shared collector:

```bash
python3 tests/cpu_paper_benchmarks/scripts/collect_benchmark_results.py \
  /path/to/build-host/benchmark_results \
  tests/cpu_paper_benchmarks/curated/all_runs.csv \
  --stats-output tests/cpu_paper_benchmarks/curated/repeat_stats.csv
```
