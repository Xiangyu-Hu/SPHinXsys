# 2D diffusion with Neumann boundary condition

Copied from `tests/tests_sycl/2d_examples/test_2d_diffusion_NeumannBC_sycl`,
including separate host CK and SYCL regression assets. The original test
snapshot is git commit `19136915145d1d0040e9878c332e23806af2fd7a`; this copy is
dedicated to paper benchmarking and remains independent of the source case.

## Benchmark options

The case accepts `--benchmark[=true|false]`, `--output[=true|false]`,
`--result-dir PATH`, `--run-id ID`, `--end-time VALUE`,
`--output-interval VALUE`, `--dp VALUE`, and `--resolution LABEL`. Benchmark
options are removed before the remaining command line is passed to
`SPHSystem`.

Resolution presets are:

- `coarse`: `H/50` (`0.02`)
- `standard` or `medium`: `H/100` (`0.01`, the original test spacing)
- `fine`: `H/200` (`0.005`)

An explicit `--dp` overrides `--resolution`. The original defaults are
`end-time=1.0` and `output-interval=0.1`.

Benchmark mode disables periodic VTP and observer writes unless
`--output=true` is supplied, and always skips regression database generation
and comparison I/O. Verification mode retains the original output and
regression/reference flow.

Each run writes `environment.json` and `summary.csv` below
`RESULT_DIR/diffusion_neumann_2d/RUN_ID`. The summary includes particle
counts, physical time, diffusion and outer loop counts, wall/init/I/O/compute
timing, whether observer data was sampled, and whether the reference was
generated, tested, skipped for restart, or not checked in benchmark mode.
When verification output is explicitly disabled, regression checking is also
reported as not checked because there are no observation snapshots to compare.
The current regression API does not expose L2/Linf values here; extracting
those metrics from regression artifacts remains a post-processing task.
