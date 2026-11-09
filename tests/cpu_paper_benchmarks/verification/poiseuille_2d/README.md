# 2D Poiseuille flow

Copied from `tests/tests_sycl/2d_examples/test_2d_mixed_poiseuille_flow_sycl`.
The source snapshot is git commit
`19136915145d1d0040e9878c332e23806af2fd7a`; this isolated copy is dedicated to
paper benchmarking and does not alter the original test.

## Benchmark options

The case accepts `--benchmark[=true|false]`, `--output[=true|false]`,
`--result-dir PATH`, `--run-id ID`, `--end-time VALUE`,
`--output-interval VALUE`, `--dp VALUE`, and `--resolution LABEL`. Benchmark
options are removed before the remaining command line is passed to
`SPHSystem`.

Resolution presets are:

- `coarse`: `DH/10` (`0.0001`)
- `standard` or `medium`: `DH/20` (`0.00005`, the original test spacing)
- `fine`: `DH/40` (`0.000025`)

An explicit `--dp` overrides `--resolution`. The original defaults are
`end-time=2.0` and `output-interval=0.25`.

Benchmark mode disables periodic VTP and observer output unless
`--output=true` is supplied. One final observer sample is still written and
timed as I/O because the analytical velocity validation requires freshly
interpolated values. Verification mode retains the original periodic output.

Each run writes `environment.json` and `summary.csv` below
`RESULT_DIR/poiseuille_2d/RUN_ID`. The summary includes particle counts,
physical time, loop counts, wall/init/I/O/compute timing, and validation
pass/fail counts. The reported analytical metric is the maximum of the
existing pointwise normalized error `abs((u_sim-u_exact)/U_f)`, tested against
the original `0.05` tolerance; it is not an L2 or Linf field.
