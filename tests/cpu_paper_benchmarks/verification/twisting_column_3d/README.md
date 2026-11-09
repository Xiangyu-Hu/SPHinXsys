# 3D twisting column

Copied from `tests/tests_sycl/3d_examples/test_3d_twisting_column_sycl`,
including its regression assets, at git commit
`19136915145d1d0040e9878c332e23806af2fd7a`. This copy is dedicated to paper
benchmarking; the original test remains unchanged.

## Benchmark options

This case accepts `--benchmark[=true|false]`, `--output[=true|false]`,
`--result-dir PATH`, `--run-id ID`, `--end-time VALUE`,
`--output-interval VALUE`, `--dp VALUE`, and `--resolution LABEL`.
Benchmark options are removed before the remaining arguments are passed to
the original `SPHSystem` command-line parser.

The original defaults are retained: `dp=PH/10` (`0.1`), `end-time=0.5`,
and `output-interval=0.002`. Resolution presets are `coarse=PH/5` (`0.2`),
`standard=PH/10` (`0.1`), and `fine=PH/20` (`0.05`); an explicit `--dp`
overrides the preset. Holder thickness and every domain/shape extent derived
from particle spacing are constructed after parsing.

Benchmark mode disables periodic VTP and Position/Velocity observer writes by
default and skips regression generation/comparison. Use `--output=true` to
retain periodic files during a benchmark. Verification mode retains the
original VTP, Position/Velocity recorders, and dynamic-time-warping regression
flow.

Each run writes `environment.json` and `summary.csv` below
`RESULT_DIR/twisting_column_3d/RUN_ID`. The summary records total particle
count (column plus observer), physical time, executed acoustic/outer steps,
wall/init/I/O/compute timing, compute time per executed step, resolution,
observer quantities, sampling state, and regression status. The environment
file creation is included in measured I/O; the final summary write occurs
after the measured simulation interval.
