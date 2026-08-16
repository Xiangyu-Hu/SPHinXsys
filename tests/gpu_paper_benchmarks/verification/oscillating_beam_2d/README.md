# 2D oscillating beam

Copied from `tests/tests_sycl/2d_examples/test_2d_oscillating_beam_sycl`,
including its regression assets, at git commit
`19136915145d1d0040e9878c332e23806af2fd7a`. This copy is dedicated to paper
benchmarking; the original test remains unchanged.

## Benchmark options

This case accepts `--benchmark[=true|false]`, `--output[=true|false]`,
`--result-dir PATH`, `--run-id ID`, `--end-time VALUE`,
`--output-interval VALUE`, `--dp VALUE`, and `--resolution LABEL`.
Benchmark options are removed before the remaining arguments are passed to
the original `SPHSystem` command-line parser.

The original defaults are retained: `dp=PH/10` (`0.002`),
`end-time=1.0`, and `output-interval=0.01`. Resolution presets are
`coarse=PH/5` (`0.004`), `standard=PH/10` (`0.002`), and
`fine=PH/20` (`0.001`); an explicit `--dp` overrides the preset.

Benchmark mode disables periodic VTP and beam-tip Position writes by default
and skips regression generation/comparison. Use `--output=true` to retain
periodic files during a benchmark. Verification mode retains the original VTP,
Position recorder, and ensemble-average regression flow.

Each run writes `environment.json` and `summary.csv` below
`RESULT_DIR/oscillating_beam_2d/RUN_ID`. The summary records total particle
count (beam plus observer), physical time, executed acoustic/outer steps,
wall/init/I/O/compute timing, compute time per executed step, resolution,
observer quantity, sampling state, and regression status. The environment
file creation is included in measured I/O; the final summary write occurs
after the measured simulation interval.
