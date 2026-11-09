# 2D fish FSI

Copied from
`tests/tests_sycl/2d_examples/test_2d_flow_stream_around_fish_sycl`. The copy
was taken at git commit `19136915145d1d0040e9878c332e23806af2fd7a` and is
isolated for paper benchmarking so the original test remains unchanged.

## Benchmark controls

This copy accepts the common benchmark options before forwarding all remaining
arguments to `SPHSystem::handleCommandlineOptions`:

- `--benchmark[=true|false]`
- `--resolution standard`
- `--dp <metres>`
- `--end-time <seconds>`
- `--output[=true|false]`
- `--output-interval <seconds>`
- `--result-dir <path>`
- `--run-id <safe-name>`

The `standard` setting preserves the original case:
`dp=0.0025`, `end-time=1.7`, and `output-interval=0.01`. Other paper
resolutions must use explicit `--dp`; no unvalidated named levels are implied.
An explicit `--dp` takes precedence over `--resolution`.

All spacing-dependent geometry (`DL_sponge`, boundary width, emitter/disposer
boxes, and fish polygon cutoff) is recomputed before `SPHSystem`, bodies, or
the solver are created. A spacing above 0.0075 m is rejected, and the code also
rejects a spacing for which the polygon has no resolvable tail-center point.
This is a geometry-safety bound, not a claim that every accepted spacing has
been numerically validated.

## Output and measurements

Without `--benchmark`, output behavior remains the original verification
behavior by default: VTP states and mechanical energy are written, and restart
files are written every 500 advection steps. Benchmark mode disables VTP by
default, disables periodic mechanical-energy records, and disables restart
writes. `--output=true` can explicitly restore VTP in benchmark mode.

One lightweight kinematic feature is always sampled at the output interval:
the interpolated fish `Position` at a tail-center reference point. The point is
the rearmost centerline location retained by the same 100-segment polygon and
`dp/2` cutoff used to construct the fish. It therefore moves inward
automatically when a coarser spacing truncates more of the thin tail. The
implementation uses the repository's existing SYCL
`addObserveRecorder<Vecd>(contact, "Position")` pattern and refreshes the
observer contact relation before sampling. The resulting observer series is
written by SPHinXsys in its normal output directory; the benchmark summary
records the exact reference coordinates and sample count. To avoid presenting
an under-resolved body as supported, explicit spacing is limited to
`dp <= 0.0075` m (at least four particles across the maximum fish thickness).

The common recorder writes `environment.json` and `summary.csv` under
`<result-dir>/fish_fsi_2d/<run-id>/`. The summary includes initial/final
particle counts, final physical time, output/advection/acoustic counts, solid
substeps, wall/initialization/I/O/compute times, and compute time per advection
step. `outer_steps` denotes periodic feature/output triggers, while
`advection_steps` counts steps executed by this run (not the restart index).
Timed simulation file and console writes are accumulated as I/O; the
final summary write itself is necessarily outside that completed timing
snapshot.
