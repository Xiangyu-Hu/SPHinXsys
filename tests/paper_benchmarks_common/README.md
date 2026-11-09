# Shared paper-benchmark recorder

Used by both `tests/gpu_paper_benchmarks` and `tests/cpu_paper_benchmarks`.

Keep the summary schema identical across GPU and CPU campaigns. Derived
throughput fields use `particle_updates = particle_count * <case step metric>`
with an explicit `particle_update_definition`. Neighbor-based MNIPS/GPIPS and
peak GPU memory remain blank/`unavailable` until a portable measurement exists.
