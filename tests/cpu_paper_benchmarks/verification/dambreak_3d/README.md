# 3D dambreak

Copied from `tests/tests_sycl/3d_examples/test_3d_dambreak_sycl`, including its
regression assets, at git commit
`19136915145d1d0040e9878c332e23806af2fd7a`. This isolated copy is the first
fully instrumented paper benchmark template.

Particle counts are read from the public `TotalRealParticles()` API. Wall time
includes initialization, compute, and explicitly measured output time; compute
time excludes initialization and measured I/O. SPHinXsys exposes no stable
public queue-wait API to this case, so overall timing relies on the existing
`exec()` dependency and host-visible reduction semantics and must be validated
on the target GPU.
