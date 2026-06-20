# test_3d_ophelie_phi_device_vector_ops

P1 scaffold: SYCL volume-weighted `dot` / `norm` / `axpy` vs host (`electromagnetic_ophelie_phi_device_vector_ops.h`).

```bash
cd ~/SPHinXsysSYCL/build
cmake ..
cmake --build . --target test_3d_ophelie_phi_device_vector_ops -j$(nproc)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_device_vector_ops/bin/test_3d_ophelie_phi_device_vector_ops
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_device_vector_ops/bin/test_3d_ophelie_phi_device_vector_ops --reload=1
```

`--gmres-parity` compares host vs device GMRES `eq_res_vol` (lattice, ~1 s).

Reload-scale GMRES parity (n≈16066, ~2× partial GMRES, several minutes):

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_device_vector_ops/bin/test_3d_ophelie_phi_device_vector_ops --reload=1 --gmres-parity
```

Production GMRES on GPU (SphinxSys `particle_reduce` + `LoopRangeCK` + `StateDynamics` CK):

```bash
--phi-gmres-device-ops=1
--phi-gmres-device-ops=1 --phi-gmres-device-krylov=1
```

Implementation: `electromagnetic_ophelie_phi_krylov_ck.h` (dot/norm/axpy/subtract CK); Krylov basis in USM shared, staged to `PhiKrylov*` particle fields for reductions.
