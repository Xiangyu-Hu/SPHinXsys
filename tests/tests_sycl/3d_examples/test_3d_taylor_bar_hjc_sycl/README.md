# HJC Taylor-bar impact with SYCL

A concrete cylinder strikes a fixed wall at 30 m/s. The geometry, material
parameters and 60 microsecond duration follow the
[CPU HJC example](../../../3d_examples/test_3d_taylor_bar_hjc/README.md).
This example uses the same HJC constitutive update on the host and device,
with the library's CK particle interactions and wall contact.

## Build and run

Install the dependencies described in the repository's build instructions.
For a SYCL build, use the Intel compiler and select a target supported by
your device. For example, `spir64` allows runtime device selection:

```sh
cmake -S . -B build-sycl -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=icpx -DSPHINXSYS_USE_SYCL=ON \
  -DSPHINXSYS_SYCL_TARGETS=spir64 -DSPHINXSYS_2D=OFF \
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=OFF \
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF
cmake --build build-sycl --target test_3d_taylor_bar_hjc_sycl test_hjc_sycl
```

SPHinXsys uses single precision for SYCL. A CPU build of the same example
can be configured with `SPHINXSYS_USE_SYCL=OFF` and
`SPHINXSYS_USE_FLOAT=ON` for a comparison at the same precision.

Run each calculation in its own directory:

```sh
mkdir impact-gpu
cd impact-gpu
ONEAPI_DEVICE_SELECTOR=level_zero:gpu \
  ../build-sycl/tests/tests_sycl/3d_examples/test_3d_taylor_bar_hjc_sycl/test_3d_taylor_bar_hjc_sycl
```

The default particle spacing is 0.5 mm. `--spacing`, `--speed`, `--cfl`
and `--end-time` control the calculation; the maximum CFL is 0.2.
`history.csv` records contact force, axial velocity, damage, kinetic
energy and the minimum deformation determinant at 31 fixed times.
VTP files contain the particle damage, stress and compaction fields.

```sh
ctest --test-dir build-sycl -R '^(test_hjc_sycl|test_3d_taylor_bar_hjc_sycl)$' --output-on-failure
python tests/tests_sycl/3d_examples/test_3d_taylor_bar_hjc_sycl/plot.py impact-cpu impact-gpu --output figures
```

The impact regression uses 1 mm spacing with the included energy and mean
damage reference data. CMake selects the single- or double-precision data
to match the build. The single-precision references include native CPU and
GPU runs with Level Zero and OpenCL, generated using the library's ensemble
regression procedure. Ordinary runs do not require the reference files.
The material tests compare prescribed loading paths with the CPU material
interface and check rigid rotation, irreversible history and invalid inputs.

## Response and particle fields

The figures use 0.5 mm spacing (12,640 concrete particles), single precision,
and the Level Zero GPU backend. The CPU comparison runs the same CK example.
The OpenCL GPU backend was checked separately. The differences below describe
these runs; they are not accuracy bounds for other devices or loading paths.

| Largest sampled CPU/GPU difference | Level Zero | OpenCL |
| --- | ---: | ---: |
| Contact force | 1.24 N | 1.05 N |
| Kinetic energy | 4.89e-5 J | 5.68e-5 J |
| Mean damage | 3.38e-4 | 1.66e-4 |

![CPU and GPU response curves](response.png)

![GPU damage and equivalent stress at 20, 40 and 60 microseconds](impact.png)

The cutaway shows particle values without smoothing. Local damage is more
sensitive than the global response: at 60 microseconds, its particlewise
CPU/GPU RMS difference is 0.010–0.015, with maximum differences of 0.20–0.31
across the two backends. Similar mean damage does not establish pointwise
agreement or convergence of a fracture pattern.

## Numerical details

`HJCIntegration1stHalfCK` retains the incremental logarithmic strain,
objective stress rotation, total Lagrangian force and pair damping used by
the CPU HJC implementation. `HJCAcousticTimeStepCK` accounts for the evolving
EOS modulus. The second half step uses `StructureIntegration2ndHalf`.

The wall uses `RepulsionFactor` and `RepulsionForceCK<Wall>` with their
standard settings. This contact force uses the wall normal and contact
damping; the original CPU example uses a different contact discretization.
CPU/device comparisons should therefore use this same CK example.

Device updates write `HJCIntegrationStatus`. The example checks it after
every first half step and stops on an invalid constitutive state. An invalid
update does not commit material history; it does not roll back the complete
particle time step.

The HJC assumptions and limitations remain those of the CPU example.
Damage is local and particles are retained after full damage. The example
demonstrates constitutive integration rather than a calibrated fracture test.
The coarse regression and fine illustration are not a mesh-convergence study.

Single precision can accumulate stress drift during many small rigid rotations.
Near complete damage, the small remaining strength makes this drift capable of
producing additional plastic strain and damage. This also occurs in the host
single-precision update; use a double-precision comparison when studying such
loading paths. In a 4096-step pure-rotation material-point check with initial
damage 0.99, GPU damage increased by about 0.0018 even though the prescribed
motion contained no strain.

## References

1. Holmquist, T. J., Johnson, G. R., and Cook, W. H. (1993).
   *A computational constitutive model for concrete subjected to large strains,
   high strain rates, and high pressures*. Proceedings of the 14th International
   Symposium on Ballistics, Quebec, pp. 591–600.
2. Meyer, C. S. (2011). [*Development of Geomaterial Parameters for Numerical
   Simulations Using the Holmquist-Johnson-Cook Constitutive Model for Concrete*](https://www.govinfo.gov/content/pkg/GOVPUB-D101-PURL-gpo10967/pdf/GOVPUB-D101-PURL-gpo10967.pdf).
   ARL-TR-5556.
