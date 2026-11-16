# Johnson–Cook necking with CPU/SYCL computing kernels

This complete example uses the SPHinXsys execution-policy and delegated-variable
interfaces. The same value-only constitutive model is used by classic CPU and
CK/SYCL integration: plastic flow, rate history, plastic-work heating, thermal
softening and cumulative dynamic failure are shared. Failed particles transmit
no stress or pair damping and retain their ballistic mass. Device integration
errors are reported through `JCIntegrationStatus` and stop the calculation.

The circular-arc rod, uniform particle generator, steel constants, symmetric
grips and command-line options are reused from the
[classic example](../../../3d_examples/test_3d_johnson_cook_necking/README.md).
Default spacing is 0.05 mm and end time is 150 microseconds. No external
particle file, mesh or additional solver installation is required.

Configure a normal SPHinXsys SYCL build with an IntelLLVM compiler and the
dependencies required by the repository. Select a target supported by the
installed SYCL toolchain, for example:

```sh
cmake -S . -B build-sycl -DCMAKE_CXX_COMPILER=icpx \
  -DSPHINXSYS_2D=OFF -DSPHINXSYS_3D=ON \
  -DSPHINXSYS_USE_SYCL=ON -DSPHINXSYS_SYCL_TARGETS=spir64 \
  -DSPHINXSYS_BUILD_TESTS=ON
cmake --build build-sycl --target test_johnson_cook_sycl test_3d_johnson_cook_necking_sycl
cd build-sycl/tests/tests_sycl/3d_examples/test_3d_johnson_cook_necking_sycl/bin
./test_3d_johnson_cook_necking_sycl
```

The selected device is printed at startup. Device selection uses the SYCL
runtime's environment settings. SPHinXsys currently selects single precision
for a SYCL build. No FP64 emulation is requested by this example. Use `--host`
to execute the same CK equations and precision on the native host, independently
of the installed SYCL CPU backend:

```sh
./test_3d_johnson_cook_necking_sycl --host
```

Run host and device calculations in separate working directories because both
write the same output names. The program transfers two grip reactions per step;
full field synchronization is performed only for diagnostics and VTP output.
It writes `case.json`, `history.csv`, `summary.json` and `output/*.vtp` using the
classic example's SI units and columns. `case.json` also records `real_bytes`.
The shared `postprocess.py` plots either SPH run using only its own
history and particle-state output.

With SYCL disabled, the target remains a native CK example. The default CTest
uses 0.2 mm uniform spacing and runs to 80 microseconds to exercise damage,
complete failure and boundary reactions. Its `RegressionTestEnsembleAverage`
checks compare total kinetic energy and mean damage at five fixed times against
the included precision-specific SPH data. The regression quantities are reduced
on the host after synchronizing output snapshots; particle-state recording is
not required. Use `--regression-test --spacing=0.0002 --end-time=0.00008
--frames=4` to request this check manually, with the copied `input` directory.
Material checks cover
mixed loading/unloading and rate histories, melting, complete failure and
invalid-input rollback on the selected execution device.

The uniform 0.2 mm functional case was run to 80 microseconds with SYCL GPU
and native CPU execution policies using the same
single-precision executable. Both used 1,648 particles and finished with 32
failed particles. Peak reactions were 1954.290 N (GPU) and 1954.343 N (host),
a difference of -0.00272%. Recorded peak temperatures were 652.147 K and
652.170 K. All three device material tests passed. These are functional checks
at coarse spacing; the default fine-resolution SYCL run has not been rerun.

CK neighbors use SPHinXsys's tabulated smoothing kernel, whereas the classic
example evaluates its kernel directly. Together with single precision and
parallel reductions, this can produce differences that grow after localization
and failure. Bitwise agreement, spatial convergence and experimental validation
are not claimed. Each run should be judged from its own recorded fields.
