# Johnson–Cook tensile necking

This three-dimensional example stretches a circular rod with a smooth waist.
It exercises strain hardening, rate dependence, adiabatic heating, thermal
softening, and complete Johnson–Cook dynamic failure. The example generates
its geometry and uniform particles in C++; no external mesh, particle file,
Python installation, or additional solver is needed to run it.

The material supports both the classic CPU solid dynamics interface and shared
CPU/SYCL computing kernels. It also supports plane-strain stress histories for
two-dimensional bodies. The companion
[`test_3d_johnson_cook_necking_sycl`](../../tests_sycl/3d_examples/test_3d_johnson_cook_necking_sycl)
example reuses this geometry, material parameters and loading.

## Build and run

From a normally configured SPHinXsys checkout:

```sh
cmake -S . -B build -DSPHINXSYS_3D=ON -DSPHINXSYS_BUILD_TESTS=ON
cmake --build build --target test_3d_johnson_cook_necking test_johnson_cook
cd build/tests/3d_examples/test_3d_johnson_cook_necking/bin
./test_3d_johnson_cook_necking
```

The default run uses 0.05 mm spacing, 150 microseconds, and 60 output intervals.
Options are `--spacing`, `--end-time`, `--frames`, `--threads`, and
`--state_recording=false`. Spacing must divide the initial diameter and length.
All units in the implementation and CSV files are SI.

The CTest case uses a small uniform lattice and compares total kinetic energy
and mean damage against the included SPH regression data at five fixed times.
Use `--regression-test --spacing=0.0002 --end-time=0.00008 --frames=4` to run
the same check manually from a directory containing the copied `input` files.
The material tests cover analytical elastic and
linear-hardening responses, heat generation, melting, the failure stress sign,
complete failure, and rotation covariance.
The [regression data](regression_test_tool/README.md) document the fixed case,
precision-specific samples and reference-generation procedure.

```sh
ctest --test-dir build -R 'JohnsonCook\.|test_3d_johnson_cook_necking' --output-on-failure
```

## Geometry and loading

The rod is 6.4 mm long, with nominal radius 0.8 mm and waist radius 0.75 mm.
Two equal circular arcs of radius 2.825 mm connect the waist to the straight
shaft over each 0.75 mm half-transition. Position and tangent are continuous.
There is no initial crack. The 0.05 mm Cartesian lattice contains 102,320
equal-volume particles; the initial surface remains discretized by particles.

Regions at `|z| >= 2.8 mm` are gripped. Their transverse displacement is fixed;
their axial velocities ramp symmetrically to ±10 m/s over 20 microseconds with
a cosine ramp. The CFL is 0.3. The acoustic time step, total Lagrangian pair
forces and Verlet updates use the existing solid dynamics implementation.
Completely failed particles transmit neither stress nor numerical damping;
their mass remains in ballistic motion. The reference neighborhood remains
fixed, so surviving particles can transmit forces across a damaged region.

The boundary reaction includes the prescribed acceleration of the grip mass.
The example writes `history.csv`, `case.json`, `summary.json`, and the standard
`output/*.vtp` particle states. Work, kinetic energy, and plastic work are
recorded; a complete energy closure is not claimed.

## Material

The verification steel has density 7830 kg/m³, Young's modulus 200 GPa and
Poisson ratio 0.29. The JC constants are A=792 MPa, B=510 MPa, n=0.26,
C=0.014, m=1.03, and reference rate 1 s⁻¹. Reference and melting temperatures
are 298 K and 1793 K; heat capacity is 477 J/(kg·K). Plastic work is converted
to heat with fraction 1.0. D1–D5 are 0.05, 3.44, −2.12, 0.002, and 0.61.

The explicit flow update uses the preceding plastic strain rate. Failure uses
the current deviatoric total strain rate and preceding temperature. Hydrostatic
pressure is `K * (1/J - 1)`. Failure is accumulated as `sum(delta_ep / ep_failure)`;
all stress components become zero at damage 1. No fracture-energy regularization
or post-initiation damage evolution is specified.

## Postprocessing

Optional postprocessing requires NumPy and Matplotlib:

```sh
python /path/to/this/example/postprocess.py /path/to/run/directory
```

This command reads only the SPH run's `case.json`, `history.csv` and any available
ASCII VTP particle states. It writes raw response plots and material sections
at exactly matching 35 and 80 microseconds. Runs without particle-state output
still produce the response plot from their history. The script uses only the
files written by the SPH example.
