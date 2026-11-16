# Tensile-necking regression data

The example CTests compare total kinetic energy and mean `JCDamage` at
0, 20, 40, 60 and 80 microseconds using `RegressionTestEnsembleAverage`.
The fixed case uses 0.2 mm uniform spacing, 1,648 particles and four output
intervals. `--regression-test` rejects other geometry/time settings.

The `double` samples come from completed classic CPU and native CK runs.
The `float` samples come from native-host and SYCL GPU CK runs.
The files store the two-run mean and the larger of the observed squared
deviation and the framework's 1% scale floor. They are deterministic regression
references, not a statistical-convergence or spatial-accuracy assessment.

The data are copied into each CTest working directory's `input` folder according
to the build precision. Verification never regenerates or relaxes references.
Device fields are synchronized at output times before the host reduces the two
regression quantities. Particle-state recording is unnecessary for verification.

To regenerate after an intentional numerical change, first run two complete
cases in separate directories, with `--spacing=0.0002 --end-time=0.00008
--frames=4 --state_recording=true`. Then run:

```sh
python generate_reference.py --run /path/to/run1 --run /path/to/run2 --output double
```

For single precision, use single-precision runs and select `--output float`.
Review the changes to both means and envelopes before replacing the baselines.
