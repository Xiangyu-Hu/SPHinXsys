# French Reduced Case

Parameterised cylindrical glass + prescribed multiloop circular filament coil.  
**French-paper-inspired reduced geometry** — not exact CAD reconstruction.

Full assumptions: [FRENCH_REDUCED_CASE_ASSUMPTIONS.md](FRENCH_REDUCED_CASE_ASSUMPTIONS.md)

Test README: [test_3d_ophelie_french_reduced/README.md](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/README.md)

## Run

```bash
cd /home/yyc/SPHinXsysSYCL/build

# Level0 smoke
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --dp=0.02 --no-phi --state_recording=0

# With phi + VTP
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --dp=0.02 --state_recording=1

# Visual bodies (coil annulus + crucible wall shell)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --dp=0.02 --coil-visual --crucible-visual --state_recording=1
```

## Defaults

| Category | Parameter | Value |
|----------|-----------|-------|
| Literature | D | 650 mm |
| Literature | f | 300 kHz |
| Literature | σ | 16 S/m |
| Literature | P_target | 50 kW |
| Assumed | glass H | 0.50 m |
| Assumed | coil R | 0.40 m |
| Assumed | loops | 8 |
| Assumed | dp | 0.02 m |

## Regression tests

```bash
for t in test_3d_ophelie_biot_savart_circular_loop_axis test_3d_ophelie_biot_savart_multiloop_axis \
         test_3d_ophelie_phi_laplace_constant test_3d_ophelie_phi_rhs_constant_A \
         test_3d_ophelie_phi_reduces_divJ test_3d_ophelie_joule_to_heat_one_way; do
  ./tests/extra_source_and_tests/3d_examples/$t/bin/$t
done
```
