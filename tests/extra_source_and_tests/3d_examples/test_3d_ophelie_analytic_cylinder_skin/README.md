# Analytic cylinder skin / self-induction (SPH)

Prescribes uniform axial \(B_0\) via rotational \(A\), runs complex edge-flux Picard,
compares mid-plane \(|J_\theta|\) and \(Q\) to the infinite-cylinder Bessel solution.

```bash
./test_3d_ophelie_analytic_cylinder_skin --dp=0.025 --radius=0.10 --half-height=0.30 \
  --b0=0.01 --sigma=16 --frequency=300000 --self-induction-max-iter=15
```

Doc: `docs/ophelie/34_LONG_CYLINDER_ANALYTIC_SELF_INDUCTION_2026-09-20.md`
