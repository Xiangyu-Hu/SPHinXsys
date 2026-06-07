# Stage 10.11 — Coulomb Projection design notes (design only, no implementation)

> **Date**: 2026-05-21  
> **Status**: design document; **not entering code mainline yet** (ChatGPT Stage 10.11-D)  
> **Trigger conditions**: see §7

---

## 1. Goal

After frequency-domain A–φ solve, project A onto approximately divergence-free subspace:

```text
solve Laplace(chi) = div A
A_new = A_old - grad chi
```

If E must remain unchanged, may synchronously update φ (phasor convention to be aligned):

```text
phi_new = phi_old + j omega chi   (sign and convention to be verified)
```

---

## 2. Open design questions

| # | Question | Current inclination |
|---|------|----------|
| 1 | Which Laplace operator for χ equation? | **pairwise uncorrected** Laplace same family as A equation |
| 2 | χ boundary conditions? | box: Dirichlet χ=0 or Neumann ∂χ/∂n=0; Contact undefined |
| 3 | Which gradient for grad χ? | same family as penalty/apply **pairwise uncorrected** |
| 4 | Update φ? | If E = -jωA - grad φ, φ correction needed after A projection to keep E |
| 5 | Does E remain unchanged? | acceptance criterion: **E_combined** change before/after projection < tolerance |
| 6 | How for Contact multibody? | needs Contact Poisson + interface χ conditions; **not designed** |
| 7 | Embed in iteration or post-solve? | first **post-solve correction** prototype; do not embed in GMRES |

---

## 3. Why not implement yet

1. Stage 10.10 shows **core region pairwise divA already correct**; large global divA mainly from boundary support artifact.  
2. projection may **over-correct boundary artifact**, introducing non-physical E/J changes.  
3. pairwise A-penalty not yet complete with boundary + exact observable full matrix validation.  
4. Contact Poisson and interface BC complexity high.

---

## 4. Relationship with A-penalty

| Method | Role | Status |
|------|------|------|
| **η_A penalty** | soft constraint divA during solve | Inner research prototype |
| **projection** | hard constraint divA after solve | design only |
| **B-corrected diagnostic** | posterior report divA | dual track adopted |

The three should not mix different discretization families.

---

## 5. Phasor convention notes

Current convention (Az2D, `φ=0`, A as real phasor):

```text
E_real = omega * A_imag - grad(phi_real) = 0
E_imag = -omega * A_real - grad(phi_imag) = -omega * A_real
```

After projection, if only A changes not φ, `E_imag` will change. If E must remain unchanged, φ correction is required.

---

## 6. Recommended prototype acceptance (if implemented in future)

```text
1. Az2D exact field, η=0, post-solve projection
2. core divA (gradDen) decreases
3. E_combined / Joule vs exact does not worsen (or improves)
4. B=curlA vs exact does not worsen
5. GMRES baseline (no projection) still pass
```

---

## 7. Conditions to enter implementation (ChatGPT 10.10-E / 10.11-D)

Open projection prototype only when all satisfied:

```text
1. core divA validated (pairwise gradDen ~1e-7)
2. boundary divA explained or ghost/buffer strategy
3. B/E/J/Joule exact validation passed
4. η=0.1 no unacceptable deviation on exact observable
5. A-penalty cannot reduce boundary divA without breaking observable
```

---

## 8. Related documents

- [`stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md`](../../../stage10_11_boundary_divA_B_E_H_projection_eta_plan_for_cursor.md)  
- [`CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md`](CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md)
