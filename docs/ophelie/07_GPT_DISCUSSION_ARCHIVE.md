# GPT Discussion Archive (Merged Summaries)

> All former `discussion_bundle/*.md` GPT packs merged here as English index + one-line conclusions.  
> Full verbatim discussion docs preserved in [`archive/discussion/`](archive/discussion/).

---

## 1. Complex edge-flux mainline (French reload)

**Primary doc:** [`archive/discussion/OPHELIE_GPT_DISCUSSION_SINCE_COMPLEX_EDGE_FLUX_PLAN.md`](archive/discussion/OPHELIE_GPT_DISCUSSION_SINCE_COMPLEX_EDGE_FLUX_PLAN.md)

**Summary (2026-06):**

- Stage 3–4: complex φ, E, J fields on French reload
- 14 structured GPT questions on production closure
- Baseline plan: [`archive/plans/OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md`](archive/plans/OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md)
- Run record: [`archive/discussion/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md`](archive/discussion/OPHELIE_COMPLEX_EDGE_FLUX_RUN_RECORD.md)

**Deferred blockers:** [`archive/discussion/OPHELIE_GPT_DISCUSSION_BLOCKERS.md`](archive/discussion/OPHELIE_GPT_DISCUSSION_BLOCKERS.md) (B1–B5)

---

## 2. Edge-flux round summary (2026-06-08)

**Doc:** [`archive/discussion/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md`](archive/discussion/OPHELIE_GPT_ROUND_SUMMARY_2026-06-08.md)

**One-line conclusion:** Edge φ/residual succeeded; graph energy was wrongly used for 50 kW calibration; switching to **`P_recon`** made H production pass.

**Stage 2 work record:** [`archive/discussion/OPHELIE_STAGE2_CURSOR_WORK_RECORD.md`](archive/discussion/OPHELIE_STAGE2_CURSOR_WORK_RECORD.md)

---

## 3. TEAM7 L2 edge-flux

| Round | Archive | Key topics |
|-------|---------|------------|
| Round 2 | `OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md` | Initial L2 probe mismatch |
| Round 3 | `OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md` | P0–P4, solver-local audit, 12 GPT questions |
| Closure | `OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` | Round 3 decisions |

**First prompts (copy-paste):**

- `discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt`
- `discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND3.txt`
- `discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND3_CLOSURE.txt`

**Upload manifests:** `TEAM7_L2_UPLOAD_MANIFEST.md`, `TEAM7_L2_UPLOAD_MANIFEST_ROUND3.md` → moved to `archive/discussion/`

---

## 4. TEAM7 P5 no-flux boundary (2026-06-14)

**Primary doc:** [`archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_P5_NO_FLUX_BOUNDARY.md`](archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_P5_NO_FLUX_BOUNDARY.md)

- P5.0–P5.5 outcome table + 8 GPT decision questions
- First prompt: `discussion_bundle/TEAM7_P5_GPT_FIRST_PROMPT.txt`
- Zip bundle: `discussion_bundle/TEAM7_P5_NO_FLUX_BOUNDARY_BUNDLE.zip` (if present locally)

Duplicate copies in `discussion_bundle/team7_p5_pack/` → markdown in [`archive/discussion_packs/team7_p5/`](archive/discussion_packs/team7_p5/)

---

## 5. TEAM7 P56 lite closure

**Doc:** [`archive/discussion_packs/team7_p56/OPHELIE_GPT_DISCUSSION_TEAM7_P56_LITE_CLOSURE.md`](archive/discussion_packs/team7_p56/)

P56 lite pack outputs: `discussion_bundle/team7_p56_lite_pack/outputs/*.csv`

---

## 6. French three-way summaries (text)

English interpretation of H/A/G comparison runs:

| File | Version |
|------|---------|
| `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v3.txt` | Latest (H passed) |
| `discussion_bundle/FRENCH_THREE_WAY_SUMMARY_v2.txt` | Post power-fix |
| `discussion_bundle/FRENCH_THREE_WAY_SUMMARY.txt` | v1 pre power-fix |

See [03_EDGE_FLUX_PRODUCTION.md](03_EDGE_FLUX_PRODUCTION.md) §5.

---

## 7. φ operator GPT package (legacy)

Former [`archive/legacy_packages/PHI_GPT_DISCUSSION_PACKAGE.md`](archive/legacy_packages/PHI_GPT_DISCUSSION_PACKAGE.md) — RHS sign, DivSigmaA vs legacy flux, reload solvability. Superseded by consolidated [02_PHI_OPERATOR_AND_BOUNDARY.md](02_PHI_OPERATOR_AND_BOUNDARY.md).

---

## 8. How to start a new GPT round

1. Read [01_MASTER_DEVELOPMENT_TIMELINE.md](01_MASTER_DEVELOPMENT_TIMELINE.md) for current status
2. Copy template from latest `OPHELIE_GPT_ROUND_SUMMARY_*.md` in archive
3. Attach CSV from `discussion_bundle/team7_l2_outputs/` or RH200 `./output/`
4. State three facts upfront (see PHI package §1 pattern): what is **settled**, what is **production**, what is **deferred**
