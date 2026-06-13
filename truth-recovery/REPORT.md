# Truth-Recovery Yardstick for meta-stats-core

meta-stats-core is validated against **R metafor to 6 decimals** for its point
outputs. That proves each formula is computed correctly. It does **not** tell you
which of the library's methods actually **recovers the truth** when a meta-
analysis is distorted by heterogeneity *and* publication selection together.

This module adds that missing evidence, reusing the allmeta truth-recovery
yardstick (a known-truth DGP that layers a parameterised publication-selection
mechanism on top of a known `mu` and `tau^2`). For each simulated meta-analysis
we measure how often each method's interval covers the **true** `mu`.

> Truth-first: every number is produced by `harness.mjs` from a seeded run
> (`seed=20260613`). Reproduce: `node truth-recovery/harness.mjs --reps 400`.

**Grid:** `mu=0.3`, `tau2=0.05`, `k in {5,10,15,25}`, 5 selection scenarios,
400 reps/cell.

## Mean coverage of the true `mu` over all selection cells

| method | mean coverage | mean \|bias\| | mean RMSE |
|---|---|---|---|
| DL + Wald (legacy) | 0.726 | 0.085 | 0.140 |
| REML + Wald | 0.723 | 0.085 | 0.140 |
| **REML + HKSJ (library default)** | **0.795** | 0.085 | 0.140 |
| PM + HKSJ | 0.788 | 0.086 | 0.141 |
| PET-PEESE (point only) | n/a | **0.036** | 0.214 |

## What this measures (two concrete results)

1. **The library's default is the right default — by +7 coverage points.**
   `metaAnalysis()` defaults to REML + HKSJ. The yardstick confirms that with a
   number: REML+HKSJ recovers the true `mu` in **0.795** of meta-analyses vs
   **0.726** for the legacy DL+Wald, at *identical* bias (0.085) — i.e. the gain
   is honest interval width (t_{k-1} HKSJ), not a point-estimate change. On clean
   (no-selection) data REML+HKSJ sits at ~0.95 (nominal) across k while DL+Wald
   under-covers at 0.90–0.93. (`test-truth-recovery.mjs` asserts this invariant.)

2. **Honest limitation: no inverse-variance interval recovers truth under strong
   selection — and the one method that reduces the bias has no CI.** Coverage of
   every IV method collapses as selection strengthens (`step_strong`: 0.85 at
   k=5 → **0.05** at k=25) because none of them model selection. PET-PEESE is the
   only method whose point estimate meaningfully cuts the selection bias
   (mean |bias| **0.036** vs 0.085), yet the library returns no interval for it.

## Actionable recommendation (measured, not yet implemented here)

Expose a **PET-PEESE confidence interval** (the WLS intercept already carries
`seIntercept` inside `petPeese`/`wls`) so users have a selection-robust inference
option, not just a bias-reduced point. That is the single highest-value method
addition the yardstick identifies for this library. Left as a follow-up so this
branch stays a pure, additive *measurement* layer with no change to the audited
estimators.

## Files

| File | Role |
|---|---|
| `dgp.mjs` | standalone seeded known-truth DGP (joint `tau^2` + selection) |
| `harness.mjs` | scores each method on coverage / bias / RMSE / width / type-I |
| `test-truth-recovery.mjs` | DGP reproducibility + measured invariants (6 tests) |

## What transferred from the allmeta estimator work, and what did not

- **Transferred:** the known-truth simulation harness and the joint
  heterogeneity-plus-selection DGP — they drop straight onto this library's own
  methods and produce a measured leaderboard it previously lacked.
- **Did not transfer:** the NPE / normalizing-flow estimator and conformal/SBC
  calibration are a *new method* (amortized SBI) that would be a large addition,
  not a tweak; they are out of scope for a deterministic formula library. The
  yardstick, however, quantifies exactly the gap such a method would fill (the
  strong-selection coverage collapse), which is the honest contribution here.
