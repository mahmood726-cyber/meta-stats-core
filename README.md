# meta-stats-core

Journal-grade meta-analysis statistics engine in pure JavaScript. Single importable module, zero dependencies.

Validated against R metafor 4.8-0 to 6 decimal places.

## What it does

- **Pooling**: REML (Fisher scoring), Paule-Mandel, DerSimonian-Laird
- **CIs**: HKSJ with t_{k-1} and floor correction (prevents narrowing below Wald when Q < k-1)
- **Prediction intervals**: t_{k-2} (undefined for k < 3)
- **Heterogeneity**: Q, I2 (tau2-based, matching metafor), H2, Q-profile CI (Viechtbauer 2007)
- **Publication bias**: Radial Egger, Peters' test, PET-PEESE (conditional), trim-and-fill (sensitivity only)
- **Convenience**: `metaAnalysis(yi, vi)` returns everything for a forest plot + report

## Quick start

```javascript
import { metaAnalysis } from './meta-stats-core.js';

const yi = [-0.89, -1.59, -1.35, -1.44, -0.22];  // log-RR
const vi = [0.33, 0.19, 0.42, 0.02, 0.05];        // variances

const result = metaAnalysis(yi, vi);
// result.theta     = pooled effect (REML)
// result.ci        = [lower, upper] (HKSJ)
// result.pi        = [lower, upper] (prediction interval)
// result.tau2      = between-study variance
// result.I2        = heterogeneity (%)
// result.I2_ci     = [lower, upper] (Q-profile)
// result.weights   = study weights (%)
```

## For single-file HTML apps

Since single-file apps cannot use ES imports, inline the functions directly into the `<script>` section. The module has no dependencies and no DOM references.

## Design rules

| Rule | Implementation |
|------|---------------|
| Never DL for k<10 | REML default, DL available as legacy |
| HKSJ uses t, not z | `tQuantile(0.025, k-1)` |
| HKSJ floor | `max(hksjVar, waldVar)` |
| PI uses t_{k-2} | `tQuantile(0.025, k-2)`, NaN for k<3 |
| I2 != magnitude | tau2 always reported alongside |
| Trim-fill = sensitivity | `isSensitivityAnalysis: true` |
| PET-PEESE conditional | PET first; PEESE only if PET p<0.10 |
| Peters for binary | Regress on 1/n, not 1/SE |

## Tests

```bash
node --test test/*.js
```

45 tests across 7 suites. Reference values from R metafor 4.8-0 (BCG vaccine dataset).

## Stats

- Engine: pure JavaScript, single module
- Tests: 45 tests across 7 suites
- Zero dependencies
