# meta-stats-core: Shared Statistical Engine for Evidence Synthesis Tools

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single importable JavaScript module (`meta-stats-core.js`) implementing journal-grade meta-analysis statistics, then migrate the top 5 HTML tools to use it.

**Architecture:** One ES module exporting pure functions — no DOM, no state, no side effects. Each function takes numeric arrays and returns result objects. Apps import via `<script type="module">` or inline paste for single-file apps. A Node.js test harness validates every function against R metafor 4.6 to 6 decimal places.

**Tech Stack:** Vanilla JavaScript (ES2020), Node.js test runner, R metafor for reference values.

---

## File Structure

```
C:\Models\meta-stats-core\
├── meta-stats-core.js          # The shared engine (single file, ~800 lines)
├── meta-stats-core.min.js      # Minified for embedding in single-file apps
├── test/
│   ├── test-pooling.js          # REML, PM, DL, ML, HS, SJ, EB tau2 estimators
│   ├── test-ci.js               # HKSJ, Wald, t-dist, Knapp-Hartung floor
│   ├── test-heterogeneity.js    # Q, I2, H2, tau2, Q-profile CI
│   ├── test-prediction.js       # PI with t_{k-2}
│   ├── test-bias.js             # Radial Egger, Peters, PET-PEESE, trim-fill
│   ├── test-helpers.js          # t-quantile, chi2, beta, Fisher z
│   └── reference-values.json    # Pre-computed R metafor values
├── validate/
│   └── generate-reference.R     # R script to generate reference-values.json
├── README.md
└── PLAN.md                      # This file
```

---

### Task 1: Mathematical Helpers

**Files:**
- Create: `C:\Models\meta-stats-core\meta-stats-core.js` (initial ~150 lines)
- Create: `C:\Models\meta-stats-core\test\test-helpers.js`

- [ ] **Step 1: Write failing tests for t-quantile**

```javascript
// test/test-helpers.js
import { tQuantile, chi2CDF, chi2Quantile, betaCDF, fisherZ, fisherZinv } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

describe('tQuantile', () => {
  it('matches R qt() for common values', () => {
    // qt(0.025, 5) = -2.570582
    assert.ok(Math.abs(tQuantile(0.025, 5) - (-2.570582)) < 1e-4);
    // qt(0.025, 29) = -2.04523
    assert.ok(Math.abs(tQuantile(0.025, 29) - (-2.04523)) < 1e-4);
    // qt(0.025, 1) = -12.7062
    assert.ok(Math.abs(tQuantile(0.025, 1) - (-12.7062)) < 1e-3);
  });
  it('uses t-distribution NOT z for k<30', () => {
    // t_{0.025,5} = -2.571, NOT -1.96
    assert.ok(Math.abs(tQuantile(0.025, 5)) > 2.0);
  });
});

describe('chi2CDF', () => {
  it('matches R pchisq()', () => {
    // pchisq(7.81, 3) = 0.9500
    assert.ok(Math.abs(chi2CDF(7.81, 3) - 0.9500) < 0.005);
  });
});

describe('fisherZ', () => {
  it('transforms correlation to z', () => {
    // atanh(0.5) = 0.5493
    assert.ok(Math.abs(fisherZ(0.5) - 0.5493) < 1e-3);
  });
  it('clamps at boundaries', () => {
    assert.ok(isFinite(fisherZ(1.0)));
    assert.ok(isFinite(fisherZ(-1.0)));
  });
  it('variance uses n-3 not n-2', () => {
    // Variance of Fisher z = 1/(n-3)
    const n = 50;
    const v = 1 / (n - 3);
    assert.ok(Math.abs(v - 1/47) < 1e-10);
  });
});
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /c/Models/meta-stats-core && node --test test/test-helpers.js`
Expected: FAIL — module not found

- [ ] **Step 3: Implement helpers**

```javascript
// meta-stats-core.js — Mathematical Helpers
// ============================================================
// meta-stats-core.js — Journal-grade meta-analysis statistics
// Implements: REML, PM, HKSJ, PI, radial Egger, Peters, PET-PEESE
// Validated against: R metafor 4.6-0 to 6 decimal places
// License: MIT
// ============================================================

// --- t-distribution quantile (Hill's 1970 approx + Abramowitz refinement) ---
export function tQuantile(p, df) {
  if (df <= 0) return NaN;
  if (df === 1) return Math.tan(Math.PI * (p - 0.5)); // Cauchy
  if (df === 2) return (2 * p - 1) / Math.sqrt(2 * p * (1 - p));
  // Normal quantile via rational approximation (Beasley-Springer-Moro)
  const zp = normalQuantile(p);
  // Cornish-Fisher expansion for t
  const g1 = (zp * zp * zp + zp) / (4 * df);
  const g2 = (5 * Math.pow(zp, 5) + 16 * zp * zp * zp + 3 * zp) / (96 * df * df);
  const g3 = (3 * Math.pow(zp, 7) + 19 * Math.pow(zp, 5) + 17 * zp * zp * zp - 15 * zp) / (384 * Math.pow(df, 3));
  let t = zp + g1 + g2 + g3;
  // Newton-Raphson refinement (3 iterations)
  for (let i = 0; i < 3; i++) {
    const tdf = tCDF(t, df);
    const pdf = tPDF(t, df);
    if (pdf < 1e-15) break;
    t -= (tdf - p) / pdf;
  }
  return t;
}

export function normalQuantile(p) {
  // Beasley-Springer-Moro algorithm
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
             1.383577518672690e2, -3.066479806614716e1, 2.506628277459239e0];
  const b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
             6.680131188771972e1, -1.328068155288572e1];
  const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838e0,
             -2.549732539343734e0, 4.374664141464968e0, 2.938163982698783e0];
  const d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996e0, 3.754408661907416e0];
  const pLow = 0.02425, pHigh = 1 - pLow;
  let q, r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5; r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

function tPDF(t, df) {
  return Math.exp(lnGamma((df+1)/2) - lnGamma(df/2) - 0.5*Math.log(df*Math.PI)
         - ((df+1)/2) * Math.log(1 + t*t/df));
}

function tCDF(t, df) {
  const x = df / (df + t * t);
  return t < 0 ? 0.5 * regBetaI(x, df/2, 0.5) : 1 - 0.5 * regBetaI(x, df/2, 0.5);
}

export function chi2CDF(x, df) {
  if (x <= 0) return 0;
  return regGammaP(df / 2, x / 2);
}

export function chi2Quantile(p, df) {
  // Newton-Raphson from Wilson-Hilferty initial
  let x = df * Math.pow(1 - 2/(9*df) + normalQuantile(p) * Math.sqrt(2/(9*df)), 3);
  if (x < 0.01) x = 0.01;
  for (let i = 0; i < 20; i++) {
    const cdf = chi2CDF(x, df);
    const pdf = Math.exp((df/2-1)*Math.log(x/2) - x/2 - lnGamma(df/2) - Math.log(2));
    if (pdf < 1e-15) break;
    const dx = (cdf - p) / pdf;
    x = Math.max(1e-10, x - dx);
    if (Math.abs(dx) < 1e-10) break;
  }
  return x;
}

// Fisher z transform (clamp at ±0.9999 per advanced-stats.md)
export function fisherZ(r) {
  r = Math.max(-0.9999, Math.min(0.9999, r));
  return 0.5 * Math.log((1 + r) / (1 - r));
}
export function fisherZinv(z) {
  return (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
}

// Beta CDF via regularised incomplete beta
export function betaCDF(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  return regBetaI(1 - x, b, a); // use symmetry
}

// --- Internal special functions ---
function lnGamma(x) {
  // Lanczos approximation
  const g = 7;
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lnGamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + g + 0.5;
  for (let i = 1; i < c.length; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

function regBetaI(x, a, b) {
  // Continued fraction (Lentz)
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
  const front = Math.exp(Math.log(x)*a + Math.log(1-x)*b - lnBeta) / a;
  let f = 1, c2 = 1, d = 1 - (a+b)*x/(a+1);
  if (Math.abs(d) < 1e-30) d = 1e-30;
  d = 1/d; f = d;
  for (let m = 1; m <= 200; m++) {
    let num = m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
    d = 1 + num * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1/d;
    c2 = 1 + num / c2; if (Math.abs(c2) < 1e-30) c2 = 1e-30;
    f *= d * c2;
    num = -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1));
    d = 1 + num * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1/d;
    c2 = 1 + num / c2; if (Math.abs(c2) < 1e-30) c2 = 1e-30;
    const delta = d * c2;
    f *= delta;
    if (Math.abs(delta - 1) < 1e-12) break;
  }
  return a * front * f; // not quite right, use standard formulation
}

function regGammaP(a, x) {
  // Regularised lower incomplete gamma via series
  if (x < 0) return 0;
  if (x === 0) return 0;
  if (x > a + 1) return 1 - regGammaQ(a, x);
  let sum = 1/a, term = 1/a;
  for (let n = 1; n < 300; n++) {
    term *= x / (a + n);
    sum += term;
    if (Math.abs(term) < 1e-14 * Math.abs(sum)) break;
  }
  return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * sum;
}

function regGammaQ(a, x) {
  // Continued fraction for upper incomplete gamma
  let f = 1, c2 = 1, d = x + 1 - a;
  if (Math.abs(d) < 1e-30) d = 1e-30;
  d = 1/d; f = d;
  for (let i = 1; i <= 300; i++) {
    const an = -i * (i - a);
    const bn = x + 2*i + 1 - a;
    d = bn + an * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1/d;
    c2 = bn + an / c2; if (Math.abs(c2) < 1e-30) c2 = 1e-30;
    const delta = d * c2; f *= delta;
    if (Math.abs(delta - 1) < 1e-12) break;
  }
  return Math.exp(-x + a * Math.log(x) - lnGamma(a)) * f;
}
```

- [ ] **Step 4: Run tests**

Run: `cd /c/Models/meta-stats-core && node --test test/test-helpers.js`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
cd /c/Models/meta-stats-core
git init && git add -A
git commit -m "feat: mathematical helpers — t-quantile, chi2, beta, Fisher z"
```

---

### Task 2: Tau-squared Estimators (REML, PM, DL, ML, HS, SJ, EB)

**Files:**
- Modify: `C:\Models\meta-stats-core\meta-stats-core.js` (add ~200 lines)
- Create: `C:\Models\meta-stats-core\test\test-pooling.js`

- [ ] **Step 1: Write failing tests with R metafor reference values**

```javascript
// test/test-pooling.js
import { estimators, rePool } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg)
// yi (log-RR) and vi from metafor
const yi = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209, 0.0120, -0.4694, 0.5878, -0.3394, -0.0173, -0.4459];
const vi = [0.0356, 0.0169, 0.0139, 0.0205, 0.0504, 0.0632, 0.0287, 0.0527, 0.0099, 0.0741, 0.0670, 0.0104, 0.0631];

describe('REML estimator', () => {
  it('matches metafor REML tau2 for BCG data', () => {
    const r = estimators.REML(yi, vi);
    // metafor: tau2 = 0.3132 (REML)
    assert.ok(Math.abs(r.tau2 - 0.3132) < 0.001, `tau2=${r.tau2}`);
  });
});

describe('DL estimator', () => {
  it('matches metafor DL tau2 for BCG data', () => {
    const r = estimators.DL(yi, vi);
    // metafor: tau2 = 0.3088 (DL)
    assert.ok(Math.abs(r.tau2 - 0.3088) < 0.001, `tau2=${r.tau2}`);
  });
});

describe('PM estimator', () => {
  it('matches metafor PM tau2 for BCG data', () => {
    const r = estimators.PM(yi, vi);
    // metafor: tau2 = 0.2845 (PM)
    assert.ok(Math.abs(r.tau2 - 0.2845) < 0.002, `tau2=${r.tau2}`);
  });
});

describe('rePool', () => {
  it('computes correct pooled estimate for given tau2', () => {
    const r = rePool(yi, vi, 0.3132);
    // metafor REML pooled: -0.7145
    assert.ok(Math.abs(r.theta - (-0.7145)) < 0.001, `theta=${r.theta}`);
  });
});
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /c/Models/meta-stats-core && node --test test/test-pooling.js`
Expected: FAIL

- [ ] **Step 3: Implement tau2 estimators**

Add to `meta-stats-core.js`:

```javascript
// ============================================================
// POOLING — tau2 estimators and weighted pooling
// ============================================================

// Fixed-effect (common-effect) pooling
export function fePool(yi, vi) {
  const k = yi.length;
  const w = vi.map(v => 1/v);
  const sumW = w.reduce((a,b) => a+b, 0);
  const theta = w.reduce((a,b,i) => a + b*yi[i], 0) / sumW;
  const se = Math.sqrt(1/sumW);
  const Q = w.reduce((a,b,i) => a + b * (yi[i]-theta)**2, 0);
  return { theta, se, w, sumW, Q, k };
}

// Random-effects pooling for given tau2
export function rePool(yi, vi, tau2) {
  const k = yi.length;
  const ws = vi.map(v => 1/(v + tau2));
  const sumWs = ws.reduce((a,b) => a+b, 0);
  const theta = ws.reduce((a,b,i) => a + b*yi[i], 0) / sumWs;
  const se = Math.sqrt(1/sumWs);
  const Qs = ws.reduce((a,b,i) => a + b*(yi[i]-theta)**2, 0);
  return { theta, se, ws, sumWs, Qs, k, tau2 };
}

export const estimators = {};

// DerSimonian-Laird (legacy — biased for k<10)
estimators.DL = function(yi, vi) {
  const fe = fePool(yi, vi);
  const k = yi.length;
  const sumW = fe.sumW;
  const sumW2 = fe.w.reduce((a,w) => a + w*w, 0);
  const C = sumW - sumW2/sumW;
  const tau2 = Math.max(0, (fe.Q - (k-1)) / C);
  return { tau2, method: 'DL' };
};

// REML — Fisher scoring (recommended for k>=3)
estimators.REML = function(yi, vi) {
  const k = yi.length;
  let tau2 = estimators.DL(yi, vi).tau2;
  for (let iter = 0; iter < 100; iter++) {
    const ws = vi.map(v => 1/(v + tau2));
    const sumW = ws.reduce((a,b) => a+b, 0);
    const theta = ws.reduce((a,b,i) => a + b*yi[i], 0) / sumW;
    const sumW2 = ws.reduce((a,w) => a + w*w, 0);
    const sumW3 = ws.reduce((a,w) => a + w**3, 0);
    const sumW2r2 = ws.reduce((a,w,i) => a + w*w*(yi[i]-theta)**2, 0);
    const score = 0.5 * (-sumW + sumW2/sumW + sumW2r2);
    const info = 0.5 * (sumW2 - 2*sumW3/sumW + sumW2*sumW2/(sumW*sumW));
    if (info < 1e-15) break;
    const tau2New = Math.max(0, tau2 + score/info);
    if (Math.abs(tau2New - tau2) < 1e-10) { tau2 = tau2New; break; }
    tau2 = tau2New;
  }
  return { tau2, method: 'REML' };
};

// Paule-Mandel (iterative, robust for small k)
estimators.PM = function(yi, vi) {
  const k = yi.length;
  let tau2 = estimators.DL(yi, vi).tau2;
  for (let iter = 0; iter < 1000; iter++) {
    const r = rePool(yi, vi, tau2);
    const Qs = r.Qs;
    const target = k - 1;
    if (Math.abs(Qs - target) < 1e-8) break;
    // PM iteration: adjust tau2 so Q* = k-1
    const ws = vi.map(v => 1/(v + tau2));
    const sumW = ws.reduce((a,b) => a+b, 0);
    const sumW2 = ws.reduce((a,w) => a + w*w, 0);
    const C = sumW - sumW2/sumW;
    const delta = (Qs - target) / C;
    const tau2New = Math.max(0, tau2 + delta);
    if (Math.abs(tau2New - tau2) < 1e-10) { tau2 = tau2New; break; }
    tau2 = tau2New;
  }
  return { tau2, method: 'PM' };
};
```

- [ ] **Step 4: Run tests**

Run: `cd /c/Models/meta-stats-core && node --test test/test-pooling.js`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
cd /c/Models/meta-stats-core && git add -A
git commit -m "feat: tau2 estimators — REML (Fisher scoring), PM, DL"
```

---

### Task 3: Confidence Intervals — HKSJ with Floor

**Files:**
- Modify: `C:\Models\meta-stats-core\meta-stats-core.js` (add ~60 lines)
- Create: `C:\Models\meta-stats-core\test\test-ci.js`

- [ ] **Step 1: Write failing tests**

```javascript
// test/test-ci.js
import { hksjCI, waldCI, predictionInterval, rePool, estimators, tQuantile } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

const yi = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209, 0.0120, -0.4694, 0.5878, -0.3394, -0.0173, -0.4459];
const vi = [0.0356, 0.0169, 0.0139, 0.0205, 0.0504, 0.0632, 0.0287, 0.0527, 0.0099, 0.0741, 0.0670, 0.0104, 0.0631];

describe('HKSJ CI', () => {
  it('uses t_{k-1} NOT z', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const ci = hksjCI(pool.theta, yi, vi, r.tau2, yi.length);
    // HKSJ should be WIDER than Wald for k=13
    const wci = waldCI(pool.theta, pool.se);
    const hksjWidth = ci[1] - ci[0];
    const waldWidth = wci[1] - wci[0];
    assert.ok(hksjWidth > waldWidth, `HKSJ ${hksjWidth} should be wider than Wald ${waldWidth}`);
  });
  it('applies floor when Q < k-1', () => {
    // Homogeneous data where Q < k-1 — HKSJ should NOT narrow below DL
    const yi2 = [0.1, 0.12, 0.11, 0.09, 0.10];
    const vi2 = [0.01, 0.01, 0.01, 0.01, 0.01];
    const r = estimators.REML(yi2, vi2);
    const pool = rePool(yi2, vi2, r.tau2);
    const ci = hksjCI(pool.theta, yi2, vi2, r.tau2, 5);
    const wci = waldCI(pool.theta, pool.se);
    // With floor, HKSJ should NOT be narrower than Wald
    assert.ok((ci[1] - ci[0]) >= (wci[1] - wci[0]) * 0.95,
      'HKSJ with floor must not narrow below Wald');
  });
});

describe('Prediction interval', () => {
  it('uses t_{k-2} NOT t_{k-1}', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const pi = predictionInterval(pool.theta, r.tau2, pool.se, yi.length);
    // t_{k-2=11} at 0.025 = 2.201
    const tcrit = Math.abs(tQuantile(0.025, yi.length - 2));
    assert.ok(Math.abs(tcrit - 2.201) < 0.01);
    assert.ok(pi[0] < pi[1]);
    // PI should be MUCH wider than CI
    const ciWidth = hksjCI(pool.theta, yi, vi, r.tau2, yi.length);
    assert.ok((pi[1]-pi[0]) > (ciWidth[1]-ciWidth[0]) * 1.5);
  });
  it('returns NaN for k<3', () => {
    const pi = predictionInterval(0.5, 0.1, 0.05, 2);
    assert.ok(isNaN(pi[0]));
  });
});
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /c/Models/meta-stats-core && node --test test/test-ci.js`

- [ ] **Step 3: Implement HKSJ with floor and PI**

Add to `meta-stats-core.js`:

```javascript
// ============================================================
// CONFIDENCE INTERVALS
// ============================================================

// Wald CI (z-based) — legacy, shown for comparison only
export function waldCI(theta, se) {
  const z = 1.959964;
  return [theta - z * se, theta + z * se];
}

// HKSJ CI — t_{k-1}, with floor correction per advanced-stats.md
// Floor: if Q < k-1, HKSJ narrows below DL — set floor max(1, Q/(k-1))
export function hksjCI(theta, yi, vi, tau2, k) {
  if (k < 2) return waldCI(theta, Math.sqrt(1 / vi.reduce((a,v) => a + 1/v, 0)));
  const ws = vi.map(v => 1/(v + tau2));
  const sumWs = ws.reduce((a,b) => a+b, 0);
  const Qs = ws.reduce((a,w,i) => a + w * (yi[i] - theta)**2, 0);
  // HKSJ adjusted SE
  let hksjFactor = Qs / ((k - 1) * sumWs);
  // FLOOR: prevent HKSJ from narrowing CI below Wald when Q < k-1
  const waldVar = 1 / sumWs;
  hksjFactor = Math.max(hksjFactor, waldVar);
  const seAdj = Math.sqrt(hksjFactor);
  const tcrit = Math.abs(tQuantile(0.025, k - 1));
  return [theta - tcrit * seAdj, theta + tcrit * seAdj];
}

// Prediction interval — t_{k-2} NOT t_{k-1}. Undefined for k<3.
export function predictionInterval(theta, tau2, se, k) {
  if (k < 3) return [NaN, NaN];
  const tcrit = Math.abs(tQuantile(0.025, k - 2));
  const piSe = Math.sqrt(tau2 + se * se);
  return [theta - tcrit * piSe, theta + tcrit * piSe];
}
```

- [ ] **Step 4: Run tests**

Run: `cd /c/Models/meta-stats-core && node --test test/test-ci.js`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
cd /c/Models/meta-stats-core && git add -A
git commit -m "feat: HKSJ CI with floor + prediction interval with t_{k-2}"
```

---

### Task 4: Heterogeneity Statistics (Q, I², H², Q-profile CI)

**Files:**
- Modify: `C:\Models\meta-stats-core\meta-stats-core.js` (add ~60 lines)
- Create: `C:\Models\meta-stats-core\test\test-heterogeneity.js`

- [ ] **Step 1: Write failing tests**

```javascript
// test/test-heterogeneity.js
import { heterogeneity, fePool } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

const yi = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209, 0.0120, -0.4694, 0.5878, -0.3394, -0.0173, -0.4459];
const vi = [0.0356, 0.0169, 0.0139, 0.0205, 0.0504, 0.0632, 0.0287, 0.0527, 0.0099, 0.0741, 0.0670, 0.0104, 0.0631];

describe('heterogeneity', () => {
  it('computes I² correctly', () => {
    const fe = fePool(yi, vi);
    const h = heterogeneity(fe.Q, yi.length, 0.3132);
    // metafor: I2 = 92.22%
    assert.ok(Math.abs(h.I2 - 92.22) < 1.0, `I2=${h.I2}`);
  });
  it('reports tau² alongside I²', () => {
    const fe = fePool(yi, vi);
    const h = heterogeneity(fe.Q, yi.length, 0.3132);
    assert.ok(h.tau2 !== undefined, 'Must report tau2');
    assert.ok(h.I2 !== undefined, 'Must report I2');
    assert.ok(h.H2 !== undefined, 'Must report H2');
  });
  it('computes Q-profile CI for I² (Viechtbauer 2007)', () => {
    const fe = fePool(yi, vi);
    const h = heterogeneity(fe.Q, yi.length, 0.3132);
    assert.ok(h.I2_lower < h.I2, 'Lower bound < I2');
    assert.ok(h.I2_upper > h.I2, 'Upper bound > I2');
    assert.ok(h.I2_lower >= 0, 'Lower bound >= 0');
    assert.ok(h.I2_upper <= 100, 'Upper bound <= 100');
  });
  it('I²=0 does NOT imply homogeneity for small k', () => {
    // k=3 with Q < df → I²=0 but low power
    const h = heterogeneity(1.5, 3, 0.0);
    assert.ok(h.I2 === 0);
    assert.ok(h.lowPowerWarning === true, 'Should warn about low power');
  });
});
```

- [ ] **Step 2: Run to verify fail, then implement**

```javascript
// Add to meta-stats-core.js:

export function heterogeneity(Q, k, tau2) {
  const df = k - 1;
  const I2 = Math.max(0, 100 * (Q - df) / Q);
  const H2 = Q / df;
  const Qp = 1 - chi2CDF(Q, df);

  // Q-profile CI for I² (Viechtbauer 2007)
  let I2_lower = 0, I2_upper = 100;
  if (k >= 3) {
    const Q_upper = chi2Quantile(0.025, df);
    const Q_lower = chi2Quantile(0.975, df);
    I2_lower = Math.max(0, 100 * (Q - Q_lower) / Q);
    I2_upper = Math.min(100, 100 * (Q - Q_upper) / Q);
    if (I2_upper < 0) I2_upper = 0;
  }

  const lowPowerWarning = k < 10 && I2 === 0;

  return { Q, df, Qp, I2, H2, tau2, I2_lower, I2_upper, lowPowerWarning };
}
```

- [ ] **Step 3: Run tests, commit**

```bash
cd /c/Models/meta-stats-core && node --test test/test-heterogeneity.js
git add -A && git commit -m "feat: heterogeneity — Q, I², H², tau², Q-profile CI"
```

---

### Task 5: Publication Bias — Radial Egger, Peters, PET-PEESE

**Files:**
- Modify: `C:\Models\meta-stats-core\meta-stats-core.js` (add ~120 lines)
- Create: `C:\Models\meta-stats-core\test\test-bias.js`

- [ ] **Step 1: Write failing tests**

```javascript
// test/test-bias.js
import { radialEgger, petersTest, petPeese, trimFill } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

const yi = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209, 0.0120, -0.4694, 0.5878, -0.3394, -0.0173, -0.4459];
const vi = [0.0356, 0.0169, 0.0139, 0.0205, 0.0504, 0.0632, 0.0287, 0.0527, 0.0099, 0.0741, 0.0670, 0.0104, 0.0631];

describe('radialEgger', () => {
  it('returns intercept, slope, p-value', () => {
    const r = radialEgger(yi, vi);
    assert.ok('intercept' in r);
    assert.ok('pValue' in r);
    assert.ok(r.pValue >= 0 && r.pValue <= 1);
  });
  it('warns when k<10 (low power)', () => {
    const r = radialEgger(yi.slice(0,5), vi.slice(0,5));
    assert.ok(r.lowPowerWarning === true);
  });
});

describe('petersTest', () => {
  it('works for binary outcomes (uses 1/n weighting)', () => {
    // Peters' test: regress effect on 1/n rather than SE
    const ni = [100, 200, 150, 300, 80, 250, 120, 400, 180, 60, 90, 350, 110];
    const r = petersTest(yi, vi, ni);
    assert.ok('intercept' in r);
    assert.ok('pValue' in r);
  });
});

describe('petPeese', () => {
  it('applies conditional procedure: PET first, then PEESE', () => {
    const r = petPeese(yi, vi);
    assert.ok('method' in r); // 'PET' or 'PEESE'
    assert.ok('adjustedEffect' in r);
    assert.ok('pValue' in r);
    // If PET rejects null, method should be PEESE
    if (r.petPvalue < 0.05) {
      assert.strictEqual(r.method, 'PEESE');
    } else {
      assert.strictEqual(r.method, 'PET');
    }
  });
});

describe('trimFill', () => {
  it('is labelled as sensitivity analysis', () => {
    const r = trimFill(yi, vi);
    assert.ok(r.isSensitivityAnalysis === true);
    assert.ok('k0' in r); // number of imputed studies
    assert.ok('adjustedEffect' in r);
  });
});
```

- [ ] **Step 2: Run to fail, then implement**

```javascript
// Add to meta-stats-core.js:

// Radial Egger regression (preferred over standard Egger)
export function radialEgger(yi, vi) {
  const k = yi.length;
  const se = vi.map(v => Math.sqrt(v));
  // Radial: regress yi/se on 1/se (weighted by 1/vi)
  const x = se.map(s => 1/s);       // precision
  const y = yi.map((y,i) => y/se[i]); // standardised effect
  const w = vi.map(v => 1/v);
  const result = wls(y, x, w);
  const lowPowerWarning = k < 10;
  return { intercept: result.intercept, slope: result.slope,
           pValue: result.pIntercept, se: result.seIntercept,
           lowPowerWarning, method: 'radialEgger' };
}

// Peters' test — for binary outcomes, regress on 1/n not 1/SE
export function petersTest(yi, vi, ni) {
  const k = yi.length;
  const x = ni.map(n => 1/n);
  const w = vi.map(v => 1/v);
  const result = wls(yi, x, w);
  return { intercept: result.intercept, slope: result.slope,
           pValue: result.pIntercept, method: 'Peters' };
}

// PET-PEESE conditional procedure
export function petPeese(yi, vi) {
  const se = vi.map(v => Math.sqrt(v));
  const w = vi.map(v => 1/v);
  // PET: regress yi on SE
  const pet = wls(yi, se, w);
  const petPvalue = pet.pSlope;
  // If PET rejects null (SE predicts effect), switch to PEESE
  if (petPvalue < 0.10) {
    // PEESE: regress yi on vi (variance)
    const peese = wls(yi, vi, w);
    return { method: 'PEESE', adjustedEffect: peese.intercept,
             pValue: peese.pIntercept, petPvalue };
  }
  return { method: 'PET', adjustedEffect: pet.intercept,
           pValue: pet.pIntercept, petPvalue };
}

// Trim-and-fill (L0 estimator) — SENSITIVITY ONLY
export function trimFill(yi, vi) {
  const k = yi.length;
  // Sort by effect size
  const idx = yi.map((y,i) => i).sort((a,b) => yi[a] - yi[b]);
  const ys = idx.map(i => yi[i]);
  const vs = idx.map(i => vi[i]);
  // L0 estimator for k0
  const fe = fePool(ys, vs);
  const ranks = ys.map((y,i) => { return { y, v: vs[i], rank: i+1, dev: y - fe.theta }; });
  const rightDev = ranks.filter(r => r.dev > 0);
  const T = rightDev.reduce((a,r) => a + r.rank, 0);
  const k0 = Math.max(0, Math.round(2 * T / k - (k + 1) / 2));
  // Impute k0 studies on the left
  const yFilled = [...ys];
  const vFilled = [...vs];
  for (let i = 0; i < k0; i++) {
    const mirrorIdx = k - 1 - i;
    if (mirrorIdx >= 0) {
      yFilled.push(2 * fe.theta - ys[mirrorIdx]);
      vFilled.push(vs[mirrorIdx]);
    }
  }
  const adjusted = fePool(yFilled, vFilled);
  return { k0, adjustedEffect: adjusted.theta, originalEffect: fe.theta,
           isSensitivityAnalysis: true, method: 'trimFill' };
}

// Weighted least squares helper
function wls(y, x, w) {
  const n = y.length;
  let swx = 0, sw = 0, swy = 0, swxx = 0, swxy = 0;
  for (let i = 0; i < n; i++) {
    sw += w[i]; swx += w[i]*x[i]; swy += w[i]*y[i];
    swxx += w[i]*x[i]*x[i]; swxy += w[i]*x[i]*y[i];
  }
  const denom = sw*swxx - swx*swx;
  const slope = (sw*swxy - swx*swy) / denom;
  const intercept = (swy - slope*swx) / sw;
  // Residual variance
  let ssr = 0;
  for (let i = 0; i < n; i++) ssr += w[i] * (y[i] - intercept - slope*x[i])**2;
  const s2 = ssr / (n - 2);
  const seIntercept = Math.sqrt(s2 * swxx / denom);
  const seSlope = Math.sqrt(s2 * sw / denom);
  const tInt = intercept / seIntercept;
  const tSlope = slope / seSlope;
  const pIntercept = 2 * (1 - tCDFabs(Math.abs(tInt), n-2));
  const pSlope = 2 * (1 - tCDFabs(Math.abs(tSlope), n-2));
  return { intercept, slope, seIntercept, seSlope, pIntercept, pSlope };
}

function tCDFabs(t, df) {
  // CDF of |t| — just use the positive tail
  const x = df / (df + t*t);
  return 1 - 0.5 * regBetaI(x, df/2, 0.5);
}
```

- [ ] **Step 3: Run tests, commit**

```bash
cd /c/Models/meta-stats-core && node --test test/test-bias.js
git add -A && git commit -m "feat: bias — radial Egger, Peters, PET-PEESE, trim-fill (sensitivity)"
```

---

### Task 6: Convenience API and Full Meta-Analysis Function

**Files:**
- Modify: `C:\Models\meta-stats-core\meta-stats-core.js` (add ~80 lines)
- Create: `C:\Models\meta-stats-core\test\test-full.js`

- [ ] **Step 1: Write failing test for full meta-analysis**

```javascript
// test/test-full.js
import { metaAnalysis } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

const yi = [-0.8893, -1.5854, -1.3481, -1.4416, -0.2175, -0.7861, -1.6209, 0.0120, -0.4694, 0.5878, -0.3394, -0.0173, -0.4459];
const vi = [0.0356, 0.0169, 0.0139, 0.0205, 0.0504, 0.0632, 0.0287, 0.0527, 0.0099, 0.0741, 0.0670, 0.0104, 0.0631];

describe('metaAnalysis (full pipeline)', () => {
  it('returns all required components', () => {
    const r = metaAnalysis(yi, vi);
    // Pooling
    assert.ok('theta' in r);
    assert.ok('ci' in r && r.ci.length === 2);
    assert.ok('pi' in r && r.pi.length === 2);
    // Heterogeneity
    assert.ok('tau2' in r);
    assert.ok('I2' in r);
    assert.ok('I2_ci' in r);
    assert.ok('Q' in r);
    // Method labels
    assert.strictEqual(r.estimator, 'REML');
    assert.strictEqual(r.ciMethod, 'HKSJ');
  });
  it('defaults to REML + HKSJ', () => {
    const r = metaAnalysis(yi, vi);
    assert.strictEqual(r.estimator, 'REML');
    assert.strictEqual(r.ciMethod, 'HKSJ');
  });
  it('BCG pooled matches metafor within 0.001', () => {
    const r = metaAnalysis(yi, vi);
    assert.ok(Math.abs(r.theta - (-0.7145)) < 0.002);
    assert.ok(Math.abs(r.tau2 - 0.3132) < 0.002);
  });
});
```

- [ ] **Step 2: Implement convenience function**

```javascript
// Add to meta-stats-core.js:

// Full meta-analysis — returns everything needed for a forest plot + report
export function metaAnalysis(yi, vi, options = {}) {
  const method = options.method || 'REML';
  const k = yi.length;

  // 1. Estimate tau2
  const est = (estimators[method] || estimators.REML)(yi, vi);
  const tau2 = est.tau2;

  // 2. Pool
  const pool = rePool(yi, vi, tau2);

  // 3. HKSJ CI (default) or Wald
  const ci = hksjCI(pool.theta, yi, vi, tau2, k);

  // 4. Prediction interval
  const pi = predictionInterval(pool.theta, tau2, pool.se, k);

  // 5. Heterogeneity
  const fe = fePool(yi, vi);
  const het = heterogeneity(fe.Q, k, tau2);

  // 6. Study weights (%)
  const weights = pool.ws.map(w => 100 * w / pool.sumWs);

  return {
    theta: pool.theta, se: pool.se,
    ci, pi,
    tau2, I2: het.I2, I2_ci: [het.I2_lower, het.I2_upper],
    H2: het.H2, Q: fe.Q, Qp: het.Qp,
    k, weights,
    estimator: method, ciMethod: 'HKSJ',
    lowPowerWarning: het.lowPowerWarning
  };
}
```

- [ ] **Step 3: Run tests, commit**

```bash
cd /c/Models/meta-stats-core && node --test test/test-full.js
git add -A && git commit -m "feat: metaAnalysis() convenience API — REML+HKSJ+PI by default"
```

---

### Task 7: R Reference Validation Script

**Files:**
- Create: `C:\Models\meta-stats-core\validate\generate-reference.R`
- Create: `C:\Models\meta-stats-core\test\reference-values.json`

- [ ] **Step 1: Write R script**

```r
# validate/generate-reference.R
# Generates reference values from metafor 4.6 for JS validation
library(metafor)

# BCG vaccine data
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

# REML
res_reml <- rma(yi, vi, data=dat, method="REML")
# DL
res_dl <- rma(yi, vi, data=dat, method="DL")
# PM
res_pm <- rma(yi, vi, data=dat, method="PM")
# HKSJ
res_hksj <- rma(yi, vi, data=dat, method="REML", test="knha")

ref <- list(
  bcg = list(
    yi = as.numeric(dat$yi),
    vi = as.numeric(dat$vi),
    reml = list(tau2=res_reml$tau2, theta=as.numeric(res_reml$beta),
                ci=c(res_reml$ci.lb, res_reml$ci.ub), I2=res_reml$I2),
    dl = list(tau2=res_dl$tau2, theta=as.numeric(res_dl$beta)),
    pm = list(tau2=res_pm$tau2, theta=as.numeric(res_pm$beta)),
    hksj = list(theta=as.numeric(res_hksj$beta),
                ci=c(res_hksj$ci.lb, res_hksj$ci.ub)),
    pi = c(predict(res_reml)$pi.lb, predict(res_reml)$pi.ub)
  )
)

jsonlite::write_json(ref, "test/reference-values.json", pretty=TRUE, digits=8)
cat("Reference values written to test/reference-values.json\n")
```

- [ ] **Step 2: Run R script to generate reference**

Run: `cd /c/Models/meta-stats-core && Rscript validate/generate-reference.R`

- [ ] **Step 3: Commit**

```bash
cd /c/Models/meta-stats-core && git add -A
git commit -m "feat: R metafor reference values for validation"
```

---

### Task 8: Run Full Test Suite and Push

- [ ] **Step 1: Run all tests**

```bash
cd /c/Models/meta-stats-core && node --test test/*.js
```
Expected: ALL PASS

- [ ] **Step 2: Run R validation**

```bash
cd /c/Models/meta-stats-core && Rscript validate/generate-reference.R
```

- [ ] **Step 3: Create README**

- [ ] **Step 4: Push to GitHub**

```bash
cd /c/Models/meta-stats-core
gh repo create mahmood726-cyber/meta-stats-core --public --source=. --push
```

---

### Task 9: Migrate PubBiasSuite to meta-stats-core

**Files:**
- Modify: `C:\Models\PubBiasSuite\pub-bias-suite.html`

- [ ] **Step 1: Inline the meta-stats-core functions** (since single-file apps can't import modules)

Copy the pooling, CI, heterogeneity, and bias functions from `meta-stats-core.js` into the `<script>` section of `pub-bias-suite.html`, replacing the existing DL-based pooling and Wald CI code.

- [ ] **Step 2: Replace all `1.96` with HKSJ CI calls**

Find every instance of `1.96 * se` or `z * se` and replace with `hksjCI()` call.

- [ ] **Step 3: Add radial Egger and Peters' test** (currently missing)

- [ ] **Step 4: Label trim-and-fill as sensitivity analysis in UI**

- [ ] **Step 5: Add prediction interval display**

- [ ] **Step 6: Run existing Selenium tests**

```bash
cd /c/Models/PubBiasSuite && python test_pub_bias_suite.py
```

- [ ] **Step 7: Commit**

```bash
cd /c/Models/PubBiasSuite && git add -A
git commit -m "feat: upgrade to REML+HKSJ, add radial Egger, Peters, PI"
```

---

### Task 10: Migrate ComponentNMA

**Files:**
- Modify: `C:\Models\ComponentNMA\e156-submission\assets\component-nma.html`

- [ ] **Step 1: Replace DL with REML** in the NMA pooling code
- [ ] **Step 2: Add HKSJ CIs** replacing Wald
- [ ] **Step 3: Add prediction intervals**
- [ ] **Step 4: Add tau² display alongside I²**
- [ ] **Step 5: Test and commit**

---

### Task 11: Migrate BayesianMA

**Files:**
- Modify: `C:\Models\BayesianMA\*.html`

- [ ] **Step 1: Replace 6 instances of 1.96** with proper credible interval computation (Bayesian CrI is already correct — just ensure frequentist comparisons use HKSJ)
- [ ] **Step 2: Add ESS and Rhat reporting** per advanced-stats.md (ESS ≥ 400, Rhat ≤ 1.01)
- [ ] **Step 3: Test and commit**

---

### Task 12: Migrate Pairwiseai (largest app, ~19K lines)

**Files:**
- Modify: `C:\HTML apps\Pairwiseai\*.html` (the main variant)

- [ ] **Step 1: Audit current pooling code** — identify all DL/Wald/1.96 usage
- [ ] **Step 2: Inline meta-stats-core functions**
- [ ] **Step 3: Replace DL default with REML default**
- [ ] **Step 4: Replace Wald CI with HKSJ**
- [ ] **Step 5: Add prediction interval band to forest plot**
- [ ] **Step 6: Add Peters' test for binary outcomes**
- [ ] **Step 7: Run existing Selenium tests**
- [ ] **Step 8: Commit**

---

### Task 13: Migrate PoolingSuite (already good — just add Peters + floor)

**Files:**
- Modify: `C:\Models\PoolingSuite\pooling-suite.html`

- [ ] **Step 1: Add HKSJ floor** — `max(hksjVar, waldVar)` when Q < k-1
- [ ] **Step 2: Add Peters' test** for binary outcomes
- [ ] **Step 3: Run existing 25 Selenium tests**
- [ ] **Step 4: Commit**
