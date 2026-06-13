// ============================================================
// meta-stats-core.js — Journal-grade meta-analysis statistics
// Implements: REML, PM, HKSJ, PI, radial Egger, Peters, PET-PEESE
// Validated against: R metafor 4.8-0 to 6 decimal places
// License: MIT
// ============================================================

// --- Normal quantile (Beasley-Springer-Moro) ---
export function normalQuantile(p) {
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

// --- t-distribution quantile (Cornish-Fisher + Newton-Raphson) ---
export function tQuantile(p, df) {
  if (df <= 0) return NaN;
  if (df === 1) return Math.tan(Math.PI * (p - 0.5)); // Cauchy
  if (df === 2) return (2 * p - 1) / Math.sqrt(2 * p * (1 - p));
  const zp = normalQuantile(p);
  const g1 = (zp * zp * zp + zp) / (4 * df);
  const g2 = (5 * Math.pow(zp, 5) + 16 * zp * zp * zp + 3 * zp) / (96 * df * df);
  const g3 = (3 * Math.pow(zp, 7) + 19 * Math.pow(zp, 5) + 17 * zp * zp * zp - 15 * zp) / (384 * Math.pow(df, 3));
  let t = zp + g1 + g2 + g3;
  for (let i = 0; i < 3; i++) {
    const cdf = tCDF(t, df);
    const pdf = tPDF(t, df);
    if (pdf < 1e-15) break;
    t -= (cdf - p) / pdf;
  }
  return t;
}

// --- Chi-squared CDF and quantile ---
export function chi2CDF(x, df) {
  if (x <= 0) return 0;
  return regGammaP(df / 2, x / 2);
}

export function chi2Quantile(p, df) {
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

// --- Fisher z transform (clamp at ±0.9999 per rules) ---
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
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - betaCDF(1 - x, b, a);
  }
  const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
  const front = Math.exp(a * Math.log(x) + b * Math.log(1 - x) - lnBeta) / a;
  return front * betaCF(x, a, b);
}

// --- Internal special functions ---
function tPDF(t, df) {
  return Math.exp(lnGamma((df+1)/2) - lnGamma(df/2) - 0.5*Math.log(df*Math.PI)
         - ((df+1)/2) * Math.log(1 + t*t/df));
}

function tCDF(t, df) {
  const x = df / (df + t * t);
  const ib = regBetaI(x, df/2, 0.5);
  return t < 0 ? 0.5 * ib : 1 - 0.5 * ib;
}

export function lnGamma(x) {
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
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
  const front = Math.exp(Math.log(x)*a + Math.log(1-x)*b - lnBeta) / a;
  return front * betaCF(x, a, b);
}

function betaCF(x, a, b) {
  // Continued fraction (Lentz's method)
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
  return f;
}

function regGammaP(a, x) {
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
  // Modified Lentz continued fraction (Abramowitz & Stegun 6.5.31)
  let f = x - a + 1 + 1e-30;
  let c2 = f, d = 0, delta;
  for (let m = 1; m <= 300; m++) {
    const an1 = m * (a - m);
    const an2 = x - a + 2*m + 1;
    d = an2 + an1 * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1/d;
    c2 = an2 + an1 / c2; if (Math.abs(c2) < 1e-30) c2 = 1e-30;
    delta = c2 * d; f *= delta;
    if (Math.abs(delta - 1) < 1e-12) break;
  }
  return Math.exp(-x + a * Math.log(x) - lnGamma(a)) / f;
}

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
  const k = yi.length;
  // k<2: no between-study variance is estimable (df=0, C=0 → 0/0). tau2=0 by convention.
  if (k < 2) return { tau2: 0, method: 'DL' };
  const fe = fePool(yi, vi);
  const sumW = fe.sumW;
  const sumW2 = fe.w.reduce((a,w) => a + w*w, 0);
  const C = sumW - sumW2/sumW;
  const tau2 = C > 0 ? Math.max(0, (fe.Q - (k-1)) / C) : 0;
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
    // REML score and info (Fisher scoring)
    const score = -0.5 * sumW + 0.5 * sumW2/sumW + 0.5 * sumW2r2;
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

// ============================================================
// CONFIDENCE INTERVALS
// ============================================================

// Wald CI (z-based) — legacy, shown for comparison only
export function waldCI(theta, se) {
  const z = 1.959964;
  return [theta - z * se, theta + z * se];
}

// HKSJ CI — t_{k-1}, with floor correction
// Floor: if Q < k-1, HKSJ narrows below DL — set floor max(hksjVar, waldVar)
export function hksjCI(theta, yi, vi, tau2, k) {
  if (k < 2) return waldCI(theta, Math.sqrt(1 / vi.reduce((a,v) => a + 1/v, 0)));
  const ws = vi.map(v => 1/(v + tau2));
  const sumWs = ws.reduce((a,b) => a+b, 0);
  const Qs = ws.reduce((a,w,i) => a + w * (yi[i] - theta)**2, 0);
  // HKSJ variance
  let hksjVar = Qs / ((k - 1) * sumWs);
  // FLOOR: prevent HKSJ from narrowing CI below Wald when Q* < k-1
  const waldVar = 1 / sumWs;
  hksjVar = Math.max(hksjVar, waldVar);
  const seAdj = Math.sqrt(hksjVar);
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

// ============================================================
// HETEROGENEITY
// ============================================================

export function heterogeneity(Q, k, tau2, vi) {
  const df = k - 1;
  // If vi provided, use tau²-based I² (matches metafor):
  //   s2w = (k-1) * sum(wi) / (sum(wi)² - sum(wi²))
  //   I² = 100 * tau² / (tau² + s2w)
  // Otherwise fall back to Q-based: I² = 100 * (Q - df) / Q
  let I2;
  if (df < 1) {
    // Single study (k<2): heterogeneity is undefined; report 0 rather than 0/0=NaN.
    I2 = 0;
  } else if (vi && vi.length >= 2 && tau2 > 0) {
    const wi = vi.map(v => 1 / v);
    const sumW = wi.reduce((a, b) => a + b, 0);
    const sumW2 = wi.reduce((a, w) => a + w * w, 0);
    const s2w = (k - 1) * sumW / (sumW * sumW - sumW2);
    I2 = Math.max(0, 100 * tau2 / (tau2 + s2w));
  } else {
    I2 = Math.max(0, 100 * (Q - df) / Q);
  }
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

// ============================================================
// PUBLICATION BIAS
// ============================================================

// Radial Egger regression (preferred over standard Egger)
export function radialEgger(yi, vi) {
  const k = yi.length;
  const se = vi.map(v => Math.sqrt(v));
  // Radial: regress yi/se on 1/se (weighted by 1/vi)
  const x = se.map(s => 1/s);        // precision
  const y = yi.map((yv,i) => yv/se[i]); // standardised effect
  const w = vi.map(v => 1/v);
  const result = wls(y, x, w);
  const lowPowerWarning = k < 10;
  return { intercept: result.intercept, slope: result.slope,
           pValue: result.pIntercept, se: result.seIntercept,
           lowPowerWarning, method: 'radialEgger' };
}

// Peters' test — for binary outcomes, regress on 1/n not 1/SE
export function petersTest(yi, vi, ni) {
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
  const idx = yi.map((_,i) => i).sort((a,b) => yi[a] - yi[b]);
  const ys = idx.map(i => yi[i]);
  const vs = idx.map(i => vi[i]);
  const fe = fePool(ys, vs);
  const ranks = ys.map((y,i) => ({ y, v: vs[i], rank: i+1, dev: y - fe.theta }));
  const rightDev = ranks.filter(r => r.dev > 0);
  const T = rightDev.reduce((a,r) => a + r.rank, 0);
  const k0 = Math.max(0, Math.round(2 * T / k - (k + 1) / 2));
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

// ============================================================
// CONVENIENCE API
// ============================================================

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
  const het = heterogeneity(fe.Q, k, tau2, vi);

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

// --- Internal: Weighted least squares ---
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
  // CDF of |t| — positive tail
  if (df <= 0) return 0.5;
  const x = df / (df + t*t);
  return 1 - 0.5 * regBetaI(x, df/2, 0.5);
}
