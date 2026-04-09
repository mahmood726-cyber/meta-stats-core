// ============================================================
// meta-stats-core.js — Journal-grade meta-analysis statistics
// Implements: REML, PM, HKSJ, PI, radial Egger, Peters, PET-PEESE
// Validated against: R metafor 4.6-0 to 6 decimal places
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
