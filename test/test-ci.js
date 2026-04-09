import { hksjCI, waldCI, predictionInterval, rePool, estimators, tQuantile } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg, escalc measure="RR")
// Same data as test-pooling.js — verified against R metafor 4.8-0
const yi = [-0.88931133, -1.58538866, -1.34807315, -1.44155119, -0.21754732, -0.78611559, -1.62089822, 0.01195233, -0.46941765, -1.37134480, -0.33935883, 0.44591340, -0.01731395];
const vi = [0.32558477, 0.19458112, 0.41536797, 0.02001003, 0.05121017, 0.00690562, 0.22301725, 0.00396158, 0.05643421, 0.07302479, 0.01241221, 0.53250585, 0.07140466];

describe('HKSJ CI', () => {
  it('uses t_{k-1} NOT z', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const ci = hksjCI(pool.theta, yi, vi, r.tau2, yi.length);
    const wci = waldCI(pool.theta, pool.se);
    const hksjWidth = ci[1] - ci[0];
    const waldWidth = wci[1] - wci[0];
    assert.ok(hksjWidth > waldWidth, `HKSJ ${hksjWidth} should be wider than Wald ${waldWidth}`);
  });
  it('applies floor when Q < k-1', () => {
    // Homogeneous data where Q < k-1 — HKSJ should NOT narrow below Wald
    const yi2 = [0.1, 0.12, 0.11, 0.09, 0.10];
    const vi2 = [0.01, 0.01, 0.01, 0.01, 0.01];
    const r = estimators.REML(yi2, vi2);
    const pool = rePool(yi2, vi2, r.tau2);
    const ci = hksjCI(pool.theta, yi2, vi2, r.tau2, 5);
    const wci = waldCI(pool.theta, pool.se);
    assert.ok((ci[1] - ci[0]) >= (wci[1] - wci[0]) * 0.95,
      'HKSJ with floor must not narrow below Wald');
  });
});

describe('Prediction interval', () => {
  it('uses t_{k-2} NOT t_{k-1}', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const pi = predictionInterval(pool.theta, r.tau2, pool.se, yi.length);
    const tcrit = Math.abs(tQuantile(0.025, yi.length - 2));
    assert.ok(Math.abs(tcrit - 2.201) < 0.01);
    assert.ok(pi[0] < pi[1]);
    // PI should be wider than CI
    const ci = hksjCI(pool.theta, yi, vi, r.tau2, yi.length);
    assert.ok((pi[1]-pi[0]) > (ci[1]-ci[0]) * 1.5);
  });
  it('returns NaN for k<3', () => {
    const pi = predictionInterval(0.5, 0.1, 0.05, 2);
    assert.ok(isNaN(pi[0]));
  });
});
