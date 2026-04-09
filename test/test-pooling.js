import { estimators, rePool, fePool } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg, escalc measure="RR")
// Verified against R metafor 4.8-0 (8 significant figures to avoid rounding artefacts)
const yi = [-0.88931133, -1.58538866, -1.34807315, -1.44155119, -0.21754732, -0.78611559, -1.62089822, 0.01195233, -0.46941765, -1.37134480, -0.33935883, 0.44591340, -0.01731395];
const vi = [0.32558477, 0.19458112, 0.41536797, 0.02001003, 0.05121017, 0.00690562, 0.22301725, 0.00396158, 0.05643421, 0.07302479, 0.01241221, 0.53250585, 0.07140466];

describe('fePool', () => {
  it('computes fixed-effect pooled estimate', () => {
    const r = fePool(yi, vi);
    // metafor FE theta: -0.4303
    assert.ok(Math.abs(r.theta - (-0.4303)) < 0.001, `theta=${r.theta}`);
    // metafor FE Q: 152.233
    assert.ok(r.Q > 100, `Q=${r.Q} should be large for heterogeneous data`);
  });
});

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
    // metafor: tau2 = 0.3181 (PM)
    assert.ok(Math.abs(r.tau2 - 0.3181) < 0.002, `tau2=${r.tau2}`);
  });
});

describe('rePool', () => {
  it('computes correct pooled estimate for given tau2', () => {
    const r = rePool(yi, vi, 0.3132);
    // metafor REML pooled: -0.7145
    assert.ok(Math.abs(r.theta - (-0.7145)) < 0.001, `theta=${r.theta}`);
  });
  it('returns study weights that sum to 100', () => {
    const r = rePool(yi, vi, 0.3132);
    const wPct = r.ws.map(w => 100 * w / r.sumWs);
    const total = wPct.reduce((a,b) => a+b, 0);
    assert.ok(Math.abs(total - 100) < 0.01);
  });
});
