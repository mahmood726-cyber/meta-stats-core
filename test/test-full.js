import { metaAnalysis } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg, escalc measure="RR")
// Verified against R metafor 4.8-0 (8 significant figures to avoid rounding artefacts)
const yi = [-0.88931133, -1.58538866, -1.34807315, -1.44155119, -0.21754732, -0.78611559, -1.62089822, 0.01195233, -0.46941765, -1.37134480, -0.33935883, 0.44591340, -0.01731395];
const vi = [0.32558477, 0.19458112, 0.41536797, 0.02001003, 0.05121017, 0.00690562, 0.22301725, 0.00396158, 0.05643421, 0.07302479, 0.01241221, 0.53250585, 0.07140466];

describe('metaAnalysis (full pipeline)', () => {
  it('returns all required components', () => {
    const r = metaAnalysis(yi, vi);
    assert.ok('theta' in r);
    assert.ok('ci' in r && r.ci.length === 2);
    assert.ok('pi' in r && r.pi.length === 2);
    assert.ok('tau2' in r);
    assert.ok('I2' in r);
    assert.ok('I2_ci' in r);
    assert.ok('Q' in r);
    assert.ok('weights' in r);
    assert.strictEqual(r.estimator, 'REML');
    assert.strictEqual(r.ciMethod, 'HKSJ');
  });
  it('defaults to REML + HKSJ', () => {
    const r = metaAnalysis(yi, vi);
    assert.strictEqual(r.estimator, 'REML');
    assert.strictEqual(r.ciMethod, 'HKSJ');
  });
  it('allows DL method override', () => {
    const r = metaAnalysis(yi, vi, { method: 'DL' });
    assert.strictEqual(r.estimator, 'DL');
  });
  it('BCG pooled matches metafor within tolerance', () => {
    const r = metaAnalysis(yi, vi);
    // REML theta ~ -0.7145, tau2 ~ 0.3132
    assert.ok(Math.abs(r.tau2 - 0.3132) < 0.002, `tau2=${r.tau2}`);
    assert.ok(r.ci[0] < r.theta && r.theta < r.ci[1], 'theta within CI');
    assert.ok(r.pi[0] < r.ci[0], 'PI wider than CI');
  });
  it('weights sum to approximately 100', () => {
    const r = metaAnalysis(yi, vi);
    const sum = r.weights.reduce((a,b) => a+b, 0);
    assert.ok(Math.abs(sum - 100) < 0.01);
  });
  it('k=2 produces valid result with no PI', () => {
    const r = metaAnalysis(yi.slice(0,2), vi.slice(0,2));
    assert.ok(isNaN(r.pi[0]), 'PI undefined for k<3');
    assert.ok(!isNaN(r.ci[0]), 'CI should still work');
  });
  it('k=1 single study does not poison theta/se/CI with NaN', () => {
    // df=0 → DL C=0 (0/0); guard sets tau2=0 so the single estimate passes through.
    const r = metaAnalysis([-0.5], [0.1]);
    assert.strictEqual(r.theta, -0.5, 'theta is the single study effect');
    assert.ok(Math.abs(r.se - Math.sqrt(0.1)) < 1e-12, 'se = sqrt(v)');
    assert.strictEqual(r.tau2, 0, 'tau2 = 0 for k=1');
    assert.ok(!isNaN(r.ci[0]) && !isNaN(r.ci[1]), 'CI is finite');
    assert.ok(r.ci[0] < r.theta && r.theta < r.ci[1], 'theta within CI');
    assert.ok(isNaN(r.pi[0]), 'PI undefined for k<3');
    assert.strictEqual(r.I2, 0, 'I2=0 (heterogeneity undefined for single study)');
  });
});
