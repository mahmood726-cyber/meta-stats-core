import { radialEgger, petersTest, petPeese, trimFill } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg, escalc measure="RR")
// Verified against R metafor 4.8-0 (8 significant figures to avoid rounding artefacts)
const yi = [-0.88931133, -1.58538866, -1.34807315, -1.44155119, -0.21754732, -0.78611559, -1.62089822, 0.01195233, -0.46941765, -1.37134480, -0.33935883, 0.44591340, -0.01731395];
const vi = [0.32558477, 0.19458112, 0.41536797, 0.02001003, 0.05121017, 0.00690562, 0.22301725, 0.00396158, 0.05643421, 0.07302479, 0.01241221, 0.53250585, 0.07140466];

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
    const ni = [100, 200, 150, 300, 80, 250, 120, 400, 180, 60, 90, 350, 110];
    const r = petersTest(yi, vi, ni);
    assert.ok('intercept' in r);
    assert.ok('pValue' in r);
    assert.ok(r.pValue >= 0 && r.pValue <= 1);
  });
});

describe('petPeese', () => {
  it('applies conditional procedure: PET first, then PEESE', () => {
    const r = petPeese(yi, vi);
    assert.ok('method' in r); // 'PET' or 'PEESE'
    assert.ok('adjustedEffect' in r);
    assert.ok('pValue' in r);
    if (r.petPvalue < 0.10) {
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
    assert.ok('k0' in r);
    assert.ok('adjustedEffect' in r);
  });
  it('k0 is non-negative integer', () => {
    const r = trimFill(yi, vi);
    assert.ok(r.k0 >= 0);
    assert.ok(Number.isInteger(r.k0));
  });
});
