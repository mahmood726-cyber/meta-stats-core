import { heterogeneity, fePool } from '../meta-stats-core.js';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

// BCG vaccine dataset (metafor::dat.bcg, escalc measure="RR")
// Verified against R metafor 4.8-0 (8 significant figures to avoid rounding artefacts)
const yi = [-0.88931133, -1.58538866, -1.34807315, -1.44155119, -0.21754732, -0.78611559, -1.62089822, 0.01195233, -0.46941765, -1.37134480, -0.33935883, 0.44591340, -0.01731395];
const vi = [0.32558477, 0.19458112, 0.41536797, 0.02001003, 0.05121017, 0.00690562, 0.22301725, 0.00396158, 0.05643421, 0.07302479, 0.01241221, 0.53250585, 0.07140466];

describe('heterogeneity', () => {
  it('computes I² correctly', () => {
    const fe = fePool(yi, vi);
    const h = heterogeneity(fe.Q, yi.length, 0.3132);
    // I2 should be high for BCG data (>80%)
    assert.ok(h.I2 > 80, `I2=${h.I2}`);
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
    const h = heterogeneity(1.5, 3, 0.0);
    assert.ok(h.I2 === 0);
    assert.ok(h.lowPowerWarning === true, 'Should warn about low power');
  });
  it('Q p-value is computed', () => {
    const fe = fePool(yi, vi);
    const h = heterogeneity(fe.Q, yi.length, 0.3132);
    assert.ok(h.Qp < 0.001, 'Q p-value should be tiny for BCG data');
  });
});
