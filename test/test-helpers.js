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
    assert.ok(Math.abs(tQuantile(0.025, 5)) > 2.0);
  });
});

describe('chi2CDF', () => {
  it('matches R pchisq()', () => {
    // pchisq(7.81, 3) = 0.9500
    assert.ok(Math.abs(chi2CDF(7.81, 3) - 0.9500) < 0.005);
  });
});

describe('chi2Quantile', () => {
  it('matches R qchisq()', () => {
    // qchisq(0.975, 12) = 23.33666
    assert.ok(Math.abs(chi2Quantile(0.975, 12) - 23.33666) < 0.01);
    // qchisq(0.025, 12) = 4.40379
    assert.ok(Math.abs(chi2Quantile(0.025, 12) - 4.40379) < 0.01);
  });
});

describe('fisherZ', () => {
  it('transforms correlation to z', () => {
    assert.ok(Math.abs(fisherZ(0.5) - 0.5493) < 1e-3);
  });
  it('clamps at boundaries', () => {
    assert.ok(isFinite(fisherZ(1.0)));
    assert.ok(isFinite(fisherZ(-1.0)));
  });
  it('inverse roundtrips', () => {
    assert.ok(Math.abs(fisherZinv(fisherZ(0.7)) - 0.7) < 1e-6);
  });
});
