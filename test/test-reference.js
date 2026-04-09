// test/test-reference.js — Validates JS engine against R metafor 4.8-0 reference values
import { metaAnalysis, fePool, estimators, rePool, hksjCI, predictionInterval } from '../meta-stats-core.js';
import { readFileSync } from 'node:fs';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

const ref = JSON.parse(readFileSync(new URL('./reference-values.json', import.meta.url)));
const bcg = ref.bcg;
const yi = bcg.yi;
const vi = bcg.vi;

describe('vs R metafor 4.8-0 (BCG vaccine data)', () => {
  it('Fixed-effect theta matches within 1e-6', () => {
    const fe = fePool(yi, vi);
    assert.ok(Math.abs(fe.theta - bcg.fe.theta[0]) < 1e-6,
      `FE theta: JS=${fe.theta} R=${bcg.fe.theta[0]}`);
  });

  it('Cochran Q matches within 1e-4', () => {
    const fe = fePool(yi, vi);
    assert.ok(Math.abs(fe.Q - bcg.fe.Q[0]) < 1e-4,
      `Q: JS=${fe.Q} R=${bcg.fe.Q[0]}`);
  });

  it('REML tau² matches within 1e-5', () => {
    const r = estimators.REML(yi, vi);
    assert.ok(Math.abs(r.tau2 - bcg.reml.tau2[0]) < 1e-5,
      `REML tau2: JS=${r.tau2} R=${bcg.reml.tau2[0]}`);
  });

  it('DL tau² matches within 1e-6', () => {
    const r = estimators.DL(yi, vi);
    assert.ok(Math.abs(r.tau2 - bcg.dl.tau2[0]) < 1e-6,
      `DL tau2: JS=${r.tau2} R=${bcg.dl.tau2[0]}`);
  });

  it('PM tau² matches within 1e-4', () => {
    const r = estimators.PM(yi, vi);
    assert.ok(Math.abs(r.tau2 - bcg.pm.tau2[0]) < 1e-4,
      `PM tau2: JS=${r.tau2} R=${bcg.pm.tau2[0]}`);
  });

  it('REML pooled theta matches within 1e-6', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    assert.ok(Math.abs(pool.theta - bcg.reml.theta[0]) < 1e-6,
      `REML theta: JS=${pool.theta} R=${bcg.reml.theta[0]}`);
  });

  it('REML SE matches within 1e-6', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    assert.ok(Math.abs(pool.se - bcg.reml.se[0]) < 1e-6,
      `REML se: JS=${pool.se} R=${bcg.reml.se[0]}`);
  });

  it('HKSJ CI matches within 1e-5', () => {
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const ci = hksjCI(pool.theta, yi, vi, r.tau2, yi.length);
    assert.ok(Math.abs(ci[0] - bcg.hksj.ci[0]) < 1e-5,
      `HKSJ CI lower: JS=${ci[0]} R=${bcg.hksj.ci[0]}`);
    assert.ok(Math.abs(ci[1] - bcg.hksj.ci[1]) < 1e-5,
      `HKSJ CI upper: JS=${ci[1]} R=${bcg.hksj.ci[1]}`);
  });

  it('PI uses t_{k-2} (deliberately wider than metafor default z-based PI)', () => {
    // NOTE: metafor default predict() uses z=1.96 for PI.
    // Our engine uses t_{k-2} per advanced-stats.md rules — MORE CONSERVATIVE.
    // So our PI is wider. This is BY DESIGN, not a bug.
    const r = estimators.REML(yi, vi);
    const pool = rePool(yi, vi, r.tau2);
    const pi = predictionInterval(pool.theta, r.tau2, pool.se, yi.length);
    // Our PI should be wider than R's z-based PI
    const jsWidth = pi[1] - pi[0];
    const rWidth = bcg.pi[1] - bcg.pi[0];
    assert.ok(jsWidth > rWidth, `JS PI width (${jsWidth.toFixed(4)}) should be wider than R z-based PI (${rWidth.toFixed(4)})`);
    // But center should match
    const jsCenter = (pi[0] + pi[1]) / 2;
    const rCenter = (bcg.pi[0] + bcg.pi[1]) / 2;
    assert.ok(Math.abs(jsCenter - rCenter) < 1e-5, 'PI center matches');
  });

  it('I² matches within 0.01%', () => {
    const ma = metaAnalysis(yi, vi);
    assert.ok(Math.abs(ma.I2 - bcg.reml.I2[0]) < 0.01,
      `I2: JS=${ma.I2} R=${bcg.reml.I2[0]}`);
  });
});
