// node --test truth-recovery/test-truth-recovery.mjs
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { generate, makeRng, SCENARIOS } from './dgp.mjs';
import { runCell, runGrid, summarize } from './harness.mjs';

describe('DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(0.3, 0.05, 12, 'step_strong', makeRng(7));
    const b = generate(0.3, 0.05, 12, 'step_strong', makeRng(7));
    assert.deepEqual(a.yi, b.yi);
    assert.deepEqual(a.vi, b.vi);
  });

  it('returns exactly k positive-variance studies for every scenario', () => {
    const rng = makeRng(1);
    for (const scen of SCENARIOS) {
      const { yi, vi } = generate(0.3, 0.05, 10, scen, rng);
      assert.equal(yi.length, 10);
      assert.equal(vi.length, 10);
      assert.ok(vi.every(v => v > 0));
    }
  });

  it('strong positive selection inflates the naive mean vs no selection', () => {
    const rng = makeRng(123);
    let none = 0, sel = 0;
    const N = 300;
    for (let i = 0; i < N; i++) {
      none += avg(generate(0.3, 0.05, 10, 'none', rng).yi);
      sel += avg(generate(0.3, 0.05, 10, 'step_strong', rng).yi);
    }
    assert.ok(sel / N > none / N, `sel=${sel / N} none=${none / N}`);
  });
});

describe('Truth-recovery yardstick (measured invariants)', () => {
  it('REML+HKSJ recovers truth at least as well as legacy DL+Wald (mean coverage)', () => {
    const grid = runGrid({ reps: 150 });
    const s = summarize(grid);
    // The library DEFAULTS to REML+HKSJ; this is the measured justification.
    assert.ok(s['REML+HKSJ'].meanCoverage >= s['DL+Wald'].meanCoverage,
      `REML+HKSJ ${s['REML+HKSJ'].meanCoverage} < DL+Wald ${s['DL+Wald'].meanCoverage}`);
  });

  it('clean data (no selection) is near nominal for REML+HKSJ; strong selection collapses coverage', () => {
    const rng = makeRng(20260613);
    const clean = runCell(0.3, 0.05, 15, 'none', 300, rng);
    const sel = runCell(0.3, 0.05, 25, 'step_strong', 300, rng);
    assert.ok(clean['REML+HKSJ'].coverage > 0.88,
      `clean coverage ${clean['REML+HKSJ'].coverage}`);
    // honest: no IV interval models selection, so strong selection wrecks it
    assert.ok(sel['REML+HKSJ'].coverage < 0.5,
      `selection coverage ${sel['REML+HKSJ'].coverage}`);
  });

  it('PET-PEESE reduces selection bias relative to RE pooling', () => {
    const rng = makeRng(42);
    const cell = runCell(0.3, 0.05, 25, 'copas_strong', 300, rng);
    assert.ok(Math.abs(cell['PET-PEESE'].bias) < Math.abs(cell['REML+HKSJ'].bias),
      `PET-PEESE |bias| ${cell['PET-PEESE'].bias} vs REML+HKSJ ${cell['REML+HKSJ'].bias}`);
  });
});

function avg(a) { return a.reduce((x, y) => x + y, 0) / a.length; }
