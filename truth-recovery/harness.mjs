// ============================================================
// harness.mjs — Truth-recovery yardstick for meta-stats-core.
//
// meta-stats-core is validated against R metafor to 6 decimals for its POINT
// outputs. That proves it computes each formula correctly; it does NOT tell you
// which of its methods actually RECOVERS THE TRUTH when a meta-analysis is
// distorted by heterogeneity AND publication selection together. This harness
// supplies that missing evidence: inject a known (mu, tau^2) + selection
// mechanism, then measure how often each method's interval covers the TRUE mu.
//
// Truth-first: every number is produced here from seeded simulation. Run:
//   node truth-recovery/harness.mjs --reps 400
// ============================================================

import {
  estimators, rePool, fePool, hksjCI, waldCI, petPeese, heterogeneity,
} from '../meta-stats-core.js';
import { generate, makeRng, SCENARIOS } from './dgp.mjs';

const BASE_SEED = 20260613;

// --- the methods under test: each (yi, vi) -> { mu, lo, hi, hasCI } ---
const METHODS = {
  'DL+Wald': (yi, vi) => {
    const tau2 = estimators.DL(yi, vi).tau2;
    const p = rePool(yi, vi, tau2);
    const [lo, hi] = waldCI(p.theta, p.se);
    return { mu: p.theta, lo, hi, hasCI: true };
  },
  'REML+Wald': (yi, vi) => {
    const tau2 = estimators.REML(yi, vi).tau2;
    const p = rePool(yi, vi, tau2);
    const [lo, hi] = waldCI(p.theta, p.se);
    return { mu: p.theta, lo, hi, hasCI: true };
  },
  'REML+HKSJ': (yi, vi) => {
    const tau2 = estimators.REML(yi, vi).tau2;
    const p = rePool(yi, vi, tau2);
    const [lo, hi] = hksjCI(p.theta, yi, vi, tau2, yi.length);
    return { mu: p.theta, lo, hi, hasCI: true };
  },
  'PM+HKSJ': (yi, vi) => {
    const tau2 = estimators.PM(yi, vi).tau2;
    const p = rePool(yi, vi, tau2);
    const [lo, hi] = hksjCI(p.theta, yi, vi, tau2, yi.length);
    return { mu: p.theta, lo, hi, hasCI: true };
  },
  // PET-PEESE: small-study-effect bias correction (point estimate only; no
  // honest CI is produced by the library, so it is scored on bias/RMSE only).
  'PET-PEESE': (yi, vi) => {
    const r = petPeese(yi, vi);
    return { mu: r.adjustedEffect, lo: NaN, hi: NaN, hasCI: false };
  },
};

function mean(a) { return a.reduce((x, y) => x + y, 0) / a.length; }

export function runCell(mu, tau2, k, scenario, reps, rng) {
  const acc = {};
  for (const name of Object.keys(METHODS)) {
    acc[name] = { covered: 0, nCI: 0, biasSum: 0, sqSum: 0, widthSum: 0,
                  reject0: 0, n: 0 };
  }
  for (let r = 0; r < reps; r++) {
    const { yi, vi } = generate(mu, tau2, k, scenario, rng);
    for (const [name, fn] of Object.entries(METHODS)) {
      let out;
      try { out = fn(yi, vi); } catch { continue; }
      if (!isFinite(out.mu)) continue;
      const a = acc[name];
      a.n++;
      a.biasSum += out.mu - mu;
      a.sqSum += (out.mu - mu) ** 2;
      if (out.hasCI && isFinite(out.lo) && isFinite(out.hi)) {
        a.nCI++;
        if (out.lo <= mu && mu <= out.hi) a.covered++;
        a.widthSum += out.hi - out.lo;
        if (!(out.lo <= 0 && 0 <= out.hi)) a.reject0++;   // 0 excluded
      }
    }
  }
  const res = {};
  for (const [name, a] of Object.entries(acc)) {
    res[name] = {
      n: a.n,
      bias: a.n ? +(a.biasSum / a.n).toFixed(4) : null,
      rmse: a.n ? +Math.sqrt(a.sqSum / a.n).toFixed(4) : null,
      coverage: a.nCI ? +(a.covered / a.nCI).toFixed(4) : null,
      meanWidth: a.nCI ? +(a.widthSum / a.nCI).toFixed(4) : null,
      reject0: a.nCI ? +(a.reject0 / a.nCI).toFixed(4) : null,
    };
  }
  return res;
}

export function runGrid({ reps = 400, ks = [5, 10, 15, 25], tau2 = 0.05,
                          scenarios = SCENARIOS, mu = 0.3 } = {}) {
  const rng = makeRng(BASE_SEED);
  const grid = [];
  for (const scen of scenarios) {
    for (const k of ks) {
      grid.push({ scen, k, results: runCell(mu, tau2, k, scen, reps, rng) });
    }
  }
  return grid;
}

// Aggregate mean coverage / |bias| per method across all selection cells.
export function summarize(grid) {
  const names = Object.keys(grid[0].results);
  const out = {};
  for (const name of names) {
    const cov = [], absBias = [], rmse = [];
    for (const cell of grid) {
      const r = cell.results[name];
      if (r.coverage != null) cov.push(r.coverage);
      if (r.bias != null) absBias.push(Math.abs(r.bias));
      if (r.rmse != null) rmse.push(r.rmse);
    }
    out[name] = {
      meanCoverage: cov.length ? +mean(cov).toFixed(4) : null,
      meanAbsBias: +mean(absBias).toFixed(4),
      meanRmse: +mean(rmse).toFixed(4),
    };
  }
  return out;
}

// --- CLI ---
const isMain = import.meta.url === `file://${process.argv[1]}` ||
               process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const arg = (f, d) => {
    const i = process.argv.indexOf(f);
    return i >= 0 ? Number(process.argv[i + 1]) : d;
  };
  const reps = arg('--reps', 400);
  const t0 = Date.now();
  const grid = runGrid({ reps });
  const summary = summarize(grid);
  console.log(`\n# Truth-recovery yardstick — meta-stats-core`);
  console.log(`reps=${reps}/cell  mu=0.3  tau2=0.05  seed=${BASE_SEED}\n`);
  console.log('## Mean over all selection cells (coverage of TRUE mu)\n');
  console.log('method        meanCov  meanAbsBias  meanRMSE');
  for (const [name, s] of Object.entries(summary)) {
    console.log(name.padEnd(13),
      String(s.meanCoverage ?? '   n/a').padStart(7),
      String(s.meanAbsBias).padStart(11),
      String(s.meanRmse).padStart(9));
  }
  console.log('\n## Per-cell coverage (REML+HKSJ vs DL+Wald)\n');
  console.log('scenario       k   DL+Wald  REML+HKSJ  PM+HKSJ');
  for (const cell of grid) {
    const d = cell.results['DL+Wald'].coverage;
    const h = cell.results['REML+HKSJ'].coverage;
    const p = cell.results['PM+HKSJ'].coverage;
    console.log(cell.scen.padEnd(13), String(cell.k).padStart(3),
      String(d).padStart(8), String(h).padStart(9), String(p).padStart(8));
  }
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
