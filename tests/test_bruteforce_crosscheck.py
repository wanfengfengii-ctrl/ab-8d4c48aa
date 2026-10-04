"""Cross-check the DP solver against an exhaustive brute-force enumerator
on random small instances (fixed seed, deterministic)."""

import itertools
import random

from nanopore_aligner.solver import align


def _compositions(n, k, lo, hi):
    """All ways to split n observations into k segments of lo..hi samples."""
    if k == 1:
        if lo <= n <= hi:
            yield (n,)
        return
    for first in range(lo, hi + 1):
        rest = n - first
        if rest < (k - 1) * lo or rest > (k - 1) * hi:
            continue
        for tail in _compositions(rest, k - 1, lo, hi):
            yield (first,) + tail


def brute_force(levels, obs, d_min, d_max, limit, dwell_min, dwell_max, max_skips):
    """Optimal objective key (skips, sum, max, drift, boundaries) by exhaustion."""
    m, n = len(levels), len(obs)
    best = None
    interior = range(1, m - 1)
    for drift in range(d_min, d_max + 1):
        for s in range(max_skips + 1):
            for skipped in itertools.combinations(interior, s):
                skipped = set(skipped)
                adopted = [j for j in range(m) if j not in skipped]
                k = len(adopted)
                for comp in _compositions(n, k, dwell_min, dwell_max):
                    starts = [0]
                    for length in comp[:-1]:
                        starts.append(starts[-1] + length)
                    total = 0.0
                    peak = 0.0
                    feasible = True
                    for seg_idx, (j, start, length) in enumerate(
                        zip(adopted, starts, comp)
                    ):
                        expected = levels[j] + drift
                        for t in range(start, start + length):
                            r = abs(obs[t] - expected)
                            if r > limit + 1e-9:
                                feasible = False
                                break
                            total += r
                            peak = max(peak, r)
                        if not feasible:
                            break
                    if not feasible:
                        continue
                    key = (s, total, peak, drift, tuple(starts))
                    if best is None or key < best:
                        best = key
    return best


def _random_case(rng):
    m = rng.randint(8, 10)
    levels = [rng.randint(-30, 30) for _ in range(m)]
    dwell_min = rng.randint(1, 2)
    dwell_max = rng.randint(dwell_min, 3)
    max_skips = rng.randint(0, 2)
    # Build observations from a random true alignment, then add noise.
    true_drift = rng.randint(-3, 3)
    adopted = sorted(
        {0, m - 1} | set(rng.sample(range(m), rng.randint(max(2, m - 3), m)))
    )
    obs = []
    for j in adopted:
        for _ in range(rng.randint(dwell_min, dwell_max)):
            obs.append(levels[j] + true_drift + rng.randint(-2, 2))
    while len(obs) < 8:
        obs.append(levels[adopted[-1]] + true_drift)
    d_lo = rng.randint(-4, 0)
    d_hi = rng.randint(0, 4)
    limit = rng.randint(0, 4)
    return levels, obs, d_lo, d_hi, limit, dwell_min, dwell_max, max_skips


def test_solver_matches_brute_force_on_random_instances():
    rng = random.Random(20261004)
    checked = 0
    for _ in range(300):
        levels, obs, d_lo, d_hi, limit, dw_lo, dw_hi, max_skips = _random_case(rng)
        expected = brute_force(levels, obs, d_lo, d_hi, limit, dw_lo, dw_hi, max_skips)
        got = align(levels, obs, d_lo, d_hi, limit, dw_lo, dw_hi, max_skips)
        if expected is None:
            assert got is None, (levels, obs, d_lo, d_hi, limit, dw_lo, dw_hi, max_skips)
            continue
        assert got is not None, (levels, obs, d_lo, d_hi, limit, dw_lo, dw_hi, max_skips)
        got_key = (
            got.skips,
            got.sum_abs_residual,
            got.max_abs_residual,
            got.drift,
            got.boundaries,
        )
        exp_key = (expected[0], expected[1], expected[2], expected[3], expected[4])
        assert len(got_key[4]) == len(exp_key[4])
        for g, e in zip(got_key, exp_key):
            assert abs(g - e) < 1e-9 if isinstance(g, float) else g == e
        # The returned alignment itself must be feasible.
        assert got.level_indices[0] == 0
        assert got.level_indices[-1] == len(levels) - 1
        assert len(got.boundaries) == len(got.level_indices)
        bounds = list(got.boundaries) + [len(obs)]
        for idx, j in enumerate(got.level_indices):
            a, b = bounds[idx], bounds[idx + 1]
            assert dw_lo <= b - a <= dw_hi
            for t in range(a, b):
                assert abs(obs[t] - (levels[j] + got.drift)) <= limit + 1e-9
        checked += 1
    assert checked > 100  # make sure we actually compared feasible cases
