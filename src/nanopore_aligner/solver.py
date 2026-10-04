"""Joint alignment of ion-current observations to reference levels.

Given 8..24 integer reference levels and 8..60 observations, the solver
chooses — jointly, without any pre-estimated baseline —

* an integer drift ``d`` inside the submitted closed interval,
* a subsequence of the reference levels that keeps the first and the last
  level and skips at most ``max_internal_skips`` interior levels, and
* a contiguous partition of the observations into one non-empty segment
  per adopted level (segment length inside the dwell range),

such that every observation belongs to exactly one adopted level and every
residual ``obs - (level + drift)`` satisfies ``abs(residual) <= limit``.

Among all feasible alignments the solver minimises, in lexicographic
order:

1. the number of skipped (non-adopted) reference levels,
2. the sum of absolute residuals,
3. the maximum absolute residual,
4. the drift value,
5. the boundary sequence (segment start indices, compared
   lexicographically).

The search is an exact dynamic program per integer drift.  State
``dp[j][i][s]`` holds the best ``(sum_abs, max_abs, boundaries, path)``
for covering ``obs[:i]`` with adopted levels ending exactly at reference
``j`` (level ``j`` owns the last segment, which ends at observation
``i``), using exactly ``s`` interior skips among references ``0..j``.
Because at most ``max_skips`` levels may be skipped, consecutive adopted
references differ by at most ``max_skips + 1``, so each state has at most
``(max_skips + 1) * dwell_choices`` predecessors.  Keeping the
lexicographically smallest boundary tuple per state is valid: all
candidates reaching one state have equally many segments, and appending a
common suffix preserves tuple order.
"""

from __future__ import annotations

from dataclasses import dataclass

_EPS = 1e-9


@dataclass(frozen=True)
class Alignment:
    """One optimal alignment for a trace."""

    drift: int
    skips: int
    sum_abs_residual: float
    max_abs_residual: float
    boundaries: tuple[int, ...]  # start index of each segment (len == #segments)
    level_indices: tuple[int, ...]  # adopted reference index per segment


def _segment_stats(obs, start, end, expected, limit):
    """Return ``(sum_abs, max_abs)`` for ``obs[start:end]`` against one shifted
    level, or ``None`` if any absolute residual exceeds ``limit``."""
    total = 0.0
    peak = 0.0
    for t in range(start, end):
        r = obs[t] - expected
        if r < 0:
            r = -r
        if r > peak:
            peak = r
            if peak > limit + _EPS:
                return None
        total += r
    return total, peak


def _best_for_drift(levels, obs, drift, limit, dwell_min, dwell_max, max_skips):
    """Optimal ``(sum, max, boundaries, path)`` for one fixed integer drift,
    minimised as ``(skips, sum, max, boundaries)``, or ``None`` if infeasible."""
    m = len(levels)
    n = len(obs)
    shifted = [lv + drift for lv in levels]

    dp = [[[None] * (max_skips + 1) for _ in range(n + 1)] for _ in range(m)]

    first = shifted[0]
    for length in range(dwell_min, dwell_max + 1):
        if length > n:
            break
        stats = _segment_stats(obs, 0, length, first, limit)
        if stats is not None:
            dp[0][length][0] = (stats[0], stats[1], (0,), (0,))

    for j in range(1, m):
        expected = shifted[j]
        dpj = dp[j]
        jp_lo = max(0, j - 1 - max_skips)
        for i in range(dwell_min, n + 1):
            for length in range(dwell_min, dwell_max + 1):
                a = i - length
                if a < 0:
                    break
                # Collect feasible predecessors before costing the segment.
                preds = []
                for jp in range(jp_lo, j):
                    skip_add = j - jp - 1
                    row = dp[jp][a]
                    sp_hi = max_skips - skip_add
                    for sp in range(sp_hi + 1):
                        prev = row[sp]
                        if prev is not None:
                            preds.append((prev, sp + skip_add))
                if not preds:
                    continue
                stats = _segment_stats(obs, a, i, expected, limit)
                if stats is None:
                    continue
                ssum, smax = stats
                for prev, s in preds:
                    psum, pmax, pbounds, ppath = prev
                    cand_sum = psum + ssum
                    cand_max = pmax if pmax >= smax else smax
                    cand_bounds = pbounds + (a,)
                    cur = dpj[i][s]
                    if cur is None or (cand_sum, cand_max, cand_bounds) < cur[:3]:
                        dpj[i][s] = (cand_sum, cand_max, cand_bounds, ppath + (j,))

    best = None
    best_key = None
    last_row = dp[m - 1][n]
    for s in range(max_skips + 1):
        v = last_row[s]
        if v is None:
            continue
        key = (s, v[0], v[1], v[2])
        if best_key is None or key < best_key:
            best_key = key
            best = v
    return best


def align(
    levels,
    observations,
    drift_min,
    drift_max,
    residual_limit,
    dwell_min,
    dwell_max,
    max_internal_skips,
):
    """Return the lexicographically optimal :class:`Alignment`, or ``None``
    when no feasible alignment exists."""
    if dwell_min < 1 or dwell_max < dwell_min:
        raise ValueError("invalid dwell range")
    if max_internal_skips < 0:
        raise ValueError("max_internal_skips must be >= 0")
    if drift_min > drift_max:
        raise ValueError("drift_min must be <= drift_max")

    m = len(levels)
    n = len(observations)

    # Structural feasibility: some adopted-level count k must be able to
    # partition n observations into segments of dwell_min..dwell_max samples.
    k_lo = m - max_internal_skips
    if k_lo < 1:
        k_lo = 1
    if not any(k * dwell_min <= n <= k * dwell_max for k in range(k_lo, m + 1)):
        return None

    best = None
    best_key = None
    for drift in range(drift_min, drift_max + 1):
        solved = _best_for_drift(
            levels, observations, drift, residual_limit,
            dwell_min, dwell_max, max_internal_skips,
        )
        if solved is None:
            continue
        ssum, smax, bounds, path = solved
        skips = m - len(path)
        key = (skips, ssum, smax, drift, bounds)
        if best_key is None or key < best_key:
            best_key = key
            best = Alignment(
                drift=drift,
                skips=skips,
                sum_abs_residual=ssum,
                max_abs_residual=smax,
                boundaries=bounds,
                level_indices=path,
            )
    return best
