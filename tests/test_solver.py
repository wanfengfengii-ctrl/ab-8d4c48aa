"""Unit tests for the joint alignment solver."""

from nanopore_aligner.solver import align

LEVELS = [10, 20, 30, 40, 50, 60, 70, 80]


def test_exact_alignment_zero_residual():
    obs = [lv for lv in LEVELS for _ in range(2)]
    res = align(LEVELS, obs, -2, 2, 0, 1, 3, 2)
    assert res is not None
    assert res.drift == 0
    assert res.skips == 0
    assert res.sum_abs_residual == 0
    assert res.max_abs_residual == 0
    assert res.boundaries == tuple(range(0, 16, 2))
    assert res.level_indices == tuple(range(8))


def test_drift_tiebreak_prefers_smaller_drift():
    # d=1 gives residuals (0,1) per level, d=2 gives (-1,0): identical
    # (skips, sum, max), so the smaller drift must win.
    obs = [x for lv in LEVELS for x in (lv + 1, lv + 2)]
    res = align(LEVELS, obs, -2, 2, 3, 1, 3, 2)
    assert res is not None
    assert res.drift == 1
    assert res.skips == 0
    assert res.sum_abs_residual == 8
    assert res.max_abs_residual == 1


def test_skip_count_dominates_residual_sum():
    # Skipping level index 1 would give sum 0, but skips=0 is feasible with
    # sum 10 — the skip count is minimised first, so skips=0 must win.
    levels = [0, 10, 20, 30, 40, 50, 60, 70]
    obs = [0, 20, 20, 30, 40, 50, 60, 70]
    res = align(levels, obs, 0, 0, 10, 1, 3, 2)
    assert res is not None
    assert res.skips == 0
    assert res.sum_abs_residual == 10
    assert res.boundaries == tuple(range(8))


def test_max_residual_tiebreak_before_drift():
    # d in {0,1,2} all give sum 8; d=1 alone gives max 1 (others give 2),
    # so the max-residual criterion must beat the smaller drifts.
    obs = [0, 0, 0, 0, 2, 2, 2, 2]
    res = align([0] * 8, obs, 0, 2, 5, 1, 3, 2)
    assert res is not None
    assert res.drift == 1
    assert res.sum_abs_residual == 8
    assert res.max_abs_residual == 1


def test_boundary_tiebreak_prefers_lexicographically_small():
    # 9 observations over 8 identical levels: exactly one segment has
    # length 2, everything ties, so the lexicographically smallest start
    # sequence (0,1,...,7) must be returned.
    res = align([0] * 8, [0] * 9, 0, 0, 0, 1, 3, 2)
    assert res is not None
    assert res.boundaries == (0, 1, 2, 3, 4, 5, 6, 7)
    assert res.level_indices == tuple(range(8))


def test_internal_skip_recovers_missing_level():
    # No observation is compatible with level 40 (index 3).
    obs = [11, 12, 21, 22, 31, 32, 51, 52, 61, 62, 71, 72, 81, 82]
    res = align(LEVELS, obs, -2, 2, 3, 1, 3, 2)
    assert res is not None
    assert res.skips == 1
    assert res.level_indices == (0, 1, 2, 4, 5, 6, 7)
    assert res.drift == 1
    assert res.boundaries == (0, 2, 4, 6, 8, 10, 12)
    assert res.sum_abs_residual == 7
    assert res.max_abs_residual == 1


def test_skip_budget_is_enforced():
    # Same trace as above but no skip allowed -> infeasible.
    obs = [11, 12, 21, 22, 31, 32, 51, 52, 61, 62, 71, 72, 81, 82]
    assert align(LEVELS, obs, -2, 2, 3, 1, 3, 0) is None


def test_dwell_range_is_respected():
    res = align([0] * 8, [0] * 16, 0, 0, 0, 2, 2, 2)
    assert res is not None
    bounds = res.boundaries + (16,)
    assert all(b - a == 2 for a, b in zip(bounds, bounds[1:]))


def test_structural_infeasibility_too_many_observations():
    # 8 levels with at most 3 samples each cover at most 24 observations.
    assert align(list(range(8)), [0] * 30, 0, 0, 100, 1, 3, 2) is None


def test_structural_infeasibility_too_few_observations():
    # 24 levels with at most 2 skips need at least 22 segments/observations.
    assert align(list(range(24)), [0] * 8, 0, 0, 100, 1, 3, 2) is None


def test_residual_infeasibility_when_drift_out_of_range():
    obs = [x for lv in LEVELS for x in (lv + 1, lv + 2)]
    assert align(LEVELS, obs, 5, 6, 3, 1, 3, 2) is None


def test_float_observations():
    res = align([0] * 8, [0.5] * 8, 0, 0, 1, 1, 1, 2)
    assert res is not None
    assert abs(res.sum_abs_residual - 4.0) < 1e-9
    assert abs(res.max_abs_residual - 0.5) < 1e-9


def test_first_and_last_levels_are_always_used():
    # Even though interior levels may be skipped, the anchors may not.
    levels = [100, 0, 0, 0, 0, 0, 0, 100]
    obs = [100, 0, 0, 0, 0, 0, 0, 100]
    res = align(levels, obs, 0, 0, 0, 1, 1, 2)
    assert res is not None
    assert res.level_indices[0] == 0
    assert res.level_indices[-1] == len(levels) - 1
