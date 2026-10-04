"""FastAPI application exposing the current-trace alignment endpoint."""

from __future__ import annotations

from fastapi import FastAPI

from .schemas import (
    AlignRequest,
    AlignResponse,
    AlignedResponse,
    InfeasibleResponse,
    Objective,
    ResidualEvidence,
    Segment,
)
from .solver import align

app = FastAPI(
    title="Nanopore Current-Trace Aligner",
    version="1.0.0",
    description=(
        "Re-aligns ion-current observations to known reference levels by "
        "jointly choosing an integer drift, a first/last-anchored reference "
        "subsequence and contiguous sample boundaries."
    ),
)


@app.get("/health", tags=["meta"])
def health() -> dict:
    return {"status": "ok"}


@app.get("/", tags=["meta"])
def root() -> dict:
    return {
        "service": "nanopore-aligner",
        "version": "1.0.0",
        "align": "/api/current-traces/align",
        "docs": "/docs",
    }


@app.post(
    "/api/current-traces/align",
    response_model=AlignResponse,
    tags=["align"],
    summary="Jointly align a current trace to reference levels",
)
def align_trace(req: AlignRequest) -> AlignResponse:
    n = len(req.observations)
    m = len(req.levels)

    result = align(
        levels=req.levels,
        observations=req.observations,
        drift_min=req.drift.min,
        drift_max=req.drift.max,
        residual_limit=req.residual_limit,
        dwell_min=req.dwell.min,
        dwell_max=req.dwell.max,
        max_internal_skips=req.max_internal_skips,
    )

    if result is None:
        return InfeasibleResponse(
            reason=_infeasible_reason(req),
            detail={
                "levels": m,
                "observations": n,
                "drift": {"min": req.drift.min, "max": req.drift.max},
                "residual_limit": req.residual_limit,
                "dwell": {"min": req.dwell.min, "max": req.dwell.max},
                "max_internal_skips": req.max_internal_skips,
            },
        )

    boundaries = list(result.boundaries)
    segments: list[Segment] = []
    residuals: list[ResidualEvidence] = []
    for idx, (start, ref_idx) in enumerate(zip(boundaries, result.level_indices)):
        end = boundaries[idx + 1] if idx + 1 < len(boundaries) else n
        level = req.levels[ref_idx]
        expected = level + result.drift
        segments.append(
            Segment(
                level_index=ref_idx,
                level=level,
                expected=expected,
                start=start,
                end=end,
                length=end - start,
            )
        )
        for t in range(start, end):
            r = req.observations[t] - expected
            residuals.append(
                ResidualEvidence(
                    index=t,
                    level_index=ref_idx,
                    observed=req.observations[t],
                    expected=expected,
                    residual=r,
                    abs_residual=abs(r),
                )
            )

    adopted = set(result.level_indices)
    return AlignedResponse(
        drift=result.drift,
        skips=result.skips,
        skipped_level_indices=[j for j in range(m) if j not in adopted],
        boundaries=boundaries,
        segments=segments,
        residuals=residuals,
        objective=Objective(
            skips=result.skips,
            sum_abs_residual=result.sum_abs_residual,
            max_abs_residual=result.max_abs_residual,
            drift=result.drift,
            boundaries=boundaries,
        ),
    )


def _infeasible_reason(req: AlignRequest) -> str:
    m = len(req.levels)
    n = len(req.observations)
    k_lo = max(1, m - req.max_internal_skips)
    ks = [
        k
        for k in range(k_lo, m + 1)
        if k * req.dwell.min <= n <= k * req.dwell.max
    ]
    if not ks:
        return (
            f"structurally infeasible: {n} observations cannot be partitioned "
            f"into segments of {req.dwell.min}..{req.dwell.max} samples for "
            f"any adopted-level count k in [{k_lo}, {m}] (first and last "
            f"reference levels must be used, at most "
            f"{req.max_internal_skips} interior skips)"
        )
    return (
        f"no integer drift in [{req.drift.min}, {req.drift.max}] admits an "
        f"alignment with all absolute residuals <= {req.residual_limit}"
    )
