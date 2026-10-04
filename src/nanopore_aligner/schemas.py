"""Request/response schemas for the current-trace alignment API."""

from __future__ import annotations

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

# Accepted field ranges. Requests violating any of them are rejected
# (HTTP 422) before any solving happens.
LEVELS_MIN, LEVELS_MAX = 8, 24
OBS_MIN, OBS_MAX = 8, 60
LEVEL_ABS_MAX = 1_000_000
OBS_ABS_MAX = 1_000_000_000
DRIFT_ABS_MAX = 1_000_000
DRIFT_WIDTH_MAX = 1_000
RESIDUAL_LIMIT_MAX = 1_000_000_000
DWELL_MIN, DWELL_MAX = 1, 3
SKIPS_MAX = 2


class DriftInterval(BaseModel):
    """Closed integer interval ``[min, max]`` for the uniform drift."""

    model_config = ConfigDict(extra="forbid")

    min: Annotated[int, Field(strict=True, ge=-DRIFT_ABS_MAX, le=DRIFT_ABS_MAX)]
    max: Annotated[int, Field(strict=True, ge=-DRIFT_ABS_MAX, le=DRIFT_ABS_MAX)]

    @model_validator(mode="after")
    def _check_interval(self) -> "DriftInterval":
        if self.min > self.max:
            raise ValueError("drift.min must be <= drift.max")
        if self.max - self.min > DRIFT_WIDTH_MAX:
            raise ValueError(
                f"drift interval width must be <= {DRIFT_WIDTH_MAX}"
            )
        return self


class DwellRange(BaseModel):
    """Allowed number of samples per adopted level (1..3)."""

    model_config = ConfigDict(extra="forbid")

    min: Annotated[int, Field(strict=True, ge=DWELL_MIN, le=DWELL_MAX)] = DWELL_MIN
    max: Annotated[int, Field(strict=True, ge=DWELL_MIN, le=DWELL_MAX)] = DWELL_MAX

    @model_validator(mode="after")
    def _check_range(self) -> "DwellRange":
        if self.min > self.max:
            raise ValueError("dwell.min must be <= dwell.max")
        return self


class AlignRequest(BaseModel):
    """Payload for ``POST /api/current-traces/align``."""

    model_config = ConfigDict(extra="forbid")

    levels: Annotated[
        list[Annotated[int, Field(strict=True, ge=-LEVEL_ABS_MAX, le=LEVEL_ABS_MAX)]],
        Field(min_length=LEVELS_MIN, max_length=LEVELS_MAX),
    ]
    observations: Annotated[
        list[
            Annotated[
                float,
                Field(allow_inf_nan=False, ge=-OBS_ABS_MAX, le=OBS_ABS_MAX),
            ]
        ],
        Field(min_length=OBS_MIN, max_length=OBS_MAX),
    ]
    drift: DriftInterval
    residual_limit: Annotated[
        float, Field(allow_inf_nan=False, ge=0, le=RESIDUAL_LIMIT_MAX)
    ]
    dwell: DwellRange = Field(default_factory=DwellRange)
    max_internal_skips: Annotated[int, Field(strict=True, ge=0, le=SKIPS_MAX)] = (
        SKIPS_MAX
    )


class Segment(BaseModel):
    """Sample interval owned by one adopted reference level."""

    level_index: int
    level: int
    expected: float  # level + drift
    start: int  # inclusive observation index
    end: int  # exclusive observation index
    length: int


class ResidualEvidence(BaseModel):
    """Per-observation residual proof."""

    index: int
    level_index: int
    observed: float
    expected: float
    residual: float  # observed - expected
    abs_residual: float


class Objective(BaseModel):
    """The lexicographic objective vector of the returned alignment."""

    skips: int
    sum_abs_residual: float
    max_abs_residual: float
    drift: int
    boundaries: list[int]


class AlignedResponse(BaseModel):
    status: Literal["aligned"] = "aligned"
    drift: int
    skips: int
    skipped_level_indices: list[int]
    boundaries: list[int]
    segments: list[Segment]
    residuals: list[ResidualEvidence]
    objective: Objective


class InfeasibleResponse(BaseModel):
    status: Literal["infeasible"] = "infeasible"
    reason: str
    detail: dict[str, Any] = Field(default_factory=dict)


AlignResponse = Annotated[
    AlignedResponse | InfeasibleResponse, Field(discriminator="status")
]
