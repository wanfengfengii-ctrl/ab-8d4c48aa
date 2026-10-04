"""API-level tests: happy path, infeasible conclusion, validation rejects."""

import copy

import pytest
from fastapi.testclient import TestClient

from nanopore_aligner.main import app

client = TestClient(app)

BASE_PAYLOAD = {
    "levels": [10, 20, 30, 40, 50, 60, 70, 80],
    "observations": [11, 12, 21, 22, 31, 32, 41, 42, 51, 52, 61, 62, 71, 72, 81, 82],
    "drift": {"min": -2, "max": 2},
    "residual_limit": 3,
    "dwell": {"min": 1, "max": 3},
    "max_internal_skips": 2,
}


def _payload(**overrides):
    p = copy.deepcopy(BASE_PAYLOAD)
    p.update(overrides)
    return p


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_align_feasible_trace():
    r = client.post("/api/current-traces/align", json=BASE_PAYLOAD)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "aligned"
    assert body["drift"] == 1
    assert body["skips"] == 0
    assert body["skipped_level_indices"] == []
    assert body["boundaries"] == [0, 2, 4, 6, 8, 10, 12, 14]
    assert body["objective"] == {
        "skips": 0,
        "sum_abs_residual": 8,
        "max_abs_residual": 1,
        "drift": 1,
        "boundaries": [0, 2, 4, 6, 8, 10, 12, 14],
    }

    segments = body["segments"]
    assert len(segments) == 8
    assert segments[0]["start"] == 0
    assert segments[-1]["end"] == 16
    for a, b in zip(segments, segments[1:]):
        assert a["end"] == b["start"]
    for seg in segments:
        assert 1 <= seg["length"] <= 3
        assert seg["expected"] == seg["level"] + body["drift"]

    residuals = body["residuals"]
    assert len(residuals) == 16
    for ev in residuals:
        assert ev["expected"] == BASE_PAYLOAD["levels"][ev["level_index"]] + body["drift"]
        assert abs(ev["observed"] - ev["expected"] - ev["residual"]) < 1e-9
        assert ev["abs_residual"] == abs(ev["residual"])
        assert ev["abs_residual"] <= 3 + 1e-9
    assert abs(sum(ev["abs_residual"] for ev in residuals) - 8) < 1e-9
    assert max(ev["abs_residual"] for ev in residuals) == 1


def test_align_infeasible_trace_returns_conclusion():
    p = _payload(drift={"min": 5, "max": 6})
    r = client.post("/api/current-traces/align", json=p)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "infeasible"
    assert "no integer drift" in body["reason"]


def test_align_structurally_infeasible_trace():
    p = _payload(observations=[11, 12] * 15)  # 30 obs > 8 levels * 3 samples
    r = client.post("/api/current-traces/align", json=p)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "infeasible"
    assert "structurally infeasible" in body["reason"]


INVALID_OVERRIDES = [
    ("too few levels", {"levels": [1, 2, 3, 4, 5, 6, 7]}),
    ("too many levels", {"levels": list(range(25))}),
    ("too few observations", {"observations": [0] * 7}),
    ("too many observations", {"observations": [0] * 61}),
    ("non-integer level", {"levels": [10, 20, 30, 40, 50, 60, 70, 80.5]}),
    ("float drift bound", {"drift": {"min": 0.5, "max": 2}}),
    ("drift min above max", {"drift": {"min": 3, "max": -3}}),
    ("drift interval too wide", {"drift": {"min": 0, "max": 5000}}),
    ("negative residual limit", {"residual_limit": -1}),
    ("dwell min below 1", {"dwell": {"min": 0, "max": 2}}),
    ("dwell max above 3", {"dwell": {"min": 1, "max": 4}}),
    ("dwell min above max", {"dwell": {"min": 3, "max": 1}}),
    ("too many internal skips", {"max_internal_skips": 3}),
    ("unknown field", {"baseline_hint": 12}),
]


@pytest.mark.parametrize("override", [o for _, o in INVALID_OVERRIDES], ids=[n for n, _ in INVALID_OVERRIDES])
def test_invalid_requests_are_rejected(override):
    r = client.post("/api/current-traces/align", json=_payload(**override))
    assert r.status_code in (400, 422)
