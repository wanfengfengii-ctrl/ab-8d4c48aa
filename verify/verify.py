#!/usr/bin/env python3
"""One-shot verification for the nanopore current-trace aligner.

Steps, each reported on stdout; the process exit code is 0 only if every
step passed:

1. package build — ``pip install --no-build-isolation --no-deps .``
2. code tests    — ``pytest``
3. HTTP smoke    — a feasible trace aligns, an infeasible trace reports a
                   no-solution conclusion, and an invalid payload is
                   rejected (4xx)

The HTTP smoke targets ``API_BASE_URL`` (default ``http://localhost:8000``;
the docker-compose ``verify`` service sets it to ``http://api:8000``).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
API_BASE = os.environ.get("API_BASE_URL", "http://localhost:8000").rstrip("/")
ALIGN_URL = f"{API_BASE}/api/current-traces/align"
HEALTH_URL = f"{API_BASE}/health"

LEVELS = [10, 20, 30, 40, 50, 60, 70, 80]

# Level 40 (index 3) is deliberately unobserved: the optimal alignment must
# skip exactly that level with drift 1.
FEASIBLE_PAYLOAD = {
    "levels": LEVELS,
    "observations": [11, 12, 21, 22, 31, 32, 51, 52, 61, 62, 71, 72, 81, 82],
    "drift": {"min": -2, "max": 2},
    "residual_limit": 3,
    "dwell": {"min": 1, "max": 3},
    "max_internal_skips": 2,
}

# The drift interval is shifted so far that no residual can stay within the
# limit: no legal alignment exists.
INFEASIBLE_PAYLOAD = {
    "levels": LEVELS,
    "observations": [11, 12, 21, 22, 31, 32, 41, 42, 51, 52, 61, 62, 71, 72, 81, 82],
    "drift": {"min": 5, "max": 6},
    "residual_limit": 3,
    "dwell": {"min": 1, "max": 3},
    "max_internal_skips": 2,
}

# Only 3 reference levels (minimum is 8): must be rejected.
INVALID_PAYLOAD = {
    "levels": [1, 2, 3],
    "observations": [0] * 16,
    "drift": {"min": 0, "max": 1},
    "residual_limit": 1,
}


def _run(cmd: list[str], cwd: Path = ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=600)


def step_package_build() -> tuple[bool, str]:
    # Build from a scratch copy so the source tree stays free of build
    # artifacts (build/, *.egg-info, ...).
    with tempfile.TemporaryDirectory(prefix="aligner-build-") as tmp:
        tmp_path = Path(tmp)
        for name in ("pyproject.toml", "README.md"):
            shutil.copy2(ROOT / name, tmp_path / name)
        shutil.copytree(ROOT / "src", tmp_path / "src")
        proc = _run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-build-isolation",
                "--no-deps",
                "--force-reinstall",
                tmp,
            ],
            cwd=tmp_path,
        )
    tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-8:])
    return proc.returncode == 0, tail


def step_code_tests() -> tuple[bool, str]:
    proc = _run([sys.executable, "-m", "pytest", "tests", "-q"])
    tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-15:])
    return proc.returncode == 0, tail


def wait_for_api(timeout_s: int = 90) -> bool:
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(HEALTH_URL, timeout=3) as resp:
                if resp.status == 200:
                    return True
        except Exception:
            pass
        time.sleep(1)
    return False


def _post(payload: dict) -> tuple[int, dict]:
    req = urllib.request.Request(
        ALIGN_URL,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        raw = exc.read().decode() or "{}"
        try:
            return exc.code, json.loads(raw)
        except json.JSONDecodeError:
            return exc.code, {"raw": raw}


def step_smoke_feasible() -> tuple[bool, str]:
    status, body = _post(FEASIBLE_PAYLOAD)
    problems: list[str] = []
    if status != 200:
        problems.append(f"HTTP {status}")
    if body.get("status") != "aligned":
        problems.append(f"status={body.get('status')!r}")
        return False, "; ".join(problems)

    n = len(FEASIBLE_PAYLOAD["observations"])
    if body.get("drift") != 1:
        problems.append(f"drift={body.get('drift')} (expected 1)")
    if body.get("skips") != 1:
        problems.append(f"skips={body.get('skips')} (expected 1)")
    if body.get("skipped_level_indices") != [3]:
        problems.append(
            f"skipped_level_indices={body.get('skipped_level_indices')} (expected [3])"
        )

    segments = body.get("segments", [])
    if not segments or segments[0].get("start") != 0 or segments[-1].get("end") != n:
        problems.append("segments do not cover the whole trace")
    elif any(a["end"] != b["start"] for a, b in zip(segments, segments[1:])):
        problems.append("segments are not contiguous")
    if any(not (1 <= s.get("length", 0) <= 3) for s in segments):
        problems.append("dwell range violated")
    if segments and (segments[0].get("level_index") != 0 or segments[-1].get("level_index") != len(LEVELS) - 1):
        problems.append("first/last reference levels not used")

    residuals = body.get("residuals", [])
    if len(residuals) != n:
        problems.append(f"residual evidence has {len(residuals)} entries (expected {n})")
    else:
        drift = body.get("drift", 0)
        limit = FEASIBLE_PAYLOAD["residual_limit"]
        for ev in residuals:
            if ev.get("expected") != LEVELS[ev["level_index"]] + drift:
                problems.append(f"expected level mismatch at index {ev.get('index')}")
                break
            if abs(ev.get("residual", 0)) > limit + 1e-9:
                problems.append(f"residual limit violated at index {ev.get('index')}")
                break
        objective = body.get("objective", {})
        total = sum(abs(ev.get("residual", 0)) for ev in residuals)
        if abs(total - objective.get("sum_abs_residual", -1)) > 1e-6:
            problems.append("objective sum_abs_residual inconsistent with evidence")
        peak = max(abs(ev.get("residual", 0)) for ev in residuals)
        if abs(peak - objective.get("max_abs_residual", -1)) > 1e-6:
            problems.append("objective max_abs_residual inconsistent with evidence")

    detail = "; ".join(problems) if problems else "aligned with drift=1, skips=1"
    return not problems, detail


def step_smoke_infeasible() -> tuple[bool, str]:
    status, body = _post(INFEASIBLE_PAYLOAD)
    problems: list[str] = []
    if status != 200:
        problems.append(f"HTTP {status}")
    if body.get("status") != "infeasible":
        problems.append(f"status={body.get('status')!r} (expected 'infeasible')")
    if not body.get("reason"):
        problems.append("missing infeasibility reason")
    detail = "; ".join(problems) if problems else "no-solution conclusion returned"
    return not problems, detail


def step_smoke_invalid() -> tuple[bool, str]:
    status, _body = _post(INVALID_PAYLOAD)
    ok = status in (400, 422)
    detail = f"HTTP {status}" if ok else f"HTTP {status} (expected 4xx)"
    return ok, detail


def main() -> int:
    report: list[tuple[str, bool, str]] = []

    def record(name: str, ok: bool, detail: str = "") -> None:
        report.append((name, ok, detail))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
        if not ok and detail:
            print(textwrap.indent(detail, "      "), flush=True)

    ok, detail = step_package_build()
    record("package build", ok, detail)

    ok, detail = step_code_tests()
    record("code tests", ok, detail)

    smoke_steps = [
        ("http smoke: feasible trace", step_smoke_feasible),
        ("http smoke: infeasible trace", step_smoke_infeasible),
        ("http smoke: invalid payload rejected", step_smoke_invalid),
    ]
    if wait_for_api():
        record("api health", True)
        for name, fn in smoke_steps:
            try:
                ok, detail = fn()
            except Exception as exc:  # noqa: BLE001 - report, never crash
                ok, detail = False, repr(exc)
            record(name, ok, detail)
    else:
        record("api health", False, f"no response at {HEALTH_URL}")
        for name, _ in smoke_steps:
            record(name, False, "api unavailable")

    failed = [name for name, ok, _ in report if not ok]
    total = len(report)
    print(
        f"VERIFY RESULT: {'PASS' if not failed else 'FAIL'} "
        f"({total - len(failed)}/{total} steps ok)",
        flush=True,
    )
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
