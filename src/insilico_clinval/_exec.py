"""The single subprocess wrapper for external tools.

Every call checks the return code, captures stdout/stderr and measures wall time and the peak
resident set size of the child process (``os.wait4``). Pipeline tools are not called directly but
through ``insilico_clinval._runner.Toolbox``, which builds the runner command line (local or
docker) and records tool version and step in the run manifest.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Mapping
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger(__name__)


class ToolError(RuntimeError):
    """An external tool is missing or exited with a non-zero status."""


def elapsed(t0: float) -> str:
    """Wall time since ``t0`` (a ``time.perf_counter()`` value) as H:MM:SS, for progress logs."""
    s = int(time.perf_counter() - t0)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


@dataclass(frozen=True)
class ExecResult:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    started_at: str
    wall_s: float
    peak_rss_mb: float | None


def _maxrss_to_mb(ru_maxrss: int) -> float:
    # ru_maxrss is in bytes on macOS and in kilobytes on Linux.
    if sys.platform == "darwin":
        return ru_maxrss / (1024 * 1024)
    return ru_maxrss / 1024


def run(
    argv: list[str],
    *,
    cwd: Path | None = None,
    check: bool = True,
    stderr_path: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> ExecResult:
    """Run ``argv`` (no shell) and return its captured result.

    With ``stderr_path``, stderr is written there as it is produced (it can be followed while the
    tool runs and is kept when it fails) and is also returned in the result. ``env`` adds
    variables to the inherited environment.
    Raises ToolError if the executable is missing or, with ``check``, on a non-zero exit.
    """
    if not argv:
        raise ValueError("empty argv")
    if shutil.which(argv[0]) is None:
        raise ToolError(f"required tool not found on PATH: {argv[0]}")
    started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    t0 = time.perf_counter()
    with ExitStack() as stack:
        out_f = stack.enter_context(tempfile.TemporaryFile())
        if stderr_path is not None:
            err_f = stack.enter_context(stderr_path.open("w+b"))
        else:
            err_f = stack.enter_context(tempfile.TemporaryFile())
        child_env = {**os.environ, **env} if env else None
        proc = subprocess.Popen(argv, cwd=cwd, stdout=out_f, stderr=err_f, env=child_env)
        _, status, rusage = os.wait4(proc.pid, 0)
        returncode = os.waitstatus_to_exitcode(status)
        proc.returncode = returncode
        wall_s = time.perf_counter() - t0
        out_f.seek(0)
        err_f.seek(0)
        stdout = out_f.read().decode(errors="replace")
        stderr = err_f.read().decode(errors="replace")
    result = ExecResult(
        argv=list(argv),
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        started_at=started_at,
        wall_s=round(wall_s, 3),
        peak_rss_mb=round(_maxrss_to_mb(rusage.ru_maxrss), 1),
    )
    log.debug("exec %s -> %d (%.3f s)", " ".join(argv), returncode, wall_s)
    if check and returncode != 0:
        raise ToolError(
            f"{argv[0]} exited with status {returncode}: {' '.join(argv)}\n{stderr.strip()}"
        )
    return result
