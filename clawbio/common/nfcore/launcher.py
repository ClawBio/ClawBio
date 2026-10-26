"""Let a pipeline wrapper clean up its own detached Nextflow process tree."""
from __future__ import annotations

import signal
import subprocess
import threading
import time
from collections.abc import Callable, Sequence

# The wrapper gives Nextflow ten seconds before escalating to SIGKILL. Leave
# enough time for that cleanup and for its captured output to drain.
_WRAPPER_SHUTDOWN_GRACE_SECONDS = 30
# How often the capture checks for a recorded SIGINT/SIGTERM.
_STOP_POLL_SECONDS = 0.2


class _ShutdownRequested(BaseException):
    def __init__(self, signum: int) -> None:
        super().__init__(signum)
        self.signum = signum


def _stop_wrapper(process: subprocess.Popen, signum: int) -> None:
    if process.poll() is not None:
        return
    process.send_signal(signum)
    try:
        process.communicate(timeout=_WRAPPER_SHUTDOWN_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        process.kill()
        process.communicate()


def _communicate_unless_stopped(
    process: subprocess.Popen, timeout: int | None, stop_requested: Callable[[], bool]
) -> tuple[str | None, str | None]:
    """The wrapper's captured output, or (None, None) once a stop is requested.

    Communicates in short slices; a retried communicate() loses no output.
    """
    deadline = None if timeout is None else time.monotonic() + timeout
    while not stop_requested():
        remaining = None if deadline is None else deadline - time.monotonic()
        if remaining is not None and remaining <= 0:
            raise subprocess.TimeoutExpired(process.args, timeout)
        slice_seconds = _STOP_POLL_SECONDS if remaining is None else min(_STOP_POLL_SECONDS, remaining)
        try:
            return process.communicate(timeout=slice_seconds)
        except subprocess.TimeoutExpired:
            continue
    return None, None


def run_wrapper(
    command: Sequence[str], *, cwd: str, timeout: int | None
) -> subprocess.CompletedProcess:
    """Capture a POSIX wrapper while forwarding stops before waiting for cleanup.

    subprocess.run kills its immediate child on interruption or timeout. A wrapper
    must instead receive SIGINT/SIGTERM and stop Nextflow in its separate session.
    Main-thread calls isolate the wrapper from terminal signals so it receives
    each stop once. Worker-thread calls retain normal process-group membership
    and never install Python signal handlers.
    """
    process = None
    pending_signum = None
    previous_handlers = {}
    main_thread = threading.current_thread() is threading.main_thread()

    def request_shutdown(signum, frame):
        # Only record the first stop. Raising here could leave Popen's wait lock held
        # (see clawbio.common.nfcore.execution._stops_recorded), and later signals
        # must not interrupt the wrapper's cleanup grace period.
        nonlocal pending_signum
        if pending_signum is None:
            pending_signum = signum

    try:
        if main_thread:
            for signum in (signal.SIGTERM, signal.SIGINT):
                previous_handlers[signum] = signal.signal(signum, request_shutdown)
        if pending_signum is not None:
            raise _ShutdownRequested(pending_signum)
        process = subprocess.Popen(
            command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, start_new_session=main_thread,
        )
        stdout, stderr = _communicate_unless_stopped(process, timeout, lambda: pending_signum is not None)
        if pending_signum is not None:
            raise _ShutdownRequested(pending_signum)
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    except (_ShutdownRequested, KeyboardInterrupt) as stop:
        signum = stop.signum if isinstance(stop, _ShutdownRequested) else signal.SIGINT
        if process is not None:
            _stop_wrapper(process, signum)
        if signum == signal.SIGINT:
            # Preserve the public run_skill API's catchable cancellation.
            raise KeyboardInterrupt from None
        raise SystemExit(128 + signum) from None
    except BaseException:
        if process is not None:
            _stop_wrapper(process, signal.SIGTERM)
        raise
    finally:
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        if process is not None:
            for stream in (process.stdout, process.stderr):
                if stream is not None:
                    stream.close()
