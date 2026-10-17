"""Let a pipeline wrapper clean up its own detached Nextflow process tree."""
from __future__ import annotations

import signal
import subprocess
import threading
from collections.abc import Sequence

# The wrapper gives Nextflow ten seconds before escalating to SIGKILL. Leave
# enough time for that cleanup and for its captured output to drain.
_WRAPPER_SHUTDOWN_GRACE_SECONDS = 30


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
    stopping = False
    previous_handlers = {}
    main_thread = threading.current_thread() is threading.main_thread()

    def request_shutdown(signum, frame):
        nonlocal pending_signum, stopping
        pending_signum = signum
        if process is not None and not stopping:
            stopping = True
            raise _ShutdownRequested(signum)
        # A signal during spawn is handled as soon as Popen returns. Further
        # signals during cleanup must not interrupt the wrapper's grace period.

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
        if pending_signum is not None:
            raise _ShutdownRequested(pending_signum)
        stdout, stderr = process.communicate(timeout=timeout)
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    except (_ShutdownRequested, KeyboardInterrupt) as stop:
        stopping = True
        signum = stop.signum if isinstance(stop, _ShutdownRequested) else signal.SIGINT
        if process is not None:
            _stop_wrapper(process, signum)
        if signum == signal.SIGINT:
            # Preserve the public run_skill API's catchable cancellation.
            raise KeyboardInterrupt from None
        raise SystemExit(128 + signum) from None
    except BaseException:
        stopping = True
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
