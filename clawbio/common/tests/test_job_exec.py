"""Tests for the shared run/submit helper.

Run with: pytest clawbio/common/tests/test_job_exec.py -v

No script is ever executed by these tests; subprocess is mocked throughout.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from clawbio.common.job_exec import (
    ExecutionRefused,
    preflight,
    run_scripts,
    submit_scripts,
)


def _script(tmp_path, name="dl.sh"):
    path = tmp_path / name
    path.write_text("#!/bin/bash\necho hi\n")
    path.chmod(0o755)
    return path


class TestPreflight:
    def test_passes_for_a_readable_script_with_its_tools_present(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/wget"):
            assert preflight([script], require=("wget",)) == []

    def test_reports_a_missing_script(self, tmp_path):
        problems = preflight([tmp_path / "absent.sh"])
        assert any("absent.sh" in p for p in problems)

    def test_reports_a_missing_required_tool(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", return_value=None):
            problems = preflight([script], require=("wget",))
        assert any("wget" in p for p in problems)

    def test_reports_missing_sbatch_when_submitting(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", side_effect=lambda n: None if n == "sbatch" else "/bin/x"):
            problems = preflight([script], submit=True)
        assert any("sbatch" in p for p in problems)

    def test_reports_an_unwritable_target_directory(self, tmp_path):
        script = _script(tmp_path)
        target = tmp_path / "ro"
        target.mkdir()
        target.chmod(0o500)
        try:
            problems = preflight([script], write_targets=[target])
            assert any("writable" in p for p in problems)
        finally:
            target.chmod(0o700)


class TestRun:
    def test_refuses_to_run_when_preflight_fails(self, tmp_path):
        """Fail before the expensive step, not on the compute node."""
        with pytest.raises(ExecutionRefused) as excinfo:
            run_scripts([tmp_path / "absent.sh"])
        assert "absent.sh" in str(excinfo.value)

    def test_runs_each_script_in_order(self, tmp_path):
        a, b = _script(tmp_path, "a.sh"), _script(tmp_path, "b.sh")
        with patch("subprocess.run", return_value=MagicMock(returncode=0)) as run:
            results = run_scripts([a, b])
        assert [c.args[0][1] for c in run.call_args_list] == [str(a), str(b)]
        assert all(r.returncode == 0 for r in results)

    def test_stops_at_the_first_failure(self, tmp_path):
        """A failed prefetch must not be followed by a dump over nothing."""
        a, b = _script(tmp_path, "a.sh"), _script(tmp_path, "b.sh")
        with patch("subprocess.run", return_value=MagicMock(returncode=1)) as run:
            results = run_scripts([a, b])
        assert len(results) == 1 and run.call_count == 1


class TestSubmit:
    def test_refuses_without_sbatch(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", return_value=None):
            with pytest.raises(ExecutionRefused, match="sbatch"):
                submit_scripts([script])

    def test_chains_jobs_with_afterok(self, tmp_path):
        """Step two must not start until step one has actually succeeded."""
        a, b = _script(tmp_path, "a.sh"), _script(tmp_path, "b.sh")
        with patch("shutil.which", return_value="/usr/bin/sbatch"), \
             patch("subprocess.run", side_effect=[
                 MagicMock(returncode=0, stdout="101\n"),
                 MagicMock(returncode=0, stdout="102\n"),
             ]) as run:
            ids = submit_scripts([a, b])
        assert ids == ["101", "102"]
        first, second = (c.args[0] for c in run.call_args_list)
        assert "--dependency" not in " ".join(first)
        assert "--dependency=afterok:101" in second

    def test_uses_parsable_so_the_job_id_can_be_read(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/sbatch"), \
             patch("subprocess.run", return_value=MagicMock(returncode=0, stdout="7\n")) as run:
            submit_scripts([script])
        assert "--parsable" in run.call_args_list[0].args[0]

    def test_raises_when_a_submission_is_rejected(self, tmp_path):
        script = _script(tmp_path)
        with patch("shutil.which", return_value="/usr/bin/sbatch"), \
             patch("subprocess.run", return_value=MagicMock(
                 returncode=1, stdout="", stderr="invalid partition specified")):
            with pytest.raises(ExecutionRefused, match="invalid partition"):
                submit_scripts([script])
