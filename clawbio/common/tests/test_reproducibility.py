"""Tests for clawbio.common.reproducibility."""

import json
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from clawbio.common.checksums import sha256_file
from clawbio.common.reproducibility import (
    write_checksums,
    write_environment_yml,
    write_commands_sh,
    write_conda_lock,
    write_ro_crate,
)


# ---------------------------------------------------------------------------
# TestWriteChecksums
# ---------------------------------------------------------------------------


class TestWriteChecksums:
    def test_creates_checksums_file(self, tmp_path):
        """write_checksums creates reproducibility/checksums.sha256."""
        f = tmp_path / "out.csv"
        f.write_bytes(b"data")
        write_checksums([f], tmp_path)
        assert (tmp_path / "reproducibility" / "checksums.sha256").exists()

    def test_returns_path(self, tmp_path):
        """Returns the path of the written checksums file."""
        f = tmp_path / "out.csv"
        f.write_bytes(b"data")
        result = write_checksums([f], tmp_path)
        assert result == tmp_path / "reproducibility" / "checksums.sha256"

    def test_creates_reproducibility_dir_if_missing(self, tmp_path):
        """reproducibility/ is created automatically."""
        f = tmp_path / "out.csv"
        f.write_bytes(b"data")
        assert not (tmp_path / "reproducibility").exists()
        write_checksums([f], tmp_path)
        assert (tmp_path / "reproducibility").exists()

    def test_format_filename_only_by_default(self, tmp_path):
        """Default (no anchor): line is '<hash>  <filename>'."""
        f = tmp_path / "masks.tif"
        f.write_bytes(b"pixels")
        write_checksums([f], tmp_path)
        line = (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines()[0]
        assert line == f"{sha256_file(f)}  {f.name}"

    def test_format_relative_path_with_anchor(self, tmp_path):
        """With anchor, line uses path relative to anchor."""
        sub = tmp_path / "figures"
        sub.mkdir()
        f = sub / "plot.png"
        f.write_bytes(b"png")
        write_checksums([f], tmp_path, anchor=tmp_path)
        line = (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines()[0]
        assert line == f"{sha256_file(f)}  figures/plot.png"

    def test_multiple_files_produce_multiple_lines(self, tmp_path):
        """Each file gets its own line."""
        paths = []
        for name in ("a.tif", "b.csv", "c.npy"):
            p = tmp_path / name
            p.write_bytes(name.encode())
            paths.append(p)
        write_checksums(paths, tmp_path)
        lines = (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines()
        assert len(lines) == 3

    def test_skips_missing_files(self, tmp_path):
        """Files that don't exist are silently skipped."""
        real = tmp_path / "real.csv"
        real.write_bytes(b"x")
        ghost = tmp_path / "ghost.tif"
        write_checksums([real, ghost], tmp_path)
        text = (tmp_path / "reproducibility" / "checksums.sha256").read_text()
        assert "real.csv" in text
        assert "ghost.tif" not in text

    def test_digest_matches_sha256_file(self, tmp_path):
        """Digest written matches sha256_file from commons."""
        f = tmp_path / "verify.bin"
        f.write_bytes(b"clawbio")
        write_checksums([f], tmp_path)
        written = (tmp_path / "reproducibility" / "checksums.sha256").read_text().strip()
        assert written.split("  ")[0] == sha256_file(f)

    def test_empty_list_produces_empty_file(self, tmp_path):
        """Empty path list writes an empty checksums file without error."""
        write_checksums([], tmp_path)
        assert (tmp_path / "reproducibility" / "checksums.sha256").read_text() == ""


# ---------------------------------------------------------------------------
# TestWriteEnvironmentYml
# ---------------------------------------------------------------------------


class TestWriteEnvironmentYml:
    def test_creates_environment_yml(self, tmp_path):
        """write_environment_yml creates reproducibility/environment.yml."""
        write_environment_yml(tmp_path, "clawbio-test", ["numpy", "scipy"])
        assert (tmp_path / "reproducibility" / "environment.yml").exists()

    def test_returns_path(self, tmp_path):
        """Returns the path of the written environment.yml."""
        result = write_environment_yml(tmp_path, "clawbio-test", ["numpy"])
        assert result == tmp_path / "reproducibility" / "environment.yml"

    def test_creates_reproducibility_dir_if_missing(self, tmp_path):
        """reproducibility/ is created automatically."""
        assert not (tmp_path / "reproducibility").exists()
        write_environment_yml(tmp_path, "clawbio-test", ["numpy"])
        assert (tmp_path / "reproducibility").exists()

    def test_name_appears_in_file(self, tmp_path):
        """The env name is written as 'name: <env_name>'."""
        write_environment_yml(tmp_path, "clawbio-cell-detection", ["numpy"])
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert "name: clawbio-cell-detection" in text

    def test_pip_deps_appear_in_file(self, tmp_path):
        """Each pip dependency appears in the file."""
        write_environment_yml(tmp_path, "clawbio-test", ["cellpose>=4.0", "tifffile"])
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert "cellpose>=4.0" in text
        assert "tifffile" in text

    def test_channels_present(self, tmp_path):
        """Standard conda channels are included."""
        write_environment_yml(tmp_path, "clawbio-test", ["numpy"])
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert "conda-forge" in text

    def test_conda_deps_separate_from_pip(self, tmp_path):
        """conda_deps kwarg lists packages outside the pip block."""
        write_environment_yml(tmp_path, "clawbio-test", ["cellpose>=4.0"],
                              conda_deps=["numpy"])
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert "numpy" in text

    def test_empty_pip_deps_produces_valid_yaml(self, tmp_path):
        """Empty pip_deps must not emit 'pip: null' — should omit the pip block."""
        import yaml
        write_environment_yml(tmp_path, "clawbio-test", pip_deps=[])
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        parsed = yaml.safe_load(text)
        deps = parsed["dependencies"]
        for item in deps:
            if isinstance(item, dict):
                assert item.get("pip") is not None, "pip key must not be null"

    def test_python_version_not_duplicated_when_in_conda_deps(self, tmp_path):
        """python= in conda_deps must not produce a duplicate python= line."""
        write_environment_yml(tmp_path, "clawbio-test", ["cellpose"],
                              conda_deps=["python=3.11", "numpy"],
                              python_version="3.11")
        text = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert text.count("python=3.11") == 1


# ---------------------------------------------------------------------------
# TestWriteCommandsSh
# ---------------------------------------------------------------------------


class TestWriteCommandsSh:
    def test_creates_commands_sh(self, tmp_path):
        """write_commands_sh creates reproducibility/commands.sh."""
        write_commands_sh(tmp_path, "python skill.py --demo")
        assert (tmp_path / "reproducibility" / "commands.sh").exists()

    def test_returns_path(self, tmp_path):
        """Returns the path of the written commands.sh."""
        result = write_commands_sh(tmp_path, "python skill.py --demo")
        assert result == tmp_path / "reproducibility" / "commands.sh"

    def test_creates_reproducibility_dir_if_missing(self, tmp_path):
        """reproducibility/ is created automatically."""
        assert not (tmp_path / "reproducibility").exists()
        write_commands_sh(tmp_path, "python skill.py")
        assert (tmp_path / "reproducibility").exists()

    def test_command_appears_in_file(self, tmp_path):
        """The command string appears verbatim in the file."""
        cmd = "python skills/cell-detection/cell_detection.py --demo --output /tmp/out"
        write_commands_sh(tmp_path, cmd)
        text = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert cmd in text

    def test_has_shebang(self, tmp_path):
        """File starts with a bash shebang."""
        write_commands_sh(tmp_path, "python skill.py")
        text = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert text.startswith("#!/")

    def test_is_executable_bash(self, tmp_path):
        """Shebang targets bash or env bash."""
        write_commands_sh(tmp_path, "python skill.py")
        first_line = (tmp_path / "reproducibility" / "commands.sh").read_text().splitlines()[0]
        assert "bash" in first_line

    def test_multiline_command(self, tmp_path):
        """Multi-line commands are written intact."""
        cmd = "python skill.py \\\n  --input data.csv \\\n  --output /tmp/out"
        write_commands_sh(tmp_path, cmd)
        text = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert cmd in text

    def test_file_is_executable(self, tmp_path):
        """commands.sh must have the executable bit set."""
        import stat
        write_commands_sh(tmp_path, "python skill.py")
        mode = (tmp_path / "reproducibility" / "commands.sh").stat().st_mode
        assert mode & stat.S_IXUSR, "owner execute bit not set"


# ---------------------------------------------------------------------------
# TestWriteCondaLock
# ---------------------------------------------------------------------------


class TestWriteCondaLock:
    def test_calls_conda_lock_with_correct_args(self, tmp_path):
        """write_conda_lock calls conda-lock with correct arguments."""
        repro_dir = tmp_path / "reproducibility"
        repro_dir.mkdir()
        (repro_dir / "environment.yml").write_text("name: test\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            write_conda_lock(tmp_path)

        mock_run.assert_called_once_with(
            ["conda-lock", "lock", "-f", "environment.yml"],
            cwd=repro_dir,
            check=True,
        )

    def test_returns_lockfile_path(self, tmp_path):
        """write_conda_lock returns the path to conda-lock.yml."""
        repro_dir = tmp_path / "reproducibility"
        repro_dir.mkdir()
        (repro_dir / "environment.yml").write_text("name: test\n")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            result = write_conda_lock(tmp_path)

        assert result == repro_dir / "conda-lock.yml"

    def test_raises_if_environment_yml_missing(self, tmp_path):
        """write_conda_lock raises FileNotFoundError if environment.yml is missing."""
        repro_dir = tmp_path / "reproducibility"
        repro_dir.mkdir()

        with pytest.raises(FileNotFoundError):
            write_conda_lock(tmp_path)

    def test_raises_clear_error_if_conda_lock_not_installed(self, tmp_path):
        """write_conda_lock raises FileNotFoundError with install instructions if conda-lock binary is missing."""
        repro_dir = tmp_path / "reproducibility"
        repro_dir.mkdir()
        (repro_dir / "environment.yml").write_text("name: test\n")

        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = FileNotFoundError("conda-lock")
            with pytest.raises(FileNotFoundError, match="pip install conda-lock"):
                write_conda_lock(tmp_path)

    def test_propagates_called_process_error(self, tmp_path):
        """write_conda_lock propagates CalledProcessError from subprocess.run."""
        repro_dir = tmp_path / "reproducibility"
        repro_dir.mkdir()
        (repro_dir / "environment.yml").write_text("name: test\n")

        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "conda-lock")
            with pytest.raises(subprocess.CalledProcessError):
                write_conda_lock(tmp_path)


# ---------------------------------------------------------------------------
# TestWriteRoCrate
# ---------------------------------------------------------------------------


class TestWriteRoCrate:
    def test_creates_metadata_file(self, tmp_path):
        out = write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
        )
        assert out == tmp_path / "ro-crate-metadata.json"
        assert out.exists()

    def test_metadata_is_valid_json(self, tmp_path):
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
        )
        data = json.loads((tmp_path / "ro-crate-metadata.json").read_text())
        assert "@context" in data
        assert "@graph" in data

    def test_root_dataset_contains_skill_version(self, tmp_path):
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="1.2.3",
            script_path="skills/test-skill/test_skill.py",
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        root = next(e for e in graph if e.get("@id") == "./")
        assert root.get("version") == "1.2.3"

    def test_create_action_present(self, tmp_path):
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        actions = [e for e in graph if e.get("@type") == "CreateAction"]
        assert len(actions) == 1
        assert actions[0]["instrument"]["@id"] == "skills/test-skill/test_skill.py"

    def test_params_included_as_property_values(self, tmp_path):
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
            params={"input": "sample.vcf", "threshold": "0.05"},
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        pv = [e for e in graph if e.get("@type") == "PropertyValue"]
        names = {e["name"] for e in pv}
        assert {"input", "threshold"} == names

    def test_output_files_registered_as_has_part(self, tmp_path):
        (tmp_path / "report.md").write_text("# Report")
        (tmp_path / "figure.png").write_bytes(b"PNG")
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        root = next(e for e in graph if e.get("@id") == "./")
        has_part_ids = {e["@id"] for e in root.get("hasPart", [])}
        assert "report.md" in has_part_ids
        assert "figure.png" in has_part_ids

    def test_description_propagated(self, tmp_path):
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
            description="Demo run for testing",
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        root = next(e for e in graph if e.get("@id") == "./")
        assert root.get("description") == "Demo run for testing"

    def test_completed_at_override(self, tmp_path):
        ts = "2026-01-01T00:00:00+00:00"
        write_ro_crate(
            tmp_path,
            skill_name="test-skill",
            skill_version="0.1.0",
            script_path="skills/test-skill/test_skill.py",
            completed_at=ts,
        )
        graph = json.loads((tmp_path / "ro-crate-metadata.json").read_text())["@graph"]
        action = next(e for e in graph if e.get("@type") == "CreateAction")
        assert action["endTime"] == ts
