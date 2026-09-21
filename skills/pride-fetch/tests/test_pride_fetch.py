"""Tests for pride-fetch.

Run with: pytest skills/pride-fetch/tests/test_pride_fetch.py -v

No network required: the vendored client's single HTTP entry point is served
from the committed fixtures in ../examples/.
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "pride_fetch.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO = "PXD084218"

sys.path.insert(0, str(SKILL_DIR))


class TestVendoredApi:
    def test_sdrf_out_path_enforces_the_extension(self):
        """quantms rejects .sdrf / .tsv / .csv; only .sdrf.tsv is accepted."""
        import pride_fetch_api as api

        assert str(api.sdrf_out_path("x.sdrf")).endswith(".sdrf.tsv")
        assert str(api.sdrf_out_path("x.tsv")).endswith(".sdrf.tsv")
        assert str(api.sdrf_out_path("x.sdrf.tsv")).endswith(".sdrf.tsv")

    def test_is_ms_file_recognises_acquisition_formats(self):
        import pride_fetch_api as api

        assert api.is_ms_file("run1.raw")
        assert api.is_ms_file("run1.mzML")
        assert not api.is_ms_file("proteins.fasta")
        assert not api.is_ms_file("README.txt")

    def test_ftp_urls_are_rewritten_to_https(self):
        import pride_fetch_api as api

        assert api._https("ftp://ftp.pride.ebi.ac.uk/a/b.raw") == \
            "https://ftp.pride.ebi.ac.uk/a/b.raw"

    def test_clean_val_replaces_control_characters(self):
        import pride_fetch_api as api

        assert api._clean_val("a\tb\nc") == "a_b_c"

    def test_minimal_sdrf_columns_are_the_documented_nineteen(self):
        import pride_fetch_api as api

        assert len(api.MINIMAL_SDRF_COLUMNS) == 19


class TestMetadataTableFromSdrf:
    """The SDRF-driven path cannot be demoed (submitter SDRFs live on
    ftp.pride.ebi.ac.uk), so it is covered here with a synthetic SDRF."""

    def test_sdrf_rows_become_harmonised_rows(self, tmp_path):
        import pride_fetch_api as api

        sdrf = (
            "source name\tcharacteristics[organism]\tcharacteristics[disease]\t"
            "comment[technical replicate]\n"
            "S1\tArabidopsis thaliana\tnormal\t1\n"
            "S2\tArabidopsis thaliana\tnormal\t2\n")
        out = tmp_path / "metadata.tsv"
        args = SimpleNamespace(accession="PXD1", out=str(out))

        with patch.object(api, "get_json", return_value=["ftp://host/x.sdrf.tsv"]), \
             patch.object(api, "http_get", return_value=sdrf.encode()), \
             patch.object(api, "merge_biosample", side_effect=lambda s, a: a):
            api.cmd_metadata_table(args)

        rows = list(csv.DictReader(out.open(), delimiter="\t"))
        assert [r["sample"] for r in rows] == ["S1", "S2"]
        assert {r["species"] for r in rows} == {"Arabidopsis thaliana"}
        assert [r["replicate"] for r in rows] == ["1", "2"]


class TestDemo:
    def test_demo_makes_no_network_call(self, tmp_path):
        import pride_fetch as app
        import pride_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            app.main(["--demo", "--output", str(tmp_path)])
            mocked.assert_not_called()

    def test_demo_writes_report_and_result(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()
        assert json.loads((tmp_path / "result.json").read_text())["skill"] == "pride-fetch"

    def test_demo_report_carries_the_disclaimer(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in text
        assert "not a medical device" in text

    def test_demo_writes_the_reproducibility_bundle(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (tmp_path / "reproducibility" / name).exists()

    def test_checksum_labels_resolve_from_the_output_dir(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for line in (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines():
            assert (tmp_path / line.split("  ", 1)[1]).exists()

    def test_demo_generates_a_minimal_sdrf(self, tmp_path):
        """The demo project has no submitter SDRF, so `auto` generates one."""
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        sdrf = tmp_path / f"{DEMO}.sdrf.tsv"
        assert sdrf.exists()
        header = sdrf.read_text().splitlines()[0].split("\t")
        assert header[0] == "source name"
        assert "characteristics[organism]" in header

    def test_demo_metadata_table_falls_back_to_project_level(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "tables" / "metadata.tsv").open(), delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["species"] == "Arabidopsis thaliana"

    def test_demo_download_script_downloads_nothing(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        script = tmp_path / "download_pride.sh"
        assert script.exists()
        body = script.read_text()
        assert "set -euo pipefail" in body
        assert "# #SBATCH --partition=<your_partition>" in body
        assert not (tmp_path / "pride_data").exists()


class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
        assert result.returncode != 0

    def test_command_flag_dispatches(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--command", "metadata", "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()

    def test_upstream_positional_form_still_works(self, tmp_path):
        import pride_fetch as app

        app._install_demo_transport()
        app.main(["metadata", DEMO, "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()


class TestSafety:
    def test_warns_before_overwriting(self, tmp_path, capsys):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwritten" in capsys.readouterr().err

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        import pride_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(tmp_path / "out")])
        assert list(cwd.iterdir()) == []

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import pride_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert str(tmp_path) not in (tmp_path / "report.md").read_text()


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree."""
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files, parents = [], {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' must be produced."""

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True)
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing))
