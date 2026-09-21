"""Tests for geo-fetch.

Run with: pytest skills/geo-fetch/tests/test_geo_fetch.py -v

No network required: the vendored client's single HTTP entry point replays a
recorded URL->body map captured from one real GSE30720 run.
"""

import csv
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "geo_fetch.py"
DEMO = "GSE30720"

sys.path.insert(0, str(SKILL_DIR))


class TestCredentialGate:
    """'Presence of an env var is not consent.' Upstream's rule, kept verbatim."""

    def test_credentials_are_not_sent_by_default(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "me@example.org")
        monkeypatch.setattr(api, "NCBI_API_KEY", "secret")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", False)
        params = api._eutil_params({"db": "gds"})
        assert "email" not in params
        assert "api_key" not in params
        assert "tool" not in params

    def test_credentials_are_sent_only_when_opted_in(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "me@example.org")
        monkeypatch.setattr(api, "NCBI_API_KEY", "secret")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", True)
        params = api._eutil_params({"db": "gds"})
        assert params["email"] == "me@example.org"
        assert params["api_key"] == "secret"
        assert params["tool"] == "geo-skill"

    def test_opting_in_with_no_env_vars_sends_nothing(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "")
        monkeypatch.setattr(api, "NCBI_API_KEY", "")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", True)
        assert api._eutil_params({"db": "gds"}) == {"db": "gds"}

    def test_demo_reports_that_no_credentials_were_sent(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["summary"]["ncbi_credentials_sent"] is False


class TestVendoredApi:
    def test_ftp_paths_shard_the_accession(self):
        """GEO's FTP layout buckets by accession: GSE30720 -> series/GSE30nnn."""
        import geo_fetch_api as api

        assert api.geo_ftp_dir("GSE30720") == ("series", "GSE30nnn")
        assert api.geo_ftp_dir("GSM12345") == ("samples", "GSM12nnn")
        assert api.geo_ftp_dir("GSE123") == ("series", "GSEnnn")

    def test_clean_val_replaces_control_characters(self):
        import geo_fetch_api as api

        assert api._clean_val("a\tb\nc") == "a_b_c"

    def test_metadata_tsv_is_lf_terminated(self, tmp_path):
        import geo_fetch_api as api

        out = tmp_path / "metadata.tsv"
        api.write_metadata_tsv([api.harmonize_row("GSM1", "1", {})], out)
        assert b"\r\n" not in out.read_bytes()

    def test_read_map_rejects_malformed_values(self):
        import geo_fetch_api as api

        for bad in ("0,1", "x,2", "1", "1,2,3"):
            with pytest.raises(SystemExit):
                api.parse_read_map(bad)


class TestDemo:
    def test_demo_makes_no_network_call(self, tmp_path):
        import geo_fetch as app
        import geo_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            app.main(["--demo", "--output", str(tmp_path)])
            mocked.assert_not_called()

    def test_demo_writes_report_and_result(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()
        assert json.loads((tmp_path / "result.json").read_text())["skill"] == "geo-fetch"

    def test_demo_report_carries_the_disclaimer(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in text
        assert "not a medical device" in text

    def test_demo_writes_the_reproducibility_bundle(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (tmp_path / "reproducibility" / name).exists()

    def test_checksum_labels_resolve_from_the_output_dir(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for line in (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines():
            assert (tmp_path / line.split("  ", 1)[1]).exists()

    def test_demo_metadata_table_has_one_row_per_sample(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "tables" / "metadata.tsv").open(), delimiter="\t"))
        assert len(rows) == 42
        assert {r["species"] for r in rows} == {"Arabidopsis thaliana"}

    def test_demo_writes_a_pipeline_ready_samplesheet(self, tmp_path):
        """GEO resolves to its SRA project, then to ENA for the FASTQ links."""
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "samplesheet.csv").open()))
        assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "strandedness"]
        assert len(rows) == 42
        assert all(r["fastq_1"].startswith("http") for r in rows)

    def test_demo_download_script_downloads_nothing(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        body = (tmp_path / "download_geo.sh").read_text()
        assert "set -euo pipefail" in body
        assert "# #SBATCH --partition=<your_partition>" in body
        assert not (tmp_path / "fastq").exists()


class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
        assert result.returncode != 0

    def test_command_flag_dispatches(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--command", "metadata", "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()

    def test_upstream_positional_form_still_works(self, tmp_path):
        import geo_fetch as app

        app._install_demo_transport()
        app.main(["metadata", DEMO, "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()

    def test_samplesheet_requires_an_assay(self, tmp_path):
        import geo_fetch as app

        with pytest.raises(SystemExit):
            app.main(["--command", "samplesheet", "--accession", DEMO,
                      "--output", str(tmp_path)])


class TestSafety:
    def test_warns_before_overwriting(self, tmp_path, capsys):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwritten" in capsys.readouterr().err

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        """Upstream's `runtable --out .` defaulted to the working directory."""
        import geo_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(tmp_path / "out")])
        assert list(cwd.iterdir()) == []

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import geo_fetch as app

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
