"""Tests for pride-fetch.

Run with: pytest skills/pride-fetch/tests/test_pride_fetch.py -v

No network required: the vendored client's single HTTP entry point is served
from the committed fixtures in ../examples/.
"""

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

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


class TestSdrfPlaceholders:
    """Defaults that are wrong for some rows must be `not available`, not a guess."""

    META = {"organisms": [{"name": "Homo sapiens (human)"}],
            "diseases": [{"name": "Alzheimer disease"}, {"name": "normal"}]}

    def test_project_disease_is_not_stamped_on_every_row(self, capsys):
        """A case/control project's first disease would label the controls too."""
        import pride_fetch_api as api

        with patch.object(api, "get_json", return_value=self.META):
            d = api.minimal_defaults("PXD0", None)
        assert d["characteristics[disease]"] == "not available"
        assert d["characteristics[organism]"] == "Homo sapiens"
        err = capsys.readouterr().err
        assert "Alzheimer disease" in err and "normal" in err

    @pytest.mark.parametrize("acq, expected", [
        (None, "not available"),
        ("dia", "data-independent acquisition"),
        ("dda", "data-dependent acquisition"),
    ])
    def test_acquisition_is_only_set_when_given(self, acq, expected):
        import pride_fetch_api as api

        with patch.object(api, "get_json", return_value={}):
            d = api.minimal_defaults("PXD0", acq)
        assert d["comment[proteomics data acquisition method]"] == expected

    def test_a_submitter_sdrf_is_not_completed_with_a_guessed_acquisition(self, tmp_path):
        import pride_fetch_api as api

        with patch.object(api, "get_json", return_value={}):
            defaults = api.minimal_defaults("PXD0", None)
        out = tmp_path / "x.sdrf.tsv"
        api.complete_existing_sdrf("source name\tcomment[data file]\nS1\ta.raw\n",
                                   str(out), defaults)
        rows = list(csv.DictReader(out.open(), delimiter="\t"))
        assert rows[0]["comment[proteomics data acquisition method]"] == "not available"

    def test_the_entry_point_does_not_default_the_acquisition(self, tmp_path):
        import pride_fetch as app

        args = app._build_parser().parse_args(
            ["--command", "samplesheet", "--accession", "PXD0"])
        assert "--acquisition" not in app._to_upstream_argv(args, tmp_path)
        args = app._build_parser().parse_args(
            ["--command", "samplesheet", "--accession", "PXD0", "--acquisition", "dda"])
        argv = app._to_upstream_argv(args, tmp_path)
        assert argv[argv.index("--acquisition") + 1] == "dda"


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

    def test_download_script_mode_honours_the_umask(self, tmp_path):
        """Execute is added only where read already is: a 077 umask gives 0700."""
        import pride_fetch as app

        old = os.umask(0o077)
        try:
            app.main(["--demo", "--output", str(tmp_path)])
        finally:
            os.umask(old)
        mode = (tmp_path / "download_pride.sh").stat().st_mode & 0o777
        assert mode == 0o700

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        import pride_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(tmp_path / "out")])
        assert list(cwd.iterdir()) == []

    @pytest.mark.parametrize("tool", ["curl", "wget"])
    def test_download_script_passes_hostile_names_literally(self, tmp_path, tool):
        """PRIDE keeps its own emitter, so the shell-quoting fix is tested here too:
        a file name with `$(...)` or backticks must reach curl/wget/unzip as text."""
        import pride_fetch_api as api
        from clawbio.common.tests.shims import run_with_shims

        names = ["a$(touch PWNED_SUBST).raw", "b`touch PWNED_TICK`.zip"]
        recs = [{"publicFileLocations": [
            {"name": "FTP Protocol", "value": f"ftp://ftp.pride.ebi.ac.uk/p/{n}"}]}
            for n in names]
        script = tmp_path / "dl.sh"
        args = api.argparse.Namespace(
            accession="PXD0", ext=None, tool=tool, outdir="pride $(touch PWNED_OUT)",
            out=str(script), unzip=True, no_slurm=True)
        with patch.object(api, "iter_files", return_value=recs):
            api.cmd_download_script(args)

        result, calls = run_with_shims(script, tmp_path)
        assert result.returncode == 0, result.stderr
        assert not list(tmp_path.glob("PWNED*"))
        transfers = [c for c in calls if "--help" not in c and c[:2] != ["-q", "-o"]]
        assert [c[-1] for c in transfers] == [
            f"https://ftp.pride.ebi.ac.uk/p/{n}" for n in names]
        unzip = [c for c in calls if c[:2] == ["-q", "-o"]]
        assert unzip == [["-q", "-o", f"pride $(touch PWNED_OUT)/{names[1]}",
                          "-d", "pride $(touch PWNED_OUT)"]]

    UNSAFE_LOCATIONS = [
        "-o/etc/passwd",
        "--config=evil.cfg",
        "file:///etc/passwd",
        "prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2026/01/PXD0/a.raw",
    ]

    @pytest.mark.parametrize("value", UNSAFE_LOCATIONS)
    def test_download_script_emits_only_https_urls(self, tmp_path, capsys, value):
        """ftp_url() falls back to the first location of any kind, so a value that
        starts with `-` would reach curl/wget as an option; shlex.quote does not
        stop that. Only https (after the ftp rewrite) may be emitted."""
        import pride_fetch_api as api

        recs = [
            {"publicFileLocations": [{"name": "Other", "value": value}]},
            {"publicFileLocations": [
                {"name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/p/ok.raw"}]},
        ]
        script = tmp_path / "dl.sh"
        args = api.argparse.Namespace(
            accession="PXD0", ext=None, tool="curl", outdir="pride",
            out=str(script), unzip=False, no_slurm=True)
        with patch.object(api, "iter_files", return_value=recs):
            api.cmd_download_script(args)
        body = script.read_text()
        assert value not in body
        assert "https://ftp.pride.ebi.ac.uk/p/ok.raw" in body
        assert "skipping" in capsys.readouterr().err

    def test_download_script_refuses_a_project_with_no_https_location(self, tmp_path):
        import pride_fetch_api as api

        recs = [{"publicFileLocations": [{"name": "Other", "value": "-o/etc/passwd"}]}]
        args = api.argparse.Namespace(
            accession="PXD0", ext=None, tool="curl", outdir="pride",
            out=str(tmp_path / "dl.sh"), unzip=False, no_slurm=True)
        with patch.object(api, "iter_files", return_value=recs):
            with pytest.raises(SystemExit):
                api.cmd_download_script(args)
        assert not (tmp_path / "dl.sh").exists()

    @pytest.mark.parametrize("url", ["file:///etc/passwd", "-x", "data:,hello"])
    def test_download_file_refuses_non_http_schemes(self, tmp_path, url):
        """urllib opens file:// and data: URLs happily; a server-supplied
        location must never read a local file into the output directory."""
        import pride_fetch_api as api

        def no_network(*a, **k):
            raise AssertionError("urlopen must not be reached")

        with patch.object(api.urllib.request, "urlopen", no_network):
            with pytest.raises(SystemExit):
                api.download_file(url, str(tmp_path))
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.parametrize("field", ["job_name", "partition", "account", "cpus",
                                       "mem", "time", "email"])
    def test_slurm_header_values_cannot_break_out_of_their_line(self, tmp_path, field):
        import pride_fetch_api as api

        recs = [{"publicFileLocations": [
            {"name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/p/ok.raw"}]}]
        opts = dict(job_name="j", partition="p", account="a", cpus="1", mem="4G",
                    time="1:00:00", email="me@example.org")
        opts[field] = "x\nrm -rf ~\r\n#SBATCH --wrap=evil"
        script = tmp_path / "dl.sh"
        args = api.argparse.Namespace(
            accession="PXD0", ext=None, tool="curl", outdir="pride",
            out=str(script), unzip=False, no_slurm=False, **opts)
        with patch.object(api, "iter_files", return_value=recs):
            api.cmd_download_script(args)
        lines = script.read_text().splitlines()
        assert not any(line.startswith(("rm", "#SBATCH --wrap")) for line in lines)

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


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/pride-fetch/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("pride-fetch", tmp_path) == snapshot({'PXD084218.sdrf.tsv': """\
source name	characteristics[organism]	characteristics[organism part]	characteristics[disease]	characteristics[biological replicate]	assay name	technology type	comment[technical replicate]	comment[data file]	comment[file uri]	comment[fraction identifier]	comment[label]	comment[instrument]	comment[proteomics data acquisition method]	comment[cleavage agent details]	comment[modification parameters]	comment[precursor mass tolerance]	comment[fragment mass tolerance]	factor value[condition]
Sample 1	Arabidopsis thaliana	not available	not available	1	run 1	proteomic profiling by mass spectrometry	1	E20260218-08.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-08.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 2	Arabidopsis thaliana	not available	not available	1	run 2	proteomic profiling by mass spectrometry	1	E20260218-16.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-16.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 3	Arabidopsis thaliana	not available	not available	1	run 3	proteomic profiling by mass spectrometry	1	E20260218-12.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-12.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 4	Arabidopsis thaliana	not available	not available	1	run 4	proteomic profiling by mass spectrometry	1	E20260218-07.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-07.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 5	Arabidopsis thaliana	not available	not available	1	run 5	proteomic profiling by mass spectrometry	1	E20260218-04.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-04.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 6	Arabidopsis thaliana	not available	not available	1	run 6	proteomic profiling by mass spectrometry	1	E20260218-11.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-11.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 7	Arabidopsis thaliana	not available	not available	1	run 7	proteomic profiling by mass spectrometry	1	E20260218-03.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-03.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 8	Arabidopsis thaliana	not available	not available	1	run 8	proteomic profiling by mass spectrometry	1	E20260218-02.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-02.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 9	Arabidopsis thaliana	not available	not available	1	run 9	proteomic profiling by mass spectrometry	1	E20260218-06.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-06.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 10	Arabidopsis thaliana	not available	not available	1	run 10	proteomic profiling by mass spectrometry	1	E20260218-10.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-10.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 11	Arabidopsis thaliana	not available	not available	1	run 11	proteomic profiling by mass spectrometry	1	E20260218-14.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-14.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
Sample 12	Arabidopsis thaliana	not available	not available	1	run 12	proteomic profiling by mass spectrometry	1	E20260218-15.raw	https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-15.raw	1	label free sample	NT=Orbitrap Eclipse;AC=MS:1003029	not available	NT=Trypsin;AC=MS:1001251	NT=Carbamidomethyl;AC=UNIMOD:4;MT=Fixed;PP=Anywhere;TA=C	10 ppm	20 ppm	not available
""", 'download_pride.sh': """\
#!/bin/bash
#
#SBATCH --job-name=pride_download
# #SBATCH --partition=<your_partition>   # set for your cluster
# #SBATCH --account=<your_account>        # set if required
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=<you@example.org>

set -euo pipefail

OUTDIR=pride_data
mkdir -p "$OUTDIR"

# curl >= 7.71 retries HTTP errors (403/429/503) too, not just network
# failures. Older curl rejects the flag outright, so probe before using it.
RETRY_ALL=""
if curl --help all 2>/dev/null | grep -q -- --retry-all-errors; then
  RETRY_ALL="--retry-all-errors"
fi

# RAW: E20260218-08.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-08.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-08.raw

# OTHER: sample_description.txt
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/sample_description.txt https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/sample_description.txt

# SEARCH: Shoot_combined_modified_peptide.tsv
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/Shoot_combined_modified_peptide.tsv https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_modified_peptide.tsv

# SEARCH: Root_combined_modified_peptide.tsv
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/Root_combined_modified_peptide.tsv https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_modified_peptide.tsv

# SEARCH: Shoot_combined_protein.tsv
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/Shoot_combined_protein.tsv https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_protein.tsv

# RAW: E20260218-16.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-16.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-16.raw

# RAW: E20260218-12.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-12.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-12.raw

# RAW: E20260218-07.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-07.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-07.raw

# RAW: E20260218-04.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-04.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-04.raw

# SEARCH: Root_WT_vs_VZ312_combined_protein.csv
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/Root_WT_vs_VZ312_combined_protein.csv https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_WT_vs_VZ312_combined_protein.csv

# RAW: E20260218-11.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-11.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-11.raw

# RAW: E20260218-03.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-03.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-03.raw

# RAW: E20260218-02.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-02.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-02.raw

# RAW: E20260218-06.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-06.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-06.raw

# RAW: E20260218-10.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-10.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-10.raw

# SEARCH: Root_combined_protein.tsv
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/Root_combined_protein.tsv https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_protein.tsv

# RAW: E20260218-14.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-14.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-14.raw

# RAW: E20260218-15.raw
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/E20260218-15.raw https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-15.raw

echo "Downloaded 18 file(s) to $OUTDIR"
""", 'report.md': """\
# PRIDE report — PXD084218

Source: [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/projects/PXD084218)

## metadata

```
accession   : PXD084218
title       : Genetic code expansion enables plant-directed control of bacterial activity
submitted   : 2026-09-15   published: 2026-09-16
type        : PARTIAL
organisms   : Arabidopsis thaliana (mouse-ear cress)
instruments : Orbitrap Eclipse
diseases    : \n\
keywords    : Bioengineering, Non-canonical amino acids, Bacillus subtilis, Plant-microbe interactions, O-methyl-l-tyrosine (omy)

description : Microbes have the potential to provide growth-promoting benefits to plants and crops, but effective plant–microbe partnerships require mechanisms to control microbial activity in a programmable and host-specific manner. Using O-methyl-L-tyrosine (OMY), a non-canonical amino acid (ncAA), as a model compound, we established a plant-based approach for controlling the activity of associated microbes. Arabidopsis thaliana plants engineered to produce OMY (proUBQ10::sfGFP-MfnG) exhibited reduced primary root length, decreased lateral root density, and impaired gravitropism compared with wild-type Col-0 plants. To investigate the basis of OMY-associated growth inhibition, we tested whether OMY phytotoxicity results from its misincorporation into endogenous proteins, a mechanism reported for other ncAAs. Specifically, we examined whether OMY-producing Arabidopsis plants accumulate O-methylated L-tyrosine residues in native proteins. Proteomic analysis was performed using shoot and root tissues from OMY-producing proUBQ10::sfGFP-MfnG plants and wild-type controls. The resulting proteomic data were analyzed to identify peptides containing O-methylated L-tyrosine residues and to determine whether their abundance was increased in OMY-producing plants relative to wild-type controls.
```

## files

```
18 file(s) for PXD084218:

E20260218-08.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-08.raw
sample_description.txt	[OTHER]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/sample_description.txt
Shoot_combined_modified_peptide.tsv	[SEARCH]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_modified_peptide.tsv
Root_combined_modified_peptide.tsv	[SEARCH]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_modified_peptide.tsv
Shoot_combined_protein.tsv	[SEARCH]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_protein.tsv
E20260218-16.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-16.raw
E20260218-12.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-12.raw
E20260218-07.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-07.raw
E20260218-04.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-04.raw
Root_WT_vs_VZ312_combined_protein.csv	[SEARCH]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_WT_vs_VZ312_combined_protein.csv
E20260218-11.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-11.raw
E20260218-03.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-03.raw
E20260218-02.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-02.raw
E20260218-06.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-06.raw
E20260218-10.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-10.raw
Root_combined_protein.tsv	[SEARCH]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_protein.tsv
E20260218-14.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-14.raw
E20260218-15.raw	[RAW]	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-15.raw
```

## metadata-table

```
Wrote 1 row(s) x 10 column(s) for PXD084218 from project metadata (no SDRF) to <output>/tables/metadata.tsv
```

## samplesheet

```
No submitter SDRF in PRIDE; generating a minimal SDRF from data files + metadata.
Generated minimal SDRF for 12 run(s) -> <output>/PXD084218.sdrf.tsv
Minimal SDRF columns present: 19/19.
Review placeholder values (acquisition method, disease, instrument, tolerances, enzyme, modifications, organism part, factor value) before running quantms.
```

## download-script

```
Wrote <output>/download_pride.sh: 18 download command(s) using curl.
Submit with: sbatch <output>/download_pride.sh   (or run directly: ./<output>/download_pride.sh)
```

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "status": "ok",
  "ok": true,
  "skill": "pride-fetch",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "accession": "PXD084218",
    "commands": [
      "metadata",
      "files",
      "metadata-table",
      "samplesheet",
      "download-script"
    ],
    "demo": true
  },
  "data": {
    "sections": {
      "metadata": "accession   : PXD084218\\ntitle       : Genetic code expansion enables plant-directed control of bacterial activity\\nsubmitted   : 2026-09-15   published: 2026-09-16\\ntype        : PARTIAL\\norganisms   : Arabidopsis thaliana (mouse-ear cress)\\ninstruments : Orbitrap Eclipse\\ndiseases    : \\nkeywords    : Bioengineering, Non-canonical amino acids, Bacillus subtilis, Plant-microbe interactions, O-methyl-l-tyrosine (omy)\\n\\ndescription : Microbes have the potential to provide growth-promoting benefits to plants and crops, but effective plant\\u2013microbe partnerships require mechanisms to control microbial activity in a programmable and host-specific manner. Using O-methyl-L-tyrosine (OMY), a non-canonical amino acid (ncAA), as a model compound, we established a plant-based approach for controlling the activity of associated microbes. Arabidopsis thaliana plants engineered to produce OMY (proUBQ10::sfGFP-MfnG) exhibited reduced primary root length, decreased lateral root density, and impaired gravitropism compared with wild-type Col-0 plants. To investigate the basis of OMY-associated growth inhibition, we tested whether OMY phytotoxicity results from its misincorporation into endogenous proteins, a mechanism reported for other ncAAs. Specifically, we examined whether OMY-producing Arabidopsis plants accumulate O-methylated L-tyrosine residues in native proteins. Proteomic analysis was performed using shoot and root tissues from OMY-producing proUBQ10::sfGFP-MfnG plants and wild-type controls. The resulting proteomic data were analyzed to identify peptides containing O-methylated L-tyrosine residues and to determine whether their abundance was increased in OMY-producing plants relative to wild-type controls.\\n",
      "files": "18 file(s) for PXD084218:\\n\\nE20260218-08.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-08.raw\\nsample_description.txt\\t[OTHER]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/sample_description.txt\\nShoot_combined_modified_peptide.tsv\\t[SEARCH]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_modified_peptide.tsv\\nRoot_combined_modified_peptide.tsv\\t[SEARCH]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_modified_peptide.tsv\\nShoot_combined_protein.tsv\\t[SEARCH]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Shoot_combined_protein.tsv\\nE20260218-16.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-16.raw\\nE20260218-12.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-12.raw\\nE20260218-07.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-07.raw\\nE20260218-04.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-04.raw\\nRoot_WT_vs_VZ312_combined_protein.csv\\t[SEARCH]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_WT_vs_VZ312_combined_protein.csv\\nE20260218-11.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-11.raw\\nE20260218-03.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-03.raw\\nE20260218-02.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-02.raw\\nE20260218-06.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-06.raw\\nE20260218-10.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-10.raw\\nRoot_combined_protein.tsv\\t[SEARCH]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/Root_combined_protein.tsv\\nE20260218-14.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-14.raw\\nE20260218-15.raw\\t[RAW]\\tftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD084218/E20260218-15.raw\\n",
      "metadata-table": "Wrote 1 row(s) x 10 column(s) for PXD084218 from project metadata (no SDRF) to <output>/tables/metadata.tsv\\n",
      "samplesheet": "No submitter SDRF in PRIDE; generating a minimal SDRF from data files + metadata.\\nGenerated minimal SDRF for 12 run(s) -> <output>/PXD084218.sdrf.tsv\\nMinimal SDRF columns present: 19/19.\\nReview placeholder values (acquisition method, disease, instrument, tolerances, enzyme, modifications, organism part, factor value) before running quantms.\\n",
      "download-script": "Wrote <output>/download_pride.sh: 18 download command(s) using curl.\\nSubmit with: sbatch <output>/download_pride.sh   (or run directly: ./<output>/download_pride.sh)\\n"
    }
  }
}\
""", 'tables/metadata.tsv': """\
sample	replicate	species	sex	age	condition	genotype	treatment	tissue	instrument
PXD084218	1	Arabidopsis thaliana	NA	NA	NA	NA	NA	NA	Orbitrap Eclipse
"""})
