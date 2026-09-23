"""Tests for biostudies-fetch.

Run with: pytest skills/biostudies-fetch/tests/test_biostudies_fetch.py -v

No network required: every HTTP call is mocked against the committed demo
fixtures in ../examples/.
"""

import json
import re
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "biostudies_fetch.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO_STUDY = EXAMPLES / "demo_S-BSST2074.json"

sys.path.insert(0, str(SKILL_DIR))


def _load(name):
    return json.loads((EXAMPLES / name).read_text())


class _FakeResponse:
    """Minimal stand-in for the object urlopen returns as a context manager."""

    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self, n=None):
        if n is None:
            payload, self._payload = self._payload, b""
            return payload
        chunk, self._payload = self._payload[:n], self._payload[n:]
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake_urlopen(payload: dict | bytes):
    """Return a urlopen replacement that records the URLs it was asked for."""
    body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
    calls = []

    def _open(req, *args, **kwargs):
        calls.append(getattr(req, "full_url", req))
        return _FakeResponse(body)

    _open.calls = calls
    return _open


# --------------------------------------------------------------------------
# Vendored API module: the upstream logic, unchanged behaviour
# --------------------------------------------------------------------------


class TestVendoredApi:
    def test_walk_files_finds_every_file_node(self):
        import biostudies_fetch_api as api

        study = _load("demo_S-BSST2074.json")
        files = list(api.walk_files(study.get("section", {})))
        assert [f["path"] for f in files] == ["GRCm38_masked_allStrains.zip"]

    def test_clean_val_replaces_control_characters(self):
        """'Clean output fields': control chars, including CR/LF/tab, become _."""
        import biostudies_fetch_api as api

        assert api._clean_val("a\tb\nc") == "a_b_c"
        assert api._clean_val("  spaced   out  ") == "spaced out"

    def test_clean_val_maps_null_markers_to_empty(self):
        import biostudies_fetch_api as api

        for null in ("", "NA", "n/a", "not applicable", "unknown", "--"):
            assert api._clean_val(null) == ""

    def test_harmonize_row_always_has_the_core_columns(self):
        import biostudies_fetch_api as api

        row = api.harmonize_row("S1", "1", {"Organism": "Mus caroli"})
        for col in api.METADATA_COLS:
            assert col in row
        assert row["species"] == "Mus caroli"
        assert row["sex"] == "NA"

    def test_harmonize_row_promotes_unknown_characteristics(self):
        """Nothing is dropped: an unmapped characteristic gets its own column."""
        import biostudies_fetch_api as api

        row = api.harmonize_row("S1", "1", {"Sequencing platform": "Illumina"})
        assert row["sequencing_platform"] == "Illumina"

    def test_metadata_tsv_is_lf_terminated_not_crlf(self, tmp_path):
        """csv.writer defaults to CRLF; the skill must write LF."""
        import biostudies_fetch_api as api

        out = tmp_path / "metadata.tsv"
        api.write_metadata_tsv([api.harmonize_row("S1", "1", {})], out)
        assert b"\r\n" not in out.read_bytes()


    def test_ebi_biosample_attrs_ignores_non_biosample_ids(self):
        """Only SAME*/SAMN*/SAMD* are BioSample ids; anything else short-circuits."""
        import biostudies_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            assert api.ebi_biosample_attrs("S-BSST2074") == {}
            mocked.assert_not_called()


# --------------------------------------------------------------------------
# Entry point: the ClawBio CLI contract
# --------------------------------------------------------------------------


class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
        assert result.returncode != 0

    def test_command_flag_dispatches_to_metadata(self, tmp_path):
        """--command is the runner-reachable form; positionals are dropped by cli.py."""
        import biostudies_fetch as app
        import biostudies_fetch_api as api

        fake = _fake_urlopen(_load("demo_S-BSST2074.json"))
        with patch.object(api.urllib.request, "urlopen", fake):
            app.main(["--command", "metadata", "--accession", "S-BSST2074",
                      "--output", str(tmp_path)])
        assert "S-BSST2074" in fake.calls[0]
        assert (tmp_path / "report.md").exists()

    def test_upstream_positional_form_still_works(self, tmp_path):
        """Anyone following upstream docs must not be stranded."""
        import biostudies_fetch as app
        import biostudies_fetch_api as api

        fake = _fake_urlopen(_load("demo_S-BSST2074.json"))
        with patch.object(api.urllib.request, "urlopen", fake):
            app.main(["metadata", "S-BSST2074", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()

    def test_unknown_command_is_rejected(self, tmp_path):
        import biostudies_fetch as app

        with pytest.raises(SystemExit):
            app.main(["--command", "not-a-command", "--accession", "S-BSST2074",
                      "--output", str(tmp_path)])


class TestDemo:
    def test_demo_makes_no_network_call(self, tmp_path):
        """--demo must be fully offline: nightly_demo_sweep runs it with egress blocked."""
        import biostudies_fetch as app
        import biostudies_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            app.main(["--demo", "--output", str(tmp_path)])
            mocked.assert_not_called()

    def test_demo_writes_report_and_result(self, tmp_path):
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["skill"] == "biostudies-fetch"

    def test_demo_report_carries_the_disclaimer(self, tmp_path):
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in text
        assert "not a medical device" in text

    def test_demo_writes_the_reproducibility_bundle(self, tmp_path):
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        repro = tmp_path / "reproducibility"
        assert (repro / "commands.sh").exists()
        assert (repro / "environment.yml").exists()
        assert (repro / "checksums.sha256").exists()

    def test_checksum_labels_resolve_from_the_output_dir(self, tmp_path):
        """anchor=output_dir, so `cd <out> && sha256sum -c` works."""
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        lines = (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines()
        assert lines
        for line in lines:
            label = line.split("  ", 1)[1]
            assert (tmp_path / label).exists(), f"unresolvable checksum label: {label}"

    def test_demo_writes_the_metadata_table(self, tmp_path):
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        tsv = tmp_path / "tables" / "metadata.tsv"
        assert tsv.exists()
        header = tsv.read_text().splitlines()[0].split("\t")
        assert header[:2] == ["sample", "replicate"]


class TestSafety:
    def test_warns_before_overwriting(self, tmp_path, capsys):
        """AGENTS.md Safety Boundary 5."""
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwritten" in capsys.readouterr().err

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        """Upstream defaulted outputs to cwd; nothing may land in the repo root."""
        import biostudies_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        out = tmp_path / "out"
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(out)])
        assert list(cwd.iterdir()) == []


# --------------------------------------------------------------------------
# Output contract — copied verbatim from scaffold_skill.py
# --------------------------------------------------------------------------


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree.

    Returns output-relative file paths. Directory lines, and any entry whose
    inline comment contains 'optional', are skipped. Returns [] when there is no
    parseable section, so skills without the section are simply not gated.
    """
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files = []
    parents = {}
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
            continue  # the root output_directory/ line
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
    """Every artifact promised in SKILL.md '## Output Structure' must be produced.

    Guards against doc/code drift: a skill documenting an output it never writes.
    Mark conditional artifacts with '(optional)' in the tree comment to exempt them.
    """

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing) + ". Write them, mark them "
            "'(optional)' in the SKILL.md tree, or remove them from the "
            "documented Output Structure."
        )


class TestDeterminism:
    """The report must describe what was produced, not where it was put."""

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import biostudies_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert str(tmp_path) not in (tmp_path / "report.md").read_text()

    def test_two_output_dirs_give_identical_reports(self, tmp_path):
        import biostudies_fetch as app

        a, b = tmp_path / "a", tmp_path / "b"
        app.main(["--demo", "--output", str(a)])
        app.main(["--demo", "--output", str(b)])
        assert (a / "report.md").read_text() == (b / "report.md").read_text()


@pytest.fixture(autouse=True)
def _clear_info_cache():
    """study_info() memoises per accession; keep that out of other tests."""
    import biostudies_fetch_api as api

    api._INFO_CACHE.clear()
    yield
    api._INFO_CACHE.clear()


class TestDownloadBase:
    """The download base is resolved from /studies/{acc}/info, not hardcoded.

    Upstream's /biostudies/files/{acc}/{path} is not broken -- it 302s to the
    same place. We resolve directly because `httpLink` names a different tree
    per collection (fire/ for E-MTAB, pub/databases/ for S-BSST), so no single
    constant is correct, and because the redirect is ~40x slower and times out
    on multi-gigabyte files. Re-verified live 2026-09-22.
    """

    def test_file_url_is_resolved_from_the_info_endpoint(self):
        import biostudies_fetch_api as api

        info = {"httpLink": "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/074/S-BSST2074",
                "relPath": "S-BSST/074/S-BSST2074"}
        with patch.object(api, "get_json", return_value=info):
            url = api.file_url("S-BSST2074", "GRCm38_masked_allStrains.zip")
        assert url == (
            "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/074/S-BSST2074"
            "/Files/GRCm38_masked_allStrains.zip")

    def test_file_url_percent_encodes_under_the_resolved_base(self):
        import biostudies_fetch_api as api

        info = {"httpLink": "https://ftp.ebi.ac.uk/biostudies/fire/S-B/001/S-B1"}
        with patch.object(api, "get_json", return_value=info):
            assert api.file_url("S-B1", "a dir/b.txt").endswith("/Files/a%20dir/b.txt")

    def test_falls_back_to_the_legacy_path_when_info_has_no_link(self):
        """Never crash on an info payload that omits httpLink; degrade instead."""
        import biostudies_fetch_api as api

        with patch.object(api, "get_json", return_value={}):
            assert api.file_url("S-B1", "b.txt") == \
                "https://www.ebi.ac.uk/biostudies/files/S-B1/b.txt"
