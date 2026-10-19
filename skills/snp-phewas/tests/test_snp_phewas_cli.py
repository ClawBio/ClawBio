"""The entry script: --demo end to end (offline), the output contract SKILL.md promises,
the reproducibility bundle, --quick, --list-datasets, --input, and that a run writes
nothing outside --output."""

from __future__ import annotations

import hashlib
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

import snp_phewas

SKILL_DIR = Path(__file__).resolve().parent.parent
ROOT = SKILL_DIR.parents[1]


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("snp_phewas_demo")
    assert snp_phewas.main(["--demo", "--output", str(out)]) == 0
    return out


def _parse_output_contract(skill_md: Path) -> list[str]:
    """Files the `## Output Structure` tree promises (lines marked optional excluded)."""
    text = skill_md.read_text(encoding="utf-8")
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files, parents = [], {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        prefix, name = re.match(r"^([\s│├└─]*)(.*)$", entry).groups()
        name = name.strip()
        depth = len(prefix) // 4
        if not name or depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(f"{rel}/{name}" if rel else name)
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced(self, demo_out):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert "reproducibility/checksums.sha256" in promised
        missing = [p for p in promised if not (demo_out / p).exists()]
        assert not missing, "SKILL.md promises artifacts the demo did not produce: " + ", ".join(missing)


def test_demo_report_states_scope_sources_and_disclaimer(demo_out):
    report = (demo_out / "report.md").read_text()
    assert "**Demo data.**" in report
    assert "We want to acknowledge the participants and investigators of the FinnGen study" in report
    assert "using the UK Biobank Resource" in report
    assert "not a medical device" in report
    assert "replicates (significant in both, same direction): 2" in report
    assert "1.58e-465" in report  # LDL p from -log10 p 464.8, not an underflowed 0
    import json

    res = json.loads((demo_out / "result.json").read_text())
    assert res["skill"] == "snp-phewas" and res["summary"]["demo"] is True
    assert res["summary"]["replication"]["verdicts"]["replicates"] == 2
    assert res["data"]["uncataloged"] == ["K11_GINGIVITIS_PERIODONTAL"]


def test_reproducibility_bundle_verifies_from_the_output_dir(demo_out):
    repro = demo_out / "reproducibility"
    assert "python skills/snp-phewas/snp_phewas.py --demo --output" in (repro / "commands.sh").read_text()
    assert "matplotlib" in (repro / "environment.yml").read_text()
    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {ln.split("  ", 1)[1] for ln in lines}
    assert {"report.md", "result.json", "tables/snp_phewas_rows.tsv",
            "reproducibility/commands.sh"} <= labels
    for ln in lines:
        digest, label = ln.split("  ", 1)
        assert hashlib.sha256((demo_out / label).read_bytes()).hexdigest() == digest, label


def test_quick_demo_reads_finngen_only_and_says_so(tmp_path):
    assert snp_phewas.main(["--demo", "--quick", "--output", str(tmp_path)]) == 0
    report = (tmp_path / "report.md").read_text()
    assert "Quick mode: FinnGen web API only" in report
    assert not (tmp_path / "tables" / "replication.tsv").exists()
    assert "pan_ukbb_eur" not in (tmp_path / "tables" / "snp_phewas_rows.tsv").read_text()


def test_bad_flag_combinations_are_refused(tmp_path):
    out = ["--output", str(tmp_path)]
    for argv in (["--quick", "--all-outcomes", "--variant", "rs1", *out],
                 ["--quick", "--dataset", "pan_ukbb_eur", "--variant", "rs1", *out],
                 ["--variant", "rs1"],
                 out):
        with pytest.raises(SystemExit):
            snp_phewas.main(argv)
    assert not any(tmp_path.iterdir())


def test_list_datasets(capsys):
    assert snp_phewas.main(["--list-datasets"]) == 0
    out = capsys.readouterr().out
    assert "pan_ukbb_eur (default)" in out and "7,160 outcomes" in out
    assert "finngen_r12 (default)" in out and "2,469 outcomes" in out
    assert "is not included" in out


def test_input_file_holds_exactly_one_variant(tmp_path):
    one = tmp_path / "one.txt"
    one.write_text("# my variant\nrs12740374\n")
    assert snp_phewas._variant_from_input(str(one)) == "rs12740374"
    two = tmp_path / "two.txt"
    two.write_text("rs1\nrs2\n")
    with pytest.raises(SystemExit, match="one variant per run"):
        snp_phewas._variant_from_input(str(two))


@given(st.floats(0, 2000, allow_nan=False))
def test_format_p_never_underflows(nlp):
    s = snp_phewas.format_p(nlp)
    mant, exp = s.split("e")
    assert float(mant) > 0
    assert math.isclose(-math.log10(float(mant)) - int(exp), nlp, abs_tol=0.01)


def test_demo_through_the_runner_writes_only_under_output(tmp_path):
    """`python clawbio.py run snp-phewas --demo` from the registered alias, and no file
    appears in the skill folder (the nightly sweep's rule)."""
    def snapshot():
        return {p for p in SKILL_DIR.rglob("*")
                if "__pycache__" not in p.parts and ".hypothesis" not in p.parts}

    before = snapshot()
    out = tmp_path / "run"
    proc = subprocess.run([sys.executable, "clawbio.py", "run", "snp-phewas", "--demo",
                           "--output", str(out)], cwd=ROOT, capture_output=True, text=True,
                          timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert (out / "report.md").exists() and (out / "figures" / "manhattan.png").exists()
    assert snapshot() == before


def test_a_live_run_warns_about_the_index_download_before_it_starts(monkeypatch, capsys, tmp_path):
    """The first default run downloads 2.08 GB of index files; the user hears it first."""
    import snp_phewas as cli

    def stop(*a, **k):
        raise SystemExit(0)

    monkeypatch.setattr(cli, "lookup_variant", stop)
    with pytest.raises(SystemExit):
        cli.main(["--variant", "rs12740374", "--output", str(tmp_path)])
    assert "2.08 GB" in capsys.readouterr().err
