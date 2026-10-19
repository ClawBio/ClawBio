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
    assert "Million Veteran Program (Verma A et al., Science 2024" in report
    assert "not a medical device" in report
    assert "1.58e-465" in report  # LDL p from -log10 p 464.8, not an underflowed 0
    assert "6.02e-587" in report  # MVP LDL: p = 0 in the file, derived from beta / SE
    # outcomes looked up per biobank: the threshold is scoped, never called phenome-wide
    summary = report.split("## Summary", 1)[1].split("## Notes", 1)[0]
    assert summary.count("(over the 3 outcomes looked up)") == 2
    assert summary.count("(over the 2 outcomes looked up)") == 1
    assert "phenome-wide" not in summary
    # an MVP study with no readable file is never counted as read
    assert ("tested in 2 of 2 outcomes read (2,070 in the catalog, 849 of them with no "
            "readable file)") in summary
    import json

    res = json.loads((demo_out / "result.json").read_text())
    assert res["skill"] == "snp-phewas" and res["summary"]["demo"] is True
    pairs = {tuple(p["pair"]): p for p in res["summary"]["replication"]["pairs"]}
    assert {k: p["verdicts"]["replicates"] for k, p in pairs.items()} == {
        ("pan_ukbb_eur", "finngen_r12"): 2, ("pan_ukbb_eur", "mvp_eur"): 2,
        ("finngen_r12", "mvp_eur"): 1}
    assert res["summary"]["replication"]["three_way"]["classes"]["replicates_in_all_three"] == 1
    assert res["summary"]["datasets"]["mvp_eur"]["n_read"] == 2
    assert res["data"]["uncataloged"] == ["K11_GINGIVITIS_PERIODONTAL"]


def test_demo_report_names_every_pair_and_all_three(demo_out, capsys):
    """With three biobanks, a report of the first pair alone would hide two comparisons
    the tables and the matrix carry."""
    report = (demo_out / "report.md").read_text()
    for a, b in (("Pan-UK Biobank (European ancestry)", "FinnGen R12 (Finnish)"),
                 ("Pan-UK Biobank (European ancestry)",
                  "Million Veteran Program (European genetic ancestry)"),
                 ("FinnGen R12 (Finnish)",
                  "Million Veteran Program (European genetic ancestry)")):
        assert f"**{a} and {b}**" in report, (a, b)
        assert f"**Replication, {a} and {b}**" in report, (a, b)
    assert "**All three biobanks** (Pan-UKBB EUR, FinnGen and MVP EUR)" in report
    assert "replicates in all three (significant in each, same direction): 1 (Ischemic Heart Disease)" in report
    assert "harmonized across at least two of the datasets" in report
    # the multi-dataset sentences use one list form: no label's own comma reads as a separator
    assert ("**Datasets**: Pan-UK Biobank (European ancestry), FinnGen R12 (Finnish) and "
            "Million Veteran Program (European genetic ancestry)") in report


def test_the_replication_lines_reach_stdout(tmp_path, capsys):
    assert snp_phewas.main(["--demo", "--output", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "**All three biobanks**" in out
    assert "**FinnGen R12 (Finnish) and Million Veteran Program (European genetic ancestry)**" in out


class _Absent:
    """A reader that answers 'absent' for every outcome asked about (no network)."""

    def read(self, outcomes, tgt, catalog_ids=None):
        import snp_phewas_schema as S
        from snp_phewas_sources import RawRecord

        recs = {o.outcome_id: RawRecord(S.ST_ABSENT, route="tabix") for o in outcomes}
        return (recs, {"route": "tabix"}, []) if catalog_ids is not None else (recs, {"route": "tabix"})


def test_a_matched_read_says_so_first_with_its_counts(offline_kwargs, tmp_path):
    """The narrowing is said where the answer is given, per biobank, counted from the rows."""
    import argparse

    from snp_phewas_lookup import lookup_variant

    offline_kwargs["readers"] = {d: _Absent() for d in ("pan_ukbb_eur", "finngen_r12", "mvp_eur")}
    offline_kwargs["outcomes"] = None
    r = lookup_variant("rs12740374", **offline_kwargs)
    args = argparse.Namespace(demo=False, quick=False)
    report = snp_phewas.write_report(r, None, tmp_path, args, 0.0).read_text()
    covers = report.split("## What this answer covers", 1)[1].split("## Summary", 1)[0]
    assert ("**This was the quicker read: Pan-UK Biobank (European ancestry), FinnGen R12 "
            "(Finnish) and Million Veteran Program (European genetic ancestry) were read only "
            "over the outcomes matched across the biobanks, the same rule in each.**") in covers
    assert ("Million Veteran Program (European genetic ancestry): read over the 727 outcomes "
            "with a readable counterpart in FinnGen or Pan-UK Biobank") in covers
    assert "494 were not looked up and 849 have no readable file" in covers
    summary = report.split("## Summary", 1)[1].split("## Notes", 1)[0]
    assert "tested in 0 of 727 outcomes read (2,070 in the catalog, 849 of them" in summary


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
    assert "mvp_eur (default)" in out and "2,070 outcomes" in out and "mvp_meta" in out
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
    """The first default run downloads 4.53 GB of index files (summed from Content-Length for
    the matched files of the three biobanks); the user hears it, and what will be read, first."""
    import snp_phewas as cli

    def stop(*a, **k):
        raise SystemExit(0)

    monkeypatch.setattr(cli, "lookup_variant", stop)
    with pytest.raises(SystemExit):
        cli.main(["--variant", "rs12740374", "--output", str(tmp_path)])
    err = capsys.readouterr().err
    assert "4.53 GB" in err and "(1,156, 622, 727 files in Pan-UK Biobank (European ancestry), " \
        "FinnGen R12 (Finnish) and Million Veteran Program (European genetic ancestry))" in err
    with pytest.raises(SystemExit):
        cli.main(["--variant", "rs12740374", "--all-outcomes", "--output", str(tmp_path)])
    err = capsys.readouterr().err
    assert "20.0 GB" in err and "(7,160, 2,469, 1,221 files" in err


def test_a_low_confidence_value_is_flagged_and_never_counted_as_a_pass(offline_kwargs, tmp_path):
    import argparse

    import snp_phewas_schema as S
    from snp_phewas_lookup import lookup_variant

    r = lookup_variant("rs12740374", **offline_kwargs)
    (ldl,) = [x for x in r.rows if x.outcome_id == "biomarkers-30780-both_sexes-irnt"]
    ldl.status = S.ST_LOW_CONFIDENCE  # p = 1.6e-465, far past any threshold
    args = argparse.Namespace(demo=False, quick=False)
    report = snp_phewas.write_report(r, None, tmp_path, args, 0.0).read_text()
    (line,) = [ln for ln in report.splitlines() if ln.startswith("- **Pan-UK Biobank")]
    assert "2 below the threshold" in line, line
    assert "(1 low-confidence value is not counted)" in line, line
    assert "tested, flagged low confidence by the source" in report  # the table's Status
