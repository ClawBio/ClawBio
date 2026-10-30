"""Tests for clawbio.privacy: sealed results and the genomic-data policy."""

import gzip
import json
import struct
import zipfile
from pathlib import Path

import pytest

from hypothesis import example, given
from hypothesis import strategies as st

from clawbio import privacy
from clawbio.privacy import SEALED_FIELDS, seal_result

# Uppercase bases and digits never occur in the sealed keys, the skill name or
# the output path below, so finding one of these strings in the sealed output
# can only mean it leaked from a field that should have been dropped.
CANARY = st.text(alphabet="ACGT0123456789", min_size=8, max_size=40)


def _result(**overrides):
    base = {
        "skill": "pharmgx",
        "success": True,
        "exit_code": 0,
        "output_dir": "/tmp/out",
        "files": ["report.md", "result.json"],
        "stdout": "",
        "stderr": "",
        "duration_seconds": 1.5,
    }
    base.update(overrides)
    return base


@given(canary=CANARY, field=st.sampled_from(["stdout", "stderr", "report_md", "chat_summary", "x_new_field"]))
@example(canary="CYP2D6000", field="stdout")
def test_sealed_result_never_carries_content_fields(canary, field):
    sealed = seal_result(_result(**{field: f"rs3892097 {canary} *4/*4 Poor Metabolizer"}))
    assert canary not in json.dumps(sealed)


def test_sealed_result_keeps_only_delivery_metadata():
    sealed = seal_result(_result(stdout="genotype data", report_md="# Report"))
    assert set(sealed) <= set(SEALED_FIELDS) | {"sealed", "note"}
    assert sealed["sealed"] is True
    assert sealed["skill"] == "pharmgx"
    assert sealed["files"] == ["report.md", "result.json"]
    assert sealed["output_dir"] == "/tmp/out"


def test_sealed_result_fails_closed_on_missing_success():
    sealed = seal_result({"skill": "pharmgx", "stdout": "data"})
    assert sealed["success"] is False


@pytest.mark.parametrize("success", [True, False, None, 0, 1, 1.0, "True", "yes", ""])
def test_success_is_true_only_for_exactly_true(success):
    assert seal_result(_result(success=success))["success"] is (success is True)


@given(canary=CANARY)
def test_sealing_is_idempotent(canary):
    once = seal_result(_result(stdout=canary))
    assert seal_result(once) == once


def test_cli_sealed_run_prints_only_delivery_metadata(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    clawbio_py = Path(__file__).resolve().parents[2] / "clawbio.py"
    out = tmp_path / "out"
    proc = subprocess.run(
        [sys.executable, str(clawbio_py), "run", "pharmgx", "--demo", "--sealed", "--output", str(out)],
        capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    sealed = json.loads(proc.stdout)
    assert sealed["sealed"] is True and sealed["success"] is True
    assert "report.md" in sealed["files"]
    # The report is on disk for the user, but none of it is on stdout for the model.
    report = (out / "report.md").read_text()
    assert "CYP2D6" in report
    assert "CYP2D6" not in proc.stdout


def test_unknown_skill_gets_a_fixed_hint_not_the_cli_error(monkeypatch, capsys):
    # Live test: the model guessed `bio-orchestrator`, got success=false with no reason, and gave up.
    from clawbio import cli

    monkeypatch.setattr("sys.argv", ["clawbio.py", "run", "bio-orchestrator", "--demo", "--sealed"])
    with pytest.raises(SystemExit):
        cli.main()
    sealed = json.loads(capsys.readouterr().out)
    assert "clawbio.py list" in sealed["error"]
    assert "Available" not in json.dumps(sealed)


def test_failed_run_never_passes_stderr_through():
    sealed = seal_result(_result(success=False, exit_code=-1, stderr="Upload failed: rs3892097 CYP2D6000"))
    assert "error" not in sealed
    assert "CYP2D6000" not in json.dumps(sealed)


def test_cli_sealed_summary_skill_writes_a_report(tmp_path, monkeypatch, capsys):
    # `compare` prints to stdout unless given an output folder; sealed drops stdout, so it delivered nothing.
    from clawbio import cli

    monkeypatch.setattr(cli, "DEFAULT_OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr("sys.argv", ["clawbio.py", "run", "compare", "--demo", "--sealed"])
    with pytest.raises(SystemExit) as exit_info:
        cli.main()
    sealed = json.loads(capsys.readouterr().out)
    assert exit_info.value.code == 0
    assert sealed["output_dir"].startswith(str(tmp_path))
    assert "report.md" in sealed["files"]


# --- The genomic-data policy ------------------------------------------------- #

GENOME_CANARY = "CANARY7F3A"
VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n"
    f"chr22\t42127941\trs3892097\tG\tA\t50\tPASS\tMARK={GENOME_CANARY}\tGT\t0/1\n"
)
TWENTY_THREE = (
    "# rsid\tchromosome\tposition\tgenotype\n"
    "rs4477212\t1\t82154\tAA\n"
    "rs3094315\t1\t752566\tAG\n"
    f"rs3131972\t1\t752721\t{GENOME_CANARY}\n"
)
_SEALED = "python clawbio.py run pharmgx --input /tmp/s.vcf --sealed"

CASES = json.loads(Path(__file__).with_name("privacy_cases.json").read_text())
CHECKS = {
    "genomic_text": privacy.is_genomic_text,
    "genomic_path": privacy.is_genomic_path,
    "sealed_run": privacy.is_sealed_run,
}


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"{c['check']}:{c['input'][:40]!r}")
def test_detection_case(case):
    assert CHECKS[case["check"]](case["input"]) is case["expected"]


def test_every_check_has_positive_and_negative_cases():
    for check in CHECKS:
        assert {c["expected"] for c in CASES if c["check"] == check} == {True, False}, check


@given(
    stem=st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-", min_size=1, max_size=12),
    ext=st.sampled_from(privacy.GENOMIC_EXTENSIONS),
    upper=st.booleans(),
)
def test_any_genomic_extension_marks_a_path_genomic(stem, ext, upper):
    name = f"/uploads/{stem}{ext}"
    assert privacy.is_genomic_path(name.upper() if upper else name)


def test_reading_a_genomic_file_is_refused_with_the_sealed_command(tmp_path):
    message = privacy.check_path(str(tmp_path / "sample.vcf"), [])
    assert message.startswith(privacy.WITHHELD)
    assert "--sealed" in message


def test_reading_an_ordinary_file_is_allowed(tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text("meeting notes\n")
    assert privacy.check_path(str(notes), []) is None


def test_reading_clawbio_output_is_refused(tmp_path):
    assert privacy.check_path(str(tmp_path / "output" / "pharmgx_20261030_111316" / "report.md"), [str(tmp_path)])


def test_the_sealed_run_is_the_only_allowed_command_on_a_genomic_file():
    assert privacy.check_command(_SEALED, []) is None
    assert privacy.check_command("cat /tmp/s.vcf", [])
    assert privacy.check_command("python clawbio.py run pharmgx --input /tmp/s.vcf", [])
    assert privacy.check_command("ls -la /tmp", []) is None


@given(caption=st.text())
def test_scrubbed_text_never_carries_genotypes(caption):
    assert GENOME_CANARY not in privacy.scrub_text(f"{caption}\n{VCF}", []) or GENOME_CANARY in caption


def test_a_saved_upload_is_swapped_for_its_path_and_the_sealed_command(tmp_path):
    path = tmp_path / "sample.vcf"
    path.write_text(VCF)
    text = privacy.replace_uploads(f"[Content of sample.vcf]:\n{VCF}\n\nhi", [str(path)], [])
    assert GENOME_CANARY not in text
    assert str(path) in text and "--sealed" in text and text.endswith("hi")


@pytest.mark.parametrize("command", [
    "python clawbio.py run pharmgx --input /tmp/x--sealed.vcf",
    "python clawbio.py run pharmgx --demo --output /tmp/o--sealed",
    "python clawbio.py run pharmgx --trait --sealed-ish --input a.vcf",
])
def test_sealed_inside_a_path_or_value_is_not_the_flag(command):
    # Such a run is unsealed and prints the report for the model.
    assert not privacy.is_sealed_run(command)


def test_quoted_input_path_with_spaces_is_still_a_sealed_run():
    assert privacy.is_sealed_run("python clawbio.py run pharmgx --input '/tmp/my sample.vcf' --sealed")


@pytest.mark.parametrize("sep", [";", "&&", "||", "|", "&", "\n", ">", "<", "$(", "`"])
def test_nothing_can_be_appended_to_a_sealed_run(sep):
    assert not privacy.is_sealed_run(f"{_SEALED} {sep}cat /tmp/s.vcf")


@pytest.mark.parametrize("flag", ["--o", "--ou", "--out", "--outp", "--outpu", "--output"])
@pytest.mark.parametrize("value", [" /tmp/x", "=/tmp/x"])
def test_no_prefix_of_output_rides_along_with_a_sealed_run(flag, value):
    # argparse accepts any unambiguous prefix (--out /tmp/x); the report would land outside output/.
    assert not privacy.is_sealed_run(f"{_SEALED} {flag}{value}")


@pytest.mark.parametrize("body", [
    VCF,
    TWENTY_THREE,
    "rsid\tchromosome\tposition\tresult\n",
    "##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
])
def test_a_leading_byte_order_mark_never_hides_a_genomic_header(body):
    # Windows and Excel exports start with a byte-order mark, which hid the header from `^`.
    assert privacy.is_genomic_text("\ufeff" + body)


def _bcf_bytes():
    # BCF: magic, header length, then the VCF header text, all BGZF (gzip) compressed.
    header = VCF.encode() + b"\x00"
    return gzip.compress(b"BCF\x02\x02" + struct.pack("<I", len(header)) + header)


def test_compressed_genomic_files_are_found_by_content(tmp_path):
    gz = tmp_path / "genome.txt.gz"
    gz.write_bytes(gzip.compress(TWENTY_THREE.encode()))
    zipped = tmp_path / "genome.zip"
    with zipfile.ZipFile(zipped, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("README.txt", "23andMe raw data export\n")
        z.writestr("genome_Test_User_v5_Full.txt", TWENTY_THREE)
    bcf = tmp_path / "sample.bcf"
    bcf.write_bytes(_bcf_bytes())
    for path in (gz, zipped, bcf):
        assert privacy.is_genomic_path(str(path)), path.name


def test_broken_compressed_files_do_not_crash_the_policy(tmp_path):
    for name, data in [("x.gz", b"\x1f\x8b\x08junk"), ("x.zip", b"PK\x03\x04junk"), ("empty.zip", b"")]:
        (tmp_path / name).write_bytes(data)
        assert privacy.is_genomic_path(str(tmp_path / name)) in (True, False)


@pytest.mark.parametrize("position", [0, 1, 2, 3])
@pytest.mark.parametrize("method", [zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED])
def test_a_genomic_member_is_found_wherever_it_sits_in_a_zip(tmp_path, position, method):
    members = [(f"notes_{i}.txt", "meeting notes\n" * 50) for i in range(3)]
    members.insert(position, ("genome_Test_User_v5_Full.txt", TWENTY_THREE))
    path = tmp_path / "export.zip"
    with zipfile.ZipFile(path, "w", method) as z:
        for name, body in members:
            z.writestr(name, body)
    assert privacy.is_genomic_path(str(path))
