"""
Tests for the ClawBio igv-validator skill.
Run with: pytest skills/igv-validator/tests/test_igv_validator.py -v

Counting tests compare against demo/demo_truth.json, written by demo/make_demo_data.py
when it planted each synthetic variant. None of these tests need IGV.
"""

import csv
import json
import time
import re
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "igv_validator.py"
DEMO = SKILL_DIR / "demo"
TRUTH = json.loads((DEMO / "demo_truth.json").read_text())["variants"]
DISCLAIMER_START = "ClawBio is a research and educational tool"

sys.path.insert(0, str(SKILL_DIR))
import igv_validator as iv  # noqa: E402


def run_cli(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo")
    r = run_cli("--demo", "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    return out


@pytest.fixture(scope="module")
def demo_result(demo_out):
    return json.loads((demo_out / "result.json").read_text())


def by_id(result):
    return {v["id"]: v for v in result["data"]["variants"]}


# ── Required ClawBio tests ─────────────────────────────────────────────────────

def test_demo_runs_without_error(demo_out):
    assert (demo_out / "report.md").exists()


def test_output_structure(demo_out):
    for rel in ["report.md", "report.html", "result.json", "tables/support_counts.tsv",
                "reproducibility/commands.sh", "reproducibility/environment.yml",
                "reproducibility/checksums.sha256"]:
        assert (demo_out / rel).exists(), rel


def test_report_contains_disclaimer(demo_out):
    assert DISCLAIMER_START in (demo_out / "report.md").read_text()
    assert DISCLAIMER_START in (demo_out / "report.html").read_text()


def test_rejects_malformed_input(tmp_path):
    bad = tmp_path / "bad.vcf"
    bad.write_text("this is not a VCF\n1 2 3\n")
    r = run_cli("--vcf", bad, "--tumor", DEMO / "demo_tumor.bam", "--normal", DEMO / "demo_normal.bam",
                "--no-igv", "--output", tmp_path / "out")
    assert r.returncode != 0
    assert "Traceback" not in r.stderr
    assert "error" in r.stderr.lower()


def test_empty_input_handled(tmp_path):
    empty = tmp_path / "empty.vcf"
    empty.write_text("")
    r = run_cli("--vcf", empty, "--tumor", DEMO / "demo_tumor.bam", "--normal", DEMO / "demo_normal.bam",
                "--no-igv", "--output", tmp_path / "out")
    assert r.returncode != 0
    assert "Traceback" not in r.stderr
    assert "empty" in r.stderr.lower() or "no variant" in r.stderr.lower()


def test_missing_bam_index_is_a_clean_error(tmp_path):
    bam = tmp_path / "noindex.bam"
    bam.write_bytes((DEMO / "demo_tumor.bam").read_bytes())
    r = run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", bam, "--normal", DEMO / "demo_normal.bam",
                "--no-igv", "--output", tmp_path / "out")
    assert r.returncode != 0 and "Traceback" not in r.stderr
    assert "index" in r.stderr.lower()


# ── Counts against the planted truth ───────────────────────────────────────────

def test_mate_breakends_become_one_variant(demo_result):
    ids = sorted(by_id(demo_result))
    assert ids == ["V1", "V2", "V3", "V4", "V5_1", "V6"]


@pytest.mark.parametrize("vid", ["V1", "V2", "V3", "V6"])
def test_snv_counts_match_truth(demo_result, vid):
    v = by_id(demo_result)[vid]
    for sample, key in [("demo_tumor", "tumor"), ("demo_normal", "normal")]:
        assert v[key]["alt"] == TRUTH[vid][sample]["alt"]
        assert v[key]["depth"] == TRUTH[vid][sample]["depth"]


def test_snv_strands(demo_result):
    v2 = by_id(demo_result)["V2"]["tumor"]
    assert (v2["alt_fwd"], v2["alt_rev"]) == (0, TRUTH["V2"]["demo_tumor"]["alt"])


def test_deletion_counts_clipped_reads_with_reference(demo_result):
    v4 = by_id(demo_result)["V4"]
    assert v4["tumor"]["alt"] == TRUTH["V4"]["demo_tumor"]["alt"]
    assert v4["normal"]["alt"] == 0


def test_deletion_without_reference_is_conservative():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf")
    v4 = next(v for v in variants if v.id == "V4")
    c = iv.count_variant(v4, DEMO / "demo_tumor.bam", fasta=None)
    assert c["alt"] == TRUTH["V4"]["demo_tumor"]["alt_without_clip_rescue"]


def test_translocation_counts(demo_result):
    v5 = by_id(demo_result)["V5_1"]
    t = TRUTH["V5"]["demo_tumor"]
    assert (v5["tumor"]["split"], v5["tumor"]["discordant"], v5["tumor"]["alt"]) == (t["split"], t["discordant"], t["union"])
    assert v5["normal"]["alt"] == 0


def test_overlapping_mates_count_once(tmp_path):
    """Regression: both reads of an overlapping pair carrying the same deletion are one molecule."""
    seq = "".join("ACGT"[(i * 7) % 4] for i in range(400))
    fa = tmp_path / "r.fa"
    fa.write_text(">c1\n" + seq + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 400}]}
    raw = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for flag, start, cig in [(99, 150, "50M3D50M"), (147, 160, "40M3D60M")]:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = "frag1", flag, 0, start
            a.cigarstring, a.mapping_quality = cig, 60
            a.query_sequence = seq[start:start + 50] + seq[start + 53:start + 103] if flag == 99 else \
                seq[start:start + 40] + seq[start + 43:start + 103]
            a.query_qualities = pysam.qualitystring_to_array("?" * 100)
            a.next_reference_id, a.next_reference_start = 0, 160 if flag == 99 else 150
            out.write(a)
    bam = tmp_path / "o.bam"
    pysam.sort("-o", str(bam), str(raw)); pysam.index(str(bam))
    v = iv.Variant(id="d", chrom="c1", pos=200, ref=seq[199:203], alt=seq[199], kind="deletion")
    assert iv.count_variant(v, bam, fasta=fa)["alt"] == 1


# ── Tumor-only mode ────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def tumor_only(tmp_path_factory):
    out = tmp_path_factory.mktemp("tonly")
    r = run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", DEMO / "demo_tumor.bam", "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    return out, json.loads((out / "result.json").read_text())


def test_tumor_only_runs_without_normal(tumor_only):
    out, res = tumor_only
    assert res["summary"]["mode"] == "tumor_only"
    v = by_id(res)
    assert all(x["normal"] is None for x in v.values())
    assert v["V1"]["tumor"]["alt"] == TRUTH["V1"]["demo_tumor"]["alt"]
    assert "strand_bias" in v["V2"]["flags"] and "low_support" in v["V3"]["flags"]


def test_tumor_only_skips_normal_based_checks(tumor_only):
    _, res = tumor_only
    v = by_id(res)
    assert not any(f in x["flags"] for x in v.values() for f in ("normal_support", "germline_site"))
    # without a normal a germline het cannot be told apart from a somatic call
    assert v["V6"]["status"] == "supported"


def test_tumor_only_report_says_so(tumor_only):
    out, _ = tumor_only
    md = (out / "report.md").read_text().lower()
    assert "tumor-only" in md and "germline" in md
    rows = (out / "tables" / "support_counts.tsv").read_text().splitlines()[1:]
    assert rows and all("\tdemo_tumor\t" in r for r in rows)


def test_demo_tumor_only_flag(tmp_path):
    r = run_cli("--demo", "--tumor-only", "--no-igv", "--output", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    assert json.loads((tmp_path / "o" / "result.json").read_text())["summary"]["mode"] == "tumor_only"


def test_igv_commands_tumor_only_load_one_track():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf")
    cmds = iv.igv_commands(variants[0], DEMO / "demo_tumor.bam", None, Path("/tmp/x"))
    assert sum(c.startswith("load ") for c in cmds) == 1


# ── Caller counts vs reads ─────────────────────────────────────────────────────

def test_caller_counts_are_reported(demo_result):
    v = by_id(demo_result)
    assert v["V1"]["caller"] == {"alt": 18, "depth": 64, "vaf_pct": 28.1, "sample": "demo_tumor"}
    assert "caller_disagrees" not in v["V1"]["flags"]
    assert v["V5_1"]["caller"] is None  # breakends carry no AD


def test_caller_disagreement_is_flagged(demo_result):
    v3 = by_id(demo_result)["V3"]
    assert v3["caller"]["alt"] == 12 and v3["tumor"]["alt"] == 2
    assert "caller_disagrees" in v3["flags"]


def test_caller_sample_can_be_chosen(tmp_path):
    out = tmp_path / "o"
    r = run_cli("--demo", "--no-igv", "--caller-sample", "demo_normal", "--output", out)
    assert r.returncode == 0, r.stderr
    v6 = by_id(json.loads((out / "result.json").read_text()))["V6"]
    assert v6["caller"]["alt"] == 25 and v6["caller"]["sample"] == "demo_normal"


def test_tumor_column_guessed_from_name_only_when_unambiguous():
    assert iv.guess_tumor_sample(["normal-ucla", "tumor-ucla-T1"]) == "tumor-ucla-T1"
    assert iv.guess_tumor_sample(["SG_NORMAL", "SG_TUMOR"]) == "SG_TUMOR"
    assert iv.guess_tumor_sample(["tumor_A", "tumor_B"]) is None   # two candidates: ask
    assert iv.guess_tumor_sample(["S1", "S2"]) is None             # no hint: ask


def test_unknown_caller_sample_is_a_clean_error(tmp_path):
    r = run_cli("--demo", "--no-igv", "--caller-sample", "nobody", "--output", tmp_path / "o")
    assert r.returncode != 0 and "Traceback" not in r.stderr and "nobody" in r.stderr


# ── Flags and status ───────────────────────────────────────────────────────────

def test_flags_and_status(demo_result):
    v = by_id(demo_result)
    assert v["V1"]["status"] == "supported" and v["V1"]["flags"] == []
    assert "strand_bias" in v["V2"]["flags"] and v["V2"]["status"] == "flagged"
    assert "low_support" in v["V3"]["flags"] and v["V3"]["status"] == "insufficient"
    assert v["V4"]["status"] == "supported"
    assert v["V5_1"]["status"] == "supported"
    assert "normal_support" in v["V6"]["flags"] and v["V6"]["status"] == "flagged"


# ── Variant selection ──────────────────────────────────────────────────────────

def test_select_by_gene():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf", genes={"DEMO1", "DEMO4"})
    assert sorted(v.id for v in variants) == ["V1", "V4"]


def test_select_by_list():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf", variants_tsv=DEMO / "demo_variants.tsv")
    assert sorted(v.id for v in variants) == ["V1", "V4"]


def test_breakend_selected_by_either_position(tmp_path):
    """Mate breakends collapse to one variant; listing the other end must still select it."""
    lst = tmp_path / "bnd.tsv"
    lst.write_text("chrom\tpos\ndemo2\t5001\n")
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf", variants_tsv=lst)
    assert [v.id for v in variants] == ["V5_1"]


def test_not_evaluated_only_lists_selected_records(tmp_path):
    """A whole-genome VCF has many records the skill cannot evaluate; only report the ones asked for."""
    vcf = tmp_path / "calls.vcf"
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    header = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    extra = ["demo1\t100\tINS1\tA\t<INS>\t.\tPASS\tSVTYPE=INS;GENE=OTHER\tGT\t0/0\t0/1",
             "demo1\t200\tMNV1\tAC\tGT\t.\tPASS\tGENE=DEMO1\tGT\t0/0\t0/1"]
    header.insert(-1, '##ALT=<ID=INS,Description="Insertion">')
    vcf.write_text("\n".join(header + sorted(extra + body, key=lambda l: int(l.split("\t")[1]))) + "\n")
    lst = tmp_path / "one.tsv"
    lst.write_text("demo1\t3000\n")
    _, skipped = iv.load_variants(vcf, variants_tsv=lst)
    assert skipped == []
    _, skipped = iv.load_variants(vcf, genes={"DEMO1"})
    assert [s["id"] for s in skipped] == ["MNV1"]


def test_high_depth_flag_marks_repeat_like_outliers():
    """A site with several times the run's median depth (collapsed repeat) gets high_depth."""
    rows = {f"v{i}": {"tumor": {"depth": d}, "normal": {"depth": 60}, "flags": [], "status": "supported"}
            for i, d in enumerate([60, 70, 65, 55, 470])}
    iv.flag_depth_outliers(rows)
    assert rows["v4"]["flags"] == ["high_depth"] and rows["v4"]["status"] == "flagged"
    assert all(rows[f"v{i}"]["flags"] == [] for i in range(4))
    few = {"a": {"tumor": {"depth": 60}, "normal": {"depth": 60}, "flags": [], "status": "supported"},
           "b": {"tumor": {"depth": 900}, "normal": {"depth": 60}, "flags": [], "status": "supported"}}
    iv.flag_depth_outliers(few)  # under 3 variants there is no meaningful median
    assert few["b"]["flags"] == []


def test_no_reads_is_reported_as_no_coverage(tmp_path):
    """A position the BAMs do not cover is not 'weak evidence': say the region has no reads."""
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    header = [l for l in lines if l.startswith("#")]
    vcf = tmp_path / "gap.vcf"
    vcf.write_text("\n".join(header + ["demo1\t100\tGAP\tA\tC\t.\tPASS\tGENE=DEMO1\tGT\t0/0\t0/1"]) + "\n")
    out = tmp_path / "out"
    r = run_cli("--vcf", vcf, "--tumor", DEMO / "demo_tumor.bam", "--normal", DEMO / "demo_normal.bam",
                "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    v = json.loads((out / "result.json").read_text())["data"]["variants"][0]
    assert v["status"] == "no_coverage" and v["flags"] == ["no_coverage"]
    md = (out / "report.md").read_text()
    assert "no reads" in md.lower() and "NA%" not in md


def test_records_without_id_get_short_readable_ids(tmp_path):
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    rows = []
    for l in lines:
        if l.startswith("#"):
            rows.append(l)
        else:
            f = l.split("\t"); f[2] = "."; rows.append("\t".join(f))
    vcf = tmp_path / "noid.vcf"
    vcf.write_text("\n".join(rows) + "\n")
    ids = {v.id for v in iv.load_variants(vcf)[0]}
    assert "demo1:9000:del12" in ids and "demo1:3000:C>A" in ids
    assert all(len(i) <= 30 for i in ids)


def test_wrong_reference_is_a_clean_error(tmp_path):
    """A FASTA that does not match the BAMs (missing contig or different length) must stop with a message."""
    fa = tmp_path / "other.fa"
    fa.write_text(">demo1\n" + "A" * 500 + "\n")
    pysam.faidx(str(fa))
    r = run_cli("--demo", "--reference", fa, "--no-igv", "--output", tmp_path / "out")
    assert r.returncode != 0 and "Traceback" not in r.stderr
    assert "reference" in r.stderr.lower() and "demo1" in r.stderr


def test_swapped_tumor_and_normal_warning(tmp_path):
    out = tmp_path / "out"
    r = run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", DEMO / "demo_normal.bam", "--normal",
                DEMO / "demo_tumor.bam", "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    assert "swapped" in (out / "report.md").read_text().lower()


def test_not_evaluated_reasons_are_summarised(tmp_path):
    out = tmp_path / "out"
    r = run_cli("--demo", "--no-igv", "--max-variants", "2", "--output", out)
    assert r.returncode == 0, r.stderr
    md = (out / "report.md").read_text()
    assert "4 over --max-variants 2" in md


def test_germline_site_flag():
    """An SNV at a position where the normal already carries another allele (ROS1-style artifact)."""
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf")
    v6 = next(v for v in variants if v.id == "V6")
    third = next(b for b in "ACGT" if b not in (v6.ref, v6.alt))
    probe = iv.Variant(id="p", chrom=v6.chrom, pos=v6.pos, ref=v6.ref, alt=third, kind="snv")
    t = iv.count_variant(probe, DEMO / "demo_tumor.bam")
    n = iv.count_variant(probe, DEMO / "demo_normal.bam")
    flags, _ = iv.flag_variant(probe, t, n)
    assert "germline_site" in flags


def _vcf_with(tmp_path, name, info_header, records):
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    header = [l for l in lines if l.startswith("##")] + info_header + [l for l in lines if l.startswith("#CHROM")]
    p = tmp_path / name
    p.write_text("\n".join(header + records) + "\n")
    return p


def test_genes_from_annovar_fields(tmp_path):
    """ANNOVAR writes gene names to Gene.refGene (several genes joined with \\x3b)."""
    vcf = _vcf_with(tmp_path, "annovar.vcf",
                    ['##INFO=<ID=Gene.refGene,Number=.,Type=String,Description="Gene.refGene annotation provided by ANNOVAR">'],
                    ["demo1\t3000\tA1\tC\tA\t.\tPASS\tGene.refGene=GENEA\tGT:AD\t0/0:60,0\t0/1:46,18",
                     "demo1\t5000\tA2\tA\tT\t.\tPASS\tGene.refGene=GENEB\\x3bGENEA\tGT:AD\t0/0:55,0\t0/1:47,6",
                     "demo1\t7000\tA3\tA\tT\t.\tPASS\tGene.refGene=GENEC\tGT:AD\t0/0:60,0\t0/1:44,2"])
    variants, _ = iv.load_variants(vcf, genes={"GENEB"})
    assert [v.id for v in variants] == ["A2"]
    variants, _ = iv.load_variants(vcf, genes={"GENEA"})
    assert sorted(v.id for v in variants) == ["A1", "A2"]


def test_genes_from_annovar_refgenewithver(tmp_path):
    """ANNOVAR run with -protocol refGeneWithVer writes Gene.refGeneWithVer."""
    vcf = _vcf_with(tmp_path, "annovar_ver.vcf",
                    ['##INFO=<ID=Gene.refGeneWithVer,Number=.,Type=String,Description="Gene.refGeneWithVer annotation provided by ANNOVAR">'],
                    ["demo1\t3000\tB1\tC\tA\t.\tPASS\tGene.refGeneWithVer=GENEE\tGT:AD\t0/0:60,0\t0/1:46,18"])
    assert [v.id for v in iv.load_variants(vcf, genes={"GENEE"})[0]] == ["B1"]


def test_select_by_regions_bed_labels_genes(tmp_path):
    """A BED of gene coordinates selects calls (incl. SVs, by either breakend) and names their gene."""
    bed = tmp_path / "genes.bed"
    bed.write_text("demo1\t2900\t3100\tGENE_A\ndemo2\t4900\t5100\tGENE_B\n")
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf", regions=bed)
    got = {v.id: v.gene for v in variants}
    assert set(got) == {"V1", "V5_1"}  # V5_1's second breakend is demo2:5001
    assert "GENE_A" in got["V1"] and "GENE_B" in got["V5_1"]


def test_survivor_tra_records_are_translocations(tmp_path):
    vcf = _vcf_with(tmp_path, "survivor.vcf",
                    ['##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END">',
                     '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
                     '##ALT=<ID=TRA,Description="Translocation">'],
                    ["demo1\t12000\tT1\tN\t<TRA>\t.\tPASS\tSVTYPE=TRA;CHR2=demo2;END=5001\tGT:AD\t0/0:.\t0/1:."])
    variants, skipped = iv.load_variants(vcf)
    assert skipped == [] and len(variants) == 1
    v = variants[0]
    assert (v.kind, v.chrom2, v.pos2) == ("bnd", "demo2", 5001)
    c = iv.count_variant(v, DEMO / "demo_tumor.bam")
    assert c["alt"] == TRUTH["V5"]["demo_tumor"]["union"]


def test_pass_only_drops_filtered_calls():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf", pass_only=True)
    assert "V3" not in {v.id for v in variants}


def test_unknown_gene_is_a_clean_error(tmp_path):
    r = run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", DEMO / "demo_tumor.bam", "--normal",
                DEMO / "demo_normal.bam", "--genes", "NOTAGENE", "--no-igv", "--output", tmp_path / "out")
    assert r.returncode != 0 and "Traceback" not in r.stderr
    assert "NOTAGENE" in r.stderr


# ── Screenshot helpers (no IGV needed) ─────────────────────────────────────────

def test_blank_image_is_detected_as_not_rendered(tmp_path):
    from PIL import Image
    Image.new("RGB", (1200, 400), "white").save(tmp_path / "blank.png")
    assert not iv.image_rendered(tmp_path / "blank.png")


def test_igv_batch_loads_genome_before_anything_else():
    variants, _ = iv.load_variants(DEMO / "demo_calls.vcf")
    cmds = iv.igv_commands(variants[0], DEMO / "demo_tumor.bam", DEMO / "demo_normal.bam", Path("/tmp/x"))
    assert cmds[0] == "new" and not any(c.startswith("genome") for c in cmds)
    assert any(c.startswith("goto") for c in cmds) and cmds[-1].startswith("snapshot")


def _fake_command(tmp_path, name):
    exe = tmp_path / name
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    return exe


def test_finds_igv_command_on_path(tmp_path, monkeypatch):
    """HPC modules usually provide an `igv` command rather than igv.sh or a .app bundle."""
    exe = _fake_command(tmp_path, "igv")
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr(iv, "MAC_APP_DIRS", [])  # ignore any IGV app installed on this machine
    argv, cwd, why = iv.find_igv(None)
    assert argv == [str(exe)] and why == ""


def test_igv_path_accepts_a_command_name(tmp_path, monkeypatch):
    exe = _fake_command(tmp_path, "igv")
    monkeypatch.setenv("PATH", str(tmp_path))
    argv, _, why = iv.find_igv("igv")
    assert argv == [str(exe)] and why == ""


def test_launcher_that_exits_early_is_not_an_error(tmp_path):
    """Module wrappers may start IGV in the background and return 0 at once; keep waiting for the port."""
    assert iv.launcher_failed(returncode=0) is False
    assert iv.launcher_failed(returncode=1) is True
    assert iv.launcher_failed(returncode=None) is False


def test_missing_igv_skips_screenshots_cleanly(tmp_path):
    r = run_cli("--demo", "--igv-path", tmp_path / "no-such-igv", "--output", tmp_path / "out")
    assert r.returncode == 0, r.stderr
    assert "screenshots skipped" in (tmp_path / "out" / "report.md").read_text().lower()


# ── Output contract: every file in SKILL.md '## Output Structure' is produced ──

def _parse_output_contract(skill_md):
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
        files.append(rel + "/" + name if rel else name)
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced(self, demo_out):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable Output Structure"
        missing = [p for p in promised if not (demo_out / p).exists()]
        assert not missing, "promised but not produced: " + ", ".join(missing)


# ── Variant lists straight from an ANNOVAR-style table ─────────────────────────

def test_variant_list_reads_annovar_table_by_column_name(tmp_path):
    """snv.functional.tsv-style table: columns found by name; ANNOVAR deletions start at the first
    deleted base with alt '-' (VCF POS is one base earlier, the anchor)."""
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    v4 = next(l.split("\t") for l in lines if "\tV4\t" in l)
    deleted = v4[3][1:]
    tbl = tmp_path / "functional.tsv"
    tbl.write_text("sample\tchr\tstart\tend\tref\talt\tgene\n"
                   "demo_tumor\tdemo1\t3000\t3000\tC\tA\tDEMO1\n"
                   f"demo_tumor\tdemo1\t9001\t{9000 + len(deleted)}\t{deleted}\t-\tDEMO4\n")
    variants, skipped = iv.load_variants(DEMO / "demo_calls.vcf", variants_tsv=tbl)
    assert sorted(v.id for v in variants) == ["V1", "V4"] and skipped == []


# ── Copy-number mode ───────────────────────────────────────────────────────────

CNV_TRUTH = TRUTH["CNV"]


@pytest.fixture(scope="module")
def cnv_result(tmp_path_factory):
    out = tmp_path_factory.mktemp("cnv")
    r = run_cli("--demo-cnv", "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    res = json.loads((out / "result.json").read_text())
    return out, {g["gene"]: g for g in res["data"]["cnv"]}


def test_cnv_homozygous_deletion_confirmed_and_low_mapq_noted(cnv_result):
    _, g = cnv_result
    h = g["HOMDEL"]
    assert h["depth_ratio"] <= CNV_TRUTH["HOMDEL"]["ratio_max"]
    assert h["status"] == "supported" and "low_mapq_remaining" in h["flags"]
    assert h["segments"][0]["cn_call"] == "del"


def test_cnv_one_copy_loss_and_neutral(cnv_result):
    _, g = cnv_result
    assert CNV_TRUTH["HEMI"]["ratio_min"] <= g["HEMI"]["depth_ratio"] <= CNV_TRUTH["HEMI"]["ratio_max"]
    assert g["HEMI"]["status"] == "supported"
    assert CNV_TRUTH["NEUTRAL"]["ratio_min"] <= g["NEUTRAL"]["depth_ratio"] <= CNV_TRUTH["NEUTRAL"]["ratio_max"]


def test_cnv_wrong_call_is_flagged(cnv_result):
    _, g = cnv_result
    assert "cnv_disagrees" in g["WRONGDEL"]["flags"] and g["WRONGDEL"]["status"] == "flagged"


def test_cnv_gene_between_segments(cnv_result):
    """Like GENEB in S5: no segment covers the gene, but its depth is still measured."""
    _, g = cnv_result
    assert g["GAPGENE"]["status"] == "no_segment" and g["GAPGENE"]["segments"] == []
    assert g["GAPGENE"]["depth_ratio"] > 0.8


def test_cnv_outputs(cnv_result):
    out, _ = cnv_result
    for rel in ["report.md", "report.html", "result.json", "tables/cnv_depth.tsv", "figures/cnv/HOMDEL.png",
                "reproducibility/commands.sh", "reproducibility/checksums.sha256"]:
        assert (out / rel).exists(), rel
    md = (out / "report.md").read_text()
    assert "Copy number" in md and DISCLAIMER_START in md


def test_cnv_requires_regions(tmp_path):
    r = run_cli("--cnv", DEMO / "demo_cnv.tsv", "--tumor", DEMO / "demo_tumor.bam", "--output", tmp_path / "o")
    assert r.returncode != 0 and "Traceback" not in r.stderr and "--regions" in r.stderr


# ── Small inversions (read-pair orientation) ───────────────────────────────────

def test_small_inversion_counted_from_same_strand_pairs(tmp_path):
    """<INV> under 1 kb: pairs with both mates on the same strand near the breakpoints, plus split reads."""
    seq = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(3000))
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 3000}]}
    raw = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        def add(name, flag, start, mstart, tags=()):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, start
            a.cigarstring, a.mapping_quality = "100M", 60
            a.query_sequence = seq[start:start + 100]
            a.query_qualities = pysam.qualitystring_to_array("?" * 100)
            a.next_reference_id, a.next_reference_start = 0, mstart
            a.set_tags(list(tags))
            out.write(a)
        for k in range(5):   # inversion between 1000 and 1400: forward-forward pairs spanning the left breakpoint
            add(f"ff{k}", 1 | 64 | 32 * 0, 850 + k * 10, 1250 + k * 10)
            add(f"ff{k}", 1 | 128, 1250 + k * 10, 850 + k * 10)
        for k in range(4):   # reverse-reverse pairs spanning the right breakpoint
            add(f"rr{k}", 1 | 64 | 16 | 32, 1150 + k * 10, 1450 + k * 10)
            add(f"rr{k}", 1 | 128 | 16 | 32, 1450 + k * 10, 1150 + k * 10)
        for k in range(20):  # normal proper pairs everywhere
            s = 300 + k * 100
            add(f"n{k}", 1 | 2 | 64 | 32, s, s + 250)
            add(f"n{k}", 1 | 2 | 128 | 16, s + 250, s)
    bam = tmp_path / "inv.bam"
    pysam.sort("-o", str(bam), str(raw)); pysam.index(str(bam))
    v = iv.Variant(id="inv", chrom="c1", pos=1000, ref="N", alt="<INV>", kind="inv", chrom2="c1", pos2=1400)
    c = iv.count_variant(v, bam)
    assert c["discordant"] == 9 and c["alt"] == 9


def test_small_inversion_is_evaluated_not_skipped(tmp_path):
    vcf = _vcf_with(tmp_path, "inv.vcf", ['##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
                                          '##ALT=<ID=INV,Description="Inversion">'],
                    ["demo1\t1000\tI1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=1400\tGT:AD\t0/0:.\t0/1:."])
    variants, skipped = iv.load_variants(vcf)
    assert [v.kind for v in variants] == ["inv"] and skipped == []


# ── IGV versions and start-up (a fake IGV stands in for the real one) ──────────

FAKE_IGV = r'''
import socket, sys, time
args = sys.argv[1:]
port = int(args[args.index("--port") + 1]); home = args[args.index("--igvDirectory") + 1]
genome = args[args.index("-g") + 1]
import os
mode = os.environ["FAKE_IGV_MODE"]   # "old" (IGV <= 2.16), "new", or "dead"
if mode == "dead":
    print("INFO fake IGV: failing to start", flush=True); sys.exit(0)
print("INFO [CommandListener] Listening on port %d" % port, flush=True)
srv = socket.socket(); srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
srv.bind(("127.0.0.1", port)); srv.listen(5)
loaded = time.time() + 1.5
while True:
    c, _ = srv.accept(); f = c.makefile("rw")
    for line in f:
        cmd = line.strip()
        if time.time() > loaded and not getattr(sys, "_said", False):
            print("INFO [GenomeManager] Loading genome: " + genome, flush=True); sys._said = True
        if cmd == "exit":
            f.write("OK\n"); f.flush(); sys.exit(0)
        if cmd == "echo":
            f.write("echo\n")
        elif cmd == "currentGenomePath":
            f.write(("UNKOWN COMMAND: currentGenomePath" if mode == "old" else (genome if time.time() > loaded else "")) + "\n")
        else:
            f.write("OK\n")
        f.flush()
'''


def _fake_igv(tmp_path, monkeypatch, mode):
    script = tmp_path / "fake_igv.py"
    script.write_text(FAKE_IGV)
    monkeypatch.setenv("FAKE_IGV_MODE", mode)
    return [sys.executable, str(script)]


def test_igv_216_without_currentgenomepath_still_starts(tmp_path, monkeypatch):
    """IGV 2.16 answers 'UNKOWN COMMAND: currentGenomePath'; fall back to its log instead of timing out."""
    argv = _fake_igv(tmp_path, monkeypatch, "old")
    t = time.time()
    s = iv.IGVSession(argv, None, DEMO / "demo_ref.fa", timeout=20)
    try:
        assert time.time() - t < 15
        assert s.can_query_genome is False
        assert s.send("goto demo1:1-100") == "OK"
    finally:
        s.close()


def test_igv_with_currentgenomepath_waits_for_the_genome(tmp_path, monkeypatch):
    argv = _fake_igv(tmp_path, monkeypatch, "new")
    s = iv.IGVSession(argv, None, DEMO / "demo_ref.fa", timeout=20)
    try:
        assert s.can_query_genome is True
        assert s.current_genome() == (DEMO / "demo_ref.fa")
    finally:
        s.close()


def test_igv_start_failure_keeps_the_log(tmp_path, monkeypatch):
    argv = _fake_igv(tmp_path, monkeypatch, "dead")
    keep = tmp_path / "kept" / "igv.log"
    with pytest.raises(RuntimeError):
        iv.IGVSession(argv, None, DEMO / "demo_ref.fa", timeout=4, log_copy=keep)
    assert keep.exists() and "failing to start" in keep.read_text()


def test_igv_timeout_flag():
    assert iv.build_parser().parse_args(["--demo", "--igv-timeout", "600"]).igv_timeout == 600


# ── Copy-number flanks sit outside the segment that covers the gene ────────────

def test_cnv_flanks_placed_outside_a_large_deletion():
    """A 4.5 kb gene inside a 212 kb deletion: flanks inside the deletion would read ~0 (GENED in a xenograft)."""
    seg = [{"contig": "c", "start": 1_000_000, "end": 1_212_000, "log2": -9.0, "cn_call": "del", "loh": False}]
    windows = iv.cnv_flank_windows(1_100_000, 1_104_500, seg, contig_len=5_000_000)
    assert windows and all(b < 1_000_000 or a > 1_212_000 for a, b in windows)
    assert any(b < 1_000_000 for a, b in windows) and any(a > 1_212_000 for a, b in windows)


# ── Gene names come from the variant table when the VCF has none ───────────────

def test_gene_label_taken_from_variant_table(tmp_path):
    tbl = tmp_path / "t.tsv"
    tbl.write_text("sample\tchr\tstart\tend\tref\talt\tgene\ndemo_tumor\tdemo1\t3000\t3000\tC\tA\tMYGENE\n")
    lines = [l for l in (DEMO / "demo_calls.vcf").read_text().splitlines()]
    vcf = tmp_path / "nogene.vcf"
    vcf.write_text("\n".join(l.replace("GENE=DEMO1", "X=1") if "\tV1\t" in l else l for l in lines) + "\n")
    v = iv.load_variants(vcf, variants_tsv=tbl)[0]
    assert [x.gene for x in v] == ["MYGENE"]


# ── Copy number is checked against the sample-wide depth, like GATK's log2 ─────

@pytest.mark.parametrize("seg_log2,obs_log2,agrees", [
    (-0.44, -0.40, True),    # broad shallow loss reproduced (GENEE-like)
    (-0.01, -2.80, False),   # neutral call, but the reads show a deep drop (GENEA in S1)
    (-7.60, -8.00, True),    # homozygous deletion, both deep
    (-7.60, -0.20, False),   # deep deletion called, normal depth
    (-1.50, 0.05, False),    # the demo's wrong call
    (0.05, 0.30, True),
])
def test_cnv_log2_agreement(seg_log2, obs_log2, agrees):
    assert iv.cnv_log2_agrees(seg_log2, obs_log2) is agrees


def test_cnv_depth_log2_reproduces_segment_calls(cnv_result):
    _, g = cnv_result
    assert abs(g["HEMI"]["depth_log2"] - (-1.0)) < 0.4
    assert abs(g["NEUTRAL"]["depth_log2"]) < iv.CNV_LOG2_TOL  # a 60 kb demo contig; real genomes are flatter
    assert g["HOMDEL"]["depth_log2"] < -3
    assert g["HOMDEL"]["baseline_depth"] > 10


def test_cnv_flanks_skip_adjacent_deleted_segments():
    """A deletion split over adjacent segments (GENED-like): flanks go past the whole run."""
    segs = [{"contig": "c", "start": 1_000_001, "end": 1_050_000, "log2": -7.6, "cn_call": "del", "loh": False},
            {"contig": "c", "start": 1_050_001, "end": 1_150_000, "log2": -9.1, "cn_call": "del", "loh": False},
            {"contig": "c", "start": 950_001, "end": 1_000_000, "log2": -6.0, "cn_call": "del", "loh": False},
            {"contig": "c", "start": 700_001, "end": 950_000, "log2": 0.0, "cn_call": "neutral", "loh": False}]
    windows = iv.cnv_flank_windows(1_010_000, 1_014_500, segs[:1], contig_len=5_000_000, all_segments=segs)
    assert windows and all(b <= 950_000 or a > 1_150_000 for a, b in windows)


def test_cnv_broad_segment_keeps_local_flanks():
    """A 65 Mb shallow loss: flanks stay next to the gene instead of jumping 65 Mb away."""
    seg = [{"contig": "c", "start": 1, "end": 65_000_000, "log2": -0.44, "cn_call": "del", "loh": False}]
    windows = iv.cnv_flank_windows(30_000_000, 30_070_000, seg, contig_len=100_000_000, all_segments=seg)
    assert all(abs(a - 30_000_000) < 2_000_000 or abs(b - 30_070_000) < 2_000_000 for a, b in windows)


# ── Summary across finished runs ───────────────────────────────────────────────

@pytest.fixture(scope="module")
def summary_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp("reports")
    for args, sub in ((["--demo"], "demo_tumor/variants"), (["--demo-cnv"], "demo_tumor/cnv")):
        r = run_cli(*args, "--no-igv", "--output", root / sub)
        assert r.returncode == 0, r.stderr
    r = run_cli("--summarize", root)
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(root / "summary.tsv"), delimiter="\t"))
    return root, rows


def test_summary_lists_every_call_with_a_concordance_label(summary_dir):
    _, rows = summary_dir
    by = {(r["gene"], r["type"]): r for r in rows}
    assert by[("HOMDEL", "CNV")]["concordance"] == "confirmed"
    assert by[("WRONGDEL", "CNV")]["concordance"] == "questioned"
    assert by[("GAPGENE", "CNV")]["concordance"] == "no call"
    assert by[("DEMO3", "SNV/indel")]["concordance"] == "weak"
    assert by[("DEMO1", "SNV/indel")]["concordance"] == "confirmed"
    assert by[("DEMO5", "SV")]["concordance"] == "confirmed"
    assert all(r["sample"] == "demo_tumor" for r in rows)


def test_summary_html_grid_links_the_evidence(summary_dir):
    root, _ = summary_dir
    page = (root / "summary.html").read_text()
    assert "HOMDEL" in page and "demo_tumor" in page and DISCLAIMER_START in page
    assert "demo_tumor/cnv/report.html" in page


def test_summarize_needs_reports(tmp_path):
    r = run_cli("--summarize", tmp_path)
    assert r.returncode != 0 and "Traceback" not in r.stderr and "result.json" in r.stderr


# ── SVs whose span covers a gene; empty selections are a result, not an error ──

def test_deletion_spanning_a_gene_is_selected(tmp_path):
    """A homozygous deletion removing a whole gene has both breakpoints outside it (GENEA/GENEB in S5)."""
    vcf = _vcf_with(tmp_path, "span.vcf", ['##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
                                           '##ALT=<ID=DEL,Description="Deletion">'],
                    ["demo1\t1000\tD1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=16000\tGT:AD\t0/0:.\t0/1:."])
    bed = tmp_path / "g.bed"
    bed.write_text("demo1\t5000\t6000\tINSIDE\n")
    variants, _ = iv.load_variants(vcf, regions=bed)
    assert [v.id for v in variants] == ["D1"] and "INSIDE" in variants[0].gene


def test_translocation_still_needs_a_breakpoint_in_the_gene(tmp_path):
    bed = tmp_path / "g.bed"
    bed.write_text("demo1\t12500\t13500\tNEAR\n")   # between the V5 breakends' positions, not at one
    with pytest.raises(iv.InputError):
        iv.load_variants(DEMO / "demo_calls.vcf", regions=bed)


def test_no_calls_in_regions_writes_an_empty_report(tmp_path):
    bed = tmp_path / "g.bed"
    bed.write_text("demo1\t100\t200\tEMPTYGENE\n")
    out = tmp_path / "o"
    r = run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", DEMO / "demo_tumor.bam", "--regions", bed,
                "--no-igv", "--output", out)
    assert r.returncode == 0, r.stderr
    res = json.loads((out / "result.json").read_text())
    assert res["data"]["variants"] == [] and res["summary"]["variants_checked"] == 0
    assert "no calls" in (out / "report.md").read_text().lower()


def test_summary_lists_runs_without_calls(tmp_path):
    bed = tmp_path / "g.bed"
    bed.write_text("demo1\t100\t200\tEMPTYGENE\n")
    root = tmp_path / "rep"
    assert run_cli("--vcf", DEMO / "demo_calls.vcf", "--tumor", DEMO / "demo_tumor.bam", "--tumor-name", "S1",
                   "--regions", bed, "--no-igv", "--output", root / "S1" / "sv").returncode == 0
    assert run_cli("--demo-cnv", "--no-igv", "--output", root / "S1" / "cnv").returncode == 0
    assert run_cli("--summarize", root).returncode == 0
    page = (root / "summary.html").read_text()
    assert "no calls" in page.lower() and "S1" in page


# ── Summary shows which copy-number state was checked ──────────────────────────

def test_summary_cnv_rows_name_the_state(summary_dir):
    root, rows = summary_dir
    by = {(r["gene"], r["type"]): r for r in rows}
    assert by[("HOMDEL", "CNV")]["state"] == "del"
    assert by[("NEUTRAL", "CNV")]["state"] == "neutral"
    # no segment: the depth decides what is shown
    assert by[("GAPGENE", "CNV")]["state"] == "no segment; depth neutral"
    page = (root / "summary.html").read_text()
    assert "CNV del" in page and "CNV neutral" in page


@pytest.mark.parametrize("log2,state", [(-10.0, "deep loss"), (-0.9, "loss"), (0.1, "neutral"), (0.7, "gain")])
def test_depth_state_words(log2, state):
    assert iv.depth_state(log2) == state


# ── Reads that are present but ambiguous (MAPQ 0, alternate contigs) ───────────

def test_cnv_ambiguous_mapping_is_not_confirmed(cnv_result):
    """GENED with alternate-haplotype copies: reads are there but map equally well to _alt contigs, so a MAPQ-filtering caller sees a
    false deletion. IGV shows those reads, so the call must not be 'supported'."""
    _, g = cnv_result
    a = g["ALTGENE"]
    assert "ambiguous_mapping" in a["flags"] and a["status"] == "flagged"
    assert a["depth_log2_all"] > -1          # counted like IGV shows them, the gene is covered
    assert a["ambiguous_fraction"] >= 0.9


def test_real_homozygous_deletion_is_not_ambiguous(cnv_result):
    _, g = cnv_result
    assert "ambiguous_mapping" not in g["HOMDEL"]["flags"] and g["HOMDEL"]["status"] == "supported"


# ── The simple answer: does IGV agree with what the callers called? ────────────

def test_agreement_table_lists_only_caller_calls(summary_dir):
    root, _ = summary_dir
    rows = list(csv.DictReader(open(root / "igv_agreement.tsv"), delimiter="\t"))
    ans = {(r["gene"], r["type"]): r for r in rows}
    assert ans[("HOMDEL", "CNV")]["igv_agrees"] == "yes"
    assert ans[("WRONGDEL", "CNV")]["igv_agrees"] == "no"
    assert ans[("ALTGENE", "CNV")]["igv_agrees"] == "no" and "alternate" in ans[("ALTGENE", "CNV")]["why"]
    assert ans[("DEMO1", "SNV/indel")]["igv_agrees"] == "yes"
    assert ans[("DEMO3", "SNV/indel")]["igv_agrees"] == "unclear"
    assert ("NEUTRAL", "CNV") not in ans and ("GAPGENE", "CNV") not in ans   # no call was made there
    assert "Does IGV agree" in (root / "summary.html").read_text()


# ── Which tool made each call, and what IGV shows, in plain words ──────────────

@pytest.mark.parametrize("header,tool", [
    ('##GATKCommandLine=<ID=Mutect2,CommandLine="Mutect2 -R ref.fa">', "Mutect2"),
    ("##source=SURVIVOR", "SURVIVOR"),
    ('##DRAGENCommandLine=<ID=dragen,Version="SW: 4.4">', "DRAGEN"),
    ("##cmdline=/usr/local/bin/configManta.py --tumorBam t.bam", "Manta"),
    ("##purpleVersion=4.3", "PURPLE"),
    ("##source=strelka", "Strelka"),
])
def test_vcf_tool_detected_from_header(tmp_path, header, tool):
    lines = (DEMO / "demo_calls.vcf").read_text().splitlines()
    vcf = tmp_path / "x.vcf"
    vcf.write_text("\n".join([lines[0], header] + lines[1:]) + "\n")
    assert iv.vcf_tool(vcf) == tool


def test_vcf_tool_falls_back_to_file_name():
    assert iv.vcf_tool(DEMO / "demo_calls.vcf") == "demo_calls.vcf"


def test_agreement_table_says_what_igv_shows(summary_dir):
    root, _ = summary_dir
    rows = {(r["gene"], r["type"]): r for r in csv.DictReader(open(root / "igv_agreement.tsv"), delimiter="\t")}
    assert rows[("DEMO1", "SNV/indel")]["tool"] == "demo_calls.vcf"
    assert "18 of 64 reads carry" in rows[("DEMO1", "SNV/indel")]["igv_shows"]
    assert "support the rearrangement" in rows[("DEMO5", "SV")]["igv_shows"]
    alt = rows[("ALTGENE", "CNV")]["igv_shows"]
    assert "counting all reads" in alt and "ambiguous" in alt
    assert "copy-number segments" in rows[("HOMDEL", "CNV")]["tool"]


# ── One overview image per gene and sample ─────────────────────────────────────

def test_reports_record_their_inputs(demo_result):
    inp = demo_result["data"]["inputs"]
    assert inp["tumor"].endswith("demo_tumor.bam") and Path(inp["tumor"]).is_absolute()


def test_overview_plan_one_entry_per_sample_gene(summary_dir, tmp_path):
    root, rows = summary_dir
    bed = tmp_path / "genes.bed"
    bed.write_text("demo1\t2900\t3100\tDEMO1\ndemo3\t20999\t24000\tHOMDEL\ndemo1\t11900\t12100\tDEMO5\n")
    plan = iv.overview_plan(iv._summary_rows(root), bed)
    by = {(p["sample"], p["gene"]): p for p in plan}
    assert set(by) == {("demo_tumor", "DEMO1"), ("demo_tumor", "HOMDEL"), ("demo_tumor", "DEMO5")}
    h = by[("demo_tumor", "HOMDEL")]
    assert h["chrom"] == "demo3" and h["start"] < 21000 and h["end"] > 24000
    assert h["bam"].endswith("demo_tumor.bam")
    assert any("del" in b for b in h["bed_lines"])                 # the GATK segment is marked
    assert any("HOMDEL" in c and "agrees" in c for c in h["caption"]) or any("IGV agrees" in c for c in h["caption"])
    d5 = by[("demo_tumor", "DEMO5")]
    assert any("12000" in b or "11999" in b for b in d5["bed_lines"])  # the SV breakpoint is marked


def test_overview_needs_regions(tmp_path, summary_dir):
    root, _ = summary_dir
    r = run_cli("--summarize", root, "--overview")
    assert r.returncode != 0 and "Traceback" not in r.stderr and "--regions" in r.stderr


# ── Captions for genes without a call must not claim agreement blindly ─────────

@pytest.fixture
def odd_cnv_rows(summary_dir, tmp_path):
    """S2 GENED: no segment and every read MAPQ 0. S1 GENEA: GATK neutral but depth shows a loss."""
    root, _ = summary_dir
    d = json.loads((root / "demo_tumor" / "cnv" / "result.json").read_text())
    for c in d["data"]["cnv"]:
        if c["gene"] == "GAPGENE":
            c.update(flags=["ambiguous_mapping"], depth_log2=-10.0, ambiguous_fraction=1.0)
        if c["gene"] == "NEUTRAL":
            c.update(flags=["cnv_disagrees"], depth_log2=-2.8, status="flagged")
    run = tmp_path / "rep" / "S" / "cnv"
    run.mkdir(parents=True)
    (run / "result.json").write_text(json.dumps(d))
    return {r["gene"]: r for r in iv._summary_rows(tmp_path / "rep")}


def _caption(rows, gene, tmp_path):
    bed = tmp_path / "g.bed"
    r = rows[gene]
    bed.write_text(f"{r['_chrom']}\t{r['_gstart'] - 1}\t{r['_gend']}\t{gene}\n")
    return " ".join(iv.overview_plan([r], bed)[0]["caption"])


def test_no_segment_with_ambiguous_reads_is_not_called_a_deep_loss(odd_cnv_rows, tmp_path):
    r = odd_cnv_rows["GAPGENE"]
    assert "deep loss" not in r["state"] and "ambiguous" in r["state"]
    cap = _caption(odd_cnv_rows, "GAPGENE", tmp_path)
    assert "the reads agree" not in cap and "deep loss" not in cap and "ambiguous" in cap


def test_neutral_call_contradicted_by_depth_is_a_disagreement(odd_cnv_rows, tmp_path):
    cap = _caption(odd_cnv_rows, "NEUTRAL", tmp_path)
    assert "the reads agree" not in cap and "IGV agrees: no" in cap
    ag = iv._agreement(odd_cnv_rows["NEUTRAL"])
    assert ag and ag["igv_agrees"] == "no" and "does not match" in ag["why"]


# ── The summary reuses the IGV the per-sample runs found ───────────────────────

def test_reports_record_the_igv_launcher(tmp_path):
    exe = _fake_command(tmp_path, "igv")
    out = tmp_path / "o"
    assert run_cli("--demo", "--no-igv", "--igv-path", exe, "--output", out).returncode == 0
    assert json.loads((out / "result.json").read_text())["data"]["inputs"]["igv"] == str(exe.resolve())


def test_overview_reuses_the_recorded_igv(tmp_path, monkeypatch):
    """`module load igv` was active for the runs but not in the terminal running --summarize."""
    exe = _fake_command(tmp_path, "igv")
    monkeypatch.setenv("PATH", str(tmp_path / "nothing"))
    monkeypatch.setattr(iv, "MAC_APP_DIRS", [])
    assert iv.overview_launcher(None, [{"igv": str(exe)}]) == str(exe)
    assert iv.overview_launcher(None, [{"igv": str(tmp_path / "gone")}]) is None
    assert iv.overview_launcher("/my/igv", [{"igv": str(exe)}]) == "/my/igv"   # an explicit path wins


def test_overview_without_igv_says_how_to_fix_it(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "nothing"))
    monkeypatch.setattr(iv, "MAC_APP_DIRS", [])
    done, note = iv.take_overviews([{"igv": None, "ref": None}], tmp_path, None)
    assert done == {} and "module load igv" in note and "--igv-path" in note


# ── The same SNV in several unrelated samples is not a tumor mutation ──────────

def _snv_row(sample, loc="chrX:1000100", call="chrX:1000100 T>G", status="supported"):
    return {"sample": sample, "gene": "GENEC", "type": "SNV/indel", "tool": "Mutect2", "state": "", "call": call,
            "location": loc, "igv_shows": "26 of 86 reads carry G", "flags": "", "status": status,
            "concordance": iv.CONCORDANCE[status], "report": "r.html"}


def test_recurrent_snv_is_marked_not_called_a_mutation():
    """GENEC chrX:1000100 T>G in all 5 samples: the reads carry it, but it is germline or artifact."""
    rows = [_snv_row(s) for s in ("A", "B", "C", "D", "E")] + [_snv_row("A", "chrX:1", "chrX:1 C>T")]
    iv.mark_recurrent(rows)
    ag = iv._agreement(rows[0])
    assert ag["igv_agrees"] == "in reads, recurrent"
    assert "5 of 5 samples" in ag["why"] and "germline or artifact" in ag["why"]
    assert rows[0]["concordance"] == "recurrent"
    assert iv._agreement(rows[-1])["igv_agrees"] == "yes"      # private to one sample: unchanged


def test_two_samples_sharing_a_variant_is_not_enough():
    rows = [_snv_row(s) for s in ("A", "B")] + [_snv_row("C", "chrX:1", "chrX:1 C>T")]
    iv.mark_recurrent(rows)
    assert all(iv._agreement(r)["igv_agrees"] == "yes" for r in rows)


def test_recurrent_flag_reaches_the_summary_page(tmp_path):
    root = tmp_path / "rep"
    for s in ("S1", "S2", "S3"):
        assert run_cli("--demo", "--no-igv", "--tumor-name", s, "--output", root / s / "snv").returncode == 0
    assert run_cli("--summarize", root).returncode == 0
    rows = list(csv.DictReader(open(root / "igv_agreement.tsv"), delimiter="\t"))
    assert any(r["igv_agrees"] == "in reads, recurrent" for r in rows)
    assert "recurrent" in (root / "summary.html").read_text()


# ── Your heatmap vs IGV: one line per heatmap cell ─────────────────────────────

def _heatmap(tmp_path, cells):
    p = tmp_path / "heatmap.csv"
    p.write_text("sample,gene,has_snv,alteration,pathway\n" +
                 "".join(f"{s},{g},NA,{a},x\n" for s, g, a in cells))
    return p


def test_heatmap_vs_igv_says_where_they_match_and_differ(summary_dir, tmp_path):
    root, _ = summary_dir
    hm = _heatmap(tmp_path, [("demo_tumor", "HOMDEL", "DEL"), ("demo_tumor", "WRONGDEL", "DEL"),
                             ("demo_tumor", "ALTGENE", "DEL"), ("demo_tumor", "NEUTRAL", "WT"),
                             ("demo_tumor", "DEMO1", "WT"), ("demo_tumor", "DEMO5", "SV"),
                             ("demo_tumor", "GENEG", "WT")])
    by = {(c["sample"], c["gene"]): c for c in iv.heatmap_vs_igv(iv._summary_rows(root), hm)}
    assert by[("demo_tumor", "HOMDEL")]["match"] == "matches"
    assert by[("demo_tumor", "WRONGDEL")]["match"] == "differs"
    alt = by[("demo_tumor", "ALTGENE")]
    assert alt["match"] == "differs" and "alternate" in alt["why"]         # the GENED case
    assert by[("demo_tumor", "NEUTRAL")]["match"] == "matches"
    d1 = by[("demo_tumor", "DEMO1")]
    assert d1["match"] == "differs" and "does not show" in d1["why"]       # a supported SNV the heatmap calls WT
    assert by[("demo_tumor", "DEMO5")]["match"] == "matches"
    assert ("demo_tumor", "GENEG") not in by                                 # never checked with IGV


def test_heatmap_wt_with_a_recurrent_variant_is_explained(tmp_path):
    """GENEC: WT in the heatmap, a G in the reads. Not a contradiction: the variant is in every sample."""
    root = tmp_path / "rep"
    for s in ("S1", "S2", "S3"):
        assert run_cli("--demo", "--no-igv", "--tumor-name", s, "--output", root / s / "snv").returncode == 0
    rows = iv._summary_rows(root)
    iv.mark_recurrent(rows)
    hm = _heatmap(tmp_path, [("S1", "DEMO1", "WT")])
    c = iv.heatmap_vs_igv(rows, hm)[0]
    assert c["match"] == "matches"
    assert "recurrent" in c["igv_found"] and "germline or artifact" in c["why"]


def test_heatmap_flag_writes_the_comparison(summary_dir, tmp_path):
    root, _ = summary_dir
    hm = _heatmap(tmp_path, [("demo_tumor", "HOMDEL", "DEL"), ("demo_tumor", "ALTGENE", "DEL")])
    out = tmp_path / "s"
    r = run_cli("--summarize", root, "--heatmap", hm, "--output", out)
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(out / "heatmap_vs_igv.tsv"), delimiter="\t"))
    assert {r["gene"]: r["match"] for r in rows} == {"HOMDEL": "matches", "ALTGENE": "differs"}
    assert "Your heatmap vs IGV" in (out / "summary.html").read_text()


def test_heatmap_needs_sample_gene_alteration(summary_dir, tmp_path):
    root, _ = summary_dir
    bad = tmp_path / "bad.csv"
    bad.write_text("a,b\n1,2\n")
    r = run_cli("--summarize", root, "--heatmap", bad)
    assert r.returncode != 0 and "Traceback" not in r.stderr and "alteration" in r.stderr


# ── The comparison follows the heatmap's own rules (heatmap table) ─────────────────

def _cn_row(sample, gene, state, log2, size, depth, status="supported", flags="", chrom="chr2"):
    seg = [] if state.startswith("no segment") else \
        [{"contig": chrom, "start": 1, "end": size, "log2": log2, "cn_call": state, "loh": False}]
    return {"sample": sample, "gene": gene, "type": "CNV", "state": state, "call": f"{state} log2 {log2}",
            "flags": flags, "status": status, "igv_shows": "", "_depth": depth, "_segs": seg, "_chrom": chrom}


def _sv_row(sample, gene):
    return {"sample": sample, "gene": gene, "type": "SV", "state": "BND", "call": "chrX:1 <-> chr6:2 breakend",
            "flags": "", "status": "supported", "igv_shows": "40 reads support the rearrangement"}


def _compare(tmp_path, rows, cells):
    return {c["gene"]: c for c in iv.heatmap_vs_igv(rows, _heatmap(tmp_path, cells))}


def test_broad_shallow_cnv_is_left_out_by_design(tmp_path):
    """GENEE: a 65 Mb segment at log2 -0.44. The heatmap shows DEL only if < 1 Mb or log2 < -2."""
    c = _compare(tmp_path, [_cn_row("S", "GENEE", "del", -0.44, 65_000_000, "loss")], [("S", "GENEE", "WT")])
    assert c["GENEE"]["match"] == "matches" and "by design" in c["GENEE"]["why"]


def test_focal_cnv_missing_from_the_heatmap_differs(tmp_path):
    c = _compare(tmp_path, [_cn_row("S", "G", "del", -1.4, 300_000, "loss")], [("S", "G", "WT")])
    assert c["G"]["match"] == "differs" and "does not show" in c["G"]["why"]


def test_sv_takes_priority_over_copy_number_in_a_cell(tmp_path):
    """S2 GENEB: SV and a broad DEL; S1 GENEA: SV label and a deep depth-only loss. One label per cell."""
    rows = [_sv_row("S", "GENEB"), _cn_row("S", "GENEB", "del", -2.6, 200_000, "deep loss")]
    c = _compare(tmp_path, rows, [("S", "GENEB", "SV")])
    assert c["GENEB"]["match"] == "matches" and "priority" in c["GENEB"]["why"]


def test_one_x_copy_is_not_a_deletion(tmp_path):
    """GENEC in a male sample: depth about half the autosomes (log2 ~ -1), no GATK segment."""
    rows = [_cn_row("S", "GENEC", "no segment; depth loss", None, 0, "loss", status="no_segment", chrom="chrX")]
    c = _compare(tmp_path, rows, [("S", "GENEC", "WT")])
    assert c["GENEC"]["match"] == "matches"


def test_depth_only_deep_loss_counts_as_a_deletion(tmp_path):
    """S5 GENEB: no segment, but the reads are gone."""
    rows = [_cn_row("S", "GENEB", "no segment; depth deep loss", None, 0, "deep loss", status="no_segment")]
    assert _compare(tmp_path, rows, [("S", "GENEB", "DEL")])["GENEB"]["match"] == "matches"


def test_wt_over_ambiguous_reads_explains_itself(tmp_path):
    """S2 GENED: no segment, all reads MAPQ 0. WT is right, but not because 'the reads show no change'."""
    rows = [_cn_row("S", "GENED", "no segment; reads ambiguous", None, 0, "deep loss", status="no_segment",
                    flags="ambiguous_mapping")]
    c = _compare(tmp_path, rows, [("S", "GENED", "WT")])["GENED"]
    assert c["match"] == "matches" and "ambiguous" in c["why"] and "no change" not in c["why"]


def test_weak_copy_of_a_recurrent_snv_is_recurrent_too():
    """S1 GENEC: 2 reads at the site that the other 4 samples carry strongly."""
    rows = [_snv_row(s) for s in ("A", "B", "C", "D")] + [_snv_row("E", status="insufficient")]
    iv.mark_recurrent(rows)
    assert "recurrent" in rows[-1]["flags"]


def test_matching_cell_gets_a_one_line_why_and_full_details(tmp_path):
    """S5 GENEC: WT, two recurrent SNVs and a mild chrX drop. The answer fits on one line."""
    rows = [_snv_row(s, loc=f"chrX:{p}", call=f"chrX:{p} T>G") for s in "ABCDE" for p in (1000100, 1000106)]
    rows.append(_cn_row("A", "GENEC", "no segment; depth loss", None, 0, "loss", status="no_segment", chrom="chrX"))
    iv.mark_recurrent(rows)
    c = _compare(tmp_path, rows, [("A", "GENEC", "WT")])["GENEC"]
    assert c["match"] == "matches" and c["why"].startswith("WT is right")
    assert len(c["why"]) < 160 and "5 of 5 samples" in c["why"]
    assert "chrX:1000100 T>G" in c["details"] and "mild loss" in c["details"]


# ── Small copy-number changes cannot be confirmed from depth ───────────────────

def test_small_change_is_not_called_real(tmp_path):
    """S4 GENEF: GATK amp at log2 0.16. Depth within 0.5 of it is also within 0.5 of normal."""
    r = _cn_row("S", "GENEF", "amp", 0.16, 11_500_000, "neutral")
    kind, verdict, text = iv._evidence(r)
    assert verdict == "not supported" and "normal range" in text and "real" not in text


def test_small_change_is_unclear_in_the_callers_table():
    r = _cn_row("S", "GENEF", "amp", 0.16, 11_500_000, "neutral")
    r.update(igv_shows="", tool="GATK", report="r.html")
    ag = iv._agreement(r)
    assert ag["igv_agrees"] == "unclear" and "too small" in ag["why"]


def test_depth_deeper_than_the_segment_is_a_focal_loss(tmp_path):
    """S2 GENEA: GATK del -0.54 over a broad segment, but the gene itself reads much lower. Not 'right to ignore'."""
    r = _cn_row("S", "GENEA", "del", -0.54, 9_000_000, "deep loss", status="flagged", flags="cnv_disagrees")
    r["_depth_log2"] = -3.1
    kind, verdict, text = iv._evidence(r)
    assert (kind, verdict) == ("DEL", "depth only") and "-3.1" in text
    c = _compare(tmp_path, [r], [("S", "GENEA", "WT")])["GENEA"]
    assert c["match"] == "differs"


def test_depth_disagreement_names_the_measured_value():
    r = _cn_row("S", "GENEA", "del", -0.54, 9_000_000, "neutral", status="flagged", flags="cnv_disagrees")
    r["_depth_log2"] = 0.05
    kind, verdict, text = iv._evidence(r)
    assert verdict == "not supported" and "0.05" in text


def test_small_segment_with_a_clear_depth_change_is_confirmed(tmp_path):
    """S2 GENEF: GATK amp 0.46 (small), but the gene reads 0.93 above the sample: the gain is visible."""
    r = _cn_row("S", "GENEF", "amp", 0.46, 3_600_000, "gain")
    r["_depth_log2"] = 0.93
    kind, verdict, text = iv._evidence(r)
    assert verdict == "hidden" and "cannot confirm" not in text and "broad" in text
    r.update(igv_shows="", tool="GATK", report="r.html")
    assert iv._agreement(r)["igv_agrees"] == "yes"


# ── Copy number inside a gene: steps and per-segment parts (real-data review) ─

def test_depth_step_found_at_a_breakpoint_inside_the_gene():
    """S1 GENEC: ~23x, then ~60x from a translocation breakpoint on."""
    st = iv.depth_step([23.0] * 30 + [60.0] * 10, start=1000, binsize=100)
    assert st and st["pos"] == 1000 + 30 * 100 and st["log2_change"] > 1


def test_one_copy_step_is_found_too():
    """S1 GENEF: ~45x then ~70x at the tandem duplication start."""
    st = iv.depth_step([44.0, 46.0] * 10 + [69.0, 71.0] * 10, start=0, binsize=10)
    assert st and 0.5 < st["log2_change"] < 0.8


def test_noise_and_a_single_dip_are_not_a_step():
    vals = [58.0, 62.0, 60.0, 57.0, 63.0] * 8
    vals[20] = 20.0    # a one-bin mappability dip (seen at chr2:2,000 kb in three samples)
    assert iv.depth_step(vals, start=0, binsize=10) is None


def test_each_segment_part_is_compared_with_its_own_depth():
    """S2 GENEA: GATK neutral 0.08 then del -0.54 from 5,000 kb; the gene average (+0.02) hides the step."""
    segs = [{"contig": "c", "start": 1, "end": 1500, "log2": 0.08, "cn_call": "neutral", "loh": False},
            {"contig": "c", "start": 1501, "end": 9000, "log2": -0.54, "cn_call": "del", "loh": False}]
    parts = iv.segment_parts([62.0] * 5 + [40.0] * 5, gstart=1001, binsize=100, segs=segs, baseline=55.7)
    assert [p["cn_call"] for p in parts] == ["neutral", "del"]
    assert all(p["agrees"] for p in parts)
    assert parts[1]["depth_log2"] < -0.3


def test_part_of_a_gene_deeply_lost_is_a_deletion(tmp_path):
    """S1 GENEB: 5' third of the gene at ~8x, the rest ~95x (sample 59x); gene average +0.14."""
    r = _cn_row("S", "GENEB", "neutral", -0.01, 50_000_000, "neutral")
    r.update(_depth_log2=0.14, _step={"pos": 2_000_000, "left": 8.0, "right": 95.0, "log2_change": 3.6,
                                      "left_log2": -2.9, "right_log2": 0.69})
    kind, verdict, text = iv._evidence(r)
    assert (kind, verdict) == ("DEL", "depth only") and "2,000,000" in text and "part of the gene" in text


def test_step_replaces_the_one_x_copy_explanation():
    r = _cn_row("S", "GENEC", "no segment; depth loss", None, 0, "loss", status="no_segment", chrom="chrX")
    r.update(_depth_log2=-0.85, _step={"pos": 1_200_000, "left": 23.0, "right": 60.0, "log2_change": 1.38,
                                       "left_log2": -1.36, "right_log2": 0.02})
    kind, verdict, text = iv._evidence(r)
    assert "X copy" not in text and "1,200,000" in text


def test_called_change_the_depth_cannot_see_is_not_confirmed():
    """S2 GENEE: GATK del -0.55, depth -0.16 (normal range)."""
    r = _cn_row("S", "GENEE", "del", -0.55, 64_800_000, "neutral")
    r["_depth_log2"] = -0.16
    kind, verdict, text = iv._evidence(r)
    assert verdict == "not supported" and "normal range" in text


def test_duplicate_sv_records_are_merged():
    """S2: Manta BND chr2:5,000,438<->chr3:3,000,319 is SURVIVOR's chr3:3,000,319<->chr2:5,000,441."""
    a = {"sample": "S", "gene": "GENEA", "type": "SV", "_chrom": "chr5", "_pos": 3000319, "_chrom2": "chr17",
         "_pos2": 5000441, "call": "a"}
    b = {"sample": "S", "gene": "GENEA", "type": "SV", "_chrom": "chr17", "_pos": 5000438, "_chrom2": "chr5",
         "_pos2": 3000319, "call": "b"}
    c = {"sample": "S", "gene": "GENEA", "type": "SV", "_chrom": "chr5", "_pos": 3001421, "_chrom2": "chr17",
         "_pos2": 4999984, "call": "c"}
    assert [r["call"] for r in iv.dedupe_svs([a, b, c])] == ["a", "c"]


def test_recurrent_snv_at_low_allele_fractions_points_to_artifact():
    """GENEC pair: 8-30% in all 5 samples. Germline would be ~50% or ~100%."""
    rows = [_snv_row(s) for s in "ABCDE"]
    for r, v in zip(rows, (30.2, 22.7, 8.0, 16.7, 13.7)):
        r["_vaf"] = v
    iv.mark_recurrent(rows)
    assert "too low for an inherited" in iv._agreement(rows[0])["why"]
