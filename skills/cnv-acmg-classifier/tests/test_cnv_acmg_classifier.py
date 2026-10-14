"""Red/green TDD suite for the CNV ACMG classifier.

Encodes the ClinGen / ACMG 2019 (Riggs et al. 2020) copy-number point framework
as deterministic expectations. Scoring is ADDITIVE across Sections 1-5 (no
terminal short-circuit on 2A/2F); Section 3 gene-count is omitted only when a
complete established call (2A/2F) is made; partial-overlap sub-calls
(2C-1/2C-2/2D/2E) are derived from breakpoint geometry, not a free-text flag.
"""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "cnv_acmg_classifier.py"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("cnv_acmg_classifier", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- reference fixtures (inline, independent of demo data) ---

# TP53: GRCh38 chr17:7,668,421-7,687,550, minus strand; CDS within the gene.
DOSAGE_MAP = [
    {"chrom": "17", "start": 7668421, "end": 7687550, "name": "TP53", "hi_score": 3,
     "ts_score": 0, "benign": False, "element_type": "gene", "strand": "-",
     "cds_start": 7669000, "cds_end": 7687000},
    {"chrom": "22", "start": 18900000, "end": 21500000, "name": "22q11.2", "hi_score": 3,
     "ts_score": 3, "benign": False, "element_type": "region", "strand": "+",
     "cds_start": None, "cds_end": None},
    {"chrom": "1", "start": 152000000, "end": 152100000, "name": "BENIGN_DEMO", "hi_score": 0,
     "ts_score": 0, "benign": True, "element_type": "region", "strand": "+",
     "cds_start": None, "cds_end": None},
]

GENE_MODEL = [
    {"chrom": "17", "start": 7668421, "end": 7687550, "gene": "TP53"},
    {"chrom": "1", "start": 152030000, "end": 152060000, "gene": "BENIGNGENE1"},
    {"chrom": "2", "start": 50100000, "end": 50200000, "gene": "GENE_X"},
]


def cnv(**kw):
    base = {"cnv_id": "x", "chrom": "1", "start": 1, "end": 2, "type": "DEL",
            "inheritance": "unknown", "case_evidence_points": 0.0}
    base.update(kw)
    return base


def codes_of(res):
    return {e["code"] for e in res["evidence"]}


# --------------------------------------------------------------------------- #
# Thresholds & type normalisation
# --------------------------------------------------------------------------- #
def test_classify_score_thresholds():
    m = load_module()
    assert m.classify_score(0.99) == "Pathogenic"
    assert m.classify_score(1.45) == "Pathogenic"
    assert m.classify_score(0.98) == "Likely pathogenic"
    assert m.classify_score(0.90) == "Likely pathogenic"
    assert m.classify_score(0.89) == "Variant of uncertain significance"
    assert m.classify_score(0.0) == "Variant of uncertain significance"
    assert m.classify_score(-0.89) == "Variant of uncertain significance"
    assert m.classify_score(-0.90) == "Likely benign"
    assert m.classify_score(-0.98) == "Likely benign"
    assert m.classify_score(-0.99) == "Benign"
    assert m.classify_score(-1.00) == "Benign"


def test_normalize_type():
    m = load_module()
    for v in ("DEL", "LOSS", "loss", "<DEL>", "CN0", "CN1"):
        assert m.normalize_type(v) == "loss"
    for v in ("DUP", "GAIN", "gain", "<DUP>", "CN3", "CN4"):
        assert m.normalize_type(v) == "gain"
    try:
        m.normalize_type("INV")
    except ValueError as exc:
        assert "type" in str(exc).lower()
    else:
        raise AssertionError("invalid SV type should raise")


def test_cn1_on_sex_chromosome_is_ambiguous():
    m = load_module()
    assert m.normalize_type("CN1", "chr2") == "loss"
    assert m.normalize_type("DEL", "chrX") == "loss"
    try:
        m.normalize_type("CN1", "chrX")
    except ValueError as exc:
        assert "sex chromosome" in str(exc).lower() or "hemizygous" in str(exc).lower()
    else:
        raise AssertionError("CN1 on chrX should be rejected as ambiguous")


# --------------------------------------------------------------------------- #
# Section 1
# --------------------------------------------------------------------------- #
def test_no_protein_coding_content_is_1B():
    m = load_module()
    res = m.score_cnv(cnv(chrom="5", start=1000, end=2000, type="DEL"), DOSAGE_MAP, GENE_MODEL)
    assert "1B" in codes_of(res)
    assert res["total_score"] == -0.60
    assert res["classification"] == "Variant of uncertain significance"


# --------------------------------------------------------------------------- #
# Section 2 — complete overlaps (additive, NOT terminal)
# --------------------------------------------------------------------------- #
def test_full_loss_of_established_HI_gene_is_pathogenic():
    m = load_module()
    res = m.score_cnv(cnv(chrom="17", start=7660000, end=7695000, type="DEL"),
                      DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2A" in c
    assert not any(x.startswith("3") for x in c)   # Section 3 omitted on a complete 2A
    assert res["total_score"] == 1.00
    assert res["classification"] == "Pathogenic"


def test_2A_loss_inherited_sums_to_vus():
    # ClinGen worked example: 2A (+1.00) + 5B inherited-unaffected (-0.30) = +0.70 = VUS.
    # Sections are additive; 2A is NOT terminal.
    m = load_module()
    res = m.score_cnv(cnv(chrom="17", start=7660000, end=7695000, type="DEL",
                          inheritance="inherited"), DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2A" in c and "5B" in c
    assert res["total_score"] == 0.70
    assert res["classification"] == "Variant of uncertain significance"


def test_2A_gain_de_novo_evidence_is_retained():
    # de novo confirmatory evidence must be summed, not discarded: 2A + 5A = 1.45.
    m = load_module()
    res = m.score_cnv(cnv(chrom="22", start=18800000, end=21600000, type="DUP",
                          inheritance="de_novo"), DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2A" in c and "5A" in c
    assert res["total_score"] == 1.45
    assert res["classification"] == "Pathogenic"


def test_containment_in_benign_region_is_benign():
    m = load_module()
    res = m.score_cnv(cnv(chrom="1", start=152030000, end=152070000, type="DEL"),
                      DOSAGE_MAP, GENE_MODEL)
    assert "2F" in codes_of(res)
    assert res["total_score"] == -1.00
    assert res["classification"] == "Benign"


# --------------------------------------------------------------------------- #
# Section 2 — partial overlaps from breakpoint geometry
# --------------------------------------------------------------------------- #
def test_partial_5prime_coding_overlap_is_2C1():
    # CNV deletes the 5' end of TP53 (minus strand: 5' end at gene.end) including CDS.
    m = load_module()
    res = m.score_cnv(cnv(chrom="17", start=7680000, end=7700000, type="DEL"),
                      DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2C-1" in c
    assert res["total_score"] == 0.90
    assert res["classification"] == "Likely pathogenic"


def test_partial_5prime_utr_only_is_2C2_zero():
    # CNV covers the 5' end but only the 5'UTR (downstream of cds_end on minus strand),
    # no coding sequence -> 2C-2 (0.00), not credited.
    m = load_module()
    res = m.score_cnv(cnv(chrom="17", start=7687100, end=7700000, type="DEL"),
                      DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2C-2" in c
    assert "2C-1" not in c
    assert res["total_score"] == 0.0
    assert res["classification"] == "Variant of uncertain significance"


def test_partial_overlap_of_region_entry_is_2B_zero():
    m = load_module()
    res = m.score_cnv(cnv(chrom="22", start=21000000, end=22000000, type="DEL"),
                      DOSAGE_MAP, GENE_MODEL)
    c = codes_of(res)
    assert "2B" in c
    assert "2C-1" not in c
    assert res["total_score"] == 0.0


# --------------------------------------------------------------------------- #
# Section 3 — gene count, omitted only on a complete established call
# --------------------------------------------------------------------------- #
def test_section3_omitted_only_on_complete_established_call():
    m = load_module()
    dosage = [{"chrom": "17", "start": 7668421, "end": 7687550, "name": "TP53",
               "hi_score": 3, "ts_score": 0, "benign": False, "element_type": "gene",
               "strand": "-", "cds_start": 7669000, "cds_end": 7687000}]
    genes = [{"chrom": "17", "start": 7000000 + i * 5000, "end": 7000000 + i * 5000 + 2000,
              "gene": f"G{i:03d}"} for i in range(40)]
    res = m.score_cnv(cnv(chrom="17", start=6990000, end=7690000, type="DEL"), dosage, genes)
    c = codes_of(res)
    assert "2A" in c
    assert not any(x.startswith("3") for x in c)
    assert res["total_score"] == 1.00


def test_zero_scoring_partial_overlap_still_gets_gene_count():
    # 2B (clip of an established region, scores 0) + >=35 genes -> 3C (+0.90) -> LP.
    m = load_module()
    dosage = [{"chrom": "5", "start": 100000, "end": 200000, "name": "REGION_X",
               "hi_score": 3, "ts_score": 0, "benign": False, "element_type": "region",
               "strand": "+", "cds_start": None, "cds_end": None}]
    genes = [{"chrom": "5", "start": 150000 + i * 10000, "end": 150000 + i * 10000 + 4000,
              "gene": f"GX{i:03d}"} for i in range(40)]
    res = m.score_cnv(cnv(chrom="5", start=150000, end=700000, type="DEL"), dosage, genes)
    c = codes_of(res)
    assert "2B" in c and "3C" in c
    assert res["total_score"] == 0.90
    assert res["classification"] == "Likely pathogenic"


def test_gene_dense_loss_triggers_3C():
    m = load_module()
    genes = [{"chrom": "19", "start": 52000000 + i * 10000, "end": 52000000 + i * 10000 + 5000,
              "gene": f"ZNFDEMO{i:03d}"} for i in range(40)]
    res = m.score_cnv(cnv(chrom="19", start=51990000, end=52410000, type="DEL"),
                      DOSAGE_MAP, genes)
    assert res["gene_count"] == 40
    assert "3C" in codes_of(res)
    assert res["total_score"] == 0.90


def test_gain_gene_count_tier_is_more_permissive_than_loss():
    m = load_module()
    genes = [{"chrom": "19", "start": 52000000 + i * 10000, "end": 52000000 + i * 10000 + 5000,
              "gene": f"ZNFDEMO{i:03d}"} for i in range(40)]
    res = m.score_cnv(cnv(chrom="19", start=51990000, end=52410000, type="DUP"),
                      DOSAGE_MAP, genes)
    assert "3B" in codes_of(res)
    assert res["total_score"] == 0.45


# --------------------------------------------------------------------------- #
# Sections 4 / 5
# --------------------------------------------------------------------------- #
def test_inherited_unaffected_parent_is_vus():
    m = load_module()
    res = m.score_cnv(cnv(chrom="2", start=50120000, end=50180000, type="DEL",
                          inheritance="inherited"), DOSAGE_MAP, GENE_MODEL)
    assert "5B" in codes_of(res)
    assert res["total_score"] == -0.30


def test_case_evidence_is_clamped():
    m = load_module()
    res = m.score_cnv(cnv(chrom="2", start=50120000, end=50180000, type="DEL",
                          case_evidence_points=5.0), DOSAGE_MAP, GENE_MODEL)
    ev = next(e for e in res["evidence"] if e["code"] == "4")
    assert ev["points"] == 0.90


# --------------------------------------------------------------------------- #
# Boundary tests on SUMMED (not hand-set) scores, incl. float accumulation
# --------------------------------------------------------------------------- #
def test_summed_float_3B_plus_case_lands_exactly_on_090():
    # 3B (0.45) + Section 4 (0.45) must sum to exactly 0.90 after round(total, 2),
    # i.e. no float drift to 0.8999... -> Likely pathogenic (the +0.90 boundary).
    m = load_module()
    genes = [{"chrom": "7", "start": 1000 + i * 1000, "end": 1000 + i * 1000 + 400,
              "gene": f"GG{i:03d}"} for i in range(30)]   # 25-34 genes -> 3B
    res = m.score_cnv(cnv(chrom="7", start=900, end=40000, type="DEL",
                          case_evidence_points=0.45), [], genes)
    assert "3B" in codes_of(res)
    assert res["total_score"] == 0.90
    assert res["classification"] == "Likely pathogenic"


def test_summed_score_crosses_pathogenic_threshold():
    m = load_module()
    # 2C-1 (0.90) + case (0.45) = 1.35 -> Pathogenic (>= 0.99)
    res = m.score_cnv(cnv(chrom="17", start=7680000, end=7700000, type="DEL",
                          case_evidence_points=0.45), DOSAGE_MAP, GENE_MODEL)
    assert res["total_score"] == 1.35
    assert res["classification"] == "Pathogenic"


def test_summed_score_likely_benign_then_benign_boundary():
    m = load_module()
    # 1A + case(-0.90) = -0.90 -> Likely benign (the -0.90 boundary, summed)
    lb = m.score_cnv(cnv(chrom="2", start=50120000, end=50180000, type="DEL",
                         case_evidence_points=-0.90), DOSAGE_MAP, GENE_MODEL)
    assert lb["total_score"] == -0.90
    assert lb["classification"] == "Likely benign"
    # add inherited (-0.30): -1.20 -> Benign (<= -0.99, summed)
    b = m.score_cnv(cnv(chrom="2", start=50120000, end=50180000, type="DEL",
                        case_evidence_points=-0.90, inheritance="inherited"),
                    DOSAGE_MAP, GENE_MODEL)
    assert b["total_score"] == -1.20
    assert b["classification"] == "Benign"


# --------------------------------------------------------------------------- #
# Loaders & CLI
# --------------------------------------------------------------------------- #
def test_load_cnvs_rejects_missing_columns(tmp_path):
    m = load_module()
    bad = tmp_path / "bad.csv"
    bad.write_text("cnv_id,chrom,start\nx,1,10\n", encoding="utf-8")
    try:
        m.load_cnvs(bad)
    except ValueError as exc:
        assert "missing" in str(exc).lower()
    else:
        raise AssertionError("missing columns should fail")


def test_load_vcf_parses_svtype_and_end(tmp_path):
    m = load_module()
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "17\t7660000\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=7695000\n"
        "22\t18800000\tsv2\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=21600000\n",
        encoding="utf-8",
    )
    rows = m.load_cnvs(vcf)
    assert len(rows) == 2
    assert rows[0]["chrom"] == "17" and int(rows[0]["end"]) == 7695000
    assert m.normalize_type(rows[1]["type"]) == "gain"


def test_cli_rejects_malformed_input_without_traceback(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("cnv_id,chrom,start\nx,1,10\n", encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(bad), "--output", str(tmp_path / "out")],
        text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_demo_cli_writes_expected_outputs(tmp_path):
    out = tmp_path / "cnv_out"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
        text=True, capture_output=True, check=True,
    )
    assert "CNV ACMG" in completed.stdout
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report
    assert "Synthetic demo data" in report
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["skill"] == "cnv-acmg-classifier"
    assert result["summary"]["cnv_count"] == 7
    tiers = {c["classification"] for c in result["classifications"]}
    assert {"Pathogenic", "Likely pathogenic", "Variant of uncertain significance",
            "Likely benign", "Benign"} <= tiers
    assert (out / "tables" / "cnv_classifications.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()
    env = (out / "reproducibility" / "environment.yml").read_text(encoding="utf-8")
    assert "conda-forge" in env and "nodefaults" in env
    checksums = (out / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8")
    assert "report.md" in checksums and len(checksums.splitlines()) >= 3


def test_demo_reproducibility_bundle_verifies(tmp_path):
    import hashlib

    out = tmp_path / "cnv_out"
    subprocess.run([sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
                   text=True, capture_output=True, check=True)
    repro = out / "reproducibility"
    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "python skills/cnv-acmg-classifier/cnv_acmg_classifier.py" in commands
    assert "--demo" in commands
    assert str(SKILL_DIR) not in commands
    assert "name: clawbio-cnv-acmg-classifier" in (repro / "environment.yml").read_text(encoding="utf-8")
    lines = (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "tables/cnv_classifications.csv",
            "reproducibility/environment.yml"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest, label
