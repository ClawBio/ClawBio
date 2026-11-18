"""
Tests for fastreer.py
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

SCRIPT = Path(__file__).parent.parent / "fastreer.py"


def run(args, **kwargs):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
        **kwargs,
    )


class TestDemoMode:
    def test_demo_exits_zero(self, tmp_path):
        result = run(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, result.stderr

    def test_demo_creates_report_md(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()

    def test_demo_creates_result_json(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        rj = tmp_path / "result.json"
        assert rj.exists()
        data = json.loads(rj.read_text())
        assert "command" in data
        assert "samples" in data

    def test_demo_creates_reproducibility_bundle(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "reproducibility" / "commands.sh").exists()

    def test_demo_replay_records_demo_not_a_generated_input(self, tmp_path):
        """The demo writes its input into the output dir and falls back to
        synthetic output when fastreeR or Java is missing. Replaying it as
        --input takes the strict path and fails in exactly that case."""
        run(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
        assert "--demo" in text
        assert "demo_samples.vcf" not in text

    def test_window_variants_is_recorded(self, tmp_path):
        run(["--demo", "--output", str(tmp_path), "--window-variants", "50"])
        text = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
        assert "--window-variants" in text and "50" in text

    def test_demo_reproducibility_bundle_is_complete(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        repro = tmp_path / "reproducibility"

        commands_text = (repro / "commands.sh").read_text(encoding="utf-8")
        assert "CLAWBIO_ROOT" in commands_text
        assert "$OUTPUT_DIR" in commands_text
        assert str(tmp_path) not in commands_text

        environment = (repro / "environment.yml").read_text(encoding="utf-8")
        assert "name: clawbio-fastreer" in environment
        assert "fastreer" in environment

        # the Java/pip snapshot stays alongside the conda recipe
        assert (repro / "environment.txt").exists()

        checksum_lines = [
            line
            for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        assert checksum_lines
        labels = set()
        for line in checksum_lines:
            digest, label = line.split("  ", 1)
            assert len(digest) == 64
            labels.add(label)
            assert (tmp_path / label).exists()
        assert {"report.md", "result.json", "tree.nwk"} <= labels

    def test_demo_default_command_is_vcf2tree(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["command"] == "VCF2TREE"

    def test_demo_vcf2dist_command(self, tmp_path):
        result = run(["--demo", "--command", "VCF2DIST", "--output", str(tmp_path)])
        assert result.returncode == 0, result.stderr
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["command"] == "VCF2DIST"

    def test_demo_fasta2dist_command(self, tmp_path):
        result = run(["--demo", "--command", "FASTA2DIST", "--output", str(tmp_path)])
        assert result.returncode == 0, result.stderr
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["command"] == "FASTA2DIST"


class TestOutputContents:
    def test_report_md_contains_sample_count(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        report = (tmp_path / "report.md").read_text()
        assert "sample" in report.lower()

    def test_report_md_contains_disclaimer(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        report = (tmp_path / "report.md").read_text()
        assert "ClawBio" in report

    def test_result_json_has_required_keys(self, tmp_path):
        run(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        for key in ("command", "samples", "input_file", "output_dir"):
            assert key in data, f"Missing key: {key}"

    def test_vcf2tree_creates_nwk_file(self, tmp_path):
        run(["--demo", "--command", "VCF2TREE", "--output", str(tmp_path)])
        assert (tmp_path / "tree.nwk").exists()

    def test_vcf2dist_creates_dist_file(self, tmp_path):
        run(["--demo", "--command", "VCF2DIST", "--output", str(tmp_path)])
        assert (tmp_path / "distances.dist").exists()

    def test_fasta2dist_creates_dist_file(self, tmp_path):
        run(["--demo", "--command", "FASTA2DIST", "--output", str(tmp_path)])
        assert (tmp_path / "distances.dist").exists()


class TestValidation:
    def test_invalid_command_exits_nonzero(self, tmp_path):
        result = run(["--command", "INVALID", "--input", "x.vcf", "--output", str(tmp_path)])
        assert result.returncode != 0

    def test_missing_input_without_demo_exits_nonzero(self, tmp_path):
        result = run(["--command", "VCF2TREE", "--output", str(tmp_path)])
        assert result.returncode != 0

    def test_missing_output_exits_nonzero(self):
        result = run(["--demo"])
        assert result.returncode != 0


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/fastreer/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("fastreer", tmp_path) == snapshot({'demo_samples.vcf': """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4	SAMPLE5
chr1	10000	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1	0/1	0/0
chr1	20000	rs2	C	T	100	PASS	.	GT	0/1	1/1	0/0	0/1	1/1
chr1	30000	rs3	G	A	100	PASS	.	GT	1/1	0/1	0/1	0/0	0/1
chr1	40000	rs4	T	C	100	PASS	.	GT	0/0	0/0	0/1	1/1	0/1
chr1	50000	rs5	A	T	100	PASS	.	GT	0/1	0/1	0/0	0/1	1/1
chr1	60000	rs6	C	G	100	PASS	.	GT	1/1	0/0	0/1	0/1	0/0
chr1	70000	rs7	G	C	100	PASS	.	GT	0/1	1/1	1/1	0/0	0/1
chr1	80000	rs8	T	A	100	PASS	.	GT	0/0	0/1	0/0	1/1	0/1
chr1	90000	rs9	A	C	100	PASS	.	GT	0/1	0/0	1/1	0/1	0/0
chr1	100000	rs10	C	A	100	PASS	.	GT	1/1	0/1	0/1	0/0	1/1
chr1	110000	rs11	G	T	100	PASS	.	GT	0/0	1/1	0/0	0/1	0/1
chr1	120000	rs12	T	G	100	PASS	.	GT	0/1	0/1	1/1	1/1	0/0
chr1	130000	rs13	A	G	100	PASS	.	GT	1/1	0/0	0/1	0/1	1/1
chr1	140000	rs14	C	T	100	PASS	.	GT	0/1	1/1	0/0	0/0	0/1
chr1	150000	rs15	G	A	100	PASS	.	GT	0/0	0/1	0/1	1/1	1/1
chr1	160000	rs16	T	C	100	PASS	.	GT	1/1	0/0	1/1	0/1	0/0
chr1	170000	rs17	A	T	100	PASS	.	GT	0/1	0/1	0/0	1/1	0/1
chr1	180000	rs18	C	G	100	PASS	.	GT	0/0	1/1	0/1	0/0	1/1
chr1	190000	rs19	G	C	100	PASS	.	GT	0/1	0/1	1/1	0/1	0/0
chr1	200000	rs20	T	A	100	PASS	.	GT	1/1	0/0	0/1	1/1	0/1
""", 'report.md': """\
# fastreer Report

**Command**: `VCF2TREE`
**Input**: `demo_samples.vcf`
**Samples**: 5
**Variants**: 20
**Date**: <TIMESTAMP>
**fastreer version**: 2.2.0

> **Note**: fastreer or Java 11+ was not found. Output shown is synthetic demo data representing the expected format.

## Samples (5)

  - SAMPLE1
  - SAMPLE2
  - SAMPLE3
  - SAMPLE4
  - SAMPLE5

## Phylogenetic Tree

**Output format**: Newick
**File**: `tree.nwk`

```newick
((SAMPLE1:0.120,SAMPLE2:0.098):0.045,(SAMPLE3:0.110,(SAMPLE4:0.087,SAMPLE5:0.132):0.062):0.038);
```

The tree was built using hierarchical clustering on a cosine distance matrix
computed from genotype data. Internal branch lengths represent evolutionary
distance; higher values indicate greater genomic divergence between samples.

## Interpretation

The Newick tree can be visualised in tools such as FigTree, iTOL, or Dendroscope. Bootstrap values (if computed) appear at internal nodes.

## Install / Requirements

```bash
pip install fastreer      # Python wrapper
# Java 11+ is also required:
sudo apt install default-jre    # Debian/Ubuntu
brew install openjdk@17          # macOS
```

## Citation

Gkanogiannis A (2016) A scalable assembly-free variable selection algorithm for
biomarker discovery from metagenomes. *BMC Bioinformatics* 17, 311.
https://doi.org/10.1186/s12859-016-1186-3

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "command": "VCF2TREE",
  "input_file": "<OUTPUT_DIR>/demo_samples.vcf",
  "output_dir": "<OUTPUT_DIR>",
  "samples": [
    "SAMPLE1",
    "SAMPLE2",
    "SAMPLE3",
    "SAMPLE4",
    "SAMPLE5"
  ],
  "n_samples": 5,
  "n_variants": 20,
  "output_file": "<OUTPUT_DIR>/tree.nwk",
  "synthetic_demo": true,
  "fastreer_version": "2.2.0",
  "timestamp": "<TIMESTAMP>"
}\
""", 'tree.nwk': '((SAMPLE1:0.120,SAMPLE2:0.098):0.045,(SAMPLE3:0.110,(SAMPLE4:0.087,SAMPLE5:0.132):0.062):0.038);'})
