"""
test_bioqc_mcp.py — Tests for the BioQC MCP skill
==================================================
All tests run offline. Tests confirm metadata, stubs,
reproducibility bundle, and custom chart generation.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml
from inline_snapshot import snapshot

from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import bioqc_mcp  # noqa: E402


def test_frontmatter_loads():
    """Verify that SKILL.md has valid frontmatter."""
    skill_md = SKILL_DIR / "SKILL.md"
    assert skill_md.exists()
    content = skill_md.read_text(encoding='utf-8')
    assert content.startswith("---")
    parts = content.split("---", 2)
    meta = yaml.safe_load(parts[1])
    assert meta["name"] == "bioqc-mcp"
    assert "metadata" in meta
    assert "openclaw" in meta["metadata"]
    assert len(meta["metadata"]["openclaw"]["trigger_keywords"]) >= 3


def test_fastq_file_finder(tmp_path):
    """find_fastq_files scans and retrieves all fastq formats."""
    # Write empty files
    (tmp_path / "read1.fastq").write_text("")
    (tmp_path / "read2.fastq.gz").write_text("")
    (tmp_path / "other.txt").write_text("")

    files = bioqc_mcp.find_fastq_files(str(tmp_path))
    assert len(files) == 2
    assert files[0]["name"] == "read1.fastq"
    assert files[1]["name"] == "read2.fastq.gz"


def test_parse_fastqc_data(tmp_path):
    """parse_fastqc_data correctly extracts metrics from summary and data files."""
    (tmp_path / "summary.txt").write_text("PASS\tBasic Statistics\tfile.fastq\nFAIL\tAdapter Content\tfile.fastq\n")
    (tmp_path / "fastqc_data.txt").write_text("##FastQC\n>>Basic Statistics\tpass\n#Measure\tValue\nTotal Sequences\t1000\n%GC\t48\n>>END_MODULE\n")

    res = bioqc_mcp.parse_fastqc_data(str(tmp_path))
    assert res["success"] is True
    assert res["summary"]["Basic Statistics"] == "PASS"
    assert res["summary"]["Adapter Content"] == "FAIL"
    assert res["metrics"]["Total Sequences"] == "1000"
    assert res["metrics"]["%GC"] == "48"


def test_write_repro_bundle(tmp_path):
    """write_repro_bundle successfully creates three reproducibility files."""
    out = tmp_path / "run"
    out.mkdir()
    (out / "report.md").write_text("# mock\n")
    (out / "multiqc_output").mkdir()
    (out / "multiqc_output" / "multiqc_report.html").write_text("<html></html>")

    bioqc_mcp.write_repro_bundle(out, is_demo=True, input_dir=None)

    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    assert (repro / "checksums.sha256").exists()
    
    cmd_text = (repro / "commands.sh").read_text()
    assert "skills/bioqc-mcp/bioqc_mcp.py" in cmd_text
    assert "--demo" in cmd_text


def test_demo_fastq_generation(tmp_path):
    """make_demo_fastq_and_fastqc successfully scaffolds stubs for SAMPLE_01 & 02."""
    reports = bioqc_mcp.make_demo_fastq_and_fastqc(tmp_path)
    assert len(reports) == 2
    for r in reports:
        assert r.exists()
        assert (r.parent / "summary.txt").exists()
        assert (r.parent / "fastqc_data.txt").exists()


def test_generate_chart_from_data():
    """generate_chart_from_data successfully renders line and bar charts."""
    if not bioqc_mcp.HAS_VIZ:
        pytest.skip("Visualization packages are missing")

    data = [
        {"pos": 1, "val": 35},
        {"pos": 2, "val": 36},
        {"pos": 3, "val": 34}
    ]
    res = bioqc_mcp.generate_chart_from_data("line", data, "Test Line Chart", "Pos", "Val")
    assert res["success"] is True
    assert res["mime_type"] == "image/png"
    assert len(res["image_data"]) > 100


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/bioqc-mcp/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("bioqc-mcp", tmp_path) == snapshot({'fastqc_output/SAMPLE_01_fastqc/SAMPLE_01_fastqc.html': '<html><head><title>FastQC Report - SAMPLE_01</title></head><body><h1>FastQC Report</h1></body></html>', 'fastqc_output/SAMPLE_01_fastqc/fastqc_data.txt': """\
##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Filename	SAMPLE_01.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	1200000
Sequences flagged as poor quality	0
Sequence length	150
%GC	48
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile
1	36	34	38
6	36	34	38
11	35	33	37
16	34	32	36
21	34	32	36
26	35	33	37
31	36	34	38
36	36	34	38
41	35	33	37
46	32	30	34
51	28	26	30
56	22	20	24
61	18	16	20
66	14	12	16
>>END_MODULE
>>Per sequence GC content	pass
#GC Content	Count
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
21	0
22	0
23	0
24	0
25	0
26	0
27	0
28	1
29	3
30	5
31	9
32	16
33	27
34	42
35	65
36	97
37	140
38	196
39	266
40	349
41	444
42	548
43	655
44	759
45	852
46	927
47	977
48	999
49	990
50	950
51	884
52	798
53	697
54	591
55	485
56	386
57	298
58	222
59	161
60	113
61	77
62	51
63	32
64	20
65	12
66	7
67	3
68	2
69	1
70	0
71	0
72	0
73	0
74	0
75	0
76	0
77	0
78	0
79	0
80	0
81	0
82	0
83	0
84	0
85	0
86	0
87	0
88	0
89	0
90	0
91	0
92	0
93	0
94	0
95	0
96	0
97	0
98	0
99	0
100	0
>>END_MODULE
""", 'fastqc_output/SAMPLE_01_fastqc/summary.txt': """\
PASS	Basic Statistics	SAMPLE_01.fastq.gz
PASS	Per base sequence quality	SAMPLE_01.fastq.gz
WARN	Per sequence quality scores	SAMPLE_01.fastq.gz
PASS	Per base sequence content	SAMPLE_01.fastq.gz
PASS	Per sequence GC content	SAMPLE_01.fastq.gz
PASS	Per base N content	SAMPLE_01.fastq.gz
PASS	Sequence Length Distribution	SAMPLE_01.fastq.gz
FAIL	Sequence Duplication Levels	SAMPLE_01.fastq.gz
PASS	Overrepresented sequences	SAMPLE_01.fastq.gz
PASS	Adapter Content	SAMPLE_01.fastq.gz
""", 'fastqc_output/SAMPLE_02_fastqc/SAMPLE_02_fastqc.html': '<html><head><title>FastQC Report - SAMPLE_02</title></head><body><h1>FastQC Report</h1></body></html>', 'fastqc_output/SAMPLE_02_fastqc/fastqc_data.txt': """\
##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Filename	SAMPLE_02.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	980000
Sequences flagged as poor quality	0
Sequence length	150
%GC	50
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile
1	35	33	37
6	35	33	37
11	34	32	36
16	32	30	34
21	30	28	32
26	28	26	30
31	26	24	28
36	24	22	26
41	20	18	22
46	16	14	18
51	12	10	14
56	10	8	12
61	8	6	10
66	6	4	8
>>END_MODULE
>>Per sequence GC content	pass
#GC Content	Count
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
21	0
22	0
23	0
24	0
25	0
26	0
27	0
28	0
29	0
30	1
31	2
32	3
33	7
34	12
35	20
36	32
37	51
38	77
39	113
40	161
41	222
42	298
43	386
44	485
45	591
46	697
47	798
48	884
49	950
50	990
51	999
52	977
53	927
54	852
55	759
56	655
57	548
58	444
59	349
60	266
61	196
62	140
63	97
64	65
65	42
66	27
67	16
68	9
69	5
70	3
71	1
72	0
73	0
74	0
75	0
76	0
77	0
78	0
79	0
80	0
81	0
82	0
83	0
84	0
85	0
86	0
87	0
88	0
89	0
90	0
91	0
92	0
93	0
94	0
95	0
96	0
97	0
98	0
99	0
100	0
>>END_MODULE
""", 'fastqc_output/SAMPLE_02_fastqc/summary.txt': """\
PASS	Basic Statistics	SAMPLE_02.fastq.gz
PASS	Per base sequence quality	SAMPLE_02.fastq.gz
WARN	Per sequence quality scores	SAMPLE_02.fastq.gz
PASS	Per base sequence content	SAMPLE_02.fastq.gz
PASS	Per sequence GC content	SAMPLE_02.fastq.gz
PASS	Per base N content	SAMPLE_02.fastq.gz
PASS	Sequence Length Distribution	SAMPLE_02.fastq.gz
FAIL	Sequence Duplication Levels	SAMPLE_02.fastq.gz
PASS	Overrepresented sequences	SAMPLE_02.fastq.gz
PASS	Adapter Content	SAMPLE_02.fastq.gz
""", 'multiqc_output/multiqc_data/multiqc_data.json': """\
{
  "report_general_stats_data": {
    "fastqc": {
      "SAMPLE_01": {
        "percent_duplicates": 12.5,
        "percent_gc": 48.0,
        "total_sequences": 1200000.0
      },
      "SAMPLE_02": {
        "percent_duplicates": 18.2,
        "percent_gc": 51.0,
        "total_sequences": 980000.0
      }
    }
  }
}\
""", 'multiqc_output/multiqc_report.html': '<html><head><title>MultiQC Report</title></head><body><h1>MultiQC Report</h1></body></html>', 'report.md': """\
# BioQC Quality Control Report

**Date**: <TIMESTAMP> UTC  \n\
**Pipeline Run**: Demo Mode (Synthetic Data)  \n\

## Per-Sample Summary Metrics

| Sample | Total Sequences | GC Content (%) | Duplication Rate (%) | Status |
|--------|-----------------|----------------|----------------------|--------|
| SAMPLE_01 | 1,200,000 | 48.0% | 12.5% | ❌ FAIL (1 checks) |
| SAMPLE_02 | 980,000 | 51.0% | 18.2% | ❌ FAIL (1 checks) |

## QC Key Visualizations

### Position Quality Scores
![Base Quality](figures/position_quality_line.png)

### GC Distribution Profile
![GC Distribution](figures/gc_distribution.png)

## Outputs & Documentation

- `report.md` — Consolidated Markdown analysis summary
- `fastqc_output/` — Per-sample FastQC folders and reports
- `multiqc_output/multiqc_report.html` — Interactive aggregate MultiQC HTML
- `figures/` — High-resolution visualization files

## Reproducibility Bundle

- `reproducibility/commands.sh` — Portable script to replay this run
- `reproducibility/environment.yml` — Python/Conda dependencies
- `reproducibility/checksums.sha256` — Integrity digests

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "chat_summary_lines": [
    "\\u2705 BioQC Quality Control pipeline execution complete.",
    "\\ud83d\\udcca Output generated successfully inside output folder: <OUTPUT_DIR>",
    "\\ud83d\\udcc8 Rich custom visualizations rendered under figures/."
  ],
  "preferred_artifacts": [
    "<OUTPUT_DIR>/report.md",
    "<OUTPUT_DIR>/multiqc_output/multiqc_report.html"
  ],
  "report_md": "# BioQC Quality Control Report\\n\\n**Date**: <TIMESTAMP> UTC  \\n**Pipeline Run**: Demo Mode (Synthetic Data)  \\n\\n## Per-Sample Summary Metrics\\n\\n| Sample | Total Sequences | GC Content (%) | Duplication Rate (%) | Status |\\n|--------|-----------------|----------------|----------------------|--------|\\n| SAMPLE_01 | 1,200,000 | 48.0% | 12.5% | \\u274c FAIL (1 checks) |\\n| SAMPLE_02 | 980,000 | 51.0% | 18.2% | \\u274c FAIL (1 checks) |\\n\\n## QC Key Visualizations\\n\\n### Position Quality Scores\\n![Base Quality](figures/position_quality_line.png)\\n\\n### GC Distribution Profile\\n![GC Distribution](figures/gc_distribution.png)\\n\\n## Outputs & Documentation\\n\\n- `report.md` \\u2014 Consolidated Markdown analysis summary\\n- `fastqc_output/` \\u2014 Per-sample FastQC folders and reports\\n- `multiqc_output/multiqc_report.html` \\u2014 Interactive aggregate MultiQC HTML\\n- `figures/` \\u2014 High-resolution visualization files\\n\\n## Reproducibility Bundle\\n\\n- `reproducibility/commands.sh` \\u2014 Portable script to replay this run\\n- `reproducibility/environment.yml` \\u2014 Python/Conda dependencies\\n- `reproducibility/checksums.sha256` \\u2014 Integrity digests\\n\\n---\\n\\n*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*\\n"
}\
"""})
