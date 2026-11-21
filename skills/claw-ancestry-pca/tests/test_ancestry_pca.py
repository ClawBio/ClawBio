"""Tests for the claw-ancestry-pca skill."""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_SKILL_DIR = Path(__file__).resolve().parent.parent

# Add both project root (for clawbio.common) and skill dir (for ancestry_pca)
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
sys.path.insert(0, str(_SKILL_DIR))

from ancestry_pca import (
    compute_pca,
    generate_report,
    load_population_map,
    run_analysis,
    run_summary,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DEMO_VCF = _PROJECT_ROOT / "examples" / "demo_populations.vcf"
DEMO_POP_MAP = _PROJECT_ROOT / "examples" / "demo_population_map.csv"


@pytest.fixture
def small_geno_matrix():
    """Deterministic 10-sample, 20-variant genotype matrix."""
    rng = np.random.RandomState(42)
    mat = rng.choice([0, 1, 2], size=(10, 20))
    return mat


@pytest.fixture
def geno_with_missing():
    """Genotype matrix with some missing values (-1)."""
    rng = np.random.RandomState(42)
    mat = rng.choice([0, 1, 2], size=(10, 20))
    # Sprinkle missing
    mat[0, 0] = -1
    mat[3, 5] = -1
    mat[7, 15] = -1
    return mat


@pytest.fixture
def pop_map_csv(tmp_path):
    """Create a temporary population map CSV."""
    csv_path = tmp_path / "pop_map.csv"
    csv_path.write_text(
        "sample_id,population\n"
        "S1,AFR\nS2,AFR\nS3,EUR\nS4,EUR\nS5,EAS\n"
    )
    return csv_path


@pytest.fixture
def pop_map_tsv(tmp_path):
    """Create a temporary population map TSV."""
    tsv_path = tmp_path / "pop_map.tsv"
    tsv_path.write_text(
        "sample_id\tpopulation\n"
        "S1\tAFR\nS2\tAFR\nS3\tEUR\nS4\tEUR\nS5\tEAS\n"
    )
    return tsv_path


# ---------------------------------------------------------------------------
# Population map parsing
# ---------------------------------------------------------------------------


class TestLoadPopulationMap:
    def test_csv_parsing(self, pop_map_csv):
        samples = ["S1", "S2", "S3", "S4", "S5"]
        result = load_population_map(pop_map_csv, samples)
        assert result["S1"] == "AFR"
        assert result["S3"] == "EUR"
        assert result["S5"] == "EAS"
        assert len(result) == 5

    def test_tsv_parsing(self, pop_map_tsv):
        samples = ["S1", "S2", "S3", "S4", "S5"]
        result = load_population_map(pop_map_tsv, samples)
        assert result["S1"] == "AFR"
        assert result["S5"] == "EAS"

    def test_infer_from_prefix(self):
        samples = ["AFR_001", "EUR_002", "EAS_003", "XYZ_004"]
        result = load_population_map(None, samples)
        assert result["AFR_001"] == "AFR"
        assert result["EUR_002"] == "EUR"
        assert result["EAS_003"] == "EAS"
        assert result["XYZ_004"] == "XYZ"

    def test_missing_file_infers(self):
        samples = ["POP_A", "POP_B"]
        result = load_population_map(Path("/nonexistent.csv"), samples)
        # Falls back to prefix inference
        assert result["POP_A"] == "POP"
        assert result["POP_B"] == "POP"

    def test_flexible_column_names(self, tmp_path):
        csv_path = tmp_path / "flex.csv"
        csv_path.write_text(
            "iid,ancestry\nSAM1,AFR\nSAM2,EUR\n"
        )
        result = load_population_map(csv_path, ["SAM1", "SAM2"])
        assert result["SAM1"] == "AFR"
        assert result["SAM2"] == "EUR"


# ---------------------------------------------------------------------------
# PCA computation
# ---------------------------------------------------------------------------


class TestComputePCA:
    def test_basic_shape(self, small_geno_matrix):
        coords, var = compute_pca(small_geno_matrix, n_components=5)
        assert coords.shape == (10, 5)
        assert len(var) == 5

    def test_deterministic(self, small_geno_matrix):
        c1, v1 = compute_pca(small_geno_matrix, n_components=3)
        c2, v2 = compute_pca(small_geno_matrix, n_components=3)
        np.testing.assert_array_almost_equal(np.abs(c1), np.abs(c2))
        np.testing.assert_array_almost_equal(v1, v2)

    def test_variance_sums_to_one_or_less(self, small_geno_matrix):
        _, var = compute_pca(small_geno_matrix, n_components=9)
        assert var.sum() <= 1.0 + 1e-6

    def test_variance_decreasing(self, small_geno_matrix):
        _, var = compute_pca(small_geno_matrix, n_components=5)
        for i in range(len(var) - 1):
            assert var[i] >= var[i + 1] - 1e-10

    def test_missing_data_handling(self, geno_with_missing):
        coords, var = compute_pca(geno_with_missing, n_components=3)
        assert coords.shape[0] == 10
        assert not np.isnan(coords).any()

    def test_component_clamping(self, small_geno_matrix):
        # Request more components than samples-1
        coords, var = compute_pca(small_geno_matrix, n_components=50)
        assert coords.shape[1] <= small_geno_matrix.shape[0] - 1

    def test_all_missing_column(self):
        mat = np.array([[0, -1], [1, -1], [2, -1]])
        coords, var = compute_pca(mat, n_components=1)
        assert coords.shape == (3, 1)
        assert not np.isnan(coords).any()

    def test_single_component(self, small_geno_matrix):
        coords, var = compute_pca(small_geno_matrix, n_components=1)
        assert coords.shape[1] == 1
        assert len(var) == 1


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


class TestReportGeneration:
    def test_report_has_header(self, tmp_path):
        var = np.array([0.5, 0.2, 0.1])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=50,
            n_variants=500,
            n_components=3,
            explained_var=var,
            pop_counts={"AFR": 20, "EUR": 30},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "Ancestry Decomposition PCA" in report
        assert "claw-ancestry-pca" in report

    def test_report_has_disclaimer(self, tmp_path):
        var = np.array([0.5, 0.2])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=2,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "not a medical device" in report

    def test_report_has_variance_table(self, tmp_path):
        var = np.array([0.4, 0.3, 0.1])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=3,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "PC1" in report
        assert "Variance Explained" in report

    def test_report_has_population_table(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"AFR": 5, "EAS": 3, "EUR": 2},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "AFR" in report
        assert "EAS" in report
        assert "EUR" in report

    def test_report_figure_section_when_generated(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=True,
        )
        assert "pca_composite.png" in report

    def test_report_references(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"EUR": 10},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "Mallick" in report
        assert "ClawBio" in report


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


class TestCSVOutput:
    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_coordinates_csv_written(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        csv_path = tmp_path / "tables" / "pc_coordinates.csv"
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "sample_id" in df.columns
        assert "population" in df.columns
        assert "PC1" in df.columns

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_variance_csv_written(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        csv_path = tmp_path / "tables" / "variance_explained.csv"
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "component" in df.columns
        assert "variance_explained" in df.columns
        assert "cumulative" in df.columns


# ---------------------------------------------------------------------------
# Demo mode (end-to-end)
# ---------------------------------------------------------------------------


class TestDemoMode:
    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_full_report(self, tmp_path):
        result = run_analysis(
            DEMO_VCF, DEMO_POP_MAP, tmp_path,
            n_components=5, no_figures=True,
        )
        assert result["n_samples"] == 50
        assert result["n_variants"] > 0
        assert result["n_components"] <= 5
        assert len(result["explained_var"]) > 0
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_result_json(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        result_json = json.loads((tmp_path / "result.json").read_text())
        assert result_json["skill"] == "claw-ancestry-pca"
        assert "summary" in result_json
        assert "n_samples" in result_json["summary"]

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_summary_mode(self):
        text = run_summary(DEMO_VCF, DEMO_POP_MAP, n_components=3)
        assert "ANCESTRY DECOMPOSITION PCA" in text
        assert "PC1" in text
        assert "Samples: 50" in text

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_populations_detected(self, tmp_path):
        result = run_analysis(
            DEMO_VCF, DEMO_POP_MAP, tmp_path,
            n_components=3, no_figures=True,
        )
        pops = set(result["pop_counts"].keys())
        assert "AFR" in pops
        assert "EUR" in pops


# ---------------------------------------------------------------------------
# Reproducibility bundle
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not DEMO_VCF.exists(), reason="Demo VCF not available")
def test_demo_writes_verifiable_repro_bundle(tmp_path):
    import hashlib
    import subprocess

    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(_SKILL_DIR / "ancestry_pca.py"),
         "--demo", "--output", str(out), "--no-figures"],
        check=True, cwd=_PROJECT_ROOT, capture_output=True,
    )
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).exists(), name

    cmd = (repro / "commands.sh").read_text()
    assert "python skills/claw-ancestry-pca/ancestry_pca.py" in cmd
    assert "--demo" in cmd
    assert str(_PROJECT_ROOT) not in cmd

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "reproducibility/environment.yml",
            "tables/pc_coordinates.csv"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest, label


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/claw-ancestry-pca/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("claw-ancestry-pca", tmp_path) == snapshot({'report.md': """\
# Ancestry Decomposition PCA

**Date**: <TIMESTAMP> UTC
**Skill**: claw-ancestry-pca
**Version**: \n\
**Samples**: 50
**Populations**: 5
**Variants analysed**: 500
**Components computed**: 10
**Input files**:
- `demo_populations.vcf`: `f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68`
- `demo_population_map.csv`: `38e694611ada17aa7d4d48013b7e35e3eec294a80192dff5455514a665080152`

---
## Variance Explained

| Component | Variance | Cumulative |
|-----------|----------|------------|
| PC1 | 7.61% | 7.61% |
| PC2 | 4.87% | 12.48% |
| PC3 | 3.82% | 16.30% |
| PC4 | 3.70% | 20.00% |
| PC5 | 3.00% | 23.00% |
| PC6 | 2.97% | 25.97% |
| PC7 | 2.80% | 28.77% |
| PC8 | 2.78% | 31.56% |
| PC9 | 2.62% | 34.17% |
| PC10 | 2.59% | 36.76% |

## Population Composition

| Population | Count | Proportion |
|------------|-------|------------|
| AFR | 8 | 16.0% |
| AMR | 5 | 10.0% |
| EAS | 7 | 14.0% |
| EUR | 22 | 44.0% |
| SAS | 8 | 16.0% |

## PCA Plots

![PCA Composite](figures/pca_composite.png)

**Panel A** (PC1 vs PC2): Main axis of population structure.
**Panel B** (PC2 vs PC3): Finer population differentiation.
**Panel C** (PC1 vs PC3): Third axis of variation.
**Panel D** (Scree plot): Variance explained per component.

## Interpretation

- **PC1** typically captures the largest axis of global differentiation
- **PC2** separates major continental groups
- **PC3** often reveals finer substructure within continental groups
- Examine the scree plot to assess how many PCs carry meaningful signal

## Methods

- **Tool**: ClawBio Ancestry PCA v0.1.0
- **PCA**: scikit-learn PCA on mean-imputed genotype matrix (0/1/2 encoding)
- **Missing data**: Per-variant mean imputation
- **Input format**: VCF with biallelic SNPs

## References

- Mallick, S. et al. (2016). The Simons Genome Diversity Project. Nature, 538, 201-206.
- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio

## Reproducibility

```bash
python skills/claw-ancestry-pca/ancestry_pca.py \\
    --input demo_populations.vcf --pop-map demo_population_map.csv --output demo_output
```

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "claw-ancestry-pca",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:f847bb32a384508cd062bda6caac0979297dbef3fe14d7e645c2b21f1f9f7c68",
  "datasets": {},
  "summary": {
    "n_samples": 50,
    "n_variants": 500,
    "n_populations": 5,
    "n_components": 10,
    "pc1_variance": 7.61,
    "pc2_variance": 4.87
  },
  "data": {
    "population_counts": {
      "AFR": 8,
      "AMR": 5,
      "EAS": 7,
      "EUR": 22,
      "SAS": 8
    },
    "variance_explained": [
      0.076102,
      0.048713,
      0.038195,
      0.036982,
      0.030041,
      0.029709,
      0.027982,
      0.02784,
      0.026182,
      0.025866
    ]
  }
}\
""", 'tables/pc_coordinates.csv': """\
sample_id,population,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
AFR_001,AFR,6.588584104475454,-2.4032581363144176,-2.2068474793289954,3.1866098376441974,2.0454337657954884,3.8635128801405574,1.7034123070236582,1.5517224073007785,4.697853495704004,2.1259158541491447
AFR_002,AFR,5.500074576751851,-1.3239001391177558,-1.7152282564758632,2.7491666501223153,4.405876490176529,0.20797273103427333,-5.3104849585621245,0.41620429371337103,-2.2677662986157516,-2.1748619481679365
AFR_003,AFR,6.622481881926115,-2.219365140592427,1.0863013276629039,4.285988572446751,1.69190063931399,-2.520552274871412,0.11362742048501404,1.5530410147476235,0.23723659218679,2.8373207906738243
AFR_004,AFR,6.752906044066396,-1.8842529026646946,1.8568062452644551,3.521393640328821,-0.10404068340067744,-2.440706074757615,-1.5920131310673347,5.010883538596307,4.027657861212642,1.5265766553801992
AFR_005,AFR,7.064368213726153,-2.748596396396258,-0.5696956941090658,1.1900391874062488,-0.7376323162481185,1.3424176388898015,3.0321010398385457,-0.688602750992322,-3.4543793327827386,-1.108611473969081
AFR_006,AFR,7.361290987534113,-1.1879065700063174,1.2413814422478051,3.271183924017567,-2.8537720676637237,5.4413404113079284,-1.1719854506230654,-2.338518623262053,-3.3737671466428045,-0.1542551996120561
AFR_007,AFR,5.509413302593154,-2.0355283632179066,1.3288569320440924,3.1395412257950874,-0.7998620374068854,-2.956476227942584,5.675543484462023,-3.5805162431123234,1.6521397509822409,1.3551543226662495
AFR_008,AFR,5.304182542963415,-2.1721167946900093,1.4507895287394856,3.7970737313434624,-3.0152166554071114,-2.626436035814368,-2.511177739610372,-2.4526074551567367,-0.5799524447829263,-4.47703435050791
AMR_001,AMR,1.9922999399294583,2.4664906132627404,6.0623506638358915,-3.7681890334846924,-1.569228349662947,1.6480843977425483,-0.4251030779476245,1.7230575197136018,0.06829281867062546,4.215350249595932
AMR_002,AMR,2.0250481395907407,0.6898415631979713,6.256563557697707,-4.928425016811298,-1.0915042882186494,-0.7226825900841357,0.7849425395879771,1.2352505481419789,0.4046240633822593,-1.5452120277365766
AMR_003,AMR,2.0157448270476315,0.8010006022901897,5.259737902186945,-2.3921073556019414,5.652125215193494,4.256418146577487,-1.7300492702937655,-1.3579882152377405,3.926843133234107,-1.2994890829409065
AMR_004,AMR,0.4191880014742338,1.2777905896843,6.150447939220487,-3.863328646288267,2.031367476181642,-0.9656896947750097,-3.1784010652307537,1.2854731321863797,-0.8993596737178122,1.3460226988235926
AMR_005,AMR,2.0834824581871527,0.6432646812366181,4.5460377119941855,-3.046594395727094,-3.977197878475847,1.406943092869682,2.3931520796229404,-1.980888538942484,0.6388995488401232,-1.7684750550243973
EAS_001,EAS,-2.4650626070106867,5.825048734806144,0.9409384211685249,1.8415770670711495,0.46725592709954594,2.8301046230270224,-3.2103462256150217,0.7829193967443473,-0.7464825706265793,-0.6174407004279078
EAS_002,EAS,-2.012843148762523,5.716532023273461,2.323689469329714,2.7982063076361174,-0.4472146917729554,3.38902958023497,3.1436853772315776,-0.17440992148980422,-1.4374771718318031,-1.5846370833657113
EAS_003,EAS,-2.3194837983285406,6.248314783484846,0.6119450685929405,2.270882403683786,0.5678170273622014,-2.367234311644643,-0.010756610962734184,-0.10206070867756942,0.4955678689341346,-1.669259734533774
EAS_004,EAS,-1.6152145289936413,4.947917500192079,-0.5033790537962072,1.7751646530957503,1.2180649810509971,-0.5963665365874592,-0.5121644176032389,0.6286919051806614,-0.690474100718058,-1.610622407411548
EAS_005,EAS,-1.9619495547207135,5.457645797757849,-1.5318756279626158,3.649364956655847,-1.8994176342998224,-0.00907161233984969,1.4903552091017946,3.1810402296860105,-0.5160426188050643,1.0463882498247872
EAS_006,EAS,-2.1784452403338235,6.958535652072063,1.2639434944024976,3.547378109916199,-0.05361018865300107,-1.3143214027821417,0.8920972441851238,-0.15289021911816292,-1.7462789607603655,-0.3175330055734937
EAS_007,EAS,-3.200352230063597,6.05914813514275,-0.24369122872562882,2.948966990691313,2.098329899387515,-3.041452196703103,2.5206632683696175,-0.8734147833598541,2.7985529524983392,2.5554004345433365
EUR_001,EUR,-3.005019793391407,-2.385174092298822,-2.0022240625172945,0.34516419918973307,2.8747801509866604,2.3448535582002807,0.43178515287109875,-1.2470268960530286,0.1605997783804465,2.9662523870142916
EUR_002,EUR,-2.9508428961943176,-3.3226599934056167,0.5324082307118398,-0.5666602475014101,-1.8427856536056602,-0.2117539527396795,1.0733407611256338,-1.5934649932351934,2.485800879532795,-1.4668660191917033
EUR_003,EUR,-3.662559393455796,-3.3674306828791583,0.9858517924418586,0.8477329865731218,-1.0969963954023643,2.6663603075643745,-0.9446763881369973,-0.040109882793275146,-1.2056614533007213,-2.9273956515320867
EUR_004,EUR,-2.6206419533121723,-0.8839423090305017,1.5176582326266466,0.5723342753885596,0.2608301860924485,1.3867423904970837,-1.8672189463843352,-4.973894275751252,0.23279753047493876,0.8238097581044438
EUR_005,EUR,-4.171600699158689,-2.4616450145777384,-1.6933434884773488,-1.1197056549293,1.4896044438001144,-0.5829468590383821,-3.1130099984366266,0.15050228388509204,3.872263599882872,1.7225611297959944
EUR_006,EUR,-2.5074912854617057,-1.3640669899869442,-0.6088850728115721,0.7444004037245336,-3.519759246405483,1.0835774637416429,1.5920459971767673,2.517094908485661,0.6968555783927926,0.5366232570485254
EUR_007,EUR,-3.3626134828898553,-1.7015491315460487,1.059421143189859,-1.6071999381368651,-2.2761094160384685,-4.608250317202594,-0.8017300987407052,1.7756995206339135,1.9613298598754765,-2.9678175849209674
EUR_008,EUR,-4.132714929769827,-1.702109606346508,-1.0239794802248272,2.2472115688022347,-0.6883587072228878,-0.009239952940869973,1.6257205938222747,1.676526396424764,0.561012721695785,0.5723117242511886
EUR_009,EUR,-2.29320622235729,-2.225279720821073,-1.913068087010586,-0.12478930716059103,0.9559722961568691,1.853744232517353,1.3949686513192268,3.770632389987607,0.060356309741397426,-3.3387404413861956
EUR_010,EUR,-2.2601118780396185,-2.543803307020824,0.30320039918932595,1.0036966478917837,-1.558403241875427,-1.130663326985731,-2.4682151380716486,-3.6299324904193004,-1.4905833363269043,3.9108305485369086
EUR_011,EUR,-4.167693640552878,-1.6917502072574144,-0.2555230458541588,-0.8016387337285942,1.2633428359265497,1.7832584646231446,0.31898403049524854,-1.329744930688891,-3.19204395044521,0.22355104789845492
EUR_012,EUR,-3.0274329728100065,-0.044896682112787266,-0.5313592363242923,0.33779341821860065,-5.346005808116709,0.7805093860818844,-4.647976463313365,0.7969496524832422,2.8100230296577102,1.725861772995509
EUR_013,EUR,-2.9875903352605375,-2.976681314621191,-1.57524218103799,0.5189689359156046,-1.655279152001147,-1.5834758579642043,-0.5366366531501585,-1.741182589475753,2.563469032308726,-1.1976192787355489
EUR_014,EUR,-3.768635938048815,-4.008616014068064,0.5262553884370307,-0.9852596016242556,2.1070442849458026,-0.9285063450083825,4.182443885824196,4.330383147465063,-1.6559600300706594,-1.1845434513667799
EUR_015,EUR,-2.050461147412664,-0.46990569097860446,-0.0190350253352877,0.5699860315011933,2.01540542728575,0.5202439132494335,-0.52783592650727,-1.4828714497725248,-0.1571103450397415,-2.728831145875962
EUR_016,EUR,-4.488670978975385,-1.0300412301843107,-1.2833272966017844,-0.649493681771606,-0.48832770573113027,0.08770471741485941,0.185033253422366,-4.978413911280321,-0.5011118318683919,3.622479306680456
EUR_017,EUR,-2.7462219392758196,-3.17853622141385,1.9175094838199191,-0.4657175398603004,1.7348473414006438,-2.4550950594499388,1.9391384680807544,-0.5941663418209022,-2.3624574860356407,-0.07728101605274396
EUR_018,EUR,-3.4181475318635552,-1.121493881761713,-1.9193486433277203,-0.10524020498555858,-1.7467971616301334,-0.25906656849360593,-3.9071203194265007,0.09426523450933495,-0.22987064622866757,1.3690718743119172
EUR_019,EUR,-2.564824675350273,-2.156469281876718,-2.1414498556129984,-2.45047876577546,-0.956291093607508,4.905205445604439,3.515026802633921,1.1669873305556604,-0.5498883057110354,-0.2174865314771929
EUR_020,EUR,-2.6237866568184267,-2.8655666794461196,-1.0678847446536919,1.8346374423555523,1.976544654654581,0.31047004097930925,-2.9412913370582943,4.35577404801082,-3.139435255806892,-0.42810620859583465
EUR_021,EUR,-2.4566442452352475,-1.9968760401118864,2.123086676251704,-0.4074857304820548,2.1843918448292996,0.6131075028548423,1.5553041033509234,1.0910581639712704,1.0025900952049125,1.292281119380368
EUR_022,EUR,-3.2194015411023065,-1.9058672010898403,-0.002924607852025605,1.1774382148884814,1.3580522542443163,-4.225656769519704,1.4363037767721998,-2.025918342611011,-1.9381424780893302,-1.124460555581554
SAS_001,SAS,2.8578186371108254,1.2351303356942944,-4.708140183597929,-4.972077174767886,-0.8806250021576648,0.005731216770895178,0.2175591207935577,2.0246169276842316,1.6974141756319634,-3.0151846979063506
SAS_002,SAS,3.507510275183635,0.11234609764862559,0.19636177390118398,-4.747257638662236,2.039174423563891,-4.56226246645032,-0.1717508243603705,-0.7202569484017027,-1.5221346763500665,-1.6054531229098221
SAS_003,SAS,2.6923698332260453,2.104587889489548,-3.4499606925227493,-3.4272834196174413,1.6262896329198084,-0.0657831291465016,0.01407323762164448,1.510254963076357,-5.205477079923621,3.658888306738136
SAS_004,SAS,2.015706479437921,2.7110994887559103,-3.0502127938155024,-1.950779052329317,-4.428722846753163,1.1318273023884993,0.4100737275645082,-0.0648926319703702,1.3268125945715254,-0.8743657684536601
SAS_005,SAS,3.9826334062726896,0.309241608298859,-1.6205977294089486,-2.930504442397768,-3.1605064898482937,-3.4326571616078905,-0.6389467060669752,1.2215400462949262,-3.9392450425011845,2.9834910075001293
SAS_006,SAS,3.119435085319239,3.157071684696992,-5.073430511974161,-1.9752730218328045,3.257579187832211,-0.4552034720222489,-0.4927911413610592,-3.8485070250195856,1.6228137967897986,-3.8741977466440214
SAS_007,SAS,2.7448187993851185,1.6826053261226213,-3.294758815107916,-3.8477528836677393,3.2321412551190827,-0.09619586765693448,2.022851785138156,-1.6265347179868441,1.8269892018923846,2.591333889777265
SAS_008,SAS,4.080307708748778,2.975672628727656,-5.536129930489844,-3.0386598951595176,-2.360506929713652,0.3085866202169739,-0.9525514293904145,-0.25175411284999194,0.9723059673031795,0.3483049042110572
""", 'tables/variance_explained.csv': """\
component,variance_explained,cumulative
PC1,0.07610241021794688,0.07610241021794688
PC2,0.04871314960803677,0.12481555982598365
PC3,0.03819483049928456,0.1630103903252682
PC4,0.036982422460285026,0.19999281278555323
PC5,0.0300408445637548,0.23003365734930803
PC6,0.029708674636452843,0.25974233198576085
PC7,0.027981539713915957,0.28772387169967684
PC8,0.02783990168429496,0.3155637733839718
PC9,0.026182144903722535,0.34174591828769435
PC10,0.025865955517337214,0.3676118738050316
"""})
