"""Tests for the rare-disease-rnaseq skill."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rare_disease_rnaseq import (
    DEFAULT_PANEL,
    call_outliers,
    cpm_log2,
    generate_demo,
    per_gene_robust_stats,
)


def test_default_panel_loads_and_has_required_columns():
    panel = pd.read_csv(DEFAULT_PANEL)
    assert {"gene", "mechanism", "phenotype"}.issubset(panel.columns)
    assert len(panel) >= 30
    assert "FBN1" in panel["gene"].values
    assert "NF1" in panel["gene"].values


def test_cpm_log2_normalises_library_size():
    counts = pd.DataFrame({"S1": [10, 20, 30], "S2": [100, 200, 300]},
                          index=["G1", "G2", "G3"])
    expr = cpm_log2(counts)
    # After CPM normalisation S1 and S2 should yield identical per-gene values
    np.testing.assert_allclose(expr["S1"].values, expr["S2"].values)


def test_per_gene_robust_stats_zero_mad_for_constant_gene():
    counts = pd.DataFrame(
        {f"C{i}": [100, 200] for i in range(20)},
        index=["constant_gene", "constant_gene_2"],
    )
    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, list(counts.columns))
    assert (stats["mad"] == 0).all()


def test_zero_mad_genes_skipped_in_outlier_call():
    # Construct expr + stats directly so library-size effects do not
    # confound the unit. constant_gene has MAD=0 (must be skipped);
    # noisy_gene has MAD>0 and the case sits 10 SDs above (must be flagged).
    samples = [f"C{i}" for i in range(20)] + ["CASE-1"]
    expr = pd.DataFrame(
        {
            "constant_gene": [10.0] * 21,
            "noisy_gene": [10.0] * 20 + [15.0],
        },
        index=samples,
    ).T
    stats = pd.DataFrame({
        "median": {"constant_gene": 10.0, "noisy_gene": 10.0},
        "mad": {"constant_gene": 0.0, "noisy_gene": 0.3},
    })
    panel = pd.DataFrame([
        {"gene": "constant_gene", "mechanism": "test", "phenotype": "test"},
        {"gene": "noisy_gene", "mechanism": "test", "phenotype": "test"},
    ])
    outliers = call_outliers(expr, ["CASE-1"], stats, panel, z_threshold=3.0)
    assert "constant_gene" not in outliers["gene"].values
    assert "noisy_gene" in outliers["gene"].values
    assert (outliers[outliers["gene"] == "noisy_gene"]["z_score"] > 0).all()


def test_demo_recovers_injected_outliers(tmp_path):
    counts_path, cases_path, controls_path = generate_demo(seed=42, output_dir=tmp_path)
    counts = pd.read_csv(counts_path, index_col=0)
    case_ids = [s.strip() for s in cases_path.read_text().splitlines() if s.strip()]
    control_ids = [s.strip() for s in controls_path.read_text().splitlines() if s.strip()]

    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, control_ids)
    panel = pd.read_csv(DEFAULT_PANEL)
    outliers = call_outliers(expr, case_ids, stats, panel, z_threshold=3.0)

    # FBN1 must be a strong DOWN outlier in CASE-001
    fbn1 = outliers[(outliers["case"] == "CASE-001") & (outliers["gene"] == "FBN1")]
    assert len(fbn1) == 1
    assert fbn1.iloc[0]["direction"] == "down"
    assert fbn1.iloc[0]["z_score"] < -5

    # NF1 must be a strong UP outlier in CASE-002
    nf1 = outliers[(outliers["case"] == "CASE-002") & (outliers["gene"] == "NF1")]
    assert len(nf1) == 1
    assert nf1.iloc[0]["direction"] == "up"
    assert nf1.iloc[0]["z_score"] > 5

    # Both flagged outliers should be in the disease panel
    panel_hits = outliers[outliers["in_panel"]]
    assert "FBN1" in panel_hits["gene"].values
    assert "NF1" in panel_hits["gene"].values


def test_threshold_controls_call_count(tmp_path):
    counts_path, cases_path, controls_path = generate_demo(seed=7, output_dir=tmp_path)
    counts = pd.read_csv(counts_path, index_col=0)
    case_ids = [s.strip() for s in cases_path.read_text().splitlines() if s.strip()]
    control_ids = [s.strip() for s in controls_path.read_text().splitlines() if s.strip()]
    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, control_ids)
    panel = pd.read_csv(DEFAULT_PANEL)
    loose = call_outliers(expr, case_ids, stats, panel, z_threshold=2.0)
    strict = call_outliers(expr, case_ids, stats, panel, z_threshold=5.0)
    assert len(strict) <= len(loose)


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `uv run pytest skills/rare-disease-rnaseq/tests --inline-snapshot=fix`; it
# rewrites the snapshot() values below, and the PR diff shows which numbers moved.

from inline_snapshot import snapshot  # noqa: E402

from clawbio.common import demo_runner  # noqa: E402


@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("rare-disease-rnaseq", tmp_path) == snapshot({'demo_data/demo_cases.txt': """\
CASE-001
CASE-002
""", 'demo_data/demo_controls.txt': """\
CTRL-000
CTRL-001
CTRL-002
CTRL-003
CTRL-004
CTRL-005
CTRL-006
CTRL-007
CTRL-008
CTRL-009
CTRL-010
CTRL-011
CTRL-012
CTRL-013
CTRL-014
CTRL-015
CTRL-016
CTRL-017
CTRL-018
CTRL-019
CTRL-020
CTRL-021
CTRL-022
CTRL-023
CTRL-024
CTRL-025
CTRL-026
CTRL-027
CTRL-028
CTRL-029
CTRL-030
CTRL-031
CTRL-032
CTRL-033
CTRL-034
CTRL-035
CTRL-036
CTRL-037
CTRL-038
CTRL-039
CTRL-040
CTRL-041
CTRL-042
CTRL-043
CTRL-044
CTRL-045
CTRL-046
CTRL-047
CTRL-048
CTRL-049
CTRL-050
CTRL-051
CTRL-052
CTRL-053
CTRL-054
CTRL-055
CTRL-056
CTRL-057
CTRL-058
CTRL-059
CTRL-060
CTRL-061
CTRL-062
CTRL-063
CTRL-064
CTRL-065
CTRL-066
CTRL-067
CTRL-068
CTRL-069
CTRL-070
CTRL-071
CTRL-072
CTRL-073
CTRL-074
CTRL-075
CTRL-076
CTRL-077
CTRL-078
CTRL-079
CTRL-080
CTRL-081
CTRL-082
CTRL-083
CTRL-084
CTRL-085
CTRL-086
CTRL-087
CTRL-088
CTRL-089
CTRL-090
CTRL-091
CTRL-092
CTRL-093
CTRL-094
CTRL-095
CTRL-096
CTRL-097
CTRL-098
CTRL-099
""", 'report.md': """\
# Rare-Disease Blood RNA-seq Outlier Report

Generated: <TIMESTAMP>
Skill: `rare-disease-rnaseq` v0.1.0

## Cohort summary

- Cases analysed: **2**
- Control reference panel: **100** samples
- Genes screened: **200**
- Outlier threshold: **|z| >= 3.0**
- Outlier events called: **5**

## Per-case candidate diagnoses

### CASE-001

- Total outlier genes: 1
- In disease panel: 1
- Off-panel (informational): 0

**Disease-panel outliers (clinical priority):**

| Gene | Direction | z-score | Mechanism | Associated phenotype |
|------|-----------|---------|-----------|----------------------|
| FBN1 | down | -18.62 | haploinsufficiency | Marfan syndrome |

**Suggested next step:** correlate flagged outlier(s) with phenotype, review WGS variant calls in the same gene, present at MDT.

### CASE-002

- Total outlier genes: 4
- In disease panel: 2
- Off-panel (informational): 2

**Disease-panel outliers (clinical priority):**

| Gene | Direction | z-score | Mechanism | Associated phenotype |
|------|-----------|---------|-----------|----------------------|
| NF1 | up | 9.03 | haploinsufficiency | Neurofibromatosis type 1 |
| SMC1A | down | -3.66 | haploinsufficiency | Cornelia de Lange syndrome |

**Suggested next step:** correlate flagged outlier(s) with phenotype, review WGS variant calls in the same gene, present at MDT.

## Method

Per-gene robust outlier scoring on log2(CPM+1). For each gene, modified z-score = 0.6745 * (x - median) / MAD computed on the control reference panel. Outliers are events with |z| above threshold; clinical priority is given to genes in the dosage-sensitive disease panel.

**Production swap:** clinical-grade calls use the DROP pipeline (gagneurlab/drop) with OUTRIDER (denoising autoencoder, NB-distributed) and FRASER2 (splicing outliers). The skill's I/O contract is unchanged.

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "rare-disease-rnaseq",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:1f779df56c7bf38a7f5660b9f4543b29fb4fb84bf662980f5831420ed243791c",
  "datasets": {},
  "summary": {
    "cases": 2,
    "controls": 100,
    "genes_screened": 200,
    "z_threshold": 3.0,
    "outlier_events": 5,
    "panel_outlier_events": 3,
    "cases_with_panel_hit": 2
  },
  "data": {
    "outliers": [
      {
        "case": "CASE-001",
        "gene": "FBN1",
        "z_score": -18.62,
        "direction": "down",
        "log2_cpm_case": 8.74,
        "log2_cpm_control_median": 12.15,
        "in_panel": true,
        "mechanism": "haploinsufficiency",
        "phenotype": "Marfan syndrome"
      },
      {
        "case": "CASE-002",
        "gene": "NF1",
        "z_score": 9.03,
        "direction": "up",
        "log2_cpm_case": 11.95,
        "log2_cpm_control_median": 10.2,
        "in_panel": true,
        "mechanism": "haploinsufficiency",
        "phenotype": "Neurofibromatosis type 1"
      },
      {
        "case": "CASE-002",
        "gene": "SMC1A",
        "z_score": -3.66,
        "direction": "down",
        "log2_cpm_case": 10.58,
        "log2_cpm_control_median": 11.15,
        "in_panel": true,
        "mechanism": "haploinsufficiency",
        "phenotype": "Cornelia de Lange syndrome"
      },
      {
        "case": "CASE-002",
        "gene": "BG0003",
        "z_score": -4.61,
        "direction": "down",
        "log2_cpm_case": 10.72,
        "log2_cpm_control_median": 11.26,
        "in_panel": false,
        "mechanism": "",
        "phenotype": ""
      },
      {
        "case": "CASE-002",
        "gene": "BG0122",
        "z_score": -3.21,
        "direction": "down",
        "log2_cpm_case": 10.61,
        "log2_cpm_control_median": 11.09,
        "in_panel": false,
        "mechanism": "",
        "phenotype": ""
      }
    ]
  }
}\
""", 'tables/outlier_calls.csv': """\
case,gene,z_score,direction,log2_cpm_case,log2_cpm_control_median,in_panel,mechanism,phenotype
CASE-001,FBN1,-18.62,down,8.74,12.15,True,haploinsufficiency,Marfan syndrome
CASE-002,NF1,9.03,up,11.95,10.2,True,haploinsufficiency,Neurofibromatosis type 1
CASE-002,SMC1A,-3.66,down,10.58,11.15,True,haploinsufficiency,Cornelia de Lange syndrome
CASE-002,BG0003,-4.61,down,10.72,11.26,False,,
CASE-002,BG0122,-3.21,down,10.61,11.09,False,,
""", 'tables/per_gene_stats.csv': """\
gene,median,mad
FBN1,12.147585365311105,0.12334711741188276
NF1,10.197374597540751,0.13070609016535517
TP53,12.817802082896652,0.090427115575892
PTEN,13.120848941992264,0.08159260891996478
BRCA1,8.922051028675543,0.11283010620032563
BRCA2,9.881746182711924,0.12882647996580854
MLH1,11.901794564659852,0.11163942260711934
MSH2,11.286208526593075,0.14039538212135572
COL1A1,11.681038605560996,0.11279987133619507
COL2A1,10.499899758186494,0.09890695075073452
COL3A1,12.98761165220716,0.09692149838467845
COL5A1,12.84974797681226,0.10264873277033892
MYH7,11.828938580691336,0.09907785067596286
LMNA,13.349580910000576,0.09110410488441811
RYR2,12.418809955474144,0.11864290713473924
KCNQ1,10.459441829345977,0.11375227959980805
RB1,12.25518725018616,0.11874491120248098
APC,10.338494287970985,0.13424444042400818
CDH1,13.00663276822652,0.08104990387553901
VHL,11.614190623259168,0.0958781061566576
NF2,11.463157799957553,0.08741100791355372
RUNX1,10.742031065722458,0.11691391074234936
SETD2,13.49187094484954,0.10426278159864744
EP300,11.501840922690057,0.10114149060654753
CBFB,11.108425220707332,0.10693609414643301
DDX3X,11.205768186620523,0.10465229773757123
EHMT1,12.531979261866088,0.09078613766945232
KMT2A,12.230758224750522,0.09251879935979268
KMT2D,12.371692977757398,0.09013712461031886
MED12,12.338402906407634,0.10041398842935045
NSD1,14.813489231004286,0.08528731312235127
SMC1A,11.14709636073676,0.10404520304769527
SOX2,10.978306069724844,0.09552751835866946
SOX9,10.547219480278603,0.12409763471289992
SOX10,12.612544284092598,0.0981745926244173
ARID1B,13.35742880737951,0.0888453704656813
CHD7,11.569337395127437,0.11866208152893964
CHD8,10.486355680834192,0.11362718099811708
CREBBP,10.53772109897944,0.1080866561723397
FOXG1,12.670537140814046,0.09992122245516111
FOXP2,12.82025699402598,0.12387967300956859
GLI3,12.526529281357234,0.09875418339856701
KAT6A,10.727731949781152,0.12752853878880188
NIPBL,12.054422419549457,0.09547474323867888
PAX3,11.904543099582146,0.08458580914803537
SATB2,12.048833983767185,0.09268229766906977
SHH,13.008065855002807,0.10163409652672684
SMARCA4,12.071810082055498,0.0975350083837947
SMARCB1,12.697648163732858,0.1148570531080857
TCF4,11.814676041859867,0.10258507461361255
BG0000,12.100980836759515,0.10592803706643483
BG0001,12.664567984591313,0.08961249606476862
BG0002,9.654464993619861,0.11176469884344797
BG0003,11.25886031021167,0.07897405208282393
BG0004,11.05948768605591,0.10171829798603049
BG0005,10.820374871872259,0.0992211038418862
BG0006,11.353635604208975,0.09062586534121753
BG0007,13.854310053569257,0.07935380918040469
BG0008,10.505401183926082,0.12402030226633265
BG0009,13.143781421441544,0.09899442741955244
BG0010,9.334026768956726,0.10619973325482057
BG0011,11.258682178280017,0.0978008031440547
BG0012,11.962204681675523,0.10346804278877375
BG0013,12.575257075631713,0.13055177056609324
BG0014,12.77907946524411,0.08054151063381987
BG0015,12.880906039508506,0.1025348728656601
BG0016,11.259775265239831,0.11470196344399408
BG0017,11.070748595113914,0.10803534436946904
BG0018,12.953803852018588,0.10569166325947066
BG0019,11.428555806030943,0.12070523194698346
BG0020,9.84666967109915,0.10607936980069077
BG0021,10.121996257102328,0.11434031031126057
BG0022,10.396916049734607,0.10425899709351061
BG0023,12.480567275695902,0.0846002731367097
BG0024,11.94999190948261,0.08978518335291064
BG0025,12.695122351408523,0.13892486968758355
BG0026,11.112329087736573,0.09831098502508517
BG0027,11.928600781095845,0.11001651295781478
BG0028,12.639828860642563,0.12201031351988068
BG0029,11.28580622923907,0.10018150462665787
BG0030,12.422620515762246,0.12378728400397065
BG0031,10.787593134632809,0.07818696490280086
BG0032,11.246692079011497,0.08393093281473263
BG0033,11.178139295350723,0.09123558515100871
BG0034,10.00440075053858,0.11585259350320332
BG0035,12.418348489544005,0.08882645881504647
BG0036,11.025693871388428,0.08434816468836104
BG0037,11.760639655178522,0.10068804664501751
BG0038,12.463996113573991,0.09911617575770837
BG0039,12.37256820899702,0.11911804189854536
BG0040,12.678238555594803,0.11220649059242138
BG0041,11.603659420980746,0.1135919617132144
BG0042,11.100304883380828,0.12065005592456757
BG0043,11.617080075862843,0.0851860351798015
BG0044,9.302482546086843,0.1268640807397574
BG0045,9.597885296559177,0.10070929041211851
BG0046,9.800123475519175,0.1384459182041473
BG0047,10.328228062245092,0.1256600492430211
BG0048,12.276143258703161,0.08718987997109195
BG0049,10.425590907841105,0.10956686796232518
BG0050,11.19850444403345,0.09262355141682743
BG0051,13.59780749355824,0.07727834330415462
BG0052,11.205511422154839,0.10328193883228565
BG0053,12.786371759993596,0.10710724846785613
BG0054,10.378929517097609,0.11514923550848799
BG0055,11.439145915045218,0.11264206890005735
BG0056,10.362689583547908,0.09391390883129613
BG0057,11.238081968571347,0.12030254192736578
BG0058,12.945372249508624,0.11752052642100885
BG0059,9.270019670702002,0.14210693956189413
BG0060,12.351260532364037,0.09184282081815098
BG0061,12.069869501180046,0.10452623899310254
BG0062,10.873485315976948,0.11554081575981634
BG0063,9.647977954857772,0.11509420244719326
BG0064,11.80432795217127,0.11165596479646833
BG0065,10.95440322543901,0.12558328966756704
BG0066,12.028945124508365,0.0933039259560351
BG0067,11.747616109780719,0.10694420387174386
BG0068,14.048009432409057,0.11075123667392361
BG0069,11.385439783060207,0.09876695857742934
BG0070,10.275634223149748,0.10754731768580772
BG0071,11.966393194494808,0.10529895060759653
BG0072,12.092375626615368,0.10535109087998507
BG0073,13.657690223562568,0.10895321462537222
BG0074,12.920242384581297,0.08867441191334269
BG0075,12.210761392695856,0.13736830544451184
BG0076,13.826086827855745,0.08118445177748779
BG0077,10.044127391876287,0.1172640487033485
BG0078,10.79878854271644,0.10522488277766406
BG0079,10.390413389859289,0.11769767437038059
BG0080,11.17329390751217,0.09826086698221648
BG0081,9.755983446771095,0.11956406245442786
BG0082,12.64027722712699,0.0855126320773758
BG0083,11.396322444879171,0.11114331826306856
BG0084,9.570221242788646,0.12923529006668577
BG0085,10.303076763567375,0.11206079188823814
BG0086,12.164357701294467,0.09809423934337413
BG0087,12.955750369279198,0.09600505392935865
BG0088,14.566291663947307,0.08355921846863712
BG0089,15.942602143783361,0.0951483385750187
BG0090,12.34361447314034,0.0983269596542673
BG0091,10.319252440655697,0.125037654281682
BG0092,8.645344206115208,0.1186645695997699
BG0093,12.101169148511186,0.10063011250900988
BG0094,10.5394414474708,0.08107876121018709
BG0095,11.144463771073626,0.12086263785914042
BG0096,10.851783194546343,0.1023150556120207
BG0097,11.51594601329888,0.12352589531533553
BG0098,13.228977520452542,0.11944684195978628
BG0099,11.948646754079396,0.09469532804618552
BG0100,11.497685877870634,0.11916901278544767
BG0101,10.211287016519565,0.12204098372340244
BG0102,9.291244598097556,0.12056011974607284
BG0103,11.002950977094539,0.10626707188460305
BG0104,11.643186473809305,0.10125900231971219
BG0105,14.291660307094581,0.09235015509431754
BG0106,11.954189937212647,0.10257220672489265
BG0107,13.154466876990062,0.12778549826396457
BG0108,11.011266275741928,0.09365910169900005
BG0109,9.994465899354456,0.10375415046142056
BG0110,10.311671994820959,0.1063178377831333
BG0111,10.65804082211821,0.11748105766793504
BG0112,14.788988761363278,0.08380279235024712
BG0113,10.573285607644745,0.1370936594434191
BG0114,12.945304182425218,0.10361859471040358
BG0115,10.406578261329837,0.1078176690125563
BG0116,13.094265059356442,0.09752177342674884
BG0117,12.328196913002948,0.08577226017647188
BG0118,11.497216257126684,0.12804384344268627
BG0119,11.64797687611048,0.12209005363308822
BG0120,10.782611044992194,0.12375189022663591
BG0121,12.368041873854544,0.09686369284955543
BG0122,11.088323709162218,0.10053523625171401
BG0123,9.968796307244329,0.12231364588627969
BG0124,9.865332893246725,0.10234579364873753
BG0125,11.995945336907038,0.0952784187636615
BG0126,13.996297768998241,0.10555394765656789
BG0127,11.972058549098268,0.1204545802726873
BG0128,11.558341442596095,0.09615913673868981
BG0129,12.16060914213973,0.09906271645963294
BG0130,13.64338296718772,0.09010833724678768
BG0131,12.012968462568153,0.0930796407916219
BG0132,11.080488883198232,0.11596702319218721
BG0133,13.323707843921479,0.10289970475342525
BG0134,12.347087204797276,0.11846693721348789
BG0135,13.974554302812013,0.08994977099271928
BG0136,11.997245674419231,0.11643488933961432
BG0137,9.965085363656474,0.14234619903954115
BG0138,9.735671014539767,0.10899936404543009
BG0139,14.093352491095386,0.10625143267764159
BG0140,14.204080451467613,0.08921435461444105
BG0141,11.50347402053696,0.11917194420401866
BG0142,11.18291150339665,0.0928709524246143
BG0143,13.827892197512199,0.08517034078295005
BG0144,10.119260554498801,0.11071292051765358
BG0145,10.362250610551786,0.11276762414596675
BG0146,12.655978798022087,0.10365745802901749
BG0147,11.14469233052196,0.08897568878221396
BG0148,11.70063932453736,0.138715019365808
BG0149,11.525044536897665,0.09386088748392218
"""})
