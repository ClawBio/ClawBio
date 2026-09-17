"""Coalescent-simulation P-values for Tajima's D, R2 and Fu's Fs.

The simulator computes each replicate's statistics directly from the simulated
genealogy rather than writing out sequences, because that is what makes 10 000
replicates affordable in pure Python. The central tests therefore check that fast
path against the skill's own sequence-based code: every simulated tree is also written
out as sequences and run through analyse_region, and the two must agree exactly.

The simulator itself is checked against exact results for the standard coalescent
rather than approximate critical values: E[S] = theta * a1 and
Var(S) = a1 * theta + a2 * theta^2 (Watterson 1975), E[k] = theta and
Var(k) = (n + 1) theta / (3 (n - 1)) + 2 (n^2 + n + 3) theta^2 / (9 n (n - 1))
(Tajima 1983).
"""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import dnasp  # noqa: E402


def _harmonic(n: int, power: int) -> float:
    return sum(1.0 / i ** power for i in range(1, n))


# ----------------------------------------------------------------------------------------------
# The refactored Fu's Fs helper must not change the existing statistic
# ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("n,H,k", [(10, 3, 0.822222), (34, 34, 16.1355), (53, 53, 61.0), (6, 2, 1.5), (4, 3, 1.166667)])
def test_fu_fs_statistic_equals_compute_fu_fs(n, H, k):
    expected = dnasp.compute_fu_fs(["A"] * n, H, k).Fs
    assert dnasp.fu_fs_statistic(n, H, k) == pytest.approx(expected, rel=0, abs=1e-12)


def test_fu_fs_statistic_undefined_cases():
    assert dnasp.fu_fs_statistic(10, 1, 2.0) is None
    assert dnasp.fu_fs_statistic(10, 3, 0.0) is None
    assert dnasp.fu_fs_statistic(1, 1, 1.0) is None


# ----------------------------------------------------------------------------------------------
# Fast tree statistics equal the sequence-based statistics on the same genealogy
# ----------------------------------------------------------------------------------------------

def _replicate(seed: int, n: int, S: int):
    rng = random.Random(seed)
    tree = dnasp._coalescent_tree(rng, n)
    counts = dnasp._place_mutations(rng, tree, n_mutations=S)
    return tree, counts


@pytest.mark.parametrize("seed", range(40))
@pytest.mark.parametrize("n,S", [(3, 4), (5, 1), (8, 12), (20, 35), (53, 60)])
def test_tree_statistics_match_analyse_region(seed, n, S):
    tree, counts = _replicate(seed, n, S)
    fast = dnasp._tree_statistics(tree, counts)

    seqs = dnasp._tree_sequences(tree, counts)
    assert len(seqs) == n and all(len(s) == S for s in seqs)
    names = [f"s{i}" for i in range(n)]
    ref = dnasp.analyse_region(seqs, names, "sim", len(seqs[0]))

    assert fast.S == ref.S == S
    assert fast.k == pytest.approx(ref.k, abs=1e-9)
    assert fast.H == ref.H
    assert fast.TajimaD == pytest.approx(ref.TajimaD, abs=1e-9) if ref.TajimaD is not None else fast.TajimaD is None
    assert fast.R2 == pytest.approx(ref.R2, abs=1e-9) if ref.R2 is not None else fast.R2 is None
    ref_fs = dnasp.compute_fu_fs(seqs, ref.H, ref.k).Fs
    assert fast.Fs == pytest.approx(ref_fs, abs=1e-9) if ref_fs is not None else fast.Fs is None


def test_every_replicate_given_s_has_exactly_s_segregating_sites():
    rng = random.Random(7)
    for _ in range(300):
        tree = dnasp._coalescent_tree(rng, 15)
        counts = dnasp._place_mutations(rng, tree, n_mutations=9)
        assert sum(counts) == 9


def test_tree_has_2n_minus_1_nodes_and_one_root():
    tree = dnasp._coalescent_tree(random.Random(3), 12)
    assert len(tree.parent) == 2 * 12 - 1
    assert sum(1 for p in tree.parent if p < 0) == 1
    root = tree.parent.index(-1)
    assert tree.desc[root] == 12


# ----------------------------------------------------------------------------------------------
# The simulator reproduces exact moments of the standard coalescent (theta conditioning)
# ----------------------------------------------------------------------------------------------

def test_theta_conditioning_reproduces_exact_moments_of_s_and_k():
    n, theta, reps = 10, 5.0, 12000
    rng = random.Random(20260917)
    S_vals, k_vals = [], []
    for _ in range(reps):
        tree = dnasp._coalescent_tree(rng, n)
        counts = dnasp._place_mutations(rng, tree, theta=theta)
        st = dnasp._tree_statistics(tree, counts)
        S_vals.append(st.S)
        k_vals.append(st.k)

    a1, a2 = _harmonic(n, 1), _harmonic(n, 2)
    ES, VS = theta * a1, a1 * theta + a2 * theta ** 2
    Ek = theta
    Vk = (n + 1) * theta / (3 * (n - 1)) + 2 * (n ** 2 + n + 3) * theta ** 2 / (9 * n * (n - 1))

    def mean(v): return sum(v) / len(v)
    def var(v):
        m = mean(v); return sum((x - m) ** 2 for x in v) / (len(v) - 1)

    # four standard errors for the means; 12 per cent relative for the variances
    assert abs(mean(S_vals) - ES) < 4 * math.sqrt(VS / reps)
    assert abs(mean(k_vals) - Ek) < 4 * math.sqrt(Vk / reps)
    assert var(S_vals) == pytest.approx(VS, rel=0.12)
    assert var(k_vals) == pytest.approx(Vk, rel=0.12)


# ----------------------------------------------------------------------------------------------
# P-value definitions
# ----------------------------------------------------------------------------------------------

def test_tail_counts_are_inclusive():
    null = [-2.0, -1.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
    assert dnasp._tail_counts(-1.0, null) == (3, 7, 8)


def test_tail_counts_ignore_undefined_replicates():
    assert dnasp._tail_counts(0.5, [None, 0.0, 1.0, None]) == (1, 1, 2)
    assert dnasp._tail_counts(0.5, [None, None]) is None
    assert dnasp._tail_counts(None, [0.0, 1.0]) is None


def test_monte_carlo_p_value_counts_the_data_as_a_replicate_and_is_never_zero():
    # Phipson and Smyth (2010): P = (b + 1) / (N + 1) for b of N replicates as extreme.
    assert dnasp._monte_carlo_p(0, 10000) == pytest.approx(1 / 10001)
    assert dnasp._monte_carlo_p(17, 10000) == pytest.approx(18 / 10001)
    assert dnasp._monte_carlo_p(10000, 10000) == 1.0
    assert dnasp._monte_carlo_p(None, 10) is None


def test_two_tailed_is_twice_the_smaller_tail_capped_at_one():
    assert dnasp._two_tailed(4 / 9, 8 / 9) == pytest.approx(8 / 9)
    assert dnasp._two_tailed(0.9, 0.95) == 1.0
    # zero replicates in the lower tail of 10 000: 2 / 10 001, not 1 / 10 001
    lower, upper = dnasp._monte_carlo_p(0, 10000), dnasp._monte_carlo_p(10000, 10000)
    assert dnasp._two_tailed(lower, upper) == pytest.approx(2 / 10001)


# ----------------------------------------------------------------------------------------------
# The test as a whole
# ----------------------------------------------------------------------------------------------

def _region(n: int = 20, seed: int = 11):
    rng = random.Random(seed)
    tree = dnasp._coalescent_tree(rng, n)
    counts = dnasp._place_mutations(rng, tree, n_mutations=25)
    seqs = dnasp._tree_sequences(tree, counts)
    return dnasp.analyse_region(seqs, [f"s{i}" for i in range(n)], "r", len(seqs[0]))


def test_coalescent_test_reports_valid_probabilities_and_its_settings():
    rs = _region()
    ct = dnasp.coalescent_test(rs, n_sim=400, given="S", seed=5, label="r")
    assert ct.n_sim == 400 and ct.given == "S" and ct.seed == 5
    assert ct.n == rs.n and ct.S == rs.S
    for p in (ct.TajimaD_p_lower, ct.TajimaD_p_upper, ct.TajimaD_p_two_tailed, ct.R2_p_lower, ct.Fs_p_lower):
        assert p is not None and 1 / 401 <= p <= 1.0
    assert ct.n_valid_TajimaD == 400
    # the stored P-values are (b + 1) / (N + 1) from the stored counts, so DnaSP's b / N is recoverable
    assert ct.TajimaD_p_lower == pytest.approx((ct.TajimaD_count_lower + 1) / 401)
    assert ct.TajimaD_p_upper == pytest.approx((ct.TajimaD_count_upper + 1) / 401)
    assert ct.R2_p_lower == pytest.approx((ct.R2_count_lower + 1) / (ct.n_valid_R2 + 1))
    assert ct.Fs_p_lower == pytest.approx((ct.Fs_count_lower + 1) / (ct.n_valid_Fs + 1))
    assert ct.TajimaD_count_lower + ct.TajimaD_count_upper >= 400
    assert ct.TajimaD_p_two_tailed == pytest.approx(min(1.0, 2 * min(ct.TajimaD_p_lower, ct.TajimaD_p_upper)))


def test_same_seed_gives_identical_results_and_the_label_separates_regions():
    rs = _region()
    a = dnasp.coalescent_test(rs, n_sim=200, given="S", seed=42, label="r")
    b = dnasp.coalescent_test(rs, n_sim=200, given="S", seed=42, label="r")
    c = dnasp.coalescent_test(rs, n_sim=200, given="S", seed=42, label="other-region")
    assert a == b
    assert (a.TajimaD_p_lower, a.R2_p_lower, a.Fs_p_lower) != (c.TajimaD_p_lower, c.R2_p_lower, c.Fs_p_lower)


def test_theta_conditioning_records_theta_and_can_lose_replicates():
    rs = _region()
    ct = dnasp.coalescent_test(rs, n_sim=300, given="theta", seed=1, label="r")
    assert ct.theta == pytest.approx(rs.ThetaW)
    assert 0 < ct.n_valid_TajimaD <= 300


def test_no_simulation_without_variation():
    seqs = ["ACGT"] * 6
    rs = dnasp.analyse_region(seqs, [f"s{i}" for i in range(6)], "r", 4)
    ct = dnasp.coalescent_test(rs, n_sim=100, given="S", seed=1, label="r")
    assert ct.TajimaD_p_two_tailed is None and ct.R2_p_lower is None and ct.Fs_p_lower is None
    assert ct.TajimaD_count_lower is None and ct.R2_count_lower is None and ct.Fs_count_lower is None
    assert "no segregating site" in ct.note


def test_invalid_settings_are_refused():
    rs = _region()
    with pytest.raises(ValueError):
        dnasp.coalescent_test(rs, n_sim=-1, given="S", seed=1, label="r")
    with pytest.raises(ValueError):
        dnasp.coalescent_test(rs, n_sim=10, given="eta", seed=1, label="r")


# ----------------------------------------------------------------------------------------------
# Wiring: run_analysis, the command line and every output
# ----------------------------------------------------------------------------------------------

def _alignment():
    return dnasp.load_alignment(SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas")


def test_run_analysis_without_n_sim_adds_nothing():
    results = dnasp.run_analysis(_alignment())
    assert results.get("coalescent") is None


def test_run_analysis_with_n_sim_adds_the_test():
    results = dnasp.run_analysis(_alignment(), n_sim=150, sim_given="S", sim_seed=9)
    ct = results["coalescent"]
    assert ct is not None and ct.n_sim == 150 and ct.seed == 9


def test_cli_writes_p_values_to_summary_and_report(tmp_path):
    src = SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas"
    out = tmp_path / "out"
    rc = dnasp.main(["--input", str(src), "--analysis", "polymorphism,fufs", "--n-sim", "120",
                     "--sim-seed", "3", "--output", str(out)])
    assert rc == 0
    import json
    summary = json.loads((out / "summary.json").read_text())
    ct = summary["coalescent"]
    assert ct["n_sim"] == 120 and ct["seed"] == 3 and ct["given"] == "S"
    assert ct["TajimaD_p_two_tailed"] is not None and ct["Fs_p_lower"] is not None
    report = (out / "report.md").read_text()
    assert "Coalescent simulation" in report
    assert "120 replicates" in report


def test_cli_records_a_generated_seed_when_none_is_given(tmp_path):
    src = SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas"
    out = tmp_path / "out"
    assert dnasp.main(["--input", str(src), "--n-sim", "50", "--output", str(out)]) == 0
    import json
    seed = json.loads((out / "summary.json").read_text())["coalescent"]["seed"]
    assert isinstance(seed, int)
    assert f"--sim-seed {seed}" in (out / "report.md").read_text()


def test_report_without_simulation_makes_no_significance_claim(tmp_path):
    src = SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas"
    out = tmp_path / "out"
    assert dnasp.main(["--input", str(src), "--output", str(out)]) == 0
    report = (out / "report.md").read_text()
    for claim in ("Consistent with neutrality", "suggests population expansion", "selective sweep or population expansion"):
        assert claim not in report
    assert "significance not assessed" in report


def test_cli_rejects_a_negative_replicate_count(tmp_path):
    src = SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas"
    with pytest.raises(SystemExit):
        dnasp.main(["--input", str(src), "--n-sim", "-5", "--output", str(tmp_path / "o")])


# ----------------------------------------------------------------------------------------------
# Small samples, report consistency and replayable provenance (added after code review)
# ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("S", [1, 3, 7])
def test_two_sequences_tree_statistics_match_analyse_region(seed, S):
    tree, counts = _replicate(seed, 2, S)
    fast = dnasp._tree_statistics(tree, counts)
    seqs = dnasp._tree_sequences(tree, counts)
    ref = dnasp.analyse_region(seqs, ["a", "b"], "sim", len(seqs[0]))
    assert fast.S == ref.S and fast.H == ref.H
    assert fast.k == pytest.approx(ref.k, abs=1e-9)
    assert fast.TajimaD is None and ref.TajimaD is None
    assert fast.R2 == pytest.approx(ref.R2, abs=1e-9)


def test_two_sequences_still_give_r2_and_fs_p_values():
    seqs = ["AAAAAA", "AGAGAA"]
    rs = dnasp.analyse_region(seqs, ["a", "b"], "r", 6)
    ct = dnasp.coalescent_test(rs, n_sim=200, given="S", seed=4, label="r")
    assert ct.TajimaD_p_two_tailed is None and ct.n_valid_TajimaD == 0
    assert ct.R2_p_lower is not None and ct.n_valid_R2 == 200
    assert ct.Fs_p_lower is not None


def _one_site_three_sequences():
    seqs = ["AAAA", "AAGA", "AAAA"]
    return dnasp.analyse_region(seqs, ["a", "b", "c"], "r", 4)


def test_report_shows_r2_p_value_even_when_tajimas_d_has_no_valid_replicate(tmp_path):
    rs = _one_site_three_sequences()
    ct = dnasp.coalescent_test(rs, n_sim=100, given="S", seed=2, label="r")
    assert ct.n_valid_R2 > 0
    src = tmp_path / "three.fas"
    src.write_text(">a\nAAAA\n>b\nAAGA\n>c\nAAAA\n")
    out = tmp_path / "out"
    assert dnasp.main(["--input", str(src), "--n-sim", "100", "--sim-seed", "2", "--output", str(out)]) == 0
    report = (out / "report.md").read_text()
    table = report[report.index("## Neutrality Tests"):report.index("## Coalescent simulation")]
    r2_row = next(line for line in table.splitlines() if line.startswith("| Ramos-Onsins"))
    assert "not assessed" not in r2_row
    assert "(lower tail)" in r2_row


def test_the_generated_seed_is_in_the_replayable_command(tmp_path):
    src = SKILL_DIR / "tests" / "fixtures" / "inputs" / "Ex_n1.fas"
    out = tmp_path / "out"
    assert dnasp.main(["--input", str(src), "--n-sim", "30", "--output", str(out)]) == 0
    import json
    seed = json.loads((out / "summary.json").read_text())["coalescent"]["seed"]
    commands = next(out.rglob("commands.sh")).read_text()
    assert f"--sim-seed {seed}" in commands


def test_vcf_chromosomes_share_the_seed_but_are_labelled_apart(tmp_path):
    import json
    vcf = SKILL_DIR / "tests" / "fixtures" / "inputs" / "vcf" / "Data_Example_DiploidPhased.vcf"
    out = tmp_path / "out"
    assert dnasp.main(["--vcf", str(vcf), "--n-sim", "40", "--sim-seed", "8", "--output", str(out)]) == 0
    records = [json.loads(p.read_text())["coalescent"] for p in sorted(out.glob("*/summary.json"))]
    assert len(records) == 3
    assert {r["seed"] for r in records} == {8}
    assert len({r["label"] for r in records}) == 3


def test_observed_statistics_come_from_the_completely_deleted_alignment():
    aln = dnasp.load_alignment(SKILL_DIR / "tests" / "fixtures" / "inputs" / "rp49_36.nex")
    results = dnasp.run_analysis(aln, n_sim=20, sim_seed=1)
    g, ct = results["global"], results["coalescent"]
    assert g.L_net < aln.L                      # gapped columns were removed
    assert (ct.n, ct.S, ct.TajimaD, ct.R2) == (g.n, g.S, g.TajimaD, g.R2)


# ---------------------------------------------------------------------------
# Display of P-values from finite replicates (added after manuscript review)
# ---------------------------------------------------------------------------

def test_p_values_are_displayed_without_ever_rounding_to_zero():
    # four decimals from 0.001 upward; below it two significant digits, trailing zero kept
    assert dnasp._fmt_p(0.1706) == "0.1706"
    assert dnasp._fmt_p(18 / 10001) == "0.0018"
    assert dnasp._fmt_p(0.001) == "0.0010"
    assert dnasp._fmt_p(0.00099995) == "0.0010"
    assert dnasp._fmt_p(1 / 10001) == "0.00010"
    assert dnasp._fmt_p(2 / 10001) == "0.00020"
    assert dnasp._fmt_p(1 / 100001) == "0.000010"
    assert dnasp._fmt_p(1 / 10 ** 13) == "0.00000000000010"
    assert dnasp._fmt_p(None) == "n.a."


def test_interpretation_states_the_p_value():
    assert "(P = 0.00010)" in dnasp._p_interp(1 / 10001, "lower")
    assert "No significant departure" in dnasp._p_interp(0.1706, "two-tailed")


def test_report_gives_counts_and_corrected_p_values(tmp_path):
    out = tmp_path / "out"
    rc = dnasp.main(["--input", str(SKILL_DIR / "tests" / "fixtures" / "inputs" / "rp49_36.nex"),
                      "--analysis", "polymorphism,fufs", "--n-sim", "50", "--sim-seed", "3",
                      "--output", str(out)])
    assert rc == 0
    ct = __import__("json").loads((out / "summary.json").read_text())["coalescent"]
    report = (out / "report.md").read_text()
    section = report[report.index("## Coalescent simulation"):]
    assert "Replicates <= observed" in section and "Replicates >= observed" in section
    assert f"| {ct['Fs_count_lower']} |" in section
    assert dnasp._fmt_p(ct["Fs_p_lower"]) in section
    assert "(b + 1)/(N + 1)" in section
    assert "0 of " not in report and "P <= " not in report
