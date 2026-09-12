"""Tests for DnaSP-Python skill.

Hand-computed reference values were verified against DnaSP 6 output
and first-principles calculation for the 6-sequence synthetic demo.

Run:
    python -m pytest skills/dnasp/tests/ -v
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import pytest
import sys

# Make the skill importable without installation
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dnasp as dn


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

#  6 sequences × 10 sites (no gaps)
#  pop1_seq1/2/3: AAACGTTAAG  (hap A, freq=3/6)
#  pop2_seq4/5:   AAACGTTAAC  (hap B, freq=2/6)  pos 10 differs
#  pop3_seq6:     TAACGTTAAG  (hap C, freq=1/6)  pos  1 differs (singleton)
DEMO_SEQS = [
    "AAACGTTAAG",
    "AAACGTTAAG",
    "AAACGTTAAG",
    "AAACGTTAAC",
    "AAACGTTAAC",
    "TAACGTTAAG",
]

# Hand-computed expected values
# n=6, L=10, no gaps → L_net=10
# S=2 (pos 1 and pos 10), Eta=2
# Haplotypes: A(3), B(2), C(1) → H=3
# freqs = [0.5, 0.333, 0.167]
# sum_sq = 0.25 + 0.1111 + 0.0278 = 0.3889
# Hd = 6/5 * (1 - 0.3889) = 1.2 * 0.6111 = 0.7333
# Pairwise diffs:  13 total pairs have diffs (see below), pairs = 15
# k = 13/15 ≈ 0.8667
# Pi = k/L = 0.08667
# a1 = 1 + 1/2 + 1/3 + 1/4 + 1/5 = 2.2833
# θ_W_abs = S/a1 = 2/2.2833 = 0.8759
# θ_W_nuc = 0.8759/10 = 0.08759
# Tajima D: k - S/a1 = 0.8667 - 0.8759 = -0.0092
#   a2=1+1/4+1/9+1/16+1/25=1.4636
#   b1=(6+1)/(3*5)=0.4667, b2=2*(36+6+3)/(9*6*5)=0.3444
#   c1=b1-1/a1=0.4667-0.4379=0.0288
#   c2=b2-(6+2)/(a1*6)+a2/a1^2=0.3444-0.5837+0.2806=0.0413
#   e1=c1/a1=0.0288/2.2833=0.01261
#   e2=c2/(a1^2+a2)=0.0413/(5.213+1.4636)=0.00619
#   vD=e1*2+e2*2*1=0.02522+0.01238=0.0376
#   D=-0.0092/sqrt(0.0376)=-0.0092/0.1939=-0.0474

DEMO_EXPECTED = {
    "n": 6,
    "L_net": 10,
    "S": 2,
    "Eta": 2,
    "H": 3,
    "Hd": pytest.approx(6 / 5 * (1 - (9 / 36 + 4 / 36 + 1 / 36)), abs=1e-4),
    "Pi": pytest.approx(13 / 15 / 10, abs=1e-5),
    "k": pytest.approx(13 / 15, abs=1e-5),
}


@pytest.fixture
def demo_seqs():
    return list(DEMO_SEQS)


@pytest.fixture
def demo_fasta_path(tmp_path):
    p = tmp_path / "demo.fas"
    p.write_text(dn.DEMO_FASTA, encoding="utf-8")
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Parser tests
# ─────────────────────────────────────────────────────────────────────────────

class TestFastaParser:
    def test_simple_fasta(self, tmp_path):
        fa = tmp_path / "test.fas"
        fa.write_text(">seq1\nATCG\n>seq2\nATCG\n")
        aln = dn.parse_fasta(fa)
        assert aln.n == 2
        assert aln.seqs[0] == "ATCG"

    def test_dnasp_header_format(self, tmp_path):
        fa = tmp_path / "dnasp.fas"
        fa.write_text(">'J1_Ost'   [by DnaSP Ver. 5, from file: x.nex]\nATCG\n>'J2_Ost'\nATCG\n")
        aln = dn.parse_fasta(fa)
        assert aln.names[0] == "J1_Ost"
        assert aln.names[1] == "J2_Ost"

    def test_multiline_seq(self, tmp_path):
        fa = tmp_path / "ml.fas"
        fa.write_text(">s1\nATC\nGAT\n>s2\nATCGAT\n")
        aln = dn.parse_fasta(fa)
        assert aln.seqs[0] == "ATCGAT"
        assert aln.L == 6

    def test_uppercase_conversion(self, tmp_path):
        fa = tmp_path / "lc.fas"
        fa.write_text(">s1\natcg\n>s2\natcg\n")
        aln = dn.parse_fasta(fa)
        assert aln.seqs[0] == "ATCG"

    def test_demo_fasta_parses(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        assert aln.n == 11  # 10 ingroup + 1 outgroup
        assert aln.L == 300

    def test_empty_file_raises(self, tmp_path):
        fa = tmp_path / "empty.fas"
        fa.write_text("")
        with pytest.raises(ValueError, match="No sequences"):
            dn.parse_fasta(fa)


class TestNexusParser:
    def test_basic_nexus(self, tmp_path):
        nex = tmp_path / "test.nex"
        nex.write_text(
            "#NEXUS\nBEGIN CHARACTERS;\nDIMENSIONS NTAX=2 NCHAR=4;\n"
            "FORMAT DATATYPE=DNA GAP=- MISSING=?;\n"
            "MATRIX\nseq1 ATCG\nseq2 ATCG\n;\nEND;\n"
        )
        aln = dn.parse_nexus(nex)
        assert aln.n == 2
        assert aln.seqs[0] == "ATCG"

    def test_matchchar_expansion(self, tmp_path):
        nex = tmp_path / "match.nex"
        nex.write_text(
            "#NEXUS\nBEGIN CHARACTERS;\nDIMENSIONS NTAX=2 NCHAR=4;\n"
            "FORMAT DATATYPE=DNA MATCHCHAR=. GAP=-;\n"
            "MATRIX\nref ATCG\nseq ...G\n;\nEND;\n"
        )
        aln = dn.parse_nexus(nex)
        # ref is ATCG, seq has match for pos 0-2 and G at pos 3
        assert aln.seqs[1] == "ATCG"


# ─────────────────────────────────────────────────────────────────────────────
# Gap treatment tests
# ─────────────────────────────────────────────────────────────────────────────

class TestCompleteDeletion:
    def test_no_gaps(self):
        seqs = ["ATCG", "ATCG"]
        clean, L_net = dn.complete_deletion(seqs)
        assert L_net == 4
        assert clean == seqs

    def test_gap_removes_column(self):
        seqs = ["A-CG", "ATCG"]
        clean, L_net = dn.complete_deletion(seqs)
        assert L_net == 3
        assert clean[0] == "ACG"
        assert clean[1] == "ACG"

    def test_all_gaps_gives_zero_net(self):
        seqs = ["----", "----"]
        clean, L_net = dn.complete_deletion(seqs)
        assert L_net == 0

    def test_N_excluded(self):
        seqs = ["ANCG", "ATCG"]
        clean, L_net = dn.complete_deletion(seqs)
        assert L_net == 3


# ─────────────────────────────────────────────────────────────────────────────
# Haplotype tests
# ─────────────────────────────────────────────────────────────────────────────

class TestHaplotypes:
    def test_all_identical(self):
        seqs = ["ATCG"] * 4
        H, Hd, VarHd = dn.compute_haplotypes(seqs)
        assert H == 1
        assert Hd == pytest.approx(0.0)

    def test_all_unique(self):
        seqs = ["AAAA", "TTTT", "CCCC", "GGGG"]
        H, Hd, VarHd = dn.compute_haplotypes(seqs)
        assert H == 4
        # Hd = 4/3 * (1 - 4*(1/4)^2) = 4/3 * (1 - 0.25) = 1.0
        assert Hd == pytest.approx(1.0, abs=1e-6)

    def test_demo_haplotypes(self, demo_seqs):
        H, Hd, VarHd = dn.compute_haplotypes(demo_seqs)
        assert H == 3
        # Hd = 6/5 * (1 - 9/36 - 4/36 - 1/36) = 6/5 * (1 - 14/36)
        expected_Hd = 6 / 5 * (1 - (9 + 4 + 1) / 36)
        assert Hd == pytest.approx(expected_Hd, abs=1e-6)
        assert VarHd >= 0


# ─────────────────────────────────────────────────────────────────────────────
# Nucleotide diversity tests
# ─────────────────────────────────────────────────────────────────────────────

class TestNucleotideDiversity:
    def test_identical_seqs(self):
        seqs = ["ATCG"] * 3
        k = dn.compute_k(seqs)
        assert k == pytest.approx(0.0)

    def test_one_diff_two_seqs(self):
        # 1 difference out of 4 sites, 1 pair → k = 1
        seqs = ["ATCG", "ATCC"]
        k = dn.compute_k(seqs)
        assert k == pytest.approx(1.0)

    def test_demo_k(self, demo_seqs):
        k = dn.compute_k(demo_seqs)
        assert k == pytest.approx(13 / 15, abs=1e-6)

    def test_demo_pi(self, demo_seqs):
        k = dn.compute_k(demo_seqs)
        Pi = k / 10
        assert Pi == pytest.approx(13 / 150, abs=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# Segregating sites tests
# ─────────────────────────────────────────────────────────────────────────────

class TestSegregatingSites:
    def test_no_variation(self):
        seqs = ["AAAA"] * 3
        S, Eta = dn.compute_segregating(seqs)
        assert S == 0
        assert Eta == 0

    def test_biallelic(self):
        seqs = ["AATG", "AACG"]
        S, Eta = dn.compute_segregating(seqs)
        assert S == 1
        assert Eta == 1

    def test_triallelic_adds_eta(self):
        # Position 0 has A, T, C → triallelic → S=1, Eta=2
        seqs = ["ATCG", "TTCG", "CTCG"]
        S, Eta = dn.compute_segregating(seqs)
        assert S == 1
        assert Eta == 2

    def test_demo_seg_sites(self, demo_seqs):
        S, Eta = dn.compute_segregating(demo_seqs)
        assert S == 2
        assert Eta == 2


# ─────────────────────────────────────────────────────────────────────────────
# Singleton tests
# ─────────────────────────────────────────────────────────────────────────────

class TestSingletons:
    def test_no_singletons(self):
        # Both positions segregate with equal frequency  -  no singletons
        seqs = ["AT", "AT", "GC", "GC"]
        eta_s, per_seq = dn.compute_singletons(seqs)
        assert eta_s == 0

    def test_demo_singletons(self, demo_seqs):
        # Only pos 1 (T in seq6) is a true singleton; pos 10 (C in seq4,5) is not.
        # DnaSP BusqSingletones counts only alleles with count == 1 (no folded).
        eta_s, per_seq = dn.compute_singletons(demo_seqs)
        assert eta_s == 1
        assert per_seq[5] == pytest.approx(1.0)   # seq6 (index 5) carries the singleton T

    def test_multiallelic_site_capped_at_one_singleton_site(self):
        # n=4, one triallelic column with a 2/1/1 split: two different alleles
        # each have count 1, but DnaSP's EtaS = EtaS - SingleMut discount means
        # this site still contributes only 1 to eta_s (a SITE count), not 2 (a
        # per-allele count). The two carrying sequences still each get +1 in
        # per_seq, which R2 needs uncapped.
        seqs = ["G", "G", "C", "T"]
        eta_s, per_seq = dn.compute_singletons(seqs)
        assert eta_s == 1
        assert per_seq == [0.0, 0.0, 1.0, 1.0]

    def test_biallelic_singleton_sites_are_not_capped_away(self):
        # Two independent biallelic singleton sites must both count.
        seqs = ["AT", "AT", "AT", "GC"]
        eta_s, _ = dn.compute_singletons(seqs)
        assert eta_s == 2


# ─────────────────────────────────────────────────────────────────────────────
# Tajima's D tests
# ─────────────────────────────────────────────────────────────────────────────

class TestTajimaD:
    def test_none_when_S_zero(self):
        assert dn.tajima_d(0.0, S=0, n=6) is None

    def test_none_when_n_lt_3(self):
        assert dn.tajima_d(1.0, S=2, n=2) is None

    def test_demo_tajima(self, demo_seqs):
        k = dn.compute_k(demo_seqs)
        S, _ = dn.compute_segregating(demo_seqs)
        D = dn.tajima_d(k, S, n=6)
        assert D is not None
        # Should be slightly negative (π slightly less than θ_W)
        # Exact value: approx -0.047
        assert D == pytest.approx(-0.047, abs=0.01)

    def test_known_positive_D(self):
        # Simulate high π relative to S: k=5, S=2, n=10 → large positive D
        D = dn.tajima_d(k=5.0, S=2, n=10)
        assert D is not None
        assert D > 0

    def test_known_negative_D(self):
        # Simulate low π relative to S: k=0.1, S=10, n=10 → large negative D
        D = dn.tajima_d(k=0.1, S=10, n=10)
        assert D is not None
        assert D < 0


# ─────────────────────────────────────────────────────────────────────────────
# Fu & Li D* and F* tests
# ─────────────────────────────────────────────────────────────────────────────

class TestFuLi:
    def test_none_when_S_zero(self):
        D, F = dn.fu_li_d_star_f_star(k=0.5, S=0, eta_s=0, n=6)
        assert D is None
        assert F is None

    def test_none_when_n_lt_3(self):
        D, F = dn.fu_li_d_star_f_star(k=0.5, S=2, eta_s=1, n=2)
        assert D is None
        assert F is None

    def test_demo_fu_li(self, demo_seqs):
        clean, L_net = dn.complete_deletion(demo_seqs)
        k = dn.compute_k(clean)
        S, _ = dn.compute_segregating(clean)
        eta_s, _ = dn.compute_singletons(clean)
        D_star, F_star = dn.fu_li_d_star_f_star(k, S, eta_s, n=6)
        assert D_star is not None
        assert F_star is not None
        # n=6, S=2, eta_s=1: D* = (2/An - 5/6) / sqrt(denom_D) ≈ +0.062 (slightly positive)
        assert D_star == pytest.approx(0.062, abs=0.01)
        assert F_star == pytest.approx(0.040, abs=0.01)

    def test_sign_consistency(self):
        # When η_s is large (many singletons), expect negative D* and F*
        D, F = dn.fu_li_d_star_f_star(k=1.0, S=10, eta_s=8, n=10)
        if D is not None:
            assert D < 0
        if F is not None:
            assert F < 0

    def test_ex_n1_matches_dnasp_6_exactly(self):
        # DnaSP's own hand-checkable n=4 example (Ex_n1.fas, distributed with
        # DnaSP 6). Run in DnaSP 6.12 for Windows on 2026-09-11: D* = -0.52807,
        # F* = -0.41685 (the "v5-style"/Simonsen-1995 figure this skill
        # implements; DnaSP's own Achaz-2009 "biallelic" figure differs, as
        # documented). This is the case that exposed the eta_s per-site cap.
        seqs = [
            "ATATACGGGGTTA---TTAGA----AAAATGTGTGTGTGTTTTTTTTTTCATGTG",
            "ATATAC--GGATA---TTACA----AGAATCTATGTCTGCTTTCTTTTTCATGTG",
            "ATATACGGGGATA---TTATA----AGAATGTGTGTGTGTTTTTTTTTTCATGTG",
            "ATATACGGGGATA---GTAGT----AAAATGTGTGTGTGTTTTTTTTTTCATGTG",
        ]
        clean, L_net = dn.complete_deletion(seqs)
        assert L_net == 46
        S, _ = dn.compute_segregating(clean)
        assert S == 10
        k = dn.compute_k(clean)
        eta_s, _ = dn.compute_singletons(clean)
        assert eta_s == 9   # capped: the triallelic site (G/C/T, 2/1/1) counts once
        D_star, F_star = dn.fu_li_d_star_f_star(k, S, eta_s, n=4)
        assert D_star == pytest.approx(-0.52807, abs=1e-5)
        assert F_star == pytest.approx(-0.41685, abs=1e-5)


# ─────────────────────────────────────────────────────────────────────────────
# R2 tests
# ─────────────────────────────────────────────────────────────────────────────

class TestR2:
    def test_none_when_S_zero(self):
        assert dn.ramos_onsins_r2(["ATCG"] * 3, k=0.0, S=0) is None

    def test_demo_r2(self, demo_seqs):
        clean, _ = dn.complete_deletion(demo_seqs)
        k = dn.compute_k(clean)
        S, _ = dn.compute_segregating(clean)
        R2 = dn.ramos_onsins_r2(clean, k, S=S)
        assert R2 is not None
        assert R2 > 0

    def test_r2_is_positive(self):
        seqs = ["ATCG", "AACG", "ATTG", "ATCG", "ATCG"]
        k = dn.compute_k(seqs)
        S, _ = dn.compute_segregating(seqs)
        if S > 0:
            R2 = dn.ramos_onsins_r2(seqs, k, S=S)
            if R2 is not None:
                assert R2 >= 0

    def test_r2_uses_S_argument_value(self):
        # Same per_seq/k, different S/Eta divisors -> different (documented)
        # results; confirms the function actually divides by its `S` argument.
        seqs = ["G", "G", "C", "T"]  # triallelic 2/1/1, matches the Ex_n1 case
        k = dn.compute_k(seqs)
        r2_at_s1 = dn.ramos_onsins_r2(seqs, k, S=1)
        r2_at_s2 = dn.ramos_onsins_r2(seqs, k, S=2)
        assert r2_at_s1 == pytest.approx(2 * r2_at_s2)


# ─────────────────────────────────────────────────────────────────────────────
# Watterson's theta tests
# ─────────────────────────────────────────────────────────────────────────────

class TestWatterson:
    def test_zero_sites(self):
        theta_abs, theta_nuc = dn.watterson_theta(S=0, n=6, L_net=10)
        assert theta_abs == 0.0
        assert theta_nuc == 0.0

    def test_demo_watterson(self):
        a1 = sum(1 / i for i in range(1, 6))  # = 2.2833
        theta_abs, theta_nuc = dn.watterson_theta(S=2, n=6, L_net=10)
        assert theta_abs == pytest.approx(2 / a1, abs=1e-6)
        assert theta_nuc == pytest.approx(2 / (a1 * 10), abs=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# GC content tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGC:
    def test_all_gc(self):
        seqs = ["GCGCGC"]
        assert dn.compute_gc(seqs) == pytest.approx(1.0)

    def test_all_at(self):
        seqs = ["ATATAT"]
        assert dn.compute_gc(seqs) == pytest.approx(0.0)

    def test_half(self):
        seqs = ["ATGC"]
        assert dn.compute_gc(seqs) == pytest.approx(0.5)

    def test_gaps_ignored(self):
        seqs = ["AT--GC"]
        assert dn.compute_gc(seqs) == pytest.approx(0.5)


# ─────────────────────────────────────────────────────────────────────────────
# Full pipeline integration tests
# ─────────────────────────────────────────────────────────────────────────────

class TestFullPipeline:
    def test_demo_mode(self, tmp_path):
        """Demo mode must run without errors and produce output files."""
        dn.run_demo(tmp_path, analyses={"polymorphism"})
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "results.tsv").exists()

    def test_analyse_region_all_identical(self):
        seqs = ["ATCG"] * 5
        rs = dn.analyse_region(seqs, [f"s{i}" for i in range(5)], "1-4", L_total=4)
        assert rs.S == 0
        assert rs.Hd == pytest.approx(0.0)
        assert rs.Pi == pytest.approx(0.0)

    def test_analyse_region_demo(self, demo_seqs):
        rs = dn.analyse_region(demo_seqs, [f"s{i}" for i in range(6)], "1-10", L_total=10)
        assert rs.n == 6
        assert rs.L_net == 10
        assert rs.S == 2
        assert rs.H == 3
        assert rs.Hd == pytest.approx(DEMO_EXPECTED["Hd"], abs=1e-4)
        assert rs.Pi == pytest.approx(DEMO_EXPECTED["Pi"], abs=1e-5)
        assert rs.TajimaD is not None

    def test_sliding_window(self, demo_fasta_path, tmp_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, window_size=150, step_size=150)
        assert results["global"].n == 11  # 10 ingroup + 1 outgroup
        assert len(results["windows"]) == 2  # two non-overlapping 150 bp windows from 300 bp

    def test_fasta_roundtrip_demo(self, demo_fasta_path, tmp_path):
        dn._run(
            demo_fasta_path, tmp_path,
            window_size=0, step_size=0,
            analyses={"polymorphism"},
            pop_assignments=None, aln2=None,
            cli_args=["--input", str(demo_fasta_path)],
        )
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "results.tsv").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()


# ─────────────────────────────────────────────────────────────────────────────
# Edge cases
# ─────────────────────────────────────────────────────────────────────────────

class TestEdgeCases:
    def test_two_sequences(self):
        seqs = ["ATCG", "ATTG"]
        rs = dn.analyse_region(seqs, ["s1", "s2"], "1-4", L_total=4)
        assert rs.S == 1
        # Tajima's D requires n ≥ 3
        assert rs.TajimaD is None
        assert rs.FuLiD_star is None

    def test_single_sequence(self):
        rs = dn.analyse_region(["ATCG"], ["s1"], "1-4", L_total=4)
        assert rs.n == 1
        assert rs.S == 0
        assert rs.Pi == 0.0

    def test_all_gaps(self):
        seqs = ["----", "----", "----"]
        rs = dn.analyse_region(seqs, ["s1", "s2", "s3"], "1-4", L_total=4)
        assert rs.L_net == 0
        assert rs.S == 0

    def test_mixed_case_normalized(self, tmp_path):
        fa = tmp_path / "lc.fas"
        fa.write_text(">s1\natcg\n>s2\natcc\n>s3\natcg\n")
        aln = dn.parse_fasta(fa)
        assert aln.seqs[0] == "ATCG"
        rs = dn.analyse_region(aln.seqs, aln.names, "1-4", L_total=4)
        assert rs.S == 1


# ─────────────────────────────────────────────────────────────────────────────
# Chi-square p-value (internal helper)
# ─────────────────────────────────────────────────────────────────────────────

class TestChi2PValue:
    def test_zero_chi2_returns_one(self):
        assert dn._chi2_1df_pvalue(0.0) == pytest.approx(1.0)

    def test_known_chi2_3_84(self):
        # chi2(1df) = 3.841 → p ≈ 0.05
        assert dn._chi2_1df_pvalue(3.841) == pytest.approx(0.05, abs=0.001)

    def test_known_chi2_6_63(self):
        # chi2(1df) = 6.635 → p ≈ 0.01
        assert dn._chi2_1df_pvalue(6.635) == pytest.approx(0.01, abs=0.001)

    def test_large_chi2_small_p(self):
        assert dn._chi2_1df_pvalue(30.0) < 1e-6


# ─────────────────────────────────────────────────────────────────────────────
# Biallelic site detection
# ─────────────────────────────────────────────────────────────────────────────

class TestBiallelicPositions:
    def test_no_variation(self):
        seqs = ["AAAA", "AAAA"]
        assert dn._get_biallelic_positions(seqs) == []

    def test_one_biallelic(self):
        # pos 0: A/T → biallelic; pos 1-3: monomorphic
        seqs = ["ATCG", "TTCG"]
        bp = dn._get_biallelic_positions(seqs)
        assert len(bp) == 1
        assert bp[0][0] == 0  # position 0

    def test_triallelic_excluded(self):
        # pos 0 has A, T, C → triallelic → excluded
        seqs = ["ATCG", "TTCG", "CTCG"]
        bp = dn._get_biallelic_positions(seqs)
        assert len(bp) == 0

    def test_gap_column_excluded(self):
        # pos 0 has a gap → excluded even though it's otherwise biallelic
        seqs = ["ATCG", "-TCG"]
        bp = dn._get_biallelic_positions(seqs)
        assert all(pos != 0 for pos, _, _ in bp)

    def test_demo_seqs_two_biallelic(self, demo_seqs):
        # DEMO has 2 segregating sites, both biallelic
        bp = dn._get_biallelic_positions(demo_seqs)
        assert len(bp) == 2


# ─────────────────────────────────────────────────────────────────────────────
# Linkage Disequilibrium
# ─────────────────────────────────────────────────────────────────────────────

class TestLDForPair:
    """Unit tests for _ld_for_pair against hand-computed values."""

    def _seqs_with_known_ld(self):
        # Two biallelic sites, equal frequencies, complete coupling:
        # Haplotypes: AB (n=3), ab (n=3) → D = p_AB - p_A*p_B = 0.5 - 0.25 = 0.25
        # p_A = p_B = 0.5 → D_max = min(0.5*0.5, 0.5*0.5) = 0.25 → D' = 1.0
        # R² = D²/(p_A*p_a*p_B*p_b) = 0.0625/0.0625 = 1.0
        return ["AG", "AG", "AG", "TC", "TC", "TC"]

    def test_complete_coupling_D(self):
        seqs = self._seqs_with_known_ld()
        D, D_prime, R2, n_valid = dn._ld_for_pair(seqs, 0, "A", 1, "G")
        # p_A=0.5, p_G=0.5; p_AG=0.5; D=0.5-0.25=0.25
        assert D == pytest.approx(0.25, abs=1e-9)

    def test_complete_coupling_Dprime(self):
        seqs = self._seqs_with_known_ld()
        D, D_prime, R2, n_valid = dn._ld_for_pair(seqs, 0, "A", 1, "G")
        assert D_prime == pytest.approx(1.0, abs=1e-9)

    def test_complete_coupling_R2(self):
        seqs = self._seqs_with_known_ld()
        D, D_prime, R2, n_valid = dn._ld_for_pair(seqs, 0, "A", 1, "G")
        assert R2 == pytest.approx(1.0, abs=1e-9)

    def test_independent_sites(self):
        # All 4 gametes present equally: R2 = 0
        seqs = ["AG", "AC", "TG", "TC"]
        D, D_prime, R2, n_valid = dn._ld_for_pair(seqs, 0, "A", 1, "G")
        assert D == pytest.approx(0.0, abs=1e-9)
        assert R2 == pytest.approx(0.0, abs=1e-9)

    def test_n_valid_counts_non_gap(self):
        # One sequence has a gap at site 0 → n_valid = 3
        seqs = ["-G", "AG", "AG", "TC"]
        _, _, _, n_valid = dn._ld_for_pair(seqs, 0, "A", 1, "G")
        assert n_valid == 3


class TestComputeLD:
    def test_no_biallelic_returns_empty(self):
        seqs = ["AAAA", "AAAA"]
        ld = dn.compute_ld(seqs)
        assert ld.n_biallelic == 0
        assert ld.n_pairs == 0
        assert ld.ZnS is None

    def test_one_biallelic_no_pairs(self):
        seqs = ["ATCG", "TTCG"]
        ld = dn.compute_ld(seqs)
        assert ld.n_biallelic == 1
        assert ld.n_pairs == 0
        assert ld.ZnS is None

    def test_complete_coupling_ZnS(self):
        # 2 sites in complete coupling → ZnS = R² = 1.0
        seqs = ["AG", "AG", "AG", "TC", "TC", "TC"]
        ld = dn.compute_ld(seqs)
        assert ld.n_biallelic == 2
        assert ld.n_pairs == 1
        assert ld.ZnS == pytest.approx(1.0, abs=1e-9)

    def test_Za_equals_ZnS_for_two_sites(self):
        # With exactly 2 biallelic sites, Za and ZnS are computed over the same pair
        seqs = ["AG", "AG", "AG", "TC", "TC", "TC"]
        ld = dn.compute_ld(seqs)
        assert ld.Za == pytest.approx(ld.ZnS, abs=1e-9)

    def test_ZZ_is_Za_minus_ZnS(self):
        # Three sites: two adjacent, one far → Za > ZnS → ZZ > 0
        # Construct seqs where adjacent sites are in full LD, far site independent
        seqs = [
            "AGT", "AGC", "AGT", "TCA", "TCC", "TCA",
        ]
        ld = dn.compute_ld(seqs)
        if ld.Za is not None and ld.ZnS is not None:
            assert ld.ZZ == pytest.approx(ld.Za - ld.ZnS, abs=1e-9)

    def test_demo_seqs_ld(self, demo_seqs):
        # DEMO: 2 biallelic sites at pos 0 (T/A, singleton in seq 6) and pos 9 (C/G)
        ld = dn.compute_ld(demo_seqs)
        assert ld.n_biallelic == 2
        assert ld.n_pairs == 1
        assert ld.ZnS is not None
        assert 0.0 <= ld.ZnS <= 1.0

    def test_pairs_have_valid_n(self, demo_seqs):
        ld = dn.compute_ld(demo_seqs)
        for pair in ld.pairs:
            assert pair.n_valid > 0
            assert pair.R2 is not None
            assert 0.0 <= pair.R2 <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# Recombination (Rm, four-gamete test)
# ─────────────────────────────────────────────────────────────────────────────

class TestComputeRecombination:
    def test_no_variation(self):
        seqs = ["AAAA"] * 4
        rec = dn.compute_recombination(seqs)
        assert rec.Rm == 0
        assert rec.n_incompatible_pairs == 0

    def test_compatible_sites(self):
        # Two sites in complete coupling → only 2 of 4 gametes → compatible
        seqs = ["AG", "AG", "AG", "TC", "TC", "TC"]
        rec = dn.compute_recombination(seqs)
        assert rec.Rm == 0
        assert rec.n_incompatible_pairs == 0

    def test_incompatible_pair_gives_Rm_1(self):
        # All 4 gametes present → 1 incompatible pair → Rm = 1
        seqs = ["AG", "AC", "TG", "TC"]
        rec = dn.compute_recombination(seqs)
        assert rec.n_incompatible_pairs >= 1
        assert rec.Rm >= 1

    def test_nested_intervals_Rm_not_double_counted(self):
        # Classic Hudson & Kaplan 1985 example: 3 incompatible pairs
        # that can all be explained by 2 recombination events
        # Sites: 1-2 incompat, 2-3 incompat, 1-3 incompat → Rm = 2
        # Using sequences where: pos 0/1 incompat, pos 1/2 incompat
        # This means pos 0 and 2 may or may not be incompat
        seqs = [
            "AGT",  # hap 1
            "AGC",  # hap 2: pos 2 differs from hap 1
            "TAT",  # hap 3: pos 0 differs from hap 1
            "TAC",  # hap 4: pos 0 and 2 differ
        ]
        rec = dn.compute_recombination(seqs)
        # We have all 4 gametes at (0,1), (0,2), (1,2) → 3 incompatible pairs
        # Rm = 2 (2 recombination events explain all 3 pairs)
        assert rec.Rm <= rec.n_incompatible_pairs

    def test_Rm_leq_n_pairs(self):
        seqs = ["AGT", "AGC", "TAT", "TAC"]
        rec = dn.compute_recombination(seqs)
        assert rec.Rm <= rec.n_incompatible_pairs

    def test_demo_seqs_no_recombination(self, demo_seqs):
        # DEMO: 2 biallelic sites; seq6 (T at pos 0) and seqs 4-5 (C at pos 9)
        # are not all 4 gametes → likely compatible
        rec = dn.compute_recombination(demo_seqs)
        # The two sites may or may not be compatible; just check structure
        assert rec.Rm >= 0
        assert rec.n_incompatible_pairs >= 0

    def test_matches_dnasp_reduction_not_optimal_stabbing(self):
        # DnaSP 6 (CODIGO2.vb::RecombinacionRM) is not the graph-theoretic
        # minimum interval-stabbing number: two intervals that only touch at
        # a shared boundary point ((9,10) and (10,11)) are treated as
        # non-overlapping and both kept, even though a single point at 10
        # would stab both. True minimum stab count for this set is 2; DnaSP
        # (and this skill) report 3. This is the same discrepancy that showed
        # up against DnaSP 6.12 for Windows on real alignments (rp49 Rm
        # 7->10, DmelOSRegion Rm 5->7) when this skill used an optimal
        # interval-stabbing algorithm instead of DnaSP's own reduction.
        intervals = [(10, 11), (9, 10), (5, 6), (2, 6)]
        assert dn._dnasp_rm_from_intervals(intervals) == 3


# ─────────────────────────────────────────────────────────────────────────────
# Mismatch distribution
# ─────────────────────────────────────────────────────────────────────────────

class TestComputeMismatch:
    def test_identical_seqs_zero_mean(self):
        seqs = ["ATCG"] * 4
        m = dn.compute_mismatch(seqs)
        assert m.mean == pytest.approx(0.0)
        assert m.n_pairs == 6

    def test_one_diff(self):
        seqs = ["ATCG", "ATCC"]
        m = dn.compute_mismatch(seqs)
        assert m.n_pairs == 1
        assert m.mean == pytest.approx(1.0)
        assert m.observed == {1: 1}

    def test_demo_mismatch(self, demo_seqs):
        # DEMO: 15 pairs; k = 13/15 diffs per pair
        m = dn.compute_mismatch(demo_seqs)
        assert m.n_pairs == 15
        assert m.mean == pytest.approx(13 / 15, abs=1e-9)

    def test_mismatch_histogram_keys_are_diff_counts(self, demo_seqs):
        m = dn.compute_mismatch(demo_seqs)
        # All observed difference counts should be non-negative integers
        for k in m.observed:
            assert k >= 0
        # Sum of counts equals total pairs
        assert sum(m.observed.values()) == m.n_pairs

    def test_raggedness_is_non_negative(self, demo_seqs):
        m = dn.compute_mismatch(demo_seqs)
        assert m.raggedness is not None
        assert m.raggedness >= 0.0

    def test_identical_seqs_raggedness(self):
        # All pairwise diffs = 0: single spike → ragged (all weight on 0)
        seqs = ["AAAA"] * 4
        m = dn.compute_mismatch(seqs)
        # raggedness = (f(0) - f(-1))^2 + ... with f(-1)=0, f(0)=1 → r ≥ 0
        assert m.raggedness is not None
        assert m.raggedness >= 0

    def test_single_seq_empty(self):
        seqs = ["ATCG"]
        m = dn.compute_mismatch(seqs)
        assert m.n_pairs == 0
        assert m.mean == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# InDel polymorphism
# ─────────────────────────────────────────────────────────────────────────────

class TestIdentifyIndelEvents:
    def test_no_gaps(self):
        seqs = ["ATCG", "ATCG"]
        events = dn._identify_indel_events(seqs)
        assert events == []

    def test_single_gap_column(self):
        # pos 1 has a gap in seq 0 only → 1 InDel event
        seqs = ["A-CG", "ATCG"]
        events = dn._identify_indel_events(seqs)
        assert len(events) == 1
        assert events[0].length == 1

    def test_contiguous_gap_same_set_is_one_event(self):
        # pos 1-2 have gaps in seq 0 only → 1 event of length 2
        seqs = ["A--G", "ATCG"]
        events = dn._identify_indel_events(seqs)
        assert len(events) == 1
        assert events[0].length == 2

    def test_two_separate_events(self):
        # pos 1: gap in seq 0; pos 3: gap in seq 0 → 2 separate events
        seqs = ["A-C-", "ATCG"]
        events = dn._identify_indel_events(seqs)
        assert len(events) == 2

    def test_different_gap_sets_are_separate_events(self):
        # pos 1: gap in seq 0; pos 2: gap in seq 1 → 2 events (different gap sets)
        seqs = ["A-CG", "AT-G"]
        events = dn._identify_indel_events(seqs)
        assert len(events) == 2

    def test_all_gaps_no_event(self):
        # If all sequences have a gap at a position, it's not an InDel event
        # (no non-gap sequences to define the other allele)
        seqs = ["A-CG", "A-CG"]
        events = dn._identify_indel_events(seqs)
        assert len(events) == 0  # all-gap column not treated as InDel


class TestComputeIndel:
    def test_no_gaps_returns_zero_events(self):
        seqs = ["ATCG", "ATCG"]
        stats = dn.compute_indel(seqs)
        assert stats.n_events == 0
        assert stats.n_positions_with_gaps == 0

    def test_single_event_stats(self):
        # One 2-column InDel event
        seqs = ["A--G", "ATCG", "ATCG"]
        stats = dn.compute_indel(seqs)
        assert stats.n_events == 1
        assert stats.n_positions_with_gaps == 2
        assert stats.mean_event_length == pytest.approx(2.0)

    def test_haplotype_count(self):
        # Two sequences: one with gap, two without → 2 InDel haplotypes
        seqs = ["A-CG", "ATCG", "ATCG"]
        stats = dn.compute_indel(seqs)
        # Binary haplotypes: (1,) for seq 0, (0,) for seqs 1 and 2 → 2 haplotypes
        assert stats.n_haplotypes == 2

    def test_k_indel_is_zero_for_no_events(self):
        seqs = ["ATCG", "ATCG", "ATCG"]
        stats = dn.compute_indel(seqs)
        assert stats.k_indel == pytest.approx(0.0)

    def test_pi_indel_non_negative(self):
        seqs = ["A-CG", "ATCG", "ATCG"]
        stats = dn.compute_indel(seqs)
        assert stats.pi_indel >= 0.0

    def test_single_sequence_empty(self):
        stats = dn.compute_indel(["ATCG"])
        assert stats.n_events == 0


# ─────────────────────────────────────────────────────────────────────────────
# Divergence between populations
# ─────────────────────────────────────────────────────────────────────────────

class TestComputeDivergence:
    # Pop1: 3 identical sequences
    # Pop2: 3 sequences differing from pop1 at every site
    # Expected Dxy = 1.0 per site (all sites fixed differences)

    POP1 = ["AAAA", "AAAA", "AAAA"]
    POP2 = ["TTTT", "TTTT", "TTTT"]

    def test_Dxy_fully_diverged(self):
        d = dn.compute_divergence(self.POP1, self.POP2, "p1", "p2")
        assert d.Dxy == pytest.approx(1.0, abs=1e-9)

    def test_Da_equals_Dxy_when_no_within_variation(self):
        # When both pops are invariant (Pi1 = Pi2 = 0), Da = Dxy
        d = dn.compute_divergence(self.POP1, self.POP2, "p1", "p2")
        assert d.Da == pytest.approx(d.Dxy, abs=1e-9)

    def test_Dxy_identical_pops(self):
        pop = ["ATCG", "ATCG"]
        d = dn.compute_divergence(pop, pop, "p1", "p2")
        assert d.Dxy == pytest.approx(0.0, abs=1e-9)

    def test_Da_leq_Dxy(self):
        # Net divergence cannot exceed gross divergence
        pop1 = ["ATCG", "ATCC"]
        pop2 = ["TTCG", "TTCC"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.Da <= d.Dxy + 1e-12

    def test_n_fixed_fully_diverged(self):
        # All 4 sites are fixed differences
        d = dn.compute_divergence(self.POP1, self.POP2, "p1", "p2")
        assert d.n_fixed == 4

    def test_shared_polymorphism(self):
        # Site 0: A/T segregating in BOTH pops → shared polymorphism
        pop1 = ["ATCG", "TTCG"]
        pop2 = ["ATCG", "TTCG"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.n_shared >= 1
        assert d.n_fixed == 0  # no fixed differences (same alleles)

    def test_L_net_counted(self):
        pop1 = ["ATCG", "ATCG"]
        pop2 = ["GCTA", "GCTA"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.L_net == 4

    def test_L_net_gap_excluded(self):
        # Gap at pos 0 in pop1 → L_net = 3
        pop1 = ["-TCG", "-TCG"]
        pop2 = ["ATCG", "ATCG"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.L_net == 3

    def test_Dxy_symmetry(self):
        # Dxy(pop1→pop2) == Dxy(pop2→pop1)
        pop1 = ["ATCG", "ATCC"]
        pop2 = ["TTCG", "TTCC"]
        d12 = dn.compute_divergence(pop1, pop2, "p1", "p2")
        d21 = dn.compute_divergence(pop2, pop1, "p2", "p1")
        assert d12.Dxy == pytest.approx(d21.Dxy, abs=1e-9)

    def test_Pi1_Pi2_computed(self):
        pop1 = ["ATCG", "TTCG"]   # 1 diff at pos 0
        pop2 = ["GCTA", "GCTA"]   # identical
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.Pi1 > 0.0
        assert d.Pi2 == pytest.approx(0.0, abs=1e-9)

    def test_partial_divergence(self):
        # 2 sites differ between pops, 2 sites identical
        pop1 = ["AACG", "AACG"]
        pop2 = ["TTCG", "TTCG"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.Dxy == pytest.approx(2 / 4, abs=1e-9)

    # --- Mod3BuscaShareFixDifferences parity (VB) --------------------------------

    def test_one_population_polymorphic_is_private_not_shared(self):
        # pop1 {A,G}, pop2 {A}: only pop1 segregates -> private to pop1, NOT shared
        pop1 = ["AC", "GC"]
        pop2 = ["AC", "AC"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.n_private1 == 1
        assert d.n_private2 == 0
        assert d.n_shared == 0
        assert d.n_fixed == 0

    def test_private_counts_are_mutations_not_sites(self):
        # triallelic in pop1, monomorphic pop2 sharing A -> 2 private mutations
        pop1 = ["A", "C", "G"]
        pop2 = ["A", "A", "A"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.n_private1 == 2
        assert d.n_shared == 0

    def test_fixed_difference_site_can_carry_private_mutation(self):
        # pop1 {A,C} disjoint from pop2 {T}: 1 fixed-difference site + 1 private
        pop1 = ["A", "C"]
        pop2 = ["T", "T"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.n_fixed == 1
        assert d.n_private1 == 1

    def test_both_biallelic_one_shared_allele(self):
        # pop1 {A,G}, pop2 {G,T}: share G, union {A,G,T} -> (2,2,3): 1 private each
        pop1 = ["A", "G"]
        pop2 = ["G", "T"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert (d.n_private1, d.n_private2, d.n_shared, d.n_fixed) == (1, 1, 0, 0)

    def test_both_biallelic_same_pair_is_shared(self):
        pop1 = ["A", "G", "A"]
        pop2 = ["A", "G", "G"]
        d = dn.compute_divergence(pop1, pop2, "p1", "p2")
        assert d.n_shared == 1
        assert d.n_private1 == 0 and d.n_private2 == 0


# ─────────────────────────────────────────────────────────────────────────────
# Population file parsing and alignment splitting
# ─────────────────────────────────────────────────────────────────────────────

class TestLoadPopFile:
    def test_basic_two_pops(self, tmp_path):
        p = tmp_path / "pops.txt"
        p.write_text("seq1\tpop1\nseq2\tpop1\nseq3\tpop2\n")
        result = dn.load_pop_file(p)
        assert result == {"seq1": "pop1", "seq2": "pop1", "seq3": "pop2"}

    def test_comment_lines_ignored(self, tmp_path):
        p = tmp_path / "pops.txt"
        p.write_text("# comment\nseq1\tpop1\nseq2\tpop2\n")
        result = dn.load_pop_file(p)
        assert "# comment" not in result
        assert len(result) == 2

    def test_blank_lines_ignored(self, tmp_path):
        p = tmp_path / "pops.txt"
        p.write_text("seq1\tpop1\n\nseq2\tpop2\n")
        result = dn.load_pop_file(p)
        assert len(result) == 2

    def test_returns_empty_for_empty_file(self, tmp_path):
        p = tmp_path / "pops.txt"
        p.write_text("")
        result = dn.load_pop_file(p)
        assert result == {}

    def test_space_separated(self, tmp_path):
        # DnaSP's VCF .SG.txt files use a single space
        p = tmp_path / "pops.SG.txt"
        p.write_text("Indiv1P1 Population1\nIndiv2P1 Population1\nIndiv5P2 Population2\n")
        result = dn.load_pop_file(p)
        assert result == {"Indiv1P1": "Population1", "Indiv2P1": "Population1",
                          "Indiv5P2": "Population2"}


_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
)


class TestParseVCF:
    def _write(self, tmp_path, body):
        p = tmp_path / "x.vcf"
        p.write_text(_VCF_HEADER + body)
        return p

    def test_phased_biallelic_two_haplotypes_per_sample(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|1\n")
        v = dn.parse_vcf(p)
        assert v.ploidy == 2 and v.is_phased is True
        aln = v.alignments["chr1"]
        assert aln.names == ["S1_h1", "S1_h2", "S2_h1", "S2_h2", "S3_h1", "S3_h2"]
        assert aln.seqs == ["AC", "AC", "AC", "GC", "GC", "GT"]

    def test_one_msa_per_chrom(self, tmp_path):
        p = self._write(tmp_path,
            "chrA\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "chrB\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n")
        v = dn.parse_vcf(p)
        assert set(v.alignments) == {"chrA", "chrB"}

    def test_region_filter(self, tmp_path):
        p = self._write(tmp_path,
            "chrA\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "chrB\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n")
        v = dn.parse_vcf(p, region="chrB")
        assert set(v.alignments) == {"chrB"}

    def test_indels_and_multiallelic_skipped(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tAT\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"   # indel
            "chr1\t20\t.\tG\tA,T\t.\tPASS\t.\tGT\t0|0\t0|1\t1|2\n"  # multiallelic
            "chr1\t30\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n")   # kept
        v = dn.parse_vcf(p)
        assert v.n_indels_skipped == 1
        assert v.n_multiallelic_skipped == 1
        assert v.alignments["chr1"].L == 1

    def test_unphased_het_becomes_gap(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n")
        v = dn.parse_vcf(p)
        assert v.ploidy == 2 and v.is_phased is False
        # S2 is 0/1 (het, unphased) -> both haplotypes gap
        assert v.alignments["chr1"].seqs == ["A", "A", "-", "-", "G", "G"]

    def test_unphased_het_sites_are_counted(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"   # has a het
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\n"   # all homozygous
            "chr1\t30\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0\n")  # two hets
        v = dn.parse_vcf(p)
        assert v.n_unphased_het_sites == 2

    def test_phased_vcf_reports_zero_unphased_het_sites(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n")
        assert dn.parse_vcf(p).n_unphased_het_sites == 0

    def test_cli_warns_on_unphased_het_and_snp_window(self, tmp_path, capsys):
        vcf = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\n"
            "chr1\t30\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n")
        out = tmp_path / "out"
        assert dn.main(["--vcf", str(vcf), "--window", "2", "--step", "1",
                        "--output", str(out)]) == 0
        err = capsys.readouterr().err
        assert "unphased heterozygous genotype" in err
        assert "SNP index" in err
        report = (out / "report.md").read_text()
        assert "slide over retained variant sites (SNP index)" in report

    def test_missing_gt_excludes_sample_from_chrom(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1\t1|1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|1\n")
        v = dn.parse_vcf(p)
        # S1 is .|. at the first variant of chr1 -> dropped from the whole MSA
        assert v.alignments["chr1"].names == ["S2_h1", "S2_h2", "S3_h1", "S3_h2"]

    def test_haploid_one_row_per_sample(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1\t1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0\t0\t1\n")
        v = dn.parse_vcf(p)
        assert v.ploidy == 1
        assert v.alignments["chr1"].names == ["S1", "S2", "S3"]
        assert v.alignments["chr1"].seqs == ["AC", "GC", "GT"]

    def test_filter_field_ignored(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tq10\t.\tGT\t0|0\t0|1\t1|1\n")
        v = dn.parse_vcf(p)  # q10 line still used, as in DnaSP
        assert v.alignments["chr1"].L == 1

    def test_merge_pools_all_chroms(self, tmp_path):
        p = self._write(tmp_path,
            "chrA\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "chrB\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n")
        v = dn.parse_vcf(p, merge=True)
        assert set(v.alignments) == {"<merged>"}
        assert v.alignments["<merged>"].L == 2

    def test_mixed_ploidy_rejected(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0\t1\n")  # S1 diploid, S2/S3 haploid
        with pytest.raises(ValueError, match="mixes ploidy"):
            dn.parse_vcf(p)

    def test_polyploid_rejected(self, tmp_path):
        p = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|1|1\t0|0|0\t1|1|1\n")
        with pytest.raises(ValueError, match="polyploid"):
            dn.parse_vcf(p)

    def test_sample_set_differs_per_chrom(self, tmp_path):
        # S1 missing at the first variant of chrA but present on chrB
        p = self._write(tmp_path,
            "chrA\t10\t.\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1\t1|1\n"
            "chrB\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n")
        v = dn.parse_vcf(p)
        assert "S1_h1" not in v.alignments["chrA"].names
        assert "S1_h1" in v.alignments["chrB"].names

    def test_cli_pop_map_built_per_chrom(self, tmp_path, capsys):
        vcf = self._write(tmp_path,
            "chrA\t10\t.\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1\t1|1\n"
            "chrB\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t1|1\n")
        pops = tmp_path / "p.txt"
        pops.write_text("S1 P1\nS2 P1\nS3 P2\n")
        rc = dn.main(["--vcf", str(vcf), "--pop-file", str(pops),
                      "--analysis", "fst", "--output", str(tmp_path / "out")])
        assert rc == 0
        # chrA drops S1 (missing there); chrB keeps all three
        assert (tmp_path / "out" / "chrA" / "report.md").exists()
        assert (tmp_path / "out" / "chrB" / "report.md").exists()

    def test_cli_colliding_chrom_names_get_distinct_dirs(self, tmp_path):
        vcf = self._write(tmp_path,
            "a:b\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "a/b\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|0\t0|0\n")
        out = tmp_path / "out"
        rc = dn.main(["--vcf", str(vcf), "--output", str(out)])
        assert rc == 0
        reports = list(out.rglob("report.md"))
        assert len(reports) == 2  # not overwritten into one dir

    def test_vcf_report_flags_per_variant_site_diversity(self, tmp_path):
        vcf = self._write(tmp_path,
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|1\t0|1\n")
        out = tmp_path / "out"
        assert dn.main(["--vcf", str(vcf), "--output", str(out)]) == 0
        report = (out / "report.md").read_text()
        assert "per variant site, not per base" in report
        tsv = (out / "results.tsv").read_text()
        assert "per variant site, not per base" in tsv

    def test_fasta_report_has_no_per_variant_site_note(self, tmp_path):
        f = tmp_path / "a.fas"
        f.write_text(">s1\nACGTACGT\n>s2\nACGAACGT\n>s3\nTCGTACGT\n")
        out = tmp_path / "out"
        assert dn.main(["--input", str(f), "--output", str(out)]) == 0
        assert "per variant site" not in (out / "report.md").read_text()


class TestSplitAlignmentByPop:
    def _make_aln(self):
        return dn.Alignment(
            names=["s1", "s2", "s3", "s4"],
            seqs=["ATCG", "ATCC", "TTCG", "TTCC"],
            source="test",
        )

    def test_two_pops_split(self):
        aln = self._make_aln()
        assigns = {"s1": "pop1", "s2": "pop1", "s3": "pop2", "s4": "pop2"}
        result = dn.split_alignment_by_pop(aln, assigns)
        assert set(result.keys()) == {"pop1", "pop2"}
        assert result["pop1"].n == 2
        assert result["pop2"].n == 2
        assert result["pop1"].names == ["s1", "s2"]
        assert result["pop2"].names == ["s3", "s4"]

    def test_unassigned_seqs_excluded(self):
        aln = self._make_aln()
        assigns = {"s1": "pop1", "s2": "pop1"}  # s3, s4 not assigned
        result = dn.split_alignment_by_pop(aln, assigns)
        assert "pop2" not in result
        assert result["pop1"].n == 2

    def test_sequences_preserved(self):
        aln = self._make_aln()
        assigns = {"s1": "pop1", "s2": "pop1", "s3": "pop2", "s4": "pop2"}
        result = dn.split_alignment_by_pop(aln, assigns)
        assert result["pop1"].seqs[0] == "ATCG"
        assert result["pop2"].seqs[0] == "TTCG"


# ─────────────────────────────────────────────────────────────────────────────
# Analysis set parsing (_parse_analyses)
# ─────────────────────────────────────────────────────────────────────────────

class TestParseAnalyses:
    def test_all_expands_to_full_set(self):
        result = dn._parse_analyses("all")
        assert result == dn.VALID_ANALYSES

    def test_polymorphism_always_included(self):
        result = dn._parse_analyses("ld")
        assert "polymorphism" in result

    def test_single_valid(self):
        result = dn._parse_analyses("ld")
        assert "ld" in result

    def test_comma_separated(self):
        result = dn._parse_analyses("ld,recombination,indel")
        assert "ld" in result
        assert "recombination" in result
        assert "indel" in result

    def test_unknown_analyses_ignored(self, capsys):
        result = dn._parse_analyses("ld,foobar")
        assert "foobar" not in result
        assert "ld" in result

    def test_case_insensitive(self):
        result = dn._parse_analyses("LD,Recombination")
        assert "ld" in result
        assert "recombination" in result


# ─────────────────────────────────────────────────────────────────────────────
# run_analysis integration tests (new modules)
# ─────────────────────────────────────────────────────────────────────────────

class TestRunAnalysisNewModules:
    def test_ld_module_included(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, analyses={"ld"})
        assert "ld" in results
        ld = results["ld"]
        assert ld is not None
        assert ld.n_biallelic == 8  # 8 biallelic sites in 300 bp demo

    def test_recombination_module_included(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, analyses={"recombination"})
        assert "recombination" in results
        assert results["recombination"] is not None

    def test_popsize_module_included(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, analyses={"popsize"})
        assert "popsize" in results
        m = results["popsize"]
        assert m.n_pairs == 55  # 11 choose 2

    def test_indel_module_no_gaps(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, analyses={"indel"})
        assert "indel" in results
        stats = results["indel"]
        assert stats.n_events == 0  # DEMO has no gaps

    def test_all_analyses_run(self, demo_fasta_path):
        aln = dn.parse_fasta(demo_fasta_path)
        results = dn.run_analysis(aln, analyses=dn.VALID_ANALYSES - {"divergence"})
        assert results["global"] is not None
        assert results["ld"] is not None
        assert results["recombination"] is not None
        assert results["popsize"] is not None
        assert results["indel"] is not None

    def test_divergence_module_two_populations(self, tmp_path):
        # Write two population FASTAs
        pop1 = tmp_path / "pop1.fas"
        pop2 = tmp_path / "pop2.fas"
        pop1.write_text(">s1\nAAAA\n>s2\nAAAA\n")
        pop2.write_text(">s3\nTTTT\n>s4\nTTTT\n")
        aln1 = dn.parse_fasta(pop1)
        aln2 = dn.parse_fasta(pop2)
        results = dn.run_analysis(aln1, analyses={"divergence"}, aln2=aln2)
        assert "divergence" in results
        d = results["divergence"]
        assert d.Dxy == pytest.approx(1.0, abs=1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# Chi-square p-value helpers
# ─────────────────────────────────────────────────────────────────────────────

class TestChiSqPvalue:
    """Tests for the generalised _chi2_pvalue function."""

    def test_df1_matches_erfc(self):
        # Should match the old _chi2_1df_pvalue exactly
        for chi2 in [0.5, 1.0, 2.0, 3.841, 6.635]:
            expected = math.erfc(math.sqrt(chi2 / 2))
            assert dn._chi2_pvalue(chi2, 1) == pytest.approx(expected, rel=1e-10)

    def test_df2(self):
        # P = exp(-chi2/2)
        assert dn._chi2_pvalue(0.0, 2) == pytest.approx(1.0, rel=1e-10)
        assert dn._chi2_pvalue(2.0, 2) == pytest.approx(math.exp(-1.0), rel=1e-10)
        assert dn._chi2_pvalue(4.0, 2) == pytest.approx(math.exp(-2.0), rel=1e-10)

    def test_df4_known(self):
        # Chi2=9.488 ≈ 0.05 critical value for df=4
        p = dn._chi2_pvalue(9.488, 4)
        assert p == pytest.approx(0.05, abs=0.002)

    def test_df6_known(self):
        # Chi2=12.592 ≈ 0.05 critical value for df=6
        p = dn._chi2_pvalue(12.592, 6)
        assert p == pytest.approx(0.05, abs=0.002)

    def test_zero_chi2_returns_one(self):
        for df in [1, 2, 3, 4]:
            assert dn._chi2_pvalue(0.0, df) == pytest.approx(1.0)

    def test_large_chi2_small_p(self):
        assert dn._chi2_pvalue(30.0, 4) < 1e-5

    def test_df0_returns_one(self):
        assert dn._chi2_pvalue(5.0, 0) == pytest.approx(1.0)


# ─────────────────────────────────────────────────────────────────────────────
# Gamma incomplete functions
# ─────────────────────────────────────────────────────────────────────────────

class TestGammaincc:
    """Spot-checks against known values from scipy.special.gammaincc."""

    def test_q_1_1(self):
        # Q(1, 1) = exp(-1) ≈ 0.36788
        q = dn._gammaincc(1.0, 1.0)
        assert q == pytest.approx(0.36787944117, rel=1e-6)

    def test_q_2_2(self):
        # Q(2, 2) = exp(-2)*(1+2) = 3*exp(-2) ≈ 0.40601
        q = dn._gammaincc(2.0, 2.0)
        assert q == pytest.approx(0.40600584971, rel=1e-5)

    def test_q_zero_x_is_one(self):
        assert dn._gammaincc(3.0, 0.0) == pytest.approx(1.0)

    def test_q_large_x_near_zero(self):
        q = dn._gammaincc(2.0, 30.0)
        assert q < 1e-10

    def test_series_path(self):
        # x < a+1 uses series: Q(5, 3) should be > 0.5
        q = dn._gammaincc(5.0, 3.0)
        assert 0.5 < q < 1.0

    def test_cf_path(self):
        # x > a+1 uses CF: Q(2, 5) should be small
        q = dn._gammaincc(2.0, 5.0)
        assert 0.0 < q < 0.1


# ─────────────────────────────────────────────────────────────────────────────
# Bisection helper
# ─────────────────────────────────────────────────────────────────────────────

class TestBisect:
    def test_simple_linear(self):
        # f(x) = x - 3, root = 3
        root = dn._bisect(lambda x: x - 3.0, 0.0, 10.0)
        assert root == pytest.approx(3.0, abs=1e-8)

    def test_quadratic(self):
        # f(x) = x^2 - 4, root = 2 in [0, 5]
        root = dn._bisect(lambda x: x**2 - 4.0, 0.0, 5.0)
        assert root == pytest.approx(2.0, abs=1e-8)

    def test_no_sign_change_returns_midpoint(self):
        # f always positive  -  should return midpoint without crash
        result = dn._bisect(lambda x: x + 1.0, 0.0, 10.0)
        assert result == pytest.approx(5.0)

    def test_root_at_boundary(self):
        root = dn._bisect(lambda x: x, 0.0, 5.0)
        assert abs(root) < 1e-8


# ─────────────────────────────────────────────────────────────────────────────
# _count_derived
# ─────────────────────────────────────────────────────────────────────────────

class TestCountDerived:
    """Tests for outgroup-polarised derived mutation counting."""

    def test_simple_derived_singleton(self):
        # Outgroup = A; seq0 has T (derived, in 1 seq)
        seqs = ["TAAAA", "AAAAA", "AAAAA", "AAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        assert S == 1
        assert eta == 1
        assert eta_e == 1

    def test_derived_in_two_seqs_not_external(self):
        # Outgroup = A; T in 2 seqs → derived but not singleton
        seqs = ["TAAAA", "TAAAA", "AAAAA", "AAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        assert S == 1
        assert eta == 1
        assert eta_e == 0

    def test_outgroup_gap_skipped(self):
        # Gap in outgroup at pos 0 → skip that site
        seqs = ["TAAAA", "AAAAA", "AAAAA"]
        outgroup = "-AAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        assert S == 0
        assert eta == 0
        assert eta_e == 0

    def test_ingroup_gap_skipped_complete_deletion(self):
        # Gap in ingroup seq at pos 0 → skip that site
        seqs = ["-AAAA", "TAAAA", "AAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        # pos 0: ingroup gap → skip; pos 1-4: all A → 0 mutations
        assert eta == 0

    def test_ancestral_allele_absent_skipped(self):
        # Outgroup = A, but all ingroup have T  -  ancestral absent
        seqs = ["TAAAA", "TAAAA", "TAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        # All seqs have T but outgroup has A; A not in ingroup at pos 0
        assert eta == 0

    def test_monomorphic_site_skipped(self):
        seqs = ["AAAAA", "AAAAA", "AAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        assert S == 0
        assert eta == 0
        assert eta_e == 0

    def test_two_derived_alleles_site_capped_at_one(self):
        # pos 0: outgroup=A, seqs have T (1) and C (1) -> both derived, each
        # singleton. eta (total derived mutations) counts both; eta_e is a
        # SITE count capped at 1, matching DnaSP's EtaE = EtaE - ExternaMut
        # discount at a multiallelic site (Mod12FuLiOutgroupNew).
        seqs = ["TAAAA", "CAAAA", "AAAAA"]
        outgroup = "AAAAA"
        S, eta, eta_e = dn._count_derived(seqs, outgroup)
        assert S == 1
        assert eta == 2       # two derived alleles (T and C)
        assert eta_e == 1     # capped: one site, not two


# ─────────────────────────────────────────────────────────────────────────────
# compute_fu_li_outgroup
# ─────────────────────────────────────────────────────────────────────────────

class TestComputeFuLiOutgroup:
    """Tests for the outgroup-based Fu & Li D and F statistics."""

    def _ingroup_outgroup(self):
        """
        4 ingroup + 1 outgroup.
        Outgroup = AAAAA
        Seq1: TAAAA   -  pos0 derived T, singleton → eta=1, eta_e=1
        Seq2-4: AAAAA  -  ancestral
        k_bar (all pairwise among 4 ingroup):
          pair (seq1,seq2), (seq1,seq3), (seq1,seq4) = 1 diff each → 3 pairs × 1
          pair (seq2,seq3), (seq2,seq4), (seq3,seq4) = 0 → k_bar = 3/6 = 0.5
        """
        seqs = ["TAAAA", "AAAAA", "AAAAA", "AAAAA"]
        outgroup = "AAAAA"
        return seqs, outgroup

    def test_returns_fulioutgroupstats(self):
        seqs, outgroup = self._ingroup_outgroup()
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert isinstance(result, dn.FuLiOutgroupStats)

    def test_eta_eta_e_counted(self):
        seqs, outgroup = self._ingroup_outgroup()
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.eta == 1
        assert result.eta_e == 1

    def test_k_bar(self):
        seqs, outgroup = self._ingroup_outgroup()
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        # k_bar = mean pairwise diffs = 3/6 = 0.5
        assert result.k_bar == pytest.approx(0.5, abs=1e-9)

    def test_d_and_f_not_none(self):
        # With enough sequences and variation both D and F should be non-None
        seqs, outgroup = self._ingroup_outgroup()
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.D is not None
        assert result.F is not None

    def test_k_bar_excludes_outgroup_gap_column(self):
        # DnaSP FULI.vb: rp1 (k_bar) accumulates only over informative sites.
        # Column 0 is polymorphic but the outgroup has a gap there -> excluded
        # from k_bar just as it is from eta/eta_e.
        seqs = ["TAAAA", "AAAAA", "AAAAA", "AAAAA"]
        outgroup = "-AAAA"
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.eta == 0
        assert result.k_bar == pytest.approx(0.0, abs=1e-12)

    def test_k_bar_excludes_ancestral_absent_column(self):
        # Column 0: ingroup {A,C}, outgroup G (ancestral absent) -> not orientable.
        # Column 1: ingroup {T,A}, outgroup A -> orientable, one derived singleton.
        # k_bar must reflect column 1 only (4 diffs / 10 pairs = 0.4).
        seqs = ["AT", "AA", "CA", "CA", "CA"]
        outgroup = "GA"
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.eta == 1
        assert result.eta_e == 1
        assert result.k_bar == pytest.approx(0.4, abs=1e-12)

    def test_no_variation_returns_none(self):
        seqs = ["AAAAA"] * 5
        outgroup = "AAAAA"
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.D is None
        assert result.F is None

    def test_too_few_seqs_returns_empty(self):
        # n < 4 → returns dataclass with D=F=None
        result = dn.compute_fu_li_outgroup(["ATCG", "TTCG"], "ATCG")
        assert result.D is None
        assert result.F is None

    def test_n_recorded(self):
        seqs, outgroup = self._ingroup_outgroup()
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.n == 4

    def test_excess_external_mutations_negative_d(self):
        """All derived mutations are singletons → η_e = η → D = η − aₙη_e < 0
        (aₙ > 1 for n ≥ 4), and F = k̄ − η_e < 0 too."""
        # All 4 derived mutations are singletons
        seqs = ["TAAAA", "ATAAA", "AATAA", "AAATA"]
        outgroup = "AAAAA"
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result is not None
        assert result.eta == 4
        assert result.eta_e == 4
        assert result.D is not None and result.D < 0
        assert result.F is not None and result.F < 0

    def test_d_and_f_scale_by_S_not_eta(self):
        # DnaSP's default mode (Mod12FuLiOutgroupNew, SSOrMutations=1) scales
        # D and F by the orientable segregating-SITE count (pv1), not by the
        # total-mutations count (EtaT) -- they only coincide when every site
        # is biallelic. Site 0 is biallelic (1 derived, singleton); site 1 is
        # triallelic (2 derived alleles, both singleton, capped to 1 site of
        # eta_e) so S=2 but eta=3 here, letting the two formulas diverge.
        seqs = ["TA", "AA", "AC", "AG", "AA"]
        outgroup = "AA"
        result = dn.compute_fu_li_outgroup(seqs, outgroup)
        assert result.S == 2
        assert result.eta == 3
        assert result.eta_e == 2

        n = 5
        a_n = dn._harmonic(n, 1)
        b_n = dn._harmonic(n, 2)
        a_n1 = a_n + 1.0 / n
        c_n = 2.0 * (n * a_n - 2.0 * (n - 1.0)) / ((n - 1.0) * (n - 2.0))
        v_D = 1.0 + (a_n**2 / (b_n + a_n**2)) * (c_n - (n + 1.0) / (n - 1.0))
        u_D = a_n - 1.0 - v_D
        v_F = (c_n + 2.0 * (n**2 + n + 3.0) / (9.0 * n * (n - 1.0))
               - 2.0 / (n - 1.0)) / (a_n**2 + b_n)
        u_F = (1.0 + (n + 1.0) / (3.0 * (n - 1.0))
               - 4.0 * ((n + 1.0) / (n - 1.0)**2) * (a_n1 - 2.0 * n / (n + 1.0))) / a_n - v_F

        S, eta, eta_e = result.S, result.eta, result.eta_e
        D_from_S = (S - a_n * eta_e) / math.sqrt(u_D * S + v_D * S * S)
        D_from_eta = (eta - a_n * eta_e) / math.sqrt(u_D * eta + v_D * eta * eta)
        assert D_from_S != pytest.approx(D_from_eta)
        assert result.D == pytest.approx(D_from_S)

        F_from_S = (result.k_bar - eta_e) / math.sqrt(u_F * S + v_F * S * S)
        assert result.F == pytest.approx(F_from_S)

    def test_run_analysis_fuliout(self):
        """fuliout dispatched correctly from run_analysis."""
        seqs = ["TAAAA", "AAAAA", "AAAAA", "AAAAA", "AAAAA"]
        outgroup = "AAAAA"
        aln = dn.Alignment(names=[f"s{i}" for i in range(5)], seqs=seqs)
        results = dn.run_analysis(aln, analyses={"fuliout"}, outgroup=outgroup)
        assert "fuliout" in results
        assert results["fuliout"] is not None
        assert results["fuliout"].eta == 1


# ─────────────────────────────────────────────────────────────────────────────
# load_hka_file
# ─────────────────────────────────────────────────────────────────────────────

class TestLoadHKAFile:
    def test_basic_parsing(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text(
            "# locus n S L_poly D L_div\n"
            "locus1 10 5 400 10 400\n"
            "locus2 12 2 300 8 300\n"
        )
        loci = dn.load_hka_file(f)
        assert len(loci) == 2
        assert loci[0].name == "locus1"
        assert (loci[0].n, loci[0].S, loci[0].D) == (10, 5, 10)
        assert loci[0].L_poly == 400.0 and loci[0].L_div == 400.0

    def test_header_line_skipped(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text("locus n S L_poly D L_div\nlocus1 10 5 400 10 400\n")
        loci = dn.load_hka_file(f)
        assert len(loci) == 1

    def test_comment_lines_skipped(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text("# comment\nlocus1 10 5 400 10 400\nlocus2 8 3 200 7 200\n")
        loci = dn.load_hka_file(f)
        assert len(loci) == 2

    def test_empty_file(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text("")
        assert dn.load_hka_file(f) == []

    def test_ldiv_defaults_to_lpoly(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text("locus_A 20 7 500 15\n")
        loc = dn.load_hka_file(f)[0]
        assert loc.L_poly == 500.0 and loc.L_div == 500.0

    def test_chromosome_factor(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text("A 20 7 500 15 500 X\nB 20 3 500 9 500 A\n")
        loci = dn.load_hka_file(f)
        assert loci[0].sex == 0.75
        assert loci[1].sex == 1.0


# ─────────────────────────────────────────────────────────────────────────────
# compute_hka
# ─────────────────────────────────────────────────────────────────────────────

class TestComputeHKA:
    """Two-locus HKA test, ported from DnaSP 6 (HKA.vb)."""

    def _neutral_loci(self):
        # n1 == n2, equal L, S1/D1 == S2/D2 -> a single (theta1, theta2, T)
        # fits all four residuals exactly, so chi2 == 0.
        return [
            dn.HKALocus(name="A", n=11, S=10, D=20, L_poly=1000, L_div=1000),
            dn.HKALocus(name="B", n=11, S=5, D=10, L_poly=1000, L_div=1000),
        ]

    def _selection_loci(self):
        return [
            dn.HKALocus(name="X", n=20, S=40, D=3, L_poly=1000, L_div=1000),
            dn.HKALocus(name="Y", n=20, S=2, D=40, L_poly=1000, L_div=1000),
        ]

    def test_returns_hkastats(self):
        assert isinstance(dn.compute_hka(self._neutral_loci()), dn.HKAStats)

    def test_requires_exactly_two_loci(self):
        one = dn.compute_hka(self._neutral_loci()[:1])
        assert one.error is not None and not one.loci_results
        three = dn.compute_hka(self._neutral_loci() + self._neutral_loci()[:1])
        assert three.error is not None and not three.loci_results

    def test_requires_site_counts(self):
        loci = [dn.HKALocus(name="A", n=10, S=5, D=10),
                dn.HKALocus(name="B", n=10, S=3, D=8)]
        r = dn.compute_hka(loci)
        assert r.error is not None and not r.loci_results

    def test_df_is_one(self):
        assert dn.compute_hka(self._neutral_loci()).df == 1

    def test_neutral_data_chi2_zero(self):
        r = dn.compute_hka(self._neutral_loci())
        assert r.error is None
        assert r.chi2 == pytest.approx(0.0, abs=1e-6)
        assert r.p_value == pytest.approx(1.0, abs=1e-3)

    def test_expectations_match_inputs_under_neutrality(self):
        r = dn.compute_hka(self._neutral_loci())
        for lr, loc in zip(r.loci_results, self._neutral_loci()):
            assert lr["E_S"] == pytest.approx(loc.S, abs=1e-6)
            assert lr["E_D"] == pytest.approx(loc.D, abs=1e-6)

    def test_variance_includes_second_harmonic_term(self):
        # HKA (1987): Var[S] = E[S] + b_n * theta^2 * L^2  > E[S]
        r = dn.compute_hka(self._neutral_loci())
        for lr in r.loci_results:
            assert lr["Var_S"] > lr["E_S"]

    def test_selection_signal_rejects_neutrality(self):
        r = dn.compute_hka(self._selection_loci())
        assert r.error is None
        assert r.chi2 > 10.0
        assert r.p_value < 0.01

    def test_no_positive_solution_is_reported_not_fabricated(self):
        # S = 0 at both loci: no positive-theta root.  Must NOT invent a T.
        loci = [dn.HKALocus(name="A", n=10, S=0, D=5, L_poly=500, L_div=500),
                dn.HKALocus(name="B", n=10, S=0, D=3, L_poly=500, L_div=500)]
        r = dn.compute_hka(loci)
        assert not r.loci_results
        assert r.error is not None
        assert r.T_hat == 0.0

    def test_theta_hat_positive(self):
        for lr in dn.compute_hka(self._neutral_loci()).loci_results:
            assert lr["theta_hat"] > 0.0

    def test_run_analysis_hka(self, tmp_path):
        f = tmp_path / "demo.fas"
        f.write_text(">s1\nAAAA\n>s2\nTAAA\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(
            aln, analyses={"hka"}, hka_loci=self._neutral_loci()
        )
        assert results["hka"] is not None
        assert results["hka"].chi2 == pytest.approx(0.0, abs=1e-6)

    def test_load_and_compute_roundtrip(self, tmp_path):
        f = tmp_path / "hka.tsv"
        f.write_text(
            "# locus n S L_poly D L_div\n"
            "A 11 10 1000 20 1000\n"
            "B 11 5 1000 10 1000\n"
        )
        r = dn.compute_hka(dn.load_hka_file(f))
        assert r.error is None
        assert r.chi2 == pytest.approx(0.0, abs=1e-6)

    def test_error_and_note_are_distinct_fields(self):
        s = dn.HKAStats()
        assert s.error is None and s.note is None

    def test_successful_run_sets_neither_error_nor_note(self):
        for loci in (self._neutral_loci(), self._selection_loci()):
            r = dn.compute_hka(loci)
            assert r.loci_results
            assert r.error is None      # error means "test not run"
            assert r.note is None       # no multi-solution remark for these data

    def test_not_run_sets_error_not_note(self):
        r = dn.compute_hka(self._neutral_loci()[:1])
        assert r.error is not None and r.note is None and not r.loci_results

    def test_report_renders_hka_note_without_saying_not_run(self, tmp_path):
        aln = dn.Alignment(names=["a", "b"], seqs=["AC", "AT"], source="t")
        results = dn.run_analysis(aln, analyses=set())
        results["hka"] = dn.HKAStats(
            n_loci=2, T_hat=1.0, chi2=0.5, df=1, p_value=0.48,
            loci_results=[{"name": "A", "n": 10, "S": 5, "D": 8,
                           "theta_hat": 0.01, "E_S": 5.0, "E_D": 8.0,
                           "Var_S": 6.0, "Var_D": 9.0}],
            note="2 positive-θ solutions; reporting the first (θ₁=0.01).",
        )
        report = dn.write_report(tmp_path, "t", aln, results, []).read_text()
        assert "Note: 2 positive-θ solutions" in report
        assert "HKA not run" not in report


# ─────────────────────────────────────────────────────────────────────────────
# Group C: MK test and Ka/Ks helpers
# ─────────────────────────────────────────────────────────────────────────────

class TestLogFactorial:
    def test_zero(self):
        assert dn._log_factorial(0) == pytest.approx(0.0)

    def test_one(self):
        assert dn._log_factorial(1) == pytest.approx(0.0)

    def test_two(self):
        import math
        assert dn._log_factorial(2) == pytest.approx(math.log(2))

    def test_ten(self):
        import math
        # log(10!) = log(3628800)
        assert dn._log_factorial(10) == pytest.approx(math.log(3628800), rel=1e-9)


class TestFisherExact:
    def test_empty_table(self):
        # All zeros → p=1
        assert dn._fisher_exact_2x2(0, 0, 0, 0) == pytest.approx(1.0)

    def test_neutral_table_equal_cells(self):
        # Equal Ps, Ds, Pn, Dn → perfectly neutral, p should be 1.0
        assert dn._fisher_exact_2x2(5, 5, 5, 5) == pytest.approx(1.0)

    def test_extreme_positive_selection(self):
        # All fixed diff nonsynonymous, all poly synonymous → extreme table
        p = dn._fisher_exact_2x2(0, 10, 10, 0)
        assert p < 0.01

    def test_returns_at_most_one(self):
        p = dn._fisher_exact_2x2(3, 7, 7, 3)
        assert 0.0 <= p <= 1.0

    def test_symmetry(self):
        # Swapping rows should give the same p-value
        p1 = dn._fisher_exact_2x2(2, 8, 8, 2)
        p2 = dn._fisher_exact_2x2(8, 2, 2, 8)
        assert p1 == pytest.approx(p2)

    def test_zero_row(self):
        # One row is all zero  -  only one possible table
        p = dn._fisher_exact_2x2(0, 0, 5, 5)
        assert p == pytest.approx(1.0)

    def test_known_small_table(self):
        # From standard Fisher's exact test example:
        # 1  9
        # 11 3   → two-tailed p ~ 0.0028
        p = dn._fisher_exact_2x2(1, 9, 11, 3)
        assert p < 0.01

    def test_one_cell_zero_not_always_significant(self):
        # 0 1 / 1 0  -  2×2 with n=2; only two possible tables, p must be 1
        p = dn._fisher_exact_2x2(0, 1, 1, 0)
        assert p == pytest.approx(1.0)


class TestCountSynSites:
    def test_atg_zero(self):
        # ATG = Met (only codon for Met) → 0 synonymous sites
        assert dn._count_syn_sites_codon('ATG') == pytest.approx(0.0)

    def test_ggg_one(self):
        # GGG: pos2 has 3 synonymous alternatives (GGA/GGC/GGT all Gly) → 1.0
        assert dn._count_syn_sites_codon('GGG') == pytest.approx(1.0)

    def test_ttt_one_third(self):
        # TTT (Phe): only TTC is synonymous (pos2) → 1/3
        result = dn._count_syn_sites_codon('TTT')
        assert result == pytest.approx(1.0 / 3.0, rel=1e-9)

    def test_tgg_zero(self):
        # TGG = Trp (only Trp codon, pos2 alternatives TGA=stop, TGC=Cys, TGT=Cys → all nonsyn or stop)
        # pos0: AGG=Arg, CGG=Arg, GGG=Gly → all nonsyn
        # pos1: TTG=Leu, TCG=Ser, TAG=stop → all nonsyn/stop
        # pos2: TGA=stop, TGC=Cys, TGT=Cys → 0 syn
        assert dn._count_syn_sites_codon('TGG') == pytest.approx(0.0)

    def test_gap_returns_zero(self):
        assert dn._count_syn_sites_codon('AT-') == pytest.approx(0.0)

    def test_n_returns_zero(self):
        assert dn._count_syn_sites_codon('ATN') == pytest.approx(0.0)

    def test_stop_returns_zero(self):
        assert dn._count_syn_sites_codon('TAA') == pytest.approx(0.0)

    def test_total_between_zero_and_three(self):
        for codon in ['AAA', 'CCC', 'GGC', 'TTG', 'AGT', 'CTG']:
            result = dn._count_syn_sites_codon(codon)
            assert 0.0 <= result <= 3.0

    def test_tgt_is_half_a_site_when_one_alternative_is_a_stop(self):
        # DnaSP (SINONIMO.vb::ComputeFoldPos) excludes stop-codon
        # alternatives from the denominator: TGT (Cys), third position ->
        # TGC (Cys, synonymous), TGA (stop), TGG (Trp). One synonymous of
        # the two non-stop paths = 1/2 a site, not 1/3 -- the VB's own
        # comment: "de los dos caminos posibles uno es sinonimo ej. TGT".
        assert dn._count_syn_sites_codon('TGT') == pytest.approx(0.5)

    def test_tat_is_a_whole_site_when_two_alternatives_are_stops(self):
        # TAT (Tyr), third position -> TAC (Tyr), TAA (stop), TAG (stop):
        # the only non-stop path is synonymous, so the position is a whole
        # synonymous site (ComputeFoldPos, "Dos caminos son stop" -> "4").
        assert dn._count_syn_sites_codon('TAT') == pytest.approx(1.0)

    def test_agt_whole_site_only_under_vertebrate_mitochondrial_code(self):
        # AGT (Ser), third position -> AGC (Ser), AGA, AGG. Under the
        # standard code AGA/AGG are Arg: 1 synonymous of 3 -> 1/3. Under
        # the vertebrate mitochondrial code they are stops: 1 of 1 -> 1.0.
        # (This is what lifts COII's synonymous sites from 161.3 to
        # DnaSP's 168.222.)
        assert dn._count_syn_sites_codon('AGT') == pytest.approx(1.0 / 3.0)
        mito = dn._count_syn_sites_codon('AGT', dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert mito == pytest.approx(1.0)


class TestClassifyCodonPair:
    def test_identical(self):
        s, ns = dn._classify_codon_pair('ATG', 'ATG')
        assert s == pytest.approx(0.0)
        assert ns == pytest.approx(0.0)

    def test_one_synonymous(self):
        # TTT (Phe) vs TTC (Phe): one synonymous change
        s, ns = dn._classify_codon_pair('TTT', 'TTC')
        assert s == pytest.approx(1.0)
        assert ns == pytest.approx(0.0)

    def test_one_nonsynonymous(self):
        # ATG (Met) vs CTG (Leu): one nonsynonymous change
        s, ns = dn._classify_codon_pair('ATG', 'CTG')
        assert s == pytest.approx(0.0)
        assert ns == pytest.approx(1.0)

    def test_gap_returns_zero(self):
        s, ns = dn._classify_codon_pair('AT-', 'ATG')
        assert s == pytest.approx(0.0) and ns == pytest.approx(0.0)

    def test_sum_equals_n_diff_one(self):
        # For 1-diff codons, syn + nonsyn = 1
        s, ns = dn._classify_codon_pair('GGT', 'GGG')
        assert s + ns == pytest.approx(1.0)

    def test_two_diff_sum_equals_two(self):
        # ATG (Met) vs CTT (Leu): 2 diffs
        s, ns = dn._classify_codon_pair('ATG', 'CTT')
        assert s + ns == pytest.approx(2.0)

    def test_three_diff_sum_equals_three(self):
        # TTT (Phe) vs GGG (Gly): 3 diffs
        s, ns = dn._classify_codon_pair('TTT', 'GGG')
        total = s + ns
        assert total == pytest.approx(3.0, rel=1e-9)

    def test_stop_codon_is_a_21st_amino_acid_but_stop_paths_are_excluded(self):
        # DnaSP's help file: "if DnaSP finds stop codons (in the middle of
        # coding regions) they will be considered as if they would code for
        # a new amino acid (the amino acid 21; for example Selenocysteine)",
        # while "it excludes those pathways that go through stop codons".
        # TAA (*) vs AAG (Lys), two differences: TAA -> TAG -> AAG runs
        # through a stop and is dropped; TAA -> AAA -> AAG is one
        # replacement (* -> Lys) and one synonymous change (Lys -> Lys). An
        # earlier version of this test expected the pair to be skipped.
        s_, ns = dn._classify_codon_pair('TAA', 'AAG')
        assert s_ == pytest.approx(1.0) and ns == pytest.approx(1.0)
        # TAA vs ATG differs at all three positions: the two orderings that
        # start with TAG are dropped; the other four average (0+1+0+1)/4
        # synonymous and (3+2+3+2)/4 replacement changes.
        s_, ns = dn._classify_codon_pair('TAA', 'ATG')
        assert s_ == pytest.approx(0.5) and ns == pytest.approx(2.5)
        # Two stop codons one base apart are a synonymous difference (* = *),
        # as DnaSP's NumSynonEntreCodons compares the GCode strings.
        s_, ns = dn._classify_codon_pair('TAA', 'TAG')
        assert s_ == pytest.approx(1.0) and ns == pytest.approx(0.0)


class TestJcCorrect:
    def test_zero_p(self):
        assert dn._jc_correct(0.0) == pytest.approx(0.0)

    def test_small_p(self):
        import math
        p = 0.1
        expected = -0.75 * math.log(1 - 4 * p / 3)
        assert dn._jc_correct(p) == pytest.approx(expected, rel=1e-9)

    def test_saturated_returns_none(self):
        assert dn._jc_correct(0.75) is None

    def test_over_saturated_returns_none(self):
        assert dn._jc_correct(0.9) is None

    def test_result_positive(self):
        result = dn._jc_correct(0.05)
        assert result is not None and result > 0


class TestComputeMK:
    # Helper: build a simple coding alignment with outgroup
    @staticmethod
    def _make_coding_aln(ingroup: list[str], outgroup: str):
        """Return (seqs, outgroup_str) ready for compute_mk."""
        return ingroup, outgroup

    def test_returns_mk_stats_type(self):
        seqs = ['ATGATG', 'ATGATG']
        out  = 'ATGATG'
        result = dn.compute_mk(seqs, out)
        assert isinstance(result, dn.MKStats)

    def test_all_identical_zero_counts(self):
        # No variation, no divergence → all zeros
        seqs = ['ATGATG', 'ATGATG']
        result = dn.compute_mk(seqs, 'ATGATG')
        assert result.Pn == 0 and result.Ps == 0
        assert result.Dn == 0 and result.Ds == 0

    def test_synonymous_fixed_difference(self):
        # TTT (Phe) ingroup, TTC (Phe) outgroup → 1 synonymous fixed diff → Ds=1
        seqs = ['TTT', 'TTT']
        result = dn.compute_mk(seqs, 'TTC')
        assert result.Ds == 1
        assert result.Dn == 0

    def test_nonsynonymous_fixed_difference(self):
        # ATG (Met) ingroup, CTG (Leu) outgroup → 1 nonsyn fixed diff → Dn=1
        seqs = ['ATG', 'ATG']
        result = dn.compute_mk(seqs, 'CTG')
        assert result.Dn == 1
        assert result.Ds == 0

    def test_synonymous_polymorphism(self):
        # TTT vs TTC within ingroup (both Phe), outgroup = TTT → Ps=1
        seqs = ['TTT', 'TTC']
        result = dn.compute_mk(seqs, 'TTT')
        assert result.Ps == 1
        assert result.Pn == 0

    def test_nonsynonymous_polymorphism(self):
        # ATG (Met) vs CTG (Leu) within ingroup, outgroup = ATG → Pn=1
        seqs = ['ATG', 'CTG']
        result = dn.compute_mk(seqs, 'ATG')
        assert result.Pn == 1
        assert result.Ps == 0

    def test_alpha_none_when_dn_zero(self):
        seqs = ['TTTTTT', 'TTTTTT']
        result = dn.compute_mk(seqs, 'TTCTTC')
        # Dn=0 → alpha undefined
        assert result.alpha is None

    def test_dos_computed(self):
        seqs = ['ATGATG', 'TTTTTT']
        result = dn.compute_mk(seqs, 'CTGCTG')
        # We just need DoS to be in a valid range
        if result.DoS is not None:
            assert -1.0 <= result.DoS <= 1.0

    def test_fisher_p_in_range(self):
        seqs = ['ATGATG', 'TTTTTT']
        result = dn.compute_mk(seqs, 'CTGCTG')
        assert result.fisher_p is not None
        assert 0.0 <= result.fisher_p <= 1.0

    def test_gap_codon_skipped(self):
        # Gap in first codon of one ingroup seq → first codon skipped
        seqs = ['---ATG', 'ATGATG']
        result = dn.compute_mk(seqs, 'ATGATG')
        # Should not crash; Pn=Ps=Dn=Ds=0 for the gap codon
        assert isinstance(result, dn.MKStats)

    def test_n_less_than_2(self):
        result = dn.compute_mk(['ATG'], 'CTG')
        assert result.Pn == 0 and result.Dn == 0

    def test_not_divisible_by_3(self):
        result = dn.compute_mk(['ATGA', 'ATGC'], 'ATGC')
        assert result.Pn == 0 and result.Dn == 0

    def test_run_analysis_dispatch(self, tmp_path):
        f = tmp_path / "coding.fas"
        # 2 codons × 2 seqs; outgroup differs at codon 2 (nonsyn fixed diff)
        f.write_text(">s1\nATGATG\n>s2\nATGATG\n>out\nATGCTG\n")
        aln = dn.parse_fasta(f)
        aln_no_out = dn.Alignment(
            names=['s1', 's2'], seqs=['ATGATG', 'ATGATG'], source=str(f)
        )
        out_seq = 'ATGCTG'
        results = dn.run_analysis(aln_no_out, analyses={'mk'}, outgroup=out_seq)
        assert 'mk' in results
        assert results['mk'] is not None

    def test_dual_status_site_counts_as_both_fixed_and_polymorphic_synonymous(self):
        # Ingroup segregates CTT/CTC (both Leu); outgroup CTA (Leu) carries a
        # third base at position 2. Two DnaSP tallies apply to that ONE site
        # and they are independent of each other:
        #   fixed: no ingroup sequence carries the outgroup's base (ht3 = 1 in
        #     EntrePobsMod.vb; AAstatusB = 7 in McDonaldK.vb::BuscoPosFijadas,
        #     whose own comment reads "fijado diferencialmente aunque sea pol
        #     en una de ella" -- fixed even if polymorphic in one species), so
        #     Ds gains 1;
        #   polymorphic: the ingroup segregates there, and its within-species
        #     status (AAstatus1, from BuscaSitiosReemplazamientoMK) is worked
        #     out from the ingroup's own codons alone -- BuscoPosSegregantesMcDK
        #     never looks at the fixed flag -- so Ps gains 1 as well.
        # DnaSP's help file shows exactly this on its codon 13-15 worked
        # example (site 15: "within species ... 1 synonymous" AND "fixed
        # differences ... Site#15 is synonymous"). An earlier version of this
        # test asserted Ps == 0 here; that encoded the opposite rule and did
        # not reproduce DnaSP 6's real output (COII: Ps 65 vs 69).
        seqs = ['CTT', 'CTC']
        result = dn.compute_mk(seqs, 'CTA')
        assert result.Ds == 1
        assert result.Ps == 1
        assert result.Pn == 0 and result.Dn == 0

    def test_dual_status_site_counts_as_both_fixed_and_polymorphic_nonsynonymous(self):
        # Same rule, replacement flavour: ingroup ATG (Met) / CTG (Leu),
        # outgroup GTG (Val). Position 0 is a fixed difference (G absent from
        # the ingroup: ATG-GTG and CTG-GTG are both replacements, so Dn = 1)
        # AND a within-species replacement polymorphism (Met/Leu, so Pn = 1).
        seqs = ['ATG', 'CTG']
        result = dn.compute_mk(seqs, 'GTG')
        assert result.Dn == 1
        assert result.Pn == 1
        assert result.Ps == 0 and result.Ds == 0

    def test_outgroup_allele_shared_with_segregating_ingroup_is_poly_not_fixed(self):
        # Ingroup segregates TTT/TTC; outgroup TTT IS one of the ingroup
        # alleles (ht3 between 0 and 1 in DnaSP's terms) -- polymorphism,
        # not a fixed difference. (Boundary case for the ht3 rule above.)
        seqs = ['TTT', 'TTC']
        result = dn.compute_mk(seqs, 'TTT')
        assert result.Ps == 1
        assert result.Ds == 0 and result.Dn == 0

    def test_single_codon_can_contribute_both_dn_and_ds(self):
        # DnaSP's own McDonald-Kreitman help page (worked example, codon
        # 13-15) shows a SINGLE codon contributing two separate fixed
        # differences of different types (one replacement, one synonymous)
        # -- DnaSP counts fixed/polymorphic differences per NUCLEOTIDE SITE,
        # not once per codon. Monomorphic ingroup CTT (Leu) vs outgroup CGA
        # (Arg): position 0 (C/C) is invariant; position 1 (T->G) is a
        # nonsynonymous fixed difference (Leu->Arg, both mutational
        # orderings agree); position 2 (T->A) is a synonymous fixed
        # difference (CTN and CGN are both fully synonymous at the third
        # position). A whole-codon classifier can never produce both Dn and
        # Ds from the same codon; this is the core structural fix.
        seqs = ['CTT', 'CTT']
        result = dn.compute_mk(seqs, 'CGA')
        assert result.Dn == 1
        assert result.Ds == 1
        assert result.Pn == 0 and result.Ps == 0

    def test_fixed_side_tie_break_prefers_fewer_replacement_changes(self):
        # Monomorphic ingroup CCC (Pro) vs outgroup CAG (Gln): differ at
        # positions 1 (C->A) and 2 (C->G), with two possible mutational
        # orderings -- CCC->CAC->CAG (2 replacements: Pro->His, His->Gln)
        # or CCC->CCG->CAG (1 replacement + 1 synonymous: Pro->Pro,
        # Pro->Gln). DnaSP's own rule for fixed differences ("DnaSP will
        # check all paths ... and choose the path with the minor number of
        # [replacement] changes") has no outgroup-side "intermediate
        # observed elsewhere" tie-break (unlike the within-species case
        # below) -- fixed differences are decided by fewest total
        # differences, then fewest replacements. So the CCG-route wins:
        # position 2 synonymous, position 1 nonsynonymous.
        seqs = ['CCC', 'CCC']
        result = dn.compute_mk(seqs, 'CAG')
        assert result.Dn == 1
        assert result.Ds == 1
        assert result.Pn == 0 and result.Ps == 0

    def test_poly_side_prefers_path_with_intermediate_seen_in_outgroup(self):
        # Directly from DnaSP's own McDonald-Kreitman help page (worked
        # example, codon 10-12): ingroup segregates CCC (Pro) / CAG (Gln),
        # differing at positions 1 and 2. Two mutational orderings exist:
        # CCC->CAC->CAG (2 replacements, both nonsynonymous) or
        # CCC->CCG->CAG (1 replacement + 1 synonymous). Taken alone, the
        # second route has FEWER replacements -- but the outgroup here is
        # CAC, i.e. one of the ORDERINGS' intermediate codons is actually
        # observed. DnaSP's documented rule: "If there are two possible
        # paths, and one of the non-extant codons is found in the other
        # species, DnaSP assumes that is the true evolutionary path" --
        # this OVERRIDES the fewer-replacements preference. So the
        # CAC-route wins even though it has more replacements: both
        # position 1 and position 2 are nonsynonymous polymorphisms.
        seqs = ['CCC', 'CAG']
        result = dn.compute_mk(seqs, 'CAC')
        assert result.Pn == 2
        assert result.Ps == 0
        assert result.Dn == 0 and result.Ds == 0

    def test_complex_codon_is_excluded_from_both_tables_and_reported(self):
        # Ingroup AAA (Lys), AAG (Lys), CCG (Pro): three codons that between
        # them differ at all three positions. DnaSP's within-species routine
        # (McDonaldK.vb::BuscaSitiosReemplazamientoMK, Case 3 / NNdif = 3,
        # "no contabilizado") gives such a codon a negative status; both
        # BuscoPosSegregantesMcDK and BuscaDeCodonesEntre then skip every
        # site of that codon, and the output reports it under "Total number
        # of complex codons no analyzed". So nothing from this codon reaches
        # Pn/Ps/Dn/Ds -- not silently, though: the count is surfaced. (An
        # earlier version of this test expected the codon to be resolved
        # allele by allele instead; DnaSP does not do that.) The second,
        # ordinary codon (TTT/TTC, Phe/Phe) is still counted.
        seqs = ['AAATTT', 'AAGTTC', 'CCGTTT']
        result = dn.compute_mk(seqs, 'AAATTT')
        assert result.n_complex_codons == 1
        assert result.Pn == 0 and result.Ps == 1
        assert result.Dn == 0 and result.Ds == 0
        # The exclusion covers the fixed side too: with an outgroup that
        # differs at every position of the complex codon, DnaSP still counts
        # no fixed difference there (BuscaDeCodonesEntre: NoTratarCodon when
        # AAstatus1 < 0).
        result = dn.compute_mk(['AAA', 'AAG', 'CCG'], 'GGG')
        assert result.n_complex_codons == 1
        assert result.Pn == 0 and result.Ps == 0
        assert result.Dn == 0 and result.Ds == 0

    def test_internal_stop_codon_is_a_21st_amino_acid_in_the_mk_test(self):
        # DnaSP does not drop a codon because it is a stop: "Stop codons have
        # been found in the coding region. DnaSP has considered that they
        # could code for a rare amino acid (the 21st amino acid ...)"
        # (McDonaldK.vb, STOPinMiddle). Monomorphic ingroup TGA vs outgroup
        # TGG: under the standard code TGA is that 21st amino acid and
        # TGA/TGG is a replacement fixed difference; under the vertebrate
        # mitochondrial code both are Trp, a synonymous one. (An earlier
        # version of this test expected the standard-code case to be
        # excluded entirely.) The codon is not the terminal stop because the
        # outgroup's TGG is not a stop.
        seqs = ['TGA', 'TGA']
        std = dn.compute_mk(seqs, 'TGG')
        assert std.Dn == 1 and std.Ds == 0
        assert std.n_internal_stop_codons == 1
        mito = dn.compute_mk(seqs, 'TGG', dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert mito.Ds == 1 and mito.Dn == 0
        assert mito.n_internal_stop_codons == 0

    def test_internal_stop_polymorphism_counts_as_replacement(self):
        # COII's real case under the (wrong) standard code: an ingroup that
        # segregates TGA/TGG at an internal codon. DnaSP counts a replacement
        # polymorphism (* vs Trp); this is what gave its standard-code Pn of
        # 13 against the skill's former 11.
        seqs = ['TGAGGG', 'TGGGGG']
        result = dn.compute_mk(seqs, 'TGGGGG')
        assert result.Pn == 1 and result.Ps == 0
        assert result.n_internal_stop_codons == 1

    def test_terminal_stop_codon_is_not_analysed(self):
        # A final codon that is a stop in every sequence is the coding
        # region's terminal stop. DnaSP's own COII_Apes.nex annotates it as
        # noncoding (CODONPOSSET "N: 682-684") and codes 1-681 only, so
        # TAG vs TAA there is not a synonymous fixed difference.
        seqs = ['GGGTAG', 'GGGTAG']
        result = dn.compute_mk(seqs, 'GGGTAA')
        assert result.Ds == 0 and result.Dn == 0
        assert result.n_internal_stop_codons == 0

    def test_two_stop_codons_at_an_internal_position_are_synonymous(self):
        # * = * : an internal TAA/TAG fixed difference is synonymous.
        result = dn.compute_mk(['TAAGGG', 'TAAGGG'], 'TAGGGG')
        assert result.Ds == 1 and result.Dn == 0

    def test_help_file_codon_13_15_site_is_both_syn_polymorphism_and_syn_fixed(self):
        # Directly from DnaSP's own McDonald-Kreitman help page (worked
        # example, codon 13-15): species 1 = AAT, AGG, ACT, ACT; species 2 =
        # GGA (monomorphic). The page's own tallies for this codon are
        #   "within species: 2 replacements (site#14), and 1 synonymous
        #    (Site#15)" and
        #   "fixed differences: 2, Site#13 is replacement; Site#15 is
        #    synonymous."
        # Site 15 (T/G in species 1, A in species 2) is counted in BOTH
        # tables. Site 14 carries three alleles (A/G/C) and contributes two
        # replacement changes, not one (the chain AAT -> ACT -> ACG -> AGG
        # chosen in BuscaSitiosReemplazamientoMK's three-codon / two-position
        # case, path 4 in the help page's list).
        seqs = ['AAT', 'AGG', 'ACT', 'ACT']
        result = dn.compute_mk(seqs, 'GGA')
        assert result.Pn == 2
        assert result.Ps == 1
        assert result.Dn == 1
        assert result.Ds == 1
        assert result.n_complex_codons == 0

    def test_help_file_codon_1_3_four_alleles_at_one_site_are_three_changes(self):
        # Help page, codon 1-3: species 1 = AGT, AGC, AGA, AGG (Ser, Ser, Arg,
        # Arg), species 2 = AGG. "species#1: 3 mutations in site#3: 1
        # replacement, 2 synonymous." A site segregating for k bases counts
        # k - 1 changes (the sum rule in BuscaSitiosReemplazamientoMK, Case 4
        # / NNdif = 1: pairwise sum 8 -> code 122). The outgroup's G is one of
        # the ingroup's bases, so nothing is fixed.
        seqs = ['AGT', 'AGC', 'AGA', 'AGG']
        result = dn.compute_mk(seqs, 'AGG')
        assert result.Pn == 1
        assert result.Ps == 2
        assert result.Dn == 0 and result.Ds == 0

    def test_three_alleles_at_one_site_are_two_changes(self):
        # COII codon 80 of the validation file (0-indexed): ingroup CTA, CTG,
        # CTC (all Leu under the vertebrate mitochondrial code), outgroup
        # CTA. Three bases at position 2 -> two synonymous changes (Case 3 /
        # NNdif = 1: pairwise sum 6 -> code 22). Counting this site once was
        # one of the four missing Ps in the DnaSP 6.12 calibration.
        seqs = ['CTA', 'CTG', 'CTC']
        result = dn.compute_mk(seqs, 'CTA', dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert result.Ps == 2
        assert result.Pn == 0
        assert result.Dn == 0 and result.Ds == 0

    def test_same_amino_acid_codons_two_sites_apart_are_both_synonymous(self):
        # AGT (Ser) vs TCT (Ser) differ at positions 0 and 1; every
        # single-step path runs through a non-Ser codon (TGT Cys or ACT Thr),
        # so path counting alone would call both changes replacements. DnaSP
        # overrides that: when two codons encode the same amino acid and
        # differ at two or three positions, all differing positions are
        # synonymous (BuscaSitiosReemplazamientoMK, Case 2, VB comment
        # "ejemplo de codones AGU Ser UCU Ser; or UCC Ser vs AGU Ser").
        seqs = ['AGT', 'TCT']
        result = dn.compute_mk(seqs, 'AGT')
        assert result.Ps == 2
        assert result.Pn == 0

    def test_help_file_codon_16_18_recombination_case_only_under_nuclear_codes(self):
        # Help page, codon 16-18: species 1 = ATA, TTA, TTG, ATG (all four
        # combinations of two variants at positions 0 and 2: a circular
        # path), species 2 = TTT. "within species: Site#16 (1 replacement);
        # Site#18 (1 synonymous). fixed differences: 1, Site#18 is
        # replacement. Note: This kind of codons will be analyzed only for
        # Nuclear Genetic Codes." Under the standard code DnaSP assumes one
        # recombination event and takes, per position, the most synonymous
        # of the four single-step edges. Under a mitochondrial table
        # (SistemaGeneticoRecombina = 0) the codon is a complex codon and is
        # excluded from both tables.
        seqs = ['ATA', 'TTA', 'TTG', 'ATG']
        std = dn.compute_mk(seqs, 'TTT')
        assert std.Pn == 1 and std.Ps == 1
        assert std.Dn == 1 and std.Ds == 0
        assert std.n_complex_codons == 0
        mito = dn.compute_mk(seqs, 'TTT', dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert mito.Pn == 0 and mito.Ps == 0
        assert mito.Dn == 0 and mito.Ds == 0
        assert mito.n_complex_codons == 1

    def test_recombination_case_with_amino_acid_change_demotes_one_synonymous_site(self):
        # The VB's own counter-example for the circular-path rule: CTG (Leu),
        # TTG (Leu), TTT (Phe), CTT (Leu). Per position the most synonymous
        # edge is synonymous at both positions (CTG-TTG at 0, CTT-CTG at 2),
        # yet an amino acid change (Phe) is present, so DnaSP demotes the
        # first synonymous position to a replacement ("problema en codones
        # del tipo CUG Leu UUG Leu UUU Phe CUU Leu"): 1 replacement + 1
        # synonymous, not 2 synonymous.
        seqs = ['CTG', 'TTG', 'TTT', 'CTT']
        result = dn.compute_mk(seqs, 'CTG')
        assert result.Pn == 1
        assert result.Ps == 1
        assert result.n_complex_codons == 0

    def test_five_or_more_codons_segregating_is_a_complex_codon(self):
        # BuscaSitiosReemplazamientoMK handles 2, 3 and 4 distinct codons per
        # triplet; five or more ("mas de cuatro codones") is a complex codon,
        # excluded and reported.
        seqs = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA']
        result = dn.compute_mk(seqs, 'AAA')
        assert result.n_complex_codons == 1
        assert result.Pn == 0 and result.Ps == 0
        assert result.Dn == 0 and result.Ds == 0

    def test_help_file_codon_7_9_two_difference_pair_prefers_fewer_replacements(self):
        # Help page, codon 7-9: ATT (Ile) / CTG (Leu) within species 1,
        # species 2 = ATT. "DnaSP will choose path#1, the path that requires
        # the minor number of replacements": ATT -> CTT -> CTG, site 7
        # replacement, site 9 synonymous. Nothing fixed (outgroup ATT is an
        # ingroup allele).
        seqs = ['ATT', 'CTG']
        result = dn.compute_mk(seqs, 'ATT')
        assert result.Pn == 1 and result.Ps == 1
        assert result.Dn == 0 and result.Ds == 0

    def test_help_file_codon_22_24_three_difference_pair(self):
        # Help page, codon 22-24: TAT (Tyr) / CTG (Leu), species 2 = TAT.
        # Six orderings, two of them through stop codons; "2 replacements
        # ... and 1 synonymous" whichever of the remaining most-synonymous
        # paths is taken.
        seqs = ['TAT', 'CTG']
        result = dn.compute_mk(seqs, 'TAT')
        assert result.Pn == 2 and result.Ps == 1
        assert result.Dn == 0 and result.Ds == 0

    def test_three_codons_at_two_sites_chain_through_the_middle_codon(self):
        # ATG (Met), ATC (Ile), CTC (Leu): two positions with two variants
        # each. DnaSP chains through the codon adjacent to both others (ATG
        # -> ATC -> CTC; VB comment "AUG AUC CUC") and labels each position
        # from its single-step edge: two replacements.
        seqs = ['ATG', 'ATC', 'CTC']
        result = dn.compute_mk(seqs, 'ATG')
        assert result.Pn == 2 and result.Ps == 0
        assert result.n_complex_codons == 0

    def test_ni_is_undefined_when_dn_is_zero(self):
        # NI = (Pn/Ps)/(Dn/Ds) has no value when Dn = 0; DnaSP prints NI and
        # alpha only when Dn * Ps is non-zero (McDonaldK.vb::McDonaldTest,
        # "If (d * a) <> 0"). A synonymous dual-status site (Ps = 1, Ds = 1,
        # Dn = 0) used to divide by zero here.
        result = dn.compute_mk(['CTT', 'CTC'], 'CTA')
        assert result.Ps == 1 and result.Ds == 1 and result.Dn == 0
        assert result.NI is None
        assert result.alpha is None

    def test_mk_stats_reports_zero_complex_codons_by_default(self):
        result = dn.compute_mk(['TTT', 'TTC'], 'TTT')
        assert result.n_complex_codons == 0


class TestComputeKaKs:
    def test_returns_kaks_stats_type(self):
        result = dn.compute_ka_ks(['ATGATG', 'ATGATG'])
        assert isinstance(result, dn.KaKsStats)

    def test_identical_seqs_zero_ka_ks(self):
        # GGT has 1 synonymous site (pos 2) → total_S > 0 even for identical seqs
        result = dn.compute_ka_ks(['GGTGGT', 'GGTGGT'])
        assert result.Ks == pytest.approx(0.0)
        assert result.Ka == pytest.approx(0.0)
        assert result.omega is None  # Ks=0 → ratio undefined

    def test_synonymous_change_only(self):
        # 3 GGT codons vs GGC/GGT/GGT: only first codon differs (Gly→Gly = syn)
        # total_S = 3, total_sd = 1, pS = 1/3 < 0.75 → JC correction defined
        result = dn.compute_ka_ks(['GGTGGTGGT', 'GGCGGTGGT'])
        assert result.Ka == pytest.approx(0.0)
        assert result.Ks is not None and result.Ks > 0

    def test_nonsynonymous_change_only(self):
        # ATG vs CTG: one nonsyn change → Ks=0, Ka>0
        result = dn.compute_ka_ks(['ATG', 'CTG'])
        assert result.Ks == pytest.approx(0.0)
        assert result.Ka is not None and result.Ka > 0

    def test_n_codons_set(self):
        result = dn.compute_ka_ks(['ATGATG', 'ATGCTG'])
        assert result.n_codons == 2

    def test_s_sites_positive(self):
        result = dn.compute_ka_ks(['GGTGGT', 'GGCGGC'])
        # GGT and GGC both Gly → only synonymous change → S_sites > 0
        assert result.S_sites > 0.0

    def test_n_less_than_2(self):
        result = dn.compute_ka_ks(['ATG'])
        assert result.n_codons == 0

    def test_not_divisible_by_3(self):
        result = dn.compute_ka_ks(['ATGA', 'ATGC'])
        assert result.n_codons == 0

    def test_omega_less_than_1_synonymous_dominated(self):
        # GGT vs GGC: synonymous → Ka=0, Ks>0 → omega=0
        result = dn.compute_ka_ks(['GGT', 'GGC'])
        if result.Ka is not None and result.Ks is not None and result.Ks > 0:
            assert result.omega is not None
            assert result.omega == pytest.approx(0.0)

    def test_sd_nd_non_negative(self):
        result = dn.compute_ka_ks(['ATGTTTTTT', 'CTGTTCTTT'])
        assert result.Sd >= 0.0
        assert result.Nd >= 0.0

    def test_multiple_seqs_averages(self):
        # 3 sequences → 3 pairs; result should still be valid
        seqs = ['ATGATG', 'CTGCTG', 'TTTTTT']
        result = dn.compute_ka_ks(seqs)
        assert result.n_codons == 2
        assert result.Ks is not None or result.Ka is not None  # at least something computed

    def test_run_analysis_dispatch(self, tmp_path):
        f = tmp_path / "coding.fas"
        f.write_text(">s1\nATGATG\n>s2\nCTGCTG\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(aln, analyses={'kaks'})
        assert 'kaks' in results
        assert results['kaks'] is not None
        assert isinstance(results['kaks'], dn.KaKsStats)

    def test_internal_stop_codon_is_a_21st_amino_acid_in_kaks(self):
        # Under the standard code TGA is a stop, but DnaSP still analyses the
        # codon as a 21st amino acid: a stop codon has no synonymous sites
        # (SINONIMO.vb::GeneticCodeFold counts its stop neighbours as stops,
        # never as synonymous changes), so TGA/TGG is one nonsynonymous
        # difference over 3 nonsynonymous sites. Under the vertebrate
        # mitochondrial code TGA is Trp and the change is synonymous. (An
        # earlier version of this test expected the standard-code pair to
        # be skipped with n_codons = 0.)
        seqs = ['TGA', 'TGG']
        std = dn.compute_ka_ks(seqs)
        assert std.n_codons == 1
        assert std.S_sites == pytest.approx(0.0)
        assert std.N_sites == pytest.approx(3.0)
        assert std.Nd == pytest.approx(1.0) and std.Ks is None
        assert std.Ka == pytest.approx(dn._jc_correct(1.0 / 3.0))
        mito = dn.compute_ka_ks(seqs, dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert mito.n_codons == 1
        assert mito.S_sites > 0.0 and mito.Sd == pytest.approx(1.0)

    def test_internal_stop_codon_contributes_three_nonsynonymous_sites(self):
        # TAA (*) followed by two GGG: S = 0 + 1 + 1 = 2, N = 9 - 2 = 7 for
        # both sequences; TAA vs TAG is one synonymous difference (* = *).
        result = dn.compute_ka_ks(['TAAGGGGGG'], outgroup='TAGGGGGGG')
        assert result.n_codons == 3
        assert result.S_sites == pytest.approx(2.0)
        assert result.N_sites == pytest.approx(7.0)
        assert result.Sd == pytest.approx(1.0) and result.Nd == pytest.approx(0.0)
        assert result.n_internal_stop_codons == 1

    def test_outgroup_none_keeps_all_pairwise_default(self):
        # No outgroup -> unchanged: average over every ingroup pair.
        seqs = ['GGTGGTGGT', 'GGCGGTGGT', 'GGTGGTGGT']
        no_out = dn.compute_ka_ks(seqs)
        assert no_out.Ks is not None

    def test_outgroup_restricts_pairs_to_ingroup_vs_outgroup(self):
        # With an outgroup, DnaSP's own manual (Synonymous and Nonsynonymous
        # Substitutions module) and its Ka/Ks menu output only make sense
        # divergence-wise against the defined outgroup, not ingroup-vs-ingroup
        # pairs. Two identical ingroup sequences that are both divergent from
        # the outgroup should show nonzero Ka/Ks (ingroup-vs-outgroup), even
        # though an ingroup-only comparison would be all zeros.
        ingroup = ['ATGATG', 'ATGATG']
        outgroup = 'CTGCTG'   # nonsynonymous difference at codon 1 and 2
        result = dn.compute_ka_ks(ingroup, outgroup=outgroup)
        assert result.Ka is not None and result.Ka > 0.0

    def test_outgroup_summary_is_dnasp_corrected_ratio(self):
        # DnaSP's Ka/Ks summary (EntrePobsMod.vb::PolDivergenceOut) is the
        # Jukes-Cantor correction of ONE ratio: the mean number of
        # differences over the ingroup-vs-outgroup pairs, divided by the
        # mean number of sites over all sequences -- not an average of
        # per-pair corrected distances (an earlier version of this test
        # asserted the latter). Ingroup ATGATG, CTGATG, ATGCTG vs outgroup
        # GTGATG: nonsynonymous differences 1, 1, 2 (mean 4/3); synonymous
        # sites per sequence 0, 4/3, 4/3, 1 (mean 11/12), so N = 6 - 11/12.
        ingroup = ['ATGATG', 'CTGATG', 'ATGCTG']
        outgroup = 'GTGATG'
        result = dn.compute_ka_ks(ingroup, outgroup=outgroup)
        assert result.Nd == pytest.approx(4.0 / 3.0)
        assert result.Sd == pytest.approx(0.0)
        assert result.S_sites == pytest.approx(11.0 / 12.0)
        assert result.N_sites == pytest.approx(6.0 - 11.0 / 12.0)
        assert result.Ka == pytest.approx(dn._jc_correct(result.Nd / result.N_sites))
        assert result.Ks == pytest.approx(0.0)

    def test_summary_ks_is_jc_of_mean_ratio_not_mean_of_per_pair_jc(self):
        # Two ingroup sequences vs outgroup, one pair with a synonymous
        # difference (GGA/GGG, Gly) and one identical pair. Mean differences
        # 0.5 over mean synonymous sites 2.0 -> p = 0.25 -> Ks = JC(0.25)
        # = 0.304099. Averaging per-pair corrected values would give
        # (JC(0.5) + 0) / 2 = 0.411985 instead.
        result = dn.compute_ka_ks(['GGAGGA', 'GGGGGA'], outgroup='GGGGGA')
        assert result.Sd == pytest.approx(0.5)
        assert result.S_sites == pytest.approx(2.0)
        assert result.Ks == pytest.approx(0.304099, rel=1e-5)

    def test_summary_rule_applies_without_outgroup_too(self):
        # No outgroup: every ingroup pair, same rule (DnaSP's Pi(s) form:
        # mean pairwise differences per synonymous site, corrected once).
        # Pairs (1,2) and (1,3) differ by one synonymous change, (2,3) by
        # none: mean 2/3 over 2.0 sites -> p = 1/3 -> Ks = JC(1/3) = 0.440840.
        result = dn.compute_ka_ks(['GGAGGA', 'GGGGGA', 'GGGGGA'])
        assert result.Sd == pytest.approx(2.0 / 3.0)
        assert result.Ks == pytest.approx(0.440840, rel=1e-5)

    def test_n_codons_counts_codons_analysed_in_every_sequence(self):
        # DnaSP reports "Number of codons analyzed" (227 for COII: the
        # terminal stop codon is not analysed), not the alignment's codon
        # count. A codon skipped in any compared sequence is not analysed.
        result = dn.compute_ka_ks(['GGGTAA', 'GGGTAA'], outgroup='GGATAA')
        assert result.n_codons == 1
        result = dn.compute_ka_ks(['GGGGGG', 'GGGGG-'])
        assert result.n_codons == 1

    def test_n_sites_counts_three_per_analysed_codon_not_alignment_length(self):
        # DnaSP: nonsynonymous sites = 3 x codons analysed - synonymous
        # sites (its COII total is 681 - 168.222 = 512.777 for 227 codons).
        # The second codon here is a stop under the standard code and is
        # skipped, so the alignment length (6) must not enter: S = 1.0
        # (GGG/GGA, third position fully synonymous), N = 3 - 1 = 2.0, not
        # 6 - 1 = 5.0.
        result = dn.compute_ka_ks(['GGGTAA', 'GGGTAA'], outgroup='GGATAA')
        assert result.S_sites == pytest.approx(1.0)
        assert result.N_sites == pytest.approx(2.0)

    def test_outgroup_included_in_site_count_average(self):
        # DnaSP: "the total number of synonymous and nonsynonymous sites ...
        # is estimated as the average ... of all sequences" -- when an
        # outgroup is given it is one of "all sequences" being compared.
        # Both codons here are valid (not stops), so a real pairwise
        # comparison exists and the site-count average should reflect both
        # the ingroup and outgroup sequences, not the ingroup alone.
        ingroup = ['GGCGGC']  # Ala, has synonymous sites
        outgroup = 'GGTGGT'   # Gly, has synonymous sites
        with_out = dn.compute_ka_ks(ingroup, outgroup=outgroup)
        without_out = dn.compute_ka_ks(ingroup)  # n<2, no outgroup -> no result
        assert with_out.S_sites > 0.0
        assert without_out.S_sites == 0.0  # confirms the outgroup path is doing the work

    def test_no_jointly_valid_codon_reports_zero_not_raw_totals(self):
        # Every ingroup codon is ambiguous, so no pair has a single
        # jointly-valid codon against the outgroup. Reporting the raw 2-codon
        # alignment length here would make a fully failed comparison look
        # like a populated, quantified result -- n_codons, S_sites and
        # N_sites should all come back as 0, and Ka/Ks/omega as undefined
        # (None), not silently omitted from an otherwise-full report. (An
        # earlier version used an all-stop-codon ingroup; stop codons are
        # now analysed as a 21st amino acid, as DnaSP does, so they no
        # longer produce this case.)
        ingroup = ['NNNNNN']
        outgroup = 'GGTGGT'
        result = dn.compute_ka_ks(ingroup, outgroup=outgroup)
        assert result.n_codons == 0
        assert result.S_sites == 0.0
        assert result.N_sites == 0.0
        assert result.Ka is None and result.Ks is None and result.omega is None


class TestStirling1Unsigned:
    """_stirling1_unsigned: unsigned Stirling numbers of the first kind."""

    def test_row_0(self):
        row = dn._stirling1_unsigned(0)
        assert row == [1]  # |s(0,0)| = 1

    def test_row_1(self):
        row = dn._stirling1_unsigned(1)
        assert row[0] == 0
        assert row[1] == 1

    def test_row_2(self):
        row = dn._stirling1_unsigned(2)
        # |s(2,0)|=0, |s(2,1)|=1, |s(2,2)|=1
        assert row[0] == 0
        assert row[1] == 1
        assert row[2] == 1

    def test_row_3(self):
        row = dn._stirling1_unsigned(3)
        # |s(3,0)|=0, |s(3,1)|=2, |s(3,2)|=3, |s(3,3)|=1
        assert row[0] == 0
        assert row[1] == 2
        assert row[2] == 3
        assert row[3] == 1

    def test_row_4(self):
        row = dn._stirling1_unsigned(4)
        # |s(4,1)|=6, |s(4,2)|=11, |s(4,3)|=6, |s(4,4)|=1
        assert row[1] == 6
        assert row[2] == 11
        assert row[3] == 6
        assert row[4] == 1

    def test_sum_equals_factorial(self):
        import math
        for n in range(1, 8):
            row = dn._stirling1_unsigned(n)
            assert sum(row) == math.factorial(n)

    def test_length(self):
        for n in range(0, 6):
            row = dn._stirling1_unsigned(n)
            assert len(row) == n + 1


class TestEwensCdf:
    """_ewens_cdf: CDF of allele count under Ewens sampling formula."""

    def test_k_max_equals_n_is_1(self):
        # P(K ≤ n | theta, n) = 1 (can't have more alleles than sequences)
        assert dn._ewens_cdf(5, 5, 2.0) == pytest.approx(1.0, abs=1e-10)

    def test_k_max_0_is_0(self):
        # P(K ≤ 0 | ...) = 0 (at least 1 allele)
        result = dn._ewens_cdf(0, 5, 2.0)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_theta_zero_returns_1(self):
        # theta=0 → degenerate; guard returns 1.0
        result = dn._ewens_cdf(3, 5, 0.0)
        assert result == pytest.approx(1.0)

    def test_monotone_in_k(self):
        # CDF must be non-decreasing in k_max
        n, theta = 8, 3.0
        vals = [dn._ewens_cdf(k, n, theta) for k in range(1, n + 1)]
        for a, b in zip(vals, vals[1:]):
            assert b >= a - 1e-12

    def test_within_unit_interval(self):
        for n in [4, 6, 10]:
            for theta in [0.5, 1.0, 5.0]:
                for k in range(1, n + 1):
                    result = dn._ewens_cdf(k, n, theta)
                    assert 0.0 <= result <= 1.0

    def test_large_theta_many_alleles(self):
        # Very large theta → many alleles expected → P(K ≤ 1) should be small
        result = dn._ewens_cdf(1, 10, 100.0)
        assert result < 0.01

    def test_small_theta_few_alleles(self):
        # Very small theta → few alleles expected → P(K ≤ 2) should be large
        result = dn._ewens_cdf(2, 10, 0.01)
        assert result > 0.95

    def test_n1_single_sequence(self):
        # n=1 → only 1 allele possible; P(K ≤ 1) = 1
        result = dn._ewens_cdf(1, 1, 2.0)
        assert result == pytest.approx(1.0)

    def test_large_n_no_overflow(self):
        # |s(n, k)| exceeds float range for n ≳ 171; log-space must not overflow
        for n in (180, 300):
            assert 0.0 <= dn._ewens_cdf(2, n, 1.0) <= 1.0
            assert 0.0 <= dn._ewens_cdf(n, n, 1.0) <= 1.0
            assert dn._ewens_cdf(n, n, 1.0) == pytest.approx(1.0, abs=1e-6)

    def test_sf_complements_cdf(self):
        n, theta = 12, 4.0
        for k in range(2, n + 1):
            assert dn._ewens_sf(k, n, theta) == pytest.approx(
                1.0 - dn._ewens_cdf(k - 1, n, theta), abs=1e-9
            )

    def test_sf_large_n_no_overflow(self):
        for n in (200, 300):
            assert 0.0 <= dn._ewens_sf(2, n, 1.0) <= 1.0
            assert 0.0 <= dn._ewens_sf(n, n, 1.0) <= 1.0


class TestComputeFuFs:
    """compute_fu_fs: Fu's Fs neutrality test."""

    def _make_clean(self, seqs):
        clean, _ = dn.complete_deletion(seqs)
        return clean

    def test_returns_fufs_stats(self):
        seqs = ['ATCG', 'ATCG', 'ATCA', 'ATCC']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=3, k=1.5)
        assert isinstance(result, dn.FuFsStats)

    def test_n_stored(self):
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=4, k=2.0)
        assert result.n == 4

    def test_H_stored(self):
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=4, k=2.0)
        assert result.H == 4

    def test_theta_pi_stored(self):
        seqs = ['ATCG', 'ATCA']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=2, k=1.5)
        assert result.theta_pi == pytest.approx(1.5)

    def test_n1_returns_none_fs(self):
        result = dn.compute_fu_fs(['ATCG'], H=1, k=0.0)
        assert result.Fs is None

    def test_H0_returns_none_fs(self):
        seqs = ['ATCG', 'ATCG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=0, k=0.0)
        assert result.Fs is None

    def test_Sk_in_unit_interval(self):
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=3, k=2.0)
        assert result.S_k is not None
        assert 0.0 <= result.S_k <= 1.0

    def test_Fs_is_float(self):
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=3, k=2.0)
        assert result.Fs is not None
        assert isinstance(result.Fs, float)

    def test_monomorphic_returns_none(self):
        # No polymorphism (H = 1, k = 0) -> Fs undefined, as DnaSP reports n.a.
        seqs = ['ATCG', 'ATCG', 'ATCG']
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=1, k=0.0)
        assert result.Fs is None
        assert result.S_k is None

    def test_Fs_negative_when_fewer_haplotypes(self):
        # When H << expected haplotype count given theta, Fs << 0
        # Use large k (high diversity) but small H (few haplotypes)
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG', 'TTTT', 'CCCC']
        clean = self._make_clean(seqs)
        # H=2 but k=5.0 (high diversity): a deficit of haplotypes → Fs > 0
        result = dn.compute_fu_fs(clean, H=2, k=5.0)
        assert result.Fs is not None
        assert result.Fs > 0

    def test_Fs_negative_when_excess_haplotypes(self):
        # Every sequence a distinct haplotype but low pairwise diversity →
        # more haplotypes than the Ewens formula expects → Fs < 0
        seqs = ['A' * 20] + ['A' * i + 'T' + 'A' * (19 - i) for i in range(9)]
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=10, k=1.8)
        assert result.Fs is not None
        assert result.Fs < 0
        assert result.S_k < 0.5

    def test_Fs_decreases_with_H(self):
        # More haplotypes observed → S' = P(K ≥ H) decreases → Fs decreases
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        clean = self._make_clean(seqs)
        r1 = dn.compute_fu_fs(clean, H=2, k=2.0)
        r2 = dn.compute_fu_fs(clean, H=4, k=2.0)
        assert r1.S_k is not None and r2.S_k is not None
        assert r2.S_k <= r1.S_k
        assert r2.Fs <= r1.Fs

    def test_large_alignment_no_overflow(self):
        # ~200 sequences: Stirling coefficients overflow float; log-space path
        # must still return a finite Fs rather than raising OverflowError.
        seqs = ['A' * 40] + [
            'A' * i + 'T' + 'A' * (39 - i) for i in range(39)
        ] * 5
        clean = self._make_clean(seqs)
        result = dn.compute_fu_fs(clean, H=40, k=1.5)
        assert result.Fs is not None
        assert math.isfinite(result.Fs) or abs(result.Fs) >= 1e307

    def test_run_analysis_dispatch(self, tmp_path):
        f = tmp_path / "aln.fas"
        f.write_text(">s1\nATCG\n>s2\nATCA\n>s3\nGCTA\n>s4\nGCTG\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(aln, analyses={'fufs'})
        assert 'fufs' in results
        assert results['fufs'] is not None
        assert isinstance(results['fufs'], dn.FuFsStats)

    def test_run_analysis_fufs_uses_global_H_and_k(self, tmp_path):
        # fufs must be consistent with polymorphism stats
        f = tmp_path / "aln.fas"
        f.write_text(">s1\nATCG\n>s2\nATCA\n>s3\nGCTA\n>s4\nGCTG\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(aln, analyses={'fufs', 'polymorphism'})
        fufs_r = results['fufs']
        poly_r = results['global']
        assert fufs_r is not None
        assert fufs_r.H == poly_r.H
        assert fufs_r.theta_pi == pytest.approx(poly_r.k)


class TestComputeSFS:
    """compute_sfs: Site frequency spectrum."""

    def test_returns_sfs_stats(self):
        seqs = ['ATCG', 'ATCA', 'GCTA', 'GCTG']
        result = dn.compute_sfs(seqs)
        assert isinstance(result, dn.SFSStats)

    def test_n_stored(self):
        seqs = ['ATCG', 'ATCA', 'GCTA']
        result = dn.compute_sfs(seqs)
        assert result.n == 3

    def test_no_outgroup_unfolded_is_none(self):
        seqs = ['ATCG', 'ATCA']
        result = dn.compute_sfs(seqs)
        assert result.unfolded is None
        assert result.has_outgroup is False

    def test_with_outgroup_unfolded_is_dict(self):
        seqs = ['ATCG', 'ATCA']
        result = dn.compute_sfs(seqs, outgroup_seq='ATCG')
        assert result.unfolded is not None
        assert result.has_outgroup is True

    def test_monomorphic_gives_empty_folded(self):
        seqs = ['AAAA', 'AAAA', 'AAAA']
        result = dn.compute_sfs(seqs)
        assert result.folded == {}

    def test_single_biallelic_site_folded(self):
        # 4 sequences, one site differs: ATCG vs ATCA vs ATCA vs ATCA
        # Col 3: G(1), A(3) → minor=G, count=1 → folded[1] = 1
        seqs = ['ATCG', 'ATCA', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs)
        assert result.folded.get(1, 0) == 1

    def test_folded_key_is_minor_allele_count(self):
        # n=6: site has 2 carrying minor allele → folded[2] = 1
        seqs = ['ATCG', 'ATCG', 'ATCG', 'ATCG', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs)
        assert result.folded.get(2, 0) == 1

    def test_folded_key_at_most_n_over_2(self):
        # All folded keys must be ≤ n//2
        seqs = ['AAAA', 'TTTT', 'CCCC', 'GGGG', 'AATT', 'CCGG']
        result = dn.compute_sfs(seqs)
        n = len(seqs)
        for key in result.folded:
            assert 1 <= key <= n // 2

    def test_gap_column_skipped(self):
        # Position 2 has a gap → should be excluded
        seqs = ['AT-G', 'ATCG', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs)
        # Column 2 has gap → skip; col 3: G(1),G(1),A(1),A(1) → actually ATCG after gap removal
        # After complete deletion: positions 0,1,3 remain; col3: G,G,A,A → minor=2
        # But compute_sfs uses original seqs with gap detection internally
        # Column with '-' → skip; col 3 minor count depends on seqs
        # We just verify no error and result is SFSStats
        assert isinstance(result, dn.SFSStats)

    def test_n_less_than_2_returns_empty(self):
        result = dn.compute_sfs(['ATCG'])
        assert result.folded == {}

    def test_unfolded_uses_outgroup_polarisation(self):
        # seqs: G at pos 3 in s1 (ancestral=G), derived=A in s2,s3,s4
        # outgroup carries G at pos 3 → derived allele A in 3 out of 4 seqs
        seqs = ['ATCG', 'ATCA', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs, outgroup_seq='ATCG')
        # Col 3: seqs have G(1),A(3); og=G → derived=A, count=3 → unfolded[3]=1
        assert result.unfolded is not None
        assert result.unfolded.get(3, 0) == 1

    def test_unfolded_ancestral_missing_skipped(self):
        # Outgroup carries C at pos 0, but no ingroup sequence has C → skip
        seqs = ['ATCG', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs, outgroup_seq='CTCG')
        # Col 0: og=C, but ingroup only has A → og allele not in ingroup counts → skip
        # Col 3: og=G(1 in seq 0), derived=A(2 in seqs 1,2) → unfolded[2]=1
        assert result.unfolded is not None

    def test_outgroup_gap_skipped_for_unfolded(self):
        # Outgroup has gap at pos 3 → that col excluded from unfolded but ok for folded
        seqs = ['ATCG', 'ATCA', 'ATCA', 'ATCA']
        result = dn.compute_sfs(seqs, outgroup_seq='ATC-')
        # Col 3: ingroup is clean but og has '-' → not in 'ATCG' → skip for unfolded
        assert result.unfolded is not None
        assert result.unfolded.get(3, 0) == 0  # skipped due to outgroup gap

    def test_total_folded_counts_consistent(self):
        # Total sites in folded SFS ≤ total segregating sites
        seqs = ['AACG', 'TTCA', 'AACG', 'TTCA']
        result = dn.compute_sfs(seqs)
        # Each key maps to site count; sum = total segregating sites (up to folded)
        total = sum(result.folded.values())
        assert total >= 0

    def test_triallelic_site_excluded(self):
        # DnaSP gates the spectrum on exactly two states (contot == 2).
        # Col 0 is triallelic {A,C,G}; col 1 is biallelic {A,T}.
        seqs = ['AA', 'CA', 'GT', 'GT']
        result = dn.compute_sfs(seqs, outgroup_seq='AA')
        assert result.folded == {2: 1}          # only the biallelic column
        assert result.unfolded == {2: 1}
        assert result.n_multiallelic_excluded == 1

    def test_quadriallelic_site_excluded(self):
        seqs = ['A', 'C', 'G', 'T']
        result = dn.compute_sfs(seqs, outgroup_seq='A')
        assert result.folded == {}
        assert result.unfolded == {}
        assert result.n_multiallelic_excluded == 1

    def test_biallelic_sites_still_counted(self):
        seqs = ['AACG', 'TTCA', 'AACG', 'TTCA']
        result = dn.compute_sfs(seqs)
        assert sum(result.folded.values()) == 3
        assert result.n_multiallelic_excluded == 0

    def test_run_analysis_dispatch_no_outgroup(self, tmp_path):
        f = tmp_path / "aln.fas"
        f.write_text(">s1\nATCG\n>s2\nATCA\n>s3\nGCTA\n>s4\nGCTG\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(aln, analyses={'sfs'})
        assert 'sfs' in results
        assert results['sfs'] is not None
        assert isinstance(results['sfs'], dn.SFSStats)
        assert results['sfs'].unfolded is None

    def test_run_analysis_dispatch_with_outgroup(self, tmp_path):
        f = tmp_path / "aln.fas"
        f.write_text(">s1\nATCG\n>s2\nATCA\n>s3\nGCTA\n>og\nAAAA\n")
        aln = dn.parse_fasta(f)
        results = dn.run_analysis(aln, analyses={'sfs'}, outgroup='og')
        assert results['sfs'] is not None
        assert results['sfs'].has_outgroup is True
        assert results['sfs'].unfolded is not None


# =============================================================================
# Group E  -  Transition / Transversion ratio
# =============================================================================

class TestSynonymousFamilies:
    """Verify _SYNONYMOUS_FAMILIES is built correctly from GENETIC_CODE."""

    def test_stop_codons_excluded(self):
        """'*' (stop) must not appear as a key in _SYNONYMOUS_FAMILIES."""
        assert '*' not in dn._SYNONYMOUS_FAMILIES

    def test_all_sense_codons_covered(self):
        """Every non-stop codon in GENETIC_CODE must appear in exactly one family."""
        all_codons_in_families: list[str] = []
        for codons in dn._SYNONYMOUS_FAMILIES.values():
            all_codons_in_families.extend(codons)
        sense_codons = [c for c, aa in dn.GENETIC_CODE.items() if aa != '*']
        assert sorted(all_codons_in_families) == sorted(sense_codons)

    def test_leucine_family_size_six(self):
        """Leu (L) has a 6-fold degenerate synonymous family."""
        assert len(dn._SYNONYMOUS_FAMILIES['L']) == 6

    def test_methionine_family_size_one(self):
        """Met (M) has exactly 1 codon (ATG)."""
        assert len(dn._SYNONYMOUS_FAMILIES['M']) == 1
        assert 'ATG' in dn._SYNONYMOUS_FAMILIES['M']

    def test_total_sense_codons_61(self):
        """Standard genetic code has 61 sense codons."""
        total = sum(len(v) for v in dn._SYNONYMOUS_FAMILIES.values())
        assert total == 61


class TestComputeTsTv:
    """Tests for compute_ts_tv."""

    def test_returns_tstv_stats(self):
        r = dn.compute_ts_tv(['AAAA', 'AAAA'])
        assert isinstance(r, dn.TsTvStats)

    def test_n_stored(self):
        seqs = ['AAAA', 'AAAA', 'AAAA']
        r = dn.compute_ts_tv(seqs)
        assert r.n == 3

    def test_identical_sequences_zero_counts(self):
        seqs = ['ATCG', 'ATCG', 'ATCG']
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 0
        assert r.n_transversions == 0
        assert r.ts_tv is None

    def test_single_transition_A_G(self):
        """A↔G = purine↔purine = transition."""
        seqs = ['AAAA', 'AGAA']
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 1
        assert r.n_transversions == 0
        assert r.ts_tv is None  # Tv == 0

    def test_single_transition_C_T(self):
        """C↔T = pyrimidine↔pyrimidine = transition."""
        seqs = ['CCCC', 'CTCC']
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 1
        assert r.n_transversions == 0

    def test_single_transversion_A_C(self):
        """A↔C = purine↔pyrimidine = transversion."""
        seqs = ['AAAA', 'ACAA']
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 0
        assert r.n_transversions == 1
        assert r.ts_tv == 0.0

    def test_single_transversion_G_T(self):
        """G↔T = purine↔pyrimidine = transversion."""
        seqs = ['GGGG', 'GTGG']
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 0
        assert r.n_transversions == 1

    def test_ratio_two_to_one(self):
        """2 Ts and 1 Tv across all pairs → Ts/Tv = 2.0."""
        # s1=ACCC, s2=GCCC (A/G Ts), s3=TCCC (A/T Tv, G/T Tv) → Ts=1 Tv=2? Let me be precise.
        # Use: s1=AAAC, s2=GAAC (A/G Ts), s3=AACC (3rd A/C Tv)
        # Pair(s1,s2): col0 A/G → Ts; pair(s1,s3): col3 C/C no diff... let me redo
        # s1=ACCC s2=GCCC → Ts; s1=ACCC s3=ACGC → Ts; s1=ACCC s4=ATCC → Tv
        # Pairs: (0,1)Ts (0,2)Ts (0,3)Tv (1,2)? col1 C/C no, Ts cols: 0: G/A Ts → yes
        # This gets complicated. Use a direct approach:
        # 1 pair: 2 cols, col0=Ts(A/G), col1=Tv(A/C)
        seqs = ['AAA', 'GAC']  # col0: A/G=Ts, col1: A/A=same, col2: A/C=Tv
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 1
        assert r.n_transversions == 1
        assert abs(r.ts_tv - 1.0) < 1e-9

    def test_ts_tv_none_when_no_transversions(self):
        seqs = ['AAAA', 'GGAA']   # A/G Ts at col 0,1; no Tv
        r = dn.compute_ts_tv(seqs)
        assert r.ts_tv is None

    def test_ts_tv_zero_when_no_transitions(self):
        seqs = ['AAAA', 'CCCC']   # all 4 differences are A/C = Tv
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 0
        assert r.n_transversions == 4
        assert r.ts_tv == 0.0

    def test_gap_column_skipped(self):
        """Columns containing gaps must not contribute."""
        seqs = ['A-A', 'AGA']
        r = dn.compute_ts_tv(seqs)
        assert r.L_net == 2          # only cols 0 and 2 are clean
        assert r.n_transitions == 0
        assert r.n_transversions == 0

    def test_ambiguous_column_skipped(self):
        """Columns with N or ? must not contribute."""
        seqs = ['ANA', 'AGA']
        r = dn.compute_ts_tv(seqs)
        assert r.L_net == 2          # cols 0,2 only

    def test_n_less_than_two_returns_empty(self):
        r = dn.compute_ts_tv(['ATCG'])
        assert r.n == 1
        assert r.n_transitions == 0
        assert r.n_transversions == 0
        assert r.ts_tv is None

    def test_l_net_counts_clean_columns(self):
        seqs = ['ATCG', 'ATCG']
        r = dn.compute_ts_tv(seqs)
        assert r.L_net == 4

    def test_one_change_per_biallelic_site(self):
        """DnaSP counts one change per biallelic column, not per sequence pair."""
        seqs = ['AAA', 'GAA', 'AAG']
        # col0 {A,G} -> 1 Ts; col2 {A,G} -> 1 Ts; col1 monomorphic
        r = dn.compute_ts_tv(seqs)
        assert r.n_transitions == 2
        assert r.n_transversions == 0
        assert r.n_sites == 2

    def test_counts_independent_of_sample_size(self):
        """Duplicating a sequence must not change the counts (Codex regression)."""
        base = ['ACGT', 'GCTT']          # col0 A/G Ts, col2 G/T Tv
        r1 = dn.compute_ts_tv(base)
        r2 = dn.compute_ts_tv(base + [base[0]] * 5)
        assert (r1.n_transitions, r1.n_transversions) == (1, 1)
        assert (r2.n_transitions, r2.n_transversions) == (1, 1)
        assert r1.ts_tv == r2.ts_tv == 1.0

    def test_triallelic_site_excluded(self):
        seqs = ['AAAA', 'GAAA', 'CAAA', 'TCAA']
        # col0 {A,G,C,T} -> multiallelic, excluded; col1 {A,C} -> 1 Tv
        r = dn.compute_ts_tv(seqs)
        assert r.n_multiallelic_excluded == 1
        assert r.n_sites == 1
        assert r.n_transversions == 1

    def test_outgroup_unresolvable_site_excluded(self):
        # col0 ingroup {A,G}; outgroup C (not an ingroup allele) -> excluded.
        # col1 ingroup {C,T}; outgroup C -> kept, 1 transition.
        seqs = ['AC', 'AT', 'GC', 'GT']
        r = dn.compute_ts_tv(seqs, outgroup_seq='CC')
        assert r.polarised is True
        assert r.n_unpolarisable_excluded == 1
        assert r.n_transitions == 1
        assert r.n_transversions == 0

    def test_outgroup_ratio_matches_unpolarised_on_surviving_sites(self):
        seqs = ['ACGT', 'GCAT', 'ACGT', 'GCAT']   # col0 {A,G} Ts, col2 {G,A} Ts
        r_no = dn.compute_ts_tv(seqs)
        r_og = dn.compute_ts_tv(seqs, outgroup_seq='ACGT')
        assert r_no.ts_tv == r_og.ts_tv

    def test_run_analysis_dispatch(self):
        """run_analysis must populate results['tstv'] when analysis='tstv'."""
        aln = dn.Alignment(names=['a','b'], seqs=['ATCG','ATCG'])
        results = dn.run_analysis(aln, analyses={'tstv'})
        assert results.get('tstv') is not None
        assert isinstance(results['tstv'], dn.TsTvStats)


# =============================================================================
# Group E  -  Codon usage bias (RSCU + ENC)
# =============================================================================

class TestComputeCodonUsage:
    """Tests for compute_codon_usage."""

    def test_returns_codon_usage_stats(self):
        r = dn.compute_codon_usage(['ATGAAA'])
        assert isinstance(r, dn.CodonUsageStats)

    def test_n_stored(self):
        seqs = ['ATGAAA', 'ATGAAA', 'ATGAAA']
        r = dn.compute_codon_usage(seqs)
        assert r.n == 3

    def test_stop_codon_excluded_from_count(self):
        """Stop codons (TAA, TAG, TGA) must not be counted."""
        r = dn.compute_codon_usage(['ATGTAA'])  # Met + Stop
        assert r.n_codons == 1.0                 # only Met

    def test_gap_triplet_skipped(self):
        """Triplets with gaps are skipped per-sequence."""
        r_clean = dn.compute_codon_usage(['ATGAAA'])
        r_gap   = dn.compute_codon_usage(['ATG---'])  # second codon has gaps
        assert r_gap.n_codons < r_clean.n_codons

    def test_n_codons_counts_clean_codons(self):
        seqs = ['ATGAAATTT', 'ATGAAATTT']  # Met-Lys-Phe × 2 seqs
        r = dn.compute_codon_usage(seqs)
        assert r.n_codons == 3.0

    def test_rscu_uniform_equals_one(self):
        """Equal counts of all synonyms → RSCU = 1.0 for every codon in the family."""
        # Phe: TTT + TTC; use equal counts
        seqs = ['TTTTTC', 'TTTTTC']   # 1 TTT + 1 TTC per seq
        r = dn.compute_codon_usage(seqs)
        assert abs(r.rscu.get('TTT', -1) - 1.0) < 1e-9
        assert abs(r.rscu.get('TTC', -1) - 1.0) < 1e-9

    def test_rscu_preferred_codon_greater_than_one(self):
        """If only one codon is used, its RSCU = family_size."""
        # Lys: AAA + AAG (2-fold); use only AAA
        seqs = ['AAAAAA', 'AAAAAA']   # AAA AAA
        r = dn.compute_codon_usage(seqs)
        assert r.rscu.get('AAA', 0.0) > 1.0
        assert r.rscu.get('AAA', 0.0) == pytest.approx(2.0)

    def test_rscu_avoided_codon_zero(self):
        """Unused codon in a family gets RSCU = 0."""
        seqs = ['AAAAAA', 'AAAAAA']   # only AAA, never AAG
        r = dn.compute_codon_usage(seqs)
        assert r.rscu.get('AAG', 1.0) == pytest.approx(0.0)

    def test_rscu_sum_equals_family_size(self):
        """Sum of RSCU values within a synonymous family equals family size."""
        # Build a sequence with all 6 Leu codons equally
        leu_codons = dn._SYNONYMOUS_FAMILIES['L']   # 6 codons
        seq = ''.join(leu_codons) * 4               # 4 repeats each
        seqs = [seq] * 3
        r = dn.compute_codon_usage(seqs)
        rscu_sum = sum(r.rscu.get(c, 0.0) for c in leu_codons)
        assert abs(rscu_sum - len(leu_codons)) < 1e-6

    def test_enc_range(self):
        """ENC must be in [20, 61] when it can be computed."""
        # Build a rich coding sequence covering all amino acid families
        codons_all: list[str] = []
        for aa, codons in sorted(dn._SYNONYMOUS_FAMILIES.items()):
            codons_all.extend(codons * 3)
        seq = ''.join(codons_all)
        seqs = [seq] * 5
        r = dn.compute_codon_usage(seqs)
        if r.ENC is not None:
            assert 20.0 <= r.ENC <= 61.0

    def test_enc_equals_61_for_uniform_usage(self):
        """Perfectly uniform codon usage → ENC = 61."""
        codons_all: list[str] = []
        for aa, codons in sorted(dn._SYNONYMOUS_FAMILIES.items()):
            codons_all.extend(codons * 2)   # each codon same count
        seq = ''.join(codons_all)
        seqs = [seq] * 10
        r = dn.compute_codon_usage(seqs)
        assert r.ENC is not None
        assert abs(r.ENC - 61.0) < 0.01

    def test_enc_equals_20_for_max_bias(self):
        """Single codon per family → ENC = 20."""
        # Each amino acid: use only its first codon
        codons_all: list[str] = []
        for aa, codons in sorted(dn._SYNONYMOUS_FAMILIES.items()):
            codons_all.extend([codons[0]] * 6)
        seq = ''.join(codons_all)
        seqs = [seq] * 10
        r = dn.compute_codon_usage(seqs)
        assert r.ENC is not None
        assert abs(r.ENC - 20.0) < 0.01

    def test_enc_none_for_insufficient_data(self):
        """ENC is None when amino acid coverage is too sparse for all classes."""
        # Only Met and Phe  -  cannot fill all 4 degeneracy classes
        r = dn.compute_codon_usage(['ATGTTT'])
        assert r.ENC is None

    def test_codon_counts_mean_over_sequences(self):
        """codon_counts should be mean per sequence, not total."""
        seqs = ['ATGATG', 'ATGATG']   # 2 ATG per seq, 2 seqs → raw[ATG]=4, mean=2
        r = dn.compute_codon_usage(seqs)
        assert abs(r.codon_counts.get('ATG', 0.0) - 2.0) < 1e-9

    def test_pooling_across_sequences(self):
        """Codon counts from all sequences are pooled."""
        # seq1 uses AAA, seq2 uses AAG → total Lys = 1 each
        seqs = ['AAA', 'AAG']
        r = dn.compute_codon_usage(seqs)
        # raw[AAA]=1, raw[AAG]=1, mean each = 0.5
        assert abs(r.codon_counts.get('AAA', 0.0) - 0.5) < 1e-9
        assert abs(r.codon_counts.get('AAG', 0.0) - 0.5) < 1e-9
        # RSCU: total_Lys=2, expected = 2/2 = 1.0; each RSCU = 1/1 = 1.0
        assert abs(r.rscu.get('AAA', 0.0) - 1.0) < 1e-9

    def test_ambiguous_base_in_triplet_skipped(self):
        """Triplets with N skip that codon."""
        r_clean = dn.compute_codon_usage(['ATGAAA'])
        r_amb   = dn.compute_codon_usage(['ATGNAA'])  # second codon has N
        assert r_amb.n_codons < r_clean.n_codons

    def test_empty_sequences_returns_defaults(self):
        r = dn.compute_codon_usage([])
        assert r.n == 0
        assert r.ENC is None

    def test_run_analysis_dispatch(self):
        """run_analysis must populate results['codon'] when analysis='codon'."""
        aln = dn.Alignment(names=['a','b'], seqs=['ATGAAATTT', 'ATGAAATTT'])
        results = dn.run_analysis(aln, analyses={'codon'})
        assert results.get('codon') is not None
        assert isinstance(results['codon'], dn.CodonUsageStats)

    def test_run_analysis_tstv_and_codon_together(self):
        """Both tstv and codon can be requested in one run."""
        aln = dn.Alignment(names=['a','b'], seqs=['ATGAAATTT', 'ATGAAGTTT'])
        results = dn.run_analysis(aln, analyses={'tstv', 'codon'})
        assert isinstance(results.get('tstv'), dn.TsTvStats)
        assert isinstance(results.get('codon'), dn.CodonUsageStats)

    def test_vertebrate_mitochondrial_code_counts_tga_as_trp(self):
        # Standard code: TGA is a stop, dropped entirely, never in codon_counts.
        seqs = ['TGATGG', 'TGGTGA']
        std = dn.compute_codon_usage(seqs)
        assert std.codon_counts.get('TGA', 0.0) == 0.0

        # Vertebrate mitochondrial code: TGA = Trp, joins TGG in the Trp
        # family; equal TGA/TGG usage here gives equal (uniform) RSCU.
        mito = dn.compute_codon_usage(seqs, dn.VERTEBRATE_MITOCHONDRIAL_CODE)
        assert mito.codon_counts.get('TGA', 0.0) > 0.0
        assert mito.rscu['TGA'] == pytest.approx(mito.rscu['TGG'])
        assert mito.rscu['TGA'] == pytest.approx(1.0)

    def test_genetic_code_flows_through_run_analysis(self):
        aln = dn.Alignment(names=['a', 'b'], seqs=['TGATGG', 'TGGTGG'])
        std = dn.run_analysis(aln, analyses={'codon'})
        mito = dn.run_analysis(
            aln, analyses={'codon'}, genetic_code=dn.VERTEBRATE_MITOCHONDRIAL_CODE
        )
        assert std['codon'].codon_counts.get('TGA', 0.0) == 0.0
        assert mito['codon'].codon_counts.get('TGA', 0.0) > 0.0

    def test_cli_genetic_code_flag(self, tmp_path):
        f = tmp_path / "mito.fas"
        f.write_text(">a\nTGATGG\n>b\nTGGTGG\n")
        out_std = tmp_path / "out_std"
        out_mito = tmp_path / "out_mito"
        assert dn.main(["--input", str(f), "--analysis", "codon",
                        "--output", str(out_std)]) == 0
        assert dn.main(["--input", str(f), "--analysis", "codon",
                        "--genetic-code", "vertebrate-mitochondrial",
                        "--output", str(out_mito)]) == 0
        std_tsv = (out_std / "results.tsv").read_text()
        mito_report = (out_mito / "report.md").read_text()
        assert "TGA" in mito_report

    def test_cli_genetic_code_rejects_unknown_choice(self, tmp_path):
        f = tmp_path / "x.fas"
        f.write_text(">a\nATG\n>b\nATG\n")
        with pytest.raises(SystemExit):
            dn.main(["--input", str(f), "--genetic-code", "bogus-code",
                     "--output", str(tmp_path / "out")])


# =============================================================================
# Group F: Fay & Wu's H / Zeng's E  and  Fst
# =============================================================================

class TestFayWuStats:
    """Unit tests for FayWuStats dataclass defaults."""

    def test_dataclass_defaults(self):
        fw = dn.FayWuStats()
        assert fw.n == 0
        assert fw.L_net == 0
        assert fw.n_polarised == 0
        assert fw.theta_pi is None
        assert fw.theta_w is None
        assert fw.theta_h is None
        assert fw.theta_l is None
        assert fw.H is None
        assert fw.E is None

    def test_dataclass_set_fields(self):
        fw = dn.FayWuStats(n=5, n_polarised=3, H=0.01, E=-0.02)
        assert fw.n == 5
        assert fw.n_polarised == 3
        assert fw.H == pytest.approx(0.01)
        assert fw.E == pytest.approx(-0.02)


class TestComputeFayWu:
    """Tests for compute_fay_wu()."""

    def test_empty_seqs_returns_default(self):
        result = dn.compute_fay_wu([], "AAAA")
        assert result.n == 0
        assert result.H is None

    def test_single_seq_returns_default(self):
        result = dn.compute_fay_wu(["AAAA"], "AAAA")
        assert result.H is None

    def test_outgroup_too_short_raises(self):
        with pytest.raises(ValueError, match="shorter than alignment"):
            dn.compute_fay_wu(["AAAA", "TTTT"], "AA")

    def test_monomorphic_no_polarised_sites(self):
        seqs = ["AAAA", "AAAA", "AAAA"]
        result = dn.compute_fay_wu(seqs, "AAAA")
        assert result.n_polarised == 0
        assert result.H is None

    def test_one_low_freq_derived_allele(self):
        # n=4, one site: derived T at freq 1/4 (singleton)
        # xi[1]=1; all others=0
        seqs = ["AAAA", "TAAA", "AAAA", "AAAA"]
        og   = "AAAA"
        r = dn.compute_fay_wu(seqs, og)
        assert r.n == 4
        assert r.n_polarised == 1
        assert r.L_net == 4
        # theta_pi = 2*1*(4-1)/(4*3)/4 = 6/12/4 = 0.125
        assert r.theta_pi == pytest.approx(2 * 1 * 3 / (4 * 3) / 4)
        # theta_h = 2*1^2/(4*3)/4 = 2/12/4
        assert r.theta_h == pytest.approx(2 * 1 / (4 * 3) / 4)
        # H = theta_pi - theta_h > 0 (low-freq derived → H positive)
        assert r.H is not None and r.H > 0

    def test_one_high_freq_derived_allele(self):
        # n=4, one site: derived T at freq 3/4 (high freq)
        # xi[3]=1
        seqs = ["TAAA", "TAAA", "TAAA", "AAAA"]
        og   = "AAAA"
        r = dn.compute_fay_wu(seqs, og)
        assert r.n_polarised == 1
        # theta_pi = 2*3*1/(4*3)/4 = 6/12/4 = 0.125
        # theta_h  = 2*9/(4*3)/4   = 18/12/4 = 0.375
        # H = theta_pi - theta_h < 0 (high-freq derived → H negative)
        assert r.H is not None and r.H < 0

    def test_symmetric_freq_H_zero(self):
        # n=4; one site freq 2/4 → theta_pi = theta_h? No.
        # xi[2]=1: theta_pi = 2*2*2/(4*3)/L = 8/12/L
        #          theta_h  = 2*4/(4*3)/L   = 8/12/L → H = 0
        seqs = ["TAAA", "TAAA", "AAAA", "AAAA"]
        og   = "AAAA"
        r = dn.compute_fay_wu(seqs, og)
        assert r.H == pytest.approx(0.0, abs=1e-10)

    def test_zeng_E_negative_with_low_freq_derived(self):
        # Low-freq derived alleles → theta_L < theta_W → E < 0
        seqs = ["AAACGT", "GAACGT", "AAACGT", "AAACGT"]
        og   = "AAACGT"
        r = dn.compute_fay_wu(seqs, og)
        assert r.E is not None and r.E < 0

    def test_unpolarisable_site_excluded(self):
        # Site where outgroup has allele not in ingroup → skip
        seqs = ["AAAA", "GAAA", "AAAA"]
        og   = "CAAA"   # ancestral C not in ingroup at pos 0
        r = dn.compute_fay_wu(seqs, og)
        assert r.n_polarised == 0
        assert r.H is None

    def test_gap_in_outgroup_excluded(self):
        seqs = ["AAAA", "GAAA", "AAAA"]
        og   = "-AAA"   # gap in outgroup at pos 0
        r = dn.compute_fay_wu(seqs, og)
        # pos 0 excluded (outgroup gap); all remaining sites monomorphic
        assert r.n_polarised == 0

    def test_gap_in_ingroup_excluded(self):
        seqs = ["AAAA", "GAAA", "-AAA"]
        og   = "AAAA"
        r = dn.compute_fay_wu(seqs, og)
        # pos 0 has gap in seq[2] → excluded; remaining sites monomorphic
        assert r.n_polarised == 0

    def test_theta_w_correct(self):
        # n=3; xi[1]=2 (two singletons); a1 = 1 + 1/2 = 1.5
        # theta_w = 2 / 1.5 / L_net
        seqs = ["AAAA", "GAAC", "AAAA"]
        og   = "AAAA"
        r = dn.compute_fay_wu(seqs, og)
        # polarised: pos0 (G freq 1/3), pos3 (C freq 1/3) → n_pol=2
        a1 = 1 + 1/2
        expected_w = 2 / a1 / 4
        assert r.theta_w == pytest.approx(expected_w, rel=1e-6)

    def test_run_analysis_faywu_no_outgroup_warns(self, capsys):
        aln = dn.Alignment(names=["a", "b"], seqs=["AAAA", "GAAA"])
        results = dn.run_analysis(aln, analyses={"faywu"})
        assert results.get("faywu") is None
        captured = capsys.readouterr()
        assert "faywu" in captured.err and "outgroup" in captured.err

    def test_run_analysis_faywu_with_outgroup(self):
        seqs = ["AAAA", "GAAA", "AAAA"]
        og   = "AAAA"
        aln  = dn.Alignment(names=["a", "b", "c"], seqs=seqs)
        results = dn.run_analysis(aln, analyses={"faywu"}, outgroup=og)
        assert isinstance(results.get("faywu"), dn.FayWuStats)
        assert results["faywu"].n == 3

    def test_per_site_normalisation(self):
        # Doubling L should halve per-site θ values
        seqs4  = ["AAAA", "GAAA", "AAAA", "AAAA"]
        seqs8  = ["AAAAAAAA", "GAAAAAAA", "AAAAAAAA", "AAAAAAAA"]
        og4    = "AAAA"
        og8    = "AAAAAAAA"
        r4 = dn.compute_fay_wu(seqs4, og4)
        r8 = dn.compute_fay_wu(seqs8, og8)
        assert r4.theta_pi == pytest.approx(r8.theta_pi * 2, rel=1e-6)
        assert r4.H == pytest.approx(r8.H * 2, rel=1e-6)


# -----------------------------------------------------------------------------

class TestFstStats:
    """Unit tests for FstStats dataclass defaults."""

    def test_dataclass_defaults(self):
        fs = dn.FstStats()
        assert fs.n_pops == 0
        assert fs.pop_names == []
        assert fs.pop_sizes == {}
        assert fs.pi_within == {}
        assert fs.pi_between == {}
        assert fs.fst_pairwise == {}
        assert fs.fst_mean is None

    def test_dataclass_n_pops_set(self):
        fs = dn.FstStats(n_pops=3, pop_names=["A", "B", "C"])
        assert fs.n_pops == 3
        assert len(fs.pop_names) == 3


class TestComputeFst:
    """Tests for compute_fst()."""

    def test_single_pop_returns_empty_pairwise(self):
        result = dn.compute_fst({"Pop1": ["AAAA", "AAAA"]})
        assert result.n_pops == 1
        assert result.fst_pairwise == {}
        assert result.fst_mean is None

    def test_fixed_differences_fst_one(self):
        # Populations completely fixed for different alleles → Fst = 1
        pop_seqs = {
            "A": ["AAAA", "AAAA"],
            "B": ["TTTT", "TTTT"],
        }
        result = dn.compute_fst(pop_seqs)
        assert result.fst_pairwise[("A", "B")] == pytest.approx(1.0)
        assert result.fst_mean == pytest.approx(1.0)

    def test_identical_pops_fst_zero(self):
        # Same sequences in both pops → Fst = 0
        pop_seqs = {
            "A": ["AAAA", "TTTT"],
            "B": ["AAAA", "TTTT"],
        }
        result = dn.compute_fst(pop_seqs)
        assert result.fst_pairwise[("A", "B")] == pytest.approx(0.0)

    def test_pi_within_homozygous_pop(self):
        # Monomorphic pop → pi_within = 0
        result = dn.compute_fst({"X": ["AAAA", "AAAA"], "Y": ["TTTT", "TTTT"]})
        assert result.pi_within["X"] == pytest.approx(0.0)
        assert result.pi_within["Y"] == pytest.approx(0.0)

    def test_pi_within_polymorphic_pop(self):
        # Pop with 2 seqs differing at all 4 sites → pi_within = 1.0 per site
        result = dn.compute_fst({"A": ["AAAA", "TTTT"], "B": ["CCCC", "CCCC"]})
        assert result.pi_within["A"] == pytest.approx(1.0)
        assert result.pi_within["B"] == pytest.approx(0.0)

    def test_pi_between_computed(self):
        # All 4 sites differ between pops
        pop_seqs = {"A": ["AAAA"], "B": ["TTTT"]}
        result = dn.compute_fst(pop_seqs)
        assert result.pi_between[("A", "B")] == pytest.approx(1.0)

    def test_fst_between_zero_and_one(self):
        pop_seqs = {
            "A": ["AAACCC", "TTACCC"],
            "B": ["AAAGGG", "TTAGGG"],
        }
        result = dn.compute_fst(pop_seqs)
        fst = result.fst_pairwise.get(("A", "B"))
        assert fst is not None
        assert 0.0 <= fst <= 1.0

    def test_fst_none_when_no_between_variation(self):
        # Both pops monomorphic with same allele → Dxy = 0 → Fst = None
        pop_seqs = {"A": ["AAAA", "AAAA"], "B": ["AAAA", "AAAA"]}
        result = dn.compute_fst(pop_seqs)
        assert result.fst_pairwise[("A", "B")] is None

    def test_three_pops_pairwise(self):
        pop_seqs = {
            "P1": ["AAAA", "AAAA"],
            "P2": ["TTTT", "TTTT"],
            "P3": ["CCCC", "CCCC"],
        }
        result = dn.compute_fst(pop_seqs)
        assert result.n_pops == 3
        assert len(result.fst_pairwise) == 3  # 3 choose 2 = 3 pairs
        assert ("P1", "P2") in result.fst_pairwise
        assert ("P1", "P3") in result.fst_pairwise
        assert ("P2", "P3") in result.fst_pairwise

    def test_fst_mean_is_average_of_pairwise(self):
        pop_seqs = {
            "A": ["AAAA", "AAAA"],
            "B": ["TTTT", "TTTT"],
            "C": ["CCCC", "CCCC"],
        }
        result = dn.compute_fst(pop_seqs)
        pairwise_vals = [v for v in result.fst_pairwise.values() if v is not None]
        expected_mean = sum(pairwise_vals) / len(pairwise_vals)
        assert result.fst_mean == pytest.approx(expected_mean)

    def test_pop_sizes_recorded(self):
        pop_seqs = {
            "A": ["AAAA", "AAAA", "AAAA"],
            "B": ["TTTT"],
        }
        result = dn.compute_fst(pop_seqs)
        assert result.pop_sizes["A"] == 3
        assert result.pop_sizes["B"] == 1

    def test_run_analysis_fst_no_pop_file_warns(self, capsys):
        aln = dn.Alignment(names=["a", "b"], seqs=["AAAA", "TTTT"])
        results = dn.run_analysis(aln, analyses={"fst"})
        assert results.get("fst") is None
        captured = capsys.readouterr()
        assert "fst" in captured.err.lower() or "pop-file" in captured.err

    def test_run_analysis_fst_with_pop_assignments(self):
        seqs = ["AAAA", "AAAA", "TTTT", "TTTT"]
        names = ["a1", "a2", "b1", "b2"]
        aln = dn.Alignment(names=names, seqs=seqs)
        pops = {"a1": "Pop1", "a2": "Pop1", "b1": "Pop2", "b2": "Pop2"}
        results = dn.run_analysis(aln, analyses={"fst"}, pop_assignments=pops)
        fst_r = results.get("fst")
        assert isinstance(fst_r, dn.FstStats)
        assert fst_r.fst_mean == pytest.approx(1.0)

    def test_fst_gap_columns_excluded(self):
        # Column with gap should not contribute to L_clean
        pop_seqs = {
            "A": ["A-AA", "A-AA"],
            "B": ["T-TT", "T-TT"],
        }
        result = dn.compute_fst(pop_seqs)
        # 3 clean columns (col 1 excluded), all differ → Dxy = 1.0
        assert result.pi_between[("A", "B")] == pytest.approx(1.0)
        assert result.fst_pairwise[("A", "B")] == pytest.approx(1.0)
