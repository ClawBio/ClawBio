"""Tests for reading the files, pairing the variants, and the run as a whole.

The arithmetic is tested in `test_coloc_abf_math.py`. What is tested here is
everything around it: the column aliases, which variant a row belongs to, what
happens to a row that cannot be used, the two demos, the output contract, and
the reproducibility bundle.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import coloc_abf_region_run as skill  # noqa: E402
from coloc_abf_math import ColocError  # noqa: E402

SCRIPT = SKILL_DIR / "coloc_abf_region_run.py"

HEADER = "variant_id\tchromosome\tposition\tref\talt\tbeta\tse\tp_value\teaf\tn"


def _write_tsv(path: Path, rows: list[str], header: str = HEADER) -> Path:
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def _row(vid, pos, ref, alt, beta, se, eaf=0.3, n=10_000, chrom="1"):
    return f"{vid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{beta}\t{se}\t0.1\t{eaf}\t{n}"


# ---------------------------------------------------------------------------
# Reading a file
# ---------------------------------------------------------------------------

class TestLoading:

    def test_reads_the_canonical_columns(self, tmp_path):
        path = _write_tsv(tmp_path / "a.tsv", [
            _row("rs1", 100, "A", "G", 0.1, 0.02),
            _row("rs2", 200, "C", "T", -0.2, 0.03),
        ])
        rows = skill.load_regional_sumstats(path)
        assert [r.variant_id for r in rows] == ["rs1", "rs2"]
        assert rows[0].beta == pytest.approx(0.1)
        assert rows[1].se == pytest.approx(0.03)
        assert rows[0].n == pytest.approx(10_000)

    def test_reads_common_column_aliases(self, tmp_path):
        """A file straight out of a region fetch or a harmoniser, unrewritten."""
        path = tmp_path / "aliased.tsv"
        path.write_text(
            "SNP\tCHR\tBP\tA2\tA1\tBETA\tSE\tP\tFREQ\tN\n"
            "rs1\t1\t100\tA\tG\t0.1\t0.02\t0.4\t0.31\t9000\n",
            encoding="utf-8",
        )
        rows = skill.load_regional_sumstats(path)
        assert rows[0].variant_id == "rs1"
        assert (rows[0].ref, rows[0].alt) == ("A", "G")
        assert rows[0].eaf == pytest.approx(0.31)
        assert rows[0].p_value == pytest.approx(0.4)

    def test_reads_a_comma_separated_file(self, tmp_path):
        path = tmp_path / "a.csv"
        path.write_text(
            "variant_id,chromosome,position,ref,alt,beta,se\n"
            "rs1,1,100,A,G,0.1,0.02\n",
            encoding="utf-8",
        )
        assert len(skill.load_regional_sumstats(path)) == 1

    def test_a_missing_required_column_names_what_is_missing(self, tmp_path):
        path = tmp_path / "a.tsv"
        path.write_text("variant_id\tchromosome\tposition\tbeta\tse\nrs1\t1\t1\t0.1\t0.02\n",
                        encoding="utf-8")
        with pytest.raises(skill.ColocInputError, match=r"ref.*alt|alt.*ref"):
            skill.load_regional_sumstats(path)

    def test_na_and_empty_values_become_none(self, tmp_path):
        path = _write_tsv(tmp_path / "a.tsv", [
            _row("rs1", 100, "A", "G", "NA", 0.02),
            _row("rs2", 200, "C", "T", 0.1, ""),
        ])
        rows = skill.load_regional_sumstats(path)
        assert rows[0].beta is None
        assert rows[1].se is None

    def test_a_missing_file_says_so(self, tmp_path):
        with pytest.raises(skill.ColocInputError, match="not found"):
            skill.load_regional_sumstats(tmp_path / "absent.tsv")


# ---------------------------------------------------------------------------
# Which variant a row belongs to
# ---------------------------------------------------------------------------

class TestVariantIdentity:

    def test_one_variant_written_in_either_allele_order_is_one_variant(self):
        forward = skill.SumstatRow("rs1", "1", 100, "A", "G", 0.1, 0.02, None, None, None)
        reverse = skill.SumstatRow("rs1", "chr1", 100, "G", "A", -0.1, 0.02, None, None, None)
        assert forward.identity == reverse.identity

    def test_two_alleles_at_one_position_are_two_variants(self):
        a = skill.SumstatRow("rs1", "1", 100, "A", "G", 0.1, 0.02, None, None, None)
        b = skill.SumstatRow("rs1", "1", 100, "A", "T", 0.1, 0.02, None, None, None)
        assert a.identity != b.identity

    def test_a_chr_prefix_does_not_make_a_different_variant(self):
        a = skill.SumstatRow("rs1", "1", 100, "A", "G", 0.1, 0.02, None, None, None)
        b = skill.SumstatRow("rs1", "CHR1", 100, "A", "G", 0.1, 0.02, None, None, None)
        assert a.identity == b.identity

    def test_lower_case_alleles_are_the_same_variant(self):
        a = skill.SumstatRow("rs1", "1", 100, "A", "G", 0.1, 0.02, None, None, None)
        b = skill.SumstatRow("rs1", "1", 100, "a", "g", 0.1, 0.02, None, None, None)
        assert a.identity == b.identity


# ---------------------------------------------------------------------------
# Pairing, and what happens to a row that cannot be used
# ---------------------------------------------------------------------------

class TestPairing:

    @staticmethod
    def _rows(spec):
        return [skill.SumstatRow(*s) for s in spec]

    def test_pairs_on_the_variant_not_on_the_identifier(self):
        """Two files naming the same variant differently still pair.

        One carries an rsID and the alleles one way round; the other carries a
        positional identifier and the alleles the other way round.
        """
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs2", "1", 200, "C", "T", 0.2, 0.02, None, 0.3, 500),
        ])
        right = self._rows([
            ("1_100_G_A", "chr1", 100, "G", "A", -0.1, 0.03, None, 0.3, 500),
            ("1_200_C_T", "chr1", 200, "C", "T", 0.2, 0.03, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, right)
        assert len(paired.keys) == 2
        assert paired.excluded == []
        # The identifiers come back in the first trait's terms.
        assert paired.display_ids == ["rs1", "rs2"]

    def test_a_variant_in_only_one_trait_is_recorded_not_dropped_silently(self):
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs_only_left", "1", 300, "A", "G", 0.1, 0.02, None, 0.3, 500),
        ])
        right = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs_only_right", "1", 400, "A", "G", 0.1, 0.02, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, right)
        assert paired.keys == ["1_100_A_G"]
        reasons = {e["variant_id"]: e["reason"] for e in paired.excluded}
        assert reasons == {
            "rs_only_left": skill.EXCLUDED_NOT_SHARED,
            "rs_only_right": skill.EXCLUDED_NOT_SHARED,
        }

    @pytest.mark.parametrize("beta,se,reason", [
        (None, 0.02, skill.EXCLUDED_MISSING_STATS),
        (0.1, None, skill.EXCLUDED_MISSING_STATS),
        (0.1, 0.0, skill.EXCLUDED_INVALID_SE),
        (0.1, -0.02, skill.EXCLUDED_INVALID_SE),
        (float("nan"), 0.02, skill.EXCLUDED_INVALID_SE),
    ])
    def test_an_unusable_row_is_recorded_with_its_reason(self, beta, se, reason):
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs_bad", "1", 200, "A", "G", beta, se, None, 0.3, 500),
        ])
        right = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs_bad", "1", 200, "A", "G", 0.1, 0.02, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, right)
        by_variant = {(e["trait"], e["variant_id"]): e["reason"] for e in paired.excluded}
        assert by_variant[("trait1", "rs_bad")] == reason

    def test_a_row_with_no_alleles_is_recorded(self):
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs_noalleles", "1", 200, "", "", 0.1, 0.02, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, left[:1])
        assert any(e["reason"] == skill.EXCLUDED_MISSING_ALLELES
                   for e in paired.excluded)

    def test_a_duplicated_row_is_dropped_and_recorded(self):
        """One variant written twice with the same numbers is a duplicated row."""
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs1_again", "1", 100, "G", "A", 0.1, 0.02, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, left[:1])
        assert paired.keys == ["1_100_A_G"]
        assert [e["reason"] for e in paired.excluded
                if e["trait"] == "trait1"] == [skill.EXCLUDED_DUPLICATE_ROW]

    def test_two_rows_for_one_variant_that_disagree_stop_the_run(self):
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs1_other", "1", 100, "A", "G", 0.9, 0.05, None, 0.3, 500),
        ])
        with pytest.raises(ColocError, match="cannot choose"):
            skill.pair_region(left, left[:1])

    def test_kept_and_excluded_account_for_every_row_read(self):
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs2", "1", 200, "C", "T", 0.2, 0.02, None, 0.3, 500),
            ("rs_bad", "1", 300, "A", "G", None, 0.02, None, 0.3, 500),
            ("rs_left_only", "1", 400, "A", "G", 0.1, 0.02, None, 0.3, 500),
        ])
        right = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.03, None, 0.3, 500),
            ("rs2", "1", 200, "C", "T", 0.2, 0.03, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, right)
        for trait, n_read in (("trait1", 4), ("trait2", 2)):
            n_excluded = sum(1 for e in paired.excluded if e["trait"] == trait)
            assert len(paired.keys) + n_excluded == n_read

    def test_a_multi_allelic_site_under_one_rsid_is_two_variants_and_keeps_pairing(self):
        """One rsID, two alleles: the alleles decide, and the key is reported.

        Pairing by identifier alone would match the wrong allele here; pairing by
        alleles gives two variants, so the reported identifier falls back to the
        positional key, which names exactly one of them.
        """
        left = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.02, None, 0.3, 500),
            ("rs1", "1", 100, "A", "T", 0.8, 0.02, None, 0.3, 500),
        ])
        right = self._rows([
            ("rs1", "1", 100, "A", "G", 0.1, 0.03, None, 0.3, 500),
            ("rs1", "1", 100, "A", "T", 0.8, 0.03, None, 0.3, 500),
        ])
        paired = skill.pair_region(left, right)
        assert sorted(paired.keys) == ["1_100_A_G", "1_100_A_T"]
        assert sorted(paired.display_ids) == ["1_100_A_G", "1_100_A_T"]


# ---------------------------------------------------------------------------
# One run, end to end
# ---------------------------------------------------------------------------

class TestRunPair:

    @staticmethod
    def _pair(peak_left=100, peak_right=100, n=120):
        """Two traits over one window, each peaking at the index given."""
        import math
        left, right = [], []
        for i in range(n):
            pos = 100_000 + 500 * i
            z_l = 9.0 * math.exp(-abs(i - peak_left) / 6.0)
            z_r = 8.0 * math.exp(-abs(i - peak_right) / 6.0)
            left.append(skill.SumstatRow(f"rs{i}", "1", pos, "A", "G",
                                         z_l * 0.03, 0.03, None, 0.3, 500))
            right.append(skill.SumstatRow(f"rs{i}", "1", pos, "A", "G",
                                          z_r * 0.02, 0.02, None, 0.3, 250_000))
        return left, right

    def test_a_shared_signal_comes_out_as_h4(self):
        left, right = self._pair(100, 100)
        run = skill.run_pair(left, right, trait1_type="quant", trait2_type="quant",
                             label1="expression", label2="a trait",
                             sdY1=1.0, sdY2=1.0)
        assert run["result"].h4 > 0.9
        assert run["lead_display_id"] == "rs100"
        assert run["sensitivity"].robust

    def test_distinct_signals_come_out_as_h3_with_a_caveat(self):
        left, right = self._pair(40, 100)
        run = skill.run_pair(left, right, trait1_type="quant", trait2_type="quant",
                             label1="expression", label2="a trait",
                             sdY1=1.0, sdY2=1.0)
        assert run["result"].h3 > 0.9
        assert any("different causal variants" in c for c in run["caveats"])

    def test_the_result_does_not_depend_on_which_allele_was_reported(self):
        """The end-to-end form of the property in test_coloc_abf_math.py.

        Flipping the effect allele of every row in one file, as a file reporting
        the other allele would, leaves the posteriors unchanged.
        """
        import dataclasses
        left, right = self._pair(100, 100)
        flipped = [
            dataclasses.replace(r, ref=r.alt, alt=r.ref, beta=-r.beta, eaf=1 - r.eaf)
            for r in right
        ]
        kwargs = dict(trait1_type="quant", trait2_type="quant",
                      label1="expression", label2="a trait", sdY1=1.0, sdY2=1.0)
        a = skill.run_pair(left, right, **kwargs)["result"].pp
        b = skill.run_pair(left, flipped, **kwargs)["result"].pp
        assert a == b

    def test_estimating_sdy_is_flagged_in_the_caveats(self):
        left, right = self._pair(100, 100)
        run = skill.run_pair(left, right, trait1_type="quant", trait2_type="quant",
                             label1="expression", label2="a trait",
                             sdY2=1.0, n1=500)
        assert run["result"].sdy_estimated["trait1"] is True
        assert any("sdY was estimated" in c for c in run["caveats"])

    def test_a_quantitative_trait_with_no_sdy_and_no_frequency_says_what_is_missing(self):
        import dataclasses
        left, right = self._pair(100, 100)
        left = [dataclasses.replace(r, eaf=None) for r in left]
        with pytest.raises(ColocError, match="allele frequency"):
            skill.run_pair(left, right, trait1_type="quant", trait2_type="quant",
                           label1="expression", label2="a trait", sdY2=1.0, n1=500)


# ---------------------------------------------------------------------------
# The demos
# ---------------------------------------------------------------------------

class TestDemos:

    def test_both_demos_are_listed(self):
        assert set(skill.available_demos()) == {
            "sort1_eqtl_x_ldl_gwas", "distinct_signals_negative"}

    def test_the_shared_demo_colocalises_robustly(self, tmp_path):
        assert skill.main(["--demo", "--output", str(tmp_path)]) == 0
        payload = json.loads((tmp_path / "result.json").read_text())
        assert payload["posteriors"]["PP.H4"] > 0.9
        assert payload["sensitivity"]["robust"] is True
        assert payload["n_shared_variants"] == 200
        assert payload["excluded_variants"]["total"] == 0

    def test_the_distinct_demo_does_not_colocalise(self, tmp_path):
        assert skill.main(["--demo", "distinct_signals_negative",
                           "--output", str(tmp_path)]) == 0
        payload = json.loads((tmp_path / "result.json").read_text())
        assert payload["posteriors"]["PP.H3"] > 0.9
        assert payload["sensitivity"]["shared_signal_at_default_p12"] is False

    def test_an_unknown_demo_name_exits_non_zero(self, tmp_path):
        assert skill.main(["--demo", "nope", "--output", str(tmp_path)]) == 2

    def test_list_demos_prints_them(self, capsys):
        assert skill.main(["--list-demos"]) == 0
        printed = capsys.readouterr().out
        assert "sort1_eqtl_x_ldl_gwas" in printed
        assert "distinct_signals_negative" in printed

    def test_the_demo_data_is_reproducible(self, tmp_path):
        """Regenerating the committed tables must reproduce them byte for byte.

        This is what makes the committed demo data auditable rather than
        something that has to be taken on trust.
        """
        import importlib.util
        maker = SKILL_DIR / "examples" / "make_demo_data.py"
        spec = importlib.util.spec_from_file_location("coloc_abf_make_demo_data", maker)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for regenerated in module.write_demo_data(tmp_path):
            committed = SKILL_DIR / "examples" / "sumstats" / regenerated.name
            assert regenerated.read_bytes() == committed.read_bytes(), (
                f"{regenerated.name} differs from the committed copy; re-run "
                "examples/make_demo_data.py and commit the result"
            )


# ---------------------------------------------------------------------------
# The CLI
# ---------------------------------------------------------------------------

class TestCli:

    def test_the_script_runs_as_a_subprocess(self, tmp_path):
        done = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True,
        )
        assert done.returncode == 0, done.stderr
        assert "PP.H4=" in done.stdout

    def test_a_config_with_no_trait_block_exits_two(self, tmp_path):
        config = tmp_path / "bad.yaml"
        config.write_text("trait1: {sumstats: a.tsv, trait_type: quant}\n", encoding="utf-8")
        assert skill.main(["--input", str(config), "--output", str(tmp_path / "out")]) == 2

    def test_a_config_with_a_bad_trait_type_exits_two(self, tmp_path):
        config = tmp_path / "bad.yaml"
        config.write_text(
            "trait1: {sumstats: a.tsv, trait_type: continuous}\n"
            "trait2: {sumstats: b.tsv, trait_type: quant}\n", encoding="utf-8")
        assert skill.main(["--input", str(config), "--output", str(tmp_path / "out")]) == 2

    def test_a_json_config_is_accepted(self, tmp_path):
        demo_data = SKILL_DIR / "examples" / "sumstats"
        config = tmp_path / "config.json"
        config.write_text(json.dumps({
            "trait1": {"sumstats": str(demo_data / "synthetic_eqtl_shared.tsv"),
                       "trait_type": "quant", "sdY": 1.0, "label": "expression"},
            "trait2": {"sumstats": str(demo_data / "synthetic_gwas_shared.tsv"),
                       "trait_type": "quant", "sdY": 1.0, "label": "a trait"},
        }), encoding="utf-8")
        out = tmp_path / "out"
        assert skill.main(["--input", str(config), "--output", str(out)]) == 0
        assert json.loads((out / "result.json").read_text())["posteriors"]["PP.H4"] > 0.9

    def test_the_rule_can_be_set_from_the_config(self, tmp_path):
        demo_data = SKILL_DIR / "examples" / "sumstats"
        config = tmp_path / "config.json"
        config.write_text(json.dumps({
            "trait1": {"sumstats": str(demo_data / "synthetic_eqtl_shared.tsv"),
                       "trait_type": "quant", "sdY": 1.0},
            "trait2": {"sumstats": str(demo_data / "synthetic_gwas_shared.tsv"),
                       "trait_type": "quant", "sdY": 1.0},
            "rule_h4_min": 0.999,
        }), encoding="utf-8")
        out = tmp_path / "out"
        assert skill.main(["--input", str(config), "--output", str(out)]) == 0
        payload = json.loads((out / "result.json").read_text())
        assert payload["sensitivity"]["rule"] == "PP.H4 >= 0.999"

    def test_neither_input_nor_demo_is_an_argument_error(self):
        with pytest.raises(SystemExit):
            skill.main(["--output", "/tmp/whatever"])


# ---------------------------------------------------------------------------
# The output contract
# ---------------------------------------------------------------------------

def _parse_output_contract(skill_md: Path) -> list[str]:
    """Extract the files promised in the SKILL.md '## Output Structure' tree.

    Returns output-relative paths. Directory lines, and any entry whose inline
    comment contains 'optional', are skipped.
    """
    if not skill_md.exists():
        return []
    match = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```",
                      skill_md.read_text(), re.S)
    if not match:
        return []
    files: list[str] = []
    parents: dict[int, str] = {}
    for raw in match.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        prefix, name = re.match(r"^([\s│├└─]*)(.*)$", entry).groups()
        name = name.strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for deeper in [k for k in parents if k > depth]:
                del parents[deeper]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(f"{rel}/{name}" if rel else name)
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md must actually be written."""

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable '## Output Structure' section"
        assert skill.main(["--demo", "--output", str(tmp_path)]) == 0
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md promises outputs the skill did not write: " + ", ".join(missing)
        )

    def test_the_report_holds_no_absolute_output_path(self, tmp_path):
        skill.main(["--demo", "--output", str(tmp_path)])
        assert str(tmp_path) not in (tmp_path / "report.md").read_text()

    def test_the_report_carries_the_disclaimer(self, tmp_path):
        skill.main(["--demo", "--output", str(tmp_path)])
        assert skill.DISCLAIMER in (tmp_path / "report.md").read_text()

    def test_two_runs_give_the_same_result_apart_from_the_timestamp(self, tmp_path):
        a, b = tmp_path / "a", tmp_path / "b"
        skill.main(["--demo", "--output", str(a)])
        skill.main(["--demo", "--output", str(b)])
        first = json.loads((a / "result.json").read_text())
        second = json.loads((b / "result.json").read_text())
        first.pop("generated_utc")
        second.pop("generated_utc")
        assert first == second

    def test_the_checksums_match_the_files_they_name(self, tmp_path):
        skill.main(["--demo", "--output", str(tmp_path)])
        lines = (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines()
        assert lines
        for line in lines:
            digest, label = line.split("  ", 1)
            target = tmp_path / label
            assert target.exists(), f"checksums name a file that is not there: {label}"
            assert hashlib.sha256(target.read_bytes()).hexdigest() == digest

    def test_the_tables_are_readable_and_complete(self, tmp_path):
        skill.main(["--demo", "--output", str(tmp_path)])
        with (tmp_path / "tables" / "posteriors.tsv").open(newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [r["hypothesis"] for r in rows] == [
            "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
        with (tmp_path / "tables" / "sensitivity_p12.tsv").open(newline="") as fh:
            sweep = list(csv.DictReader(fh, delimiter="\t"))
        assert len(sweep) >= 9
        assert all(r["rule_passes"] in ("yes", "no") for r in sweep)

    def test_the_excluded_table_has_a_header_even_when_nothing_was_excluded(self, tmp_path):
        skill.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "tables" / "excluded_variants.tsv").read_text()
        assert text.splitlines()[0].split("\t") == [
            "trait", "variant_id", "chromosome", "position", "ref", "alt", "reason"]

    def test_every_excluded_row_reaches_the_table(self, tmp_path):
        """A dropped variant is a row with a reason, not a number in a log."""
        data = tmp_path / "inputs"
        data.mkdir()
        left = _write_tsv(data / "left.tsv", [
            _row("rs1", 100, "A", "G", 0.1, 0.02),
            _row("rs_left_only", 300, "A", "G", 0.1, 0.02),
            _row("rs_bad", 400, "A", "G", "NA", 0.02),
        ])
        right = _write_tsv(data / "right.tsv", [
            _row("rs1", 100, "A", "G", 0.1, 0.02),
            _row("rs_right_only", 500, "A", "G", 0.1, 0.02),
        ])
        config = tmp_path / "config.json"
        config.write_text(json.dumps({
            "trait1": {"sumstats": left.name, "trait_type": "quant", "sdY": 1.0},
            "trait2": {"sumstats": right.name, "trait_type": "quant", "sdY": 1.0},
        }), encoding="utf-8")
        (data / "config.json").write_text(config.read_text(), encoding="utf-8")
        out = tmp_path / "out"
        assert skill.main(["--input", str(data / "config.json"), "--output", str(out)]) == 0
        with (out / "tables" / "excluded_variants.tsv").open(newline="") as fh:
            excluded = list(csv.DictReader(fh, delimiter="\t"))
        assert len(excluded) == 3
        assert {r["reason"] for r in excluded} == {
            skill.EXCLUDED_NOT_SHARED, skill.EXCLUDED_MISSING_STATS}
        payload = json.loads((out / "result.json").read_text())
        assert payload["excluded_variants"]["total"] == 3
        assert sum(payload["excluded_variants"]["by_reason"].values()) == 3
