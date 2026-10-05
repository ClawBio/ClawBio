"""Tests for isolate-amr-typing. Written first (red/green TDD)."""

from __future__ import annotations

import gzip
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "isolate_amr_typing.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO_FASTA = EXAMPLES / "demo_isolate.fasta"
DEMO_AMR = EXAMPLES / "demo_amrfinder.tsv"
DEMO_MLST = EXAMPLES / "demo_mlst.tsv"
DEMO_PLASMID = EXAMPLES / "demo_plasmidfinder.tsv"

AMR_V4_HEADER = (
    "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\t"
    "Subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t"
    "% Coverage of reference\t% Identity to reference\tAlignment length\t"
    "Closest reference accession\tClosest reference name\tHMM accession\tHMM description"
)
AMR_V3_HEADER = (
    "Protein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\t"
    "Element type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\t"
    "Reference sequence length\t% Coverage of reference sequence\t"
    "% Identity to reference sequence\tAlignment length\tAccession of closest sequence\t"
    "Name of closest sequence\tHMM id\tHMM description"
)


def _load():
    spec = importlib.util.spec_from_file_location("isolate_amr_typing", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


def _amr_row(symbol, method, *, contig="c1", typ="AMR", subtype="AMR", cls="BETA-LACTAM",
             subclass="CEPHALOSPORIN", scope="core", cov="100.00", ident="100.00"):
    return "\t".join([
        "NA", contig, "10", "900", "+", symbol, f"{symbol} product", scope, typ, subtype, cls,
        subclass, method, "297", "297", cov, ident, "297", "SYNTHETIC", symbol, "NA", "NA",
    ])


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                          capture_output=True, text=True)


# ── Assembly ───────────────────────────────────────────────────────────────────


class TestAssembly:
    def test_demo_assembly_stats(self, mod):
        stats = mod.read_assembly(DEMO_FASTA)
        assert stats["n_contigs"] == 4
        assert stats["total_length"] == 20300
        assert stats["n50"] == 5600  # 9000 alone is under half of 20300
        assert stats["contigs"]["contig_2"] == 4200
        assert 45 < stats["gc_percent"] < 55

    def test_contig_names_stop_at_whitespace(self, mod):
        assert "contig_1" in mod.read_assembly(DEMO_FASTA)["contigs"]

    def test_gzipped_assembly(self, mod, tmp_path):
        gz = tmp_path / "asm.fasta.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(DEMO_FASTA.read_text())
        assert mod.read_assembly(gz)["total_length"] == 20300

    def test_n50_is_length_weighted(self, mod, tmp_path):
        fa = _write(tmp_path / "a.fa", ">a\n" + "A" * 10 + "\n>b\n" + "C" * 6 + "\n>c\n" + "G" * 4 + "\n")
        assert mod.read_assembly(fa)["n50"] == 10

    def test_empty_file_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="no sequences"):
            mod.read_assembly(_write(tmp_path / "e.fa", ""))

    def test_non_fasta_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="FASTA"):
            mod.read_assembly(_write(tmp_path / "x.fa", "@read1\nACGT\n+\nIIII\n"))

    def test_protein_fasta_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="nucleotide"):
            mod.read_assembly(_write(tmp_path / "p.faa", ">p1\nMKLVFWQEDHRPSTYNI\n"))

    def test_duplicate_contig_names_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            mod.read_assembly(_write(tmp_path / "d.fa", ">a\nACGT\n>a\nACGT\n"))


# ── AMRFinderPlus ──────────────────────────────────────────────────────────────


class TestAmrfinderParser:
    def test_demo_table(self, mod):
        hits = mod.parse_amrfinder(DEMO_AMR)
        assert len(hits) == 16
        ctx = next(h for h in hits if h["symbol"] == "blaCTX-M-15")
        assert ctx["contig"] == "contig_2"
        assert ctx["call"] == "allele"
        assert ctx["element_type"] == "AMR"
        assert ctx["drug_class"] == "BETA-LACTAM"
        assert ctx["identity"] == 100.0 and ctx["coverage"] == 100.0

    def test_v3_and_v4_headers_parse_identically(self, mod, tmp_path):
        row = _amr_row("blaKPC-2", "ALLELEX")
        v4 = mod.parse_amrfinder(_write(tmp_path / "v4.tsv", AMR_V4_HEADER + "\n" + row + "\n"))
        v3 = mod.parse_amrfinder(_write(tmp_path / "v3.tsv", AMR_V3_HEADER + "\n" + row + "\n"))
        assert v3 == v4
        assert v4[0]["symbol"] == "blaKPC-2"

    def test_leading_name_column_is_tolerated(self, mod, tmp_path):
        text = "Name\t" + AMR_V4_HEADER + "\nisoA\t" + _amr_row("sul2", "EXACTX") + "\n"
        hits = mod.parse_amrfinder(_write(tmp_path / "n.tsv", text))
        assert hits[0]["symbol"] == "sul2"

    def test_multiple_samples_rejected(self, mod, tmp_path):
        text = ("Name\t" + AMR_V4_HEADER + "\nisoA\t" + _amr_row("sul2", "EXACTX")
                + "\nisoB\t" + _amr_row("sul1", "EXACTX") + "\n")
        with pytest.raises(ValueError, match="one isolate"):
            mod.parse_amrfinder(_write(tmp_path / "m.tsv", text))

    def test_header_only_means_zero_hits(self, mod, tmp_path):
        assert mod.parse_amrfinder(_write(tmp_path / "h.tsv", AMR_V4_HEADER + "\n")) == []

    def test_unrecognised_header_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="AMRFinderPlus"):
            mod.parse_amrfinder(_write(tmp_path / "bad.tsv", "gene\tcontig\nblaTEM\tc1\n"))

    def test_empty_file_rejected(self, mod, tmp_path):
        # An empty file is a failed run, not "no resistance genes".
        with pytest.raises(ValueError, match="empty"):
            mod.parse_amrfinder(_write(tmp_path / "empty.tsv", ""))

    @pytest.mark.parametrize("method,call", [
        ("ALLELEX", "allele"), ("ALLELEP", "allele"),
        ("EXACTX", "exact"), ("EXACTP", "exact"),
        ("BLASTX", "blast"), ("BLASTP", "blast"),
        ("PARTIALX", "partial"), ("PARTIALP", "partial"),
        ("PARTIAL_CONTIG_ENDX", "partial_contig_end"), ("PARTIAL_CONTIG_ENDP", "partial_contig_end"),
        ("INTERNAL_STOP", "internal_stop"),
        ("HMM", "hmm"),
        ("POINTX", "point"), ("POINTN", "point"), ("POINTP", "point"),
        ("SOMETHING_NEW", "other"),
    ])
    def test_method_to_call(self, mod, method, call):
        assert mod.classify_method(method) == call

    @pytest.mark.parametrize("call,tier", [
        ("allele", "confident"), ("exact", "confident"), ("blast", "confident"),
        ("point", "confident"),
        ("partial", "review"), ("partial_contig_end", "review"), ("hmm", "review"),
        ("other", "review"),
        ("internal_stop", "disrupted"),
    ])
    def test_call_to_tier(self, mod, call, tier):
        assert mod.CALL_TIER[call] == tier


# ── MLST ───────────────────────────────────────────────────────────────────────


class TestMlstParser:
    def test_demo(self, mod):
        m = mod.parse_mlst(DEMO_MLST)
        assert m["scheme"] == "ecoli_achtman_4"
        assert m["st"] == "131"
        assert m["call"] == "assigned"
        assert [a["locus"] for a in m["alleles"]] == ["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]
        assert all(a["flag"] == "exact" for a in m["alleles"])

    @pytest.mark.parametrize("token,allele,flag", [
        ("adk(53)", "53", "exact"),
        ("adk(~53)", "53", "novel"),
        ("adk(53?)", "53", "partial"),
        ("adk(-)", None, "missing"),
        ("adk(53,54)", "53,54", "multiple"),
    ])
    def test_allele_notation(self, mod, token, allele, flag):
        parsed = mod.parse_mlst_allele(token)
        assert parsed == {"locus": "adk", "allele": allele, "flag": flag}

    def test_novel_allele_blocks_st(self, mod, tmp_path):
        m = mod.parse_mlst(_write(tmp_path / "m.tsv", "a.fa\tecoli_achtman_4\t-\tadk(~53)\tfumC(40)\n"))
        assert m["st"] is None
        assert m["call"] == "unassigned"
        assert "adk" in m["note"]

    def test_novel_combination(self, mod, tmp_path):
        m = mod.parse_mlst(_write(tmp_path / "m.tsv", "a.fa\tsaureus\t-\tarcC(1)\taroE(2)\n"))
        assert m["call"] == "unassigned"
        assert "combination" in m["note"]

    def test_no_scheme(self, mod, tmp_path):
        m = mod.parse_mlst(_write(tmp_path / "m.tsv", "a.fa\t-\t-\n"))
        assert m["scheme"] is None and m["st"] is None
        assert m["call"] == "no_scheme"

    def test_header_line_is_skipped(self, mod, tmp_path):
        m = mod.parse_mlst(_write(tmp_path / "m.tsv", "FILE\tSCHEME\tST\tadk\na.fa\tecoli\t10\tadk(10)\n"))
        assert m["st"] == "10"

    def test_multiple_isolates_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="one isolate"):
            mod.parse_mlst(_write(tmp_path / "m.tsv", "a.fa\tecoli\t10\tadk(10)\nb.fa\tecoli\t11\tadk(11)\n"))

    def test_other_tools_table_is_named_as_wrong_format(self, mod):
        # Passing the AMRFinderPlus table as --mlst must say so, not "several assemblies".
        with pytest.raises(ValueError, match="not mlst output"):
            mod.parse_mlst(DEMO_AMR)

    def test_empty_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="empty"):
            mod.parse_mlst(_write(tmp_path / "m.tsv", "\n"))


# ── Organism inference ─────────────────────────────────────────────────────────


class TestOrganism:
    @pytest.mark.parametrize("scheme,organism", [
        ("ecoli_achtman_4", "Escherichia"),
        ("ecoli", "Escherichia"),
        ("senterica_achtman_2", "Salmonella"),
        ("klebsiella", "Klebsiella_pneumoniae"),
        ("saureus", "Staphylococcus_aureus"),
        ("paeruginosa", "Pseudomonas_aeruginosa"),
        ("abaumannii_2", "Acinetobacter_baumannii"),
    ])
    def test_known_schemes(self, mod, scheme, organism):
        assert mod.infer_organism(scheme) == organism

    @pytest.mark.parametrize("scheme", ["neisseria", "some_future_scheme", None, ""])
    def test_ambiguous_or_unknown_is_none(self, mod, scheme):
        # A wrong --organism silently applies the wrong point-mutation panel; never guess.
        assert mod.infer_organism(scheme) is None

    def test_organism_listing_ignores_the_surrounding_log(self, mod):
        # Real `amrfinder --list_organisms` output: other log lines must not become organisms.
        text = ("Running: amrfinder -l\nSoftware directory: '/env/bin/'\nSoftware version: 4.2.7\n"
                "Database directory: '/env/share/amrfinderplus/data/2026-08-07.1'\n"
                "Database version: 2026-08-07.1\n\n"
                "Available --organism options: Acinetobacter_baumannii, Escherichia, Vibrio_cholerae\n")
        assert mod.parse_organism_listing(text) == {"Acinetobacter_baumannii", "Escherichia", "Vibrio_cholerae"}
        assert mod.parse_organism_listing("Software version: 4.2.7\n") == set()

    def test_parse_organism_listing(self, mod):
        text = "Available --organism options: Acinetobacter_baumannii, Escherichia, Salmonella\n"
        assert mod.parse_organism_listing(text) == {"Acinetobacter_baumannii", "Escherichia", "Salmonella"}


# ── PlasmidFinder ──────────────────────────────────────────────────────────────


class TestPlasmidParser:
    def test_demo_results_tab(self, mod):
        reps = mod.parse_plasmidfinder(DEMO_PLASMID)
        assert [r["replicon"] for r in reps] == ["IncFII", "IncFIA", "Col156"]
        assert reps[0]["contig"] == "contig_2"  # description after whitespace dropped
        assert reps[0]["start"] == 3800 and reps[0]["stop"] == 4060
        assert reps[0]["identity"] == 100.0
        assert reps[0]["coverage"] == 100.0

    def test_coverage_from_query_template_length(self, mod, tmp_path):
        text = ("Database\tPlasmid\tIdentity\tQuery / Template length\tContig\tPosition in contig\t"
                "Note\tAccession number\nenterobacteriales\tIncI1\t97.5\t120 / 160\tc9\t5..124\t\tX\n")
        rep = mod.parse_plasmidfinder(_write(tmp_path / "p.tsv", text))[0]
        assert rep["coverage"] == 75.0

    def test_abricate_layout(self, mod, tmp_path):
        text = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t"
                "%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
                "a.fa\tc7\t100\t360\t+\tIncFII_1\t1-261/261\t===============\t0/0\t100.00\t99.62\t"
                "plasmidfinder\tAY458016\tIncFII\t\n")
        rep = mod.parse_plasmidfinder(_write(tmp_path / "a.tsv", text))[0]
        assert rep == {"replicon": "IncFII_1", "contig": "c7", "start": 100, "stop": 360,
                       "identity": 99.62, "coverage": 100.0, "accession": "AY458016"}

    def test_header_only_means_zero_replicons(self, mod, tmp_path):
        text = ("Database\tPlasmid\tIdentity\tQuery / Template length\tContig\tPosition in contig\t"
                "Note\tAccession number\n")
        assert mod.parse_plasmidfinder(_write(tmp_path / "p.tsv", text)) == []

    def test_unrecognised_header_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="PlasmidFinder"):
            mod.parse_plasmidfinder(_write(tmp_path / "p.tsv", "a\tb\n1\t2\n"))


# ── Integration of the three results ───────────────────────────────────────────


@pytest.fixture(scope="module")
def demo_result(mod):
    return mod.summarise(
        sample="demo",
        assembly=mod.read_assembly(DEMO_FASTA),
        amr_hits=mod.parse_amrfinder(DEMO_AMR),
        mlst=mod.parse_mlst(DEMO_MLST),
        replicons=mod.parse_plasmidfinder(DEMO_PLASMID),
        organism=None,
        mode="demo",
    )


class TestSummarise:
    def test_identity(self, demo_result):
        assert demo_result["skill"] == "isolate-amr-typing"
        assert demo_result["mlst"]["st"] == "131"
        assert "not a medical device" in demo_result["disclaimer"]

    def test_amr_and_point_split(self, demo_result):
        amr = demo_result["amr"]
        assert len(amr["genes"]) == 10
        assert len(amr["point_mutations"]) == 5
        assert {p["symbol"] for p in amr["point_mutations"]} == {
            "gyrA_S83L", "gyrA_D87N", "parC_S80I", "parC_E84V", "parE_I529L"}

    def test_non_amr_elements_kept_apart(self, demo_result):
        amr = demo_result["amr"]
        assert [e["symbol"] for e in amr["other_elements"]] == ["qacEdelta1"]
        assert "qacEdelta1" not in {g["symbol"] for g in amr["genes"]}

    def test_partial_hit_is_review_tier(self, demo_result):
        cat = next(g for g in demo_result["amr"]["genes"] if g["symbol"] == "catB3")
        assert cat["tier"] == "review"

    def test_class_summary_splits_multi_class(self, demo_result):
        classes = {c["drug_class"]: c for c in demo_result["amr"]["class_summary"]}
        # aac(6')-Ib-cr5 is AMINOGLYCOSIDE/QUINOLONE and must count under both.
        assert "aac(6')-Ib-cr5" in classes["AMINOGLYCOSIDE"]["confident"]
        assert "aac(6')-Ib-cr5" in classes["QUINOLONE"]["confident"]
        assert "gyrA_S83L" in classes["QUINOLONE"]["confident"]
        assert classes["PHENICOL"]["confident"] == []
        assert classes["PHENICOL"]["review"] == ["catB3"]
        assert "QUATERNARY AMMONIUM" not in classes

    def test_colocation_with_replicons(self, demo_result):
        genes = {g["symbol"]: g for g in demo_result["amr"]["genes"]}
        assert genes["blaCTX-M-15"]["replicons_on_contig"] == ["IncFII"]
        assert genes["sul1"]["replicons_on_contig"] == ["IncFIA"]
        assert genes["blaEC"]["replicons_on_contig"] == []

    def test_supplied_mlst_never_implies_an_organism(self, demo_result):
        # Only an mlst run made by the skill (stderr seen, no tie) may drive the organism.
        assert demo_result["amr"]["organism"] is None
        assert demo_result["amr"]["organism_source"] is None
        assert demo_result["amr"]["point_mutation_screening"] == "detected"

    def test_disrupted_gene_excluded_from_class_summary(self, mod, tmp_path):
        text = AMR_V4_HEADER + "\n" + _amr_row("blaTEM-1", "INTERNAL_STOP", subclass="BETA-LACTAM") + "\n"
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", text)),
                            mlst=None, replicons=None, organism=None, mode="precomputed")
        assert res["amr"]["genes"][0]["tier"] == "disrupted"
        beta = next(c for c in res["amr"]["class_summary"] if c["drug_class"] == "BETA-LACTAM")
        assert beta["confident"] == [] and beta["review"] == [] and beta["disrupted"] == ["blaTEM-1"]

    def test_point_screening_unknown_without_organism(self, mod, tmp_path):
        text = AMR_V4_HEADER + "\n" + _amr_row("sul2", "EXACTX", cls="SULFONAMIDE", subclass="SULFONAMIDE") + "\n"
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", text)),
                            mlst=None, replicons=None, organism=None, mode="precomputed")
        assert res["amr"]["point_mutation_screening"] == "unknown"
        assert any("point mutation" in w.lower() for w in res["warnings"])

    def test_own_run_without_organism_is_not_screened(self, mod, tmp_path):
        # We ran AMRFinderPlus ourselves with no organism, so we know nothing was screened.
        text = AMR_V4_HEADER + "\n" + _amr_row("sul2", "EXACTX", cls="SULFONAMIDE", subclass="SULFONAMIDE") + "\n"
        mlst = mod.parse_mlst(_write(tmp_path / "m.tsv", "a.fa\tneisseria\t1\tabcZ(1)\n"))
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", text)),
                            mlst=mlst, replicons=None, organism=None, mode="live", amr_run_here=True,
                            mlst_run_here=True)
        assert res["amr"]["point_mutation_screening"] == "not_screened"
        assert res["amr"]["organism"] is None
        assert res["amr"]["organism_applied"] is False
        assert any("Neisseria_gonorrhoeae" in w for w in res["warnings"])

    def test_applied_organism_is_recorded(self, mod, tmp_path):
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", AMR_V4_HEADER + "\n")),
                            mlst=None, replicons=None, organism=None, mode="live",
                            applied_organism="Salmonella", amr_run_here=True)
        assert res["amr"]["organism_applied"] is True
        assert res["amr"]["point_mutation_screening"] == "screened"

    def test_inferred_organism_is_not_claimed_as_applied(self, demo_result):
        assert demo_result["amr"]["organism_applied"] is False

    def test_user_organism_marks_screened(self, mod, tmp_path):
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", AMR_V4_HEADER + "\n")),
                            mlst=None, replicons=None, organism="Salmonella", mode="precomputed")
        assert res["amr"]["organism_source"] == "user"
        assert res["amr"]["point_mutation_screening"] == "screened"

    def test_missing_components_are_not_run_not_empty(self, mod):
        res = mod.summarise(sample="s", assembly=None, amr_hits=None,
                            mlst=mod.parse_mlst(DEMO_MLST), replicons=None, organism=None,
                            mode="precomputed")
        assert res["amr"]["status"] == "not_run"
        assert res["plasmids"]["status"] == "not_run"
        assert res["mlst"]["status"] == "ok"

    def test_contig_mismatch_warns(self, mod, tmp_path):
        text = AMR_V4_HEADER + "\n" + _amr_row("sul2", "EXACTX", contig="NODE_99") + "\n"
        res = mod.summarise(sample="s", assembly=mod.read_assembly(DEMO_FASTA),
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", text)),
                            mlst=None, replicons=None, organism=None, mode="precomputed")
        assert any("NODE_99" in w for w in res["warnings"])


# ── Tool invocation (never executed here) ──────────────────────────────────────


class TestCommands:
    def test_amrfinder_command(self, mod):
        cmd = mod.build_amrfinder_cmd(Path("a.fa"), Path("o.tsv"), organism="Escherichia", threads=4)
        assert cmd[0] == "amrfinder"
        assert cmd[cmd.index("--nucleotide") + 1] == "a.fa"
        assert cmd[cmd.index("--output") + 1] == "o.tsv"
        assert cmd[cmd.index("--organism") + 1] == "Escherichia"
        assert cmd[cmd.index("--threads") + 1] == "4"
        assert "--plus" in cmd

    def test_amrfinder_command_without_organism(self, mod):
        assert "--organism" not in mod.build_amrfinder_cmd(Path("a.fa"), Path("o.tsv"), organism=None, threads=1)

    def test_mlst_command(self, mod):
        assert mod.build_mlst_cmd(Path("a.fa")) == ["mlst", "a.fa"]

    def test_plasmidfinder2_command_pins_thresholds(self, mod):
        cmd = mod.build_plasmid_cmd("plasmidfinder2", ["plasmidfinder.py"], Path("a.fa"), Path("pf"), db=Path("/db"))
        assert cmd[0] == "plasmidfinder.py"
        assert cmd[cmd.index("-i") + 1] == "a.fa"
        assert cmd[cmd.index("-o") + 1] == "pf"
        assert cmd[cmd.index("-p") + 1] == "/db"
        assert cmd[cmd.index("-l") + 1] == "0.60"
        assert cmd[cmd.index("-t") + 1] == "0.95"
        assert "-x" in cmd

    def test_plasmidfinder3_command_writes_json(self, mod):
        # PlasmidFinder 3 ships as a module, has no plasmidfinder.py and writes no results_tab.tsv.
        cmd = mod.build_plasmid_cmd("plasmidfinder3", ["/env/python", "-m", "plasmidfinder"],
                                    Path("a.fa"), Path("pf"), db=None)
        assert cmd[:3] == ["/env/python", "-m", "plasmidfinder"]
        assert cmd[cmd.index("-j") + 1] == str(Path("pf") / "plasmidfinder.json")
        assert cmd[cmd.index("-l") + 1] == "0.60"
        assert cmd[cmd.index("-t") + 1] == "0.95"
        assert "-p" not in cmd

    def test_abricate_command_matches_plasmidfinder_thresholds(self, mod):
        # abricate defaults to 80/80; without these flags it reports hits PlasmidFinder would not.
        cmd = mod.build_plasmid_cmd("abricate", ["abricate"], Path("a.fa"), Path("pf"), db=None)
        assert cmd[:3] == ["abricate", "--db", "plasmidfinder"]
        assert cmd[cmd.index("--minid") + 1] == "95"
        assert cmd[cmd.index("--mincov") + 1] == "60"
        assert cmd[-1] == "a.fa"

    @pytest.mark.parametrize("on_path,module,expected", [
        ({"plasmidfinder.py", "python3", "abricate"}, True, "plasmidfinder2"),
        ({"python3", "abricate"}, True, "plasmidfinder3"),
        ({"python3", "abricate"}, False, "abricate"),
        ({"python3"}, False, None),
    ])
    def test_plasmid_backend_preference(self, mod, monkeypatch, on_path, module, expected):
        monkeypatch.setattr(mod.shutil, "which", lambda n: f"/bin/{n}" if n in on_path else None)
        monkeypatch.setattr(mod, "_module_available", lambda python, name: module)
        backend = mod.find_plasmid_backend()
        assert (backend[0] if backend else None) == expected

    def test_missing_tools_reported_together(self, mod, monkeypatch):
        monkeypatch.setattr(mod.shutil, "which", lambda name: None)
        with pytest.raises(mod.MissingToolError) as exc:
            mod.require_tools(["amrfinder", "mlst"])
        assert "amrfinder" in str(exc.value) and "mlst" in str(exc.value)


# ── CLI ────────────────────────────────────────────────────────────────────────


class TestCLI:
    def test_no_args_exits_nonzero(self):
        assert _run().returncode != 0

    def test_demo_writes_core_outputs(self, tmp_path):
        r = _run("--demo", "--output", tmp_path)
        assert r.returncode == 0, r.stderr
        for rel in ("report.md", "result.json", "tables/amr_determinants.csv",
                    "tables/mlst.csv", "tables/plasmid_replicons.csv",
                    "tables/drug_class_summary.csv", "figures/drug_class_determinants.png",
                    "reproducibility/commands.sh", "reproducibility/environment.yml",
                    "reproducibility/checksums.sha256"):
            assert (tmp_path / rel).exists(), rel

    def test_demo_report_content(self, tmp_path):
        _run("--demo", "--output", tmp_path)
        report = (tmp_path / "report.md").read_text()
        assert "ST131" in report
        assert "blaCTX-M-15" in report
        assert "gyrA_S83L" in report
        assert "IncFII" in report
        assert "synthetic" in report.lower()
        assert "not a medical device" in report.lower()
        assert "phenotyp" in report.lower()  # genotype is not a susceptibility result

    def test_demo_is_deterministic(self, tmp_path):
        _run("--demo", "--output", tmp_path / "a")
        _run("--demo", "--output", tmp_path / "b")
        assert (tmp_path / "a/result.json").read_text() == (tmp_path / "b/result.json").read_text()

    def test_checksums_verify_from_output_dir(self, tmp_path):
        _run("--demo", "--output", tmp_path)
        lines = (tmp_path / "reproducibility/checksums.sha256").read_text().splitlines()
        assert lines
        import hashlib
        for line in lines:
            digest, label = line.split("  ", 1)
            assert hashlib.sha256((tmp_path / label).read_bytes()).hexdigest() == digest

    def test_precomputed_mode_without_assembly(self, tmp_path):
        r = _run("--amrfinder", DEMO_AMR, "--mlst", DEMO_MLST, "--output", tmp_path)
        assert r.returncode == 0, r.stderr
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["mode"] == "precomputed"
        assert result["plasmids"]["status"] == "not_run"
        assert result["assembly"] is None
        assert "NOT RUN" in (tmp_path / "report.md").read_text()

    def test_assembly_plus_all_precomputed_needs_no_tools(self, tmp_path):
        r = _run("--input", DEMO_FASTA, "--amrfinder", DEMO_AMR, "--mlst", DEMO_MLST,
                 "--plasmidfinder", DEMO_PLASMID, "--output", tmp_path)
        assert r.returncode == 0, r.stderr
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["assembly"]["n_contigs"] == 4
        assert result["sample"] == "demo_isolate"

    def test_missing_input_file(self, tmp_path):
        r = _run("--input", tmp_path / "nope.fasta", "--output", tmp_path / "o")
        assert r.returncode != 0
        assert "not found" in r.stderr.lower()

    def test_invalid_assembly_fails_without_traceback(self, tmp_path):
        bad = _write(tmp_path / "p.faa", ">p1\nMKLVFWQEDHRPSTYNI\n")
        r = _run("--input", bad, "--amrfinder", DEMO_AMR, "--mlst", DEMO_MLST,
                 "--plasmidfinder", DEMO_PLASMID, "--output", tmp_path / "o")
        assert r.returncode != 0
        assert "Traceback" not in r.stderr
        assert not (tmp_path / "o" / "report.md").exists()

    def test_live_run_without_tools_names_them(self, tmp_path, monkeypatch):
        mod = _load()
        monkeypatch.setattr(mod.shutil, "which", lambda name: None)
        rc = mod.main(["--input", str(DEMO_FASTA), "--output", str(tmp_path / "o")])
        assert rc != 0
        assert not (tmp_path / "o" / "report.md").exists()

    def test_existing_output_warns(self, tmp_path):
        _run("--demo", "--output", tmp_path)
        r = _run("--demo", "--output", tmp_path)
        assert r.returncode == 0
        assert "already exists" in r.stderr


# ── mlst scheme ties (seen with mlst 2.35.0 on real genomes) ───────────────────

TIE_ECOLI = ("WARNING: salmonella(3529)==ecoli_achtman_4(131) score=100 genomes/jj1886.fna\n"
             "Remember that --minscore is only used when using automatic scheme detection.\n")
TIE_KLEB = "WARNING: ecoli_achtman_4(14464)==klebsiella(11) score=100 genomes/hs11286.fna\n"


class TestMlstSchemeTies:
    def test_parse_ties(self, mod):
        assert mod.parse_mlst_ties(TIE_ECOLI) == [("salmonella", "3529"), ("ecoli_achtman_4", "131")]
        assert mod.parse_mlst_ties(TIE_KLEB) == [("ecoli_achtman_4", "14464"), ("klebsiella", "11")]
        assert mod.parse_mlst_ties("This is mlst 2.35.0\nDone.\n") == []

    def test_new_salmonella_scheme_name_is_mapped(self, mod):
        assert mod.infer_organism("salmonella") == "Salmonella"

    def test_forced_scheme_command(self, mod):
        assert mod.build_mlst_cmd(Path("a.fa"), scheme="klebsiella") == ["mlst", "--scheme", "klebsiella", "a.fa"]

    def test_tie_resolved_only_by_a_stated_organism(self, mod):
        ties = mod.parse_mlst_ties(TIE_KLEB)
        assert mod.resolve_tie(ties, "Klebsiella_pneumoniae") == "klebsiella"
        assert mod.resolve_tie(ties, "Escherichia") == "ecoli_achtman_4"
        assert mod.resolve_tie(ties, None) is None
        assert mod.resolve_tie(ties, "Staphylococcus_aureus") is None

    def test_tied_scheme_never_drives_the_organism(self, mod, tmp_path):
        # mlst picked ecoli_achtman_4 for a Klebsiella genome; inferring Escherichia from
        # that would apply the wrong point-mutation panel.
        mlst = mod.parse_mlst(_write(tmp_path / "m.tsv", "hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n"))
        res = mod.summarise(sample="s", assembly=None,
                            amr_hits=mod.parse_amrfinder(_write(tmp_path / "a.tsv", AMR_V4_HEADER + "\n")),
                            mlst=mlst, replicons=None, organism=None, mode="live", amr_run_here=True,
                            mlst_ties=mod.parse_mlst_ties(TIE_KLEB))
        assert res["mlst"]["call"] == "ambiguous"
        assert res["mlst"]["tied_schemes"] == [{"scheme": "ecoli_achtman_4", "st": "14464"},
                                               {"scheme": "klebsiella", "st": "11"}]
        assert res["amr"]["organism"] is None
        assert res["amr"]["point_mutation_screening"] == "not_screened"
        assert any("klebsiella" in w and "ecoli_achtman_4" in w for w in res["warnings"])

    def test_ambiguous_st_is_not_headlined(self, mod, tmp_path):
        mlst = mod.parse_mlst(_write(tmp_path / "m.tsv", "hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n"))
        res = mod.summarise(sample="s", assembly=None, amr_hits=None, mlst=mlst, replicons=None,
                            organism=None, mode="live", mlst_ties=mod.parse_mlst_ties(TIE_KLEB))
        summary = mod.build_report(res).split("## MLST")[0]
        assert "AMBIGUOUS" in summary
        assert "ST14464 (scheme" not in summary

    def test_run_tools_reruns_mlst_with_the_scheme_for_the_stated_organism(self, mod, tmp_path, monkeypatch):
        calls = []

        class Proc:
            def __init__(self, stderr=""):
                self.stderr, self.stdout = stderr, ""

        def fake_run(cmd, stdout_path=None):
            calls.append(cmd)
            if "--scheme" in cmd:
                stdout_path.write_text("hs.fna\tklebsiella\t11\tgapA(3)\n")
                return Proc()
            stdout_path.write_text("hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n")
            return Proc(TIE_KLEB)

        monkeypatch.setattr(mod, "_run_cmd", fake_run)
        monkeypatch.setattr(mod, "_version", lambda cmd: "x")
        paths = {}
        info = mod.run_tools(DEMO_FASTA, tmp_path / "raw", need={"mlst"}, organism="Klebsiella_pneumoniae",
                             threads=1, plasmidfinder_db=None, paths=paths, mlst_scheme=None)
        assert calls[-1] == ["mlst", "--scheme", "klebsiella", str(DEMO_FASTA)]
        assert mod.parse_mlst(paths["mlst"])["st"] == "11"
        assert info["mlst_ties"] == []
        assert info["mlst_scheme_forced"] == "klebsiella"

    def test_run_tools_keeps_tie_when_nothing_resolves_it(self, mod, tmp_path, monkeypatch):
        class Proc:
            stderr, stdout = TIE_KLEB, ""

        def fake_run(cmd, stdout_path=None):
            stdout_path.write_text("hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n")
            return Proc()

        monkeypatch.setattr(mod, "_run_cmd", fake_run)
        monkeypatch.setattr(mod, "_version", lambda cmd: "x")
        info = mod.run_tools(DEMO_FASTA, tmp_path / "raw", need={"mlst"}, organism=None, threads=1,
                             plasmidfinder_db=None, paths={}, mlst_scheme=None)
        assert info["mlst_ties"] == [("ecoli_achtman_4", "14464"), ("klebsiella", "11")]
        assert info["applied_organism"] is None


# ── Preflight checks added after live runs ─────────────────────────────────────


class TestPreflight:
    def test_unknown_organism_is_rejected_with_the_valid_names(self, mod):
        with pytest.raises(ValueError, match="Escherichia"):
            mod.check_organism("Ecoli", {"Escherichia", "Salmonella"})

    def test_known_organism_passes(self, mod):
        mod.check_organism("Salmonella", {"Escherichia", "Salmonella"})

    def test_organism_not_checked_when_listing_unavailable(self, mod):
        mod.check_organism("Anything", set())

    def test_plasmidfinder3_needs_a_database_before_anything_runs(self, mod):
        with pytest.raises(ValueError, match="--plasmidfinder-db"):
            mod.check_plasmid_db("plasmidfinder3", None, {})
        mod.check_plasmid_db("plasmidfinder3", Path("/db"), {})
        mod.check_plasmid_db("plasmidfinder3", None, {"CGE_PLASMIDFINDER_DB": "/db"})
        mod.check_plasmid_db("plasmidfinder2", None, {})
        mod.check_plasmid_db("abricate", None, {})

    def test_decompressed_assembly_is_not_left_in_the_output(self, mod, tmp_path, monkeypatch):
        gz = tmp_path / "asm.fasta.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(DEMO_FASTA.read_text())
        seen = {}

        class Proc:
            stderr, stdout = "", ""

        def fake_run(cmd, stdout_path=None):
            seen["existed_during_run"] = Path(cmd[-1]).exists()
            stdout_path.write_text("asm.fasta\tsaureus\t36\tarcC(2)\n")
            return Proc()

        monkeypatch.setattr(mod, "_run_cmd", fake_run)
        monkeypatch.setattr(mod, "_version", lambda cmd: "mlst 2.35.0")
        info = mod.run_tools(gz, tmp_path / "raw", need={"mlst"}, organism=None, threads=1,
                             plasmidfinder_db=None, paths={}, mlst_scheme=None)
        assert seen["existed_during_run"]
        assert not (tmp_path / "raw" / "assembly.fasta").exists()
        assert info["versions"]["mlst"] == "2.35.0"


# ── Rerun hints for ties, replicon distance, split environments ────────────────


class TestTieRerunHints:
    def _tied(self, mod, tmp_path, stderr):
        mlst = mod.parse_mlst(_write(tmp_path / "m.tsv", "hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n"))
        return mod.summarise(sample="s", assembly=None, amr_hits=None, mlst=mlst, replicons=None,
                             organism=None, mode="live", mlst_run_here=True,
                             mlst_ties=mod.parse_mlst_ties(stderr))

    def test_each_tied_scheme_gets_exact_flags(self, mod, tmp_path):
        res = self._tied(mod, tmp_path, TIE_KLEB)
        assert res["mlst"]["rerun_options"] == [
            {"scheme": "ecoli_achtman_4", "st": "14464", "organism": "Escherichia",
             "flags": "--mlst-scheme ecoli_achtman_4 --organism Escherichia"},
            {"scheme": "klebsiella", "st": "11", "organism": "Klebsiella_pneumoniae",
             "flags": "--mlst-scheme klebsiella --organism Klebsiella_pneumoniae"},
        ]

    def test_scheme_without_an_organism_gets_scheme_flag_only(self, mod, tmp_path):
        res = self._tied(mod, tmp_path, "WARNING: bsubtilis(1)==neisseria(2) score=100 x.fna\n")
        assert [o["flags"] for o in res["mlst"]["rerun_options"]] == [
            "--mlst-scheme bsubtilis", "--mlst-scheme neisseria"]
        assert all(o["organism"] is None for o in res["mlst"]["rerun_options"])

    def test_report_shows_the_flags(self, mod, tmp_path):
        report = mod.build_report(self._tied(mod, tmp_path, TIE_KLEB))
        assert "`--mlst-scheme klebsiella --organism Klebsiella_pneumoniae`" in report
        assert "`--mlst-scheme ecoli_achtman_4 --organism Escherichia`" in report

    def test_no_hints_without_a_tie(self, demo_result):
        assert "rerun_options" not in demo_result["mlst"]

    def test_cli_prints_full_rerun_commands(self, tmp_path, monkeypatch, capsys):
        mod = _load()
        out = tmp_path / "o"

        def fake_run_tools(assembly, raw_dir, *, need, organism, threads, plasmidfinder_db, paths, mlst_scheme):
            raw_dir.mkdir(parents=True, exist_ok=True)
            paths["mlst"] = _write(raw_dir / "mlst.tsv", "hs.fna\tecoli_achtman_4\t14464\tadk(1769)\n")
            paths["amrfinder"] = _write(raw_dir / "amrfinder.tsv", AMR_V4_HEADER + "\n")
            paths["plasmidfinder"] = raw_dir / "plasmidfinder.tsv"
            paths["plasmidfinder"].write_text(DEMO_PLASMID.read_text())
            return {"applied_organism": None, "versions": {}, "mlst_scheme_forced": None,
                    "mlst_ties": mod.parse_mlst_ties(TIE_KLEB)}

        monkeypatch.setattr(mod, "run_tools", fake_run_tools)
        monkeypatch.setattr(mod, "require_tools", lambda names: None)
        monkeypatch.setattr(mod, "find_plasmid_backend", lambda: ("abricate", ["abricate"]))
        assert mod.main(["--input", str(DEMO_FASTA), "--output", str(out)]) == 0
        err = capsys.readouterr().err
        assert "mlst could not choose" in err
        line = next(l for l in err.splitlines() if "--mlst-scheme klebsiella" in l)
        assert "--organism Klebsiella_pneumoniae" in line
        assert f"--input {DEMO_FASTA}" in line
        # the plasmid result does not depend on the organism, so the rerun reuses it
        assert f"--plasmidfinder {out / 'tool_outputs' / 'plasmidfinder.tsv'}" in line


class TestReplicon_Distance:
    @pytest.mark.parametrize("a,b,expected", [
        ((150, 1022), (3800, 4060), 2778),
        ((3800, 4060), (150, 1022), 2778),
        ((100, 500), (400, 900), 0),
        ((100, 500), (501, 900), 1),
        ((None, 500), (501, 900), None),
    ])
    def test_interval_distance(self, mod, a, b, expected):
        assert mod.interval_distance(*a, *b) == expected

    def test_nearest_replicon_and_contig_length(self, demo_result):
        genes = {g["symbol"]: g for g in demo_result["amr"]["genes"]}
        assert genes["blaCTX-M-15"]["nearest_replicon"] == {"replicon": "IncFII", "distance_bp": 2778}
        assert genes["blaCTX-M-15"]["contig_length"] == 4200
        assert genes["blaEC"]["nearest_replicon"] is None
        assert genes["blaEC"]["contig_length"] == 9000

    def test_nearest_of_several(self, mod):
        res = mod.summarise(
            sample="JJ1886", assembly=None,
            amr_hits=mod.parse_amrfinder(FIXTURES / "jj1886_amrfinder_v4.tsv"),
            mlst=None, replicons=mod.parse_plasmidfinder(FIXTURES / "jj1886_plasmidfinder_v3.json"),
            organism="Escherichia", mode="precomputed")
        oxa = next(g for g in res["amr"]["genes"] if g["symbol"] == "blaOXA-1")
        near = oxa["nearest_replicon"]
        gaps = {r["replicon"]: mod.interval_distance(oxa["start"], oxa["stop"], r["start"], r["stop"])
                for r in res["plasmids"]["replicons"]}
        assert near["distance_bp"] == min(gaps.values())
        assert gaps[near["replicon"]] == near["distance_bp"]
        assert oxa["contig_length"] is None  # no assembly supplied

    def test_report_and_table_show_distance(self, tmp_path):
        _run("--demo", "--output", tmp_path)
        assert "IncFII (2.8 kb)" in (tmp_path / "report.md").read_text()
        header = (tmp_path / "tables/amr_determinants.csv").read_text().splitlines()[0]
        assert "nearest_replicon,nearest_replicon_distance_bp,contig_length" in header


class TestEnvironmentsAndOverrides:
    def test_two_environment_files(self, tmp_path):
        # AMRFinderPlus 4.2 and mlst 2.35 cannot be solved into one conda environment.
        _run("--demo", "--output", tmp_path)
        main_env = (tmp_path / "reproducibility/environment.yml").read_text()
        mlst_env = (tmp_path / "reproducibility/environment-mlst.yml").read_text()
        assert "ncbi-amrfinderplus" in main_env and "plasmidfinder" in main_env
        assert "mlst" not in main_env.replace("clawbio-isolate-amr-typing", "")
        assert "- mlst" in mlst_env and "amrfinder" not in mlst_env
        assert "name: clawbio-isolate-amr-typing-mlst" in mlst_env
        assert "bioconda" in mlst_env

    def test_default_tool_commands(self, mod, monkeypatch):
        for var in ("CLAWBIO_MLST_CMD", "CLAWBIO_AMRFINDER_CMD", "CLAWBIO_PLASMIDFINDER_CMD"):
            monkeypatch.delenv(var, raising=False)
        assert mod.tool_cmd("mlst") == ["mlst"]
        assert mod.tool_cmd("amrfinder") == ["amrfinder"]

    def test_override_commands_are_used(self, mod, monkeypatch):
        monkeypatch.setenv("CLAWBIO_MLST_CMD", "micromamba run -n typing mlst")
        monkeypatch.setenv("CLAWBIO_AMRFINDER_CMD", "/opt/amr/bin/amrfinder")
        assert mod.build_mlst_cmd(Path("a.fa"), "saureus") == [
            "micromamba", "run", "-n", "typing", "mlst", "--scheme", "saureus", "a.fa"]
        assert mod.build_amrfinder_cmd(Path("a.fa"), Path("o.tsv"), organism=None, threads=1)[0] == "/opt/amr/bin/amrfinder"

    @pytest.mark.parametrize("command,backend", [
        ("micromamba run -n pf plasmidfinder.py", "plasmidfinder2"),
        ("micromamba run -n pf python -m plasmidfinder", "plasmidfinder3"),
        ("micromamba run -n ab abricate", "abricate"),
    ])
    def test_plasmid_override_selects_backend(self, mod, monkeypatch, command, backend):
        monkeypatch.setenv("CLAWBIO_PLASMIDFINDER_CMD", command)
        monkeypatch.setattr(mod.shutil, "which", lambda n: None)
        assert mod.find_plasmid_backend() == (backend, command.split())

    def test_required_tool_check_follows_the_override(self, mod, monkeypatch):
        monkeypatch.setenv("CLAWBIO_MLST_CMD", "micromamba run -n typing mlst")
        monkeypatch.setattr(mod.shutil, "which", lambda n: "/bin/micromamba" if n == "micromamba" else None)
        mod.require_tools(["mlst"])
        with pytest.raises(mod.MissingToolError, match="amrfinder"):
            mod.require_tools(["mlst", "amrfinder"])

    def test_overrides_are_recorded_for_replay(self, mod):
        lines = mod.override_preflight({"CLAWBIO_MLST_CMD": "micromamba run -n typing mlst", "HOME": "/x"})
        assert lines == ["export CLAWBIO_MLST_CMD=\"${CLAWBIO_MLST_CMD:-micromamba run -n typing mlst}\""]
        assert mod.override_preflight({}) == []


# ── Real tool outputs (see tests/fixtures/README.md) ───────────────────────────

FIXTURES = Path(__file__).resolve().parent / "fixtures"


class TestRealToolOutputs:
    def test_amrfinder_v3(self, mod):
        hits = mod.parse_amrfinder(FIXTURES / "jj1886_amrfinder_v3.tsv")
        assert len(hits) == 33
        ctx = next(h for h in hits if h["symbol"] == "blaCTX-M-15")
        assert (ctx["call"], ctx["contig"], ctx["drug_class"]) == ("allele", "NC_022648.1", "BETA-LACTAM")

    def test_amrfinder_v4(self, mod):
        hits = mod.parse_amrfinder(FIXTURES / "jj1886_amrfinder_v4.tsv")
        assert len(hits) == 31
        assert all(h["symbol"] and h["element_type"] and h["method"] for h in hits)
        assert all(isinstance(h["start"], int) and isinstance(h["identity"], float) for h in hits)
        assert {h["call"] for h in hits} == {"allele", "exact", "blast", "partial", "point"}

    def test_mlst(self, mod):
        m = mod.parse_mlst(FIXTURES / "jj1886_mlst.tsv")
        assert (m["scheme"], m["st"], m["call"]) == ("ecoli_achtman_4", "131", "assigned")

    def test_plasmidfinder_v2_and_v3_agree(self, mod):
        v2 = mod.parse_plasmidfinder(FIXTURES / "jj1886_plasmidfinder_v2.tsv")
        v3 = mod.parse_plasmidfinder(FIXTURES / "jj1886_plasmidfinder_v3.json")
        assert v2 == v3
        assert v3[0] == {"replicon": "IncFIA", "contig": "NC_022651.1", "start": 6717, "stop": 7104,
                         "identity": 99.74, "coverage": 100.0, "accession": "AP001918"}
        assert [r["replicon"] for r in v3] == ["IncFIA", "IncFII"]

    def test_plasmidfinder_v3_json_with_no_hits(self, mod, tmp_path):
        empty = _write(tmp_path / "pf.json", json.dumps({"software_name": "PlasmidFinder", "seq_regions": {}}))
        assert mod.parse_plasmidfinder(empty) == []

    def test_other_json_rejected(self, mod, tmp_path):
        with pytest.raises(ValueError, match="PlasmidFinder"):
            mod.parse_plasmidfinder(_write(tmp_path / "x.json", json.dumps({"hello": 1})))

    def test_abricate(self, mod):
        reps = mod.parse_plasmidfinder(FIXTURES / "jj1886_abricate_plasmidfinder.tsv")
        assert len(reps) == 5
        assert reps[1]["replicon"] == "IncFIA_1" and reps[1]["contig"] == "NC_022651.1"

    def test_summary_of_real_isolate(self, mod):
        res = mod.summarise(
            sample="JJ1886", assembly=None,
            amr_hits=mod.parse_amrfinder(FIXTURES / "jj1886_amrfinder_v4.tsv"),
            mlst=mod.parse_mlst(FIXTURES / "jj1886_mlst.tsv"),
            replicons=mod.parse_plasmidfinder(FIXTURES / "jj1886_plasmidfinder_v3.json"),
            organism="Escherichia", mode="precomputed",
        )
        amr = res["amr"]
        genes = {g["symbol"]: g for g in amr["genes"]}
        # JJ1886 carries blaCTX-M-15 on the chromosome and blaOXA-1 on plasmid pJJ1886_5.
        assert genes["blaCTX-M-15"]["replicons_on_contig"] == []
        assert genes["blaOXA-1"]["replicons_on_contig"] == ["IncFIA", "IncFII"]
        assert genes["catB3"]["tier"] == "review"
        assert {"gyrA_S83L", "gyrA_D87N", "parC_S80I", "parC_E84V", "parE_I529L"} <= {
            p["symbol"] for p in amr["point_mutations"]}
        assert all(e["element_type"] in ("VIRULENCE", "STRESS") for e in amr["other_elements"])
        assert len(amr["genes"]) + len(amr["point_mutations"]) + len(amr["other_elements"]) == 31
        assert amr["point_mutation_screening"] == "detected"

    def test_cli_accepts_plasmidfinder_json(self, tmp_path):
        r = _run("--amrfinder", FIXTURES / "jj1886_amrfinder_v3.tsv", "--mlst", FIXTURES / "jj1886_mlst.tsv",
                 "--plasmidfinder", FIXTURES / "jj1886_plasmidfinder_v3.json", "--output", tmp_path)
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "tool_outputs" / "plasmidfinder.json").exists()
        assert "IncFIA" in (tmp_path / "report.md").read_text()


# ── SKILL.md output contract ───────────────────────────────────────────────────


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree."""
    if not skill_md.exists():
        return []
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", skill_md.read_text(), re.S)
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
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
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
    """Every artifact promised in SKILL.md '## Output Structure' must be produced."""

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md needs a parseable '## Output Structure' tree"
        assert any(p.startswith("reproducibility/") for p in promised)
        r = _run("--demo", "--output", tmp_path)
        assert r.returncode == 0, r.stderr
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, "SKILL.md promises outputs the skill did not write: " + ", ".join(missing)
