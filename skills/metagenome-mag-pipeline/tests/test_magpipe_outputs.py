"""Parsing the pipeline's final tables, and deciding which tables must exist.

Two rules dominate this module:

* Column names are upstream's, read from the pinned commit's `workflow/scripts/`.
  A guess here would silently mislabel a column.
* An empty cell means the axis was not measured. It is never zero, because "we
  could not measure this" and "we measured zero" are different claims about the
  data, and upstream says so in its own docs.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_config import RunOptions, build_config
from magpipe_errors import ErrorCode, SkillError
from magpipe_outputs import expected_outputs, parse_outputs
from magpipe_samplesheet import Sample

BASE = yaml.safe_load((FIXTURES / "upstream_config.yaml").read_text())


def config_for(preset: str, tmp_path: Path, **kwargs) -> dict:
    return build_config(
        copy.deepcopy(BASE),
        [Sample("S1", (Path("/r/S1_R1.fastq.gz"),), (Path("/r/S1_R2.fastq.gz"),), "S1", None, True)],
        RunOptions(preset=preset, **kwargs),
        samplesheet=tmp_path / "pipeline" / "samplesheet.tsv",
        results_dir=tmp_path / "results",
        database_dir=tmp_path / "db",
    )


def full_config(fake_results: Path, tmp_path: Path) -> dict:
    db = fake_results.parent / "dbs"
    (db / "gtdbtk").mkdir(parents=True, exist_ok=True)
    (db / "taxdump").mkdir(parents=True, exist_ok=True)
    for name in ("checkm2.dmnd", "nr.dmnd", "kegg.dmnd", "cog.dmnd", "Pfam-A.hmm"):
        (db / name).write_text("x\n")
    return config_for(
        "full",
        tmp_path,
        checkm2_db=db / "checkm2.dmnd",
        gtdbtk_db=db / "gtdbtk",
        taxonomy_db=db / "nr.dmnd",
        taxdump=db / "taxdump",
        kegg_db=db / "kegg.dmnd",
        cog_db=db / "cog.dmnd",
        pfam_db=db / "Pfam-A.hmm",
    )


class TestExpectedOutputs:
    def test_always_present(self, tmp_path: Path) -> None:
        for preset in ("qc", "profile", "assembly", "full"):
            paths = expected_outputs(config_for(preset, tmp_path))
            for required in (
                "final/qc/read_counts.tsv",
                "final/diagnostic/bottleneck.tsv",
                "final/benchmarks/runtime.tsv",
                "provenance/software_versions.tsv",
                "provenance/pipeline_revision.tsv",
                "provenance/reference_db.md5",
                "provenance/run_manifest.json",
            ):
                assert required in paths, f"{preset} must require {required}"

    def test_qc_preset(self, tmp_path: Path) -> None:
        paths = expected_outputs(config_for("qc", tmp_path))
        assert "final/profile/profile.tsv" not in paths
        assert "final/assembly/quast.tsv" not in paths
        assert "final/binning/bins.tsv" not in paths
        assert "final/heterogeneity/strain_heterogeneity.tsv" not in paths

    def test_profile_preset(self, tmp_path: Path) -> None:
        paths = expected_outputs(config_for("profile", tmp_path))
        assert "final/profile/profile.tsv" in paths
        assert "final/assembly/quast.tsv" not in paths

    def test_assembly_preset(self, tmp_path: Path) -> None:
        paths = expected_outputs(config_for("assembly", tmp_path))
        assert "final/assembly/quast.tsv" in paths
        assert "final/profile/profile.tsv" not in paths
        assert "final/binning/bins.tsv" not in paths

    def test_full_preset(self, fake_results: Path, tmp_path: Path) -> None:
        paths = expected_outputs(full_config(fake_results, tmp_path))
        for required in (
            "final/profile/profile.tsv",
            "final/assembly/quast.tsv",
            "final/binning/bins.tsv",
            "final/binning/checkm2/quality_summary.tsv",
            "final/binning/drep/GenomeInfo.csv",
            "final/binning/gtdb/gtdb_summary.tsv",
            "final/taxonomy/contig_taxonomy.tsv",
            "final/annotation/function_abundance.tsv",
            "final/heterogeneity/strain_heterogeneity.tsv",
        ):
            assert required in paths

    def test_binning_tables_per_assembly_unit_are_a_glob(self, fake_results: Path, tmp_path: Path) -> None:
        paths = expected_outputs(full_config(fake_results, tmp_path))
        assert "intermediate/binning/dastool/*/contig2bin.tsv" in paths

    def test_multiqc_only_when_enabled(self, tmp_path: Path) -> None:
        config = config_for("qc", tmp_path)
        assert "reports/multiqc/multiqc_report.html" in expected_outputs(config)
        config["report"]["multiqc"] = False
        assert "reports/multiqc/multiqc_report.html" not in expected_outputs(config)

    def test_sourmash_only_when_enabled(self, tmp_path: Path) -> None:
        config = config_for("qc", tmp_path)
        config["qc"]["sourmash"]["enabled"] = False
        assert "final/qc/sourmash_similarities.csv" not in expected_outputs(config)

    def test_taxonomy_needs_assembly(self, fake_results: Path, tmp_path: Path) -> None:
        config = full_config(fake_results, tmp_path)
        config["assembly"]["enabled"] = False
        config["binning"]["enabled"] = False
        config["mag"]["enabled"] = False
        paths = expected_outputs(config)
        assert "final/taxonomy/contig_taxonomy.tsv" not in paths
        assert "final/annotation/function_abundance.tsv" not in paths

    def test_heterogeneity_needs_binning(self, fake_results: Path, tmp_path: Path) -> None:
        config = full_config(fake_results, tmp_path)
        config["binning"]["enabled"] = False
        config["mag"]["enabled"] = False
        assert "final/heterogeneity/strain_heterogeneity.tsv" not in expected_outputs(config)

    def test_minpath_outputs_only_with_a_data_dir(self, fake_results: Path, tmp_path: Path) -> None:
        config = full_config(fake_results, tmp_path)
        assert not any("minpath" in path for path in expected_outputs(config))
        config["annotation"]["minpath"]["data_dir"] = str(tmp_path / "minpath")
        assert any("minpath" in path for path in expected_outputs(config))

    def test_paths_are_relative_to_the_results_dir(self, tmp_path: Path) -> None:
        for path in expected_outputs(config_for("full", tmp_path)):
            assert not path.startswith("/")


class TestParseFull:
    @pytest.fixture
    def parsed(self, fake_results: Path, tmp_path: Path) -> dict:
        return parse_outputs(fake_results, full_config(fake_results, tmp_path))

    def test_top_level_keys(self, parsed: dict) -> None:
        assert set(parsed) == {
            "qc",
            "profile",
            "assembly",
            "bins",
            "diagnostic",
            "taxonomy",
            "annotation",
            "heterogeneity",
            "provenance",
        }

    def test_read_counts(self, parsed: dict) -> None:
        rows = parsed["qc"]["read_counts"]
        assert [row["sample"] for row in rows] == ["S1", "S2"]
        assert rows[0]["retained_reads"] == 4000
        assert rows[0]["retained_bases"] == 1200000
        assert rows[0]["retained_pairs"] == 2000
        assert rows[0]["layout"] == "paired"

    def test_blank_cell_is_none_not_zero(self, parsed: dict) -> None:
        rows = parsed["qc"]["read_counts"]
        assert rows[1]["retained_pairs"] is None
        assert rows[1]["retained_reads"] == 2500

    def test_sourmash_matrix(self, parsed: dict) -> None:
        matrix = parsed["qc"]["sourmash"]
        assert matrix["samples"] == ["S1", "S2"]
        assert matrix["values"]["S1"]["S1"] == 1.0
        assert matrix["values"]["S2"]["S1"] == 0.31

    def test_profile_is_relative_abundance(self, parsed: dict) -> None:
        profile = parsed["profile"]
        assert profile["version"] == "mpa_vJan26_CHOCOPhlAnSGB_202605"
        assert profile["samples"] == ["S1", "S2"]
        top = {row["clade_name"]: row for row in profile["clades"]}
        assert top["k__Bacteria"]["S1"] == 98.5
        top_taxa = {row["clade_name"]: row["relative_abundance_pct"] for row in profile["top_taxa"]["S2"]}
        assert top_taxa["k__Bacteria|p__Proteobacteria"] == 70.5
        assert profile["top_taxa"]["S2"][0]["clade_name"] == "k__Bacteria"

    def test_quast_is_transposed_per_assembly_unit(self, parsed: dict) -> None:
        quast = parsed["assembly"]["quast"]
        assert quast["S1"]["# contigs"] == 331
        assert quast["S1"]["GC (%)"] == 51.2
        assert parsed["assembly"]["units"] == ["S1"]

    def test_bins_summary(self, parsed: dict) -> None:
        bins = parsed["bins"]
        assert bins["count"] == 3
        assert bins["passing_qa"] == 2
        assert bins["genomes"] == ["S1.bin.1", "S1.bin.2", "S2.bin.1"]
        assert bins["checkm2"][0]["completeness"] == 94.53
        assert bins["checkm2"][0]["passes_qa"] is True
        assert bins["checkm2"][2]["passes_qa"] is False

    def test_bin_thresholds_come_from_the_config(self, parsed: dict, fake_results: Path, tmp_path: Path) -> None:
        config = full_config(fake_results, tmp_path)
        assert parsed["bins"]["thresholds"] == {
            "checkm2_min_completeness": config["mag"]["checkm2"]["min_completeness"],
            "checkm2_max_contamination": config["mag"]["checkm2"]["max_contamination"],
            "drep_completeness": config["mag"]["drep"]["completeness"],
            "drep_contamination": config["mag"]["drep"]["contamination"],
        }

    def test_bottleneck_blank_capture_is_none(self, parsed: dict) -> None:
        rows = parsed["diagnostic"]["bottleneck"]
        assert rows[0]["assembly_capture_pct"] == 77.5
        assert rows[1]["assembly_capture_pct"] is None
        assert rows[1]["catalogue_capture_pct"] == 60.0

    def test_taxonomy_rows(self, parsed: dict) -> None:
        rows = parsed["taxonomy"]["contigs"]
        assert rows[0]["taxonomy"] == "k__Bacteria;p__Bacillota"
        assert rows[0]["disparity"] == 0.05
        assert rows[1]["taxonomy"] == "unclassified"

    def test_annotation_rows(self, parsed: dict) -> None:
        rows = parsed["annotation"]["functions"]
        assert rows[0] == {"bin": "unassigned", "source": "kegg", "function": "ko:K00001", "genes": 31}

    def test_heterogeneity_blank_rate_is_none(self, parsed: dict) -> None:
        rows = parsed["heterogeneity"]["rows"]
        assert rows[0]["heterogeneity_pct"] == 0.1493
        assert rows[2]["heterogeneity_pct"] is None
        assert rows[2]["covered_positions"] == 0

    def test_provenance(self, parsed: dict) -> None:
        provenance = parsed["provenance"]
        assert provenance["pipeline"]["commit"] == "7351a702d29857801963cae7aeff6600d59ebe77"
        assert provenance["pipeline"]["dirty"] == "false"
        assert {row["tool"] for row in provenance["tools"]} == {
            "snakemake",
            "bbtools",
            "metaphlan",
        }
        assert provenance["manifest"]["samplesheet_sha256"] == "def"

    def test_reference_checksums_are_listed(self, parsed: dict) -> None:
        resources = {row["resource"] for row in parsed["provenance"]["references"]}
        assert "host_grch38.fa" in resources

    def test_runtime_is_available(self, parsed: dict) -> None:
        rows = parsed["provenance"]["runtime"]
        assert {row["rule"] for row in rows} == {"bbtools", "spades"}


class TestParseStagesOff:
    def test_qc_only_config_maps_stages_to_none(self, fake_results: Path, tmp_path: Path) -> None:
        config = config_for("qc", tmp_path)
        config["qc"]["sourmash"]["enabled"] = False
        for name in (
            "final/profile/profile.tsv",
            "final/assembly/quast.tsv",
            "final/binning/bins.tsv",
            "final/taxonomy/contig_taxonomy.tsv",
            "final/annotation/function_abundance.tsv",
            "final/heterogeneity/strain_heterogeneity.tsv",
            "final/qc/sourmash_similarities.csv",
        ):
            (fake_results / name).unlink(missing_ok=True)
        (fake_results / "provenance" / "run_manifest.json").write_text("{}\n")
        (fake_results / "provenance" / "reference_db.md5").write_text("")
        parsed = parse_outputs(fake_results, config)
        assert parsed["profile"] is None
        assert parsed["assembly"] is None
        assert parsed["bins"] is None
        assert parsed["taxonomy"] is None
        assert parsed["annotation"] is None
        assert parsed["heterogeneity"] is None
        assert parsed["qc"]["sourmash"] is None
        assert parsed["diagnostic"] is not None


class TestMissingOutputs:
    def test_missing_file_raises_listing_it(self, fake_results: Path, tmp_path: Path) -> None:
        (fake_results / "final" / "qc" / "read_counts.tsv").unlink()
        with pytest.raises(SkillError) as excinfo:
            parse_outputs(fake_results, full_config(fake_results, tmp_path))
        error = excinfo.value
        assert error.error_code == ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND
        assert "final/qc/read_counts.tsv" in error.details["missing"]

    def test_all_missing_files_are_listed_together(self, fake_results: Path, tmp_path: Path) -> None:
        (fake_results / "final" / "qc" / "read_counts.tsv").unlink()
        (fake_results / "provenance" / "software_versions.tsv").unlink()
        with pytest.raises(SkillError) as excinfo:
            parse_outputs(fake_results, full_config(fake_results, tmp_path))
        assert len(excinfo.value.details["missing"]) == 2

    def test_empty_glob_is_reported(self, fake_results: Path, tmp_path: Path) -> None:
        for path in (fake_results / "intermediate" / "binning" / "dastool").glob("*/contig2bin.tsv"):
            path.unlink()
        with pytest.raises(SkillError) as excinfo:
            parse_outputs(fake_results, full_config(fake_results, tmp_path))
        assert "contig2bin.tsv" in " ".join(excinfo.value.details["missing"])

    def test_results_dir_missing_entirely(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            parse_outputs(tmp_path / "nope", config_for("qc", tmp_path))
        assert excinfo.value.error_code == ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND


class TestNumberCoercion:
    def test_integers_floats_and_booleans(self, fake_results: Path, tmp_path: Path) -> None:
        parsed = parse_outputs(fake_results, full_config(fake_results, tmp_path))
        row = parsed["bins"]["checkm2"][0]
        assert isinstance(row["completeness"], float)
        assert isinstance(row["passes_qa"], bool)
        assert isinstance(parsed["qc"]["read_counts"][0]["retained_reads"], int)

    def test_identifiers_stay_strings(self, fake_results: Path, tmp_path: Path) -> None:
        parsed = parse_outputs(fake_results, full_config(fake_results, tmp_path))
        assert parsed["bins"]["genomes"][0] == "S1.bin.1"
        assert isinstance(parsed["bins"]["genomes"][0], str)

    def test_empty_strings_become_none(self, fake_results: Path, tmp_path: Path) -> None:
        parsed = parse_outputs(fake_results, full_config(fake_results, tmp_path))
        unclassified = next(row for row in parsed["bins"]["gtdb"] if row["user_genome"] == "S1.bin.2")
        assert unclassified["classification"] is None
        assert unclassified["closest_genome_reference"] is None
