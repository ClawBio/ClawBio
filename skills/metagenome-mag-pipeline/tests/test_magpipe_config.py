"""Preset switches, upstream's config guards and resource capping.

The guards tested here are upstream's own, from the `Snakefile` and
`workflow/rules/helpers.smk` at the pinned commit. The wrapper repeats them so a
run fails in a second with a named flag rather than hours into a cohort, after
the tools are built and the reads are processed.
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

from magpipe_config import (
    RunOptions,
    build_config,
    validate_config,
    write_config,
)
from magpipe_errors import ErrorCode, SkillError
from magpipe_samplesheet import Sample
from magpipe_schemas import BINNERS, DIAMOND_SENSITIVITIES

BASE = yaml.safe_load((FIXTURES / "upstream_config.yaml").read_text())
BASE_TEST = yaml.safe_load((FIXTURES / "upstream_config_test.yaml").read_text())


def paired(sample: str = "S1") -> Sample:
    return Sample(
        sample=sample,
        fastq_1=(Path(f"/reads/{sample}_R1.fastq.gz"),),
        fastq_2=(Path(f"/reads/{sample}_R2.fastq.gz"),),
        group=sample,
        sra_run=None,
        paired=True,
    )


def single(sample: str = "S1") -> Sample:
    return Sample(
        sample=sample,
        fastq_1=(Path(f"/reads/{sample}.fastq.gz"),),
        fastq_2=(),
        group=sample,
        sra_run=None,
        paired=False,
    )


def remote(sample: str = "SRR1761673") -> Sample:
    return Sample(
        sample=sample, fastq_1=(), fastq_2=(), group=sample, sra_run=sample, paired=True
    )


@pytest.fixture
def databases(tmp_path: Path) -> dict[str, Path]:
    """Existing database paths of the right kind for every database flag."""
    checkm2 = tmp_path / "db" / "checkm2.dmnd"
    checkm2.parent.mkdir(parents=True, exist_ok=True)
    checkm2.write_text("dmnd\n")
    gtdbtk = tmp_path / "db" / "gtdbtk"
    gtdbtk.mkdir()
    (gtdbtk / "marker").write_text("x\n")
    others = {}
    for name in ("nr.dmnd", "kegg.dmnd", "cog.dmnd"):
        path = tmp_path / "db" / name
        path.write_text("dmnd\n")
        others[name] = path
    taxdump = tmp_path / "db" / "taxdump"
    taxdump.mkdir()
    (taxdump / "nodes.dmp").write_text("x\n")
    pfam = tmp_path / "db" / "Pfam-A.hmm"
    pfam.write_text("HMMER3/f\n")
    return {
        "checkm2_db": checkm2,
        "gtdbtk_db": gtdbtk,
        "taxonomy_db": others["nr.dmnd"],
        "taxdump": taxdump,
        "kegg_db": others["kegg.dmnd"],
        "cog_db": others["cog.dmnd"],
        "pfam_db": pfam,
    }


def full_options(databases: dict[str, Path], **kwargs) -> RunOptions:
    return RunOptions(preset="full", **databases, **kwargs)


def render(options: RunOptions, samples=None, base=None, tmp_path: Path | None = None):
    root = tmp_path or Path("/out")
    samples = [paired()] if samples is None else samples
    return build_config(
        copy.deepcopy(base or BASE),
        samples,
        options,
        samplesheet=root / "pipeline" / "samplesheet.tsv",
        results_dir=root / "results",
        database_dir=root / "db",
    )


class TestPresetSwitches:
    @pytest.mark.parametrize(
        ("preset", "profiler", "assembly", "binning", "heterogeneity"),
        [
            ("qc", "none", False, False, False),
            ("profile", "metaphlan4", False, False, False),
            ("assembly", "none", True, False, False),
            ("full", "metaphlan4", True, True, True),
        ],
    )
    def test_preset_switches(self, preset, profiler, assembly, binning, heterogeneity) -> None:
        config = render(RunOptions(preset=preset))
        assert config["profile"]["profiler"] == profiler
        assert config["assembly"]["enabled"] is assembly
        assert config["binning"]["enabled"] is binning
        assert config["heterogeneity"]["enabled"] is heterogeneity

    def test_qc_preset_turns_off_database_stages(self) -> None:
        config = render(RunOptions(preset="qc"))
        assert config["mag"]["enabled"] is False
        assert config["taxonomy"]["enabled"] is False
        assert config["annotation"]["enabled"] is False

    def test_unknown_preset_rejected(self) -> None:
        with pytest.raises(SkillError) as excinfo:
            render(RunOptions(preset="everything"))
        assert excinfo.value.error_code == ErrorCode.INVALID_PRESET

    def test_full_enables_mag_taxonomy_annotation_when_databases_given(
        self, databases: dict[str, Path], tmp_path: Path
    ) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        assert config["mag"]["enabled"] is True
        assert config["taxonomy"]["enabled"] is True
        assert config["annotation"]["enabled"] is True

    def test_full_disables_stages_when_no_databases_given(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="full"), tmp_path=tmp_path)
        assert config["binning"]["enabled"] is True
        assert config["mag"]["enabled"] is False
        assert config["taxonomy"]["enabled"] is False
        assert config["annotation"]["enabled"] is False


class TestPartialDatabases:
    def test_mag_needs_both_databases(self, databases: dict[str, Path], tmp_path: Path) -> None:
        options = full_options({**databases, "gtdbtk_db": None})
        with pytest.raises(SkillError) as excinfo:
            render(options, tmp_path=tmp_path)
        assert excinfo.value.error_code == ErrorCode.MISSING_DATABASE
        assert "--gtdbtk-db" in excinfo.value.fix

    def test_gtdbtk_error_names_the_unreachable_outputs(
        self, databases: dict[str, Path], tmp_path: Path
    ) -> None:
        options = full_options({**databases, "gtdbtk_db": None})
        with pytest.raises(SkillError) as excinfo:
            render(options, tmp_path=tmp_path)
        for path in (
            "final/binning/bins.tsv",
            "final/diagnostic/bottleneck.tsv",
            "final/benchmarks/runtime.tsv",
        ):
            assert path in excinfo.value.message, f"message does not name {path}"

    def test_taxonomy_needs_both_databases(
        self, databases: dict[str, Path], tmp_path: Path
    ) -> None:
        options = full_options({**databases, "taxdump": None})
        with pytest.raises(SkillError) as excinfo:
            render(options, tmp_path=tmp_path)
        assert excinfo.value.error_code == ErrorCode.MISSING_DATABASE
        assert "--taxdump" in excinfo.value.fix

    @pytest.mark.parametrize("missing", ["kegg_db", "cog_db", "pfam_db"])
    def test_annotation_needs_all_three(
        self, databases: dict[str, Path], tmp_path: Path, missing: str
    ) -> None:
        options = full_options({**databases, missing: None})
        with pytest.raises(SkillError) as excinfo:
            render(options, tmp_path=tmp_path)
        assert excinfo.value.error_code == ErrorCode.MISSING_DATABASE
        assert f"--{missing.replace('_db', '')}-db" in excinfo.value.fix

    def test_nothing_is_silently_disabled(self, databases: dict[str, Path], tmp_path: Path) -> None:
        """A partial set is an error, never a stage quietly turned off."""
        options = full_options({**databases, "kegg_db": None})
        with pytest.raises(SkillError) as excinfo:
            render(options, tmp_path=tmp_path)
        assert excinfo.value.details.get("stage") == "annotation"
        assert excinfo.value.details.get("missing") == ["--kegg-db"]


class TestDatabaseFlagsAgainstDisabledStages:
    """A database passed to a preset that disables its stage does nothing.

    The same honesty rule the demo follows: accepting the flag and leaving the
    stage off would have the user believe a database was in use when no rule read
    it. The stage is named, and the preset that would enable it.
    """

    @pytest.mark.parametrize(
        ("preset", "flag", "attribute", "stage"),
        [
            ("assembly", "--checkm2-db", "checkm2_db", "mag"),
            ("assembly", "--gtdbtk-db", "gtdbtk_db", "mag"),
            ("qc", "--taxonomy-db", "taxonomy_db", "taxonomy"),
            ("qc", "--taxdump", "taxdump", "taxonomy"),
            ("profile", "--kegg-db", "kegg_db", "annotation"),
            ("profile", "--cog-db", "cog_db", "annotation"),
            ("profile", "--pfam-db", "pfam_db", "annotation"),
        ],
    )
    def test_refused(
        self, tmp_path: Path, preset: str, flag: str, attribute: str, stage: str
    ) -> None:
        database = tmp_path / "db.dmnd"
        database.write_text("x\n")
        with pytest.raises(SkillError) as excinfo:
            render(RunOptions(preset=preset, **{attribute: database}), tmp_path=tmp_path)
        error = excinfo.value
        assert error.error_code == ErrorCode.INVALID_CONFIG
        assert flag in error.message
        # The config key is machine-readable, the message names the stage.
        assert error.details["stage"] == stage
        assert error.details["databases"]
        assert "--preset" in error.fix

    def test_assembly_preset_without_databases_is_unaffected(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="assembly"), tmp_path=tmp_path)
        assert config["mag"]["checkm2"]["db"] is None
        assert config["mag"]["enabled"] is False

    def test_full_still_accepts_the_databases(self, databases: dict[str, Path], tmp_path: Path) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        assert config["mag"]["enabled"] is True

    def test_qc_preset_with_no_database_flags_is_fine(self, tmp_path: Path) -> None:
        validate_config(render(RunOptions(preset="qc"), tmp_path=tmp_path), [paired()])


class TestConfigGuards:
    def test_binning_requires_assembly(self) -> None:
        config = render(RunOptions(preset="full"))
        config["assembly"]["enabled"] = False
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert excinfo.value.error_code == ErrorCode.INVALID_CONFIG
        assert "assembly" in excinfo.value.message

    def test_mag_requires_binning(self, databases: dict[str, Path], tmp_path: Path) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        config["binning"]["enabled"] = False
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert "binning" in excinfo.value.message

    def test_assembly_requires_paired_end(self) -> None:
        config = render(RunOptions(preset="assembly"), samples=[paired(), single("S2")])
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired(), single("S2")])
        assert excinfo.value.error_code == ErrorCode.ASSEMBLY_REQUIRES_PAIRED
        assert "S2" in excinfo.value.message

    def test_isolate_mode_rejected(self) -> None:
        config = render(RunOptions(preset="assembly"))
        config["assembly"]["spades"]["isolate"] = True
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert excinfo.value.error_code == ErrorCode.INVALID_CONFIG
        assert "isolate" in excinfo.value.message.lower()

    def test_unknown_diamond_sensitivity(self, databases: dict[str, Path], tmp_path: Path) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        config["taxonomy"]["sensitivity"] = "extremely-sensitive"
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert "mid-sensitive" in excinfo.value.message

    def test_every_diamond_sensitivity_is_accepted(
        self, databases: dict[str, Path], tmp_path: Path
    ) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        for sensitivity in DIAMOND_SENSITIVITIES:
            config["taxonomy"]["sensitivity"] = sensitivity
            validate_config(config, [paired()])

    def test_taxonomy_databases_required(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="full"))
        config["taxonomy"]["enabled"] = True
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert excinfo.value.error_code == ErrorCode.MISSING_DATABASE
        assert "taxonomy.db" in excinfo.value.message

    def test_annotation_databases_required(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="full"))
        config["annotation"]["enabled"] = True
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert "annotation.kegg_db" in excinfo.value.message

    def test_checkm2_db_must_be_a_file(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="full"))
        config["mag"]["enabled"] = True
        config["mag"]["checkm2"]["db"] = str(tmp_path)
        config["mag"]["gtdbtk"]["db"] = str(tmp_path)
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert excinfo.value.error_code == ErrorCode.MISSING_DATABASE
        assert "checkm2" in excinfo.value.message

    def test_gtdbtk_db_must_be_a_directory(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="full"))
        config["mag"]["enabled"] = True
        config["mag"]["checkm2"]["db"] = str(tmp_path / "x.dmnd")
        (tmp_path / "x.dmnd").write_text("dmnd\n")
        config["mag"]["gtdbtk"]["db"] = str(tmp_path / "x.dmnd")
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert "gtdbtk" in excinfo.value.message

    def test_unknown_binner(self) -> None:
        config = render(RunOptions(preset="full", binners=("metabat", "concoct")))
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert excinfo.value.error_code == ErrorCode.INVALID_CONFIG
        assert "concoct" in excinfo.value.message

    def test_consensus_needs_two_binners(self) -> None:
        config = render(RunOptions(preset="full", binners=("metabat",)))
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [paired()])
        assert "at least two" in excinfo.value.message

    def test_binner_allowlist(self) -> None:
        assert BINNERS == ("metabat", "maxbin", "vamb")
        for binner in BINNERS:
            config = render(RunOptions(preset="full", binners=BINNERS[:2] + (binner,)))
            if len(set(BINNERS[:2] + (binner,))) < 2:
                continue
            validate_config(config, [paired()])


class TestRemoteInputs:
    def test_fetch_enabled_only_with_the_flag(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="qc", allow_remote_inputs=True), samples=[remote()])
        assert config["fetch"]["sra"]["enabled"] is True

    def test_fetch_disabled_without_the_flag(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="qc"), samples=[remote()])
        assert config["fetch"]["sra"]["enabled"] is False

    def test_remote_sample_without_the_flag_is_refused(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="qc"), samples=[remote()])
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [remote()])
        assert excinfo.value.error_code == ErrorCode.REMOTE_INPUT_NOT_ALLOWED
        assert "--allow-remote-inputs" in excinfo.value.fix

    def test_remote_sample_with_the_flag_validates(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="qc", allow_remote_inputs=True), samples=[remote()])
        validate_config(config, [remote()])

    def test_local_samples_never_enable_fetching(self) -> None:
        config = render(RunOptions(preset="qc", allow_remote_inputs=True))
        assert config["fetch"]["sra"]["enabled"] is False


class TestPaths:
    def test_core_paths_are_absolute(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="assembly"))
        for key in ("samplesheet", "results_dir", "database_dir"):
            assert Path(config[key]).is_absolute(), key

    def test_host_fasta_is_absolute(self, tmp_path: Path) -> None:
        host = tmp_path / "host.fa"
        host.write_text(">h\nACGT\n")
        config = render(RunOptions(preset="qc", host_fasta=host))
        assert Path(config["qc"]["decontaminate"]["fasta"]).is_absolute()
        assert config["qc"]["decontaminate"]["enabled"] is True

    def test_adapter_fasta_stays_relative_to_the_checkout(self, tmp_path: Path) -> None:
        """qc.adapter.fasta lives in the checkout; absolutising it would break it."""
        config = render(RunOptions(preset="qc"))
        assert config["qc"]["adapter"]["fasta"] == "workflow/resources/illumina_adapters.fa"

    def test_databases_are_absolute(self, databases: dict[str, Path], tmp_path: Path) -> None:
        config = render(full_options(databases), tmp_path=tmp_path)
        for value in (
            config["mag"]["checkm2"]["db"],
            config["mag"]["gtdbtk"]["db"],
            config["taxonomy"]["db"],
            config["taxonomy"]["taxdump"],
            config["annotation"]["kegg_db"],
            config["annotation"]["cog_db"],
            config["annotation"]["pfam_db"],
        ):
            assert Path(value).is_absolute(), value

    def test_skip_host_removal(self) -> None:
        config = render(RunOptions(preset="qc", skip_host_removal=True))
        assert config["qc"]["decontaminate"]["enabled"] is False

    def test_coassemble_and_binners(self) -> None:
        config = render(RunOptions(preset="full", coassemble=True, binners=("metabat", "vamb")))
        assert config["assembly"]["coassemble"] is True
        assert config["binning"]["binners"] == ["metabat", "vamb"]


class TestResourceCapping:
    def test_threads_capped_at_cores(self) -> None:
        config = render(RunOptions(preset="full", cores=2))
        assert max(config["threads"].values()) <= 2
        assert config["threads"]["spades"] == 2

    def test_memory_capped_at_mem_mb(self) -> None:
        config = render(RunOptions(preset="full", mem_mb=4000))
        assert max(config["mem_mb"].values()) <= 4000

    def test_capping_never_raises_a_base_value(self) -> None:
        generous = render(RunOptions(preset="full", cores=64, mem_mb=999999))
        for key, value in BASE["threads"].items():
            assert generous["threads"][key] == value
        for key, value in BASE["mem_mb"].items():
            assert generous["mem_mb"][key] == value

    def test_resources_stay_positive(self) -> None:
        config = render(RunOptions(preset="full", cores=1, mem_mb=1))
        assert all(value >= 1 for value in config["threads"].values())
        assert all(value >= 1 for value in config["mem_mb"].values())


class TestImmutabilityAndIO:
    def test_base_config_not_mutated(self, tmp_path: Path) -> None:
        base = copy.deepcopy(BASE)
        snapshot = copy.deepcopy(base)
        build_config(
            base,
            [paired()],
            RunOptions(preset="full", cores=1),
            samplesheet=tmp_path / "s.tsv",
            results_dir=tmp_path / "r",
            database_dir=tmp_path / "db",
        )
        assert base == snapshot

    def test_demo_config_differs_only_in_the_three_path_keys(self, tmp_path: Path) -> None:
        config = build_config(
            copy.deepcopy(BASE_TEST),
            [paired()],
            RunOptions(preset="assembly", demo=True),
            samplesheet=tmp_path / "s.tsv",
            results_dir=tmp_path / "results",
            database_dir=tmp_path / "db",
        )
        differing = {key for key in set(config) | set(BASE_TEST) if config.get(key) != BASE_TEST.get(key)}
        assert differing == {"samplesheet", "results_dir", "database_dir"}
        assert config["assembly"]["enabled"] is True
        assert config["profile"]["profiler"] == "none"

    def test_demo_config_refuses_a_preset_it_cannot_apply(self, tmp_path: Path) -> None:
        """`--demo --preset full` is a contradiction, not a silent no-op."""
        with pytest.raises(SkillError) as excinfo:
            build_config(
                copy.deepcopy(BASE_TEST),
                [paired()],
                RunOptions(preset="full", demo=True),
                samplesheet=tmp_path / "s.tsv",
                results_dir=tmp_path / "results",
                database_dir=tmp_path / "db",
            )
        assert excinfo.value.error_code == ErrorCode.INVALID_PRESET
        assert "--demo" in excinfo.value.fix
        assert excinfo.value.details["demo_base_config"] == "workflow/config/config_test.yaml"

    def test_demo_config_keeps_the_test_databases_off(self, tmp_path: Path) -> None:
        """Database flags cannot turn a stage on in the demo: it has none enabled."""
        config = build_config(
            copy.deepcopy(BASE_TEST),
            [paired()],
            RunOptions(preset="assembly", demo=True, checkm2_db=tmp_path / "c.dmnd"),
            samplesheet=tmp_path / "s.tsv",
            results_dir=tmp_path / "results",
            database_dir=tmp_path / "db",
        )
        assert config["mag"]["enabled"] is False
        assert config["mag"]["checkm2"]["db"] is None

    def test_write_config_round_trips(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="assembly"))
        dest = write_config(config, tmp_path / "pipeline" / "config.yaml")
        reloaded = yaml.safe_load(dest.read_text())
        assert reloaded == config

    def test_written_yaml_keeps_upstream_keys(self, tmp_path: Path) -> None:
        config = render(RunOptions(preset="assembly"))
        dest = write_config(config, tmp_path / "config.yaml")
        reloaded = yaml.safe_load(dest.read_text())
        assert set(reloaded) == set(BASE)
        assert reloaded["qc"]["adapter"]["fasta"] == BASE["qc"]["adapter"]["fasta"]

    def test_validate_accepts_every_preset(self, tmp_path: Path) -> None:
        for preset in ("qc", "profile", "assembly"):
            config = render(RunOptions(preset=preset))
            validate_config(config, [paired()])

    def test_validate_rejects_an_empty_samplesheet(self) -> None:
        config = render(RunOptions(preset="qc"))
        with pytest.raises(SkillError) as excinfo:
            validate_config(config, [])
        assert excinfo.value.error_code == ErrorCode.INVALID_CONFIG
