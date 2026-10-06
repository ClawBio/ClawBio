"""Render the upstream config from a preset, and re-check upstream's own guards.

`build_config` never mutates the base config it is given: it renders a new dict
with the preset's stage switches, absolute paths and capped resources.
`validate_config` repeats every guard upstream's `Snakefile` and `helpers.smk`
apply, so an invalid combination fails here — in a second, naming the flag —
rather than after the tools are built and the reads are processed.
"""

from __future__ import annotations

import copy
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_schemas import BINNERS, DIAMOND_SENSITIVITIES, PRESETS

STAGE = "config"

# What each preset switches on. `full` enables the database-dependent stages only
# when their databases are supplied; see `PRESET_SWITCHES["full"]` handling in
# build_config for the partial-set rule.
PRESET_SWITCHES: dict[str, dict[str, Any]] = {
    "qc": {
        "profiler": "none",
        "assembly": False,
        "binning": False,
        "heterogeneity": False,
        "mag": False,
        "taxonomy": False,
        "annotation": False,
    },
    "profile": {
        "profiler": "metaphlan4",
        "assembly": False,
        "binning": False,
        "heterogeneity": False,
        "mag": False,
        "taxonomy": False,
        "annotation": False,
    },
    "assembly": {
        "profiler": "none",
        "assembly": True,
        "binning": False,
        "heterogeneity": False,
        "mag": False,
        "taxonomy": False,
        "annotation": False,
    },
    "full": {
        "profiler": "metaphlan4",
        "assembly": True,
        "binning": True,
        "heterogeneity": True,
        "mag": True,
        "taxonomy": True,
        "annotation": True,
    },
}

# stage -> (RunOptions attribute holding its database, CLI flag) pairs. Every one
# of them is required when the stage is on: upstream rejects a null at
# config-parse time, and a half-configured database is worse than none.
STAGE_DATABASES: dict[str, tuple[tuple[str, str], ...]] = {
    "mag": (("checkm2_db", "--checkm2-db"), ("gtdbtk_db", "--gtdbtk-db")),
    "taxonomy": (("taxonomy_db", "--taxonomy-db"), ("taxdump", "--taxdump")),
    "annotation": (
        ("kegg_db", "--kegg-db"),
        ("cog_db", "--cog-db"),
        ("pfam_db", "--pfam-db"),
    ),
}


# How each stage is named in an error message. The keys are the same stages
# STAGE_DATABASES covers.
_STAGE_LABELS = {
    "mag": "MAG assessment (CheckM2, dRep, GTDB-Tk)",
    "taxonomy": "contig taxonomy",
    "annotation": "functional annotation",
}


@dataclass
class RunOptions:
    """Everything the CLI can vary that lands in the pipeline config."""

    preset: str = "assembly"
    cores: int = 4
    mem_mb: int = 8000
    coassemble: bool = False
    host_fasta: Path | None = None
    skip_host_removal: bool = False
    allow_remote_inputs: bool = False
    binners: tuple[str, ...] = ("metabat", "maxbin", "vamb")
    checkm2_db: Path | None = None
    gtdbtk_db: Path | None = None
    taxonomy_db: Path | None = None
    taxdump: Path | None = None
    kegg_db: Path | None = None
    cog_db: Path | None = None
    pfam_db: Path | None = None
    demo: bool = False


def _fail(code: str, message: str, fix: str, **details: object) -> SkillError:
    return SkillError(
        stage=STAGE, error_code=code, message=message, fix=fix, details=dict(details)
    )


def _absolute(path: Path | str | None) -> str | None:
    """Absolute, symlink-resolved form of a user-supplied path."""
    if path is None or path == "":
        return None
    candidate = Path(path).expanduser()
    return str(candidate if candidate.is_absolute() else candidate.resolve())


def _missing_flags(stage: str, options: RunOptions) -> list[str]:
    return [flag for attribute, flag in STAGE_DATABASES[stage] if getattr(options, attribute) is None]


def _cap(config: dict, section: str, limit: int) -> None:
    """Cap a resource table without ever raising a value above the base config."""
    table = config.get(section)
    if not isinstance(table, dict):
        return
    for key, value in list(table.items()):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        table[key] = max(1, min(int(value), limit))


def build_config(
    base_config: dict,
    samples: list,
    options: RunOptions,
    *,
    samplesheet: Path,
    results_dir: Path,
    database_dir: Path,
) -> dict:
    """Return a new config for this run; `base_config` is left untouched."""
    config = copy.deepcopy(base_config)
    config["samplesheet"] = str(Path(samplesheet).expanduser().resolve())
    config["results_dir"] = str(Path(results_dir).expanduser().resolve())
    config["database_dir"] = str(Path(database_dir).expanduser().resolve())

    if options.demo:
        # The demo is upstream's own test configuration, which is deliberately a
        # reduced instance of the schema: QC, assembly, QUAST, MultiQC and
        # provenance, and nothing that needs a reference database. Only the three
        # paths move, so a demo run is comparable to upstream's run_test.sh.
        #
        # That makes the preset meaningless here. Applying preset switches on top
        # would stop it being the upstream test configuration, and ignoring the
        # preset would tell the user they ran something they did not, so the
        # contradiction is refused here as well as in the CLI.
        if options.preset != "assembly":
            raise _fail(
                ErrorCode.INVALID_PRESET,
                f"--demo cannot apply --preset {options.preset}: a demo run uses "
                "upstream's test configuration, which enables assembly and nothing "
                "else",
                "Drop --demo to run that preset against your own reads, or drop "
                "--preset to run the demo.",
                preset=options.preset,
                demo_base_config="workflow/config/config_test.yaml",
            )
        return config

    switches = PRESET_SWITCHES.get(options.preset)
    if switches is None:
        raise _fail(
            ErrorCode.INVALID_PRESET,
            f"unknown preset {options.preset!r}",
            f"Choose one of {', '.join(PRESETS)}.",
            preset=str(options.preset),
        )

    config.setdefault("profile", {})["profiler"] = switches["profiler"]
    for section in ("assembly", "binning", "heterogeneity", "mag", "taxonomy", "annotation"):
        config.setdefault(section, {})["enabled"] = switches[section]

    # A database-dependent stage is on only when every one of its databases is
    # supplied. Half a set is an error: silently dropping the stage would report
    # a run as complete while a requested analysis never happened.
    for stage in ("mag", "taxonomy", "annotation"):
        if not switches[stage]:
            continue
        missing = _missing_flags(stage, options)
        if not missing:
            continue
        supplied = [flag for attr, flag in STAGE_DATABASES[stage] if getattr(options, attr)]
        if not supplied:
            config[stage]["enabled"] = False
            continue
        message = (
            f"preset 'full' enables {stage}, which needs "
            + ", ".join(missing)
            + " as well as the "
            + ", ".join(supplied)
            + " already given"
        )
        if stage == "mag":
            # Without the GTDB-Tk reference the run fails at classification and
            # these tables are never written; the user must hear that here, not
            # after the tools are built and the reads are processed.
            message += (
                "; without it the run fails at the GTDB-Tk classification step "
                "and never writes final/binning/bins.tsv, "
                "final/diagnostic/bottleneck.tsv or final/benchmarks/runtime.tsv"
            )
        raise _fail(
            ErrorCode.MISSING_DATABASE,
            message,
            "Pass " + ", ".join(missing) + " for the " + stage + " stage, or narrow the "
            "run with --preset assembly so the stage is off and needs no database.",
            stage=stage,
            missing=missing,
        )

    # A database flag is only meaningful for a stage this preset enables. The
    # path would otherwise be written into a config that switches the stage off,
    # and the user would believe a reference was in use when no rule reads it.
    for stage, label in _STAGE_LABELS.items():
        if switches[stage]:
            continue
        flags = [flag for _, flag in STAGE_DATABASES[stage]]
        supplied = [
            f"{attribute}={_absolute(getattr(options, attribute))}"
            for attribute, _ in STAGE_DATABASES[stage]
            if getattr(options, attribute) is not None
        ]
        if not supplied:
            continue
        raise _fail(
            ErrorCode.INVALID_CONFIG,
            f"preset {options.preset!r} leaves {label} off, so "
            + " ".join(flags)
            + " cannot be used (given: "
            + ", ".join(supplied)
            + ")",
            f"Use --preset full to enable {label}, or drop {flags[0]} and the rest "
            "of its databases.",
            stage=stage,
            preset=options.preset,
            databases=supplied,
        )

    config.setdefault("assembly", {})["coassemble"] = bool(options.coassemble)
    config["assembly"].setdefault("spades", {})["isolate"] = False
    config.setdefault("binning", {})["binners"] = list(options.binners)

    decontaminate = config.setdefault("qc", {}).setdefault("decontaminate", {})
    decontaminate["enabled"] = not options.skip_host_removal
    host = _absolute(options.host_fasta)
    if host is not None:
        decontaminate["fasta"] = host

    remote = any(getattr(sample, "sra_run", None) for sample in samples)
    config.setdefault("fetch", {}).setdefault("sra", {})["enabled"] = bool(
        remote and options.allow_remote_inputs
    )

    config["mag"]["checkm2"]["db"] = _absolute(options.checkm2_db)
    config["mag"]["gtdbtk"]["db"] = _absolute(options.gtdbtk_db)
    config["taxonomy"]["db"] = _absolute(options.taxonomy_db)
    config["taxonomy"]["taxdump"] = _absolute(options.taxdump)
    config["annotation"]["kegg_db"] = _absolute(options.kegg_db)
    config["annotation"]["cog_db"] = _absolute(options.cog_db)
    config["annotation"]["pfam_db"] = _absolute(options.pfam_db)

    _cap(config, "threads", options.cores)
    _cap(config, "mem_mb", options.mem_mb)
    return config


def _require(condition: bool, code: str, message: str, fix: str, **details: object) -> None:
    if not condition:
        raise _fail(code, message, fix, **details)


def validate_config(config: dict, samples: list) -> None:
    """Every guard upstream applies, plus the two MAG databases it does not check."""
    _require(
        bool(samples),
        ErrorCode.INVALID_CONFIG,
        "no samples to run",
        "Check the samplesheet, or run with --demo.",
    )

    remote = [sample.sample for sample in samples if getattr(sample, "sra_run", None)]
    fetch_enabled = bool(
        config.get("fetch", {}).get("sra", {}).get("enabled", False)
    )
    if remote and not fetch_enabled:
        raise _fail(
            ErrorCode.REMOTE_INPUT_NOT_ALLOWED,
            f"samplesheet names public runs ({', '.join(remote)}) but fetch.sra.enabled "
            "is false",
            "Re-run with --allow-remote-inputs to download those runs, or supply local "
            "FASTQ paths instead of accessions.",
            samples=remote,
        )

    assembly = bool(config.get("assembly", {}).get("enabled", True))
    binning = bool(config.get("binning", {}).get("enabled", True))
    mag = bool(config.get("mag", {}).get("enabled", True))

    _require(
        not binning or assembly,
        ErrorCode.INVALID_CONFIG,
        "binning.enabled requires assembly.enabled=true",
        "Enable assembly (--preset assembly or --preset full), or turn binning off.",
    )
    _require(
        not mag or binning,
        ErrorCode.INVALID_CONFIG,
        "mag.enabled requires binning.enabled=true",
        "Enable binning (--preset full), or turn MAG assessment off.",
    )

    if assembly:
        single = [sample.sample for sample in samples if not getattr(sample, "paired", False)]
        _require(
            not single,
            ErrorCode.ASSEMBLY_REQUIRES_PAIRED,
            "metaSPAdes assembly requires paired-end samples; single-end library: "
            + ", ".join(single),
            "Use --preset qc or --preset profile for single-end libraries, or supply "
            "both mates.",
            samples=single,
        )
        _require(
            not bool(config.get("assembly", {}).get("spades", {}).get("isolate", False)),
            ErrorCode.INVALID_CONFIG,
            "assembly.spades.isolate is incompatible with this metagenome assembly "
            "workflow",
            "Set assembly.spades.isolate to false.",
        )

    if binners := config.get("binning", {}).get("binners", list(BINNERS)):
        unknown = [binner for binner in binners if binner not in BINNERS]
        _require(
            not unknown,
            ErrorCode.INVALID_CONFIG,
            f"binning.binners names unknown binner(s) {unknown}",
            f"Choose from {', '.join(sorted(BINNERS))}, or pass --binners.",
            unknown=unknown,
        )
        _require(
            len(binners) >= 2,
            ErrorCode.INVALID_CONFIG,
            f"binning.binners selects {list(binners)}; consensus binning needs at "
            "least two",
            "Pass at least two of --binners metabat,maxbin,vamb.",
            binners=list(binners),
        )

    if config.get("taxonomy", {}).get("enabled", True):
        sensitivity = config["taxonomy"].get("sensitivity")
        _require(
            sensitivity in DIAMOND_SENSITIVITIES,
            ErrorCode.INVALID_CONFIG,
            f"taxonomy.sensitivity is {sensitivity!r}; choose from "
            + ", ".join(DIAMOND_SENSITIVITIES),
            "Set taxonomy.sensitivity to a DIAMOND setting such as mid-sensitive.",
            sensitivity=sensitivity,
        )
        for key, flag in (("db", "--taxonomy-db"), ("taxdump", "--taxdump")):
            _require(
                bool(config["taxonomy"].get(key)),
                ErrorCode.MISSING_DATABASE,
                f"taxonomy is enabled but taxonomy.{key} is unset; both databases are "
                "external and are not bundled",
                f"Pass {flag}.",
                database=key,
            )

    if config.get("annotation", {}).get("enabled", True):
        missing = [
            f"annotation.{key}"
            for key in ("kegg_db", "cog_db", "pfam_db")
            if not config["annotation"].get(key)
        ]
        _require(
            not missing,
            ErrorCode.MISSING_DATABASE,
            "annotation is enabled but " + ", ".join(missing) + " is unset; all three "
            "databases are external and are not bundled",
            "Pass --kegg-db, --cog-db and --pfam-db, or narrow the run.",
            missing=missing,
        )

    if mag:
        checkm2 = config.get("mag", {}).get("checkm2", {}).get("db")
        gtdbtk = config.get("mag", {}).get("gtdbtk", {}).get("db")
        _require(
            bool(checkm2) and Path(str(checkm2)).is_file(),
            ErrorCode.MISSING_DATABASE,
            "MAG assessment needs mag.checkm2.db to point at a CheckM2 1.1.0 DIAMOND "
            f".dmnd file; {checkm2!r} is not one",
            "Pass --checkm2-db /path/to/checkm2.dmnd.",
            database="checkm2",
        )
        _require(
            bool(gtdbtk) and Path(str(gtdbtk)).is_dir(),
            ErrorCode.MISSING_DATABASE,
            "MAG assessment needs mag.gtdbtk.db to point at an extracted GTDB-Tk "
            f"reference directory; {gtdbtk!r} is not one. Without it the run fails "
            "at the GTDB-Tk classification step and never writes "
            "final/binning/bins.tsv, final/diagnostic/bottleneck.tsv or "
            "final/benchmarks/runtime.tsv",
            "Pass --gtdbtk-db /path/to/gtdbtk.",
            database="gtdbtk",
            unreachable_outputs=[
                "final/binning/bins.tsv",
                "final/diagnostic/bottleneck.tsv",
                "final/benchmarks/runtime.tsv",
            ],
        )


def write_config(config: dict, dest: Path) -> Path:
    """Write the effective config as YAML, preserving the upstream key order."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Effective configuration rendered by the ClawBio metagenome-mag-pipeline\n"
        "# skill. Derived from the pinned upstream workflow/config/config.yaml;\n"
        "# edit the skill's flags, not this file.\n"
    )
    dest.write_text(header + yaml.safe_dump(config, sort_keys=False, default_flow_style=False))
    return dest
