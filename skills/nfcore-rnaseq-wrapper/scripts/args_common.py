"""Shared CLI-arg-namespace construction for the nfcore_rnaseq_wrapper DAG steps.

ClawBio 1.x's ``nfcore_rnaseq_wrapper.py`` builds one big ``argparse.Namespace``
(~100 flags) in ``main()`` and threads it, mutated in place by a handful of
``_apply_*_overrides`` functions, through ``run_preflight`` /
``build_effective_params`` / ``build_nextflow_command``. Those three real
functions are duck-typed over that namespace (almost everything is read via
``getattr(args, field, default)``), so a much smaller namespace -- containing
only the flags this 2.0 port exposes as ``pipeline.yaml`` inputs -- drives the
same real code paths correctly; unset fields simply fall back to the same
defaults ``getattr`` already encodes.

Every ClawBio 2.0 step script imports ``build_args`` from here so the override
logic (``_apply_demo_overrides``, ``_sync_profile_flags``, ...) is ported once
and applied identically at every DAG step, matching how 1.x's ``main()``
applies it once before dispatching to preflight/samplesheet/params.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from schemas import DEFAULT_PIPELINE_VERSION, DEFAULT_PROFILE  # noqa: E402

SENTINEL = "none"


def none_or(value: str | None) -> str | None:
    """Map the pipeline.yaml sentinel string 'none' (and '', None) to None."""
    if value is None:
        return None
    v = value.strip()
    return None if v == "" or v.lower() == "none" else v


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


# Reference/index fields --demo and self-contained nf-core test profiles clear
# before params.yaml is written (ported from nfcore_rnaseq_wrapper.py's
# _DEMO_CLEARED_REFERENCE_FIELDS, trimmed to the reference fields this port
# exposes as pipeline.yaml inputs: genome/fasta/gtf/gff. 1.x additionally
# clears --igenomes-base/--transcript-fasta/--additional-fasta/--gene-bed/
# --splicesites/--*-index/--sortmerna-index, which this port does not expose.)
_DEMO_CLEARED_REFERENCE_FIELDS = ("genome", "fasta", "gtf", "gff")

_PROKARYOTIC_PROFILE_TOKENS = frozenset({"prokaryotic", "test_prokaryotic"})
_SELF_CONTAINED_TEST_PROFILES = frozenset(
    {"test", "test_full", "test_prokaryotic", "test_full_aws", "test_full_gcp", "test_full_azure", "test_gpu"}
)


def _apply_demo_overrides(args: SimpleNamespace) -> None:
    """Port of nfcore_rnaseq_wrapper.py's ``_apply_demo_overrides``."""
    if not getattr(args, "demo", False):
        return
    requested_aligner = getattr(args, "aligner", None)
    if requested_aligner and requested_aligner != "star_salmon":
        print(f"WARNING: --demo forces --aligner star_salmon (requested: {requested_aligner!r}).", file=sys.stderr)
    args.aligner = "star_salmon"
    if getattr(args, "input", None):
        print("WARNING: --demo ignores --input; the upstream test profile provides its own samplesheet.", file=sys.stderr)
        args.input = None
    cleared = [field for field in _DEMO_CLEARED_REFERENCE_FIELDS if getattr(args, field, None)]
    if cleared:
        flags = ", ".join("--" + field for field in cleared)
        print(
            f"WARNING: --demo ignores reference flags ({flags}); "
            "the upstream test profile provides its own bundled references.",
            file=sys.stderr,
        )
        for field in cleared:
            setattr(args, field, None)


def _sync_profile_flags(args: SimpleNamespace) -> None:
    """Port of nfcore_rnaseq_wrapper.py's ``_sync_profile_flags``."""
    profile_parts = {p.strip() for p in (getattr(args, "profile", "") or "").split(",") if p.strip()}
    if profile_parts & _PROKARYOTIC_PROFILE_TOKENS and not getattr(args, "prokaryotic", False):
        args.prokaryotic = True
    if "arm64" in profile_parts and not getattr(args, "arm", False):
        args.arm = True
    if profile_parts & _SELF_CONTAINED_TEST_PROFILES and not getattr(args, "_noinput", False):
        args._noinput = True


def _apply_noinput_overrides(args: SimpleNamespace) -> None:
    """Port of nfcore_rnaseq_wrapper.py's ``_apply_noinput_overrides``."""
    if not getattr(args, "_noinput", False) or getattr(args, "demo", False):
        return
    cleared = [field for field in _DEMO_CLEARED_REFERENCE_FIELDS if getattr(args, field, None)]
    if cleared:
        flags = ", ".join("--" + field for field in cleared)
        print(
            f"WARNING: a self-contained nf-core test profile owns its own references; "
            f"ignoring reference flags ({flags}).",
            file=sys.stderr,
        )
        for field in cleared:
            setattr(args, field, None)


def _apply_prokaryotic_overrides(args: SimpleNamespace) -> None:
    """Port of nfcore_rnaseq_wrapper.py's ``_apply_prokaryotic_overrides``.

    1.x distinguishes "user never passed --aligner" (None sentinel) from "user
    explicitly passed --aligner star_salmon" so --prokaryotic can safely default
    to bowtie2_salmon only in the former case. This 2.0 port's ``aligner`` input
    is ``required: true`` (pipeline.yaml has no per-run defaults -- see
    phylogenetics_builder for the precedent), so every run supplies an explicit
    value and that distinction cannot be reconstructed here. Consequently
    --prokaryotic does NOT auto-coerce the aligner in this port; pass
    ``aligner=bowtie2_salmon`` explicitly for prokaryotic runs (documented in
    SKILL.md Gotchas).
    """
    return


def _apply_aligner_default(args: SimpleNamespace) -> None:
    if getattr(args, "aligner", None) is None:
        args.aligner = "star_salmon"


def build_args(
    *,
    output: str,
    demo: bool = False,
    check: bool = True,
    input: str | None = None,
    aligner: str = "star_salmon",
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    gff: str | None = None,
    profile: str = DEFAULT_PROFILE,
    prokaryotic: bool = False,
    allow_remote_inputs: bool = False,
    pipeline_version: str = DEFAULT_PIPELINE_VERSION,
    allow_pipeline_version_override: bool = False,
    skip_alignment: bool = False,
    resume: bool = False,
    work_dir: str | None = None,
    timeout_hours: int = 12,
) -> SimpleNamespace:
    """Build the duck-typed args namespace shared by every DAG step's script.

    Every caller passes the exact same pipeline.yaml inputs (see SKILL.md), so
    every step re-derives an identical, deterministic namespace independently
    -- no step depends on another step's in-memory Python state, only on file
    outputs (matching the clawbio contract: state lives in files, not in a
    live process).
    """
    args = SimpleNamespace()
    args.output = output
    args.demo = as_bool(demo)
    args.check = as_bool(check)
    args.input = none_or(input)
    args.aligner = none_or(aligner)
    args._aligner_explicit = args.aligner is not None
    args.genome = none_or(genome)
    args.fasta = none_or(fasta)
    args.gtf = none_or(gtf)
    args.gff = none_or(gff)
    args.profile = profile
    args.prokaryotic = as_bool(prokaryotic)
    args.rapid_quant = False
    args.arm = False
    args.allow_remote_inputs = as_bool(allow_remote_inputs)
    args.pipeline_version = pipeline_version
    args.allow_pipeline_version_override = as_bool(allow_pipeline_version_override)
    args.skip_alignment = as_bool(skip_alignment)
    args.resume = as_bool(resume)
    args.work_dir = none_or(work_dir)
    args.timeout_hours = timeout_hours
    args._noinput = False

    _apply_demo_overrides(args)
    _sync_profile_flags(args)
    _apply_noinput_overrides(args)
    _apply_prokaryotic_overrides(args)
    _apply_aligner_default(args)
    return args
