"""Shared argparse.Namespace-compatible defaults for the ported nf-core/scrnaseq
wrapper modules (preflight.py, params_builder.py, samplesheet_builder.py,
command_builder.py, pipeline_source.py).

Those modules were written against ClawBio 1.x's ``nfcore_scrnaseq_wrapper.py``
CLI, whose ~60-flag argparse.Namespace they read from (mostly via
``getattr(args, name, None)``, a handful of fields directly). ClawBio 2.0's
pipeline.yaml only exposes a curated subset of that surface as DAG inputs (see
../pipeline.yaml and ../SKILL.md's "Curated input surface" note) — mirroring
the source wrapper's own stated policy that it "is not a free-form
passthrough" of the upstream nf-core parameter set.

``build_args`` returns a ``types.SimpleNamespace`` pre-populated with every
attribute the ported modules read, set to the same defaults
``nfcore_scrnaseq_wrapper.py``'s ``build_parser()`` used, then overridden by
whatever keyword values the calling step script passes in. This keeps the
ported library modules byte-for-byte faithful (no internal call site needed
to change) while each step script only has to know about the handful of
fields it actually varies.
"""
from __future__ import annotations

from types import SimpleNamespace


def build_args(**overrides: object) -> SimpleNamespace:
    defaults: dict[str, object] = dict(
        verbose=False,
        no_banner=False,
        input=None,
        output=None,
        demo=False,
        check=False,
        profile="docker",
        pipeline_version="4.1.0",
        allow_dirty_pipeline=False,
        require_local_pipeline=False,
        allow_pipeline_version_override=False,
        trust_config_params=False,
        extra_config=[],
        preset="standard",
        preset_explicit=False,
        aligner=None,
        protocol=None,
        email=None,
        email_on_fail=None,
        multiqc_title=None,
        multiqc_config=None,
        multiqc_logo=None,
        multiqc_methods_description=None,
        publish_dir_mode=None,
        trace_report_suffix=None,
        monochrome_logs=False,
        expected_cells=None,
        timeout_hours=12.0,
        work_dir=None,
        allow_remote_inputs=False,
        allow_conda_cellranger=False,
        resume=False,
        skip_cellbender=False,
        skip_fastqc=False,
        skip_emptydrops=False,
        skip_multiqc=False,
        skip_cellranger_renaming=False,
        skip_cellrangermulti_vdjref=False,
        genome=None,
        igenomes_base=None,
        igenomes_ignore=False,
        save_reference=False,
        save_align_intermeds=None,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        star_feature=None,
        star_ignore_sjdbgtf=False,
        seq_center=None,
        simpleaf_umi_resolution=None,
        kb_workflow=None,
        kb_t1c=None,
        kb_t2c=None,
        motifs=None,
        cellrangerarc_config=None,
        cellrangerarc_reference=None,
        cellranger_vdj_index=None,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
        run_downstream=False,
        skip_downstream=False,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def parse_bool(value: str | None, *, default: bool = False) -> bool:
    if value is None:
        return default
    token = str(value).strip().lower()
    if token in ("true", "1", "yes", "y"):
        return True
    if token in ("false", "0", "no", "n", "", "none"):
        return False
    raise ValueError(f"expected a boolean-ish value, got {value!r}")


def parse_optional_str(value: str | None) -> str | None:
    if value is None:
        return None
    token = str(value).strip()
    if not token or token.lower() == "none":
        return None
    return token


def parse_optional_int(value: str | None) -> int | None:
    token = parse_optional_str(value)
    if token is None:
        return None
    return int(token)
