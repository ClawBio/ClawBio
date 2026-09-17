"""ClawBio common utilities — shared parsers, profiles, reports, checksums, reproducibility.

Attributes are resolved lazily (PEP 562). Importing a submodule directly, as in
`from clawbio.common.reproducibility import write_checksums`, must not drag numpy,
pandas and opentelemetry in through this file: skills that are otherwise
standard-library-only declare `pip_deps=[]` in their reproducibility recipe, and a
recipe that cannot import its own helper cannot replay the run it documents.
"""

from importlib import import_module

_EXPORTS = {
    "detect_format": "parsers",
    "parse_genetic_file": "parsers",
    "GenotypeRecord": "parsers",
    "sha256_file": "checksums",
    "sha256_hex": "checksums",
    "generate_report_header": "report",
    "generate_report_footer": "report",
    "DISCLAIMER": "report",
    "PatientProfile": "profile",
    "HtmlReportBuilder": "html_report",
    "write_html_report": "html_report",
    "compute_input_checksum": "scrna_io",
    "detect_processed_input_reason": "scrna_io",
    "load_count_adata": "scrna_io",
    "load_10x_mtx_data": "scrna_io",
    "resolve_input_source": "scrna_io",
    "write_checksums": "reproducibility",
    "write_environment_yml": "reproducibility",
    "write_commands_sh": "reproducibility",
    "write_conda_lock": "reproducibility",
    "SarekConfig": "sarek",
    "SarekSample": "sarek",
    "SarekWrapper": "sarek",
    "build_samplesheet": "sarek",
    "QcConfig": "vcf_qc",
    "QcResult": "vcf_qc",
    "VcfQC": "vcf_qc",
}


def __getattr__(name: str):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(import_module(f"{__name__}.{module}"), name)


__all__ = list(_EXPORTS)
