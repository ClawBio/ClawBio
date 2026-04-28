"""ClawBio common utilities — shared parsers, profiles, reports, checksums, reproducibility."""

from clawbio.common.parsers import (
    detect_format,
    parse_genetic_file,
    GenotypeRecord,
)
from clawbio.common.checksums import sha256_file, sha256_hex
from clawbio.common.report import (
    generate_report_header,
    generate_report_footer,
    DISCLAIMER,
)
from clawbio.common.profile import PatientProfile
from clawbio.common.html_report import HtmlReportBuilder, write_html_report
from clawbio.common.scrna_io import (
    compute_input_checksum,
    detect_processed_input_reason,
    load_count_adata,
    load_10x_mtx_data,
    resolve_input_source,
)
from clawbio.common.reproducibility import (
    write_checksums,
    write_environment_yml,
    write_commands_sh,
    write_conda_lock,
)

__all__ = [
    "detect_format",
    "parse_genetic_file",
    "GenotypeRecord",
    "sha256_file",
    "sha256_hex",
    "generate_report_header",
    "generate_report_footer",
    "DISCLAIMER",
    "PatientProfile",
    "HtmlReportBuilder",
    "write_html_report",
    "compute_input_checksum",
    "detect_processed_input_reason",
    "load_count_adata",
    "load_10x_mtx_data",
    "resolve_input_source",
    "write_checksums",
    "write_environment_yml",
    "write_commands_sh",
    "write_conda_lock",
]
