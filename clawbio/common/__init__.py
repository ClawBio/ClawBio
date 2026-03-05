"""ClawBio common utilities — shared parsers, profiles, reports, checksums."""

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
]
