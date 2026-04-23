#!/usr/bin/env python3
"""
WES Clinical Report Generator (English) - ClawBio Skill
Generates professional clinical PDF reports from WES (Whole Exome Sequencing)
markdown reports. English language output.

Project: X202SC26016276-Z01-F001 (Novogene WES, 7 samples, hg38)
"""

import re
from pathlib import Path
from datetime import datetime

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
    from reportlab.platypus import (
        Paragraph, Spacer, Table, TableStyle,
        HRFlowable, NextPageTemplate, Image, Frame, PageTemplate,
        BaseDocTemplate,
    )
    from reportlab.platypus.doctemplate import _doNothing
    _REPORTLAB_AVAILABLE = True
except ImportError:
    colors = None
    A4 = None
    getSampleStyleSheet = None
    ParagraphStyle = None
    mm = None
    TA_LEFT = None
    TA_CENTER = None
    TA_JUSTIFY = None
    Paragraph = None
    Spacer = None
    Table = None
    TableStyle = None
    HRFlowable = None
    NextPageTemplate = None
    Image = None
    Frame = None
    PageTemplate = None
    BaseDocTemplate = None
    _doNothing = None
    _REPORTLAB_AVAILABLE = False

# ── Paths (defaults, overridable via CLI) ───────────────────────────────
REPORT_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/REPORTS")
OUTPUT_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/REPORTS/PDF-EN")
LOGO_LEFT = str(REPORT_DIR / "logo_predice.jpg")
LOGO_RIGHT = str(REPORT_DIR / "logo_inbiomedic.jpg")
ANCESTRY_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/ANCESTRY/OUTPUT")
SAMPLES = [f"Sample{i}" for i in range(1, 8)]

if _REPORTLAB_AVAILABLE:
    # ── Colour palette ──────────────────────────────────────────────────
    NAVY       = colors.HexColor("#1B2A4A")
    DARK_BLUE  = colors.HexColor("#2C4A7C")
    MID_BLUE   = colors.HexColor("#4A90D9")
    LIGHT_BLUE = colors.HexColor("#D6E8F7")
    PALE_BLUE  = colors.HexColor("#EBF3FB")
    ACCENT_RED = colors.HexColor("#C0392B")
    ACCENT_AMBER = colors.HexColor("#E67E22")
    ACCENT_GREEN = colors.HexColor("#27AE60")
    WARM_GREY  = colors.HexColor("#7F8C8D")
    LIGHT_GREY = colors.HexColor("#F4F6F8")
    WHITE      = colors.white
    BLACK      = colors.black
    ROW_ALT    = colors.HexColor("#F8FAFC")

    PAGE_W, PAGE_H = A4
    MARGIN = 20 * mm
else:
    NAVY = DARK_BLUE = MID_BLUE = LIGHT_BLUE = PALE_BLUE = None
    ACCENT_RED = ACCENT_AMBER = ACCENT_GREEN = None
    WARM_GREY = LIGHT_GREY = WHITE = BLACK = ROW_ALT = None
    PAGE_W = PAGE_H = MARGIN = None


def _require_reportlab():
    """Raise a helpful error when PDF generation is requested without reportlab."""
    if not _REPORTLAB_AVAILABLE:
        raise ImportError(
            "reportlab is required for PDF generation. Install it with "
            "`pip install -r skills/wes-clinical-report-en/requirements.txt`."
        )


# ── Markdown parser ─────────────────────────────────────────────────────

def parse_markdown_report(filepath):
    """Parse a WES markdown report into structured sections."""
    text = filepath.read_text(encoding="utf-8")
    report = {"sample_id": "", "metadata": {}, "sections": []}

    m = re.search(r"Report:\s*(Sample\d+)", text)
    if m:
        report["sample_id"] = m.group(1)

    meta_block = re.search(r">\s*\*\*Project\*\*(.+?)(?=\n---|\n\n##)", text, re.DOTALL)
    if meta_block:
        raw = meta_block.group(0)
        for pair in re.findall(r"\*\*(\w[\w\s]*?)\*\*\s*(.+?)(?=\s*\||\s*$)", raw, re.MULTILINE):
            report["metadata"][pair[0].strip()] = pair[1].strip()

    section_splits = re.split(r"\n## (\d+\.\s+.+)", text)
    for i in range(1, len(section_splits), 2):
        heading = section_splits[i].strip()
        body = section_splits[i + 1] if i + 1 < len(section_splits) else ""
        report["sections"].append({"heading": heading, "body": body.strip()})

    return report


def parse_table(text):
    """Parse a markdown table into headers and rows."""
    lines = [l.strip() for l in text.strip().split("\n") if l.strip()]
    if len(lines) < 2:
        return [], []
    headers = [c.strip().replace("**", "") for c in lines[0].split("|") if c.strip()]
    rows = []
    for line in lines[2:]:
        if "|" not in line:
            break
        cells = [c.strip().replace("**", "") for c in line.split("|") if c.strip() != ""]
        if cells:
            rows.append(cells)
    return headers, rows


def extract_tables_and_text(body):
    """Split section body into alternating text and table blocks."""
    parts = []
    lines = body.split("\n")
    buf = []
    in_table = False
    table_buf = []

    for line in lines:
        stripped = line.strip()
        is_table_line = stripped.startswith("|") and "|" in stripped[1:]

        if is_table_line:
            if not in_table:
                if buf:
                    parts.append(("text", "\n".join(buf)))
                    buf = []
                in_table = True
                table_buf = []
            table_buf.append(stripped)
        else:
            if in_table:
                parts.append(("table", "\n".join(table_buf)))
                table_buf = []
                in_table = False
            buf.append(line)

    if in_table and table_buf:
        parts.append(("table", "\n".join(table_buf)))
    if buf:
        joined = "\n".join(buf).strip()
        if joined:
            parts.append(("text", joined))

    return parts


# ── PDF styles ──────────────────────────────────────────────────────────

def build_styles():
    ss = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle(
            "ClinTitle", parent=ss["Title"],
            fontName="Helvetica-Bold", fontSize=22,
            textColor=NAVY, spaceAfter=4 * mm, alignment=TA_LEFT,
        ),
        "subtitle": ParagraphStyle(
            "ClinSubtitle", parent=ss["Normal"],
            fontName="Helvetica", fontSize=10,
            textColor=WARM_GREY, spaceAfter=6 * mm,
        ),
        "heading1": ParagraphStyle(
            "ClinH1", parent=ss["Heading1"],
            fontName="Helvetica-Bold", fontSize=14,
            textColor=NAVY, spaceBefore=8 * mm, spaceAfter=3 * mm,
            borderWidth=0, borderPadding=0,
        ),
        "heading2": ParagraphStyle(
            "ClinH2", parent=ss["Heading2"],
            fontName="Helvetica-Bold", fontSize=11,
            textColor=DARK_BLUE, spaceBefore=5 * mm, spaceAfter=2 * mm,
        ),
        "body": ParagraphStyle(
            "ClinBody", parent=ss["Normal"],
            fontName="Helvetica", fontSize=9,
            textColor=BLACK, leading=13, alignment=TA_JUSTIFY,
            spaceAfter=2 * mm,
        ),
        "body_bold": ParagraphStyle(
            "ClinBodyBold", parent=ss["Normal"],
            fontName="Helvetica-Bold", fontSize=9,
            textColor=BLACK, leading=13, spaceAfter=2 * mm,
        ),
        "bullet": ParagraphStyle(
            "ClinBullet", parent=ss["Normal"],
            fontName="Helvetica", fontSize=9,
            textColor=BLACK, leading=13,
            leftIndent=12, bulletIndent=4,
            spaceBefore=1, spaceAfter=1,
        ),
        "small": ParagraphStyle(
            "ClinSmall", parent=ss["Normal"],
            fontName="Helvetica", fontSize=7.5,
            textColor=WARM_GREY, leading=10, spaceAfter=1 * mm,
        ),
        "disclaimer": ParagraphStyle(
            "ClinDisclaimer", parent=ss["Normal"],
            fontName="Helvetica-Oblique", fontSize=7.5,
            textColor=ACCENT_RED, leading=10,
            spaceBefore=4 * mm, spaceAfter=2 * mm,
        ),
        "table_header": ParagraphStyle(
            "ClinTH", fontName="Helvetica-Bold", fontSize=7.5,
            textColor=WHITE, leading=10, alignment=TA_LEFT,
        ),
        "table_cell": ParagraphStyle(
            "ClinTD", fontName="Helvetica", fontSize=7.5,
            textColor=BLACK, leading=10, alignment=TA_LEFT,
        ),
        "table_cell_bold": ParagraphStyle(
            "ClinTDBold", fontName="Helvetica-Bold", fontSize=7.5,
            textColor=BLACK, leading=10, alignment=TA_LEFT,
        ),
        "kpi_number": ParagraphStyle(
            "KPINum", fontName="Helvetica-Bold", fontSize=18,
            textColor=NAVY, alignment=TA_CENTER, leading=22,
        ),
        "kpi_label": ParagraphStyle(
            "KPILabel", fontName="Helvetica", fontSize=7.5,
            textColor=WARM_GREY, alignment=TA_CENTER, leading=10,
        ),
        "footer": ParagraphStyle(
            "ClinFooter", fontName="Helvetica", fontSize=7,
            textColor=WARM_GREY, alignment=TA_CENTER,
        ),
    }
    return styles


# ── Page templates ──────────────────────────────────────────────────────

def _draw_logos(canvas, y_center):
    """Draw institutional logos (left and right) centred on y_center."""
    pred_h = 16 * mm
    pred_w = pred_h * (166 / 100)
    pred_x = MARGIN
    pred_y = y_center - pred_h / 2
    try:
        canvas.drawImage(LOGO_LEFT, pred_x, pred_y, pred_w, pred_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    inb_h = 14 * mm
    inb_w = inb_h * (431 / 133)
    inb_x = PAGE_W - MARGIN - inb_w
    inb_y = y_center - inb_h / 2
    try:
        canvas.drawImage(LOGO_RIGHT, inb_x, inb_y, inb_w, inb_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass


def header_footer_en(canvas, doc):
    """Header and footer for body pages (English)."""
    canvas.saveState()

    # Small logo strip
    logo_bar_h = 12 * mm
    logo_bar_y = PAGE_H - logo_bar_h
    canvas.setFillColor(WHITE)
    canvas.rect(0, logo_bar_y, PAGE_W, logo_bar_h, fill=1, stroke=0)

    pred_h = 8 * mm
    pred_w = pred_h * (166 / 100)
    try:
        canvas.drawImage(LOGO_LEFT, MARGIN, logo_bar_y + 2 * mm,
                         pred_w, pred_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    inb_h = 7 * mm
    inb_w = inb_h * (431 / 133)
    try:
        canvas.drawImage(LOGO_RIGHT,
                         PAGE_W - MARGIN - inb_w, logo_bar_y + 2.5 * mm,
                         inb_w, inb_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    # Navy header bar
    navy_bar_y = logo_bar_y - 11 * mm
    canvas.setFillColor(NAVY)
    canvas.rect(0, navy_bar_y, PAGE_W, 11 * mm, fill=1, stroke=0)

    canvas.setFillColor(WHITE)
    canvas.setFont("Helvetica-Bold", 8)
    canvas.drawString(MARGIN, navy_bar_y + 3.5 * mm,
                      "WHOLE EXOME SEQUENCING CLINICAL REPORT")

    canvas.setFont("Helvetica", 8)
    sample_id = getattr(doc, '_sample_id', '')
    canvas.drawRightString(PAGE_W - MARGIN, navy_bar_y + 3.5 * mm,
                           sample_id)

    # Accent line
    canvas.setStrokeColor(MID_BLUE)
    canvas.setLineWidth(1.5)
    canvas.line(0, navy_bar_y - 0.5 * mm, PAGE_W, navy_bar_y - 0.5 * mm)

    # Footer
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(WARM_GREY)
    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 12 * mm, PAGE_W - MARGIN, 12 * mm)

    canvas.drawString(MARGIN, 8 * mm, "ClawBio WES Analysis Pipeline v1.0")
    canvas.drawCentredString(PAGE_W / 2, 8 * mm,
                             "Project X202SC26016276-Z01-F001")
    canvas.drawRightString(PAGE_W - MARGIN, 8 * mm,
                           f"Page {doc.page}")

    canvas.setFont("Helvetica-Oblique", 6)
    canvas.setFillColor(colors.HexColor("#BDC3C7"))
    canvas.drawCentredString(PAGE_W / 2, 4 * mm,
                             "CONFIDENTIAL - FOR RESEARCH USE ONLY"
                             " - NOT FOR CLINICAL DIAGNOSIS")
    canvas.restoreState()


def first_page_header_en(canvas, doc):
    """Cover page header (English)."""
    canvas.saveState()

    # Logo strip
    logo_strip_h = 22 * mm
    logo_strip_y = PAGE_H - logo_strip_h
    canvas.setFillColor(WHITE)
    canvas.rect(0, logo_strip_y, PAGE_W, logo_strip_h, fill=1, stroke=0)
    _draw_logos(canvas, logo_strip_y + logo_strip_h / 2)

    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.8)
    canvas.line(MARGIN, logo_strip_y, PAGE_W - MARGIN, logo_strip_y)

    # Navy header block
    header_h = 48 * mm
    header_y = logo_strip_y - header_h
    canvas.setFillColor(NAVY)
    canvas.rect(0, header_y, PAGE_W, header_h, fill=1, stroke=0)

    # Accent stripe
    canvas.setFillColor(MID_BLUE)
    canvas.rect(0, header_y - 2 * mm, PAGE_W, 2 * mm, fill=1, stroke=0)

    # Title
    canvas.setFillColor(WHITE)
    canvas.setFont("Helvetica-Bold", 24)
    canvas.drawString(MARGIN, header_y + header_h - 16 * mm,
                      "Whole Exome")
    canvas.drawString(MARGIN, header_y + header_h - 26 * mm,
                      "Sequencing Report")

    # Sample ID badge
    sample_id = getattr(doc, '_sample_id', '')
    canvas.setFont("Helvetica-Bold", 14)
    canvas.setFillColor(MID_BLUE)
    badge_w = 50 * mm
    badge_h = 12 * mm
    badge_x = PAGE_W - MARGIN - badge_w
    badge_y = header_y + header_h - 29 * mm
    canvas.roundRect(badge_x, badge_y, badge_w, badge_h, 3 * mm,
                     fill=1, stroke=0)
    canvas.setFillColor(WHITE)
    canvas.drawCentredString(badge_x + badge_w / 2, badge_y + 3 * mm,
                             sample_id)

    # Metadata
    canvas.setFont("Helvetica", 9)
    canvas.setFillColor(colors.HexColor("#A0B4CC"))
    canvas.drawString(MARGIN, header_y + 9 * mm,
                      "Illumina NovaSeq PE150  |  Xplus WES (60.5 Mb)  |  "
                      "GRCh38/hg38  |  GATK HaplotypeCaller 4.3.0")
    now = datetime.now()
    date_str = now.strftime("%d %B %Y")
    canvas.drawString(MARGIN, header_y + 3 * mm,
                      f"Report generated: {date_str}  |  "
                      "Analysis: ClawBio + ANNOVAR + ClinPGx")

    # Footer
    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 12 * mm, PAGE_W - MARGIN, 12 * mm)

    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(WARM_GREY)
    canvas.drawString(MARGIN, 8 * mm, "ClawBio WES Analysis Pipeline v1.0")
    canvas.drawCentredString(PAGE_W / 2, 8 * mm,
                             "Project X202SC26016276-Z01-F001")
    canvas.drawRightString(PAGE_W - MARGIN, 8 * mm, f"Page {doc.page}")

    canvas.setFont("Helvetica-Oblique", 6)
    canvas.setFillColor(colors.HexColor("#BDC3C7"))
    canvas.drawCentredString(PAGE_W / 2, 4 * mm,
                             "CONFIDENTIAL - FOR RESEARCH USE ONLY"
                             " - NOT FOR CLINICAL DIAGNOSIS")
    canvas.restoreState()


# ── KPI cards ───────────────────────────────────────────────────────────

def build_kpi_row(metrics, styles):
    cells = []
    for label, value in metrics.items():
        cell_content = [
            Paragraph(str(value), styles["kpi_number"]),
            Paragraph(label, styles["kpi_label"]),
        ]
        cells.append(cell_content)

    n = len(cells)
    col_w = (PAGE_W - 2 * MARGIN) / n
    data = [cells]

    t = Table(data, colWidths=[col_w] * n, rowHeights=[22 * mm])
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("BOX", (0, 0), (-1, -1), 0.5, LIGHT_BLUE),
        ("INNERGRID", (0, 0), (-1, -1), 0.5, LIGHT_BLUE),
        ("BACKGROUND", (0, 0), (-1, -1), PALE_BLUE),
        ("TOPPADDING", (0, 0), (-1, -1), 3 * mm),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2 * mm),
        ("ROUNDEDCORNERS", [3, 3, 3, 3]),
    ]))
    return t


# ── Table builder ───────────────────────────────────────────────────────

def build_clinical_table(headers, rows, styles):
    usable_w = PAGE_W - 2 * MARGIN

    header_paras = [Paragraph(h, styles["table_header"]) for h in headers]

    data_paras = []
    for row in rows:
        para_row = []
        for j, cell in enumerate(row):
            cell_text = cell[:120] + "..." if len(cell) > 120 else cell
            if j == 0:
                para_row.append(Paragraph(cell_text, styles["table_cell_bold"]))
            else:
                para_row.append(Paragraph(cell_text, styles["table_cell"]))
        data_paras.append(para_row)

    table_data = [header_paras] + data_paras
    ncols = len(headers)

    if ncols <= 2:
        col_widths = [usable_w / ncols] * ncols
    else:
        base = usable_w / (ncols + 1)
        col_widths = [base] * (ncols - 1) + [base * 2]

    total = sum(col_widths)
    if total > usable_w:
        scale = usable_w / total
        col_widths = [w * scale for w in col_widths]

    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), DARK_BLUE),
        ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 7.5),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 3 * mm),
        ("TOPPADDING", (0, 0), (-1, 0), 3 * mm),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 1), (-1, -1), 7.5),
        ("TOPPADDING", (0, 1), (-1, -1), 2 * mm),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 2 * mm),
        ("LEFTPADDING", (0, 0), (-1, -1), 3 * mm),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2 * mm),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#DEE2E6")),
        ("LINEBELOW", (0, 0), (-1, 0), 1.2, NAVY),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, ROW_ALT]),
    ]

    t = Table(table_data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle(style_cmds))
    return t


# ── KPI extraction ──────────────────────────────────────────────────────

def extract_exome_kpis(section_body):
    """Extract key performance indicators from the Exome Summary section."""
    kpis = {}
    parts = extract_tables_and_text(section_body)
    for kind, content in parts:
        if kind == "table":
            headers, rows = parse_table(content)
            break
    else:
        return kpis

    for row in rows:
        if len(row) >= 2:
            label, value = row[0].strip(), row[1].strip().replace(",", "")
            if "Total SNP" in label:
                kpis["Total SNPs"] = f"{int(value):,}"
            elif "Missense" == label:
                kpis["Missense"] = f"{int(value):,}"
            elif "ClinVar" in label:
                kpis["ClinVar P/LP"] = value
            elif "Rare + comp" in label or "Rare + comput" in label:
                kpis["Rare Damaging"] = value
            elif "Stopgain" in label:
                kpis["Stopgain"] = value

    return kpis


# ── Data extraction helpers ─────────────────────────────────────────────

def _extract_pathogenic_variants(report):
    """Extract P/LP variants from section 2."""
    variants = []
    for sec in report["sections"]:
        if "Clinically Significant" in sec["heading"]:
            parts = extract_tables_and_text(sec["body"])
            for kind, content in parts:
                if kind == "table":
                    headers, rows = parse_table(content)
                    if "Classification" in headers:
                        cls_idx = headers.index("Classification")
                        gene_idx = headers.index("Gene") if "Gene" in headers else 0
                        var_idx = headers.index("Variant") if "Variant" in headers else 1
                        zyg_idx = headers.index("Zygosity") if "Zygosity" in headers else 2
                        cond_idx = headers.index("Associated Condition") if "Associated Condition" in headers else -1
                        for row in rows:
                            if len(row) > cls_idx:
                                cls = row[cls_idx].strip()
                                if "Pathogenic" in cls and "Conflicting" not in cls:
                                    variants.append({
                                        "gene": row[gene_idx].strip(),
                                        "variant": row[var_idx].strip() if len(row) > var_idx else "",
                                        "zygosity": row[zyg_idx].strip() if len(row) > zyg_idx else "",
                                        "condition": row[cond_idx].strip() if cond_idx >= 0 and len(row) > cond_idx else "",
                                    })
            break
    return variants


def _extract_pgx_alerts(report):
    """Extract PGx findings from section 3."""
    alerts = []
    for sec in report["sections"]:
        if "Pharmacogenomics" in sec["heading"]:
            parts = extract_tables_and_text(sec["body"])
            for kind, content in parts:
                if kind == "table":
                    headers, rows = parse_table(content)
                    if "Clinical Effect" in headers:
                        gene_idx = headers.index("Gene") if "Gene" in headers else 0
                        allele_idx = headers.index("Allele") if "Allele" in headers else 2
                        zyg_idx = headers.index("Zygosity") if "Zygosity" in headers else 3
                        eff_idx = headers.index("Clinical Effect")
                        med_idx = headers.index("Affected Medications") if "Affected Medications" in headers else -1
                        for row in rows:
                            if len(row) > eff_idx:
                                alerts.append({
                                    "gene": row[gene_idx].strip(),
                                    "allele": row[allele_idx].strip() if len(row) > allele_idx else "",
                                    "zygosity": row[zyg_idx].strip() if len(row) > zyg_idx else "",
                                    "effect": row[eff_idx].strip(),
                                    "meds": row[med_idx].strip() if med_idx >= 0 and len(row) > med_idx else "",
                                })
            break
    return alerts


def _extract_rare_damaging(report):
    """Extract top rare damaging variants from section 5."""
    variants = []
    for sec in report["sections"]:
        if "Rare Damaging" in sec["heading"] or "Prioritised" in sec["heading"]:
            parts = extract_tables_and_text(sec["body"])
            for kind, content in parts:
                if kind == "table":
                    headers, rows = parse_table(content)
                    gene_idx = headers.index("Gene") if "Gene" in headers else 0
                    revel_idx = headers.index("REVEL") if "REVEL" in headers else -1
                    omim_idx = -1
                    for i, h in enumerate(headers):
                        if "OMIM" in h:
                            omim_idx = i
                            break
                    for row in rows[:5]:
                        if len(row) > gene_idx:
                            revel = row[revel_idx].strip() if revel_idx >= 0 and len(row) > revel_idx else ""
                            omim = row[omim_idx].strip() if omim_idx >= 0 and len(row) > omim_idx else ""
                            variants.append({
                                "gene": row[gene_idx].strip(),
                                "revel": revel,
                                "omim": omim[:80],
                            })
            break
    return variants


def _extract_metrics(report):
    """Extract key metrics from section 1."""
    metrics = {}
    for sec in report["sections"]:
        if "Exome Summary" in sec["heading"]:
            parts = extract_tables_and_text(sec["body"])
            for kind, content in parts:
                if kind == "table":
                    headers, rows = parse_table(content)
                    for row in rows:
                        if len(row) >= 2:
                            metrics[row[0].strip()] = row[1].strip()
                    break

            m = re.search(r"ratio\s+is\s+\*?\*?(\d+\.\d+)", sec["body"])
            if m:
                metrics["hom_het_ratio"] = float(m.group(1))
            break
    return metrics


def _extract_acmg_cancer(report):
    """Extract ACMG and cancer panel counts from section 2."""
    info = {"acmg_plp": 0, "cancer_plp": 0}
    for sec in report["sections"]:
        if "Clinically Significant" in sec["heading"]:
            m = re.search(r"of which (\d+) have ClinVar P/LP", sec["body"])
            if m:
                info["acmg_plp"] = int(m.group(1))
            m = re.search(r"Cancer predisposition panel:\s*(\d+)", sec["body"])
            if m:
                info["cancer_plp"] = int(m.group(1))
            break
    return info


# ── Interpretation paragraph (English) ──────────────────────────────────

def build_interpretation_paragraph(report, styles):
    """Build an interpretive summary paragraph in English."""
    elements = []

    pathogenic = _extract_pathogenic_variants(report)
    pgx = _extract_pgx_alerts(report)
    rare_dmg = _extract_rare_damaging(report)
    metrics = _extract_metrics(report)
    acmg = _extract_acmg_cancer(report)

    hom_het = metrics.get("hom_het_ratio", 0)
    rare_count = metrics.get("Rare + computationally damaging", "0")
    sample_id = report["sample_id"]

    sentences = []

    # Pathogenic variants
    if pathogenic:
        gene_list = ", ".join(v["gene"] for v in pathogenic)
        hom_genes = [v["gene"] for v in pathogenic if v["zygosity"] == "Hom"]
        het_genes = [v["gene"] for v in pathogenic if v["zygosity"] == "Het"]
        sentences.append(
            f"The exome analysis of {sample_id} identified "
            f"{len(pathogenic)} variant(s) classified as pathogenic "
            f"in ClinVar, located in the genes {gene_list}."
        )
        conditions_seen = set()
        for v in pathogenic:
            cond = v["condition"]
            if cond and cond not in conditions_seen:
                conditions_seen.add(cond)
        if conditions_seen:
            cond_text = "; ".join(c[:80] for c in sorted(conditions_seen))
            sentences.append(
                f"These variants are associated with the following conditions: "
                f"{cond_text}."
            )
        if hom_genes:
            sentences.append(
                f"Variants in {', '.join(hom_genes)} are homozygous, "
                f"indicating that both copies of the gene are affected "
                f"and the clinical effect may be more pronounced."
            )
    else:
        sentences.append(
            f"The exome analysis of {sample_id} did not identify any "
            f"variants classified as confirmed pathogenic in ClinVar."
        )

    # ACMG / Cancer
    if acmg["acmg_plp"] > 0 or acmg["cancer_plp"] > 0:
        parts = []
        if acmg["acmg_plp"] > 0:
            parts.append(
                f"{acmg['acmg_plp']} variant(s) in medically actionable "
                f"genes according to ACMG SF v3.2 guidelines"
            )
        if acmg["cancer_plp"] > 0:
            parts.append(
                f"{acmg['cancer_plp']} variant(s) in cancer "
                f"predisposition genes"
            )
        sentences.append(
            "Additionally, " + ", and ".join(parts) +
            " were identified, which may require specialised clinical follow-up."
        )

    # PGx
    if pgx:
        pgx_parts = []
        for p in pgx:
            gene = p["gene"]
            eff = p["effect"].lower()
            meds = p["meds"]
            if meds:
                pgx_parts.append(f"{gene} ({eff}, affecting {meds})")
            else:
                pgx_parts.append(f"{gene} ({eff})")
        sentences.append(
            f"From a pharmacogenomic perspective, relevant variants were "
            f"detected in {', '.join(pgx_parts)}. "
            f"These findings may influence individual drug response and "
            f"should be considered before prescribing the indicated medications."
        )

    # Homozygosity
    if hom_het > 2.0:
        sentences.append(
            f"An elevated homozygosity pattern was observed "
            f"(Hom/Het ratio = {hom_het:.2f}, compared to the expected ~1.5 "
            f"in outbred populations), which may reflect specific population "
            f"structure or parental consanguinity, and should be taken into "
            f"account when interpreting recessive disease risk."
        )

    # Rare damaging
    if rare_count and int(rare_count) > 0:
        top_genes = ", ".join(v["gene"] for v in rare_dmg[:3]) if rare_dmg else ""
        sent = (
            f"{rare_count} rare variants with computational pathogenicity "
            f"predictions (REVEL > 0.5 or CADD > 20) were prioritised"
        )
        if top_genes:
            sent += f", including variants in {top_genes}"
        sent += (
            ". These variants require correlation with the patient's clinical "
            "phenotype to determine their diagnostic relevance."
        )
        sentences.append(sent)

    # Closing
    sentences.append(
        "It is recommended that clinically relevant variants be validated "
        "by Sanger sequencing or another independent methodology, and that "
        "the patient be referred to a clinical genetics service for "
        "counselling and family evaluation as appropriate."
    )

    # Render
    elements.append(Paragraph(
        "<b>Results Interpretation</b>",
        styles["heading2"]
    ))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=MID_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    full_text = " ".join(sentences)
    elements.append(Paragraph(full_text, styles["body"]))
    elements.append(Spacer(1, 4 * mm))

    return elements


# ── Ancestry section ────────────────────────────────────────────────────

def _build_ancestry_section(sample_id, styles):
    """Build ancestry section showing admixture percentages."""
    import json

    elements = []
    elements.append(Paragraph(
        f'<font color="{MID_BLUE.hexval()}">8.</font> '
        f'Ancestry Estimation',
        styles["heading1"]
    ))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=LIGHT_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    json_path = ANCESTRY_DIR / "ancestry_results.json"
    ancestry_data = None
    if json_path.exists():
        with open(json_path) as f:
            all_data = json.load(f)
            ancestry_data = all_data.get(sample_id)

    if ancestry_data:
        props = ancestry_data["proportions"]
        n_markers = ancestry_data["n_markers"]
        sorted_pops = sorted(props.items(), key=lambda x: -x[1])

        top_pop, top_pct = sorted_pops[0]
        components = [f"{pop} ({pct:.1f}%)" for pop, pct in sorted_pops if pct >= 0.5]

        interp = (
            f"The ancestral composition of {sample_id} was estimated using "
            f"non-negative least squares (NNLS) decomposition of observed "
            f"allele frequencies against continental reference populations "
            f"from the 1000 Genomes Project, using "
            f"{n_markers:,} ancestry-informative markers (AIMs). "
            f"The predominant component is <b>{top_pop} ({top_pct:.1f}%)</b>, "
            f"followed by {', '.join(components[1:])}."
        )
        elements.append(Paragraph(interp, styles["body"]))
        elements.append(Spacer(1, 3 * mm))

        headers = ["Ancestral Category", "Proportion (%)"]
        rows = [[pop, f"{pct:.1f}%"] for pop, pct in sorted_pops]
        t = build_clinical_table(headers, rows, styles)
        elements.append(t)
        elements.append(Spacer(1, 2 * mm))

        pie_path = ANCESTRY_DIR / f"{sample_id}_ancestry_pie.png"
        if pie_path.exists():
            elements.append(Spacer(1, 2 * mm))
            usable_w = PAGE_W - 2 * MARGIN
            img_w = usable_w * 0.65
            img_h = img_w * 0.78
            img = Image(str(pie_path), width=img_w, height=img_h)
            img.hAlign = "CENTER"
            elements.append(img)
            elements.append(Spacer(1, 2 * mm))

        elements.append(Paragraph(
            f"<i>Estimation based on {n_markers:,} ancestry-informative markers "
            f"(AIMs) with known population frequencies from the 1000 Genomes "
            f"Project. Method: non-negative least squares (NNLS).</i>",
            styles["small"]
        ))
    else:
        elements.append(Paragraph(
            "No ancestry data available for this sample.",
            styles["body"]
        ))

    elements.append(Spacer(1, 2 * mm))
    return elements


# ── Limitations section ─────────────────────────────────────────────────

def _build_limitations_section(styles):
    """Build limitations section in English."""
    elements = []
    elements.append(Paragraph(
        f'<font color="{MID_BLUE.hexval()}">9.</font> Limitations',
        styles["heading1"]
    ))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=LIGHT_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    limitations = [
        (
            "Exome coverage",
            "Whole exome sequencing (WES) captures approximately 1-2% of the "
            "genome, corresponding to protein-coding regions. Variants in deep "
            "intronic, intergenic, regulatory, or non-coding DNA regions are "
            "not detected by this methodology. Large structural variants "
            "(deletions, duplications, inversions) and repeat expansions are "
            "also not reliably detected by WES."
        ),
        (
            "Variant classification",
            "Pathogenicity classifications are derived from ClinVar and "
            "reflect the consensus of submitting laboratories at the time "
            "of annotation. These classifications may change as new evidence "
            "accumulates. Variants of uncertain significance (VUS) should "
            "not be interpreted as pathogenic or benign."
        ),
        (
            "Pharmacogenomics",
            "Pharmacogenomic star alleles are inferred from individual "
            "variants (SNVs) and do not include copy number analysis (CNV), "
            "which is particularly relevant for CYP2D6. Metaboliser phenotype "
            "assignment requires clinical validation through certified "
            "pharmacogenomic assays before modifying any prescription."
        ),
        (
            "Ancestry estimation",
            "Ancestry is estimated by correlating allele frequencies with "
            "continental reference populations (1000 Genomes, gnomAD). This "
            "method provides a continental-level approximation but does not "
            "replace a formal principal component analysis (PCA) with a dense "
            "reference panel such as the Simons Genome Diversity Project. "
            "Low variant overlap between WES and genomic reference panels "
            "limits subcontinental resolution."
        ),
        (
            "Fitness and nutrition traits",
            "Genotype-phenotype associations for nutrition and fitness traits "
            "are based on observational studies with generally modest effect "
            "sizes. These findings should not be used as the sole basis for "
            "dietary or exercise recommendations. Evidence grades (A, B, C) "
            "reflect the level of replication in the literature, not the "
            "magnitude of the effect."
        ),
        (
            "Homozygosity and consanguinity",
            "Elevated homozygosity ratios observed across samples may be due "
            "to population structure, consanguinity, or inherent exome capture "
            "bias. A formal analysis with whole-genome data and inbreeding "
            "coefficients (F-statistics) would be required to distinguish "
            "between these causes."
        ),
        (
            "Clinical validation",
            "This report is generated for research purposes and has not been "
            "validated according to accredited clinical laboratory standards "
            "(ISO 15189). All variants with potential clinical relevance must "
            "be confirmed by Sanger sequencing or another independent "
            "methodology in a certified laboratory before making medical "
            "decisions."
        ),
    ]

    for title, text in limitations:
        elements.append(Paragraph(
            f"\u2022  <b>{title}</b>: {text}",
            styles["bullet"]
        ))
        elements.append(Spacer(1, 1.5 * mm))

    elements.append(Spacer(1, 2 * mm))
    return elements


# ── Section rendering ───────────────────────────────────────────────────

def clean_md_text(text):
    """Convert markdown bold/italic to reportlab tags."""
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"\*(.+?)\*", r"<i>\1</i>", text)
    text = text.replace("\u2014", " - ").replace("\u2013", " - ")
    return text


def render_section(section, styles):
    """Convert a parsed section into reportlab flowables (English)."""
    elements = []
    heading = section["heading"]
    body = section["body"]

    num_match = re.match(r"(\d+)\.\s+(.*)", heading)
    if num_match:
        num, title = num_match.groups()
        h_text = f'<font color="{MID_BLUE.hexval()}">{num}.</font> {title}'
    else:
        h_text = heading

    elements.append(Paragraph(h_text, styles["heading1"]))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=LIGHT_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    parts = extract_tables_and_text(body)

    for kind, content in parts:
        if kind == "table":
            headers, rows = parse_table(content)
            if headers and rows:
                max_rows = 20
                truncated = False
                if len(rows) > max_rows:
                    rows = rows[:max_rows]
                    truncated = True

                t = build_clinical_table(headers, rows, styles)
                elements.append(t)
                if truncated:
                    elements.append(Paragraph(
                        f"<i>(showing top {max_rows}; "
                        f"see complete data in TSV files)</i>",
                        styles["small"]
                    ))
                elements.append(Spacer(1, 3 * mm))

        elif kind == "text":
            lines = content.split("\n")
            for line in lines:
                stripped = line.strip()
                if not stripped:
                    continue

                if stripped.startswith("### "):
                    sub = stripped[4:].replace("**", "")
                    elements.append(Paragraph(sub, styles["heading2"]))

                elif stripped.startswith("- "):
                    bullet_text = clean_md_text(stripped[2:])
                    elements.append(Paragraph(
                        f"\u2022  {bullet_text}", styles["bullet"]
                    ))

                elif "Disclaimer" in stripped or "not a clinical" in stripped.lower():
                    en_text = (
                        "This report is generated for research and educational "
                        "purposes only. It is not a clinical diagnostic report "
                        "and should not be used for making medical decisions "
                        "without consulting a qualified healthcare professional. "
                        "Variant classifications reflect ClinVar submissions at "
                        "the time of annotation and may change as evidence "
                        "accumulates."
                    )
                    elements.append(Paragraph(en_text, styles["disclaimer"]))

                elif stripped.startswith("*Report prepared"):
                    elements.append(Spacer(1, 2 * mm))
                    sig = stripped.strip("*")
                    elements.append(Paragraph(
                        clean_md_text(sig), styles["small"]
                    ))

                elif stripped.startswith("**") and stripped.endswith("**"):
                    elements.append(Paragraph(
                        clean_md_text(stripped), styles["body_bold"]
                    ))

                else:
                    cleaned = clean_md_text(stripped)
                    if cleaned and not stripped.startswith("|") and not stripped.startswith("---"):
                        elements.append(Paragraph(cleaned, styles["body"]))

    return elements


# ── Build PDF (English) ─────────────────────────────────────────────────

def build_sample_pdf_en(md_path, output_path):
    """Generate an English clinical PDF from a WES markdown report."""
    _require_reportlab()
    report = parse_markdown_report(md_path)
    sample_id = report["sample_id"] or md_path.stem.split("_")[0]
    styles = build_styles()

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=75 * mm,
        bottomMargin=18 * mm,
        title=f"WES Clinical Report - {sample_id}",
        author="M. Corpas - ClawBio WES Analysis Pipeline",
        subject="Whole Exome Sequencing Clinical Report",
    )
    doc._sample_id = sample_id

    frame_cover = Frame(
        MARGIN, 18 * mm,
        PAGE_W - 2 * MARGIN, PAGE_H - 75 * mm - 18 * mm,
        id="cover_frame"
    )
    frame_body = Frame(
        MARGIN, 18 * mm,
        PAGE_W - 2 * MARGIN, PAGE_H - 28 * mm - 18 * mm,
        id="body_frame"
    )

    doc.addPageTemplates([
        PageTemplate(id="cover", frames=frame_cover, onPage=first_page_header_en),
        PageTemplate(id="body", frames=frame_body, onPage=header_footer_en),
    ])

    elements = []
    elements.append(NextPageTemplate("body"))

    # KPI cards
    if report["sections"]:
        kpis = extract_exome_kpis(report["sections"][0]["body"])
        if kpis:
            elements.append(Spacer(1, 2 * mm))
            elements.append(build_kpi_row(kpis, styles))
            elements.append(Spacer(1, 4 * mm))

    # Interpretive paragraph
    interp = build_interpretation_paragraph(report, styles)
    elements.extend(interp)

    # Render all sections
    for section in report["sections"]:
        section_elements = render_section(section, styles)
        elements.extend(section_elements)
        elements.append(Spacer(1, 2 * mm))

    # Section 8: Ancestry
    elements.extend(_build_ancestry_section(sample_id, styles))

    # Section 9: Limitations
    elements.extend(_build_limitations_section(styles))

    # Final disclaimer
    elements.append(Spacer(1, 6 * mm))
    elements.append(HRFlowable(width="100%", thickness=1, color=ACCENT_RED,
                               spaceAfter=3 * mm))

    disclaimer_en = (
        "IMPORTANT NOTICE: This report is generated for research and "
        "educational purposes only. It is not a clinical diagnostic report, "
        "has not been validated for clinical use, and should not be used for "
        "making medical decisions without consulting a qualified healthcare "
        "professional and performing appropriate confirmatory testing. Variant "
        "classifications reflect ClinVar submissions at the time of annotation "
        "and may change as evidence accumulates. Pharmacogenomic predictions "
        "require clinical validation before modifying drug prescriptions."
    )
    elements.append(Paragraph(disclaimer_en, styles["disclaimer"]))

    now = datetime.now()
    date_str = now.strftime("%d %B %Y")
    elements.append(Spacer(1, 4 * mm))
    elements.append(Paragraph(
        f"Report prepared by M. Corpas using ClawBio WES Analysis "
        f"Pipeline v1.0 on {date_str}.",
        styles["small"]
    ))
    elements.append(Paragraph(
        "University of Westminster, London, United Kingdom  |  github.com/ClawBio",
        styles["small"]
    ))

    doc.build(elements)
    return output_path


# ── CLI ─────────────────────────────────────────────────────────────────

def _build_argparser():
    import argparse
    parser = argparse.ArgumentParser(
        description="WES Clinical Report Generator - English (ClawBio Skill)"
    )
    parser.add_argument("--report-dir", type=str, default=str(REPORT_DIR),
                        help="Directory containing WES markdown reports")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_DIR),
                        help="Output directory for PDF reports")
    parser.add_argument("--logo-left", type=str, default=LOGO_LEFT,
                        help="Left institutional logo (JPG/PNG)")
    parser.add_argument("--logo-right", type=str, default=LOGO_RIGHT,
                        help="Right institutional logo (JPG/PNG)")
    parser.add_argument("--samples", type=str, default=None,
                        help="Comma-separated list of samples (default: Sample1-7)")
    parser.add_argument("--demo", action="store_true",
                        help="Run with default Novogene project data")
    return parser


def main():
    import sys

    parser = _build_argparser()
    args = parser.parse_args()

    _mod = sys.modules[__name__]
    report_dir = Path(args.report_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _mod.LOGO_LEFT = args.logo_left
    _mod.LOGO_RIGHT = args.logo_right

    if args.samples:
        samples = [s.strip() for s in args.samples.split(",")]
    else:
        samples = SAMPLES

    if args.demo:
        report_dir = REPORT_DIR
        output_dir = OUTPUT_DIR
        samples = SAMPLES

    print("WES Clinical Report Generator (English)")
    print(f"Input directory:  {report_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Samples: {', '.join(samples)}")
    print(f"{'=' * 60}")

    generated = []
    for sample in samples:
        md_file = report_dir / f"{sample}_WES_Report.md"
        if not md_file.exists():
            print(f"  SKIPPED  {sample} - markdown report not found")
            continue

        pdf_file = output_dir / f"{sample}_WES_Clinical_Report.pdf"
        print(f"  Generating {sample}...", end="", flush=True)
        try:
            build_sample_pdf_en(md_file, pdf_file)
            size_kb = pdf_file.stat().st_size / 1024
            print(f"  OK ({size_kb:.0f} KB)")
            generated.append(pdf_file)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'=' * 60}")
    print(f"Generated {len(generated)}/{len(samples)} PDF reports")
    for f in generated:
        print(f"  {f}")


if __name__ == "__main__":
    main()
