"""Render the prs-abstain markdown reports to PDF.

Kept separate from prs_abstain.py so the skill still runs, and still writes every
text artefact, when reportlab is not installed. Import failures degrade to "no
PDF" rather than "no report".

Note on glyphs: ReportLab's built-in fonts have no Unicode subscript/superscript,
so those render as black boxes. Formulae are written with ASCII markers and
converted to ReportLab's <super>/<sub> XML tags here.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

BRAND = "#1F4E5F"
ACCENT = "#C1443C"
OK = "#2E8B7A"
MUTED = "#6B6B6B"


def _available() -> bool:
    try:
        import reportlab  # noqa: F401
        return True
    except Exception:
        return False


def _inline(text: str) -> str:
    """Markdown inline -> ReportLab mini-XML, with escaping done first."""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", text)             # images handled separately
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)        # links -> plain text
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*(?!\s)(.+?)(?<!\s)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`([^`]+)`", r'<font face="Courier">\1</font>', text)
    text = re.sub(r"\^(\d+)", r"<super>\1</super>", text)          # r^2 -> superscript
    text = text.replace("|", "&#124;")
    return text


def _split_row(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def markdown_to_pdf(md_path: Path, pdf_path: Path, title: str,
                    subtitle: str = "", figures_dir: Path | None = None) -> bool:
    """Render one markdown report to PDF. Returns False if reportlab is missing."""
    if not _available():
        return False

    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph,
                                    SimpleDocTemplate, Spacer, Table, TableStyle)

    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["BodyText"], fontSize=9.5, leading=14,
                          spaceAfter=6, alignment=TA_LEFT)
    h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontSize=17, leading=21,
                        textColor=colors.HexColor(BRAND), spaceBefore=4, spaceAfter=10)
    h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontSize=12.5, leading=16,
                        textColor=colors.HexColor(BRAND), spaceBefore=13, spaceAfter=6)
    h3 = ParagraphStyle("h3", parent=ss["Heading3"], fontSize=10.5, leading=14,
                        textColor=colors.HexColor("#333333"), spaceBefore=10, spaceAfter=4)
    quote = ParagraphStyle("quote", parent=body, leftIndent=8, borderPadding=6,
                           backColor=colors.HexColor("#FBF3F2"),
                           textColor=colors.HexColor(ACCENT), spaceBefore=6, spaceAfter=8)
    code = ParagraphStyle("code", parent=body, fontName="Courier", fontSize=8.2, leading=11,
                          leftIndent=8, backColor=colors.HexColor("#F4F4F4"), spaceAfter=8)
    small = ParagraphStyle("small", parent=body, fontSize=7.8, leading=10,
                           textColor=colors.HexColor(MUTED))
    cell = ParagraphStyle("cell", parent=body, fontSize=7.6, leading=9.6, spaceAfter=0)
    cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")

    story: list[Any] = [Paragraph(_inline(title), h1)]
    if subtitle:
        story.append(Paragraph(_inline(subtitle), small))
    story.append(Spacer(1, 5 * mm))

    lines = md_path.read_text().splitlines()
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i < len(lines) and lines[i].startswith("# "):
        i += 1  # the cover block above already carries this title
    in_code, code_buf = False, []

    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip()

        if line.startswith("```"):
            if in_code:
                story.append(Paragraph("<br/>".join(
                    _inline(c) or "&nbsp;" for c in code_buf), code))
                code_buf, in_code = [], False
            else:
                in_code = True
            i += 1
            continue
        if in_code:
            code_buf.append(line)
            i += 1
            continue

        if not line.strip():
            i += 1
            continue

        # Tables
        if line.lstrip().startswith("|") and i + 1 < len(lines) and \
                set(lines[i + 1].replace("|", "").replace(":", "").strip()) <= {"-", " "} and \
                "-" in lines[i + 1]:
            header = _split_row(line)
            i += 2
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append(_split_row(lines[i]))
                i += 1
            ncol = len(header)
            data = [[Paragraph(_inline(c), cellb) for c in header]]
            for r in rows:
                r = (r + [""] * ncol)[:ncol]
                data.append([Paragraph(_inline(c), cell) for c in r])
            avail = 170 * mm
            widths = [avail / ncol] * ncol
            if ncol >= 3:
                widths[0] = avail * (0.22 if ncol <= 5 else 0.14)
                rest = (avail - widths[0]) / (ncol - 1)
                widths[1:] = [rest] * (ncol - 1)
            t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
            t.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(BRAND)),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#DDDDDD")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                 [colors.white, colors.HexColor("#F7F9FA")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]))
            story.extend([t, Spacer(1, 4 * mm)])
            continue

        # Images
        m = re.match(r"!\[[^\]]*\]\(([^)]+)\)", line.strip())
        if m and figures_dir is not None:
            src = (figures_dir.parent / m.group(1)).resolve()
            if src.exists():
                try:
                    from PIL import Image as PILImage
                    w, h = PILImage.open(src).size
                    disp_w = 165 * mm
                    story.append(KeepTogether([
                        Image(str(src), width=disp_w, height=disp_w * h / w),
                        Spacer(1, 4 * mm)]))
                except Exception:
                    pass
            i += 1
            continue

        if line.startswith("#### "):
            story.append(Paragraph(_inline(line[5:]), h3))
        elif line.startswith("### "):
            story.append(Paragraph(_inline(line[4:]), h3))
        elif line.startswith("## "):
            story.append(Paragraph(_inline(line[3:]), h2))
        elif line.startswith("# "):
            story.append(Paragraph(_inline(line[2:]), h1))
        elif line.startswith("> "):
            story.append(Paragraph(_inline(line[2:]), quote))
        elif line.startswith("---"):
            story.append(Spacer(1, 3 * mm))
        elif re.match(r"^\s*[-*] ", line):
            story.append(Paragraph("• " + _inline(re.sub(r"^\s*[-*] ", "", line)),
                                   ParagraphStyle("li", parent=body, leftIndent=10)))
        elif re.match(r"^\s*\d+\. ", line):
            story.append(Paragraph(_inline(line.strip()),
                                   ParagraphStyle("ol", parent=body, leftIndent=10)))
        else:
            story.append(Paragraph(_inline(line), body))
        i += 1

    def _decorate(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor(BRAND))
        canvas.setLineWidth(2)
        canvas.line(20 * mm, 285 * mm, 190 * mm, 285 * mm)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor(MUTED))
        canvas.drawString(20 * mm, 288 * mm, "ClawBio prs-abstain")
        canvas.drawRightString(190 * mm, 288 * mm,
                               "Research tool. Not a medical device.")
        canvas.drawCentredString(105 * mm, 12 * mm, f"Page {doc.page}")
        canvas.restoreState()

    SimpleDocTemplate(
        str(pdf_path), pagesize=A4,
        leftMargin=20 * mm, rightMargin=20 * mm,
        topMargin=25 * mm, bottomMargin=18 * mm,
        title=title, author="ClawBio prs-abstain",
    ).build(story, onFirstPage=_decorate, onLaterPages=_decorate)
    return True
