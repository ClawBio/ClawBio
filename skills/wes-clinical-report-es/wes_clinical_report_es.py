#!/usr/bin/env python3
"""
Generador de Informes Clinicos PDF en Espanol - Exomas Novogene WES
Analiza los informes WES en markdown y genera PDFs clinicos profesionales
traducidos al espanol.

Proyecto: X202SC26016276-Z01-F001 (Novogene WES, 7 muestras, hg38)
"""

import re
import textwrap
from pathlib import Path
from datetime import datetime

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm, cm
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        PageBreak, HRFlowable, KeepTogether, Frame, PageTemplate,
        BaseDocTemplate, NextPageTemplate, Image
    )
    from reportlab.platypus.doctemplate import _doNothing
    from reportlab.graphics.shapes import Drawing, Rect, String, Line
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    _REPORTLAB_AVAILABLE = True
except ImportError:
    colors = None
    A4 = None
    getSampleStyleSheet = None
    ParagraphStyle = None
    mm = None
    cm = None
    TA_LEFT = None
    TA_CENTER = None
    TA_RIGHT = None
    TA_JUSTIFY = None
    SimpleDocTemplate = None
    Paragraph = None
    Spacer = None
    Table = None
    TableStyle = None
    PageBreak = None
    HRFlowable = None
    KeepTogether = None
    Frame = None
    PageTemplate = None
    BaseDocTemplate = None
    NextPageTemplate = None
    Image = None
    _doNothing = None
    Drawing = None
    Rect = None
    String = None
    Line = None
    pdfmetrics = None
    TTFont = None
    _REPORTLAB_AVAILABLE = False

# ── Paths ────────────────────────────────────────────────────────────────
REPORT_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/REPORTS")
OUTPUT_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/REPORTS/PDF-ES")

LOGO_PREDICE = str(REPORT_DIR / "logo_predice.jpg")
LOGO_INBIOMEDIC = str(REPORT_DIR / "logo_inbiomedic.jpg")
ANCESTRY_DIR = Path("/Volumes/CPM-16Tb/NOVOGENE/ANALYSIS/ANCESTRY/OUTPUT")

SAMPLES = [f"Sample{i}" for i in range(1, 8)]

# ── Spanish month names ──────────────────────────────────────────────────
MESES_ES = {
    1: "enero", 2: "febrero", 3: "marzo", 4: "abril",
    5: "mayo", 6: "junio", 7: "julio", 8: "agosto",
    9: "septiembre", 10: "octubre", 11: "noviembre", 12: "diciembre",
}

def fecha_es():
    now = datetime.now()
    return f"{now.day} de {MESES_ES[now.month]} de {now.year}"

# ── Translation dictionaries ─────────────────────────────────────────────

# Section headings
HEADING_ES = {
    "Exome Summary": "Resumen del Exoma",
    "Clinically Significant Variants": "Variantes de Significancia Clinica",
    "Pharmacogenomics": "Farmacogenomica",
    "Fitness and Nutrition Traits": "Rasgos de Aptitud Fisica y Nutricion",
    "Prioritised Rare Damaging Variants": "Variantes Raras Patogenicas Priorizadas",
    "Disease and Pathway Context": "Contexto de Enfermedad y Vias Metabolicas",
    "Methods": "Metodologia",
}

# Table headers
HEADER_ES = {
    "Gene": "Gen",
    "Variant": "Variante",
    "Zygosity": "Cigosidad",
    "Classification": "Clasificacion",
    "Consequence": "Consecuencia",
    "Associated Condition": "Condicion Asociada",
    "Metric": "Metrica",
    "Count": "Cantidad",
    "Allele": "Alelo",
    "Clinical Effect": "Efecto Clinico",
    "Affected Medications": "Medicamentos Afectados",
    "Trait": "Rasgo",
    "Interpretation": "Interpretacion",
    "Ev.": "Ev.",
    "REVEL": "REVEL",
    "CADD": "CADD",
    "gnomAD AF": "gnomAD AF",
    "OMIM Disease": "Enfermedad OMIM",
}

# Table cell values / metric labels
METRIC_ES = {
    "Total SNP variants": "Variantes SNP totales",
    "Missense": "Missense (cambio de sentido)",
    "Synonymous": "Sinonimas",
    "Stopgain": "Ganancia de codon de parada",
    "Frameshift": "Cambio de marco de lectura",
    "Splicing": "Variantes de empalme",
    "Loss-of-function (stopgain + frameshift)": "Perdida de funcion (stopgain + frameshift)",
    "Rare coding (gnomAD < 1%)": "Codificantes raras (gnomAD < 1%)",
    "Rare + computationally damaging": "Raras + computacionalmente patogenicas",
    "ClinVar Pathogenic / Likely Pathogenic": "ClinVar Patogenica / Prob. Patogenica",
}

# Zygosity
ZYGOSITY_ES = {
    "Hom": "Hom",
    "Het": "Het",
}

# Classification
CLASSIF_ES = {
    "Pathogenic": "Patogenica",
    "Pathogenic/LP": "Patogenica/PP",
    "Likely Pathogenic": "Probablemente Patogenica",
    "Benign": "Benigna",
    "Likely Benign": "Probablemente Benigna",
    "VUS": "VUS",
}

# Consequence
CONSEQ_ES = {
    "missense SNV": "SNV de cambio de sentido",
    "intronic": "intronica",
    "synonymous SNV": "SNV sinonima",
    "stopgain": "ganancia de parada",
    "frameshift deletion": "delecion con cambio de marco",
    "frameshift insertion": "insercion con cambio de marco",
    "splicing": "empalme",
    "UTR3": "UTR3",
    "UTR5": "UTR5",
}

# Clinical effects (PGx)
EFFECT_ES = {
    "Non-expressor": "No expresor",
    "Decreased transport": "Transporte reducido",
    "Slow metaboliser": "Metabolizador lento",
    "Slow acetylator": "Acetilador lento",
    "Normal/Increased function": "Funcion normal/aumentada",
    "Variable": "Variable",
    "Low-dose warfarin": "Warfarina dosis baja",
    "Increased transport": "Transporte aumentado",
    "Ultra-rapid (inducible)": "Ultra-rapido (inducible)",
    "Decreased function": "Funcion reducida",
    "Decreased transport": "Transporte reducido",
}

# Body text patterns to translate
BODY_ES = {
    "Key implications:": "Implicaciones clinicas principales:",
    "variant(s) classified as Pathogenic or Likely Pathogenic in ClinVar:":
        "variante(s) clasificadas como Patogenicas o Probablemente Patogenicas en ClinVar:",
    "An additional": "Adicionalmente,",
    "variant(s) have conflicting or uncertain classifications.":
        "variante(s) tienen clasificaciones conflictivas o inciertas.",
    "The most clinically relevant are listed below (coding variants in disease-associated genes):":
        "Las mas relevantes clinicamente se enumeran a continuacion (variantes codificantes en genes asociados a enfermedad):",
    "ACMG SF v3.2 actionable genes": "Genes accionables ACMG SF v3.2",
    "coding variants identified across": "variantes codificantes identificadas en",
    "medically actionable genes, of which": "genes medicamente accionables, de los cuales",
    "have ClinVar P/LP classification.": "tienen clasificacion P/PP en ClinVar.",
    "Cancer predisposition panel": "Panel de predisposicion al cancer",
    "with P/LP classification across": "con clasificacion P/PP en",
    "cancer predisposition genes.": "genes de predisposicion al cancer.",
    "The following pharmacogenomic markers were identified from CPIC-defined star-allele positions. Variants are reported where the genotype differs from the reference allele.":
        "Los siguientes marcadores farmacogenomicos se identificaron a partir de posiciones de alelos estrella definidas por CPIC. Se informan las variantes donde el genotipo difiere del alelo de referencia.",
    "Genotypes at positions associated with fitness and nutrition traits":
        "Genotipos en posiciones asociadas con rasgos de aptitud fisica y nutricion",
    "Corpas et al. 2021, Tables 3": "Corpas et al. 2021, Tablas 3",
    "Only markers captured by the WES panel and with non-reference genotypes are shown.":
        "Solo se muestran los marcadores capturados por el panel WES con genotipos no referencia.",
    "Evidence grades: A = strong replication, B = moderate, C = preliminary.":
        "Grados de evidencia: A = replicacion fuerte, B = moderada, C = preliminar.",
    "variants pass all filters: coding, rare (gnomAD AF < 0.01), and computationally predicted damaging (CADD > 20 or REVEL > 0.5). Top 15 ranked by pathogenicity prediction score:":
        "variantes pasan todos los filtros: codificantes, raras (gnomAD AF < 0.01), y predichas computacionalmente como patogenicas (CADD > 20 o REVEL > 0.5). Las 15 principales ordenadas por puntuacion de prediccion de patogenicidad:",
    "Across the full variant set:": "En el conjunto completo de variantes:",
    "variants map to OMIM disease entries,": "variantes mapean a entradas de enfermedad OMIM,",
    "overlap GWAS Catalog associations, and": "se solapan con asociaciones del Catalogo GWAS, y",
    "have COSMIC somatic mutation records.": "tienen registros de mutaciones somaticas en COSMIC.",
    "KEGG pathways enriched in rare coding variants:":
        "Vias KEGG enriquecidas en variantes codificantes raras:",
    "Whole exome sequencing was performed on an Illumina NovaSeq 6000 platform using 150 bp paired-end reads with the Xplus capture kit (60.5 Mb target region).":
        "La secuenciacion del exoma completo se realizo en una plataforma Illumina NovaSeq 6000 utilizando lecturas pareadas de 150 pb con el kit de captura Xplus (region objetivo de 60,5 Mb).",
    "Reads were aligned to the GRCh38/hg38 reference genome using BWA-MEM.":
        "Las lecturas se alinearon al genoma de referencia GRCh38/hg38 utilizando BWA-MEM.",
    "Variant calling was performed with GATK HaplotypeCaller v4.3.0 following GATK Best Practices.":
        "La identificacion de variantes se realizo con GATK HaplotypeCaller v4.3.0 siguiendo las Mejores Practicas de GATK.",
    "Functional annotation was performed with ANNOVAR, incorporating ClinVar (2024), gnomAD v3.1.2 (9 population groups), COSMIC, OMIM, SIFT, PolyPhen-2, CADD, REVEL, and 15 additional databases.":
        "La anotacion funcional se realizo con ANNOVAR, incorporando ClinVar (2024), gnomAD v3.1.2 (9 grupos poblacionales), COSMIC, OMIM, SIFT, PolyPhen-2, CADD, REVEL y 15 bases de datos adicionales.",
    "Pharmacogenomic analysis used CPIC star-allele definitions with evidence enrichment from the ClinPGx API (PharmGKB).":
        "El analisis farmacogenomico utilizo definiciones de alelos estrella de CPIC con enriquecimiento de evidencia de la API ClinPGx (PharmGKB).",
    "Fitness and nutrition trait interpretation followed the evidence framework of Corpas et al. (2021)":
        "La interpretacion de rasgos de aptitud fisica y nutricion siguio el marco de evidencia de Corpas et al. (2021)",
    "Frontiers in Genetics 12:535123.": "Frontiers in Genetics 12:535123.",
    "Variant prioritisation applied sequential filters: coding consequence":
        "La priorizacion de variantes aplico filtros secuenciales: consecuencia codificante",
    "population frequency (gnomAD AF < 0.01)": "frecuencia poblacional (gnomAD AF < 0,01)",
    "computational pathogenicity (CADD > 20 or REVEL > 0.5).":
        "patogenicidad computacional (CADD > 20 o REVEL > 0,5).",
}

# Fitness/nutrition trait translations
TRAIT_ES = {
    "Muscle fibre type (power vs endurance)": "Tipo de fibra muscular (potencia vs resistencia)",
    "Fat metabolism during exercise": "Metabolismo de grasas durante ejercicio",
    "Folate metabolism (C677T)": "Metabolismo del folato (C677T)",
    "Vitamin D binding protein": "Proteina transportadora de vitamina D",
    "Vitamin D synthesis": "Sintesis de vitamina D",
    "Omega-3 conversion": "Conversion de Omega-3",
    "Fat sensitivity & T2D risk": "Sensibilidad a grasas y riesgo de DT2",
    "Dietary fat absorption": "Absorcion de grasa dietaria",
    "Type 2 diabetes risk": "Riesgo de diabetes tipo 2",
    "Salt sensitivity / hypertension": "Sensibilidad a la sal / hipertension",
    "Salt sensitivity": "Sensibilidad a la sal",
    "Alcohol metabolism speed": "Velocidad del metabolismo del alcohol",
    "Bitter taste perception": "Percepcion del sabor amargo",
    "Bitter taste (linked)": "Sabor amargo (asociado)",
    "Bitter taste (3rd SNP)": "Sabor amargo (3er SNP)",
    "Beta-carotene conversion": "Conversion de beta-caroteno",
    "Vitamin E & triglycerides": "Vitamina E y trigliceridos",
    "Caffeine metabolism": "Metabolismo de la cafeina",
    "Triglyceride response to fat": "Respuesta de trigliceridos a la grasa",
}

INTERP_ES = {
    "XX — endurance phenotype": "XX: fenotipo de resistencia",
    "XX - endurance phenotype": "XX: fenotipo de resistencia",
    "Moderate effect": "Efecto moderado",
    "Altered fat metabolism": "Metabolismo de grasas alterado",
    "CT — 35% reduced": "CT: 35% reducido",
    "CT - 35% reduced": "CT: 35% reducido",
    "Lower vitamin D": "Vitamina D reducida",
    "Reduced synthesis": "Sintesis reducida",
    "Poor converter": "Conversor deficiente",
    "Moderate benefit": "Beneficio moderado",
    "Moderate": "Moderado",
    "Increased absorption": "Absorcion aumentada",
    "Salt-sensitive": "Sensible a la sal",
    "Ultra-rapid": "Ultra-rapido",
    "Medium taster": "Percepcion media",
    "Medium": "Media",
    "Reduced": "Reducida",
    "Higher triglycerides": "Trigliceridos elevados",
    "Slow metaboliser": "Metabolizador lento",
}


if _REPORTLAB_AVAILABLE:
    # ── Colour palette ───────────────────────────────────────────────────
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
            "`pip install -r skills/wes-clinical-report-es/requirements.txt`."
        )


# ── Translation helper ───────────────────────────────────────────────────

def translate_text(text):
    """Apply all body text translations."""
    result = text
    for en, es in BODY_ES.items():
        result = result.replace(en, es)
    return result


def translate_header(header):
    """Translate a single table header."""
    return HEADER_ES.get(header, header)


def translate_metric(label):
    """Translate a metric table row label."""
    return METRIC_ES.get(label, label)


def translate_heading(heading):
    """Translate a section heading, preserving the number prefix."""
    m = re.match(r"(\d+\.\s+)(.*)", heading)
    if m:
        prefix, title = m.groups()
        return prefix + HEADING_ES.get(title.strip(), title)
    return HEADING_ES.get(heading.strip(), heading)


def translate_cell(header, cell):
    """Translate a table cell value based on its column header."""
    h = header.strip()
    c = cell.strip()

    if h in ("Zygosity", "Cigosidad"):
        return ZYGOSITY_ES.get(c, c)
    if h in ("Classification", "Clasificacion"):
        return CLASSIF_ES.get(c, c)
    if h in ("Consequence", "Consecuencia"):
        return CONSEQ_ES.get(c, c)
    if h in ("Clinical Effect", "Efecto Clinico"):
        return EFFECT_ES.get(c, c)
    if h in ("Trait", "Rasgo"):
        return TRAIT_ES.get(c, c)
    if h in ("Interpretation", "Interpretacion"):
        return INTERP_ES.get(c, c)
    if h in ("Metric", "Metrica"):
        return METRIC_ES.get(c, c)
    return cell


# ── Markdown parser (same as English version) ────────────────────────────

def parse_markdown_report(filepath):
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


# ── PDF styles ───────────────────────────────────────────────────────────

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
            textColor=BLACK, leading=13,
            spaceAfter=2 * mm,
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
            textColor=WARM_GREY, leading=10,
            spaceAfter=1 * mm,
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


# ── Page templates (Spanish) ─────────────────────────────────────────────

def header_footer_es(canvas, doc):
    canvas.saveState()

    # ── Small logo strip above navy bar ──────────────────────────────
    logo_bar_h = 12 * mm
    logo_bar_y = PAGE_H - logo_bar_h
    canvas.setFillColor(WHITE)
    canvas.rect(0, logo_bar_y, PAGE_W, logo_bar_h, fill=1, stroke=0)

    # Predice logo (small)
    pred_h = 8 * mm
    pred_w = pred_h * (166 / 100)
    try:
        canvas.drawImage(LOGO_PREDICE, MARGIN, logo_bar_y + 2 * mm,
                         pred_w, pred_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    # Inbiomedic logo (small)
    inb_h = 7 * mm
    inb_w = inb_h * (431 / 133)
    try:
        canvas.drawImage(LOGO_INBIOMEDIC,
                         PAGE_W - MARGIN - inb_w, logo_bar_y + 2.5 * mm,
                         inb_w, inb_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    # ── Navy header bar ──────────────────────────────────────────────
    navy_bar_y = logo_bar_y - 11 * mm
    canvas.setFillColor(NAVY)
    canvas.rect(0, navy_bar_y, PAGE_W, 11 * mm, fill=1, stroke=0)

    canvas.setFillColor(WHITE)
    canvas.setFont("Helvetica-Bold", 8)
    canvas.drawString(MARGIN, navy_bar_y + 3.5 * mm,
                      "INFORME DE SECUENCIACION DEL EXOMA COMPLETO")

    canvas.setFont("Helvetica", 8)
    sample_id = getattr(doc, '_sample_id', '')
    canvas.drawRightString(PAGE_W - MARGIN, navy_bar_y + 3.5 * mm,
                           sample_id)

    # Accent line below navy bar
    canvas.setStrokeColor(MID_BLUE)
    canvas.setLineWidth(1.5)
    canvas.line(0, navy_bar_y - 0.5 * mm, PAGE_W, navy_bar_y - 0.5 * mm)

    # ── Footer ───────────────────────────────────────────────────────
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(WARM_GREY)

    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 12 * mm, PAGE_W - MARGIN, 12 * mm)

    canvas.drawString(MARGIN, 8 * mm, "ClawBio WES Analysis Pipeline v1.0")
    canvas.drawCentredString(PAGE_W / 2, 8 * mm,
                             "Proyecto X202SC26016276-Z01-F001")
    canvas.drawRightString(PAGE_W - MARGIN, 8 * mm,
                           f"Pagina {doc.page}")

    canvas.setFont("Helvetica-Oblique", 6)
    canvas.setFillColor(colors.HexColor("#BDC3C7"))
    canvas.drawCentredString(PAGE_W / 2, 4 * mm,
                             "CONFIDENCIAL - SOLO PARA USO EN INVESTIGACION"
                             " - NO PARA DIAGNOSTICO CLINICO")

    canvas.restoreState()


def _draw_logos(canvas, y_center):
    """Draw Predice (left) and Inbiomedic (right) logos centred on y_center."""
    from reportlab.lib.utils import ImageReader

    # Predice logo - left side
    # Original: 166x100 px. Render at 20mm height.
    pred_h = 16 * mm
    pred_w = pred_h * (166 / 100)
    pred_x = MARGIN
    pred_y = y_center - pred_h / 2
    try:
        canvas.drawImage(LOGO_PREDICE, pred_x, pred_y, pred_w, pred_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass

    # Inbiomedic logo - right side
    # Original: 431x133 px. Render at 16mm height.
    inb_h = 14 * mm
    inb_w = inb_h * (431 / 133)
    inb_x = PAGE_W - MARGIN - inb_w
    inb_y = y_center - inb_h / 2
    try:
        canvas.drawImage(LOGO_INBIOMEDIC, inb_x, inb_y, inb_w, inb_h,
                         preserveAspectRatio=True, mask='auto')
    except Exception:
        pass


def first_page_header_es(canvas, doc):
    canvas.saveState()

    # ── Logo strip (white background above the navy header) ──────────
    logo_strip_h = 22 * mm
    logo_strip_y = PAGE_H - logo_strip_h
    canvas.setFillColor(WHITE)
    canvas.rect(0, logo_strip_y, PAGE_W, logo_strip_h, fill=1, stroke=0)
    _draw_logos(canvas, logo_strip_y + logo_strip_h / 2)

    # Thin separator below logos
    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.8)
    canvas.line(MARGIN, logo_strip_y, PAGE_W - MARGIN, logo_strip_y)

    # ── Navy header block (shifted down by logo strip) ───────────────
    header_h = 48 * mm
    header_y = logo_strip_y - header_h
    canvas.setFillColor(NAVY)
    canvas.rect(0, header_y, PAGE_W, header_h, fill=1, stroke=0)

    # Accent stripe
    canvas.setFillColor(MID_BLUE)
    canvas.rect(0, header_y - 2 * mm, PAGE_W, 2 * mm, fill=1, stroke=0)

    # Title text
    canvas.setFillColor(WHITE)
    canvas.setFont("Helvetica-Bold", 24)
    canvas.drawString(MARGIN, header_y + header_h - 16 * mm,
                      "Secuenciacion del")
    canvas.drawString(MARGIN, header_y + header_h - 26 * mm,
                      "Exoma Completo")

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

    # Metadata lines
    canvas.setFont("Helvetica", 9)
    canvas.setFillColor(colors.HexColor("#A0B4CC"))
    canvas.drawString(MARGIN, header_y + 9 * mm,
                      "Illumina NovaSeq PE150  |  Xplus WES (60,5 Mb)  |  "
                      "GRCh38/hg38  |  GATK HaplotypeCaller 4.3.0")
    canvas.drawString(MARGIN, header_y + 3 * mm,
                      f"Informe generado: {fecha_es()}  |  "
                      "Analisis: ClawBio + ANNOVAR + ClinPGx")

    # ── Footer ───────────────────────────────────────────────────────
    canvas.setStrokeColor(LIGHT_BLUE)
    canvas.setLineWidth(0.5)
    canvas.line(MARGIN, 12 * mm, PAGE_W - MARGIN, 12 * mm)

    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(WARM_GREY)
    canvas.drawString(MARGIN, 8 * mm, "ClawBio WES Analysis Pipeline v1.0")
    canvas.drawCentredString(PAGE_W / 2, 8 * mm,
                             "Proyecto X202SC26016276-Z01-F001")
    canvas.drawRightString(PAGE_W - MARGIN, 8 * mm, f"Pagina {doc.page}")

    canvas.setFont("Helvetica-Oblique", 6)
    canvas.setFillColor(colors.HexColor("#BDC3C7"))
    canvas.drawCentredString(PAGE_W / 2, 4 * mm,
                             "CONFIDENCIAL - SOLO PARA USO EN INVESTIGACION"
                             " - NO PARA DIAGNOSTICO CLINICO")

    canvas.restoreState()


# ── KPI cards ────────────────────────────────────────────────────────────

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


# ── Table builder ────────────────────────────────────────────────────────

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


# ── Section rendering (Spanish) ──────────────────────────────────────────

def clean_md_text(text):
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"\*(.+?)\*", r"<i>\1</i>", text)
    text = text.replace("—", " - ").replace("–", " - ")
    return text


# Sub-heading translations
SUBHEADING_ES = {
    "Fitness": "Aptitud Fisica",
    "Nutrition": "Nutricion",
}


def render_section_es(section, styles):
    """Convert a parsed section into Spanish reportlab flowables."""
    elements = []
    heading = translate_heading(section["heading"])
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
            headers_en, rows = parse_table(content)
            if headers_en and rows:
                # Translate headers
                headers_es = [translate_header(h) for h in headers_en]

                # Translate cell values based on column
                rows_es = []
                for row in rows:
                    translated_row = []
                    for j, cell in enumerate(row):
                        if j < len(headers_en):
                            translated_row.append(translate_cell(headers_en[j], cell))
                        else:
                            translated_row.append(cell)
                    rows_es.append(translated_row)

                max_rows = 20
                truncated = False
                if len(rows_es) > max_rows:
                    rows_es = rows_es[:max_rows]
                    truncated = True

                t = build_clinical_table(headers_es, rows_es, styles)
                elements.append(t)
                if truncated:
                    elements.append(Paragraph(
                        f"<i>(mostrando las {max_rows} principales; "
                        f"ver datos completos en archivos TSV)</i>",
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
                    sub_es = SUBHEADING_ES.get(sub, sub)
                    elements.append(Paragraph(sub_es, styles["heading2"]))

                elif stripped.startswith("- "):
                    bullet_text = clean_md_text(translate_text(stripped[2:]))
                    elements.append(Paragraph(
                        f"\u2022  {bullet_text}", styles["bullet"]
                    ))

                elif "Disclaimer" in stripped or "not a clinical" in stripped.lower():
                    es_text = (
                        "Este informe se genera exclusivamente con fines de "
                        "investigacion y educacion. No es un informe de "
                        "diagnostico clinico y no debe utilizarse para tomar "
                        "decisiones medicas sin consultar a un profesional "
                        "sanitario cualificado. Las clasificaciones de variantes "
                        "reflejan las presentaciones de ClinVar en el momento de "
                        "la anotacion y pueden cambiar conforme se acumule evidencia."
                    )
                    elements.append(Paragraph(es_text, styles["disclaimer"]))

                elif stripped.startswith("*Report prepared"):
                    elements.append(Spacer(1, 2 * mm))
                    sig = stripped.strip("*").replace(
                        "Report prepared by", "Informe preparado por"
                    ).replace("using", "utilizando").replace("on ", "el ")
                    elements.append(Paragraph(
                        clean_md_text(sig), styles["small"]
                    ))

                elif stripped.startswith("**") and stripped.endswith("**"):
                    elements.append(Paragraph(
                        clean_md_text(translate_text(stripped)),
                        styles["body_bold"]
                    ))

                else:
                    cleaned = clean_md_text(translate_text(stripped))
                    if cleaned and not stripped.startswith("|") and not stripped.startswith("---"):
                        elements.append(Paragraph(cleaned, styles["body"]))

    return elements


# ── Extract KPI metrics ──────────────────────────────────────────────────

KPI_LABELS_ES = {
    "Total SNPs": "SNPs Totales",
    "Missense": "Missense",
    "Stopgain": "Stopgain",
    "Rare Damaging": "Raras Patog.",
    "ClinVar P/LP": "ClinVar P/PP",
}

def extract_exome_kpis_es(section_body):
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
                kpis["SNPs Totales"] = f"{int(value):,}"
            elif "Missense" == label:
                kpis["Missense"] = f"{int(value):,}"
            elif "ClinVar" in label:
                kpis["ClinVar P/PP"] = value
            elif "Rare + comp" in label or "Rare + comput" in label:
                kpis["Raras Patog."] = value
            elif "Stopgain" in label:
                kpis["Stopgain"] = value

    return kpis


# ── Clinical interpretation summary builder ──────────────────────────────

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
                    for row in rows[:5]:  # top 5
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

            # Extract Hom/Het ratio from text
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


def build_interpretation_paragraph(report, styles):
    """Build a concise interpretive paragraph explaining the significant findings."""
    elements = []

    pathogenic = _extract_pathogenic_variants(report)
    pgx = _extract_pgx_alerts(report)
    rare_dmg = _extract_rare_damaging(report)
    metrics = _extract_metrics(report)
    acmg = _extract_acmg_cancer(report)

    hom_het = metrics.get("hom_het_ratio", 0)
    rare_count = metrics.get("Rare + computationally damaging", "0")
    sample_id = report["sample_id"]

    # ── Build the paragraph ─────────────────────────────────────────
    sentences = []

    # Pathogenic variants
    if pathogenic:
        gene_list = ", ".join(v["gene"] for v in pathogenic)
        hom_genes = [v["gene"] for v in pathogenic if v["zygosity"] == "Hom"]
        het_genes = [v["gene"] for v in pathogenic if v["zygosity"] == "Het"]
        sentences.append(
            f"El analisis del exoma de {sample_id} identifico "
            f"{len(pathogenic)} variante(s) clasificadas como patogenicas "
            f"en ClinVar, localizadas en los genes {gene_list}."
        )
        # Describe conditions
        conditions_seen = set()
        for v in pathogenic:
            cond = v["condition"]
            if cond and cond not in conditions_seen:
                conditions_seen.add(cond)
        if conditions_seen:
            cond_text = "; ".join(
                c[:80] for c in sorted(conditions_seen)
            )
            sentences.append(
                f"Estas variantes estan asociadas a las siguientes condiciones: "
                f"{cond_text}."
            )
        if hom_genes:
            sentences.append(
                f"Las variantes en {', '.join(hom_genes)} se encuentran "
                f"en estado homocigoto, lo que indica que ambas copias del "
                f"gen estan afectadas y el efecto clinico puede ser mas "
                f"pronunciado."
            )
    else:
        sentences.append(
            f"El analisis del exoma de {sample_id} no identifico variantes "
            f"clasificadas como patogenicas confirmadas en ClinVar."
        )

    # ACMG / Cancer
    if acmg["acmg_plp"] > 0 or acmg["cancer_plp"] > 0:
        parts = []
        if acmg["acmg_plp"] > 0:
            parts.append(
                f"{acmg['acmg_plp']} variante(s) en genes medicamente "
                f"accionables segun las guias ACMG SF v3.2"
            )
        if acmg["cancer_plp"] > 0:
            parts.append(
                f"{acmg['cancer_plp']} variante(s) en genes de "
                f"predisposicion al cancer"
            )
        sentences.append(
            "Adicionalmente se identificaron " + ", y ".join(parts) +
            ", que pueden requerir seguimiento clinico especializado."
        )

    # PGx interpretation
    if pgx:
        pgx_parts = []
        for p in pgx:
            gene = p["gene"]
            eff = EFFECT_ES.get(p["effect"], p["effect"]).lower()
            meds = p["meds"]
            if meds:
                pgx_parts.append(f"{gene} ({eff}, afecta a {meds})")
            else:
                pgx_parts.append(f"{gene} ({eff})")
        sentences.append(
            f"Desde el punto de vista farmacogenomico, se detectaron "
            f"variantes relevantes en {', '.join(pgx_parts)}. "
            f"Estos hallazgos pueden condicionar la respuesta individual "
            f"a determinados farmacos y deben considerarse antes de "
            f"prescribir los medicamentos indicados."
        )

    # Homozygosity
    if hom_het > 2.0:
        sentences.append(
            f"Se observa un patron de homocigosidad elevada "
            f"(ratio Hom/Het = {hom_het:.2f}, frente al ~1,5 esperado en "
            f"poblaciones no consanguineas), lo que puede reflejar "
            f"estructura poblacional especifica o consanguinidad parental, "
            f"y debe tenerse en cuenta al interpretar el riesgo de "
            f"enfermedades recesivas."
        )

    # Rare damaging
    if rare_count and int(rare_count) > 0:
        top_genes = ", ".join(v["gene"] for v in rare_dmg[:3]) if rare_dmg else ""
        sent = (
            f"Se priorizaron {rare_count} variantes raras con prediccion "
            f"computacional de patogenicidad (REVEL > 0,5 o CADD > 20)"
        )
        if top_genes:
            sent += f", entre las que destacan variantes en {top_genes}"
        sent += (
            ". Estas variantes requieren correlacion con el fenotipo "
            "clinico del paciente para determinar su relevancia diagnostica."
        )
        sentences.append(sent)

    # Closing
    sentences.append(
        "Se recomienda validar las variantes de mayor relevancia clinica "
        "mediante secuenciacion Sanger u otra metodologia independiente, "
        "y derivar al paciente a un servicio de genetica clinica para "
        "asesoramiento y evaluacion familiar si procede."
    )

    # ── Render as styled paragraph block ─────────────────────────────
    elements.append(Paragraph(
        "<b>Interpretacion de resultados</b>",
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


# ── Ancestry section ─────────────────────────────────────────────────────

def _build_ancestry_section(sample_id, styles):
    """Build ancestry section showing admixture percentages for this sample only."""
    import json

    elements = []
    elements.append(Paragraph(
        f'<font color="{MID_BLUE.hexval()}">8.</font> '
        f'Estimacion de Ancestria',
        styles["heading1"]
    ))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=LIGHT_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    # Load ancestry results
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

        # Build interpretation paragraph
        top_pop, top_pct = sorted_pops[0]
        components = [f"{pop} ({pct:.1f}%)" for pop, pct in sorted_pops if pct >= 0.5]

        interp = (
            f"La composicion ancestral de {sample_id} se estimo mediante "
            f"descomposicion por minimos cuadrados no negativos (NNLS) de las "
            f"frecuencias alelicas observadas frente a las poblaciones de "
            f"referencia continentales del proyecto 1000 Genomes, utilizando "
            f"{n_markers:,} marcadores informativos de ancestria (AIMs). "
            f"El componente predominante es <b>{top_pop} ({top_pct:.1f}%)</b>, "
            f"seguido de {', '.join(components[1:])}."
        )

        if "Americano" in top_pop:
            interp += (
                " Este perfil es consistente con una ancestria latinoamericana "
                "admixta, con un componente principal amerindio y contribuciones "
                "menores de origen asiatico oriental y africano, como es tipico "
                "de poblaciones de America Latina."
            )

        elements.append(Paragraph(interp, styles["body"]))
        elements.append(Spacer(1, 3 * mm))

        # Admixture table
        headers = ["Categoria Ancestral", "Proporcion (%)"]
        rows = []
        for pop, pct in sorted_pops:
            rows.append([pop, f"{pct:.1f}%"])
        t = build_clinical_table(headers, rows, styles)
        elements.append(t)
        elements.append(Spacer(1, 2 * mm))

        # Embed per-sample pie chart
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
            f"<i>Estimacion basada en {n_markers:,} marcadores informativos de "
            f"ancestria (AIMs) con frecuencias poblacionales conocidas en el "
            f"proyecto 1000 Genomes. Metodo: minimos cuadrados no negativos (NNLS).</i>",
            styles["small"]
        ))
    else:
        elements.append(Paragraph(
            "No se dispone de datos de ancestria para esta muestra.",
            styles["body"]
        ))

    elements.append(Spacer(1, 2 * mm))
    return elements


# ── Limitations section ──────────────────────────────────────────────────

def _build_limitations_section(styles):
    """Build limitations section in Spanish."""
    elements = []
    elements.append(Paragraph(
        f'<font color="{MID_BLUE.hexval()}">9.</font> Limitaciones',
        styles["heading1"]
    ))
    elements.append(HRFlowable(
        width="100%", thickness=0.8, color=LIGHT_BLUE,
        spaceAfter=3 * mm, spaceBefore=0
    ))

    limitations = [
        (
            "Cobertura del exoma",
            "La secuenciacion del exoma completo (WES) captura aproximadamente "
            "el 1-2% del genoma, correspondiente a las regiones codificantes de "
            "proteinas. Variantes en regiones intronicas profundas, intergenicas, "
            "reguladoras o en ADN no codificante no son detectadas por esta "
            "metodologia. Variantes estructurales grandes (deleciones, "
            "duplicaciones, inversiones) y expansiones de repeticiones tampoco "
            "se detectan de forma fiable mediante WES."
        ),
        (
            "Clasificacion de variantes",
            "Las clasificaciones de patogenicidad provienen de ClinVar y "
            "reflejan el consenso de los laboratorios que han enviado "
            "interpretaciones hasta la fecha de la anotacion. Estas "
            "clasificaciones pueden cambiar a medida que se acumule nueva "
            "evidencia. Las variantes de significado incierto (VUS) no deben "
            "interpretarse como patogenicas ni como benignas."
        ),
        (
            "Farmacogenomica",
            "Los alelos estrella farmacogenomicos se infieren a partir de "
            "variantes individuales (SNVs) y no incluyen analisis de numero "
            "de copias (CNV), que es particularmente relevante para CYP2D6. "
            "La asignacion de fenotipos metabolizadores requiere validacion "
            "clinica mediante ensayos farmacogenomicos certificados antes de "
            "modificar cualquier prescripcion."
        ),
        (
            "Estimacion de ancestria",
            "La ancestria se estima mediante correlacion de frecuencias "
            "alelicas con poblaciones de referencia continentales (1000 Genomes, "
            "gnomAD). Este metodo proporciona una aproximacion a nivel "
            "continental, pero no sustituye un analisis de componentes "
            "principales (PCA) formal con un panel de referencia denso como el "
            "Simons Genome Diversity Project. La baja superposicion de variantes "
            "entre WES y paneles genomicos de referencia limita la resolucion "
            "subcontinental."
        ),
        (
            "Rasgos de aptitud fisica y nutricion",
            "Las asociaciones genotipo-fenotipo para rasgos de nutricion y "
            "aptitud fisica se basan en estudios observacionales con tamanos "
            "de efecto generalmente modestos. Estos hallazgos no deben "
            "utilizarse como unica base para recomendaciones dieteticas o de "
            "ejercicio. Los grados de evidencia (A, B, C) reflejan el nivel "
            "de replicacion en la literatura, no la magnitud del efecto."
        ),
        (
            "Homocigosidad y consanguinidad",
            "La elevada ratio de homocigosidad observada en todas las muestras "
            "puede deberse a estructura poblacional, consanguinidad o al sesgo "
            "inherente de la captura exomica. Un analisis formal con datos "
            "genomicos completos y coeficientes de consanguinidad (F-statistics) "
            "seria necesario para distinguir entre estas causas."
        ),
        (
            "Validacion clinica",
            "Este informe se genera con fines de investigacion y no ha sido "
            "validado segun los estandares de un laboratorio clinico acreditado "
            "(ISO 15189). Todas las variantes con potencial relevancia clinica "
            "deben confirmarse mediante secuenciacion Sanger u otra metodologia "
            "independiente en un laboratorio certificado antes de tomar "
            "decisiones medicas."
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


# ── Build PDF (Spanish) ──────────────────────────────────────────────────

def build_sample_pdf_es(md_path, output_path):
    _require_reportlab()
    report = parse_markdown_report(md_path)
    sample_id = report["sample_id"] or md_path.stem.split("_")[0]
    styles = build_styles()

    doc = BaseDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=75 * mm,   # cover: logo strip + navy header
        bottomMargin=18 * mm,
        title=f"Informe Clinico WES - {sample_id}",
        author="M. Corpas - ClawBio WES Analysis Pipeline",
        subject="Informe de Secuenciacion del Exoma Completo",
    )
    doc._sample_id = sample_id

    # Cover page: logo strip (22mm) + navy block (48mm) + accent (2mm) = 72mm
    frame_cover = Frame(
        MARGIN, 18 * mm,
        PAGE_W - 2 * MARGIN, PAGE_H - 75 * mm - 18 * mm,
        id="cover_frame"
    )
    # Body pages: logo bar (12mm) + navy bar (11mm) + accent (1.5mm) ~ 25mm
    frame_body = Frame(
        MARGIN, 18 * mm,
        PAGE_W - 2 * MARGIN, PAGE_H - 28 * mm - 18 * mm,
        id="body_frame"
    )

    doc.addPageTemplates([
        PageTemplate(id="cover", frames=frame_cover, onPage=first_page_header_es),
        PageTemplate(id="body", frames=frame_body, onPage=header_footer_es),
    ])

    elements = []
    elements.append(NextPageTemplate("body"))

    # KPI cards
    if report["sections"]:
        kpis = extract_exome_kpis_es(report["sections"][0]["body"])
        if kpis:
            elements.append(Spacer(1, 2 * mm))
            elements.append(build_kpi_row(kpis, styles))
            elements.append(Spacer(1, 4 * mm))

    # ── Interpretive paragraph ───────────────────────────────────
    interp = build_interpretation_paragraph(report, styles)
    elements.extend(interp)

    # Render all sections in Spanish
    for section in report["sections"]:
        section_elements = render_section_es(section, styles)
        elements.extend(section_elements)
        elements.append(Spacer(1, 2 * mm))

    # ── Section 8: Ancestry ──────────────────────────────────────────
    elements.extend(_build_ancestry_section(sample_id, styles))

    # ── Section 9: Limitations ───────────────────────────────────────
    elements.extend(_build_limitations_section(styles))

    # Final disclaimer
    elements.append(Spacer(1, 6 * mm))
    elements.append(HRFlowable(width="100%", thickness=1, color=ACCENT_RED,
                               spaceAfter=3 * mm))

    disclaimer_es = (
        "AVISO IMPORTANTE: Este informe se genera exclusivamente con fines de "
        "investigacion y educacion. No es un informe de diagnostico clinico, no ha "
        "sido validado para uso clinico y no debe utilizarse para tomar decisiones "
        "medicas sin consultar a un profesional sanitario cualificado y realizar las "
        "pruebas confirmatorias apropiadas. Las clasificaciones de variantes reflejan "
        "las presentaciones de ClinVar en el momento de la anotacion y pueden cambiar "
        "conforme se acumule evidencia. Las predicciones farmacogenomicas requieren "
        "validacion clinica antes de modificar prescripciones farmacologicas."
    )
    elements.append(Paragraph(disclaimer_es, styles["disclaimer"]))

    elements.append(Spacer(1, 4 * mm))
    elements.append(Paragraph(
        f"Informe preparado por M. Corpas utilizando ClawBio WES Analysis "
        f"Pipeline v1.0 el {fecha_es()}.",
        styles["small"]
    ))
    elements.append(Paragraph(
        "University of Westminster, Londres, Reino Unido  |  github.com/ClawBio",
        styles["small"]
    ))

    doc.build(elements)
    return output_path


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generador de Informes Clinicos WES en Espanol (ClawBio Skill)"
    )
    parser.add_argument("--report-dir", type=str, default=str(REPORT_DIR),
                        help="Directorio con los informes WES en markdown")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_DIR),
                        help="Directorio de salida para los PDFs")
    parser.add_argument("--logo-left", type=str, default=LOGO_PREDICE,
                        help="Logo institucional izquierdo (JPG/PNG)")
    parser.add_argument("--logo-right", type=str, default=LOGO_INBIOMEDIC,
                        help="Logo institucional derecho (JPG/PNG)")
    parser.add_argument("--samples", type=str, default=None,
                        help="Lista de muestras separadas por coma (default: Sample1-7)")
    parser.add_argument("--demo", action="store_true",
                        help="Ejecutar con datos de ejemplo del proyecto Novogene")

    args = parser.parse_args()

    # Override module-level logo paths with CLI args
    import sys
    _mod = sys.modules[__name__]
    report_dir = Path(args.report_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _mod.LOGO_PREDICE = args.logo_left
    _mod.LOGO_INBIOMEDIC = args.logo_right

    if args.samples:
        samples = [s.strip() for s in args.samples.split(",")]
    else:
        samples = SAMPLES

    if args.demo:
        # Use default Novogene paths
        report_dir = REPORT_DIR
        output_dir = OUTPUT_DIR
        samples = SAMPLES

    print("Generador de Informes Clinicos PDF (Espanol)")
    print(f"Directorio de entrada: {report_dir}")
    print(f"Directorio de salida:  {output_dir}")
    print(f"Muestras: {', '.join(samples)}")
    print(f"{'=' * 60}")

    generated = []
    for sample in samples:
        md_file = report_dir / f"{sample}_WES_Report.md"
        if not md_file.exists():
            print(f"  OMITIDO  {sample} - informe markdown no encontrado")
            continue

        pdf_file = output_dir / f"{sample}_Informe_Clinico_WES.pdf"
        print(f"  Generando {sample}...", end="", flush=True)
        try:
            build_sample_pdf_es(md_file, pdf_file)
            size_kb = pdf_file.stat().st_size / 1024
            print(f"  OK ({size_kb:.0f} KB)")
            generated.append(pdf_file)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'=' * 60}")
    print(f"Generados {len(generated)}/{len(samples)} informes PDF")
    for f in generated:
        print(f"  {f}")


if __name__ == "__main__":
    main()
