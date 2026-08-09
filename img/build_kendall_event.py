"""Build the ClawBio @ Kendall event graphics, matching docs.clawbio.ai.

Same Material-for-MkDocs "clawbio" dark theme as build_hero.py:
  bg #0d1117, cards #161b22 with #30363d hairline borders (10px radius),
  green #3fb950 bold headings (letter-spacing -0.02em), link blue #58a6ff,
  muted #8b949e, system sans + SF Mono for code-like labels.

Event details live in the card grid, with the capacity card carrying the
green hover-glow accent. Three formats: square (Luma cover / IG / X),
landscape (LinkedIn), story (vertical). Render to PNG with rsvg-convert.
"""
import base64
import io
import os
import subprocess
from pathlib import Path

HERE = Path(os.path.dirname(os.path.abspath(__file__)))
RAPTOR = HERE / "clawbio-logo.jpeg"

SANS = "Helvetica Neue, Helvetica, Arial, sans-serif"   # closest to the docs system sans
MONO = "SF Mono, Menlo, Consolas, monospace"            # docs --md-code-font

GREEN = "#3fb950"
FG = "#e6edf3"
MUTED = "#8b949e"
DIM = "#6e7681"
BLUE = "#58a6ff"
CARD = "#161b22"
BORDER = "#30363d"
BG = "#0d1117"

EYEBROW = "BOSTON  ·  AUGUST 5, 2026"
TITLE = "ClawBio @ Kendall"
SUB1 = "Agentic AI for Biology. A networking evening."
SUB2 = "Trustworthy. Auditable. Agentic genomics."
URL = "luma.com/plvedw0h"

_cache = {}


def raptor_b64() -> str:
    """Recolour the raptor to off-white with a transparent background."""
    if "r" in _cache:
        return _cache["r"]
    from PIL import Image
    img = Image.open(RAPTOR).convert("L")
    px = img.load()
    w, h = img.size
    out = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    opx = out.load()
    r, g, b = 0xE6, 0xED, 0xF3
    for y in range(h):
        for x in range(w):
            a = 255 - px[x, y]
            if a > 18:
                opx[x, y] = (r, g, b, a)
    buf = io.BytesIO()
    out.save(buf, format="PNG")
    _cache["r"] = base64.b64encode(buf.getvalue()).decode()
    return _cache["r"]


def raptor_cropped():
    """Raptor trimmed to its ink bounds. Returns (b64, width, height)."""
    if "c" in _cache:
        return _cache["c"]
    import base64 as _b64
    from PIL import Image
    raw = _b64.b64decode(raptor_b64())
    img = Image.open(io.BytesIO(raw))
    img = img.crop(img.getbbox())
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    _cache["c"] = (_b64.b64encode(buf.getvalue()).decode(), img.width, img.height)
    return _cache["c"]


def stat_card(x, y, w, h, number, label, sublabel="", highlight=False, size=52):
    stroke = GREEN if highlight else BORDER
    glow = ('<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="%s" opacity="0.10"/>'
            % (x - 6, y - 6, w + 12, h + 12, GREEN)) if highlight else ""
    sub = (f'<text x="{x + 30}" y="{y + 132}" font-family="{MONO}" font-size="14" '
           f'font-weight="500" letter-spacing="2" fill="{DIM}">{sublabel}</text>') if sublabel else ""
    return f'''{glow}
  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{stroke}" stroke-width="1"/>
  <text x="{x + 28}" y="{y + 78}" font-family="{SANS}" font-size="{size}" font-weight="800" fill="{GREEN}" letter-spacing="-1">{number}</text>
  <text x="{x + 30}" y="{y + 112}" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="2.5" fill="{MUTED}">{label}</text>
  {sub}'''


TERM = [
    ("$ claw run variant-classify --sample NA12878", FG),
    ("  plan: 4 steps  ·  tools: 3  ·  ref: GRCh38", MUTED),
    ("✓ ACMG PS3, PM2  ·  provenance: 4 sources", GREEN),
    ("⏸ abstained: 1 of 12 calls, evidence &lt; 2", BLUE),
]


def term_card(x, y, w, h, fs=17, lh=34):
    """Docs-style code block: the agent showing its work."""
    rows = "\n".join(
        f'  <text x="{x + 26}" y="{y + 52 + i * lh}" font-family="{MONO}" '
        f'font-size="{fs}" font-weight="500" fill="{c}" '
        f'xml:space="preserve">{t}</text>'
        for i, (t, c) in enumerate(TERM)
    )
    return f'''<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{BORDER}" stroke-width="1"/>
  <line x1="{x}" y1="{y + 34}" x2="{x + w}" y2="{y + 34}" stroke="{BORDER}" stroke-width="1"/>
  <circle cx="{x + 22}" cy="{y + 17}" r="5" fill="#f85149"/>
  <circle cx="{x + 40}" cy="{y + 17}" r="5" fill="#d29922"/>
  <circle cx="{x + 58}" cy="{y + 17}" r="5" fill="{GREEN}"/>
  <text x="{x + w - 20}" y="{y + 23}" text-anchor="end" font-family="{MONO}" font-size="13" letter-spacing="1.5" fill="{DIM}">AUDITABLE BY DEFAULT</text>
{rows}'''


def watermark(cx, cy, size, opacity=0.05):
    return (f'<image x="{cx}" y="{cy}" width="{size}" height="{size}" opacity="{opacity}" '
            f'xlink:href="data:image/png;base64,{raptor_b64()}"/>')


def glow(cx, cy, r):
    return f'''<defs><radialGradient id="gl{int(cx)}{int(cy)}">
    <stop offset="0" stop-color="{GREEN}" stop-opacity="0.16"/>
    <stop offset="1" stop-color="{GREEN}" stop-opacity="0"/>
  </radialGradient></defs>
  <circle cx="{cx}" cy="{cy}" r="{r}" fill="url(#gl{int(cx)}{int(cy)})"/>'''


def header(w, logo=72, y=20, name=46):
    return f'''<image x="48" y="{y}" width="{logo}" height="{logo}" xlink:href="data:image/png;base64,{raptor_b64()}"/>
  <text x="{48 + logo + 12}" y="{y + logo * 0.75:.0f}" font-family="{SANS}" font-size="{name}" font-weight="700" fill="{FG}">ClawBio</text>
  <text x="{w-50}" y="{y + 46}" text-anchor="end" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1.5" fill="{MUTED}">clawbio.ai</text>
  <line x1="0" y1="{y + logo + 20}" x2="{w}" y2="{y + logo + 20}" stroke="{BORDER}" stroke-width="1"/>'''


def cta(x, y, w=330, label="Register on Luma"):
    return f'''<g transform="translate({x},{y})">
    <rect x="0" y="0" width="{w}" height="44" rx="8" fill="{GREEN}"/>
    <path d="M22,8 L26,17.5 L36,18.3 L28.4,24.9 L30.8,34.6 L22,29.4 L13.2,34.6 L15.6,24.9 L8,18.3 L18,17.5 Z" fill="{BG}"/>
    <text x="46" y="29" font-family="{SANS}" font-size="19" font-weight="700" fill="#0b2a13">{label}</text>
  </g>'''


def svg_open(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{w}" height="{h}" viewBox="0 0 {w} {h}">\n  <rect width="{w}" height="{h}" fill="{BG}"/>')


CARDS = [
    dict(number="6:00 PM", label="WED, AUG 5", sublabel="UNTIL 8:30 PM", size=44),
    dict(number="Catalyst", label="TECHNOLOGY SQ", sublabel="CAMBRIDGE, MA", size=44),
    dict(number="Free", label="TO ATTEND", sublabel="RSVP REQUIRED", highlight=True, size=44),
]


def landscape():
    W, H = 1200, 627
    cw, gap, x0, cy, ch = 326, 30, 60, 380, 156
    cards = "\n".join(
        stat_card(x0 + i * (cw + gap), cy, cw, ch, **c) for i, c in enumerate(CARDS)
    )
    return f'''{svg_open(W, H)}
  {watermark(880, 96, 330, 0.055)}
  {header(W)}
  <text x="62" y="190" font-family="{MONO}" font-size="16" font-weight="600" letter-spacing="4" fill="{GREEN}">{EYEBROW}</text>
  <text x="60" y="256" font-family="{SANS}" font-size="62" font-weight="800" fill="{GREEN}" letter-spacing="-1.5">{TITLE}</text>
  <text x="62" y="310" font-family="{SANS}" font-size="27" font-weight="400" fill="{MUTED}">{SUB1}</text>
  <text x="62" y="350" font-family="{SANS}" font-size="27" font-weight="400" fill="{FG}">{SUB2}</text>
  {cards}
  {cta(60, 560)}
  <text x="410" y="589" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1" fill="{MUTED}">{URL}  &#183;  <tspan fill="{BLUE}">free to attend</tspan></text>
</svg>'''


def square():
    W = H = 1200
    cw, gap, x0, cy, ch = 344, 24, 60, 730, 168
    cards = "\n".join(
        stat_card(x0 + i * (cw + gap), cy, cw, ch, **c) for i, c in enumerate(CARDS)
    )
    return f'''{svg_open(W, H)}
  {watermark(820, 210, 420, 0.05)}
  {header(W)}
  <text x="62" y="250" font-family="{MONO}" font-size="17" font-weight="600" letter-spacing="4" fill="{GREEN}">{EYEBROW}</text>
  <text x="60" y="340" font-family="{SANS}" font-size="76" font-weight="800" fill="{GREEN}" letter-spacing="-2">{TITLE}</text>
  <text x="62" y="404" font-family="{SANS}" font-size="30" font-weight="400" fill="{MUTED}">{SUB1}</text>
  <text x="62" y="448" font-family="{SANS}" font-size="30" font-weight="400" fill="{FG}">{SUB2}</text>
  {term_card(60, 500, 1080, 190)}
  {cards}
  {cta(60, 1064)}
  <text x="410" y="1093" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1" fill="{MUTED}">{URL}  &#183;  <tspan fill="{BLUE}">free to attend</tspan></text>
</svg>'''


def story():
    W, H = 1080, 1920
    cw, gap, x0, cy, ch = 960, 26, 60, 860, 150
    cards = "\n".join(
        stat_card(x0, cy + i * (ch + gap), cw, ch, **c) for i, c in enumerate(CARDS)
    )
    return f'''{svg_open(W, H)}
  {watermark(600, 1180, 520, 0.05)}
  {header(W, logo=80, y=300, name=48)}
  <text x="62" y="500" font-family="{MONO}" font-size="18" font-weight="600" letter-spacing="4" fill="{GREEN}">{EYEBROW}</text>
  <text x="60" y="584" font-family="{SANS}" font-size="78" font-weight="800" fill="{GREEN}" letter-spacing="-2">ClawBio</text>
  <text x="60" y="668" font-family="{SANS}" font-size="78" font-weight="800" fill="{GREEN}" letter-spacing="-2">@ Kendall</text>
  <text x="62" y="726" font-family="{SANS}" font-size="30" font-weight="400" fill="{MUTED}">{SUB1}</text>
  <text x="62" y="770" font-family="{SANS}" font-size="30" font-weight="400" fill="{FG}">{SUB2}</text>
  {cards}
  {cta(60, 1470, w=360)}
  <text x="62" y="1560" font-family="{MONO}" font-size="18" font-weight="500" letter-spacing="1" fill="{MUTED}">{URL}  &#183;  <tspan fill="{BLUE}">free to attend</tspan></text>
</svg>'''


def tile():
    """Logo-forward brand tile: raptor dominant, event lockup beneath."""
    W = H = 1200
    cx = W / 2
    rb64, rw, rh = raptor_cropped()
    ih = 470
    iw = ih * rw / rh
    return f'''{svg_open(W, H)}
  <rect x="32" y="32" width="{W-64}" height="{H-64}" rx="16" fill="none" stroke="{BORDER}" stroke-width="1"/>
  <image x="{cx-iw/2:.1f}" y="200" width="{iw:.1f}" height="{ih}" xlink:href="data:image/png;base64,{rb64}"/>
  <text x="{cx}" y="132" text-anchor="middle" font-family="{MONO}" font-size="24" font-weight="600" letter-spacing="6" fill="{GREEN}">AGENTIC AI FOR BIOLOGY</text>
  <text x="{cx}" y="782" text-anchor="middle" font-family="{SANS}" font-size="112" font-weight="800" fill="{FG}" letter-spacing="-3">ClawBio</text>
  <text x="{cx}" y="896" text-anchor="middle" font-family="{SANS}" font-size="112" font-weight="800" fill="{GREEN}" letter-spacing="-3">@ Kendall</text>
  <line x1="{cx-300}" y1="954" x2="{cx+300}" y2="954" stroke="{BORDER}" stroke-width="1"/>
  <text x="{cx}" y="1022" text-anchor="middle" font-family="{MONO}" font-size="32" font-weight="600" letter-spacing="2" fill="{FG}">WED, AUGUST 5, 2026 · 6:00 PM</text>
  <text x="{cx}" y="1072" text-anchor="middle" font-family="{MONO}" font-size="25" font-weight="500" letter-spacing="1.5" fill="{MUTED}">CATALYST · TECHNOLOGY SQ, CAMBRIDGE MA</text>
  <text x="{cx}" y="1124" text-anchor="middle" font-family="{MONO}" font-size="25" font-weight="600" letter-spacing="2" fill="{GREEN}">FREE TO ATTEND · RSVP REQUIRED</text>
</svg>'''


def main() -> None:
    jobs = [
        ("clawbio-kendall-tile", tile(), 1200, 1200),
        ("clawbio-kendall-square", square(), 1200, 1200),
        ("clawbio-kendall-linkedin", landscape(), 1200, 627),
        ("clawbio-kendall-story", story(), 1080, 1920),
    ]
    for name, svg, w, h in jobs:
        svg_path = HERE / f"{name}.svg"
        svg_path.write_text(svg)
        for scale, suffix in ((1, ""), (2, "@2x")):
            png = HERE / f"{name}{suffix}.png"
            subprocess.run(
                ["rsvg-convert", "-w", str(w * scale), "-h", str(h * scale),
                 "-o", str(png), str(svg_path)],
                check=True,
            )
            print(f"wrote {png.name}")


if __name__ == "__main__":
    main()
