"""Build the ClawBio Hackathon Berlin graphics, co-branded with Nebius.

Same layout as build_kendall_event.py so the two events read as one series:
header bar with hairline, mono eyebrow, green title, two subtitles, a docs-style
code block, a three-card grid, and a green CTA. The only structural difference is
the Nebius badge in the header, top right, where the Kendall tile puts clawbio.ai.

Theme is the Material-for-MkDocs "clawbio" dark theme, as build_hero.py:
bg #0d1117, cards #161b22 with #30363d hairlines at 10px radius, green #3fb950
headings, muted #8b949e, system sans + SF Mono.

Nebius mark is the official one from nebius.com/logo.svg, rasterised to
nebius-logo.png. Replace that file to update it.
"""
import base64
import io
import os
import subprocess
from pathlib import Path

HERE = Path(os.path.dirname(os.path.abspath(__file__)))
RAPTOR = HERE / "clawbio-logo.jpeg"
NEBIUS_LOGO = HERE / "nebius-logo.png"

SANS = "Helvetica Neue, Helvetica, Arial, sans-serif"
MONO = "SF Mono, Menlo, Consolas, monospace"

GREEN = "#3fb950"
FG = "#e6edf3"
MUTED = "#8b949e"
DIM = "#6e7681"
BLUE = "#58a6ff"
CARD = "#161b22"
BORDER = "#30363d"
BG = "#0d1117"

EYEBROW = "BERLIN  ·  AUGUST 18, 2026"
TITLE = "ClawBio Hackathon"
SUB1 = "Agentic AI for Genomics. One day, teams of five."
SUB2 = "Build a ClawBio skill. Demo it. Compute provided."
URL = "luma.com/clawbio-q8pw"

CARDS = [
    dict(number="12:00 PM", label="TUE, AUG 18", sublabel="UNTIL 6:00 PM", size=44),
    dict(number="Impact Hub", label="THE LOOP", sublabel="ROLLBERGSTR. 28A, BERLIN", size=38),
    dict(number="Free", label="TO ATTEND", sublabel="COMPUTE + CREDITS", highlight=True, size=44),
]

TERM = [
    ("$ claw run rnaseq-de --counts salmon/ --design design.tsv", FG),
    ("  plan: 6 steps  ·  tools: 4  ·  compute: nebius", MUTED),
    ("✓ 312 DE genes  ·  provenance: 6 sources, 0 inferred", GREEN),
    ("⏸ abstained: low-count filter, evidence &lt; 2", BLUE),
]

_c = {}


def raptor_b64() -> str:
    if "r" in _c:
        return _c["r"]
    from PIL import Image
    img = Image.open(RAPTOR).convert("L")
    px = img.load()
    w, h = img.size
    out = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    opx = out.load()
    for y in range(h):
        for x in range(w):
            a = 255 - px[x, y]
            if a > 18:
                opx[x, y] = (0xE6, 0xED, 0xF3, a)
    buf = io.BytesIO()
    out.save(buf, format="PNG")
    _c["r"] = base64.b64encode(buf.getvalue()).decode()
    return _c["r"]


def raptor_cropped():
    """Raptor trimmed to its ink bounds. Returns (b64, width, height)."""
    if "c" in _c:
        return _c["c"]
    from PIL import Image
    img = Image.open(io.BytesIO(base64.b64decode(raptor_b64())))
    img = img.crop(img.getbbox())
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    _c["c"] = (base64.b64encode(buf.getvalue()).decode(), img.width, img.height)
    return _c["c"]


def nebius_centred(cx, cy, h=44):
    """Official Nebius badge, centred on cx."""
    if not NEBIUS_LOGO.exists():
        return (f'<text x="{cx}" y="{cy+h*0.35:.0f}" text-anchor="middle" font-family="{SANS}" '
                f'font-size="{h*0.8:.0f}" font-weight="700" letter-spacing="3" fill="{FG}">NEBIUS</text>')
    from PIL import Image
    w0, h0 = Image.open(NEBIUS_LOGO).size
    b64 = base64.b64encode(NEBIUS_LOGO.read_bytes()).decode()
    w = h * w0 / h0
    return (f'<image x="{cx-w/2:.0f}" y="{cy-h/2:.0f}" width="{w:.0f}" height="{h}" '
            f'xlink:href="data:image/png;base64,{b64}"/>')


def tile():
    """Logo-forward brand tile, matching clawbio-kendall-tile, plus Nebius."""
    W = H = 1200
    cx = W / 2
    rb64, rw, rh = raptor_cropped()
    ih = 380
    iw = ih * rw / rh
    return f'''{head(W, H)}
  <rect x="32" y="32" width="{W-64}" height="{H-64}" rx="16" fill="none" stroke="{BORDER}" stroke-width="1"/>
  {nebius_centred(cx, 128, 84)}
  <text x="{cx}" y="228" text-anchor="middle" font-family="{MONO}" font-size="26" font-weight="600" letter-spacing="6" fill="{GREEN}">AGENTIC AI FOR GENOMICS</text>
  <image x="{cx-iw/2:.1f}" y="268" width="{iw:.1f}" height="{ih}" xlink:href="data:image/png;base64,{rb64}"/>
  <text x="{cx}" y="740" text-anchor="middle" font-family="{SANS}" font-size="108" font-weight="800" fill="{FG}" letter-spacing="-3">ClawBio</text>
  <text x="{cx}" y="850" text-anchor="middle" font-family="{SANS}" font-size="108" font-weight="800" fill="{GREEN}" letter-spacing="-3">@ Berlin</text>
  <line x1="{cx-320}" y1="908" x2="{cx+320}" y2="908" stroke="{BORDER}" stroke-width="1"/>
  <text x="{cx}" y="972" text-anchor="middle" font-family="{MONO}" font-size="34" font-weight="600" letter-spacing="2" fill="{FG}">TUE, AUG 18, 2026 · 12:00</text>
  <text x="{cx}" y="1026" text-anchor="middle" font-family="{MONO}" font-size="28" font-weight="500" letter-spacing="1.5" fill="{MUTED}">IMPACT HUB · ROLLBERGSTR. 28A</text>
  <text x="{cx}" y="1082" text-anchor="middle" font-family="{MONO}" font-size="28" font-weight="600" letter-spacing="2" fill="{GREEN}">FREE · COMPUTE PROVIDED</text>
</svg>'''


def nebius_mark(right_x, cy, h=38):
    """Official Nebius badge, right-aligned to right_x, vertically centred on cy."""
    if not NEBIUS_LOGO.exists():
        return (f'<text x="{right_x}" y="{cy+8}" text-anchor="end" font-family="{SANS}" '
                f'font-size="{h*0.8:.0f}" font-weight="700" letter-spacing="3" fill="{FG}">NEBIUS</text>')
    from PIL import Image
    w0, h0 = Image.open(NEBIUS_LOGO).size
    b64 = base64.b64encode(NEBIUS_LOGO.read_bytes()).decode()
    w = h * w0 / h0
    return (f'<image x="{right_x-w:.0f}" y="{cy-h/2:.0f}" width="{w:.0f}" height="{h}" '
            f'xlink:href="data:image/png;base64,{b64}"/>')


def stat_card(x, y, w, h, number, label, sublabel="", highlight=False, size=52):
    stroke = GREEN if highlight else BORDER
    glow = ('<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="%s" opacity="0.10"/>'
            % (x - 6, y - 6, w + 12, h + 12, GREEN)) if highlight else ""
    sub = (f'<text x="{x+30}" y="{y+132}" font-family="{MONO}" font-size="14" font-weight="500" '
           f'letter-spacing="2" fill="{DIM}">{sublabel}</text>') if sublabel else ""
    return f'''{glow}
  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{stroke}" stroke-width="1"/>
  <text x="{x+28}" y="{y+78}" font-family="{SANS}" font-size="{size}" font-weight="800" fill="{GREEN}" letter-spacing="-1">{number}</text>
  <text x="{x+30}" y="{y+112}" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="2.5" fill="{MUTED}">{label}</text>
  {sub}'''


def term_card(x, y, w, h, fs=17, lh=34):
    rows = "\n".join(
        f'  <text x="{x+26}" y="{y+52+i*lh}" font-family="{MONO}" font-size="{fs}" '
        f'font-weight="500" fill="{c}" xml:space="preserve">{t}</text>'
        for i, (t, c) in enumerate(TERM))
    return f'''<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{BORDER}" stroke-width="1"/>
  <line x1="{x}" y1="{y+34}" x2="{x+w}" y2="{y+34}" stroke="{BORDER}" stroke-width="1"/>
  <circle cx="{x+22}" cy="{y+17}" r="5" fill="#f85149"/>
  <circle cx="{x+40}" cy="{y+17}" r="5" fill="#d29922"/>
  <circle cx="{x+58}" cy="{y+17}" r="5" fill="{GREEN}"/>
  <text x="{x+w-20}" y="{y+23}" text-anchor="end" font-family="{MONO}" font-size="13" letter-spacing="1.5" fill="{DIM}">AUDITABLE BY DEFAULT</text>
{rows}'''


def watermark(x, y, size, opacity=0.05):
    return (f'<image x="{x}" y="{y}" width="{size}" height="{size}" opacity="{opacity}" '
            f'xlink:href="data:image/png;base64,{raptor_b64()}"/>')


def header(w, logo=72, y=20, name=46):
    """ClawBio left, Nebius badge right, hairline under. Mirrors the Kendall header."""
    return f'''<image x="48" y="{y}" width="{logo}" height="{logo}" xlink:href="data:image/png;base64,{raptor_b64()}"/>
  <text x="{48+logo+12}" y="{y+logo*0.75:.0f}" font-family="{SANS}" font-size="{name}" font-weight="700" fill="{FG}">ClawBio</text>
  <text x="{w-260}" y="{y+logo*0.62:.0f}" text-anchor="end" font-family="{SANS}" font-size="26" font-weight="400" fill="{DIM}">&#215;</text>
  {nebius_mark(w-50, y+logo*0.5, 40)}
  <line x1="0" y1="{y+logo+20}" x2="{w}" y2="{y+logo+20}" stroke="{BORDER}" stroke-width="1"/>'''


def cta(x, y, w=330, label="Register on Luma"):
    return f'''<g transform="translate({x},{y})">
    <rect x="0" y="0" width="{w}" height="44" rx="8" fill="{GREEN}"/>
    <path d="M22,8 L26,17.5 L36,18.3 L28.4,24.9 L30.8,34.6 L22,29.4 L13.2,34.6 L15.6,24.9 L8,18.3 L18,17.5 Z" fill="{BG}"/>
    <text x="46" y="29" font-family="{SANS}" font-size="19" font-weight="700" fill="#0b2a13">{label}</text>
  </g>'''


def head(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{w}" height="{h}" viewBox="0 0 {w} {h}">\n'
            f'  <rect width="{w}" height="{h}" fill="{BG}"/>')


def square():
    W = H = 1200
    cw, gap, x0, cy, ch = 344, 24, 60, 730, 168
    cards = "\n".join(stat_card(x0 + i*(cw+gap), cy, cw, ch, **c) for i, c in enumerate(CARDS))
    return f'''{head(W, H)}
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


def landscape():
    W, H = 1200, 627
    cw, gap, x0, cy, ch = 326, 30, 60, 380, 156
    cards = "\n".join(stat_card(x0 + i*(cw+gap), cy, cw, ch, **c) for i, c in enumerate(CARDS))
    return f'''{head(W, H)}
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


def main() -> None:
    for name, svg, w, h in [("clawbio-berlin-tile", tile(), 1200, 1200),
                            ("clawbio-berlin-square", square(), 1200, 1200),
                            ("clawbio-berlin-linkedin", landscape(), 1200, 627)]:
        p = HERE / f"{name}.svg"
        p.write_text(svg)
        for scale, suffix in ((1, ""), (2, "@2x")):
            png = HERE / f"{name}{suffix}.png"
            subprocess.run(["rsvg-convert", "-w", str(w*scale), "-h", str(h*scale),
                            "-o", str(png), str(p)], check=True)
            print(f"wrote {png.name}")


if __name__ == "__main__":
    main()
