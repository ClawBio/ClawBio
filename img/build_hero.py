"""Build the ClawBio 1,000-stars hero image, matching docs.clawbio.ai.

Faithful to the Material-for-MkDocs "clawbio" dark theme:
  bg #0d1117, cards #161b22 with #30363d hairline borders (10px radius),
  green #3fb950 bold headings (letter-spacing -0.02em), link blue #58a6ff,
  muted #8b949e, system sans + SF Mono for code-like labels.

Live metrics row (forks/contributors/commits) styled as the docs card grid,
with the stars hero card carrying the green hover-glow accent.
Render to PNG with rsvg-convert.
"""
import base64
import io
import os
from pathlib import Path

HERE = Path(os.path.dirname(os.path.abspath(__file__)))
RAPTOR = HERE / "clawbio-logo.jpeg"
OUT_SVG = HERE / "clawbio-1k-stars-hero.svg"

W, H = 1200, 627

SANS = "Helvetica Neue, Helvetica, Arial, sans-serif"   # closest to the docs system sans
MONO = "SF Mono, Menlo, Consolas, monospace"            # docs --md-code-font

GREEN = "#3fb950"
FG = "#e6edf3"
MUTED = "#8b949e"
BLUE = "#58a6ff"
CARD = "#161b22"
BORDER = "#30363d"
BG = "#0d1117"


def raptor_b64() -> str:
    """Recolour the raptor to off-white with a transparent background."""
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
    return base64.b64encode(buf.getvalue()).decode()


def stat_card(x, y, w, h, number, label, sublabel="", highlight=False):
    stroke = GREEN if highlight else BORDER
    glow = ('<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="%s" opacity="0.10"/>'
            % (x - 6, y - 6, w + 12, h + 12, GREEN)) if highlight else ""
    sub = (f'<text x="{x + 30}" y="{y + 132}" font-family="{MONO}" font-size="14" '
           f'font-weight="500" letter-spacing="2" fill="#6e7681">{sublabel}</text>') if sublabel else ""
    return f'''{glow}
  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{stroke}" stroke-width="1"/>
  <text x="{x + 28}" y="{y + 78}" font-family="{SANS}" font-size="52" font-weight="800" fill="{GREEN}" letter-spacing="-1">{number}</text>
  <text x="{x + 30}" y="{y + 112}" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="2.5" fill="{MUTED}">{label}</text>
  {sub}'''


def main() -> None:
    raptor = raptor_b64()

    cw, gap, x0, cy, ch = 326, 30, 60, 380, 156
    cards = "\n".join([
        stat_card(x0, cy, cw, ch, "1,007", "STARS", "GITHUB MILESTONE", highlight=True),
        stat_card(x0 + (cw + gap), cy, cw, ch, "6,416", "CLONES", "LAST 14 DAYS"),
        stat_card(x0 + 2 * (cw + gap), cy, cw, ch, "41", "CONTRIBUTORS", "AND GROWING"),
    ])

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <rect width="{W}" height="{H}" fill="{BG}"/>

  <!-- header bar, like the docs top nav -->
  <image x="48" y="20" width="72" height="72" xlink:href="data:image/png;base64,{raptor}"/>
  <text x="132" y="74" font-family="{SANS}" font-size="46" font-weight="700" fill="{FG}">ClawBio</text>
  <text x="{W-50}" y="66" text-anchor="end" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1.5" fill="{MUTED}">clawbio.github.io/ClawBio</text>
  <line x1="0" y1="112" x2="{W}" y2="112" stroke="{BORDER}" stroke-width="1"/>

  <!-- eyebrow -->
  <text x="62" y="190" font-family="{MONO}" font-size="16" font-weight="600" letter-spacing="4" fill="{GREEN}">MILESTONE</text>

  <!-- hero title, green bold like docs h1 -->
  <text x="60" y="256" font-family="{SANS}" font-size="62" font-weight="800" fill="{GREEN}" letter-spacing="-1.5">1,000 GitHub stars</text>

  <!-- subtitle, muted -->
  <text x="62" y="310" font-family="{SANS}" font-size="27" font-weight="400" fill="{MUTED}">The bioinformatics-native AI agent skill library.</text>
  <text x="62" y="350" font-family="{SANS}" font-size="27" font-weight="400" fill="{FG}">Trustworthy. Auditable. Agentic genomics.</text>

  <!-- live-metrics card grid -->
  {cards}

  <!-- footer: star call-to-action -->
  <g transform="translate(60,560)">
    <rect x="0" y="0" width="312" height="44" rx="8" fill="{GREEN}"/>
    <path d="M22,8 L26,17.5 L36,18.3 L28.4,24.9 L30.8,34.6 L22,29.4 L13.2,34.6 L15.6,24.9 L8,18.3 L18,17.5 Z" fill="#0d1117"/>
    <text x="46" y="29" font-family="{SANS}" font-size="19" font-weight="700" fill="#0b2a13">Star us on GitHub</text>
  </g>
  <text x="392" y="589" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1" fill="{MUTED}">github.com/ClawBio/ClawBio  &#183;  <tspan fill="{BLUE}">MIT licensed</tspan></text>
</svg>'''
    OUT_SVG.write_text(svg)
    print(f"wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
