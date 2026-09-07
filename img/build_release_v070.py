"""Build the ClawBio v0.7.0 release tile (Agent Skills, MCP server deprecated).

Sibling of build_mcp_release.py; same palette, raptor treatment and card system.
Frozen release graphic: the numbers are the release claims and must not drift after
publication.

    python3 build_release_v070.py

PNG rendering uses rsvg-convert (brew install librsvg) when present, else the
cairosvg Python package (needs brew's cairo). QuickLook's qlmanage is deliberately
not used: it renders at the display's 2x factor into a square canvas and silently
ships a cropped, zoomed tile.
"""
import base64
import io
import os
import shutil
import subprocess
from pathlib import Path

HERE = Path(os.path.dirname(os.path.abspath(__file__)))
RAPTOR = HERE / "clawbio-logo.jpeg"
VERSION = "0.7.0"
SLUG = "v" + VERSION.replace(".", "")
OUT_SVG = HERE / f"clawbio-release-{SLUG}.svg"
REPO = "ClawBio/ClawBio"

W, H = 1200, 627

SANS = "Helvetica Neue, Helvetica, Arial, sans-serif"
MONO = "SF Mono, Menlo, Consolas, monospace"

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
    sub = (f'<text x="{x + 26}" y="{y + 124}" font-family="{MONO}" font-size="12.5" '
           f'font-weight="500" letter-spacing="1.6" fill="#6e7681">{sublabel}</text>') if sublabel else ""
    return f'''{glow}
  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{CARD}" stroke="{stroke}" stroke-width="1"/>
  <text x="{x + 24}" y="{y + 72}" font-family="{SANS}" font-size="48" font-weight="800" fill="{GREEN}" letter-spacing="-1">{number}</text>
  <text x="{x + 26}" y="{y + 104}" font-family="{MONO}" font-size="15" font-weight="500" letter-spacing="2.2" fill="{MUTED}">{label}</text>
  {sub}'''


def render_png(svg: Path, png: Path, scale: int) -> bool:
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-w", str(W * scale), "-h", str(H * scale),
                        "-o", str(png), str(svg)], check=True)
        return True
    # cairocffi locates libcairo through ctypes.util.find_library, which on macOS
    # honours DYLD_FALLBACK_LIBRARY_PATH but does not look in Homebrew's prefix.
    for prefix in ("/opt/homebrew/lib", "/usr/local/lib"):
        if Path(prefix, "libcairo.2.dylib").is_file():
            os.environ.setdefault("DYLD_FALLBACK_LIBRARY_PATH", prefix)
    try:
        import cairosvg
    except (ImportError, OSError):
        return False
    cairosvg.svg2png(url=str(svg), write_to=str(png), output_width=W * scale, output_height=H * scale)
    return True


def main() -> None:
    raptor = raptor_b64()

    cw, gap, x0, cy, ch = 326, 30, 60, 384, 148
    cards = "\n".join([
        stat_card(x0, cy, cw, ch, "97", "SKILLS", "PLAIN AGENT SKILLS FOLDERS", highlight=True),
        stat_card(x0 + (cw + gap), cy, cw, ch, "4", "EDITORS", "CURSOR / VS CODE / CODEX / ZED"),
        stat_card(x0 + 2 * (cw + gap), cy, cw, ch, "0", "SERVERS", "NOTHING TO CONFIGURE"),
    ])

    cmd = "/plugin marketplace add ClawBio/ClawBio"

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <rect width="{W}" height="{H}" fill="{BG}"/>

  <image x="48" y="20" width="72" height="72" xlink:href="data:image/png;base64,{raptor}"/>
  <text x="132" y="74" font-family="{SANS}" font-size="46" font-weight="700" fill="{FG}">ClawBio</text>
  <text x="{W-50}" y="66" text-anchor="end" font-family="{MONO}" font-size="17" font-weight="500" letter-spacing="1.5" fill="{MUTED}">clawbio.ai</text>
  <line x1="0" y1="112" x2="{W}" y2="112" stroke="{BORDER}" stroke-width="1"/>

  <text x="62" y="176" font-family="{MONO}" font-size="16" font-weight="600" letter-spacing="4" fill="{GREEN}">RELEASE v{VERSION} &#183; AGENT SKILLS</text>

  <text x="60" y="244" font-family="{SANS}" font-size="60" font-weight="800" fill="{GREEN}" letter-spacing="-1.5">Genomics skills, in your editor</text>

  <text x="62" y="292" font-family="{SANS}" font-size="25" font-weight="400" fill="{MUTED}">No bridge required. The editors now read the skill folders directly,</text>
  <text x="62" y="326" font-family="{SANS}" font-size="25" font-weight="400" fill="{FG}">so the MCP server is retiring and 41 pull requests of hardening ship.</text>

  {cards}

  <rect x="60" y="556" width="700" height="46" rx="8" fill="{CARD}" stroke="{BORDER}" stroke-width="1"/>
  <text x="80" y="585" font-family="{MONO}" font-size="18" font-weight="500" fill="{GREEN}"><tspan fill="{MUTED}">&gt;</tspan><tspan dx="11">{cmd}</tspan></text>

  <text x="{W-50}" y="585" text-anchor="end" font-family="{MONO}" font-size="16" font-weight="500" letter-spacing="1" fill="{MUTED}">github.com/{REPO}  &#183;  MIT</text>
</svg>'''

    OUT_SVG.write_text(svg)
    print(f"wrote {OUT_SVG.name}")
    for scale, suffix in ((1, ""), (2, "@2x")):
        png = HERE / f"clawbio-release-{SLUG}{suffix}.png"
        if render_png(OUT_SVG, png, scale):
            print(f"wrote {png.name}")
        else:
            raise SystemExit("no SVG renderer: brew install librsvg, or pip install cairosvg")


if __name__ == "__main__":
    main()
