"""
viewer.py — Self-contained 3Dmol.js HTML viewer coloured by pLDDT.

3Dmol.js is embedded inline so the file opens without a network connection.
On first use the library is downloaded from cdnjs and cached locally at
~/.clawbio/3dmol-min.js.  If the download fails the viewer falls back to a
CDN <script> tag with a visible warning in the page.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_CDN_URL   = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"
_CACHE_PATH = Path.home() / ".clawbio" / "3dmol-min.js"


def _get_3dmol_script() -> str:
    """Return an HTML script block with 3Dmol.js embedded inline.

    Tries (in order):
    1. Local cache at ~/.clawbio/3dmol-min.js
    2. Download from CDN and save to cache
    3. Fall back to a <script src=...> CDN tag (requires internet)
    """
    if _CACHE_PATH.exists():
        return f"<script>{_CACHE_PATH.read_text()}</script>"

    try:
        import urllib.request
        js = urllib.request.urlopen(_CDN_URL, timeout=15).read().decode()
        _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        _CACHE_PATH.write_text(js)
        return f"<script>{js}</script>"
    except Exception:
        pass

    # Offline fallback — viewer will not work without internet
    return (
        f'<script src="{_CDN_URL}"></script>\n'
        '  <p style="color:#f66;padding:8px">'
        "⚠ 3Dmol.js could not be loaded offline. "
        "Run with internet access once to cache the library.</p>"
    )


# AlphaFold / Boltz pLDDT colour bands (hex, same as report.py BANDS)
def _plddt_colour(score: float) -> str:
    if score >= 90:
        return "#0053D6"
    if score >= 70:
        return "#65CBF3"
    if score >= 50:
        return "#FFDB13"
    return "#FF7D45"


def generate_viewer_html(
    output_path: Path,
    cif_path: Path,
    plddt: np.ndarray,
    chain_boundaries: list[dict],
) -> None:
    """Write a self-contained HTML file with an embedded 3Dmol.js viewer.

    The CIF structure is inlined as a JavaScript string so the file opens
    directly from the filesystem with no server required.

    Args:
        output_path: Destination HTML file.
        cif_path: Path to the predicted CIF structure.
        plddt: Per-residue pLDDT array (0–100).
        chain_boundaries: List of {chain_id, start, end} dicts from confidence.py.
    """
    cif_text = cif_path.read_text()
    # Escape backticks so the string is safe inside a JS template literal
    cif_escaped = cif_text.replace("\\", "\\\\").replace("`", "\\`")

    # Build per-residue colour spec for 3Dmol setStyle calls
    # Group consecutive residues of the same colour to reduce JS payload
    colour_calls = _build_colour_calls(plddt, chain_boundaries)
    legend_html  = _build_legend_html()
    script_block = _get_3dmol_script()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <title>Struct Predictor — 3D Viewer</title>
  {script_block}
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: system-ui, sans-serif; background: #1a1a2e; color: #e0e0e0; }}
    #header {{ padding: 14px 20px; background: #16213e; border-bottom: 1px solid #0f3460; }}
    #header h1 {{ font-size: 1.1rem; font-weight: 600; color: #e0e0e0; }}
    #header p  {{ font-size: 0.75rem; color: #888; margin-top: 2px; }}
    #main {{ display: flex; height: calc(100vh - 56px); }}
    #viewer {{ flex: 1; position: relative; }}
    #sidebar {{
      width: 220px; background: #16213e; padding: 16px 14px;
      border-left: 1px solid #0f3460; overflow-y: auto;
      display: flex; flex-direction: column; gap: 16px;
    }}
    #sidebar h2 {{ font-size: 0.8rem; font-weight: 600; text-transform: uppercase;
                   letter-spacing: 0.05em; color: #aaa; }}
    .legend-item {{ display: flex; align-items: center; gap: 8px; font-size: 0.78rem; }}
    .swatch {{ width: 14px; height: 14px; border-radius: 3px; flex-shrink: 0; }}
    .stat {{ font-size: 0.78rem; color: #ccc; }}
    .stat span {{ font-weight: 600; color: #fff; }}
    #controls {{ display: flex; flex-direction: column; gap: 6px; }}
    button {{
      background: #0f3460; border: 1px solid #1a5276; color: #e0e0e0;
      border-radius: 4px; padding: 5px 8px; font-size: 0.75rem; cursor: pointer;
    }}
    button:hover {{ background: #1a5276; }}
    #disclaimer {{
      font-size: 0.65rem; color: #666; line-height: 1.4;
      border-top: 1px solid #0f3460; padding-top: 10px;
    }}
  </style>
</head>
<body>
  <div id="header">
    <h1>Struct Predictor — 3D Structure Viewer</h1>
    <p>Coloured by pLDDT confidence &nbsp;|&nbsp; Boltz-2</p>
  </div>
  <div id="main">
    <div id="viewer"></div>
    <div id="sidebar">
      <div>
        <h2>pLDDT Legend</h2>
        <div style="display:flex;flex-direction:column;gap:6px;margin-top:8px;">
          {legend_html}
        </div>
      </div>
      <div>
        <h2>Controls</h2>
        <div id="controls">
          <button onclick="setStyle('cartoon')">Cartoon</button>
          <button onclick="setStyle('stick')">Stick</button>
          <button onclick="setStyle('sphere')">Sphere</button>
          <button onclick="viewer.zoomTo()">Reset view</button>
          <button onclick="viewer.spin(!spinning); spinning=!spinning;">Toggle spin</button>
        </div>
      </div>
      <div id="disclaimer">
        ClawBio is a research and educational tool. It is not a medical device
        and does not provide clinical diagnoses. Consult a healthcare
        professional before making any medical decisions.
      </div>
    </div>
  </div>

  <script>
    const CIF_DATA = `{cif_escaped}`;

    const COLOUR_CALLS = {colour_calls};

    let currentStyle = 'cartoon';
    let spinning = false;

    const viewer = $3Dmol.createViewer('viewer', {{
      backgroundColor: '#1a1a2e',
    }});

    viewer.addModel(CIF_DATA, 'cif');

    function applyColours(style) {{
      viewer.setStyle({{}}, {{}});
      for (const [sel, colour] of COLOUR_CALLS) {{
        const styleSpec = style === 'cartoon'
          ? {{ cartoon: {{ color: colour }} }}
          : style === 'stick'
          ? {{ stick: {{ color: colour, radius: 0.15 }} }}
          : {{ sphere: {{ color: colour, scale: 0.4 }} }};
        viewer.setStyle(sel, styleSpec);
      }}
    }}

    function setStyle(s) {{
      currentStyle = s;
      applyColours(s);
      viewer.render();
    }}

    applyColours('cartoon');
    viewer.zoomTo();
    viewer.render();
  </script>
</body>
</html>
"""
    output_path.write_text(html)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_colour_calls(plddt: np.ndarray, chain_boundaries: list[dict]) -> str:
    """Build a compact JS array of [selector, colour] pairs."""
    # Map residue index → chain_id and 1-based resi
    chain_map: list[tuple[str, int]] = []
    for cb in chain_boundaries:
        chain_id = cb["chain_id"]
        for res_idx in range(cb["start"], cb["end"] + 1):
            resi = res_idx - cb["start"] + 1
            chain_map.append((chain_id, resi))

    # Group consecutive residues with the same colour
    groups: list[tuple[str, str, list[int]]] = []  # (chain_id, colour, [resi, ...])
    for i, score in enumerate(plddt):
        colour = _plddt_colour(float(score))
        chain_id, resi = chain_map[i] if i < len(chain_map) else ("A", i + 1)
        if groups and groups[-1][0] == chain_id and groups[-1][1] == colour:
            groups[-1][2].append(resi)
        else:
            groups.append((chain_id, colour, [resi]))

    lines = []
    for chain_id, colour, residues in groups:
        resi_str = ",".join(str(r) for r in residues)
        lines.append(f'  [{{chain: "{chain_id}", resi: [{resi_str}]}}, "{colour}"]')

    return "[\n" + ",\n".join(lines) + "\n]"


def _build_legend_html() -> str:
    items = [
        ("#0053D6", "Very high  ≥ 90"),
        ("#65CBF3", "High  70–90"),
        ("#FFDB13", "Low  50–70"),
        ("#FF7D45", "Very low  &lt; 50"),
    ]
    parts = []
    for colour, label in items:
        parts.append(
            f'<div class="legend-item">'
            f'<div class="swatch" style="background:{colour}"></div>'
            f'<span>{label}</span></div>'
        )
    return "\n          ".join(parts)
