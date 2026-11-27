"""
14-pca.py -- PCA of Genomebook agents based on 26 genetic traits

Computes PCA from trait_scores across all genomes, exports an interactive
scatter plot. Color by generation, founder lineage, sex, or health.

Usage:
    python 14-pca.py
    python 14-pca.py --max-gen 5
"""

import argparse
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
GENOMES_DIR = DATA / "GENOMES"
DEFAULT_OUTPUT = BASE.parent / "slides" / "genomebook" / "pca.html"


def load_genomes(max_gen=None):
    genomes = []
    for gf in sorted(GENOMES_DIR.glob("*.genome.json")):
        g = json.load(open(gf))
        if max_gen is not None and g.get("generation", 0) > max_gen:
            continue
        genomes.append(g)
    return genomes


def get_founder_lineage(genome, all_genomes_map):
    """Trace back to generation-0 ancestors."""
    founders = set()
    visited = set()
    queue = [genome["id"]]
    while queue:
        gid = queue.pop(0)
        if gid in visited:
            continue
        visited.add(gid)
        g = all_genomes_map.get(gid)
        if not g:
            continue
        if g.get("generation", 0) == 0:
            founders.add(g.get("name", gid))
        else:
            for pid in g.get("parents", []):
                if pid:
                    queue.append(pid)
    return sorted(founders)


def pca_2d(data_matrix):
    """Pure Python PCA. Returns 2D coordinates.

    Args:
        data_matrix: list of lists (N samples x M features)
    Returns:
        list of [x, y] coordinates
    """
    n = len(data_matrix)
    m = len(data_matrix[0])

    # 1. Center the data
    means = [0.0] * m
    for row in data_matrix:
        for j in range(m):
            means[j] += row[j]
    means = [x / n for x in means]

    centered = []
    for row in data_matrix:
        centered.append([row[j] - means[j] for j in range(m)])

    # 2. Covariance matrix (M x M)
    cov = [[0.0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i, m):
            s = sum(centered[k][i] * centered[k][j] for k in range(n))
            cov[i][j] = s / (n - 1) if n > 1 else 0
            cov[j][i] = cov[i][j]

    # 3. Power iteration for top 2 eigenvectors
    def power_iteration(matrix, num_iter=200):
        dim = len(matrix)
        import random
        random.seed(42)
        v = [random.gauss(0, 1) for _ in range(dim)]
        norm = math.sqrt(sum(x * x for x in v))
        v = [x / norm for x in v]

        for _ in range(num_iter):
            new_v = [0.0] * dim
            for i in range(dim):
                for j in range(dim):
                    new_v[i] += matrix[i][j] * v[j]
            norm = math.sqrt(sum(x * x for x in new_v))
            if norm < 1e-10:
                break
            v = [x / norm for x in new_v]
            eigenvalue = norm
        return v, eigenvalue

    def deflate(matrix, eigvec, eigval):
        dim = len(matrix)
        new_m = [row[:] for row in matrix]
        for i in range(dim):
            for j in range(dim):
                new_m[i][j] -= eigval * eigvec[i] * eigvec[j]
        return new_m

    # PC1
    pc1_vec, ev1 = power_iteration(cov)
    cov_deflated = deflate(cov, pc1_vec, ev1)

    # PC2
    pc2_vec, ev2 = power_iteration(cov_deflated)

    # 4. Project data
    coords = []
    for row in centered:
        x = sum(row[j] * pc1_vec[j] for j in range(m))
        y = sum(row[j] * pc2_vec[j] for j in range(m))
        coords.append([round(x, 4), round(y, 4)])

    total_var = sum(cov[i][i] for i in range(m))
    var_explained = [
        round(ev1 / total_var * 100, 1) if total_var > 0 else 0,
        round(ev2 / total_var * 100, 1) if total_var > 0 else 0,
    ]

    return coords, var_explained


def main():
    parser = argparse.ArgumentParser(description="PCA of Genomebook agents")
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT))
    parser.add_argument("--max-gen", type=int, default=None)
    args = parser.parse_args()

    genomes = load_genomes(max_gen=args.max_gen)
    print(f"Loaded {len(genomes)} genomes")

    if len(genomes) < 3:
        print("Need at least 3 genomes for PCA")
        return

    # Get consistent trait order
    all_traits = set()
    for g in genomes:
        all_traits.update(g.get("trait_scores", {}).keys())
    trait_names = sorted(all_traits)
    print(f"Traits: {len(trait_names)}")

    # Build data matrix
    data_matrix = []
    for g in genomes:
        traits = g.get("trait_scores", {})
        row = [traits.get(t, 0.5) for t in trait_names]
        data_matrix.append(row)

    # Compute PCA
    coords, var_explained = pca_2d(data_matrix)
    print(f"Variance explained: PC1={var_explained[0]}%, PC2={var_explained[1]}%")

    # Build genome map for lineage tracing
    genome_map = {g["id"]: g for g in genomes}

    # Build point data
    gen_colors = ['#3fb950', '#58a6ff', '#bc8cff', '#e3b341', '#f85149',
                  '#f778ba', '#8b949e', '#39d353', '#d2a8ff', '#ffa657']

    points = []
    for i, g in enumerate(genomes):
        lineage = get_founder_lineage(g, genome_map)
        name = g.get("name", g["id"])
        if len(name) > 45:
            name = name[:42] + "..."

        top3 = sorted(g.get("trait_scores", {}).items(), key=lambda x: x[1], reverse=True)[:3]

        points.append({
            "x": coords[i][0],
            "y": coords[i][1],
            "id": g["id"],
            "name": name,
            "gen": g.get("generation", 0),
            "sex": g.get("sex", "?"),
            "health": round(g.get("health_score", 1.0), 2),
            "lineage": ", ".join(lineage[:3]) + ("..." if len(lineage) > 3 else ""),
            "primary_ancestor": lineage[0] if lineage else "?",
            "top_traits": ", ".join(f"{t.replace('_',' ')}: {s:.2f}" for t, s in top3),
        })

    # Get unique ancestors for coloring
    ancestors = sorted(set(p["primary_ancestor"] for p in points))

    points_json = json.dumps(points, default=str)
    traits_json = json.dumps(trait_names)
    ancestors_json = json.dumps(ancestors)

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genomebook PCA -- Genetic Clustering</title>
<style>
:root {{
  --bg: #0d1117; --bg2: #161b22; --bg3: #21262d;
  --border: #30363d; --text: #e6edf3; --muted: #8b949e;
  --green: #3fb950; --blue: #58a6ff; --red: #f85149;
  --purple: #bc8cff; --orange: #e3b341; --pink: #f778ba;
}}
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ background: var(--bg); color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; height: 100vh; display: flex; flex-direction: column; }}

.header {{
  background: var(--bg2); border-bottom: 1px solid var(--border);
  padding: 0.7rem 1.5rem; display: flex; justify-content: space-between; align-items: center; flex-shrink: 0;
}}
.header h1 {{ font-size: 1.1rem; font-weight: 800; }}
.header h1 span {{ color: var(--green); }}
.header a {{ color: var(--blue); text-decoration: none; font-size: 0.75rem; }}

.controls {{
  background: var(--bg2); border-bottom: 1px solid var(--border);
  padding: 0.5rem 1.5rem; display: flex; gap: 1rem; align-items: center; font-size: 0.78rem; flex-shrink: 0;
}}
.controls label {{ color: var(--muted); }}
.controls select {{
  background: var(--bg3); color: var(--text); border: 1px solid var(--border);
  border-radius: 4px; padding: 0.3rem 0.5rem; font-size: 0.78rem;
}}
.info {{ font-size: 0.75rem; color: var(--muted); }}

.chart-area {{ flex: 1; position: relative; }}
canvas {{ position: absolute; top: 0; left: 0; width: 100%; height: 100%; }}

.tooltip {{
  position: fixed; background: var(--bg2); border: 1px solid var(--green);
  border-radius: 8px; padding: 0.8rem; font-size: 0.78rem; max-width: 320px;
  pointer-events: none; display: none; z-index: 100; line-height: 1.5;
  box-shadow: 0 4px 20px rgba(0,0,0,0.5);
}}
.tt-name {{ font-weight: 700; font-size: 0.88rem; margin-bottom: 0.3rem; }}
.tt-row {{ margin-top: 0.15rem; font-size: 0.75rem; }}
.tt-label {{ color: var(--muted); font-size: 0.68rem; text-transform: uppercase; }}

.legend {{
  position: absolute; bottom: 1rem; right: 1rem; background: var(--bg2);
  border: 1px solid var(--border); border-radius: 8px; padding: 0.6rem 0.8rem;
  font-size: 0.7rem; max-height: 40vh; overflow-y: auto;
}}
.legend-item {{ display: flex; align-items: center; gap: 0.4rem; margin: 0.15rem 0; cursor: pointer; }}
.legend-dot {{ width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }}
.legend-label {{ color: var(--muted); }}
.legend-label:hover {{ color: var(--text); }}
</style>
</head>
<body>

<div class="header">
  <h1>&#x1F9EC; <span>Genomebook</span> PCA</h1>
  <div>
    <span class="info">PC1: {var_explained[0]}% &middot; PC2: {var_explained[1]}% variance explained &middot; {len(points)} agents &middot; {len(trait_names)} traits</span>
    &nbsp;&middot;&nbsp;
    <a href="phylogeny.html">Phylogeny</a> &middot;
    <a href="demo.html">Observatory</a> &middot;
    <a href="https://clawbio.github.io/ClawBio/slides/genomebook/">Slides</a>
  </div>
</div>

<div class="controls">
  <label>Color by:</label>
  <select id="color-mode" onchange="draw()">
    <option value="generation">Generation</option>
    <option value="lineage">Founder Lineage</option>
    <option value="sex">Sex</option>
    <option value="health">Health Score</option>
  </select>
  <label>Size by:</label>
  <select id="size-mode" onchange="draw()">
    <option value="fixed">Fixed</option>
    <option value="health">Health</option>
    <option value="generation">Generation</option>
  </select>
  <label>Highlight gen:</label>
  <select id="gen-filter" onchange="draw()">
    <option value="all">All</option>
    {"".join(f'<option value="{i}">{i}</option>' for i in range(max(p["gen"] for p in points) + 1))}
  </select>
</div>

<div class="chart-area">
  <canvas id="canvas"></canvas>
  <div class="legend" id="legend"></div>
</div>

<div class="tooltip" id="tooltip"></div>

<script>
const POINTS = {points_json};
const TRAITS = {traits_json};
const ANCESTORS = {ancestors_json};

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');

const GEN_COLORS = ['#3fb950','#58a6ff','#bc8cff','#e3b341','#f85149','#f778ba','#8b949e','#39d353','#d2a8ff','#ffa657'];
const LINEAGE_COLORS = {{}};
const lColors = ['#3fb950','#58a6ff','#bc8cff','#e3b341','#f85149','#f778ba','#39d353','#d2a8ff','#ffa657','#8b949e',
                 '#79c0ff','#7ee787','#ffc680','#ff9bce','#b392f0','#56d364','#a371f7','#ffa198','#7ee787','#d2a8ff'];
ANCESTORS.forEach((a, i) => LINEAGE_COLORS[a] = lColors[i % lColors.length]);

function getColor(p, mode) {{
  if (mode === 'generation') return GEN_COLORS[p.gen % GEN_COLORS.length];
  if (mode === 'lineage') return LINEAGE_COLORS[p.primary_ancestor] || '#8b949e';
  if (mode === 'sex') return p.sex === 'Male' ? '#58a6ff' : '#f778ba';
  if (mode === 'health') {{
    if (p.health >= 0.85) return '#3fb950';
    if (p.health >= 0.7) return '#58a6ff';
    if (p.health >= 0.5) return '#e3b341';
    return '#f85149';
  }}
  return '#58a6ff';
}}

function getSize(p, mode) {{
  if (mode === 'health') return 3 + p.health * 6;
  if (mode === 'generation') return 8 - p.gen * 0.5;
  return 5;
}}

function draw() {{
  const W = canvas.width = canvas.offsetWidth * 2;
  const H = canvas.height = canvas.offsetHeight * 2;
  ctx.scale(2, 2);
  const w = W / 2, h = H / 2;
  ctx.clearRect(0, 0, w, h);

  const colorMode = document.getElementById('color-mode').value;
  const sizeMode = document.getElementById('size-mode').value;
  const genFilter = document.getElementById('gen-filter').value;

  // Compute bounds
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (const p of POINTS) {{
    xMin = Math.min(xMin, p.x); xMax = Math.max(xMax, p.x);
    yMin = Math.min(yMin, p.y); yMax = Math.max(yMax, p.y);
  }}
  const xPad = (xMax - xMin) * 0.08 || 1;
  const yPad = (yMax - yMin) * 0.08 || 1;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

  const pad = {{ top: 20, right: 20, bottom: 30, left: 40 }};
  const pw = w - pad.left - pad.right;
  const ph = h - pad.top - pad.bottom;

  function toScreen(px, py) {{
    return [
      pad.left + (px - xMin) / (xMax - xMin) * pw,
      pad.top + (1 - (py - yMin) / (yMax - yMin)) * ph
    ];
  }}

  // Grid
  ctx.strokeStyle = '#21262d'; ctx.lineWidth = 0.5;
  for (let i = 0; i <= 4; i++) {{
    const y = pad.top + ph * i / 4;
    ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(w - pad.right, y); ctx.stroke();
    const x = pad.left + pw * i / 4;
    ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, pad.top + ph); ctx.stroke();
  }}

  // Axis labels
  ctx.fillStyle = '#8b949e'; ctx.font = '10px system-ui';
  ctx.textAlign = 'center';
  ctx.fillText('PC1 ({var_explained[0]}%)', w / 2, h - 5);
  ctx.save(); ctx.translate(12, h / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('PC2 ({var_explained[1]}%)', 0, 0);
  ctx.restore();

  // Draw points (dimmed first, then highlighted)
  const layers = [[], []];
  for (const p of POINTS) {{
    const active = genFilter === 'all' || p.gen === parseInt(genFilter);
    layers[active ? 1 : 0].push(p);
  }}

  for (let layer = 0; layer < 2; layer++) {{
    for (const p of layers[layer]) {{
      const [sx, sy] = toScreen(p.x, p.y);
      const r = getSize(p, sizeMode);
      const color = getColor(p, colorMode);
      const alpha = layer === 0 ? 0.15 : 0.85;

      ctx.beginPath();
      ctx.arc(sx, sy, r, 0, Math.PI * 2);
      ctx.fillStyle = color.replace(')', ',' + alpha + ')').replace('rgb', 'rgba').replace('#', '');

      // Convert hex to rgba
      const hex = color.replace('#', '');
      const cr = parseInt(hex.substr(0, 2), 16);
      const cg = parseInt(hex.substr(2, 2), 16);
      const cb = parseInt(hex.substr(4, 2), 16);
      ctx.fillStyle = 'rgba(' + cr + ',' + cg + ',' + cb + ',' + alpha + ')';

      ctx.fill();

      // Founder labels
      if (p.gen === 0 && layer === 1) {{
        ctx.fillStyle = 'rgba(' + cr + ',' + cg + ',' + cb + ',0.9)';
        ctx.font = '8px system-ui';
        ctx.textAlign = 'left';
        const label = p.name.split(' ').pop();
        ctx.fillText(label, sx + r + 3, sy + 3);
      }}
    }}
  }}

  ctx.setTransform(1, 0, 0, 1, 0, 0);

  // Legend
  renderLegend(colorMode);
}}

function renderLegend(mode) {{
  const legend = document.getElementById('legend');
  let items = [];
  if (mode === 'generation') {{
    const gens = [...new Set(POINTS.map(p => p.gen))].sort((a,b) => a - b);
    items = gens.map(g => ({{ color: GEN_COLORS[g % GEN_COLORS.length], label: 'Gen ' + g, count: POINTS.filter(p => p.gen === g).length }}));
  }} else if (mode === 'lineage') {{
    items = ANCESTORS.map(a => ({{ color: LINEAGE_COLORS[a], label: a, count: POINTS.filter(p => p.primary_ancestor === a).length }}));
  }} else if (mode === 'sex') {{
    items = [{{ color: '#58a6ff', label: 'Male', count: POINTS.filter(p => p.sex === 'Male').length }},
             {{ color: '#f778ba', label: 'Female', count: POINTS.filter(p => p.sex === 'Female').length }}];
  }} else {{
    items = [{{ color: '#3fb950', label: '>= 0.85' }}, {{ color: '#58a6ff', label: '0.7-0.85' }},
             {{ color: '#e3b341', label: '0.5-0.7' }}, {{ color: '#f85149', label: '< 0.5' }}];
  }}
  legend.innerHTML = items.map(it =>
    '<div class="legend-item"><div class="legend-dot" style="background:' + it.color + '"></div>' +
    '<span class="legend-label">' + it.label + (it.count ? ' (' + it.count + ')' : '') + '</span></div>'
  ).join('');
}}

// Tooltip on hover
canvas.addEventListener('mousemove', (e) => {{
  const rect = canvas.getBoundingClientRect();
  const mx = e.clientX - rect.left;
  const my = e.clientY - rect.top;

  // Convert to data space
  const w = rect.width, h = rect.height;
  const pad = {{ top: 20, right: 20, bottom: 30, left: 40 }};
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (const p of POINTS) {{
    xMin = Math.min(xMin, p.x); xMax = Math.max(xMax, p.x);
    yMin = Math.min(yMin, p.y); yMax = Math.max(yMax, p.y);
  }}
  const xPad = (xMax - xMin) * 0.08 || 1;
  const yPad = (yMax - yMin) * 0.08 || 1;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

  const pw = w - pad.left - pad.right;
  const ph = h - pad.top - pad.bottom;

  let closest = null, minDist = 15;
  for (const p of POINTS) {{
    const sx = pad.left + (p.x - xMin) / (xMax - xMin) * pw;
    const sy = pad.top + (1 - (p.y - yMin) / (yMax - yMin)) * ph;
    const d = Math.sqrt((mx - sx) ** 2 + (my - sy) ** 2);
    if (d < minDist) {{ minDist = d; closest = p; }}
  }}

  if (closest) {{
    tooltip.style.display = 'block';
    tooltip.style.left = Math.min(e.clientX + 15, window.innerWidth - 340) + 'px';
    tooltip.style.top = Math.min(e.clientY + 15, window.innerHeight - 180) + 'px';
    const hc = closest.health >= 0.7 ? 'color:#3fb950' : closest.health >= 0.5 ? 'color:#e3b341' : 'color:#f85149';
    tooltip.innerHTML =
      '<div class="tt-name">' + closest.name + '</div>' +
      '<div class="tt-row"><span class="tt-label">Gen: </span>' + closest.gen + ' &middot; ' + closest.sex +
      ' &middot; <span style="' + hc + '">Health: ' + closest.health + '</span></div>' +
      '<div class="tt-row"><span class="tt-label">Lineage: </span>' + closest.lineage + '</div>' +
      '<div class="tt-row"><span class="tt-label">Top traits: </span>' + closest.top_traits + '</div>';
  }} else {{
    tooltip.style.display = 'none';
  }}
}});

canvas.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; }});
window.addEventListener('resize', draw);
draw();
</script>
</body>
</html>'''

    out = Path(args.output)
    out.write_text(html)
    print(f"Written: {out}")
    print(f"Size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
