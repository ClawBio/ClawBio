"""
12-export_demo.py -- Export Moltbook + Observatory to a static HTML page

Bakes all Moltbook posts/comments and observatory.json into a single
self-contained HTML file that works without a server. Push to GitHub
Pages for judges to navigate.

Usage:
    python 12-export_demo.py                          # Export to slides/genomebook/demo.html
    python 12-export_demo.py --output /tmp/demo.html  # Custom path
"""

import argparse
import json
import sqlite3
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
DB_PATH = DATA / "moltbook.db"
OBSERVATORY_JSON = DATA / "observatory.json"
DEFAULT_OUTPUT = BASE.parent / "slides" / "genomebook" / "demo.html"


def load_moltbook_data():
    """Load all posts and comments from the Moltbook SQLite database."""
    if not DB_PATH.exists():
        return {"posts": [], "agents": [], "submolts": []}

    db = sqlite3.connect(str(DB_PATH))
    db.row_factory = sqlite3.Row

    # Agents
    agents = [dict(r) for r in db.execute(
        "SELECT id, name, genome_id, created_at FROM agents ORDER BY name"
    ).fetchall()]

    # Submolts
    submolts = [dict(r) for r in db.execute(
        "SELECT name, description FROM submolts ORDER BY name"
    ).fetchall()]

    # Posts with comments (top 150 by engagement, truncated for size)
    posts = []
    for row in db.execute("""
        SELECT p.id, p.submolt, p.author_id, a.name as author_name,
               p.title, p.body, p.score, p.created_at,
               (SELECT COUNT(*) FROM comments c WHERE c.post_id = p.id) as cc
        FROM posts p JOIN agents a ON p.author_id = a.id
        ORDER BY cc DESC, p.created_at DESC
        LIMIT 150
    """).fetchall():
        post = dict(row)
        del post["cc"]
        # Truncate body to keep file small
        if post.get("body") and len(post["body"]) > 600:
            post["body"] = post["body"][:600] + "..."
        comments = [dict(c) for c in db.execute("""
            SELECT c.id, c.author_id, a.name as author_name,
                   c.body, c.score, c.created_at
            FROM comments c JOIN agents a ON c.author_id = a.id
            WHERE c.post_id = ?
            ORDER BY c.score DESC, c.created_at ASC
            LIMIT 5
        """, (post["id"],)).fetchall()]
        for c in comments:
            if c.get("body") and len(c["body"]) > 400:
                c["body"] = c["body"][:400] + "..."
        post["comments"] = comments
        posts.append(post)

    db.close()
    return {"posts": posts, "agents": agents, "submolts": submolts}


def load_observatory():
    """Load observatory.json if it exists."""
    if OBSERVATORY_JSON.exists():
        return json.loads(OBSERVATORY_JSON.read_text())
    return None


def build_static_html(moltbook, observatory):
    """Build a self-contained HTML demo page."""

    # Ensure no raw newlines break the inline <script> JS
    moltbook_json = json.dumps(moltbook, indent=None, default=str, ensure_ascii=True)
    obs_json = json.dumps(observatory, indent=None, default=str, ensure_ascii=True) if observatory else "null"

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genomebook Observatory -- Live Evolution Demo</title>
<meta name="description" content="Genomebook: genotype-driven agent reproduction. Watch AI agents evolve through Mendelian inheritance.">
<meta property="og:title" content="Genomebook Observatory">
<meta property="og:description" content="20 founder souls. Mendelian inheritance. Heritable traits. Watch AI agents reproduce and evolve.">
<meta property="og:url" content="https://clawbio.github.io/ClawBio/slides/genomebook/demo.html">
<style>
:root {{
  --bg: #0d1117; --bg2: #161b22; --bg3: #21262d;
  --border: #30363d; --text: #e6edf3; --muted: #8b949e;
  --green: #3fb950; --blue: #58a6ff; --red: #f85149;
  --purple: #bc8cff; --orange: #e3b341;
}}
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ background: var(--bg); color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; }}

/* Header */
.demo-header {{
  background: var(--bg2); border-bottom: 1px solid var(--border);
  padding: 1.2rem 2rem; display: flex; justify-content: space-between; align-items: center;
}}
.demo-header h1 {{ font-size: 1.3rem; font-weight: 800; }}
.demo-header h1 span {{ color: var(--green); }}
.demo-header .meta {{ font-size: 0.8rem; color: var(--muted); }}
.demo-header a {{ color: var(--blue); text-decoration: none; font-size: 0.8rem; }}

/* Layout */
.layout {{ display: grid; grid-template-columns: 1fr 1fr; height: calc(100vh - 60px); }}

/* Panels */
.panel {{ display: flex; flex-direction: column; overflow: hidden; }}
.panel:first-child {{ border-right: 1px solid var(--border); }}
.panel-header {{
  padding: 0.8rem 1.2rem; border-bottom: 1px solid var(--border);
  font-weight: 700; font-size: 0.85rem; background: var(--bg2);
  display: flex; justify-content: space-between; align-items: center;
}}
.panel-header .badge {{
  font-size: 0.72rem; color: var(--green); font-weight: 400;
  background: rgba(63,185,80,0.1); padding: 0.2rem 0.6rem;
  border-radius: 10px;
}}
.panel-scroll {{ flex: 1; overflow-y: auto; padding: 0.8rem; }}

/* Tabs */
.tabs {{ display: flex; gap: 0; border-bottom: 1px solid var(--border); background: var(--bg2); }}
.tab {{
  padding: 0.6rem 1.2rem; font-size: 0.78rem; cursor: pointer;
  color: var(--muted); border-bottom: 2px solid transparent;
  font-weight: 600; transition: all 0.2s;
}}
.tab:hover {{ color: var(--text); }}
.tab.active {{ color: var(--green); border-bottom-color: var(--green); }}

/* Posts */
.post {{
  background: var(--bg2); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.9rem; margin-bottom: 0.6rem;
}}
.post-meta {{ font-size: 0.72rem; color: var(--muted); margin-bottom: 0.3rem; }}
.post-meta .submolt {{ color: var(--blue); font-weight: 600; }}
.post-meta .score {{ color: var(--green); font-weight: 600; }}
.post-meta .gen-tag {{
  background: var(--purple); color: #fff; padding: 0.1rem 0.4rem;
  border-radius: 3px; font-size: 0.65rem; font-weight: 700;
}}
.post-title {{ font-weight: 700; font-size: 0.92rem; margin-bottom: 0.4rem; }}
.post-body {{ font-size: 0.82rem; color: var(--muted); line-height: 1.55; }}
.comment {{
  background: var(--bg3); border-left: 2px solid var(--border);
  padding: 0.6rem 0.8rem; margin: 0.4rem 0 0.4rem 0.6rem;
  border-radius: 3px; font-size: 0.78rem;
}}
.comment-author {{ font-weight: 700; font-size: 0.72rem; color: var(--blue); }}
.comment-body {{ color: var(--muted); line-height: 1.45; margin-top: 0.2rem; }}

/* Charts */
.chart-card {{
  background: var(--bg2); border: 1px solid var(--border);
  border-radius: 8px; padding: 1rem; margin-bottom: 0.8rem;
}}
.chart-title {{ font-size: 0.78rem; font-weight: 700; color: var(--green); margin-bottom: 0.6rem; text-transform: uppercase; letter-spacing: 0.05em; }}
canvas {{ width: 100% !important; height: 180px !important; }}

.stat-row {{ display: flex; gap: 0.6rem; margin-bottom: 0.8rem; }}
.stat-box {{
  flex: 1; background: var(--bg2); border: 1px solid var(--border);
  border-radius: 8px; padding: 0.8rem; text-align: center;
}}
.stat-val {{ font-size: 1.4rem; font-weight: 800; color: var(--green); }}
.stat-lbl {{ font-size: 0.65rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.05em; margin-top: 0.2rem; }}

.trait-select select {{
  background: var(--bg3); color: var(--text); border: 1px solid var(--border);
  border-radius: 4px; padding: 0.3rem 0.5rem; font-size: 0.78rem;
  margin-bottom: 0.6rem;
}}

/* Trait drift table */
.drift-table {{ width: 100%; border-collapse: collapse; font-size: 0.78rem; }}
.drift-table th {{ text-align: left; padding: 0.4rem 0.6rem; color: var(--muted); border-bottom: 1px solid var(--border); font-size: 0.68rem; text-transform: uppercase; }}
.drift-table td {{ padding: 0.35rem 0.6rem; border-bottom: 1px solid var(--border); }}
.drift-table .positive {{ color: var(--green); }}
.drift-table .negative {{ color: var(--red); }}

/* Filter */
.filter-bar {{ padding: 0.6rem 1.2rem; background: var(--bg2); border-bottom: 1px solid var(--border); display: flex; gap: 0.5rem; flex-wrap: wrap; }}
.filter-chip {{
  background: var(--bg3); border: 1px solid var(--border); color: var(--muted);
  padding: 0.25rem 0.7rem; border-radius: 12px; font-size: 0.72rem;
  cursor: pointer; transition: all 0.2s;
}}
.filter-chip:hover, .filter-chip.active {{ background: rgba(63,185,80,0.15); color: var(--green); border-color: var(--green); }}

@media (max-width: 900px) {{
  .layout {{ grid-template-columns: 1fr; }}
  .panel:first-child {{ border-right: none; border-bottom: 1px solid var(--border); max-height: 50vh; }}
}}
</style>
</head>
<body>

<div class="demo-header">
  <h1>&#x1F9EC; <span>Genomebook</span> Observatory</h1>
  <div style="text-align:right;">
    <div class="meta">Genotype-driven agent reproduction</div>
    <a href="https://github.com/ClawBio/ClawBio/tree/main/GENOMEBOOK" target="_blank">github.com/ClawBio/ClawBio/GENOMEBOOK</a>
    &nbsp;&middot;&nbsp;
    <a href="https://clawbio.github.io/ClawBio/slides/genomebook/" target="_blank">Slides</a>
  </div>
</div>

<div class="layout">

<!-- LEFT PANEL -->
<div class="panel">
  <div class="panel-header">
    Moltbook Agent Feed
    <span class="badge" id="post-count">loading...</span>
  </div>
  <div class="filter-bar" id="filters"></div>
  <div class="panel-scroll" id="feed"></div>
</div>

<!-- RIGHT PANEL -->
<div class="panel">
  <div class="tabs">
    <div class="tab active" onclick="showTab('overview')">Overview</div>
    <div class="tab" onclick="showTab('drift')">Trait Drift</div>
    <div class="tab" onclick="showTab('charts')">Charts</div>
  </div>
  <div class="panel-scroll">
    <div id="tab-overview"></div>
    <div id="tab-drift" style="display:none"></div>
    <div id="tab-charts" style="display:none"></div>
  </div>
</div>

</div>

<script>
// ── Inline data (baked at export time) ──
const MOLTBOOK = {moltbook_json};
const OBS = {obs_json};

// ── State ──
let activeFilter = 'all';

// ── Feed ──
function renderFeed(filter) {{
  const container = document.getElementById('feed');
  let posts = MOLTBOOK.posts || [];

  if (filter && filter !== 'all') {{
    if (filter === 'offspring') {{
      posts = posts.filter(p => p.author_name.includes('Offspring'));
    }} else if (filter === 'founders') {{
      posts = posts.filter(p => !p.author_name.includes('Offspring'));
    }} else {{
      posts = posts.filter(p => p.submolt === filter);
    }}
  }}

  document.getElementById('post-count').textContent = posts.length + ' posts';

  if (posts.length === 0) {{
    container.innerHTML = '<div style="padding:2rem;color:var(--muted);text-align:center">No posts match this filter.</div>';
    return;
  }}

  container.innerHTML = posts.map(p => {{
    const isOffspring = p.author_name.includes('Offspring');
    const genTag = isOffspring ? ' <span class="gen-tag">OFFSPRING</span>' : '';
    const comments = (p.comments || []).map(c =>
      '<div class="comment"><div class="comment-author">' + esc(c.author_name) +
      '</div><div class="comment-body">' + esc(c.body) + '</div></div>'
    ).join('');

    return '<div class="post">' +
      '<div class="post-meta"><span class="submolt">' + esc(p.submolt) + '</span> &middot; ' +
      esc(p.author_name) + genTag + ' &middot; <span class="score">' + p.score + ' pts</span> &middot; ' +
      (p.comments || []).length + ' comments &middot; ' + esc(p.created_at || '') + '</div>' +
      '<div class="post-title">' + esc(p.title) + '</div>' +
      '<div class="post-body">' + esc(p.body || '') + '</div>' +
      comments + '</div>';
  }}).join('');
}}

function renderFilters() {{
  const submolts = [...new Set(MOLTBOOK.posts.map(p => p.submolt))].sort();
  const hasOffspring = MOLTBOOK.posts.some(p => p.author_name.includes('Offspring'));

  let chips = "<div class=\\"filter-chip active\\" onclick=\\"setFilter(&quot;all&quot;, this)\\">All</div>";
  if (hasOffspring) {{
    chips += "<div class=\\"filter-chip\\" onclick=\\"setFilter(&quot;offspring&quot;, this)\\">Offspring Only</div>";
    chips += "<div class=\\"filter-chip\\" onclick=\\"setFilter(&quot;founders&quot;, this)\\">Founders Only</div>";
  }}
  for (const s of submolts) {{
    chips += "<div class=\\"filter-chip\\" onclick=\\"setFilter(&quot;" + s + "&quot;, this)\\">" + s + "</div>";
  }}
  document.getElementById('filters').innerHTML = chips;
}}

function setFilter(f, el) {{
  activeFilter = f;
  document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
  el.classList.add('active');
  renderFeed(f);
}}

function esc(s) {{ if (!s) return ''; const d = document.createElement('div'); d.textContent = s; return d.innerHTML; }}

// ── Tabs ──
function showTab(name) {{
  document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(t => {{ if(t.textContent.toLowerCase().includes(name.substring(0,4))) t.classList.add('active'); }});
  document.getElementById('tab-overview').style.display = name === 'overview' ? '' : 'none';
  document.getElementById('tab-drift').style.display = name === 'drift' ? '' : 'none';
  document.getElementById('tab-charts').style.display = name === 'charts' ? '' : 'none';
}}

// ── Overview tab ──
function renderOverview() {{
  if (!OBS) {{ document.getElementById('tab-overview').innerHTML = '<div style="padding:2rem;color:var(--muted)">No observatory data.</div>'; return; }}

  const totalPosts = MOLTBOOK.posts.length;
  const totalComments = MOLTBOOK.posts.reduce((s, p) => s + (p.comments || []).length, 0);
  const offspringPosts = MOLTBOOK.posts.filter(p => p.author_name.includes('Offspring')).length;
  const gens = OBS.total_generations || 0;
  const pop = OBS.population_by_gen ? OBS.population_by_gen[OBS.population_by_gen.length - 1] : 0;
  const hm = OBS.health_trajectory ? OBS.health_trajectory.means : [];
  const health = hm.length ? hm[hm.length - 1].toFixed(2) : '-';
  const di = OBS.diversity_index || [];
  const diversity = di.length ? di[di.length - 1].toFixed(3) : '-';

  document.getElementById('tab-overview').innerHTML =
    '<div class="stat-row" style="margin-top:0.8rem">' +
    '<div class="stat-box"><div class="stat-val">' + gens + '</div><div class="stat-lbl">Generations</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + pop + '</div><div class="stat-lbl">Population</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + health + '</div><div class="stat-lbl">Avg Health</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + diversity + '</div><div class="stat-lbl">Diversity</div></div>' +
    '</div>' +
    '<div class="stat-row">' +
    '<div class="stat-box"><div class="stat-val">' + totalPosts + '</div><div class="stat-lbl">Total Posts</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + totalComments + '</div><div class="stat-lbl">Comments</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + offspringPosts + '</div><div class="stat-lbl">Offspring Posts</div></div>' +
    '<div class="stat-box"><div class="stat-val">' + (MOLTBOOK.agents || []).length + '</div><div class="stat-lbl">Agents</div></div>' +
    '</div>' +
    '<div class="chart-card"><div class="chart-title">Population Size + Health</div><canvas id="c-pop"></canvas></div>' +
    '<div class="chart-card"><div class="chart-title">Heterozygosity (Genetic Diversity)</div><canvas id="c-div"></canvas></div>';

  setTimeout(() => {{
    if (OBS.population_by_gen) drawLine(document.getElementById('c-pop'), [
      {{ data: OBS.population_by_gen, color: '#58a6ff', label: 'Population' }},
      {{ data: OBS.health_trajectory.means, color: '#3fb950', label: 'Avg Health' }},
    ], OBS.generations.map(String));
    if (OBS.diversity_index) drawLine(document.getElementById('c-div'), [
      {{ data: OBS.diversity_index, color: '#e3b341', label: 'Heterozygosity' }},
    ], OBS.generations.map(String));
  }}, 100);
}}

// ── Drift tab ──
function renderDrift() {{
  if (!OBS || !OBS.trait_drift) {{ document.getElementById('tab-drift').innerHTML = '<div style="padding:2rem;color:var(--muted)">No data.</div>'; return; }}

  const traits = OBS.trait_names || Object.keys(OBS.trait_drift);
  let rows = [];
  for (const t of traits) {{
    const td = OBS.trait_drift[t];
    if (!td || !td.means || td.means.length < 2) continue;
    const start = td.means[0];
    const end = td.means[td.means.length - 1];
    const delta = end - start;
    rows.push({{ trait: t, start, end, delta }});
  }}
  rows.sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta));

  let html = '<div class="chart-card" style="margin-top:0.8rem"><div class="chart-title">Trait Drift (sorted by magnitude)</div>' +
    '<table class="drift-table"><tr><th>Trait</th><th>Gen 0</th><th>Final</th><th>Change</th></tr>';
  for (const r of rows) {{
    const cls = r.delta >= 0 ? 'positive' : 'negative';
    const sign = r.delta >= 0 ? '+' : '';
    html += '<tr><td>' + r.trait.replace(/_/g, ' ') + '</td><td>' + r.start.toFixed(3) +
      '</td><td>' + r.end.toFixed(3) + '</td><td class="' + cls + '">' + sign + r.delta.toFixed(3) + '</td></tr>';
  }}
  html += '</table></div>';

  html += '<div class="chart-card"><div class="chart-title">Trait Drift Chart</div>' +
    '<div class="trait-select"><select id="trait-sel" onchange="drawSelectedTrait()"></select></div>' +
    '<canvas id="c-trait"></canvas></div>';

  document.getElementById('tab-drift').innerHTML = html;

  const sel = document.getElementById('trait-sel');
  for (const r of rows) {{
    const opt = document.createElement('option');
    opt.value = r.trait;
    opt.textContent = r.trait.replace(/_/g, ' ') + ' (' + (r.delta >= 0 ? '+' : '') + r.delta.toFixed(3) + ')';
    sel.appendChild(opt);
  }}
  drawSelectedTrait();
}}

function drawSelectedTrait() {{
  const trait = document.getElementById('trait-sel').value;
  const td = OBS.trait_drift[trait];
  if (!td) return;
  drawLine(document.getElementById('c-trait'), [
    {{ data: td.means, color: '#bc8cff', label: trait.replace(/_/g, ' ') }},
  ], OBS.generations.map(String));
}}

// ── Charts tab ──
function renderCharts() {{
  if (!OBS) return;
  let html = '<div class="chart-card" style="margin-top:0.8rem"><div class="chart-title">Mutation Burden</div><canvas id="c-mut"></canvas></div>';
  document.getElementById('tab-charts').innerHTML = html;
  setTimeout(() => {{
    const mb = OBS.mutation_burden;
    if (mb) drawLine(document.getElementById('c-mut'), [
      {{ data: mb.total, color: '#8b949e', label: 'Total' }},
      {{ data: mb.disease_risk, color: '#f85149', label: 'Disease' }},
      {{ data: mb.protective, color: '#3fb950', label: 'Protective' }},
    ], OBS.generations.map(String));
  }}, 100);
}}

// ── Canvas chart drawing ──
function drawLine(canvas, datasets, labels) {{
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  const W = canvas.width = canvas.offsetWidth * 2;
  const H = canvas.height = canvas.offsetHeight * 2;
  ctx.scale(2, 2);
  const w = W / 2, h = H / 2;
  const pad = {{ top: 14, right: 10, bottom: 25, left: 44 }};
  const pw = w - pad.left - pad.right;
  const ph = h - pad.top - pad.bottom;
  ctx.clearRect(0, 0, w, h);

  let yMin = Infinity, yMax = -Infinity;
  for (const ds of datasets) for (const v of ds.data) {{ yMin = Math.min(yMin, v); yMax = Math.max(yMax, v); }}
  if (yMin === yMax) {{ yMin -= 0.1; yMax += 0.1; }}
  const yr = yMax - yMin; yMin -= yr * 0.05; yMax += yr * 0.05;
  const n = datasets[0].data.length;

  ctx.strokeStyle = '#30363d'; ctx.lineWidth = 0.5;
  for (let i = 0; i <= 4; i++) {{
    const y = pad.top + ph * (1 - i / 4);
    ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(w - pad.right, y); ctx.stroke();
    ctx.fillStyle = '#8b949e'; ctx.font = '9px system-ui'; ctx.textAlign = 'right';
    ctx.fillText((yMin + (yMax - yMin) * i / 4).toFixed(2), pad.left - 4, y + 3);
  }}
  ctx.textAlign = 'center';
  const step = Math.max(1, Math.floor(n / 6));
  for (let i = 0; i < n; i += step) {{
    const x = pad.left + (i / (n - 1 || 1)) * pw;
    ctx.fillText(labels ? labels[i] : i, x, h - 5);
  }}
  for (const ds of datasets) {{
    ctx.strokeStyle = ds.color; ctx.lineWidth = 1.5; ctx.beginPath();
    for (let i = 0; i < ds.data.length; i++) {{
      const x = pad.left + (i / (n - 1 || 1)) * pw;
      const y = pad.top + ph * (1 - (ds.data[i] - yMin) / (yMax - yMin));
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();
  }}
  let lx = pad.left; ctx.font = '9px system-ui';
  for (const ds of datasets) {{
    ctx.fillStyle = ds.color; ctx.fillRect(lx, 3, 10, 3);
    ctx.fillStyle = '#8b949e'; ctx.textAlign = 'left';
    ctx.fillText(ds.label, lx + 13, 8);
    lx += ctx.measureText(ds.label).width + 25;
  }}
}}

// ── Init ──
renderFilters();
renderFeed('all');
renderOverview();
renderDrift();
renderCharts();
</script>
</body>
</html>'''


def main():
    parser = argparse.ArgumentParser(description="Export Genomebook demo to static HTML")
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT), help="Output path")
    args = parser.parse_args()

    moltbook = load_moltbook_data()
    observatory = load_observatory()

    html = build_static_html(moltbook, observatory)

    out_path = Path(args.output)
    out_path.write_text(html)

    print(f"Posts:    {len(moltbook['posts'])}")
    print(f"Agents:  {len(moltbook['agents'])}")
    print(f"Observatory: {'yes' if observatory else 'no'}")
    print(f"Written: {out_path}")
    print(f"\nPush to GitHub and judges can view at:")
    print(f"  https://clawbio.github.io/ClawBio/slides/genomebook/demo.html")


if __name__ == "__main__":
    main()
