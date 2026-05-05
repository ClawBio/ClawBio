#!/usr/bin/env python3
"""
analyze_fasta.py - Analisis general de archivos FASTA.
Detecta automaticamente si es nucleotido o proteina y genera un reporte completo.

Uso:
    python3 analyze_fasta.py archivo.fasta
    python3 analyze_fasta.py archivo.fasta --json
    python3 analyze_fasta.py archivo.fasta --html reporte.html
"""

import sys
import json
import html as html_mod
import argparse
from pathlib import Path
from collections import Counter
from datetime import datetime

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ──────────────────────────────────────────────
# Deteccion de tipo
# ──────────────────────────────────────────────

def detect_sequence_type(seq_str):
    """Detecta si la secuencia es nucleotido o proteina."""
    nuc_chars = set("ATCGUNatcgun")
    sample = seq_str[:500].replace("-", "").replace(".", "")
    nuc_ratio = sum(1 for c in sample if c in nuc_chars) / max(len(sample), 1)
    return "nucleotide" if nuc_ratio > 0.85 else "protein"


# ──────────────────────────────────────────────
# Analisis de nucleotidos
# ──────────────────────────────────────────────

def analyze_nucleotide(record):
    seq = record.seq
    seq_str = str(seq).upper()
    length = len(seq)
    comp = Counter(seq_str)

    result = {
        "id": record.id,
        "description": record.description,
        "length_bp": length,
        "gc_content": round(gc_fraction(seq) * 100, 2),
        "composition": {
            "A": comp.get("A", 0),
            "T": comp.get("T", 0),
            "G": comp.get("G", 0),
            "C": comp.get("C", 0),
            "N": comp.get("N", 0),
        },
        "composition_pct": {
            "A": round(comp.get("A", 0) / max(length, 1) * 100, 2),
            "T": round(comp.get("T", 0) / max(length, 1) * 100, 2),
            "G": round(comp.get("G", 0) / max(length, 1) * 100, 2),
            "C": round(comp.get("C", 0) / max(length, 1) * 100, 2),
        },
        "at_content": round((comp.get("A", 0) + comp.get("T", 0)) / max(length, 1) * 100, 2),
    }

    # Dinucleotidos
    dinuc = Counter(seq_str[i:i+2] for i in range(len(seq_str) - 1))
    top_dinuc = dinuc.most_common(5)
    result["top_dinucleotides"] = {d: c for d, c in top_dinuc}

    # Buscar ORFs simples (ATG ... stop)
    stops = {"TAA", "TAG", "TGA"}
    orfs = []
    for frame in range(3):
        i = frame
        while i < length - 2:
            codon = seq_str[i:i+3]
            if codon == "ATG":
                start = i
                j = i + 3
                while j < length - 2:
                    c = seq_str[j:j+3]
                    if c in stops:
                        orf_len = j - start + 3
                        if orf_len >= 300:
                            orfs.append({
                                "frame": f"+{frame+1}",
                                "start": start + 1,
                                "end": j + 3,
                                "length_bp": orf_len,
                                "length_aa": orf_len // 3,
                            })
                        break
                    j += 3
                i = j + 3
            else:
                i += 3

    result["orfs_found"] = len(orfs)
    result["orfs"] = orfs[:10]

    try:
        mw = molecular_weight(seq, "DNA")
        result["molecular_weight_da"] = round(mw, 1)
    except Exception:
        pass

    return result


# ──────────────────────────────────────────────
# Analisis de proteinas
# ──────────────────────────────────────────────

def analyze_protein(record):
    seq_str = str(record.seq).upper().replace("X", "").replace("*", "")
    length = len(seq_str)
    comp = Counter(seq_str)

    result = {
        "id": record.id,
        "description": record.description,
        "length_aa": length,
        "composition": {aa: comp.get(aa, 0) for aa in sorted(comp.keys())},
        "composition_pct": {
            aa: round(count / max(length, 1) * 100, 2)
            for aa, count in sorted(comp.items())
        },
    }

    try:
        pa = ProteinAnalysis(seq_str)
        result["molecular_weight_da"] = round(pa.molecular_weight(), 1)
        result["isoelectric_point"] = round(pa.isoelectric_point(), 2)
        result["aromaticity"] = round(pa.aromaticity(), 4)
        result["instability_index"] = round(pa.instability_index(), 2)
        result["stability"] = "estable" if pa.instability_index() < 40 else "inestable"
        result["gravy"] = round(pa.gravy(), 4)
        result["hydrophobicity"] = "hidrofilica" if pa.gravy() < 0 else "hidrofobica"

        helix, turn, sheet = pa.secondary_structure_fraction()
        result["secondary_structure_pct"] = {
            "helix": round(helix * 100, 1),
            "turn": round(turn * 100, 1),
            "sheet": round(sheet * 100, 1),
        }

        aa_pct = pa.amino_acids_percent
        charged = aa_pct.get("R", 0) + aa_pct.get("K", 0) + aa_pct.get("D", 0) + aa_pct.get("E", 0)
        result["charged_residues_pct"] = round(charged * 100, 1)

        aromatic = aa_pct.get("F", 0) + aa_pct.get("W", 0) + aa_pct.get("Y", 0)
        result["aromatic_residues_pct"] = round(aromatic * 100, 1)

    except Exception as e:
        result["analysis_error"] = str(e)

    return result


# ──────────────────────────────────────────────
# Analisis principal
# ──────────────────────────────────────────────

def analyze_fasta(filepath):
    path = Path(filepath)
    if not path.exists():
        return {"error": f"Archivo no encontrado: {filepath}"}

    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        return {"error": "No se encontraron secuencias en el archivo FASTA"}

    # Loud-failure validations (no silent degradation)
    first_seq = str(records[0].seq)
    if len(first_seq) < 10:
        return {"error": f"Secuencia demasiado corta ({len(first_seq)} caracteres). Minimo 10."}
    n_ratio = first_seq.upper().count("N") / max(len(first_seq), 1)
    if n_ratio > 0.5:
        return {"error": f"Secuencia con {round(n_ratio*100,1)}% de Ns - datos insuficientes para analisis confiable."}

    seq_type = detect_sequence_type(first_seq)

    report = {
        "file": path.name,
        "file_path": str(path.resolve()),
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_sequences": len(records),
        "sequence_type": seq_type,
        "sequences": [],
        "summary": {},
    }

    lengths = []
    for record in records:
        if seq_type == "nucleotide":
            analysis = analyze_nucleotide(record)
        else:
            analysis = analyze_protein(record)
        report["sequences"].append(analysis)
        lengths.append(len(record.seq))

    report["summary"] = {
        "total_sequences": len(records),
        "total_residues": sum(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "avg_length": round(sum(lengths) / len(lengths), 1),
        "n50": calculate_n50(lengths),
    }

    if seq_type == "nucleotide":
        gc_values = [s.get("gc_content", 0) for s in report["sequences"]]
        report["summary"]["avg_gc_content"] = round(sum(gc_values) / len(gc_values), 2)
        report["summary"]["total_orfs"] = sum(s.get("orfs_found", 0) for s in report["sequences"])

    return report


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running = 0
    for length in sorted_lengths:
        running += length
        if running >= total / 2:
            return length
    return 0


# ──────────────────────────────────────────────
# Formato texto
# ──────────────────────────────────────────────

def format_text_report(report):
    lines = []
    lines.append("=" * 60)
    lines.append(f"  ANALISIS FASTA: {report['file']}")
    lines.append("=" * 60)
    lines.append("")

    s = report["summary"]
    lines.append(f"Tipo de secuencia:   {report['sequence_type']}")
    lines.append(f"Total secuencias:    {s['total_sequences']}")
    lines.append(f"Total residuos:      {s['total_residues']:,}")
    lines.append(f"Longitud min:        {s['min_length']:,}")
    lines.append(f"Longitud max:        {s['max_length']:,}")
    lines.append(f"Longitud promedio:   {s['avg_length']:,}")
    lines.append(f"N50:                 {s['n50']:,}")

    if report["sequence_type"] == "nucleotide":
        lines.append(f"GC promedio:         {s.get('avg_gc_content', 'N/A')}%")
        lines.append(f"ORFs encontrados:    {s.get('total_orfs', 0)}")

    lines.append("")
    lines.append("-" * 60)

    for i, seq in enumerate(report["sequences"][:20]):
        lines.append("")
        lines.append(f"  Secuencia {i+1}: {seq.get('id', 'unknown')}")
        lines.append(f"  {seq.get('description', '')}")
        lines.append("")

        if report["sequence_type"] == "nucleotide":
            lines.append(f"    Longitud:       {seq['length_bp']:,} bp")
            lines.append(f"    GC content:     {seq['gc_content']}%")
            lines.append(f"    AT content:     {seq['at_content']}%")
            comp = seq["composition_pct"]
            lines.append(f"    Composicion:    A={comp['A']}%  T={comp['T']}%  G={comp['G']}%  C={comp['C']}%")
            if seq.get("molecular_weight_da"):
                lines.append(f"    Peso molecular: {seq['molecular_weight_da']:,.1f} Da")
            if seq["orfs_found"] > 0:
                lines.append(f"    ORFs (>=100aa): {seq['orfs_found']}")
                for orf in seq["orfs"][:5]:
                    lines.append(f"      {orf['frame']} pos {orf['start']}-{orf['end']} ({orf['length_aa']} aa)")
        else:
            lines.append(f"    Longitud:          {seq['length_aa']} aa")
            if seq.get("molecular_weight_da"):
                lines.append(f"    Peso molecular:    {seq['molecular_weight_da']:,.1f} Da")
            if seq.get("isoelectric_point"):
                lines.append(f"    Punto isoelectrico:{seq['isoelectric_point']}")
            if seq.get("instability_index"):
                lines.append(f"    Indice estabilidad: {seq['instability_index']} ({seq.get('stability', '')})")
            if seq.get("gravy"):
                lines.append(f"    GRAVY:             {seq['gravy']} ({seq.get('hydrophobicity', '')})")
            if seq.get("aromaticity"):
                lines.append(f"    Aromaticidad:      {seq['aromaticity']}")
            if seq.get("secondary_structure_pct"):
                ss = seq["secondary_structure_pct"]
                lines.append(f"    Estr. secundaria:  helice={ss['helix']}%  lamina={ss['sheet']}%  giro={ss['turn']}%")
            if seq.get("charged_residues_pct"):
                lines.append(f"    Residuos cargados: {seq['charged_residues_pct']}%")

        lines.append("    " + "-" * 40)

    if len(report["sequences"]) > 20:
        lines.append(f"\n  ... y {len(report['sequences']) - 20} secuencias mas.")

    lines.append("")
    lines.append("=" * 60)
    return "\n".join(lines)


# ──────────────────────────────────────────────
# Formato HTML
# ──────────────────────────────────────────────

def generate_html_report(report, json_str):
    """Genera un reporte HTML completo con analisis, JSON y seccion 'Como funciona'."""

    s = report["summary"]
    seq_type = report["sequence_type"]
    is_nuc = seq_type == "nucleotide"

    # ── Summary table rows
    summary_rows = f"""
        <tr><td>Tipo de secuencia</td><td><span class="badge {'nuc' if is_nuc else 'prot'}">{seq_type}</span></td></tr>
        <tr><td>Total secuencias</td><td>{s['total_sequences']}</td></tr>
        <tr><td>Total residuos</td><td>{s['total_residues']:,} {'bp' if is_nuc else 'aa'}</td></tr>
        <tr><td>Longitud minima</td><td>{s['min_length']:,}</td></tr>
        <tr><td>Longitud maxima</td><td>{s['max_length']:,}</td></tr>
        <tr><td>Longitud promedio</td><td>{s['avg_length']:,}</td></tr>
        <tr><td>N50</td><td>{s['n50']:,}</td></tr>"""

    if is_nuc:
        summary_rows += f"""
        <tr><td>GC promedio</td><td>{s.get('avg_gc_content', 'N/A')}%</td></tr>
        <tr><td>ORFs encontrados (&ge;100 aa)</td><td>{s.get('total_orfs', 0)}</td></tr>"""

    # ── Per-sequence detail cards
    seq_cards = ""
    for i, seq in enumerate(report["sequences"][:20]):
        sid = html_mod.escape(seq.get("id", "unknown"))
        sdesc = html_mod.escape(seq.get("description", ""))

        if is_nuc:
            comp = seq["composition_pct"]
            gc = seq["gc_content"]
            at = seq["at_content"]

            # Composition bar widths
            bar_a = comp["A"]
            bar_t = comp["T"]
            bar_g = comp["G"]
            bar_c = comp["C"]

            orfs_html = ""
            if seq["orfs_found"] > 0:
                orf_rows = ""
                for orf in seq["orfs"][:8]:
                    orf_rows += f'<tr><td>{orf["frame"]}</td><td>{orf["start"]:,} - {orf["end"]:,}</td><td>{orf["length_bp"]:,} bp</td><td>{orf["length_aa"]:,} aa</td></tr>'
                orfs_html = f"""
                <h4>ORFs detectados ({seq['orfs_found']})</h4>
                <table class="inner-table">
                    <thead><tr><th>Frame</th><th>Posicion</th><th>Longitud</th><th>Proteina</th></tr></thead>
                    <tbody>{orf_rows}</tbody>
                </table>"""

            dinuc_bars = ""
            for dn, count in seq.get("top_dinucleotides", {}).items():
                max_dn = max(seq["top_dinucleotides"].values()) if seq["top_dinucleotides"] else 1
                w = count / max_dn * 100
                dinuc_bars += f'<div class="mini-bar-row"><span class="mini-label">{dn}</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:{w}%"></div></div><span class="mini-val">{count:,}</span></div>'

            seq_cards += f"""
            <div class="seq-card">
                <div class="seq-header">
                    <span class="seq-num">#{i+1}</span>
                    <div>
                        <h3>{sid}</h3>
                        <p class="seq-desc">{sdesc}</p>
                    </div>
                </div>
                <div class="seq-grid">
                    <div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Longitud</span><span class="metric-value">{seq['length_bp']:,} bp</span></div>
                            <div class="metric"><span class="metric-label">GC content</span><span class="metric-value">{gc}%</span></div>
                            <div class="metric"><span class="metric-label">AT content</span><span class="metric-value">{at}%</span></div>
                        </div>
                        <h4>Composicion de bases</h4>
                        <div class="comp-bar">
                            <div class="comp-seg seg-a" style="width:{bar_a}%" title="A: {bar_a}%">A</div>
                            <div class="comp-seg seg-t" style="width:{bar_t}%" title="T: {bar_t}%">T</div>
                            <div class="comp-seg seg-g" style="width:{bar_g}%" title="G: {bar_g}%">G</div>
                            <div class="comp-seg seg-c" style="width:{bar_c}%" title="C: {bar_c}%">C</div>
                        </div>
                        <div class="comp-legend">
                            <span><i class="dot dot-a"></i>A {bar_a}%</span>
                            <span><i class="dot dot-t"></i>T {bar_t}%</span>
                            <span><i class="dot dot-g"></i>G {bar_g}%</span>
                            <span><i class="dot dot-c"></i>C {bar_c}%</span>
                        </div>
                    </div>
                    <div>
                        <h4>Top dinucleotidos</h4>
                        {dinuc_bars}
                    </div>
                </div>
                {orfs_html}
            </div>"""

        else:  # protein
            mw = seq.get("molecular_weight_da", "N/A")
            mw_str = f"{mw:,.1f} Da" if isinstance(mw, (int, float)) else mw
            pi = seq.get("isoelectric_point", "N/A")
            stab = seq.get("instability_index", "N/A")
            stab_label = seq.get("stability", "")
            gravy = seq.get("gravy", "N/A")
            hydro = seq.get("hydrophobicity", "")
            arom = seq.get("aromaticity", "N/A")
            ss = seq.get("secondary_structure_pct", {})
            charged = seq.get("charged_residues_pct", "N/A")
            aromatic_pct = seq.get("aromatic_residues_pct", "N/A")

            ss_html = ""
            if ss:
                ss_html = f"""
                <h4>Estructura secundaria estimada</h4>
                <div class="ss-bars">
                    <div class="ss-row"><span class="ss-label">Helice</span><div class="mini-track"><div class="mini-fill helix-fill" style="width:{ss.get('helix',0)}%"></div></div><span class="mini-val">{ss.get('helix',0)}%</span></div>
                    <div class="ss-row"><span class="ss-label">Lamina</span><div class="mini-track"><div class="mini-fill sheet-fill" style="width:{ss.get('sheet',0)}%"></div></div><span class="mini-val">{ss.get('sheet',0)}%</span></div>
                    <div class="ss-row"><span class="ss-label">Giro</span><div class="mini-track"><div class="mini-fill turn-fill" style="width:{ss.get('turn',0)}%"></div></div><span class="mini-val">{ss.get('turn',0)}%</span></div>
                </div>"""

            seq_cards += f"""
            <div class="seq-card">
                <div class="seq-header">
                    <span class="seq-num">#{i+1}</span>
                    <div>
                        <h3>{sid}</h3>
                        <p class="seq-desc">{sdesc}</p>
                    </div>
                </div>
                <div class="seq-grid">
                    <div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Longitud</span><span class="metric-value">{seq['length_aa']:,} aa</span></div>
                            <div class="metric"><span class="metric-label">Peso molecular</span><span class="metric-value">{mw_str}</span></div>
                            <div class="metric"><span class="metric-label">Punto isoelectrico</span><span class="metric-value">{pi}</span></div>
                        </div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Estabilidad</span><span class="metric-value">{stab} <small>({stab_label})</small></span></div>
                            <div class="metric"><span class="metric-label">GRAVY</span><span class="metric-value">{gravy} <small>({hydro})</small></span></div>
                            <div class="metric"><span class="metric-label">Aromaticidad</span><span class="metric-value">{arom}</span></div>
                        </div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Residuos cargados</span><span class="metric-value">{charged}%</span></div>
                            <div class="metric"><span class="metric-label">Residuos aromaticos</span><span class="metric-value">{aromatic_pct}%</span></div>
                        </div>
                    </div>
                    <div>{ss_html}</div>
                </div>
            </div>"""

    overflow_note = ""
    if len(report["sequences"]) > 20:
        overflow_note = f'<p class="overflow-note">Mostrando 20 de {len(report["sequences"])} secuencias.</p>'

    # ── Escaped JSON for display
    json_escaped = html_mod.escape(json_str)

    # ── Full HTML
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analisis FASTA - {html_mod.escape(report['file'])}</title>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Merriweather:wght@400;700;900&family=Source+Sans+3:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
        :root {{
            --bg:#fff; --surface:#fafafa; --surface-2:#f3f3f0; --border:#ddd9d0; --border-lt:#e8e5de;
            --text:#1a1a1a; --text-2:#4a4a4a; --text-dim:#6b6b6b; --text-cap:#888;
            --uy:#FFC600; --uyd:#e0ad00; --uys:#fffbeb; --ub:#1a1a1a;
            --blue:#2c5282; --green:#276749; --amber:#b7791f; --red:#9b2c2c;
            --blue-bg:#ebf0f7; --green-bg:#e6f4ed; --amber-bg:#fef6e6; --red-bg:#fde8e8;
        }}
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:'Source Sans 3',system-ui,sans-serif; background:var(--bg); color:var(--text); line-height:1.7; }}
        .top-bar {{ background:var(--ub); padding:0.6rem 2rem; display:flex; align-items:center; justify-content:space-between; }}
        .top-bar-logo {{ font-family:'Merriweather',serif; font-weight:900; font-size:1.1rem; color:var(--uy); }}
        .top-bar-info {{ font-size:0.75rem; color:rgba(255,255,255,0.5); }}
        .yellow-strip {{ height:4px; background:var(--uy); }}
        .container {{ max-width:940px; margin:0 auto; padding:2rem 2rem 4rem; }}

        .hero {{ text-align:center; padding:2.5rem 0 2rem; border-bottom:2px solid var(--uy); margin-bottom:2.5rem; }}
        .hero h1 {{ font-family:'Merriweather',serif; font-size:1.6rem; font-weight:700; color:var(--ub); margin-bottom:0.5rem; }}
        .hero-meta {{ font-size:0.82rem; color:var(--text-cap); margin-top:0.5rem; }}
        .hero-file {{ font-family:'JetBrains Mono',monospace; font-size:0.9rem; color:var(--blue); background:var(--blue-bg); padding:0.2rem 0.7rem; border-radius:3px; }}

        .sn {{ display:inline-flex; align-items:center; gap:0.5rem; font-family:'Merriweather',serif; font-size:0.72rem; font-weight:700; color:var(--uyd); text-transform:uppercase; letter-spacing:0.1em; margin-bottom:0.6rem; }}
        .sn::before {{ content:''; width:24px; height:3px; background:var(--uy); border-radius:2px; }}
        h2 {{ font-family:'Merriweather',serif; font-size:1.35rem; font-weight:700; color:var(--ub); margin-bottom:0.8rem; }}
        h3 {{ font-family:'Merriweather',serif; font-size:1rem; font-weight:700; color:var(--ub); margin-bottom:0.3rem; }}
        h4 {{ font-family:'Merriweather',serif; font-size:0.82rem; font-weight:700; color:var(--text-dim); text-transform:uppercase; letter-spacing:0.06em; margin:1rem 0 0.5rem; }}
        p {{ color:var(--text-2); font-size:0.93rem; margin-bottom:0.7rem; }}
        section {{ margin-bottom:3rem; }}
        small {{ color:var(--text-dim); }}

        .badge {{ display:inline-block; padding:0.12rem 0.5rem; border-radius:3px; font-family:'JetBrains Mono',monospace; font-size:0.72rem; font-weight:600; }}
        .badge.nuc {{ background:var(--blue-bg); color:var(--blue); }}
        .badge.prot {{ background:var(--amber-bg); color:var(--amber); }}

        /* Summary table */
        .summary-table {{ width:100%; border-collapse:collapse; margin-top:1rem; }}
        .summary-table td {{ padding:0.55rem 0.8rem; font-size:0.88rem; border-bottom:1px solid var(--border-lt); }}
        .summary-table td:first-child {{ font-weight:600; color:var(--text); width:45%; }}
        .summary-table td:last-child {{ color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.82rem; }}
        .summary-table tr:hover td {{ background:var(--uys); }}

        /* Sequence cards */
        .seq-card {{ background:var(--surface); border:1px solid var(--border); border-radius:4px; padding:1.5rem 1.8rem; margin-bottom:1.2rem; }}
        .seq-header {{ display:flex; align-items:baseline; gap:0.8rem; margin-bottom:1rem; padding-bottom:0.8rem; border-bottom:1px solid var(--border-lt); }}
        .seq-num {{ font-family:'Merriweather',serif; font-size:1.5rem; font-weight:900; color:var(--surface-2); }}
        .seq-desc {{ font-size:0.8rem; color:var(--text-cap); margin:0; word-break:break-all; }}
        .seq-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; }}
        @media(max-width:640px) {{ .seq-grid {{ grid-template-columns:1fr; }} }}

        .metric-group {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:0.5rem; margin-bottom:0.8rem; }}
        .metric {{ background:var(--bg); border:1px solid var(--border-lt); border-radius:4px; padding:0.5rem 0.7rem; }}
        .metric-label {{ display:block; font-size:0.68rem; font-weight:600; text-transform:uppercase; letter-spacing:0.05em; color:var(--text-cap); }}
        .metric-value {{ font-family:'JetBrains Mono',monospace; font-size:0.88rem; font-weight:600; color:var(--text); }}

        /* Composition bar */
        .comp-bar {{ display:flex; height:28px; border-radius:4px; overflow:hidden; margin:0.5rem 0; }}
        .comp-seg {{ display:flex; align-items:center; justify-content:center; font-size:0.7rem; font-weight:700; color:#fff; min-width:20px; }}
        .seg-a {{ background:#4299e1; }}
        .seg-t {{ background:#f56565; }}
        .seg-g {{ background:#48bb78; }}
        .seg-c {{ background:#ed8936; }}
        .comp-legend {{ display:flex; gap:1rem; font-size:0.75rem; color:var(--text-dim); flex-wrap:wrap; }}
        .dot {{ display:inline-block; width:8px; height:8px; border-radius:2px; margin-right:3px; vertical-align:middle; }}
        .dot-a {{ background:#4299e1; }} .dot-t {{ background:#f56565; }} .dot-g {{ background:#48bb78; }} .dot-c {{ background:#ed8936; }}

        /* Mini bar charts */
        .mini-bar-row {{ display:flex; align-items:center; gap:0.5rem; margin-bottom:0.3rem; }}
        .mini-label {{ font-family:'JetBrains Mono',monospace; font-size:0.72rem; width:28px; text-align:right; color:var(--text-dim); }}
        .mini-track {{ flex:1; height:16px; background:var(--surface-2); border-radius:3px; overflow:hidden; }}
        .mini-fill {{ height:100%; border-radius:3px; }}
        .nuc-fill {{ background:var(--blue); opacity:0.6; }}
        .mini-val {{ font-family:'JetBrains Mono',monospace; font-size:0.68rem; color:var(--text-cap); width:50px; }}

        /* Inner table (ORFs) */
        .inner-table {{ width:100%; border-collapse:collapse; margin-top:0.5rem; font-size:0.82rem; }}
        .inner-table th {{ text-align:left; padding:0.4rem 0.6rem; font-size:0.68rem; font-weight:700; letter-spacing:0.06em; text-transform:uppercase; color:var(--text-dim); background:var(--surface-2); border-bottom:2px solid var(--uy); }}
        .inner-table td {{ padding:0.35rem 0.6rem; border-bottom:1px solid var(--border-lt); color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.78rem; }}

        /* SS bars */
        .ss-bars {{ margin:0.3rem 0; }}
        .ss-row {{ display:flex; align-items:center; gap:0.5rem; margin-bottom:0.3rem; }}
        .ss-label {{ font-size:0.75rem; width:55px; text-align:right; color:var(--text-dim); }}
        .helix-fill {{ background:var(--red); opacity:0.6; }}
        .sheet-fill {{ background:var(--blue); opacity:0.6; }}
        .turn-fill {{ background:var(--green); opacity:0.6; }}

        .overflow-note {{ font-size:0.82rem; color:var(--text-cap); font-style:italic; text-align:center; margin-top:1rem; }}

        /* JSON block */
        .json-toggle {{ background:var(--ub); color:var(--uy); border:none; padding:0.45rem 1rem; border-radius:4px; font-family:'Source Sans 3',sans-serif; font-size:0.82rem; font-weight:600; cursor:pointer; margin-bottom:0.8rem; }}
        .json-toggle:hover {{ background:#333; }}
        .json-block {{ background:var(--ub); border-radius:4px; padding:1.2rem 1.5rem; overflow-x:auto; display:none; }}
        .json-block pre {{ font-family:'JetBrains Mono',monospace; font-size:0.75rem; line-height:1.6; color:#e4e4ed; white-space:pre-wrap; word-break:break-word; }}

        /* Como funciona */
        .how-section {{ background:var(--uys); border:1px solid var(--uy); border-radius:4px; padding:2rem 2.2rem; position:relative; }}
        .how-section::before {{ content:''; position:absolute; top:0; left:0; width:100%; height:4px; background:var(--uy); border-radius:4px 4px 0 0; }}
        .how-section h2 {{ margin-bottom:1rem; }}
        .how-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; margin-top:1rem; }}
        @media(max-width:640px) {{ .how-grid {{ grid-template-columns:1fr; }} }}
        .how-card {{ background:var(--bg); border:1px solid var(--border); border-radius:4px; padding:1.2rem; }}
        .how-card h3 {{ font-size:0.88rem; margin-bottom:0.4rem; }}
        .how-card p {{ font-size:0.85rem; margin-bottom:0.3rem; }}
        .how-card code {{ font-family:'JetBrains Mono',monospace; font-size:0.78em; background:var(--surface-2); padding:0.1em 0.35em; border-radius:2px; }}
        .how-card ul {{ list-style:none; padding:0; margin:0.3rem 0 0; }}
        .how-card ul li {{ font-size:0.82rem; padding:0.2rem 0 0.2rem 0.9rem; position:relative; color:var(--text-2); }}
        .how-card ul li::before {{ content:''; position:absolute; left:0; top:0.6rem; width:5px; height:5px; border-radius:50%; background:var(--uyd); }}

        .pipe-flow {{ display:flex; align-items:center; justify-content:center; gap:0.4rem; flex-wrap:wrap; margin:1rem 0; }}
        .pipe-node {{ padding:0.5rem 0.8rem; border-radius:4px; font-size:0.78rem; font-weight:600; text-align:center; }}
        .pipe-arrow {{ color:var(--text-cap); font-size:1rem; }}
        .pipe-node.n-input {{ background:var(--blue-bg); color:var(--blue); border:1px solid rgba(44,82,130,0.2); }}
        .pipe-node.n-python {{ background:var(--green-bg); color:var(--green); border:1px solid rgba(39,103,73,0.2); }}
        .pipe-node.n-claude {{ background:var(--amber-bg); color:var(--amber); border:1px solid rgba(183,121,31,0.2); }}
        .pipe-node.n-output {{ background:var(--red-bg); color:var(--red); border:1px solid rgba(155,44,44,0.2); }}

        .divider {{ width:50px; height:3px; background:var(--uy); margin:2.5rem auto; }}
        .footer {{ margin-top:2.5rem; padding:1.2rem 0; border-top:3px solid var(--uy); text-align:center; }}
        .footer p {{ font-size:0.78rem; color:var(--text-cap); margin:0; }}
    </style>
</head>
<body>

<div class="top-bar">
    <div class="top-bar-logo">FASTA Analyzer</div>
    <div class="top-bar-info">analyze_fasta.py &middot; Biopython</div>
</div>
<div class="yellow-strip"></div>

<div class="container">

    <header class="hero">
        <h1>Reporte de analisis FASTA</h1>
        <div><span class="hero-file">{html_mod.escape(report['file'])}</span></div>
        <div class="hero-meta">Generado el {report.get('analysis_date', '')} &middot; analyze_fasta.py</div>
    </header>

    <!-- ===== RESUMEN ===== -->
    <section>
        <div class="sn">1. Resumen general</div>
        <h2>Estadisticas del archivo</h2>
        <table class="summary-table">
            <tbody>{summary_rows}</tbody>
        </table>
    </section>

    <!-- ===== DETALLE POR SECUENCIA ===== -->
    <section>
        <div class="sn">2. Analisis por secuencia</div>
        <h2>Detalle de cada secuencia</h2>
        {seq_cards}
        {overflow_note}
    </section>

    <div class="divider"></div>

    <!-- ===== JSON OUTPUT ===== -->
    <section>
        <div class="sn">3. JSON output</div>
        <h2>Datos crudos del analisis</h2>
        <p>Output completo en formato JSON, tal como lo consume Claude para la interpretacion biologica.</p>
        <button class="json-toggle" onclick="var b=document.getElementById('json-block');b.style.display=b.style.display==='block'?'none':'block';this.textContent=b.style.display==='block'?'Ocultar JSON':'Mostrar JSON'">Mostrar JSON</button>
        <div class="json-block" id="json-block">
            <pre>{json_escaped}</pre>
        </div>
    </section>

    <div class="divider"></div>

    <!-- ===== COMO FUNCIONA ===== -->
    <section>
        <div class="sn">4. Como funciona</div>
        <div class="how-section">
            <h2>Skill: analyze-fasta</h2>
            <p>
                Esta herramienta combina un <strong>script de Python</strong> que realiza el analisis bioinformatico
                con <strong>Claude (IA de Anthropic)</strong> que interpreta los resultados en contexto biologico.
                El script es el musculo (calculos); Claude es el cerebro (interpretacion).
            </p>

            <div class="pipe-flow">
                <div class="pipe-node n-input">Archivo FASTA</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-python">analyze_fasta.py<br><small>Biopython + NumPy</small></div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-output">JSON + HTML</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-claude">Claude (skill)<br><small>Interpretacion biologica</small></div>
            </div>

            <div class="how-grid">
                <div class="how-card">
                    <h3>Librerias utilizadas</h3>
                    <ul>
                        <li><strong>Biopython</strong> (<code>Bio.SeqIO</code>) &mdash; parsing de FASTA, acceso a secuencias</li>
                        <li><strong>Bio.SeqUtils</strong> &mdash; calculo de GC%, peso molecular</li>
                        <li><strong>Bio.SeqUtils.ProtParam</strong> &mdash; analisis fisicoquimico de proteinas (pI, GRAVY, estabilidad, estructura secundaria)</li>
                        <li><strong>NumPy</strong> &mdash; dependencia de Biopython para calculos numericos</li>
                        <li><strong>Python stdlib</strong> &mdash; Counter, json, argparse, pathlib</li>
                    </ul>
                </div>

                <div class="how-card">
                    <h3>Que analiza</h3>
                    <ul>
                        <li><strong>Deteccion automatica</strong> de tipo (nucleotido vs proteina) por composicion</li>
                        <li><strong>Nucleotidos:</strong> GC%, composicion de bases, dinucleotidos, ORFs en 3 frames, peso molecular, N50</li>
                        <li><strong>Proteinas:</strong> peso molecular, punto isoelectrico, indice de estabilidad, GRAVY, aromaticidad, estructura secundaria, residuos cargados</li>
                        <li><strong>Multi-secuencia:</strong> estadisticas agregadas, N50 para evaluar calidad de ensamblaje</li>
                    </ul>
                </div>

                <div class="how-card">
                    <h3>Deteccion de ORFs</h3>
                    <p>
                        Busca marcos de lectura abiertos en los 3 frames forward (+1, +2, +3).
                        Un ORF se define como una secuencia que empieza con <code>ATG</code> y
                        termina en un codon stop (<code>TAA</code>, <code>TAG</code>, <code>TGA</code>)
                        con un minimo de 100 aminoacidos (300 bp).
                    </p>
                    <p>
                        <strong>Limitacion:</strong> no busca en cadena reversa complementaria ni modela intrones.
                        Para prediccion genica completa se recomienda Augustus o GenScan.
                    </p>
                </div>

                <div class="how-card">
                    <h3>Rol de Claude (skill)</h3>
                    <p>
                        La skill <code>analyze-fasta.md</code> instruye a Claude para:
                    </p>
                    <ul>
                        <li>Ejecutar el script con <code>--json</code></li>
                        <li>Interpretar GC%, ORFs, pI, GRAVY en contexto biologico</li>
                        <li>Inferir organismo, funcion, o tipo de secuencia probable</li>
                        <li>Comparar con valores de referencia</li>
                        <li>Sugerir analisis adicionales (BLAST, filogenia, etc.)</li>
                    </ul>
                </div>
            </div>

            <h4 style="margin-top:1.5rem">Modos de uso</h4>
            <table class="summary-table" style="margin-top:0.5rem">
                <tr><td><code>python3 analyze_fasta.py archivo.fasta</code></td><td>Reporte en texto por terminal</td></tr>
                <tr><td><code>python3 analyze_fasta.py archivo.fasta --json</code></td><td>Output JSON (para Claude o pipelines)</td></tr>
                <tr><td><code>python3 analyze_fasta.py archivo.fasta --html reporte.html</code></td><td>Genera este reporte HTML</td></tr>
            </table>
        </div>
    </section>

    <footer class="footer">
        <p>Generado por analyze_fasta.py &middot; Biopython {get_biopython_version()} &middot; Python {sys.version.split()[0]}</p>
    </footer>
</div>

</body>
</html>"""


def get_biopython_version():
    try:
        import Bio
        return Bio.__version__
    except Exception:
        return "?"


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

CLAWBIO_DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def format_markdown_report(report):
    """Generate a clean Markdown report with ClawBio disclaimer."""
    lines = []
    lines.append(f"# analyze-fasta Report")
    lines.append("")
    lines.append(f"**Input file:** `{report['file']}`  ")
    lines.append(f"**Analysis date:** {report['analysis_date']}  ")
    lines.append(f"**Sequence type:** `{report['sequence_type']}`  ")
    lines.append(f"**Total sequences:** {report['total_sequences']}  ")
    lines.append("")

    summary = report.get("summary", {})
    if summary:
        lines.append("## Summary")
        lines.append("")
        lines.append("| Metric | Value |")
        lines.append("|---|---|")
        for k, v in summary.items():
            lines.append(f"| {k} | {v} |")
        lines.append("")

    lines.append("## Per-sequence metrics")
    lines.append("")
    for i, seq in enumerate(report.get("sequences", [])[:20], 1):
        lines.append(f"### {i}. {seq.get('id', 'unknown')}")
        lines.append("")
        lines.append(f"- **Description:** {seq.get('description', '')}")
        if report["sequence_type"] == "nucleotide":
            lines.append(f"- **Length:** {seq.get('length_bp', 0)} bp")
            lines.append(f"- **GC content:** {seq.get('gc_content', 0)}%")
            lines.append(f"- **AT content:** {seq.get('at_content', 0)}%")
            orfs = seq.get("orfs", [])
            lines.append(f"- **ORFs (>=100 aa):** {len(orfs)}")
            if seq.get("molecular_weight_da"):
                lines.append(f"- **Molecular weight:** {seq['molecular_weight_da']} Da")
        else:
            lines.append(f"- **Length:** {seq.get('length_aa', 0)} aa")
            lines.append(f"- **Molecular weight:** {seq.get('molecular_weight_da', 0)} Da")
            lines.append(f"- **Isoelectric point (pI):** {seq.get('isoelectric_point', 0)}")
            lines.append(f"- **GRAVY:** {seq.get('gravy', 0)}")
            lines.append(f"- **Instability index:** {seq.get('instability_index', 0)}")
        lines.append("")

    if len(report.get("sequences", [])) > 20:
        lines.append(f"_Showing first 20 of {len(report['sequences'])} sequences. Full data in `result.json`._")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(f"_{CLAWBIO_DISCLAIMER}_")
    return "\n".join(lines) + "\n"


def write_clawbio_output(report, output_dir, input_file):
    """Write ClawBio-convention output: report.md + result.json + reproducibility/."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    json_str = json.dumps(report, indent=2, ensure_ascii=False)
    (out / "result.json").write_text(json_str, encoding="utf-8")
    (out / "report.md").write_text(format_markdown_report(report), encoding="utf-8")

    # Reproducibility bundle
    repro = out / "reproducibility"
    repro.mkdir(exist_ok=True)
    (repro / "commands.sh").write_text(
        f"#!/bin/bash\n"
        f"# Reproduce this analyze-fasta run:\n"
        f"python3 {Path(__file__).name} --input {input_file} --output {out}\n",
        encoding="utf-8",
    )
    run_meta = {
        "skill": "analyze-fasta",
        "version": "0.1.0",
        "input_file": str(input_file),
        "input_size_bytes": Path(input_file).stat().st_size if Path(input_file).exists() else 0,
        "ran_at": datetime.now().isoformat(),
        "biopython_version": get_biopython_version(),
        "python_version": sys.version.split()[0],
    }
    (repro / "run.json").write_text(json.dumps(run_meta, indent=2), encoding="utf-8")

    # Optional HTML alongside the markdown for visual inspection
    (out / "report.html").write_text(generate_html_report(report, json_str), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Analisis general de archivos FASTA (ClawBio skill)")
    parser.add_argument("fasta", nargs="?", help="(legacy) Archivo FASTA posicional")
    parser.add_argument("--input", "-i", help="Archivo FASTA a analizar (ClawBio convention)")
    parser.add_argument("--output", "-o", help="Directorio de salida (genera report.md + result.json + reproducibility/)")
    parser.add_argument("--demo", action="store_true", help="Modo demo con datos sinteticos bundleados")
    parser.add_argument("--json", action="store_true", help="(legacy) imprime JSON a stdout")
    parser.add_argument("--html", nargs="?", const="auto", default=None, metavar="OUTPUT",
                        help="(legacy) Genera reporte HTML standalone")
    args = parser.parse_args()

    # Resolve input
    if args.demo:
        input_file = Path(__file__).resolve().parent / "example_data" / "demo_nucleotide.fasta"
    elif args.input:
        input_file = Path(args.input)
    elif args.fasta:
        input_file = Path(args.fasta)
    else:
        parser.error("Specify --input FILE, --demo, or pass a positional FASTA path")

    report = analyze_fasta(str(input_file))

    if "error" in report:
        print(f"ERROR: {report['error']}", file=sys.stderr)
        sys.exit(1)

    json_str = json.dumps(report, indent=2, ensure_ascii=False)

    # ClawBio convention: --output DIR writes report.md + result.json
    if args.output:
        write_clawbio_output(report, args.output, input_file)
        print(f"Wrote {args.output}/report.md")
        print(f"Wrote {args.output}/result.json")
        return

    # Legacy paths
    if args.html is not None:
        if args.html == "auto":
            fasta_stem = Path(input_file).stem
            project_root = Path(input_file).resolve().parent
            output_dir = project_root / "output"
            output_dir.mkdir(exist_ok=True)
            output_path = output_dir / f"output_{fasta_stem}.html"
        else:
            output_path = Path(args.html)
        html_content = generate_html_report(report, json_str)
        output_path.write_text(html_content, encoding="utf-8")
        print(f"Reporte HTML generado: {output_path}")
    elif args.json:
        print(json_str)
    else:
        print(format_text_report(report))


if __name__ == "__main__":
    main()
