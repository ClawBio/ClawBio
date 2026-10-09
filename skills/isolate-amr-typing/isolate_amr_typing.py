#!/usr/bin/env python3
"""Isolate AMR Typing - genotypic characterisation of one bacterial isolate assembly.

Combines three established tools into a single report:
  - mlst (Seemann)            sequence type from the PubMLST schemes
  - AMRFinderPlus (NCBI)      AMR genes and, when the organism is known, point mutations
  - PlasmidFinder (CGE)       plasmid replicons

The skill either runs the tools on an assembly (--input) or reads results you already
have (--mlst / --amrfinder / --plasmidfinder). It does not re-filter the tools' calls by
identity or coverage; it classifies, cross-references and reports them.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SKILL_DIR = Path(__file__).resolve().parent
SKILL_NAME = "isolate-amr-typing"
VERSION = "0.1.0"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

DEMO_FASTA = SKILL_DIR / "examples" / "demo_isolate.fasta"
DEMO_AMR = SKILL_DIR / "examples" / "demo_amrfinder.tsv"
DEMO_MLST = SKILL_DIR / "examples" / "demo_mlst.tsv"
DEMO_PLASMID = SKILL_DIR / "examples" / "demo_plasmidfinder.tsv"

# PlasmidFinder thresholds passed explicitly on live runs (Carattoli et al. 2014 defaults).
PLASMIDFINDER_MIN_COVERAGE = "0.60"
PLASMIDFINDER_MIN_IDENTITY = "0.95"

# mlst scheme -> AMRFinderPlus --organism. Only unambiguous, one-to-one cases belong
# here: a wrong organism applies the wrong point-mutation panel without any error.
SCHEME_TO_ORGANISM = {
    "ecoli": "Escherichia",
    "ecoli_achtman_4": "Escherichia",
    "salmonella": "Salmonella",
    "senterica": "Salmonella",
    "senterica_achtman_2": "Salmonella",
    "klebsiella": "Klebsiella_pneumoniae",
    "koxytoca": "Klebsiella_oxytoca",
    "saureus": "Staphylococcus_aureus",
    "spseudintermedius": "Staphylococcus_pseudintermedius",
    "paeruginosa": "Pseudomonas_aeruginosa",
    "abaumannii": "Acinetobacter_baumannii",
    "abaumannii_2": "Acinetobacter_baumannii",
    "efaecalis": "Enterococcus_faecalis",
    "efaecium": "Enterococcus_faecium",
    "campylobacter": "Campylobacter",
    "spneumoniae": "Streptococcus_pneumoniae",
    "spyogenes": "Streptococcus_pyogenes",
    "sagalactiae": "Streptococcus_agalactiae",
    "cdifficile": "Clostridioides_difficile",
    "vcholerae": "Vibrio_cholerae",
    "ecloacae": "Enterobacter_cloacae",
    "cfreundii": "Citrobacter_freundii",
    "bpseudomallei": "Burkholderia_pseudomallei",
}
# Schemes spanning several AMRFinderPlus organisms: the user must choose.
AMBIGUOUS_SCHEMES = {
    "neisseria": "Neisseria_gonorrhoeae or Neisseria_meningitidis",
}

# AMRFinderPlus "Method" prefixes, longest first so PARTIAL_CONTIG_END wins over PARTIAL.
_METHOD_PREFIXES = (
    ("PARTIAL_CONTIG_END", "partial_contig_end"),
    ("INTERNAL_STOP", "internal_stop"),
    ("PARTIAL", "partial"),
    ("ALLELE", "allele"),
    ("EXACT", "exact"),
    ("BLAST", "blast"),
    ("POINT", "point"),
    ("HMM", "hmm"),
)
# confident: counted as a determinant. review: reported, needs a human look.
# disrupted: gene has a premature stop; listed but never counted.
CALL_TIER = {
    "allele": "confident",
    "exact": "confident",
    "blast": "confident",
    "point": "confident",
    "partial": "review",
    "partial_contig_end": "review",
    "hmm": "review",
    "other": "review",
    "internal_stop": "disrupted",
}

# Canonical field -> header names used by AMRFinderPlus v4 and v3.
_AMR_COLUMNS = {
    "contig": ("Contig id",),
    "start": ("Start",),
    "stop": ("Stop",),
    "strand": ("Strand",),
    "symbol": ("Element symbol", "Gene symbol"),
    "name": ("Element name", "Sequence name"),
    "scope": ("Scope",),
    "element_type": ("Type", "Element type"),
    "subtype": ("Subtype", "Element subtype"),
    "drug_class": ("Class",),
    "subclass": ("Subclass",),
    "method": ("Method",),
    "coverage": ("% Coverage of reference", "% Coverage of reference sequence"),
    "identity": ("% Identity to reference", "% Identity to reference sequence"),
}
_AMR_REQUIRED = ("symbol", "element_type", "drug_class", "method")

_INSTALL_HINTS = {
    "amrfinder": "conda install -c bioconda ncbi-amrfinderplus && amrfinder -u",
    "mlst": "conda install -c bioconda mlst",
    "plasmidfinder.py": "conda install -c bioconda plasmidfinder (v2 gives plasmidfinder.py, v3 gives "
                        "`python -m plasmidfinder`), or abricate with its plasmidfinder database",
}
_NUCLEOTIDES = set("ACGTUNRYKMSWBDHV")


class MissingToolError(RuntimeError):
    """One or more external tools needed for a live run are not on PATH."""


# ── Small helpers ──────────────────────────────────────────────────────────────


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def _to_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _rows(path: Path) -> list[list[str]]:
    with open(path, encoding="utf-8") as fh:
        return [line.rstrip("\r\n").split("\t") for line in fh if line.strip()]


def sample_name_from(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".fasta", ".fna", ".fa", ".fas"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name or "isolate"


# ── Assembly ───────────────────────────────────────────────────────────────────


def read_assembly(path: Path) -> dict:
    """Validate a nucleotide FASTA and return contig lengths and summary statistics."""
    contigs: dict[str, int] = {}
    counts = dict.fromkeys("ACGT", 0)
    non_nucleotide = 0
    current = None
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0] if line[1:].split() else ""
                if not current:
                    raise ValueError(f"{path.name}: FASTA record with an empty name")
                if current in contigs:
                    raise ValueError(f"{path.name}: duplicate contig name '{current}'")
                contigs[current] = 0
                continue
            if current is None:
                raise ValueError(f"{path.name}: not a FASTA file (first record does not start with '>')")
            seq = line.upper()
            contigs[current] += len(seq)
            for base in "ACGT":
                counts[base] += seq.count(base)
            non_nucleotide += sum(1 for ch in seq if ch not in _NUCLEOTIDES)
    if not contigs or sum(contigs.values()) == 0:
        raise ValueError(f"{path.name}: no sequences found")
    if non_nucleotide:
        raise ValueError(
            f"{path.name}: {non_nucleotide} non-nucleotide characters; this skill needs a "
            "nucleotide assembly, not a protein FASTA"
        )
    total = sum(contigs.values())
    running, n50 = 0, 0
    for length in sorted(contigs.values(), reverse=True):
        running += length
        if running * 2 >= total:
            n50 = length
            break
    acgt = sum(counts.values())
    return {
        "file": path.name,
        "n_contigs": len(contigs),
        "total_length": total,
        "n50": n50,
        "longest_contig": max(contigs.values()),
        "gc_percent": round(100 * (counts["G"] + counts["C"]) / acgt, 2) if acgt else None,
        "contigs": contigs,
    }


# ── AMRFinderPlus ──────────────────────────────────────────────────────────────


def classify_method(method: str) -> str:
    method = (method or "").upper()
    for prefix, call in _METHOD_PREFIXES:
        if method.startswith(prefix):
            return call
    return "other"


def parse_amrfinder(path: Path) -> list[dict]:
    """Parse an AMRFinderPlus TSV (v3 or v4 column names) for a single isolate."""
    rows = _rows(path)
    if not rows:
        raise ValueError(
            f"{path.name}: file is empty. AMRFinderPlus writes a header even with no hits, "
            "so an empty file is a failed run, not a negative result"
        )
    header = [h.strip() for h in rows[0]]
    index = {}
    for field, names in _AMR_COLUMNS.items():
        for name in names:
            if name in header:
                index[field] = header.index(name)
                break
    missing = [f for f in _AMR_REQUIRED if f not in index]
    if missing:
        raise ValueError(
            f"{path.name}: not an AMRFinderPlus table (no column for: {', '.join(missing)})"
        )
    name_col = header.index("Name") if "Name" in header else None
    if name_col is not None and len({r[name_col] for r in rows[1:] if len(r) > name_col}) > 1:
        raise ValueError(f"{path.name}: results for several samples; this skill reports one isolate per run")

    def cell(row: list[str], field: str) -> str | None:
        i = index.get(field)
        if i is None or i >= len(row):
            return None
        value = row[i].strip()
        return None if value in ("", "NA") else value

    hits = []
    for row in rows[1:]:
        method = cell(row, "method") or ""
        hits.append({
            "symbol": cell(row, "symbol"),
            "name": cell(row, "name"),
            "element_type": cell(row, "element_type"),
            "subtype": cell(row, "subtype"),
            "drug_class": cell(row, "drug_class"),
            "subclass": cell(row, "subclass"),
            "scope": cell(row, "scope"),
            "method": method,
            "call": classify_method(method),
            "contig": cell(row, "contig"),
            "start": _to_int(cell(row, "start")),
            "stop": _to_int(cell(row, "stop")),
            "strand": cell(row, "strand"),
            "coverage": _to_float(cell(row, "coverage")),
            "identity": _to_float(cell(row, "identity")),
        })
    return hits


# ── MLST ───────────────────────────────────────────────────────────────────────


def parse_mlst_allele(token: str) -> dict:
    """Decode one mlst allele field: adk(53), adk(~53), adk(53?), adk(-), adk(53,54)."""
    match = re.match(r"^(.+?)\((.*)\)$", token.strip())
    if not match:
        raise ValueError(f"unrecognised mlst allele field '{token}'")
    locus, value = match.group(1), match.group(2)
    if value in ("-", ""):
        return {"locus": locus, "allele": None, "flag": "missing"}
    if "," in value:
        return {"locus": locus, "allele": value, "flag": "multiple"}
    if value.startswith("~"):
        return {"locus": locus, "allele": value[1:], "flag": "novel"}
    if value.endswith("?"):
        return {"locus": locus, "allele": value[:-1], "flag": "partial"}
    return {"locus": locus, "allele": value, "flag": "exact"}


def parse_mlst(path: Path) -> dict:
    """Parse default `mlst` output (FILE, SCHEME, ST, alleles...) for a single isolate."""
    rows = [r for r in _rows(path) if not (len(r) > 1 and r[1].strip().upper() == "SCHEME")]
    if not rows:
        raise ValueError(f"{path.name}: file is empty; mlst always prints one line per assembly")
    allele = re.compile(r"^.+\(.*\)$")
    for r in rows:
        if len(r) < 3 or not all(allele.match(tok.strip()) for tok in r[3:] if tok.strip()):
            raise ValueError(f"{path.name}: not mlst output (expected FILE, SCHEME, ST, then locus(allele) fields)")
    if len(rows) > 1:
        raise ValueError(f"{path.name}: results for several assemblies; this skill reports one isolate per run")
    row = [c.strip() for c in rows[0]]
    scheme = None if row[1] in ("-", "") else row[1]
    st = row[2] if row[2].isdigit() else None
    alleles = [parse_mlst_allele(tok) for tok in row[3:] if tok]

    if scheme is None:
        call, note = "no_scheme", "No MLST scheme matched this assembly."
    elif st is not None:
        call, note = "assigned", ""
    else:
        call = "unassigned"
        imperfect = [f"{a['locus']} ({a['flag']})" for a in alleles if a["flag"] != "exact"]
        if imperfect:
            note = "No ST: " + ", ".join(imperfect) + "."
        else:
            note = ("No ST: every allele is an exact match but the combination is not in the "
                    "scheme (possible novel ST).")
    return {"scheme": scheme, "st": st, "call": call, "note": note, "alleles": alleles}


def parse_mlst_ties(stderr: str) -> list[tuple[str, str]]:
    """Schemes mlst scored equally, from its stderr warning. Empty when there was no tie.

    mlst prints e.g. `WARNING: salmonella(3529)==ecoli_achtman_4(131) score=100 <file>` and
    then reports one of the two, not always the same one across runs.
    """
    match = re.search(r"WARNING:\s*(\S+?)\(([^)]*)\)==(\S+?)\(([^)]*)\)\s+score=", stderr or "")
    if not match:
        return []
    return [(match.group(1), match.group(2)), (match.group(3), match.group(4))]


def resolve_tie(ties: list[tuple[str, str]], organism: str | None) -> str | None:
    """The tied scheme belonging to the organism the user stated, if exactly one does."""
    if not organism:
        return None
    matching = [scheme for scheme, _ in ties if SCHEME_TO_ORGANISM.get(scheme) == organism]
    return matching[0] if len(matching) == 1 else None


# ── Organism ───────────────────────────────────────────────────────────────────


def infer_organism(scheme: str | None) -> str | None:
    """AMRFinderPlus organism for an mlst scheme, or None when not one-to-one."""
    return SCHEME_TO_ORGANISM.get(scheme or "")


def parse_organism_listing(text: str) -> set[str]:
    """Organism names from `amrfinder --list_organisms` output."""
    match = re.search(r"Available --organism options:\s*(.+)", text)
    if not match:
        return set()
    return {name.strip() for name in match.group(1).split(",") if name.strip()}


# ── PlasmidFinder ──────────────────────────────────────────────────────────────


def parse_plasmidfinder(path: Path) -> list[dict]:
    """Parse PlasmidFinder 2 `results_tab.tsv`, PlasmidFinder 3 `-j` JSON, or abricate output."""
    text = path.read_text(encoding="utf-8")
    if text.lstrip().startswith("{"):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path.name}: not valid PlasmidFinder JSON ({exc.msg})") from exc
        if not isinstance(data.get("seq_regions"), dict):
            raise ValueError(f"{path.name}: not PlasmidFinder 3 JSON (no 'seq_regions' block)")
        replicons = []
        for region in data["seq_regions"].values():
            contig = str(region.get("query_id") or "").split()
            replicons.append({
                "replicon": region.get("name"),
                "contig": contig[0] if contig else None,
                "start": _to_int(region.get("query_start_pos")),
                "stop": _to_int(region.get("query_end_pos")),
                "identity": _to_float(region.get("identity")),
                "coverage": _to_float(region.get("coverage")),
                "accession": region.get("ref_acc") or None,
            })
        return replicons
    rows = _rows(path)
    if not rows:
        raise ValueError(f"{path.name}: file is empty; expected at least a header line")
    header = [h.strip().lstrip("#") for h in rows[0]]
    col = {name: i for i, name in enumerate(header)}

    def get(row: list[str], name: str) -> str:
        i = col[name]
        return row[i].strip() if i < len(row) else ""

    replicons = []
    if {"Plasmid", "Query / Template length", "Contig"} <= col.keys():
        for row in rows[1:]:
            lengths = [_to_float(x) for x in get(row, "Query / Template length").split("/")]
            coverage = None
            if len(lengths) == 2 and lengths[0] is not None and lengths[1]:
                coverage = round(100 * lengths[0] / lengths[1], 2)
            position = get(row, "Position in contig").split("..") if "Position in contig" in col else []
            contig = get(row, "Contig").split()
            replicons.append({
                "replicon": get(row, "Plasmid"),
                "contig": contig[0] if contig else None,
                "start": _to_int(position[0]) if len(position) == 2 else None,
                "stop": _to_int(position[1]) if len(position) == 2 else None,
                "identity": _to_float(get(row, "Identity")),
                "coverage": coverage,
                "accession": (get(row, "Accession number") if "Accession number" in col else "") or None,
            })
    elif {"SEQUENCE", "GENE", "%IDENTITY", "%COVERAGE"} <= col.keys():
        for row in rows[1:]:
            replicons.append({
                "replicon": get(row, "GENE"),
                "contig": get(row, "SEQUENCE") or None,
                "start": _to_int(get(row, "START")),
                "stop": _to_int(get(row, "END")),
                "identity": _to_float(get(row, "%IDENTITY")),
                "coverage": _to_float(get(row, "%COVERAGE")),
                "accession": (get(row, "ACCESSION") if "ACCESSION" in col else "") or None,
            })
    else:
        raise ValueError(
            f"{path.name}: not PlasmidFinder results_tab.tsv, PlasmidFinder 3 JSON or abricate "
            "output (unrecognised header)"
        )
    return replicons


# ── Integration ────────────────────────────────────────────────────────────────


def interval_distance(a_start: int | None, a_stop: int | None,
                      b_start: int | None, b_stop: int | None) -> int | None:
    """Gap in bp between two features on one contig; 0 if they overlap. Contigs are linear."""
    if None in (a_start, a_stop, b_start, b_stop):
        return None
    a_lo, a_hi = sorted((a_start, a_stop))
    b_lo, b_hi = sorted((b_start, b_stop))
    return max(0, b_lo - a_hi, a_lo - b_hi)


def _unique(items: list[str]) -> list[str]:
    return list(dict.fromkeys(items))


def summarise(
    *,
    sample: str,
    assembly: dict | None,
    amr_hits: list[dict] | None,
    mlst: dict | None,
    replicons: list[dict] | None,
    organism: str | None,
    mode: str,
    applied_organism: str | None = None,
    amr_run_here: bool = False,
    mlst_ties: list[tuple[str, str]] | None = None,
    mlst_run_here: bool = False,
    tool_versions: dict | None = None,
) -> dict:
    """Merge the three tool results into one structured finding.

    `organism` is what the user stated; `applied_organism` is what a live run actually
    passed to AMRFinderPlus, and `amr_run_here` says this run executed AMRFinderPlus
    itself (so the organism used is known). None for a component means it was not run or provided,
    which is reported as such and never as an empty (negative) result.
    """
    warnings: list[str] = []
    scheme = mlst["scheme"] if mlst else None
    if mlst and mlst_ties:
        tied = " and ".join(f"{s} (ST{st})" for s, st in mlst_ties)
        mlst = {**mlst, "call": "ambiguous",
                "tied_schemes": [{"scheme": s, "st": st} for s, st in mlst_ties],
                "rerun_options": [
                    {"scheme": s, "st": st, "organism": SCHEME_TO_ORGANISM.get(s),
                     "flags": f"--mlst-scheme {s}"
                              + (f" --organism {SCHEME_TO_ORGANISM[s]}" if s in SCHEME_TO_ORGANISM else "")}
                    for s, st in mlst_ties
                ],
                "note": (f"mlst scored {tied} equally and reported one of them arbitrarily. "
                         "Rerun with the flags for the scheme that fits this isolate.")}

    by_contig: dict[str, list[dict]] = {}
    for rep in replicons or []:
        if rep["contig"]:
            by_contig.setdefault(rep["contig"], []).append(rep)

    genes, points, others = [], [], []
    for hit in amr_hits or []:
        hit = dict(hit)
        hit["tier"] = CALL_TIER[hit["call"]]
        neighbours = by_contig.get(hit["contig"], [])
        hit["replicons_on_contig"] = _unique([r["replicon"] for r in neighbours])
        gaps = [(interval_distance(hit["start"], hit["stop"], r["start"], r["stop"]), r["replicon"])
                for r in neighbours]
        gaps = sorted((g, name) for g, name in gaps if g is not None)
        hit["nearest_replicon"] = {"replicon": gaps[0][1], "distance_bp": gaps[0][0]} if gaps else None
        hit["contig_length"] = assembly["contigs"].get(hit["contig"]) if assembly else None
        if hit["element_type"] != "AMR":
            others.append(hit)
        elif hit["call"] == "point" or (hit["subtype"] or "").upper().startswith("POINT"):
            points.append(hit)
        else:
            genes.append(hit)

    classes: dict[str, dict[str, list[str]]] = {}
    for hit in genes + points:
        for drug_class in (hit["drug_class"] or "UNCLASSIFIED").split("/"):
            bucket = classes.setdefault(drug_class.strip(), {"confident": [], "review": [], "disrupted": []})
            bucket[hit["tier"]].append(hit["symbol"])
    class_summary = [
        {"drug_class": name, **{tier: _unique(symbols) for tier, symbols in buckets.items()}}
        for name, buckets in sorted(classes.items())
    ]

    # Organism and point-mutation screening status.
    if organism:
        resolved, source = organism, "user"
    elif applied_organism:
        resolved, source = applied_organism, "mlst_scheme"
    else:
        resolved, source = None, None

    if amr_hits is None:
        screening = "not_run"
    elif points:
        screening = "detected"
    elif organism or applied_organism:
        screening = "screened"
    elif amr_run_here:
        screening = "not_screened"
        reason = ("mlst could not choose between schemes" if mlst_ties
                  else "the MLST result was supplied, and a supplied result cannot show scheme "
                       "ties, so no organism was inferred from it" if mlst and not mlst_run_here
                  else f"MLST scheme '{scheme}' covers {AMBIGUOUS_SCHEMES[scheme]}" if scheme in AMBIGUOUS_SCHEMES
                  else f"MLST scheme '{scheme}' has no unambiguous AMRFinderPlus organism" if scheme
                  else "no MLST scheme was assigned")
        warnings.append(
            f"Point mutations were NOT screened: {reason}. Rerun with --organism to screen them."
        )
    else:
        screening = "unknown"
        warnings.append(
            "Point mutation screening status is unknown: AMRFinderPlus only reports point "
            "mutations when it was run with --organism, and the supplied table has none. "
            "Absence here is not evidence of absence; pass --organism if it was used."
        )

    if assembly is not None:
        known = assembly["contigs"]
        strangers = _unique([
            x["contig"] for x in (amr_hits or []) + (replicons or [])
            if x["contig"] and x["contig"] not in known
        ])
        if strangers:
            shown = ", ".join(strangers[:5]) + (" ..." if len(strangers) > 5 else "")
            warnings.append(
                f"Tool results name contigs that are not in {assembly['file']} ({shown}). "
                "The result files may belong to a different assembly."
            )
    review = _unique([h["symbol"] for h in genes if h["tier"] == "review"])
    if review:
        warnings.append(
            "Partial, contig-end or HMM-only hits need manual review before being treated as "
            "functional genes: " + ", ".join(review) + "."
        )
    if mlst and mlst["call"] != "assigned":
        warnings.append("MLST: " + mlst["note"])
    if mlst and not mlst_run_here and mode != "demo" and mlst["scheme"]:
        warnings.append(
            "The MLST result was supplied, not run here. mlst reports scheme ties only on "
            "stderr, so check its log: with a tie the scheme and ST shown may be the wrong one."
        )

    return {
        "skill": SKILL_NAME,
        "version": VERSION,
        "sample": sample,
        "mode": mode,
        "assembly": assembly,
        "mlst": {"status": "ok", **mlst} if mlst else {"status": "not_run"},
        "amr": {
            "status": "ok" if amr_hits is not None else "not_run",
            "organism": resolved,
            "organism_source": source,
            "organism_applied": bool(organism or applied_organism) if amr_hits is not None else False,
            "point_mutation_screening": screening,
            "genes": genes,
            "point_mutations": points,
            "other_elements": others,
            "class_summary": class_summary,
        },
        "plasmids": (
            {"status": "ok", "replicons": replicons} if replicons is not None else {"status": "not_run"}
        ),
        "warnings": warnings,
        "tool_versions": tool_versions or {},
        "disclaimer": DISCLAIMER,
    }


# ── Running the tools ──────────────────────────────────────────────────────────


_TOOL_ENV = {
    "mlst": "CLAWBIO_MLST_CMD",
    "amrfinder": "CLAWBIO_AMRFINDER_CMD",
    "plasmidfinder": "CLAWBIO_PLASMIDFINDER_CMD",
}


def tool_cmd(name: str) -> list[str]:
    """Command prefix for a tool, honouring its CLAWBIO_*_CMD override.

    The overrides exist because the tools do not all install into one conda
    environment; e.g. CLAWBIO_MLST_CMD="micromamba run -n mlst mlst".
    """
    override = os.environ.get(_TOOL_ENV.get(name, ""), "").strip()
    return shlex.split(override) if override else [name]


def override_preflight(env: dict) -> list[str]:
    """Shell lines that carry the tool-command overrides into commands.sh."""
    return [f'export {var}="${{{var}:-{env[var]}}}"' for var in _TOOL_ENV.values() if env.get(var)]


def build_mlst_cmd(assembly: Path, scheme: str | None = None) -> list[str]:
    return [*tool_cmd("mlst"), *(["--scheme", scheme] if scheme else []), str(assembly)]


def build_amrfinder_cmd(assembly: Path, out_tsv: Path, *, organism: str | None, threads: int) -> list[str]:
    cmd = [*tool_cmd("amrfinder"), "--nucleotide", str(assembly), "--output", str(out_tsv),
           "--plus", "--threads", str(threads)]
    if organism:
        cmd += ["--organism", organism]
    return cmd


def _module_available(python: str, name: str) -> bool:
    try:
        return subprocess.run([python, "-c", f"import {name}"], capture_output=True, timeout=60).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def find_plasmid_backend() -> tuple[str, list[str]] | None:
    """Pick how to run PlasmidFinder: (backend, command prefix), or None if unavailable.

    PlasmidFinder 2 installs `plasmidfinder.py` and writes results_tab.tsv. PlasmidFinder 3
    is a Python module only (`python -m plasmidfinder`) and writes JSON. abricate with its
    bundled plasmidfinder database is the last resort.
    """
    override = tool_cmd("plasmidfinder")
    if override != ["plasmidfinder"]:
        names = [Path(tok).name for tok in override]
        if "abricate" in names:
            return "abricate", override
        return ("plasmidfinder2" if names[-1] == "plasmidfinder.py" else "plasmidfinder3"), override
    if shutil.which("plasmidfinder.py"):
        return "plasmidfinder2", ["plasmidfinder.py"]
    for name in ("python3", "python"):
        python = shutil.which(name)
        if python and _module_available(python, "plasmidfinder"):
            return "plasmidfinder3", [python, "-m", "plasmidfinder"]
    if shutil.which("abricate"):
        return "abricate", ["abricate"]
    return None


def build_plasmid_cmd(backend: str, prefix: list[str], assembly: Path, out_dir: Path, *,
                      db: Path | None) -> list[str]:
    if backend == "abricate":
        # abricate defaults to 80/80; match the PlasmidFinder thresholds instead.
        return [*prefix, "--db", "plasmidfinder",
                "--minid", str(round(float(PLASMIDFINDER_MIN_IDENTITY) * 100)),
                "--mincov", str(round(float(PLASMIDFINDER_MIN_COVERAGE) * 100)), str(assembly)]
    cmd = [*prefix, "-i", str(assembly), "-o", str(out_dir),
           "-l", PLASMIDFINDER_MIN_COVERAGE, "-t", PLASMIDFINDER_MIN_IDENTITY]
    cmd += ["-x"] if backend == "plasmidfinder2" else ["-q", "-j", str(out_dir / "plasmidfinder.json")]
    if db:
        cmd += ["-p", str(db)]
    return cmd


def check_organism(organism: str, available: set[str]) -> None:
    """Reject an --organism the installed AMRFinderPlus database does not know."""
    if available and organism not in available:
        raise ValueError(
            f"'{organism}' is not an AMRFinderPlus organism in the installed database. "
            "Valid values: " + ", ".join(sorted(available))
        )


def check_plasmid_db(backend: str, db: Path | None, env: dict) -> None:
    """PlasmidFinder 3 has no default database location; fail before any tool runs."""
    if backend == "plasmidfinder3" and not db and not env.get("CGE_PLASMIDFINDER_DB"):
        raise ValueError(
            "PlasmidFinder 3 needs its database location: pass --plasmidfinder-db <dir> "
            "or set CGE_PLASMIDFINDER_DB"
        )


def _amrfinder_organisms() -> set[str]:
    try:
        listing = subprocess.run([*tool_cmd("amrfinder"), "--list_organisms"], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return set()
    return parse_organism_listing(listing.stdout + listing.stderr)


def require_tools(names: list[str]) -> None:
    missing = [n for n in names if shutil.which(tool_cmd(n)[0]) is None]
    if missing:
        lines = [f"  {n}: {_INSTALL_HINTS.get(n, 'install it and put it on PATH')}" for n in missing]
        raise MissingToolError(
            "Required tools not found on PATH:\n" + "\n".join(lines) + "\n"
            "Or supply results you already have with --mlst / --amrfinder / --plasmidfinder."
        )


def _run_cmd(cmd: list[str], stdout_path: Path | None = None) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr.strip() or proc.stdout.strip()).splitlines()[-8:])
        raise RuntimeError(f"{cmd[0]} failed (exit {proc.returncode}):\n{tail}")
    if stdout_path is not None:
        stdout_path.write_text(proc.stdout, encoding="utf-8")
    return proc


def _version(cmd: list[str]) -> str | None:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (proc.stdout or proc.stderr).strip()
    return text.splitlines()[0] if text else None


def run_tools(
    assembly: Path,
    raw_dir: Path,
    *,
    need: set[str],
    organism: str | None,
    threads: int,
    plasmidfinder_db: Path | None,
    paths: dict[str, Path],
    mlst_scheme: str | None = None,
) -> dict:
    """Run whichever of mlst / amrfinder / plasmidfinder is in `need`.

    Fills `paths` with the raw result file for each tool run. Returns the organism
    actually passed to AMRFinderPlus, tool versions, any unresolved mlst scheme tie and
    the scheme forced on mlst. mlst runs first because its scheme decides the organism.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    work = assembly
    if str(assembly).endswith(".gz"):
        work = raw_dir / "assembly.fasta"
        with gzip.open(assembly, "rb") as src, open(work, "wb") as dst:
            shutil.copyfileobj(src, dst)
    try:
        return _run_tools(work, raw_dir, need=need, organism=organism, threads=threads,
                          plasmidfinder_db=plasmidfinder_db, paths=paths, mlst_scheme=mlst_scheme)
    finally:
        if work != assembly:
            work.unlink(missing_ok=True)


def _run_tools(work: Path, raw_dir: Path, *, need: set[str], organism: str | None, threads: int,
               plasmidfinder_db: Path | None, paths: dict[str, Path], mlst_scheme: str | None) -> dict:
    versions: dict[str, str | None] = {}
    ties: list[tuple[str, str]] = []
    forced = mlst_scheme
    backend = find_plasmid_backend() if "plasmidfinder" in need else None
    if backend:
        check_plasmid_db(backend[0], plasmidfinder_db, os.environ)
    if "amrfinder" in need and organism:
        check_organism(organism, _amrfinder_organisms())

    if "mlst" in need:
        paths["mlst"] = raw_dir / "mlst.tsv"
        proc = _run_cmd(build_mlst_cmd(work, mlst_scheme), stdout_path=paths["mlst"])
        ties = [] if mlst_scheme else parse_mlst_ties(proc.stderr)
        pick = resolve_tie(ties, organism)
        if pick:
            _run_cmd(build_mlst_cmd(work, pick), stdout_path=paths["mlst"])
            forced, ties = pick, []
        versions["mlst"] = (_version([*tool_cmd("mlst"), "--version"]) or "").removeprefix("mlst ").strip() or None

    applied = organism
    if "amrfinder" in need:
        # Infer only from an mlst run made here with no scheme tie: a supplied mlst file
        # cannot show ties, and a tied pick may belong to another genus.
        if applied is None and "mlst" in need and not ties:
            applied = infer_organism(parse_mlst(paths["mlst"])["scheme"])
            if applied:
                available = _amrfinder_organisms()
                if available and applied not in available:
                    print(f"WARNING: this AMRFinderPlus database has no organism '{applied}'; "
                          "running without --organism", file=sys.stderr)
                    applied = None
        paths["amrfinder"] = raw_dir / "amrfinder.tsv"
        proc = _run_cmd(build_amrfinder_cmd(work, paths["amrfinder"], organism=applied, threads=threads))
        versions["amrfinder"] = _version([*tool_cmd("amrfinder"), "--version"])
        db = re.search(r"Database version:\s*(\S+)", proc.stderr)
        versions["amrfinder_database"] = db.group(1) if db else None

    if "plasmidfinder" in need:
        backend, prefix = backend
        pf_dir = raw_dir / "plasmidfinder"
        pf_dir.mkdir(exist_ok=True)
        cmd = build_plasmid_cmd(backend, prefix, work, pf_dir, db=plasmidfinder_db)
        if backend == "abricate":
            paths["plasmidfinder"] = raw_dir / "plasmidfinder.tsv"
            _run_cmd(cmd, stdout_path=paths["plasmidfinder"])
            versions["plasmidfinder"] = _version([*prefix, "--version"])
        elif backend == "plasmidfinder3":
            paths["plasmidfinder"] = raw_dir / "plasmidfinder.json"
            _run_cmd(cmd)
            shutil.copyfile(pf_dir / "plasmidfinder.json", paths["plasmidfinder"])
            versions["plasmidfinder"] = "PlasmidFinder " + (_version([*prefix, "-v"]) or "3.x")
        else:
            paths["plasmidfinder"] = raw_dir / "plasmidfinder.tsv"
            _run_cmd(cmd)
            shutil.copyfile(pf_dir / "results_tab.tsv", paths["plasmidfinder"])
            versions["plasmidfinder"] = "PlasmidFinder 2.x (plasmidfinder.py reports no version)"
        shutil.rmtree(pf_dir, ignore_errors=True)
    return {"applied_organism": applied if "amrfinder" in need else None, "versions": versions,
            "mlst_ties": ties, "mlst_scheme_forced": forced}


# ── Outputs ────────────────────────────────────────────────────────────────────


def _write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(["" if v is None else v for v in row] for row in rows)


def _write_tables(result: dict, tables: Path) -> list[Path]:
    amr = result["amr"]
    hits = amr.get("genes", []) + amr.get("point_mutations", []) + amr.get("other_elements", [])
    fields = ["symbol", "name", "element_type", "subtype", "drug_class", "subclass", "method",
              "call", "tier", "scope", "contig", "start", "stop", "strand", "coverage", "identity"]
    _write_csv(tables / "amr_determinants.csv",
               fields + ["replicons_on_contig", "nearest_replicon", "nearest_replicon_distance_bp", "contig_length"],
               [[h[f] for f in fields] + [";".join(h["replicons_on_contig"]),
                                          (h["nearest_replicon"] or {}).get("replicon"),
                                          (h["nearest_replicon"] or {}).get("distance_bp"),
                                          h["contig_length"]] for h in hits])
    _write_csv(tables / "drug_class_summary.csv", ["drug_class", "confident", "review", "disrupted"],
               [[c["drug_class"], ";".join(c["confident"]), ";".join(c["review"]), ";".join(c["disrupted"])]
                for c in amr.get("class_summary", [])])
    mlst = result["mlst"]
    _write_csv(tables / "mlst.csv", ["sample", "scheme", "st", "call", "locus", "allele", "flag"],
               [[result["sample"], mlst.get("scheme"), mlst.get("st"), mlst.get("call"),
                 a["locus"], a["allele"], a["flag"]] for a in mlst.get("alleles", [])])
    rep_fields = ["replicon", "contig", "start", "stop", "identity", "coverage", "accession"]
    _write_csv(tables / "plasmid_replicons.csv", rep_fields,
               [[r[f] for f in rep_fields] for r in result["plasmids"].get("replicons", [])])
    return [tables / n for n in ("amr_determinants.csv", "drug_class_summary.csv",
                                 "mlst.csv", "plasmid_replicons.csv")]


def _write_figure(result: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = result["amr"].get("class_summary", [])
    fig, ax = plt.subplots(figsize=(7, max(2.2, 0.42 * len(summary) + 1.2)))
    if summary:
        labels = [c["drug_class"].title() for c in summary][::-1]
        left = [0] * len(summary)
        for tier, colour in (("confident", "#2a6f97"), ("review", "#e9c46a"), ("disrupted", "#b0b0b0")):
            values = [len(c[tier]) for c in summary][::-1]
            ax.barh(labels, values, left=left, color=colour, label=tier)
            left = [a + b for a, b in zip(left, values)]
        ax.set_xlabel("AMR determinants (genes and point mutations)")
        ax.xaxis.get_major_locator().set_params(integer=True)
        ax.legend(frameon=False, fontsize=8)
    else:
        message = ("AMRFinderPlus not run" if result["amr"]["status"] == "not_run"
                   else "No AMR determinants reported")
        ax.text(0.5, 0.5, message, ha="center", va="center")
        ax.set_axis_off()
    ax.set_title(f"{result['sample']}: AMR determinants by drug class", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _fmt(value) -> str:
    return "-" if value in (None, "") else str(value)


def _fmt_bp(n: int | None) -> str:
    if n is None:
        return "-"
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f} Mb"
    return f"{n / 1000:.1f} kb" if n >= 1000 else f"{n} bp"


def _hit_table(hits: list[dict]) -> list[str]:
    lines = ["| Determinant | Class | Subclass | Call | Tier | % id | % cov | Contig | Contig length "
             "| Nearest replicon on contig (distance) |",
             "|---|---|---|---|---|---:|---:|---|---:|---|"]
    for h in hits:
        lines.append(
            f"| {h['symbol']} | {_fmt(h['drug_class'])} | {_fmt(h['subclass'])} | {h['method']} | "
            f"{h['tier']} | {_fmt(h['identity'])} | {_fmt(h['coverage'])} | {_fmt(h['contig'])} | "
            f"{_fmt_bp(h['contig_length'])} | "
            + (f"{h['nearest_replicon']['replicon']} ({_fmt_bp(h['nearest_replicon']['distance_bp'])})"
               if h["nearest_replicon"] else ", ".join(h["replicons_on_contig"]) or "-") + " |"
        )
    return lines


def build_report(result: dict) -> str:
    amr, mlst, plasmids, asm = result["amr"], result["mlst"], result["plasmids"], result["assembly"]
    mode_text = {
        "demo": "Demo: synthetic, hand-written tool outputs (not a real isolate)",
        "live": "Tools run by this skill on the supplied assembly",
        "precomputed": "Pre-computed tool outputs supplied by the user",
        "mixed": "Some tools run by this skill, some outputs supplied by the user",
    }[result["mode"]]
    out = [f"# Isolate AMR and Typing Report: {result['sample']}", "", f"**Mode**: {mode_text}", ""]

    out += ["## Summary", ""]
    if mlst["status"] == "not_run":
        out.append("- **Sequence type**: NOT RUN")
    elif mlst["call"] == "assigned":
        out.append(f"- **Sequence type**: ST{mlst['st']} (scheme `{mlst['scheme']}`)")
    elif mlst["call"] == "ambiguous":
        tied = " or ".join(f"{t['scheme']} ST{t['st']}" for t in mlst["tied_schemes"])
        out.append(f"- **Sequence type**: AMBIGUOUS, {tied} (equal mlst score; see MLST)")
    elif mlst["call"] == "unassigned":
        out.append(f"- **Sequence type**: not assigned (scheme `{mlst['scheme']}`). {mlst['note']}")
    else:
        out.append(f"- **Sequence type**: {mlst['note']}")
    if amr["status"] == "not_run":
        out.append("- **AMR determinants**: NOT RUN")
    else:
        n_conf = sum(1 for h in amr["genes"] if h["tier"] == "confident")
        n_other = len(amr["genes"]) - n_conf
        classes = [c["drug_class"].title() for c in amr["class_summary"] if c["confident"]]
        out.append(f"- **AMR genes**: {n_conf} confident, {n_other} needing review or disrupted")
        point_text = {
            "detected": f"{len(amr['point_mutations'])} detected",
            "screened": "screened, none detected",
            "not_screened": "NOT SCREENED (no organism; see Limits)",
            "unknown": "screening status unknown (see Limits)",
        }[amr["point_mutation_screening"]]
        out.append(f"- **Resistance point mutations**: {point_text}")
        out.append(f"- **Drug classes with a confident determinant**: {', '.join(classes) or 'none'}")
    if plasmids["status"] == "not_run":
        out.append("- **Plasmid replicons**: NOT RUN")
    else:
        names = [r["replicon"] for r in plasmids["replicons"]]
        out.append(f"- **Plasmid replicons**: {', '.join(names) or 'none detected'}")
    out.append("")

    if asm:
        out += ["## Assembly", "",
                "| File | Contigs | Total length (bp) | N50 (bp) | Longest contig (bp) | GC % |",
                "|---|---:|---:|---:|---:|---:|",
                f"| {asm['file']} | {asm['n_contigs']} | {asm['total_length']} | {asm['n50']} | "
                f"{asm['longest_contig']} | {_fmt(asm['gc_percent'])} |", ""]

    out += ["## MLST", ""]
    if mlst["status"] == "not_run":
        out += ["NOT RUN: no mlst result was produced or supplied.", ""]
    else:
        out += [f"Scheme: `{_fmt(mlst['scheme'])}`. ST: {_fmt(mlst['st'])}. {mlst['note']}".rstrip(), ""]
        if mlst.get("rerun_options"):
            out += ["The scheme and ST above are mlst's arbitrary pick. Rerun with the flags that fit this isolate:", ""]
            out += [f"- `{o['flags']}` for scheme `{o['scheme']}` (ST{o['st']})"
                    + (f", organism {o['organism']}" if o["organism"] else "")
                    for o in mlst["rerun_options"]]
            out.append("")
        if mlst["alleles"]:
            out += ["| Locus | Allele | Match |", "|---|---|---|"]
            out += [f"| {a['locus']} | {_fmt(a['allele'])} | {a['flag']} |" for a in mlst["alleles"]]
            out.append("")

    out += ["## AMR genes", ""]
    if amr["status"] == "not_run":
        out += ["NOT RUN: no AMRFinderPlus result was produced or supplied.", ""]
    else:
        if not amr["organism"]:
            org = "none"
        elif amr["organism_source"] == "user":
            org = f"`{amr['organism']}` (stated by user)"
        else:
            org = f"`{amr['organism']}` (inferred from the MLST scheme and applied)"
        out += [f"AMRFinderPlus organism: {org}.", ""]
        out += (_hit_table(amr["genes"]) if amr["genes"] else ["No AMR genes reported."]) + [""]
        out += ["## Resistance point mutations", ""]
        out += (_hit_table(amr["point_mutations"]) if amr["point_mutations"]
                else ["None reported. Read the screening status in the Summary before treating this as negative."])
        out.append("")
        out += ["## Drug class summary", ""]
        if amr["class_summary"]:
            out += ["| Drug class | Confident | Needs review | Disrupted |", "|---|---|---|---|"]
            out += [f"| {c['drug_class']} | {', '.join(c['confident']) or '-'} | "
                    f"{', '.join(c['review']) or '-'} | {', '.join(c['disrupted']) or '-'} |"
                    for c in amr["class_summary"]]
        else:
            out.append("No AMR determinants reported.")
        out.append("")
        if amr["other_elements"]:
            out += ["## Stress and virulence elements", "",
                    "Reported by AMRFinderPlus `--plus`. These are not antimicrobial resistance "
                    "determinants and are excluded from the drug class summary.", "",
                    "| Element | Type | Subtype | Class | Call | Contig |", "|---|---|---|---|---|---|"]
            out += [f"| {h['symbol']} | {_fmt(h['element_type'])} | {_fmt(h['subtype'])} | "
                    f"{_fmt(h['drug_class'])} | {h['method']} | {_fmt(h['contig'])} |"
                    for h in amr["other_elements"]]
            out.append("")

    out += ["## Plasmid replicons", ""]
    if plasmids["status"] == "not_run":
        out += ["NOT RUN: no PlasmidFinder result was produced or supplied.", ""]
    elif plasmids["replicons"]:
        out += ["| Replicon | % id | % cov | Contig | Position |", "|---|---:|---:|---|---|"]
        out += [f"| {r['replicon']} | {_fmt(r['identity'])} | {_fmt(r['coverage'])} | {_fmt(r['contig'])} | "
                f"{_fmt(r['start'])}..{_fmt(r['stop'])} |" for r in plasmids["replicons"]]
        out.append("")
    else:
        out += ["No replicons detected.", ""]

    out += ["## Limits", ""]
    out += [f"- {w}" for w in result["warnings"]]
    out += [
        "- This is a genotype. It is not a phenotypic susceptibility result: a detected "
        "determinant does not prove resistance, and no detected determinant does not prove "
        "susceptibility.",
        "- The gene list includes any species-intrinsic genes AMRFinderPlus reports (for example "
        "chromosomal blaEC in E. coli); the tools do not flag them as intrinsic.",
        "- \"Nearest replicon on contig\" is the closest replicon hit on the same contig and the "
        "gap to it, measured along the contig as if linear. It is a distance, not a verdict: a "
        "blank cell does not place a gene on the chromosome (plasmids fragment in short-read "
        "assemblies), and a replicon megabases away on a closed chromosome does not place a gene "
        "on a plasmid.",
        "- Calls are the tools' own; this skill does not re-filter by identity or coverage.",
        "",
        "## Methods", "",
        "- MLST: `mlst` (Seemann) against PubMLST schemes.",
        "- AMR: NCBI AMRFinderPlus with `--plus`; point mutations require `--organism`.",
        f"- Plasmid replicons: PlasmidFinder database (live runs: identity >= {PLASMIDFINDER_MIN_IDENTITY}, "
        f"coverage >= {PLASMIDFINDER_MIN_COVERAGE}; supplied results keep the thresholds they were run with).",
        "- Tiers: confident = ALLELE, EXACT, BLAST, POINT; review = PARTIAL, PARTIAL_CONTIG_END, HMM; "
        "disrupted = INTERNAL_STOP.",
    ]
    versions = {k: v for k, v in result["tool_versions"].items() if v}
    if versions:
        out.append("- Versions: " + "; ".join(f"{k} {v}" for k, v in sorted(versions.items())) + ".")
    out += ["", f"*{DISCLAIMER}*", ""]
    return "\n".join(out)


def write_outputs(result: dict, output_dir: Path, raw: dict[str, Path], repro_args: list,
                  preflight: list[str]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for sub in ("tables", "figures", "tool_outputs"):
        (output_dir / sub).mkdir(exist_ok=True)

    written = [output_dir / "report.md", output_dir / "result.json"]
    for tool, src in sorted(raw.items()):
        dest = output_dir / "tool_outputs" / f"{tool}{'.json' if src.suffix == '.json' else '.tsv'}"
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
        written.append(dest)

    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (output_dir / "report.md").write_text(build_report(result), encoding="utf-8")
    written += _write_tables(result, output_dir / "tables")
    figure = output_dir / "figures" / "drug_class_determinants.png"
    _write_figure(result, figure)
    written.append(figure)

    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills") / SKILL_NAME / "isolate_amr_typing.py",
            args=repro_args + ["--output", ReproPath(output_dir, anchor="output_dir")],
            comment=f"Reproduce this {SKILL_NAME} run",
            preflight=preflight,
        ),
        repo_root=_PROJECT_ROOT,
    )
    # Two environments: current AMRFinderPlus and mlst do not solve into one.
    channels = ["conda-forge", "bioconda"]
    write_environment_yml(
        output_dir,
        f"clawbio-{SKILL_NAME}",
        pip_deps=[],
        conda_deps=["matplotlib", "ncbi-amrfinderplus", "plasmidfinder"],
        channels=channels,
    )
    with tempfile.TemporaryDirectory() as tmp:
        mlst_env = write_environment_yml(tmp, f"clawbio-{SKILL_NAME}-mlst", pip_deps=[],
                                         conda_deps=["mlst"], channels=channels)
        shutil.copyfile(mlst_env, output_dir / "reproducibility" / "environment-mlst.yml")
    write_checksums(written, output_dir, anchor=output_dir)


# ── CLI ────────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="MLST, AMR genes, resistance point mutations and plasmid replicons for one "
                    "bacterial isolate assembly."
    )
    parser.add_argument("--input", type=Path, help="Assembled genome, nucleotide FASTA (.gz accepted)")
    parser.add_argument("--output", type=Path, default=Path("isolate_amr_typing_out"), help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run on bundled synthetic data; no tools needed")
    parser.add_argument("--mlst", type=Path, help="Existing `mlst` output; skips running mlst")
    parser.add_argument("--amrfinder", type=Path, help="Existing AMRFinderPlus TSV; skips running amrfinder")
    parser.add_argument("--plasmidfinder", type=Path,
                        help="Existing PlasmidFinder results_tab.tsv (v2), -j JSON (v3) or abricate TSV; skips running it")
    parser.add_argument("--organism",
                        help="AMRFinderPlus --organism value (default: inferred from the MLST scheme)")
    parser.add_argument("--mlst-scheme",
                        help="Force this mlst scheme instead of auto-detection (resolves scheme ties)")
    parser.add_argument("--sample-name", help="Label for the report (default: assembly file name)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for AMRFinderPlus (default 4)")
    parser.add_argument("--plasmidfinder-db", type=Path, help="PlasmidFinder database directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.demo:
        args.input, args.mlst, args.amrfinder, args.plasmidfinder = DEMO_FASTA, DEMO_MLST, DEMO_AMR, DEMO_PLASMID
    supplied = {k: v for k, v in (("mlst", args.mlst), ("amrfinder", args.amrfinder),
                                  ("plasmidfinder", args.plasmidfinder)) if v}
    if args.input is None and not supplied:
        parser.error("provide --input <assembly.fasta>, pre-computed tool outputs, or --demo")

    for path in [args.input, *supplied.values()]:
        if path is not None and not path.is_file():
            print(f"ERROR: input file not found: {path}", file=sys.stderr)
            return 2

    need = set() if args.input is None else {"mlst", "amrfinder", "plasmidfinder"} - supplied.keys()
    if args.demo:
        mode = "demo"
    elif not need:
        mode = "precomputed"
    else:
        mode = "mixed" if supplied else "live"

    output_dir = args.output
    try:
        assembly = read_assembly(args.input) if args.input else None
        if need:
            tools = sorted(need - {"plasmidfinder"})
            if "plasmidfinder" in need and find_plasmid_backend() is None:
                tools.append("plasmidfinder.py")
            require_tools(tools)

        if output_dir.exists() and any(output_dir.iterdir()):
            print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}",
                  file=sys.stderr)

        raw = dict(supplied)
        run = {"applied_organism": None, "versions": {}, "mlst_ties": [], "mlst_scheme_forced": None}
        if need:
            run = run_tools(
                args.input, output_dir / "tool_outputs", need=need, organism=args.organism,
                threads=args.threads, plasmidfinder_db=args.plasmidfinder_db, paths=raw,
                mlst_scheme=args.mlst_scheme,
            )
        result = summarise(
            sample=args.sample_name or (sample_name_from(args.input) if args.input else "isolate"),
            assembly=assembly,
            amr_hits=parse_amrfinder(raw["amrfinder"]) if "amrfinder" in raw else None,
            mlst=parse_mlst(raw["mlst"]) if "mlst" in raw else None,
            replicons=parse_plasmidfinder(raw["plasmidfinder"]) if "plasmidfinder" in raw else None,
            organism=args.organism,
            mode=mode,
            applied_organism=run["applied_organism"],
            amr_run_here="amrfinder" in need,
            mlst_ties=run["mlst_ties"],
            mlst_run_here="mlst" in need,
            tool_versions=run["versions"],
        )
    except MissingToolError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.demo:
        repro_args: list = ["--demo"]
    else:
        repro_args = []
        for flag, path in (("--input", args.input), ("--mlst", args.mlst), ("--amrfinder", args.amrfinder),
                           ("--plasmidfinder", args.plasmidfinder), ("--plasmidfinder-db", args.plasmidfinder_db)):
            if path is not None:
                repro_args += [flag, ReproPath(path.resolve())]
        if args.organism:
            repro_args += ["--organism", args.organism]
        if run["mlst_scheme_forced"]:
            repro_args += ["--mlst-scheme", run["mlst_scheme_forced"]]
        if args.sample_name:
            repro_args += ["--sample-name", f'"{args.sample_name}"']
        if need:
            repro_args += ["--threads", str(args.threads)]
    preflight = override_preflight(os.environ) if need else []
    for tool in sorted(need - {"plasmidfinder"}):
        exe = tool_cmd(tool)[0]
        preflight.append(f'command -v {exe} >/dev/null || {{ echo "{exe} not on PATH" >&2; exit 1; }}')

    write_outputs(result, output_dir, raw, repro_args, preflight)
    options = result["mlst"].get("rerun_options")
    if options and args.input:
        def kept(tool: str) -> str:
            suffix = ".json" if raw[tool].suffix == ".json" else ".tsv"
            return f" --{tool} {output_dir / 'tool_outputs' / (tool + suffix)}"

        print(f"WARNING: mlst could not choose between schemes for {result['sample']}. Rerun with one of:",
              file=sys.stderr)
        for opt in options:
            # The plasmid result never depends on the organism; the AMR result does only
            # when the scheme has one to apply.
            reuse = (kept("plasmidfinder") if "plasmidfinder" in raw else "") + (
                kept("amrfinder") if "amrfinder" in raw and not opt["organism"] else "")
            print(f"  python {Path(__file__).resolve()} --input {args.input} {opt['flags']}{reuse} "
                  f"--output {output_dir}_{opt['scheme']}", file=sys.stderr)
    print(f"Isolate AMR Typing wrote {output_dir / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
