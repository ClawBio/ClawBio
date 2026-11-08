#!/usr/bin/env python3
"""Stage 1 of the gwas_prs DAG: resolve PGS Catalog scoring file(s).

Port of ClawBio 1.x `gwas-prs` (skills/gwas-prs/gwas_prs.py) for the
ClawBio 2.0 `clawbio` CLI. Given a comma-separated list of PGS IDs (or the
literal "all" for the 6 bundled curated scores), locates a local scoring
file for each — from this skill's bundled data/ directory first, falling
back to a live PGS Catalog fetch (stdlib `urllib`, no `requests` dependency)
for any ID not bundled. Writes a manifest JSON consumed by
scripts/calculate_prs.py.

This is the "score-file resolution" stage; genotype parsing
(scripts/parse_genotype.py) is a separate, independent stage — neither
depends on the other's output, so they run as DAG siblings (see
pipeline.yaml). This mirrors the resolve->fetch split in gwas_lookup and
the independent-stage split in claw_metagenomics.

Network fetches degrade gracefully: on any failure (no network, unknown
PGS ID, HTTP error) the manifest records an "error" entry for that ID
instead of crashing the run, so `calculate_prs` can skip it with a warning
— the same "cache-first, graceful stub on failure" pattern gwas_lookup uses.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BUILD = "GRCh37"
USER_AGENT = "ClawBio-GWAS-PRS/2.0 (clawbio2.0 port)"
PGS_API_BASE = "https://www.pgscatalog.org/rest"
PGS_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"


CURATED_PANEL_MARKER = "#clawbio_panel=curated_demo"


def is_curated_demo_panel(filepath: str | Path) -> bool:
    """True if the file's leading comment block carries the curated-demo marker.

    Same contract as 1.x `gwas_prs.is_curated_demo_panel`: reads only the
    leading `#` block and is safe on a missing or unreadable path.
    """
    path = Path(filepath)
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:  # type: ignore[operator]
            for line in handle:
                if not line.startswith("#"):
                    return False
                if line.strip() == CURATED_PANEL_MARKER:
                    return True
    except (OSError, EOFError, UnicodeDecodeError):
        return False
    return False


def _load_curated(skill_dir: Path) -> dict:
    path = skill_dir / "curated_scores.json"
    return json.loads(path.read_text())


def _http_get_json(url: str, timeout: int = 30) -> dict:
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def _fetch_live(pgs_id: str, cache_dir: Path) -> Path:
    """Fetch metadata + harmonized scoring file from the live PGS Catalog.

    Downloads the gzipped harmonized file from EBI FTP and decompresses it
    into cache_dir. Raises on any network/HTTP error — caller catches.
    """
    # Confirm the score exists (also surfaces a clean 404 for bad IDs).
    _http_get_json(f"{PGS_API_BASE}/score/{pgs_id}")

    url = f"{PGS_FTP_BASE}/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_{BUILD}.txt.gz"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    cache_dir.mkdir(parents=True, exist_ok=True)
    gz_path = cache_dir / f"{pgs_id}_hmPOS_{BUILD}.txt.gz"
    with urllib.request.urlopen(req, timeout=60) as resp:
        gz_path.write_bytes(resp.read())

    txt_path = cache_dir / f"{pgs_id}_hmPOS_{BUILD}.txt"
    with gzip.open(gz_path, "rt", errors="replace") as fh:
        txt_path.write_text(fh.read())
    return txt_path


def resolve(pgs_ids: list[str], skill_dir: Path, cache_dir: Path) -> list[dict]:
    curated = _load_curated(skill_dir)
    data_dir = skill_dir / "data"

    manifest: list[dict] = []
    for pgs_id in pgs_ids:
        pgs_id = pgs_id.strip().upper()
        if not pgs_id:
            continue

        meta = curated.get(pgs_id, {})
        local_path = data_dir / f"{pgs_id}_hmPOS_{BUILD}.txt"
        if not local_path.exists():
            # 1.x curated panels are named <PANEL-ID>_GRCh37.txt
            local_path = data_dir / f"{pgs_id}_{BUILD}.txt"

        if local_path.exists():
            curated_flag = is_curated_demo_panel(local_path)
            manifest.append({
                "pgs_id": pgs_id,
                "trait": meta.get("trait", "Unknown"),
                "filepath": str(local_path),
                "metadata": {
                    "publication": meta.get("publication", ""),
                    "variants_count": meta.get("variants_count", 0),
                },
                "reference_distribution": meta.get("reference_distribution"),
                "source": "bundled",
                "curated_demo_panel": curated_flag,
                "curated_panel_id": (
                    meta.get("curated_panel_id") or pgs_id if curated_flag else None
                ),
                "legacy_pgs_id": meta.get("legacy_pgs_id") if curated_flag else None,
            })
            print(f"  {pgs_id}: bundled scoring file ({meta.get('trait', 'Unknown')})")
            continue

        print(f"  {pgs_id}: not bundled, attempting live PGS Catalog fetch...")
        try:
            fetched_path = _fetch_live(pgs_id, cache_dir)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            print(f"  {pgs_id}: WARNING fetch failed ({e}) — will be skipped downstream")
            manifest.append({
                "pgs_id": pgs_id,
                "trait": meta.get("trait", "Unknown"),
                "filepath": None,
                "metadata": {},
                "reference_distribution": meta.get("reference_distribution"),
                "source": "error",
                "curated_demo_panel": False,
                "curated_panel_id": None,
                "legacy_pgs_id": None,
                "error": str(e),
            })
            continue

        manifest.append({
            "pgs_id": pgs_id,
            "trait": meta.get("trait", "Unknown"),
            "filepath": str(fetched_path),
            "metadata": {
                "publication": meta.get("publication", ""),
                "variants_count": meta.get("variants_count", 0),
            },
            "reference_distribution": meta.get("reference_distribution"),
            "source": "fetched",
            "curated_demo_panel": False,
            "curated_panel_id": None,
            "legacy_pgs_id": None,
        })
        print(f"  {pgs_id}: fetched live scoring file")

    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Resolve PGS Catalog scoring file(s)")
    parser.add_argument(
        "--pgs-ids", required=True,
        help='Comma-separated PGS IDs (e.g. "PGS000013,PGS000011"), or "all" '
             "for the 6 bundled curated scores",
    )
    parser.add_argument("--skill-dir", required=True, help="Path to the gwas_prs skill directory")
    parser.add_argument("--out", required=True, help="Output path for scoring_manifest.json")
    args = parser.parse_args()

    skill_dir = Path(args.skill_dir).resolve()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = out_path.parent / "pgs_cache"

    curated = _load_curated(skill_dir)
    if args.pgs_ids.strip().lower() == "all":
        pgs_ids = sorted(curated.keys())
    else:
        pgs_ids = [p.strip() for p in args.pgs_ids.split(",") if p.strip()]

    if not pgs_ids:
        print("Error: no PGS IDs given (use --pgs-ids all or a comma-separated list)")
        sys.exit(1)

    print(f"Resolving {len(pgs_ids)} scoring file(s): {', '.join(pgs_ids)}")
    manifest = resolve(pgs_ids, skill_dir, cache_dir)

    resolved = [m for m in manifest if m["source"] != "error"]
    if not resolved:
        print("Error: no scoring files could be resolved (all bundled lookups and live fetches failed)")
        sys.exit(1)

    out_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nResolved {len(resolved)}/{len(manifest)} scoring file(s) -> {out_path}")


if __name__ == "__main__":
    main()
