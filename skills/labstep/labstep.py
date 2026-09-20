#!/usr/bin/env python3
"""
labstep.py — Labstep ELN bridge for ClawBio
============================================
Query experiments, protocols, resources, and inventory via the Labstep API
(labstepPy).  Run with --demo for an offline showcase using synthetic data.

Usage:
    python skills/labstep/labstep.py --demo --output /tmp/labstep
    python skills/labstep/labstep.py --experiments [--search QUERY] [--count N] [--output DIR]
    python skills/labstep/labstep.py --experiment-id ID [--output DIR]
    python skills/labstep/labstep.py --protocols [--search QUERY] [--count N] [--output DIR]
    python skills/labstep/labstep.py --protocol-id ID [--output DIR]
    python skills/labstep/labstep.py --inventory [--search QUERY] [--output DIR]
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DIR = SKILL_DIR / "demo"
_PROJECT_ROOT = SKILL_DIR.parent.parent

if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions.*"
)


# ---------------------------------------------------------------------------
# Authentication
# ---------------------------------------------------------------------------


def get_labstep_user():
    """Authenticate with Labstep; returns a labstep User object."""
    try:
        import labstep  # type: ignore
    except ImportError:
        print(
            "ERROR: labstepPy not installed. Run: pip install labstep",
            file=sys.stderr,
        )
        sys.exit(1)

    key = os.environ.get("LABSTEP_API_KEY")
    if not key:
        settings_path = Path(".claude/settings.json")
        if settings_path.exists():
            cfg = json.loads(settings_path.read_text(encoding="utf-8"))
            key = cfg.get("skillsConfig", {}).get("labstep", {}).get("apiKey")

    if not key:
        print(
            "ERROR: No Labstep API key found.\n"
            "  Set the LABSTEP_API_KEY environment variable, or configure\n"
            "  .claude/settings.json → skillsConfig.labstep.apiKey",
            file=sys.stderr,
        )
        sys.exit(1)

    return labstep.authenticate(apikey=key)


# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------


def _fmt_date(iso: str) -> str:
    """Convert ISO datetime string to YYYY-MM-DD."""
    try:
        dt = datetime.fromisoformat(iso.replace("Z", "+00:00"))
        return dt.strftime("%Y-%m-%d")
    except Exception:
        return iso[:10] if iso else "—"


def format_experiments(experiments: list[dict], title: str = "Experiments") -> str:
    lines = [f"# 🔬 Labstep — {title}\n", f"**{len(experiments)} experiment(s)**\n"]
    for exp in experiments:
        lines.append(f"## [{exp['id']}] {exp['name']}\n")
        lines.append(f"- **Created**: {_fmt_date(exp.get('created_at', ''))}")
        lines.append(f"- **Updated**: {_fmt_date(exp.get('updated_at', ''))}")

        tags = exp.get("tags", [])
        if tags:
            tag_str = ", ".join(f"`{t['name']}`" for t in tags)
            lines.append(f"- **Tags**: {tag_str}")

        data_fields = exp.get("data_fields", [])
        if data_fields:
            lines.append("\n**Data Fields:**\n")
            lines.append("| Field | Value |")
            lines.append("|---|---|")
            for df in data_fields:
                lines.append(f"| {df['label']} | {df.get('value', '—')} |")

        protocols = exp.get("protocols", [])
        if protocols:
            proto_str = ", ".join(f"{p['name']} (#{p['id']})" for p in protocols)
            lines.append(f"\n**Linked Protocols**: {proto_str}")

        comments = exp.get("comments", [])
        if comments:
            lines.append("\n**Comments:**\n")
            for c in comments:
                lines.append(f"- *{_fmt_date(c.get('created_at', ''))}* — {c['text']}")

        lines.append("")

    lines.append(f"---\n{DISCLAIMER}")
    return "\n".join(lines)


def format_protocols(protocols: list[dict], title: str = "Protocols") -> str:
    lines = [f"# 📋 Labstep — {title}\n", f"**{len(protocols)} protocol(s)**\n"]
    for proto in protocols:
        lines.append(f"## [{proto['id']}] {proto['name']}  (v{proto.get('version', 1)})\n")
        lines.append(f"- **Created**: {_fmt_date(proto.get('created_at', ''))}")
        lines.append(f"- **Updated**: {_fmt_date(proto.get('updated_at', ''))}")

        steps = proto.get("steps", [])
        if steps:
            lines.append(f"\n**Steps ({len(steps)}):**\n")
            for s in steps:
                lines.append(f"**{s['position']}. {s.get('title', 'Step')}**\n")
                lines.append(f"{s.get('body', '')}\n")

        inv = proto.get("inventory_fields", [])
        if inv:
            inv_str = ", ".join(f"{f['label']} (#{f['resource_id']})" for f in inv)
            lines.append(f"**Inventory Fields**: {inv_str}")

        lines.append("")

    lines.append(f"---\n{DISCLAIMER}")
    return "\n".join(lines)


def format_inventory(data: dict, search: str | None = None) -> str:
    resources = data.get("resources", [])
    locations = data.get("locations", [])

    if search:
        q = search.lower()
        resources = [r for r in resources if q in r["name"].lower() or q in r.get("category", "").lower()]

    title = f"Inventory — \"{search}\"" if search else "Inventory"
    lines = [
        f"# 🧪 Labstep — {title}\n",
        f"**{len(resources)} resource(s)** across {len(locations)} location(s)\n",
    ]

    by_category: dict[str, list] = {}
    for r in resources:
        cat = r.get("category", "Uncategorised")
        by_category.setdefault(cat, []).append(r)

    for cat, items in sorted(by_category.items()):
        lines.append(f"## {cat}\n")
        for r in items:
            meta = r.get("metadata", {})
            supplier = meta.get("supplier", "")
            lot = meta.get("lot", "")
            expiry = meta.get("expiry", "")
            hazard = meta.get("hazard", "")
            n_items = len(r.get("items", []))
            lines.append(f"### [{r['id']}] {r['name']}\n")
            if supplier:
                lines.append(f"- **Supplier/Cat#**: {supplier}")
            if lot:
                lines.append(f"- **Lot**: {lot}")
            if expiry:
                lines.append(f"- **Expiry**: {expiry}")
            if hazard:
                lines.append(f"- **Hazard**: {hazard}")
            lines.append(f"- **Stock items**: {n_items}")

            for item in r.get("items", []):
                lines.append(
                    f"  - Item #{item['id']}: {item.get('name', '—')} | "
                    f"{item.get('amount', '—')} | 📍 {item.get('location', '—')}"
                )
            lines.append("")

    if locations and not search:
        lines.append("## Storage Locations\n")
        lines.append("| Location | Items |")
        lines.append("|---|---|")
        for loc in locations:
            lines.append(f"| {loc['name']} | {loc['item_count']} |")
        lines.append("")

    lines.append(f"---\n{DISCLAIMER}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _write_reproducibility(output_dir: Path, repro_args: list[str]) -> None:
    """Write the reproducibility bundle for a run into output_dir."""
    write_environment_yml(
        output_dir,
        env_name="clawbio-labstep",
        pip_deps=["labstep>=3.0"],
        python_version="3.11",
    )
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/labstep/labstep.py"),
            args=[*repro_args, "--output", ReproPath(output_dir, "output_dir")],
            comment="Reproduce this labstep run",
        ),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [output_dir / "report.md", output_dir / "result.json"],
        output_dir,
        anchor=output_dir,
    )


def _write_output(output_dir: Path, content: str, label: str, repro_args: list[str]) -> None:
    """Write markdown content to output_dir/report.md and a result.json envelope."""
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.md"
    report_path.write_text(content, encoding="utf-8")
    print(f"Report written to {report_path}")

    # Standardised result.json — use common helper when available
    try:
        if str(_PROJECT_ROOT) not in sys.path:
            sys.path.insert(0, str(_PROJECT_ROOT))
        from clawbio.common.report import write_result_json  # type: ignore
        write_result_json(
            output_dir=output_dir,
            skill="labstep",
            version="0.2.0",
            summary={"label": label, "output": str(report_path)},
            data={},
        )
    except ImportError:
        pass  # common helper not available — report.md is sufficient

    _write_reproducibility(output_dir, repro_args)
    print(f"Full output in {output_dir}/")


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------


def _load_demo(filename: str) -> dict | list:
    path = DEMO_DIR / filename
    if not path.exists():
        print(f"ERROR: Demo file missing: {path}", file=sys.stderr)
        sys.exit(1)
    return json.loads(path.read_text(encoding="utf-8"))


def run_demo(output_dir: Path | None = None) -> None:
    print("\nLabstep ELN Bridge — Demo Mode (offline synthetic data)")
    print("=" * 58)

    experiments: list[dict] = _load_demo("demo_experiments.json")  # type: ignore
    protocols: list[dict] = _load_demo("demo_protocols.json")  # type: ignore
    inventory: dict = _load_demo("demo_inventory.json")  # type: ignore

    exp_md = format_experiments(experiments, title="Recent Experiments (demo)")
    proto_md = format_protocols([protocols[0]], title="Protocol Detail (demo)")
    inv_md = format_inventory(inventory)
    inv_search_md = format_inventory(inventory, search="RNA")

    print("\n\n--- EXPERIMENTS ---\n")
    print(exp_md)

    print("\n\n--- PROTOCOL DETAIL: Lentiviral sgRNA Library Transduction ---\n")
    print(proto_md)

    print("\n\n--- INVENTORY SNAPSHOT ---\n")
    print(inv_md)

    print("\n\n--- INVENTORY SEARCH: \"RNA\" ---\n")
    print(inv_search_md)

    if output_dir is not None:
        sections = [
            "# 🔬 Labstep ELN — Demo Report\n",
            "## Experiments\n",
            exp_md,
            "\n## Protocol Detail\n",
            proto_md,
            "\n## Inventory Snapshot\n",
            inv_md,
            "\n## Inventory Search: \"RNA\"\n",
            inv_search_md,
        ]
        _write_output(output_dir, "\n".join(sections), label="demo", repro_args=["--demo"])


# ---------------------------------------------------------------------------
# Live API helpers
# ---------------------------------------------------------------------------


def live_experiments(user, search: str | None, count: int) -> list[dict]:
    kwargs: dict = {"count": count}
    if search:
        kwargs["search_query"] = search
    exps = user.getExperiments(**kwargs)
    results = []
    for e in exps:
        tags = []
        try:
            tags = [{"name": t.name} for t in (e.getTags() or [])]
        except Exception:
            pass
        data_fields = []
        try:
            for df in (e.getDataFields() or []):
                data_fields.append({
                    "label": getattr(df, "label", ""),
                    "field_type": getattr(df, "field_type", "default"),
                    "value": str(getattr(df, "value", "") or ""),
                })
        except Exception:
            pass
        linked_protocols = []
        try:
            for p in (e.getProtocols() or []):
                linked_protocols.append({"id": p.id, "name": p.name})
        except Exception:
            pass
        results.append({
            "id": e.id,
            "name": e.name,
            "created_at": str(getattr(e, "created_at", "")),
            "updated_at": str(getattr(e, "updated_at", "")),
            "tags": tags,
            "data_fields": data_fields,
            "protocols": linked_protocols,
        })
    return results


def live_experiment_detail(user, exp_id: int) -> list[dict]:
    e = user.getExperiment(exp_id)
    tags = []
    try:
        tags = [{"name": t.name} for t in (e.getTags() or [])]
    except Exception:
        pass
    data_fields = []
    try:
        for df in (e.getDataFields() or []):
            data_fields.append({
                "label": getattr(df, "label", ""),
                "field_type": getattr(df, "field_type", "default"),
                "value": str(getattr(df, "value", "") or ""),
            })
    except Exception:
        pass
    linked_protocols = []
    try:
        for p in (e.getProtocols() or []):
            linked_protocols.append({"id": p.id, "name": p.name})
    except Exception:
        pass
    comments = []
    try:
        for c in (e.getComments() or []):
            comments.append({
                "text": getattr(c, "message", str(c)),
                "created_at": str(getattr(c, "created_at", "")),
            })
    except Exception:
        pass
    return [{
        "id": e.id,
        "name": e.name,
        "created_at": str(getattr(e, "created_at", "")),
        "updated_at": str(getattr(e, "updated_at", "")),
        "tags": tags,
        "data_fields": data_fields,
        "protocols": linked_protocols,
        "comments": comments,
    }]


def live_protocols(user, search: str | None, count: int) -> list[dict]:
    kwargs: dict = {"count": count}
    if search:
        kwargs["search_query"] = search
    protos = user.getProtocols(**kwargs)
    results = []
    for p in protos:
        results.append({
            "id": p.id,
            "name": p.name,
            "created_at": str(getattr(p, "created_at", "")),
            "updated_at": str(getattr(p, "updated_at", "")),
            "version": getattr(p, "version", 1),
            "steps": [],
        })
    return results


def live_protocol_detail(user, proto_id: int) -> list[dict]:
    p = user.getProtocol(proto_id)
    steps = []
    try:
        for i, s in enumerate(p.getSteps() or [], 1):
            steps.append({
                "position": i,
                "title": getattr(s, "title", f"Step {i}"),
                "body": getattr(s, "body", ""),
            })
    except Exception:
        pass
    inv_fields = []
    try:
        for f in (p.getInventoryFields() or []):
            inv_fields.append({
                "label": getattr(f, "label", ""),
                "resource_id": getattr(f, "resource_id", 0),
            })
    except Exception:
        pass
    return [{
        "id": p.id,
        "name": p.name,
        "created_at": str(getattr(p, "created_at", "")),
        "updated_at": str(getattr(p, "updated_at", "")),
        "version": getattr(p, "version", 1),
        "steps": steps,
        "inventory_fields": inv_fields,
    }]


def live_inventory(user, search: str | None, count: int) -> dict:
    kwargs: dict = {"count": count}
    if search:
        kwargs["search_query"] = search
    raw_resources = user.getResources(**kwargs)
    resources = []
    for r in raw_resources:
        items = []
        try:
            for item in (r.getItems() or []):
                loc = ""
                try:
                    loc_obj = item.getLocation()
                    loc = getattr(loc_obj, "name", "") if loc_obj else ""
                except Exception:
                    pass
                items.append({
                    "id": item.id,
                    "name": getattr(item, "name", ""),
                    "amount": str(getattr(item, "amount", "") or ""),
                    "location": loc,
                })
        except Exception:
            pass
        resources.append({
            "id": r.id,
            "name": r.name,
            "category": (
                getattr(r, "resource_category", {}).get("name", "")
                if isinstance(getattr(r, "resource_category", None), dict)
                else ""
            ),
            "items": items,
            "metadata": {},
        })

    raw_locs = user.getResourceLocations(count=50)
    locations = [
        {"guid": getattr(loc, "guid", ""), "name": loc.name, "item_count": 0}
        for loc in raw_locs
    ]
    return {"resources": resources, "locations": locations}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio Labstep ELN bridge — query experiments, protocols, inventory"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--demo", action="store_true", help="Run offline demo with synthetic data")
    mode.add_argument("--experiments", action="store_true", help="List experiments")
    mode.add_argument("--experiment-id", type=int, metavar="ID", help="Get full detail for one experiment")
    mode.add_argument("--protocols", action="store_true", help="List protocols")
    mode.add_argument("--protocol-id", type=int, metavar="ID", help="Get full detail for one protocol (with steps)")
    mode.add_argument("--inventory", action="store_true", help="List resources / inventory")

    parser.add_argument("--search", metavar="QUERY", help="Filter by keyword")
    parser.add_argument("--count", type=int, default=20, help="Max items to return (default: 20)")
    parser.add_argument("--output", metavar="DIR", help="Write report.md and result.json to this directory")

    args = parser.parse_args()
    out = Path(args.output) if args.output else None

    # Only non-default flags: the recipe should show what this run chose.
    filters: list[str] = []
    if args.search:
        filters += ["--search", shlex.quote(args.search)]
    if args.count != 20:
        filters += ["--count", str(args.count)]

    if args.demo:
        run_demo(output_dir=out)
        return

    user = get_labstep_user()

    if args.experiments:
        data = live_experiments(user, args.search, args.count)
        title = f"Experiments — \"{args.search}\"" if args.search else "Recent Experiments"
        md = format_experiments(data, title=title)
        print(md)
        if out:
            _write_output(out, md, label="experiments", repro_args=["--experiments", *filters])

    elif args.experiment_id:
        data = live_experiment_detail(user, args.experiment_id)
        md = format_experiments(data, title=f"Experiment #{args.experiment_id}")
        print(md)
        if out:
            _write_output(out, md, label=f"experiment-{args.experiment_id}", repro_args=["--experiment-id", str(args.experiment_id)])

    elif args.protocols:
        data = live_protocols(user, args.search, args.count)
        title = f"Protocols — \"{args.search}\"" if args.search else "Recent Protocols"
        md = format_protocols(data, title=title)
        print(md)
        if out:
            _write_output(out, md, label="protocols", repro_args=["--protocols", *filters])

    elif args.protocol_id:
        data = live_protocol_detail(user, args.protocol_id)
        md = format_protocols(data, title=f"Protocol #{args.protocol_id}")
        print(md)
        if out:
            _write_output(out, md, label=f"protocol-{args.protocol_id}", repro_args=["--protocol-id", str(args.protocol_id)])

    elif args.inventory:
        data = live_inventory(user, args.search, args.count)
        md = format_inventory(data, search=args.search)
        print(md)
        if out:
            _write_output(out, md, label="inventory", repro_args=["--inventory", *filters])


if __name__ == "__main__":
    main()
