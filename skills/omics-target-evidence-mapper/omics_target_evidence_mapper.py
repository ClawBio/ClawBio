#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import requests

from clawbio.common.reproducibility import write_checksums, write_environment_yml, write_ro_crate


DEMO_GENE = "IL6R"
DEMO_DISEASE = "coronary artery disease"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate public target-level evidence across omics and translational sources."
    )
    parser.add_argument("--gene", type=str, help="Gene or target symbol")
    parser.add_argument("--disease", type=str, help="Optional disease term")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument("--max-papers", type=int, default=5, help="Maximum number of PubMed hits")
    parser.add_argument("--max-trials", type=int, default=5, help="Maximum number of trial hits")
    parser.add_argument("--demo", action="store_true", help="Run the built-in demo query")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.demo and not args.gene:
        raise ValueError("Provide --gene or use --demo.")
    if args.max_papers < 1 or args.max_papers > 20:
        raise ValueError("--max-papers must be between 1 and 20.")
    if args.max_trials < 1 or args.max_trials > 20:
        raise ValueError("--max-trials must be between 1 and 20.")


def safe_request_json(
    method: str,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    json_body: dict[str, Any] | None = None,
    timeout: int = 30,
) -> dict[str, Any] | None:
    try:
        response = requests.request(method, url, params=params, json=json_body, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as exc:
        print(f"[WARN] {method} {url} failed: {exc}", file=sys.stderr)
        return None


def fetch_uniprot_summary(gene: str) -> dict[str, Any]:
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene} AND reviewed:true AND organism_id:9606",
        "format": "json",
        "size": 1,
    }
    data = safe_request_json("GET", url, params=params)

    if data is None:
        return {"status": "unavailable", "gene": gene}
    if not data.get("results"):
        return {"status": "no_result", "gene": gene}

    entry = data["results"][0]
    protein_desc = (
        entry.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
    )

    organism = entry.get("organism", {}).get("scientificName")
    primary_accession = entry.get("primaryAccession")
    uni_name = entry.get("uniProtkbId")

    return {
        "status": "ok",
        "gene": gene,
        "primary_accession": primary_accession,
        "uniprot_id": uni_name,
        "protein_name": protein_desc,
        "organism": organism,
    }


def fetch_pubmed_hits(gene: str, disease: str | None, max_papers: int) -> dict[str, Any]:
    term = gene if not disease else f"{gene} AND {disease}"
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": max_papers,
        "sort": "pub date",
    }
    search_data = safe_request_json("GET", search_url, params=search_params)

    if search_data is None:
        return {"status": "unavailable", "items": []}

    pmids = search_data.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return {"status": "no_result", "items": []}

    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "json",
    }
    summary_data = safe_request_json("GET", summary_url, params=summary_params)
    if summary_data is None:
        return {"status": "unavailable", "items": []}

    results = []
    for pmid in pmids:
        item = summary_data.get("result", {}).get(pmid, {})
        if not item:
            continue
        results.append(
            {
                "pmid": pmid,
                "title": item.get("title"),
                "pubdate": item.get("pubdate"),
                "source": item.get("source"),
            }
        )
    return {"status": "ok" if results else "no_result", "items": results}


def fetch_open_targets_evidence(gene: str, disease: str | None) -> dict[str, Any]:
    if not disease:
        return {"status": "skipped", "reason": "No disease term provided."}

    url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Step 1: resolve IDs via search
    search_query = """
    query Search($geneText: String!, $diseaseText: String!) {
      geneResult: search(queryString: $geneText, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits { id }
      }
      diseaseResult: search(queryString: $diseaseText, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits { id }
      }
    }
    """
    search_data = safe_request_json("POST", url, json_body={"query": search_query, "variables": {"geneText": gene, "diseaseText": disease}})

    if not search_data or "data" not in search_data:
        return {"status": "unavailable", "gene": gene, "disease": disease}

    target_id = ((search_data["data"].get("geneResult") or {}).get("hits") or [{}])[0].get("id")
    disease_id = ((search_data["data"].get("diseaseResult") or {}).get("hits") or [{}])[0].get("id")

    if not target_id and not disease_id:
        return {"status": "no_result", "gene_query": gene, "disease_query": disease, "matched_target": None, "matched_disease": None}

    # Step 2: enrich with structured target/disease data
    enrich_query = """
    query Enrich($ensemblId: String!, $efoId: String!) {
      target(ensemblId: $ensemblId) { id approvedSymbol approvedName biotype }
      disease(efoId: $efoId) { id name description }
    }
    """
    enrich_data = safe_request_json("POST", url, json_body={"query": enrich_query, "variables": {"ensemblId": target_id or "", "efoId": disease_id or ""}})

    enrich = (enrich_data or {}).get("data", {})
    target_hit = enrich.get("target")
    disease_hit = enrich.get("disease")

    return {
        "status": "ok" if target_hit or disease_hit else "no_result",
        "gene_query": gene,
        "disease_query": disease,
        "matched_target": target_hit,
        "matched_disease": disease_hit,
    }


def fetch_trials(gene: str, disease: str | None, max_trials: int) -> dict[str, Any]:
    query = gene if not disease else f"{gene} {disease}"
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": query,
        "pageSize": max_trials,
        "format": "json",
    }
    data = safe_request_json("GET", url, params=params)

    if data is None:
        return {"status": "unavailable", "items": []}

    results = []
    for study in data.get("studies", []):
        proto = study.get("protocolSection", {})
        id_mod = proto.get("identificationModule", {})
        status_mod = proto.get("statusModule", {})
        design_mod = proto.get("designModule", {})
        results.append(
            {
                "nct_id": id_mod.get("nctId"),
                "title": id_mod.get("briefTitle"),
                "status": status_mod.get("overallStatus"),
                "phase": (design_mod.get("phases") or [None])[0],
            }
        )
    return {"status": "ok" if results else "no_result", "items": results}


def build_evidence(args: argparse.Namespace) -> dict[str, Any]:
    gene = DEMO_GENE if args.demo else args.gene
    disease = DEMO_DISEASE if args.demo else args.disease

    target_summary = fetch_uniprot_summary(gene)
    disease_association = fetch_open_targets_evidence(gene, disease)
    literature = fetch_pubmed_hits(gene, disease, args.max_papers)
    trials = fetch_trials(gene, disease, args.max_trials)

    limitations = [
        "This tool aggregates public evidence for research triage only.",
        "It does not infer causality from association evidence.",
        "It does not provide clinical recommendations.",
        "Public API availability may affect completeness.",
    ]
    if target_summary.get("status") == "unavailable":
        limitations.append("UniProt was unavailable; target summary was not assessed.")
    if disease_association.get("status") == "unavailable":
        limitations.append("Open Targets was unavailable; disease association was not assessed.")
    if literature.get("status") == "unavailable":
        limitations.append("PubMed was unavailable; literature was not assessed.")
    if trials.get("status") == "unavailable":
        limitations.append("ClinicalTrials.gov was unavailable; trials were not assessed.")

    return {
        "query": {
            "gene": gene,
            "disease": disease,
            "demo_mode": args.demo,
        },
        "target_summary": target_summary,
        "disease_association": disease_association,
        "literature": literature,
        "trials": trials,
        "limitations": limitations,
        "provenance": {
            "sources": ["UniProt", "Open Targets", "PubMed", "ClinicalTrials.gov"],
            "generated_at_utc": datetime.now(UTC).isoformat(),
            "version": "0.1.0",
        },
    }


def build_report(evidence: dict[str, Any]) -> str:
    query = evidence["query"]
    target = evidence["target_summary"]
    disease_assoc = evidence["disease_association"]
    literature = evidence["literature"]
    if not isinstance(literature, dict):
        literature = {
            "status": "ok" if literature else "no_result",
            "items": literature or [],
        }
    trials = evidence["trials"]
    if not isinstance(trials, dict):
        trials = {
            "status": "ok" if trials else "no_result",
            "items": trials or [],
        }
    limitations = evidence["limitations"]
    literature_items = literature.get("items") or []
    trial_items = trials.get("items") or []

    lines = [
        "# Omics-to-Target Evidence Mapper Report",
        "",
        "## Query",
        f"- Gene: {query.get('gene')}",
        f"- Disease: {query.get('disease')}",
        f"- Demo mode: {query.get('demo_mode')}",
        "",
        "## Target Summary",
        f"- Status: {target.get('status')}",
        f"- UniProt accession: {target.get('primary_accession')}",
        f"- UniProt ID: {target.get('uniprot_id')}",
        f"- Protein name: {target.get('protein_name')}",
        f"- Organism: {target.get('organism')}",
        "",
        "## Disease Association",
        f"- Status: {disease_assoc.get('status')}",
        f"- Gene query: {disease_assoc.get('gene_query')}",
        f"- Disease query: {disease_assoc.get('disease_query')}",
        f"- Matched target: {json.dumps(disease_assoc.get('matched_target'), ensure_ascii=False)}",
        f"- Matched disease: {json.dumps(disease_assoc.get('matched_disease'), ensure_ascii=False)}",
        "",
        "## Literature Snapshot",
    ]

    if literature.get("status") == "unavailable":
        lines.append("- Literature unavailable; not assessed.")
    elif literature_items:
        for paper in literature_items:
            lines.extend(
                [
                    f"- PMID: {paper.get('pmid')}",
                    f"  - Title: {paper.get('title')}",
                    f"  - Date: {paper.get('pubdate')}",
                    f"  - Source: {paper.get('source')}",
                ]
            )
    else:
        lines.append("- No literature hits found.")

    lines.extend(["", "## Trial Landscape"])
    if trials.get("status") == "unavailable":
        lines.append("- Trials unavailable; not assessed.")
    elif trial_items:
        for trial in trial_items:
            lines.extend(
                [
                    f"- NCT ID: {trial.get('nct_id')}",
                    f"  - Title: {trial.get('title')}",
                    f"  - Status: {trial.get('status')}",
                    f"  - Phase: {trial.get('phase')}",
                ]
            )
    else:
        lines.append("- No trial hits found.")

    lines.extend(["", "## Limitations and Safety Notes"])
    for item in limitations:
        lines.append(f"- {item}")

    lines.extend(
        [
            "",
            "## Provenance",
            "- Sources: UniProt, Open Targets, PubMed, ClinicalTrials.gov",
        ]
    )

    return "\n".join(lines) + "\n"


def write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def write_json(path: Path, content: dict[str, Any]) -> None:
    path.write_text(json.dumps(content, indent=2, ensure_ascii=False), encoding="utf-8")


def _source_items(block: Any) -> list[Any]:
    if isinstance(block, dict):
        items = block.get("items")
        return items if isinstance(items, list) else []
    if isinstance(block, list):
        return block
    return []


def _source_status(block: Any) -> str:
    if isinstance(block, dict) and isinstance(block.get("status"), str):
        return block["status"]
    return "ok" if block else "no_result"


def main() -> None:
    args = parse_args()
    validate_args(args)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    evidence = build_evidence(args)
    report = build_report(evidence)

    evidence_path = output_dir / "evidence.json"
    report_path = output_dir / "report.md"
    metadata_path = output_dir / "metadata.json"

    metadata = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "query": evidence["query"],
        "sources": evidence["provenance"]["sources"],
        "counts": {
            "literature": len(_source_items(evidence["literature"])),
            "trials": len(_source_items(evidence["trials"])),
        },
        "status": {
            "literature": _source_status(evidence["literature"]),
            "trials": _source_status(evidence["trials"]),
        },
    }

    write_json(evidence_path, evidence)
    write_text(report_path, report)
    write_json(metadata_path, metadata)
    environment_path = write_environment_yml(
        output_dir,
        env_name="clawbio-omics-target-evidence-mapper",
        pip_deps=["requests", "rocrate"],
        python_version="3.11",
    )
    ro_crate_path = write_ro_crate(
        output_dir,
        skill_name="omics-target-evidence-mapper",
        skill_version="0.1.0",
        script_path="skills/omics-target-evidence-mapper/omics_target_evidence_mapper.py",
        description="Aggregate public target-level evidence across omics and translational sources",
        params={
            "gene": evidence["query"]["gene"],
            "disease": evidence["query"].get("disease", ""),
            "max_papers": args.max_papers,
            "max_trials": args.max_trials,
            "demo": args.demo,
        },
    )
    write_checksums(
        [evidence_path, report_path, metadata_path, environment_path, ro_crate_path],
        output_dir,
        anchor=output_dir,
    )

    print(f"Done. Output written to: {output_dir}")


if __name__ == "__main__":
    main()
