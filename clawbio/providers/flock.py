"""FLock API provider for ClawBio intelligent routing.

FLock.io provides open-source model inference through an OpenAI-compatible API.
This module uses FLock as the LLM backend for the Bio Orchestrator's
intelligent routing when keyword matching fails.

Usage:
    from clawbio.providers.flock import FlockRouter
    router = FlockRouter()
    skill = router.route_query("What drugs should I avoid with my genotype?")
"""

from __future__ import annotations

import json
import os

try:
    from openai import OpenAI
except ImportError:
    OpenAI = None  # type: ignore[assignment,misc]

FLOCK_BASE_URL = "https://api.flock.io/v1"
FLOCK_DEFAULT_MODEL = "gemini-3-flash-preview"

# All available skills with descriptions for the LLM to choose from
SKILL_DESCRIPTIONS = {
    "pharmgx-reporter": "Pharmacogenomics report from consumer genetic data (23andMe/AncestryDNA). Drug-gene interactions, star alleles, CPIC guidelines, metaboliser phenotypes.",
    "drug-photo": "Identify a medication from a photo and return a personalised dosage card against the user's genotype.",
    "clinpgx": "Gene-drug lookup from ClinPGx, PharmGKB, CPIC, and FDA drug label databases.",
    "gwas-lookup": "Federated variant query across 9 genomic databases (GWAS Catalog, gnomAD, ClinVar, GTEx, etc.).",
    "gwas-prs": "Polygenic risk score calculation from consumer genetic data using PGS Catalog.",
    "profile-report": "Unified personal genomic profile report aggregating results from multiple skills.",
    "equity-scorer": "HEIM diversity metrics (FST, heterozygosity, representation) from VCF or ancestry data.",
    "nutrigx": "Personalised nutrigenomics recommendations from genetic data — diet, vitamins, caffeine, lactose.",
    "claw-ancestry-pca": "Ancestry decomposition PCA against the Simons Genome Diversity Project.",
    "claw-semantic-sim": "Semantic Isolation Index for disease research equity using PubMed abstracts.",
    "claw-metagenomics": "Shotgun metagenomics profiling — taxonomy, resistome, and functional pathways.",
    "genome-compare": "Pairwise IBS comparison vs George Church (PGP-1) plus ancestry estimation.",
    "ukb-navigator": "Semantic search across 22,000+ UK Biobank fields.",
    "scrna-orchestrator": "Scanpy single-cell RNA-seq pipeline — QC, clustering, marker DE analysis.",
    "bio-orchestrator": "Meta-agent that routes bioinformatics requests to the right specialist skill.",
}

ROUTING_SYSTEM_PROMPT = """You are the ClawBio Bio Orchestrator routing agent. Given a user's bioinformatics query, determine which skill should handle it.

Available skills:
{skills}

Respond with ONLY a JSON object: {{"skill": "<skill-name>", "confidence": <0.0-1.0>, "reasoning": "<one sentence>"}}

If no skill matches, respond: {{"skill": null, "confidence": 0.0, "reasoning": "<why>"}}"""


class FlockRouter:
    """Routes bioinformatics queries to skills using FLock open-source models."""

    def __init__(
        self,
        api_key: str | None = None,
        model: str = FLOCK_DEFAULT_MODEL,
        base_url: str = FLOCK_BASE_URL,
    ):
        if OpenAI is None:
            raise ImportError(
                "openai package required for FLock provider. "
                "Install with: pip install openai>=1.0"
            )
        self.api_key = api_key or os.environ.get("FLOCK_API_KEY", "")
        if not self.api_key:
            raise ValueError(
                "FLOCK_API_KEY not set. Get one at https://platform.flock.io"
            )
        self.model = model
        self.client = OpenAI(
            base_url=base_url,
            api_key=self.api_key,
            default_headers={"x-litellm-api-key": self.api_key},
        )

    def route_query(self, query: str) -> dict:
        """Route a natural language query to a ClawBio skill.

        Returns:
            dict with keys: skill (str|None), confidence (float), reasoning (str)
        """
        skills_text = "\n".join(
            f"- **{name}**: {desc}" for name, desc in SKILL_DESCRIPTIONS.items()
        )
        system_prompt = ROUTING_SYSTEM_PROMPT.format(skills=skills_text)

        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": query},
            ],
            temperature=0.1,
            max_tokens=256,
        )

        content = response.choices[0].message.content.strip()

        # Parse JSON response — handle markdown code fences
        if content.startswith("```"):
            content = content.split("\n", 1)[1].rsplit("```", 1)[0].strip()

        try:
            result = json.loads(content)
        except json.JSONDecodeError:
            result = {"skill": None, "confidence": 0.0, "reasoning": f"Failed to parse LLM response: {content[:100]}"}

        return result

    def route_query_safe(self, query: str) -> dict:
        """Route with error handling — never raises, returns null skill on failure."""
        try:
            return self.route_query(query)
        except Exception as e:
            return {"skill": None, "confidence": 0.0, "reasoning": f"FLock API error: {e}"}
