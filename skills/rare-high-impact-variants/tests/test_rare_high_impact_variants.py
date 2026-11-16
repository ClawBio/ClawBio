"""Tests for rare-high-impact-variants. Red/green TDD: these should fail until implementation is complete."""

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "rare_high_impact_variants.py"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"


class TestCLI:
    """CLI interface tests."""

    def test_no_args_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            capture_output=True, text=True
        )
        assert result.returncode != 0

    def test_demo_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    def test_input_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(DEMO_INPUT),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()

    @pytest.mark.parametrize("report", ["report.md", "result.json"])
    def test_reports_do_not_record_the_input_directory(self, tmp_path, report):
        """Reports are shared; an absolute input path leaks the user's name and folder layout."""
        private_dir = tmp_path / "private_user_dir"
        private_dir.mkdir()
        vcf = private_dir / "sample.vcf"
        vcf.write_bytes(DEMO_INPUT.read_bytes())
        out = tmp_path / "out"
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--input", str(vcf), "--output", str(out)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        text = (out / report).read_text()
        assert "sample.vcf" in text
        assert str(private_dir) not in text and str(private_dir.resolve()) not in text

    def test_missing_input_exits_nonzero(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(tmp_path / "nonexistent.txt"),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode != 0


class TestOutputFormat:
    """Output format validation."""

    def test_result_json_is_valid(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert isinstance(result, dict)
        assert "skill" in result
        assert result["skill"] == "rare-high-impact-variants"

    def test_report_contains_disclaimer(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        report = (tmp_path / "report.md").read_text()
        assert "not a medical device" in report.lower()

    def test_result_has_variants_count(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert "variants_processed" in result
        assert result["variants_processed"] > 0


class TestDemoData:
    """Demo data integrity."""

    def test_demo_input_exists(self):
        assert DEMO_INPUT.exists(), f"Demo data missing: {DEMO_INPUT}"

    def test_demo_input_has_content(self):
        content = DEMO_INPUT.read_text()
        lines = [l for l in content.splitlines() if l.strip() and not l.startswith("#")]
        assert len(lines) > 0, "Demo input has no data lines"


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree.

    Returns output-relative file paths. Directory lines, and any entry whose
    inline comment contains 'optional', are skipped. Returns [] when there is no
    parseable section, so skills without the section are simply not gated.
    """
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files = []
    parents = {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue  # the root output_directory/ line
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' must be produced.

    Guards against doc/code drift: a skill documenting an output it never writes.
    Mark conditional artifacts with '(optional)' in the tree comment to exempt them.
    """

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing) + ". Write them, mark them "
            "'(optional)' in the SKILL.md tree, or remove them from the "
            "documented Output Structure."
        )


class TestReproducibilityBundle:
    """The demo writes a reproducibility bundle whose checksums verify."""

    def test_demo_bundle_verifies(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        repro = tmp_path / "reproducibility"
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (repro / name).exists(), f"missing reproducibility/{name}"

        commands = (repro / "commands.sh").read_text()
        assert "skills/rare-high-impact-variants/rare_high_impact_variants.py" in commands
        assert "--demo" in commands
        assert str(SKILL_DIR) not in commands

        lines = (repro / "checksums.sha256").read_text().splitlines()
        labels = {line.split("  ", 1)[1] for line in lines}
        assert {"report.md", "result.json", "tables/results.csv",
                "reproducibility/environment.yml"} <= labels
        for line in lines:
            digest, label = line.split("  ", 1)
            path = tmp_path / label
            assert path.is_file(), f"checksum label does not resolve: {label}"
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, label


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/rare-high-impact-variants/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("rare-high-impact-variants", tmp_path) == snapshot({'report.md': """\
# Rare High-Impact Variants Report

**Input**: demo_input.txt
**Rarity threshold**: population AF < 0.01

## 3 rare high-impact variants carried

Of 6 carried, annotated variants, 5 are high-impact (loss-of-function). Of those:

- **3 rare** with documented population frequency below 0.01 (ultra-rare AF < 0.001: 1; rare: 2)
- 1 common (documented AF at or above the threshold)
- 1 with no population-frequency data, so they cannot be confirmed rare (absence of a frequency is not evidence of rarity; many are common LoF polymorphisms)

## Variants

| Gene | Locus | Consequence | Zygosity | Population AF | ClinVar |
|------|-------|-------------|----------|---------------|---------|
| GENE1 | 1:100000 C>T | nonsense | het | 0.0002 | Pathogenic |
| GENE7 | 7:700000 C>T | splice_acceptor | het | 0.002 | - |
| GENE5 | 5:500000 C>G | nonsense | het | 0.004 | - |

## Scope

Counts high-impact (loss-of-function) variants annotated with molecular consequence and population frequency in the input VCF. 'Rare' requires a documented frequency below the threshold; variants with no frequency are reported separately and not called rare. Genome-wide novel LoF calling (VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD) are out of scope for v0.1.0.

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "rare-high-impact-variants",
  "version": "0.1.0",
  "source": "demo_input.txt",
  "rarity_threshold": 0.01,
  "variants_processed": 7,
  "carried_variants": 6,
  "high_impact_carried": 5,
  "rare_high_impact_count": 3,
  "high_impact_common": 1,
  "high_impact_frequency_unknown": 1,
  "by_rarity": {
    "ultra_rare": 1,
    "rare": 2
  },
  "by_consequence": {
    "nonsense": 2,
    "splice_acceptor": 1
  },
  "findings": [
    {
      "gene": "GENE1",
      "rsid": "rs1",
      "locus": "1:100000 C>T",
      "consequence": "nonsense",
      "zygosity": "het",
      "population_af": 0.0002,
      "clnsig": "Pathogenic",
      "rarity": "ultra_rare"
    },
    {
      "gene": "GENE7",
      "rsid": "rs7",
      "locus": "7:700000 C>T",
      "consequence": "splice_acceptor",
      "zygosity": "het",
      "population_af": 0.002,
      "clnsig": "",
      "rarity": "rare"
    },
    {
      "gene": "GENE5",
      "rsid": "rs5",
      "locus": "5:500000 C>G",
      "consequence": "nonsense",
      "zygosity": "het",
      "population_af": 0.004,
      "clnsig": "",
      "rarity": "rare"
    }
  ],
  "frequency_unknown_genes": [
    "GENE2"
  ],
  "scope_note": "Counts high-impact (loss-of-function) variants annotated with molecular consequence and population frequency in the input VCF. 'Rare' requires a documented frequency below the threshold; variants with no frequency are reported separately and not called rare. Genome-wide novel LoF calling (VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD) are out of scope for v0.1.0.",
  "disclaimer": "ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."
}\
""", 'tables/results.csv': """\
gene,locus,consequence,zygosity,population_af,rarity,clnsig
GENE1,1:100000 C>T,nonsense,het,0.0002,ultra_rare,Pathogenic
GENE7,7:700000 C>T,splice_acceptor,het,0.002,rare,
GENE5,5:500000 C>G,nonsense,het,0.004,rare,
"""})
