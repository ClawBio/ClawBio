"""
Tests for phylogenetics_builder.py — Tech Lead quality TDD suite.

Red/Green TDD cycle:
  1. Run:  pytest skills/phylogenetics-builder/tests/ -v  → expect failures (Red)
  2. Implement phylogenetics_builder.py to make all tests pass (Green)
  3. Refactor with tests always passing

Coverage areas:
  • FASTA parsing (aligned + unaligned modes)
  • MSA command generation (all 6 aligners)
  • trimAl command generation
  • ModelFinder output parsing
  • IQ-TREE command generation (ufboot / standard / all bootstrap modes)
  • RAxML-NG command generation (check + run)
  • Midpoint rooting (ETE3, with graceful fallback)
  • Demo mode end-to-end (no external binaries)
  • CLI argument validation
  • Output file contract
  • result.json schema
  • Disclaimer in report

Markers:
  integration — requires external binaries (iqtree2, raxml-ng, mafft, trimal)
                skip with:  pytest -m "not integration" -v
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "phylogenetics_builder.py"
DEMO_INPUT = SKILL_DIR / "demo_alignment.fasta"
EXAMPLES_DIR = SKILL_DIR / "examples"

VALID_ALIGNERS = ("mafft", "muscle", "clustalw", "kalign", "tcoffee", "prank")
VALID_ENGINES = ("iqtree2", "raxml-ng")
VALID_BOOTSTRAP = ("ufboot", "standard", "all")


# ── Module loader ──────────────────────────────────────────────────────────────


def get_module():
    if not SCRIPT.exists():
        raise ImportError(f"Script {SCRIPT} does not exist yet.")
    spec = spec_from_file_location("phylogenetics_builder", SCRIPT)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_cli(args, **kwargs):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
        **kwargs,
    )


# ── FASTA parsing ──────────────────────────────────────────────────────────────


class TestParseFASTA:
    def test_valid_aligned_fasta(self, tmp_path):
        fasta = ">s1\nACGT\n>s2\nACGT\n>s3\nACGT\n"
        f = tmp_path / "a.fasta"
        f.write_text(fasta)
        m = get_module()
        seqs = m.parse_fasta(f, require_aligned=True)
        assert len(seqs) == 3
        assert seqs["s1"] == "ACGT"

    def test_unaligned_fasta_passes_without_length_check(self, tmp_path):
        fasta = ">s1\nACGTT\n>s2\nACGT\n>s3\nACGTACGT\n"
        f = tmp_path / "u.fasta"
        f.write_text(fasta)
        m = get_module()
        seqs = m.parse_fasta(f, require_aligned=False)
        assert len(seqs) == 3

    def test_aligned_rejects_unequal_lengths(self, tmp_path):
        fasta = ">s1\nACGT\n>s2\nACGTT\n>s3\nACGT\n"
        f = tmp_path / "u.fasta"
        f.write_text(fasta)
        m = get_module()
        with pytest.raises(ValueError, match="same length"):
            m.parse_fasta(f, require_aligned=True)

    def test_too_few_sequences(self, tmp_path):
        fasta = ">s1\nACGT\n>s2\nACGT\n"
        f = tmp_path / "few.fasta"
        f.write_text(fasta)
        m = get_module()
        with pytest.raises(ValueError, match="at least 3"):
            m.parse_fasta(f)

    def test_empty_header_raises(self, tmp_path):
        fasta = ">s1\nACGT\n>\nACGT\n>s3\nACGT\n"
        f = tmp_path / "eh.fasta"
        f.write_text(fasta)
        m = get_module()
        with pytest.raises(ValueError, match="Empty sequence header"):
            m.parse_fasta(f)

    def test_duplicate_header_raises(self, tmp_path):
        fasta = ">s1\nACGT\n>s1\nACGG\n>s3\nACGT\n"
        f = tmp_path / "dup.fasta"
        f.write_text(fasta)
        m = get_module()
        with pytest.raises(ValueError, match="Duplicate sequence header"):
            m.parse_fasta(f)

    def test_multiline_sequences_concatenated(self, tmp_path):
        fasta = ">s1\nAC\nGT\n>s2\nACGT\n>s3\nACGT\n"
        f = tmp_path / "ml.fasta"
        f.write_text(fasta)
        m = get_module()
        seqs = m.parse_fasta(f)
        assert seqs["s1"] == "ACGT"

    def test_file_not_found_raises(self, tmp_path):
        m = get_module()
        with pytest.raises(FileNotFoundError):
            m.parse_fasta(tmp_path / "ghost.fasta")


# ── MSA command generation ─────────────────────────────────────────────────────


class TestMSACommands:
    """build_msa_command() returns a dict with keys: cmd (list), stdin_from, stdout_to."""

    def _cmd(self, aligner, tmp_path):
        m = get_module()
        return m.build_msa_command(aligner, tmp_path / "in.fa", tmp_path / "out.fa")

    def test_mafft_uses_auto_flag(self, tmp_path):
        r = self._cmd("mafft", tmp_path)
        assert "mafft" in r["cmd"][0]
        assert "--auto" in r["cmd"]
        assert r.get("stdout_to") is not None  # mafft writes to stdout

    def test_muscle_uses_output_flag(self, tmp_path):
        r = self._cmd("muscle", tmp_path)
        assert "muscle" in r["cmd"][0]
        assert "-align" in r["cmd"]
        assert "-output" in r["cmd"]

    def test_clustalw_uses_outfile_flag(self, tmp_path):
        r = self._cmd("clustalw", tmp_path)
        assert "clustalw" in r["cmd"][0] or "clustal" in r["cmd"][0]
        assert any(
            "-OUTFILE" in x or "-outfile" in x or str(tmp_path / "out.fa") in x
            for x in r["cmd"]
        )

    def test_kalign_uses_stdin_stdout(self, tmp_path):
        r = self._cmd("kalign", tmp_path)
        assert "kalign" in r["cmd"][0]
        assert r.get("stdin_from") is not None
        assert r.get("stdout_to") is not None

    def test_tcoffee_uses_infile_outfile(self, tmp_path):
        r = self._cmd("tcoffee", tmp_path)
        assert "t_coffee" in r["cmd"][0] or "tcoffee" in r["cmd"][0]
        cmd_str = " ".join(r["cmd"])
        assert "infile" in cmd_str or str(tmp_path / "in.fa") in cmd_str

    def test_prank_uses_d_flag(self, tmp_path):
        r = self._cmd("prank", tmp_path)
        assert "prank" in r["cmd"][0]
        cmd_str = " ".join(r["cmd"])
        assert "-d=" in cmd_str or "-d" in cmd_str

    def test_invalid_aligner_raises(self, tmp_path):
        m = get_module()
        with pytest.raises(ValueError, match="Unknown aligner"):
            m.build_msa_command("fakealigner", tmp_path / "in.fa", tmp_path / "out.fa")


# ── trimAl command generation ──────────────────────────────────────────────────


class TestTrimAlCommand:
    def test_default_strategy_is_automated1(self, tmp_path):
        m = get_module()
        cmd = m.build_trimal_command(tmp_path / "in.fa", tmp_path / "out.fa")
        assert "trimal" in cmd[0]
        assert "-automated1" in cmd

    def test_nogaps_strategy(self, tmp_path):
        m = get_module()
        cmd = m.build_trimal_command(
            tmp_path / "in.fa", tmp_path / "out.fa", strategy="-nogaps"
        )
        assert "-nogaps" in cmd
        assert "-automated1" not in cmd

    def test_input_output_flags(self, tmp_path):
        m = get_module()
        cmd = m.build_trimal_command(tmp_path / "in.fa", tmp_path / "out.fa")
        assert "-in" in cmd
        assert "-out" in cmd


# ── ModelFinder parsing ────────────────────────────────────────────────────────


class TestModelFinderParsing:
    def test_parses_bic_model(self, tmp_path):
        iqtree_content = """
...
Best-fit model according to BIC: TIM3+F+G4

List of models sorted by BIC scores:
...
"""
        f = tmp_path / "run.iqtree"
        f.write_text(iqtree_content)
        m = get_module()
        assert m.parse_iqtree_model(f) == "TIM3+F+G4"

    def test_parses_tvme_model(self, tmp_path):
        iqtree_content = "Best-fit model according to BIC: TVMe+I+R2\n"
        f = tmp_path / "run.iqtree"
        f.write_text(iqtree_content)
        m = get_module()
        assert m.parse_iqtree_model(f) == "TVMe+I+R2"

    def test_fallback_when_no_model_line(self, tmp_path):
        f = tmp_path / "run.iqtree"
        f.write_text("No model info here.\n")
        m = get_module()
        model = m.parse_iqtree_model(f)
        assert "GTR" in model  # fallback is GTR-based

    def test_handles_legacy_format(self, tmp_path):
        iqtree_content = "Best-fit model: HKY+G\n"
        f = tmp_path / "run.iqtree"
        f.write_text(iqtree_content)
        m = get_module()
        assert m.parse_iqtree_model(f) == "HKY+G"


# ── IQ-TREE command generation ─────────────────────────────────────────────────


class TestIQTREECommands:
    def _make_cmd(self, bootstrap="ufboot", outgroup=None, threads=2, seed=42):
        m = get_module()
        return m.build_iqtree_command(
            "iqtree2",
            Path("/tmp/aln.fa"),
            "GTR+F+G4",
            Path("/tmp/prefix"),
            bootstrap,
            outgroup,
            threads,
            seed,
        )

    def test_ufboot_uses_bb_flag(self):
        cmd = self._make_cmd("ufboot")
        assert "-bb" in cmd
        idx = cmd.index("-bb")
        assert cmd[idx + 1] == "1000"
        assert "-b" not in cmd
        assert "-alrt" not in cmd

    def test_standard_bootstrap_uses_b_flag(self):
        cmd = self._make_cmd("standard")
        assert "-b" in cmd
        idx = cmd.index("-b")
        assert cmd[idx + 1] == "100"
        assert "-bb" not in cmd

    def test_all_bootstrap_has_triple_flags(self):
        cmd = self._make_cmd("all")
        assert "-bb" in cmd
        assert "-alrt" in cmd
        assert "-abayes" in cmd

    def test_outgroup_is_appended(self):
        cmd = self._make_cmd(outgroup="Mus_musculus")
        assert "-o" in cmd
        idx = cmd.index("-o")
        assert cmd[idx + 1] == "Mus_musculus"

    def test_multiple_outgroups_comma_separated(self):
        cmd = self._make_cmd(outgroup="Mus_musculus,Rattus_norvegicus")
        assert "-o" in cmd
        idx = cmd.index("-o")
        assert "Mus_musculus" in cmd[idx + 1]

    def test_no_outgroup_by_default(self):
        cmd = self._make_cmd(outgroup=None)
        assert "-o" not in cmd

    def test_input_file_in_command(self):
        cmd = self._make_cmd()
        cmd_str = " ".join(cmd)
        assert "/tmp/aln.fa" in cmd_str

    def test_model_in_command(self):
        cmd = self._make_cmd()
        assert "-m" in cmd
        idx = cmd.index("-m")
        assert cmd[idx + 1] == "GTR+F+G4"

    def test_prefix_in_command(self):
        cmd = self._make_cmd()
        assert "--prefix" in cmd


# ── RAxML-NG command generation ────────────────────────────────────────────────


class TestRAxMLCommands:
    def test_check_command_structure(self, tmp_path):
        m = get_module()
        cmd = m.build_raxml_check_command(
            "raxml-ng", tmp_path / "aln.fa", "TIM3+G4", tmp_path / "prefix"
        )
        assert "raxml-ng" in cmd[0]
        assert "--check" in cmd
        assert "--msa" in cmd
        assert "--model" in cmd

    def test_run_command_has_msa_model_prefix(self, tmp_path):
        m = get_module()
        cmd = m.build_raxml_run_command(
            "raxml-ng",
            tmp_path / "aln.fa",
            "TIM3+G4",
            tmp_path / "prefix",
            "standard",
            None,
            2,
            42,
        )
        assert "--msa" in cmd
        assert "--model" in cmd
        assert "--prefix" in cmd
        assert "--threads" in cmd
        assert "--seed" in cmd

    def test_run_command_bootstrap_trees(self, tmp_path):
        m = get_module()
        cmd = m.build_raxml_run_command(
            "raxml-ng",
            tmp_path / "aln.fa",
            "TIM3+G4",
            tmp_path / "prefix",
            "standard",
            None,
            2,
            42,
        )
        assert "--bs-trees" in cmd

    def test_outgroup_added_when_specified(self, tmp_path):
        m = get_module()
        cmd = m.build_raxml_run_command(
            "raxml-ng",
            tmp_path / "aln.fa",
            "TIM3+G4",
            tmp_path / "prefix",
            "ufboot",
            "Mus_musculus",
            2,
            42,
        )
        assert "--outgroup" in cmd
        idx = cmd.index("--outgroup")
        assert "Mus_musculus" in cmd[idx + 1]

    def test_check_not_in_run_command(self, tmp_path):
        m = get_module()
        cmd = m.build_raxml_run_command(
            "raxml-ng",
            tmp_path / "aln.fa",
            "GTR+G4",
            tmp_path / "prefix",
            "ufboot",
            None,
            2,
            42,
        )
        assert "--check" not in cmd


# ── Model adaptation ───────────────────────────────────────────────────────────


class TestModelAdaptation:
    def test_strips_plus_f_for_raxml(self):
        m = get_module()
        adapted = m.adapt_model_for_engine("TIM3+F+G4", "raxml-ng")
        assert "+F" not in adapted
        assert "TIM3" in adapted
        assert "G4" in adapted

    def test_gtr_f_g4_stripped(self):
        m = get_module()
        adapted = m.adapt_model_for_engine("GTR+F+G4", "raxml-ng")
        assert "+F" not in adapted

    def test_iqtree2_model_unchanged(self):
        m = get_module()
        adapted = m.adapt_model_for_engine("TIM3+F+G4", "iqtree2")
        assert adapted == "TIM3+F+G4"

    def test_model_without_f_unchanged(self):
        m = get_module()
        adapted = m.adapt_model_for_engine("HKY+G4", "raxml-ng")
        assert adapted == "HKY+G4"


# ── Newick table parsing ───────────────────────────────────────────────────────


class TestNewickToTable:
    def test_parses_support_values(self):
        newick = "(A:0.02,B:0.03)100:0.05;"
        m = get_module()
        table = m.newick_to_table(newick)
        assert len(table) > 0
        assert any(row["support"] == 100 for row in table)

    def test_branch_length_parsed(self):
        newick = "(A:0.02,B:0.03)100:0.05;"
        m = get_module()
        table = m.newick_to_table(newick)
        lengths = [row["length"] for row in table]
        assert any(abs(l - 0.02) < 0.001 for l in lengths)

    def test_handles_triple_support_labels(self):
        newick = "(A:0.02,B:0.03)80/0.8/95:0.05;"
        m = get_module()
        table = m.newick_to_table(newick)
        assert len(table) > 0


# ── Midpoint rooting ───────────────────────────────────────────────────────────


class TestMidpointRooting:
    def test_returns_newick_string(self):
        newick = "(A:0.02,B:0.03,C:0.05);"
        m = get_module()
        result = m.root_midpoint(newick)
        assert isinstance(result, str)
        assert len(result) > 5

    def test_graceful_fallback_without_ete3(self):
        newick = "(A:0.02,B:0.03,C:0.05);"
        m = get_module()
        with patch.dict("sys.modules", {"ete3": None}):
            result = m.root_midpoint(newick)
        assert isinstance(result, str)
        assert len(result) > 0


# ── CLI argument validation ────────────────────────────────────────────────────


class TestCLIArgs:
    def test_demo_mode_requires_output_or_uses_default(self, tmp_path):
        result = run_cli(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, f"stderr: {result.stderr}"

    def test_missing_input_nonzero_exit(self, tmp_path):
        result = run_cli(
            ["--input", str(tmp_path / "nope.fa"), "--output", str(tmp_path)]
        )
        assert result.returncode != 0

    def test_no_args_nonzero_exit(self):
        result = run_cli([])
        assert result.returncode != 0

    def test_invalid_aligner_rejected(self, tmp_path):
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--aligner",
                "bogus_aligner",
                "--aligned",
            ]
        )
        assert result.returncode != 0

    def test_invalid_engine_rejected(self, tmp_path):
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--engine",
                "bogus_engine",
                "--aligned",
            ]
        )
        assert result.returncode != 0

    def test_invalid_bootstrap_rejected(self, tmp_path):
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--bootstrap",
                "bogus_bootstrap",
                "--aligned",
            ]
        )
        assert result.returncode != 0

    def test_aligned_flag_accepted(self, tmp_path):
        result = run_cli(["--demo", "--output", str(tmp_path), "--aligned"])
        assert result.returncode == 0, f"stderr: {result.stderr}"

    def test_no_trim_flag_accepted(self, tmp_path):
        result = run_cli(["--demo", "--output", str(tmp_path), "--no-trim"])
        assert result.returncode == 0, f"stderr: {result.stderr}"


# ── Demo mode end-to-end ───────────────────────────────────────────────────────


class TestDemoMode:
    def test_demo_creates_required_outputs(self, tmp_path):
        result = run_cli(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "phylo_tree.nwk").exists()
        assert (tmp_path / "figures" / "phylogram.png").exists()
        assert (tmp_path / "tables" / "branch_support.csv").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "environment.yml").exists()

    def test_demo_no_binaries_in_path(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PATH", "")
        result = run_cli(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    def test_demo_newick_is_valid(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        nwk = (tmp_path / "phylo_tree.nwk").read_text().strip()
        assert nwk.startswith("(")
        assert nwk.endswith(";")

    def test_demo_taxa_count_matches_fasta(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        result_data = json.loads((tmp_path / "result.json").read_text())
        m = get_module()
        seqs = m.parse_fasta(DEMO_INPUT)
        assert result_data["num_taxa"] == len(seqs)


# ── result.json schema ─────────────────────────────────────────────────────────


class TestResultJSON:
    def test_required_fields_present(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        required = [
            "engine",
            "aligner",
            "aligned",
            "trimmed",
            "model",
            "bootstrap_mode",
            "num_taxa",
            "tree",
            "input_file",
            "output_dir",
            "pipeline_steps",
            "status",
            "chat_summary_lines",
            "preferred_artifacts",
            "workflow_state",
            "suggested_actions",
        ]
        for field in required:
            assert field in data, f"Missing field: {field}"

    def test_status_is_success(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["status"] == "success"

    def test_engine_is_known_value(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["engine"] in ("iqtree2", "raxml-ng", "precomputed")

    def test_workflow_state_is_completed(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["workflow_state"] == "completed"

    def test_pipeline_steps_is_list(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        assert isinstance(data["pipeline_steps"], list)
        assert len(data["pipeline_steps"]) >= 1

    def test_preferred_artifacts_list_with_paths(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        data = json.loads((tmp_path / "result.json").read_text())
        artifacts = data["preferred_artifacts"]
        assert isinstance(artifacts, list)
        assert len(artifacts) >= 1
        for a in artifacts:
            assert "path" in a


# ── Report content ─────────────────────────────────────────────────────────────


class TestReportContent:
    def test_disclaimer_in_report(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        content = (tmp_path / "report.md").read_text().lower()
        assert "not a medical device" in content

    def test_model_name_in_report(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        content = (tmp_path / "report.md").read_text()
        assert "Model" in content or "model" in content

    def test_pipeline_section_in_report(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        content = (tmp_path / "report.md").read_text()
        assert "Pipeline" in content or "pipeline" in content or "Step" in content


# ── Reproducibility bundle ─────────────────────────────────────────────────────


class TestReproducibility:
    def test_checksums_file_created(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "reproducibility" / "checksums.sha256").exists()

    def test_environment_yml_has_dependencies(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        content = (tmp_path / "reproducibility" / "environment.yml").read_text()
        assert "biopython" in content or "iqtree" in content

    def test_commands_sh_references_script(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        content = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert "phylogenetics_builder" in content


# ── Output contract (SKILL.md → actual files) ─────────────────────────────────


def _parse_output_contract(skill_md: Path) -> list[str]:
    """Extract files promised in SKILL.md '## Output Structure' tree."""
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files = []
    parents: dict[int, str] = {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry = parts[0]
        comment = parts[1] if len(parts) > 1 else ""
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
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
    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = run_cli(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md promises artifacts the skill did not produce: "
            + ", ".join(missing)
        )


# ── Integration tests (require actual binaries) ────────────────────────────────


@pytest.mark.integration
class TestIntegrationIQTREE:
    """Require iqtree2 in PATH.  Skip with: pytest -m 'not integration'."""

    def test_full_pipeline_aligned_input(self, tmp_path):
        import shutil

        if not shutil.which("iqtree2") and not shutil.which("iqtree"):
            pytest.skip("iqtree2 not in PATH")
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--aligned",
                "--no-trim",
                "--bootstrap",
                "ufboot",
            ]
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "phylo_tree.nwk").exists()
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["engine"] == "iqtree2"


@pytest.mark.integration
class TestIntegrationRAxML:
    """Require raxml-ng in PATH."""

    def test_raxml_pipeline(self, tmp_path):
        import shutil

        if not shutil.which("raxml-ng"):
            pytest.skip("raxml-ng not in PATH")
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--aligned",
                "--no-trim",
                "--engine",
                "raxml-ng",
                "--bootstrap",
                "standard",
            ]
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["engine"] == "raxml-ng"


@pytest.mark.integration
class TestIntegrationMSA:
    """Require mafft in PATH."""

    def test_mafft_msa_then_iqtree(self, tmp_path):
        import shutil

        if not shutil.which("mafft"):
            pytest.skip("mafft not in PATH")
        if not shutil.which("iqtree2") and not shutil.which("iqtree"):
            pytest.skip("iqtree2 not in PATH")
        result = run_cli(
            [
                "--input",
                str(DEMO_INPUT),
                "--output",
                str(tmp_path),
                "--aligner",
                "mafft",
            ]
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["aligner"] == "mafft"
        assert data["aligned"] is False


# ── Alignment export ───────────────────────────────────────────────────────────


class TestAlignmentExport:
    """Alignments produced in the temporary working directory are kept."""

    def _run(self, tmp_path, extra_args, trimal_fails=False):
        m = get_module()
        newick, _ = m.get_demo_fallback()

        def fake_msa(input_fasta, aligner, output_fasta):
            output_fasta.write_text(DEMO_INPUT.read_text())

        def fake_trimal(input_fasta, output_fasta, strategy="-automated1"):
            if trimal_fails:
                raise RuntimeError("trimAl failed")
            output_fasta.write_text(DEMO_INPUT.read_text().replace("-", ""))

        out = tmp_path / "out"
        with patch.object(m.shutil, "which", return_value="/usr/bin/true"), \
                patch.object(m, "run_msa", side_effect=fake_msa), \
                patch.object(m, "run_trimal", side_effect=fake_trimal), \
                patch.object(m, "run_modelfinder", return_value="GTR+G"), \
                patch.object(m, "run_iqtree_main", return_value=newick):
            m.main(["--input", str(DEMO_INPUT), "--output", str(out)] + extra_args)
        return out, json.loads((out / "result.json").read_text())

    def test_msa_and_trimmed_alignments_are_exported(self, tmp_path):
        out, result = self._run(tmp_path, [])
        assert result["alignment_files"] == {
            "aligned": "alignment/aligned.fasta",
            "trimmed": "alignment/trimmed.fasta",
        }
        assert (out / "alignment" / "aligned.fasta").read_text() == DEMO_INPUT.read_text()
        assert (out / "alignment" / "trimmed.fasta").is_file()
        artifact_paths = {a["path"] for a in result["preferred_artifacts"]}
        assert {"alignment/aligned.fasta", "alignment/trimmed.fasta"} <= artifact_paths
        checksums = (out / "reproducibility" / "checksums.sha256").read_text()
        assert "alignment/aligned.fasta" in checksums
        assert "alignment/trimmed.fasta" in checksums
        assert "`alignment/aligned.fasta`" in (out / "report.md").read_text()

    def test_pre_aligned_input_without_trimming_exports_nothing(self, tmp_path):
        out, result = self._run(tmp_path, ["--aligned", "--no-trim"])
        assert result["alignment_files"] == {}
        assert not (out / "alignment").exists()

    def test_rerun_without_alignment_stages_removes_earlier_alignments(self, tmp_path):
        self._run(tmp_path, [])
        out, result = self._run(tmp_path, ["--aligned", "--no-trim"])
        assert result["alignment_files"] == {}
        assert not (out / "alignment").exists()
        assert "alignment/" not in (out / "reproducibility" / "checksums.sha256").read_text()

    def test_rerun_with_failed_trimming_keeps_only_the_new_alignment(self, tmp_path):
        self._run(tmp_path, [])
        out, result = self._run(tmp_path, [], trimal_fails=True)
        assert result["alignment_files"] == {"aligned": "alignment/aligned.fasta"}
        assert not (out / "alignment" / "trimmed.fasta").exists()

    def test_rerun_on_an_exported_alignment_does_not_delete_its_own_input(self, tmp_path):
        """Reusing alignment/aligned.fasta as the input of a rerun into the same
        output directory must not delete the file before it is read."""
        out, _ = self._run(tmp_path, [])
        exported = out / "alignment" / "aligned.fasta"
        before = exported.read_text()
        m = get_module()
        newick, _ = m.get_demo_fallback()
        with patch.object(m.shutil, "which", return_value="/usr/bin/true"), \
                patch.object(m, "run_trimal", side_effect=lambda i, o, strategy="-automated1": o.write_text(before)), \
                patch.object(m, "run_modelfinder", return_value="GTR+G"), \
                patch.object(m, "run_iqtree_main", return_value=newick):
            m.main(["--input", str(exported), "--aligned", "--output", str(out)])
        assert exported.is_file() and exported.read_text() == before

