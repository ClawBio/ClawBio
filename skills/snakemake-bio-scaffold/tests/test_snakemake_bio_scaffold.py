"""Tests for snakemake-bio-scaffold.

The skill generates a Snakemake project in the config-by-concern layout
(data/analysis/software.yaml + config_loader.py, one rules/<stage>.smk and one
standalone scripts/<stage>.py per stage, results/done/ sentinels). These tests
check the generator's own contract and that the *generated* project is
internally consistent and runnable.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import py_compile
import re
import shutil
import subprocess
import sys
import tokenize
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "snakemake_bio_scaffold.py"
DEMO_SPEC = SKILL_DIR / "examples" / "demo_spec.yaml"
DEMO_PROJECT = "demo_sumstats_qc"
DEMO_STAGES = ["format_sumstats", "filter_maf", "filter_pvalue"]


def run_cli(args, **kwargs):
    return subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        capture_output=True,
        text=True,
        **kwargs,
    )


def _load_module():
    spec = importlib.util.spec_from_file_location("snakemake_bio_scaffold", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo")
    result = run_cli(["--demo", "--output", str(out)])
    assert result.returncode == 0, result.stderr
    return out


@pytest.fixture(scope="module")
def demo_project(demo_out):
    return demo_out / DEMO_PROJECT


def _tokens(path: Path) -> list[tokenize.TokenInfo]:
    with open(path, "rb") as fh:
        return list(tokenize.tokenize(fh.readline))


def _load_generated_config(project: Path, config_dir: Path | None = None):
    loader_path = project / "config" / "config_loader.py"
    spec = importlib.util.spec_from_file_location(f"cfg_loader_{id(project)}", loader_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, mod.load_config(str(config_dir) if config_dir else None)


# ── Generator: top-level outputs ───────────────────────────────────────────────


class TestDemoRun:
    def test_report_and_result_written(self, demo_out):
        assert (demo_out / "report.md").exists()
        data = json.loads((demo_out / "result.json").read_text())
        assert data["skill"] == "snakemake-bio-scaffold"
        assert data["project"] == DEMO_PROJECT
        assert [s["name"] for s in data["stages"]] == DEMO_STAGES
        assert data["targets"] == ["cohort_a", "cohort_b"]
        assert data["rule_count"] == len(DEMO_STAGES)
        assert "Snakefile" in data["files"]

    def test_reproducibility_bundle(self, demo_out):
        repro = demo_out / "reproducibility"
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (repro / name).exists(), name
        assert "--demo" in (repro / "commands.sh").read_text()

    def test_report_has_disclaimer_and_next_steps(self, demo_out):
        text = (demo_out / "report.md").read_text(encoding="utf-8")
        assert "not a medical device" in text
        assert "snakemake -n" in text
        for stage in DEMO_STAGES:
            assert stage in text


# ── Generated project: layout ─────────────────────────────────────────────────


class TestProjectLayout:
    def test_top_level_layout(self, demo_project):
        for rel in (
            "Snakefile",
            "README.md",
            ".gitignore",
            "config/data.yaml",
            "config/analysis.yaml",
            "config/software.yaml",
            "config/config_loader.py",
            "scripts/lib/pipeline_io.py",
            "utilities/README.md",
            "resources/README.md",
        ):
            assert (demo_project / rel).exists(), rel

    def test_one_rule_and_one_script_per_stage(self, demo_project):
        rules = sorted(p.stem for p in (demo_project / "rules").glob("*.smk"))
        scripts = sorted(p.stem for p in (demo_project / "scripts").glob("*.py"))
        assert rules == sorted(DEMO_STAGES)
        assert scripts == sorted(DEMO_STAGES)

    def test_demo_inputs_are_synthetic_tsvs(self, demo_project):
        for item in ("cohort_a", "cohort_b"):
            path = demo_project / "input" / f"{item}.tsv"
            lines = path.read_text().splitlines()
            assert lines[0].split("\t")[:3] == ["SNP", "CHR", "BP"]
            assert len(lines) > 50

    def test_gitignore_excludes_results(self, demo_project):
        text = (demo_project / ".gitignore").read_text()
        assert "results/" in text
        assert ".snakemake/" in text

    def test_all_python_compiles(self, demo_project):
        for py in demo_project.rglob("*.py"):
            py_compile.compile(str(py), doraise=True)

    def test_snakefile_and_rules_tokenize(self, demo_project):
        """Snakemake syntax is Python-tokenizable; this catches broken literals without Snakemake."""
        for path in [demo_project / "Snakefile", *(demo_project / "rules").glob("*.smk")]:
            _tokens(path)

    def test_all_yaml_parses(self, demo_project):
        for y in (demo_project / "config").glob("*.yaml"):
            yaml.safe_load(y.read_text())


# ── Generated project: conventions ────────────────────────────────────────────


class TestConventions:
    def test_every_rule_touches_a_done_sentinel(self, demo_project):
        for smk in (demo_project / "rules").glob("*.smk"):
            text = smk.read_text()
            assert re.search(r'done=touch\(f"\{OUT\}/done/' + smk.stem + r"_", text), smk.name

    def test_rules_use_sys_executable_not_bare_python(self, demo_project):
        snakefile = (demo_project / "Snakefile").read_text()
        assert "PYTHON = sys.executable" in snakefile
        for smk in (demo_project / "rules").glob("*.smk"):
            text = smk.read_text()
            assert "'{PYTHON:q} scripts/" in text
            assert "'python " not in text and '"python ' not in text

    def test_downstream_stage_depends_on_upstream_sentinel(self, demo_project):
        text = (demo_project / "rules" / "filter_maf.smk").read_text()
        assert 'f"{OUT}/done/format_sumstats_{{dataset}}.done"' in text

    def test_rule_all_targets_last_stage(self, demo_project):
        text = (demo_project / "Snakefile").read_text()
        assert 'f"{OUT}/done/filter_pvalue_{n}.done" for n in TARGETS' in text

    def test_wildcard_constraint_is_closed_set(self, demo_project):
        text = (demo_project / "Snakefile").read_text()
        assert "wildcard_constraints:" in text
        assert 'dataset="|".join(re.escape(' in text

    def test_scripts_do_not_import_snakemake(self, demo_project):
        for py in (demo_project / "scripts").rglob("*.py"):
            text = py.read_text()
            assert "import snakemake" not in text
            assert "snakemake." not in text

    def test_params_come_from_analysis_yaml(self, demo_project):
        analysis = yaml.safe_load((demo_project / "config" / "analysis.yaml").read_text())
        assert analysis["filter_maf"] == {"maf": 0.01}
        assert analysis["filter_pvalue"] == {"p_threshold": 5e-08}
        smk = (demo_project / "rules" / "filter_maf.smk").read_text()
        assert '_item_analysis(wc.dataset, "filter_maf")["maf"]' in smk
        assert "--maf {params.maf:q}" in smk


# ── Generated project: config_loader.py ───────────────────────────────────────


class TestGeneratedConfigLoader:
    def test_loads_and_resolves_paths(self, demo_project):
        _, cfg = _load_generated_config(demo_project)
        items = cfg["data"]["datasets"]
        assert set(items) == {"cohort_a", "cohort_b"}
        path = items["cohort_a"]["path"]
        assert Path(path).is_absolute()
        assert "\\" not in path
        assert "\\" not in cfg["analysis"]["output_dir"]

    def test_per_item_overrides_merge(self, demo_project):
        _, cfg = _load_generated_config(demo_project)
        a = cfg["data"]["datasets"]["cohort_a"]["analysis"]
        b = cfg["data"]["datasets"]["cohort_b"]["analysis"]
        assert a["filter_maf"]["maf"] == 0.01
        assert b["filter_maf"]["maf"] == 0.05  # demo spec overrides cohort_b

    def _copy_config(self, demo_project, tmp_path):
        dst = tmp_path / "proj"
        shutil.copytree(demo_project, dst)
        return dst

    def test_missing_input_file_raises_config_error(self, demo_project, tmp_path):
        proj = self._copy_config(demo_project, tmp_path)
        (proj / "input" / "cohort_a.tsv").unlink()
        mod_path = proj / "config" / "config_loader.py"
        spec = importlib.util.spec_from_file_location("cfg_missing", mod_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with pytest.raises(mod.ConfigError, match="not found"):
            mod.load_config()

    def test_unknown_target_raises_config_error(self, demo_project, tmp_path):
        proj = self._copy_config(demo_project, tmp_path)
        analysis_path = proj / "config" / "analysis.yaml"
        analysis = yaml.safe_load(analysis_path.read_text())
        analysis["targets"].append("no_such_item")
        analysis_path.write_text(yaml.safe_dump(analysis))
        spec = importlib.util.spec_from_file_location("cfg_target", proj / "config" / "config_loader.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with pytest.raises(mod.ConfigError, match="no_such_item"):
            mod.load_config()

    def test_typo_in_override_raises_config_error(self, demo_project, tmp_path):
        proj = self._copy_config(demo_project, tmp_path)
        data_path = proj / "config" / "data.yaml"
        data = yaml.safe_load(data_path.read_text())
        data["datasets"]["cohort_a"]["overrides"] = {"filter_maf": {"mfa": 0.2}}
        data_path.write_text(yaml.safe_dump(data))
        spec = importlib.util.spec_from_file_location("cfg_typo", proj / "config" / "config_loader.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with pytest.raises(mod.ConfigError, match="mfa"):
            mod.load_config()


# ── Generated project: stage scripts are standalone CLIs ──────────────────────


class TestStageScripts:
    def test_stub_runs_standalone_and_passes_rows_through(self, demo_project, tmp_path):
        src = demo_project / "input" / "cohort_a.tsv"
        out = tmp_path / "x.tsv"
        summary = tmp_path / "x.summary.json"
        result = subprocess.run(
            [
                sys.executable,
                str(demo_project / "scripts" / "filter_maf.py"),
                "--input", str(src),
                "--out", str(out),
                "--summary-json", str(summary),
                "--maf", "0.01",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert out.read_text() == src.read_text()
        s = json.loads(summary.read_text())
        assert s["stage"] == "filter_maf"
        assert s["rows_in"] == s["rows_out"] > 0
        assert s["params"] == {"maf": 0.01}

    def test_stub_has_docstring_example(self, demo_project):
        text = (demo_project / "scripts" / "filter_pvalue.py").read_text()
        assert "Example:" in text
        assert "--p-threshold" in text


# ── Spec validation and quick mode ────────────────────────────────────────────


class TestSpecValidation:
    def _write_spec(self, tmp_path, spec):
        p = tmp_path / "spec.yaml"
        p.write_text(yaml.safe_dump(spec))
        return p

    @pytest.mark.parametrize(
        "stages, message",
        [
            ([], "at least one stage"),
            ([{"name": "bad-name"}], "bad-name"),
            ([{"name": "a"}, {"name": "a"}], "duplicate"),
            ([{"name": "all"}], "reserved"),
            ([{"name": "a", "params": {"x": None}}], "null"),
            ([{"name": "a", "params": {"bad-key": 1}}], "bad-key"),
        ],
    )
    def test_invalid_specs_rejected(self, tmp_path, stages, message):
        spec = self._write_spec(tmp_path, {"project": "p", "stages": stages})
        result = run_cli(["--input", str(spec), "--output", str(tmp_path / "out")])
        assert result.returncode != 0
        assert message in result.stderr

    def test_invalid_project_name_rejected(self, tmp_path):
        spec = self._write_spec(tmp_path, {"project": "../escape", "stages": [{"name": "a"}]})
        result = run_cli(["--input", str(spec), "--output", str(tmp_path / "out")])
        assert result.returncode != 0
        assert "project" in result.stderr

    def test_custom_wildcard_and_items_key(self, tmp_path):
        spec = self._write_spec(
            tmp_path,
            {
                "project": "rna",
                "items": {"key": "samples", "wildcard": "sample",
                          "entries": {"s1": {"path": "input/s1.tsv"}}},
                "stages": [{"name": "count"}],
            },
        )
        result = run_cli(["--input", str(spec), "--output", str(tmp_path / "out")])
        assert result.returncode == 0, result.stderr
        proj = tmp_path / "out" / "rna"
        assert "samples" in yaml.safe_load((proj / "config" / "data.yaml").read_text())
        assert "{{sample}}" in (proj / "rules" / "count.smk").read_text()

    def test_refuses_non_empty_project_dir(self, tmp_path):
        out = tmp_path / "out"
        (out / DEMO_PROJECT).mkdir(parents=True)
        (out / DEMO_PROJECT / "keep.txt").write_text("user file")
        result = run_cli(["--demo", "--output", str(out)])
        assert result.returncode != 0
        assert "--force" in result.stderr
        assert (out / DEMO_PROJECT / "keep.txt").read_text() == "user file"

    def test_force_overwrites_scaffold_but_keeps_other_files(self, tmp_path):
        out = tmp_path / "out"
        (out / DEMO_PROJECT).mkdir(parents=True)
        (out / DEMO_PROJECT / "keep.txt").write_text("user file")
        result = run_cli(["--demo", "--output", str(out), "--force"])
        assert result.returncode == 0, result.stderr
        assert (out / DEMO_PROJECT / "Snakefile").exists()
        assert (out / DEMO_PROJECT / "keep.txt").exists()

    def test_quick_mode_from_flags(self, tmp_path):
        result = run_cli(
            ["--name", "quick", "--stages", "align,count,report", "--output", str(tmp_path)]
        )
        assert result.returncode == 0, result.stderr
        rules = sorted(p.stem for p in (tmp_path / "quick" / "rules").glob("*.smk"))
        assert rules == ["align", "count", "report"]

    def test_no_input_and_no_demo_errors(self, tmp_path):
        result = run_cli(["--output", str(tmp_path)])
        assert result.returncode != 0


# ── Hostile specs: spec text must never become code (PR #525 review) ──────────

# Every payload creates a PWNED_* file if it ever runs. Each payload keeps
# "open(" and "PWNED" together, so a payload that stays inert sits inside one
# STRING or COMMENT token.
_PAYLOAD = "open('PWNED_{}', 'w')"
_BREAKOUTS = ['"""\n{p}\n"""', "'''\n{p}\n'''", '\\"); {p}; ("', '\\', '$({p})', '`{p}`', '@@ARGPARSE@@ @@INCLUDES@@']


def _hostile(tag: str) -> str:
    payload = _PAYLOAD.format(tag)
    return " ".join(b.format(p=payload) for b in _BREAKOUTS) + " " + payload + " \\"


def _hostile_line(tag: str) -> str:
    """Single-line variant, for values that reach a command line (newlines are rejected there)."""
    return _hostile(tag).replace("\n", " ")


def _hostile_spec(path: str | None = None) -> dict:
    return {
        "project": "hostile",
        "description": _hostile("desc"),
        "items": {"entries": {"a": {"path": path or "input/a " + _hostile_line("path") + ".tsv",
                                    "label": _hostile("label"), "description": _hostile("itemdesc")}}},
        "stages": [
            {"name": "first", "description": _hostile("stage1"),
             "params": {"txt": _hostile_line("param"), "trail": "ends with backslash \\", "n": 3}},
            {"name": "second", "description": "trailing backslash \\", "params": {"flag": True}},
        ],
    }


def _scaffold(tmp_path, spec, *extra):
    p = tmp_path / "spec.yaml"
    p.write_text(yaml.safe_dump(spec), encoding="utf-8")
    out = tmp_path / "out"
    result = run_cli(["--input", str(p), "--output", str(out), *extra])
    assert result.returncode == 0, result.stderr
    return out / spec["project"]


def _pwned(root: Path) -> list[Path]:
    return [p for p in root.rglob("PWNED_*")]


class TestHostileSpec:
    def test_generated_sources_keep_spec_text_inert(self, tmp_path):
        proj = _scaffold(tmp_path, _hostile_spec())
        sources = [proj / "Snakefile", *proj.glob("rules/*.smk"), *proj.rglob("*.py")]
        assert len(sources) >= 6
        for path in sources:
            for tok in _tokens(path):
                if "PWNED" in tok.string:
                    assert tok.type in (tokenize.STRING, tokenize.COMMENT), (path.name, tok)
                    assert "open(" in tok.string, (path.name, tok)
        for py in proj.rglob("*.py"):
            ast.parse(py.read_text(encoding="utf-8"), filename=str(py))

    def test_stage_scripts_help_runs_without_side_effects(self, tmp_path):
        proj = _scaffold(tmp_path, _hostile_spec())
        for script in proj.glob("scripts/*.py"):
            r = subprocess.run([sys.executable, str(script), "--help"], cwd=proj, capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
        assert _pwned(tmp_path) == []

    def test_spec_text_round_trips_into_yaml(self, tmp_path):
        spec = _hostile_spec()
        proj = _scaffold(tmp_path, spec)
        analysis = yaml.safe_load((proj / "config" / "analysis.yaml").read_text(encoding="utf-8"))
        assert analysis["first"]["txt"] == spec["stages"][0]["params"]["txt"]
        assert analysis["first"]["trail"].endswith("\\")
        data = yaml.safe_load((proj / "config" / "data.yaml").read_text(encoding="utf-8"))
        # Item paths are normalised to forward slashes; otherwise verbatim.
        assert data["datasets"]["a"]["path"] == spec["items"]["entries"]["a"]["path"].replace("\\", "/")

    @pytest.mark.parametrize("value", ["/abs/results", "C:/results", "C:\\results", "../results",
                                       "results/../../x", "", "\\\\server\\share"])
    def test_output_dir_must_stay_inside_project(self, tmp_path, value):
        p = tmp_path / "spec.yaml"
        p.write_text(yaml.safe_dump({"project": "p", "stages": ["a"], "output_dir": value}))
        r = run_cli(["--input", str(p), "--output", str(tmp_path / "out")])
        assert r.returncode != 0
        assert "output_dir" in r.stderr

    @pytest.mark.parametrize("where", ["param", "path"])
    @pytest.mark.parametrize("char", ["\n", "\r", "\x00", "\x1b"])
    def test_control_characters_rejected_in_command_line_values(self, tmp_path, where, char):
        """Params and item paths become command-line arguments; on Windows cmd a newline ends the command."""
        spec = {"project": "p", "stages": [{"name": "a", "params": {"x": "ok"}}],
                "items": {"entries": {"i": {"path": "input/i.tsv"}}}}
        if where == "param":
            spec["stages"][0]["params"]["x"] = f"line1{char}line2"
        else:
            spec["items"]["entries"]["i"]["path"] = f"input/i{char}.tsv"
        p = tmp_path / "spec.yaml"
        p.write_text(yaml.safe_dump(spec))
        r = run_cli(["--input", str(p), "--output", str(tmp_path / "out")])
        assert r.returncode != 0
        assert "control character" in r.stderr

    def test_relative_output_dir_accepted(self, tmp_path):
        p = tmp_path / "spec.yaml"
        p.write_text(yaml.safe_dump({"project": "p", "stages": ["a"], "output_dir": "results/run1"}))
        assert run_cli(["--input", str(p), "--output", str(tmp_path / "out")]).returncode == 0

    def test_shell_arguments_use_snakemake_quoting(self, tmp_path):
        proj = _scaffold(tmp_path, _hostile_spec())
        text = (proj / "rules" / "first.smk").read_text(encoding="utf-8")
        for field in ("{input.data:q}", "{output.result:q}", "{output.summary:q}", "{params.txt:q}",
                      "{params.trail:q}", "{params.n:q}", "{log:q}"):
            assert field in text, field
        assert '"{input.data}"' not in text and '"{params.txt}"' not in text


class TestOutputGuard:
    def test_refuses_to_overwrite_existing_report(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "report.md").write_text("earlier report")
        r = run_cli(["--name", "q", "--stages", "a", "--output", str(out)])
        assert r.returncode != 0
        assert "--force" in r.stderr
        assert (out / "report.md").read_text() == "earlier report"
        assert not (out / "q").exists()

    @pytest.mark.parametrize("existing", ["result.json", "reproducibility"])
    def test_refuses_other_existing_outputs(self, tmp_path, existing):
        out = tmp_path / "out"
        (out / existing).mkdir(parents=True)
        r = run_cli(["--name", "q", "--stages", "a", "--output", str(out)])
        assert r.returncode != 0

    def test_force_overwrites_report(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "report.md").write_text("earlier report")
        r = run_cli(["--name", "q", "--stages", "a", "--output", str(out), "--force"])
        assert r.returncode == 0, r.stderr
        assert "Snakemake Bio Scaffold Report" in (out / "report.md").read_text(encoding="utf-8")

    def test_unrelated_files_in_output_are_fine(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "notes.txt").write_text("mine")
        assert run_cli(["--name", "q", "--stages", "a", "--output", str(out)]).returncode == 0


class TestModuleApi:
    def test_validate_spec_fills_defaults(self):
        mod = _load_module()
        spec = mod.validate_spec({"project": "p", "stages": [{"name": "a"}]})
        assert spec["items"]["key"] == "datasets"
        assert spec["items"]["wildcard"] == "dataset"
        assert spec["stages"][0]["ext"] == "tsv"
        assert spec["stages"][0]["params"] == {}


# ── Output contract ───────────────────────────────────────────────────────────


def _parse_output_contract(skill_md: Path) -> list[str]:
    """Extract files promised in SKILL.md '## Output Structure' tree."""
    text = skill_md.read_text(encoding="utf-8")
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
        files.append(f"{rel}/{name}" if rel else name)
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md must document an Output Structure tree"
        result = run_cli(["--demo", "--output", str(tmp_path)])
        assert result.returncode == 0, result.stderr
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, "SKILL.md promises artifacts the skill did not produce: " + ", ".join(missing)


# ── Integration: the generated project actually runs under Snakemake ──────────


def _snakemake_available() -> bool:
    return importlib.util.find_spec("snakemake") is not None


@pytest.mark.integration
@pytest.mark.skipif(not _snakemake_available(), reason="snakemake not installed")
class TestSnakemakeIntegration:
    def test_dry_run(self, demo_project):
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "-n", "--cores", "1"],
            cwd=demo_project,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    def test_full_run_produces_sentinels(self, demo_project, tmp_path):
        proj = tmp_path / "proj"
        shutil.copytree(demo_project, proj)
        result = subprocess.run(
            # --drop-metadata: .snakemake/metadata filenames are base64 of the full
            # output path, which overflows Windows MAX_PATH under pytest's deep tmp dirs.
            [sys.executable, "-m", "snakemake", "--cores", "1", "--drop-metadata"],
            cwd=proj,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        for stage in DEMO_STAGES:
            for item in ("cohort_a", "cohort_b"):
                assert (proj / "results" / "done" / f"{stage}_{item}.done").exists()
                assert (proj / "results" / stage / f"{item}.summary.json").exists()

    def test_hostile_spec_runs_without_side_effects(self, tmp_path):
        """Full run of a hostile spec: nothing executes, and odd param values reach the script verbatim."""
        spec = _hostile_spec(path="input/a.tsv")
        odd = "a b $(echo X) `echo Y` \"q\" 'r' ; & | > x"
        spec["stages"][1]["params"]["odd"] = odd
        proj = _scaffold(tmp_path, spec)
        (proj / "input" / "a.tsv").write_text("SNP\tP\nrs1\t0.5\n")
        dry = subprocess.run([sys.executable, "-m", "snakemake", "-n", "--cores", "1"],
                             cwd=proj, capture_output=True, text=True)
        assert dry.returncode == 0, dry.stdout + dry.stderr
        assert _pwned(tmp_path) == []
        result = subprocess.run(
            [sys.executable, "-m", "snakemake", "--cores", "1", "--drop-metadata"],
            cwd=proj, capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert _pwned(tmp_path) == []
        first = json.loads((proj / "results" / "first" / "a.summary.json").read_text())
        assert first["params"]["txt"] == spec["stages"][0]["params"]["txt"]
        assert first["params"]["trail"] == "ends with backslash \\"
        second = json.loads((proj / "results" / "second" / "a.summary.json").read_text())
        assert second["params"]["odd"] == odd

    def test_check_flag_records_dry_run(self, tmp_path):
        result = run_cli(["--demo", "--output", str(tmp_path), "--check"])
        assert result.returncode == 0, result.stderr
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["dry_run"]["status"] == "passed"
