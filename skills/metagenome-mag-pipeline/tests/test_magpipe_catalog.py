"""The skills/catalog.json entry must mirror SKILL.md and the CLI registration.

The catalog is the machine-readable index an agent routes on, so it must not
drift from the skill's own description/version/tags/trigger_keywords, and the
CLI alias the catalog advertises must be the alias `clawbio.py run` accepts.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parent.parent
CATALOG = SKILL_DIR.parent / "catalog.json"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from clawbio.cli import SKILLS

ALIAS = "mag-pipeline"
FOLDER = "metagenome-mag-pipeline"


def _frontmatter() -> dict:
    text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert text.startswith("---"), "SKILL.md must open with a YAML frontmatter block"
    _, block, _ = text.split("---", 2)
    return yaml.safe_load(block)


def _entry() -> dict:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    items = catalog if isinstance(catalog, list) else catalog.get("skills", [])
    return next(entry for entry in items if entry.get("name") == FOLDER)


def test_catalog_entry_mirrors_skill_frontmatter() -> None:
    frontmatter = _frontmatter()
    metadata = frontmatter.get("metadata", {})
    entry = _entry()
    assert entry["description"] == frontmatter["description"]
    assert entry["version"] == metadata["version"]
    assert entry["tags"] == metadata["tags"]
    assert entry["trigger_keywords"] == metadata["openclaw"]["trigger_keywords"]
    assert entry["license"] == frontmatter["license"]


def test_catalog_entry_advertises_the_cli_alias() -> None:
    entry = _entry()
    assert entry["cli_alias"] == ALIAS
    assert entry["has_script"] is True
    assert entry["has_tests"] is True
    # The demo command must carry the fetch-consent flag: the demo clones the
    # pinned third-party workflow when the cache is empty.
    assert entry["demo_command"] == f"python clawbio.py run {ALIAS} --demo --allow-fetch"


def test_alias_is_registered_with_the_same_script() -> None:
    info = SKILLS[ALIAS]
    assert Path(info["script"]) == SKILL_DIR / "metagenome_mag_pipeline.py"
    # --allow-fetch: the demo clones the pinned third-party workflow, so the
    # registered demo command carries the consent flag explicitly.
    assert info["demo_args"] == ["--demo", "--allow-fetch"]
    assert info["accepts_genotypes"] is False
    assert Path(info["script"]).is_file()


def test_chaining_partners_are_registered_skills() -> None:
    entry = _entry()
    assert entry["chaining_partners"], "the skill chains with claw-metagenomics and others"
    for partner in entry["chaining_partners"]:
        assert (SKILL_DIR.parent / partner / "SKILL.md").is_file(), (
            f"{partner} is advertised as a chaining partner but has no SKILL.md"
        )


def test_catalog_counts_include_this_skill() -> None:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8"))
    assert catalog["skill_count"] == len(catalog["skills"])
    assert any(entry["name"] == FOLDER for entry in catalog["skills"])


# Flags whose value is a filesystem path. clawbio.py runs the script from the skill
# folder, so any of these missing from extra_path_flags turns `--host-fasta ref.fa`
# into skills/metagenome-mag-pipeline/ref.fa. (--input is resolved by the runner.)
PATH_FLAGS = {
    "--host-fasta", "--checkm2-db", "--gtdbtk-db", "--taxonomy-db", "--taxdump",
    "--kegg-db", "--cog-db", "--pfam-db", "--pipeline-dir", "--conda-prefix",
}


def test_every_path_flag_is_resolved_by_the_runner() -> None:
    sys.path.insert(0, str(SKILL_DIR))
    try:
        import importlib

        parser = importlib.import_module("metagenome_mag_pipeline").build_parser()
    finally:
        sys.path.remove(str(SKILL_DIR))
    flags = {opt for action in parser._actions for opt in action.option_strings if opt.startswith("--")}
    looks_like_a_path = {f for f in flags if f.endswith(("-db", "-dir", "-fasta", "-prefix")) or f == "--taxdump"}
    assert looks_like_a_path == PATH_FLAGS
    assert SKILLS[ALIAS].get("extra_path_flags") == PATH_FLAGS
    assert PATH_FLAGS <= SKILLS[ALIAS]["allowed_extra_flags"]


def test_relative_path_flag_through_clawbio_run_resolves_against_cwd(tmp_path: Path) -> None:
    import subprocess

    (tmp_path / "ref.fa").write_text(">h\nACGT\n")
    probe = tmp_path / "probe.py"
    # Stand-in skill script that reports the --host-fasta value it was given.
    probe.write_text("import sys\nargs = sys.argv\nprint('HOST=' + args[args.index('--host-fasta') + 1])\n")
    code = (
        "import sys; from pathlib import Path; import clawbio.cli as cli\n"
        f"cli.SKILLS['{ALIAS}'] = dict(cli.SKILLS['{ALIAS}'], script=Path({str(probe)!r}))\n"
        f"r = cli.run_skill('{ALIAS}', demo=True, extra_args=['--host-fasta', 'ref.fa'], "
        f"output_dir={str(tmp_path / 'out')!r})\n"
        "print(r.get('stdout', ''))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=tmp_path, capture_output=True, text=True,
                            env={"PYTHONPATH": str(REPO_ROOT), "PATH": "/usr/bin:/bin"})
    assert f"HOST={tmp_path / 'ref.fa'}" in result.stdout, result.stdout + result.stderr
