"""Licence policy for skills, enforced against the generated catalogue.

Why this exists: a PROPRIETARY skill and a GPL-3.0 skill inside an MIT wheel is
the first thing a procurement reader hits, and a blank licence is worse because
nobody can tell what they are allowed to do. The policy gates the CODE licence
only. Data licences are reported, not gated, because they are declared by 97
contributors and a gate none of them can pass today is a gate that gets
disabled tomorrow.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = ROOT / "skills"
CATALOG = SKILLS_DIR / "catalog.json"

sys.path.insert(0, str(ROOT / "scripts"))
from generate_catalog import parse_yaml_frontmatter  # noqa: E402

# Code licences that may ship inside the MIT wheel.
WHEEL_LICENCES = {"MIT", "Apache-2.0"}

# Skills whose code licence is incompatible with the wheel. They remain in the
# repository as folders labelled by licence and are excluded from the wheel.
# Every entry needs a reason so the exception cannot grow silently.
FOLDER_ONLY = {
    "fastreer": "GPL-3.0: copyleft cannot be redistributed inside an MIT wheel",
    "wes-clinical-report-en": "PROPRIETARY: declared by the contributor; not redistributable",
    "wes-clinical-report-es": "PROPRIETARY: declared by the contributor; not redistributable",
}


def _licences_on_disk() -> dict[str, str]:
    """Licences read from disk, not from the catalogue: generate_catalog.py
    drops folders listed in EXCLUDED_FOLDERS, so a catalogue-only check cannot
    see the skills most likely to be unredistributable."""
    return {
        d.name: str(parse_yaml_frontmatter((d / "SKILL.md").read_text(encoding="utf-8")).get("license", ""))
        for d in sorted(SKILLS_DIR.iterdir())
        if d.is_dir() and (d / "SKILL.md").is_file()
    }


def _wheel_excluded_skills() -> set[str]:
    """Read WHEEL_EXCLUDED_SKILLS from the AST rather than importing the module:
    hatch_build imports hatchling, a build-time dependency that need not be
    installed wherever the tests run."""
    import ast
    tree = ast.parse((ROOT / "hatch_build.py").read_text())
    values = [
        ast.literal_eval(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "WHEEL_EXCLUDED_SKILLS" for t in node.targets)
    ]
    assert values, "WHEEL_EXCLUDED_SKILLS not defined in hatch_build.py"
    return set(values[0])


def _skills() -> list[dict]:
    return json.loads(CATALOG.read_text())["skills"]


def test_every_skill_declares_a_code_licence():
    blank = [s["name"] for s in _skills() if not (s.get("license") or "").strip()]
    assert blank == [], f"skills with no code licence: {blank}"


def test_wheel_skills_carry_a_wheel_compatible_licence():
    offenders = [
        (s["name"], s.get("license"))
        for s in _skills()
        if s["name"] not in FOLDER_ONLY and (s.get("license") or "") not in WHEEL_LICENCES
    ]
    assert offenders == [], f"licences that cannot ship in the MIT wheel: {offenders}"


def test_folder_only_exceptions_are_real_and_still_needed():
    """Anchored to the folder on disk, not to the catalogue: an exception may
    legitimately be absent from catalog.json (generate_catalog.py has its own
    EXCLUDED_FOLDERS), and that is precisely when the wheel needs gating."""
    on_disk = _licences_on_disk()
    for name, reason in FOLDER_ONLY.items():
        assert name in on_disk, f"{name} is listed as an exception but has no skills/{name}/SKILL.md"
        lic = on_disk[name]
        assert lic not in WHEEL_LICENCES, (
            f"{name} is now {lic!r}; remove it from FOLDER_ONLY, the exception is stale")
        assert lic.split(":")[0].strip() in reason, f"{name}: reason must name the licence"


def test_wheel_exclusion_matches_the_policy():
    """Three views of one policy must agree: the declared exceptions, the build
    hook's constant, and what the folders on disk actually say. The build hook
    also derives exclusions from disk at build time, so this is what stops the
    declaration drifting silently behind that gate."""
    declared = _wheel_excluded_skills()
    assert declared == set(FOLDER_ONLY)
    on_disk = {name for name, lic in _licences_on_disk().items()
               if lic and lic not in WHEEL_LICENCES}
    assert declared == on_disk, (
        f"declared exceptions {sorted(declared)} do not match the folders on "
        f"disk {sorted(on_disk)}; update FOLDER_ONLY and WHEEL_EXCLUDED_SKILLS")


def test_data_licence_is_reported_not_gated(capsys):
    """Published count, so the gap is visible on every run without blocking anyone."""
    undeclared = [s["name"] for s in _skills()
                  if (s.get("data_license") or "none").lower() in {"none", ""}]
    print(f"\n[licence policy] data_license undeclared on {len(undeclared)} of {len(_skills())} skills")
    assert isinstance(undeclared, list)


def test_wheel_excludes_every_unredistributable_folder_on_disk():
    """The wheel is built by walking skills/ on disk, so the policy has to be
    checked against disk too. A folder kept out of catalog.json by
    EXCLUDED_FOLDERS is still bundled by hatch_build.py."""
    excluded = _wheel_excluded_skills()
    offenders = [
        (name, lic)
        for name, lic in _licences_on_disk().items()
        if lic and lic not in WHEEL_LICENCES and name not in excluded
    ]
    assert offenders == [], f"unredistributable skills that would ship in the MIT wheel: {offenders}"


def test_every_skill_folder_on_disk_declares_a_code_licence():
    """A blank licence on disk is the same hole as a non-wheel licence: the
    build hook has nothing to gate on."""
    blank = [name for name, lic in _licences_on_disk().items() if not lic]
    assert blank == [], f"skill folders with no code licence in SKILL.md: {blank}"


def test_cli_registered_skills_survive_the_wheel_filter():
    """A skill can be registered as a CLI action or excluded from the wheel, but
    if it is both, `clawbio run <skill>` is broken for every pip user. Any such
    skill must be declared in cli.UNBUNDLED_SKILLS so the CLI can say why."""
    sys.path.insert(0, str(ROOT))
    from clawbio import cli

    excluded = _wheel_excluded_skills()
    registered = {
        alias
        for alias, info in cli.SKILLS.items()
        if "script" in info and SKILLS_DIR in Path(info["script"]).parents
    }
    broken = {
        alias
        for alias in registered
        if Path(cli.SKILLS[alias]["script"]).parent.name in excluded
    }
    assert broken <= set(cli.UNBUNDLED_SKILLS), (
        f"CLI-registered skills excluded from the wheel but not declared in "
        f"UNBUNDLED_SKILLS: {sorted(broken - set(cli.UNBUNDLED_SKILLS))}")
    assert set(cli.UNBUNDLED_SKILLS) == broken, (
        f"UNBUNDLED_SKILLS is stale; it names {sorted(set(cli.UNBUNDLED_SKILLS) - broken)} "
        f"which are not excluded from the wheel")
    for alias, reason in cli.UNBUNDLED_SKILLS.items():
        lic = _licences_on_disk()[Path(cli.SKILLS[alias]["script"]).parent.name]
        assert lic.split("-")[0] in reason, f"{alias}: reason must name the licence"


def test_unbundled_skill_reports_why_it_is_missing(monkeypatch, tmp_path):
    """The pip user gets the licence reason, not a bare path that does not exist."""
    sys.path.insert(0, str(ROOT))
    from clawbio import cli

    alias = next(iter(cli.UNBUNDLED_SKILLS))
    missing = tmp_path / "not-in-the-wheel" / f"{alias}.py"
    monkeypatch.setitem(cli.SKILLS[alias], "script", missing)

    result = cli.run_skill(alias, demo=True)
    assert result["success"] is False
    assert "not bundled" in result["stderr"]
    assert cli.UNBUNDLED_SKILLS[alias].split(":")[0] in result["stderr"]
    assert "github.com" in result["stderr"] or "checkout" in result["stderr"]
