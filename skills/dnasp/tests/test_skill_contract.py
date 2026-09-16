"""The SKILL.md contract must be usable by an agent without the help text.

An agent that reads only the registered metadata has to turn a logical input name
into a command line. Nothing in the metadata used to say which flag an input maps
to, so both operator models in the publication benchmark translated `window_size`
literally and produced `--window-size`, which does not exist. These tests hold the
mapping in place: every declared input names a real flag, and no flag the metadata
claims can drift away from the parser.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import dnasp  # noqa: E402


def _front_matter() -> str:
    text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
    assert match, "SKILL.md has no YAML front matter"
    return match.group(1)


def _declared_inputs() -> list[dict[str, str]]:
    """Read the inputs block the way the catalog generator does.

    PyYAML is a declared dependency, so this imports it rather than skipping:
    a skipped contract test is indistinguishable from a passing one in CI.
    """
    import yaml

    data = yaml.safe_load(_front_matter())
    return data["metadata"]["inputs"]


def _parser_flags() -> set[str]:
    flags: set[str] = set()
    for action in dnasp.build_parser()._actions:
        flags.update(action.option_strings)
    return flags


def test_every_declared_input_names_its_flag() -> None:
    missing = [i["name"] for i in _declared_inputs() if not i.get("cli_flag")]
    assert not missing, (
        "inputs in SKILL.md with no cli_flag, so an agent must guess the flag: "
        + ", ".join(missing)
    )


def test_declared_flags_exist_in_the_parser() -> None:
    flags = _parser_flags()
    wrong = [
        (i["name"], i["cli_flag"])
        for i in _declared_inputs()
        if i.get("cli_flag") and i["cli_flag"] not in flags
    ]
    assert not wrong, f"SKILL.md names flags the parser does not accept: {wrong}"


def test_the_guessable_mistake_is_not_a_real_flag() -> None:
    """Guards the specific failure the benchmark found."""
    flags = _parser_flags()
    assert "--window-size" not in flags and "--step-size" not in flags
    assert "--window" in flags and "--step" in flags


def test_no_flag_is_claimed_by_two_inputs() -> None:
    claimed = [i["cli_flag"] for i in _declared_inputs() if i.get("cli_flag")]
    assert len(claimed) == len(set(claimed)), f"duplicate cli_flag: {claimed}"


def test_version_is_consistent_across_the_skill() -> None:
    """A partial bump is the recurring mistake: metadata, any top-level version
    key, and the script constant must agree."""
    import yaml

    data = yaml.safe_load(_front_matter())
    versions = {"metadata.version": str(data["metadata"]["version"]),
                "dnasp.__version__": dnasp.__version__}
    if "version" in data:
        versions["frontmatter.version"] = str(data["version"])
    assert len(set(versions.values())) == 1, f"versions disagree: {versions}"


# The mapping is asserted explicitly, not merely checked for existence. Four tests
# that only ask "does this flag exist" would pass if window_size and step_size were
# swapped, or alignment mapped to --vcf. This table is the contract.
EXPECTED_FLAGS = {
    "alignment": "--input",
    "vcf": "--vcf",
    "alignment2": "--input2",
    "pop_file": "--pop-file",
    "outgroup": "--outgroup",
    "hka_file": "--hka-file",
    "analyses": "--analysis",
    "window_size": "--window",
    "step_size": "--step",
    "genetic_code": "--genetic-code",
}


def test_each_input_maps_to_the_right_flag() -> None:
    declared = {i["name"]: i.get("cli_flag") for i in _declared_inputs()}
    assert declared == EXPECTED_FLAGS, (
        "the metadata's input-to-flag mapping is not what the CLI means; "
        f"declared {declared}"
    )


def test_the_flag_a_mapping_names_takes_the_value_the_input_describes() -> None:
    """A swap between two real flags would pass every other test here. Check that
    each declared flag's argparse metavar or type is consistent with its input."""
    actions = {opt: a for a in dnasp.build_parser()._actions for opt in a.option_strings}
    for name, flag in EXPECTED_FLAGS.items():
        action = actions[flag]
        if name in {"window_size", "step_size"}:
            assert action.type is int, f"{flag} should take an integer for {name}"
        if name == "genetic_code":
            assert action.choices, f"{flag} should be a fixed set of codes for {name}"
