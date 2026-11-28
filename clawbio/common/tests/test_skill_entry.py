"""Which file runs a skill, and which standard flags it takes, read from the skill folder."""

from hypothesis import given, settings, strategies as st

from clawbio.common.skill_entry import STANDARD_RUN_FLAGS, run_flags


def _script(tmp_path, body, name="s.py"):
    path = tmp_path / name
    path.write_text(body)
    return path


def test_lists_the_standard_flags_the_script_defines(tmp_path):
    script = _script(tmp_path, (
        "import argparse\n"
        "p = argparse.ArgumentParser()\n"
        "p.add_argument('--output')\n"
        "p.add_argument('--demo', action='store_true')\n"
        "p.add_argument('--genes')\n"
    ))
    assert run_flags(script) == ["--demo", "--output"]


def test_flags_added_to_an_argument_group_count(tmp_path):
    script = _script(tmp_path, (
        "import argparse\n"
        "p = argparse.ArgumentParser()\n"
        "g = p.add_argument_group('io')\n"
        "g.add_argument('-i', '--input')\n"
    ))
    assert run_flags(script) == ["--input"]


def test_a_non_python_script_has_no_run_flags(tmp_path):
    assert run_flags(_script(tmp_path, "# --input --output --demo\n", "s.R")) == []


def test_an_unparseable_script_has_no_run_flags(tmp_path):
    assert run_flags(_script(tmp_path, "def (:\n")) == []


def test_no_script_has_no_run_flags():
    assert run_flags(None) == []


_OTHER = ["--genes", "--top", "--inputs", "--output-format", "--demo-size"]


@settings(deadline=None, max_examples=50)
@given(
    flags=st.lists(st.sampled_from([*STANDARD_RUN_FLAGS, *_OTHER]), unique=True),
    grouped=st.booleans(),
)
def test_run_flags_is_exactly_the_standard_subset_in_any_order(tmp_path_factory, flags, grouped):
    # Order of add_argument calls and look-alike flags (--inputs, --output-format)
    # never change which standard flags are reported.
    target = "g" if grouped else "p"
    lines = ["import argparse", "p = argparse.ArgumentParser()", "g = p.add_argument_group('x')"]
    lines += [f"{target}.add_argument({flag!r})" for flag in flags]
    script = tmp_path_factory.mktemp("s") / "s.py"
    script.write_text("\n".join(lines) + "\n")
    assert run_flags(script) == sorted(f for f in flags if f in STANDARD_RUN_FLAGS)
