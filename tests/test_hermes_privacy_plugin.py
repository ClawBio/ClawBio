"""Tests for the Hermes Agent privacy plugin in integrations/hermes/clawbio-privacy.

The plugin keeps genomic file contents out of the model's context. It registers
four hooks and one middleware; each is tested here with plain values, no Hermes install needed.
"""

import gzip
import importlib.util
import json
import struct
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

_PLUGIN = Path(__file__).resolve().parents[1] / "integrations" / "hermes" / "clawbio-privacy" / "__init__.py"
_spec = importlib.util.spec_from_file_location("clawbio_hermes_privacy", _PLUGIN)
plugin = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plugin)

CANARY = "CANARY7F3A"
VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n"
    f"chr22\t42127941\trs3892097\tG\tA\t50\tPASS\tMARK={CANARY}\tGT\t0/1\n"
)
TWENTY_THREE = (
    "# rsid\tchromosome\tposition\tgenotype\n"
    "rs4477212\t1\t82154\tAA\n"
    "rs3094315\t1\t752566\tAG\n"
    f"rs3131972\t1\t752721\t{CANARY}\n"
)


@dataclass
class Event:
    """The MessageEvent fields the hook reads."""
    text: str
    media_urls: list = field(default_factory=list)


def _upload(tmp_path, name, content):
    path = tmp_path / f"doc_0123456789ab_{name}"
    path.write_text(content)
    return path


def _inlined(path, name, content, caption=""):
    # Mirrors Hermes's Telegram adapter: injection, blank line, then the caption.
    injection = f"[Content of {name}]:\n{content}"
    return Event(text=f"{injection}\n\n{caption}" if caption else injection, media_urls=[str(path)])


# --- pre_gateway_dispatch: pasted file contents ---------------------------- #

def test_inlined_vcf_is_withheld_and_caption_kept(tmp_path):
    path = _upload(tmp_path, "sample.vcf", VCF)
    result = plugin.on_pre_gateway_dispatch(event=_inlined(path, "sample.vcf", VCF, "What's in this file?"))
    assert result["action"] == "rewrite"
    assert CANARY not in result["text"]
    assert "What's in this file?" in result["text"]
    assert str(path) in result["text"]
    assert "--sealed" in result["text"]
    assert "attached automatically" in result["text"]


def test_inlined_23andme_txt_is_withheld(tmp_path):
    path = _upload(tmp_path, "genome.txt", TWENTY_THREE)
    result = plugin.on_pre_gateway_dispatch(event=_inlined(path, "genome.txt", TWENTY_THREE))
    assert CANARY not in result["text"]


def test_non_genomic_upload_is_untouched(tmp_path):
    notes = "shopping list\nmilk\neggs\n"
    path = _upload(tmp_path, "notes.txt", notes)
    assert plugin.on_pre_gateway_dispatch(event=_inlined(path, "notes.txt", notes)) is None


def test_message_without_media_is_untouched():
    assert plugin.on_pre_gateway_dispatch(event=Event(text="hello")) is None


# Hypothesis reruns a test per example and warns about per-test fixtures; tmp_path and monkeypatch are
# safe to share here because every example writes the same file and sets the same variable.
_SHARED_FIXTURES = settings(suppress_health_check=[HealthCheck.function_scoped_fixture])


@_SHARED_FIXTURES
@given(caption=st.text())
def test_canary_never_survives_any_caption(tmp_path, caption):
    path = _upload(tmp_path, "sample.vcf", VCF)
    result = plugin.on_pre_gateway_dispatch(event=_inlined(path, "sample.vcf", VCF, caption))
    assert CANARY not in result["text"] or CANARY in caption


def test_dispatch_fails_closed_on_error(tmp_path, monkeypatch):
    path = _upload(tmp_path, "genome.txt", TWENTY_THREE)

    def boom(*_a, **_k):
        raise RuntimeError("bug")

    monkeypatch.setattr(plugin._policy(), "is_genomic_text", boom)
    result = plugin.on_pre_gateway_dispatch(event=_inlined(path, "genome.txt", TWENTY_THREE))
    assert result["action"] == "rewrite"
    assert CANARY not in result["text"]


# --- pre_tool_call: reading genomic files ---------------------------------- #

def _blocked(tool_name, **args):
    result = plugin.on_pre_tool_call(tool_name=tool_name, args=args)
    return isinstance(result, dict) and result.get("action") == "block"


def test_read_file_on_vcf_is_blocked(tmp_path):
    assert _blocked("read_file", path=str(_upload(tmp_path, "s.vcf", VCF)))


def test_read_file_on_23andme_txt_is_blocked_by_content(tmp_path):
    assert _blocked("read_file", path=str(_upload(tmp_path, "genome.txt", TWENTY_THREE)))


def test_read_file_on_ordinary_file_is_allowed(tmp_path):
    assert not _blocked("read_file", path=str(_upload(tmp_path, "notes.md", "# notes\n")))


def test_terminal_cat_of_vcf_is_blocked():
    assert _blocked("terminal", command="cat /tmp/doc_x_sample.vcf")


def test_terminal_sealed_clawbio_run_is_allowed():
    assert not _blocked("terminal", command="python clawbio.py run pharmgx --input /tmp/doc_x_sample.vcf --sealed")


def test_terminal_sealed_run_cannot_be_chained():
    assert _blocked("terminal", command="python clawbio.py run pharmgx --input /tmp/s.vcf --sealed; cat /tmp/s.vcf")


def test_unsealed_clawbio_run_is_blocked():
    assert _blocked("terminal", command="python clawbio.py run pharmgx --input /tmp/s.vcf")


def test_execute_code_touching_vcf_is_blocked():
    assert _blocked("execute_code", code="print(open('/tmp/s.vcf').read())")


def test_unrelated_tools_are_ignored():
    assert plugin.on_pre_tool_call(tool_name="web_search", args={"query": "rs3892097 CYP2D6"}) is None


# --- pre_tool_call: connector tools (Hermes 0.22) -------------------------- #
# A connector tool reads any local path it is given and uploads the bytes to the connector's
# storage (tools/connectors/gateway/files.py). Uploads in cache/documents are sendable.

def test_connector_sending_a_vcf_is_blocked(tmp_path):
    vcf = _upload(tmp_path, "sample.vcf", VCF)
    assert _blocked("connectors__gmail__send_email", to="me@example.org", attachment=str(vcf))


def test_connector_sending_raw_data_found_by_content_is_blocked(tmp_path):
    raw = _upload(tmp_path, "genome.txt", TWENTY_THREE)
    assert _blocked("connectors__googledrive__upload_file", file={"path": str(raw), "name": "notes"})


def test_connector_batch_with_a_genomic_file_is_blocked(tmp_path):
    vcf = _upload(tmp_path, "sample.vcf", VCF)
    calls = [{"name": "connectors__slack__send_message", "arguments": {"text": "hi"}},
             {"name": "connectors__gmail__send_email", "arguments": {"attachments": [str(vcf)]}}]
    assert _blocked("connectors__execute", calls=calls)


def test_connector_sending_a_report_is_blocked(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert _blocked("connectors__gmail__send_email", attachment=f"{LIVE_DIR}/output/pharmgx_x/report.md")


def test_connector_text_that_looks_like_a_home_path_does_not_crash():
    # Every string argument is checked, so prose reaches Path.expanduser: "~approx" names no user.
    assert plugin.on_pre_tool_call(tool_name="connectors__gmail__send_email", args={"body": "~approx 3pm"}) is None


# Tool arguments are arbitrary JSON the model wrote: text that looks like paths, nesting, lists.
_ARG_TEXT = st.text(st.sampled_from("~/._-aZ09 \x00"), max_size=30) | st.sampled_from(["~", "~approx", "/", "x" * 5000])
_ARGS = st.recursive(_ARG_TEXT | st.integers() | st.none(),
                     lambda inner: st.lists(inner, max_size=3) | st.dictionaries(st.text(max_size=5), inner, max_size=3),
                     max_leaves=8)


@given(args=st.dictionaries(st.text(max_size=5), _ARGS, max_size=4))
@example(args={"body": "~approx 3pm"})
@example(args={"": "\x00"})  # Path() refuses a NUL with ValueError
def test_connector_guard_never_raises(args):
    result = plugin.on_pre_tool_call(tool_name="connectors__gmail__send_email", args=args)
    assert result is None or result["action"] == "block"


@given(args=st.dictionaries(st.text(max_size=5), _ARGS, max_size=4), where=st.lists(st.sampled_from(["a", "b"]), max_size=3))
def test_connector_with_a_genomic_path_anywhere_is_blocked(args, where):
    nested = "/home/u/doc_0123456789ab_sample.vcf"
    for key in reversed(where):
        nested = {key: [nested]}
    assert _blocked("connectors__execute", **args, calls=nested)


def test_connector_call_without_genomic_files_is_allowed(tmp_path):
    notes = _upload(tmp_path, "notes.txt", "meeting at 3pm")
    assert plugin.on_pre_tool_call(tool_name="connectors__gmail__send_email",
                                   args={"to": "me@example.org", "body": "hello", "attachment": str(notes)}) is None


# --- transform_tool_result: backstop --------------------------------------- #

def test_genomic_tool_output_is_withheld():
    result = plugin.on_transform_tool_result(tool_name="terminal", args={}, result=VCF)
    assert isinstance(result, str) and CANARY not in result


def test_sealed_json_passes_through():
    sealed = '{"skill": "pharmgx", "success": true, "sealed": true, "files": ["report.md"]}'
    assert plugin.on_transform_tool_result(tool_name="terminal", args={}, result=sealed) is None


def test_web_results_are_not_filtered():
    assert plugin.on_transform_tool_result(tool_name="web_search", args={}, result=VCF) is None


# --- registration ---------------------------------------------------------- #

def test_register_wires_all_hooks_and_the_request_middleware():
    registered = {}

    class Ctx:
        def register_hook(self, name, fn):
            registered[name] = fn

        def register_middleware(self, kind, fn):
            registered[kind] = fn

    plugin.register(Ctx())
    assert set(registered) == {
        "pre_gateway_dispatch", "pre_tool_call", "transform_tool_result", "transform_llm_output", "llm_request",
    }


def test_note_uses_configured_clawbio_dir_and_is_itself_allowed(tmp_path, monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", "/opt/ClawBio")
    path = _upload(tmp_path, "sample.vcf", VCF)
    text = plugin.on_pre_gateway_dispatch(event=_inlined(path, "sample.vcf", VCF))["text"]
    command = f"/opt/ClawBio/.venv/bin/python /opt/ClawBio/clawbio.py run pharmgx --input {path} --sealed"
    assert "/opt/ClawBio/.venv/bin/python /opt/ClawBio/clawbio.py run <skill>" in text
    assert not _blocked("terminal", command=command)


# --- llm_request middleware: uploads that skip pre_gateway_dispatch ---------- #
# Live test 2026-10-11: a caption-less VCF sent while a turn was running was queued
# by Hermes's busy path (gateway/run_busy.py), which never calls pre_gateway_dispatch,
# and became the next user turn with the file pasted in (gateway/run_turn.py).

def _busy_turn_text(path, name, content, caption=""):
    note = (f"[The user sent a text document: '{name}'. Its content has been included below. "
            f"The file is also saved at: {path}]")
    return f"{note}\n\n{_inlined(path, name, content, caption).text}"


def _chat_request(user_content):
    return {"model": "m", "messages": [{"role": "system", "content": "You are helpful."},
                                       {"role": "user", "content": user_content}]}


def _middleware(request):
    result = plugin.on_llm_request(request=request, original_request=request)
    return result["request"] if result else request


def test_busy_path_upload_is_withheld_from_the_model_request(tmp_path):
    path = _upload(tmp_path, "sample.vcf", VCF)
    sent = _middleware(_chat_request(_busy_turn_text(path, "sample.vcf", VCF, "What's in this file?")))
    text = sent["messages"][1]["content"]
    assert CANARY not in text
    assert "What's in this file?" in text
    assert str(path) in text and "--sealed" in text
    assert sent["messages"][0]["content"] == "You are helpful."


def test_content_parts_are_scrubbed_too(tmp_path):
    # Anthropic messages and the Responses API carry text in content parts.
    path = _upload(tmp_path, "genome.txt", TWENTY_THREE)
    parts = [{"type": "input_text", "text": _busy_turn_text(path, "genome.txt", TWENTY_THREE)}]
    request = {"model": "m", "input": [{"role": "user", "content": parts}]}
    assert CANARY not in str(_middleware(request))


def test_pasted_genotypes_without_a_saved_file_are_withheld():
    assert CANARY not in str(_middleware(_chat_request(f"look at this\n{VCF}")))


def test_ordinary_request_is_untouched():
    assert plugin.on_llm_request(request=_chat_request("hello"), original_request=_chat_request("hello")) is None


def test_middleware_fails_closed_on_error(tmp_path, monkeypatch):
    # Hermes runs middleware fail-open: an exception would send the original request.
    path = _upload(tmp_path, "sample.vcf", VCF)

    def boom(*_a, **_k):
        raise RuntimeError("bug")

    monkeypatch.setattr(plugin._policy(), "scrub_text", boom)
    request = _chat_request(_busy_turn_text(path, "sample.vcf", VCF))
    assert CANARY not in str(_middleware(request))


@_SHARED_FIXTURES
@given(caption=st.text(), parts=st.booleans())
def test_canary_never_reaches_the_model_on_the_busy_path(tmp_path, caption, parts):
    path = _upload(tmp_path, "sample.vcf", VCF)
    text = _busy_turn_text(path, "sample.vcf", VCF, caption)
    content = [{"type": "text", "text": text}] if parts else text
    assert CANARY not in str(_middleware(_chat_request(content))) or CANARY in caption


# --- transform_llm_output: attach the report without trusting the model's formatting --- #
# Only report.md is sent: Telegram renders it in its document viewer (live test 2026-10-14).

def _sealed_run(tmp_path, report=True):
    """A sealed run's output dir and the terminal result Hermes hands the hook for it."""
    out = tmp_path / "output" / "pharmgx_20261014_111134"
    out.mkdir(parents=True)
    if report:
        (out / "report.md").write_text("# ClawBio PharmGx Report\n")
    sealed = {"skill": "pharmgx", "success": True, "sealed": True, "output_dir": str(out)}
    # Hermes wraps terminal output as {"output": "<stdout>"}, so the sealed JSON arrives escaped.
    return out, json.dumps({"output": json.dumps(sealed), "exit_code": 0})


def test_only_report_md_is_attached_and_the_models_own_output_tags_are_dropped(tmp_path):
    out, result = _sealed_run(tmp_path)
    assert plugin.on_transform_tool_result(tool_name="terminal", result=result, session_id="s1") is None
    reply = plugin.on_transform_llm_output(
        response_text=f"Ready. MEDIA:{out}/report.html\nMEDIA:{out}/result.json", session_id="s1")
    assert reply.splitlines()[-1] == f"MEDIA:{out}/report.md"
    assert reply.count("MEDIA:") == 1


def test_report_is_attached_even_when_the_model_mangles_the_tag(tmp_path):
    out, result = _sealed_run(tmp_path)
    plugin.on_transform_tool_result(tool_name="terminal", result=result, session_id="s2")
    reply = plugin.on_transform_llm_output(response_text=f"Done. ||MEDIA:{out}/report.md||", session_id="s2")
    assert reply.splitlines()[-1] == f"MEDIA:{out}/report.md"
    assert "||" not in reply


def test_report_is_attached_when_the_model_forgets_the_tag(tmp_path):
    out, result = _sealed_run(tmp_path)
    plugin.on_transform_tool_result(tool_name="terminal", result=result, session_id="s3")
    assert plugin.on_transform_llm_output(response_text="Analysis done.", session_id="s3").endswith(
        f"\nMEDIA:{out}/report.md")


def test_attachment_is_per_session_and_once(tmp_path):
    _, result = _sealed_run(tmp_path)
    plugin.on_transform_tool_result(tool_name="terminal", result=result, session_id="s4")
    assert plugin.on_transform_llm_output(response_text="hi", session_id="other") is None
    assert plugin.on_transform_llm_output(response_text="hi", session_id="s4") is not None
    assert plugin.on_transform_llm_output(response_text="hi again", session_id="s4") is None


def test_reply_without_a_run_or_output_tags_is_untouched():
    assert plugin.on_transform_llm_output(response_text="MEDIA:/tmp/chart.png", session_id="s5") is None


# --- pre_tool_call: ClawBio output (reports) ------------------------------- #
# Live test 2026-10-11: a routed profile's agent ran this grep on another profile's
# report and got the text back; report sections without rsIDs pass the content check.

LIVE_DIR = "/Users/assistant/clawbio-privacy-test"
LIVE_GREP = f"grep -A3 'Gene Profiles' {LIVE_DIR}/output/pharmgx_20261010_224608/report.md"


def test_live_grep_of_another_profiles_report_is_blocked(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert _blocked("terminal", command=LIVE_GREP)


def test_read_file_under_output_is_blocked(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert _blocked("read_file", path=f"{LIVE_DIR}/output/pharmgx_x/report.md")


def test_report_block_message_names_reports_not_only_genomic_files(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    report = f"{LIVE_DIR}/output/pharmgx_x/report.md"
    for tool, args in [("terminal", {"command": LIVE_GREP}), ("read_file", {"path": report})]:
        message = plugin.on_pre_tool_call(tool_name=tool, args=args)["message"]
        assert "ClawBio output" in message
        assert f"--input {report}" not in message


def test_read_file_through_a_symlinked_checkout_is_blocked(tmp_path, monkeypatch):
    real = tmp_path / "clawbio"
    (real / "output" / "run1").mkdir(parents=True)
    (real / "output" / "run1" / "report.md").write_text("| CYP2C19 | *1/*1 |\n")
    alias = tmp_path / "alias"
    alias.symlink_to(real)
    monkeypatch.setenv("CLAWBIO_DIR", str(real))
    assert _blocked("read_file", path=str(alias / "output" / "run1" / "report.md"))
    assert _blocked("terminal", command=f"cat {alias}/output/run1/report.md")


def test_cd_into_checkout_then_relative_output_is_blocked(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert _blocked("terminal", command=f"cd {LIVE_DIR} && cat output/pharmgx_x/report.md")


def test_clawbio_discovery_commands_stay_allowed(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    py = f"{LIVE_DIR}/.venv/bin/python {LIVE_DIR}/clawbio.py"
    assert not _blocked("terminal", command=f"{py} list")
    assert not _blocked("terminal", command=f"{py} run --help")
    assert not _blocked("terminal", command=f"{py} run pharmgx --input /tmp/s.vcf --sealed")


def test_profile_env_clawbio_dir_wins_and_both_outputs_are_blocked(tmp_path, monkeypatch):
    home = tmp_path / "profiles" / "alice"
    home.mkdir(parents=True)
    (home / ".env").write_text('# comment\nCLAWBIO_DIR="/srv/alice-clawbio"\n')
    monkeypatch.setattr(plugin, "_hermes_home", lambda: home)
    monkeypatch.setenv("CLAWBIO_DIR", "/srv/shared-clawbio")
    path = _upload(tmp_path, "sample.vcf", VCF)
    text = plugin.on_pre_gateway_dispatch(event=_inlined(path, "sample.vcf", VCF))["text"]
    assert "/srv/alice-clawbio/clawbio.py run <skill>" in text
    assert _blocked("read_file", path="/srv/alice-clawbio/output/r/report.md")
    assert _blocked("read_file", path="/srv/shared-clawbio/output/r/report.md")


_segment = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_.-", min_size=1, max_size=12).filter(
    lambda s: s not in (".", ".."))


@_SHARED_FIXTURES
@given(parts=st.lists(_segment, min_size=1, max_size=4))
@example(parts=["pharmgx_20261010_224608", "report.md"])
def test_anything_under_output_is_blocked_for_every_reader(monkeypatch, parts):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    target = "/".join([LIVE_DIR, "output", *parts])
    assert _blocked("read_file", path=target)
    assert _blocked("search_files", path=target)
    assert _blocked("terminal", command=f"head -50 {target}")
    assert _blocked("execute_code", code=f"print(open('{target}').read())")


# --- Delivery: only the report, whatever the model asks for -------------- #

@_SHARED_FIXTURES
@given(
    run=_segment,
    name=st.sampled_from(["report.html", "report.md", "result.json", "figures/plot.png"]),
    wrap=st.sampled_from(["{}", "||{}||", "`{}`", "**{}**"]),
)
# Live Discord test 2026-10-14, turn 2: no run this turn, the model re-attached an earlier run's output.
@example(run="pharmgx_20261014_223854", name="report.html", wrap="{}")
def test_no_tag_for_clawbio_output_survives_any_turn(monkeypatch, run, name, wrap):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    tag = wrap.format(f"MEDIA:{LIVE_DIR}/output/{run}/{name}")
    reply = plugin.on_transform_llm_output(response_text=f"Here you go:\n{tag}", session_id="no-run")
    assert reply is not None and f"{LIVE_DIR}/output" not in reply


# --- Found by a Hypothesis review, 2026-10-14 -------------------------------- #


def test_unsealed_run_on_a_genomic_file_is_blocked_despite_sealed_in_the_path():
    assert _blocked("terminal", command="python clawbio.py run pharmgx --input /tmp/x--sealed.vcf")


def test_bom_prefixed_txt_upload_is_withheld(tmp_path):
    body = "\ufeff##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\tPASS\tM=" + CANARY + "\tGT\t0/1\n"
    path = _upload(tmp_path, "sample.txt", body)
    result = plugin.on_pre_gateway_dispatch(event=_inlined(path, "sample.txt", body, "hi"))
    assert result and CANARY not in result["text"]


# --- Found by the live Hermes test on the Mac mini, 2026-10-15 ---------------- #

@pytest.mark.parametrize("command", [
    # Hermes keeps the terminal's folder between calls, so `cd $CLAWBIO_DIR` then a relative path.
    "sed -n 54,80p output/pharmgx_20261015_093246/report.md",
    "cat ./output/pharmgx_20261015_093246/report.md",
    # `cd output` first, then the run folder alone.
    "grep -A12 'Gene Profiles' pharmgx_20261015_093246/report.md",
])
def test_relative_clawbio_output_is_blocked_with_or_without_clawbio_dir(command, monkeypatch):
    monkeypatch.delenv("CLAWBIO_DIR", raising=False)
    assert _blocked("terminal", command=command)
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert _blocked("terminal", command=command)


def test_relative_read_file_under_output_is_blocked(monkeypatch):
    monkeypatch.delenv("CLAWBIO_DIR", raising=False)
    assert _blocked("read_file", path="output/pharmgx_20261015_093246/report.md")


def test_listing_output_and_skills_stays_allowed(monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", LIVE_DIR)
    assert not _blocked("terminal", command="ls -t output | head -1")
    assert not _blocked("terminal", command=f"{LIVE_DIR}/.venv/bin/python {LIVE_DIR}/clawbio.py list")


def test_report_with_a_skill_prefix_is_attached(tmp_path):
    # nutrigx writes nutrigx_report.md; live Telegram test delivered nothing for it.
    out = tmp_path / "output" / "nutrigx_20261015_094022"
    out.mkdir(parents=True)
    (out / "nutrigx_report.md").write_text("# NutriGx\n")
    sealed = {"skill": "nutrigx", "success": True, "sealed": True, "output_dir": str(out)}
    plugin.on_transform_tool_result(tool_name="terminal", result=json.dumps({"output": json.dumps(sealed)}),
                                    session_id="nx")
    reply = plugin.on_transform_llm_output(response_text="Done.", session_id="nx")
    assert reply.splitlines()[-1] == f"MEDIA:{out}/nutrigx_report.md"


def test_manifest_declares_the_request_middleware():
    # `hermes plugins validate` fails a plugin that registers middleware it does not declare.
    import yaml

    manifest = yaml.safe_load(_PLUGIN.with_name("plugin.yaml").read_text())
    assert manifest["provides_middleware"] == ["llm_request"]


@pytest.mark.parametrize("body, kind", [(VCF, "VCF"), (TWENTY_THREE, "raw genotype")])
def test_upload_note_says_what_the_file_is_and_where_to_find_skills(tmp_path, body, kind):
    # Live test: the model saw only `run <skill>` for a .txt upload and guessed bio-orchestrator.
    path = _upload(tmp_path, "genome.txt", body)
    text = plugin.on_pre_gateway_dispatch(event=_inlined(path, "genome.txt", body))["text"]
    assert kind in text
    assert "clawbio.py list" in text


# --- Found by the file-format matrix on the Mac mini, 2026-10-15 ------------- #
# Synthetic files read through Hermes's real tool path (model_tools.handle_function_call).

def _hermes_wrapped(tool, text):
    # What transform_tool_result actually receives: terminal output and read_file content as JSON.
    if tool == "read_file":
        numbered = "\n".join(f"{i}|{line}" for i, line in enumerate(text.splitlines(), 1))
        return json.dumps({"content": numbered, "total_lines": len(text.splitlines())})
    return json.dumps({"output": text, "exit_code": 0, "error": None})


@given(body=st.sampled_from([VCF, TWENTY_THREE]),
       tool=st.sampled_from(["terminal", "read_file"]), prefix=st.text(max_size=20))
@example(body=VCF, tool="terminal", prefix="")  # live: full `bcftools view` output passed through
def test_genomic_output_is_withheld_inside_hermes_json(body, tool, prefix):
    result = plugin.on_transform_tool_result(tool_name=tool, result=_hermes_wrapped(tool, prefix + "\n" + body))
    assert result == plugin.WITHHELD


def test_glob_naming_a_bam_is_blocked(tmp_path):
    (tmp_path / "sample.bam").write_bytes(b"\x1f\x8b\x08\x04")
    assert _blocked("terminal", command=f"gzip -dc < {tmp_path}/sample.b?m | strings | head -40")
    assert _blocked("terminal", command=f"for f in {tmp_path}/sample.b?m; do samtools view $f; done")


def _bcf_bytes():
    # BCF: magic, header length, then the VCF header text, all BGZF (gzip) compressed.
    header = VCF.encode() + b"\x00"
    return gzip.compress(b"BCF\x02\x02" + struct.pack("<I", len(header)) + header)


def test_compressed_genomic_files_are_found_by_content(tmp_path):
    gz = tmp_path / "genome.txt.gz"
    gz.write_bytes(gzip.compress(TWENTY_THREE.encode()))
    zipped = tmp_path / "genome.zip"
    with zipfile.ZipFile(zipped, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("README.txt", "23andMe raw data export\n")
        z.writestr("genome_Test_User_v5_Full.txt", TWENTY_THREE)
    bcf = tmp_path / "sample.bcf"
    bcf.write_bytes(_bcf_bytes())
    for path in (gz, zipped, bcf):
        assert _blocked("terminal", command=f"gzip -dc < {path} | cut -f4")


def test_sam_is_found_by_its_header(tmp_path):
    sam = tmp_path / "sample.sam"
    sam.write_text(f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\nread_{CANARY}\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n")
    assert _blocked("read_file", path=str(sam))
    assert _blocked("terminal", command=f"cat {sam}")


@pytest.mark.parametrize("photo", ["/Users/u/.hermes/cache/images/img_0123456789ab_IMG_20261015_093000.jpg",
                                   "/Users/u/Pictures/Screenshot_20261015_093000.png"])
def test_phone_photos_with_timestamps_are_not_clawbio_runs(photo):
    # Android names photos IMG_YYYYMMDD_HHMMSS; only ClawBio's lowercase <skill>_<timestamp> folders are runs.
    assert not _blocked("read_file", path=photo)
    assert not _blocked("terminal", command=f"ls -la {photo}")
    assert plugin.on_transform_llm_output(response_text=f"Here it is\nMEDIA:{photo}", session_id="photo") is None


# --- Sealed runs, run folders and archives ---------------------------------- #

_SEALED = "python clawbio.py run pharmgx --input /tmp/s.vcf --sealed"


@pytest.mark.parametrize("skill", ["pharmgx", "nutrigx", "full_profile", "genome-compare", "a"])
def test_run_folders_are_blocked_and_photos_with_the_same_stamp_are_not(skill):
    stamp = "20261015_093246"
    assert _blocked("terminal", command=f"cat {skill}_{stamp}/report.md")
    assert not _blocked("terminal", command=f"ls -la /photos/IMG_{stamp}.jpg")


# --- Finding the policy: a thin adapter over clawbio/privacy.py -------------- #

def _copied_plugin(tmp_path):
    # A copied install (~/.hermes/plugins/clawbio-privacy/) cannot reach the checkout through its own path.
    copy = tmp_path / "plugins" / "clawbio-privacy" / "__init__.py"
    copy.parent.mkdir(parents=True)
    copy.write_text(_PLUGIN.read_text())
    spec = importlib.util.spec_from_file_location("clawbio_hermes_privacy_copy", copy)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_copied_plugin_without_clawbio_dir_fails_closed(tmp_path, monkeypatch):
    monkeypatch.delenv("CLAWBIO_DIR", raising=False)
    copied = _copied_plugin(tmp_path)
    notes = tmp_path / "notes.txt"
    notes.write_text("meeting notes\n")
    blocked = copied.on_pre_tool_call(tool_name="read_file", args={"path": str(notes)})
    assert blocked["action"] == "block" and "CLAWBIO_DIR" in blocked["message"]
    assert copied.on_transform_tool_result(tool_name="terminal", result="hello") == copied.WITHHELD


def test_copied_plugin_finds_the_policy_through_clawbio_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("CLAWBIO_DIR", str(_PLUGIN.parents[3]))
    copied = _copied_plugin(tmp_path)
    notes = tmp_path / "notes.txt"
    notes.write_text("meeting notes\n")
    assert copied.on_pre_tool_call(tool_name="read_file", args={"path": str(notes)}) is None
    assert copied.on_pre_tool_call(tool_name="read_file", args={"path": str(tmp_path / "s.vcf")})["action"] == "block"


def test_withheld_text_matches_the_policy():
    # The adapter keeps its own copy for when the policy is missing; the two must not drift.
    assert plugin.WITHHELD == plugin._policy().WITHHELD


@pytest.mark.parametrize("missing", ["check_path", "check_command"])
def test_an_older_policy_without_the_checks_still_blocks(monkeypatch, missing):
    # A checkout older than the adapter loads, but lacks the check; the hook must block, not let the tool run.
    monkeypatch.delattr(plugin._policy(), missing)
    tool, args = ("read_file", {"path": "/etc/hosts"}) if missing == "check_path" else ("terminal", {"command": "ls"})
    assert _blocked(tool, **args)
