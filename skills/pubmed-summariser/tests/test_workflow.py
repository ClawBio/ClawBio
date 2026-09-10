"""Functional contract: real CLI, XML parser and output writers, fake services."""
import hashlib
import io
import json
import shlex
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))
import pubmed_summariser as app
from clawbio.providers import GenerationResult, ProviderError
from clawbio.common.report import DISCLAIMER

BACKGROUND = "Opening sentence. " + "Background detail. " * 25
RESULTS = "The BRCA1 signal was uncertain & requires replication."
XML = f"""<PubmedArticleSet>
<PubmedArticle><MedlineCitation><PMID>1001</PMID><Article>
<ArticleTitle>Synthetic &lt;script&gt; study</ArticleTitle>
<Journal><Title>Test journal</Title><JournalIssue><PubDate><Year>2025</Year></PubDate></JournalIssue></Journal>
<Abstract><AbstractText Label="Background" NlmCategory="BACKGROUND">{BACKGROUND}</AbstractText>
<AbstractText Label="Findings" NlmCategory="RESULTS">The <i>BRCA1</i> signal was uncertain &amp; requires replication.</AbstractText></Abstract>
</Article></MedlineCitation></PubmedArticle>
<PubmedArticle><MedlineCitation><PMID>1002</PMID><Article>
<ArticleTitle>Second synthetic study</ArticleTitle>
<Abstract><AbstractText>Another opening. Another finding.</AbstractText></Abstract>
</Article></MedlineCitation></PubmedArticle>
<PubmedArticle><MedlineCitation><PMID>1003</PMID><Article>
<ArticleTitle>No abstract study</ArticleTitle>
</Article></MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


@pytest.fixture
def services(monkeypatch):
    search = Mock()
    search.json.return_value = {"esearchresult": {"idlist": ["1003", "1002", "1001"]}}
    fetch = Mock(content=XML.encode())
    get = Mock(side_effect=[search, fetch])
    monkeypatch.setattr(app.pubmed_api.requests, "get", get)
    return get


def run(tmp_path, *options):
    app.main(["--query", "BRCA1", "--output", str(tmp_path), *options])
    return json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))


def test_search_preserves_existing_query_limit_and_order(services, tmp_path, capsys):
    app.main(["--query", "BRCA1", "--max-results", "20", "--output", str(tmp_path)])
    params = services.call_args_list[0].kwargs["params"]
    assert params["term"] == "BRCA1 AND english[la]"
    assert params["retmax"] == 20
    assert params["sort"] == "date"
    output = capsys.readouterr().out
    assert output.index("Synthetic <script> study") < output.index("Second synthetic study")


def test_default_preserves_source_and_saves_excerpt(services, tmp_path, capsys):
    result = run(tmp_path)
    paper = result["data"]["papers"][0]
    assert paper["abstract"] == BACKGROUND.strip() + "\n\n" + RESULTS
    assert paper["abstract_sections"] == [
        {"label": "Background", "category": "BACKGROUND", "text": BACKGROUND.strip()},
        {"label": "Findings", "category": "RESULTS", "text": RESULTS},
    ]
    assert paper["summary"]["text"] == "Opening sentence."
    assert paper["summary"]["method"] == "first-sentence"
    assert paper["summary"]["provider"] is None
    assert result["data"]["summary_config"]["method"] == "first-sentence"
    assert result["data"]["retrieved_at"]
    assert result["data"]["papers"][2]["abstract_sections"] == []
    assert result["data"]["papers"][2]["summary"]["method"] == "unavailable"
    terminal = capsys.readouterr().out
    assert "Summary method: first-sentence" in terminal
    assert "Abstract opening — no AI" in terminal
    assert BACKGROUND.strip() not in terminal
    html = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<details>" in html and "BRCA1 signal" in html
    assert "Abstract opening — no AI" in html
    assert "&lt;script&gt;" in html and "<script> study" not in html
    for file in ("report.md", "report.html"):
        assert DISCLAIMER in (tmp_path / file).read_text(encoding="utf-8")
    assert DISCLAIMER in terminal


@pytest.mark.parametrize("failure", [False, True])
def test_llm_summary_and_per_paper_fallback(services, tmp_path, monkeypatch, failure):
    provider = Mock(name="provider")
    provider.name = "ollama"
    provider.model = "test-model"
    provider.base_url = "http://localhost:11434/v1"
    generated = GenerationResult("The finding needs replication.", "ollama", "test-model", "stop", 200, 12)
    provider.generate.side_effect = [ProviderError("Request failed.") if failure else generated, generated]
    factory = Mock(return_value=provider)
    monkeypatch.setattr(app, "create_provider", factory, raising=False)
    result = run(tmp_path, "--summary-method", "llm", "--provider", "ollama", "--model", "test-model",
                 "--model-params", '{"reasoning_effort":"none"}', "--summary-max-tokens", "200")
    assert factory.call_count == 1
    assert provider.generate.call_count == 2  # skip missing abstract
    call = provider.generate.call_args_list[0]
    assert RESULTS in call.args[0]
    assert BACKGROUND.strip() in call.args[0]
    assert call.kwargs["model_params"] == {"reasoning_effort": "none"}
    assert call.kwargs["max_output_tokens"] == 200
    assert "1-2 short sentences" in call.kwargs["system"]
    assert "Attribute findings or interpretations to the paper or its authors" in call.kwargs["system"]
    assert "Preserve the strength of the evidence" in call.kwargs["system"]
    papers = result["data"]["papers"]
    assert papers[0]["abstract"].endswith(RESULTS)
    first = papers[0]["summary"]
    assert first["method"] == ("first-sentence" if failure else "llm")
    assert first["fallback_reason"] == ("Request failed." if failure else None)
    assert papers[1]["summary"]["text"] == generated.text
    assert papers[1]["summary"]["provider"] == "ollama"
    assert papers[1]["summary"]["prompt_version"] == "3"
    assert result["data"]["summary_config"]["prompt_version"] == "3"
    assert papers[1]["summary"]["usage"]["completion_tokens"] == 12
    assert bool(result["data"]["warnings"]) == failure
    provider.close.assert_called_once()
    assert "AI summary — ollama / test-model" in (tmp_path / "report.html").read_text(encoding="utf-8")


def test_explicit_first_sentence_method(services, tmp_path):
    result = run(tmp_path, "--summary-method", "first-sentence")
    assert result["data"]["papers"][0]["summary"]["method"] == "first-sentence"


def test_help_explains_default_without_ai(capsys):
    with pytest.raises(SystemExit) as exc:
        app.main(["--help"])
    assert exc.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "Default: first-sentence" in help_text
    assert "no AI model is used" in help_text


@pytest.mark.parametrize("options", [
    ["--summary-method", "llm"],
    ["--provider", "ollama"],
    ["--model", "unused"],
    ["--max-results", "0"],
    ["--max-results", "-1"],
    ["--summary-method", "llm", "--provider", "ollama", "--model", "x", "--model-params", "[]"],
    ["--summary-method", "llm", "--provider", "ollama", "--model", "x", "--model-params", '{"api_key":"secret"}'],
    ["--summary-method", "llm", "--provider", "ollama", "--model", "x", "--llm-timeout", "nan"],
    ["--summary-method", "llm", "--provider", "ollama", "--model", "x", "--summary-max-tokens", "0"],
])
def test_invalid_configuration_stops_before_search(services, tmp_path, options):
    with pytest.raises(SystemExit) as exc:
        run(tmp_path, *options)
    assert exc.value.code == 2
    services.assert_not_called()
    assert not (tmp_path / "result.json").exists()


def test_empty_search_has_valid_outputs(services, tmp_path):
    services.side_effect = [Mock(json=lambda: {"esearchresult": {"idlist": []}})]
    result = run(tmp_path)
    assert result["data"]["papers"] == []
    assert result["summary"]["paper_count"] == 0
    assert services.call_count == 1
    assert (tmp_path / "report.html").exists()
    assert (tmp_path / "report.md").exists()


def test_reproducibility_bundle_and_overwrite_warning(services, tmp_path, capsys):
    (tmp_path / "report.html").write_text("existing report")
    run(tmp_path)
    assert "overwrit" in capsys.readouterr().err.lower()
    repro = tmp_path / "reproducibility"
    assert "--query BRCA1" in (repro / "commands.sh").read_text()
    assert "requests" in (repro / "environment.yml").read_text()
    lines = (repro / "checksums.sha256").read_text().splitlines()
    assert len(lines) >= 3
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute()
        assert hashlib.sha256((tmp_path / label).read_bytes()).hexdigest() == digest


def test_fetch_failure_does_not_write_success_report(services, tmp_path):
    services.side_effect = requests.ConnectionError("offline")
    with pytest.raises(SystemExit) as exc:
        run(tmp_path)
    assert exc.value.code == 1
    assert not (tmp_path / "result.json").exists()


def test_query_file_and_replay_command(services, tmp_path):
    query_file = tmp_path / "query.txt"
    query_file.write_text("type 2 diabetes\n", encoding="utf-8")
    app.main(["--input", str(query_file), "--output", str(tmp_path)])
    assert services.call_args_list[0].kwargs["params"]["term"] == "type 2 diabetes AND english[la]"
    command = (tmp_path / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    args = shlex.split(command)
    assert args[args.index("--query") + 1] == "type 2 diabetes"
    assert Path(args[args.index("--output") + 1]).is_absolute()


def test_blank_query_stops_before_search(services, tmp_path):
    with pytest.raises(SystemExit) as exc:
        app.main(["--query", "   ", "--output", str(tmp_path)])
    assert exc.value.code == 2
    services.assert_not_called()


def test_missing_llm_model_stops_before_search(services, tmp_path, monkeypatch):
    monkeypatch.delenv("OLLAMA_MODEL", raising=False)
    monkeypatch.delenv("CLAWBIO_MODEL", raising=False)
    with pytest.raises(SystemExit) as exc:
        run(tmp_path, "--summary-method", "llm", "--provider", "ollama")
    assert exc.value.code == 2
    services.assert_not_called()


def test_replay_records_resolved_model(services, tmp_path, monkeypatch):
    provider = Mock()
    provider.model = "resolved-model"
    provider.base_url = "http://localhost:11434/v1"
    provider.generate.return_value = GenerationResult("Short summary.", "ollama", "resolved-model", "stop")
    monkeypatch.setattr(app, "create_provider", Mock(return_value=provider))
    run(tmp_path, "--summary-method", "llm", "--provider", "ollama")
    command = (tmp_path / "reproducibility" / "commands.sh").read_text().splitlines()[-1]
    args = shlex.split(command)
    assert args[args.index("--model") + 1] == "resolved-model"


def test_provider_closed_when_pubmed_fails(services, tmp_path, monkeypatch):
    services.side_effect = requests.ConnectionError("offline")
    provider = Mock()
    monkeypatch.setattr(app, "create_provider", Mock(return_value=provider))
    with pytest.raises(SystemExit):
        run(tmp_path, "--summary-method", "llm", "--provider", "ollama", "--model", "test")
    provider.close.assert_called_once()


def test_unicode_papers_survive_legacy_windows_stdout(services, tmp_path, monkeypatch):
    search = Mock(json=lambda: {"esearchresult": {"idlist": ["1001", "1002", "1003"]}})
    services.side_effect = [search, Mock(content=XML.replace("Synthetic", "Synthetic α").encode())]
    output = io.TextIOWrapper(io.BytesIO(), encoding="cp1252", errors="strict")
    monkeypatch.setattr(sys, "stdout", output)
    try:
        result = run(tmp_path)
        assert "α" in result["data"]["papers"][0]["title"]
    finally:
        output.close()
