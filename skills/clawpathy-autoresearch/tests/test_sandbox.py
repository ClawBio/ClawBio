from pathlib import Path
import pytest
from skills.clawpathy_autoresearch.sandbox import Sandbox, load_sandbox


def test_load_sandbox_full(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text(
        "allowed_tools: [Read, Bash]\n"
        "network:\n  allowed: false\n"
        "data_dir: ./data\n"
        "python_packages: [pandas]\n"
        "timeout_seconds: 300\n"
    )
    sb = load_sandbox(path)
    assert sb.allowed_tools == ["Read", "Bash"]
    assert sb.network_allowed is False
    assert sb.network_hosts == []
    assert sb.data_dir == "./data"
    assert sb.python_packages == ["pandas"]
    assert sb.timeout_seconds == 300


def test_load_sandbox_network_hostlist(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text(
        "allowed_tools: [WebFetch]\n"
        "network:\n  allowed: [ebi.ac.uk, ensembl.org]\n"
    )
    sb = load_sandbox(path)
    assert sb.network_allowed is True
    assert sb.network_hosts == ["ebi.ac.uk", "ensembl.org"]


def test_load_sandbox_defaults(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text("allowed_tools: [Read]\n")
    sb = load_sandbox(path)
    assert sb.network_allowed is False
    assert sb.timeout_seconds == 600
    assert sb.python_packages == []


def test_missing_allowed_tools_raises(tmp_path: Path):
    path = tmp_path / "sandbox.yaml"
    path.write_text("network:\n  allowed: false\n")
    with pytest.raises(ValueError, match="allowed_tools"):
        load_sandbox(path)


def test_sandbox_as_contract_text():
    sb = Sandbox(
        allowed_tools=["Read"],
        network_allowed=False,
        network_hosts=[],
        data_dir="./data",
        python_packages=["pandas"],
        timeout_seconds=300,
    )
    text = sb.as_contract_text()
    assert "allowed_tools" in text
    assert "Read" in text
    assert "network" in text
