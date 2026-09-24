import argparse
import hashlib
import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "omics_target_evidence_mapper.py"
SPEC = importlib.util.spec_from_file_location("omics_target_evidence_mapper", MODULE_PATH)
mapper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mapper)


def test_checksum_manifest_covers_every_advertised_output(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "output"
    monkeypatch.setattr(
        mapper,
        "parse_args",
        lambda: argparse.Namespace(
            demo=True,
            gene=None,
            disease=None,
            output=str(output_dir),
            max_papers=5,
            max_trials=5,
        ),
    )
    monkeypatch.setattr(
        mapper,
        "build_evidence",
        lambda _args: {
            "query": {"gene": "IL6R", "disease": "coronary artery disease", "demo_mode": True},
            "provenance": {"sources": ["UniProt"], "version": "0.1.0"},
            "literature": [],
            "trials": [],
            "target_summary": {"status": "no_result"},
            "disease_association": {"status": "skipped"},
            "limitations": [],
        },
    )
    monkeypatch.setattr(mapper, "build_report", lambda _evidence: "# Report\n")

    def write_ro_crate(output_dir, **_kwargs):
        path = Path(output_dir) / "ro-crate-metadata.json"
        path.write_text("{}\n", encoding="utf-8")
        return path

    monkeypatch.setattr(mapper, "write_ro_crate", write_ro_crate)

    mapper.main()

    manifest = output_dir / "reproducibility" / "checksums.sha256"
    entries = dict(line.split("  ", 1) for line in manifest.read_text().splitlines())
    expected = {
        "evidence.json",
        "report.md",
        "metadata.json",
        "reproducibility/environment.yml",
        "ro-crate-metadata.json",
    }
    assert set(entries.values()) == expected
    for relative_path in expected:
        content = (output_dir / relative_path).read_bytes()
        assert entries[hashlib.sha256(content).hexdigest()] == relative_path
    assert "reproducibility/checksums.sha256" not in entries.values()
