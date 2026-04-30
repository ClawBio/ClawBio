"""Tests for the struct-predictor skill."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

# Add skill dir to path
SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from struct_predictor_core.io import validate_and_prepare, VALID_AA

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

VALID_YAML = """\
sequences:
  - id: proteinA
    entity_type: protein
    sequence: ACDEFGHIKLM
  - id: proteinB
    entity_type: protein
    sequence: NPQRSTVWY
"""

# Native Boltz format (nested protein: key)
NATIVE_BOLTZ_YAML = """\
version: 1
sequences:
  - protein:
      id: A
      sequence: NLYIQWLKDGGPSSGRPPPS
"""

# Native Boltz format with multiple chain ids
NATIVE_BOLTZ_MULTIMER_YAML = """\
version: 1
sequences:
  - protein:
      id: A
      sequence: MAHHHHHHVAVDAVS
  - protein:
      id: B
      sequence: MRYAFAAEATTCNAF
"""

BAD_YAML_MISSING_SEQ = """\
sequences:
  - id: proteinA
    entity_type: protein
"""

BAD_YAML_INVALID_AA = """\
sequences:
  - id: Bad
    entity_type: protein
    sequence: ACDEFGHIKLMBXZ
"""

EMPTY_YAML = """\
sequences:
  - id: EmptyProt
    entity_type: protein
    sequence: ""
"""

# ---------------------------------------------------------------------------
# TestValidateAndPrepare
# ---------------------------------------------------------------------------

class TestValidateAndPrepare:
    def test_fasta_raises(self, tmp_path):
        f = tmp_path / "prot.fasta"
        f.write_text(">HP35\nLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF\n")
        with pytest.raises(ValueError, match="no longer supported"):
            validate_and_prepare(f, tmp_path / "work")

    def test_yaml_valid(self, tmp_path):
        f = tmp_path / "complex.yaml"
        f.write_text(VALID_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        assert result["input_type"] == "yaml"
        assert len(result["sequences"]) == 2
        assert result["boltz_input_path"].exists()

    def test_yaml_missing_sequence_raises(self, tmp_path):
        f = tmp_path / "bad.yaml"
        f.write_text(BAD_YAML_MISSING_SEQ)
        work = tmp_path / "work"
        with pytest.raises(ValueError, match="missing 'sequence'"):
            validate_and_prepare(f, work)

    def test_invalid_aa_raises(self, tmp_path):
        f = tmp_path / "bad.yaml"
        f.write_text(BAD_YAML_INVALID_AA)
        work = tmp_path / "work"
        with pytest.raises(ValueError, match="Invalid amino acid"):
            validate_and_prepare(f, work)

    def test_native_boltz_yaml_parsed(self, tmp_path):
        """Native Boltz format {protein: {id: A, sequence: ...}} must be accepted."""
        f = tmp_path / "trpcage.yaml"
        f.write_text(NATIVE_BOLTZ_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        assert result["input_type"] == "yaml"
        assert len(result["sequences"]) == 1
        assert result["sequences"][0]["name"] == "A"
        assert result["sequences"][0]["sequence"] == "NLYIQWLKDGGPSSGRPPPS"

    def test_native_boltz_yaml_msa_injected(self, tmp_path):
        """Output YAML must contain msa: empty so Boltz runs offline."""
        f = tmp_path / "trpcage.yaml"
        f.write_text(NATIVE_BOLTZ_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        output_text = result["boltz_input_path"].read_text()
        assert "msa" in output_text and "empty" in output_text

    def test_native_boltz_multimer_parsed(self, tmp_path):
        f = tmp_path / "multimer.yaml"
        f.write_text(NATIVE_BOLTZ_MULTIMER_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        assert len(result["sequences"]) == 2

    def test_work_dir_created(self, tmp_path):
        f = tmp_path / "complex.yaml"
        f.write_text(VALID_YAML)
        work = tmp_path / "new_work_dir"
        assert not work.exists()
        validate_and_prepare(f, work)
        assert work.exists()

    def test_written_yaml_contains_sequence(self, tmp_path):
        f = tmp_path / "trpcage.yaml"
        f.write_text(NATIVE_BOLTZ_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        content = result["boltz_input_path"].read_text()
        assert "NLYIQWLKDGGPSSGRPPPS" in content

    def test_empty_sequence_raises(self, tmp_path):
        f = tmp_path / "empty.yaml"
        f.write_text(EMPTY_YAML)
        work = tmp_path / "work"
        with pytest.raises(ValueError, match="is empty"):
            validate_and_prepare(f, work)

    def test_chain_ids_assigned(self, tmp_path):
        f = tmp_path / "multi.yaml"
        f.write_text(VALID_YAML)
        work = tmp_path / "work"
        result = validate_and_prepare(f, work)
        assert result["sequences"][0]["chain_id"] == "A"
        assert result["sequences"][1]["chain_id"] == "B"


from struct_predictor_core.predict import run_boltz, _find_cif, _build_boltz_cmd


class TestBuildBoltzCmd:
    def test_basic_structure(self, tmp_path):
        cmd = _build_boltz_cmd(tmp_path / "input.yaml", tmp_path / "out")
        assert cmd[0] == "boltz"
        assert cmd[1] == "predict"
        assert str(tmp_path / "input.yaml") in cmd

    def test_out_dir_flag(self, tmp_path):
        out = tmp_path / "out"
        cmd = _build_boltz_cmd(tmp_path / "input.yaml", out)
        assert "--out_dir" in cmd
        assert str(out) in cmd

    def test_no_msa_server_flag(self, tmp_path):
        """Boltz always runs offline — --use_msa_server must never appear."""
        cmd = _build_boltz_cmd(tmp_path / "input.yaml", tmp_path / "out")
        assert "--use_msa_server" not in cmd
        assert "--use-msa-server" not in cmd


class TestFindCif:
    def test_finds_cif(self, tmp_path):
        pred_dir = tmp_path / "predictions" / "Trpcage"
        pred_dir.mkdir(parents=True)
        cif = pred_dir / "Trpcage_model_0.cif"
        cif.write_text("data_Trpcage\n")
        conf = pred_dir / "confidence_Trpcage_model_0.json"
        conf.write_text('{"plddt": [95.0]}')

        result = _find_cif(tmp_path)
        assert result["cif_path"] == cif
        assert result["confidence_json_path"] == conf

    def test_raises_if_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="No CIF file found"):
            _find_cif(tmp_path)


class TestRunBoltz:
    def test_success(self, tmp_path):
        pred_dir = tmp_path / "boltz_out" / "predictions" / "Trpcage"
        pred_dir.mkdir(parents=True)
        cif = pred_dir / "Trpcage_model_0.cif"
        cif.write_text("data_Trpcage\n")
        conf = pred_dir / "confidence_Trpcage_model_0.json"
        conf.write_text('{"plddt": [95.0], "pae": [[0.5]]}')

        mock_proc = MagicMock()
        mock_proc.returncode = 0

        with patch("struct_predictor_core.predict.subprocess.run", return_value=mock_proc):
            result = run_boltz(
                input_path=tmp_path / "input.yaml",
                boltz_output_dir=tmp_path / "boltz_out",
            )

        assert result["cif_path"] == cif

    def test_nonzero_exit_raises(self, tmp_path):
        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stderr = "CUDA out of memory"

        with patch("struct_predictor_core.predict.subprocess.run", return_value=mock_proc):
            with pytest.raises(RuntimeError, match="Boltz exited with code 1"):
                run_boltz(
                    input_path=tmp_path / "input.yaml",
                    boltz_output_dir=tmp_path / "boltz_out",
                )


# ---------------------------------------------------------------------------
# TestConfidence
# ---------------------------------------------------------------------------

from struct_predictor_core.confidence import extract_confidence, _parse_plddt_from_cif, _parse_pae_from_json, _read_atom_site_columns


def _make_synthetic_cif(tmp_path: Path, n_residues: int = 5, n_chains: int = 1) -> Path:
    """Write a minimal mmCIF with CA and CB atoms."""
    lines = [
        "data_TEST", "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.label_atom_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.B_iso_or_equiv",
    ]
    atom_id = 1
    per_chain = n_residues // n_chains
    for c_idx in range(n_chains):
        chain_id = chr(ord("A") + c_idx)
        for r in range(per_chain):
            bfac = 70.0 + r * 2.0
            lines.append(f"ATOM {atom_id} CA {chain_id} {r+1} {r*1.0:.3f} 0.000 0.000 {bfac:.2f}")
            atom_id += 1
            lines.append(f"ATOM {atom_id} CB {chain_id} {r+1} {r*1.0+0.5:.3f} 0.000 0.000 99.99")
            atom_id += 1
    cif_path = tmp_path / "structure_model_0.cif"
    cif_path.write_text("\n".join(lines) + "\n")
    return cif_path


def _make_synthetic_conf_json(tmp_path: Path, n: int) -> Path:
    plddt = [80.0 + i for i in range(n)]
    pae = [[float(abs(i - j)) for j in range(n)] for i in range(n)]
    data = {"plddt": plddt, "pae": pae, "ptm": 0.85, "iptm": 0.72}
    p = tmp_path / "confidence_TEST_model_0.json"
    p.write_text(json.dumps(data))
    return p


class TestReadAtomSiteColumns:
    """Unit tests for the pure-Python mmCIF loop reader."""

    def _write(self, tmp_path, content, name="test.cif"):
        p = tmp_path / name
        p.write_text(content)
        return p

    def test_happy_path_returns_four_columns(self, tmp_path):
        cif = _make_synthetic_cif(tmp_path, n_residues=3)
        cols = _read_atom_site_columns(cif)
        assert len(cols["_atom_site.label_atom_id"]) == 6   # CA + CB x 3
        assert len(cols["_atom_site.label_asym_id"]) == 6
        assert len(cols["_atom_site.label_seq_id"]) == 6
        assert len(cols["_atom_site.B_iso_or_equiv"]) == 6

    def test_no_loop_returns_empty(self, tmp_path):
        cif = self._write(tmp_path, "data_EMPTY\n_some.field value\n")
        cols = _read_atom_site_columns(cif)
        assert cols["_atom_site.label_atom_id"] == []

    def test_wrong_loop_returns_empty(self, tmp_path):
        content = "data_WRONG\nloop_\n_other.field_a\n_other.field_b\nFOO BAR\n"
        cif = self._write(tmp_path, content)
        cols = _read_atom_site_columns(cif)
        assert cols["_atom_site.label_atom_id"] == []

    def test_comment_lines_skipped(self, tmp_path):
        content = (
            "data_TEST\n"
            "# This is a comment\n"
            "loop_\n"
            "_atom_site.label_atom_id\n"
            "_atom_site.label_asym_id\n"
            "_atom_site.label_seq_id\n"
            "_atom_site.B_iso_or_equiv\n"
            "# another comment\n"
            "CA A 1 85.0\n"
            "CB A 1 99.9\n"
        )
        cif = self._write(tmp_path, content)
        cols = _read_atom_site_columns(cif)
        assert cols["_atom_site.label_atom_id"] == ["CA", "CB"]
        assert cols["_atom_site.B_iso_or_equiv"] == ["85.0", "99.9"]

    def test_multi_loop_picks_atom_site(self, tmp_path):
        content = (
            "data_TEST\n"
            "loop_\n"
            "_chem_comp.id\n"
            "_chem_comp.name\n"
            "ALA ALANINE\n"
            "loop_\n"
            "_atom_site.label_atom_id\n"
            "_atom_site.label_asym_id\n"
            "_atom_site.label_seq_id\n"
            "_atom_site.B_iso_or_equiv\n"
            "CA A 1 92.5\n"
        )
        cif = self._write(tmp_path, content)
        cols = _read_atom_site_columns(cif)
        assert cols["_atom_site.label_atom_id"] == ["CA"]
        assert cols["_atom_site.B_iso_or_equiv"] == ["92.5"]

    def test_ca_values_match_bfactors(self, tmp_path):
        cif = _make_synthetic_cif(tmp_path, n_residues=3)
        cols = _read_atom_site_columns(cif)
        ca_indices = [i for i, a in enumerate(cols["_atom_site.label_atom_id"]) if a == "CA"]
        bfacs = [float(cols["_atom_site.B_iso_or_equiv"][i]) for i in ca_indices]
        assert bfacs == [70.0, 72.0, 74.0]


class TestConfidence:
    def test_plddt_shape(self, tmp_path):
        cif = _make_synthetic_cif(tmp_path, n_residues=5)
        plddt = _parse_plddt_from_cif(cif)
        assert plddt.shape == (5,)

    def test_plddt_only_ca(self, tmp_path):
        cif = _make_synthetic_cif(tmp_path, n_residues=3)
        plddt = _parse_plddt_from_cif(cif)
        assert not np.any(np.isclose(plddt, 99.99))

    def test_pae_shape(self, tmp_path):
        n = 5
        conf_json = _make_synthetic_conf_json(tmp_path, n)
        pae = _parse_pae_from_json(conf_json, n)
        assert pae.shape == (n, n)

    def test_pae_missing_json_returns_zeros(self, tmp_path):
        n = 4
        pae = _parse_pae_from_json(None, n)
        assert pae.shape == (n, n)
        assert np.all(pae == 0.0)

    def test_extract_confidence_single_chain(self, tmp_path):
        n = 4
        cif = _make_synthetic_cif(tmp_path, n_residues=n, n_chains=1)
        conf = _make_synthetic_conf_json(tmp_path, n)
        result = extract_confidence(cif, conf)
        assert result["plddt"].shape == (n,)
        assert result["pae"].shape == (n, n)
        assert len(result["chain_boundaries"]) == 1
        assert result["chain_boundaries"][0]["chain_id"] == "A"
        assert result["chain_boundaries"][0]["start"] == 0
        assert result["chain_boundaries"][0]["end"] == n - 1

    def test_extract_confidence_multi_chain(self, tmp_path):
        n = 4  # 2 residues per chain
        cif = _make_synthetic_cif(tmp_path, n_residues=n, n_chains=2)
        conf = _make_synthetic_conf_json(tmp_path, n)
        result = extract_confidence(cif, conf)
        assert len(result["chain_boundaries"]) == 2
        assert result["chain_boundaries"][0]["chain_id"] == "A"
        assert result["chain_boundaries"][1]["chain_id"] == "B"


# ---------------------------------------------------------------------------
# TestReport
# ---------------------------------------------------------------------------

from struct_predictor_core.report import generate_report, _confidence_band_breakdown, DISCLAIMER


class TestReport:
    def _make_inputs(self, n=10):
        plddt = np.array([50.0, 55.0, 65.0, 72.0, 80.0,
                          88.0, 91.0, 95.0, 85.0, 45.0], dtype=np.float32)[:n]
        pae = np.zeros((n, n), dtype=np.float32)
        chain_boundaries = [{"chain_id": "A", "start": 0, "end": n - 1}]
        sequences_info = [{"name": "TestProt", "sequence": "A" * n, "chain_id": "A"}]
        return plddt, pae, chain_boundaries, sequences_info

    def test_report_md_written(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict input.yaml --out_dir /tmp/boltz",
            input_label="test_protein.yaml", demo=False,
        )
        report = (tmp_path / "out" / "report.md").read_text()
        assert "# Struct Predictor Report" in report
        assert "ClawBio" in report

    def test_disclaimer_in_report(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict input.yaml", input_label="test.yaml", demo=False,
        )
        report = (tmp_path / "out" / "report.md").read_text()
        assert DISCLAIMER in report

    def test_result_json_written(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict input.yaml", input_label="test.yaml", demo=False,
        )
        result = json.loads((tmp_path / "out" / "result.json").read_text())
        assert "mean_plddt" in result
        assert "n_residues" in result
        assert "band_breakdown" in result

    def test_figures_created(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict input.yaml", input_label="test.yaml", demo=False,
        )
        assert (tmp_path / "out" / "figures" / "plddt.png").exists()
        assert (tmp_path / "out" / "figures" / "pae.png").exists()

    def test_viewer_html_written(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict input.yaml", input_label="test.yaml", demo=False,
        )
        html_path = tmp_path / "out" / "viewer.html"
        assert html_path.exists()
        html = html_path.read_text()
        assert "3Dmol" in html or "3dmol" in html
        assert "data_TEST" in html  # CIF content inlined


    def test_reproducibility_files(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        cmd_str = "boltz predict input.yaml --out_dir /tmp/boltz_raw"
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd=cmd_str, input_label="test.yaml", demo=False,
        )
        commands_sh = (tmp_path / "out" / "reproducibility" / "commands.sh").read_text()
        assert cmd_str in commands_sh
        assert (tmp_path / "out" / "reproducibility" / "environment.txt").exists()

    def test_band_breakdown_sums_to_100(self):
        plddt = np.array([45.0, 55.0, 65.0, 75.0, 85.0, 95.0])
        bd = _confidence_band_breakdown(plddt)
        assert abs(sum(bd.values()) - 100.0) < 0.01

    def test_band_breakdown_includes_100(self):
        """pLDDT = 100.0 must be counted in 'Very high', not lost."""
        plddt = np.array([100.0, 95.0, 85.0])
        bd = _confidence_band_breakdown(plddt)
        assert abs(sum(bd.values()) - 100.0) < 0.01
        assert bd["Very high"] > 0

    def test_demo_label_in_report(self, tmp_path):
        plddt, pae, cb, seqs = self._make_inputs()
        cif_src = tmp_path / "fake.cif"
        cif_src.write_text("data_TEST\n")
        generate_report(
            output_dir=tmp_path / "out", sequences_info=seqs, plddt=plddt,
            pae=pae, chain_boundaries=cb, cif_path=cif_src,
            cmd="boltz predict demo.yaml",
            input_label="Trp-cage miniprotein (demo)", demo=True,
        )
        report = (tmp_path / "out" / "report.md").read_text()
        assert "demo" in report.lower()


# ---------------------------------------------------------------------------
# TestPipeline + TestCLI
# ---------------------------------------------------------------------------

import struct_predictor
from struct_predictor import run_struct_prediction, DEMO_NAME


class TestDemoConstants:
    def test_demo_name(self):
        assert DEMO_NAME == "Trpcage"

    def test_demo_yaml_exists(self):
        from struct_predictor import _DEMO_YAML
        assert _DEMO_YAML.exists(), f"Demo YAML not found at {_DEMO_YAML}"

    def test_demo_yaml_has_sequence(self):
        from struct_predictor import _DEMO_YAML
        text = _DEMO_YAML.read_text()
        assert "NLYIQWLKDGGPSSGRPPPS" in text


class TestPipeline:
    def _make_fake_boltz_output(self, boltz_out_dir: Path, n: int = 20) -> None:
        pred_dir = boltz_out_dir / "predictions" / "Trpcage"
        pred_dir.mkdir(parents=True, exist_ok=True)
        cif_lines = [
            "data_Trpcage", "loop_",
            "_atom_site.group_PDB", "_atom_site.id",
            "_atom_site.label_atom_id", "_atom_site.label_asym_id",
            "_atom_site.label_seq_id", "_atom_site.Cartn_x",
            "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "_atom_site.B_iso_or_equiv",
        ]
        for i in range(n):
            cif_lines.append(f"ATOM {i+1} CA A {i+1} {i:.1f} 0.0 0.0 {90.0+i*0.5:.1f}")
        (pred_dir / "Trpcage_model_0.cif").write_text("\n".join(cif_lines))
        pae = [[float(abs(r - c)) for c in range(n)] for r in range(n)]
        (pred_dir / "confidence_Trpcage_model_0.json").write_text(
            json.dumps({"plddt": [90.0 + i * 0.5 for i in range(n)], "pae": pae})
        )

    def test_demo_end_to_end(self, tmp_path):
        mock_proc = MagicMock()
        mock_proc.returncode = 0

        def fake_run(cmd, **kwargs):
            if "predict" in cmd:
                out_dir_idx = cmd.index("--out_dir") + 1
                boltz_out = Path(cmd[out_dir_idx])
                self._make_fake_boltz_output(boltz_out, n=20)
            return mock_proc

        with patch("struct_predictor_core.predict.subprocess.run", side_effect=fake_run):
            result = run_struct_prediction(
                input_path=None, output_dir=tmp_path / "out", demo=True,
            )

        out = tmp_path / "out"
        # ClawBio artifacts at root
        assert (out / "report.md").exists()
        assert (out / "result.json").exists()
        assert (out / "figures" / "plddt.png").exists()
        assert (out / "figures" / "pae.png").exists()
        assert (out / "viewer.html").exists()
        # Boltz native layout preserved — CIF lives in predictions/<name>/
        assert (out / "predictions" / "Trpcage" / "Trpcage_model_0.cif").exists()
        assert result["demo"] is True
        assert result["n_residues"] == 20

    def test_no_input_no_demo_raises(self, tmp_path):
        with pytest.raises(ValueError, match="Provide --input or --demo"):
            run_struct_prediction(
                input_path=None, output_dir=tmp_path / "out", demo=False,
            )


class TestCLI:
    def test_demo_flag_parsed(self):
        parser = struct_predictor._build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/test"])
        assert args.demo is True
        assert args.output == "/tmp/test"

    def test_input_flag(self):
        parser = struct_predictor._build_parser()
        args = parser.parse_args(["--input", "prot.yaml", "--output", "/tmp/test"])
        assert args.input == "prot.yaml"

    def test_no_msa_flag_not_present(self):
        """--no-msa is removed; all runs are offline by default."""
        parser = struct_predictor._build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--demo", "--output", "/tmp/test", "--no-msa"])
