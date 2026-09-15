"""Execute the public first-skill notebook, including its pinned bootstrap."""
import json
import os
import re
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "docs" / "tutorial-first-skill.ipynb"


class FirstSkillNotebookTest(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("CLAWBIO_TEST_COLAB") == "1", "Opt-in notebook execution downloads the pinned environment")
    def test_clean_notebook_runs_and_checks_missing_evidence(self):
        import nbformat
        from nbclient import NotebookClient
        notebook = nbformat.read(NOTEBOOK, as_version=4)
        code = "\n".join(cell.source for cell in notebook.cells if cell.cell_type == "code")
        self.assertRegex(code, r'CLAWBIO_COMMIT = "[0-9a-f]{40}"')
        self.assertNotIn("summary.json", code)
        self.assertIn("--no-enrich", code)
        self.assertIn("--locked", code)
        with tempfile.TemporaryDirectory(prefix="clawbio-notebook-test-") as workspace:
            # Verification runs inside the same kernel, after all public cells.
            notebook.cells.append(nbformat.v4.new_code_cell('''
assert baseline["data"]["gene_profiles"]["CYP2C19"]["diplotype"].startswith("Indeterminate")
assert changed["data"]["gene_profiles"]["CYP2C19"]["diplotype"] == "NOT_TESTED"
assert baseline["summary"]["clinpgx_enriched"] == 0
assert changed["summary"]["clinpgx_enriched"] == 0
assert baseline["data"]["gene_profiles"]["DPYD"] == changed["data"]["gene_profiles"]["DPYD"]
assert bundle.is_file()
assert json.loads((run_root / "tutorial-provenance.json").read_text())["commit"] == CLAWBIO_COMMIT
assert "not a medical device" in (baseline_dir / "report.md").read_text()
'''))
            NotebookClient(notebook, timeout=900, kernel_name="python3", resources={"metadata": {"path": workspace}}).execute()

    def test_notebook_is_shareable_without_saved_results(self):
        notebook = json.loads(NOTEBOOK.read_text())
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                self.assertEqual(cell.get("outputs"), [])
                self.assertIsNone(cell.get("execution_count"))
            source = "".join(cell["source"])
            self.assertIsNone(re.search(r"(?m)^\s*(---|\*\*\*|___)\s*$", source))


if __name__ == "__main__":
    unittest.main()
