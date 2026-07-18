"""API entry point for skill registry and orchestrator integration."""

import sys
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

_mod_path = Path(__file__).parent / "population_equity_auditor.py"
_spec = spec_from_file_location("population_equity_auditor", _mod_path)
_mod = module_from_spec(_spec)
# Register before exec so the module's dataclasses (string annotations) resolve their fields.
sys.modules["population_equity_auditor"] = _mod
_spec.loader.exec_module(_mod)


def run(input_path: str, output_dir: str = "/tmp/population-equity-auditor") -> dict:
    """Run the skill programmatically. Returns result dict."""
    data = _mod.validate_input(Path(input_path))
    result = _mod.run_analysis(data)
    _mod.write_report(result, Path(output_dir))
    return result
