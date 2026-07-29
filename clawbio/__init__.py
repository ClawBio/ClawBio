"""ClawBio: Bioinformatics AI Agent shared library."""

__version__ = "0.6.0"

from .runner import list_skills, run_skill, upload_profile

__all__ = ["__version__", "run_skill", "list_skills", "upload_profile"]
