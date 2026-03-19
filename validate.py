import yaml
from pathlib import Path

skill_path = Path("skills/marcelo-emr-fhir-drug-analyser/SKILL.md")

if not skill_path.exists():
    print(f"File not found: {skill_path}")
    exit(1)

text = skill_path.read_text(encoding="utf-8")

# Split frontmatter safely
parts = text.split("---")
if len(parts) < 3:
    print("Invalid format: missing YAML frontmatter")
    exit(1)

front = parts[1]

try:
    meta = yaml.safe_load(front)
except Exception as e:
    print(f"YAML parsing error: {e}")
    exit(1)

required = ["name", "version", "author", "domain", "description", "inputs", "outputs"]
missing = [f for f in required if f not in meta]

if missing:
    print(f"Missing fields: {missing}")
else:
    print(f"SKILL.md valid: {meta['name']} v{meta['version']}")