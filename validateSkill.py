#!/usr/bin/env python3
"""
Validate the structure of a SKILL.md file.

Checks for required front matter fields and mandatory sections.
Usage: python validate_skill.py <path_to_SKILL.md>
"""

import sys
import yaml
import argparse
from pathlib import Path

# Required front matter fields
REQUIRED_FIELDS = ['name', 'version', 'author', 'domain', 'inputs', 'outputs', 'demo_data']

# Required section headings (must appear in the Markdown after front matter)
REQUIRED_SECTIONS = ['## Domain Decisions', '## Safety Rules', '## Agent Boundary']

def extract_front_matter(content):
    """Extract YAML front matter from the markdown content."""
    parts = content.split('---')
    if len(parts) < 3:
        return None, "No valid front matter found (need at least two '---' separators)."
    # The first part is empty or leading text, second is the front matter
    front_matter = parts[1].strip()
    if not front_matter:
        return None, "Front matter is empty."
    return front_matter, None

def validate_skill_file(file_path):
    """Run all validations on the given SKILL.md file."""
    path = Path(file_path)
    if not path.is_file():
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

    content = path.read_text(encoding='utf-8')

    # Extract and parse front matter
    front_yaml, err = extract_front_matter(content)
    if err:
        print(f"Error: {err}")
        sys.exit(1)

    try:
        meta = yaml.safe_load(front_yaml)
    except yaml.YAMLError as e:
        print(f"Error: Invalid YAML in front matter:\n{e}")
        sys.exit(1)

    if not isinstance(meta, dict):
        print("Error: Front matter is not a YAML dictionary.")
        sys.exit(1)

    # Perform checks
    checks = []

    # Check front matter fields
    for field in REQUIRED_FIELDS:
        ok = field in meta
        checks.append((ok, f"field '{field}'"))

    # Check required sections
    # Remove front matter part to search in the rest
    rest_of_doc = '---'.join(content.split('---')[2:])  # after second '---'
    for section in REQUIRED_SECTIONS:
        ok = section in rest_of_doc
        checks.append((ok, f"section '{section}'"))

    # Print results
    failed = []
    for ok, label in checks:
        status = 'PASS' if ok else 'FAIL'
        print(f'  [{status}] {label}')
        if not ok:
            failed.append(label)

    if failed:
        print(f'\n{len(failed)} check(s) failed. Fix these before submitting.')
        sys.exit(1)
    else:
        print('\nAll checks passed!')
        sys.exit(0)

def main():
    parser = argparse.ArgumentParser(description='Validate a SKILL.md file structure.')
    parser.add_argument('file', help='Path to SKILL.md')
    args = parser.parse_args()
    validate_skill_file(args.file)

if __name__ == '__main__':
    main()