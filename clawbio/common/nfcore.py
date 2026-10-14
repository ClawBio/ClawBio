"""Input and runtime helpers shared by the pinned nf-core wrappers."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path


def read_delimited_samplesheet(
    path: Path, *, delimiter: str = ","
) -> tuple[list[str], list[dict[str, str]]]:
    """Reject ambiguous headers and ragged records before creating row maps."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter, strict=True)
        header = next(reader, [])
        if not header or any(not name.strip() for name in header):
            raise ValueError("Samplesheet contains an empty column name.")
        if len(set(header)) != len(header):
            raise ValueError("Samplesheet contains duplicate column names.")
        rows = []
        for record in reader:
            if not record:
                continue
            if len(record) != len(header):
                raise ValueError(
                    f"Record ending at line {reader.line_num} has {len(record)} fields; expected {len(header)}."
                )
            rows.append(dict(zip(header, record)))
    return header, rows


def read_log_tail(path: Path, limit: int = 65536) -> str:
    """Read at most limit bytes even when a run produced very large logs."""
    try:
        with path.open("rb") as handle:
            handle.seek(0, 2)
            handle.seek(max(0, handle.tell() - limit))
            return handle.read(limit).decode("utf-8", errors="replace")
    except OSError:
        return ""


def finite_nonnegative_hours(value: str) -> float:
    hours = float(value)
    if not math.isfinite(hours) or hours < 0:
        raise ValueError("Hours must be finite and >= 0 (0 disables the timeout).")
    return hours


_CONFIG_TOKEN_RE = re.compile(
    r'/\*.*?\*/|//[^\n]*|\'(?:\\.|[^\'\\])*\'|"(?:\\.|[^"\\])*"|[A-Za-z_][A-Za-z_0-9]*|[^\s]',
    re.DOTALL,
)


def read_manifest_version(path: Path) -> str:
    """Read a literal top-level manifest version without interpreting Groovy."""
    try:
        tokens = [
            token
            for token in _CONFIG_TOKEN_RE.findall(path.read_text(encoding="utf-8"))
            if not token.startswith(("//", "/*"))
        ]
    except (OSError, UnicodeError):
        return ""
    depth, in_manifest = 0, False
    for index, token in enumerate(tokens):
        if (
            token == "manifest"
            and depth == 0
            and tokens[index + 1 : index + 2] == ["{"]
        ):
            in_manifest = True
        elif token == "{":
            depth += 1
        elif token == "}":
            depth -= 1
            if in_manifest and depth == 0:
                return ""
        elif in_manifest and depth == 1 and token == "version":
            if tokens[index - 1 : index] == ["."]:
                continue
            assignment = tokens[index + 1 : index + 3]
            if (
                len(assignment) == 2
                and assignment[0] == "="
                and assignment[1][0] in "'\""
            ):
                return assignment[1][1:-1].strip()
    return ""


def load_unique_yaml(text: str):
    """Load safe YAML with explicit duplicate mapping keys rejected."""
    import yaml

    class UniqueLoader(yaml.SafeLoader):
        def construct_mapping(self, node, deep=False):
            seen = set()
            for key_node, _ in node.value:
                if key_node.tag == "tag:yaml.org,2002:merge":
                    continue
                key = self.construct_object(key_node, deep=deep)
                if not isinstance(key, str) or key in seen:
                    raise yaml.constructor.ConstructorError(
                        "while reading a mapping",
                        node.start_mark,
                        f"Expected a unique string key, got {key!r}",
                        key_node.start_mark,
                    )
                seen.add(key)
            return super().construct_mapping(node, deep=deep)

    return yaml.load(text, Loader=UniqueLoader)


def validate_parameter_scalars(
    params: dict, contract: dict, *, patterns: bool = True
) -> None:
    """Check pinned scalar constraints; leave relationships to preflight.

    nf-schema ignores null/false unset values, including reference sentinels.
    Quoted booleans stay strings so they cannot silently become truthy flags.
    """
    for name, value in params.items():
        spec = contract.get(name)
        if not spec or value is None or value is False or value == "":
            continue
        kind = spec.get("type")
        allowed = [kind] if isinstance(kind, str) else (kind or [])
        numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
        matches = {
            "string": isinstance(value, str),
            "boolean": isinstance(value, bool),
            "integer": numeric
            and (not isinstance(value, float) or math.isfinite(value))
            and int(value) == value,
            "number": numeric
            and (not isinstance(value, float) or math.isfinite(value)),
            "object": isinstance(value, dict),
            "array": isinstance(value, list),
        }
        if allowed and not any(matches.get(option, False) for option in allowed):
            raise ValueError(
                f"{name} must have type {kind}; got {type(value).__name__}. Numeric values must be finite."
            )
        if patterns:
            if "enum" in spec and value not in spec["enum"]:
                raise ValueError(f"{name} must be one of {spec['enum']}.")
            if isinstance(value, str):
                if "pattern" in spec and not re.search(spec["pattern"], value):
                    raise ValueError(
                        f"{name} does not match the upstream pattern {spec['pattern']!r}."
                    )
                if "minLength" in spec and len(value) < spec["minLength"]:
                    raise ValueError(
                        f"{name} is shorter than its upstream minimum length."
                    )
                if "maxLength" in spec and len(value) > spec["maxLength"]:
                    raise ValueError(f"{name} exceeds its upstream maximum length.")
        if numeric:
            if "minimum" in spec and value < spec["minimum"]:
                raise ValueError(f"{name} must be >= {spec['minimum']}.")
            if "maximum" in spec and value > spec["maximum"]:
                raise ValueError(f"{name} must be <= {spec['maximum']}.")
            branches = spec.get("oneOf", [])
            if branches and not any(
                ("const" in branch and value == branch["const"])
                or ("minimum" in branch and value >= branch["minimum"])
                for branch in branches
            ):
                raise ValueError(
                    f"{name} does not satisfy its numeric alternatives: {branches}."
                )


def canonical_work_dir(args, output_dir: Path) -> str:
    raw = getattr(args, "work_dir", None)
    if not raw:
        return "upstream/work"
    if "://" in str(raw):
        return str(raw).strip()
    return Path(raw).expanduser().resolve().as_posix()


def write_manifest_checkpoint(output_dir: Path, manifest: dict) -> Path:
    """Commit a launch record before execution so failures can safely resume."""
    import json
    from clawbio.common.textio import write_text_lf_atomic

    path = output_dir / "reproducibility/manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf_atomic(
        path, json.dumps({**manifest, "run_state": "prepared"}, indent=2) + "\n"
    )
    return path


def load_unique_json(text: str):
    """Reject repeated JSON keys instead of silently taking the final value."""
    import json

    def mapping(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON mapping key: {key!r}")
            result[key] = value
        return result

    return json.loads(text, object_pairs_hook=mapping)


def docker_cpu_count() -> int | None:
    """Best-effort CPU capacity of the active Docker daemon."""
    import shutil
    import subprocess

    if not shutil.which("docker"):
        return None
    try:
        result = subprocess.run(
            ["docker", "info", "--format", "{{.NCPU}}"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0 and result.stdout.strip().isdigit():
            return int(result.stdout.strip()) or None
    except (OSError, subprocess.SubprocessError):
        # Capacity probing is optional; callers fall back to host capacity.
        pass
    return None


def read_json_mapping(path: Path) -> dict:
    value = load_unique_json(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("The manifest must contain a JSON object.")
    return value


def uses_local_docker_capacity(args, supported_profiles: set[str]) -> bool:
    """Use launcher capacity only for recognized local Docker profiles.

    Site/cloud profiles may submit to larger compute nodes. Their resource
    settings belong to the site configuration, rather than the launch machine.
    """
    parts = {
        part.strip()
        for part in str(getattr(args, "profile", "")).split(",")
        if part.strip()
    }
    cloud = any(
        any(token in part for token in ("aws", "azure", "gcp")) for part in parts
    )
    remote_work = "://" in str(getattr(args, "work_dir", "") or "")
    return (
        "docker" in parts
        and parts <= supported_profiles
        and not cloud
        and not remote_work
    )


def build_runtime_resource_guard() -> str:
    """Standalone Bash 3.2 guard that bounds recorded limits to the replay VM.

    Only a temporary generated config is changed; the checksummed original stays
    intact. This fragment is embedded in Nextflow-direct replay scripts.
    """
    return r"""# Bound recorded resource limits to the current Docker runtime.
RESOURCE_CONFIG="-c reproducibility/resource_limits.config"
SOURCE_RESOURCE_CONFIG="reproducibility/resource_limits.config"
REPLAY_RESOURCE_CONFIG="reproducibility/resource_limits.replay.config"
RUNTIME_CAPACITY="$(docker info --format '{{.NCPU}} {{.MemTotal}}' 2>/dev/null || true)"
RUNTIME_CPUS="$(printf '%s\n' "$RUNTIME_CAPACITY" | awk '{print $1}')"
RUNTIME_BYTES="$(printf '%s\n' "$RUNTIME_CAPACITY" | awk '{print $2}')"
SOURCE_CPUS="$(awk '/cpus[[:space:]]*[:=]/ {gsub(/[^0-9]/, ""); print; exit}' "$SOURCE_RESOURCE_CONFIG")"
SOURCE_MEMORY="$(awk '/memory[[:space:]]*[:=]/ {gsub(/[^0-9]/, ""); print; exit}' "$SOURCE_RESOURCE_CONFIG")"
if [[ -z "$RUNTIME_CPUS" || -z "$RUNTIME_BYTES" || -z "$SOURCE_CPUS" || -z "$SOURCE_MEMORY" ]]; then
  echo "ERROR: Cannot determine Docker capacity or recorded resource limits. Start Docker and check resource_limits.config." >&2
  exit 1
fi
case "$RUNTIME_CPUS $RUNTIME_BYTES $SOURCE_CPUS $SOURCE_MEMORY" in
  *[!0-9\ ]*) echo "ERROR: Invalid Docker capacity or recorded resource limits." >&2; exit 1 ;;
esac
BOUNDED_CPUS="$(awk -v original="$SOURCE_CPUS" -v runtime="$RUNTIME_CPUS" 'BEGIN {print (original < runtime ? original : runtime)}')"
BOUNDED_MEMORY="$(awk -v original="$SOURCE_MEMORY" -v bytes="$RUNTIME_BYTES" 'BEGIN {limit=int(bytes/1073741824*0.9); if(limit<1) limit=1; print (original < limit ? original : limit)}')"
awk -v cpus="$BOUNDED_CPUS" -v memory="$BOUNDED_MEMORY" -v quote="'" '
  /^executor\.cpus[[:space:]]*=/ {print "executor.cpus = " cpus; next}
  /^executor\.memory[[:space:]]*=/ {print "executor.memory = " quote memory ".GB" quote; next}
  /^[[:space:]]*cpus[[:space:]]*:/ {print "        cpus: " cpus ","; next}
  /^[[:space:]]*memory[[:space:]]*:/ {print "        memory: " quote memory ".GB" quote ","; next}
  {print}
' "$SOURCE_RESOURCE_CONFIG" > "$REPLAY_RESOURCE_CONFIG"
RESOURCE_CONFIG="-c $REPLAY_RESOURCE_CONFIG"
trap 'rm -f "$REPLAY_RESOURCE_CONFIG"' EXIT

"""


def json_safe_details(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe_details(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe_details(item) for item in value]
    return value
