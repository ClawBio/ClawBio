"""Local Docker capacity selection and portable replay resource guards."""

from collections.abc import Mapping


def cap_to_profile_limits(cpus: int, memory_gb: int, limits: Mapping[str, int] | None) -> tuple[int, int]:
    """Bound launcher capacity by a test profile's own resourceLimits.

    --demo runs nf-core's ``test`` profile, whose resourceLimits define the
    canonical test (4 CPUs and 15 GB in the pinned releases). The launcher may
    lower them on a smaller machine, but never raises them. ``None`` leaves the
    capacity unchanged (a real run).
    """
    if not limits:
        return cpus, memory_gb
    return min(cpus, limits["cpus"]), min(memory_gb, limits["memory_gb"])


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
# The value of the first line matching an awk pattern, after its = or :.
_clawbio_config_value() {
  awk -v pattern="$1" '$0 ~ pattern {sub(/^[^=:]*[=:][[:space:]]*/, ""); sub(/[,[:space:]]+$/, ""); print; exit}' "$SOURCE_RESOURCE_CONFIG"
}
# Gigabytes in a Nextflow memory value such as '15.5.GB' or '16384 MB', or nothing.
_clawbio_gigabytes() {
  printf '%s\n' "$1" | awk '{
    value = $0; gsub(/["'\''[:space:]]/, "", value)
    if (!match(value, /[.]?[KMGT]B$/)) exit
    number = substr(value, 1, RSTART - 1); unit = substr(value, RSTART); sub(/^[.]/, "", unit)
    if (number !~ /^[0-9]+([.][0-9]+)?$/) exit
    scale = (unit == "KB") ? 1 / 1048576 : (unit == "MB") ? 1 / 1024 : (unit == "TB") ? 1024 : 1
    printf "%.6f\n", number * scale
  }'
}
_clawbio_bound_cpus() { awk -v recorded="$1" -v runtime="$RUNTIME_CPUS" 'BEGIN {print (recorded < runtime ? recorded : runtime)}'; }
# Whole gigabytes keep the N.GB form; a fraction is written as exact megabytes.
_clawbio_bound_memory() {
  awk -v recorded="$1" -v bytes="$RUNTIME_BYTES" 'BEGIN {
    limit = int(bytes / 1073741824 * 0.9); if (limit < 1) limit = 1
    gb = (recorded + 0 < limit) ? recorded + 0 : limit
    if (gb == int(gb)) printf "%d.GB\n", gb; else printf "%d.MB\n", int(gb * 1024 + 0.5)
  }'
}
# The executor lines are the run's aggregate budget; resourceLimits are per process
# (a demo keeps its test profile's smaller ones). Each is bounded on its own.
SOURCE_CPUS="$(_clawbio_config_value '^executor[.]cpus[[:space:]]*=')"
SOURCE_MEMORY="$(_clawbio_gigabytes "$(_clawbio_config_value '^executor[.]memory[[:space:]]*=')")"
LIMIT_CPUS="$(_clawbio_config_value '^[[:space:]]*cpus[[:space:]]*:')"
LIMIT_CPUS="${LIMIT_CPUS:-$SOURCE_CPUS}"
LIMIT_MEMORY="$(_clawbio_gigabytes "$(_clawbio_config_value '^[[:space:]]*memory[[:space:]]*:')")"
LIMIT_MEMORY="${LIMIT_MEMORY:-$SOURCE_MEMORY}"
if [[ -z "$RUNTIME_CPUS" || -z "$RUNTIME_BYTES" || -z "$SOURCE_CPUS" || -z "$SOURCE_MEMORY" ]]; then
  echo "ERROR: Cannot determine Docker capacity or recorded resource limits. Start Docker and check resource_limits.config." >&2
  exit 1
fi
case "$RUNTIME_CPUS $RUNTIME_BYTES $SOURCE_CPUS $LIMIT_CPUS" in
  *[!0-9\ ]*) echo "ERROR: Invalid Docker capacity or recorded resource limits." >&2; exit 1 ;;
esac
case "$SOURCE_MEMORY $LIMIT_MEMORY" in
  *[!0-9.\ ]*) echo "ERROR: Invalid recorded resource limits." >&2; exit 1 ;;
esac
# One resourceLimits block is what the wrapper writes; a second one (a hand-added
# withLabel, say) would be rewritten with the first block's limits, so refuse it.
_clawbio_count() { awk -v pattern="$1" '$0 ~ pattern {n++} END {print n + 0}' "$SOURCE_RESOURCE_CONFIG"; }
if [[ "$(_clawbio_count '^[[:space:]]*cpus[[:space:]]*:')" -gt 1 || "$(_clawbio_count '^[[:space:]]*memory[[:space:]]*:')" -gt 1 ]]; then
  echo "ERROR: resource_limits.config has more than one resourceLimits block; this replay can bound only one. Keep a single block." >&2
  exit 1
fi
BOUNDED_CPUS="$(_clawbio_bound_cpus "$SOURCE_CPUS")"
BOUNDED_MEMORY="$(_clawbio_bound_memory "$SOURCE_MEMORY")"
BOUNDED_LIMIT_CPUS="$(_clawbio_bound_cpus "$LIMIT_CPUS")"
BOUNDED_LIMIT_MEMORY="$(_clawbio_bound_memory "$LIMIT_MEMORY")"
awk -v cpus="$BOUNDED_CPUS" -v memory="$BOUNDED_MEMORY" -v limit_cpus="$BOUNDED_LIMIT_CPUS" \
    -v limit_memory="$BOUNDED_LIMIT_MEMORY" -v quote="'" '
  /^executor[.]cpus[[:space:]]*=/ {print "executor.cpus = " cpus; next}
  /^executor[.]memory[[:space:]]*=/ {print "executor.memory = " quote memory quote; next}
  /^[[:space:]]*cpus[[:space:]]*:/ {print "        cpus: " limit_cpus ","; next}
  /^[[:space:]]*memory[[:space:]]*:/ {print "        memory: " quote limit_memory quote ","; next}
  {print}
' "$SOURCE_RESOURCE_CONFIG" > "$REPLAY_RESOURCE_CONFIG"
RESOURCE_CONFIG="-c $REPLAY_RESOURCE_CONFIG"
trap 'rm -f "$REPLAY_RESOURCE_CONFIG"' EXIT

"""
