"""Local Docker capacity selection and portable replay resource guards."""


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
# The recorded limits are the executor.cpus and executor.memory lines the wrapper
# wrote; memory keeps its decimal point and unit (a hand-edited 15.5.GB is 15.5 GB).
SOURCE_CPUS="$(awk '/^executor\.cpus[[:space:]]*=/ {sub(/^[^=]*=[[:space:]]*/, ""); sub(/[[:space:]]+$/, ""); print; exit}' "$SOURCE_RESOURCE_CONFIG")"
SOURCE_MEMORY="$(awk '/^executor\.memory[[:space:]]*=/ {
  value = $0; sub(/^[^=]*=[[:space:]]*/, "", value); gsub(/["'\''[:space:]]/, "", value)
  if (!match(value, /[.]?[KMGT]B$/)) exit
  number = substr(value, 1, RSTART - 1); unit = substr(value, RSTART); sub(/^[.]/, "", unit)
  if (number !~ /^[0-9]+([.][0-9]+)?$/) exit
  scale = (unit == "KB") ? 1 / 1048576 : (unit == "MB") ? 1 / 1024 : (unit == "TB") ? 1024 : 1
  printf "%.6f\n", number * scale; exit
}' "$SOURCE_RESOURCE_CONFIG")"
if [[ -z "$RUNTIME_CPUS" || -z "$RUNTIME_BYTES" || -z "$SOURCE_CPUS" || -z "$SOURCE_MEMORY" ]]; then
  echo "ERROR: Cannot determine Docker capacity or recorded resource limits. Start Docker and check resource_limits.config." >&2
  exit 1
fi
case "$RUNTIME_CPUS $RUNTIME_BYTES $SOURCE_CPUS" in
  *[!0-9\ ]*) echo "ERROR: Invalid Docker capacity or recorded resource limits." >&2; exit 1 ;;
esac
case "$SOURCE_MEMORY" in
  *[!0-9.]*) echo "ERROR: Invalid recorded resource limits." >&2; exit 1 ;;
esac
BOUNDED_CPUS="$(awk -v original="$SOURCE_CPUS" -v runtime="$RUNTIME_CPUS" 'BEGIN {print (original < runtime ? original : runtime)}')"
# Whole gigabytes keep the N.GB form; a fraction is written as exact megabytes.
BOUNDED_MEMORY="$(awk -v original="$SOURCE_MEMORY" -v bytes="$RUNTIME_BYTES" 'BEGIN {
  limit = int(bytes / 1073741824 * 0.9); if (limit < 1) limit = 1
  gb = (original + 0 < limit) ? original + 0 : limit
  if (gb == int(gb)) printf "%d.GB\n", gb; else printf "%d.MB\n", int(gb * 1024 + 0.5)
}')"
awk -v cpus="$BOUNDED_CPUS" -v memory="$BOUNDED_MEMORY" -v quote="'" '
  /^executor\.cpus[[:space:]]*=/ {print "executor.cpus = " cpus; next}
  /^executor\.memory[[:space:]]*=/ {print "executor.memory = " quote memory quote; next}
  /^[[:space:]]*cpus[[:space:]]*:/ {print "        cpus: " cpus ","; next}
  /^[[:space:]]*memory[[:space:]]*:/ {print "        memory: " quote memory quote ","; next}
  {print}
' "$SOURCE_RESOURCE_CONFIG" > "$REPLAY_RESOURCE_CONFIG"
RESOURCE_CONFIG="-c $REPLAY_RESOURCE_CONFIG"
trap 'rm -f "$REPLAY_RESOURCE_CONFIG"' EXIT

"""
