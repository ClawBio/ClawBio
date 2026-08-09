#!/usr/bin/env bash
################################################################################
# claw-amplicon-qc end-to-end test
#
# Runs the skill in --demo mode (bundled tiny FASTQ fixture at
# tests/fixtures/) and verifies:
#   1. The script exits 0
#   2. All expected output files are produced
#   3. Output files have non-zero content
#   4. Overall retention is in a plausible range (>50%)
#   5. Number of samples matches the fixture (1 pair)
#
# Usage:
#   cd skills/claw-amplicon-qc/
#   bash tests/run_test.sh
#
# Requires: the amplicon-qc conda environment activated, or the tools
# (Rscript, cutadapt, seqkit) otherwise on PATH.
################################################################################

set -euo pipefail

# Locate the skill root regardless of where the test was invoked from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(dirname "$SCRIPT_DIR")"
SCRIPT="$SKILL_DIR/amplicon_qc.R"

# Isolated temp output — cleaned up on exit whether we pass or fail
OUTPUT_DIR="$(mktemp -d -t claw_amplicon_qc_test.XXXXXX)"
trap 'rm -rf "$OUTPUT_DIR"' EXIT

# Colours if the terminal supports them
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'
  RED='\033[0;31m'
  YELLOW='\033[0;33m'
  NC='\033[0m'
else
  GREEN=''; RED=''; YELLOW=''; NC=''
fi

pass_count=0
fail_count=0

pass() {
  echo -e "  ${GREEN}✓${NC} $1"
  pass_count=$((pass_count + 1))
}

fail() {
  echo -e "  ${RED}✗${NC} $1"
  fail_count=$((fail_count + 1))
}

echo "─────────────────────────────────────────────────────────"
echo "  claw-amplicon-qc — end-to-end test"
echo "─────────────────────────────────────────────────────────"
echo "  Skill:       $SKILL_DIR"
echo "  Test output: $OUTPUT_DIR"
echo "─────────────────────────────────────────────────────────"
echo ""

# --------------------------------------------------------------------------
# Prerequisites
# --------------------------------------------------------------------------
echo "Checking prerequisites..."
for tool in Rscript cutadapt seqkit; do
  if command -v "$tool" >/dev/null 2>&1; then
    pass "$tool on PATH"
  else
    fail "$tool NOT found on PATH — activate the amplicon-qc conda env?"
  fi
done

if [[ -f "$SCRIPT" ]]; then
  pass "amplicon_qc.R found"
else
  fail "amplicon_qc.R not found at $SCRIPT"
fi

if [[ -f "$SKILL_DIR/tests/fixtures/demo_R1.fastq.gz" ]] && \
   [[ -f "$SKILL_DIR/tests/fixtures/demo_R2.fastq.gz" ]]; then
  pass "bundled fixture present"
else
  fail "bundled fixture missing — expected demo_R1.fastq.gz and demo_R2.fastq.gz in tests/fixtures/"
fi

# Bail early if any prerequisite failed
if [[ $fail_count -gt 0 ]]; then
  echo ""
  echo -e "${RED}Prerequisite checks failed — cannot run pipeline test.${NC}"
  exit 1
fi

echo ""

# --------------------------------------------------------------------------
# Run the pipeline in --demo mode
# --------------------------------------------------------------------------
echo "Running pipeline in --demo mode..."
echo ""

if Rscript "$SCRIPT" --demo --output "$OUTPUT_DIR" 2>&1 | sed 's/^/    /'; then
  pass "pipeline exited 0"
else
  fail "pipeline exited non-zero"
  exit 1
fi

echo ""
echo "Checking outputs..."

# --------------------------------------------------------------------------
# Output file checks
# --------------------------------------------------------------------------
expected_files=(
  "01_raw_stats.txt"
  "02_filtN_stats.txt"
  "03_trimmed_stats.txt"
  "04_cutadapt_log.txt"
  "report.md"
  "qc_summary.json"
)

for f in "${expected_files[@]}"; do
  if [[ -s "$OUTPUT_DIR/$f" ]]; then
    pass "$f present and non-empty"
  else
    fail "$f missing or empty"
  fi
done

expected_dirs=(
  "filtN"
  "cutadapt_trimmed"
)

for d in "${expected_dirs[@]}"; do
  if [[ -d "$OUTPUT_DIR/$d" ]] && [[ -n "$(ls -A "$OUTPUT_DIR/$d")" ]]; then
    pass "$d/ present and non-empty"
  else
    fail "$d/ missing or empty"
  fi
done

# --------------------------------------------------------------------------
# Sanity-check the actual retention numbers
# --------------------------------------------------------------------------

# Fixture is one sample pair; report should say "successfully: 1"
if grep -q "Samples processed successfully: \*\*1\*\*" "$OUTPUT_DIR/report.md" 2>/dev/null; then
  pass "report shows 1 sample processed successfully"
else
  fail "report doesn't show 1 successful sample"
fi

if grep -q "Samples failed (cutadapt): \*\*0\*\*" "$OUTPUT_DIR/report.md" 2>/dev/null; then
  pass "report shows 0 failed samples"
else
  fail "report shows failures or missing summary"
fi

# Overall retention in qc_summary.json should be >50% for the demo fixture
# (real 515F/806R data — retention ~90% expected; we set the bar low to
# tolerate small parameter changes in the pipeline over time)
if command -v python3 >/dev/null 2>&1; then
  pct=$(python3 -c "
import json, sys
try:
    with open('$OUTPUT_DIR/qc_summary.json') as f:
        d = json.load(f)
    total_in = d['totals']['raw_reads_total']
    total_out = d['totals']['after_cutadapt_total']
    if total_in == 0:
        print(0)
    else:
        print(int(100 * total_out / total_in))
except Exception as e:
    print(-1, file=sys.stderr)
    sys.exit(1)
" 2>/dev/null || echo "-1")

  if [[ "$pct" -ge 50 ]]; then
    pass "overall retention $pct% (>= 50% threshold)"
  elif [[ "$pct" -eq -1 ]]; then
    fail "could not parse qc_summary.json"
  else
    fail "overall retention $pct% is below 50% — pipeline may be misconfigured"
  fi
else
  echo -e "  ${YELLOW}~${NC} python3 not available — skipping retention sanity check"
fi

# --------------------------------------------------------------------------
# Result
# --------------------------------------------------------------------------
echo ""
echo "─────────────────────────────────────────────────────────"
if [[ $fail_count -eq 0 ]]; then
  echo -e "  ${GREEN}All $pass_count checks passed.${NC}"
  echo "─────────────────────────────────────────────────────────"
  exit 0
else
  echo -e "  ${RED}$fail_count check(s) failed${NC} ($pass_count passed)"
  echo "─────────────────────────────────────────────────────────"
  exit 1
fi
