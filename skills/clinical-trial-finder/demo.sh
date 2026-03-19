#!/usr/bin/env bash
# demo.sh -- End-to-end clinical-trial-finder demonstration
#
# Runs the full pipeline: demo, gene mode, rsID mode, profile-report integration.
# No API keys needed -- all data sources are public and free.
#
# Usage:  bash skills/clinical-trial-finder/demo.sh
# From:   ClawBio repo root
#
# Windows: This script requires bash (Git Bash or WSL).
#          Alternatively, run the Python commands directly:
#            python3 skills/clinical-trial-finder/clinical_trial_finder.py --demo --fhir --output out/demo
#            python3 skills/clinical-trial-finder/clinical_trial_finder.py --gene TP53 --fhir --output out/gene
#            python3 skills/clinical-trial-finder/clinical_trial_finder.py --rsid rs3798220 --fhir --output out/rsid

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT="/tmp/ctf_demo_$(date +%s)"

cd "$REPO_DIR"

echo "============================================================"
echo "  ClawBio Clinical Trial Finder -- End-to-End Demo"
echo "============================================================"
echo ""
echo "Output directory: $OUT"
echo ""

# --- 1. Demo mode (built-in BRCA1 data) ---
echo "[1/5] Demo mode: BRCA1 breast cancer"
echo "--------------------------------------------------------------"
python3 skills/clinical-trial-finder/clinical_trial_finder.py \
  --demo --fhir --output "$OUT/demo"
echo ""

# --- 2. Gene mode (OpenTargets enrichment) ---
echo "[2/5] Gene mode: TP53 via OpenTargets"
echo "--------------------------------------------------------------"
python3 skills/clinical-trial-finder/clinical_trial_finder.py \
  --gene TP53 --fhir --output "$OUT/gene_tp53"
echo ""

# --- 3. rsID mode (GWAS Catalog variant lookup) ---
echo "[3/5] rsID mode: rs3798220 (LPA / coronary artery disease)"
echo "--------------------------------------------------------------"
python3 skills/clinical-trial-finder/clinical_trial_finder.py \
  --rsid rs3798220 --fhir --output "$OUT/rsid_rs3798220"
echo ""

# --- 4. Country filter ---
echo "[4/5] Country filter: BRCA1 trials in United States only"
echo "--------------------------------------------------------------"
python3 skills/clinical-trial-finder/clinical_trial_finder.py \
  --demo --country "United States" --output "$OUT/country_us"
echo ""

# --- 5. Profile-report integration (if available) ---
echo "[5/5] Profile report integration (PRS -> trials)"
echo "--------------------------------------------------------------"
if python3 -c "import numpy" 2>/dev/null; then
  python3 skills/profile-report/profile_report.py \
    --demo --output "$OUT/profile"
  echo ""
  echo "Profile report sections:"
  grep "^## " "$OUT/profile/profile_report.md" | head -10
else
  echo "  Skipped (numpy not installed -- run: pip install numpy)"
fi

echo ""
echo "============================================================"
echo "  Demo complete!"
echo "============================================================"
echo ""
echo "Outputs:"
find "$OUT" -type f | sort | while read -r f; do
  size=$(du -h "$f" | cut -f1)
  echo "  $size  ${f#$OUT/}"
done
echo ""
echo "Verify checksums:"
for dir in "$OUT"/*/; do
  if [ -f "$dir/checksums.sha256" ]; then
    name=$(basename "$dir")
    (cd "$dir" && sha256sum -c checksums.sha256 2>&1 | tail -1) && echo "    $name: ALL OK" || echo "    $name: FAIL"
  fi
done
echo ""
echo "Open HTML reports in browser:"
find "$OUT" -name "report.html" | while read -r f; do
  echo "  file://$f"
done
