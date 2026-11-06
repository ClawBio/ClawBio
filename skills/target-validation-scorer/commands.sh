#!/bin/bash
# Target Validation Scorer - Reproducibility commands

# Demo mode (no API calls, uses cached TGFBR1/IPF data)
python skills/target-validation-scorer/target_validation_scorer.py --demo --output /tmp/target_validation_output

# Check output
cat /tmp/target_validation_output/report.md
ls /tmp/target_validation_output/figures/scoring_summary.png
