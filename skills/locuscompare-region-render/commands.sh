#!/usr/bin/env bash
# Example invocations for locuscompare-region-render.
#
# Each run writes its own replay script to `<dir>/reproducibility/commands.sh`
# (with environment.yml and checksums.sha256) via clawbio.common.reproducibility.
# This file is only a hand-written usage example.

set -euo pipefail

# Run the bundled offline demo:
python skills/locuscompare-region-render/cli.py --demo --output /tmp/locuscompare_demo

# Or run with a real config:
# python skills/locuscompare-region-render/cli.py \
#     --input <CONFIG_PATH> \
#     --output <OUTPUT_DIR>
