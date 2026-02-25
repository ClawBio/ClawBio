# Live Demo Script — Equity Scorer

Run this during Slide 14 (Tip 9). Total time: 2-3 minutes.

## Pre-Demo Checklist (do before the talk)

```bash
# Terminal setup
cd ~/ClawBio/examples
clear

# Verify everything works
python3 ../skills/equity-scorer/equity_scorer.py \
    --input demo_populations.vcf \
    --pop-map demo_population_map.csv \
    --output /tmp/test_demo
# Should complete in <5 seconds

# Clean up
rm -rf /tmp/test_demo

# Pre-open the report and figures
# Have examples/demo_report/report.md open in a browser tab (Markdown preview)
# Have examples/demo_report/figures/ ready to show
```

## Font Size

Set terminal font to 24pt+ so the back row can read it.

```bash
# macOS Terminal: Cmd+Shift+. to increase font size
# iTerm2: Cmd+= to zoom in
```

## The Demo (what to type and say)

### Step 1: Show the input

```bash
head -3 demo_population_map.csv
```

Say: "Here's our input. 50 samples across 5 populations. Notice the European overrepresentation: 22 out of 50 samples."

### Step 2: Run the scorer

```bash
python3 ../skills/equity-scorer/equity_scorer.py \
    --input demo_populations.vcf \
    --pop-map demo_population_map.csv \
    --output demo_report
```

Say: "It's parsing the VCF, computing real heterozygosity from genotype data, calculating pairwise FST between every population pair, running PCA on the genotype matrix, and computing the HEIM Equity Score."

Wait for output. Point to key numbers as they appear:
- "AFR has the highest heterozygosity, as expected from the out-of-Africa model"
- "AFR vs EAS has the highest FST at 0.10"
- "HEIM Score: 76 out of 100. Good, but not excellent."

### Step 3: Show the report

```bash
# Open the report in the browser
open demo_report/report.md
# Or if no markdown viewer:
cat demo_report/report.md | head -40
```

Say: "It generates a full markdown report with methods, results, and a reproducibility block. Every figure is saved as PNG. Every table as CSV."

### Step 4: Show figures (switch to slides)

Go to Slide 16 (PCA), 17 (FST heatmap), 18 (ancestry + het).

Say: "Five publication-quality figures. The PCA shows clear population structure. The FST heatmap quantifies divergence. And the ancestry chart shows exactly where the bias is."

## If the Demo Fails

Fallback: skip to Slide 15 (terminal output screenshot) and say "Let me show you what the output looks like" while showing the pre-generated figures on slides 16-19.

The pre-generated output is already in `examples/demo_report/`. It will always be there as a safety net.

## After the Demo

Return to Slide 19 (HEIM Score breakdown), then continue to Slide 20 (Why This Matters).
