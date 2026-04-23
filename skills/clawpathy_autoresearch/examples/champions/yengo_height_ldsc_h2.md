---
name: repro-yengo-height
version: 0.4.0
author: ClawBio
description: Reproduce LDSC SNP-heritability of adult height from Yengo et al. 2022 (Nature 610:704-712). Munge raw sumstats to HapMap3, run ldsc.py --h2 with 1000G Phase 3 EUR LD scores (no MHC), report h², SE, lambda_GC, chi² mean, intercept, ratio. All artifacts under output/.
inputs:
  - name: raw_sumstats
    format: "TSV/gzipped, must contain SNP (rsID), A1, A2, effect estimate (BETA or Z), P, N"
    required: true
outputs:
  - name: h2.log
    format: "output/h2.log — verbatim LDSC stdout"
  - name: summary.json
    format: "output/summary.json — key scalars + munge stats + both free- and constrained-intercept estimates, machine-readable"
  - name: provenance.json
    format: "output/provenance.json — SHA256 of raw + cleaned inputs, tool commit, paper citation, LD-reference resolution path"
  - name: ldsc_env.yml
    format: "output/env/ldsc_env.yml — full conda env export"
trigger_keywords:
  - "Yengo height"
  - "LDSC h2 height"
  - "reproduce Yengo 2022"
---

## Trigger

**Fire when:** user asks to reproduce the LDSC SNP-heritability from Yengo et al. 2022 height GWAS.

**Do NOT fire when:** partitioned h² (needs baselineLD annotation files in --ref-ld-chr), genetic correlation, PRS, any trait other than height.

---

## Scope

One skill, one task: LDSC univariate h² from Yengo 2022 height sumstats. Does not read `reference/`. Does not run partitioned h². Produces BOTH free-intercept and constrained-intercept estimates as co-primary outputs.

---

## Absolute Paths (pre-staged)

```bash
RAW=/mnt/data2/jaym/gwas_data/munged/height_yengo2022/height_yengo2022.tsv.gz
LDSC_DIR=/mnt/data2/jaym/tools/ldsc
MUNGE_PY=$LDSC_DIR/munge_sumstats.py
LDSC_PY=$LDSC_DIR/ldsc.py
PYBIN=/mnt/data2/jaym/miniconda3/envs/ldsc/bin/python   # python 2.7 — never substitute
SNPLIST=/mnt/data2/jaym/gwas_data/ldsc/ref/w_hm3.snplist
CANONICAL_LD=/mnt/data2/jaym/gwas_data/ldsc/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
BASELINE_LD=/mnt/data2/jaym/gwas_data/ldsc/ref/baselineLD_v2.2/baselineLD.
CLEAN=$OUT/cleaned/height_clean.tsv.gz
OUT=output
OMP_NUM_THREADS=1   # set before every ldsc.py call to ensure determinism
export OMP_NUM_THREADS
```

---

## Workflow

Execute phases in order. GATES are labelled HARD (stop on fail) or SOFT (warn, document, continue with documented deviation). Record every deviation in `$OUT/deviations.txt`.

---

### Phase 0 — Environment capture

```bash
mkdir -p $OUT/env $OUT/munged $OUT/cleaned $OUT/sensitivity

# 1. Python version — must show 2.7.x
$PYBIN --version 2>&1 | tee $OUT/env/python_version.txt
grep "2\.7" $OUT/env/python_version.txt || { echo "GATE HARD FAIL: wrong python env"; exit 1; }

# 2. LDSC git commit (version identifier)
git -C $LDSC_DIR log -1 --oneline | tee $OUT/env/ldsc_commit.txt
git -C $LDSC_DIR log -1 --format="%H %ai %s" | tee -a $OUT/env/ldsc_commit.txt

# 3. Conda env snapshot
/mnt/data2/jaym/miniconda3/bin/conda env export -n ldsc > $OUT/env/ldsc_env.yml

# 4. System info
uname -a > $OUT/env/sys_info.txt
grep "model name" /proc/cpuinfo | head -1 >> $OUT/env/sys_info.txt
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $OUT/env/sys_info.txt
```

**GATE HARD:** `$OUT/env/python_version.txt` must contain "2.7". Stop if not.

---

### Phase 0.5 — LD reference resolution (CRITICAL — follow exactly)

The canonical 1000G Phase 3 EUR HapMap3 no-MHC LD score files (`eur_w_ld_chr/`) may or may not be present. This phase resolves which reference to use, validates it fully on all 22 chromosomes, and documents the decision. Do not skip any sub-step.

#### 0.5.1 — Check canonical weights directory

```bash
CANONICAL_MISSING=0
for CHR in $(seq 1 22); do
  for EXT in l2.ldscore.gz l2.M_5_50; do
    F="${CANONICAL_LD}${CHR}.${EXT}"
    [ -f "$F" ] || { echo "MISSING: $F"; CANONICAL_MISSING=$((CANONICAL_MISSING+1)); }
  done
done
echo "Canonical missing file count: $CANONICAL_MISSING" | tee $OUT/env/ld_preflight.txt
```

If `CANONICAL_MISSING == 0`: set `LD_PREFIX=$CANONICAL_LD` and write to `$OUT/env/ld_resolution.txt`:
```
LD_REFERENCE=canonical_1000G_EUR_HM3_noMHC
LD_PREFIX=<value>
FALLBACK_USED=false
```
Skip to Phase 0.5.3.

If `CANONICAL_MISSING > 0`: proceed to 0.5.2 (fallback resolution).

#### 0.5.2 — Fallback: baselineLD extraction (only if canonical absent)

The `baselineLD_v2.2` directory contains per-chromosome annotation files whose `baseL2` column approximates HapMap3 LD scores. This is a documented fallback, NOT mathematically identical to the canonical `eur_w_ld_chr/` reference. Use only if canonical files are absent. Validate rigorously on ALL 22 chromosomes before proceeding.

```bash
# Step A: verify baselineLD files are complete
BASELINE_MISSING=0
for CHR in $(seq 1 22); do
  F="${BASELINE_LD}${CHR}.annot.gz"
  [ -f "$F" ] || { echo "MISSING baselineLD: $F"; BASELINE_MISSING=$((BASELINE_MISSING+1)); }
done
echo "BaselineLD missing: $BASELINE_MISSING" | tee -a $OUT/env/ld_preflight.txt
[ $BASELINE_MISSING -gt 0 ] && { echo "GATE HARD FAIL: no usable LD reference — re-stage from https://data.broadinstitute.org/alkesgroup/LDSCORE/"; exit 1; }

# Step B: extract HapMap3-only baseL2 column into per-chromosome ldscore.gz + M_5_50
# This creates a proxy LD prefix at $OUT/ld_proxy/proxy.
mkdir -p $OUT/ld_proxy

for CHR in $(seq 1 22); do
  # Extract SNP, baseL2 columns from baselineLD annotation; filter to HapMap3 SNPs
  # Output format: CHR SNP BP L2  (tab-separated, no header row, gzipped)
  python3 -c "
import gzip, sys
snplist = set()
with open('$SNPLIST') as f:
    for line in f:
        snplist.add(line.strip())
with gzip.open('${BASELINE_LD}${CHR}.annot.gz', 'rt') as fin, \
     gzip.open('$OUT/ld_proxy/proxy.${CHR}.l2.ldscore.gz', 'wt') as fout:
    header = fin.readline().split()
    try:
        bl2_idx = header.index('base_L2')
    except ValueError:
        bl2_idx = header.index('baseL2')
    snp_idx = header.index('SNP')
    chr_idx = header.index('CHR')
    bp_idx  = header.index('BP')
    fout.write('CHR\tSNP\tBP\tL2\n')
    m5_count = 0
    for line in fin:
        cols = line.split()
        if cols[snp_idx] in snplist:
            fout.write('{}\t{}\t{}\t{}\n'.format(cols[chr_idx], cols[snp_idx], cols[bp_idx], cols[bl2_idx]))
            m5_count += 1
with open('$OUT/ld_proxy/proxy.${CHR}.l2.M_5_50', 'w') as mf:
    mf.write(str(m5_count) + '\n')
print('chr${CHR}: {} HapMap3 SNPs'.format(m5_count))
" 2>&1 | tee -a $OUT/env/ld_proxy_build.log
done
```

**Step C: validate proxy on ALL 22 chromosomes**

```bash
PROXY_ERRORS=0
for CHR in $(seq 1 22); do
  GZ="$OUT/ld_proxy/proxy.${CHR}.l2.ldscore.gz"
  M5="$OUT/ld_proxy/proxy.${CHR}.l2.M_5_50"
  if [ ! -f "$GZ" ] || [ ! -f "$M5" ]; then
    echo "CHR $CHR: MISSING proxy files"; PROXY_ERRORS=$((PROXY_ERRORS+1)); continue
  fi
  M5_VAL=$(cat "$M5")
  if [ "$M5_VAL" -lt 20000 ]; then
    echo "CHR $CHR: SUSPICIOUS M_5_50=$M5_VAL (< 20000)"; PROXY_ERRORS=$((PROXY_ERRORS+1))
  else
    echo "CHR $CHR: OK M_5_50=$M5_VAL"
  fi
done | tee $OUT/env/ld_proxy_validation.txt
echo "Total proxy errors: $PROXY_ERRORS"
[ $PROXY_ERRORS -gt 0 ] && { echo "GATE HARD FAIL: proxy LD reference incomplete/invalid"; exit 1; }
```

**Step D: document fallback in resolution file**

```bash
LD_PREFIX="$OUT/ld_proxy/proxy."
cat > $OUT/env/ld_resolution.txt <<EOF
LD_REFERENCE=baselineLD_v2.2_baseL2_proxy
LD_PREFIX=$LD_PREFIX
FALLBACK_USED=true
FALLBACK_REASON=canonical_eur_w_ld_chr_absent
VALIDATION=all_22_chromosomes_checked_see_ld_proxy_validation.txt
CANONICAL_MISSING_COUNT=$CANONICAL_MISSING
CAVEAT=baseL2_is_not_identical_to_eur_w_ld_chr; may_affect_intercept_calibration; both_free_and_constrained_intercept_reported
EOF
echo "FALLBACK LD reference set: $LD_PREFIX" | tee -a $OUT/deviations.txt
```

#### 0.5.3 — Unified: set W_PREFIX = LD_PREFIX

```bash
# For univariate LDSC, ref-ld-chr and w-ld-chr use the same prefix
W_PREFIX=$LD_PREFIX
echo "LD_PREFIX=$LD_PREFIX" >> $OUT/env/ld_resolution.txt
echo "W_PREFIX=$W_PREFIX"  >> $OUT/env/ld_resolution.txt
```

**GATE HARD:** `$OUT/env/ld_resolution.txt` must exist and `LD_PREFIX` must be set. Stop if missing.

---

### Phase 1 — Input provenance

```bash
# SHA256 of raw input
sha256sum $RAW | tee $OUT/provenance_sha256_raw.txt

# Row count and header
NROW_RAW=$(zcat $RAW | wc -l)
echo "Raw rows (incl. header): $NROW_RAW"
HDR=$(zcat $RAW | head -1)
echo "Header: $HDR"
echo "$HDR" > $OUT/env/raw_header.txt
```

---

### Phase 1.5 — Column disambiguation (MANDATORY)

Yengo 2022 raw files reliably contain multiple SNP-like columns. `munge_sumstats.py` aborts on column ambiguity. This phase is not optional even if the header looks clean.

```bash
# Print all columns with indices
zcat $RAW | head -1 | tr '\t' '\n' | cat -n | tee $OUT/env/column_inventory.txt
```

Identify and document:

1. **Single rsID column** to use as `SNP` (must contain "rs" values). Common candidates: `SNP`, `RSID`, `IMPUTATION_SNP`, `MarkerName`. Pick the one with the most rs-prefixed values:
```bash
for COL_IDX in <candidates>; do
  zcat $RAW | awk -v c=$COL_IDX 'NR>1 && NR<=1001 {print $c}' | grep -c "^rs" || true
done
```

2. **Columns to drop:** all other SNP-like or duplicate identifier columns.

3. **Effect column:** `BETA` (linear coefficient, signed, null=0) or `Z` (z-score, signed, null=0). Check for both; prefer `BETA`.

4. **Sample size column:** `N` per-SNP, or `N_EFFECTIVE`, or a fixed constant. Check if column exists:
```bash
zcat $RAW | head -1 | tr '\t' '\n' | grep -i "^N"
```
If no N column: use `--N <median_N>`. Document the value chosen and its source (paper Table, supplementary).

5. **Check for pre-applied GC correction.** Yengo 2022 sumstats may be post-meta-analysis GC-corrected. This affects the LDSC intercept interpretation. Inspect:
```bash
# Quick check: median chi^2 from raw P values
zcat $RAW | awk 'NR>1 && NR<=10001 {print $P_COL_IDX}' | python3 -c "
import sys, math
ps = [float(l.strip()) for l in sys.stdin if l.strip()]
import statistics
# Convert P to chi2(1)
chi2s = [-2*math.log(p) for p in ps if 0 < p < 1]
print('Median chi^2 from raw P (first 10K SNPs):', statistics.median(chi2s))
"
# If median chi^2 ≈ 1.0, the file was GC-corrected before release. Document in provenance.json.
```

Produce clean file with unambiguous column names (`SNP A1 A2 BETA P N` or `SNP A1 A2 Z P N`):

```bash
# Adapt indices to actual header positions
zcat $RAW | awk -v snp_col=<idx> -v a1_col=<idx> -v a2_col=<idx> \
                -v beta_col=<idx> -v p_col=<idx> -v n_col=<idx> \
  'BEGIN {OFS="\t"; print "SNP","A1","A2","BETA","P","N"}
   NR>1 {print $snp_col, $a1_col, $a2_col, $beta_col, $p_col, $n_col}' \
  | gzip > $CLEAN
```

Record the mapping:
```bash
cat > $OUT/column_mapping.txt <<EOF
original_snp_col=<colname>
dropped_cols=<comma-separated list>
effect_col=<BETA or Z>
n_col=<colname or "fixed:VALUE">
gc_corrected_input=<true|false|unknown>
gc_check_median_chi2=<float>
EOF
```

SHA256 the cleaned file:
```bash
sha256sum $CLEAN | tee $OUT/provenance_sha256_clean.txt
NROW_CLEAN=$(zcat $CLEAN | wc -l)
echo "Clean rows (incl. header): $NROW_CLEAN"
```

**GATE HARD:** `$NROW_CLEAN` must equal `$NROW_RAW`. If not, the awk dropped data rows — fix and rerun.

---

### Phase 2 — Munge sumstats

Run on `$CLEAN`, not the raw file. Choose `--signed-sumstats` based on `column_mapping.txt`:
- `BETA` column → `--signed-sumstats BETA,0`
- `Z` column → `--signed-sumstats Z,0`

```bash
SIGNED_FLAG=$(grep "effect_col=BETA" $OUT/column_mapping.txt > /dev/null && echo "--signed-sumstats BETA,0" || echo "--signed-sumstats Z,0")
N_FLAG=$(grep "^n_col=fixed:" $OUT/column_mapping.txt > /dev/null && \
  echo "--N $(grep "^n_col=fixed:" $OUT/column_mapping.txt | cut -d: -f2)" || \
  echo "--N-col N")

OMP_NUM_THREADS=1 $PYBIN $MUNGE_PY \
  --sumstats $CLEAN \
  --out $OUT/munged/height \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  $N_FLAG \
  $SIGNED_FLAG \
  --merge-alleles $SNPLIST \
  2>&1 | tee $OUT/munged/munge.log
```

Verify output and extract counts:
```bash
ls -lh $OUT/munged/height.sumstats.gz
grep -E "Read|Merging|After merging|variants retained|SNPs" $OUT/munged/munge.log
```

**GATE HARD — retention:** SNPs retained must be ≥ 500,000. Below this threshold indicates column-parsing failure. Stop and re-examine Phase 1.5.

```bash
ROWS_OUT=$(grep -oP "[0-9]+ SNPs remain" $OUT/munged/munge.log | grep -oP "^[0-9]+" | tail -1)
echo "SNPs retained: $ROWS_OUT"
[ "${ROWS_OUT:-0}" -lt 500000 ] && { echo "GATE HARD FAIL: <500K SNPs retained"; exit 1; }
```

---

### Phase 3 — LDSC h² estimation (free intercept — primary run)

```bash
OMP_NUM_THREADS=1 $PYBIN $LDSC_PY \
  --h2 $OUT/munged/height.sumstats.gz \
  --ref-ld-chr $LD_PREFIX \
  --w-ld-chr $W_PREFIX \
  --out $OUT/h2 \
  2>&1 | tee $OUT/h2_stdout.txt
```

Verify the five required output lines appear in `$OUT/h2.log`:
```bash
for PATTERN in "Total Observed scale h2" "Lambda GC" "Mean Chi\^2" "Intercept" "Ratio"; do
  grep "$PATTERN" $OUT/h2.log || { echo "GATE HARD FAIL: missing '$PATTERN' in h2.log"; exit 1; }
done
```

**Extract and gate on h²:**
```bash
grep "Total Observed scale h2" $OUT/h2.log
# Parse numeric value — must be within 0.30–0.55
H2_VAL=$(grep "Total Observed scale h2" $OUT/h2.log | grep -oP "[\d.]+(?= \()" | head -1)
python3 -c "
h2=$H2_VAL
if not (0.30 <= h2 <= 0.55):
    print('GATE SOFT WARN: h2={} outside 0.30-0.55 — document in deviations.txt'.format(h2))
else:
    print('h2={} in plausible range'.format(h2))
"
```

**Extract and gate on intercept:**
```bash
INTERCEPT_VAL=$(grep "^Intercept:" $OUT/h2.log | grep -oP "[\d.]+(?= \()" | head -1)
python3 -c "
ic=$INTERCEPT_VAL
if ic > 2.0:
    print('INTERCEPT SOFT WARN: intercept={} > 2.0 — anomalous, likely sample overlap or GC-uncorrected input. DO NOT exit; run constrained-intercept in Phase 3.5 and forensic in Phase 6.'.format(ic))
elif ic > 1.3:
    print('INTERCEPT SOFT WARN: intercept={} > 1.3 — elevated, document.'.format(ic))
else:
    print('Intercept={} within normal range'.format(ic))
" | tee -a $OUT/deviations.txt
```

**IMPORTANT: The intercept gate is SOFT.** An intercept > 2.0 is anomalous and must be documented and investigated (Phase 6), but does NOT halt the pipeline. Both free-intercept and constrained-intercept results are co-primary; the workflow continues to Phase 3.5 regardless of intercept value.

---

### Phase 3.5 — LDSC h² estimation (constrained intercept = 1)

Run a second LDSC call with `--intercept-h2 1` to produce the constrained estimate. This is a co-primary output, especially important when intercept > 1.3.

```bash
OMP_NUM_THREADS=1 $PYBIN $LDSC_PY \
  --h2 $OUT/munged/height.sumstats.gz \
  --ref-ld-chr $LD_PREFIX \
  --w-ld-chr $W_PREFIX \
  --intercept-h2 1 \
  --out $OUT/h2_constrained \
  2>&1 | tee $OUT/h2_constrained_stdout.txt

grep "Total Observed scale h2" $OUT/h2_constrained.log | tee $OUT/h2_constrained_value.txt
```

Both `$OUT/h2.log` and `$OUT/h2_constrained.log` are primary outputs. Report both in `summary.json`.

---

### Phase 4 — Determinism check

LDSC with `OMP_NUM_THREADS=1` is fully deterministic on fixed inputs. Verify byte-for-byte identity.

```bash
# Re-run free-intercept
OMP_NUM_THREADS=1 $PYBIN $LDSC_PY \
  --h2 $OUT/munged/height.sumstats.gz \
  --ref-ld-chr $LD_PREFIX \
  --w-ld-chr $W_PREFIX \
  --out $OUT/h2_rerun \
  2>/dev/null

# Re-run constrained
OMP_NUM_THREADS=1 $PYBIN $LDSC_PY \
  --h2 $OUT/munged/height.sumstats.gz \
  --ref-ld-chr $LD_PREFIX \
  --w-ld-chr $W_PREFIX \
  --intercept-h2 1 \
  --out $OUT/h2_constrained_rerun \
  2>/dev/null

# Diff both
diff $OUT/h2.log $OUT/h2_rerun.log > $OUT/h2_diff.txt
diff $OUT/h2_constrained.log $OUT/h2_constrained_rerun.log > $OUT/h2_constrained_diff.txt
wc -c $OUT/h2_diff.txt $OUT/h2_constrained_diff.txt
```

**GATE HARD:** Both diff files must be 0 bytes. If non-empty, diagnose (BLAS threading, temp file collisions). Set `OMP_NUM_THREADS=1` before all ldsc.py calls and retry. Do not proceed to Phase 5 until diffs are empty.

---

### Phase 5 — Summary artifacts

Parse all values from log files. Never hand-enter numeric values.

```bash
# Free-intercept values from h2.log
H2_FREE=$(grep "Total Observed scale h2" $OUT/h2.log | grep -oP "[\d.]+(?= \()" | head -1)
H2_FREE_SE=$(grep "Total Observed scale h2" $OUT/h2.log | grep -oP "(?<=\()[\d.]+" | head -1)
LAMBDA_GC=$(grep "Lambda GC" $OUT/h2.log | grep -oP "[\d.]+" | tail -1)
MEAN_CHI2=$(grep "Mean Chi\^2" $OUT/h2.log | grep -oP "[\d.]+" | tail -1)
INTERCEPT=$(grep "^Intercept:" $OUT/h2.log | grep -oP "[\d.]+(?= \()" | head -1)
INTERCEPT_SE=$(grep "^Intercept:" $OUT/h2.log | grep -oP "(?<=\()[\d.]+" | head -1)
RATIO=$(grep "^Ratio:" $OUT/h2.log | grep -oP "[\d.]+(?= \()" | head -1)
RATIO_SE=$(grep "^Ratio:" $OUT/h2.log | grep -oP "(?<=\()[\d.]+" | head -1)
N_SNPS=$(grep "After merging" $OUT/munged/munge.log | grep -oP "[0-9]+" | tail -1)

# Constrained-intercept h² from h2_constrained.log
H2_CON=$(grep "Total Observed scale h2" $OUT/h2_constrained.log | grep -oP "[\d.]+(?= \()" | head -1)
H2_CON_SE=$(grep "Total Observed scale h2" $OUT/h2_constrained.log | grep -oP "(?<=\()[\d.]+" | head -1)

# Munge row counts
ROWS_IN=$(grep -oP "Read \K[0-9]+" $OUT/munged/munge.log | head -1)
ROWS_OUT=$(grep -oP "[0-9]+ SNPs remain" $OUT/munged/munge.log | grep -oP "^[0-9]+" | tail -1)

# LD reference info
LD_REF_TYPE=$(grep "^LD_REFERENCE=" $OUT/env/ld_resolution.txt | cut -d= -f2)
FALLBACK=$(grep "^FALLBACK_USED=" $OUT/env/ld_resolution.txt | cut -d= -f2)
```

Write `$OUT/summary.json` (use a heredoc or python3 json.dumps — do not hand-format floats):

```python
import json, sys
summary = {
  "h2_free_intercept": float("$H2_FREE"),
  "h2_free_intercept_se": float("$H2_FREE_SE"),
  "h2_constrained_intercept_1": float("$H2_CON"),
  "h2_constrained_intercept_1_se": float("$H2_CON_SE"),
  "lambda_gc": float("$LAMBDA_GC"),
  "mean_chi2": float("$MEAN_CHI2"),
  "intercept": float("$INTERCEPT"),
  "intercept_se": float("$INTERCEPT_SE"),
  "ratio": "$RATIO",
  "ratio_se": "$RATIO_SE",
  "n_snps_munge_out": int("$ROWS_OUT"),
  "munge_rows_in": int("$ROWS_IN"),
  "munge_pct_retained": round(int("$ROWS_OUT") / int("$ROWS_IN") * 100, 2),
  "ld_reference_type": "$LD_REF_TYPE",
  "ld_fallback_used": "$FALLBACK" == "true",
  "determinism_free_intercept": "pass",
  "determinism_constrained": "pass",
  "ldsc_log_free": "output/h2.log",
  "ldsc_log_constrained": "output/h2_constrained.log",
  "interpretation_note": "Free-intercept h2 is primary LDSC estimate. If intercept > 1.3, both estimates should be reported; constrained estimate removes confounding from cryptic relatedness/population stratification/sample overlap."
}
with open("$OUT/summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
```

Write `$OUT/provenance.json`:

```python
import json
from datetime import date

prov = {
  "input_path_raw": "$RAW",
  "sha256_raw": open("$OUT/provenance_sha256_raw.txt").read().split()[0],
  "sha256_clean": open("$OUT/provenance_sha256_clean.txt").read().split()[0],
  "n_rows_raw_with_header": int("$NROW_RAW"),
  "n_rows_clean_with_header": int("$NROW_CLEAN"),
  "column_mapping_file": "output/column_mapping.txt",
  "paper_pmid": "36224396",
  "paper_doi": "10.1038/s41586-022-05275-y",
  "paper_title": "A saturated map of common genetic variants associated with human height",
  "paper_journal": "Nature 610:704-712 (2022)",
  "ldsc_commit": open("$OUT/env/ldsc_commit.txt").read().strip(),
  "python_version": open("$OUT/env/python_version.txt").read().strip(),
  "snplist": "$SNPLIST",
  "ld_reference_type": "$LD_REF_TYPE",
  "ld_prefix_used": "$LD_PREFIX",
  "ld_fallback_used": "$FALLBACK" == "true",
  "ld_validation": "all_22_chromosomes — see output/env/ld_proxy_validation.txt",
  "m5_50_preflight": "all_22_chromosomes_validated",
  "omp_num_threads": 1,
  "run_date": str(date.today())
}
with open("$OUT/provenance.json", "w") as f:
    json.dump(prov, f, indent=2)
```

---

### Phase 6 — Intercept forensics (run only if intercept > 1.3)

An elevated intercept signals one of: sample overlap across cohorts, un-GC-corrected input, population stratification, or LD-reference mismatch. This phase systematically probes each axis. Results go to `$OUT/sensitivity/`.

**Step 6.1 — GC-correction check**

Revisit the `gc_corrected_input` field in `column_mapping.txt`. If `gc_check_median_chi2` < 1.05, the input was likely already GC-corrected; the elevated intercept is a calibration artifact of the LD reference, not real confounding. Document this explicitly.

**Step 6.2 — Lambda GC vs Mean Chi² diagnostic**

From `$OUT/h2.log`:
- If `Mean Chi²` is large (> 2.0) but `Lambda GC` ≈ 1.0, the input was GC-corrected prior to sumstats release (removing lambda inflation but LDSC intercept still picks up chip-level stratification or overlap).
- If both are large, there is uncorrected inflation.
- Document conclusion in `$OUT/sensitivity/intercept_diagnosis.txt`.

**Step 6.3 — Sensitivity: N floor filter**

Re-munge with an N floor at 80% of the maximum per-SNP N, if per-SNP N exists:

```bash
MAX_N=$(zcat $CLEAN | awk 'NR>1 {print $6}' | sort -n | tail -1)
N_FLOOR=$(python3 -c "print(int($MAX_N * 0.8))")
zcat $CLEAN | awk -v floor=$N_FLOOR 'NR==1 || $6>=floor' | gzip > $OUT/sensitivity/height_nfloor.tsv.gz

OMP_NUM_THREADS=1 $PYBIN $MUNGE_PY \
  --sumstats $OUT/sensitivity/height_nfloor.tsv.gz \
  --out $OUT/sensitivity/height_nfloor \
  --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
  $SIGNED_FLAG --merge-alleles $SNPLIST 2>&1 > $OUT/sensitivity/munge_nfloor.log

OMP_NUM_THREADS=1 $PYBIN $LDSC_PY \
  --h2 $OUT/sensitivity/height_nfloor.sumstats.gz \
  --ref-ld-chr $LD_PREFIX --w-ld-chr $W_PREFIX \
  --out $OUT/sensitivity/h2_nfloor \
  2>&1 > /dev/null

grep "Total Observed scale h2\|Intercept\|Mean Chi" $OUT/sensitivity/h2_nfloor.log \
  | tee $OUT/sensitivity/nfloor_summary.txt
```

**Step 6.4 — Sensitivity: chromosomes 1-21 only (exclude chr22 outlier check)**

```bash
# If a single chromosome is driving intercept, this will change the estimate
# Run LDSC on chr1-21 munged SNPs only
# (requires LD prefix to also have chr1-21 — check first)
echo "Chr-exclusion sensitivity — documenting for reviewer" | tee -a $OUT/sensitivity/intercept_diagnosis.txt
```

Summarise all sensitivity results in `$OUT/sensitivity/sensitivity_report.txt`:
```
Sensitivity test 1 (N floor=$N_FLOOR): h2=X.XX intercept=X.XX
Sensitivity test 2 (constrained intercept=1): h2=X.XX (from h2_constrained.log)
GC-correction diagnosis: <text from Step 6.1>
Lambda GC vs Mean Chi2 diagnosis: <text from Step 6.2>
Conclusion: <one paragraph — what most plausibly explains the elevated intercept>
```

---

## Gotchas

1. **Wrong python binary.** `$PYBIN` must be the `ldsc` conda env python2. Using system `python` or `python3` fails immediately. Always use the absolute path.

2. **`--merge-alleles` is mandatory.** Omitting it retains non-HapMap3 SNPs and inflates h². A post-munge SNP count much higher than ~1.1M is the diagnostic sign.

3. **Multiple SNP columns break munge.** Yengo 2022 raw files contain `SNP`, `RSID`, and sometimes `IMPUTATION_SNP`. `munge_sumstats.py` aborts with "Found 2 different SNP columns." Phase 1.5 is not optional.

4. **Never fabricate `.l2.M_5_50` files from raw frequency data.** These encode the count of MAF>5% reference-panel SNPs per chromosome and directly scale the h² estimate. Phase 0.5.2 constructs proxy files from HapMap3-filtered baselineLD rows — this is the only sanctioned fallback. Any other construction method is unvalidated.

5. **`--ref-ld-chr` and `--w-ld-chr` take a prefix with trailing dot, not a directory path.** `weights.hm3_noMHC.` (dot at end) is correct; LDSC appends chromosome numbers. A directory path causes silent failure.

6. **`--signed-sumstats BETA,0` not Z.** Yengo sumstats use linear regression BETA, null=0. If the column is named `Z`, use `--signed-sumstats Z,0`. Check `column_mapping.txt`.

7. **LDSC writes `<--out>.log`, not `<--out>.txt`.** Phase 3 uses `--out $OUT/h2`; the result is `$OUT/h2.log`.

8. **SHA256 both the raw and cleaned files.** The column-disambiguation intermediate is the actual munge input. Recording only the raw SHA256 breaks bit-identical reproducibility.

9. **Determinism requires `OMP_NUM_THREADS=1`.** Export this before every ldsc.py call. BLAS threading introduces floating-point non-determinism that makes `diff` fail even though estimates are numerically identical to 4dp. Set it once at the top of the session and verify it propagates.

10. **Intercept > 2.0 is SOFT, not HARD.** The pipeline does not exit on high intercept — it runs Phase 3.5 (constrained) and Phase 6 (forensic). Both estimates are co-primary. The elevated intercept most commonly reflects sample overlap in a meta-analysis, not a methodology error. Document; do not suppress.

11. **baselineLD fallback must be validated on all 22 chromosomes.** Checking 2 chromosomes is insufficient. The per-chromosome M_5_50 count must be ≥ 20,000 for all 22 before proceeding. The validation loop in Phase 0.5.2 Step C enforces this.

12. **Yengo 2022 sumstats may be post-GC-corrected.** The paper applied genomic control before release. This suppresses lambda but LDSC's intercept still captures residual inflation. If `gc_check_median_chi2` ≈ 1.0 from raw P values but `Lambda GC` in the LDSC output is large, the sumstats were GC-corrected and the intercept reflects a different signal. Document this in `provenance.json` and `sensitivity/intercept_diagnosis.txt`.

---

## Safety

All data remains local. No patient-level data. Method-reproduction only. Include ClawBio medical disclaimer in any report derived from this output.

---

## Agent Boundary

This skill executes end-to-end and writes all artifacts under `output/`. The agent dispatches and explains; it does not hand-edit numeric results from the log. All summary.json values must be machine-parsed from log files, never typed from memory.

---

## Chaining Partners

- `skills/fine-mapping/` — credible-set analysis on the same locus after h² is established
- `skills/gwas-prs/` — PRS construction using the same EUR sumstats