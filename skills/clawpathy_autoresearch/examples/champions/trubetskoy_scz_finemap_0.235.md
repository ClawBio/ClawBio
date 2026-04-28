---
name: repro-trubetskoy-scz
version: 1.3.0
author: ClawBio
description: Reproduce Trubetskoy et al. 2022 PGC3 SCZ fine-mapping — FINEMAP v1.4, 1000G phase 3 EUR LD, credible sets with PIPs, 10-locus tuning subset, scalable to all 287.
inputs:
  - name: sumstats
    format: "TSV daner format, GRCh37 (PGC3 SCZ public release)"
    required: true
  - name: loci_list
    format: "TSV from paper Supplementary Table 3 (chr, start, end, locus_id, lead_rsid)"
    required: true
  - name: ld_reference
    format: "1000G phase 3 per-chromosome VCF.gz + .tbi (downloaded in Phase 0)"
    required: false
outputs:
  - name: credible_sets
    format: "output/credible_sets/<locus_id>_credible_sets.tsv (per-variant)"
  - name: provenance
    format: "output/provenance.json, output/provenance/*.txt"
  - name: environment
    format: "output/env/tool_versions.txt, output/env/environment.yml"
trigger_keywords:
  - "Trubetskoy"
  - "PGC3 SCZ fine-mapping"
  - "FINEMAP SCZ"
  - "schizophrenia credible sets"
  - "reproduce Trubetskoy 2022"
---

## Trigger

**Fire when:** user asks to reproduce Trubetskoy et al. 2022 fine-mapping; build FINEMAP v1.4 credible sets on PGC3 SCZ; replicate PGC3 SCZ PIPs.

**Do NOT fire when:** generic SuSiE fine-mapping (→ `skills/fine-mapping/`); other GWAS traits; variant annotation only (→ `skills/vcf-annotator/`).

---

## Scope

One skill, one task: FINEMAP v1.4 fine-mapping of PGC3 SCZ GWAS loci with 1000G EUR LD, producing per-locus credible sets with PIPs. Does not read `reference/` under any circumstances.

---

## Four Decisions That Must Be Made Before Writing Code

Run these four checks and record results before Phase 1. They determine branching points later.

### Decision 1 — Does FINEMAP v1.4 have `--seed`?

```bash
$FINEMAP --help 2>&1 | grep -i seed \
  && echo "SEED_FLAG_PRESENT" \
  || echo "NO_SEED_FLAG"
```

**If NO_SEED_FLAG** (expected for v1.4): use 5-replicate median strategy throughout. Record in `output/provenance/determinism_note.txt`. Rerun verdicts are exactly **PASS / WARN / FAIL** — no other labels. PASS: corr ≥ 0.999 AND max_diff ≤ 0.01. WARN: corr ≥ 0.95. FAIL: corr < 0.95 → exclude locus, do not relabel.

### Decision 2 — What is Neff?

```python
# scripts/derive_neff.py  (run as: python scripts/derive_neff.py raw_sumstats.tsv)
import pandas as pd, numpy as np, pathlib, sys
df = pd.read_csv(sys.argv[1], sep="\t", nrows=200000)
nca = next((c for c in df.columns if c.lower() in ["nca","ncas","ncases"]), None)
nco = next((c for c in df.columns if c.lower() in ["nco","ncon","ncontrols"]), None)
if nca and nco:
    neff = float(np.median(4 / (1/df[nca].astype(float) + 1/df[nco].astype(float))))
    method = "4/(1/Nca+1/Nco) per-variant median"
elif "Neff_half" in df.columns:
    neff = float(df["Neff_half"].median()) * 2; method = "Neff_half*2"
else:
    n_col = next((c for c in df.columns if c.upper() in ["N","NEFF","N_EFF"]), None)
    neff = float(df[n_col].median()) if n_col else 126000
    method = f"column_{n_col}" if n_col else "paper_fallback"
pathlib.Path("output/provenance").mkdir(parents=True, exist_ok=True)
pathlib.Path("output/provenance/neff_derivation.txt").write_text(
    f"neff_used={neff:.0f}\nderivation_method={method}\n"
    f"paper_full_meta_neff=126000\n"
    f"public_subset_limitation={'YES - public file is subset; PIPs lower, CS wider than paper' if neff < 100000 else 'No'}\n"
)
print(f"Neff={neff:.0f} via {method}. Paper full meta: ~126000.")
```

### Decision 3 — Confirm `prior_std = 0.0224` everywhere

Benner 2016 (PLoS Genet 12:e1005942 Table 1) binary-trait default for case-control log(OR). **Never 0.05.** Every FINEMAP invocation must pass `--prior-std 0.0224`.

### Decision 4 — Confirm PLINK binary is v1.9

```bash
plink --version 2>&1 | grep -q "1\.9" \
  || { echo "ERROR: plink 1.9 required — found: $(plink --version 2>&1 | head -1)"; exit 1; }
```

Store the confirmed PLINK path in `output/config.sh`. Do not mix plink2 pgen commands with plink1.9 r-matrix commands — they are format-incompatible.

---

## Workflow

### Phase 0 — Environment (run once)

**Step 0.1 — Download FINEMAP v1.4 and record all tool versions to STDOUT and file.**

```bash
mkdir -p output/env output/provenance output/credible_sets output/loci output/sensitivity scripts data/1kg_eur

# FINEMAP v1.4
wget http://www.christianbenner.com/finemap_v1.4_x86_64.tgz -O output/env/finemap_v1.4_x86_64.tgz
sha256sum output/env/finemap_v1.4_x86_64.tgz | tee output/provenance/finemap_sha256.txt
tar -xzf output/env/finemap_v1.4_x86_64.tgz -C output/env/
FINEMAP=output/env/finemap_v1.4_x86_64/finemap_v1.4_x86_64
chmod +x $FINEMAP
sha256sum $FINEMAP | tee -a output/provenance/finemap_sha256.txt

# Version strings — MUST appear in transcript
echo "=== TOOL VERSIONS ===" | tee output/env/tool_versions.txt
$FINEMAP --version 2>&1 | tee -a output/env/tool_versions.txt
plink --version 2>&1 | head -1 | tee -a output/env/tool_versions.txt
bcftools --version 2>&1 | head -1 | tee -a output/env/tool_versions.txt
python --version 2>&1 | tee -a output/env/tool_versions.txt
uname -srm | tee -a output/env/tool_versions.txt

# Fail loudly if not plink 1.9
grep -q "1\.9" output/env/tool_versions.txt \
  || { echo "ABORT: plink 1.9 not found in tool_versions.txt"; exit 1; }

# Record seed status
$FINEMAP --help 2>&1 | grep -iq seed \
  && echo "SEED_FLAG_PRESENT" \
  || echo "NO_SEED_FLAG: 5-replicate median strategy in use" \
  | tee output/provenance/determinism_note.txt

pip freeze > output/env/pip_freeze.txt
conda env export > output/env/environment.yml 2>/dev/null \
  || echo "conda not active; pip_freeze.txt is environment record"
```

**Step 0.2 — Write config.**

```bash
python scripts/derive_neff.py raw_sumstats.tsv

cat > output/config.sh <<'EOF'
export FINEMAP=output/env/finemap_v1.4_x86_64/finemap_v1.4_x86_64
export PLINK=plink
export NEFF=$(grep neff_used output/provenance/neff_derivation.txt | cut -d= -f2)
export PRIOR_STD=0.0224
export N_CAUSAL_MAX=5
export N_REPLICATES=5
export CRED_COVERAGE=0.95
EOF
source output/config.sh
echo "Config: NEFF=${NEFF} PRIOR_STD=${PRIOR_STD} N_CAUSAL_MAX=${N_CAUSAL_MAX} N_REPLICATES=${N_REPLICATES}"
```

**Step 0.3 — Record input checksums.**

```python
# scripts/record_provenance.py <sumstats_path> <retrieval_date_ISO8601> <loci_tsv_path>
import hashlib, json, datetime, sys, pathlib

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

arts = {
    "pgc3_scz_sumstats": {
        "path": sys.argv[1], "sha256": sha256(sys.argv[1]),
        "source_url": "https://pgc.unc.edu/for-researchers/download-results/",
        "retrieval_date": sys.argv[2], "build": "GRCh37",
        "description": "PGC3 SCZ GWAS summary statistics EUR meta-analysis (public release)",
    },
    "loci_list": {
        "path": sys.argv[3], "sha256": sha256(sys.argv[3]),
        "description": "Supplementary Table 3, Trubetskoy et al. 2022 (Nature 604:502-508, DOI 10.1038/s41586-022-04434-5)",
    },
    "finemap_binary": {
        "path": "output/env/finemap_v1.4_x86_64/finemap_v1.4_x86_64",
        "sha256": sha256("output/env/finemap_v1.4_x86_64/finemap_v1.4_x86_64"),
        "source_url": "http://www.christianbenner.com/finemap_v1.4_x86_64.tgz",
        "version": "1.4",
    },
    "record_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
}
json.dump(arts, open("output/provenance.json", "w"), indent=2)
print("Provenance written.")
```

Run: `python scripts/record_provenance.py raw_sumstats.tsv <ISO8601_date> data/Supplementary_Table3.tsv`

**CHECKPOINT 0 — Verify before Phase 1:**

```bash
for f in output/provenance.json output/provenance/neff_derivation.txt \
          output/provenance/determinism_note.txt output/env/tool_versions.txt \
          output/config.sh; do
  [ -f "$f" ] && echo "OK: $f" || echo "MISSING: $f — stop and fix"
done
grep -q "1\.9" output/env/tool_versions.txt && echo "PLINK v1.9 confirmed" || echo "ABORT: plink 1.9 not confirmed"
```

---

### Phase 1 — Sumstats QC

**Step 1.1 — QC and standardize.**

Inspect the header first: `head -1 raw_sumstats.tsv`. PGC3 daner format columns include `CHR BP SNP A1 A2 FRQ_A_<N> FRQ_U_<N> INFO OR SE P Nca Nco Neff_half`. Column names vary by release; the script normalizes them.

```python
# scripts/qc_sumstats.py <infile> <outfile>
import pandas as pd, numpy as np, re, pathlib, sys

infile, outfile = sys.argv[1], sys.argv[2]
df = pd.read_csv(infile, sep="\t", low_memory=False)
n0 = len(df)
print(f"Input: {n0} variants. Columns: {list(df.columns)}")

# Normalize column names
cm = {}
for c in df.columns:
    cl = c.lower()
    if cl in ["chr","chrom","chromosome"]: cm[c]="CHROM"
    elif cl in ["bp","pos","position"]: cm[c]="POS"
    elif cl in ["snp","rsid","markername"]: cm[c]="SNP"
    elif cl=="a1": cm[c]="EA"
    elif cl=="a2": cm[c]="NEA"
    elif re.match(r"frq_a",cl) or cl in ["eaf","freq","frq"]: cm[c]="EAF"
    elif cl=="info": cm[c]="INFO"
    elif cl in ["or","odds_ratio"]: cm[c]="OR"
    elif cl in ["beta","effect"]: cm[c]="BETA"
    elif cl in ["se","stderr"]: cm[c]="SE"
    elif cl=="p" or cl.startswith("p_"): cm[c]="P"
df.rename(columns=cm, inplace=True)

# BETA from OR — PGC3 daner uses OR, not BETA
if "BETA" not in df.columns:
    if "OR" in df.columns:
        df["BETA"] = np.log(df["OR"].astype(float))
        print("BETA = log(OR)")
    else:
        raise SystemExit("ERROR: no OR or BETA column")

# Effective N
if "Nca" in df.columns and "Nco" in df.columns:
    df["N"] = 4 / (1/df["Nca"].astype(float) + 1/df["Nco"].astype(float))
elif "Neff_half" in df.columns:
    df["N"] = df["Neff_half"].astype(float) * 2
else:
    df["N"] = np.nan

# Filters (log counts at each step)
df["CHROM"] = pd.to_numeric(df["CHROM"], errors="coerce")
df = df[df["CHROM"].between(1,22)]; print(f"After autosome: {len(df)}")
if "INFO" in df.columns:
    n=len(df); df=df[df["INFO"].astype(float)>=0.6]; print(f"INFO<0.6 dropped: {n-len(df)}")
if "EAF" in df.columns:
    maf=df["EAF"].astype(float).clip(0,1); maf=maf.where(maf<=0.5,1-maf)
    n=len(df); df=df[maf>=0.01]; print(f"MAF<0.01 dropped: {n-len(df)}")
df["SE"]=pd.to_numeric(df["SE"],errors="coerce"); df["BETA"]=pd.to_numeric(df["BETA"],errors="coerce")
df=df[df["SE"].notna()&(df["SE"]>0)&df["BETA"].notna()]
df["Z"]=df["BETA"]/df["SE"]
n=len(df); df["CHROM"]=df["CHROM"].astype(int); df["POS"]=df["POS"].astype(int)
df=df.drop_duplicates(subset=["CHROM","POS","EA","NEA"]); print(f"Duplicates: {n-len(df)}")

lambda_gc = float(np.median(df["Z"]**2)/0.4549)
print(f"Lambda_GC: {lambda_gc:.4f}")
if lambda_gc > 2.5: print("WARNING: Lambda_GC > 2.5")

pathlib.Path("output/provenance/qc_metrics.txt").write_text(
    f"lambda_gc={lambda_gc:.6f}\nn_input={n0}\nn_output={len(df)}\n"
)
df.to_csv(outfile, sep="\t", index=False)
print(f"QC: {n0} → {len(df)}")
```

Run: `python scripts/qc_sumstats.py raw_sumstats.tsv output/sumstats_qc.tsv`

**Step 1.2 — LDSC intercept (best effort; record outcome either way).**

```bash
mkdir -p output/provenance/ldsc
if command -v ldsc.py &>/dev/null; then
    munge_sumstats.py --sumstats output/sumstats_qc.tsv --N-col N \
        --snp SNP --a1 EA --a2 NEA --p P --signed-sumstats Z,0 \
        --out output/provenance/ldsc/pgc3_scz \
        --merge-alleles data/w_hm3.snplist 2>&1 | tee output/provenance/ldsc/munge.log
    ldsc.py --h2 output/provenance/ldsc/pgc3_scz.sumstats.gz \
        --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ \
        --out output/provenance/ldsc/h2 2>&1 | tee output/provenance/ldsc/ldsc.log
    grep "Intercept" output/provenance/ldsc/h2.log \
        | tee output/provenance/ldsc/ldsc_intercept.txt
else
    echo "LDSC UNAVAILABLE: ldsc not installed; lambda_GC=$(grep lambda_gc output/provenance/qc_metrics.txt)" \
        | tee output/provenance/ldsc/ldsc_intercept.txt
fi
# Gate: file must exist (content "UNAVAILABLE" is acceptable)
[ -f output/provenance/ldsc/ldsc_intercept.txt ] && echo "LDSC gate OK" || { echo "ABORT: ldsc gate"; exit 1; }
```

---

### Phase 2 — Locus Definition and Tuning Subset

**Step 2.1 — Parse Supplementary Table 3.**

```python
# scripts/parse_st3.py  →  output/loci_raw.tsv
import pandas as pd, sys, pathlib

try:
    df = pd.read_excel("data/Supplementary_Table3.xlsx", header=0)
except Exception:
    df = pd.read_csv("data/Supplementary_Table3.tsv", sep="\t")
print("Columns:", list(df.columns))

# Map whatever column names the file uses
cm = {}
for c in df.columns:
    cl = c.lower().strip()
    if "locus" in cl and "id" in cl: cm[c]="locus_id"
    elif cl in ["chr","chrom","chromosome"]: cm[c]="chr"
    elif cl in ["start","bp_start","pos_start"]: cm[c]="start"
    elif cl in ["end","bp_end","pos_end"]: cm[c]="end"
    elif "lead" in cl and ("rsid" in cl or "snp" in cl): cm[c]="lead_rsid"
df.rename(columns=cm, inplace=True)

missing = [c for c in ["locus_id","chr","start","end"] if c not in df.columns]
if missing: raise SystemExit(f"ERROR: missing columns {missing}. Inspect file headers and update cm.")

df = df[["locus_id","chr","start","end","lead_rsid"]].dropna(subset=["chr","start","end"])
df["chr"] = df["chr"].astype(str).str.replace("chr","",regex=False)
df = df[~df["chr"].str.upper().isin(["X","Y","MT"])]
df["chr"]=df["chr"].astype(int); df["start"]=df["start"].astype(int); df["end"]=df["end"].astype(int)

# Exclude MHC (chr6:25–35 Mb) — LD structure is too complex for standard fine-mapping
n=len(df)
df = df[~((df["chr"]==6)&(df["end"]>=25000000)&(df["start"]<=35000000))]
print(f"MHC excluded: {n-len(df)}. Total non-MHC autosomal loci: {len(df)}")
df.to_csv("output/loci_raw.tsv", sep="\t", index=False)
```

**Step 2.2 — Count sumstats variants per locus.**

```python
# scripts/count_locus_variants.py  output/sumstats_qc.tsv  output/loci_raw.tsv  output/loci_with_counts.tsv
import pandas as pd, sys
ss = pd.read_csv(sys.argv[1], sep="\t")
loci = pd.read_csv(sys.argv[2], sep="\t")
loci["n_variants"] = [
    int(((ss["CHROM"]==int(r["chr"]))&(ss["POS"]>=int(r["start"]))&(ss["POS"]<=int(r["end"]))).sum())
    for _, r in loci.iterrows()
]
loci.to_csv(sys.argv[3], sep="\t", index=False)
print(f"{len(loci)} loci written. n_variants summary:\n{loci['n_variants'].describe().to_string()}")
```

**Step 2.3 — Select 10-locus tuning subset (stratified by n_variants, fixed seed=42).**

```python
# scripts/select_tuning_subset.py  output/loci_with_counts.tsv  output/tuning_subset.tsv
import pandas as pd, numpy as np, sys

SEED = 42; N = 10
loci = pd.read_csv(sys.argv[1], sep="\t")
loci = loci[loci["n_variants"]>=10].sort_values("n_variants").reset_index(drop=True)

q33, q67 = np.percentile(loci["n_variants"], [33.3, 66.7])
loci["tertile"] = pd.cut(loci["n_variants"], bins=[-1,q33,q67,float("inf")],
                          labels=["small","medium","large"])
print(f"Tertile bounds: small≤{q33:.0f}, medium≤{q67:.0f}, large>{q67:.0f}")

chosen = pd.concat([
    loci[loci["tertile"]==t].sample(n=min(k,len(loci[loci["tertile"]==t])), random_state=SEED)
    for t, k in [("small",3),("medium",3),("large",4)]
]).sort_values("chr").reset_index(drop=True)

print(f"Tuning subset ({len(chosen)} loci):")
print(chosen[["locus_id","chr","start","end","n_variants","tertile"]].to_string())
chosen.to_csv(sys.argv[2], sep="\t", index=False)
# Gate
assert len(chosen) == N, f"Expected 10 loci, got {len(chosen)}"
```

Run: `python scripts/select_tuning_subset.py output/loci_with_counts.tsv output/tuning_subset.tsv`

**GATE 2:** `output/tuning_subset.tsv` must have exactly 10 rows. If a locus has < 10 sumstats variants, log `SKIP` in `output/run_status.txt` and continue with N-1. **Do not substitute.**

---

### Phase 3 — LD Reference (build once, cached per locus)

**Step 3.1 — Download 1000G phase 3 and extract EUR samples.**

```bash
# scripts/download_1kg.sh
set -euo pipefail
BASE=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
for CHR in $(seq 1 22); do
    F="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    wget -c "${BASE}/${F}" -O "data/1kg_eur/${F}"
    wget -c "${BASE}/${F}.tbi" -O "data/1kg_eur/${F}.tbi"
    sha256sum "data/1kg_eur/${F}" >> output/provenance/1kg_checksums.txt
done
PANEL=data/1kg_eur/integrated_call_samples_v3.20130502.ALL.panel
wget -c "${BASE}/$(basename $PANEL)" -O $PANEL
sha256sum $PANEL >> output/provenance/1kg_checksums.txt

# EUR superpopulations: CEU, TSI, FIN, GBR, IBS
awk '$3=="EUR"{print $1}' $PANEL > data/1kg_eur/EUR_samples.txt
EUR_N=$(wc -l < data/1kg_eur/EUR_samples.txt)
echo "EUR samples: ${EUR_N} (expected 503)" | tee -a output/provenance/ld_reference_log.txt
[ "${EUR_N}" -eq 503 ] || echo "WARNING: expected 503, got ${EUR_N}"
awk '{print $1, $1}' data/1kg_eur/EUR_samples.txt > data/1kg_eur/EUR_keep.txt
echo "EUR subpopulations included: CEU TSI FIN GBR IBS" >> output/provenance/ld_reference_log.txt
```

**Step 3.2 — Per-locus LD matrix (cached; Pearson r, NOT r²).**

```bash
# scripts/compute_ld_locus.sh <locus_id> <chr> <start> <end>
set -euo pipefail
source output/config.sh
LOCUS_ID=$1; CHR=$2; START=$3; END=$4
OUTDIR="output/loci/${LOCUS_ID}"
mkdir -p "${OUTDIR}"

# Cache: skip if both files exist and are non-empty
if [ -s "${OUTDIR}/ld.ld" ] && [ -s "${OUTDIR}/ld.snplist" ]; then
    echo "LD cached for ${LOCUS_ID} ($(wc -l < "${OUTDIR}/ld.snplist") variants)"
    exit 0
fi

VCF="data/1kg_eur/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# VCF → BED (plink 1.9; filters: biallelic SNPs, MAF≥1%, missingness≤5%)
$PLINK \
    --vcf "${VCF}" \
    --keep data/1kg_eur/EUR_keep.txt \
    --chr "${CHR}" --from-bp "${START}" --to-bp "${END}" \
    --snps-only just-acgt --biallelic-only \
    --maf 0.01 --geno 0.05 \
    --make-bed --out "${OUTDIR}/ld_plink" \
    2>&1 | tee "${OUTDIR}/plink_extract.log"

N_VAR=$(wc -l < "${OUTDIR}/ld_plink.bim" 2>/dev/null || echo 0)
echo "Extracted ${N_VAR} variants for ${LOCUS_ID}"

if [ "${N_VAR}" -lt 10 ]; then
    echo "SKIP_LD: ${N_VAR} variants" | tee "${OUTDIR}/status.txt"; exit 0
fi

# Pearson r matrix (FINEMAP v1.4 requires r, NOT r²)
$PLINK --bfile "${OUTDIR}/ld_plink" --r square --out "${OUTDIR}/ld" \
    2>&1 | tee -a "${OUTDIR}/plink_extract.log"

awk '{print $2}' "${OUTDIR}/ld_plink.bim" > "${OUTDIR}/ld.snplist"
echo "LD_N=$(wc -l < "${OUTDIR}/ld.snplist") FORMAT=pearson_r" | tee -a "${OUTDIR}/status.txt"

# LD spot-check (diagonal = 1.0; print to STDOUT for first locus processed)
python3 - "${OUTDIR}" <<'PYEOF'
import sys, numpy as np, pandas as pd, pathlib
d = pathlib.Path(sys.argv[1])
ld = pd.read_csv(d/"ld.ld", sep="\t", header=None).values
snps = pd.read_csv(d/"ld.snplist", header=None, names=["rsid"])
print(f"LD shape: {ld.shape}. Diagonal range: {np.diag(ld).min():.4f}–{np.diag(ld).max():.4f}")
print("Top-left 4x4 (Pearson r):")
print(pd.DataFrame(ld[:4,:4], index=snps["rsid"][:4], columns=snps["rsid"][:4]).round(4).to_string())
if not np.allclose(np.diag(ld), 1.0, atol=0.02):
    print("WARNING: diagonal ≠ 1.0 — LD matrix may be wrong format or corrupted")
PYEOF
```

---

### Phase 4 — Allele Alignment and Z-file Construction

**Step 4.1 — Harmonize sumstats alleles against 1000G LD panel per locus.**

PGC3 convention: A1 = effect allele. 1000G PLINK BIM A1 = minor allele (arbitrary). Expect ~50% direct, ~50% flipped across a locus. **If n_direct == 0 across all aligned variants, abort the locus — this is a systematic encoding error.**

```python
# scripts/harmonize_locus.py <locus_id> <chr> <start> <end> <sumstats_path>
import pandas as pd, numpy as np, sys, pathlib, json

locus_id, chr_, start, end, ss_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
outdir = pathlib.Path(f"output/loci/{locus_id}")
status_path = outdir / "status.txt"

neff = float(next(l.split("=")[1] for l in open("output/provenance/neff_derivation.txt") if l.startswith("neff_used")))

ss = pd.read_csv(ss_path, sep="\t")
lss = ss[(ss["CHROM"]==int(chr_))&(ss["POS"]>=start)&(ss["POS"]<=end)].copy()
if len(lss) < 10:
    status_path.write_text(f"SKIP_HARMONIZE: {len(lss)} sumstats variants\n"); sys.exit(0)

bim = pd.read_csv(outdir/"ld_plink.bim", sep="\t", header=None,
                  names=["bim_chr","rsid","cm","pos","a1","a2"])
ld_snps = pd.read_csv(outdir/"ld.snplist", header=None, names=["rsid"])

def comp(s): return s.translate(str.maketrans("ACGTacgt","TGCAtgca"))
def ambig(a,b): return {a.upper(),b.upper()} in [{"A","T"},{"C","G"}]

merged = lss.merge(bim[["rsid","pos","a1","a2"]], left_on="POS", right_on="pos", how="inner")
records, n_dir=[], n_flip=[], n_strand_dir=[], n_strand_flip=[], n_ambig=0, n_mismatch=0
# (counters initialized as ints)
n_dir=n_flip=n_sd=n_sf=n_ambig=n_mm=0
recs=[]
for _,row in merged.iterrows():
    ea=str(row["EA"]).upper(); nea=str(row["NEA"]).upper()
    la1=str(row["a1"]).upper(); la2=str(row["a2"]).upper()
    if ambig(ea,nea): n_ambig+=1; continue
    z,eaf,st=np.nan,np.nan,"mismatch"
    if ea==la1 and nea==la2:   z,eaf,st=row["Z"],row["EAF"],"direct";        n_dir+=1
    elif ea==la2 and nea==la1: z,eaf,st=-row["Z"],1-row["EAF"],"flipped";    n_flip+=1
    else:
        cea,cnea=comp(ea),comp(nea)
        if cea==la1 and cnea==la2:   z,eaf,st=row["Z"],row["EAF"],"strand";      n_sd+=1
        elif cea==la2 and cnea==la1: z,eaf,st=-row["Z"],1-row["EAF"],"strand_flip"; n_sf+=1
        else: n_mm+=1; continue
    recs.append({"rsid":row["rsid"],"POS":int(row["POS"]),"a1":la1,"a2":la2,
                 "maf":float(min(eaf,1-eaf)) if not np.isnan(eaf) else np.nan,
                 "BETA":float(row["BETA"]),"SE":float(row["SE"]),"Z_aligned":float(z),
                 "N":float(row["N"]) if "N" in row.index and not np.isnan(row["N"]) else neff,
                 "original_rsid":str(row.get("SNP",f"chr{chr_}:{int(row['POS'])}"))})

n_total=n_dir+n_flip+n_sd+n_sf
log={"locus_id":locus_id,"n_sumstats_in_window":len(lss),"n_ambiguous_excluded":n_ambig,
     "n_mismatch_excluded":n_mm,"n_direct":n_dir,"n_flipped":n_flip,
     "n_strand":n_sd,"n_strand_flip":n_sf,"n_aligned_total":n_total,
     "flip_rate_pct":round(100*(n_flip+n_sf)/max(1,n_total),1)}
json.dump(log,open(outdir/"harmonization_summary.json","w"),indent=2)
print(json.dumps(log,indent=2))

# HARD ABORT: n_direct == 0 is a systematic allele encoding error
if n_dir==0 and n_total>0:
    msg=(f"ABORT_ALLELE_ERROR {locus_id}: n_direct=0 across {n_total} aligned variants. "
         f"Check: (1) correct EA column, (2) log(OR) applied, (3) 1000G BIM A1 convention.")
    with open(status_path,"a") as f: f.write(msg+"\n")
    raise SystemExit(msg)

aligned=pd.DataFrame(recs)
aligned=aligned[aligned["rsid"].isin(ld_snps["rsid"])]
ld_snps["order"]=range(len(ld_snps))
aligned=aligned.merge(ld_snps,on="rsid",how="inner").sort_values("order").drop("order",axis=1)
if len(aligned)<10:
    with open(status_path,"a") as f: f.write(f"SKIP_FINEMAP: {len(aligned)} aligned variants\n")
    print(f"SKIP {locus_id}: {len(aligned)} aligned"); sys.exit(0)

# Write FINEMAP z-file (space-separated; columns must match FINEMAP v1.4 spec)
zf=aligned[["rsid","POS","a1","a2","maf","BETA","SE","Z_aligned","N"]].copy()
zf.insert(1,"chromosome",int(chr_))
zf.rename(columns={"POS":"position","a1":"allele1","a2":"allele2",
                    "BETA":"beta","SE":"se","Z_aligned":"z","N":"n"},inplace=True)
zf["n"]=zf["n"].fillna(neff).round().astype(int)
zf["maf"]=zf["maf"].clip(0.01,0.49)
zf.to_csv(outdir/"finemap_input.z",sep=" ",index=False)

aligned[["rsid","original_rsid","POS"]].to_csv(outdir/"rsid_map.tsv",sep="\t",index=False)

# Print first 20 lines to TRANSCRIPT (makes allele errors visible without opening files)
print(f"\n=== Z-FILE PREVIEW: {locus_id} (first 20 lines) ===")
print(zf.head(20).to_string())
print(f"\n{locus_id}: {len(aligned)} variants in z-file")
```

**Step 4.2 — Lead-SNP |Z| sanity check (run after each z-file is written).**

A lead SNP with |Z| below the 25th percentile of its locus indicates a sign error or wrong-column bug.

```python
# scripts/check_lead_z.py <locus_id>
import pandas as pd, numpy as np, pathlib, sys, json

lid=sys.argv[1]
outdir=pathlib.Path(f"output/loci/{lid}")
loci=pd.read_csv("output/loci_with_counts.tsv",sep="\t")
row=loci[loci["locus_id"]==lid].iloc[0]
zf=pd.read_csv(outdir/"finemap_input.z",sep=" ")
z_abs=zf["z"].abs().sort_values(ascending=False)
p25=float(np.percentile(z_abs,25))
result={"locus_id":lid,"z_max":round(z_abs.max(),3),"z_p25":round(p25,3)}
lead=str(row.get("lead_rsid",""))
if lead and lead!="nan":
    m=zf[zf["rsid"]==lead]
    if not m.empty:
        lz=float(m["z"].abs().iloc[0])
        result.update({"lead_rsid":lead,"lead_z_abs":round(lz,3)})
        print(f"{lid}: lead {lead} |Z|={lz:.2f} vs locus p25={p25:.2f}")
        if lz<p25:
            msg=(f"SANITY_FAIL {lid}: lead rsID |Z|={lz:.2f} < p25={p25:.2f}. "
                 f"Likely sign error or wrong effect column. Fix allele alignment.")
            with open(outdir/"status.txt","a") as f: f.write(msg+"\n")
            raise SystemExit(msg)
    else:
        result["lead_note"]=f"{lead} not in z-file"
else:
    result["lead_note"]="no lead_rsid"
json.dump(result,open(outdir/"lead_z_check.json","w"),indent=2)
print(f"{lid}: z_max={result['z_max']:.2f} p25={p25:.2f}")
```

---

### Phase 5 — FINEMAP v1.4 (5-replicate median)

**Step 5.1 — Write per-locus config files.**

```python
# scripts/write_finemap_config.py <locus_id>
import pathlib, json, sys

lid=sys.argv[1]
outdir=pathlib.Path(f"output/loci/{lid}")
neff=int(float(next(l.split("=")[1] for l in open("output/provenance/neff_derivation.txt") if l.startswith("neff_used"))))

params={"locus_id":lid,"tool":"FINEMAP","version":"1.4",
        "n_causal_max":5,
        "prior_std":0.0224,"prior_std_justification":"Benner 2016 (PLoS Genet 12:e1005942 Table 1) binary-trait default; calibrated for case-control log(OR). Not 0.05.",
        "prob_csnp":1e-4,"n_iterations":100000,"n_convergence":1000,
        "n_replicates":5,"seed_strategy":"no --seed in v1.4; 5-replicate median used",
        "cred_coverage":0.95,"n_samples":neff}
json.dump(params,open(outdir/"finemap_params.json","w"),indent=2)

# Standalone rerun script for this locus
cmd=(f"#!/bin/bash\n# Standalone rerun: {lid}\n"
     f"FINEMAP=output/env/finemap_v1.4_x86_64/finemap_v1.4_x86_64\n"
     f"NEFF={neff}\n"
     f"mkdir -p output/loci/{lid}/rerun\n"
     f"printf 'z;ld;snp;config;cred;log;n_samples\\n"
     f"output/loci/{lid}/finemap_input.z;output/loci/{lid}/ld.ld;"
     f"output/loci/{lid}/rerun/finemap.snp;output/loci/{lid}/rerun/finemap.config;"
     f"output/loci/{lid}/rerun/finemap.cred;output/loci/{lid}/rerun/finemap.log;{neff}\\n' "
     f"> output/loci/{lid}/rerun/master\n"
     f"$FINEMAP --shotgun --in-files output/loci/{lid}/rerun/master "
     f"--n-causal-snps 5 --prior-std 0.0224 --prob-csnp 1e-4 "
     f"--n-iterations 100000 --n-convergence 1000 --cred 0.95\n")
p=outdir/"finemap_cmd.sh"; p.write_text(cmd); p.chmod(0o755)
print(f"Config written: {lid}")
```

**Step 5.2 — Run 5 replicates per locus.**

```bash
# scripts/run_finemap_locus.sh <locus_id>
set -euo pipefail
source output/config.sh
LID=$1; OUTDIR="output/loci/${LID}"
[ -f "${OUTDIR}/finemap_input.z" ] || { echo "SKIP_NO_ZFILE ${LID}"; echo "SKIP_FINEMAP_NO_ZFILE" >> "${OUTDIR}/status.txt"; exit 0; }
NEFF_VAL=$(grep neff_used output/provenance/neff_derivation.txt | cut -d= -f2 | xargs printf '%.0f')
echo "[$(date -u)] FINEMAP ${N_REPLICATES} reps: ${LID}"
for REP in $(seq 1 ${N_REPLICATES}); do
    RD="${OUTDIR}/rep${REP}"; mkdir -p "${RD}"
    printf "z;ld;snp;config;cred;log;n_samples\n${OUTDIR}/finemap_input.z;${OUTDIR}/ld.ld;${RD}/finemap.snp;${RD}/finemap.config;${RD}/finemap.cred;${RD}/finemap.log;${NEFF_VAL}\n" > "${RD}/master"
    $FINEMAP --shotgun --in-files "${RD}/master" \
        --n-causal-snps ${N_CAUSAL_MAX} --prior-std ${PRIOR_STD} \
        --prob-csnp 1e-4 --n-iterations 100000 --n-convergence 1000 \
        --cred ${CRED_COVERAGE} 2>&1 | tee "${RD}/finemap_run.log"
    echo "  rep${REP} done"
done
echo "FINEMAP_DONE: ${LID}" >> output/run_status.txt
```

**Step 5.3 — Aggregate replicates (median PIP) and QC.**

```python
# scripts/aggregate_replicates.py <locus_id>
import pandas as pd, numpy as np, pathlib, json, sys

lid=sys.argv[1]; outdir=pathlib.Path(f"output/loci/{lid}"); N_REPS=5
dfs=[]
for r in range(1,N_REPS+1):
    f=outdir/f"rep{r}/finemap.snp"
    if f.exists():
        d=pd.read_csv(f,sep=" "); d=d.rename(columns={"prob":f"pip_r{r}"}); dfs.append(d.set_index("rsid"))
if not dfs: raise SystemExit(f"No replicate outputs for {lid}")

pip_cols=[c for d in dfs for c in d.columns if c.startswith("pip_")]
comb=dfs[0]
for d in dfs[1:]: comb=comb.join(d[[c for c in d.columns if c.startswith("pip_")]],how="outer")
pc=[c for c in comb.columns if c.startswith("pip_")]
comb["pip_median"]=comb[pc].median(axis=1)
comb["pip_std"]=comb[pc].std(axis=1)
comb["pip_min"]=comb[pc].min(axis=1); comb["pip_max"]=comb[pc].max(axis=1)

corrs=[]
for i in range(len(pc)):
    for j in range(i+1,len(pc)):
        m=comb[pc[i]].notna()&comb[pc[j]].notna()
        if m.sum()>10: corrs.append(float(np.corrcoef(comb.loc[m,pc[i]],comb.loc[m,pc[j]])[0,1]))
mean_corr=float(np.mean(corrs)) if corrs else np.nan
max_diff=float((comb["pip_max"]-comb["pip_min"]).max())

n_sat=int((comb["pip_median"]>=0.99).sum()); frac_sat=n_sat/max(1,len(comb))
sat_warn=frac_sat>0.05
if sat_warn:
    ld_mat=pd.read_csv(outdir/"ld.ld",sep="\t",header=None).values
    cond=float(np.linalg.cond(ld_mat)); min_eig=float(np.linalg.eigvalsh(ld_mat).min())
    print(f"SANITY_WARN_PIP_SATURATED {lid}: {n_sat}/{len(comb)} variants PIP≥0.99 "
          f"({frac_sat*100:.1f}%). LD cond={cond:.2e}, min_eig={min_eig:.4f}. "
          f"Check ld.ld is Pearson r not r².")
    with open(outdir/"status.txt","a") as f:
        f.write(f"PIP_SATURATED: {n_sat}/{len(comb)} variants. cond={cond:.2e}\n")
    cond_out,mineig_out=round(cond,2),round(min_eig,6)
else: cond_out,mineig_out=None,None

qc={"locus_id":lid,"n_replicates":len(pc),"mean_pairwise_pip_corr":round(mean_corr,6) if not np.isnan(mean_corr) else None,
    "max_abs_pip_diff":round(max_diff,6),"n_pip_sat_ge99":n_sat,
    "pip_saturation_warn":sat_warn,"ld_condition_number":cond_out,"ld_min_eig":mineig_out}
json.dump(qc,open(outdir/"replicate_qc.json","w"),indent=2)
print(json.dumps(qc,indent=2))

out=comb.reset_index()
out["pip"]=out["pip_median"]
keep=["rsid","chromosome","position","allele1","allele2","maf","beta","se","z","pip","pip_std","pip_min","pip_max"]
keep=[c for c in keep if c in out.columns]
out[keep].sort_values("pip",ascending=False).to_csv(outdir/"finemap.snp",sep=" ",index=False)
print(f"Consensus PIPs written: {len(out)} variants, lead_pip={out['pip'].max():.4f}")
```

---

### Phase 6 — Credible Set Extraction

```python
# scripts/extract_credible_sets.py <locus_id>
import pandas as pd, numpy as np, pathlib, json, sys

lid=sys.argv[1]; outdir=pathlib.Path(f"output/loci/{lid}")
cs_outdir=pathlib.Path("output/credible_sets"); cs_outdir.mkdir(exist_ok=True)
COVERAGE=0.95; PURITY_THR=0.1

if not (outdir/"finemap.snp").exists(): print(f"No .snp for {lid}"); sys.exit(0)
snps=pd.read_csv(outdir/"finemap.snp",sep=" ").sort_values("pip",ascending=False).reset_index(drop=True)

snplist=pd.read_csv(outdir/"ld.snplist",header=None,names=["rsid"])
snplist["idx"]=range(len(snplist)); ridx=dict(zip(snplist["rsid"],snplist["idx"]))
ld_mat=pd.read_csv(outdir/"ld.ld",sep="\t",header=None).values
sat_flag="PIP_SATURATED" in ((outdir/"status.txt").read_text() if (outdir/"status.txt").exists() else "")

# Parse FINEMAP .cred file (rep1) for CS membership
cs_records={}
cred_f=outdir/"rep1/finemap.cred"
if cred_f.exists():
    cur=None
    for line in open(cred_f):
        parts=line.strip().split()
        if not parts or parts[0].startswith("prob"): continue
        try: float(parts[0]); float(parts[1]); int(parts[2]); cur=f"CS{len(cs_records)+1}"; cs_records[cur]=[]; continue
        except (ValueError,IndexError): pass
        if cur and len(parts)>=2:
            try: float(parts[0]); cs_records[cur].append(parts[1])
            except ValueError: pass
if not cs_records:  # Fallback: build one CS by PIP cumsum
    cumsum=0.0; vs=[]
    for _,r in snps.iterrows():
        vs.append(r["rsid"]); cumsum+=r["pip"]
        if cumsum>=COVERAGE: break
    cs_records={"CS1":vs}

def purity(rsids):
    idxs=[ridx[r] for r in rsids if r in ridx]
    if len(idxs)<2: return np.nan
    sub=ld_mat[np.ix_(idxs,idxs)]**2; np.fill_diagonal(sub,np.nan); return float(np.nanmin(sub))

pur={cid:purity(rs) for cid,rs in cs_records.items()}
cov_ach={cid:min(1.0,float(snps[snps["rsid"].isin(rs)]["pip"].sum())) for cid,rs in cs_records.items()}
r2cs={r:cid for cid,rs in cs_records.items() for r in rs}

rsid_map={}
if (outdir/"rsid_map.tsv").exists():
    rm=pd.read_csv(outdir/"rsid_map.tsv",sep="\t"); rsid_map=dict(zip(rm["rsid"],rm["original_rsid"]))

snps["cs_id"]=snps["rsid"].map(r2cs).fillna("none")
snps["cs_coverage"]=snps["cs_id"].map(cov_ach)
snps["cs_purity_min_r2"]=snps["cs_id"].map(pur)
def conf(row):
    if row["cs_id"]=="none": return "not_in_cs"
    if sat_flag: return "pip_saturated"
    if not np.isnan(row.get("cs_purity_min_r2",np.nan)) and row["cs_purity_min_r2"]<PURITY_THR: return "low_confidence"
    return "in_cs"
snps["cs_confidence"]=snps.apply(conf,axis=1)
snps["locus_id"]=lid
snps["variant_id"]=snps["rsid"].map(rsid_map).fillna(snps["rsid"])  # original sumstats rsID
snps["ld_rsid"]=snps["rsid"]  # 1000G rsID

out=snps[["variant_id","ld_rsid","chromosome","position","allele1","allele2","maf",
          "pip","pip_std","cs_id","cs_coverage","cs_purity_min_r2","cs_confidence","locus_id"]].copy()
out.rename(columns={"chromosome":"chr","position":"pos"},inplace=True)

out.to_csv(cs_outdir/f"{lid}_credible_sets.tsv",sep="\t",index=False)
out.to_csv(outdir/"credible_sets.tsv",sep="\t",index=False)

cs_sum=[{"locus_id":lid,"cs_id":cid,"n_variants":len(rs),
          "cs_coverage":cov_ach.get(cid),"cs_purity_min_r2":pur.get(cid),
          "lead_pip":float(snps[snps["rsid"].isin(rs)]["pip"].max()) if rs else None}
         for cid,rs in cs_records.items()]
pd.DataFrame(cs_sum).to_csv(outdir/"cs_summary.tsv",sep="\t",index=False)
print(f"{lid}: {len(cs_records)} CS, lead_pip={out['pip'].max():.3f}")
```

**Step 6.2 — Rerun check (run on one representative locus after all 10 complete).**

```python
# scripts/rerun_check.py <locus_id>
import pandas as pd, numpy as np, subprocess, pathlib, json, sys

lid=sys.argv[1]; outdir=pathlib.Path(f"output/loci/{lid}")
orig=pd.read_csv(outdir/"finemap.snp",sep=" ").set_index("rsid")["pip"]

subprocess.run(["bash","scripts/run_finemap_locus.sh",lid],check=True)
subprocess.run(["python","scripts/aggregate_replicates.py",lid],check=True)

new=pd.read_csv(outdir/"finemap.snp",sep=" ").set_index("rsid")["pip"]
com=orig.index.intersection(new.index)
corr=float(np.corrcoef(orig[com],new[com])[0,1])
max_diff=float((orig[com]-new[com]).abs().max())

if corr>=0.999 and max_diff<=0.01: verdict="PASS"
elif corr>=0.95: verdict="WARN"
else:
    verdict="FAIL"
    with open(outdir/"status.txt","a") as f: f.write(f"RERUN_FAIL: corr={corr:.6f}\n")

result={"locus_id":lid,"pip_pearson_corr":round(corr,6),"max_abs_pip_diff":round(max_diff,6),"verdict":verdict}
json.dump(result,open("output/rerun_check.json","w"),indent=2)
with open("output/rerun_check.txt","a") as f: f.write(f"{lid}\tcorr={corr:.6f}\tmax_diff={max_diff:.6f}\t{verdict}\n")
print(json.dumps(result,indent=2))
if verdict=="WARN":
    print("WARN: corr < 0.999. Consider increasing N_REPLICATES to 10.")
```

---

### Phase 7 — Parallelized Orchestration with Scale-Out Knob

```python
#!/usr/bin/env python3
# scripts/run_pipeline.py
"""
PGC3 SCZ FINEMAP v1.4 reproduction pipeline.

Usage:
  python scripts/run_pipeline.py --loci-subset tuning10 --n-workers 16
  python scripts/run_pipeline.py --loci-subset all287 --dry-run
  python scripts/run_pipeline.py --loci-subset all287 --n-workers 32
"""
import argparse, subprocess, concurrent.futures, pathlib, pandas as pd, sys
from datetime import datetime

def parse_args():
    p=argparse.ArgumentParser(description="PGC3 SCZ FINEMAP v1.4 pipeline (Trubetskoy 2022)")
    p.add_argument("--loci-subset",choices=["tuning10","all287"],required=True,
                   help="tuning10=10-locus subset; all287=full locus set (same code path)")
    p.add_argument("--sumstats",default="output/sumstats_qc.tsv")
    p.add_argument("--n-workers",type=int,default=16)
    p.add_argument("--dry-run",action="store_true",help="List loci and exit")
    return p.parse_args()

STEPS=[
    (["bash","scripts/compute_ld_locus.sh","{lid}","{chr}","{start}","{end}"],"ld"),
    (["python","scripts/harmonize_locus.py","{lid}","{chr}","{start}","{end}","{ss}"],"harmonize"),
    (["python","scripts/check_lead_z.py","{lid}"],"lead_z"),
    (["python","scripts/write_finemap_config.py","{lid}"],"config"),
    (["bash","scripts/run_finemap_locus.sh","{lid}"],"finemap"),
    (["python","scripts/aggregate_replicates.py","{lid}"],"aggregate"),
    (["python","scripts/extract_credible_sets.py","{lid}"],"extract"),
]

def run_locus(row,ss):
    lid=str(row["locus_id"]); t0=datetime.utcnow()
    subs={"lid":lid,"chr":str(int(row["chr"])),"start":str(int(row["start"])),
          "end":str(int(row["end"])),"ss":ss}
    for tmpl,step in STEPS:
        cmd=[t.format(**subs) for t in tmpl]
        r=subprocess.run(cmd,capture_output=True,text=True)
        if r.returncode!=0: return(lid,f"FAILED_{step.upper()}",r.stderr[-300:])
        sf=pathlib.Path(f"output/loci/{lid}/status.txt")
        if sf.exists():
            st=sf.read_text()
            if "SKIP" in st: return(lid,"SKIP",st.strip()[:150])
            if "ABORT" in st: return(lid,"ABORT_ALLELE_ERROR",st.strip()[:150])
    return(lid,"OK",f"{(datetime.utcnow()-t0).total_seconds():.0f}s")

def main():
    args=parse_args()
    lf="output/tuning_subset.tsv" if args.loci_subset=="tuning10" else "output/loci_with_counts.tsv"
    if not pathlib.Path(lf).exists(): sys.exit(f"ERROR: {lf} not found. Run Phases 1-2 first.")
    loci=pd.read_csv(lf,sep="\t")
    print(f"\nLoci ({args.loci_subset}): {len(loci)}")
    print(loci[["locus_id","chr","start","end","n_variants"]].to_string())
    if args.dry_run: print(f"\n[DRY RUN] {len(loci)} loci listed."); sys.exit(0)
    pathlib.Path("output/run_status.txt").unlink(missing_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as ex:
        futs={ex.submit(run_locus,row,args.sumstats):row["locus_id"] for _,row in loci.iterrows()}
        for fut in concurrent.futures.as_completed(futs):
            lid,status,info=fut.result()
            print(f"  {lid}: {status} ({info[:80]})")
            with open("output/run_status.txt","a") as f: f.write(f"{lid}\t{status}\t{info[:150]}\n")
    df=pd.read_csv("output/run_status.txt",sep="\t",header=None,names=["locus_id","status","info"])
    print("\n--- Run Summary ---"); print(df["status"].value_counts().to_string())

if __name__=="__main__": main()
```

**Run these two commands and show output in transcript (mandatory for scale-out gate):**

```bash
python scripts/run_pipeline.py --help
python scripts/run_pipeline.py --loci-subset all287 --dry-run
```

---

### Phase 8 — Sensitivity Analysis

Run on one locus per tertile from the tuning subset. Vary k ∈ {1, 3, 5, 10}. Print results to transcript.

```bash
# scripts/sensitivity_k.sh <locus_id> <k>
set -euo pipefail
source output/config.sh
LID=$1; K=$2
SD="output/loci/${LID}/sensitivity_k${K}"; mkdir -p "${SD}"
cp "output/loci/${LID}/finemap_input.z" "${SD}/"
cp "output/loci/${LID}/ld.ld" "${SD}/"
NEFF_VAL=$(grep neff_used output/provenance/neff_derivation.txt | cut -d= -f2 | xargs printf '%.0f')
printf "z;ld;snp;config;cred;log;n_samples\n${SD}/finemap_input.z;${SD}/ld.ld;${SD}/finemap.snp;${SD}/finemap.config;${SD}/finemap.cred;${SD}/finemap.log;${NEFF_VAL}\n" > "${SD}/master"
$FINEMAP --shotgun --in-files "${SD}/master" \
    --n-causal-snps ${K} --prior-std ${PRIOR_STD} --prob-csnp 1e-4 \
    --n-iterations 100000 --n-convergence 1000 --cred ${CRED_COVERAGE} \
    2>&1 | tee "${SD}/finemap_run.log"
```

```python
# scripts/run_sensitivity.py
import pandas as pd, numpy as np, subprocess, pathlib, json

tuning=pd.read_csv("output/tuning_subset.tsv",sep="\t")
sens_loci=tuning.groupby("tertile").apply(lambda g:g.iloc[0]).reset_index(drop=True)["locus_id"].tolist()
print(f"Sensitivity loci (one per tertile): {sens_loci}")
results=[]
for lid in sens_loci:
    base=pathlib.Path(f"output/loci/{lid}/finemap.snp")
    if not base.exists(): continue
    bp=pd.read_csv(base,sep=" ").set_index("rsid")["pip"]
    for k in [1,3,5,10]:
        subprocess.run(["bash","scripts/sensitivity_k.sh",lid,str(k)],check=False)
        af=pathlib.Path(f"output/loci/{lid}/sensitivity_k{k}/finemap.snp")
        if not af.exists(): continue
        ap=pd.read_csv(af,sep=" ").set_index("rsid")["pip"]
        com=bp.index.intersection(ap.index)
        corr=float(np.corrcoef(bp[com],ap[com])[0,1]) if len(com)>1 else np.nan
        results.append({"locus_id":lid,"base_k":5,"alt_k":k,"pip_corr":round(corr,4)})

df=pd.DataFrame(results)
df.to_csv("output/sensitivity/k_sensitivity.tsv",sep="\t",index=False)
print("\n=== K SENSITIVITY (pip corr vs base k=5) ===")
print(df.pivot(index="locus_id",columns="alt_k",values="pip_corr").to_string())
print("k=1 corr < 0.9 → multi-causal locus; k≥10 adds negligible CS change → k=5 is stable.")
```

---

### Phase 9 — Final Audit

```python
# scripts/final_audit.py
import pathlib, json, pandas as pd, sys

errors, warnings = [], []
loci=pd.read_csv("output/tuning_subset.tsv",sep="\t")

# Root-level files
for f in ["output/provenance.json","output/provenance/neff_derivation.txt",
          "output/provenance/determinism_note.txt","output/provenance/ldsc/ldsc_intercept.txt",
          "output/provenance/qc_metrics.txt","output/config.sh",
          "output/env/tool_versions.txt","output/env/environment.yml",
          "output/rerun_check.txt","output/run_status.txt",
          "output/tuning_subset.tsv","output/loci_with_counts.tsv",
          "output/sensitivity/k_sensitivity.tsv"]:
    if not pathlib.Path(f).exists(): errors.append(f"MISSING: {f}")

# Tool version strings
if pathlib.Path("output/env/tool_versions.txt").exists():
    tv=pathlib.Path("output/env/tool_versions.txt").read_text()
    if "1.4" not in tv: errors.append("tool_versions.txt missing '1.4'")
    if "1.9" not in tv.lower(): errors.append("tool_versions.txt missing 'v1.9' (plink)")

# Per-locus checks
for _,row in loci.iterrows():
    lid=str(row["locus_id"]); ld=pathlib.Path(f"output/loci/{lid}")
    st=(ld/"status.txt").read_text() if (ld/"status.txt").exists() else ""
    if any(x in st for x in ["SKIP","RERUN_FAIL","ABORT"]): print(f"  {lid}: skipped/failed"); continue
    for f in ["finemap_input.z","finemap.snp","finemap_params.json",
              "harmonization_summary.json","replicate_qc.json","lead_z_check.json",
              "finemap_cmd.sh","ld.ld","ld.snplist","rsid_map.tsv",
              "rep1/finemap.snp","rep1/finemap.cred","rep1/finemap.log"]:
        if not (ld/f).exists(): errors.append(f"MISSING: output/loci/{lid}/{f}")
    cs=pathlib.Path(f"output/credible_sets/{lid}_credible_sets.tsv")
    if not cs.exists(): errors.append(f"MISSING: {cs}")
    else:
        df=pd.read_csv(cs,sep="\t")
        for col in ["variant_id","ld_rsid","chr","pos","allele1","allele2","pip","cs_id","cs_purity_min_r2","cs_confidence"]:
            if col not in df.columns: errors.append(f"SCHEMA: {cs} missing '{col}'")
    if (ld/"finemap_params.json").exists():
        p=json.load(open(ld/"finemap_params.json"))
        if abs(p.get("prior_std",0)-0.0224)>1e-6: errors.append(f"prior_std wrong for {lid}: {p.get('prior_std')}")
    if (ld/"harmonization_summary.json").exists():
        h=json.load(open(ld/"harmonization_summary.json"))
        if h.get("n_direct",1)==0: errors.append(f"ALLELE_ERROR: {lid} n_direct=0")

# Rerun verdicts — only PASS/WARN/FAIL allowed
if pathlib.Path("output/rerun_check.txt").exists():
    for line in pathlib.Path("output/rerun_check.txt").read_text().splitlines():
        if not line: continue
        if "FAIL" in line: warnings.append(f"Rerun FAIL: {line}")
        parts=line.split("\t")
        if len(parts)>=3 and parts[-1].strip() not in ["PASS","WARN","FAIL"]:
            errors.append(f"Invalid rerun verdict: {line[:80]}")

for w in warnings: print(f"WARN: {w}")
if errors:
    print("\nAUDIT FAILED:")
    for e in errors: print(f"  {e}")
    sys.exit(1)
else: print("\nAUDIT PASSED")
```

---

## Expected Output Layout

```
output/
├── provenance.json                     # SHA256 for sumstats, loci TSV, FINEMAP binary
├── provenance/
│   ├── neff_derivation.txt             # neff_used=N, method, subset note
│   ├── determinism_note.txt            # NO_SEED_FLAG or SEED_FLAG_PRESENT
│   ├── qc_metrics.txt                  # lambda_gc=X, n_input, n_output
│   ├── ldsc/ldsc_intercept.txt         # intercept value or "LDSC UNAVAILABLE"
│   ├── 1kg_checksums.txt
│   └── finemap_sha256.txt
├── config.sh                           # NEFF, FINEMAP, PLINK, PRIOR_STD, N_CAUSAL_MAX
├── env/
│   ├── tool_versions.txt               # FINEMAP 1.4, PLINK v1.9, bcftools, python, OS
│   ├── environment.yml
│   └── pip_freeze.txt
├── sumstats_qc.tsv
├── loci_with_counts.tsv                # all non-MHC autosomal loci with n_variants
├── tuning_subset.tsv                   # 10 loci (stratified, seed=42)
├── run_status.txt                      # locus_id\tOK/SKIP/FAILED_X\tinfo
├── rerun_check.txt                     # locus_id\tcorr=X\tmax_diff=Y\tPASS|WARN|FAIL
├── credible_sets/
│   └── <locus_id>_credible_sets.tsv   # per-variant, all 10 tuning loci
├── sensitivity/
│   └── k_sensitivity.tsv
└── loci/<locus_id>/
    ├── finemap_input.z                 # space-sep; rsid chrom pos a1 a2 maf beta se z n
    ├── ld.ld, ld.snplist               # Pearson r matrix + SNP order
    ├── ld_plink.bim                    # allele reference
    ├── rsid_map.tsv                    # 1000G rsid → original sumstats rsid
    ├── finemap_params.json             # prior_std=0.0224, k=5, justification
    ├── finemap_cmd.sh                  # standalone rerun script
    ├── harmonization_summary.json      # n_direct, n_flipped, flip_rate_pct
    ├── lead_z_check.json
    ├── replicate_qc.json               # pip_corr, pip_saturation, ld_cond
    ├── rep1/ … rep5/                   # individual FINEMAP outputs
    ├── finemap.snp                     # consensus (median) PIPs
    ├── credible_sets.tsv
    ├── cs_summary.tsv                  # one row per CS with purity and coverage
    └── sensitivity_k{1,3,10}/
```

`credible_sets/<locus_id>_credible_sets.tsv` required columns:
```
variant_id  ld_rsid  chr  pos  allele1  allele2  maf  pip  pip_std
cs_id  cs_coverage  cs_purity_min_r2  cs_confidence  locus_id
```
`cs_confidence` ∈ `{in_cs, low_confidence, pip_saturated, not_in_cs}`

---

## Gotchas

1. **FINEMAP v1.4 has no `--seed` flag.** Verify immediately post-install. Use 5-replicate median strategy. Rerun verdicts are exactly **PASS / WARN / FAIL** — no other strings.

2. **`n_direct == 0` is a hard abort.** If all aligned variants are flipped with zero direct matches, the allele encoding is wrong. Fix the EA column or OR→BETA conversion before proceeding. Never skip this check.

3. **LD matrix is Pearson r, NOT r².** `plink --r square` returns r. FINEMAP v1.4 expects r. Purity computation squares it to r² internally. Do not pre-square `ld.ld`.

4. **PGC3 daner uses OR, not BETA.** Always compute `BETA = log(OR)`. Passing raw OR as BETA silently produces wrong PIPs.

5. **`prior_std = 0.0224` everywhere.** In every FINEMAP invocation, `finemap_params.json`, and `config.sh`. The generic default 0.05 is incorrect for case-control binary traits.

6. **Neff from the sumstats file, never hardcoded.** The public PGC3 release is a subset (~58K vs ~126K). Document the gap; PIPs will be lower and CS wider — that is correct behavior, not a bug.

7. **PIP saturation (>5% of variants at PIP ≥ 0.99) signals LD rank deficiency.** Compute condition number of `ld.ld`. Likely cause: matrix is r² not r, or wrong genomic window. Log and flag as `pip_saturated` in output.

8. **plink 1.9 only.** plink2 uses pgen format and different `--r` semantics. Mixing them produces silent format errors. Check `plink --version` says "v1.9" before every LD step.

9. **Print z-file preview and LD spot-check to STDOUT.** For at least the first locus, print the first 20 lines of `finemap_input.z` and the 4×4 top-left block of `ld.ld`. Allele errors are immediately visible without file inspection.

10. **`variant_id` must be the original sumstats rsID.** Store the 1000G rsID separately as `ld_rsid`. Positional IDs cannot be joined against published credible-set tables.

11. **Never read from `reference/`.** Not for debugging, validation, or any other purpose. All comparison to published results is judge-only.

12. **Do not substitute loci on failure.** Log SKIP/ABORT in `run_status.txt` and continue with N-1 loci.

---

## Safety

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

This pipeline processes public GWAS summary statistics only. The 1000G reference panel is used for LD estimation.

---

## Agent Boundary

The agent reads this SKILL.md, writes the scripts under `scripts/`, monitors gate checkpoints, and never reads `reference/`. Pipeline scripts invoke FINEMAP, PLINK, and Python. The agent verifies gate conditions by running the checkpoint commands shown, not by examining output files directly.

---

## Chaining Partners

- `skills/fine-mapping/`: SuSiE cross-validation on the same 10 loci (CS overlap as sensitivity)
- `skills/vcf-annotator/`: downstream annotation of top-PIP variants
- `skills/gwas-lookup/`: federated database lookup of credible-set variants

---

## Maintenance

Review when: PGC releases updated SCZ sumstats; FINEMAP v1.4.x adds `--seed` (remove replicate strategy); 1000G FTP paths change; christianbenner.com FINEMAP binary unavailable. Deprecation: PGC4 SCZ supersedes PGC3 with incompatible locus boundaries.