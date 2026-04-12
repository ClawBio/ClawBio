# Fine-mapping for tight credible sets — parameter-optimized

## Workflow

1. **Prepare output directory:**
   ```bash
   mkdir -p /tmp/finemapping_output
   ```

2. **Run fine-mapping on locus L1** with balanced SuSiE parameters (L=2, coverage=0.70, min_purity=0.97, prior_variance=25):
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L1.tsv \
     --ld data/locus_L1_ld.npy \
     --output /tmp/finemapping_output/L1 \
     --L 2 \
     --coverage 0.70 \
     --min-purity 0.97 \
     --prior-variance 25
   ```

3. **Run fine-mapping on locus L2** with aggressive tightening (L=2, coverage=0.50, min_purity=0.99, prior_variance=25):
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L2.tsv \
     --ld data/locus_L2_ld.npy \
     --output /tmp/finemapping_output/L2 \
     --L 2 \
     --coverage 0.50 \
     --min-purity 0.99 \
     --prior-variance 25
   ```

4. **Run fine-mapping on locus L3** with balanced SuSiE parameters (L=2, coverage=0.70, min_purity=0.97, prior_variance=25):
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L3.tsv \
     --ld data/locus_L3_ld.npy \
     --output /tmp/finemapping_output/L3 \
     --L 2 \
     --coverage 0.70 \
     --min-purity 0.97 \
     --prior-variance 25
   ```

5. **Extract credible set rsID lists and L from each result.json:**
   - Read `/tmp/finemapping_output/L1/result.json`, copy field `credible_set` (list of rsIDs) and field `L` (integer)
   - Read `/tmp/finemapping_output/L2/result.json`, copy field `credible_set` (list of rsIDs) and field `L` (integer)
   - Read `/tmp/finemapping_output/L3/result.json`, copy field `credible_set` (list of rsIDs) and field `L` (integer)

6. **Assemble and return output JSON:**
   ```json
   {
     "L1": {
       "credible_set": [<rsID_list from L1 result.json>],
       "L": <integer from L1 result.json>
     },
     "L2": {
       "credible_set": [<rsID_list from L2 result.json>],
       "L": <integer from L2 result.json>
     },
     "L3": {
       "credible_set": [<rsID_list from L3 result.json>],
       "L": <integer from L3 result.json>
     }
   }
   ```

## Output JSON Shape

Executor must return exactly this structure: