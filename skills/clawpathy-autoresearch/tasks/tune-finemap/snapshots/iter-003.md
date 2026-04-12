# Fine-mapping for tight credible sets

## Workflow

1. **Prepare output directory:**
   ```bash
   mkdir -p /tmp/finemapping_output
   ```

2. **Run fine-mapping on locus L1** with conservative SuSiE parameters (L=2, coverage=0.80, min_purity=0.95, prior_variance=25):
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L1.tsv \
     --ld data/locus_L1_ld.npy \
     --output /tmp/finemapping_output/L1 \
     --L 2 \
     --coverage 0.80 \
     --min-purity 0.95 \
     --prior-variance 25
   ```

3. **Run fine-mapping on locus L2** with same parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L2.tsv \
     --ld data/locus_L2_ld.npy \
     --output /tmp/finemapping_output/L2 \
     --L 2 \
     --coverage 0.80 \
     --min-purity 0.95 \
     --prior-variance 25
   ```

4. **Run fine-mapping on locus L3** with same parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L3.tsv \
     --ld data/locus_L3_ld.npy \
     --output /tmp/finemapping_output/L3 \
     --L 2 \
     --coverage 0.80 \
     --min-purity 0.95 \
     --prior-variance 25
   ```

5. **Extract credible set data** from each result JSON:
   - Read `/tmp/finemapping_output/L1/result.json` and extract fields `credible_set` and `L`
   - Read `/tmp/finemapping_output/L2/result.json` and extract fields `credible_set` and `L`
   - Read `/tmp/finemapping_output/L3/result.json` and extract fields `credible_set` and `L`

6. **Return JSON output** combining all three loci:
   ```json
   {
     "L1": {
       "credible_set": [rsid list from L1 result.json],
       "L": integer from L1 result.json
     },
     "L2": {
       "credible_set": [rsid list from L2 result.json],
       "L": integer from L2 result.json
     },
     "L3": {
       "credible_set": [rsid list from L3 result.json],
       "L": integer from L3 result.json
     }
   }
   ```

## Output JSON Shape