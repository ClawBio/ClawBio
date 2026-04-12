# Fine-mapping for tight credible sets

## Workflow

1. **Prepare output directory:**
   ```bash
   mkdir -p /tmp/finemapping_output
   ```

2. **Run fine-mapping on locus L1** with SuSiE parameters: L=3, coverage=0.95, min_purity=0.8, prior_variance=50:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L1.tsv \
     --ld data/locus_L1_ld.npy \
     --output /tmp/finemapping_output/L1 \
     --L 3 \
     --coverage 0.95 \
     --min-purity 0.8 \
     --prior-variance 50
   ```
   Extract `credible_set` (array of rsIDs) and `L` (integer) from `/tmp/finemapping_output/L1/result.json`.

3. **Run fine-mapping on locus L2** with same parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L2.tsv \
     --ld data/locus_L2_ld.npy \
     --output /tmp/finemapping_output/L2 \
     --L 3 \
     --coverage 0.95 \
     --min-purity 0.8 \
     --prior-variance 50
   ```
   Extract `credible_set` and `L` from `/tmp/finemapping_output/L2/result.json`.

4. **Run fine-mapping on locus L3** with same parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L3.tsv \
     --ld data/locus_L3_ld.npy \
     --output /tmp/finemapping_output/L3 \
     --L 3 \
     --coverage 0.95 \
     --min-purity 0.8 \
     --prior-variance 50
   ```
   Extract `credible_set` and `L` from `/tmp/finemapping_output/L3/result.json`.

5. **Assemble output JSON** combining all three loci:
   ```json
   {
     "L1": {
       "credible_set": [<rsid list from step 2>],
       "L": <L value from step 2>
     },
     "L2": {
       "credible_set": [<rsid list from step 3>],
       "L": <L value from step 3>
     },
     "L3": {
       "credible_set": [<rsid list from step 4>],
       "L": <L value from step 4>
     }
   }
   ```

## Output JSON