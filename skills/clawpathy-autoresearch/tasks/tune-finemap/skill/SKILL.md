# Fine-mapping for tight credible sets

## Workflow

1. **Run fine-mapping on locus L1** with single-signal SuSiE parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L1.tsv \
     --ld data/locus_L1_ld.npy \
     --output /tmp/finemapping_L1 \
     --L 1 \
     --coverage 0.50 \
     --min-purity 0.99 \
     --prior-variance 5
   ```

2. **Run fine-mapping on locus L2** with single-signal SuSiE parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L2.tsv \
     --ld data/locus_L2_ld.npy \
     --output /tmp/finemapping_L2 \
     --L 1 \
     --coverage 0.50 \
     --min-purity 0.99 \
     --prior-variance 5
   ```

3. **Run fine-mapping on locus L3** with single-signal SuSiE parameters:
   ```bash
   python /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio/skills/fine-mapping/fine_mapping.py \
     --sumstats data/locus_L3.tsv \
     --ld data/locus_L3_ld.npy \
     --output /tmp/finemapping_L3 \
     --L 1 \
     --coverage 0.50 \
     --min-purity 0.99 \
     --prior-variance 5
   ```

4. **Extract credible sets from result files:**
   - Read `/tmp/finemapping_L1/result.json` and extract the `credible_set` array and the `L` integer.
   - Read `/tmp/finemapping_L2/result.json` and extract the `credible_set` array and the `L` integer.
   - Read `/tmp/finemapping_L3/result.json` and extract the `credible_set` array and the `L` integer.

5. **Assemble and return output JSON with this exact structure:**
   ```json
   {
     "L1": {
       "credible_set": [<array of rsID strings from L1 result.json>],
       "L": <integer L value from L1 result.json>
     },
     "L2": {
       "credible_set": [<array of rsID strings from L2 result.json>],
       "L": <integer L value from L2 result.json>
     },
     "L3": {
       "credible_set": [<array of rsID strings from L3 result.json>],
       "L": <integer L value from L3 result.json>
     }
   }
   ```

## Output JSON Shape