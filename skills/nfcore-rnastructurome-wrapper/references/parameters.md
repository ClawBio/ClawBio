# nf-core/rnastructurome parameter reference

Generated from `nextflow_schema.json` in [nf-core/rnastructurome](https://github.com/nf-core/rnastructurome) v1.0.0.
Regenerate after a pipeline version bump rather than hand-editing: see the skill's SKILL.md `## Maintenance` section.
This is the full audited surface; see the skill's own CLI Reference section for the subset used in everyday runs.

## Input/output options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--input` (required) | string | `—` | Path to comma-separated file containing information about the samples in the experiment. |
| `--outdir` (required) | string | `—` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. |
| `--email` | string | `—` | Email address for completion summary. |
| `--multiqc_title` | string | `—` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. |

## Reference genome options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | string | `—` | Path to transcript FASTA file. |
| `--gtf` | string | `—` | Path to transcript annotation GTF file. |
| `--genomes` | object | `—` | Optional custom reference map keyed by organism/reference name. |
| `--ensembl_base_url` | string | `https://ftp.ensembl.org/pub` | Ensembl FTP base URL used when auto-resolving transcript FASTA. |
| `--ensembl_release` | string | `current` | Ensembl release for transcript FASTA auto-resolution (`current`, `latest`, or release number like `114`). |
| `--ensembl_species_map` | object | `—` | Optional map from organism keys to Ensembl species names used for transcript FASTA auto-resolution. |
| `--ncbi_accessions_map` | object | `—` | Internal map from normalized viral organism keys to NCBI accession lists. |

## Sample metadata options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--sample_id` | string | `—` | Fallback sample ID for rows missing `sample_id`. |
| `--method` | string | `—` | Fallback probing method for rows missing `method`. |
| `--principle` | string | `—` | Fallback probing principle for rows missing `principle`. Allowed: RT-stop, MaP, rt-stop, map. |
| `--organism` | string | `—` | Fallback organism for rows missing `organism`, typically a Latin binomial such as `Homo sapiens`. |
| `--chemical` | string | `—` | Fallback chemical probing reagent for rows missing `chemical` (e.g. `1M7`). |
| `--RT_enzyme` | string | `—` | Fallback reverse transcriptase enzyme for rows missing `RT_enzyme` (e.g. `M-MLV`). |
| `--pH` | number | `—` | Fallback DMS reaction pH for rows missing `pH`. |
| `--umi_pattern` | string | `—` | Fallback UMI pattern for rows missing `umi_pattern`. Supplying a pattern enables umi_tools extract before cutadapt. |
| `--fuzzy_untreated_pairing` | boolean | `True` | When true (default), a treated group with no exact `sample_group`+`replicate` untreated match falls back to the untreated sample sharing the same sample_group base token (portion before the first `_`) at the same replicate — e.g. `MDA-MB-231_untreated_r1` can serve as control for `MDA-MB-231_MTX_treated_r1`. If a treated group still has no untreated after that, and exactly one untreated control exists anywhere else in the same reference, that single control is reused for it (with a warning) — e.g. one shared untreated backing several treated replicates. Set false to require exact matches; unmatched treated groups then proceed without an untreated control (scoring method 2 or 4). |

## Read trimming options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--cutadapt_adapter_5p` | string | `AGATCGGAAGAGC` | Global fallback 5' adapter sequence(s), comma-separated. If omitted and no per-sample adapter is provided, the pipeline defaults to AGATCGGAAGAGC. |
| `--cutadapt_adapter_3p` | string | `AGATCGGAAGAGC` | Global fallback 3' adapter sequence(s), comma-separated. If omitted and no per-sample adapter is provided, the pipeline defaults to AGATCGGAAGAGC. |
| `--cutadapt_quality_only` | boolean | `False` | If true, skip adapter trimming and apply quality/length trimming only. |
| `--cutadapt_5quality` | integer | `20` | Base quality threshold for 5' trimming (used for MaP; RT-stop is forced to 0). |
| `--cutadapt_3quality` | integer | `20` | Base quality threshold for 3' trimming. |
| `--cutadapt_len` | integer | `25` | Minimum read length retained after trimming. |
| `--cutadapt_min_align` | integer | `3` | Minimum adapter overlap required for trimming. |
| `--cutadapt_trim_n` | boolean | `True` | Trim terminal N bases. |

## Alignment options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--transcriptome` | boolean | `False` | Align against a transcript-level FASTA using Bowtie (RT-stop) / Bowtie2 (MaP). Default (false) downloads a genome FASTA and aligns with STAR. |
| `--count_genome` | boolean | `False` | Genome route (STAR) quantification method. Default (false) uses STAR --quantMode TranscriptomeSAM projected onto transcripts, tag-corrected with samtools calmd, counted with rf-count. Set true for the legacy rf-count-genome + rf-rctools extract path. |
| `--star_multimap_nmax` | integer | `10` | Maximum number of loci a read is allowed to map to in STAR (`--outFilterMultimapNmax`). Reads mapping to more loci than this are discarded. Increase when studying repetitive non-coding RNAs. |
| `--bowtie_k` | integer | `1` | Report up to this number of alignments per read (`-k`). |
| `--bowtie_all` | boolean/string | `True` | Report all valid alignments (`-a`). |
| `--bowtie_trim5` | integer | `0` | Trim this many bases from 5' read end (`--trim5`). |
| `--bowtie_trim3` | integer | `0` | Trim this many bases from 3' read end (`--trim3`). |
| `--bowtie_index` | string | `—` | Optional path to a prebuilt Bowtie index (currently not consumed by workflow logic). |
| `--bowtie_n` | integer | `2` | Bowtie v1 `-n` seed mismatch mode value. |
| `--bowtie_v` | integer | `—` | Bowtie v1 `-v` mismatch mode value. |
| `--bowtie_max` | integer | `1` | Bowtie v1 `-m` maximum reportable alignments per read before suppression. |
| `--bowtie_chunkmbs` | integer | `512` | Bowtie v1 `--chunkmbs` memory chunk size. |
| `--bowtie2_preset` | string | `--very-sensitive-local` | Bowtie2 alignment preset. Default is equivalent to `--local -N 0 -D 20 -R 3 -L 20 -i S,1,0.50`. |
| `--bowtie2_mp` | string | `6,2` | Bowtie2 `--mp` mismatch penalties as `max,min`. |
| `--bowtie2_dpad` | integer | `15` | Bowtie2 `--dpad` padding around DP table. |
| `--bowtie2_rdg` | string | `5,3` | Bowtie2 `--rdg` read gap penalties as `open,extend`. |
| `--bowtie2_rfg` | string | `5,3` | Bowtie2 `--rfg` reference gap penalties as `open,extend`. |
| `--bowtie2_softclip` | boolean | `—` | Enable Bowtie2 local alignment mode (`--local`). |
| `--bowtie2_ma` | integer | `2` | Bowtie2 `--ma` match bonus used with `--local`. |
| `--bowtie2_dovetail` | boolean | `True` | Allow dovetailing mates (`--dovetail`). |
| `--skip_markdup` | boolean | `True` | Skip SAMtools markdup deduplication (default: true). Position-based duplicate removal is not valid for chemical-probing libraries without UMIs, where reads sharing a 5' start are independent RT events rather than PCR duplicates. Leave enabled; UMI libraries are still deduplicated via UMI-tools. Set false only if you have a specific reason. |
| `--star_map_sjdb_overhang` | integer | `200` | STAR `--sjdbOverhang` for MaP alignment. Ideally `readLength - 1`; default 200 suits 201 bp reads. |

## RNAframework options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--rfcount_img` | boolean | `False` | Deprecated: rf-count plot generation (`-g`) is always enabled. |
| `--rfcount_strandedness` | string | `unstranded` | Library strandedness fallback for rf-count-genome when RSeQC inference is unavailable (e.g. viral/bacterial references without a usable BED annotation). Allowed: unstranded, first, second. |
| `--rfcount_genome_map_quality` | integer | `0` | Minimum MAPQ score for rf-count-genome (`--map-quality`). Useful for filtering multi-mappers; 0 keeps all alignments. |
| `--rfcount_genome_block_size` | integer | `100000` | Block size for per-chromosome memory allocation in rf-count-genome (`-bs`); increase for large genomes. |
| `--rfcount_genome_samtools_threads` | integer | `1` | SAMtools working threads per processor instance in rf-count-genome (`-wt`); total threads = `-p` × `-wt`. |
| `--rfcount_map_median_quality` | integer | `20` | MaP only: discard reads with median Phred+33 base quality below this value (`-eq`). |
| `--rfcount_map_no_deletions` | boolean | `False` | MaP only: ignore all deletion events during mutation counting (`-nd`). |
| `--rfcount_map_no_ambiguous` | boolean | `—` | MaP only: ignore ambiguously aligned deletion events (`-na`). Auto: off for M-MLV, on for Group II Intron (e.g. TGIRT); leave unset to use the auto default. |
| `--rfcount_map_no_insertions` | boolean | `—` | MaP only: ignore insertion events during mutation counting (`-ni`). Auto: off for M-MLV, on for Group II Intron (e.g. TGIRT); leave unset to use the auto default. |
| `--rfcount_map_right_deletion` | boolean | `—` | MaP only: mark only the right-most base of a deletion as mutated (`--right-deletion`). Auto: on for M-MLV, off for Group II Intron (e.g. TGIRT); leave unset to use the auto default. |
| `--rfcount_map_discard_consecutive` | integer | `—` | MaP only: discard mutations within N nt of each other (`-dc N`). Auto: 3 when `RT_enzyme` is a Group II Intron RT (e.g. TGIRT), off for M-MLV; leave unset to use the auto default. |
| `--rfcount_map_eval_surrounding` | boolean | `True` | MaP only: also evaluate quality of ±1 nt surrounding each mutation (`-es`). |
| `--rfcount_trim_5prime` | integer | `0` | Trim this many bases from the 5' read end for rf-count (`-t5`). |
| `--rfcount_mask_file` | string | `—` | Path to rf-count mask file (`--mask-file`). |
| `--rfcount_primary_only` | boolean | `False` | Use only primary alignments in rf-count (`--primary-only`). |
| `--rfcount_discard_clipped_reads` | boolean | `False` | Discard reads with excessive 5' or 3' end clipping in rf-count / rf-count-genome (`--discard-clipped-reads`). |
| `--rfcount_paired_only` | boolean | `False` | For paired-end samples, use only reads with both mates mapped (`--paired-only`). |
| `--rfcount_properly_paired` | boolean | `False` | For paired-end samples, use only properly paired reads (`--properly-paired`). |
| `--rfcount_map_sort_by_read_name` | boolean | `True` | In mutation mode (`-m`), pre-sort paired-end reads by read name (`--sort-by-read-name`). No effect on single-end data. |
| `--rfcount_map_discard_shorter` | integer | `1` | In mutation mode (`-m`), discard reads shorter than this value (`--discard-shorter`). |
| `--rfcount_map_min_quality` | integer | `20` | In mutation mode (`-m`), minimum Phred+33 base quality for mutations (`--min-quality`). |
| `--rfcount_map_collapse_consecutive` | boolean | `—` | In mutation mode (`-m`), collapse consecutive mutations/indels (`--collapse-consecutive`). Auto: on unless `RT_enzyme` is a Group II Intron RT (e.g. TGIRT); leave unset to use the auto default. |
| `--rfcount_map_max_collapse_distance` | integer | `2` | In mutation mode with collapsing enabled, max distance to collapse (`--max-collapse-distance`). |
| `--rfnorm_remap_reactivities` | boolean | `False` | Remap normalized reactivities to the 0-1 range according to Zarringhalam et al. (`--remap-reactivities`). |
| `--rfnorm_reactive_bases` | string | `—` | Reactive bases used for normalization window selection (`--reactive-bases`), e.g. `AC` for DMS. |
| `--rfnorm_norm_window` | integer | `—` | Normalization window size (`--norm-window`). Auto: 50 for DMS or RT-stop, unset otherwise; leave unset to use the auto default. |
| `--rfnorm_window_offset` | integer | `—` | Normalization window offset (`--window-offset`). |
| `--rfnorm_dynamic_window` | integer/boolean | `—` | Enable dynamic normalization window (`--dynamic-window`). Set to a positive integer to enable (value is passed to `--norm-window`); set to `false` or `0` to explicitly disable even for DMS (which enables it by default). Leave unset to use the per-method default. |
| `--rfnorm_norm_independent` | boolean | `False` | Normalize each reactive base independently (`--norm-independent`). |
| `--rfnorm_use_normfactor` | boolean | `—` | Cross-experiment normalization via rf-normfactor (one factor set per reference, fed to rf-norm via `-nf`), putting reactivities on a common scale. Leave unset for auto (enabled only for a reference with more than one treated sample to cross-normalise); set `true` to force on, `false` to force off (per-sample box-plot). |
| `--rfnorm_normfactor_min_coverage` | integer | `—` | rf-normfactor minimum coverage (`-mc`): bases below this coverage are excluded from the normalization factor calculation. Auto: 1000 for MaP, 50 for RT-stop; leave unset to use the auto default. |
| `--rfnorm_score_method` | integer | `—` | Override rf-norm scoring method (`-sm`). 1=Ding, 2=Rouskin, 3=Siegfried, 4=Zubradt. Allowed: 1, 2, 3, 4. |
| `--rfnorm_norm_method` | integer | `—` | Override rf-norm normalization method (`-nm`). 2=90% Winsorizing, 3=Box-plot, 4=Mitchell (MaP only). Allowed: 2, 3, 4. |
| `--rfnorm_raw` | boolean | `False` | Score raw reactivities without applying normalization (`--raw`). |
| `--rfnorm_pseudocount` | number | `—` | Pseudocount used by Ding scoring (`--pseudocount`). |
| `--rfnorm_ignore_lower_than_untreated` | boolean | `False` | Set reactivities lower than untreated to zero for Ding or Siegfried scoring (`--ignore-lower-than-untreated`). |
| `--rfnorm_mean_coverage` | number | `0` | Discard transcripts with mean coverage below this threshold (`--mean-coverage`). |
| `--rfnorm_median_coverage` | number | `0` | Discard transcripts with median coverage below this threshold (`--median-coverage`). |
| `--rfnorm_nan` |  | `—` | Positions with read coverage below this threshold are reported as NaN (`--nan`). Auto: 1000 for MaP, 50 for RT-stop. Set to 0 to disable NaN masking. |
| `--rfnorm_prefilter_min_coverage` | integer | `1` | Genome route only: after rf-rctools extract, keep transcripts with at least one position at this coverage or higher before rf-norm. Set to 0 to keep the full annotation RC. |
| `--rfnorm_img` | boolean | `False` | Deprecated: rf-norm plot generation (`--img`) is always enabled. |
| `--rffold_img` | boolean | `False` | Generate ViennaRNA RNAplot structure diagrams via rf-fold (`-g`). Disabled by default; slow on large transcriptomes. Use `--r2dt` for template-based diagrams instead. |
| `--jackknife_reference` | string | `—` | Path to a reference `.db` structure file for rf-jackknife normalisation assessment. When provided, rf-jackknife runs between rf-norm and rf-fold and the output CSV reports optimal slope/intercept values. When omitted, rf-fold runs directly using whatever slope/intercept params are set. |
| `--rfjackknife_slope` | string | `0,5` | Comma-separated slope range to search in rf-jackknife (`-sl`), e.g. `0,5`. |
| `--rfjackknife_intercept` | string | `-3,0` | Comma-separated intercept range to search in rf-jackknife (`-in`), e.g. `-3,0`. |
| `--rfjackknife_slope_step` | number | `0.2` | Step size for slope grid search in rf-jackknife (`-ss`). |
| `--rfjackknife_intercept_step` | number | `0.2` | Step size for intercept grid search in rf-jackknife (`-is`). |
| `--rfjackknife_mfmi` | boolean | `True` | Use modified FMI (Lan et al. 2022) instead of standard FMI in rf-jackknife (`-m`). |
| `--rfjackknife_relaxed` | boolean | `True` | Use relaxed FMI criteria (Deigan et al. 2009) in rf-jackknife (`-x`). |
| `--rfjackknife_keep_lonelypairs` | boolean | `True` | Retain lonely base-pairs (1 bp helices) in the reference structure before comparison (`-kl`). |
| `--rfjackknife_keep_pseudoknots` | boolean | `—` | Retain pseudoknotted base-pairs in the reference structure during FMI comparison (`-kp`). Disabled by default because rf-fold cannot predict pseudoknots, so including them systematically lowers FMI. |
| `--rfjackknife_only_common` | boolean | `—` | Only use transcripts present across all replicates in rf-jackknife (`-oc`). |
| `--rfjackknife_img` | boolean | `—` | Generate R heatmap of grid-search results in rf-jackknife (`-g`). |
| `--rfjackknife_rf_fold_params` | string | `-md 600` | Additional rf-fold parameters passed inside rf-jackknife (`-rp`), e.g. `-md 500`. |
| `--rfjackknife_pool_all` | boolean | `True` | Pool XMLs from all sample groups into a single rf-jackknife run (output: `jackknife/all_groups/`) instead of running one jackknife per group. Default: true. |
| `--stop_after_jackknife` | boolean | `False` | Stop the pipeline after rf-jackknife completes, skipping rf-fold, structure visualisation, browser track generation, and RDAT output. Useful for calibration runs where you only want jackknife statistics. |
| `--rfrctools_gtf_feature` | string | `exon` | GTF feature type to extract with rf-rctools (`-f`); override for non-standard annotations. |
| `--rfrctools_gtf_attribute` | string | `transcript_id` | GTF attribute used as the output RC entry ID by rf-rctools (`-b`). |
| `--rfeval_reference` | string | `—` | Path to a `.db` file of known RNA secondary structures. Enables rf-eval when provided. |
| `--rfeval_windows` | string | `—` | Optional manifest of sub-region references, one per line: `ref_seq_id start end structure_id [strand]` (1-based inclusive, in the XML's own coordinate space). When set, rf-norm reactivities are sliced to each window before rf-eval, so a sub-region structure (e.g. an Rfam element on a whole chromosome) is scored against a matching windowed XML instead of the diluted full transcript. `structure_id` must match an entry in `--rfeval_reference`. |
| `--rfeval_reactivity_cutoff` | number | `0.7` | Cutoff for classifying a base as highly-reactive in rf-eval unpaired-coefficient calculation (`-c`). |
| `--rfeval_ignore_terminal` | boolean | `True` | Exclude terminal base-pairs from rf-eval calculations (`-it`). |
| `--rfeval_terminal_as_unpaired` | boolean | `False` | Treat terminal base-pairs as unpaired in rf-eval (`-tu`). |
| `--rfeval_keep_pseudoknots` | boolean | `True` | Retain pseudoknotted base-pairs in rf-eval (`-kp`). |
| `--rfeval_keep_lonelypairs` | boolean | `True` | Retain lonely/isolated base-pairs in rf-eval (`-kl`). |
| `--rfeval_img` | boolean | `False` | Generate R metric plots in rf-eval (`-g`). |
| `--structextract` | boolean | `False` | Run rf-structextract after rf-fold to extract high-confidence, low-reactivity / low-Shannon structural motifs from the folded structures. |
| `--structextract_win_size` | integer | `50` | rf-structextract window size in nt for median reactivity/Shannon calculations (`-w`). |
| `--structextract_min_transcript_len` | integer | `500` | rf-structextract skips low-reactivity/low-Shannon evaluation for transcripts below this length (`-ml`). |
| `--structextract_min_value_frac` | number | `0.4` | rf-structextract windows with less than this fraction of bases covered are set to NaN (`-mv`). |
| `--structextract_min_below_median` | number | `0.7` | rf-structextract minimum fraction of bases whose Shannon and reactivity are below the transcript median (`-mb`). |
| `--structextract_min_paired_frac` | number | `0.45` | rf-structextract discards elements with less than this fraction of paired bases (`-mp`). |
| `--structextract_min_motif_len` | integer | `50` | rf-structextract discards structure elements shorter than this (`-mm`). |
| `--structextract_max_motif_len` | integer | `—` | rf-structextract discards structure elements longer than this (`-xm`). Unset for no limit. |
| `--structextract_max_loop_size` | integer | `—` | rf-structextract discards elements with a loop larger than this (`-xl`). Unset for no limit. |
| `--structextract_ignore_react` | boolean | `False` | rf-structextract skips low-reactivity evaluation (`-ir`). When false (default), reactivity is evaluated so only regions with probing support are extracted. |
| `--structextract_ignore_shannon` | boolean | `False` | rf-structextract skips low-Shannon evaluation (`-is`). When `false`, the Shannon-entropy test is applied so only high-confidence regions are extracted. |
| `--structextract_multiway_only` | boolean | `False` | rf-structextract only reports elements encompassing multiway junctions (`-mo`). |
| `--structextract_one_per_file` | boolean | `False` | rf-structextract reports each extracted element in a separate file (`-opf`). |
| `--structextract_eval_energy` | boolean | `False` | rf-structextract only reports elements with a free energy significantly lower than expected by chance (`-ee`). |
| `--structextract_pvalue` | number | `0.05` | rf-structextract p-value threshold for energy significance (`-v`); requires `--structextract_eval_energy`. |
| `--structextract_n_shufflings` | integer | `100` | rf-structextract number of sequence shufflings for energy evaluation (`-ns`). |
| `--structextract_dinucl_shuffle` | boolean | `False` | rf-structextract preserves dinucleotide frequencies when shuffling for energy evaluation (`-ds`). |
| `--structextract_plot` | boolean | `True` | Render an SVG diagram per extracted rf-structextract motif via ViennaRNA RNAplot. |
| `--correlate_replicates` | boolean | `True` | Run rf-correlate to compute pairwise reactivity-profile correlations between replicates of each sample group (a replicate-reproducibility QC surfaced in MultiQC). Only runs for sample groups with more than one replicate; a no-op otherwise. |
| `--rfcorrelate_cap_react` | number | `1.5` | Cap reactivities to this value before Pearson correlation in rf-correlate (`--cap-react`). Not applied to the Spearman run, which is rank-based. |
| `--correlate_min_values` | number | `—` | rf-correlate minimum number of values to calculate a correlation (`-m`); a value between 0 and 1 is interpreted as a fraction of transcript length. Unset uses the tool default (off). |
| `--correlate_ignore_sequence` | boolean | `False` | rf-correlate ignores sequence differences (e.g. SNVs) between compared transcripts (`-I`). |
| `--correlate_img` | boolean | `False` | Generate the rf-correlate correlation heatmap PDF (`-g`; requires R). |
| `--rffold_ct` | boolean | `False` | Write CT format structures with rf-fold (`-ct`). Disabled by default to keep dot-bracket output. |
| `--rffold_window` | integer | `1000` | Enable windowed MFE folding in rf-fold with this window size (`-w -fw N`). Set null/unset to fold the whole transcript instead. |
| `--rffold_partition_window` | integer | `1000` | Partition-function window size in rf-fold (`-pw N`). |
| `--rffold_unconstrained` | boolean | `False` | Fold without reactivity constraints in rf-fold (`-u`). |
| `--rffold_vienna_no_lonely_pairs` | boolean | `False` | Pass ViennaRNA no-lonely-pairs mode to rf-fold (`-vnlp`). |
| `--rffold_vienna_constrained` | boolean | `False` | Use ViennaRNA hard constraints with rf-fold (`-vc`). |
| `--rffold_vienna_max_bp_span` | integer | `600` | Maximal base-pair span for ViennaRNA in rf-fold (`-vmd`). |
| `--rffold_vienna_bp_span` | integer | `—` | Minimal base-pair span for ViennaRNA in rf-fold (`-vms`). |
| `--rffold_only_common` | integer | `—` | Only fold transcripts covered in at least this number of XML experiments (`-oc`). |
| `--rffold_fold_constraint_file` | string | `—` | Constraint file for allowed base-pairing positions (`-fc`). |
| `--rffold_unpaired_constraint_file` | string | `—` | Constraint file for required unpaired positions (`-uc`). |
| `--rffold_dotplot` | boolean | `True` | Generate dot plots from rf-fold (`-d`). |
| `--rffold_shannon_entropy` | boolean | `True` | Compute and report Shannon entropy in rf-fold (`-sh`). |
| `--rffold_slope` | number | `—` | Slope for reactivity-to-folding-constraint conversion in rf-fold (`-sl`). Auto by `chemical`: DMS 4.6, NAI 2.2, 2A3 1, other/unset 1.8. Overridden by jackknife calibration when available. |
| `--rffold_intercept` | number | `—` | Intercept for reactivity-to-folding-constraint conversion in rf-fold (`-in`). Auto by `chemical`: DMS -2, NAI -0.8, 2A3 -0.4, other/unset -0.6. Overridden by jackknife calibration when available. |
| `--rffold_vienna_rnaplot` | string | `RNAplot` | Path to ViennaRNA `RNAplot` binary for SVG structure plots in rf-fold (`-vrp`). |
| `--r2dt` | boolean | `True` | Generate template-based 2D RNA structure diagrams using R2DT, coloured by normalised SHAPE/chemical-probing reactivity. |
| `--r2dt_templatable_biotypes` | string | `rRNA,Mt_rRNA,tRNA,Mt_tRNA,snoRNA,scaRNA,snRNA,SRP_RNA,RNase_P_RNA,RNase_MRP_RNA,tmRNA` | Comma-separated GTF biotypes sent to R2DT. R2DT only has templates for structured ncRNA classes; transcripts of other biotypes (mRNA, lncRNA, ...) are excluded so R2DT never force-fits them to the wrong template, and are drawn by ViennaRNA instead. Matched case-insensitively against transcript/gene biotype. |
| `--rfwiggle_report_zeroes` | boolean | `True` | Report positions with zero reactivity in the WIG/BigWig output (`-z`). Enabled by default so genome browsers receive complete tracks without gaps. |
| `--rfwiggle_keep_bases` | string | `—` | Restrict WIG output to specific bases using an IUPAC code or combination (e.g. `AC` for DMS probing) (`-kb`). Default: all bases (`N`). |
| `--rnaframework_r_path` | string | `/usr/bin/R` | Path to R executable used by rf-count, rf-norm, and rf-fold plotting. |

## Institutional config options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--custom_config_version` | string | `master` | Git commit id for Institutional configs. |
| `--custom_config_base` | string | `https://raw.githubusercontent.com/nf-core/configs/master` | Base directory for Institutional configs. |
| `--config_profile_name` | string | `—` | Institutional config name. |
| `--config_profile_description` | string | `—` | Institutional config description. |
| `--config_profile_contact` | string | `—` | Institutional config contact information. |
| `--config_profile_url` | string | `—` | Institutional config URL link. |

## Executor options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--slurm_queue` | string | `—` | SLURM partition / queue name passed via `--partition`. |
| `--slurm_account` | string | `—` | SLURM account / allocation passed via `--account`. |
| `--slurm_clusterOptions` | string | `—` | Extra sbatch flags appended verbatim to `process.clusterOptions`. |

## Generic options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--version` | boolean | `—` | Display version and exit. |
| `--publish_dir_mode` | string | `copy` | Method used to save pipeline results to output directory. Allowed: symlink, rellink, link, copy, copyNoFollow, move. |
| `--email_on_fail` | string | `—` | Email address for completion summary, only when pipeline fails. |
| `--plaintext_email` | boolean | `—` | Send plain-text email instead of HTML. |
| `--max_multiqc_email_size` | string | `25.MB` | File size limit when attaching MultiQC reports to summary emails. |
| `--monochrome_logs` | boolean | `—` | Do not use coloured log outputs. |
| `--hook_url` | string | `—` | Incoming hook URL for messaging service |
| `--multiqc_config` | string | `—` | Custom config file to supply to MultiQC. |
| `--multiqc_logo` | string | `—` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file |
| `--multiqc_methods_description` | string | `—` | Custom MultiQC yaml file containing HTML including a methods description. |
| `--validate_params` | boolean | `True` | Boolean whether to validate parameters against the schema at runtime |
| `--pipelines_testdata_base_path` | string | `https://raw.githubusercontent.com/nf-core/test-datasets/` | Base URL or local path to location of pipeline test dataset files |
| `--trace_report_suffix` | string | `—` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. |
| `--help` | boolean/string | `—` | Display the help message. |
| `--help_full` | boolean | `—` | Display the full detailed help message. |
| `--show_hidden` | boolean | `—` | Display hidden parameters in the help message (only works when --help or --help_full are provided). |

