include { readStructure } from 'plugin/nf-fgbio'

workflow SAMPLESHEET_TO_CHANNEL {
    take:
    ch_from_samplesheet           // samplesheet
    aligner                       // String: aligner
    ascat_alleles                 // Path: ascat alleles
    ascat_loci                    // Path: ascat loci
    ascat_loci_gc                 // Path: ascat loci gc
    ascat_loci_rt                 // Path: ascat loci rt
    bcftools_annotations          // Path: bcftools annotations
    bcftools_annotations_tbi      // Path: bcftools annotations tbi
    bcftools_columns              // Path: bcftools columns
    bcftools_header_lines         // Path: bcftools header lines
    build_only_index              // Boolean: build only index
    dbsnp                         // Path: dbsnp
    fasta                         // Path: fasta
    germline_resource             // Path: germline resource
    intervals                     // Path: intervals
    joint_germline                // Boolean: joint_germline
    joint_mutect2                 // Boolean: joint_mutect2
    known_indels                  // Path: known indels
    known_snps                    // Path: known snps
    no_intervals                  // Boolean: no intervals
    pon                           // Path: pon
    sentieon_dnascope_emit_mode   // String: sentieon dnascope emit mode
    sentieon_haplotyper_emit_mode // String: sentieon haplotyper emit mode
    seq_center                    // String: seq center
    seq_platform                  // String: seq platform
    skip_tools                    // Array: skip tools
    snpeff_cache                  // Path: snpeff cache
    snpeff_db                     // String: snpeff db
    step                          // String: step
    tools                         // Array: tools
    umi_length                    // Integer: umi length for fastp extraction
    umi_location                  // String: umi location for fastp extraction
    umi_in_read_header            // Boolean: umi in read header
    umi_read_structure            // String: umi read structure for fgbio consensus
    wes                           // wes

    main:
    ch_from_samplesheet.dump(tag: "ch_from_samplesheet")

    ch_from_samplesheet
        .map { meta, _fastq_1, _fastq_2, _spring_1, _spring_2, _table, _cram, _crai, _bam, _bai, _contamination, _vcf, _variantcaller ->
            // Get only the patient, sample and status fields from the meta map
            [meta.patient, meta.subMap('sample', 'status')]
        }
        .unique()
        .groupTuple()
        .map { patient, samples ->
            // Count samples with status 0 and status 1
            def status0_count = samples.count { it.status == 0 }
            def status1_count = samples.count { it.status == 1 }

            // Check the condition and exit with an error if met
            if (status1_count == 1 && status0_count > 1) {
                System.err.println("Patient [${patient}] has more than one sample [${status0_count}] with normal status [0] and one sample with tumor status [1].")
                error("Execution halted due to sample status inconsistency.")
            }
        }

    ch_from_samplesheet
        .map { meta, _fastq_1, _fastq_2, _spring_1, _spring_2, _table, _cram, _crai, _bam, _bai, _contamination, _vcf, _variantcaller ->
            // Create a unique key for patient-sample-status-lane combination
            def combination_key = "${meta.patient}-${meta.sample}-${meta.status}-${meta.lane}"
            [combination_key, [meta.patient, meta.sample, meta.status, meta.lane]]
        }
        .groupTuple()
        .map { combination_key, combination_list ->
            if (combination_list.size() > 1) {
                def patient = combination_list[0][0]
                def sample = combination_list[0][1]
                def status = combination_list[0][2]
                def lane = combination_list[0][3]
                System.err.println("Duplicate patient-sample-status-lane combination found: Patient '${patient}', Sample '${sample}', Status '${status}', Lane '${lane}' appears ${combination_list.size()} times. Please ensure each combination is unique.")
                error("Execution halted due to duplicate patient-sample-status-lane combination.")
            }
        }

    ch_from_samplesheet
        .map { meta, _fastq_1, _fastq_2, _spring_1, _spring_2, _table, _cram, _crai, _bam, _bai, _contamination, _vcf, _variantcaller ->
            // Get only the patient, sample and status fields from the meta map
            [meta.patient, meta.subMap('sample', 'status')]
        }
        .unique()
        .groupTuple()
        .map { patient, samples ->
            // Return the patient and the list of sample ids
            [patient, samples.collect { it.sample }]
        }
        // Flatten to [sample_id, patient] pairs
        .flatMap { patient, sample_ids -> sample_ids.collect { sample_id -> [sample_id, patient] } }
        // Group by sample_id to collect all patient ids per sample
        .groupTuple()
        .map { sample_id, patient_ids ->
            def unique_patients = patient_ids.unique()
            if (unique_patients.size() > 1) {
                System.err.println("Sample ID '${sample_id}' is associated with multiple patient IDs: ${unique_patients.join(', ')}. Please ensure each sample ID is unique to a single patient.")
                error("Execution halted due to sample status inconsistency.")
            }
        }

    // Process the input channel to group lanes by patient and sample
    // Generate patient_sample key to group lanes together
    // Save the channel ch_with_patient_sample for later use
    // Group by patient_sample to get all lanes
    // Count number of lanes per sample
    // Combine with channel ch_with_patient_sample to add numLanes information
    input_sample = ch_from_samplesheet
        .map { meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, contamination, vcf, variantcaller ->
            [meta.patient + meta.sample, [meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, contamination, vcf, variantcaller]]
        }
        .tap { ch_with_patient_sample }
        .groupTuple()
        .map { patient_sample, ch_items ->
            [patient_sample, ch_items.size()]
        }
        .combine(ch_with_patient_sample, by: 0)
        .map { _patient_sample, num_lanes, ch_items ->
            def (meta, fastq_1, fastq_2, spring_1, spring_2, table, cram, crai, bam, bai, contamination, vcf, variantcaller) = ch_items

            if (contamination) {
                meta = meta + [ contamination: contamination]
            }
            if ((meta.lane || meta.lane == 0) && fastq_2) {
                // mapping from fastq files
                meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "fastq_gz", num_lanes: num_lanes.toInteger(), size: 1]

                if (step == 'mapping') {
                    return [meta, [fastq_1, fastq_2]]
                }
                else {
                    error("Samplesheet contains fastq files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if ((meta.lane || meta.lane == 0) && spring_1 && spring_2) {
                // mapping from TWO spring-files - one with R1 and one with R2
                meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "two_fastq_gz_spring", num_lanes: num_lanes.toInteger(), size: 1]

                if (step == 'mapping') {
                    return [meta, [spring_1, spring_2]]
                }
                else {
                    error("Samplesheet contains spring files (in columns `spring_1` and `spring_2`) but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if ((meta.lane || meta.lane == 0) && spring_1 && !spring_2) {
                // mapping from ONE spring-file containing both R1 and R2
                meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), data_type: "one_fastq_gz_spring", num_lanes: num_lanes.toInteger(), size: 1]

                if (step == 'mapping') {
                    return [meta, [spring_1]]
                }
                else {
                    error("Samplesheet contains a spring file (in columns `spring_1`) but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if ((meta.lane || meta.lane == 0) && bam) {
                // Any step from BAM
                if (step != 'mapping' && !bai) {
                    error("BAM index (bai) should be provided.")
                }
                // Only use lane-based ID for mapping step where lanes are processed separately;
                // for all other steps (variant_calling, recalibrate, etc.) use sample-based ID
                // to ensure consistent ID matching in downstream channel joins
                meta = meta + [id: step == 'mapping' ? "${meta.sample}-${meta.lane}".toString() : meta.sample.toString()]
                def CN = seq_center ? "CN:${seq_center}\\t" : ''
                def read_group = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${fasta}\\tPL:${seq_platform}\""

                meta = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

                if (step != 'annotate') {
                    return [meta - meta.subMap('lane'), bam, bai]
                }
                else {
                    error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if (table && cram) {
                // recalibration from CRAM
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(step == 'mapping' || step == 'annotate')) {
                    return [meta - meta.subMap('lane'), cram, crai, table]
                }
                else {
                    error("Samplesheet contains cram files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if (table && bam) {
                // recalibration when skipping MarkDuplicates
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(step == 'mapping' || step == 'annotate')) {
                    return [meta - meta.subMap('lane'), bam, bai, table]
                }
                else {
                    error("Samplesheet contains bam files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if (cram) {
                // prepare_recalibration or variantcalling from CRAM
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(step == 'mapping' || step == 'annotate')) {
                    return [meta - meta.subMap('lane'), cram, crai]
                }
                else {
                    error("Samplesheet contains cram files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if (bam) {
                // prepare_recalibration when skipping MarkDuplicates or markduplicates
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(step == 'mapping' || step == 'annotate')) {
                    return [meta - meta.subMap('lane'), bam, bai]
                }
                else {
                    error("Samplesheet contains bam files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else if (vcf) {
                // annotation
                meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

                if (step == 'annotate') {
                    return [meta - meta.subMap('lane'), vcf]
                }
                else {
                    error("Samplesheet contains vcf files but step is `${step}`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            }
            else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }

    if (step != 'annotate' && tools && !build_only_index) {
        // Two checks for ensuring that the pipeline stops with a meaningful error message if
        // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
        // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
        input_sample
            .filter { it[0].status == 1 }
            .ifEmpty {
                // In this case, the sample-sheet contains no tumor-samples
                if (!build_only_index) {
                    def tools_tumor = ['ascat', 'controlfreec', 'mutect2', 'msisensorpro']
                    def tools_tumor_asked = []
                    tools_tumor.each { tool ->
                        if (tools && tools.split(',').contains(tool)) {
                            tools_tumor_asked.add(tool)
                        }
                    }
                    if (!tools_tumor_asked.isEmpty()) {
                        error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
                    }
                }
            }

        input_sample
            .filter { it[0].status == 0 }
            .ifEmpty {
                // In this case, the sample-sheet contains no normal/germline-samples
                def tools_requiring_normal_samples = ['ascat', 'deepvariant', 'haplotypecaller', 'msisensorpro']
                def requested_tools_requiring_normal_samples = []
                tools_requiring_normal_samples.each { tool_requiring_normal_samples ->
                    if (tools && tools.split(',').contains(tool_requiring_normal_samples)) {
                        requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
                    }
                }
                if (!requested_tools_requiring_normal_samples.isEmpty()) {
                    error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
                }
            }
    }

    // Fails when wrong extension for intervals file
    if (wes && !step == 'annotate') {
        if (intervals && !intervals.endsWith("bed")) {
            error("Target file specified with `--intervals` must be in BED format for targeted data")
        }
        else {
            log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
        }
    }
    else if (intervals && !intervals.endsWith("bed") && !intervals.endsWith("list")) {
        error("Intervals file must end with .bed, .list, or .interval_list")
    }

    if (step == 'mapping' && aligner.contains("dragmap") && !(skip_tools && skip_tools.split(',').contains("baserecalibrator"))) {
        log.warn("DragMap was specified as aligner. Base recalibration is not contained in --skip_tools. It is recommended to skip baserecalibration when using DragMap\nhttps://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode")
    }

    if (step == 'mapping' && (aligner.contains("parabricks") || aligner.contains("sentieon-bwamem")) && umi_read_structure) {
        error("${aligner} is currently not compatible with FGBio UMI handling. Please choose a different aligner.")
    }

    if (step == 'mapping' && aligner.contains("parabricks" ) && umi_in_read_header) {
        error("${aligner} is currently not compatible with extracting UMIs from read headers. Please choose a different aligner.")
    }

    if (step == 'mapping' && aligner.contains("parabricks")  && umi_location) {
        error("${aligner} is currently not compatible with UMI extraction from reads through fastp. Please choose a different aligner.")
    }

    if (step == 'mapping' && umi_read_structure && umi_location) {
        error("UMI extraction from reads through fastp (umi_location) and fgbio consensus read generation (umi_read_structure) cannot be used together. Please choose one of the two options.")
    }

    if (step == 'mapping' && umi_read_structure && umi_in_read_header) {
        // If UMIs are in read header, then we cannot use umi_read_structure separately, and instead the UMIs will be taken directly from the header
        // This requires us to set umi_read_structure to "+T +T" to indicate that UMIs are in the read header
        if( umi_read_structure != "+T +T" ) {
            error("UMI extraction from read headers (`umi_in_read_header`) will override  `umi_read_structure` when using fgbio consensus generation. Please set `umi_read_structure` to '+T +T'.")
        }
    }

    if (step == 'mapping' && umi_in_read_header && umi_location) {
        error("UMI extraction from read headers (umi_in_read_header) and UMI extraction (umi_location) from reads through fastp cannot be used together. Please choose one of the two options.")
    }

    if (step == 'mapping' && umi_location && !umi_length) {
        error("UMI extraction (umi_location) from reads through fastp requires a UMI length to be specified.")
    }

    if (tools && tools.split(',').contains("sentieon_haplotyper") && joint_germline && (!sentieon_haplotyper_emit_mode || !sentieon_haplotyper_emit_mode.contains('gvcf'))) {
        error("When setting the option `--joint_germline` and including `sentieon_haplotyper` among the requested tools, please set `--sentieon_haplotyper_emit_mode` to include `gvcf`.")
    }

    // MarkDuplicatesSpark can't do UMI based deduplication
    if (tools && tools.split(',').contains('markduplicates_spark') && (umi_in_read_header || umi_location )) {
        error("UMI based deduplication is not supported by MarkDuplicatesSpark. Please choose a different tool for deduplication.")
    }

    // Check the UMI read structure is correct using fgbio plugin
    // Copied from fastquorum
    if( umi_read_structure) {
        umi_read_structure.tokenize(" ").each { rs->
        // If parsing the read structure fails, then a java.lang.reflect.InvocationTargetException will be thrown, with
        // the cause containing the exception produced by fgbio.
        try {
            readStructure(rs)
        } catch (java.lang.reflect.InvocationTargetException ex) {
            def message = """
                |Please check the input UMI Read structure`${umi_read_structure}` invalid
                |
                |   ${ex.getCause().getMessage()}
                |
                |   For more information on read structures, visit: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
                |
                |   Validate your read structures here: https://fulcrumgenomics.github.io/fgbio/validate-read-structure.html
                |""".stripMargin()
            error(message)
            throw ex
        }
    }
    }

    // Fails or warns when missing files or params for ascat
    if (tools && tools.split(',').contains('ascat')) {
        if (!ascat_alleles) {
            error("No allele files were provided for running ASCAT. Please provide a zip folder with allele files.")
        }
        if (!ascat_loci) {
            error("No loci files were provided for running ASCAT. Please provide a zip folder with loci files.")
        }
        if (!ascat_loci_gc && !ascat_loci_rt) {
            log.warn("No LogRCorrection performed in ASCAT. For LogRCorrection to run, please provide either loci gc files or both loci gc files and loci rt files.")
        }
        if (wes) {
            log.warn("Default reference files not suited for running ASCAT on WES data. It's recommended to use the reference files provided here: https://github.com/Wedge-lab/battenberg#required-reference-files")
        }
    }

    // Warns when missing files or params for mutect2
    if (tools && tools.split(',').contains('mutect2')) {
        if (!pon) {
            log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
        }
        if (!germline_resource) {
            log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
        }
        if (pon && pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")) {
            log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
        }
    }

    // Fails when missing resources for baserecalibrator
    // Warns when missing resources for haplotypecaller
    if (!dbsnp && !known_indels) {
        if (step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!skip_tools || (skip_tools && !skip_tools.split(',').contains('baserecalibrator')))) {
            error("Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command.")
        }
        if (tools && (tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope'))) {
            log.warn("If GATK's Haplotypecaller, Sentieon's Dnascope or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-")
        }
    }
    if (joint_germline && (!tools || !(tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope')))) {
        error("The GATK's Haplotypecaller, Sentieon's Dnascope or Sentieon's Haplotyper should be specified as one of the tools when doing joint germline variant calling.) ")
    }

    if (tools && (tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope')) && joint_germline && (!dbsnp || !known_indels || !known_snps || no_intervals)) {
        log.warn(
            """If GATK's Haplotypecaller, Sentieon's Dnascope and/or Sentieon's Haplotyper is specified, but without `--dbsnp`, `--known_snps`, `--known_indels` or the associated resource labels (ie `known_snps_vqsr`), no variant recalibration will be done. For recalibration you must provide all of these resources.\nFor more information see VariantRecalibration: https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator \n\
Joint germline variant calling also requires intervals in order to genotype the samples. As a result, if `--no_intervals` is set to `true` the joint germline variant calling will not be performed."""
        )
    }

    if (tools && tools.split(',').contains('sentieon_dnascope') && joint_germline && (!sentieon_dnascope_emit_mode || !sentieon_dnascope_emit_mode.split(',').contains('gvcf'))) {
        error("When using Sentieon Dnascope for joint-germline variant-calling the option `--sentieon_dnascope_emit_mode` has to include `gvcf`.")
    }

    if (tools && tools.split(',').contains('sentieon_haplotyper') && joint_germline && (!sentieon_haplotyper_emit_mode || !sentieon_haplotyper_emit_mode.split(',').contains('gvcf'))) {
        error("When using Sentieon Haplotyper for joint-germline variant-calling the option `--sentieon_haplotyper_emit_mode` has to include `gvcf`.")
    }


    // Fails when --joint_mutect2 is used without enabling mutect2
    if (joint_mutect2 && (!tools || !tools.split(',').contains('mutect2'))) {
        error("The mutect2 should be specified as one of the tools when doing joint somatic variant calling with Mutect2. (The mutect2 could be specified by adding `--tools mutect2` to the nextflow command.)")
    }

    // Fails when missing tools for variant_calling or annotate
    if ((step == 'variant_calling' || step == 'annotate') && !tools) {
        error("Please specify at least one tool when using `--step ${step}`.\nhttps://nf-co.re/sarek/parameters#tools")
    }

    // Fails when missing sex information for CNV tools or varlociraptor
    if (tools && (tools.split(',').contains('ascat') || tools.split(',').contains('controlfreec') || tools.split(',').contains('varlociraptor'))) {
        input_sample.map {
            if (it[0].sex == 'NA') {
                error("Please specify sex information for each sample in your samplesheet when using '--tools' with 'ascat' or 'controlfreec' or 'varlociraptor'.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }
        }
    }

    // Fails when varlociraptor is enable for tumor samples but no contamination is provided
    if (tools && tools.split(',').contains('varlociraptor')) {
        input_sample.map {
            if (it[0].status == 1 && !it[0].containsKey('contamination')) {
                error("Please specify contamination information for each tumor sample in your samplesheet when using '--tools' with 'varlociraptor'.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
            }
        }
    }

    // Fails when bcftools annotate is used but no files are supplied
    if (tools && tools.split(',').contains('bcfann') && !(bcftools_annotations && bcftools_annotations_tbi && bcftools_header_lines)) {
        error("Please specify --bcftools_annotations, --bcftools_annotations_tbi, and --bcftools_header_lines, when using BCFTools annotations")
    }

    // Fails when snpeff annotation is enabled but snpeff_db is not specified
    if ((snpeff_cache && tools && (tools.split(',').contains("snpeff") || tools.split(',').contains('merge'))) && !snpeff_db) {
        error("Please specify --snpeff_db")
    }

    emit:
    input_sample
}
