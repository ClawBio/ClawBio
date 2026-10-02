rule qc_filter:
    """Stage 3: drop invalid / palindromic / duplicate variants (counted by reason), sort."""
    input:
        data=f"{OUT}/derive_effects/{{dataset}}.tsv",
        done=f"{OUT}/done/derive_effects_{{dataset}}.done",
    output:
        result=f"{OUT}/qc_filter/{{dataset}}.tsv",
        summary=f"{OUT}/qc_filter/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/qc_filter_{{dataset}}.done"),
    params:
        palindromic=_stage("qc_filter")["palindromic"],
        min_maf=_stage("qc_filter")["min_maf"],
        keep_indels=_flag(_stage("qc_filter")["keep_indels"]),
    log:
        f"{OUT}/logs/qc_filter/{{dataset}}.log",
    shell:
        '"{PYTHON}" "{SCRIPTS}/qc_filter.py" '
        '--input "{input.data}" --out "{output.result}" --summary-json "{output.summary}" '
        '--palindromic "{params.palindromic}" --min-maf "{params.min_maf}" '
        '--keep-indels "{params.keep_indels}" '
        '> "{log}" 2>&1'
