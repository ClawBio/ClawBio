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
        python=PYTHON,
        script=f"{SCRIPTS}/qc_filter.py",
        palindromic=_stage("qc_filter")["palindromic"],
        min_maf=_stage("qc_filter")["min_maf"],
        keep_indels=_flag(_stage("qc_filter")["keep_indels"]),
    conda:
        "../envs/python.yaml"
    log:
        f"{OUT}/logs/qc_filter/{{dataset}}.log",
    shell:
        # Only params/input/output/log reach the shell (snakemake --lint), each :q-quoted.
        '{params.python:q} {params.script:q} '
        '--input {input.data:q} --out {output.result:q} --summary-json {output.summary:q} '
        '--palindromic {params.palindromic:q} --min-maf {params.min_maf:q} '
        '--keep-indels {params.keep_indels:q} '
        '> {log:q} 2>&1'
