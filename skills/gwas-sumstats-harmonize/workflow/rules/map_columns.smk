rule map_columns:
    """Stage 1: map raw column names onto canonical names (values copied verbatim)."""
    input:
        data=lambda wc: DATASETS[wc.dataset]["path"],
        columns=f"{OUT}/tmp/{{dataset}}.columns.json",
    output:
        result=f"{OUT}/map_columns/{{dataset}}.tsv",
        summary=f"{OUT}/map_columns/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/map_columns_{{dataset}}.done"),
    params:
        python=PYTHON,
        script=f"{SCRIPTS}/map_columns.py",
    conda:
        "../envs/python.yaml"
    log:
        f"{OUT}/logs/map_columns/{{dataset}}.log",
    shell:
        # Only params/input/output/log reach the shell (snakemake --lint), each :q-quoted.
        '{params.python:q} {params.script:q} '
        '--input {input.data:q} --columns-json {input.columns:q} '
        '--out {output.result:q} --summary-json {output.summary:q} '
        '> {log:q} 2>&1'
