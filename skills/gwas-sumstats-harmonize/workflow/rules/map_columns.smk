rule map_columns:
    """Stage 1: map raw column names onto canonical names (values copied verbatim)."""
    input:
        data=lambda wc: DATASETS[wc.dataset]["path"],
        columns=f"{OUT}/tmp/{{dataset}}.columns.json",
    output:
        result=f"{OUT}/map_columns/{{dataset}}.tsv",
        summary=f"{OUT}/map_columns/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/map_columns_{{dataset}}.done"),
    log:
        f"{OUT}/logs/map_columns/{{dataset}}.log",
    shell:
        '"{PYTHON}" "{SCRIPTS}/map_columns.py" '
        '--input "{input.data}" --columns-json "{input.columns}" '
        '--out "{output.result}" --summary-json "{output.summary}" '
        '> "{log}" 2>&1'
