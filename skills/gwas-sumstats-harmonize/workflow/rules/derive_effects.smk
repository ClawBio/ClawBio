rule derive_effects:
    """Stage 2: normalise values; derive BETA from OR, SE and P where missing."""
    input:
        data=f"{OUT}/map_columns/{{dataset}}.tsv",
        done=f"{OUT}/done/map_columns_{{dataset}}.done",
    output:
        result=f"{OUT}/derive_effects/{{dataset}}.tsv",
        summary=f"{OUT}/derive_effects/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/derive_effects_{{dataset}}.done"),
    params:
        python=PYTHON,
        script=f"{SCRIPTS}/derive_effects.py",
    conda:
        "../envs/python.yaml"
    log:
        f"{OUT}/logs/derive_effects/{{dataset}}.log",
    shell:
        # Only params/input/output/log reach the shell (snakemake --lint), each :q-quoted.
        '{params.python:q} {params.script:q} '
        '--input {input.data:q} --out {output.result:q} --summary-json {output.summary:q} '
        '> {log:q} 2>&1'
