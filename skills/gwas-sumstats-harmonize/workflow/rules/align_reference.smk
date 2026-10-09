rule align_reference:
    """Stage 4: align alleles so EA = reference ALT; fill EAF allele-matched."""
    input:
        data=f"{OUT}/qc_filter/{{dataset}}.tsv",
        done=f"{OUT}/done/qc_filter_{{dataset}}.done",
        reference=[REFERENCE] if REFERENCE else [],
    output:
        result=f"{OUT}/align_reference/{{dataset}}.tsv",
        summary=f"{OUT}/align_reference/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/align_reference_{{dataset}}.done"),
    params:
        python=PYTHON,
        script=f"{SCRIPTS}/align_reference.py",
        drop_unmatched=_flag(_stage("align_reference")["drop_unmatched"]),
    conda:
        "../envs/python.yaml"
    log:
        f"{OUT}/logs/align_reference/{{dataset}}.log",
    shell:
        # Only params/input/output/log reach the shell (snakemake --lint), each :q-quoted.
        # With no reference, input.reference formats to nothing (not ''), so the
        # --reference=<value> form is needed: it gives "" and the script passes rows through.
        '{params.python:q} {params.script:q} '
        '--input {input.data:q} --reference={input.reference:q} '
        '--out {output.result:q} --summary-json {output.summary:q} '
        '--drop-unmatched {params.drop_unmatched:q} '
        '> {log:q} 2>&1'
