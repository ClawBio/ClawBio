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
        reference_flag=f'--reference "{REFERENCE}"' if REFERENCE else "",
        drop_unmatched=_flag(_stage("align_reference")["drop_unmatched"]),
    log:
        f"{OUT}/logs/align_reference/{{dataset}}.log",
    shell:
        '"{PYTHON}" "{SCRIPTS}/align_reference.py" '
        '--input "{input.data}" {params.reference_flag} '
        '--out "{output.result}" --summary-json "{output.summary}" '
        '--drop-unmatched "{params.drop_unmatched}" '
        '> "{log}" 2>&1'
