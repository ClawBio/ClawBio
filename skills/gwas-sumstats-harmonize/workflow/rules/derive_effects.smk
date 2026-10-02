rule derive_effects:
    """Stage 2: normalise values; derive BETA from OR, SE and P where missing."""
    input:
        data=f"{OUT}/map_columns/{{dataset}}.tsv",
        done=f"{OUT}/done/map_columns_{{dataset}}.done",
    output:
        result=f"{OUT}/derive_effects/{{dataset}}.tsv",
        summary=f"{OUT}/derive_effects/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/derive_effects_{{dataset}}.done"),
    log:
        f"{OUT}/logs/derive_effects/{{dataset}}.log",
    shell:
        '"{PYTHON}" "{SCRIPTS}/derive_effects.py" '
        '--input "{input.data}" --out "{output.result}" --summary-json "{output.summary}" '
        '> "{log}" 2>&1'
