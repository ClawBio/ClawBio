# rules/common.smk — helper functions shared by the rules.
# Kept out of the Snakefile so rules and functions are not mixed (snakemake --lint).


def _stage(name):
    return CFG["analysis"][name]


def _flag(value):
    return str(value).lower() if isinstance(value, bool) else str(value)
