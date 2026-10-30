# Real tool outputs (regression fixtures)

Unlike `examples/`, these files are real tool output. Each was produced by running the
named tool on the public RefSeq assembly GCF_000493755.1 (*Escherichia coli* JJ1886).
They pin the column layouts the parsers depend on.

| File | Tool and version | Command |
|---|---|---|
| `jj1886_amrfinder_v3.tsv` | AMRFinderPlus 3.12.8, database 2024-07-22.1 | `amrfinder --nucleotide <fna> --plus --organism Escherichia` |
| `jj1886_amrfinder_v4.tsv` | AMRFinderPlus 4.2.7, database 2026-08-07.1 | same |
| `jj1886_mlst.tsv` | mlst 2.35.0 | `mlst <fna>` |
| `jj1886_plasmidfinder_v2.tsv` | PlasmidFinder 2.1.6 `results_tab.tsv` | `plasmidfinder.py -i <fna> -o <dir> -x -l 0.60 -t 0.95` |
| `jj1886_plasmidfinder_v3.json` | PlasmidFinder 3.0.3 `-j` JSON | `python -m plasmidfinder -i <fna> -o <dir> -l 0.60 -t 0.95 -p <db> -j <json>` |
| `jj1886_abricate_plasmidfinder.tsv` | abricate 1.4.0, default thresholds | `abricate --db plasmidfinder <fna>` |

Edits: the `genomes/` path prefix was removed from the FILE column of the mlst and abricate
files, and the alignment strings and the `software_executions` block (which held local
paths) were removed from the PlasmidFinder 3 JSON.
