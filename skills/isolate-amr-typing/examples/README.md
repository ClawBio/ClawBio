# Demo data (synthetic)

These files drive `--demo` and the test suite. None of them comes from a real isolate.

| File | What it is |
|---|---|
| `demo_isolate.fasta` | Four contigs of seeded random sequence (20.3 kb). Used only for assembly statistics and contig-name checks. |
| `demo_amrfinder.tsv` | Hand-written table in AMRFinderPlus v4 column layout. |
| `demo_mlst.tsv` | Hand-written line in `mlst` default output layout. |
| `demo_plasmidfinder.tsv` | Hand-written table in PlasmidFinder `results_tab.tsv` layout. |

The tool-output files were written by hand to resemble an ST131-like *Escherichia coli*
genotype. They were **not** produced by running the tools on `demo_isolate.fasta`, which is
random sequence and contains no genes. Coordinates, lengths, identities and accessions
(`SYNTHETIC`) are placeholders.
