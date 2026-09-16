# Round-5 DnaSP 6.12.03 evidence used by the tests

These files record what DnaSP 6.12.03 for Windows printed during the round-5
comparison (15 September 2026). Tests read them as recorded evidence; no test
runs DnaSP.

| File | What it is | Used by |
|------|------------|---------|
| `A2_rp49_win100_25.out` | DnaSP output file saved by DnaSP: rp49 ingroup, sliding window 100 / step 25 | `test_validation5_regressions.py` window placement and midpoint tests |
| `C3_COII_win100_50.out` | DnaSP output file saved by DnaSP: COII, sliding window 100 / step 50 | same |
| `E2_Dmel_win500_250.out` | DnaSP output file saved by DnaSP: *D. melanogaster* Os region, sliding window 500 / step 250 | same |
| `M_ENC_standard_GUI_tables_transcribed.txt` | **Transcription, not a DnaSP export.** Per-sequence ENC table under the standard code, transcribed by Claude Code from the operator's screenshots | `test_per_sequence_enc_transcriptions` |
| `M_ENC_vertebrate-mitochondrial_GUI_tables_transcribed.txt` | **Transcription, not a DnaSP export.** Per-sequence ENC table under the vertebrate mitochondrial code, transcribed by Claude Code from the operator's pasted text and screenshots | same |

DnaSP 6.12.03 could not export the per-sequence codon-usage table when the
coding assignment has no noncoding region: the output file stops after the table
header and DnaSP reports "Error in Codon Usage module". The two ENC tables were
therefore transcribed from the screen. They are labelled as transcriptions in
their headers, and the comparisons that rely on them are flagged as transcribed
in Supplementary Table S2 of the DnaSP-ClawBio manuscript.

DnaSP abbreviates local input paths in its output as `C:\...`; no other local
information is recorded. The input alignments are described in
`../inputs/README.md`.
