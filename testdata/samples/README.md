# Ready-made test samples

Four paired-end samples, gzipped (1.3 MB total), 5 000 pairs each. Drop a pair straight into the
app — `_R1`/`_R2` are grouped automatically and `.gz` is read natively.

| Sample | Content | Human index | rRNA index |
|---|---|---|---|
| `sample_with_human` | real human pairs (GIAB HG002, 125 bp) | **94.1 %** | 0.3 % |
| `sample_without_human` | random DNA + *E. coli* 16S | **0.0 %** | 50.0 % |
| `sample_with_rrna` | *E. coli* 16S + human 18S | 50.0 % | **100.0 %** |
| `sample_without_rrna` | random DNA only | **0.0 %** | **0.0 %** |

Percentages are reads removed at threshold 0.5, measured against the published indexes.

**The off-diagonal figures are correct, not bugs.** `with_rrna` loses half to the human index
because half of it is the human 18S gene, which genuinely belongs to the human genome.
`without_human` loses half to the rRNA index because its *E. coli* 16S half is genuinely ribosomal.
`with_human` loses 0.3 % to the rRNA index because real human reads sometimes fall in ribosomal
regions — biology, not a false positive.

Negatives are random DNA on purpose: a real bacterial genome carries rRNA operons and would quietly
break the `without_rrna` guarantee.

## Quickest check

Load the **Ribosomal RNA** index (282 MB, downloads in seconds), then run `sample_with_rrna`
(expect ~100 % removed) and `sample_without_rrna` (expect 0 %). If both land, the whole chain works.

## Provenance

- Human reads: [GIAB](https://www.nist.gov/programs-projects/genome-bottle) HG002 (NIST, public
  reference data), NIST_Stanford_Illumina_6kb_matepair.
- *E. coli* 16S: GenBank `J01859.1` · human 18S: RefSeq `NR_003286.4`.
- Regenerate with `python3 ../gen_testsets.py <dir with hg002_R1/R2.fastq>`.
