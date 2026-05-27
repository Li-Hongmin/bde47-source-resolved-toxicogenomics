# Reproducibility notes

This file records the command-level route used for the manuscript-facing
outputs. It is not a raw-data mirror.

## Manuscript and figures

From the parent AlphaScience workspace:

```bash
python3 scripts/build_bde47_latex_commsbio_v1_6.py
```

The build produces:

- `manuscript/bde47_commsbio_manuscript_and_supplement_v1_6.pdf`
- `manuscript/communications_biology_cover_letter_v1_6.pdf`
- `figures/*.png`
- `figures/*.pdf`
- `figures/*.svg`
- `figures/*.tiff`
- `provenance/latex_build_status_commsbio_v1_6.json`

## CTD sign-tension scaffold

The CTD scaffold used the local run identifier `20260523T105644Z` and was
bounded to the first 2,000,000 non-comment CTD Chemical-Gene Interaction rows.
Rows were retained only when they had:

- a PubMed identifier;
- direct chemical-gene evidence;
- an increase or decrease direction in expression/activity action classes;
- a source-readiness pass flag;
- a chemical-to-gene orientation sentence.

The orientation gate was deliberately fail-closed. Obvious gene/protein-to-
chemical metabolism or activation sentences were excluded. The BDE-47
background audit used a 20-chemical row-count-nearest background and a
10,000-iteration BDE-47 sign-label permutation with seed 47.

## PMID 31675489 exosome-cargo mapping

The 31675489 block was mapped from CTD gene symbols to the PBDE47 supplement
table. Listed zero fold-change values were retained as reported and floored
only for log-scale display. Direction concordance was computed from the
reported fold-change direction, not from the display-floor values.

Source table used by the v1.6 build:

```text
/Volumes/DevWork/AlphaScience_external/contradiction-driven-discovery/outputs/supplement_tables/20260523T105644Z/31675489_pbde47_ctd_mapped_fold_changes.csv
```

The same source-output table is mirrored in this release at:

```text
source_data/31675489/31675489_pbde47_ctd_mapped_fold_changes.csv
```

## PMID 37385330 Table S1 contrast

Cellular proteome contrasts were computed from positive Table S1 relative
abundance values as:

```text
log2(mean BDE47 condition / mean matched-time vehicle)
```

No pseudocount was added for retained proteins. The six manuscript-facing
BDE-47 dose-time conditions were:

- 1 uM 15 h
- 5 uM 15 h
- 1 uM 24 h
- 5 uM 24 h
- 1 uM 39 h
- 5 uM 39 h

## DIA-NN raw pilot

The bounded raw pilot used MassIVE `MSV000087870`, BDE-47 5 uM 39 h versus
vehicle 39 h. SCIEX WIFF/.scan files were converted in a Windows-compatible
route to seven `.wiff.dia` inputs: three BDE-47 5 uM 39 h files and four
vehicle 39 h files.

DIA-NN 1.8.1 was run in library-free mode using a human UniProt SwissProt FASTA
downloaded on 2026-05-24. The command structure was:

```bash
DiaNN.exe \
  --f <bde47_rep1.wiff.dia> \
  --f <bde47_rep2.wiff.dia> \
  --f <bde47_rep3.wiff.dia> \
  --f <vehicle_rep1.wiff.dia> \
  --f <vehicle_rep2.wiff.dia> \
  --f <vehicle_rep3.wiff.dia> \
  --f <vehicle_rep4.wiff.dia> \
  --fasta <human_uniprot_swissprot_2026-05-24.fasta> \
  --fasta-search \
  --predictor \
  --matrices \
  --out report.tsv
```

Target extraction used `report.pg_matrix.tsv`. Missing or zero protein-group
quantities were left missing for log2 summaries and were not imputed. Target
p-values are descriptive exact target-only permutations and are not global
proteome FDR or author-level limma/FDR.

## GSE104896 mRNA boundary analysis

GSE104896 was used as an independent primary villous cytotrophoblast mRNA
boundary check. Processed GEO expression records were mapped to gene symbols
using platform-based probe annotation. Donor-paired BDE-47 versus vehicle
deltas were used where pairing metadata were available.

Matched-background variables:

- selected-probe mean expression;
- expression standard deviation;
- paired-delta standard deviation;
- probe count.

The gene-level boundary used donor-paired exact sign-flip testing with
all-background Benjamini-Hochberg FDR.

The GSE104896 source-output files used for manuscript-facing checks are
mirrored in this release under:

```text
source_data/GSE104896/
```

Key files:

- `GSE104896_external_validation.json`
- `GSE104896_ctb_mrna_validation_rows.csv`
- `GSE104896_matched_background_permutation.csv`
- `GSE104896_paired_model_summary.csv`
- `GSE104896_gene_background_log2_delta.csv`
