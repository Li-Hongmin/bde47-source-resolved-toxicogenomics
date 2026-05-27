# BDE-47 source-resolved toxicogenomics

This repository contains the Communications Biology pre-submission package for:

**Source-resolved toxicogenomics reveals compartment-conditioned BDE-47
responses in human placental evidence**

The manuscript tests a bounded biological-evidence claim: BDE-47-associated
database sign conflicts can reveal source-resolvable placental readout
structure. CTD is used as a discovery scaffold, PMID 31675489 is reconstructed
as placenta-derived exosomal protein cargo, PMID 37385330 is used as a bounded
primary villous cytotrophoblast cellular proteome contrast, and GSE104896 is
used as a cytotrophoblast mRNA boundary check.

## Main outputs

- `manuscript/bde47_commsbio_manuscript_and_supplement_v1_6.pdf`:
  combined main manuscript and Supplementary Information.
- `manuscript/bde47_commsbio_manuscript_and_supplement_v1_6.tex`:
  native LaTeX source for the combined manuscript.
- `manuscript/communications_biology_cover_letter_v1_6.pdf`:
  draft Communications Biology cover letter.
- `figures/`: submission-facing PNG, PDF, SVG and TIFF figure exports.
- `tables/`: manuscript-facing derived tables and claim-boundary audits.
- `source_data/`: audited source-output tables mirrored from the external
  analysis tree for the 31675489 exosome-cargo extraction and GSE104896 mRNA
  boundary analysis.
- `scripts/build_bde47_latex_commsbio_v1_6.py`: manuscript/figure build script.
- `provenance/latex_build_status_commsbio_v1_6.json`: build validation status.
- `docs/reproducibility.md`: command-level reproduction notes.

## Evidence boundary

The compact proteome contrast is bounded directional evidence. It is not
author-level limma/FDR differential abundance, target-protein significance,
magnitude replication or full raw-reconstructed proteome reanalysis. SNX17
retains the replicate-missingness and magnitude caveats stated in the
manuscript.

## Public data used

- CTD Chemical-Gene Interaction records: <https://ctdbase.org>
- PMID 31675489 source article and exosome proteomics supplement
- PMID 37385330 source article, Table S1 and MassIVE MSV000087870 raw route
- GSE104896 from NCBI GEO

Large raw mass-spectrometry files are not included in this repository. The
repository contains derived manuscript-facing tables, claim-boundary audits,
figure source data, audited source-output tables and build scripts.

## Rebuild

From the parent AlphaScience workspace, rebuild the manuscript with:

```bash
python3 scripts/build_bde47_latex_commsbio_v1_6.py
```

The release itself is intended as a reproducibility and review package for the
manuscript-facing outputs, not as a raw-data mirror. See
`docs/reproducibility.md` for the CTD, Table S1, DIA-NN and GSE104896
reproducibility details.

## Audited source-output mirrors

The package includes the source-output files needed to check the manuscript's
two most direct derived claims:

- `source_data/31675489/31675489_pbde47_ctd_mapped_fold_changes.csv`
  supports the 35/35 CTD-mapped PBDE47 exosome-cargo direction concordance.
- `source_data/GSE104896/` contains the external mRNA validation JSON, paired
  model summaries, matched-background permutation outputs and gene-level
  background table used for the mRNA boundary analysis.

## Release and DOI status

Target GitHub release tag:

```text
v1.0-commsbio
```

Zenodo metadata is prepared in `.zenodo.json`. A DOI-bearing Zenodo archive
still requires the GitHub release to be connected to Zenodo or manually
uploaded by an authenticated user. No DOI is claimed until Zenodo returns one.
