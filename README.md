# BDE-47 source-resolved toxicogenomics

This repository contains the pre-submission Communications Biology manuscript
package for:

**Source-resolved toxicogenomics reveals compartment-conditioned BDE-47
responses in human placental evidence**

The manuscript tests a bounded claim: BDE-47-associated database sign conflicts
can reveal source-resolvable placental readout structure. In this package, CTD
is used as a discovery scaffold, PMID 31675489 is reconstructed as
placenta-derived exosomal protein cargo, PMID 37385330 is used as a bounded CTB
cellular proteome contrast, and GSE104896 is used as an mRNA boundary check.

## Main outputs

- `manuscript/bde47_commsbio_manuscript_and_supplement_v1_5.pdf`:
  combined main manuscript and Supplementary Information.
- `manuscript/bde47_commsbio_manuscript_and_supplement_v1_5.tex`:
  native LaTeX source for the combined manuscript.
- `manuscript/communications_biology_cover_letter_v1_5.pdf`:
  draft cover letter.
- `figures/`: PNG/PDF/SVG figure files.
- `tables/`: manuscript-facing derived tables and claim-boundary audits.
- `scripts/build_bde47_latex_commsbio_v1_5.py`: manuscript/figure build script.
- `provenance/latex_build_status_commsbio_v1_5.json`: build validation status.

## Evidence boundary

The compact proteome contrast is bounded directional evidence. It is not
author-level limma/FDR differential abundance, target-protein significance,
magnitude replication, or full raw-reconstructed proteome reanalysis. SNX17
retains the replicate-missingness and magnitude caveats stated in the
manuscript.

## Public data used

- CTD Chemical-Gene Interaction records: <https://ctdbase.org>
- PMID 31675489 source article and supplement
- PMID 37385330 source article, Table S1 and MassIVE MSV000087870 raw route
- GSE104896 from NCBI GEO

Large raw mass-spectrometry files are not included in this repository. The
repository contains derived manuscript-facing tables and figure source data.

## Rebuild

From the parent AlphaScience workspace, rebuild the manuscript with:

```bash
python3 scripts/build_bde47_latex_commsbio_v1_5.py
```

The script expects the local AlphaScience manuscript workspace layout used to
generate this release. The release itself is intended as a reproducibility and
review package for the manuscript-facing outputs, not as a raw-data mirror.

## DOI status

This GitHub release is the pre-submission public route referenced in the
manuscript. A DOI-bearing Zenodo archive should be minted from the release
before final journal submission.

