# GSE104896 External CTB mRNA Validation

Run ID: `20260524T034100Z`
Generated UTC: `2026-05-24T06:09:29.722900+00:00`

## Boundary

This is an independent context/readout check, not an exact replication of the 31675489 exosome-cargo experiment. GSE104896 measures RMA-normalized mRNA in primary villous cytotrophoblasts exposed to 1 uM BDE-47 for 24 h. PMID 31675489 measures placenta-derived exosomal protein cargo after PBDE-47 exposure, and PMID 37385330 measures CTB proteome abundance. The analysis below asks whether the 31675489 exosome decrease block is mirrored by a uniform decrease in an independent cellular mRNA context.

## Core Decision

The external CTB mRNA data do not show a uniform cellular mRNA decrease for the 31675489 exosome-cargo decrease block. All 35 CTD-mapped genes were mapped from Affymetrix transcript-cluster IDs to gene symbols. Among the 32 genes that decrease in exosome cargo, 21 have positive mean paired CTB mRNA deltas and 10 have negative deltas after excluding one unchanged gene. Against the GSE104896 background positive fraction of 0.497, this positive-skew check has a descriptive binomial upper-tail p value of 0.0330.

This supports the compartment/readout-divergence model at the level of directional consistency: exosome cargo depletion should not be interpreted as equivalent to intracellular or cellular mRNA depletion. It does not yet prove a universal BDE-47 sign law or replace the need for author-level differential-expression statistics.

The paired-model upgrade uses a donor-paired exact sign-flip test on BDE-47 minus vehicle deltas for each selected gene probe, followed by Benjamini-Hochberg FDR correction across the mapped GSE104896 gene background. This is a Python-native paired model, not limma; it is included to formalize the external mRNA check without adding an unavailable R/Bioconductor runtime.

The matched-background permutation adds a stricter background check using selected-probe mean expression, selected-probe expression SD, paired-delta SD, and probe count. This still supports the same bounded conclusion: the exosome-cargo decrease block is not a simple uniform cellular mRNA-decrease block.

## Summary

- CTD-mapped genes checked: `35`
- GSE104896 background genes with mapped probes: `26932`
- Exosome decrease genes: `32`
- Exosome decrease genes with positive CTB mRNA mean delta: `21`
- Exosome decrease genes with negative CTB mRNA mean delta: `10`
- Exosome decrease genes with unchanged CTB mRNA mean delta: `1`
- Correlation between exosome log2FC and CTB mRNA mean delta (zero FC floored for plotting): `-0.016`
- Paired model: `donor_paired_exact_signflip_mean_delta`
- Paired model BH universe: `26932` genes
- CTD-mapped genes with paired-model FDR < 0.05: `0`
- Exosome-decrease genes with paired-model FDR < 0.05: `0`

## Key Genes

- `ILK`: exosome decrease (FC=0.7), CTB mRNA increase (mean log2 delta=0.073)
- `NFKB1`: exosome decrease (FC=0.0), CTB mRNA decrease (mean log2 delta=-0.100)
- `SNX17`: exosome decrease (FC=0.3), CTB mRNA increase (mean log2 delta=0.012)

## Source Tables

| Artifact | URL | Local path | Bytes |
|---|---|---|---:|
| GSE104896_series_matrix.txt.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104896/matrix/GSE104896_series_matrix.txt.gz | outputs/external_omics/GSE104896/raw/GSE104896_series_matrix.txt.gz | 1170614 |
| GPL16686_full.txt | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16686&targ=self&form=text&view=full | outputs/external_omics/GSE104896/raw/GPL16686_full.txt | 4043207 |
| hugene20sttranscriptcluster.db_8.8.0.tar.gz | https://bioconductor.org/packages/3.23/data/annotation/src/contrib/hugene20sttranscriptcluster.db_8.8.0.tar.gz | outputs/external_omics/GSE104896/raw/hugene20sttranscriptcluster.db_8.8.0.tar.gz | 4928528 |
| Homo_sapiens.gene_info.gz | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | outputs/external_omics/GSE104896/raw/Homo_sapiens.gene_info.gz | 5147322 |

## CTD-Mapped Gene Deltas

| Gene | Exosome FC | Exosome sign | CTB mRNA mean log2 delta | CTB mRNA sign | Paired p | BH FDR | Positive pairs | Negative pairs |
|---|---:|---|---:|---|---:|---:|---:|---:|
| AARS1 | 0.7 | decrease | 0.072 | increase | 0.0938 | 0.7357 | 5 | 1 |
| ARHGEF7 | 0.5 | decrease | 0.000 | unchanged | 1.0000 | 1.0000 | 2 | 4 |
| CASP1 | 0 | decrease | 0.267 | increase | 0.2500 | 0.8449 | 3 | 2 |
| CDK2 | 0.4 | decrease | -0.073 | decrease | 0.2188 | 0.8409 | 1 | 5 |
| CYP11A1 | 0.6 | decrease | 0.033 | increase | 0.4688 | 0.9271 | 4 | 2 |
| FCN1 | 2 | increase | -0.128 | decrease | 0.3125 | 0.8734 | 2 | 4 |
| FCN2 | 2.5 | increase | 0.098 | increase | 0.6250 | 0.9605 | 3 | 3 |
| FIS1 | 0.5 | decrease | 0.055 | increase | 0.2188 | 0.8409 | 4 | 2 |
| FNIP2 | 0 | decrease | 0.315 | increase | 0.0312 | 0.6653 | 6 | 0 |
| FUS | 0 | decrease | -0.085 | decrease | 0.1250 | 0.7501 | 2 | 4 |
| GFPT1 | 0.7 | decrease | 0.013 | increase | 0.8125 | 0.9896 | 3 | 2 |
| GMPPB | 0.7 | decrease | 0.083 | increase | 0.2500 | 0.8449 | 4 | 2 |
| IFIT1 | 0.5 | decrease | 0.250 | increase | 0.3438 | 0.8901 | 4 | 2 |
| ILK | 0.7 | decrease | 0.073 | increase | 0.1250 | 0.7501 | 4 | 1 |
| LGMN | 0.7 | decrease | -0.052 | decrease | 0.0625 | 0.6793 | 0 | 5 |
| MARS1 | 0.7 | decrease | 0.110 | increase | 0.0625 | 0.6793 | 5 | 1 |
| MGLL | 0.3 | decrease | -0.013 | decrease | 0.8750 | 0.9966 | 3 | 3 |
| MMP19 | 0.7 | decrease | 0.098 | increase | 0.1562 | 0.7974 | 5 | 1 |
| MMRN1 | 2 | increase | -0.077 | decrease | 0.5312 | 0.9430 | 2 | 4 |
| MX1 | 0.6 | decrease | 0.300 | increase | 0.0938 | 0.7357 | 5 | 1 |
| MYO1E | 0 | decrease | -0.162 | decrease | 0.0312 | 0.6653 | 0 | 6 |
| NFKB1 | 0 | decrease | -0.100 | decrease | 0.2500 | 0.8449 | 1 | 5 |
| OAS3 | 0.3 | decrease | 0.222 | increase | 0.1562 | 0.7974 | 5 | 1 |
| PARP14 | 0 | decrease | 0.087 | increase | 0.7500 | 0.9794 | 2 | 3 |
| PCNA | 0.3 | decrease | 0.037 | increase | 0.5625 | 0.9500 | 3 | 3 |
| PGM3 | 0.6 | decrease | 0.097 | increase | 0.0625 | 0.6793 | 5 | 1 |
| PRKDC | 0.7 | decrease | 0.023 | increase | 0.7812 | 0.9891 | 4 | 2 |
| SEC24A | 0.7 | decrease | -0.010 | decrease | 0.8125 | 0.9896 | 2 | 4 |
| SEC24D | 0.6 | decrease | 0.018 | increase | 0.5938 | 0.9597 | 3 | 3 |
| SNRPA1 | 0.7 | decrease | -0.033 | decrease | 0.5938 | 0.9597 | 3 | 3 |
| SNX17 | 0.3 | decrease | 0.012 | increase | 0.9062 | 1.0000 | 3 | 3 |
| TGFB1 | 0 | decrease | -0.122 | decrease | 0.0625 | 0.6793 | 1 | 5 |
| VWA5A | 0.5 | decrease | -0.147 | decrease | 0.4688 | 0.9271 | 3 | 3 |
| WIPI1 | 0.7 | decrease | 0.060 | increase | 0.3125 | 0.8734 | 4 | 2 |
| YARS1 | 0.5 | decrease | 0.122 | increase | 0.0312 | 0.6653 | 6 | 0 |

## Module Summary

| Module | Genes | Mean CTB mRNA log2 delta | Positive | Negative |
|---|---|---:|---:|---:|
| immune_antiviral_inflammatory | `CASP1|IFIT1|MX1|NFKB1|OAS3|PARP14|TGFB1` | 0.129 | 5 | 2 |
| vesicle_trafficking_sorting | `ARHGEF7|ILK|MYO1E|SEC24A|SEC24D|SNX17` | -0.011 | 3 | 2 |
| cell_cycle_dna_repair | `CDK2|PCNA|PRKDC` | -0.004 | 2 | 1 |
| metabolism_translation_glycosylation | `AARS1|CYP11A1|GFPT1|GMPPB|MARS1|MGLL|PGM3|YARS1` | 0.065 | 7 | 1 |
| autophagy_stress_organelle | `FIS1|FNIP2|FUS|LGMN|WIPI1` | 0.059 | 3 | 2 |
| matrix_invasion_signaling | `MMP19|VWA5A` | -0.024 | 1 | 1 |
| rna_processing | `SNRPA1` | -0.033 | 0 | 1 |

## Matched-Background Permutation

| Metric | Observed | Null mean | Null SD | P(null >= obs) | P(null <= obs) |
|---|---:|---:|---:|---:|---:|
| positive_fraction | 0.6774 | 0.5045 | 0.0880 | 0.0301 | 0.9720 |
| negative_fraction | 0.3226 | 0.4955 | 0.0880 | 0.9720 | 0.0301 |
| mean_log2_delta | 0.0484 | 0.0053 | 0.0190 | 0.0122 | 0.9879 |

## Paired Model Summary

| Metric | Value | Interpretation |
|---|---:|---|
| model | donor_paired_exact_signflip_mean_delta | Exact paired sign-flip test on BDE-47 minus vehicle donor deltas |
| background_genes_bh_universe | 26932 | All mapped GSE104896 genes with selected probes |
| ctd_mapped_genes | 35 | 31675489 CTD-mapped genes checked in GSE104896 |
| exosome_decrease_genes | 32 | 31675489 exosome-cargo decrease block |
| exosome_decrease_positive_nonzero_fraction | 0.6774193548387096 | Positive CTB mRNA paired-model logFC fraction |
| ctd_mapped_genes_fdr_lt_0_05 | 0 | CTD-mapped genes passing all-background BH FDR |
| exosome_decrease_genes_fdr_lt_0_05 | 0 | Exosome decrease genes passing all-background BH FDR |
| background_genes_fdr_lt_0_05 | 0 | All-background genes passing BH FDR |

## Interpretation Boundary

The useful result is negative in a good way: an independent CTB mRNA dataset fails to support the simplistic interpretation that the 31675489 exosome-cargo decrease block is just cellular expression decrease. The external dataset is directionally compatible with a compartment/readout split, but exact mechanism still rests on the 31675489 exosome table and the 37385330 CTB proteome table.
