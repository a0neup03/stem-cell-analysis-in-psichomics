# Using psichomics to analyse human isogenic stem cells against fibroblasts

1. SRP063867 data loading
2. Gene expression data filtering and normalisation
3. Alternative splicing quantification
4. Creation of groups of interest
5. Exploratory analysis on gene expression and alternative splicing data using
Principal Component Analysis (PCA)
6. Differential gene expression
7. Differential splicing
8. Correlation of alternative splicing events with gene expression of RBPs:
    - Test correlation using [GTEx]() data
    - Test correlation using [TCGA]() data and test for differences in survival

* [analysis.R](): main script using psichomics to analyse human isogenic stem
cells against fibroblasts ([SRP063867][SRA]), where the aforementioned steps are
performed step-by-step
* [helper_functions.R](): auxiliary functions used in the main script
* [wilcox_ties.R](): explain the second stratum obtained during the differential
splicing analysis (Fig5b)

For any doubt/issues, contact [nunoagostinho@medicina.ulisboa.pt][mail]

[SRA]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063867
[TCGA]: https://tcga-data.nci.nih.gov
[GTEx]: http://www.gtexportal.org
[mail]: mailto:nunoagostinho@medicina.ulisboa.pt