# Analysing isogenic human stem cells vs fibroblasts using psichomics

## Pipeline

1. [SRP063867][] loading containing biological data for human isogenic stem cells
and fibroblasts
2. Gene expression data filtering and normalisation
3. Alternative splicing quantification
4. Creation of groups of interest
5. Exploratory analysis on gene expression and alternative splicing data using
[Principal Component Analysis (PCA)][PCA]
6. Differential gene expression
7. Differential splicing
8. Correlation of alternative splicing events with gene expression of
RNA-binding proteins:
    - Test correlation using [GTEx][] data
    - Test correlation using [TCGA][] data and test for differences in survival

## File description

* [analysis.R](analysis.R): main script using psichomics to perform all the
steps of the aforementioned pipeline
* [helper_functions.R](helper_functions.R): auxiliary functions used in the main
script
* [wilcox_ties.R](wilcox_ties.R): explain the second stratum obtained during the
differential splicing analysis (Fig5b)

## Contact

For any doubt and issues, contact [nunoagostinho@medicina.ulisboa.pt][mail]

[SRP063867]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063867
[TCGA]: https://tcga-data.nci.nih.gov
[GTEx]: http://www.gtexportal.org
[mail]: mailto:nunoagostinho@medicina.ulisboa.pt
[PCA]: https://en.wikipedia.org/wiki/Principal_component_analysis