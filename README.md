# Analysing isogenic human stem cells vs fibroblasts using psichomics

This project contains the code for the analysis performed in the following
book chapter:

> Nuno Saraiva-Agostinho and Nuno L. Barbosa-Morais (2020). 
**[Interactive Alternative Splicing Analysis of Human Stem Cells Using psichomics][chapter]**. In: Kidder B. (eds) Stem Cell Transcriptional Networks. *Methods in Molecular Biology*, vol 2117. Humana, New York, NY

## Pipeline

1. Loading [SRP063867][] data, containing data for human isogenic stem cells
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

| File                                     | Description |
| ---------------------------------------- | ----------- |
| [analysis.R](analysis.R)                 | Main script using psichomics to perform the aforementioned pipeline |
| [helper_functions.R](helper_functions.R) | Auxiliary functions used in the main script |
| [wilcox_ties.R](wilcox_ties.R)           | Explanation of the second stratum obtained during the differential splicing analysis (Fig5b) |

## Contact

In case of doubts and issues, please contact [nunoagostinho@medicina.ulisboa.pt][mail]

[chapter]: https://doi.org/10.1007/978-1-0716-0301-7_10
[SRP063867]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063867
[TCGA]: https://tcga-data.nci.nih.gov
[GTEx]: http://www.gtexportal.org
[mail]: mailto:nunoagostinho@medicina.ulisboa.pt
[PCA]: https://en.wikipedia.org/wiki/Principal_component_analysis
