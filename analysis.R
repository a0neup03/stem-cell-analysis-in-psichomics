# Isogenic stem cell vs isogenic fibroblast
# Nuno Agostinho, 19 Feb 2019

source("helper_functions.R")
library(psichomics)

# Load SRA data ----------------------------------------------------------------
data <- loadSRAproject("SRP063867")
clinical      <- data[[1]]$`Clinical data`
sampleInfo    <- data[[1]]$`Sample metadata`
junctionQuant <- data[[1]]$`Junction quantification`
geneExpr      <- data[[1]]$`Gene expression`

# Filter and normalise gene expression data ------------------------------------
plotGeneExprPerSample(geneExpr)
plotRowStats(geneExpr, "mean", "log10(var)")

fig3a <- plotDistribution(log10(colSums(geneExpr)), rugLabels=TRUE,
                          vLine=FALSE) %>%
    hc_xAxis(title=list(text="log10(Library sizes)")) %>%
    hc_yAxis(title=list(text="Density")) %>%
    hc_legend(enabled=FALSE) %>%
    hc_exporting(enabled=FALSE)

# Filtering using function edgeR::filterByExpr()
geneExprNorm <- list()
geneExprNorm$filterByExpr <- prepareGeneExpr(geneExpr)
plotDistribution(colSums(geneExprNorm$filterByExpr), rugLabels=TRUE)

fig3b <- plotGeneExprPerSample(
    geneExprNorm$filterByExpr, sortByMedian=FALSE, showXlabels=FALSE,
    title="Low read count filtering and raw library scaling") %>%
    hc_yAxis(title=list(text="log2CPM")) %>%
    hc_exporting(enabled=FALSE)
# Bad normalisation...

# What if we use voom?
geneExprNorm$voom <- prepareGeneExpr(geneExpr, performVoom=TRUE)
plotDistribution(colSums(geneExprNorm$voom), rugLabels=TRUE)

fig3c <- plotGeneExprPerSample(
    geneExprNorm$voom, sortByMedian=FALSE, showXlabels=FALSE,
    title="Low read count filtering, raw library scaling and voom") %>%
    hc_yAxis(title=list(text="log2CPM")) %>%
    hc_exporting(enabled=FALSE)

# Applying Voom with quantile normalisation
geneExprNorm$voomQuantile <- prepareGeneExpr(geneExpr, performVoom=TRUE,
                                             method="quantile")
plotDistribution(colSums(geneExprNorm$voomQuantile), rugLabels=TRUE)

fig3d <- plotGeneExprPerSample(
    geneExprNorm$voomQuantile, sortByMedian=FALSE, showXlabels=FALSE,
    title="Low read count filtering and voom using quantile normalisation") %>%
    hc_yAxis(title=list(text="log2CPM")) %>%
    hc_exporting(enabled=FALSE)

# Can we improve it by removing bad samples?
plotDistribution(log10(colSums(geneExpr)), rugLabels=TRUE)
samplesToDiscard <- log10(colSums(geneExpr)) < 8
samplesToDiscard <- names(samplesToDiscard[samplesToDiscard])

geneExprFiltered <- geneExpr[ , !colnames(geneExpr) %in% samplesToDiscard]
geneExprNorm$voomQuantileDiscardBadSamples <- prepareGeneExpr(
    geneExprFiltered, performVoom=TRUE, method="quantile")

plotDistribution(colSums(geneExprNorm$voomQuantileDiscardBadSamples),
                 rugLabels=TRUE)

fig3e <- plotGeneExprPerSample(
    geneExprNorm$voomQuantileDiscardBadSamples, sortByMedian=FALSE,
    showXlabels=FALSE,
    title=paste("Discarded SRR2453313, followed by low read count filtering",
                "and voom using quantile normalisation")) %>%
    hc_yAxis(title=list(text="log2CPM")) %>%
    hc_exporting(enabled=FALSE)

fig3 <- browsable(hw_grid(fig3a, fig3b, fig3c, fig3d, fig3e,
                          ncol=1, rowheight="150px"))
fig3 # Missing labels/annotation of SRR2453313
save_html(fig3, "Fig3.html")

# Alternative splicing quantification ------------------------------------------
# Load Human (hg38 assembly) annotation
human <- grep("hg38", listSplicingAnnotations(), value=TRUE)
annotation <- loadAnnotation(human)

junctionQuantFiltered <- junctionQuant[
    , !colnames(junctionQuant) %in% samplesToDiscard]
psi <- list()
psi$psi <- quantifySplicing(annotation, junctionQuantFiltered)

library(miscTools)
fig4a <- plotRowStats(psi$psi, x="median", y="var", xmin=0.05, xmax=0.95) +
    xlab("PSI median") +
    ylab("PSI variance") +
    theme_classic()

plotRangeVarWithGreyPoints <- function(allPSI, filteredPSI) {
    rangeVarPSIplot <- plotRowStats(
        allPSI, x="range", y="log10(var)", xmin=0.15, ymin=-3,
        ylim=c(-15, 0))
    rangeVarPSIplot <- rangeVarPSIplot +
        geom_point(aes(x=range, y=log10(var)), colour="grey")
    return(rangeVarPSIplot)
}
plotRangeVarWithGreyPoints(psi$psi, psi$filteredByMedian)

psi$filteredByMedian <- psi$psi[filterPSI(
    psi$psi, minMedian=0.05, maxMedian=0.95), ]
# fig4b <- plotRangeVarWithGreyPoints(psi$psi, psi$filteredByMedian) +
fig4b <- plotRowStats(psi$filteredByMedian, "range", "log10(var)",
                      ylim=c(-15, 0), xmin=0.15, ymin=-3) +
    xlab("PSI range") +
    ylab("log10(PSI variance)") +
    theme_classic()

fig4 <- plot_grid(fig4a, fig4b, labels=c("a", "b"))
fig4
saveToTIFF(fig4, width=10, filepath="Fig4.tiff")

filteredPSI <- psi$psi[filterPSI(psi$psi, minMedian=0.05, maxMedian=0.95,
                                 minLogVar=-3, minRange=0.15), ]

# Create groups ----------------------------------------------------------------
groups <- prepareGroups(sampleInfo)

# PCA on gene expression -------------------------------------------------------
round_perc_in_string <- function(str) {
    perc <- gsub(".*?([0-9]{1,}\\.[0-9]*)%.*", "\\1", str)
    gsub("[0-9]{1,}\\.[0-9]*", round(as.numeric(perc)), str)
}

pca <- list()
pca$geneExprNorm <- performPCA(t(geneExprNorm$voomQuantileDiscardBadSamples$E),
                               scale.=TRUE)
fig5a <- plotVariance(pca$geneExprNorm) %>%
    hc_yAxis(title=list(text="Cumulative percentage of variance")) %>%
    hc_plotOptions(series=list(dataLabels=list(enabled=FALSE)))
    
fig5b <- plotPCA(pca$geneExprNorm, groups=groups$cellType)
fig5c <- plotPCA(pca$geneExprNorm, loadings=TRUE, individuals=FALSE,
                 nLoadings=100) %>% hc_subtitle(text="")
    # hc_subtitle(text="Bubble size ~ relative contribution to PC1 and PC2")
fig5d <- plotPCA(pca$geneExprNorm, groups=groups$dataSet)

head(calculateLoadingsContribution(pca$geneExprNorm))

# PCA on alternative splicing quantification -----------------------------------
pca$filteredPSI <- performPCA(t(filteredPSI))
# plotVariance(pca)
fig5e <- plotPCA(pca$filteredPSI, groups=groups$cellType)
fig5f <- plotPCA(pca$filteredPSI, groups=groups$dataSet)
# plotPCA(pca, loadings=TRUE, individuals=FALSE)
fig5g <- calculateLoadingsContribution(pca$filteredPSI)

prepareChartForPublication <- function(
    hc, filename="chart", size=16, title=size, subtitle=size,
    axisTitle=size, axisLabel=size, dataLabel=size, legend=size) {
    
    # Increase font size
    hc <- hc %>%
        hc_title(text="") %>%
        hc_subtitle(style=list(fontSize=subtitle)) %>%
        hc_xAxis(title=list(style=list(fontSize=axisTitle)),
                 labels=list(style=list(fontSize=axisLabel))) %>%
        hc_yAxis(title=list(style=list(fontSize=axisTitle)),
                 labels=list(style=list(fontSize=axisLabel))) %>%
        hc_plotOptions(series=list(dataLabels=list(
            style=list(fontSize=dataLabel)))) %>%
        hc_legend(itemStyle=list(fontSize=legend)) %>%
        hc_chart(width=700/2.1, height=700/2.1) %>%
        hc_exporting(filename=filename)
    
    # Round percante of variance explained
    roundPercentageInStr <- function(str) {
        perc <- gsub(".*?([0-9]{1,}\\.[0-9]*)%.*", "\\1", str)
        if (str != perc)
            res <- gsub("[0-9]{1,}\\.[0-9]*", round(as.numeric(perc)), str)
        else
            res <- str
        return(res)
    }
    hc$x$hc_opts$xAxis$title$text <- roundPercentageInStr(
        hc$x$hc_opts$xAxis$title$text)
    hc$x$hc_opts$yAxis$title$text <- roundPercentageInStr(
        hc$x$hc_opts$yAxis$title$text)
    save_html(hc, paste0(filename, ".html"))
    return(hc)
}

# fig5af <- browsable(hw_grid(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f,
#                            ncol=2, rowheight="300px"))
# save_html(fig5af, "Fig5a-f.html")
prepareChartForPublication(fig5a, "Fig5a")
prepareChartForPublication(fig5b, "Fig5b")
prepareChartForPublication(fig5c, "Fig5c")
prepareChartForPublication(fig5d, "Fig5d")
prepareChartForPublication(fig5e, "Fig5e")
prepareChartForPublication(fig5f, "Fig5f")
fig5g

# Differential gene expression between isogenic stem cells and fibroblasts -----
diffExpr <- performDiffExpr(
    geneExprNorm$voomQuantileDiscardBadSamples,
    groups$stemCellsVSfib_isogenic)
fig6a <- plotDiffExpr(diffExpr) +
    xlab("log2-Fold change (isoFib - isoSC)") +
    ylab("-log10(FDR)")

# Differential splicing between isogenic stem cells and fibroblasts ------------
diffSplicing <- diffAnalyses(filteredPSI, groups$stemCellsVSfib_isogenic)
diffSplicing <- diffSplicing[order(
    diffSplicing$`Wilcoxon p-value (BH adjusted)`), ]

diffSplicing$label <- ifelse(
    rownames(diffSplicing) ==
        "SE_1_+_207785682_207790253_207790345_207793519_CD46", "CD46", NA)
fig6b <- plotDiffSplicing(diffSplicing, deltaPSIthreshold=0.1,
                          pvalueThreshold=0.01) +
    xlab("∆ Median PSI (isoFib - isoSC)") +
    ylab("-log10(FDR)") +
    geom_label_repel(aes(label=label))

fig6ab <- plot_grid(fig6a, fig6b, labels=c("a", "b"))
fig6ab

# fig6
head(diffSplicing)
saveToTIFF(fig6ab, "Fig6a,b.tiff", width=10)

### What is up with those two strata? Check wilcox_ties.R ###
topDiffSplicingEvents <- head(rownames(diffSplicing), 100)

# Alternative splicing of CD46 is mediated by ESRP: exon inclusion introduces
# a premature termination codon (PTC); transcripts may be NMD-targeted
CD46_PTC_event_hg38 <- grep("CD46", topDiffSplicingEvents, value=TRUE)

groups$stemCellsVSfib_isogenic2 <- groups$stemCellsVSfib_isogenic
names(groups$stemCellsVSfib_isogenic2) <- c("Isogenic stem cells",
                                            "Isogenic fibroblasts")
fig7a <- plotDistribution(filteredPSI[CD46_PTC_event_hg38, ],
                          groups$stemCellsVSfib_isogenic2) %>%
    hc_yAxis(title=list(text="Density"))
fig7a

CD46_ESRP_corr <- correlateGEandAS(geneExprNorm$voomQuantileDiscardBadSamples$E,
                                   filteredPSI, paste0("ESRP", 1:2),
                                   CD46_PTC_event_hg38, method="pearson")
fig7b <- plot(CD46_ESRP_corr)[[1]][[1]]
fig7c <- plot(CD46_ESRP_corr)[[1]][[2]]

fig7bc <- plot_grid(fig7b, fig7c, labels=c("b", "c"))
fig7bc

info  <- queryEnsemblByEvent(CD46_PTC_event_hg38, species="human",
                             assembly="hg38")
fig7d <- plotTranscripts(info, event=CD46_PTC_event_hg38)
fig7d

save_html(fig7a, "Fig7a.html")
saveToTIFF(fig7bc, filepath="Fig7b,c.tiff", width=10)
save_html(fig7d, "Fig7d.html")

# Correlate CD46 penultimate exon inclusion with RBP gene expression -----------
rbps <- getGeneList()$`Sebestyen et al. 2016`$`RNA-binding proteins`
corr <- correlateGEandAS(geneExprNorm$voomQuantileDiscardBadSamples$E,
                         filteredPSI, rbps, CD46_PTC_event_hg38)

corrTable <- as.table(corr)
# plot(density(corrTable$`Pearson's product-moment correlation`))
# plot(density(abs(corrTable$`Pearson's product-moment correlation`)))
# boxplot(corrTable$`Pearson's product-moment correlation`)

# Scatterplots for genes which the correlation was more significant
plot(corr[order(corrTable$`p-value (BH adjusted)`)[1:3]])

# Check correlation based on GTEx data -----------------------------------------
CD46_PTC_event_hg19 <- "SE_1_+_207959027_207963598_207963690_207966864_CD46"

# Check GTEx tissues available based on the sample attributes
GTEx_data <- loadGtexData()[[1]]
GTEx_data$`Gene expression normalised` <- prepareGeneExpr(
    GTEx_data$`Gene expression`, method="quantile", performVoom=TRUE)
rownames(GTEx_data$`Gene expression normalised`) <- convertGeneIdentifiers(
    org.Hs.eg.db, rownames(GTEx_data$`Gene expression normalised`))

GTEx_data$psi <- quantifyPSI(
    GTEx_data$`Junction quantification`,
    genes=parseSplicingEvent(CD46_PTC_event_hg19)$gene[[1]])

GTEx_tissues <- createGroupByAttribute("Tissue Type (area of retrieval)",
                                       GTEx_data$`Sample metadata`)
attr(GTEx_tissues, "Colour") <- setNames(rainbow(length(GTEx_tissues)),
                                         names(GTEx_tissues))

GTEx_corr <- correlateGEandAS(GTEx_data$`Gene expression normalised`,
                              GTEx_data$psi, paste0("ESRP", 1:2),
                              CD46_PTC_event_hg19)
fig8a <- plot(GTEx_corr, colourGroups=GTEx_tissues, legend=TRUE)[[1]][[1]] +
    guides(colour=guide_legend(ncol=3)) +
    theme(legend.text=element_text(size=12),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14, face="bold")) +
    xlab("CD46 penultimate exon inclusion (PSI)") +
    ylab("ESRP2 gene expression (log2CPM)")
fig8a

GTExTissueAnalysis <- lapply(names(GTEx_tissues), corrPerGTExTissue,
                             GTEx_data, gene=paste0("ESRP", 1:2),
                             ASevent=CD46_PTC_event_hg19)
names(GTExTissueAnalysis) <- names(GTEx_tissues)
lapply(names(GTExTissueAnalysis)[-6], function(i)
    plot(GTExTissueAnalysis[[i]])[[1]][[1]] + ggtitle(i))

# Display a higher p-value to make the plot prettier
GTExTissueAnalysis$Esophagus[[1]]$ESRP2$cor$p.value <- 4e-300
  
GTExTissueAnalysis_plots_ESRP2 <- lapply(
    names(GTExTissueAnalysis), prepareCorrPlot, GTExTissueAnalysis, "ESRP2",
    singlePlot=FALSE)
names(GTExTissueAnalysis_plots_ESRP2) <- names(GTExTissueAnalysis)

fig8b <- plot_grid(plotlist=Filter(length, GTExTissueAnalysis_plots_ESRP2))
fig8b <- arrangeGrob(
    fig8b,
    left=textGrob("ESRP2 gene expression (log2CPM)",
                  gp=gpar(fontface="bold", fontsize=14), rot=90),
    bottom=textGrob("CD46 penultimate exon inclusion (PSI)",
                    gp=gpar(fontface="bold", fontsize=14)))

fig8 <- plot_grid(fig8a, fig8b, labels=c("a", "b"), nrow=2, rel_heights=c(1, 2),
                  label_size=14)
saveToTIFF(fig8, "Fig8.tiff", width=10, height=12)
fig8

# Check correlation based on TCGA data -----------------------------------------
library(survival)
panTCGAanalysis <- lapply(names(getFirebrowseCohorts()), function(cohort) {
    try(performCorrAndSurvPerTCGACohort(cohort, CD46_PTC_event_hg19,
                                        c(paste0("ESRP", 1:2), "CD46")))
})
names(panTCGAanalysis) <- names(getFirebrowseCohorts())
panTCGAanalysis_clean <- Filter(function(x) !is(x, "try-error"),
                                panTCGAanalysis)

panTCGAanalysis_corr <- lapply(panTCGAanalysis_clean, "[[", "corr")

# Ignore multi-cohorts
multiCohorts <- c("COADREAD", "GBMLGG", "KIPAN", "STES")
for (cohort in multiCohorts) panTCGAanalysis_corr[[cohort]] <- NULL

panTCGAanalysis_corr_plots <- list()
panTCGAanalysis_corr_plots$ESRP1 <- lapply(
    names(panTCGAanalysis_corr), prepareCorrPlot, panTCGAanalysis_corr,
    "ESRP1|54845")
panTCGAanalysis_corr_plots$ESRP2 <- lapply(
    names(panTCGAanalysis_corr), prepareCorrPlot, panTCGAanalysis_corr,
    "ESRP2|80004")

fig9 <- plot_grid(plotlist=Filter(length, panTCGAanalysis_corr_plots$ESRP1))
fig9 <- arrangeGrob(
    fig9,
    left=textGrob("ESRP1 gene expression (log2CPM)",
                  gp=gpar(fontface="bold", fontsize=14), rot=90),
    bottom=textGrob("CD46 penultimate exon inclusion (PSI)",
                    gp=gpar(fontface="bold", fontsize=14)))
ggsave("Fig9.tiff", dpi=600, height=10, width=12, fig9, compression="lzw")
plot(fig9)

fig10 <- plot_grid(plotlist=Filter(length, panTCGAanalysis_corr_plots$ESRP2))
fig10 <- arrangeGrob(
    fig10,
    left=textGrob("ESRP2 gene expression (log2CPM)",
                  gp=gpar(fontface="bold", fontsize=14), rot=90),
    bottom=textGrob("CD46 penultimate exon inclusion (PSI)",
                    gp=gpar(fontface="bold", fontsize=14)))
ggsave("Fig10.tiff", dpi=600, height=10, width=12, fig10, compression="lzw")
fig10 <- plot(fig10)
fig10

# Check survival on TCGA data --------------------------------------------------
plotSurvs <- function(item, panTCGAanalysis_clean) {
    surv <- panTCGAanalysis_clean[[item]]$surv
    pval <- panTCGAanalysis_clean[[item]]$pvalue
    plotSurvivalCurves(surv, pvalue=pval, mark=FALSE, title=item)
}
lapply(names(panTCGAanalysis_clean), plotSurvs, panTCGAanalysis_clean)
fig11a <- plotSurvs("LGG", panTCGAanalysis_clean)
fig11b <- plotSurvs("LUAD", panTCGAanalysis_clean)

# Plot p-values by cutoff
tcgaSurv <- loadFirebrowseData(folder=getDownloadsFolder(),
                               cohort=c("LGG", "LUAD"),
                               data=c("clinical", "junction_quantification"),
                               date="2016-01-28")

survPvaluesByCutoff <- function(data) {
    clinical <- data[["Clinical data"]]
    junctionQuant <- data[["Junction quantification (Illumina HiSeq)"]]
    sampleInfo <- data[["Sample metadata"]]
    
    message("Quantifying alternative splicing...")
    psi <- quantifyPSI(junctionQuant, genes=parseSplicingEvent(
        CD46_PTC_event_hg19)$gene[[1]])
    eventPSI <- psi[CD46_PTC_event_hg19, ]
    
    message("Plotting survival p-value plots by cutoff...")
    match <- getSubjectFromSample(colnames(psi), clinical, sampleInfo=sampleInfo)
    
    types  <- createGroupByAttribute("Sample types", sampleInfo)
    tumour <- types$`Primary solid Tumor`
    
    subjectPSI <- assignValuePerPatient(eventPSI, match, clinical,
                                        samples=unlist(tumour))
    plotSurvivalPvaluesByCutoff(clinical, subjectPSI, censoring="right",
                                event="days_to_death",
                                timeStart="days_to_death")
}
fig11c <- survPvaluesByCutoff(tcgaSurv$`Brain Lower Grade Glioma 2016-01-28`)
fig11d <- survPvaluesByCutoff(tcgaSurv$`Lung adenocarcinoma 2016-01-28`)

fig11 <- browsable(hw_grid(fig11a, fig11b, fig11c, fig11d, ncol=2))
fig11
save_html(hw_grid(fig11a, fig11b, fig11c, fig11d, ncol=2), "Fig11.html")
