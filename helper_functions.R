# Stem cell analysis using recount dataset SRP063867 (functions)
# Nuno Agostinho, 19 Feb 2019

library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
library(htmltools)
library(highcharter)
library(scales)

# Groups -----------------------------------------------------------------------
prepareGroups <- function(sampleInfo) {
    groups <- list()
    groups$cellType <- createGroupByAttribute("cell type", sampleInfo)
    groups$dataSet  <- createGroupByAttribute("data set",  sampleInfo)
    groups$stemCellsVSfib$SC <- unlist(groups$cellType[c("iPSC", "ESC")])
    groups$stemCellsVSfib$Fib <- groups$cellType$Fibroblast
    groups$stemCellsVSfib_isogenic <- groups$stemCellsVSfib
    for (group in names(groups$stemCellsVSfib_isogenic)) {
        groups$stemCellsVSfib_isogenic[[group]] <- intersect(
            groups$stemCellsVSfib_isogenic[[group]], groups$dataSet$isogenic)
    }
    
    title   <- createGroupByAttribute("title", sampleInfo)
    groups$iPSCfib_vs_ESCfib$iPSCfib <- intersect(
        groups$cellType$Fibroblast, unlist(title[grepl("iPSC", names(title))]))
    groups$iPSCfib_vs_ESCfib$ESCfib  <- intersect(
        groups$cellType$Fibroblast, unlist(title[grepl("ESC", names(title))]))
    
    return(groups)
}

# Gene expression filtering, normalisation and analysis ------------------------
# Convert gene names from ENSEMBL to gene symbol (if available)
library(data.table)
library(edgeR)
library(org.Hs.eg.db)

prepareGeneExpr <- function(geneExpr, minCounts=10, minTotalCounts=15, ...) {
    keep <- filterGeneExpr(geneExpr, minCounts=minCounts,
                           minTotalCounts=minTotalCounts)
    geneExprNorm <- normaliseGeneExpression(geneExpr, keep, ...)
    rownames(geneExprNorm) <- convertGeneIdentifiers(org.Hs.eg.db,
                                                     rownames(geneExprNorm))
    return(geneExprNorm)
}

# Differential expression analysis ---------------------------------------------
library(limma)
performDiffExpr <- function(data, groups, coef=2) {
    # Prepare groups of samples to analyse and further filter unavailable
    # samples in selected groups for gene expression
    ge <- data
    groups <- psichomics:::discardOutsideSamplesFromGroups(groups, colnames(ge))
    
    ge           <- ge[ , unlist(groups)]
    isFromGroup1 <- colnames(ge) %in% groups[[1]]
    design       <- cbind(1, ifelse(isFromGroup1, 0, 1))
    
    # Fit a gene-wise linear model based on selected groups
    fit <- lmFit(ge, design)
    
    # Calculate moderated t-statistics and DE log-odds using limma::eBayes
    ebayesFit <- eBayes(fit, trend=!is(data, "EList"))
    
    # Prepare data summary
    pvalueAdjust <- "BH" # Benjamini-Hochberg p-value adjustment (FDR)
    summary <- topTable(ebayesFit, number=nrow(fit), coef=coef, sort.by="none",
                        adjust.method=pvalueAdjust, confint=TRUE)
    summary$ID <- NULL
    names(summary) <- c("log2 Fold-Change", "CI (low)", "CI (high)",
                        "Average expression", "moderated t-statistics",
                        "p-value",
                        paste0("p-value (", pvalueAdjust, " adjusted)"),
                        "B-statistics")
    attr(summary, "groups") <- groups
    
    # Calculate basic statistics
    stats <- diffAnalyses(ge, groups, "basicStats", pvalueAdjust=NULL)
    final <- cbind(stats, summary)
    
    final$minusLog10qvalue <- -log10(final$`p-value (BH adjusted)`)
    final <- final[order(final$`p-value (BH adjusted)`), ]
    final <- final[ , -c(1:3)]
    return(final)
}

library(ggrepel)
plotDiffExpr <- function(data) {
    # cognateGenes <- unlist(parseSplicingEvent(events)$gene)
    logFCthreshold  <- abs(data$`log2 Fold-Change`) > 1
    pvalueThreshold <- data$`p-value (BH adjusted)` < 0.01
    
    data$genes <- gsub("\\|.*$", "\\1", rownames(data))
    
    ggplot(data, aes(`log2 Fold-Change`,
                     -log10(`p-value (BH adjusted)`))) +
        geom_point(data=data[logFCthreshold & pvalueThreshold, ],
                   colour="orange", alpha=0.5, size=3) +
        geom_point(data=data[!logFCthreshold | !pvalueThreshold, ],
                   colour="gray", alpha=0.5, size=3) +
        # geom_text_repel(data=data[cognateGenes, ], aes(label=genes),
        # box.padding=0.4, size=5) +
        theme_light(16) +
        ylab("-log10(q-value)")
}

# PSI filtering ----------------------------------------------------------------
# Filter splicing quantification based on variance
samplesPerRow <- function(df) rowSums(!is.na(df))

# Differential splicing --------------------------------------------------------
library(ggplot2)
plotDiffSplicing <- function(diffSplicing,
                             pvalCol="Wilcoxon p-value (BH adjusted)",
                             deltaPSIthreshold=0.1, pvalueThreshold=0.01) {
    deltaPSIpoints <- abs(diffSplicing$`∆ Median`) > deltaPSIthreshold
    pvaluePoints   <- diffSplicing[[pvalCol]] <= pvalueThreshold
    
    diffSplicing$sampleNumber <- rowSums(
        diffSplicing[ , grep("Samples", colnames(diffSplicing))[[1]],
                      drop=FALSE])
    
    subset <- nrow(diffSplicing[deltaPSIpoints & pvaluePoints, ])
    total  <- nrow(diffSplicing)
    perc   <- subset/total * 100
    message(sprintf("%s (%s%%) out of %s splicing events highlighted",
                    subset, round(perc), total))
    
    highlighted <- diffSplicing[which(deltaPSIpoints & pvaluePoints), ]
    notHighlighted <- diffSplicing[which(!deltaPSIpoints | !pvaluePoints), ]
    
    ggplot(diffSplicing, aes(diffSplicing$`∆ Median`,
                             -log10(diffSplicing[[pvalCol]]))) +
        geom_point(
            data=notHighlighted,
            mapping=aes(notHighlighted$`∆ Median`,
                        -log10(notHighlighted[[pvalCol]])),
            colour="grey", alpha=0.5, size=2, na.rm=TRUE) +
        geom_point(
            data=highlighted,
            mapping=aes(highlighted$`∆ Median`, -log10(highlighted[[pvalCol]])),
            alpha=0.5, size=2, na.rm=TRUE, colour="orange") +
        # geom_vline(xintercept=-deltaPSIthreshold) +
        # geom_vline(xintercept=deltaPSIthreshold) +
        # geom_hline(yintercept=-log10(pvalueThreshold)) +
        # guides(colour=guide_legend(title="Samples")) +
        theme_light(16) +
        ylab("-log10(FDR)")
}

# GTEx and TCGA ----------------------------------------------------------------
quantifyPSI <- function(junctionQuant, genes=NULL) {
    human <- listSplicingAnnotations()[[1]]
    annotation <- loadAnnotation(human)
    psi <- quantifySplicing(annotation, junctionQuant, genes=genes)
}

performCorrAndSurvPerTCGACohort <- function(cohort, ASevent, gene) {
    data <- loadFirebrowseData(
        folder=getDownloadsFolder(), date="2016-01-28", cohort=cohort,
        data=c("clinical", "junction_quantification", "RSEM_genes"))[[1]]
    
    geneExpr <- data$`Gene expression`
    if (is.null(geneExpr))
        geneExpr <- data$`Gene expression (Illumina HiSeq)`
    data[["Gene expression normalised"]] <- prepareGeneExpr(geneExpr,
                                                            performVoom=TRUE,
                                                            method="quantile")
    
    junctionQuant <- data$`Gene expression`
    if (is.null(junctionQuant))
        junctionQuant <- data$`Junction quantification (Illumina HiSeq)`
    data$psi <- quantifyPSI(junctionQuant,
                            parseSplicingEvent(ASevent)$gene[[1]])
    
    corr <- correlateGEandAS(data[["Gene expression normalised"]]$E, data$psi,
                             gene, ASevent)
    
    # Check survival differences for AS event based on tumour samples
    # Find optimal cutoff for the event
    sampleTypes <- createGroupByAttribute(
        "Sample types", data$`Sample metadata`)
    if ("Primary solid Tumor" %in% names(sampleTypes)) {
        tumour <- unlist(sampleTypes[["Primary solid Tumor"]])
    } else {
        tumour <- unlist(sampleTypes[[grep("Tumor", names(sampleTypes))]])
    }
    match <- getSubjectFromSample(colnames(data$psi),
                                  data$`Clinical data`,
                                  sampleInfo=data$`Sample metadata`)
    eventPSI <- assignValuePerPatient(data$psi[ASevent, ], match,
                                      data$`Clinical data`, samples=tumour)
    opt <- optimalSurvivalCutoff(data$`Clinical data`, eventPSI,
                                 censoring="right", event="days_to_death",
                                 timeStart="days_to_death")
    (optimalCutoff <- opt$par)    # Optimal exon inclusion level
    (optimalPvalue <- opt$value)  # Respective p-value
    
    label     <- labelBasedOnCutoff(eventPSI, round(optimalCutoff, 2),
                                    label="PSI values")
    survTerms <- processSurvTerms(data$`Clinical data`, censoring="right",
                                  event="days_to_death",
                                  timeStart="days_to_death",
                                  group=label, scale="years")
    surv <- survfit(survTerms)
    pvalue <- testSurvival(survTerms)
    plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
    
    res <- list(corr=corr, surv=surv, pvalue=pvalue)
    return(res)
}

prepareCorrPlot <- function(type, corr, gene, singlePlot=TRUE) {
    if (is.null(corr[[type]][[1]][[gene]]$cor)) return(NULL)
    p <- plot(corr[[type]])[[1]][[gene]]
    r <- round(corr[[type]][[1]]$ESRP2$cor$estimate, 1)
    pval <- signif(corr[[type]][[1]]$ESRP2$cor$p.value, 1)
    if (!is.null(p)) {
        p <- p + ggtitle(type, sprintf("r: %s (P: %s)", r, pval)) +
            xlab(NULL) + ylab(NULL) +
            scale_x_continuous(breaks=pretty_breaks(3), expand=c(0, 0),
                               labels=c("0"="0", "0.5"="0.5", "1"="1")) +
            scale_y_continuous(breaks=pretty_breaks(3, min.n=3))
        
        if (type %in% c("Spleen", "Nerve")) {
            p <- p + suppressWarnings(
                coord_cartesian(expand=FALSE, xlim=c(0, 1)))
        } else if (type == "TGCT" && gene == "ESRP1|54845") {
            p <- p + suppressWarnings(
                coord_cartesian(expand=FALSE, xlim=c(0, 1), ylim=c(0, 10)))
        }
        
        if (singlePlot) {
            p <- p + theme(axis.text=element_text(size=14),
                           plot.title=element_text(size=14),
                           plot.subtitle=element_text(size=14))
        } else {
            p <- p + theme(axis.text=element_text(size=14),
                           plot.title=element_text(size=13),
                           plot.subtitle=element_text(size=12))
        }
    }
    return(p)
}

corrPerGTExTissue <- function(tissue, GTEx_data, gene, ASevent) {
    GTEx_tissues <- createGroupByAttribute("Tissue Type (area of retrieval)",
                                           GTEx_data$`Sample metadata`)
    samples <- GTEx_tissues[[tissue]]
    
    geneExpr <- GTEx_data[["Gene expression normalised"]]$E
    geneExpr <- geneExpr[ , colnames(geneExpr) %in% samples]
    
    psi2 <- GTEx_data$psi
    psi2 <- psi2[ , colnames(psi2) %in% samples]
    
    corr <- correlateGEandAS(geneExpr, psi2, gene=gene, ASevents=ASevent)
    return(corr)
}

# Save figures -----------------------------------------------------------------

saveToPDF <- function(object, filepath, width=12/2.2, height=8/2.2) {
    cairo_pdf(filepath, width=width, height=height)
    print(object)
    dev.off()
}

saveToTIFF <- function(object, filepath, width=5, height=5, res=600) {
    tiff(filepath, units="in", width=width, height=height, res=res,
         compression="lzw")
    print(object)
    dev.off()
}
