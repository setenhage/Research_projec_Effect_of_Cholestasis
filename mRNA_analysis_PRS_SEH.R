#==============================================================================
# Analysis of Microarray data of liver tissue samples
# mRNA
#
# Version 1.0
# Course: PRO4002 Research project 1 MSB
# Names: S.E. ten Hage, P.R. Stolpe, MaCSBio,
#        Maastricht University
# ID: i6201596, i6195641
#
#
# History: 
#    1.0: Creation
#==============================================================================
# 0. load required packages
#==============================================================================
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "biomaRt",
                     "clusterProfiler", "rJava",
                     "RDAVIDWebService", "org.Hs.eg.db"))
library(gplots)
library(limma)
library(pcaMethods)
library(biomaRt)
library(dplyr)
library(xlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(rJava)
library(RDAVIDWebService)
options(stringsAsFactors = F)
#
#
#==========================================================
# 1. import data files
#==========================================================
# define directories according to where the files are in and
# results are supposed to be found after analysis
DATA.DIR <- "~/Documents/Master/Maastricht/Research project 1/data"
RESULTS.DIR <- "~/Documents/Master/Maastricht/Research project 1/results/results mRNA"
# define working directory
setwd(DATA.DIR)

# 1.1 import the normalized gene expression data into an object: Hepa.norm
#------------------------------------------------------------------------------
Hepa.norm <- read.delim("GeneExpressionNormalized.txt",
                        row.names = 1)

# 1.2 reshuffle the data according to experimental design
#------------------------------------------------------------------------------
# create three vectors with three groups of sampples: drained, cholestatic and control
drained.norm <- cbind("FGS_02_410978_1_2", "FGS_03_410978_1_3", 
                      "FGS_08_410978_2_2", "FGS_11_410979_1_1", "FGS_14_410979_1_4", 
                      "FGS_17_410979_2_3", "FGS_20_410980_1_1", "FGS_23_410980_1_4",
                      "FGS_26_410980_2_3", "FGS_30_412287_1_3")

cholestatic.norm <- cbind("FGS_01_410978_1_1", "FGS_04_410978_1_4", 
                          "FGS_10_422569_2_3_H", "FGS_13_410979_1_3", "FGS_16_410979_2_2",
                          "FGS_19_410979_2_4", "FGS_22_410980_1_3", "FGS_25_410980_2_2",
                          "FGS_28_412287_1_1")

control.norm <- cbind("FGS_06_410978_2_1", "FGS_09_410978_2_3", 
                      "FGS_12_410979_1_2", "FGS_15_410979_2_1", "FGS_21_410980_1_2",
                      "FGS_24_410980_2_1", "FGS_27_410980_2_4", "FGS_29_412287_1_2",
                      "FGS_32_412287_1_4")
# combine reordered sample with Genesymbol column: Hepa.norm
Hepa.norm <- cbind(Genesymbol = Hepa.norm$Genesymbol, Hepa.norm[, cholestatic.norm], 
                        Hepa.norm[, drained.norm], Hepa.norm[, control.norm])
# use dimension check as check for correct combining
dim(Hepa.norm)

#==========================================================
# 2. Inspect data object: Hepa.norm
#==========================================================
class(Hepa.norm)
str(Hepa.norm)
max(Hepa.norm[, 2:29])
min(Hepa.norm[, 2:29])
mean(as.matrix(Hepa.norm[, 2:29]))

#==============================================================================
# 3. Visualization of diagnostic plots
#==============================================================================

# 3.1 Boxplot of normalized data 
#------------------------------------------------------------------------------
setwd(RESULTS.DIR)
pdf("visualization of normalized data mRNA.pdf", 7, 5)
boxplot(Hepa.norm[, 2:29], col = c(rep("lightblue", 9), rep("orange", 10), 
                              rep("darkgreen", 9)), 
        ylab = "gene expression levels",
        outline = F)
legend("topright", legend = c("cholestatic", "drained", "control"),
       fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)
title("Boxplot of normalized gene expression data")

# 3.2 PCA of normalized data and remove batch effects
#------------------------------------------------------------------------------
#make diagnostic plot before the removal of batch effects 
plotMDS(Hepa.norm[, 2:29], 
        pch = c(rep(0, 9), rep(1, 10), rep(2, 9)),
        col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9))
        , gene.selection = "common",
        main = "PCA before removal of batch effects")
legend("topright",
       legend = c("cholestatic", "drained", "control"),
       pch = c(0,1,2), cex = 0.75,
       col = c("lightblue", "orange", "darkgreen"))
# match chip IDs to samples: A,B,C,D -> chips 1-4
batch <- c("A", "A", "A", "B", "B", "B", "C", "C", "D",
           "A", "A", "A", "B", "B", "B", "C", "C", "C", "D",
           "A", "A", "B", "B", "C", "C", "C", "D", "D")
# produce design matrix according to experimental design
design <- model.matrix(~0+factor(c(1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3)))
# remove batch effects
Hepa.norm[, 2:29] <- removeBatchEffect(Hepa.norm[, 2:29], batch = batch, 
                                       design = design)
# make diagnostic plot after removal of batch effects for confirmation
# also to see if batch effects were true 
plotMDS(Hepa.norm[, 2:29], 
        pch = c(rep(0, 9), rep(1, 10), rep(2, 9)),
        col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9))
        , gene.selection = "common",
        main = "PCA after removal of batch effects")
legend("topleft",
       legend = c("cholestatic", "drained", "control"),
       pch = c(0,1,2), cex = 0.75,
       col = c("lightblue", "orange", "darkgreen"))

# 3.3 produce demonstration of more PCs using pca function
#------------------------------------------------------------------------------
# # replace NaN with NA to be readable for PCA function
 Hepa.norm[is.na(Hepa.norm)] <- NA
# # again check for correct dimensions
# dim(Hepa.norm.rBE)
# # pca of Hepa.norm.shuf using ppca: pcaRes
pcaRes <- pca(t(Hepa.norm[, 2:29]), method = "ppca", center = TRUE,  nPcs = 10)
# # plot Pca Results 
 plot(pcaRes, main = "PCA of normalized data")
# # make plots of principal components agains each other
# plotPcs(pcaRes, c(1,2), type = "scores", sl = NULL, 
#         col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
# legend("topright", legend = c("cholestatic", "drained", "control"),
#        fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)
# 
plotPcs(pcaRes, c(1,3), sl = NULL, 
         col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
# legend("topright", legend = c("cholestatic", "drained", "control"),
#        fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)
# 
plotPcs(pcaRes, c(2,3), sl = NULL, 
         col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
# legend("topright", legend = c("cholestatic", "drained", "control"),
#        fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)

dev.off()



#==============================================================================
# 4. statistical analysis
#==============================================================================
# Design matrix of the experimental design
design <- model.matrix(~0+factor(c(1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3)))
colnames(design) <- c("cholestatic", 
                             "drained", 
                             "control")
# use lmFit to create linear regression model of gene expression according to
# design matrix
Hepa.norm.fit <- lmFit(Hepa.norm[, 2:29], design = design)

# build contrast matrix of the desired comparisons
cm <- makeContrasts(cholestatic-control, 
                    drained-control, 
                    cholestatic-drained,
                    levels = design)
# make contrast fit to submit for further analysis of significant
# differences in gene expression levels
Hepa.norm.fit2 <- contrasts.fit(Hepa.norm.fit, cm)
# compute t-statistics, F-statistics and log-odds of linear model fit
Hepa.norm.fit2.eBayes <- eBayes(Hepa.norm.fit2)
# create object with decisions based on default options BH, p<0.05 and lfc = 1.5
Hepa.norm.fit2.decide <- decideTests(Hepa.norm.fit2.eBayes, lfc = 1.5)
# adjust for multiple testing using topTable: BH, p<0.05, lfc = 1.5
statscc <- topTable(Hepa.norm.fit2.eBayes,
                  number = nrow(Hepa.norm.fit2.eBayes),
                  adjust.method = "BH", 
                  coef = 1, p.value = 0.05, lfc = 1.5,
                  sort.by = "none")
# controlvdrained: 1 significant gene: no further analysis
statscd <- topTable(Hepa.norm.fit2.eBayes,
                    number = nrow(Hepa.norm.fit2.eBayes),
                    adjust.method = "BH", 
                    coef = 2, p.value = 0.05, lfc = 1.5,
                    sort.by = "none")
# drainedvcholestatic
statsdc <- topTable(Hepa.norm.fit2.eBayes,
                    number = nrow(Hepa.norm.fit2.eBayes),
                    adjust.method = "BH", 
                    coef = 3, p.value = 0.05, lfc = 1.5,
                    sort.by = "none")


#dim(stats)
# produce histogram of p values for controlvcholestatic (unadjusted)
pdf("Histgram Volcanoplot mRNA.pdf", 7, 5)
# hist(statscc[, "P.Value"], xlab = "P.Value", 
#      main = "Histogram of P-Values")
# produce volcanoplot and highlight 5 most significant genes with names
volcanoplot(Hepa.norm.fit2.eBayes, highlight = 5,
            coef = 1,
            names = Hepa.norm$Genesymbol,
            xlab = "Log2 FC", 
            ylab = "-Log10(pValue)",
            main = "Volcanoplot of P.Values cholestatic v control")
volcanoplot(Hepa.norm.fit2.eBayes, highlight = 5,
            coef = 3,
            names = Hepa.norm$Genesymbol,
            xlab = "Log2 FC", 
            ylab = "-Log10(pValue)",
            main = "Volcanoplot of P.Values cholestatic v drained")
dev.off()

# limma::plotMA(Hepa.norm.fit2.eBayes,
#               coef = 1,
#        xlab = "Average log-expression", 
#        ylab = "sqrt(sigma)",
#        zero.weights = FALSE,
#        pch = 16, 
#        cex = 0.3)
# limma::plotMA(Hepa.norm.fit2.eBayes,
#               coef = 3,
#               xlab = "Average log-expression", 
#               ylab = "sqrt(sigma)",
#               zero.weights = FALSE,
#               pch = 16, 
#               cex = 0.3)
# limma::plotMA(Hepa.norm.fit2.eBayes,
#               coef = 2,
#               xlab = "Average log-expression", 
#               ylab = "sqrt(sigma)",
#               zero.weights = FALSE,
#               pch = 16, 
#               cex = 0.3)
# creating venn diagrams of statistical analysis
pdf("Venn Diagram Heatmaps mRNA.pdf", 10, 9)   
# upregulated genes
vennDiagram(Hepa.norm.fit2.decide, include = c("up","down"), 
            counts.col = c("red", "green"),
            show.include = TRUE)
title("Venn Diagram of up-/downregulated genes after statistical analysis")
# # produce heatmap of differential expressed genes from cholestatic v control
# heatmap(as.matrix(Hepa.norm[rownames(statscc), 2:29]), distfun = dist,
#         ColSideColors = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)),
#         labRow = Hepa.norm$Genesymbol)
# title("Heatmap of differentially expressed genes in cholestatic v control")
# # produce heatmap of differential expressed genes from cholestatic v drained
# heatmap(as.matrix(Hepa.norm[rownames(statsdc), 2:29]), distfun = dist,
#         ColSideColors = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)),
#         labRow = Hepa.norm$Genesymbol)
# title("Heatmap of differentially expressed genes in cholestatic v drained")
dev.off()


#==============================================================================
# 5. Adding annotation to the normalized data obejct:
# Hepa.norm (gene symbols, names), statscc, statsdc
#==============================================================================
# connect to ensembl database
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
# retrieve gene ids, and description from ensembl database
ensembl.query <- getBM(attributes = c("ensembl_gene_id", "description", "hgnc_symbol",
                                      "entrezgene"), 
                       filters = "hgnc_symbol", 
                       values = Hepa.norm$Genesymbol, 
                       mart = ensembl)

#match gene annotation with normalized data 
ensembl.match <- ensembl.query[match(table = ensembl.query$hgnc_symbol, 
                                             x = Hepa.norm$Genesymbol), ]
# attach annotation: Hepa.norm
Hepa.norm <- cbind(Hepa.norm, ensembl.match)

#match gene annotation with significant gene lists
statscc.match <- ensembl.query[match(table = ensembl.query$entrezgene,
                                     x = rownames(statscc)), ]
statsdc.match <- ensembl.query[match(table = ensembl.query$entrezgene,
                                     x = rownames(statsdc)), ]
statscc <- cbind(statscc, statscc.match)
statsdc <- cbind(statsdc, statsdc.match)
setwd(RESULTS.DIR)
# export significantly changed genes to excel file: SignificantGenes
xlsx::write.xlsx(statscc, "SignificantGenes.xlsx", sheetName = "controlvcholestatic",
           row.names = TRUE, col.names = TRUE)
xlsx::write.xlsx(statsdc, "SignificantGenes.xlsx", sheetName = "drainedvcholestatic",
           row.names = TRUE, col.names = TRUE, append = TRUE)

#==============================================================================
# 6. Pathway and Gene ontology enrichment
#==============================================================================
# create vector with entrez gene ids
gene.symb <- Hepa.norm$entrezgene
# perform GO enrichment analysis on controlvcholestatic
# and drainedvcholestatic using goana
go.enrichcc <- goana(Hepa.norm.fit2.eBayes, geneid = gene.symb, coef = 1,
                     species = "Hs")
go.enrichdc <- goana(Hepa.norm.fit2.eBayes, geneid = gene.symb, coef = 3,
                   species = "Hs")
# print results to objects
GOcc <- topGO(go.enrichcc, ontology = "BP", number = 20)
GOdc <- topGO(go.enrichdc, ontology = "BP", number = 20)
# create xlsx file with two sheets
xlsx::write.xlsx(GOcc, "EnrichmentAnalysis.xlsx", sheetName = "GO_controlvcholestatic")
xlsx::write.xlsx(GOdc, "EnrichmentAnalysis.xlsx", sheetName = "GO_drainedvcholestatic", 
           append = TRUE)

#------------------------------------------------------------------------------
# using DAVID GO enrichment function (adjusted for MT)
ccdavid <- enrichDAVID(rownames(statscc), idType = "ENTREZ_GENE_ID" , 
                       pvalueCutoff = 0.05, minGSSize = 10, pAdjustMethod = "BH", 
                       david.user = "r.stolpe@student.maastrichtuniversity.nl")
# enrichment is negative with multiple testing: no proceedings
dcdavid <- enrichDAVID(rownames(statsdc), idType = "ENTREZ_GENE_ID" , pvalueCutoff = 0.05,
            pAdjustMethod = "BH", david.user = "r.stolpe@student.maastrichtuniversity.nl")
# dcdavid has been found empty and will not be printed to excel file
xlsx::write.xlsx(ccdavid, "EnrichmentAnalysis.xlsx", sheetName = "enrichDAVID_controlvcholestatic",
           col.names = TRUE, append = TRUE)

#------------------------------------------------------------------------------
# perform KEGG database enerichment analysis on controlvcholestatic
# and drainedvcholestatic <- no adjustment for multiple testing
kegg.enrichdc <- kegga(Hepa.norm.fit2.eBayes, geneid = gene.symb, 
                                      coef = 3, species = "Hs")
kegg.enrichcc <- kegga(Hepa.norm.fit2.eBayes, geneid = gene.symb, 
                                      coef = 1, species = "Hs")
# print results to objects
KEGGdc <- topKEGG(kegg.enrichdc, number = 20)
KEGGcc <- topKEGG(kegg.enrichcc, number = 20)
# create xlsx file with two sheets
xlsx::write.xlsx(KEGGcc, "EnrichmentAnalysis.xlsx", sheetName = "kegga_controlvcholestatic", 
           append = TRUE)
xlsx::write.xlsx(KEGGdc, "EnrichmentAnalysis.xlsx", sheetName = "kegga_drainedvcholestatic", 
           append = TRUE)

# perform KEGG pathway enrichment analysis: cc.kegg (multiple testing adjusted)
# this will produce enriched kegg pathways to visualize using online 
# pathway viewer
cc.kegg <- enrichKEGG(rownames(statscc), organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
dc.kegg <- enrichKEGG(rownames(statsdc), organism = "hsa", keyType = "kegg", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
# results are now looked at. Genes will be fed into online pathway viewer and
# visualized according to enrichment. See also report

# # import of file which contains confirmed effects of miRNA on mRNA: 21 genes 
# miRNAcheckcc <- read.delim("miRNA mRNA gene pairs. ConChol", header = TRUE)
# # performing enrichment analysis on effected genes
statscc.go <- enrichGO(rownames(statscc), OrgDb = org.Hs.eg.db, ont = "BP", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH")
statsdc.go <- enrichGO(rownames(statsdc), OrgDb = org.Hs.eg.db, ont = "BP",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")
# only control v cholestatic produced results 
setwd(RESULTS.DIR)
pdf("Gene Ontology plots.pdf", 14,11)
dotplot(statscc.go, title = "Gene Ontology Analysis using enrichGO")
dotplot(ccdavid, title = "Gene Ontology Analysis using enrichDAVID")
dev.off()
#==============================================================================
# # import of phenotype data
# setwd(DATA.DIR)
# pheno.data <- read.xlsx("Copy of metadata tbv arrays.xlsx")



# #==============================================================================
# # 7. WGCNA analysis controlvcholestatic
# #==============================================================================
# NETWORKS.DIR <- "~/Documents/Master/Maastricht/Research project 1/results/networks"
# setwd(NETWORKS.DIR)
# # check for the quality of samples and genes -> if last command returns 
# # true the dataset is complete and doesnt contain too many missing values
# gsg <- goodSamplesGenes(Hepa.norm[, 2:29], verbose = 3)
# gsg$allOK
# #create subset of significantly changed genes: signcc
# signcc <- Hepa.norm[rownames(statscc), ]
# # check for outliers and also look at the distance of samples to each other
# samplesTreecc <- hclust(dist(t(signcc[, 2:29])), method = "average")
# sizeGrWindow(12,9)
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# # plot clustering of samples 
# plot(samplesTreecc, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# 
# # Choose a set of soft-thresholding powers: powers
# powers = c(1:20)
# # Call the network topology analysis function
# sft = pickSoftThreshold(t(signcc[, 2:29]), powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90, col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# 
# #wgcna
# # pick power that has scale independence above threshold and perform wgcna
# # plot results of different modules that have been created 
# ## here is where optimization needs to happen because no good scale independence could be 
# ## achieved yet: better statistical analysis meaning remove batch effects and different 
# ## testing method. 
# powercc <- 13
# wgcnaRescc <- blockwiseModules(t(signcc[, 2:29]), power = powercc, TOMType = "unsigned", 
#                              minModuleSize = 10, verbose = 3)
# mergedColorscc <- wgcnaRescc$colors
# plotDendroAndColors(wgcnaRescc$dendrograms[[1]],
#                     mergedColorscc[wgcnaRescc$blockGenes[[1]]],
#                     "Module colors", dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# # Create list with all genes from the different modules (for easy lookup)
# moduleGenescc <- split(signcc[, 1], wgcnaRescc$colors)
# 
# # Define modules of interest. 
# moduleSelectioncc <- c("blue", "turquoise", "brown", "yellow")
# 
# 
# #-----------------------------------------------------------------------------#
# # Export data for Cytoscape
# #-----------------------------------------------------------------------------#
# # Export the original input data with module membership
# # pick threshold so to cut off modules connected at a high position .. the higher
# # the threshold the better for the network formation and interpretation 
# cytoScapeAnnotation <- cbind(signcc[, 2:29], "WGCNA-Module" = wgcnaRescc$colors)
# write.table(cytoScapeAnnotation, "DataForCytoscapecc.txt", sep = "\t", quote = F, row.names = F)
# 
# TOM <- TOMsimilarityFromExpr(t(signcc[, 2:29]), power = powercc)
# probes <- signcc$Genesymbol
# moduleColors <- labels2colors(wgcnaRescc$colors)
# moduleSelection <- c("blue", "turquoise", "brown", "yellow")
# inModule = is.finite(match(moduleColors, moduleSelection))
# modProbes <- probes[inModule]
# modTOM <- TOM[inModule, inModule]
# dimnames(modTOM) <- list(modProbes, modProbes)
# cyt <- exportNetworkToCytoscape(modTOM,
#                                 edgeFile=paste("CytoscapeInput-edges-", paste(moduleSelection, collapse="-"), ".txt", sep=""),
#                                 nodeFile=paste("CytoscapeInput-nodes-", paste(moduleSelection, collapse="-"), ".txt", sep=""),
#                                 weighted=TRUE,
#                                 threshold=0.2,
#                                 nodeNames=modProbes,
#                                 altNodeNames=modProbes,
#                                 nodeAttr=moduleColors[inModule])
# 
# 
# 
# ##### add more comments to repitions of wgcna analysis and export 
# 
# 
# 
# #==============================================================================
# # 8. WGCNA analysis drainedvcholestatic
# #==============================================================================
# # check for the quality of samples and genes -> if last command returns 
# # true the dataset is complete and doesnt contain too many missing values
# gsg <- goodSamplesGenes(Hepa.norm[, 2:29], verbose = 3)
# gsg$allOK
# 
# signdc <- Hepa.norm[row.names(statsdc), ]
# # check for outliers 
# samplesTreedc <- hclust(dist(t(signdc[, 2:29])), method = "average")
# sizeGrWindow(12,9)
# #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(samplesTreedc, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# # pheno.col <- numbers2colors(design, signed = FALSE)
# # 
# # plotDendroAndColors(samplesTree, pheno.col, 
# #                     groupLabels = colnames(design),
# #                     main = "Sample dendrogram and phenotype heatmap")
# 
# # Choose a set of soft-thresholding powers
# powers = c(1:20)
# # Call the network topology analysis function
# sft <- pickSoftThreshold(t(signdc[, 2:29]), powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 <- 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90, col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# 
# #wgcna
# powerdc <- 19
# wgcnaResdc <- blockwiseModules(t(signdc[, 2:29]), power = powerdc, TOMType = "unsigned", 
#                              minModuleSize = 5, verbose = 3)
# mergedColorsdc <- wgcnaResdc$colors
# plotDendroAndColors(wgcnaResdc$dendrograms[[1]],
#                     mergedColorsdc[wgcnaResdc$blockGenes[[1]]],
#                     "Module colors", dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# 
# 
# # Create list with all genes from the different modules (for easy lookup)
# moduleGenesdc <- split(signdc[, 1], wgcnaResdc$colors)
# 
# # Define modules of interest. For the sake of time, let's look at only three
# #  modules: blue, brown, pink
# moduleSelection <- c("blue", "turquoise")
# 
# #-----------------------------------------------------------------------------#
# # Export data for Cytoscape tutorial
# #-----------------------------------------------------------------------------#
# # Export the original input data with module membership
# cytoScapeAnnotation <- cbind(signdc[, 2:29], "WGCNA-Module" = wgcnaResdc$colors)
# write.table(cytoScapeAnnotation, "DataForCytoscapedc.txt", sep = "\t", quote = F, row.names = F)
# 
# TOM <- TOMsimilarityFromExpr(t(signdc[, 2:29]), power = powerdc)
# probes <- signdc$Genesymbol
# moduleColors <- labels2colors(wgcnaResdc$colors)
# moduleSelection <- c("blue", "turquoise")
# inModule = is.finite(match(moduleColors, moduleSelection))
# modProbes <- probes[inModule]
# modTOM <- TOM[inModule, inModule]
# dimnames(modTOM) <- list(modProbes, modProbes)
# cyt <- exportNetworkToCytoscape(modTOM,
#                                 edgeFile=paste("CytoscapeInput-edges-", paste(moduleSelection, collapse="-"), ".txt", sep=""),
#                                 nodeFile=paste("CytoscapeInput-nodes-", paste(moduleSelection, collapse="-"), ".txt", sep=""),
#                                 weighted=TRUE,
#                                 threshold=0.2,
#                                 nodeNames=modProbes,
#                                 altNodeNames=modProbes,
#                                 nodeAttr=moduleColors[inModule])
# 
