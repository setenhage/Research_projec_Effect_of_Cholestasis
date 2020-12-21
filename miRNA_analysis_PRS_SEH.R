#=========================================================
# Analysis of miRNA Microarray data of liver tissue samples
# 
#
# Version 1.0
# Course: PRO4002 Research project 1 MSB
# Names: S.E. ten Hage, P.R. Stolpe, 
#        MaCSBio, Maastricht University
# ID: i6201596, i6195641
#
#
# History: 
#    1.0: Creation
#==========================================================
# 0. load required packages
#==========================================================
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "biomaRt",
                       "microRNA", "mirbase.db",
                       "multiMiR", "org.Hs.eg.db"))
require(ggplot2)
require(gplots)
require(limma)
require(pcaMethods)
require(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(microRNA)
library(mirbase.db)
library(multiMiR)
options(stringsAsFactors = F)
#
#
#==========================================================
# 1. Import data files
#==========================================================
# define directories according to where the files are in and
# results are supposed to be found after analysis
DATA.DIR <- "C:\\Users\\Suzanne\\Documents\\Studie\\Systems Biology\\PRO4002 - Cholestasis project\\Data"
RESULTS.DIR <- "C:\\Users\\Suzanne\\Documents\\Studie\\Systems Biology\\PRO4002 - Cholestasis project\\Results"

# define working directory
setwd(DATA.DIR)

# 1.1 import the normalized gene expression and miRNA data into an object: Hepa.norm
#----------------------------------------------------------
#norm.miRNA <- read.delim("")
norm.miRNA <- read.delim("miRNAexpression.txt", header = TRUE,
                         row.names = 1)
# import differentially expressed mRNA gens cholestatic v control
sig.mRNA.ConChol<- read.delim("convcholsignificant.txt", header = TRUE,
                       row.names = 1)
# import differentially expressed mRNA genes cholestatic v drained
sig.mRNA.DrainChol<- read.delim("Gdrainedvcholsignificant.txt", header = TRUE,
                               row.names = 1)

#norm.mRNA.ConDrain<- read.delim("GeneExpressionLimmaprocessed.txt", header = TRUE,
#                               row.names = 1)

#==========================================================
# 2. inspect the miRNA data object & remove batch effect
#==========================================================

# replace NaN with NA to be readable for PCA function
norm.miRNA[is.na(norm.miRNA)] <- NA

# Design matrix of the experimental design
design.miRNA <- model.matrix(~0+factor(c(1,1,1,1,1,1,1,1,1,
                                   2,2,2,2,2,2,2,2,2,2,
                                   3,3,3,3,3,3,3,3,3)))

colnames(design.miRNA) <- c("cholestatic", 
                      "drained",
                      "control")
#remove batch effect
#plotMDS(norm.miRNA, labels = c(1,1,1,1,1,1,1,1,1,
#                               2,2,2,2,2,2,2,2,2,2,
#                               3,3,3,3,3,3,3,3,3),
#       gene.selection = "common")

#title("plot before removal of batch effect")

#batch.miRNA <- c(rep("A", 8), rep("B", 8), rep("C", 8), rep("D", 4))
#batch.miRNA <- c("A", "A", "A", "B", "B", "B", "C", "C", "D",
#           "A", "A", "A", "B", "B", "B", "C", "C", "C", "D",
#           "A", "A", "B", "B", "C", "C", "C", "D", "D")
#
#norm.miRNA <- removeBatchEffect(norm.miRNA,
#                                batch = batch.miRNA,
#                                 design = design.miRNA)

#plotMDS(norm.miRNA, labels = c(1,1,1,1,1,1,1,1,1,
#                                          2,2,2,2,2,2,2,2,2,2,
#                                          3,3,3,3,3,3,3,3,3),
#        gene.selection = "common")
#title("plot after removal of batch effect")

#==============================================================================
# 3. create diagnostic plots of the data 
#==============================================================================

# 3.1 create boxplot of microRNA data
#------------------------------------------------------------------------------
setwd(RESULTS.DIR)

pdf("visualization_of_normalized_data miRNA.pdf", 7, 5)
dev.off()
# create boxplot of normalized gene expression data
boxplot(norm.miRNA, ylab = "gene expression levels", 
        col = c(rep("lightblue", 9), rep("orange", 10), 
                rep("darkgreen", 9)), 
        outline = F)
legend("topright", legend = c("cholestatic", "drained", "control"),
       fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)
title("Boxplot of normalized gene expression data")

# 3.2 PCA of normalized data 
#----------------------------------------------------------

# again check for correct dimensions
dim(norm.miRNA)

# pca of Hepa.norm.shuf using ppca: pcaRes
pcaRes.miRNA <- pca(t(norm.miRNA), 
              method = "ppca", 
              center = TRUE,  
              nPcs = 10)

# plot Pca Results 
plot(pcaRes.miRNA, main = "PCA of normalized data")

# make plots of principal components agains each other
plotPcs(pcaRes.miRNA, c(1,2), type = "scores", sl = NULL, 
        col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
legend("topright", legend = c("cholestatic", "drained", "control"),
       fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)

plotPcs(pcaRes.miRNA, c(1,3), sl = NULL, 
        col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
legend("topright", legend = c("cholestatic", "drained", "control"),
       fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)

plotPcs(pcaRes.miRNA, c(2,3), sl = NULL, 
        col = c(rep("lightblue", 9), rep("orange", 10), rep("darkgreen", 9)))
legend("topright", legend = c("cholestatic", "drained", "control"),
       fill = c("lightblue", "orange", "darkgreen"), ncol = 1, cex = 0.65)

dev.off()


#==========================================================
# 4. statistical analysis using limma package
#==========================================================
# 4.1 Statistical analysis
#---------------------------------------------------------------------

# use lmFit to create statistical model of microRNA expression according to
# design matrix
norm.fit.miRNA <- lmFit(norm.miRNA,
                        design = design.miRNA)

# build contrast matrix of the desired comparisons
contrast.matrix <- makeContrasts(cholestatic-control, 
                    cholestatic-drained, 
                    drained-control, 
                    levels = design.miRNA)
  


# make contrast fit to submit for further analysis of significant
# differences in gene expression levels
norm.fit2.contrast.miRNA <- contrasts.fit(norm.fit.miRNA, contrast.matrix)
# compute t-statistics, F-statistics and log-odds using eBayes
norm.eBayes.miRNA <- eBayes(norm.fit2.contrast.miRNA)
# create an MA-plot 
limma::plotMA(norm.eBayes.miRNA,
              xlab = "Average log-expression", 
              ylab = "sqrt(sigma)", 
              zero.weights = FALSE,
              pch = 16, 
              cex = 0.3 
              )
# produce object to identify differential expressed genes unsing decideTests
norm.decide.miRNA <- decideTests(norm.eBayes.miRNA)

#which <- which(norm.decide.miRNA != 0, arr.ind = T)
#signif.miRNA <- rownames(norm.decide.miRNA[which[, 1], ])

# assign statistical analysis: stats
stats.miRNA <- topTable(norm.eBayes.miRNA,
                        number = nrow(norm.eBayes.miRNA),
                        adjust.method = "BH",
                        sort.by = "none")
dim(stats.miRNA)

stats.miRNA.conchol <- topTable(norm.eBayes.miRNA,
                                coef = 1,
                                number = nrow(norm.eBayes.miRNA),
                                adjust.method = "BH",
                                sort.by = "none")

stats.miRNA.CholDrain <- topTable(norm.eBayes.miRNA,
                                coef = 2,
                                number = nrow(norm.eBayes.miRNA),
                                adjust.method = "BH",
                                sort.by = "none")

stats.miRNA.ConDrain <- topTable(norm.eBayes.miRNA,
                                  coef = 3,
                                  number = nrow(norm.eBayes.miRNA),
                                  adjust.method = "BH",
                                  sort.by = "none")

# produce histogram of p values 
pdf("Histogram and volcano plot miRNA.pdf", 7,5)
dev.off()

hist(stats.miRNA[, "P.Value"], xlab = "P.Value", 
     main = "Histogram of P-Values")

# produce volcanoplot and highlight 5 most significant genes with names
volcanoplot(norm.eBayes.miRNA, highlight = 5,
            names = rownames(norm.miRNA),
            xlab = "Log2 FC", 
            ylab = "-Log10(pValue)", 
            main = "Volcanoplot of P.Values")
dev.off()


# combine statistical analysis with data frame
norm.miRNA <- data.frame(norm.miRNA, 
                         adj.P.Values = stats.miRNA$adj.P.Val,
                         FC = norm.fit2.contrast.miRNA$coefficients)

# creating venn diagrams of statistical analysis
pdf("Venn Diagram miRNA.pdf", 10, 9)   
dev.off()

# VennDiagram to visualize number of up and downregulated 
# genes in each comparison
vennDiagram(norm.decide.miRNA, 
            include = c("up", "down"), 
            counts.col = c("red","green"), 
            cex = 1.25)
title("Venn Diagram of up and downregulated genes")

dev.off()

#==========================================================
# 5. Selecting significantly changed miRNAs
#==========================================================
#select significantly changed miRNAs 

#control vs. cholestatic -> 4
sig.miRNA.ConChol <- stats.miRNA.conchol[rownames(stats.miRNA.conchol)[
  which(stats.miRNA.conchol$adj.P.Val < 0.05 & 
          (stats.miRNA.conchol$logFC > 1.5 |
             stats.miRNA.conchol$logFC < -1.5))], ]

#drained vs. cholestatic -> 0 -> not further analyzed. 
sig.miRNA.CholDrain <- stats.miRNA.CholDrain[rownames(stats.miRNA.CholDrain)[
  which(stats.miRNA.CholDrain$adj.P.Val < 0.05 & 
          (stats.miRNA.CholDrain$logFC > 1.5 |
             stats.miRNA.CholDrain$logFC < -1.5))], ]

#control vs. drained -> 0 -> not further analyzed. 
sig.miRNA.ConDrain <- stats.miRNA.ConDrain[rownames(stats.miRNA.ConDrain)[
  which(stats.miRNA.ConDrain$adj.P.Val < 0.05 & 
          (stats.miRNA.ConDrain$logFC > 1.5 |
             stats.miRNA.ConDrain$logFC < -1.5))], ]

#==========================================================
# 6. Finding targets
#==========================================================
#Find miRNA - gene pairs of significantly changed miRNAs and genes.

#First make a vector of entrez symbols of all significantly changed mRNAs

sig.mRNA.ConChol.entrez <- row.names(sig.mRNA.ConChol)

#Control vs. Cholestatic: 5.  using a predicted cutoff of 20 = same result. 
gene.miRNA.int.ConChol <- get_multimir(org = "hsa",
                             mirna = rownames(sig.miRNA.ConChol), 
                             target = sig.mRNA.ConChol.entrez,
                             table = 'all', 
                             summary = TRUE,
                             predicted.cutoff.type = "p",
                             predicted.cutoff = 10,
                             use.tibble = TRUE
                           )

table(gene.miRNA.int.ConChol@data$type)

miRNA.gene.pairs.ConChol <- select(gene.miRNA.int.ConChol, 
                 keytype = "type", 
                 keys = "validated", 
                 columns = columns(gene.miRNA.int.ConChol))

#Select unique pairs from result 
unique.pairs.ConChol <- miRNA.gene.pairs.ConChol[!duplicated(
  miRNA.gene.pairs.ConChol[, c("mature_mirna_id", "target_entrez")]), ]

setwd(RESULTS.DIR)
write.table(unique.pairs.ConChol, "miRNA mRNA gene pairs. ConChol", 
          append = FALSE, sep = "\t", dec = ".",
          row.names = TRUE, col.names = TRUE)


#==========================================================
# 7. Analyse long and short drainage time seperately. 
#==========================================================
design.miRNA.2 <- model.matrix(~0+factor(c(1,1,1,1,1,1,1,1,1,
                                         2,2,2,2,2,2,
                                         3,3,3,
                                         2,
                                         4,4,4,4,4,4,4,4,4)
                                         )
                               )

colnames(design.miRNA.2) <- c("cholestatic", 
                            "drained.low",
                            "drained.high",
                            "control")

norm.fit.miRNA.2 <- lmFit(norm.miRNA[ , 1:28],
                        design = design.miRNA.2)

# build contrast matrix of the desired comparisons
contrast.matrix.2 <- makeContrasts(cholestatic-control, 
                                 cholestatic-drained.low,
                                 cholestatic-drained.high,
                                 drained.low-control,
                                 drained.high-control,
                                 levels = design.miRNA.2)

# make contrast fit to submit for further analysis of significant
# differences in gene expression levels
norm.fit2.contrast.miRNA.2 <- contrasts.fit(norm.fit.miRNA.2, contrast.matrix.2)
# compute t-statistics, F-statistics and log-odds using eBayes
norm.eBayes.miRNA.2 <- eBayes(norm.fit2.contrast.miRNA.2)
# produce object to identify differential expressed genes unsing decideTests
norm.decide.miRNA.2 <- decideTests(norm.eBayes.miRNA.2)

# attempt to select significant genes 
#which <- which(norm.decide.miRNA != 0, arr.ind = T)
#signif.miRNA <- rownames(norm.decide.miRNA[which[, 1], ])

# assign statistical analysis: stats
stats.miRNA.2 <- topTable(norm.eBayes.miRNA.2, 
                        number = nrow(norm.eBayes.miRNA.2),
                        adjust.method = "bonferroni",
                        sort.by = "none")
dim(stats.miRNA)

# produce histogram of p values 
pdf("Histogram and volcano plot miRNA.pdf", 7,5)
dev.off()

hist(stats.miRNA.2[, "P.Value"], xlab = "P.Value", 
     main = "Histogram of P-Values")

#find the significant miRNA
sig.miRNA.ConChol.2 <- stats.miRNA.2[rownames(stats.miRNA.2)[
  which(stats.miRNA.2$adj.P.Val < 0.05 & 
          (stats.miRNA.2$cholestatic...control > 1.7 |
             stats.miRNA.2$cholestatic...control < -1.7))], ]

View(sig.miRNA.ConChol.2)

#No clear difference in expression between long (>15) and short
#drainage time for this miRNA. Power is to low to analyse
#anything else (n = 3 in long drainage time). 

#==========================================================
# 8. Export files
#==========================================================

write.csv(sig.mRNA.genes, "significant mRNA genes", 
          append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE)




