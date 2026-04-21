#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/DESeq2/treatment_genotype")
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/DESeq2/treatment_genotype")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/DESeq2/treatment_genotype")

##
# Packages
##

# un-comment to install packages, if necessary
#install.packages("ggplot2")
#install.packages("ggplotify")
#install.packages("pheatmap")
#install.packages("RColorBrewer")
#install.packages("ghibli")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")

# import libraries
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(ghibli)
library(DESeq2)


##
# Data Setup
##

# import gene count data
#gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/counts_merged.csv", row.names="gene"))
gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/counts_merged.csv", row.names="gene"))
#gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/counts_merged.csv", row.names="gene"))

# trim the data table to remove lines with counting statistics (htseq)
removeList <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
gene_counts <- gene_counts[!row.names(gene_counts) %in% removeList,]

# check out the number of imported genes
nrow(gene_counts)

# import grouping factors
#colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/data/DESeq2/study_design_genotype_treatment.csv", row.names="sample")
colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment.csv", row.names="sample")

# verify that the order of the samples in the counts and groupings files match
#colnames(gene_counts)
#rownames(colData)

# convert the grouping data into factors 
colData$treatment <- factor(colData$treatment)
colData$genotype <- factor(colData$genotype)

# verify that the data is now of the factor type
#is.factor(colData$treatment)
#is.factor(colData$genotype)

# create DESeqDataSet list object
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = colData,
                              design = ~ treatment + genotype)

# inspect the list object
#dds

# specify the reference level
dds$treatment <- relevel(dds$treatment, ref = "Control")
dds$genotype <- relevel(dds$genotype, ref = "CON_6")

# verify the re-leveling
#dds$treatment
#dds$genotype


##
# Pre-Filtering
##

# here there are 3 treated samples
smallestGroupSize <- 3

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

# filter the list object
dds <- dds[keep,]

# check out the number of kept genes
nrow(dds)


##
# Plotting Palettes
##

# change the graphical parameters
#par(mfrow=c(9,3))

# view all available ghibli palettes
#for(i in names(ghibli_palettes)) print(ghibli_palette(i))

# close the plot and return the display to the default graphical parameters
#dev.off()

# retrieve the vector of colors associated with PonyoMedium
ponyo_colors <- ghibli_palette("PonyoMedium", type = "discrete")

# retrieve the vector of colors associated with YesterdayMedium
yest_colors <- ghibli_palette("YesterdayMedium", type = "discrete")

# retrieve the vector of colors associated with KikiMedium
kiki_colors <- ghibli_palette("KikiMedium", type = "discrete")


##
# Data Exploration
##

# vst the data
vsd <- vst(dds, blind=FALSE)

# visualize the overall effect of the experimental treatments or any batch effects
plotPCA(vsd, intgroup=c("treatment", "genotype"))

# save the PCA
pcaData <- plotPCA(vsd, intgroup=c("treatment", "genotype"), returnData=FALSE)

# customize the PCA
ggplot(pcaData@data, aes(PC1, PC2, color=treatment, shape=genotype)) +
  geom_point(size=3) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_colour_manual(values = c(ponyo_colors[4], ponyo_colors[3], ponyo_colors[2], ponyo_colors[1]))

# store the PCA plot
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=treatment, shape=genotype)) +
  geom_point(size=3) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_colour_manual(values = c(ponyo_colors[4], ponyo_colors[3], ponyo_colors[2], ponyo_colors[1]))

# save the PCA plot
ggsave("sample_pca.png", plot = sample_pca, device = "png")

# transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

# convert the distances to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# update the column names
colnames(sampleDistMatrix) <- NULL

# create a list of continuous colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# view the similarities and dissimilarities between samples
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# store the clustering plot as a ggplot object
sample_clust <- as.ggplot(pheatmap(sampleDistMatrix,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         col=colors))

# save the plot to a png file
ggsave("sample_clustering.png", plot = sample_clust, bg = "white", device = "png")


##
# DE Analysis Contrasts
##

# standard differential expression analysis steps are wrapped into a single function
dds <- DESeq(dds)

# extract a results table with log2 fold changes, p values and adjusted p values
res <- results(dds)

# check out the results
#res

# directly specify the comparison
res <- results(dds, contrast=c("treatment","UV","Control"))

# check out the results
#res

# summarize some basic tallies
summary(res)

# retrieve information about which variables and tests were used
mcols(res)$description

# order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

# save the ordered results to a csv file
write.csv(as.data.frame(resOrdered), file="UV_Control_results.csv")

# check the number of adjusted p-values were less than 0.5
sum(res$padj < 0.5, na.rm=TRUE)

# set the adjusted p-value cut off to 0.5
res05 <- results(dds, alpha=0.5)

# summarize the results
summary(res05)

# save the filtered results to a csv file
write.csv(as.data.frame(res05), file="UV_Control_results_0.5.csv")

##
# Results Exploration
##

# retrieve the ordered row means for the top 20 most abundant genes
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=FALSE)[1:20]

# setup list of colors associated with treatments
anno_colors = list(
  treatment = c(Control = ponyo_colors[4], UV = ponyo_colors[3]),
  genotype = c("CON_6" = yest_colors[8], "E2_1" = yest_colors[7], 
               "GRO_3" = yest_colors[6], "NGD_1" = yest_colors[5], 
               "Y002_3_2" = yest_colors[4], "Y019_2_1" = yest_colors[3]))

# explore the count matrix using a heatmap of the vst data
pheatmap(assay(vsd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=colData,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(kiki_colors[4], "white", kiki_colors[3]))(100))

# store the clustering plot as a ggplot object
vst_clust <- as.ggplot(pheatmap(assay(vsd)[select,], 
                                   cluster_rows=FALSE, 
                                   show_rownames=FALSE,
                                   cluster_cols=FALSE, 
                                   annotation_col=colData,
                                   annotation_colors = anno_colors,
                                   color = colorRampPalette(c(kiki_colors[4], "white", kiki_colors[3]))(100)))

# save the plot to a png file
ggsave("vst_clustering.png", plot = vst_clust, bg = "white", device = "png")

# show the log2 fold changes attributable to a given variable over the mean 
plotMA(res, ylim=c(-2,2))

# save the plot to a jpg file
jpeg("samples_log2fc.jpg")
plotMA(res, ylim=c(-2,2))
dev.off()

# shrink the log2 fold changes to remove the noise 
resLFC <- lfcShrink(dds, coef="treatment_UV_vs_Control", type="apeglm")

# inspect the shrunken log2 fold changes
#resLFC

# it is more useful to visualize the MA-plot for the shrunken log2 fold changes
plotMA(resLFC, ylim=c(-2,2))

# save the plot to a jpg file
jpeg("shrunken_log2fc.jpg")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# it can also be useful to examine the counts of reads for a single gene across the groups
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")

# save the plot of counts for a single gene across the groups
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", returnData=TRUE)

# customize the plot of counts for a single gene
ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme_bw()

# store the plot of counts
gene_counts <- ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme_bw()

# save the plot to a png file
ggsave(paste(rownames(dds[which.min(res$padj)]), "gene_counts.png", sep = "_"), 
       plot = gene_counts, bg = "white", device = "png")
