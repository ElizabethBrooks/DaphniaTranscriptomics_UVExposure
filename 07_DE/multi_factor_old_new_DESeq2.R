#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype_noPA")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/FDR0.1_LFC0.1")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/Olympic")
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/Sierra")

##
# Packages
##

# un-comment to install packages, if necessary
#install.packages("ggplot2")
#install.packages("ggplotify")
#install.packages("pheatmap")
#install.packages("RColorBrewer")
#install.packages("rcartocolor")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")

# import libraries
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(rcartocolor)
library(DESeq2)
library(dplyr)

##
# Plotting Palettes
##

# TO-DO: double check
# color blind safe plotting palettes
defaultColors <- palette.colors(palette = "Okabe-Ito")
blindColors <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safeColors <- carto_pal(12, "Safe")
plotColors <- c(safeColors, blindColors, defaultColors)


##
# Data Setup
##

# import gene count data for new UV data
old_gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/data/OLYM_dMel_UV/cleaned.csv", row.names="gene"))

# import gene count data for new UV data
new_gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/counts_merged.csv", row.names="gene"))

# merge the old and new data
gene_counts <- merge(new_gene_counts, old_gene_counts, by = "row.names")
rownames(gene_counts) <- gene_counts$Row.names
gene_counts <- gene_counts[,-1]

# trim the data table to remove lines with counting statistics (htseq)
removeList <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
gene_counts <- gene_counts[!row.names(gene_counts) %in% removeList,]

# check out the number of imported genes
nrow(gene_counts)

# import grouping factors
#design <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment.csv", row.names="sample")
#design <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment_batch.csv", row.names="sample")
design <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design.csv", row.names="sample")

# verify that the order of the samples in the counts and groupings files match
#targets(gene_counts)
#rownames(design)

# convert the grouping data into factors 
targets <- as.data.frame(lapply(design, as.factor))
rownames(targets) <- rownames(design)

# remove data
`%ni%` <- Negate(`%in%`)
remove_list <- row.names(targets[grepl("Olympic", targets$group),])
gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
targets <- targets[!grepl("Olympic", targets$group),]
#remove_list <- row.names(targets[grepl("Sierra", targets$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#targets <- targets[!grepl("Sierra", targets$group),]
gene_counts <- dplyr::select(gene_counts, -contains("PA"))
targets <- targets[!grepl("PA", targets$group),]
targets <- dplyr::select(targets, -contains("group"))
targets <- dplyr::select(targets, -contains("batch"))
targets <- dplyr::select(targets, -contains("tolerance"))

# create DESeqDataSet list object
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = targets,
                              design = ~ treatment + genotype)
                              #design = ~ 0 + treatment + genotype)
                              #design = ~ group + treatment + genotype + batch + tolerance)

# inspect the list object
#dds

# design to test for interactions as well
#dds$group <- factor(paste0(dds$treatment, dds$genotype))
#design(dds) <- ~ group

# specify the reference level
#dds$treatment <- relevel(dds$treatment, ref = "VIS")
#dds$genotype <- relevel(dds$genotype, ref = "PA")

# verify the re-leveling
#dds$treatment
#dds$genotype

# output the input data
gene_counts <- cbind(gene = row.names(gene_counts), gene_counts)
write.csv(as.data.frame(gene_counts), file="gene_counts.csv", quote = FALSE, row.names = FALSE)
targets <- cbind(sample = row.names(targets), targets)
write.csv(as.data.frame(targets), file="sample_data.csv", quote = FALSE, row.names = FALSE)


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
# Data Exploration
##

# vst the data
vsd <- vst(dds, blind=FALSE)

# save the PCA
targets <- dplyr::select(targets, -contains("sample"))
pcaData <- plotPCA(vsd, intgroup=colnames(targets), returnData=FALSE)
                      
# store the PCA plot
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=pcaData@data[,5], shape=pcaData@data[,6])) +
  geom_point(size=3) + 
  scale_colour_manual(name = colnames(pcaData@data)[5], values = c(plotColors[seq(1, length(levels(pcaData@data[,5])))])) +
  scale_shape_manual(name = colnames(pcaData@data)[6], values = seq(0, length(levels(pcaData@data[,6]))-1)) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed()
# save the PCA plot
ggsave("sample_pca.png", plot = sample_pca, device = "png", width = 9, height = 8, units = "in")

# transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

# convert the distances to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# update the column names
colnames(sampleDistMatrix) <- NULL

# create a list of continuous colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# store the clustering plot as a ggplot object
sample_clust <- as.ggplot(pheatmap(sampleDistMatrix,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         col=colors))
# save the plot to a png file
ggsave("sample_clustering.png", plot = sample_clust, bg = "white", device = "png", width = 9, height = 8, units = "in")


##
# DE Analysis Contrasts
##

# standard differential expression analysis steps are wrapped into a single function
dds <- DESeq(dds)

# extract a results table with log2 fold changes, p values and adjusted p values
#res <- results(dds)

# view the group names
#resultsNames(dds)

# directly specify the comparison
#res <- results(dds, contrast=c(colnames(targets)[2],levels(targets[,2])))

# summarize some basic tallies
#summary(res)

# retrieve information about which variables and tests were used
#mcols(res)$description

# order our results table by the smallest p value
#resOrdered <- res[order(res$pvalue),]

# save the ordered results to a csv file
#resOrdered_out <- as.data.frame(resOrdered)
#resOrdered_out <- cbind(gene = row.names(resOrdered_out), resOrdered_out)
#write.csv(resOrdered_out, file="UV_VIS_results.csv", quote = FALSE, row.names = FALSE)

# check the number of adjusted p-values were less than 0.05
#sum(res$padj < 0.05, na.rm=TRUE)

# thresholds
cutFDR=0.05
cutLFC=0

# set the adjusted p-value cut off to 0.05
#res05 <- results(dds, contrast=c(colnames(targets)[2],levels(targets[,2])), alpha=cutFDR)
# set the adjusted p-value cut off to 0.05 and LFC to 1.2
# If lfcThreshold is specified, the results are for Wald tests, and LRT p-values will be overwritten.
res05 <- results(dds, contrast=c(colnames(targets)[2],levels(targets[,2])), alpha=cutFDR, lfcThreshold=cutLFC)

# order our results table by the smallest p value
res05Ordered <- res05[order(res05$pvalue),]

# summarize the results
summary(res05Ordered)
# number of up expressed
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[4], " ")[[1]][9])
gsub(",", "", strsplit(capture.output(summary(res05Ordered))[4], " ")[[1]][12])
# number of down expressed
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[5], " ")[[1]][6])
gsub(",", "", strsplit(capture.output(summary(res05Ordered))[5], " ")[[1]][10])

# save the filtered results to a csv file
res05_out <- as.data.frame(res05Ordered)
res05_out <- cbind(gene = row.names(res05_out), res05_out)
write.csv(res05_out, file="UV_VIS_results.csv", quote = FALSE, row.names = FALSE)


##
# Results Exploration
##

# identify significantly DE genes
DGESubset <- na.omit(res05_out[res05_out$padj <= cutFDR,])
DGESubset <- na.omit(DGESubset[DGESubset$log2FoldChange >= cutLFC | DGESubset$log2FoldChange <= (-1*cutLFC),])
# extract vst counts
vsdCounts <- assay(vsd)

# subset the vst data by the DGE set
DGESubset.keep <- rownames(vsdCounts) %in% rownames(DGESubset)
vsdSubset <- vsdCounts[DGESubset.keep, ]

# combine all columns into one period separated
#exp_factor <- data.frame(Sample = unlist(targets, use.names = FALSE))
#rownames(exp_factor) <- colnames(vsdSubset)

# TO-DO: use color blind safe palette for sample dendrogram
# create heatmap for DGE
vst_dge <- as.ggplot(
  pheatmap(vsdSubset, scale="row", #annotation_col = exp_factor, 
           main="Heatmap of GLM DE Genes", show_rownames = FALSE,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[6]))(100))
)
# save the plot to a png file
ggsave("vst_dge.png", plot = vst_dge, bg = "white", device = "png", width = 9, height = 8, units = "in")

# retrieve the ordered row means for the top 20 most abundant genes
#select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=FALSE)[1:20]

# setup list of colors associated with treatments
#anno_colors = list(
#  treatment = c(VIS = plotColors[1], UV = plotColors[2]),
#  genotype = c("CON_6" = plotColors[12], "E2_1" = plotColors[11], 
#               "GRO_3" = plotColors[10], "NGD_1" = plotColors[9], 
#               "Y002_3_2" = plotColors[8], "Y019_2_1" = plotColors[7],
#               "Y05" = plotColors[6], "Y023_5" = plotColors[5],
#               "E05" = plotColors[4], "R2" = plotColors[3],
#               "PA" = plotColors[2], 
#               "Sierra" = plotColors[1]))

# store the clustering plot as a ggplot object
#vst_clust <- as.ggplot(pheatmap(assay(vsd)[select,],
#                                   cluster_rows=FALSE, 
                                   #show_rownames=FALSE,
#                                   cluster_cols=FALSE, 
#                                   annotation_col=targets,
#                                   annotation_colors = anno_colors,
#                                   color = colorRampPalette(c(plotColors[4], "white", plotColors[3]))(100)))
# save the plot to a png file
#ggsave("vst_clustering.png", plot = vst_clust, bg = "white", device = "png", width = 9, height = 8, units = "in")

# show the log2 fold changes attributable to a given variable over the mean 
# save the plot to a png file
png("samples_log2fc.png")
plotMA(res05, ylim=c(-2,2))
dev.off()

# shrink the log2 fold changes to remove the noise 
resLFC <- lfcShrink(dds, coef="treatment_VIS_vs_UV", type="apeglm")
#resLFC <- lfcShrink(dds, coef="treatment_UV_vs_VIS", type="apeglm")
# inspect the shrunken log2 fold changes
#resLFC

# it is more useful to visualize the MA-plot for the shrunken log2 fold changes
# save the plot to a png file
png("shrunken_log2fc.png")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# it can also be useful to examine the counts of reads for a single gene across the groups
plotCounts(dds, gene=which.min(res05$padj), intgroup="treatment")

# save the plot of counts for a single gene across the groups
d <- plotCounts(dds, gene=which.min(res05$padj), intgroup="treatment", returnData=TRUE)
#d <- plotCounts(dds, gene="gene-LOC124188748", intgroup="treatment", returnData=TRUE)

# store the plot of counts
gene_counts <- ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme_bw()

# save the plot to a png file
ggsave(paste(rownames(dds[which.min(res05$padj)]), "gene_counts.png", sep = "_"), 
       plot = gene_counts, bg = "white", device = "png")
#ggsave(paste(rownames(dds["gene-LOC124188748"]), "LOC124188748_gene_counts.png", sep = "_"), 
#       plot = gene_counts, bg = "white", device = "png")
