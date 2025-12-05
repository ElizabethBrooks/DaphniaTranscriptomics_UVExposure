#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory for treatment analysis
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment/PA")

# set the working directory for treatment_genotype analysis
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype_noPA")
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype_TG_noPA")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/FDR0.1_LFC0.1")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/Olympic")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype/Sierra")

# set the working directory for treatment_group_batch analysis
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_group_batch")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_group_batch_TG")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_group_batch_TG_noPA")

# set the working directory for treatment_group analysis
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_group")

# set the working directory for treatment_batch analysis
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_batch")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_batch/Olympic")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_batch/Sierra")

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
blindColors <- c("#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", 
                 "#CC79A7", "#000000", "#999999")
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

# remove specific count data
#`%ni%` <- Negate(`%in%`)
#remove_list <- row.names(targets[grepl("Olympic", targets$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#targets <- targets[!grepl("Olympic", targets$group),]
#remove_list <- row.names(targets[grepl("Sierra", targets$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#targets <- targets[!grepl("Sierra", targets$group),]
gene_counts <- dplyr::select(gene_counts, -contains("PA"))
targets <- targets[!grepl("PA", targets$group),]

# remove unnecessary factors
targets <- dplyr::select(targets, -contains("group"))
targets <- dplyr::select(targets, -contains("batch"))
#targets <- dplyr::select(targets, -contains("tolerance"))
#targets <- dplyr::select(targets, -contains("genotype"))

# create DESeqDataSet list object
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = targets,
                              #design = ~ treatment + genotype)
                              design = ~ treatment + genotype + treatment:genotype)
                              #design = ~ treatment + group + batch)
                              #design = ~ treatment + group + batch + treatment:group)
                              #design = ~ treatment + group)
                              #design = ~ treatment + batch)
                              #design = ~ treatment)

# inspect the list object
#dds

# design to test for interactions as well
#dds$group <- factor(paste0(dds$treatment, dds$genotype))
#design(dds) <- ~ group

# specify the reference level
dds$treatment <- relevel(dds$treatment, ref = "VIS")
#dds$genotype <- relevel(dds$genotype, ref = "PA")
#dds$group <- relevel(dds$group, ref = "PA")
dds$genotype <- relevel(dds$genotype, ref = "Sierra")
#dds$group <- relevel(dds$group, ref = "Sierra")

# verify the re-leveling
#dds$treatment
#dds$genotype

# output the input data
gene_counts <- cbind(gene = row.names(gene_counts), gene_counts)
write.csv(as.data.frame(gene_counts), file="gene_counts.csv", quote = FALSE, row.names = FALSE)
targets <- cbind(sample = row.names(targets), targets)
write.csv(as.data.frame(targets), file="sample_data.csv", quote = FALSE, row.names = FALSE)

# clean up the targets
targets <- dplyr::select(targets, -contains("sample"))


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
pcaData <- plotPCA(vsd, intgroup=colnames(targets), returnData=FALSE)
# store the PCA plot of PC1 and PC2
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=pcaData@data[,6], shape=pcaData@data[,5])) +
  geom_point(size=3) + 
  scale_colour_manual(name = colnames(pcaData@data)[6], values = c(plotColors[seq(1, length(levels(pcaData@data[,6])))])) +
  scale_shape_manual(name = colnames(pcaData@data)[5], values = seq(0, length(levels(pcaData@data[,5]))-1)) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() +
  #stat_ellipse(aes(group = batch)) +
  theme_bw()
# save the PCA plot
ggsave("sample_pca.png", plot = sample_pca, device = "png", width = 9, height = 8, units = "in")

# save the PCA with PC1 and PC3
pcaData_pc1_pc3 <- plotPCA(vsd, intgroup=colnames(targets), returnData=FALSE, pcsToUse = c(1,3))
# store the PCA plot of PC1 and PC3
sample_pc1_pc3 <- ggplot(pcaData_pc1_pc3@data, aes(PC1, PC3, color=pcaData_pc1_pc3@data[,6], shape=pcaData_pc1_pc3@data[,5])) +
  geom_point(size=3) + 
  scale_colour_manual(name = colnames(pcaData_pc1_pc3@data)[6], values = c(plotColors[seq(1, length(levels(pcaData_pc1_pc3@data[,6])))])) +
  scale_shape_manual(name = colnames(pcaData_pc1_pc3@data)[5], values = seq(0, length(levels(pcaData_pc1_pc3@data[,5]))-1)) +
  xlab(pcaData_pc1_pc3@labels$x) +
  ylab(pcaData_pc1_pc3@labels$y) + 
  coord_fixed() +
  #stat_ellipse(aes(group = batch)) +
  theme_bw()
# save the PCA plot
ggsave("sample_pc1_pc3.png", plot = sample_pc1_pc3, device = "png", width = 9, height = 8, units = "in")

# save the PCA with PC2 and PC3
pcaData_pc2_pc3 <- plotPCA(vsd, intgroup=colnames(targets), returnData=FALSE, pcsToUse = c(2,3))
# store the PCA plot of PC2 and PC3
sample_pc2_pc3 <- ggplot(pcaData_pc2_pc3@data, aes(PC2, PC3, color=pcaData_pc2_pc3@data[,6], shape=pcaData_pc2_pc3@data[,5])) +
  geom_point(size=3) + 
  scale_colour_manual(name = colnames(pcaData_pc2_pc3@data)[6], values = c(plotColors[seq(1, length(levels(pcaData_pc2_pc3@data[,6])))])) +
  scale_shape_manual(name = colnames(pcaData_pc2_pc3@data)[5], values = seq(0, length(levels(pcaData_pc2_pc3@data[,5]))-1)) +
  xlab(pcaData_pc2_pc3@labels$x) +
  ylab(pcaData_pc2_pc3@labels$y) + 
  coord_fixed() +
  #stat_ellipse(aes(group = batch)) +
  theme_bw()
# save the PCA plot
ggsave("sample_pc2_pc3.png", plot = sample_pc2_pc3, device = "png", width = 9, height = 8, units = "in")

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

# LRT
#dds <- DESeq(dds, test="LRT", reduced=~batch)

# extract a results table with log2 fold changes, p values and adjusted p values
#res <- results(dds)

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

# view the group names
resultsNames(dds)

# set the adjusted p-value and LFC cut offs
# if lfcThreshold is specified, the results are for Wald tests, and LRT p-values will be overwritten
#res05 <- results(dds, contrast=c(colnames(targets)[2],levels(targets[,2])), alpha=cutFDR)
#res05 <- results(dds, contrast=c("treatment", "UV", "VIS"), alpha=cutFDR, lfcThreshold=cutLFC)

# testing comparisons
# https://rpubs.com/ge600/deseq2
# the effect of treatment in Sierra (the main effect)
res05 <- results(dds, contrast=c("treatment", "VIS", "UV"), alpha=cutFDR, lfcThreshold=cutLFC)
# the effect of treatment in Olympic
#res05 <- results(dds, contrast=list(c("treatment_UV_vs_VIS", "treatmentUV.groupOlympic")), alpha=cutFDR, lfcThreshold=cutLFC)
# what is the difference between Sierra and Olympic without treatment?
#res05 <- results(dds, contrast=c("group", "Sierra", "Olympic"), alpha=cutFDR, lfcThreshold=cutLFC)
# with treatment, what is the difference between Sierra and Olympic?
#res05 <- results(dds, contrast=list(c("group_Olympic_vs_Sierra", "treatmentUV.groupOlympic")), alpha=cutFDR, lfcThreshold=cutLFC)
# the different response in groups (interaction term)
#res05 <- results(dds, name="treatmentUV.groupOlympic")
# what is the difference between Old and New without treatment?
#res05 <- results(dds, contrast=c("batch", "Old", "New"), alpha=cutFDR, lfcThreshold=cutLFC)
# with treatment, what is the difference between Old and New?
#res05 <- results(dds, contrast=list(c("batch_Old_vs_New", "treatmentUV.groupOlympic")), alpha=cutFDR, lfcThreshold=cutLFC)

# LRT results
#res05 <- results(dds, alpha=cutFDR, lfcThreshold=cutLFC)

# order our results table by the smallest p value
res05Ordered <- res05[order(res05$pvalue),]

# summarize the results
summary(res05Ordered)
# number of up expressed
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[4], " ")[[1]][9])
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[4], " ")[[1]][12])
# number of down expressed
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[5], " ")[[1]][6])
#gsub(",", "", strsplit(capture.output(summary(res05Ordered))[5], " ")[[1]][10])

# format results
res05_out <- as.data.frame(res05Ordered)
res05_out <- cbind(gene = row.names(res05_out), res05_out)

# save the filtered results to a csv file
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
#exp_factor <- data.frame(
#                  Sample = unlist(paste(targets$genotype, 
#                                        targets$treatment, 
#                                        sep = "."), 
#                  use.names = FALSE))
#rownames(exp_factor) <- colnames(vsdSubset)

# setup annotation names
ann_groups <- data.frame(
                treatment = as.factor(targets$treatment),
                genotype = as.factor(targets$genotype),
                #group = as.factor(targets$group),
                #batch = as.factor(targets$batch),
                row.names = colnames(vsdSubset))

# setup annotation colors
ann_colors <- list(
  genotype = plotColors[seq(1, length(unique(targets$genotype)))],
  #batch = blindColors[seq(1, length(unique(targets$batch)))],
  #group = safeColors[seq(1, length(unique(targets$group)))],
  treatment = c("black", "white"))

# setup annotation color names
names(ann_colors[["treatment"]]) <- unique(targets$treatment)
names(ann_colors[["genotype"]]) <- unique(targets$genotype)
#names(ann_colors[["group"]]) <- unique(targets$group)
#names(ann_colors[["batch"]]) <- unique(targets$batch)


# create heatmap for DGE
vst_dge <- as.ggplot(
  pheatmap(vsdSubset, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           main="Heatmap of DE Genes", 
           show_rownames = FALSE,
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
#resLFC <- lfcShrink(dds, coef="treatment_VIS_vs_UV", type="apeglm")
resLFC <- lfcShrink(dds, coef="treatment_UV_vs_VIS", type="apeglm")

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
