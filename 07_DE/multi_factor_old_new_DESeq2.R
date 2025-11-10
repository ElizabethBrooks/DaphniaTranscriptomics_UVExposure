#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype_noPA")
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DESeq2/treatment_genotype")

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
library(dplyr)

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
#targets <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment.csv", row.names="sample")
#targets <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment_batch.csv", row.names="sample")
targets <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design.csv", row.names="sample")

# verify that the order of the samples in the counts and groupings files match
#colnames(gene_counts)
#rownames(targets)

# convert the grouping data into factors 
colData <- as.data.frame(lapply(targets, as.factor))
rownames(colData) <- rownames(targets)

# remove data
#`%ni%` <- Negate(`%in%`)
#remove_list <- row.names(colData[grepl("E2_1", colData$genotype),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#colData <- colData[!grepl("E2_1", colData$genotype),]
#remove_list <- row.names(colData[grepl("Sierra", colData$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#colData <- colData[!grepl("Sierra", colData$group),]
#gene_counts <- select(gene_counts, -contains("PA"))
#colData <- colData[!grepl("PA", colData$group),]
colData <- dplyr::select(colData, -contains("group"))
colData <- dplyr::select(colData, -contains("batch"))
colData <- dplyr::select(colData, -contains("tolerance"))

# create DESeqDataSet list object
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = colData,
                              design = ~ treatment + genotype)
                              #design = ~ 0 + treatment + genotype)
                              #design = ~ group + treatment + genotype + batch + tolerance)

# inspect the list object
dds

# specify the reference level
dds$treatment <- relevel(dds$treatment, ref = "VIS")
#dds$genotype <- relevel(dds$genotype, ref = "Sierra")
dds$genotype <- relevel(dds$genotype, ref = "PA")

# verify the re-leveling
#dds$treatment
#dds$genotype

# output the input data
gene_counts <- cbind(gene = row.names(gene_counts), gene_counts)
write.csv(as.data.frame(gene_counts), file="gene_counts.csv", quote = FALSE, row.names = FALSE)
colData <- cbind(sample = row.names(colData), colData)
write.csv(as.data.frame(colData), file="sample_data.csv", quote = FALSE, row.names = FALSE)


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

# visualize the overall effect of the experimental treatments or any batch effects
#plotPCA(vsd, intgroup=c("treatment", "genotype"))

# save the PCA
pcaData <- plotPCA(vsd, intgroup=c("treatment", "genotype"), returnData=FALSE)
                      
# store the PCA plot
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=treatment, shape=genotype)) +
  geom_point(size=3) +
  #scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_colour_manual(values = c(ponyo_colors[4], ponyo_colors[3]))

# view the pca
sample_pca

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
ggsave("sample_clustering.png", plot = sample_clust, bg = "white", device = "png", width = 9, height = 8, units = "in")


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
res <- results(dds, contrast=c("treatment","UV","VIS"))

# check out the results
#res

# summarize some basic tallies
summary(res)

# retrieve information about which variables and tests were used
mcols(res)$description

# order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

# save the ordered results to a csv file
resOrdered_out <- as.data.frame(resOrdered)
resOrdered_out <- cbind(gene = row.names(resOrdered_out), resOrdered_out)
write.csv(resOrdered_out, file="UV_VIS_results.csv", quote = FALSE, row.names = FALSE)

# check the number of adjusted p-values were less than 0.05
sum(res$padj < 0.05, na.rm=TRUE)

# set the adjusted p-value cut off to 0.05 and LFC to 1.2
res05 <- results(dds, alpha=0.05, lfcThreshold=1.2)

# summarize the results
summary(res05)
#summary(as.data.frame(res05))

# save the filtered results to a csv file
res05_out <- as.data.frame(res05)
res05_out <- cbind(gene = row.names(res05_out), res05_out)
write.csv(res05_out, file="UV_VIS_results_FDR0.05_LFC1.2.csv", quote = FALSE, row.names = FALSE)


##
# Results Exploration
##

# retrieve the ordered row means for the top 20 most abundant genes
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=FALSE)[1:20]

# setup list of colors associated with treatments
anno_colors = list(
  treatment = c(VIS = ponyo_colors[4], UV = ponyo_colors[3]),
  genotype = c("CON_6" = "white", "E2_1" = yest_colors[7], 
               "GRO_3" = yest_colors[6], "NGD_1" = yest_colors[5], 
               "Y002_3_2" = yest_colors[4], "Y019_2_1" = yest_colors[3],
               "Y05" = yest_colors[2], "Y023_5" = yest_colors[1],
               "E05" = ponyo_colors[7], "R2" = ponyo_colors[6],
               "PA" = ponyo_colors[5], 
               "Sierra" = ponyo_colors[2]))

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
                                   #show_rownames=FALSE,
                                   cluster_cols=FALSE, 
                                   annotation_col=colData,
                                   annotation_colors = anno_colors,
                                   color = colorRampPalette(c(kiki_colors[4], "white", kiki_colors[3]))(100)))

# save the plot to a png file
ggsave("vst_clustering.png", plot = vst_clust, bg = "white", device = "png", width = 9, height = 8, units = "in")

# show the log2 fold changes attributable to a given variable over the mean 
plotMA(res, ylim=c(-2,2))

# save the plot to a jpg file
jpeg("samples_log2fc.jpg")
plotMA(res, ylim=c(-2,2))
dev.off()

# shrink the log2 fold changes to remove the noise 
resLFC <- lfcShrink(dds, coef="treatment_UV_vs_VIS", type="apeglm")

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
#d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", returnData=TRUE)
d <- plotCounts(dds, gene="gene-LOC124188748", intgroup="treatment", returnData=TRUE)

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
#ggsave(paste(rownames(dds[which.min(res$padj)]), "gene_counts.png", sep = "_"), 
#       plot = gene_counts, bg = "white", device = "png")
ggsave(paste(rownames(dds["gene-LOC124188748"]), "LOC124188748_gene_counts.png", sep = "_"), 
       plot = gene_counts, bg = "white", device = "png")
