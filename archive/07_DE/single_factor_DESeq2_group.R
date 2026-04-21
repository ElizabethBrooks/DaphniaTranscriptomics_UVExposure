#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/DESeq2/group")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/DESeq2/group")
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/DESeq2/group")

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
# the values in the input matrix should be un-normalized counts or  
# estimated counts of sequencing reads (SE) or fragments (PE)

# import gene count data
#gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx/counts_merged.csv", row.names="gene"))
#gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_EGAPx_test/counts_merged.csv", row.names="gene"))
gene_counts <- as.matrix(read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/short_read_data_processed_test/counts_merged.csv", row.names="gene"))

# trim the data table to remove lines with counting statistics (htseq)
removeList <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
gene_counts <- gene_counts[!row.names(gene_counts) %in% removeList,]

# check out the number of imported genes
nrow(gene_counts)

# import grouping factors
colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/data/DESeq2/study_design_group.csv", row.names="sample")

# verify that the order of the samples in the counts and groupings files match
#colnames(gene_counts)
#rownames(colData)

# convert the grouping data into factors 
colData$group <- factor(colData$group)

# verify that the data is now of the factor type
#is.factor(colData$group)

# create DESeqDataSet list object
# the design formula expresses the variables which will be used in modeling
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = colData,
                              design = ~ group)

# inspect the list object
#dds

# by default, the reference level for factors is based on alphabetical order
# specify the reference level
dds$group <- relevel(dds$group, ref = "Olympic")

# verify the re-leveling
#dds$group


##
# Pre-Filtering
##
# pre-filter low count genes to improve performance and visualizations

# here there are 3 treated samples
smallestGroupSize <- 3

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

# filter the list object
dds <- dds[keep,]

# check out the number of kept genes
nrow(dds)


##
# Normalization
##
# the DESeq2 model internally corrects for library size, so transformed or 
# normalized values such as counts scaled by library size should not be used


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


##
# Data Exploration
##
# for example, variance stabilizing transformations (vst) produces
# transformed data on the log2 scale which has been normalized with respect 
# to library size or other normalization factors

# vst the data
vsd <- vst(dds, blind=FALSE)

# visualize the overall effect of the experimental groups or any batch effects
plotPCA(vsd, intgroup=c("group"))

# save the PCA
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=FALSE)

# customize the PCA
ggplot(pcaData@data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_colour_manual(values = c(ponyo_colors[4], ponyo_colors[3]))

# store the PCA plot
sample_pca <- ggplot(pcaData@data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(pcaData@labels$x) +
  ylab(pcaData@labels$y) + 
  coord_fixed() + 
  scale_colour_manual(values = c(ponyo_colors[4], ponyo_colors[3]))

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
# note that we have to provide a hierarchical clustering to the heatmap  
# function based on the sample distances, or else it would be based on the 
# distances between the rows/columns of the distance matrix
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
# with no additional arguments to results, the log2 fold change and Wald test
# p value will be for the last last level of the factored variable over the reference level

# standard differential expression analysis steps are wrapped into a single function
dds <- DESeq(dds)

# extract a results table with log2 fold changes, p values and adjusted p values
res <- results(dds)

# check out the results
# group Olympic vs Sierra indicates that the estimates are of the 
# logarithmic fold change log2(treated/untreated)
#res

# note that the order of the variables of the design do not matter 
# so long as the we directly specify the comparison
res <- results(dds, contrast=c("group","Sierra","Olympic"))

# check out the results
#res

# summarize some basic tallies
# notice that the default adjusted p-value cut off is 0.1
summary(res)

# retrieve information about which variables and tests were used
mcols(res)$description

# order our results table by the smallest p value
resOrdered <- res[order(res$pvalue),]

# save the ordered results to a csv file
write.csv(as.data.frame(resOrdered), file="Sierra_Olympic_results.csv")

# check the number of adjusted p-values were less than 0.05
sum(res$padj < 0.05, na.rm=TRUE)

# set the adjusted p-value cut off to 0.05
res05 <- results(dds, alpha=0.05)

# summarize the results
summary(res05)

# save the filtered results to a csv file
write.csv(as.data.frame(res05), file="Sierra_Olympic_results_0.5.csv")


##
# Results Exploration
##

# retrieve the ordered row means for the top 20 most abundant genes
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=FALSE)[1:20]

# setup list of colors associated with groups
anno_colors = list(group = c(Olympic = ponyo_colors[4], Sierra = ponyo_colors[3]))

# explore the count matrix using a heatmap of the vst data
pheatmap(assay(vsd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=colData,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(ponyo_colors[2], "white", ponyo_colors[5]))(100))

# store the clustering plot as a ggplot object
vst_clust <- as.ggplot(pheatmap(assay(vsd)[select,], 
                                cluster_rows=FALSE, 
                                show_rownames=FALSE,
                                cluster_cols=FALSE, 
                                annotation_col=colData,
                                annotation_colors = anno_colors,
                                color = colorRampPalette(c(ponyo_colors[2], "white", ponyo_colors[5]))(100)))

# save the plot to a png file
ggsave("vst_clustering.png", plot = vst_clust, bg = "white", device = "png")

# show the log2 fold changes attributable to a given variable over the mean 
# of normalized counts for all the samples
plotMA(res, ylim=c(-2,2))

# save the plot to a jpg file
jpeg("samples_log2fc.jpg")
plotMA(res, ylim=c(-2,2))
dev.off()

# shrink the log2 fold changes to remove the noise associated with log2 fold 
# changes from low count genes without requiring arbitrary filtering thresholds
resLFC <- lfcShrink(dds, coef="group_Sierra_vs_Olympic", type="apeglm")

# inspect the shrunken log2 fold changes
resLFC

# it is more useful to visualize the MA-plot for the shrunken log2 fold changes
plotMA(resLFC, ylim=c(-2,2))

# save the plot to a jpg file
jpeg("shrunken_log2fc.jpg")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# it can also be useful to examine the counts of reads for a single gene across the groups
plotCounts(dds, gene=which.min(res$padj), intgroup="group")

# save the plot of counts for a single gene across the groups
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="group", returnData=TRUE)

# customize the plot of counts for a single gene
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme_bw()

# store the plot of counts
gene_counts <- ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  theme_bw()

# save the plot to a png file
ggsave(paste(rownames(dds[which.min(res$padj)]), "gene_counts.png", sep = "_"), 
       plot = gene_counts, bg = "white", device = "png")
