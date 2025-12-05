#!/usr/bin/env Rscript

# edgeR user guide
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

##
# Working Directory
##

# set the working directory
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype/Olympic")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype/Sierra")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype/PA")
#setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype_noPA")


##
# Packages
##

# install packages, if necessary
#install.packages("ggplot2")
#install.packages("ghibli")
#install.packages("ggVennDiagram")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

# import libraries
library(ggplot2)
library(ghibli)
library(ggVennDiagram)
library(edgeR)
library(dplyr)
library(rcartocolor)
library(ggplotify)
library(pheatmap)


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
# Data
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

# remove specific count data
#`%ni%` <- Negate(`%in%`)
#remove_list <- row.names(targets[grepl("Olympic", targets$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#targets <- targets[!grepl("Olympic", targets$group),]
#remove_list <- row.names(targets[grepl("Sierra", targets$group),])
#gene_counts <- subset(gene_counts,select = names(gene_counts) %ni% remove_list)
#targets <- targets[!grepl("Sierra", targets$group),]
#gene_counts <- dplyr::select(gene_counts, -contains("PA"))
#targets <- targets[!grepl("PA", targets$group),]

# remove unnecessary factors
targets <- dplyr::select(targets, -contains("group"))
targets <- dplyr::select(targets, -contains("batch"))
targets <- dplyr::select(targets, -contains("tolerance"))


##
# GLM Design
##

# import grouping factor
glm_targets <- targets

# setup a design matrix
glm_group <- factor(paste(glm_targets$genotype, glm_targets$treatment, sep="."))

# begin to construct the DGE list object
glm_list <- DGEList(counts=gene_counts, group=glm_group)

# add the sample names
colnames(glm_list) <- rownames(glm_targets)

# parametrize the experimental design with a one-way layout 
#glm_design <- model.matrix(~ glm_group)
glm_design <- model.matrix(~ 0 + glm_group)
#glm_design <- model.matrix(~ 0 + genotype + treatment, data=glm_targets)

# add group names
colnames(glm_design) <- levels(glm_group)

##
# GLM Normalization
##

# filter the list of gene counts based on expression levels
glm_keep <- filterByExpr(glm_list)

# view the number of filtered genes
table(glm_keep)

# remove genes that are not expressed in either experimental condition
glm_list <- glm_list[glm_keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
glm_list <- calcNormFactors(glm_list)

# compute counts per million (CPM) using normalized library sizes
norm_glm_list <- cpm(glm_list, normalized.lib.sizes=TRUE)

# retrieve the number of grouping levels
stringLevels <- gsub("\\..*","", (levels(glm_group)))
# setup colors and points
colors <- plotColors[1:length(stringLevels)]
#points <- c(0:length(unique(stringLevels)))
png("sample_pca.png")
# add extra space to right of plot area and change clipping to figure
par(mar=c(6.5, 5.5, 5.5, 9.5), xpd=TRUE)
# PCA plot with distances approximating log2 fold changes
#plotMDS(list, col=colors[group], pch=points[group], gene.selection="common", main = "Principal Component Analysis")
plotMDS(glm_list, col=colors[glm_group], gene.selection="common", main = "Principal Component Analysis")
# place the legend outside the right side of the plot
#legend("topright", inset=c(-0.5,0), legend=levels(group), pch=points, col=colors)
legend("topright", inset=c(-0.5,0), legend=levels(glm_group), fill=colors)
dev.off()

##
# GLM Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
glm_list <- estimateDisp(glm_list, glm_design, robust=TRUE)

# estimate the QL dispersions
glm_fit <- glmQLFit(glm_list, glm_design, robust=TRUE)

# plot the QL dispersions
png("QL_dispersions.png")
plotQLDisp(glm_fit)
dev.off()

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
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")

# view the selected color palette
#ghibli_colors

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])

##
# GLM Contrasts
##

###
## condition
###

# examine the overall effect of condition
con.condition <- makeContrasts(set.condition = 
                               (CON_6.UV + E05.UV + E2_1.UV + GRO_3.UV + NGD_1.UV + PA.UV + R2.UV + Sierra.UV + Y002_3_2.UV + Y019_2_1.UV + Y023_5.UV + Y05.UV) -
                               #(E05.UV + E2_1.UV + R2.UV + Y002_3_2.UV + Y019_2_1.UV + Y023_5.UV + Y05.UV) -
                               #(CON_6.UV + GRO_3.UV + NGD_1.UV + Sierra.UV) -
                               #(PA.UV) -
                               #(CON_6.UV + E05.UV + E2_1.UV + GRO_3.UV + NGD_1.UV + R2.UV + Sierra.UV + Y002_3_2.UV + Y019_2_1.UV + Y023_5.UV + Y05.UV) -
                               (CON_6.VIS + E05.VIS + E2_1.VIS + GRO_3.VIS + NGD_1.VIS + PA.VIS + R2.VIS + Sierra.VIS + Y002_3_2.VIS + Y019_2_1.VIS + Y023_5.VIS + Y05.VIS),
                               #(E05.VIS + E2_1.VIS + R2.VIS + Y002_3_2.VIS + Y019_2_1.VIS + Y023_5.VIS + Y05.VIS),
                               #(CON_6.VIS + GRO_3.VIS + NGD_1.VIS + Sierra.VIS),
                               #(PA.VIS),
                               #(CON_6.VIS + E05.VIS + E2_1.VIS + GRO_3.VIS + NGD_1.VIS + R2.VIS + Sierra.VIS + Y002_3_2.VIS + Y019_2_1.VIS + Y023_5.VIS + Y05.VIS),
                               levels = glm_design)

# conduct gene wise statistical tests
anov.condition <- glmTreat(glm_fit, contrast=con.condition)

# view summary of DE genes
summary(decideTests(anov.condition))

# create MD plot of DE genes
png("UV_VIS_MD.png")
plotMD(anov.condition)
# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")
dev.off()

# generate table of DE genes
tagsTbl_condition <- topTags(anov.condition, n=nrow(anov.condition$table), adjust.method="fdr")$table
# add gene row name tag
resultsTbl_condition <- as_tibble(tagsTbl_condition, rownames = "gene")
# output table
write.table(resultsTbl_condition, "UV_VIS_results.csv", sep=",", row.names=FALSE, quote=FALSE)

# add column for identifying direction of DE gene expression
tagsTbl_condition$topDE <- "NA"

# thresholds
cutFDR=0.05
cutLFC=0

# identify significantly up DE genes
tagsTbl_condition$topDE[tagsTbl_condition$logFC > cutLFC & tagsTbl_condition$FDR < cutFDR] <- "UP"

# identify significantly down DE genes
tagsTbl_condition$topDE[tagsTbl_condition$logFC < (-1*cutLFC) & tagsTbl_condition$FDR < cutFDR] <- "DOWN"

# create volcano plot
volcano_condition <- ggplot(data=tagsTbl_condition, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))
# save the plot to a png file
ggsave("UV_VIS_volcano.png", plot = volcano_condition, bg = "white", device = "png", width = 9, height = 8, units = "in")

# identify significantly DE genes by FDR
tagsTbl_condition.glm_keep <- tagsTbl_condition$FDR < cutFDR

# create filtered results table of DE genes
tagsTbl_condition_filtered <- tagsTbl_condition[tagsTbl_condition.glm_keep,]
# add gene row name tag
resultsTbl_condition_filtered <- as_tibble(tagsTbl_condition_filtered, rownames = "gene")
# output table
write.table(resultsTbl_condition_filtered, "UV_VIS_results_sig.csv", sep=",", row.names=FALSE, quote=FALSE)

# identify significantly DE genes
DGESubset_condition <- tagsTbl_condition_filtered#[tagsTbl_condition_filtered$logFC > cutLFC | tagsTbl_condition_filtered$logFC < (-1*cutLFC),]
# subset the log2 CPM by the DGE set
DGESubset_condition.keep <- rownames(norm_glm_list) %in% rownames(DGESubset_condition)
logcpmSubset_condition <- norm_glm_list[DGESubset_condition.keep, ]

# setup annotation names
ann_groups <- data.frame(genotype = as.factor(targets$genotype),
                         treatment = as.factor(targets$treatment),
                         row.names = colnames(logcpmSubset_condition))
# setup annotation colors
ann_colors <- list(
  genotype = plotColors[seq(1, length(unique(targets$genotype)))],
  treatment = c("black", "white"))
#treatment = plotColors[c(length(plotColors)-2,length(plotColors)-1)])

# setup annotation color names
names(ann_colors[[1]]) <- unique(targets$genotyp)
names(ann_colors[[2]]) <- unique(targets$treatment)

# create heatmap for DGE
pheatmap_condition <- as.ggplot(
  pheatmap(logcpmSubset_condition, 
           scale="row", 
           annotation_col = ann_groups, 
           annotation_colors = ann_colors,
           main="Heatmap of DE Genes", 
           show_rownames = FALSE,
           color = colorRampPalette(c(plotColors[5], "white", plotColors[6]))(100))
)
# save the plot to a png file
ggsave("UV_VIS_pheatmap.png", plot = pheatmap_condition, bg = "white", device = "png", width = 9, height = 8, units = "in")
