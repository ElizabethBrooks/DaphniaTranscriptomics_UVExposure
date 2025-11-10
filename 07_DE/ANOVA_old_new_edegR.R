#!/usr/bin/env Rscript

# edgeR user guide
# https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

##
# Working Directory
##

# set the working directory
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/edgeR/treatment_genotype")

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
#colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment.csv", row.names="sample")
#colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design_genotype_treatment_batch.csv", row.names="sample")
colData <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/study_design.csv", row.names="sample")

# verify that the order of the samples in the counts and groupings files match
#colnames(gene_counts)
#rownames(colData)

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
colData <- select(colData, -contains("group"))
colData <- select(colData, -contains("batch"))
colData <- select(colData, -contains("tolerance"))


##
# GLM Design
##

# import grouping factor
glm_targets <- colData

# setup a design matrix
glm_group <- factor(paste(glm_targets$genotype, glm_targets$treatment, sep="."))

# begin to construct the DGE list object
glm_list <- DGEList(counts=gene_counts, group=glm_group)

# add the sample names
colnames(glm_list) <- rownames(glm_targets)

# parametrize the experimental design with a one-way layout 
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

##
# GLM Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
glm_list <- estimateDisp(glm_list, glm_design, robust=TRUE)

# estimate the QL dispersions
glm_fit <- glmQLFit(glm_list, glm_design, robust=TRUE)

# plot the QL dispersions
plotQLDisp(glm_fit)

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
                               (CON_6.VIS + E05.VIS + E2_1.VIS + GRO_3.VIS + NGD_1.VIS + PA.VIS + R2.VIS + Sierra.VIS + Y002_3_2.VIS + Y019_2_1.VIS + Y023_5.VIS + Y05.VIS),
                               levels = glm_design)

# conduct gene wise statistical tests
anov.condition <- glmTreat(glm_fit, contrast=con.condition)

# view summary of DE genes
summary(decideTests(anov.condition))

# create MD plot of DE genes
plotMD(anov.condition)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_condition <- topTags(anov.condition, n=nrow(anov.condition$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_condition$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_condition$topDE[tagsTbl_condition$logFC > 1 & tagsTbl_condition$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_condition$topDE[tagsTbl_condition$logFC < -1 & tagsTbl_condition$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_condition, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_condition.glm_keep <- tagsTbl_condition$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_condition_filtered <- tagsTbl_condition[tagsTbl_condition.glm_keep,]

###
## group
###

# examine the overall effect of group
con.group <- makeContrasts(set.group = 
                           (CON_6.UV + CON_6.VIS + GRO_3.UV + GRO_3.VIS + NGD_1.VIS + NGD_1.UV + Sierra.UV + Sierra.VIS) -
                           (E05.UV + E2_1.UV + E05.VIS + E2_1.VIS + R2.UV + R2.VIS + Y002_3_2.VIS + Y019_2_1.VIS + Y023_5.VIS + Y05.VIS + Y002_3_2.UV + Y019_2_1.UV + Y023_5.UV + Y05.UV),
                           levels = glm_design)

# conduct gene wise statistical tests
anov.group <- glmTreat(glm_fit, contrast=con.group)

# view summary of DE genes
summary(decideTests(anov.group))

# create MD plot of DE genes
plotMD(anov.group)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_group <- topTags(anov.group, n=nrow(anov.group$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_group$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_group$topDE[tagsTbl_group$logFC > 1 & tagsTbl_group$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_group$topDE[tagsTbl_group$logFC < -1 & tagsTbl_group$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_group, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_group.glm_keep <- tagsTbl_group$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_group_filtered <- tagsTbl_group[tagsTbl_group.glm_keep,]


##
# GLM Results Exploration
##

# retrieve set of DE gene names for group contrast
geneSet_group <- rownames(tagsTbl_group_filtered)

# retrieve set of DE gene names for interaction contrast
geneSet_condition <- rownames(tagsTbl_condition_filtered)

# create combined glm_list of DE gene names
glm_list_venn <- list(group = geneSet_group, 
                          interaction = geneSet_condition)

# create venn diagram
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("group","condition")) +
  scale_color_brewer(palette = "Paired")
