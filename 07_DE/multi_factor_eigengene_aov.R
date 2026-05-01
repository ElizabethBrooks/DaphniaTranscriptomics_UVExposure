#!/usr/bin/env Rscript

# usage: Rscript ANOVA_OLYMGenotypes_aov.r workingDir countsFile factorGroupingFile
# usage Ex: Rscript ANOVA_OLYMGenotypes_aov.r /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_genotypes_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression_line.csv /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVGenotypes/InputData/expDesign_OlympicsGenotypes.csv
# R script to perform statistical analysis of gene count tables using aov
# note: https://www.r-bloggers.com/2022/05/two-way-anova-example-in-r-quick-guide/
# https://conjugateprior.org/2013/01/formulae-in-r-anova/
# https://environmentalcomputing.net/statistics/linear-models/anova/anova-nested/

# turn off scientific notation
options(scipen=999)

# install packages
#install.packages("ggpubr")
#install.packages("multcomp")

# load libraries
library(ggpubr)
library(multcomp)
library(rcartocolor)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[11], plotColors[6], plotColors[4], plotColors[5])

# retrieve input file name of gene counts
args = commandArgs(trailingOnly=TRUE)

# set working directory
setwd("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DA/edgeR_LFC0.01_FDR0.05_normalized/genotype_noE2/NA/DA/DESeq2_LFC0.1_FDR0.05_normalized/treatment_group_batch_TG_noE2")

# import expression data
inputTable <- read.csv("/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/DA/edgeR_LFC0.01_FDR0.05_normalized/genotype_noE2/NA/eigengeneExpression.csv", row.names="gene")

# transpose expression data
inputTable <- as.data.frame(t(inputTable))

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/PfrenderLab/melanica_UV_exposure/old_new_merged/data/study_design_noE2.csv", row.names="sample")

# retrieve input design samples
designSamples <- data.frame(ID1 = rownames(targets))
# retrieve input gene counts samples
countsSamples <- data.frame(ID2 = rownames(inputTable))
# first simply check if the number of samples matches
if(nrow(designSamples) != nrow(countsSamples)) return(NULL)
# find samples in counts, but not in design
mismatch_counts <- designSamples %>% 
  filter(!designSamples$ID1 %in% countsSamples$ID2)
# find samples in design, but not in counts
mismatch_design <- countsSamples %>% 
  filter(!countsSamples$ID2 %in% designSamples$ID1)
# check total non matches
totalMismatches <- nrow(mismatch_counts) + nrow(mismatch_design)

# setup data frame
expData <- merge(targets, inputTable, by = 'row.names') 

# set row names
rownames(expData) <- expData$Row.names

# remove the Row.names column
expData <- expData[,-1]

# retrieve columns to factor
targets_cols <- colnames(targets)

# loop over each factor column
for (col in 1:length(targets_cols)) {
  expData[,targets_cols[col]] <- factor(expData[,targets_cols[col]])
}

# check the structure of the expression data
#str(expData)

# make frequency tables
#table(expData$treatment, expData$genotype, expData$tolerance, expData$replicate)

# retrieve module column name
#modName <- colnames(expData)[6]

# update module column name
#colnames(expData)[colnames(expData) == modName] <- "expression"

# view data
#print(expData)

# input expression
in_exp <- "treatment + group + batch + treatment:group"
#in_exp <- "batch - group"

# initialize data
start <- length(targets_cols)+1
stop <- length(expData)-1
aov_exp <- NA

# loop over each sample
for (col in start:stop) {
  # get sample name
  samp <- colnames(expData)[col]
  # setup expression
  aov_exp <- paste(samp, "~", in_exp, sep=" ")
  # compute anova
  affy.aov <- aov(as.formula(aov_exp), data = expData)
  # retrieve summary stats
  stats <- summary(affy.aov)
  # separate summary stats
  df <- stats[[1]][1]
  sum_sq <- stats[[1]][2]
  mean_sq <- stats[[1]][3]
  f_val <- stats[[1]][4]
  pr <- data.frame(Pr = stats[[1]][5][1:nrow(df)-1,])
  # add sample name to columns
  colnames(df) <- samp
  colnames(sum_sq) <- samp
  colnames(mean_sq) <- samp
  colnames(f_val) <- samp
  colnames(pr) <- samp
  # fix row names
  rownames(pr) <- rownames(df)[1-nrow(df)-1]
  # initialize data frames
  if (col == start) {
    df_table <- data.frame(Df = df)
    sum_sq_table <- data.frame(SumSq = sum_sq)
    mean_sq_table <- data.frame(MeanSq = mean_sq)
    f_val_table <- data.frame(FValue = f_val)
    pr_table <- data.frame(Pr = pr)
  } else {
    # add current sample data to each results table
    df_table <- cbind(df_table, df)
    sum_sq_table <- cbind(sum_sq_table, sum_sq)
    mean_sq_table <- cbind(mean_sq_table, mean_sq)
    f_val_table <- cbind(f_val_table, f_val)
    pr_table <- cbind(pr_table, pr)
  }
}




# create a colored box plot
exportFile <- paste(modName, "coloredBoxPlot.png", sep="_")
png(file=exportFile)
ggboxplot(data=expData, x="treatment", y="expression", color="group",
          palette = plotColorSubset)
dev.off()

# box plot with two variable factors
#png(file=exportFile)
#boxplot(expression ~ treatment * (tolerance/genotype), data=expData, frame=FALSE,
#        col = plotColorSubset)
#dev.off()

# two way interaction plot
exportFile <- paste(modName, "twoWayInteractionPlot.png", sep="_")
png(file=exportFile)
interaction.plot(x.factor = expData$treatment, 
                trace.factor = expData$group,
                response = expData$expression, 
                fun = mean,
                type = "b", 
                legend = TRUE,
                xlab = "Treatment", 
                ylab="Expression",
                pch=c(1,1,1), 
                col = plotColorSubset
                )
dev.off()


## check the validity of ANOVA assumptions
# The data must be regularly distributed, and the variation between groups must be homogeneous

## examine the assumption of homogeneity of variance
# the residuals versus fits graphic are used to assess for variance homogeneity

# plot homogeneity of variances
exportFile <- paste(modName, "homogeneityPlot.png", sep="_")
png(file=exportFile)
plot(affy.aov, 1)
dev.off()

# examine the assumption of normality
# in a residuals’ normality plot the residuals quantiles are displayed against the normal distribution quantiles
# the residuals’ normal probability plot is used to confirm that the residuals are normally distributed
# the residuals’ normal probability plot should roughly follow a straight line

# normality plot
exportFile <- paste(modName, "normalityPlot.png", sep="_")
png(file=exportFile)
plot(affy.aov, 2)
dev.off()

# extract the residuals
aovRes <- residuals(object = affy.aov)

# run Shapiro-Wilk test
swTest <- shapiro.test(x = aovRes)

# write test statistics to a file
exportFile <- paste(modName, "shapiroTest_summary.csv", sep="_")
capture.output(swTest, file=exportFile)
