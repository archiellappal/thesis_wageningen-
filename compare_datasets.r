##########

#Code to process the CDF data files, edition 2
#Name - Architha Ellappalayam
#Reg No - 940704222120

###########

#Clear the workspace in the global environment
rm(list = ls())

#Set the working directory
setwd("/Users/ellap001/Documents/thesis/data/gse_29008/")

#Installing the packages
install.packages("statmod")
install.packages("arrayQualityMetrics")

#Package for exploring oligonucleotide array analysis
source("https://bioconductor.org/biocLite.R")
biocLite()

#Downloading the packages from bioclite
biocLite("affy")
biocLite("simpleaffy")
biocLite("limma")
biocLite("genefilter")
biocLite("annotate")
biocLite("arrayQualityMetrics")
biocLite("GenomicRanges")
biocLite("hgu133plus2hsentrezgcdf")
biocLite("hgu133plus2.db")
biocLite("hgu133plus2hsentrezgcdf")
biocLite("mouse4302.db")
biocLite("mouse4302cdf")
biocLite("mouse4302probe")

#Loading the libraries in the environment
library(affy)
library(simpleaffy)
library(affyPLM)
library(limma)
library(genefilter)
library(annotate)
library(GenomicRanges)          #package to perfrom differential expression analysis
library(statmod)
library(RColorBrewer)         #load the color libraries
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(mouse4302.db)
library(mouse4302cdf)
library(mouse4302probe)

#Read the data files
hct_eight <- ReadAffy(cdfname = "hgu133plus2cdf" )

#Processing the data in RMA menthod with background correction
hct_eight_rma <- rma(hct_eight)

#creating the expression matrix(genes in rows and condiitons in columns)
hct_eight_normal = exprs(hct_eight_rma)

#Setting the column names in the matrix
colnames(hct_eight_normal) = c("Control_2h_R1","Control_2h_R2","ToxA_2h_R1","ToxA_2h_R2","ToxB_2h_R1","ToxB_2h_R2","Control_6h_R2","Control_6h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1","ToxB_6h_R2","ToxB_6h_R3","Control_24h_R1","Control_24h_R3", "ToxA_24h_R1", "ToxA_24h_R3","ToxB_24h_R1", "Toxb_24h_R3")

#log transforming the expression values for normal distribution
hct_eight_normal = log(hct_eight_normal, 2)

#Set the working directory
setwd("../gse_44091/")

#Read the data files
mouse_cells <- ReadAffy(cdfname = "mouse4302cdf")

#Processing the data in RMA menthod with background correction
mouse_cells_rma <- rma(mouse_cells)

#creating the expression matrix(genes in rows and conditons in columns)
mouse_cells_normal = exprs(mouse_cells_rma)

#Setting the column names in the matrix
colnames(mouse_cells_normal) = c("ToxA_2h_R1","ToxA_2h_R2","ToxA_2h_R3","ToxB_2h_R1","ToxB_2h_R2","ToxB_2h_R3","ToxAB_2h_R1","ToxAB_2h_R2","ToxAB_2h_R3","Sham_2h_R1","Sham_2h_R2","Sham_2h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1", "ToxB_6h_R2", "ToxB_6h_R3","ToxAB_6h_R1", "Sham_6h_R1", "Sham_6h_R2", "Sham_6h_R3","ToxA_16h_R1","ToxA_16h_R2", "ToxA_16h_R3","ToxB_16h_R1", "ToxB_16h_R2", "ToxB_16h_R3","Sham_16h_R1","Sham_16h_R2", "Sham_16h_R3","Sham_16h_R4")

#log transforming the expression values for normal distribution
mouse_cells_normal = log(mouse_cells_normal, 2)


#####Quality control checks for the expression data######

#box plot fo the original data without any RMA normalization
boxplot(hct_eight, col= "chocolate1",main="HCT8 Probe intensities")

#box plot of the normlaized data
boxplot(hct_eight_normal, col="cadetblue1",main="HCT8 RMA expression values")

#creating hierarchical clutsering to see effect of genes in conditions in HCT cells
distance <- dist(t(hct_eight_normal), method="maximum")
clusters <- hclust(distance)
plot(clusters)

#creating hierarchical clutsering to see effect of genes in conditions in mouse cells
distance <- dist(t(mouse_cells_normal), method="maximum")
clusters <- hclust(distance)
plot(clusters)

#####Finding differentially expressed genes for hct data ########
samples.hct <- c("Control_2h","Control_2h","ToxA_2h","ToxA_2h","ToxB_2h","ToxB_2h","Control_6h","Control_6h","ToxA_6h","ToxA_6h","ToxA_6h","ToxB_6h","ToxB_6h","ToxB_6h","Control_24h","Control_24h", "ToxA_24h", "ToxA_24h","ToxB_24h", "ToxB_24h")
samples.hct <- as.factor(samples.hct)
design.hct <- model.matrix(~0 + samples.hct)
colnames(design.hct) <- c("Control_24h", "Control_2h", "Control_6h", "ToxA_24h", "ToxA_2h","ToxA_6h","ToxB_24h", "ToxB_2h", "ToxB_6h")

#providing the data for limma analysis

#fit the linear model to the expression set
fit.hct <- lmFit(hct_eight_normal ,design.hct)

#setting up the contrast matrix for the design matrix
contrast.matrix.hct <- makeContrasts(ToxA_24h_ToxB_2h <- ToxA_24h - ToxB_2h, ToxA_24h_ToxB_6h <- ToxA_24h - ToxB_6h, ToxA_24h_ToxB_24h <- ToxA_24h - ToxB_24h, levels = design.hct)

#contrast matrix is now combined with the per-probeset linear model fit
hct_fits <- contrasts.fit(fit.hct, contrast.matrix.hct)
hct_ebFit <- eBayes(hct_fits)

#return the top results of the contrast fit
probeset.list.hct <- toptable(hct_ebFit, coef=1, number = 100) ##coef = 1( 2 hours), coef = 2(6 hours), coef = 3(24 hours)###

#####Finding ifferentially expressed genes for the mouse data#######
samples.mouse <- c("ToxA_2h","ToxA_2h","ToxA_2h","ToxB_2h","ToxB_2h","ToxB_2h","ToxAB_2h","ToxAB_2h","ToxAB_2h","Sham_2h","Sham_2h","Sham_2h","ToxA_6h","ToxA_6h","ToxA_6h","ToxB_6h", "ToxB_6h", "ToxB_6h","ToxAB_6h", "Sham_6h", "Sham_6h", "Sham_6h","ToxA_16h","ToxA_16h", "ToxA_16h","ToxB_16h", "ToxB_16h", "ToxB_16h","Sham_16h","Sham_16h", "Sham_16h","Sham_16h")
samples.mouse <- as.factor(samples.mouse)
design.mouse <- model.matrix(~0 + samples.mouse)
colnames(design.mouse) <- c("Sham_16h", "Sham_2h", "Sham_6h", "ToxAB_2h", "ToxAB_6h", "ToxA_16h", "ToxA_2h", "ToxA_6h", "ToxB_16h", "ToxB_2h", "ToxAB_6h")

#providing the data for limma analysis

#fit the linear model to the expression set
fit.mouse <- lmFit(mouse_cells_normal ,design.mouse)

#setting up the contrast matrix for the design matrix
contrast.matrix.mouse <- makeContrasts(ToxA_2h_ToxB_2h = ToxA_2h - ToxB_2h,ToxA_16h_ToxB_16h = ToxA_16h - ToxB_16h,  levels = design.mouse)

#contrast matrix is now combined with the per-probeset linear model fit
mouse_fits <- contrasts.fit(fit.mouse, contrast.matrix.mouse)
mouse_ebFit <- eBayes(mouse_fits)

#return the top results of the contrast fit for ToxA and Tox B at 2 hours
probeset.list.mouse <- toptable(mouse_ebFit, coef=2, number = 5000 )  ### coef = 1( 2hours), coef = 2(6 hours)###
