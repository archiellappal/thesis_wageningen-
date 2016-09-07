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
