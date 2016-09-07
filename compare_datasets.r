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