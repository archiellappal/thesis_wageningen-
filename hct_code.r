##########

#Code to process the CDF data files, edition 2
#Name - Architha Ellappalayam
#Reg No - 940704222120

###########

#Set the working directory
setwd("/Users/ellap001/Desktop/thesis/gse_29008/")

#Installing the packages
install.packages("statmod")
install.packages("arrayQualityMetrics")

#installing the packages from BrainArray 
install.packages("~/packages/hgu133plus2hsentrezgcdf_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/hgu133plus2hsentrezg.db_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/hgu133plus2hsentrezgprobe_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/AnnBuilder_1.16.0.tar.gz", repos = NULL, type = "source")

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
biocLite("GOSummaries")


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
library(mouse4302mmentrezgcdf)
library(mouse4302mmentrezg.db)
library(mouse4302mmentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(hgu133plus2hsentrezgprobe)
library(AnnBuilder)
library(GOsummaries)
#Read the data files
hct_eight <- ReadAffy(cdfname = "hgu133plus2hsentrezgcdf" )

#Processing the data in RMA menthod with background correction
hct_eight_rma <- rma(hct_eight)

#creating the expression matrix(genes in rows and condiitons in columns)
hct_eight_normal = exprs(hct_eight_rma)

#Setting the column names in the matrix
colnames(hct_eight_normal) = c("Control_2h_R1","Control_2h_R2","ToxA_2h_R1","ToxA_2h_R2","ToxB_2h_R1","ToxB_2h_R2","Control_6h_R2","Control_6h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1","ToxB_6h_R2","ToxB_6h_R3","Control_24h_R1","Control_24h_R3", "ToxA_24h_R1", "ToxA_24h_R3","ToxB_24h_R1", "Toxb_24h_R3")

#log transforming the expression values for normal distribution
hct_eight_normal = log(hct_eight_normal, 2)

#####Quality control checks for the expression data######

#box plot fo the original data without any RMA normalization
boxplot(hct_eight, col= "chocolate1",main="HCT8 Probe intensities")

#box plot of the normlaized data
boxplot(hct_eight_normal, col="cadetblue1",main="HCT8 RMA expression values")

#creating hierarchical clutsering to see effect of genes in conditions in HCT cells
distance <- dist(t(hct_eight_normal), method="maximum")
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
contrast.matrix.hct <- makeContrasts(ToxA_2h_ToxB_2h = ToxA_2h - ToxB_2h, ToxA_6h_ToxB_6h = ToxA_6h - ToxB_6h, ToxA_24h_ToxB_24h <- ToxA_24h - ToxB_24h, levels=design.hct)

#contrast matrix is now combined with the per-probeset linear model fit
hct_fits <- contrasts.fit(fit.hct, contrast.matrix.hct)
hct_ebFit <- eBayes(hct_fits)

#write the differentially expressed genes to separate tables 
probeset_one <- topTable(hct_ebFit, coef = 1, number = 100000)
probeset_two <- topTable(hct_ebFit, coef = 2, number = 100000)
probeset_three <- topTable(hct_ebFit, coef = 3, number = 100000)

