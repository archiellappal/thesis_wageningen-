#########

#Name - Architha Ellappalayam 
#Reg no - 9407074222120 
#Code description - Code to load the CaCo2 dataset and normalize the values present 

#########

#Installing the CDF files from BRAINARRAY 
install.packages("C:/Users/ellap001/Dropbox/thesis/brain_array_packages/hugene11sthsentrezg.db_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/ellap001/Dropbox/thesis/brain_array_packages/hugene11sthsentrezgcdf_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/ellap001/Dropbox/thesis/brain_array_packages/hugene11sthsentrezgprobe_19.0.0.tar.gz", repos = NULL, type = "source")

#Packages for microarray analysis are to be downloaded from BiocLite 
source("https://bioconductor.org/biocLite.R")
biocLite()

#Downloading the required packages
biocLite("affy")
biocLite("limma")
biocLite("annotate")
biocLite("GEOquery")
biocLite("affyPLM")
biocLite("genefilter")


#Loading the libraries in the environment 
library(affy)
library(limma)
library(annotate)
library(GEOquery)
library(affyPLM)
library(genefilter)
library(hugene11sthsentrezg.db)
library(hugene11sthsentrezgcdf)
library(hugene11sthsentrezgprobe)

#Set the working directory 
setwd("/Users/ellap001/Dropbox/thesis/data/caco2_files/")

#Read the CEL files from the folder 
caco_raw <- ReadAffy(cdfname = "hugene11sthsentrezgcdf")

#Performing RMA normalization on the raw data 
caco_rma <- rma(caco_raw)

#Creating the expression matrix for the normalized data 
caco_rma_exprs <- exprs(caco_rma)
colnames(caco_rma_exprs) <- c("Control_7d_R1","ToxA_7d_R1", "ToxB_7d_R1", "ToxAB_7d_R1", "Control_7d_R2", "ToxA_7d_R2", "ToxB_7d_R2", "ToxAB_7d_R2", "Control_7d_R3","ToxA_7d_R3", "ToxB_7d_R3", "ToxAB_7d_R3", "Control_21d_R1","ToxA_21d_R1", "ToxB_21d_R1", "ToxAB_21d_R1", "Control_21d_R2","ToxA_21d_R2", "ToxB_21d_R2", "ToxAB_21d_R2", "Control_21d_R3","ToxA_21d_R3", "ToxB_21d_R3", "ToxAB_21d_R3") 

#log transforming the values in the matrix for a more normal distribution 
caco_normal <- log2(caco_rma_exprs)

#Printing the expression matrix 
write.csv(caco_normal, "../../code_output/caco_data.csv")


#####Quality control checks for the expression data######

#box plot fo the original data without any RMA normalization
boxplot(caco_raw, col= "violetred1",main="HCT8- Before normlaization")

#box plot of the normlaized data
boxplot(caco_normal, col="thistle1",main="HCT8 RMA expression values")

#Histograms for the data before and after normalisation 
hist_raw_caco <- hist(caco_raw, main = "Histogram before normalisation")
hist_rma_caco <- hist(caco_rma, main = "Histogram after normalisation")

#creating hierarchical clutsering to see effect of genes in conditions in HCT cells
distance.caco <- dist(t(caco_normal), method="maximum")
clusters.caco <- hclust(distance.caco)
plot(clusters.caco)

#Perform metric calculations on the CEL files 
caco_raw.qc <- fitPLM(caco_raw)

#Using affyPLM to provide more informative boxplots by Relative Log Expression 
#The values here should be close to zero
rle_image_caco <- RLE(caco_raw.qc, main = "RLE", col = "blueviolet")

#Using NUSE ( Normalised Unscaled Standard Errors)
#The median standard error should be 1 for most genes
nuse_image_caco <- NUSE(caco_raw.qc, main = "NUSE", col = "yellow")

#####Finding differentially expressed genes for caco2 data ########
samples.caco <- c("Control_7d","ToxA_7d", "ToxB_7d", "ToxAB_7d", "Control_7d", "ToxA_7d", "ToxB_7d", "ToxAB_7d", "Control_7d","ToxA_7d", "ToxB_7d", "ToxAB_7d", "Control_21d","ToxA_21d", "ToxB_21d", "ToxAB_21d", "Control_21d","ToxA_21d","ToxB_21d", "ToxAB_21d", "Control_21d","ToxA_21d", "ToxB_21d", "ToxAB_21d") 
samples.caco <- as.factor(samples.caco)

#The levels in the design matrix are checked by looking at samples 
design.caco <- model.matrix(~0 + samples.caco)
colnames(design.caco) <- c("Control_21d", "Control_7d", "ToxA_21d", "ToxA_7d", "ToxAB_21d", "ToxAB_7d", "ToxB_21d", "ToxB_7d")
  
#providing the data for limma analysis

#fit the linear model to the expression set
fit.caco <- lmFit(caco_normal ,design.caco)

#setting up the contrast matrix for the design matrix
contrast.matrix.caco <- makeContrasts(ToxA_7d_ToxB_7d = ToxA_7d - ToxB_7d,ToxA_21d_ToxB_21d = ToxA_21d - ToxB_21d,  levels = design.caco)

#contrast matrix is now combined with the per-probeset linear model fit
caco_fits <- contrasts.fit(fit.caco, contrast.matrix.caco)
caco_ebFit <- eBayes(caco_fits)

#return the top results of the contrast fit
genelist.caco <- toptable(caco_ebFit, coef= 2, number = 20000) ##coef = 1( 2 hours), coef = 2(6 hours), coef = 3(24 hours)###

#Condition to check for differentially expressed genes with p-value less than 0.01
sum(genelist.caco$adj.P.Val < 0.01)

#Converting the log fold changes to regular Fold Changes 
genelist.caco$FC <- 2^genelist.caco$logFC

#write the data onto a separate file 
write.csv(genelist.caco, "../../code_output/genelist_caco.csv")






