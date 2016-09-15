#########

#Name - Architha Ellappalayam 
#Reg no - 9407074222120 
#Code description - Code to load the GSE29008 dataset and normalize the values present 

#########

#Installing the CDF files from BRAINARRAY 
install.packages("~/packages/hgu133plus2hsentrezgcdf_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/hgu133plus2hsentrezg.db_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/hgu133plus2hsentrezgprobe_19.0.0.tar.gz", repos = NULL, type = "source")

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

#Set the working directory 
setwd("/Users/ellap001/Dropbox/thesis/gse_29008/")

#Read the CEL files from the folder 
hct_raw <- ReadAffy(cdfname = "hgu133plus2hsentrezgcdf")

#Performing RMA normalization on the raw data 
hct_rma <- rma(hct_raw)

#Creating the expression matrix for the normalized data 
hct_rma_exprs <- exprs(hct_rma)
colnames(hct_rma_exprs) <-  c("Control_2h_R1","Control_2h_R2","ToxA_2h_R1","ToxA_2h_R2","ToxB_2h_R1","ToxB_2h_R2","Control_6h_R2","Control_6h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1","ToxB_6h_R2","ToxB_6h_R3","Control_24h_R1","Control_24h_R3", "ToxA_24h_R1", "ToxA_24h_R3","ToxB_24h_R1", "Toxb_24h_R3")

#log transforming the values in the matrix for a more normal distribution 
hct_normal <- log2(hct_rma_exprs)

#Printing the expression matrix 
write.csv(hct_normal, "../hct_data.csv")
 

#####Quality control checks for the expression data######

#box plot fo the original data without any RMA normalization
boxplot(hct_raw, col= "chocolate1",main="HCT8- Before normlaization")

#box plot of the normlaized data
boxplot(hct_normal, col="cadetblue1",main="HCT8 RMA expression values")

#Histograms for the data before and after normalisation 
hist_raw <- hist(hct_raw, main = "Histogram before normalisation")
hist_rma <- hist(hct_rma, main = "Histogram after normlaisation")

#creating hierarchical clutsering to see effect of genes in conditions in HCT cells
distance.hct <- dist(t(hct_normal), method="maximum")
clusters.hct <- hclust(distance.hct)
plot(clusters.hct)

#Perform metric calculations on the CEL files 
hct_raw.qc <- fitPLM(hct_raw)

#Using affyPLM to provide more informative boxplots by Relative Log Expression 
#The values here should be close to zero
rle_image <- RLE(hct_raw.qc, main = "RLE", col = "cadetblue1")

#Using NUSE ( Normalised Unscaled Standard Errors)
#The median standard error should be 1 for most genes
nuse_image <- NUSE(hct_raw.qc, main = "NUSE", col = "brown1")


#####Finding differentially expressed genes for hct data ########
samples.hct <- c("Control_2h","Control_2h","ToxA_2h","ToxA_2h","ToxB_2h","ToxB_2h","Control_6h","Control_6h","ToxA_6h","ToxA_6h","ToxA_6h","ToxB_6h","ToxB_6h","ToxB_6h","Control_24h","Control_24h", "ToxA_24h", "ToxA_24h","ToxB_24h", "ToxB_24h")
samples.hct <- as.factor(samples.hct)

#The levels in the design matrix are checked by looking at samples 
design.hct <- model.matrix(~0 + samples.hct)
colnames(design.hct) <- c("Control_24h", "Control_2h", "Control_6h", "ToxA_24h", "ToxA_2h","ToxA_6h","ToxB_24h", "ToxB_2h", "ToxB_6h")

#providing the data for limma analysis

#fit the linear model to the expression set
fit.hct <- lmFit(hct_normal ,design.hct)

#setting up the contrast matrix for the design matrix
contrast.matrix.hct <- makeContrasts(ToxA_2h_ToxB_6h = ToxA_2h - ToxB_6h, ToxA_2h_ToxB_24h = ToxA_2h - ToxB_24h, ToxA_6h_ToxB_24h <- ToxA_6h - ToxB_24h, levels=design.hct)

#contrast matrix is now combined with the per-probeset linear model fit
hct_fits <- contrasts.fit(fit.hct, contrast.matrix.hct)
hct_ebFit <- eBayes(hct_fits)

#return the top results of the contrast fit
genelist.hct<- toptable(hct_ebFit, coef= 1, number = 20000) ##coef = 1( 2 hours), coef = 2(6 hours), coef = 3(24 hours)###

#Condition to check for differentially expressed genes with p-value less than 0.01
sum(genelist.hct$adj.P.Val < 0.01)

#Converting the log fold changes to regular Fold Changes 
genelist.hct$FC <- 2^genelist.hct$logFC

#write the data onto a separate file 
write.csv(genelist.hct, "genelist.csv")
