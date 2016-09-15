#########

#Name - Architha Ellappalayam 
#Reg no - 9407074222120 
#Code description - Code to load the GSE244091 dataset and normalize the values present 

#########

#Installing the CDF files from BRAINARRAY 
install.packages("~/packages/mouse4302mmentrezgcdf_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/mouse4302mmentrezg.db_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/mouse4302mmentrezgprobe_19.0.0.tar.gz", repos = NULL, type = "source")

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
setwd("/Users/ellap001/Dropbox/thesis/gse_44091/")

#Read the CEL files from the folder 
mouse_raw <- ReadAffy(cdfname = "mouse4302mmentrezgcdf")

#Performing RMA normalization on the raw data 
mouse_rma <- rma(mouse_raw)

#Creating the expression matrix for the normalized data 
mouse_rma_exprs <- exprs(mouse_rma)
colnames(mouse_rma_exprs) <- c("ToxA_2h_R1","ToxA_2h_R2","ToxA_2h_R3","ToxB_2h_R1","ToxB_2h_R2","ToxB_2h_R3","ToxAB_2h_R1","ToxAB_2h_R2","ToxAB_2h_R3","Sham_2h_R1","Sham_2h_R2","Sham_2h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1", "ToxB_6h_R2", "ToxB_6h_R3","ToxAB_6h_R1", "Sham_6h_R1", "Sham_6h_R2", "Sham_6h_R3","ToxA_16h_R1","ToxA_16h_R2", "ToxA_16h_R3","ToxB_16h_R1", "ToxB_16h_R2", "ToxB_16h_R3","Sham_16h_R1","Sham_16h_R2", "Sham_16h_R3","Sham_16h_R4")

#log transforming the values in the matrix for a more normal distribution 
mouse_normal <- log2(mouse_rma_exprs)

#Printing the expression matrix 
write.csv(mouse_normal, "../hct_data.csv")


#####Quality control checks for the expression data######

#box plot fo the original data without any RMA normalization
boxplot(mouse_raw, col= "chocolate1",main="HCT8- Before normlaization")

#box plot of the normlaized data
boxplot(mouse_normal, col="cadetblue1",main="HCT8 RMA expression values")

#Histograms for the data before and after normalisation 
hist_raw_mouse <- hist(hct_raw, main = "Histogram before normalisation")
hist_rma_mouse <- hist(hct_rma, main = "Histogram after normlaisation")

#creating hierarchical clutsering to see effect of genes in conditions in HCT cells
distance.mouse <- dist(t(mouse_normal), method="maximum")
clusters.mouse <- hclust(distance.mouse)
plot(clusters.mouse)

#Perform metric calculations on the CEL files 
mouse_raw.qc <- fitPLM(mouse_raw)

#Using affyPLM to provide more informative boxplots by Relative Log Expression 
#The values here should be close to zero
rle_image_mouse <- RLE(mouse_raw.qc, main = "RLE", col = "cadetblue1")

#Using NUSE ( Normalised Unscaled Standard Errors)
#The median standard error should be 1 for most genes
nuse_image_mouse <- NUSE(mouse_raw.qc, main = "NUSE", col = "brown1")

#After examining the NUSE image, the 24th value seems out of proportion and hence removed 
mouse_normal <- mouse_normal[, -24]

#####Finding differentially expressed genes for hct data ########
samples.mouse <- c("ToxA_2h_R1","ToxA_2h_R2","ToxA_2h_R3","ToxB_2h_R1","ToxB_2h_R2","ToxB_2h_R3","ToxAB_2h_R1","ToxAB_2h_R2","ToxAB_2h_R3","Sham_2h_R1","Sham_2h_R2","Sham_2h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1", "ToxB_6h_R2", "ToxB_6h_R3","ToxAB_6h_R1", "Sham_6h_R1", "Sham_6h_R2", "Sham_6h_R3","ToxA_16h_R1","ToxA_16h_R2", "ToxA_16h_R3","ToxB_16h_R1", "ToxB_16h_R2", "ToxB_16h_R3","Sham_16h_R1","Sham_16h_R2", "Sham_16h_R3","Sham_16h_R4")
samples.mouse <- as.factor(samples.mouse)

#The levels in the design matrix are checked by looking at samples 
design.mouse <- model.matrix(~0 + samples.mouse)
colnames(design.mouse) <- c("Sham_16h", "Sham_2h", "Sham_6h", "ToxAB_2h", "ToxAB_6h", "ToxA_16h", "ToxA_2h", "ToxA_6h", "ToxB_16h", "ToxB_2h", "ToxAB_6h")

#providing the data for limma analysis

#fit the linear model to the expression set
fit.mouse <- lmFit(mouse_normal ,design.mouse)

#setting up the contrast matrix for the design matrix
contrast.matrix.mouse <- makeContrasts(ToxA_2h_ToxB_2h = ToxA_2h - ToxB_2h,ToxA_16h_ToxB_16h = ToxA_16h - ToxB_16h,  levels = design.mouse)

#contrast matrix is now combined with the per-probeset linear model fit
mouse_fits <- contrasts.fit(fit.mouse, contrast.matrix.mouse)
mouse_ebFit <- eBayes(mouse_fits)

#return the top results of the contrast fit
genelist.mouse<- toptable(mouse_ebFit, coef= 2, number = 20000) ##coef = 1( 2 hours), coef = 2(6 hours), coef = 3(24 hours)###

#Condition to check for differentially expressed genes with p-value less than 0.01
sum(genelist.mouse$adj.P.Val < 0.01)

#Converting the log fold changes to regular Fold Changes 
genelist.mouse$FC <- 2^genelist.mouse$logFC

#write the data onto a separate file 
write.csv(genelist.mouse, "genelist_mouse.csv")
