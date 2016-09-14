##########

#Code to process the CDF data files, edition 2
#Name - Architha Ellappalayam
#Reg No - 940704222120

###########

#Set the working directory
setwd("/Users/ellap001/Desktop/thesis/gse_44091//")

#Installing the packages
install.packages("statmod")
install.packages("arrayQualityMetrics")
install.packages("compare")

#installing the packages from BrainArray 
install.packages("~/packages/mouse4302mmentrezgcdf_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/mouse4302mmentrezg.db_19.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/packages/mouse4302mmentrezgprobe_19.0.0.tar.gz", repos = NULL, type = "source")

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
library(compare)

#Read the data files
mouse_cells <- ReadAffy(cdfname = "mouse4302mmentrezgcdf")

#Processing the data in RMA menthod with background correction
mouse_cells_rma <- rma(mouse_cells)

#creating the expression matrix(genes in rows and conditons in columns)
mouse_cells_normal = exprs(mouse_cells_rma)

#Setting the column names in the matrix
colnames(mouse_cells_normal) = c("ToxA_2h_R1","ToxA_2h_R2","ToxA_2h_R3","ToxB_2h_R1","ToxB_2h_R2","ToxB_2h_R3","ToxAB_2h_R1","ToxAB_2h_R2","ToxAB_2h_R3","Sham_2h_R1","Sham_2h_R2","Sham_2h_R3","ToxA_6h_R1","ToxA_6h_R2","ToxA_6h_R3","ToxB_6h_R1", "ToxB_6h_R2", "ToxB_6h_R3","ToxAB_6h_R1", "Sham_6h_R1", "Sham_6h_R2", "Sham_6h_R3","ToxA_16h_R1","ToxA_16h_R2", "ToxA_16h_R3","ToxB_16h_R1", "ToxB_16h_R2", "ToxB_16h_R3","Sham_16h_R1","Sham_16h_R2", "Sham_16h_R3","Sham_16h_R4")

#log transforming the expression values for normal distribution
mouse_cells_normal = log(mouse_cells_normal, 2)

#bind these data to the originial table 
complete_probeset <- cbind(mouse_twohours_toxA,mouse_twohours_toxB, mouse_sixhours_toxA, mouse_sixhours_toxB,mouse_teenhours_toxA, mouse_teenhours_toxB)

#creating hierarchical clutsering to see effect of genes in conditions in mouse cells
distance <- dist(t(mouse_cells_normal), method="maximum") 
clusters <- hclust(distance)
plot(clusters)

#####Finding ifferentially expressed genes for the mouse data#######
samples.mouse <- c("ToxA_2h","ToxA_2h","ToxA_2h","ToxB_2h","ToxB_2h","ToxB_2h","ToxAB_2h","ToxAB_2h","ToxAB_2h","Sham_2h","Sham_2h","Sham_2h","ToxA_6h","ToxA_6h","ToxA_6h","ToxB_6h", "ToxB_6h", "ToxB_6h","ToxAB_6h", "Sham_6h", "Sham_6h", "Sham_6h","ToxA_16h","ToxA_16h", "ToxA_16h","ToxB_16h", "ToxB_16h", "ToxB_16h","Sham_16h","Sham_16h", "Sham_16h","Sham_16h")
samples.mouse <- as.factor(samples.mouse)
design.mouse <- model.matrix(~0 + samples.mouse)
colnames(design.mouse) <- c("Sham_16h", "Sham_2h", "Sham_6h","ToxA_16h", "ToxA_2h", "ToxA_6h","ToxAB_2h", "ToxAB_6h",  "ToxB_16h", "ToxB_2h", "ToxB_6h")

#providing the data for limma analysis

#fit the linear model to the expression set
fit.mouse <- lmFit(mouse_cells_normal ,design.mouse)

#setting up the contrast matrix for the design matrix
contrast.matrix.mouse <- makeContrasts(ToxA_6h_ToxA_6h <- ToxA_6h - ToxA_6h, ToxA_16h_ToxA_16h = ToxA_16h - ToxA_16h,  levels = design.mouse)

#contrast matrix is now combined with the per-probeset linear model fit
mouse_fits <- contrasts.fit(fit.mouse, contrast.matrix.mouse)
mouse_ebFit <- eBayes(mouse_fits)

# create two lists for differentially expressed genes at 6 and sixteen hours 
sixteen <- topTable(mouse_ebFit, coef = 2, number= 2500)
six  <- topTable(mouse_ebFit, coef = 1, number = 2000)

#Create a combined list of differentially expressed genes 
de_genes <- unique(append((row.names(six)), row.names(sixteen)))


