##########

#Code to compare the differentially expressed genes of mouse and HCT8 over a time course using maSigPro
#Name - Architha Ellappalayam
#Reg No - 940704222120

#Description - For the analysis of single and multi series time course microarray experiments. It uses the 
# regression strategy to find the genes with temporal expression changes and differences between 
# experimental groups. 
############

source("https://bioconductor.org/biocLite.R")
biocLite()

#Package installation 
install.packages("maSigPro")
biocLite("maSigPro")
biocLite("GEOquery")

#loading the libraries 
library(maSigPro)
library(GEOquery)
require(stringr)
require(maSigPro)
require(biomaRt)

#Downloading the GSE29008 datset 
dataset_hct <- getGEO("GSE29008")

#Get the expression values of the dataset 
exp_hct <- exprs(dataset_hct$GSE29008_series_matrix.txt.gz)

#Get the details of the experimental data- the PhenoData 
pd_hct <- pData(dataset_hct$GSE29008_series_matrix.txt.gz)

#Creating the design matrix for the experiment
pd_hct.des <- as.data.frame(str_match(pd_hct$title, "^(.*?)-(.*?)-(.*?)-(.*?)$")[, 2:5]) 
colnames(pd_hct.des) <- c("cell", "agent", "Time", "Replicate")
pd_hct.des$hct <- ifelse(pd_hct.des$cell == "Hct8", 1, 0)
pd_hct.des$Control <- ifelse(pd_hct.des$agent == "Control", 1, 0)
pd_hct.des$TcdA <- ifelse(pd_hct.des$agent == "TcdA", 1, 0)
pd_hct.des$TcdB <- ifelse(pd_hct.des$agent == "TcdB", 1, 0)
pd_hct.des$Time <- gsub("tp", "", pd_hct.des$Time)
pd_hct.des$Time <- gsub("hr", "", pd_hct.des$Time)
pd_hct.des$Time <- as.numeric(pd_hct.des$Time)


#Setting up the rownames for the table 
rownames(pd_hct.des) <- pd_hct$geo_accession

# now we can make the design matrix from the appropriate columns
design_hct <- make.design.matrix(pd_hct.des[, c(3, 4, 6:8)], degree = 3)

exp_hct.des <- exp_hct[, c(1:20)]

#fitting a regression model to discover probesets with significant differential expression over time. 
#The functions p.vector() and T.fit() use print() to report progress, so we’re hiding
#that output here using capture.output().

hide_hct <- capture.output(fit_hct <- p.vector(exp_hct.des, design_hct))
hide_hct <- capture.output(tstep_hct <- T.fit(fit_hct, step.method = "backward", alfa = 0.05))
sigs_hct <- get.siggenes(tstep_hct, rsq = 0.6, vars = "groups")

#Stroing the signififcant p-values for different conditions 
Control_hct <- sigs_hct$sig.genes$Control$sig.pvalues
ToxinA_hct <- sigs_hct$sig.genes$TcdAvsControl$sig.pvalues
ToxinB_hct <- sigs_hct$sig.genes$TcdBvsControl$sig.pvalues

#Matching probesets to genes using BioMart 
getGenes <- function(sig, bm) {
  genes <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"), filters = "affy_hg_u133_plus_2", 
                 values = rownames(sig), mart = bm)
  m <- match(rownames(sig), genes$affy_hg_u133_plus_2)
  sig$gene <- genes[m, "hgnc_symbol"]
  return(sig)
}

mart.hs.hct <- useMart("ensembl", "hsapiens_gene_ensembl")

#Identifying the names of genes using BioMart 
Control_hct <- getGenes(Control_hct, mart.hs.hct)
ToxinA_hct <- getGenes(ToxinA_hct, mart.hs.hct)
ToxinB_hct <- getGenes(ToxinB_hct, mart.hs.hct)

#Creating a function in ggplot2 to view the time course of differentially expressed genes 
plotGenes <- function(e, probe, g, md) {
  require(ggplot2)
  d <- as.data.frame(e[p, ])
  colnames(d) <- "value"
  d$Rep <- md$Replicate
  d$time <- md$Time
  d$agent <- md$agent
  gg <- ggplot(d, aes(time, value)) + geom_boxplot(aes(position = factor(time)), 
                                                   outlier.shape = NA) + theme_bw() + scale_x_continuous(breaks = unique(d$time)) + 
    geom_jitter(aes(color = factor(agent))) + geom_smooth() + labs(title = paste(g, 
                                                                                 probe, sep = "/"), x = "time (hours)", y = "RMA value") + scale_color_discrete(name = "treatment")
  return(gg)
}

#Viewing the data of the top view genes with the significant p-value for Control_hct 
head(Control_hct[order(Control_hct$`p-value`, decreasing = FALSE), ], n = 20)

#plotting the graph for the significant genes in Control_hct ( here, it is gene1(TACSTD))
p <- rownames(Control_hct[order(Control_hct$`p-value`, decreasing = FALSE), ])[7] #The number signifies the gene 
gene <- ifelse(is.na(subset(Control_hct, rownames(Control_hct) == p)$gene), p, subset(Control_hct, rownames(Control_hct) == p)$gene)
plotGenes(exp_hct.des, p, gene, pd_hct.des)

#plotting the graph for the significant gene in Control_hct(NCAPG2)                               
p    <- rownames(Control_hct[order(Control_hct$`p-value`, decreasing = FALSE), ])[20]
gene <- ifelse(is.na(subset(Control_hct, rownames(Control_hct) == p)$gene), p, subset(Control_hct, rownames(Control_hct) == p)$gene)
plotGenes(exp_hct.des, p, gene, pd_hct.des)

#Viewing the top differntial genes of Toxin A with significant p-value 
head(ToxinA_hct[order(ToxinA_hct$`p-value`, decreasing = FALSE),], n= 20)

#Plotting the graph for the gene TACSTD2 (ToxA vs Control)
p    <- rownames(ToxinA_hct[order(ToxinA_hct$p.valor_TcdAvsControl, decreasing = FALSE), ])[3]
gene <- ifelse(is.na(subset(ToxinA_hct, rownames(ToxinA_hct) == p)$gene), p, subset(ToxinA_hct, rownames(ToxinA_hct) == p)$gene)
plotGenes(exp_hct.des, p, gene, pd_hct.des)

#Plotting the graph for the gene JUN (ToxB vs Control)
p    <- rownames(ToxinA_hct[order(ToxinA_hct$p.valor_TcdBvsControl, decreasing = FALSE), ])[3]
gene <- ifelse(is.na(subset(ToxinA_hct, rownames(ToxinA_hct) == p)$gene), p, subset(ToxinA_hct, rownames(ToxinA_hct) == p)$gene)
plotGenes(exp_hct.des, p, gene, pd_hct.des)


#Viewing the top differntial genes of Toxin B with significant p-value 
head(ToxinB_hct[order(ToxinB_hct$`p-value`, decreasing = FALSE), ], n= 20)

#Plotting the graph for the gene TACSTD2 (ToxA vs Control)
p    <- rownames(ToxinB_hct[order(ToxinB_hct$p.valor_TcdAvsControl, decreasing = FALSE), ])[1]
gene <- ifelse(is.na(subset(ToxinB_hct, rownames(ToxinB_hct) == p)$gene), p, subset(ToxinB_hct, rownames(ToxinB_hct) == p)$gene)
plotGenes(exp_hct.des, p, gene, pd_hct.des)


#Visualising the data using Venn diagrams 
venn_diagram <- suma2Venn(sigs_hct$summary[,c(1:3)])

