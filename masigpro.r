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
dataset <- getGEO("GSE29008")

#Get the expression values of the dataset 
exp <- exprs(dataset$GSE29008_series_matrix.txt.gz)

#Get the details of the experimental data- the PhenoData 
pd <- pData(dataset$GSE29008_series_matrix.txt.gz)

#Creating the design matrix for the experiment
pd.des <- as.data.frame(str_match(pd$title, "^(.*?)-(.*?)-(.*?)-(.*?)$")[, 2:5]) 
colnames(pd.des) <- c("cell", "agent", "Time", "Replicate")
pd.des$hct <- ifelse(pd.des$cell == "Hct8", 1, 0)
pd.des$Control <- ifelse(pd.des$agent == "Control", 1, 0)
pd.des$TcdA <- ifelse(pd.des$agent == "TcdA", 1, 0)
pd.des$TcdB <- ifelse(pd.des$agent == "TcdB", 1, 0)
pd.des$Time <- gsub("tp", "", pd.des$Time)
pd.des$Time <- gsub("hr", "", pd.des$Time)
pd.des$Time <- as.numeric(pd.des$Time)


#Setting up the rownames for the table 
rownames(pd.des) <- pd$geo_accession

# now we can make the design matrix from the appropriate columns
design <- make.design.matrix(pd.des[, c(3, 4, 6:8)], degree = 3)

exp.des <- exp[, c(1:20)]

#fitting a regression model to discover probesets with significant differential expression over time. 
#The functions p.vector() and T.fit() use print() to report progress, so we're hiding
#that output here using capture.output().

hide <- capture.output(fit <- p.vector(exp.des, design))
hide <- capture.output(tstep <- T.fit(fit, step.method = "backward", alfa = 0.05))
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

#Stroing the signififcant p-values for different conditions 
control <- sigs$sig.genes$Control$sig.pvalues
ToxinA <- sigs$sig.genes$TcdAvsControl$sig.pvalues
ToxinB <- sigs$sig.genes$TcdBvsControl$sig.pvalues

#Matching probesets to genes using BioMart 
getGenes <- function(sig, bm) {
  genes <- getBM(attributes = c("affy_hg_u133a_2", "hgnc_symbol"), filters = "affy_hg_u133a_2", 
                 values = rownames(sig), mart = bm)
  m <- match(rownames(sig), genes$affy_hg_u133a_2)
  sig$gene <- genes[m, "hgnc_symbol"]
  return(sig)
}

mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")

#Identifying the names of genes using BioMart 
control <- getGenes(control, mart.hs)
ToxinA <- getGenes(ToxinA, mart.hs)
ToxinB <- getGenes(ToxinB, mart.hs)

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

#Viewing the data of the top view genes with the significant p-value for COntrol 
head(control[order(control$`p-value`, decreasing = FALSE), ])

#plotting the graph for the significant genes in Control ( here, it is gene1(TACSTD))
p <- rownames(control[order(control$`p-value`, decreasing = FALSE), ])[1] #The number signifies the gene 
gene <- ifelse(is.na(subset(control, rownames(control) == p)$gene), p, subset(control, rownames(control) == p)$gene)
plotGenes(exp.des, p, gene, pd.des)

#plotting the graph for the significant gene in Control(NCAPG2)                               
p    <- rownames(control[order(control$`p-value`, decreasing = FALSE), ])[6]
gene <- ifelse(is.na(subset(control, rownames(control) == p)$gene), p, subset(control, rownames(control) == p)$gene)
plotGenes(exp.des, p, gene, pd.des)

#Viewing the top differntial genes of Toxin A with significant p-value 
head(ToxinA[order(ToxinA$`p-value`, decreasing = FALSE), ])

#Plotting the graph for the gene TACSTD2 (ToxA vs Control)
p    <- rownames(ToxinA[order(ToxinA$p.valor_TcdAvsControl, decreasing = FALSE), ])[3]
gene <- ifelse(is.na(subset(ToxinA, rownames(ToxinA) == p)$gene), p, subset(ToxinA, rownames(ToxinA) == p)$gene)
plotGenes(exp.des, p, gene, pd.des)

#Plotting the graph for the gene JUN (ToxB vs Control)
p    <- rownames(ToxinA[order(ToxinA$p.valor_TcdBvsControl, decreasing = FALSE), ])[3]
gene <- ifelse(is.na(subset(ToxinA, rownames(ToxinA) == p)$gene), p, subset(ToxinA, rownames(ToxinA) == p)$gene)
plotGenes(exp.des, p, gene, pd.des)


#Viewing the top differntial genes of Toxin B with significant p-value 
head(ToxinB[order(ToxinB$`p-value`, decreasing = FALSE), ])

#Plotting the graph for the gene TACSTD2 (ToxA vs Control)
p    <- rownames(ToxinB[order(ToxinB$p.valor_TcdAvsControl, decreasing = FALSE), ])[1]
gene <- ifelse(is.na(subset(ToxinB, rownames(ToxinB) == p)$gene), p, subset(ToxinB, rownames(ToxinB) == p)$gene)
plotGenes(exp.des, p, gene, pd.des)


#Visualising the data using Venn diagrams 
venn_diagram <- suma2Venn(sigs$summary[,c(1:3)])


