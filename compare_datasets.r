##########

#Code to compare the differentially expressed genes of CaCo2, mouse and HCT8 using XGSA 
#Name - Architha Ellappalayam
#Reg No - 940704222120

#Description - XGSA is a statistical method for cross species gene set analysis, to 
# determine if a given set of differentially expressed genes are statistically enriched for 
# genes in other pre-defined sets. This method avoids false positive discoveries while maintaining 
# good statistical power. 

############

#Set the working directory
setwd("/Users/ellap001/Dropbox/thesis/code_output/")

#Package for exploring oligonucleotide array analysis
source("https://bioconductor.org/biocLite.R")
biocLite()

#Installing the packages for graphs 
install.packages("igraph")
install.packages("devtools")
install.packages("slam", "igraph")

#Downloading the packages from bioclite
biocLite("biomaRt")
biocLite("AnnotationDbi")
biocLite("GO.db")
biocLite("graph")

#Downlaoding the XGSA tool from GITHUB
devtools::install_github('VCCRI/XGSA')

#Loading the libraries in the environment
library(biomaRt)
library(AnnotationDbi)
library(GO.db)
library(graph)
library(igraph)
library(devtools)
library(xgsa)
library(annotate)
library(mouse4302mmentrezg.db)  

####### Setting up the mouse data##############

#Extract the mouse_genes value 
mouse_genes <- read.table("genelist_mouse.csv", header = TRUE ,sep = ",")

#Renaming the title
mouse_genes$gene_id <- mouse_genes$X

#Deleting the other extra columns and reatinging only the gene ids 
mouse_genes <- mouse_genes[, -c(1,2,3,4,5,6)]
mouse_genes$FC <- NULL

#creating a mouse table with their Entrez Id's 
PROBE <- as.character(mouse_genes$gene_id)

#Extracting their Entrez ID and removing their probe ID
mouse_table <- select(mouse4302mmentrezg.db,PROBE , c("ENTREZID"))
mouse_table$PROBEID <- NULL

#Extracting the whole list of Entrez IDs from Mus Musculus and removing the external Gene Name 
mouse.ensembl.symbol.map <- get_ENSEMBL_symbol_map(species = 'mmusculus')
mouse.ensembl.symbol.map$external_gene_name <- NULL

#Adding another additonal column with gene ids for easier identification 
mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

#Converting the entire list of ENSEMBLE ids to gene IDS
ensemble_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene"),values= mouse.ensembl.symbol.map,mart= mart)

#Finding the ENSEMBLE ID for all the gene_ids present in the mouse_table 
mouse.cecum.ensembl.symbols <- ensemble_table$ensembl_gene_id[ensemble_table$entrezgene %in% mouse_table$ENTREZID]

#Creating a proper XGSA dataset with the ENSEMBLE IDS of the Mus Musculus data 
mouse.data <- new_XGSA_dataset(species = 'mmusculus', data = list(mouseCecumGenes = mouse.cecum.ensembl.symbols), type = 'genesetlist', name = 'MouseCecumGenes', universe = unique(ensemble_table$ensembl_gene_id))

##########Setting up the human data###############

#In the gene universe, we use all the genes that are present in homo sapiens for our analysis
human.GO <- get_GO('hsapiens', ontologies = "biological_process")
#> [1] "retrieved GO"
human.GO <- human.GO[lapply(human.GO, length) > 10 & lapply(human.GO, length) < 500]
human.GO.data <- new_XGSA_dataset(species = "hsapiens", data = human.GO, type = 'genesetlist', name = "humanGO", universe = unique(unlist(human.GO)))

#Comparing the mouse_genes against the whole set of human chromosomes 
mouse_vs_humans <- run_XGSA_test(mouse.data, human.GO.data)

# We need to separate the pvalues and the overlapping gene IDs, because XGSA returns both.
resulting.pvals <- lapply(mouse_vs_humans, function(X){ X[["pvals"]] })
resulting.overlap.genes <- lapply(mouse_vs_humans, function(X){ X[["genes"]] })

# The Benjamini Hochberg multiple hypothesis testing correction to the pvalues is performed 
adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")

#The names of the results are made interpretable
names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ,"\\."), function(X){return(X[[2]])}))

# We can use another XGSA helper function to find out the GO term names.
human.GO.names <- get_GO_names('hsapiens')

# And finally we get interpretable names
names(adjusted.pvals) <- human.GO.names[match( names(adjusted.pvals), human.GO.names$go_id),"name_1006"]
#

significant.GO.Terms <- adjusted.pvals
print(head(sort(significant.GO.Terms),10))

par(mar=c(10,5,4,2))
barplot(-log10(head(sort(significant.GO.Terms),10)), ylab = "- log 10 p-value", las=2)

