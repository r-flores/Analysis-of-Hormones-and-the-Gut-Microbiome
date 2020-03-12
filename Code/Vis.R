install.packages("BiocManager")
BiocManager::install("WGCNA")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

library(htmltools)
library('tidyverse')
library(ggplot2)
library('DESeq2')

GSE114400_raw_exon_male_liver <- read.delim("C:/Users/Ricky/Desktop/Hormone and Gut Microbiome/Data/GSE114400_raw_exon_male_liver.txt")


convR_geneEX <- data.frame(select(GSE114400_raw_exon_male_liver,1,3:14))
# Filter Values that are less than 10
convR_geneEX_10 <- convR_geneEX[apply(convR_geneEX[,-1], 1, function(x) !all(x<=10)),]

GF_geneEX <- data.frame(select(GSE114400_raw_exon_male_liver,1,15:26))
# Filter Values that are less than 10
GF_geneEX_10 <- GF_geneEX[apply(GF_geneEX[,-1], 1, function(x) !all(x<=10)),]

write.csv(convR_geneEX_10, "convR.csv")
write.csv(convR_geneEX_10, "GF.csv")

sampleTree = hclust(dist(convR_geneEX_10), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Clustering to detect outlier", sub = '', xlab = '', cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#Network Construction
