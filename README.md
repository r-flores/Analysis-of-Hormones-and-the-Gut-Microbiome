# Analysis-of-Hormones-and-the-Gut-Microbiome
Repo for BIOI 4890 Senior Project Spring 2020
Clone the repository or download the Zip file,
this folder will conatin all the necessary data
and R script which when run will generate the network
that can be imported to Cytoscape for Analysis.
This tutorial is for the generation of the conventually raised mouse gene co-expression network.<br>
<br>
Install the following packages if required.
```R
BiocManager::install("WGCNA")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")
```
Set the working directory to file location, import the necessary libraries, and set options
```R
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(htmltools)
library('WGCNA')
library('tidyverse')
library(ggplot2)
library('DESeq2')
options(stringsAsFactors = FALSE);
workDir <- getwd()
```
Read in the Data
```R
Raw_Count <- read.delim(paste(workDir, "/GSE114400_raw_exon_male_liver.txt", sep = ''))
Female_Raw_Count <- read.delim(paste(workDir, "/GSE114400_raw_exon_female.txt", sep = ''))
GH <- read.delim(paste(workDir, "/GSE114610_raw_exon_GH.txt", sep = ''))
```
Seperate each sample type
```R
convR <- data.frame(select(Raw_Count,1,2,3:14))
GF <- data.frame(select(Raw_Count,1,2,15:26))
Female_convR <- data.frame(select(Female_Raw_Count,1,2,3:14))
Female_GF <- data.frame(select(Female_Raw_Count,1,2,15:26))
```
Data cleaning and normalization
```R
## Take Average of every row for both experiments
convR <- data.frame(gene=convR[,1], ExonLength = convR[,2], select(convR,3:14), Average=rowMeans(convR[3:14]))
GF <- data.frame(gene=GF[,1], ExonLength = GF[,2], select(GF,3:14), Average=rowMeans(GF[3:14]))
GH <- data.frame(gene=GH[,1], ExonLength = GH[,2], select(GH,3:13), Average=rowMeans(GH[3:13]))

Female_convR <- data.frame(gene=Female_convR[,1], ExonLength = Female_convR[,2], select(Female_convR,3:14), Average=rowMeans(Female_convR[3:14]))
Female_GF <- data.frame(gene=Female_GF[,1], ExonLength = Female_GF[,2], select(Female_GF,3:14), Average=rowMeans(Female_GF[3:14]))

## Remove values with an average count less than 500
convR <- convR[!(convR$Average < 500),]
GF <- GF[!(GF$Average < 500),]
GH <- GH[!(GH$Average < 500),]
Female_convR <- Female_convR[!(Female_convR$Average < 500),]
Female_GF <- Female_GF[!(Female_GF$Average < 500),]

##Expression Value
convR <- data.frame(gene=convR[,1], convR[,3:14] / convR[,2])
GF <- data.frame(gene=GF[,1], GF[,3:14] / GF[,2])
GH <- data.frame(gene=GH[,1], GH[,3:13] / GH[,2])
Female_convR <- data.frame(gene=Female_convR[,1], Female_convR[,3:14] / Female_convR[,2])
Female_GF <- data.frame(gene=Female_GF[,1], Female_GF[,3:14] / Female_GF[,2])

##Data formatting

convR <- separate(convR, "gene", c("probe", "gene_symbol"), sep = '\\|')
GF <- separate(GF, "gene", c("probe", "gene_symbol"), sep = '\\|')
GH <- separate(GH, "gene", c("probe", "gene_symbol"), sep = '\\|')
Female_convR <- separate(Female_convR, "gene", c("probe", "gene_symbol"), sep = '\\|')
Female_GF <- separate(Female_GF, "gene", c("probe", "gene_symbol"), sep = '\\|')

convR_samples_genes <- convR[c(1,2)]
gf_samples_genes <- GF[c(1,2)]
GH_samples_genes <- GH[c(1,2)]
Female_convR_samples_genes <- Female_convR[c(1,2)]
Female_gf_samples_genes <- Female_GF[c(1,2)]

rownames(convR) <- convR[,1]
rownames(GF) <- GF[,1]
rownames(GH) <- GH[,1]
rownames(Female_convR) <- Female_convR[,1]
rownames(Female_GF) <- Female_GF[,1]

convR[,1] <- NULL
convR[,1] <- NULL
GF[,1] <- NULL
GF[,1] <- NULL
GH[,1] <- NULL
GH[,1] <- NULL
Female_convR[,1] <- NULL
Female_convR[,1] <- NULL
Female_GF[,1] <- NULL
Female_GF[,1] <- NULL
##Transpose the data for analysis
convR <- t(convR)
GF <- t(GF)
GH <- t(GH)
Female_convR <- t(Female_convR)
Female_GF <- t(Female_GF)
```
