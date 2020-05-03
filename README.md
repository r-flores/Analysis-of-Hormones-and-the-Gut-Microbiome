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
Check for missing values (Should return True)
```R
gsg = goodSamplesGenes(convR, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(GF, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(GH, verbose = 3);
gsg$allOK
```
#Network Generation
First select which network you wish to create.<br>
[1]Conventionally Rasied Male Mice = convR<br>
[2]Germ Free Rasied Male Mice = GF<br>
[3]Growth Hormone Treated Male Mice = GH<br>
[4]Conventionally Rasied Female Mice = Female_convR<br>
[5]Germ Free Rasied Female Mice = Female_GF<br>
For this walkthrough we will be using Conventionally Rasied Male Mice
```R
datExpr0 = convR
```
Cluster to detect and remove outliers
```R
sampleTree1 = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# We detect one outlier in the data for ConvR mice
# Plot a line to show the cut
abline(h = 3000, col = "red");
# Determine cluster under the line
clust1 = cutreeStatic(sampleTree1, cutHeight = 3000, minSize = 10)
table(clust1)
# clust 1 contains the samples we want to keep.
keepSamples = (clust1==1)

datExpr = datExpr0[keepSamples, ]
```
Uncomment the following if you wish to save all samples.<br>
A count of the genes and samples will be recorded
```R
###IF No Samples need to be removed
#datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
```
The rest of the walkthrough will incorperate multithreading.
```R
enableWGCNAThreads()
```
Determine a soft-threshold, the following will generate a recommendation based on the toplogical analysis
```R
powers = c(c(1:20), seq(from = 20, to=80, by=10))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.6,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
Based on the topological Analysis we chose a soft-power of 16
and generate the adjacency matrix.
```R
softpower = 16
adjacency = adjacency(datExpr, power = softpower)
```
