setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
install.packages("BiocManager")
BiocManager::install("WGCNA")
library('tidyverse')
library('WGCNA')

Raw_Count <- read.delim("C:/Users/Ricky.Flores/Desktop/WS/GSE114400_raw_exon_male_liver.txt")
convR <- data.frame(select(Raw_Count,1,2,3:14))
GF <- data.frame(select(Raw_Count,1,2,15:26))

#Data normaliztion
## Take Average of every row for both experiments
AvgconvR <- data.frame(gene=convR[,1], ExonLength = convR[,2], select(convR,3:14), Average=rowMeans(convR[3:14]))
AvgGF <- data.frame(gene=GF[,1], ExonLength = GF[,2], select(convR,3:14), Average=rowMeans(GF[3:14]))

## Remove values with an average count less than 500
convR_clean <- AvgconvR[!(AvgconvR$Average < 500),]
GF_clean <- AvgGF[!(AvgGF$Average < 500),]

##Expression Value
convRExpression <- data.frame(gene=convR_clean[,1], convR_clean[,3:14] / convR_clean[,2])
GFExpression <- data.frame(gene=GF_clean[,1], GF_clean[,3:14] / GF_clean[,2])