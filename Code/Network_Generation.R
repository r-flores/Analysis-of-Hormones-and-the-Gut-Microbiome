BiocManager::install("WGCNA")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(htmltools)
library('WGCNA')
library('tidyverse')
library(ggplot2)
library('DESeq2')
options(stringsAsFactors = FALSE);
workDir <- getwd()

Raw_Count <- read.delim(paste(workDir, "/GSE114400_raw_exon_male_liver.txt", sep = ''))
Female_Raw_Count <- read.delim(paste(workDir, "/GSE114400_raw_exon_female.txt", sep = ''))
GH <- read.delim(paste(workDir, "/GSE114610_raw_exon_GH.txt", sep = ''))

convR <- data.frame(select(Raw_Count,1,2,3:14))
GF <- data.frame(select(Raw_Count,1,2,15:26))
Female_convR <- data.frame(select(Female_Raw_Count,1,2,3:14))
Female_GF <- data.frame(select(Female_Raw_Count,1,2,15:26))

###########Data Normalization and cleaning###############
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
##Check For population(Should return TRUE)
gsg = goodSamplesGenes(convR, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(GF, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(GH, verbose = 3);
gsg$allOK


####################################################################################################
####################################################################################################
################Change This value to the variable for which network you want to Create##############
datExpr0 = convR ###################################################################################
####################################################################################################


####CLUSTERING SAMPLES TO DETECT OUTLIERS####
sampleTree1 = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# We detect one outlier in the data for ConvR mice  and will remove it (Sample_CTR2_Liver_ZT10)
# Plot a line to show the cut
abline(h = 3000, col = "red");
# Determine cluster under the line
clust1 = cutreeStatic(sampleTree1, cutHeight = 3000, minSize = 10)
table(clust1)
# clust 1 contains the samples we want to keep.
keepSamples = (clust1==1)

datExpr = datExpr0[keepSamples, ]

###IF No Samples need to be removed###
#datExpr = datExpr0

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#######################
#######################
#End of data cleaning##
#######################
#######################
#################################################Network Creation#################################################
###enable multithreading###
enableWGCNAThreads()

####We select a soft-threshold Power by analzing Network Topology for ConvR_Mice (NETWORK CONSTURCTION FOR convR)####
# Choose a set of soft-thresholding powers
powers = c(c(1:20), seq(from = 20, to=80, by=10))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.6,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##We determine the Power to be 16##
###Co-Expression adjacencies and similarity for ConvRMice###
############################################################################
############################################################################
softpower = 16
adjacency = adjacency(datExpr, power = softpower)

##TOM
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.06);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.06,
                    addGuide = TRUE, guideHang = 0.08,
                    main = "Gene dendrogram and module colors")

###merging expression profiles of similar genes for ConvR##
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.32
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.06,
                    addGuide = TRUE, guideHang = 0.08,
                    main = "Cluster Dendrogram")
#dev.off()
##rename Variables for furthur analysis
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
##########################
##########################
#End of Network Creation##
##########################
##########################

#######Export to Cytoscape###############
# Select modules
modules = c("darkorange", "midnightblue", "blue", "white", "green", "yellow", "grey")
# Select module probes
#####################
##################### replace _sample_genes with the sample currently in use.
#####################
probes = Female_gf_samples_genes$probe
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = Female_gf_samples_genes$gene_symbol[match(modProbes, Female_gf_samples_genes$probe)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("Results/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("Results/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
