install.packages("BiocManager")
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

Raw_Count <- read.delim("C:/Users/Ricky/Desktop/Hormone and Gut Microbiome/Data/GSE114400_raw_exon_male_liver.txt")

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

##make sure no gene values are empty
convRExpression_Values_Only <- data.frame(convRExpression[,2:13])
gsg = goodSamplesGenes(test, verbose = 3);
gsg$allOK

GFExpression_Values_Only <- data.frame(GFExpression[,2:13])
gsg = goodSamplesGenes(test, verbose = 3);
gsg$allOK

####multiExperiment####
nSets = 2
setLabels = c("Conventually_Raised", "Germ-Free_Raised")
shortLabels = c("convR", "GF")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(convRExpression[c(2:13)])));
names(multiExpr[[1]]$data) = convRExpression$gene;
rownames(multiExpr[[1]]$data) = names(convRExpression)[c(2:13)];
multiExpr[[2]] = list(data = as.data.frame(t(GFExpression[c(2:13)])));
names(multiExpr[[2]]$data) = GFExpression$gene;
rownames(multiExpr[[2]]$data) = names(GFExpression)[c(2:13)];


####plot a tree to detect SAMPLE outliers####
sampleTree = hclust(dist(convRExpression), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut H is where we want to cut off
abline(h = 6000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = test[keepSamples, ]
nGenes = ncol(convRExpression)
nSamples = nrow(convRExpression)

####Nerwork Construction####

options(stringsAsFactors = FALSE)
# Choose a set of soft-thresholding powers
powers = c(c(1:20), seq(from = 20, to=100, by=10))
# Call the network topology analysis function
sft = pickSoftThreshold(convRExpression_Values_Only, powerVector = powers, verbose = 5)
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
abline(h=-0.13,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##12 is the lowest highest value on scale of independance###
net = blockwiseModules(convRExpression_Values_Only, power = 12,
                       TOMType = "unsigned", minModuleSize = 1,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ConvrMiceTOM",
                       verbose = 3)
table(net$colors)
###Clustering Dendrogram###
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "ConVRMice-02-networkConstruction-auto.RData")

####Visualizations####
nGenes = ncol(convRExpression)
nSamples = ncol(convRExpression)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(convRExpression_Values_Only, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


####Export####

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(convRExpression_Values_Only, power = 12);
# Select module
module = "turquoise";
# Select module probes
probes = names(convRExpression_Values_Only)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(convRExpression$gene))


#####TestingCode#####
##Adjancency###
convRExpression_Values_Only <- data.frame(convRExpression[,2:13])
softPower = 1
adjacency = adjacency(convRExpression_Values_Only, power = softPower)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10;
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
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
