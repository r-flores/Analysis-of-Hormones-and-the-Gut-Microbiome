# Analysis-of-Hormones-and-the-Gut-Microbiome
Repo for BIOI 4890 Senior Project Spring 2020
Clone the repository or download the Zip file,
this folder will conatin all the necessary data
and R script which when run will generate the network
that can be imported to Cytoscape for Analysis.
This tutorial is for the generation of the conventually raised mouse gene co-expression network.<br>
Install the following packages if required.
```R
BiocManager::install("WGCNA")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")
```
