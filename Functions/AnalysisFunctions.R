###
###Include all required packages
###
# Install BioConductor package manager
if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
BiocManager::install()
# General plotting requirements
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
library("ggplot2")
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}
library('scales')
if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel")}
library("ggrepel")
if("heatmap3" %in% rownames(installed.packages()) == FALSE) {install.packages('heatmap3')}
library("heatmap3")
if('RColorBrewer' %in% rownames(installed.packages()) == FALSE) {install.packages('RColorBrewer')}
library('RColorBrewer')
if('viridis' %in% rownames(installed.packages()) == FALSE) {install.packages('viridis')}
library('viridis')
if('randomcoloR' %in% rownames(installed.packages()) == FALSE) {install.packages('randomcoloR')}
library('randomcoloR')
# General data manipulation packages
if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages('reshape2')}
library("reshape2")
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages('stringr')}
library("stringr")
if("openxlsx" %in% rownames(installed.packages()) == FALSE) {install.packages('openxlsx')}
library("openxlsx")
# Required packages for regulon analysis
if("cluster" %in% rownames(installed.packages()) == FALSE) {install.packages('cluster')}
library("cluster")
#if('factoextra' %in% rownames(installed.packages()) == FALSE) {install.packages('factoextra')}
#library('factoextra')
# Required packages for gene networks
if('sna' %in% rownames(installed.packages()) == FALSE) {install.packages('sna')}
library('sna')
if('network' %in% rownames(installed.packages()) == FALSE) {install.packages('network')}
library('network')
if('ggnetwork' %in% rownames(installed.packages()) == FALSE) {install.packages('ggnetwork')}
library('ggnetwork')
# Required for gene ontology
if("GOstats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install('GOstats')}
library("GOstats")
if("org.Mm.eg.db" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("org.Mm.eg.db")}
library("org.Mm.eg.db")
if("org.Hs.eg.db" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("org.Hs.eg.db")}
library("org.Hs.eg.db")
if("GSEABase" %in% rownames(installed.packages()) == FALSE) {BiocManager::install('GSEABase')}
library("GSEABase")

###
###Include all supporting functions
###
source(paste0(scriptsDirectory,"/Functions/generateColumnNames.R"))
source(paste0(scriptsDirectory,"/Functions/loadData.R"))
source(paste0(scriptsDirectory,"/Functions/applyExclusions.R"))
source(paste0(scriptsDirectory,"/Functions/loadMetaData.R"))

###
###Include all analysis functions
###
source(paste0(scriptsDirectory,"/Functions/rawDataQualityControl.R"))
source(paste0(scriptsDirectory,"/Functions/cutoffControl.R"))
source(paste0(scriptsDirectory,"/Functions/geneQuality.R"))
source(paste0(scriptsDirectory,"/Functions/readDistribution.R"))
source(paste0(scriptsDirectory,"/Functions/buildSCSObject.R"))
source(paste0(scriptsDirectory,"/Functions/processedDataQualityControl.R"))
source(paste0(scriptsDirectory,"/Functions/analysePlateContribution.R"))
source(paste0(scriptsDirectory,"/Functions/analyseClusterDistribution.R"))
source(paste0(scriptsDirectory,"/Functions/analyseClusterComposition.R"))
source(paste0(scriptsDirectory,"/Functions/analyseClusterDifferentialExpression.R"))
source(paste0(scriptsDirectory,"/Functions/analyseMarkerGenesPerCluster.R"))
source(paste0(scriptsDirectory,"/Functions/plotGeneExpressionTSNE.R"))
source(paste0(scriptsDirectory,"/Functions/plotGeneExpressionCutoffTSNE.R"))
source(paste0(scriptsDirectory,"/Functions/plotMetaDataTSNE.R"))
source(paste0(scriptsDirectory,"/Functions/plotClusterTSNE.R"))
source(paste0(scriptsDirectory,"/Functions/plotClusterHeatmap.R"))
source(paste0(scriptsDirectory,"/Functions/analyseCorrelationAB.R"))
source(paste0(scriptsDirectory,"/Functions/analyseCorrelationMetaDataGene.R"))
source(paste0(scriptsDirectory,"/Functions/analyseCorrelation.R"))
source(paste0(scriptsDirectory,"/Functions/analyseMetaDataCorrelation.R"))
source(paste0(scriptsDirectory,"/Functions/plotCorrelationVolcano.R"))
source(paste0(scriptsDirectory,"/Functions/compareClustersAcrossGroups.R"))
source(paste0(scriptsDirectory,"/Functions/analyseClusterTyping.R"))
source(paste0(scriptsDirectory,"/Functions/returnTopDifferentialGenesPerCluster.R"))
source(paste0(scriptsDirectory,"/Functions/analyseClusterGeneOntology.R"))
source(paste0(scriptsDirectory,"/Functions/plotClusterGeneOntology.R"))
source(paste0(scriptsDirectory,"/Functions/analyseCorrelationGeneOntology.R"))
source(paste0(scriptsDirectory,"/Functions/plotCorrelationGeneOntology.R"))
source(paste0(scriptsDirectory,"/Functions/generateCorrelationMatrix.R"))
source(paste0(scriptsDirectory,"/Functions/analyzeRegulons.R"))
source(paste0(scriptsDirectory,"/Functions/plotRegulons.R"))
source(paste0(scriptsDirectory,"/Functions/analyzeRegulonGeneOntology.R"))
source(paste0(scriptsDirectory,"/Functions/plotRegulonGeneOntology.R"))
