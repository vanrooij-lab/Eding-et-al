# Load the scripts that add in the functions
# Determine output folder
{
  #Set scriptsDirectory
  scriptsDirectory = '/Users/j.eding/Work/Scripts/SingleCellSeq/'
  
  #Load required script files. #TODO: make finding these automatically easier (based on enforcing the folder structure?)
  source(paste0(scriptsDirectory,"/Functions/RaceID2_StemID_class.R"))
  source(paste0(scriptsDirectory,"/Functions/AnalysisFunctions.R"))
  
  #Generate output dir, include date & time to prevent overwriting. 
  outputDir <- paste0('analysis_', format(Sys.time(), "%Y-%m-%d-%Hh_%Mm"))
}

# Make default configuration
{
  #Default options, here for convenience
  options(stringsAsFactors = FALSE)

  #Set working directory (this directory contains all data files)
  setwd("/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/AllDataFiles")
  
  #Initialize a config variable. This variable will hold _all_ configuration for the rest of the analysis.
  config <- list()
    
  #Add scriptsDirectory to config so it's available within function too
  config$scriptsDirectory = scriptsDirectory
  
  #Make a list of samples and which plates should be in which sample
  ###Example: config$samples[['sampleName']] <- c('namePlate1', 'namePlate2', namePlateX')
  ###In the example, the filename for namePlate1 is namePlate1_TranscriptCounts.tsv
  config$samples <- list()
  config$samples[['patient1']] <- c('JE1', 'JE2', 'JE3', 'JE4')
  config$samples[['patient2']] <- c('JE5', 'JE6', 'JE7', 'JE8')
  config$samples[['patient3']] <- c('MW5', 'MW6', 'MW7', 'MW8')
  config$samples[['patient4']] <- c('JE10', 'JE11')
  
  #Make a list of groups and define by which samples they are made up.
  config$groups <- list()
  config$groups[['patientAll']] <- c('patient1', 'patient2', 'patient3')
  config$groups[['patientAllWithIndex']] <- c('patient1', 'patient2', 'patient3', 'patient4')
  config$groups[['patient1']] <- c('patient1')
  config$groups[['patient2']] <- c('patient2')
  config$groups[['patient3']] <- c('patient3')
  config$groups[['patient4']] <- c('patient4')

  
  #Make a list of colors, some function will use colors matching to a group or sample in graphing
  config$colors <- list()
  config$colors[['patient1']] <- "#67AD5D"
  config$colors[['patient2']] <- "#E88B3F"
  config$colors[['patient3']] <- "#D76A9D"
  config$colors[['patient4']] <- "#CCCCCC"
  
  #Make a list of text_elements used for ggPlot
  #TODO: Consider removing this (is it really used?)
  config$plotOptions <- list()
  config$plotOptions[['title']] <- element_text()
  config$plotOptions[['subtitle']] <-element_text(size=7, face="italic", color="black")
  config$plotOptions[['axis.title']] <- element_text(size=10)
  config$plotOptions[['axis.text']] <- element_text(size=8)
  config$plotOptions[['legend.title']] <- element_text(size=rel(0.75))
  config$plotOptions[['legend.text']] <- element_text(size=rel(0.5))
  config$plotOptions[['axis.text.x']] <- element_text(angle = 90, hjust = 1, vjust=1)
  config$plotOptions[['axis.text.y']] <- element_text()
  config$plotOptions[['aspect.ratio']] <- 1
  config$plotOptions[['panel.border']] <- element_rect(colour = "black", fill=NA, size=1)
  config$plotOptions[['axis.line']] <- element_blank()
  
  # Enter gene identifier type, pick one from:
  #   - "ensembl_gene_id"               (identifier looks like 'ENSMUSG00000000001')
  #   - "external_gene_name"            (identifier looks like 'Gnai3')
  #   - "external_gene_name_with_chr"   (identifier looks like 'Gnai3__chr1')
  #   - "entrezgene"                    (identifier looks like '14679')
  config$geneIdentifierType <- "external_gene_name_with_chr"
  
  # Specify species, pick one from:
  #   - "human"
  #   - "mouse"
  config$species <- "human"
  
  #Exclude mitochondrial genes? Set to T(RUE) or F(ALSE)
  config$excludeMitochondrialGenes = T
  
  ###
  ### Configure data selection and clustering
  ###
  #Define minimum total number of unique reads per cell (default = 1000, check read distribution over cells to pick a definite value)
  config$minTotalReads = 1000
  #Define p-value at which a cell is considered an outlier within its own cluster (default = 1e-7)
  config$outlierPValue = 1e-7
  
  #Define lists of cell-type markers
  ###Case-sensitive (pay attention with human vs murine capitalization)
  ###Matches specifically to only the full gene name you put in. So: 'MYH7' matches with 'MYH7' but not with 'MYH7B'
  config$markerGenes <- list(
    'cardiomyocytes' = c('MYH6','MYH7','MYL2','MYL3', 'ATP2A2', 'RYR2', 'TTN', 'ACTN2', 'ACTC1', 'TNNT2', 'MYBPC3'),
    'fibroblast' = c('COL1A2', 'COL3A1', 'FSTL1', 'DCN', 'COL6A1', 'FBN1', 'POSTN', 'MMP2', 'VIM'),
    'vascular' = c('PECAM1', 'SDPR', 'CAV1','VEGFA', 'ACTA2', 'IGFBP5', 'EGFL7', 'LY6E','AQP1', 'ENG'),
    'macrophage' = c('FTL','LYZ','SPP2','MPEG1','CD68','CSTB','LGALS3','APOE','MRC2','CCL2')
  )
  config$markerGenes <- list(
    'cardiomyocytes' = c('MYH6','MYH7', 'TTN', 'TNNT2', 'MYBPC3', 'NPPA'),
    'fibroblast' = c('COL1A2', 'COL3A1', 'FSTL1', 'POSTN', 'VIM'),
    'vascular' = c('PECAM1', 'CAV1','VEGFA'),
    'macrophage' = c('FTL','LYZ', 'CD68', 'LGALS3','APOE')
  )
  if(!config$excludeMitochondrialGenes) {
    config$markerGenes <- c(config$markerGenes, list('mitochondrial' = c('MT-')))
  }
  config$markerGenes <- append(config$markerGenes$cardiomyocytes, append(config$markerGenes$fibroblast, append(config$markerGenes$vascular,config$markerGenes$macrophage)))
  config$markerGenes <- config$markerGenes[length(config$markerGenes):1]
}

# Load data
groupedData <- loadData(config, c('patientAll', 'patientAllWithIndex', 'patient1', 'patient2', 'patient3', 'patient4'))

# Apply exclusions
groupedData <- applyExclusions(config, groupedData)

# Perform quality control on the raw data
rawDataQualityControl(config, groupedData, groupsToAnalyse = NULL, 'save')
  
# Build SCSObjects for the groupedData
groupedSCS <- buildSCSObject(config, groupedData)

# Perform quality control on the processed data
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')

# Look at clustering
plotClusterTSNE(config, groupedSCS, groupNames=NULL, outputMode='pdf')
analyseClusterComposition(config, groupedSCS, groupNames = NULL, outputMode = 'pdf')
analyseClusterDistribution(config, groupedSCS, groupNames = NULL, outputMode = 'save')
plotClusterHeatmap(config,groupedSCS,groupNames=NULL,includeOutliers=T,hmethod='single',outputMode='pdf')
clusterDifferentialExpression <- analyseClusterDifferentialExpression(config, groupedSCS, groupNames = NULL, outputMode = 'save')

# Analyse marker gene expression per cluster
markerGeneNames = c('^TTN','MYH7_','VIM_','COL3A1','PECAM1','CAV1','FTL','CSTB','GAPDH')
markerGeneNames = c('MYH7_','NPPA','MTRNR2L1_','TCAP_','KCNQ1OT1','MALAT1')
markerGeneNames = c('MYH7_', 'TNNT1', 'TNNT2', '^TNNI3_', 'ACTN2', 'RYR2', 'ATP2A2')
analyseMarkerGenesPerCluster(config, clusterDifferentialExpression, geneNames=markerGeneNames, groupNames='patientAllMod', outputMode='pdf')
analyseMarkerGenesPerCluster(config, clusterDifferentialExpression, geneNames=markerGeneNames, groupNames=NULL, outputMode='save')
analyseClusterTyping(config, groupedSCS, clusterDifferentialExpression, geneNames = NULL, groupNames = 'patientAllMod', numGenesPerCluster = 1, normalize='no', displayLog = 'yes', overrideDir=F, outputMode = 'show')

# Compare clusters accross groups
compareClustersAcrossGroups(config, groupedSCS, referenceGroup='patientAll', groupNames=NULL, outputMode = 'pdf')

# Draw t-SNE maps for a set of interesting genes
geneNames <- c('^TTN_','CMYA5','NPPA','NPPB','MYH7_','VIM_','COL3A1','PECAM1','CAV1','FTL','CSTB','GAPDH','ZFP106')
geneNames <- c('GATA4', 'GNAS_','GNAS-', 'PAX2','PAX7','TBX5_','TBX5-','^DSP_','NKX2-5','LMX1B','^EN2_') #From TRIAGE
plotGeneExpressionTSNE(config, groupedSCS, geneNames=geneNames, groupNames = 'patientAll', includeOutlierCells = T, colorScheme='Spectral', expressionScale='linear', outputMode='pdf')
plotGeneExpressionTSNE(config, groupedSCS, geneNames=geneNames, groupNames = NULL, includeOutlierCells = T, colorScheme='Spectral', expressionScale='linear', outputMode='pdf')

# Analyse correlation for one gene with all others (Volcano)
correlationResultTTN = analyseCorrelation(config, groupedSCS, geneName='^TTN_', groupName='patientAll', removeNoExpressCells = 'yes', percentage=20)
plotCorrelationVolcano(config, correlationResult=correlationResultTTN, coefficientCutoff = c(-0.25,0.25), pValueCutoff=1/10^5, outputMode='pdf')
correlationResultKCNQ1OT1 = analyseCorrelation(config, groupedSCS, geneName='KCNQ1OT', groupName='patientAll', removeNoExpressCells = 'yes', percentage=20)
plotCorrelationVolcano(config, correlationResult=correlationResultKCNQ1OT1, coefficientCutoff=c(-0.25,0.25), pValueCutoff=1/10^8, outputMode='pdf')
correlationResultNPPA = analyseCorrelation(config, groupedSCS, geneName='NPPA', groupName='patientAll', removeNoExpressCells = 'yes', percentage=20)
plotCorrelationVolcano(config, correlationResult=correlationResultNPPA, coefficientCutoff=c(-0.15,0.15), pValueCutoff=1/10^8, outputMode='pdf')
correlationResultNPPB = analyseCorrelation(config, groupedSCS, geneName='NPPB', groupName='patientAll', removeNoExpressCells = 'yes', percentage=20)
plotCorrelationVolcano(config, correlationResult=correlationResultNPPB, coefficientCutoff=c(-0.15,0.15), pValueCutoff=1/10^8, outputMode='pdf')
correlationResultRPL3L = analyseCorrelation(config, groupedSCS, geneName='RPL3L', groupName='patientAll', removeNoExpressCells = 'yes', percentage=25)
plotCorrelationVolcano(config, correlationResult=correlationResultRPL3L, coefficientCutoff=c(-0.25,0.25), pValueCutoff=1/10^8, outputMode='pdf')

# Analyse correlation for two genes
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','NPPB'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','^TTN_'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','RYR2'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','MYOM1'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','IGFBP2'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('NPPA','ACTC1'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','CMYA5'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('RPL3L','CALU'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','NPPA'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','NPPB'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','XIRP2'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','MYH7_'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
analyseCorrelationAB(config, groupedSCS, geneNames=c('^TTN_','MYH6'), correlationMethod = 'pearson', groupNames=NULL, outputMode='pdf')
for(geneOfInterest in correlationResultKCNQ1OT1$correlationData[which(correlationResultKCNQ1OT1$correlationData$correlation<=-coefficientCutoff & correlationResultKCNQ1OT1$correlationData$pValue<=pValueCutoff),'genes']) {  #Only functional if correlationCutoff & pValueCutoff manually defined
  analyseCorrelationAB(config, groupedSCS, geneNames=c('KCNQ1OT1',geneOfInterest), correlationMethod='pearson', groupNames=NULL, outputMode = 'save')
}

# Exploratory analysis of index data from FACS
# Load data
plateJE10 = read.csv('/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
# Fix well names
plateJE10$Well = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Well = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
# Merge into one dataframe
indexData = rbind(plateJE10, plateJE11)
# Fix column names
colnames(indexData) <- c('Cell','Events','Parent','FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')
# Convert measurement columns to numeric
indexData$FSC_A = as.numeric(unlist(lapply(indexData$FSC_A, function(x) {str_replace(x, ',','.')})))
indexData$FSC_W = as.numeric(unlist(lapply(indexData$FSC_W, function(x) {str_replace(x, ',','.')})))
indexData$FSC_H = as.numeric(unlist(lapply(indexData$FSC_H, function(x) {str_replace(x, ',','.')})))
indexData$SSC_A = as.numeric(unlist(lapply(indexData$SSC_A, function(x) {str_replace(x, ',','.')})))
indexData$SSC_W = as.numeric(unlist(lapply(indexData$SSC_W, function(x) {str_replace(x, ',','.')})))
indexData$SSC_H = as.numeric(unlist(lapply(indexData$SSC_H, function(x) {str_replace(x, ',','.')})))
indexData$BV421_A = as.numeric(unlist(lapply(indexData$BV421_A, function(x) {str_replace(x, ',','.')})))

# Load indexData as metaData
loadMetaData(config, groupedSCS, indexData) 

# Plot tSNE-maps
plotMetaDataTSNE(config, groupedSCS, metadataColumns = 'FSC_A', groupNames = NULL, colorScheme = 'Spectral', expressionScale = 'linear', overrideDir = FALSE, overrideTheme = NULL, outputMode = 'pdf')

# Get correlations with gene expression
FSCACorrelation <- analyseMetaDataCorrelation(config, groupedSCS, parameterName = 'FSC_A', groupName ='patientAllWithIndex', removeNoExpressCells = 'yes', percentage=10)
plotCorrelationVolcano(config, FSCACorrelation, coefficientCutoff = c(-0.4,0.2), pValueCutoff = 5/10^3, overrideDir = F, outputMode = 'show')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYH7_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndex', excludeNoExprCells = T, overrideDir = F, outputMode = 'show')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYH7B_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndex', excludeNoExprCells = T, overrideDir = F, outputMode = 'show')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYL2_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndex', excludeNoExprCells = T, overrideDir = F, outputMode = 'show')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYL7_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndex', excludeNoExprCells = T, overrideDir = F, outputMode = 'show')

# 
KCNQ1OT1Data = as.double(groupedSCS$patientAll@ndata['KCNQ1OT1__chr11',])
KCNQ1OT1Data[which(KCNQ1OT1Data > 15)] = 15
ggplot(
) + geom_histogram(
  data=data.frame('KCNQ1OT1' = KCNQ1OT1Data),
  aes(
    x=KCNQ1OT1
  )
)


configSets = list(
  #'minExp 0.1 - 100 cells - 20 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 20,'minCellExp' = 0.1), #11 clusters
  #'minExp 0.1 - 100 cells - 25 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 25,'minCellExp' = 0.1), #13 clusters
  'minExp 0.1 - 100 cells - 30 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 30,'minCellExp' = 0.1), #6 clusters 
  'minExp 0.1 - 100 cells - 35 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 35,'minCellExp' = 0.1), #6 clusters
  'minExp 0.1 - 100 cells - 40 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 40,'minCellExp' = 0.1) #6 clusters
)
configSets = list(
  #'minExp 0.1 - 5 per cent - 20 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 20,'minCellExp' = 0.1), # >10 clusters
  #'minExp 0.1 - 5 per cent - 30 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 30,'minCellExp' = 0.1), # 9 clusters, 4 very small ones
  'minExp 0.1 - 5 per cent - 40 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 40,'minCellExp' = 0.1), # 5 clusters
  'minExp 0.1 - 5 per cent - 50 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 50,'minCellExp' = 0.1) # 6 clusters
)
configSets = list(
  'minExp 0.1 - 10 per cent - 20 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 20,'minCellExp' = 0.1), # 6 clusters
  'minExp 0.1 - 10 per cent - 30 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 30,'minCellExp' = 0.1), # 6 clusters
  'minExp 0.1 - 10 per cent - 40 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 40,'minCellExp' = 0.1) # 6 clusters
) 
configSets = list(
  #'minExp 0.2 - 100 cells - 20 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 20,'minCellExp' = 0.2), # 12 clusters
  #'minExp 0.2 - 100 cells - 25 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 25,'minCellExp' = 0.2), # 12 clusters
  'minExp 0.2 - 100 cells - 30 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 30,'minCellExp' = 0.2), # 6 clusters
  #'minExp 0.2 - 100 cells - 35 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 35,'minCellExp' = 0.2), # 14 clusters
  'minExp 0.2 - 100 cells - 40 correlations' = list('minCellFraction' = 100/1615,'minCorrelations' = 40,'minCellExp' = 0.2) # 6 clusters
)
configSets = list(
  #'minExp 0.2 - 5 per cent - 20 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 20,'minCellExp' = 0.2), #  13 clusters
  'minExp 0.2 - 5 per cent - 30 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 30,'minCellExp' = 0.2), # 6 clusters
  'minExp 0.2 - 5 per cent - 40 correlations' = list('minCellFraction' = 0.05,'minCorrelations' = 40,'minCellExp' = 0.2), # 6 clusters
  #'minExp 0.2 - 5 per cent - 50 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 50,'minCellExp' = 0.2) # 9 Clusters
)
configSets = list(
  #'minExp 0.2 - 10 per cent - 20 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 20,'minCellExp' = 0.2), # 7 clusters
  #'minExp 0.2 - 10 per cent - 30 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 30,'minCellExp' = 0.2), # 7 clusters
  'minExp 0.2 - 10 per cent - 40 correlations' = list('minCellFraction' = 0.10,'minCorrelations' = 40,'minCellExp' = 0.2) # 6 clusters
) 

for(configSetName in names(configSets)) {
  outputDir = configSetName
  configSet = configSets[[configSetName]]
  # Generate correlation matrix
  correlationMatrices = generateCorrelationMatrix(config, groupedSCS, groupNames=c('patientAllMod'), excludeOutlierCells=T, minCellFraction=configSet$minCellFraction, minCellExpression=configSet$minCellExp, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf')
  # Plot number of genes at different cutoffs
  genesPerCutoffPlot = ggplot(data=correlationMatrices$patientAllMod$genesPerCutoff) + geom_line(aes(x=cutoff, y=nGenes)) + ylim(c(0,800)) + xlim(c(0, 80))
  ggsave(paste0(outputDir,'/genesPerCutoffPlot.pdf'), plot=genesPerCutoffPlot)
  # Determine regulons from the correlation matrix
  regulons <- analyzeRegulons(config, correlationMatrices, minCorrelations=configSet$minCorrelations, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir=F, outputMode='pdf')
  # Plot the clustering gap statistic
  gapStatPlot = ggplot()+geom_point(data=data.frame(x=1:nrow(regulons$patientAllMod$gap_stat[[1]]),gap=regulons$patientAllMod$gap_stat[[1]][,'gap']), aes(x=x, y=gap))
  ggsave(paste0(outputDir,'/gapStatPlot.pdf'), plot=gapStatPlot)
  # Plot regulon gene networks
  plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', plotLabels=F, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  #plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='seham', colorBy='Regulon', edgeCategory='direction', plotLabels=F, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  #plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='eigen', colorBy='Regulon', edgeCategory='direction', plotLabels=F, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  #plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Cluster', edgeCategory='direction', plotLabels=F, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  #plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='TF', edgeCategory='direction', plotLabels=F, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', plotLabels=T, plotFullMatrix=F, overrideDir=F, outputMode='pdf')
  #plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', plotLabels=F, plotFullMatrix=T, overrideDir=F, outputMode='pdf')
  # Analyze regulon gene ontology
  #regulonGO <- analyzeRegulonGeneOntology(config, groupedSCS, correlationMatrices, regulons, pathwayPCutoff=0.05, GOKegg='GO', includeChildTerms=F)
  #plotRegulonGeneOntology(config, groupedSCS, regulonGO, regulons, topX=10, plotRegulatedGenes=FALSE, overrideDir=FALSE, outputMode = 'pdf')
}



TTNUpLabels = c('RYR2','CMYA5','XIRP2','DST','SVIL','ZFP106','DMD','MYOM1','DSP','SORBS1','SYNM','NRAP','NEXN','NEBL','MYH7','SYNM', 'EIF3A','EIF5B', 'TRDN')
TTNDownLabels = c('ACTA1','TNNI3','CRYAB','ACTC1','MYL2','HSPB1','MYL7','TPM1','TNNC1','LDB3','CSRP3','PLN','MYL3','LDB3', 'CSRP3')
TTNLabels = unlist(c(TTNUpLabels,TTNDownLabels, 'TTN'))
regulonLabels = unlist(returnTopDifferentialGenesPerCluster(config, groupedSCS, clusterDifferentialExpression, numGenes = 20)$patientAllMod[1:5])
regulonLabels = str_replace_all(regulonLabels, "__chr(\\d+|[MXY])",'')
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', plotLabels=F, plotFullMatrix=T, overrideDir=F, outputMode='show')
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', plotLabels=regulonLabels, plotFullMatrix=F, overrideDir=F, outputMode='show')
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='TF', edgeCategory='direction', plotLabels=humanTranscriptionFactors, plotFullMatrix=F, overrideDir=F, outputMode='show')
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Cluster', edgeCategory='direction', plotLabels=regulonLabels, plotFullMatrix=F, overrideDir=F, outputMode='show')






### Look into some better cutoffs for discarding the KCNQ1OT1/MALAT1-high cells
# Plot original clustering
plotClusterTSNE(config, groupedSCS, groupNames='patientAll', includeOutlierClusters = T, includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS1/B', overrideTheme = TSNEImprovementTheme, outputMode = 'show')
# Just exclude cluster 7&8
groupedData$patientAllNo78 <- groupedData$patientAll[,which((colnames(groupedData$patientAll) %in% names(groupedSCS$patientAll@cpart[which(groupedSCS$patientAll@cpart < 7)])))]
# Make a sum of KNCQ1OT1 & MALAT1
groupedSCS$patientAll@ndata['MALAT1OT1',] = groupedSCS$patientAll@ndata['MALAT1_',] + groupedSCS$patientAll@ndata['KCNQ1OT1',]
plotGeneExpressionTSNE(config, groupedSCS, geneNames=c('MALAT1_','KCNQ1OT1','MALAT1OT1'), groupNames='patientAll', includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', overrideDir='Fig2/C', pointSizeModifier=1, overrideTheme = TSNEImprovementTheme, outputMode='show')
plotGeneExpressionCutoffTSNE(config, groupedSCS, geneNames = 'MALAT1OT1', overUnder='Over', cutoff=10, groupNames='patientAllMALAT1OT1', includeOutlierCells = T, pointSizeModifier=1, overrideDir = 'FigS1/C',overrideTheme = TSNEImprovementTheme, outputMode = 'show')
# Only take main cluster cells with low MALAT1/KCNQ1OT1
groupedData$patientAllMALAT1OT1 <- groupedData$patientAll[,which((colnames(groupedData$patientAll) %in% names(groupedSCS$patientAll@cpart[which(groupedSCS$patientAll@cpart < 9)])))]
groupedData$patientAllMALAT1OT1 = groupedData$patientAllMALAT1OT1[,which(!(colnames(groupedData$patientAllMALAT1OT1) %in% colnames(groupedSCS$patientAll@ndata)[which(as.double(groupedSCS$patientAll@ndata['MALAT1OT1',])>15)]))]
# Cluster new data
groupedSCSMod <- buildSCSObject(config, groupedData, groupNames=c('patientAllNo78','patientAllMALAT1OT1'))
groupedSCS$patientAllNo78 = groupedSCSMod$patientAllNo78
groupedSCS$patientAllMALAT1OT1 = groupedSCSMod$patientAllMALAT1OT1
plotClusterTSNE(config, groupedSCS, groupNames=NULL, includeOutlierClusters = T, includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS1/B', overrideTheme = TSNEImprovementTheme, outputMode = 'show')
groupedSCS$patientAllMALAT1OT1@ndata['MALAT1OT1',] = groupedSCS$patientAllMALAT1OT1@ndata['MALAT1',] + groupedSCS$patientAllMALAT1OT1@ndata['KCNQ1OT1',]
plotGeneExpressionTSNE(config, groupedSCS, geneNames=c('MALAT1_','KCNQ1OT1','MALAT1OT1'), groupNames='patientAllMALAT1OT1', includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', overrideDir='Fig2/C', pointSizeModifier=1, overrideTheme = TSNEImprovementTheme, outputMode='show')
