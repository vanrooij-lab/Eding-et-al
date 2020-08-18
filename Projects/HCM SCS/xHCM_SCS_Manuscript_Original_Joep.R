
######################################################################
# Original version of manscript figures by Joep.
######################################################################

######################################################################

# Used for backwards compatibilty with the older R-versions (random number sampling is changed in R 3.6)
RNGkind(sample.kind="Rounding")

# Load the scripts that add in the functions
# Determine output folder
{
  #Set scriptsDirectory
  scriptsDirectory = '/Users/j.eding/Workg/Scripts/SingleCellSeq/'
  
  #Load required script files. #TODO: make finding these automatically easier (based on enforcing the folder structure?)
  source(paste0(scriptsDirectory,"/Functions/RaceID2_StemID_class.R"))
  source(paste0(scriptsDirectory,"/Functions/AnalysisFunctions.R"))
  
  #Generate output dir, include date & time to prevent overwriting. 
  outputDir <- paste0('analysis_', format(Sys.time(), "%Y-%m-%d-%Hh_%Mm"))
}

# Make default configuration
{
  #Default options, here for convenience
  options(stringsAsFactors = FALSE) # Changing this may break a lot of things

  #Set working directory (this directory contains all data files)
  setwd("/Users/j.eding/Workg/Projects/SCS in HCM/SCS-data/AllDataFiles")
  
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
  config$samples[['patient5']] <- c('AL1', 'AL2')
  
  #Make a list of groups and define by which samples they are made up.
  config$groups <- list()
  config$groups[['patientAll']] <- c('patient1', 'patient2', 'patient3', 'patient4', 'patient5')
  config$groups[['patientAllWithIndex']] <- c('patient1', 'patient2', 'patient3', 'patient5', 'patient4')
  config$groups[['patient1']] <- c('patient1')
  config$groups[['patient2']] <- c('patient2')
  config$groups[['patient3']] <- c('patient3')
  config$groups[['patient4']] <- c('patient4')
  config$groups[['patient5']] <- c('patient5')
  
  #Make a list of colors, some function will use colors matching to a group or sample in graphing
  config$colors <- list()
  config$colors[['patient1']] <- "#67AD5D"
  config$colors[['patient2']] <- "#E88B3F"
  config$colors[['patient3']] <- "#D76A9D"
  config$colors[['patient4']] <- "#CCCCCC"
  config$colors[['patient5']] <- "#333333"
  config$colors[['patient1Mod']] <- "#67AD5D"
  config$colors[['patient2Mod']] <- "#E88B3F"
  config$colors[['patient3Mod']] <- "#D76A9D"
  config$colors[['patient4Mod']] <- "#CCCCCC"
  config$colors[['patient5Mod']] <- "#333333"
  
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
groupedData <- loadData(config, c('patientAll', 'patientAllWithIndex', 'patient1', 'patient2', 'patient3', 'patient4', 'patient5'))
groupedDataBackup <- groupedData
# Apply exclusions
groupedData <- applyExclusions(config, groupedData)

# Build SCSObjects for the groupedData
groupedSCS <- buildSCSObject(config, groupedData, groupNames = 'patientAll')

# Perform quality control on the processed data
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')

#### General figure settings ####
TSNEImprovementTheme = theme(
  axis.text = element_text(size=14, colour = 'black'),
  axis.text.x = element_text(angle = 0, hjust = 0.5)
)

#### Figure 1 ####
# Anne makes this. Figure consists of schematics and histology, no R graphs as of yet.

#### Supplementary Figure 1 ####
# Figure S1A (Show distribution of read counts per cell, show cutoff for cell inclusion)
cutoffControl(config, groupedData, 'patientAll', numBins = 100, cutoff=1000, overrideDir = 'FigS1/A2', overrideTheme = NULL, outputMode = 'pdf')

# Figure S1B (Show t-SNE with all cells included)
plotClusterTSNE(config, groupedSCS, groupNames='patientAll', includeOutlierClusters = T, includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS1/B', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure S1C (exclude all cells with KCNQ1OT1 > 10)
plotGeneExpressionCutoffTSNE(config, groupedSCS, geneNames = 'KCNQ1OT1', overUnder='Over', cutoff=10, groupNames='patientAll', includeOutlierCells = T, pointSizeModifier=1, overrideDir = 'FigS1/C',overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')
cellsToExclude = colnames(groupedSCS$patientAll@ndata)[which(as.double(groupedSCS$patientAll@ndata['KCNQ1OT1__chr11',])>10)]

# Re-analyse data without these cells
groupedData$patientAllMod = groupedData$patientAll[,which(!(colnames(groupedData$patientAll) %in% cellsToExclude))]
groupedData$patient1Mod = groupedData$patient1[,which(!(colnames(groupedData$patient1) %in% cellsToExclude))]
groupedData$patient2Mod = groupedData$patient2[,which(!(colnames(groupedData$patient2) %in% cellsToExclude))]
groupedData$patient3Mod = groupedData$patient3[,which(!(colnames(groupedData$patient3) %in% cellsToExclude))]
groupedData$patient4Mod = groupedData$patient4[,which(!(colnames(groupedData$patient4) %in% cellsToExclude))]
groupedData$patient5Mod = groupedData$patient5[,which(!(colnames(groupedData$patient5) %in% cellsToExclude))]
groupedData$patientAllWithIndexMod = groupedData$patientAllWithIndex[,which(!(colnames(groupedData$patientAllWithIndex) %in% cellsToExclude))]
# Build all SCS-objects
groupedSCSMod <- buildSCSObject(config, groupedData, groupNames = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllMod','patientAllWithIndexMod'))
groupedSCSMod$patientAll <- groupedSCS$patientAll
groupedSCS = groupedSCSMod
remove(groupedSCSMod)

#Figure S1D: clustering t-SNE with outlier clusters
plotClusterTSNE(config, groupedSCS, groupNames='patientAllMod', includeOutlierClusters = T, includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS1/D', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

#Figure S1E: 
#See 2B

# Figure S1F:
dir.create(paste0(outputDir,'/FigS1/F'), recursive=T)
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')  # This method doesn't currently support the overrideDir parameter, hence the symlinking.
file.symlink(paste0(getwd(),'/',outputDir,'/patientAllMod/Quality control/processed/patientAllMod.pdf'), paste0(outputDir,'/FigS1/F/patientAllModQualityControl.pdf'))


#### Figure 1 ####
# Figure consists of histology and FACS and is made by Anne.


#### Figure 2 #### 
# Figure 2A
plotClusterTSNE(config, groupedSCS, groupNames='patientAllMod', includeOutlierClusters = F, includeOutlierCells = F, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='Fig2/A', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 2B
clusterDifferentialExpression <- analyseClusterDifferentialExpression(config, groupedSCS, groupNames = NULL, overrideDir='Fig2/D', outputMode = 'save')
compareClustersAcrossGroups(config, groupedSCS, referenceGroup='patientAllMod', groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllMod'), orientation='vertical', includeOutliers=FALSE, overrideDir='Fig2/B', outputMode='pdf')
compareClustersAcrossGroups(config, groupedSCS, referenceGroup='patientAllMod', groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllMod'), orientation='vertical', includeOutliers=TRUE, overrideDir='FigS1/F', outputMode='pdf')

# Figure 2C
plotGeneExpressionTSNE(config, groupedSCS, geneNames='MYH7_', groupNames='patientAllMod', includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', overrideDir='Fig2/C', pointSizeModifier=1, overrideTheme = TSNEImprovementTheme, outputMode='pdf')

# Figure 2D
# Generated in 2B

# Figure 2E
plotGeneExpressionTSNE(config, groupedSCS, geneNames=c('NPPB','MYL2','^TTN_'), groupNames='patientAllMod', includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', pointSizeModifier=1, overrideDir='Fig2/E', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 2F/G
# See figure S2

#### Supplementary Database 1 ####
# Clusters-specific expression
dir.create(paste0(outputDir,'/DatabaseS1/'), recursive=T)
file.symlink(paste0(getwd(),'/',outputDir,'/Fig2/D/patientAllMod/clusterSpecificExpression.xlsx'), paste0(outputDir,'/DatabaseS1/clusterSpecificExpression_shortcut_to_2-D.xlsx'))
# TODO: Get gene lists based on a cutoff for supplementary table 1.


#### Supplementary Figure 2 ####
clusterGO = analyseClusterGeneOntology(config, groupedSCS, clusterDifferentialExpression, groupNames = 'patientAllMod', includeOutlierClusters = F, minExpression=0.20, FCCutoff = 1.2, PCutoff = NULL, upDown = 'Up', pathwayPCutoff = 0.05, GOKegg='GO', includeChildTerms=F)
plotClusterGeneOntology(config, groupedSCS, clusterGO, clusterDifferentialExpression, topX=10, plotRegulatedGenes = F, overrideDir = 'FigS2/GO minExp0.2 - FC1.2Up - px', outputMode = 'pdf')

#### Figure 3 ####
# Figure 3A:
plotClusterTSNE(config, groupedSCS, groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod'), includeOutlierClusters = F, includeOutlierCells = F, colorScheme = 'Spectral', pointSizeModifier=1.5, overrideDir='Fig3/A', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 3B:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = '^TTN_', groupNames = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1.5, overrideDir='Fig3/B', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 3C: Venn diagrams for the overlap between TTN clusters
library('UpSetR')
dir.create(paste0(outputDir,'/Fig3/C/'), recursive=T)
pdf(paste0(outputDir,'/Fig3/C/TTN_cluster_overlap.pdf'))
TTNClusterGenes = list(
  p1.4 = (clusterDifferentialExpression$patient1Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::filter(fc > 1.5))$rowname,
  p2.4 = (clusterDifferentialExpression$patient2Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::filter(fc > 1.5))$rowname,
  p3.8 = (clusterDifferentialExpression$patient3Mod$cl.8 %>% tibble::rownames_to_column() %>% dplyr::filter(fc > 1.5))$rowname,
  p4.4 = (clusterDifferentialExpression$patient4Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::filter(fc > 1.5))$rowname,
  p5.1 = (clusterDifferentialExpression$patient5Mod$cl.1 %>% tibble::rownames_to_column() %>% dplyr::filter(fc > 1.5))$rowname
)
UpSetData <- fromList(TTNClusterGenes[seq(length(TTNClusterGenes),1)])
UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  UpSetData,
  empty.intersections = NULL,
  order.by = "freq",
  group.by = "degree",
  nsets = 36,
  nintersects = 252,
  mb.ratio = c(0.5,0.5),
  queries = list(
    list(query = intersects, params = list('p1.4','p2.4','p3.8','p4.4','p5.1'), color = "green", active = T)
  )
)
dev.off()

# Figure 3D: Histology/Immunochemistry
# Anne takes care of this

#### Supplementary Database 2 ####
p1.4 = clusterDifferentialExpression$patient1Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::arrange(desc(fc))
p2.4 = clusterDifferentialExpression$patient2Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::arrange(desc(fc))
p3.8 = clusterDifferentialExpression$patient3Mod$cl.8 %>% tibble::rownames_to_column() %>% dplyr::arrange(desc(fc))
p4.4 = clusterDifferentialExpression$patient4Mod$cl.4 %>% tibble::rownames_to_column() %>% dplyr::arrange(desc(fc))
p5.1 = clusterDifferentialExpression$patient5Mod$cl.1 %>% tibble::rownames_to_column() %>% dplyr::arrange(desc(fc))
dir.create(paste0(outputDir,'/DatabaseS2/'), recursive=T)
ClusterGenesWB = createWorkbook()
addWorksheet(ClusterGenesWB, 'patient 1 - cluster 4')
writeDataTable(ClusterGenesWB, 'patient 1 - cluster 4', p1.4)
addWorksheet(ClusterGenesWB, 'patient 2 - cluster 4')
writeDataTable(ClusterGenesWB, 'patient 2 - cluster 4', p2.4)
addWorksheet(ClusterGenesWB, 'patient 3 - cluster 8')
writeDataTable(ClusterGenesWB, 'patient 3 - cluster 8', p3.8)
addWorksheet(ClusterGenesWB, 'patient 4 - cluster 4')
writeDataTable(ClusterGenesWB, 'patient 4 - cluster 4', p4.4)
addWorksheet(ClusterGenesWB, 'patient 5 - cluster 1')
writeDataTable(ClusterGenesWB, 'patient 5 - cluster 1', p5.1)
saveWorkbook(ClusterGenesWB, paste0(outputDir,'/DatabaseS2/','TTNclusters.xlsx')) 


#### Figure 4 ####
# Figure 4A
plotGeneExpressionTSNE(config, groupedSCS, geneNames = 'NPPA', groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig4/A', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 4B: Histology/Immuno

# Figure 4C:
correlationResultNPPA= analyseCorrelation(config, groupedSCS, geneName='^NPPA_', groupName='patientAllMod', removeNoExpressCells = 'no', removeOutlierCells=T, percentage=10)
NPPAUpLabels = c('NPPB','ACTC1','IGFBP2','XIRP1','RTN4','ACTR1A','CYB5R3','GATA4','MEF2A','CRYAB')
NPPADownLabels = c('RYR2','TTN','MYOM1')
NPPALabels = unlist(c(NPPAUpLabels,NPPADownLabels))
plotCorrelationVolcano(config, correlationResult=correlationResultNPPA, coefficientCutoff=c(-0.1,0.1), pValueCutoff=1/10^2, overrideLabels=NPPALabels ,overrideDir='Fig4/C', outputMode='pdf')

# Figure 4D:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = c('NPPB', 'ACTC1', 'CRYAB'), groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig4/D', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 4E:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = c('RYR2'), groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig4/E', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 4F:
NPPAposGO <- analyseCorrelationGeneOntology(config, groupedSCS, correlationResultNPPA, minExpression=0.2, corrCutoff=0.1, pCutoff=1/10^2, posNeg='Pos', pathwayPCutoff=0.05)
plotCorrelationGeneOntology(config, groupedSCS, GOResults = NPPAposGO, corrResults=correlationResultNPPA, topX=10, plotRegulatedGenes=FALSE, overrideDir=paste0('Fig4/F/minExp','0.2'), outputMode = 'pdf') 

# Figure 4G:
NPPAnegGO <- analyseCorrelationGeneOntology(config, groupedSCS, correlationResultNPPA, minExpression=0.2, corrCutoff=0.1, pCutoff=1/10^2, posNeg='Neg', pathwayPCutoff=0.05)
plotCorrelationGeneOntology(config, groupedSCS, GOResults = NPPAnegGO, corrResults=correlationResultNPPA, topX=10, plotRegulatedGenes=FALSE, overrideDir=paste0('Fig4/G/minExp','0.2'), outputMode = 'pdf') 


#### Supplementary Database 3 ####
# Correlation analysis with NPPA
dir.create(paste0(outputDir,'/DatabaseS3/'), recursive=T)
NPPACorrelationWB = createWorkbook()
addWorksheet(
  NPPACorrelationWB,
  'NPPA'
)
writeDataTable(
    NPPACorrelationWB,
    sheet='NPPA',
    correlationResultNPPA$correlationData
)
saveWorkbook(NPPACorrelationWB, paste0(outputDir,'/DatabaseS3/','NPPA_correlation.xlsx')) 


#### Figure 5 #### 
# Fig 5A: Calculate correlation matrix
correlationMatrices = generateCorrelationMatrix(config, groupedSCS, groupNames=c('patientAllMod'), excludeOutlierCells=T, minCellFraction=0.05, minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf')
regulons <- analyzeRegulons(config, correlationMatrices, minCorrelations=40, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig5/A', outputMode='pdf')
gapStatPlot = ggplot()+geom_point(data=data.frame(x=1:nrow(regulons$patientAllMod$gap_stat[[1]]),gap=regulons$patientAllMod$gap_stat[[1]][,'gap']), aes(x=x, y=gap))
ggsave(paste0(outputDir,'/Fig5/A/gapStatPlot.pdf'), plot=gapStatPlot)

# Fig 5B: Plot gene correlation network
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', whichCorrelations = "Pos", plotLabels=F, plotFullMatrix=F, overrideDir='Fig5/B', outputMode='pdf')

# Fig 5C: Top correlated genes per regulon
# See excel file in the Fig5A folder
# Print transcription factors per cluster to 
{
  cat("Transcription factors in: ")
  for(i in names(regulons$patientAllMod$regulons)) {
      cat(paste0("\n - ",i,": "))
      regulons$patientAllMod$regulons[[i]][which(gsub("__chr(\\d+|[MXY])", '', regulons$patientAllMod$regulons[[i]]) %in% humanTranscriptionFactors)] %>% cat
  }
}
# Transcription factor motif enrichment analysis: Done by Kees Boogerd
# Gene-ontology for regulons
regulonGO <- analyzeRegulonGeneOntology(config, groupedSCS, correlationMatrices, regulons, pathwayPCutoff=0.05, GOKegg='GO', includeChildTerms=F)
plotRegulonGeneOntology(config, groupedSCS, regulonGO, regulons, topX=10, plotRegulatedGenes=FALSE, overrideDir='Fig5/C', outputMode = 'pdf')


#### Supplementary Database 4 ####
dir.create(paste0(outputDir,'/DatabaseS4/'), recursive=T)
file.symlink(paste0(getwd(),'/',outputDir,'/Fig5/A/regulonHeatmap_patientAllMod_minCell0.05_minExp0.1_rCutoff_minCorr40_squish0.01'), paste0(outputDir,'/DatabaseS4/regulonHeatmap_patientAllMod_minCell0.05_minExp0.1_rCutoff_minCorr40_squish0.01.xlsx'))


#### Figure 6 ####
# Figure 6A: Pie-chart with genetic make-up of VUmc cohort
# Figure 6B: Bulk qPCR for phenotype validation controls vs HCM
# Figure 6C: Collagens in bulk RNA

# Load qPCR data
qpcrData <- read.xlsx(
  '/Volumes/dataniob/Group VanRooij/Group members/Joep - Maya - Anne/RNA bulk data/RNA samples VUmc/myectomy samples VU qPCR.xlsx',
  sheet=1, 
  startRow=2, 
  skipEmptyCols = T,
  na.strings = c('NA','N/A',"N/A","NA")
)
rownames(qpcrData) <- qpcrData$FHCM # Fix rownames
qpcrData$group = 'HCM'
qpcrData[which(grepl('Control', rownames(qpcrData))),'group'] = 'Control'
# Force all columns to numeric
for(n in colnames(qpcrData)[2:(ncol(qpcrData)-1)]) {
  qpcrData[[n]] <- as.numeric(qpcrData[[n]])
}
# Iterate over values to remove
for(y in colnames(qpcrData)[2:(ncol(qpcrData)-1)]) {
  for(x in 1:nrow(qpcrData)) {
    if(is.na(qpcrData[x,y])) {next}
    if(qpcrData[x,y] > 50) {qpcrData[x,y] = NA}
  }
}
geneList <- list()
for(n in colnames(qpcrData)[2:ncol(qpcrData)-1]) {
  if(n %in% c('FHCM','GAPDH','HPRT','avgHG','avgHG2','KCNJ')) { next() }
  qpcrData[[n]] <- as.numeric(qpcrData[[n]]) - as.numeric(qpcrData$GAPDH)
  geneList = unlist(c(geneList,n))
}

# Figure 6DE: Correlate all qPCR genes with all other qPCR genes
dir.create(paste0(outputDir,'/Fig6/qPCR/'), recursive=T, showWarnings = F)
geneList = geneList[which(!(geneList %in% c('TNNI3','PTGDS')))] #Exclude genes with no amplification, otherwise these cause errors
for(stableGene in geneList) {
  for(goi in geneList) {
    ABCorrelation <-  ggplot(   #Start ggplot with data and grouping
      qpcrData,
      aes(y=get(stableGene), x=get(goi), color=group)
    ) + geom_point(    #Print points
      na.rm = T
    ) + scale_colour_manual(
      values=c('Control'='red','HCM'='black')
    ) + theme_classic(         #Get rid of default ggplot theme
    ) + theme(
      plot.subtitle=element_text(size=7, face="italic", color="black")
    ) + geom_smooth(
      method = 'lm',
      se=T,
      color='black',
      na.rm=T
    ) + labs(
      y=stableGene, 
      x=goi,     #Get rid of axis labels
      title=paste0(stableGene,'  vs. ',goi),
      subtitle = paste0(
        'R: ',
        round(summary(lm(get(stableGene) ~ get(goi), qpcrData))$coefficients['get(goi)',1],3),
        '     p: ',
        round(summary(lm(get(stableGene) ~ get(goi), qpcrData))$coefficients['get(goi)',4],7)
      )
    ) + TSNEImprovementTheme
    ggsave(paste0(stableGene,'_vs_',goi,'.pdf'),path = paste0(outputDir,'/Fig6/qPCR/'))
  }
}


#### Figure 7 ####
# Figure 7A:
# Load index data
plateJE10 = read.csv('/Users/j.eding/Workg/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/j.eding/Workg/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
plateAL1 = read.csv('/Users/j.eding/Workg/Projects/SCS in HCM/SCS-data/AL1.csv', header=T ,sep=",", row.names=NULL, skip=15)
plateAL2 = read.csv('/Users/j.eding/Workg/Projects/SCS in HCM/SCS-data/AL2.csv', header=T ,sep=",", row.names=NULL, skip=15)
# Fix well names
plateJE10$Well = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Well = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1$Well = unlist(lapply(plateAL1$Well, function(x) {return(paste0('AL1_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1 = plateAL1[colnames(plateJE10)] # Throw away ridiculous column overload
plateAL2$Well = unlist(lapply(plateAL2$Well, function(x) {return(paste0('AL2_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL2 = plateAL2[colnames(plateJE10)] # Throw away ridiculous column overload
# Merge into one dataframe
indexData = rbind(plateJE10, plateJE11, plateAL1, plateAL2)
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
# Get correlations with gene expression
FSCACorrelation <- analyseMetaDataCorrelation(config, groupedSCS, parameterName = 'FSC_A', groupName ='patientAllWithIndexMod', removeNoExpressCells = 'parameter', removeOutlierCells=T, percentage=20)
FSCAUpLabels = c('ACTA1','FUNDC2','HSPB1','MYH7','TNNT2','NDUF4A','MYL2','MYL7','CRYAB','ERLIN2','CRYAB')
FSCADownLabels = c('MALAT1','LMNA','GDI2','TUBB4B','PRRC2B')
FSCALabels = unlist(c(FSCAUpLabels,FSCADownLabels))
plotCorrelationVolcano(config, FSCACorrelation, coefficientCutoff = c(-0.15,0.15), pValueCutoff = 0.05, overrideLabels = regulons$patientAllMod$regulons$regulon.2, overrideDir = 'Fig7/A', outputMode = 'pdf')
FSCA2Labels = regulons$patientAllMod$regulons$regulon.2
# Figure 7B: Dot plots for positive correlations
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'ACTA1_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig7/B', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYL2_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig7/B', outputMode = 'pdf')
# Figure 7C: Dot plots for negative correlations
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'LMNA_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig7/C', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = '^PRRC2B_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig7/C', outputMode = 'pdf')

# Figure 7D&F: WGA/ACTA1 immunofluorescence co-staining
# Figure 7E&G: Correlation of cell size to ACTA1-intensity

#### Supplementary Table 5 ####
dir.create(paste0(outputDir,'/DatabaseS5/'), recursive=T)
FSCACorrelationWB = createWorkbook()
addWorksheet(
  FSCACorrelationWB,
  'FSCA'
)
writeDataTable(
    FSCACorrelationWB,
    sheet='FSCA',
    FSCACorrelation$correlationData
)
saveWorkbook(FSCACorrelationWB, paste0(outputDir,'/DatabaseS5/','FSCA_correlation.xlsx')) 


#### Determine some numbers for use in the text ####
# Mean number of mitochondrial transcripts per cell
{
cat("\n\nSpread of mitochondrial read percentage\n")
mitoReads <- groupedDataBackup$patientAll %>% 
  dplyr::select(colnames(groupedSCS$patientAll@ndata)) %>% 
  tibble::rownames_to_column('gene') %>%
  dplyr::filter(stringr::str_detect(gene, '__chrM', negate=F)) %>%
  tibble::column_to_rownames('gene') %>%
  colSums()
geneReads <- groupedDataBackup$patientAll %>% 
  dplyr::select(colnames(groupedSCS$patientAll@ndata)) %>% 
  tibble::rownames_to_column('gene') %>%
  dplyr::filter(stringr::str_detect(gene, '(__chrM)|(ERCC-)', negate=T)) %>%
  tibble::column_to_rownames('gene') %>%
  colSums() 
qualityTable <- data.frame(
  'mito' = mitoReads,
  'gene' = geneReads,
  'fraction' = mitoReads/(mitoReads+geneReads)
)
qualityTable$group='All'
qualityTable %>% summary %>% print
mean(qualityTable$fraction) %>% print
sd(qualityTable$fraction) %>% print
ggplot2::ggplot(data=qualityTable) + ggplot2::geom_violin(ggplot2::aes(y=fraction, group=group, x=group), draw_quantiles = c(0.5))

cat("\n\nSpread of mitochondrial read percentage after KCNQ1OT1-high cell exclusion\n")
mitoReadsMod <- groupedDataBackup$patientAll %>% 
  dplyr::select(colnames(groupedSCS$patientAllMod@ndata)) %>% 
  tibble::rownames_to_column('gene') %>%
  dplyr::filter(stringr::str_detect(gene, '__chrM', negate=F)) %>%
  tibble::column_to_rownames('gene') %>%
  colSums()
geneReadsMod <- groupedDataBackup$patientAll %>% 
  dplyr::select(colnames(groupedSCS$patientAllMod@ndata)) %>% 
  tibble::rownames_to_column('gene') %>%
  dplyr::filter(stringr::str_detect(gene, '(__chrM)|(ERCC-)', negate=T)) %>%
  tibble::column_to_rownames('gene') %>%
  colSums() 
qualityTableMod <- data.frame(
  'mito' = mitoReadsMod,
  'gene' = geneReadsMod,
  'fraction' = mitoReadsMod/(mitoReadsMod+geneReadsMod)
)
qualityTableMod$group='Mod'
qualityTableMod %>% summary %>% print
mean(qualityTableMod$fraction) %>% print
sd(qualityTableMod$fraction) %>% print
ggplot2::ggplot(data=qualityTableMod) + ggplot2::geom_violin(ggplot2::aes(y=fraction, group=group, x=group), draw_quantiles = c(0.5))

cat("\n\nNumber of cells per patient before exclusion of KCNQ1OT1-high cells")
cat(paste0(
      "\n Pooled cells: ",
      ncol(groupedSCS$patientAll@ndata),
      "\n Patient 1: ",
      length(grep(paste0('(',config$samples$patient1,')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 2: ",
      length(grep(paste0('(',config$samples$patient2,')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 3: ",
      length(grep(paste0('(',config$samples$patient3,')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 4: ",
      length(grep(paste0('(',config$samples$patient4,')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 5: ",
      length(grep(paste0('(',config$samples$patient5,')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T))
 ))

cat("\n\nNumber of cells per patient after exclusion of KCNQ1OT1-high cells")
cat(paste0(
      "\n Pooled cells: ",
      ncol(groupedSCS$patientAllMod@ndata),
      "\n Patient 1: ",
      length(grep(paste0('(',config$samples$patient1,')', collapse = '|'), colnames(groupedSCS$patient1Mod@ndata), value=T)),
      "\n Patient 2: ",
      length(grep(paste0('(',config$samples$patient2,')', collapse = '|'), colnames(groupedSCS$patient2Mod@ndata), value=T)),
      "\n Patient 3: ",
      length(grep(paste0('(',config$samples$patient3,')', collapse = '|'), colnames(groupedSCS$patient3Mod@ndata), value=T)),
      "\n Patient 4: ",
      length(grep(paste0('(',config$samples$patient4,')', collapse = '|'), colnames(groupedSCS$patient4Mod@ndata), value=T)),
      "\n Patient 5: ",
      length(grep(paste0('(',config$samples$patient5,')', collapse = '|'), colnames(groupedSCS$patient5Mod@ndata), value=T))
 ))

}


# Calculate some regulon expression statistics
testData <- data.frame(
    averageExp = unlist(c(
        as.double(rowMeans(groupedSCS$patientAllMod@ndata[regulons$patientAllMod$regulons$regulon.1,])),
        as.double(rowMeans(groupedSCS$patientAllMod@ndata[regulons$patientAllMod$regulons$regulon.2,])),
        as.double(rowMeans(groupedSCS$patientAllMod@ndata[regulons$patientAllMod$regulons$regulon.3,])),
        as.double(rowMeans(groupedSCS$patientAllMod@ndata[regulons$patientAllMod$regulons$regulon.4,])),
        as.double(rowMeans(groupedSCS$patientAllMod@ndata[regulons$patientAllMod$regulons$regulon.5,]))
    )),
    group = unlist(c(
        rep('Reg1', length(regulons$patientAllMod$regulons$regulon.1)),
        rep('Reg2', length(regulons$patientAllMod$regulons$regulon.2)),
        rep('Reg3', length(regulons$patientAllMod$regulons$regulon.3)),
        rep('Reg4', length(regulons$patientAllMod$regulons$regulon.4)),
        rep('Reg5', length(regulons$patientAllMod$regulons$regulon.5))
    ))
)
ggplot(
  data=testData
) + geom_boxplot(
  aes(
    y=averageExp,
    group=group,
    x=group
  )
)
