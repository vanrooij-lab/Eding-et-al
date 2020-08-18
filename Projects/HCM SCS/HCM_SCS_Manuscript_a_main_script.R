
# load("~/Data/_2019_02_HCM_SCS/_sessions/HCM_2020may_re-analysis_tempsave-v4.Rdata")

# This is a slightly edited version of the original script that generated figuers for the "SCS HCM" manuscript from Joep.
# See Git for a change history.
#
# More extensive additional analysis are now included in additional scripts, most
# of these require (parts of the) main script to be run first (also, for some scripts
# it might be required to run parts of other scripts, I'm working to streamline this
# better):
# - HCM_SCS_Manuscript_a_regulons_per_patient.R: regulon analysis performed per patient
# - HCM_SCS_Manuscript_b_FSCA_part1.R: analysis FSCA-data and correlations
# - HCM_SCS_Manuscript_b_FSCA_part2.R: analysis FSCA-data and correlations (part 2)
# - HCM_SCS_Manuscript_b_plotting_regulons_on_tSNE.R: composite regulon expression on tSNE
# - HCM_SCS_Manuscript_c_comparing_gene_detection_patients.R: comparing gene detection between patients
# - HCM_SCS_Manuscript_c_comparing_gene_expression_patients.R: comparing gene expression between patients
#
# Some analyses that I also did but which will likely not be incorporated in the manuscript:
# - HCM_SCS_Manuscript_y_extra_check_plots.R: pre-processing validation and TTN checks.
# - HCM_SCS_Manuscript_y_extra_preprocessing_plots.R: pre-processing validation and TTN checks.
#
# I also made some scrpits earlier that I don't use any more, they can be found in the
# "old_scripts" folder (HCM SCS Manuscript_extensions_MW_fsca.R, )


# Used only for backwards compatibilty with the older R-versions (random number sampling is changed in R 3.6):
# RNGkind(sample.kind="Rounding") # not used here 

# To load a saved analysis: load libraries, load workspace, set workdirectory, e.g.
# - Execute lines which say "source(...)", l. 16-21 below
# load("~/Data_notbacked/temp_sessions/HCM_2020may_re-analysis_tempsave.Rdata")
# setwd("/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/Datafiles/")

# Load the scripts that add in the functions
# Determine output folder
{
  #Set scriptsDirectory
  scriptsDirectory = '/Users/m.wehrens/Documents/git_repos/SCS_Joep/'
  
  #Load required script files. #TODO: make finding these automatically easier (based on enforcing the folder structure?)
  # Note: don't select all these lines and then execute them with cmd+enter, that'll confuse R
  # E.g. seperately executing them works 
  source(paste0(scriptsDirectory,"/Functions/RaceID2_StemID_class.R"))
  source(paste0(scriptsDirectory,"/Functions/AnalysisFunctions.R"))
  # Some additional functions that are needed (these file locations are not convenient -- todo)
  source(paste0(scriptsDirectory,"/Functions/functions_mw_copy/MW_standard_plots.R"))
  source(paste0(scriptsDirectory,"/Functions/functions_mw_copy/MW_analysis_functions.R"))
  source(paste0(scriptsDirectory,"/Functions/functions_mw_copy/MW_some_colors.R"))
  source(paste0(scriptsDirectory,"/Functions/functions_mw_copy/MW_GOKEGG_analysis.R"))
  # original path:
  #source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_standard_plots.R")
  #source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_analysis_functions.R")
  #source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_some_colors.R")
  #source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_GOKEGG_analysis.R")
  source(paste0(scriptsDirectory,"/Functions/shorthandFunctionsMW.R"))
  
  #Set working directory (this directory contains all data files)
  setwd("/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/Datafiles/")
  
  #Generate output dir, include date & time to prevent overwriting. 
  outputDir <- paste0('analysis_', format(Sys.time(), "%Y-%m-%d-%Hh_%Mm"))
}

# Make default configuration
{
  #Default options, here for convenience
  options(stringsAsFactors = FALSE) # Changing this may break a lot of things
  
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
  config$groups[['patientAllWithIndex']] <- c('patient4', 'patient5')
  config$groups[['patient1']] <- c('patient1')
  config$groups[['patient2']] <- c('patient2')
  config$groups[['patient3']] <- c('patient3')
  config$groups[['patient4']] <- c('patient4')
  config$groups[['patient5']] <- c('patient5')
  
  #Make a list of colors, some function will use colors matching to a group or sample in graphing
  config$colors <- list()
  config$colors[['patient1']] <- "deepskyblue4"
  config$colors[['patient2']] <- "darksalmon"
  config$colors[['patient3']] <- "darkgray"
  config$colors[['patient4']] <- "darkred"
  config$colors[['patient5']] <- "chocolate2"
  config$colors[['patient1Mod']] <- "deepskyblue4"
  config$colors[['patient2Mod']] <- "darksalmon"
  config$colors[['patient3Mod']] <- "darkgray"
  config$colors[['patient4Mod']] <- "darkred"
  config$colors[['patient5Mod']] <- "chocolate2"
  
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
groupedData <- loadData(config, c('patientAll', 'patientAllWithIndex', 'patient1', 'patient2', 'patient3','patient4', 'patient5'))
groupedDataBackup <- groupedData
# Apply exclusions
groupedData <- applyExclusions(config, groupedData)

# Build SCSObjects for the groupedData
groupedSCS <- buildSCSObject(config, groupedData, groupNames = 'patientAll')

# Perform quality control on the processed data
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')

# See end of this file for some stats on the scRNA-seq data (nr of cells, nr of reads, etc)

#### General figure settings ####
TSNEImprovementTheme = theme(
  axis.text = element_text(size=14, colour = 'black'),
  axis.text.x = element_text(angle = 0, hjust = 0.5)
)

# (MW) Let's contemplate which clusters or genes we should remove from analysis
# before continuing analysis
clusterDifferentialExpression <- analyseClusterDifferentialExpression(config, groupedSCS, groupNames = NULL, overrideDir='FigMW/original_diff', outputMode = 'save')

#### Supplementary Figure 1 ####
# Figure S1A gating strategy
# Figure S1B DRAQ5 assay
# Figure S1C pictures pre and after sort
# Figure S1D RIN values of patient 1-3

#### Supplementary Figure 2 ####
# Figure S2A: Martijn added this code?


# Figure S2B (Show distribution of read counts per cell, show cutoff for cell inclusion)
cutoffControl(config, groupedData, 'patientAll', numBins = 100, cutoff=1000, overrideDir = 'FigS2/AB2', overrideTheme = NULL, outputMode = 'pdf')

# Figure S2C (Show t-SNE with all cells included)
plotClusterTSNE(config, groupedSCS, groupNames='patientAll', includeOutlierClusters = T, includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS2/C', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# ---
# (Modified by MW)
# Initially, cells with KCNQ1OT1 > 10 were removed; now the aim is to remove cells
# with KCNQ1OT1 > 10 OR cells that are member of a KCNQ1OT1 enriched cluster.
# Also, the KCNQ1OT1 and MALAT1 genes will be removed from the analysis beforehand.
# ---

# Figure S2D (exclude all cells with KCNQ1OT1 > 10) ---> THIS NEEDS TO BE RE-MADE
#plotGeneExpressionCutoffTSNE(config, groupedSCS, geneNames = 'KCNQ1OT1', overUnder='Over', cutoff=10, 
#                             groupNames='patientAll', includeOutlierCells = T, pointSizeModifier=1, overrideDir = 'FigS2/D',
#                             overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')
KCNQ1OT1_positive_cells =  groupedSCS$patientAll@ndata['KCNQ1OT1__chr11',]>10

# Now look at fold-change for KCNQ1OT1 per cluster
fc_K = sapply(clusterDifferentialExpression$patientAll, function(X) {  X['KCNQ1OT1__chr11',]$fc })
# Get names of clusters with higher than 2 fold change
KCNQ1OT1_cluster_cells = paste0('cl.',groupedSCS$patientAll@cpart) %in% names(fc_K[fc_K>2])
#KCNQ1OT1_cluster_cells =  groupedSCS$patientAll@cpart %in% c(10,11,13,16)

cellsToExclude = colnames(groupedSCS$patientAll@ndata)[KCNQ1OT1_positive_cells|KCNQ1OT1_cluster_cells]

# print info
print(paste0('# Cells excluded: ',length(cellsToExclude),' (',
  round(100*length(cellsToExclude)/length(KCNQ1OT1_positive_cells),2),'%)',
  ', leaving ', length(KCNQ1OT1_positive_cells)-length(cellsToExclude),' for analysis'))

# make plot
p=shorthand_selection(the_selection = as.logical(KCNQ1OT1_positive_cells|KCNQ1OT1_cluster_cells), current_SCS = groupedSCS$patientAll,
                    mytitle = 'Original',sel_title = 'Exclusion',mypointsize = 1.5)+
                    shorthand_tsne_joeptheme()+xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
p
shorthand_save(p,thedir = paste0(outputDir,'/FigS2/D'),filename = 'selection.pdf',mywidth = 7*2.54, myheight = 5.3*2.54)

genesToExclude = c('KCNQ1OT1__chr11','MALAT1__chr11')

# Remove those cells from the data
groupedData$patientAllMod = groupedData$patientAll[,!(colnames(groupedData$patientAll) %in% cellsToExclude)]
groupedData$patientAllWithIndexMod = groupedData$patientAllWithIndex[,!(colnames(groupedData$patientAllWithIndex) %in% cellsToExclude)]
groupedData$patient1Mod = groupedData$patient1[,!(colnames(groupedData$patient1) %in% cellsToExclude)]
groupedData$patient2Mod = groupedData$patient2[,!(colnames(groupedData$patient2) %in% cellsToExclude)]
groupedData$patient3Mod = groupedData$patient3[,!(colnames(groupedData$patient3) %in% cellsToExclude)]
groupedData$patient4Mod = groupedData$patient4[,!(colnames(groupedData$patient4) %in% cellsToExclude)]
groupedData$patient5Mod = groupedData$patient5[,!(colnames(groupedData$patient5) %in% cellsToExclude)]
# Remove the genes from the data
config$custom_gene_blacklist = genesToExclude
groupedData <- applyExclusions(config, groupedData)
# Re-analyse data without these cells
# Build all SCS-objects
groupedSCSMod <- buildSCSObject(config, groupedData, groupNames = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllWithIndexMod', 'patientAllMod'))
groupedSCSMod$patientAll <- groupedSCS$patientAll
groupedSCS = groupedSCSMod
remove(groupedSCSMod)

#Figure S2E: clustering t-SNE with outlier clusters
plotClusterTSNE(config, groupedSCS, groupNames='patientAllMod', includeOutlierClusters = T, 
                includeOutlierCells = T, colorScheme = 'Spectral', pointSizeModifier=1, overrideDir='FigS2/E', 
                overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

#Figure S2F: 
#See 1D

# Figure S2G:
dir.create(paste0(outputDir,'/FigS2/G'), recursive=T)
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')  # This method doesn't currently support the overrideDir parameter, hence the symlinking.
file.symlink(paste0(getwd(),'/',outputDir,'/patientAllMod/Quality control/processed/patientAllMod.pdf'), paste0(outputDir,'/FigS2/G/patientAllModQualityControl.pdf'))


#### Figure 1 ####
# Figure 1A schematic myectomy
# Figure 1B Histology myectomy tissue

# Figure 1C
plotClusterTSNE(config, groupedSCS, groupNames='patientAllMod', includeOutlierClusters = F, includeOutlierCells = F, colorScheme = 'Spectral', pointSizeModifier=1, 
                overrideDir='Fig1/C', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 1D


clusterDifferentialExpression <- analyseClusterDifferentialExpression(config, groupedSCS, groupNames = NULL, overrideDir='Fig1/D', outputMode = 'save')
compareClustersAcrossGroups(config, groupedSCS, referenceGroup='patientAllMod', 
                            groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllMod'), 
                            orientation='vertical', includeOutliers=FALSE, overrideDir='Fig1/D', outputMode='pdf')
compareClustersAcrossGroups(config, groupedSCS, referenceGroup='patientAllMod', 
                            groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllMod'), 
                            orientation='vertical', includeOutliers=TRUE, overrideDir='FigS2/F', outputMode='pdf')

# Figure 1E
plotGeneExpressionTSNE(config, groupedSCS, geneNames='MYH7_', groupNames='patientAllMod', includeOutlierCells = F, colorScheme='Spectral', 
                       useSquish = 0.01, expressionScale='linear', overrideDir='Fig1/E', pointSizeModifier=1, overrideTheme = TSNEImprovementTheme, outputMode='pdf')

# Note: one could also use the alternative function
shorthand_plot_patient_contr_clusters(config, groupedSCS$patientAllMod)

# Figure 1F
# Generated in 1C

###########check top expressed genes in clusters############

# Figure 1G
plotGeneExpressionTSNE(config, groupedSCS, geneNames=c('^RPS12', 'FHL1','NPPB','MYL2','^TTN_','ACTA1', 'KRT5_', 'HSPB6', 'NPPA','ZFP106','NEAT1'), groupNames='patientAllMod', 
                       includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', pointSizeModifier=1, overrideDir='Fig1/G', 
                       overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')
# Also with alt color scheme
plotGeneExpressionTSNE(config, groupedSCS, geneNames=c('^RPS12', 'FHL1','NPPB','MYL2','^TTN_','ACTA1', 'KRT5_', 'HSPB6', 'NPPA','ZFP106','NEAT1'), groupNames='patientAllMod', 
                       includeOutlierCells = F, colorScheme=col_viridius_inferno_inv, useSquish = 0.01, expressionScale='linear', pointSizeModifier=1, overrideDir='Fig1/G/alternativeColorScheme/', 
                       overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')
# Plus some additional plots
plotGeneExpressionTSNE(config, groupedSCS, geneNames=paste0('^',c('ZNF106', 'CHGA', 'HSPB6', 'TCAP', 'RYR2'),'_'), groupNames='patientAllMod', 
                       includeOutlierCells = F, colorScheme='Spectral', useSquish = 0.01, expressionScale='linear', pointSizeModifier=1, overrideDir='FigX/', 
                       overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')
# And alternative color scheme
plotGeneExpressionTSNE(config, groupedSCS, geneNames=paste0('^',c('ZNF106', 'CHGA', 'HSPB6', 'TCAP', 'RYR2'),'_'), groupNames='patientAllMod', 
                       includeOutlierCells = F, colorScheme=col_viridius_inferno_inv, useSquish = 0.01, expressionScale='linear', pointSizeModifier=1, overrideDir='FigX/otherColorScheme/', 
                       overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

for (gene in c('ZNF106', 'CHGA', 'HSPB6', 'TCAP', 'RYR2')) {
  p=shorthand_expression(gene, groupedSCS$patientAllMod, paste0(gene, ' expression'))
  shorthand_save(p, 'analysis_2020-05-04-18h_17m/FigX/', paste0('mw_tSNE_expr_',gene,'.pdf'))
}

# Also some numbers of min-max changes in gene expression
TTN_expr=get_expression_gene2(expression_matrix = as.matrix(groupedSCS$patientAllMod@ndata),gene_query = 'TTN')
TTN_percentiles0199 = calc_limits(TTN_expr)
TTN_percentiles0199[2]/TTN_percentiles0199[1]

#### Supplementary Database 1 ####
# Clusters-specific expression
dir.create(paste0(outputDir,'/DatabaseS1/'), recursive=T)
file.symlink(paste0(getwd(),'/',outputDir,'/Fig1/F/patientAllMod/clusterSpecificExpression.xlsx'), paste0(outputDir,'/DatabaseS1/clusterSpecificExpression_shortcut_to_2-D.xlsx'))
# TODO: Get gene lists based on a cutoff for supplementary table 1.


#### Supplementary Figure 3 ####
clusterGO = analyseClusterGeneOntology(config, groupedSCS, clusterDifferentialExpression, groupNames = 'patientAllMod', 
                                       includeOutlierClusters = F, minExpression=0.20, FCCutoff = 1.2, PCutoff = NULL, upDown = 'Up', pathwayPCutoff = 0.05, GOKegg='GO', includeChildTerms=F)

plotClusterGeneOntology(config, groupedSCS, clusterGO, clusterDifferentialExpression, topX=10, plotRegulatedGenes = F, overrideDir = 'FigS3/GO minExp0.2 - FC1.2Up - px', outputMode = 'pdf')

#### Figure 2 ####
# Figure 2A:
plotClusterTSNE(config, groupedSCS, groupNames=c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod'), includeOutlierClusters = F, 
                includeOutlierCells = F, colorScheme = 'Spectral', pointSizeModifier=1.5, overrideDir='Fig2/A', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 2B:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = '^TTN_', groupNames = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod'), 
                       includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1.5, overrideDir='Fig2/B', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 2C: Venn diagrams for the overlap between TTN clusters
# Now let's get the top-expressed genes per cluster, per patient
clusterDifferentialExpression[[current_patient]][[X]]$fc
rownames(clusterDifferentialExpression$patient1Mod$cl.1[order(clusterDifferentialExpression$patient1Mod$cl.1$fc),])[1]
# Note loops over patients in following list:
groupNames_toCheck = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod')
cluster_top_DE_pp=list()
for (current_patient in groupNames_toCheck) {
  # gather top-expressed genes and their fold-change
  cluster_top_DE_pp[[current_patient]] = sapply(
          names(clusterDifferentialExpression[[current_patient]]),
          function(X) {c(rownames(clusterDifferentialExpression[[current_patient]][[X]][order(clusterDifferentialExpression[[current_patient]][[X]]$fc, decreasing = T),])[1],
                         clusterDifferentialExpression[[current_patient]][[X]][order(clusterDifferentialExpression[[current_patient]][[X]]$fc, decreasing = T),]$fc[1])})
}

# Note that if there are more than one TTN-enriched clusters, you might need
# to manually inspect cluster_top_DE_pp
for (current_patient in groupNames_toCheck) {
  print(paste0(current_patient, '-->', names(cluster_top_DE_pp[[current_patient]][1,
    cluster_top_DE_pp[[current_patient]][1,]=='TTN__chr2'])))
}
 
library('UpSetR')
dir.create(paste0(outputDir,'/Fig2/C/'), recursive=T)
pdf(paste0(outputDir,'/Fig2/C/TTN_cluster_overlap.pdf'))
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

# Figure 2D: Histology/Immunochemistry
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
saveWorkbook(ClusterGenesWB, paste0(outputDir,'/DatabaseS2/','TTNclusters.xlsx'), overwrite = T) 


#### Figure 3 ####
# Figure 3A
plotGeneExpressionTSNE(config, groupedSCS, geneNames = 'NPPA', groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, 
                       expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig3/A', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 3B: Histology/Immuno

# Figure 3C:
correlationResultNPPA= analyseCorrelation(config, groupedSCS, geneName='^NPPA_', groupName='patientAllMod', removeNoExpressCells = 'no', removeOutlierCells=T, percentage=10)
NPPAUpLabels = c('NPPB','ACTC1','IGFBP2','XIRP1','RTN4','ACTR1A','CYB5R3','GATA4','MEF2A','CRYAB')
NPPADownLabels = c('RYR2','TTN','MYOM1')
NPPALabels = unlist(c(NPPAUpLabels,NPPADownLabels))
# Copy original parameter
correlationResultNPPA_corrected = correlationResultNPPA
# Now perform adjustment here
# save old values
correlationResultNPPA_corrected$correlationData$pValue.old=correlationResultNPPA_corrected$correlationData$pValue
# put corrected values in pValue slot
correlationResultNPPA_corrected$correlationData$pValue=p.adjust(correlationResultNPPA_corrected$correlationData$pValue.old, method="BH")
plotCorrelationVolcano(config, correlationResult=correlationResultNPPA_corrected, coefficientCutoff=c(-0.1,0.1), pValueCutoff=1/10^2, overrideLabels=NPPALabels ,overrideDir='Fig3/C', outputMode='pdf')

# Some info printed:
NPPA_significant_corrs = correlationResultNPPA_corrected$correlationData[correlationResultNPPA_corrected$correlationData$pValue<.01,]
print(paste0(sum(NPPA_significant_corrs$correlation>0), ' positive and ',
        sum(NPPA_significant_corrs$correlation<0), ' negative genes had p<0.01 for NPPA correlations'))

# Figure 3D:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = c('NPPB', 'ACTC1', 'CRYAB', 'RTN4_', 'IGFBP2'), groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig3/D', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 3E:
plotGeneExpressionTSNE(config, groupedSCS, geneNames = c('RYR2'), groupNames = c('patientAllMod'), includeOutlierCells = F, colorScheme = 'Spectral', useSquish = 0.01, expressionScale = 'linear', pointSizeModifier=1, overrideDir='Fig3/E', overrideTheme = TSNEImprovementTheme, outputMode = 'pdf')

# Figure 3F:
# Note: I changed the correlation minimum to 0 (no cutoff), because the cutoff can be based on the p-values
NPPAposGO <- analyseCorrelationGeneOntology(config, groupedSCS, correlationResultNPPA_corrected, minExpression=0.2, corrCutoff=0, pCutoff=1/10^2, posNeg='Pos', pathwayPCutoff=0.05)
plotCorrelationGeneOntology(config, groupedSCS, GOResults = NPPAposGO, corrResults=correlationResultNPPA_corrected, topX=10, plotRegulatedGenes=FALSE, overrideDir=paste0('Fig3/F/minExp','0.2'), outputMode = 'pdf') 
#plot_volcano(correlationResultNPPA$correlationData)

#NPPAposGO_corrected <- analyseCorrelationGeneOntology(config, groupedSCS, correlationResultNPPA_corrected, minExpression=0.2, corrCutoff=0.1, pCutoff=1/10^2, posNeg='Pos', pathwayPCutoff=0.05)
#plotCorrelationGeneOntology(config, groupedSCS, GOResults = NPPAposGO, corrResults=correlationResultNPPA, topX=10, plotRegulatedGenes=FALSE, overrideDir=paste0('Fig3/F/minExp','0.2_p.adjust'), outputMode = 'pdf') 

# Addition MW: With selection for GO terms that are related to min. 3 genes
p1=shorthand_GOplot(NPPAposGO$patientAllMod$NPPA__chr1, minGeneCount = 1, mycolor = '#E8442A')+
  theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())+
  xlab(element_blank())+give_better_textsize_plot(8)
p1
shorthand_save(p1,thedir = paste0(outputDir,'/Fig3/F'),filename = 'GOterms_NPPA_pos_notselgenes_mw.pdf',mywidth = 12.5, myheight = 7.5) # mywidth = 7*2.54, myheight = 5.3*2.54
p2=shorthand_GOplot(NPPAposGO$patientAllMod$NPPA__chr1, minGeneCount = 3, mycolor = '#E8442A')+
  theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())+
  xlab(element_blank())+give_better_textsize_plot(8)
p2
shorthand_save(p2,thedir = paste0(outputDir,'/Fig3/F'),filename = 'GOterms_NPPA_pos_sel3genes_mw.pdf',mywidth = 12.5, myheight = 7.5) # mywidth = 7*2.54, myheight = 5.3*2.54

# Figure 3G:
# Note: I changed the correlation minimum to 0 (no cutoff), because the cutoff can be based on the p-values
NPPAnegGO <- analyseCorrelationGeneOntology(config, groupedSCS, correlationResultNPPA_corrected, minExpression=0.2, corrCutoff=0, pCutoff=1/10^2, posNeg='Neg', pathwayPCutoff=0.05)
plotCorrelationGeneOntology(config, groupedSCS, GOResults = NPPAnegGO, corrResults=correlationResultNPPA_corrected, topX=10, plotRegulatedGenes=FALSE, overrideDir=paste0('Fig3/G/minExp','0.2'), outputMode = 'pdf') 

# Addition MW: With selection for GO terms that are related to min. 3 genes
p=shorthand_GOplot(NPPAnegGO$patientAllMod$NPPA__chr1, minGeneCount = 3, mycolor = 'darkblue')+
  theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())+
  xlab(element_blank())+give_better_textsize_plot(8)
p
shorthand_save(p,thedir = paste0(outputDir,'/Fig3/G'),filename = 'GOterms_NPPA_neg_sel3genes_mw.pdf',mywidth = 12.5, myheight = 7.5) # mywidth = 7*2.54, myheight = 5.3*2.54
# also mw plot w/o selection
p2=shorthand_GOplot(NPPAnegGO$patientAllMod$NPPA__chr1, minGeneCount = 1, mycolor = 'darkblue')+
  theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())+
  xlab(element_blank())+give_better_textsize_plot(8)
p2
shorthand_save(p2,thedir = paste0(outputDir,'/Fig3/G'),filename = 'GOterms_NPPA_neg_not_sel_genes_mw.pdf',mywidth = 12.5, myheight = 7.5) # mywidth = 7*2.54, myheight = 5.3*2.54

#### Supplementary Database 3 ####

# Let's do this cleanly, by adding p.adjust parameter to the original data again (to have better param names)
correlationResultNPPA$correlationData$pVal.adj = p.adjust(correlationResultNPPA$correlationData$pValue, method="BH")

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
saveWorkbook(NPPACorrelationWB, paste0(outputDir,'/DatabaseS3/','NPPA_correlation.xlsx'), overwrite = T) 


#### Figure 4 #### 
# Fig 4A: Calculate correlation matrix
correlationMatrices = generateCorrelationMatrix(config, groupedSCS, groupNames=c('patientAllMod'), excludeOutlierCells=T, minCellFraction=0.05, 
                                                minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf')
regulons <- analyzeRegulons(config, correlationMatrices, minCorrelations=40, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig4/A', outputMode='pdf')
gapStatPlot = ggplot()+geom_point(data=data.frame(x=1:nrow(regulons$patientAllMod$gap_stat[[1]]),gap=regulons$patientAllMod$gap_stat[[1]][,'gap']), aes(x=x, y=gap))
ggsave(paste0(outputDir,'/Fig4/A/gapStatPlot.pdf'), plot=gapStatPlot)

# Fig 4B: Plot gene correlation network
plotRegulons(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', whichCorrelations = "Pos", 
             plotLabels=F, plotFullMatrix=F, overrideDir='Fig4/B', outputMode='pdf')

# Fig 4C: Top correlated genes per regulon
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

# The code below would generate GO terms for all of the separate regulons;
# but we would like to look at GO terms for the shared regulons
if (F) { # so we don't need the code here:
  regulonGO <- analyzeRegulonGeneOntology(config, groupedSCS, correlationMatrices, regulons, pathwayPCutoff=0.05, GOKegg='GO', includeChildTerms=F)
  plotRegulonGeneOntology(config, groupedSCS, regulonGO, regulons, topX=10, plotRegulatedGenes=FALSE, overrideDir='Fig4/C', outputMode = 'pdf')
} # instead, look at "HCM_SCS_Manuscript_a_regulons_per_patient.R"


#### Supplementary Database 4 ####
dir.create(paste0(outputDir,'/DatabaseS4/'), recursive=T)
file.symlink(paste0(getwd(),'/',outputDir,'/Fig4/A/regulonHeatmap_patientAllMod_minCell0.05_minExp0.1_rCutoff_minCorr40_squish0.01'), 
             paste0(outputDir,'/DatabaseS4/regulonHeatmap_patientAllMod_minCell0.05_minExp0.1_rCutoff_minCorr40_squish0.01.xlsx'))


#### Figure 5 ####
# Figure 5A: Pie-chart with genetic make-up of VUmc cohort
# Figure 5B: Bulk qPCR for phenotype validation controls vs HCM
# Figure 5C: Collagens in bulk RNA

# Load qPCR data
qpcrData <- read.xlsx(
  #'/Users/a.leeuw/Documents/Single cell sequencing/20200326R_Analysis/Datafiles/myectomy samples VU qPCR.xlsx',
  'myectomy samples VU qPCR.xlsx',
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

# Figure 5DE: Correlate all qPCR genes with all other qPCR genes
dir.create(paste0(outputDir,'/Fig5/qPCR/'), recursive=T, showWarnings = F)
geneList = geneList[which(!(geneList %in% c('TNNI3','PTGDS')))] #Exclude genes with no amplification, otherwise these cause errors
#yGeneList = c('ACTC1', 'XIRP2')
for (DontPlotControl in c(F, T)) {
  if (DontPlotControl){
    data_selection = !(qpcrData$group == 'Control')
    controlString= '_noControl_'
  } else { data_selection = T; controlString='' }
  # plot corr plots for all correlations
  for(stableGene in geneList) {
    for(goi in geneList) {
      ABCorrelation <-  ggplot(   #Start ggplot with data and grouping
        qpcrData[data_selection,],
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
          round(summary(lm(get(stableGene) ~ get(goi), qpcrData[data_selection,]))$coefficients['get(goi)',1],3),
          '     p: ',
          round(summary(lm(get(stableGene) ~ get(goi), qpcrData[data_selection,]))$coefficients['get(goi)',4],7)
        )
      ) + TSNEImprovementTheme
      if (DontPlotControl) {ABCorrelation=ABCorrelation+theme(legend.position = 'none')}
      ggsave(paste0(controlString,stableGene,'_vs_',goi,'.pdf'),path = paste0(outputDir,'/Fig5/qPCR/'))
    }
  }
}

#### Figure 6 ####
# Figure 6A:
# Load index data
plateJE10 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
plateAL1 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL1.csv', header=T ,sep=",", row.names=NULL, skip=15)
plateAL2 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL2.csv', header=T ,sep=",", row.names=NULL, skip=15)
# Fix well names
plateJE10$Well = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Well = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1$Well = unlist(lapply(plateAL1$Well, function(x) {return(paste0('AL1_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1 = plateAL1[colnames(plateJE10)] # Throw away ridiculous column overload
plateAL2$Well = unlist(lapply(plateAL2$Well, function(x) {return(paste0('AL2_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL2 = plateAL2[colnames(plateJE10)] # Throw away ridiculous column overload
# Merge into one dataframe
indexData = rbind(plateAL1, plateAL2)
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
FSCACorrelation <- analyseMetaDataCorrelation(config, groupedSCS, parameterName = 'FSC_A', groupName ='patient5Mod', 
                                              removeNoExpressCells = 'parameter', removeOutlierCells=T, percentage=20)
FSCAUpLabels = c('ACTA1','FUNDC2','HSPB1','MYH7','TNNT2','NDUF4A','MYL2','MYL7','CRYAB','ERLIN2','CRYAB','UQCR10',
                 'MYL12A','TMEM261','MDH1_','NDUFB4_', 'PLN_', 'ENAH_')
FSCADownLabels = c('MALAT1','LMNA','GDI2','TUBB4B','PRRC2B')
FSCALabels = unlist(c(FSCAUpLabels,FSCADownLabels))
# Copy original parameter
correlationResultFSCA_corrected = FSCACorrelation
# Now perform adjustment here
# save old values
correlationResultFSCA_corrected$correlationData$pValue.old=correlationResultFSCA_corrected$correlationData$pValue
# put corrected values in pValue slot
correlationResultFSCA_corrected$correlationData$pValue=p.adjust(correlationResultFSCA_corrected$correlationData$pValue.old, method="BH")
plotCorrelationVolcano(config, FSCACorrelation, coefficientCutoff = c(-0.15,0.15), pValueCutoff = 0.05, overrideLabels = regulons$patientAllMod$regulons$regulon.2, overrideDir = 'Fig6/A', outputMode = 'pdf')
#plotCorrelationVolcano(config, FSCACorrelation, coefficientCutoff = c(-0.15,0.15), pValueCutoff = 0.05, overrideDir = 'Fig6/A', outputMode = 'pdf')
FSCA2Labels = regulons$patientAllMod$regulons$regulon.2

# Figure 6B: Dot plots for positive correlations
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'ACTA1_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/B', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYL2_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/B', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'MYL12A_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/B', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'UQCR10_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/B', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = 'LMNA_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/C', outputMode = 'pdf')
analyseCorrelationMetaDataGene(config, groupedSCS, parameterName = 'FSC_A', geneName = '^PRRC2B_', correlationMethod = 'pearson', groupNames = 'patientAllWithIndexMod', excludeNoExprCells = F, removeOutlierCells=T, plotLine = T, pointSizeModifier=1.5, overrideDir = 'Fig6/C', outputMode = 'pdf')

# Figure 6C&D: WGA/ACTA1 immunofluorescence co-staining
# Figure 6E&F: Correlation of cell size to ACTA1-intensity

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
      length(grep(paste0('(',paste0(config$samples$patient1,'_'),')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 2: ",
      length(grep(paste0('(',paste0(config$samples$patient2,'_'),')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 3: ",
      length(grep(paste0('(',paste0(config$samples$patient3,'_'),')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 4: ",
      length(grep(paste0('(',paste0(config$samples$patient4,'_'),')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T)),
      "\n Patient 5: ",
      length(grep(paste0('(',paste0(config$samples$patient5,'_'),')', collapse = '|'), colnames(groupedSCS$patientAll@ndata), value=T))
 ))

cat("\n\nNumber of cells per patient after exclusion of KCNQ1OT1-high cells")
cat(paste0(
      "\n Pooled cells: ",
      ncol(groupedSCS$patientAllMod@ndata),
      "\n Patient 1: ",
      length(grep(paste0('(',paste0(config$samples$patient1,'_'),')', collapse = '|'), colnames(groupedSCS$patient1Mod@ndata), value=T)),
      "\n Patient 2: ",
      length(grep(paste0('(',paste0(config$samples$patient2,'_'),')', collapse = '|'), colnames(groupedSCS$patient2Mod@ndata), value=T)),
      "\n Patient 3: ",
      length(grep(paste0('(',paste0(config$samples$patient3,'_'),')', collapse = '|'), colnames(groupedSCS$patient3Mod@ndata), value=T)),
      "\n Patient 4: ",
      length(grep(paste0('(',paste0(config$samples$patient4,'_'),')', collapse = '|'), colnames(groupedSCS$patient4Mod@ndata), value=T)),
      "\n Patient 5: ",
      length(grep(paste0('(',paste0(config$samples$patient5,'_'),')', collapse = '|'), colnames(groupedSCS$patient5Mod@ndata), value=T))
 ))

}

# nr of unique genes
nrow(groupedSCS$patientAllMod@ndata)

# Show some stats to show in the manuscript:
cells_kept = colnames(groupedSCS$patientAll@expdata) %in% colnames(groupedSCS$patientAll@ndata) # needed for some calculations
# First numbers of cells that were kept:
p_source_all = shorthand_give_p_source(config, current_SCS = groupedSCS$patientAll)
total_cell_count_pt = aggregate(data.frame(p_source_all, count=1)$count, by=list(p_source_all), FUN=sum)
# Then number of unique genes detected
dim(groupedSCS$patientAll@expdata)
# Total transcript counts per patient
p_source_all_raw = shorthand_give_p_source_forcellnames(config, colnames(groupedSCS$patientAll@expdata))
total_transcript_pt = aggregate(data.frame(p_source_all_raw[cells_kept], count=apply(groupedSCS$patientAll@expdata[,cells_kept],2,sum))$count, by=list(p_source_all_raw[cells_kept]), FUN=sum)
round(total_transcript_pt$x)
# Transcript count after cell-to-cell normalization
# (This is meaningless)
#total_transcript_pt_ndata = aggregate(data.frame(p_source_all, count=apply(groupedSCS$patientAll@ndata,2,sum))$count, by=list(p_source_all), FUN=sum)
#round(total_transcript_pt_ndata$x)
# Mean cell transcript total per patient
transcript_totals_pt = aggregate(data.frame(p_source_all_raw[cells_kept], count=apply(groupedSCS$patientAll@expdata[,cells_kept],2,sum))$count, by=list(p_source_all_raw[cells_kept]), FUN=mean)
round(transcript_totals_pt$x)

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
