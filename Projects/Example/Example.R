# Load the scripts that add in the functions
# Determine output folder
{
  #Set scriptsDirectory
  scriptsDirectory = '/Users/j.eding/Work/Scripts/SingleCellSeq/'     # Change this to where you have the R-scripts on your computer
  
  #Load required script files. #TODO: make finding these automatically easier (based on enforcing the folder structure?)
  source(paste0(scriptsDirectory,"/Functions/RaceID2_StemID_class.R")) # Don't change this
  source(paste0(scriptsDirectory,"/Functions/AnalysisFunctions.R")) # Don't change this
  
  #Generate output dir, include date & time to prevent overwriting. 
  outputDir <- paste0('analysis_', format(Sys.time(), "%Y-%m-%d-%Hh_%Mm")) # Don't change this unless you know what you're doing.
}


# Make default configuration
{
  #Default options, here for convenience
  options(stringsAsFactors = FALSE) # Changing this may break a lot of things

  #Set working directory (this directory contains all data files)
  setwd("/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/AllDataFiles") # Change this to where your data files are
  
  #Initialize a config variable. This variable will hold _all_ configuration for the rest of the analysis.
  config <- list() # Don't change this
  
  #Add scriptsDirectory to config so it's available within function too
  config$scriptsDirectory = scriptsDirectory # Don't change this
  
  #Make a list of samples and which plates should be in which sample
  ###Example: config$samples[['sampleName']] <- c('namePlate1', 'namePlate2', namePlateX')
  ###In the example, the filename for namePlate1 is namePlate1_TranscriptCounts.tsv
  config$samples <- list() # Don't change this
  config$samples[['patient1']] <- c('JE1', 'JE2', 'JE3', 'JE4') # Make a list of your samples ('patient1') and which plates have the cells from this samplje ('JE1' until 'JE4)
  config$samples[['patient2']] <- c('JE5', 'JE6', 'JE7', 'JE8') # Make a list of your samples ('patient2') and which plates have the cells from this samplje ('JE5' until 'JE8)
  config$samples[['patient3']] <- c('MW5', 'MW6', 'MW7', 'MW8') # Make a list of your samples ('patient3') and which plates have the cells from this samplje ('MW5' until 'MW8)
  config$samples[['patient4']] <- c('JE10', 'JE11') # Make a list of your samples ('patient4') and which plates have the cells from this samplje ('JE10' until 'JE11)
  
  #Make a list of groups and define by which samples they are made up.
  config$groups <- list() # Don't change this
  config$groups[['patientAll']] <- c('patient1', 'patient2', 'patient3') # Make a list of your groups. Your groups are what is being analysed. Here we're pooling patient1 to 3 in one group for analysis of everything together.
  config$groups[['patient1']] <- c('patient1')   # It is possible to make multiple groups. Groups can also consist of just one sample. Samples can be in multiple groups
  config$groups[['patient2']] <- c('patient2')
  config$groups[['patient3']] <- c('patient3')
  config$groups[['patient4']] <- c('patient4')
  config$groups[['patientAllWithIndex']] <- c('patient1', 'patient2', 'patient3', 'patient4')
  
  # Enter gene identifier type, pick one from:
  #   - "ensembl_gene_id"               (identifier looks like 'ENSMUSG00000000001')
  #   - "external_gene_name"            (identifier looks like 'Gnai3')
  #   - "external_gene_name_with_chr"   (identifier looks like 'Gnai3__chr1')
  #   - "entrezgene"                    (identifier looks like '14679')
  config$geneIdentifierType <- "external_gene_name_with_chr" # Change this to match one of the 4 options above
  
  # Specify species, pick one from:
  #   - "human"
  #   - "mouse"
  config$species <- "human" # Change this to match one of the 2 options above
  
  #Exclude mitochondrial genes? Set to T(RUE) or F(ALSE)
  config$excludeMitochondrialGenes = T
  
  ###
  ### Configure data selection and clustering
  ###
  #Define minimum total number of unique reads per cell (default = 1000, check read distribution over cells to pick a definite value)
  config$minTotalReads = 1000
  #Define p-value at which a cell is considered an outlier within its own cluster (default = 1e-7)
  config$outlierPValue = 1e-7
}


# Load data (This function takes some time, the amount of time depends on the number of groups you load and the number of plates that are in each group)
# This command loads all of your data into a variable that's called 'groupedData'
groupedData <- loadData(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  c('patientAll','patient1','patient2','patient3','patient4','patientAllWithIndex') # Make a list here of all the groups that you want to analyze. Any groups you list here need to be specified in the config variable in the configuration part of this script.
)

# Apply exclusions
# This function removes some doubled gene names that arise from loading additional gene identifier names,
groupedData <- applyExclusions(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedData # This needs to exactly match the variable name of the variable that has the data from loadData (generated in the step above).
)

# Build SCSObjects for the groupedData
# This function makes the single cell sequencing models that you need for all further analysis functions
groupedSCS <- buildSCSObject(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedData, # This needs to exactly match the variable name of the variable that has the data from loadData (generated in the step above).
  groupNames = c('patientAll','patient1','patient2','patient3','patient4','patientAllWithIndex') # Make a list here of all the groups that you want to analyze. Any groups you list here need to also be loaded with loadData.
)

# Perform quality control on the processed data
processedDataQualityControl(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS, # This needs to exactly match the variable name of the variable that has the the generated SCS models created in the function above.
  groupNames = NULL, # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  outputMode = 'save' # Leave as 'save', will make a pdf file with your quality control graphs.
)

# General figure settings. This is an optional variable to add to the styling of your graphs. Google 'ggplot2 theme()' to see what's possible.
# Changing this allows you to easily influence graph text size etc.
TSNEImprovementTheme = theme(
  axis.text = element_text(size=14, colour = 'black'),
  axis.text.x = element_text(angle = 0, hjust = 0.5)
)


###
###
###   PLOT SOME t-SNE MAPS
###
###


# Plot clustering on a tSNE map
plotClusterTSNE(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS, # This needs to exactly match the variable name of the variable that has the the generated SCS models created in the function above.
  groupNames = c('patientAll'), # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  includeOutlierClusters = T, # Set as T if you want to include outlier clusters. Set to F if you want to exclude outlier clusters, then outlier cells are plotted in the cluster they originated from.
  includeOutlierCells = T, # Set as T if you want to plot outlier cells. Set to F if you don't want to plot the outlier cells.
  colorScheme = 'Spectral', # Sets the color scheme. See this link for possible values: https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
  pointSizeModifier=1, # Changes the size of the dots in the tSNE. Bigger numbers give bigger dots. Can have decimals (e.g. 0.01 is a valid value)
  overrideDir=FALSE, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  overrideTheme = TSNEImprovementTheme, # Include the graphics improvement theme defined above.
  outputMode = 'show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)

# Analyze cluster differential expression
clusterDifferentialExpression <- analyseClusterDifferentialExpression(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  groupNames = NULL, # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  overrideDir=FALSE, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode = 'save'
)

# Plot gene expression on a t-SNE-map
plotGeneExpressionTSNE(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  geneNames='MYH7_', # Specify gene name to draw tSNE for. Will automatically try to find the gene name you mean, but it needs to be specific for a single gene. Use ^ to specify where the gene name begins, use _ to specify end of gene name if chromosomes are included in gene names, use $ to specify end of gene name if chromosomes are not included. Can also accept a list of gene names to plot multiple t-SNE maps, e.g: c('MYH7_','^TTN_') 
  groupNames='patientAll', # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  includeOutlierCells = F, # Set to F to exclude cells in outlier clusters, set to T to include cells in outlier clusters
  colorScheme='Spectral', # Sets the color scheme. See this link for possible values: https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/ ; Also takes the values from this package https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html#the-color-scales 
  useSquish = 0.01, # Needs a number, is interpreted as percentage (so 0.01 is 1%). Orders all values by level of expression. Sets the cutoff for the color scale at this percentage from the top and bottom. So 0.01 makes makes the lowest 1% the lowest color and the highest 1% all the highest color. This makes the outliers less visible, but gives a bigger color range for the range of values in which 98% of the samples are. Can be turned of by specifying 0 as a value.
  expressionScale='linear', # Specify scaling for expression values. Default is linear. Also accepts log, log2, log10.
  overrideDir=FALSE, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  pointSizeModifier=1, # Changes the size of the dots in the tSNE. Bigger numbers give bigger dots. Can have decimals (e.g. 0.01 is a valid value)
  overrideTheme = TSNEImprovementTheme, # Include the graphics improvement theme defined above.
  outputMode='show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)

# Plot gene expression on a t-SNE-map
plotGeneExpressionCutoffTSNE(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  geneNames = 'KCNQ1OT1', # Specify gene name to draw tSNE for. Will automatically try to find the gene name you mean, but it needs to be specific for a single gene. Use ^ to specify where the gene name begins, use _ to specify end of gene name if chromosomes are included in gene names, use $ to specify end of gene name if chromosomes are not included. Can also accept a list of gene names to plot multiple t-SNE maps, e.g: c('MYH7_','^TTN_') 
  overUnder='Under', # Specify whether genes over or under the threshold should be colored
  cutoff=10, # Specify the cutoff of the expression level
  groupNames='patientAll', # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  includeOutlierCells = T, # Set as T if you want to plot outlier cells. Set to F if you don't want to plot the outlier cells.
  pointSizeModifier=1, # Changes the size of the dots in the tSNE. Bigger numbers give bigger dots. Can have decimals (e.g. 0.01 is a valid value)
  overrideDir = FALSE, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  overrideTheme = TSNEImprovementTheme, # Include the graphics improvement theme defined above.
  outputMode = 'show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)



###
###
###   SIMPLE CORRELATION ANALYSES
###
###


# Correlate two genes with eachother
analyseCorrelationAB(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  geneNames=c('^TTN_','CMYA5'), # Specify a list of exactly 2 genes to correlate with eachother.
  correlationMethod = 'pearson', # Specify correlation method, accepts "pearson", "kendall" or "spearman"
  groupNames=NULL, # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do analysis for all groups.
  removeOutlierCells = T, # Set as F if you want to exclude  cells in outlier clusters. Set to T if you don't want to exclude the outlier cells.
  plotLine=F, # Use T to plot a line for the correlation coefficient and a 95% confidence interval. Use F to not plot a line.
  overrideDir=F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode='show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)

# Correlate 1 gene with all other genes
# This function outputs the correlation data, so you need to assign that to a new variable
correlationResultTTN = analyseCorrelation(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  geneName='^TTN_', # Specify gene name to draw tSNE for. Will automatically try to find the gene name you mean, but it needs to be specific for a single gene. Use ^ to specify where the gene name begins, use _ to specify end of gene name if chromosomes are included in gene names, use $ to specify end of gene name if chromosomes are not included. 
  groupName='patientAll', # Specify a single group name to perform this analysis for. Needs to be a groupname for which the SCS model is made. 
  cellNames=NULL, # Optional. Manually provide a list of cell names for which you want to do this analysis. 
  removeNoExpressCells = 'no', # Whether to remove some cells. 3 options: 'no' never removes any cells, 'goi' removes cells that don't express the gene specified in geneName, 'yes' removes cells that don't express either the specified gene or (determined separately for every comparison made) the gene it's being compared with. 
  removeOutlierCells=T, # Set as F if you want to exclude  cells in outlier clusters. Set to T if you don't want to exclude the outlier cells.
  percentage=20, # Minimum percentage of all cells that should remain after excluding cells. If it's less, the gene that is being compared with the specified gene is skipped.
  correlationMethod='pearson' # Specify correlation method, accepts "pearson", "kendall" or "spearman"
)
# Plot correlation analysis in a Volcano
TTNLabels = c('^TTN_', 'CMYA5', 'MYL12A')
plotCorrelationVolcano(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  correlationResult=correlationResultTTN, # Specify the variable that you saved the correlation analysis in.
  coefficientCutoff=c(-0.20,0.20), # Specify a list of 2 cutoffs for the correlation coefficient (one negative, one positive). This plot a vertical line and determines which genes get a label
  pValueCutoff=1/10^5, # Specify a single cutoff for the p-value. This draws a horizontal line on the volcano plot and determines which genes get an automatic label
  overrideLabels=TTNLabels, # Optionally specify a list of gene names (like for the t-SNE functions). If specified, only the listed genes will get a label. Leave as FALSE to plot default labels
  overrideDir=F,  # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode='show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)



###
###
###   ADVANCED CORRELATION ANALYSIS
###
###
# Correlate all genes to eachother
correlationMatrices = generateCorrelationMatrix(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  groupNames=c('patientAll'), # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  excludeOutlierCells=T, # Set to T to exclude cells in outlier clusters, set to F to include cells in outlier clusters
  minCellFraction=0.1, # Specify minimum number of cells that need to express a gene for it to be taken along in the analysis.
  minCellExpression=0.1, # Specify the expression level a cell needs to _exceed_ to be considered epxressing a gene. 0.1 (not 0!) is the value for cells that have no read for a certain gene.
  desiredPValue=0.00001, # Enter a desired p-value for correlations to be considered significant. 
  adjustP=T, # Must be either T or F. If T the function p.adjust is used to adjust all p-values for the correlations. 
  saveMatrices=F, # Whether to save the correlation matrix as a file. Best to leave it as F. Has a high chance of crashing RStudio if you set it to T.
  overrideDir=F,  # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode='pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)
# Plot number of genes at different cutoffs
genesPerCutoffPlot = ggplot(data=correlationMatrices$patientAll$genesPerCutoff) + geom_line(aes(x=cutoff, y=nGenes)) + ylim(c(0,800)) + xlim(c(0, 80))
ggsave(paste0(outputDir,'/genesPerCutoffPlot.pdf'), plot=genesPerCutoffPlot)
# Determine regulons from the correlation matrix
regulons <- analyzeRegulons(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  correlationMatrices, # Supply the variable that has the correlation matrices generated in the function above.
  minCorrelations=40, # Set the minimum number of significant correlations that a gene needs to have with other genes to be taken along in this analysis
  clusteringMethod='ward.D2', # Set the clustering method to use. Can be 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'
  overrideClusterNum=9, # Set to a number to force that number of clusters to be made. If left to F, it will run k-means clustering for a range of cluster numbers and then calculate the gap-statistic to determine the optimal number of clusters.
  useSquish=0.01,  # Needs a number, is interpreted as percentage (so 0.01 is 1%). Orders all values by level of expression. Sets the cutoff for the color scale at this percentage from the top and bottom. So 0.01 makes makes the lowest 1% the lowest color and the highest 1% all the highest color. This makes the outliers less visible, but gives a bigger color range for the range of values in which 98% of the samples are. Can be turned of by specifying 0 as a value.
  overrideDir=F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode='show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)
# Plot the clustering gap statistic
gapStatPlot = ggplot()+geom_point(data=data.frame(x=1:nrow(regulons$patientAll$gap_stat[[1]]),gap=regulons$patientAll$gap_stat[[1]][,'gap']), aes(x=x, y=gap))
ggsave(paste0(outputDir,'/gapStatPlot.pdf'), plot=gapStatPlot)
# Plot regulon gene networks
plotRegulons(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  clusterDifferentialExpression, # Specify the variable that contains the cluster differential expression analysis (used here if you want to color nodes not by which regulon they are in but for the cluster in which they are most highly expressed)
  correlationMatrices, # Specify the variable in which the correlation matrices made by generateCorrelationMatrices are saved
  regulons, # Specify the variable in which the regulons generated by analyzeRegulons are saved
  plotMode='geodist', # Specify the plotting method, this has a big influence on how your plots look. Specify 'bla' to generate an error that shows you the available plotting methods. Some plotting methods are less useful than others, this function generates a warning when you use the less useful methods.
  colorBy='Regulon', # Specify what determines the color of the nodes(=genes) in the network. Accepts 'Regulon', 'Cluster', 'TF'
  edgeCategory='direction', # Specify what determines edge category (currently only accepts 'direction', which makes positive and negative correlations get different lines)
  plotLabels=F, # Specify whether gene names should be plotted, accepts T (plots labels), F (plots no labels), or a list of gene names which causes only those gene names to be plotted.
  plotFullMatrix=F, # Specify whether the full correlation matrix should be plotted (T) or only the genes that were used for the regulon analysis (F). Plotting the full matrix is much slower and not always results in nice plots.
  overrideDir=F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode='show' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)






###
###
###   INDEX DATA
###
###
# EXPLANATION:
# To be able to couple index data from the FACSorting to your transciptomes, you first need to load this data.
# The structure of this data is probably dependent on the FACS machine/software version. 
# You need to load this data and then make sure that the script can understand which data belongs to which cell, that's what we're doing here.
# Load index data (For each plate, load the document that contains the data)
plateJE10 = read.csv('/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/j.eding/Work/Projects/SCS in HCM/SCS-data/Patient 6 (JE10&11)/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
# Fix well names (In my data, the wells are named A212, B300, C050. The lines below add the plate-specific label)
plateJE10$Well = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Well = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
# Merge into one dataframe (Now all metadata tables are merged into one table (this is possible because now there's a plate-specific prefix for each well))
indexData = rbind(plateJE10, plateJE11)
# Fix column names (Rename the columns to make sure they're called something that is easy to work with (you need to know the order of your columns to be able to do this))
colnames(indexData) <- c('Cell','Events','Parent','FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')
# Sometimes the numeric values are accidentally imported as text. The lines below need to be run for every column that should contain numeric data (so not the cell names!) and converts this column to numeric data
indexData$FSC_A = as.numeric(unlist(lapply(indexData$FSC_A, function(x) {str_replace(x, ',','.')})))
indexData$FSC_W = as.numeric(unlist(lapply(indexData$FSC_W, function(x) {str_replace(x, ',','.')})))
indexData$FSC_H = as.numeric(unlist(lapply(indexData$FSC_H, function(x) {str_replace(x, ',','.')})))
indexData$SSC_A = as.numeric(unlist(lapply(indexData$SSC_A, function(x) {str_replace(x, ',','.')})))
indexData$SSC_W = as.numeric(unlist(lapply(indexData$SSC_W, function(x) {str_replace(x, ',','.')})))
indexData$SSC_H = as.numeric(unlist(lapply(indexData$SSC_H, function(x) {str_replace(x, ',','.')})))
indexData$BV421_A = as.numeric(unlist(lapply(indexData$BV421_A, function(x) {str_replace(x, ',','.')})))
# This function loads the new metadata and combines it with the single cell data
# This function creates a new object (scsMetadata) that you won't see but is used by subsequent metadata functions
# Don't run this function multiple times with the same metadata
loadMetaData(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  indexData
) 

# Plot metadata on a t-SNE map
plotMetaDataTSNE(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  metadataColumns='FSC_A', # Specify the name of the metadata column you want to compare.
  groupNames='patientAllWithIndex', # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  colorScheme='Spectral', # Sets the color scheme. See this link for possible values: https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/ ; Also takes the values from this package https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html#the-color-scales 
  expressionScale='linear', # Specify scaling for expression values. Default is linear. Also accepts log, log2, log10.
  pointSizeModifier=1, # Changes the size of the dots in the tSNE. Bigger numbers give bigger dots. Can have decimals (e.g. 0.01 is a valid value)
  overrideDir=F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  overrideTheme = NULL, # Include the graphics improvement theme defined above.
  outputMode='pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)
# Analyse correlation between metadata and a single gene
analyseCorrelationMetaDataGene(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  parameterName = 'FSC_A', # Specify the name of the metadata column you want to compare.
  geneName = 'ACTA1_', # Specify the name of the gene that you want to compare the metadata paremeter with. Does some automatic matching, see the t-SNE functions for help with gene name matching.
  correlationMethod = 'pearson', # Specify correlation method, accepts "pearson", "kendall" or "spearman"
  groupNames = 'patientAllWithIndex', # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  excludeNoExprCells = F, # Exclude cells that don't express the gene that the metadata paremeter is being compared with.  Can be TRUE or FALSE.
  removeOutlierCells=T, # Set to T to exclude cells in outlier clusters, set to F to include cells in outlier clusters
  plotLine = T, # Use T to plot a line for the correlation coefficient and a 95% confidence interval. Use F to not plot a line.
  pointSizeModifier=1.5, # Changes the size of the dots in the tSNE. Bigger numbers give bigger dots. Can have decimals (e.g. 0.01 is a valid value)
  overrideDir = F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode = 'pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)
# Analyse correlation between metadata and all genes
FSCACorrelation <- analyseMetaDataCorrelation(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  parameterName = 'FSC_A', # Specify the name of the metadata column you want to compare.
  groupName ='patientAllWithIndex', # Specify a singe groupname to do this analysis for. 
  cellNames=NULL, # Optional. Manually provide a list of cell names for which you want to do this analysis. 
  removeNoExpressCells = 'parameter', # Either 'yes' or 'parameter'. If parameter, exclude only cells that have no value for the metadata parameter. If 'yes' excludes all cells that have no value for the metadata parameter, and for every comparison also excludes the cells that have no expression of that gene.
  removeOutlierCells=T, # Set to T to exclude cells in outlier clusters, set to F to include cells in outlier clusters
  percentage=20, # Minimum percentage of all cells that should remain after excluding cells. If it's less, the gene that is being compared with the specified gene is skipped.
  correlationMethod='pearson' # Specify correlation method, accepts "pearson", "kendall" or "spearman"
)
# Plot metadata correlation in a volcano plot (This is the same function as for plotting a volcano for the regular correlation analyses)
plotCorrelationVolcano(
  config, 
  FSCACorrelation, 
  coefficientCutoff = c(-0.2,0.2), 
  pValueCutoff = 1/10^2, 
  overrideLabels = F, 
  overrideDir = F, 
  outputMode = 'pdf' 
)




###
###
###   GENE ONTOLOGY
###
###


# Analyze gene ontology for clusters
clusterGO = analyseClusterGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  clusterDifferentialExpression, # Supply the variable name of the variable that contains the clusterdifferentialexpression data generated by analyseClusterDifferentialExpression (done further up in this script)
  groupNames = c('patientAll'), # Make a list of groupnames for which you want to do the quality control analysis. List should look like the 'groupNames' parameter of functions above. Leave as NULL to do quality control for all groups.
  includeOutlierClusters = F, # Set to T to include cells in outlier clusters, set to F to exclude cells in outlier clusters
  minExpression=0.20, # Set the minimum expression level for genes to be taken along for the background
  FCCutoff = 1.2, # Set the minimum fold change (enrichment) for a gene to be taken as specific for that cluster
  PCutoff = NULL, # Set a maximum pvalue for the enrichment of a gene to be considered significant.
  upDown = 'Up', # Specify wether to look only for upregulated ('Up'), downregulated ('Down') or both ('Both') genes.
  pathwayPCutoff = 0.05, # Specify the maximum p-value for a pathway to be considered significant
  GOKegg='GO', # Specify whether to do GO ('GO') or KEGG ('KEGG') analysis
  includeChildTerms=F # Only matters for GO. If T it counts hits for a certain term also als a hit for all parent terms.
)
# Plot gene ontology for clusters
plotClusterGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  clusterGO, # Specify the variable that has the cluster gene ontology results
  clusterDifferentialExpression, # Supply the variable name of the variable that contains the clusterdifferentialexpression data generated by analyseClusterDifferentialExpression (done further up in this script)
  topX=10, # For plotting, specify the topX number of pathways that should be plotted.
  plotRegulatedGenes = F, # Should the function plot t-SNE maps of all the regulated genes in each regulated pathway (T if you want this, F if you don't.). WARNING: MAY TAKE A LOT OF TIME TO PLOT
  overrideDir = F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode = 'pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)



# Analyze gene ontology for simple correlation analysis
TTNposGO <- analyseCorrelationGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  correlationResultTTN, # Specify the variable that contains the result of your correlation analysis
  minExpression=0.2, # Specify the minimum epxression level that is required to consider a gene expressed in a cell
  corrCutoff=0.2, # Specify the correlation cutoff value for correlation to be significant
  pCutoff=1/10^5, # Specify the pvalue cutoff value for correlation to be significant
  posNeg='Pos', # Specify which correlations should be looked at, only the positive ones ('Pos'), only the negative ones ('Neg') or both ('Both')
  pathwayPCutoff=0.05, # Specify the maximum p-value for a pathway to be considered significant
  GOKegg='GO', # Specify whether to do GO ('GO') or KEGG ('KEGG') analysis
  includeChildTerms=F # Only matters for GO. If T it counts hits for a certain term also als a hit for all parent terms.
)
# Plot gene ontology for simple correlation analysis
plotCorrelationGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS, # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  GOResults = TTNposGO, # Specify the variable that has the correlation gene ontology results
  corrResults=correlationResultTTN, # Specify the variable that has the correlation analysis results
  topX=10, # For plotting, specify the topX number of pathways that should be plotted.
  plotRegulatedGenes=FALSE, # Should the function plot t-SNE maps of all the regulated genes in each regulated pathway (T if you want this, F if you don't.). WARNING: MAY TAKE A LOT OF TIME TO PLOT
  overrideDir=F, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode = 'pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
) 



# Analyze gene ontology for regulons
regulonGO <- analyzeRegulonGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  correlationMatrices, # Specify the variable that contains the correlation matrices
  regulons, # Specify the variable that contains your regulons
  pathwayPCutoff=0.05, # Specify the maximum p-value for a pathway to be considered significant
  GOKegg='GO', # Specify whether to do GO ('GO') or KEGG ('KEGG') analysis
  includeChildTerms=F # Only matters for GO. If T it counts hits for a certain term also als a hit for all parent terms.
)
# Plot gene ontology for regulons
plotRegulonGeneOntology(
  config, # Don't change this. It needs the configuration variables you set up in the configuration part of the script.
  groupedSCS,  # This needs to exactly match the variable name of the variable that has the the generated SCS models.
  regulonGO, # Specify the variable that has the regulon gene ontology results
  regulons, # Specify the variable that has the regulons
  topX=10, # For plotting, specify the topX number of pathways that should be plotted.
  plotRegulatedGenes=FALSE, # Should the function plot t-SNE maps of all the regulated genes in each regulated pathway (T if you want this, F if you don't.). WARNING: MAY TAKE A LOT OF TIME TO PLOT
  overrideDir=FALSE, # Override the standard saving directory for the generated graphs. Can be used to plot graphs in figure-specific folders like 'Fig2/B'. Leave as FALSE to plot in standard directories.
  outputMode = 'pdf' # Use 'show' to show the graph in RStudio. Change to 'pdf' or 'png' to save the graphs.
)
