### analyseClusterComposition <- function(config, groupedSCS, groupNames=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does analysis for all.
### outputMode: Pass 'save' (saves output to Excel file). Other values won't be accepted.
###
### Return: void. Saves a set of PDF files with t-SNE maps (1 per cluster) coloured according to sample origin.
###               Additionally saves an excel file containing the sample contribution to each cluster
### 
### Authors
###   - Joep Eding
### 
### TODO
###   -
analyseClusterComposition <- function(config, groupedSCS, groupNames=NULL, overrideDir=F, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop("You want to analyse cluster size for SCSObject '",groupName,"' that is not in groupedSCS")}
    }
  }
  
  #Check outputMode is set to save. Not saving makes no sense for this function.
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }

  #Iterate over groups
  for(groupName in groupNames) {
    #Create directory if necessary
    if(outputMode != 'show') {
      if(overrideDir == FALSE) {
        saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Clusters')
      } else {
        saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
      }
      if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Get TSNE-location for all cells in this group
    tsneLocation = groupedSCS[[groupName]]@tsne
    
    #Find out how much clusters
    numClusters = length(groupedSCS[[groupName]]@cluster$jaccard)
    
    #Make a dataframe to save statistics in
    stats <- data.frame(
      'clusterNum' = seq(1, numClusters)
    )
    for(sampleName in unlist(config$groups[[groupName]])) {stats[,sampleName] = 0}
    
    #For each cluster, plot new cell location on map of all cells
    for(i in 1:numClusters) {
      #Select names of all cells belonging to this cluster
      clusterCells <- names(groupedSCS[[groupName]]@cluster$kpart[which(groupedSCS[[groupName]]@cluster$kpart == i)])
      
      #Create a copy of tsneLocation, add a column for the grouping
      tsneGrouping <- tsneLocation
      tsneGrouping['sample'] = 'nonCluster'
      
      #Add sample origin for all cells in the cluster being analysed
      for(cellName in rownames(tsneGrouping)) {
        #Check if sample is in cluster, skip if it is not.
        if(!(cellName %in% clusterCells)) {next}
        
        #Iterate over samples to find sample this cell belongs to
        for(sampleName in names(config$samples)) {
          if(strsplit(cellName, '_')[[1]][1] %in% config$samples[[sampleName]]) {
            tsneGrouping[cellName, 'sample'] = sampleName
            break
          }
        }
      }
      
      #Prepare a named color for the cells that don't get a sample color. 
      nonClusterColor = c('#000000')
      names(nonClusterColor) <- c('nonCluster')
      
      #Generate graph
      clusterComposition <- ggplot(   #Start ggplot with data and grouping
        tsneGrouping,
        aes(x=V1, y=V2, group=sample)
      ) + geom_point(    #Print points
        aes(color=sample)
      ) + scale_colour_manual(    #Define coloring per sample based on colors from config$colors
        values=unlist(append(nonClusterColor, config$colors))
      ) + theme_classic(         #Get rid of default ggplot theme
      ) + labs(
        x='', y='',     #Get rid axis labels
        color='Patient',
        title = paste0(groupName,' sample contribution for cluster ',i)
      )
      
      #Show or save graph
      if(outputMode == 'show') {
        print(clusterComposition)
      } else {
        fileName = paste0(groupName,'_clusterComposition_cluster',i,'.',outputMode)
        ggsave(fileName, path=saveDir)
        print(paste0("t-SNE map saved as ",saveDir,.Platform$file.sep,fileName))
      }
      
      #Add data to stats variable
      for(sampleName in unlist(config$groups[[groupName]])) {
        stats[which(stats$clusterNum==i),sampleName] = nrow(tsneGrouping[which(tsneGrouping$sample==sampleName),])
      }
    }
    
    #Store 'stats' in Excel workbook
    statsWB <- createWorkbook()
    addWorksheet(statsWB, sheetName='Sample contribution to cluster')
    writeDataTable(
      statsWB,
      sheet='Sample contribution to cluster',
      stats
    )
    if(outputMode != 'show') {
      saveWorkbook(statsWB, paste0(saveDir,.Platform$file.sep,'sampleContribution.xlsx'), overwrite = TRUE) 
    }
  }
}