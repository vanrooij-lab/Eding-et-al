### plotClusterTSNE(config, groupedSCS, groupNames=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### groupNames: Optional. Pass the name of a single group to plot a clustering tSNE for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does analysis for all groups in groupedSCS.
### outputMode: Save will save the graphs to the tSNE-directory. Show will just print the graph. (default: 'save')  
### 
### Authors
###   - Joep Eding
### 
### TODO: Update comment
### TODO: Change filename and/or subtitle based on settings for outlierExclusion?
### TODO: Add checks for pointSizeModifier
### TODO: Add the viridis color package to colorScheme as well (or: make gene expression and clustering TSNE's all use the same color picking logic?)
### 
### Return: void. Just saves graphs or prints them to viewer.          
plotClusterTSNE <- function(config, groupedSCS, groupNames=NULL, includeOutlierClusters=T, includeOutlierCells=T, colorScheme='Spectral', pointSizeModifier=1, overrideDir=F, overrideTheme = NULL, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }
    
  #Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  #Sanitize input includeOutlierClusters
  if(!is.logical(includeOutlierClusters)) {stop('Parameter includeOutlierClusters must be either TRUE or FALSE. Pay attention to NOT using quote')}
  
  #Sanitize input includeOutlierCells
  if(!is.logical(includeOutlierCells)) {stop('Parameter includeOutlierCells must be either TRUE or FALSE. Pay attention to NOT using quote')}
  if(includeOutlierClusters && !includeOutlierCells) {stop("Can't include outlierClusters but exclude outlierCells, check your settings")}
  
  #Process input colors
  if(is.null(colorScheme) || !is.character(colorScheme)) {
    stop("Parameter colorScheme needs a character input compatible with the ColorBrewer package.")
  } else if(length(colorScheme) == 1) {
    # TODO Implement a color check for existing palettes here.
    tSNEColoring = scale_color_distiller(
      palette=colorScheme
    )
  } else {
    stop("Parameter colorScheme needs a character input compatible with the ColorBrewer package.")
  } 
  
  #Iterate over provided groups
  for(groupName in groupNames) {
    #Create directory if graphs need to be saved
    if(outputMode != 'show') {
      if(overrideDir == FALSE) {
          saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'t-SNE')
      } else {
        saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
      }
      if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Figure out which clusterNumbers are 'major' and which are 'outlier' clusters
    majorClusterNums = 1:max(groupedSCS[[groupName]]@cluster$kpart)
    outlierClusterNums = 1:max(groupedSCS[[groupName]]@cpart)
    outlierClusterNums = outlierClusterNums[which(!(outlierClusterNums %in% majorClusterNums))]
    
    #Prepare colors for the major clusters
    TSNEMainColors = colorRampPalette(brewer.pal(9, colorScheme))(max(majorClusterNums))
    #Prepare colors for the outlier clusters
    TSNEOutlierColors = colorRampPalette(brewer.pal(9, colorScheme))(length(outlierClusterNums))
    #Now repeat the same colors for the outlier clusters (which will get a different shape)
    TSNEColors = append(TSNEMainColors, TSNEOutlierColors)
    
    #Determine shapes for the clusters
    TSNEShapes = append(rep(16, length(majorClusterNums)), rep(21, length(outlierClusterNums)))
    
    #Prepare dataframe for ggplot
    cellsForTSNE = data.frame(
      cellNames = names(groupedSCS[[groupName]]@cpart),
      clusterNum = as.factor(groupedSCS[[groupName]]@cpart),
      V1 = groupedSCS[[groupName]]@tsne$V1,
      V2 = groupedSCS[[groupName]]@tsne$V2
    )
    
    # Clean up the dataframe according to outlierSettings
    if(includeOutlierClusters && includeOutlierCells) {
      # No need to do anything, dataframe is great as is.
    } else if(includeOutlierClusters && !includeOutlierCells) {
      # This makes no sense, but also can't happen since it's caught at the start of the function
    } else if(!includeOutlierClusters && !includeOutlierCells) {
      # Don't include the outlier clusters, don't include the outlier cells.
      cellsForTSNE = cellsForTSNE[which(cellsForTSNE$clusterNum %in% majorClusterNums),]
    } else if(!includeOutlierClusters && includeOutlierCells) {
      # Don't inlude the outlier clusters but do include the cells.
      # Cells in minor clusters need to be replaced to their original clusters
      for(cellName in names(groupedSCS[[groupName]]@cpart[which(groupedSCS[[groupName]]@cpart %in% outlierClusterNums)])) {
        cellsForTSNE[match(cellName, cellsForTSNE$cellNames),'clusterNum'] = groupedSCS[[groupName]]@cluster$kpart[[cellName]]
      }
    }
    
    #Plot graph
    clusterTSNE <- ggplot(
      data=cellsForTSNE,
      aes(x=V1,y=V2)
    ) + geom_point(
      aes(color=clusterNum, shape=clusterNum),
      size=1.5*pointSizeModifier
    ) + scale_color_manual(
      values=TSNEColors
    ) + scale_shape_manual(
      values=TSNEShapes
    ) + theme(
      aspect.ratio = 1
    ) + theme_classic(
    ) + theme(
      aspect.ratio = config$plotOptions$aspect.ratio,
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.line = element_blank()
    ) + labs(
      x=element_blank(),
      y=element_blank(),
      color='Cluster',
      shape='Cluster',
      title=paste0(groupName," - t-SNE - Clustering")
    ) + overrideTheme
    if(outputMode != 'show') {
      ggsave(paste0(groupName,'_clustering.',outputMode),path = saveDir, useDingbats=F)
    } else {
      print(clusterTSNE)
    }
  }
}