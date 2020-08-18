### compareClustersAcrossGroups()
###
### TODO Upgrade this comment
### TODO Consider renaming this function
### TODO Consider splitting this function into 3 functions
### 
### Authors
###   - Joep Eding

### Return: 
### --- Outputs 3 graphs:
###    --- clusterOriginBarplot_clusterResolution: Shows which sample clusters the clusters in referenceGroup originate from.
###    --- clusterDestinationBarplot: Shows which referenceGroup clusters each sampleCluster contributes to.
###    --- clusterOriginBarplot_sampleResolution: Similar to the clusterResolution version, now shows only from which samples the clusters originate.
### --- Outputs 1 excel file: (only when outputMode is not 'show')
###    --- clusterComparison.xlsx: contains the raw data the graphs are based upon.
compareClustersAcrossGroups <- function(config, groupedSCS, referenceGroup=NULL, groupNames=NULL, orientation='vertical', includeOutliers=TRUE, overrideDir=F, outputMode='pdf') {
  # Sanitize input groupNames & input referenceGroup
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to compare clusters across groups for '",groupName,"' that is not in groupedSCS."))}
    }
  }
  if(is.null(referenceGroup)) {stop(paste0("You want to compare clusters accross groups without setting a reference group."))}
  if(!(referenceGroup %in% groupNames)) {stop(paste0("You want to compare clusters across groups for reference group '", referenceGroup,"' that is not in groupNames."))}
  if(length(referenceGroup) > 1) {stop("You have selected more than 1 referenceGroup.")}
  if(all(groupNames %in% referenceGroup)) {stop("You have not selected any other groups than the referenceGroup.")}
  
  # Sanitize input orientation
  if(!(orientation %in% c('horizontal','vertical'))) {stop("Parameter 'orientation' needs to be one of: 'horizontal', 'vertical'")}
  
  # Sanitize input includeOutliers
  if(!is.logical(includeOutliers)) {stop("Parameter 'includeOutliers' needs to be either TRUE or FALSE. Pay attention: no quotes!")}
  
  # Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }

  # Generate list of all clusters  
  clusterList = list()
  for(groupName in groupNames) {
    for(clusterName in names(clusterDifferentialExpression[[groupName]])) {
      #Check if whether this isn't an outlierCluster
      if(includeOutliers) {
        if(!(as.integer(strsplit(clusterName, 'cl.')[[1]][2]) %in% 1:max(groupedSCS[[groupName]]@cpart))) { next }
      } else {
        if(!(as.integer(strsplit(clusterName, 'cl.')[[1]][2]) %in% 1:max(groupedSCS[[groupName]]@cluster$kpart))) { next }
      }
      
      
      newClusterName = paste0(groupName,'_',clusterName)
      clusterList[[newClusterName]] = list()
      clusterList[[newClusterName]]['groupName'] = groupName
      clusterList[[newClusterName]]['clusterNum'] = clusterName
    }
  }
  
  # Calculate comparisonTable: contains the number of overlapping cells between all of the clusters (compares individual samples vs groups too)
  comparisonTable <- data.frame()
  for(clusterAName in names(clusterList)) {
    groupAName = clusterList[[clusterAName]]$groupName
    clusterANum = clusterList[[clusterAName]]$clusterNum
    
    if(includeOutliers) {
      clusterAcells <- names(groupedSCS[[groupAName]]@cpart[groupedSCS[[groupAName]]@cpart==as.integer(strsplit(clusterANum, 'cl.')[[1]][2])])
    } else {
      clusterAcells <- names(groupedSCS[[groupAName]]@cluster$kpart[groupedSCS[[groupAName]]@cluster$kpart==as.integer(strsplit(clusterANum, 'cl.')[[1]][2])])  
    }
    
    
    for(clusterBName in names(clusterList)) {
      groupBName = clusterList[[clusterBName]]$groupName
      clusterBNum = clusterList[[clusterBName]]$clusterNum
      
      if(includeOutliers) {
        clusterBcells <- names(groupedSCS[[groupBName]]@cpart[groupedSCS[[groupBName]]@cpart==as.integer(strsplit(clusterBNum, 'cl.')[[1]][2])])  
      } else {
        clusterBcells <- names(groupedSCS[[groupBName]]@cluster$kpart[groupedSCS[[groupBName]]@cluster$kpart==as.integer(strsplit(clusterBNum, 'cl.')[[1]][2])])
      }
      
      
      comparisonTable[clusterAName,clusterBName] = length(which(clusterAcells %in% clusterBcells))
    }
  }
  
  if(outputMode != 'show') {
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,referenceGroup,.Platform$file.sep,'Clusters')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
    }
    # Create directory if necessary
    if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
  }
  
  # Generate source data for the stacked bar graph by melting the comparisonTable
  stackedBarSource = melt(as.matrix(comparisonTable)) #Uses as.matrix to keep the rowNames from 
  colnames(stackedBarSource) <- c('destinationCluster','sourceCluster','number')
  # Remove rows that don't have a cluster in the referenceGroup as destinationCluster
  stackedBarSource = stackedBarSource[which(grepl(referenceGroup,stackedBarSource$destinationCluster)),]
  # Remove rows that have a cluster in the referenceGroup as sourceCluster
  stackedBarSource = stackedBarSource[which(!grepl(referenceGroup,stackedBarSource$sourceCluster)),]
  
  
  # Generate lists of all cluster names per sample
  sourceGroupNames = names(config$groups)[which(names(config$groups) != referenceGroup)]
  sourceGroups = list()
  for(sourceName in unique(stackedBarSource$sourceCluster)) {
    sourceGroup = strsplit(sourceName, '_cl.')[[1]][1]
    sourceNum = strsplit(sourceName, '_cl.')[[1]][2]
    
    sourceGroups[[sourceGroup]] = append(sourceName, sourceGroups[[sourceGroup]])
  }
  
  # Generate different colors for each sample, and then variations (lighter/darker) on that color for cluster in that sample
  colorList = list()
  for(sampleNum in 1:length(sourceGroups)) {
    #Count number of clusters for this sample
    numColors = length(sourceGroups[[sampleNum]])
    
    #Test whether a colour is specified for this group
    if(!(nchar(config$colors[names(sourceGroups)[sampleNum]]) > 0)) {
      stop("No colour is provided for sample '",names(sourceGroups)[sampleNum],"'")
    }
    
    #Genereate colors
    colorList[[sampleNum]] = colorRampPalette(c(
      config$colors[names(sourceGroups)[sampleNum]],
      #rainbow(length(sourceGroups))[sampleNum], #Generate distinct colors for each sample
      '#000000'))(numColors*2)[1:numColors] #Generate a colorRampPalette ranging from the sample color to white
  }
  
  
  colorList = unlist(colorList)
  # Name colors so they mach the clusters
  names(colorList) <- unlist(sourceGroups)
  
  #Draw a plot that shows which clusters the cells in each cluster of referenceGroup came from. 
  clusterOriginGraph = ggplot(
    data=stackedBarSource[which(stackedBarSource$number>5),],
    aes(x=destinationCluster, y=number, fill=sourceCluster)
  ) + geom_bar(
    stat='identity'
  ) + labs(
    x="Cluster",
    y="Cells (#)",
    title="Cluster origin",
    subtitle="Shows from which sample-specific clusters the final clusters derive.\nExcludes contributions <5 cells per sample cluster"
  ) + scale_fill_manual(    #Define coloring per sample based on colors from config$colors
    values=colorList
  ) + geom_text(
    aes(label=lapply(sourceCluster, function(x){strsplit(as.character(x),'_')[[1]][2]})),
    position = position_stack(vjust = 0.5),
    angle=if(orientation=='horizontal') 90 else 0, #Tilt x-axis labels only when numbers go above 9
    size=5/.pt
  ) + scale_x_discrete(
    breaks=unique(stackedBarSource$destinationCluster), 
    labels=lapply(unique(stackedBarSource$destinationCluster), function(x){as.character(strsplit(as.character(x),'_cl.')[[1]][2])})
  ) + theme_classic(
  ) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank(),
    legend.title = element_text(size=rel(0.75)),
    legend.text = element_text(size=rel(0.5)),
    plot.subtitle = element_text(size=rel(0.5))
  )
  if(orientation == 'horizontal') {
    clusterOriginGraph = clusterOriginGraph+coord_flip()
  }
  if(outputMode != 'show') {
    ggsave(paste0('clusterOriginBarplot_clusterResolution.',outputMode), path=saveDir)
  } else {
    print(clusterOriginGraph)
  }
  
  #Draw a plot that shows which cluster in referenceGroup each of the cells from all group-specific cluster go to.
  clusterDestinationGraph = ggplot(
    data=stackedBarSource[which(stackedBarSource$number>5),],
    aes(x=sourceCluster, y=number, fill=as.character(lapply(destinationCluster, function(x){strsplit(as.character(x),'_cl.')[[1]][2]})))
  ) + geom_bar(
    stat='identity'
  ) + labs(
    x="Sample clusters",
    y="Cells (#)",
    fill="Destination cluster",
    title="Cluster destination",
    subtitle="Shows which final clusters the sample-specific clusters contribute to.\nExcludes contributions <5 cells"
  ) + geom_text(
    aes(label=lapply(destinationCluster, function(x){strsplit(as.character(x),'_cl.')[[1]][2]})),
    position = position_stack(vjust = 0.5)
  ) + theme_classic(
  ) + theme(
    plot.subtitle = element_text(size=rel(0.5)),
    axis.text.x = element_text(
      angle = 90, 
      vjust = 0.5,
      hjust = 0
    )
  )
  if(orientation == 'horizontal') {
    clusterDestinationGraph = clusterDestinationGraph+coord_flip()
  }
  if(outputMode != 'show') {
    ggsave(paste0('clusterDestinationBarplot.',outputMode), path=saveDir)
  } else {
    print(clusterDestinationGraph)
  }
  
  #Generate new dataframe for a sample-resolution clusterOriginBarPlot
  clusterOriginSource = data.frame(
    destinationCluster = grep(referenceGroup,names(clusterList),value=T),
    stringsAsFactors = F
  )
  rownames(clusterOriginSource) <- grep(referenceGroup,names(clusterList),value=T)
  for(destinationCluster in clusterOriginSource$destinationCluster) {
    for(sourceSample in names(sourceGroups)) {
      clusterOriginSource[destinationCluster,sourceSample] = sum(comparisonTable[destinationCluster,grepl(sourceSample,colnames(comparisonTable))])
    }
  }
  clusterOriginSource <- melt(clusterOriginSource,id.vars='destinationCluster')
  colnames(clusterOriginSource) <- c('destinationCluster','sourceSample','number')
  clusterOriginSource$destinationCluster <- unlist(lapply(clusterOriginSource$destinationCluster, function(x){as.character(str_pad(strsplit(x,'_cl.')[[1]][2],2,side='left',pad=' '))}))
  #clusterOriginSource$destinationCluster <- as.integer(clusterOriginSource$destinationCluster)
  
  #Draw a plot that shows which samples the cells in each cluster of referenceGroup came from. 
  clusterOriginGraph = ggplot(
    data=clusterOriginSource,
    aes(x=destinationCluster, y=number, fill=sourceSample)
  ) + geom_bar(
    stat='identity'
  ) + labs(
    x="Cluster",
    y="Cells (#)",
    fill="Sample",
    title="Cluster origin",
    subtitle="Shows from which samples the final clusters derive."
  ) + scale_fill_manual(    #Define coloring per sample based on colors from config$colors
    values=unlist(config$colors)
  ) + theme_classic(
  ) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank(),
    plot.subtitle = element_text(size=rel(0.5))
  )
  if(orientation == 'horizontal') {
    clusterOriginGraph = clusterOriginGraph+coord_flip()
  }
  if(outputMode != 'show') {
    ggsave(paste0('clusterOriginBarplot_sampleResolution.',outputMode), path=saveDir)
  } else {
    print(clusterOriginGraph)
  }
  
  clusterComparisonWB = createWorkbook()
  #Write comparisonTable to the new sheet
  addWorksheet(clusterComparisonWB, sheetName='comparisonTable')
  writeDataTable(
    clusterComparisonWB,
    sheet = 'comparisonTable',
    comparisonTable
  )
  #Write stackedBarSource to the new sheet
  addWorksheet(clusterComparisonWB, sheetName = 'stackedBarSource')
  writeDataTable(
    clusterComparisonWB,
    sheet='stackedBarSource',
    stackedBarSource
  )
  #Save workbook
  if(outputMode != 'show') {
    saveWorkbook(clusterComparisonWB, paste0(saveDir,.Platform$file.sep,'clusterComparison.xlsx'))
  }
}