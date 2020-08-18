### analyseClusterTyping()
### Attempts to make a graph like Figure 6D in this paper: http://dev.biologists.org.proxy.library.uu.nl/content/145/18/dev168609
### 
### Parameters
###   -
### 
### Return
###   -
###   
### Authors
###   - Joep Eding
### 
### TODO
###   - 
analyseClusterTyping <- function(config, groupedSCS, clusterDifferentialExpression, geneNames=NULL, groupNames=NULL, numGenesPerCluster=NULL, normalize='no', displayLog='yes', overrideDir=F, outputMode='pdf') {
  # Sanitize input normalize & displayLog
  if(!(normalize %in% c('cluster', 'gene', 'no'))) {stop("Value for 'normalize' can only be 'cluster', 'gene', or 'no'.")}
  if(displayLog=='yes') {
    displayLog=T
  } else if(displayLog=='no') {
    displayLog=F
  } else {
    stop("Value of 'displayLog' must be either 'yes' or 'no'.")
  }
  
  # Sanitize input groupNames & input referenceGroup
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to compare clusters across groups for '",groupName,"' that is not in groupedSCS."))}
    }
  }
  
  # Sanitize input geneNames
  if(is.null(geneNames) & is.null(numGenesPerCluster)) {stop("Specify a value for at least one of 'geneNames' or 'numGenesPerCluster'")}
  #If no geneNames are set (but numGenesPerCluster is provided), automatically determine top enriched genes per cluster to use as geneNames
  if(is.null(geneNames)) {
    if(!is.numeric(numGenesPerCluster)) {stop("'numGenesPerCluster' needs to be a number")} else {numGenesPerCluster=as.integer(numGenesPerCluster)}
    geneNamesTemp = returnTopDifferentialGenesPerCluster(config, groupedSCS, clusterDifferentialExpression, numGenes = numGenesPerCluster)
    geneNames = list()
    for(groupName in groupNames) {
      geneNames[[groupName]] = list()
      for(clusterNum in 1:length(geneNamesTemp[[groupName]])) {
        geneNames[[groupName]] = unlist(append(geneNames[[groupName]], geneNamesTemp[[groupName]][[clusterNum]][1:numGenesPerCluster]))
      }
      geneNames[[groupName]] = unique(geneNames[[groupName]])
    }
  } else {
    #Check whether groupNames are already defined (then marker genes should be specified per group)
    if(any(names(geneNames) %in% groupNames)) {
      #Check if all groups for which marker genes are specified are also in groupNames
      for(groupName in names(geneNames)) {
        if(!(groupName %in% groupNames)) {
          stop("You specified marker genes for group '",groupName,"' which is not in groupNames")
        }
      }
      
      #Check if all groups in groupNames have marker genes specified
      for(groupName in groupNames) {
        if(!(groupName %in% names(geneNames))) {
          stop("You want to analyze cluster typing for group '",groupName," but have not specified marker genes for this group")
        }
      }
    } else { #If no groups are specified, assume there's one list of markerGenes to use for all groups
      #Check that list is not empty
      if(length(geneNames) < 1) {error("No genes specified")}
      
      #Throw warning to notify user that this assumption is being made
      warning("No group distinction made in geneNames, the same list of genes will be used for every group", immediate. = T)
      
      #Reformat geneList to fit for every group
      geneNamesTemp = geneNames
      geneNames = list()
      for(groupName in groupNames) {
        geneNames[[groupName]] = geneNamesTemp
      }
    }
  }
  #Now check whether all geneNames are possible
  for(groupName in names(geneNames)) {
    for(geneName in geneNames[[groupName]]) {
      #Find geneName in correct format
      newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@ndata), value=T)
      
      if(length(newGeneName) == 0) {
        #Gene not found
        stop(paste0("You want to analyze marker gene ",geneName," which is not expressed in group ",groupName))
      } else if(length(newGeneName) > 1) {
        #Gene name matches multiple
        stop(paste0("You want to analyze marker gene ",geneName," which is ambiguous. Did you mean: ", paste(newGeneName, collapse = " OR "),"?"))
      } else {
        #Gene is unique, update identifier
        geneNames[[groupName]][match(geneName, geneNames[[groupName]])] = newGeneName
      }
    }
  }
  
  # Sanitize input outputMode
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
    
    #Further process data
    if(any(grepl('graphingData',ls()))) {remove('graphingData')} #Remove previous versions of graphingData to clear factor levels
    graphingData = data.frame(
      'cluster' = character(),
      'gene' = character(),
      'mean' = double(),
      'percentage' = double()
    )
    
    #Extract clustering information
    clustering = groupedSCS[[groupName]]@cpart
    
    #Fill graphingData by cluster and by gene
    for(clusterNum in 1:max(clustering)) {
      for(geneName in geneNames[[groupName]]) {
        #Select data
        clusterGeneData = groupedSCS[[groupName]]@ndata[geneName, names(clustering[which(clustering == clusterNum)])]
        
        #Then process
        if(length(clusterGeneData) == 1) {
          meanExpression=clusterGeneData
          if(meanExpression==0.1) {percentageExpression=0} else {percentageExpression=100};
        } else {
          meanExpression = rowMeans(clusterGeneData[,which(clusterGeneData[geneName,] >0.0)])[[1]]
          percentageExpression = length(which(clusterGeneData[geneName,] > 0.1))/ncol(clusterGeneData)*100
        }
        
        #Then add to graphingData
        rowNum = nrow(graphingData) + 1
        graphingData[rowNum, 'cluster'] = clusterNum
        graphingData[rowNum, 'gene'] = geneName
        graphingData[rowNum, 'mean'] = meanExpression
        graphingData[rowNum, 'percentage'] = percentageExpression
        
      }
    }
    
    #Normalize expression
    if(normalize=='cluster') {
      for(clusterNum in 1:max(clustering)) {
        graphingData[which(graphingData$cluster==clusterNum),]$mean = scale(graphingData[which(graphingData$cluster==clusterNum),]$mean)
      }
      normalizationTag="Normalized expression within cluster"
    } else if(normalize=='gene') {
      for(geneName in geneNames[[groupName]]) {
        graphingData[which(graphingData$gene==geneName),]$mean = scale(graphingData[which(graphingData$gene==geneName),]$mean)
      }
      normalizationTag="Normalized expression per gene"
    } else if(normalize=='no') {
      normalizationTag="No normalization of expression values"
    }

    
    #Make sure X-axis order & Y-axis order are preserverd
    graphingData$gene <- factor(graphingData$gene, levels=geneNames[[groupName]])
    graphingData$cluster <- factor(graphingData$cluster, levels=lapply(1:max(clustering), function(x) {as.character(x)}))
    
    xAxisLabels = unlist(lapply(graphingData$gene, function(x) {strsplit(as.character(x), '__')[[1]][1]}))
    names(xAxisLabels) = graphingData$gene
    
    clusterTyping <- ggplot(
      graphingData, 
      aes(y=cluster, x=gene)
    ) + geom_count(
      aes(
        size=percentage, 
        color=if(displayLog) log(mean) else mean
        )
    ) + scale_color_distiller(
      palette="Spectral"
    ) + scale_x_discrete(
      labels=xAxisLabels
    ) + theme_classic(
    ) + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)
    ) + labs(
      x="Genes",
      y="Clusters",
      title=paste0(groupName," - Marker genes per cluster"),
      subtitle=normalizationTag,
      color=if(displayLog) 'log(mean expression)' else 'mean expression',
      size="% expression"
    )
    
    if(outputMode != 'show') {
      ggsave(paste0(groupName,'_clusterTyping_',normalize,'Normalized.',outputMode), path=saveDir)
    } else {
      print(clusterTyping)
    }
  }
}