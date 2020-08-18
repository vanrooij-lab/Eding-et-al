### plotGeneExpressionCutoffTSNE()
### 
### Dichotomizes gene expression based on a provided cutoff. Outputs a black/white t-SNE
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
plotGeneExpressionCutoffTSNE <- function(config, groupedSCS, geneNames=NULL, overUnder='Under', cutoff=NULL, groupNames=NULL, includeOutlierCells=T, colorScheme='Spectral', pointSizeModifier = 1, overrideDir=F, overrideTheme = NULL, outputMode='pdf') {
  # Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }
  
  #Sanitize input includeOutlierCells
  if(!is.logical(includeOutlierCells)) {stop('Parameter includeOutlierCells must be either TRUE or FALSE. Pay attention to NOT using quote')}

  if(!(overUnder %in% c('Over','Under'))) {
    stop("overUnder needs to be either 'Over' or 'Under'.")
  }
  
  # Sanitize input cutoff
  if(!is.numeric(cutoff)) {
    stop("Please input a number for parameter 'cutoff'.")
  } else {
    cutoff = as.double(cutoff)
  }
  
  # Sanitize input geneNames
  if(is.null(geneNames)) {stop("You did not provide (a) gene name(s) to draw (a) t-SNE map(s) for")}
  for(geneName in geneNames) {
    for(groupName in groupNames) {
      # Find geneName in correct format
      newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@ndata), value=T)
      
      if(length(newGeneName) == 0) {
        # Gene not found
        stop(paste0("You want to analyze marker gene ",geneName," which is not expressed in group ",groupName))
      } else if(length(newGeneName) > 1) {
        # Gene name matches multiple
        stop(paste0("You want to generate a t-SNE for gene ",geneName," which is ambiguous. Did you mean: ", paste(newGeneName, collapse = " OR "),"?"))
      } else {
        # Gene is unique, update identifier
        geneNames[match(geneName, geneNames)] = newGeneName
      }
    }
  }

  # Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  # Iterate over provided groups
  for(groupName in groupNames) {
    # Create directory if necessary
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
    
    # Iterate over provided geneNames
    for(geneName in geneNames) {
      #Prepare dataframe for ggplot
      cellsForTSNE = data.frame(
        cellNames = names(groupedSCS[[groupName]]@cpart),
        clusterNum = as.factor(groupedSCS[[groupName]]@cpart),
        V1 = groupedSCS[[groupName]]@tsne$V1,
        V2 = groupedSCS[[groupName]]@tsne$V2,
        expression = as.double(groupedSCS[[groupName]]@ndata[geneName,])
      )
      
      # If cells from outlier clusters should be removed, do so now
      if(!includeOutlierCells) {
        cellsForTSNE = cellsForTSNE[which(cellsForTSNE$clusterNum %in% majorClusterNums),]
      } 
      
      if(overUnder == 'Over') {
        cellsForTSNE$status = 'Not over'
        cellsForTSNE[which(cellsForTSNE$expression > cutoff),'status'] = 'Over'
      } else if(overUnder == 'Under') {
        cellsForTSNE$status = 'Not under'
        cellsForTSNE[which(cellsForTSNE$expression < cutoff),'status'] = 'Under'
      }
      
      # Plot graph
      expressionTSNE = ggplot(
        data=cellsForTSNE
      ) + geom_point(
        aes(
          x=V1, 
          y=V2,
          color=status
          ),
        size=1.5*pointSizeModifier
      ) + theme_classic(
      ) + theme(
        aspect.ratio = config$plotOptions$aspect.ratio,
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()
      ) + labs(
        x=element_blank(),
        y=element_blank(),
        title=paste0(groupName," - t-SNE - ",strsplit(geneName,'__')[[1]][1]," - cutoff ",cutoff),
        subtitle = paste0('Cutoff: ',cutoff)
      ) + overrideTheme
      if(outputMode == 'show') {
        print(expressionTSNE)
      } else {
        fileName = paste0(groupName,'_',strsplit(geneName,'__')[[1]][1],'_cutoff',cutoff,'.',outputMode)
        ggsave(fileName, path=saveDir, useDingbats=FALSE)
        print(paste0("t-SNE map saved as ",saveDir,.Platform$file.sep,fileName))
      }
    }
  }
}