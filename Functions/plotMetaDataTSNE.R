### plotMetaDataTSNE(config, groupedSCS, genenames=NULL, groupNames=NULL, outputMode='save')
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
plotMetaDataTSNE <- function(config, groupedSCS, metadataColumns=NULL, groupNames=NULL, colorScheme='Spectral', expressionScale='linear', pointSizeModifier=1, overrideDir=F, overrideTheme = NULL, outputMode='pdf') {
  # Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }

  # Check whether a metadata object exists
  if(!exists('scsMetadata')) {
    stop("Load metadata first by running loadMetaData().")
  }
  
  # Sanitize input metaDataColumns
  if(is.null(metadataColumns)) {stop("You did not provide (a) gene name(s) to draw (a) t-SNE map(s) for")}
  for(metadataName in metadataColumns) {
    if(!(metadataName %in% colnames(scsMetadata))) {
      print(paste0("Metadata variable '",metadataName,"' not loaded, skipping."))
      next() # This achieves nothing. 
    }
  }
  
  # Sanitize input expressionScale
  if(!(expressionScale %in% c('linear','log','log2','log10'))) {stop(paste0("Value '",expressionScale,"' is not an accepted value for expressionScale. Use 'linear' (default), 'log', 'log2' or 'log10'"))}
    
  # Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  # Process input colors
  if(is.null(colorScheme)) {
    tSNEColoring = NULL
  } else if(length(colorScheme) == 1) {
    # TODO Implement a color check for existing palettes here.
    tSNEColoring = scale_color_distiller(
      palette=colorScheme
    )
  } else if(length(colorScheme) == 2) {
    tSNEColoring = scale_color_gradient(
      low=colorScheme[1],
      high=colorScheme[2],
      guide=guide_colorbar(barheight=15)
    )
  } else if(length(colorScheme) == 3) {
    tSNEColoring = scale_color_gradient2(
      low=colorScheme[1],
      mid=colorScheme[2],
      high=colorScheme[3],
      guide=guide_colorbar(barheight=15)
    )
  } else if(length(colorScheme) > 3) {
    tSNEColoring = scale_color_gradient(
      colours=colorScheme[1],
      guide=guide_colorbar(barheight=15)
    )
  } else {
    stop("Did not understand your input to 'colorScheme', needs to be a list of colors with a minimum length of 2.")
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
    
    # Iterate over provided geneNames
    for(metadataName in metadataColumns) {
      # Find data
      tsneLoc = groupedSCS[[groupName]]@tsne
      tsneLoc$Cell = rownames(tsneLoc)
      
      # Figure out whether there's enough metadata for this group and metadata parameter to continue
      if(length(unique(scsMetadata[which(scsMetadata$Cell %in% tsneLoc$Cell), metadataName])) < 2) {
        print(paste0("Skipping parameter '",metadataName,"' for group '",groupName,"'; not enough data."))
        next
      }
      
      # Merge metadata into tSNEData
      mergeData = merge(tsneLoc, scsMetadata[,c('Cell',metadataName)], by='Cell', all.x=T, all.y=F)
      colnames(mergeData) <- c('Cell','V1','V2','metadataParameter')
      
      # Process expression data according to expressionScale
      expressionScaleTag = ''
      if(expressionScale != 'linear') {
        mergeData$metadataParameter = get(expressionScale)(mergeData$metadataParameter)
        expressionScaleTag = expressionScale
      }

      # Plot graph
      expressionTSNE = ggplot(
      ) + geom_point( # Plot 'background' of cells without a value for this metadata parameter first
        data=mergeData[which(is.na(mergeData$metadataParameter)),], 
        aes(
          x=V1, y=V2
        ),
        size=1.5,
        color='#CCCCCC'
      ) + geom_point( # Then plot 'background' of cells without a value for this metadata parameter first
        data=mergeData[which(!is.na(mergeData$metadataParameter)),], 
        aes(
          x=V1, y=V2,
          color=mergeData[which(!is.na(mergeData$metadataParameter)),]$metadataParameter
        ),
        size=1.5*pointSizeModifier
      ) + tSNEColoring + theme_classic(
      ) + theme(
        aspect.ratio = config$plotOptions$aspect.ratio,
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()
      ) + labs(
        x=element_blank(),
        y=element_blank(),
        title=paste0(groupName," - t-SNE - ",strsplit(metadataName,'__')[[1]][1]),
        color=paste0(expressionScaleTag,' ',metadataName)
      ) + overrideTheme
      if(outputMode == 'show') {
        print(expressionTSNE)
      } else {
        fileName = paste0(groupName,'_',strsplit(metadataName,'__')[[1]][1],'_',expressionScale,'.',outputMode)
        ggsave(fileName, path=saveDir, useDingbats = FALSE)
        print(paste0("t-SNE map saved as ",saveDir,.Platform$file.sep,fileName))
      }
    }
  }
}
