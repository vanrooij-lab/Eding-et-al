### plotGeneExpressionTSNE(config, groupedSCS, genenames=NULL, groupNames=NULL, outputMode='save')
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
plotGeneExpressionTSNE <- function(config, groupedSCS, geneNames=NULL, groupNames=NULL, includeOutlierCells=T, colorScheme='Spectral', useSquish=0, expressionScale='linear', pointSizeModifier = 1, overrideDir=F, overrideTheme = NULL, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }

  #Sanitize input includeOutlierCells
  if(!is.logical(includeOutlierCells)) {stop('Parameter includeOutlierCells must be either TRUE or FALSE. Pay attention to NOT using quote')}

  #Sanitize input geneNames
  if(is.null(geneNames)) {stop("You did not provide (a) gene name(s) to draw (a) t-SNE map(s) for")}
  for(geneName in geneNames) {
    for(groupName in groupNames) {
      #Find geneName in correct format
      newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@ndata), value=T)
      
      if(length(newGeneName) == 0) {
        #Gene not found
        stop(paste0("You want to analyze marker gene ",geneName," which is not expressed in group ",groupName))
      } else if(length(newGeneName) > 1) {
        #Gene name matches multiple
        stop(paste0("You want to generate a t-SNE for gene ",geneName," which is ambiguous. Did you mean: ", paste(newGeneName, collapse = " OR "),"?"))
      } else {
        #Gene is unique, update identifier
        geneNames[match(geneName, geneNames)] = newGeneName
      }
    }
  }
  
  #Sanitize input expressionScale
  if(!(expressionScale %in% c('linear','log','log2','log10'))) {stop(paste0("Value '",expressionScale,"' is not an accepted value for expressionScale. Use 'linear' (default), 'log', 'log2' or 'log10'"))}
    
  #Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  #Iterate over provided groups
  for(groupName in groupNames) {
    #Create directory if necessary
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
    
    #Iterate over provided geneNames
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
      
      #Process expression data according to expressionScale
      expressionScaleTag = ''
      if(expressionScale != 'linear') {
        cellsForTSNE$expression = get(expressionScale)(cellsForTSNE$expression)
        expressionScaleTag = expressionScale
      }
      
      #Build color scheme
      if(is.null(colorScheme)) {
        tSNEColoring = NULL
      } else if(colorScheme %in% c('magma','inferno','plasma','viridis','cividis')) {
        tSNEColoring = scale_color_gradientn(
          colors = viridis_pal(option=colorScheme)(100),
          limits=c(
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*useSquish))],
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*(1-useSquish)))] 
          ),
          oob=squish
        )
      } else if(length(colorScheme) == 1) {
        # TODO Implement a color check for existing palettes here.
        tSNEColoring = scale_color_distiller(
          palette=colorScheme,
          limits=c(
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*useSquish))],
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*(1-useSquish)))] 
          ),
          oob=squish
        )
      } else if(length(colorScheme) == 2) {
        tSNEColoring = scale_color_gradient(
          low=colorScheme[1],
          high=colorScheme[2],
          guide=guide_colorbar(barheight=15),
          limits=c(
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*useSquish))],
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*(1-useSquish)))] 
          ),
          oob=squish
        )
      } else if(length(colorScheme) == 3) {
        tSNEColoring = scale_color_gradient2(
          low=colorScheme[1],
          mid=colorScheme[2],
          high=colorScheme[3],
          guide=guide_colorbar(barheight=15),
          limits=c(
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*useSquish))],
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*(1-useSquish)))] 
          ),
          oob=squish
        )
      } else if(length(colorScheme) > 3) {
        tSNEColoring = scale_color_gradientn( # change by mw to allow custom gradient
          
          colours=colorScheme,
          guide=guide_colorbar(barheight=15),
          limits=c(
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*useSquish))],
              cellsForTSNE$expression[order(cellsForTSNE$expression, decreasing=F)][max(1,ceiling(nrow(cellsForTSNE)*(1-useSquish)))] 
          ),
          oob=squish
        )
      } else {
        stop("Did not understand your input to 'colorScheme', needs to be a list of colors with a minimum length of 2.")
      }
  

      #Plot graph
      expressionTSNE = ggplot(
        cellsForTSNE, 
        aes(x=V1, y=V2)
      ) + geom_point(
        aes(
          color=expression
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
        title=paste0(groupName," - t-SNE - ",strsplit(geneName,'__')[[1]][1]),
        color=paste0(expressionScaleTag,"\nExp")
      ) + overrideTheme
      if(outputMode == 'show') {
        print(expressionTSNE)
      } else {
        fileName = paste0(groupName,'_',strsplit(geneName,'__')[[1]][1],'_',expressionScale,'.',outputMode)
        ggsave(fileName, path=saveDir, useDingbats=FALSE)
        print(paste0("t-SNE map saved as ",saveDir,.Platform$file.sep,fileName))
      }
    }
  }
}