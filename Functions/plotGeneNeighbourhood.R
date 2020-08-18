### Plot gene neighbourhood
### 
### Parameters
###   - 
### 
### Return
###   - 
### 
### Authors
### - Martijn Wehrens - Original idea and script
### - Joep Eding - Adaptation for accessibility and uniformity with other functions
### 
### TODO
### - Add a title to the graph
### - Add some of the gene selection parameters (and gene count) in the subtitle
plotGeneNeighbourhood <- function(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, groupName=NULL, geneOrRegulon=NULL, nameOrNum=NULL, extraStepsOut=0, plotMode='geodist', colorBy='Regulon', edgeCategory='direction', whichCorrelations="Both", plotLabels=T, overrideDir=F, outputMode='pdf') {
    # Sanitize input outputMode & overrideDir
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  } else {
    # Create directory if graphs need to be saved
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,'regulons')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir)
    }
    if(!dir.exists(saveDir)) {
      dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print(paste0("Created directory: ",saveDir))
    }
  }
  
  # Sanitize input parameter whichCorrelations
  if(!(whichCorrelations %in% c('Pos','Neg','Both'))) {
    stop("Invalid value for parameter whichCorrelations: '",whichCorrelations,"'. Should be one of 'Pos','Neg' or 'Both'.")
  }

  # Sanitize input correlationMatrices & regulons & groupedSCS & clusterDifferentialExpression & groupName
  # - Check that groupName is in all
  if(!(groupName %in% names(correlationMatrices))) {
    stop("Specified groupName '",groupName,"' is not in correlationMatrices. Options are: '",paste(names(correlationMatrices), collapse = "','"),"'.")
  }
  if(!(groupName %in% names(regulons))) {
    stop("Specified groupName '",groupName,"' is not in regulons. Options are: '",paste(names(regulons), collapse = "','"),"'.")
  }
  if(!(groupName %in% names(groupedSCS))) {
    stop("Specified groupName '",groupName,"' is not in groupedSCS. Options are: '",paste(names(groupedSCS), collapse = "','"),"'.")
  }
  if(!(groupName %in% names(clusterDifferentialExpression))) {
    stop("Specified groupName '",groupName,"' is not in clusterDifferentialExpression. Options are: '",paste(names(clusterDifferentialExpression), collapse = "','"),"'.")
  }
  
  # Sanitize geneOrRegulon and nameOrNum
  if(is.null(geneOrRegulon)) {
    stop("No value provided for parameter 'geneOrRegulon'. Needs to be either 'Gene' or 'Regulon'.")
  }
  if(is.null(nameOrNum)) {
    stop("No value provided for parameter 'nameOrNum'. Needs to be either a gene name or regulon number, depending on the value of parameter 'geneOrRegulon'.")
  }
  if(geneOrRegulon == 'Gene') {
    # Find original gene name
    newGeneName = grep(nameOrNum, rownames(groupedSCS[[1]]@ndata), value = T)
    if(length(newGeneName) == 0) {
      stop("Gene '",nameOrNum,"' is not expressed")
    }
    if(length(newGeneName) > 1) {
      stop("Gene '",nameOrNum,"' is ambiguous. Did you mean: '",paste(newGeneName,collapse = "','"),"'?")
    }
    nameOrNum = newGeneName
  } else if(geneOrRegulon == 'Regulon') {
    if(!is.numeric(nameOrNum)) {
      stop("If parameter geneOrRegulon is 'Regulon', then nameOrNum needs to be a number. Options are: ",paste(seq(1:length(names(regulons[[groupName]]$regulons))),collapse = ', '),".")
    }
    if(!(nameOrNum %in% seq(1, length(names(regulons[[groupName]]$regulons))))) {
      stop("Invalid value for parameter nameOrNum. Needs to be the number of a regulon. Options are: ",paste(seq(1:length(names(regulons[[groupName]]$regulons))),collapse = ', '),".")
    }
    nameOrNum = paste0('regulon.',nameOrNum)
  }
  
  # Sanitize input plotMode
  if(plotMode %in% c('seham', 'adj', 'mds', 'princoord', 'rmds', 'hall', 'geodist', 'segeo', 'eigen')) {
    # Do nothing. These values will result in the same plotting with every re-run.
  } else if(plotMode %in% c('circle', 'circrand', 'fruchtermanreingold', 'kamadakawai', 'random', 'spring', 'springrepulse', 'target')) {
    warning(paste0("Network plot mode '",plotMode,"' will not result in plotting on the same coordinates with a second run."))
  } else {
    stop("Invalid value for parameter 'plotMode'. Valid values are: 'seham', 'adj', 'mds', 'princoord', 'rmds', 'hall', 'geodist', 'segeo', 'eigen', 'circle', 'circrand', 'fruchtermanreingold', 'kamadakawai', 'random', 'spring', 'springrepulse' and 'target'")
  }
  
  # Sanitize input colorBy
  
  
  # Sanitize input plotLabels
  if(is.character(plotLabels)) {
    plotLabels = grep(paste0('(^',paste(plotLabels, collapse='_)|(^'),'_)'), rownames(groupedSCS[[1]]@ndata), value=T)
  } else if(is.logical(plotLabels)) {
    # Do nothing, T/F can be passed to this function
  } else {
    stop("Incorrect value for parameter 'plotLabels'. Can be either TRUE, FALSE of a list of gene names.")
  }
  
  # Make a matrix for plotting
  significanceMatrix <- correlationMatrices[[groupName]]$filteredPValueMatrix < correlationMatrices$patientAllMod$config$desiredPValue
  # If direction of the correlation is important for the plotting, then filter here
  if(whichCorrelations == 'Pos') {
    directionalityMatrix <- correlationMatrices[[groupName]]$filteredCorrelationMatrix > 0
  } else if(whichCorrelations == 'Neg') {
    directionalityMatrix <- correlationMatrices[[groupName]]$filteredCorrelationMatrix < 0
  } else if(whichCorrelations == 'Both') {
    directionalityMatrix <- significanceMatrix # Only used to have a matrix with T/F on the same places as the significanceMatrix
  }
  # Combine significance and directionality to arrive at a final matrix of edges
  binaryCorrelationMatrix = significanceMatrix & directionalityMatrix
  
  # Set colnames and rownames
  rownames(binaryCorrelationMatrix) <- rownames(correlationMatrices$patientAllMod$correlationMatrix)
  colnames(binaryCorrelationMatrix) <- rownames(correlationMatrices$patientAllMod$correlationMatrix)
  
  # Select which genes to plot
  cat("-Selecting genes to plot: ")
  if(geneOrRegulon == 'Gene') {
    # Only select genes that have a correlation with the selected gene
    selectedGenes = names(which(binaryCorrelationMatrix[,nameOrNum] | is.na(binaryCorrelationMatrix[,nameOrNum])))
  } else if(geneOrRegulon == 'Regulon') {
    # Get names of genes that are in the regulon
    selectedGenes = regulons[[groupName]]$regulons[[nameOrNum]]
    # Only select genes that have a correlation with any of the genes in the selected regulon
    selectedGenes = names(which(rowSums(binaryCorrelationMatrix[,selectedGenes]*1, na.rm = T) > 0) )
  }
  
  # Optionally look at genes that are correlated with with any of the selected genes
  if(extraStepsOut > 0) {
    # Select all genes that are correlated with the currently selected genes
    # Repeat this step for the amount of times that extraStepsOut specifies
    cat("\n---Started with ",length(selectedGenes)," genes. Expanding ",extraStepsOut," steps outward.\n")
    for(stepNum in seq(1, extraStepsOut)) {
      selectedGenes = names(which(rowSums(binaryCorrelationMatrix[,selectedGenes]*1, na.rm = T) > 0) )
      cat(paste0("-----Step ",stepNum,": ",length(selectedGenes)," genes selected.\n"))
    }
  } else {
    cat(paste0(length(selectedGenes)," genes selected.\n"))
  }
  
  # Select matrix to plot
  binaryCorrelationMatrix <- binaryCorrelationMatrix[selectedGenes,selectedGenes]
  
  cat("-Generate network.\n")
  # Generate gene network
  geneNetwork <- network(binaryCorrelationMatrix, directed = F)
  geneNetworkData <- ggnetwork(geneNetwork, layout=plotMode)
  
  # Add vertex coloring information to network
  cat("-Determine node coloring.\n")
  if(colorBy == 'Regulon') {
    # Add regulon column to networkData
    geneNetworkData$Regulon = rep(as.character(NA), nrow(geneNetworkData))
    # Iterate over vertices
    for(vertexName in names(regulons[[groupName]]$clustering)) {
      # Add regulon number to each vertex
      geneNetworkData[which(
          #geneNetworkData$x == geneNetworkData$xend &  # Used to identify a vertex instead of an edge
          #geneNetworkData$y == geneNetworkData$yend & # Used to identify a vertex instead of an edge
          geneNetworkData$vertex.names == vertexName
        ), 'Regulon'] <- as.integer(regulons[[groupName]]$clustering[vertexName])
    }
  } else if(colorBy == 'Cluster') {
    # Add cluster column to networkData
    geneNetworkData$Cluster = rep(as.character(NA), nrow(geneNetworkData))
    # Iterate over vertices
    for(vertexName in regulons[[groupName]]$config$selectedGenes) {
      highestCluster = 1
      highestExpr = 0
      # Find cluster in which this gene is expressed highest
      for(clusterNum in unique(groupedSCS[[groupName]]@cluster$kpart)) {
        if(clusterDifferentialExpression[[groupName]][[paste0('cl.',clusterNum)]][vertexName,'mean.cl'] > highestExpr) {
          highestCluster = clusterNum
          highestExpr = clusterDifferentialExpression[[groupName]][[paste0('cl.',clusterNum)]][vertexName,'mean.cl']
        }
      }
      # Add cluster number to each vertex
      geneNetworkData[which(
        #geneNetworkData$x == geneNetworkData$xend &  # Used to identify a vertex instead of an edge
        #geneNetworkData$y == geneNetworkData$yend &    # Used to identify a vertex instead of an edge
        geneNetworkData$vertex.names == vertexName
      ), 'Cluster'] = highestCluster
    }
  } else if(colorBy == 'TF') {
    # Load transcription factor list first
    if(config$species == 'human') {
      source(paste0(config$scriptsDirectory,'Resources/humanTFs.R'))
      listOfTranscriptionFactors = humanTranscriptionFactors
    } else if(config$species == 'mouse') {
      stop("Transcription factor list not available for organism '",config$species,"'.")
    } else {
      stop("Transcription factor list not available for organism '",config$species,"'.")
    }
    # This part currently only works with config$geneIdentifierType == external_gene_name_with_chr
    # TODO: Fix this
    if(config$geneIdentifierType != 'external_gene_name_with_chr') {
      stop("This function can currently only do colorBy='TF' when config$geneIdentifierType == 'external_gene_name_with_chr'. Find a bio-informatician to fix this, or use a different vertex coloring option.")
    }
    # Add TF column to networkData
    geneNetworkData$TF = 'No'
    # Iterate over vertices
    for(vertexName in unique(geneNetworkData$vertex.names)) {
      # Check whether this vertex is in the transcription factor list
      if(gsub("__chr(\\d+|[MXY])", '', vertexName) %in% listOfTranscriptionFactors) {
        # Add TF status to vertex
        geneNetworkData[which(
          #geneNetworkData$x == geneNetworkData$xend &  # Used to identify a vertex instead of an edge
          #geneNetworkData$y == geneNetworkData$yend &   # Used to identify a vertex instead of an edge
          geneNetworkData$vertex.names == vertexName
        ), 'TF']= 'Yes'
      }
    }
  }
  
  # Add edge coloring information to network
  cat("-Determine edge category\n")
  if(edgeCategory == 'direction') {
    #set.edge.attribute(geneNetwork, edgeColoring, value=NA) # This probably isn't doing much..
    # Iterate over binaryCorrelationMatrix
    for(rowNum in 1:nrow(binaryCorrelationMatrix)) {
      for(colNum in rowNum:nrow(binaryCorrelationMatrix)) {
        # Only proceed if this is a significant correlation (otherwise there isn't even an edge to color)
        if(!is.na(binaryCorrelationMatrix[rowNum, colNum] == FALSE)) {
          if(binaryCorrelationMatrix[rowNum, colNum] == FALSE) {next}
        }
        
        # Find coordinates for the two genes
        rowGene = rownames(binaryCorrelationMatrix)[rowNum]
        colGene = colnames(binaryCorrelationMatrix)[colNum]
        rowGeneCoord = geneNetworkData[which(
          geneNetworkData$vertex.names == rowGene &
          geneNetworkData$x == geneNetworkData$xend & 
          geneNetworkData$y == geneNetworkData$yend
        ), c('x','y')]
        if(nrow(rowGeneCoord) > 1) {
          rowGeneCoordX = as.double(rowGeneCoord[1,'x'])
          rowGeneCoordY = as.double(rowGeneCoord[1,'y'])
        } else {
          rowGeneCoordX = as.double(rowGeneCoord['x'])
          rowGeneCoordY = as.double(rowGeneCoord['y'])
        }
        colGeneCoord = geneNetworkData[which(
          geneNetworkData$vertex.names == colGene &
          geneNetworkData$x == geneNetworkData$xend & 
          geneNetworkData$y == geneNetworkData$yend
        ), c('x','y')]
        if(nrow(colGeneCoord) > 1) {
          colGeneCoordX = as.double(colGeneCoord[1,'x'])
          colGeneCoordY = as.double(colGeneCoord[1,'y'])
        } else {
          colGeneCoordX = as.double(colGeneCoord['x'])
          colGeneCoordY = as.double(colGeneCoord['y'])
        }

        #rowGeneCoord = unlist(geneNetworkData[rowGene, c('x','y')])
        #colGeneCoord = unlist(geneNetworkData[colGene, c('x','y')])
        
        # Find corresponding edge (edge can go from either rowGene -> colGene or the other way around, but not both)
        edgeRow = which(
          (
            geneNetworkData$x == rowGeneCoordX & 
            geneNetworkData$y == rowGeneCoordY & 
            geneNetworkData$xend == colGeneCoordX & 
            geneNetworkData$yend == colGeneCoordY
          )
          |
          (
            geneNetworkData$x == colGeneCoordX & 
            geneNetworkData$y == colGeneCoordY & 
            geneNetworkData$xend == rowGeneCoordX & 
            geneNetworkData$yend == rowGeneCoordY
          )
        )
        
        # Find correlation direction 
        if(correlationMatrices[[groupName]]$correlationMatrix[rowGene, colGene] >= 0) {
          correlationDirection = 'Pos'
        } else {
          correlationDirection = 'Neg'
        }
        
        # Change value for edge
        geneNetworkData[edgeRow, edgeCategory] = correlationDirection
      }
    }
  }

  # Start network plot, plot edges and nodes
  cat("-Plot network.\n")
  geneNetworkPlot <- ggplot(
    geneNetworkData, 
    aes(x = x, y = y, xend = xend, yend = yend)
  ) + geom_edges(
    aes_string(linetype=edgeCategory),
    size=0.1,
    color='grey80'
  ) + geom_nodes( # Plot the uncolored nodes first
    aes_string(color=colorBy),
    size=1,
    data=function(x) { x[is.na(x[,colorBy]), ]}
  ) + geom_nodes( # Then the colored nodes go on top to make them stand out more
    aes_string(color=colorBy),
    size=1,
    data=function(x) { x[!is.na(x[,colorBy]), ]}
  ) + theme_blank() 
  # Optionally add text labels
  if(is.character(plotLabels) | (is.logical(plotLabels) && plotLabels)) {
    if(is.logical(plotLabels)) {
      # Plot all vertex names
      geneNetworkPlot = geneNetworkPlot + geom_nodetext(
        aes_string(
          color = colorBy, 
          label = 'vertex.names'
        ),
        size = 2
      )  
    } else {
      # Plot only the specified labels
      geneNetworkPlot = geneNetworkPlot + geom_nodetext(
        aes_string(
          color = colorBy, 
          label = 'vertex.names'
        ),
        size = 2,
        data = function(x) { x[ x$vertex.names %in% plotLabels, ]}
      )  
    }
  }
  #print(geneNetworkPlot)

  #geneNetworkPlot = ggnet2(geneNetwork, mode=plotMode, node.size = 1, color = colorBy, palette = colorsForNetwork, label=plotLabels, label.size = 2)
  if(outputMode != 'show') {
    plotLabelTag = ifelse(is.logical(plotLabels), plotLabels, substr(paste(plotLabels, collapse = '-'),0,10))
    fileName = paste0("geneNeighbourHood",groupName,"_selected",nameOrNum,"and",extraStepsOut,"extraStepsOut_corDir",whichCorrelations,"_plotMode",plotMode,"_colorBy",colorBy,'_edgeCategory',edgeCategory,'_plotLabels',plotLabelTag,'.',outputMode)
    ggsave(fileName,path = saveDir, useDingbats=FALSE)
    print(paste0("Saved plot as: ",saveDir,'/',fileName))
  } else {
    print(geneNetworkPlot)
  }
}