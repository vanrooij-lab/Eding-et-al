### Plot regulon networks
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
### - Spectral colorscheme can't handle more than 10/11 colors. So function crashes with more regulons than that
### - Add option to limit coloring of edges to those connected to a specified (set of) gene(s) (vertices)
### - If colorBy=='TF', this currently only works when config$geneIdentifierType == 'external_gene_name_with_chr'. Fix this
### - Add which edges (pos, neg, both) are being plotted in the output file name as well
plotRegulons <- function(config, groupedSCS, clusterDifferentialExpression, correlationMatrices, regulons, plotMode='geodist', colorBy='Regulon', edgeCategory='Direction', whichCorrelations="Both", plotLabels=T, plotFullMatrix=F, overrideDir=F, outputMode='pdf') {
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
  
  # Sanitize input correlationMatrices & regulons & groupedSCS & clusterDifferentialExpression
  # - Check that groupnames in both are the same
  
  # Sanitize input plotFullMatrix
  # - Check that it's LOGICAL
  
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
  
  # Iterate over groups
  for(groupName in names(regulons)) {
    cat(paste0("Plotting regulons for '",groupName,"':\n"))
    
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
    
    # Optionally select only the genes in regulons
    if(!plotFullMatrix) {
      # Only select selected genes from the full matrix
      binaryCorrelationMatrix <- binaryCorrelationMatrix[regulons[[groupName]]$config$selectedGenes,regulons[[groupName]]$config$selectedGenes]
    }
    
    # Remove 'singletons' from binarycorrelationMatrix
    cat(paste0("-Removing 'single' genes':\n"))
    removedGenes = c('start_loop') # Make sure loop gets started
    while(length(removedGenes) > 0) {
      oldRowNames <- binaryCorrelationMatrix %>% rownames
      binaryCorrelationMatrix <- binaryCorrelationMatrix[apply(binaryCorrelationMatrix, 1, function(x) {length(which(x))}) > 1, apply(binaryCorrelationMatrix, 1, function(x) {length(which(x))}) > 1]
      newRowNames <- binaryCorrelationMatrix %>% rownames
      removedGenes <- oldRowNames[which(!(oldRowNames %in% newRowNames))]
      cat(paste0("---Removed ",length(removedGenes)," 'single' genes': '",paste(removedGenes, collapse = "', '"),"'\n"))
    }
    
    cat("-Generate network.\n")
    # Generate gene network
    geneNetwork <- network(binaryCorrelationMatrix, directed = F)
    geneNetworkData <- ggnetwork(geneNetwork, layout=plotMode)
    
    # Add vertex coloring information to network
    cat("-Determine node coloring.\n")
    nodeColoring = NULL
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
      nodeColoring = regulons[[groupName]]$geneClusteringPalette
      names(nodeColoring) <- as.character(seq(1,length(nodeColoring)))
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
      data=function(x) { x[is.na(x[,colorBy]), ]},
    ) 
    # Optionally specify the colors for the nodes
    if(!is.null(nodeColoring)) {
      geneNetworkPlot <- geneNetworkPlot + scale_colour_manual(
        values=nodeColoring
      )
    }
    geneNetworkPlot <- geneNetworkPlot + geom_nodes( # Then the colored nodes go on top to make them stand out more
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
      plotFullMatrixTag = ifelse(plotFullMatrix, 'FullMatrix','SubMatrix')
      fileName = paste0("RegulonNetwork_",groupName,"_plot",plotFullMatrixTag,"_plotMode",plotMode,"_colorBy",colorBy,'_edgeCategory',edgeCategory,'_plotLabels',plotLabelTag,'.',outputMode)
      ggsave(fileName,path = saveDir, useDingbats=FALSE)
      print(paste0("Saved plot as: ",saveDir,'/',fileName))
    } else {
      print(geneNetworkPlot)
    }
  }
}
