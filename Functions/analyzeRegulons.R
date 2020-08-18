### Function that draws a heatmap from a correlation plot and returns a list of regulons
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
### - Add proper labels to graphs
### - Add thorough sanitization for input parameters
### - Data sanitization for combinations of cluster cutoff values (give warnings regarding priority)
### - Make clusGap parameters available to user as function parameters (and/or provide alternative to gap method)
### - Visualize gap statistics results
### - Actually respect the settings to save graphs
### - Input sanitization for overrideClusterColors
analyzeRegulons <- function(config, correlationMatrices, minCorrelations=0, clusteringMethod='ward.D2', overrideClusterNum=F, overrideRegulonColors=F, useSquish=0, overrideDir=F, outputMode='pdf', K.max=25, filename_precursor='') {
  # Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf'))) {
    stop("'outputMode' needs to be one of 'show','png' or 'pdf'. Pay attention to capitalisation.")
  } else {
    # Create directory if graphs need to be saved
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,'regulons','/')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,'/')
    }
    if(!dir.exists(saveDir)) {
      dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print(paste0("Created directory: ",saveDir))
    }
  }
  
  # Sanitize input parameter clusteringMethod
  if(!(clusteringMethod %in% c('ward.D','ward.D2','single','complete','average','mcquitty','median','centroid'))) {
    stop("Invalid entry '",clusteringMethod,"' for parameter clusteringMethod. Should be one of 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  
  # Initialize return variable
  returnData = list()
  
  # Iterate over groupnames
  for(groupName in names(correlationMatrices)) {
    # Print status
    print(paste0("Analyzing regulons for '",groupName,"'."))
    
    # Select genes for 'regulon' matrix
    selectedGenes <- correlationMatrices[[groupName]]$correlationsPerGene[which(correlationMatrices[[groupName]]$correlationsPerGene > minCorrelations)] %>% names
    
    # Select subsection of correlationMatrix for these genes
    subCorrelationMatrix = correlationMatrices[[groupName]]$correlationMatrix[selectedGenes,selectedGenes]
    K.max = min(K.max, nrow(subCorrelationMatrix)-1) # this prevents errors when there are less genes than clusters
    
    # Hierarchically cluster genes based on correlations
    geneClusteringTree <- hclust(as.dist(1-subCorrelationMatrix), method=clusteringMethod)
    # Try to determine optimal point of cutting dendogram to determine clusters
    if(overrideClusterNum == FALSE) {
      print(paste0(" - Automatically determining optimal cluster number: "))
      # Automatically determine cluster cutoff height
      gap_stat <- clusGap(subCorrelationMatrix, FUN = kmeans, nstart = 10, K.max = K.max, B = 10)
      nCluster = maxSE(f=gap_stat$Tab[,'gap'], SE.f = gap_stat$Tab[,'SE.sim'])
      #fviz_gap_stat(gap_stat)
      cat(paste0(nCluster," clusters is optimal.\n"))
    } else {
      nCluster = overrideClusterNum
    }
    
    # Use determined number of clusters to cut dendogram into clusters
    geneClustering = cutree(geneClusteringTree, k=nCluster)
    
    # Determine cluster colors
    if(length(overrideRegulonColors) > 1) {
      # An override list of cluster colors has been specified, use this instead of random colors
      # But check whether enough colors have been specified first
      if(length(overrideRegulonColors) < nCluster) {
        stop("You have specified a list of ",length(overrideRegulonColors)," colors. There are ",nCluster," regulons. Please specify at least ",nCluster," colors.")
      }
      geneClusteringPalette = overrideRegulonColors
    } else {
      geneClusteringPalette = distinctColorPalette(max(geneClustering))
    }
    geneClusteringColors <- geneClusteringPalette[geneClustering]      
    
    # Calculate silhouette score to get a measure of how well each gene fits in its regulon
    clusteringSilhouette <- silhouette(geneClustering, as.dist(1-subCorrelationMatrix))
    silhouetteData <- as.data.frame(clusteringSilhouette[,])
    rownames(silhouetteData) <- names(geneClustering)
    
    # Generate heatmap with clusters
    #   - First, determine 98% interval for color boundaries (if required)
    subFilteredCorrelationMatrix <- subCorrelationMatrix
    subFilteredCorrelationMatrix[row(subFilteredCorrelationMatrix)==col(subFilteredCorrelationMatrix)] <- NA
    subFilteredCorrelationVector <- as.vector(subFilteredCorrelationMatrix)
    lowerBoundary <- subFilteredCorrelationVector[order(subFilteredCorrelationVector,decreasing = F)][max(1,round(useSquish*length(which((subFilteredCorrelationVector<1)))))]
    upperBoundary <- subFilteredCorrelationVector[order(subFilteredCorrelationVector,decreasing = F)][max(1,round((1-useSquish)*length(which((subFilteredCorrelationVector<1)))))]
    stepSize   <- (upperBoundary-lowerBoundary)/100
    
    # Open appropriate device to save plot (if required)
    filename = paste0(filename_precursor,'regulonHeatmap_',groupName,'_minCell',correlationMatrices[[groupName]]$config$minCellFraction,'_minExp',correlationMatrices[[groupName]]$config$minCellExpression,'_rCutoff',correlationMatrices[[groupName]]$rCutoff,'_minCorr',minCorrelations,'_squish',useSquish)
    if(outputMode == 'png') {
      png(
        filename = paste0(saveDir,filename,'.png'),
        width = 1000,
        height = 1000,
        units = 'px'
      )
    } else if(outputMode == 'pdf') {
      pdf(
        file = paste0(saveDir,filename,'.pdf'),
        width = 10,
        height = 10
      )
    }
    
    #   - Then draw heatmap
    heatmap3(
        #subCorrelationMatrix, 
        squish(subCorrelationMatrix, range = c(lowerBoundary, upperBoundary)),
        ColSideColors = geneClusteringColors, 
        Colv=as.dendrogram(geneClusteringTree), 
        Rowv=as.dendrogram(geneClusteringTree), 
        col = viridis_pal(option = "D")(100),
        labRow = NA,
        labCol = NA,
        ColSideLabs = NA,
        scale='none'
        #breaks = seq(lowerBoundary,upperBoundary,stepSize),
        #tracecol=NA
    )    
    legend(
      'right',
      title = "Clusters",
      legend=seq(1, length(geneClusteringPalette)), 
      fill=geneClusteringPalette, 
      cex=0.8, 
      box.lty=0,
      inset=-0.03,
      xpd=T
    )
    
    # Close device if plot is being saved
    if(outputMode %in% c('png','pdf')) {
      dev.off()
      print(paste0("Saved sample distance plot as: ",saveDir,filename,outputMode))
    }
    
    # Save regulons to excel sheet
    # First create a workbook
    regulonBook = createWorkbook()
    # Iterate over regulons
    for(regulonNum in 1:max(geneClustering)) {
      # Add a worksheet for this regulon
      addWorksheet(regulonBook, paste0('Regulon.',regulonNum))
      # Get correlation matrix of the genes in this regulon
      regulonGeneNames = names(geneClustering[geneClustering==regulonNum])
      regulonGenes = as.data.frame(subCorrelationMatrix[regulonGeneNames,regulonGeneNames])
      # Add column with silhouette score
      regulonGenes$silhouette <- rep(NA, nrow(regulonGenes))
      for(geneName in regulonGeneNames) {
        regulonGenes[geneName,'silhouette'] = silhouetteData[geneName,'sil_width']
      }
      # Order data based on the silhouette column
      regulonGenes <- regulonGenes[order(regulonGenes$silhouette, decreasing = T),]
      # Remove chromosome names from row and colnames if present
      rownames(regulonGenes) <- gsub("__chr(\\d+|[MXY])", '', rownames(regulonGenes))
      colnames(regulonGenes) <- gsub("__chr(\\d+|[MXY])", '', colnames(regulonGenes))
      # Write data to the worksheet
      writeDataTable(regulonBook, paste0('Regulon.',regulonNum), regulonGenes, rowNames = T)
    }
    addWorksheet(regulonBook, paste0('Background'))
    writeData(regulonBook, 'Background',  gsub("__chr(\\d+|[MXY])", '', rownames(correlationMatrices[[groupName]]$correlationMatrix)))
    saveWorkbook(regulonBook, paste0(saveDir,filename,'.xlsx'), overwrite = T)
    
    # Store cluster composition in return variable
    returnData[[groupName]]$config$minCorrelations = minCorrelations
    returnData[[groupName]]$config$nCluster = nCluster
    returnData[[groupName]]$config$selectedGenes = selectedGenes
    returnData[[groupName]]$config$clusteringMethod = clusteringMethod
    returnData[[groupName]]$gap_stat = ifelse(overrideClusterNum, F, gap_stat)
    returnData[[groupName]]$clustering = geneClustering
    returnData[[groupName]]$geneClusteringPalette = geneClusteringPalette
    
    # Store two other important params
    returnData[[groupName]]$subCorrelationMatrix = subCorrelationMatrix
    returnData[[groupName]]$silhouetteData = silhouetteData
    
    for(clusterNum in 1:max(geneClustering)) {
      returnData[[groupName]]$regulons[[paste0('regulon.',clusterNum)]] = names(geneClustering[geneClustering==clusterNum])
    }
  }
  return(returnData)
}
