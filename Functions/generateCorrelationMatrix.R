### Function that analyses all correlations
### 
### Parameters
###   -
### 
### Return
###   -
### 
### Authors
###  - Martijn Wehrens - Original idea and script
###  - Joep Eding - Adaptation for accessibility and uniformity with other functions
### 
### Warning
###   -RRREEEAAALLLYYY slow function with large datasets (when scrambling expressionMatrix)
### 
### TODO
###   - Make randomization of matrix optional?
###   - Convert histograms to ggplot2
###   - Add saving options to the histograms
###   - Add option to exclude cells in outlier clusters
###   - Simplify histogram generation, presumably:
###     - The combination of hist() and geom_line can be replaced by geom_density (Seems much much slower like that)
###     - Color generation can be automatic
###     - Either 1 or 2 datasets (depending on value of scrambeExpressionMatrix)
###       can be combined in 1 dataframe with a grouping variable
###   - Uncouple rCutoff-generation from correlationMatrix generation
###     - The calculation for corrected pCutoff (and rCutoff) based on desiredPValue
###       seems to be based on just the number of cells.
###   - Make sure that pValue-generation takes into account the remaining number of genes/columns after exclusions!
generateCorrelationMatrix <- function(config, groupedSCS, groupNames, excludeOutlierCells=T, minCellFraction=0, minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf', which_genes_to_select=T, filename_precursor='') {
  # Sanitize input outputMode
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
  
  # Sanitize input parameter minCellFraction.
  if(!is.numeric(minCellFraction)) {
    stop("Invalid input '",minCellFraction,"' for parameter minCellFraction. Must be a number between 0 and 1.")
  } else if((minCellFraction < 0) | (minCellFraction > 1)) {
    stop("Invalid input '",minCellFraction,"' for parameter minCellFraction. Must be a number between 0 and 1.")
  }
  
  # Sanitize input parameter adjustP
  if(is.logical(adjustP)) {
    # TRUE and FALSE are valid values for adjustP but need to be changed to something that p.adjust understands
    if(adjustP == TRUE) {
      # In case of TRUE, change value to default method 'fdr' for p.adjust
      adjustP = 'fdr'
    } else {
      # In case of FALSE, change value to 'none', to make p.adjust a pass-through
      adjustP = 'none'
    }
  } else if(adjustP %in% p.adjust.methods) {
    # Do nothing, these are valid values.
  } else {
    stop("Invalid value '",adjustP,"' for adjustP, use one of TRUE, FALSE, '",paste(p.adjust.methods, collapse="', '"),"'")
  }
  
  # Initialize return variable
  returnData = list()
  
  # Iterate over groups
  for(groupName in groupNames) {
    # Print progress
    print(paste0("Preparing expression matrices for '",groupName,"'"))
    
    # Prepare expression matrix (optionally select genes)
    expressionMatrix = groupedSCS[[groupName]]@ndata[which_genes_to_select,]
    
    # Exclude cells in outlier clusters if required
    if(excludeOutlierCells) {
      cellsToSelect <- names(groupedSCS[[groupName]]@cpart[groupedSCS[[groupName]]@cpart <= max(groupedSCS[[groupName]]@cluster$kpart)])
      expressionMatrix <- expressionMatrix
    }
    
    # Select genes based on cutoffs, per gene:
    #   -count number of cells expressing gene above threshold minCellExpression
    #   -keep genes with number of cells above threshold minCellFraction
    rowSelector <- rowSums((expressionMatrix > minCellExpression)*1)/ncol(expressionMatrix) >= minCellFraction
    expressionMatrix <- expressionMatrix[rowSelector,]
    expressionMatrixDimensions <- dim(expressionMatrix)
    
    # Generate correlation matrix for expressionMatrix
    print(paste0("Generating correlation matrix for actual data"))
    startTime <- Sys.time()
    correlationMatrix <- cor(t(expressionMatrix))
    endTime <- Sys.time(); 
    print(paste0('This took from ',startTime,' until ',endTime))
    endTime - startTime   
    
    # Generate a pValue matrix
    # Generate correlation matrix for expressionMatrix
    print(paste0("Generating pValue matrix for actual data"))
    startTime <- Sys.time()
    degreesOfFreedom <- ncol(expressionMatrix)-2
    pValueMatrix = matrix(
      sapply(
        as.vector(correlationMatrix),
        function(R, df) {
          t1<-sqrt(df) * R/sqrt(1 - R^2)
          pValue <- 2*min(pt(t1,df),pt(t1,df,lower.tail = FALSE))
          return(pValue)          
        },
        df=degreesOfFreedom
      ),
      nrow=expressionMatrixDimensions[1]
    )
    endTime <- Sys.time(); 
    print(paste0('This took from ',startTime,' until ',endTime))
    endTime - startTime  
    # Now adjust the pValue matrix according to the method selected by the user
    print(paste0("Adjusting pValue matrix for actual data"))
    startTime <- Sys.time()
    degreesOfFreedom <- ncol(expressionMatrix)-2
    pValueMatrix <- matrix(
      p.adjust(
        as.vector(pValueMatrix),
        method = adjustP
      ), 
      nrow=expressionMatrixDimensions[1]
    )
    endTime <- Sys.time(); 
    print(paste0('This took from ',startTime,' until ',endTime))
    endTime - startTime  
    
    # Filter 'perfect correlations' on the diagonal of the matrix (change to NA)
    # Do the same for the pValueMatrix
    filteredCorrelationMatrix = correlationMatrix
    filteredCorrelationMatrix[row(filteredCorrelationMatrix)==col(filteredCorrelationMatrix)] <- NA
    filteredPValueMatrix = pValueMatrix
    filteredPValueMatrix[row(filteredPValueMatrix)==col(filteredPValueMatrix)] <- NA
    
    # Count number of significant correlations per gene
    correlationsPerGene <- apply(
      1*(filteredPValueMatrix < desiredPValue),
      2, 
      sum, 
      na.rm=T
    )
    names(correlationsPerGene) <- colnames(correlationMatrix)
    
    # Graph the number of genes that have more than X correlations
    # This can be used to determine the 'connectedness' cutoff for the regulon identification
    genesPerCutoff = data.frame(
      'cutoff' = 1:(0.5*max(correlationsPerGene))
    )
    genesPerCutoff$nGenes = NA
    for(rowNum in 1:nrow(genesPerCutoff)) {
      genesPerCutoff[rowNum,'nGenes'] = length(which(correlationsPerGene > rowNum))
    }
    ggplot(
      data = genesPerCutoff
    ) + geom_line(
      aes(x=cutoff,y=nGenes)
    )
    
    # Add calculated data to return variable
    returnData[[groupName]] = list(
      'correlationMatrix' = correlationMatrix,
      'filteredCorrelationMatrix' = filteredCorrelationMatrix,
      'pValueMatrix' = pValueMatrix,
      'filteredPValueMatrix' = filteredPValueMatrix,
      'correlationsPerGene' = correlationsPerGene,
      'genesPerCutoff' = genesPerCutoff,
      'config' = list(
        'minCellFraction' = minCellFraction,
        'minCellExpression' = minCellExpression,
        'desiredPValue' = desiredPValue,
        'adjustP' = adjustP
      )
    )
    
    # Save the generated correlationMatrix
    if(saveMatrices) {
      write.csv(correlationMatrix, paste0(saveDir,'/',filename_precursor,'_',groupName,'_correlationMatrix.csv'))  
    }
  }
  return(returnData)
}