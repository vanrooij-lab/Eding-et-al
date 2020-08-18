### analyseMetaDataCorrelation(config, groupedSCS, geneName, groupName, cellNames, removeNoExpressCells, percentage, correlationMethod)
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### geneName: Pass a single gene name (the 'gene of interest') to calculate the correlation with all other genes for
### groupName: Pass a single groupName.
### cellNames: Pass a list of names of individual cells or leave as NULL to use all cells in this group
### removeNoExpressCells: Three options:
###                         -no: no cells will be excluded, all will be used for correlation analysis
###                         -goi: cells not expressing gene-of-interest (goi) will be excluded, rest will be used for correlation analysis
###                         -yes: cells not expressing goi or the gene it's being compared with are excluded, rest used for corr analysis
### percentage: If cells are removed due to removeNoExpressCells settings, at least this percentage of cells must remain for the gene to 
###             be taken along for correlation analysis. If less cells remain, then this gene is skipped.
### correlationMethod: the correlationMethod to use in cor.test
### 
### Return: List containig groupName, geneName, cellNames, removeNoExpressCells, percentage, correlationMethod and the correlationData. 
# 
# Authors
#   - Joep Eding
# 
### TODO:
###   - Update comment
analyseMetaDataCorrelation <- function(config, groupedSCS, parameterName=NULL, groupName=NULL, cellNames=NULL, removeNoExpressCells='parameter', removeOutlierCells=F, percentage=0, correlationMethod='pearson') {
  #Sanitize input groupName
  if(is.null(groupName)) {stop("You did not provide the name of the group in which you want to analyze expression correlation")} 
  if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to analyze expression correlation in ",groupName,"' that is not in groupedSCS"))}

  # Check whether a metadata object exists
  if(!exists('scsMetadata')) {
    stop("Load metadata first by running loadMetaData().")
  }
  
  #Sanitize input parameterName
  if(is.null(parameterName)) {stop("You did not provide a metadata parameter name to analyse correlation for.")}
  if(!(parameterName %in% colnames(scsMetadata))) {
    print(paste0("Metadata variable '",parameterName,"' not loaded, skipping."))
    next()
  }
  
  #Sanitize input cellNames
  if(is.null(cellNames)) {
    #Use all cells
    cells = colnames(groupedSCS[[groupName]]@ndata)
  } else if(length(cellNames) < 3) {
    stop("You need to define at least 3 cell names in order to be able to analyze correlation")
  } else {
    for(cellName in cellNames) {
      if(!(cellName %in% colnames(groupedSCS[[groupName]]@ndata))) {stop(paste0("Cell '",cellName,"' does not exist in group '",groupName,"'"))}
    }
  }
    
  # Sanitize input removeNoExpressCells
  if(!(removeNoExpressCells %in% c('parameter', 'yes'))) {stop("Only valid values for 'removeNoExpressCells' are 'parameter' & 'yes'")}
  
  # Sanitize input percentage
  if(!( (is.double(percentage)) || (is.integer(percentage)) )) {stop("Input 'percentage' should be a number")}
  if((percentage < 0) || percentage > 100) {stop("Input 'percentage' should be between 0 and 100")}
  
  #Prepare some variables
  inputData <- as.data.frame(t(groupedSCS[[groupName]]@ndata))
  inputData$Cell = rownames(inputData)
  inputData <- merge(inputData, scsMetadata[,c('Cell',parameterName)], by='Cell', all.x=T, all.y=F)
  inputData <- inputData[,which(colnames(inputData) != 'Cell')]
  genes <- NULL
  allCorrelations<- NULL
  allPValues <- NULL
  averageExpression <- NULL
  nCellsExpressGenes <- NULL
  
  # Remove outlier cells if required
  if(removeOutlierCells) {
    inputData <- inputData[!(rownames(inputData) %in% names(groupedSCS[[groupName]]@cpart[which(groupedSCS[[groupName]]@cpart>max(groupedSCS[[groupName]]@cluster$kpart))])),]
  }
  
  #Process inputData according to removeNoExpressCells setting
  if(removeNoExpressCells=="no"){
    # Can't happen. It also makes no sense to take along all cells, because
    # assuming a missing parameter is zero isn't a viable option (like it is for 
    # expression levels).
  } else if(removeNoExpressCells=="parameter") {
    #Remove all cells that don't express the gene of interest
    for(i in 1:ncol(inputData)){
      #Generate matrix with expression values of geneName and comparison gene
      correlationMatrix <- data.frame(parameter=inputData[,parameterName],gene=inputData[,i])
      
      #Exclude al cells that have NA for the metadata parameter
      correlationMatrix <- correlationMatrix[!is.na(correlationMatrix$parameter),]
      
      #Only proceed with adding this gene comparison to the list if more than 3 cells express it.
      if(nrow(correlationMatrix)<3) {next}
      
      #Only proceed with adding this gene comparison to the list if more than the required percentage of cells express GOI
      if(nrow(correlationMatrix[which(correlationMatrix$gene > 0.1),]) < nrow(correlationMatrix)*(percentage/100)) {next}
      
      #Add this comparison to the list
      genes <- c(genes,colnames(inputData)[i])
      #Calculate correlation
      correlationTest <- cor.test(correlationMatrix[,1],correlationMatrix[,2],method= correlationMethod)
      #Add correlation coefficient and pValue to lists
      allCorrelations <- c(allCorrelations,as.numeric(correlationTest["estimate"]))
      allPValues <- c(allPValues,as.numeric(correlationTest["p.value"]))

      #Add mean expression and number of cells with expression to the list
      expression <- mean(correlationMatrix[,"gene"])
      averageExpression <- c(averageExpression,expression)
      nCells <- nrow(correlationMatrix[correlationMatrix[,"gene"]>0.1,])
      nCellsExpressGenes <- c(nCellsExpressGenes,nCells)
    }
  } else if(removeNoExpressCells=="yes") {
    #Remove all cells that don't express either gene of interest or comparison gene
    for(i in 1:ncol(inputData)){
      correlationMatrix <- data.frame(parameter=inputData[,parameterName],gene=inputData[,i])
      correlationMatrix <- correlationMatrix[which(
          !is.na(correlationMatrix$parameter)
          &
          correlationMatrix$gene>0.1),]
      
      #Only proceed with adding this gene comparison to the list if more than 3 cells express it.
      if(nrow(correlationMatrix)<3) {next}
      
      #Only proceed with adding this gene comparison to the list if more than the required percentage of cells express GOI
      if(nrow(correlationMatrix) < (nrow(scsMetadata[which(!is.na(scsMetadata[[parameterName]])),])*(percentage/100))) {next}
      
      #Add this comparison to the list
      genes <- c(genes,colnames(inputData)[i])

      #Calculate correlation
      correlationTest <- cor.test(correlationMatrix[,1],correlationMatrix[,2],method= correlationMethod)
      #Add correlation coefficient and pValue to lists
      allCorrelations <- c(allCorrelations,as.numeric(correlationTest["estimate"]))
      allPValues <- c(allPValues,as.numeric(correlationTest["p.value"]))
      
      #Add mean expression and number of cells with expression to the list
      expression <- mean(correlationMatrix[,"gene"])
      averageExpression <- c(averageExpression,expression)
      nCells <- nrow(correlationMatrix)
      nCellsExpressGenes <- c(nCellsExpressGenes,nCells)
    }
  }
  
  correlationDataFrame <- data.frame(
    genes=genes, 
    correlation=allCorrelations, 
    pValue = allPValues, 
    averageExpression = averageExpression,
    nCellsExpressGenes= nCellsExpressGenes
  )  
  
  returnList = list(
    'geneName' = parameterName,
    'groupName' = groupName,
    'cellNames' = cellNames,
    'removeNoExpressCells' = removeNoExpressCells,
    'percentage' = percentage,
    'correlationMethod' = correlationMethod,
    'correlationData' = correlationDataFrame
  )
  return(returnList)
}
