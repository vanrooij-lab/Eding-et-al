### analyseCorrelation(config, groupedSCS, geneName, groupName, cellNames, removeNoExpressCells, percentage, correlationMethod)
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
### 
### Authors
###   - Joep Eding
### 
### TODO
###   -
analyseCorrelation <- function(config, groupedSCS, geneName=NULL, groupName=NULL, cellNames=NULL, removeNoExpressCells='no', removeOutlierCells=F, percentage=0, correlationMethod='pearson') {
  #Sanitize input groupName
  if(is.null(groupName)) {stop("You did not provide the name of the group in which you want to analyze expression correlation")} 
  if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to analyze expression correlation in ",groupName,"' that is not in groupedSCS"))}

  #Sanitize input geneName
  if(is.null(geneName)) {stop("You did not provide a gene name to analyse correlation for.")}
  #Find geneName in correct format
  newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@ndata), value=T)
  if(length(newGeneName) == 0) {
    #Gene not found
    stop(paste0("You want to analyze correlation for gene ",geneName," which is not expressed in group ",groupName))
  } else if(length(newGeneName) > 1) {
    #Gene name matches multiple
    stop(paste0("You want to analyse correlation for gene ",geneName," which is ambiguous. Did you mean: ", paste(newGeneName, collapse = " OR "),"?"))
  } else {
    #Gene is unique, update identifier
    oldGeneName <- geneName
    geneName = newGeneName
  }
  
  #Sanitize input cellNames
  if(is.null(cellNames)) {
    #Use all cells
    cells = colnames(groupedSCS[[groupName]]@ndata)
  } else if(length(cellNames) < 3) {
    stop("You need to define at least 3 cell names in order to be able to analyze correlation")
  } else {
    stop("You're trying to analyze correlation in a specific subset of cell. This functionality is not currently working.")
    # When adding in this functionality, make sure that clashes with manually specifying cell names are adequately resolved.
    for(cellName in cellNames) {
      if(!(cellName %in% colnames(groupedSCS[[groupName]]@ndata))) {stop(paste0("Cell '",cellName,"' does not exist in group '",groupName,"'"))}
    }
  }
    
  # Sanitize input removeNoExpressCells
  if(!(removeNoExpressCells %in% c('no', 'goi', 'yes'))) {stop("Only valid values for 'removeNoExpressCells' are 'no', 'goi' & 'yes'")}
  
  # Sanitize input percentage
  if(!( (is.double(percentage)) || (is.integer(percentage)) )) {stop("Input 'percentage' should be a number")}
  if((percentage < 0) || percentage > 100) {stop("Input 'percentage' should be between 0 and 100")}
  
  # Prepare some variables
  inputData <- as.data.frame(t(groupedSCS[[groupName]]@ndata))
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
    #Don't remove any cells
    for(i in 1:ncol(inputData)){
      #Generate matrix with expression values of geneName and comparison gene
      correlationMatrix <- data.frame(gene1=inputData[,grepl(geneName,colnames(inputData))],gene2=inputData[,i])
      
      #Only proceed with adding this gene comparison to the list if more than 3 cells express it.
      if(nrow(correlationMatrix)<3) {next}
      
      #Only proceed with adding this gene comparison to the list if more than the required percentage of cells express gene2
      if(nrow(correlationMatrix[correlationMatrix$gene2>0.1,]) < (nrow(correlationMatrix)*(percentage/100))) {next}
      
      #Add this comparison to the list
      genes <- c(genes,colnames(inputData)[i])
      
      #Calculate correlation between these two genes and append to list
      correlation <- as.numeric(cor.test(correlationMatrix[,1],correlationMatrix[,2],method=correlationMethod)["estimate"])
      allCorrelations <- c(allCorrelations,correlation)
      
      #Calculate p-value for this correlation and append to list
      pValue <- as.numeric(cor.test(correlationMatrix[,1],correlationMatrix[,2],method=correlationMethod)["p.value"])
      allPValues <- c(allPValues,pValue)
      
      #Calculate mean expreesion for this correlation and append to list
      expression <- mean(correlationMatrix[,"gene2"])
      averageExpression <- c(averageExpression,expression)
      
      #Calculate number of cells expressing this gene and append to list
      nCells <- nrow(correlationMatrix[correlationMatrix[,1]>0.1&
                                                correlationMatrix[,2]>0.1,])
      nCellsExpressGenes <- c(nCellsExpressGenes,nCells)
    } 
  } else if(removeNoExpressCells=="goi") {
    #Remove all cells that don't express the gene of interest
    for(i in 1:ncol(inputData)){
      #Generate matrix with expression values of geneName and comparison gene
      correlationMatrix <- data.frame(gene1=inputData[,grepl(geneName,colnames(inputData))],gene2=inputData[,i])
      #Exclude al cells that have an expression of 0.1 for the genoe of interest
      correlationMatrixGOI <- correlationMatrix[correlationMatrix[,"gene1"]>0.1,]
      
      #Only proceed with adding this gene comparison to the list if more than 3 cells express it.
      if(nrow(correlationMatrixGOI)<3) {next}
      
      #Only proceed with adding this gene comparison to the list if more than the required percentage of cells express GOI
      if(nrow(correlationMatrixGOI) < (nrow(correlationMatrix)*(percentage/100))) {next}
      
      #Add this comparison to the list
      genes <- c(genes,colnames(inputData)[i])
      #Calculate & Add correlation
      correlation <- as.numeric(cor.test(correlationMatrixGOI[,1],correlationMatrixGOI[,2],method= correlationMethod)["estimate"])
      allCorrelations <- c(allCorrelations,correlation)
      #Calculate & add pValue
      pValue <- as.numeric(cor.test(correlationMatrixGOI[,1],correlationMatrixGOI[,2],method= correlationMethod)["p.value"])
      allPValues <- c(allPValues,pValue)
      #Add mean expression and number of cells with expression to the list
      expression <- mean(correlationMatrixGOI[,"gene2"])
      averageExpression <- c(averageExpression,expression)
      nCells <- nrow(correlationMatrix[correlationMatrix[,"gene1"]>0.1&
                                               correlationMatrix[,"gene2"]>0.1,])
      nCellsExpressGenes <- c(nCellsExpressGenes,nCells)
    }
  } else if(removeNoExpressCells=="yes") {
    #Remove all cells that don't express either gene of interest or comparison gene
    for(i in 1:ncol(inputData)){
      correlationMatrix <- data.frame(gene1=inputData[,grepl(geneName,colnames(inputData))],gene2=inputData[,i])
      correlationMatrix <- correlationMatrix[correlationMatrix[,"gene1"]>0.1&
                                                 correlationMatrix[,"gene2"]>0.1,]
      
      #Only proceed with adding this gene comparison to the list if more than 3 cells express it.
      if(nrow(correlationMatrix)<3) {next}
      
      #Only proceed with adding this gene comparison to the list if more than the required percentage of cells express GOI
      if(nrow(correlationMatrix) < (ncol(groupedSCS[[groupName]]@ndata)*(percentage/100))) {next}
      
      #Add this comparison to the list
      genes <- c(genes,colnames(inputData)[i])

      #Calculate & Add correlation
      correlation <- as.numeric(cor.test(correlationMatrix[,1],correlationMatrix[,2],method= correlationMethod)["estimate"])
      allCorrelations <- c(allCorrelations,correlation)
      
      #Calculate & add pValue
      pValue <- as.numeric(cor.test(correlationMatrix[,1],correlationMatrix[,2],method= correlationMethod)["p.value"])
      allPValues <- c(allPValues,pValue)
      
      #Add mean expression and number of cells with expression to the list
      expression <- mean(correlationMatrix[,"gene2"])
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
  
  print(paste0(
    nrow(inputData[inputData[,grepl(geneName, colnames(inputData))]>0.1,]) / length(cells)*100,
    " % of cells express gene '",geneName,"'"
  ))
  
  returnList = list(
    'geneName' = geneName,
    'groupName' = groupName,
    'cellNames' = cellNames,
    'removeNoExpressCells' = removeNoExpressCells,
    'percentage' = percentage,
    'correlationMethod' = correlationMethod,
    'correlationData' = correlationDataFrame
  )
  return(returnList)
}