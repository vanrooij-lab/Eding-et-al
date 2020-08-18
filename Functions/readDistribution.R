
### readDistribution(config, groupedData, groupNames=NULL)
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedData: Pass the full groupedData-variable obtained from loadData.
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does QC for all.
###
### Return: list containing number of spike-ins, chromosomal reads and gene reads per group.
### 
### Authors
###   - Joep Eding
### 
### TODO
readDistribution <- function(config, groupedData, groupNames=NULL) {
  #Sanitize input groupName
  if(is.null(groupNames)) {
    groupNames = names(groupedData)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedData))) {stop("You want to analyze readDistribution for group '",groupName,"' that is not in groupedData")}
    }
  }
  
  #Calculate read distribution for each group
  readDistr <- list()
  for(groupName in groupNames) {
    readDistr[[groupName]]['totalReads'] <- sum(rowSums(groupedData[[groupName]]))
    readDistr[[groupName]]['ERCCReads'] <- sum(rowSums(groupedData[[groupName]][grepl('ERCC-',rownames(groupedData[[groupName]])),]))
    readDistr[[groupName]]['chrMReads'] <- sum(rowSums(groupedData[[groupName]][grepl('__chrM',rownames(groupedData[[groupName]])),]))
    readDistr[[groupName]]['geneReads'] <- sum(rowSums(groupedData[[groupName]][!grepl('__chrM|ERCC-',rownames(groupedData[[groupName]])),]))
    readDistr[[groupName]]['mitoRatio'] <- readDistr[[groupName]]['chrMReads']/readDistr[[groupName]]['geneReads']
  }
  
  #Return data
  return(readDistr)
}
