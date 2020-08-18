### returnTopDifferentialGenePerCluster(config, groupedSCS)
### TODO Make a sensible comment here.
### 
### Parameters
###   -
### 
### Return
###   -
### 
### Authors
###   - Joep Eding
### 
### TODO
###   - 
returnTopDifferentialGenesPerCluster <- function(config, groupedSCS, clusterDifferentialExpression, numGenes=5) {
  # TODO Add input sanitization here
  
  # Initialize return variable
  returnVariable = list()
  
  # Loop through groups to fill returnVariable
  for(groupName in names(groupedSCS)) {
    #Find number of clusters:
    numClusters = max(groupedSCS[[groupName]]@cpart)
  
    geneList = list()
    for(clusterNum in 1:numClusters) {
      clustDiff = clusterDifferentialExpression[[groupName]][[paste0('cl.',clusterNum)]]
      clustDiff = clustDiff[order(clustDiff$fc, decreasing = T),]
      geneList[[clusterNum]] = rownames(clustDiff)[1:numGenes]
    }
    
    returnVariable[[groupName]] = geneList
  }
  
  return(returnVariable)
}