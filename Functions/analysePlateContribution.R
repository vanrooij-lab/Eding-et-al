
### analysePlateContribution(config, groupedSCS)
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### 
### Return: Named list containing the number of cells included from each plate
### 
### Authors
###   - Joep Eding
### 
### TODO
analysePlateContribution <- function(config, groupedSCS) {
  #Initialize return variable
  plateContribution <- list()
  
  #Iterate over all groups in groupedSCS
  for(groupName in names(groupedSCS)) {
    #Iterate over all cells in this group, count how many cells from each plate contribute to this group.
    for(cellName in colnames(groupedSCS[[groupName]]@fdata)) {
      if(is.null(plateContribution[[groupName]][substr(cellName,1,3)][[1]])) {plateContribution[[groupName]][substr(cellName,1,3)]=0}
      plateContribution[[groupName]][substr(cellName,1,3)] = as.numeric(plateContribution[[groupName]][substr(cellName,1,3)])+1
    }
  }
  
  #Return plateContribution
  return(plateContribution)
}