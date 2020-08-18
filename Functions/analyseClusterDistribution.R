
### analyseClusterDistribution(config, groupedSCS, groupNames=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does analysis for all.
### outputMode: Pass 'save' (saves output to Excel file). Other values won't be accepted.
### 
### Return: void. Saves an Excel-file for each group containing which cluster each cell belongs to as well as the size of each cluster.
### 
### Authors
###   - Joep Eding
### 
### TODO
###   -
analyseClusterDistribution <- function(config, groupedSCS, groupNames=NULL, overrideDir=F, outputMode='save') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop("You want to analyse cluster size for SCSObject '",groupName,"' that is not in groupedSCS")}
    }
  }
  
  #Check outputMode is set to save. Not saving makes no sense for this function.
  if(outputMode != 'save') {stop("You tried to analyze cluster size with outputMode set to something else than 'save', that won't work for this function.")}
  
  #Analyze clusterDistribution for every groupName
  for(groupName in groupNames) {
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Clusters')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
    }
    #Create directory if necessary
    if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    
    #Initialize Excel workbook to store differential expression in. 
    clusterDistribution <- createWorkbook()
    
    #Generate list of which cell goes in which cluster
    cellClusters <- data.frame(CellID=names(groupedSCS[[groupName]]@cpart),cluster=groupedSCS[[groupName]]@cpart)
    cellClusters <- cellClusters[order(cellClusters$cluster, decreasing=F),]
    
    #Add a sheet containing which cell goes in which cluster
    addWorksheet(clusterDistribution, sheetName='Cells')
    #Write cell destination to the new sheet
    writeDataTable(
      clusterDistribution,
      sheet = 'Cells',
      cellClusters
    )
    
    #Determine how much cells are in each cluster
    clusterSize <- NULL
    for(i in seq(1,max(cellClusters[,2]))){
      x <- NULL
      x <- length(cellClusters[cellClusters[,2]==i,2])
      clusterSize <- c(clusterSize,x)
    }
    clusterSize <- data.frame(seq(1,max(cellClusters[,2])),clusterSize)
    
    #Add a sheet containing which cell goes in which cluster
    addWorksheet(clusterDistribution, sheetName='ClusterSize')
    #Write clusterSize to this new sheet
    writeDataTable(
      clusterDistribution,
      sheet = 'ClusterSize',
      clusterSize
    )
    
    #Save workbook
    saveWorkbook(clusterDistribution, paste0(saveDir,.Platform$file.sep,'clusterDistribution.xlsx'), overwrite = TRUE)
  }
}