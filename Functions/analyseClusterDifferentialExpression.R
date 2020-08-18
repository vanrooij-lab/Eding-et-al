### analyseClusterDifferentialExpression
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
### 
### Return: List with differential expression per cluster per group
analyseClusterDifferentialExpression <- function(config, groupedSCS, groupNames=NULL, overrideDir=F, outputMode='save') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop("You want to analyse cluster differential expression for SCSObject '",groupName,"' that is not in groupedSCS")}
    }
  }
  
  #Check outputMode is set to save. Not saving makes no sense for this function.
  if(!(outputMode %in% c('save','show'))) {stop("You tried to analyze cluster differential expression with outputMode set to something else than 'save', that won't work for this function.")}
  
  #Initialize a return variable
  diffExpReturn <- list()
  
  #Analyse clusterDifferentialExpression for all groups
  for(groupName in groupNames) {
    #Create directory if necessary
    if(outputMode != 'show') {
      if(overrideDir == FALSE) {
        saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Clusters')
      } else {
        saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
      }
      if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Initialize Excel workbook to store differential expression in. 
    diffExpWB <- createWorkbook()
    
    #Generate differential expression lists
    diffExp <- clustdiffgenes(groupedSCS[[groupName]],pvalue=1) 
    
    #Process differential expression in each cluster
    for(clusterName in names(diffExp)) {
      #Create new sheet in the workbook for this cluster
      addWorksheet(diffExpWB, sheetName=paste0("Cluster ",substring(clusterName,4)))
      
      #Reorder differential expression data based on fold change
      diffExp[[clusterName]] <- diffExp[[clusterName]] %>% tibble::rownames_to_column(var='geneName') %>% dplyr::arrange(desc(fc)) %>% tibble::column_to_rownames(var='geneName')
      #Write dataframe to the new sheet
      writeDataTable(
        diffExpWB,
        sheet = paste0("Cluster ",substring(clusterName,4)),
        data.frame(
          GeneID=rownames(diffExp[[clusterName]]),
          diffExp[[clusterName]]
        )
      )
    }
    
    #Save workbook if outputMode is not 'show'
    if(outputMode != 'show') {
      saveWorkbook(diffExpWB, paste0(saveDir,.Platform$file.sep,'clusterSpecificExpression.xlsx'), overwrite = TRUE) 
    }
    
    #Add differential expression to return variable
    diffExpReturn[[groupName]] <- diffExp
  }
  
  #Return differential expression
  return(diffExpReturn)
}