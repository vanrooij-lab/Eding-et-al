### loadData(config, groups)
### Returns a named list where each element is named according to the groups to analyse and contains a merged data.frame of all transcriptCounts files for that group.
### 
### Parameters
###   - config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
###   - groupsToAnalyse: Pass a list with the names of the groups that you want to analyse. Pass names(config$groups) to analyse all groups. 
### 
### Return
###   - 
### 
### Authors
###   - Joep Eding
###
### TODO
###   - 
loadData <- function(config, groupsToAnalyse) {
  #Check whether groups contains valid group names.
  for(groupName in groupsToAnalyse) {
    if(!(groupName %in% names(config$groups))) {
      stop(paste0("Trying to load data for group ",groupName," that's not defined in config$groups."))
    }
    #Check whether all samples defined in this group are defined in config$samples
    for(sampleName in config$groups[[groupName]]) {
      if(!(sampleName %in% names(config$samples))) {
        stop(paste0("Trying to load data for sample ",sampleName," in group ",groupName," that is not defined in config$samples"))
      }
    }
  }
  
  #Load all data files that will be used during this analysis
  datafiles <- list()
  for(fileName in dir()) {
    #Check whether it's a TranscriptCounts-file, skip to next file otherwise
    if(!grepl('TranscriptCounts.tsv',fileName)) { next }
    
    #Get the plate-specific part of the name
    plateName <- strsplit(fileName,'_TranscriptCounts.tsv')[[1]]
    
    #Check whether this plate is required to be loaded, skip to next file otherwise
    if(!(plateName %in% unlist(config$sample[unlist(config$groups[groupsToAnalyse])]))) { next }
    
    #Read file and save to a temporary list
    datafiles[[plateName]] <- read.table(fileName, header=T, row.names=1, sep="\t")
    
    #Make the column names nice
    colnames(datafiles[[plateName]]) <- generateColumnNames(plateName)
  }
  
  #Check whether all requested plates have been loaded, stop otherwise
  for(groupName in groupsToAnalyse) {
    for(plateName in unlist(config$sample[unlist(config$groups[[groupName]])])) {
      if(!(plateName %in% names(datafiles))) {stop(paste0("Plate with name ",plateName," in group ",groupName," was not loaded"))}
    }
  }
  
  #For every group to analyse, merge all relevant plates into a single dataframe
  groupedData <- list()
  for(groupName in groupsToAnalyse) {
    #Merge all relevant plates for this group into a single dataframe
    tempData = Reduce(
      function(x,y) {
        x <- merge(x,y,by.x=0,by.y=0,all=TRUE)
        rownames(x) <- x[,1]
        x[,"Row.names"] <- NULL
        return(x)
      },
      datafiles[which(names(datafiles) %in% unlist(config$sample[unlist(config$groups[[groupName]])]))]
    )
    
    #Set NAs to 0 as the further processing logic can't handle NAs.
    tempData[is.na(tempData)] <- 0
    
    #Add this group to the groupedData list
    groupedData[[groupName]] = tempData
  }
  
  #Clear temporary data to free up memory
  remove('tempData', 'datafiles')
  
  #Return groupedData
  return(groupedData)
}