### rawDataQualityControl(config, groupedData, groupToAnalyse=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedData: Pass the full groupedData-variable obtained from loadData.
### groupsToAnalyse: Optionally, pass the name of a single group to perform quality control for or pass a list of names. By default (=NULL) does QC for all.
### outputMode: Pass either 'save' (saves output to PDF) or 'show' (prints output to Viewer). 'save' is default.
### 
### Return: void. Prints some output to console, additionally generates graphs that are either saved to PDF or printed to Viewer.
### 
### 
### Authors
###   - Joep Eding
### 
rawDataQualityControl <- function(config, groupedData, groupsToAnalyse=NULL, outputMode='save') {
  #If groupsToAnalyse is not defined, set it to contain all groups.
  if(is.null(groupsToAnalyse)) {
    groupsToAnalyse=names(groupedData)
  } 
  
  #Check if all groupsToAnalyse exist
  for(groupName in groupsToAnalyse) {
    if(!(groupName %in% names(groupedData))) {stop(paste0("Asked to do quality control for group ",groupName," for which no data was loaded"))}
  }
  
  #Set saveGraphs to a boolean based on outputMode
  if(outputMode == 'save') {
    saveGraphs <- TRUE
  } else if(outputMode == 'show') {
    saveGraphs <- FALSE
  } else {
    stop("outputMode needs to be either 'save' or 'show', no other values accepted, pay attention to capitalisation.")
  }
  
  #Perform quality control for all the groups
  for(groupName in groupsToAnalyse) {
    #Create directory if graphs need to be saved
    if(saveGraphs) {
      QCDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Quality control',.Platform$file.sep,'raw')
      if(!dir.exists(QCDir)) {dir.create(QCDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Run cutoffControl
    cutoffControl(config, groupedData, groupName, QCDir, saveGraphs)
  }
}