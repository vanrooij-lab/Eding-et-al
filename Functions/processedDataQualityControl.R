### processedDataQualityControl(config, groupedSCS, groupNames=NULL, outputMode='save')
### 
### Authors
###   - Joep Eding
### 
### TODO Improve this comment
### - Switch to ggplot2
processedDataQualityControl <- function(config, groupedSCS, groupNames=NULL, outputMode='save') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop("You want to perform quality control for SCSObject '",groupName,"' that is not in groupedSCS")}
    }
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
  for(groupName in groupNames) {
    #Create directory if graphs need to be saved
    if(saveGraphs) {
      QCDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Quality control',.Platform$file.sep,'processed')
      if(!dir.exists(QCDir)) {dir.create(QCDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Open PDF file if required to save
    if(saveGraphs) {pdf(file=paste0(QCDir,.Platform$file.sep,groupName,".pdf"), width=12, height=10, pointsize=12)}
    
    #Plot within-cluster dispersion as a function of cluster number (only if sat==TRUE)
    plotsaturation(groupedSCS[[groupName]],disp=TRUE)
    
    #Plot change of within-cluster dispersion as a function of cluster number (only if sat==TRUE)
    plotsaturation(groupedSCS[[groupName]])
    
    #Plot silhouette of k-medoids clusters
    plotsilhouette(groupedSCS[[groupName]])
    
    #Plot Jaccard's similarity of k-medoids clusters
    plotjaccard(groupedSCS[[groupName]])
    
    #Plot barchart of outlier probabilities
    plotoutlierprobs(groupedSCS[[groupName]])
    
    #Plot regression of background model
    plotbackground(groupedSCS[[groupName]])
    
    #Plot dependence of outlier number on probability threshold
    plotsensitivity(groupedSCS[[groupName]])
    
    #Plot heatmap of k-medoids cluster
    clustheatmap(groupedSCS[[groupName]],final=FALSE,hmethod="single")
    
    #Plot heatmap of final cluster
    clustheatmap(groupedSCS[[groupName]],final=TRUE,hmethod="single")
    
    #Plot highlight of k-medoids clusters in t-SNE map
    plottsne(groupedSCS[[groupName]],final=FALSE)
    
    #Plot highlight of final clusters in t-SNE map
    plottsne(groupedSCS[[groupName]],final=TRUE)
    
    #Highlight cell labels in t-SNE map
    plotlabelstsne(groupedSCS[[groupName]],labels=sub("(\\_\\d+)","",names(groupedSCS[[groupName]]@ndata)))
    
    #Hightlight groups of cells by symbols in t-SNE map
    plotsymbolstsne(groupedSCS[[groupName]],types=sub("(\\_\\d+)$","", names(groupedSCS[[groupName]]@ndata)))  
    
    #Close PDF file if it was opened
    if(saveGraphs) {dev.off()}
  }
}
