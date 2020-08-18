### analyseMarkerGenesPerCluster(config, clusterDifferentialExpression, geneNames=NULL, groupNames=NULL, outputMode='save')
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
analyseMarkerGenesPerCluster <- function(config, clusterDifferentialExpression, geneNames=NULL, groupNames=NULL, overrideDir=F, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to analyse cluster differential expression for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }

  #Sanitize input geneNames
  if(is.null(geneNames)) {
    geneNames = config$markerGenes
  }
  #Check whether the list of markergenes is present in every group
  for(geneName in geneNames) {
    for(groupName in groupNames) {
      #Find geneName in correct format
      newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@fdata), value=T)
      
      if(length(newGeneName) == 0) {
        #Gene not found
        stop(paste0("You want to analyze marker gene ",geneName," which is not expressed in group ",groupName))
      } else if(length(newGeneName) > 1) {
        #Gene name matches multiple
        stop(paste0("You want to generate a t-SNE for gene ",geneName," which is ambiguous. Did you mean: ", paste(newGeneName, collapse = " OR "),"?"))
      } else {
        #Gene is unique, update identifier
        geneNames[match(geneName, geneNames)] = newGeneName
      }
    }
  }
    
  #Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  #Analyze marker genes for each group
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
    
    #Make a report of marker gene expression for each cluster
    for (clusterName in names(clusterDifferentialExpression[[groupName]])) {
      #Make graph
      markerGenes <- ggplot(
      ) + geom_col(
        data=clusterDifferentialExpression[[groupName]][[clusterName]][grepl(paste(geneNames, collapse = '|'), row.names(clusterDifferentialExpression[[groupName]][[clusterName]])),],
        aes(
          x=row.names(clusterDifferentialExpression[[groupName]][[clusterName]][grepl(paste(geneNames, collapse = '|'), row.names(clusterDifferentialExpression[[groupName]][[clusterName]])),]),
          y=mean.cl
        ),
        color='red'
      ) + theme_classic(
      ) + theme(
        title = config$plotOptions[['title']],
        #subtitle = config$plotOptions[['subtitle']],
        axis.title = config$plotOptions[['axis.title']],
        axis.text = config$plotOptions[['axis.text']],
        axis.text.x = config$plotOptions[['axis.text.x']],
        axis.text.y = config$plotOptions[['axis.text.y']],
        legend.position = 'right',
        legend.title = config$plotOptions[['legend.title']],
        legend.text = config$plotOptions[['legend.text']],
        aspect.ratio = config$plotOptions[['aspect.ratio']],
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = config$plotOptions[['axis.line']]
      ) + labs(
        x="Genes",
        y="Mean transcript count in cluster",
        title=paste0(groupName,' - ',clusterName," - marker genes")
      )
      
      #Show or save graph
      if(outputMode == 'show') {
        print(markerGenes)
      } else {
        fileName = paste0(groupName,"_markerGenes_",clusterName,'.',outputMode)
        ggsave(fileName, path=saveDir)
        print(paste0("t-SNE map saved as ",saveDir,.Platform$file.sep,fileName))
      }
    }
  }
}
