### analyseCorrelationAB(config, groupedSCS, geneNameA, geneNameB, groupNames=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### geneNames: Pass a list of exactly two gene names to plot the correlation for. First gene will go on X-axis, second on Y-axis.
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does analysis for all groups in groupedSCS.
### outputMode: Save will save the graphs to the correlation-directory. Show will just print the graph. (default: 'save')
### plotLine 
### Return: void. Either saves graphs or prints them to viewer. 
### 
### Authors
###   - Joep Eding
### 
### TODO
###  - Update comment.
analyseCorrelationAB <- function(config, groupedSCS, geneNames=NULL, correlationMethod = 'pearson', groupNames=NULL, removeOutlierCells=F, plotLine=T, overrideDir=F, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }
  
  #Sanitize input correlationMethod
  if(!(correlationMethod %in% c("pearson", "kendall", "spearman"))) {
    stop(paste0("You want tot analyze correlation using ",correlationMethod," which is not either 'pearson', 'kendall', or 'spearman'. Mind capitalisation."))
  }
  
  #Sanitize input geneNames
  inputGeneNames <- geneNames
  if(length(geneNames) != 2) {stop("You did not provide the correct number (exactly 2) of gene names to analyse correlation for")}
  for(geneName in geneNames) {
    for(groupName in groupNames) {
      #Find geneName in correct format
      newGeneName <- grep(geneName, rownames(groupedSCS[[groupName]]@ndata), value=T)
      
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
  
  #Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  #Rename geneName variables to be more convenient in use. 
  geneX = geneNames[1]
  geneY = geneNames[2]
  
  for(groupName in groupNames) {
    #Create directory if graphs need to be saved
    if(outputMode != 'show') {
      if(overrideDir == FALSE) {
        saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Correlation')
      } else {
        saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
      }
      if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    # Fetch expression data for the specified two genes
    cellData <- groupedSCS[[groupName]]@ndata[grepl(paste0('(',geneX,'|',geneY,')'),rownames(groupedSCS[[groupName]]@ndata)),]
    
    # Remove outlier cells if required
    if(removeOutlierCells) {
      cellData <- cellData[,!(colnames(cellData) %in% names(groupedSCS[[groupName]]@cpart[which(groupedSCS[[groupName]]@cpart>max(groupedSCS[[groupName]]@cluster$kpart))]))]
    }
    
    # Reform data for ggplot
    cellDataForGraph <- data.frame(
      X = unlist(cellData[geneX,]),
      Y = unlist(cellData[geneY,])
    )

    #Perform statistical correlation analysis
    correlationResult = cor.test(cellDataForGraph$X, cellDataForGraph$Y, alternative='two.sided', method=correlationMethod)
    expressionPercentageX = round(nrow(cellDataForGraph[cellDataForGraph$X>0.1,])/nrow(cellDataForGraph)*100,1)
    expressionPercentageY = round(nrow(cellDataForGraph[cellDataForGraph$Y>0.1,])/nrow(cellDataForGraph)*100,1)
    
    #Generate graph
    ABCorrelation <-  ggplot(   #Start ggplot with data and grouping
          cellDataForGraph,
          aes(x=X, y=Y)
        ) + geom_point(    #Print points
        ) + theme_classic(         #Get rid of default ggplot theme
        ) + theme(
          plot.subtitle=element_text(size=7, face="italic", color="black")
        ) + labs(
          x=inputGeneNames[1], 
          y=inputGeneNames[2],     #Get rid axis labels
          title=paste0(inputGeneNames[1],' vs. ',inputGeneNames[2]),
          subtitle=paste0(
            "group: ",groupName,'   ',
            "||   R: ",signif(correlationResult$estimate,3),'   ',
            "||   R-sq: ",signif(correlationResult$estimate^2,3),'   ',
            "||   P: ",signif(correlationResult$p.value,3),'   ',
            #"\n",
            '||   ',
            "Expr(",inputGeneNames[1],"): ",expressionPercentageX,'%   ',
            "||   Expr(",inputGeneNames[2],"): ",expressionPercentageY,'%   '
          )
        )
    if(plotLine) {
      ABCorrelation <- ABCorrelation + geom_smooth(
        method = 'lm',
        se=T,
        color='black'
      )
    }
    #Show or save graph
    if(outputMode != 'show') {
      ggsave(paste0(groupName,'_',inputGeneNames[1],'_vs_',inputGeneNames[2],'.',outputMode),path = saveDir)
    } else {
      print(ABCorrelation)
    }
  }
}