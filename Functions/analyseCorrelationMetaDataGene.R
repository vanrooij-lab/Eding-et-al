### analyseCorrelationMetaDataGene(config, groupedSCS, geneNameA, geneNameB, groupNames=NULL, outputMode='save')
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedSCS: Pass the full groupedSCS-variable obtained from buildSCSObject
### parameterName: Name of the metadata parameter to analyse
### geneNames: Pass a list of exactly two gene names to plot the correlation for. First gene will go on X-axis, second on Y-axis.
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) does analysis for all groups in groupedSCS.
### excludeNoExprCells: 
### overrideDir: 
### outputMode: Save will save the graphs to the correlation-directory. Show will just print the graph. (default: 'save')
### 
### Return: void. Either saves graphs or prints them to viewer. 
### 
### Authors
###   - Joep Eding
### 
### TODO: 
###   - Make excluding no expression cells optional
###   - Fix parameters in comment
analyseCorrelationMetaDataGene <- function(config, groupedSCS, parameterName=NULL, geneName=NULL, correlationMethod = 'pearson', groupNames=NULL, excludeNoExprCells=T, removeOutlierCells=F, plotLine=T, pointSizeModifier = 1, overrideDir=F, outputMode='pdf') {
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
  
  #Sanitize input geneName
  inputGeneName <- geneName
  if(length(geneName) != 1) {stop("You did not provide the correct number (exactly 1) of gene names to analyse correlation for")}
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
      geneName = newGeneName
    }
  }
  
  # Check whether a metadata object exists
  if(!exists('scsMetadata')) {
    stop("Load metadata first by running loadMetaData().")
  }
  
  #Sanitize input parameterName
  if(is.null(parameterName)) {stop("You did not provide a metadata parameter name to analyse correlation for.")}
  if(!(parameterName %in% colnames(scsMetadata))) {stop("Parameter '",parameterName,"' is not loaded")}
  
  #Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
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
    
    #Fetch expression data 
    cellData <- data.frame(
      stringsAsFactors = F,
      Cell = colnames(groupedSCS[[groupName]]@ndata),
      Gene = as.numeric(groupedSCS[[groupName]]@ndata[geneName,])
    )
    
    # Merge in metadata parameter data
    cellData <- merge(
      cellData,
      scsMetadata[,c('Cell',parameterName)],
      by='Cell',
      all.x=T,
      all.y=F
    )
    colnames(cellData) <- c('Cell', 'Gene', 'Parameter')
    
    # Filter for rows that have a value for Parameter
    cellData <- cellData %>% dplyr::filter(!is.na(Parameter))
    expressionPercentageGene = round(nrow(cellData[cellData$Gene>0.1,])/nrow(cellData)*100,1)
    if(excludeNoExprCells) {
      cellData = cellData %>% dplyr::filter(Gene > 0.1)
    }
    
    # Remove outlier cells if required
    if(removeOutlierCells) {
      cellData <- cellData[!(cellData$Cell %in% names(groupedSCS[[groupName]]@cpart[which(groupedSCS[[groupName]]@cpart>max(groupedSCS[[groupName]]@cluster$kpart))])),]
    }
    
    # Skip this group if there's not enough cells that have metadata
    if(nrow(cellData) < 3) {
      print(paste0("Not enough metadata available for parameter '",parameterName,"' and group '",groupName,"', skipping."))
      next
    }
    
    #Perform statistical correlation analysis
    correlationResult = cor.test(cellData$Parameter, cellData$Gene, alternative='two.sided', method=correlationMethod)
    
    #Generate graph
    ABCorrelation <-  ggplot(   #Start ggplot with data and grouping
          cellData,
          aes(x=Parameter, y=Gene)
        ) + geom_point(    #Print points
          size=1.5*pointSizeModifier
        ) + theme_classic(         #Get rid of default ggplot theme
        ) + theme(
          plot.subtitle=element_text(size=7, face="italic", color="black")
        ) + labs(
          x=parameterName, 
          y=inputGeneName,     #Get rid of axis labels
          title=paste0(parameterName,' vs. ',inputGeneName),
          subtitle=paste0(
            "group: ",groupName,'   ',
            "||   R: ",signif(correlationResult$estimate,3),'   ',
            "||   R-sq: ",signif(correlationResult$estimate^2,3),'   ',
            "||   P: ",signif(correlationResult$p.value,3),'   ',
            #"\n",
            '||   ',
            "Expr(",inputGeneName,"): ",expressionPercentageGene,'%   '
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
      ggsave(paste0(groupName,'_',parameterName,'_vs_',inputGeneName,'.',outputMode),path = saveDir)
    } else {
      print(ABCorrelation)
    }
  }
}
