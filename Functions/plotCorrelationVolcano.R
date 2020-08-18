### plotCorrelationVolcano
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters
### correlationResult: Pass the output from the function analyseCorrelation or analyseMetaDataCorrelation, this function will automatically extract relevant parameters from that.
### coefficientCutoff: Pass a list of two doubles (between 0 & 1), correlation coefficient needs to be outside these two values to be considered significant
### pValueCutoff: Pass a double (between 0 & 1), pValue for correlation needs to be smaller than cutoff to be considered significant.
### overrideDir: 
### outputMode: Save will save the graphs to the correlation-directory. Show will just print the graph. (default: 'save')
### 
### Authors
###   - Joep Eding
### 
### TODO 
### - Improve this comment
### - Improve label overriding: Make it possible to also show labels for genes that don't match significance cutoffs.
### Return: 
plotCorrelationVolcano <- function(config, correlationResult=NULL, coefficientCutoff=c(-0.3, 0.3), pValueCutoff=0.05, overrideLabels=FALSE, overrideDir=F, outputMode='pdf') {
  # Sanitize input correlationResult
  if(is.null(correlationResult)) {stop("You did not provide a data.frame with a correlation result.")}
  
  # Exclude genes for which analyseCorrelation failed
  volcanoInput = correlationResult$correlationData[!is.na(correlationResult$correlationData[,'correlation']),]
  
  # Exclude the gene of interest
  volcanoInput = volcanoInput[volcanoInput$genes != correlationResult$geneName, ]
  
  # Sanitize input overrideLabels
  if(!isFALSE(overrideLabels)) {
    if(!is.character(overrideLabels)) {
      stop("Inappropriate input for parameter 'overrideLabels'. Use a list of regex-ready gene names.")
    }
  }
  
  # Sanitize input coefficientCutoff
  if(length(coefficientCutoff) != 2) {
    stop("Parameter 'coefficientCutoff' needs to be a list of exactly 2 numerical values. Default: 'c(-0.3, 0.3)'.")
  } else if(!is.numeric(coefficientCutoff[1]) | !is.numeric(coefficientCutoff[2])) {
    stop("Parameter 'coefficientCutoff' needs to be a list of exactly 2 numerical values. Default: 'c(-0.3, 0.3)'.")
  } else if(any(coefficientCutoff > 1) | any(coefficientCutoff < -1)) {
    stop("Both values of 'coefficientCutoff' need to be between -1 & 1. Default: 'c(-0.3, 0.3)'.")
  } else if(coefficientCutoff[1] == coefficientCutoff[2]) {
    stop("The two values in parameter 'coefficientCutoff' need to be different from each other. Default: 'c(-0.3, 0.3)'.")
  } else if(coefficientCutoff[1] > coefficientCutoff[2]) {
    stop("The first value in parameter 'coefficientCutoff' needs to be smaller than the second value. Default: 'c(-0.3, 0.3)'.")
  }
  
  #Calculate boundaries for the graph
  if(max(volcanoInput[,"correlation"])>abs(min(volcanoInput[,"correlation"]))){
    limit <- max(volcanoInput[,"correlation"])
  } else {
    limit <- abs(min(volcanoInput[,"correlation"]))
  }
  if(limit < max(abs(coefficientCutoff))) {
    limit = max(abs(coefficientCutoff))
  }
  limit <- limit + 0.05
  
  #Collect positively correlated genes
  positivelyCorrelatedGenes  <- volcanoInput[volcanoInput[,"correlation"]>=coefficientCutoff[2]&
                                            volcanoInput[,"pValue"]<=pValueCutoff,]
  #Collect negatively correlated genes
  negativelyCorrelatedGenes  <- volcanoInput[volcanoInput[,"correlation"]<=coefficientCutoff[1]&
                                            volcanoInput[,"pValue"]<=pValueCutoff,]
  
  #Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  #Create directory if graphs need to be saved
  if(outputMode != 'show') {
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,correlationResult$groupName,.Platform$file.sep,'Correlation')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,correlationResult$groupName)
    }
    if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
  }
  
  #Add grouping to volcanoInput based on whether genes make the cutoff for strength & significance of correlation
  volcanoInput$group = 'Neutral'
  volcanoInput[(volcanoInput$correlation>=coefficientCutoff[2] & volcanoInput$pValue<=pValueCutoff), 'group'] = 'Up'
  volcanoInput[(volcanoInput$correlation<=coefficientCutoff[1] & volcanoInput$pValue<=pValueCutoff), 'group'] = 'Down'
  #Add labels for cells that make the cutoff
  volcanoInput$label = ''
  volcanoInput[volcanoInput$group != 'Neutral','label']=volcanoInput[volcanoInput$group != 'Neutral','genes']
  #Generate a named list of colors
  namedColors = c('black', 'red3', 'royalblue3')
  names(namedColors) = c('Neutral', 'Up', 'Down')
  #Make sure gene names are in character type
  # TODO Figure out why this is necessary, maybe this information is being supplied wrong?
  volcanoInput$genes <- as.character(volcanoInput$genes)
  
  #Set tag for exclusion of cells based on removeNoExpressCells in the correlationResult
  if(correlationResult$removeNoExpressCells == 'no') {
    removeNoExpressCellsTag = 'All cells & all genes were used for this analysis.'
  } else if(correlationResult$removeNoExpressCells == 'yes') {
    removeNoExpressCellsTag = paste0("Only cells expressing both ",strsplit(correlationResult$geneName,'__')[[1]][1],
                                    " and the comparison gene were included.\nGenes excluded if less than ",
                                    correlationResult$percentage,"% of cells express them.")
  } else if(correlationResult$removeNoExpressCells == 'goi') {
    removeNoExpressCellsTag = paste0("Only cells expressing '",strsplit(correlationResult$geneName,'__')[[1]][1],
                                "' were included.\nGenes excluded if less than ",
                                correlationResult$percentage,"% of cells express it.")
  } else if(correlationResult$removeNoExpressCells == 'parameter') {
    removeNoExpressCellsTag = paste0("Only cells for which the metadata parameter was unavailable were excluded")
  }
  
  if(isFALSE(overrideLabels)) {
    volcanoLabelInput <- volcanoInput
  } else {
    volcanoLabelInput <- volcanoInput %>% dplyr::filter(stringr::str_detect(genes, paste0('(',paste(overrideLabels, collapse=')|('),')')))
  }
  
  # Generate graph
  correlationVolcano <- ggplot( #Load data
    data=volcanoInput,
    aes(
      x=correlation,
      y=-log10(pValue)
    )
  ) + geom_point( #Plot points
    aes(color=group)
  ) + scale_color_manual( #Color according to significantly correlated
    values=namedColors
  ) + xlim( #Set plot limits
    -limit,
    limit
  ) + geom_vline( #Draw lines for cutoff for positive or negative correlation
    aes(xintercept=coefficientCutoff[1]),color='royalblue3'
  ) + geom_vline(
    aes(xintercept=coefficientCutoff[2]),color='red3'
  ) + geom_hline(
    aes(yintercept=-log10(pValueCutoff)), color='black'
  ) + geom_text_repel( #Draw labels for positively correlated genes
    aes(label=lapply(genes, function(x){strsplit(x,'__')[[1]][1]})),
    data=volcanoLabelInput[which(volcanoLabelInput$correlation>=coefficientCutoff[2] & volcanoLabelInput$pValue<=pValueCutoff),],
    #data=volcanoLabelInput[which(volcanoLabelInput$correlation>=0),],
    size=3,
    #label.padding=0.1,
    segment.size = 0.3,
    min.segment.length = if(isFALSE(overrideLabels)) 0.5 else 0,
    xlim=c(coefficientCutoff[2], NA)
  ) + geom_text_repel( #Draw labels for negatively correlated genes
    aes(label=lapply(genes, function(x){strsplit(x,'__')[[1]][1]})),
    data=volcanoLabelInput[which(volcanoLabelInput$correlation<=coefficientCutoff[1] & volcanoLabelInput$pValue<=pValueCutoff),],
    #data=volcanoLabelInput[which(volcanoLabelInput$correlation<=0),],
    size=3,
    #label.padding=0.1,
    segment.size = 0.3,
    min.segment.length = if(isFALSE(overrideLabels)) 0.5 else 0,
    xlim=c(NA, coefficientCutoff[1])
  ) + theme_classic(
  ) + theme( #Remove legend, make subtitle smaller
    legend.position = 'none',
    plot.subtitle=element_text(size=7, face="italic", color="black"),
    axis.title=element_text(size=10),
    axis.text=element_text(size=8)
    #aspect.ratio = 0.40
  ) + labs( #Add titles and axis labels.
    x="Correlation coefficient",
    y="-log10(pValue)",
    title = paste0("Correlation analysis for ",strsplit(correlationResult$geneName,'__')[[1]][1]),
    subtitle = paste0("group: ",correlationResult$groupName,"        p-value cutoff: ",pValueCutoff,"        correlation cutoffs: ",coefficientCutoff[1]," & ",coefficientCutoff[2],"\n",removeNoExpressCellsTag)
  )
  #Show or save graph
  if(outputMode != 'show') {
      ggsave(paste0(correlationResult$groupName,'_',strsplit(correlationResult$geneName,'__')[[1]][1],'_',correlationResult$removeNoExpressCells,correlationResult$percentage,'_correlationVolcano.',outputMode),path = saveDir, useDingbats=FALSE)
      #ggsave(paste0(correlationResult$groupName,'_',strsplit(correlationResult$geneName,'__')[[1]][1],'_',correlationResult$removeNoExpressCells,correlationResult$percentage,'_correlationVolcano.',outputMode),path = saveDir, useDingbats=FALSE, width=15, height=25, units='cm')
  } else {
      print(correlationVolcano)
  }
  
  #return(volcanoInput)
}
