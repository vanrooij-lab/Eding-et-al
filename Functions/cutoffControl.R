### cutoffControl(config, groupedData, groupName, QCDir, saveGraphs)
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedData: Pass the full groupedData-variable obtained from loadData.
### groupName: Pass the name of the group to analyse the cutoff value for.
### TODO: Update parameters in comment for documentation
### 
### Authors
###   - Joep Eding
### 
### Return: void. Prints some output to console, additionally generates graphs that are either saved to PDF or printed to Viewer.
cutoffControl <- function(config, groupedData, groupName, numBins=100, cutoff=1000, overrideDir=NULL, overrideTheme=NULL, outputMode='pdf') {
  # Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf'. Pay attention to capitalisation.")
  } else if(outputMode != 'show') {
    # Create directory if graphs need to be saved
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,'QA',.Platform$file.sep)
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep)
    }
    if(!dir.exists(saveDir)) {
      dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print(paste0("Created directory: ",saveDir))
    }
  }
  
  # Sanitize input groupName
  if(!(groupName %in% names(groupedData))) {stop(paste0("You want to analyze cutoffControl for group '",groupName,"' that is not in groupedData"))}
  
  # Sanitize input cutoff
  if(!is.numeric(cutoff)) {
    stop("You need to supply a number for parameter 'cutoff'")
  }
  
  # Get readcounts and process for coloring
  readCounts = data.frame(
    counts = log10(colSums(groupedData[[groupName]])+1),
    color = unlist(lapply(as.numeric(colSums(groupedData[[groupName]])), function(x) {if(x >= cutoff) {'Over'} else {'Under'}}))
  )
  
  #Plot graph
  readsHistogram = ggplot(
    readCounts, 
    aes(x=counts, fill=color)
  ) + geom_histogram(
    bins = 200,
    position = 'stack',
    show.legend = F
  ) + theme_classic(
  ) + theme(
    aspect.ratio = config$plotOptions$aspect.ratio,
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank()
  ) + labs(
    y='Number of cells',
    x='Log10(transcript count per section + 1)',
    title=paste0("Transcriptome quality - ",groupName)
    #color=paste0(expressionScaleTag,"\nExp")
  ) + geom_vline(
    xintercept = log10(cutoff)
  ) + overrideTheme
  
  if(outputMode == 'show') {
    print(readsHistogram)
  } else {
    fileName = paste0(groupName,'_cutoffControl.',outputMode)
    ggsave(fileName, path=saveDir, useDingbats=FALSE)
    print(paste0("Histogram saved as ",saveDir,.Platform$file.sep,fileName))
  }
}