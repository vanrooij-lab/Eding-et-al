
### geneQuality
### parameters
### 
### Return: void. Generates a single graph that is either saved or showed on screen.
### 
### Authors
###   - Joep Eding
### 
### TODO:
###  - Make this comment useful?
###  - Think about coloring options?
geneQuality <- function(config, groupedData, groupedSCS, groupName, overrideDir=NULL, overrideTheme=NULL, outputMode='pdf') {
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
  
  # Get readcounts and process for coloring
  readCounts <- data.frame(
    meanCounts = log2(rowMeans(groupedData[[groupName]][,colnames(groupedSCS[[groupName]]@ndata)]+1)),
    numCells = as.numeric(apply(groupedData[[groupName]][,colnames(groupedSCS[[groupName]]@ndata)], 1, function(x) {length(which(x > 1))}))
  )
   
  #Plot graph
  geneQuality = ggplot(
    readCounts, 
    aes(x=meanCounts, y=numCells)
  ) + geom_point(
  ) + theme_classic(
  ) + theme(
    aspect.ratio = config$plotOptions$aspect.ratio,
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank()
  ) + labs(
    y="Number of cells \nwith expression",
    x='Log2(mRNA per cell + 1)',
    title=paste0("Gene quality - ",groupName)
    #color=paste0(expressionScaleTag,"\nExp")
  ) + overrideTheme
  
  if(outputMode == 'show') {
    print(geneQuality)
  } else {
    fileName = paste0(groupName,'_geneQuality.',outputMode)
    ggsave(fileName, path=saveDir, useDingbats=FALSE)
    print(paste0("Plot saved as ",saveDir,.Platform$file.sep,fileName))
  }
  
}