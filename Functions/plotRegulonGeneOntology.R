### Plot regulon gene ontology
### 
### Parameters
###   - 
### 
### Return
###   -
### 
### Authors
###   - Martijn Wehrens - Original idea and script
###   - Joep Eding - Adaptation for accessibility and uniformity with other functions
### 
### TODO
###   - 
plotRegulonGeneOntology <- function(config, groupedSCS, GOResults, corrResults, topX=NULL, plotRegulatedGenes=FALSE, overrideDir=FALSE, outputMode = 'pdf') {
  # Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  } else if(outputMode != 'show') {
    # Create directory if graphs need to be saved
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,'regulons/ontology/')
    } else {
      saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir)
    }
    if(!dir.exists(saveDir)) {
      dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print(paste0("Created directory: ",saveDir))
    }
  }
  
  # Iterate over groups & clusters
  for(groupName in names(GOResults)) {
    regulatedTermWB = createWorkbook()
    for(regulonName in names(GOResults[[groupName]])) {
      # Get data and make sure it's ordered with low->high pValue
      ontologyData = GOResults[[groupName]][[regulonName]]
      ontologyData = ontologyData[order(ontologyData$Pvalue, decreasing = FALSE),]
      
      # Write regulated terms to sheet
      addWorksheet(regulatedTermWB, regulonName)
      writeDataTable(regulatedTermWB, regulonName, ontologyData)
      
      # Take only topX enriched terms
      if(!is.null(topX)) {
        if(topX < nrow(ontologyData)) {
          ontologyData=ontologyData[1:topX,]
        }
      }
      
      # Plot gene ontology
      ontologyPlot <- ggplot(
        data=ontologyData,
        aes(x=reorder(Term,-Pvalue),y=-log10(Pvalue))
      ) + geom_bar(
        stat='identity'
      ) + coord_flip(
      ) + labs( #Add titles and axis labels.
        x="Term",
        y="-log10(pValue)",
        title = paste0("Gene ontology ",groupName,': ',regulonName,' correlated genes')
      )
      
      if(outputMode != 'show') {
        fileName = paste0("GO_",groupName,"_",regulonName,"_regulatedTerms.",outputMode)
        ggsave(fileName,path = saveDir)
        print(paste0("Saved plot as: ",saveDir,fileName))
      } else {
        print(ontologyPlot)
      }
    }
    saveWorkbook(regulatedTermWB, paste0(saveDir, "/","GO","_",groupName,"_RegulonGeneOntology.xlsx"))
  }
} 