### Function that plots cluster gene ontology results
### 
### Parameters
###   - 
### 
### Return
###   -
### 
### Authors
###   - Joep Eding
### 
### TODO: 
###   - Make this comment useful.
###   - If overrideDir is set, pass it along to plotGeneExpression
###   - Make an option to pass along includeOutlierCells for the tSNE. (now stuck on FALSE by default)
###   - Make it possible to influence tSNE-map properties more when plotRegulatedGenes=T
plotClusterGeneOntology <- function(config, groupedSCS, GOResults, clusterSpecificExpression, topX=NULL, plotRegulatedGenes=FALSE, overrideDir=FALSE, outputMode = 'pdf') {
  # Sanitize input outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  } else if(outputMode != 'show') {
    # Create directory if graphs need to be saved
    if(overrideDir == FALSE) {
      saveDir <- paste0(outputDir,.Platform$file.sep,'geneOntology')
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
    for(clusterName in names(GOResults[[groupName]])) {
      # Create workbook to save differential expression data of regulated genes per term
      regulatedGeneWB = createWorkbook()
      
      # Get data and make sure it's ordered with low->high pValue
      ontologyData = GOResults[[groupName]][[clusterName]]
      ontologyData = ontologyData[order(ontologyData$Pvalue, decreasing = FALSE),]
      
      # Write regulated terms to sheet
      addWorksheet(regulatedTermWB, clusterName)
      writeDataTable(regulatedTermWB, clusterName, ontologyData)
      
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
        title = paste0("Gene ontology ",groupName,': ',clusterName)
      )
      
      if(outputMode != 'show') {
        fileName = paste0("GO_",groupName,"_",clusterName,"_regulatedTerms.",outputMode)
        ggsave(fileName,path = saveDir)
        print(paste0("Saved plot as: ",saveDir,fileName))
      } else {
        print(ontologyPlot)
      }
      
      for(i in 1:nrow(ontologyData)) {
        # Reinstate list of genes regulated in this ontology term and make grep-ready
        listOfGenes = unique(unlist(str_split(ontologyData[i,'genes'], ';')))
        listOfGenes = unlist(lapply(listOfGenes, function(x){return(paste0('^',x,'__'))}))
                  
        # Plot tSNE map for regulated genes
        if(plotRegulatedGenes){ 
          plotGeneExpressionTSNE(config, groupedSCS, geneNames = listOfGenes, groupNames = groupName, includeOutlierCells = F, colorScheme = 'Spectral', expressionScale = 'linear', overrideDir = paste0(overrideDir,'/',groupName,'/',clusterName,'/',ontologyData[i,'Term']), outputMode = outputMode) 
        }
      
        # Fetch differential expression data of the genes regulated for this term
        DEdata = clusterSpecificExpression[[groupName]][[clusterName]] %>% tibble::rownames_to_column()
        DEdata = DEdata[which(grepl(paste0('(',paste(listOfGenes, collapse=')|('),')'), DEdata$rowname)),]
        DEdata = DEdata[order(DEdata$fc, decreasing = T),]
        
        # Add table with regulated genes to workbook
        sheetName = substr(ontologyData[i,'Term'], 1, 30)
        addWorksheet(regulatedGeneWB, sheetName) # Can't handle sheetnames > 31 characters
        writeDataTable(regulatedGeneWB, sheetName, DEdata)
      }
      regulatedGeneWBFilename = paste0(saveDir, '/GO_',groupName,'_',clusterName,'_regulatedGenes.xlsx')
      print(paste0("Saving ",regulatedGeneWBFilename))
      saveWorkbook(regulatedGeneWB, regulatedGeneWBFilename)
    }
    saveWorkbook(regulatedTermWB, paste0(saveDir, "/GO_",groupName,"_clusterGeneOntology.xlsx"))
  }
}