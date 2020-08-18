### Function that performs gene ontology analysis based on correlation analysis results
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
### TODO
###   - Make this comment useful
###   - Consider storing the GOFrame globally to reduce computing time for repeated function calls.
###     - Save to a global variable, and try to see if it is defined before computing it. 
###     - Also allow option to override loading the existing set for whatever reason it may be necessary.
###   - Make this work seamlessly with KEGG analysis too. 
###   - Save xlsx of the ontology analysis too
###   - Test with data where config$geneIdentifierType != 'external_gene_name(_with_chr)'
###   - Remove requirement for corrCutoff to be doubly specified with posNeg == 'Both'
###   - Consider changing background gene selection to use number of cells with expression rather than mean expression over all cells to better match with the selection in the correlation analysis function. 
###   - Make a more flexible way of setting geneNames (now uses rownames from the first group's ndata)
analyseCorrelationGeneOntology <- function(config, groupedSCS, correlationResult, minExpression=0, corrCutoff=NULL, pCutoff=NULL, posNeg='Pos', pathwayPCutoff=0.05, GOKegg='GO', includeChildTerms=F) {
  # Sanitize input upDown and cutoffs.
  if(!(posNeg %in% c('Pos','Neg','Both'))) {
    stop("Parameter 'upDown' should be one of 'Pos', 'Neg', or 'Both'; Pay attention to capitalisation.")
  }
  if(all(is.null(corrCutoff), is.null(pCutoff))) {
    stop("Set a cutoff to determine significance in at least one of 'corrCutoff' and 'pCutoff'")
  }
  if(posNeg == 'Both') {
    if(!is.null(corrCutoff)) {
      if(length(corrCutoff) != 2) {
        stop("When posNeg is 'Both', please specify FCCutoff as either 'NULL' or 'c(x,y)' where 'x' is the cutoff for downregulated genes and 'y' the cutoff for upregulated genes")
      }
    }
  }
  
  # Sanitize input groupNames
  # TODO: Actually sanitize this input.
  
  print("Preparing gene associations")
  if(config$species == 'mouse') {
    organismToStudy = 'Mus musculus'
    if(GOKegg == 'GO') {
      if(includeChildTerms) {
        Terms = toTable(org.Mm.egGO2ALLEGS)
        includeChildTerms = 'WithChildren'
      } else {
        Terms = toTable(org.Mm.egGO2EG)
        includeChildTerms = 'WithoutChildren'
      }
      termFrameData = data.frame(Terms$go_id, Terms$Evidence, Terms$gene_id)
      termFrame = GOFrame(termFrameData,organism=organismToStudy)
      termFrame = GOAllFrame(termFrame)
      termGSC = GeneSetCollection(termFrame, setType = GOCollection())
    } else {
      Terms = toTable(org.Mm.egPATH)
      termFrameData = data.frame(Terms$path_id, Terms$gene_id)
      termFrame = KEGGFrame(termFrameData,organism=organismToStudy)
      termGSC = GeneSetCollection(termFrame, setType = KEGGCollection())
    }
  } else if(config$species == 'human') {
    organismToStudy = 'Homo sapiens'
    if(GOKegg == 'GO') {
      if(includeChildTerms) {
        Terms = toTable(org.Hs.egGO2ALLEGS)
        includeChildTerms = 'WithChildren'
      } else {
        Terms = toTable(org.Hs.egGO2EG)
        includeChildTerms = 'WithoutChildren'
      }
      termFrameData = data.frame(Terms$go_id, Terms$Evidence, Terms$gene_id)
      termFrame = GOFrame(termFrameData,organism=organismToStudy)
      termFrame = GOAllFrame(termFrame)
      termGSC = GeneSetCollection(termFrame, setType = GOCollection())
    } else {
      Terms = toTable(org.Hs.egPATH)
      termFrameData = data.frame(Terms$path_id, Terms$gene_id)
      termFrame = KEGGFrame(termFrameData,organism=organismToStudy)
      termGSC = GeneSetCollection(termFrame, setType = KEGGCollection())
    }
  }
  
  # Get external_gene_name for all ensembl_gene_id (where available)
  if(exists('geneIdentifierConversionTable')) {
    # Use previously defined geneIdentifierConversionTable
    print('Loaded existing geneIdentifierConversionTable. Use remove(\'geneIdentifierConversionTable\') to regenerate this table.')
    geneNames <- geneIdentifierConversionTable
  } else {
    # Define geneIdentifierConversionTable
    print('Generating geneIdentifierConversionTable.')

    if(config$geneIdentifierType %in% c('external_gene_name','external_gene_name_with_chr')) {
      # Start dataframe with data that's available
      geneNames = data.frame(
        'ensembl_gene_id' = NA,
        'external_gene_name' = if(config$geneIdentifierType == 'external_gene_name_with_chr') gsub("__chr(\\d+|[MXY])", '', rownames(groupedSCS[[1]]@ndata)) else rownames(groupedSCS[[1]]@ndata),
        'entrezgene' = NA,
        stringsAsFactors = T
      )
      # Then add entrez genes
      geneNames = merge(geneNames, toTable(org.Hs.egSYMBOL), by.x='external_gene_name', by.y='symbol', all.x=T)
      geneNames$entrezgene = geneNames$gene_id
      # Then add ensembl gene ids
      geneNames = merge(geneNames, toTable(org.Hs.egENSEMBL), by.x='entrezgene', by.y='gene_id')
      geneNames$ensembl_gene_id = geneNames$ensembl_id
      # Throw away extra columns
      geneNames = geneNames[,c('ensembl_gene_id', 'external_gene_name','entrezgene')]
      
    } else if(config$geneIdentifierType == 'entrezgene') {
      # Start dataframe with available data
      geneNames = data.frame(
        'ensembl_gene_id' = NA,
        'external_gene_name' = NA,
        'entrezgene' = rownames(groupedSCS[[1]]@ndata),
        stringsAsFactors = T
      )
      # Then add ensembl gene ids
      geneNames = merge(geneNames, toTable(org.Hs.egENSEMBL), by.x='entrezgene', by.y='gene_id')
      geneNames$ensembl_gene_id = geneNames$ensembl_id
      # Then add gene symbols
      geneNames = merge(geneNames, toTable(org.Hs.egSYMBOL), by.x='entrezgene', by.y='gene_id')
      geneNames$external_gene_name = geneNames$symbol
      # Throw away extra columns
      geneNames = geneNames[,c('ensembl_gene_id', 'external_gene_name','entrezgene')]
      
    } else if(config$geneIdentifierType == 'ensembl_gene_id') {
      # Start dataframe with available data
      geneNames = data.frame(
        'ensembl_gene_id' = rownames(groupedSCS[[1]]@ndata),
        'external_gene_name' = NA,
        'entrezgene' = NA,
        stringsAsFactors = T
      )
      # Then add entrez genes
      geneNames = merge(geneNames, toTable(org.Hs.egENSEMBL), by.x='ensembl_gene_id', by.y='gene_id', all.x=T)
      geneNames$entrezgene = geneNames$gene_id
      # Then add gene symbols
      geneNames = merge(geneNames, toTable(org.Hs.egSYMBOL), by.x='entrezgene', by.y='gene_id')
      geneNames$external_gene_name = geneNames$symbol
      # Throw away extra columns
      geneNames = geneNames[,c('ensembl_gene_id', 'external_gene_name','entrezgene')]
      
    }
    # Save geneIdentifierConversionTable so it doesn't need to be generated for the next group (and is available outside this function)
    assign('geneIdentifierConversionTable', geneNames, .GlobalEnv)
  }
  
  # Prepare return object
  returnData = list()
  
  # Iterate over the available groups
  for(groupName in correlationResult$groupName) {
    # Determine the background
    meanExpression = rowMeans(groupedSCS[[groupName]]@ndata)
    backgroundGenes = names(meanExpression)[which(meanExpression > minExpression)]
    if(config$geneIdentifierType == 'external_gene_name_with_chr') {
      backgroundGenes = geneNames$entrezgene[which(geneNames$external_gene_name %in% gsub("__chr(\\d+|[MXY])", '', backgroundGenes))]
    } else {
      backgroundGenes = geneNames$entrezgene[which(geneNames[[config$geneIdentifierType]] %in% backgroundGenes)]
    }

    filteredGenes <- correlationResult$correlationData
    if(posNeg == 'Pos') {
      filteredGenes <- filteredGenes %>% dplyr::filter(correlation >= corrCutoff)
    } else if(posNeg == 'Neg') {
      filteredGenes <- filteredGenes %>% dplyr::filter(correlation <= corrCutoff)
    } else if(posNeg == 'Both') {
      filteredGenes <- filteredGenes %>% dplyr::filter(correlation >= corrCutoff[1] | correlation <= 1/corrCutoff[2])
    } else if(is.null(posNeg)) {
      #No filtering on fold change
      filteredGenes <- correlationResult$correlationData
    }
      
    # Further filter genes for pValue (if required) 
    if(!is.null(pCutoff)) {
      filteredGenes <- filteredGenes %>% dplyr::filter(pValue <= pCutoff)
    }
    
    # Just continue with the gene names
    filteredGenes <- filteredGenes$genes
    
    # Get entrez gene IDs
    if(config$geneIdentifierType == 'external_gene_name_with_chr') {
      filteredGenes = geneNames$entrezgene[which(geneNames$external_gene_name %in% gsub("__chr(\\d+|[MXY])", '', filteredGenes))]
    } else {
      filteredGenes = geneNames$entrezgene[which(geneNames[[config$geneIdentifierType]] %in% filteredGenes)]
    }
  
    # Filter out genes that are not in the background (min expression filtering is not applied to DE-list)
    filteredGenes <- filteredGenes[filteredGenes %in% backgroundGenes]
    
    # Perform gene ontology analysis. 
    if(length(filteredGenes) > 0) {
      # Set up analysis
      print(paste0("Analyzing '",groupName))
      if(GOKegg == 'GO') {
        params <- GSEAGOHyperGParams(
          name="GO GSEA",
          geneSetCollection=termGSC,
          geneIds = filteredGenes,
          universeGeneIds = backgroundGenes,
          ontology = 'BP',
          pvalueCutoff = pathwayPCutoff,
          conditional = T,
          testDirection = 'over'
        )          
      } else {
        params <- GSEAKEGGHyperGParams(
          name="KEGG GSEA",
          geneSetCollection=termGSC,
          geneIds = filteredGenes,
          universeGeneIds = backgroundGenes,
          pvalueCutoff = pathwayPCutoff,
          testDirection = 'over'
        )          
      }
      
      # Run enrichment analysis
      hyperGResult <- hyperGTest(params)
      termResult <- summary(hyperGResult)
      
      # Process results
      if(nrow(termResult) > 0) {
        # Prepare collection table for the regulated genes per GO Term
        if(GOKegg == 'GO') {
          genesInTerm = data.frame(
            'TermID' = termResult$GOBPID,
            'genes' = '',
            stringsAsFactors = F
          )
        } else {
          genesInTerm = data.frame(
            'TermID' = termResult$KEGGID,
            'genes' = '',
            stringsAsFactors = F
          )
        }
        
        # Rename ID columns in termFrame for GO/KEGG so they can be parsed the same way
        if(GOKegg == 'GO') {
          colnames(termFrame@data) <- c('term_id', 'evidence', 'gene_id')
        } else if(GOKegg == 'KEGG') {
          colnames(termFrame@data) <- c('term_id', 'gene_id')
        }
        
        # For each GO Term, find all regulated genes association with it, store them in the previously made collection table.
        for(i in 1:nrow(termResult)) {
          #GOid to look up
          termID = termResult[i,1]
          
          #Genes in this GO-Term (each gene can be represented multiple times due to multiple evidence levels, hence the call to unique() )
          termGenes = unique(termFrame@data[which((termFrame@data$term_id == termID) & (termFrame@data$gene_id %in% backgroundGenes)), 'gene_id'])
          #Regulated genes in this GO-Term
          regulatedTermGenes = termGenes[which(termGenes %in% filteredGenes)]
          #Convert entrez IDs back to external_gene_names
          regulatedTermGeneNames = geneNames %>% dplyr::filter(entrezgene %in% regulatedTermGenes)
          regulatedTermGeneNames = unique(regulatedTermGeneNames$external_gene_name)
          #Add names of regulated genes to genesInGo table
          genesInTerm[i,] = c(termID, paste0(regulatedTermGeneNames, collapse=';'))
        } 
        
        # Merge list of regulated genes per GO Term with result of GO analysis
        if(GOKegg == 'GO'){
          termResult <- merge(termResult, genesInTerm, by.x='GOBPID', by.y='TermID')
        } else {
          termResult <- merge(termResult, genesInTerm, by.x='KEGGID', by.y='TermID')
        }
        
        # Order goResult by pValue
        termResult <- termResult[order(termResult$Pvalue),]
      } else {
        termResult$genes = character()
      } 
    
      #Add analysis to returnData variable
      returnData[[groupName]][[correlationResult$geneName]] = termResult
    } else {
      print(paste0("No significantly regulated genes for '",groupName,"': ",correlationResult$geneName,"."))
    }
  }
    
  #Return analysis data
  return(returnData)
}