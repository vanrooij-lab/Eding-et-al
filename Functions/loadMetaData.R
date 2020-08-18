### loadMetaData(config, groupedSCS, metadata)
### 
### Parameters
###   - config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
###   - groupedSCS: 
###   - metadata: 
### 
### Return
###   - 
###   
### Authors
###   - 
### 
### TODO:
###  - Add mechanism to prevent double column names. Currently double column names result in addition of suffixes which may break analysis.
loadMetaData <- function(config, groupedSCS, metaData) {
  # Check whether a metadata variable already exists
  if(exists('scsMetadata')) {
    # metadata variable already exists:
  } else {
    # metadata variable does not exist yet, create it:
    print("Creating metadata variable 'scsMetadata'")
    assign(
      'scsMetadata', 
      data.frame(
        stringsAsFactors = F,
        Cell = character()
      ),
      envir = .GlobalEnv
    )
  }
  
  # Check whether metadata refers to existing cells
  allCellNames = list()
  for(groupName in names(groupedSCS)) {
    allCellNames <- unlist(append(allCellNames, colnames(groupedSCS[[groupName]]@ndata)))
  }
  if(!all(metaData$Cell %in% allCellNames)) {
    print(paste0("Ignoring metadata for cells that don't exist in any group: '",paste(metaData$Cell[which(!(metaData$Cell %in% allCellNames))], collapse = "', '"),"'."))
    metaData = metaData[which(metaData$Cell %in% allCellNames),]
  }
  
  # Add metadata to metadata variable
  assign(
    'scsMetadata',
    merge(scsMetadata, metaData, by='Cell', all=T),
    envir = .GlobalEnv
  )
}