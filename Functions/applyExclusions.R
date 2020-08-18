### applyExclusions(config, groupedData)
### 
### Parameters
###   - config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
###   - groupedData: Pass the full groupedData-variable obtained from loadData.
###
### Return: 
###   - groupedData with the relevant exclusions applied.
### 
### Authors
###   - 
###   
### TODO
###   - 
applyExclusions <- function(config, groupedData) {
  #Apply exclusions
  for(groupName in names(groupedData)) {
    groupedData[[groupName]] <- groupedData[[groupName]][(
      #Exclude spike-ins
      !grepl('ERCC', row.names(groupedData[[groupName]]))
      &
      #Exclude mitochondrial genes if required
      (
        !grepl('__chrM',row.names(groupedData[[groupName]]))
        |
        !config$excludeMitochondrialGenes
      ) &
      # Now add custom exclusions if applicable
      (
          if (is.null(config$custom_gene_blacklist)) { 
              T 
          } else {
                !grepl(paste(config$custom_gene_blacklist,collapse="|"), 
                        row.names(groupedData[[groupName]]))
          }
      )
    ),]
  }
  
  #return data after exclusions were applied
  return(groupedData)
}
