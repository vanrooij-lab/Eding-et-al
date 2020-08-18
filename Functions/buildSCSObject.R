
### buildSCSObject(config, groupedData, groupNames=NULL)
### config: Pass the full config variable from the analysis pipeline, this function automatically extracts the correct parameters.
### groupedData: Pass the full groupedData-variable obtained from loadData.
### groupNames: Optional. Pass the name of a single group to perform read distribution analysis for a single group, or a list of 
###             names to analyze multiple groups.. By default (=NULL) builds for all.
###             
### Return: list containing an SCSObject per group
### 
### Authors
###   - Joep Eding
### 
### TODO
buildSCSObject <- function(config, groupedData, groupNames=NULL, groupedSCS= list()) {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedData)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedData))) {stop("You want to build SCSObject for group '",groupName,"' that is not in groupedData")}
    }
  }
  
  # The return variable "groupedSCS" is initialized in the 
  # function declaration; allowing you to also provide an 
  # already existing "groupedSCS" object to the function
  # to which new analyses can be added.
  
  for(groupName in groupNames) {
    #Make SingleCellSequencing object
    print(paste0("Making SCSObject for ",groupName))
    groupedSCS[[groupName]] <- SCseq(groupedData[[groupName]])
    
    #Filter expression data
    print("---Filtering expression data")
    groupedSCS[[groupName]] <- filterdata(
      groupedSCS[[groupName]], 
      mintotal = config$minTotalReads, 
      minexpr = 3, 
      minnumber = 1, 
      maxexpr = Inf, 
      downsample = T, 
      dsn = 1, 
      rseed = 17000
    )
    
    #Do k-medoids clustering
    print("---Clustering")
    groupedSCS[[groupName]] <- clustexp(
      groupedSCS[[groupName]],
      clustnr=30,
      bootnr=50,
      metric="pearson",
      do.gap=FALSE,
      sat=TRUE,
      SE.method="Tibs2001SEmax",
      SE.factor=.25,
      B.gap=50,
      cln=0,
      rseed=17000,
      FUNcluster="kmedoids"
    )
    
    #Compute t-SNE map
    print("---Compute tSNE-map")
    groupedSCS[[groupName]] <- comptsne(groupedSCS[[groupName]], rseed=15555)    
    
    #Detect outliers and redefine cluster
    print("---Detect outliers")
    groupedSCS[[groupName]] <- findoutliers(
      groupedSCS[[groupName]], 
      outminc=5,
      outlg=4,
      probthr=config$outlierPValue,
      thr=2**-(1:40),
      outdistquant=.95
    )
  }
  
  #Return SCSObjects
  return(groupedSCS)
}