TRIAGEData <- read.table('/Users/j.eding/Work/Projects/Nathan Palpant/Single-Cell in HCM/discordance.csv', header=T, row.names=1, sep=";", dec=',')

oldCluster = list()
oldClusterTriage = createWorkbook()
for(clusterNum in 1:max(groupedSCS$patientAll@cpart)) {
  #Copy cluster diff exp
  oldCluster[[paste0('cl.',clusterNum)]] = as.data.frame(clusterDifferentialExpression$patientAll[[paste0('cl.',clusterNum)]])
  
  #Add a column containing the gene name without chromosome information
  oldCluster[[paste0('cl.',clusterNum)]]$colForMerge = as.character(lapply(rownames(oldCluster[[paste0('cl.',clusterNum)]]), function(x) {strsplit(x, '__')[[1]][1]}))
  oldCluster[[paste0('cl.',clusterNum)]]$rowNames = rownames(oldCluster[[paste0('cl.',clusterNum)]])
  
  #Get cells from cluster
  clusterCells = names(groupedSCS$patientAll@cpart[groupedSCS$patientAll@cpart == clusterNum])
  TRIAGESubset = TRIAGEData[,clusterCells]
  if(ncol(as.data.frame(TRIAGESubset)) == 1) {
    TRIAGEMEANS = TRIAGESubset
  } else {
    TRIAGEMeans = rowMeans(TRIAGESubset, na.rm=F)
  }
  #TRIAGEMEANS[which(is.na(TRIAGEMEANS))] = 0
  TRIAGEClusterMean = data.frame(
    colForMerge = as.character(rownames(TRIAGEData)),
    TRIAGE = TRIAGEMeans
  )
  
  #Merge TRIAGE data into diff exp
  oldCluster[[paste0('cl.',clusterNum)]] = merge(as.data.frame(oldCluster[[paste0('cl.',clusterNum)]]), TRIAGEClusterMean, by.x='colForMerge', all=T)
  
  #Write to xlsx-file
  addWorksheet(oldClusterTriage, sheetName=paste0('cl.',clusterNum))
  writeDataTable(
    oldClusterTriage,
    sheet = paste0('cl.',clusterNum),
    oldCluster[[paste0('cl.',clusterNum)]]
  )
}
saveWorkbook(oldClusterTriage, paste0(clusterDir,.Platform$file.sep,'clustDiffExpTRIAGE.xlsx'))

geneNames = c('TTN','MYH7','MYBPC3','KCNQ1OT1','NPPA','MTRNR2L1')
for(geneName in geneNames) {
  ggplot(
    data=data.frame(
      cell = rownames(TRIAGEData[geneName,]),
      discrepancy = as.numeric(TRIAGEData[geneName,]),
      cluster = as.character(groupedSCS$patientAll@cluster$kpart),
      stringsAsFactors = F
    ),
    aes(
      x=discrepancy
    )
  ) + geom_density(
    aes(
      group=cluster,
      color=cluster
    )
  ) + geom_density(
    aes(
    )
  ) + labs(
    title=geneName,
    x = paste0(geneName,' discrepancy')
  )
  ggsave(paste0('discr_',geneName,'.pdf'),path='/Users/j.eding/Work/Projects/Nathan Palpant/Single-Cell in HCM/')
  
  ggplot(
    data=data.frame(
      cell = rownames(TRIAGEData[geneName,]),
      expression = as.numeric(groupedSCS$patientAll@ndata[grepl(paste0('^',geneName,'_'),rownames(groupedSCS$patientAll@ndata)),]),
      cluster = as.character(groupedSCS$patientAll@cluster$kpart),
      stringsAsFactors = F
    ),
    aes(
      x=expression
    )
  ) + geom_density(
    aes(
      group=cluster,
      color=cluster
    )
  ) + geom_density(
    aes(
    )
  ) + labs(
    title=geneName,
    x = paste0(geneName,' expression')
  )
  ggsave(paste0('expr_',geneName,'.pdf'),path='/Users/j.eding/Work/Projects/Nathan Palpant/Single-Cell in HCM/')
}



TRIAGEData2 = list()
TRIAGEData2$TRIAGE = TRIAGEData
for(colN in 1:ncol(TRIAGEData2$TRIAGE)) {
  TRIAGEData2$TRIAGE[which(is.na(TRIAGEData2$TRIAGE[,colN])), colN] = 0  
}
TRIAGEData2$TRIAGE = floor(TRIAGEData2$TRIAGE * (1000))
TRIAGEconfig = config
#TRIAGEconfig$minTotalReads = 1
TRIAGEconfig$outlierPValue = 1e-13
TRIAGESCS = buildSCSObject(TRIAGEconfig, TRIAGEData2, 'TRIAGE')

#Look at data
processedDataQualityControl(TRIAGEconfig, TRIAGESCS, groupNames = NULL, outputMode = 'save')
analyseClusterComposition(TRIAGEconfig, TRIAGESCS, groupNames = NULL, outputMode = 'save')
analyseClusterDistribution(TRIAGEconfig, TRIAGESCS, groupNames = NULL, outputMode = 'save')
TRIAGEClusterDifferentialExpression <- analyseClusterDifferentialExpression(TRIAGEconfig, TRIAGESCS, groupNames = NULL, outputMode = 'save')
plotClusterTSNE(config, TRIAGESCS, groupNames = 'TRIAGE', includeOutlierClusters = F, includeOutlierCells = T, outputMode='pdf')
plotGeneExpressionTSNE(config, TRIAGESCS, geneNames = 'TTN',groupNames = 'TRIAGE', colorScheme='Spectral', expressionScale = 'linear', outputMode = 'pdf')
