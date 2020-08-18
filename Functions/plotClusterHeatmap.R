### plotClusterHeatmap(config, groupedSCS, groupNames=NULL, outputMode='pdf') 
### 
### Authors
###   - Joep Eding
### 
### TODO Make this a useful comment
plotClusterHeatmap <- function(config, groupedSCS, groupNames=NULL, includeOutliers=T, hmethod='single', overrideDir=F, outputMode='pdf') {
  #Sanitize input groupNames
  if(is.null(groupNames)) {
    groupNames = names(groupedSCS)
  } else {
    for(groupName in groupNames) {
      if(!(groupName %in% names(groupedSCS))) {stop(paste0("You want to draw gene expression t-SNE map(s) for SCSObject '",groupName,"' that is not in groupedSCS"))}
    }
  }
  
  #Set saveGraphs based on outputMode
  if(!(outputMode %in% c('show','png','pdf','eps','ps','tiff','svg'))) {
    stop("'outputMode' needs to be one of 'show','png','pdf','eps','ps','tiff','svg'. Pay attention to capitalisation.")
  }
  
  for(groupName in groupNames) {
    #Create directory if necessary
    if(outputMode!='show') {
      if(overrideDir == FALSE) {
        saveDir <- paste0(outputDir,.Platform$file.sep,groupName,.Platform$file.sep,'Quality control',.Platform$file.sep,'processed')
      } else {
        saveDir <- paste0(outputDir,.Platform$file.sep,overrideDir,.Platform$file.sep,groupName)
      }
      if(!dir.exists(saveDir)) {dir.create(saveDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")}
    }
    
    #Get the clustering information
    clustering <- if(includeOutliers) groupedSCS[[groupName]]@cpart else groupedSCS[[groupName]]$kpart
    
    #Calculate mean expression level of every gene per cluster
    na <- c() #List for which clusters were included
    j <- 0 #Counter for number of clusters included  (Only used once? Could be replaced by length(na)?)
    for(clusterNum in 1:max(clustering)) {
      if(sum(clustering == clusterNum) == 0) { next }
      j=j+1
      na <- append(na,clusterNum)
      clusterCells = groupedSCS[[groupName]]@fdata[,clustering==clusterNum]
      
      #Calculate row means ('center')
      if(sum(clustering==clusterNum) == 1) { 
        center <- clusterCells 
      } else { 
        center = apply(clusterCells,1,mean) #Shouldn't this just be rowMeans for clarity? There's probably a speed difference.
      }
      
      #Put new cluster as column into tmp dataframe
      if(j==1) {
        tmp <- data.frame(center)
      } else {
        tmp <- cbind(tmp,center)
      }
    }
    # Make appropriate column names
    names(tmp) <- paste('cl',na,sep='.')
    
    #Calculate (local???) distances?
    ld <- if(groupedSCS[[groupName]]@clusterpar$FUNcluster == 'kmedoids') dist.gen(t(tmp), method=groupedSCS[[groupName]]@clusterpar$metric) else dist.gen(as.matrix(dist.gen(t(tmp),method=groupedSCS[[groupName]]@clusterpar$metric)))
    if(max(clustering) > 1) clusterOrdering <- hclust(ld, method=hmethod)$order else clusterOrdering <- 1

    # ???    
    originalClustering <- clustering
    tmpClustering <- clustering
    for(clusterNum in 1:max(clustering)) {
      tmpClustering[clustering == na[clusterOrdering[clusterNum]]] <- clusterNum
    }
    clustering <- tmpClustering
         
    #Calculate distances?
    distances <- if(groupedSCS[[groupName]]@clusterpar$FUNcluster == 'kmedoids') groupedSCS[[groupName]]@distances else as.data.frame(as.matrix(dist.gen(t(groupedSCS[[groupName]]@distances))))
    pto <- clustering[order(clustering, decreasing=FALSE)]
    ptn <- c()
    for(i in 1:max(pto)) {
      pt <- names(pto)[pto==i]; 
      z <- if(length(pt)==1) pt else pt[hclust(as.dist(t(distances[pt,pt])),method=hmethod)$order]; 
      ptn <- append(ptn,z) 
    }
    
    #Calculate cluster centers
    clusterCenters = list()
    for(clusterNum in 1:max(clustering)) {
      #Get all cells in cluster
      clusterCellNames = names(originalClustering[which(originalClustering==clusterNum)])
      
      #Get ordered cells in cluster
      orderedClusterCellNames = ptn[which(ptn %in% clusterCellNames)]
      
      #Get cluster center
      newCenter = orderedClusterCellNames[ceiling(length(orderedClusterCellNames)/2)] #Use 'ceiling' to get an integer and always have at least '1'
      
      clusterCenters <- append(clusterCenters, newCenter)
    }
    clusterCenters = unlist(clusterCenters)
    
    ptnNumbers = seq(0,1,length.out = length(ptn))
    names(ptnNumbers) <- ptn
    
    heatmapData = melt(as.matrix(distances[ptn,ptn]))
    heatmapData$Var1 = unlist(lapply(heatmapData$Var1, function(x) {as.double(ptnNumbers[x])}))
    heatmapData$Var2 = unlist(lapply(heatmapData$Var2, function(x) {as.double(ptnNumbers[x])}))
    
    colorRampNames = seq(0,max(heatmapData$value+0.001),0.001)
    heatmapData$value = round(heatmapData$value,3)
    
    colorRampResult <- colorRampPalette(brewer.pal(7, 'RdYlBu'))(length(colorRampNames))
    names(colorRampResult) <- colorRampNames
    heatmapData$xLabel = NA
    heatmapData$yLabel = NA
    
    #Add labels for x-axis
    # TODO Make this faster by making a call to the right row straigh away instead of looping over all.
    for(clusterNum in 1:length(clusterCenters)) {
      #print(heatmapData[which((heatmapData$Var1 == clusterCenter) & (heatmapData$Var2 == ptn[1])),])
      heatmapData[which((heatmapData$Var1 == ptnNumbers[clusterCenters[clusterNum]]) & (heatmapData$Var2 == ptnNumbers[ptn[1]])),'xLabel'] = clusterNum
    }
    #Add labels for y-axis
    # TODO Make this faster by making a call to the right row straigh away instead of looping over all.
    for(clusterNum in 1:length(clusterCenters)) {
      heatmapData[which((heatmapData$Var2 == ptnNumbers[clusterCenters[clusterNum]]) & (heatmapData$Var1 == ptnNumbers[ptn[1]])),'yLabel'] = clusterNum
    }
    heatmapData2 <- heatmapData[1:(200*1809),]
    heatmapDataBackup <- heatmapData
    
    #Plot heatmap as geom_points (this is kinda hacky, but essentially the same as the native RaceID2 functions does it)
    # TODO Consider putting cluster labels on axes as major ticks instead of using geom_text_repel
    heatmap <- ggplot(
      data=heatmapData,
      aes(x=as.double(Var1),y=as.double(Var2),color=as.character(value))
    ) + geom_point(
      shape=15, # Using circles for points results in weird curves at the upper edge of the heatmap.
      size=0.01 # Use very small points, otherwise the top row and right column appear too large
    ) + scale_color_manual(
      values = colorRampResult,
      guide=FALSE
    ) + theme(
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = "left"
    ) + labs(
      title='Cluster heatmap'
    ) + scale_y_continuous(
      limits = c(-0.2,1.0),
      expand = c(0,0)
    ) + scale_x_continuous(
      limits = c(-0.2,1.0),
      expand = c(0,0)
    ) + geom_text_repel(
      aes(label=xLabel),
      ylim = c(-0.02,-0.18),
      xlim = c(-0,1),
      min.segment.length = grid::unit(0,'pt'),
      color='grey30',
      size=5/.pt,
      hjust=0.5,
      box.padding = 0.15,
      segment.size = 0.01,
      force = 1,
      angle=if(max(clustering)>9) 270 else 0, #Tilt x-axis labels only when numbers go above 9.
      na.rm = T #Makes removal of missing values silent
    ) + geom_text_repel(
      aes(label=yLabel),
      xlim = c(-0.02,-0.09),
      ylim = c(0,1),
      min.segment.length = grid::unit(0,'pt'),
      color='grey30',
      size=5/.pt,
      hjust=0.5,
      box.padding = 0.15,
      segment.size = 0.01,
      force = 1,
      na.rm = T #Makes removal of missing values silent
    )
    # TODO Add cluster colors along the axes.

    #Save or print
    if(outputMode != 'show') {
      ggsave(paste0(groupName,'_clusterHeatmap.',outputMode),path = saveDir)
    } else {
      print(heatmap)
    }
  }
}