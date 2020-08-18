 for(groupName in groupNames) {
    #Get the clustering information
    clustering <- if(includeOutliers) groupedSCS[[groupName]]@cpart else groupedSCS[[groupName]]$kpart
    
    #Calculate mean expression level of every gene per cluster
    na <- c() #List for which clusters were included
    j <- 0 #Counter for number of clusters included  (Only used once? Could be replaced by length(na)?)
    for(clusterNum in 1:max(clustering)) {
      if(sum(clustering == clusterNum) == 0) { next }
      j=j+1
      na <- append(na,clusterNum)
      clusterCells = groupedSCS[[groupName]]@ndata[,clustering==clusterNum]
      
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
    names(tmp) <- paste(groupName,'cl',na,sep='')
    assign(paste0(groupName,'Tmp'),tmp)
 }

#Merge pt1 & 2
testTmp <- merge(patient1Tmp,patient2Tmp,by=0)
rownames(testTmp) <- testTmp$Row.names
testTmp <- testTmp[2:ncol(testTmp)]
#Merge in pt4
testTmp <- merge(testTmp,patient4Tmp,by=0)
rownames(testTmp) <- testTmp$Row.names
testTmp <- testTmp[2:ncol(testTmp)]
#Merge in ptAll
testTmp <- merge(testTmp,patientAllTmp,by=0)
rownames(testTmp) <- testTmp$Row.names
testTmp <- testTmp[2:ncol(testTmp)]

#Prepare for DESeq
testTmp <- floor((testTmp-0.1)*10)
metadataThing <- data.frame(
  cluster=colnames(testTmp),
  sample = unlist(lapply(colnames(testTmp),function(x){str_split(x,'cl')[[1]][1]}))
)
rownames(metadataThing) <- metadataThing$cluster

#Uncomment line below to exclude clusters
#metadataThing <- metadataThing[!grepl('1cl9|2cl11',rownames(metadataThing2)),]

#Make deseq dataset
dds <- DESeqDataSetFromMatrix(
  countData=testTmp[,rownames(metadataThing)], 
  colData=metadataThing, 
  design= ~ sample
)
dds2 <- DESeq(dds)
rld <- rlog(dds2, blind=FALSE)
pcadata <- plotPCA(rld, intgroup=colnames(metadataThing), returnData=T)
print(
  ggplot(
    pcadata, 
    aes(PC1, PC2)
  ) + geom_point(
    aes(fill=sample), 
    color='black', 
    shape=21, 
    size=4, 
    alpha=.9
  ) + geom_text_repel(
    aes(label=cluster), 
    vjust=3, 
    size=3
  )
)
