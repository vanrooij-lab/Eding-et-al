

# THIS ENTIRE SCRIPT IS OLD AND CAN BE REMOVED


# First redo some stuff
regulons_redone <- analyzeRegulons(config, correlationMatrices, minCorrelations=40, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig4/A/extra', outputMode='pdf')
regulons_selectedGenes_redone <- analyzeRegulons(config, correlationMatrices_selectedGenes, minCorrelations=10, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig4/A', outputMode='pdf', K.max = 50, filename_precursor='genes_allPs')

# Now make some little plots
length(regulons_selectedGenes_redone$patient1Mod$regulons$regulon.1)
length(regulons_selectedGenes$patient1Mod$regulons$regulon.1)


#########





# Old stuff
"

ggplot(melt(silhouette_tbl[[1]]))+
    geom_bar(aes(x=Var1,y=value,fill=Var2),position = 'dodge',stat='identity')

ggplot(melt(silhouette_tbl[[1]]))+
    geom_boxplot(aes(x=Var1,y=value))+
    coord_flip()

ggplot(melt(silhouette_tbl[[1]]))+
    geom_boxplot(aes(x=Var1,y=value))+
    coord_flip()
    
# same but sorted by median
ggplot(melt(silhouette_tbl_md[[1]]))+
    geom_jitter(aes(x=Var1,y=value, color=Var1))+
    coord_flip()+
    scale_color_manual(values=rep(sample(col_Dark2),3))+ #col_vector_60
    theme_bw()+theme(legend.position = 'none')



regulons_selectedGenes$patient2Mod$silhouetteData[current_genenames_chrXX,]$sil_width


current_geneClustering = regulons$patient1Mod$clustering

clusteringSilhouette <- silhouette(current_geneClustering, as.dist(1-subCorrelationMatrix))
    silhouetteData <- as.data.frame(clusteringSilhouette[,])
    rownames(silhouetteData) <- names(current_geneClustering)
"
