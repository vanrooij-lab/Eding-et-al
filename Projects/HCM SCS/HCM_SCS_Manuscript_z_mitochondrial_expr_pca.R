

# Script that looks into cell characterization using only mitochondrial reads.
# It's interesting to see that the patients separate pretty well, so 
# this is another indication that maybe more sophisticated methods are needed
# to correct for batch effects.

# Note: this is generated in the script
# HCM_SCS_Manuscript_d_mitochondrial_and_metaparameters_correlations.R
mitoCountTable_filt

# Get some info
cellsnames_patientAllMod = names(groupedSCS$patientAllMod@cluster$kpart)
clustering_patientAllMod = as.factor(groupedSCS$patientAllMod@cluster$kpart)
patient_source_ptAllMod=shorthand_give_p_source_forcellnames(config = config,cellnames = cellsnames_patientAllMod)

# Perform PCA
pca_out = prcomp(t(scale(mitoCountTable_filt[,cellsnames_patientAllMod]))) # note: selection of correct cells is redundant
summary(pca_out)

# Let's see if we can see some pattern in the PCA (not really)
ggplot(data.frame(pca1=pca_out$x[,1],pca2=pca_out$x[,2]))+
    geom_point(aes(x=pca1, y=pca2))+theme_bw()
# Let's see if patient source shows pattern by mitochondrial expression (it does)
ggplot(data.frame(pca1=pca_out$x[,1],pca2=pca_out$x[,2], patient=patient_source_ptAllMod))+
    geom_point(aes(x=pca1, y=pca2, color=patient))+theme_bw()+theme_bw()
# pca3, pca4 (patient coloring)
ggplot(data.frame(pca3=pca_out$x[,3],pca4=pca_out$x[,4],
                patient = patient_source_ptAllMod))+
    geom_point(aes(x=pca3, y=pca4, color=patient))+theme_bw()

# See if the clustering from Race is related to mitochondrial gene expression
# pca1, pca2
ggplot(data.frame(pca1=pca_out$x[,1],pca2=pca_out$x[,2],
                cluster = clustering_patientAllMod))+
    geom_point(aes(x=pca1, y=pca2, color=cluster))+theme_bw()
# pca3, pca4
ggplot(data.frame(pca3=pca_out$x[,3],pca4=pca_out$x[,4],
                cluster = clustering_patientAllMod))+
    geom_point(aes(x=pca3, y=pca4, color=cluster))+theme_bw()


# Also use package to check out relation original genes to pcs

#install.packages('devtools')
#library(devtools)
#install_github("vqv/ggbiplot")

library(ggbiplot)

ggbiplot(pca_out)+theme_bw()
ggbiplot(pca_out, ellipse = T, alpha=0.2)+theme_bw()

################################################################################

# While we're add it, also perform tsne
# Perform tsne
rtsne_out = Rtsne::Rtsne(t(scale(mitoCountTable_filt[,cellsnames_patientAllMod])))
summary(rtsne_out)

# Show result
ggplot(data.frame(tsne1=rtsne_out$Y[,1],tsne2=rtsne_out$Y[,2]))+
    geom_point(aes(x=tsne1, y=tsne2))+theme_bw()
# Clusters
ggplot(data.frame(tsne1=rtsne_out$Y[,1],tsne2=rtsne_out$Y[,2],
                cluster = clustering_patientAllMod))+
    geom_point(aes(x=tsne1, y=tsne2, color=cluster))+theme_bw()
# Patients
ggplot(data.frame(tsne1=rtsne_out$Y[,1],tsne2=rtsne_out$Y[,2],
                patient = patient_source_ptAllMod))+
    geom_point(aes(x=tsne1, y=tsne2, color=patient))+theme_bw()

