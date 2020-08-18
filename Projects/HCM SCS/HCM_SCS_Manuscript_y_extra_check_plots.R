

################################################################################
# This script tries several pre-processing choices, and compares
# their results.
# Importantly, it also includes checks to see whether the patterns 
# identified from the data don't rely too heavily on highly-expressed
# genes, specifically TTN.
################################################################################


# Some repeated analysis as Joep did them, with some variations
# First run "HCM SCS Manuscript_AE.R" up to "groupedDataBackup <- groupedData"
# -----

# Apply exclusions
# Only default exclusions like mito genes
config$custom_gene_blacklist=NULL
groupedData_all      <- applyExclusions(config, groupedDataBackup)
# Also excluse KCNQ1OT & MALAT1 ("no K & M")
config$custom_gene_blacklist=c("^MALAT1_","^KCNQ1OT1_")
groupedData_noKM     <- applyExclusions(config, groupedDataBackup)
# Also excluse KCNQ1OT, MALAT1, TTN ("no K, M & T")
config$custom_gene_blacklist=c("^MALAT1_","^KCNQ1OT1_","^TTN_")
groupedData_noKMT    <- applyExclusions(config, groupedDataBackup)

# -----
# Build SCSObjects for the groupedData
groupedSCS_all   <- buildSCSObject(config, groupedData_all,   groupNames = 'patientAll')
groupedSCS_noKM  <- buildSCSObject(config, groupedData_noKM,  groupNames = 'patientAll')
groupedSCS_noKMT <- buildSCSObject(config, groupedData_noKMT, groupNames = 'patientAll')


# -----
clusterDifferentialExpression_all <- analyseClusterDifferentialExpression(config, groupedSCS_all, groupNames = NULL, overrideDir='DE_all/', outputMode = 'save')
clusterDifferentialExpression_noKM <- analyseClusterDifferentialExpression(config, groupedSCS_noKM, groupNames = NULL, overrideDir='DE_noKM/', outputMode = 'save')
clusterDifferentialExpression_noKMT <- analyseClusterDifferentialExpression(config, groupedSCS_noKMT, groupNames = NULL, overrideDir='DE_noKMT/', outputMode = 'save')


# -----
# Now also calculate "original" clustering etc with KCNQ+ cells removed
# This was referred to as "mod"

cellsToExclude = colnames(groupedSCS_all$patientAll@ndata)[which(as.double(groupedSCS_all$patientAll@ndata['KCNQ1OT1__chr11',])>10)]

# Re-analyse data without these cells
groupedData_all$patientAllMod = groupedData_all$patientAll[,which(!(colnames(groupedData_all$patientAll) %in% cellsToExclude))]
groupedData_all$patientAllWithIndexMod = groupedData_all$patientAllWithIndex[,which(!(colnames(groupedData_all$patientAllWithIndex) %in% cellsToExclude))]
groupedData_all$patient1Mod = groupedData_all$patient1[,which(!(colnames(groupedData_all$patient1) %in% cellsToExclude))]
groupedData_all$patient2Mod = groupedData_all$patient2[,which(!(colnames(groupedData_all$patient2) %in% cellsToExclude))]
groupedData_all$patient3Mod = groupedData_all$patient3[,which(!(colnames(groupedData_all$patient3) %in% cellsToExclude))]
groupedData_all$patient4Mod = groupedData_all$patient4[,which(!(colnames(groupedData_all$patient4) %in% cellsToExclude))]
groupedData_all$patient5Mod = groupedData_all$patient5[,which(!(colnames(groupedData_all$patient5) %in% cellsToExclude))]

# Build SCS-object
#groupedSCS_allMod <- buildSCSObject(config, groupedData_all, groupNames = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod','patientAllWithIndexMod', 'patientAllMod'))
groupedSCS_allMod <- buildSCSObject(config, groupedData_all, groupNames = c('patientAllMod'))
# Again differential expression
clusterDifferentialExpression_allMod <- analyseClusterDifferentialExpression(config, groupedSCS_allMod, groupNames = c('patientAllMod'), overrideDir='DE_patientAllMod/', outputMode = 'save')

######################################################################
# Some functionality to allow easy plotting

source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_standard_plots.R")
source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_analysis_functions.R")
source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_some_colors.R")

shorthand_XY_plot_patient = function(current_SCS,mytitle) {
    p1=plot_sth_in_XY(comp_1 = current_SCS@tsne[,1],
             comp_2 = current_SCS@tsne[,2],
             color_by = as.factor(current_SCS@cluster$kpart),
             name_for_color = 'cluster', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('r-tsne: ',mytitle))#+
        #theme(legend.position = 'none')
    return(p1)
}
shorthand_expression = function(some_gene,current_SCS,mytitle) {
    Expr=get_expression_gene2(
            expression_matrix = as.matrix(current_SCS@fdata),
            gene_query = some_gene)
    p1=plot_sth_in_XY(comp_1 = current_SCS@tsne[,1],
             comp_2 = current_SCS@tsne[,2],
             color_by = Expr,
             name_for_color = some_gene, print_yes = F)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne: ',mytitle))+
             scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(Expr),oob=squish)#+
             #theme(legend.position = 'none')
    print(p1)
    
    return(p1)
}

################################################################################
# Now just plot the tSNE maps

# Including all genes
shorthand_XY_plot_patient(groupedSCS_all$patientAll,'Original')
shorthand_expression('TTN',groupedSCS_all$patientAll,'Original')

# Including all genes, KCNQ1OT1+ cells excluded
shorthand_XY_plot_patient(groupedSCS_allMod$patientAllMod,'KCNQ1OT1+ excluded')
shorthand_expression('TTN',groupedSCS_allMod$patientAllMod,'KCNQ1OT1+ excluded')

# Ignoring K+M
shorthand_XY_plot_patient(groupedSCS_noKM$patientAll,'ignore KCNQ1OT1/MALAT')
shorthand_expression('TTN',groupedSCS_noKM$patientAll,'ignore KCNQ1OT1/MALAT')

# Now without TTN
shorthand_XY_plot_patient(groupedSCS_noKMT$patientAll,'ignore KCNQ1OT1/MALAT/TTN')
shorthand_expression('RYR2',groupedSCS_noKMT$patientAll,'ignore KCNQ1OT1/MALAT/TTN')

# Now retrieve TTN expression
# First check if colnames (cells) match
cells_to_fetch = colnames(groupedSCS_noKMT$patientAll@fdata)
cells_source   = colnames(groupedSCS_noKM$patientAll@fdata)
cells_found    = cells_to_fetch[cells_to_fetch %in% cells_source]

Expr=get_expression_gene2(
            expression_matrix = as.matrix(groupedSCS_noKM$patientAll@fdata[,cells_found]),
            gene_query = 'TTN')
p1=plot_sth_in_XY(   comp_1 = groupedSCS_noKMT$patientAll@tsne[cells_found,1],
                     comp_2 = groupedSCS_noKMT$patientAll@tsne[cells_found,2],
                     color_by = Expr,
                     name_for_color = some_gene, print_yes = F)+
                     theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne: ','ignoring TTN'))+
                     scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(Expr),oob=squish)
print(p1)



################################################################################
# Gradients are in tact, now see how the clustering matches

#-----
# First print clustering one on the other map
p2=plot_sth_in_XY(comp_1 = groupedSCS_noKM$patientAll@tsne[cells_found,1],
             comp_2 = groupedSCS_noKM$patientAll@tsne[cells_found,2],
             color_by = as.factor(groupedSCS_all$patientAll@cluster$kpart[cells_found]),
             name_for_color = 'cluster', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne(no-KM) with cluster(all-genes) on top'))
p2

#-----
# Repeat this for comparing removal of MALAT1/KCNQ
p1=plot_sth_in_XY(comp_1 = groupedSCS_noKM$patientAll@tsne[cells_found,1],
             comp_2 = groupedSCS_noKM$patientAll@tsne[cells_found,2],
             color_by = as.factor(groupedSCS_noKM$patientAll@cluster$kpart[cells_found]),
             name_for_color = 'cluster', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne(no-KMT) with cluster(no-KM) on top'))
p1

#-----
# Then alluvial



# Create alluvial plot
comparing_cluster_df = data.frame(clustering_noKM = groupedSCS_noKM$patientAll@cluster$kpart[cells_found],
        clustering_noKM_noTTN = groupedSCS_noKMT$patientAll@cluster$kpart[cells_found],
        freq=1)

alluvial_comparing_cluster_df = aggregate(comparing_cluster_df$freq, 
    by=list(clustering_noKM=comparing_cluster_df$clustering_noKM,
            clustering_noKM_noTTN=comparing_cluster_df$clustering_noKM_noTTN),
    FUN=sum)

library(ggalluvial)
ggplot(alluvial_comparing_cluster_df,
       aes(y = x, axis1 = clustering_noKM, axis2 = clustering_noKM_noTTN)) +
  geom_alluvium(aes(fill = as.factor(clustering_noKM)), width = 1/12)+
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_label(stat = "stratum", infer.label = TRUE)+theme_bw()+
    ggtitle('Alluvial to compare')
  #scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  #ggtitle("UC Berkeley admissions and rejections, by sex and department")

# -----
# Now difference in clustering between removing KM or not
cells1 = colnames(groupedSCS_all$patientAll@fdata)
cells2 = colnames(groupedSCS_noKM$patientAll@fdata)
cells_overlap12   = cells1[cells1 %in% cells2]

comparing_cluster_df = data.frame(
        clustering_all  = groupedSCS_all$patientAll@cluster$kpart[cells_overlap12],
        clustering_noKM = groupedSCS_noKM$patientAll@cluster$kpart[cells_overlap12],
        freq=1)

alluvial_comparing_cluster_df = aggregate(comparing_cluster_df$freq, 
    by=list(clustering_all=comparing_cluster_df$clustering_all,
            clustering_noKM=comparing_cluster_df$clustering_noKM),
    FUN=sum)

library(ggalluvial)
ggplot(alluvial_comparing_cluster_df,
       aes(y = x, axis1 = clustering_all, axis2 = clustering_noKM)) +
  geom_alluvium(aes(fill = as.factor(clustering_all)), width = 1/12)+
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_label(stat = "stratum", infer.label = TRUE)+theme_bw()+
    ggtitle('Alluvial to compare')


# See where KCN expression was:

Expr=get_expression_gene2(
            expression_matrix = as.matrix(groupedSCS_all$patientAll@fdata[,cells_overlap12]),
            gene_query = 'KCNQ1OT1')
p1=plot_sth_in_XY(   comp_1 = groupedSCS_noKM$patientAll@tsne[cells_overlap12,1],
                     comp_2 = groupedSCS_noKM$patientAll@tsne[cells_overlap12,2],
                     color_by = Expr,
                     name_for_color = 'KCNQ1OT1', print_yes = F)+
                     theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne: ','ignoring K/M'))+
                     scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(Expr),oob=squish)
                     #theme(legend.position = 'none')
print(p1)

################################################################################
# -----
# Perform quality control on the processed data
processedDataQualityControl(config, groupedSCS, groupNames = NULL, outputMode = 'save')
