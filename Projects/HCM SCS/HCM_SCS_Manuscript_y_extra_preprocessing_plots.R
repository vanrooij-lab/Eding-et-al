

######################################################################
# Please see notes in "HCM SCS Manuscript_AE.R" first

######################################################################
# Script May 2020
# ---
# This script provides some extra plots to facilitate
# choosing the right pre-processing steps, regarding filtering of 
# cells and genes that cause artifacts.
#


######################################################################
# Some functionality to allow easy plotting

library(stats)
library(fdrtool) # not sure required for mutoss?
library(mutoss) # (1) BiocManager::install("multtest") (2) install.packages('mutoss')
library('patchwork') # install.packages('patchwork')
source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_standard_plots.R")
source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_analysis_functions.R")
source("/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_some_colors.R")
source("/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/shorthandFunctionsMW.R")


######################################################################

# Show original tSNE (Including all genes)
shorthand_XY_plot_clusters(groupedSCS$patientAll,'Original')
shorthand_XY_plot_clusters(groupedSCS$patientAll,'Original', cl2 = T)
# Show expression on top of tSNE of a few genes
pK = shorthand_expression('KCNQ1OT1',groupedSCS$patientAll,'Original')
pM = shorthand_expression('MALAT1',groupedSCS$patientAll,'Original')
pNA = shorthand_expression('NPPA',groupedSCS$patientAll,'Original')
pNB = shorthand_expression('NPPB',groupedSCS$patientAll,'Original')
pM7 = shorthand_expression('MYH7',groupedSCS$patientAll,'Original')
shorthand_expression('NEAT1',groupedSCS$patientAll,'Original')
shorthand_expression('TTN',groupedSCS$patientAll,'Original')

# Exporting some plots
p_prettyK = pK+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(p_prettyK, thedir = 'MW_custom_plots/', filename = 'tSNE_beforefilter_KCNQ.png', mywidth = 7.5, myheight = 7.5)
p_prettyM = pM+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(p_prettyM, thedir = 'MW_custom_plots/', filename = 'tSNE_beforefilter_MALAT.png', mywidth = 7.5, myheight = 7.5)

# perhaps also quickly look at Volcano of these genes
volcano_df_MALAT1 = get_volcano_df(expression_matrix = as.matrix(groupedSCS$patientAll@ndata), my_gene = 'MALAT1', calc_qvals = F)
#p=plot_volcano(volcano_df_MALAT1, custom_highlight_group = genes_set_of_interest_sel_1, NRLABELED_hlgroup = 20)
p=plot_volcano(volcano_df_MALAT1, custom_highlight_group = shared_regulon_core_list[[1]], NRLABELED_hlgroup = 20)+give_better_textsize_plot(10)
p; shorthand_save(p,thedir = 'MW_custom_plots/', filename = 'volcano_MALAT1_custom.png', mywidth = 7, myheight = 7)

volcano_df_KCNQ = get_volcano_df(expression_matrix = as.matrix(groupedSCS$patientAll@ndata), my_gene = 'KCNQ1OT1', calc_qvals = F)
#p=plot_volcano(volcano_df_KCNQ, custom_highlight_group = genes_set_of_interest_sel_1, NRLABELED_hlgroup = 20)
p=plot_volcano(volcano_df_KCNQ, custom_highlight_group = shared_regulon_core_list[[1]], NRLABELED_hlgroup = 20, manual_gene_name = 'KCNQ1OT1')+give_better_textsize_plot(10)
p; shorthand_save(p,thedir = 'MW_custom_plots/', filename = 'volcano_KCNQ_custom.png', mywidth = 7, myheight = 7)

# Show expression per cluster
shorthand_expr_per_cluster('KCNQ1OT1', groupedSCS$patientAll, cl2 = T)
shorthand_expr_per_cluster('TTN', groupedSCS$patientAll, cl2 = T)
shorthand_expr_per_cluster('MYH7', groupedSCS$patientAll, cl2 = T)

# So let's add qvals and adjust pvals for differential expression of the clusters
for (clname in names(clusterDifferentialExpression$patientAll)) {
    clusterDifferentialExpression$patientAll[[clname]]$qv = pval2qval(clusterDifferentialExpression$patientAll[[clname]]$pv)[[1]]
    clusterDifferentialExpression$patientAll[[clname]]$pv.adj = p.adjust(clusterDifferentialExpression$patientAll[[clname]]$pv, method='BH')
}

# If you wanna view cluster stast on a gene:
#clusterDifferentialExpression$patientAll$cl.10['KCNQ1OT1__chr11',]

# Now let's make an overview plot
fc_K = sapply(clusterDifferentialExpression$patientAll, function(X) {  X['KCNQ1OT1__chr11',]$fc })
#pv_K = sapply(clusterDifferentialExpression$patientAll, function(X) {  X['KCNQ1OT1__chr11',]$pv })
qv_K = sapply(clusterDifferentialExpression$patientAll, function(X) {  X['KCNQ1OT1__chr11',]$qv })
pv_K = sapply(clusterDifferentialExpression$patientAll, function(X) {  X['KCNQ1OT1__chr11',]$pv.adj })
ggplot(data.frame(fold_change=fc_K, pv=pv_K, qv=qv_K, cluster=names(fc_K)))+
    geom_point(aes(x=fold_change, y=-log10(pv), color=cluster))+
    geom_text_repel(aes(x=fold_change, y=-log10(pv), label=cluster))+theme_bw()+
    geom_hline(yintercept = -log10(.05))+
    geom_vline(xintercept = 1)+give_better_textsize_plot(14)

CUTOFF=2
df_toplot=data.frame(fold_change=fc_K, pv=pv_K, qv=qv_K, cluster=names(fc_K))
ggplot(df_toplot)+
    geom_point(aes(x=log10(fold_change), y=-log10(pv), color=cluster))+
    #geom_text_repel(data=df_toplot[df_toplot$fold_change>CUTOFF,], aes(x=log10(fold_change), y=-log10(pv), label=cluster))+theme_bw()+
    geom_text_repel(aes(x=log10(fold_change), y=-log10(pv), label=cluster))+theme_bw()+
    geom_hline(yintercept = -log10(.05))+
    geom_vline(xintercept = log10(1))+
    geom_vline(xintercept = log10(CUTOFF),color='red')+give_better_textsize_plot(14)

ggplot(data.frame(fold_change=fc_K, pv=pv_K, qv=qv_K, cluster=names(fc_K)))+
    geom_point(aes(x=qv, y=-log10(pv), color=cluster))+
    geom_text_repel(aes(x=qv, y=-log10(pv), label=cluster))+theme_bw()+
    geom_hline(yintercept = -log10(.05))+
    geom_vline(xintercept = .1)

clusterDifferentialExpression

# Show histogram KCNQ1OT1 again
shorthand_histogram('KCNQ1OT1',groupedSCS$patientAll,'Original')
# Show histogram MALAT1 (MALAT1__chr11)
shorthand_histogram('MALAT1',groupedSCS$patientAll,'Original')

# Now for those genes listed in Gehart
#Gehart_genes = c('RN45s', 'MALAT1', 'KCNQ1OT1', 'A630089N07RIK', 'GM17821')
#Gehart_genes %in% strip__chrXX(rownames(groupedSCS$patientAll@ndata))
#Gehart_genes[Gehart_genes %in% strip__chrXX(rownames(groupedSCS$patientAll@ndata))]
# These are not 
shorthand_histogram('MALAT1',groupedSCS$patientAll,'Original')





# Make a seleciton based on KCNQ1OT1
gene_names=rownames(groupedSCS$patientAll@ndata)
KCNQ1OT1_sel = as.double(groupedSCS$patientAll@ndata['KCNQ1OT1__chr11',])>10
shorthand_selection(KCNQ1OT1_sel,groupedSCS$patientAll,'Original','KCNQ1OT1')

# Make a selection based on KCNQ1OT1 enriched clusters
cl_sel = groupedSCS$patientAll@cpart %in% c(10,11,13,16)
shorthand_selection(cl_sel,groupedSCS$patientAll,'Original','bad_clusters')

# MYH7 low clusters
cl_sel_M7 = groupedSCS$patientAll@cpart %in% c(9,14,15,16,17,18)
shorthand_selection(cl_sel_M7,groupedSCS$patientAll,'Original','low_MYH7')

# Combination of those two
shorthand_selection(KCNQ1OT1_sel|cl_sel,groupedSCS$patientAll,'Original','combined')

# Show some MALAT1 stats
print(find__chrXX('MALAT1', gene_names))
MALAT1_sel = as.double(groupedSCS$patientAll@ndata['MALAT1__chr11',])>10
shorthand_selection(MALAT1_sel,groupedSCS$patientAll,'Original','MALAT1')

# Combination of KCNQ1OT1 & MALAT1
shorthand_selection(MALAT1_sel|KCNQ1OT1_sel,groupedSCS$patientAll,'MALAT1+ or KCNQ1OT1+','MALAT_KCN')

################################################################################
# Now show the result of the clustering after selection

shorthand_XY_plot_clusters(groupedSCS$patientAllMod,'Original')
shorthand_XY_plot_clusters(groupedSCS$patientAllMod,'Original', cl2 = T)
 
# Show were patients are   
plot_sth_in_XY(comp_1 = groupedSCS$patientAllMod@tsne[,1],
             comp_2 = groupedSCS$patientAllMod@tsne[,2],
             color_by = shorthand_give_p_source(config, groupedSCS$patientAllMod),
             name_for_color = 'patient', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)
shorthand_expression('TTN',groupedSCS$patientAllMod,'Filtered')
shorthand_expression('MYH7',groupedSCS$patientAllMod,'Filtered')
shorthand_expression('NPPA',groupedSCS$patientAllMod,'Filtered')
shorthand_expression('NPPB',groupedSCS$patientAllMod,'Filtered')


my_genes = c('ZNF106', 'CHGA', 'MTRNR2L1', 'FILIP1', 'ANKRD2', 'TCAP', 'NPPA', 'VIM')
p=list()
for (idx in 1:length(my_genes)) {
    current_gene=my_genes[idx]
    p[[idx]]=shorthand_expression(current_gene,groupedSCS$patientAllMod,'Filtered')
    print(p[[idx]])
}
wrap_plots(p)


# Highlight clusters 9-14
cl_9_14 = groupedSCS$patientAllMod@cpart %in% 9:14
shorthand_selection(cl_9_14,groupedSCS$patientAllMod,'Original','cl_9_14')

# Some notes about the clustering:
groupedSCS$patientAllMod@cluster$kpart[groupedSCS$patientAllMod@cpart %in% 9:14]
# "groupedSCS$patientAllMod@cluster$kpart" is identical to "groupedSCS$patientAllMod@cpart"
# except that the clusters with few cells in groupedSCS$patientAllMod@cpart are re-assigned
# to larger clusters, such that cluster$kpart has fewer clusters. The cluster membership
# of cells belonging to large clusters in both is conserved between these two.
any(!groupedSCS$patientAllMod@cluster$kpart[!(groupedSCS$patientAllMod@cpart %in% 9:14)]==
        groupedSCS$patientAllMod@cpart[!(groupedSCS$patientAllMod@cpart %in% 9:14)])

# Highlight cluster 8 cells in original map
cl8_names = names(groupedSCS$patientAllMod@cluster$kpart)[groupedSCS$patientAllMod@cluster$kpart %in% 8]
cl8_in_orig = colnames(groupedSCS$patientAll@ndata) %in% cl8_names
shorthand_selection(cl8_in_orig,groupedSCS$patientAll,'Original','in_new_cl8')

################################################################################
# Now look at cluster of patients with index data

shorthand_XY_plot_clusters(groupedSCS$patientAllWithIndexMod,'Original')
#shorthand_XY_plot_clusters(groupedSCS$patientAllWithIndexMod,'Original', cl2 = T)

# Show were patients are   
plot_sth_in_XY(comp_1 = groupedSCS$patientAllWithIndexMod@tsne[,1],
             comp_2 = groupedSCS$patientAllWithIndexMod@tsne[,2],
             color_by = shorthand_give_p_source(config, groupedSCS$patientAllWithIndexMod),
             name_for_color = 'patient', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)


