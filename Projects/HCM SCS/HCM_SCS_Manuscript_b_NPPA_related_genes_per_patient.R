

# Since we're doing TTN per patient, perhaps also NPPA should be done per patient.
# This since a caveat of the pooled data might be that patient-to-patient 
# (biological) variations dominate the correlations.
#
# Unfortunately indepedent analysis per patient yields inconsistent results.
#
# There is however some consolidation in the fact that when we first determine
# average gene expression per patient, and then determine correlations (ie correlations
# purely based on patient-to-patient variations), those are not the same as the ones
# from the pooled expressoin data.
#
# This seems to suggest that determining correlations from the pooled patient data does offer
# some additional resolution. 

NPPA_volcano_df_pts = list()
for (current_patient in c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod')) {
    NPPA_volcano_df_pts[[current_patient]] = 
        get_volcano_df(expression_matrix = as.matrix(groupedSCS[[current_patient]]@ndata), my_gene = 'NPPA',calc_qvals = F,no_expression_val = 0.1)
}
    
# Filter out correlations that are based on observations in <10 cells
NPPA_volcano_df_pts_sel=list()
for (pt_name in names(NPPA_volcano_df_pts)){
  NPPA_volcano_df_pts_sel[[pt_name]] = NPPA_volcano_df_pts[[pt_name]][NPPA_volcano_df_pts[[pt_name]]$nrCellsExpressGene>=10,]
}

pNPPA=
(plot_volcano(NPPA_volcano_df_pts_sel[['patient1Mod']])+ggtitle('Patient 1'))+give_better_textsize_plot(8)+
(plot_volcano(NPPA_volcano_df_pts_sel[['patient2Mod']])+ggtitle('Patient 2'))+give_better_textsize_plot(8)+
(plot_volcano(NPPA_volcano_df_pts_sel[['patient3Mod']])+ggtitle('Patient 3'))+give_better_textsize_plot(8)+
(plot_volcano(NPPA_volcano_df_pts_sel[['patient4Mod']])+ggtitle('Patient 4'))+give_better_textsize_plot(8)+
(plot_volcano(NPPA_volcano_df_pts_sel[['patient5Mod']])+ggtitle('Patient 5'))+give_better_textsize_plot(8)+
plot_layout(nrow = 2)

pNPPA
shorthand_save(pNPPA, 'analysis_2020-05-04-18h_17m/FigX/', paste0('mw_NPPA_Vulcano_',current_patient,'_all.pdf'), mywidth = 20, myheight = 20)

#######

# Top X
NPPA_topX_pt1 = NPPA_volcano_df_pts_sel[['patient1Mod']][order(NPPA_volcano_df_pts_sel[['patient1Mod']]$corr, decreasing = T),]$gene_name_short[2:51]
NPPA_topX_pt2 = NPPA_volcano_df_pts_sel[['patient2Mod']][order(NPPA_volcano_df_pts_sel[['patient2Mod']]$corr, decreasing = T),]$gene_name_short[2:51]
NPPA_topX_pt3 = NPPA_volcano_df_pts_sel[['patient3Mod']][order(NPPA_volcano_df_pts_sel[['patient3Mod']]$corr, decreasing = T),]$gene_name_short[2:51]
NPPA_topX_pt4 = NPPA_volcano_df_pts_sel[['patient4Mod']][order(NPPA_volcano_df_pts_sel[['patient4Mod']]$corr, decreasing = T),]$gene_name_short[2:51]
NPPA_topX_pt5 = NPPA_volcano_df_pts_sel[['patient5Mod']][order(NPPA_volcano_df_pts_sel[['patient5Mod']]$corr, decreasing = T),]$gene_name_short[2:51]

# Create a manual matrix
#all_medium_expressed_genes = unique(c(genes_pt1, genes_pt2, genes_pt3, genes_pt4, genes_pt5))
p_common = unique(c(NPPA_topX_pt1,NPPA_topX_pt2,NPPA_topX_pt3,NPPA_topX_pt4,NPPA_topX_pt5))
p1_yes = as.double(p_common %in% NPPA_topX_pt1)
p2_yes = as.double(p_common %in% NPPA_topX_pt2)
p3_yes = as.double(p_common %in% NPPA_topX_pt3)
p4_yes = as.double(p_common %in% NPPA_topX_pt4)
p5_yes = as.double(p_common %in% NPPA_topX_pt5)

manual_upset_matrix_NPPAtop = data.frame(p1=p1_yes, p2=p2_yes, p3=p3_yes, p4=p4_yes, p5=p5_yes)
rownames(manual_upset_matrix_NPPAtop) = p_common # all_medium_expressed_genes

library('UpSetR')
#dir.create(paste0(outputDir,'/Fig2/C/'), recursive=T)
#pdf(paste0(outputDir,'/XXX/C/NPPA_cluster_overlap.pdf'))
#all_expressed_genes = list(
#  p1 = genes_pt1, p2 = genes_pt2, p3 = genes_pt3, p4 = genes_pt4, p5 = genes_pt5
#)
#UpSetData <- UpSetR::fromList(all_expressed_genes)
#UpSetData$n <- sample(1:nrow(UpSetData))
p=upset(
  manual_upset_matrix_NPPAtop,
  #empty.intersections = NULL,
  #order.by = "freq",
  #group.by = "degree",
  #nsets = 36,
  #nintersects = 252,
  #mb.ratio = c(0.5,0.5), 
  #queries = list(
  #  list(query = intersects, params = list('p1','p2','p3','p4','p5'), color = "green", active = T)
  #)
)

# save it
pdf(file = paste0(getwd(), '/analysis_2020-05-04-18h_17m/FigX/NPPA_TopX_shared.pdf'), width = .39*10, height=.39*10)
p
dev.off()

# doesn't work
#shorthand_save(p, 
#    '/analysis_2020-05-04-18h_17m/FigX/', paste0('NPPA_TopX_shared.pdf'),mywidth = 10, myheight = 10)

# genes correlated with NPPA in all patients
#toString(p_common[apply(manual_upset_matrix_NPPAtop, 1, all)])
toString(p_common[apply(manual_upset_matrix_NPPAtop, 1, sum)==5])
toString(p_common[apply(manual_upset_matrix_NPPAtop, 1, sum)>=4])
toString(p_common[apply(manual_upset_matrix_NPPAtop, 1, sum)>=3])
toString(p_common[apply(manual_upset_matrix_NPPAtop, 1, sum)<3])

################################################################################

# Unfortunately, there doesn't seem to be much overlap from one patient to the next..
# So, let's see what happens if we would determine correlations based on patient-to-patient
# variation

#
all_gene_names = rownames(groupedSCS$patientAllMod@ndata)
# or another criterium, at least in 10 cells
genes_detected_min10Cells = rownames(groupedSCS$patientAllMod@ndata[apply(groupedSCS$patientAllMod@ndata>.1, 1, sum)>=10,])
# or use this one (testing showed this resulted in the highest overlap with pooled result):
all_genes_minExpr_chrXX = find__chrXX(genes_detected_all_patients_mincelltreshold, all_gene_names)

patient_average_expr=list()
patient_average_expr_expanded=list()
for (patient_name in c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod')) {
    patient_average_expr[[patient_name]] = apply(groupedSCS[[patient_name]]@ndata, 1, mean)
    patient_average_expr_expanded[[patient_name]] = patient_average_expr[[patient_name]][all_genes_minExpr_chrXX]
    names(patient_average_expr_expanded[[patient_name]]) = all_genes_minExpr_chrXX
}

for (patient_name in c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod')) {
    patient_average_expr_expanded[[patient_name]][is.na(patient_average_expr_expanded[[patient_name]])] = 
        0.1
}

patient_avg_expr_df = data.frame(patient_average_expr_expanded)
View(patient_avg_expr_df)

NPPA_avg_pt_volcano = get_volcano_df(expression_matrix = as.matrix(patient_avg_expr_df), 
    my_gene = 'NPPA',calc_qvals = F,no_expression_val = 0.1)
plot_volcano(NPPA_avg_pt_volcano)
View(NPPA_avg_pt_volcano)

NPP_top30_avg_pt = NPPA_avg_pt_volcano[order(NPPA_avg_pt_volcano$corr, decreasing = T),][1:30,]$gene_name
NPP_top30_pooled = correlationResultNPPA$correlationData[order(correlationResultNPPA$correlationData$correlation, decreasing = T),][1:30,]$genes

# so the overlapping genes found in 
sum(NPP_top30_avg_pt %in% NPP_top30_pooled)
