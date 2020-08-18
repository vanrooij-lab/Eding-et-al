library(mutoss)
library(patchwork)

TTN_volcano_df_pts = list()
for (current_patient in c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod')) {
    TTN_volcano_df_pts[[current_patient]] = 
        get_volcano_df(expression_matrix = as.matrix(groupedSCS[[current_patient]]@ndata), my_gene = 'TTN',calc_qvals = F,no_expression_val = 0.1)
}
    
# Filter out correlations that are based on observations in <10 cells
# This isn't necessary any more since I increased the treshold to a fraction of 10% detected
#TTN_volcano_df_pts_sel=list()
#for (pt_name in names(TTN_volcano_df_pts)){
#  TTN_volcano_df_pts_sel[[pt_name]] = TTN_volcano_df_pts[[pt_name]][TTN_volcano_df_pts[[pt_name]]$nrCellsExpressGene>=10,]
#}


(plot_volcano(TTN_volcano_df_pts[['patient1Mod']])+ggtitle('Patient 1'))+
(plot_volcano(TTN_volcano_df_pts[['patient2Mod']])+ggtitle('Patient 2'))+
(plot_volcano(TTN_volcano_df_pts[['patient3Mod']])+ggtitle('Patient 3'))+
(plot_volcano(TTN_volcano_df_pts[['patient4Mod']])+ggtitle('Patient 4'))+
(plot_volcano(TTN_volcano_df_pts[['patient5Mod']])+ggtitle('Patient 5'))+
plot_layout(nrow = 2)

View(expressionData_FSCA)


# Save extra plots with "Joep-look"
patients_of_interest = c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod')
for (idx in 1:length(patients_of_interest)) {
  current_patient= patients_of_interest[idx]
  p=plot_volcano(TTN_volcano_df_pts[[current_patient]], mycex = 7, mypointsize = 1,mylinesize = 1)+ggtitle(paste0('Patient ',idx))+theme_classic()+give_better_textsize_plot(24)+xlim(c(-.8,.8))+shorthand_tsne_joeptheme()
  shorthand_save(p, 'analysis_2020-05-04-18h_17m/FigX/', paste0('mw_TTN_Vulcano_',current_patient,'.pdf'), mywidth = 17.86, myheight = 14.27)
}
# Save extra plots with "Joep-look" (small)
patients_of_interest = c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod')
for (idx in 1:length(patients_of_interest)) {
  current_patient= patients_of_interest[idx]
  p=plot_volcano(TTN_volcano_df_pts[[current_patient]], mycex = 1.5, mypointsize = .1,mylinesize = .25)+ggtitle(paste0('Patient ',idx))+theme_classic()+give_better_textsize_plot(7)+xlim(c(-.8,.8))+shorthand_tsne_joeptheme()
  shorthand_save(p, 'analysis_2020-05-04-18h_17m/FigX/', paste0('mw_TTN_Vulcano_',current_patient,'_small.pdf'), mywidth = 5, myheight = 4.33)
}
# note, US letter = 215.9 x 279.4 mm

# Also save the correlations with TTN to an Excel sheet
excel_export_df_TTN = TTN_volcano_df_pts # create a new var just to have nice names
names(excel_export_df_TTN)=gsub('Mod','',names(excel_export_df_TTN))
openxlsx::write.xlsx(excel_export_df_TTN,paste0(getwd(),'/analysis_2020-05-04-18h_17m/FigX/TTN_correlated_genes.xlsx'))
rm('excel_export_df_TTN') # remove extra paramater 

########

# Top X
TTN_topX_pt1 = TTN_volcano_df_pts[['patient1Mod']][order(TTN_volcano_df_pts[['patient1Mod']]$corr, decreasing = T),]$gene_name_short[2:31]
TTN_topX_pt2 = TTN_volcano_df_pts[['patient2Mod']][order(TTN_volcano_df_pts[['patient2Mod']]$corr, decreasing = T),]$gene_name_short[2:31]
TTN_topX_pt3 = TTN_volcano_df_pts[['patient3Mod']][order(TTN_volcano_df_pts[['patient3Mod']]$corr, decreasing = T),]$gene_name_short[2:31]
TTN_topX_pt4 = TTN_volcano_df_pts[['patient4Mod']][order(TTN_volcano_df_pts[['patient4Mod']]$corr, decreasing = T),]$gene_name_short[2:31]
TTN_topX_pt5 = TTN_volcano_df_pts[['patient5Mod']][order(TTN_volcano_df_pts[['patient5Mod']]$corr, decreasing = T),]$gene_name_short[2:31]

# Create a manual matrix
#all_medium_expressed_genes = unique(c(genes_pt1, genes_pt2, genes_pt3, genes_pt4, genes_pt5))
p_common = unique(c(TTN_topX_pt1,TTN_topX_pt2,TTN_topX_pt3,TTN_topX_pt4,TTN_topX_pt5))
p1_yes = as.double(p_common %in% TTN_topX_pt1)
p2_yes = as.double(p_common %in% TTN_topX_pt2)
p3_yes = as.double(p_common %in% TTN_topX_pt3)
p4_yes = as.double(p_common %in% TTN_topX_pt4)
p5_yes = as.double(p_common %in% TTN_topX_pt5)

manual_upset_matrix_TTNtop = data.frame(p1=p1_yes, p2=p2_yes, p3=p3_yes, p4=p4_yes, p5=p5_yes)
rownames(manual_upset_matrix_TTNtop) = p_common # all_medium_expressed_genes

library('UpSetR')
#dir.create(paste0(outputDir,'/Fig2/C/'), recursive=T)
#pdf(paste0(outputDir,'/XXX/C/TTN_cluster_overlap.pdf'))
#all_expressed_genes = list(
#  p1 = genes_pt1, p2 = genes_pt2, p3 = genes_pt3, p4 = genes_pt4, p5 = genes_pt5
#)
#UpSetData <- UpSetR::fromList(all_expressed_genes)
#UpSetData$n <- sample(1:nrow(UpSetData))
p=upset(
  manual_upset_matrix_TTNtop,
  empty.intersections = NULL,
  order.by = "freq",
  group.by = "degree",
  nsets = 36,
  nintersects = 252,
  mb.ratio = c(0.5,0.5), 
  queries = list(
    list(query = intersects, params = list('p1','p2','p3','p4','p5'), color = "green", active = T)
  )
)

# save it
pdf(file = paste0(getwd(), '/analysis_2020-05-04-18h_17m/FigX/TTN_TopX_shared.pdf'), width = .39*10, height=.39*10)
p
dev.off()

# doesn't work
#shorthand_save(p, 
#    '/analysis_2020-05-04-18h_17m/FigX/', paste0('TTN_TopX_shared.pdf'),mywidth = 10, myheight = 10)

# genes correlated with TTN in all patients
#toString(p_common[apply(manual_upset_matrix_TTNtop, 1, all)])
toString(p_common[apply(manual_upset_matrix_TTNtop, 1, sum)==5])
toString(p_common[apply(manual_upset_matrix_TTNtop, 1, sum)>=4])
toString(p_common[apply(manual_upset_matrix_TTNtop, 1, sum)>=3])
toString(p_common[apply(manual_upset_matrix_TTNtop, 1, sum)<3])
#toString(p_common[apply(manual_upset_matrix_TTNtop[, -4], 1, all)])

print(paste0('nr of top-30 genes overlapping between 4+ patients: ',sum(apply(manual_upset_matrix_TTNtop, 1, sum)>=4)))




