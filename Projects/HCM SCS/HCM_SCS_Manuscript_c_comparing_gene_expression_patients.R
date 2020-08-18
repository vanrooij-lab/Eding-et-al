
# Analysis of differential gene expression per patient


# This script uses this param:
# groupedSCS$patientAllMod@ndata

# Get patient source and names
patient_source_Mod = shorthand_give_p_source(config, groupedSCS$patientAllMod)
patient_names = names(config$samples)

# determine in how many cells each gene is detected
gene_detected_counts = apply(groupedSCS$patientAllMod@ndata>.1, 1, sum)

# treshold for taking genes along (expression in X cells)
MIN_DETECTED = 100

# core of the script: determine mean expression per patient
# set up empty matrix, sized nr_genes*nr_patients
mean_gene_expression_per_patient = matrix(rep(NA, sum(gene_detected_counts>MIN_DETECTED)*length(patient_names)),ncol = length(patient_names))
# set appropriate row and column names
rownames(mean_gene_expression_per_patient) = names(gene_detected_counts[gene_detected_counts>MIN_DETECTED])
colnames(mean_gene_expression_per_patient) = patient_names
# Now determine the mean values of each gene per patient, loop over patients
for (idx in 1:length(patient_names)) {
    # mean gene expression (gene=row) on applicable sub-matrix
    mean_gene_expression_per_patient[,idx] = apply(
        groupedSCS$patientAllMod@ndata[gene_detected_counts>MIN_DETECTED,patient_source_Mod==paste0('patient',idx)], 1, mean)
}
# mean_gene_expression_per_patient

# Now scale to compare genes, divide each gene by their mean expression among patients
mean_gene_expression_per_patient_scaled = t(apply(mean_gene_expression_per_patient, 1, function(X) {X/mean(X)}))

# plot histogram of determined fold change values (scaled mean values = fold change)
ggplot(data.frame(fold_change=as.vector(mean_gene_expression_per_patient_scaled)))+
    geom_histogram(aes(x=fold_change))+theme_bw()

# determine max fold change per gene
max_fc = apply(mean_gene_expression_per_patient_scaled, 1, max)
# also determine mean as sanity check; should be ±1
mean_fc = apply(mean_gene_expression_per_patient_scaled, 2, mean) 

# View fc-table, ordered by max fold changes
View(mean_gene_expression_per_patient_scaled[order(max_fc, decreasing = T), ])

# Export to xls
library(openxlsx)
#wb <- createWorkbook()
my_list=list()
for (idx in 1:length(patient_names)) {
    
    top10_fc = as.data.frame(mean_gene_expression_per_patient_scaled[order(mean_gene_expression_per_patient_scaled[,idx], decreasing = T), ][1:10,])
    top10_fc$gene_name = rownames(top10_fc)
    my_list[[paste0('patient',idx)]] = top10_fc
    
}

my_list$all_data = as.data.frame(mean_gene_expression_per_patient_scaled[order(max_fc, decreasing = T), ])
my_list$all_data$gene_name = strip__chrXX(rownames(my_list$all_data))
write.xlsx(file=paste0(outputDir, '/patient_top_fc.xlsx'), x = my_list, overwrite = TRUE)





