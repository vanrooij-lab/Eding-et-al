
# I first wanted to calculate the overlap between detected genes in patients by comparing highly expressed genes;
# However, this might not be the best way, as around the cutoff, there might be many genes that have an expression
# just above the cutoff in one patient, and just below it in another. 
# 
# It might be better to calculate the differential expression of genes between patients.

# Let's look at the differential gene expression between patients.

# Filter genes that are found in X% of cells (note that mean of {1,0} is fraction expressed)
MINPERCENTAGCELLSEEXPRESSED=.2
genes_pt1 = strip__chrXX(rownames(groupedSCS$patient1Mod@ndata[apply(groupedSCS$patient1Mod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_pt2 = strip__chrXX(rownames(groupedSCS$patient2Mod@ndata[apply(groupedSCS$patient2Mod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_pt3 = strip__chrXX(rownames(groupedSCS$patient3Mod@ndata[apply(groupedSCS$patient3Mod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_pt4 = strip__chrXX(rownames(groupedSCS$patient4Mod@ndata[apply(groupedSCS$patient4Mod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_pt5 = strip__chrXX(rownames(groupedSCS$patient5Mod@ndata[apply(groupedSCS$patient5Mod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
# Filter genes with a minimal average expression (also done for GO terms)
MIN_EXPRESSION=.2 # I find it hard to interpret what this cutoff really means; but the #genes selected with mean>.2
                  # is roughly similar to the above cutoff of genes expressed in >20% of cells
ggplot(data.frame(mean_expression=apply(groupedSCS$patientAllMod@ndata,1,mean)))+
    geom_histogram(aes(x=log10(.1+mean_expression)), bins=30)+theme_bw()+geom_vline(xintercept = log10(.1+MIN_EXPRESSION))+
    ggtitle('Gene selection')
genes_pt1 = strip__chrXX(rownames(groupedSCS$patient1Mod@ndata[apply(groupedSCS$patient1Mod@ndata, 1, mean)>MIN_EXPRESSION,]))
genes_pt2 = strip__chrXX(rownames(groupedSCS$patient2Mod@ndata[apply(groupedSCS$patient2Mod@ndata, 1, mean)>MIN_EXPRESSION,]))
genes_pt3 = strip__chrXX(rownames(groupedSCS$patient3Mod@ndata[apply(groupedSCS$patient3Mod@ndata, 1, mean)>MIN_EXPRESSION,]))
genes_pt4 = strip__chrXX(rownames(groupedSCS$patient4Mod@ndata[apply(groupedSCS$patient4Mod@ndata, 1, mean)>MIN_EXPRESSION,]))
genes_pt5 = strip__chrXX(rownames(groupedSCS$patient5Mod@ndata[apply(groupedSCS$patient5Mod@ndata, 1, mean)>MIN_EXPRESSION,]))
# What is the aim?
# (a) avoid artifacts from genes just above detection limit
# (b) 
# So I think the way below is the way to go, first make a general selection, and then take along
# only those that are found in all patients
MINPERCENTAGCELLSEEXPRESSED=.2
genes_detected = strip__chrXX(rownames(groupedSCS$patientAllMod@ndata[apply(groupedSCS$patientAllMod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_detected_all_patients_mincelltreshold = genes_detected

# Create a manual matrix
#all_medium_expressed_genes = unique(c(genes_pt1, genes_pt2, genes_pt3, genes_pt4, genes_pt5))
p1_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS$patient1Mod@ndata)))
p2_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS$patient2Mod@ndata)))
p3_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS$patient3Mod@ndata)))
p4_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS$patient4Mod@ndata)))
p5_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS$patient5Mod@ndata)))

manual_upset_matrix = data.frame(p1=p1_yes, p2=p2_yes, p3=p3_yes, p4=p4_yes, p5=p5_yes)
rownames(manual_upset_matrix) = genes_detected # all_medium_expressed_genes
manual_upset_matrix_genesexpressedpatients = manual_upset_matrix

library('UpSetR')
#dir.create(paste0(outputDir,'/Fig2/C/'), recursive=T)
#pdf(paste0(outputDir,'/XXX/C/TTN_cluster_overlap.pdf'))
#all_expressed_genes = list(
#  p1 = genes_pt1, p2 = genes_pt2, p3 = genes_pt3, p4 = genes_pt4, p5 = genes_pt5
#)
#UpSetData <- UpSetR::fromList(all_expressed_genes)
#UpSetData$n <- sample(1:nrow(UpSetData))
upset(
  manual_upset_matrix,
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
#dev.off()

# So, which genes are so patient-specific?
for (patient_idx in 1:5) {
    print(paste0('patient ', patient_idx))
    print(toString(rownames(manual_upset_matrix)[manual_upset_matrix[,patient_idx] & !apply(manual_upset_matrix[,-patient_idx],1,any)]))
}







