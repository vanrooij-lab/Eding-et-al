
################################################################################
# 
# This scripts generates regulons per patients, using Joep's functions.
# 
# It basically does it twice, the second time i think in a better way,
# taking along the same gene set, and with lower (but more appropriate) tresholds
# (The second part is also labeled as such, i.e. "Part 2")


################################################################################
# ggnet2 / network plotting stuff
library(sna)
library(network)
library(GGally) # needed for ggnet2

# Patient specific regulons

home_outputDir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/'

################################################################################
# Calculate correlation matrices for all patients
correlationMatrices = generateCorrelationMatrix(config, groupedSCS, groupNames=c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod'), excludeOutlierCells=T, minCellFraction=0.05, 
                                            minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf')
# Determine the regulons    
regulons <- analyzeRegulons(config, correlationMatrices, minCorrelations=40, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig4/A', outputMode='pdf', K.max = 50)
# View some stuff
View(regulons$patient4Mod$regulons)

# now create a parameter which holds all regulolons in groups
joined_regulons_ = unlist(lapply(regulons, function(X) {X$regulons}), recursive = F)
joined_regulons = lapply(joined_regulons_, strip__chrXX)

# make the names shorter
names(joined_regulons)=sub(pattern = 'patient', replacement = 'p', x=names(joined_regulons))
names(joined_regulons)=sub(pattern = 'Mod.regulon.', replacement = '.R', x=names(joined_regulons))

# create distance matrix
dist_matrix = matrix(rep(NA, length(joined_regulons)*length(joined_regulons)), nrow=length(joined_regulons))
for (reg_idx in 1:length(joined_regulons)) {
    dist_matrix[reg_idx, ] = 
        sapply(joined_regulons, function(X) {sum(joined_regulons[[reg_idx]] %in% X)/min(length(X), length(joined_regulons[[reg_idx]]))})
}
# give distance matrix names
rownames(dist_matrix) = names(joined_regulons)
colnames(dist_matrix) = names(joined_regulons)

# show heatmap of distance matrix
pheatmap(dist_matrix, cluster_cols = F, cluster_rows = F)
pheatmap(dist_matrix)

overlap_other_regulons = apply(dist_matrix,1,sum)
pheatmap(dist_matrix[overlap_other_regulons>1.1,overlap_other_regulons>1.1])

# show network
regulon_network = network(dist_matrix, directed = FALSE)
ggnet2(regulon_network, node.size = 1, label=T, label.size = 5) # simple plotting
# show network using cutoff
regulon_network = network(dist_matrix>.25, directed = FALSE)
ggnet2(regulon_network, node.size = 1, label=T, label.size = 5) # simple plotting

# So it appears that especially p1.R6/p2.R3/p3.R8/p5.R4 are very similar
similar_regulon_set1 = c('p1.R6','p2.R3','p3.R8','p5.R4')
pheatmap(dist_matrix[similar_regulon_set1,similar_regulon_set1], cluster_cols = F, cluster_rows = F)
dist_matrix[similar_regulon_set1,similar_regulon_set1]
lapply(similar_regulon_set1, function(X) {joined_regulons[[X]]})
genes_from_similar_regulons_set1 = unique(unlist(lapply(similar_regulon_set1, function(X) {joined_regulons[[X]]})))
toString(genes_from_similar_regulons_set1)
# Now create a little table

library("org.Hs.eg.db")
library("AnnotationDbi") # mapIds
#library(Biobase)
#getEG(x = 'ANKRD1',data = 'org.Hs.eg')
#lookUp(x = 'ANKRD1',data = 'org.Hs.eg', "EG")

geneIDs <- mapIds(org.Hs.eg.db, keys='ANKRD1', column="ENTREZID", keytype="SYMBOL", multiVals="first")
geneID_Ensemble <- mapIds(org.Hs.eg.db, keys='ANKRD1', column="ENSEMBL", keytype="SYMBOL", multiVals="first")
# ENTREZID, SYMBOL
# Use columns(org.Hs.eg.db) to get a list of possible, ENTREZID, ENSEMBL and SYMBOL are probably most relevant

################################################################################
# Homer analysis

# To execute a Homer analysis, (1) install Homer:
# See: http://homer.ucsd.edu/homer/introduction/install.html
# (2) Create a text file with the gene names on each row
# (See instructions also on http://homer.ucsd.edu/homer/microarray/index.html)
# (3) Execute: findMotifs.pl input_file_path human output_dir_path
# the value "human" here indicates (obviously) the organism, but the parameter is also referred to as "promoter set"

################################################################################
################################################################################
# Part 2 || Part II
################################################################################
################################################################################

################################################################################
# Let's do another thing, and take genes that are expressed in all patients,
# and make a low-treshold correlation matrix based on those
# Note 1: the analysis below was finally used for the manuscript
# Note 2: the idea here is that we again generate a corr matrix for each
# of the patients separately, but we only look at genes expressed in all
# of the patients; this slightly lowers the amount of genes considered.
# Therefor, but also because in general that might be better, we decrease
# the thresholds used previously to take into account more genes.

######

# In the script "HCM_SCS_Manuscript_c_comparing_gene_detection_patients.R" it is determined
# which genes that are in the pooled count tables are found in 20% of cells, are also 
# found (=very minimal demands from config parameter selection criteria for taking along
# genes) in which patients.
# 
# Now let's identify genes that are expressed in all patients
all_genes_patients = rownames(manual_upset_matrix_genesexpressedpatients)[apply(manual_upset_matrix_genesexpressedpatients,1,all)]
sum(unique(unlist(joined_regulons)) %in% all_genes_patients)
# As a side/more general note, one can also wonder how many genes are found in e.g. at least 3 patients
genes_at_least_3_pts = rownames(manual_upset_matrix_genesexpressedpatients)[apply(manual_upset_matrix_genesexpressedpatients,1,sum)>=3]
sum(unique(unlist(joined_regulons)) %in% genes_at_least_3_pts)

# Now take those genes (this line just renames them w. the chromosome suffix)
all_genes_patients_chrXX = find__chrXX(all_genes_patients, rownames(groupedSCS$patientAll@ndata))

# Calculate corr matrix based on lower p-val and only genes expressed in all patients
correlationMatrices_selectedGenes = generateCorrelationMatrix(config, groupedSCS, which_genes_to_select = all_genes_patients_chrXX, groupNames=c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod'), excludeOutlierCells=T, minCellFraction=0.05, 
                                            minCellExpression=0.1, desiredPValue=0.001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf', filename_precursor='genes_allPs')
regulons_selectedGenes <- analyzeRegulons(config, correlationMatrices_selectedGenes, minCorrelations=10, clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, overrideDir='Fig4/A', outputMode='pdf', K.max = 50, filename_precursor='genes_allPs')

# Some stats that are mentioned in the manuscript:
# Nr of genes that are selected based on the minCellExpression criterium

# Nr of genes that are selected based on the minCorrelations criterium
sapply(names(regulons_selectedGenes), function(X) {dim(regulons_selectedGenes[[X]]$subCorrelationMatrix)[1]})

# Print the nr of regulons per patient
nr_regs_pt = sapply(names(regulons_selectedGenes), function(X) {length(regulons_selectedGenes[[X]]$regulons)})
toString(nr_regs_pt)

# So now repeat what was done above, ie create an overview of which regulons are similar between patients
# ---

# now create a parameter which holds all regulolons in groups
joined_regulons_sel_ = unlist(lapply(regulons_selectedGenes, function(X) {X$regulons}), recursive = F)
joined_regulons_sel = lapply(joined_regulons_sel_, strip__chrXX)

# make the names shorter
names(joined_regulons_sel)=sub(pattern = 'patient', replacement = 'p', x=names(joined_regulons_sel))
names(joined_regulons_sel)=sub(pattern = 'Mod.regulon.', replacement = '.R', x=names(joined_regulons_sel))

# create distance matrix
dist_matrix_sel = matrix(rep(NA, length(joined_regulons_sel)*length(joined_regulons_sel)), nrow=length(joined_regulons_sel))
for (reg_idx in 1:length(joined_regulons_sel)) {
    dist_matrix_sel[reg_idx, ] = 
        sapply(joined_regulons_sel, function(X) {sum(joined_regulons_sel[[reg_idx]] %in% X)/min(length(X), length(joined_regulons_sel[[reg_idx]]))})
}
# give distance matrix names
rownames(dist_matrix_sel) = names(joined_regulons_sel)
colnames(dist_matrix_sel) = names(joined_regulons_sel)

# show heatmap of distance matrix
pheatmap(dist_matrix_sel, cluster_cols = F, cluster_rows = F)
# cluster
hclust_out = hclust(dist(dist_matrix_sel), method='ward.D2')
plot(as.dendrogram(hclust_out))
cutree_out = cutree(hclust_out, h = 3)
cutree_df  = as.data.frame(as.factor(cutree_out)); colnames(cutree_df) = c('group')
#annotation_colors = col_Dark2[1:max(cutree_out)]
# Create a little heatmap
annotation_colors = col_vector_60[1:max(cutree_out)]
names(annotation_colors) = unique(cutree_out)
annotation_colors=list(group=annotation_colors)
p=pheatmap(dist_matrix_sel, cluster_rows = hclust_out,cluster_cols = hclust_out, 
    annotation_col = cutree_df, annotation_row = cutree_df, annotation_colors = annotation_colors)
    #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
shorthand_save(p, 'MW_custom_plots/', 'regulons_heatmap_shared_clusters.png', mywidth = 15, myheight = 12)

# 
for (i in 1:5) {
  overlaps=as.vector(dist_matrix_sel[cutree_out==i,cutree_out==i])
  print(paste0('overlap reg.', i,'= ', 
    round(100*median(overlaps[overlaps<1])),'%'))
}


# an example of overlap
toString(joined_regulons_sel$p1.R1)
toString(joined_regulons_sel$p2.R2)
joined_regulons_sel$p1.R1[joined_regulons_sel$p1.R1 %in% joined_regulons_sel$p2.R2] # overlap
sum(joined_regulons_sel$p1.R1 %in% joined_regulons_sel$p2.R2)/min(length(joined_regulons_sel$p2.R2), length(joined_regulons_sel$p1.R1))
# Showing all members of shared regulon 1
toString(names(cutree_out[cutree_out==1]))

# show network
regulon_network_sel = network(dist_matrix_sel, directed = FALSE)# , names.eval = "weights", matrix.type = "bipartite")
#regulon_network_sel = network(dist_matrix_sel, directed = FALSE, names.eval = "weights")
ggnet2(regulon_network_sel, node.size = 1, label=T, label.size = 5) # simple plotting
ggnet2(regulon_network_sel, node.size = 1, label=T, label.size = 5, edge.size = "weights") # simple plotting
# show network using cutoff
regulon_network_sel = network(dist_matrix_sel>.5, directed = FALSE)
ggnet2(regulon_network_sel, node.size = 1, label=T, label.size = 5) # simple plotting

# We can plot a graph
#plot(regulon_network_sel, vertex.cex = 3, mode = "circle")
test_edges = melt(dist_matrix_sel)
label_to_int=1:dim(dist_matrix_sel)[1]; names(label_to_int)=rownames(dist_matrix_sel)
colnames(test_edges) = c('from', 'to', 'weight')
nodes=data.frame(id=label_to_int, label=names(label_to_int))
test_edges$from=label_to_int[test_edges$from]
test_edges$to=label_to_int[test_edges$to]
#regulon_network_sel = network(test_edges[test_edges$weight>0,], matrix.type = "edgelist", vertex.attr = label_to_int, ignore.eval = FALSE)
regulon_network_sel = network(test_edges[test_edges$weight>0,], matrix.type = "edgelist", vertex.attr = nodes, ignore.eval = FALSE)
library(ggraph)
ggraph(regulon_network_sel, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8)+
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE, color='red') +
  #labs(edge_width = "Letters") +
  theme_graph()

cutree_out[cutree_out==1]

# The heatmap clearly shows a few groups
regulon_set_of_interest_sel_1 = names(cutree_out)[cutree_out==1] # c('p4.R4','p4.R2', 'p3.R1', 'p5.R1', 'p1.R2', 'p2.R3')
regulon_set_of_interest_sel_2 = names(cutree_out)[cutree_out==2] # c('p3.R4','p1.R3','p2.R4','p4.R3','p5.R4','p5.R3','p5.R5')
regulon_set_of_interest_sel_3 = names(cutree_out)[cutree_out==3] # c('p1.R1', 'p3.R3', 'p5.R2', 'p2.R2', 'p4.R1')
regulon_set_of_interest_sel_4 = names(cutree_out)[cutree_out==4] # c('p5.R6', 'p2.R5', 'p3.R5', 'p2.R1', 'p3.R2') # here the group is a bit less clear
# Gather data for export
regulon_data_for_export = list()
regulon_data_for_export$shared_regulon1 = as.data.frame(regulon_set_of_interest_sel_1)
regulon_data_for_export$shared_regulon2 = as.data.frame(regulon_set_of_interest_sel_2)
regulon_data_for_export$shared_regulon3 = as.data.frame(regulon_set_of_interest_sel_3)
regulon_data_for_export$shared_regulon4 = as.data.frame(regulon_set_of_interest_sel_4)

# Let's create lists of genes that belong to each of the shared groups
genes_set_of_interest_sel_1 = unique(unlist(lapply(regulon_set_of_interest_sel_1, function(X) {joined_regulons_sel[[X]]})))
genes_set_of_interest_sel_2 = unique(unlist(lapply(regulon_set_of_interest_sel_2, function(X) {joined_regulons_sel[[X]]})))
genes_set_of_interest_sel_3 = unique(unlist(lapply(regulon_set_of_interest_sel_3, function(X) {joined_regulons_sel[[X]]})))
genes_set_of_interest_sel_4 = unique(unlist(lapply(regulon_set_of_interest_sel_4, function(X) {joined_regulons_sel[[X]]})))
# gather for export
for (elem in names(joined_regulons_sel)) {
  regulon_data_for_export[[elem]] = as.data.frame(joined_regulons_sel[[elem]])
  }
regulon_data_for_export$shared_regulon1_genes = as.data.frame(genes_set_of_interest_sel_1)
regulon_data_for_export$shared_regulon2_genes = as.data.frame(genes_set_of_interest_sel_2)
regulon_data_for_export$shared_regulon3_genes = as.data.frame(genes_set_of_interest_sel_3)
regulon_data_for_export$shared_regulon4_genes = as.data.frame(genes_set_of_interest_sel_4)

# save these gene sets to input files for Homer
write.table(as.data.frame(genes_set_of_interest_sel_1), file = paste0(home_outputDir,'sel_set1.txt'),quote = F,col.names = F, row.names = F)
write.table(as.data.frame(genes_set_of_interest_sel_2), file = paste0(home_outputDir,'sel_set2.txt'),quote = F,col.names = F, row.names = F)
write.table(as.data.frame(genes_set_of_interest_sel_3), file = paste0(home_outputDir,'sel_set3.txt'),quote = F,col.names = F, row.names = F)
write.table(as.data.frame(genes_set_of_interest_sel_4), file = paste0(home_outputDir,'sel_set4.txt'),quote = F,col.names = F, row.names = F)
# Also write the background
write.table(as.data.frame(genes_detected_all_patients_mincelltreshold), file = paste0(home_outputDir,'background_mincelltreshold.txt'),quote = F,col.names = F, row.names = F)

toString(genes_set_of_interest_sel_1)
toString(genes_set_of_interest_sel_2)
toString(genes_set_of_interest_sel_3)
toString(genes_set_of_interest_sel_4)

# note that there is 95% overlap between genes in in genes_set_of_interest_sel_3 and genes from genes_from_similar_regulons_set1
# (the latter was calculated earlier in this script)
sum(genes_set_of_interest_sel_3 %in% genes_from_similar_regulons_set1) / min(length(genes_from_similar_regulons_set1), length(genes_set_of_interest_sel_3))

# Then run Homer analysis in command prompt
"
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_set1.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_1 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_set2.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_2 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_set3.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_3 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_set4.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_4 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
"

# Now check whether resulting transcription factors are actually in the dataset
# (use \(.* as grepl search and replace to quickly get gene names from homer output)
all_gene_symbols = rownames(groupedSCS$patientAllMod@ndata)

genes1 = c('Mef2b','Mef2d','CArG','Mef2c','PRDM14','TBP','Mef2a')
genes2 = c('SpiB', 'PU.1', 'ETS1', 'Etv2', 'IRF8', 'ERG', 'GABPA', 'ETS', 'ELF3', 'EWS:FLI1-fusion', 'Fli1', 'ETV1', 'ELF1', 'Elk4', 'Elk1', 'NFY', 'EWS:ERG-fusion', 'Ets1-distal', 'EHF')
genes3 = c('GATA:SCL', 'SCL')
genes4 = c('JunB', 'Atf3', 'Fra2', 'Fra1', 'BATF', 'Fosl2', 'AP-1', 'ZNF264', 'Jun-AP1', 'Usf2')

for (gene_set in list(genes1,genes2,genes3,genes4)){
  print('========')
  for (gene in gene_set) {
   
      if (  any(grepl(paste0('^',toupper(gene),'_'), all_gene_symbols))  ) {
          print(paste0('X ',toupper(gene),' (found)'))
      } else {
          print(paste0('  ',toupper(gene),' (not found)'))
      }
         
  }
}


# Another "upset" matrix, to see gene memberships
regulon_set_of_interest_sel_list = 
  list(regulon_set_of_interest_sel_1,regulon_set_of_interest_sel_2,regulon_set_of_interest_sel_3,regulon_set_of_interest_sel_4)
gene_set_of_interest_sel_list = 
  list(genes_set_of_interest_sel_1,genes_set_of_interest_sel_2,genes_set_of_interest_sel_3,genes_set_of_interest_sel_4)
shared_regulon_core_list=list()
for (idx in 1:length(regulon_set_of_interest_sel_list)) {
  
  # get current list with all genes connected to this shared regulon
  current_list = regulon_set_of_interest_sel_list[[idx]]
  # determine for each of the regulons whether each of the gene is a member
  memberships1 = lapply(current_list, function(X) {gene_set_of_interest_sel_list[[idx]] %in% joined_regulons_sel[[X]]})
  # convert outcome to dataframe
  memberships1_df = as.data.frame(do.call(cbind, memberships1))
  
  # set row/colnames
  rownames(memberships1_df) = gene_set_of_interest_sel_list[[idx]]
  colnames(memberships1_df) = current_list
  
  # now take only those genes that are present in at least 3 patients
  shared_regulon_core_list[[idx]]=gene_set_of_interest_sel_list[[idx]][apply(memberships1_df, 1, sum)>=3]
    # note that in principle we're looking for genes at least present in X patients, 
    # however, since two regulons from the same patient cannot contain the same genes
    # we can use this method of counting despite the fact that there are multiple
    # regulons from 1 patient present sometimes
  
  # print for user
  print(paste0("For shared regulon ", idx))
  print(toString(shared_regulon_core_list[[idx]]))
}

# export to xlsx
names(shared_regulon_core_list) = paste0('SharedRegulon',1:length(regulon_set_of_interest_sel_list))
write.xlsx(shared_regulon_core_list, 'MW_custom_plots/SharedRegulons.xlsx')

# more elaborate xlsx 
for (elem in names(shared_regulon_core_list)) {
  regulon_data_for_export[[paste0(elem,'core')]] = as.data.frame(shared_regulon_core_list[[elem]])
}
write.xlsx(regulon_data_for_export, 'MW_custom_plots/SharedRegulons_extended.xlsx')

# another export for homer
for (idx in 1:4) {
  write.table(as.data.frame(shared_regulon_core_list[[idx]]), file = paste0(home_outputDir,'sel_core_set',idx,'.txt'),quote = F,col.names = F, row.names = F)
}
# Then run Homer analysis in command prompt
"
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_core_set1.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_core_1 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_core_set2.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_core_2 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_core_set3.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_core_3 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/sel_core_set4.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_sel_core_4 -bg /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/background_mincelltreshold.txt
"

##############################

# Now look at GO terms for these shared regulons; 
# Take the "core" genes, which are associated with this regulon in only
# three patients.
# The query genes are in: 
# ---> shared_regulon_core_list
# The background is a bit harder to determine; we want to take along only genes that
# are detected in the patients because we're not interested in finding out 
# genes that we sequenced are more likely to be heart-related genes.
# We also want to have an equal selection on expression; so if we for the 
# regulon matrix only used genes with expression >20% cells (which we did)
# I think that's also the criterion for taking along genes now.
# I don't think additional criteria are necessary (you could e.g. also
# take as background the genes that are considered in the correlatino matrix,
# but those might already be a gene set with special properties, and we might 
# be interested in those properties already)
# So then are background set would be: 
# ---> genes_detected_all_patients_mincelltreshold
# (determined in "HCM_SCS_Manuscript_c_comparing_gene_detection_patients.R")

# Code from previous version:
#regulonGO <- analyzeRegulonGeneOntology(config, groupedSCS, correlationMatrices, regulons, pathwayPCutoff=0.05, GOKegg='GO', includeChildTerms=F)
#plotRegulonGeneOntology(config, groupedSCS, regulonGO, regulons, topX=10, plotRegulatedGenes=FALSE, overrideDir='Fig4/C', outputMode = 'pdf')

# Remove gene conversion table if not all genes of interests are in it
# This will result in regeneration fo the table
# Note: for the genes detected, not all are in the conversion table
# even if regenerated
if (!all(genes_detected_all_patients_mincelltreshold %in% geneIdentifierConversionTable$external_gene_name)) {
  rm('geneIdentifierConversionTable')
}
# Now generate the GO terms
regulonGO_MW = list()
for (regulon_name in names(shared_regulon_core_list)) {
  regulonGO_MW[[regulon_name]] = analyzeGeneOntology_MW(config = config, 
    all_genes = rownames(groupedSCS$patientAll@ndata), 
    background_genes = genes_detected_all_patients_mincelltreshold, 
    genes_query = shared_regulon_core_list[[regulon_name]],
    pathwayPCutoff =0.05 ,GOKegg = 'GO',includeChildTerms = F)
}
# Then export to xls
openxlsx::write.xlsx(regulonGO_MW, file = paste0(getwd(), '/analysis_2020-05-04-18h_17m/FigX/regulon_GO_mw.xlsx'))
# Then create some plots
plotRegulonGeneOntology(config, groupedSCS, regulonGO_MW[[1]], 
  regulons, topX=10, plotRegulatedGenes=FALSE, overrideDir='FigX/', outputMode = 'pdf')
# Plot gene ontology
for (regulon_idx in 1:length(regulonGO_MW)) {
  shorthand_save(shorthand_GOplot(regulonGO_MW[[regulon_idx]]), 
    '/analysis_2020-05-04-18h_17m/FigX/', paste0('regulon_sharedregulon',regulon_idx,'_GO_terms.pdf'))
}


################################################################################
# Which genes are most characteristic per regulon?
# Let's find out using the silhouette scores

# silhouette scores were calculated per regulon, now do some juggling to get them in tables
patient_names_mod = c('patient1Mod','patient2Mod','patient3Mod','patient4Mod','patient5Mod')
silhouette_tbl=list()
average_silh=list()
median_silh=list()
for (reg_idx in 1:length(shared_regulon_core_list)) {

    reg_name=names(shared_regulon_core_list)[reg_idx]

    current_genenames_chrXX=find__chrXX(shared_regulon_core_list[[reg_name]], rownames(groupedSCS$patientAll@expdata))

    silhouette_tbl[[reg_idx]]=matrix(rep(NA,length(current_genenames_chrXX)*5),nrow=length(current_genenames_chrXX))
    rownames(silhouette_tbl[[reg_idx]]) = shared_regulon_core_list[[reg_name]]
    colnames(silhouette_tbl[[reg_idx]]) = patient_names_mod

    # put silhouette scores in table for this regulon
    for (pt_idx in 1:5) {
        silhouette_tbl[[reg_idx]][,pt_idx] = 
            regulons_selectedGenes[[patient_names_mod[pt_idx]]]$silhouetteData[current_genenames_chrXX,]$sil_width
    }
    
    # also order it
    average_silh[[reg_idx]] = apply(silhouette_tbl[[reg_idx]], 1, mean, na.rm=T)
    median_silh[[reg_idx]] = apply(silhouette_tbl[[reg_idx]], 1, median, na.rm=T)
    silhouette_tbl[[reg_idx]]=silhouette_tbl[[reg_idx]][order(average_silh[[reg_idx]], decreasing = F),]
    #median_silh = apply(silhouette_tbl[[1]], 1, median, na.rm=T)
    #silhouette_tbl_md=list()
    #silhouette_tbl_md[[1]]=silhouette_tbl[[1]][order(median_silh, decreasing = F),]
    
}

# Now plot
listp=list()
for (reg_idx in 1:4) {
    # create plot of silhouette scores
    listp[[reg_idx]]=ggplot(melt(silhouette_tbl[[reg_idx]]))+
        geom_boxplot(aes(x=Var1,y=value), color='gray')+
        geom_jitter(aes(x=Var1,y=value, color=Var1))+
        coord_flip()+
        scale_color_manual(values=rep(sample(col_Dark2),10))+ #col_vector_60
        theme_bw()+theme(legend.position = 'none')+ggtitle(paste0('Silhouette scores reg ', reg_idx))+
        give_better_textsize_plot(10)
    # show
    print(listp[[reg_idx]])
    # save
    shorthand_save(listp[[reg_idx]], '/analysis_2020-05-04-18h_17m/FigX/', 
        paste0('regulon_silhouette_sharedregulon',reg_idx,'.pdf'),
        mywidth = 10, myheight = length(shared_regulon_core_list[[reg_idx]])*.5+3)
}

# create slightly extended table:
extended_silhouette_tbl=list()
for (reg_idx in 1:4) {
  #print(paste0(reg_idx))
  extended_silhouette_tbl[[reg_idx]]=cbind(data.frame(gene_name=rownames(silhouette_tbl[[reg_idx]])),as.data.frame(silhouette_tbl[[reg_idx]]))
  extended_silhouette_tbl[[reg_idx]]$average_silh = average_silh[[reg_idx]][rownames(silhouette_tbl[[reg_idx]])]
  extended_silhouette_tbl[[reg_idx]]$median_silh = median_silh[[reg_idx]][rownames(silhouette_tbl[[reg_idx]])]
  extended_silhouette_tbl[[reg_idx]] = 
    extended_silhouette_tbl[[reg_idx]][order(extended_silhouette_tbl[[reg_idx]]$average_silh, decreasing = T),]
}
names(extended_silhouette_tbl) = names(shared_regulon_core_list)

# save to excel
openxlsx::write.xlsx(extended_silhouette_tbl, 
    file = paste0(getwd(), '/analysis_2020-05-04-18h_17m/FigX/regulon_silhouette.xlsx'))

################################################################################
# I wonder whether TEAD4 was MALAT1-related?
# --> don't think so (see below); but MALAT does have some interesting genes it's correlated with

library(mutoss)
gene_counts=apply(groupedSCS$patientAll@ndata>.1,1,sum)
volcano_TEAD4_df = get_volcano_df(expression_matrix = as.matrix(groupedSCS$patientAll@ndata[gene_counts>3,]), my_gene = 'TEAD4')
volcano_MALAT1_df = get_volcano_df(expression_matrix = as.matrix(groupedSCS$patientAll@ndata[gene_counts>3,]), my_gene = 'MALAT1')
volcano_KCNQ1OT1_df = get_volcano_df(expression_matrix = as.matrix(groupedSCS$patientAll@ndata[gene_counts>3,]), my_gene = 'KCNQ1OT1')

plot_volcano(volcano_TEAD4_df)
plot_volcano(volcano_MALAT1_df)
plot_volcano(volcano_KCNQ1OT1_df)

################################################################################
# Also, is there a correlation between FSC-A and regulon 1?
# Yes, there is, but it's not so visible on tSNE maps

# Plot regulon 1 on tSNE (patients 4+5 only)
shorthand_composite_expression(current_SCS = groupedSCS$patientAllWithIndexMod, current_gene_set = genes_set_of_interest_sel_1)
# Plot regulon 1 on tSNE (all patients)
shorthand_composite_expression(current_SCS = groupedSCS$patientAllMod, current_gene_set = genes_set_of_interest_sel_1)
p=shorthand_composite_expression(current_SCS = groupedSCS$patientAllMod, current_gene_set = shared_regulon_core_list[[1]])
pp = p + shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
pp; shorthand_save(pp, thedir = 'MW_custom_plots/', filename = 'tSNE_f_all_shared-reg-1.png', mywidth = 7.5, myheight = 7.5)

# Get correlation with FSC-A
regulon1_composite_expression = shorthand_composite_expression_values(current_SCS = groupedSCS$patientAllWithIndexMod, current_gene_set = genes_set_of_interest_sel_1)
regulon2_composite_expression = shorthand_composite_expression_values(current_SCS = groupedSCS$patientAllWithIndexMod, current_gene_set = genes_set_of_interest_sel_2)
regulon3_composite_expression = shorthand_composite_expression_values(current_SCS = groupedSCS$patientAllWithIndexMod, current_gene_set = genes_set_of_interest_sel_3)
regulon4_composite_expression = shorthand_composite_expression_values(current_SCS = groupedSCS$patientAllWithIndexMod, current_gene_set = genes_set_of_interest_sel_4)

cellnames.p4 = colnames(groupedSCS$patient4Mod@ndata)
cellnames.p5 = colnames(groupedSCS$patient5Mod@ndata)

ggplot(data.frame(regulon1_composite=regulon1_composite_expression[cellnames.p4], FSC_A=indexData_pt4pt5[cellnames.p4,]$FSC_A))+
    geom_point(aes(x=regulon1_composite, y=FSC_A))

toplot_df = 
    rbind(
    data.frame(regulon_composite=regulon1_composite_expression[cellnames.p4], FSC_A=indexData_pt4pt5[cellnames.p4,]$FSC_A, regulon='R1', patient='p4'),
    data.frame(regulon_composite=regulon2_composite_expression[cellnames.p4], FSC_A=indexData_pt4pt5[cellnames.p4,]$FSC_A, regulon='R2', patient='p4'),
    data.frame(regulon_composite=regulon3_composite_expression[cellnames.p4], FSC_A=indexData_pt4pt5[cellnames.p4,]$FSC_A, regulon='R3', patient='p4'),
    data.frame(regulon_composite=regulon4_composite_expression[cellnames.p4], FSC_A=indexData_pt4pt5[cellnames.p4,]$FSC_A, regulon='R4', patient='p4'),
    data.frame(regulon_composite=regulon1_composite_expression[cellnames.p5], FSC_A=indexData_pt4pt5[cellnames.p5,]$FSC_A, regulon='R1', patient='p5'),
    data.frame(regulon_composite=regulon2_composite_expression[cellnames.p5], FSC_A=indexData_pt4pt5[cellnames.p5,]$FSC_A, regulon='R2', patient='p5'),
    data.frame(regulon_composite=regulon3_composite_expression[cellnames.p5], FSC_A=indexData_pt4pt5[cellnames.p5,]$FSC_A, regulon='R3', patient='p5'),
    data.frame(regulon_composite=regulon4_composite_expression[cellnames.p5], FSC_A=indexData_pt4pt5[cellnames.p5,]$FSC_A, regulon='R4', patient='p5')
    )
ggplot(toplot_df,aes(x=regulon_composite, y=FSC_A))+
    geom_point()+
    geom_smooth(method='lm')+
    facet_grid(regulon ~ patient)+theme_bw()
ggplot(toplot_df[toplot_df$regulon=='R1',],aes(x=regulon_composite, y=FSC_A))+
    geom_point()+
    geom_smooth(method='lm')+
    facet_grid(regulon ~ patient)+theme_bw()+give_better_textsize_plot(12)

# Calculate correlations (perhaps this is a good parameter for the overall claim!)
cor.test(regulon1_composite_expression[cellnames],indexData_pt4pt5[cellnames,]$FSC_A)
cor.test(regulon2_composite_expression[cellnames],indexData_pt4pt5[cellnames,]$FSC_A)
cor.test(regulon3_composite_expression[cellnames],indexData_pt4pt5[cellnames,]$FSC_A)
cor.test(regulon4_composite_expression[cellnames],indexData_pt4pt5[cellnames,]$FSC_A)


# and with outlier removal
toplot_df.p4 = toplot_df[toplot_df$regulon=='R1'&toplot_df$patient=='p4',]
toplot_df.p5 = toplot_df[toplot_df$regulon=='R1'&toplot_df$patient=='p5',]
mylims.p4 = calc_limits(toplot_df.p4$FSC_A)
mylims.p5 = calc_limits(toplot_df.p5$FSC_A)
toplot_df.p4.sel=toplot_df.p4[toplot_df.p4$FSC_A>mylims.p4[1]&toplot_df.p4$FSC_A<mylims.p4[2],]
toplot_df.p5.sel=toplot_df.p5[toplot_df.p5$FSC_A>mylims.p5[1]&toplot_df.p5$FSC_A<mylims.p5[2],]
ggplot(toplot_df.p4,aes(x=regulon_composite,y=FSC_A))+
        geom_point()+
        geom_point(data=toplot_df.p4.sel, color='red')+
        geom_smooth(method='lm')

ggplot(toplot_df.p5,aes(x=regulon_composite,y=FSC_A))+
        geom_point()+
        geom_point(data=toplot_df.p5.sel, color='red')+
        geom_smooth(method='lm')    

# Also test whether the correlation is not determined by outliers
cor.test(toplot_df.p4$regulon_composite, toplot_df.p4$FSC_A)
cor.test(toplot_df.p5$regulon_composite, toplot_df.p5$FSC_A)
        
cor.test(toplot_df.p4.sel$regulon_composite, toplot_df.p4.sel$FSC_A)
cor.test(toplot_df.p5.sel$regulon_composite, toplot_df.p5.sel$FSC_A)

# Perhaps further illustrate by showing FSC-A for categories of low, medium and high reg expr
ggplot(toplot_df.p5)+
  geom_histogram(aes(x=FSC_A))
ggplot(toplot_df.p5)+
  geom_histogram(aes(x=regulon_composite))
toplot_df.p5$cat_hlm = NA
lims=calc_limits(toplot_df.p5$regulon_composite, percentile = .33)
toplot_df.p5$cat_hlm[toplot_df.p5$regulon_composite>=lims[2]] = 'high'
toplot_df.p5$cat_hlm[toplot_df.p5$regulon_composite<=lims[1]] = 'low'
toplot_df.p5$cat_hlm[toplot_df.p5$regulon_composite<lims[2]&toplot_df.p5$regulon_composite>lims[1]] = 'medium'
toplot_df.p5$cat_hlm = factor(toplot_df.p5$cat_hlm, levels=c('low','medium','high'))
ggplot(toplot_df.p5,aes(x=cat_hlm,y=FSC_A))+
        geom_boxplot()+
        geom_jitter()
p=ggplot(toplot_df.p5,aes(x=cat_hlm,y=FSC_A))+
        geom_violin()+xlab('Regulon 1 expression')+theme_minimal()+give_better_textsize_plot(12)
p; shorthand_save(p, 'MW_custom_plots/', 'FSCA_violins_regulon1_3cats.png', mywidth = 7.5, myheight = 7.5)
#ggplot(toplot_df.p5,aes(x=cat_hlm,y=regulon_composite))+
#        geom_boxplot()+
#        geom_jitter()
ggplot(toplot_df.p5,aes(x=cat_hlm,y=FSC_A))+
        geom_bar(aes(x=cat_hlm,y=mean(FSC_A, na.rm=T)), stat='identity')


