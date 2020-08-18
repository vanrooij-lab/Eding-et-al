
######################################################################
# Please see notes in "HCM SCS Manuscript_AE.R" first

################################################################################
# RE-ANALYSIS OF FSC-A DATA INCLUDING BOTH PATIENTS
#
# Some things I did here:
# - Analyze FSCA-gene correlations per patient
# - Also pooled patient FSCA-gene correlations, but wouldn't recommend that
# - Scatter plots & corrs with outlier removal
# - Combining patient data by normalization beforehand, then determine corrs
# - Determining corrs from randomized data as reference, to see how high corrs you get from random data

################################################################################
current_outputdir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/analysis_2020-03/fsca-custom/'

# Analysis by RaceID2 analysis for patients 4+5 combined should be available
# in "groupedSCS$patientAllWithIndexMod"

################################################################################

#----

# First construct a complete indexData parameter
# Merge into one dataframe
indexData_MW = rbind(plateAL1, plateAL2, plateJE10, plateJE11)
# Fix column names
colnames(indexData_MW) <- c('Cell','Events','Parent','FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')
# Convert measurement columns to numeric
for (column_name in c('FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')) {
    indexData_MW[[column_name]] = as.numeric(unlist(lapply(indexData_MW[[column_name]], function(x) {str_replace(x, ',','.')})))    
}

# Create patient annotation, based on combined RaceID2 count table
# Create empty param with size of RaceID2 table
patient_source_pt4.pt5=rep(NA,dim((groupedSCS$patientAllWithIndexMod@ndata))[2])
# Assign column (cell) names to the annotation for later reference
names(patient_source_pt4.pt5)=colnames(groupedSCS$patientAllWithIndexMod@ndata)
for (current_name in names(config$samples)) {
    # go over patients and look for their identifiers in the column names
    patient_source_pt4.pt5[grepl(paste(config$samples[[current_name]],collapse = '|'),
                        colnames(groupedSCS$patientAllWithIndexMod@ndata))] = current_name
}
# NOTE: I MADE THE ABOVE LINES IN A FUNCTION ALSO: "shorthand_give_p_source()"

# Let's create convenient parameters to handle the indexdata
pt4.pt5.sel.idxs = patient_source_pt4.pt5 %in% c('patient4','patient5')
pt4.pt5.cellnames = colnames(groupedSCS$patientAllWithIndexMod@ndata)[pt4.pt5.sel.idxs]
rownames(indexData_MW) = indexData_MW$Cell
indexData_pt4 = indexData_MW[colnames(groupedSCS$patient4Mod@ndata),]
indexData_pt5 = indexData_MW[colnames(groupedSCS$patient5Mod@ndata),]
indexData_pt4pt5 = indexData_MW[colnames(groupedSCS$patientAllWithIndexMod@ndata),]
expressionData_FSCA = groupedSCS$patientAllWithIndexMod@ndata
FSCA_matrix = matrix(indexData_pt4pt5$FSC_A,nrow=1)
rownames(FSCA_matrix) = c('FSCA')
colnames(FSCA_matrix) = colnames(expressionData_FSCA)
expressionData_FSCA = rbind(expressionData_FSCA, as.data.frame(FSCA_matrix))
    # sanity check (outcome should be false)
    any(colnames(groupedSCS$patientAllWithIndexMod@ndata) != indexData_pt4pt5$Cell)
# also add patient annotation
indexData_pt4pt5$Patient = patient_source_pt4.pt5

# just a little test, correlate gene 1 with FSC-A
thecor = cor.test(as.matrix(groupedSCS$patientAllWithIndexMod@ndata)[1,], indexData_pt4pt5$FSC_A)

# create some selection criteria
gene_in_cells = apply(groupedSCS$patientAllWithIndexMod@ndata>0.1,1,sum)
gene_in_cells_pt4 = apply(groupedSCS$patientAllWithIndexMod@ndata[,patient_source_pt4.pt5=='patient4']>0.1,1,sum)
gene_in_cells_pt5 = apply(groupedSCS$patientAllWithIndexMod@ndata[,patient_source_pt4.pt5=='patient5']>0.1,1,sum)
gene_in_cells_f = apply(groupedSCS$patientAllWithIndexMod@fdata>0.1,1,sum)
COUNT_CUTOFF_perc = .33
COUNT_CUTOFF.p45=COUNT_CUTOFF_perc*length(patient_source_pt4.pt5)
COUNT_CUTOFF.p4=COUNT_CUTOFF_perc*sum(patient_source_pt4.pt5=='patient4')
COUNT_CUTOFF.p5=COUNT_CUTOFF_perc*sum(patient_source_pt4.pt5=='patient5')

# Determine correlations, first for pooled patients, then for patient4, then patient5

# Corrs Pooled
# using ndata (groupedSCS$patientAllWithIndexMod@ndata)
cor.test.all.out = as.data.frame(t(apply(as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells>COUNT_CUTOFF.p45,]),1,function(X) {thecor=cor.test(X, indexData_pt4pt5$FSC_A); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.all.out$pval.adj = p.adjust(cor.test.all.out$p.value, method="BH")

# Corrs P4 (this is where the actual corrs are calculated)
cor.test.pt4.out = as.data.frame(t(apply(as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt4>COUNT_CUTOFF.p4,patient_source_pt4.pt5=='patient4']),1,function(X) {thecor=cor.test(X, indexData_pt4pt5$FSC_A[patient_source_pt4.pt5=='patient4']); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.pt4.out$pval.adj = p.adjust(cor.test.pt4.out$p.value, method="BH")

# Corrs P5
cor.test.pt5.out = as.data.frame(t(apply(as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt5>COUNT_CUTOFF.p5,patient_source_pt4.pt5=='patient5']),1,function(X) {thecor=cor.test(X, indexData_pt4pt5$FSC_A[patient_source_pt4.pt5=='patient5']); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.pt5.out$pval.adj = p.adjust(cor.test.pt5.out$p.value, method="BH")    

# Now let's check which genes are significantly correlated
# Spoiler: almost none, ie only 1
cor.test.pt4.out[cor.test.pt4.out$pval.adj<0.05,]
cor.test.pt5.out[cor.test.pt5.out$pval.adj<0.05,]

# top 10 correlated genes for the two patients
cor.test.pt4.out[order(cor.test.pt4.out$pval.adj),][1:10,]$gene_name_short
cor.test.pt5.out[order(cor.test.pt5.out$pval.adj),][1:10,]$gene_name_short

# Can compare in Venn-diagrams top 30 of correlated genes for both patients
venn1=cor.test.pt4.out[order(cor.test.pt4.out$corr, decreasing = T),][1:30,]$gene_name_short
venn2=cor.test.pt5.out[order(cor.test.pt5.out$corr, decreasing = T),][1:30,]$gene_name_short

venn1[venn1 %in% venn2]

library("limma")
venn_genes = union(venn1, venn2)
gene_overlap_df=as.matrix(data.frame(pt4=venn_genes %in% venn1, pt5=venn_genes %in% venn2))
my_ven <- vennCounts(gene_overlap_df)
vennDiagram(my_ven, cex=1) 

# Create volcano plots
# for pat 4
cor.test.pt4.out$gene_name = row.names(cor.test.pt4.out)
cor.test.pt4.out$gene_name_short = strip__chrXX(row.names(cor.test.pt4.out))
cor.test.pt4.out$corr=cor.test.pt4.out$cor.cor
p.p4=plot_volcano(cor.test.pt4.out, manual_gene_name = 'FSC-A (pt.4)', NRLABELED = 10, mypvaltreshold = 0.05)
p.p4
# for pat 5
cor.test.pt5.out$gene_name = row.names(cor.test.pt5.out)
cor.test.pt5.out$gene_name_short = strip__chrXX(row.names(cor.test.pt5.out))
cor.test.pt5.out$corr=cor.test.pt5.out$cor.cor
p.p5=plot_volcano(cor.test.pt5.out, manual_gene_name = 'FSC-A (pt.5)', NRLABELED = 10, mypvaltreshold = 0.05)
p.p5
# Combined plot
(p.p4+ggtitle('Pt.4')+ylim(c(0,3.5)))+(p.p5+ggtitle('Pt.5')+ylim(c(0,3.5)))+plot_annotation(title='Gene correlations with FSC-A')

# Volcano for the combined patients (should be used! pooling patients leads to artefacts!)
cor.test.all.out$gene_name = row.names(cor.test.all.out)
cor.test.all.out$gene_name_short = strip__chrXX(row.names(cor.test.all.out))
cor.test.all.out$corr=cor.test.all.out$cor.cor
p.p45=plot_volcano(cor.test.all.out, manual_gene_name = 'FSC-A (pt.4&pt.5)', NRLABELED = 10, mypvaltreshold = 0.05)
p.p45

########################################################################

# Let's create a scatter plots of the correlation coefficients
overlapping_genes = cor.test.pt4.out$gene_name[cor.test.pt4.out$gene_name %in% cor.test.pt5.out$gene_name]
correlations_pt4.pt5 = 
    data.frame(cor.pt4=cor.test.pt4.out[overlapping_genes,]$cor.cor, cor.pt5=cor.test.pt5.out[overlapping_genes,]$cor.cor,
        pval.pt4=cor.test.pt4.out[overlapping_genes,]$pval.adj, pval.pt5=cor.test.pt5.out[overlapping_genes,]$pval.adj,
        gene=strip__chrXX(overlapping_genes))
correlations_pt4.pt5$distance = sqrt((correlations_pt4.pt5$cor.pt4-correlations_pt4.pt5$cor.pt5)^2)

# create a fit line
the_fit = lm(cor.pt4~cor.pt5, correlations_pt4.pt5)

#df_toplot = rbind(correlations_pt4.pt5, c(-.25, rep(NA,5)))
#df_toplot = rbind(df_toplot, c(.25, rep(NA,5)))
p=ggplot(correlations_pt4.pt5,aes(x=cor.pt4,y=cor.pt5))+
    geom_hline(yintercept = 0, colour='black')+geom_vline(xintercept = 0, colour='black')+ #x=0, y=0 axes lines
    geom_point(shape=1, color='grey')+
    geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,],color='red',shape=1)+
    geom_abline(intercept = the_fit$coefficients[1], slope = the_fit$coefficients[2])+
    geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,],color='red',
        aes(label=gene), cex=3)+
    #geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.25 | correlations_pt4.pt5$cor.pt5>.2,],color='blue',
    #    aes(label=gene), cex=4)+
    theme_bw()+give_better_textsize_plot(12)+xlab('Gene-cellsize correlation patient 4')+ylab('Gene-cellsize correlation patient 5')
    #(method='lm', color='black', se=F, size=.5)+
p

# let's also print 'red cloud' genes to user as string
toString(correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,]$gene)
# let's also export to a little text file for homer analysis
write.table(as.data.frame(correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,]$gene), file = paste0(home_outputDir,'set_fsca_corr_high.txt'),quote = F,col.names = F, row.names = F)
# use this command in terminal:
"findMotifs.pl /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/set_fsca_corr_high.txt human /Users/m.wehrens/Data/_2019_02_HCM_SCS/Homer_analyses/output_fsca_corr_high"
# export the plots
p; shorthand_save(p, 'MW_custom_plots/', 'fsca_corr_scatter_redcloud.png')
shorthand_save(p, 'MW_custom_plots/', 'fsca_corr_scatter_redcloud_small.png',mywidth = 10, myheight = 7.5)
# and also a slightly modified version
p=p+shorthand_tsne_joeptheme()+give_better_textsize_plot(10)#+coord_cartesian(xlim=c(-.25, .25), ylim=c(-.25, .25))
p; shorthand_save(p, 'analysis_2020-05-04-18h_17m/FigX/', 'fsca_corr_scatter_redcloud_submitstyle.pdf',mywidth = 10, myheight = 10)

ggplot(correlations_pt4.pt5,aes(x=cor.pt4,y=cor.pt5))+
    geom_point()+
    geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,],color='red')+
    geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.2&correlations_pt4.pt5$cor.pt5<=.1 | 
            correlations_pt4.pt5$cor.pt5>.2&correlations_pt4.pt5$cor.pt4<=.1,],color='blue')+
    geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.2&correlations_pt4.pt5$cor.pt5<=.1 | 
            correlations_pt4.pt5$cor.pt5>.2&correlations_pt4.pt5$cor.pt4<=.1,],color='blue',
        aes(label=gene), cex=4)+
    theme_bw()+give_better_textsize_plot(14)
# Plot with regulon 1 genes highlighted, if available
# (You need to have determined the regulons before doing this with the regulon script)
if (exists('genes_set_of_interest_sel_1')){
    ggplot(correlations_pt4.pt5,aes(x=cor.pt4,y=cor.pt5))+
        geom_point(shape=1)+
        geom_density_2d(color='black', bins=5)+
        geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% genes_set_of_interest_sel_1,],color='red')+
        geom_density_2d(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% genes_set_of_interest_sel_1,],color='red', bins=5)+
        #geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% genes_set_of_interest_sel_1,],color='green',
        #    aes(label=gene), cex=4)+
        theme_bw()+give_better_textsize_plot(14)
}
# Same plot, but only with those shared regulon 1 genes seen in 3 patients
if (exists('shared_regulon_core_list')){
    
    p1=ggplot(correlations_pt4.pt5,aes(x=cor.pt4,y=cor.pt5))+
        geom_point(shape=1)+
        geom_density_2d(color='black', bins=5)+
        geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='red')+
        geom_density_2d(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='red', bins=5)+
        #geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='black',
        #    aes(label=gene), cex=4)+
        theme_bw()+give_better_textsize_plot(12)
    shorthand_save(p1, 'MW_custom_plots/', 'regulon1_cellsize_scatter.png', mywidth = 10, myheight = 7.5)
    
    p2=ggplot(correlations_pt4.pt5,aes(x=cor.pt4,y=cor.pt5))+
        geom_hline(yintercept = 0, colour='black')+geom_vline(xintercept = 0, colour='black')+ #x=0, y=0 axes lines
        geom_point(shape=1, color='grey')+
        geom_density_2d(color='grey', bins=5)+
        geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='red')+
        geom_density_2d(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='red', bins=5)+
        geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$gene %in% shared_regulon_core_list[[1]],],color='black',
            aes(label=gene), cex=3)+
        theme_bw()+give_better_textsize_plot(12)
    shorthand_save(p2, 'MW_custom_plots/', 'regulon1_cellsize_scatter_names.png', mywidth = 10, myheight = 7.5)
    
    # alt. style
    p2=p2+shorthand_joeptheme()+xlab('Gene-cellsize correlation patient 4')+ylab('Gene-cellsize correlation patient 5')+
        give_better_textsize_plot(10)
    p2; shorthand_save(p2, 'analysis_2020-05-04-18h_17m/FigX/', 'fsca_corr_scatter_regulon1_submitstyle.pdf',mywidth = 10, myheight = 10)

}

# correlations of correlations :)
cor.test(correlations_pt4.pt5$cor.pt4,correlations_pt4.pt5$cor.pt5)

# Let's create an additional scatter of the p-values
ggplot(correlations_pt4.pt5,aes(x=pval.pt4,y=pval.pt5))+
    geom_point()
    #geom_point(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1 & correlations_pt4.pt5$distance<.05,],color='red')+
    #geom_text_repel(data=correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1 & correlations_pt4.pt5$distance<.05,],color='red',
    #    aes(label=gene), cex=4)+theme_bw()+give_better_textsize_plot(14)

# Naive way of determining composite p-value
correlations_pt4.pt5$pval.p4p5.mult=correlations_pt4.pt5$pval.pt4*correlations_pt4.pt5$pval.pt5

# Export this to Excel
library(openxlsx)
corr_pt4pt5_sel_cols = c("gene", "cor.pt4", "cor.pt5", "pval.pt4", "pval.pt5")
correlations_pt4.pt5$mean_corr = apply(correlations_pt4.pt5[,c('cor.pt4','cor.pt5')], 1, mean)
#openxlsx::write.xlsx(correlations_pt4.pt5[order(correlations_pt4.pt5$pval.p4p5.mult),corr_pt4pt5_sel_cols], 
#    paste0('analysis_2020-05-04-18h_17m/FigX/','FSCA_correlations_two_patients.xlsx'))
openxlsx::write.xlsx(correlations_pt4.pt5[order(correlations_pt4.pt5$mean_corr, decreasing = T),corr_pt4pt5_sel_cols], 
    paste0('analysis_2020-05-04-18h_17m/FigX/','FSCA_correlations_two_patients.xlsx'))

################################################################################
################################################################################
# This part is a bit more unpolished
# I don't think using fdata is necessarily the way to go

# using fdata (groupedSCS$patient4Mod@fdata)
COUNT_CUTOFF2 = 50
cellnames_pt4 = colnames(groupedSCS$patient4Mod@fdata)
gene_in_cells_p4 = apply(groupedSCS$patient4Mod@fdata>0.1,1,sum)
cor.test.pt4.out.f = as.data.frame(t(apply(as.matrix(groupedSCS$patient4Mod@fdata[gene_in_cells_p4>COUNT_CUTOFF2,]),1,function(X) {thecor=cor.test(X, indexData_pt4pt5[cellnames_pt4,]$FSC_A); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.pt4.out.f$pval.adj = p.adjust(cor.test.pt4.out.f$p.value, method="BH")
cor.test.pt4.out.f$corr = cor.test.pt4.out.f$cor.cor
cor.test.pt4.out.f$gene_name_short = strip__chrXX(rownames(cor.test.pt4.out.f))
p.f=plot_volcano(cor.test.pt4.out.f, manual_gene_name = 'FSC-A (pt.4, fdata)', NRLABELED = 10, mypvaltreshold = 0.05)+give_better_textsize_plot(14)
p.f

# Now filter the data
boundaries_FSC_A = calc_limits(indexData_MW$FSC_A,percentile = .02)

# Show histogram of two patients combined
ggplot(data=indexData_MW, aes(x=FSC_A))+
    geom_histogram()+
    geom_vline(xintercept = boundaries_FSC_A[1])+
    geom_vline(xintercept = boundaries_FSC_A[2])+
    theme_bw()

# Histogram separated per patient
ggplot(data=indexData_pt4pt5, aes(x=FSC_A, color=as.factor(Patient)))+
    geom_freqpoly(stat='density')+
    geom_vline(xintercept = boundaries_FSC_A[1])+
    geom_vline(xintercept = boundaries_FSC_A[2])+
    theme_bw()+give_better_textsize_plot(12)

# calculate correlations from combined data
FSCA_volcano_df = get_volcano_df(expression_matrix = as.matrix(expressionData_FSCA), my_gene = 'FSCA')
plot_volcano(my_corrs_df_current)
View(expressionData_FSCA)

################################################################################
# Testing outlier removal

# -----
# Shorthand function simple scatter with outlier removal indicated
shorthand_scatter_fsca_gene_sel = function(gene_expression, FSC_A,mytitle,mylims=c(-Inf,Inf),gene_name) {
    df1=data.frame(gene_expression=gene_expression, FSC_A=FSC_A)
    df1_sel = df1[df1$FSC_A>mylims[1]&df1$FSC_A<mylims[2],]
    cor_out_sel = cor.test(df1_sel$FSC_A,df1_sel$gene_expression)
    ggplot(df1,aes(x=FSC_A,y=gene_expression))+
            geom_point()+
            geom_point(data=df1_sel, color='red')+
            geom_smooth(method='lm')+
            geom_smooth(data=df1_sel, method='lm', color='red')+
            theme_bw()+ggtitle(paste0(mytitle,' (cor = ',round(cor_out_sel$estimate,2),')'))+
            ylab(gene_name)
}
shorthand_scatter_fsca_gene_ptcol = function(gene_expression, FSC_A,mytitle,gene_name,patient_annotation) {
    df1=data.frame(gene_expression=gene_expression, FSC_A=FSC_A, patient=patient_annotation)
    cor_out_sel = cor.test(df1_sel$FSC_A,df1_sel$gene_expression)
    ggplot(df1,aes(x=FSC_A,y=gene_expression))+
            geom_point(aes(color=patient))+
            geom_smooth(method='lm', color='black')+
            theme_bw()+ggtitle(paste0(mytitle,' (cor = ',round(cor_out_sel$estimate,2),')'))+
            ylab(gene_name)
}
# -----
# now for patient 4, take out outliers
current_indexdata=indexData_pt4pt5$FSC_A[patient_source_pt4.pt5=='patient4']
names(current_indexdata) = indexData_pt4pt5$Cell[patient_source_pt4.pt5=='patient4']
current_lims=calc_limits(current_indexdata,percentile = .02)
current_indexdata_sel = current_indexdata[current_indexdata>current_lims[1]&
        current_indexdata<current_lims[2]]
cor.test.pt4.sel.out = as.data.frame(t(apply(as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt4>COUNT_CUTOFF.p4,names(current_indexdata_sel)]),1,function(X) {thecor=cor.test(X, current_indexdata_sel); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.pt4.sel.out$pval.adj = p.adjust(cor.test.pt4.sel.out$p.value, method="BH")
cor.test.pt4.sel.out$corr = cor.test.pt4.sel.out$cor.cor
cor.test.pt4.sel.out$gene_name = rownames(cor.test.pt4.sel.out)
cor.test.pt4.sel.out$gene_name_short = strip__chrXX(cor.test.pt4.sel.out$gene_name)
# show vulcano
p4v=plot_volcano(cor.test.pt4.sel.out,mypvaltreshold = .05, NRLABELED = 10, manual_gene_name = 'FSC-A (pt.4)')
# show scatter 
p4s=shorthand_scatter_fsca_gene_sel(gene_expression = as.vector(get_expression_gene2(expression_matrix = as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt4>COUNT_CUTOFF.p4,names(current_indexdata)]),
                                    gene_query = 'ACTA1')),
                                FSC_A  = current_indexdata, mytitle = 'pt.4',mylims=current_lims, gene_name = 'ACTA1')
# repeat for patient 5
current_indexdata=indexData_pt4pt5$FSC_A[patient_source_pt4.pt5=='patient5'&!is.na(indexData_pt4pt5$FSC_A)]
names(current_indexdata) = indexData_pt4pt5$Cell[patient_source_pt4.pt5=='patient5'&!is.na(indexData_pt4pt5$FSC_A)]
current_lims=calc_limits(current_indexdata,percentile = .02)
current_indexdata_sel = current_indexdata[current_indexdata>current_lims[1]&current_indexdata<current_lims[2]]
cor.test.pt5.sel.out = as.data.frame(t(apply(as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt5>COUNT_CUTOFF.p5,names(current_indexdata_sel)]),1,function(X) {thecor=cor.test(X, current_indexdata_sel); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
cor.test.pt5.sel.out$pval.adj = p.adjust(cor.test.pt5.sel.out$p.value, method="BH")
cor.test.pt5.sel.out$corr = cor.test.pt5.sel.out$cor.cor
cor.test.pt5.sel.out$gene_name = rownames(cor.test.pt5.sel.out)
cor.test.pt5.sel.out$gene_name_short = strip__chrXX(cor.test.pt5.sel.out$gene_name)
# show vulcano
p5v=plot_volcano(cor.test.pt5.sel.out,mypvaltreshold = .05, NRLABELED = 10, manual_gene_name = 'FSC-A (pt.5)')
# show scatter 
p5s=shorthand_scatter_fsca_gene_sel(gene_expression = as.vector(get_expression_gene2(expression_matrix = as.matrix(groupedSCS$patientAllWithIndexMod@ndata[gene_in_cells_pt5>COUNT_CUTOFF.p5,names(current_indexdata)]),
                                    gene_query = 'ACTA1')),
                                FSC_A  = current_indexdata, mytitle = 'pt.5',mylims=current_lims,
                                gene_name = 'ACTA1')
(p4s+p5s)/(p4v+p5v)

#----
# Now also compare original correlations with outlier-removed correlations
gene_sel = rownames(cor.test.pt4.sel.out)
df_compare_sel = data.frame(cor.sel = cor.test.pt4.sel.out[gene_sel,]$corr,
                            cor.all = cor.test.pt4.out[gene_sel,]$corr,
                            gene_name = cor.test.pt4.out[gene_sel,]$gene_name,
                            gene_name_short = cor.test.pt4.out[gene_sel,]$gene_name_short)
ggplot(df_compare_sel, aes(x=cor.sel, y=cor.all))+
    geom_point()+
    geom_point(data=df_compare_sel[df_compare_sel$gene_name_short=='ACTA1',],color='red')+
    theme_bw()+ggtitle('Pt 4')
# Same for patient 5
gene_sel = rownames(cor.test.pt5.sel.out)
df_compare_sel_p5 = data.frame(cor.sel = cor.test.pt5.sel.out[gene_sel,]$corr,
                            cor.all = cor.test.pt5.out[gene_sel,]$corr,
                            gene_name = cor.test.pt5.out[gene_sel,]$gene_name,
                            gene_name_short = cor.test.pt5.out[gene_sel,]$gene_name_short)
ggplot(df_compare_sel_p5, aes(x=cor.sel, y=cor.all))+
    geom_point()+
    geom_point(data=df_compare_sel_p5[df_compare_sel_p5$gene_name_short=='ACTA1',],color='red')+
    theme_bw()+ggtitle('Pt 5')

################################################################################
# Some more less polished stuff
# -----

expression_data_df = as.data.frame(t(expressionData_FSCA))
expression_data_df$patient = patient_source_pt4.pt5

ggplot(expression_data_df)+
    geom_point(aes(x=FSCA, y=UQCR10__chr22))

library(patchwork)
scatter_pt4pt5 = function(expression_data_df,feature1, feature2) {
    
    p1=ggplot(expression_data_df, aes_string(x=feature1, y=feature2))+
        geom_point()+
        geom_point(aes_string(color='patient'))+
        geom_smooth(method='lm',aes_string(color='patient'))+
        geom_smooth(method='lm',color='black')+
        theme_bw()+ggtitle('Combined')
    p1
    p2=ggplot(dplyr::filter(expression_data_df,patient=='patient4'),
        aes_string(x=feature1, y=feature2))+
        geom_point()+
        geom_smooth(method='lm')+
        theme_bw()+ggtitle('Pt.4')
    p3=ggplot(dplyr::filter(expression_data_df,patient=='patient5'),
        aes_string(x=feature1, y=feature2))+
        geom_point()+
        geom_smooth(method='lm')+
        theme_bw()+ggtitle('Pt.5')
    p1+p2+p3
}
scatter_pt4pt5(expression_data_df,feature1='FSCA', feature2='ACTA1__chr1')

FSCA_selection=
    expression_data_df$FSCA>boundaries_FSC_A[1]&expression_data_df$FSCA<boundaries_FSC_A[2]
scatter_pt4pt5(expression_data_df[FSCA_selection,],feature1='FSCA', feature2='ACTA1__chr1')

ggplot(expression_data_df)+
    geom_point(aes(x=FSCA, y=ACTA1__chr1, color=patient))+theme_bw()

##########

"scatter_pt4pt5 = function(expression_data_df, percentile=0.05, feature1, feature2) {
    
    boundaries_FSC_A = calc_limits(expression_data_df$FSCA, percentile = percentile)
    selection = expression_data_df$FSCA>boundaries_FSC_A[1]&expression_data_df$FSCA<boundaries_FSC_A[2]
    
    expression_data_df_sel = expression_data_df[selection,]
    
    totalCor = cor.test(expression_data_df_sel[,feature1],expression_data_df_sel[,feature2])
    p1=ggplot(expression_data_df_sel, aes_string(x=feature1, y=feature2))+
        geom_point(data=expression_data_df)+
        geom_point()+
        geom_point(aes_string(color='patient'))+
        geom_smooth(method='lm',aes_string(color='patient'))+
        geom_smooth(method='lm',color='black')+
        ggtitle(paste0('Correlation: ',round(totalCor$estimate,2),'(non-adjusted p-value: ',totalCor$p.value,')'))+
        theme_bw()
    p1
    p2=ggplot(dplyr::filter(expression_data_df,patient=='patient4'),
        aes_string(x=feature1, y=feature2))+
        geom_point()+
        geom_smooth(method='lm')+
        theme_bw()
    p3=ggplot(dplyr::filter(expression_data_df,patient=='patient5'),
        aes_string(x=feature1, y=feature2))+
        geom_point()+
        geom_smooth(method='lm')+
        theme_bw()
    p1+p2+p3
}"


################################################################################
################################################################################
# Yet another method, combining the two count tables by
# beforehand normalizing per patient

# Let's create convenient parameters to handle the indexdata

# Create z-score params for patient 4
sc.gene.in.cells.pt4 = apply(sc.countdata.pt4>0.1,1,sum)
sc.gene.sel.pt4 = sc.gene.in.cells.pt4>round(dim(sc.countdata.pt4)[2]*.1) &
                    !is.na(sc.gene.in.cells.pt4)
sc.countdata.pt4 = t(scale(t(groupedSCS$patient4Mod@ndata[sc.gene.sel.pt4,])))
sc.indexData.pt4 = scale(indexData_MW[colnames(groupedSCS$patient4Mod@ndata),]$FSC_A)

# Create z-score params for patient 5
sc.gene.in.cells.pt5 = apply(sc.countdata.pt5>0.1,1,sum)
sc.gene.sel.pt5 = sc.gene.in.cells.pt5>round(dim(sc.countdata.pt5)[2]*.1) &
                    !is.na(sc.gene.in.cells.pt5)
sc.countdata.pt5 = t(scale(t(groupedSCS$patient5Mod@ndata[sc.gene.sel.pt5,])))
sc.indexData.pt5 = scale(indexData_MW[colnames(groupedSCS$patient5Mod@ndata),]$FSC_A)

# Combine the patients
shared_genes = rownames(sc.countdata.pt4)[rownames(sc.countdata.pt4) %in% rownames(sc.countdata.pt5)]
sc.countdata.p4p5 = cbind(sc.countdata.pt4[shared_genes,],sc.countdata.pt5[shared_genes,])
sc.indexData.p4p5 = c(as.vector(sc.indexData.pt4), as.vector(sc.indexData.pt5))


#sc.fsca.cormat = cor.test(sc.countdata.p4p5[2,], sc.indexData.p4p5)
#sc.gene.in.cells = apply(sc.countdata.p4p5>0,1,sum)
#sc.gene.sel = sc.gene.in.cells>100
sc.cor.test.p4.p5.out = as.data.frame(t(apply(as.matrix(sc.countdata.p4p5),1,function(X) {thecor=cor.test(X, sc.indexData.p4p5); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
sc.cor.test.p4.p5.out$pval.adj = p.adjust(sc.cor.test.p4.p5.out$p.value, method="BH")
sc.cor.test.p4.p5.out$corr = sc.cor.test.p4.p5.out$cor.cor
sc.cor.test.p4.p5.out$gene_name = row.names(sc.cor.test.p4.p5.out)
sc.cor.test.p4.p5.out$gene_name_short = strip__chrXX(row.names(sc.cor.test.p4.p5.out))

plot_volcano(sc.cor.test.p4.p5.out, manual_gene_name = 'FSC-A (pt.4+5) [n]', NRLABELED = 10, mypvaltreshold = 0.05)

sc.cor.test.p4.p5.out[sc.cor.test.p4.p5.out$gene_name_short=='ACTA1',]

shorthand_scatter_fsca_gene_ptcol( gene_expression = get_expression_gene2(sc.countdata.p4p5, 'REV1'),
    FSC_A = sc.indexData.p4p5, mytitle = 'Pt4+Pt5 [n]',gene_name = 'REV1',patient_annotation = patient_source_pt4.pt5[colnames(sc.countdata.p4p5)])
        
################################################################################
# Look at a randomized matrix as reference
# This shows that you can get pretty high correlations in random data

# Create shuffled matrix
sc.countdata.p4p5.shuffled = t(apply(sc.countdata.p4p5, 1, sample))
sh.sc.cor.test.p4.p5.out = as.data.frame(t(apply(as.matrix(sc.countdata.p4p5.shuffled),1,function(X) {thecor=cor.test(X, sc.indexData.p4p5); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
sh.sc.cor.test.p4.p5.out$pval.adj = p.adjust(sh.sc.cor.test.p4.p5.out$p.value, method="BH")
sh.sc.cor.test.p4.p5.out$corr = sh.sc.cor.test.p4.p5.out$cor.cor
sh.sc.cor.test.p4.p5.out$gene_name = row.names(sh.sc.cor.test.p4.p5.out)
sh.sc.cor.test.p4.p5.out$gene_name_short = strip__chrXX(row.names(sh.sc.cor.test.p4.p5.out))

ggplot(sh.sc.cor.test.p4.p5.out)+
    give_better_textsize_plot(12)+
    geom_histogram(aes(x=corr),stat='density')+theme_bw()+
    ggtitle('Randomized expression')

max.cor=c()
for (trial in 1:10) {
    print(paste('trial = ', trial))
    sc.countdata.p4p5.shuffled = t(apply(sc.countdata.p4p5, 1, sample))
    sh.sc.cor.test.p4.p5.out = as.data.frame(t(apply(as.matrix(sc.countdata.p4p5.shuffled),1,function(X) {thecor=cor.test(X, sc.indexData.p4p5); return(c(cor=thecor$estimate, p.value=thecor$p.value))})))
    max.cor[[trial]]=max(sh.sc.cor.test.p4.p5.out$cor.cor)
}

"
expressionData_FSCA = rbind(expressionData_FSCA, as.data.frame(FSCA_matrix))
    # sanity check (outcome should be false)
    any(colnames(groupedSCS_pt4.pt5$pt4.pt5Mod@ndata) != indexData_pt4pt5$Cell)
# also add patient annotation
indexData_pt4pt5$Patient = patient_source_pt4.pt5
"
################################################################################
