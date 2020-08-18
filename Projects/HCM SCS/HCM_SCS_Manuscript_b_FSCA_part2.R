

######################################################################
# Please see notes in "HCM SCS Manuscript_AE.R" first
#
# This script performs some additional analysis:
# a) total transcript count vs. FSCA analysis
# b) project FSC-A on tSNE map
# c) check how "red" genes are correlated to each other

################################################################################
# Perform total transcript count vs. FSCA analysis

# Load the applicable data for patient 4 and 5
config2 = config
config2$custom_gene_blacklist = NULL
config2$excludeMitochondrialGenes = F
groupedData_raw_p4p5 <- loadData(config2, c('patient4', 'patient5'))
groupedData_raw_p4p5 <- applyExclusions(config2, groupedData_raw_p4p5) # only removes spike-ins with settings above

# First calculate cell totals
raw_cell_count_totals_p4 = apply(groupedData_raw_p4p5$patient4,2,sum)
raw_cell_count_totals_p5 = apply(groupedData_raw_p4p5$patient5,2,sum)
# Then get FSC-A for those cells
#FSCA_p4 = indexData_pt4pt5[names(raw_cell_count_totals_p4),]$FSC_A # Ooops, this only contains filtered cells
#FSCA_p5 = indexData_pt4pt5[names(raw_cell_count_totals_p5),]$FSC_A
FSCA_p4 = indexData_MW[names(raw_cell_count_totals_p4),]$FSC_A
FSCA_p5 = indexData_MW[names(raw_cell_count_totals_p5),]$FSC_A

# Remove the raw data
rm('groupedData_raw_p4p5')

# A few outliers are messing up the corrs, remove those
cor.test.out4=cor.test(raw_cell_count_totals_p4[raw_cell_count_totals_p4<30000&raw_cell_count_totals_p4>500],FSCA_p4[raw_cell_count_totals_p4<30000&raw_cell_count_totals_p4>500])
cor.test.out5=cor.test(raw_cell_count_totals_p5[raw_cell_count_totals_p5<30000&raw_cell_count_totals_p5>500],FSCA_p5[raw_cell_count_totals_p5<30000&raw_cell_count_totals_p5>500])
cor.test.out4
cor.test.out5

cor.test(raw_cell_count_totals_p4,FSCA_p4)
cor.test(raw_cell_count_totals_p5,FSCA_p5)

# Create plots that scatter the two
# first create dataframes
df_toplot_pt4=data.frame(total_count=raw_cell_count_totals_p4, FSC_A=FSCA_p4, patient=4)
df_toplot_pt5=data.frame(total_count=raw_cell_count_totals_p5, FSC_A=FSCA_p5, patient=5)

#ggplot(data = data.frame(x = c(min(df_toplot_pt4$total_count),max(df_toplot_pt4$total_count))), aes(x)) +
#  stat_function(fun = dnorm, n = 101, args = list(mean = mean(df_toplot_pt4$total_count), sd = sd(df_toplot_pt4$total_count))) + ylab("") +
#  scale_y_continuous(breaks = NULL)

# look at some stats
mean_pt4=mean(df_toplot_pt4$total_count)
sd_pt4=sd(df_toplot_pt4$total_count)
mean_pt5=mean(df_toplot_pt5$total_count)
sd_pt5=sd(df_toplot_pt5$total_count)
mean_pt4+5*sd_pt4 # could be used as cutoffs
mean_pt5+5*sd_pt5

OUTLIER_CUTOFF = 30000 # this seams a realistic cutoff
MIN_READCOUNT = 500

# Plot the histogram
df_histplot = rbind(df_toplot_pt4,df_toplot_pt5)
df_histplot$patient=as.factor(df_histplot$patient)
pH=ggplot(df_histplot) +
    geom_jitter(data=df_histplot[is.na(df_histplot$FSC_A),],aes(x=patient, y=log10(1+total_count)),size=.1,color='red')+
    geom_jitter(data=df_histplot[!is.na(df_histplot$FSC_A),],aes(x=patient, y=log10(1+total_count)),size=.1,color='black')+
    #geom_jitter(aes(x=patient, y=log10(1+total_count)),size=.1)+
    geom_violin(data=df_histplot[!is.na(df_histplot$FSC_A),],aes(x=patient, y=log10(1+total_count)),alpha=.8,fill='gray',color='black')+
    #geom_boxplot(aes(x=patient, y=log10(1+total_count),fill=as.factor(patient)), alpha=.5)+
    geom_hline(yintercept = log10(1+OUTLIER_CUTOFF))+theme_bw()+shorthand_tsne_joeptheme()+
    geom_hline(yintercept = log10(1+MIN_READCOUNT))+theme_bw()+shorthand_tsne_joeptheme()+
    ylab('log10(1+UMI count)')+xlab('Patient')+give_better_textsize_plot(8)
pH
shorthand_save(pH, 
    thedir = 'MW_custom_plots/', filename = 'fsca_totalreadcount_histogramForLims.pdf', mywidth = 7.5, myheight = 7.5)

sel1=df_toplot_pt4$total_count<=OUTLIER_CUTOFF&df_toplot_pt4$total_count>MIN_READCOUNT
p1=ggplot(df_toplot_pt4[sel1,],aes(x=FSC_A, y=total_count))+
    geom_point(data=df_toplot_pt4[!sel1,],color='red',size=.25)+
    geom_point(color='black',size=.25)+
    geom_smooth(method='lm')+
    theme_bw()+give_better_textsize_plot(14)+ggtitle('Patient 4')+xlab('Cell size (FSCA)')+ylab('mRNA count \n(UMI count)')
sel2=df_toplot_pt5$total_count<=OUTLIER_CUTOFF&df_toplot_pt5$total_count>MIN_READCOUNT
p2=ggplot(df_toplot_pt5[sel2,],aes(x=FSC_A, y=total_count))+
    geom_point(data=df_toplot_pt5[!sel2,],color='red',size=.25)+
    geom_point(color='black',size=.25)+
    geom_smooth(method='lm')+
    theme_bw()+give_better_textsize_plot(14)+ggtitle('Patient 5')+xlab('Cell size (FSCA)')+ylab('mRNA count \n(UMI count)')
p1+p2+p1+ylim(c(0,30000))+p2+ylim(c(0,30000))

# Make plot with correlations also
format(round(cor.test.out4$estimate, 2), nsmall = 2)
formatC(cor.test.out4$p.value, format = "e", digits = 2)
format(round(cor.test.out5$estimate, 2), nsmall = 2)
formatC(cor.test.out5$p.value, format = "e", digits = 2)
pO = p1+ylim(c(0,30000))+give_better_textsize_plot(8)+
        ggtitle(paste0('Patient 4 \nCorr = ',round(cor.test.out4$estimate, 2) ,', p=',formatC(cor.test.out4$p.value, format = "e", digits = 2)))+
     p2+ylim(c(0,30000))+give_better_textsize_plot(8)+
        ggtitle(paste0('Patient 5 \nCorr = ',round(cor.test.out5$estimate, 2) ,', p=',formatC(cor.test.out5$p.value, format = "e", digits = 2)))
pO
# export scatter plots
shorthand_save(pO, 
    thedir = 'MW_custom_plots/', filename = 'fsca_totalreadcount.pdf', mywidth = 15, myheight = 7.5)

# export more formats to file
shorthand_save(p1+p2, 
    thedir = 'MW_custom_plots/', filename = 'fsca_totalreadcount_x.pdf', mywidth = 15, myheight = 7.5)
shorthand_save(p1+ylim(c(0,30000))+p2+ylim(c(0,30000)), 
    thedir = 'MW_custom_plots/', filename = 'fsca_totalreadcount_lims.pdf', mywidth = 15, myheight = 7.5)



################################################################################
# Now project FSC-A on tSNE map

# Patient 4
FSCA_p4=indexData_pt4pt5[rownames(groupedSCS$patient4Mod@tsne),]$FSC_A
p1=plot_sth_in_XY(comp_1 = groupedSCS$patient4Mod@tsne[,1],
         comp_2 = groupedSCS$patient4Mod@tsne[,2],
         color_by = FSCA_p4,
         name_for_color = 'FSCA', print_yes = F)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne patient 4'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p4),oob=squish)#+
         #theme(legend.position = 'none')
print(p1)
# Patient 5
FSCA_p5=indexData_pt4pt5[rownames(groupedSCS$patient5Mod@tsne),]$FSC_A
p2=plot_sth_in_XY(comp_1 = groupedSCS$patient5Mod@tsne[,1],
         comp_2 = groupedSCS$patient5Mod@tsne[,2],
         color_by = FSCA_p5,
         name_for_color = 'FSCA', print_yes = F)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne patient 5'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p5),oob=squish)#+
         #theme(legend.position = 'none')
print(p2)

library(patchwork)
p1+p2

# Both
FSCA_p4p5=indexData_pt4pt5[rownames(groupedSCS$patientAllWithIndexMod@tsne),]$FSC_A
p3=plot_sth_in_XY(comp_1 = groupedSCS$patientAllWithIndexMod@tsne[,1],
         comp_2 = groupedSCS$patientAllWithIndexMod@tsne[,2],
         color_by = FSCA_p4p5,
         name_for_color = 'FSCA', print_yes = F)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne patient 4+5'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p4p5),oob=squish)#+
         #theme(legend.position = 'none')
print(p3)

################################################################################
# Just project on big thing, both before and after filtering

# First without filtering
#mycells = rownames(groupedSCS$patientAll@tsne)[rownames(groupedSCS$patientAll@tsne) %in% indexData_pt4pt5$Cell]
FSCA_all=indexData_pt4pt5[rownames(groupedSCS$patientAll@tsne),]$FSC_A
p_all=plot_sth_in_XY(comp_1 = groupedSCS$patientAll@tsne[,1],
         comp_2 = groupedSCS$patientAll@tsne[,2],
         color_by = FSCA_all,
         name_for_color = 'FSCA', print_yes = F, mypointsize = .5)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p5),oob=squish)#+
         #theme(legend.position = 'none')
print(p_all)
# export
p_all = p_all+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(p_all, thedir = 'MW_custom_plots/', filename = 'tSNE_unf_cellsize_all.png', mywidth = 7.5, myheight = 7.5)

# then including filtering
#mycells = rownames(groupedSCS$patientAll@tsne)[rownames(groupedSCS$patientAll@tsne) %in% indexData_pt4pt5$Cell]
FSCA_all=indexData_pt4pt5[rownames(groupedSCS$patientAllMod@tsne),]$FSC_A
p_all=plot_sth_in_XY(comp_1 = groupedSCS$patientAllMod@tsne[,1],
         comp_2 = groupedSCS$patientAllMod@tsne[,2],
         color_by = FSCA_all,
         name_for_color = 'FSCA', print_yes = F, mypointsize = .5)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p5),oob=squish)#+
         #theme(legend.position = 'none')
print(p_all)
# export
p_all = p_all+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(p_all, thedir = 'MW_custom_plots/', filename = 'tSNE_f_cellsize_all.png', mywidth = 7.5, myheight = 7.5)


# Now the same but don't show the grey cells
cells_to_get=rownames(groupedSCS$patientAllWithIndexMod@tsne)
FSCA_p4p5=indexData_pt4pt5[cells_to_get,]$FSC_A
p_all2=plot_sth_in_XY(comp_1 = groupedSCS$patientAllMod@tsne[cells_to_get,1],
         comp_2 = groupedSCS$patientAllMod@tsne[cells_to_get,2],
         color_by = FSCA_p4p5,
         name_for_color = 'FSCA', print_yes = F, mypointsize = .5)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(FSCA_p5),oob=squish)+
    xlim(c(min(groupedSCS$patientAllMod@tsne[,1]),max(groupedSCS$patientAllMod@tsne[,1])))+
    ylim(c(min(groupedSCS$patientAllMod@tsne[,2]),max(groupedSCS$patientAllMod@tsne[,2])))#+
         #theme(legend.position = 'none')
print(p_all2)
# export
p_all2 = p_all2+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(p_all2, thedir = 'MW_custom_plots/', filename = 'tSNE_f_cellsize_all_hideNA.png', mywidth = 7.5, myheight = 7.5)

################################################################################
# Now check how "red" genes are correlated to each other

correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1 & correlations_pt4.pt5$distance<.05,]$gene

red_genes_names_short = correlations_pt4.pt5[correlations_pt4.pt5$cor.pt4>.1 & correlations_pt4.pt5$cor.pt5>.1,]$gene

all_gene_names = rownames(groupedSCS$patientAllWithIndexMod@ndata)

red_genes_names_long = find__chrXX(gene_query = red_genes_names_short,  gene_names = all_gene_names)

red_genes_expression = groupedSCS$patientAllWithIndexMod@ndata[red_genes_names_long,]

red_cor_matrix = t(cor(t(red_genes_expression)))

library(pheatmap)
max_corr = max(red_cor_matrix[red_cor_matrix!=1])
pheatmap(red_cor_matrix, breaks = seq(0, max_corr, .01))

# in all data
red_genes_expression = groupedSCS$patientAllMod@ndata[red_genes_names_long,]
red_cor_matrix = t(cor(t(red_genes_expression)))
library(pheatmap)
max_corr = max(red_cor_matrix[red_cor_matrix!=1])
pheatmap(red_cor_matrix, breaks = seq(0, max_corr, .01))

# in patient 5
red_genes_expression = groupedSCS$patient5Mod@ndata[red_genes_names_long,]
red_cor_matrix = t(cor(t(red_genes_expression)))
library(pheatmap)
max_corr = max(red_cor_matrix[red_cor_matrix!=1])
pheatmap(red_cor_matrix, breaks = seq(0, max_corr, .01))


