

# let's look some more in mitochondrial stuff
config2=config
config2$excludeMitochondrialGenes=F
groupedData_InclMito <- loadData(config2, c('patientAll'))
filtered_cells = colnames(groupedSCS$patientAllMod@ndata)

mitoGenes = grepl('^MT-',rownames(groupedData_InclMito$patientAll))
mitoCountTable = groupedData_InclMito$patientAll[mitoGenes,]
mitoCountTable_filt = groupedData_InclMito$patientAll[mitoGenes,filtered_cells]

totalReads = apply(groupedData_InclMito$patientAll,2,sum)
totalMitoReads = apply(mitoCountTable,2,sum)

totalReads_filt = apply(groupedData_InclMito$patientAll[,filtered_cells],2,sum)
totalMitoReads_filt = apply(mitoCountTable_filt,2,sum)

fractionMitoReads = totalMitoReads/totalReads
fractionMitoReads_filt = totalMitoReads_filt/totalReads_filt

rm('groupedData_InclMito')

# Now create violin plot per patient with mito fractions / cell
patient_annotation_mito_fracs = shorthand_give_p_source_forcellnames(config, names(fractionMitoReads_filt))
patient_annotation_mito_fracs = gsub('patient','',patient_annotation_mito_fracs)
mito_frac_df = data.frame(fractionMitoReads_filt = fractionMitoReads_filt, patient = patient_annotation_mito_fracs)
p=ggplot(mito_frac_df)+
    geom_violin(aes(x=patient, y=fractionMitoReads_filt))+
    ylim(c(0,1))+
    theme_bw()+
    give_better_textsize_plot(10)+
    ylab('Fraction mitochondrial reads')+xlab('Patient')+
    shorthand_tsne_joeptheme()
p
shorthand_save(p, thedir = 'analysis_2020-05-04-18h_17m/FigX/', filename = 'MitochondrialFractReads_Violins.pdf', mywidth = 7.5, myheight = 7.5)

# Also show stat regarding this data
#mito_frac_df$patient_str=paste0('p', mito_frac_df$patient)
#aggregate(mito_frac_df, by=list(mito_frac_df$patient_str), mean)$fractionMitoReads_filt
mean(mito_frac_df$fractionMitoReads_filt)
sd(mito_frac_df$fractionMitoReads_filt)

# =========
# Now print this on the tSNE map
MitoFract_forTsne=fractionMitoReads_filt[rownames(groupedSCS$patientAllMod@tsne)]
p_MtFract=plot_sth_in_XY(comp_1 = groupedSCS$patientAllMod@tsne[,1],
         comp_2 = groupedSCS$patientAllMod@tsne[,2],
         color_by = MitoFract_forTsne,
         name_for_color = 'Mito_Fract', print_yes = F, mypointsize = .5)+
         theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne'))+
         scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(MitoFract_forTsne),oob=squish)#+
         #theme(legend.position = 'none')
print(p_MtFract)
# export
p_MtFract = p_MtFract+shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())

p_MtFract = p_MtFract+shorthand_tsne_joeptheme()+ xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
p_MtFract
shorthand_save(p_MtFract, thedir = 'MW_custom_plots/', filename = 'tSNE_MitoFract.png', mywidth = 7.5, myheight = 7.5)

# now obtain regulon composite expression again
regulon_composite_expression=list()
for (reg_idx in 1:4){
    regulon_composite_expression[[reg_idx]] = 
        shorthand_composite_expression_values(current_SCS = groupedSCS$patientAllMod, current_gene_set = shared_regulon_core_list[[reg_idx]])
}

# also retrieve FSCA per cells accordingly
FSCA_forCorr = indexData_pt4pt5[names(regulon_composite_expression[[1]]),]$FSC_A
names(FSCA_forCorr) = names(regulon_composite_expression[[1]])
# mito fraction accordingly 
MitoFract_forCorr = fractionMitoReads[names(regulon_composite_expression[[1]])]

# check some correlations
cor.test(MitoFract_forCorr,regulon_composite_expression[[1]])
cor.test(MitoFract_forCorr,regulon_composite_expression[[2]])
cor.test(MitoFract_forCorr,regulon_composite_expression[[3]])
cor.test(MitoFract_forCorr,regulon_composite_expression[[4]])

cor.test(regulon_composite_expression[[1]],regulon_composite_expression[[2]])
cor.test(regulon_composite_expression[[1]],regulon_composite_expression[[3]])
cor.test(regulon_composite_expression[[1]],regulon_composite_expression[[4]])

FSCA_forCorr_noNA = FSCA_forCorr[!is.na(FSCA_forCorr)]
cor(FSCA_forCorr_noNA, multiple_Params_Comparison[names(FSCA_forCorr_noNA),])

# systematically compare correlations
multiple_Params_Comparison = data.frame( FSCA_forCorr,
            MitoFract_forCorr, 
            regulon_composite_expression[[1]],
            regulon_composite_expression[[2]], 
            regulon_composite_expression[[3]], 
            regulon_composite_expression[[4]] )
    
corParams = cor(multiple_Params_Comparison)

pheatmap(corParams, cluster_rows = F, cluster_cols = F, breaks=seq(-.4,.4,.01))

# let's put in FSCA-correlations, based only on cells where we have data
# Note: this perhaps combines things that shouldn't be combined, because
# the correlations in this matrix will now be calculated based of 
# different cellular subsets
corParams2 = corParams
corParams2[1,] = cor(FSCA_forCorr_noNA, multiple_Params_Comparison[names(FSCA_forCorr_noNA),])
corParams2[,1] = cor(FSCA_forCorr_noNA, multiple_Params_Comparison[names(FSCA_forCorr_noNA),])

pheatmap(corParams2, cluster_rows = F, cluster_cols = F, breaks=seq(-.4,.4,.01))
pheatmap(corParams2, breaks=seq(-.4,.4,.01))


