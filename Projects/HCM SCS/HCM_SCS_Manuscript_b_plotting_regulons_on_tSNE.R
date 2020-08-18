

library(patchwork)


# make list of regulon genes
shared_regulon_genes_list = list(genes_set_of_interest_sel_1, genes_set_of_interest_sel_2, genes_set_of_interest_sel_3, genes_set_of_interest_sel_4)

# plot 'm & save 'm
listP1=list()
for (idx in 1:length(shared_regulon_genes_list)) {
    
    listP1[[idx]]=shorthand_composite_expression(current_SCS = groupedSCS$patientAllMod, current_gene_set = shared_regulon_genes_list[[idx]],mypointsize = 1.5)+
        ggtitle(paste0('Shared regulon ', idx))+
            labs(x=element_blank(),y=element_blank())+shorthand_joeptheme()
    #if (idx<length(shared_regulon_genes_list)) {listP1[[idx]]=listP1[[idx]]+theme(legend.position = 'none')}
    
    shorthand_save(listP1[[idx]], 'analysis_2020-05-04-18h_17m/FigX/', paste0('tSNE_allPts_sharedRegulon',idx,'_expression.pdf'), mywidth = 7*2.54, myheight = 5.3*2.54)

}
wrap_plots(listP1, nrow = 2)


# Now for separate patients
uberListP = list()
for (reg_idx in 1:4) { # loop over shared regulons
    listP = list()
    for (patientNr in 1:5) { # loop over patients
        listP[[patientNr]]=shorthand_composite_expression(current_SCS = groupedSCS[[paste0('patient',patientNr,'Mod')]], 
            current_gene_set = shared_regulon_genes_list[[reg_idx]],
            mypointsize = .5)+
            ggtitle(paste0('Reg-',reg_idx, ', pt-',patientNr))+
            labs(x=element_blank(),y=element_blank())+shorthand_joeptheme()
        # remove title for all but last
        if (patientNr<5) {listP[[patientNr]]=listP[[patientNr]]+theme(legend.position = 'none')}
    }
    
    print(wrap_plots(listP, nrow = 1))
    uberListP[[reg_idx]]=listP
}

wrap_plots(listP, nrow = 1)

# save all
for (reg_idx in 1:4) {
    p=wrap_plots(uberListP[[reg_idx]], nrow = 1)
    shorthand_save(p, 'analysis_2020-05-04-18h_17m/FigX/', paste0('tSNE_regulon',reg_idx,'_expression.pdf'), mywidth = 20, myheight = 4)
}



 


