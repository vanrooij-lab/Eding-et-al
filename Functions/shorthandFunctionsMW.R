shorthand_XY_plot_clusters = function(current_SCS,mytitle,cl2=F) {
    clusters = as.factor(current_SCS@cluster$kpart)
    if (cl2) { clusters = as.factor(current_SCS@cpart) }
    cl_centers = get_cluster_centers(df_tsne = data.frame(
        tsne1=current_SCS@tsne[,1],tsne2=current_SCS@tsne[,2],cluster=clusters
        ),varname1 = 'tsne1',varname2 = 'tsne2',clustervarname = 'cluster')
    p1=plot_sth_in_XY(comp_1 = current_SCS@tsne[,1],
             comp_2 = current_SCS@tsne[,2],
             color_by = clusters,
             name_for_color = 'cluster', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('r-tsne: ',mytitle))+
             geom_text(data=cl_centers,aes(x=coord1,y=coord2,label=clName),fontface = 'bold')
        #+
        #theme(legend.position = 'none')
    return(p1)
}
shorthand_expression = function(some_gene,current_SCS,mytitle,textsize=10) {
    Expr=get_expression_gene2(
            expression_matrix = as.matrix(current_SCS@ndata),
            gene_query = some_gene)
    p1=plot_sth_in_XY(comp_1 = current_SCS@tsne[,1],
             comp_2 = current_SCS@tsne[,2],
             color_by = Expr,
             name_for_color = some_gene, print_yes = F)+
             theme_bw()+give_better_textsize_plot(textsize)+ggtitle(paste0('tsne: ',mytitle))+
             scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(Expr),oob=squish)+
             xlab('tSNE-1')+ylab('tSNE-2')#+
             #theme(legend.position = 'none')
    print(p1)
    
    return(p1)
}
shorthand_selection = function(the_selection,current_SCS,mytitle,sel_title,mypointsize=1) {
    plot_sth_in_XY(comp_1 = current_SCS@tsne[,1],
             comp_2 = current_SCS@tsne[,2],
             color_by = the_selection,
             name_for_color = sel_title, print_yes = F, mypointsize=mypointsize)+
             theme_bw()+give_better_textsize_plot(15)+ggtitle(paste0('tsne: ',mytitle))   
}
shorthand_histogram = function(some_gene,current_SCS,mytitle) {
    Expr=get_expression_gene2(
            expression_matrix = as.matrix(current_SCS@ndata),
            gene_query = some_gene)
    p1=ggplot(data.frame(expression=Expr))+
        geom_histogram(aes(x=expression), bins=50)+theme_bw()+ggtitle(paste0('Expression of ',some_gene))+
        give_better_textsize_plot(14)
    p2=ggplot(data.frame(expression=Expr))+
        geom_histogram(aes(x=expression), bins=50)+theme_bw()+
        give_better_textsize_plot(14)+xlim(c(0,calc_limits(Expr,0.05)[2]))#+ggtitle(paste0('Expression of ',some_gene))
    p1+p2+plot_layout(ncol = 1)
    
}

# Violin expr / cluster
shorthand_expr_per_cluster = function(some_gene, current_SCS, cl2=F) {
    
    clusters = as.factor(current_SCS@cluster$kpart)
    if (cl2) { clusters = as.factor(current_SCS@cpart) }

    Expr=get_expression_gene2(
            expression_matrix = as.matrix(current_SCS@ndata),
            gene_query = some_gene)

    df_toplot = data.frame(expr=Expr,cluster=clusters)

    ggplot(df_toplot)+
        geom_jitter(aes(x=cluster, y=log10(.1+expr)))+
        #geom_violin(aes(x=cluster, y=log10(.1+expr)), alpha=.5)+
        geom_boxplot(aes(x=cluster, y=log10(.1+expr)), alpha=.5)+
        theme_bw()+ggtitle(paste0(some_gene, ' expression per cluster'))
}

shorthand_give_p_source = function(config, current_SCS, additional_suffix = '_') {

    print('Set additional_suffix correctly!')
  
    cellnames = colnames(current_SCS@ndata)
    patient_source = rep(NA, length(cellnames))
    
    
    for (current_patient in names(config$samples)) {
        current_patient_ids=paste0(config$samples[[current_patient]], additional_suffix)
        patient_source[grepl(paste(current_patient_ids, collapse = '|'), cellnames)] = current_patient
    }
    
    return(patient_source)
    
}
shorthand_give_p_source_forcellnames = function(config, cellnames, additional_suffix = '_') {

    patient_source = rep(NA, length(cellnames))
    
    for (current_patient in names(config$samples)) {
        current_patient_ids=paste0(config$samples[[current_patient]], additional_suffix)
        patient_source[grepl(paste(current_patient_ids, collapse = '|'), cellnames)] = current_patient
    }
    
    return(patient_source)
    
}

shorthand_save = function(p,thedir,filename, myunits='cm', mywidth=15, myheight=10){
    
    print("Saving to workdir..")
    if (!dir.exists(paste0(getwd(),'/',thedir))) {dir.create(paste0(getwd(),'/',thedir),recursive = TRUE)}
    ggsave(filename = paste0(getwd(),'/',thedir,'/',filename), plot=p, width = mywidth, height = myheight, units = 'cm')
    
}

shorthand_pretty_tsne = function(mylegendsize=NULL){
    
        returnvalue = theme_minimal()+theme(
            axis.text.x = element_blank(),axis.text.y = element_blank(),
            axis.ticks = element_blank(),panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
        if (!is.null(mylegendsize)) {
            returnvalue = returnvalue+theme(legend.key.size = unit(mylegendsize,'cm')) }
        return(returnvalue)
        
        #+
        #xlab('tSNE1')+ylab('tSNE2')+
}

# used in "shorthand_composite_expression()"
scale_non_zero = function(X, threshold=0) {
    X[X>threshold] = scale(X[X>threshold], center = F, scale = T)
    return(as.vector(X))
}

shorthand_composite_expression_values = function(current_SCS, current_gene_set, show10violins=F) {

    # get full gene names
    current_gene_set_chrXX = find__chrXX(gene_query = current_gene_set, gene_names = rownames(current_SCS@ndata))
    # Get desired subset of expression matrix
    current_composite_matrix = current_SCS@ndata[current_gene_set_chrXX, ]
    
    # Now scale the matrix to make genes "comparable" in expression
    current_composite_matrix_scaled = t(apply(current_composite_matrix-.1,1,scale_non_zero))
    colnames(current_composite_matrix_scaled) = colnames(current_composite_matrix)
    
    # show effect of normalization of genes visually
    if (show10violins) {
        p=ggplot(melt(t(current_composite_matrix_scaled[1:10,])))+
            geom_violin(aes(x=Var2,y=value))
        print(p)
    }
    
    # Alternative scaling
    #current_composite_matrix_scaled = t(scale(t(current_composite_matrix), center = T, scale = T))
    #ggplot(melt(t(current_composite_matrix_scaled[1:10,])))+
    #    geom_violin(aes(x=Var2,y=value))
    
    # Calculate composite expression (mean of normalized expression)
    composite_expression = apply(current_composite_matrix_scaled, 2, mean)
    
    return(composite_expression)
}

shorthand_composite_expression = function(current_SCS, current_gene_set, show10violins=F, mypointsize=1) {

    # get full gene names
    current_gene_set_chrXX = find__chrXX(gene_query = current_gene_set, gene_names = rownames(current_SCS@ndata))
    # Get desired subset of expression matrix
    current_composite_matrix = current_SCS@ndata[current_gene_set_chrXX, ]
    
    # Now scale the matrix to make genes "comparable" in expression
    current_composite_matrix_scaled = t(apply(current_composite_matrix-.1,1,scale_non_zero))
    
    # show effect of normalization of genes visually
    if (show10violins) {
        p=ggplot(melt(t(current_composite_matrix_scaled[1:10,])))+
            geom_violin(aes(x=Var2,y=value))
        print(p)
    }
    
    # Alternative scaling
    #current_composite_matrix_scaled = t(scale(t(current_composite_matrix), center = T, scale = T))
    #ggplot(melt(t(current_composite_matrix_scaled[1:10,])))+
    #    geom_violin(aes(x=Var2,y=value))
    
    # Calculate composite expression (mean of normalized expression)
    composite_expression = apply(current_composite_matrix_scaled, 2, mean)
    
    # Plot composite expression
    plot_sth_in_XY(current_SCS@tsne[,1],current_SCS@tsne[,2],
        color_by = composite_expression, name_for_color = 'comp_expr', print_yes = F, mypointsize = mypointsize)+
        theme_bw()+give_better_textsize_plot(10)+
        scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(composite_expression),oob=squish)

}
# alias
shorthand_plot_composite_expression = shorthand_composite_expression


# "Joep's theme" for tSNE maps
# Note: point size should be 1.5
# figure is typically saved to sizes mywidth = 7*2.54, myheight = 5.3*2.54 (cms)
# 
shorthand_joeptheme = function() {
    theme(
        aspect.ratio = config$plotOptions$aspect.ratio,
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
      ) 
}


shorthand_GOplot = function(GO_df, topX=10, minGeneCount=0, mycolor='black') {
    
    topX_ = min(topX, dim(GO_df)[1])
    selGC = GO_df$Count>=minGeneCount
    
    p <- ggplot(data=GO_df[selGC,][order(GO_df[selGC,]$Pvalue),][1:topX_,],
                aes(x=reorder(Term,-Pvalue),y=-log10(Pvalue))) + 
            geom_bar(stat='identity', fill=mycolor) + coord_flip() + 
            labs( x="Term", y="-log10(p-value)")+theme_bw()

    return(p)
    
}


# patient contributions to clusters
# current_SCS = groupedSCS$patientAllMod
shorthand_plot_patient_contr_clusters = function(config, current_SCS) {
  misc_cluster_assignments = current_SCS@cluster$kpart
  misc_patient_assignments = shorthand_give_p_source(config, current_SCS)
  
  p=ggplot(data.frame(  cluster_assignments=as.factor(misc_cluster_assignments), 
                      patient_source=as.factor(misc_patient_assignments)))+
      geom_histogram(aes(x=cluster_assignments, fill=patient_source), stat='count')+theme_bw()
  
  return(p)
}




