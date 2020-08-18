
library(ggalt)
library(scales) # Needed for squish in oob setting

# Input is freq_frame, which has the structure of a melt of multiple histogram lines
# Should have the fields centers, counts, condition
plot_multiple_hists<-function(freq_frame) {
  
  # Plot the last one
  TEXTSIZE=15
  p<-ggplot(data=freq_frame, mapping=aes(x=centers, y=counts,fill=condition)) +
    #geom_bar(stat="identity")+
    geom_line(aes(color=condition))+
    #geom_point()+
    xlab("Transcript count")+
    ylab("Number of times observed")+
    ggtitle(paste("Distribution of transcript counts of ", gene_oi_realname,'\nZero count = ',toString(zero_count))) +
    theme(text = element_text(size=TEXTSIZE),
          axis.text = element_text(size=TEXTSIZE),
          plot.title = element_text(size=TEXTSIZE))+
    scale_x_continuous(labels = comma)+
    scale_y_continuous(labels = comma)
  
  return(p)
  
}


# takes dataframe with raw transcript count observations as input (field name count)
# and condition field, takes additional parameter that holds list of names matching the
# different condition values
plot_violins <- function(gene_counts_df, labelnames) {
  
  TEXTSIZE=15
  ggplot(gene_counts_df, aes(factor(condition), counts)) + 
    #geom_violin(aes(fill = condition))+
    geom_boxplot(aes(fill = condition))+
    scale_x_discrete(breaks=seq(1,length(labelnames)),
                     labels=labelnames)+
    ggtitle(paste("Distribution of transcript counts of ", gene_oi_realname,'\nZero count = ',toString(zero_count))) +
    xlab(element_blank())+ylab('Transcript count')+
    theme(legend.position="none",
          text = element_text(size=TEXTSIZE),
          axis.text = element_text(size=TEXTSIZE),
          plot.title = element_text(size=TEXTSIZE),
          legend.text = element_text(size=TEXTSIZE),
          axis.text.x = element_text(angle = 90, hjust = 1, size=TEXTSIZE))
  
}


plot_scatter_w_highlighted_clusters<-function(df_toplot,x,y,cluster_assignments,
                                              myxlabel,myylabel,mytitle,col_vector,mylegendposition="right",
                                              thexlim=NULL, theylim=NULL) {
  
  TEXTSIZE=15
  p<-ggplot(data=df_toplot)+
    geom_point(aes_string(x=x,y=y,color=cluster_assignments))+#, color=cell_123_array)+
    ggtitle(mytitle)+
    xlab(myxlabel)+ylab(myylabel)+
    scale_color_manual(values=col_vector)+
    theme(legend.position=mylegendposition,
      text = element_text(size=TEXTSIZE),
      axis.text = element_text(size=TEXTSIZE),
      plot.title = element_text(size=TEXTSIZE),
      legend.text = element_text(size=TEXTSIZE))
  
  if (!is.null(thexlim)) { p<-p+xlim(thexlim) }
  if (!is.null(theylim)) { p<-p+ylim(theylim) }
      
  return(p)
  
}


  

plot_scatter_w_highlighted_clusters_condition<-function(df_toplot,x,y,cluster_assignments,conditions,condition_names,condition_markers,
                                              myxlabel,myylabel,mytitle,col_vector) {
  
  TEXTSIZE=15
  p<-ggplot(data=df_toplot)+
    geom_point(aes_string(x=x,y=y,color=cluster_assignments,shape=conditions),size=2,alpha=.8)+#, color=cell_123_array)+
    ggtitle(mytitle)+
    xlab(myxlabel)+ylab(myylabel)+
    scale_color_manual(values=col_vector)+
    scale_shape_manual(values=condition_markers,labels=condition_names)+
    theme(#legend.position="none",
      text = element_text(size=TEXTSIZE),
      axis.text = element_text(size=TEXTSIZE),
      plot.title = element_text(size=TEXTSIZE),
      legend.text = element_text(size=TEXTSIZE))
  
  return(p)
  
}

plot_scatter_w_highlighted_clusters_condition_exprgrad<-function(df_toplot,x,y,cluster_assignments_varname,conditions_varname,condition_names,condition_markers,
                                                        myxlabel,myylabel,mytitle,col_vector,selected_gene_expression_varname,savelocation,mylimits=NULL) {
  
  TEXTSIZE=15
  
  # The gene expression plot
  p1<-ggplot(data=df_toplot)+
    geom_point(aes_string(x=x,y=y,size=selected_gene_expression_varname,color=selected_gene_expression_varname))+
    scale_color_gradient(low="grey70", high="black",limits=mylimits)+
    scale_size_continuous(range = c(1, 10))+#mylimits
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )+
    ggtitle(mytitle)+
    xlab(myxlabel)+ylab(myylabel)+
    theme(#legend.position="none",
      text = element_text(size=TEXTSIZE),
      axis.text = element_text(size=TEXTSIZE),
      plot.title = element_text(size=TEXTSIZE),
      legend.text = element_text(size=TEXTSIZE))
  
  # The usual clustered scatter plot
  p2<-ggplot(data=df_toplot)+
    geom_point(aes_string(x=x,y=y,color=cluster_assignments_varname,shape=conditions_varname),size=2,alpha=.8)+#, color=cell_123_array)+
    scale_shape_manual(values=condition_markers,labels=condition_names)+
    scale_color_manual(values=col_vector)+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )+
    ggtitle(mytitle)+
    xlab(myxlabel)+ylab(myylabel)+
    theme(#legend.position="none",
      text = element_text(size=TEXTSIZE),
      axis.text = element_text(size=TEXTSIZE),
      plot.title = element_text(size=TEXTSIZE),
      legend.text = element_text(size=TEXTSIZE))
  
  # Combined plot needs to be saved while writing it
  # https://stackoverflow.com/questions/29708821/how-to-save-a-grid-plot-in-r
  pdf(savelocation, height = 6, width = 6, paper = "special")
  
  # Create the combined plot
  # https://stackoverflow.com/questions/11508902/plotting-discrete-and-continuous-scales-in-same-ggplot
  grid.newpage()
  pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
  print( p1 + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
  print( p2 + theme(legend.position="none") , vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
  
  # end of saving
  dev.off()
  
  return(list(p1,p2))
  
}

plot_scatter_gene_expression<-function(df_toplot,x,y,conditions_varname,condition_names,condition_markers,
                                       myxlabel,myylabel,mytitle,col_vector_grad=c('#91F3C8','#242122'),selected_gene_expression_varname,
                                       mylegendposition="right", mylimits=NULL,pointsize=2,
                                       cluster_outline=F,cl_str=NULL,mytextsize=15) {
  
  #colorScheme<-brewer.pal(3,"Spectral")

  # The gene expression plot
  p1<-ggplot(data=df_toplot)+
    geom_point(aes_string(x=x,y=y,
                          color=selected_gene_expression_varname, 
                          shape=conditions_varname
                          ),
               size=pointsize,alpha=0.8)+
    scale_color_gradientn(colours = col_vector_grad, limits=mylimits,oob=squish)+
    #scale_color_distiller(palette="Spectral", limits=mylimits)+
    #scale_color_gradient2(low=colorScheme[1], mid=colorScheme[2],  high=colorScheme[3], guide=guide_colorbar(barheight=15))+
    ggtitle(mytitle)+
    xlab(myxlabel)+ylab(myylabel)+
    theme(legend.position=mylegendposition,
      text = element_text(size=mytextsize),
      axis.text = element_text(size=mytextsize),
      plot.title = element_text(size=mytextsize),
      legend.text = element_text(size=mytextsize))
  
  if (!is.null(condition_markers)){
      p1<-p1+scale_shape_manual(values=condition_markers,labels=condition_names)
    }
  
  # Create cluster outline if desired
  if (cluster_outline) {
    for (ii in levels(df_toplot[[cl_str]])) {
        p1<-p1+geom_encircle(data=df_toplot[df_toplot[[cl_str]]==ii,],
            aes_string(x=x,y=y),expand=0)
    }
  }
  
  return(p1)
  
}

pump_out_freq_df<-function(p,gene_expression, color, shape) {
  NR_BINS=30 # 40
  
  # determine 80% percentile
  ordered_expression<-gene_expression[gene_expression>0] # remove zero-values
  ordered_expression<-ordered_expression[order(ordered_expression, decreasing=TRUE)]
  percentile_value<-ordered_expression[round(length(ordered_expression)*.8)] # currently not used
  determined_max_value<-max(ordered_expression)
  
  #mybinwidth=ceiling((max(gene_expression)+1)/NR_BINS) # raceid3 has much lower normalized values
  mybinwidth=(determined_max_value)/NR_BINS
  
  freq <- hist(x=as.numeric(gene_expression),
               breaks=seq(0,max(gene_expression)+mybinwidth,mybinwidth),
               plot=FALSE)
  
  freq_df <- data.frame(centers = freq$mids, counts = freq$counts, my_gene_nr=as.factor(ii))
  
  #geom_line(data=freq_df,
  #             stat="identity", 
  #             mapping=aes(x=centers, y=counts),
  #             color=color)
  #shape
  
}

my_title_row<-function(mytitle) {
  par(mar = c(0,0,0,0))
  p<-plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(mytitle), 
       cex = 1.6, col = "black", srt=90)
  return(p)
}

my_title_col<-function(mytitle) {
  par(mar = c(0,0,0,0))
  p<-plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("My title"), 
       cex = 1.6, col = "black", srt=0)
  return(p)
}

barplot_differential_expression<-function(selected_data_df,
                                          centers_varname,
                                          differential_expression_varname,
                                          all_gene_names,
                                          lowcol,highcol,
                                          ylabtext, mytitle) {
  TEXTSIZE=15
  
  mybreaks<-selected_data_df[[centers_varname]]
  
  ggplot(data=selected_data_df, mapping=aes_string(x=centers_varname, y=differential_expression_varname))+#,fill=dataset_id)) +
    geom_bar(stat="identity", mapping=aes_string(fill=differential_expression_varname))+
    scale_x_discrete(breaks=mybreaks,
                     labels=all_gene_names[selected_data_df$original_nr])+
    coord_flip()+
    xlab("Genes")+
    ylab(ylabtext)+
    scale_fill_gradient(low=lowcol, high=highcol)+
    ggtitle(mytitle)+
    theme(legend.position="none",
          text = element_text(size=TEXTSIZE),
          axis.text = element_text(size=TEXTSIZE),
          plot.title = element_text(size=TEXTSIZE))
}

barplot_differential_expression_v2<-function(selected_data_df,
                                          differential_expression_varname,
                                          center_varname,
                                          gene_name_varname,
                                          lowcol,highcol,
                                          ylabtext, mytitle) {
  TEXTSIZE=15
  
  mybreaks <- selected_data_df[[center_varname]]
  mylabels <- selected_data_df[[gene_name_varname]]
  
  ggplot(data=selected_data_df, mapping=aes_string(x=center_varname, y=differential_expression_varname))+#,fill=dataset_id)) +
    geom_bar(stat="identity", mapping=aes_string(fill=differential_expression_varname))+
    scale_x_discrete(breaks=mybreaks,
                     labels=mylabels)+
    coord_flip()+
    xlab("Genes")+
    ylab(ylabtext)+
    scale_fill_gradient(low=lowcol, high=highcol)+
    ggtitle(mytitle)+
    theme(legend.position="none",
          text = element_text(size=TEXTSIZE),
          axis.text = element_text(size=TEXTSIZE),
          plot.title = element_text(size=TEXTSIZE))
}

give_better_textsize_plot <- function(TEXTSIZE){
  theme(#legend.position="none",
        text = element_text(size=TEXTSIZE),
        axis.text = element_text(size=TEXTSIZE),
        plot.title = element_text(size=TEXTSIZE))
}


# TODO: this function needs to be finalized!!!!!!
plot_tsne_with_cluster_labels <- function(dataframe_cells, x_str, y_str, 
    cl_str, a_title, textsize=15, mylabelsize=3, pointsize=3,
    cluster_outline=F, custom_labels=NULL, myalpha=1, mycolors=NULL) {
  
  # determine cluster locations
  coord1_clusters <- numeric()
  coord2_clusters <- numeric()
  for (ii in seq(1,max(as.integer(levels(dataframe_cells[[cl_str]]))))) {
    coord1_clusters[[ii]] <- mean(dataframe_cells[[x_str]][dataframe_cells[[cl_str]]==ii])
    coord2_clusters[[ii]] <- mean(dataframe_cells[[y_str]][dataframe_cells[[cl_str]]==ii])
  }
  
  # Create dataframe with coordinates and labels
  if (!is.null(custom_labels)) { df_cluster_labels <- data.frame(coord1=coord1_clusters, coord2=coord2_clusters,clName=custom_labels) 
  } else { df_cluster_labels <- data.frame(coord1=coord1_clusters, coord2=coord2_clusters,clName=seq(1,max(as.integer(levels(dataframe_cells[[cl_str]]))))) }
        
  # Make the plot
  p <- ggplot(data=dataframe_cells)+
    geom_point(aes_string(x=x_str,y=y_str,color=cl_str),size = pointsize,alpha=myalpha)+#,shape=condition))+
    ggtitle(a_title)+
    give_better_textsize_plot(textsize)+
    #geom_label(data=df_cluster_labels,aes(x=coord1,y=coord2,label=clName),size=mylabelsize)+
    geom_text(data=df_cluster_labels,aes(x=coord1,y=coord2,label=clName),size=mylabelsize,fontface = 'bold')+
    theme(legend.position = 'none')
  
  if (!is.null(mycolors)) {
   p<-p+scale_color_manual(values = mycolors)   
  }
  
  # Create cluster outline if desired
  if (cluster_outline) {
    for (ii in levels(dataframe_cells[[cl_str]])) {
        p<-p+geom_encircle(data=dataframe_cells[dataframe_cells[[cl_str]]==ii,],
            aes_string(x=x_str,y=y_str),expand=0)
    }
  }
  
  return(p)
    
}

plot_tsne_with_cluster_labels2 <- function(dataframe_cells, x_str, y_str, 
    cl_str, a_title, textsize=15, mylabelsize=3, pointsize=3,
    cluster_outline=F, custom_labels=NULL, myalpha=1, mycolors=NULL) {
  
  # determine cluster locations
  coord1_clusters <- numeric()
  coord2_clusters <- numeric()

  count=0
  for (lvl in levels(dataframe_cells[[cl_str]])) {
    count=count+1
    coord1_clusters[[count]] <- mean(dataframe_cells[[x_str]][dataframe_cells[[cl_str]]==lvl])
    coord2_clusters[[count]] <- mean(dataframe_cells[[y_str]][dataframe_cells[[cl_str]]==lvl])
  }
  
  # Create dataframe with coordinates and labels
  df_cluster_labels <- data.frame(coord1=coord1_clusters, coord2=coord2_clusters,clName=levels(dataframe_cells[[cl_str]]))
        
  # Make the plot
  p <- ggplot(data=dataframe_cells)+
    geom_point(aes_string(x=x_str,y=y_str,color=cl_str),size = pointsize,alpha=myalpha)+#,shape=condition))+
    ggtitle(a_title)+
    give_better_textsize_plot(textsize)+
    #geom_label(data=df_cluster_labels,aes(x=coord1,y=coord2,label=clName),size=mylabelsize)+
    geom_text(data=df_cluster_labels,aes(x=coord1,y=coord2,label=clName),size=mylabelsize,fontface = 'bold')+
    theme(legend.position = 'none')
  
  if (!is.null(mycolors)) {
   p<-p+scale_color_manual(values = mycolors)   
  }
  
  # Create cluster outline if desired
  if (cluster_outline) {
    for (ii in levels(dataframe_cells[[cl_str]])) {
        p<-p+geom_encircle(data=dataframe_cells[dataframe_cells[[cl_str]]==ii,],
            aes_string(x=x_str,y=y_str),expand=0)
    }
  }
  
  return(p)
    
}

get_cluster_centers <- function(df_tsne, varname1='V1', varname2='V2', clustervarname='cluster_assignments') {

  # determine cluster centers by looping over cluster indices, and taking means
  coord1_clusters <- numeric()
  coord2_clusters <- numeric()
  #for (ii in seq(1,max(as.integer(df_tsne[[clustervarname]])))) {
  for (ii in seq(1,max(as.integer(levels(as.factor(df_tsne[[clustervarname]])))))) {
    coord1_clusters[[ii]] <- mean(df_tsne[[varname1]][df_tsne[[clustervarname]]==ii])
    coord2_clusters[[ii]] <- mean(df_tsne[[varname2]][df_tsne[[clustervarname]]==ii])
  }
  
  # store in dataframe
  df_cluster_labels <- data.frame(coord1=coord1_clusters, 
                                  coord2=coord2_clusters,
                                  clName=seq(1,max(as.integer(levels(as.factor(df_tsne[[clustervarname]]))))))
  
  # return the locations in a dataframe
  return(df_cluster_labels)

}

# function to convert cut rownames to centers of bins
# not sure if this is the best solution but it works
# code from
# https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
# example:
# get_midpoints_from_cut(levels(bin_assignment))
get_midpoints_from_cut <- function(x){#, dp=2){
    lower <- as.numeric(gsub(',.*','',gsub('\\(|\\[|\\)|\\]','', x)))
    upper <- as.numeric(gsub('.*,','',gsub('\\(|\\[|\\)|\\]','', x)))
    #return(round(lower+(upper-lower)/2, dp))
    return(lower+(upper-lower)/2)
}


return_short_name <- function(gene_name) {
    
    substr(gene_name,1,(strfind(gene_name,'_')[1]-1))
}
    







plot_vulcano_correlations_for_gene <- function() {
    
    stop('Function script not finished yet.. Script needs to be edited.')
    
    # Make a vulcano plot for desired genes
    
    GOI_NAMES<-c('NGF', 'GDNF', 'NTF3', 'SEMA3A')
    corr<-list()
    corr_filt<-list()
    corr_filt_30<-list()
    
    PVAL_SHOW <- 1/10^8
    PVAL_CUTOFF <- 1/100
    
    i=1
    for (i in 1:length(GOI_NAMES)){
    
        CURRENT_GENE_NAME <- GOI_NAMES[[i]]
        
        # gather statistics
        corr[[i]] <- get_gene_correlations(expression_matrix_sel,CURRENT_GENE_NAME,show_hist=F)
            # note that genes expressed in few cells can create weird correlations
        # add nr cells in which each gene is found
        corr[[i]] <- cbind(corr[[i]],
                        data.frame(
                            nr_cells=apply(1*(expression_matrix_sel>0),1,sum),
                            short_name = gsub("__chr(\\d+|[MXY])", '', corr[[i]]$gene_name),
                            is_TF = TF_search_out$boolean_names,
                            is_L  = ligand_search_out$boolean_names,
                            R_associated = R_associated_$R_associated_filled_str))
        
        # identify in how many cells goi is found
        nr_cells_goi <- corr[[i]][corr[[i]]$self,]$nr_cells
        
        # now make subselection of significant ones, only take top 30
        #corr_filt[[i]]<-corr[[i]][corr[[i]]$pval_adjusted<PVAL_SHOW,]
        corr_filt[[i]]<-corr[[i]][corr[[i]]$pval_adjusted<PVAL_CUTOFF,]
        corr_filt_30[[i]]<-corr_filt[[i]][order(corr_filt[[i]]$pval,decreasing = F),][1:min(dim(corr_filt[[i]])[1],30),]
        
        # Now create a vulcano plot
        thenames<-gsub("__chr(\\d+|[MXY])", '', corr[[i]]$gene_name)
        p<-ggplot(corr[[i]][!corr[[i]]$self,])+
            geom_point(aes(x=correlations, y=-log10(pval_adjusted)))+
            geom_point(data=corr_filt[[i]][!corr_filt[[i]]$self,],aes(x=correlations, y=-log10(pval_adjusted)),color='red')+
            geom_text_repel(data=corr_filt_30[[i]][!corr_filt_30[[i]]$self,],
                aes(x=correlations, y=-log10(pval_adjusted)),label=corr_filt_30[[i]][!corr_filt_30[[i]]$self,]$short_name,
                size=5,color='black')+
            xlab('Correlation coefficient')+ylab('-log10(adjusted p_value)')+
            ggtitle(paste0('Correlation with ',corr[[i]]$short_name[corr[[i]]$self]))+
            give_better_textsize_plot(15)#+
            #xlim(-.5,.5)
        p
        #ggsave(paste0(outputDir_Iliana,'correlations/correlations_',CURRENT_GENE_NAME,'.png'),plot=p,
        #            units='mm',width=100,height=100,dpi=300)
        
        #===
        # Create another vulcano plot with ligands highlighted
        
        # Now create a vulcano plot
        thenames<-gsub("__chr(\\d+|[MXY])", '', corr[[i]]$gene_name)
        p<-ggplot(corr[[i]][!corr[[i]]$self,])+
            geom_point(aes(x=correlations, y=-log10(pval_adjusted)))+
            geom_point(data=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_L,],aes(x=correlations, y=-log10(pval_adjusted)),color='pink')+
            geom_text_repel(data=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_L,],
                aes(x=correlations, y=-log10(pval_adjusted)),label=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_L,]$short_name,
                size=3,color='purple')+
            xlab('Correlation coefficient')+ylab('-log10(adjusted p_value)')+
            ggtitle(paste0('Correlation with ',corr[[i]]$short_name[corr[[i]]$self]),' (ligands highlighted)')+
            give_better_textsize_plot(15)#+
            #xlim(-.5,.5)
        p
        #ggsave(paste0(outputDir_Iliana,'correlations/correlations_',CURRENT_GENE_NAME,'_ligands.png'),plot=p,
        #            units='mm',width=100,height=100,dpi=300)
        
        
        #===
        # Create another vulcano plot with transcription factors highlighted
        
        # Now create a vulcano plot
        thenames<-gsub("__chr(\\d+|[MXY])", '', corr[[i]]$gene_name)
        p<-ggplot(corr[[i]][!corr[[i]]$self,])+
            geom_point(aes(x=correlations, y=-log10(pval_adjusted)))+
            geom_point(data=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_TF,],aes(x=correlations, y=-log10(pval_adjusted)),color='orange')+
            geom_text_repel(data=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_TF,],
                aes(x=correlations, y=-log10(pval_adjusted)),label=corr_filt[[i]][!corr_filt[[i]]$self&corr_filt[[i]]$is_TF,]$short_name,
                size=3,color='red')+
            xlab('Correlation coefficient')+ylab('-log10(adjusted p_value)')+
            ggtitle(paste0('Correlation with ',corr[[i]]$short_name[corr[[i]]$self]),' (TFs highlighted)')+
            give_better_textsize_plot(15)#+
            #xlim(-.5,.5)
        p
        #ggsave(paste0(outputDir_Iliana,'correlations/correlations_',CURRENT_GENE_NAME,'_TFs.png'),plot=p,
        #            units='mm',width=100,height=100,dpi=300)
        
        #===
        # export xls files with this data
        #write.xlsx(corr[[i]], file=paste0(outputDir_Iliana,'correlations/correlations_',CURRENT_GENE_NAME,'.xlsx'))#, sheetName=paste0('cluster',cl_idx))
    
    }
}









# gene_list: list of genes you want to have the expression plotted for on tSNE
# HSMM_filtered: Monocle dataframe 
plot_tSNE_w_expression <- function(config, HSMM_filtered, gene_list, 
                            pointsize=.35, mytextsize=8, mywidth=32, myheight=32) {
    
    expression_genes_normalized <- Matrix::t(Matrix::t(exprs(HSMM_filtered))/sizeFactors(HSMM_filtered))
    
    df_Mcle_cells_exprs <- data.frame(
        tSNE1=reducedDimA(HSMM_filtered)[1,],
        tSNE2=reducedDimA(HSMM_filtered)[2,],
        condition=HSMM_filtered$condition,
        cluster=HSMM_filtered$Cluster)
    
    # Now repeat previous code, which adds the gene expression to the 
    # cell dataframe, and then plots it
    for (idx in 1:length(gene_list)) {
        
        # CHECK WHETHER THIS WORKS
        goi_idx <- give_idxs_on_gene_list(gene_list[idx],rownames(HSMM_filtered))$boolean_list
        expr_goi<-expression_genes_normalized[goi_idx,]
        df_Mcle_cells_exprs[gene_list[idx]]<-expr_goi
    
        min<-expr_goi[order(expr_goi,decreasing = F)][ceiling(length(expr_goi)*0.01)]
        max<-expr_goi[order(expr_goi,decreasing = F)][ceiling(length(expr_goi)*0.99)]
        mylimits<-c(min,max)
        
        p<-plot_scatter_gene_expression(df_Mcle_cells_exprs,'tSNE1','tSNE2',
                                            'condition',config$shortdatasetnames,condition_markers = c(16,16),
                                           element_blank(),element_blank(),gene_list[idx],
                                           col_vector_grad=col_viridius_inferno_inv,
                                           selected_gene_expression_varname=gene_list[idx],
                                           mylegendposition="none", mylimits=mylimits,
                                           pointsize=pointsize,mytextsize=mytextsize,
                                            cluster_outline=T,cl_str='cluster')
        
        #p <- plot_tsne_with_cluster_labels(df_Mcle_cells_exprs,'tSNE1','tSNE2',marker_list[idx],'Cluster assignment',
        #            textsize=8,mylabelsize=3,pointsize=.5,cluster_outline = F)
        p<-p+theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
                    axis.text.x = element_blank(),axis.text.y = element_blank(),
                    plot.title=element_text(margin = margin(b=0, t = 0), hjust=0))+
                    monocle_theme_opts()
                    # inside = b=-20, t = 10
        print(p)
        
        # Save it
        ggsave(paste0(config$outputDir,'tSNE_expression_',gene_list[idx],'.svg'),plot=p,
            units='mm',width=mywidth,height=myheight,dpi=300)
        ggsave(paste0(config$outputDir,'tSNE_expression_',gene_list[idx],'.png'),plot=p,
            units='mm',width=mywidth,height=myheight,dpi=300)
    }

}

calc_limits <- function(values, percentile=.01) {
    
        min<-values[order(values,decreasing = F)][ceiling(length(values)*percentile)]
        max<-values[order(values,decreasing = F)][ceiling(length(values)*(1-percentile))]
        mylimits<-c(min,max)
    
}








# Here is a convenient way to plot expression on tSNE maps
plotgeneinxy <- function(name_goi, expr_normalized, comp_1, comp_2, print_yes=T) {
    #name_goi = 'KCNQ1OT1'
    goi_idx <- give_idxs_on_gene_list(name_goi,rownames(expr_normalized))$boolean_list
    if (!any(goi_idx)) { warning(paste0('Gene ',name_goi,' not found.')); return(list(p=NA,expr_goi=NA)) }
    expr_goi<-expr_normalized[goi_idx,]
    #expr_goi<-1*(expr_goi>50)
    #
    #MYPCA1=1; MYPCA2=2
    plot_df=data.frame(comp_1=comp_1, comp_2=comp_2,
            #cluster=HSMM_filtered_current$Cluster, condition=HSMM_filtered_current$condition,
            expr_goi=expr_goi)
    p<-ggplot(plot_df)+
        geom_point(aes(x=comp_1,y=comp_2,color=expr_goi))+
        #geom_point(data=filter(plot_df,expr_goi>75),aes(x=x,y=y),color='turquoise1')+
        xlim(calc_limits(comp_1))+ylim(calc_limits(comp_2))+
        #xlab(paste0('PCA-',MYPCA1))+ylab(paste0('PCA-',MYPCA2))+
        scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(expr_goi),oob=squish)+
        give_better_textsize_plot(10)+ggtitle(name_goi)
    if (print_yes) { print(p) }
    # hist(expr_goi,1000,xlim=c(0,500))
    return(list(p,expr_goi))
}

plot_sth_in_XY = function(comp_1, comp_2, color_by, name_for_color, print_yes=T,mypointsize=1,stringent_lims=F,custom_colors=NULL) {
    plot_df=data.frame(comp_1=comp_1, comp_2=comp_2)
    plot_df[[name_for_color]]=color_by
    p<-ggplot(plot_df)+
        geom_point(aes_string(x='comp_1',y='comp_2',color=name_for_color),size=mypointsize)+
        #geom_point(data=filter(plot_df,expr_goi>75),aes(x=x,y=y),color='turquoise1')+
        give_better_textsize_plot(20)
    if (!is.null(custom_colors)){
        p=p+scale_color_manual(breaks = custom_colors$breaks,
                               values= custom_colors$values)
    }
    if (stringent_lims) {p=p+xlim(calc_limits(comp_1))+ylim(calc_limits(comp_2))}
        #xlab(paste0('PCA-',MYPCA1))+ylab(paste0('PCA-',MYPCA2))+
        #scale_color_gradientn(colours = col_viridius_inferno_inv, limits=calc_limits(expr_goi),oob=squish)
    if (print_yes) {print(p)}
    return(p)
}



# Function to plot tSNE maps (and color by condition or cluster)
plot_Race_tSNE = function(output_object, config, color_by = 'Race_cluster', sel_source=1:99,
                        Expr=NULL) {
    Race_tsne1 = output_object$RaceID_object@tsne[,1]
    Race_tsne2 = output_object$RaceID_object@tsne[,2]
    Race_cluster   = as.factor(output_object$RaceID_object@cluster$kpart)
    Race_cluster2  = as.factor(output_object$RaceID_object@cpart)

    Plate_Origin = dt_pooled$origin[names(output_object$RaceID_object@cpart)]
    Condition = as.factor(config$conversion_numbers_int[Plate_Origin])

    df_toplot = data.frame(tSNE_1=Race_tsne1, tSNE_2=Race_tsne2,
                        Race_cluster=Race_cluster,
                        Condition=Condition,
                        Plate_Origin=Plate_Origin)
    
    if (!is.null(Expr)) { df_toplot$Expr = Expr; if (!is.factor(Expr)) {df_toplot$Expr_Log10 = log10(Expr)}}
    
    p=ggplot(df_toplot)+
        geom_point(aes_string(x='tSNE_1',y='tSNE_2',color=color_by),size=.5)+
        theme_bw()+
        give_better_textsize_plot(10)
        
    #print(p)
    return(p)
}




    ############################################################


createVennTwoLists = function(list1, list2, name1='list1', name2='list2') {
    
    # Create overlap matrix
    all_genes = unique(c(list1, list2))
    gene_overlap_df=as.matrix(1*data.frame(list1=all_genes %in% list1, 
                                           list2=all_genes %in% list2))
    
    # Now plot this with Venn
    
    # library 1
    # library("limma")
    my_ven <- vennCounts(gene_overlap_df)
    vennDiagram(my_ven, cex=1, names = c(name1,name2)) # note: cex is general text scaling parameter
    
    # library 2
    #library(venneuler)
    #v<-venneuler::venneuler(gene_overlap_df)
    #v$labels=c(paste0(name1,'\n(',sum(gene_overlap_df[,1]),')'),
    #           paste0(name2,'\n(',sum(gene_overlap_df[,2]),')'))
    #par(cex = 1.4) 
    #p=plot(v, main = paste0("Overlap\n",sum(gene_overlap_df[,1]&gene_overlap_df[,2])), cex.main = 1)

    #return(p)
    
}
