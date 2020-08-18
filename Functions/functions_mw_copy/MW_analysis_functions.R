
######################################################################
# libraries
######################################################################

# Libraries required for GO analysis
"
# Install if necessary
if('GOstats' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('GOstats')}
if('org.Mm.eg.db' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('org.Mm.eg.db')}
if('org.Hs.eg.db' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('org.Hs.eg.db')}
if('GSEABase' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('GSEABase')}
if('biomaRt' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('biomaRt')}
"
library("GOstats")
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("GSEABase")
library("biomaRt")  # BiocManager::install("biomaRt")

library(gridExtra) # required for grid.arrange

library(DESeq2) # BiocManager::install("DESeq2")

# org.Hs.egGO

######################################################################

# gene_names<-rownames(groupedSCS$patientAll@ndata)
get_idx_for_gene_name<-function(gene_names,name_str,showmessages=T) {
  
  # find it
  gene_idx<-which(grepl(name_str,gene_names))
  
  # do some checks
  if (length(gene_idx)==0){
    if (showmessages) { print(paste(name_str, ": No gene found with that name")) }
    gene_idx=NaN
    # next
  } else if (length(gene_idx)>1) {
    if (showmessages) { print(paste(name_str, ": Multiple genes found with that name")) } 
    gene_idx=NaN
    # next
  } else {
    if (showmessages) { print(paste(gene_names[gene_idx], ' found gene with that name; index= ', gene_idx, '.'),sep="") }
  }
  
  # return index
  return(gene_idx)

}

# all_gene_expression needs to have rownames that give the gene names
get_expression_gene <- function(all_gene_expression, name_str) {
  
  # get index for gene of interest based search string
  gene_idx<-get_idx_for_gene_name(gene_names=rownames(all_gene_expression),name_str)
  
  # if found
  if (!is.nan(gene_idx)){
    
    # select the expression of the gene from the table
    the_gene_expression <- all_gene_expression[gene_idx,]
    
  } else {
   
    # throw error if gene not found
    stop('Error, couldn\'t find gene.')
     
  }
    
  return(the_gene_expression)
  
}

# very simple expression lookup
get_expression_gene2 = function(expression_matrix, gene_query) {
    
    sel_matrix = 
        expression_matrix[ gsub("__chr(\\d+|[MXY])", '', rownames(expression_matrix)) %in% 
                            gsub("__chr(\\d+|[MXY])", '', gene_query),]
    
    return(sel_matrix)
    
}

# shorthand
strip__chrXX = function(gene_query) {
    return(gsub("__chr(\\d+|[MXY])", '', gene_query))
}

find__chrXX = function(gene_query, gene_names) {

    # perform the search (do one-by-one to preserve order) -- TODO might be better way    
    search_result = unlist(lapply(gene_query,function(X) {which(give_idxs_on_gene_list(X, gene_names)$boolean_list)}))

    # generate output
    output_gene_names = gene_names[search_result]
    
}



## get_differential_gene_expression()
# calculates fold-change as set1/set2
# indices_set1 : indices of columns (cells) that should be selected for set 1
# indices_set2 : indices of columns (cells) that should be selected for set 2
# all_gene_expression : matrix with all gene expression, normalized as grun/oudenaarden (raceID2)
# all_gene_expression_raw : same as previous, but raw transcript counts
# method='race_fixed' : pertains to rescaling of expression values; see comments in function
# pcutoff=0.01 : cutoff value for p values, entries above will be filtered out
# -- written by MW, 03-2019
get_differential_gene_expression <- function(indices_set1,indices_set2,all_gene_expression,all_gene_expression_raw,method='race_fixed',pcutoff=0.01) {
  
  # get names of genes
  gene_names <- rownames(all_gene_expression)
  
  # Get the gene expression of the two subsets
  gene_expression_set1  <- all_gene_expression[,indices_set1]
  gene_expression_set2  <- all_gene_expression[,indices_set2]
  
  # Calculate mean and median expression values
  mean_gene_expression_set1   <- rowMeans(as.matrix(gene_expression_set1))
  mean_gene_expression_set2   <- rowMeans(as.matrix(gene_expression_set2))
  median_gene_expression_set1 <- rowMedians(as.matrix(gene_expression_set1))
  median_gene_expression_set2 <- rowMedians(as.matrix(gene_expression_set2))
  
  # Calculate standard deviations
  sd_gene_expression_set1 <- rowSds(as.matrix(gene_expression_set1))
  sd_gene_expression_set2 <- rowSds(as.matrix(gene_expression_set2))
  
  # Now also calculate counts of positive in cells (similar to done before)
  detectioncount_gene_expression_set1<-rowSums(1*(gene_expression_set1>0.1))
  detectioncount_gene_expression_set2<-rowSums(1*(gene_expression_set2>0.1))
  
  # Calculate rescaling factor to rescale data to get cell total transcript counts to be rescaled
  # to the set median or minimum.
  # Calculate rescaling factors as done Grun & Van Oudenaarden (no pro'lly stands for normalization)
  # IMPORTANT NOTE: this value critically depends on whether "min" or "median" is used for normalization;
  # this is done in the raceid code, where ndata is calculated
  #
  # chose rescaling factor to continue with
  if (identical(method,'min')) {
    no_set1 <- (median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (min(apply(all_gene_expression_raw,2,sum)))
    no_set2 <- (median(apply(all_gene_expression_raw[,indices_set2],2,sum))) / (min(apply(all_gene_expression_raw,2,sum)))
  } else if (identical(method,'min_same')) {
    no_set1 <- (median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (min(apply(all_gene_expression_raw,2,sum)))
    no_set2 <- no_set1
  } else if (identical(method,'median')) {
    no_set1 <- (median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (median(apply(all_gene_expression_raw,2,sum)))
    no_set2 <- (median(apply(all_gene_expression_raw[,indices_set2],2,sum))) / (median(apply(all_gene_expression_raw,2,sum)))
  } else if (identical(method,'race_12')){
    # Rescaling both with separate scale
    no_set1<-(median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (min(apply(all_gene_expression_raw+.1,2,sum)))
    no_set2<-(median(apply(all_gene_expression_raw[,indices_set2],2,sum))) / (min(apply(all_gene_expression_raw+.1,2,sum)))
  } else if (identical(method,'race_bug')){
    # This is what raceID does, but calculted in a more straightforward way
    # But also with 0.1 bug
    no_set1<-(median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (min(apply(all_gene_expression_raw+.1,2,sum)))
    no_set2<-no_set1
  } else if (identical(method,'race_fixed')){
    # This is what raceID does, but calculted in a more straightforward way
    # And without 0.1 bug
    no_set1<-(median(apply(all_gene_expression_raw[,indices_set1],2,sum))) / (min(apply(all_gene_expression_raw,2,sum)))
    no_set2<-no_set1
  } else if (identical(method,'race_identical')){
    # entirely identical to raceID2 (including +0.1 bug)
    no_set1<-if ( sum(indices_set1) > 1 ) median(apply(all_gene_expression_raw[,indices_set1],2,sum))/median(apply(all_gene_expression[,indices_set1],2,sum)) else sum(all_gene_expression_raw[,indices_set1])/sum(all_gene_expression[,indices_set1])
    no_set2<-no_set1
  } else if (identical(method,'none')) {
    no_set1<-1
    no_set2<-1
  } else {
    stop('Invalid method supplied.')
  }
  ratio_rescaling_factors <- no_set1/no_set2 # interesting to see whole-set expression difference (not used)
  
  # Now rescale them accordingly
  mean_gene_expression_set1_rescaled <- no_set1*mean_gene_expression_set1
  mean_gene_expression_set2_rescaled <- no_set2*mean_gene_expression_set2
  
  # Now count the differential expression
  differential_gene_expression <- mean_gene_expression_set1_rescaled / mean_gene_expression_set2_rescaled
  
  # Now calculate standard deviations
  differential_gene_expression_stdev <- differential_gene_expression*
    sqrt((sd_gene_expression_set1*abs(no_set1)/mean_gene_expression_set1_rescaled)^2+
           (sd_gene_expression_set2*abs(no_set2)/mean_gene_expression_set2_rescaled)^2)
  
  # now calculate p-values
  pv <- binompval(mean_gene_expression_set2_rescaled/sum(mean_gene_expression_set2_rescaled),
                  sum(mean_gene_expression_set1_rescaled),mean_gene_expression_set1_rescaled)
  
  # create dataframe similar grun & van oudenaarden
  diff_expr_df <- data.frame(mean.set2=mean_gene_expression_set2_rescaled,mean.set1=mean_gene_expression_set1_rescaled,
                             fc=differential_gene_expression,
                             fc_inv=1/differential_gene_expression,
                             pv=pv,
                             stdev=differential_gene_expression_stdev,
                             cellcount.set1=detectioncount_gene_expression_set1,
                             cellcount.set2=detectioncount_gene_expression_set2,
                             row.names = gene_names,
                             gene_name = gene_names,
                             stringsAsFactors = FALSE)
  
  # filter based on p values, order by fold-change, add row line number (useful bar plot)
  diff_expr_df_filterpv <- diff_expr_df[diff_expr_df$pv<pcutoff,]
  diff_expr_df_filterpv <- diff_expr_df_filterpv[order(diff_expr_df_filterpv$fc,decreasing=T),]
  diff_expr_df_filterpv <- mutate(diff_expr_df_filterpv,n123=factor(1:nrow(diff_expr_df_filterpv)))
  
  return(list(diff_expr_df,diff_expr_df_filterpv))  
  
}


# calculates a factor vector that maps the combined cellular data to the conditions
# required input: config (needed for the conversion numbers)
get_condition_factors <- function(config,thedata) {
  
  names_to_convert_factors = colnames(thedata)
  condition_factors <- factor(rep(0,length(names_to_convert_factors)), 
                              levels= seq(0,NR_CONDITIONS))
  for (ii in seq(1,length(config$conversion_numbers))) {
    hit_idxs<-which(grepl(config$conversion_searchterms[ii],names_to_convert_factors))
    condition_factors[hit_idxs]<-config$conversion_numbers[ii]
  }
  
  return(condition_factors)
}


# this function converts a list of wellnames (generally the col names of a transcript
# count matrix) to a list of factors that match the condition
# example:
# my_conditions <- convert_wellname_to_factor(colnames(count_matrix), config$conversion_searchterms, config$conversion_numbers) 
convert_wellname_to_factor <- function(wellnames, searchterms, conditions_as_factors) {
    
    well_conditions <- factor(rep(0, length(wellnames)), levels=levels(conditions_as_factors))
    
    for (ii in seq(length(searchterms))) {
        current_searchterm <- searchterms[[ii]]
        well_conditions[grepl(current_searchterm,wellnames)] <- conditions_as_factors[[ii]]
    }
    
    return(well_conditions)
    
}

give_partialgenenames_idxs <- function(partial_name_list, gene_names) {
    
    hit_list_idxs <- integer()
    for (idx in seq(1,length(partial_name_list))) {
      
      hit_idx <- which(grepl(partial_name_list[[idx]], gene_names))
      
      if (length(hit_idx)==1) {
        hit_list_idxs[idx]<-hit_idx
      } else {
        hit_list_idxs[idx]<-NA
      }
      
    }
    hits_list_idxs_clean<-hit_list_idxs[!is.na(hit_list_idxs)]
    
    return(hits_list_idxs_clean)
    
}


# e.g.
# df_cells <- return_convenient_dataframe_cells(groupedSCS,  group_oi = 'patientAll', config)
return_convenient_dataframe_cells <- function(groupedSCS, group_oi, config) {
    
    dataset_annotations <- convert_wellname_to_factor(colnames(groupedSCS[[group_oi]]@ndata), config$conversion_searchterms, config$conversion_numbers)
    
    df_cells <-
    data.frame(V1=groupedSCS[[group_oi]]@tsne$V1,
               V2=groupedSCS[[group_oi]]@tsne$V2,
               dataset_annotations=dataset_annotations,
               cluster_assignment=as.factor(groupedSCS[[group_oi]]@cpart))

}


search_gene_for_me <- function(searchterm, gene_names){
    
    found_names<-gene_names[which(grepl(searchterm,gene_names))]
    found_gene_idxs<-which(grepl(searchterm,gene_names))

    return(list(found_names,found_gene_idxs))
       
}


add_goi_expression_to_df<-function(df_cells, gene_str, all_gene_expression){
 
    goi_expression<-get_expression_gene(all_gene_expression = all_gene_expression, name_str = gene_str)
    
    df_cells[gene_str]<-goi_expression
       
    return(df_cells)
    
}







# Look for genes with high differential expression by looking for
# non-uniformity of expression
#
# Important arguments:
# expression_minimum: the # of cells where the 
#
# example:
# find_pattern_in_cells2dplot(gene_expression_matrix=normalized_Mcle_expr, 
#                            cellpositionsx=df_tsne_3sets$tSNE1_Mcle_he,
#                            cellpositionsy=df_tsne_3sets$tSNE2_Mcle_he,
#                            xbins=15,
#                            ybins=15,
#                            expression_minimum=.2)
find_pattern_in_cells2dplot <- function(gene_expression_matrix, cellpositionsx, cellpositionsy,
                                        xbins, ybins, expression_minimum) {

    
    # per gene in how many cell there is any signal >0
    gene_expr_nr_cells <- apply(gene_expression_matrix>0,1,sum)
    # now select for genes with signal in at least 20% of cells
    my_gene_selection<-which(gene_expr_nr_cells>round(dim(gene_expression_matrix)[2]*expression_minimum))
    
    # now loop over that selection
    xgroup <- cut(cellpositionsx, xbins, include.lowest=TRUE, labels=F)
    ygroup <- cut(cellpositionsy, ybins, include.lowest=TRUE, labels=F)
    # First calculate "heatmap" for cells
    cell_map=matrix(rep(0,max(xgroup)*max(ygroup)),nrow=max(ygroup))
    for (i in 1:length(goi_expression_2)) {
        cell_map[xgroup[i],ygroup[i]]<-cell_map[xgroup[i],ygroup[i]]+1
    }
    final_cell_map<-cell_map/sum(as.vector(cell_map),na.rm=T)
    #image(final_cell_map)
    #image(log(final_cell_map),main='log cells')
    
    # now for all genes
    df_diff_unbiased <- data.frame(gene_idf=my_gene_selection,diff_expr=rep(0,length(my_gene_selection)),
                                   idx_tmp=rep(0,length(my_gene_selection)))
    for (jj in 1:length(my_gene_selection)) {
    #for (jj in 1:10) {
    
        gene_idx <- my_gene_selection[jj]
        goi_expression <- gene_expression_matrix[gene_idx,]
        goi_expression_2 <- goi_expression-min(goi_expression)
        # now remove highest 1% of values to prevent outliers to dominate the pattern
        cutoff<-goi_expression_2[order(goi_expression_2,decreasing = F)][round(length(goi_expression_2)*.99)]
        goi_expression_2[goi_expression_2>cutoff]<-0
    
        # Now calculate for genes
        gene_expr_map=matrix(rep(0,max(xgroup)*max(ygroup)),nrow=max(ygroup))
        for (i in 1:length(goi_expression_2)) {
            #if (!is.na(goi_expression_2[i])){
            gene_expr_map[xgroup[i],ygroup[i]]<-gene_expr_map[xgroup[i],ygroup[i]]+goi_expression_2[i]
            #}
        }
        final_gene_expr_map<-gene_expr_map/sum(as.vector(gene_expr_map),na.rm=T)
        #image(gene_expr_map)
        #image(final_gene_expr_map,main='gene expr')
        
        # Now calculate some distance measure of the two:
        #image(abs(final_cell_map-final_gene_expr_map))
        #image(final_cell_map-final_gene_expr_map,main='difference')
        
        df_diff_unbiased$diff_expr[jj]<-sum(as.vector((final_cell_map-final_gene_expr_map)^2))
        df_diff_unbiased$idx_tmp[jj]<-jj
        
    }
    
    df_diff_unbiased_ordered<-df_diff_unbiased[order(df_diff_unbiased$diff_expr,decreasing = T),]
    
    # Convenient commands to proceed with this are:
    #df_diff_unbiased_ordered[1:10,]
    #gene_names[df_diff_unbiased_ordered[1:100,]$gene_idf]
    
    return(df_diff_unbiased_ordered)

}

# This functions returns indices of genes in gene_names that occur
# in the list gene_list
# 
# gene_names should be geneName__chrX
# gene_list should be list of geneName only
#
# TODO: this function is not quite optimally written
# using gsub, list could be converted first
# ---> NOTE this function is deprecated, use 
# give_idxs_on_gene_list instead
give_transcription_factors_idxs <- function(gene_names, gene_list) {
    
        
    # Now look at transcription factors only
    tf_list_idxs <- double()
    for (idx in seq(1,length(gene_list))) {
      
      tf_idx <- which(grepl(paste0('^',gene_list[[idx]],'_'), gene_names))
      #tf_idx <- which(grepl(gene_list[[idx]], gene_names))
      
      if (length(tf_idx)==1) {
        tf_list_idxs[idx]<-tf_idx
      } else {
        tf_list_idxs[idx]<-NA
      }
      
    }
    tf_list_idxs_clean<-tf_list_idxs[!is.na(tf_list_idxs)]
    
    if (exists('mw_data')) {
        percentage_hit <- sum(mw_data[tf_list_idxs_clean,])/sum(mw_data)
        print(paste0('Fraction genes on list = ',toString(percentage_hit)))
    }
    
    return(tf_list_idxs_clean)
    
}

# gene_list: e.g. list of transcription factors or ligands
give_idxs_on_gene_list <- function(gene_query, gene_list) {

    # remove suffixes (won't affect gene names without suffixes)
    gene_query_short <- gsub("__chr(\\d+|[MXY])", '', gene_query)
    gene_list_short  <- gsub("__chr(\\d+|[MXY])", '', gene_list)
    
    # now find the list of genes 
    boolean_query   <- !is.na(match(gene_query_short, gene_list_short))
    boolean_list    <- !is.na(match(gene_list_short, gene_query_short))
    
    # return both ways
    out<-list(boolean_query=boolean_query,boolean_list=boolean_list)
    return(out)
}
# just create a synonym because I can never find this fn
get_idxs_on_gene_list<-give_idxs_on_gene_list

# for ligand-receptor pairs, this is convenient:
get_R_associated <- function(ligand_receptor_pairs_df, ligand_search_out) {
    
    warning('It is more convenient to use the function get_X_associated instead of get_R_associated.')
    
    R_associated <- lapply(unique(ligand_receptor_pairs_df$ligand[ligand_search_out$boolean_list]),
            function(X){ as.character( ligand_receptor_pairs_df$receptor[X==ligand_receptor_pairs_df$ligand] ) })
    
    R_associated_filled <- rep(NA,length(ligand_search_out$boolean_query))
    R_associated_filled[which(ligand_search_out$boolean_query)]<-R_associated
    
    R_associated_filled_str <- sapply(R_associated_filled, function(X) {toString(X)})
    
    return(list(R_associated=R_associated,R_associated_filled=R_associated_filled,R_associated_filled_str=R_associated_filled_str))
    
}

# more general function, which finds receptors associated with ligands, or ligands associated with receptors
get_X_associated <- function(ligand_receptor_pairs_df, X_search_out, search_for='R') {
    
    if (search_for == 'L') {
        X_associated <- lapply(unique(ligand_receptor_pairs_df$ligand[X_search_out$boolean_list]),
            function(X){ as.character( ligand_receptor_pairs_df$receptor[X==ligand_receptor_pairs_df$ligand] ) })
    } else if (search_for == 'R') {
        X_associated <- lapply(unique(ligand_receptor_pairs_df$receptor[X_search_out$boolean_list]),
            function(X){ as.character( ligand_receptor_pairs_df$ligand[X==ligand_receptor_pairs_df$receptor] ) })
    } else { stop('Unknown option for search_for parameter..') }
        
    
    X_associated_filled <- rep(NA,length(X_search_out$boolean_query))
    X_associated_filled[which(X_search_out$boolean_query)]<-X_associated
    
    X_associated_filled_str <- sapply(X_associated_filled, function(X) {toString(X)})
    
    return(list(X_associated=X_associated,X_associated_filled=X_associated_filled,X_associated_filled_str=X_associated_filled_str))
    
}





# used by MW_correlate1vsall function below    
cor.test.wrapper<-function(x,y){
    out<-cor.test(x,y); return(c(cor=out$estimate,p.value=out$p.value))
}
# extended version of  cor.test.wrapper   
cor.test.wrapper.ext<-function(x,y,treshold=0){
    out<-cor.test(x,y)
    pos_vals1<-sum(x>treshold)
    pos_vals2<-sum(y>treshold)
    return(c(cor=out$estimate,p.value=out$p.value,pos_vals1=pos_vals1,pos_vals2=pos_vals2))
}

# This is a refactored version of earlier written correlation function
# it correlates one row i of a matrix against all other rows j.
# 
# It also reports for all pairs of rows i and j how many datapoints were non-zero
MW_correlate1vsall <- function(expression_matrix, gene_of_interest,treshold=0) {
    
    # convert to index if a string was given
    if (is.character(gene_of_interest)) {
        # convert
        gene_of_interest<-which(row.names(expression_matrix)==gene_of_interest)
        # but if no match, throw error
        if (length(gene_of_interest)==0) { stop('Gene name not found.') }
    }

    # Correlate goi vs. all using apply function
    out<-apply(
        as.matrix(expression_matrix)[-gene_of_interest,], 
        1, cor.test.wrapper.ext, 
        y=as.matrix(expression_matrix)[gene_of_interest,],
        treshold=treshold)
    
    # convert to dataframe, add gene names and return
    out<-as.data.frame(t(out))
    out['gene_name']<-rownames(as.matrix(expression_matrix)[-gene_of_interest,])
    return(out)
    
    # # Only correlations:
    # out<-apply(as.matrix(expression_matrix)[-gene_of_interest,], 1, cor, y=as.matrix(expression_matrix)[gene_of_interest,])
    
    # # Should be equivalent, but doesn't work???
    # expression_goi    <- as.matrix(expression_matrix[gene_of_interest,])
    # expression_other  <- as.matrix(expression_matrix[-gene_of_interest,])
    # out2<-apply(expression_other, 1, cor.test.wrapper, y=expression_goi)
    
    # # These were selection criteria previously used
    # I don't think it is a good idea to select based on gene expression
    # if (cutoff==-1) {
    #      idx_sel<-which(expression_goi>0.1 & expression_other>0.1)
    # } else if(cutoff<=1 & cutoff>0) {
    #      idx_sel<-which(expression_goi>0.1 & sum(as.numeric((expression_other>0.1)))/length(expression_other)>cutoff )
    # } else if (cutoff==0){
    #      idx_sel<-which(expression_goi>0.1)
    # } else {  stop()  }
}


# Gene Ontology function based on Joep Eding's function for GO terms enriched in clusters of cells
# See: https://github.com/vanrooij-lab/SingleCellSeq/blob/development/AnalysisFunctions.R
# Function name: analyseClusterGeneOntology
#
#
# Required libraries
# library(semtree) # not sure whether this should be loaded to use toTable
#   AnnotationDbi seems to be the package with toTable, which comes from bioconductor
library(GO.db)

# Requires the following params:
# config$species <- 'human'
# config$geneIdentifierType <- "external_gene_name_with_chr"
# Function to perform go analysis on a given set of gene names
give_GO_terms <- function(config, gene_names_SCS, genes_query, genes_background=NULL) {

    warning('I recommend using an updated versino of this function, MW_GOKEGG_analysis.R > analyzeGeneOntology_MW()')
    
    print('Obtaining GO terms')
    
 if(config$species == 'mouse') {
    GOTerms = AnnotationDbi::toTable(org.Mm.egGO)
    goFrameOrganism = "Mus musculus"
  } else if(config$species == 'human') {
    GOTerms = AnnotationDbi::toTable(org.Hs.egGO)
    goFrameOrganism = "Homo sapiens"
  }

    if(config$species == 'human') {
        mitoGeneIdentifier <- '^MT-'
    } else if(config$species == 'mouse') {
        mitoGeneIdentifier <- '^mt-'
    }
    
  if(exists('geneIdentifierConversionTable')) {
      print('Conversion table is present')
    geneNames <- geneIdentifierConversionTable
  } else {
      
      print('Conversion table not present; Initializing bioMart')
      
      # Initialize bioMart
      if(!(config$species %in% c('human','mouse'))) {
        stop("Incorrect setting for 'config$species', enter valid species name")
      } else if(config$species == 'human') {
        mart <- useDataset('hsapiens_gene_ensembl', useMart("ENSEMBL_MART_ENSEMBL",host="uswest.ensembl.org", ensemblRedirect = FALSE))
            # Note that the parameters host="uswest.ensembl.org", ensemblRedirect = FALSE are not necessary,
            # but apparently servers are sometimes not functioning properly, so this redirects to a server
            # that seems to work (most of the time)
            # (The non-descript error thrown is "biomaRt expected a character string of length 1.")
        
      } else if(config$species == 'mouse') {
        mart <- useDataset("mmusculus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="uswest.ensembl.org", ensemblRedirect = FALSE))  
            # Note that the parameters host="uswest.ensembl.org", ensemblRedirect = FALSE are not necessary,
            # but apparently servers are sometimes not functioning properly, so this redirects to a server
            # that seems to work (most of the time)
            # (The non-descript error thrown is "biomaRt expected a character string of length 1.")    
        
      }
      
    print('Getting gene names from Ensembl. This might take a while...')
    geneNames <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), # note: entrezgene_id used to be entrezgene
      filters = ifelse((config$geneIdentifierType == 'external_gene_name_with_chr'), 'external_gene_name', config$geneIdentifierType),
      values = if(config$geneIdentifierType == 'external_gene_name_with_chr') gsub("__chr(\\d+|[MXY])", '', gene_names_SCS) else gene_names_SCS,
      mart <- mart
    )
    assign('geneIdentifierConversionTable', geneNames, .GlobalEnv)
  }
    
  # Prepare goFrame containing associations between genes and GO terms
  if(!exists('goGSC')) {
      print("Preparing gene associations")
      goframeData = data.frame(GOTerms$go_id, GOTerms$Evidence, GOTerms$gene_id)
      goFrame=GOFrame(goframeData,organism=goFrameOrganism)
      goAllFrame=GOAllFrame(goFrame)
      goGSC <- GeneSetCollection(goAllFrame, setType = GOCollection())
      assign('goGSC', goGSC, .GlobalEnv)
  } else {
      print("Using goGSC that is already present")
  }
    
  # Now convert the gene names    
  genes_query_entrez      = geneNames$entrezgene[which(geneNames$external_gene_name %in% gsub("__chr(\\d+|[MXY])", '', genes_query))]
  if (!is.null(genes_background)) {
    genes_background_entrez = geneNames$entrezgene[which(geneNames$external_gene_name %in% gsub("__chr(\\d+|[MXY])", '', genes_background))]
  } else { genes_background_entrez <- NULL }
  
  print('Setting up analysis')
    # Set up analysis
    params <- GSEAGOHyperGParams(name="custom GSEA",
                                 geneSetCollection=goGSC,
                                 geneIds = genes_query_entrez,
                                 universeGeneIds = genes_background_entrez, # defines background
                                 ontology = 'BP',
                                 pvalueCutoff = 0.001,
                                 conditional = T,
                                 testDirection = 'over'
    )
    
    print('Running analysis')
    # Run enrichment analysis
    hyperGResult <- hyperGTest(params)
    goResult <- summary(hyperGResult)
  
    # Joep now also checks which other genes are regulated with these GO terms
    
    # But let's just use the summary
    return(goResult)
}



# example:
# hclust_out      <- hclust(as.dist(1-correlation_matrix))
# plot_clust_join_density(hclust_out)
plot_clust_join_density <- function(hclust_out) {
    
    # Let's do this a bit rigourously
    # create data frame of joining data
    df_hclust <- data.frame(x1=hclust_out$merge[,1],x2=hclust_out$merge[,2],h=hclust_out$height)
        # note that for each time two clusters are joined, this is stored
        # as a row in merge with points i,j where i is the index of one point and j the
        # the index of the other. Singletons are given negative numbers, 
        # and conglomerate clusters are given positive numbers.
        # height gives the corresponding distance of points being joints.

    # now look at increase of height    
    df_height  <- data.frame(#nr_points=seq(length(heights),1,-1),
                    nr_points=seq(1,length(df_hclust$h)),
                    heights=df_hclust$h)
    
    ################################################################################
    # 1st method to define a cutoff distance

    # In an ideal situation, there would be clearly distinct clusters in the data
    # In this case, all points first are joined to clusters, and then only at the 
    # end are the clusters joined together. At that point, the clusters are separated
    # much farther than the typical point-to-point distance, or point-to-cluster distance.
    # 
    # Thus, we can look at the distance between the last points that were merged into
    # clusters to define a cutoff distance.
    
    # get distances where single points are joined at 
    singleton_heights <- hclust_out$height[(hclust_out$merge[,1]<0)|(hclust_out$merge[,2]<0)]
    # now take max of that as a distance where "clusters" have formed
    # (or, actually, take the 98% percentile to avoid outliers)
    cutoff_distance1 <-
        sort(singleton_heights,decreasing = F)[round(length(singleton_heights)*.98)]
    cutoff_distance2 <- max(singleton_heights)
    
    # let's look at how connectedness increases
    p1<-ggplot(data=df_hclust)+
        geom_point(aes(x=x1,y=h))+
        geom_point(aes(x=x2,y=h),color='red')+
        geom_hline(yintercept = cutoff_distance1)+
        geom_hline(yintercept = cutoff_distance2)+
        ggtitle(paste0('Positive x: conglomerate joining\nNegative x: singleton joining\nHeights of lines: ',
                        round(cutoff_distance1,4),' and ', round(cutoff_distance2,4) ,'.'))+
        give_better_textsize_plot(15)+xlab('Element index')+ylab('Joining distance')
    
    ################################################################################
    # 2nd method to detect a good cutoff distance
    
    # This method tries to find the 'tipping point', where the curve goes from
    # slowly rising (short distances being added) to steeply rising (large distances
    # being added). The idea would be similar to the first method, which is that
    # the typical distance between points in clusters is smaller than the typical
    # distance between clusters of points. 
    
    # give shorthand names to heights and points
    x=df_height$nr_points
    y=df_height$heights
    
    # Extrapolate data between points
    approx_out <- approx(x=x, y=y, method='linear', n = 100)
    
    # Put the extrapolated line into convenient parameters.
    x_ = approx_out$x
    y_ = approx_out$y
    
    # calculate the derivative line of the extrapolated line
    d_y = (y_[2:length(y_)]-y_[1:(length(y_)-1)]) # just the point-to-point difference (dy)
    x_mid = (x_[2:length(y_)]+x_[1:(length(y_)-1)])/2 # the x-locations
    dy_div_dx = (y_[2:length(y_)]-y_[1:(length(y_)-1)])/(x_[2]-x_[1]) # actual derivative; dy divided by dx
    
    # calculate the 2nd derivative
    dd_y = (d_y[2:length(d_y)]-d_y[1:(length(d_y)-1)])
    x_mm = (x_mid[2:length(x_mid)]+x_mid[1:(length(x_mid)-1)])/2
    
    # determine tipping point, which should be reflected by a bump in the 2nd derivative
    # (in an ideal situation)
    # we define a bump by the cumulative sum reaching a certain percentage (e.g. 1/100) of the total area
    x_cutoff1<-x_mid[which(cumsum(dd_y)>sum(dd_y)/100)[1]]
    y_cutoff1<-y_[1+which(cumsum(dd_y)>sum(dd_y)/100)[1]]
    x_cutoff2<-x_mid[which(cumsum(dd_y)>sum(dd_y)/10)[1]]
    y_cutoff2<-y_[1+which(cumsum(dd_y)>sum(dd_y)/10)[1]]
    
    # now plot all of this in a figure
    p2<-ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))+
            geom_point(data=data.frame(x=approx_out$x, y=approx_out$y), aes(x=x,y=y,color='extrapolated'))+
            geom_point(aes(color='data'))+
            geom_vline(xintercept = x_cutoff1)+geom_hline(yintercept = y_cutoff1)+
            geom_vline(xintercept = x_cutoff2)+geom_hline(yintercept = y_cutoff2)+
            geom_text(aes(x=x_cutoff1,y=max(y),label=round(x_cutoff1,3)),hjust=1,color='red')+
            geom_text(aes(x=x_cutoff2,y=max(y)*.9,label=round(x_cutoff2,3)),hjust=1,color='red')+
            geom_text(aes(x=0,y=y_cutoff1*1.07,label=round(y_cutoff1,3)),hjust=0,color='red')+
            geom_text(aes(x=0,y=y_cutoff2*1.07,label=round(y_cutoff2,3)),hjust=0,color='red')+
            scale_color_manual(values=c('black','gray'))+
            xlab('Number of merges')+ylab('Distance of joining')+
            give_better_textsize_plot(8)
    #print(p2)
    
            #theme(legend.position = c(0.8, 0.2))
    p3<-ggplot(data=data.frame(x=x,y=y), aes(x=x,y=y))+
            geom_point(data=data.frame(x=x_mid, y=d_y), aes(x=x,y=y,color='derivative'))+
            geom_line(data=data.frame(x=x_mid, y=d_y), aes(x=x,y=y,color='derivative'))+
            geom_point(data=data.frame(x=x_mm, y=dd_y), aes(x=x,y=y,color='2nd_derivative'))+
            geom_line(data=data.frame(x=x_mm, y=dd_y), aes(x=x,y=y,color='2nd_derivative'))+
            geom_vline(xintercept = x_cutoff1)+geom_text(aes(x=x_cutoff1,y=max(dd_y),label=round(x_cutoff1,3)),hjust=1)+
            geom_vline(xintercept = x_cutoff2)+geom_text(aes(x=x_cutoff2,y=max(dd_y*.9),label=round(x_cutoff2,3)),hjust=1)+
            scale_color_manual(values=c('red','blue'))+
            #theme(legend.position = c(0.8, 0.2))+
            xlab('Number of merges')+ylab('derivative')+
            give_better_textsize_plot(8)
    #p3
    
    g1<-grid.arrange(p2,p3,nrow=2)

    return(list(p1=p1,
        p2=p2,
        p3=p3,g1=g1,
        distance=c(cutoff_distance1,cutoff_distance2),
        y_cutoff=c(y_cutoff1,y_cutoff2)))
    
}




# Function to run DESeq2 differential expression analysis
# Note that setting up a DESeq2 is rather straightforward, so this
# function is not very complicated
run_DESeq2_MW <- function(raw_expression_data, indices_set1, indices_set2, pvalue_sel=0.01, PARRALEL_YESNO=F) {
    
    # Set up a dataframe which combines the two observations
    #
    # Note this data now might some redundancy, but this setup
    # originates from the fact that DESeq2 was designed for
    # 'non single cell' RNA sequencing, where you have two
    # separate sets of samples.
    df_set1_and_set2 <- raw_expression_data[,c(indices_set1, indices_set2)]
    # now create a matching parameter to annotate the sets
    mycondition <-factor(c(rep(1,length(indices_set1)),rep(2,length(indices_set2))))
    # colnames
    mycolnames  <- colnames(df_set1_and_set2)
    # libType
    my_libType  <- rep("single-end", dim(df_set1_and_set2)[2])
    
    # run on sc@expdata
    des <- data.frame( row.names = mycolnames, 
                       condition = mycondition, 
                       libType = my_libType)
    
    cds <- DESeqDataSetFromMatrix(countData=df_set1_and_set2,
                                  colData=des,
                                  design =~ condition)#,...) 
      
    # start timer
    start_time <- Sys.time() 
    
    # now run DESeq
    cds <- DESeq(cds, fitType='local', parallel=PARRALEL_YESNO)
    res <- results(cds)
    
    # end timer
    end_time <- Sys.time()
    runtime<-end_time-start_time; runtime

    
    # convert to dataframe
    df_res<-as.data.frame(res)
    
    # Add the number of cells each gene is found in
    df_res$num_cells = apply(raw_expression_data>0,1,sum)
    
    # now select & process
    # establish selection indices (ones with a p-value<0.1 and not NA)
    sel_idxs<-(df_res$padj<pvalue_sel)&(!is.na(df_res$padj))
        # According to Michael Love (DESeq author) "The genes with adjusted 
        # p-value of NA have less mean normalized counts than the optimal threshold." 
        # "(..)  it's a sliding scale of "counts high enough to have good power 
        # to detect differential expression".
        # see https://support.bioconductor.org/p/76144/
        # therefor i think it's best not to take those genes into account

    # apply selection
    df_res_sel <- df_res[sel_idxs,]
    # obtain gene names
    gene_names_sel <- rownames(df_res)[sel_idxs]
    # add fold-change columns and gene name column
    df_res_sel<-mutate(df_res_sel,
                       fc=2^df_res_sel$log2FoldChange,
                       fc_inv=2^-df_res_sel$log2FoldChange,
                       gene_name=gene_names_sel)
    # re-establish rownames
    rownames(df_res_sel)<-gene_names_sel
    # order by FC
    df_res_sel<-df_res_sel[order(df_res_sel$fc,decreasing = T),]

    return(list(df_res,df_res_sel))
    
}


# calculate p-values for correlation coefficients
# requires stats library
# See also: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
# df is degrees of freedom, so #observations minus 2 
# (for genes, it is #cells-2: dim(expression_matrix)[2]-2)
calculate_pvalue <- function(R, df) {
    t1<-sqrt(df) * R/sqrt(1 - R^2)
    p_val <- 2*min(pt(t1,df),pt(t1,df,lower.tail = FALSE))
    return(p_val)
}


# Function to calculate correlations with other genes for a particular gene
# note that genes expressed in few cells can create weird correlations
# example of use:
# get_gene_correlations(expression_matrix_sel,'NGF')
get_gene_correlations <- function(expression_matrix, gene_name, show_hist=T) {
 
    # get index of row with gene expression
    gene_idx<-grepl(paste0('^',gene_name,'_'),rownames(expression_matrix))
       
    # throw error if a single gene wasn't identified
    if (!(sum(gene_idx*1)==1)) {
        stop('Gene not found. Give full gene name without chromosome suffix.')
    }
    
    # calculate correlations
    cors_goi <- apply(expression_matrix,1,cor,y=expression_matrix[gene_idx,])
    
    # calculate p-values
    my_freedom <- dim(expression_matrix)[2]-2
    pvals_goi <- sapply(cors_goi, calculate_pvalue, df=my_freedom)
    # show histogram of p-values (to eyeball false discovery/positive rate)
    if (show_hist) {
        print(hist(pvals_goi)) }
    # calculate adjusted pvalues
    pvals_adjusted <- p.adjust(pvals_goi,method='BH')
    
    # return the correlations
    return_df <- data.frame(gene_name=names(cors_goi), correlations=cors_goi,
               pval=pvals_goi, pval_adjusted=pvals_adjusted,
               self=gene_idx)
    return(return_df)
}













# This function in principle is just a wrapper around the
# Monocle function differentialGeneTest; it executes
# this function where it tests for each of the clusters
# whether being inside/outside the cluster has predictive
# value for the expression of each of the genes.
# I.e. it loops and executes differentialGeneTest for each
# of the clusters.
# It additionally marks transcription factors and ligands
# in this analysis.
#
# TODO: check whether an adjusted p-value should be calculated 
#       they use likelihood-ratio test, but i'm not sure about the details of that
differential_gene_expr_Monocle_clusters <- function(config, HSMM_filtered, mycores=11, PVALCUTOFF = 10^-5) {
    
    all_gene_expression <- Matrix::t(Matrix::t(exprs(HSMM_filtered))/sizeFactors(HSMM_filtered))
    
    # differential expression specific to state
    differential_genes_by_cluster_nr<-list()
    differential_genes_by_cluster_ordered<-list()
    differential_genes_by_cluster_ordered_up<-list()
    differential_genes_by_cluster_ordered_down<-list()
    goResult_percluster_up<-list()
    goResult_percluster_down<-list()
    goResult_percluster_seen<-list()
    goResult_percluster_up_100<-list()
    goResult_percluster_down_100<-list()
    for (cluster_idx in 1:max(as.numeric(levels(HSMM_filtered$Cluster)))) {
        
        print(paste0("Analyzing cluster ", cluster_idx))
        
        # Create boolean flag that gives cluster membership for cluster X
        HSMM_filtered[[paste0('cluster',cluster_idx)]] <- HSMM_filtered$Cluster==cluster_idx
        
        # calculate genes that are differentially expression
        differential_genes_by_cluster_nr[[cluster_idx]] <-  
            differentialGeneTest(HSMM_filtered, 
                                 fullModelFormulaStr = paste0("~cluster",cluster_idx),
                                 cores=mycores)
        # now calculate median:median ratio
        median_expression_inside  <- apply(all_gene_expression[,HSMM_filtered[[paste0("cluster",cluster_idx)]]],1,median)
        median_expression_outside <- apply(all_gene_expression[,!HSMM_filtered[[paste0("cluster",cluster_idx)]]],1,median)
        mean_expression_inside  <- apply(all_gene_expression[,HSMM_filtered[[paste0("cluster",cluster_idx)]]],1,mean)
        mean_expression_outside <- apply(all_gene_expression[,!HSMM_filtered[[paste0("cluster",cluster_idx)]]],1,mean)
        expression_ratios<-median_expression_inside/median_expression_outside
        expression_ratios_mean<-mean_expression_inside/mean_expression_outside
        differential_genes_by_cluster_nr[[cluster_idx]]$ratio<-expression_ratios_mean[as.character(differential_genes_by_cluster_nr[[cluster_idx]]$gene_short_name)]    
        # select significant ones for up and down regulation
        differential_genes_by_cluster_ordered[[cluster_idx]]<-differential_genes_by_cluster_nr[[cluster_idx]][order(differential_genes_by_cluster_nr[[cluster_idx]]$pval,decreasing = F),]

        # Now add information about whether the genes are transcription factors
        differential_genes_by_cluster_ordered[[cluster_idx]][['TranscriptionFactor_Flag']]<-F
        differential_genes_by_cluster_ordered[[cluster_idx]][['TranscriptionFactor_Flag']][give_transcription_factors_idxs(rownames(differential_genes_by_cluster_ordered[[cluster_idx]]), human_TFs)] <- T
        # Now add information about whether the genes are ligands
        ligands_out    <- give_idxs_on_gene_list(gene_query = rownames(differential_genes_by_cluster_ordered[[cluster_idx]]), gene_list = ligand_receptor_pairs_df$ligand)
        receptors_assoc <-get_R_associated(ligand_receptor_pairs_df = ligand_receptor_pairs_df, ligand_search_out = ligands_out)
        differential_genes_by_cluster_ordered[[cluster_idx]][['Ligand_Flag']]<-F
        differential_genes_by_cluster_ordered[[cluster_idx]][['Ligand_Flag']][ligands_out$boolean_query] <- T
        differential_genes_by_cluster_ordered[[cluster_idx]][['R_associated']]<-receptors_assoc$R_associated_filled_str
        
        # make selections
        differential_genes_by_cluster_ordered_up[[cluster_idx]]<-differential_genes_by_cluster_ordered[[cluster_idx]][differential_genes_by_cluster_ordered[[cluster_idx]]$pval<PVALCUTOFF&differential_genes_by_cluster_ordered[[cluster_idx]]$ratio>1,]
        differential_genes_by_cluster_ordered_down[[cluster_idx]]<-differential_genes_by_cluster_ordered[[cluster_idx]][differential_genes_by_cluster_ordered[[cluster_idx]]$pval<PVALCUTOFF&differential_genes_by_cluster_ordered[[cluster_idx]]$ratio<1,]
    
        # store the differentially regulated genes
        df_export_up   <- as.data.frame(differential_genes_by_cluster_ordered_up[[cluster_idx]])
        df_export_down <- as.data.frame(differential_genes_by_cluster_ordered_down[[cluster_idx]])
        
        # Export data to Excel sheet
        #Per sheet would be better, but for some reason this is not working.. Bug??
        "if (cluster_idx==1) {
            xlsx::write.xlsx(df_export_up, file=paste0(config$directory_with_data,'out_mw_analyses/diff_genes_per_cluster.xlsx'), sheetName=paste0('cluster',cluster_idx,'up'))
            xlsx::write.xlsx(df_export_down, file=paste0(config$directory_with_data,'out_mw_analyses/diff_genes_per_cluster.xlsx'), sheetName=paste0('cluster',cluster_idx,'down'),append=TRUE)
        } else {
            xlsx::write.xlsx(df_export_up,   file=paste0(config$directory_with_data,'out_mw_analyses/diff_genes_per_cluster.xlsx'), sheetName=paste0('cluster',cluster_idx,'up'),   append=TRUE)
            xlsx::write.xlsx(df_export_down, file=paste0(config$directory_with_data,'out_mw_analyses/diff_genes_per_cluster.xlsx'), sheetName=paste0('cluster',cluster_idx,'down'), append=TRUE)
        }"
        # Export each to separate file
        outputDir_sub = paste0(config$outputDir, '/Monocle_Differential_Expression/')
        if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
        xlsx::write.xlsx(df_export_up,   file=paste0(outputDir_sub, 'Monocle_diff_genes_',config$current_run_name,'_cluster',cluster_idx,'up.xlsx'), sheetName=paste0('cluster',cluster_idx,'up'))
        xlsx::write.xlsx(df_export_down, file=paste0(outputDir_sub, 'Monocle_diff_genes_',config$current_run_name,'_cluster',cluster_idx,'down.xlsx'), sheetName=paste0('cluster',cluster_idx,'down'))
        
        # This is somethin that perhaps should be added later (GO analysis)
        "
        # Now perform GO analysis for those
        goResult_percluster_up[[cluster_idx]] <- give_GO_terms(config = config, gene_names_SCS = rownames(all_gene_expression), 
                                  genes_query = as.character(differential_genes_by_cluster_ordered_up[[cluster_idx]]$gene_short_name),
                                  genes_background = rownames(all_gene_expression))
        # Now perform GO analysis for those
        goResult_percluster_down[[cluster_idx]] <- give_GO_terms(config = config, gene_names_SCS = rownames(all_gene_expression), 
                              genes_query = as.character(differential_genes_by_cluster_ordered_down[[cluster_idx]]$gene_short_name),
                              genes_background = rownames(all_gene_expression))
        
        # Now repeat but take top #100 --- note that p-value cutoff is not applied here..
        # Now perform GO analysis for those
        goResult_percluster_up_100[[cluster_idx]] <- give_GO_terms(config = config, gene_names_SCS = rownames(all_gene_expression), 
                                  genes_query = as.character(differential_genes_by_cluster_ordered_up[[cluster_idx]][1:100,]$gene_short_name),
                                  genes_background = rownames(all_gene_expression))
        # Now perform GO analysis for those
        goResult_percluster_down_100[[cluster_idx]] <- give_GO_terms(config = config, gene_names_SCS = rownames(all_gene_expression), 
                              genes_query = as.character(differential_genes_by_cluster_ordered_down[[cluster_idx]][1:100,]$gene_short_name),
                              genes_background = rownames(all_gene_expression))
        "
        
    }
    
    # and finally export background
    outputDir_sub = paste0(config$outputDir, '/Monocle_Differential_Expression/')
    if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    df_export_allnames   <- as.data.frame(rownames(all_gene_expression))
    write.xlsx(df_export_allnames, file=paste0(outputDir_sub, 'diff_genes_cluster_background_all_genes_detected_',config$current_run_name,'.xlsx'), sheetName=paste0('background_all_genes_detected'),append=F)

    return(list(differential_genes_by_cluster_ordered,differential_genes_by_cluster_ordered_up,differential_genes_by_cluster_ordered_down,df_export_allnames))
    
}

# Note: For RaceID, already a p-value cutoff has been made, so PVAL_CUTOFF is just here to offer
# additional selection
process_and_export_race_clusters_DE <- function(RaceID_object, RaceID_diff_genes_all_cls, PVAL_CUTOFF=0.01) {
    
    # Now let's add a field to identify TFs and ligands, and export the result to an Excel file
    # Loop over the different clusters again
    RaceID_diff_genes_all_cls_selected <- list()
    for (cl_idx in unique(RaceID_object@cpart)) {
        print(paste0('Post-processing cluster ', cl_idx, '..'))
        
        # Obtain data frame with significantly differentially expressed genes, for current cluster
        RaceID_diff_genes_all_cls_selected[[cl_idx]] <- 
            RaceID_diff_genes_all_cls[[cl_idx]][RaceID_diff_genes_all_cls[[cl_idx]]$padj<PVAL_CUTOFF,]
    
        # For export purposes, also add the gene names
        RaceID_diff_genes_all_cls_selected[[cl_idx]]$gene_names =                  
            rownames(RaceID_diff_genes_all_cls_selected[[cl_idx]])
                
        # determine which of those are transcription factors
        RaceID_diff_genes_all_cls_selected[[cl_idx]][['TF_flag']] <- 
            get_idxs_on_gene_list(rownames(RaceID_diff_genes_all_cls_selected[[cl_idx]]), human_TFs)$boolean_query
        
        # determine which of those are ligands
        Ligand_search_result <- 
            get_idxs_on_gene_list(rownames(RaceID_diff_genes_all_cls_selected[[cl_idx]]), 
                ligand_receptor_pairs_df$ligand)
        RaceID_diff_genes_all_cls_selected[[cl_idx]][['Ligand_flag']] <- Ligand_search_result$boolean_query
        # get associated receptors
        associated_R <- get_R_associated(ligand_receptor_pairs_df, Ligand_search_result)[[3]]
        RaceID_diff_genes_all_cls_selected[[cl_idx]][['Associated_Receptors']] <- associated_R
        
        # Now the other way around, determine which of those are receptors
        receptor_search_result <- 
            get_idxs_on_gene_list(rownames(RaceID_diff_genes_all_cls_selected[[cl_idx]]), 
                ligand_receptor_pairs_df$receptor)
        RaceID_diff_genes_all_cls_selected[[cl_idx]][['Receptor_flag']] <- receptor_search_result$boolean_query
        # get associated ligands
        associated_L <- get_X_associated(ligand_receptor_pairs_df, receptor_search_result)[[3]]
        RaceID_diff_genes_all_cls_selected[[cl_idx]][['Associated_Ligands']] <- associated_L
        
        # Now export this excel file
        # store the differentially regulated genes
        df_export_up   <- filter(RaceID_diff_genes_all_cls_selected[[cl_idx]], fc>1)
        df_export_down <- filter(RaceID_diff_genes_all_cls_selected[[cl_idx]], fc<1)
        
        # If empty, add row with NA placeholders
        if (dim(df_export_up)[1] == 0) { df_export_up[1,]=NA }
        if (dim(df_export_down)[1] == 0) { df_export_down[1,]=NA }
        
        # Export each to separate file
        outputDir_sub = paste0(config$outputDir, '/Race_Differential_Expression/')
        if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
        xlsx::write.xlsx(df_export_up,   file=paste0(outputDir_sub, 'Race_diff_genes_',config$current_run_name,'_cluster',cl_idx,'up.xlsx'), sheetName=paste0('cluster',cl_idx,'up'))
        xlsx::write.xlsx(df_export_down, file=paste0(outputDir_sub, 'Race_diff_genes_',config$current_run_name,'_cluster',cl_idx,'down.xlsx'), sheetName=paste0('cluster',cl_idx,'down'))
        
    }
    
    return(RaceID_diff_genes_all_cls_selected)
}

column1_to_rownames <- function(my_data) {
 
    row.names(my_data) <- my_data[,1]
    my_data <- my_data[,2:length(my_data)]
       
}



################################################################################
# Differential gene expression calculated by Monocle for two arbitrary sets of genes



# =====
# To calculate DE for two arbitrary sets in RaceID, use:
# A <- cellnames_set_of_interest
# B <- cellnames_set_of_reference
# diffexpnb_out <- RaceID::diffexpnb(getfdata(RaceID_analysis_out_allconditions$RaceID_object,n=c(A,B)), 
#     A=A, 
#     B=B )
# diffexpnb_out$res holds the results
# =====
# The function below further processes the output data
process_and_export_race_customset_DE <- function(RaceID_object, diffexpnb_out_res, comparison_name, PVAL_CUTOFF=0.01) {
    
    # Now let's add a field to identify TFs and ligands, and export the result to an Excel file

    print(paste0('Post-processing custom set ', comparison_name, '..'))
    
    # Obtain data frame with significantly differentially expressed genes, for current cluster
    diffexpnb_out_selected <- 
        diffexpnb_out_res[diffexpnb_out_res$padj<PVAL_CUTOFF,]

    # For export purposes, also add the gene names
    diffexpnb_out_selected$gene_names =                  
        rownames(diffexpnb_out_selected)
            
    # determine which of those are transcription factors
    diffexpnb_out_selected[['TF_flag']] <- 
        get_idxs_on_gene_list(rownames(diffexpnb_out_selected), human_TFs)$boolean_query
    
    # determine which of those are ligands
    Ligand_search_result <- 
        get_idxs_on_gene_list(rownames(diffexpnb_out_selected), 
            ligand_receptor_pairs_df$ligand)
    diffexpnb_out_selected[['Ligand_flag']] <- Ligand_search_result$boolean_query
    # get associated receptors
    associated_R <- get_R_associated(ligand_receptor_pairs_df, Ligand_search_result)[[3]]
    diffexpnb_out_selected[['Associated_Receptors']] <- associated_R
    
    # Now the other way around, determine which of those are receptors
    receptor_search_result <- 
        get_idxs_on_gene_list(rownames(diffexpnb_out_selected), 
            ligand_receptor_pairs_df$receptor)
    diffexpnb_out_selected[['Receptor_flag']] <- receptor_search_result$boolean_query
    # get associated ligands
    associated_L <- get_X_associated(ligand_receptor_pairs_df, receptor_search_result)[[3]]
    diffexpnb_out_selected[['Associated_Ligands']] <- associated_L
    
    # Now export this excel file
    # store the differentially regulated genes
    df_export_up   <- filter(diffexpnb_out_selected, foldChange>1)
    df_export_down <- filter(diffexpnb_out_selected, foldChange<1)
    
    # Export each to separate file
    outputDir_sub = paste0(config$outputDir, '/Race_Differential_Expression/')
    if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    xlsx::write.xlsx(df_export_up,   file=paste0(outputDir_sub, 'Race_diff_genes_',config$current_run_name,'_customSet_',comparison_name,'_up.xlsx'), sheetName=paste0(comparison_name,'_up'))
    xlsx::write.xlsx(df_export_down, file=paste0(outputDir_sub, 'Race_diff_genes_',config$current_run_name,'_customSet_',comparison_name,'_down.xlsx'), sheetName=paste0(comparison_name,'_down'))
    

    
    return(diffexpnb_out_res)
}


################################################################################
# Differential gene expression calculated by Monocle for two arbitrary sets of genes
differential_gene_expr_Monocle_custom_sets <- function(config, HSMM_filtered, set_interest, set_reference, comparison_name_str, mycores=2, PVALCUTOFF = 10^-5) {
    
    all_gene_expression <- Matrix::t(Matrix::t(exprs(HSMM_filtered))/sizeFactors(HSMM_filtered))
        
    print(paste0('Analyzing differential expression'))
    
    # Create boolean flag that gives cluster membership for cluster X
    HSMM_filtered[[paste0('set_of_interest')]]  <- set_interest
    HSMM_filtered[[paste0('set_of_reference')]] <- set_reference
    # Only keep cells that are in one of the sets
    current_HSMM_filtered = HSMM_filtered[,HSMM_filtered$set_of_interest | HSMM_filtered$set_of_reference]
    current_gene_expression = all_gene_expression[,HSMM_filtered$set_of_interest | HSMM_filtered$set_of_reference,]
    
    # calculate genes that are differentially expression
    differential_genes <-  
        differentialGeneTest(current_HSMM_filtered, 
                             fullModelFormulaStr = paste0("~set_of_interest"),
                             cores=mycores)
    
    # now calculate median:median ratio and mean:mean ratio
    median_expression_inside  <- apply(current_gene_expression[,current_HSMM_filtered$set_of_interest],1,median)
    median_expression_outside <- apply(current_gene_expression[,current_HSMM_filtered$set_of_reference],1,median)
    mean_expression_inside  <- apply(current_gene_expression[,current_HSMM_filtered$set_of_interest],1,mean)
    mean_expression_outside <- apply(current_gene_expression[,current_HSMM_filtered$set_of_reference],1,mean)
    expression_ratios<-median_expression_inside/median_expression_outside
    expression_ratios_mean<-mean_expression_inside/mean_expression_outside
    # save the mean ratio
    differential_genes$ratio<-expression_ratios_mean[as.character(differential_genes$gene_short_name)]    
    
    # select significant ones for up and down regulation
    differential_genes_ordered<-differential_genes[order(differential_genes$pval,decreasing = F),]

    # Now add information about whether the genes are transcription factors
    differential_genes_ordered[['TranscriptionFactor_Flag']]<-F
    differential_genes_ordered[['TranscriptionFactor_Flag']][give_transcription_factors_idxs(rownames(differential_genes_ordered), human_TFs)] <- T
    # Now add information about whether the genes are ligands
    ligands_out    <- give_idxs_on_gene_list(gene_query = rownames(differential_genes_ordered), gene_list = ligand_receptor_pairs_df$ligand)
    receptors_assoc <-get_R_associated(ligand_receptor_pairs_df = ligand_receptor_pairs_df, ligand_search_out = ligands_out)
    differential_genes_ordered[['Ligand_Flag']]<-F
    differential_genes_ordered[['Ligand_Flag']][ligands_out$boolean_query] <- T
    differential_genes_ordered[['R_associated']]<-receptors_assoc$R_associated_filled_str
    
    # Now the other way around, determine which of those are receptors
    receptor_search_result <- 
        get_idxs_on_gene_list(rownames(differential_genes_ordered), 
            ligand_receptor_pairs_df$receptor)
    differential_genes_ordered[['Receptor_flag']] <- receptor_search_result$boolean_query
    # get associated ligands
    associated_L <- get_X_associated(ligand_receptor_pairs_df, receptor_search_result)[[3]]
    differential_genes_ordered[['Associated_Ligands']] <- associated_L

    # make selections
    differential_genes_ordered_up<-differential_genes_ordered[differential_genes_ordered$pval<PVALCUTOFF&differential_genes_ordered$ratio>1,]
    differential_genes_ordered_down<-differential_genes_ordered[differential_genes_ordered$pval<PVALCUTOFF&differential_genes_ordered$ratio<1,]

    # store the differentially regulated genes
    df_export_up   <- as.data.frame(differential_genes_ordered_up)
    df_export_down <- as.data.frame(differential_genes_ordered_down)
    
    # Export each to separate file
    outputDir_sub = paste0(config$outputDir, '/Monocle_Differential_Expression/')
    if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    xlsx::write.xlsx(df_export_up,   file=paste0(outputDir_sub, 'Monocle_diff_genes_',config$current_run_name,'_customset_',comparison_name_str,'_up.xlsx'), sheetName=paste0(comparison_name_str,'_up'))
    xlsx::write.xlsx(df_export_down, file=paste0(outputDir_sub, 'Monocle_diff_genes_',config$current_run_name,'_customset_',comparison_name_str,'_down.xlsx'), sheetName=paste0(comparison_name_str,'_down'))
    
    # This is somethin that perhaps should be added later (GO analysis)
    ## Now perform GO analysis for those
    #goResult_percluster_up <- give_GO_terms(config = config, gene_names_SCS = rownames(current_gene_expression), 
    #                          genes_query = as.character(differential_genes_ordered_up$gene_short_name),
    #                          genes_background = rownames(current_gene_expression))
    ## Now perform GO analysis for those
    #goResult_percluster_down <- give_GO_terms(config = config, gene_names_SCS = rownames(current_gene_expression), 
    #                      genes_query = as.character(differential_genes_ordered_down$gene_short_name),
    #                      genes_background = rownames(current_gene_expression))
    #
    ## Now repeat but take top #100 --- note that p-value cutoff is not applied here..
    ## Now perform GO analysis for those
    #goResult_percluster_up_100 <- give_GO_terms(config = config, gene_names_SCS = rownames(current_gene_expression), 
    #                          genes_query = as.character(differential_genes_ordered_up[1:100,]$gene_short_name),
    #                          genes_background = rownames(current_gene_expression))
    ## Now perform GO analysis for those
    #goResult_percluster_down_100 <- give_GO_terms(config = config, gene_names_SCS = rownames(current_gene_expression), 
    #                      genes_query = as.character(differential_genes_ordered_down[1:100,]$gene_short_name),
    #                      genes_background = rownames(current_gene_expression))
    
    return(list(differential_genes_ordered=differential_genes_ordered,
        differential_genes_ordered_up=differential_genes_ordered_up,
        differential_genes_ordered_down=differential_genes_ordered_down))
    
}




################################################################################
# New volcano stuff

# Volcano plot of NPPA
#GENE = 'NPPA__chr1'
get_volcano_df = function(expression_matrix,my_gene,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1) {
    
    # get the expression of the gene of interest (goi)
    goi_expr = get_expression_gene2(expression_matrix, my_gene)
    if (mean(goi_expr>no_expression_val)<min_cell_expressed) {
      warning(paste0('Your gene of interest is expressed in only ',round(mean(goi_expr>no_expression_val)*100,1), '% cells (below min treshold))'))
    }
    
    # Determine which genes to take along
    seen_in_frac_cells = apply(expression_matrix>no_expression_val, 1, mean)
    
    # calculate the correlations again, but now for cell_counts vs. genes
    start_time <- Sys.time()
    my_corrs_current = lapply(as.data.frame(t(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])), function(X) {cor.test(X,goi_expr)})
    end_time <- Sys.time(); elapsed_time <- end_time - start_time; elapsed_time
    
    # Put result in dataframe
    my_corrs_df_current = data.frame(corr=unlist(lapply(my_corrs_current,function(X) {X$estimate})),
                  pval=unlist(lapply(my_corrs_current,function(X) {X$p.value})))
    
    # Add some more params
    my_corrs_df_current$pval.adj = p.adjust(my_corrs_df_current$pval, method = 'BH') # BH
    if (calc_qvals) {
      my_corrs_df_current$qval = pval2qval(my_corrs_df_current$pval)[[1]]
    }    
    
    rownames(my_corrs_df_current) = rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])
    my_corrs_df_current$gene_name = rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])
    my_corrs_df_current$gene_name_short = strip__chrXX(rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,]))
    
    # add some more parameter to export to excel file
    my_corrs_df_current$nrCellsExpressGene  = apply(expression_matrix[seen_in_frac_cells>=min_cell_expressed,]>no_expression_val,1,sum)
    #my_corrs_df_current$averageExpression   = averageExpression_current
    
    # export correlations to excel file
    #outputDir_sub = paste0(outputDir2, '/Excel/')
    #if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    #xlsx::write.xlsx(my_corrs_df_current,   file=paste0(outputDir_sub, 'correlations_',gsub("__chr(\\d+|[MXY])", '',GENE),'.xlsx'), sheetName=paste0(gsub("__chr(\\d+|[MXY])", '',GENE)))
    
    # look at the p-val distribution to get an impression of the FDR
    #hist(my_corrs_df_current$pval,100)
    
    # Order the output
    my_corrs_df_current = my_corrs_df_current[order(my_corrs_df_current$pval.adj),]
    
    return(my_corrs_df_current)
}

# now plot the volcano
plot_volcano = function(my_corrs_df_current,mytextsize=15,mycex=3,manual_gene_name=NULL,
        NRLABELED=20,mypvaltreshold=0.01,custom_highlight_group=NULL,NRLABELED_hlgroup=NULL,
                  mypointsize=.1, mylinesize=.25) {

    if (!is.null(manual_gene_name)){
        current_gene=manual_gene_name
    } else {
        current_gene = my_corrs_df_current$gene_name_short[my_corrs_df_current$corr==1&!is.na(my_corrs_df_current$corr)]
    }
    
    pval_idx = order(my_corrs_df_current$pval.adj, decreasing = F)
    
    p = ggplot(my_corrs_df_current)+
        geom_point(aes(x=corr,y=-log10(pval.adj)),size=mypointsize)+
        geom_point(data= dplyr::filter(my_corrs_df_current,pval.adj<0.01&corr<0),aes(x=corr,y=-log10(pval.adj)),size=mypointsize, color='blue')+
        geom_point(data= dplyr::filter(my_corrs_df_current,pval.adj<0.01&corr>0),aes(x=corr,y=-log10(pval.adj)),size=mypointsize, color='red')+
        xlab('Correlation')+ylab('-log10(p-value)')+
        #geom_text_repel(data=my_corrs_df_current[my_corrs_df_current$pval<.001,],
        #    mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name),cex=4,color='red',segment.color='gray')+
        give_better_textsize_plot(mytextsize)+
        theme_bw()+
        geom_hline(yintercept = -log10(mypvaltreshold), size=mylinesize)+
        #xlim(c(-.6,.6))+
        ggtitle(paste0('Genes correlated with ', current_gene))
    if (is.null(custom_highlight_group)) {
      # add labels to 1st N genes (determined by p-val)
      p=p+geom_text_repel(data=my_corrs_df_current[pval_idx[1:NRLABELED],],
            mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name_short),color='black',cex=mycex,segment.size	=mylinesize)#,segment.color='gray')+#,cex=.3
    } else {
      # select a custom group of genes
      toshow_df = my_corrs_df_current[my_corrs_df_current$gene_name_short %in% strip__chrXX(custom_highlight_group),]
      # only display 1st N labels if desired
      if (!is.null(NRLABELED_hlgroup)) {toshow_df = toshow_df[1:NRLABELED_hlgroup,]}
      # add labels
      p=p+geom_text_repel(data=toshow_df,
            mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name_short),color='black',cex=mycex)#,segment.color='gray')+#,cex=.3
        
    }
    p
    
}







##########

# used in "shorthand_composite_expression()"
scale_non_zero_2 = function(X, threshold=0) {
    X[X>threshold] = scale(X[X>threshold], center = F, scale = T)
    return(as.vector(X))
}

calculate_composite_expression_values = function(expression_matrix, current_gene_set, show10violins=F, doScale=T) {

    # get full gene names
    current_gene_set_chrXX = find__chrXX(gene_query = current_gene_set, gene_names = rownames(expression_matrix))
    # Get desired subset of expression matrix
    current_composite_matrix = expression_matrix[current_gene_set_chrXX, ]
    
    # Now scale the matrix to make genes "comparable" in expression
    if (doScale) {
      current_composite_matrix_final = t(apply(current_composite_matrix-.1,1,scale_non_zero_2))
      colnames(current_composite_matrix_final) = colnames(current_composite_matrix)
    } else {
      current_composite_matrix_final = current_composite_matrix
    }
      
    # show effect of normalization of genes visually
    if (show10violins) {
        p=ggplot(melt(t(current_composite_matrix_final[1:min(10,length(current_gene_set_chrXX)),])))+
            geom_violin(aes(x=X2,y=value))
        print(p)
    }
    
    # Alternative scaling
    #current_composite_matrix_final = t(scale(t(current_composite_matrix), center = T, scale = T))
    #ggplot(melt(t(current_composite_matrix_final[1:10,])))+
    #    geom_violin(aes(x=Var2,y=value))
    
    # Calculate composite expression (mean of normalized expression)
    composite_expression = apply(current_composite_matrix_final, 2, mean)
    
    return(composite_expression)
}

# searchterms can be list of groups
give_source_forcellnames = function(searchterms, cellnames) {

    # create output parameter
    cell_source = rep(NA, length(cellnames))
    
    idx=0
    for (current_terms in searchterms) {
      
        idx = idx+1
        
        # test if search term map 1:1 to cells (ie cells don't come up multiple times and identities are overwritten)
        if (any(!is.na((cell_source[grepl(paste(current_terms, collapse = '|'), cellnames)])))) {
          stop('Search terms don\'t map uniquely')
        }
        
        # assign identity to cells
        cell_source[grepl(paste(current_terms, collapse = '|'), cellnames)] = idx
    }
    
    return(cell_source)
    
}




