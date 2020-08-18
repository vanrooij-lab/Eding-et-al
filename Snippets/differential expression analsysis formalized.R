# functions for permutation-based statistics to deterimene
# cell population size dymanics across different conditions
# strategy of Ferhabi et al (2019)* was adapted for this script
# this sctipt requires a data.frame with teh first column being cell names which include the condition they are derived from
# the second column should depict the cluster number of each cell.

# SCript by: Bas Molenaar 15-07-2019

# * Nona Farbehi, et al. (2019). Single-cell expression profiling reveals dynamic flux of cardiac stromal, vascular 
#   and immune cells in health and injury. E-life,8


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
##### function to get contigency table of raw number of cells                                                     #####
##### per condition in all clusters                                                                               #####  
##### all_clus = data frame showing the cluster number for each cell that                                         #####
##### was included in the RaceID clustering (after filtering) first column should be cell names,                  #####
##### second should column should be cluster number                                                               #####
##### ... = depict here which condition the cells from a group originates (e.g. here "sham_1d","Sham_14d","IR_1d) #####
##### output is contigency table with every row being different clusters                                          #####
##### and every column being different conditions                                                                 #####
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

create_contigence_table <- function(all_clus,...){
  for(condition in list(...)){
    if (any(grepl(condition,all_clus[,1]))==FALSE){
      stop("Warning, one or more conditions are not in the dataset")
    }
  }
  output <- data.frame()
  clusters_nr <- NULL
  for(i in 1:max(all_clus[,2])){
    clusters_nr <- c(clusters_nr,i)
    num_cell <- NULL
    for(y in seq(1,length(list(...)))){
      cells <- nrow(subset(all_clus,all_clus[,2]==i & grepl(unlist(list(...)[y]),all_clus[,1])))
      num_cell <- c(num_cell,cells) 
    }
    output<- rbind(output,num_cell)
  }
  output <- cbind(clusters_nr,output)
  colnames_data1 <- "cluster number"
  for(colna in list(...)){
    name1 <- paste("nr of cells in cluster from condition",colna)
    colnames_data1<- c(colnames_data1,name1)
  }
  colnames(output)<- colnames_data1
  rownames(output) <- output[,"cluster number"]
  output[,"cluster number"] <- NULL
  return(output)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### function to create the condition and cluster information list for each cell.                          ###
### Creates list with 2 vectors from a contigency table containing nuber of cells from each               ###
### condition per cluster ,with colnamen being different conditions and rownames being different cluster. ###
### Numbers of each element in each vector represents a seperate cell, with information about             ###
### the condition from that cell in the first vector (condition vector) and information about the cluster ###
### in the second row (cluster vector). Used for downstream bootstrapping permutation procedure           ###
### input =  contigency table with every row being different clusters                                     ###
### and every column being different conditions                                                           ###                                                                        ###
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


makevariablesforcells<-function(input){
  cond_names<-colnames(input)
  clus_names<-rownames(input)
  # create condition vector with each condition repeated x times, with x being the total number of cells
  # in a specific condition. Follows the order of the rownames of the contigency table 
  cond<-rep(cond_names, apply(input,2,sum)) 
  # creates cluster vector for each each condition, per condition each cluster is repeated x times, with x being the total
  # number of cells in a a cluster of the specific condition. 
  clus_list<-apply(input,2,function(x){
    rep(clus_names, x )
  })
  # catonate the cluster vectors of each condtion, follows the order of the rownames of the contigency table 
  clus <- NULL
  for(i in seq(1:length(clus_list))){
    clus <- c(clus,clus_list[[i]])
  }
  #table(cond,clus);
  return(list(condition=cond, cluster=clus));
}

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
### function to perform permutation based null distribution for each condition across all clusters           ###
### Creates multi-way array with 3 dimensions. 1 dimension (rownames) depict each cluster                    ### 
### 1 dimension (colnames) depict each condition clusters. last dimension depict number of iterations for    ###
### purmutation based null distribution.                                                                     ###
### IMPORTANT, MAKE SURE TO HAVE THE FUNCTION makevariablesforcells IN THE GENERAL ENVIROMENT                ###
### input =  contigency table with every row being different clusters                                        ###
### and every column being different conditions                                                              ###                                             ###
### n = iteration number for perfoming permutations in background model                                      ###
### p = number of cells for which clusters will be permuted by SRSWOR from the observed cluster destribution ###
### (equation --> number_cells_SRSWOR = total cells*p)                                                       ###
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################


generatebackground_model <-function(input, n=10^5, p=0.1){
  # making sleketal frame for background distribution, with dimnames same as the input contigency table
  background_mod <- array(NA,dim=c(nrow(input),ncol(input),n))
  dimnames(background_mod)[[3]] <- 1:n
  # grenerating lists showing for each cell the condition/cluster combination
  observed_sample <- makevariablesforcells(input)
  # permuting the clusters in a random group of cells by simple random sampling without replacement from the observed
  # cluster distribution
  for(i in 1:n)
  {
    if(i %% 5000 == 0){
      print(paste(i, "of ", n, "permutations"))
    }
    permuted_sample <- observed_sample
    sampled_cells_index <-sample(1:length(permuted_sample$cluster),round(length(permuted_sample$cluster)*p))
    permuted_sample$cluster[sampled_cells_index]<-sample(observed_sample$cluster,length(sampled_cells_index),replace=F)
    permuted_sample_cont_table <- table(permuted_sample$cluster,permuted_sample$condition)
    background_mod[,,i]<-permuted_sample_cont_table
  }
  # manipulation the colnames (the conditions) of the background model in teh array in the correct orde
  dimnames(background_mod)[[1]] <- rownames(permuted_sample_cont_table)
  dimnames(background_mod)[[2]] <- colnames(permuted_sample_cont_table)
  return(background_mod)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
### perform signficance test between all conditions                                                       ###
### returns a 3 dimensional array, with the x and y dimension being all the conditions in from the        ###
### contigency table, and the z dimension being the different clusters                                    ###
### input =  contigency table with every row being different clusters                                     ###
### and every column being different conditions                                                           ###
### background_mod = the background model array from the generatebackground_model function                ###
### output is an array with 3 dimensions, first and second dimension from the P-value matrix in a cluster ###
### between all possible conditions. Third dimension represent the different clusters                     ###
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


DPA_test <- function(input,background_mod){
  p_value_array <- array(-1,dim=c(length(colnames(input)),length(colnames(input)),length(rownames(input))))
  dimnames(p_value_array)[1] <- list(colnames(input))
  dimnames(p_value_array)[2] <- list(colnames(input))
  dimnames(p_value_array)[3] <- list(rownames(input))
  counter <- 0
  for(i in colnames(input)){
    for(y in colnames(input)){
      for(z in rownames(input)){
        observed_difference <- input[z,y]/sum(input[,y])-
          input[z,i]/sum(input[,i])
        background_differences <- background_mod[z,y,]/apply(background_mod[,y,],2,sum)-
          background_mod[z,i,]/apply(background_mod[,i,],2,sum)
        increase_observed <- length(observed_difference[observed_difference>background_differences])/length(background_differences)
        decrease_observed <- length(observed_difference[observed_difference<background_differences])/length(background_differences)
        if(increase_observed==0&decrease_observed==0){
          p_value_array[i,y,z] <- 1
        } else{
          p_value_array[i,y,z] <- min(increase_observed,decrease_observed)
        }
      }
    }
    counter=counter+1
   perc_done <- counter/length(colnames(input))*100
   print(paste("currently",perc_done,"% done"))
  }
  return(p_value_array)
}

