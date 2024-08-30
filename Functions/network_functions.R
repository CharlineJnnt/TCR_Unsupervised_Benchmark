#network functions

#create edges between CDR3s belong to a cluster. 
#cluster_list : list of df with cluster and CDR3/CDR3a_CDR3b columns
make_edges <- function(cluster_list, column){
  all_cluster <- data.frame()
  for(cluster in names(cluster_list)){
    #print(cluster)
    subset <- cluster_list[[cluster]] %>% mutate(score = 0) %>% unique() %>% reshape2::dcast(get(column)~get(column))
    subset[is.na(subset)] <- 1
    subset_melt <- melt(subset) %>% mutate(cluster = cluster) %>% as.data.frame()
    all_cluster <- rbind(all_cluster, subset_melt) 
  }
  return(all_cluster)
}


color_specificity <- function(vertices_list, epi_one, color){
  row_toKeep_list <- list()
  for(i in names(vertices_list)){
    #print(i)
    if(nrow(vertices_list[[i]])>1){#if multiple specificities
      if(vertices_list[[i]]$Epitope==epi_one){
        row_toKeep_df <- vertices_list[[i]] %>% filter(Epitope == epi_one)#keep specific epitope
      }else{row_toKeep_df <- vertices_list[[i]][1,]}#otherwise keep the first row
      
    }else{row_toKeep_df <- vertices_list[[i]][1,]}#otherwise keep the first row
    
    row_toKeep_list[[i]] <- row_toKeep_df 
  }
  
  vertices <- rbindlist(row_toKeep_list) %>% unique() %>% drop_na()
  
  vertices <- vertices  %>% mutate(color_type_epi = ifelse(is.na(Epitope), "grey",#non labeled cdr3
                                                           ifelse(Epitope == epi_one, color, "black")))
  
  vertices <- as.data.frame(vertices)
  return(vertices)
}


make_edges_gliph2 <- function(cluster_list){
  all_cluster <- data.frame()
  for(cluster in names(cluster_list)){
    #print(cluster)
    subset <- cluster_list[[cluster]] %>% mutate(score = 0) %>% unique() 
    subset$CDR3b_CDR3a <- paste(subset$CDR3b, subset$CDR3a, sep = "")
    subset <- subset %>% reshape2::dcast(CDR3b_CDR3a~CDR3b_CDR3a, value.var = "score")
    subset[is.na(subset)] <- 1
    subset_melt <- melt(subset) %>% mutate(cluster = cluster) %>% as.data.frame()
    
    colnames(subset_melt) <- c("CDR3b_CDR3a", "variable", "value", "cluster")
    subset_melt$CDR3b_CDR3a <- paste(subset_melt$CDR3b_CDR3a, subset_melt$cluster, sep="_")#add the cluster index to edge list in order to separate cluster
    subset_melt$variable <- paste(subset_melt$variable, subset_melt$cluster, sep="_")
    
    
    all_cluster <- rbind(all_cluster, subset_melt) 
  }
  return(all_cluster)
}