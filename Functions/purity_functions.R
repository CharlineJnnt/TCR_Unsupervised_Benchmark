# Purity functions

purity_function <- function(output, th, colonne_pair){
  output <- output %>% filter(size >th)#filter on cluster size
  if(missing(colonne_pair)){
    output_spe <- merge(output, dataframe_alpha_beta[,c("CDR3", "Epitope", "chain")], by = "CDR3") %>% unique()
  }else{output_spe <- merge(output[,c(colonne_pair, "clusters", "size")], data_human[,c("CDR3b_CDR3a", "Epitope")], by.x = colonne_pair, by.y = "CDR3b_CDR3a") %>% unique()}
  output_spe <- output_spe %>% group_by(clusters, Epitope) %>% mutate(nb_epi_cluster = n())
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_majoritaire = max(nb_epi_cluster))
  output_spe_sum <- output_spe %>% select(c(clusters, size, nb_majoritaire)) %>% unique()
  
  return(sum(output_spe_sum$nb_majoritaire)/sum(output_spe_sum$size))
}

all_purity_onechain <- function(output, th){
  output <- output %>% filter(size >th)#big clusters only
  output_spe <- merge(output, dataframe_alpha_beta[,c("CDR3", "Epitope", "chain")], by = "CDR3") %>% unique()
  output_spe <- output_spe %>% group_by(clusters, Epitope) %>% mutate(nb_epi_cluster = n())
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_majoritaire = max(nb_epi_cluster))
  output_spe_sum <- output_spe %>% select(c(clusters, size, nb_majoritaire)) %>% unique()
  output_spe_sum <- output_spe_sum %>% mutate(purity = nb_majoritaire/size)
  
  return(output_spe_sum)
}

all_purity_pairing <- function(output, colonne_pair){
  output <- output %>% filter(size >th)
  output_spe <- merge(output[,c(colonne_pair, "clusters", "size")], data_human[,c("CDR3b_CDR3a", "Epitope")], by.x = colonne_pair, by.y = "CDR3b_CDR3a") %>% unique()
  output_spe <- output_spe %>% group_by(clusters, Epitope) %>% mutate(nb_epi_cluster = n())
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_majoritaire = max(nb_epi_cluster))
  output_spe_sum <- output_spe %>% select(c(clusters, size, nb_majoritaire)) %>% unique()
  output_spe_sum <- output_spe_sum %>% mutate(purity = nb_majoritaire/size)
  
  return(output_spe_sum)
}

purity_chain_function <- function(output, colonne_pair){
  output_spe <- merge(output, dataframe_alpha_beta[,c("CDR3", "chain")], by = "CDR3") %>% unique()
  output_spe <- output_spe %>% group_by(clusters, chain) %>% mutate(nb_chain_cluster = n())
  output_spe <- output_spe %>% group_by(clusters) %>% mutate(nb_majoritaire = max(nb_chain_cluster))
  output_spe_sum <- output_spe %>% select(c(clusters, size, nb_majoritaire)) %>% unique()
  
  return(sum(output_spe_sum$nb_majoritaire)/sum(output_spe_sum$size))
}