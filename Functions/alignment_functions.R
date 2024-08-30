#alignment functions

#create alignment list
list_alignement_fun <- function(method, output, chain){
  if(method == "Hamming_Distance" | method == "Levenshtein_Distance" | method == "TCRMatch" | method == "iSMART" | method == "GIANA"){
    if(method == "GIANA"){output <- output %>% group_by(clusters) %>% mutate(nb_chain = n_distinct(chain)) %>% filter(nb_chain == 1) %>% select(-nb_chain)} #remove clusters with alpha and beta chain
    output <- merge(output, liste_sequence_TRA_TRB_chain, by= "CDR3")
    output_beta <- output %>% filter(chain =="beta")
    output_alpha <- output %>% filter(chain =="alpha")
    if(chain == "a"){list_alignement <- split(output_alpha[c("clusters", "CDR3")], output_alpha$clusters)  
    }else if(chain =="b"){list_alignement <- split(output_beta[c("clusters", "CDR3")], output_beta$clusters)}
  }
  else if(method == "DeepTCR" | method == "GLIPH2" | method == "TCRdist" | method == "clusTCR"){
    if(chain == "a"){
      output <- select(output, clusters, CDR3a) %>% dplyr::rename(., "CDR3" = "CDR3a")
      list_alignement <- split(output, output$clusters)   
    }else if (chain =="b"){
      output <- select(output, clusters, CDR3b) %>% dplyr::rename(., "CDR3" = "CDR3b")
      list_alignement <- split(output, output$clusters)}
  }
  return(list_alignement)
}
