#Import and process output results from clustering methods
library(stringdist)
library(igraph)
library(readr)

directory <- "../Data/Output_methods/"

## Hamming Distance ##
HD_matrix <- as.matrix(stringdistmatrix(liste_sequence_TRA_TRB$CDR3, method = "hamming",useNames = "strings", useBytes=TRUE))
edge_list_HD <- melt(HD_matrix)
colnames(edge_list_HD) <- c("input_sequence", "match_sequence", "score")
edge_list_HD <- filter(edge_list_HD, score == 1) %>% unique()

graph_HD <- graph_from_data_frame(edge_list_HD, directed = F)
component <- components(graph_HD)
output_HD <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
  group_by(clusters) %>%dplyr:: mutate(size = n())

output_HD_alpha <- merge(output_HD, data.frame(CDR3= liste_sequence_TRA$CDR3), by = "CDR3")
output_HD_beta <- merge(output_HD, data.frame(CDR3= liste_sequence_TRB$CDR3), by = "CDR3")

## Levenshtein Distance ##
LD_matrix <- as.matrix(stringdistmatrix(liste_sequence_TRA_TRB$CDR3, method = "lv",useNames = "strings", useBytes=TRUE))
edge_list_LD <- melt(LD_matrix)
colnames(edge_list_LD) <- c("input_sequence", "match_sequence", "score")
edge_list_LD <- filter(edge_list_LD, score == 1) %>% unique()

graph_LD <- graph_from_data_frame(edge_list_LD, directed = F)
component <- components(graph_LD)
output_LD <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
  group_by(clusters) %>% mutate(size = n())

output_LD_alpha <- merge(output_LD, data.frame(CDR3= liste_sequence_TRA$CDR3), by = "CDR3")
output_LD_beta <- merge(output_LD, data.frame(CDR3= liste_sequence_TRB$CDR3), by = "CDR3")

## TCRMatch ##
edge_list_tcrm <- data.table::fread(paste0(directory,"TCRMatch_output.txt"), header = F)
colnames(edge_list_tcrm) <- c("input_sequence", "match_sequence", "score")
edge_list_tcrm <- edge_list_tcrm %>% filter(score != 1)

graph_tcrm <- graph_from_data_frame(edge_list_tcrm, directed = F)
component <- components(graph_tcrm)
output_tcrm <- data.frame(clusters = component$membership) %>% rownames_as_first_column(., "CDR3") %>% 
  group_by(clusters) %>% mutate(size = n())

output_tcrm_alpha <- merge(output_tcrm, data.frame(CDR3= liste_sequence_TRA$CDR3), by = "CDR3")
output_tcrm_beta <- merge(output_tcrm, data.frame(CDR3= liste_sequence_TRB$CDR3), by = "CDR3")

## iSMART ##
output_ismart <- read.table(paste0(directory, "iSMART_output.txt"), header=F,sep='\t',stringsAsFactors = F)
colnames(output_ismart) <- c('CDR3', 'clusters')
output_ismart <- output_ismart %>% group_by(clusters) %>% mutate(size = n())

output_ismart_alpha <- merge(output_ismart, data.frame(CDR3= liste_sequence_TRA$CDR3), by = "CDR3")
output_ismart_beta <- merge(output_ismart, data.frame(CDR3= liste_sequence_TRB$CDR3), by = "CDR3")

## GIANA ##
output_giana <- read.table(paste0(directory, "GIANA_output.txt"), header=F,sep='\t',stringsAsFactors = F)
colnames(output_giana) <- c('CDR3', 'clusters', 'TRBV')
output_giana <- output_giana %>% group_by(clusters) %>% mutate(size = n()) %>% select(-TRBV)

output_giana_alpha <- merge(output_giana, data.frame(CDR3= liste_sequence_TRA$CDR3), by = "CDR3")
output_giana_beta <- merge(output_giana, data.frame(CDR3= liste_sequence_TRB$CDR3), by = "CDR3")

## clusTCR ##
output_clustcr <- read_delim(file=paste0(directory,"clusTCR_output.csv"), delim = ",")
colnames(output_clustcr) <- c("CDR3b_CDR3a","clusters")
output_clustcr_summary <- read_delim(file=paste0(directory, "clusTCR_output_summary.csv"), delim = ",")
colnames(output_clustcr_summary)[1] <- "clusters"
edge_list_clustcr <- read_delim(file=paste0(directory, "clusTCR_output_edgelist.txt"), delim = "\t", col_names = c("From","To"))
output_clustcr_merge <- merge(output_clustcr, output_clustcr_summary[, c("clusters", "size")], by="clusters")

pairing <- data.frame(CDR3b = data_human$CDR3_beta, CDR3a = data_human$CDR3_alpha) %>% unique()
pairing$CDR3b_CDR3a <- paste(pairing$CDR3b, pairing$CDR3a, sep = "")
output_clustcr_pairing <- merge(output_clustcr_merge, pairing, by.x = "CDR3b_CDR3a", "CDR3b_CDR3a")
colnames(output_clustcr_pairing)[1] <- "CDR3b_CDR3a"

## GLIPH2 ##
output_gliph <- read_delim(file=paste0(directory, "GLIPH2_output_cluster.csv"), delim = ",")
output_gliph <- output_gliph[1:18]
output_gliph <- na.omit(output_gliph)
output_gliph <- filter(output_gliph, pattern != "single") 
output_gliph <- output_gliph %>% group_by(index) %>% mutate(size = n())
output_gliph <- dplyr::rename(output_gliph, "CDR3b" = "TcRb",
                              "CDR3a" = "TcRa", 
                              "clusters" = "index")
output_gliph$CDR3b_CDR3a <- paste(output_gliph$CDR3b, output_gliph$CDR3a, sep = "")
output_gliph <- output_gliph %>% separate(type, c('type_cluster', 'motif'), sep = "-") %>% select(-motif)

## DeepTCR ##
output_deeptcr <- read_delim(file=paste0(directory,"DeepTCR_output.csv"), delim = ",")
output_deeptcr <- dplyr::rename(output_deeptcr,  "CDR3b" = "Beta_Sequences",
                                "CDR3a" ="Alpha_Sequences", 
                                "clusters" = "cluster")
output_deeptcr$CDR3b_CDR3a <- paste(output_deeptcr$CDR3b, output_deeptcr$CDR3a, sep = "")


## TCRdist3 ##
output_tcrdist <- read_delim(file=paste0(directory,"TCRdist3_hclust_output.csv"), delim = ";")
output_tcrdist <- output_tcrdist %>% group_by(clusters) %>% dplyr::mutate(size= n()) %>% filter(size != 1)
output_tcrdist$CDR3b_CDR3a <- paste(output_tcrdist$CDR3b, output_tcrdist$CDR3a, sep = "")


