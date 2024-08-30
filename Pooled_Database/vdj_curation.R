#Curation of VDJ database with alpha_beta pairing
library(dplyr)
library(purrr)
library(jsonlite)

directory = "../Data/Database/"


data_vdj <- data.table::fread(paste0(directory, "vdjdb.txt"), header = T)
data_vdj <- filter(data_vdj, species == "HomoSapiens" )

#single chain
data_vdj_alone <- data_vdj %>% filter(complex.id == 0)

data_vdj_alone_alpha <- data_vdj_alone %>% filter(gene == "TRA") %>% select(c(complex.id, cdr3, v.segm, j.segm, antigen.epitope, antigen.gene, antigen.species, mhc.class, method, vdjdb.score, web.method, web.method.seq, reference.id))
colnames(data_vdj_alone_alpha) <- c("complex.id", "CDR3a", "TRAV", "TRAJ", "Epitope", "Gene", "Species", "Cell_subset", "Identification_method", "vdjdb.score", "web.method", "web.method.seq", "PubMed_ID")
data_vdj_alone_alpha <- data_vdj_alone_alpha %>% mutate(CDR3b=NA, TRBV = NA, TRBJ = NA) %>% relocate(c(complex.id, Epitope,Gene,Species, Cell_subset, PubMed_ID, Identification_method, vdjdb.score, web.method, web.method.seq), .after = TRBJ) %>% select(., -complex.id)

data_vdj_alone_beta <- data_vdj_alone %>% filter(gene == "TRB") %>% select(c(complex.id, cdr3, v.segm, j.segm, antigen.epitope, antigen.gene, antigen.species, mhc.class, reference.id,  method, vdjdb.score, web.method, web.method.seq))
colnames(data_vdj_alone_beta) <- c("complex.id", "CDR3b", "TRBV", "TRBJ", "Epitope", "Gene", "Species", "Cell_subset","PubMed_ID", "Identification_method", "vdjdb.score", "web.method", "web.method.seq")
data_vdj_alone_beta <- data_vdj_alone_beta %>% mutate(CDR3a=NA, TRAV = NA, TRAJ = NA) %>% relocate(c(complex.id, CDR3a,TRAV,TRAJ), .before = CDR3b) %>% select(., -complex.id)

#alpha beta pairing
data_vdj <- data_vdj %>% filter(complex.id != 0)
data_vdj_beta <- data_vdj %>%  filter(gene == "TRB") %>% select(c(complex.id, cdr3, v.segm, j.segm, antigen.epitope, antigen.gene, antigen.species, mhc.class,reference.id, method, vdjdb.score, web.method, web.method.seq))
colnames(data_vdj_beta) <- c("complex.id", "CDR3b", "TRBV", "TRBJ", "Epitope", "Gene", "Species", "Cell_subset","PubMed_ID", "Identification_method", "vdjdb.score", "web.method", "web.method.seq")

data_vdj_alpha <- data_vdj %>%  filter(gene == "TRA") %>% select(c(complex.id, cdr3, v.segm, j.segm, antigen.epitope, antigen.gene, antigen.species, mhc.class,reference.id, method, vdjdb.score, web.method, web.method.seq))
colnames(data_vdj_alpha) <- c("complex.id", "CDR3a", "TRAV", "TRAJ", "Epitope", "Gene", "Species", "Cell_subset", "PubMed_ID","Identification_method", "vdjdb.score", "web.method", "web.method.seq")

data_vdj_alphabeta <- merge(data_vdj_alpha, data_vdj_beta, by = c("complex.id","Epitope", "Gene", "Species", "Cell_subset", "PubMed_ID","Identification_method", "vdjdb.score", "web.method", "web.method.seq"))
data_vdj_alphabeta <- data_vdj_alphabeta %>% relocate(c(complex.id,Epitope,Gene,Species, Cell_subset, PubMed_ID,Identification_method, vdjdb.score, web.method, web.method.seq), .after = TRBJ) %>% select(., -complex.id)

#merge single chain and pairing
data_vdj_final <- rbind(data_vdj_alphabeta, data_vdj_alone_alpha) %>% rbind(., data_vdj_alone_beta) %>% unique()
data_vdj_final <- subset(data_vdj_final, nchar(as.character(data_vdj_final$CDR3a))<=23 & nchar(as.character(data_vdj_final$CDR3a))>=6 | nchar(as.character(data_vdj_final$CDR3b))<=23 & nchar(as.character(data_vdj_final$CDR3b))>=6)
colnames(data_vdj_final) <- c("CDR3_alpha", "V_alpha", "J_alpha", "CDR3_beta", "V_beta", "J_beta", "Epitope", "Antigen", "Antigen_organism", "Cell_subset", "PubMed_ID","Identification_method", "vdjdb.score", "web.method", "web.method.seq")
data_vdj_final$Cell_subset[data_vdj_final$Cell_subset =="MHCI"] <- "CD8"
data_vdj_final$Cell_subset[data_vdj_final$Cell_subset =="MHCII"] <- "CD4"
data_vdj_final$V_alpha <- gsub("\\*0.*","\\1", data_vdj_final$V_alpha)
data_vdj_final$J_alpha <- gsub("\\*0.*","\\1", data_vdj_final$J_alpha)
data_vdj_final$V_beta <- gsub("\\*0.*","\\1", data_vdj_final$V_beta)
data_vdj_final$J_beta <- gsub("\\*0.*","\\1", data_vdj_final$J_beta)

data_vdj_final <- data_vdj_final %>% #split the Identification_method column
  select(-Identification_method) %>%
  bind_cols(map_df(data_vdj_final$Identification_method, fromJSON))

#sequencing method
data_vdj_final <- data_vdj_final %>% mutate(Sequencing_method = ifelse(web.method.seq == sequencing, sequencing, ifelse(sequencing == "", web.method.seq, ifelse(web.method.seq == "amplicon", sequencing, paste0(web.method.seq, '_', sequencing)))))
data_vdj_final <- data_vdj_final %>% select(-c(web.method, web.method.seq, sequencing, singlecell))

#antigen identification
data_vdj_final <- data_vdj_final %>% mutate(Antigen_type = ifelse(grepl(pattern = "antigen-expressing", x = identification)==T & grepl(pattern = "sort", x = identification)==T, "antigenic organism, protein or peptide, Linear peptide",
                                                                  ifelse(grepl(pattern = "antigen-expressing", x = identification)==T, "antigenic organism, protein or peptide",
                                                                         ifelse(grepl(pattern = "antigen-loaded", x = identification)==T & grepl(pattern = "sort", x = identification)==T, "antigenic peptide, Linear peptide",
                                                                                ifelse(grepl(pattern = "antigen-loaded-target", x = identification)==T, "antigenic peptide",
                                                                                       ifelse(grepl(pattern = "sort", x = identification)==T & grepl(pattern = "peptide", x = identification)==T, "antigenic peptide, Linear peptide",
                                                                                              ifelse(grepl(pattern = "peptide", x = identification)==T, "antigenic peptide", 
                                                                                                     ifelse(grepl(pattern = "sort", x = identification)==T, "Linear peptide",
                                                                                                            ifelse(identification == "", NA, "other")))))))))

data_vdj_final$identification[data_vdj_final$identification == "cla"] <- "CLA-pos"
data_vdj_final <- data_vdj_final %>% select(-c(frequency, vdjdb.score, verification))
data_vdj_final <- data_vdj_final %>% relocate(Sequencing_method, .before = identification) 
data_vdj_final <- data_vdj_final %>% mutate(Verified_score = ifelse(!is.na(CDR3_alpha) & !is.na(CDR3_beta), 2,
                                                                    ifelse(!is.na(CDR3_alpha) & is.na(CDR3_beta), 1.1, 
                                                                           ifelse(is.na(CDR3_alpha) & !is.na(CDR3_beta), 1.2, 0))))
colnames(data_vdj_final)[13] <- "Antigen_identification"
data_vdj_final <- data_vdj_final %>% mutate(Database = "VDJdb")
data_vdj_final$PubMed_ID <- gsub("PMID:*", "",data_vdj_final$PubMed_ID)
data_vdj_final <- unique(data_vdj_final)


data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://github.com/antigenomics/vdjdb-db/issues/291"] <- "32668259"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://github.com/antigenomics/vdjdb-db/issues/311"] <- "28636589, 32341563"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://github.com/antigenomics/vdjdb-db/issues/313"] <- "33101265"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://github.com/antigenomics/vdjdb-db/issues/315"] <- "33326767"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://github.com/antigenomics/vdjdb-db/issues/326"] <- "35383307"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID=="https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#"] <- "https://pages.10xgenomics.com/rs/446-PBO-704/images/10x_AN047_IP_A_New_Way_of_Exploring_Immunity_Digital.pdf"
data_vdj_final$PubMed_ID[data_vdj_final$PubMed_ID==""] <- NA


#####################################################################################################################################################
rm(data_vdj, data_vdj_alone,  data_vdj_alone_alpha,  data_vdj_alone_beta, data_vdj_alpha, data_vdj_beta, data_vdj_alphabeta)#remove intermediate df #
#####################################################################################################################################################
