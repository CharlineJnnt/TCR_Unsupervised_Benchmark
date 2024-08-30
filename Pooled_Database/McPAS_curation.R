#Curation of Mc-PAS database with alpha_beta pairing
library(dplyr)

directory = "../Data/Database/"


data_mcpas <- read.csv(paste0(directory, "McPAS-TCR.csv"), sep =",", encoding= "UTF-8")
data_mcpas <- data_mcpas%>% filter(Species == "Human") %>% as.data.frame()

line_test <- filter(data_mcpas, CDR3.alpha.aa == "TRAJ37")#incorrect row
new_line <- line_test
new_line[1,] <- c("CALDGPSNTGKLIF", "CATSESSGQTYEQYF", "Human", "Cancer", "Tumor associated antigen (TAA)", "D018290", NA, 1, "Yes", "Yes", NA, NA, "FLCMKALLL", NA, "HLA-A2:01",
                  "PBMC", "CD8", NA, NA, "TRAV16", "TRAJ37", "TRBV15", NA, "TRBJ2-7", NA, NA, NA, 33562731, NA)
data_mcpas <- data_mcpas %>% filter(CDR3.alpha.aa != "TRAJ37")#add corrected line
data_mcpas <- rbind(data_mcpas, new_line) %>% as.data.frame()

#keep NGS for TCR identification and Antigen.identification.method of antigen identification
data_mcpas <- data_mcpas %>% select(c(CDR3.alpha.aa, TRAV, TRAJ, CDR3.beta.aa, TRBV, TRBJ, Category, Antigen.protein, Epitope.peptide, Pathology, MHC, Tissue, T.Cell.Type, Remarks, NGS, Single.cell, Antigen.identification.method, PubMed.ID)) 
data_mcpas$CDR3.alpha.aa <- ifelse(data_mcpas$CDR3.alpha.aa == "", NA, data_mcpas$CDR3.alpha.aa)

#if Remarks != NA, we modify CDR3
data_mcpas <- data_mcpas %>% mutate(CDR3_alpha = ifelse(Remarks == "No final F" & CDR3.alpha.aa != "", paste0(CDR3.alpha.aa, "F"), 
                                                        ifelse(Remarks == "No first Cysteine" & CDR3.alpha.aa != "", paste0("C", CDR3.alpha.aa), 
                                                               ifelse(Remarks == "Brackets represent optional amino acid" & CDR3.alpha.aa != "", NA, CDR3.alpha.aa))), 
                                    CDR3_beta = ifelse(Remarks == "No final F" & CDR3.beta.aa != "", paste0(CDR3.beta.aa, "F"), 
                                                       ifelse(Remarks == "No first Cysteine" & CDR3.beta.aa != "", paste0("C", CDR3.beta.aa), 
                                                              ifelse(Remarks == "Brackets represent optional amino acid" & CDR3.beta.aa != "", NA ,CDR3.beta.aa))))
#if Remarks == NA, we keep original CDR3
data_mcpas$CDR3_alpha <- ifelse(is.na(data_mcpas$CDR3_alpha) & is.na(data_mcpas$Remarks), data_mcpas$CDR3.alpha.aa, data_mcpas$CDR3_alpha)
data_mcpas$CDR3_beta <- ifelse(is.na(data_mcpas$CDR3_beta) & is.na(data_mcpas$Remarks), data_mcpas$CDR3.beta.aa, data_mcpas$CDR3_beta)

#remove column with NA in both CDR3 alpha & beta
data_mcpas <- data_mcpas %>% filter(!(is.na(CDR3_alpha)) | !(is.na(CDR3_beta)))
data_mcpas <- data_mcpas[grepl(pattern = "\\*", data_mcpas$CDR3_alpha)==FALSE,]#remove sequence with "*"
data_mcpas <- data_mcpas[grepl(pattern = "\\*", data_mcpas$CDR3_beta)==FALSE,]#remove sequence with "*"

data_mcpas <- data_mcpas[grepl(pattern = "\xa0", data_mcpas$CDR3_alpha)==FALSE,]#remove non ASCII characters
data_mcpas$CDR3_alpha <- stringr::str_conv(data_mcpas$CDR3_alpha, "UTF-8")
data_mcpas <- subset(data_mcpas, nchar(as.character(data_mcpas$CDR3_alpha))<=23 & nchar(as.character(data_mcpas$CDR3_alpha))>=6 | nchar(as.character(data_mcpas$CDR3_beta))<=23 & nchar(as.character(data_mcpas$CDR3_beta))>=6)

data_mcpas <- data_mcpas %>% select(CDR3_alpha, TRAV, TRAJ, CDR3_beta, TRBV, TRBJ, Epitope.peptide, Antigen.protein, Pathology, T.Cell.Type, NGS, Single.cell, Antigen.identification.method, PubMed.ID)
colnames(data_mcpas) <- c("CDR3_alpha", "V_alpha", "J_alpha", "CDR3_beta", "V_beta", "J_beta", "Epitope", "Antigen", "Antigen_organism", "Cell_subset", "Next_generation_sequencing", "Single_Cell", "Antigen_identification_score", "PubMed_ID")
data_mcpas <- unique(data_mcpas)

data_mcpas <- data_mcpas %>% mutate(Ag_identification = ifelse(Antigen_identification_score == 1.0 , "Peptide-MHC tetramers", 
                                                               ifelse(Antigen_identification_score == 2.0 , "In vitro stimulation", 
                                                                      ifelse(Antigen_identification_score == 2.1 , "In vitro stimulation with a peptide",
                                                                             ifelse(Antigen_identification_score == 2.2 , "In vitro stimulation with a protein", 
                                                                                    ifelse(Antigen_identification_score == 2.3 , "In vitro stimulation with a pathogen",
                                                                                           ifelse(Antigen_identification_score == 2.4 , "In vitro stimulation with tumor cells",
                                                                                                  ifelse(Antigen_identification_score == 2.5 , "Other types of in vitro stimulation",
                                                                                                         ifelse(Antigen_identification_score == 3 , "Isolated from a specific tissue under a specific pathology/condition", 
                                                                                                                ifelse(Antigen_identification_score == 4, "cell line",NA))))))))))


data_mcpas <- data_mcpas %>% mutate(Sequencing_method = ifelse(Next_generation_sequencing == "Yes" & Single_Cell == "Yes", "NGS-SingleCell", 
                                                               ifelse(Next_generation_sequencing == "Yes" & Single_Cell == "No", "NGS-noSingleCell", 
                                                                      ifelse(Next_generation_sequencing == "No" & Single_Cell == "No", "noNGS-noSingleCell", 
                                                                             ifelse(Next_generation_sequencing == "No" & Single_Cell == "Yes", "noNGS-SingleCell", NA)))))
data_mcpas <- data_mcpas %>% select(-c(Next_generation_sequencing, Single_Cell))


data_mcpas <- data_mcpas %>% mutate(Antigen_identification = ifelse(grepl(x= Ag_identification, pattern = "vitro stimulation"), "In vitro stimulation",
                                                                    ifelse(grepl(x= Ag_identification, pattern = "Peptide-MHC tetramers"), "tetramer-sort",
                                                                           ifelse(Ag_identification == "Isolated from a specific tissue under a specific pathology/condition", "Isolated from a specific tissue under a specific pathology/condition", NA))))

data_mcpas <- data_mcpas %>% mutate(Antigen_type = ifelse(grepl(x= Ag_identification, pattern = "peptide"), "antigenic peptide",
                                                          ifelse(grepl(x= Ag_identification, pattern = "protein"), "protein",
                                                                 ifelse(grepl(x= Ag_identification, pattern = "pathogen"), "pathogen", 
                                                                        ifelse(grepl(x= Ag_identification, pattern = "tumor cells"), "tumor cells", 
                                                                               ifelse(grepl(x= Ag_identification, pattern = "Other"), "other", 
                                                                                      ifelse(Ag_identification == "Peptide-MHC tetramers", "Linear peptide",
                                                                                             ifelse(Ag_identification == "cell line", "cell line",NA))))))))

data_mcpas <- data_mcpas %>% select(-c(Ag_identification, Antigen_identification_score))
data_mcpas <- data_mcpas %>% mutate(Verified_alpha = ifelse(is.na(CDR3_alpha) | V_alpha == "mTRDV2-2" , "No", "alpha"))
data_mcpas$Verified_alpha[is.na(data_mcpas$Verified_alpha)] <- "alpha"

data_mcpas <- data_mcpas %>% mutate(Verified_beta = ifelse(is.na(CDR3_beta) , "No", "beta"))
data_mcpas$Verified_alphabeta <- paste(data_mcpas$Verified_alpha, data_mcpas$Verified_beta)

data_mcpas  <- data_mcpas %>% mutate(Verified_score = ifelse(Verified_alphabeta == "alpha beta", 2,
                                                             ifelse(Verified_alphabeta == "alpha No", 1.1, 
                                                                    ifelse(Verified_alphabeta == "No beta", 1.2, 0))))

data_mcpas <- data_mcpas %>% select(-c(Verified_alpha, Verified_beta, Verified_alphabeta))
data_mcpas <- data_mcpas %>% mutate(Database = "Mc-PAS")

data_mcpas$V_alpha <- gsub("\\:.*","", data_mcpas$V_alpha)
data_mcpas$V_beta <- gsub("\\:.*","", data_mcpas$V_beta)
data_mcpas$J_alpha <- gsub("\\:.*","", data_mcpas$J_alpha)
data_mcpas$J_beta <- gsub("\\:.*","", data_mcpas$J_beta)

data_mcpas$V_alpha <- gsub("\\-0","\\-", data_mcpas$V_alpha)
data_mcpas$J_alpha <- gsub("\\-0","\\-", data_mcpas$J_alpha)
data_mcpas$V_beta <- gsub("\\-0","\\-", data_mcpas$V_beta)
data_mcpas$J_beta <- gsub("\\-0","\\-", data_mcpas$J_beta)

data_mcpas$J_beta <- ifelse(data_mcpas$J_beta == "Donor 13", NA, data_mcpas$J_beta)

data_mcpas <- unique(data_mcpas)

########################
rm(line_test, new_line)#
########################