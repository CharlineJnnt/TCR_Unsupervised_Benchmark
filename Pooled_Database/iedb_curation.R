#Curation of IEDB database with alpha_beta pairing

#I. creation of database pooled
#1) IEDB MHCI & MHCII
#Starting from the IEDB homepage (iedb.org), the following filters were applied:
#- Epitope: " Any Epitopes";
#- Assay: " Positive Assays Only"  and " T Cell Assays" ;
#- Antigen: no filters for Organism or Antigen Name;
#- MHC Restriction: "Class I or II";
#- Host: " Human"; 
#- Disease: " Any Disease" ;
#- TCR type: TCR alpha/beta

library(dplyr)
library(readr)
library(data.table)
library(tidyr)
library(purrr)
library(jsonlite)
library(bio3d)

directory = "../Data/Database/"

#data CD4/CD8 file
print("IEDB database loading")

data_iedb_CD8 <- read_delim(paste0(directory, "IEDB_tcell_receptor_table_MHCI.csv"), delim = ",")
data_iedb_CD8 <- data_iedb_CD8 %>% mutate(Cell_subset = "CD8") %>% unique()

data_iedb_CD4 <- read_delim(paste0(directory,  "IEDB_tcell_receptor_table_MHCII.csv"), delim = ",")
data_iedb_CD4 <- data_iedb_CD4 %>% mutate(Cell_subset = "CD4") %>% unique()

data_iedb_CD48 <- rbind(data_iedb_CD4, data_iedb_CD8)
data_iedb_CD48 <- data_iedb_CD48 %>% select("Epitope ID", "Reference ID", "Chain 1 Type", "Chain 1 Full Sequence",
                                  "Chain 1 CDR3 Start Curated", "Chain 1 CDR3 End Curated",
                                  "Chain 1 CDR3 Curated", "Curated Chain 1 V Gene", "Curated Chain 1 J Gene",
                                  "Chain 1 CDR3 Start Calculated", "Chain 1 CDR3 End Calculated",
                                  "Chain 1 CDR3 Calculated", "Calculated Chain 1 V Gene", "Calculated Chain 1 J Gene",
                                  "Chain 2 Type", "Chain 2 Full Sequence",
                                  "Chain 2 CDR3 Start Curated", "Chain 2 CDR3 End Curated",
                                  "Chain 2 CDR3 Curated", "Curated Chain 2 V Gene", "Curated Chain 2 J Gene",
                                  "Chain 2 CDR3 Start Calculated", "Chain 2 CDR3 End Calculated",
                                  "Chain 2 CDR3 Calculated", "Calculated Chain 2 V Gene", "Calculated Chain 2 J Gene",
                                  "Description", "Antigen", "Organism", "Cell_subset")


data_iedb_CD48 <- data_iedb_CD48 %>% unique()

colnames(data_iedb_CD48) <- c("Epitope ID", "Reference ID", "Chain_alpha", "fullSeq_alpha", "startCur_alpha", "endCur_alpha","CDR3cur_alpha", "Vcur_alpha", "Jcur_alpha", 
                         "startCal_alpha", "endCal_alpha","CDR3cal_alpha", "Vcal_alpha", "Jcal_alpha",
                         "Chain_beta", "fullSeq_beta", "startCur_beta", "endCur_beta","CDR3cur_beta", "Vcur_beta", "Jcur_beta", 
                         "startCal_beta", "endCal_beta","CDR3cal_beta", "Vcal_beta", "Jcal_beta",
                         "Epitope", "Antigen", "Organism", "Cell_subset")

data_iedb_CD48$pairing_id <- 1:nrow(data_iedb_CD48)


#assay file: import assay table to recup information about identification of TCR or Epitope
print("Assay table from IEDB loading")

assay_iedb_CD8 <- read_delim(paste0(directory, "IEDB_tcell_assay_table_MHCI.csv"), delim = ",", col_names = F)
assay_iedb_CD4 <- read_delim(paste0(directory, "IEDB_tcell_assay_table_MHCII.csv"), delim = ",", col_names = F)

col_df_assay_CD8 <- data.frame(col = unlist(assay_iedb_CD8[1,]), subcol = unlist(assay_iedb_CD8[2,]))#change colnames
col_df_assay_CD8 <- col_df_assay_CD8 %>% mutate(new_col = paste0(col,"_", subcol))
colnames(assay_iedb_CD8) <- col_df_assay_CD8$new_col
assay_iedb_CD8 <- assay_iedb_CD8[-c(1,2),]

col_df_assay_CD4 <- data.frame(col = unlist(assay_iedb_CD4[1,]), subcol = unlist(assay_iedb_CD4[2,]))
col_df_assay_CD4 <- col_df_assay_CD4 %>% mutate(new_col = paste0(col,"_", subcol))
colnames(assay_iedb_CD4) <- col_df_assay_CD4$new_col
assay_iedb_CD4 <- assay_iedb_CD4[-c(1,2),]

assay_iedb <- rbind(assay_iedb_CD8, assay_iedb_CD4)

assay_iedb <- assay_iedb %>% select(c('Epitope_Epitope ID', 'Reference_Reference ID', 'Reference_PubMed ID',"Assay_Method/Technique", "Assay Antigen_Antigen Object Type")) %>% unique()
colnames(assay_iedb) <- c('Epitope ID', 'Reference ID', 'PubMed_ID', "Antigen_identification", "Antigen_type")

#merge data and assay files
merge_iedb <- merge(x= data_iedb_CD48, y = assay_iedb, by= c('Epitope ID', 'Reference ID'), all = T) %>% unique()

#separate alpha /beta chain to avoid duplication in filters
merge_iedb <- merge_iedb %>% select(-`Epitope ID`) %>% unique()

data_iedb_alpha <- merge_iedb %>% select(c("Chain_alpha", "fullSeq_alpha", "startCur_alpha", "endCur_alpha","CDR3cur_alpha", "Vcur_alpha", "Jcur_alpha", 
                                           "startCal_alpha", "endCal_alpha","CDR3cal_alpha", "Vcal_alpha", "Jcal_alpha", 
                                           "Epitope", "Antigen", "Organism", "Cell_subset","PubMed_ID", "pairing_id", "Antigen_identification", "Antigen_type"))
colnames(data_iedb_alpha) <- c("Chain", "fullSeq", "startCur", "endCur","CDR3cur", "Vcur", "Jcur", "startCal", "endCal","CDR3cal", "Vcal", "Jcal", 
                               "Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID","pairing_id", "Antigen_identification", "Antigen_type")

data_iedb_beta <- merge_iedb %>% select(c("Chain_beta", "fullSeq_beta", "startCur_beta", "endCur_beta","CDR3cur_beta", "Vcur_beta", "Jcur_beta", 
                                          "startCal_beta", "endCal_beta","CDR3cal_beta", "Vcal_beta", "Jcal_beta",
                                          "Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID","pairing_id", "Antigen_identification", "Antigen_type"))

colnames(data_iedb_beta) <- c("Chain", "fullSeq", "startCur", "endCur","CDR3cur", "Vcur", "Jcur", "startCal", "endCal","CDR3cal", "Vcal", "Jcal", 
                               "Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID","pairing_id", "Antigen_identification", "Antigen_type")


data_iedb <- rbind(data_iedb_alpha, data_iedb_beta) #rbind alpha and beta


### Curation of IEDB accross different cases
print("processing filters")
# a- Keep pairings where we have at least the full sequence and the calculated sequence
data_iedb_curFull <- data_iedb %>% filter(CDR3cal != "" & fullSeq != "")
#Assign a "New curated" sequence from the full sequence and the start and end position of calculated sequence
data_iedb_curFull <- data_iedb_curFull %>% mutate(CDR3curMe = ifelse(fullSeq != "" & startCal != "" & endCal != "", substr(fullSeq, startCal - 1, endCal + 1), 
                                                                           ifelse(fullSeq != "" & startCur != "" & endCur != "", substr(fullSeq, startCur - 1, endCur + 1), "")))
data_iedb_curFull <- data_iedb_curFull %>% mutate(verified = ifelse(Chain == "alpha", "verified_alpha",
                                                                    ifelse(Chain == "beta", "verified_beta", NA)))


# b- For sequence that don't have the full sequence but a curated and a calculated sequence
data_iedb_curCal <- data_iedb %>% filter(is.na(fullSeq) & CDR3cur != "" & CDR3cal != "")
data_iedb_curCal$aaSeq <- paste("C", data_iedb_curCal$CDR3cal, "F", sep = "")
#Check if aaSeq == CDR3cur to assign the "New curated" sequence otherwise indicate that is no verified
data_iedb_curCal <- data_iedb_curCal %>% mutate(CDR3curMe = CDR3cur)
data_iedb_curCal <- data_iedb_curCal %>% mutate(verified = ifelse(aaSeq == CDR3cur & Chain =="beta", "verified_beta", NA))
data_iedb_curCal <- data_iedb_curCal %>% select(-aaSeq)
#for sequence that C-CDR3cal-F : check if CDR3cur == C-CDR3cal-W -> no sequences


# c- For sequence that have just curated sequence
data_iedb_cur <- data_iedb %>% filter(is.na(fullSeq) & is.na(CDR3cal) & CDR3cur != "")
data_iedb_cur <- data_iedb_cur %>% mutate(CDR3curMe = CDR3cur)
data_iedb_cur <- data_iedb_cur %>% mutate(verified = NA)


# d- For sequence with Full seq and curated sequence
data_iedb_full <- data_iedb %>% filter(!is.na(fullSeq) & is.na(CDR3cal) & CDR3cur != "")

CDR3cur_list <- unlist(data_iedb_full$CDR3cur)
fullSeq_list <- unlist(data_iedb_full$fullSeq)
pos_start_list <- c()
pos_end_list <- c()
for(i in 1:length(CDR3cur_list)){       # for-loop over rows
  if(!is.na(CDR3cur_list[i]) & !is.na(fullSeq_list[i])){
    pos <- motif.find(motif = CDR3cur_list[i], sequence = fullSeq_list[i])
    pos_start_list  <- append(pos_start_list, pos[1])
    pos_end_list  <- append(pos_end_list, pos[nchar(CDR3cur_list[i])])
    
  }else{
    pos_start_list  <- append(pos_start_list, NA)
    pos_end_list  <- append(pos_end_list, NA)
  }
}

data_iedb_full <- data_iedb_full %>% mutate(CDR3_start = pos_start_list,
                                            CDR3_end = pos_end_list)
data_iedb_full <- data_iedb_full %>% mutate(CDR3curMe = ifelse(CDR3_start != -1 & fullSeq != CDR3cur & !is.na(CDR3_start) & !is.na(CDR3_end) , substr(fullSeq, CDR3_start - 1, CDR3_end + 1), CDR3cur))
data_iedb_full <- data_iedb_full %>% mutate(verified = ifelse(CDR3_start != -1 & fullSeq != CDR3cur & !is.na(CDR3_start) & !is.na(CDR3_end) & Chain =="alpha", "verified_alpha", 
                                                              ifelse(CDR3_start != -1 & fullSeq != CDR3cur & !is.na(CDR3_start) & !is.na(CDR3_end) & Chain =="beta", "verified_beta", NA)))
data_iedb_full <- data_iedb_full %>% select(-c(CDR3_start, CDR3_end, CDR3_start, CDR3_end))

print("merge all filters")

## merge 4 dfs
data_iedb_filtre <- rbind(data_iedb_curFull, data_iedb_curCal, data_iedb_cur, data_iedb_full) %>% unique()
data_iedb_filtre_alpha <- data_iedb_filtre %>% filter(Chain == "alpha")
data_iedb_filtre_beta <- data_iedb_filtre %>% filter(Chain == "beta")

colnames(data_iedb_filtre_alpha) <- c("Chain_alpha", "fullSeq_alpha", "startCur_alpha", "endCur_alpha","CDR3cur_alpha", "Vcur_alpha", "Jcur_alpha", 
                                           "startCal_alpha", "endCal_alpha","CDR3cal_alpha", "Vcal_alpha", "Jcal_alpha", 
                                           "Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID","pairing_id", "Antigen_identification", "Antigen_type", "CDR3cur_alphaMe","verified_alpha")


colnames(data_iedb_filtre_beta) <- c("Chain_beta", "fullSeq_beta", "startCur_beta", "endCur_beta","CDR3cur_beta", "Vcur_beta", "Jcur_beta", 
                                          "startCal_beta", "endCal_beta","CDR3cal_beta", "Vcal_beta", "Jcal_beta",
                                          "Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID","pairing_id", "Antigen_identification", "Antigen_type", "CDR3cur_betaMe","verified_beta")

data_iedb_curated <- merge(data_iedb_filtre_alpha, data_iedb_filtre_beta, by = c("Epitope", "Antigen", "Organism", "Cell_subset", "PubMed_ID", "pairing_id", "Antigen_identification", "Antigen_type"), all = T)

#Genes harmonization
data_iedb_curated <- data_iedb_curated %>% mutate(Vcur_alphaMe = ifelse(Vcal_alpha != "", Vcal_alpha, Vcur_alpha),
                                                  Jcur_alphaMe = ifelse(Jcal_alpha != "", Jcal_alpha, Jcur_alpha),
                                                  Vcur_betaMe = ifelse(Vcal_beta != "", Vcal_beta, Vcur_beta),
                                                  Jcur_betaMe = ifelse(Jcal_beta != "", Jcal_beta, Jcur_beta))

data_iedb_curated$Vcur_alphaMe <- gsub(x= data_iedb_curated$Vcur_alphaMe, pattern ="TCR", replacement ="TR")
data_iedb_curated$Jcur_alphaMe <- gsub(x= data_iedb_curated$Jcur_alphaMe, pattern ="TCR", replacement ="TR")
data_iedb_curated$Vcur_betaMe <- gsub(x= data_iedb_curated$Vcur_betaMe, pattern ="TCR", replacement ="TR")
data_iedb_curated$Jcur_betaMe <- gsub(x= data_iedb_curated$Jcur_betaMe, pattern ="TCR", replacement ="TR")


data_iedb_curated <- data_iedb_curated %>% mutate(V_alpha= gsub("\\*0.*","\\1",Vcur_alphaMe),
                                                  J_alpha= gsub("\\*0.*","\\1",Jcur_alphaMe),
                                                  V_beta= gsub("\\*0.*","\\1",Vcur_betaMe),
                                                  J_beta= gsub("\\*0.*","\\1",Jcur_betaMe))
data_iedb_curated <- data_iedb_curated %>% select(CDR3cur_alphaMe, V_alpha, J_alpha, CDR3cur_betaMe, V_beta, J_beta, Epitope, Antigen, Organism, Cell_subset, PubMed_ID, Antigen_identification, Antigen_type, pairing_id, verified_alpha, verified_beta)
colnames(data_iedb_curated) <- c("CDR3_alpha", "V_alpha", "J_alpha", "CDR3_beta", "V_beta", "J_beta", "Epitope", "Antigen", "Antigen_organism", "Cell_subset", "PubMed_ID","Antigen_identification", "Antigen_type", "pairing_id", "verified_alpha", "verified_beta")


##verified score
data_iedb_curated$verified_alpha_beta <- paste(data_iedb_curated$verified_alpha, data_iedb_curated$verified_beta)
data_iedb_curated <- data_iedb_curated %>% mutate(Verified_score = ifelse(is.na(CDR3_alpha) & is.na(CDR3_beta) | verified_alpha_beta == "NA NA",0,
                                                                          ifelse(is.na(CDR3_alpha) & !is.na(CDR3_beta) | verified_alpha_beta == "NA verified_beta",1.2,
                                                                                 ifelse(!is.na(CDR3_alpha) & is.na(CDR3_beta) | verified_alpha_beta == "verified_alpha NA",1.1,
                                                                                        ifelse(!is.na(CDR3_alpha) & !is.na(CDR3_beta) | verified_alpha_beta == "verified_alpha verified_beta",2, NA)))))

data_iedb_curated <- data_iedb_curated %>% select(-c(verified_alpha, verified_beta, verified_alpha_beta))
data_iedb_curated <- data_iedb_curated %>% relocate(c(CDR3_alpha, V_alpha, J_alpha, CDR3_beta, V_beta, J_beta), .before = Epitope) 
data_iedb_curated <- subset(data_iedb_curated, nchar(as.character(data_iedb_curated$CDR3_alpha))<=23 & nchar(as.character(data_iedb_curated$CDR3_alpha))>=6 | nchar(as.character(data_iedb_curated$CDR3_beta))<=23 & nchar(as.character(data_iedb_curated$CDR3_beta))>=6)

data_iedb_curated <- data_iedb_curated %>% mutate(Database = "IEDB") %>% unique()
data_iedb_curated$Antigen_identification[data_iedb_curated$Antigen_identification == "ICS"] <- "Intracellular Cytokine Staining (ICS)" 
data_iedb_curated <- data_iedb_curated %>% mutate(Sequencing_method = NA)
data_iedb_curated <- data_iedb_curated %>% relocate(Sequencing_method, .before = Antigen_identification)

data_iedb_curated <- data_iedb_curated %>% select(-pairing_id)
data_iedb_curated <- unique(data_iedb_curated)

##########################################################################################################################################################
rm(data_iedb_CD4, data_iedb_CD48, data_iedb_CD8, col_df_assay_CD8, col_df_assay_CD4, assay_iedb_CD8, assay_iedb_CD4, assay_iedb,                         #
   merge_iedb,data_iedb_alpha, data_iedb_beta, data_iedb, data_iedb_curFull, data_iedb_curCal, data_iedb_cur, data_iedb_full, CDR3cur_list, fullSeq_list,#
   pos_start_list, pos_end_list, data_iedb_filtre, data_iedb_filtre_alpha, data_iedb_filtre_beta, i, pos)#remove intermediate df                         #
##########################################################################################################################################################
