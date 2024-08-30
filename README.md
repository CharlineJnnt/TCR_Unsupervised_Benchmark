# Benchmarking unsupervised methods for inferring TCR specificity
This repository contains scripts and data to accompagny the study "Benchmarking unsupervised methods for inferring TCR specificity".

## Introduction
Identifying T cell receptor (TCR) specificity is crucial for advancing the understanding of adaptive immunity. Despite the development of computational methods to infer TCR specificity, their performance has not been thoroughly compared. We addressed this by curating a unified database of 190,670 human TCRs with known specificities for 2,313 epitopes across 121 organisms, combining data from IEDB ([iedb](https://www.iedb.org/)), McPAS-TCR ([McPAS](https://friedmanlab.weizmann.ac.il/McPAS-TCR/)), and VDJdb ([VDJdb](https://vdjdb.cdr3.net/)). Nine methods for clustering TCRs based on similarity were benchmarked against this dataset. DeepTCR demonstrated the best retention, while ClusTCR, TCRMatch, and GLIPH2 excelled in cluster purity but had lower retention. DeepTCR, Levenshtein distance, and TCRdist3 generated large clusters, with DeepTCR showing high purity. Smaller, antigen-specific clusters were produced by Hamming distance, GIANA, and iSMART. GLIPH2 and DeepTCR were the most sensitive in capturing antigen-specific TCRs. This study offers a unified TCR database and a benchmark of specificity inference methods, guiding researchers in selecting appropriate tools.

## Worflow
![](Workflow_study.png)

## Analysis
This repository as the paper is composed of two parts: the process of the curated pooled database and the benchmark analysis of clustering methods. This repository is organized as follow: 


```bash
C:.
│   .RData
│   .Rhistory
│   README.md
│   Workflow_study.png
│
├───Analysis
│       .Rhistory
│       benchmark_tools_pooled_database_CD8.Rmd
│       output_processing.R
│
├───Data
│   ├───Database
│   │       database_pooled_human_2023_03_15.txt
│   │       IEDB_tcell_assay_table_MHCI.csv
│   │       IEDB_tcell_assay_table_MHCII.csv
│   │       IEDB_tcell_receptor_table_MHCI.csv
│   │       IEDB_tcell_receptor_table_MHCII.csv
│   │       McPAS-TCR.csv
│   │       vdjdb.txt
│   │
│   └───Output_methods
│           clusTCR_output.csv
│           clusTCR_output_edgelist.txt
│           clusTCR_output_summary.csv
│           DeepTCR_output.csv
│           GIANA_output.txt
│           GLIPH2_output_cluster.csv
│           iSMART_output.txt
│           TCRdist3_hclust_output.csv
│           TCRMatch_output.txt
│
├───Functions
│       .Rhistory
│       alignment_functions.R
│       database_processing.R
│       general_functions.R
│       network_functions.R
│       purity_functions.R
│       similarity_functions.R
│
├───Plot
│       6EpitopeSelection.pdf
│       AllMetricsTh_1.pdf
│       AllMetricsTh_1_byChain.pdf
│       AllSpecificitiesInSpeClusters_GIL.pdf
│       clusterSizeDistribution.pdf
│       correlationPurityRetention.pdf
│       NumberSpecificities_GIL.pdf
│       PercClusteredSeq_GIL.pdf
│       PercClusters_GIL.pdf
│       purityChain.pdf
│       SensibilityAllEpitopes.pdf
│       Sensibility_GIL.pdf
│
└───Pooled_Database
        database_processing_human.R
        description_database.R
        iedb_curation.R
        McPAS_curation.R
        vdj_curation.R
```

### 1. Process the pooled database
Scripts to process the 3 public databases are available in the /Pooled_Database/ folder.
- iedb_curation.R: Curation of the IEDB database.
- McPAS_curation.R: Curation of the McPAS database.
- vdj_curation.R: Curation of the VDJdb database.
- database_processing_human.R: Imports the 3 curated databases and assemble them into a single pooled database.
All original database and final pooled human database are provided in the /Data/Database/ folder.
- description_database.R: Create descriptive plots for the pooled human database saved in the /Plot/ folder (referring to the Figure 1 of the study).

### 2. Perform the benchmarking analysis
- (main) /Analysis/benchmark_tools_pooled_database_CD8.Rmd: Rmarkdown file which prove all analyses to compare the clustering methods. Results plots are stored in the /Plot/ folder (referring to other Figures of the paper).
- /Analysis/output_processing.R: Reads and processes the output of clustering methods stored in the /Data/Output_methods/ folder. 
- /Functions/ folder : Contains the database_processing.R script that processes and filter the pooled human database for the analysis and all other functions useful for the benchmark.

## Data acquisition
WARNING: Data from public databases were downloaded between January and March 2023. So please keep in mind that if you use a more recent version of these databases, the final pooled human database resulting from curation may be different. The original data used in this analysis is provided in this repository, but some versions of the R package may have been modified since the database_processing_human.R script was created and may also lead to a different human grouped database. 

Steps to download public databases: 
- Data from McPAS-TCR database were directly downloaded from the home page of the web site at https://friedmanlab.weizmann.ac.il/McPAS-TCR/
- Data from VDJdb were uploaded from the original web site at https://vdjdb.cdr3.net/about
- For IEDB data, we have started from the home page https://www.iedb.org/ and selected the following criteria:
                - Epitope:"Any"
                - Assay:"T Cell" and Outcome:"Positive"
                - Epitope source:no filter
                - MHC Restriction:"Class I"
                - Host:"human"
                - Disease:"Any"
Once the page was loaded, the filter TCR_type: "has TCR sequence" was added. Then, we have downloaded the Assays and Receptors tables. Finally, we have re-performed this manip using MHC-Restriction:"Class II" resulting in 4 files.

## Figures
All plots of the paper are available in the /Plot/ folder of this repository.

## Support
For any questions regarding the analysis, bug reports or feature requests, please submit an issue. 