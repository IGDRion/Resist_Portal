# Project plan for Resist_Portal Rshiny tool

11/07/2024
Victor Le Bars
Aurore Besson

## 11/07/2024: First meeting :D

First tasks:

    - [x] initiate tool on github
    - [ ] visual plan for the tool with options/input data/output data/modules to create
    - [x] create module template 
    - [x] filter summary table with a search button + save a variable of the search for further data manipulation on other future files
    - [ ] input files modification ? g
        * gene/tx summary to merge ? + NA values not visible + geneID/geneName on the file 
        * RDS file ? 
        * check isoformSwitchAnalyzer files to use it for DGE/DTE/DTU (no DTU for glioblastoma ?!)
    - [ ] start with only one gene as input + only on cdna data + only expression by sample


**Objectives:**

* create Rshiny app to visualise gene/isoform expression/differential expression based on lncrna_resist project data
* input: gene name/ID or list of gene name/ID
* output: expression in 48 samples
    * 4 cancers: melanoma, glioblastoma, lung cancer, prostate cancer
    * 2 conditions: control vs cancer
    * 2 protocols: cdna, drna
    * gene expression and transcript expression
    * DGE/DTE/DTU
* options: 
    * possibility to create a summary in pdf with selected results
    * create UCSC link to visualise selected gene position on browser
* Important: focus on biotype (not only PCG) and isoforms to add value to the tool (notre outil est mieux!)
