# Resist_Portal

## Description

Resist Portal is a Rshiny portal for the exploration of long read sequencing data generated during the lncrna_resist project. The portal includes modules for explorating Differential Gene Expression (DGE), Differential Transcript Expression (DTE), and Differential Transcript Usage (DTU).

The app is separated in multiple menus:

- Summary: Information on genes/transcripts: positions on the genome, IDs & names, biotype, new (discovered by Bambu) or known (present in the annotation file) genes, filters.
- Count: Read count per transcript for each sample. Counts are normalized based on TPM (transcripts per million).
- DGE: Differential Gene Expression data.
- DTE: Differential Transcript Expression data.  
DGE and DTE were performed with DESeq2 per cancer using a negative binomial generalized linear model. Significance of differential gene/transcript expression is defined by p-value adjusted (`padj`) and log 2 fold change (`log2FoldChange`).  
- DTU: Differential Transcript Usage data performed with isoformSwitchAnalyzeR. A significant isoform switch is defined by `alpha` corresponding to FDR corrected P-value cut-off and `dIF` reflecting the effect size.

On every menu, a search bar is here to search a gene (by its name or its gene ID (Ensembl)) get more detail on it. When searching for a gene, a "Query" tab appears to get more specific information on the searched gene.


## Access to the application & Installation

A ready-to-use version of the application is available on [here](https://shiny-dog.univ-rennes.fr/Resist_Portal/). It can also be installed from gitHub and launched with the command: runApp('path/to/Resist_Portal'), or with the graphical interface of R studio through the app.R file. 

If you choose to install from gitHub, download the following data file in the `DATA`repository: 
[DTUall.qs](https://shiny-dog.univ-rennes.fr/dwnd_data/DTUall.qs) 
