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
[DTUall.qs](https://homeandco.genouest.org/api/download/groups/igdrion/lncrna_resist_cgo/secondary/cdna_SQK-DCS109/INPUT_RESIST_PORTAL/DTUall.qs?token=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJ1c2VyIjp7IlVpZCI6InZsZWJhcnMiLCJQYXNzd29yZCI6IiIsIlVpZE51bWJlciI6NTkyNjUsIkdpZE51bWJlciI6NDAzOTksIkhvbWUiOiIvaG9tZS9nZW5vdWVzdC9jbnJzX3VtcjYyOTAvdmxlYmFycyIsIkdyb3VwcyI6W3siR2lkTnVtYmVyIjo0MDM5OSwiTmFtZSI6ImNucnNfdW1yNjI5MCJ9LHsiR2lkTnVtYmVyIjo0MDk0OSwiTmFtZSI6InByal9pZ2RyaW9uIn1dLCJBZG1pbiI6ZmFsc2UsIlNoYXJlcyI6W3siUGF0aCI6Ii9ncm91cHMiLCJSZWFkV3JpdGUiOnRydWUsIlNoYXJlZEJ5IjpudWxsLCJGb3JiaWRzIjpudWxsfV19LCJuYW1lIjoiIiwiZXhwIjoxNzI5NTM4NDAzLCJpc3MiOiJob21lYW5kY28ifQ.EEBAvG_S1DafVWg6Jdi1xo622Bw7vt9UMsqYn9uA52I) 
