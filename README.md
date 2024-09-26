# Resist Portal

## Description
Resist Portal is a portal made with Rshiny that allows to explore data generated during the lncrna_resist project. 
It allows to explore Differential Gene Expression (DGE), Differential Transcript Expression (DTE), and Differential Transcript Usage (DTU).

The app is separated in multiple menus:
- Summary: Information on genes/transcripts: positions on the genome, IDs & names, biotype, and more.
- Count: Read count per transcript for each sample. Counts are normalized based on TPM (transcripts per million).
- DGE: Differential Gene Expression data. 
- DTE: Differential Transcript Expression data.
- DTU: Differential Transcript Usage data.

On every menu, a search bar is here to search a gene (by its name or its gene ID (Ensembl)) get more detail on it.
When searching for a gene, a "Query" tab appears to get more specific information on the searched gene.

## Access to the application & Installation
A ready-to-use version of the application is available on here.
It can also be installed from gitHub and launched with the command: `runApp('path/to/Resist_Portal')`, or with the graphical interface of R studio through the app.R file.
