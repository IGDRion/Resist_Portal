library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(IsoformSwitchAnalyzeR)
library(patchwork)
library(plotly)
library(purrr)
library(reshape2)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(stringr)
library(summaryBox)
library(tidyr)

# Bootstrap 4
theme <- bslib::bs_theme(version = 4)

# Set directory to the src directory of the app so that the path
# is correctly resolved later (otherwise they won't be found)
setwd(getSrcDirectory(function(){})[1])
# Source all modules
source("./module_search.R")
source("./module_boxplot.R")
source("./module_barplot.R")
source("./module_filters.R")
source("./module_volcano.R")
source("./module_statbox.R")
source("./module_switchplot.R")

# Ui
ui <- page_navbar(
  
  # Set up shinyjs to use its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  
  title = "Resist Portal",
  bg = "#70b684",
  inverse = TRUE,
  nav_panel(title = "Summary", 
            p(
              HTML(
                paste(
                  "Resist Portal is an R shiny application developed to visualize data from the lncRNA Resist Project. Samples, coming from 4 cancers, have been analysed by long-read sequencing (nanopore technology). Resist Portal provides an overview of lncRNA resist analysis in terms of expression (tpm) and differential analysis at gene and transcript levels for a specific gene.<br><br>",
                  "Conditions are the following:",
                  gsub("- ", "<br>- ", "- 4 cell lines: Melanoma, Glioblastoma, Lung Cancer, Prostate Cancer - 2 condtions: resistant and sensitive to treatment - triplicates - cdna and drna protocoles (only cdna results shown)"),
                  sep = ""
                )
              )
            ), 
            # Image of experimental design
            HTML('<center><img(src = "experimental_design.png", alt = "Experimental Design", width = "50%", height = "auto"></center>'),
            # Search bar 
            searchBarUI("searchBar1", "submit_btn", "reset_btn"),
            # Main summary table
            DTOutput(outputId = "SummaryTable") %>% withSpinner(), # withSpinner() display a loading spinner while the dataframe is not on screen
            # download table
            downloadButton("downloadSummary", "Download Summary Table"),
            # UCSC link
            uiOutput("UCSClink")
            ),
  
  nav_panel(title = "Count",
            p(""),
            # Write which gene is currently selected at the top of the page
            uiOutput("SelectedGeneTextCount"),
            
            # Search bar 
            searchBarUI("searchBar2", "submit_btn", "reset_btn"),
            p("Expressions are normalized based on TPM (transcripts per million)."),
            
            tabsetPanel(id = "TabsetCount",
                        
                        tabPanel(title = "All",
                                 # Count table
                                 DTOutput(outputId = "CountTableFull") %>% withSpinner()
                                 ),
                        tabPanel(title = "Query Gene Level",
                                 # Count table
                                 DTOutput(outputId = "CountTable") %>% withSpinner(),
                                 # Boxplot
                                 boxplotUI(id = "boxplotGene")),
                        tabPanel(title = "Query Transcript Level",
                                 # Count table + Barplot
                                 DTOutput(outputId = "CountTableTx") %>% withSpinner(),
                                 barplotUI(id = "barplotTX") %>% withSpinner(),
                                 # Boxplot
                                 boxplotUI(id = "boxplotTX")
                                 )
                        )
  ),
  nav_panel(title = "DGE",
            p(""),
            # Write which gene is currently selected at the top of the page
            uiOutput("SelectedGeneTextDGE"),
            
            tabsetPanel(id= "TabsetDGE",
                        tabPanel(title = "All",
                                 p(""),
                                 
                                 # Search bar + Filter box
                                 fluidRow(
                                   column(width = 6, align = "center",
                                          filtersBoxUI("filtersDGE", Dtype = "DE")
                                   ), 
                                   column(width = 6, align = "center",
                                          searchBarUI("searchBar3", "submit_btn", "reset_btn"),
                                          # Box test
                                          statboxUI("statboxDGEAll")
                                   )
                                 ),
                                 
                                 # DGE table All
                                 DTOutput(outputId = "DGETableAll") %>% withSpinner(),
                                 
                                 # download table
                                 downloadButton("downloadDGEall", "Download all DGE Table")
                                        
                        ),
                        tabPanel(title = "Query",
                                 p(""),
                                 
                                 # Search bar + a slider for padj threshold
                                 fluidRow(
                                   column(6, align = "right",
                                     sliderTextInput(inputId = "padj_threshold_DGEQuery",
                                                     label = "padj Threshold:",
                                                     choices = c(0.01, 0.05, "NONE"),
                                                     selected = "NONE",
                                                     grid = TRUE)
                                   ),
                                   column(6, align = "left",
                                          searchBarUI("searchBar4", "submit_btn", "reset_btn"))
                                 ),
                                 
                                 # Box test
                                 statboxUI("statboxDGEQuery"),
                                 
                                 # DGE table Query
                                 DTOutput(outputId = "DGETableQuery") %>% withSpinner(),
                                 
                                 # Volcano plot
                                 volcanoUI("volcanoDGE"),
                                 
                                 # Text to explain that volcano plot points are stopped when they are above a certain limit
                                 uiOutput("volcanoText")
                        )
            )
  ),
  
  nav_panel(title = "DTE",
            p(""),
            # Write which gene is currently selected at the top of the page
            uiOutput("SelectedGeneTextDTE"),
            
            tabsetPanel(id = "TabsetDTE",
                        tabPanel(title = "All",
                                 # Search bar + Filter box
                                 fluidRow(
                                   column(width = 6, align = "center",
                                          filtersBoxUI("filtersDTE", Dtype = "DE")
                                   ), 
                                   column(width = 6, align = "center",
                                          searchBarUI("searchBar5", "submit_btn", "reset_btn"))
                                   
                                 ),
                                 # DTE table
                                 DTOutput(outputId = "DTETableAll") %>% withSpinner(),
                                 
                                 # download table
                                 downloadButton("downloadDTEall", "Download all DTE Table")
                        ),
                        tabPanel(title = "Query",
                                 p(""),
                                 
                                 # Search bar + a slider for padj threshold
                                 fluidRow(
                                   column(6, align = "right",
                                     sliderTextInput(inputId = "padj_threshold_DTEQuery",
                                                     label = "padj Threshold:",
                                                     choices = c(0.01, 0.05, "NONE"),
                                                     selected = "NONE",
                                                     grid = TRUE)
                                   ),
                                   column(6, align = "left",
                                          searchBarUI("searchBar6", "submit_btn", "reset_btn"))
                                 ),
                                 
                                 # Box test
                                 #statboxUI("statboxDTEQuery"),
                                 
                                 # DGE table Query
                                 DTOutput(outputId = "DTETableQuery") %>% withSpinner(),
                                 
                                 # download table
                                 downloadButton("downloadDTEquery", "Download query DTE Table"),
                                 
                                 # Volcano plot
                                 volcanoUI("volcanoDTE"),
                                 
                                 # Text to explain that volcano plot points are stopped when they are above a certain limit
                                 uiOutput("volcanoText")
                        )
            )
            
  ),
  
  nav_panel(title = "DTU",
            p(""),
            # Write which gene is currently selected at the top of the page
            uiOutput("SelectedGeneTextDTU"),
            
            tabsetPanel(id= "TabsetDTU",
                        tabPanel(title = "All",
                                 # Search bar + Filter box
                                 fluidRow(
                                   column(width = 6, align = "center",
                                          filtersBoxUI("filtersDTU", Dtype = "DU")
                                   ),
                                   column(width = 6, align = "center",
                                          searchBarUI("searchBar7", "submit_btn", "reset_btn"))
                                 ),
                                 # DTU table
                                 DTOutput(outputId = "DTUTableAll") %>% withSpinner(),
                                 # download table
                                 downloadButton("downloadDTUall", "Download all DTU Table")
                        ),
                        tabPanel(title = "Query",
                                 p(""),
                                 fluidRow(
                                   column(width = 12, align = "center",
                                          searchBarUI("searchBar8", "submit_btn", "reset_btn"))
                                 ),
                                 # DTU table
                                 DTOutput(outputId = "DTUTableQuery") %>% withSpinner(),
                                 # download table
                                 downloadButton("downloadDTUquery", "Download query DTU Table"),
                                 # switch plot
                                 switchPlotUI(id = "switchplot") %>% withSpinner()
                                 )
            )
  ),
  
  nav_spacer(),
  nav_item(
    actionButton(
      inputId = "infoButton",
      label = "Info",
      icon = icon("info-circle"),
      style = "background-color: white; color: #70b684; border: none;"
    )
  ),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("GitHub", href = "https://github.com/IGDRion/Resist_Portal", target = "_blank")),
    nav_item(tags$a("Shiny", href = "https://shiny.posit.co", target = "_blank")),
    nav_item(tags$a("DESeq2", href = "https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html", target = "_blank")),
    nav_item(tags$a("ISA", href = "https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html", target = "_blank"))
  )
)