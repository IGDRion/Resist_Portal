if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("IsoformSwitchAnalyzeR")) BiocManager::install("IsoformSwitchAnalyzeR")
if (!require("patchwork")) devtools::install_github("thomasp85/patchwork")
if (!require("plotly")) install.packages("plotly")
if (!require("purrr")) install.packages("purrr")
if (!require("remotes")) install.packages("remotes")
if (!require("reshape2")) install.packages("reshape2")
if (!require("shiny")) install.packages("shiny")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("stringr")) install.packages("stringr")
if (!require("summaryBox")) remotes::install_github("deepanshu88/summaryBox")
if (!require("tidyr")) install.packages("tidyr")


library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(ggpubr)
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
                                 DTOutput(outputId = "DGETableAll") %>% withSpinner()
                                        
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
                                 DTOutput(outputId = "DTETableAll") %>% withSpinner()
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
                                 DTOutput(outputId = "DTUTableAll") %>% withSpinner()),
                        tabPanel(title = "Query",
                                 # DTU table
                                 DTOutput(outputId = "DTUTableQuery") %>% withSpinner(),
                                 # switch plot
                                 switchPlotUI(id = "switchplot1") %>% withSpinner(),
                                 switchPlotUI(id = "switchplot2") %>% withSpinner(),
                                 switchPlotUI(id = "switchplot3") %>% withSpinner(),
                                 switchPlotUI(id = "switchplot4") %>% withSpinner())
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

# Server
server <- function(input, output, session) {
  
  # Import datasets
  summary_data <- read.table("./DATA/summary_lncRNAresist.txt", header = TRUE)
  autocomplete_list <- unique(c(summary_data$gene_id, summary_data$gene_name))

  count_data <- readRDS("./DATA/resist_transcript_expression.rds")
  
  load("./DATA/Differential_analysis.RData") #Directly load into the environment DGEall, DTEall
  DGEall <- DGEall
  DTEall <- DTEall

  switch_data <- readRDS("./DATA/All_switchlist_DEXSeq.Rds")
  DTUall <- switch_data$isoformFeatures %>%
    select(isoform_id, gene_id, condition_1, condition_2, gene_name, 
         gene_biotype, iso_biotype, IF_overall, IF1,IF2, dIF, isoform_switch_q_value) %>%
    rename(transcript_id = isoform_id, transcript_biotype = iso_biotype) %>%
    mutate(cancer = case_when(grepl("glioblastoma", condition_1)~"Glioblastoma",
                             grepl("lung", condition_1)~ "Lung",                        
                             grepl("melanoma", condition_1)~ "Melanoma",
                             grepl("prostate", condition_1)~ "Prostate"))
  
  ############################
  ###  SEARCH TERM VALUE ###
  ############################
  
  # Initialize all search bars
  searchBarServer("searchBar1", autocomplete_list)
  searchBarServer("searchBar2", autocomplete_list)
  searchBarServer("searchBar3", autocomplete_list)
  searchBarServer("searchBar4", autocomplete_list)
  searchBarServer("searchBar5", autocomplete_list)
  searchBarServer("searchBar6", autocomplete_list)
  searchBarServer("searchBar7", autocomplete_list)
  
  # Create a reactive value to store the gene searched through the researched bar, so it can be used everywhere in the app
  search_term <- reactiveVal("")
  
  # Add this observeEvent so when the search bar module is used, search_term changes
  observeEvent(input$submit_btn, {
    req(any(
      input$searchBar1 != "",
      input$searchBar2 != "",
      input$searchBar3 != "",
      input$searchBar4 != "",
      input$searchBar5 != "",
      input$searchBar6 != "",
      input$searchBar7 != ""
    ))
    
    search_inputs <- c(input$searchBar1, input$searchBar2, input$searchBar3,
                       input$searchBar4, input$searchBar5, input$searchBar6,
                       input$searchBar7)
    new_value <- search_inputs[search_inputs != ""][1]
    
    # Search term becomes the value searched in whichever search bar
    search_term(new_value)
    
    # Reset search bars to prepare next research
    searchBarServer("searchBar1", autocomplete_list)
    searchBarServer("searchBar2", autocomplete_list)
    searchBarServer("searchBar3", autocomplete_list)
    searchBarServer("searchBar4", autocomplete_list)
    searchBarServer("searchBar5", autocomplete_list)
    searchBarServer("searchBar6", autocomplete_list)
    searchBarServer("searchBar7", autocomplete_list)
    
  })
  
  
  # Hide all Query tabs upon initialization while no gene as been searched
  hideTab(inputId = "TabsetCount", target = "Query Gene Level")
  hideTab(inputId = "TabsetCount", target = "Query Transcript Level")
  
  hideTab(inputId = "TabsetDGE", target = "Query")
  
  hideTab(inputId = "TabsetDTE", target = "Query")
  
  hideTab(inputId = "TabsetDTU", target = "Query")
  
  
  # Add this observeEvent so when the reset button of the search bar module is used, search_term becomes NULL and filtered_data becomes the full dataset again
  observeEvent(input$reset_btn, {
    updateTabsetPanel(session, "TabsetCount", selected = "All")
    updateTabsetPanel(session, "TabsetDGE", selected = "All")
    updateTabsetPanel(session, "TabsetDTE", selected = "All")
    updateTabsetPanel(session, "TabsetDTU", selected = "All")

    isolate({
      # Reset search term
      search_term("")
      
      # Remove all query tabs when reset button pressed
      hideTab(inputId = "TabsetCount", target = "Query Gene Level")
      hideTab(inputId = "TabsetCount", target = "Query Transcript Level")
      
      hideTab(inputId = "TabsetDGE", target = "Query")
      
      hideTab(inputId = "TabsetDTE", target = "Query")
      
      hideTab(inputId = "TabsetDTU", target = "Query")
    })
  })
  
  observe({
    if (search_term() != "") {
      showTab(inputId = "TabsetCount", target = "Query Gene Level")
      showTab(inputId = "TabsetCount", target = "Query Transcript Level")
      
      showTab(inputId = "TabsetDGE", target = "Query")
      
      showTab(inputId = "TabsetDTE", target = "Query")
      
      showTab(inputId = "TabsetDTU", target = "Query")
    }
  })
  
  
  ###########################
  ###  TABSET 1 : SUMMARY ###
  ###########################
  
  
  # Create a reactive value to store the filtered summary_data, changing if search_term is not NULL
  filtered_summary_data <- reactive({
    if (search_term() != "") {
      summary_data %>% 
        filter(gene_id == search_term() | gene_name == search_term())
    } else {
      summary_data  # Return full summary_data if search is empty
    }
  })
  
  # Render the table, with filtered_summary_data (by default full table, filtered only if a gene has been searched (e.g search_term is not NULL))
  output$SummaryTable <- renderDT({
    datatable(filtered_summary_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  
  UCSC_url <- reactive({
    if (search_term() != ""){
      term <- search_term()
      chr <- filtered_summary_data()$seqnames
      start <- filtered_summary_data()$start
      end <- filtered_summary_data()$end
      a("View on UCSC with cDNA and dRNA data", href="https://genome-euro.ucsc.edu/s/abesson/lncRNAresist_hg38_Ensembl108_cDNA_dRNA", target = "_blank")
    } else {
      ""
    }
  })
  
  output$UCSClink <- renderUI({
    UCSC_url()
  })
  
  #########################
  ###  TABSET 2 : COUNT ###
  #########################

  # Text to display which gene has been searched and giving info to click on the query sub-tab.
  output$SelectedGeneTextCount <- renderUI({
    if (search_term() != ""){
      paste0("Currently selected gene: ", search_term(), ". Please click on the \"Query\" tabs (at gene or transcript level) to get detailed information on it.")
    } else {
      ""
    }
  })
  
  
  ### PREPARE DATAFRAME (FILTERING WITH SEARCH TERM) AND DISPLAY THEM IN THE OUTPUT ###
  
  # Create a reactive value to store the count data, filtering it if a gene has been searched
  filtered_count_data <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) %>% # If a gene has been searched, find and take only the row with the name of the gene
        group_by(gene_id, gene_name) %>% # Group by gene
        summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop') # Sum all numeric columns
    } else {
      NULL # If no gene has been searched, becomes NULL
    }
  })
  
  filtered_count_data_tx <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) %>% # If a gene has been searched, find and take only the row with the name of the gene
        as_tibble() # to match filtered_count_data() that went tibble when using group_by and summarize
    } else {
      NULL  # If no gene has been searched, becomes NULL
    }
  })
  
  # renaming columns for homogeneity
  output$CountTableFull <- renderDT({
    old_names <- c("501Mel_1_S", "501Mel_2_S", "501Mel_3_S",
                   "501Mel_1_R", "501Mel_2_R", "501Mel_3_R",
                   "ADCA72_1_S", "ADCA72_2_S", "ADCA72_3_S",
                   "ADCA72_1_R", "ADCA72_2_R", "ADCA72_3_R",
                   "PC3_1_S", "PC3_2_S", "PC3_3_S",
                   "PC3_1_R", "PC3_2_R", "PC3_3_R",
                   "U251_1_S", "U251_2_S", "U251_3_S",
                   "U251_1_R", "U251_2_R", "U251_3_R")
    new_names <- c("Melanoma_1_S", "Melanoma_2_S", "Melanoma_3_S",
                   "Melanoma_1_R", "Melanoma_2_R", "Melanoma_3_R",
                   "Lung_1_S", "Lung_2_S", "Lung_3_S",
                   "Lung_1_R", "Lung_2_R", "Lung_3_R",
                   "Prostate_1_S", "Prostate_2_S", "Prostate_3_S",
                   "Prostate_1_R", "Prostate_2_R", "Prostate_3_R",
                   "Glioblastoma_1_S", "Glioblastoma_2_S", "Glioblastoma_3_S",
                   "Glioblastoma_1_R", "Glioblastoma_2_R", "Glioblastoma_3_R")
    
    DT::datatable(count_data %>%
                    mutate_if(is.numeric, ~round(., digits = 3)) %>%
                    dplyr::rename_with(~ new_names[match(., old_names)],
                                       all_of(old_names))
                  ) 
  })
  
  # Render the table with gene, with filtered_count_data
  # Taking median of the 3 replicates for each cancer/condition, and renaming cell lines to cancer name
  output$CountTable <- renderDT({
    DT::datatable(filtered_count_data() %>%
                    group_by(gene_id, gene_name) %>%
                    summarise(
                      Melanoma_Sensitive = median(c(`501Mel_1_S`, `501Mel_2_S`, `501Mel_3_S`), na.rm = TRUE),
                      Melanoma_Resistant = median(c(`501Mel_1_R`, `501Mel_2_R`, `501Mel_3_R`), na.rm = TRUE),
                      Lung_Sensitive = median(c(ADCA72_1_S, ADCA72_2_S, ADCA72_3_S), na.rm = TRUE),
                      Lung_Resistant = median(c(ADCA72_1_R, ADCA72_2_R, ADCA72_3_R), na.rm = TRUE),
                      Prostate_Sensitive = median(c(PC3_1_S, PC3_2_S, PC3_3_S), na.rm = TRUE),
                      Prostate_Resistant = median(c(PC3_1_R, PC3_2_R, PC3_3_R), na.rm = TRUE),
                      Glioblastoma_Sensitive = median(c(U251_1_S, U251_2_S, U251_3_S), na.rm = TRUE),
                      Glioblastoma_Resistant = median(c(U251_1_R, U251_2_R, U251_3_R), na.rm = TRUE),
                      .groups = 'drop'
                    ) %>%
                    mutate_if(is.numeric, ~round(., digits = 3)),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  
  # Render the table with transcripts, with filtered_count_data_tx
  # Taking median of the 3 replicates for each cancer/condition, and renaming cell lines to cancer name
  output$CountTableTx <- renderDT({
    DT::datatable(filtered_count_data_tx() %>%
                    group_by(transcript_id, gene_id, gene_name) %>%
                    summarise(
                      Melanoma_Sensitive = median(c(`501Mel_1_S`, `501Mel_2_S`, `501Mel_3_S`), na.rm = TRUE),
                      Melanoma_Resistant = median(c(`501Mel_1_R`, `501Mel_2_R`, `501Mel_3_R`), na.rm = TRUE),
                      Lung_Sensitive = median(c(ADCA72_1_S, ADCA72_2_S, ADCA72_3_S), na.rm = TRUE),
                      Lung_Resistant = median(c(ADCA72_1_R, ADCA72_2_R, ADCA72_3_R), na.rm = TRUE),
                      Prostate_Sensitive = median(c(PC3_1_S, PC3_2_S, PC3_3_S), na.rm = TRUE),
                      Prostate_Resistant = median(c(PC3_1_R, PC3_2_R, PC3_3_R), na.rm = TRUE),
                      Glioblastoma_Sensitive = median(c(U251_1_S, U251_2_S, U251_3_S), na.rm = TRUE),
                      Glioblastoma_Resistant = median(c(U251_1_R, U251_2_R, U251_3_R), na.rm = TRUE),
                      .groups = 'drop'
                    ) %>%
                    mutate_if(is.numeric, ~round(., digits = 3)),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  
  
  ### PREPARE DATAFRAME USED IN BOXPLOTS (FROM FILTERED DATAFRAME ABOVE) AND DISPLAY THEM IN THE OUTPUT ###
  
  prepare_table_boxplot <- function(data, type) {
    
    if (type == "gene") {
      # Prepare data for the boxplot with gene_id
      data <- melt(data)
    } else if (type == "tx") {
      # Prepare data for the boxplot with transcript_id
      data <- data %>%
        dplyr::select(-c(gene_id, gene_name))
      data <- melt(data)
    }
    
    data <- data %>% 
      mutate(variable = as.character(variable)) %>%
      # Add new column with cancer type
      mutate(cancer_type = case_when(
        startsWith(variable, "ADCA72") ~ "Lung",
        startsWith(variable, "PC3") ~ "Prostate",
        startsWith(variable, "U251") ~ "Glioblastoma",
        startsWith(variable, "501Mel") ~ "Melanoma"
      )) %>%
      # Add new column with condition information (sensitive/resistant)
      mutate(condition = case_when(
        endsWith(variable, "S") ~ "Sensitive",
        endsWith(variable, "R") ~ "Resistant"
      )) %>%
      mutate(cancer_condition = paste0(cancer_type, "_", condition))
    
    if (type == "gene") {
      data <- data %>%
        group_by(gene_id, variable, cancer_type, condition, cancer_condition) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = 'drop')
    }
    
    print(data)
    return(data)
  }
  
  
  # If search_term() is not NULL (e.g a gene has been searched)
  observe({

    if (search_term() != "") {
      
      shinyjs::show("boxplotGene")
      shinyjs::show("boxplotTX")
      shinyjs::show("barplotTX")
      
      count_boxplot_data_gene <- prepare_table_boxplot(filtered_count_data(), "gene")
      
      count_boxplot_data_tx <- prepare_table_boxplot(filtered_count_data_tx(), "tx")
      
      boxplotServer(id = "boxplotGene",
                    type = "gene",
                    count_boxplot_data_gene,
                    search_term())
      
      boxplotServer(id = "boxplotTX",
                    type = "tx",
                    count_boxplot_data_tx,
                    search_term())
      
      barplotServer(id = "barplotTX",
                    filtered_count_data_tx())
      
    } else {
      
      output$boxplotGene <- renderUI({ NULL })
      output$boxplotTX <- renderUI({ NULL })
      output$barplotTX <- renderUI({ NULL })
    }
    
  })

  #########################
  ###  TABSET 3 : DGE ###
  #########################
  
  # Text to display which gene has been searched and giving info to click on the query sub-tab.
  output$SelectedGeneTextDGE <- renderUI({
    if (search_term() != ""){
      paste0("Currently selected gene: ", search_term(), ". Please click on the \"Query\" tab to get detailed information on it.")
    } else {
      ""
    }
  })
  
  
  filtersDGE <- filtersBoxServer("filtersDGE")
  
  # Create a reactive value to store the DGE data,filtering with the sliders log2fc / padj (on all data)
  filtered_DGE_data_All <- reactive({
    data <- DGEall
    
    # Apply log2FoldChange filter based on DEside
    data <- data %>% 
      filter(case_when(
        filtersDGE$DEside() == "both" ~ abs(log2FoldChange) >= filtersDGE$log2fc_threshold(),
        filtersDGE$DEside() == "up" ~ log2FoldChange >= filtersDGE$log2fc_threshold(),
        filtersDGE$DEside() == "down" ~ log2FoldChange <= -filtersDGE$log2fc_threshold()
      ))
    
    # Apply padj filter if not "NONE"
    if (filtersDGE$padj_threshold() != "NONE") {
      data <- data %>% filter(padj <= as.numeric(filtersDGE$padj_threshold()))
    }
    
    # Apply cancer type filter
    if (!is.null(filtersDGE$cancer_types())) {
      data <- data %>% filter(cancer %in% filtersDGE$cancer_types())
    }
 
    return(data)
  })
  
  
  # Create a reactive value to store the DGE data,filtering to keep only the searched gene
  filtered_DGE_data_Query <- reactive({
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- DGEall
      data <- data %>% filter(gene_id == search_term() | gene_name == search_term())
      return(data)
    } else {
      NULL
    }
    
  })
  
  # Render the filtered DGE table of All genes
  output$DGETableAll <- renderDT({
    datatable(filtered_DGE_data_All() %>%
                mutate(across(where(is.numeric) & !padj, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatStyle(columns = 'cancer',
                  target = 'cell',
                  backgroundColor = styleEqual(
                    c("Melanoma", "Lung", "Prostate", "Glioblastoma"),
                    c("#f5ab05", "#00cc94", "#0084cf", "#eb8dc1")
                  ))
  })
  
  # Render the DGE table for the query gene (4 rows, one per cancer)
  output$DGETableQuery <- renderDT({
    datatable(filtered_DGE_data_Query() %>% 
                mutate(across(where(is.numeric) & !padj, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatStyle(columns = 'cancer',
                  target = 'cell',
                  backgroundColor = styleEqual(
                    c("Melanoma", "Lung", "Prostate", "Glioblastoma"),
                    c("#f5ab05", "#00cc94", "#0084cf", "#eb8dc1")
                  ))
  })
  
  
  # Creating list required for statbox server module + launching it
  statboxDataDGE <- reactive({
    if(search_term() != ""){
      list(
        search_term = search_term(),
        DGEall = DGEall,
        padj_threshold = input$padj_threshold_DGEQuery
      )
    } else {
      list(
        DGEall = DGEall,
        padj_threshold = input$padj_threshold_DGEQuery
      )
    }
  })
  
  statboxServer("statboxDGEAll", statboxDataDGE, "All")
  statboxServer("statboxDGEQuery", statboxDataDGE, "Query")
  
  # Creating list required for volcano server module + launching it
  volcanoDataDGE <- reactive({
    req(search_term() != "")
    list(
      search_term = search_term(),
      dataAll = DGEall,
      cancer_types = c("Melanoma", "Lung", "Prostate", "Glioblastoma")
    )
  })
  
  volcanoServer("volcanoDGE", volcanoDataDGE, "DGE")
  
  output$volcanoText <- renderUI({
    if (search_term() != ""){
      HTML("<i>Note: Outlier points beyond the axes limits have been constrained to improve plot readability.</i>")
    } else {
      ""
    }
  })
  
  #########################
  ###  TABSET 4 : DTE ###
  #########################
  
  # Text to display which gene has been searched and giving info to click on the query sub-tab.
  output$SelectedGeneTextDTE <- renderUI({
    if (search_term() != ""){
      paste0("Currently selected gene: ", search_term(), ". Please click on the \"Query\" tab to get detailed information on it.")
    } else {
      ""
    }
  })
  
  filtersDTE <- filtersBoxServer("filtersDTE")
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTE_data_All <- reactive({
    data <- DTEall

    # Apply log2FoldChange filter based on DEside
    data <- data %>% 
      filter(case_when(
        filtersDTE$DEside() == "both" ~ abs(log2FoldChange) >= filtersDTE$log2fc_threshold(),
        filtersDTE$DEside() == "up" ~ log2FoldChange >= filtersDTE$log2fc_threshold(),
        filtersDTE$DEside() == "down" ~ log2FoldChange <= -filtersDTE$log2fc_threshold()
      ))
    
    # Apply padj filter if not "NONE"
    if (filtersDTE$padj_threshold() != "NONE") {
      data <- data %>% filter(padj <= as.numeric(filtersDTE$padj_threshold()))
    }
    
    # Apply cancer type filter
    if (!is.null(filtersDTE$cancer_types())) {
      data <- data %>% filter(cancer %in% filtersDTE$cancer_types())
    }
    
    return(data)
  })

    # Create a reactive value to store the DTE data,filtering to keep only the searched gene
  filtered_DTE_data_Query <- reactive({
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- DTEall
      data <- data %>% filter(gene_id == search_term() | gene_name == search_term()) %>%
        # Add unique ids to query transcripts to print it on volcano plot
        mutate(nbr = match(transcript_id, unique(transcript_id))) %>%
        select(nbr, everything())
      
      return(data)
    } else {
      NULL
    }
    
  })
  
  output$DTETableAll <- renderDT({
    datatable(filtered_DTE_data_All() %>% 
                mutate(across(where(is.numeric) & !padj, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatStyle(columns = 'cancer',
                  target = 'cell',
                  backgroundColor = styleEqual(
                    c("Melanoma", "Lung", "Prostate", "Glioblastoma"),
                    c("#f5ab05", "#00cc94", "#0084cf", "#eb8dc1")
                  ))
  })

    # Render the DTE table for the query gene (4 rows, one per cancer)
  output$DTETableQuery <- renderDT({
    datatable(filtered_DTE_data_Query() %>% 
                mutate(across(where(is.numeric) & !padj, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatStyle(columns = 'cancer',
                  target = 'cell',
                  backgroundColor = styleEqual(
                    c("Melanoma", "Lung", "Prostate", "Glioblastoma"),
                    c("#f5ab05", "#00cc94", "#0084cf", "#eb8dc1")
                  ))
  })
  
  
  # Creating list required for volcano server module + launching it
  volcanoDataDTE <- reactive({
    req(search_term() != "")
    list(
      search_term = search_term(),
      dataAll = DTEall,
      cancer_types = c("Melanoma", "Lung", "Prostate", "Glioblastoma")
    )
  })
  
  volcanoServer("volcanoDTE", volcanoDataDTE, "DTE")
  
  #########################
  ###  TABSET 5 : DTU ###
  #########################
  
  # Text to display which gene has been searched and giving info to click on the query sub-tab.
  output$SelectedGeneTextDTU <- renderUI({
    if (search_term() != ""){
      paste0("Currently selected gene: ", search_term(), ". Please click on the \"Query\" tab to get detailed information on it.")
    } else {
      ""
    }
  })
  
  ### PREPARE TABLE (FILTERING WITH SEARCH TERM) AND DISPLAY IT IN THE OUTPUT ###
  
  filtersDTU <- filtersBoxServer("filtersDTU")
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTU_data_All <- reactive({
    data <- DTUall
    
    # Apply log2FoldChange filter based on DEside
    data <- data %>% 
      filter(case_when(
        filtersDTU$DEside() == "both" ~ abs(dIF) >= filtersDTU$log2fc_threshold(),
        filtersDTU$DEside() == "up" ~ dIF >= filtersDTU$log2fc_threshold(),
        filtersDTU$DEside() == "down" ~ dIF <= -filtersDTU$log2fc_threshold()
      ))
    
    # Apply padj filter if not "NONE"
    if (filtersDTU$padj_threshold() != "NONE") {
      data <- data %>% filter(isoform_switch_q_value <= as.numeric(filtersDTU$padj_threshold()))
    }
    
    # Apply cancer type filter
    if (!is.null(filtersDTU$cancer_types())) {
      data <- data %>% filter(cancer %in% filtersDTU$cancer_types())
    }
    
    return(data)
  })
  
  # Create a reactive value to store the DTU data, filtering it if a gene as been searched
  filtered_DTU_data_Query <- reactive({
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- DTUall
      data <- data %>% filter(gene_id == search_term() | gene_name == search_term())
      return(data)
    } else {
      NULL
    }
    
  })
  
  
  output$DTUTableAll <- renderDT({
    datatable(filtered_DTU_data_All() %>%
                mutate(across(where(is.numeric) & !isoform_switch_q_value, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  output$DTUTableQuery <- renderDT({
    datatable(filtered_DTU_data_Query()%>%
                mutate(across(where(is.numeric) & !isoform_switch_q_value, ~round(., digits = 3))) %>%
                mutate(across(everything(), ~ifelse(is.na(.), "NA", as.character(.)))),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  ### PREPARE SWITCHPLOT ###
  
  observe({
    
    if (search_term() != "") {
      
      if (search_term() %in% DTUall$gene_id | search_term() %in% DTUall$gene_name){
        
        shinyjs::show("switchplot1")
        shinyjs::show("switchplot2")
        shinyjs::show("switchplot3")
        shinyjs::show("switchplot4")
        
        switchPlotServer(id = "switchplot1",
                         switch_data,
                         "Glioblastoma",
                         search_term())
        switchPlotServer(id = "switchplot2",
                         switch_data,
                         "Lung",
                         search_term())
        switchPlotServer(id = "switchplot3",
                         switch_data,
                         "Melanoma",
                         search_term())
        switchPlotServer(id = "switchplot4",
                         switch_data,
                         "Prostate",
                         search_term())
        
      } else {
        shinyjs::hide("switchplot1")
        shinyjs::hide("switchplot2")
        shinyjs::hide("switchplot3")
        shinyjs::hide("switchplot4")

      }
      
    } else {
      shinyjs::hide("switchplot1")
      shinyjs::hide("switchplot2")
      shinyjs::hide("switchplot3")
      shinyjs::hide("switchplot4")
    }
    
  })
  
  #########################
  ###  INFO BUTTON  ###
  #########################
  
  observeEvent(input$infoButton, {
    showModal(modalDialog(
      title = "Information",
      markdown("
               #### MAIN TAB
               
               Information on genes/transcripts analysed including position and biotype.
               When a gene/transcript is present in the annotation file, the type is annoted `known`, the `new` type corresponds to ones discovered by bambu.
               A robust CAGE (Cap Analysis of Gene Expression) dataset is used to validate transcripts at 5' end, especially for new transcripts/genes (`Y` for CAGE validation).
               The filter category gives another level of confidence for new genes/transcripts and corresponds to the filters included in the pipeline:
               - bambu cut-off: NDR < 0.2
               - transforKmers cut-off: TFK < 0.2
               
               #### COUNT TAB
               
               Count Data corresponds to the number of sequence fragments that have been assigned to each gene or transcript.
               These data are normalized based on TPM (transcripts per million) corresponding to a normalization by gene length followed by a normalization for sequencing depth.
               In the table, expression corresponds to the triplicates median for each condition.
               
               #### DGE/DTE TAB
               
               Differential expression analysis was performed with DESeq2 per cancer using a negative binomial generalized linear models.
               The significance of differential Gene/transcript expression is defined by parameters:
               - p-value adjusted (padj)
               - log 2 fold change corresponding to the change in expression between the 2 conditions (log2FoldChange)
               
               #### DTU TAB
               
               Differential usage analysis was performed with isoformSwitchAnalyzeR. 
               This R package enables statistical identification of isoform switches based on DEXSeq tool.
               A significant isoform switch is defined by parameters:
               - `alpha` corresponding to FDR corrected P-value cut-off (q-value) -> statistical certainty of the difference between 2 conditions
               - `dIF` indicating the minimum absolute change in isoforme usage (dIF) -> reflect the effect size
               "),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
}

# Run the application
shinyApp(ui, server)