if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("plotly")) install.packages("plotly")
if (!require("remotes")) install.packages("remotes")
if (!require("reshape2")) install.packages("reshape2")
if (!require("shiny")) install.packages("shiny")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("stringr")) install.packages("stringr")
if (!require("summaryBox")) remotes::install_github("deepanshu88/summaryBox")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("IsoformSwitchAnalyzeR")) BiocManager::install("IsoformSwitchAnalyzeR")

library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(reshape2)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(stringr)
library(summaryBox)
library(IsoformSwitchAnalyzeR)

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
  nav_panel(title = "Main", 
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
            searchBarUI("searchBar", "submit_btn", "reset_btn"),
            # Main summary table
            DTOutput(outputId = "SummaryTable") %>% withSpinner(), # withSpinner() display a loading spinner while the dataframe is not on screen
            # UCSC link
            uiOutput("UCSClink")
            ),
  
  nav_panel(title = "Count",
            p("Expressions are normalized based on TPM (transcripts per million)."),
            # Search bar 
            searchBarUI("searchBar", "submit_btn", "reset_btn"),
            tabsetPanel(id = "TabsetCount",
                        tabPanel(title = "Gene",
                                 # Count table
                                 DTOutput(outputId = "CountTable") %>% withSpinner(),
                                 # Boxplot
                                 boxplotUI(id = "boxplot1")),
                        tabPanel(title = "Transcript",
                                 # Count table + Barplot
                                  DTOutput(outputId = "CountTableTx") %>% withSpinner(),
                                  barplotUI(id = "barplotTX") %>% withSpinner(),
                                 # Boxplot
                                  boxplotUI(id = "boxplot2"))
            )
  ),
  nav_panel(title = "DGE",
            p(""),
            # Search bar + Filter box
            fluidRow(
              column(width = 6, align = "center",
                     filtersBoxUI("filtersDGE", Dtype = "DE")
              ), 
              column(width = 6, align = "center",
                     searchBarUI("searchBar", "submit_btn", "reset_btn"),
                     # Box test
                     statboxUI("statboxDGE")
                     )
            ),
            
            # DGE table
            DTOutput(outputId = "DGETable") %>% withSpinner(),
            
            # Volcano plot
            volcanoUI("volcanoDGE"),
            
            # Text to explain that volcano plot points are stopped when they are above a certain limit
            uiOutput("volcanoText")
            ),
  
  nav_panel(title = "DTE",
            p(""),
            # Write which gene is currently selected at the top of the page
            uiOutput("SelectedGeneText"),
            # Search bar + Filter box
            fluidRow(
              column(width = 6, align = "center",
                     filtersBoxUI("filtersDTE", Dtype = "DE")
              ), 
              column(width = 6, align = "center",
                     searchBarUI("searchBar", "submit_btn", "reset_btn"))
              
            ),
            # DTE table
            DTOutput(outputId = "DTETable") %>% withSpinner()
            ),
  
  nav_panel(title = "DTU",
            p(""),
            # Search bar + Filter box
            fluidRow(
              column(width = 6, align = "center",
                     filtersBoxUI("filtersDTU", Dtype = "DU")
              ), 
              column(width = 6, align = "center",
                     searchBarUI("searchBar", "submit_btn", "reset_btn"))
              
            ),
            # DTU table
            DTOutput(outputId = "DTUTable") %>% withSpinner(),
            # switch plot
            switchPlotUI(id = "switchplot1") %>% withSpinner(),
            switchPlotUI(id = "switchplot2") %>% withSpinner(),
            switchPlotUI(id = "switchplot3") %>% withSpinner(),
            switchPlotUI(id = "switchplot4")
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
  count_data <- readRDS("./DATA/resist_transcript_expression.rds") %>%
    mutate_if(is.numeric, ~round(., digits = 3))
  load("./DATA//Differential_analysis.RData") #Directly load into the environment DGEall, DTEall, and DTUall
  DGEall <- DGEall %>%
    mutate(across(where(is.numeric) & !padj, ~round(., digits = 3)))
  DTEall <- DTEall %>%
    mutate(across(where(is.numeric) & !padj, ~round(., digits = 3)))
  DTUall <- DTUall %>%
    mutate_if(is.numeric, ~round(., digits = 3))
  switch_data <- readRDS("./DATA/All_switchlist_DEXSeq.Rds")
  
  ############################
  ###  SEARCH BUTTON VALUE ###
  ############################
  
  # Create a reactive value to store the gene searched through the researched bar, so it can be used everywhere in the app
  search_term <- reactiveVal("")
  
  # Add this observeEvent so when the search bar module is used, search_term changes
  observeEvent(input$submit_btn, {
    search_term(input$searchBar)
  })
  
  # Add this observeEvent so when the reset button of the search bar module is used, search_term becomes NULL and filtered_data becomes the full dataset again
  observeEvent(input$reset_btn, {
    search_term("")
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
      a(paste0("View ", term, " on UCSC"), href=paste0("https://genome.ucsc.edu/cgi-bin/hgSearch?search=",term,"&db=hg38"), target = "_blank")
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
  
  ### PREPARE DATAFRAME (FILTERING WITH SEARCH TERM) AND DISPLAY THEM IN THE OUTPUT ###
  
  # Create a reactive value to store the count data, filtering it if a gene as been searched
  filtered_count_data <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) %>% # If a gene has been searched, find and take only the row with the name of the gene
        group_by(gene_id, gene_name) %>% # Group by gene
        summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop') %>% # Sum all numeric columns
        mutate_if(is.numeric, ~round(., digits = 3)) # Round it again because the sum made it with more than 3 digits again
    } else {
      count_data # If no gene has been searched, take the full count data
    }
  })
  
  filtered_count_data_tx <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) %>% # If a gene has been searched, find and take only the row with the name of the gene
        as_tibble() # to match filtered_count_data() that went tibble when using group_by and summarize
    } else {
      count_data  # If no gene has been searched, take the full count data
    }
  })
  
  # Render the table with gene, with filtered_count_data
  output$CountTable <- renderDT({
    DT::datatable(filtered_count_data() %>%
                    group_by(gene_id, gene_name) %>%
                    summarise(
                      X501Mel_S = median(c(X501Mel_1_S, X501Mel_2_S, X501Mel_3_S), na.rm = TRUE),
                      X501Mel_R = median(c(X501Mel_1_R, X501Mel_2_R, X501Mel_3_R), na.rm = TRUE),
                      ADCA72_S = median(c(ADCA72_1_S, ADCA72_2_S, ADCA72_3_S), na.rm = TRUE),
                      ADCA72_R = median(c(ADCA72_1_R, ADCA72_2_R, ADCA72_3_R), na.rm = TRUE),
                      PC3_S = median(c(PC3_1_S, PC3_2_S, PC3_3_S), na.rm = TRUE),
                      PC3_R = median(c(PC3_1_R, PC3_2_R, PC3_3_R), na.rm = TRUE),
                      U251_S = median(c(U251_1_S, U251_2_S, U251_3_S), na.rm = TRUE),
                      U251_R = median(c(U251_1_R, U251_2_R, U251_3_R), na.rm = TRUE),
                      .groups = 'drop'
                    ),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  
  # Render the table with transcripts, with filtered_count_data_tx
  output$CountTableTx <- renderDT({
    DT::datatable(filtered_count_data_tx() %>%
                    group_by(transcript_id, gene_id, gene_name) %>%
                    summarise(
                      X501Mel_S = median(c(X501Mel_1_S, X501Mel_2_S, X501Mel_3_S), na.rm = TRUE),
                      X501Mel_R = median(c(X501Mel_1_R, X501Mel_2_R, X501Mel_3_R), na.rm = TRUE),
                      ADCA72_S = median(c(ADCA72_1_S, ADCA72_2_S, ADCA72_3_S), na.rm = TRUE),
                      ADCA72_R = median(c(ADCA72_1_R, ADCA72_2_R, ADCA72_3_R), na.rm = TRUE),
                      PC3_S = median(c(PC3_1_S, PC3_2_S, PC3_3_S), na.rm = TRUE),
                      PC3_R = median(c(PC3_1_R, PC3_2_R, PC3_3_R), na.rm = TRUE),
                      U251_S = median(c(U251_1_S, U251_2_S, U251_3_S), na.rm = TRUE),
                      U251_R = median(c(U251_1_R, U251_2_R, U251_3_R), na.rm = TRUE),
                      .groups = 'drop'
                    ),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  
  
  ### PREPARE DATAFRAME USED IN BOXPLOTS (FROM FILTERED DATAFRAME ABOVE) AND DISPLAY THEM IN THE OUTPUT ###
  
  prepare_table <- function(data, type) {
    
    if (type == "gene") {
      # Prepare data for the barplot with gene_id
      data <- melt(data)
    } else if (type == "tx") {
      # Prepare data for the barplot with transcript_id
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
        startsWith(variable, "X501Mel") ~ "Melanoma"
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
    
    return(data)
  }
  
  
  # If search_term() is not NULL (e.g a gene has been searched)
  observe({

    if (search_term() != "") {
      
      shinyjs::show("boxplot1")
      shinyjs::show("boxplot2")
      shinyjs::show("barplotTX")
      
      count_barplot_data <- prepare_table(filtered_count_data(), "gene")
      
      count_barplot_data_tx <- prepare_table(filtered_count_data_tx(), "tx")
      
      boxplotServer(id = "boxplot1",
                    type = "gene",
                    count_barplot_data,
                    search_term())
      
      boxplotServer(id = "boxplot2",
                    type = "tx",
                    count_barplot_data_tx,
                    search_term())
      
      barplotServer(id = "barplotTX",
                    filtered_count_data_tx())
      
    } else {
      output$boxplot1 <- renderUI({ NULL })
      output$boxplot2 <- renderUI({ NULL })
      output$barplotTX <- renderUI({ NULL })
    }
    
  })

  #########################
  ###  TABSET 3 : DGE ###
  #########################

  filtersDGE <- filtersBoxServer("filtersDGE")
  
  # Create a reactive value to store the DGE data, filtering it if a gene as been searched, and filtering with the sliders log2fc / padj
  filtered_DGE_data <- reactive({
    data <- DGEall
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- data %>% filter(geneID == search_term() | gene_name == search_term())
    }
    
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

  # Render the filtered DGE table
  output$DGETable <- renderDT({
    datatable(filtered_DGE_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatStyle(columns = 'cancer',
                  target = 'row',
                  backgroundColor = styleEqual(
                    c("Melanoma", "Lung", "Prostate", "Glioblastoma"),
                    c("#f5ab05", "#00cc94", "#0084cf", "#eb8dc1")
                  )) # FormatStyle not working
  })
  
  # Creating list required for statbox server module + launching it
  statboxDataDGE <- reactive({
    req(search_term() != "")
    list(
      search_term = search_term(),
      DGEall = DGEall,
      cancer_types = filtersDGE$cancer_types()
    )
  })
  
  statboxServer("statboxDGE", statboxDataDGE)
  
  # Creating list required for volcano server module + launching it
  volcanoData <- reactive({
    req(search_term() != "")
    list(
      search_term = search_term(),
      DGEall = DGEall,
      filtered_data = filtered_DGE_data(),
      log2fc_threshold = filtersDGE$log2fc_threshold(),
      padj_threshold = filtersDGE$padj_threshold(),
      cancer_types = filtersDGE$cancer_types()
    )
  })
  
  volcanoServer("volcanoDGE", volcanoData)
  
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
  
  output$SelectedGeneText <- renderUI({
    if (search_term() != ""){
      paste0("Currently selected gene: ", search_term())
    } else {
      ""
    }
  })
  
  
  filtersDTE <- filtersBoxServer("filtersDTE")
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTE_data <- reactive({
    data <- DTEall
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- data %>% filter(gene_id == search_term() | gene_name == search_term())
    }
    
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
  
  output$DTETable <- renderDT({
    datatable(filtered_DTE_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  
  #########################
  ###  TABSET 5 : DTU ###
  #########################
  
  ### PREPARE TABLE (FILTERING WITH SEARCH TERM) AND DISPLAY IT IN THE OUTPUT ###
  
  filtersDTU <- filtersBoxServer("filtersDTU")
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTU_data <- reactive({
    data <- DTUall
    
    # Apply gene search filter if a search term is provided
    if (search_term() != "") {
      data <- data %>% filter(gene_id == search_term() | gene_name == search_term())
    }
    
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
  
  output$DTUTable <- renderDT({
    datatable(filtered_DTU_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  ### PREPARE SWITCHPLOT ###
  
  observe({
    
    if (search_term() != "") {
      
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