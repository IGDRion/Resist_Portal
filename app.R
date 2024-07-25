if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("shiny")) install.packages("shiny")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("stringr")) install.packages("stringr")

library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(reshape2)
library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinyWidgets)
library(stringr)


# Set directory to the src directory of the app so that the path
# is correctly resolved later (otherwise they won't be found)
setwd(getSrcDirectory(function(){})[1])
# Source all modules
source("./module_search.R")
source("./module_boxplot.R")
source("./module_barplot.R")

# Color code of cancers
color_mapping <- list(
  "X501Mel_S" = "#f5ab05",
  "X501Mel_R" = "#8c6100",
  "ADCA72_S" = "#00cc94",
  "ADCA72_R" = "#00664a",
  "PC3_S" = "#0084cf",
  "PC3_R" = "#014c75",
  "U251_S" = "#eb8dc1",
  "U251_R" = "#824a69"
)

# Ui
ui <- page_navbar(
  
  # Set up shinyjs to use its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  
  title = "Resist Portal",
  bg = "#70b684",
  inverse = TRUE,
  nav_panel(title = "Main", 
            p("ECRIRE ICI AURORE ;)"), # Put text here if you need a text on top of the page
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
            # Main summary table
            DTOutput(outputId = "SummaryTable") %>% withSpinner() # withSpinner() display a loading spinner while the dataframe is not on screen
            ),
  
  nav_panel(title = "Count",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
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
                     value_box(value = "",
                               title = "Filters",
                               sliderInput(inputId = "log2fc_threshold",
                                           label = "Log2FC Threshold:",
                                           min = 0, max = 5, value = 0, step = 0.5),
                               sliderTextInput(inputId = "padj_threshold",
                                               label = "padj Threshold: ",
                                               choices = c(0.01, 0.05, "NONE"),
                                               selected = "NONE",
                                               grid = TRUE))
              ), 
              column(width = 6, align = "center",
                     searchBarUI("searchBar", "submit_btn"))
              
            ),
            # DGE table
            DTOutput(outputId = "DGETable") %>% withSpinner()
            ),
  
  nav_panel(title = "DTE",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
            # DTE table
            DTOutput(outputId = "DTETable") %>% withSpinner()
            ),
  
  nav_panel(title = "DTU",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
            # DTE table
            DTOutput(outputId = "DTUTable") %>% withSpinner()
            ),
  
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("GitHub", href = "https://github.com/IGDRion/Resist_Portal")),
    nav_item(tags$a("Shiny", href = "https://shiny.posit.co"))
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
    mutate_if(is.numeric, ~round(., digits = 3))
  DTEall <- DTEall %>%
    mutate_if(is.numeric, ~round(., digits = 3))
  DTUall <- DTUall %>%
    mutate_if(is.numeric, ~round(., digits = 3))
  
  ############################
  ###  SEARCH BUTTON VALUE ###
  ############################
  
  # Create a reactive value to store the gene searched through the researched bar, so it can be used everywhere in the app
  search_term <- reactiveVal("")
  
  # Add this observeEvent so when the search bar module is used, search_term changes
  observeEvent(input$submit_btn, {
    search_term(input$searchBar)
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
        select(-c(gene_id, gene_name))
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
      shinyjs::hide("boxplot1")
      shinyjs::hide("boxplot2")
      shinyjs::hide("barplotTX")
    }
    
  })

  #########################
  ###  TABSET 3 : DGE ###
  #########################

  
  # Create a reactive value to store the DGE data, filtering it if a gene as been searched
  filtered_DGE_data <- reactive({
    if (search_term() != "") {
      DGEall %>% 
        filter(geneID == search_term() | gene_name == search_term()) # If a gene has been searched, find and take only the row with the name of the gene
    } else {
      DGEall # If no gene has been searched, take the full count data
    }
  })

  # Render the filtered DGE table
  output$DGETable <- renderDT({
    datatable(filtered_DGE_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  
  #########################
  ###  TABSET 4 : DTE ###
  #########################
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTE_data <- reactive({
    if (search_term() != "") {
      DTEall %>% 
        filter(gene_id == search_term() | gene_name == search_term()) # If a gene has been searched, find and take only the row with the name of the gene
    } else {
      DTEall # If no gene has been searched, take the full count data
    }
  })
  
  output$DTETable <- renderDT({
    datatable(filtered_DTE_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  
  #########################
  ###  TABSET 5 : DTU ###
  #########################
  
  # Create a reactive value to store the DTE data, filtering it if a gene as been searched
  filtered_DTU_data <- reactive({
    if (search_term() != "") {
      DTUall %>% 
        filter(gene_id == search_term() | gene_name == search_term()) # If a gene has been searched, find and take only the row with the name of the gene
    } else {
      DTUall # If no gene has been searched, take the full count data
    }
  })
  
  output$DTUTable <- renderDT({
    datatable(filtered_DTU_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  
}

# Run the application
shinyApp(ui, server)