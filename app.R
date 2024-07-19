if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("shiny")) install.packages("shiny")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("stringr")) install.packages("stringr")

library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(reshape2)
library(shiny)
library(shinyjs)
library(stringr)


# Set directory to the src directory of the app so that the path
# is correctly resolved later (otherwise they won't be found)
setwd(getSrcDirectory(function(){})[1])
# Source all modules
source("./module_search.R")
source("./module_boxplot.R")

# Ui
ui <- page_navbar(
  
  # Set up shinyjs to use its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  
  title = "Resist Portal",
  bg = "#70b684",
  inverse = TRUE,
  nav_panel(title = "Main", 
            p(""), # Put text here if you need a text on top of the page
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
            # Main summary table
            DTOutput(outputId = "SummaryTable")),
  
  nav_panel(title = "Count",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn"),
            tabsetPanel(id = "TabsetCount",
                        tabPanel(title = "Gene",
                                 # Count table
                                 DTOutput(outputId = "CountTable"),
                                 # Boxplot
                                 boxplotUI(id = "boxplot1")),
                        tabPanel(title = "Transcript",
                                 # Count table
                                 DTOutput(outputId = "CountTableTx"))
                        )
            ),
  
  nav_panel(title = "DGE",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn")),
  
  nav_panel(title = "DTE",
            p(""),
            # Search bar 
            searchBarUI("searchBar", "submit_btn")),
  
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
    mutate_if(is.numeric, ~round(., digits = 2))
  
  
  ###########################
  ###  TABSET 1 : SUMMARY ###
  ###########################
  
  # Create a reactive value to store the gene searched through the researched bar, so it can be used everywhere in the app
  search_term <- reactiveVal("")
  
  # Add this observeEvent so when the search bar module is used, search_term changes
  observeEvent(input$submit_btn, {
    search_term(input$searchBar)
  })
  
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
  
  # Create a reactive value to store the count data, NULL by default to not load anything before filtering
  filtered_count_data <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) %>% # If a gene has been searched, find and take only the row with the name of the gene
        group_by(gene_id, gene_name) %>% # Group by gene
        summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop') # Sum all numeric columns
    } else {
      count_data # If no gene has been searched, take the full count data
    }
  })
  
  filtered_count_data_tx <- reactive({
    if (search_term() != "") {
      count_data %>% 
        filter(gene_id == search_term() | gene_name == search_term()) # If a gene has been searched, find and take only the row with the name of the gene
    } else {
      count_data  # If no gene has been searched, take the full count data
    }
  })
  
  # Render the table with gene, with filtered_count_data
  output$CountTable <- renderDT({
    DT::datatable(filtered_count_data(),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  
  # Render the table with transcripts, with filtered_count_data_tx
  output$CountTableTx <- renderDT({
    DT::datatable(filtered_count_data_tx(),
                  options = list(ordering = TRUE, pageLength = 10),
                  rownames = TRUE)
  })
  

  # If search_term() is not NULL (e.g a gene has been searched)
  observe({
    
    if (search_term() != "") {
      # Prepare data for the barplot
      count_barplot_data <- melt(filtered_count_data())
      count_barplot_data <- count_barplot_data %>%
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
        mutate(cancer_condition = paste0(cancer_type, "_", condition)) %>%
        group_by(gene_id, variable, cancer_type, condition, cancer_condition) %>%
        summarise(value = sum(value, na.rm = TRUE), .groups = 'drop')
      
      
      
      boxplotServer(id = "boxplot1",
                    count_barplot_data,
                    search_term())
      
    } else {
      output$boxplot1 <- NULL
    }
  })

}

# Run the application
shinyApp(ui, server)