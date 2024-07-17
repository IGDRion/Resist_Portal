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
source("./module_boxplot.R")

# Ui
ui <- page_navbar(
  
  # Set up shinyjs to use its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  
  title = "Resist_Portal",
  bg = "#70b684",
  inverse = TRUE,
  nav_panel(title = "Main", 
            p(""), # Put text here if you need a text on top of the page
            
            # Search bar 
            fluidRow(
              column(6, align = "right",
                     textInput(inputId = "searchBar",
                               label = "Enter a gene ID/name to explore",
                               placeholder = "Gene ID/name")
              ),
              column(6, align = "left",
                     actionButton(inputId = "submit_btn", label = "Submit", style = "margin-top: 25px;")
              )
            ),
            
            # Main summary table
            DTOutput(outputId = "SummaryTable")),
  
  nav_panel(title = "Count",
            p(""),
            boxplotUI(id = "boxplot1")),
  
  nav_panel(title = "DGE",
            p("Differential gene expression")),
  
  nav_panel(title = "DTE",
            p("Differential transcript expression")),
  
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
  summary_data <- read.table("./DATA/summary_lncRNAresist_annotFull_byGene.txt", header = TRUE)
  count_data <- read.table("./DATA/countMtx_rlogNorm_lncRNAresist.Rdata", header = TRUE)
  
  
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
        filter(gene_id == search_term())
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
        filter(row.names(.) == search_term()) # If a gene has been searched, find and take only the row with the name of the gene
    } else {
      NULL  # If no gene has been searched, becomes NULL
    }
  })

  # If filtered_count_data() is not NULL (e.g a gene has been searched)
  observe({
    
    if (!(is.null(filtered_count_data()))) {
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
          str_detect(variable, "S[1-3]?$") ~ "Sensitive",
          str_detect(variable, "R[1-3]?$") | str_detect(variable, "TMZ") ~ "Resistant"
        ))
      
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