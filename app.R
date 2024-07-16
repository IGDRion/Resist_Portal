if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("htmlwidgets")) install.packages("htmlwidgets")
if (!require("plotly")) install.packages("plotly")
if (!require("shiny")) install.packages("shiny")
if (!require("shinyjs")) install.packages("shinyjs")

library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(shiny)
library(shinyjs)


# Ui
ui <- page_navbar(
  
  # Set up shinyjs to use its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  
  title = "Resist_Portal",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "One", 
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
            DTOutput(outputId = "MainTable")
            ),
  nav_panel(title = "Two", p("Second page content.")),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Posit", href = "https://posit.co")),
    nav_item(tags$a("Shiny", href = "https://shiny.posit.co"))
  )
)

# Server
server <- function(input, output, session) {
  
  data <- read.table("DATA/summary_lncRNAresist_annotFull_byGene.txt", header = TRUE)
  
  # Create a reactive value to store the filtered data
  filtered_data <- reactiveVal(data)
  
  # Render the table, with filtered_data if a gene has been searched
  output$MainTable <- renderDT({
    datatable(filtered_data(),
              options = list(ordering = TRUE, pageLength = 10),
              rownames = FALSE)
  })
  
  # Add this observeEvent for the search bar module
  observeEvent(input$submit_btn, {
    search_term <- input$searchBar
    if (search_term != "") {
      filtered <- data %>%
        filter(gene_id == search_term)
      filtered_data(filtered)
    } else {
      filtered_data(data)  # Reset to full data if search is empty
    }
  })
}

# Run the application
shinyApp(ui, server)