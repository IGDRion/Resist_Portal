if (!require("bslib")) install.packages("bslib")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("plotly")) install.packages("plotly")
if (!require("shiny")) install.packages("shiny")
if (!require("shinyjs")) install.packages("shinyjs")

library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(shiny)
library(shinyjs)

# Ui
ui <- page_navbar(
  title = "Resist_Portal",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "One", 
            p("First page content."),
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
  
  data <- read.table("/Users/Victor/Downloads/summary_lncRNAresist_annotFull_byGene.txt", header = TRUE)
  print(data)
  
  # Render the filtered table
  output$MainTable <- renderDT({
    datatable(data, options = list(ordering = TRUE, pageLength = 10), rownames = FALSE)
  })
  
}

# Run the application
shinyApp(ui, server)