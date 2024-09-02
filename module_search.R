library(shiny)

# UI module
searchBarUI <- function(id_bar, id_button, id_reset){

  # Search bar 
  fluidRow(
    column(6, align = "right",
           textInput(inputId = id_bar,
                     label = "Enter a gene ID/name to explore",
                     placeholder = "Gene ID/name")
    ),
    column(6, align = "left",
           actionButton(inputId = id_button, 
                        label = "Submit", 
                        style = "margin-top: 25px;"),
           actionButton(inputId = id_reset, 
                        label = "Reset", 
                        style = "margin-top: 25px;")
    )
  )
  
}
