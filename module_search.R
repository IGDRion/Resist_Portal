library(shiny)

# UI module
searchBarUI <- function(id_bar, id_button, id_reset){
  tagList(
    tags$style(HTML("
      .selectize-control.single .selectize-input::after {
        content: none !important;
      }
    ")),
    fluidRow(
      column(6, align = "right",
             selectizeInput(inputId = id_bar,
                            label = "Enter a gene ID/name to explore",
                            choices = NULL,
                            selected = "",  # Set to empty string
                            multiple = FALSE,
                            options = list(
                              placeholder = 'e.g ENSG00000141510 / TP53',
                              create = FALSE,
                              maxOptions = 10,
                              loadThrottle = 100
                            )
             )
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
  )
}

# Server function
searchBarServer <- function(id_bar, autocomplete_list) {
  updateSelectizeInput(inputId = id_bar, 
                       server = TRUE,
                       choices = autocomplete_list,
                       selected = "",  # Set to empty string when updating the search bars
                       session = getDefaultReactiveDomain(),
                       options = list(
                         placeholder = 'e.g ENSG00000141510 / TP53',
                         create = FALSE,
                         maxOptions = 10,
                         loadThrottle = 100
                       )
  )
}