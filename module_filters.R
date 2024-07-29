library(shiny)

# UI module
filtersBoxUI <- function(id){
  ns <- NS(id)
  
  # Filter box
  value_box(value = "",
            title = "Filters",
            fluidRow(
              column(width = 6, align = "center",
                     sliderInput(inputId = ns("log2fc_threshold"),
                                 label = "Log2FC Threshold:",
                                 min = 0, max = 3, value = 0, step = 0.5)),
              column(width = 6, align = "center",
                     sliderTextInput(inputId = ns("padj_threshold"),
                                     label = "padj Threshold:",
                                     choices = c(0.01, 0.05, "NONE"),
                                     selected = "NONE",
                                     grid = TRUE))
            ),
            fluidRow(
              column(width = 12, align = "center",
                     radioButtons(inputId = ns("DEside"),
                                  label = "DE type:",
                                  choices = c("Down Regulated" = "down",
                                              "Up Regulated" = "up",
                                              "Both" = "both"),
                                  selected = "both"))
            ),
  )
  
}

filtersBoxServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    list(
      DEside = reactive(input$DEside),
      log2fc_threshold = reactive(input$log2fc_threshold),
      padj_threshold = reactive(input$padj_threshold)
    )
  })
}