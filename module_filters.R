library(shiny)

# UI module
filtersBoxUI <- function(id, Dtype){
  ns <- NS(id)
  
  # Filter box
  if (Dtype == "DE"){
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
              br(),
              fluidRow(
                column(width = 6, align = "center",
                       radioButtons(inputId = ns("DEside"),
                                    label = "DE type:",
                                    choices = c("Down Regulated" = "down",
                                                "Up Regulated  " = "up",
                                                "Both          " = "both"),
                                    selected = "both")),
                column(width = 6, align = "center",
                       checkboxGroupInput(inputId = ns("cancer_types"),
                                    label = "Cancer:",
                                    choices = c("Melanoma",
                                                "Lung",
                                                "Prostate",
                                                "Glioblastoma"),
                                    selected = c("Melanoma",
                                                "Lung",
                                                "Prostate",
                                                "Glioblastoma")
                                    ))
              )
    )
  } else if (Dtype == "DU"){
    value_box(value = "",
              title = "Filters",
              fluidRow(
                column(width = 6, align = "center",
                       sliderInput(inputId = ns("log2fc_threshold"),
                                   label = "dIF Threshold:",
                                   min = 0, max = 1, value = 0, step = 0.1)),
                column(width = 6, align = "center",
                       sliderTextInput(inputId = ns("padj_threshold"),
                                       label = "qvalue Threshold:",
                                       choices = c(0.01, 0.05, "NONE"),
                                       selected = "NONE",
                                       grid = TRUE))
              ),
              br(),
              fluidRow(
                column(width = 6, align = "center",
                       radioButtons(inputId = ns("DEside"),
                                    label = "DE type:",
                                    choices = c("Down Regulated" = "down",
                                                "Up Regulated  " = "up",
                                                "Both          " = "both"),
                                    selected = "both")),
                column(width = 6, align = "center",
                       checkboxGroupInput(inputId = ns("cancer_types"),
                                          label = "Cancer:",
                                          choices = c("Melanoma",
                                                      "Lung",
                                                      "Prostate",
                                                      "Glioblastoma"),
                                          selected = c("Melanoma",
                                                       "Lung",
                                                       "Prostate",
                                                       "Glioblastoma")
                       ))
              )
    )
  }
  
  
}

filtersBoxServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    list(
      DEside = reactive(input$DEside),
      log2fc_threshold = reactive(input$log2fc_threshold),
      padj_threshold = reactive(input$padj_threshold),
      cancer_types = reactive(input$cancer_types)
    )
  })
}