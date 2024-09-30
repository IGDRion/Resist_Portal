library(shiny)
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(ggpubr)

# UI module example
switchPlotUI <- function(id){ # ADD ARGUMENTS IF NEEDED
  ns <- NS(id)
 
  # UI HERE
  fluidRow(
    column(12, align = "center",
           plotOutput(outputId = ns("switchPlot"), height = "auto"))
  )

}


# Server module example
switchPlotServer <- function(id, switch_data, cancerName, search_term) { # without alpha and IF arguments
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      cond1 <- case_when(
        cancerName == "Prostate" ~ "prostate_cancer_sensitive",
        cancerName == "Glioblastoma" ~ "glioblastoma_sensitive",
        cancerName == "Melanoma" ~ "melanoma_sensitive",
        cancerName == "Lung"~ "lung_cancer_sensitive"
      )
      cond2 <- case_when(
        cancerName == "Prostate" ~ "prostate_cancer_resistant",
        cancerName == "Glioblastoma" ~ "glioblastoma_resistant",
        cancerName == "Melanoma" ~ "melanoma_resistant",
        cancerName == "Lung"~ "lung_cancer_resistant"
      )
      swPlot_tx <- switchPlotTranscript(
        switch_data,
        gene = search_term,
        condition1 = cond1,
        condition2 = cond2,
        plotTopology = FALSE,
        localTheme = theme_bw(base_size = 16) # making text sightly larger for vignette
      )
      swPlot_usage <- switchPlotIsoUsage(
        switch_data, 
        gene = search_term,
        condition1 = cond1,
        condition2 = cond2,
        localTheme = theme_bw(base_size = 16) 
      )
      switchPlots <- ggarrange(swPlot_tx, swPlot_usage,
                               ncol = 2, nrow = 1)
      
      output$switchPlot <- renderPlot({
        switchPlots
      }, width = 1920, height = function() {
        300  # Adjust height based on the number of plots
      })
      
    }
  )
}
