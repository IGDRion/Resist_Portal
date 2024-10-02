# UI module
switchPlotUI <- function(id){ # ADD ARGUMENTS IF NEEDED
  ns <- NS(id)
 
  # UI HERE
  fluidRow(
    column(12, align = "center",
           plotOutput(outputId = ns("switchPlot"), height = "auto"))
  )

}


# Server module
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
      # rename plot title to fit the search term
      swPlot_tx_name <-  paste0("Isoform switch in ", search_term, " (",cond1," vs ",cond2,")")
      swPlot_tx$labels$title <- swPlot_tx_name
      swPlot_tx <- swPlot_tx + theme(plot.title = element_text(size = 17)) # adjust title size
      
      swPlot_usage <- switchPlotIsoUsage(
        switch_data, 
        gene = search_term,
        condition1 = cond1,
        condition2 = cond2,
        localTheme = theme_bw(base_size = 16) 
      )

      # rename plot title
      swPlot_usage$labels$title <- "Isoform usage"
      
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
