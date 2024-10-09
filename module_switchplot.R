# UI module
switchPlotUI <- function(id){ # ADD ARGUMENTS IF NEEDED
  ns <- NS(id)
 
  # UI HERE
  fluidRow(
    column(12, align = "center",
           uiOutput(outputId = ns("switchPlot")))
  )

}


# Server module
switchPlotServer <- function(id, switch_data, search_term, validity) { # without alpha and IF arguments
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      if (validity == "valid"){
        
        all_outputs <- list()
        
        for(cancerName in c("Prostate","Glioblastoma","Melanoma","Lung")){
          
          
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
          
          # check table if data available
          checkTable <- switch_data$isoformFeatures %>% 
            mutate(maxIFval = max(IF1,IF2)) %>%
            dplyr::filter(condition_1 == cond1, 
                          (gene_id == search_term | gene_name == search_term), 
                          maxIFval > 0.05)
          
          # plot if data are available in checkTable
          if (nrow(checkTable) > 0){
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
            # merge plots if any
            all_outputs[[cancerName]] <- list(type = "plot", content = switchPlots)
            
          } else {
            printText <- paste0("No significant isoform switch for ",search_term," in ",cancerName,".")
            
            # merge no value text if any
            all_outputs[[cancerName]] <- list(type = "text", content = printText)
          }
        } 
        
        
        output$switchPlot <- renderUI({
          plot_outputs <- list()
          text_outputs <- list()
          
          for (output in all_outputs) {
            if (output$type == "plot") {
              plot_outputs[[length(plot_outputs) + 1]] <- output$content
            } else {
              text_outputs[[length(text_outputs) + 1]] <- output$content
            }
          }
          
          tagList(
            if (length(plot_outputs) > 0) {
              renderPlot({
                ggarrange(plotlist = plot_outputs, ncol = 1, nrow = length(plot_outputs))
              }, width = 1920, height = 300 * length(plot_outputs))
            },
            if (length(text_outputs) > 0) {
              div(
                style = "margin-top: 20px; font-size: 16px;",
                lapply(text_outputs, function(text) p(text))
              )
            }
          )
        })
        
      } else if (validity == "not-valid"){
        
        output$switchPlot <- renderUI({
          
          p(paste0("No significant isoform switch for ",search_term," in all cancers"))
          
        })
        
      }
      
    }
  )
}
