library(shiny)
library(ggplot2)

# UI module example
boxplotUI <- function(id){
  ns <- NS(id)
  
  # UI HERE
  column(width = 12, align = "center",
         plotOutput(outputId = ns("countBarplot")))
  
}


# Server module example
boxplotServer <- function(id, count_barplot_data, search_term) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      # LOGIC HERE
      count_boxplot <- ggplot(data = count_barplot_data, aes(x = cancer_type, y = value, fill = cancer_type)) +
        geom_violin(trim = FALSE, alpha = 0.5) +
        geom_boxplot(width = 0.1, fill = "grey") +
        geom_jitter(aes(shape = condition), size = 3) +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = c("Melanoma" = '#E69F00' ,"Lung" = '#009E73' ,
                                     "Prostate" = '#0072B2' ,"Glioblastoma" = '#CC79A7'), 
                          breaks=c("Melanoma", "Lung", "Prostate","Glioblastoma")) + # Change the colors of the fill aestethics
        labs(fill = "Cancer type", shape = "Condition") + # Rename the legend titles
        guides(fill = guide_legend(override.aes = list(shape = NA))) + # This removes the point from the fill legend
        ylab("Normalized count") +
        xlab(" ") +
        ggtitle(paste0("Count data of ", search_term, " across cancers")) +
        theme_bw(base_size = 20)
      
      
      output$countBarplot <- renderPlot({
        count_boxplot
      }, width = 1000, height = 400)
      
      
    }
  )
}