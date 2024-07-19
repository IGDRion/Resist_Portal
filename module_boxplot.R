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
      
      print(count_barplot_data)
      
      # Reorder table for the plot
      count_barplot_data$cancer_condition <- factor(count_barplot_data$cancer_condition, 
                                                    levels = c("Melanoma_Sensitive", "Melanoma_Resistant",
                                                               "Lung_Sensitive", "Lung_Resistant",
                                                               "Prostate_Sensitive", "Prostate_Resistant",
                                                               "Glioblastoma_Sensitive", "Glioblastoma_Resistant"))
      
      print(count_barplot_data)
      
      # LOGIC HERE
      count_boxplot <- ggplot(data = count_barplot_data, aes(x = cancer_condition, y = value, fill = cancer_condition)) +
        geom_violin(trim = FALSE, alpha = 0.5) +
        geom_boxplot(width = 0.1, fill = "grey") +
        geom_jitter(aes(shape = condition), size = 3) +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = c("Melanoma_Sensitive" = '#f5ab05' , "Melanoma_Resistant" = '#8c6100',
                                     "Lung_Sensitive" = '#00cc94' , "Lung_Resistant" = '#00664a',
                                     "Prostate_Sensitive" = '#0084cf' , "Prostate_Resistant" = '#014c75',
                                     "Glioblastoma_Sensitive" = '#eb8dc1', "Glioblastoma_Resistant" = '#824a69'), 
                          breaks=c("Melanoma_Sensitive" , "Melanoma_Resistant",
                                   "Lung_Sensitive", "Lung_Resistant",
                                   "Prostate_Sensitive", "Prostate_Resistant",
                                   "Glioblastoma_Sensitive", "Glioblastoma_Resistant")) + # Change the colors of the fill aesthetics + order in legend
        labs(fill = "Cancer type", shape = "Condition") + # Rename the legend titles
        guides(fill = guide_legend(override.aes = list(shape = NA))) + # This removes the point from the fill legend
        ylab("Normalized count") +
        xlab(" ") +
        ggtitle(paste0("Count data of ", search_term, " across cancers")) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      
      output$countBarplot <- renderPlot({
        count_boxplot
      }, width = 1000, height = 500)
      
      
    }
  )
}