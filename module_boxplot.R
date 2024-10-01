# UI module example
boxplotUI <- function(id){
  ns <- NS(id)
  
  # UI HERE
  column(width = 12, align = "center",
         plotOutput(outputId = ns("countBoxplot"), height = "auto"))
  
}


# Server module example
boxplotServer <- function(id, type, count_boxplot_data, search_term) {
  moduleServer(
    id = id,
    module = function(input, output, session) {

      # Reorder table for the plot
      count_boxplot_data$cancer_condition <- factor(count_boxplot_data$cancer_condition, 
                                                    levels = c("Melanoma_Sensitive", "Melanoma_Resistant",
                                                               "Lung_Sensitive", "Lung_Resistant",
                                                               "Prostate_Sensitive", "Prostate_Resistant",
                                                               "Glioblastoma_Sensitive", "Glioblastoma_Resistant"))
      
      plot_count <- function(data, id) {
        
        ggplot(data = data, aes(x = condition, y = value, fill = cancer_condition)) +
          #geom_violin(trim = FALSE, alpha = 0.5) +
          geom_boxplot(width = 0.2) +
          geom_jitter(aes(shape = condition), size = 3) +
          scale_y_continuous(labels = scales::comma) +
          scale_fill_manual(values = c("Melanoma_Sensitive" = '#f5ab05' , "Melanoma_Resistant" = '#8c6100',
                                       "Lung_Sensitive" = '#00cc94' , "Lung_Resistant" = '#00664a',
                                       "Prostate_Sensitive" = '#0084cf' , "Prostate_Resistant" = '#014c75',
                                       "Glioblastoma_Sensitive" = '#eb8dc1', "Glioblastoma_Resistant" = '#824a69'), 
                            breaks=c("Melanoma_Sensitive" , "Melanoma_Resistant",
                                     "Lung_Sensitive", "Lung_Resistant",
                                     "Prostate_Sensitive", "Prostate_Resistant",
                                     "Glioblastoma_Sensitive", "Glioblastoma_Resistant")) + 
          labs(fill = "Cancer type", shape = "Condition") + 
          guides(fill = guide_legend(override.aes = list(shape = NA))) + 
          ylab("Normalized count (tpm)") +
          xlab(" ") +
          ggtitle(paste0("Count data of ", id, " across cancers")) +
          theme_bw(base_size = 20) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          facet_wrap(~ cancer_type, scales = "fixed", ncol = 4)
        
      }
      
      # depending of the type variable when the module was called
      # Do one plot for the whole gene
      if (type == "gene") {
        
        boxplot <- plot_count(count_boxplot_data, search_term)
        
        output$countBoxplot <- renderPlot({
          boxplot
        }, width = 1000, height = 500)
      
      # Or one plot per transcript, on top of each other    
      } else if (type == "tx"){
        
        boxplots <- count_boxplot_data %>%
          group_by(transcript_id) %>%
          group_map(~ plot_count(.x, .y$transcript_id[1]))
        
        stacked_plot <- purrr::reduce(boxplots, `+`) + plot_layout(ncol = 1) # Stack all boxplot on top of each other
        
        output$countBoxplot <- renderPlot({
          stacked_plot
        }, width = 1000, height = function() {
          500 * length(boxplots)  # Adjust height based on the number of plots
        })
      }
      
    }
  )
}