# UI module example
boxplotUI <- function(id){
  ns <- NS(id)
  
  # UI HERE
  column(width = 12, align = "center",
         uiOutput(outputId = ns("countBoxplot"))
         )
  
}


plot_count <- function(data, id) {
  
  if (!(sum(data$value >= 0.5))) {
    
    return(list(
      type = "text",
      content = paste0("No significant data for ", id)
    ))
    
  } else {
    
    plot <- ggplot(data = data, aes(x = condition, y = value, fill = cancer_condition)) +
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
    
    return(list(
      type = "plot",
      content = plot
    ))
    
  }
}

# Server part
boxplotServer <- function(id, type, count_boxplot_data, DGEall, search_term) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      # Reorder table for the plot
      # Reorder the cancers
      count_boxplot_data$cancer_condition <- factor(count_boxplot_data$cancer_condition, 
                                                    levels = c("Melanoma_Sensitive", "Melanoma_Resistant",
                                                               "Lung_Sensitive", "Lung_Resistant",
                                                               "Prostate_Sensitive", "Prostate_Resistant",
                                                               "Glioblastoma_Sensitive", "Glioblastoma_Resistant"))
      # Reorder the conditions
      count_boxplot_data$condition <- factor(count_boxplot_data$condition, 
                                             levels = c("Sensitive", "Resistant"))
      
      if (type == "gene") { # For the count tab at gene level, make one boxplot for the 4 cancers
        boxplot_result <- plot_count(count_boxplot_data, search_term)
        
        output$countBoxplot <- renderUI({
          if (boxplot_result$type == "plot") {
            renderPlot({
              boxplot_result$content
            }, width = 1000, height = 500)
          } else {
            div(
              style = "margin-top: 20px; font-size: 16px;",
              p(boxplot_result$content)
            )
          }
        })
        
      } else if (type == "tx") { # For the count tab at transcript level, make as many boxplot for the 4 cancers as there are transcript for the query gene
        boxplot_results <- count_boxplot_data %>%
          group_by(transcript_id) %>%
          group_split() %>%
          map(~ plot_count(.x, .x$transcript_id[1]))
        
        output$countBoxplot <- renderUI({
          plot_outputs <- list()
          text_outputs <- list()
          
          for (result in boxplot_results) {
            if (result$type == "plot") {
              plot_outputs[[length(plot_outputs) + 1]] <- result$content
            } else {
              text_outputs[[length(text_outputs) + 1]] <- result$content
            }
          }
          
          tagList(
            if (length(plot_outputs) > 0) {
              renderPlot({
                wrap_plots(plot_outputs, ncol = 1)
              }, width = 1000, height = 500 * length(plot_outputs))
            },
            if (length(text_outputs) > 0) {
              div(
                style = "margin-top: 20px; font-size: 16px;",
                lapply(text_outputs, function(text) p(text))
              )
            }
          )
        })
        
      } else if (type == "DGE") { # For the DGE tab, make a boxplot per cancer, but separated with patchwork as for volcano plots, to align them under the volcano plots
        
        plot_list <- list()
        
        for (cancer_name in c("Melanoma", "Lung", "Prostate", "Glioblastoma")){
          log2FC <- DGEall %>%
            filter(gene_id == search_term | gene_name == search_term) %>%
            filter(cancer == cancer_name) %>%
            select(log2FoldChange)
            
          padj <- DGEall %>%
            filter(gene_id == search_term | gene_name == search_term) %>%
            filter(cancer == cancer_name) %>%
            select(padj)
          
          data <- count_boxplot_data %>%
            filter(cancer_condition == paste0(cancer_name, "_Sensitive") | cancer_condition == paste0(cancer_name, "_Resistant"))
          
          plot <- ggplot(data = data, aes(x = condition, y = value, fill = cancer_condition)) +
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
            ggtitle(paste0("Log2FC ", round(log2FC, digits = 3), " / padj ", signif(padj, digits = 3))) +
            theme_bw(base_size = 20) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1))
          
            plot_list[[cancer_name]] <- plot
        }
        
        # Make a plot but just to have complete legend
        plot_for_legend <- ggplot(data = count_boxplot_data, aes(x = condition, y = value, fill = cancer_condition)) +
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
          ggtitle(paste0("Log2FC ", round(log2FC, digits = 3), " / padj ", signif(padj, digits = 3))) +
          theme_bw(base_size = 20) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1))
        
        
        
        # Extract legend from one plot to put it in the wrap
        legend_plot <- plot_for_legend + 
          theme(legend.position = "bottom") +
          guides(fill = guide_legend(override.aes = list(size = 6, shape = 21), nrow = 1),
                 size = guide_legend(nrow = 1),
                 shape = guide_legend(override.aes = list(size = 6), nrow = 1))
        legend <- get_legend(legend_plot)
        
        # Combine plots and legend using patchwork
        combined_plot <- wrap_plots(plot_list, ncol = 4) +
          plot_layout(guides = "collect") &
          theme(legend.position = "none")
          
        final_plot <- combined_plot / legend +
          plot_layout(heights = c(4, 1))  # Adjust the ratio as needed
        
        output$countBoxplot <- renderUI({
          renderPlot({final_plot})
        })
        
      }
    }
  )
}