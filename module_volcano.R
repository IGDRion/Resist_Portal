library(shiny)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(shinycssloaders)

volcanoUI <- function(id){
  ns <- NS(id)
  plotOutput(outputId = ns("plotVolcano")) %>% withSpinner()
}


volcanoServer <- function(id, data) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      filtered_data <- reactive({
        req(data())
        thresholdLOG2FC <- data()$log2fc_threshold
        thresholdPADJ <- data()$padj_threshold
        if (thresholdPADJ == "NONE") {
          thresholdPADJ <- 1
        } else {
          thresholdPADJ <- as.numeric(thresholdPADJ)
        }
        
        filtered <- data()$DGEall %>%
          mutate(
            # Change padj column to -log 
            padj = ifelse((padj == 0 | -log(padj) > 30), exp(-30), padj),
            log2FoldChange = ifelse(abs(log2FoldChange) > 10, 10, log2FoldChange),
            
            # If biotype not lncRNA or protein_coding, it becomes "other"
            gene_biotype = ifelse(is.na(gene_biotype), "other", gene_biotype),
            
            # Add diff column to tag each gene and know in which color they should be on the plot
            diff = case_when(
              geneID == data()$search_term | gene_name == data()$search_term ~ "SEARCHED",
              log2FoldChange >= thresholdLOG2FC & padj <= thresholdPADJ ~ "UPPER",
              log2FoldChange <= -thresholdLOG2FC & padj <= thresholdPADJ ~ "UNDER",
              TRUE ~ "NONE"
            )
          )

          # Keep only upper and under regulated genes to reduce lag when plotting
          filtered <- filtered %>% 
            filter(diff != "NONE") %>%
          # Remove rows with NA in specific columns (columns needed for volcano)
          filter(complete.cases(padj, log2FoldChange, gene_biotype))
        
        return(filtered)
      })

      
      ### Create the volcano plot ###     
      VolcanoPlot <- reactive({
        req(data()$filtered_data)
        
        plot_list <- list()
        
        for(cancer_name in data()$cancer_types) {
          cancer_data <- filtered_data() %>% filter(cancer == cancer_name)
          
          plot <- ggplot(cancer_data, aes(x = log2FoldChange, y = -log10(padj),
                                          fill = factor(diff),
                                          size = gene_annot,
                                          shape = gene_biotype)) +
            # First layer: UNDER and UPPER points
            geom_point(data = subset(cancer_data, diff != "SEARCHED"),
                       aes(stroke = .2),
                       alpha = 0.5) +
            # Second layer: SEARCHED points
            geom_point(data = subset(cancer_data, diff == "SEARCHED"),
                       aes(stroke = .2),
                       size = 10) +
            # add some lines
            geom_vline(xintercept = data()$log2fc_threshold, linetype = "dashed", color = "grey") +
            geom_vline(xintercept = -data()$log2fc_threshold, linetype = "dashed", color = "grey") +
            geom_vline(xintercept = 0, color = "black") +
            # Colors for Under/Overexpressed genes + shapes + x axis limits
            scale_fill_manual(values = c("UNDER" = "#56B4E9", "UPPER" = "#D55E00", "SEARCHED" = "#6ccf41")) +
            scale_shape_manual(values = c("lncRNA" = 21, "protein_coding" = 23, "other" = 22)) +
            scale_x_continuous(limits = c(min(cancer_data$log2FoldChange), 
                                          max(cancer_data$log2FoldChange))) +
            # Text & axis style
            labs(x = "Log2FC", y = "-Log10(p.adj)", title = cancer_name,
                 fill = "Differential Expression", size = "Origin", shape = "Gene Biotype") +
            theme_minimal() +
            theme(legend.position = "none", # Remove individual legend
                  plot.title = element_text(size = 20),
                  axis.title = element_text(size = 15),
                  axis.text = element_text(size = 12),
                  legend.title = element_text(size = 15), # Element text for legend text because of legend extraction just after
                  legend.text = element_text(size = 12))
          
          plot_list[[cancer_name]] <- plot
        }
        
        # Extract legend from one plot to put it in th wrap
        legend_plot <- plot_list[[1]] + 
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
        
        return(final_plot)
      })
      
      ## OLD PLOTLY VERSION TEST ##      

      # VolcanoPlot <- reactive({
      #   req(filtered_data())
      #   
      #   # Separate data for interactive and non-interactive points
      #   interactive_points <- filtered_data() %>% filter(diff == "SEARCHED")
      #   static_points <- filtered_data() %>% filter(diff != "SEARCHED")
      #   
      #   # One per cancer selected by the user
      #   for(cancer in data()$cancer_types){
      #     
      #   }
      #   
      #   # Create the base plot
      #   plot <- plot_ly()
      #   
      #   # Add non-interactive (static) points for UPPER and UNDER
      #   plot <- plot %>% add_trace(
      #     data = static_points,
      #     x = ~log2FoldChange, 
      #     y = ~-log10(padj),
      #     type = 'scatter',
      #     mode = 'markers',
      #     marker = list(
      #       color = ~ifelse(diff == "UPPER", "#D55E00", "#56B4E9"),
      #       size = 8,
      #       opacity = 0.5
      #     ),
      #     showlegend = FALSE,
      #     hoverinfo = 'skip'
      #   )
      #   
      #   # Add layout details
      #   plot <- plot %>% layout(
      #     xaxis = list(title = "Log2FC"),
      #     yaxis = list(title = "-Log10(p.adj)"),
      #     shapes = list(
      #       list(type = "line", x0 = data()$log2fc_threshold, x1 = data()$log2fc_threshold, y0 = 0, y1 = 1, yref = "paper", line = list(color = "grey", dash = "dash")),
      #       list(type = "line", x0 = -data()$log2fc_threshold, x1 = -data()$log2fc_threshold, y0 = 0, y1 = 1, yref = "paper", line = list(color = "grey", dash = "dash")),
      #       list(type = "line", x0 = 0, x1 = 0, y0 = 0, y1 = 1, yref = "paper", line = list(color = "black"))
      #     )
      #   )
      #   
      #   # Add interactive points for SEARCHED
      #   plot <- plot %>% add_trace(
      #     data = interactive_points,
      #     x = ~log2FoldChange, 
      #     y = ~-log10(padj),
      #     type = 'scatter',
      #     mode = 'markers',
      #     marker = list(
      #       color = '#6ccf41',
      #       size = 10,
      #       symbol = 'circle',
      #       line = list(color = 'black', width = 1)  # Add a black border
      #     ),
      #     text = ~paste("GeneID:", geneID, "<br>Gene name:", gene_name, 
      #                   "<br>Log2FoldChange:", log2FoldChange, 
      #                   "<br>p-adjusted:", padj),
      #     hoverinfo = 'text'
      #   )
      #   
      #   return(plot)
      # })

      # Print plot
      
      output$plotVolcano <- renderPlot({
        VolcanoPlot()
      })
      
    }
  )
}
