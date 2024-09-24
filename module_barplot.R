library(ggplot2)
library(shiny)
library(tidyverse)

# UI module example
barplotUI <- function(id) {
  ns <- NS(id)
  
  # UI HERE
  plotOutput(outputId = ns("countBarplotTX"))
}

# Server module example
barplotServer <- function(id, count_barplot_data) {
  moduleServer(
    id = id,
    module = function(input, output, session) {
      
      # Reshape the data from wide to long format
      long_data <- count_barplot_data %>%
        pivot_longer(cols = -c(transcript_id, gene_id, gene_name), names_to = "sample", values_to = "value")
      
      # Summarize the total value for each transcript_id and each group of samples
      summarized_data <- long_data %>%
        mutate(group = case_when(
          str_detect(sample, "^501Mel") ~ "Melanoma",
          str_detect(sample, "^ADCA72") ~ "Lung",
          str_detect(sample, "^PC3") ~ "Prostate",
          str_detect(sample, "^U251") ~ "Glioblastoma"
        )) %>%
        group_by(transcript_id, group) %>%
        summarize(group_sum = sum(value), .groups = 'drop') %>%
        group_by(transcript_id) %>%
        mutate(total_sum = sum(group_sum),
               proportion = group_sum / total_sum) %>%
        ungroup()
      
      # Set the order of the groups explicitly
      summarized_data$group <- factor(summarized_data$group, levels = c("Melanoma", "Lung", "Prostate", "Glioblastoma"))
      
      num_bars <- length(unique(summarized_data$transcript_id))
      
      barplot <- ggplot(summarized_data, aes(x = transcript_id, y = group_sum, fill = group, group = transcript_id)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = c("Melanoma" = '#f5ab05',
                                     "Lung" = '#00cc94',
                                     "Prostate" = '#0084cf',
                                     "Glioblastoma" = '#eb8dc1'),
                          breaks = c("Melanoma", "Lung", "Prostate", "Glioblastoma")) +
        theme_minimal() +
        theme_bw(base_size = 20) +
        coord_flip() +
        scale_x_discrete(limits = rev) + # Inverse order of the bars
        labs(title = "",
             x = "",
             y = "Normalized count (tpm)",
             fill = "Cancer") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      output$countBarplotTX <- renderPlot({
        barplot
      }, height = 300)
    }
  )
}