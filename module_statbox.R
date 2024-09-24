statboxUI <- function(id){
  ns <- NS(id)
  
  tagList(
    fluidRow(
      summaryBox(
        title = HTML('<span style="color: #f5ab05;">Melanoma</span>'),  # Yellow
        HTML(paste(tags$span(style = "font-size: 0.8em;", textOutput(ns("melanoma_1st_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("melanoma_2nd_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("melanoma_3rd_line"))))), 
        width = 6, icon = "fas fa-dna", style = "success", border = "bottom"
      ),
      summaryBox(
        title = HTML('<span style="color: #00cc94;">Lung</span>'),  # Green
        HTML(paste(tags$span(style = "font-size: 0.8em;", textOutput(ns("lung_1st_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("lung_2nd_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("lung_3rd_line"))))), 
        width = 6, icon = "fas fa-dna", style = "success", border = "bottom"
      )
    ),
    br(),
    fluidRow(
      summaryBox(
        title = HTML('<span style="color: #0084cf;">Prostate</span>'),  # Blue
        HTML(paste(tags$span(style = "font-size: 0.8em;", textOutput(ns("prostate_1st_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("prostate_2nd_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("prostate_3rd_line"))))), 
        width = 6, icon = "fas fa-dna", style = "success", border = "bottom"
      ),
      summaryBox(
        title = HTML('<span style="color: #eb8dc1;">Glioblastoma</span>'),  # Pink
        HTML(paste(tags$span(style = "font-size: 0.8em;", textOutput(ns("glioblastoma_1st_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("glioblastoma_2nd_line"))),
                   #"<br>",
                   tags$span(style = "font-size: 0.8em;", textOutput(ns("glioblastoma_3rd_line"))))), 
        width = 6, icon = "fas fa-dna", style = "success", border = "bottom"
      )
    )
  )
}

statboxServer <- function(id, data, target) {
  moduleServer(
    id,
    function(input, output, session) {
      
      if (target == "All") {
        
        # the mean of log2FC of all gene for each cancer (separating over and under expressed genes)
        cancer_means <- reactive({
          req(data()$DGEall)
          data()$DGEall %>%
            group_by(cancer) %>%
            summarise(
              mean_log2FC_up = mean(log2FoldChange[log2FoldChange >= 0], na.rm = TRUE),
              mean_log2FC_down = mean(log2FoldChange[log2FoldChange < 0], na.rm = TRUE)
            )
        })

        # Melanoma : upper and under mean
        output$melanoma_1st_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Melanoma") %>% 
            dplyr::pull(mean_log2FC_up)
          
          sprintf("Mean upper Log2FC: %.3f", mean_value)
        })

        output$melanoma_2nd_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Melanoma") %>% 
            dplyr::pull(mean_log2FC_down)
          
          sprintf("Mean under Log2FC: %.3f", mean_value)
        })
        
        # Lung : upper and under mean
        output$lung_1st_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Lung") %>% 
            dplyr::pull(mean_log2FC_up)
          
          sprintf("Mean upper Log2FC: %.3f", mean_value)
        })
        
        output$lung_2nd_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Lung") %>% 
            dplyr::pull(mean_log2FC_down)
          
          sprintf("Mean under Log2FC: %.3f", mean_value)
        })
        
        # Prostate : upper and under mean
        output$prostate_1st_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Prostate") %>% 
            dplyr::pull(mean_log2FC_up)
          
          sprintf("Mean upper Log2FC: %.3f", mean_value)
        })
        
        output$prostate_2nd_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Prostate") %>% 
            dplyr::pull(mean_log2FC_down)
          
          sprintf("Mean under Log2FC: %.3f", mean_value)
        })
        
        # Glioblastoma : upper and under mean
        output$glioblastoma_1st_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Glioblastoma") %>% 
            dplyr::pull(mean_log2FC_up)
          
          sprintf("Mean upper Log2FC: %.3f", mean_value)
        })
        
        output$glioblastoma_2nd_line <- renderText({
          mean_value <- cancer_means() %>% 
            dplyr::filter(cancer == "Glioblastoma") %>% 
            dplyr::pull(mean_log2FC_down)
          
          sprintf("Mean under Log2FC: %.3f", mean_value)
        })
        
      }
      
      if (target == "Query"){
        
        # The log2FC of the searched gene
        log2fc_values <- reactive({
          req(data()$DGEall, data()$search_term)
          data()$DGEall %>%
            filter(geneID == data()$search_term | gene_name == data()$search_term) %>%
            dplyr::select(cancer, log2FoldChange)
        })
        
        # The padj of the searched gene
        padj_values <- reactive({
          req(data()$DGEall, data()$search_term)
          data()$DGEall %>%
            filter(geneID == data()$search_term | gene_name == data()$search_term) %>%
            dplyr::select(cancer, padj)
        })
        
        # padj threshold value
        padj_threshold <- reactive({
          req(data()$padj_threshold)
          if (data()$padj_threshold == "NONE") {
            return(1)
          } else {
            return(data()$padj_threshold)
          }
        })
        
        # Calculate ranks of the searched gene for all cancer types, among genes with the same orientation of differential expression as the searched gene
        rank_values <- reactive({
          req(data()$DGEall)
          req(data()$padj_threshold)
          
          # Separating up and down regulated genes
          DGEall_up <- data()$DGEall %>%
            filter(log2FoldChange > 0 & padj < padj_threshold())
          DGEall_down <- data()$DGEall %>%
            filter(log2FoldChange < 0 & padj < padj_threshold())
          
          # Calculating rank
          rank_up <- DGEall_up %>%
            mutate(orientation = "up") %>%
            group_by(cancer) %>%
            arrange(desc(log2FoldChange)) %>%
            mutate(rank = row_number()) %>%
            ungroup() %>%
            filter(geneID == data()$search_term | gene_name == data()$search_term) %>%
            dplyr::select(cancer, rank, orientation)
          
          rank_down <- DGEall_down %>%
            mutate(orientation = "down") %>%
            group_by(cancer) %>%
            arrange(log2FoldChange) %>%
            mutate(rank = row_number()) %>%
            ungroup() %>%
            filter(geneID == data()$search_term | gene_name == data()$search_term) %>%
            dplyr::select(cancer, rank, orientation)
          
          # Making final table with ranks of the searched genes
          rank <- bind_rows(rank_up, rank_down)
          
          return(rank)
        })
        
        
        # Melanoma: log2FC and rank
        output$melanoma_1st_line <- renderText({
          log2fc <- log2fc_values() %>% 
            dplyr::filter(cancer == "Melanoma") %>% 
            dplyr::pull(log2FoldChange)
          sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
        })
        
        output$melanoma_2nd_line <- renderText({
          padj <- padj_values() %>% 
            dplyr::filter(cancer == "Melanoma") %>% 
            dplyr::pull(padj)
          paste0(data()$search_term, " padj: ", padj)
        })
        
        output$melanoma_3rd_line <- renderText({
          rank <- rank_values() %>% 
            dplyr::filter(cancer == "Melanoma") %>% 
            dplyr::pull(rank)
          paste0(data()$search_term, " rank among genes with padj < ", padj_threshold(), ": ", rank)
        })
        
        # Lung: log2FC and rank
        output$lung_1st_line <- renderText({
          log2fc <- log2fc_values() %>% 
            dplyr::filter(cancer == "Lung") %>% 
            dplyr::pull(log2FoldChange)
          sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
        })
        
        output$lung_2nd_line <- renderText({
          padj <- padj_values() %>% 
            dplyr::filter(cancer == "Lung") %>% 
            dplyr::pull(padj)
          paste0(data()$search_term, " padj: ", padj)
        })
        
        output$lung_3rd_line <- renderText({
          rank <- rank_values() %>% 
            dplyr::filter(cancer == "Lung") %>% 
            dplyr::pull(rank)
          paste0(data()$search_term, " rank among genes with padj < ", padj_threshold(), ": ", rank)
        })
        
        # Prostate: log2FC and rank
        output$prostate_1st_line <- renderText({
          log2fc <- log2fc_values() %>% 
            dplyr::filter(cancer == "Prostate") %>% 
            dplyr::pull(log2FoldChange)
          sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
        })
        
        output$prostate_2nd_line <- renderText({
          padj <- padj_values() %>% 
            dplyr::filter(cancer == "Prostate") %>% 
            dplyr::pull(padj)
          paste0(data()$search_term, " padj: ", padj)
        })
        
        output$prostate_3rd_line <- renderText({
          rank <- rank_values() %>% 
            dplyr::filter(cancer == "Prostate") %>% 
            dplyr::pull(rank)
          paste0(data()$search_term, " rank among genes with padj < ", padj_threshold(), ": ", rank)
        })
        
        # Glioblastoma: log2FC and rank
        output$glioblastoma_1st_line <- renderText({
          log2fc <- log2fc_values() %>% 
            dplyr::filter(cancer == "Glioblastoma") %>% 
            dplyr::pull(log2FoldChange)
          sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
        })
        
        output$glioblastoma_2nd_line <- renderText({
          padj <- padj_values() %>% 
            dplyr::filter(cancer == "Glioblastoma") %>% 
            dplyr::pull(padj)
          paste0(data()$search_term, " padj: ", padj)
        })
        
        output$glioblastoma_3rd_line <- renderText({
          rank <- rank_values() %>% 
            dplyr::filter(cancer == "Glioblastoma") %>% 
            dplyr::pull(rank)
          paste0(data()$search_term, " rank among genes with padj < ", padj_threshold(), ": ", rank)
        })
        
      }
      
    }
  )
}