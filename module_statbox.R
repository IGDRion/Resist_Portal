statboxUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      summaryBox("Melanoma", 
                 HTML(paste(textOutput(ns("melanoma_log2fc")), "<br>", textOutput(ns("melanoma_mean")))), 
                 width = 6, icon = "fas fa-dna", style = "success", border = "bottom"),
      summaryBox("Lung", 
                 HTML(paste(textOutput(ns("lung_log2fc")), "<br>", textOutput(ns("lung_mean")))), 
                 width = 6, icon = "fas fa-dna", style = "success", border = "bottom")
    ),
    br(),
    fluidRow(
      summaryBox("Prostate", 
                 HTML(paste(textOutput(ns("prostate_log2fc")), "<br>", textOutput(ns("prostate_mean")))), 
                 width = 6, icon = "fas fa-dna", style = "success", border = "bottom"),
      summaryBox("Glioblastoma", 
                 HTML(paste(textOutput(ns("glioblastoma_log2fc")), "<br>", textOutput(ns("glioblastoma_mean")))), 
                 width = 6, icon = "fas fa-dna", style = "success", border = "bottom")
    )
  )
}

statboxServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      ### Calculating the data needed to go inside the box
      # The log2FC of the searched gene
      log2fc_values <- reactive({
        req(data()$DGEall, data()$search_term)
        data()$DGEall %>%
          filter(geneID == data()$search_term | gene_name == data()$search_term) %>%
          select(cancer, log2FoldChange)
      })
      
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
      
      
      ### Writing the text inside every box
      output$melanoma_log2fc <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Melanoma") %>% 
          pull(log2FoldChange)
        sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
      })
      
      output$melanoma_mean <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Melanoma") %>% 
          pull(log2FoldChange)
        
        mean_value <- cancer_means() %>% 
          filter(cancer == "Melanoma") %>% 
          pull(if(log2fc >= 0) mean_log2FC_up else mean_log2FC_down)
        
        sprintf("Mean Log2FC: %.3f", mean_value)
      })
      
      output$lung_log2fc <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Lung") %>% 
          pull(log2FoldChange)
        sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
      })
      
      output$lung_mean <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Lung") %>% 
          pull(log2FoldChange)
        
        mean_value <- cancer_means() %>% 
          filter(cancer == "Lung") %>% 
          pull(if(log2fc >= 0) mean_log2FC_up else mean_log2FC_down)
        
        sprintf("Mean Log2FC: %.3f", mean_value)
      })
      
      output$prostate_log2fc <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Prostate") %>% 
          pull(log2FoldChange)
        sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
      })
      
      output$prostate_mean <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Prostate") %>% 
          pull(log2FoldChange)
        
        mean_value <- cancer_means() %>% 
          filter(cancer == "Prostate") %>% 
          pull(if(log2fc >= 0) mean_log2FC_up else mean_log2FC_down)
        
        sprintf("Mean Log2FC: %.3f", mean_value)
      })
      
      output$glioblastoma_log2fc <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Glioblastoma") %>% 
          pull(log2FoldChange)
        sprintf(paste0(data()$search_term," Log2FC: %.3f"), log2fc)
      })
      
      output$glioblastoma_mean <- renderText({
        log2fc <- log2fc_values() %>% 
          filter(cancer == "Glioblastoma") %>% 
          pull(log2FoldChange)
        
        mean_value <- cancer_means() %>% 
          filter(cancer == "Glioblastoma") %>% 
          pull(if(log2fc >= 0) mean_log2FC_up else mean_log2FC_down)
        
        sprintf("Mean Log2FC: %.3f", mean_value)
      })
      
    }
  )
}