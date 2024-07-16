library(shiny)

# UI module example
searchBarUI <- function(id){ # ADD ARGUMENTS IF NEEDED
  ns <- NS(id)
  
  #plotOutput(outputId = ns("plot"))
  # UI HERE
  
}


# Server module example
searchBarServer <- function(id) { # ADD ARGUMENTS IF NEEDED
  moduleServer(
    id = id,
    module = function(input, output, session) {
        
      # myplot <- ggplot(...)
      # LOGIC HERE
      
    }
  )
}
