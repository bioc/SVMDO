#' @title SVMDO
#' @name gene_list_table_visualization_ui
#' @param id connection input
#' @return Providing table form of discriminative gene sets in GUI

test_data_table_ui  <- function(id){
  ns <- NS(id)
  tagList(
    dataTableOutput(outputId = ns("table")),
    textOutput(outputId = ns("my_text"))  
  )  
}
