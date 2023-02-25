#' @title SVMDO
#' @name top_val_based_deg_filtration_ui
#' @param id connection input
#' @return UI section of selecting differentially expressed genes based on top gene value

innerUI_top_gene_selection <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(1,      
           actionButton(ns("initiate_top_gene_selection"),"Top Gene Number Selection")))
}
