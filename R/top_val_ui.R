#' @title SVMDO
#' @name top_val_ui
#' @param id connection input
#' @return UI section of entering top gene value


innerUI_top_gene_val <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,      
           numericInput(ns("num_val"),"Input Size (Fixed for Test Datasets)", value = 50,min = 50)))
}