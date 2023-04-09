#' @title SVMDO
#' @name gene_list_name_ui
#' @param id connection input
#' @return UI section of entering top gene value


innerUI_gene_names <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,      
           textInput(ns("gene_names"),"Enter gene names")))
}