#' @title SVMDO
#' @name deg_ui
#' @param id connection input

#' @return UI section of differential gene expression analysis

innerUI_deg_analysis <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(1,      
           actionButton(ns("initiate_deg_analysis"), "DEG Analysis")))
}