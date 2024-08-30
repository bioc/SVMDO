#' @title SVMDO
#' @name test_data_selection_ui
#' @param id connection input
#' @return UI section of providing information about selected radio button 

innerUI_test_data <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,radioButtons(ns("test_datasets"),"Load Test Dataset (DEG + Survival)",
                          choices = c("None","COAD"))))
}