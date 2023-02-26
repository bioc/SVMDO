#' @title SVMDO
#' @name innerServer_exp_ui
#' @param id connection input
#' @return UI section of providing expression dataset

innerUI_exp_data <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,
           fileInput(ns("file1"),
                     "Choose Your Expression Dataset",
                     accept = ".txt")))
}
