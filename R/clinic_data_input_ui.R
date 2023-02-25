#' @title SVMDO
#' @name clinic_data_input_ui
#' @param id connection input
#' @return UI section of loading clinical data 

innerUI_clinic_data <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,
           fileInput(ns("file2"),
                     "Choose Clinical Data",
                     accept = ".txt")))
}