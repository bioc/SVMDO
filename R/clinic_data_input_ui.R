#' @title SVMDO
#' @name clinic_data_input_ui
#' @param id connection input
#' @return UI section of loading clinical data 

# Clinical dataset to be used in survival analysis consists of 3 main columns:
# TCGA id
# Patient Days to Death Information
# Patient Vital Status (Death/Alive)


innerUI_clinic_data <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,
           fileInput(ns("file2"),
                     "Choose Clinical Data",
                     accept = ".txt")))
}