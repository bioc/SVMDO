#' @title SVMDO
#' @name gene_directory_selection_ui
#' @param id connection input
#' @return UI section of entering output/working for gene list directory

innerUI_path <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(5,
           shinyDirButton(ns("dir"), "Choose Output Directory", ns("Upload")),
           verbatimTextOutput(ns("dir"), placeholder = TRUE)))
}