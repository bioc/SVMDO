#' @title SVMDO
#' @name directory_selection_ui
#' @param id connection input
#' @return UI section of entering output/working directory

innerUI_path <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(5,
           shinyDirButton(ns("dir"), "Choose Directory", ns("Upload")),
           verbatimTextOutput(ns("dir"), placeholder = TRUE)))
}