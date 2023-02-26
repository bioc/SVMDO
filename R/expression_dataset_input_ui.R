#' @title SVMDO
#' @name innerServer_exp_ui
#' @param id connection input
#' @return UI section of providing expression dataset into GUI

# Expression dataset to be used in DEG analysis consists of 3 main columns:
# TCGA id
# Tissue type (Normal/Tumour)
# Expression values (Named with Gene Symbols)

# If there is not any requirements for survival analysis, TCGA id is not required

innerUI_exp_data <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(4,
           fileInput(ns("file1"),
                     "Choose Your Expression Dataset",
                     accept = ".txt")))
}

