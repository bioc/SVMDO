#' @title SVMDO
#' @name table_ui
#' @param id connection input
#' @return UI section of providing discriminative gene set for preparing table


innerUI_table_show <- function(id) {
  ns <- NS(id)
  actionButton(ns("file"),"Show Gene Results")}