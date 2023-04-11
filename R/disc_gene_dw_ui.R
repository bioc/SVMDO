#' @title SVMDO
#' @name disc_gene_download_ui
#' @param id connection input
#' @return UI section of discriminative gene set download button

disc_gene_download_ui<- function(id){
  ns <- NS(id)
  actionButton(ns("dw_genes"),"Download Gene List")
}