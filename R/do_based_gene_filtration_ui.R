#' @title SVMDO
#' @name do_based_gene_filtration_ui
#' @param id connection input
#' @return UI section of disease ontology based filtration of differentially expressed genes

innerUI_disease_ont_class<- function(id){
  ns <- NS(id)
  actionButton(ns("initiate_do_analysis"),"DO Analysis")
}
