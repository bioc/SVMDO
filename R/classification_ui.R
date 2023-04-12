#' @title SVMDO
#' @name classification_ui
#' @param id connection input
#' @return UI section of wilks lambda filtration and SVM classification of disease filtered differentially expressed gene set


innerUI_classification<- function(id){
  ns <- NS(id)
  actionButton(ns("initiate_wlsvmdo_analysis"),"Classification Analysis")
}
