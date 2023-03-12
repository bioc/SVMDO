#' @title SVMDO
#' @name plot_push_ui
#' @param id connection input
#' @return UI section of providing information about total number of survival plots for visualization


innerUI_plot_inject<-function(id){
  ns <- NS(id)
  actionButton(ns("plot_add"),"Show Survival Plot List")  
}