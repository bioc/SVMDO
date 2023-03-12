#' @title SVMDO
#' @name surv_plot_dw_ui
#' @param id connection input
#' @return UI section of downloading survival plots of discriminative gene set

surv_plots_download_ui<-function(id){
  ns<-NS(id)
  actionButton(ns("dw_plots"),"Download Survival Plot List")
}