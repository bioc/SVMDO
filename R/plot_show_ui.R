#' @title SVMDO
#' @name plot_show_ui
#' @param id connection input
#' @return UI section of providing information about total number of survival plots for visualization


innerUI_plot_show<-function(id){
  ns <- NS(id)
  uiOutput(outputId = ns("plots"))
}