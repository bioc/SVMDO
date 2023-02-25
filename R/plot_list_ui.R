#' @title SVMDO
#' @name plot_list_ui
#' @param id connection output
#' @return UI section of preparing plot list to be visualized in GUI page

innerUI_collect_plot_data<-function(id){
  ns <- NS(id)
  actionButton(ns("list_data"),"Prepare Plot Lists")
}