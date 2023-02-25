#' @title SVMDO
#' @name gui_obj_removal_ui
#' @param id connection input
#' @return UI section of workspace clearance

innerUI_clear_env<- function(id){
  ns <- NS(id)
  actionButton(ns("clean_workspace"),"Clear Environment")
}