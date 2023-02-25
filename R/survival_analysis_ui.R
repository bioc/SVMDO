#' @title SVMDO
#' @name survival_analysis_ui
#' @param id connection input
#' @return UI section of survival analysis of final discriminative gene set


innerUI_surv<-function(id){
  ns <- NS(id)
  actionButton(ns("initiate_surv_analysis"),"Survival Analysis")
}
