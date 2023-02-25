#' @title SVMDO
#' @name test_data_selection_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of providing information about selected radio button 

innerServer_rad<-function(input, output, session) {
  rval <- reactive({input$test_datasets})
  return(rval)
}
