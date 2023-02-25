#' @title SVMDO
#' @name top_val_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of entering top gene value

innerServer_4<-function(input, output, session) {
  top_val <- reactive({input$num_val})
  return(top_val)  
}