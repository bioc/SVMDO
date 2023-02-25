#' @title SVMDO
#' @name plot_push_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of providing information about total number of survival plots for visualization

plot_push_server<-function(input,output,session){
  max_data<-eventReactive(input$plot_add,{
    max_plots
  })
  return(max_data)
}