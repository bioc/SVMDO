#' @title SVMDO
#' @name gene_list_name_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of entering final gene list name 

innerServer_10<-function(input, output, session) {
  gene_list_val <- reactive({input$gene_names})
  return(gene_list_val)
}