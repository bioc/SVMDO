#' @title SVMDO
#' @name expression_dataset_input_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of providing expression dataset

innerServer_exp <- function(input, output, session) {
  rawData <- eventReactive(input$file1, {fread(input$file1$datapath,sep = "\t")})
  return(rawData)
}