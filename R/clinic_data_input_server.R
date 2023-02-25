#' @title SVMDO
#' @name clinic_data_input_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of loading clinical data 


innerServer_clinic <- function(input, output, session) {
  rawData_2 <- eventReactive(input$file2, {read.table(input$file2$datapath,sep = "\t",header = TRUE)})
  return(rawData_2)
}
