#' @title SVMDO
#' @name expression_dataset_input_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of providing expression dataset

# Expression dataset to be used in DEG analysis consists of 3 main columns:
# TCGA id
# Tissue type (Normal/Tumour)
# Expression values (Named with Gene Symbols)

# If there is not any requirements for survival analysis, TCGA id is not required

innerServer_exp <- function(input, output, session) {
  rawData <- eventReactive(input$file1, {
    fread(input$file1$datapath,sep = "\t")
    })
  return(rawData)
  
}