#' @title SVMDO
#' @name table_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of providing discriminative gene set for preparing table

table_server<-function(input,output,session){
  getData<-eventReactive(input$file,{
    if (exists("final_discriminative_gene_set")) {
      final_discriminative_gene_set
    }else{
      showModal(
        modalDialog(
          title = "Error in Gene Set Visualization",
          "Classification Step is Missing",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
}