#' @title SVMDO
#' @name disc_gene_dw_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of final discriminative gene set download button

disc_gene_dw_server<-function(input,output,session){
  observeEvent(input$dw_genes, {
    if(exists("final_discriminative_gene_set")){
      write.table(final_discriminative_gene_set,"final_discriminative_gene_set.txt",sep = "\t")
      showModal(
        modalDialog(
          title = "Gene List Download Completed",
          paste("Save Location:",sep=" ",getwd()),
          easyClose = TRUE,
          footer = NULL
        )
      )
    }else{
      showModal(
        modalDialog(
          title = "Gene List Download Failed",
          "Gene List Not Found",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
}