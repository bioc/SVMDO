#' @title SVMDO
#' @name disc_gene_dw_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of final discriminative gene set download button

disc_gene_dw_server<-function(input,output,session,gene_list_val,global){
  observeEvent(input$dw_genes, {
    if(exists("final_discriminative_gene_set")){
      gene_name_variable<-gene_list_val()
      if ((gene_name_variable)!="") {
        path_loc<-file.path(direct_val_gene_list,"/",gene_name_variable,".txt",fsep="")
        write.table(final_discriminative_gene_set,file = path_loc,sep = "\t")
        showModal(
          modalDialog(
            title = "Gene List Download Completed",
            paste("Save Location:",sep=" ",direct_val_gene_list),
            easyClose = TRUE,
            footer = NULL
          )
        )
      }else{
        showModal(
          modalDialog(
            title = "No variable name entered",
            "Please enter gene list name",
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    }else{
      showModal(
        modalDialog(
          title = "FAILED",
          "MISSING STEP",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
}