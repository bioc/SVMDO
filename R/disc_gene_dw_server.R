#' @title SVMDO
#' @name disc_gene_dw_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @param gene_list_val discriminative gene set list variable
#' @return Server section of discriminative gene set download button

disc_gene_dw_server<-function(input,output,session,gene_list_val){
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
            title = "Filename Error",
            "Please Enter Filename",
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    }else{
      showModal(
        modalDialog(
          title = "Gene List Download Failed",
          "Missing Gene List Data",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
  
}