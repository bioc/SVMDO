#' @title SVMDO
#' @name top_val_based_deg_filtration
#' @param input server input
#' @param output server output
#' @param session server session
#' @param top_val top gene number value provided from top_val_server
#' @return Server section of selecting differentially expressed genes based on top gene value

innerServer_5<- function(input,output,session,top_val){
  observeEvent(input$initiate_top_gene_selection, {
    top_gene_selection<-NULL
    top_gene_number<-top_val()

    if (exists("sorted_new_bound_form_A")) {
      
      if (nrow(sorted_new_bound_form_A) < top_gene_number | nrow(sorted_new_bound_form_B) < top_gene_number) {
        max_down_genes<-sorted_new_bound_form_A
        max_up_genes<-sorted_new_bound_form_B
        message("Insufficient gene number All up/downregulated genes were selected")
        message("Enter a lower value of input size")
      }else{
        max_down_genes<-sorted_new_bound_form_A[seq.int((nrow(sorted_new_bound_form_A)-(top_gene_number-1)),nrow(sorted_new_bound_form_A)),]
        max_up_genes<-sorted_new_bound_form_B[seq.int((nrow(sorted_new_bound_form_B)-(top_gene_number-1)),nrow(sorted_new_bound_form_B)),]
      }
      top_combined_genes<-rbind(max_down_genes,max_up_genes)
      rownames(top_combined_genes)<-top_combined_genes[,1]
      changed_whole_data<-subset(total_exp_dataset,select=top_combined_genes$Genes)
      assign("top_genes",changed_whole_data,envir =.GlobalEnv)
      top_gene_selection<-1
    }else{
      showModal(
        modalDialog(
          title = "Error in Top Gene Selection",
          "Process Failed",
          easyClose = TRUE,footer = NULL))}
    
    if (!is.null(top_gene_selection)) {
      showModal(
        modalDialog(
          title = "Top Gene Selection Result",
          "Process Completed",
          easyClose = TRUE,footer = NULL))}
  })
}