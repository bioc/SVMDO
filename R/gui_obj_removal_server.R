#' @title SVMDO
#' @name gui_obj_removal_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of workspace clearance


innerServer_9<-function(input,output,session){
  observeEvent(input$clean_workspace, {
    string_names_1<-c("complete_deg_gene_list","complete_deg_gene_list_test","disease_filtered_gene_data",
                    "final_discriminative_gene_set","new_tissue_type_list",
                    "sorted_new_bound_form_A","sorted_new_bound_form_A_test",
                    "sorted_new_bound_form_B", "sorted_new_bound_form_B_test",
                    "tcga_id_list","tcga_sample_comb","tissue_type_list",
                    "top_genes_test","top_genes","total_exp_dataset",
                    "total_exp_dataset_test","max_plots")
    
    n_list<-NULL
    y_list<-NULL
    sig_val<-NULL
    
    rm_var_names_1<-lapply(seq.int(string_names_1), function(i){
      if (exists(string_names_1[i],envir = .GlobalEnv )) {
        rm(list = c(string_names_1[i]), envir = .GlobalEnv)
        y_list<-c(0,y_list)
      }else{
        n_list<-c(1,n_list)
      }
      return(list(y_list=y_list,n_list=n_list))
      })
    
    if (all(names(table(unlist(rm_var_names_1), useNA="ifany")) %in% 1)) {
      y_list<-0
      n_list<-length(string_names_1)
    }else{
      y_list<-as.numeric(table(unlist(rm_var_names_1))[1])
      
      n_list<-as.numeric(table(unlist(rm_var_names_1))[2])

    }

    if(y_list==0){
      sig_val<-0
    }
    if(y_list>0){
      if ((n_list+y_list)==length(string_names_1)) {
        sig_val<-1
      }
    }
    
    if (length(ls(,pattern = "^fit1_",envir = .GlobalEnv))>0) {
      a<-length(ls(,pattern = "^fit1_",envir = .GlobalEnv))
      b<-ls(,pattern = "^fit1_",envir = .GlobalEnv)
      c<-ls(,pattern= "^modulename_",envir = .GlobalEnv)
      d<-ls(,pattern= "^p_",envir = .GlobalEnv)
      e<-ls(,pattern= "^hr_",envir = .GlobalEnv)
      rm_var_names_fit<-lapply(seq.int(a), function(i){
        rm(list = b[i], envir = .GlobalEnv)})
      rm_var_names_mod<-lapply(seq.int(a), function(i){
        rm(list = c[i], envir = .GlobalEnv)})
      rm_var_names_pval<-lapply(seq.int(a), function(i){
        rm(list = d[i], envir = .GlobalEnv)})
      rm_var_names_hr<-lapply(seq.int(a), function(i){
        rm(list = e[i], envir = .GlobalEnv)})
      }
      
      

    
    if (exists("plot_prep_sign",envir = .GlobalEnv )) {
      rm(list = c("fit_list","hr_list","p_list","modulename_list","plot_prep_sign"), envir = .GlobalEnv)
    }
    

    
    if (sig_val==1) {
      showModal(
        modalDialog(
          title = "Clearing SVMDO Objects",
          "Process Completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    if(sig_val==0){
      showModal(
        modalDialog(
          title = "Error in SVMDO Objects",
          "Objects were not found",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
}