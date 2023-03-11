#' @title SVMDO
#' @name gui_obj_removal_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of workspace clearance


innerServer_9<-function(input,output,session){
  observeEvent(input$clean_workspace, {
    string_names_1<-c("complete_deg_gene_list","disease_filtered_gene_data",
                    "final_discriminative_gene_set","new_tissue_type_list",
                    "sorted_new_bound_form_A","sorted_new_bound_form_B",
                    "tcga_id_list","tcga_sample_comb","tissue_type_list",
                    "top_genes_test","top_genes","total_exp_dataset",
                    "max_plots")
    
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

    if (length(ls(pattern = "fit1_",envir = .GlobalEnv))>0) {
      
      pat_list<-c("fit1_","modulename_","p_","hr_")
      
      rm_pat_list<-lapply(seq.int(pat_list),function(i){
        rm(list = pat_list[grep(file.path("^",pat_list[i],fsep = ""), pat_list)],envir = .GlobalEnv)
        })
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