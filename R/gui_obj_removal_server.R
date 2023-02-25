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
    
    string_names_2<-c("fit_list","hr_list","p_list","modulename_list")
    
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

    if (!is.null(ls(pattern = "fit1_",envir = .GlobalEnv))) {
      fit_list<-ls(pattern = "fit1_",envir = .GlobalEnv)
      modulename_list<-ls(pattern = "modulename_",envir = .GlobalEnv)
      p_list<-ls(pattern = "^p_",envir = .GlobalEnv)
      hr_list<-ls(pattern = "hr_",envir = .GlobalEnv)
      
      assign("fit_list",fit_list,envir = .GlobalEnv)
      assign("modulename_list",modulename_list,envir = .GlobalEnv)
      assign("hr_list",hr_list,envir = .GlobalEnv)
      assign("p_list",p_list,envir = .GlobalEnv)
    }
    
    rm_fit_names<-lapply(seq.int(fit_list), function(i){
      if (exists(fit_list[i],envir = .GlobalEnv )) {
        rm(list = c(fit_list[i]), envir = .GlobalEnv)
        }
      })
    
    rm_p_names<-lapply(seq.int(p_list), function(i){
      if (exists(p_list[i],envir = .GlobalEnv )) {
        rm(list = c(p_list[i]), envir = .GlobalEnv)
      }
    })
    
    rm_hr_names<-lapply(seq.int(hr_list), function(i){
      if (exists(hr_list[i],envir = .GlobalEnv )) {
        rm(list = c(hr_list[i]), envir = .GlobalEnv)
      }
    })
    
    rm_mod_names<-lapply(seq.int(modulename_list), function(i){
      if (exists(modulename_list[i],envir = .GlobalEnv )) {
        rm(list = c(modulename_list[i]), envir = .GlobalEnv)
      }
    })
    
    rm_var_names_2<-lapply(seq.int(string_names_2), function(i){
      if (exists(string_names_2[i],envir = .GlobalEnv )) {
        rm(list = c(string_names_2[i]), envir = .GlobalEnv)
      }})
    
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