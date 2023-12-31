#' @title SVMDO
#' @name plot_list_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of preparing plot list to be visualized in GUI page

plot_list_server<-function(input,output,session){
  observeEvent(input$list_data, {
    
    if (length(ls(pattern = "^fit1_",envir = .GlobalEnv))>0) {
      fit_list<-ls(pattern = "^fit1_",envir = .GlobalEnv)

      modulename_list<-ls(pattern = "modulename_",envir = .GlobalEnv)
      p_list<-ls(pattern = "^p_",envir = .GlobalEnv)
      hr_list<-ls(pattern = "^hr_",envir = .GlobalEnv)
      
      assign("fit_list",fit_list,envir = .GlobalEnv)
      assign("modulename_list",modulename_list,envir = .GlobalEnv)
      assign("hr_list",hr_list,envir = .GlobalEnv)
      assign("p_list",p_list,envir = .GlobalEnv)
      assign("max_plots",length(fit_list),envir = .GlobalEnv)
      assign("plot_prep_sign",1,envir = .GlobalEnv)
      
      showModal(
        modalDialog(
          title = "Survival Plot Preparation",
          "Process Completed",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }else{
      showModal(
        modalDialog(
          title = "Survival Plot Preparation",
          "Missing Survival Plots",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
    })}