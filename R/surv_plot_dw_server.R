#' @title SVMDO
#' @name surv_plot_dw_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of downloading survival plots of discriminative gene set


surv_plot_dw_server<-function(input,output,session){
  observeEvent(input$dw_plots,{
    if(exists("final_discriminative_gene_set")){
      lapply(seq.int(length(fit_list)),function(a){
        png(filename=paste0(get(modulename_list[a]),".png"),width = 800,height = 600)
        plot(get(fit_list[a]),main=get(modulename_list[a]),col = c(2,4),frame.plot=0,ylab="Survival Probability",xlab="Time (day)",cex.axis=1.5,cex.lab=1.5,cex.main=2)
        legend("top",legend = c("High risk group", "Low risk group"),text.col = c(2,4),bty = "n",horiz = TRUE,cex = 1.5)
        mtext(text = paste0("p = ",get(p_list[a])),side=3,adj = 1,line = -6,cex=1.5)
        mtext(text = paste0("HR = ",get(hr_list[a])),side=3,adj = 1,line = -8,cex = 1.5)
        dev.off()})
      showModal(
        modalDialog(
          title = "Plot List Download Completed",
          paste("Save Location:",sep=" ",getwd()),
          easyClose = TRUE,
          footer = NULL
        )
      )
      
    }else{
      showModal(
        modalDialog(
          title = "Survival Plot List Download Failed",
          "Plot List Not Found",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })
}