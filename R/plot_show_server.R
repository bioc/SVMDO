#' @title SVMDO
#' @name plot_show_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @param max_data Information of total number of survival plots prepared with discriminative gene set
#' @return Server section of providing information about total number of survival plots for visualization

plot_show_server<-function(input,output,session,max_data){
  output$plots<-renderUI({
    max_plot_val<-max_data()
    lapply(seq.int(max_plot_val),function(a) {
      local({
        plotname <- paste("plot", a, sep="")
        output[[plotname]] <- renderPlot({
          plot(get(fit_list[a]),main=get(modulename_list[a]),col = c(2,4),frame.plot=0,ylab="Survival Probability",xlab="Time (day)",cex.axis=1.5,cex.lab=1.5,cex.main=2)
          legend("top",legend = c("High risk group", "Low risk group"),text.col = c(2,4),bty = "n",horiz = TRUE,cex = 1.5)
          mtext(text = paste0("p = ",get(p_list[a])),side=3,adj = 1,line = -6,cex=1.5)
          mtext(text = paste0("HR = ",get(hr_list[a])),side=3,adj = 1,line = -8,cex = 1.5)
        },height = 500,width=500        
        )
      })
    })    
  })
}