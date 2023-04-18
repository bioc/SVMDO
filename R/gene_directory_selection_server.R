#' @title SVMDO
#' @name gene_directory_selection_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of entering output/working for gene list directory


innerServer <- function(input, output,session) {
  if (grepl("linux-gnu", R.version$os)) {
    osSystem<-"Linux"
  }else{
    osSystem<-R.version$os
  }
  if (osSystem == "Linux") {
    def_roots <- c(home = "~")
  }else {
    def_roots <- getVolumes()()
  }
  shinyDirChoose(
    input,
    'dir',
    roots=def_roots,
    filetypes = c("txt")
  )
  
  global <- reactiveValues(datapath=getwd())
  dir <- reactive(input$dir)
  output$dir <- renderText({
    global$datapath
  })
  
  assign("direct_val_gene_list",getwd(), envir =.GlobalEnv)
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                 if (osSystem == "Linux") {
                   filtered<-gsub("^[^/]*/", "", global$datapath)
                   assign("direct_val_gene_list",global$datapath, envir =.GlobalEnv)
                 }else{
                   filtered<-gsub("^[^/]*/", "", global$datapath)
                   selected_dir_val<-file.path(as.character(def_roots),filtered,fsep = "/")
                   global$datapath<-selected_dir_val
                   assign("direct_val_gene_list",global$datapath, envir =.GlobalEnv)

                   
                 }})

  
}
