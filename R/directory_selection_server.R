#' @title SVMDO
#' @name directory_selection_server
#' @param input server input
#' @param output server output
#' @param session server session
#' @return Server section of entering output/working directory


innerServer <- function(input, output,session) {
  osSystem <- Sys.info()["sysname"]
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
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                 if (osSystem == "Linux") {
                   setwd(global$datapath)
                 }else{
                   filtered<-gsub("^[^/]*/", "", global$datapath)
                   selected_dir_val<-file.path(as.character(def_roots),filtered,fsep = "/")
                   setwd(selected_dir_val)
                   global$datapath<-selected_dir_val
                 }})
}
