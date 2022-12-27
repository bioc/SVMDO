linebreaks <- function(n){HTML(strrep(br(), n))}
ui <- fluidPage(
  
  titlePanel("SVMDO Analysis"),
  fluidRow(
    column(3,
           actionButton("dir_select","Choose directory",
                        icon = icon("folder")
           ),htmlOutput("dir_select_ready"))),
  
  linebreaks(1),
  
  fluidRow(
    column(4,
           fileInput("file1", "Choose Your Expression Dataset",accept = ".txt"))),
  fluidRow(
    column(4,radioButtons("test_datasets","Load Test Datasets (DEG + Survival)",choices = c("None","COAD","LUSC")))),
  
  fluidRow(
    column(1,      
           actionButton("initiate_deg_analysis", "DEG Analysis"))),
  
  linebreaks(1),
  
  fluidRow(
    column(5,
           numericInput("num_val","Input Size (Fixed for Test Datasets)", value = 50,min = 50))),
  
  
  fluidRow(
    column(5,
           actionButton("initiate_top_gene_selection","Top Gene Number Selection"))),
  
  linebreaks(1),
  
  fluidRow(
    column(1,       
           actionButton("initiate_do_analysis","DO Analysis"))),
  
  linebreaks(1),
  
  fluidRow(
    column(1,       
           actionButton("initiate_wlsvmdo_analysis","Classification"))),
  
  linebreaks(1),
  
  fluidRow(
    column(4,
           fileInput("file2", "Choose Clinical Data",accept = ".txt"))),
  fluidRow(
    column(1,
           actionButton("initiate_surv_analysis","Survival Analysis"))),
  mainPanel(
  )
)
