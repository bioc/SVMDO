#' @title SVMDO
#' @name runGUI
#' @param n linebreak function variable
#' @return Returning GUI window screen
#' 
#' @importFrom shinyFiles shinyDirChoose shinyDirButton getVolumes
#' @importFrom golem with_golem_options
#' @import shiny
#' @importFrom shinytitle use_shiny_title
#' @importFrom nortest ad.test
#' @importFrom SummarizedExperiment assay
#' @importFrom e1071 svm
#' @importFrom BSDA z.test
#' @importFrom data.table fread fwrite first
#' @importFrom sjmisc str_contains
#' @importFrom caTools sample.split
#' @importFrom klaR greedy.wilks
#' @importFrom caret confusionMatrix
#' @importFrom survival Surv survfit coxph survdiff
#' @importFrom DOSE enrichDO
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr mutate
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend mtext
#' @importFrom stats median na.omit p.adjust predict qnorm sd wilcox.test
#' @importFrom utils read.table write.table

#' @examples
#' #SVMDO::runGUI() Calling GUI without activating library
#' #runGUI() Calling GUI after activating library
#' # Disease Ontology Enrichment of a differentially expresed gene (entrez id):
#' a_1<-DOSE::enrichDO(2981,ont="DO")

linebreaks <- function(n){HTML(strrep(br(), n))}

outerUI <- function(id) {
  
  ns <- NS(id)
  
  navbarPage(title="SVMDO",
             tabPanel(title = "Analysis",
                      column(6,
                      fluidRow(
                        column(12,
                               wellPanel(
                                 linebreaks(1),
                                 helpText(a("Click here to examine recent version of SVMDO manual in Github", 
                                            target="_blank",
                                            href="https://github.com/robogeno/SVMDO/blob/master/vignettes/SVMDO_guide.Rmd")),
                                 p("(Wait until 'Upload complete' message appears after choosing your expression dataset)"),
                                 innerUI_exp_data(ns("inner2")),
                                 p("(In a test study, clinical and expression datasets are pre-selected)"),
                                 innerUI_test_data(ns("inner3")),
                                 p("(In an input-based differential expression analysis, an error message will appear if expression dataset is incompatible)"),
                                 innerUI_deg_analysis(ns("inner4")),
                                 linebreaks(1),
                                 innerUI_top_gene_val(ns("inner5")),
                                 innerUI_top_gene_selection(ns("inner6")),
                                 linebreaks(1),
                                 div(style="display:inline-block",
                                     innerUI_disease_ont_class(ns("inner7")),
                                     innerUI_classification(ns("inner8")
                                                            )),
                                 linebreaks(1),
                                 innerUI_clinic_data(ns("inner9")),
                                 p("(In an input-based survival analysis, incompatible clinical dataset will result in an error message)"),
                                 div(style="display:inline-block",
                                     innerUI_surv(ns("inner10")),
                                     innerUI_clear_env(ns("inner11")))
                                 ))
                                 
                               )),
                      ),
             tabPanel(title = "Results",
                      fluidRow(column(6,
                                      p("(Visualizing the gene set acquired from classification step)"),
                                      innerUI_table_show(ns("inner12")),
                                      div(dataTableOutput(ns("table")), 
                                      style = "font-size: 60%;width: 50%"),
                                      div(style="display:inline-block",
                                          p("(To download or show plots, click on Prepare Survival Plot List button)"),
                                          innerUI_collect_plot_data(ns("inner13")),
                                          div(linebreaks(1),
                                              p("(Output directory of final gene set and survival plots)"),
                                              p("(To download them in different directories, it can be used sequentially)"),
                                              innerUI_path(ns("inner1"))),
                                          innerUI_gene_names(ns("inner18")),
                                          disc_gene_download_ui(ns("inner14")),
                                          surv_plots_download_ui(ns("inner15")),
                                          innerUI_plot_inject(ns("inner16")),
                                          innerUI_plot_show(ns("inner17")))
                                      ))))}

ui <- fluidPage(
  outerUI("mod1")
)

outerServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      innerResult <- callModule(innerServer,"inner1")
      innerResult_2 <-callModule(innerServer_exp, "inner2")
      innerResult_x <- callModule(innerServer_rad,"inner3")
      innerResult_3 <- callModule(innerServer_3, "inner4",innerResult_2,innerResult_x)
      innerResult_4 <- callModule(innerServer_4, "inner5")
      innerResult_5 <- callModule(innerServer_5, "inner6",innerResult_4)
      innerResult_6 <- callModule(innerServer_6, "inner7")
      innerResult_7 <- callModule(innerServer_7, "inner8")
      innerResult_clinic <- callModule(innerServer_clinic,"inner9")
      innerResult_8 <- callModule(innerServer_8, "inner10",innerResult_clinic,innerResult_x)
      innerResult_9 <- callModule(innerServer_9,"inner11")
      innerResult_10<-callModule(table_server,"inner12")
      output$table<-renderDataTable(return(innerResult_10()),options = list(pageLength = 5)
                                           ,escape=FALSE)
      innerResult_11<-callModule(plot_list_server,"inner13")
      innerResult_16<-callModule(innerServer_10,"inner18")    
      innerResult_12<-callModule(disc_gene_dw_server,"inner14",innerResult_16)
      innerResult_13<-callModule(surv_plot_dw_server,"inner15")
      innerResult_14<-callModule(plot_push_server,"inner16")
      innerResult_15<-callModule(plot_show_server,"inner17",innerResult_14)
      }
  )
}

server <- function(input, output, session) {
  outerServer("mod1")
  session$onSessionEnded(function() {
    stopApp()
  })
}

#' @export

runGUI <- function(
    reproducible = TRUE
) {
  options(shiny.maxRequestSize = 700*1024^2)
  with_golem_options(
    app = shinyApp(
      ui = ui, 
      server = server,
      options = list("launch.browser"=TRUE)
    ),
    golem_opts = list(reproducible = reproducible)
  )
}

