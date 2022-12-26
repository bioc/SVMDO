globalVariables(c("final_discriminative_gene_set","top_genes","tcga_id_list", 
                  "..elected_val","tissue_type_list","disease_filtered_gene_data",
                  "tcga_sample_comb","new_tissue_type_list","tissue_type"))

if (!requireNamespace("nortest", quietly = TRUE))
  stop("Install 'nortest' to use this package.")

if (!requireNamespace("e1071", quietly = TRUE))
  stop("Install 'e1071' to use this package.")

if (!requireNamespace("BSDA", quietly = TRUE))
  stop("Install 'BSDA' to use this package.")

if (!requireNamespace("data.table", quietly = TRUE))
  stop("Install 'data.table' to use this package.")

if (!requireNamespace("sjmisc", quietly = TRUE))
  stop("Install 'sjmisc'  to use this package.")

if (!requireNamespace("klaR", quietly = TRUE))
  stop("Install 'klaR' to use this package.")

if (!requireNamespace("caTools", quietly = TRUE))
  stop("Install 'caTools' to use this package.")

if (!requireNamespace("caret", quietly = TRUE))
  stop("Install 'caret' to use this package.")

if (!requireNamespace("survival", quietly = TRUE))
  stop("Install 'survival' to use this package.")

if (!requireNamespace("DOSE", quietly = TRUE))
  stop("Install 'DOSE' to use this package.")

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  stop("Install 'AnnotationDbi' to use this package.")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  stop("Install 'org.Hs.eg.db' to use this package.")

if (!requireNamespace("dplyr", quietly = TRUE))
  stop("Install 'dplyr' to use this package.")

if (!requireNamespace("shiny", quietly = TRUE))
  stop("Install 'shiny' to use this package.")

if (!requireNamespace("rstudioapi", quietly = TRUE))
  stop("Install 'rstudioapi' to use this package.")

#' Initiating GUI screen
#' 
#' Function of SVMDO initiating GUI screen of main dialog box including
#' analysis steps
#' 
#' @param reproducible Visualizing GUI screen

#' @return Returning GUI window screen
#' @export

#' @importFrom rstudioapi selectDirectory
#' @importFrom golem with_golem_options 
#' @importFrom nortest ad.test
#' @importFrom e1071 svm
#' @importFrom BSDA z.test
#' @importFrom data.table fread fwrite
#' @importFrom sjmisc str_contains
#' @importFrom caTools sample.split
#' @importFrom klaR greedy.wilks
#' @importFrom caret confusionMatrix
#' @importFrom survival Surv survfit coxph
#' @importFrom DOSE enrichDO
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr mutate
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend mtext
#' @importFrom stats median na.omit p.adjust predict qnorm sd wilcox.test
#' @importFrom utils read.table write.table


#' @import shiny

#' @examples
#' #SVMDO::runGUI() Calling GUI without activating library
#' #runGUI() Calling GUI after activating library
#' # Disease Ontology Enrichment of a differentially expresed gene (entrez id):
#' a_1<-DOSE::enrichDO(2981,ont="DO", maxGSSize=Inf)


runGUI <- function(
    reproducible = TRUE
) {
  options(shiny.maxRequestSize = 700*1024^2)
  with_golem_options(
    app = shinyApp(
      ui = ui, 
      server = server
    ),
    golem_opts = list(reproducible = reproducible)
  )
}