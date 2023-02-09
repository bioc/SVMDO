globalVariables(c("final_discriminative_gene_set","top_genes","tcga_id_list", 
                  "..elected_val","..num_data", "..type_num_data","tissue_type_list","disease_filtered_gene_data",
                  "tcga_sample_comb","new_tissue_type_list","tissue_type",
                  "top_genes_test","total_exp_dataset"))

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

if (!requireNamespace("shinytitle", quietly = TRUE))
  stop("Install 'shiny' to use this package.")

if (!requireNamespace("shinyFiles", quietly = TRUE))
  stop("Install 'shinyFiles' to use this package.")

if (!requireNamespace("golem", quietly = TRUE))
  stop("Install 'shiny' to use this package.")

#' Initiating GUI screen
#' 
#' Function of SVMDO initiating GUI screen of main dialog box including
#' analysis steps
#' 
#' @param reproducible Visualizing GUI screen

#' @return Returning GUI window screen
#' @export

#' @importFrom golem with_golem_options
#' @import shiny

#' @examples
#' #SVMDO::runGUI() Calling GUI without activating library
#' #runGUI() Calling GUI after activating library
#' # Disease Ontology Enrichment of a differentially expresed gene (entrez id):
#' a_1<-DOSE::enrichDO(2981,ont="DO")


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
    golem_opts = list(reproducible = reproducible),
  )
}